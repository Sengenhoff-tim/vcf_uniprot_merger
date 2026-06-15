use regex::Regex;
use std::collections::{HashMap, HashSet};
use std::io::{BufRead, Write};

use crate::variant::Variant;

use anyhow::{Context, Result};

use std::sync::LazyLock;

static REGEX_ENST: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r"^ENST\d{11}(?:\.\d+)?$").unwrap());

static REGEX_AA_CHANGE: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r"^([A-Z]+) -> ([A-Z]+)").unwrap());

/// A feature line block read from the input, kept as raw text so it can be
/// written back out unchanged (or only lightly modified).
enum FtBlock {
    /// Any non-variant feature (SIGNAL, CHAIN, MOD_RES, ...). Just text.
    NonVariant(String),
    /// A VARIANT / VAR_SEQ feature. We keep the original text *and* a parsed
    /// `Variant` purely for comparison against insert candidates. `matched_ids`
    /// collects the ids of any candidates that matched this feature.
    Variant {
        variant: Variant,
        lines: String,
        matched_ids: Vec<String>,
    },
}

/// Scratch state while a single FT feature block is being read across lines.
#[derive(Default)]
struct PendingFt {
    active: bool,
    is_variant: bool,
    /// set when the position could not be parsed (e.g. crosslinks); the block
    /// is discarded on flush, matching the previous silent-skip behaviour.
    drop: bool,
    /// raw text of every line in this block, verbatim.
    lines: String,
    /// accumulated /note= text, used only to parse aa_ref/aa_new for comparison.
    note: String,
    pos: Option<(u32, Option<u32>)>,
}

pub fn add_variants<R, W>(
    reader: R,
    writer: &mut W,
    features: &HashMap<u32, Vec<Variant>>,
    confirmed_only: bool,
) -> Result<HashSet<u32>>
where
    R: BufRead,
    W: Write,
{
    let mut insert_candidates: Vec<&[Variant]> = Vec::new();
    let mut global_features_inserted: HashSet<u32> = HashSet::new();

    // Existing FT feature blocks of the current entry, in file order.
    let mut blocks: Vec<FtBlock> = Vec::new();
    let mut pending = PendingFt::default();

    let mut seq = String::new();
    let mut seq_len: u32 = 0;
    let mut collect_seq = false;

    for (line_num, line_result) in reader.lines().enumerate() {
        let cur_line = line_result.context(format!("Failed to read cur_line {}", line_num + 1))?;

        if cur_line.starts_with("ID")
            || cur_line.starts_with("AC")
            || cur_line.starts_with("CC")
            || cur_line.starts_with("DT")
            || cur_line.starts_with("DE")
            || cur_line.starts_with("GN")
            || cur_line.starts_with("OS")
            || cur_line.starts_with("OC")
            || cur_line.starts_with("OX")
            || cur_line.starts_with("RN")
            || cur_line.starts_with("RP")
            || cur_line.starts_with("RG")
            || cur_line.starts_with("RA")
            || cur_line.starts_with("RT")
            || cur_line.starts_with("RL")
            || cur_line.starts_with("PE")
            || cur_line.starts_with("KW")
        {
            writeln!(writer, "{}", cur_line)?;
            continue;
        }

        if collect_seq && !cur_line.starts_with("//") {
            seq.push_str(&format!("{}\n", &cur_line));
        }

        if cur_line.starts_with("SQ") {
            // FT section is finished; finalise any pending FT block.
            flush_pending_ft(&mut blocks, &mut pending);

            seq_len = cur_line[13..]
                .split_whitespace()
                .next()
                .context(format!("Malforemd SQ cur_line {}", line_num + 1))?
                .parse::<u32>()
                .context(format!("Malforemd SQ cur_line {}", line_num + 1))?;

            seq.push_str(&format!("{}\n", &cur_line));

            collect_seq = true;

            continue;
        }

        if cur_line.starts_with("//") {
            // Defensive: flush in case there was no SQ line.
            flush_pending_ft(&mut blocks, &mut pending);

            let last_aa = seq.trim_end().chars().last();

            if let Some(last_aa_valid) = last_aa {
                // Index existing variant blocks for lookup. Equality/hash on
                // Variant ignores the id, so candidates match on position + aa.
                let mut index: HashMap<Variant, usize> = HashMap::new();
                for (i, block) in blocks.iter().enumerate() {
                    if let FtBlock::Variant { variant, .. } = block {
                        index.insert(variant.clone(), i);
                    }
                }

                for candidate_slice in insert_candidates.drain(..) {
                    for mut candidate in candidate_slice.iter().cloned() {
                        candidate.normalize(seq_len, last_aa_valid);

                        if let Some(&i) = index.get(&candidate) {
                            // Matches an existing feature: keep that feature's
                            // text and just record the candidate id to merge in.
                            if let FtBlock::Variant { matched_ids, .. } = &mut blocks[i] {
                                if !candidate.id.is_empty() {
                                    matched_ids.push(candidate.id.clone());
                                }
                            }
                        } else {
                            // Brand-new variant from the external source.
                            write!(writer, "{}", candidate.to_uniprot(seq_len)?)?;
                        }
                    }
                }
            } else {
                insert_candidates.clear();
            }

            // Write existing feature blocks back out in their original order.
            for block in blocks.drain(..) {
                match block {
                    FtBlock::NonVariant(text) => {
                        if !confirmed_only {
                            write!(writer, "{}", text)?;
                        }
                    }
                    FtBlock::Variant {
                        lines, matched_ids, ..
                    } => {
                        if matched_ids.is_empty() {
                            // Existing-only variant: kept unless confirmed_only.
                            if !confirmed_only {
                                write!(writer, "{}", lines)?;
                            }
                        } else {
                            // Confirmed by a candidate: keep the text, merge ids.
                            write!(writer, "{}", merge_ids_into_block(&lines, &matched_ids))?;
                        }
                    }
                }
            }

            //TODO add failure logging
            writeln!(writer, "{}{}", seq, cur_line)?;

            seq.clear();
            collect_seq = false;

            continue;
        }

        if cur_line.starts_with("DR") {
            if cur_line.starts_with("DR   Ensembl;") {
                collect_enst_variants(
                    &cur_line,
                    &REGEX_ENST,
                    features,
                    &mut global_features_inserted,
                    &mut insert_candidates,
                )
                .context(format!("Failed to parse DR cur_line {}", line_num + 1))?;
            }
            writeln!(writer, "{}", cur_line)?;
            continue;
        }

        if cur_line.starts_with("FT") {
            let content = cur_line.get(21..).unwrap_or("").trim_start();

            let ft_field = cur_line
                .get(2..21)
                .context(format!("Failed to read FT cur_line {}", line_num + 1))?
                .trim();

            if !ft_field.is_empty() {
                // A new feature header begins: finalise the previous block.
                flush_pending_ft(&mut blocks, &mut pending);

                pending.active = true;
                pending.lines.push_str(&cur_line);
                pending.lines.push('\n');

                if matches!(ft_field, "VARIANT" | "VAR_SEQ") {
                    pending.is_variant = true;

                    if content.starts_with('<') {
                        anyhow::bail!(format!("Failed to read FT cur_line {}", line_num + 1));
                    }

                    let mut parts = content.split("..");
                    let pos_start = parts.next().unwrap_or("").trim();
                    let pos_end = parts.next();

                    match pos_start.parse::<u32>() {
                        Ok(pos_start_u32) => {
                            let pos_end_u32 = match pos_end.map(|p| p.trim().parse::<u32>()) {
                                Some(Ok(val)) => Some(val),
                                //TODO maybe add logging
                                Some(Err(_)) => {
                                    pending.drop = true;
                                    None
                                }
                                None => None,
                            };
                            if !pending.drop {
                                pending.pos = Some((pos_start_u32, pos_end_u32));
                            }
                        }
                        //TODO maybe add logging. this currently skips crosslinks silently
                        Err(_) => {
                            pending.drop = true;
                        }
                    }
                }

                continue;
            }

            // Continuation line of the current feature (ft_field is blank).
            if pending.active {
                pending.lines.push_str(&cur_line);
                pending.lines.push('\n');

                if pending.is_variant {
                    if content.starts_with("/note=") {
                        pending.note = cur_line.clone();
                    } else if !content.starts_with('/') {
                        // wrapped continuation of the previous qualifier (note)
                        pending.note.push_str(cur_line.get(21..).unwrap_or(""));
                    }
                    // /evidence= and /id= are preserved in `lines`, not parsed.
                }
            }

            continue;
        }
    }

    Ok(global_features_inserted)
}

/// Finalise the in-progress FT block (if any) and push it onto `blocks`.
fn flush_pending_ft(blocks: &mut Vec<FtBlock>, pending: &mut PendingFt) {
    if pending.active {
        if pending.is_variant {
            if !pending.drop {
                if let Some(pos) = pending.pos {
                    // Build a Variant for comparison; fall back to keeping the
                    // text as a non-comparable block if there is no parseable note.
                    let parsed = extract_note(&pending.note).map(|note| build_variant(note, pos));
                    match parsed {
                        Some(variant) => blocks.push(FtBlock::Variant {
                            variant,
                            lines: std::mem::take(&mut pending.lines),
                            matched_ids: Vec::new(),
                        }),
                        None => {
                            blocks.push(FtBlock::NonVariant(std::mem::take(&mut pending.lines)))
                        }
                    }
                }
                // pos == None && !drop should not happen; discard if it does.
            } else {
                // Isoform-coordinate or otherwise unparseable variant (e.g.
                // "A4UGR9-4:757"): can't compare to candidates, but must be
                // kept verbatim so the output is not missing features.
                blocks.push(FtBlock::NonVariant(std::mem::take(&mut pending.lines)));
            }
        } else {
            blocks.push(FtBlock::NonVariant(std::mem::take(&mut pending.lines)));
        }
    }

    *pending = PendingFt::default();
}

/// Extract the text inside `/note="..."` from an accumulated note buffer.
fn extract_note(buf: &str) -> Option<&str> {
    let start = buf.find("/note=\"")?;
    let after_quote = &buf[start + 7..]; // 7 = len of /note="
    Some(match after_quote.find('"') {
        Some(end) => &after_quote[..end],
        None => after_quote,
    })
}

/// Parse a comparable Variant from a note string and position. Mirrors the old
/// `push_entry` logic; the id is left empty because the original text already
/// carries it and Variant equality ignores the id anyway.
fn build_variant(note: &str, pos: (u32, Option<u32>)) -> Variant {
    if let Some(caps) = REGEX_AA_CHANGE.captures(note) {
        Variant {
            pos_start: pos.0,
            pos_end: pos.1,
            aa_ref: caps[1].to_string(),
            aa_new: caps[2].to_string(),
            id: String::new(),
        }
    } else {
        Variant {
            pos_start: pos.0,
            pos_end: pos.1,
            aa_ref: note.to_string(),
            aa_new: String::new(),
            id: String::new(),
        }
    }
}

/// Merge candidate ids into an existing feature's raw text. The ids are appended
/// to the existing `/id="..."` value (pipe-separated, continuous string). If the
/// feature has no `/id=` line, one is added.
fn merge_ids_into_block(lines: &str, ids: &[String]) -> String {
    let addition = ids.join("|");
    let mut out = String::with_capacity(lines.len() + addition.len() + 24);
    let mut inserted = false;

    for line in lines.split_inclusive('\n') {
        if !inserted {
            if let Some(id_pos) = line.find("/id=\"") {
                let after = id_pos + 5; // len of /id="
                if let Some(rel_end) = line[after..].find('"') {
                    let end = after + rel_end;
                    out.push_str(&line[..end]);
                    if end != after {
                        // existing id value is non-empty -> separate with a pipe
                        out.push('|');
                    }
                    out.push_str(&addition);
                    out.push_str(&line[end..]);
                    inserted = true;
                    continue;
                }
            }
        }
        out.push_str(line);
    }

    if !inserted {
        if !out.ends_with('\n') {
            out.push('\n');
        }
        out.push_str(&format!("FT                   /id=\"{}\"\n", addition));
    }

    out
}

fn collect_enst_variants<'a>(
    cur_line: &str,
    regex_enst: &Regex,
    features: &'a HashMap<u32, Vec<Variant>>,
    global_features_inserted: &mut HashSet<u32>,
    insert_candidates: &mut Vec<&'a [Variant]>,
) -> Result<()> {
    for enst in cur_line
        .split_whitespace()
        .map(|t| t.trim_end_matches(';'))
        .filter(|t| regex_enst.is_match(t))
    {
        let key = enst
            .strip_prefix("ENST")
            .context("missing ENST prefix")?
            .split('.')
            .next()
            .context("missing ENST numeric part")?
            .parse::<u32>()
            .context("invalid ENST numeric part")?;

        if let Some(entry) = features.get(&key) {
            insert_candidates.push(entry.as_slice());

            global_features_inserted.insert(key);
        }
    }

    Ok(())
}