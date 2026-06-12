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
    let mut features_in_entry: HashSet<Variant> = HashSet::new();

    let mut seq = String::new();
    let mut seq_len: u32 = 0;
    let mut collect_seq = false;

    let mut aa_change_pos: Option<(u32, Option<u32>)> = None;
    let mut aa_change_line = String::new();

    let mut other_ft_lines = String::new();

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
        || cur_line.starts_with("DR")
        || cur_line.starts_with("PE")
        || cur_line.starts_with("KW")
        {
            writeln!(writer, "{}", cur_line)?;
        }

        if collect_seq && !cur_line.starts_with("//") {
            seq.push_str(&format!("{}\n", &cur_line));
        }

        if cur_line.starts_with("SQ") {
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
            let last_aa = seq.trim_end().chars().last();

            if let Some(last_aa_valid) = last_aa {
                for candidate_slice in insert_candidates.drain(..) {
                    for mut candidate in candidate_slice.iter().cloned() {
                        candidate.normalize(seq_len, last_aa_valid);
                        let match_in_features = features_in_entry.take(&candidate);
                        if let Some(feature) = match_in_features {
                            if !feature.id.is_empty() {
                                candidate.id.push_str(&format!("|{}", feature.id));
                                write!(writer, "{}", candidate.to_uniprot(seq_len)?)?;
                            }
                        } else {
                            write!(writer, "{}", candidate.to_uniprot(seq_len)?)?
                        }
                    }
                }
            }

            if !confirmed_only {
                write!(writer, "{}", other_ft_lines)?;
                other_ft_lines.clear();
                for feature in features_in_entry.drain() {
                    write!(writer, "{}", feature.to_uniprot(seq_len)?)?;
                }
            } else {
                features_in_entry.clear();
            }

            //TODO add failure logging
            writeln!(writer, "{}{}", seq, cur_line)?;

            features_in_entry.clear();

            seq.clear();

            collect_seq = false;

            continue;
        }

        if cur_line.starts_with("DR   Ensembl;") {
            collect_enst_variants(
                &cur_line,
                &REGEX_ENST,
                features,
                &mut global_features_inserted,
                &mut insert_candidates,
            )
            .context(format!("Failed to parse DR cur_line {}", line_num + 1))?;

            continue;
        }

        if cur_line.starts_with("FT") {
            let content = cur_line.get(21..).unwrap_or("").trim_start();

            let ft_field = cur_line
                .get(2..21)
                .context(format!("Failed to read FT cur_line {}", line_num + 1))?
                .trim();

            if !ft_field.is_empty() {
                if let Some(pos) = aa_change_pos {
                    push_entry(&mut aa_change_line, &pos, &mut features_in_entry, "");

                    aa_change_pos = None;
                    aa_change_line.clear();
                }

                if matches!(ft_field, "VARIANT" | "VAR_SEQ") {
                    if content.starts_with("<") {
                        anyhow::bail!(format!("Failed to read FT cur_line {}", line_num + 1));
                    } else {
                        let mut parts = content.split("..");

                        let pos_start = parts.next().unwrap();
                        //.context(format!("Failed to parse FT cur_line {}", line_num + 1))?;

                        let pos_end = parts.next();

                        let Ok(pos_start_u32) = pos_start.parse::<u32>() else {
                            //TODO maybe add logging. this currently skips crosslinks silently
                            continue;
                        };

                        let pos_end_u32 = match pos_end.map(|p| p.parse::<u32>()) {
                            Some(Ok(val)) => Some(val),
                            Some(Err(_)) => {
                                continue;
                            }
                            None => None,
                        };

                        aa_change_pos = Some((pos_start_u32, pos_end_u32));

                        continue;
                    }
                } else {
                    // Non-variant feature: buffer the header line
                    if !confirmed_only {
                        other_ft_lines.push_str(&format!("{}\n", &cur_line));
                    }
                    continue;
                }
            }

            if let Some(pos) = aa_change_pos {
                //collect note cur_line
                if content.starts_with("/note=") {
                    aa_change_line = cur_line;
                    continue;
                }
                //entry finished with id
                else if content.starts_with("/id=") {
                    let feature_id = content
                        .strip_prefix("/id=\"")
                        .and_then(|s| s.strip_suffix('"'))
                        .unwrap_or("");

                    push_entry(
                        &mut aa_change_line,
                        &pos,
                        &mut features_in_entry,
                        feature_id,
                    );

                    aa_change_pos = None;
                    aa_change_line.clear();
                    continue;
                } else if !content.starts_with("/") {
                    aa_change_line.push_str(cur_line.get(21..).unwrap());
                }
            } else if !confirmed_only && !other_ft_lines.is_empty() {
                // continuation of a non-variant feature
                other_ft_lines.push_str(&format!("{}\n", &cur_line));
                continue;
            }
        }
    }
    Ok(global_features_inserted)
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

fn push_entry(
    aa_change_line: &mut str,
    change_pos: &(u32, Option<u32>),
    features_in_entry: &mut HashSet<Variant>,
    id: &str,
) {
    let note = aa_change_line.get(28..).unwrap();
    if let Some(caps) = REGEX_AA_CHANGE.captures(note) {
        let aa_ref = &caps[1];
        let aa_new = &caps[2];

        features_in_entry.insert(Variant {
            pos_start: change_pos.0,
            pos_end: change_pos.1,
            aa_ref: aa_ref.to_string(),
            aa_new: aa_new.to_string(),
            id: id.to_string(),
        });
    } else {
        features_in_entry.insert(Variant {
            pos_start: change_pos.0,
            pos_end: change_pos.1,
            aa_ref: note
                .to_string()
                .strip_suffix('"')
                .unwrap_or(note)
                .to_string(),
            aa_new: String::new(),
            id: id.to_string(),
        });
    }
}

/*
// New feature begins before /id=
if !content.starts_with('/') && !aa_change_line.is_empty() {
    push_entry(
        &mut aa_change_line,
        &pos,
        &mut variants_in_entry,
        &String::new(),
    );

    aa_change_pos = None;
    aa_change_line.clear();
    // fall through and parse the new feature
} else {
    aa_change_line.push_str(content);
    continue;
}
*/
/*
    if !insert_candidates.is_empty() {
    //CASE: FT VARIANT or VAR_SEQ cur_line with position found in previous iteration
    if let Some(pos) = aa_change_pos {
        let segment = cur_line
            .get(21..)
            .and_then(|s| s.split_whitespace().next())
            .context(format!("Failed to read FT cur_line {}", line_num + 1))?;

        //CASE: in FT cur_line after VARIANT or VAR_SEQ
        if segment.starts_with("/note=") {
            aa_change_line = cur_line;
            continue;
        }
        else if segment.starts_with("/id=") {
            let variant_id = segment
                .strip_prefix("/id=\"")
                .and_then(|s| s.strip_suffix('"'))
                .unwrap()
                .to_string();

            //CASE: variant entry complete -> insert complete entry
            push_entry(&mut aa_change_line, &pos, &mut variants_in_entry, &variant_id);
            aa_change_pos = None;
            aa_change_line.clear();
            continue;
        }
        else if segment.starts_with("/") && aa_change_line > 1{
            //CASE: variant entry complete -> insert complete entry
            push_entry(&mut aa_change_line, &pos, &mut variants_in_entry, &variant_id);
            aa_change_pos = None;
            aa_change_line.clear();
            continue;
        }
    }
}

     */
