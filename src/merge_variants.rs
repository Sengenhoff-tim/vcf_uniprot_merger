use regex::Regex;
use std::collections::{HashMap, HashSet};
use std::io::{BufRead, Write};

use crate::variant::Variant;

use anyhow::{Context, Result};

pub fn add_variants<R, W>(
    reader: R,
    writer: &mut W,
    variants: &HashMap<u32, Vec<Variant>>,
) -> Result<HashSet<u32>>
where
    R: BufRead,
    W: Write,
{
    let regex_enst =
        Regex::new(r"^ENST\d{11}(?:\.\d+)?$").context("Failed to compile ENST regex")?;
    let regex_aa_change =
        Regex::new(r"^([A-Z]+) -> ([A-Z]+)").context("Failed to compile AA change regex")?;

    let mut insert_candidates: Vec<&[Variant]> = Vec::new();
    let mut global_variants_inserted: HashSet<u32> = HashSet::new();
    let mut variants_in_entry: HashSet<Variant> = HashSet::new();

    let mut write_seq = false;

    let mut aa_change_pos: Option<u32> = None;
    let mut aa_change_line: Option<String> = None;

    for (line_num, line_result) in reader.lines().enumerate() {
        let line = line_result.context(format!("Failed to read line {}", line_num + 1))?;

        if line.starts_with("ID")
            || line.starts_with("AC")
            || line.starts_with("CC")
            || line.starts_with("DT")
            || line.starts_with("DE")
            || line.starts_with("OS")
            || line.starts_with("OC")
            || line.starts_with("OX")
            || line.starts_with("RN")
            || line.starts_with("RP")
            || line.starts_with("RG")
            || line.starts_with("RA")
            || line.starts_with("RL")
            || line.starts_with("PE")
            || (write_seq && !line.starts_with("//"))
        {
            writeln!(writer, "{}", line)?;
            continue;
        }

        if line.starts_with("SQ") {
            for candidate_slice in insert_candidates.drain(..) {
                for candidate in candidate_slice {
                    if !variants_in_entry.contains(candidate) {
                        write!(writer, "{}", candidate.to_str_unchecked())?;
                    }
                }
            }

            writeln!(writer, "{}", line)?;

            variants_in_entry.clear();

            write_seq = true;

            continue;
        }

        if line.starts_with("//") {
            writeln!(writer, "{}", line)?;

            write_seq = false;

            continue;
        }

        if line.starts_with("DR   Ensembl;") {
            collect_enst_variants(
                &line,
                &regex_enst,
                variants,
                &mut global_variants_inserted,
                &mut insert_candidates,
            )
            .context(format!("Failed to parse DR line {}", line_num + 1))?;

            continue;
        }

        if line.starts_with("FT") {
            writeln!(writer, "{}", line)?;

            if !insert_candidates.is_empty() {
                if let Some(pos) = aa_change_pos {
                    let segment = line
                        .get(21..)
                        .and_then(|s| s.split_whitespace().next())
                        .context(format!("Failed to read FT line {}", line_num + 1))?;

                    if let Some(ref mut aa_change) = aa_change_line {
                        if !segment.starts_with("/") {
                            aa_change.push_str(line.get(21..).unwrap());
                        } else {
                            if let Some(caps) =
                                regex_aa_change.captures(aa_change.get(28..).unwrap())
                            {
                                let g1_bytes = caps[1].as_bytes();
                                let g2_bytes = caps[2].as_bytes();

                                variants_in_entry.insert(Variant {
                                    pos,
                                    aa_ref: g1_bytes.to_vec(),
                                    aa_new: g2_bytes.to_vec(),
                                });
                            }
                            aa_change_pos = None;
                            aa_change_line = None;
                        }
                    } else if segment.starts_with("/note=") {
                        aa_change_line = Some(line);
                    }
                } else {
                    let mut parts = line
                        .get(2..)
                        .map(|s| s.split_whitespace())
                        .context(format!("Failed to read FT line {}", line_num + 1))?;
                    let segment = parts
                        .next()
                        .context(format!("Failed to read FT line {}", line_num + 1))?;

                    if segment.starts_with("VARIANT") || segment.starts_with("VAR_SEQ") {
                        let position = parts
                            .next()
                            .context(format!("Failed to read FT line {}", line_num + 1))?;

                        if position.starts_with("<") {
                            anyhow::bail!(format!("Failed to read FT line {}", line_num + 1));
                        } else {
                            let pos = position
                                .split(".")
                                .next()
                                .context(format!("Failed to parse FT line {}", line_num + 1))?;

                            let Ok(pos_u32) = pos.parse::<u32>() else {
                                //TODO maybe add logging
                                continue;
                            };

                            aa_change_pos = Some(pos_u32);
                        }
                    }
                }
            }
            continue;
        }
    }

    Ok(global_variants_inserted)
}

fn collect_enst_variants<'a>(
    line: &str,
    regex_enst: &Regex,
    variants: &'a HashMap<u32, Vec<Variant>>,
    global_variants_inserted: &mut HashSet<u32>,
    insert_candidates: &mut Vec<&'a [Variant]>,
) -> Result<()> {
    for enst in line
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

        if let Some(entry) = variants.get(&key) {
            insert_candidates.push(entry.as_slice());

            global_variants_inserted.insert(key);
        }
    }

    Ok(())
}
