use anyhow::{Context, Result};
use std::{
    collections::HashMap,
    io::{BufRead, BufReader},
};

use regex::Regex;

use crate::variant::Variant;

pub fn build_variant_dict<R: std::io::Read>(
    reader: BufReader<R>,
) -> Result<HashMap<u32, Vec<Variant>>> {
    let re_id = Regex::new(r"ENST\d{11}").unwrap();
    let re_aa = Regex::new(r"(\d+)([A-Z*]+)>(\d+)([A-Z*]+)").unwrap();

    let mut variants: HashMap<u32, Vec<Variant>> = HashMap::new();

    for (line_idx, line) in reader.lines().enumerate() {
        let line = line.context(format!("Failed to read line {}", line_idx))?;

        let (id, aa_change) = line
            .split_once(' ')
            .context(format!("Missing ' ' delimiter in line {}", line_idx))?;

        //TODO handle malformed input
        if aa_change != "."
            && re_id.is_match(id)
            && let Some(caps) = re_aa.captures(aa_change)
            && let Ok(key) = id[4..].parse::<u32>()
        {
            let pos1 = caps[1].parse::<u32>()?;
            let _pos2 = caps[3].parse::<u32>()?;

            /*
            if pos1 != pos2 {
                bail!("pos old {} != pos new {} in line {}", pos1, pos2, line_idx);
            }
             */

            let aa_ref = &caps[2];
            let aa_new = &caps[4];

            let mut variant = Variant::from_match(pos1, aa_ref, aa_new);
            let variants_for_id = variants.entry(key).or_default();

            if !variants_for_id.contains(&variant) {
                variant.id = format!("{}:{}", variants_for_id.len(), id);
                variants_for_id.push(variant);
            }
        }
    }

    Ok(variants)
}
