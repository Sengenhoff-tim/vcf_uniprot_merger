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
    let re_aa = Regex::new(r"\d+[A-Z]+>\d+[A-Z]+").unwrap();

    let mut variants: HashMap<u32, Vec<Variant>> = HashMap::new();

    for line in reader.lines() {
        let line = line.context("Failed to read line")?;

        let (id, aa_change) = line.split_once(' ').context("Missing ' '")?;

        //TODO handle malformed input
        if aa_change != "."
            && re_id.is_match(id)
            && re_aa.is_match(aa_change)
            && let Ok(key) = id[4..].parse::<u32>()
        {
            variants
                .entry(key)
                .or_default()
                .push(Variant::from_string_unchecked(aa_change));
        }
    }

    Ok(variants)
}

