use anyhow::{Context, Result};
use std::collections::HashMap;

pub struct AminoMasses {
    masses: HashMap<char, f64>,
}

impl AminoMasses {
    pub fn new() -> Self {
        let masses = HashMap::from([
            ('G', 57.0519),
            ('A', 71.0788),
            ('S', 87.0782),
            ('P', 97.1167),
            ('V', 99.1326),
            ('T', 101.1051),
            ('C', 103.1388),
            ('L', 113.1594),
            ('I', 113.1594),
            ('N', 114.1038),
            ('D', 115.0886),
            ('Q', 128.1307),
            ('K', 128.1741),
            ('E', 129.1155),
            ('M', 131.1926),
            ('H', 137.1411),
            ('F', 147.1766),
            ('U', 150.0388),
            ('R', 156.1875),
            ('Y', 163.1760),
            ('W', 186.2132),
            ('O', 237.3018),
            ('J', 113.1594),
            ('X', 0.0),
            ('Z', 128.6231),
            ('B', 114.5962),
        ]);

        Self { masses }
    }

    pub fn get(&self, aa: char) -> Result<f64> {
        self.masses
            .get(&aa.to_ascii_uppercase())
            .copied()
            .context(format!("unknown amino acid: {}", aa))
    }

    pub fn sequence_mass(&self, seq: &str) -> Result<f64> {
        let mut total = 0.0;

        for aa in seq.chars() {
            total += self.get(aa)?;
        }

        total += 18.01524;

        Ok(total)
    }
}
