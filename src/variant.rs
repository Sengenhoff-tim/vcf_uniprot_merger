use anyhow::Result;

#[derive(Eq, PartialEq, Hash, Clone)]
pub struct Variant {
    pub pos_start: u32,
    pub pos_end: Option<u32>,
    pub aa_ref: String,
    pub aa_new: String,
}

impl Variant {
    pub fn from_match(pos: u32, aa_ref: &str, aa_new: &str) -> Variant {
        let pos_end = if aa_ref.len() > 1 {
            Some(pos + aa_ref.len() as u32 - 1)
        } else {
            None
        };

        Variant {
            pos_start: pos,
            pos_end,
            aa_ref: aa_ref.to_string(),
            aa_new: aa_new.to_string(),
        }
    }

    pub fn normalize(&mut self, seq_len: u32, last_seq_char: char) {
        if self.aa_ref.ends_with("*") {
            //CASE *->AA*
            let combined = format!("{}{}", last_seq_char, self.aa_new.trim_end_matches('*'));
            self.pos_start = seq_len;
            self.pos_end = Some(seq_len);
            self.aa_ref = last_seq_char.to_string();
            self.aa_new = combined;
        }

        if self.aa_new.find("*").is_some() {
            self.pos_end = Some(seq_len);
            //CASE A->* & AA->*
            if self.aa_new.len() == 1 {
                self.aa_new = String::new();
                self.aa_ref = String::new();

                return;
            }
            //CASE A->A*
            //CASE AA->A*
            self.aa_new = self.aa_new.trim_end_matches('*').to_string();
        }
        //CASE A->A
        //CASE A->AA
        //CASE AA->A
        //CASE AA->AA
    }

    pub fn to_uniprot(&self, seq_len: u32) -> Result<String> {
        if self.aa_ref.is_empty() && self.aa_new.is_empty() {
            return Ok(format!(
                "FT   VAR_SEQ         {}..{}\n\
                FT                   /note=\"Missing in sample\"\n",
                self.pos_start, seq_len
            ));
        }

        Ok(format_variant(
            self.pos_start,
            self.pos_end,
            &self.aa_ref,
            &self.aa_new,
        ))
    }
}

fn format_variant(pos_start: u32, pos_end: Option<u32>, aa_ref: &str, aa_new: &str) -> String {
    let mut result = pos_end
        .filter(|end| *end != pos_start)
        .map(|end| format!("FT   VAR_SEQ         {}..{}\n", pos_start, end))
        .unwrap_or_else(|| format!("FT   VARIANT         {}\n", pos_start));

    let note = format!("{aa_ref} -> {aa_new}");

    const LINE_LIMIT: usize = 80;
    const PREFIX: &str = "FT                   /note=\"";
    const CONTINUATION: &str = "FT                   ";

    let mut remaining = note.as_str();
    let mut first_line = true;

    while !remaining.is_empty() {
        let prefix = if first_line { PREFIX } else { CONTINUATION };

        let available = if remaining.len() + prefix.len() < LINE_LIMIT {
            LINE_LIMIT - prefix.len() - 1 // reserve room for closing quote
        } else {
            LINE_LIMIT - prefix.len()
        };

        let split_at = if remaining.len() > available {
            remaining[..available].rfind(' ').unwrap_or(available)
        } else {
            remaining.len()
        };

        let chunk = &remaining[..split_at];
        remaining = remaining[split_at..].trim_start();

        result.push_str(prefix);
        result.push_str(chunk);

        if remaining.is_empty() {
            result.push('"');
        }

        result.push('\n');
        first_line = false;
    }

    result
}
