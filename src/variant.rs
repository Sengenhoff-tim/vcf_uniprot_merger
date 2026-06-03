use anyhow::{Result, bail};

#[derive(Eq, PartialEq, Hash)]
pub struct Variant {
    pub pos_start: u32,
    pub pos_end: Option<u32>,
    pub aa_ref: String,
    pub aa_new: String,
}

impl Variant {
    pub fn from_match(pos: u32, aa_ref: &str, aa_new: &str) -> Variant {

        let pos_end = if aa_ref.len() > 1 {
            Some(pos + aa_ref.len() as u32)
        } else {
            None
        };
        
        Variant {
            pos_start: pos,
            pos_end: pos_end,
            aa_ref: aa_ref.to_string(),
            aa_new: aa_new.to_string(),
        }
    }

    pub fn to_uniprot(&self, seq_len: u32) -> Result<String> {

    if let Some(stop_pos_ref) = self.aa_ref.find("*") {
        if self.pos_start == seq_len + 1{
            return Ok(
                format!(
                    "FT   VAR_SEQ         {}\n\
                    FT                   /note=\"{} -> {}\"\n",
                    self.pos_start,
                    self.aa_ref,
                    self.aa_new,
                )
            )
        }
        else if stop_pos_ref > seq_len as usize + 1{
            bail!(format!("Suspicious sequence length {} for reference stop position {}; Old seq: {}, new seq: {}", seq_len, stop_pos_ref, self.aa_ref, self.aa_new))
        }
    }

    if let Some(stop_pos_new) = self.aa_new.find("*") {
        if stop_pos_new as u32 == self.pos_start {
            return Ok(
                format!(
                    "FT   VAR_SEQ         {}\n\
                    FT                   /note=\"Missing in sample\"\n",
                    self.pos_start
                )
            )

        }
        else {
            return Ok(
                format!(
                    "{}FT   VAR_SEQ         {}..{}\n\
                    FT                   /note=\"Missing in sample\"\n",
                    format_variant(self.pos_start, self.pos_end.map(|x| x -1), &self.aa_ref, &self.aa_new),
                    self.pos_end.unwrap(),
                    seq_len
                )
            )
        }
    }

    Ok(format_variant(self.pos_start, self.pos_end, &self.aa_ref, &self.aa_new))
}


}


fn format_variant(
    pos_start: u32,
    pos_end: Option<u32>,
    aa_ref: &str,
    aa_new: &str,
) -> String {
    let pos = pos_end
        .map(|end| format!("{pos_start}..{end}"))
        .unwrap_or_else(|| pos_start.to_string());

    let note = format!("{aa_ref} -> {aa_new}");

    const LINE_LIMIT: usize = 80;
    const PREFIX: &str = "FT                   /note=\"";
    const CONTINUATION: &str = "FT                   ";

    let mut result = format!("FT   VARIANT         {pos}\n");

    let mut remaining = note.as_str();
    let mut first_line = true;

    while !remaining.is_empty() {
        let prefix = if first_line {
            PREFIX
        } else {
            CONTINUATION
        };

        let available = if remaining.len() + prefix.len() + 1 <= LINE_LIMIT {
            LINE_LIMIT - prefix.len() - 1 // reserve room for closing quote
        } else {
            LINE_LIMIT - prefix.len()
        };

        let split_at = if remaining.len() > available {
            remaining[..available]
                .rfind(' ')
                .unwrap_or(available)
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

