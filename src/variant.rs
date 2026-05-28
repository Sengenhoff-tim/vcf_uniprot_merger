#[derive(Eq, PartialEq, Hash)]
pub struct Variant {
    pub pos: u32,
    pub aa_ref: Vec<u8>,
    pub aa_new: Vec<u8>,
}

impl Variant {
    pub fn from_string_unchecked(aa_str: &str) -> Variant {
        let (left, right) = aa_str.split_once('>').unwrap();

        let l = left.bytes().position(|b| b.is_ascii_alphabetic()).unwrap();

        let r = right.bytes().position(|b| b.is_ascii_alphabetic()).unwrap();

        let (pos, aa_ref) = left.split_at(l);

        Variant {
            pos: pos.parse::<u32>().unwrap(),
            aa_ref: aa_ref.as_bytes().to_vec(),
            aa_new: right.as_bytes()[r..].to_vec(),
        }
    }

    pub fn to_str_unchecked(&self) -> String {
        format!(
            "FT   VARIANT         {}\n\
            FT                   /note=\"{} -> {}\"\n",
            self.pos,
            std::str::from_utf8(&self.aa_ref).unwrap(),
            std::str::from_utf8(&self.aa_new).unwrap(),
        )
    }
}
