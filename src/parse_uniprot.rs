use anyhow::Result;

use crate::variant::Variant;
use crate::aa_masses::AminoMasses;

pub fn build_dummy_entry<W: std::io::Write>(
    enst: &str,
    var: Vec<Variant>,
    seq: &str,
    masses_dict: &AminoMasses,
    mut writer: W,
) -> Result<()> {
    
    let mut ft_lines = String::new();
    for variant in var {
        ft_lines.push_str(&variant.to_str_unchecked());
    }
    
    let seqs = format_seq(seq, masses_dict)?;
    
    let entry = format!(
        "ID   {}              Unreviewed;          {} AA.\n\
        AC   {};\n\
        DE   RecName: Full={}.\n\
        {}\
        {}\
        //\n",
        enst,
        seq.len(),
        enst,
        enst,
        ft_lines,
        seqs
    );
    
    writer.write_all(entry.as_bytes())?;
    Ok(())
}

    pub fn format_seq(seq: &str, masses_dict: &AminoMasses) -> Result<String> {

        let mass = masses_dict.sequence_mass(seq)?.round() as u64;

        Ok(format!(
            "SQ   SEQUENCE {} AA; {} MW; {} CRC64;\n{}",
            seq.len(),
            mass,
            "0000000000000000",
            format_sequence(seq)?
        ))
    }

    fn format_sequence(seq: &str) -> Result<String> {
        let mut result = String::new();
        let chars: Vec<char> = seq.chars().collect();
        
        for chunk in chars.chunks(60) {
            result.push_str("     ");
            
            for (j, group) in chunk.chunks(10).enumerate() {
                if j > 0 {
                    result.push(' ');
                }
                result.push_str(&group.iter().collect::<String>());
            }
            result.push('\n');
        }
        
        Ok(result)
    }
