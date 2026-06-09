use anyhow::Result;

use crate::aa_masses::AminoMasses;
use crate::variant::Variant;

pub fn build_dummy_entry<W: std::io::Write>(
    enst: &str,
    var: Vec<Variant>,
    seq: &str,
    masses_dict: &AminoMasses,
    mut writer: W,
) -> Result<()> {
    let mut ft_lines = String::new();

    let last_aa = seq.trim_end().chars().last();

    let seqs = format_seq(seq, masses_dict)?;

    if let Some(last_aa_valid) = last_aa {
        for mut variant in var {
            variant.normalize(seq.len() as u32, last_aa_valid);
            ft_lines.push_str(&variant.to_uniprot(seq.len() as u32)?);
        }
    }

    let entry = format!(
        "ID   {}              Unreviewed;          {} AA.\n\
        AC   {};\n\
        DT   1-JAN-1970, integrated into DUMMY.\n\
        DT   1-JAN-1970, sequence version 1.\n\
        DT   1-JAN-1970, entry version 1.\n\
        DE   RecName: Full={}.\n\
        OS   dummy organism (DUMMY).\n\
        OC   Unclassified; dummy (DUMMY).\n\
        OX   NCBI_TaxID=0; dummy (DUMMY).\n\
        RN   [1]\n\
        RP   SEQUENCE.\n\
        RG   dummy group;\n\
        RA   dummy authors;\n\
        RL   Unpublished.\n\
        PE   4: Predicted;\n\
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
