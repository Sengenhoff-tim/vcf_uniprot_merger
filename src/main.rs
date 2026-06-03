use anyhow::Result;
use flate2::{Compression, write::GzEncoder};
use std::fs::File;
use std::io::{BufReader, BufWriter};

mod aa_masses;
mod cli;
mod ensembl_client;
mod merge_variants;
mod parse_uniprot;
mod read_bcftools;
mod variant;
mod writer_wrapper;

use crate::merge_variants::add_variants;
use crate::read_bcftools::build_variant_dict;
use crate::writer_wrapper::WriterWrapper;

use crate::ensembl_client::fetch_sequences;

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    let config = cli::parse_args()?;

    let bcftools_file = File::open(config.bcftools_input_path)?;
    let bcftool_reader = BufReader::new(bcftools_file);

    let mut variants = build_variant_dict(bcftool_reader)?;

    let uniprot_file = File::open(config.uniprot_path)?;
    let uniprot_reader = BufReader::new(uniprot_file);

    let output_file = File::create(config.output_path)?;

    let encoder = if config.zip {
        WriterWrapper::Compressed(GzEncoder::new(output_file, Compression::default()))
    } else {
        WriterWrapper::Uncompressed(output_file)
    };

    let mut output_writer = BufWriter::new(encoder);

    let inserted = add_variants(uniprot_reader, &mut output_writer, &variants, config.confirmed_only)?;

    if config.ensembl_fallback {
        variants.retain(|key, _value| !inserted.contains(key));

        fetch_sequences(variants, &mut output_writer).await?;
    } 

    Ok(())
}
