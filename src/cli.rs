// cli.rs
use anyhow::Result;
use clap::Parser;

#[derive(Parser, Debug)]
pub struct Config {
    #[arg(short = 'o', long = "output_path", value_name = "PATH")]
    pub output_path: String,

    #[arg(short = 'i', long = "bcftools_input_path", value_name = "PATH")]
    pub bcftools_input_path: String,

    #[arg(short = 'u', long = "uniprot_input_path", value_name = "PATH")]
    pub uniprot_path: String,

    #[arg(short = 'f', long = "ensembl_fallback", value_name = "BOOLEAN")]
    pub ensembl_fallback: bool,

    #[arg(short = 'c', long = "confirmed_only", value_name = "BOOLEAN")]
    pub confirmed_only: bool,

    #[arg(short = 'z', long = "zip", value_name = "BOOLEAN")]
    pub zip: bool,
}

pub fn parse_args() -> Result<Config> {
    Ok(Config::parse())
}
