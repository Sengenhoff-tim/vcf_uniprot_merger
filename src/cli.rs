// cli.rs
use clap::Parser;
use anyhow::Result;

#[derive(Parser, Debug)]
pub struct Config {
    #[arg(short = 'o', long = "output_path", value_name = "PATH")]
    pub output_path: String,

    #[arg(short = 'i', long = "bcftools_input_path", value_name = "PATH")]
    pub bcftools_input_path: String,

    #[arg(short = 'u', long = "uniprot_input_path", value_name = "PATH")]
    pub uniprot_path: String,

    #[arg(short = 'z', long = "zip", value_name = "BOOLEAN")]
    pub zip: bool,
}

pub fn parse_args() -> Result<Config> {
    Ok(Config::parse())
}
