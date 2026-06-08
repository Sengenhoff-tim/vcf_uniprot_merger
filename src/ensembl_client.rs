use anyhow::{Context, Result, anyhow};
use std::collections::HashMap;
use std::io::Write;
use std::time::Duration;
use tokio::time::sleep;

use futures::{StreamExt, stream};
use reqwest::Client;
use serde::{Deserialize, Serialize};

use crate::{aa_masses::AminoMasses, parse_uniprot::build_dummy_entry, variant::Variant};

const MAX_POST_SIZE: usize = 50;
const CONCURRENCY: usize = 10;
const URL: &str = "https://rest.ensembl.org/sequence/id?type=protein";

#[derive(Serialize)]
struct RequestBody<'a> {
    ids: &'a [String],
}

#[derive(Debug, Deserialize)]
pub struct SequenceResponse {
    pub molecule: String,
    pub seq: String,
    pub query: String,
    pub id: String,
}

pub async fn fetch_sequences<W: Write>(
    mut variants: HashMap<u32, Vec<Variant>>,
    mut writer: W,
) -> Result<()> {
    let ids: Vec<String> = variants.keys().map(|&n| format!("ENST{:011}", n)).collect();

    let requests = ids.chunks(MAX_POST_SIZE);

    let masses_dict = AminoMasses::new();

    let client = Client::new();

    let mut responses = stream::iter(requests)
        .map({
            move |chunk| {
                let client = client.clone();

                async move {
                    let mut attempt = 0;
                    let max_attempts = 3;

                    loop {
                        attempt += 1;

                        let result = async {
                            let resp = client
                                .post(URL)
                                .json(&RequestBody { ids: chunk })
                                .send()
                                .await?
                                .error_for_status()?
                                .json::<Vec<SequenceResponse>>()
                                .await?;

                            Ok::<_, reqwest::Error>(resp)
                        }
                        .await;

                        match result {
                            Ok(resp) => return Ok(resp),

                            Err(_) if attempt < max_attempts => {
                                let backoff = Duration::from_millis(100 * 2u64.pow(attempt - 1));
                                sleep(backoff).await;
                                continue;
                            }

                            Err(err) => return Err(err),
                        }
                    }
                }
            }
        })
        .buffer_unordered(CONCURRENCY);

    while let Some(result) = responses.next().await {
        let batch = result.context("request failed")?;

        for record in batch {
            if record.molecule != "protein" {
                return Err(anyhow!(
                    "Expected molecule type 'protein', got '{}'",
                    record.molecule
                ));
            }
            let key = record
                .query
                .strip_prefix("ENST")
                .context("missing ENST prefix")?
                .parse::<u32>()
                .context("invalid ENST numeric part")?;

            let var = variants.remove(&key).context("missing variant for key")?;

            build_dummy_entry(&record.id, var, &record.seq, &masses_dict, &mut writer)?;
        }
    }

    Ok(())
}
