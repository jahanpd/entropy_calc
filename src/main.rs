// use anndata;
use clap::arg;
use database::EntropyData;
use indicatif::{MultiProgress, ProgressBar, ProgressStyle};
use std::sync::{Mutex, mpsc};
use std::collections::HashMap;
use noodles::{bam, bam::bai, sam::{Header, record::cigar::op::kind::Kind}};
use noodles::sam::record::sequence::Base;
use noodles::sam::record::data::field::tag::Tag;
use rayon::prelude::*;
mod import;
mod database;
// goal here is to create binary that achieves the following
// takes 3 arguments
//      1. path to sqlite db of gene annotations indexed by gene id eg WBGene00000000 (from gffutils)
//      2. path to newline delimited file of gene ids for which to calculate entropy
//      3. output sqlite database location - will be overwritten


#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {

    rayon::ThreadPoolBuilder::new()
        .num_threads(12)
        .build_global().unwrap();


    let matches = clap::Command::new("RNA/Genomic Entropy Calculation")
        .version("0.0.1")
        .author("Jahizzle PDizzle")
        .about("Calculating the entropy of reads in RNAseq datasets")
        .args([
            arg!(--genes <PATH> "location of the csv with gene names, contig, start, stop eg WBGene0000000, I, 120, 500 \n...")
                .default_value("/mnt/RNASeq/data/genes.csv"),
            arg!(--bams <PATH> "folder where bam and bai files are stored")
                .default_value("/mnt/RNASeq/data/scRNAseq/"),
            arg!(--output <PATH> "location to write output sqlite database")
                .default_value("/mnt/RNASeq/data/entropy.db"),
        ])
        .get_matches();
    let genes_file = matches.get_one::<String>("genes").unwrap();
    let bam_folder = matches.get_one::<String>("bams").unwrap();
    let out_db = matches.get_one::<String>("output").unwrap();
    // TODO path validation
    let bam_paths = std::fs::read_dir(bam_folder).unwrap()
        .map(|x| x.unwrap().path().as_path().to_str().unwrap().to_string())
        .filter(|x| !(x.contains("bai") || x.contains("h5ad")))
        .filter(|x| x.contains("d8_1.bam.1"))
        .collect::<Vec<String>>();
    dbg!(&bam_paths);

    // import csv
    let csv: String = std::fs::read_to_string(genes_file)?.parse()?;
    let gene_store = import::StoreGenes::new(csv);

    let m = MultiProgress::new();
    let sty = ProgressStyle::with_template(
        "[{elapsed_precise} | {eta_precise}] {bar:40.cyan/blue} {pos:>9}/{len:9} {msg}",
    )
    .unwrap()
    .progress_chars("##-");

    // initialise message passing channel to communicate between tasks
    let (tx, rx) = mpsc::sync_channel::<Option<database::EntropyData>>(10000);

    let edb = database::EntropyDB::new(out_db, 8, 60, rx).await.unwrap();
    let handle = std::thread::spawn(move || {
        edb.listen();
    });

    bam_paths.par_iter().for_each_with(tx.clone(), |tx, path| {
        let pbout = m.add(ProgressBar::new(gene_store.gene_id.len() as u64));
        pbout.set_style(sty.clone());

        // let tc = tx.clone();
        let _ = gene_store.genes().par_bridge().for_each_with(tx.clone(), |tx, g| {

            let region = format!("{}:{}-{}", g.contig, g.start, g.end).parse().unwrap();
            let mut reader = std::fs::File::open(path).map(bam::Reader::new).unwrap();
            let header: Header = reader.read_header().unwrap();
            let index = bai::read(format!("{}{}", path, ".bai")).unwrap();
            let q = reader.query(&header, &index, &region).unwrap();

            let store: Mutex<HashMap<String, HashMap<usize, Vec<f64>>>> = Mutex::new(HashMap::new());

            q.enumerate().par_bridge().for_each(|(idx, r)| {
                let record = r.unwrap();
                let seq = record.sequence();
                let astart = record.alignment_start().unwrap().get();
                let gstart = g.start as usize;
                let gend = g.end as usize;
                let mut idx = 1; // for accessing store
                let mut acurr = astart - gstart + 1; // for tracking alignment
                let maxread = gend - gstart + 1;
                let cb = record.data().get(Tag::CellBarcodeSequence).unwrap();

                let mut innerstore: HashMap<usize, usize> = HashMap::new();

                for op in record.cigar().iter() {

                    if op.kind() == Kind::Match {

                        (idx.clone()..idx.clone()+op.len()).for_each(|pos| {

                            let i = noodles::core::Position::new(pos).unwrap();
                            let base = seq.get(i).unwrap_or(&Base::L);
                            if base == &Base::A && acurr >= 0 as usize && pos <= maxread {
                                let _ = innerstore.insert(pos, 0);
                            } else if base == &Base::C && acurr >= 0 as usize && pos <= maxread {
                                let _ = innerstore.insert(pos, 1);
                            } else if base == &Base::G && acurr >= 0 as usize && pos <= maxread  {
                                let _ = innerstore.insert(pos, 2);
                            } else if base == &Base::T && acurr >= 0 as usize && pos <= maxread {
                                let _ = innerstore.insert(pos, 3);
                            } else {
                                // dbg!(base);
                            }
                            if acurr >= 0 as usize {
                                idx += 1;
                            }
                            acurr += 1;
                        });
                    } else if op.kind() == Kind::SequenceMismatch {
                        dbg!(&op);
                        if acurr >= 0 as usize {
                            idx += op.len();
                        }
                        acurr += op.len();
                    } else {
                        if acurr >= 0 as usize {
                            idx += op.len();
                        }
                        acurr += op.len();
                    }
                }
                let mut storelock = store.lock().unwrap();
                if !storelock.contains_key(&cb.to_string()) {
                    storelock.insert(cb.to_string(), HashMap::new());
                }
                innerstore.iter().for_each(|(k, v)| {
                    let ss = storelock.get_mut(&cb.to_string()).unwrap();
                    if !ss.contains_key(k) {
                        let mut ve = vec![0f64, 0f64, 0f64, 0f64, 0f64, 0f64];
                        ve[*v] += 1.;
                        let _ = ss.insert(*k, ve);
                    } else {
                        let ve = ss.get_mut(k).unwrap();
                        ve[*v] += 1.;
                    }

                });

            });

            let storelock = store.lock().unwrap();
            let entropy_celltype = storelock.par_iter().map(|(celltype, bases)| {
                // let mut seq = vec![];
                let e: f64 = bases.iter().map(|(k,v)| {
                    let mut total: f64 = v.iter().sum();
                    total += 1e-8;
                    // seq.push(v.clone().iter().map(|x| *x as i64).collect::<Vec<i64>>());
                    v.iter().map(|p| {
                        -((p/total) + 1e-8).ln() * (p/(total + 1e-8))
                    }).sum::<f64>()
                }).sum::<f64>() / bases.iter().len() as f64;
                (celltype, e)
            });

            // let entropy = entropy_celltype.clone().fold(0f64,|acc, (c, e)| acc + e);

            pbout.set_message(
                format!("{} {}",
                        path.clone().split("/").collect::<Vec<&str>>()[5],
                        g.gene_id.clone(),
                        // entropy,
                )
            );
            pbout.inc(1);

            let time = path.split("_").collect::<Vec<&str>>()[4];

            // TODO send data to loop
            entropy_celltype.for_each_with(tx.clone(), |tx, (celltype, entropy)| {
                let _ = tx.send(Some(EntropyData {
                    ent: entropy,
                    celltype: celltype.to_string(),
                    gene: g.gene_id.clone(),
                    time: time.to_string(),
                    // sequence:seq
                }));
            });

            drop(storelock);
            drop(reader);
            drop(index);

        });

    });
    let _ = tx.send(None);
    let _ = handle.join();

    Ok(())
}
