use std::{fs::File, path::PathBuf};
use structopt::StructOpt;
use std::io::Write;

mod kmer;
use crate::kmer::count_read_base_kmers;

#[derive(StructOpt, Debug)]
#[structopt(name = "basic")]
struct Opt {
    #[structopt(short, long, default_value = "6")]
    k: u8,
    #[structopt(short, long, parse(from_os_str), default_value="SRR12458125_1.fastq.gz")]
    input: PathBuf,
    #[structopt(short, long, parse(from_os_str), default_value="debug_out.tsv")]
    debugout: PathBuf,
    //#[structopt(short, long, parse(from_os_str))]
    //output: PathBuf,
    #[structopt(short, long, default_value = "0.01")]
    error_rate: f64,
}

fn main()  -> std::io::Result<()> {
    let opt = Opt::from_args();
    let kmer_counts = count_read_base_kmers(opt.k, opt.input);

    let mut file = File::create(&opt.debugout)?;
    for i in 0..kmer_counts.len() {
        let counts = &kmer_counts[i].count;
        for kmer_offset in 0..counts.len() {
            let kmer_seq = needletail::bitkmer::bitmer_to_bytes(crate::kmer::MutableKmer::kmer_low_bits_to_native_encoding(kmer_offset as u64, opt.k));
            write!(file, "{}\t{}\t{}\t{}\t{:?}\n", i, kmer_offset, std::str::from_utf8(&kmer_seq).unwrap(), counts[kmer_offset], crate::kmer::calculate_state(kmer_offset, opt.k, &kmer_counts[i], opt.error_rate))?;
        }
    }
    file.flush()?;
    Ok(())
}
