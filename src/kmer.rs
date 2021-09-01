//! TODO: crate-level documentation
use std::path::PathBuf;

use needletail::{bitkmer::BitKmer, parse_fastx_file};

static NUCLEOTIDE_BASES: &'static [u8] = &['A' as u8, 'C' as u8, 'G' as u8, 'T' as u8];

const NEEDLETAIL_ENCODING_NUCLEOTIDE_BASES: [Option<u8>; 256] = {
    let mut lookup = [None; 256];

    lookup[b'A' as usize] = Some(0);
    lookup[b'C' as usize] = Some(1);
    lookup[b'G' as usize] = Some(2);
    lookup[b'T' as usize] = Some(3);
    lookup[b'a' as usize] = Some(0);
    lookup[b'c' as usize] = Some(1);
    lookup[b'g' as usize] = Some(2);
    lookup[b't' as usize] = Some(3);

    lookup
};

pub trait MutableKmer {
    fn neighbours(&self) -> Vec<Self> where Self: Sized;
    fn next_kmers(&self) -> Vec<Self> where Self: Sized;
    fn set_position(&self, position: usize, base: u8) -> Self;
    // Ensures we can use the kmer itself as an array index
    fn to_kmer_low_bits(&self) -> u64;
    // Converts back to native kmer encoding for the library
    fn kmer_low_bits_to_native_encoding(low_bit_encoded : u64, k : u8) -> Self;
}

impl MutableKmer for BitKmer {
    fn neighbours(&self) -> Vec<BitKmer> where Self: Sized {
        let kmer = self.1;
        let mut result = Vec::new();
        for i in 0..kmer {
            for base in NUCLEOTIDE_BASES {
                let new_kmer = self.set_position(i as usize, *base);
                if new_kmer.0 != self.0 {
                    result.push(new_kmer);
                }
            }
        }
        return result;
    }
    fn next_kmers(&self) -> Vec<Self> {
        let mut next : Vec<Self> = Vec::with_capacity(4);
        // quick and dirty implementation for now
        for base in NUCLEOTIDE_BASES {
            let mut unencoded = needletail::bitkmer::bitmer_to_bytes(*self);
            unencoded.push(*base);
            let next_seq = &unencoded[1..=(self.1 as usize)];
            let next_kmer = needletail::bitkmer::BitNuclKmer::new(next_seq, self.1, false).next().unwrap().1;
            next.push(next_kmer);
        }
        return next;
    }

    fn set_position(&self, position: usize, base: u8) -> Self {
        if position >= self.1 as usize {
            panic!("position outside of kmer bounds");
        }
        let mut encoded = self.0;
        // needletail::bitkmer::nuc2bti_lookup_nocheck(base); // package private - urgh!
        let encoded_base_raw = NEEDLETAIL_ENCODING_NUCLEOTIDE_BASES[base as usize].unwrap() as u64;
        // needletail encodes the first base in the LSBs.
        let encoded_position = self.1 as usize - position - 1;
        encoded &= !(3 << (2 * encoded_position));
        // set bit
        encoded |= encoded_base_raw << (2 * encoded_position);
        return (encoded, self.1);
    }

    fn to_kmer_low_bits(&self) -> u64 {
        // just a pass-through since the first base is LSB
        // (c.f. GRIDSS 2-bit encoding scheme)
        return self.0;
    }

    fn kmer_low_bits_to_native_encoding(low_bit_encoded : u64, k : u8) -> Self {
        return (low_bit_encoded, k);
    }
}

fn add_kmer_offset_counts(kmer_counts: &mut Vec<KmerCounts>, k : u8, read : &[u8]) {
    let bnk = needletail::bitkmer::BitNuclKmer::new(read, k, false);
    for kmer in bnk {
        let bit_kmer : BitKmer = kmer.1;
        let base_index = kmer.0;
        let kmer_array_offset = bit_kmer.to_kmer_low_bits() as usize;
        kmer_counts[base_index].count[kmer_array_offset] += 1;
    }
}

/// Counts kmers in each read position for the input fastq/fasta
/// Output is a vector of all possible read kmer offsets
/// each containing the full kmer distribution of all possible kmers in that position
pub(crate) fn count_read_base_kmers(k: u8, input: PathBuf) -> Vec<KmerCounts> {
    let possible_kmers = 1 << (2 * k);
    let mut counts : Vec<KmerCounts> = Vec::with_capacity(301);
    let mut reader = parse_fastx_file(&input).expect("invalid file");
    while let Some(record) = reader.next() {
        let seqrec = record.expect("invalid record");
        let kmer_len = seqrec.num_bases() - k as usize + 1;
        while counts.len() < kmer_len {
            counts.push(KmerCounts { count : vec![0; possible_kmers] } );
        }
        add_kmer_offset_counts(&mut counts, k, seqrec.raw_seq());
    }
    return counts;
}

pub struct KmerCounts {
    pub count : Vec<u64>,
}
#[derive(Debug)]
pub enum KmerState {
    /// Kmer is not found
    NotFound,
    /// Could be a sequencing error
    PossibleSequencingError,
    TODO,
}
pub struct KmerCountStatistics {
    dominant_kmers : Vec<BitKmer>
}

pub fn calculate_state(kmer_offset : usize, k: u8, lookup : &KmerCounts, max_error_rate : f64) -> KmerState {
    let count = lookup.count[kmer_offset];
    if count == 0 || max_error_rate <= 0.0 {
        return KmerState::NotFound;
    }
    let kmer : BitKmer = MutableKmer::kmer_low_bits_to_native_encoding(kmer_offset as u64, k);
    let min_collapse_count : u64 = (count as f64 / max_error_rate).ceil() as u64;
    let neighbourhood_count : u64 = kmer.neighbours()
        .iter()
        .map(|n| lookup.count[MutableKmer::to_kmer_low_bits(n) as usize])
        .sum();
    if neighbourhood_count > min_collapse_count {
        return KmerState::PossibleSequencingError;
    }
    return KmerState::TODO;
}

#[cfg(test)]
mod tests {
    use crate::kmer::MutableKmer;
    #[test]
    fn formatting_assumptions() {
        let rawseq = b"ACGTT";
        let bnk = needletail::bitkmer::BitNuclKmer::new(rawseq, 4, false);
        let result : Vec<Vec<u8>> = bnk.map(|f| needletail::bitkmer::bitmer_to_bytes(f.1)).collect();
        assert_eq!(2, result.len());
        assert_eq!(&rawseq[0..4], result[0]);
        assert_eq!(&rawseq[1..5], result[1]);
    }
    fn test_neighbour(offset : usize, kmer : u8, seq : &[u8], expected : &Vec<&[u8]>) {
        let mut bnk = needletail::bitkmer::BitNuclKmer::new(&seq[offset..(offset + kmer as usize)], kmer, false);
        let n1 = bnk.next().unwrap().1.neighbours();
        let expected_vec : Vec<Vec<u8>> = expected.iter().map(|a| a.iter().cloned().collect()).collect();
        let actual_vec  : Vec<Vec<u8>> = n1.iter().map(|n| needletail::bitkmer::bitmer_to_bytes(*n)).collect();
        assert_eq!(expected_vec, actual_vec);
    }
    fn test_next_kmers(offset : usize, kmer : u8, seq : &[u8], expected : &Vec<&[u8]>) {
        let mut bnk = needletail::bitkmer::BitNuclKmer::new(&seq[offset..(offset + kmer as usize)], kmer, false);
        let n1 = bnk.next().unwrap().1.next_kmers();
        let expected_vec : Vec<Vec<u8>> = expected.iter().map(|a| a.iter().cloned().collect()).collect();
        let actual_vec  : Vec<Vec<u8>> = n1.iter().map(|n| needletail::bitkmer::bitmer_to_bytes(*n)).collect();
        assert_eq!(expected_vec, actual_vec);
    }
    #[test]
    fn neighbours() {
        test_neighbour(0, 4, b"ACGTT",&vec!(
            b"CCGT",
            b"GCGT",
            b"TCGT",
            b"AAGT",
            b"AGGT",
            b"ATGT",
            b"ACAT",
            b"ACCT",
            b"ACTT",
            b"ACGA",
            b"ACGC",
            b"ACGG",
        ));
        test_neighbour(2, 3, b"ACGTT",&vec!(
            b"ATT",
            b"CTT",
            b"TTT",
            b"GAT",
            b"GCT",
            b"GGT",
            b"GTA",
            b"GTC",
            b"GTG",
        ));
    }
    #[test]
    fn next_kmer() {
        test_next_kmers(0, 4, b"ACGTT",&vec!(
            b"CGTA",
            b"CGTC",
            b"CGTG",
            b"CGTT",
        ));
        test_next_kmers(1, 4, b"ACGTT",&vec!(
            b"GTTA",
            b"GTTC",
            b"GTTG",
            b"GTTT",
        ));
    }
}