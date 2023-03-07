use crate::sbh_assembler::NodeType;
use std::io::{ BufRead, BufReader, BufWriter, Write };
use std::fs::File;

// Read a fasta file
// TODO: lazy format checking. Ensure input is a fasta file
pub fn fasta_reader(fname: &str) -> Vec<Vec<u8>> {
    let file = match File::open(fname) {
        Ok(f) => f,
        Err(_) => {
            println!("\x1b[31mFATAL: Failed to open file: '{}'.\x1b[0m", fname);
            std::process::exit(1);
        }
    };
    let reader = BufReader::new(file);
    let mut reads: Vec<Vec<u8>> = vec![];
    for line in reader.lines() {
        let line = line.unwrap();
        if line.starts_with('>') { continue; }
        if line.len() != 30 {
            continue;
        }
        reads.push(line.as_bytes().to_vec());
    }
    reads
}

// Convert a sequence vec to an index
pub fn vec2idx(read: &Vec<u8>, t: NodeType) -> u32 {
    let iter: Box<dyn Iterator<Item=&u8>> = match t {
        NodeType::Prefix => Box::new(read.iter().take(15)),
        NodeType::Suffix => Box::new(read.iter().skip(read.len()-15)),
    };
    let mut idx = 0;
    for (i, c) in iter.enumerate() {
        let v = match c {
            65 => 0,
            67 => 1,
            71 => 2,
            84 => 3,
            _ => panic!("VALS CAN ONLY BE A, C, G, or T"),
        };
        idx += v * 4u32.pow(i as u32);
    }
    idx
}

// Convert a index to a sequence vector
pub fn idx2vec(mut idx: u32, length: usize) -> Vec<u8> {
    let mut result = Vec::with_capacity(length);
    for _ in 0..length {
        let c = match idx % 4 {
            0 => 65_u8,
            1 => 67_u8,
            2 => 71_u8,
            3 => 84_u8,
            _ => panic!("IMPOSSIBLE VALUE FOR MOD 4"),
        };
        result.push(c);
        idx /= 4;
    }
    result
}

// Write contigs to file in fasta format
pub fn cont2file(fname: &str, contigs: Vec<Vec<u8>>) -> std::io::Result<()> {
    let file = File::create(fname)?;
    let mut writer = BufWriter::new(file);
    for (i, cont) in contigs.iter().enumerate() {
        writeln!(writer, ">sequence{}", i + 1)?;
        let cont_str = String::from_utf8(cont.to_vec())
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
        writeln!(writer, "{}", cont_str)?;
    }

    Ok(())
}
