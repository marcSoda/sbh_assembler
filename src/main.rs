mod sbh_assembler;
mod utils;
use sbh_assembler::{ Assembler, PathType };
use std::{ time::Instant, env };

fn main() {
    let time_start = Instant::now();
    let infile = match env::args().nth(1) {
        Some(a) => a,
        None => "data/YeastReads.fasta".to_string(),
    };
    let outfile = match env::args().nth(2) {
        Some(a) => a,
        None => "cont.fasta".to_string(),
    };

    println!("\nRunning the assembler with infile: \x1b[32m{}\x1b[0m and outfile: \x1b[32m{}\x1b[0m", infile, outfile);
    println!("If you would like to run with different files, use the program like this");
    println!("\t\x1b[32msbh <infile_path> <outfile_path>\x1b[0m");
    println!("\tor, if you do not have an executable, you will have to recompile:");
    println!("\t\x1b[32mcargo run -- <infile_path> <outfile_path>\x1b[0m\n");

    let reads = utils::fasta_reader(&infile);
    let mut ass = Assembler::new(reads);

    println!("Populating Paths................................");
    ass.populate_paths_or_cycles(PathType::Path);
    let lpath = ass.paths.iter()
        .max_by(|a, b| a.len().cmp(&b.len()))
        .cloned()
        .unwrap_or_default();
    println!("Generated \x1b[32m{}\x1b[0m total paths.", ass.paths.len());
    println!("\tLongest generated path was \x1b[32m{}\x1b[0m nodes.", lpath.len());

    println!("Populating Cycles...............................");
    ass.populate_paths_or_cycles(PathType::Cycle);
    let lcycle = ass.cycles.iter()
        .max_by_key(|a| a.len())
        .cloned()
        .unwrap_or_default();
    println!("Generated \x1b[32m{}\x1b[0m total cycles.", ass.cycles.len());
    println!("\tLongest generated cycle was \x1b[32m{}\x1b[0m nodes.", lcycle.len());

    println!("Converting the paths and cycles to contigs......");
    ass.paths_cycles_to_contigs();
    println!("Generated \x1b[32m{}\x1b[0m contigs.", ass.contigs.len());

    println!("Condensing contigs..............................");
    let mut prev = usize::MAX;
    loop {
        println!("\tRemoving Contained Contigs..............");
        let removed = ass.remove_contained_contigs();
        println!("\t\tRemoved \x1b[32m{}\x1b[0m contained contigs.", removed);
        println!("\tMerging contigs. May take some time.....");
        let merged = ass.merge_contigs(15);
        println!("\t\tMerged \x1b[32m{}\x1b[0m contigs.", merged);
        if prev == ass.contigs.len() { break; }
        prev = ass.contigs.len();
    }
    println!("Successfully condensed to \x1b[32m{}\x1b[0m contigs.", ass.contigs.len());

    let lcont = ass.contigs.iter()
        .max_by(|a, b| a.len().cmp(&b.len()))
        .cloned()
        .unwrap_or_default();
    println!("Longest generated contig was \x1b[32m{}\x1b[0m nucleotides.", lcont.len());

    println!("Writing contigs to \x1b[32m{}\x1b[0m...................", outfile);
    match utils::cont2file(&outfile, ass.contigs) {
        Ok(_) => println!("Successfully wrote contigs to \x1b[32m{}\x1b[0m", outfile),
        Err(_) => {
            eprintln!("\x1b[31mThere was an error writing to {}.\n
                Please email me at masa20@lehigh.edu before you give me a 0.\x1b[0m", outfile);
            std::process::exit(1);
        }
    }
    let duration = time_start.elapsed();
    println!("Assembly took \x1b[32m{}\x1b[0m seconds.", duration.as_secs_f32());
}
