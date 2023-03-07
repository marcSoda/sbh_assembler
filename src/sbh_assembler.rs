use crate::utils;
use rayon::prelude::*;
use std::{
    collections::HashMap,
    cell::RefCell,
    rc::Rc,
    sync::atomic::{ AtomicBool, Ordering }
};

pub struct Assembler {
    pub graph: HashMap<u32, HashMap<u32, Vec<Edge>>>,
    pub nodes: HashMap<u32, Rc<RefCell<Node>>>,
    pub paths: Vec<Vec<Rc<RefCell<Node>>>>,
    pub cycles: Vec<Vec<Rc<RefCell<Node>>>>,
    pub contigs: Vec<Vec<u8>>,
}

impl Assembler {
    // Build the graph
    pub fn new(reads: Vec<Vec<u8>>) -> Self {
        let mut nodes: HashMap<u32, Rc<RefCell<Node>>> = HashMap::new();
        let mut graph: HashMap<u32, HashMap<u32, Vec<Edge>>> = HashMap::new();
        for read in reads.iter() {
            // Get indices from strings
            let pidx = utils::vec2idx(read, NodeType::Prefix);
            let sidx = utils::vec2idx(read, NodeType::Suffix);
            // Get nodes from prefixes or create them, setting odeg and ideg accordingly
            let prefix = nodes.entry(pidx)
                .and_modify(|n| { n.borrow_mut().odeg+=1; })
                .or_insert_with(|| { Node::new(pidx, 0, 1) })
                .clone();
            let suffix = nodes.entry(sidx)
                .and_modify(|n| { n.borrow_mut().ideg+=1; })
                .or_insert_with(|| { Node::new(sidx, 1, 0) })
                .clone();
            // Insert edge
            graph.entry(pidx)
                .or_insert_with(HashMap::new)
                .entry(sidx)
                .or_insert_with(Vec::new)
                .push(Edge::new(prefix, suffix));
        }
        Assembler {
            graph,
            nodes,
            paths: Vec::default(),
            cycles: Vec::default(),
            contigs: Vec::default(),
        }
    }

    // Find all paths or cycles depending on the type requested
    pub fn populate_paths_or_cycles(&mut self, typ: PathType) {
        // Get all valid start nodes depending on the type requested
        let starts: Vec<Rc<RefCell<Node>>> = self.nodes.values().filter_map(|n| {
            let node = n.borrow();
            match typ {
                // Paths only start where outdegree > indegree
                PathType::Path if node.odeg > node.ideg => Some(n.clone()),
                // Cycles only start where outdegree and indegree is positive
                PathType::Cycle if node.odeg > 0 && node.ideg > 0 => Some(n.clone()),
                _ => None,
            }
        }).collect();
        // Get all paths or cycles and populate their respective vector
        for start in starts {
            let p = self.find_path_or_cycle(start.clone(), typ);
            match typ {
                PathType::Path => if p.len() >= 5 { self.paths.push(p); },
                PathType::Cycle => if p.len() >= 3 { self.cycles.push(p); },
            }
        }
    }

    // Find the path or cycle that starts at the start node if it exists
    fn find_path_or_cycle(&mut self, start: Rc<RefCell<Node>>, typ: PathType) -> Vec<Rc<RefCell<Node>>> {
        let mut path = vec![start.clone()];
        let mut current = start.clone();
        loop {
            // Get suffix graph associated with the start node
            let sufs = match self.graph.get_mut(&current.borrow().idx) {
                Some(m) => m,
                None => break,
            };
            // Get the index of the next suffix
            let next = sufs.keys().find(|idx| {
                match sufs.get(idx) {
                    Some(idx) => idx,
                    None => return false,
                }.iter().all(|e| !e.used)
            });
            match next {
                Some(&idx) => {
                    let edge = sufs.get_mut(&idx).unwrap().iter_mut().find(|e| !e.used).unwrap();
                    edge.mark_used();
                    let node = self.nodes.get(&idx).unwrap().clone();
                    path.push(node.clone());
                    // The only difference between a path and a cycle is a cycle stops when we get back to the start node
                    if matches!(typ, PathType::Cycle) && node == start { break; }
                    current = node;
                } None => break,
            }
        }
        path
    }

    // Convert paths and cycles to contigs
    pub fn paths_cycles_to_contigs(&mut self) {
        // Chain the paths and cycles into one vector
        for path_or_cycle in self.paths.iter().chain(self.cycles.iter()) {
            let mut contig = Vec::new();
            for node in path_or_cycle {
                let seq = utils::idx2vec(node.borrow().idx, 15);
                contig.extend_from_slice(&seq);
            }
            self.contigs.push(contig);
        }
    }

    // Remove all contigs that are completely encompassed by another contig
    // This method is parallalized, making it orders of magnitudes faster for large conig arrays
    // Returns the number of contigs that was removes
    pub fn remove_contained_contigs(&mut self) -> usize {
        // Faster if sorted
        self.contigs.sort_unstable_by_key(|c| std::cmp::Reverse(c.len()));
        // Describes which contigs should be removed
        let to_remove = Vec::from_iter((0..self.contigs.len()).map(|_| AtomicBool::new(false)));
        self.contigs.par_iter().enumerate().for_each(|(i, contig_i)| {
            if to_remove[i].load(Ordering::SeqCst) { return }
            for j in (i + 1)..self.contigs.len() {
                if to_remove[j].load(Ordering::SeqCst) { continue; }
                if contig_i.len() > self.contigs[j].len() && Self::is_contig_contains(contig_i.as_slice(), self.contigs[j].as_slice()) {
                    to_remove[j].store(true, Ordering::SeqCst);
                } else if Self::is_contig_contains(self.contigs[j].as_slice(), contig_i.as_slice()) {
                    to_remove[i].store(true, Ordering::SeqCst);
                    break;
                }
            }
        });
        // Remove the contigs that were marked in to_remove
        self.contigs = self
            .contigs
            .iter()
            .enumerate()
            .filter_map(|(i, v)| if !to_remove[i].load(Ordering::SeqCst) { Some(v.clone()) } else { None })
            .collect();
        // return number of removed contigs
        to_remove.iter().filter(|&b| b.load(Ordering::SeqCst)).count()
    }

    // Returns true if sub_contig is contained in contig
    fn is_contig_contains(contig: &[u8], sub_contig: &[u8]) -> bool {
        if contig.len() < sub_contig.len() { return false; }
        let mut pos = 0;
        while pos + sub_contig.len() <= contig.len() {
            if contig[pos..pos + sub_contig.len()] == sub_contig[..] { return true; }
            pos += 1;
        }
        false
    }

    // Merges contigs if they overlap
    // This method is parallalized, making it orders of magnitudes faster for large conig arrays
    // Returns the number of contigs that were merged
    pub fn merge_contigs(&mut self, min_overlap: usize) -> usize {
        let mut merged = 0;
        self.contigs.sort_unstable_by_key(|contig| std::cmp::Reverse(contig.len()));
        let mut i = 0;
        while i < self.contigs.len() {
            let j_range = i + 1..self.contigs.len();
            let overlaps = j_range.clone().into_par_iter().map(|j| {
                Self::merge_if_overlap(&self.contigs[i], &self.contigs[j], min_overlap)
            }).collect::<Vec<_>>();
            for (j, overlap) in (j_range).zip(overlaps.into_iter()) {
                if let Some((_, new_contig)) = overlap {
                    self.contigs.swap_remove(j);
                    self.contigs.swap_remove(i);
                    self.contigs.push(new_contig);
                    i = 0;
                    merged+=2;
                    break;
                }
            }
            i += 1;
        }
        merged
    }

    // merges c1 and c2 if they overlap
    fn merge_if_overlap(c1: &[u8], c2: &[u8], min_overlap_len: usize) -> Option<(usize, Vec<u8>)> {
        let overlap_range = c1.len().min(c2.len());
        // Check for overlap at beginning of c1 and end of c2
        for i in 0..overlap_range {
            if c1.starts_with(&c2[c2.len() - i - 1..]) {
                let overlap_len = i + 1;
                if overlap_len >= min_overlap_len {
                    let new_contig = [&c2[..c2.len() - overlap_len], c1].concat();
                    return Some((overlap_len, new_contig));
                }
            }
            if c2.starts_with(&c1[c1.len() - i - 1..]) {
                let overlap_len = i + 1;
                if overlap_len >= min_overlap_len {
                    let new_contig = [&c1[..c1.len() - overlap_len], c2].concat();
                    return Some((overlap_len, new_contig));
                }
            }
        }
        None
    }
}

#[derive(PartialEq)]
pub struct Node {
    pub idx: u32,
    pub ideg: usize,
    pub odeg: usize,
}

impl Node {
    pub fn new(idx: u32, ideg: usize, odeg: usize) -> Rc<RefCell<Self>> {
        Rc::new(RefCell::new(Node {
            idx,
            ideg,
            odeg,
        }))
    }
}

pub struct Edge {
    pub prefix: Rc<RefCell<Node>>,
    pub suffix: Rc<RefCell<Node>>,
    pub used: bool,
}

impl Edge {
    pub fn new(prefix: Rc<RefCell<Node>>, suffix: Rc<RefCell<Node>>) -> Self {
        Edge {
            prefix,
            suffix,
            used: false,
        }
    }

    pub fn mark_used(&mut self) {
        self.prefix.borrow_mut().odeg -= 1;
        self.suffix.borrow_mut().ideg -= 1;
        self.used = true;
    }
}

#[derive(Copy, Clone)]
pub enum PathType {
    Path,
    Cycle,
}

pub enum NodeType {
    Prefix,
    Suffix,
}
