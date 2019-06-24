extern crate mesh_simplification;

use mesh_simplification::{Flt, Mesh};

use std::env;

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() != 4 {
        println!("Invalid arguments.\n./mesh_simplication [in.obj] [out.obj] 0.3");
    } else {
        Mesh::new(&args[1])
            .simplify(args[3].parse::<Flt>().unwrap())
            .save(&args[2]);
    }
}
