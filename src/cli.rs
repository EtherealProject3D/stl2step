use crate::converter::{convert, Options};
use std::path::PathBuf;

const HELP: &str = "stl2step - convert STL to STEP (AP214)

USAGE:
    stl2step <input.stl> [-o output.step] [--no-merge] [--scale F] [--angle-tol DEG]

OPTIONS:
    -o, --output <file>   output path (default: input name with .step)
    --no-merge            keep one face per triangle
    --scale <factor>      multiply coordinates (25.4 for inch STLs; output is mm)
    --angle-tol <deg>     coplanarity tolerance in degrees (default 0.1)
    -h, --help            show this help";

pub fn run() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if args.is_empty() || args.iter().any(|a| a == "-h" || a == "--help") {
        println!("{HELP}");
        return;
    }
    let mut input: Option<PathBuf> = None;
    let mut output: Option<PathBuf> = None;
    let mut opts = Options::default();
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "-o" | "--output" => {
                i += 1;
                output = args.get(i).map(PathBuf::from);
            }
            "--no-merge" => opts.merge = false,
            "--scale" => {
                i += 1;
                opts.scale = args.get(i).and_then(|s| s.parse().ok()).unwrap_or_else(|| {
                    eprintln!("--scale needs a number");
                    std::process::exit(2)
                });
            }
            "--angle-tol" => {
                i += 1;
                opts.angle_tol = args.get(i).and_then(|s| s.parse().ok()).unwrap_or_else(|| {
                    eprintln!("--angle-tol needs a number");
                    std::process::exit(2)
                });
            }
            other => input = Some(PathBuf::from(other)),
        }
        i += 1;
    }
    let input = input.unwrap_or_else(|| {
        eprintln!("No input file given.\n\n{HELP}");
        std::process::exit(2)
    });
    let output = output.unwrap_or_else(|| input.with_extension("step"));
    let mut log = |s: String| println!("{s}");
    if let Err(e) = convert(&input, &output, &opts, &mut log) {
        eprintln!("ERROR: {e}");
        std::process::exit(1);
    }
}
