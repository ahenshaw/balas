use anyhow::Result;
use argh::FromArgs;
use balas::Balas;
// use lp_parser_rs::model::lp_problem::LPProblem;
use lp_parser_rs::parse::parse_lp_file;
use std::fs;
use std::path::PathBuf;
use std::time::Instant;

#[derive(FromArgs)]
/// Solve a Binary-Variable Linear Program
struct Args {
    /// input file in LP format
    #[argh(positional)]
    infile: PathBuf,
}

fn main() -> Result<()> {
    let args: Args = argh::from_env();
    let code = fs::read_to_string(args.infile)?;
    let lp = parse_lp_file(&code)?;

    let mut balas = Balas::from_lp(&lp)?;

    let start = Instant::now();
    // for _ in 0..100_000 {
    balas.solve();
    // }
    println!("Elapsed time: {:?}", Instant::now() - start);
    balas.report();

    Ok(())
}
