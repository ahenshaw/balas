use anyhow::Result;
use argh::FromArgs;
use balas::Balas;
use std::{path::PathBuf, time::Instant};

#[derive(FromArgs)]
/// Solve a Binary-Variable Linear Program
struct Args {
    /// input file in LP format
    #[argh(positional)]
    infile: PathBuf,
    /// how many repetitions (for timing)
    #[argh(option, short = 'r', default = "1")]
    reps: usize,
}

fn main() -> Result<()> {
    let args: Args = argh::from_env();

    let mut balas = Balas::from_lp(&args.infile)?;

    let start = Instant::now();
    for _ in 0..args.reps {
        balas.reset();
        balas.solve();
    }
    println!("Elapsed time: {:?} (repetitions: {})", Instant::now() - start, args.reps);
    balas.report();

    Ok(())
}
