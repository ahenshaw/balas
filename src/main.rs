use anyhow::Result;
use argh::FromArgs;
use balas::Balas;
use std::fs::File;
use std::io::Write;
use std::path::PathBuf;
use std::time::Instant;

#[derive(FromArgs)]
/// Solve a Binary-Variable Linear Program
struct Args {
    /// input file in LP format
    #[argh(positional)]
    infile: PathBuf,

    /// how many repetitions (for timing)
    #[argh(option, short = 'r', default = "1")]
    reps: usize,

    /// optional recording file
    #[argh(option)]
    outfile: Option<PathBuf>,

    /// use this heuristic pre-solve
    #[argh(option)]
    heuristic: Option<f64>,

    /// use the original recursive code
    #[argh(switch)]
    recursive: bool,
    /// single-threaded mode
    #[argh(switch, short = 's')]
    single: bool,
    /// num threads
    #[argh(option, short = 't')]
    threads: Option<usize>,

}

fn main() -> Result<()> {
    let args: Args = argh::from_env();

    let mut balas = Balas::from_lp(&args.infile)?;

    let start = Instant::now();
    for _ in 0..args.reps {
        balas.reset();
        if let Some(heuristic) = args.heuristic {
            balas.best = heuristic.into();
        }
        if args.recursive {
            balas.solve_recursively();
        } else {
            if args.single {
                balas.solve()
            } else {
                let num_threads = match args.threads {
                    Some(n) => n,
                    None => 2,
                };
                balas.solve_mt(num_threads);
            }
        }
    }
    println!(
        "Elapsed time: {:?} (repetitions: {})",
        Instant::now() - start,
        args.reps
    );
    balas.report();
    if let Some(outfile) = args.outfile {
        let mut out = File::create(outfile)?;
        let buf = serde_json::to_string(&balas)?;
        out.write(buf.as_bytes())?;
    }

    Ok(())
}
