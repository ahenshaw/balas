use anyhow::Result;
use argh::FromArgs;
use balas::Balas;
use std::fs::File;
use std::io::Write;
use std::path::PathBuf;
use std::time::Instant;

#[derive(FromArgs)]
#[argh(help_triggers("-h", "--help"))]
/// Solve a Binary-Variable Linear Program
struct Args {
    /// input file in LP format
    #[argh(positional)]
    infile: PathBuf,

    /// number of threads to use
    #[argh(option, short = 't')]
    threads: Option<usize>,

    /// optional recording file
    #[argh(option, short = 'o')]
    outfile: Option<PathBuf>,

    /// use this heuristic pre-solve
    #[argh(option)]
    heuristic: Option<f64>,

    /// use the original recursive code
    #[argh(switch)]
    recursive: bool,
}

fn main() -> Result<()> {
    let args: Args = argh::from_env();

    let mut balas = Balas::from_lp(&args.infile)?;

    let start = Instant::now();
        if let Some(heuristic) = args.heuristic {
            balas.best = heuristic;
        }
        if args.recursive {
            balas.solve_recursively();
        } else {
                let num_threads = match args.threads {
                    Some(n) => n.max(1),
                    None => match std::thread::available_parallelism() {
                        Ok(n) => n.into(),
                        _ => 1usize,
                    },
                };
                // num_threads needs to be a power of two
                let used_threads = 1usize << num_threads.ilog2();
                println!("Using {used_threads} thread{}", if used_threads != 1 {"s"} else {""});
                balas.solve(used_threads);
        }
    println!(
        "Elapsed time: {:?}",
        Instant::now() - start,
    );

    balas.report();

    if let Some(outfile) = args.outfile {
        let mut out = File::create(outfile)?;
        let buf = serde_json::to_string(&balas)?;
        let _ = out.write(buf.as_bytes())?;
    }

    Ok(())
}
