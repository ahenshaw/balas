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

// fn display(lp: &LPProblem) {
//     println!("Subject To");
//     for constraint in &lp.constraints {
//         println!("{constraint:?}");
//         println!();
//     }
// }
fn main() -> Result<()> {
    let args: Args = argh::from_env();
    let code = fs::read_to_string(args.infile)?;
    let lp = parse_lp_file(&code)?;
    // display(&lp);

    // let z: Vec<f32> = vec![3.0, 5.0, 6.0, 9.0, 10.0, 10.0];
    // let constraints: Vec<Vec<f32>> = vec![
    //     vec![-2.0, -5.0, 5.0],
    //     vec![6.0, -3.0, -1.0],
    //     vec![-3.0, 1.0, 4.0],
    //     vec![4.0, 3.0, -2.0],
    //     vec![1.0, -2.0, 2.0],
    //     vec![-2.0, 1.0, -1.0],
    // ];
    // let b: Vec<f32> = vec![2.0, -2.0, 3.0];
    // let mut balas = Balas::new(&z, &constraints, &b);

    let mut balas = Balas::from_lp(&lp)?;

    let start = Instant::now();
    for _ in 0..100_000 {
        balas.solve();
    }
    println!("Elapsed time: {:?}", Instant::now() - start);
    balas.report();

    Ok(())
}
