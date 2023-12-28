use balas::Balas;
use std::time::Instant;

fn main() {
    let z: Vec<f32> = vec![3.0, 5.0, 6.0, 9.0, 10.0, 10.0];
    let constraints: Vec<Vec<f32>> = vec![
        vec![-2.0, -5.0, 5.0],
        vec![6.0, -3.0, -1.0],
        vec![-3.0, 1.0, 4.0],
        vec![4.0, 3.0, -2.0],
        vec![1.0, -2.0, 2.0],
        vec![-2.0, 1.0, -1.0],
    ];
    let b: Vec<f32> = vec![2.0, -2.0, 3.0];
    let mut balas = Balas::new(&z, &constraints, &b);

    let start = Instant::now();
    for _ in 0..1000 {
        balas.solve();
    }
    println!("Elapsed time: {:?}", Instant::now() - start);
    balas.report();
}
