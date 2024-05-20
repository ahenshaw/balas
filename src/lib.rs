mod lp_errors;
mod lp_reader;
mod recursive_solver;

use num::Bounded;
use serde::{Deserialize, Serialize};
use std::sync::{Arc, RwLock};
use std::{fmt::Display, ops::Neg};
use rayon::prelude::*;


type Array<T> = Vec<Vec<T>>;

#[derive(Deserialize, Serialize)]
pub struct Balas<T> {
    pub best: T,
    pub solution: Vec<u8>,
    pub count: usize,
    vars: Vec<String>,
    pub recording: Vec<Record>,
    fixed: Fixed<T>,
}

/// The unchanging, fixed variables representing the BIP
#[derive(Deserialize, Serialize, Clone)]
struct Fixed<T> {
    num_vars: usize,
    coefficients: Vec<T>,
    constraints: Array<T>,
    rhs: Vec<T>,
    cumulative: Array<T>,
}

impl<T> Balas<T>
where
    T: Bounded
        + std::convert::Into<f64>
        + Neg
        + Copy
        + Display
        + num::Zero
        + for<'a> std::ops::AddAssign<&'a T>
        + for<'a> std::ops::SubAssign<&'a T>
        + std::cmp::PartialOrd
        + std::fmt::Debug,
    Vec<T>: FromIterator<<T as Neg>::Output>,
{
    pub fn new(coeff: &[T], constraints: &Array<T>, b: &[T], vars: &Vec<String>) -> Balas<T> {
        let cumulative = Self::make_cumulative(constraints);
        let fixed = Fixed {
            num_vars: coeff.len(),
            coefficients: coeff.to_vec(),
            constraints: constraints.clone(),
            rhs: b.to_vec(),
            cumulative,
        };
        Balas {
            best: T::max_value(),
            solution: Vec::new(),
            count: 0,
            vars: vars.to_owned(),
            recording: vec![],
            fixed,
        }
    }

    /// multi-threaded solver
    pub fn solve(&mut self, num_threads: usize)
    where
        T: std::marker::Send + std::marker::Sync + 'static,
    {
        let fixed = Arc::new(self.fixed.clone());
        let global_best = Arc::new(RwLock::new(self.best));
        let start_index = num_threads.ilog2() as usize;

        let handles: Vec<(T, usize, Vec<u8>)> = (0..num_threads).into_par_iter().map(|i| {
            let f = Arc::clone(&fixed);
            let gb = Arc::clone(&global_best);
            Self::solve_subtree(start_index, i, gb, &f)
        }).collect();

        for handle in handles {
            let (best, count, solution) = handle;
            self.count += count;
            if best < self.best {
                self.solution = solution;
                self.best = best;
            }
        }
    }

    fn init_subtree(start_var_index: usize, tree_index: usize, fixed: &Fixed<T>) -> (Vec<T>, T, Vec<u8>) {
        // Initialize the constraint accumulator with the negation of the b vector (the
        // right-hand side of the constraints).  This way, we can just compare against 0
        // later on.
        let mut accumulator: Vec<T> = fixed.rhs.iter().map(|&b| -b).collect();
        let mut objective = T::zero();

        let mut branches = tree_index;
        let mut vars = vec![];

        // to init a subtree, we need to calculate the state for all of the ancestor nodes (if any)

        for var_index in 0..start_var_index {
            let branch = (branches & 1) as u8;
            let constraints = &fixed.constraints[var_index];
            let coefficient = fixed.coefficients[var_index];

            if branch == 1 {
                objective += &coefficient;
                accumulator
                    .iter_mut()
                    .zip(constraints)
                    .for_each(|(a, b)| *a += b);
            }
            vars.push(branch);
            branches >>= 1;
        }

        (accumulator, objective, vars)
    }

    fn solve_subtree(
        start_var_index: usize,
        tree_index: usize,
        global_best: Arc<RwLock<T>>,
        fixed: &Fixed<T>,
    ) -> (T, usize, Vec<u8>) {
        let mut var_index = start_var_index;
        let mut state = Flow::Normal;
        let mut count = 0usize;
        let mut best = T::max_value();
        let mut solution = Vec::<u8>::new();

        let (mut accumulator, mut objective, mut vars) =
            Self::init_subtree(start_var_index, tree_index, fixed);
        vars.resize(fixed.num_vars, 0);

        let mut branch = 0u8;

        let mut min_index = usize::MAX;

        loop {
            min_index = min_index.min(var_index);
            // Alias the current column of the constraints and grab the coefficients value
            let constraints = &fixed.constraints[var_index];
            let coefficient = fixed.coefficients[var_index];

            match state {
                Flow::Terminate => break,
                Flow::Backtrack => {
                    if vars[var_index] == 1 {
                        if var_index == start_var_index {
                            state = Flow::Terminate;
                        } else {
                            // we have to reverse what we did before we leave
                            accumulator
                                .iter_mut()
                                .zip(constraints)
                                .for_each(|(a, b)| *a -= b);
                            objective -= &coefficient;
                            vars[var_index] = 0;
                            var_index -= 1;
                        }
                    } else {
                        state = Flow::Normal;
                        branch = 1;
                    }
                }
                Flow::Normal => {
                    count += 1;

                    if branch == 1 {
                        vars[var_index] = 1;
                        // Update the accumulator.  This only needs to be done in the ones branch
                        accumulator
                            .iter_mut()
                            .zip(constraints)
                            .for_each(|(a, b)| *a += b);

                        // Update the current value of the objective
                        objective += &coefficient;

                        if (objective >= best) || (objective >= *global_best.read().unwrap())
                        {
                            state = Flow::Backtrack;
                            continue;
                        } else {
                            // Check if constraints satisfied, while updating the accumulator.
                            // We do not have to check the 0 branch, as the accumulator is not changed there.
                            // If all of constraints are satisfied, then we are fathomed and we can't do any better.
                            if accumulator.iter().all(|x| *x >= T::zero()) {
                                best = objective;
                                loop {
                                    if let Ok(mut gl) = global_best.try_write() {
                                        if best < * gl {
                                            *gl = best;
                                            println!("{best}");
                                            break;
                                        }
                                    }
                                }
                                // println!("{objective} {:?}", &vars[..=index]);
                                solution.clone_from(&vars);
                                state = Flow::Backtrack;
                                continue;
                            }
                        }
                    }

                    // This is the same behavior for either a 0 or 1 branch
                    // If there is a potentially feasible descendant, then keep descending the tree
                    if let Some(ccons) = fixed.cumulative.get(var_index) {
                        if accumulator
                            .iter()
                            .zip(ccons)
                            .all(|(&a, &b)| a + b >= T::zero())
                        {
                            var_index += 1;
                            branch = 0;
                        } else {
                            state = Flow::Backtrack;
                        }
                    } else {
                        state = Flow::Backtrack;
                    }
                }
            }
        }
        (best, count, solution)
    }

    fn record(&mut self, label: &str, state: NodeState) {
        self.recording.push(Record {
            node: label.to_string(),
            state,
        });
    }

    pub fn report(&self) {
        if self.best != T::max_value() {
            println!("Optimal value: {}", self.best);
            println!("Solution:");
            for (i, value) in self.solution.iter().enumerate() {
                print!("{value}");
                if i % 4 == 3 {print!(" ")}
            }
            println!();
        } else {
            println!("No solution");
        }

        println!("Examined {:?} nodes", self.count);
    }

    pub fn make_cumulative(constraints: &Array<T>) -> Array<T> {
        let num_cols = constraints[0].len();
        let mut running_total = vec![T::zero(); num_cols];
        let mut cumulative = vec![];

        for row in constraints.iter().skip(1).rev() {
            for (i, val) in row.iter().enumerate() {
                if *val > T::zero() {
                    running_total[i] += val;
                }
            }
            cumulative.push(running_total.clone())
        }
        cumulative.reverse();
        cumulative
    }
}

#[derive(Serialize, Deserialize, Clone, Debug)]
enum Flow {
    Terminate,
    Backtrack,
    Normal,
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub enum NodeState {
    Default,
    Active,
    Visited,
    Fathomed,
    Infeasible,
    Suboptimal,
    ImpossibleChildren,
    Skipped,
}

#[derive(Serialize, Deserialize, Clone)]
pub struct Record {
    pub node: String,
    pub state: NodeState,
}
