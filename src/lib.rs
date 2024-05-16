mod lp_errors;
mod lp_reader;
mod recursive_solver;

use num::Bounded;
use serde::{Deserialize, Serialize};
use std::sync::{Arc, RwLock};
use std::{fmt::Display, ops::Neg};

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
    pub fn reset(&mut self) {
        self.count = 0;
        self.best = T::max_value();
        self.solution = Vec::new();
    }

    /// single-threaded solver
    pub fn solve(&mut self) {
        let global_best = Arc::new(RwLock::new(T::max_value()));

        let (best, count, solution) = Self::solve_subtree(0, 0, global_best, &self.fixed);
        self.best = best;
        self.count = count;
        self.solution = solution;
    }

    /// multi-threaded solver
    pub fn solve_mt(&mut self, _num_threads: usize)
    where
        T: std::marker::Send + std::marker::Sync + 'static,
    {
        let fixed = Arc::new(self.fixed.clone());
        let fixed1 = Arc::clone(&fixed);
        let fixed2 = Arc::clone(&fixed);
        let global_best = Arc::new(RwLock::new(T::max_value()));
        let gb1 = Arc::clone(&global_best);
        let gb2 = Arc::clone(&global_best);
        let h1 = std::thread::spawn(move || Self::solve_subtree(1, 0, gb1, &fixed1));
        let h2 = std::thread::spawn(move || Self::solve_subtree(1, 1, gb2, &fixed2));
        let (b1, c1, s1) = h1.join().unwrap();
        let (b2, c2, s2) = h2.join().unwrap();
        // let (b1, c1, s1) = Self::solve_subtree(1, 0, gb1, &fixed1);
        // let (b2, c2, s2) = Self::solve_subtree(1, 1, gb2, &fixed2);

        self.count = c1 + c2;
        if b1 <= b2 {
            self.best = b1;
            self.solution = s1;
        } else {
            self.best = b2;
            self.solution = s2;
        };
    }

    fn init_subtree(start_var_index: usize, tree_index: usize, fixed: &Fixed<T>) -> (Vec<T>, T) {
        // Initialize the constraint accumulator with the negation of the b vector (the
        // right-hand side of the constraints).  This way, we can just compare against 0
        // later on.
        let mut accumulator: Vec<T> = fixed.rhs.iter().map(|&b| -b).collect();
        let mut objective = T::zero();

        // to init a subtree, we need to calculate the state for all of the ancestor nodes (if any)
        for var_index in 0..start_var_index {
            let branch = tree_index >> (start_var_index - var_index - 1) & 1;
            let constraints = &fixed.constraints[var_index];
            let coefficient = fixed.coefficients[var_index];

            if branch == 1 {
                objective += &coefficient;
                accumulator
                    .iter_mut()
                    .zip(constraints)
                    .for_each(|(a, b)| *a += &b);
            }
        }

        (accumulator, objective)
    }

    fn solve_subtree(
        start_var_index: usize,
        tree_index: usize,
        global_best: Arc<RwLock<T>>,
        fixed: &Fixed<T>,
    ) -> (T, usize, Vec<u8>) {
        let mut vars: Vec<u8> = vec![0; fixed.num_vars];
        let mut var_index = start_var_index;
        let mut state = Flow::Normal;
        let mut count = 0usize;
        let mut best: T = *global_best.read().unwrap();
        let mut solution = Vec::<u8>::new();

        let (mut accumulator, mut objective) =
            Self::init_subtree(start_var_index, tree_index, &fixed);

        let mut branch = 0u8;

        loop {
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
                                .for_each(|(a, b)| *a -= &b);
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
                        vars[var_index] = branch;
                        // Update the accumulator.  This only needs to be done in the ones branch
                        accumulator
                            .iter_mut()
                            .zip(constraints)
                            .for_each(|(a, b)| *a += &b);

                        // Update the current value of the objective
                        objective += &coefficient;

                        // If we're already not better than the current best objective, then
                        // we can prune this entire branch.
                        if (objective >= best)
                            && (objective
                                >= match global_best.try_read() {
                                    Ok(gl) => {
                                        best = *gl;
                                        best
                                    }
                                    _ => best,
                                })
                        {
                            state = Flow::Backtrack;
                            continue;
                        } else {
                            // Check if constraints satisfied, while updating the accumulator.
                            // We do not have to check the 0 branch, as the accumulator is not changed there.
                            // If all of constraints are satisfied, then we are fathomed and we can't do any better.
                            if accumulator.iter().all(|x| *x >= T::zero()) {
                                best = objective;
                                if let Ok(mut gl) = global_best.try_write() {
                                    *gl = best;
                                }

                                // println!("{objective} {:?}", &vars[..=index]);
                                solution = vars.clone();
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
            for (var, value) in self.vars.iter().zip(self.solution.iter()) {
                println!("  {var}: {}", value);
            }
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
