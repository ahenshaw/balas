mod lp_errors;
mod lp_reader;
mod recursive_solver;

use num::Bounded;
use serde::{Deserialize, Serialize};
use std::{fmt::Display, ops::Neg};

type Array<T> = Vec<Vec<T>>;

#[derive(Deserialize, Serialize)]
pub struct Balas<T> {
    pub best: T,
    pub solution: Vec<u8>,
    pub count: usize,
    vars: Vec<String>,
    pub recording: Vec<Record>,
    pub fixed: Fixed<T>,
}

/// The unchanging, fixed variables representing the BIP
#[derive(Deserialize, Serialize)]
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
        let fixed = Fixed{
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

    pub fn solve(&mut self) {

    }

    fn tree(start: usize, fixed: &Fixed<T>) -> (T, usize, Vec<u8>) {

        let mut vars: Vec<u8> = vec![0; fixed.num_vars];
        let mut branch = 0u8;
        let mut index = start;
        let mut objective = T::zero();
        let mut state = Flow::Normal;
        let mut count = 0usize;
        let mut best = T::max_value();
        let mut solution = Vec::<u8>::new();

        // Initialize the constraint accumulator with the negation of the b vector (the
        // right-hand side of the constraints).  This way, we can just compare against 0
        // later on.
        let mut accumulator: Vec<T> = fixed.rhs.iter().map(|&b| -b).collect();

        loop {
            // Alias the current column of the constraints and grab the coefficients value
            let cons = &fixed.constraints[index];
            let coeff = fixed.coefficients[index];

            match state {
                Flow::Terminate => break,
                Flow::Backtrack => {
                    if vars[index] == 1 {
                        if index == 0 {
                            state = Flow::Terminate;
                        } else {
                            // we have to reverse what we did before we leave
                            accumulator.iter_mut().zip(cons).for_each(|(a, b)| *a -= &b);
                            objective -= &coeff;
                            vars[index] = 0;
                            index -= 1;
                        }
                    } else {
                        state = Flow::Normal;
                        branch = 1;
                    }
                }
                Flow::Normal => {
                    count += 1;

                    if branch == 1 {
                        vars[index] = branch;
                        // Update the accumulator.  This only needs to be done in the ones branch
                        accumulator.iter_mut().zip(cons).for_each(|(a, b)| *a += &b);

                        // Update the current value of the objective
                        objective += &coeff;

                        // If we're already not better than the current best objective, then
                        // we can prune this entire branch.
                        if objective >= best {
                            state = Flow::Backtrack;
                            continue;
                        } else {
                            // Check if constraints satisfied, while updating the accumulator.
                            // We do not have to check the 0 branch, as the accumulator is not changed there.
                            // If all of constraints are satisfied, then we are fathomed and we can't do any better.
                            if accumulator.iter().all(|x| *x >= T::zero()) {
                                best = objective;
                                // println!("{objective} {:?}", &vars[..=index]);
                                solution = vars.clone();
                                state = Flow::Backtrack;
                                continue;
                            }
                        }
                    }

                    // This is the same behavior for either a 0 or 1 branch
                    // If there is a potentially feasible descendant, then keep descending the tree
                    if let Some(ccons) = fixed.cumulative.get(index) {
                        if accumulator
                            .iter()
                            .zip(ccons)
                            .all(|(&a, &b)| a + b >= T::zero())
                        {
                            index += 1;
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
        (objective, count, solution)
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
