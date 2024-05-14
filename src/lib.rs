mod lp_errors;
mod lp_reader;

use num::Bounded;
use serde::{Deserialize, Serialize};
use std::{fmt::Display, io::Write, ops::Neg};

type Array<T> = Vec<Vec<T>>;

#[derive(Deserialize, Serialize)]
pub struct Balas<T> {
    pub coefficients: Vec<T>,
    pub constraints: Array<T>,
    pub rhs: Vec<T>,
    cumulative: Array<T>,
    pub best: T,
    #[serde(skip_serializing, skip_deserializing)]
    pub solution: Vec<u8>,
    pub count: usize,
    vars: Vec<String>,
    pub recording: Vec<Record>,
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
        Balas {
            coefficients: coeff.to_vec(),
            constraints: constraints.clone(),
            rhs: b.to_vec(),
            cumulative,
            best: T::max_value(),
            solution: Vec::new(),
            count: 0,
            vars: vars.to_owned(),
            recording: vec![],
        }
    }
    pub fn reset(&mut self) {
        self.count = 0;
        self.best = T::max_value();
        self.solution = Vec::new();
    }


    pub fn solve_non_recursive(&mut self) {

        let num_vars = self.coefficients.len();
        let mut vars: Vec<u8> = vec![0;num_vars];
        let mut branch = 0u8;
        let mut index: usize = 0;
        let mut objective = T::zero();
        let mut state = Flow::Normal;

        // Initialize the constraint accumulator with the negation of the b vector (the
        // right-hand side of the constraints).  This way, we can just compare against 0
        // later on.
        let mut accumulator: Vec<T> = self.rhs.iter().map(|&b| -b).collect();

        while index < num_vars || branch == 0 {

            // Alias the current column of the constraints and grab the coefficients value
            let cons = &self.constraints[index];
            let coeff = &self.coefficients[index];

            match state {
                Flow::Backtrack => {
                    if vars[index] == 1 {
                        // we have to reverse what we did before we leave
                        accumulator.iter_mut().zip(cons).for_each(|(a, b)| *a -= b);
                        objective -= coeff;
                        index -= 1;
                    } else {
                        state = Flow::Normal;
                        branch = 0;
                    }
                },
                Flow::Normal => {
                    self.count += 1;
                    vars[index] = branch;

                    if branch == 1 {
                        // Update the accumulator.  This only needs to be done in the ones branch
                        accumulator.iter_mut().zip(cons).for_each(|(a, b)| *a += b);

                        // Update the current value of the objective
                        objective += coeff;

                        // If we're already not better than the current best objective, then
                        // we can prune this entire branch.
                        if objective >= self.best {
                            state = Flow::Backtrack;
                        } else {
                            // Check if constraints satisfied, while updating the accumulator.
                            // We do not have to check the 0 branch, as the accumulator is not changed there.
                            // If all of constraints are satisfied, then we are fathomed and we can't do any better.
                            if accumulator.iter().all(|x| *x >= T::zero()) {
                                self.best = objective;
                                // self.solution = vars;
                                state = Flow::Backtrack;
                            }
                        }
                    }
                }
            }
        }
    }


    pub fn solve(&mut self) {
        // Initialize the constraint accumulator with the negation of the b vector (the
        // right-hand side of the constraints).  This way, we can just compare against 0
        // later on.
        let accumulator: Vec<T> = self.rhs.iter().map(|&a| -a).collect();
        // let vars = BitVec::from_elem(self.coefficients.len(), false);
        let num_vars = self.coefficients.len();
        let vars = vec![0u8; num_vars];
        self.record("", NodeState::Active);
        self.record("", NodeState::Visited);

        // self.node(0, 0, &accumulator, &T::zero(), &vars, "0".to_string());
        // self.node(1, 0, &accumulator, &T::zero(), &vars, "1".to_string());
        self.node(0, 0, &accumulator, &T::zero(), &vars);
        self.node(1, 0, &accumulator, &T::zero(), &vars);
    }

    fn node(
        &mut self,
        branch: u8,
        index: usize,
        accumulator: &[T],
        objective: &T,
        vars: &Vec<u8>,
        // label: String,
    ) {
        // self.record(&label, NodeState::Active);
        let mut objective = *objective;
        let mut vars = vars.to_owned();
        let mut accumulator = accumulator.to_owned();
        // Alias the current column of the cumulative constraints
        // let ccons = &self.cumulative[index];

        self.count += 1;
        // println!("count:{}  branch:{branch}  index:{index}  objective:{objective}  accumulator:{accumulator:?}", self.count);

        if branch == 1 {
            vars[index] = 1;
            // Alias the current column of the constraints
            let cons = &self.constraints[index];

            // Update the current value of the objective
            objective += &self.coefficients[index];

            // If we're already not better than the current best objective, then
            // we can prune this entire branch.
            if objective >= self.best {
                // self.record(&label, NodeState::Suboptimal);
                return;
            }

            // Check if constraints satisfied, while updating the accumulator.
            // We do not have to check the 0 branch, as the accumulator is not changed there.
            // If all of constraints are satisfied, then we are fathomed and we can't do any better.
            accumulator.iter_mut().zip(cons).for_each(|(a, b)| *a += b);
            if accumulator.iter().all(|a| *a >= T::zero()) {
                // println!("New best objective: {} {:?}", objective, vars);
                self.best = objective;
                // print!("{objective} ");
                std::io::stdout().flush().unwrap();
                self.solution = vars;
                // self.record(&label, NodeState::Fathomed);
                return;
            }
        }
        // self.record(&label, NodeState::Visited);
        // If there is a potentially feasible descendant, then spawn 0 and 1 child nodes
        let Some(ccons) = self.cumulative.get(index) else {
            // println!("run out of vars with index: {index}");
            // self.record(&label, NodeState::Infeasible);
            return;
        };

        if accumulator
            .iter()
            .zip(ccons)
            .all(|(&a, &b)| a + b >= T::zero())
        {
            self.node(
                0,
                index + 1,
                &accumulator,
                &objective,
                &vars,
                // label.clone() + "0",
            );
            self.node(
                1,
                index + 1,
                &accumulator,
                &objective,
                &vars,
                // label.clone() + "1",
            );
        } else {
            // self.record(&label, NodeState::ImpossibleChildren);
        }
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

enum Flow {
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
