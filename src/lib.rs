mod lp_errors;
mod lp_reader;

use bit_vec::BitVec;
use num::Bounded;
use serde::{Deserialize, Serialize};
use std::{fmt::Display, ops::Neg};

type Array<T> = Vec<Vec<T>>;

#[derive(Deserialize, Serialize)]
pub struct Balas<T> {
    pub coefficients: Vec<T>,
    pub constraints: Array<T>,
    pub rhs: Vec<T>,
    cumulative: Array<T>,
    pub best: T,
    #[serde(skip_serializing, skip_deserializing)]
    pub solution: BitVec,
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
            solution: BitVec::new(),
            count: 0,
            vars: vars.to_owned(),
            recording: vec![],
        }
    }
    pub fn reset(&mut self) {
        self.count = 0;
        self.best = T::max_value();
        self.solution = BitVec::new();
    }

    pub fn solve(&mut self) {
        // Initialize the constraint accumulator with the negation of the b vector (the
        // right-hand side of the constraints).  This way, we can just compare against 0
        // later on.
        let accumulator: Vec<T> = self.rhs.iter().map(|&a| -a).collect();
        let vars = BitVec::from_elem(self.coefficients.len(), false);
        self.record("", NodeState::Active);
        self.record("", NodeState::Visited);

        // self.node(0, 0, &accumulator, &T::zero(), &vars, "0".to_string());
        // self.node(1, 0, &accumulator, &T::zero(), &vars, "1".to_string());
        self.node(false, 0, &accumulator, &T::zero(), &vars);
        self.node(true, 0, &accumulator, &T::zero(), &vars);
    }

    fn node(
        &mut self,
        branch: bool,
        index: usize,
        accumulator_in: &[T],
        objective: &T,
        vars: &BitVec,
        // label: String,
    ) {
        // self.record(&label, NodeState::Active);
        let mut objective = *objective;
        let mut vars = vars.to_owned();
        // let mut accumulator = accumulator.to_owned();
        let accumulator: Vec<T>;

        self.count += 1;

        if branch {
            // Update the current value of the objective
            objective += &self.coefficients[index];

            // If we're already not better than the current best objective, then
            // we can prune this entire branch.
            if objective >= self.best {
                // self.record(&label, NodeState::Suboptimal);
                return;
            }

            vars.set(index, true);

            // Alias the current column of the constraints
            let cons = &self.constraints[index];

            // Update the accumulator
            accumulator = accumulator_in
                .iter()
                .zip(cons)
                .map(|(a, b)| *a + *b)
                .collect();

            // Check if constraints are satisfied.
            // We do not have to check the 0 branch, as the accumulator is not changed there.
            // If all of constraints are satisfied, then we are fathomed and we can't do any better.
            if accumulator.iter().all(|a| *a >= T::zero()) {
                self.best = objective;
                self.solution = vars;
                // self.record(&label, NodeState::Fathomed);
                return;
            }
        } else {
            accumulator = accumulator_in.to_owned();
        }
        // self.record(&label, NodeState::Visited);

        // If there is a potentially feasible descendant, then spawn 0 and 1 child nodes
        let Some(ccons) = self.cumulative.get(index) else {
            // no more nodes to check in this branch
            // self.record(&label, NodeState::Infeasible);
            return;
        };

        if accumulator
            .iter()
            .zip(ccons)
            .all(|(&a, &b)| a + b >= T::zero())
        {
            self.node(
                false,
                index + 1,
                &accumulator,
                &objective,
                &vars,
                // label.clone() + "0",
            );
            self.node(
                true,
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
                println!("  {var}: {}", value as u8);
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
