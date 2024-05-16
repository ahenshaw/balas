use crate::Balas;
use crate::NodeState;
use num::Bounded;
use std::{fmt::Display, io::Write, ops::Neg};

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
    pub fn solve_recursively(&mut self) {
        // Initialize the constraint accumulator with the negation of the b vector (the
        // right-hand side of the constraints).  This way, we can just compare against 0
        // later on.
        let accumulator: Vec<T> = self.fixed.rhs.iter().map(|&a| -a).collect();
        // let vars = BitVec::from_elem(self.coefficients.len(), false);
        let num_vars = self.fixed.coefficients.len();
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
            let cons = &self.fixed.constraints[index];

            // Update the current value of the objective
            objective += &self.fixed.coefficients[index];

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
        let Some(ccons) = self.fixed.cumulative.get(index) else {
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
}
