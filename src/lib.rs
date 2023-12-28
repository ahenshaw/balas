use num::Bounded;
use std::{fmt::Display, ops::Neg};

type Array<T> = Vec<Vec<T>>;

pub struct Balas<T> {
    coefficients: Vec<T>,
    c: Array<T>,
    b: Vec<T>,
    cumulative: Array<T>,
    best: T,
    count: usize,
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
    pub fn new(coeff: &[T], c: &Array<T>, b: &[T]) -> Balas<T> {
        let cumulative = Self::make_cumulative(c);
        Balas {
            coefficients: coeff.to_vec(),
            c: c.clone(),
            b: b.to_vec(),
            cumulative,
            best: T::max_value(),
            count: 0,
        }
    }

    pub fn solve(&mut self) {
        // Initialize the constraint accumulator with the negation of the b vector (the
        // right-hand side of the constraints).  This way, we can just compare against 0
        // later on.
        let accumulator: Vec<T> = self.b.clone().into_iter().map(|a| -a).collect();

        self.node(0, 0, &accumulator, &T::zero());
        self.node(1, 0, &accumulator, &T::zero());
    }

    fn node(&mut self, branch: i32, index: usize, accumulator: &[T], objective: &T) {
        // Get aliases for the current "column" of the cumulative constraints and
        // the constraints.
        let ccons = &self.cumulative[index];
        let cons = &self.c[index];
        let coeff = &self.coefficients[index];
        let mut objective = *objective;
        let mut accumulator = accumulator.to_owned();

        if branch == 1 {
            // Update the current value of the objective
            objective += coeff;
            // If we're not better than the current best objective, then we can prune this entire branch,
            if objective >= self.best {
                return;
            }

            // Check if constraints satisfied, while updating the accumulator.
            // We do not have to check the 0 branch, as the accumulator is not changed there.
            // If all of constraints are satisfied, then we are fathomed and we can't do any better.
            accumulator.iter_mut().zip(cons).for_each(|(a, b)| *a += b);
            if accumulator.iter().all(|a| *a >= T::zero()) {
                // println!("New best objective: {}", z);
                self.best = objective;
                return;
            }
        }
        // If there is a potentially feasible descendant, then spawn 0 and 1 child nodes
        if accumulator
            .iter()
            .zip(ccons)
            .all(|(&a, &b)| a + b >= T::zero())
        {
            self.node(0, index + 1, &accumulator, &objective);
            self.node(1, index + 1, &accumulator, &objective);
        }
    }

    pub fn report(&self) {
        println!("Minimum value: {:?}", self.best);
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
