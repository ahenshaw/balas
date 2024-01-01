use crate::Balas;
use lp_parser_rs::model::constraint::Constraint;
use lp_parser_rs::model::lp_problem::LPProblem;
use lp_parser_rs::model::sense::Sense;
use lp_parser_rs::model::variable::VariableType;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum Errors {
    #[error("All variables must be Binary")]
    VarNotBinary,
    #[error("No variables found")]
    NoVars,
    #[error("All constraints must be greater-than-or equal (TO-DO: remove this restriction)")]
    AllSenseGE,
    #[error("Must be a minimization problem (TO-DO: remove this restriction)")]
    ProblemSenseNotMinimize,
    #[error("Expected objective")]
    NoObjective,
    #[error("Coefficients must be positive (TO-DO: remove this restriction)")]
    NegativeCoefficient,
    #[error("Can only handle Standard constraints")]
    UnexpectedConstraintType,
}

impl Balas<f64> {
    pub fn from_lp(lp: &LPProblem) -> Result<Balas<f64>, Errors> {
        // LP binary problem normalization checks
        if lp.problem_sense != Sense::Minimize {
            return Err(Errors::ProblemSenseNotMinimize);
        }
        if lp
            .variables
            .iter()
            .any(|(_, vtype)| *vtype != VariableType::Binary)
        {
            return Err(Errors::VarNotBinary);
        }

        let coefficients: Vec<f64>;
        let vars: Vec<String>;
        let index: std::collections::HashMap<String, usize>;
        if let Some(objective) = lp.objectives.first() {
            // sort the variables by coefficient
            let mut obj: Vec<_> = objective.coefficients.iter().collect();
            obj.sort_by(|a, b| a.coefficient.partial_cmp(&b.coefficient).unwrap());
            vars = obj.iter().map(|v| v.var_name.to_owned()).collect();

            // Create mapping from variable name to constraints column (visually)
            // For the solver, the constraints are transposed (for efficiency), so the index
            // maps to a row index.
            index = obj
                .iter()
                .enumerate()
                .map(|(i, c)| (c.var_name.to_owned(), i))
                .collect();

            coefficients = objective
                .coefficients
                .iter()
                .map(|c| c.coefficient)
                .collect()
        } else {
            return Err(Errors::NoObjective);
        }
        if coefficients.iter().any(|&x| x < 0.0) {
            return Err(Errors::NegativeCoefficient);
        }
        let num_vars = lp.variables.len();
        let num_constraints = lp.constraints.len();
        let mut constraints: Vec<_> = (0..num_vars).map(|_| vec![0.0; num_constraints]).collect();
        let mut rhs = vec![];
        for (col, (_, constraint)) in lp.constraints.iter().enumerate() {
            match constraint {
                Constraint::Standard {
                    name: _,
                    coefficients,
                    sense: _,
                    rhs: this_rhs,
                } => {
                    rhs.push(*this_rhs);
                    for coeff in coefficients {
                        if let Some(row) = index.get(&coeff.var_name) {
                            constraints[*row][col] = coeff.coefficient;
                        }
                    }
                }
                _ => return Err(Errors::UnexpectedConstraintType),
            }
        }
        Ok(Balas::new(&coefficients, &constraints, &rhs, &vars))
    }
}
