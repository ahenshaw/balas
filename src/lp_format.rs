use crate::{lp_errors::LpErrors, Balas};
use lp_parser_rs::{
    model::{constraint::Constraint, sense::Sense, variable::VariableType},
    parse::parse_lp_file,
};
use std::{fs, path::Path};

impl Balas<f64> {
    pub fn from_lp(lp_path: &Path) -> Result<Balas<f64>, LpErrors> {
        let code = fs::read_to_string(lp_path).map_err(LpErrors::FileReadError)?;
        let lp = parse_lp_file(&code).map_err(LpErrors::LPParseError)?;
        dbg!(&lp);

        // LP binary problem normalization checks
        if lp.problem_sense != Sense::Minimize {
            return Err(LpErrors::ProblemSenseNotMinimize);
        }
        if lp.variables.iter().any(|(_, vtype)| *vtype != VariableType::Binary) {
            return Err(LpErrors::VarNotBinary);
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

            index = obj.iter().enumerate().map(|(i, c)| (c.var_name.to_owned(), i)).collect();

            coefficients = obj.iter().map(|c| c.coefficient).collect()
        } else {
            return Err(LpErrors::NoObjective);
        }
        if coefficients.iter().any(|&x| x < 0.0) {
            return Err(LpErrors::NegativeCoefficient);
        }
        let num_vars = lp.variables.len();
        let num_constraints = lp.constraints.len();
        let mut constraints: Vec<_> = (0..num_vars).map(|_| vec![0.0; num_constraints]).collect();
        let mut rhs = vec![];
        for (col, (_, constraint)) in lp.constraints.iter().enumerate() {
            match constraint {
                Constraint::Standard { name: _, coefficients, sense: _, rhs: this_rhs } => {
                    rhs.push(*this_rhs);
                    for coeff in coefficients {
                        if let Some(row) = index.get(&coeff.var_name) {
                            constraints[*row][col] = coeff.coefficient;
                        }
                    }
                }
                _ => return Err(LpErrors::UnexpectedConstraintType),
            }
        }
        Ok(Balas::new(&coefficients, &constraints, &rhs, &vars))
    }
}
