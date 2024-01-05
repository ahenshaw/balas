use crate::lp_errors::LpErrors;
use crate::Balas;
use lp_parser_rs::model::coefficient::Coefficient;
use lp_parser_rs::model::constraint::Constraint;
use lp_parser_rs::model::lp_problem::LPProblem;
use lp_parser_rs::model::objective::Objective;
use lp_parser_rs::model::sense::Sense;
use lp_parser_rs::model::variable::VariableType;
use lp_parser_rs::parse::parse_lp_file;
use std::collections::HashMap;
use std::fs;
use std::path::Path;

type Constraints = HashMap<String, Constraint>;

impl Balas<f64> {
    pub fn from_lp(lp_path: &Path) -> Result<Balas<f64>, LpErrors> {
        let code = fs::read_to_string(lp_path).map_err(LpErrors::FileReadError)?;
        let lp = parse_lp_file(&code).map_err(LpErrors::LPParseError)?;

        let lp = normalize_for_balas(&lp)?;

        // dbg!(&lp);

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

            coefficients = obj.iter().map(|c| c.coefficient).collect()
        } else {
            return Err(LpErrors::NoObjective);
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
                _ => return Err(LpErrors::UnexpectedConstraintType),
            }
        }
        // dbg!(&coefficients);
        // dbg!(&constraints);
        // dbg!(&rhs);
        // dbg!(&vars);
        Ok(Balas::new(&coefficients, &constraints, &rhs, &vars))
    }
}

/// The Balas algorithm requires that:
/// - the problem sense must be "minimize".  A "maximize" sense
///     will be converted by negating the objective coefficients.
/// - all constraints be of the >= sense.  So, this function
///     will convert <= sense constraints to >= by negating
///     the coefficients and the rhs. Equality constraints
///     need to have a negated constraint added.
/// - all objective coefficients must be positive.  Negative
///     coefficients will be converted by replacing "x"
///     with "y = 1 - x"
///
fn normalize_for_balas(lp: &LPProblem) -> Result<LPProblem, LpErrors> {
    let problem_name = format!("{}_balas", lp.problem_name);
    let objective = create_min_objective(lp)?;
    let constraints = create_ge_constraints(lp)?;
    let (objective, constraints) = fix_neg_variables(&objective, &constraints);

    // copy variables while making sure they all are binary
    let mut variables = HashMap::new();
    for (s, vtype) in &lp.variables {
        if *vtype != VariableType::Binary {
            return Err(LpErrors::VarNotBinary);
        }
        variables.insert(s.clone(), VariableType::Binary);
    }

    Ok(LPProblem {
        problem_name,
        problem_sense: Sense::Minimize,
        variables,
        objectives: vec![objective],
        constraints,
    })
}

fn fix_neg_variables(objective: &Objective, constraints: &Constraints) -> (Objective, Constraints) {
    let mut to_change = Vec::<&str>::new();

    let mut coeff: Vec<Coefficient> = vec![];
    for coeff_var in &objective.coefficients {
        let mut c = coeff_var.coefficient;
        if coeff_var.coefficient < 0.0 {
            c = -c;
            to_change.push(&coeff_var.var_name);
        }
        coeff.push(Coefficient {
            var_name: coeff_var.var_name.clone(),
            coefficient: c,
        });
    }
    let objective = Objective {
        name: objective.name.clone(),
        coefficients: coeff,
    };

    let mut new_constraints = Constraints::new();
    for (label, constraint) in constraints {
        if let Constraint::Standard {
            name,
            coefficients,
            sense,
            rhs,
        } = constraint
        {
            let mut rhs = *rhs;
            let mut coefficients: Vec<Coefficient> = coefficients
                .iter()
                .map(|c| Coefficient {
                    var_name: c.var_name.clone(),
                    coefficient: c.coefficient,
                })
                .collect();

            for &var_name in &to_change {
                coefficients.iter_mut().for_each(|c| {
                    if c.var_name == var_name {
                        rhs -= c.coefficient;
                        c.coefficient = -c.coefficient;
                    }
                });
            }
            let new_constraint = Constraint::Standard {
                name: name.to_owned(),
                coefficients,
                sense: sense.to_owned(),
                rhs,
            };
            new_constraints.insert(label.clone(), new_constraint);
        }
    }
    (objective, new_constraints)
}

fn create_ge_constraints(lp: &LPProblem) -> Result<Constraints, LpErrors> {
    // make all constraints be of the >= sense
    let mut additional: Vec<(String, Constraint)> = vec![];
    let mut constraints: Constraints = lp
        .constraints
        .iter()
        .map(|(label, constraint)| {
            match constraint {
                Constraint::Standard {
                    name,
                    coefficients,
                    sense,
                    rhs,
                } => {
                    let mut my_sense = sense.to_owned();
                    let mut my_coefficients: Vec<Coefficient> = coefficients
                        .iter()
                        .map(|c| Coefficient {
                            var_name: c.var_name.clone(),
                            coefficient: c.coefficient,
                        })
                        .collect();
                    let mut my_rhs = *rhs;
                    match sense.as_str() {
                        "<=" => {
                            my_sense = ">=".to_owned();
                            my_rhs = -rhs;
                            my_coefficients
                                .iter_mut()
                                .for_each(|c| c.coefficient = -c.coefficient);
                        }
                        "=" => {
                            // change over "=" to ">="
                            my_sense = ">=".to_owned();

                            // create negated copy of constraint
                            let new_coeff: Vec<Coefficient> = my_coefficients
                                .iter()
                                .map(|c| Coefficient {
                                    var_name: c.var_name.clone(),
                                    coefficient: -c.coefficient,
                                })
                                .collect();
                            additional.push((
                                format!("{}_balas", label),
                                Constraint::Standard {
                                    name: format!("{}_balas", name),
                                    coefficients: new_coeff,
                                    sense: ">=".to_owned(),
                                    rhs: -rhs,
                                },
                            ));
                        }
                        _ => {}
                    }
                    (
                        label.to_owned(),
                        Constraint::Standard {
                            name: name.to_owned(),
                            coefficients: my_coefficients,
                            sense: my_sense,
                            rhs: my_rhs,
                        },
                    )
                }
                _ => unimplemented!(),
            }
        })
        .collect();
    for (label, constraint) in additional {
        constraints.insert(label, constraint);
    }
    Ok(constraints)
}

fn create_min_objective(lp: &LPProblem) -> Result<Objective, LpErrors> {
    let negate = lp.problem_sense == Sense::Maximize;
    if let Some(objective) = lp.objectives.first() {
        let outer: Vec<Coefficient> = objective
            .coefficients
            .iter()
            .map(|c| {
                let coeff = if negate {
                    -c.coefficient
                } else {
                    c.coefficient
                };
                Coefficient {
                    var_name: c.var_name.clone(),
                    coefficient: coeff,
                }
            })
            .collect();
        Ok(Objective {
            name: objective.name.clone(),
            coefficients: outer,
        })
    } else {
        Err(LpErrors::NoObjective)
    }
}
