use thiserror::Error;

#[derive(Debug, Error)]
pub enum LpErrors {
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

    #[error("failed to read the LP file")]
    FileReadError(#[source] std::io::Error),

    #[error("failed to parse the LP file")]
    LPParseError(#[source] anyhow::Error),
}
