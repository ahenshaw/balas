use thiserror::Error;

#[derive(Debug, Error)]
pub enum LpErrors {
    #[error("All variables must be Binary")]
    VarNotBinary,

    #[error("No variables found")]
    NoVars,

    #[error("Expected objective")]
    NoObjective,

    #[error("Can only handle Standard constraints")]
    UnexpectedConstraintType,

    #[error("failed to read the LP file")]
    FileReadError(#[source] std::io::Error),

    #[error("failed to parse the LP file")]
    LPParseError(#[source] anyhow::Error),
}
