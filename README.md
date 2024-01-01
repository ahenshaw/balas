Balas Additive algorithm implemented in Rust.  

This is not yet ready for release.

The Balas Additive algorithm solves binary-variable linear programming problems.  In addition to requiring that
all of the variables must be binary, the algorithm places other limitations on the formulation of the problem (which can 
all be worked around):
- the objective must be *minimized*
- the coefficients for the variables in the objective must be positive
- the variables in the objective must be arranged in ascending order (by coefficient)
- the inequality used in all constraints must be >= (greater-than-or-equal)

All of these limitations (except binary variables only) are on the roadmap to be removed.
