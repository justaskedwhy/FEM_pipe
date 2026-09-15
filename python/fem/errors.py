"""Exception classes for the python/fem orchestration layer.

Exit codes follow doc 01 section 7.1: user input = 1, model semantics = 2.
Numerical (3) and solver (4) failures surface from femcore as RuntimeError or
via SolverOutput.converged=False.
"""


class FEMError(Exception):
    exit_code = 5


class InputError(FEMError):
    exit_code = 1


class ModelError(FEMError):
    exit_code = 2
