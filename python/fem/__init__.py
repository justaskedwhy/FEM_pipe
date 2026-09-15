"""fem: Python orchestration layer of the hybrid Python-C++ FE solver.

Modules: M1 parser, M2 model builder, M3 preprocessor, M10 VTK writers,
plus the meta cache (invariant I7) and the CLI (Layer 5).
The numerical core (M4-M9) lives in the compiled ``femcore`` extension.
"""

from .errors import FEMError, InputError, ModelError
from .parser import RawModel, load_inp, parse_inp
from .model_builder import FEMModel, build_model
from .preprocessor import PreprocessedInput, preprocess

__all__ = [
    "FEMError",
    "InputError",
    "ModelError",
    "RawModel",
    "load_inp",
    "parse_inp",
    "FEMModel",
    "build_model",
    "PreprocessedInput",
    "preprocess",
]

__version__ = "0.1.0"


def solve_cli(inp_path: str, out_path: str, verbose: bool = False) -> int:
    from .cli import main

    return main(
        ["--in", inp_path, "--out", out_path] + (["--verbose"] if verbose else [])
    )
