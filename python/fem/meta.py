"""ELEMENT_META cache (invariant I7).

The canonical source of element metadata is the C++ ElementRegistry,
queried through femcore.element_meta(name). This module is an ephemeral
per-name cache and a test-injectable indirection so the Python layer never
hardcodes element metadata (invariant I2).
"""

from typing import Callable, Optional

from .errors import ModelError

ELEMENT_META: dict[str, dict] = {}

_provider: Optional[Callable[[str], dict]] = None


def _default_provider(name: str) -> dict:
    try:
        import femcore
    except ImportError:
        raise ModelError(
            f"Cannot validate element '{name}': femcore extension is not built/unavailable. "
            "Run `pip install . -v` first."
        )
    return femcore.element_meta(name)


def configure(provider: Optional[Callable[[str], dict]]) -> None:
    global _provider
    _provider = provider


def element_meta(name: str) -> dict:
    if name in ELEMENT_META:
        return ELEMENT_META[name]
    provider = _provider if _provider is not None else _default_provider
    meta = provider(name)
    ELEMENT_META[name] = dict(meta)
    return meta


def registered_elements() -> list[str]:
    try:
        import femcore
    except ImportError:
        return []
    return list(femcore.registered_elements())
