"""Module M1: ABAQUS .inp subset parser.

Lexical and syntactic parsing of the supported keyword subset defined in
doc 04 section 1. Unsupported keywords are skipped with a logged warning.
Raises InputError (with line number) on malformed input.
"""

import logging
from dataclasses import dataclass, field
from typing import Optional

from .errors import InputError

logger = logging.getLogger("fem.parser")

SUPPORTED_KEYWORDS = frozenset(
    {
        "NODE",
        "ELEMENT",
        "NSET",
        "ELSET",
        "MATERIAL",
        "ELASTIC",
        "SOLID SECTION",
        "BOUNDARY",
        "CLOAD",
        "DLOAD",
        "STEP",
        "STATIC",
        "END STEP",
    }
)

UNSUPPORTED_KEYWORDS = frozenset(
    {
        "RESTART",
        "OUTPUT",
        "EL PRINT",
        "NODE PRINT",
        "EL FILE",
        "NODE FILE",
        "SURFACE",
        "TIE",
        "CONTACT",
        "EQUATION",
        "MPC",
        "KINEMATIC",
        "SHELL SECTION",
        "BEAM SECTION",
        "AMPLITUDE",
        "TEMPERATURE",
    }
)


@dataclass
class RawNode:
    id: int
    x: float
    y: float
    z: float = 0.0


@dataclass
class RawElement:
    id: int
    type_name: str
    node_ids: list[int]


@dataclass
class RawMaterial:
    name: str
    E: float
    nu: float


@dataclass
class RawSection:
    elset_name: str
    material_name: str
    thickness: float


@dataclass
class RawBoundary:
    node_id: int
    dof_first: int
    dof_last: int
    value: float
    node_set: Optional[str] = None


@dataclass
class RawPointLoad:
    node_id: int
    dof: int
    magnitude: float
    node_set: Optional[str] = None


@dataclass
class RawPressureLoad:
    elem_id: int
    face_label: str
    magnitude: float


@dataclass
class RawModel:
    nodes: dict[int, RawNode] = field(default_factory=dict)
    elements: dict[int, RawElement] = field(default_factory=dict)
    materials: dict[str, RawMaterial] = field(default_factory=dict)
    sections: list[RawSection] = field(default_factory=list)
    boundaries: list[RawBoundary] = field(default_factory=list)
    point_loads: list[RawPointLoad] = field(default_factory=list)
    pressure_loads: list[RawPressureLoad] = field(default_factory=list)
    nsets: dict[str, list[int]] = field(default_factory=dict)
    elsets: dict[str, list[int]] = field(default_factory=dict)


def _parse_float(text: str, line_no: int) -> float:
    text = text.strip()
    try:
        return float(text)
    except ValueError:
        raise InputError(f"Malformed float at line {line_no}: {text!r}")


def _parse_int(text: str, line_no: int) -> int:
    text = text.strip()
    try:
        value = int(text)
    except ValueError:
        try:
            value = int(float(text))
        except ValueError:
            raise InputError(f"Malformed integer at line {line_no}: {text!r}")
    if value <= 0:
        raise InputError(f"Identifier must be a positive integer at line {line_no}: {text!r}")
    return value


def _split_keyword(line: str) -> tuple[str, dict[str, Optional[str]]]:
    body = line[1:].strip()
    parts = [p.strip() for p in body.split(",")]
    keyword_name = parts[0].upper()
    params: dict[str, Optional[str]] = {}
    for p in parts[1:]:
        if not p:
            continue
        if "=" in p:
            key, _, value = p.partition("=")
            params[key.strip().upper()] = value.strip()
        else:
            params[p.upper()] = None
    return keyword_name, params


def parse_inp(text: str) -> RawModel:
    model = RawModel()
    current_keyword: Optional[str] = None
    current_params: dict[str, Optional[str]] = {}
    current_material: Optional[str] = None

    lines = text.splitlines()
    line_no = 0
    while line_no < len(lines):
        line_no += 1
        raw = lines[line_no - 1]
        stripped = raw.strip()
        if not stripped:
            continue
        if stripped.startswith("**"):
            continue

        if stripped.startswith("*"):
            keyword_name, params = _split_keyword(stripped)
            current_keyword = keyword_name
            current_params = params
            if keyword_name in UNSUPPORTED_KEYWORDS:
                logger.warning(
                    "[parser] Ignoring unsupported keyword %s (data block skipped)", keyword_name
                )
                current_keyword = "__skip__"
            elif keyword_name not in SUPPORTED_KEYWORDS:
                logger.warning(
                    "[parser] Ignoring unknown keyword %s (data block skipped)", keyword_name
                )
                current_keyword = "__skip__"
            elif keyword_name == "MATERIAL":
                name = params.get("NAME")
                if not name:
                    raise InputError(f"*MATERIAL requires NAME= at line {line_no}")
                current_material = name
                model.materials[name] = RawMaterial(name=name, E=0.0, nu=0.0)
            elif keyword_name == "STEP":
                if params.get("NLGEOM") not in (None, "NO"):
                    raise InputError(f"*STEP NLGEOM must be NO at line {line_no}")
            elif keyword_name == "SOLID SECTION":
                elset_name = params.get("ELSET")
                material_name = params.get("MATERIAL")
                if not elset_name or not material_name:
                    raise InputError(
                        f"*SOLID SECTION requires ELSET= and MATERIAL= at line {line_no}"
                    )
                thickness = 1.0
                while line_no < len(lines):
                    next_raw = lines[line_no].strip()
                    if not next_raw or next_raw.startswith("**"):
                        line_no += 1
                        continue
                    if next_raw.startswith("*"):
                        break
                    fields = [f.strip() for f in next_raw.split(",")]
                    if fields and fields[0]:
                        thickness = _parse_float(fields[0], line_no + 1)
                    line_no += 1
                    break
                model.sections.append(
                    RawSection(
                        elset_name=elset_name,
                        material_name=material_name,
                        thickness=thickness,
                    )
                )
                current_keyword = "__skip__"
            continue

        if current_keyword is None:
            raise InputError(f"Data line without preceding keyword at line {line_no}")

        fields = [f.strip() for f in stripped.split(",")]

        if current_keyword == "__skip__":
            continue

        if current_keyword == "NODE":
            if len(fields) < 3:
                raise InputError(f"*NODE entry needs id,x,y at line {line_no}")
            node_id = _parse_int(fields[0], line_no)
            if node_id in model.nodes:
                raise InputError(f"Duplicate node id {node_id} at line {line_no}")
            x = _parse_float(fields[1], line_no)
            y = _parse_float(fields[2], line_no)
            z = 0.0
            if len(fields) >= 4 and fields[3]:
                z = _parse_float(fields[3], line_no)
            model.nodes[node_id] = RawNode(id=node_id, x=x, y=y, z=z)
            nset_name = current_params.get("NSET")
            if nset_name:
                model.nsets.setdefault(nset_name, []).append(node_id)

        elif current_keyword in ("NSET", "ELSET"):
            is_nset = current_keyword == "NSET"
            set_name = current_params.get("NSET") if is_nset else current_params.get("ELSET")
            if not set_name:
                param = "NSET" if is_nset else "ELSET"
                raise InputError(f"*{current_keyword} requires {param}= at line {line_no}")
            target = model.nsets if is_nset else model.elsets
            for token in fields:
                if not token:
                    continue
                if token.upper() in ("GENERATE", "UNSORTED"):
                    raise InputError(
                        f"*{current_keyword} GENERATE form is not supported "
                        f"at line {line_no}"
                    )
                target.setdefault(set_name, []).append(_parse_int(token, line_no))

        elif current_keyword == "ELEMENT":
            if "TYPE" not in current_params:
                raise InputError(f"*ELEMENT requires TYPE= at line {line_no}")
            elem_id = _parse_int(fields[0], line_no)
            if elem_id in model.elements:
                raise InputError(f"Duplicate element id {elem_id} at line {line_no}")
            node_ids = [_parse_int(f, line_no) for f in fields[1:] if f]
            if not node_ids:
                raise InputError(f"*ELEMENT entry needs at least one node at line {line_no}")
            model.elements[elem_id] = RawElement(
                id=elem_id, type_name=current_params["TYPE"].upper(), node_ids=node_ids
            )
            elset_name = current_params.get("ELSET")
            if elset_name:
                model.elsets.setdefault(elset_name, []).append(elem_id)

        elif current_keyword == "MATERIAL":
            pass

        elif current_keyword == "ELASTIC":
            if current_material is None:
                raise InputError(f"*ELASTIC outside *MATERIAL at line {line_no}")
            if len(fields) < 2:
                raise InputError(f"*ELASTIC entry needs E,nu at line {line_no}")
            E = _parse_float(fields[0], line_no)
            nu = _parse_float(fields[1], line_no)
            model.materials[current_material] = RawMaterial(name=current_material, E=E, nu=nu)

        elif current_keyword == "BOUNDARY":
            if len(fields) < 3:
                raise InputError(f"*BOUNDARY entry needs node,dof_first,dof_last at line {line_no}")
            node_token = fields[0].strip()
            if _is_int(node_token):
                node_id = _parse_int(node_token, line_no)
                node_set = None
            else:
                node_id = 0
                node_set = node_token
            dof_first = _parse_int(fields[1], line_no)
            dof_last = _parse_int(fields[2], line_no)
            if not (1 <= dof_first <= dof_last <= 3):
                raise InputError(
                    f"*BOUNDARY invalid DOF range {dof_first}..{dof_last} at line {line_no}"
                )
            value = 0.0
            if len(fields) >= 4 and fields[3]:
                value = _parse_float(fields[3], line_no)
            model.boundaries.append(
                RawBoundary(
                    node_id=node_id,
                    dof_first=dof_first,
                    dof_last=dof_last,
                    value=value,
                    node_set=node_set,
                )
            )

        elif current_keyword == "CLOAD":
            if len(fields) < 3:
                raise InputError(f"*CLOAD entry needs node,dof,magnitude at line {line_no}")
            node_token = fields[0].strip()
            if _is_int(node_token):
                node_id = _parse_int(node_token, line_no)
                node_set = None
            else:
                node_id = 0
                node_set = node_token
            dof = _parse_int(fields[1], line_no)
            if not 1 <= dof <= 3:
                raise InputError(f"*CLOAD invalid dof {dof} at line {line_no}")
            magnitude = _parse_float(fields[2], line_no)
            model.point_loads.append(
                RawPointLoad(node_id=node_id, dof=dof, magnitude=magnitude, node_set=node_set)
            )

        elif current_keyword == "DLOAD":
            if len(fields) < 3:
                raise InputError(f"*DLOAD entry needs element,face,magnitude at line {line_no}")
            elem_id = _parse_int(fields[0], line_no)
            face_label = fields[1].upper()
            if not (face_label.startswith("P") and face_label[1:].isdigit()):
                raise InputError(f"*DLOAD invalid face label {face_label!r} at line {line_no}")
            magnitude = _parse_float(fields[2], line_no)
            model.pressure_loads.append(
                RawPressureLoad(elem_id=elem_id, face_label=face_label, magnitude=magnitude)
            )

        elif current_keyword in ("STEP", "STATIC", "END STEP"):
            continue

        else:
            logger.warning("[parser] Ignoring unexpected data under %s at line %d",
                           current_keyword, line_no)

    return model


def load_inp(path: str) -> RawModel:
    with open(path, "r", encoding="utf-8") as f:
        text = f.read()
    return parse_inp(text)


def _is_int(token: str) -> bool:
    token = token.strip()
    try:
        int(token)
        return True
    except ValueError:
        return False