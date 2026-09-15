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
        "PART",
        "END PART",
        "ASSEMBLY",
        "INSTANCE",
        "END INSTANCE",
        "END ASSEMBLY",
        "SURFACE",
        "MATERIAL",
        "ELASTIC",
        "SOLID SECTION",
        "BOUNDARY",
        "CLOAD",
        "DLOAD",
        "DSLOAD",
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
        "HEADING",
        "PREPRINT",
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
    surface_name: Optional[str] = None


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
    surfaces: dict[str, list[tuple[str, str]]] = field(default_factory=dict)


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
        raise InputError(
            f"Identifier must be a positive integer at line {line_no}: {text!r}"
        )
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
    current_part: Optional[str] = None
    last_instance: Optional[str] = None
    current_surface_name: Optional[str] = None

    def store_name(name: str) -> str:
        instance = current_params.get("INSTANCE")
        if instance:
            return f"instance:{instance}:{name}"
        if current_part:
            return f"part:{current_part}:{name}"
        return name

    def ref_name(name: str, target: dict) -> str:
        if current_part:
            return f"part:{current_part}:{name}"
        if last_instance:
            key = f"instance:{last_instance}:{name}"
            if key in target:
                return key
        return name

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
                    "[parser] Ignoring unsupported keyword %s (data block skipped)",
                    keyword_name,
                )
                current_keyword = "__skip__"
            elif keyword_name not in SUPPORTED_KEYWORDS:
                logger.warning(
                    "[parser] Ignoring unknown keyword %s (data block skipped)",
                    keyword_name,
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
            elif keyword_name in ("PART", "INSTANCE"):
                instance = params.get("NAME")
                if keyword_name == "INSTANCE":
                    last_instance = instance
                    current_keyword = "__skip__"
                else:
                    current_part = instance
            elif keyword_name in (
                "END PART",
                "END INSTANCE",
                "END ASSEMBLY",
                "ASSEMBLY",
            ):
                if keyword_name == "END PART":
                    current_part = None
                current_keyword = "__skip__"
            elif keyword_name == "SURFACE":
                current_surface_name = params.get("NAME")
                surface_type = params.get("TYPE")
                if not current_surface_name:
                    raise InputError(f"*SURFACE requires NAME= at line {line_no}")
                if surface_type and surface_type.upper() != "ELEMENT":
                    raise InputError(f"*SURFACE type must be ELEMENT at line {line_no}")
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
                        elset_name=ref_name(elset_name, model.elsets),
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
                model.nsets.setdefault(store_name(nset_name), []).append(node_id)

        elif current_keyword in ("NSET", "ELSET"):
            is_nset = current_keyword == "NSET"
            set_name = (
                current_params.get("NSET") if is_nset else current_params.get("ELSET")
            )
            if not set_name:
                param = "NSET" if is_nset else "ELSET"
                raise InputError(
                    f"*{current_keyword} requires {param}= at line {line_no}"
                )
            target = model.nsets if is_nset else model.elsets
            key = store_name(set_name)
            utils = list(target.setdefault(key, []))
            if "GENERATE" in current_params:
                if len(fields) < 2:
                    raise InputError(
                        f"*{current_keyword} GENERATE needs first,last[,step] "
                        f"at line {line_no}"
                    )
                first = _parse_int(fields[0], line_no)
                last = _parse_int(fields[1], line_no)
                step = 1
                if len(fields) >= 3 and fields[2]:
                    step = _parse_int(fields[2], line_no)
                if step <= 0:
                    raise InputError(
                        f"*{current_keyword} GENERATE step must be > 0 at line {line_no}"
                    )
                for value in range(first, last + 1, step):
                    utils.append(value)
                current_keyword = "__skip__"
            else:
                for token in fields:
                    if not token:
                        continue
                    token = token.strip()
                    if token in ("UNSORTED", "INTERNAL"):
                        continue
                    if token.isdigit():
                        utils.append(_parse_int(token, line_no))
                    else:
                        if is_nset:
                            if token not in model.nsets:
                                raise InputError(
                                    f"*NSET references undefined set '{token}' "
                                    f"at line {line_no}"
                                )
                            utils.extend(model.nsets[token])
                        else:
                            if token not in model.elsets:
                                raise InputError(
                                    f"*ELSET references undefined set '{token}' "
                                    f"at line {line_no}"
                                )
                            utils.extend(model.elsets[token])
            target[key] = list(dict.fromkeys(utils))

        elif current_keyword == "ELEMENT":
            if "TYPE" not in current_params:
                raise InputError(f"*ELEMENT requires TYPE= at line {line_no}")
            elem_id = _parse_int(fields[0], line_no)
            if elem_id in model.elements:
                raise InputError(f"Duplicate element id {elem_id} at line {line_no}")
            node_ids = [_parse_int(f, line_no) for f in fields[1:] if f]
            if not node_ids:
                raise InputError(
                    f"*ELEMENT entry needs at least one node at line {line_no}"
                )
            model.elements[elem_id] = RawElement(
                id=elem_id, type_name=current_params["TYPE"].upper(), node_ids=node_ids
            )
            elset_name = current_params.get("ELSET")
            if elset_name:
                key = store_name(elset_name)
                model.elsets.setdefault(key, []).append(elem_id)

        elif current_keyword == "MATERIAL":
            pass

        elif current_keyword == "ELASTIC":
            if current_material is None:
                raise InputError(f"*ELASTIC outside *MATERIAL at line {line_no}")
            if len(fields) < 2:
                raise InputError(f"*ELASTIC entry needs E,nu at line {line_no}")
            E = _parse_float(fields[0], line_no)
            nu = _parse_float(fields[1], line_no)
            model.materials[current_material] = RawMaterial(
                name=current_material, E=E, nu=nu
            )

        elif current_keyword == "BOUNDARY":
            if len(fields) < 3:
                raise InputError(
                    f"*BOUNDARY entry needs node,dof_first,dof_last at line {line_no}"
                )
            node_token = fields[0].strip()
            if _is_int(node_token):
                node_id = _parse_int(node_token, line_no)
                node_set = None
            else:
                node_id = 0
                node_set = ref_name(node_token, model.nsets)
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
                raise InputError(
                    f"*CLOAD entry needs node,dof,magnitude at line {line_no}"
                )
            node_token = fields[0].strip()
            if _is_int(node_token):
                node_id = _parse_int(node_token, line_no)
                node_set = None
            else:
                node_id = 0
                node_set = ref_name(node_token, model.nsets)
            dof = _parse_int(fields[1], line_no)
            if not 1 <= dof <= 3:
                raise InputError(f"*CLOAD invalid dof {dof} at line {line_no}")
            magnitude = _parse_float(fields[2], line_no)
            model.point_loads.append(
                RawPointLoad(
                    node_id=node_id, dof=dof, magnitude=magnitude, node_set=node_set
                )
            )

        elif current_keyword == "DLOAD":
            if len(fields) < 3:
                raise InputError(
                    f"*DLOAD entry needs element,face,magnitude at line {line_no}"
                )
            elem_id = _parse_int(fields[0], line_no)
            face_label = fields[1].upper()
            if not (face_label.startswith("P") and face_label[1:].isdigit()):
                raise InputError(
                    f"*DLOAD invalid face label {face_label!r} at line {line_no}"
                )
            magnitude = _parse_float(fields[2], line_no)
            model.pressure_loads.append(
                RawPressureLoad(
                    elem_id=elem_id, face_label=face_label, magnitude=magnitude
                )
            )

        elif current_keyword == "SURFACE":
            if len(fields) < 2:
                raise InputError(f"*SURFACE entry needs set,face at line {line_no}")
            set_token = fields[0].strip()
            face_label = fields[1].upper()
            if not (face_label.startswith("S") and face_label[1:].isdigit()):
                raise InputError(
                    f"*SURFACE invalid face label {face_label!r} at line {line_no}"
                )
            if current_surface_name is None:
                raise InputError(f"*SURFACE data without NAME= at line {line_no}")
            resolved_set = (
                () if set_token.isdigit() else ref_name(set_token, model.elsets)
            )
            if set_token.isdigit():
                resolved_set = set_token
            model.surfaces.setdefault(current_surface_name, []).append(
                (resolved_set, face_label)
            )

        elif current_keyword == "DSLOAD":
            if len(fields) < 3:
                raise InputError(
                    f"*DSLOAD entry needs surface,P,magnitude at line {line_no}"
                )
            surface_name = fields[0].strip()
            label = fields[1].upper()
            if label != "P":
                raise InputError(
                    f"*DSLOAD invalid load label {label!r} at line {line_no}"
                )
            magnitude = _parse_float(fields[2], line_no)
            model.pressure_loads.append(
                RawPressureLoad(
                    elem_id=0,
                    face_label="",
                    magnitude=magnitude,
                    surface_name=surface_name,
                )
            )

        elif current_keyword in ("STEP", "STATIC", "END STEP"):
            continue

        else:
            logger.warning(
                "[parser] Ignoring unexpected data under %s at line %d",
                current_keyword,
                line_no,
            )

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
