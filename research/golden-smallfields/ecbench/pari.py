"""Optional independent PARI/GP oracle for canonical ecbench records.

PARI is not needed by the package's standard-library verifier. This module
checks point counts, abstract groups, named-point orders, and all successful
curve arithmetic vectors. Binary polynomial-basis coordinates are interpreted
in their exact recorded fields, and all point comparisons happen inside GP.
No generated dataset is modified.
"""
from __future__ import annotations

import json
import os
import shutil
import subprocess
from typing import Any


def _binary_element(value: int) -> str:
    if value < 0:
        raise ValueError("negative binary field coordinate")
    terms = ["u^"+str(i) for i in range(value.bit_length()) if value >> i & 1]
    return "(0*u"+("+"+"+".join(terms) if terms else "")+")"


def _element(value: str, field: dict[str, Any]) -> str:
    number = int(value, 16)
    cardinality = int(field["cardinality"])
    if not 0 <= number < cardinality:
        raise ValueError("noncanonical coordinate supplied to PARI oracle")
    if field["kind"] == "prime":
        return f"Mod({number},{cardinality})"
    if field["kind"] == "binary":
        return _binary_element(number)
    raise ValueError("unsupported field kind")


def _point(point: Any, field: dict[str, Any]) -> str:
    if point is None:
        return "[0]"
    return "["+_element(point["x"],field)+","+_element(point["y"],field)+"]"


def _field_definition(field: dict[str, Any]) -> str:
    if field["kind"] == "prime":
        cardinality = int(field["cardinality"])
        if int(field["representation"]["modulus"],16) != cardinality:
            raise ValueError("prime field cardinality/modulus mismatch")
        return ""
    if field["kind"] != "binary":
        raise ValueError("unsupported field kind")
    exponents = list(map(int,field["representation"]["polynomial_exponents"]))
    if (not exponents or min(exponents)<0 or len(set(exponents))!=len(exponents)
            or max(exponents)!=int(field["element_bits"])):
        raise ValueError("invalid binary defining-polynomial exponents")
    modulus = sum(1 << exponent for exponent in exponents)
    if modulus != int(field["representation"]["modulus"],16):
        raise ValueError("binary defining-polynomial mask/exponent mismatch")
    if int(field["cardinality"]) != 1 << int(field["element_bits"]):
        raise ValueError("binary field cardinality/degree mismatch")
    polynomial = "+".join("z^"+str(exponent) for exponent in exponents)
    return f"u=ffgen(Mod({polynomial},2),'u);"


def verify_with_pari(fields, curves, vectors, gp_path=None):
    """Return check counts, or raise on the first independent GP mismatch.

    ``fields``, ``curves``, and ``vectors`` are canonical record lists. The
    oracle intentionally leaves field, encoding, and membership vectors to
    the standard-library reference verifier; it independently checks every
    add/double/neg/mul curve answer and every stored point order. GP lookup is
    explicit path, CURVE_REPORT_GP, PATH, then the local session fallback.
    """
    executable = (gp_path or os.environ.get("CURVE_REPORT_GP") or
                  shutil.which("gp") or "/tmp/pari-local/usr/bin/gp")
    field_by_id = {field["id"]:field for field in fields}
    curve_by_id = {curve["id"]:curve for curve in curves}
    if len(field_by_id)!=len(fields) or len(curve_by_id)!=len(curves):
        raise ValueError("duplicate field or curve IDs")
    vector_by_curve = {cid:[] for cid in curve_by_id}
    ignored_vector_count = 0
    for vector in vectors:
        if vector["domain"] != "curve":
            ignored_vector_count+=1
            continue
        if vector["curve_id"] not in vector_by_curve:
            raise ValueError("curve vector refers to unknown curve")
        if vector["operation"] not in {"add","double","neg","mul"}:
            raise ValueError("unsupported curve arithmetic operation")
        if "point" not in vector["expected"]:
            # Rejections are exercised by the canonical reference validator.
            ignored_vector_count+=1
            continue
        vector_by_curve[vector["curve_id"]].append(vector)
    lines = ["checks=0;"]
    counts = {"cardinality_checks":0,"group_structure_checks":0,
              "named_point_order_checks":0,"curve_arithmetic_checks":0,
              "subgroup_generator_order_checks":0}
    def check(expression, label, counter):
        message = json.dumps("FAIL "+label,ensure_ascii=True)
        lines.append(f"if({expression},checks++,print({message});quit(1));")
        counts[counter]+=1
    for curve in curves:
        cid = curve["id"]
        field = field_by_id[curve["field_id"]]
        definition = _field_definition(field)
        if definition:
            lines.append(definition)
        coefficients = [_element(curve["coefficients"][name],field)
                        for name in ("a1","a2","a3","a4","a6")]
        lines.append("E=ellinit(["+",".join(coefficients)+"]);")
        order = int(curve["invariants"]["order"])
        d1,d2 = map(int,curve["invariants"]["group_invariant_factors"])
        if d1*d2 != order or d2%d1:
            raise ValueError("invalid point-group invariant factors")
        check(f"ellcard(E)=={order}",cid+" cardinality","cardinality_checks")
        gp_group = f"[{d2}]" if d1==1 else f"[{d2},{d1}]"
        check(f"ellgroup(E)=={gp_group}",cid+" point group","group_structure_checks")
        for name,record in curve["named_points"].items():
            point = _point(record["point"],field)
            point_order = int(record["order"])
            check(f"ellorder(E,{point},{order})=={point_order}",
                  cid+" point order "+name,"named_point_order_checks")
        subgroup = curve["subgroup"]
        r = int(subgroup["order"])
        generator = _point(subgroup["generator"],field)
        check(f"ellorder(E,{generator},{order})=={r}",
              cid+" subgroup generator order","subgroup_generator_order_checks")
        for vector in vector_by_curve[cid]:
            operation = vector["operation"]
            inputs = vector["inputs"]
            P = _point(inputs["P"],field)
            if operation=="add":
                Q = _point(inputs["Q"],field)
                expression = f"elladd(E,{P},{Q})"
            elif operation=="double":
                expression = f"elladd(E,{P},{P})"
            elif operation=="neg":
                expression = f"ellneg(E,{P})"
            else:
                scalar = int(inputs["scalar"])
                expression = f"ellmul(E,{P},{scalar})"
            expected = _point(vector["expected"]["point"],field)
            check(f"{expression}=={expected}",str(vector["id"]),"curve_arithmetic_checks")
    lines.append("print([checks,version()]);")
    try:
        result = subprocess.run([str(executable),"-f","-q","-s","256M"],
                                input="\n".join(lines)+"\n",text=True,
                                capture_output=True,check=False)
    except FileNotFoundError as exc:
        raise RuntimeError("PARI/GP executable was not found; set CURVE_REPORT_GP or pass gp_path") from exc
    output = result.stdout.strip()
    if result.returncode or result.stderr.strip():
        details = output or result.stderr.strip()
        raise RuntimeError("PARI verification failed: "+details[:2000])
    try:
        actual_checks, version = json.loads(output)
    except (ValueError,TypeError) as exc:
        raise RuntimeError("PARI oracle returned an unexpected summary: "+output[:500]) from exc
    if actual_checks != sum(counts.values()):
        raise RuntimeError("PARI oracle executed an unexpected number of checks")
    return {"status":"passed","oracle":"PARI/GP",
            "pari_version":".".join(map(str,version)),
            "fields":len(fields),"curves":len(curves),
            **counts,"total_checks":actual_checks,
            "vectors_outside_oracle_scope":ignored_vector_count}
