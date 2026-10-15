#!/usr/bin/env python3
"""Independent prime-field arithmetic and deterministic ecbench answers.

This implementation uses only Python's standard library. It deliberately uses
affine coordinates, integer modular inverses, and Tonelli-Shanks, independent
of the PARI point-count computations and CM calculation implementation.

Raw output contains integer coordinates and null for infinity. SEC1 vectors
use fixed-width, big-endian coordinates and enforce canonical field values.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import shutil
import subprocess
from pathlib import Path
from typing import Optional

Point = Optional[tuple[int, int]]


def sqrt_mod(a: int, p: int) -> Optional[int]:
    """Return the smaller square root, or None if a is a nonresidue."""
    a %= p
    if a == 0:
        return 0
    if pow(a, (p - 1) // 2, p) != 1:
        return None
    if p % 4 == 3:
        x = pow(a, (p + 1) // 4, p)
        return min(x, p - x)
    q = p - 1
    s = 0
    while q % 2 == 0:
        q //= 2
        s += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1:
        z += 1
    c = pow(z, q, p)
    x = pow(a, (q + 1) // 2, p)
    t = pow(a, q, p)
    m = s
    while t != 1:
        i = 1
        tt = t * t % p
        while tt != 1:
            tt = tt * tt % p
            i += 1
            if i >= m:
                raise ArithmeticError("Tonelli-Shanks invariant failed")
        b = pow(c, 1 << (m - i - 1), p)
        x = x * b % p
        t = t * b * b % p
        c = b * b % p
        m = i
    assert x * x % p == a
    return min(x, p - x)


class PrimeCurve:
    def __init__(self, p: int, A: int, B: int):
        self.p, self.A, self.B = p, A, B
        self.width = (p.bit_length() + 7) // 8
        assert p > 3 and 0 <= A < p and 0 <= B < p
        assert (4 * A**3 + 27 * B**2) % p != 0

    def on_curve(self, P: Point) -> bool:
        if P is None:
            return True
        x, y = P
        return (0 <= x < self.p and 0 <= y < self.p and
                (y * y - x * x * x - self.A * x - self.B) % self.p == 0)

    def neg(self, P: Point) -> Point:
        return None if P is None else (P[0], (-P[1]) % self.p)

    def add(self, P: Point, Q: Point) -> Point:
        if P is None:
            return Q
        if Q is None:
            return P
        x1, y1 = P
        x2, y2 = Q
        p = self.p
        if x1 == x2:
            if (y1 + y2) % p == 0:
                return None
            slope = (3 * x1 * x1 + self.A) * pow(2 * y1, -1, p) % p
        else:
            slope = (y2 - y1) * pow(x2 - x1, -1, p) % p
        x3 = (slope * slope - x1 - x2) % p
        result = x3, (slope * (x1 - x3) - y1) % p
        assert self.on_curve(result)
        return result

    def mul(self, n: int, P: Point) -> Point:
        if n < 0:
            return self.mul(-n, self.neg(P))
        result = None
        while n:
            if n & 1:
                result = self.add(result, P)
            P = self.add(P, P)
            n >>= 1
        return result

    def scan(self, start: int = 0):
        for x in range(start, min(self.p, 100000)):
            y = sqrt_mod((x**3 + self.A * x + self.B) % self.p, self.p)
            if y is not None:
                yield (x, y)

    def encode(self, P: Point, compressed: bool = True) -> bytes:
        if P is None:
            return b"\x00"
        if not self.on_curve(P):
            raise ValueError("point_not_on_curve_or_noncanonical")
        x, y = P
        xb = x.to_bytes(self.width, "big")
        if compressed:
            return bytes([2 + (y & 1)]) + xb
        return b"\x04" + xb + y.to_bytes(self.width, "big")

    def decode(self, data: bytes) -> Point:
        if data == b"\x00":
            return None
        if not data:
            raise ValueError("invalid_encoding")
        tag = data[0]
        if tag in (2, 3):
            if len(data) != self.width + 1:
                raise ValueError("invalid_encoding")
            x = int.from_bytes(data[1:], "big")
            if x >= self.p:
                raise ValueError("noncanonical_coordinate")
            y = sqrt_mod((x**3 + self.A * x + self.B) % self.p, self.p)
            if y is None:
                raise ValueError("point_not_on_curve")
            if (y & 1) != (tag & 1):
                y = (-y) % self.p
            if (y & 1) != (tag & 1):
                raise ValueError("invalid_y_selector")
            return x, y
        if tag == 4:
            if len(data) != 2 * self.width + 1:
                raise ValueError("invalid_encoding")
            x = int.from_bytes(data[1:self.width+1], "big")
            y = int.from_bytes(data[self.width+1:], "big")
            if x >= self.p or y >= self.p:
                raise ValueError("noncanonical_coordinate")
            if not self.on_curve((x, y)):
                raise ValueError("point_not_on_curve")
            return x, y
        raise ValueError("invalid_encoding")


def field_vectors(p: int) -> list[dict]:
    v = []
    def put(op, args, result):
        v.append({"op": op, "args": args, "result": result})
    for a, b in [(0,0),(1,p-1),(p-1,p-1),(p//2,p//3)]:
        put("add", [a,b], (a+b)%p)
    for a, b in [(0,1),(1,1),(p-1,p//2)]:
        put("sub", [a,b], (a-b)%p)
    for a, b in [(0,p-1),(1,p-1),(p-1,p-1),(p//2,p//3)]:
        put("mul", [a,b], a*b%p)
    for a in [0,1,p-1,p//2]:
        put("square", [a], a*a%p)
    for a in [1,2,p-1,p//3]:
        put("inverse", [a], pow(a,-1,p))
    v.append({"op":"inverse","args":[0],"error":"inverse_zero"})
    for a in [0,p-1,p,p+1,2*p-1,(1<<(2*p.bit_length()))+12345,-1]:
        put("reduce", [a], a%p)
    return v


def cyclic_contains(E: PrimeCurve, G: Point, r: int, H: Point) -> bool:
    """Exact small-order BSGS membership for the rank-two torsion traps."""
    assert r <= 100000
    m = math.isqrt(r) + 1
    babies = {}
    point = None
    for j in range(m):
        babies[point] = j
        point = E.add(point,G)
    step = E.neg(E.mul(m,G))
    gamma = H
    for i in range(m):
        j = babies.get(gamma)
        if j is not None and i*m+j < r:
            return True
        gamma = E.add(gamma,step)
    return False


def exact_point_order(E: PrimeCurve, P: Point, N: int, factors: dict) -> int:
    order = N
    assert E.mul(order,P) is None
    for prime in map(int,factors):
        while order%prime==0 and E.mul(order//prime,P) is None:
            order//=prime
    assert E.mul(order,P) is None
    for prime in map(int,factors):
        if order%prime==0:
            assert E.mul(order//prime,P) is not None
    return order


def curve_vectors(record: dict) -> dict:
    p, A, B = record["q"], record["A"], record["B"]
    N, r = record["N"], record["r"]
    d1, d2 = record["group"]
    E = PrimeCurve(p,A,B)
    assert d1*d2 == N and d2%d1 == 0 and d2%r == 0
    iterator = E.scan()
    P = next(iterator)
    Q = next(point for point in iterator if point not in (P,E.neg(P)))
    projector = d2//r
    G = None
    for source in E.scan():
        candidate = E.mul(projector,source)
        if candidate is not None:
            assert E.mul(r,candidate) is None
            G = candidate
            G_source = source
            break
    if G is None:
        raise RuntimeError(f"No subgroup generator found for {record['id']}")
    assert E.mul(N,P) is None and E.mul(N,Q) is None
    assert E.mul(d2,P) is None and E.mul(d2,Q) is None
    # The largest r is already a certified prime in the source bundle. Nonzero
    # G with [r]G=O proves exact order r, including the smooth r=2 examples.
    points = {"O":None,"P":P,"Q":Q,"G":G,"neg_P":E.neg(P),
              "neg_G":E.neg(G),"double_P":E.add(P,P),"double_G":E.add(G,G),
              "P_plus_Q":E.add(P,Q),"mul_17_P":E.mul(17,P),"mul_7_G":E.mul(7,G)}
    # A maximum-order whole-group point is useful for benchmarking even when
    # the first scan points happen to have order2. This is one point of order
    # d2; it is not a claim of a basis for a noncyclic group.
    exponent_primes = [int(prime) for prime in record["factors"] if d2%int(prime)==0]
    R = next(point for point in E.scan()
             if E.mul(d2,point) is None and
                all(E.mul(d2//prime,point) is not None for prime in exponent_primes))
    points["R"] = R
    points["neg_R"] = E.neg(R)
    points["double_R"] = E.add(R,R)
    points["R_plus_Q"] = E.add(R,Q)
    rank = 2 if d1%r==0 else 1
    H = None
    H_source = None
    H_derivation = None
    if rank==2:
        assert r<=100000
        for source in E.scan():
            # A single d2/r map can kill the smaller r-primary factor when
            # its exponent is below that of d2 (M61 E4 C2 x C2^60). Also
            # inspect directly encountered r-torsion, e.g. its order2 roots.
            candidates = []
            if E.mul(r,source) is None:
                candidates.append((source,"direct_r_torsion"))
            candidates.append((E.mul(projector,source),"projector_d2_over_r"))
            for candidate,derivation in candidates:
                if candidate is not None and not cyclic_contains(E,G,r,candidate):
                    H,H_source = candidate,source
                    H_derivation = derivation
                    assert E.mul(r,H) is None
                    points["H"] = H
                    break
            if H is not None:
                break
        assert H is not None
    T2 = None
    if d2%2==0:
        for source in E.scan():
            candidate = E.mul(d2//2,source)
            if candidate is not None:
                assert candidate[1]==0 and E.mul(2,candidate) is None
                T2 = candidate
                points["T2"] = candidate
                break
        assert T2 is not None
    cv = []
    def put(op,args,result):
        cv.append({"op":op,"args":args,"result":result})
    for args in [(P,Q),(P,None),(None,Q),(None,None),(P,E.neg(P)),(G,G)]:
        put("add",list(args),E.add(*args))
    for point in [P,G,None]:
        put("double",[point],E.add(point,point))
        put("neg",[point],E.neg(point))
    wide_scalar = (1 << (p.bit_length()-1)) + 0x12345
    for n,point in [(0,P),(1,P),(2,P),(3,P),(17,P),(-1,P),(r,P),(N,P),
                    (0,G),(1,G),(2,G),(3,G),(7,G),(-1,G),
                    (r-1,G),(r,G),(r+1,G),(wide_scalar,G),((1<<128)+1,G),
                    (wide_scalar,P),((1<<128)+1,P)]:
        put("mul",[n,point],E.mul(n,point))
    for n in [0,1,2,r,N,wide_scalar]:
        put("mul",[n,R],E.mul(n,R))
    put("mul",[-1,R],E.neg(R))
    put("double",[R],E.add(R,R))
    put("neg",[R],E.neg(R))
    put("add",[R,Q],E.add(R,Q))
    if T2 is not None:
        put("double",[T2],None)
        put("add",[T2,T2],None)
        put("neg",[T2],T2)
        put("mul",[2,T2],None)
        put("mul",[3,T2],T2)
    if H is not None:
        put("mul",[r,H],None)
    enc = []
    for name,point in points.items():
        compressed = E.encode(point,True)
        uncompressed = E.encode(point,False)
        assert E.decode(compressed)==point and E.decode(uncompressed)==point
        enc.append({"name":name,"point":point,
                    "sec1_compressed_hex":compressed.hex(),
                    "sec1_uncompressed_hex":uncompressed.hex()})
    invalid_point = next((0,y) for y in range(3) if not E.on_curve((0,y)))
    invalid = [
        {"op":"validate_point","point":invalid_point,"result":False,"reason":"off_curve"},
        {"op":"validate_point","point":(p,0),"result":False,"reason":"noncanonical_x"},
        {"op":"validate_point","point":(0,p),"result":False,"reason":"noncanonical_y"},
        {"op":"validate_point","point":(-1,0),"result":False,"reason":"negative_coordinate"},
    ]
    malformed = [b"",b"\x00\x00",b"\x02",b"\x05"+bytes(E.width),
                 b"\x02"+p.to_bytes(E.width,"big"),
                 b"\x02"+b"\xff"*E.width,
                 b"\x04"+invalid_point[0].to_bytes(E.width,"big")+
                              invalid_point[1].to_bytes(E.width,"big")]
    # An x whose RHS is a nonresidue is a distinct compressed-decoding failure.
    nonresidue_x = next(x for x in range(100000)
                       if sqrt_mod((x**3+A*x+B)%p,p) is None)
    malformed.append(b"\x02"+nonresidue_x.to_bytes(E.width,"big"))
    # Odd selector on a zero-y point is impossible and must be rejected.
    zero_y = T2
    if zero_y is not None:
        malformed.append(b"\x03"+zero_y[0].to_bytes(E.width,"big"))
    for data in malformed:
        try:
            E.decode(data)
        except ValueError as exc:
            invalid.append({"op":"decode_point","encoding_hex":data.hex(),"error":str(exc)})
        else:
            raise AssertionError("Malformed point accepted")
    subgroup = [
        {"op":"public_key_valid","point":G,"allow_infinity":False,"result":True},
        {"op":"public_key_valid","point":None,"allow_infinity":False,"result":False},
        {"op":"subgroup_contains","point":None,"allow_infinity":True,"result":True},
        {"op":"subgroup_contains","point":G,"allow_infinity":True,"result":True},
        {"op":"public_key_valid","point":invalid_point,"allow_infinity":False,"result":False},
    ]
    if H is not None:
        subgroup.extend([
            {"op":"is_r_torsion","point":H,"allow_infinity":True,"result":True},
            {"op":"subgroup_contains","point":H,"allow_infinity":True,"result":False},
            {"op":"public_key_valid","point":H,"allow_infinity":False,"result":False},
        ])
    if N != r:
        outside = next((point for point in E.scan() if E.mul(r,point) is not None),None)
        if outside is not None:
            subgroup.append({"op":"public_key_valid","point":outside,"allow_infinity":False,"result":False})
    if d2>r:
        subgroup.append({"op":"public_key_valid","point":R,"allow_infinity":False,"result":False})
    result = {"id":record["id"],"field_id":record["field"],"p":p,"A":A,"B":B,
              "order":N,"group_invariant_factors":[d1,d2],"r":r,
              "coordinate_bytes":E.width,"subgroup_generator":G,
              "subgroup_r_torsion_rank":rank,"subgroup_nonunique":rank==2,
              "alternate_r_torsion_point":H,
              "deterministic_scan":{"rule":"x ascending from zero; minimum integer square root",
                                    "P_scan_x":P[0],"Q_scan_x":Q[0],
                                    "R_scan_x":R[0],"R_rule":"first scanned point of exact order d2",
                                    "G_scan_x":G_source[0],"G_source":G_source,
                                    "H_scan_x":H_source[0] if H_source else None,
                                    "H_source":H_source,
                                    "projector":projector,"projector_formula":"d2/r"},
              "points":points,"field_vectors":field_vectors(p),"curve_vectors":cv,
              "encoding_vectors":enc,"invalid_vectors":invalid,
              "subgroup_validation_vectors":subgroup}
    result["point_orders"] = {name:exact_point_order(E,point,N,record["factors"])
                              for name,point in points.items()}
    result["point_order_factors"] = {}
    for name,order in result["point_orders"].items():
        fs = {}
        rem = order
        for prime in map(int,record["factors"]):
            e = 0
            while rem%prime==0:
                rem//=prime
                e+=1
            if e:
                fs[str(prime)] = e
        assert rem==1
        result["point_order_factors"][name] = fs
    assert result["point_orders"]["R"]==d2
    result["deterministic_scan"]["H_derivation"] = H_derivation
    return result


def crosscheck_pari(data: dict, output: Path) -> None:
    """Independent PARI known-answer and point-order crosscheck, when available."""
    gp_executable = (os.environ.get("CURVE_REPORT_GP") or shutil.which("gp") or
                     "/tmp/pari-local/usr/bin/gp")
    lines = []
    expected = []
    def point_text(P):
        return "[0]" if P is None else "["+",".join(map(str,P))+"]"
    arithmetic_count = 0
    order_count = 0
    for curve in data["curves"]:
        p = curve["p"]
        lines.append(f'E=ellinit([Mod({curve["A"]},{p}),Mod({curve["B"]},{p})]);')
        lines.append('print([ellcard(E),ellgroup(E)]);')
        expected.append(("group",curve["id"],[curve["order"],curve["group_invariant_factors"]]))
        for vector in curve["curve_vectors"]:
            op,args = vector["op"],vector["args"]
            if op=="add":
                expression=f"elladd(E,{point_text(args[0])},{point_text(args[1])})"
            elif op=="double":
                expression=f"elladd(E,{point_text(args[0])},{point_text(args[0])})"
            elif op=="neg":
                expression=f"ellneg(E,{point_text(args[0])})"
            elif op=="mul":
                expression=f"ellmul(E,{point_text(args[1])},{args[0]})"
            else:
                raise ValueError(op)
            lines.append(f"print(lift({expression}));")
            expected.append(("arithmetic",curve["id"],vector["result"]))
            arithmetic_count+=1
        for name,point in curve["points"].items():
            lines.append(f'print(ellorder(E,{point_text(point)},{curve["order"]}));')
            expected.append(("order",curve["id"],curve["point_orders"][name]))
            order_count+=1
    result=subprocess.run([gp_executable,"-f","-q","-s","256M"],
                          input="\n".join(lines)+"\n",text=True,
                          capture_output=True,check=True)
    assert not result.stderr,result.stderr
    rows=[json.loads(line) for line in result.stdout.splitlines()]
    assert len(rows)==len(expected)
    for got,(kind,cid,want) in zip(rows,expected):
        if kind=="arithmetic":
            got=None if got==[0] else got
            want=None if want is None else list(want)
        elif kind=="group":
            N,group=got
            got=[N,[group[1] if len(group)>1 else 1,group[0]]]
        assert got==want,(kind,cid,got,want)
    summary={"status":"passed","pari_version":"2.17.2","curves":len(data["curves"]),
             "arithmetic_vectors":arithmetic_count,"point_order_vectors":order_count,
             "cardinality_and_group_crosschecks":len(data["curves"]),
             "operations":["add","double","neg","mul","ellorder","ellcard","ellgroup"]}
    output.write_text(json.dumps(summary,indent=2)+"\n")
    print(f"PARI independently agrees with {arithmetic_count} arithmetic vectors and {order_count} point orders")


def main():
    project_root = Path(__file__).resolve().parents[1]
    default_source = project_root/"provenance"/"normalized_data.json"
    if not default_source.is_file():
        workspace_source = Path("/workspace/curve_report/normalized_data.json")
        if workspace_source.is_file():
            default_source = workspace_source
    parser = argparse.ArgumentParser()
    parser.add_argument("--source","--input",dest="input",type=Path,default=default_source)
    parser.add_argument("--output",type=Path,default=project_root/"build"/"prime_vectors.json")
    parser.add_argument("--pari-crosscheck",action="store_true")
    args = parser.parse_args()
    data = json.loads(args.input.read_text())
    curves = [record for record in data["curves"] if record["kind"]=="prime"]
    answer = {"schema":"ecbench-prime-raw-v1","arithmetic":"independent affine Python integers and Tonelli-Shanks",
              "infinity":None,"coordinate_encoding":"SEC1 fixed-width big-endian, strict canonical coordinates",
              "curves":[]}
    for record in curves:
        vector = curve_vectors(record)
        answer["curves"].append(vector)
        print(record["id"],"Gscan",vector["deterministic_scan"]["G_scan_x"],flush=True)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(answer,indent=2)+"\n")
    print(f"Wrote {len(curves)} prime curves to {args.output}")
    if args.pari_crosscheck:
        crosscheck_pari(answer,args.output.parent/"prime_pari_crosscheck.json")


if __name__=="__main__":
    main()
