#!/usr/bin/env python3
"""Independent known answers for rational maps on the prime-field examples."""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess

from prime_vectors import PrimeCurve, sqrt_mod


def evaluate(operation, point, parameter, p):
    if point is None:
        return None
    x,y = point
    if operation=="degree2_isogeny":
        if x==0:
            assert y==0
            return None
        inverse = pow(x,-1,p)
        return (x-inverse)%p,y*(1+inverse*inverse)%p
    if operation=="gaussian_automorphism":
        return (-x)%p,parameter*y%p
    if operation=="eisenstein_automorphism":
        return parameter*x%p,y
    raise ValueError(operation)


def build(source, raw):
    records = {record["id"]:record for record in source["curves"]}
    result = {"schema":"ecbench-prime-maps-raw-v1","vectors":[],
              "scope":"Explicit auxiliary degree-two targets and rational CM unit automorphisms; no added primary curves"}
    identities = 0
    for curve in raw["curves"]:
        cid = curve["id"]
        record = records[cid]
        p = curve["p"]
        E = PrimeCurve(p,curve["A"],curve["B"])
        specifications = []
        if record["family"]=="E4" and record["index"]==0:
            assert curve["A"]==p-1 and curve["B"]==0
            specifications.append(("degree2_isogeny",None,4,0))
        if record["family"]=="E4" and record["ordinary"]:
            i = sqrt_mod(-1,p)
            assert i is not None and i*i%p==p-1
            specifications.append(("gaussian_automorphism",i,curve["A"],curve["B"]))
        if record["family"]=="E3" and record["ordinary"]:
            zeta = pow(record["g"],(p-1)//3,p)
            assert zeta!=1 and pow(zeta,3,p)==1
            specifications.append(("eisenstein_automorphism",zeta,curve["A"],curve["B"]))
        for operation,parameter,A,B in specifications:
            target = PrimeCurve(p,A,B)
            for name in ["O","G","R","P","T2"]:
                if name not in curve["points"]:
                    continue
                point = curve["points"][name]
                point = None if point is None else tuple(point)
                image = evaluate(operation,point,parameter,p)
                assert target.on_curve(image)
                if operation=="gaussian_automorphism":
                    image2=evaluate(operation,image,parameter,p)
                    assert image2==E.neg(point)
                    assert evaluate(operation,evaluate(operation,image2,parameter,p),parameter,p)==point
                    identities+=2
                elif operation=="eisenstein_automorphism":
                    assert evaluate(operation,evaluate(operation,image,parameter,p),parameter,p)==point
                    identities+=1
                result["vectors"].append({"curve_id":cid,"operation":operation,
                    "point_name":name,"point":point,"parameter":parameter,
                    "target_coefficients":{"a1":0,"a2":0,"a3":0,"a4":A,"a6":B},
                    "result":image})
            R=tuple(curve["points"]["R"]);G=tuple(curve["points"]["G"])
            assert evaluate(operation,E.add(R,G),parameter,p)==target.add(
                evaluate(operation,R,parameter,p),evaluate(operation,G,parameter,p))
            identities+=1
    result["verification"]={"target_on_curve_checks":len(result["vectors"]),
                            "unit_order_and_homomorphism_checks":identities}
    return result


def pari_crosscheck(data, raw):
    executable=(os.environ.get("CURVE_REPORT_GP") or shutil.which("gp") or
                "/tmp/pari-local/usr/bin/gp")
    curves={curve["id"]:curve for curve in raw["curves"]}
    lines=[]
    expected=[]
    for vector in data["vectors"]:
        curve=curves[vector["curve_id"]];p=curve["p"]
        operation=vector["operation"]
        point=vector["point"]
        image=vector["result"]
        lines.append(f'E=ellinit([Mod({curve["A"]},{p}),Mod({curve["B"]},{p})]);')
        gp_point="[0]" if point is None else "["+",".join(map(str,point))+"]"
        if operation=="degree2_isogeny":
            lines.append(f'iso=ellisogeny(E,[Mod(0,{p}),Mod(0,{p})]);')
            lines.append('print(lift(iso[1]));')
            expected.append([0,0,0,4,0])
            expression=f"ellisogenyapply(iso[2],{gp_point})"
        elif point is None:
            expression="[0]"
        elif operation=="gaussian_automorphism":
            expression=f'[Mod({-point[0]},{p}),Mod({vector["parameter"]*point[1]},{p})]'
        else:
            expression=f'[Mod({vector["parameter"]*point[0]},{p}),Mod({point[1]},{p})]'
        lines.append(f"print(lift({expression}));")
        expected.append([0] if image is None else list(image))
    result=subprocess.run([executable,"-f","-q","-s","256M"],input="\n".join(lines)+"\n",
                          text=True,capture_output=True,check=True)
    assert not result.stderr,result.stderr
    rows=[json.loads(line) for line in result.stdout.splitlines()]
    assert rows==expected,(rows,expected)
    data["verification"]["pari_version"]="2.17.2"
    data["verification"]["pari_map_answer_checks"]=len(data["vectors"])
    data["verification"]["pari_degree2_target_checks"]=sum(
        vector["operation"]=="degree2_isogeny" for vector in data["vectors"])


def main():
    root=Path(__file__).resolve().parents[1]
    parser=argparse.ArgumentParser()
    parser.add_argument("--source",type=Path,default=root/"provenance"/"normalized_data.json")
    parser.add_argument("--raw-vectors",type=Path,default=root/"build"/"prime_vectors.json")
    parser.add_argument("--output",type=Path,default=root/"build"/"prime_maps.json")
    parser.add_argument("--pari-crosscheck",action="store_true")
    args=parser.parse_args()
    source=json.loads(args.source.read_text());raw=json.loads(args.raw_vectors.read_text())
    answer=build(source,raw)
    if args.pari_crosscheck:
        pari_crosscheck(answer,raw)
    args.output.write_text(json.dumps(answer,indent=2)+"\n")
    print(f"Wrote {len(answer['vectors'])} explicit prime-field map vectors")


if __name__=="__main__":
    main()
