#!/usr/bin/env python3
"""Package source, golden data, actual measurements, report, and wheel."""
import hashlib
from pathlib import Path
import zipfile

ROOT=Path(__file__).resolve().parents[1]
ARCHIVE=ROOT.parent/'ecbench-golden-v1.0.0.zip'
BUILD_FILES={
    'prime_vectors.json','binary_vectors.json','prime_maps.json','binary_maps.json',
    'prime_pari_crosscheck.json','prime_algebra_cases.json','verification.json','optimized_verification.json',
    'reference_benchmark.json','endomorphism_cases.json','endomorphism_verification.json','endomorphism_benchmark.json',
    'glv4.json','glv4_verification.json','glv4_benchmark.json','algebra_cases.json','algebra_coverage_cases.json',
    'algebra_buchberger.json','algebra_f5b.json','algebra_coverage_buchberger.json','algebra_coverage_f5b.json',
    'algebra_summary.csv','algebra_generation_benchmark.json','algebra_fixture_verification.json',
    'battery_verification.json','reproducibility.json','schema_verification.json','algebra_schema_verification.json','unit_tests.log',
    'installed_wheel_verification.json',
}


def main():
    paths=[]
    for directory in ('ecbench','tools','schemas','tests','docs','provenance'):
        paths.extend(path for path in (ROOT/directory).rglob('*') if path.is_file()
            and '__pycache__' not in path.parts and path.suffix!='.pyc')
    for name in ('README.md','pyproject.toml','golden_ecbench.json','golden_ecbench.jsonl','ecbench_golden-1.0.0-py3-none-any.whl'):
        paths.append(ROOT/name)
    paths.extend(path for path in (ROOT/'build').iterdir() if path.is_file() and path.name in BUILD_FILES)
    external=ROOT/'build'/'external'/'example'
    if external.exists():paths.extend(path for path in external.iterdir() if path.is_file())
    paths=sorted(set(paths))
    checksum=ROOT/'SHA256SUMS'
    checksum.write_text(''.join(f'{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.relative_to(ROOT)}\n' for path in paths))
    paths.append(checksum)
    with zipfile.ZipFile(ARCHIVE,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=9) as archive:
        for path in paths:
            info=zipfile.ZipInfo('ecbench-golden-v1.0.0/'+str(path.relative_to(ROOT)),date_time=(2026,10,12,0,0,0))
            info.compress_type=zipfile.ZIP_DEFLATED
            info.external_attr=(0o100644<<16)
            archive.writestr(info,path.read_bytes())
    with zipfile.ZipFile(ARCHIVE) as archive:
        if archive.testzip() is not None:raise RuntimeError('Release ZIP integrity failure')
    digest=hashlib.sha256(ARCHIVE.read_bytes()).hexdigest()
    ARCHIVE.with_suffix('.zip.sha256').write_text(f'{digest}  {ARCHIVE.name}\n')
    print(f'{ARCHIVE} ({ARCHIVE.stat().st_size} bytes; {len(paths)} files)')
    print(f'SHA256 {digest}')


if __name__=='__main__':main()
