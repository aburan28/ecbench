#!/usr/bin/env python3
"""Package the research source snapshot without overwriting the frozen release."""
import hashlib
from pathlib import Path
import zipfile

from build_release import BUILD_FILES

ROOT = Path(__file__).resolve().parents[1]
DESTINATION = ROOT.parent / 'ecbench-research-followup-2026-10-15.zip'
PREFIX = 'ecbench-research-followup-2026-10-15/'


def main():
    files = []
    for directory in ('ecbench', 'tools', 'schemas', 'tests', 'docs', 'provenance', 'followups'):
        files.extend(path for path in (ROOT / directory).rglob('*') if path.is_file()
            and '__pycache__' not in path.parts and path.suffix != '.pyc')
    for name in ('README.md', 'FOLLOWUP.md', 'pyproject.toml', 'golden_ecbench.json', 'golden_ecbench.jsonl',
                 'ecbench_golden-1.0.0-py3-none-any.whl', 'SHA256SUMS'):
        files.append(ROOT / name)
    files.extend(path for path in (ROOT / 'build').iterdir() if path.is_file() and path.name in BUILD_FILES)
    example = ROOT / 'build' / 'external' / 'example'
    files.extend(path for path in example.iterdir() if path.is_file())
    files = sorted(set(files))
    sums = ROOT / 'FOLLOWUP_SHA256SUMS'
    sums.write_text(''.join(f'{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.relative_to(ROOT)}\n' for path in files))
    files.append(sums)
    with zipfile.ZipFile(DESTINATION, 'w', compression=zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
        for path in files:
            item = zipfile.ZipInfo(PREFIX + str(path.relative_to(ROOT)), date_time=(2026,10,15,0,0,0))
            item.compress_type = zipfile.ZIP_DEFLATED
            item.external_attr = 0o100644 << 16
            archive.writestr(item, path.read_bytes())
    with zipfile.ZipFile(DESTINATION) as archive:
        if archive.testzip() is not None: raise RuntimeError('ZIP integrity failure')
        for line in archive.read(PREFIX + 'FOLLOWUP_SHA256SUMS').decode().splitlines():
            expected, name = line.split('  ', 1)
            if hashlib.sha256(archive.read(PREFIX + name)).hexdigest() != expected:
                raise RuntimeError('Archive checksum mismatch: ' + name)
    actual = hashlib.sha256(DESTINATION.read_bytes()).hexdigest()
    DESTINATION.with_suffix('.zip.sha256').write_text(f'{actual}  {DESTINATION.name}\n')
    print(f'{DESTINATION} ({DESTINATION.stat().st_size} bytes; {len(files)} files)')
    print('SHA256 ' + actual)


if __name__ == '__main__': main()
