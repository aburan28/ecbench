#!/usr/bin/env python3
"""Build and independently verify supplementary GLV2/Frobenius cases."""
import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from ecbench.endomorphisms import build_cases, verify_cases
from ecbench.validate import load_dataset


def build(output=ROOT / 'build' / 'endomorphism_cases.json', data_dir=None):
    data = load_dataset(data_dir)
    provenance = json.loads((ROOT / 'provenance' / 'binary_data.json').read_text())
    artifact = build_cases(data, provenance)
    report = verify_cases(data, artifact)
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(artifact, sort_keys=True, indent=2) + '\n')
    return artifact, report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT / 'build' / 'endomorphism_cases.json')
    parser.add_argument('--data-dir', type=Path)
    args = parser.parse_args()
    artifact, report = build(args.output, args.data_dir)
    print(json.dumps({'counts': artifact['counts'], 'verification': report}, sort_keys=True))
    print(f'Wrote {args.output}')


if __name__ == '__main__':
    main()
