#!/usr/bin/env python3
"""Validate the published golden corpus and measured results against their schemas."""
import json
from pathlib import Path

import jsonschema

SUITE = Path(__file__).resolve().parents[1] / 'research' / 'golden-smallfields'
CHECKS = {
    'schemas/dataset.schema.json': ['golden_ecbench.json'],
    'schemas/algebra_cases.schema.json': [
        'build/algebra_cases.json', 'build/algebra_coverage_cases.json'],
    'schemas/algebra_results.schema.json': [
        'build/algebra_buchberger.json', 'build/algebra_f5b.json',
        'build/algebra_coverage_buchberger.json', 'build/algebra_coverage_f5b.json'],
}


def main():
    checked = 0
    for schema_name, instances in CHECKS.items():
        schema = json.loads((SUITE / schema_name).read_text())
        jsonschema.Draft202012Validator.check_schema(schema)
        validator = jsonschema.Draft202012Validator(schema, format_checker=jsonschema.FormatChecker())
        for name in instances:
            validator.validate(json.loads((SUITE / name).read_text()))
            checked += 1
            print(f'{name}: pass')
    print(json.dumps({'status': 'pass', 'schemas': len(CHECKS), 'documents_checked': checked}))


if __name__ == '__main__': main()
