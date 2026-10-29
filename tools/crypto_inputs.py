#!/usr/bin/env python3
"""Print the sparse-checkout patterns the pinned crypto crate needs to build.

Usage: crypto_inputs.py DEPS_CRYPTO_DIR

The crate is a path dependency, so Cargo needs every path crate it names
(followed through each crate's own manifest) and every file a library or
path-crate source embeds with include_str!/include_bytes!.  Both are read
from the checked-out source, so the list follows the pin instead of going
stale.  Paths come from the source tree of `HEAD`; nothing is executed.
"""

import re
import subprocess
import sys
from pathlib import PurePosixPath

PATH_DEP = re.compile(r'=\s*\{[^}\n]*\bpath\s*=\s*"([^"\n]+)"')
EMBED = re.compile(r'include_(?:str|bytes)!\(\s*"([^"\n]+)"')
EMBED_ROOT = re.compile(
    r'include_(?:str|bytes)!\(\s*concat!\(\s*env!\("CARGO_MANIFEST_DIR"\),\s*"([^"\n]+)"'
)


def git(repo, *args):
    return subprocess.run(
        ["git", "-C", repo, *args], check=True, capture_output=True, text=True
    ).stdout


def show(repo, path):
    try:
        return git(repo, "show", f"HEAD:{path}")
    except subprocess.CalledProcessError:
        return None


def norm(path):
    parts = []
    for part in PurePosixPath(path).parts:
        if part == "..":
            if not parts:
                raise SystemExit(f"path escapes the crypto checkout: {path}")
            parts.pop()
        elif part != ".":
            parts.append(part)
    return "/".join(parts)


def main():
    repo = sys.argv[1]
    tracked = set(git(repo, "ls-tree", "-r", "--name-only", "HEAD").split("\n"))

    crates, queue = [], ["."]
    while queue:
        crate = queue.pop(0)
        manifest = show(repo, "Cargo.toml" if crate == "." else f"{crate}/Cargo.toml")
        for rel in PATH_DEP.findall(manifest or ""):
            dep = norm(f"{crate}/{rel}")
            if dep not in crates:
                crates.append(dep)
                queue.append(dep)

    patterns = []
    for crate in crates:
        patterns += [f"/{crate}/Cargo.toml", f"/{crate}/Cargo.lock", f"/{crate}/build.rs", f"/{crate}/src/"]

    roots = ["src/"] + [f"{crate}/src/" for crate in crates]
    inputs = set()
    for file in sorted(tracked):
        if not file.endswith(".rs") or not any(file.startswith(r) for r in roots):
            continue
        text = show(repo, file) or ""
        crate = next((c for c in crates if file.startswith(f"{c}/")), "")
        for rel in EMBED.findall(text):
            inputs.add(norm(f"{PurePosixPath(file).parent}/{rel}"))
        for rel in EMBED_ROOT.findall(text):
            inputs.add(norm(f"{crate}/{rel.lstrip('/')}" if crate else rel.lstrip("/")))
    patterns += [f"/{p}" for p in sorted(inputs) if p in tracked and not p.startswith("src/")]
    print("\n".join(patterns))


if __name__ == "__main__":
    main()
