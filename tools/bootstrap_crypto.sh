#!/usr/bin/env bash
set -euo pipefail

readonly CRYPTO_SHA=e0378e02f75a035f3104d552693b3a750d3b265d
readonly REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
readonly CRYPTO_DIR="$REPO_ROOT/deps/crypto"

if [[ -e "$CRYPTO_DIR" && ! -d "$CRYPTO_DIR/.git" ]]; then
  echo "Refusing to replace existing $CRYPTO_DIR" >&2
  exit 1
fi

if [[ ! -d "$CRYPTO_DIR/.git" ]]; then
  mkdir -p "$REPO_ROOT/deps"
  git init -q "$CRYPTO_DIR"
  git -C "$CRYPTO_DIR" remote add origin https://github.com/aburan28/crypto.git
  git -C "$CRYPTO_DIR" config core.sparseCheckout true
  git -C "$CRYPTO_DIR" config core.sparseCheckoutCone false
  mkdir -p "$CRYPTO_DIR/.git/info"
fi

if [[ -n "$(git -C "$CRYPTO_DIR" status --porcelain --untracked-files=no)" ]]; then
  echo "Refusing to update dirty $CRYPTO_DIR" >&2
  exit 1
fi

cat > "$CRYPTO_DIR/.git/info/sparse-checkout" <<'PATTERNS'
/Cargo.toml
/Cargo.lock
/src/
/docs/curves/registry.json
/docs/ic/calibration.json
/docs/ic/boundary_targets.json
/docs/ecbench/schema.sql
/docs/ecbench/specs/
/docs/bounds/
PATTERNS

git -C "$CRYPTO_DIR" fetch --no-tags --depth=1 --filter=blob:none origin "$CRYPTO_SHA"
git -C "$CRYPTO_DIR" checkout --detach -q "$CRYPTO_SHA"
# The path crates Cargo resolves and the files the library's sources embed,
# read from the pinned source itself (tools/crypto_inputs.py).
python3 -I "$REPO_ROOT/tools/crypto_inputs.py" "$CRYPTO_DIR" >> "$CRYPTO_DIR/.git/info/sparse-checkout"
git -C "$CRYPTO_DIR" sparse-checkout reapply
actual="$(git -C "$CRYPTO_DIR" rev-parse HEAD)"
if [[ "$actual" != "$CRYPTO_SHA" ]]; then
  echo "Expected crypto $CRYPTO_SHA, got $actual" >&2
  exit 1
fi
echo "crypto source pinned at $actual"
