#!/usr/bin/env bash

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

BIN="${TDKC_BIN:-./target/release/tdkc}"
DB="data/influenza_test_db"
READS="data/test_reads/mock.fastq"
EXPECTED_DIR="tests/expected"
ACTUAL_DIR="tests/_actual"

if [ ! -x "$BIN" ]; then
  echo "ERROR: tdkc binary not found at $BIN" >&2
  echo "Build it first with: cargo build --release" >&2
  exit 2
fi
if [ ! -d "$DB" ]; then
  echo "ERROR: test DB not found at $DB" >&2; exit 2
fi
if [ ! -f "$READS" ]; then
  echo "ERROR: test reads not found at $READS" >&2; exit 2
fi
if [ ! -d "$EXPECTED_DIR" ]; then
  echo "ERROR: expected-output dir not found at $EXPECTED_DIR" >&2; exit 2
fi

rm -rf "$ACTUAL_DIR"
mkdir -p "$ACTUAL_DIR"

PASS=0
FAIL=0
FAILED_TESTS=()

normalize() {
  sed -E 's/[[:space:]]+$//' "$1" | sed -e :a -e '/^$/{$d;N;ba' -e '}'
}

diff_normalized() {
  local label="$1" actual="$2" expected="$3"
  if [ ! -f "$expected" ]; then
    echo "  [SKIP] No fixture for $label (expected at $expected) — actual was:"
    sed 's/^/        /' "$actual" || true
    return 0
  fi
  if diff -u <(normalize "$expected") <(normalize "$actual") > "$ACTUAL_DIR/${label}.diff" 2>&1; then
    echo "  [OK]   $label matches $expected"
    rm -f "$ACTUAL_DIR/${label}.diff"
    return 0
  else
    echo "  [FAIL] $label differs from $expected:"
    sed 's/^/        /' "$ACTUAL_DIR/${label}.diff"
    return 1
  fi
}

run_case() {
  local name="$1"; shift
  echo
  echo "Test: $name"
  echo "    cmd: $BIN $*"
  if ! "$BIN" "$@" > "$ACTUAL_DIR/${name}.stdout" 2> "$ACTUAL_DIR/${name}.stderr"; then
    echo "  [FAIL] tdkc exited non-zero. stderr:"
    sed 's/^/        /' "$ACTUAL_DIR/${name}.stderr"
    FAIL=$((FAIL+1)); FAILED_TESTS+=("$name (non-zero exit)")
    return 1
  fi
  return 0
}

check_case() {
  local name="$1"; shift
  local ok=1
  for f in "$@"; do
    local base; base="$(basename "$f")"
    if ! diff_normalized "${name}__${base}" "$f" "$EXPECTED_DIR/${name}/${base}"; then
      ok=0
    fi
  done
  if [ "$ok" -eq 1 ]; then
    PASS=$((PASS+1))
    echo "  [PASS] $name"
  else
    FAIL=$((FAIL+1)); FAILED_TESTS+=("$name")
    echo "  [FAIL] $name"
  fi
}


# Test 1: default query (single-end)

T1_OUT="$ACTUAL_DIR/default"
mkdir -p "$(dirname "$T1_OUT")"
if run_case default \
    query -d "$DB" -1 "$READS" -j 1 -o "$T1_OUT"; then
  check_case default "$T1_OUT.report" "$T1_OUT.output"
fi

# Test 2: query with -a (accession tracking in hit pattern)

T2_OUT="$ACTUAL_DIR/accession"
if run_case accession \
    query -d "$DB" -1 "$READS" -j 1 -a -o "$T2_OUT"; then
  check_case accession "$T2_OUT.report" "$T2_OUT.output"
fi


# Test 3: query with -a -c (accession counts .strain.txt)
# 
T3_OUT="$ACTUAL_DIR/strain"
if run_case strain \
    query -d "$DB" -1 "$READS" -j 1 -a -c -o "$T3_OUT"; then
  check_case strain "$T3_OUT.report" "$T3_OUT.output" "$T3_OUT.strain.txt"
fi

# ----------------------------------------------------------------------
echo
echo "Query integration test summary: $PASS passed, $FAIL failed"
if [ "$FAIL" -gt 0 ]; then
  echo "Failed tests:"
  for t in "${FAILED_TESTS[@]}"; do
    echo "  - $t"
  done
  exit 1
fi
echo "All query tests passed."