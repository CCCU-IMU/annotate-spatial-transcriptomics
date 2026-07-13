#!/usr/bin/env bash
set -euo pipefail

REPO="${ANNOTATION_SKILL_REPO:-CCCU-IMU/annotate-spatial-transcriptomics}"
REF="${ANNOTATION_SKILL_REF:-main}"
DEST_ROOT="${CODEX_HOME:-$HOME/.codex}/skills"
SOURCE_DIR=""
SKILL_NAME="annotate-spatial-transcriptomics"

usage() {
  printf '%s\n' \
    "Install the annotate-spatial-transcriptomics Codex Skill." \
    "" \
    "Usage: bash install.sh [--ref TAG_OR_BRANCH] [--dest SKILLS_DIR] [--repo OWNER/REPO] [--source DIR]" \
    "" \
    "Examples:" \
    "  bash install.sh" \
    "  bash install.sh --ref v1.2.0" \
    "  bash install.sh --dest /path/to/.codex/skills"
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --ref) REF="$2"; shift 2 ;;
    --dest) DEST_ROOT="$2"; shift 2 ;;
    --repo) REPO="$2"; shift 2 ;;
    --source) SOURCE_DIR="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) printf 'Unknown argument: %s\n' "$1" >&2; usage >&2; exit 2 ;;
  esac
done

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]:-/dev/null}")" 2>/dev/null && pwd || true)"
TMP_DIR=""

cleanup() {
  if [[ -n "$TMP_DIR" && -d "$TMP_DIR" ]]; then
    rm -rf "$TMP_DIR"
  fi
}
trap cleanup EXIT

if [[ -z "$SOURCE_DIR" && -n "$SCRIPT_DIR" && -f "$SCRIPT_DIR/$SKILL_NAME/SKILL.md" ]]; then
  SOURCE_DIR="$SCRIPT_DIR/$SKILL_NAME"
fi

if [[ -z "$SOURCE_DIR" ]]; then
  for cmd in curl tar; do
    command -v "$cmd" >/dev/null 2>&1 || { printf 'Required command not found: %s\n' "$cmd" >&2; exit 1; }
  done
  TMP_DIR="$(mktemp -d)"
  BRANCH_URL="https://github.com/${REPO}/archive/refs/heads/${REF}.tar.gz"
  if [[ "$REF" == v* ]]; then
    TAG_URL="https://github.com/${REPO}/archive/refs/tags/${REF}.tar.gz"
  fi
  printf 'Downloading %s@%s ...\n' "$REPO" "$REF"
  if [[ -n "${TAG_URL:-}" ]]; then
    curl --fail --location --silent "$TAG_URL" -o "$TMP_DIR/source.tar.gz" \
      || curl --fail --location --silent --show-error "$BRANCH_URL" -o "$TMP_DIR/source.tar.gz"
  else
    curl --fail --location --silent --show-error "$BRANCH_URL" -o "$TMP_DIR/source.tar.gz"
  fi
  TAR_HELP="$(tar --help 2>&1 || true)"
  if [[ "$TAR_HELP" == *"--warning=KEYWORD"* ]]; then
    tar --warning=no-timestamp -xzf "$TMP_DIR/source.tar.gz" -C "$TMP_DIR"
  else
    tar -xzf "$TMP_DIR/source.tar.gz" -C "$TMP_DIR"
  fi
  SOURCE_DIR="$(find "$TMP_DIR" -mindepth 2 -maxdepth 2 -type d -name "$SKILL_NAME" -print -quit)"
fi

if [[ -z "$SOURCE_DIR" || ! -f "$SOURCE_DIR/SKILL.md" ]]; then
  printf 'Could not locate %s/SKILL.md\n' "$SKILL_NAME" >&2
  exit 1
fi

TARGET="$DEST_ROOT/$SKILL_NAME"
mkdir -p "$DEST_ROOT"
STAGE="${TARGET}.installing.$$"
rm -rf "$STAGE"
mkdir -p "$STAGE"
cp -R "$SOURCE_DIR"/. "$STAGE"/

if [[ -e "$TARGET" ]]; then
  BACKUP="${TARGET}.backup.$(date +%Y%m%d%H%M%S)"
  mv "$TARGET" "$BACKUP"
  printf 'Existing installation backed up to %s\n' "$BACKUP"
fi
mv "$STAGE" "$TARGET"

VERIFY_SCRIPT="${SCRIPT_DIR}/scripts/verify_install.py"
if [[ -f "$VERIFY_SCRIPT" ]] && command -v python >/dev/null 2>&1; then
  python "$VERIFY_SCRIPT" "$TARGET"
else
  for required in SKILL.md references/iterative-controller.md references/multi-route-controller.md references/report-contract.md scripts/init_annotation_project.py scripts/check_completion_gate.py scripts/build_report.py; do
    [[ -f "$TARGET/$required" ]] || { printf 'Missing installed file: %s\n' "$required" >&2; exit 1; }
  done
fi

printf '\nInstalled %s to:\n  %s\n' "$SKILL_NAME" "$TARGET"
printf 'Restart Codex or open a new task, then invoke: $%s\n' "$SKILL_NAME"
