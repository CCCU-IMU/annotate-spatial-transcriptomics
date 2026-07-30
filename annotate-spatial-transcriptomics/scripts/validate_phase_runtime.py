#!/usr/bin/env python3
"""Validate the exact Python/R runtime before a formal controller phase."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
import subprocess
from pathlib import Path


def run(command: list[str]) -> tuple[bool, str]:
    try:
        result = subprocess.run(
            command, capture_output=True, text=True, check=False, timeout=120
        )
    except subprocess.TimeoutExpired as exc:
        captured = "".join(
            value.decode(errors="replace") if isinstance(value, bytes) else value or ""
            for value in (exc.stdout, exc.stderr)
        ).strip()
        detail = f"; partial output: {captured}" if captured else ""
        return False, f"runtime probe timed out after 120 seconds{detail}"
    message = (result.stdout + result.stderr).strip()
    return result.returncode == 0, message


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validate_semantic_input(specification: str) -> dict[str, object]:
    try:
        name, media_type, raw_path = specification.split(":", 2)
    except ValueError as exc:
        raise ValueError("semantic input must be name:json|tsv|file:path") from exc
    path = Path(raw_path).resolve()
    if not path.is_file():
        raise ValueError(f"semantic input is missing: {name}")
    if media_type == "json":
        json.loads(path.read_text(encoding="utf-8"))
    elif media_type == "tsv":
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.reader(handle, delimiter="\t")
            header = next(reader, [])
        if not header or any(not value for value in header):
            raise ValueError(f"TSV semantic input lacks a valid header: {name}")
    elif media_type != "file":
        raise ValueError(f"unsupported semantic input media type: {media_type}")
    return {
        "name": name, "media_type": media_type,
        "path": str(path), "sha256": sha256(path),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--python", required=True)
    ap.add_argument("--rscript")
    ap.add_argument("--python-script", action="append", type=Path, default=[])
    ap.add_argument("--r-script", action="append", type=Path, default=[])
    ap.add_argument("--python-import", action="append", default=[])
    ap.add_argument("--r-package", action="append", default=[])
    ap.add_argument("--semantic-input", action="append", default=[])
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    errors: list[str] = []
    python_path = shutil.which(args.python) or args.python
    python_ok, python_version = run([python_path, "--version"])
    if not python_ok:
        errors.append("exact Python executable is unavailable")
    script_records: list[dict[str, object]] = []
    for script in args.python_script:
        script = script.resolve()
        if not script.is_file():
            errors.append(f"Python script is missing: {script}")
            continue
        ok, message = run([python_path, "-m", "py_compile", str(script)])
        if not ok:
            errors.append(f"Python cannot compile {script.name}: {message}")
        script_records.append({
            "path": str(script), "sha256": sha256(script), "compiled": ok,
        })
    import_records: list[dict[str, object]] = []
    for module in sorted(set(args.python_import)):
        ok, message = run([
            python_path, "-c",
            f"import {module}; print(getattr({module}, '__version__', 'present'))",
        ])
        if not ok:
            errors.append(f"exact Python lacks {module}: {message}")
        import_records.append({"module": module, "available": ok, "version": message})

    rscript_path = ""
    r_version = ""
    r_package_records: list[dict[str, object]] = []
    r_script_records: list[dict[str, object]] = []
    if args.rscript:
        rscript_path = shutil.which(args.rscript) or args.rscript
        ok, r_version = run([rscript_path, "--version"])
        if not ok:
            errors.append("exact Rscript executable is unavailable")
        for script in args.r_script:
            script = script.resolve()
            if not script.is_file():
                errors.append(f"R script is missing: {script}")
                continue
            ok, message = run([
                rscript_path, "-e", f"parse(file={json.dumps(str(script))})",
            ])
            if not ok:
                errors.append(f"R cannot parse {script.name}: {message}")
            r_script_records.append({
                "path": str(script), "sha256": sha256(script), "parsed": ok,
            })
        for package in sorted(set(args.r_package)):
            ok, message = run([
                rscript_path, "-e",
                (
                    f"if (!requireNamespace('{package}', quietly=TRUE)) "
                    f"quit(status=2); cat(as.character(packageVersion('{package}')))"
                ),
            ])
            if not ok:
                errors.append(f"exact R runtime lacks {package}: {message}")
            r_package_records.append({
                "package": package, "available": ok, "version": message,
            })

    input_records: list[dict[str, object]] = []
    for specification in args.semantic_input:
        try:
            input_records.append(validate_semantic_input(specification))
        except (ValueError, json.JSONDecodeError, UnicodeDecodeError) as exc:
            errors.append(str(exc))

    args.out.mkdir(parents=True, exist_ok=True)
    manifest = {
        "schema_version": "2.2",
        "status": "PASS" if not errors else "BLOCKED",
        "artifact_role": "phase_runtime_preflight",
        "python": {"path": str(Path(python_path).resolve()), "version": python_version},
        "python_scripts": script_records,
        "python_imports": import_records,
        "rscript": {"path": str(Path(rscript_path).resolve()), "version": r_version}
        if rscript_path else None,
        "r_packages": r_package_records,
        "r_scripts": r_script_records,
        "semantic_inputs": input_records,
        "errors": errors,
    }
    manifest_path = args.out / "phase_runtime_preflight.json"
    manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(manifest, ensure_ascii=False, indent=2))
    return 0 if not errors else 2


if __name__ == "__main__":
    raise SystemExit(main())
