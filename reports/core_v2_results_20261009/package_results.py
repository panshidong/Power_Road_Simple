#!/usr/bin/env python3
"""Publishable, lossless result bundle; reads the completed run without solvers."""
from pathlib import Path
import csv
import hashlib
import io
import json
import shutil
import tarfile

DEST = Path(__file__).resolve().parent
ROOT = DEST.parents[1]
RUN = ROOT / "runtime/output/core-v2-stallfix-20261008"
ARCHIVE = DEST / "Austin_model_v2_complete_results_20261009.tar.xz"


def sha(path):
    with Path(path).open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n")


def regular(path):
    return (path.is_file() and not path.is_symlink() and "__pycache__" not in path.parts
            and path.suffix not in (".pyc", ".lock", ".tmp", ".o", ".d"))


def main():
    state = json.loads((RUN / "control/status.json").read_text())
    audit = json.loads((RUN / "control/completion_acceptance.json").read_text())
    analysis = json.loads((RUN / "analysis/status.json").read_text())
    assert state["status"] == "completed" and audit["passed"] and analysis["complete"]
    assert audit["total_tasks"] == 336 and audit["evaluation_trajectories"] == 264
    assert audit["incomplete_repairs"] == 0 and not analysis["failures"] and not analysis["missing_jobs"]
    assert len(audit["source_sha256"]) == 336
    for name, expected in audit["source_sha256"].items():
        assert sha(RUN / name) == expected, name
    assert not list((RUN / "results").rglob("*.error.json"))
    with (RUN / "analysis/scenario_rows.csv").open() as handle:
        assert len(list(csv.DictReader(handle))) == 264

    # Every non-lock run artifact, including raw events and native diagnostics.
    selected = {}
    for p in sorted(RUN.rglob("*")):
        rel = p.relative_to(RUN)
        if regular(p) and rel.parts[0] != "snapshot":
            selected["original_run/" + rel.as_posix()] = p

    # Exact executed Python/config/native objects. Bulky public input copies are
    # indexed below instead of being duplicated inside the results download.
    frozen = RUN / "snapshot/runtime"
    for directory in ("austin_runtime", "configs", "build"):
        for p in sorted((frozen / directory).rglob("*")):
            if regular(p):
                selected["code_snapshot/runtime/" + p.relative_to(frozen).as_posix()] = p
    selected["code_snapshot/runtime/background_run.py"] = frozen / "background_run.py"
    selected["code_snapshot/runtime/audit_completed_run.py"] = RUN / "control/audit_completed_run.py"
    for name in ("catalog.json", "load_catalog.csv.gz"):
        selected["prepared_metadata/" + name] = frozen / "prepared" / name
    for directory in ("src", "include"):
        for p in sorted((ROOT / "vendor/tap-b" / directory).rglob("*")):
            if regular(p):
                selected["code_snapshot/vendor/tap-b/" + p.relative_to(ROOT / "vendor/tap-b").as_posix()] = p
    for name in ("Makefile", "LICENSE.md", "SOURCE.json"):
        selected["code_snapshot/vendor/tap-b/" + name] = ROOT / "vendor/tap-b" / name
    vendor = json.loads((ROOT / "vendor/tap-b/SOURCE.json").read_text())
    for name, expected in vendor["source_sha256"].items():
        if name.startswith(("src/", "include/")):
            assert sha(ROOT / "vendor/tap-b" / name) == expected, name
    assert sha(frozen / "build/tap-b/bin/tap") == vendor["local_patches"][-1]["native_build_sha256"]
    for p in sorted((ROOT / "runtime/tests_native").rglob("*")):
        if regular(p):
            selected["code_snapshot/runtime/tests_native/" + p.relative_to(ROOT / "runtime/tests_native").as_posix()] = p
    for name in ("pyproject.toml", "requirements.txt", "MODEL_FIXES_20261007.md"):
        p = ROOT / "runtime" / name
        if p.exists():
            selected["code_snapshot/runtime/" + name] = p
    for name in ("sources.lock.json", "requirements.txt"):
        if (ROOT / name).exists():
            selected["code_snapshot/" + name] = ROOT / name

    provenance = DEST / "provenance"
    provenance.mkdir(exist_ok=True)
    report_names = [
        "model_v2_completion_20261009.json", "model_v2_completion_process_check_20261009.json",
        "tapb_stagnation_fix_20261008.json", "tapb_stagnation_full_network_audit_20261008.json",
        "tapb_stagnation_20261008.patch", "tapb_stagnation_native_tests_20261008.log",
        "tapb_stagnation_unit_tests_20261008.log", "model_v2_local_power_native_20261007.log",
    ]
    for name in report_names:
        src = ROOT / "reports" / name
        shutil.copy2(src, provenance / name)
        selected["provenance/" + name] = provenance / name
    # The standalone before/after reproductions document the solver fix.
    diagnostics = ROOT / "runtime/diagnostics/tapb-stall-20261008"
    for case in ("original-replay", "no-aec-replay", "stallfix-replay"):
        for name in ("network.tntp", "tapb.log", "process.json"):
            p = diagnostics / case / name
            selected["solver_fix_reproduction/" + case + "/" + name] = p

    inventory = []
    for directory in (RUN / "snapshot/data", frozen / "prepared"):
        for p in sorted(directory.rglob("*")):
            if regular(p):
                inventory.append(dict(path=p.relative_to(RUN / "snapshot").as_posix(), bytes=p.stat().st_size, sha256=sha(p)))
    write(provenance / "input_inventory.json", dict(
        files=inventory, source_commit="c625594ebf15e54b8d2151111d73a46309f0b91b",
        repository="panshidong/Power_Road_Simple", branch="austin-model-fixes",
        public_input_sources="See code_snapshot/sources.lock.json; full data and circuit copies remain local.",
    ))
    selected["provenance/input_inventory.json"] = provenance / "input_inventory.json"

    # Browser-friendly copies alongside the complete lossless download.
    shutil.copytree(RUN / "analysis", DEST / "analysis", dirs_exist_ok=True)
    meta = DEST / "run_metadata"; meta.mkdir(exist_ok=True)
    for name in ("run_manifest.json", "validation.json", "tables.json", "RUN_SCOPE.txt"):
        shutil.copy2(RUN / name, meta / name)
    for name in ("status.json", "launch.json", "completion_acceptance.json", "health_check.json", "supervisor.log"):
        shutil.copy2(RUN / "control" / name, meta / name)

    hashes = {name: sha(path) for name, path in sorted(selected.items())}
    sums = "".join(f"{value}  {name}\n" for name, value in hashes.items()).encode()
    uncompressed = sum(p.stat().st_size for p in selected.values()) + len(sums)
    print(f"Archiving {len(selected)} original files ({uncompressed / 2**20:.1f} MiB)", flush=True)
    with tarfile.open(ARCHIVE, "w:xz", preset=6, format=tarfile.PAX_FORMAT) as tar:
        for name, path in sorted(selected.items()):
            info = tar.gettarinfo(str(path), arcname=name)
            info.uid = info.gid = 0; info.uname = info.gname = ""; info.mtime = 0
            with path.open("rb") as handle:
                tar.addfile(info, handle)
        info = tarfile.TarInfo("SHA256SUMS.txt"); info.size = len(sums); info.mtime = 0
        tar.addfile(info, io.BytesIO(sums))
    print(f"Archive: {ARCHIVE.stat().st_size / 2**20:.1f} MiB; verifying every member", flush=True)
    seen = set()
    with tarfile.open(ARCHIVE, "r:xz") as tar:
        for member in tar:
            assert member.isfile() and member.name not in seen
            seen.add(member.name)
            data = tar.extractfile(member)
            if member.name == "SHA256SUMS.txt":
                assert data.read() == sums
            else:
                assert hashlib.file_digest(data, "sha256").hexdigest() == hashes[member.name], member.name
    assert seen == set(selected) | {"SHA256SUMS.txt"}
    result_paths = [name for name in selected if name.startswith("original_run/results/")]
    assert len(result_paths) == 336
    write(DEST / "archive_verification.json", dict(
        archive=ARCHIVE.name, bytes=ARCHIVE.stat().st_size, sha256=sha(ARCHIVE),
        members=len(seen), original_members=len(selected), uncompressed_bytes=uncompressed,
        all_members_byte_identical_to_source=True, original_results=336,
        result_hashes_match_completion_acceptance=True, evaluation_trajectories=264,
        unfinished_repairs=0, final_audit_passed=True, fingerprint=audit["fingerprint"],
        solvers_started=0, excluded="Mutable lock files, Python bytecode, bulk public input copies, and reusable solver caches; all remain local.",
    ))
    print(json.dumps(json.loads((DEST / "archive_verification.json").read_text()), indent=2), flush=True)


if __name__ == "__main__":
    main()
