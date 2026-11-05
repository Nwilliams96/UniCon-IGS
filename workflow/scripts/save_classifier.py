"""Validate and atomically save a completed SILVA or PR2 classifier.

The scientific construction steps remain visible as individual Snakemake
rules. This helper is intentionally limited to the infrastructure operation
that independent Snakemake processes cannot coordinate themselves: locking a
shared destination and replacing it only with a fully validated artifact.
"""


import fcntl
import os
import shutil
import subprocess
import re
import zipfile
import hashlib
from pathlib import Path


def classifier_is_valid(path, log_handle):
    path = Path(path)
    if not path.is_file() or path.stat().st_size == 0:
        return False
    try:
        with zipfile.ZipFile(path) as archive:
            metadata = [n for n in archive.namelist()
                        if len(n.split("/")) == 2 and n.endswith("/metadata.yaml")]
            if len(metadata) != 1 or not re.search(
                r"^type:\s*TaxonomicClassifier\s*$",
                archive.read(metadata[0]).decode("utf-8"), re.M
            ):
                return False
    except (OSError, ValueError, zipfile.BadZipFile, UnicodeError):
        return False
    result = subprocess.run(
        ["qiime", "tools", "validate", str(path)],
        stdout=log_handle,
        stderr=subprocess.STDOUT,
    )
    return result.returncode == 0


def write_marker(marker, classifier):
    marker = Path(marker)
    marker.parent.mkdir(parents=True, exist_ok=True)
    temporary = marker.with_name(marker.name + f".tmp.{os.getpid()}")
    temporary.write_text(
        f"classifier\t{Path(classifier).resolve()}\n"
        f"size_bytes\t{Path(classifier).stat().st_size}\n",
        encoding="utf-8",
    )
    os.replace(temporary, marker)


def save_classifier(source, destination, lock_path, marker, log_path,
                       read_only=False, label="PR2", expected_md5=None):
    source = Path(source)
    destination = Path(destination)
    lock_path = Path(lock_path)
    marker = Path(marker)
    log_path = Path(log_path)

    if not read_only:
        destination.parent.mkdir(parents=True, exist_ok=True)
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.parent.mkdir(parents=True, exist_ok=True)

    with log_path.open("a", encoding="utf-8") as log_handle:
        # A user-provided classifier is read-only. Only its project-local
        # marker/log are written; do not place locks beside that source file.
        if read_only:
            if source.resolve() != destination.resolve() or not classifier_is_valid(source, log_handle):
                raise RuntimeError(f"Configured pre-existing {label} classifier is invalid: {source}")
            write_marker(marker, source)
            return
        if source.resolve() != destination.resolve() and not classifier_is_valid(
            source, log_handle
        ):
            raise RuntimeError(f"Staged {label} classifier failed QIIME validation: {source}")

        log_handle.write(f"Waiting for shared {label} save lock: {lock_path}\n")
        log_handle.flush()
        with lock_path.open("a+b") as lock_handle:
            fcntl.flock(lock_handle.fileno(), fcntl.LOCK_EX)
            log_handle.write(f"Acquired shared {label} save lock\n")
            log_handle.flush()

            checksum_agrees = True
            if expected_md5 and destination.is_file():
                digest = hashlib.md5()
                with destination.open("rb") as cached:
                    for chunk in iter(lambda: cached.read(1024 * 1024), b""):
                        digest.update(chunk)
                checksum_agrees = digest.hexdigest().lower() == expected_md5.lower()
            if checksum_agrees and classifier_is_valid(destination, log_handle):
                log_handle.write(f"Reusing validated {label} classifier: {destination}\n")
                write_marker(marker, destination)
                return

            if source.resolve() == destination.resolve():
                raise RuntimeError(
                    f"The configured pre-existing {label} classifier is invalid: "
                    f"{destination}"
                )

            temporary = destination.with_name(
                destination.name + f".saving.{os.getpid()}"
            )
            try:
                with source.open("rb") as input_handle, temporary.open("wb") as output:
                    shutil.copyfileobj(input_handle, output)
                    output.flush()
                    os.fsync(output.fileno())
                if not classifier_is_valid(temporary, log_handle):
                    raise RuntimeError(
                        f"Copied {label} classifier failed validation before saving"
                    )
                if expected_md5:
                    digest = hashlib.md5()
                    with temporary.open("rb") as copied:
                        for chunk in iter(lambda: copied.read(1024 * 1024), b""):
                            digest.update(chunk)
                    if digest.hexdigest().lower() != expected_md5.lower():
                        raise RuntimeError(f"Copied {label} classifier failed MD5 verification before saving")
                os.replace(temporary, destination)
            finally:
                temporary.unlink(missing_ok=True)

            if not classifier_is_valid(destination, log_handle):
                raise RuntimeError(f"Saved {label} classifier failed QIIME validation")
            write_marker(marker, destination)
            log_handle.write(f"Saved {label} classifier atomically: {destination}\n")


if "snakemake" in globals():
    save_classifier(
        snakemake.input.classifier,
        snakemake.params.classifier,
        snakemake.params.lock,
        snakemake.output.marker,
        snakemake.log[0],
        read_only=getattr(snakemake.params, "read_only", False),
        label=getattr(snakemake.params, "label", "PR2"),
        expected_md5=getattr(snakemake.params, "md5", None),
    )
