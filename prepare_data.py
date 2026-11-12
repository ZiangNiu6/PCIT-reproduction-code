#!/usr/bin/env python3
"""Fetch the immutable public inputs; never start an analysis or replace inputs."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import tempfile
import urllib.request
import zipfile

SOURCES = {
    "ccle": [
        {
            "name": "mutation.txt",
            "url": "https://data.broadinstitute.org/ccle_legacy_data/binary_calls_for_copy_number_and_mutation_data/CCLE_MUT_CNA_AMP_DEL_binary_Revealer.gct",
            "sha256": "937430ab6ce71e0a71bba581a6ccc41b797d45d6292830c0787faf3efef6b145",
            "strip_gct_headers": True,
        },
        {
            "name": "response.csv",
            "url": "https://data.broadinstitute.org/ccle_legacy_data/pharmacological_profiling/CCLE_NP24.2009_Drug_data_2015.02.24.csv",
            "sha256": "9dce1a7d4c27627bcb3b5c2e6d14873955aac2906e5d44ad344c265a491404e5",
        },
    ],
    "frangieh": [{
        "name": "frangieh-2021-processed.zip",
        "url": "https://zenodo.org/records/10976417/files/frangieh-2021-processed.zip?download=1",
        "md5": "1ad2adc53fde94fa17cad9ac6ebca006",
    }],
}
FRANGIEH_FILES = [
    "processed/control/gene/gene_expression_matrix.odm",
    "processed/control/gene/gene_expression_metadata.rds",
    "processed/control/grna_assignment/grna_assignments_ungrouped.odm",
    "processed/control/grna_assignment/grna_assignments_ungrouped_metadata.rds",
]


def checksum(path, algorithm="sha256"):
    result = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            result.update(chunk)
    return result.hexdigest()


def verify_download(path, source):
    for algorithm in ("sha256", "md5"):
        if algorithm in source and checksum(path, algorithm) != source[algorithm]:
            raise ValueError(f"{algorithm} mismatch for {source['name']}; input was not published")


def stage_download(source, directory, source_directory):
    target = directory / source["name"]
    if source_directory:
        # --source-dir accepts original public download filenames, not renamed
        # or transformed files. This also supports an offline transfer.
        name = source["url"].split("?")[0].rsplit("/", 1)[1]
        original = source_directory / name
        if not original.is_file():
            raise FileNotFoundError(original)
        shutil.copyfile(original, target)
    else:
        print(f"Downloading {source['name']} ...", flush=True)
        with urllib.request.urlopen(source["url"], timeout=120) as response, target.open("wb") as output:
            shutil.copyfileobj(response, output, 1024 * 1024)
    verify_download(target, source)
    return target


def verify_receipt(destination, domain):
    receipt_path = destination / "input-provenance.json"
    if not receipt_path.is_file():
        raise ValueError(f"No input-provenance.json in {destination}; prepare inputs first")
    receipt = json.loads(receipt_path.read_text())
    if receipt.get("domain") != domain or receipt.get("schema") != "pcit_public_inputs_v1":
        raise ValueError("Input provenance schema/domain mismatch")
    if receipt.get("sources") != SOURCES[domain]:
        raise ValueError("Input provenance does not identify the pinned public sources")
    expected = {"mutation.txt", "response.csv"} if domain == "ccle" else set(FRANGIEH_FILES)
    records = receipt.get("files", [])
    if {record.get("path") for record in records} != expected or len(records) != len(expected):
        raise ValueError("Input provenance file set mismatch")
    for record in records:
        path = destination / record["path"]
        if not path.is_file() or path.stat().st_size != record["bytes"] or checksum(path) != record["sha256"]:
            raise ValueError(f"Missing or modified input: {path}")
    print(f"Verified {len(records)} {domain} input files in {destination}")


def prepare(domain, destination, source_directory):
    if destination.exists():
        # A valid existing receipt is a resumable success; incomplete or altered
        # inputs are never silently replaced.
        verify_receipt(destination, domain)
        return
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".pcit-inputs-", dir=destination.parent) as temporary:
        staging = Path(temporary)
        output = staging / "inputs"
        output.mkdir()
        if domain == "ccle":
            for source in SOURCES[domain]:
                original = stage_download(source, staging, source_directory)
                if source.get("strip_gct_headers"):
                    with original.open("rb") as reader, (output / source["name"]).open("wb") as writer:
                        if reader.readline().strip() != b"#1.2" or reader.readline().strip() != b"48270\t1030":
                            raise ValueError("Unexpected CCLE GCT version/dimensions")
                        shutil.copyfileobj(reader, writer)
                else:
                    shutil.move(original, output / source["name"])
            paths = ["mutation.txt", "response.csv"]
        else:
            original = stage_download(SOURCES[domain][0], staging, source_directory)
            with zipfile.ZipFile(original) as archive:
                for relative in FRANGIEH_FILES:
                    member = "frangieh-2021-processed/" + relative
                    path = output / relative
                    path.parent.mkdir(parents=True, exist_ok=True)
                    with archive.open(member) as reader, path.open("wb") as writer:
                        shutil.copyfileobj(reader, writer, 1024 * 1024)
            paths = FRANGIEH_FILES
        receipt = {
            "schema": "pcit_public_inputs_v1", "domain": domain,
            "sources": SOURCES[domain],
            "transformation": "remove exactly the first two GCT lines" if domain == "ccle" else "extract four control-condition files; no cell/guide filtering",
            "files": [{"path": name, "bytes": (output / name).stat().st_size,
                       "sha256": checksum(output / name)} for name in paths],
        }
        (output / "input-provenance.json").write_text(json.dumps(receipt, indent=2) + "\n")
        if destination.exists():
            raise FileExistsError(f"Destination appeared during download: {destination}")
        os.rename(output, destination)
    verify_receipt(destination, domain)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("domain", choices=SOURCES)
    parser.add_argument("--destination", type=Path, required=True,
                        help="new CCLE directory or Frangieh root containing processed/control")
    parser.add_argument("--source-dir", type=Path, help="directory of original downloads for offline preparation")
    parser.add_argument("--check", action="store_true", help="verify a previously prepared directory without downloading")
    args = parser.parse_args()
    try:
        if args.check:
            verify_receipt(args.destination, args.domain)
        else:
            prepare(args.domain, args.destination, args.source_dir)
    except (OSError, ValueError, KeyError, zipfile.BadZipFile) as error:
        parser.exit(1, f"Data preparation failed: {error}\n")


if __name__ == "__main__":
    main()
