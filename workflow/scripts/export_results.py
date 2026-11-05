"""Copy the small, user-facing deliverables into Results-Export."""

import csv
from collections import defaultdict
from pathlib import Path
import shutil


def write_taxonomy_summary(long_data_path, output_path, rank):
    totals = defaultdict(float)
    sample_totals = defaultdict(float)
    with Path(long_data_path).open(encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            sample = (row.get("SampleID") or "").strip()
            taxon = (row.get(rank) or "Unassigned").strip() or "Unassigned"
            try:
                abundance = float(row.get("Corrected_Sequence_Counts") or 0)
            except ValueError:
                abundance = 0
            if sample and abundance > 0:
                totals[(sample, taxon)] += abundance
                sample_totals[sample] += abundance

    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["SampleID", rank, "Corrected_Sequence_Counts", "Relative_Abundance"])
        for (sample, taxon), abundance in sorted(totals.items()):
            relative = abundance / sample_totals[sample] if sample_totals[sample] else 0
            writer.writerow([sample, taxon, f"{abundance:.10g}", f"{relative:.10g}"])


def export_results(sources, destinations, readme_path, project_name, summaries=None):
    if len(sources) != len(destinations):
        raise ValueError("Every export source must have one destination")

    copied_names = []
    for source, destination in zip(sources, destinations):
        source_path = Path(source)
        destination_path = Path(destination)
        destination_path.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source_path, destination_path)
        copied_names.append(destination_path.name)

    if summaries:
        long_data = next(
            source
            for source, destination in zip(sources, destinations)
            if str(destination).endswith(".long_data.tsv")
            and ".ISD_corrected_" not in str(destination)
        )
        for rank, destination in zip(("Phylum", "Order"), summaries):
            write_taxonomy_summary(long_data, destination, rank)
            copied_names.append(Path(destination).name)

    readme = Path(readme_path)
    readme.parent.mkdir(parents=True, exist_ok=True)
    contents = [
        f"Results export for {project_name}",
        "",
        "This folder contains the small, user-facing outputs from the completed pipeline run.",
        "The HTML report is self-contained and can be opened directly in a web browser.",
        "",
        "Files:",
        *[f"- {name}" for name in copied_names],
        "- README.txt",
        "",
        "Keep the complete results/ directory if detailed QIIME 2 artifacts or intermediate files are needed.",
        "",
    ]
    readme.write_text("\n".join(contents), encoding="utf-8")


if "snakemake" in globals():
    export_results(
        list(snakemake.input.data_files),
        list(snakemake.output.data_files),
        str(snakemake.output.readme),
        str(snakemake.config["projectName"]),
        list(snakemake.output.summaries),
    )
