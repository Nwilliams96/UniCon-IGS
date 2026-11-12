import base64
import csv
import importlib.util
import json
import re
import tempfile
import unittest
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path
from types import SimpleNamespace
from unittest import mock


SCRIPT = Path(__file__).parents[1] / "workflow" / "scripts" / "generate_pipeline_report.py"
SPEC = importlib.util.spec_from_file_location("generate_pipeline_report", SCRIPT)
REPORT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(REPORT)


def render_interaction_fixture(root):
    """Write a synthetic, offline report for Python and optional browser tests."""
    root = Path(root)
    root.mkdir(parents=True, exist_ok=True)
    sample_rows = [
        ("001", "Environmental", "-39.5", "100.5", "20.5"),
        ("002", "Environmental", "-5.25", "2.5", "21"),
        ("010", "Mock", "", "", ""),
        ("020", "", "", "", ""),
        ("030", "Environmental", "12.25", "2.05", "20.5"),
    ]
    paths = {}

    def write_table(name, headers, rows, filename=None):
        destination = root / (filename or f"{name}.tsv")
        with destination.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(headers)
            writer.writerows(rows)
        paths[name] = destination

    write_table("samples", ("sample", "condition", "latitude", "depth", "Temperature"),
                sample_rows, "custom-samples.tsv")
    write_table("split_summary", ("sample", "prok_seqs_split", "euk_seqs_split"),
                ((row[0], 1200, 300) for row in sample_rows))
    write_table("stats16s", ("sample-id", "input", "filtered", "denoised", "merged", "non-chimeric"),
                ((row[0], 1200, 1100, 1000, 900, 800) for row in sample_rows))
    write_table("stats18s", ("sample-id", "input", "filtered", "denoised", "non-chimeric"),
                ((row[0], 300, 260, 220, 180) for row in sample_rows))
    long_rows = []
    for sample, *_ in sample_rows:
        for index in range(1, 41):
            # Removing the mock reverses this taxon's abundance rank, testing colour stability.
            abundance = 5000 if sample == "010" and index == 40 else 41 - index
            long_rows.append((sample, f"ASV-{index}", "Bacteria", "SILVA kingdom",
                              "", "", "", f"Taxon {index:02d}", "",
                              "Prokaryotic_16S", abundance))
        long_rows.append((sample, "PR2-ASV", "Eukaryota", "", "PR2 supergroup",
                          "PR2 division", "PR2 subdivision", "", "",
                          "Eukaryote_18S", 10))
    write_table("long_data", ("SampleID", "ASV_hash", "Domain", "Kingdom", "Supergroup",
                             "Division", "Subdivision", "Phylum", "ProPortal_ASV_Ecotype",
                             "Sequence_Type", "Corrected_Sequence_Counts"), long_rows)
    paths["cutadapt_qc"] = []
    for sample, *_ in sample_rows:
        path = root / f"{sample}.qc.txt"
        path.write_text("Total read pairs processed: 2,000\n"
                        "Pairs written (passing filters): 1,800 (90.0%)\n")
        paths["cutadapt_qc"].append(path)
    paths["logo"] = SCRIPT.parents[1] / "assets/unicon-logo.png"
    output = root / "report.html"
    REPORT.render_report({"projectName": "Report controls — synthetic fixture"}, paths, output)
    return output


class GeneratePipelineReportTests(unittest.TestCase):
    def test_read_bars_are_solid_forest_green_with_thin_black_borders(self):
        chart = REPORT.svg_sample_read_counts({"001": 120, "002": 60}, ["001", "002"], "Reads")
        bars = ET.fromstring("<root>" + chart + "</root>").findall(".//rect")
        self.assertEqual(len(bars), 2)
        self.assertEqual([bar.attrib["height"] for bar in bars], ["285.0", "142.5"])
        for bar in bars:
            self.assertEqual(bar.attrib["fill"], "#24543d")
            self.assertEqual(bar.attrib["stroke"], "#111111")
            self.assertEqual(bar.attrib["stroke-width"], "0.8")
            self.assertNotIn("opacity", bar.attrib)

    def test_domain_display_hides_other_labels_without_changing_denominators_or_inputs(self):
        values = {"001": {"Bacteria": 60, "Eukaryota": 20, "Archaea": 10,
                          "Unassigned": 5, "Eukaryota:nucl": 5}}
        original = json.loads(json.dumps(values))
        chart = REPORT.svg_composition(values, "Domains")
        root = ET.fromstring("<root>" + chart + "</root>")
        bars = root.findall(".//rect")
        self.assertEqual([bar.attrib["height"] for bar in bars], ["180.0", "60.0", "30.0"])
        self.assertEqual([bar.attrib["fill"] for bar in bars], list(REPORT.DOMAIN_COLOURS.values()))
        self.assertNotIn("Unassigned", chart)
        self.assertNotIn("Eukaryota:nucl", chart)
        self.assertIn("60.0%", chart)
        self.assertEqual(values, original)
        hidden_only = REPORT.svg_composition({"001": {"Unassigned": 10}}, "Domains")
        self.assertNotIn("<rect", hidden_only)
        self.assertIn("001", hidden_only)  # Keep the sample even when no domain is displayed.

    def test_assignment_colours_do_not_shift_when_other_categories_are_absent(self):
        chart = REPORT.svg_assignment_counts({"Eukaryotic 18S": 12, "Unassigned": 3})
        bars = ET.fromstring(chart).findall("rect")
        self.assertEqual([bar.attrib["fill"] for bar in bars], ["#346b8b", "#9b9c95"])
        self.assertEqual([bar.attrib["width"] for bar in bars], ["570.0", "142.5"])

    def test_quality_bands_keep_their_coordinates_and_use_forest_green(self):
        points = [{"position": 1, "9%": 20, "25%": 25, "50%": 30, "75%": 35, "91%": 40},
                  {"position": 100, "9%": 10, "25%": 15, "50%": 20, "75%": 25, "91%": 30}]
        root = ET.fromstring(REPORT.svg_quality_profile(points, "Quality"))
        polygons = root.findall("polygon")
        self.assertEqual(polygons[0].attrib["points"], "58.0,54.4 800.0,115.3 800.0,237.1 58.0,176.2")
        self.assertEqual([p.attrib["fill"] for p in polygons], ["#24543d", "#24543d"])
        self.assertEqual([p.attrib["opacity"] for p in polygons], [".16", ".40"])
        self.assertEqual(root.find("path").attrib["stroke"], "#24543d")

    def test_internal_standard_markers_remain_circles_with_consistent_domain_colours(self):
        chart = REPORT.svg_log_series(
            ["001"], {"001": {"Archaea": 10, "Bacteria": 100}},
            ["Archaea", "Bacteria"], "Domain copies", category_colours=REPORT.DOMAIN_COLOURS,
        )
        root = ET.fromstring("<root>" + chart + "</root>")
        circles = root.findall(".//circle")
        self.assertEqual([point.attrib["fill"] for point in circles], ["#d59c32", "#24543d"])
        self.assertEqual([point.find("title").text for point in circles],
                         ["001 — Archaea: 10", "001 — Bacteria: 100"])
        self.assertEqual(root.findall(".//rect"), [])
        self.assertEqual(root.findall(".//polygon"), [])

    def test_plot_background_opacity_and_domain_omissions_are_explicit(self):
        with tempfile.TemporaryDirectory() as tmp:
            rendered = self.render_branding_fixture(Path(tmp))
        self.assertIn("--plot-background:rgba(255,255,255,0.60)", rendered)
        self.assertIn("Only Bacteria, Eukaryota and Archaea are displayed", rendered)
        self.assertIn("Percentages use all assignments", rendered)
        self.assertIn("other assignments remain in the exported data", rendered)

    def render_branding_fixture(self, root, logo=None):
        """Render synthetic QC records without running sequence processing."""
        inputs = {
            "samples": "sample\n001\n",
            "split_summary": "sample\tprok_seqs_split\teuk_seqs_split\n001\t100\t20\n",
            "stats16s": "sample-id\tinput\tfiltered\tdenoised\tmerged\tnon-chimeric\n001\t100\t95\t90\t85\t80\n",
            "stats18s": "sample-id\tinput\tfiltered\tdenoised\tnon-chimeric\n001\t20\t18\t15\t10\n",
            "long_data": (
                "SampleID\tASV_hash\tDomain\tPhylum\tCorrected_Sequence_Counts\n"
                "001\tA1\tArchaea\tThermoproteota\t80\n"
                "001\tA2\tEukaryota\tDinoflagellata\t10\n"
            ),
        }
        paths = {}
        for name, content in inputs.items():
            paths[name] = root / f"{name}.tsv"
            paths[name].write_text(content, encoding="utf-8")
        if logo is not None:
            paths["logo"] = logo
        output = root / "report.html"
        REPORT.render_report({"projectName": "Branding fixture"}, paths, output)
        return output.read_text(encoding="utf-8")

    def test_branding_embeds_the_bundled_logo_for_offline_reports(self):
        logo = SCRIPT.parents[1] / "assets" / "unicon-logo.png"
        with tempfile.TemporaryDirectory() as tmp:
            rendered = self.render_branding_fixture(Path(tmp), logo)
        image = re.search(r'<img class="brand-logo" src="data:image/png;base64,([^\"]+)"', rendered)
        self.assertIsNotNone(image)
        self.assertEqual(base64.b64decode(image.group(1)), logo.read_bytes())
        self.assertTrue(logo.read_bytes().startswith(b"\x89PNG\r\n\x1a\n"))
        self.assertIn('<span>UniCon-IGS</span>', rendered)
        self.assertIn('<h1>Branding fixture</h1>', rendered)
        self.assertIn('alt="UniCon-IGS unicorn"', rendered)
        self.assertIn("width:56px;height:56px;object-fit:contain;flex:none", rendered)
        self.assertNotRegex(rendered, r'<(?:img|script)[^>]+src="https?://')
        self.assertNotIn("linear-gradient(130deg,#12213c,#1d4ed8)", rendered)

    def test_standalone_report_can_omit_a_logo_but_an_explicit_bad_path_fails(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            rendered = self.render_branding_fixture(root)
            self.assertIn('<span>UniCon-IGS</span>', rendered)
            self.assertNotIn('<img class="brand-logo"', rendered)
            with self.assertRaises(FileNotFoundError):
                self.render_branding_fixture(root, root / "missing-logo.png")

    def test_report_theme_matches_parchment_and_retains_accessible_text(self):
        with tempfile.TemporaryDirectory() as tmp:
            rendered = self.render_branding_fixture(Path(tmp))
        theme = re.search(r":root\{([^}]+)\}", rendered).group(1)
        colours = dict(re.findall(r"--([a-z-]+):\s*(#[0-9a-f]{6})", theme))
        self.assertEqual(colours, {
            "ink": "#101010", "muted": "#636363", "line": "#d6cbb6",
            "paper": "#fcf8ef", "wash": "#f3ead7", "accent": "#101010",
        })

        def luminance(colour):
            channels = [int(colour[index:index + 2], 16) / 255 for index in (1, 3, 5)]
            linear = [value / 12.92 if value <= 0.04045 else ((value + 0.055) / 1.055) ** 2.4
                      for value in channels]
            return sum(value * weight for value, weight in zip(linear, (0.2126, 0.7152, 0.0722)))

        for background in ("paper", "wash"):
            for text in ("ink", "muted"):
                ratio = (luminance(colours[background]) + 0.05) / (luminance(colours[text]) + 0.05)
                self.assertGreaterEqual(ratio, 4.5, (background, text, ratio))
        self.assertEqual(REPORT.PALETTE, [
            "#332288", "#88CCEE", "#44AA99", "#117733", "#999933",
            "#DDCC77", "#CC6677", "#882255", "#AA4499", "#0072B2",
            "#E69F00", "#009E73", "#D55E00", "#56B4E9", "#64748B",
        ])

    def test_snakemake_tracks_and_passes_the_bundled_logo(self):
        rule = (SCRIPT.parents[1] / "rules" / "06-generate-report.smk").read_text()
        self.assertIn('logo=workflow.source_path("../assets/unicon-logo.png")', rule)
        scalar_names = (
            "logo", "samples", "split_summary", "stats16s", "stats18s",
            "correction_factors", "species_assignment", "chloroplast_audit",
            "chloroplast_summary", "quality_16s", "quality_18s_paired",
            "quality_18s_concatenated", "long_data",
        )
        inputs = {name: Path("inputs") / name for name in scalar_names}
        inputs["cutadapt_qc"] = []
        snakemake = SimpleNamespace(
            input=SimpleNamespace(**inputs), config={"projectName": "Test"},
            output=SimpleNamespace(html=Path("report.html")),
        )
        with mock.patch.object(REPORT, "render_report") as renderer:
            REPORT.run_from_snakemake(snakemake)
        self.assertEqual(renderer.call_args.args[1]["logo"], str(inputs["logo"]))

    def test_report_tracks_the_configured_sample_sheet_not_a_hardcoded_copy(self):
        rule = (SCRIPT.parents[1] / "rules/06-generate-report.smk").read_text()
        self.assertIn('samples=config["samplesheet"]', rule)
        self.assertNotIn('samples="config/samples.tsv"', rule)

    def test_database_rank_labels_do_not_rename_taxonomy_or_drop_missing_values(self):
        with tempfile.TemporaryDirectory() as tmp:
            rendered = render_interaction_fixture(tmp).read_text()
        payload = re.search(r'id="taxonomy-explorer-data">(.*?)</script>', rendered, re.S)
        data = json.loads(payload.group(1))
        self.assertEqual(data["rankLabels"]["Phylum"], "Phylum (SILVA only)")
        self.assertEqual(data["rankLabels"]["Kingdom"], "Kingdom (SILVA only)")
        for rank in ("Supergroup", "Division", "Subdivision"):
            self.assertEqual(data["rankLabels"][rank], f"{rank} (PR2 only)")
        self.assertEqual(data["rankLabels"]["Domain"], "Domain")
        self.assertIn("ProPortal only", data["rankLabels"]["ProPortal_ASV_Ecotype"])
        self.assertEqual(data["samples"]["001"]["metadata"]["SampleID"], "001")
        self.assertEqual(data["samples"]["001"]["metadata"]["Latitude"], "-39.5")
        self.assertEqual(data["samples"]["010"]["metadata"]["Latitude"], "")
        self.assertEqual(data["samples"]["001"]["metadata"]["Temperature"], "20.5")
        self.assertEqual(len(data["samples"]["001"]["taxonomy"]["Phylum"]), 41)
        self.assertEqual(data["samples"]["001"]["taxonomy"]["Phylum"]["Taxon 40"], 1)
        self.assertEqual(data["samples"]["001"]["taxonomy"]["Phylum"]["Unassigned at Phylum"], 10)
        self.assertEqual(sum(data["samples"]["001"]["taxonomy"]["Phylum"].values()), 830)
        self.assertNotIn("Other", data["samples"]["001"]["taxonomy"]["Phylum"])

    def test_pooled_losses_remain_visible_and_sample_tables_start_collapsed(self):
        with tempfile.TemporaryDirectory() as tmp:
            rendered = render_interaction_fixture(tmp).read_text()
        quality = re.search(r'<section id="quality">(.*?)</section>', rendered, re.S).group(1)
        bbsplit = re.search(r'<details id="bbsplit-sample-details"[^>]*>(.*?)</details>', quality, re.S)
        dada2 = re.search(r'<details id="dada2-sample-details"[^>]*>(.*?)</details>', quality, re.S)
        for match in (bbsplit, dada2):
            self.assertIsNotNone(match)
            self.assertNotIn(" open", match.group(0).split(">", 1)[0])
        self.assertEqual(bbsplit.group(1).count("<tr>"), 6)  # Header plus five samples.
        self.assertEqual(dada2.group(1).count("<tr>"), 11)  # Header plus two paths/sample.
        visible = quality.replace(bbsplit.group(0), "").replace(dada2.group(0), "")
        self.assertIn("Whole-run BBsplit summary", visible)
        self.assertIn("Whole-run DADA2 summary", visible)
        self.assertIn("<th>All samples</th><td>10,000</td><td>1,000 (10.0%)</td>", visible)
        self.assertIn("1,500 (16.7%)", visible)  # BBsplit loss from 9,000 trimmed pairs.
        self.assertIn("<td>16S paired</td><td>6,000</td>", visible)
        self.assertIn("<td>18S concatenated</td><td>1,500</td>", visible)
        self.assertIn("66.7%", visible)
        self.assertIn("60.0%", visible)

    def test_static_fallback_keeps_all_taxa_in_bounded_panels(self):
        totals = REPORT.Counter({f"Taxon-{index}": index + 1 for index in range(45)})
        self.assertEqual(REPORT.svg_top_taxa(totals).count("<rect"), 45)
        self.assertEqual(REPORT.svg_top_taxa(totals, limit=12).count("<rect"), 12)
        with tempfile.TemporaryDirectory() as tmp:
            rendered = render_interaction_fixture(tmp).read_text()
        fallback = re.search(r"<noscript>(.*?)</noscript>", rendered, re.S).group(1)
        self.assertIn("Taxon 40", fallback)
        self.assertIn('class="taxonomy-static-chart"', fallback)
        self.assertIn("max-height:240px;overflow:auto", rendered)
        self.assertIn("min-height:0;flex-shrink:0", rendered)
        self.assertNotIn("min-height:1px", rendered)
        self.assertNotIn("rows.slice(0, 20)", REPORT.TAXONOMY_EXPLORER_JS)
        self.assertNotIn("rows.slice(0, palette.length)", REPORT.TAXONOMY_EXPLORER_JS)

    def test_qiime_environment_identity_and_name_are_in_parameter_rows(self):
        rows = dict(REPORT.flatten_config({"qiime2_environment": {
            "kind": "named", "requested": "qiime-2026.7",
            "specification": "/shared/envs/qiime", "identity": "named-abc123",
            "software_versions": {"scikit-learn": "1.7.1"},
        }}))
        self.assertEqual(rows["qiime2_environment.requested"], "qiime-2026.7")
        self.assertEqual(rows["qiime2_environment.specification"], "/shared/envs/qiime")
        self.assertEqual(rows["qiime2_environment.software_versions.scikit-learn"], "1.7.1")

    def test_database_sources_are_reported_as_declared_and_escaped(self):
        rendered = REPORT.database_resources_table({
            "database_resources": {
                "SILVA": {
                    "mode": "download", "version": "144",
                    "classifier": "/refs/silva.qza", "directory": "/refs",
                    "url": "https://example.org/silva.qza?a=1&b=<test>",
                    "checksum": "abc123", "classifier_region": "full_length",
                },
                "PR2": {
                    "mode": "train", "version": "5.1.1",
                    "classifier": "/refs/pr2.qza", "reference_mode": "amplicon",
                    "training_settings": {"Forward primer": "GTGYCAGCMGCCGCGGTAA",
                                          "Dereplication": "uniq"},
                },
            },
        })
        self.assertIn("Download pre-trained database", rendered)
        self.assertIn("Train database with RESCRIPt", rendered)
        self.assertIn("Declared release", rendered)
        self.assertIn("Full-length SSU reference", rendered)
        self.assertIn("Configured primer region", rendered)
        self.assertIn("Training — Forward primer", rendered)
        self.assertIn("GTGYCAGCMGCCGCGGTAA", rendered)
        self.assertIn("uniq", rendered)
        self.assertIn("a=1&amp;b=&lt;test&gt;", rendered)
        self.assertIn("abc123", rendered)
        self.assertIn("does not verify the reference release", rendered)

        existing = REPORT.database_resources_table({
            "database_resources": {
                "SILVA": {"mode": "existing", "version": "138.2", "classifier": "/old.qza"}
            }
        })
        self.assertIn("Use existing databases", existing)
        self.assertIn("138.2", existing)
        self.assertIn("/old.qza", existing)
        self.assertIn("As supplied by the classifier", existing)

    def test_skipped_species_matching_explains_release_without_zero_cards(self):
        rendered = REPORT.species_assignment_table([{
            "status": "skipped",
            "reason": "Exact species matching is currently supported only for SILVA 144",
            "silva_version": "138.2",
        }])
        self.assertIn("Species matching skipped", rendered)
        self.assertIn("supported only for SILVA 144", rendered)
        self.assertNotIn('class="cards"', rendered)
        self.assertNotIn("species labels added", rendered)

    def test_older_release_report_has_matching_release_and_skip_reason(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            inputs = {
                "samples": "sample\n001\n",
                "split_summary": "sample\tprok_seqs_split\teuk_seqs_split\n001\t10\t2\n",
                "stats16s": "sample-id\tinput\tfiltered\tdenoised\tmerged\tnon-chimeric\n001\t10\t9\t8\t7\t6\n",
                "stats18s": "sample-id\tinput\tfiltered\tdenoised\tnon-chimeric\n001\t2\t2\t2\t2\n",
                "long_data": "SampleID\tASV_hash\tDomain\tCorrected_Sequence_Counts\n001\tASV1\tBacteria\t6\n",
                "species_assignment": (
                    "status\treason\tsilva_version\n"
                    "skipped\tExact species matching is currently supported only for SILVA 144\t138.2\n"
                ),
            }
            paths = {}
            for name, content in inputs.items():
                paths[name] = root / f"{name}.tsv"
                paths[name].write_text(content, encoding="utf-8")
            output = root / "report.html"
            REPORT.render_report(
                {"projectName": "older-release", "SILVAversion": "138.2"},
                paths,
                output,
            )
            rendered = output.read_text(encoding="utf-8")
            self.assertIn("SILVA 138.2 species matching", rendered)
            self.assertIn("Species matching skipped", rendered)
            self.assertNotIn("<h3>SILVA 144 species matching</h3>", rendered)
            self.assertNotIn("species labels added", rendered)

    def test_compiles_after_snakemake_generated_preamble(self):
        source = "snakemake = None\n" + SCRIPT.read_text(encoding="utf-8")
        compile(source, str(SCRIPT), "exec")

    def test_non_finite_values_are_treated_as_missing(self):
        for value in ("NaN", "nan", "Inf", "-Inf", float("nan"), float("inf")):
            with self.subTest(value=value):
                self.assertIsNone(REPORT.number(value))
                self.assertEqual(REPORT.fmt_count(value), "—")
                self.assertEqual(REPORT.fmt_percent(value), "—")

        chart = REPORT.svg_log_series(
            ["S-1", "S-2"],
            {
                "S-1": {"Bacteria": float("nan"), "Eukaryota": 100},
                "S-2": {"Bacteria": float("inf"), "Eukaryota": 200},
            },
            ["Bacteria", "Eukaryota"],
            "Non-finite regression test",
        )
        self.assertIn("Non-finite regression test", chart)
        self.assertNotIn(">nan<", chart.lower())
        self.assertNotIn(">inf<", chart.lower())

    def test_log_axis_extends_below_smallest_positive_value(self):
        chart = REPORT.svg_log_series(
            ["S-1", "S-2"],
            {"S-1": {"Archaea": 100}, "S-2": {"Archaea": 10000}},
            ["Archaea"],
            "Absolute Archaea",
        )
        self.assertIn(">50.1</text>", chart)

    def test_renders_summary_and_optional_figure(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "samples.tsv").write_text(
                "sample\tcondition\tLatitude [degrees_north]\tLongitude [degrees_east]\tDepth (m)\n"
                "S_1\ttest\t-39.49\t99.39\t2\n",
                encoding="utf-8",
            )
            (root / "split.tsv").write_text(
                "sample\tprok_seqs_split\teuk_seqs_split\ttotal_seqs_split\teukfrac_split\n"
                "S_1\t700\t300\t1000\t0.3\n",
                encoding="utf-8",
            )
            stats = (
                "sample-id\tinput\tfiltered\tdenoised\tmerged\tnon-chimeric\n"
                "S-1\t{input}\t{filtered}\t{denoised}\t{merged}\t{final}\n"
            )
            (root / "16.tsv").write_text(
                stats.format(input=700, filtered=650, denoised=640, merged=610, final=600),
                encoding="utf-8",
            )
            (root / "18.tsv").write_text(
                stats.format(input=300, filtered=250, denoised=220, merged="", final=200),
                encoding="utf-8",
            )
            (root / "correction-factors.tsv").write_text(
                "path\tamplicon_amount_pM\tstarting_molar_fraction\t"
                "bbsplit_assigned_reads\tobserved_read_fraction\tcorrection_factor\n"
                "16S\t1000\t0.5\t700\t0.7\t0.7142857143\n"
                "18S\t1000\t0.5\t300\t0.3\t1.666666667\n",
                encoding="utf-8",
            )
            (root / "species-assignment.tsv").write_text(
                "total_16S_ASVs\teligible_SILVA_ASVs\t"
                "exact_unambiguous_species_matches\tspecies_labels_added\t"
                "exact_matches_withheld_by_genus_check\t"
                "unmatched_or_ambiguous\n"
                "5\t4\t3\t2\t1\t1\n",
                encoding="utf-8",
            )
            (root / "chloroplast-summary.tsv").write_text(
                "PR2_min_confidence\ttotal_16S_ASVs\tSILVA_chloroplast_ASVs\t"
                "PR2_plastid_ASVs_any_confidence\tPR2_plastid_ASVs_accepted\t"
                "PR2_confirmed_SILVA_ASVs\tPR2_rescued_ASVs\t"
                "SILVA_chloroplast_unconfirmed_ASVs\t"
                "PR2_low_confidence_not_accepted_ASVs\tSILVA_retained_ASVs\n"
                "0.7\t5\t1\t3\t2\t1\t1\t0\t1\t2\n",
                encoding="utf-8",
            )
            (root / "chloroplast-audit.tsv").write_text(
                "ASV_hash\tChloroplast_detection_source\n"
                "A3\tPR2_confirmed_SILVA\n"
                "A1\tPR2_rescue\n",
                encoding="utf-8",
            )
            (root / "unmerged-summary.tsv").write_text(
                "standard_non_chimeric_16S_reads\tretained_run_non_chimeric_16S_reads\t"
                "supplemental_unmerged_16S_features\tsupplemental_unmerged_16S_reads\t"
                "samples_with_supplemental_reads\tsupplemental_chloroplast_features\t"
                "supplemental_chloroplast_reads\tinterpretation\n"
                "600\t675\t3\t75\t1\t1\t20\tExperimental supplement\n",
                encoding="utf-8",
            )
            (root / "S_1.qc.txt").write_text(
                "Total read pairs processed:              1,250\n"
                "Pairs written (passing filters):          1,100 (88.0%)\n",
                encoding="utf-8",
            )
            (root / "long.tsv").write_text(
                "SampleID\tDomain\tKingdom\tPhylum\tDivision\tClass\tSequence_Type\t"
                "plastid_16S_rRNA\tASV_hash\tCorrected_Sequence_Counts\n"
                "S-1\tBacteria\tPseudomonadati\tPseudomonadota\t\t\tProkaryotic_16S\tno\tA1\t600\n"
                "S-1\tEukaryota\t\t\tDinoflagellata\t\tEukaryote_18S\tno\tA2\t200\n"
                "S-1\tEukaryota\t\t\tChloroplastida\t\tChloroplast_16S\tyes\tA3\t100\n"
                "S-1\tEukaryota\t\t\t\tMitochondria\tEukaryote_18S\tno\tA4\t50\n"
                "S-1\tUnassigned\t\t\t\t\tUnassigned\tno\tA5\t25\n"
                "S-1\tBacteria\t\tInternal-standard\t\t\tProkaryotic_16S\tno\tA-ISD\t1000\n",
                encoding="utf-8",
            )
            quality_summary = (
                "\t1\t2\t3\t4\n"
                "count\t100\t100\t100\t100\n"
                "2%\t10\t11\t12\t13\n"
                "9%\t20\t21\t22\t23\n"
                "25%\t30\t31\t32\t33\n"
                "50%\t35\t36\t37\t38\n"
                "75%\t40\t41\t42\t43\n"
                "91%\t42\t43\t44\t45\n"
                "98%\t45\t45\t45\t45\n"
            )
            for directory, paired in (
                (root / "quality-16s", True),
                (root / "quality-18s-paired", True),
                (root / "quality-18s-concat", False),
            ):
                directory.mkdir()
                with zipfile.ZipFile(directory / "visualization.qzv", "w") as qzv:
                    qzv.writestr(
                        "fixture/data/forward-seven-number-summaries.tsv",
                        quality_summary,
                    )
                    if paired:
                        qzv.writestr(
                            "fixture/data/reverse-seven-number-summaries.tsv",
                            quality_summary,
                        )
            # A tiny valid PNG is sufficient to verify base64 embedding.
            png = root / "figure.png"
            png.write_bytes(
                b"\x89PNG\r\n\x1a\n\x00\x00\x00\rIHDR\x00\x00\x00\x01\x00\x00\x00\x01"
                b"\x08\x06\x00\x00\x00\x1f\x15\xc4\x89"
            )
            (root / "corrected.tsv").write_text(
                "SampleID\tASV_hash\tDomain\tisd_1_recovery_ratio\trecovery_mean\t"
                "Copies_BP_recovery_ratio\n"
                "S-1\tA1\tBacteria\t0.25\t0.25\t2400\n"
                "S-1\tA2\tEukaryota\t0.25\t0.25\t400\n"
                "S-1\tA3\tEukaryota\t0.25\t0.25\t200\n"
                "S-1\tA4\tEukaryota\t0.25\t0.25\t100\n"
                "S-1\tA5\tUnassigned\t0.25\t0.25\t100\n",
                encoding="utf-8",
            )
            output = root / "report.html"
            REPORT.render_report(
                {
                    "projectName": "test-run",
                    "intstds": ["BP"],
                    "trunclens": {"truncR1": 220, "truncR2": 180},
                    "dada2": {
                        "prokaryotes": {"max_ee_f": 3.5},
                        "eukaryotes": {"max_ee": 4.5},
                    },
                },
                {
                    "samples": root / "samples.tsv",
                    "split_summary": root / "split.tsv",
                    "stats16s": root / "16.tsv",
                    "stats18s": root / "18.tsv",
                    "correction_factors": root / "correction-factors.tsv",
                    "species_assignment": root / "species-assignment.tsv",
                    "chloroplast_audit": root / "chloroplast-audit.tsv",
                    "chloroplast_summary": root / "chloroplast-summary.tsv",
                    "cutadapt_qc": [root / "S_1.qc.txt"],
                    "quality_16s": root / "quality-16s",
                    "quality_18s_paired": root / "quality-18s-paired",
                    "quality_18s_concatenated": root / "quality-18s-concat",
                    "long_data": root / "long.tsv",
                    "internal_standard_figures": [png],
                    "internal_standard_table": [root / "corrected.tsv"],
                    "unmerged_16s_summary": [root / "unmerged-summary.tsv"],
                },
                output,
            )
            rendered = output.read_text(encoding="utf-8")
            self.assertIn("test-run", rendered)
            self.assertIn("<title>test-run pipeline report</title>", rendered)
            self.assertIn("S-1", rendered)
            self.assertIn("Filtering loss", rendered)
            self.assertIn("Denoising loss", rendered)
            self.assertIn("Pair-merging loss", rendered)
            self.assertIn("Chimera-removal loss", rendered)
            self.assertIn("Reads before filtering and quality control", rendered)
            self.assertIn("S-1: 1,250", rendered)
            self.assertIn("Primer trimming and BBsplit assignment", rendered)
            self.assertIn("BBsplit unassigned", rendered)
            self.assertIn("150 (12.0%)", rendered)
            self.assertIn("100 (9.1%)", rendered)
            self.assertIn("Reads retained after DADA2", rendered)
            self.assertIn("QIIME 2 sequence-quality profiles", rendered)
            self.assertIn("16S forward reads (R1)", rendered)
            self.assertIn("16S reverse reads (R2)", rendered)
            self.assertIn("18S forward reads before concatenation", rendered)
            self.assertIn("18S reverse reads before concatenation", rendered)
            self.assertIn("Concatenated 18S reads entering DADA2", rendered)
            self.assertEqual(rendered.count('class="quality-card"'), 5)
            self.assertIn('class="quality-profile"', rendered)
            self.assertIn("S-1: 800", rendered)
            self.assertIn("Not applicable", rendered)
            self.assertIn("Pseudomonadota", rendered)
            self.assertIn("Interactive taxonomy bar plot", rendered)
            self.assertIn('id="taxonomy-explorer-data"', rendered)
            self.assertIn('id="taxonomy-plot-field"', rendered)
            self.assertIn('id="taxonomy-rank"', rendered)
            self.assertIn('"Kingdom":{"Pseudomonadati":600.0,', rendered)
            self.assertIn("justify-content:flex-start", rendered)
            self.assertIn("top:308px", rendered)
            self.assertIn('class="frozen-y-axis"', rendered)
            self.assertIn('axis.className = "taxonomy-y-axis"', rendered)
            self.assertIn("#332288", rendered)
            self.assertIn(
                '"metadataFields":["SampleID","Condition","Latitude","Longitude","Depth"]',
                rendered,
            )
            self.assertIn('"Phylum":{"Pseudomonadota":600.0,', rendered)
            self.assertIn("Sequence assignments", rendered)
            self.assertIn("total 16S", rendered)
            self.assertIn("750</strong><span>total 16S", rendered)
            self.assertIn("100</strong><span>chloroplast 16S", rendered)
            self.assertIn("50</strong><span>mitochondrial 16S", rendered)
            self.assertIn("25</strong><span>unassigned", rendered)
            self.assertIn("Sequence-assignment counts", rendered)
            self.assertIn("Recovery ratios by sample (log10 scale)", rendered)
            self.assertIn("Domain abundance by correction method", rendered)
            self.assertIn("S-1 — BP: 0.25", rendered)
            self.assertNotIn("A-ISD", rendered)
            self.assertIn("trunclens.truncR1", rendered)
            self.assertIn("Effective DADA2 settings used", rendered)
            self.assertIn("16S and 18S correction factors", rendered)
            self.assertIn("Correction factor used", rendered)
            self.assertIn("SILVA 144 species matching", rendered)
            self.assertIn("unambiguous exact matches", rendered)
            self.assertIn("2</strong><span>species labels added", rendered)
            self.assertIn("PR2 plastid cross-check", rendered)
            self.assertIn("plastid ASVs rescued by PR2", rendered)
            self.assertIn("600</strong><span>abundance assigned to rescued ASVs", rendered)
            self.assertIn("PR2 plastid calls accepted at confidence ≥ 0.7", rendered)
            self.assertIn("0.714286", rendered)
            self.assertIn("1.66667", rendered)
            self.assertIn("max_ee_f", rendered)
            self.assertIn("Maximum expected errors in a forward read", rendered)
            self.assertIn("<code>3.5</code>", rendered)
            self.assertIn("<code>4.5</code>", rendered)
            self.assertIn("min_overlap", rendered)
            self.assertIn("max_merge_mismatch", rendered)
            self.assertIn("Unmerged 16S linked-read rescue", rendered)
            self.assertIn("automatically runs a second DADA2 pass", rendered)
            self.assertIn("3</strong><span>supplemental linked features", rendered)
            self.assertIn("not</strong> included in the canonical ASV table", rendered)
            self.assertIn("n_reads_learn", rendered)
            self.assertNotIn("Median base-quality profiles", rendered)
            self.assertNotIn("Where reads were retained or lost", rendered)

            explorer = REPORT.build_taxonomy_explorer_data(
                REPORT.read_tsv(root / "samples.tsv"),
                REPORT.read_tsv(root / "long.tsv"),
                "Corrected_Sequence_Counts",
            )
            self.assertEqual(
                explorer["metadataFields"],
                ["SampleID", "Condition", "Latitude", "Longitude", "Depth"],
            )
            self.assertIn("Phylum", explorer["ranks"])
            self.assertIn("Kingdom", explorer["ranks"])
            self.assertEqual(explorer["samples"]["S-1"]["metadata"]["SampleID"], "S-1")
            self.assertEqual(explorer["samples"]["S-1"]["metadata"]["Condition"], "test")
            self.assertEqual(explorer["samples"]["S-1"]["metadata"]["Latitude"], "-39.49")
            self.assertEqual(
                explorer["samples"]["S-1"]["taxonomy"]["Phylum"]["Pseudomonadota"],
                600,
            )
            self.assertEqual(
                explorer["samples"]["S-1"]["taxonomy"]["Kingdom"]["Pseudomonadati"],
                600,
            )

            domain_chart = REPORT.svg_composition(
                {"S-1": {"Bacteria": 60, "Eukaryota": 40}},
                "Domain composition by sample",
            )
            self.assertGreater(
                domain_chart.index('class="chart-legend"'),
                domain_chart.index("</svg>"),
            )


if __name__ == "__main__":
    unittest.main()
