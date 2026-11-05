"""Exercise routing and cache identity without downloading scientific data."""

import ast
import hashlib
import importlib.util
import itertools
import json
from pathlib import Path
import runpy
import tempfile
import types
import unittest
from unittest import mock
import zipfile


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "database_config", ROOT / "workflow/scripts/database_config.py"
)
DATABASE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(DATABASE)
CAPABILITY_SPEC = importlib.util.spec_from_file_location(
    "check_reference_capability", ROOT / "workflow/scripts/check_reference_capability.py"
)
CAPABILITY = importlib.util.module_from_spec(CAPABILITY_SPEC)
CAPABILITY_SPEC.loader.exec_module(CAPABILITY)


def write_classifier(path, semantic_type="TaxonomicClassifier", provenance_type="FeatureData[Sequence]"):
    """Root artifact type differs from the deliberately misleading provenance."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(path, "w") as archive:
        archive.writestr("artifact/metadata.yaml", f"uuid: artifact\ntype: {semantic_type}\n")
        archive.writestr("artifact/VERSION", "QIIME 2\narchive: 7\nframework: 2026.7.0\n")
        archive.writestr("artifact/data/classifier.pkl", b"synthetic test fixture")
        archive.writestr(
            "artifact/provenance/artifacts/ancestor/metadata.yaml",
            f"type: {provenance_type}\n",
        )
    return path


class DatabaseChoiceTests(unittest.TestCase):
    def test_bbsplit_reuse_uses_only_the_specific_setting(self):
        source = (ROOT / "workflow/rules/common.smk").read_text(encoding="utf-8")
        preamble = ast.parse(source.split("# SILVA 144 classifiers", 1)[0])
        assignment = next(
            node for node in preamble.body
            if isinstance(node, ast.Assign)
            and any(isinstance(target, ast.Name)
                    and target.id == "USE_PREEXISTING_BBSPLIT_DATABASE"
                    for target in node.targets)
        )
        expression = compile(ast.Expression(assignment.value), "common.smk", "eval")
        for config, expected in (
            ({}, False),
            ({"use_preexisting_databases": True}, False),
            ({"use_preexisting_bbsplit_database": True}, True),
            ({"use_preexisting_bbsplit_database": False,
              "use_preexisting_databases": True}, False),
        ):
            with self.subTest(config=config):
                self.assertIs(eval(expression, {"config": config}), expected)
        self.assertNotIn('"use_preexisting_databases"', source)

    def base_config(self, root):
        return dict(
            database_dir=str(Path(root) / "shared databases"),
            fwdPrimer="GTGYCAGCMGCCGCGGTAA",
            revPrimer="CCGYCAATTYMTTTRAGTTT",
        )

    def test_all_nine_independent_database_combinations(self):
        with tempfile.TemporaryDirectory() as tmp:
            for silva_mode, pr2_mode in itertools.product(
                ("existing", "download", "train"), repeat=2
            ):
                with self.subTest(SILVA=silva_mode, PR2=pr2_mode):
                    config = self.base_config(tmp)
                    for family, mode in (("SILVA", silva_mode), ("PR2", pr2_mode)):
                        lower = family.lower()
                        config["use_preexisting_" + lower + "_database"] = mode == "existing"
                        if mode == "existing":
                            config[lower + "_classifier"] = str(write_classifier(
                                Path(tmp) / "provided" / (family + ".qza")
                            ))
                        elif mode == "download":
                            config[family + "_download_link"] = (
                                f"https://references.example/{family}.qza"
                            )
                    self.assertEqual(DATABASE.resolve_database(config, "SILVA")["mode"], silva_mode)
                    self.assertEqual(DATABASE.resolve_database(config, "PR2")["mode"], pr2_mode)

    def test_existing_mode_uses_exact_file_and_ignores_inactive_url(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = write_classifier(Path(tmp) / "read only reference" / "provided.qza")
            before = source.read_bytes()
            config = self.base_config(tmp)
            config.update(
                use_preexisting_silva_database=True,
                silva_classifier=str(source),
                SILVA_download_link="not an active URL",
                SILVA_download_md5="not an active checksum",
            )
            record = DATABASE.resolve_database(config, "SILVA")
            self.assertEqual(record["classifier"], str(source))
            self.assertEqual(record["mode"], "existing")
            self.assertEqual(record["url"], "")
            self.assertIsNone(record["checksum"])
            self.assertTrue(record["lock"].startswith(".snakemake/"))
            self.assertEqual(source.read_bytes(), before)

    def test_existing_missing_or_invalid_file_fails_without_fallback(self):
        with tempfile.TemporaryDirectory() as tmp:
            config = self.base_config(tmp)
            config.update(use_preexisting_pr2_database=True,
                          pr2_classifier=str(Path(tmp) / "missing.qza"),
                          PR2_download_link="https://references.example/backup.qza")
            with self.assertRaisesRegex(ValueError, "not found"):
                DATABASE.resolve_database(config, "PR2")
            invalid = Path(tmp) / "invalid.qza"
            invalid.write_text("an HTML error page", encoding="utf-8")
            config["pr2_classifier"] = str(invalid)
            with self.assertRaisesRegex(ValueError, "TaxonomicClassifier"):
                DATABASE.resolve_database(config, "PR2")

    def test_blank_null_and_absent_urls_select_training(self):
        with tempfile.TemporaryDirectory() as tmp:
            for family in ("SILVA", "PR2"):
                for value in (None, "", " \t "):
                    with self.subTest(family=family, value=value):
                        config = self.base_config(tmp)
                        config[family + "_download_link"] = value
                        self.assertEqual(DATABASE.resolve_database(config, family)["mode"], "train")
                self.assertEqual(DATABASE.resolve_database(self.base_config(tmp), family)["mode"], "train")

    def test_absolute_per_database_destinations_and_default_locations(self):
        with tempfile.TemporaryDirectory() as tmp:
            config = self.base_config(tmp)
            for family in ("SILVA", "PR2"):
                record = DATABASE.resolve_database(config, family)
                expected = Path(config["database_dir"]) / "classification" / family
                self.assertEqual(Path(record["classifier"]).parent, expected)
                config[family.lower() + "_database_dir"] = str(Path(tmp) / "custom" / family)
                record = DATABASE.resolve_database(config, family)
                self.assertEqual(Path(record["classifier"]).parent, Path(tmp) / "custom" / family)
                self.assertTrue(Path(record["classifier"]).is_absolute())
                self.assertFalse(Path(record["build_root"]).is_absolute())

    def test_invalid_active_urls_checksums_versions_and_regions_fail(self):
        with tempfile.TemporaryDirectory() as tmp:
            for url in ("/local/classifier.qza", "file:///tmp/example.qza", "ftp://host/file", "https:///missing-host"):
                config = self.base_config(tmp)
                config["PR2_download_link"] = url
                with self.subTest(url=url), self.assertRaisesRegex(ValueError, "HTTP"):
                    DATABASE.resolve_database(config, "PR2")
            for field, value in (("PR2_download_md5", "x" * 32), ("PR2version", "../../5.1.1"), ("pr2_reference_mode", "unknown")):
                config = self.base_config(tmp)
                config["PR2_download_link"] = "https://references.example/pr2.qza"
                config[field] = value
                with self.subTest(field=field), self.assertRaises(ValueError):
                    DATABASE.resolve_database(config, "PR2")

    def test_custom_download_cache_identity_includes_url_checksum_and_release(self):
        with tempfile.TemporaryDirectory() as tmp:
            config = self.base_config(tmp)
            config["PR2_download_link"] = "https://references.example/pr2.qza"
            base = DATABASE.resolve_database(config, "PR2")
            for field, value in (
                ("PR2_download_link", "https://references.example/other.qza"),
                ("PR2_download_md5", "a" * 32),
                ("PR2version", "5.1.0"),
            ):
                changed = DATABASE.resolve_database(dict(config, **{field: value}), "PR2")
                with self.subTest(field=field):
                    self.assertNotEqual(base["classifier"], changed["classifier"])
                    self.assertNotEqual(base["ready"], changed["ready"])
                    self.assertNotEqual(base["build_root"], changed["build_root"])

    def test_official_silva_download_retains_filename_and_checksum(self):
        with tempfile.TemporaryDirectory() as tmp:
            for url, (name, checksum) in DATABASE.SILVA_PRESETS.items():
                config = self.base_config(tmp)
                config["SILVA_download_link"] = url
                record = DATABASE.resolve_database(config, "SILVA")
                self.assertEqual(Path(record["classifier"]).name, name)
                self.assertEqual(record["checksum"], checksum)
                config["SILVA_download_md5"] = "0" * 32
                with self.assertRaisesRegex(ValueError, "disagrees"):
                    DATABASE.resolve_database(config, "SILVA")

    def test_known_silva_preset_cannot_be_declared_as_an_older_release(self):
        with tempfile.TemporaryDirectory() as tmp:
            config = self.base_config(tmp)
            config.update(SILVAversion="138.2", SILVA_download_link=next(iter(DATABASE.SILVA_PRESETS)))
            with self.assertRaisesRegex(ValueError, "requires SILVAversion"):
                DATABASE.resolve_database(config, "SILVA")

    def test_effective_training_recipe_is_recorded_only_for_training(self):
        with tempfile.TemporaryDirectory() as tmp:
            config = self.base_config(tmp)
            trained = DATABASE.resolve_database(config, "SILVA")
            self.assertIn("Bacteria 1200", trained["training_settings"]["Minimum lengths (nt)"])
            self.assertEqual(trained["training_settings"]["Forward primer"], config["fwdPrimer"])
            config["SILVA_download_link"] = next(iter(DATABASE.SILVA_PRESETS))
            downloaded = DATABASE.resolve_database(config, "SILVA")
            self.assertEqual(downloaded["training_settings"], {})

    def test_training_cache_identity_separates_recipe_release_environment_and_region(self):
        with tempfile.TemporaryDirectory() as tmp:
            for family in ("SILVA", "PR2"):
                config = self.base_config(tmp)
                base = DATABASE.resolve_database(config, family)
                self.assertIn("rescript-v1", base["classifier"])
                cases = (("fwdPrimer", "AAAA"),
                         (family.lower() + "_reference_mode", "full_length"),
                         ("_database_environment_id", "different-env-definition"),
                         ("SILVAversion" if family == "SILVA" else "PR2version", "138.2" if family == "SILVA" else "5.1.0"))
                for key, value in cases:
                    with self.subTest(family=family, key=key):
                        changed = DATABASE.resolve_database(dict(config, **{key: value}), family)
                        self.assertNotEqual(base["classifier"], changed["classifier"])
                        self.assertNotEqual(base["build_root"], changed["build_root"])
                full = dict(config, **{family.lower() + "_reference_mode": "full_length"})
                self.assertEqual(
                    DATABASE.resolve_database(full, family)["classifier"],
                    DATABASE.resolve_database(dict(full, fwdPrimer="AAAA"), family)["classifier"],
                )

    def test_download_and_training_cannot_share_cache_identity(self):
        with tempfile.TemporaryDirectory() as tmp:
            for family in ("SILVA", "PR2"):
                config = self.base_config(tmp)
                trained = DATABASE.resolve_database(config, family)
                downloaded = DATABASE.resolve_database(dict(config, **{family + "_download_link": "https://references.example/model.qza"}), family)
                self.assertNotEqual(trained["classifier"], downloaded["classifier"])

    def test_existing_mode_recognizes_legacy_locations(self):
        with tempfile.TemporaryDirectory() as tmp:
            config = self.base_config(tmp)
            primers = config["fwdPrimer"] + "_" + config["revPrimer"]
            names = {
                "SILVA": "silva-ssu-nr99-tax-dereplicated-sliced_" + primers + "_dereplicated_final_classifier_USE_ME.qza",
                "PR2": "pr2_version_5.1.1_SSU_dada2.clean.culled.derep-sliced_" + primers + "_dereplicated_final_classifier_qiime2-2026.7.qza",
            }
            for family, name in names.items():
                provided = write_classifier(Path(config["database_dir"]) / "classification" / family / name)
                config["use_preexisting_" + family.lower() + "_database"] = True
                self.assertEqual(DATABASE.resolve_database(config, family)["classifier"], str(provided))

    def test_only_root_metadata_can_identify_a_classifier(self):
        with tempfile.TemporaryDirectory() as tmp:
            correct = write_classifier(Path(tmp) / "correct.qza")
            self.assertTrue(DATABASE.reusable_classifier(correct))
            wrong = write_classifier(Path(tmp) / "wrong.qza", "FeatureData[Sequence]", "TaxonomicClassifier")
            self.assertFalse(DATABASE.reusable_classifier(wrong))
            checksum = hashlib.md5(correct.read_bytes()).hexdigest()
            self.assertTrue(DATABASE.reusable_classifier(correct, checksum))
            self.assertFalse(DATABASE.reusable_classifier(correct, "0" * 32))
            with zipfile.ZipFile(correct, "a") as archive:
                archive.writestr("other/metadata.yaml", "type: TaxonomicClassifier\n")
            self.assertFalse(DATABASE.reusable_classifier(correct))

    def test_capability_check_reports_supported_release_and_fails_unsupported_release(self):
        action = types.SimpleNamespace(signature=types.SimpleNamespace(parameters={
            "version": types.SimpleNamespace(qiime_type={"138.2", "138.1"})
        }))
        manager = types.SimpleNamespace(plugins={"rescript": types.SimpleNamespace(
            actions={"get_silva_data": action}
        )})
        modules = {
            "qiime2": types.SimpleNamespace(__version__="2026.7.0"),
            "qiime2.sdk": types.SimpleNamespace(PluginManager=lambda: manager),
            "rescript": types.SimpleNamespace(__version__="2026.7.0"),
            "sklearn": types.SimpleNamespace(__version__="1.7.1"),
        }
        with tempfile.TemporaryDirectory() as tmp, mock.patch.dict("sys.modules", modules):
            output = Path(tmp) / "capability.json"
            with self.assertRaisesRegex(ValueError, "cannot download SILVA 144"):
                CAPABILITY.check_capability("SILVA", "144", output)
            self.assertFalse(output.exists())
            CAPABILITY.check_capability("SILVA", "138.2", output)
            result = json.loads(output.read_text())
            self.assertEqual(result["version"], "138.2")
            self.assertEqual(result["rescript"], "2026.7.0")

    def test_capability_script_reads_database_parameter_for_both_references(self):
        releases = {"SILVA": "138.2", "PR2": "5.1.1"}
        actions = {
            action: types.SimpleNamespace(signature=types.SimpleNamespace(parameters={
                "version": types.SimpleNamespace(qiime_type={releases[database]})
            }))
            for database, action in (("SILVA", "get_silva_data"), ("PR2", "get_pr2_data"))
        }
        manager = types.SimpleNamespace(plugins={"rescript": types.SimpleNamespace(actions=actions)})
        modules = {
            "qiime2": types.SimpleNamespace(__version__="2026.7.0"),
            "qiime2.sdk": types.SimpleNamespace(PluginManager=lambda: manager),
            "rescript": types.SimpleNamespace(__version__="2026.7.0"),
            "sklearn": types.SimpleNamespace(__version__="1.7.1"),
        }
        rules = (ROOT / "workflow/rules/prepare_classification_dbs.smk").read_text()
        with tempfile.TemporaryDirectory() as tmp, mock.patch.dict("sys.modules", modules):
            for database, version in releases.items():
                with self.subTest(database=database):
                    block = rules.split(f"rule check_{database.lower()}_reference_capability:", 1)[1].split("\nif ", 1)[0]
                    self.assertIn(f'database="{database}"', block)
                    output = Path(tmp) / database / "capability.json"
                    snakemake = types.SimpleNamespace(
                        params=types.SimpleNamespace(database=database, version=version),
                        output=[str(output)],
                    )
                    runpy.run_path(str(ROOT / "workflow/scripts/check_reference_capability.py"),
                                   init_globals={"snakemake": snakemake})
                    result = json.loads(output.read_text())
                    self.assertEqual(result["database"], database)
                    self.assertEqual(result["version"], version)


if __name__ == "__main__":
    unittest.main()
