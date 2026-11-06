"""Keep BLAST environments minimal without upgrading the search engine."""
from pathlib import Path
import unittest

import yaml


ROOT = Path(__file__).resolve().parents[1]


class BlastEnvironmentTests(unittest.TestCase):
    def read_environment(self, filename):
        return yaml.safe_load((ROOT / 'workflow/envs' / filename).read_text())

    def test_channels_prioritize_conda_forge_and_exclude_inherited_defaults(self):
        for filename in ('proportal.yml', 'blast-env.yaml'):
            with self.subTest(filename=filename):
                environment = self.read_environment(filename)
                self.assertEqual(environment['channels'], ['conda-forge', 'bioconda', 'nodefaults'])
                self.assertEqual(environment['name'], 'blast-env')

    def test_proportal_declares_only_the_tools_its_rule_needs(self):
        self.assertEqual(self.read_environment('proportal.yml')['dependencies'],
                         ['blast=2.12.0', 'seqtk', 'python=3.10', 'git'])

    def test_internal_standard_blast_keeps_the_existing_search_version(self):
        self.assertEqual(self.read_environment('blast-env.yaml')['dependencies'],
                         ['blast=2.12.0'])

    def test_rules_still_select_the_same_environment_files(self):
        prok = (ROOT / 'workflow/rules/02-denoise-and-export-prok.smk').read_text()
        standards = (ROOT / 'workflow/rules/05-internal-standard-correction.smk').read_text()
        self.assertIn('"../envs/proportal.yml"', prok)
        self.assertIn('"../envs/blast-env.yaml"', standards)


if __name__ == '__main__':
    unittest.main()
