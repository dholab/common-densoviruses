import tempfile
import unittest
from pathlib import Path

from build_site import rewrite_repository_links, split_manuscript


REPO = "https://github.com/example/common-densoviruses"


class ManuscriptParsingTests(unittest.TestCase):
    def setUp(self):
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.fixture_root = Path(self.temporary_directory.name)
        (self.fixture_root / "notes.txt").write_text("notes\n", encoding="utf-8")
        (self.fixture_root / "results").mkdir()

    def tearDown(self):
        self.temporary_directory.cleanup()

    def test_split_manuscript_extracts_masthead(self):
        source = """# Test title

Alice Example<sup>1</sup>, Bob Example<sup>2*</sup>

<sup>1</sup> First Institute<br>
<sup>2*</sup> Corresponding author: <bob@example.org>

## Introduction

Opening paragraph.
"""
        manuscript = split_manuscript(source)
        self.assertEqual(manuscript.title, "Test title")
        self.assertIn("Alice Example", manuscript.authors_html)
        self.assertIn("First Institute", manuscript.affiliations_html)
        self.assertTrue(manuscript.body_markdown.startswith("## Introduction"))

    def test_rewrite_repository_links_distinguishes_files_and_directories(self):
        rendered = (
            '<a href="notes.txt">file</a> <a href="results">directory</a> '
            '<a href="https://example.org">external</a> '
            '<a href="mailto:a@example.org">mail</a> <a href="#figure-1">anchor</a>'
        )
        rewritten = rewrite_repository_links(rendered, self.fixture_root, REPO)
        self.assertIn(f'{REPO}/blob/main/notes.txt', rewritten)
        self.assertIn(f'{REPO}/tree/main/results', rewritten)
        self.assertIn('href="https://example.org"', rewritten)
        self.assertIn('href="mailto:a@example.org"', rewritten)
        self.assertIn('href="#figure-1"', rewritten)
