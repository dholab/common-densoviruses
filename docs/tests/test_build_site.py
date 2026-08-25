import tempfile
import unittest
from pathlib import Path
from xml.etree import ElementTree

from build_site import (
    BuildError,
    collect_figures,
    make_svg_transparent,
    rewrite_repository_links,
    split_manuscript,
    stage_figure_assets,
)


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


class TransparentSvgTests(unittest.TestCase):
    def setUp(self):
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.fixture_root = Path(self.temporary_directory.name)
        self.source = self.fixture_root / "source.svg"
        self.destination = self.fixture_root / "output.svg"

    def tearDown(self):
        self.temporary_directory.cleanup()

    def write_svg(self, contents: str) -> Path:
        self.source.write_text(contents, encoding="utf-8")
        return self.source

    def test_removes_only_direct_full_canvas_white_rect(self):
        """Removing the direct canvas background must retain scientific marks."""
        source = self.write_svg(
            '<svg xmlns="http://www.w3.org/2000/svg" width="100" height="60">'
            '<rect x="0" y="0" width="100" height="60" fill="white"/>'
            '<rect x="10" y="10" width="20" height="20" fill="#163139"/>'
            "</svg>"
        )

        make_svg_transparent(source, self.destination)

        root = ElementTree.parse(self.destination).getroot()
        self.assertEqual(len(list(root)), 1)
        self.assertEqual(list(root)[0].attrib["fill"], "#163139")

    def test_removes_background_group_with_inherited_white_fill(self):
        """A white background-only group must be removed as one unit."""
        source = self.write_svg(
            '<svg xmlns="http://www.w3.org/2000/svg" width="100px" height="60px">'
            '<g style="fill: #fff"><rect x="0" y="0" width="100" height="60"/></g>'
            '<rect x="10" y="10" width="20" height="20" fill="#163139"/>'
            "</svg>"
        )

        make_svg_transparent(source, self.destination)

        root = ElementTree.parse(self.destination).getroot()
        self.assertEqual(len(list(root)), 1)
        self.assertEqual(list(root)[0].attrib["fill"], "#163139")

    def test_removes_igv_percent_canvas_backdrop_with_default_origin(self):
        """The IGV root backdrop is the sole safe percentage-dimension exception."""
        source = self.write_svg(
            '<svg xmlns="http://www.w3.org/2000/svg" width="100" height="60">'
            '<rect id="svg_output_backdrop" width="100%" height="100%" fill="white"/>'
            '<rect x="10" y="10" width="20" height="20" fill="#163139"/>'
            "</svg>"
        )

        make_svg_transparent(source, self.destination)

        root = ElementTree.parse(self.destination).getroot()
        self.assertEqual(len(list(root)), 1)
        self.assertEqual(list(root)[0].attrib["fill"], "#163139")

    def test_rejects_svg_without_expected_canvas_background(self):
        """A small white annotation must never be mistaken for a background."""
        source = self.write_svg(
            '<svg xmlns="http://www.w3.org/2000/svg" width="100" height="60">'
            '<rect x="10" y="10" width="20" height="20" fill="white"/>'
            "</svg>"
        )

        with self.assertRaises(BuildError):
            make_svg_transparent(source, self.destination)

        self.assertFalse(self.destination.exists())


class FigureCollectionTests(unittest.TestCase):
    def setUp(self):
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.fixture_root = Path(self.temporary_directory.name)
        (self.fixture_root / "figures").mkdir()
        (self.fixture_root / "figures" / "one.svg").write_text(
            '<svg xmlns="http://www.w3.org/2000/svg" width="1" height="1">'
            '<rect x="0" y="0" width="1" height="1" fill="white"/></svg>',
            encoding="utf-8",
        )

    def tearDown(self):
        self.temporary_directory.cleanup()

    def test_collects_allowed_labels_and_stages_transparent_assets(self):
        """The manuscript figure syntax yields stable published asset paths."""
        markdown_text = "![Figure 5A](figures/one.svg)"

        figures = collect_figures(markdown_text, self.fixture_root)
        staged = stage_figure_assets(figures, self.fixture_root / "site")

        self.assertEqual(figures[0].label, "Figure 5A")
        self.assertEqual(figures[0].asset_name, "figure-5a.svg")
        self.assertEqual(staged, {"figures/one.svg": "assets/figures/figure-5a.svg"})
        published = self.fixture_root / "site" / "assets" / "figures" / "figure-5a.svg"
        self.assertTrue(published.is_file())
        self.assertEqual(len(list(ElementTree.parse(published).getroot())), 0)

    def test_rejects_duplicate_figure_sources(self):
        """One scientific SVG cannot silently represent two manuscript figures."""
        markdown_text = "![Figure 1](figures/one.svg)\n![Figure 2](figures/one.svg)"

        with self.assertRaises(BuildError):
            collect_figures(markdown_text, self.fixture_root)
