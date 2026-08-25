import tempfile
import unittest
from pathlib import Path
import struct
from xml.etree import ElementTree
import zlib

from build_site import (
    BuildError,
    build,
    collect_figures,
    make_svg_transparent,
    render_markdown,
    render_page,
    rewrite_repository_links,
    split_manuscript,
    stage_figure_assets,
    wrap_figures,
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

    def test_rejects_non_svg_xml_without_publishing_destination(self):
        """A matching rectangle in arbitrary XML must not become a figure asset."""
        source = self.write_svg(
            '<document width="100" height="60">'
            '<rect x="0" y="0" width="100" height="60" fill="white"/>'
            "</document>"
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


class RepositoryContractTests(unittest.TestCase):
    """Protect the committed GitHub Pages manuscript surface."""

    def test_real_manuscript_publishes_all_figures_and_preserves_viewers(self):
        """Missing or stale committed output must fail before a rebuild can hide it."""
        repository_root = Path(__file__).resolve().parents[2]
        docs = repository_root / "docs"
        expected_labels = [
            "Figure 1", "Figure 2", "Figure 3", "Figure 4", "Figure 5A",
            "Figure 5B", "Figure 6", "Figure 7A", "Figure 7B",
        ]
        expected_assets = [
            "figure-1.svg", "figure-2.svg", "figure-3.svg", "figure-4.svg",
            "figure-5a.svg", "figure-5b.svg", "figure-6.svg", "figure-7a.svg",
            "figure-7b.svg",
        ]
        viewer_routes = [
            "01-air-samples-tgs",
            "01-air-samples-vsp",
            "02-sra-mining",
            "04-skin-microbiome",
            "05-cameroonian-plasma",
            "06-belgian-air",
        ]

        figures = collect_figures(
            (repository_root / "README.md").read_text(encoding="utf-8"), repository_root
        )
        self.assertEqual([figure.label for figure in figures], expected_labels)
        self.assertTrue(all(figure.source.is_file() for figure in figures))
        self.assertTrue(all((docs / route / "index.html").is_file() for route in viewer_routes))

        # Check the committed artifact first: this deliberately precedes build().
        published_assets = docs / "assets" / "figures"
        self.assertEqual(
            sorted(asset.name for asset in published_assets.glob("*.svg")), expected_assets
        )
        committed_page = (docs / "index.html").read_text(encoding="utf-8")
        self.assertEqual(committed_page.count('<figure class="display" id="figure-'), 9)

        output = build(repository_root, docs)
        rebuilt_page = output.read_text(encoding="utf-8")
        self.assertEqual(output, docs / "index.html")
        self.assertEqual(rebuilt_page.count('<figure class="display" id="figure-'), 9)


class PageRenderingTests(unittest.TestCase):
    def setUp(self):
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.fixture_root = Path(self.temporary_directory.name)
        self.docs = self.fixture_root / "docs"
        self.docs.mkdir()

    def tearDown(self):
        self.temporary_directory.cleanup()

    def write_complete_fixture(self, template: str) -> None:
        labels = [
            "Figure 1", "Figure 2", "Figure 3", "Figure 4", "Figure 5A",
            "Figure 5B", "Figure 6", "Figure 7A", "Figure 7B",
        ]
        manuscript = [
            "# Fixture manuscript",
            "",
            "Alice Example<sup>1*</sup>",
            "",
            "<sup>1</sup> Example Institute",
            "",
            "## Introduction",
            "",
            "Opening evidence paragraph.",
        ]
        for number, label in enumerate(labels, start=1):
            source = f"figures/{number}.svg"
            svg = self.fixture_root / source
            svg.parent.mkdir(parents=True, exist_ok=True)
            svg.write_text(
                '<svg xmlns="http://www.w3.org/2000/svg" width="1" height="1">'
                '<rect x="0" y="0" width="1" height="1" fill="white"/></svg>',
                encoding="utf-8",
            )
            manuscript.extend([
                "", f"![{label}]({source})", "", f"> **{label}.** Caption {number}.",
            ])
        (self.fixture_root / "README.md").write_text("\n".join(manuscript), encoding="utf-8")
        (self.docs / "template.html").write_text(template, encoding="utf-8")
        (self.docs / "assets").mkdir()
        (self.docs / "assets" / "og.png").write_bytes(b"placeholder social preview")

    def test_wraps_image_and_caption_as_labeled_figure(self):
        """A figure image and its next caption must remain a single display item."""
        body = (
            '<p><img alt="Figure 5A" src="docs/source.svg" /></p>'
            '<blockquote><p><strong>Figure 5A.</strong> Evidence caption.</p></blockquote>'
        )

        wrapped = wrap_figures(body, {"docs/source.svg": "assets/figures/figure-5a.svg"})

        self.assertIn('<figure class="display" id="figure-5a">', wrapped)
        self.assertIn('src="assets/figures/figure-5a.svg"', wrapped)
        self.assertIn('<span class="fig-label">Figure 5A</span>', wrapped)
        self.assertNotIn("<blockquote>", wrapped)

    def test_rejects_image_without_a_matching_caption(self):
        """A missing caption must fail rather than publish an unlabeled scientific figure."""
        body = '<p><img alt="Figure 5A" src="docs/source.svg" /></p>'

        with self.assertRaises(BuildError):
            wrap_figures(body, {"docs/source.svg": "assets/figures/figure-5a.svg"})

    def test_rejects_caption_without_a_matching_image(self):
        """An orphan caption must fail rather than silently detach its evidence."""
        body = '<blockquote><p><strong>Figure 5A.</strong> Evidence caption.</p></blockquote>'

        with self.assertRaises(BuildError):
            wrap_figures(body, {})

    def test_builds_complete_house_style_manuscript_page(self):
        """A complete fixture must render the masthead and all nine linked figures."""
        self.write_complete_fixture(
            "<!doctype html><html><head><title>{{TITLE_TEXT}}</title>"
            '<meta name="description" content="{{DESCRIPTION}}">'
            '<meta name="theme-color" content="#F8F4E9">'
            '<link rel="canonical" href="{{CANONICAL_URL}}"></head><body>'
            '<header class="masthead"><h1>{{TITLE_HTML}}</h1>{{AUTHORS}}{{AFFILIATIONS}}'
            '<p class="standfirst">{{STANDFIRST}}</p></header><main>{{BODY}}</main>'
            '<footer class="colophon">{{REPOSITORY_URL}}</footer></body></html>',
        )

        output = build(self.fixture_root, self.docs)
        page = output.read_text(encoding="utf-8")

        self.assertEqual(output, self.docs / "index.html")
        self.assertIn("Fixture manuscript", page)
        self.assertIn('class="masthead"', page)
        self.assertIn("<h2>Introduction</h2>", page)
        self.assertIn("This web version is the preferred way to read this evolving manuscript", page)
        self.assertIn('class="colophon"', page)
        self.assertIn('name="theme-color" content="#F8F4E9"', page)
        self.assertIn(
            'name="description" content="An evolving manuscript on densoviruses in the human &amp; mammalian virospheres."',
            page,
        )
        self.assertIn('rel="canonical" href="https://dholab.github.io/common-densoviruses/"', page)
        self.assertEqual(page.count('<figure class="display" id="figure-'), 9)
        self.assertNotRegex(page, r"{{[^}]+}}")
        self.assertEqual((self.docs / "assets" / "og.png").read_bytes(), b"placeholder social preview")

    def test_social_metadata_is_complete(self):
        """Shared manuscript links must expose the canonical card metadata and image."""
        template = (Path(__file__).resolve().parents[1] / "template.html").read_text(
            encoding="utf-8"
        )
        self.write_complete_fixture(template)

        page = build(self.fixture_root, self.docs).read_text(encoding="utf-8")
        image_url = "https://dholab.github.io/common-densoviruses/assets/og.png"

        self.assertTrue((self.docs / "assets" / "og.png").is_file())
        self.assertIn('rel="canonical" href="https://dholab.github.io/common-densoviruses/"', page)
        self.assertIn('property="og:title" content="Fixture manuscript"', page)
        self.assertIn(
            'property="og:description" content="An evolving manuscript on densoviruses in the human &amp; mammalian virospheres."',
            page,
        )
        self.assertIn('property="og:type" content="website"', page)
        self.assertIn('property="og:url" content="https://dholab.github.io/common-densoviruses/"', page)
        self.assertIn(f'property="og:image" content="{image_url}"', page)
        self.assertIn('name="twitter:card" content="summary_large_image"', page)
        self.assertIn('name="twitter:title" content="Fixture manuscript"', page)
        self.assertIn(
            'name="twitter:description" content="An evolving manuscript on densoviruses in the human &amp; mammalian virospheres."',
            page,
        )
        self.assertIn(f'name="twitter:image" content="{image_url}"', page)

    def test_social_preview_uses_only_approved_house_palette(self):
        """The social card must retain only the three exact house-palette colors."""
        preview = Path(__file__).resolve().parents[1] / "assets" / "og.png"
        data = preview.read_bytes()
        self.assertEqual(data[:8], b"\x89PNG\r\n\x1a\n")

        offset = 8
        header = None
        compressed = bytearray()
        while offset < len(data):
            length = struct.unpack(">I", data[offset : offset + 4])[0]
            chunk_type = data[offset + 4 : offset + 8]
            chunk = data[offset + 8 : offset + 8 + length]
            offset += 12 + length
            if chunk_type == b"IHDR":
                header = struct.unpack(">IIBBBBB", chunk)
            elif chunk_type == b"IDAT":
                compressed.extend(chunk)
            elif chunk_type == b"IEND":
                break

        self.assertIsNotNone(header)
        width, height, bit_depth, color_type, compression, filtering, interlace = header
        self.assertEqual((width, height), (1200, 630))
        self.assertEqual((bit_depth, color_type, compression, filtering, interlace), (8, 2, 0, 0, 0))

        row_bytes = width * 3
        inflated = zlib.decompress(compressed)
        self.assertEqual(len(inflated), height * (row_bytes + 1))
        previous = bytearray(row_bytes)
        colors = set()
        cursor = 0

        def paeth(left: int, above: int, upper_left: int) -> int:
            estimate = left + above - upper_left
            left_distance = abs(estimate - left)
            above_distance = abs(estimate - above)
            upper_left_distance = abs(estimate - upper_left)
            if left_distance <= above_distance and left_distance <= upper_left_distance:
                return left
            if above_distance <= upper_left_distance:
                return above
            return upper_left

        for _ in range(height):
            filter_type = inflated[cursor]
            cursor += 1
            filtered = inflated[cursor : cursor + row_bytes]
            cursor += row_bytes
            row = bytearray(row_bytes)
            for index, value in enumerate(filtered):
                left = row[index - 3] if index >= 3 else 0
                above = previous[index]
                upper_left = previous[index - 3] if index >= 3 else 0
                if filter_type == 0:
                    row[index] = value
                elif filter_type == 1:
                    row[index] = (value + left) & 0xFF
                elif filter_type == 2:
                    row[index] = (value + above) & 0xFF
                elif filter_type == 3:
                    row[index] = (value + ((left + above) // 2)) & 0xFF
                elif filter_type == 4:
                    row[index] = (value + paeth(left, above, upper_left)) & 0xFF
                else:
                    self.fail(f"Unsupported PNG filter type: {filter_type}")
            colors.update(tuple(row[index : index + 3]) for index in range(0, row_bytes, 3))
            previous = row

        self.assertSetEqual(colors, {(248, 244, 233), (22, 49, 57), (193, 106, 60)})

    def test_rejects_completed_page_with_missing_local_asset(self):
        """A local page asset must exist before a build replaces the published page."""
        self.write_complete_fixture(
            "<!doctype html><html><head><title>{{TITLE_TEXT}}</title></head><body>"
            '<header class="masthead"><h1>{{TITLE_HTML}}</h1>{{AUTHORS}}{{AFFILIATIONS}}'
            '<p class="standfirst">{{STANDFIRST}}</p></header><main>{{BODY}}'
            '<img src="assets/missing.svg"></main>'
            '<footer class="colophon">{{REPOSITORY_URL}}</footer></body></html>'
        )

        with self.assertRaisesRegex(BuildError, "local asset"):
            build(self.fixture_root, self.docs)

    def test_rejects_repository_only_single_quoted_asset(self):
        """A source-only repository file must not be accepted as a Pages asset."""
        self.write_complete_fixture(
            "<!doctype html><html><head><title>{{TITLE_TEXT}}</title></head><body>"
            '<header class="masthead"><h1>{{TITLE_HTML}}</h1>{{AUTHORS}}{{AFFILIATIONS}}'
            '<p class="standfirst">{{STANDFIRST}}</p></header><main>{{BODY}}'
            "<img src='README.md'></main>"
            '<footer class="colophon">{{REPOSITORY_URL}}</footer></body></html>'
        )

        with self.assertRaisesRegex(BuildError, "local asset"):
            build(self.fixture_root, self.docs)

    def test_rejects_stale_single_quoted_generated_figure(self):
        """A figure that publication removes must not validate from the old output."""
        self.write_complete_fixture(
            "<!doctype html><html><head><title>{{TITLE_TEXT}}</title></head><body>"
            '<header class="masthead"><h1>{{TITLE_HTML}}</h1>{{AUTHORS}}{{AFFILIATIONS}}'
            '<p class="standfirst">{{STANDFIRST}}</p></header><main>{{BODY}}'
            "<img src='assets/figures/stale.svg'></main>"
            '<footer class="colophon">{{REPOSITORY_URL}}</footer></body></html>'
        )
        stale = self.docs / "assets" / "figures" / "stale.svg"
        stale.parent.mkdir()
        stale.write_text("<svg/>", encoding="utf-8")

        with self.assertRaisesRegex(BuildError, "local asset"):
            build(self.fixture_root, self.docs)
