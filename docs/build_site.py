"""Helpers for building the manuscript GitHub Pages site from Markdown."""

from dataclasses import dataclass
from html import escape, unescape
import os
from pathlib import Path
import re
import shutil
import tempfile
from urllib.parse import urlsplit
from xml.etree import ElementTree

import markdown


DOCS = Path(__file__).resolve().parent
ROOT = DOCS.parent
REPOSITORY_URL = "https://github.com/dholab/common-densoviruses"
CANONICAL_URL = "https://dholab.github.io/common-densoviruses/"
STANDFIRST = (
    "This web version is the preferred way to read this evolving manuscript "
    "because each figure links directly to its interactive genomic viewer, "
    "methods, scripts, and source data in the repository."
)


class BuildError(ValueError):
    """Raised when a manuscript asset cannot be built without ambiguity."""


@dataclass(frozen=True)
class Figure:
    label: str
    source: Path
    asset_name: str


_FIGURE_LABELS = {
    "Figure 1",
    "Figure 2",
    "Figure 3",
    "Figure 4",
    "Figure 5A",
    "Figure 5B",
    "Figure 6",
    "Figure 7A",
    "Figure 7B",
}
_SVG_NUMBER = re.compile(r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:px)?\Z", re.IGNORECASE)
_WHITE_FILLS = {"white", "#fff", "#ffffff", "rgb(255,255,255)"}


def _local_name(element: ElementTree.Element) -> str:
    return element.tag.rsplit("}", 1)[-1]


def _numeric_svg_value(value: str | None) -> float | None:
    if value is None or _SVG_NUMBER.fullmatch(value.strip()) is None:
        return None
    return float(value.strip()[:-2] if value.strip().lower().endswith("px") else value.strip())


def _style_value(style: str | None, name: str) -> str | None:
    if style is None:
        return None
    for declaration in style.split(";"):
        property_name, separator, value = declaration.partition(":")
        if separator and property_name.strip().lower() == name:
            return value.strip()
    return None


def _is_white(fill: str | None) -> bool:
    return fill is not None and re.sub(r"\s+", "", fill.lower()) in _WHITE_FILLS


def _effective_fills(root: ElementTree.Element) -> dict[ElementTree.Element, str | None]:
    fills: dict[ElementTree.Element, str | None] = {}

    def visit(element: ElementTree.Element, inherited_fill: str | None) -> None:
        fill = _style_value(element.attrib.get("style"), "fill")
        if fill is None:
            fill = element.attrib.get("fill", inherited_fill)
        fills[element] = fill
        for child in element:
            visit(child, fill)

    visit(root, None)
    return fills


def _canvas_size(root: ElementTree.Element) -> tuple[float, float]:
    width = _numeric_svg_value(root.attrib.get("width"))
    height = _numeric_svg_value(root.attrib.get("height"))
    if width is None or height is None or width <= 0 or height <= 0:
        raise BuildError("SVG must declare positive numeric width and height")
    return width, height


def make_svg_transparent(source: Path, destination: Path) -> None:
    """Publish a validated SVG copy with its single white canvas removed."""
    try:
        tree = ElementTree.parse(source)
    except (OSError, ElementTree.ParseError) as error:
        raise BuildError(f"Cannot parse SVG source {source}") from error

    root = tree.getroot()
    if root.tag != "{http://www.w3.org/2000/svg}svg":
        raise BuildError(f"SVG source must have an SVG root element: {source}")
    canvas_width, canvas_height = _canvas_size(root)
    parents = {child: parent for parent in root.iter() for child in parent}
    fills = _effective_fills(root)
    candidates: list[ElementTree.Element] = []
    for element in root.iter():
        if _local_name(element) != "rect":
            continue
        x = _numeric_svg_value(element.attrib.get("x"))
        y = _numeric_svg_value(element.attrib.get("y"))
        width = _numeric_svg_value(element.attrib.get("width"))
        height = _numeric_svg_value(element.attrib.get("height"))
        numeric_canvas = (
            x == 0
            and y == 0
            and width == canvas_width
            and height == canvas_height
        )
        igv_root_backdrop = (
            parents.get(element) is root
            and element.attrib.get("id") == "svg_output_backdrop"
            and "x" not in element.attrib
            and "y" not in element.attrib
            and element.attrib.get("width", "").strip() == "100%"
            and element.attrib.get("height", "").strip() == "100%"
        )
        if (numeric_canvas or igv_root_backdrop) and _is_white(fills[element]):
            candidates.append(element)

    if len(candidates) != 1:
        raise BuildError(f"Expected exactly one white canvas background in {source}")

    background = candidates[0]
    parent = parents.get(background)
    if parent is None:
        raise BuildError(f"Cannot remove SVG root from {source}")
    if _local_name(parent) == "g" and len(parent) == 1:
        grandparent = parents.get(parent)
        if grandparent is None:
            raise BuildError(f"Cannot remove SVG background group from {source}")
        grandparent.remove(parent)
    else:
        parent.remove(background)

    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="wb", suffix=".svg", prefix=f".{destination.name}.",
            dir=destination.parent, delete=False
        ) as temporary_file:
            temporary_path = Path(temporary_file.name)
            tree.write(temporary_file, encoding="utf-8", xml_declaration=True)
        ElementTree.parse(temporary_path)
        temporary_path.replace(destination)
    except (OSError, ElementTree.ParseError) as error:
        raise BuildError(f"Cannot publish transparent SVG {destination}") from error
    finally:
        if temporary_path is not None and temporary_path.exists():
            temporary_path.unlink()


def collect_figures(markdown_text: str, root: Path) -> list[Figure]:
    """Collect the manuscript's fixed, local SVG figure references."""
    figures: list[Figure] = []
    labels: set[str] = set()
    sources: set[Path] = set()
    image_references = re.finditer(r"!\[([^]]*)\]\(([^\s)]+)(?:\s+[^)]*)?\)", markdown_text)
    for match in image_references:
        label, source_text = match.groups()
        if not label.startswith("Figure "):
            continue
        if label not in _FIGURE_LABELS:
            raise BuildError(f"Unexpected manuscript figure label: {label}")
        source = (root / source_text).resolve()
        try:
            source.relative_to(root.resolve())
        except ValueError as error:
            raise BuildError(f"Figure source is outside the repository: {source_text}") from error
        if source.suffix.lower() != ".svg" or not source.is_file():
            raise BuildError(f"Figure source is not an available SVG: {source_text}")
        if label in labels or source in sources:
            raise BuildError(f"Duplicate manuscript figure label or source: {label}")
        labels.add(label)
        sources.add(source)
        figures.append(Figure(label, source, f"{label.lower().replace(' ', '-')}.svg"))
    return figures


def stage_figure_assets(
    figures: list[Figure], output_dir: Path, root: Path | None = None
) -> dict[str, str]:
    """Stage validated transparent figure copies and return their published paths."""
    output_dir = output_dir.resolve()
    repository_root = (root or output_dir.parent).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    target_directory = output_dir / "assets" / "figures"
    staging_directory = Path(tempfile.mkdtemp(prefix=".figure-assets-", dir=output_dir))
    assets: dict[str, str] = {}
    try:
        for figure in figures:
            try:
                source_path = figure.source.relative_to(repository_root)
            except ValueError as error:
                raise BuildError(f"Figure source is outside the repository: {figure.source}") from error
            asset_path = Path("assets") / "figures" / figure.asset_name
            make_svg_transparent(figure.source, staging_directory / figure.asset_name)
            assets[source_path.as_posix()] = asset_path.as_posix()
        target_directory.mkdir(parents=True, exist_ok=True)
        for staged_asset in staging_directory.iterdir():
            staged_asset.replace(target_directory / staged_asset.name)
        return assets
    finally:
        shutil.rmtree(staging_directory, ignore_errors=True)


@dataclass(frozen=True)
class Manuscript:
    title: str
    authors_html: str
    affiliations_html: str
    body_markdown: str


def split_manuscript(text: str) -> Manuscript:
    """Split a manuscript Markdown document into masthead and body fields."""
    heading = re.match(r"^#\s+(.+?)\s*$", text, re.MULTILINE)
    if heading is None:
        raise ValueError("Manuscript must begin with an H1 title")

    first_section = re.search(r"^##\s+", text, re.MULTILINE)
    if first_section is None:
        raise ValueError("Manuscript must include a level-two section")

    masthead = text[heading.end() : first_section.start()].strip()
    blocks = [block.strip() for block in re.split(r"\n[ \t]*\n", masthead) if block.strip()]
    if len(blocks) != 2:
        raise ValueError("Manuscript must include nonempty author and affiliation blocks")

    return Manuscript(
        title=heading.group(1),
        authors_html=markdown.markdown(blocks[0], extensions=["nl2br"]),
        affiliations_html=markdown.markdown(blocks[1], extensions=["nl2br"]),
        body_markdown=text[first_section.start() :],
    )


def render_markdown(text: str) -> str:
    """Render manuscript body Markdown using the site extensions."""
    return markdown.markdown(
        text,
        extensions=["tables", "fenced_code", "attr_list", "sane_lists", "nl2br"],
    )


def rewrite_repository_links(
    html_text: str, root: Path, repo_url: str, branch: str = "main"
) -> str:
    """Rewrite relative links to their file or directory location in the repository."""
    def replace_link(match: re.Match[str]) -> str:
        quote, href = match.groups()
        decoded_href = unescape(href)
        parsed = urlsplit(decoded_href)
        if parsed.scheme or decoded_href.startswith(("#", "//")):
            return match.group(0)

        repository_path = decoded_href.lstrip("./")
        kind = "tree" if (root / repository_path).is_dir() else "blob"
        return f'href={quote}{repo_url}/{kind}/{branch}/{repository_path}{quote}'

    return re.sub(r"href=(['\"])(.*?)\1", replace_link, html_text)


_FIGURE_IMAGE = re.compile(
    r'<p>(<img alt="(?P<label>Figure [^"]+)" src="(?P<source>[^"]+)"\s*/?>)</p>'
)
_FIGURE_CAPTION = re.compile(
    r'<blockquote>\s*<p><strong>(?P<label>Figure [^.]+)\.\s*</strong>(?P<caption>.*?)</p>\s*</blockquote>',
    re.DOTALL,
)


def wrap_figures(body_html: str, figure_assets: dict[str, str]) -> str:
    """Turn each rendered manuscript image and its following caption into a figure."""
    images = list(_FIGURE_IMAGE.finditer(body_html))
    captions = list(_FIGURE_CAPTION.finditer(body_html))
    if not images or len(images) != len(captions):
        raise BuildError("Every manuscript figure image must have one following caption")
    if len(images) != len(figure_assets):
        raise BuildError("Every staged manuscript figure must appear in the rendered body")

    expected_labels = [image.group("label") for image in images]
    caption_labels = [caption.group("label") for caption in captions]
    if expected_labels != caption_labels:
        raise BuildError("Figure captions must immediately follow their matching images")
    if set(image.group("source") for image in images) != set(figure_assets):
        raise BuildError("Rendered manuscript images do not match staged figure assets")

    figure_pattern = re.compile(
        r'<p>(<img alt="(?P<label>Figure [^"]+)" src="(?P<source>[^"]+)"\s*/?>)</p>'
        r'\s*<blockquote>\s*<p><strong>(?P=label)\.\s*</strong>(?P<caption>.*?)</p>\s*</blockquote>',
        re.DOTALL,
    )

    def replace(match: re.Match[str]) -> str:
        label = match.group("label")
        source = match.group("source")
        asset = figure_assets.get(source)
        if asset is None:
            raise BuildError(f"No staged asset exists for {source}")
        anchor = label.lower().replace(" ", "-")
        image = match.group(1).replace(f'src="{source}"', f'src="{asset}"')
        return (
            f'<figure class="display" id="{anchor}">\n'
            f"  {image}\n"
            f"  <figcaption><span class=\"fig-label\">{label}</span>{match.group('caption')}</figcaption>\n"
            "</figure>"
        )

    wrapped, count = figure_pattern.subn(replace, body_html)
    if count != len(images):
        raise BuildError("Figure image and caption pairs could not be assembled")
    return wrapped


def render_page(manuscript: Manuscript, template: str, body_html: str) -> str:
    """Fill the house template with an escaped masthead and rendered manuscript."""
    title_html = render_markdown(manuscript.title).removeprefix("<p>").removesuffix("</p>")
    replacements = {
        "{{TITLE_TEXT}}": escape(re.sub(r"<[^>]+>", "", title_html)),
        "{{TITLE_HTML}}": title_html,
        "{{AUTHORS}}": manuscript.authors_html,
        "{{AFFILIATIONS}}": manuscript.affiliations_html,
        "{{STANDFIRST}}": STANDFIRST,
        "{{BODY}}": body_html,
        "{{REPOSITORY_URL}}": REPOSITORY_URL,
        "{{CANONICAL_URL}}": CANONICAL_URL,
        "{{SOCIAL_PREVIEW_PATH}}": "assets/og.png",
    }
    page = template
    for token, value in replacements.items():
        page = page.replace(token, value)
    if re.search(r"{{[^}]+}}", page):
        raise BuildError("Template contains an unresolved token")
    return page


def _validate_rendered_page(page: str) -> None:
    """Reject incomplete output before it can replace the published manuscript."""
    if re.search(r"{{[^}]+}}", page):
        raise BuildError("Rendered page contains an unresolved template token")
    required = ("<!doctype html>", "<html", "<head>", "<body", "</html>")
    if any(token not in page.lower() for token in required):
        raise BuildError("Rendered page is not a complete HTML document")


def build(root: Path = ROOT, docs: Path = DOCS) -> Path:
    """Build validated manuscript HTML and transparent SVG assets into ``docs``."""
    output = docs / "index.html"
    root = root.resolve()
    docs = docs.resolve()
    manuscript_text = (root / "README.md").read_text(encoding="utf-8")
    manuscript = split_manuscript(manuscript_text)
    figures = collect_figures(manuscript_text, root)
    expected_labels = [
        "Figure 1", "Figure 2", "Figure 3", "Figure 4", "Figure 5A",
        "Figure 5B", "Figure 6", "Figure 7A", "Figure 7B",
    ]
    if [figure.label for figure in figures] != expected_labels:
        raise BuildError("Manuscript figures must use the required labels in order")

    staging = Path(tempfile.mkdtemp(prefix=".manuscript-build-", dir=docs))
    try:
        assets = stage_figure_assets(figures, staging, root)
        staged_figures = staging / "assets" / "figures"
        for figure in figures:
            ElementTree.parse(staged_figures / figure.asset_name)

        body_html = render_markdown(manuscript.body_markdown)
        body_html = wrap_figures(body_html, assets)
        body_html = rewrite_repository_links(body_html, root, REPOSITORY_URL)
        template = (docs / "template.html").read_text(encoding="utf-8")
        page = render_page(manuscript, template, body_html)
        _validate_rendered_page(page)
        staged_page = staging / "index.html"
        staged_page.write_text(page, encoding="utf-8")

        target_assets = docs / "assets" / "figures"
        target_assets.mkdir(parents=True, exist_ok=True)
        expected_names = {figure.asset_name for figure in figures}
        for old_asset in target_assets.glob("*.svg"):
            if old_asset.name not in expected_names:
                old_asset.unlink()
        for figure in figures:
            os.replace(staged_figures / figure.asset_name, target_assets / figure.asset_name)
        os.replace(staged_page, docs / "index.html")
        return output
    finally:
        shutil.rmtree(staging, ignore_errors=True)


if __name__ == "__main__":
    build()
