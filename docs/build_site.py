"""Helpers for building the manuscript GitHub Pages site from Markdown."""

from dataclasses import dataclass
from html import unescape
from pathlib import Path
import re
from urllib.parse import urlsplit

import markdown


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
