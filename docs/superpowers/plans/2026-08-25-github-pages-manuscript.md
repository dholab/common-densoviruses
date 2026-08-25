# GitHub Pages Manuscript Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Generate and commit a GitHub Pages manuscript from `README.md` that matches the Cyclospora house style, preserves existing IGV routes, and embeds transparent copies of all nine scientific SVG figures.

**Architecture:** A Python build script reads the manuscript Markdown, extracts its masthead, renders the body, wraps figure/caption pairs, rewrites repository-relative links, and stages validated transparent SVG copies before atomically publishing `docs/index.html` and `docs/assets/figures/`. A standalone HTML template owns the shared visual system and metadata. Unit tests exercise parsing, link resolution, safe background removal, and full-page generation with synthetic fixtures before the real manuscript is generated.

**Tech Stack:** Python 3.12, Python-Markdown 3.x, standard-library `unittest`, `xml.etree.ElementTree`, Pixi, static HTML/CSS, GitHub Pages.

**Spec:** `docs/superpowers/specs/2026-08-25-github-pages-manuscript-design.md`

## Global Constraints

- Background is exactly `#F8F4E9`; primary text is `#163139`; accent links and figure labels are `#C16A3C`.
- Headings use `Georgia, "Iowan Old Style", "Times New Roman", serif`; body text uses `Arial, Helvetica, "Liberation Sans", sans-serif`.
- Reading measure is exactly `42rem`; body typography is `16px/1.62` before the narrow-screen adjustment.
- `README.md` is the canonical manuscript source, and `docs/index.html` is committed generated output.
- Original SVGs, BAM/BAI files, reference data, existing IGV pages, and `analysis/XX-sra-mining-acheta-domestica/` must not be modified.
- A transparent SVG transformation may remove only an explicitly white full-canvas background matching the SVG's declared width and height.
- GitHub Pages continues to deploy the `main` branch's `/docs` directory; no alternate hosting configuration is added.
- The build must fail on missing manuscript structure, unmatched figure captions, missing SVGs, unsafe background detection, invalid generated SVG, or unresolved local asset links.

---

### Task 1: Establish the tested manuscript builder

**Files:**
- Create: `docs/pixi.toml`
- Create: `docs/build_site.py`
- Create: `docs/tests/__init__.py`
- Create: `docs/tests/test_build_site.py`
- Generate: `docs/pixi.lock`

**Interfaces:**
- Consumes: Markdown text, repository root `Path`, and repository URL.
- Produces: `Manuscript(title: str, authors_html: str, affiliations_html: str, body_markdown: str)`, `split_manuscript(text: str) -> Manuscript`, `render_markdown(text: str) -> str`, `rewrite_repository_links(html_text: str, root: Path, repo_url: str, branch: str = "main") -> str`.

- [ ] **Step 1: Add the failing parser and link-rewrite tests**

Create synthetic manuscript text with a title, authors paragraph, affiliation paragraph, `## Introduction`, one local file link, one local directory link, one external link, one mail link, and one anchor. Assert that `split_manuscript()` extracts the three masthead fields, retains `## Introduction` in the body, rewrites file and directory links to `/blob/main/` and `/tree/main/`, and leaves the external/mail/anchor links unchanged.

```python
class ManuscriptParsingTests(unittest.TestCase):
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
```

- [ ] **Step 2: Run the focused tests and confirm the expected import failure**

Run: `cd docs && pixi run python -m unittest tests.test_build_site.ManuscriptParsingTests -v`

Expected: FAIL because `build_site.py` and its public interfaces do not exist.

- [ ] **Step 3: Add the Pixi environment and minimal parsing implementation**

Use this project contract in `docs/pixi.toml`:

```toml
[project]
name = "common-densoviruses-manuscript-site"
version = "0.1.0"
description = "Builds the GitHub Pages manuscript from README.md"
channels = ["conda-forge"]
platforms = ["osx-arm64", "linux-64"]

[tasks]
build = "python build_site.py"
test = "python -m unittest discover -s tests -p 'test_*.py'"

[dependencies]
python = ">=3.12,<3.13"
markdown = ">=3.6,<4"
```

Implement `Manuscript` as a frozen dataclass. Require an H1, two nonempty front-matter blocks, and a first `##` section. Render front-matter blocks with Python-Markdown and the `nl2br` extension, and render body Markdown with `tables`, `fenced_code`, `attr_list`, `sane_lists`, and `nl2br`. Decode HTML character references before classifying link schemes; use `Path.is_dir()` to choose `tree`, otherwise `blob`.

- [ ] **Step 4: Run the focused parser tests**

Run: `cd docs && pixi run python -m unittest tests.test_build_site.ManuscriptParsingTests -v`

Expected: PASS.

- [ ] **Step 5: Generate and inspect the lock file**

Run: `cd docs && pixi lock`

Expected: `docs/pixi.lock` resolves Python 3.12 and Python-Markdown for both declared platforms without dependency errors.

- [ ] **Step 6: Commit the tested builder foundation**

```bash
git add docs/pixi.toml docs/pixi.lock docs/build_site.py docs/tests/__init__.py docs/tests/test_build_site.py
git commit -m "Add tested manuscript site builder"
```

---

### Task 2: Add safe transparent SVG generation

**Files:**
- Modify: `docs/build_site.py`
- Modify: `docs/tests/test_build_site.py`

**Interfaces:**
- Consumes: SVG source/destination `Path` objects and Markdown figure references.
- Produces: `Figure(label: str, source: Path, asset_name: str)`, `collect_figures(markdown_text: str, root: Path) -> list[Figure]`, `make_svg_transparent(source: Path, destination: Path) -> None`, `stage_figure_assets(figures: list[Figure], output_dir: Path) -> dict[str, str]` mapping repository-relative source paths to `assets/figures/<asset_name>`.

- [ ] **Step 1: Add failing SVG safety tests**

Create three fixture SVGs in temporary directories: one with a direct white full-canvas rectangle, one with a white full-canvas rectangle inside a white background group, and one with only a smaller interior white rectangle. Assert that the first two transformations remove exactly one background, keep colored shapes, and parse successfully; assert that the third transformation raises `BuildError` and does not publish a destination file.

```python
def test_removes_only_direct_full_canvas_white_rect(self):
    source = self.write_svg(
        '<svg xmlns="http://www.w3.org/2000/svg" width="100" height="60">'
        '<rect x="0" y="0" width="100" height="60" fill="white"/>'
        '<rect x="10" y="10" width="20" height="20" fill="#163139"/>'
        '</svg>'
    )
    make_svg_transparent(source, self.destination)
    root = ElementTree.parse(self.destination).getroot()
    self.assertEqual(len(list(root)), 1)
    self.assertEqual(list(root)[0].attrib["fill"], "#163139")

def test_rejects_svg_without_expected_canvas_background(self):
    source = self.write_svg(
        '<svg xmlns="http://www.w3.org/2000/svg" width="100" height="60">'
        '<rect x="10" y="10" width="20" height="20" fill="white"/>'
        '</svg>'
    )
    with self.assertRaises(BuildError):
        make_svg_transparent(source, self.destination)
    self.assertFalse(self.destination.exists())
```

- [ ] **Step 2: Run the SVG tests and verify failure**

Run: `cd docs && pixi run python -m unittest tests.test_build_site.TransparentSvgTests -v`

Expected: FAIL because `make_svg_transparent()` and related types do not exist.

- [ ] **Step 3: Implement exact background detection and staged publishing**

Parse numeric `width` and `height`, accepting a trailing `px` but rejecting percentages or missing dimensions. Traverse with parent references. A candidate must be a `<rect>` with `x=0`, `y=0`, width and height numerically equal to the canvas, and effective white fill from either `fill="white"`, `fill="#fff"`, `fill="#ffffff"`, `fill="rgb(255,255,255)"`, or an inherited parent `style` declaration. Reject zero or multiple candidates. Remove the rectangle, or remove its parent when that parent is a background-only group. Serialize to a temporary file, parse the result again, then rename it to the destination.

`collect_figures()` must recognize the labels `Figure 1`, `Figure 2`, `Figure 3`, `Figure 4`, `Figure 5A`, `Figure 5B`, `Figure 6`, `Figure 7A`, and `Figure 7B`, reject duplicate labels or sources, and generate lowercase hyphenated asset names such as `figure-5a.svg`.

- [ ] **Step 4: Run SVG and full unit tests**

Run: `cd docs && pixi run test`

Expected: PASS, including successful parsing of every generated SVG fixture.

- [ ] **Step 5: Commit transparent asset processing**

```bash
git add docs/build_site.py docs/tests/test_build_site.py
git commit -m "Generate safe transparent figure assets"
```

---

### Task 3: Build the house-style manuscript page

**Files:**
- Modify: `README.md`
- Create: `docs/template.html`
- Modify: `docs/build_site.py`
- Modify: `docs/tests/test_build_site.py`
- Create: `docs/.nojekyll`

**Interfaces:**
- Consumes: `Manuscript`, rendered body HTML, figure source-to-asset mapping, and the template.
- Produces: `wrap_figures(body_html: str, figure_assets: dict[str, str]) -> str`, `render_page(manuscript: Manuscript, template: str, body_html: str) -> str`, and `build(root: Path = ROOT, docs: Path = DOCS) -> Path` returning the generated `docs/index.html` path.

- [ ] **Step 1: Add failing figure-wrapper and page-integration tests**

Assert that a rendered image immediately followed by a blockquote caption becomes one `<figure class="display" id="figure-5a">`, that the source changes to `assets/figures/figure-5a.svg`, that the caption label is `Figure 5A`, and that unmatched images or captions raise `BuildError`. Add a full fixture build test that asserts the title, masthead, `## Introduction`, standfirst, colophon, theme color, canonical URL, nine figure anchors, and zero unresolved `{{...}}` template tokens.

```python
def test_wraps_image_and_caption_as_labeled_figure(self):
    body = (
        '<p><img alt="Figure 5A" src="docs/source.svg" /></p>'
        '<blockquote><p><strong>Figure 5A.</strong> Evidence caption.</p></blockquote>'
    )
    wrapped = wrap_figures(body, {"docs/source.svg": "assets/figures/figure-5a.svg"})
    self.assertIn('<figure class="display" id="figure-5a">', wrapped)
    self.assertIn('src="assets/figures/figure-5a.svg"', wrapped)
    self.assertIn('<span class="fig-label">Figure 5A</span>', wrapped)
    self.assertNotIn("<blockquote>", wrapped)
```

- [ ] **Step 2: Run the page tests and verify failure**

Run: `cd docs && pixi run python -m unittest tests.test_build_site.PageRenderingTests -v`

Expected: FAIL because the template and page assembly interfaces do not exist.

- [ ] **Step 3: Reorder the README into manuscript form**

Move the existing author paragraph and affiliation lines directly below the H1, delete their former end section and the future-coauthor status sentence, and insert `## Introduction` before the opening densovirus background paragraph. Do not edit scientific prose, citations, captions, or figure paths in this step.

- [ ] **Step 4: Add the exact house-style template**

Port the reference site's CSS tokens, masthead, typography, links, code blocks, display-item rules, colophon, and narrow-screen media query. Remove the Cyclospora-specific PDF control and interactive Vega scripts. Add placeholders for escaped title/description, title HTML, authors, affiliations, standfirst, body, repository URL, canonical URL, and social-preview path.

The standfirst copy is:

> This web version is the preferred way to read this evolving manuscript because each figure links directly to its interactive genomic viewer, methods, scripts, and source data in the repository.

- [ ] **Step 5: Implement figure wrapping and atomic full-site generation**

Render all outputs inside a temporary directory under `docs/`. Validate figure count and label order, reparse every SVG, render and validate the complete HTML, then replace only `docs/index.html` and `docs/assets/figures/`. Preserve `docs/assets/og.png` when replacing the generated figures directory. Use `os.replace()` for the HTML file and generated SVG files after validation.

- [ ] **Step 6: Run focused and complete tests**

Run: `cd docs && pixi run python -m unittest tests.test_build_site.PageRenderingTests -v`

Expected: PASS.

Run: `cd docs && pixi run test`

Expected: PASS.

- [ ] **Step 7: Commit the manuscript source and template**

```bash
git add README.md docs/.nojekyll docs/template.html docs/build_site.py docs/tests/test_build_site.py
git commit -m "Add house-style manuscript rendering"
```

---

### Task 4: Add the branded social preview and metadata

**Files:**
- Create: `docs/assets/og.png`
- Modify: `docs/template.html`
- Modify: `docs/tests/test_build_site.py`

**Interfaces:**
- Consumes: the approved house palette, exact manuscript title, canonical URL, and generated `docs/assets/og.png`.
- Produces: a 1200×630 landscape preview plus matching Open Graph and X metadata in every generated root page.

- [ ] **Step 1: Add a failing metadata test**

Assert that generated HTML includes canonical, Open Graph title/description/type/url/image, and X card/title/description/image fields; assert that the image URL is the trusted absolute URL `https://dholab.github.io/common-densoviruses/assets/og.png` and the local file exists.

- [ ] **Step 2: Run the metadata test and verify failure**

Run: `cd docs && pixi run python -m unittest tests.test_build_site.PageRenderingTests.test_social_metadata_is_complete -v`

Expected: FAIL because the social card and complete metadata are absent.

- [ ] **Step 3: Generate and inspect one branded card**

Use one image-generation request for a 1200×630 editorial card with the exact title `Common densoviruses in the human and mammalian virospheres`, cream background, deep-teal serif title, restrained terracotta rule/accent, generous whitespace, and no invented scientific imagery, logos, people, credentials, or private data. Save the inspected result as `docs/assets/og.png`; retry once only if title text is incorrect or unreadable.

- [ ] **Step 4: Wire metadata and pass the focused test**

Add the tested canonical and sharing tags to `docs/template.html`, rerun the focused metadata test, and confirm PASS.

- [ ] **Step 5: Commit the preview and metadata**

```bash
git add docs/assets/og.png docs/template.html docs/tests/test_build_site.py
git commit -m "Add manuscript social preview"
```

---

### Task 5: Generate the real Pages artifact and documentation

**Files:**
- Replace: `docs/index.html`
- Create: `docs/assets/figures/figure-1.svg`
- Create: `docs/assets/figures/figure-2.svg`
- Create: `docs/assets/figures/figure-3.svg`
- Create: `docs/assets/figures/figure-4.svg`
- Create: `docs/assets/figures/figure-5a.svg`
- Create: `docs/assets/figures/figure-5b.svg`
- Create: `docs/assets/figures/figure-6.svg`
- Create: `docs/assets/figures/figure-7a.svg`
- Create: `docs/assets/figures/figure-7b.svg`
- Replace: `docs/README.md`
- Modify: `docs/tests/test_build_site.py`

**Interfaces:**
- Consumes: real repository README, SVG sources, template, social card, and builder.
- Produces: committed GitHub Pages files and accurate local build/deployment instructions.

- [ ] **Step 1: Add a real-repository contract test**

Assert that `collect_figures(README.md)` returns the exact nine labels in order, every source exists, every expected existing IGV route still contains `index.html`, and `build()` creates a root page containing exactly nine figure displays.

- [ ] **Step 2: Run the contract test before generation**

Run: `cd docs && pixi run python -m unittest tests.test_build_site.RepositoryContractTests -v`

Expected: FAIL until the final build output and generated figure assets exist.

- [ ] **Step 3: Generate the real site twice**

Run: `cd docs && pixi run build`

Expected: the command reports nine transparent SVG assets and writes `docs/index.html`.

Run from the repository root:

```bash
FIRST_BUILD_CHECKSUMS=$(shasum docs/index.html docs/assets/figures/*.svg)
(cd docs && pixi run build)
SECOND_BUILD_CHECKSUMS=$(shasum docs/index.html docs/assets/figures/*.svg)
test "$FIRST_BUILD_CHECKSUMS" = "$SECOND_BUILD_CHECKSUMS"
```

Expected: `test` exits 0, proving the second build is identical.

- [ ] **Step 4: Replace the Pages README with accurate maintenance instructions**

Document `cd docs && pixi run test`, `cd docs && pixi run build`, the generated-file policy, transparent-copy behavior, local `python -m http.server 8000`, preservation of the interactive viewer routes, and deployment from `main:/docs`. Remove obsolete dataset-card descriptions and advice to rewrite paths for local testing.

- [ ] **Step 5: Run the complete automated suite**

Run: `cd docs && pixi run test`

Expected: PASS for parser, links, SVG safety, page rendering, metadata, and real-repository contracts.

- [ ] **Step 6: Commit the generated Pages site**

```bash
git add docs/index.html docs/assets/figures docs/README.md docs/tests/test_build_site.py
git commit -m "Generate GitHub Pages manuscript"
```

---

### Task 6: Validate the rendered site and repository integrity

**Files:**
- Modify only if validation reveals a defect: `README.md`, `docs/template.html`, `docs/build_site.py`, `docs/tests/test_build_site.py`, `docs/index.html`, `docs/assets/figures/*`, `docs/assets/og.png`, `docs/README.md`

**Interfaces:**
- Consumes: final committed build.
- Produces: evidence that the root manuscript, representative interactive viewers, responsive layout, and repository boundaries satisfy the specification.

- [ ] **Step 1: Run clean automated verification**

Run: `cd docs && pixi run test && pixi run build && cd .. && git diff --exit-code -- docs/index.html docs/assets/figures`

Expected: tests pass, build succeeds, and generated files have no post-build diff.

- [ ] **Step 2: Serve the `docs/` directory and check HTTP routes**

Start `python3 -m http.server 8000 --directory docs` in a retained session. Confirm HTTP 200 for `/`, `/assets/og.png`, all nine `/assets/figures/*.svg` files, `/01-air-samples-vsp/`, and `/06-belgian-air/`.

- [ ] **Step 3: Perform browser layout validation**

Inspect the local root page at desktop and narrow mobile widths. Compare the background, palette, font families, 42rem measure, masthead spacing, heading hierarchy, link treatment, figure rules, caption labels, and responsive padding against the Cyclospora reference. Confirm the cream background is visible through every figure canvas, the tall IGV exports are not cropped, and no horizontal page overflow occurs.

- [ ] **Step 4: Check HTML, SVG, and link integrity**

Parse `docs/index.html` with an HTML parser; parse every generated SVG with `ElementTree`; scan `href` and `src` attributes for missing local targets; confirm there are no unresolved template tokens, `file:` URLs, temporary paths, or references to source SVG paths outside `docs/`.

- [ ] **Step 5: Audit the final Git scope**

Run: `git status --short`, `git diff --check HEAD~5..HEAD`, and `git diff --name-only 367852e..HEAD`.

Expected: only the planned README and `docs/` files changed after the design commit; `analysis/XX-sra-mining-acheta-domestica/` remains untracked and untouched.

- [ ] **Step 6: Commit only validation fixes, if any**

If validation required source changes, regenerate the site, rerun the complete verification command, and commit the exact affected files with:

```bash
git commit -m "Fix GitHub Pages validation issues" -- README.md docs/
```

If validation required no changes, do not create an empty commit.
