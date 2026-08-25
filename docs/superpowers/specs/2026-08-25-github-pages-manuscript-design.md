# GitHub Pages manuscript design

## Purpose

Replace the current dataset-card landing page with a manuscript reading version of `README.md` that matches the visual system used by the Cyclospora wastewater metagenomics GitHub Pages site. The repository README remains the canonical manuscript source, while the generated Pages document provides a more deliberate reading experience and preserves access to the existing interactive genomic viewers.

## Goals

- Match the reference site's house style: cream background, deep-teal text, terracotta accents, Georgia-family display type, Arial-family body type, narrow reading measure, restrained rules, and compact figure captions.
- Present the existing authors and affiliations in the masthead immediately below the title.
- Keep manuscript content synchronized by generating `docs/index.html` from `README.md`.
- Give every embedded scientific figure a transparent canvas so the cream page background shows through.
- Preserve the scientific content of every figure and leave the original figure files unchanged.
- Preserve all existing interactive IGV routes and make them reachable from the manuscript links.
- Produce a static site that GitHub Pages can deploy from the existing `docs/` directory without a runtime application or database.

## Non-goals

- Rewriting the scientific argument, changing figure data, or redesigning the interactive IGV viewers.
- Adding manuscript functionality that is absent from the source repository, such as an interactive chart or downloadable PDF.
- Moving or modifying BAM, BAI, reference, or other scientific data files.
- Changing GitHub Pages hosting settings or publishing through another hosting provider.

## Source structure

`README.md` will be adjusted to follow manuscript order:

1. Title
2. Authors
3. Affiliations and corresponding-author line
4. Introductory manuscript text
5. Results and subordinate sections
6. Discussion
7. Data availability
8. Acknowledgements

The existing `Authors and affiliations` section near the end will therefore move into the front matter. The temporary note about future co-authors will be removed because it is project-status text rather than manuscript content.

## Build architecture

The site will use the same generator-and-template pattern as the reference repository:

- `docs/build_site.py` reads `README.md`, renders Markdown to HTML, structures the masthead and figures, rewrites repository links, and writes `docs/index.html`.
- `docs/template.html` contains the document shell, metadata, and house-style CSS.
- `docs/assets/` contains generated, Pages-ready figure copies and any site-wide social-preview asset.
- `docs/index.html` is generated output and is committed so GitHub Pages can serve it directly.
- `docs/.nojekyll` prevents Jekyll from modifying the generated static output.

The build remains deterministic: identical README, template, generator, and source SVGs must produce identical Pages assets and HTML.

## Manuscript rendering

The generator will:

- extract the H1 title and front-matter author blocks for the masthead;
- add a short standfirst explaining that the web manuscript links directly to its source repository and interactive viewers;
- render `##` and `###` headings as manuscript sections and subsections;
- recognize each Markdown image followed by its blockquote caption as one figure display;
- preserve inline emphasis, code, links, lists, and other Markdown constructs;
- assign stable figure anchors derived from figure labels such as `Figure 5A`;
- rewrite repository-relative links that cannot be served from `docs/` to their corresponding GitHub `blob` or `tree` URLs;
- keep external, mail, anchor, and existing GitHub Pages links unchanged; and
- include a colophon linking back to `README.md` and the generator.

The current dataset-card content in `docs/index.html` will be replaced. The existing viewer directories under `docs/01-air-samples-vsp/`, `docs/01-air-samples-tgs/`, `docs/02-sra-mining/`, `docs/04-skin-microbiome/`, `docs/05-cameroonian-plasma/`, and `docs/06-belgian-air/` remain unchanged and continue to be linked from figure captions.

## Visual system

The template will reuse the reference site's core tokens and layout:

- Background: `#F8F4E9`
- Text and rules: deep teal `#163139` with transparent tints
- Links and figure labels: terracotta `#C16A3C`
- Heading type: Georgia, Iowan Old Style, or Times New Roman fallback
- Body type: Arial, Helvetica, or Liberation Sans fallback
- Reading measure: `42rem`
- Body size and leading: `16px/1.62`

Figures will be bounded by subtle horizontal rules. Captions will use a compact sans-serif style with terracotta uppercase figure labels. The layout will remain single-column and responsive, reducing typography and padding on narrow screens without introducing a separate mobile information architecture.

## Transparent figure processing

All nine embedded figures are SVGs. The generator will create Pages-ready copies rather than altering their scientific source files.

For each SVG, the transformation will remove only a full-canvas background shape that is explicitly white and covers the declared SVG width and height. This includes the white root rectangle used by IGV exports and the equivalent white background group used by the two tree exports. Interior white shapes, labels, read tracks, and other graphic elements will remain intact.

The build will fail if it cannot identify the expected full-canvas background in a source figure or if the transformed SVG is invalid. This prevents a silent broad color replacement that could damage scientific content. Generated figure copies will live under `docs/assets/figures/` and the generated HTML will reference those copies.

## Metadata and sharing

The document head will include the manuscript title, a concise description based on the README's subject, the cream theme color, canonical Pages URL, and standard Open Graph/X fields. A branded landscape social-preview image will be added using the same palette and typography; it will not substitute for or modify any scientific figure.

## Accessibility

- Existing Markdown image alt text will be preserved and improved only where it is too generic to identify the figure.
- Heading hierarchy will remain sequential.
- Links and controls will have visible keyboard focus.
- Text and accent colors will maintain readable contrast against the cream background.
- Wide figures will scale within the reading column; unusually tall IGV exports will remain fully available without cropping.
- The static page will remain usable when scripts are disabled.

## Error handling

The build exits unsuccessfully when required source sections are missing, a figure-caption pair cannot be recognized, an expected image is unavailable, SVG background processing is unsafe, or generated HTML contains unresolved repository-relative asset paths. Existing files outside the generated site surface are not changed on a failed build.

## Verification

Verification will include:

1. Run the generator twice and confirm the second run produces no diff.
2. Parse the generated HTML and all generated SVGs successfully.
3. Confirm that all nine figures appear once, use generated transparent assets, and retain their captions.
4. Confirm that each transparent SVG has no full-canvas white background while retaining interior white elements.
5. Check all local links and referenced assets in `docs/index.html`.
6. Serve `docs/` locally and verify that the manuscript and representative IGV routes return successfully.
7. Compare the rendered manuscript with the reference site's typography, width, palette, spacing, figure treatment, and responsive layout at desktop and mobile widths.
8. Run the repository's existing relevant checks, if present, and inspect the final Git diff to ensure unrelated data and the existing untracked analysis directory are untouched.

## Deployment

The repository will continue using GitHub Pages from the `main` branch's `/docs` directory. No deployment will be triggered from the local implementation. Once the changes are committed and pushed by the repository owner, GitHub Pages will publish the generated manuscript at `https://dholab.github.io/common-densoviruses/`.
