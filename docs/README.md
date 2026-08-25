# Common Densoviruses GitHub Pages manuscript

`docs/index.html` is the static, GitHub Pages version of the repository's
manuscript. The repository-root `README.md` remains the source of truth for
the manuscript text, captions, and figure references.

## Build and test

Run the automated tests from the repository root:

```bash
cd docs && pixi run test
```

Regenerate the committed Pages artifact after changing the manuscript,
template, site builder, or a source figure:

```bash
cd docs && pixi run build
```

The build writes `docs/index.html` and the nine SVG files in
`docs/assets/figures/`. These are generated files and are committed for
GitHub Pages to serve directly; do not edit them by hand. Regenerate and
commit them with the source change instead.

Each published figure is a validated transparent copy of its source SVG. The
builder removes only the single full-canvas white backdrop, leaving the
scientific source SVGs and their interior content unchanged. A build fails if
that safe transformation cannot be verified.

## Local preview

After building, serve the `docs/` directory:

```bash
cd docs && python -m http.server 8000
```

Open <http://localhost:8000>. No path rewriting is needed for local testing.

The interactive IGV viewers remain published at their existing routes:
`01-air-samples-tgs/`, `01-air-samples-vsp/`, `02-sra-mining/`,
`04-skin-microbiome/`, `05-cameroonian-plasma/`, and `06-belgian-air/`.
Keep these directories and their `index.html` files intact when maintaining
the manuscript landing page.

## Deployment

GitHub Pages deploys the committed `docs/` directory from `main` (`main:/docs`).
Build and commit the generated artifact before merging a Pages change; no
separate runtime application or deployment build step is required.
