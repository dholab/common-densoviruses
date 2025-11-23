# Common Densoviruses - GitHub Pages Setup

This directory contains the GitHub Pages website for viewing genomic alignment data using IGV.js.

## Setup Instructions

### 1. Enable GitHub Pages

1. Go to your repository on GitHub: https://github.com/dholab/common-densoviruses
2. Click on "Settings" tab
3. Scroll down to "Pages" in the left sidebar
4. Under "Source", select "Deploy from a branch"
5. Under "Branch", select "main" and "/docs" folder
6. Click "Save"

Your site will be available at: https://dholab.github.io/common-densoviruses/

### 2. Index BAM Files

Before the IGV viewer can display BAM files, they need to be indexed:

```bash
# Run from the repository root
./scripts/index_bams.sh
```

This creates `.bai` index files for all BAM files.

### 3. Commit and Push Changes

```bash
git add .
git commit -m "Add GitHub Pages with IGV viewer"
git push origin main
```

### 4. Wait for Deployment

GitHub Pages typically takes 5-10 minutes to deploy. You can check the deployment status in the "Actions" tab of your repository.

## File Structure

- `index.html` - Main IGV viewer page
- `_config.yml` - GitHub Pages configuration
- `../air-samples/vsp/bam/` - BAM alignment files
- `../air-samples/vsp/ref/` - Reference genome and annotations

## Troubleshooting

### CORS Issues

If you encounter CORS (Cross-Origin Resource Sharing) errors:
1. Make sure all files are committed and pushed to GitHub
2. Access files through the raw.githubusercontent.com URL
3. Consider using GitHub's Large File Storage (LFS) for BAM files over 100MB

### BAM Files Not Loading

1. Ensure BAM files have corresponding `.bai` index files
2. Check that file paths in the HTML match your repository structure
3. Verify files are accessible at: https://raw.githubusercontent.com/dholab/common-densoviruses/main/[path-to-file]

### Large File Issues

GitHub has a 100MB file size limit. For larger BAM files:
1. Use Git LFS (Large File Storage)
2. Host files on a separate data server
3. Consider down-sampling BAM files for web viewing

## Local Development

To test locally:
```bash
# Start a local web server in the docs directory
python3 -m http.server 8000

# Open in browser
# http://localhost:8000
```

Note: When testing locally, you'll need to adjust the base URL in the JavaScript code to point to local files.