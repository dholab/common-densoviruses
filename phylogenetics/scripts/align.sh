#!/usr/bin/env bash
set -euo pipefail

# Align NS1 amino acid sequences and build a rooted ML tree using the Parvovirinae set as an outgroup.
# Inputs:
#   - phylogenetics/parvovirinae_ns1_outgroup.fasta
#   - phylogenetics/vsp_ns1_aa_sequences.fasta
#   - phylogenetics/logan_ns1_aa_sequences.fasta
#   - phylogenetics/known densoviruses NCBI genbank 2025-11-19_ns1_aa_sequences.fasta
# Outputs:
#   - phylogenetics/alignments/ns1_alignment.fasta            (MAFFT alignment)
#   - phylogenetics/trees/ns1_tree.treefile                   (IQ-TREE ML tree, rooted on the outgroup)
#   - phylogenetics/trees/ns1_tree.iqtree/log files           (run metadata and supports)

THREADS=${THREADS:-12}
ALIGN_DIR="phylogenetics/alignments"
TREE_DIR="phylogenetics/trees"

mkdir -p "$ALIGN_DIR" "$TREE_DIR"

COMBINED="$ALIGN_DIR/ns1_combined.fasta"
FILTERED="$ALIGN_DIR/ns1_combined_filtered.fasta"
ALIGNMENT="$ALIGN_DIR/ns1_alignment.fasta"
PREFIX="$TREE_DIR/ns1_tree"

OUTGROUP="Erythroparvovirus_NS1|XZS47276.1,Dependoparvovirus_NS1|QKN88780.1,Protoparvovirus_NS1|YBV31529.1,Bocaparvovirus_NS1|YBQ66646.1"

# 1) Concatenate inputs.
cat \
  phylogenetics/parvovirinae_ns1_outgroup.fasta \
  phylogenetics/vsp_ns1_aa_sequences.fasta \
  phylogenetics/logan_ns1_aa_sequences.fasta \
  "phylogenetics/known densoviruses NCBI genbank 2025-11-19_ns1_aa_sequences.fasta" \
  > "$COMBINED"

# 2) Drop records that are all X/N/ambiguous to prevent downstream tree builder failures.
python - <<'PY'
from pathlib import Path

inp = Path("phylogenetics/alignments/ns1_combined.fasta")
outp = Path("phylogenetics/alignments/ns1_combined_filtered.fasta")
allowed = set("ACDEFGHIKLMNPQRSTVWY")
keep = []
with inp.open() as fh:
    header = None
    seq = []
    for line in fh:
        line = line.rstrip()
        if line.startswith(">"):
            if header is not None:
                s = "".join(seq).upper()
                if s and sum(c in allowed for c in s) / len(s) >= 0.1:
                    keep.append((header, s))
            header, seq = line, []
        else:
            seq.append(line.strip())
    if header is not None:
        s = "".join(seq).upper()
        if s and sum(c in allowed for c in s) / len(s) >= 0.1:
            keep.append((header, s))
with outp.open("w") as out:
    for h, s in keep:
        out.write(f"{h}\n")
        for i in range(0, len(s), 80):
            out.write(s[i:i+80] + "\n")
PY

# 3) Align with MAFFT.
mafft --thread "$THREADS" --auto "$FILTERED" > "$ALIGNMENT"

# 4) Infer ML tree with IQ-TREE, rooting on the Parvovirinae outgroup.
iqtree \
  -s "$ALIGNMENT" \
  -m LG+G4+F \
  -bb 1000 \
  -alrt 1000 \
  -nt "$THREADS" \
  -safe \
  -o "$OUTGROUP" \
  -pre "$PREFIX"

echo "Done. Alignment: $ALIGNMENT"
echo "Tree: ${PREFIX}.treefile"
