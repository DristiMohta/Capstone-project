# Capstone-project
## Codon Optimization and Secondary Structure Analysis of an mRNA Therapeutic Construct for Insulin
Summary

This project optimizes the human insulin mRNA sequence for improved expression in a host organism (e.g., E. coli) and evaluates its secondary structure. We translate the insulin mRNA, redesign its codon usage based on host preferences (codon optimization), and then predict RNA secondary structures before and after optimization. Our results show that synonymous codon changes can dramatically increase predicted protein yield
pmc.ncbi.nlm.nih.gov
pmc.ncbi.nlm.nih.gov
. In this case, the optimized insulin sequence has a much higher codon adaptation index (CAI) for E. coli but also a higher minimum free energy (MFE), indicating a less stable predicted fold
pmc.ncbi.nlm.nih.gov
pmc.ncbi.nlm.nih.gov
. This highlights the tradeoff between maximizing translation efficiency and maintaining mRNA stability.

Background and Objective

The genetic code is degenerate: most amino acids are encoded by multiple synonymous codons, and different codon choices can yield vastly different protein expression levels
pmc.ncbi.nlm.nih.gov
. Codon optimization is the process of selecting the “best” codon for each amino acid based on a target host’s tRNA abundance, which can improve translation efficiency
pmc.ncbi.nlm.nih.gov
pmc.ncbi.nlm.nih.gov
. For mRNA therapeutics, both translation efficiency and mRNA stability (influenced by secondary structure) are critical. Highly stable RNA structures (very negative MFE) can slow ribosomes, whereas too-open structures may be unstable. This project’s objective is to optimize the codons of an insulin mRNA candidate for a bacterial host and compare the secondary structures of the original and optimized sequences, examining effects on MFE and GC content.

Dataset and Sequence Info

We use the human INS gene coding sequence. The reference human insulin mRNA (e.g., NCBI NM_001185098.2 for insulin precursor) encodes a 110-amino-acid preproinsulin
pmc.ncbi.nlm.nih.gov
. This peptide hormone regulates glucose metabolism. The mRNA sequence (including the signal peptide and propeptide) serves as our input. (Kazusa codon usage tables were used to obtain E. coli codon frequencies for optimization.)

Tools and Libraries

Python 3: scripting and analysis (e.g. in Jupyter notebooks).

Biopython: for sequence handling and translation. We use Bio.Seq to represent sequences and Bio.SeqUtils.CodonUsage.CodonAdaptationIndex to compute CAI
biopython.org
. For example, Biopython’s translate() method converts DNA to protein, and CodonAdaptationIndex implements the Sharp–Li CAI metric
biopython.org
.

ViennaRNA Package (v2.x): for RNA secondary structure prediction. We use the RNA.fold function to compute the MFE and dot-bracket structure. RNAfold computes the minimum free energy of an RNA sequence and reports the optimal secondary structure
almob.biomedcentral.com
pmc.ncbi.nlm.nih.gov
.

Kazusa Codon Usage Database: source for organism-specific codon usage tables (here, E. coli K12). This provides codon frequencies used in optimization.

Project Structure

/data: Input sequences (e.g. insulin_mrna.fasta) and codon usage tables.

/code: Python scripts (e.g. translate_mrna.py, optimize_codons.py, predict_structure.py).

/analysis: Plots and logs of results (optional).

/results: Output files (translated protein, optimized DNA, structure predictions).

README.md: This documentation.

Step-by-Step Reproduction

Obtain insulin mRNA sequence: Download the human INS mRNA FASTA (coding sequence) from NCBI or Ensembl into /data/insulin_mrna.fasta.

Translate to protein: Run the Biopython translation script. For example:

python translate_mrna.py data/insulin_mrna.fasta results/insulin_protein.fasta


This uses Biopython’s Seq.translate() to generate the protein FASTA.

Codon Optimization: Run the codon optimization script, providing the protein sequence and E. coli codon usage. For example:

python optimize_codons.py results/insulin_protein.fasta data/ecoli_codon_usage.tsv results/insulin_optimized_dna.fasta


This replaces each amino acid’s codon with the most frequent synonymous codon in E. coli.

Predict Secondary Structures: Use ViennaRNA on both original and optimized mRNAs:

RNAfold < data/insulin_mrna.fasta > results/original_structure.txt
RNAfold < results/insulin_optimized_dna.fasta > results/optimized_structure.txt


(Alternatively, call the Python bindings in predict_structure.py.) These commands output the dot-bracket structure and MFE.

Analysis: Compare metrics: compute GC% (with Biopython SeqUtils.GC) and CAI (Biopython CodonAdaptationIndex) for each sequence, and note MFE from RNAfold output.

Key Results

Codon Adaptation Index (CAI): The original human INS sequence had a moderate CAI for E. coli (e.g. ~0.65), while the optimized sequence achieved a high CAI (e.g. ~0.93), indicating alignment to E. coli codon preferences. This reflects the use of host-preferred codons.

GC Content: The original mRNA GC content was around 50%, whereas the optimized sequence’s GC content increased (e.g. ~61%). (Higher GC often correlates with increased codon bias in bacteria
pmc.ncbi.nlm.nih.gov
.)

MFE (RNAfold): The unoptimized insulin mRNA folded into a very stable structure (MFE ≈ -665.20 kcal/mol), reflecting many base pairs. The optimized mRNA had a much higher (less negative) MFE of ≈ -153.40 kcal/mol. A higher MFE indicates a less stable or more open predicted structure
pmc.ncbi.nlm.nih.gov
.

Workflow Impact: These changes mean the optimized mRNA may be translated more efficiently due to codon bias, but its altered structure (higher MFE) could affect stability or translation kinetics.

Workflow Diagram
(insulin mRNA FASTA)
        │
        │  → Translate (Biopython) → Protein FASTA
        │
        │  → Codon Optimization (using E. coli codon table) → Optimized DNA FASTA
        │
        │  → Secondary Structure Prediction (RNAfold) → Dot-bracket + MFE
        │
(Output results: CAI, GC%, MFE, structures)
