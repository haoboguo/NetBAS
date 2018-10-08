# NetBAS
Both Human and Yeast PIN are used.
R mark down examples (/MarkDown) include:

Examples: using 50 permutations from yeast PIN (examples 1-6) and 100 permutations from human PIN (example 7), respectively.
1-3 can be combined in to a same mark down for the GO-GO matrix; 4-6 are for the GO enrichment for a gene set (the life-extending single-gene deletions, CM15, from Cell Metabolism 2015 paper); 7 is for the GO enrichment (or annotation refinement) for a single human gene PD-1 (2018 Nobel Prize in medicine).

Although only 50 or 100 permutations are used, the GO enrichment/suppression results are broadly consistent with those obtained from 1000 or 10k permutations.

All script below are used for the biological process (BP) GO terms. For cellular composition (CC) or molecular function (MF) terms, simply replace "bp" by "cc" or "mf" in the scripts.

Run the ms02star script for permutations. The success rate for yeast PIN is c.a. 85%, and that for human PIN is close to 100%.
Run the script ms02star/ms02.star.1k.R for 1k permutations, e.g.;
$ R -f ms02star.1k.R --args human.pin.csv human
will construct 1k permutations based on human PIN into the folder "./human"
  
1. yeast.bp.matrix.Rmd
This script calculate the frequencies between each pair of BP terms from the yeast PIN (the original BioGrid PIN)
The heatmap was plotted agains log10(freq+1) (the frequencies vary from 0 to > 7000).

2. yeast.ms02.bp.matrix.Rmd
This script do the same calculations as for 1, for permutations. 50 permutations are used here.
All results are saved in the directory "/ms02.yeast/bp-bp/"

3. yeast.bp.z.heatmap.Rmd
This script calculate the Z-scores for each BP term, resulting in a n*n matrix, where n is the total number of BP terms.
Heatmap is plotted based on this matrix, positive Z-scores are in red, and negative Z-scores are in blue, respectively.

4. yeast.cm15.bp.matrix.Rmd
This script calculate the BP terms from proteins that interact with the CM15 set (225 genes).
The result is a spectrum-like, 1D array.

5. yeast.cm15.ms02.bp.matrix.Rmd
This script do the same calculations as for 4, for permutations. 50 permutations are used here.
All results are saved in the directory "/ms02.yeast/cm15/"

6. yeast.cm15.bp.z.heatmap.Rmd
This script calculate z-scores of BP terms from CM15 set.
The enriched (Z > 3) and suppressed (Z < -3) BP terms are recorded in the files yeast.cm15.bp.enriched.csv and yeast.cm15.bp.suppressed.csv, respectively.
The top 10 enriched BP terms are plotted as heatmap.

7. human.pd1.annotation.Rmd
This script using human PIN for GO annotations of the protein PD-1 (id: PDCD1). This gene is one of the genes that function with the T-cells to block the immune reaction and is the target for immuno-therapy in cancer (Nobel Prize in Medicine 2018). 
Because only 100 permutations have been used, the sampling resulted to an Inf value term, which shall be discarded here.
Note that no suppressed GO terms have been found using the NetBAS method for the annotation of PD-1 gene.
Significantly enriched GO annotations for PDCD1 include "T cell receptor signaling pathway" (GO:0050852, z=40.757), "T cell differenciation" (GO:0030217, z=18.73), "T cell costimulation (GO:0031295, z=18.216), "T cell activation" (GO:0042110, z=9.496), in line with the role of this gene for regulating the T-cell functions (for which GO:0031295 is already in the known GO terms of the PD-1 gene).

HBG updated 10-8-2018
