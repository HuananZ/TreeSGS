Detecting Population-differentiation Copy Number Variants in Human Population Tree by Sparse Group Selection

Huanan Zhang1, David Roe2,3 and Rui Kuang1,3

1. Department of Computer Science and Engineering, University of Minnesota Twin Cities, Minneapolis, MN
2. Bioinformatics Research, National Marrow Donor Program, Minneapolis, MN
3. Biomedical Informatics and Computational Biology, Minneapolis, MN


ABSTRACT

Copy-number variants (CNVs) account for a substantial proportion of human genetic variations. Understanding the CNV diversities across populations is a computational challenge because CNV patterns are often present in several related populations and only occur in a subgroup of individuals within each of the population. This paper introduces a tree-guided sparse group selection algorithm (treeSGS) to detect population-differentiation CNV markers of subgroups across populations organized by a phylogenetic tree of human populations. The treeSGS algorithm detects CNV markers of populations associated with nodes from all levels of the tree such that the evolutionary relations among the populations are incorporated for more accurate detection of population-differentiation CNV. We applied treeSGS algorithm to study the 1179 samples from the 11 populations in Hapmap3 CNV genotype data. The treeSGS algorithm accurately identifies CNV markers of each population and the collection of populations organized under the branches of the human population tree, validated by consistency among family trios and SNP characterizations of the CNV regions. Further comparison between the detected CNV markers and other population-differentiation CNVs reported in 1000 genome data and other recent studies also shows that treeSGS can significantly improve the current annotations of population-differentiation CNV markers. Finally, a list of population-differentiation CNVs that are potential Ancestry Informative Markers (AIMs) is reported.


Contact: 
huanan@cs.umn.edu
kuang@cs.umn.edu

Funding: NSF-IIS 1149697


How to use:

run_TBL_b.m	the main script to run treeSGS_b algorithm on Hapmap 3 data.
run_TBL_g.m	the main script to run treeSGS_g algorithm on Hapmap 3 data.
In both codes, there are a few hyper parameters which can be changed by user:
k: 		number of latent components
ratio:	\theta in the paper.
pcut: 	\tau in the paper.
binnum:	a specifc paramter for treeSGS_b method, which control the bin size of histogram but have limited impact.
sigma:	a specifc paramter for treeSGS_g method, which control the smoothness of Gaussian kernel but have limited impact.


We use SLEP to solve Lasso problem on U. Information about SLEP can be found at http://www.yelab.net/software/SLEP/