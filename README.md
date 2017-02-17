run_TBL_b.m	the main script to run treeSGS_b algorithm on Hapmap 3 data.
run_TBL_g.m	the main script to run treeSGS_g algorithm on Hapmap 3 data.
In both codes, there are a few hyper parameters which can be changed by user:
k: 		number of latent components
ratio:	\theta in the paper.
pcut: 	\tau in the paper.
binnum:	a specifc paramter for treeSGS_b method, which control the bin size of histogram but have limited impact.
sigma:	a specifc paramter for treeSGS_g method, which control the smoothness of Gaussian kernel but have limited impact.


We use SLEP to solve Lasso problem on U. Information about SLEP can be found at http://www.yelab.net/software/SLEP/