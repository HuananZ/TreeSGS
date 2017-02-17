clc, clear;
close all;
format compact;

% add SLEP package for solve Lasso problem on U
addpath(genpath('../SLEP_package_4.1/SLEP'));
%% load and process data
%%%%%%%%%%%%%%%%%%%%%%%%%%
% data: CNV feature by samples matrix, 0 means no CNV change (2 copies), -1
% means 1 CNV del (1 copy), 2 means 2 CNV amp (4 copies)
% cn_pos:   location of each CNV
% ThePop: population infor mation. 3rd column: population; 4th column super population
%%%%%%%%%%%%%%%%%%%%%%%%%%
load Hapmap_CNV_processed2;

X = data;
samples = ThePop;
[m, n] = size(X);

%% phylogenetic info
% pop seqence in data
poplist = {'ASW', 'LWK', 'MKK', 'YRI', 'MXL', 'CHB', 'CHD', 'JPT', 'CEU', 'TSI', 'GIH'};
pop_num = length(poplist);

pop_ix = cell(pop_num,1);
for i = 1:pop_num
    pop_ix{i} = find(strcmp(ThePop(:,3), poplist{i}));
end;

% group
group = zeros(pop_num, n);
for i = 1:pop_num
    group(i, pop_ix{i}) = 1;
end;

% tree relations
mx = zeros(pop_num*2-1);
mx([2 3], 12) = 1;
mx([4 12], 13) = 1;
mx([1 13], 14) = 1;
mx([6 7], 15) = 1;
mx([8 15], 16) = 1;
mx([9 10], 17) = 1;
mx([5 17], 18) = 1;
mx([11 18], 19) = 1;
mx([16 19], 20) = 1;
mx([14 20], 21) = 1;

%% treeSGS
% hyper-parameters: k, ratio, pcut, sigma
for k = 60
for ratio = 0.8
for pcut = [0.0001 0.0005 0.001 0.005]
for sigma = [20 40]
disp([k ratio pcut sigma]);
% Lasso parameter
lambda = 1;
% Fused Lasso parameter, not active in this version
gamma = 0;
Vnorm = 2;

maxiter = 100;
[m, n] = size(X);

% for further use on weighted Fused Lasso, not active in this version
C = zeros(m-1, m);

[U_sps, V_sps, Cov_sps, Obj_sps, U_scale_his_sps] = ...
    TreeSGS_g(X, gamma, lambda, k, C, maxiter, [], [], group, mx, Vnorm, ratio, pcut, sigma);

save(['TreeSGS_b_k' num2str(k) '_ratio' num2str(ratio) ...
    '_pcut' num2str(pcut) '_sigma' num2str(sigma) ...
    '_lambda' num2str(lambda) '_gamma' num2str(gamma) ...
    '_Vnorm' num2str(Vnorm) '.mat'],...
    'U_sps', 'V_sps', 'Cov_sps', 'Obj_sps', 'U_scale_his_sps');
end;
end;
end;
end;