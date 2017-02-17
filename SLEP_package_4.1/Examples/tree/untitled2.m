clear, clc;



cd ..
cd ..

root=cd;
addpath(genpath([root '/SLEP']));
                     % add the functions in the folder SLEP to the path
                   
% change to the original folder
cd Examples/tree;

n=10000;
v=randn(n,1);

lambda = general_findLambdaMax(v, n, G, ind, nodes);
a = 0.999;
% Group Property
ind=[[1, 2000, a*lambda*sqrt(2000)]', [2001, 4000, a*lambda*sqrt(2000)]',... % the layer above the leaf
    [4001, 5000, a*lambda*sqrt(1000)]', [5001, 7000, a*lambda*sqrt(2000)]', [7001,10000, a*lambda*sqrt(3000)]',...
    [1, 5000, a*lambda*sqrt(5000)]', [5001, 10000, a*lambda*sqrt(5000)]']; % the higher layer



G=[1:10000];
nodes=size(ind,2);

x=general_altra(v, n, G, ind, nodes);

figure;
plot(v);
hold on;
plot(x, 'r');
title(['lambda = ' num2str(a) '*lambda\_max']);


