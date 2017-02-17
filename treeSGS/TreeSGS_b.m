function [U, V, Cov, Obj, U_scale_his] = TreeSGS_b(X, gamma, lambda, k, C, maxiter, U_ini, V_ini, group, mx, Vnorm, ratio, pcut, binnum)
%% histogram version of treeSGS algorithm

%% initialization
[m, n] = size(X);
% for cost function
Obj = zeros(maxiter,4);
% for covergency, test by change of U and V
Cov = zeros(maxiter,2); 

% U's initial, pca or given
if ~exist('U_ini', 'var') || isempty(U_ini)
    [~, SCORE] = pca(X);
    U = SCORE(:,1:k);
else
    U = U_ini;
end;

% V's initial, all zeros or given
if  ~exist('V_ini', 'var') || isempty(V_ini)
    V = zeros(k, n);
else
    V = V_ini;
    norm_row = sum(V.^Vnorm, 2);
    norm_row(norm_row == 0) = 1;
    V = V ./repmat(norm_row.^(1/Vnorm), 1, n);
end;

U_scale_his = zeros(maxiter,1);

%%
for iter = 1:maxiter
    U_old = U;
    V_old = V;
    [Obj(iter,1), Obj(iter,2), Obj(iter,3), Obj(iter,4)]...
        = CostFunction(U, V, X, gamma, lambda, C);
    %%%%%%%%%%
    % solve V
    %%%%%%%%%%
    V = SolveV_SparseGroup(X, V, U, group, mx, Vnorm, ratio, pcut, binnum);
    
    disp(['Iter:' num2str(iter) ' V solved']);
    %%%%%%%%%%
    % solve U
    %%%%%%%%%%
    % Lasso or Fused Lasso, use SLEP package
    Beta = FusedLasso_SLEP(X', V', U_old', gamma, lambda);   
    U = Beta';
    disp(['Iter:' num2str(iter) ' U solved']);
    
    % failed
    if norm(U,'fro') == 0
        U = rand(m,k)*2-1;
        disp('Warning! U=0');
        pause;
        break;
    end;
    U_scale_his(iter) = max(max(abs(U)));
    %%%%%%%%%%
    % convergency test
    %%%%%%%%%%
    re1 = norm(U - U_old, 'fro')/norm(U_old, 'fro');
    re2 = norm(V - V_old, 'fro')/norm(V_old, 'fro');
    Cov(iter,:) = [re1 re2];
    disp([re1 re2]);
    
    if max([re1,re2]) < 5e-3
        break;
    end
end

Obj = Obj(1:iter,:);
Cov = Cov(1:iter,:);
if iter == maxiter
    disp(['algo didn''t converge in ' int2str(maxiter) ' iterations.']);
else
    disp(['algo converge in ' int2str(iter) ' iterations.']);
end;

end

%%
function V = SolveV_SparseGroup(X, V, U, group, mx, Vnorm, ratio, pcut, binnum)
iter_max = 100;
k = size(U, 2);
n_nodes = size(mx,1);
n_leaves = size(group,1);

[~, children, alldes, ~, range, allanc] = gettreeinfo(mx, group);

% each group's index and leftover index
ind = cell(size(group, 1), 1);
ind_other = 1 : size(X, 1);
for i = 1 : length(ind)
    ind{i} = find(group(i, :) == 1);
    ind_other = setdiff(ind_other, ind{i});
end

for iter = 1 : iter_max
    V_old = V;
    
    parts = zeros(k,1);
    for kk = 1 : k
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % X = U*V => V = U'X/(U'U), U,V are vectors
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tmp = U(:,kk)' * U(:,kk);
        if tmp ~= 0
            vv = U(:,kk)' * (X - U(:, [1:kk-1 kk+1:end]) * V([1:kk-1 kk+1:end], :)) / tmp;
            
            % use binnum to get histogram bins
            maxvv = max(vv);
            minvv = min(vv);
            edges = minvv:((maxvv-minvv)/(binnum-1)):maxvv;

            % calculate entropy
            etpy = zeros(n_nodes,1);
            etpy_gain = zeros(n_nodes,1);
            
            for i = 1:n_nodes
%                 ccc = histcounts(vv(range{i}), edges);
                ccc = histc(vv(range{i}), edges);
                ccc(ccc==0) = [];
                etpy(i) = -sum(ccc/sum(ccc).*log(ccc/sum(ccc)))/log(binnum);
            end;

            % calculate informatin gain
            for i = (n_leaves+1):n_nodes
                if ~isempty(children{i})
                    ch = children{i};
                    range1 = range{ch(1)};
                    range2 = range{ch(2)};

                    etpy_gain(i) = etpy(i) ...
                        - length(range1)/(length(range1) + length(range2))*etpy(ch(1))...
                        - length(range2)/(length(range1) + length(range2))*etpy(ch(2));
                end;     
            end;
            
            % decide split based on pcut
            split = zeros(n_nodes,1);
            % large infor gain -> split -> split on all ancestor
            split(etpy_gain > pcut) = 1;
            for i = n_nodes:-1:(n_leaves+1)
                if split(i) == 1
                    split(allanc{i}) = 1;
                end;     
            end;            
            
            merge = 1 - split;           
            for i = (n_leaves+1):n_nodes
                if merge(i) == 1
                    merge(alldes{i}) = 0;
                end;     
            end;            

            parts(kk) = sum(merge);
            % use merge nodes as seperation on data
            % each non-overlapping group, a value 
            comb_pt = find(merge);
            comb_value = zeros(size(comb_pt));
            for i = 1:length(comb_pt)
                sel_rge = range{comb_pt(i)};
                comb_value(i) = norm(vv(sel_rge)) / sqrt(length(sel_rge));
            end;
            
            [B, IX] = sort(comb_value, 'descend');
            num_selected = find(cumsum(B) / sum(B) >= ratio, 1);
            
            ind_selected = [];
            for l = 1 : num_selected
                ind_selected = union(ind_selected, range{comb_pt(IX(l))});
            end

            % solve V by sparse group selection
            V(kk, ind_selected) = max(vv(ind_selected), 0);
            if norm(V(kk, ind_selected)) > 0
                V(kk, ind_selected) = V(kk, ind_selected) / norm(V(kk, ind_selected), Vnorm);
            else
                [~, I] = max(vv(ind_selected));
                V(kk, ind_selected(I)) = 1;
            end
            V(kk, setdiff(1: size(V, 2), ind_selected)) = 0;
        else
            V(kk, :) = 0;
        end
    end
    
    if norm(V - V_old, 'fro') < 1e-5 * norm(V_old, 'fro');
        break;
    end
end

if iter == iter_max
    disp(['GroupLasso_Ratio didn''t converge in ' int2str(iter) ' iterations.']);
end
end


%%
function [TotalObj, RecErr, L1norm, FLnorm] = CostFunction(U, V, X, gamma, lambda, C)
RecErr = 1/2 * norm(X - U * V, 'fro')^2;
FLnorm = gamma * sum(sum(abs(C*U)));
L1norm = lambda * sum(sum(abs(U)));
TotalObj = RecErr + L1norm + FLnorm;
end




