function [parent, children, alldes, sibling, range, allanc] = gettreeinfo(mx, group)
n_node = size(mx,1);

% parent
parent = zeros(n_node,1);
for i = 1:n_node
    IX = find(mx(i,:));
    if ~isempty(IX)
        parent(i) = IX;
    end;
end;

% children
children = cell(n_node,1);
for i = 1:n_node
    IX = find(mx(:,i));
    if ~isempty(IX)
        children{i} = IX;
    end;
end;

% sibling
sibling = zeros(n_node,1);
for i = 1:n_node
    p = parent(i);
    if p ~= 0
        c = children{p};
        if ~isempty(c)
            sibling(i) = setdiff(c,i);
        end;
    end;
end;

% all descent
alldes = cell(n_node, 1);
for i = 1:n_node
    if i <= (n_node+1)/2
        alldes{i} = i;
    end;
    if parent(i) ~= 0
        alldes{parent(i)} = union(alldes{parent(i)}, i);
        alldes{parent(i)} = union(alldes{parent(i)}, alldes{i});
    end;
end;

% range
range = cell(n_node,1);
for i = 1:((n_node+1)/2)
    range{i} = find(group(i,:));    
end;
for i = ((n_node+1)/2)+1:n_node
    range{i} = union(range{children{i}(1)}, range{children{i}(2)});
end;

% all ancestor
allanc = cell(n_node, 1);
for i = 1:n_node
    if ~isempty(alldes{i})
        tt = alldes{i};
        for j = 1:length(tt)
            allanc{tt(j)} = [allanc{tt(j)} i]; 
        end;
    end;
end;

end