
function Adj = get_adj(periods)
% column with same periods means a group, thus adjacency 

u  = unique(periods);
num_fre  = histc(periods,u);

num_u = length(u);

Ac =  cell(num_u);
for i =1:num_u
    tmp   = num_fre(i);
    Ac{i}  =  ones(tmp);
end

Adj = blkdiag(Ac{:});
