function [ lmap ] = fisher_net_graphcut( fnet, fnorm )
%Summary of this function goes here
%   Detailed explanation goes here

segclass = zeros(size(fnet, 2), 1);
unary = [fnet; -fnet];
pairwise = sparse(fnorm);
labelcost = [0 1; 1 0];

[ lmap ] = GCMex( segclass, single(unary), pairwise, single(labelcost), 0 );

end

