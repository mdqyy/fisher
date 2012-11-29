function [ ids, nodes ] = computeIrregularGridAssign(frames, height, width, rdivs, cdivs)
% 

if rdivs(1) ~= 1
    rdivs = [1, rdivs] ;
end

if cdivs(1) ~= 1
    cdivs = [1, cdivs] ;
end

if rdivs(end) ~= height+1
    rdivs = [rdivs, height+1] ;
end

if cdivs(end) ~= width+1
    cdivs = [cdivs, width+1] ;
end

rnd = numel(rdivs) - 1 ;
cnd = numel(cdivs) - 1 ;

binsx = vl_binsearch(cdivs, frames(1, :)) ;
binsy = vl_binsearch(rdivs, frames(2, :)) ;

ids = int32((binsx - 1) .* rnd + binsy) ;

nodes = zeros(rnd*cnd, 1) ;
nodes = vl_binsum(nodes, ones(size(ids)), ids) ;
nodes = reshape(nodes, rnd, cnd) ;

end