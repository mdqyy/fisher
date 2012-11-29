function [ ids, nodes ] = computeGridAssign(frames, height, width, cell_size)

cell_height = cell_size(1) ;
cell_width = cell_size(2) ;

rnd = ceil(height / cell_height) ;
cnd = ceil(width / cell_width) ;

binsx = vl_binsearch(1:cell_width:width+1, frames(1, :)) ;
binsy = vl_binsearch(1:cell_height:height+1, frames(2, :)) ;

ids = int32((binsx - 1) .* rnd + binsy) ;

nodes = zeros(rnd*cnd, 1) ;
nodes = vl_binsum(nodes, ones(size(ids)), ids) ;
nodes = reshape(nodes, rnd, cnd) ;

end