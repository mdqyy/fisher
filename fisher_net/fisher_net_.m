function [fnet, fnorm, nodes] = fisher_net_(filename, sift, frames, im_size, node_size, w, confs)

if ~isa(sift, 'double')
    error('The features should be stored as doubles.') ;
end

if nargin > 6
    
    fprintf('spatial pyramids not support now\n') ;
    %with spatial pyramids
    return ;
    
elseif nargin > 5
    
    [ids nodes] = computeNodeAssign(frames, im_size(1), im_size(2), node_size) ;
    ids = ids - 1 ;
    fish = mexFisherNetAssign(sift, ids, numel(nodes), filename) ;
    
    fnorm = fish' * fish ;
    fnet = reshape(w' * fish, size(nodes, 1), size(nodes, 2)) ;
    %fnet = squeeze(reshape(w' * fish, size(w, 2), size(nodes, 1), size(nodes, 2))) ;
    
elseif nargin > 4
    
    [ids nodes] = computeNodeAssign(frames, im_size(1), im_size(2), node_size) ;
    ids = ids - 1 ;
    fish = mexFisherNetAssign(sift, ids, numel(nodes), filename) ;
    
    %square-root on each cell
    %sgn = sign(fish) ;
    %fish = sgn .* sqrt(abs(fish)) ;
    
    %tt = tic ;
    fnorm = fish' * fish ;
    %fprintf('time cost: %f\n', toc(tt)) ;
    
    fnet = fish ;
    
end
