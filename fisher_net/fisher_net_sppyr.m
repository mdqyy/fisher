function [fvec ff] = fisher_net_sppyr(filename, sift, frames, im_size, confs)

if ~isa(sift, 'double')
    error('The features should be stored as doubles.') ;
end

height = im_size(1) ;
width = im_size(2) ;
% confs = [ones(1, size(confs, 2)); confs] ;

div_r = lcm(1, confs(1, 1)) ;
div_c = lcm(1, confs(1, 2)) ;
for x = 2 : size(confs, 1)
    div_r = lcm(div_r, confs(x, 1)) ;
    div_c = lcm(div_c, confs(x, 2)) ;
end

ixs = computeSpatialConfigurations(frames, height, width, confs) ;

ids = computeSpatialConfigurations(frames, height, width, [div_r, div_c]) ;
[ids, dummy] = find(ids) ;
ids = int32(ids) ;

% [ids nodes] = computeGridAssign(frames, im_size(1), im_size(2), [10, 10]) ;

fnet = mexFisherNetAssign(sift, ids, div_r * div_c, filename) ;

fvecall = sum(fnet, 2) ./ size(sift, 2) ;
fvecall = normalizeColsL2(fvecall) ;

fvec = zeros(size(fnet, 1), size(ixs, 1)) ; % Normalize it w.r.t the number of features T

for sp = 1 : size(ixs, 1)
   fi = find(ixs(sp, :)) ;
   spids = ids(fi) ;
   uqspids = unique(spids) ;
   
   fv = sum(fnet(:, uqspids), 2) ./ size(spids, 2) ; % Normalize it w.r.t the number of features T
   fvec(:, sp) = normalizeColsL2(fv) ;
end

ff = [fvecall(:), fvec] ;
fvec = normalizeColsL2(ff(:)) ;

%square-root on each cell
%sgn = sign(fish) ;
%fish = sgn .* sqrt(abs(fish)) ;
