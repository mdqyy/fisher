function [ scores ] = fisher_net_deepspm( filename, sift, frames, im_size, spm, w )
%FISHER_NET_DEEPSPM Summary of this function goes here
%   Detailed explanation goes here

if ~isa(sift, 'double')
    error('The features should be stored as doubles.') ;
end

height = im_size(1) ;
width = im_size(2) ;

divs_row = [] ;
divs_col = [] ;
for i = 1 : size(spm, 1)
    divs_row = unique([ceil(linspace(1, height+1, spm(i, 1)+1)), divs_row]) ;
    divs_col = unique([ceil(linspace(1, width+1, spm(i, 2)+1)), divs_col]) ;
end

[ids, nodes] = computeIrregularGridAssign(frames, height, width, divs_row, divs_col) ;

%c++ array starts from 0
ids = ids - 1 ;

%tic
fish = mexFisherNetAssign(sift, ids, numel(nodes), filename) ;
%toc

%tic
fnorm = fish' * fish ;
%toc

if ~isempty(w)
    
    fkdim = size(fish, 1) ;
    nmodels = size(w, 2) ;
    npyrs = sum(prod(spm, 2)) + 1 ;
    score = zeros(npyrs, nmodels) ;
    
    gbox = [1, 1, width, height] ;
    mask = computeIrregularGridMask(height, width, divs_row, divs_col, gbox) ;
    fs = sum(fish(:, mask==1)' * w(1:fkdim, :), 1) ;
    ll = sum(sum(fnorm(mask==1, mask==1))) ;
    score(1, :) = fs ./ sqrt(ll) ;
    
    ss = 2 ;
    for i = 1:size(spm, 1)
       pyr = spm(i, :) ;
       stepx = (gbox(3) - gbox(1) + 1) / pyr(2) ;
       stepy = (gbox(4) - gbox(2) + 1) / pyr(1) ;

       for p = 1:pyr(2)
           for q = 1:pyr(1)
               gbox_pyr(1) = ceil(gbox(1) + (p-1)*stepx) ;
               gbox_pyr(2) = ceil(gbox(2) + (q-1)*stepy) ;
               gbox_pyr(3) = ceil(gbox(1) + p*stepx - 1) ;
               gbox_pyr(4) = ceil(gbox(2) + q*stepy - 1) ;

               gbox_pyr(3) = min([gbox(3), gbox_pyr(3)]) ;
               gbox_pyr(4) = min([gbox(4), gbox_pyr(4)]) ;

               mask_pyr = computeIrregularGridMask(height, width, divs_row, divs_col, gbox_pyr) ;
               fs = sum(fish(:, mask_pyr==1)' * w((ss-1)*fkdim+1:ss*fkdim, :), 1) ;
               ll = sum(sum(fnorm(mask_pyr==1, mask_pyr==1))) ;
               score(ss, :) = fs ./ sqrt(ll) ;
               ss = ss + 1 ;
           end
       end

    end
    
    scores = sum(score, 1) ./ sqrt(npyrs) ;
end
