function [ fs, scores, ovlps ] = computeNMSFisherScore( fnet, fnorm, box, maskboxes, masks )
%COMPUTEAPPROXFISHERSCORE Summary of this function goes here
%   Detailed explanation goes here

nmodels = size(fnet, 2) ;
nmasks = size(masks, 1) ;
scores = zeros(nmasks, nmodels) ;

for i = 1 : nmasks
    scores(i, :) = computeFisherScore(fnet, fnorm, masks(i, :)) ;
end

fs = max(scores) ;

ovlps = boxOverlap(maskboxes', box') ;
%weighted mean
[ovlps, si] = sort(ovlps, 'descend') ;
scores = scores(si, :) ;

end

