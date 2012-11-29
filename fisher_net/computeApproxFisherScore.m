function [ fs ] = computeApproxFisherScore( fnet, fnorm, box, maskboxes, masks )
%COMPUTEAPPROXFISHERSCORE Summary of this function goes here
%   Detailed explanation goes here

expmu = 0.013 ;% 0.013 ;
nmodels = size(fnet, 2) ;
nmasks = size(masks, 1) ;
scores = zeros(nmasks, nmodels) ;

for i = 1 : nmasks
    scores(i, :) = computeFisherScore(fnet, fnorm, masks(i, :)) ;
end

if nmasks == 1
    fs = scores ;
else
    
    olaps = boxOverlap(box', maskboxes') ;
    
    %weighted mean
    [olaps, si] = sort(olaps, 'descend') ;
    scores = scores(si, :) ;
    weights = exppdf((1-olaps), expmu) ;
    weights = repmat(weights', 1, nmodels) ;
    fs = wmean(scores, weights, 1) ;

    %linear fitting
%     pw_olaps = boxOverlap(maskboxes', maskboxes') ;
%     olap = boxOverlap(box', maskboxes') ;
%     fs = zeros(1, nmodels) ;
%     for j = 1 : nmodels
%         b = scores(:, j) ;
%         A = pw_olaps ;
%         [X, ~] = regress(b, A) ;
%         fs(1, j) = olap * X ;
%         %[X, ~] = regress(b, [A, ones(nmasks, 1)]) ;
%         %fs(1, j) = [olap, 1] * X ;
%     end
    
end

end

