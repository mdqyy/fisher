function [ mask, lc ] = computeIrregularGridMask(height, width, rdivs, cdivs, roi)
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

binsx = vl_binsearch(cdivs, roi([1 3])) ;
binsy = vl_binsearch(rdivs, roi([2 4])) ;

mask = zeros(rnd, cnd) ;
mask(binsy(1):binsy(2), binsx(1):binsx(2)) = 1 ;
mask = reshape(mask, 1, numel(mask)) ;

lc = [binsx(1), binsy(1), binsx(2), binsy(2)] ;

end