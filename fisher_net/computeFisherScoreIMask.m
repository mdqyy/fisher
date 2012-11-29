function [ fs ] = computeFisherScoreIMask( fnet, fnorm, imask1, imask2 )
%
%
if nargin < 4
    imask2 = imask1;
end

fs = imask1 * fnet ;

mask2 = imask2 > 0 ;
immask2 = imask2(mask2==1) ;
mask2_norm = immask2' * immask2 ;
ll = sum(sum(fnorm(mask2==1, mask2==1) .* mask2_norm)) ;

if ll < 1e-5
    fs = zeros(1, size(fs, 2)) ;
    return ;
end

fs = fs ./ sqrt(ll) ;
        
end
