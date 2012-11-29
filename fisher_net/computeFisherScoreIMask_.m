function [ fs ] = computeFisherScoreIMask_( fnet, fnorm, imask1, imask2 )
%
%
if nargin < 4
    imask2 = imask1;
end

fs = imask1 * fnet ;

imask2_norm = imask2' * imask2 ;
ll = sum(sum(fnorm .* imask2_norm)) ;

if ll < 1e-5
    fs = zeros(1, size(fs, 2)) ;
    return ;
end

fs = fs ./ sqrt(ll) ;
        
end


