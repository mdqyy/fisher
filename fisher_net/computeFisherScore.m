function [ fs ] = computeFisherScore( fnet, fnorm, mask )

fs = sum(fnet(mask==1, :), 1) ;
ll = sum(sum(fnorm(mask==1, mask==1))) ;

if ll < 1e-5
    fs = zeros(1, size(fs, 2)) ;
    return ;
end

fs = fs ./ sqrt(ll) ;

end