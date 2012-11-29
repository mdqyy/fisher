function fish = normalizeFisherColsL2(fish, pyramids)

if nargin == 1

    %sgn = sign(fish) ;
    %fish = sgn .* sqrt(abs(fish)) ;
    fish = normalizeColsL2(fish) ;
    
else
    
    %sgn = sign(fish) ;
    %fish = sgn .* sqrt(abs(fish)) ;
    
    npyrs = sum(prod(pyramids, 2)) + 1 ;
    perpyr = size(fish, 1) / npyrs ;
    
    for i = 1 : npyrs
        
        %fish((i - 1) * perpyr + 1 : i * perpyr) = normalizeColsL2(fish((i - 1) * perpyr + 1 : i * perpyr)) ;
        fish((i - 1) * perpyr + 1 : i * perpyr, :) = normalizeColsL2(fish((i - 1) * perpyr + 1 : i * perpyr, :)) ;
        
    end
    
    fish = normalizeColsL2(fish) ;
    
end