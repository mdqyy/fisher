function fish = fisher(filename, sift, frames, height, width, spm)

if ~isa(sift, 'double')
    error('The features should be stored as doubles.') ;
end

if nargin > 5
    ixs = computeSpatialConfigurations(frames, height, width, spm) ;
    fish = cell(size(ixs, 1) + 1, 1) ;
    fish{1} = mexFisherAssign(sift, filename) ;
    
    for i = 1 : size(ixs, 1)
        fish{i + 1} = mexFisherAssign(sift(:, ixs(i, :) == 1), filename) ;
    end
    fish = cat(1, fish{:}) ;
elseif nargin > 1
    fish = mexFisherAssign(sift, filename) ;
else
    fish = fishervec(wgt, mu, sigma, sift) ;
end

