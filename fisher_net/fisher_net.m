function [ fish fnet fnorm nodes ] = fisher_net( filename, sift, frames, im_size, grid_configs, w, spm )

if ~isa(sift, 'double')
    error('The features should be stored as doubles.') ;
end

if numel(grid_configs) == 1
    cell_size = grid_configs ;
    [ids nodes] = computeGridAssign(frames, im_size(1), im_size(2), [cell_size, cell_size]) ;
elseif numel(grid_configs) == 2 && ~iscell(grid_configs)
    cell_height = grid_configs(1) ;
    cell_width = grid_configs(2) ;
    [ids nodes] = computeGridAssign(frames, im_size(1), im_size(2), [cell_height, cell_width]) ;
elseif iscell(grid_configs)
    divs_row = grid_configs{1} ;
    divs_col = grid_configs{2} ;
    [ids, nodes] = computeIrregularGridAssign(frames, im_size(1), im_size(2), divs_row, divs_col) ;
end

%c++ array starts from 0
ids = ids - 1 ;

%tic
fish = mexFisherNetAssign(sift, ids, numel(nodes), filename) ;
%toc

%tic
fnorm = fish' * fish ;
%toc

fnet = [] ;
if nargin == 6 && ~isempty(w)
    fnet = fish' * w  ;
elseif nargin == 7
    % ... needs to be done
    fkdim = size(fish, 1) ;
    nnodes = numel(nodes) ;
    nmodels = size(w, 2) ;
    npyrs = sum(prod(spm, 2)) + 1 ;
    fnet = zeros(nnodes, npyrs, nmodels) ;
    
    for i = 1:nmodels
        model = w(:, i) ;
        model = reshape(model, fkdim, npyrs) ;
        fnet(:, :, i) = fish' * model ;
    end
end

end

%square-root on each cell
%sgn = sign(fish) ;
%fish = sgn .* sqrt(abs(fish)) ;
