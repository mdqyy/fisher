function [ ff ] = aggreg_fisher( fish, fish_size, agg_factors )
%AGGREG_FISHER Summary of this function goes here
%   Detailed explanation goes here

rows = fish_size(1) ;
cols = fish_size(2) ;

fdim = size(fish,1) ;
fish = reshape(fish, fdim, rows, cols) ;

agg_rows = ceil(rows / agg_factors(1)) ;
agg_cols = ceil(cols / agg_factors(2)) ;

binsrow = vl_binsearch(1:agg_factors(1):rows+1, 1:rows) ;
binscol = vl_binsearch(1:agg_factors(2):cols+1, 1:cols) ;

ff = zeros(size(fish,1), agg_rows, agg_cols) ;

for i = 1:agg_rows
    for j = 1:agg_cols
        fish_bin = fish(:, binsrow==i, binscol==j) ; 
        fish_bin = reshape(fish_bin, fdim, size(fish_bin, 2)*size(fish_bin, 3)) ;
        ff(:, i, j) = sum(fish_bin, 2) ;
    end
end

ff = reshape(ff, fdim, agg_rows*agg_cols) ;

end

