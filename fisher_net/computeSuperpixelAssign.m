function [ ids, nsppixels, sphist ] = computeSuperpixelAssign( frames, spmap )
%COMPUTESUPERPIXELASSIGN Summary of this function goes here
%   Detailed explanation goes here

remap = reshape(spmap, 1, numel(spmap)) ;
nsppixels = numel(unique(spmap)) ;
ids = int32(remap(((frames(1, :)-1) .* size(spmap, 1)) + frames(2, :))) ;

sphist = zeros(nsppixels, 1) ;
sphist = vl_binsum(sphist, ones(size(ids)), ids) ;

end
