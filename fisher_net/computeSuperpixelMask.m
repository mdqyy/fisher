function [ mask ] = computeSuperpixelMask( spmap, roi )
%COMPUTESUPERPIXELMASK Summary of this function goes here
%   Detailed explanation goes here

roimap = spmap(roi(2):roi(4), roi(1):roi(3));
idx = unique(roimap);

mask = zeros(1, max(max(spmap)));
mask(idx) = 1;

end

