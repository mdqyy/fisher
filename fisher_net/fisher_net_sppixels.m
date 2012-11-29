function [ fnet, fnorm, spix ] = fisher_net_sppixels( img, spix, gmmModel, w )
%FISHER_NET_SUPERPIXELS Summary of this function goes here
%   Detailed explanation goes here

%generate superpixels
%spix struct generated outside

%generate fisher net
img_gray = rgb2gray(img);

[frames, sift] = vl_phow(single(img_gray), 'Step', gmmModel.ss);
nzix = sum(sift, 1) ~= 0;     % OPTIONAL
sift = double(sift(:, nzix)); % OPTIONAL, although if there are too many 0 vectors, performance decreases significantly
sift = normalizeColsL2(sift); % OPTIONAL
sift = gmmModel.pcamap' * sift;
frames = frames(:, nzix);

[ids, nsppixels, spix.hist] = computeSuperpixelAssign(frames, spix.spmap);
ids = ids - 1;

fish = mexFisherNetAssign(sift, ids, nsppixels, gmmModel.gmmfilename);
fnorm = fish' * fish;
    
if nargin > 3
    fnet = w' * fish;
elseif nargin > 2
    fnet = fish;
end

