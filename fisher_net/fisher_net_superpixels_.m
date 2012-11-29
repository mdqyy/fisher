function [ fnet, fnorm, spix ] = fisher_net_superpixels( img, spix, gmmModel, w )
%FISHER_NET_SUPERPIXELS Summary of this function goes here
%   Detailed explanation goes here

%generate superpixels
if strcmp(spix.superpixelmethod, 'turbopixels')
    [phi, boundary, spmap] = superpixels(im2double(img), spix.numSuperpixels);
    spix.spmap = spmap ;
    spix.phi = phi ;
    spix.boundary = boundary ;
    spix.nsppixels = numel(unique(spix.spmap)) ;
elseif strcmp(spix.superpixelmethod, 'arbelaez')
    [spmap, boundary, ucm2] = arbelaezMalikRegions(spix.impath, spix.arbelaezpath, spix.malikthres) ;
    spix.spmap = spmap ;
    spix.boundary = boundary ;
    spix.ucm2 = ucm2 ;
    spix.nsppixels = numel(unique(spix.spmap)) ;
else
    errormessage = sprintf('There are supported two types of superpixels: "turbopixels" [1] or "arbelaez" [2]\n') ;
    errormessage = [errormessage, sprintf('[1] http://www.cs.toronto.edu/~babalex/research.html \n')] ;
    errormessage = [errormessage, sprintf('[2] http://www.eecs.berkeley.edu/Research/Projects/CS/vision/grouping/resources.html\n')] ;
    error(errormessage) ;
end

%generate fisher net
img_gray = rgb2gray(img);

[frames, sift] = vl_phow(single(img_gray), 'Step', gmmModel.ss);
nzix = sum(sift, 1) ~= 0;     % OPTIONAL
sift = double(sift(:, nzix)); % OPTIONAL, although if there are too many 0 vectors, performance decreases significantly
sift = normalizeColsL2(sift); % OPTIONAL
sift = gmmModel.pcamap' * sift;
frames = frames(:, nzix);

[ids, nsppixels, spix.hist] = computeSuperpixelAssign(frames, spmap);
ids = ids - 1;

fish = mexFisherNetAssign(sift, ids, nsppixels, gmmModel.gmmfilename);
fnorm = fish' * fish;
    
if nargin > 3
    fnet = w' * fish;
elseif nargin > 2
    fnet = fish;
end

