function [ fnet, fnorm, spix, fish, nfeatspersp ] = fisher_net_superpixels( img, spix, gmmModel, varargin)
%FISHER_NET_SUPERPIXELS Summary of this function goes here
%   Detailed explanation goes here
fnet = [] ;
fnorm = [] ;
nonorm = false ;
w = [] ;

for q = 1 : numel(varargin)
    if strcmp(varargin{q}, 'model')
        w = varargin{q + 1} ;
    elseif strcmp(varargin{q}, 'nonorm')
        nonorm = true ;
    end
end

%generate superpixels
if isfield(spix, 'loadfrom')
    ldd = load(spix.loadfrom) ;
    spix.spmap = ldd.spix.spmap ;
    spix.phi = ldd.spix.phi ;
    spix.boundary = ldd.spix.boundary ;
    spix.nsppixels = ldd.spix.nsppixels ;
else
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
    elseif strcmp(spix.superpixelmethod, 'CPMC')
        % Nothing to be done, alreday loaded
%         1+1 ;
%         ids = spix.mask_data.sp ;
    else
        errormessage = sprintf('There are supported two types of superpixels: "turbopixels" [1] or "arbelaez" [2]\n') ;
        errormessage = [errormessage, sprintf('[1] http://www.cs.toronto.edu/~babalex/research.html \n')] ;
        errormessage = [errormessage, sprintf('[2] http://www.eecs.berkeley.edu/Research/Projects/CS/vision/grouping/resources.html\n')] ;
        error(errormessage) ;
    end
end

%generate fisher net
img_gray = rgb2gray(img);

[frames, sift] = vl_phow(single(img_gray), 'Step', gmmModel.ss);
nzix = sum(sift, 1) ~= 0;     % OPTIONAL
sift = double(sift(:, nzix)); % OPTIONAL, although if there are too many 0 vectors, performance decreases significantly
sift = normalizeColsL2(sift); % OPTIONAL
sift = gmmModel.pcamap' * sift;
frames = frames(:, nzix);

[ids, nsppixels, nfeatspersp] = computeSuperpixelAssign(frames, spix.spmap);
ids = ids - 1;

fish = mexFisherNetAssign(sift, ids, nsppixels, gmmModel.gmmfilename);

if ~nonorm
    fnorm = fish' * fish;
end
    
if ~isempty(w)
    fnet = w' * fish;
end


