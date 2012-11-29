function computeGmmForFisher(filename, features, centres, clusterIdx, maxiter)
% ClusteredFeatures2Gmm(filename, features, centres, clusterIdx, variances, weights)
% Creates GMM model from the features. Features ans centres are column based.
%
% filename:             Filename where gmm model will be saved
% features:             D x N matrix of N column based features of size D.
% centres:              D x K matrix with K centres with dimensionality D.
% clusterIdx:           N x 1 vector denoting assignment of features to centers.
% variances (optional): Variances of the clusters
% weights (optional):   Weights of the clusters
%
% GMM model needs kmeans initialisation for computational efficiency. When
% using the matlab native kmeans, centres and clusterIdx correspond to C
% and IDX respectively (see help of kmeans). Note that this kmeans uses
% row-based features while this function needs column-based features.
% 
% After training the gmm, fisher vectors are created as follows:
%   fisherVector = mexFisherAssign(newFeatures, filename)
%
% Matlab wrapper written by Jasper Uijlings
% Original GMM code by Jorge Sanchez

if ~isa(features, 'double')
    error('The features should be stored as doubles.') ;
end

if size(features,1) ~= size(centres,1)
    error('Features and cluster centres should have the same dimensionality in col direction');
end

if size(features,2) ~= length(clusterIdx)
    error('There should be as many cluster assignments as features');
end

% Calculate variances per cluster
variances = zeros(size(centres));
for i=1:max(clusterIdx)
    variances(:,i) = var(features(:,clusterIdx == i), 0, 2);
end

% Calculate weights of clusters (i.e. priors)
numCentres = max(clusterIdx);
weights = histc(clusterIdx, 1:numCentres) ./ size(features,2);

if nargin < 5
    maxiter = 100 ;
end

%mexGmm(features, centres, variances, weights, filename, maxiter);
mexGmm(features, centres, variances, weights, filename);
