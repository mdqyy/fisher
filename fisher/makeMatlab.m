% Compile Matlab wrapper around GMM code of Jorge Sanchez
%
% Learn the gmm model using ClusteredFeatures2Gmm. See this file for more
% details about the usage

% Creating mixture model 
%mex -g mexGmm.cpp gmm.cxx stat.cxx simd_math.cxx
mex mexGmm.cpp gmm.cxx stat.cxx simd_math.cxx

% For assigning fisher vectors
%mex -g mexFisherAssign.cpp fisher.cxx gmm.cxx stat.cxx simd_math.cxx
mex mexFisherAssign.cpp fisher.cxx gmm.cxx stat.cxx simd_math.cxx