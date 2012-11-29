% Compile Matlab wrapper around GMM code of Jorge Sanchez
%
% Learn the gmm model using ClusteredFeatures2Gmm. See this file for more
% details about the usage

% Creating mixture model 
%mex mexGmm.cpp gmm.cxx stat.cxx simd_math.cxx

% For assigning fisher vectors
mex mexFisherNetAssign.cpp fisher_net.cxx gmm.cxx stat.cxx simd_math.cxx
% mex mexFisherApproxDetector.cpp fisher_approx_detector.cxx fisher_net.cxx gmm.cxx stat.cxx simd_math.cxx
% mex mexFisherApproxDetectorEq.cpp fisher_approx_detector.cxx fisher_net.cxx gmm.cxx stat.cxx simd_math.cxx