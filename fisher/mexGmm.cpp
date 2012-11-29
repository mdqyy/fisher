// Matlab wrapper around GMM code of Jorge Sanchez
// Wrapper made by Jasper Uijlings
#include <cmath>
#include "mex.h"
#include "gmm.h"

void mexFunction(int nlhs, mxArray *out[], int nrhs, const mxArray *input[])
{
    // Checking input/output arguments
    if (nlhs != 0){
        mexErrMsgTxt("Error: function has one output parameter");
        return;
    }

    if (nrhs != 5){
        mexErrMsgTxt("Error: mexGmm needs five input parameters. Features, centers, variances, weights, and filename.");
        return;
    }

    if (mxIsChar(input[4]) != 1){
        mexErrMsgTxt("Fifth input argument should be a string.");
        return;
    }

    // Load in arrays
    double* features = mxGetPr( input[0] );
    double* centres = mxGetPr(input[1]);
    double* variances = mxGetPr(input[2]);
    double* priors = mxGetPr(input[3]);
    char* filename = mxArrayToString(input[4]);

    // The features are column wise. First value is the feature dimension
    int* featSize = (int*) mxGetDimensions(input[0]);
    int dimensionality = featSize[0];
    int nrFeatures = featSize[1];
    int* centreSize = (int*) mxGetDimensions(input[1]);
    int nrGauss = centreSize[1]; // # of centres for gaussian mixture
    
    // Check if centers and features have the same dimensionality
    if (featSize[0] != centreSize[0])
        mexErrMsgTxt("Cluster centers and features have a different dimensionality");

    // Initialize gaussian mixture model
    gaussian_mixture<double>* gmmModel;
    gmmModel = new gaussian_mixture<double>(nrGauss, dimensionality);
    
    // Test if model can be saved before doing any calculations
    if (gmmModel->save(filename) == -1){
        mexErrMsgTxt("Incorrect filename: cannot save.");
    }

    // Now make vector arrays from input arguments (which are double arrays)
    // Note that no copying of the values takes place
    //std::vector<double*> featVecs(nrFeatures);
    //std::vector<double*> centreVecs(nrFeatures);
    //std::vector<double*> varVecs(nrFeatures);
    //for(int i=0; i < nrFeatures; i++){
    //    featVecs[i] = features + (i * dimensionality);
    //    centreVecs[i] = centres + (i * dimensionality);
    //    varVecs[i] = variances + (i * dimensionality);
    //}
	std::vector<double*> featVecs(nrFeatures);
    for(int i=0; i < nrFeatures; i++){
        featVecs[i] = features + (i * dimensionality);
    }
	std::vector<double*> centreVecs(nrGauss);
    std::vector<double*> varVecs(nrGauss);
    for(int i=0; i < nrGauss; i++){
        centreVecs[i] = centres + (i * dimensionality);
        varVecs[i] = variances + (i * dimensionality);
    }

    std::vector<double> priorVec(priors, priors + nrGauss);

    // Set the initial parameters for the gmm model
    gmmModel->set(centreVecs, varVecs, priorVec);
    
    // Now train the gmm model 
    gmmModel->em(featVecs);

    gmmModel->save(filename);
    // gmmModel->print(true, true, true);

    // Free memory
    delete gmmModel;
    //delete filename;
    return;
}
 
