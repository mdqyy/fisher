// Matlab wrapper around GMM code of Jorge Sanchez
// Wrapper made by Jasper Uijlings
#include <cmath>
#include "mex.h"
#include "fisher.h"
#include "gmm.h"

void mexFunction(int nlhs, mxArray *out[], int nrhs, const mxArray *input[])
{
    // Checking number of arguments
    if (nlhs > 1){
        mexErrMsgTxt("Function has a single output: the fisher vector");
        return;
    }

    if (nrhs != 2){
        mexErrMsgTxt("Usage: fisherVector = mexFisherAssign(features, gmmModelFilename)");
        return;
    }
    
    if (mxIsChar(input[1]) != 1){
        mexErrMsgTxt("second input argument should be a string.");
        return;
    }

    // Load in arrays
    double* features = mxGetPr( input[0] );

    // The features are column wise. First value is the feature dimension
    int* featSize = (int*) mxGetDimensions(input[0]);
    int dimensionality = featSize[0];
    int nrFeatures = featSize[1];
    char* filename = mxArrayToString(input[1]);

    // Initialize gaussian mixture model
    gaussian_mixture<double>* gmmModel;
    gmmModel = new gaussian_mixture<double>(filename);
    
    // Check if features and gmm model correspond
    if (gmmModel->n_dim() != dimensionality){
        mexPrintf("Features: %d dimensions, gmm: %d dimensions", dimensionality, gmmModel->n_dim());
        mexErrMsgTxt("Dimensionality of features and dimensions inconsistent.");
        return;
    }   
    
    // Initialize fisher class
    fisher_param* param = new fisher_param;
    // !!Remove Power and L2 Normalization
    param->alpha = 1.0f;
    param->pnorm = 0.0;
    //param->grad_weights = true;
    //
    fisher<double>* fisherProjector;
    fisherProjector = new fisher<double>(*param);

    // Assign gmm model to fisher class
    fisherProjector->set_model(*gmmModel);

    // Now make vector arrays from input arguments (which are double arrays)
    std::vector<double*> featVecs(nrFeatures);
    for(int i=0; i < nrFeatures; i++){
        featVecs[i] = features + (i * dimensionality);
    }

    // Allocate memory for returning fisher vector
    int fisherDim = fisherProjector->dim();
    out[0] = mxCreateDoubleMatrix(fisherDim, 1, mxREAL);
    double* fisherVector = mxGetPr(out[0]);

    // Finally, compute the fisher vector
    fisherProjector->compute(featVecs, fisherVector);

    // Free memory
    delete fisherProjector;
    delete param;
    delete gmmModel;
    return;
}
 

