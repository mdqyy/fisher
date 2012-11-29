// Wrapper made by Zhenyang Li
#include <cmath>
#include "mex.h"
#include "fisher_net.h"
#include "gmm.h"

void mexFunction(int nlhs, mxArray *out[], int nrhs, const mxArray *input[])
{
    // Checking number of arguments
    if (nlhs > 1){
        mexErrMsgTxt("Function has a single output: the fisher net");
        return;
    }

    if (nrhs != 4){
        mexErrMsgTxt("Usage: fisherVector = mexFisherNetAssign(features, ids, nnd, gmmModelFilename)");
        return;
    }
    
    if (mxIsChar(input[3]) != 1){
        mexErrMsgTxt("second input argument should be a string.");
        return;
    }

    // Load in arrays
    double* features = mxGetPr(input[0]);
    int* nodes = (int*) mxGetData(input[1]);
            
    // The features are column wise. First value is the feature dimension
    double* pnrNodes = mxGetPr(input[2]);
	int nrNodes = (int) (*pnrNodes);
    int* featSize = (int*) mxGetDimensions(input[0]);
    int dimensionality = featSize[0];
    int nrFeatures = featSize[1];
    char* filename = mxArrayToString(input[3]);

    // Initialize gaussian mixture model
    gaussian_mixture<double>* gmmModel;
    gmmModel = new gaussian_mixture<double>(filename);
    
    // Check if features and gmm model correspond
    if (gmmModel->n_dim() != dimensionality){
        mexPrintf("Features: %d dimensions, gmm: %d dimensions", dimensionality, gmmModel->n_dim());
        mexErrMsgTxt("Dimensionality of features and dimensions inconsistent.");
        return;
    }   
    
    // Initialize fisher_net class
    fisher_net_param* param = new fisher_net_param;
    // !!Remove Power and L2 Normalization
    param->alpha = 1.0f;
    param->pnorm = 0.0;
    //param->grad_weights = true;
    //
    fisher_net<double>* fisherProjector;
    fisherProjector = new fisher_net<double>(*param);

    // Assign gmm model to fisher_net class
    fisherProjector->set_model(*gmmModel);

    // Now make vector arrays from input arguments (which are double arrays)
    std::vector<double*> featVecs(nrFeatures);
    for(int i=0; i < nrFeatures; i++){
        featVecs[i] = features + (i * dimensionality);
    }

    // Allocate memory for returning fisher vector
    int fisherDim = fisherProjector->dim();
    //out[0] = mxCreateDoubleMatrix(fisherDim, 1, mxREAL);
    //double* fisherVector = mxGetPr(out[0]);
    
	out[0] = mxCreateDoubleMatrix(fisherDim, nrNodes, mxREAL);
    double* fisherVectorPr = mxGetPr(out[0]);
    std::vector<double*> fisherVector(nrNodes);
    for(int i=0; i<nrNodes; i++){
        fisherVector[i] = fisherVectorPr + (i * fisherDim);
    }

	//out[0] = mxCreateNumericMatrix(fisherDim, nrFeatures, mxSINGLE_CLASS, mxREAL);
    //float* fisherVectorPr = (float*) mxGetData(out[0]);
    //std::vector<float*> fisherVector(nrFeatures);
    //for(int i=0; i<nrFeatures; i++){
    //    fisherVector[i] = fisherVectorPr + (i * fisherDim);
    //}

    // Finally, compute the fisher vector
    fisherProjector->compute_fisher_net(featVecs, nodes, fisherVector);

    // Free memory
    delete fisherProjector;
    delete param;
    delete gmmModel;
    return;
}
 
