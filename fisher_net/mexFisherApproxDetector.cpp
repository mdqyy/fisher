// Wrapper made by Zhenyang Li
#include <cmath>
#include "mex.h"
#include "fisher_approx_detector.h"

void mexFunction(int nlhs, mxArray *out[], int nrhs, const mxArray *input[])
{
    // Checking number of arguments
    if (nlhs > 1){
        mexErrMsgTxt("Function has a single output");
        return;
    }

    if (nrhs != 9){
        mexErrMsgTxt("Usage: ");
        return;
    }
    
    if (mxIsChar(input[7]) != 1){
        mexErrMsgTxt("last input argument should be a string.");
        return;
    }

    // Load in arrays
    double* features = mxGetPr(input[0]);
    int* frames = (int*) mxGetData(input[1]);
    int* im_size = (int*) mxGetData(input[2]);
    int* cell_size1 = (int*) mxGetData(input[3]);
    int* sppyrs = (int*) mxGetData(input[4]);
    double* svm_coeffs = mxGetPr(input[5]);
    int* boxes = (int*) mxGetData(input[6]);
    int* cell_size2 = (int*) mxGetData(input[8]);
     
    char* filename = mxArrayToString(input[7]);
    // The features are column wise. First value is the feature dimension
    int* featSize = (int*) mxGetDimensions(input[0]);
    int dimensionality = featSize[0];
    int nrFeatures = featSize[1];
    int* pyrSize = (int*) mxGetDimensions(input[4]);
    int nrPyrs = pyrSize[1];
    int* boxSize = (int*) mxGetDimensions(input[6]);
    int nrBoxes = boxSize[1];
    int* modelSize = (int*) mxGetDimensions(input[5]);
    int nrCoeffs = modelSize[0];
    int nrModels = modelSize[1];
       
    // initialize gmm model class
    gaussian_mixture<double> *gmm;
    gmm = new gaussian_mixture<double>(filename);
    
     // initialize fisher_net class
    fisher_net_param param;
    // remove Power and L2 Normalization
    param.alpha = 1.0f;
    param.pnorm = 0.0;
    
    fisher_net<double>* fisherNet;
    fisherNet = new fisher_net<double>(param);
    fisherNet->set_model(*gmm);

    // Now make vector arrays from input arguments (which are double arrays)
    std::vector<double*> featVecs(nrFeatures);
    for(int i=0; i < nrFeatures; i++){
        featVecs[i] = features + (i * dimensionality);
    }
    std::vector<int*> frameVecs(nrFeatures);
    for(int i=0; i < nrFeatures; i++){
        frameVecs[i] = frames + (i * 2);
    }
    std::vector<int*> pyrVecs(nrPyrs);
    for(int i=0; i <nrPyrs; i++){
        pyrVecs[i] = sppyrs + (i * 2);
    }
    std::vector<int*> boxVecs(nrBoxes);
    for(int i=0; i <nrBoxes; i++){
        boxVecs[i] = boxes + (i * 4);
    }
    std::vector<double*> modelVecs(nrModels);
    for(int i=0; i <nrModels; i++){
        modelVecs[i] = svm_coeffs + (i * nrCoeffs);
    }

    // Allocate memory for returning fisher vector
	out[0] = mxCreateDoubleMatrix(nrBoxes, nrModels, mxREAL);
    double* scores = mxGetPr(out[0]);
    std::vector<double*> scoreVecs(nrModels);
    for(int i=0; i<nrModels; ++i){
        scoreVecs[i] = scores + (i * nrBoxes);   
    }

    //fisherDetector->compute(featVecs, frameVecs, im_size, pyrVecs, svm_coeffs, boxVecs, scores, cell_size);
    fisher_approx_detector_compute(fisherNet, featVecs, frameVecs, im_size, pyrVecs, modelVecs, boxVecs, scoreVecs, cell_size1, cell_size2);
    
    // Free memory
    //delete fisherDetector;
    return;
}
 
