// Wrapper made by Zhenyang Li
#include <cmath>
#include "mex.h"
#include "fisher_detector.h"

void mexFunction(int nlhs, mxArray *out[], int nrhs, const mxArray *input[])
{
    // Checking number of arguments
    if (nlhs > 1){
        mexErrMsgTxt("Function has a single output");
        return;
    }

    if (nrhs != 8){
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
    int* cell_size = (int*) mxGetData(input[3]);
    int* sppyrs = (int*) mxGetData(input[4]);
    double* svm_coefs = mxGetPr(input[5]);
    int* boxes = (int*) mxGetData(input[6]);
    
    char* filename = mxArrayToString(input[7]);
    // The features are column wise. First value is the feature dimension
    int* featSize = (int*) mxGetDimensions(input[0]);
    int dimensionality = featSize[0];
    int nrFeatures = featSize[1];
    int* pyrSize = (int*) mxGetDimensions(input[4]);
    int nrPyrs = pyrSize[1];
    int* boxSize = (int*) mxGetDimensions(input[6]);
    int nrBoxes = boxSize[1];
       
    // Initialize fisher_detector class
    fisher_detector<double>* fisherDetector;
    fisherDetector = new fisher_detector<double>(filename);

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
    for(int i=0; i <nrPyrs ; i++){
        pyrVecs[i] = sppyrs + (i * 2);
    }
    std::vector<int*> boxVecs(nrBoxes);
    for(int i=0; i <nrBoxes ; i++){
        boxVecs[i] = boxes + (i * 4);
    }

    // Allocate memory for returning fisher vector
	out[0] = mxCreateDoubleMatrix(nrBoxes, 1, mxREAL);
    double* scores = mxGetPr(out[0]);

    fisherDetector->compute(featVecs, frameVecs, im_size, pyrVecs, svm_coefs, boxVecs, scores, cell_size);

    // Free memory
    delete fisherDetector;
    return;
}
 
