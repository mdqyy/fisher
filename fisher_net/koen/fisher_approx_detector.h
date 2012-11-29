
/// \Class    fisher_approx_detector "fisher_approx_detector.h"
///  
/// \brief
/// \

#ifndef __FISHER_APPROX_DETECTOR_H
#define __FISHER_APPROX_DETECTOR_H


#include "gmm.h"
#include "fisher_net.h"



template<class T>
int fisher_approx_detector_compute( fisher_net<T> *fnet, std::vector<T*> &x, std::vector<int*> &lc, int* im_size, std::vector<int*> &sppyr,
                       std::vector<T*>& models, std::vector<int*> &boxes, std::vector<T*> &scores, int* cell_size );

template<class T>
int fisher_approx_detector_compute( fisher_net<T> *fnet, std::vector<T*> &x, std::vector<int*> &lc, int* im_size, std::vector<int*> &sppyr,
                       std::vector<T*>& models, std::vector<int*> &boxes, std::vector<T*> &scores, int* cell_size1, int* cell_size2 );



#endif
