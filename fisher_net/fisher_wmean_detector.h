
/// \Class    fisher_detector "fisher_detector.h"
///  
/// \brief
/// \

#ifndef __FISHER_DETECTOR_H
#define __FISHER_DETECTOR_H


#include "gmm.h"
#include "fisher_net.h"


template<class T>
int fisher_detector_compute(fisher_net<T> *fnet, std::vector<T*> &x, std::vector<int*> &lc, int* im_size, std::vector<int*> &sppyr,
							std::vector<T*>& models, std::vector<int*> &box, std::vector<T*> &scores, int* cell_size );

template<class T>
T fisher_detector_compute_detector_score( int cols, int rows, T* map, T* norm, int* mask );
 
template<class T>
T fisher_detector_compute_map_score( int cols, int rows, T* map, int* mask );

void fisher_detector_compute_map_mask( int cell_width, int cell_height, int* roi, int* mask );

template<class T>
void fisher_detector_compute_fisher_map( fisher_net<T> *fnet, std::vector<T*> &x, std::vector<int*> &lc, int* im_size, int* cell_size, 
        int npyrs, std::vector<T*> &models, std::vector<T*> &fmap, T* fnorm );

template<class T>
void fisher_detector_compute_integral_map( int cols, int rows, T* map );

int fisher_detector_compute_node_assign( std::vector<int*> &lc, int width, int height,
        int cell_width, int cell_height, int* idx );


#endif
