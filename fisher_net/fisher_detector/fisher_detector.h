
/// \Class    fisher_detector "fisher_detector.h"
///  
/// \brief
/// \

#ifndef __FISHER_DETECTOR_H
#define __FISHER_DETECTOR_H


#include "gmm.h"
#include "fisher_net.h"



template<class T>
class fisher_detector
{
public:

  fisher_detector( char* gmmfile );
  ~fisher_detector( );
  
  int compute( std::vector<T*> &x, std::vector<int*> &lc, int* im_size, std::vector<int*> &sppyr,
        T* w, std::vector<int*> &box, T* score, int* cell_size );

private:

  T compute_detector_score( int cols, int rows, T* map, T* norm, int* mask );
  
  T compute_map_score( int cols, int rows, T* map, int* mask );

  void compute_map_mask( int cell_width, int cell_height, int* roi, int* mask );

  void compute_fisher_map( std::vector<T*> &x, std::vector<int*> &lc, int* im_size, int* cell_size, 
         std::vector<int*> &sppyr, T* w, std::vector<T*> &fmap, T* fnorm );

  void compute_integral_map( int cols, int rows, T* map );

  int compute_node_assign( std::vector<int*> &lc, int width, int height,
	    int cell_width, int cell_height, int* idx );

protected:
  
  gaussian_mixture<T> *gmm;
  fisher_net<T> *fnet;

};


#endif