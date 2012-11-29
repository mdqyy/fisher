
#include "fisher_detector.h"
#include "simd_math.h"
#include <time.h>
#include "mex.h"


double diffclock(clock_t clock1,clock_t clock2)
{
  double diffticks=clock1-clock2;
  double diffms=(diffticks)/CLOCKS_PER_SEC;
  return diffms;
}


template<class T>
fisher_detector<T>::fisher_detector( char* gmmfile )
{
  // initialize gmm model class
  gmm = new gaussian_mixture<T>(gmmfile);

  // initialize fisher net class
  fisher_net_param param;
  // remove Power and L2 Normalization
  param.alpha = 1.0f;
  param.pnorm = 0.0;
	
  fnet = new fisher_net<T>(param);
  fnet->set_model(*gmm);
}


template<class T>
fisher_detector<T>::~fisher_detector()
{
  delete gmm;
  gmm = 0;

  delete fnet;
  fnet = 0;
}


/////////////////////////////////////////////////////
//description of the parameters:
//	x = vector({sift1}, {sift2}, ...)
//	lc = vector({x1, y1}, {x2, y2}, ...)
//	im_size = {width, height}
//	sppyr = vector({x_partitions_1 y_partitions_1}, {x_partitions_2 y_partitions_2}, ...)
//	w = {w1, w2, ...}
//	box = vector({x1, y1, x2, y2}, {}, ...)
//	score = {score1, score2, ...}
//	cell_size = {cell_width, cell_height}
/////////////////////////////////////////////////////
template<class T>
int 
fisher_detector<T>::compute(std::vector<T*> &x, std::vector<int*> &lc, int* im_size, std::vector<int*> &sppyr,
			          T* w, std::vector<int*> &box, T* score, int* cell_size)
{
  assert(fnet);

  // to do..
  int cols = (int)ceil(im_size[0] / (float)cell_size[0]);
  int rows = (int)ceil(im_size[1] / (float)cell_size[1]);
  int ncells = rows*cols;
  
  int npyrs = 1;
  int lpyrs = (int)sppyr.size();
  for(int i=0; i<lpyrs; ++i){
    npyrs += sppyr[i][0] * sppyr[i][1];
  }

  T* fmap_pr = new T[ncells * npyrs];
  std::vector<T*> fmap(npyrs);
  for(int i=0; i<npyrs; ++i){
    fmap[i] = fmap_pr + (i * ncells);
  }
  T* fnorm = new T[ncells * ncells];

  // compute fisher map
  compute_fisher_map(x, lc, im_size, cell_size, sppyr, w, fmap, fnorm);

  // compute detection score per box
  int nboxes = (int)box.size();

  #pragma omp parallel for
  for(int i=0; i<nboxes; ++i)
  {
    int* box_i = box[i];
	int mask[4];
	// compute the whole
	compute_map_mask(cell_size[0], cell_size[1], box_i, mask);
    T score_i = compute_detector_score(cols, rows, fmap[0], fnorm, mask);
    
	int ss = 1;
	// #pragma omp parallel for
    for(int j=0; j<lpyrs; ++j)
	{
	  int* pyr_j = sppyr[j];
      T stepx = (box_i[2] - box_i[0]) / (T)pyr_j[0];
      T stepy = (box_i[3] - box_i[1]) / (T)pyr_j[1];

	  for(int p=0; p<pyr_j[0]; ++p)
	  {
        for(int q=0; q<pyr_j[1]; ++q)
		{
		  int roi[4];
		  int roi_mask[4];
		  // compute region of interest
	      roi[0] = (int)ceil(box_i[0] + p*stepx);
	      roi[1] = (int)ceil(box_i[1] + q*stepy);
          roi[2] = (int)floor(box_i[0] + (p+1)*stepx);
	      roi[3] = (int)floor(box_i[1] + (q+1)*stepy);
         
		  // compute map mask
          compute_map_mask(cell_size[0], cell_size[1], roi, roi_mask);
  
          // compute map score
          T score_i_ss = compute_detector_score(cols, rows, fmap[ss++], fnorm, roi_mask);
          score_i += score_i_ss;
		}
	  }
	}
    
	// normalize score
	score[i] = score_i / sqrt((T)npyrs);
  }

  delete[] fnorm;
  delete[] fmap_pr;

  return 0;
}


template<class T>
void
fisher_detector<T>::compute_map_mask( int cell_width, int cell_height, int* roi, int* mask )
{
  mask[0] = (int)ceil(roi[0] / (float)cell_width) - 1;
  mask[2] = (int)ceil(roi[2] / (float)cell_width) - 1;
  mask[1] = (int)ceil(roi[1] / (float)cell_height) - 1;
  mask[3] = (int)ceil(roi[3] / (float)cell_height) - 1;
}


template<class T>
T 
fisher_detector<T>::compute_detector_score( int cols, int rows, T* map, T* norm, int* mask )
{
  T score = 0;
  T scale = 0;
  int ncells = rows*cols;

  //#pragma omp parallel for reduction(+:score)
  //for(int i=mask[0]; i<=mask[2]; ++i)
  //  for(int j=mask[1]; j<=mask[3]; ++j)
	 // score += map[i*rows + j];

  score = compute_map_score(cols, rows, map, mask);

  //#pragma omp parallel for reduction(+:scale)
  //for(int i=mask[0]; i<=mask[2]; ++i)
  //  for(int j=mask[1]; j<=mask[3]; ++j)
	 // for(int p=mask[0]; p<=mask[2]; ++p)
  //      for(int q=mask[1]; q<=mask[3]; ++q)
	 //     scale += norm[(i*rows+j)*ncells + (p*rows+q)];

  #pragma omp parallel for reduction(+:scale)
  for(int i=mask[0]; i<=mask[2]; ++i)
    for(int j=i; j<=mask[2]; ++j){
      int norm_mask[4];
	  norm_mask[0] = i*rows + mask[1];
	  norm_mask[1] = j*rows + mask[1];
	  norm_mask[2] = i*rows + mask[3];
	  norm_mask[3] = j*rows + mask[3];

	  T scale_ij = compute_map_score(ncells, ncells, norm, norm_mask);
	  scale += scale_ij;

	  if(j != i)
	    scale += scale_ij;
    }

  return (score / sqrt((T)scale));
}


template<class T>
T 
fisher_detector<T>::compute_map_score( int cols, int rows, T* map, int* mask )
{
  T score = 0;

  score += map[mask[2]*rows + mask[3]];
  score -= (mask[0] ? map[(mask[0]-1)*rows + mask[3]] : 0);
  score -= (mask[1] ? map[mask[2]*rows + (mask[1]-1)] : 0);
  score += ((mask[0] && mask[1]) ? map[(mask[0]-1)*rows + (mask[1]-1)] : 0);

  return score;
}


template<class T>
void 
fisher_detector<T>::compute_fisher_map(std::vector<T*> &x, std::vector<int*> &lc, int* im_size, int* cell_size,
                      std::vector<int*> &sppyr, T* w, std::vector<T*> &fmap, T* fnorm)
{
  // compute cell assign per descriptor
  int cols = (int)ceil(im_size[0] / (float)cell_size[0]);
  int rows = (int)ceil(im_size[1] / (float)cell_size[1]);
  int ndescrs = (int)x.size();
  int* idx = new int[ndescrs];
  int ncells = compute_node_assign(lc, im_size[0], im_size[1], cell_size[0], cell_size[1], idx);

  // computer fisher vectors
  int fkdim = fnet->dim();
  T* fisher_pr = new T[fkdim * ncells];
  std::vector<T*> fisher(ncells);
  for(int i=0; i<ncells; ++i){
    fisher[i] = fisher_pr + (i * fkdim);
  }

  clock_t begin=clock();
  fnet->compute_fisher_net(x, idx, fisher);
  clock_t end=clock();
  mexPrintf("Time elapsed: %f ms\n", 
                 double(diffclock(end, begin)));

  // computer fisher norms
  begin=clock();
  vector_mult(fisher, fkdim, fnorm);
  end=clock();
  mexPrintf("Time elapsed: %f ms\n", 
                 double(diffclock(end, begin)));

  // compute integral norm map
  compute_integral_map( ncells, ncells, fnorm );

  // compute fisher maps with spatial pyramids
  int npyrs = (int)fmap.size();

  #pragma omp parallel for
  for(int i=0; i<npyrs; ++i)
  {
    T* w_i = w + i*fkdim;
	T* fmap_i = fmap[i];

	// compute cell-wise contribution to classification
	for(int j=0; j<ncells; ++j)
	{
      fmap_i[j] = simd::dot_product( fkdim, w_i, fisher[j] );
	}

	// compute integral score map
    compute_integral_map( cols, rows, fmap_i );
  }
  
  delete[] fisher_pr;
  delete[] idx;
}


template<class T>
void
fisher_detector<T>::compute_integral_map( int cols, int rows, T* map )
{
  // first col only
  T rs = 0;
  for(int i=0; i<rows; ++i)
  {
    rs+= map[i];
	map[i] = rs;
  }

  // remaining cells are sum above and to the left
  for(int i=1; i<cols; ++i)
  {
    rs = 0;
    for(int j=0; j<rows; ++j)
    {
      rs += map[i*rows + j];
	  map[i*rows + j] = rs + map[(i-1)*rows + j];
    }
  }
}


template<class T>
int
fisher_detector<T>::compute_node_assign(std::vector<int*> &lc, int width, int height,
                      int cell_width, int cell_height, int* idx)
{
  int ndescrs = (int)lc.size();
  int cols = (int)ceil(width / (float)cell_width);
  int rows = (int)ceil(height / (float)cell_height);

  #pragma omp parallel for
  for(int i=0; i<ndescrs; ++i)
  {
    int* lc_i = lc[i];
	int binx = (int)ceil(lc_i[0] / (float)cell_width) - 1;
    int biny = (int)ceil(lc_i[1] / (float)cell_height) - 1;
	idx[i] = binx*rows + biny;
  }

  return rows*cols;
}


template<class T>
void
vector_mult(std::vector<T*> &vec, int dim, T* result)
{
  int nvecs = (int)vec.size();

  #pragma omp parallel for
  for(int i=0; i<nvecs; ++i)
  {
	for(int j=i; j<nvecs; ++j)
	{
	  T dp = simd::dot_product( dim, vec[i], vec[j] );

	  // symmetric matrix
	  result[i*nvecs + j] = dp;
	  result[j*nvecs + i] = dp;
	}
  }
}



template void vector_mult<float>( std::vector<float*> &vec, int dim, float* result );
template void vector_mult<double>( std::vector<double*> &vec, int dim, double* result );

template class fisher_detector<float>;
template class fisher_detector<double>;
