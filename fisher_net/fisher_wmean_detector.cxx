
#include "fisher_detector.h"
#include "simd_math.h"

#ifndef EXP_MU
#define EXP_MU 0.013
#endif

/////////////////////////////////////////////////////
//description of the parameters:
//    x = vector({sift1}, {sift2}, ...)
//    lc = vector({x1, y1}, {x2, y2}, ...)
//    im_size = {width, height}
//    sppyr = vector({x_partitions_1 y_partitions_1}, {x_partitions_2 y_partitions_2}, ...)
//    w = {w1, w2, ...}
//    boxes = vector({x1, y1, x2, y2}, {}, ...)
//    scores = {score1, score2, ...}
//    cell_size = {cell_width, cell_height}
/////////////////////////////////////////////////////
template<class T>
int 
fisher_detector_compute(fisher_net<T> *fnet, std::vector<T*> &x, std::vector<int*> &lc, int* im_size, std::vector<int*> &sppyr,
                      std::vector<T*>& models, std::vector<int*> &boxes, std::vector<T*> &scores, int* cell_size)
{
  assert(fnet);
  
  T* w = models[0];
  T* score = scores[0];

  // to do..
  int cols = (int)ceil(im_size[0] / (float)cell_size[0]);
  int rows = (int)ceil(im_size[1] / (float)cell_size[1]);
  int ncells = rows*cols;
  std::cout << ncells << std::endl;
  
  int npyrs = 1;
  int lpyrs = (int)sppyr.size();
  for(int i=0; i<lpyrs; ++i){
    npyrs += sppyr[i][0] * sppyr[i][1];
  }

  int nboxes = (int)boxes.size();
  int nmodels = (int)models.size();
  T* fmap_pr = new T[nmodels * npyrs * ncells];
  std::vector<T*> fmap(nmodels);
  for(int i=0; i<nmodels; ++i){
    fmap[i] = fmap_pr + (i * npyrs * ncells);

    // initialize scores
    //memset( scores[i], 0, nboxes*sizeof(T));
  }
  T* fnorm = new T[ncells * ncells];

  // compute fisher map
  fisher_detector_compute_fisher_map(fnet, x, lc, im_size, cell_size, npyrs, models, fmap, fnorm);

  // compute detection score per box
  #pragma omp parallel for
  for(int i=0; i<nboxes; ++i)
  {
    int* box_i = boxes[i];
    int mask[4];
    // compute the whole
    fisher_detector_compute_map_mask(cell_size[0], cell_size[1], box_i, mask);

    // compute map score for each model
    for(int k=0; k<nmodels; ++k)
    {
      scores[k][i] = fisher_detector_compute_detector_score(cols, rows, fmap[k], fnorm, mask);
    }
    
    int ss = 1;
    // #pragma omp parallel for
    for(int j=0; j<lpyrs; ++j)
    {
      int* pyr_j = sppyr[j];
      T stepx = (box_i[2] - box_i[0] + 1) / (T)pyr_j[0];
      T stepy = (box_i[3] - box_i[1] + 1) / (T)pyr_j[1];


      //for(int q=0; q<pyr_j[1]; ++q)
      for(int p=0; p<pyr_j[0]; ++p)
      {
        //for(int p=0; p<pyr_j[0]; ++p)
        for(int q=0; q<pyr_j[1]; ++q)
        {
          int roi[4];
          int roi_mask[4];
          // compute region of interest
          //roi[0] = (int)ceil(box_i[0] + p*stepx);
          //roi[1] = (int)ceil(box_i[1] + q*stepy);
          //roi[2] = (int)floor(box_i[0] + (p+1)*stepx);
          //roi[3] = (int)floor(box_i[1] + (q+1)*stepy);
          roi[0] = (int)(box_i[0] + p*stepx);
          roi[1] = (int)(box_i[1] + q*stepy);
          roi[2] = (int)(box_i[0] + (p+1)*stepx - 1);
          roi[3] = (int)(box_i[1] + (q+1)*stepy - 1);
         
          // compute map mask
          fisher_detector_compute_map_mask(cell_size[0], cell_size[1], , roi_mask);
  
          // compute map score for each model
          for(int k=0; k<nmodels; ++k)
          {
            T* fmap_k_ss = fmap[k] + ss*ncells;
            T score_i_ss = fisher_detector_compute_detector_score(cols, rows, fmap_k_ss, fnorm, roi_mask);
            scores[k][i] += score_i_ss;
          }

          ss++;
        }
      }
    }
    
    // normalize score
    for(int k=0; k<nmodels; ++k){
      scores[k][i] /= sqrt((T)npyrs);
    }
  }

  delete[] fnorm;
  delete[] fmap_pr;

  return 0;
}


void
fisher_detector_compute_map_multi_masks( int width, int height, int cell_width, int cell_height, int* roi, std::vector<int*>& masks, std::vector<int*>& maskboxes )
{
  int binsx[2], binsy[2];
  int near_neighbors[4];

  binsx[0] = (int)ceil(roi[0] / (float)cell_width) - 1;
  binsx[1] = (int)ceil(roi[2] / (float)cell_width) - 1;
  binsy[0] = (int)ceil(roi[1] / (float)cell_height) - 1;
  binsy[1] = (int)ceil(roi[3] / (float)cell_height) - 1;

  near_neighbors[0] = (roi[0]%cell_width == 1 || (roi[0]+cell_width) > (width+1) ? 0 : 1;
  near_neighbors[1] = (roi[1]%cell_height == 1 || (roi[1]+cell_height) > (height+1) ? 0 : 1;
  near_neighbors[2] = (roi[2]%cell_width == 0 || (roi[2]-cell_width) < 0) ? 0 : 1;
  near_neighbors[3] = (roi[3]%cell_height == 0 || (roi[3]-cell_height) < 0) ? 0 : 1;

  for(int x0 = 0; x0<=near_neighbors(0); ++x0)
    for(int x1 = 0; x1<=near_neighbors(2); ++x1)
      for(int y0 = 0; y0<=near_neighbors(1); ++y0)
        for(int y1 = 0; y1<=near_neighbors(3); ++y1){
          int mask[4];
		  int maskbox[4];
          mask[0] = binsx[0] + x0;
          mask[1] = binsy[0] + y0;
          mask[2] = binsx[1] - x1;
          mask[3] = binsy[1] - y1;
          masks.push_back(mask);

		  maskbox[0] = mask[0]*cell_width + 1;
		  maskbox[1] = mask[0]*cell_height + 1;
        }
}


template<class T>
T 
fisher_detector_compute_detector_approx_score( int cols, int rows, T* map, T* norm, int* roi, std::vector<int*>& masks, std::vector<int*>& maskboxes )
{
  int nmasks = (int)masks.size();
  if( nmasks == 1 ){
    fisher_detector_compute_detector_score(cols, rows, map, norm, masks[0]);
  }else{
    //ovlps
	for(int i=0; i<nmasks; ++i)
	{
	  int* mask = masks[i];
	  
	  
	}
  }
}


void
fisher_detector_compute_map_mask( int cell_width, int cell_height, int* roi, int* mask )
{
  mask[0] = (int)ceil(roi[0] / (float)cell_width) - 1;
  mask[2] = (int)ceil(roi[2] / (float)cell_width) - 1;
  mask[1] = (int)ceil(roi[1] / (float)cell_height) - 1;
  mask[3] = (int)ceil(roi[3] / (float)cell_height) - 1;
}


template<class T>
T 
fisher_detector_compute_detector_score( int cols, int rows, T* map, T* norm, int* mask )
{
  T score = 0;
  T scale = 0;
  int ncells = rows*cols;

  // compute score
  score = fisher_detector_compute_map_score(cols, rows, map, mask);

  // compute norm scalar
  #pragma omp parallel for reduction(+:scale)
  for(int i=mask[0]; i<=mask[2]; ++i)
    for(int j=i; j<=mask[2]; ++j){
      int norm_mask[4];
      norm_mask[0] = i*rows + mask[1];
      norm_mask[1] = j*rows + mask[1];
      norm_mask[2] = i*rows + mask[3];
      norm_mask[3] = j*rows + mask[3];

      T scale_ij = fisher_detector_compute_map_score(ncells, ncells, norm, norm_mask);
      scale += scale_ij;

      if(j != i)
        scale += scale_ij;
    }

  if(scale < 1e-5) return 0.0;
  return (score / sqrt((T)scale));
}


template<class T>
T 
fisher_detector_compute_map_score( int cols, int rows, T* map, int* mask )
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
fisher_detector_compute_fisher_map( fisher_net<T> *fnet, std::vector<T*> &x, std::vector<int*> &lc, int* im_size, int* cell_size,
                      int npyrs, std::vector<T*> &models, std::vector<T*> &fmap, T* fnorm )
{
  // compute cell assign per descriptor
  int cols = (int)ceil(im_size[0] / (float)cell_size[0]);
  int rows = (int)ceil(im_size[1] / (float)cell_size[1]);
  int ndescrs = (int)x.size();
  int* idx = new int[ndescrs];
  int ncells = fisher_detector_compute_node_assign(lc, im_size[0], im_size[1], cell_size[0], cell_size[1], idx);

  // computer fisher vectors
  int fkdim = fnet->dim();
  T* fisher_pr = new T[fkdim * ncells];
  std::vector<T*> fisher(ncells);
  for(int i=0; i<ncells; ++i){
    fisher[i] = fisher_pr + (i * fkdim);
  }
  fnet->compute_fisher_net(x, idx, fisher);

  // computer fisher norms
  vector_mult(fisher, fkdim, fnorm);

  // compute integral norm map
  fisher_detector_compute_integral_map( ncells, ncells, fnorm );

  // compute fisher maps with spatial pyramids
  //int npyrs = (int)fmap.size();
  int nmodels = (int)models.size();
  
  #pragma omp parallel for
  for(int k=0; k<nmodels; ++k)
  {
    T* model_k = models[k];
    T* fmap_k = fmap[k];
    
    for(int i=0; i<npyrs; ++i)
    {
      T* w_k_i = model_k + i*fkdim;
      T* fmap_k_i = fmap_k + i*ncells;

      // compute cell-wise contribution to classification
      for(int j=0; j<ncells; ++j)
      {
        fmap_k_i[j] = simd::dot_product( fkdim, w_k_i, fisher[j] );
      }

      // compute integral score map
      fisher_detector_compute_integral_map( cols, rows, fmap_k_i );
    }    
  }

  delete[] fisher_pr;
  delete[] idx;
}


template<class T>
void
fisher_detector_compute_integral_map( int cols, int rows, T* map )
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


int
fisher_detector_compute_node_assign(std::vector<int*> &lc, int width, int height,
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
  //using namespace Impala::Core::Matrix;
  int nvecs = (int)vec.size();

  std::cout << "vecmul" << dim << " " << nvecs << std::endl;
/*  Mat* m = MatCreate<Mat>(nvecs, dim);
  for(int i=0; i<nvecs; ++i)
  {
      for(int j=0; j<dim; ++j)
      {
          *MatE(m, i, j) = vec[i][j];
      }
  }
  Mat* mT = MatTranspose(m);
  std::cout << "vecmulbe" << dim << " " << nvecs << std::endl;
  Mat* res = MatMul(m, mT);
  delete m;
  delete mT;
  std::cout << "vecmulaf" << dim << " " << nvecs << std::endl;
*/
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
  std::cout << "vecmul" << dim << std::endl;
}


template int fisher_detector_compute( fisher_net<float> *fnet, std::vector<float*> &x, std::vector<int*> &lc, int* im_size, std::vector<int*> &sppyr,
                            std::vector<float*>& models, std::vector<int*> &boxes, std::vector<float*> &scores, int* cell_size );

template int fisher_detector_compute( fisher_net<double> *fnet, std::vector<double*> &x, std::vector<int*> &lc, int* im_size, std::vector<int*> &sppyr,
                            std::vector<double*>& models, std::vector<int*> &boxes, std::vector<double*> &scores, int* cell_size );

template void vector_mult<float>( std::vector<float*> &vec, int dim, float* result );
template void vector_mult<double>( std::vector<double*> &vec, int dim, double* result );
