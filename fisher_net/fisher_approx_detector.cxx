
#include "fisher_approx_detector.h"
#include "simd_math.h"

#define DIFF_EPSILON 1.19209290E-07F

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
fisher_approx_detector_compute( fisher_net<T> *fnet, std::vector<T*> &x, std::vector<int*> &lc, int* im_size, std::vector<int*> &sppyr,
                   std::vector<T*>& models, std::vector<int*> &boxes, std::vector<T*> &scores, int* cell_size )
{
  assert(fnet);
  
  T* w = models[0];
  T* score = scores[0];

  // to do ...
  int cols = (int)ceil(im_size[0] / (float)cell_size[0]);
  int rows = (int)ceil(im_size[1] / (float)cell_size[1]);
  int ncells = rows*cols;
  
  int npyrs = 1;
  int lpyrs = (int)sppyr.size();
  for(int i=0; i<lpyrs; ++i){
    npyrs += sppyr[i][0] * sppyr[i][1];
  }

  int nboxes = (int)boxes.size();
  int nmodels = (int)models.size();
  T* fmap_pr = new T[nmodels * npyrs * ncells];
  T* fnorm_pr = new T[ncells * ncells];
  std::vector<T*> fmap(nmodels);
  std::vector<T*> fnorm(ncells);
  for(int i=0; i<nmodels; ++i){
    fmap[i] = fmap_pr + (i * npyrs * ncells);

	// initialize scores
	//memset(scores[i], 0, nboxes*nmodels*sizeof(T));
  }
  for(int i=0; i<ncells; ++i){
	fnorm[i] = fnorm_pr + (i * ncells);
  }

  // compute fisher map
  compute_fisher_map(fnet, x, lc, im_size, cell_size, npyrs, models, fmap, fnorm);

  // compute detection score per box
  #pragma omp parallel for
  for(int i=0; i<nboxes; ++i)
  {
    int* box_i = boxes[i];
    int mask[4];
    T norm_scalar, wspacing[4];

    // compute the whole
    compute_map_mask(cell_size[0], cell_size[1], box_i, mask, wspacing);
    norm_scalar = compute_map_norm(cols, rows, fnorm, mask, wspacing);

	// compute map score for each model
	for(int k=0; k<nmodels; ++k)
    {
      scores[k][i] = compute_detector_score(cols, rows, fmap[k], norm_scalar, mask, wspacing);
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
		  int roi[4], roi_mask[4];
          T roi_norm_scalar, roi_wspacing[4];
          
          // compute region of interest
		  roi[0] = (int)ceil(box_i[0] + p*stepx);
          roi[1] = (int)ceil(box_i[1] + q*stepy);
          roi[2] = (int)ceil(box_i[0] + (p+1)*stepx - 1);
          roi[3] = (int)ceil(box_i[1] + (q+1)*stepy - 1);
		  roi[2] = (roi[2] < box_i[2]) ? roi[2] : box_i[2];
		  roi[3] = (roi[3] < box_i[3]) ? roi[3] : box_i[3];
         
          // compute map mask
          compute_map_mask(cell_size[0], cell_size[1], roi, roi_mask, roi_wspacing);
          roi_norm_scalar = compute_map_norm(cols, rows, fnorm, roi_mask, roi_wspacing);
  
          // compute map score for each model
		  for(int k=0; k<nmodels; ++k)
	      {
			T* fmap_k_ss = fmap[k] + ss*ncells;
            T score_i_ss = compute_detector_score(cols, rows, fmap_k_ss, roi_norm_scalar, roi_mask, roi_wspacing);
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

  delete[] fnorm_pr;
  delete[] fmap_pr;
  return 0;
}


template<class T>
int 
fisher_approx_detector_compute( fisher_net<T> *fnet, std::vector<T*> &x, std::vector<int*> &lc, int* im_size, std::vector<int*> &sppyr,
                   std::vector<T*>& models, std::vector<int*> &boxes, std::vector<T*> &scores, int* cell_size1, int* cell_size2 )
{
  assert(fnet);

  T* w = models[0];
  T* score = scores[0];

  // to do ...
  int cols1 = (int)ceil(im_size[0] / (float)cell_size1[0]);
  int rows1 = (int)ceil(im_size[1] / (float)cell_size1[1]);
  int cols2 = (int)ceil(im_size[0] / (float)cell_size2[0]);
  int rows2 = (int)ceil(im_size[1] / (float)cell_size2[1]);
  int ncells1 = rows1*cols1;
  int ncells2 = rows2*cols2;

  int npyrs = 1;
  int lpyrs = (int)sppyr.size();
  for(int i=0; i<lpyrs; ++i){
    npyrs += sppyr[i][0] * sppyr[i][1];
  }

  int nboxes = (int)boxes.size();
  int nmodels = (int)models.size();
  T* fmap_pr = new T[nmodels * npyrs * ncells1];
  T* fnorm_pr = new T[ncells2 * ncells2];
  std::vector<T*> fmap(nmodels);
  std::vector<T*> fnorm(ncells2);
  for(int i=0; i<nmodels; ++i){
    fmap[i] = fmap_pr + (i * npyrs * ncells1);

    // initialize scores
    //memset(scores[i], 0, nboxes*nmodels*sizeof(T));
  }
  for(int i=0; i<ncells2; ++i){
    fnorm[i] = fnorm_pr + (i * ncells2);
  }

  // compute fisher map
  compute_fisher_map(fnet, x, lc, im_size, cell_size1, cell_size2, npyrs, models, fmap, fnorm);

  // compute detection score per box
  #pragma omp parallel for
  for(int i=0; i<nboxes; ++i)
  {
    int* box_i = boxes[i];
    int mask1[4], mask2[4];
    T norm_scalar, wspacing1[4], wspacing2[4];

    // compute the whole
    compute_map_mask(cell_size1[0], cell_size1[1], box_i, mask1, wspacing1);
    compute_map_mask(cell_size2[0], cell_size2[1], box_i, mask2, wspacing2);
    norm_scalar = compute_map_norm(cols2, rows2, fnorm, mask2, wspacing2);

    // compute map score for each model
    for(int k=0; k<nmodels; ++k)
    {
      scores[k][i] = compute_detector_score(cols1, rows1, fmap[k], norm_scalar, mask1, wspacing1);
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
          int roi[4], roi_mask1[4], roi_mask2[4];
          T roi_norm_scalar, roi_wspacing1[4], roi_wspacing2[4];

          // compute region of interest
          roi[0] = (int)ceil(box_i[0] + p*stepx);
          roi[1] = (int)ceil(box_i[1] + q*stepy);
          roi[2] = (int)ceil(box_i[0] + (p+1)*stepx - 1);
          roi[3] = (int)ceil(box_i[1] + (q+1)*stepy - 1);
          roi[2] = (roi[2] < box_i[2]) ? roi[2] : box_i[2];
          roi[3] = (roi[3] < box_i[3]) ? roi[3] : box_i[3];

          // compute map mask
          compute_map_mask(cell_size1[0], cell_size1[1], roi, roi_mask1, roi_wspacing1);
          compute_map_mask(cell_size2[0], cell_size2[1], roi, roi_mask2, roi_wspacing2);
          roi_norm_scalar = compute_map_norm(cols2, rows2, fnorm, roi_mask2, roi_wspacing2);

          // compute map score for each model
          for(int k=0; k<nmodels; ++k)
          {
            T* fmap_k_ss = fmap[k] + ss*ncells1;
            T score_i_ss = compute_detector_score(cols1, rows1, fmap_k_ss, roi_norm_scalar, roi_mask1, roi_wspacing1);
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

  delete[] fnorm_pr;
  delete[] fmap_pr;
  return 0;
}


inline
void
compute_map_mask( int cell_width, int cell_height, int* roi, int* mask )
{
  mask[0] = (int)ceil(roi[0] / (float)cell_width) - 1;
  mask[2] = (int)ceil(roi[2] / (float)cell_width) - 1;
  mask[1] = (int)ceil(roi[1] / (float)cell_height) - 1;
  mask[3] = (int)ceil(roi[3] / (float)cell_height) - 1;
}


template<class T>
inline
void
compute_map_mask( int cell_width, int cell_height, int* roi, int* mask, T* wspacing )
{
  mask[0] = (int)ceil(roi[0] / (float)cell_width) - 1;
  mask[2] = (int)ceil(roi[2] / (float)cell_width) - 1;
  mask[1] = (int)ceil(roi[1] / (float)cell_height) - 1;
  mask[3] = (int)ceil(roi[3] / (float)cell_height) - 1;

  wspacing[0] = (roi[0] - 1) / (T)cell_width - mask[0];
  wspacing[2] = (mask[2] + 1) - roi[2] / (T)cell_width;
  wspacing[1] = (roi[1] - 1) / (T)cell_height - mask[1];
  wspacing[3] = (mask[3] + 1) - roi[3] / (T)cell_height;
}


template<class T>
inline
T
compute_detector_score( int cols, int rows, T* map, T norm_scalar, int* mask )
{
  // compute score
  T score = compute_map_score(cols, rows, map, mask);

  // divided by norm scalar
  if(norm_scalar < 1e-5) return 0.0;
  return (score / sqrt((T)norm_scalar));
}

 
template<class T>
inline
T
compute_detector_score( int cols, int rows, T* map, T norm_scalar, int* mask, T* wspacing )
{
  // compute score
  T score = compute_map_score(cols, rows, map, mask, wspacing);

  // divided by norm scalar
  if(norm_scalar < 1e-5) return 0.0;
  return (score / sqrt((T)norm_scalar));
}


template<class T>
T 
compute_map_norm( int cols, int rows, std::vector<T*> &map, int* mask )
{
  T norm_scalar = 0;

  #pragma omp parallel for reduction(+:norm_scalar)
  for(int i=mask[0]; i<=mask[2]; ++i){
    for(int j=mask[1]; j<=mask[3]; ++j){
      T scalar_ij = compute_map_score(cols, rows, map[i*rows + j], mask);
      norm_scalar += scalar_ij;
    }
  }

  return norm_scalar;
}


template<class T>
T 
compute_map_norm( int cols, int rows, std::vector<T*> &map, int* mask, T* wspacing )
{
  T norm_scalar = 0;

  if(fabs(wspacing[0]) < DIFF_EPSILON && fabs(wspacing[1]) < DIFF_EPSILON
    &&  fabs(wspacing[2]) < DIFF_EPSILON && fabs(wspacing[3]) < DIFF_EPSILON)
  {
    return compute_map_norm(cols, rows, map, mask);;
  }

  #pragma omp parallel for reduction(+:norm_scalar)
  for(int i=mask[0]; i<=mask[2]; ++i){
    for(int j=mask[1]; j<=mask[3]; ++j){
      
      T scalar_i, scalar_j;
      scalar_i = scalar_j = 1.0;

      if(i == mask[0])
        scalar_i -= wspacing[0];
      if(i == mask[2])
        scalar_i -= wspacing[2];
      if(j == mask[1])
        scalar_j -= wspacing[1];
      if(j == mask[3])
        scalar_j -= wspacing[3];

      T scalar_ij = compute_map_score(cols, rows, map[i*rows + j], mask, wspacing);
      norm_scalar += scalar_i * scalar_j * scalar_ij;
    }
  }

  return norm_scalar;
}


template<class T>
inline
T 
compute_map_score( int cols, int rows, T* map, int* mask )
{
  T score = 0;

  score += map[mask[2]*rows + mask[3]];
  score -= (mask[0] ? map[(mask[0]-1)*rows + mask[3]] : 0);
  score -= (mask[1] ? map[mask[2]*rows + (mask[1]-1)] : 0);
  score += ((mask[0] && mask[1]) ? map[(mask[0]-1)*rows + (mask[1]-1)] : 0);

  return score;
}


template<class T>
T 
compute_map_score( int cols, int rows, T* map, int* mask, T* wspacing )
{
  T score = 0;
  int wmask[4];
  
  score = compute_map_score(cols, rows, map, mask);

  if(fabs(wspacing[0]) < DIFF_EPSILON && fabs(wspacing[1]) < DIFF_EPSILON
	  &&  fabs(wspacing[2]) < DIFF_EPSILON && fabs(wspacing[3]) < DIFF_EPSILON)
  {
    return score;
  }
  
  // minus scores
  wmask[0] = mask[0];
  wmask[1] = mask[1];
  wmask[2] = mask[0];
  wmask[3] = mask[3];
  score -=  wspacing[0] * compute_map_score(cols, rows, map, wmask);

  wmask[0] = mask[0];
  wmask[1] = mask[1];
  wmask[2] = mask[2];
  wmask[3] = mask[1];
  score -=  wspacing[1] * compute_map_score(cols, rows, map, wmask);

  wmask[0] = mask[2];
  wmask[1] = mask[1];
  wmask[2] = mask[2];
  wmask[3] = mask[3];
  score -=  wspacing[2] * compute_map_score(cols, rows, map, wmask);

  wmask[0] = mask[0];
  wmask[1] = mask[3];
  wmask[2] = mask[2];
  wmask[3] = mask[3];
  score -=  wspacing[3] * compute_map_score(cols, rows, map, wmask);

  // plus scores
  wmask[0] = mask[0];
  wmask[1] = mask[1];
  wmask[2] = mask[0];
  wmask[3] = mask[1];
  score +=  wspacing[0] * wspacing[1] * compute_map_score(cols, rows, map, wmask);

  wmask[0] = mask[2];
  wmask[1] = mask[1];
  wmask[2] = mask[2];
  wmask[3] = mask[1];
  score +=  wspacing[1] * wspacing[2] * compute_map_score(cols, rows, map, wmask);

  wmask[0] = mask[2];
  wmask[1] = mask[3];
  wmask[2] = mask[2];
  wmask[3] = mask[3];
  score +=  wspacing[2] * wspacing[3] * compute_map_score(cols, rows, map, wmask);

  wmask[0] = mask[0];
  wmask[1] = mask[3];
  wmask[2] = mask[0];
  wmask[3] = mask[3];
  score +=  wspacing[0] * wspacing[3] * compute_map_score(cols, rows, map, wmask);

  return score;
}


template<class T>
void 
compute_fisher_map( fisher_net<T> *fnet, std::vector<T*> &x, std::vector<int*> &lc, int* im_size,
                  int* cell_size, int npyrs, std::vector<T*> &models, std::vector<T*> &fmap, std::vector<T*> &fnorm )
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
  fnet->compute_fisher_net(x, idx, fisher);

  // computer fisher norms
  vector_mult(fisher, fkdim, fnorm[0]);

  // compute integral norm maps
  for(int i=0; i<ncells; ++i){
	compute_integral_map(cols, rows, fnorm[i]);
  }
  
  // compute fisher maps with spatial pyramids
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
		fmap_k_i[j] = simd::dot_product(fkdim, w_k_i, fisher[j]);
      }

      // compute integral score map
      compute_integral_map(cols, rows, fmap_k_i);
    }	
  }

  delete[] fisher_pr;
  delete[] idx;
}


template<class T>
void 
compute_fisher_map( fisher_net<T> *fnet, std::vector<T*> &x, std::vector<int*> &lc, int* im_size,
                  int* cell_size1, int* cell_size2, int npyrs, std::vector<T*> &models, std::vector<T*> &fmap, std::vector<T*> &fnorm )
{
  // compute cell assign per descriptor
  int cols1 = (int)ceil(im_size[0] / (float)cell_size1[0]);
  int rows1 = (int)ceil(im_size[1] / (float)cell_size1[1]);
  int cols2 = (int)ceil(im_size[0] / (float)cell_size2[0]);
  int rows2 = (int)ceil(im_size[1] / (float)cell_size2[1]);
  int ndescrs = (int)x.size();
  int* idx = new int[ndescrs];
  int ncells1 = compute_node_assign(lc, im_size[0], im_size[1], cell_size1[0], cell_size1[1], idx);
  int ncells2 = rows2 * cols2;

  // computer fisher vectors
  int fkdim = fnet->dim();
  T* fisher_pr1 = new T[fkdim * ncells1];
  std::vector<T*> fisher1(ncells1);
  for(int i=0; i<ncells1; ++i){
    fisher1[i] = fisher_pr1 + (i * fkdim);
  }
  fnet->compute_fisher_net(x, idx, fisher1);

  // compute fisher maps with spatial pyramids
  int nmodels = (int)models.size();

  #pragma omp parallel for
  for(int k=0; k<nmodels; ++k)
  {
    T* model_k = models[k];
    T* fmap_k = fmap[k];

    for(int i=0; i<npyrs; ++i)
    {
      T* w_k_i = model_k + i*fkdim;
      T* fmap_k_i = fmap_k + i*ncells1;

      // compute cell-wise contribution to classification
      for(int j=0; j<ncells1; ++j)
      {
        fmap_k_i[j] = simd::dot_product(fkdim, w_k_i, fisher1[j]);
      }

      // compute integral score map
      compute_integral_map(cols1, rows1, fmap_k_i);
    }	
  }

  //
  T* fisher_pr2 = new T[fkdim * ncells2];
  std::vector<T*> fisher2(ncells2);
  for(int i=0; i<ncells2; ++i){
    fisher2[i] = fisher_pr2 + (i * fkdim);
  }

  if(cell_size2[0]%cell_size1[0]==0 && cell_size2[1]%cell_size1[1]==0)
  {
    int aggreg_factor[2];
    aggreg_factor[0] = (int)(cell_size2[0]/cell_size1[0]);
    aggreg_factor[1] = (int)(cell_size2[1]/cell_size1[1]);

    aggreg_fisher(cols1, rows1, fkdim, fisher1, aggreg_factor, fisher2);

    // computer fisher norms
    vector_mult(fisher2, fkdim, fnorm[0]);

    // compute integral norm maps
    for(int i=0; i<ncells2; ++i){
      compute_integral_map(cols2, rows2, fnorm[i]);
    }

  }else{
    
    compute_node_assign(lc, im_size[0], im_size[1], cell_size2[0], cell_size2[1], idx);
    fnet->compute_fisher_net(x, idx, fisher2);

    // computer fisher norms
    vector_mult(fisher2, fkdim, fnorm[0]);

    // compute integral norm maps
    for(int i=0; i<ncells2; ++i){
      compute_integral_map(cols2, rows2, fnorm[i]);
    }
 
  }

  delete[] fisher_pr1;
  delete[] fisher_pr2;
  delete[] idx;
}


template<class T>
void
aggreg_fisher( int cols, int rows, int fkdim, std::vector<T*> &fisher, int* agg_factor, std::vector<T*> &agg_fisher )
{
  int ss = 0;
  memset(agg_fisher[0], 0, fkdim*agg_fisher.size()*sizeof(T));

  for(int i=0; i<cols; i=i+agg_factor[0]){
    for(int j=0; j<rows; j=j+agg_factor[1]){

      int ncols = std::min(cols-i, agg_factor[0]);
      int nrows = std::min(rows-j, agg_factor[1]);
      for(int p=0; p<ncols; ++p){
        for(int q=0; q<nrows; ++q){
          simd::add(fkdim, agg_fisher[ss], fisher[(i+p)*rows+(j+q)]);
        }
      }

      ss++;
    }
  }

}


template<class T>
void
compute_integral_map( int cols, int rows, T* map )
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
compute_node_assign( std::vector<int*> &lc, int width, int height,
                      int cell_width, int cell_height, int* idx )
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
vector_mult( std::vector<T*> &vec, int dim, T* result )
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


template int fisher_approx_detector_compute( fisher_net<float> *fnet, std::vector<float*> &x, std::vector<int*> &lc, int* im_size, std::vector<int*> &sppyr,
                           std::vector<float*>& models, std::vector<int*> &boxes, std::vector<float*> &scores, int* cell_size );

template int fisher_approx_detector_compute( fisher_net<double> *fnet, std::vector<double*> &x, std::vector<int*> &lc, int* im_size, std::vector<int*> &sppyr,
                           std::vector<double*>& models, std::vector<int*> &boxes, std::vector<double*> &scores, int* cell_size );

template int fisher_approx_detector_compute( fisher_net<float> *fnet, std::vector<float*> &x, std::vector<int*> &lc, int* im_size, std::vector<int*> &sppyr,
                           std::vector<float*>& models, std::vector<int*> &boxes, std::vector<float*> &scores, int* cell_size1, int* cell_size2 );

template int fisher_approx_detector_compute( fisher_net<double> *fnet, std::vector<double*> &x, std::vector<int*> &lc, int* im_size, std::vector<int*> &sppyr,
                           std::vector<double*>& models, std::vector<int*> &boxes, std::vector<double*> &scores, int* cell_size1, int* cell_size2 );

template void vector_mult<float>( std::vector<float*> &vec, int dim, float* result );
template void vector_mult<double>( std::vector<double*> &vec, int dim, double* result );
