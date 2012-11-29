
/// \Class    fisher-net fisher_net.h "fisher_net.h"
///  
/// \brief
/// \


#ifndef __FISHER_NET_H
#define __FISHER_NET_H

#include <limits>

#include "simd_math.h"
#include "gmm.h"


// -------------------------------------------------------------------------
// Fisher Vector

struct fisher_net_param {
  fisher_net_param() :
    grad_weights(false), 
    grad_means(true), 
    grad_variances(true),
    alpha(0.5), 
    pnorm(2.0) { }
  bool grad_weights;
  bool grad_means;
  bool grad_variances;
  float alpha;
  float pnorm;
  float gamma_eps;
  void print();
};

// -------------------------------------------------------------------------

template<class T>
class fisher_net
{

public:
  
  fisher_net( fisher_net_param &_param );
  ~fisher_net( );
  
  void set_model( gaussian_mixture<T> &_gmm );

  // unweighted
  int compute_fisher( std::vector<T*> &x, T *fk );

  // weighted
  int compute_fisher( std::vector<T*> &x, std::vector<T> &wgh, T *fk);
  
  ////////////////////////////////////////////////////
  // single node
  int compute_fisher_net( std::vector<T*> &x, T *fx );

  // with node index
  int compute_fisher_net( std::vector<T*> &x, int *idx, std::vector<T*> &fx );
  
  // node-wised computation
  int compute_fisher_node( std::vector<T*> &x, T *fn );
  ///////////////////////////////////////////////////

  int dim(){ return fkdim; }

private:

  bool equal( T a, T b )
  { 
    if( fabs((T)a-(T)b)<std::numeric_limits<T>::epsilon() )
      return true;
    return false;
  }

  void alpha_and_lp_normalization( T *fk );

protected:

  fisher_net_param param;

  int ndim, ngauss, ngrad, fkdim;

  gaussian_mixture<T> *gmm;

  T *iwgh, *istd;
};

#endif
