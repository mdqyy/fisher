
#include "fisher_net.h"

template<class T>
fisher_net<T>::fisher_net( fisher_net_param &_param )
  : param(_param), gmm(0), iwgh(0), istd(0)
{
  ngrad = (int)param.grad_weights + (int)param.grad_means + (int)param.grad_variances;
  assert( (param.alpha>0.0) && (param.alpha<=1.0) ); 
}

template<class T>
fisher_net<T>::~fisher_net()
{
  gmm=0;

  delete[] iwgh;
  iwgh=0;

  delete[] istd;
  istd = 0;
}

template<class T>
void
fisher_net<T>::set_model( gaussian_mixture<T> &_gmm )
{
  gmm = &_gmm;
  ndim = gmm->n_dim();
  ngauss = gmm->n_gauss();

  fkdim = 0;
  if( param.grad_weights )
  {
    fkdim += ngauss;
  }
  if( param.grad_means )
  {
    fkdim += ngauss*ndim;
  }
  if( param.grad_variances )
  {
    fkdim += ngauss*ndim;
  }

  delete[] iwgh;

  // precompute inverse weights
  iwgh = new T[ngauss];
  for( int j=0; j<ngauss; ++j )
  {
    assert( gmm->coef[j]>0.0 );
    iwgh[j] = 1.0/gmm->coef[j];
  } 

  // precompute inverse standard deviations
  if( param.grad_means || param.grad_variances )
  {
    delete[] istd;
    istd = new T[ngauss*ndim];

    for( int j=0; j<ngauss; ++j ) 
    {
      T *var_j = gmm->var[j];
      T *istd_j = istd+j*ndim;
      for( int k=ndim; k--; ) 
      {
        assert( var_j[k]>0.0 );
        istd_j[k] = (T)1.0/(T)sqrtf( (float)var_j[k] );
      }
    }    
  }
}


template<class T>
int
fisher_net<T>::compute_fisher( std::vector<T*> &x, T *fk )
{
  std::vector<T> wghx( x.size(), 1.0 );  
  return compute_fisher( x, wghx, fk );
}


template<class T>
int 
fisher_net<T>::compute_fisher( std::vector<T*> &x, std::vector<T> &wghx, T *fk )
{  

  assert(gmm);

  assert( x.size()==wghx.size() );

  int nsamples = (int)wghx.size();

  T wghsum=0.0;
#pragma omp parallel for reduction(+:wghsum)
  for( int i=0; i<nsamples; ++i )
  {
    wghsum += wghx[i];
  }
  assert( wghsum>0 );

  // accumulate statistics
  gmm->reset_stat_acc();
  for( int i=0; i<nsamples; ++i )
  {
    gmm->accumulate_statistics( x[i], true, param.grad_means||param.grad_variances, param.grad_variances );
  }

  T *p=fk;

  // Gradient w.r.t. the mixing weights
  // without the constraint \sum_i pi_i=1 => Soft-BoV
  if( param.grad_weights )
  {
    for( int j=ngauss; j--; )
    {        
      p[j] = gmm->s0[j] / ( wghsum*(T)sqrtf((float)iwgh[j]) );
	  //p[j] = gmm->s0[j] / ( (T)sqrtf((float)iwgh[j]) );
    } 
    p += ngauss;
  }

  // Gradient w.r.t. the means
  if( param.grad_means )
  {
#pragma omp parallel for
    for( int j=0; j<ngauss; j++ )
    {
      T *s1_j = gmm->s1[j];
      T *mean_j = gmm->mean[j];
      T *istd_j = istd+j*ndim;
      T *p_j = p+j*ndim;
      T mc = (T)sqrtf((float)iwgh[j])/wghsum;
	  //T mc = (T)sqrtf((float)iwgh[j]);

      for( int k=ndim; k--; )
      {
        p_j[k] = mc * ( s1_j[k] - mean_j[k] * gmm->s0[j] ) * istd_j[k];
      }      
    }
    p += ngauss*ndim; 
  }

  // Gradient w.r.t. the variances
  if( param.grad_variances )
  {
#pragma omp parallel for
    for( int j=0; j<ngauss; j++ )
    {
      T *s1_j = gmm->s1[j];
      T *s2_j = gmm->s2[j];
      T *mean_j = gmm->mean[j];
      T *var_j = gmm->var[j];
      T *p_j = p+j*ndim;
      T vc = (T)sqrtf(0.5f*(float)iwgh[j])/wghsum;
	  //T vc = (T)sqrtf(0.5f*(float)iwgh[j]);

      for( int k=ndim; k--; )
      {
        p_j[k] = vc * ( ( s2_j[k] + mean_j[k] * ( mean_j[k]*gmm->s0[j] - (T)2.0*s1_j[k] ) ) / var_j[k] - gmm->s0[j] );
      }   
    }
  } 
  
  alpha_and_lp_normalization(fk);
  
  return 0;
}


////////////////////////////////////////////////////
template<class T>
int 
fisher_net<T>::compute_fisher_net( std::vector<T*> &x, T *fx )
{
  return compute_fisher_node(x, fx);
}


template<class T>
int 
fisher_net<T>::compute_fisher_net( std::vector<T*> &x, int *idx, std::vector<T*> &fx )
{
  
  int nsamples = (int)x.size();
  int nnodes = (int)fx.size();

#pragma omp parallel for
  for( int k=0; k<nnodes; ++k )
  {
	std::vector<T*> xx;

    for( int i=0; i<nsamples; ++i )
    {
	  if(idx[i] == k)
	    xx.push_back(x[i]);
    }

	compute_fisher_node(xx, fx[k]);
  }
  
  return 0;
}


template<class T>
int 
fisher_net<T>::compute_fisher_node( std::vector<T*> &x, T *fn )
{  

  assert(gmm);

  int nsamples = (int)x.size();

  // accumulate statistics
  gmm->reset_stat_acc();
  for( int i=0; i<nsamples; ++i )
  {
    gmm->accumulate_statistics( x[i], true, param.grad_means||param.grad_variances, param.grad_variances );
  }

  T *p=fn;

  // Gradient w.r.t. the mixing weights
  // without the constraint \sum_i pi_i=1 => Soft-BoV
  if( param.grad_weights )
  {
    for( int j=ngauss; j--; )
    {
	  p[j] = gmm->s0[j] / ( (T)sqrtf((float)iwgh[j]) );
    } 
    p += ngauss;
  }

  // Gradient w.r.t. the means
  if( param.grad_means )
  {
#pragma omp parallel for
    for( int j=0; j<ngauss; j++ )
    {
      T *s1_j = gmm->s1[j];
      T *mean_j = gmm->mean[j];
      T *istd_j = istd+j*ndim;
      T *p_j = p+j*ndim;
	  T mc = (T)sqrtf((float)iwgh[j]);

      for( int k=ndim; k--; )
      {
        p_j[k] = mc * ( s1_j[k] - mean_j[k] * gmm->s0[j] ) * istd_j[k];
      }      
    }
    p += ngauss*ndim; 
  }

  // Gradient w.r.t. the variances
  if( param.grad_variances )
  {
#pragma omp parallel for
    for( int j=0; j<ngauss; j++ )
    {
      T *s1_j = gmm->s1[j];
      T *s2_j = gmm->s2[j];
      T *mean_j = gmm->mean[j];
      T *var_j = gmm->var[j];
      T *p_j = p+j*ndim;
	  T vc = (T)sqrtf(0.5f*(float)iwgh[j]);

      for( int k=ndim; k--; )
      {
        p_j[k] = vc * ( ( s2_j[k] + mean_j[k] * ( mean_j[k]*gmm->s0[j] - (T)2.0*s1_j[k] ) ) / var_j[k] - gmm->s0[j] );
      }   
    }
  }
  
  return 0;
}
////////////////////////////////////////////////


template<class T>
void
fisher_net<T>::alpha_and_lp_normalization( T *fk )
{
  // alpha normalization
  if( !equal(param.alpha,1.0f) )
  {
    if( equal(param.alpha,0.5f) )
    {
#pragma omp parallel for
      for( int i=0; i<fkdim; i++ )
      {
        if( fk[i]<0.0 )
          fk[i] = -std::sqrt(-fk[i]);
        else
          fk[i] = std::sqrt(fk[i]);
      }
    }
    else
    {
#pragma omp parallel for
      for( int i=0; i<fkdim; i++ )
      {
        if( fk[i]<0.0 )
          fk[i] = -std::pow(-fk[i],(T)param.alpha);
        else
          fk[i] = std::pow(fk[i],(T)param.alpha);
      }
    }
  }

  // Lp normalization
  if( !equal(param.pnorm,(float)0.0) )
  {
    T pnorm=0;
    if( equal(param.pnorm,(float)1.0) )
    {
#pragma omp parallel for reduction(+:pnorm)
      for( int i=0; i<fkdim; ++i )
      {
        pnorm += std::fabs(fk[i]);
      }
    }
    else if( equal(param.pnorm,2.0) )
    {
      pnorm = sqrt( simd::dot_product( fkdim, fk, fk ) );
    }
    else
    {
#pragma omp parallel for reduction(+:pnorm)
      for( int i=0; i<fkdim; ++i )
      {
        pnorm += std::pow( std::fabs(fk[i]), (T)param.pnorm );
      }
      //pnorm = std::pow( pnorm, 1.0/(T)param.pnorm );
      pnorm = std::pow( pnorm, (T)(1.0/param.pnorm) );
    }

    if( pnorm>0.0 )
    {
      simd::scale( fkdim, fk, (T)(1.0/pnorm) );
    }
  }
}

/// \bief print
/// 
/// \param none
///
/// \return none
///
/// \author Jorge Sanchez
/// \date    August 2009

void
fisher_net_param::print()
{
  std::cout << "  grad_weights = " << grad_weights << std::endl;
  std::cout << "  grad_means = " << grad_means << std::endl;
  std::cout << "  grad_variances = " << grad_variances << std::endl;
  std::cout << "  alpha = " << alpha << std::endl;
  std::cout << "  pnorm = " << pnorm << std::endl;
}

template class fisher_net<float>;
template class fisher_net<double>;
