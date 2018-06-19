#include "../aflow.h"
#include "aflow_apl.h"

namespace apl {

#ifndef QH_OP_H
#define QH_OP_H

#define _iszero(a) (std::abs(a) < _AFLOW_APL_EPS_) ? true : false
#define _isequal(a, b) (std::abs(a - b) < _AFLOW_APL_EPS_) ? true : false

#define APL_DBL_EPSILON 2.2204460492503131e-16
#define APL_DBL_MIN 2.2250738585072014e-308
#undef APL_RANGE_CHECK
#define APL_RANGE_CHECK 1
#undef APL_RANGE_COND
#define NULL_VECTOR_VIEW \
  {                      \
    { 0, 0, 0, 0, 0 }    \
  }
#define NULL_VECTOR \
  { 0, 0, 0, 0, 0 }
#define NULL_MATRIX \
  { 0, 0, 0, 0, 0, 0 }
#define NULL_MATRIX_VIEW \
  {                      \
    { 0, 0, 0, 0, 0, 0 } \
  }
#define RETURN_IF_NULL(x) \
  if (!x) { return; }
#define APL_SQRT_DBL_MAX 1.3407807929942596e+154
#define APL_SQRT_DBL_MIN 1.4916681462400413e-154
#define APL_SQRT_DBL_EPSILON 1.4901161193847656e-08
#define APL_MULTIFIT_FN_EVAL_F(F, x, y) ((*((F)->f))(x, (F)->params, (y)))
#define APL_MULTIFIT_FN_EVAL_DF(F, x, dy) ((*((F)->df))(x, (F)->params, (dy)))
#define APL_MULTIFIT_FN_EVAL_F_DF(F, x, y, dy) ((*((F)->fdf))(x, (F)->params, (y), (dy)))
#define CONST_REAL(a, i) (((const double *)a)[2 * (i)])
#define CONST_IMAG(a, i) (((const double *)a)[2 * (i) + 1])
#define CONST_REAL0(a) (((const double *)a)[0])
#define CONST_IMAG0(a) (((const double *)a)[1])
#define REAL0(a) (((double *)a)[0])
#define IMAG0(a) (((double *)a)[1])
#define REAL(a, i) (((double *)a)[2 * (i)])
#define IMAG(a, i) (((double *)a)[2 * (i) + 1])
#define APL_MV_COMPLEX_P(zp) ((zp)->dat)
#define APL_MV_COMPLEX_P_REAL(zp) ((zp)->dat[0])
#define APL_MV_COMPLEX_P_IMAG(zp) ((zp)->dat[1])
#define APL_MV_COMPLEX_EQ(z1, z2) (((z1).dat[0] == (z2).dat[0]) && ((z1).dat[1] == (z2).dat[1]))

typedef std::complex<double> _CD_;
typedef vector<_CD_> _CVEC_;
typedef vector<_CVEC_> _CMAT_;
typedef vector<double> _VEC_;
typedef vector<_VEC_> _MAT_;

typedef enum {
  APL_MV_EIGEN_SORT_VAL_ASC,
  APL_MV_EIGEN_SORT_VAL_DESC,
  APL_MV_EIGEN_SORT_ABS_ASC,
  APL_MV_EIGEN_SORT_ABS_DESC
} apl_eigen_sort_t;

enum APL_ORDER { aplRowMajor = 101,
                 aplColMajor = 102 };
enum APL_TRANSPOSE { aplNoTrans = 111,
                     aplTrans = 112,
                     aplConjTrans = 113 };
enum APL_UPLO { aplUpper = 121,
                aplLower = 122 };

typedef enum APL_TRANSPOSE APL_TRANSPOSE_t;
typedef enum APL_UPLO APL_UPLO_t;  //when UPLO is ‘U’ the upper triangle (‘L’ for lower triangle)

template <class T>
struct apl_complex {
  T dat[2];
};

#define APL_MV_REAL(z) ((z).dat[0])
#define APL_MV_IMAG(z) ((z).dat[1])
#define APL_MV_SET_COMPLEX(zp, x, y) \
  do {                               \
    (zp)->dat[0] = (x);              \
    (zp)->dat[1] = (y);              \
  } while (0)
#define APL_MV_COMPLEX_AT(zv, i) ((apl_complex<double> *)&((zv)->data[2 * (i) * (zv)->stride]))
#define APL_MV_COMPLEX_ONE (apl_complex_rect(1.0, 0.0))
#define APL_MV_COMPLEX_ZERO (apl_complex_rect(0.0, 0.0))

template <class T>
struct apl_block {
  size_t size;
  T *data;
};

template <class T>
struct apl_vector {
  size_t size;
  size_t stride;
  T *data;
  apl_block<T> *block;
  int owner;
};

//TEST
template <class T>
struct apl_block_complex {
  size_t size;
  T *data;
};

template <class T>
struct apl_vector_complex {
  size_t size;
  size_t stride;
  double *data;
  apl_block_complex<T> *block;
  int owner;
};

template <class T>
struct apl_vector_complex_view {
  apl_vector_complex<T> vector;
};

template <class T>
struct apl_matrix_complex {
  size_t size1;
  size_t size2;
  size_t tda;
  T *data;
  apl_block_complex<T> *block;
  int owner;
};

template <class T>
struct apl_matrix_complex_view {
  apl_matrix_complex<T> matrix;
};

typedef struct {
  size_t size;
  double *d;
  double *sd;
  double *tau;
  double *gc;
  double *gs;
} apl_eigen_hermv_workspace;

//TEST

template <class T>
struct apl_vector_view {
  apl_vector<T> vector;
};

template <class T>
struct apl_matrix {
  size_t size1;
  size_t size2;
  size_t tda;
  T *data;
  apl_block<T> *block;
  int owner;
};

template <class T>
struct apl_matrix_view {
  apl_matrix<T> matrix;
};

struct apl_matrix_const_view {
  apl_matrix<double> matrix;
};

//typedef const apl_matrix_const_view apl_matrix_const_view;

template <class T>
struct apl_multifit_linear_workspace {
  size_t n;  // number of observations //
  size_t p;  // number of parameters //
  apl_matrix<T> *A;
  apl_matrix<T> *Q;
  apl_matrix<T> *QSI;
  apl_vector<T> *S;
  apl_vector<T> *t;
  apl_vector<T> *xt;
  apl_vector<T> *D;
};

// Definition of vector-valued functions with parameters based on apl_vector //
template <class T>
struct apl_multifit_function_struct {
  int (*f)(const apl_vector<T> *x, void *params, apl_vector<T> *f);
  size_t n; /* number of functions */
  size_t p; /* number of independent variables */
  void *params;
};

typedef struct apl_multifit_function_struct<double> apl_multifit_function;

template <class T>
class apl_multifit_function_fdf_struct {
 public:
  int (*f)(const apl_vector<T> *x, void *params, apl_vector<T> *f);
  int (*df)(const apl_vector<T> *x, void *params, apl_matrix<T> *df);
  int (*fdf)(const apl_vector<T> *x, void *params, apl_vector<T> *f, apl_matrix<T> *df);
  size_t n; /* number of functions */
  size_t p; /* number of independent variables */
  void *params;
};

typedef class apl_multifit_function_fdf_struct<double> apl_multifit_function_fdf;

#define APL_MULTIFIT_FN_EVAL(F, x, y) (*((F)->f))(x, (F)->params, (y))

template <class T>
struct apl_multifit_fsolver_type {
  const char *name;
  size_t size;
  int (*alloc)(void *state, size_t n, size_t p);
  int (*set)(void *state, apl_multifit_function *function, apl_vector<T> *x, apl_vector<T> *f, apl_vector<T> *dx);
  int (*iterate)(void *state, apl_multifit_function *function, apl_vector<T> *x, apl_vector<T> *f, apl_vector<T> *dx);
  void (*free)(void *state);
};

template <class T>
struct apl_multifit_fsolver {
  const apl_multifit_fsolver_type<T> *type;
  apl_multifit_function *function;
  apl_vector<T> *x;
  apl_vector<T> *f;
  apl_vector<T> *dx;
  void *state;
};

template <class T>
struct apl_multifit_fdfsolver_type {
  const char *name;
  size_t size;
  int (*alloc)(void *state, size_t n, size_t p);
  int (*set)(void *state, apl_multifit_function_fdf *fdf, apl_vector<T> *x, apl_vector<T> *f, apl_matrix<T> *J, apl_vector<T> *dx);
  int (*iterate)(void *state, apl_multifit_function_fdf *fdf, apl_vector<T> *x, apl_vector<T> *f, apl_matrix<T> *J, apl_vector<T> *dx);
  void (*free)(void *state);
};

//const apl_multifit_fdfsolver_type<double> * apl_multifit_fdfsolver_lmsder;

template <class T>
struct apl_multifit_fdfsolver {
  const apl_multifit_fdfsolver_type<T> *type;
  apl_multifit_function_fdf *fdf;
  apl_vector<T> *x;
  apl_vector<T> *f;
  apl_matrix<T> *J;
  apl_vector<T> *dx;
  void *state;
};

struct apl_permutation_struct {
  size_t size;
  size_t *data;
};

typedef struct apl_permutation_struct apl_permutation;

template <class utype>
struct lmder_state_t {
  size_t iter;
  double xnorm;
  double fnorm;
  double delta;
  double par;
  apl_matrix<utype> *r;
  apl_vector<utype> *tau;
  apl_vector<utype> *diag;
  apl_vector<utype> *qtf;
  apl_vector<utype> *newton;
  apl_vector<utype> *gradient;
  apl_vector<utype> *x_trial;
  apl_vector<utype> *f_trial;
  apl_vector<utype> *df;
  apl_vector<utype> *sdiag;
  apl_vector<utype> *rptdx;
  apl_vector<utype> *w;
  apl_vector<utype> *work1;
  apl_permutation *perm;
};

class MVops {
 public:
  MVops() {}
  ~MVops() {}
  void clear() { } //this->clear(); }

 public:
  //vector function declarations
  template <class T>
  apl_block<T> *apl_block_alloc(const size_t n);
  template <class T>
  apl_vector<T> *apl_vector_alloc(const size_t n);
  template <class T>
  apl_vector<T> *apl_vector_calloc(const size_t n);
  template <class T>
  void apl_vector_free(apl_vector<T> *v);
  template <class T>
  void apl_block_free(apl_block<T> *b);
  template <class T>
  inline void apl_vector_set(apl_vector<T> *v, const size_t i, T x);
  template <class T>
  inline T apl_vector_get(const apl_vector<T> *v, const size_t i);
  template <class T>
  apl_vector_view<T> apl_vector_subvector(apl_vector<T> *v, size_t offset, size_t n);
  template <class T>
  int apl_vector_scale(apl_vector<T> *a, const double x);
  template <class T>
  void
  apl_vector_set_all(apl_vector<T> *v, T x);
  template <class T>
  void
  apl_vector_set_zero(apl_vector<T> *v);
  template <class T>
  int apl_vector_div(apl_vector<T> *a, const apl_vector<T> *b);
  template <class T>
  int apl_vector_memcpy(apl_vector<T> *dest,
                        const apl_vector<T> *src);
  template <class T>
  const apl_vector_view<T> apl_vector_const_subvector(const apl_vector<T> *v, size_t offset, size_t n);
  template <class T>
  int apl_vector_swap_elements(apl_vector<T> *v, const size_t i, const size_t j);
  template <class T>
  apl_vector_view<T>
  apl_vector_view_array(T *base, size_t n);
  template <class utype>
  static size_t
  count_nsing(const apl_matrix<utype> *r);

  //matrix function declarations
  template <class T>
  apl_matrix<T> *
  apl_matrix_alloc(const size_t n1, const size_t n2);
  template <class T>
  void apl_matrix_free(apl_matrix<T> *m);

  template <class T>
  inline T apl_matrix_get(const apl_matrix<T> *m, const size_t i, const size_t j);
  template <class T>
  inline void apl_matrix_set(apl_matrix<T> *m, const size_t i, const size_t j, const T x);

  template <class T>
  apl_matrix_view<T>
  apl_matrix_submatrix(apl_matrix<T> *m,
                       const size_t i, const size_t j,
                       const size_t n1, const size_t n2);
  template <class T>
  int apl_matrix_memcpy(apl_matrix<T> *dest,
                        const apl_matrix<T> *src);

  template <class T>
  apl_vector_view<T>
  apl_matrix_row(apl_matrix<T> *m, const size_t i);

  template <class T>
  apl_vector_view<T>
  apl_matrix_column(apl_matrix<T> *m, const size_t j);

  template <class T>
  const apl_vector_view<T>
  apl_matrix_const_column(const apl_matrix<T> *m, const size_t j);

  template <class T>
  void
  apl_matrix_set_identity(apl_matrix<T> *m);
  template <class T>
  int apl_matrix_swap_columns(apl_matrix<T> *m,
                              const size_t i, const size_t j);
  template <class T>
  const apl_vector_view<T>
  apl_matrix_const_row(const apl_matrix<T> *m, const size_t i);
  template <class T>
  apl_matrix<T> *
  apl_matrix_calloc(const size_t n1, const size_t n2);

  //apl fit
  template <class T>
  apl_multifit_linear_workspace<T> *
  apl_multifit_linear_alloc(size_t n, size_t p);
  template <class T>
  void
  apl_multifit_linear_free(apl_multifit_linear_workspace<T> *work);
  template <class T>
  int apl_linalg_balance_columns(apl_matrix<T> *A, apl_vector<T> *D);
  template <class T>
  T apl_blas_dasum(const apl_vector<T> *v);
  template <class T>
  void
  apl_blas_dscal(T alpha, apl_vector<T> *v);
  template <class T>
  T apl_blas_dnrm2(const apl_vector<T> *v);
  template <class T>
  size_t
  apl_blas_idamax(const apl_vector<T> *X);
  size_t
  cblas_idamax(const int N, const double *X, const int incX);

  template <class T>
  int apl_linalg_bidiag_decomp(apl_matrix<T> *A, apl_vector<T> *tau_U, apl_vector<T> *tau_V);

  template <class T>
  int apl_linalg_bidiag_unpack2(apl_matrix<T> *A,
                                apl_vector<T> *tau_U,
                                apl_vector<T> *tau_V,
                                apl_matrix<T> *V);

  template <class T>
  void
  svd2(apl_vector<T> *d, apl_vector<T> *f, apl_matrix<T> *U, apl_matrix<T> *V);

  template <class T>
  void
  qrstep(apl_vector<T> *d, apl_vector<T> *f, apl_matrix<T> *U, apl_matrix<T> *V);

  template <class T>
  int apl_linalg_SV_decomp(apl_matrix<T> *A, apl_matrix<T> *V, apl_vector<T> *S,
                           apl_vector<T> *work);

  template <class T>
  int apl_linalg_SV_decomp_mod(apl_matrix<T> *A,
                               apl_matrix<T> *X,
                               apl_matrix<T> *V, apl_vector<T> *S, apl_vector<T> *work);

  template <class T>
  static void
  chop_small_elements(apl_vector<T> *d, apl_vector<T> *f);
  inline static void
  create_givens(const double a, const double b, double *c, double *s);
  template <class T>
  static void
  chase_out_trailing_zero(apl_vector<T> *d, apl_vector<T> *f, apl_matrix<T> *V);
  static void
  create_schur(double d0, double f0, double d1, double *c, double *s);
  template <class T>
  static void
  chase_out_intermediate_zero(apl_vector<T> *d, apl_vector<T> *f, apl_matrix<T> *U, size_t k0);
  template <class T>
  static double
  trailing_eigenvalue(const apl_vector<T> *d, const apl_vector<T> *f);
  template <class T>
  int apl_blas_daxpy(double alpha, const apl_vector<T> *u, apl_vector<T> *v);
  void
  cblas_dgemv(const enum APL_ORDER order, const enum APL_TRANSPOSE TransA,
              const int M, const int N, const double alpha, const double *A,
              const int lda, const double *X, const int incX,
              const double beta, double *Y, const int incY);
  template <class T>
  int apl_blas_ddot(const apl_vector<T> *X, const apl_vector<T> *Y, T *result);
  double
  cblas_ddot(const int N, const double *X, const int incX, const double *Y,
             const int incY);
  template <class T>
  int apl_blas_dgemv(APL_TRANSPOSE_t TransA,
                     double alpha,
                     const apl_matrix<T> *A,
                     const apl_vector<T> *X,
                     double beta,
                     apl_vector<T> *Y);
  template <class T>
  int apl_multifit_wlinear(const apl_matrix<T> *X,
                           const apl_vector<T> *w,
                           const apl_vector<T> *y,
                           apl_vector<T> *c,
                           apl_matrix<T> *cov,
                           double *chisq,
                           apl_multifit_linear_workspace<T> *work);
  template <class T>
  T apl_linalg_householder_transform(apl_vector<T> *v);
  template <class T>
  int apl_linalg_householder_mh(double tau, const apl_vector<T> *v, apl_matrix<T> *A);
  template <class T>
  int apl_linalg_householder_hm(T tau, const apl_vector<T> *v, apl_matrix<T> *A);
  template <class T>
  int apl_linalg_householder_hm1(double tau, apl_matrix<T> *A);
  template <class T>
  static int
  multifit_wlinear_svd(const apl_matrix<T> *X,
                       const apl_vector<T> *w,
                       const apl_vector<T> *y,
                       double tol,
                       int balance,
                       size_t *rank,
                       apl_vector<T> *c,
                       apl_matrix<T> *cov,
                       double *chisq, apl_multifit_linear_workspace<T> *work);
  apl_permutation *apl_permutation_alloc(const size_t n);
  apl_permutation *apl_permutation_calloc(const size_t n);
  inline size_t apl_permutation_get(const apl_permutation *p, const size_t i);
  void apl_permutation_free(apl_permutation *p);
  void
  apl_permutation_init(apl_permutation *p);
  int apl_permutation_swap(apl_permutation *p, const size_t i, const size_t j);
  int apl_permute_inverse(const size_t *p, double *data, const size_t stride, const size_t n);
  template <class utype>
  int apl_permute_vector_inverse(const apl_permutation *p, apl_vector<utype> *v);

  static double
  scaled_enorm(const apl_vector<double> *d, const apl_vector<double> *f);
  template <class utype>
  static double
  enorm(const apl_vector<utype> *f);

  template <class utype>
  static void
  compute_gradient_direction(const apl_matrix<utype> *r, const apl_permutation *p,
                             const apl_vector<utype> *qtf, const apl_vector<utype> *diag,
                             apl_vector<utype> *g);

  //nonlinear functions
  template <class utype>
  apl_multifit_fsolver<utype> *
  apl_multifit_fsolver_alloc(const apl_multifit_fsolver_type<utype> *T,
                             size_t n, size_t p);

  template <class utype>
  void apl_multifit_fsolver_free(apl_multifit_fsolver<utype> *s);
  template <class utype>
  int apl_multifit_fsolver_set(apl_multifit_fsolver<utype> *s,
                               apl_multifit_function *f,
                               const apl_vector<utype> *x);
  template <class utype>
  int apl_multifit_fdfsolver_set(apl_multifit_fdfsolver<utype> *s,
                                 apl_multifit_function_fdf *f,
                                 const apl_vector<utype> *x);
  template <class utype>
  int apl_multifit_fdfsolver_iterate(apl_multifit_fdfsolver<utype> *s);
  template <class utype>
  int apl_linalg_QRPT_decomp(apl_matrix<utype> *A,
                             apl_vector<utype> *tau,
                             apl_permutation *p,
                             int *signum,
                             apl_vector<utype> *norm);
  template <class utype>
  static int iterate(void *vstate, apl_multifit_function_fdf *fdf, apl_vector<utype> *x, apl_vector<utype> *f, apl_matrix<utype> *J, apl_vector<utype> *dx, int scale);
  template <class utype>
  int apl_linalg_householder_hv(double tau, const apl_vector<utype> *v, apl_vector<utype> *w);
  template <class utype>
  static int
  lmpar(apl_matrix<utype> *r, const apl_permutation *perm, const apl_vector<utype> *qtf,
        const apl_vector<utype> *diag, double delta, double *par_inout,
        apl_vector<utype> *newton, apl_vector<utype> *gradient, apl_vector<utype> *sdiag,
        apl_vector<utype> *x, apl_vector<utype> *w);
  template <class utype>
  static void
  compute_newton_bound(const apl_matrix<utype> *r, const apl_vector<utype> *x,
                       double dxnorm, const apl_permutation *perm,
                       const apl_vector<utype> *diag, apl_vector<utype> *w);
  template <class utype>
  static int
  qrsolv(apl_matrix<utype> *r, const apl_permutation *p, const double lambda,
         const apl_vector<utype> *diag, const apl_vector<utype> *qtb,
         apl_vector<utype> *x, apl_vector<utype> *sdiag, apl_vector<utype> *wa);
  template <class utype>
  static void
  compute_newton_correction(const apl_matrix<utype> *r, const apl_vector<utype> *sdiag,
                            const apl_permutation *p, apl_vector<utype> *x,
                            double dxnorm,
                            const apl_vector<utype> *diag, apl_vector<utype> *w);
  template <class utype>
  int apl_linalg_QR_QTvec(const apl_matrix<utype> *QR, const apl_vector<utype> *tau, apl_vector<utype> *v);
  template <class utype>
  static void
  compute_trial_step(apl_vector<utype> *x, apl_vector<utype> *dx, apl_vector<utype> *x_trial);
  static double
  compute_actual_reduction(double fnorm, double fnorm1);
  template <class utype>
  static void
  compute_rptdx(const apl_matrix<utype> *r, const apl_permutation *p,
                const apl_vector<utype> *dx, apl_vector<utype> *rptdx);
  template <class utype>
  apl_multifit_fdfsolver<utype> *
  apl_multifit_fdfsolver_alloc(const apl_multifit_fdfsolver_type<utype> *T,
                               size_t n, size_t p);
  template <class utype>
  int apl_multifit_test_delta(const apl_vector<utype> *dx, const apl_vector<utype> *x,
                              double epsabs, double epsrel);
  template <class utype>
  int apl_multifit_covar(const apl_matrix<utype> *J, double epsrel, apl_matrix<utype> *covar);
  template <class utype>
  void
  apl_multifit_fdfsolver_free(apl_multifit_fdfsolver<utype> *s);
  template <class utype>
  static int
  lmder_alloc(void *vstate, size_t n, size_t p);

  template <class utype>
  static int lmsder_set(void *vstate, apl_multifit_function_fdf *fdf, apl_vector<utype> *x, apl_vector<utype> *f, apl_matrix<utype> *J, apl_vector<utype> *dx);
  template <class utype>
  static int lmder_set(void *vstate, apl_multifit_function_fdf *fdf, apl_vector<utype> *x, apl_vector<utype> *f, apl_matrix<utype> *J, apl_vector<utype> *dx);
  template <class utype>
  static int set(void *vstate, apl_multifit_function_fdf *fdf, apl_vector<utype> *x, apl_vector<utype> *f, apl_matrix<utype> *J, apl_vector<utype> *dx, int scale);
  template <class utype>
  static int lmder_iterate(void *vstate, apl_multifit_function_fdf *fdf, apl_vector<utype> *x, apl_vector<utype> *f, apl_matrix<utype> *J, apl_vector<utype> *dx);
  static void lmder_free(void *vstate);
  template <class utype>
  static int
  lmsder_iterate(void *vstate, apl_multifit_function_fdf *fdf, apl_vector<utype> *x, apl_vector<utype> *f, apl_matrix<utype> *J, apl_vector<utype> *dx);

  template <class utype>
  int apl_multifit_fdfsolver_dif_fdf(const apl_vector<utype> *x, apl_multifit_function_fdf *fdf,
                                     apl_vector<utype> *f, apl_matrix<utype> *J);

  template <class utype>
  static int
  fdjac(const apl_vector<utype> *x, apl_multifit_function_fdf *fdf,
        const apl_vector<utype> *f, apl_matrix<utype> *J);
  template <class utype>
  static void
  compute_diag(const apl_matrix<utype> *J, apl_vector<utype> *diag);
  template <class utype>
  static double
  compute_delta(apl_vector<utype> *diag, apl_vector<utype> *x);
  template <class utype>
  static void
  update_diag(const apl_matrix<utype> *J, apl_vector<utype> *diag);

  template <class utype>
  int apl_multifit_fdfsolver_dif_df(const apl_vector<utype> *x, apl_multifit_function_fdf *fdf,
                                    const apl_vector<utype> *f, apl_matrix<utype> *J);
  template <class utype>
  static void
  lmder_free(void *vstate);

  template <class utype>
  static void
  compute_newton_direction(const apl_matrix<utype> *r, const apl_permutation *perm,
                           const apl_vector<utype> *qtf, apl_vector<utype> *x);

  template <class utype>
  const char *
  apl_multifit_fdfsolver_name(const apl_multifit_fdfsolver<utype> *s) {
    return s->type->name;
  }

  //TEST
  template <class T>
  apl_block_complex<T> *apl_block_complex_alloc(const size_t n);
  template <class T>
  apl_vector_complex<T> *apl_vector_complex_alloc(const size_t n);
  template <class T>
  apl_vector_complex<T> *apl_vector_complex_calloc(const size_t n);
  template <class T>
  void apl_vector_complex_free(apl_vector_complex<T> *v);
  template <class T>
  void apl_block_complex_free(apl_block_complex<T> *b);
  template <class T>
  apl_vector_complex_view<T>
  apl_vector_complex_subvector(apl_vector_complex<T> *v, size_t offset, size_t n);
  template <class T>
  inline apl_complex<T> apl_vector_complex_get(const apl_vector_complex<T> *v, const size_t i);
  template <class T>
  inline void apl_vector_complex_set(apl_vector_complex<T> *v, const size_t i, apl_complex<T> z);
  template <class T>
  const apl_vector_complex_view<T>
  apl_vector_complex_const_subvector(const apl_vector_complex<T> *v, size_t offset, size_t n);
  template <class T>
  int apl_vector_complex_scale(apl_vector_complex<T> *a, const apl_complex<T> x);
  template <class T>
  int apl_vector_complex_div(apl_vector_complex<T> *a, const apl_vector_complex<T> *b);
  template <class T>
  void
  apl_vector_complex_set_all(apl_vector_complex<T> *v, T x);
  template <class T>
  int apl_vector_complex_memcpy(apl_vector_complex<T> *dest,
                                const apl_vector_complex<T> *src);
  template <class T>
  int apl_vector_complex_swap_elements(apl_vector_complex<T> *v, const size_t i, const size_t j);
  template <class T>
  apl_vector_complex_view<T>
  apl_vector_complex_view_array(T *base, size_t n);
  template <class T>
  void
  apl_vector_complex_set_zero(apl_vector_complex<T> *v);
  template <class T>
  apl_matrix_complex<T> *
  apl_matrix_complex_alloc(const size_t n1, const size_t n2);
  template <class T>
  void
  apl_matrix_complex_free(apl_matrix_complex<T> *m);
  template <class T>
  inline apl_complex<T>
  apl_matrix_complex_get(const apl_matrix_complex<T> *m, const size_t i, const size_t j);
  template <class T>
  inline void
  apl_matrix_complex_set(apl_matrix_complex<T> *m, const size_t i, const size_t j, const apl_complex<T> x);
  template <class T>
  int apl_matrix_complex_memcpy(apl_matrix_complex<T> *dest,
                                const apl_matrix_complex<T> *src);
  template <class T>
  const apl_vector_complex_view<T>
  apl_matrix_complex_const_row(const apl_matrix_complex<T> *m, const size_t i);
  template <class T>
  apl_vector_complex_view<T>
  apl_matrix_complex_row(apl_matrix_complex<T> *m, const size_t i);
  template <class T>
  apl_vector_complex_view<T>
  apl_matrix_complex_column(apl_matrix_complex<T> *m, const size_t j);
  template <class T>
  const apl_vector_complex_view<T>
  apl_matrix_complex_const_column(const apl_matrix_complex<T> *m, const size_t j);
  template <class T>
  apl_matrix_complex_view<T>
  apl_matrix_complex_submatrix(apl_matrix_complex<T> *m,
                               const size_t i, const size_t j,
                               const size_t n1, const size_t n2);
  template <class T>
  void
  apl_matrix_complex_set_identity(apl_matrix_complex<T> *m);
  template <class T>
  int apl_matrix_complex_swap_columns(apl_matrix_complex<T> *m,
                                      const size_t i, const size_t j);
  template <class T>
  apl_matrix_complex<T> *
  apl_matrix_complex_calloc(const size_t n1, const size_t n2);
  void apl_eigen_hermv_free(apl_eigen_hermv_workspace *w);
  apl_eigen_hermv_workspace *apl_eigen_hermv_alloc(const size_t n);
  template <class T>
  int apl_eigen_hermv(apl_matrix_complex<T> *A, apl_vector<T> *eval,
                      apl_matrix_complex<T> *evec,
                      apl_eigen_hermv_workspace *w);
  inline apl_complex<double> apl_complex_rect(double x, double y);
  template <class T>
  int apl_linalg_hermtd_decomp(apl_matrix_complex<T> *A, apl_vector_complex<T> *tau);
  template <class T>
  double apl_complex_abs(apl_complex<T> z);
  template <class T>
  apl_complex<T>
  apl_linalg_complex_householder_transform(apl_vector_complex<T> *v);
  template <class T>
  double apl_blas_dznrm2(const apl_vector_complex<T> *X);
  double cblas_dznrm2(const int N, const void *X, const int incX);
  template <class T>
  apl_complex<T>
  apl_complex_sub_real(apl_complex<T> a, double x);
  template <class T>
  apl_complex<T>
  apl_complex_inverse(apl_complex<T> a);
  template <class T>
  void
  apl_blas_zscal(const apl_complex<T> alpha, apl_vector_complex<T> *X);
  void cblas_zscal(const int N, const void *alpha, void *X, const int incX);
  template <class T>
  int apl_blas_zhemv(APL_UPLO Uplo, const apl_complex<T> alpha,
                     const apl_matrix_complex<T> *A, const apl_vector_complex<T> *X,
                     const apl_complex<T> beta, apl_vector_complex<T> *Y);
  void
  cblas_zhemv(const enum APL_ORDER order, const enum APL_UPLO Uplo,
              const int N, const void *alpha, const void *A, const int lda,
              const void *X, const int incX, const void *beta, void *Y,
              const int incY);

  template <class T>
  apl_complex<T>
  apl_complex_mul(apl_complex<T> a, apl_complex<T> b);
  template <class T>
  int apl_blas_zdotc(const apl_vector_complex<T> *X, const apl_vector_complex<T> *Y,
                     apl_complex<T> *dotc);
  void
  cblas_zdotc_sub(const int N, const void *X, const int incX, const void *Y,
                  const int incY, void *result);

  template <class T>
  apl_complex<T>
  apl_complex_mul_real(apl_complex<T> a, double x);
  template <class T>
  int apl_blas_zaxpy(const apl_complex<T> alpha, const apl_vector_complex<T> *X,
                     apl_vector_complex<T> *Y);
  void
  cblas_zaxpy(const int N, const void *alpha, const void *X, const int incX,
              void *Y, const int incY);

  template <class T>
  int apl_blas_zher2(APL_UPLO_t Uplo, const apl_complex<T> alpha,
                     const apl_vector_complex<T> *X, const apl_vector_complex<T> *Y,
                     apl_matrix_complex<T> *A);

  void
  cblas_zher2(const enum APL_ORDER order, const enum APL_UPLO Uplo,
              const int N, const void *alpha, const void *X, const int incX,
              const void *Y, const int incY, void *A, const int lda);

  template <class T>
  int apl_linalg_hermtd_unpack(const apl_matrix_complex<T> *A,
                               const apl_vector_complex<T> *tau,
                               apl_matrix_complex<T> *U,
                               apl_vector<T> *diag,
                               apl_vector<T> *sdiag);

  template <class T>
  int apl_linalg_complex_householder_hm(apl_complex<T> tau, const apl_vector_complex<T> *v, apl_matrix_complex<T> *A);

  template <class T>
  apl_complex<T>
  apl_complex_conjugate(apl_complex<T> a);

  template <class T>
  apl_complex<T>
  apl_complex_add(apl_complex<T> a, apl_complex<T> b);
  template <class T>
  apl_complex<T>
  apl_complex_sub(apl_complex<T> a, apl_complex<T> b);
  static void
  chop_small_elements1(const size_t N, const double d[], double sd[]);

  static void
  qrstep1(const size_t n, double d[], double sd[], double gc[], double gs[]);
  inline static double
  trailing_eigenvalue1(const size_t n, const double d[], const double sd[]);
  void PRINTM(apl_matrix_complex<double> *A);
  void PRINTV(apl_vector_complex<double> *V);
  void PRINTVD(apl_vector<double> *V);
  bool isorthogonal(apl_matrix_complex<double> *m);
  int apl_eigen_hermv_sort(apl_vector<double> *eval, apl_matrix_complex<double> *evec,
                           apl_eigen_sort_t sort_type);
  xvector<double>
  calculate_gruneisen(xmatrix<xcomplex<double> > &m0, xmatrix<xcomplex<double> > &mP, xmatrix<xcomplex<double> > &mM, double delV, double V_0);
  int apl_blas_zgemv(APL_TRANSPOSE_t TransA, const apl_complex<double> alpha,
                     const apl_matrix_complex<double> *A, const apl_vector_complex<double> *X,
                     const apl_complex<double> beta, apl_vector_complex<double> *Y);

  void
  cblas_zgemv(const enum APL_ORDER order, const enum APL_TRANSPOSE TransA,
              const int M, const int N, const void *alpha, const void *A,
              const int lda, const void *X, const int incX, const void *beta,
              void *Y, const int incY);
  xvector<int>
  trace_acoustic_mode(const xmatrix<xcomplex<double> > &mk1, const xmatrix<xcomplex<double> > &mk2);
};
#endif
}
