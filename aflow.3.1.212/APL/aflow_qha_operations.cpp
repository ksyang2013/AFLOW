// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *                Aflow PINKU NATH - Duke University 2014-2017             *
// *                                                                         *
// ***************************************************************************
// Written by Pinku Nath
// pn49@duke.edu
//this cpp file perform vector and matrix operations

#include "aflow_qha_operations.h"

#define OFFSET(N, incX) ((incX) > 0 ? 0 : ((N)-1) * (-(incX)))

namespace apl {
// ******************************************************
static const apl_multifit_fdfsolver_type<double> lmsder_type =
    {              //This is a robust and efficient version of the Levenberg-Marquardt algorithm
        "lmsder",  // name //
        sizeof(lmder_state_t<double>),
        &MVops::lmder_alloc<double>,
        &MVops::lmsder_set<double>,
        &MVops::lmsder_iterate<double>,
        &MVops::lmder_free<double>};
// ******************************************************

//static const apl_multifit_fdfsolver_type<double> lmder_type =
//    {             //Levenberg-Marquardt solver
//        "lmder",  // name //
//        sizeof(lmder_state_t<double>),
//        &MVops::lmder_alloc<double>,
//        &MVops::lmder_set<double>,
//        &MVops::lmder_iterate<double>,
//        &MVops::lmder_free<double>};

// ******************************************************
struct data {
  size_t n;
  double *y;
  double *sigma;
  double *xdata;
};
// ******************************************************
int Birch_Murnaghan_4th_order(const apl_vector<double> *x, void *data,
                              apl_vector<double> *f) {
  //Birch-Murnaghan 4th order function
  //(Birch F, Phys. Rev. 71, p809 (1947))
  MVops fit;
  size_t n = ((struct data *)data)->n;
  double *y = ((struct data *)data)->y;
  double *sigma = ((struct data *)data)->sigma;
  double *xdata = ((struct data *)data)->xdata;

  double V0 = fit.apl_vector_get<double>(x, 0);
  double E0 = fit.apl_vector_get<double>(x, 1);
  double B0 = fit.apl_vector_get<double>(x, 2);
  double B0p = fit.apl_vector_get<double>(x, 3);
  double B0pp = fit.apl_vector_get<double>(x, 4);

  size_t i;
  for (i = 0; i < n; i++) {
    double V = xdata[i];
    double power1 = pow((V0 / V), (2.0 / 3.0));
    double f1_V = (1.0 / 2.0) * (power1 - 1.0);
    double T1_V = (9.0 / 2.0) * B0 * V0 * f1_V * f1_V;
    double T2_V = 1.0 + (B0p - 4.0) * f1_V;
    double T3_V = (3.0 / 4.0) * (B0 / B0pp + B0p * (B0p - 7.0) + 143.0 / 9.0) * f1_V * f1_V;
    double A = E0 + T1_V * (T2_V + T3_V);

    fit.apl_vector_set<double>(f, i, (A - y[i]) / sigma[i]);
  }
  return 0;
}
// ******************************************************
int Birch_Murnaghan_4th_order_df(const apl_vector<double> *x, void *data,
                                 apl_matrix<double> *J) {
  //Birch-Murnaghan function derivative w.r.t parameters
  MVops fit;
  size_t n = ((struct data *)data)->n;
  double *sigma = ((struct data *)data)->sigma;
  double *xdata = ((struct data *)data)->xdata;

  double V0 = fit.apl_vector_get<double>(x, 0);
  double E0 = fit.apl_vector_get<double>(x, 1);
  (void)E0;  //to avoid compiler warning: [unused variable ‘E0’]
  double B0 = fit.apl_vector_get<double>(x, 2);
  double B0p = fit.apl_vector_get<double>(x, 3);
  double B0pp = fit.apl_vector_get<double>(x, 4);

  size_t i;

  for (i = 0; i < n; i++) {
    // Jacobian matrix J(i,j) = dfi / dxj, //
    double V = xdata[i];
    double s = sigma[i];
    double D_myf_E0 = 1.0;
    double powV1 = pow(V0 / V, 0.666667);
    double powV2 = pow(V0 / V, 1.33333);
    double powV3 = pow(V0 / V, 2);
    double D_myf_V0 = (1. / B0pp) * 1.19531 * B0 * (-1. + powV1) * (-0.176471 * B0 + (-5.62745 + (1.70588 - 0.176471 * B0p) * B0p) * B0pp + (1. * B0 + (25.6144 + B0p * (-8.88235 + 1. * B0p)) * B0pp) * powV1 + (-1.47059 * B0 + (-29.0131 + (11.7059 - 1.47059 * B0p) * B0p) * B0pp) * powV2 + (0.647059 * B0 + (10.281 + (-4.52941 + 0.647059 * B0p) * B0p) * B0pp) * powV3);
    double sq1 = (-1. + powV1) * (-1. + powV1);
    double sq2 = (-1. + powV1) * (-1. + powV1) * (-1. + powV1) * (-1. + powV1);
    double D_myf_B0 = 1.125 * V0 * (1. + 0.1875 * (15.8889 + (-7. + B0p) * B0p + B0 / B0pp) * sq1 + 0.5 * (-4. + B0p) * (-1. + powV1)) * sq1 + (0.210938 * B0 * V0 * sq2) / B0pp;
    double D_myf_B0p = 1.125 * B0 * V0 * (0.5 + 0.375 * (-3.5 + B0p) * (-1. + powV1)) * (-1.0 + powV1) * sq1;
    double D_myf_B0pp = -((0.210938 * B0 * B0 * V0 * sq2) / (B0pp * B0pp));

    fit.apl_matrix_set<double>(J, i, 0, D_myf_E0 / s);
    fit.apl_matrix_set<double>(J, i, 1, D_myf_B0 / s);
    fit.apl_matrix_set<double>(J, i, 2, D_myf_V0 / s);
    fit.apl_matrix_set<double>(J, i, 3, D_myf_B0p / s);
    fit.apl_matrix_set<double>(J, i, 4, D_myf_B0pp / s);
  }
  return 0;
}
// ******************************************************
int Birch_Murnaghan_4th_order_fdf(const apl_vector<double> *x, void *data,
                                  apl_vector<double> *f, apl_matrix<double> *J) {
  //Birch_Murnaghan_4th_order function total derivative
  Birch_Murnaghan_4th_order(x, data, f);
  Birch_Murnaghan_4th_order_df(x, data, J);
  return 0;
}
// ******************************************************
int Birch_Murnaghan_3rd_order(const apl_vector<double> *x, void *data,
                              apl_vector<double> *f) {
  //Birch-Murnaghan 3rd order function
  //(Birch F, Phys. Rev. 71, p809 (1947))
  MVops fit;
  size_t n = ((struct data *)data)->n;
  double *y = ((struct data *)data)->y;
  double *sigma = ((struct data *)data)->sigma;
  double *xdata = ((struct data *)data)->xdata;

  double V0 = fit.apl_vector_get<double>(x, 0);
  double E0 = fit.apl_vector_get<double>(x, 1);
  double B0 = fit.apl_vector_get<double>(x, 2);
  double B0p = fit.apl_vector_get<double>(x, 3);

  size_t i;
  for (i = 0; i < n; i++) {
    double V = xdata[i];
    double t1 = pow((V0 / V), (1.0 / 3.0));
    double t2 = t1 * t1;
    double t3 = t2 - 1.0;
    double A = 9.0 / 8.0 * B0 * V0 * t3 * t3 * (B0p * t3 / 2.0 - 2.0 * t2 + 3.0) + E0;
    fit.apl_vector_set<double>(f, i, (A - y[i]) / sigma[i]);
  }

  return 0;
}
// ******************************************************
int Birch_Murnaghan_3rd_order_df(const apl_vector<double> *x, void *data,
                                 apl_matrix<double> *J) {
  //Birch-Murnaghan function derivative w.r.t parameters
  MVops fit;
  size_t n = ((struct data *)data)->n;
  double *sigma = ((struct data *)data)->sigma;
  double *xdata = ((struct data *)data)->xdata;

  double V0 = fit.apl_vector_get<double>(x, 0);
  double E0 = fit.apl_vector_get<double>(x, 1);
  (void)E0;  //to avoid compiler warning: [unused variable ‘E0’]
  double B0 = fit.apl_vector_get<double>(x, 2);
  double B0p = fit.apl_vector_get<double>(x, 3);

  size_t i;

  for (i = 0; i < n; i++) {
    // Jacobian matrix J(i,j) = dfi / dxj, //
    double V = xdata[i];
    double s = sigma[i];
    double D_myf_E0 = 1.0;
    double pow1 = pow(V0 / V, 0.333333);
    double pow2 = pow(V0 / V, 0.666667);
    double pow3 = pow(V0 / V, 1.33333);
    double pow4 = pow(V0 / V, 1.66667);
    double pow5 = pow(V0 / V, 2.33333);
    double D_myf_V0 = 1. / (V * pow1) * B0 * (V0 * (-15. + 2.8125 * B0p + (10.5 - 2.25 * B0p) * pow2 + (-4.5 + 1.125 * B0p) * pow3) + V * ((3.375 - 0.5625 * B0p) * pow1 + (7.875 - 1.6875 * B0p) * pow4 + (-2.25 + 0.5625 * B0p) * pow5));
    double D_myf_B0 = 1.125 * V0 * (-1. + pow2) * (-1. + pow2) * (3. - 0.5 * B0p + (-2. + 0.5 * B0p) * pow2);
    double D_myf_B0p = 0.5625 * B0 * V0 * (-1. + pow2) * (-1. + pow2) * (-1. + pow2);

    fit.apl_matrix_set<double>(J, i, 0, D_myf_E0 / s);
    fit.apl_matrix_set<double>(J, i, 1, D_myf_B0 / s);
    fit.apl_matrix_set<double>(J, i, 2, D_myf_V0 / s);
    fit.apl_matrix_set<double>(J, i, 3, D_myf_B0p / s);
  }
  return 0;
}
// ******************************************************
int Birch_Murnaghan_3rd_order_fdf(const apl_vector<double> *x, void *data,
                                  apl_vector<double> *f, apl_matrix<double> *J) {
  //Birch_Murnaghan_3rd_order function total derivative
  Birch_Murnaghan_3rd_order(x, data, f);
  Birch_Murnaghan_3rd_order_df(x, data, J);

  return 0;
}
// ******************************************************
int Birch_Murnaghan(const apl_vector<double> *x, void *data,
                    apl_vector<double> *f) {
  //Birch-Murnaghan function
  MVops fit;
  size_t n = ((struct data *)data)->n;
  double *y = ((struct data *)data)->y;
  double *sigma = ((struct data *)data)->sigma;
  double *xdata = ((struct data *)data)->xdata;

  double Vo = fit.apl_vector_get<double>(x, 0);
  double Eo = fit.apl_vector_get<double>(x, 1);
  double Bo = fit.apl_vector_get<double>(x, 2);
  double Bp = fit.apl_vector_get<double>(x, 3);

  size_t i;

  for (i = 0; i < n; i++) {
    double V = xdata[i];
    double Yi = Eo + (Bo * V) / Bp * (pow((Vo / V), Bp) / (Bp - 1.) + 1.) - Vo * Bo / (Bp - 1.);
    fit.apl_vector_set<double>(f, i, (Yi - y[i]) / sigma[i]);
  }
  return 0;
}
// ******************************************************
int Birch_Murnaghan_df(const apl_vector<double> *x, void *data,
                       apl_matrix<double> *J) {
  //Birch-Murnaghan function derivative w.r.t parameters
  MVops fit;
  size_t n = ((struct data *)data)->n;
  double *sigma = ((struct data *)data)->sigma;
  double *xdata = ((struct data *)data)->xdata;

  double Vo = fit.apl_vector_get<double>(x, 0);
  double Eo = fit.apl_vector_get<double>(x, 1);
  (void)Eo;  //to avoid compiler warning: [unused variable ‘Eo’]
  double Bo = fit.apl_vector_get<double>(x, 2);
  double Bp = fit.apl_vector_get<double>(x, 3);

  size_t i;

  for (i = 0; i < n; i++) {
    // Jacobian matrix J(i,j) = dfi / dxj, //
    double V = xdata[i];
    double s = sigma[i];
    double D_myf_Eo = 1.0;
    double D_myf_Bo = -(Vo / (-1. + Bp)) + (V * (1. + pow((Vo / V), Bp) / (-1. + Bp))) / Bp;
    double D_myf_Vo = -(Bo / (-1. + Bp)) + (Bo * pow((Vo / V), (-1.0 + Bp))) / (-1.0 + Bp);
    double D_myf_Bp = (Bo * Vo) / ((-1. + Bp) * (-1. + Bp)) - (Bo * V * (1. + pow((Vo / V), Bp) / (-1. + Bp))) / (Bp * Bp) +
                      (Bo * V * (-(pow((Vo / V), Bp) / (-1. + Bp) * (-1. + Bp)) + (pow((Vo / V), Bp) * log(Vo / V)) / (-1.0 + Bp))) / Bp;

    fit.apl_matrix_set<double>(J, i, 0, D_myf_Eo / s);
    fit.apl_matrix_set<double>(J, i, 1, D_myf_Bo / s);
    fit.apl_matrix_set<double>(J, i, 2, D_myf_Vo / s);
    fit.apl_matrix_set<double>(J, i, 3, D_myf_Bp / s);
  }
  return 0;
}
// ******************************************************
int Birch_Murnaghan_fdf(const apl_vector<double> *x, void *data,
                        apl_vector<double> *f, apl_matrix<double> *J) {
  //Birch_Murnaghan function total derivative
  Birch_Murnaghan(x, data, f);
  Birch_Murnaghan_df(x, data, J);

  return 0;
}
// ******************************************************
template <class T>
apl_block<T> *MVops::apl_block_alloc(const size_t n) {
#define MULTIPLICITY 1
  apl_block<T> *b;

  if (n == 0) {
    cerr << "block length n must be positive integer" << std::endl;
    exit(0);
  }

  b = (apl_block<T> *)malloc(sizeof(apl_block<T>));

  if (b == 0) {
    cerr << "failed to allocate space for block struct" << std::endl;
    exit(0);
  }

  b->data = (T *)calloc(1, MULTIPLICITY * n * sizeof(T));

  if (b->data == 0) {
    free(b);  // exception in constructor, avoid memory leak //

    cerr << "failed to allocate space for block data" << std::endl;
  }

  b->size = n;
#undef MULTIPLICITY
  return b;
}
// ******************************************************
void md_lsquares::clear() {
  {
    data_read_error = false;  //no error
    nl_success_status = 0;    //success
    nl_err_msg = "";
  }

  {
    Xdata.clear();
    Ydata.clear();
    guess.clear();
  }

  {
    luncertanity_V0 = 0.0;
    luncertanity_E0 = 0.0;
    luncertanity_B0 = 0.0;
    lchisq = 0.0;
    leqmV0 = 0.0;
    leqmE0 = 0.0;
    leqmB0 = 0.0;
  }

  {
    uncertanity_V0 = 0.0;
    uncertanity_E0 = 0.0;
    uncertanity_B0 = 0.0;
    uncertanity_Bp = 0.0;
    chisq_dof = 0.0;  //chi square per degrees of freedom
    nleqmV0 = 0.0;
    nleqmE0 = 0.0;
    nleqmB0 = 0.0;
    nleqmBp = 0.0;
  }
}
// ******************************************************
void md_lsquares::cubic_polynomial_fit() {
  if (Xdata.size() != Ydata.size()) {
    data_read_error = true;
    return;
  } else if ((Xdata.size() == 0) || Ydata.size() == 0) {
    data_read_error = true;
    return;
  }

  int i, n;
  double chisq;
  apl_matrix<double> *X, *cov;
  apl_vector<double> *y, *w, *c;

  n = Xdata.size();

  X = apl_matrix_alloc<double>(n, 3);
  y = apl_vector_alloc<double>(n);
  w = apl_vector_alloc<double>(n);

  c = apl_vector_alloc<double>(3);
  cov = apl_matrix_alloc<double>(3, 3);

  for (i = 0; i != (int)Xdata.size(); i++) {
    apl_matrix_set<double>(X, i, 0, 1.0);
    apl_matrix_set<double>(X, i, 1, Xdata[i]);
    apl_matrix_set<double>(X, i, 2, Xdata[i] * Xdata[i]);
    apl_vector_set<double>(y, i, Ydata[i]);
    apl_vector_set<double>(w, i, 1.0 / (spread * spread));
  }

  {
    apl_multifit_linear_workspace<double> *work = apl_multifit_linear_alloc<double>(n, 3);
    apl_multifit_wlinear<double>(X, w, y, c, cov, &chisq, work);
    apl_multifit_linear_free<double>(work);
  }

#define C(i) (apl_vector_get<double>(c, (i)))
#define COV(i, j) (apl_matrix_get<double>(cov, (i), (j)))
#define ERR(i) sqrt(apl_matrix_get<double>(cov, i, i))
  {
    lchisq = chisq;
    luncertanity_V0 = ERR(0);
    luncertanity_E0 = ERR(1);
    luncertanity_B0 = ERR(2);
  }

  double V0 = -C(1) / (2. * C(2));
  leqmV0 = V0;
  leqmE0 = C(2) * V0 * V0 + C(1) * V0 + C(0);
  leqmB0 = 2. * C(2) * V0;

  guess.resize(4, 0.0);
  guess[0] = V0;                                 //Volulme
  guess[1] = C(2) * V0 * V0 + C(1) * V0 + C(0);  //Energy
  guess[2] = 2. * C(2) * V0;                     //Bulk Modulus
  guess[3] = 4.0;                                //from experience

  apl_matrix_free<double>(X);
  apl_vector_free<double>(y);
  apl_vector_free<double>(w);
  apl_vector_free<double>(c);
  apl_matrix_free<double>(cov);

  birch_murnaghan_fit();
#undef C
#undef COV
#undef ERR
}
// ******************************************************
void md_lsquares::birch_murnaghan_fit() {
  //Birch Murnaghan data fit
  int datapoints = (int)Xdata.size();

  int status;
  unsigned int i, iter = 0;
  const size_t n = datapoints;
  const size_t p = 4;

  apl_matrix<double> *covar = apl_matrix_alloc<double>(p, p);
  double *y, *sigma, *xdata;
  y = (double *)malloc(sizeof(double) * datapoints);
  sigma = (double *)malloc(sizeof(double) * datapoints);
  xdata = (double *)malloc(sizeof(double) * datapoints);
  if ((y == NULL) || (sigma == NULL) || (xdata == NULL)) {
    cerr << " not able to allocate the memories \n";
    exit(0);
  }
  struct data d = {n, y, sigma, xdata};
  apl_multifit_function_fdf f;

  //guess value
  double x_init[4];

  x_init[0] = guess[0];
  x_init[1] = guess[1];
  x_init[2] = guess[2];
  x_init[3] = guess[3];
  apl_vector_view<double> x = apl_vector_view_array<double>(x_init, p);

  f.f = &Birch_Murnaghan;
  f.df = &Birch_Murnaghan_df;
  f.fdf = &Birch_Murnaghan_fdf;
  f.n = n;
  f.p = p;
  f.params = &d;

  // This is the data to be fitted
  for (i = 0; i != Xdata.size(); i++) {
    xdata[i] = Xdata[i];
    y[i] = Ydata[i];
    sigma[i] = spread;
  }

#ifdef APL_multifit_fdfsolver_lmsder
  const apl_multifit_fdfsolver_type<double> *apl_multifit_fdfsolver_lmsder = &lmsder_type;
#endif
#ifdef APL_multifit_fdfsolver_lmder
  const apl_multifit_fdfsolver_type<double> *apl_multifit_fdfsolver_lmsder = &lmder_type;
#endif

  const apl_multifit_fdfsolver_type<double> *T;
  T = apl_multifit_fdfsolver_lmsder;
  apl_multifit_fdfsolver<double> *s;

  s = apl_multifit_fdfsolver_alloc<double>(T, n, p);
  apl_multifit_fdfsolver_set<double>(s, &f, &x.vector);

  do {
    iter++;
    status = apl_multifit_fdfsolver_iterate(s);

    if (status)
      break;
    status = apl_multifit_test_delta<double>(s->dx, s->x,
                                             Absolute_Error, Relative_Error);
  } while (status == -2 && iter < 500);
  apl_multifit_covar<double>(s->J, 0.0, covar);

  nl_success_status = status;

  if (status == 0)
    nl_err_msg = "status = success";
  else if (status == 29)
    nl_err_msg = "cannot reach the specified tolerance in F";
  else if (status == 30)
    nl_err_msg = "cannot reach the specified tolerance in X";
  else if (status == 31)
    nl_err_msg = "cannot reach the specified tolerance in gradient";
  else if (status == 37)
    nl_err_msg = "iteration is not making progress towards solution";
  else
    nl_err_msg = "unknown error";

  fdfsolver_name = apl_multifit_fdfsolver_name<double>(s);

#define FIT(i) apl_vector_get<double>(s->x, i)
#define ERR(i) sqrt(apl_matrix_get<double>(covar, i, i))

  {
    double chi = apl_blas_dnrm2(s->f);
    double dof = n - p;
    chisq_dof = pow(chi, 2.0) / dof;
  }

  uncertanity_V0 = ERR(0);
  uncertanity_E0 = ERR(1);
  uncertanity_B0 = ERR(2);
  uncertanity_Bp = ERR(3);
  nleqmV0 = FIT(0);
  nleqmE0 = FIT(1);
  nleqmB0 = FIT(2);
  nleqmBp = FIT(3);
#undef FIT
#undef ERR

  free(y);
  free(sigma);
  free(xdata);
  apl_multifit_fdfsolver_free<double>(s);
  apl_matrix_free<double>(covar);
}
// ******************************************************
void md_lsquares::birch_murnaghan_4th_order_fit(const xvector<double> &user_guess) {
  //Birch Murnaghan 4th_order data fit
  if (user_guess.rows != 4) {
    cerr << "birch_murnaghan_4th_order_fit can not be performed user_guess.rows!=4 \n";
    return;
  }
  if ((Xdata.size() == 0) || (Ydata.size() == 0)) {
    cerr << "birch_murnaghan_4th_order_fit can not be performed, data sizes are zero \n";
    return;
  }

  int datapoints = (int)Xdata.size();

  int status;
  unsigned int i, iter = 0;
  const size_t n = datapoints;
  const size_t p = 5;

  apl_matrix<double> *covar = apl_matrix_alloc<double>(p, p);
  double *y, *sigma, *xdata;
  y = (double *)malloc(sizeof(double) * datapoints);
  sigma = (double *)malloc(sizeof(double) * datapoints);
  xdata = (double *)malloc(sizeof(double) * datapoints);
  if ((y == NULL) || (sigma == NULL) || (xdata == NULL)) {
    cerr << " not able to allocate the memories \n";
    exit(0);
  }
  struct data d = {n, y, sigma, xdata};
  apl_multifit_function_fdf f;

  //guess value
  double x_init[5];

  //guess values from Birch-Murnaghan fits
  x_init[0] = user_guess[1];  //V0
  x_init[1] = user_guess[2];  //E0
  x_init[2] = user_guess[3];  //B0
  x_init[3] = user_guess[4];  //Bp
  x_init[4] = -1.0;
  apl_vector_view<double> x = apl_vector_view_array<double>(x_init, p);

  f.f = &Birch_Murnaghan_4th_order;
  f.df = &Birch_Murnaghan_4th_order_df;
  f.fdf = &Birch_Murnaghan_4th_order_fdf;
  f.n = n;
  f.p = p;
  f.params = &d;

  // This is the data to be fitted
  for (i = 0; i != Xdata.size(); i++) {
    xdata[i] = Xdata[i];
    y[i] = Ydata[i];
    sigma[i] = spread;
  }

#ifdef APL_multifit_fdfsolver_lmsder
  const apl_multifit_fdfsolver_type<double> *apl_multifit_fdfsolver_lmsder = &lmsder_type;
#endif
#ifdef APL_multifit_fdfsolver_lmder
  const apl_multifit_fdfsolver_type<double> *apl_multifit_fdfsolver_lmsder = &lmder_type;
#endif

  const apl_multifit_fdfsolver_type<double> *T;
  T = apl_multifit_fdfsolver_lmsder;
  apl_multifit_fdfsolver<double> *s;

  s = apl_multifit_fdfsolver_alloc<double>(T, n, p);
  apl_multifit_fdfsolver_set<double>(s, &f, &x.vector);

  do {
    iter++;
    status = apl_multifit_fdfsolver_iterate(s);

    if (status)
      break;
    status = apl_multifit_test_delta<double>(s->dx, s->x,
                                             Absolute_Error, Relative_Error);
  } while (status == -2 && iter < 500);
  apl_multifit_covar<double>(s->J, 0.0, covar);

  nl_success_status = status;

  if (status == 0)
    nl_err_msg = "status = success";
  else if (status == 29)
    nl_err_msg = "cannot reach the specified tolerance in F";
  else if (status == 30)
    nl_err_msg = "cannot reach the specified tolerance in X";
  else if (status == 31)
    nl_err_msg = "cannot reach the specified tolerance in gradient";
  else if (status == 37)
    nl_err_msg = "iteration is not making progress towards solution";
  else
    nl_err_msg = "unknown error";

  fdfsolver_name = apl_multifit_fdfsolver_name<double>(s);

#define FIT(i) apl_vector_get<double>(s->x, i)
#define ERR(i) sqrt(apl_matrix_get<double>(covar, i, i))

  {
    double chi = apl_blas_dnrm2(s->f);
    double dof = n - p;
    chisq_dof = pow(chi, 2.0) / dof;
  }

  uncertanity_V0 = ERR(0);
  uncertanity_E0 = ERR(1);
  uncertanity_B0 = ERR(2);
  uncertanity_Bp = ERR(3);
  uncertanity_Bpp = ERR(4);
  nleqmV0 = FIT(0);
  nleqmE0 = FIT(1);
  nleqmB0 = FIT(2);
  nleqmBp = FIT(3);
  nleqmBpp = FIT(4);
#undef FIT
#undef ERR
  free(y);
  free(sigma);
  free(xdata);
  apl_multifit_fdfsolver_free<double>(s);
  apl_matrix_free<double>(covar);
}
// ******************************************************
void md_lsquares::birch_murnaghan_3rd_order_fit(const xvector<double> &user_guess) {
  //Birch Murnaghan 4th_order data fit
  if (user_guess.rows != 4) {
    cerr << "birch_murnaghan_4th_order_fit can not be performed user_guess.rows!=4 \n";
    return;
  }
  if ((Xdata.size() == 0) || (Ydata.size() == 0)) {
    cerr << "birch_murnaghan_3rd_order_fit can not be performed, data sizes are zero \n";
    return;
  }

  int datapoints = (int)Xdata.size();

  int status;
  unsigned int i, iter = 0;
  const size_t n = datapoints;
  const size_t p = 4;

  apl_matrix<double> *covar = apl_matrix_alloc<double>(p, p);
  double *y, *sigma, *xdata;
  y = (double *)malloc(sizeof(double) * datapoints);
  sigma = (double *)malloc(sizeof(double) * datapoints);
  xdata = (double *)malloc(sizeof(double) * datapoints);
  if ((y == NULL) || (sigma == NULL) || (xdata == NULL)) {
    cerr << " not able to allocate the memories \n";
    exit(0);
  }
  struct data d = {n, y, sigma, xdata};
  apl_multifit_function_fdf f;

  //guess value
  double x_init[4];

  //guess values from Birch-Murnaghan fits
  x_init[0] = user_guess[1];  //V0
  x_init[1] = user_guess[2];  //E0
  x_init[2] = user_guess[3];  //B0
  x_init[3] = user_guess[4];  //Bp
  apl_vector_view<double> x = apl_vector_view_array<double>(x_init, p);

  f.f = &Birch_Murnaghan_3rd_order;
  f.df = &Birch_Murnaghan_3rd_order_df;
  f.fdf = &Birch_Murnaghan_3rd_order_fdf;
  f.n = n;
  f.p = p;
  f.params = &d;

  // This is the data to be fitted
  for (i = 0; i != Xdata.size(); i++) {
    xdata[i] = Xdata[i];
    y[i] = Ydata[i];
    sigma[i] = spread;
  }

#ifdef APL_multifit_fdfsolver_lmsder
  const apl_multifit_fdfsolver_type<double> *apl_multifit_fdfsolver_lmsder = &lmsder_type;
#endif
#ifdef APL_multifit_fdfsolver_lmder
  const apl_multifit_fdfsolver_type<double> *apl_multifit_fdfsolver_lmsder = &lmder_type;
#endif

  const apl_multifit_fdfsolver_type<double> *T;
  T = apl_multifit_fdfsolver_lmsder;
  apl_multifit_fdfsolver<double> *s;

  s = apl_multifit_fdfsolver_alloc<double>(T, n, p);
  apl_multifit_fdfsolver_set<double>(s, &f, &x.vector);

  do {
    iter++;
    status = apl_multifit_fdfsolver_iterate(s);

    if (status)
      break;
    status = apl_multifit_test_delta<double>(s->dx, s->x,
                                             Absolute_Error, Relative_Error);
  } while (status == -2 && iter < 500);
  apl_multifit_covar<double>(s->J, 0.0, covar);

  nl_success_status = status;

  if (status == 0)
    nl_err_msg = "status = success";
  else if (status == 29)
    nl_err_msg = "cannot reach the specified tolerance in F";
  else if (status == 30)
    nl_err_msg = "cannot reach the specified tolerance in X";
  else if (status == 31)
    nl_err_msg = "cannot reach the specified tolerance in gradient";
  else if (status == 37)
    nl_err_msg = "iteration is not making progress towards solution";
  else
    nl_err_msg = "unknown error";

  fdfsolver_name = apl_multifit_fdfsolver_name<double>(s);

#define FIT(i) apl_vector_get<double>(s->x, i)
#define ERR(i) sqrt(apl_matrix_get<double>(covar, i, i))

  {
    double chi = apl_blas_dnrm2(s->f);
    double dof = n - p;
    chisq_dof = pow(chi, 2.0) / dof;
  }

  uncertanity_V0 = ERR(0);
  uncertanity_E0 = ERR(1);
  uncertanity_B0 = ERR(2);
  uncertanity_Bp = ERR(3);
  nleqmV0 = FIT(0);
  nleqmE0 = FIT(1);
  nleqmB0 = FIT(2);
  nleqmBp = FIT(3);
#undef FIT
#undef ERR
  free(y);
  free(sigma);
  free(xdata);
  apl_multifit_fdfsolver_free<double>(s);
  apl_matrix_free<double>(covar);
}
// ******************************************************
template <class T>
apl_vector<T> *MVops::apl_vector_alloc(const size_t n) {
  apl_block<T> *block;
  apl_vector<T> *v;

  if (n == 0) {
    cerr << "vector length n must be positive integer" << std::endl;
    exit(0);
  }

  v = (apl_vector<T> *)malloc(sizeof(apl_vector<T>));

  if (v == 0) {
    cerr << "failed to allocate space for vector struct" << std::endl;
    exit(0);
  }

  block = apl_block_alloc<T>(n);

  if (block == 0) {
    free(v);

    cerr << "failed to allocate space for block" << std::endl;
    exit(0);
  }

  v->data = block->data;
  v->size = n;
  v->stride = 1;
  v->block = block;
  v->owner = 1;

  return v;
}
// ******************************************************
template <class T>
apl_vector<T> *
MVops::apl_vector_calloc(const size_t n) {
#define MULTIPLICITY 1
  size_t i;

  apl_vector<T> *v = apl_vector_alloc<T>(n);

  if (v == 0)
    return 0;

  // initialize vector to zero //

  for (i = 0; i < MULTIPLICITY * n; i++) {
    v->data[i] = 0;
  }
#undef MULTIPLICITY
  return v;
}
// ******************************************************
template <class T>
void MVops::apl_vector_free(apl_vector<T> *v) {
  RETURN_IF_NULL(v);

  if (v->owner) {
    apl_block_free(v->block);
  }
  free(v);
}
// ******************************************************
template <class T>
void MVops::apl_block_free(apl_block<T> *b) {
  if (b == NULL) return;
  free(b->data);
  free(b);
}
// ******************************************************
template <class T>
inline void MVops::apl_vector_set(apl_vector<T> *v, const size_t i, T x) {
#if APL_RANGE_CHECK
  if (i >= v->size) {
    cerr << "index out of range" << std::endl;
    exit(0);
  }
#endif
  v->data[i * v->stride] = x;
}
// ******************************************************
template <class T>
inline T MVops::apl_vector_get(const apl_vector<T> *v, const size_t i) {
#if APL_RANGE_CHECK
  if (i >= v->size) {
    cerr << "index out of range" << std::endl;
    exit(0);
  }
#endif
  return v->data[i * v->stride];
}
// ******************************************************
template <class T>
apl_vector_view<T> MVops::apl_vector_subvector(apl_vector<T> *v, size_t offset, size_t n) {
#define MULTIPLICITY 1
  apl_vector_view<T> view = NULL_VECTOR_VIEW;

  if (n == 0) {
    cerr << "vector length n must be positive integer" << std::endl;
    exit(0);
  }

  if (offset + (n - 1) >= v->size) {
    cerr << "view would extend past end of vector" << std::endl;
    exit(0);
  }

  {
    apl_vector<T> s = NULL_VECTOR;

    s.data = v->data + MULTIPLICITY * v->stride * offset;
    s.size = n;
    s.stride = v->stride;
    s.block = v->block;
    s.owner = 0;

    view.vector = s;
#undef MULTIPLICITY
    return view;
  }
}
// ******************************************************
template <class T>
const apl_vector_view<T> MVops::apl_vector_const_subvector(const apl_vector<T> *v, size_t offset, size_t n) {
#define MULTIPLICITY 1
  apl_vector_view<T> view = NULL_VECTOR_VIEW;

  if (n == 0) {
    cerr << "vector length n must be positive integer" << std::endl;
    exit(0);
  }

  if (offset + (n - 1) >= v->size) {
    cerr << "view would extend past end of vector" << std::endl;
    exit(0);
  }

  {
    apl_vector<T> s = NULL_VECTOR;

    s.data = v->data + MULTIPLICITY * v->stride * offset;
    s.size = n;
    s.stride = v->stride;
    s.block = v->block;
    s.owner = 0;

    view.vector = s;
#undef MULTIPLICITY
    return view;
  }
}
// ******************************************************
template <class T>
int MVops::apl_vector_scale(apl_vector<T> *a, const double x) {
  const size_t N = a->size;
  const size_t stride = a->stride;

  size_t i;

  for (i = 0; i < N; i++) {
    a->data[i * stride] *= x;
  }

  return 0;
}
// ******************************************************
template <class T>
int MVops::apl_vector_div(apl_vector<T> *a, const apl_vector<T> *b) {
  const size_t N = a->size;

  if (b->size != N) {
    cerr << "vectors must have same length" << std::endl;
    exit(0);
  } else {
    const size_t stride_a = a->stride;
    const size_t stride_b = b->stride;

    size_t i;

    for (i = 0; i < N; i++) {
      a->data[i * stride_a] /= b->data[i * stride_b];
    }

    return 0;
  }
}
// ******************************************************
template <class T>
void MVops::apl_vector_set_all(apl_vector<T> *v, T x) {
#define MULTIPLICITY 1
  T *const data = v->data;
  const size_t n = v->size;
  const size_t stride = v->stride;

  size_t i;

  for (i = 0; i < n; i++) {
    *(T *)(data + MULTIPLICITY * i * stride) = x;
  }
#undef MULTIPLICITY
}
// ******************************************************
template <class T>
int MVops::apl_vector_memcpy(apl_vector<T> *dest,
                             const apl_vector<T> *src) {
#define MULTIPLICITY 1
  const size_t src_size = src->size;
  const size_t dest_size = dest->size;

  if (src_size != dest_size) {
    cerr << "vector lengths are not equal" << std::endl;
    exit(0);
  }

  {
    const size_t src_stride = src->stride;
    const size_t dest_stride = dest->stride;
    size_t j;

    for (j = 0; j < src_size; j++) {
      size_t k;

      for (k = 0; k < MULTIPLICITY; k++) {
        dest->data[MULTIPLICITY * dest_stride * j + k] = src->data[MULTIPLICITY * src_stride * j + k];
      }
    }
  }

#undef MULTIPLICITY
  return 0;
}
// ******************************************************
template <class T>
int MVops::apl_vector_swap_elements(apl_vector<T> *v, const size_t i, const size_t j) {
#define MULTIPLICITY 1
  T *data = v->data;
  const size_t size = v->size;
  const size_t stride = v->stride;

  if (i >= size) {
    cerr << "first index is out of range" << std::endl;
    exit(0);
  }

  if (j >= size) {
    cerr << "second index is out of range" << std::endl;
    exit(0);
  }

  if (i != j) {
    const size_t s = MULTIPLICITY * stride;
    size_t k;

    for (k = 0; k < MULTIPLICITY; k++) {
      T tmp = data[j * s + k];
      data[j * s + k] = data[i * s + k];
      data[i * s + k] = tmp;
    }
  }

#undef MULTIPLICITY
  return 0;
}
// ******************************************************
template <class T>
apl_vector_view<T>
MVops::apl_vector_view_array(T *base, size_t n) {
  apl_vector_view<T> view = NULL_VECTOR_VIEW;

  if (n == 0) {
    cerr << "vector length n must be positive integer" << std::endl;
    exit(0);
  }

  {
    apl_vector<T> v = NULL_VECTOR;

    v.data = (T *)base;
    v.size = n;
    v.stride = 1;
    v.block = 0;
    v.owner = 0;

    view.vector = v;
    return view;
  }
}
// ******************************************************
template <class T>
void MVops::apl_vector_set_zero(apl_vector<T> *v) {
#define MULTIPLICITY 1
  T *const data = v->data;
  const size_t n = v->size;
  const size_t stride = v->stride;
  const T zero = 0.0;

  size_t i;

  for (i = 0; i < n; i++) {
    *(T *)(data + MULTIPLICITY * i * stride) = zero;
  }
#undef MULTIPLICITY
}
// ******************************************************
template <class T>
apl_matrix<T> *
MVops::apl_matrix_alloc(const size_t n1, const size_t n2) {
  apl_block<T> *block;
  apl_matrix<T> *m;

  if (n1 == 0) {
    cerr << "matrix dimension n1 must be positive integer" << std::endl;
    exit(0);
  } else if (n2 == 0) {
    cerr << "matrix dimension n2 must be positive integer" << std::endl;
    exit(0);
  }

  m = (apl_matrix<T> *)malloc(sizeof(apl_matrix<T>));

  if (m == 0) {
    cerr << "failed to allocate space for matrix struct" << std::endl;
    exit(0);
  }

  // FIXME: n1*n2 could overflow for large dimensions //

  block = apl_block_alloc<T>(n1 * n2);

  if (block == 0) {
    cerr << "failed to allocate space for block" << std::endl;
    exit(0);
  }

  m->data = block->data;
  m->size1 = n1;
  m->size2 = n2;
  m->tda = n2;
  m->block = block;
  m->owner = 1;

  return m;
}
// ******************************************************
template <class T>
void MVops::apl_matrix_free(apl_matrix<T> *m) {
  RETURN_IF_NULL(m);

  if (m->owner) {
    apl_block_free(m->block);
  }

  free(m);
}
// ******************************************************
template <class T>
inline T
MVops::apl_matrix_get(const apl_matrix<T> *m, const size_t i, const size_t j) {
  if (i >= m->size1) {
    cerr << "first index out of range" << std::endl;
    exit(0);
  } else if (j >= m->size2) {
    cerr << "second index out of range" << std::endl;
    exit(0);
  }
  return m->data[i * m->tda + j];
}
// ******************************************************
template <class T>
inline void
MVops::apl_matrix_set(apl_matrix<T> *m, const size_t i, const size_t j, const T x) {
  if (i >= m->size1) {
    cerr << "first index out of range" << std::endl;
    exit(0);
  } else if (j >= m->size2) {
    cerr << "second index out of range" << std::endl;
    exit(0);
  }
  m->data[i * m->tda + j] = x;
}
// ******************************************************
template <class T>
int MVops::apl_matrix_memcpy(apl_matrix<T> *dest,
                             const apl_matrix<T> *src) {
#define MULTIPLICITY 1
  const size_t src_size1 = src->size1;
  const size_t src_size2 = src->size2;
  const size_t dest_size1 = dest->size1;
  const size_t dest_size2 = dest->size2;

  if (src_size1 != dest_size1 || src_size2 != dest_size2) {
    cerr << "matrix sizes are different" << std::endl;
    exit(0);
  }

  {
    const size_t src_tda = src->tda;
    const size_t dest_tda = dest->tda;
    size_t i, j;

    for (i = 0; i < src_size1; i++) {
      for (j = 0; j < MULTIPLICITY * src_size2; j++) {
        dest->data[MULTIPLICITY * dest_tda * i + j] = src->data[MULTIPLICITY * src_tda * i + j];
      }
    }
  }

#undef MULTIPLICITY
  return 0;
}
// ******************************************************
template <class T>
const apl_vector_view<T>
MVops::apl_matrix_const_row(const apl_matrix<T> *m, const size_t i) {
#define MULTIPLICITY 1
  apl_vector_view<T> view = NULL_VECTOR_VIEW;

  if (i >= m->size1) {
    cerr << "row index is out of range" << std::endl;
    exit(0);
  }

  {
    apl_vector<T> v = NULL_VECTOR;

    v.data = m->data + i * MULTIPLICITY * m->tda;
    v.size = m->size2;
    v.stride = 1;
    v.block = m->block;
    v.owner = 0;

    view.vector = v;
#undef MULTIPLICITY
    return view;
  }
}
// ******************************************************
template <class T>
apl_vector_view<T>
MVops::apl_matrix_row(apl_matrix<T> *m, const size_t i) {
#define MULTIPLICITY 1
  apl_vector_view<T> view = NULL_VECTOR_VIEW;

  if (i >= m->size1) {
    cerr << "row index is out of range" << std::endl;
    exit(0);
  }

  {
    apl_vector<T> v = NULL_VECTOR;

    v.data = m->data + i * MULTIPLICITY * m->tda;
    v.size = m->size2;
    v.stride = 1;
    v.block = m->block;
    v.owner = 0;

    view.vector = v;
#undef MULTIPLICITY
    return view;
  }
}
// ******************************************************
template <class T>
apl_vector_view<T>
MVops::apl_matrix_column(apl_matrix<T> *m, const size_t j) {
#define MULTIPLICITY 1
  apl_vector_view<T> view = NULL_VECTOR_VIEW;

  if (j >= m->size2) {
    cerr << "column index is out of range" << std::endl;
    exit(0);
  }

  {
    apl_vector<T> v = NULL_VECTOR;

    v.data = m->data + j * MULTIPLICITY;
    v.size = m->size1;
    v.stride = m->tda;
    v.block = m->block;
    v.owner = 0;

    view.vector = v;
#undef MULTIPLICITY
    return view;
  }
}
// ******************************************************
template <class T>
const apl_vector_view<T>
MVops::apl_matrix_const_column(const apl_matrix<T> *m, const size_t j) {
#define MULTIPLICITY 1
  apl_vector_view<T> view = NULL_VECTOR_VIEW;

  if (j >= m->size2) {
    cerr << "column index is out of range" << std::endl;
    exit(0);
  }

  {
    apl_vector<T> v = NULL_VECTOR;

    v.data = m->data + j * MULTIPLICITY;
    v.size = m->size1;
    v.stride = m->tda;
    v.block = m->block;
    v.owner = 0;

    view.vector = v;
#undef MULTIPLICITY
    return view;
  }
}
// ******************************************************
template <class T>
apl_matrix_view<T>
MVops::apl_matrix_submatrix(apl_matrix<T> *m,
                            const size_t i, const size_t j,
                            const size_t n1, const size_t n2) {
#define MULTIPLICITY 1
  apl_matrix_view<T> view = NULL_MATRIX_VIEW;

  if (i >= m->size1) {
    cerr << "row index is out of range" << std::endl;
    exit(0);
  } else if (j >= m->size2) {
    cerr << "column index is out of range" << std::endl;
    exit(0);
  } else if (n1 == 0) {
    cerr << "first dimension must be non-zero" << std::endl;
    exit(0);
  } else if (n2 == 0) {
    cerr << "second dimension must be non-zero" << std::endl;
    exit(0);
  } else if (i + n1 > m->size1) {
    cerr << "first dimension overflows matrix" << std::endl;
    exit(0);
  } else if (j + n2 > m->size2) {
    cerr << "second dimension overflows matrix" << std::endl;
    exit(0);
  }

  {
    apl_matrix<T> s = NULL_MATRIX;

    s.data = m->data + MULTIPLICITY * (i * m->tda + j);
    s.size1 = n1;
    s.size2 = n2;
    s.tda = m->tda;
    s.block = m->block;
    s.owner = 0;

    view.matrix = s;
#undef MULTIPLICITY
    return view;
  }
}
// ******************************************************
template <class T>
void MVops::apl_matrix_set_identity(apl_matrix<T> *m) {
#define MULTIPLICITY 1
  size_t i, j;
  T *const data = m->data;
  const size_t p = m->size1;
  const size_t q = m->size2;
  const size_t tda = m->tda;

  for (i = 0; i < p; i++) {
    for (j = 0; j < q; j++) {
      *(T *)(data + MULTIPLICITY * (i * tda + j)) = ((i == j) ? 1. : 0.);
    }
  }
#undef MULTIPLICITY
}
// ******************************************************
template <class T>
int MVops::apl_matrix_swap_columns(apl_matrix<T> *m,
                                   const size_t i, const size_t j) {
#define MULTIPLICITY 1
  const size_t size1 = m->size1;
  const size_t size2 = m->size2;

  if (i >= size2) {
    cerr << "first column index is out of range" << std::endl;
    exit(0);
  }

  if (j >= size2) {
    cerr << "second column index is out of range" << std::endl;
    exit(0);
  }

  if (i != j) {
    T *col1 = m->data + MULTIPLICITY * i;
    T *col2 = m->data + MULTIPLICITY * j;

    size_t p;

    for (p = 0; p < size1; p++) {
      size_t k;
      size_t n = p * MULTIPLICITY * m->tda;

      for (k = 0; k < MULTIPLICITY; k++) {
        T tmp = col1[n + k];
        col1[n + k] = col2[n + k];
        col2[n + k] = tmp;
      }
    }
  }

#undef MULTIPLICITY
  return 0;
}
// ******************************************************
template <class T>
apl_matrix<T> *
MVops::apl_matrix_calloc(const size_t n1, const size_t n2) {
#define MULTIPLICITY 1
  size_t i;

  apl_matrix<T> *m = apl_matrix_alloc<T>(n1, n2);

  if (m == 0)
    return 0;

  // initialize matrix to zero //

  for (i = 0; i < MULTIPLICITY * n1 * n2; i++) {
    m->data[i] = 0;
  }

#undef MULTIPLICITY
  return m;
}
// ******************************************************
int MVops::apl_permute_inverse(const size_t *p, double *data, const size_t stride, const size_t n) {
#define MULTIPLICITY 1
  size_t i, k, pk;

  for (i = 0; i < n; i++) {
    k = p[i];

    while (k > i)
      k = p[k];

    if (k < i)
      continue;

    // Now have k == i, i.e the least in its cycle //

    pk = p[k];

    if (pk == i)
      continue;

    // shuffle the elements of the cycle in the inverse direction //

    {
      unsigned int a;

      double t[MULTIPLICITY];

      for (a = 0; a < MULTIPLICITY; a++)
        t[a] = data[k * stride * MULTIPLICITY + a];

      while (pk != i) {
        for (a = 0; a < MULTIPLICITY; a++) {
          double r1 = data[pk * stride * MULTIPLICITY + a];
          data[pk * stride * MULTIPLICITY + a] = t[a];
          t[a] = r1;
        }

        k = pk;
        pk = p[k];
      };

      for (a = 0; a < MULTIPLICITY; a++)
        data[pk * stride * MULTIPLICITY + a] = t[a];
    }
  }

#undef MULTIPLICITY
  return 0;
}
// ******************************************************

template <class T>
apl_block_complex<T> *MVops::apl_block_complex_alloc(const size_t n) {
#define MULTIPLICITY 2
  apl_block_complex<T> *b;

  if (n == 0) {
    cerr << "block length n must be positive integer" << std::endl;
    exit(0);
  }

  b = (apl_block_complex<T> *)malloc(sizeof(apl_block_complex<T>));

  if (b == 0) {
    cerr << "failed to allocate space for block struct" << std::endl;
    exit(0);
  }

  b->data = (T *)calloc(1, MULTIPLICITY * n * sizeof(T));

  if (b->data == 0) {
    free(b);  // exception in constructor, avoid memory leak //

    cerr << "failed to allocate space for block data" << std::endl;
  }

  b->size = n;

#undef MULTIPLICITY
  return b;
}
// ******************************************************
template <class T>
apl_vector_complex<T> *MVops::apl_vector_complex_alloc(const size_t n) {
  apl_block_complex<T> *block;
  apl_vector_complex<T> *v;

  if (n == 0) {
    cerr << "vector length n must be positive integer" << std::endl;
    exit(0);
  }

  v = (apl_vector_complex<T> *)malloc(sizeof(apl_vector_complex<T>));

  if (v == 0) {
    cerr << "failed to allocate space for vector struct" << std::endl;
    exit(0);
  }

  block = apl_block_complex_alloc<T>(n);

  if (block == 0) {
    free(v);

    cerr << "failed to allocate space for block" << std::endl;
    exit(0);
  }

  v->data = block->data;
  v->size = n;
  v->stride = 1;
  v->block = block;
  v->owner = 1;

  return v;
}
// ******************************************************
template <class T>
apl_vector_complex<T> *
MVops::apl_vector_complex_calloc(const size_t n) {
#define MULTIPLICITY 2
  size_t i;

  apl_vector_complex<T> *v = apl_vector_complex_alloc<T>(n);

  if (v == 0)
    return 0;

  // initialize vector to zero //

  for (i = 0; i < MULTIPLICITY * n; i++) {
    v->data[i] = 0;
  }

#undef MULTIPLICITY
  return v;
}
// ******************************************************
template <class T>
void MVops::apl_vector_complex_free(apl_vector_complex<T> *v) {
  RETURN_IF_NULL(v);

  if (v->owner) {
    apl_block_complex_free(v->block);
  }
  free(v);
}
// ******************************************************
template <class T>
void MVops::apl_block_complex_free(apl_block_complex<T> *b) {
  if (b == NULL) return;
  free(b->data);
  free(b);
}
// ******************************************************
template <class T>
inline void MVops::apl_vector_complex_set(apl_vector_complex<T> *v, const size_t i, apl_complex<T> z) {
#if APL_RANGE_CHECK
  if (i >= v->size) {
    cerr << "index out of range" << std::endl;
    exit(0);
  }
#endif
  *APL_MV_COMPLEX_AT(v, i) = z;
}
// ******************************************************
template <class T>
inline apl_complex<T> MVops::apl_vector_complex_get(const apl_vector_complex<T> *v, const size_t i) {
#if APL_RANGE_CHECK
  if (i >= v->size) {
    cerr << "index out of range" << std::endl;
    exit(0);
  }
#endif
  return *APL_MV_COMPLEX_AT(v, i);
}
// ******************************************************
template <class T>
apl_vector_complex_view<T> MVops::apl_vector_complex_subvector(apl_vector_complex<T> *v, size_t offset, size_t n) {
#define MULTIPLICITY 2
  apl_vector_complex_view<T> view = NULL_VECTOR_VIEW;

  if (n == 0) {
    cerr << "vector length n must be positive integer" << std::endl;
    exit(0);
  }

  if (offset + (n - 1) >= v->size) {
    cerr << "view would extend past end of vector" << std::endl;
    exit(0);
  }

  {
    apl_vector_complex<T> s = NULL_VECTOR;

    s.data = v->data + MULTIPLICITY * v->stride * offset;
    s.size = n;
    s.stride = v->stride;
    s.block = v->block;
    s.owner = 0;

    view.vector = s;
#undef MULTIPLICITY
    return view;
  }
}
// ******************************************************
template <class T>
const apl_vector_complex_view<T> MVops::apl_vector_complex_const_subvector(const apl_vector_complex<T> *v, size_t offset, size_t n) {
#define MULTIPLICITY 2
  apl_vector_complex_view<T> view = NULL_VECTOR_VIEW;

  if (n == 0) {
    cerr << "vector length n must be positive integer" << std::endl;
    exit(0);
  }

  if (offset + (n - 1) >= v->size) {
    cerr << "view would extend past end of vector" << std::endl;
    exit(0);
  }

  {
    apl_vector_complex<T> s = NULL_VECTOR;

    s.data = v->data + MULTIPLICITY * v->stride * offset;
    s.size = n;
    s.stride = v->stride;
    s.block = v->block;
    s.owner = 0;

    view.vector = s;
#undef MULTIPLICITY
    return view;
  }
}
// ******************************************************
template <class T>
int MVops::apl_vector_complex_scale(apl_vector_complex<T> *a, const apl_complex<T> x) {
  const size_t N = a->size;
  const size_t stride = a->stride;

  size_t i;

  T xr = APL_MV_REAL(x);
  T xi = APL_MV_IMAG(x);

  for (i = 0; i < N; i++) {
    T ar = a->data[2 * i * stride];
    T ai = a->data[2 * i * stride + 1];

    a->data[2 * i * stride] = ar * xr - ai * xi;
    a->data[2 * i * stride + 1] = ar * xi + ai * xr;
  }

  return 0;
}
// ******************************************************
template <class T>
int MVops::apl_vector_complex_div(apl_vector_complex<T> *a, const apl_vector_complex<T> *b) {
  const size_t N = a->size;

  if (b->size != N) {
    cerr << "\n vectors must have same length \n";
    exit(0);
  } else {
    const size_t stride_a = a->stride;
    const size_t stride_b = b->stride;

    size_t i;

    for (i = 0; i < N; i++) {
      T ar = a->data[2 * i * stride_a];
      T ai = a->data[2 * i * stride_a + 1];

      T br = b->data[2 * i * stride_b];
      T bi = b->data[2 * i * stride_b + 1];

      T s = 1.0 / hypot(br, bi);

      T sbr = s * br;
      T sbi = s * bi;

      a->data[2 * i * stride_a] = (ar * sbr + ai * sbi) * s;
      a->data[2 * i * stride_a + 1] = (ai * sbr - ar * sbi) * s;
    }

    return 0;
  }
}
// ******************************************************
template <class T>
void MVops::apl_vector_complex_set_all(apl_vector_complex<T> *v, T x) {
#define MULTIPLICITY 2
  T *const data = v->data;
  const size_t n = v->size;
  const size_t stride = v->stride;

  size_t i;

  for (i = 0; i < n; i++) {
    *(T *)(data + MULTIPLICITY * i * stride) = x;
  }
#undef MULTIPLICITY
}
// ******************************************************
template <class T>
int MVops::apl_vector_complex_memcpy(apl_vector_complex<T> *dest,
                                     const apl_vector_complex<T> *src) {
#define MULTIPLICITY 2
  const size_t src_size = src->size;
  const size_t dest_size = dest->size;

  if (src_size != dest_size) {
    cerr << "vector lengths are not equal" << std::endl;
    exit(0);
  }

  {
    const size_t src_stride = src->stride;
    const size_t dest_stride = dest->stride;
    size_t j;

    for (j = 0; j < src_size; j++) {
      size_t k;

      for (k = 0; k < MULTIPLICITY; k++) {
        dest->data[MULTIPLICITY * dest_stride * j + k] = src->data[MULTIPLICITY * src_stride * j + k];
      }
    }
  }

#undef MULTIPLICITY
  return 0;
}
// ******************************************************
template <class T>
int MVops::apl_vector_complex_swap_elements(apl_vector_complex<T> *v, const size_t i, const size_t j) {
#define MULTIPLICITY 2
  T *data = v->data;
  const size_t size = v->size;
  const size_t stride = v->stride;

  if (i >= size) {
    cerr << "first index is out of range" << std::endl;
    exit(0);
  }

  if (j >= size) {
    cerr << "second index is out of range" << std::endl;
    exit(0);
  }

  if (i != j) {
    const size_t s = MULTIPLICITY * stride;
    size_t k;

    for (k = 0; k < MULTIPLICITY; k++) {
      T tmp = data[j * s + k];
      data[j * s + k] = data[i * s + k];
      data[i * s + k] = tmp;
    }
  }

#undef MULTIPLICITY
  return 0;
}
// ******************************************************
template <class T>
apl_vector_complex_view<T>
MVops::apl_vector_complex_view_array(T *base, size_t n) {
  apl_vector_complex_view<T> view = NULL_VECTOR_VIEW;

  if (n == 0) {
    cerr << "vector length n must be positive integer" << std::endl;
    exit(0);
  }

  {
    apl_vector_complex<T> v = NULL_VECTOR;

    v.data = (T *)base;
    v.size = n;
    v.stride = 1;
    v.block = 0;
    v.owner = 0;

    view.vector = v;
    return view;
  }
}
// ******************************************************
template <class T>
void MVops::apl_vector_complex_set_zero(apl_vector_complex<T> *v) {
#define MULTIPLICITY 2
  T *const data = v->data;
  const size_t n = v->size;
  const size_t stride = v->stride;
  const T zero = 0.0;

  size_t i;

  for (i = 0; i < n; i++) {
    *(T *)(data + MULTIPLICITY * i * stride) = zero;
  }
#undef MULTIPLICITY
}
// ******************************************************
template <class T>
apl_matrix_complex<T> *
MVops::apl_matrix_complex_alloc(const size_t n1, const size_t n2) {
  apl_block_complex<T> *block;
  apl_matrix_complex<T> *m;

  if (n1 == 0) {
    cerr << "matrix dimension n1 must be positive integer" << std::endl;
    exit(0);
  } else if (n2 == 0) {
    cerr << "matrix dimension n2 must be positive integer" << std::endl;
    exit(0);
  }

  m = (apl_matrix_complex<T> *)malloc(sizeof(apl_matrix_complex<T>));

  if (m == 0) {
    cerr << "failed to allocate space for matrix struct" << std::endl;
    exit(0);
  }

  // FIXME: n1*n2 could overflow for large dimensions //

  block = apl_block_complex_alloc<T>(n1 * n2);

  if (block == 0) {
    cerr << "failed to allocate space for block" << std::endl;
    exit(0);
  }

  m->data = block->data;
  m->size1 = n1;
  m->size2 = n2;
  m->tda = n2;
  m->block = block;
  m->owner = 1;

  return m;
}
// ******************************************************
template <class T>
void MVops::apl_matrix_complex_free(apl_matrix_complex<T> *m) {
  RETURN_IF_NULL(m);

  if (m->owner) {
    apl_block_complex_free(m->block);
  }

  free(m);
}
// ******************************************************
template <class T>
apl_complex<T>
MVops::apl_matrix_complex_get(const apl_matrix_complex<T> *m, const size_t i, const size_t j) {
  if (i >= m->size1) {
    cerr << "first index out of range" << std::endl;
    exit(0);
  } else if (j >= m->size2) {
    cerr << "second index out of range" << std::endl;
    exit(0);
  }
  return *(apl_complex<T> *)(m->data + 2 * (i * m->tda + j));
}
// ******************************************************
template <class T>
inline void
MVops::apl_matrix_complex_set(apl_matrix_complex<T> *m, const size_t i, const size_t j, const apl_complex<T> x) {
  if (i >= m->size1) {
    cerr << "first index out of range" << std::endl;
    exit(0);
  } else if (j >= m->size2) {
    cerr << "second index out of range" << std::endl;
    exit(0);
  }
  *(apl_complex<T> *)(m->data + 2 * (i * m->tda + j)) = x;
}
// ******************************************************
template <class T>
int MVops::apl_matrix_complex_memcpy(apl_matrix_complex<T> *dest,
                                     const apl_matrix_complex<T> *src) {
#define MULTIPLICITY 2
  const size_t src_size1 = src->size1;
  const size_t src_size2 = src->size2;
  const size_t dest_size1 = dest->size1;
  const size_t dest_size2 = dest->size2;

  if (src_size1 != dest_size1 || src_size2 != dest_size2) {
    cerr << "matrix sizes are different" << std::endl;
    exit(0);
  }

  {
    const size_t src_tda = src->tda;
    const size_t dest_tda = dest->tda;
    size_t i, j;

    for (i = 0; i < src_size1; i++) {
      for (j = 0; j < MULTIPLICITY * src_size2; j++) {
        dest->data[MULTIPLICITY * dest_tda * i + j] = src->data[MULTIPLICITY * src_tda * i + j];
      }
    }
  }
#undef MULTIPLICITY

  return 0;
}
// ******************************************************
template <class T>
const apl_vector_complex_view<T>
MVops::apl_matrix_complex_const_row(const apl_matrix_complex<T> *m, const size_t i) {
#define MULTIPLICITY 2
  apl_vector_complex_view<T> view = NULL_VECTOR_VIEW;

  if (i >= m->size1) {
    cerr << "row index is out of range" << std::endl;
    exit(0);
  }

  {
    apl_vector_complex<T> v = NULL_VECTOR;

    v.data = m->data + i * MULTIPLICITY * m->tda;
    v.size = m->size2;
    v.stride = 1;
    v.block = m->block;
    v.owner = 0;

    view.vector = v;
#undef MULTIPLICITY
    return view;
  }
}
// ******************************************************
template <class T>
apl_vector_complex_view<T>
MVops::apl_matrix_complex_row(apl_matrix_complex<T> *m, const size_t i) {
#define MULTIPLICITY 2
  apl_vector_complex_view<T> view = NULL_VECTOR_VIEW;

  if (i >= m->size1) {
    cerr << "row index is out of range" << std::endl;
    exit(0);
  }

  {
    apl_vector_complex<T> v = NULL_VECTOR;

    v.data = m->data + i * MULTIPLICITY * m->tda;
    v.size = m->size2;
    v.stride = 1;
    v.block = m->block;
    v.owner = 0;

    view.vector = v;
#undef MULTIPLICITY
    return view;
  }
}
// ******************************************************
template <class T>
apl_vector_complex_view<T>
MVops::apl_matrix_complex_column(apl_matrix_complex<T> *m, const size_t j) {
#define MULTIPLICITY 2
  apl_vector_complex_view<T> view = NULL_VECTOR_VIEW;

  if (j >= m->size2) {
    cerr << "column index is out of range" << std::endl;
    exit(0);
  }

  {
    apl_vector_complex<T> v = NULL_VECTOR;

    v.data = m->data + j * MULTIPLICITY;
    v.size = m->size1;
    v.stride = m->tda;
    v.block = m->block;
    v.owner = 0;

    view.vector = v;
#undef MULTIPLICITY
    return view;
  }
}
// ******************************************************
template <class T>
const apl_vector_complex_view<T>
MVops::apl_matrix_complex_const_column(const apl_matrix_complex<T> *m, const size_t j) {
#define MULTIPLICITY 2
  apl_vector_complex_view<T> view = NULL_VECTOR_VIEW;

  if (j >= m->size2) {
    cerr << "column index is out of range" << std::endl;
    exit(0);
  }

  {
    apl_vector_complex<T> v = NULL_VECTOR;

    v.data = m->data + j * MULTIPLICITY;
    v.size = m->size1;
    v.stride = m->tda;
    v.block = m->block;
    v.owner = 0;

    view.vector = v;
#undef MULTIPLICITY
    return view;
  }
}
// ******************************************************
template <class T>
apl_matrix_complex_view<T>
MVops::apl_matrix_complex_submatrix(apl_matrix_complex<T> *m,
                                    const size_t i, const size_t j,
                                    const size_t n1, const size_t n2) {
#define MULTIPLICITY 2
  apl_matrix_complex_view<T> view = NULL_MATRIX_VIEW;

  if (i >= m->size1) {
    cerr << "row index is out of range" << std::endl;
    exit(0);
  } else if (j >= m->size2) {
    cerr << "column index is out of range" << std::endl;
    exit(0);
  } else if (n1 == 0) {
    cerr << "first dimension must be non-zero" << std::endl;
    exit(0);
  } else if (n2 == 0) {
    cerr << "second dimension must be non-zero" << std::endl;
    exit(0);
  } else if (i + n1 > m->size1) {
    cerr << "first dimension overflows matrix" << std::endl;
    exit(0);
  } else if (j + n2 > m->size2) {
    cerr << "second dimension overflows matrix" << std::endl;
    exit(0);
  }

  {
    apl_matrix_complex<T> s = NULL_MATRIX;

    s.data = m->data + MULTIPLICITY * (i * m->tda + j);
    s.size1 = n1;
    s.size2 = n2;
    s.tda = m->tda;
    s.block = m->block;
    s.owner = 0;

    view.matrix = s;
#undef MULTIPLICITY
    return view;
  }
}
// ******************************************************
template <class T>
void MVops::apl_matrix_complex_set_identity(apl_matrix_complex<T> *m) {
#define MULTIPLICITY 2
  size_t i, j;
  T *const data = m->data;
  const size_t p = m->size1;
  const size_t q = m->size2;
  const size_t tda = m->tda;

  for (i = 0; i < p; i++) {
    for (j = 0; j < q; j++) {
      *(T *)(data + MULTIPLICITY * (i * tda + j)) = ((i == j) ? 1. : 0.);
    }
  }
#undef MULTIPLICITY
}
// ******************************************************
template <class T>
int MVops::apl_matrix_complex_swap_columns(apl_matrix_complex<T> *m,
                                           const size_t i, const size_t j) {
#define MULTIPLICITY 2
  const size_t size1 = m->size1;
  const size_t size2 = m->size2;

  if (i >= size2) {
    cerr << "first column index is out of range" << std::endl;
    exit(0);
  }

  if (j >= size2) {
    cerr << "second column index is out of range" << std::endl;
    exit(0);
  }

  if (i != j) {
    T *col1 = m->data + MULTIPLICITY * i;
    T *col2 = m->data + MULTIPLICITY * j;

    size_t p;

    for (p = 0; p < size1; p++) {
      size_t k;
      size_t n = p * MULTIPLICITY * m->tda;

      for (k = 0; k < MULTIPLICITY; k++) {
        T tmp = col1[n + k];
        col1[n + k] = col2[n + k];
        col2[n + k] = tmp;
      }
    }
  }
#undef MULTIPLICITY

  return 0;
}
// ******************************************************
template <class T>
apl_matrix_complex<T> *
MVops::apl_matrix_complex_calloc(const size_t n1, const size_t n2) {
#define MULTIPLICITY 2
  size_t i;

  apl_matrix_complex<T> *m = apl_matrix_complex_alloc<T>(n1, n2);

  if (m == 0)
    return 0;

  // initialize matrix to zero //

  for (i = 0; i < MULTIPLICITY * n1 * n2; i++) {
    m->data[i] = 0;
  }

#undef MULTIPLICITY
  return m;
}
// ******************************************************
apl_complex<double> MVops::apl_complex_rect(double x, double y) {  // return z = x + i y //
  apl_complex<double> z;
  APL_MV_SET_COMPLEX(&z, x, y);
  return z;
}
// ******************************************************
template <class T>
double MVops::apl_complex_abs(apl_complex<T> z) {  // return |z| //
  return hypot(APL_MV_REAL(z), APL_MV_IMAG(z));
}
// ******************************************************
template <class T>
apl_complex<T>
MVops::apl_complex_sub_real(apl_complex<T> a, double x) {  // z=a-x //
  apl_complex<T> z;
  APL_MV_SET_COMPLEX(&z, APL_MV_REAL(a) - x, APL_MV_IMAG(a));
  return z;
}
// ******************************************************
template <class T>
apl_complex<T>
MVops::apl_complex_inverse(apl_complex<T> a) {  // z=1/a //
  double s = 1.0 / apl_complex_abs<T>(a);

  apl_complex<T> z;
  APL_MV_SET_COMPLEX(&z, (APL_MV_REAL(a) * s) * s, -(APL_MV_IMAG(a) * s) * s);
  return z;
}
// ******************************************************
template <class T>
apl_complex<T>
MVops::apl_complex_mul(apl_complex<T> a, apl_complex<T> b) {  // z=a*b //
  double ar = APL_MV_REAL(a), ai = APL_MV_IMAG(a);
  double br = APL_MV_REAL(b), bi = APL_MV_IMAG(b);

  apl_complex<T> z;
  APL_MV_SET_COMPLEX(&z, ar * br - ai * bi, ar * bi + ai * br);
  return z;
}
// ******************************************************
template <class T>
apl_complex<T>
MVops::apl_complex_mul_real(apl_complex<T> a, double x) {  // z=a*x //
  apl_complex<T> z;
  APL_MV_SET_COMPLEX(&z, x * APL_MV_REAL(a), x * APL_MV_IMAG(a));
  return z;
}
// ******************************************************
template <class T>
apl_complex<T>
MVops::apl_complex_conjugate(apl_complex<T> a) {  // z=conj(a) //
  apl_complex<T> z;
  APL_MV_SET_COMPLEX(&z, APL_MV_REAL(a), -APL_MV_IMAG(a));
  return z;
}
// ******************************************************
template <class T>
apl_complex<T>
MVops::apl_complex_add(apl_complex<T> a, apl_complex<T> b) {  // z=a+b //
  double ar = APL_MV_REAL(a), ai = APL_MV_IMAG(a);
  double br = APL_MV_REAL(b), bi = APL_MV_IMAG(b);

  apl_complex<T> z;
  APL_MV_SET_COMPLEX(&z, ar + br, ai + bi);
  return z;
}
// ******************************************************
template <class T>
apl_complex<T>
MVops::apl_complex_sub(apl_complex<T> a, apl_complex<T> b) {  // z=a-b //
  double ar = APL_MV_REAL(a), ai = APL_MV_IMAG(a);
  double br = APL_MV_REAL(b), bi = APL_MV_IMAG(b);

  apl_complex<T> z;
  APL_MV_SET_COMPLEX(&z, ar - br, ai - bi);
  return z;
}
// ******************************************************

template <class T>
apl_multifit_linear_workspace<T> *
MVops::apl_multifit_linear_alloc(size_t n, size_t p) {
  apl_multifit_linear_workspace<T> *w;

  w = (apl_multifit_linear_workspace<T> *)
      malloc(sizeof(apl_multifit_linear_workspace<T>));

  if (w == 0) {
    cerr << "failed to allocate space for multifit_linear struct" << std::endl;
    exit(0);
  }

  w->n = n;  // number of observations //
  w->p = p;  // number of parameters //

  w->A = apl_matrix_alloc<T>(n, p);

  if (w->A == 0) {
    free(w);
    cerr << "failed to allocate space for A" << std::endl;
    exit(0);
  }

  w->Q = apl_matrix_alloc<T>(p, p);

  if (w->Q == 0) {
    apl_matrix_free<T>(w->A);
    free(w);
    cerr << "failed to allocate space for Q" << std::endl;
    exit(0);
  }

  w->QSI = apl_matrix_alloc<T>(p, p);

  if (w->QSI == 0) {
    apl_matrix_free<T>(w->Q);
    apl_matrix_free<T>(w->A);
    free(w);
    cerr << "failed to allocate space for QSI" << std::endl;
    exit(0);
  }

  w->S = apl_vector_alloc<T>(p);

  if (w->S == 0) {
    apl_matrix_free<T>(w->QSI);
    apl_matrix_free<T>(w->Q);
    apl_matrix_free<T>(w->A);
    free(w);
    cerr << "failed to allocate space for S" << std::endl;
    exit(0);
  }

  w->t = apl_vector_alloc<T>(n);

  if (w->t == 0) {
    apl_vector_free<T>(w->S);
    apl_matrix_free<T>(w->QSI);
    apl_matrix_free<T>(w->Q);
    apl_matrix_free<T>(w->A);
    free(w);
    cerr << "failed to allocate space for t" << std::endl;
    exit(0);
  }

  w->xt = apl_vector_calloc<T>(p);

  if (w->xt == 0) {
    apl_vector_free<T>(w->t);
    apl_vector_free<T>(w->S);
    apl_matrix_free<T>(w->QSI);
    apl_matrix_free<T>(w->Q);
    apl_matrix_free<T>(w->A);
    free(w);
    cerr << "failed to allocate space for xt" << std::endl;
    exit(0);
  }

  w->D = apl_vector_calloc<T>(p);

  if (w->D == 0) {
    apl_vector_free<T>(w->D);
    apl_vector_free<T>(w->t);
    apl_vector_free<T>(w->S);
    apl_matrix_free<T>(w->QSI);
    apl_matrix_free<T>(w->Q);
    apl_matrix_free<T>(w->A);
    free(w);
    cerr << "failed to allocate space for xt" << std::endl;
    exit(0);
  }
  return w;
}
// ******************************************************
template <class T>
void MVops::apl_multifit_linear_free(apl_multifit_linear_workspace<T> *work) {
  RETURN_IF_NULL(work);
  apl_matrix_free<T>(work->A);
  apl_matrix_free<T>(work->Q);
  apl_matrix_free<T>(work->QSI);
  apl_vector_free<T>(work->S);
  apl_vector_free<T>(work->t);
  apl_vector_free<T>(work->xt);
  apl_vector_free<T>(work->D);
  free(work);
}
// ******************************************************
template <class T>
int MVops::apl_multifit_wlinear(const apl_matrix<T> *X,
                                const apl_vector<T> *w,
                                const apl_vector<T> *y,
                                apl_vector<T> *c,
                                apl_matrix<T> *cov,
                                double *chisq, apl_multifit_linear_workspace<T> *work) {
  size_t rank;
  int status = multifit_wlinear_svd<T>(X, w, y, APL_DBL_EPSILON, 1, &rank, c, cov, chisq, work);
  return status;
}
// ******************************************************
template <class T>
T MVops::apl_blas_dasum(const apl_vector<T> *v) {
  const int N = (int)v->size;
  const T *X = v->data;
  const int incX = (int)v->stride;
  T r = 0.0;
  int i;
  int ix = 0;

  if (incX <= 0) {
    return 0;
  }

  for (i = 0; i < N; i++) {
    r += fabs(X[ix]);
    ix += incX;
  }
  return r;
}
// ******************************************************
template <class T>
void MVops::apl_blas_dscal(T alpha, apl_vector<T> *v) {
  const int N = (int)v->size;
  T *X = v->data;
  const int incX = (int)v->stride;

  int i;
  int ix = 0;

  if (incX <= 0) {
    return;
  }

  for (i = 0; i < N; i++) {
    X[ix] *= alpha;
    ix += incX;
  }
}
// ******************************************************
template <class T>
T MVops::apl_blas_dnrm2(const apl_vector<T> *v) {
  const int N = (int)v->size;
  const T *X = v->data;
  const int incX = (int)v->stride;

  T scale = 0.0;
  T ssq = 1.0;
  int i;
  int ix = 0;

  if (N <= 0 || incX <= 0) {
    return 0;
  } else if (N == 1) {
    return fabs(X[0]);
  }

  for (i = 0; i < N; i++) {
    const T x = X[ix];

    if (x != 0.0) {
      const T ax = fabs(x);

      if (scale < ax) {
        ssq = 1.0 + ssq * (scale / ax) * (scale / ax);
        scale = ax;
      } else {
        ssq += (ax / scale) * (ax / scale);
      }
    }

    ix += incX;
  }

  return scale * sqrt(ssq);
}
// ******************************************************
template <class T>
size_t
MVops::apl_blas_idamax(const apl_vector<T> *X) {
  return cblas_idamax(int(X->size), X->data, int(X->stride));
}
// ******************************************************
size_t
MVops::cblas_idamax(const int N, const double *X, const int incX) {
  double max = 0.0;
  int ix = 0;
  int i;
  size_t result = 0;

  if (incX <= 0) {
    return 0;
  }

  for (i = 0; i < N; i++) {
    if (fabs(X[ix]) > max) {
      max = fabs(X[ix]);
      result = i;
    }
    ix += incX;
  }

  return result;
}
// ******************************************************
template <class T>
int MVops::apl_linalg_balance_columns(apl_matrix<T> *A, apl_vector<T> *D) {
  const size_t N = A->size2;
  size_t j;

  if (D->size != A->size2) {
    cerr << "length of D must match second dimension of A" << std::endl;
    exit(0);
  }

  apl_vector_set_all<T>(D, 1.0);

  for (j = 0; j < N; j++) {
    apl_vector_view<T> A_j = apl_matrix_column<T>(A, j);

    double s = apl_blas_dasum<double>(&A_j.vector);

    double f = 1.0;

    if ((_iszero(s)) || (!std::isfinite(s))) {
      apl_vector_set<T>(D, j, f);
      continue;
    }

    // FIXME: we could use frexp() here //

    while (s > 1.0) {
      s /= 2.0;
      f *= 2.0;
    }

    while (s < 0.5) {
      s *= 2.0;
      f /= 2.0;
    }

    apl_vector_set<T>(D, j, f);

    if (f != 1.0) {
      apl_blas_dscal<T>(1.0 / f, &A_j.vector);
    }
  }

  return 0;
}
// ******************************************************
template <class T>
T MVops::apl_linalg_householder_transform(apl_vector<T> *v) {
  // replace v[0:n-1] with a householder vector (v[0:n-1]) and
  //   coefficient tau that annihilate v[1:n-1] //

  const size_t n = v->size;

  if (n == 1) {
    return 0.0;  // tau = 0 //
  } else {
    double alpha, beta, tau;

    apl_vector_view<T> x = apl_vector_subvector<T>(v, 1, n - 1);

    double xnorm = apl_blas_dnrm2<T>(&x.vector);

    if (_iszero(xnorm)) {
      return 0.0;  // tau = 0 //
    }

    alpha = apl_vector_get<T>(v, 0);
    if (_iszero(alpha)) alpha = 0.00;
    beta = -(alpha >= 0.0 ? +1.0 : -1.0) * hypot(alpha, xnorm);
    tau = (beta - alpha) / beta;

    {
      double s = (alpha - beta);

      if (fabs(s) > APL_DBL_MIN) {
        apl_blas_dscal<T>(1.0 / s, &x.vector);
        apl_vector_set<T>(v, 0, beta);
      } else {
        apl_blas_dscal<T>(APL_DBL_EPSILON / s, &x.vector);
        apl_blas_dscal<T>(1.0 / APL_DBL_EPSILON, &x.vector);
        apl_vector_set<T>(v, 0, beta);
      }
    }

    return tau;
  }
}
// ******************************************************
template <class T>
int MVops::apl_linalg_householder_hm(T tau, const apl_vector<T> *v, apl_matrix<T> *A) {
  // applies a householder transformation v,tau to matrix m //

  if (_iszero(tau)) {
    return 0;
  }

  {
    size_t i, j;

    for (j = 0; j < A->size2; j++) {
      // Compute wj = Akj vk //

      double wj = apl_matrix_get<T>(A, 0, j);

      for (i = 1; i < A->size1; i++)  // note, computed for v(0) = 1 above //
      {
        wj += apl_matrix_get<T>(A, i, j) * apl_vector_get<T>(v, i);
      }

      // Aij = Aij - tau vi wj //

      // i = 0 //
      {
        double A0j = apl_matrix_get<T>(A, 0, j);
        apl_matrix_set<T>(A, 0, j, A0j - tau * wj);
      }

      // i = 1 .. M-1 //

      for (i = 1; i < A->size1; i++) {
        double Aij = apl_matrix_get<T>(A, i, j);
        double vi = apl_vector_get<T>(v, i);
        apl_matrix_set<T>(A, i, j, Aij - tau * vi * wj);
      }
    }
  }

  return 0;
}
// ******************************************************
template <class T>
int MVops::apl_linalg_householder_hm1(double tau, apl_matrix<T> *A) {
  //   applies a householder transformation v,tau to a matrix being
  //   build up from the identity matrix, using the first column of A as
  //   a householder vector

  if (_iszero(tau)) {
    size_t i, j;

    apl_matrix_set<T>(A, 0, 0, 1.0);

    for (j = 1; j < A->size2; j++) {
      apl_matrix_set<T>(A, 0, j, 0.0);
    }

    for (i = 1; i < A->size1; i++) {
      apl_matrix_set<T>(A, i, 0, 0.0);
    }

    return 0;
  }

  // w = A' v //
  {
    size_t i, j;

    for (j = 1; j < A->size2; j++) {
      double wj = 0.0;  // A0j * v0 //

      for (i = 1; i < A->size1; i++) {
        double vi = apl_matrix_get<T>(A, i, 0);
        wj += apl_matrix_get<T>(A, i, j) * vi;
      }

      // A = A - tau v w' //

      apl_matrix_set<T>(A, 0, j, -tau * wj);

      for (i = 1; i < A->size1; i++) {
        double vi = apl_matrix_get<T>(A, i, 0);
        double Aij = apl_matrix_get<T>(A, i, j);
        apl_matrix_set<T>(A, i, j, Aij - tau * vi * wj);
      }
    }

    for (i = 1; i < A->size1; i++) {
      double vi = apl_matrix_get<T>(A, i, 0);
      apl_matrix_set<T>(A, i, 0, -tau * vi);
    }

    apl_matrix_set<T>(A, 0, 0, 1.0 - tau);
  }

  return 0;
}
// ******************************************************
template <class T>
int MVops::apl_linalg_householder_mh(double tau, const apl_vector<T> *v, apl_matrix<T> *A) {
  //   applies a householder transformation v,tau to matrix m from the
  //   right hand side in order to zero out rows //

  if (_iszero(tau))
    return 0;

  // A = A - tau w v' //
  {
    size_t i, j;

    for (i = 0; i < A->size1; i++) {
      double wi = apl_matrix_get<T>(A, i, 0);

      for (j = 1; j < A->size2; j++)  // note, computed for v(0) = 1 above //
      {
        wi += apl_matrix_get<T>(A, i, j) * apl_vector_get<T>(v, j);
      }

      // j = 0 //

      {
        double Ai0 = apl_matrix_get<T>(A, i, 0);
        apl_matrix_set<T>(A, i, 0, Ai0 - tau * wi);
      }

      // j = 1 .. N-1 //

      for (j = 1; j < A->size2; j++) {
        double vj = apl_vector_get<T>(v, j);
        double Aij = apl_matrix_get<T>(A, i, j);
        apl_matrix_set<T>(A, i, j, Aij - tau * wi * vj);
      }
    }
  }
  return 0;
}
// ******************************************************
template <class T>
int MVops::apl_linalg_bidiag_decomp(apl_matrix<T> *A, apl_vector<T> *tau_U, apl_vector<T> *tau_V) {
  if (A->size1 < A->size2) {
    cerr << "bidiagonal decomposition requires M>=N" << std::endl;
    exit(0);
  } else if (tau_U->size != A->size2) {
    cerr << "size of tau_U must be N" << std::endl;
    exit(0);
  } else if (tau_V->size + 1 != A->size2) {
    cerr << "size of tau_V must be (N - 1)" << std::endl;
    exit(0);
  } else {
    const size_t M = A->size1;
    const size_t N = A->size2;
    size_t i;

    for (i = 0; i < N; i++) {
      // Apply Householder transformation to current column //

      {
        apl_vector_view<T> c = apl_matrix_column<T>(A, i);
        apl_vector_view<T> v = apl_vector_subvector<T>(&c.vector, i, M - i);
        double tau_i = apl_linalg_householder_transform<T>(&v.vector);

        // Apply the transformation to the remaining columns //

        if (i + 1 < N) {
          apl_matrix_view<T> m =
              apl_matrix_submatrix<T>(A, i, i + 1, M - i, N - (i + 1));
          apl_linalg_householder_hm<T>(tau_i, &v.vector, &m.matrix);
        }

        apl_vector_set<T>(tau_U, i, tau_i);
      }

      // Apply Householder transformation to current row //

      if (i + 1 < N) {
        apl_vector_view<T> r = apl_matrix_row<T>(A, i);
        apl_vector_view<T> v = apl_vector_subvector<T>(&r.vector, i + 1, N - (i + 1));
        double tau_i = apl_linalg_householder_transform<T>(&v.vector);

        // Apply the transformation to the remaining rows //

        if (i + 1 < M) {
          apl_matrix_view<T> m =
              apl_matrix_submatrix<T>(A, i + 1, i + 1, M - (i + 1), N - (i + 1));
          apl_linalg_householder_mh<T>(tau_i, &v.vector, &m.matrix);
        }

        apl_vector_set<T>(tau_V, i, tau_i);
      }
    }
  }

  return 0;
}
// ******************************************************
template <class T>
int MVops::apl_linalg_bidiag_unpack2(apl_matrix<T> *A,
                                     apl_vector<T> *tau_U,
                                     apl_vector<T> *tau_V,
                                     apl_matrix<T> *V) {
  const size_t M = A->size1;
  const size_t N = A->size2;

  const size_t K = std::min(M, N);

  if (M < N) {
    cerr << "matrix A must have M >= N" << std::endl;
    exit(0);
  } else if (tau_U->size != K) {
    cerr << "size of tau must be MIN(M,N)" << std::endl;
    exit(0);
  } else if (tau_V->size + 1 != K) {
    cerr << "size of tau must be MIN(M,N) - 1" << std::endl;
    exit(0);
  } else if (V->size1 != N || V->size2 != N) {
    cerr << "size of V must be N x N" << std::endl;
    exit(0);
  } else {
    size_t i, j;

    // Initialize V to the identity //

    apl_matrix_set_identity<T>(V);

    for (i = N - 1; i-- > 0;) {
      // Householder row transformation to accumulate V //
      const apl_vector_view<T> r = apl_matrix_const_row(A, i);
      const apl_vector_view<T> h =
          apl_vector_const_subvector<T>(&r.vector, i + 1, N - (i + 1));

      T ti = apl_vector_get<T>(tau_V, i);

      apl_matrix_view<T> m =
          apl_matrix_submatrix<T>(V, i + 1, i + 1, N - (i + 1), N - (i + 1));

      apl_linalg_householder_hm<T>(ti, &h.vector, &m.matrix);
    }

    // Copy superdiagonal into tau_v //

    for (i = 0; i < N - 1; i++) {
      T Aij = apl_matrix_get<T>(A, i, i + 1);
      apl_vector_set<T>(tau_V, i, Aij);
    }

    // Allow U to be unpacked into the same memory as A, copy
    //  diagonal into tau_U //

    for (j = N; j-- > 0;) {
      // Householder column transformation to accumulate U //
      T tj = apl_vector_get<T>(tau_U, j);
      T Ajj = apl_matrix_get<T>(A, j, j);
      apl_matrix_view<T> m = apl_matrix_submatrix<T>(A, j, j, M - j, N - j);

      apl_vector_set<T>(tau_U, j, Ajj);
      apl_linalg_householder_hm1<T>(tj, &m.matrix);
    }

    return 0;
  }
}
// ******************************************************
template <class T>
void MVops::chop_small_elements(apl_vector<T> *d, apl_vector<T> *f) {
  MVops fit;
  const size_t N = d->size;
  double d_i = fit.apl_vector_get(d, 0);

  size_t i;

  for (i = 0; i < N - 1; i++) {
    double f_i = fit.apl_vector_get<T>(f, i);
    double d_ip1 = fit.apl_vector_get<T>(d, i + 1);

    if (fabs(f_i) < APL_DBL_EPSILON * (fabs(d_i) + fabs(d_ip1))) {
      fit.apl_vector_set<T>(f, i, 0.0);
    }

    d_i = d_ip1;
  }
}
// ******************************************************
inline void
MVops::create_givens(const double a, const double b, double *c, double *s) {
  if (_iszero(b)) {
    *c = 1;
    *s = 0;
  } else if (fabs(b) > fabs(a)) {
    double t = -a / b;
    double s1 = 1.0 / sqrt(1 + t * t);
    *s = s1;
    *c = s1 * t;
  } else {
    double t = -b / a;
    double c1 = 1.0 / sqrt(1 + t * t);
    *c = c1;
    *s = c1 * t;
  }
}
// ******************************************************
template <class T>
void MVops::chase_out_trailing_zero(apl_vector<T> *d, apl_vector<T> *f, apl_matrix<T> *V) {
  MVops fit;
  const size_t N = V->size1;
  const size_t n = d->size;
  double c, s;
  double x, y;
  size_t k;

  x = fit.apl_vector_get<T>(d, n - 2);
  y = fit.apl_vector_get<T>(f, n - 2);

  for (k = n - 1; k-- > 0;) {
    fit.create_givens(x, y, &c, &s);

    // Compute V <= V G where G = [c, s ; -s, c] //

    {
      size_t i;

      for (i = 0; i < N; i++) {
        double Vip = fit.apl_matrix_get<T>(V, i, k);
        double Viq = fit.apl_matrix_get<T>(V, i, n - 1);
        fit.apl_matrix_set<T>(V, i, k, c * Vip - s * Viq);
        fit.apl_matrix_set<T>(V, i, n - 1, s * Vip + c * Viq);
      }
    }

    // compute B <= B G //

    fit.apl_vector_set<T>(d, k, c * x - s * y);

    if (k == n - 2)
      fit.apl_vector_set<T>(f, k, s * x + c * y);

    if (k > 0) {
      double z = fit.apl_vector_get<T>(f, k - 1);
      fit.apl_vector_set<T>(f, k - 1, c * z);

      x = fit.apl_vector_get<T>(d, k - 1);
      y = s * z;
    }
  }
}
// ******************************************************
void MVops::create_schur(double d0, double f0, double d1, double *c, double *s) {
  double apq = 2.0 * d0 * f0;

  if ((_iszero(d0)) || (_iszero(f0))) {
    *c = 1.0;
    *s = 0.0;
    return;
  }

  // Check if we need to rescale to avoid underflow/overflow //
  if (fabs(d0) < APL_SQRT_DBL_MIN || fabs(d0) > APL_SQRT_DBL_MAX || fabs(f0) < APL_SQRT_DBL_MIN || fabs(f0) > APL_SQRT_DBL_MAX || fabs(d1) < APL_SQRT_DBL_MIN || fabs(d1) > APL_SQRT_DBL_MAX) {
    double scale;
    int d0_exp, f0_exp;
    frexp(d0, &d0_exp);
    frexp(f0, &f0_exp);
    // Bring |d0*f0| into the range APL_DBL_MIN to APL_DBL_MAX //
    scale = ldexp(1.0, -(d0_exp + f0_exp) / 4);
    d0 *= scale;
    f0 *= scale;
    d1 *= scale;
    apq = 2.0 * d0 * f0;
  }

  if (apq != 0.0) {
    double t;
    double tau = (f0 * f0 + (d1 + d0) * (d1 - d0)) / apq;

    if (tau >= 0.0) {
      t = 1.0 / (tau + hypot(1.0, tau));
    } else {
      t = -1.0 / (-tau + hypot(1.0, tau));
    }

    *c = 1.0 / hypot(1.0, t);
    *s = t * (*c);
  } else {
    *c = 1.0;
    *s = 0.0;
  }
}
// ******************************************************
template <class T>
void MVops::svd2(apl_vector<T> *d, apl_vector<T> *f, apl_matrix<T> *U, apl_matrix<T> *V) {
  size_t i;
  double c, s, a11, a12, a21, a22;

  const size_t M = U->size1;
  const size_t N = V->size1;

  double d0 = apl_vector_get<T>(d, 0);
  double f0 = apl_vector_get<T>(f, 0);

  double d1 = apl_vector_get<T>(d, 1);

  if (_iszero(d0)) {
    // Eliminate off-diagonal element in [0,f0;0,d1] to make [d,0;0,0] //

    create_givens(f0, d1, &c, &s);

    // compute B <= G^T B X,  where X = [0,1;1,0] //

    apl_vector_set<T>(d, 0, c * f0 - s * d1);
    apl_vector_set<T>(f, 0, s * f0 + c * d1);
    apl_vector_set<T>(d, 1, 0.0);

    // Compute U <= U G //

    for (i = 0; i < M; i++) {
      double Uip = apl_matrix_get<T>(U, i, 0);
      double Uiq = apl_matrix_get<T>(U, i, 1);
      apl_matrix_set<T>(U, i, 0, c * Uip - s * Uiq);
      apl_matrix_set<T>(U, i, 1, s * Uip + c * Uiq);
    }

    // Compute V <= V X //

    apl_matrix_swap_columns<T>(V, 0, 1);

    return;
  } else if (_iszero(d1)) {
    // Eliminate off-diagonal element in [d0,f0;0,0] //

    create_givens(d0, f0, &c, &s);

    // compute B <= B G //

    apl_vector_set<T>(d, 0, d0 * c - f0 * s);
    apl_vector_set<T>(f, 0, 0.0);

    // Compute V <= V G //

    for (i = 0; i < N; i++) {
      double Vip = apl_matrix_get(V, i, 0);
      double Viq = apl_matrix_get<T>(V, i, 1);
      apl_matrix_set<T>(V, i, 0, c * Vip - s * Viq);
      apl_matrix_set<T>(V, i, 1, s * Vip + c * Viq);
    }

    return;
  } else {
    // Make columns orthogonal, A = [d0, f0; 0, d1] * G //

    create_schur(d0, f0, d1, &c, &s);

    // compute B <= B G //

    a11 = c * d0 - s * f0;
    a21 = -s * d1;

    a12 = s * d0 + c * f0;
    a22 = c * d1;

    // Compute V <= V G //

    for (i = 0; i < N; i++) {
      double Vip = apl_matrix_get<T>(V, i, 0);
      double Viq = apl_matrix_get<T>(V, i, 1);
      apl_matrix_set<T>(V, i, 0, c * Vip - s * Viq);
      apl_matrix_set<T>(V, i, 1, s * Vip + c * Viq);
    }

    // Eliminate off-diagonal elements, bring column with largest
    //   norm to first column

    if (hypot(a11, a21) < hypot(a12, a22)) {
      double t1, t2;

      // B <= B X //

      t1 = a11;
      a11 = a12;
      a12 = t1;
      t2 = a21;
      a21 = a22;
      a22 = t2;

      // V <= V X //

      apl_matrix_swap_columns<T>(V, 0, 1);
    }
    create_givens(a11, a21, &c, &s);

    // compute B <= G^T B //

    apl_vector_set<T>(d, 0, c * a11 - s * a21);
    apl_vector_set<T>(f, 0, c * a12 - s * a22);
    apl_vector_set<T>(d, 1, s * a12 + c * a22);

    // Compute U <= U G //

    for (i = 0; i < M; i++) {
      double Uip = apl_matrix_get<T>(U, i, 0);
      double Uiq = apl_matrix_get<T>(U, i, 1);
      apl_matrix_set<T>(U, i, 0, c * Uip - s * Uiq);
      apl_matrix_set<T>(U, i, 1, s * Uip + c * Uiq);
    }

    return;
  }
}
// ******************************************************
template <class T>
void MVops::chase_out_intermediate_zero(apl_vector<T> *d, apl_vector<T> *f, apl_matrix<T> *U, size_t k0) {
  MVops fit;
  const size_t M = U->size1;
  const size_t n = d->size;
  double c, s;
  double x, y;
  size_t k;

  x = fit.apl_vector_get<T>(f, k0);
  y = fit.apl_vector_get<T>(d, k0 + 1);

  for (k = k0; k < n - 1; k++) {
    fit.create_givens(y, -x, &c, &s);

    // Compute U <= U G //

    {
      size_t i;

      for (i = 0; i < M; i++) {
        double Uip = fit.apl_matrix_get<T>(U, i, k0);
        double Uiq = fit.apl_matrix_get<T>(U, i, k + 1);
        fit.apl_matrix_set<T>(U, i, k0, c * Uip - s * Uiq);
        fit.apl_matrix_set<T>(U, i, k + 1, s * Uip + c * Uiq);
      }
    }

    // compute B <= G^T B //

    fit.apl_vector_set<T>(d, k + 1, s * x + c * y);

    if (k == k0)
      fit.apl_vector_set<T>(f, k, c * x - s * y);

    if (k < n - 2) {
      double z = fit.apl_vector_get<T>(f, k + 1);
      fit.apl_vector_set<T>(f, k + 1, c * z);

      x = -s * z;
      y = fit.apl_vector_get<T>(d, k + 2);
    }
  }
}
// ******************************************************
template <class T>
double
MVops::trailing_eigenvalue(const apl_vector<T> *d, const apl_vector<T> *f) {
  MVops fit;
  const size_t n = d->size;

  double da = fit.apl_vector_get<T>(d, n - 2);
  double db = fit.apl_vector_get<T>(d, n - 1);
  double fa = (n > 2) ? fit.apl_vector_get<T>(f, n - 3) : 0.0;
  double fb = fit.apl_vector_get<T>(f, n - 2);

  double mu;

  {
    // We can compute mu more accurately than using the formula above
    //   since we know the roots cannot be negative.  This also avoids
    //   the possibility of NaNs in the formula above.

    //    The matrix is [ da^2 + fa^2,  da fb      ;
    //                    da fb      , db^2 + fb^2 ]
    //    and mu is the eigenvalue closest to the bottom right element.

    double ta = da * da + fa * fa;
    double tb = db * db + fb * fb;
    double tab = da * fb;

    double dt = (ta - tb) / 2.0;

    double S = ta + tb;
    double da2 = da * da, db2 = db * db;
    double fa2 = fa * fa, fb2 = fb * fb;
    double P = (da2 * db2) + (fa2 * db2) + (fa2 * fb2);
    double D = hypot(dt, tab);
    double r1 = S / 2 + D;

    if (dt >= 0) {
      // tb < ta, choose smaller root //
      mu = (r1 > 0) ? P / r1 : 0.0;
    } else {
      // tb > ta, choose larger root //
      mu = r1;
    }
  }

  return mu;
}
// ******************************************************
template <class T>
void MVops::qrstep(apl_vector<T> *d, apl_vector<T> *f, apl_matrix<T> *U, apl_matrix<T> *V) {
  MVops fit;
  const size_t M = U->size1;
  const size_t N = V->size1;
  const size_t n = d->size;
  double y, z;
  double ak, bk, zk, ap, bp, aq;  //bq;
  size_t i, k;

  if (n == 1)
    return;  // shouldn't happen //

  // Compute 2x2 svd directly //

  if (n == 2) {
    svd2<T>(d, f, U, V);
    return;
  }

  // Chase out any zeroes on the diagonal //

  for (i = 0; i < n - 1; i++) {
    double d_i = apl_vector_get<T>(d, i);

    if (_iszero(d_i)) {
      fit.chase_out_intermediate_zero<T>(d, f, U, i);
      return;
    }
  }
  // Chase out any zero at the end of the diagonal //

  {
    double d_nm1 = apl_vector_get<T>(d, n - 1);

    if (_iszero(d_nm1)) {
      fit.chase_out_trailing_zero<T>(d, f, V);
      return;
    }
  }
  // Apply QR reduction steps to the diagonal and offdiagonal //
  {
    double d0 = apl_vector_get<T>(d, 0);
    double f0 = apl_vector_get<T>(f, 0);

    double d1 = apl_vector_get<T>(d, 1);
    //double f1 = apl_vector_get<T> (f, 1);

    {
      double mu = fit.trailing_eigenvalue<T>(d, f);

      y = d0 * d0 - mu;
      z = d0 * f0;
    }

    // Set up the recurrence for Givens rotations on a bidiagonal matrix //

    ak = 0;
    bk = 0;

    ap = d0;
    bp = f0;

    aq = d1;
    //bq = f1;
  }

  for (k = 0; k < n - 1; k++) {
    double c, s;
    fit.create_givens(y, z, &c, &s);

    // Compute V <= V G //

    for (i = 0; i < N; i++) {
      double Vip = apl_matrix_get<T>(V, i, k);
      double Viq = apl_matrix_get<T>(V, i, k + 1);
      apl_matrix_set<T>(V, i, k, c * Vip - s * Viq);
      apl_matrix_set(V, i, k + 1, s * Vip + c * Viq);
    }

    // compute B <= B G //

    {
      double bk1 = c * bk - s * z;

      double ap1 = c * ap - s * bp;
      double bp1 = s * ap + c * bp;
      double zp1 = -s * aq;

      double aq1 = c * aq;

      if (k > 0) {
        apl_vector_set<T>(f, k - 1, bk1);
      }

      ak = ap1;
      bk = bp1;
      zk = zp1;

      ap = aq1;

      if (k < n - 2) {
        bp = apl_vector_get<T>(f, k + 1);
      } else {
        bp = 0.0;
      }

      y = ak;
      z = zk;
    }

    fit.create_givens(y, z, &c, &s);

    // Compute U <= U G //

    for (i = 0; i < M; i++) {
      double Uip = apl_matrix_get<T>(U, i, k);
      double Uiq = apl_matrix_get<T>(U, i, k + 1);
      apl_matrix_set<T>(U, i, k, c * Uip - s * Uiq);
      apl_matrix_set<T>(U, i, k + 1, s * Uip + c * Uiq);
    }

    // compute B <= G^T B //

    {
      double ak1 = c * ak - s * zk;
      double bk1 = c * bk - s * ap;
      double zk1 = -s * bp;

      double ap1 = s * bk + c * ap;
      double bp1 = c * bp;

      apl_vector_set<T>(d, k, ak1);

      ak = ak1;
      bk = bk1;
      zk = zk1;

      ap = ap1;
      bp = bp1;

      if (k < n - 2) {
        aq = apl_vector_get<T>(d, k + 2);
      } else {
        aq = 0.0;
      }

      y = bk;
      z = zk;
    }
  }
  apl_vector_set<T>(f, n - 2, bk);
  apl_vector_set<T>(d, n - 1, ap);
}
// ******************************************************
template <class T>
int MVops::apl_linalg_SV_decomp(apl_matrix<T> *A, apl_matrix<T> *V, apl_vector<T> *S,
                                apl_vector<T> *work) {
  size_t a, b, i, j, iter;

  const size_t M = A->size1;
  const size_t N = A->size2;
  const size_t K = std::min(M, N);

  if (M < N) {
    cerr << "svd of MxN matrix, M<N, is not implemented" << std::endl;
    exit(0);
  } else if (V->size1 != N) {
    cerr << "square matrix V must match second dimension of matrix A" << std::endl;
    exit(0);
  } else if (V->size1 != V->size2) {
    cerr << "matrix V must be square" << std::endl;
    exit(0);
  } else if (S->size != N) {
    cerr << "length of vector S must match second dimension of matrix A" << std::endl;
    exit(0);
  } else if (work->size != N) {
    cerr << "length of workspace must match second dimension of matrix A" << std::endl;
    exit(0);
  }

  // Handle the case of N = 1 (SVD of a column vector) //

  if (N == 1) {
    apl_vector_view<T> column = apl_matrix_column<T>(A, 0);
    double norm = apl_blas_dnrm2<T>(&column.vector);

    apl_vector_set<T>(S, 0, norm);
    apl_matrix_set<T>(V, 0, 0, 1.0);

    if (norm != 0.0) {
      apl_blas_dscal<T>(1.0 / norm, &column.vector);
    }

    return 0;
  }

  {
    apl_vector_view<T> f = apl_vector_subvector<T>(work, 0, K - 1);
    // bidiagonalize matrix A, unpack A into U S V //

    apl_linalg_bidiag_decomp<T>(A, S, &f.vector);
    apl_linalg_bidiag_unpack2<T>(A, S, &f.vector, V);
    // apply reduction steps to B=(S,Sd) //
    MVops fit;
    fit.chop_small_elements<T>(S, &f.vector);

    // Progressively reduce the matrix until it is diagonal //
    b = N - 1;
    iter = 0;

    while (b > 0) {
      double fbm1 = apl_vector_get<T>(&f.vector, b - 1);
      if ((_iszero(fbm1)) || (std::isnan(fbm1))) {
        b--;
        continue;
      }

      // Find the largest unreduced block (a,b) starting from b
      //   and working backwards

      a = b - 1;

      while (a > 0) {
        double fam1 = apl_vector_get<T>(&f.vector, a - 1);

        if ((_iszero(fam1)) || (std::isnan(fam1))) {
          break;
        }

        a--;
      }

      iter++;

      if (iter > 100 * N) {
        cerr << "SVD decomposition failed to converge" << std::endl;
        break;
      }

      {
        const size_t n_block = b - a + 1;
        apl_vector_view<T> S_block = apl_vector_subvector<T>(S, a, n_block);
        apl_vector_view<T> f_block = apl_vector_subvector<T>(&f.vector, a, n_block - 1);

        apl_matrix_view<T> U_block =
            apl_matrix_submatrix<T>(A, 0, a, A->size1, n_block);
        apl_matrix_view<T> V_block =
            apl_matrix_submatrix<T>(V, 0, a, V->size1, n_block);

        int rescale = 0;
        double scale = 1;
        double norm = 0;

        // Find the maximum absolute values of the diagonal and subdiagonal //

        for (i = 0; i < n_block; i++) {
          double s_i = apl_vector_get<T>(&S_block.vector, i);
          double a = fabs(s_i);
          if (a > norm) norm = a;
        }

        for (i = 0; i < n_block - 1; i++) {
          double f_i = apl_vector_get<T>(&f_block.vector, i);
          double a = fabs(f_i);
          if (a > norm) norm = a;
        }

        // Temporarily scale the submatrix if necessary //

        if (norm > APL_SQRT_DBL_MAX) {
          scale = (norm / APL_SQRT_DBL_MAX);
          rescale = 1;
        } else if (norm < APL_SQRT_DBL_MIN && norm > 0) {
          scale = (norm / APL_SQRT_DBL_MIN);
          rescale = 1;
        }

        if (rescale) {
          apl_blas_dscal<T>(1.0 / scale, &S_block.vector);
          apl_blas_dscal<T>(1.0 / scale, &f_block.vector);
        }

        // Perform the implicit QR step //

        fit.qrstep<T>(&S_block.vector, &f_block.vector, &U_block.matrix, &V_block.matrix);
        //remove any small off-diagonal elements

        chop_small_elements<T>(&S_block.vector, &f_block.vector);

        // Undo the scaling if needed //

        if (rescale) {
          apl_blas_dscal<T>(scale, &S_block.vector);
          apl_blas_dscal<T>(scale, &f_block.vector);
        }
      }
    }
  }
  // Make singular values positive by reflections if necessary //
  for (j = 0; j < K; j++) {
    double Sj = apl_vector_get<double>(S, j);

    if (Sj < 0.0) {
      for (i = 0; i < N; i++) {
        double Vij = apl_matrix_get<double>(V, i, j);
        apl_matrix_set<double>(V, i, j, -Vij);
      }

      apl_vector_set<double>(S, j, -Sj);
    }
  }

  // Sort singular values into decreasing order //

  for (i = 0; i < K; i++) {
    double S_max = apl_vector_get<double>(S, i);
    size_t i_max = i;

    for (j = i + 1; j < K; j++) {
      double Sj = apl_vector_get<double>(S, j);

      if (Sj > S_max) {
        S_max = Sj;
        i_max = j;
      }
    }

    if (i_max != i) {
      // swap eigenvalues //
      apl_vector_swap_elements<T>(S, i, i_max);

      // swap eigenvectors //
      apl_matrix_swap_columns<T>(A, i, i_max);
      apl_matrix_swap_columns<T>(V, i, i_max);
    }
  }
  return 0;
}
// ******************************************************
template <class T>
int MVops::apl_blas_daxpy(double alpha, const apl_vector<T> *u, apl_vector<T> *v) {
  if (u->size == v->size) {
    const int N = (int)u->size;
    const T *X = u->data;
    const int incX = (int)u->stride;
    T *Y = v->data;
    const int incY = v->stride;

    int i;

    if (_iszero(alpha)) {
      return 0;
    }

    if (incX == 1 && incY == 1) {
      const int m = N % 4;

      for (i = 0; i < m; i++) {
        Y[i] += alpha * X[i];
      }

      for (i = m; i + 3 < N; i += 4) {
        Y[i] += alpha * X[i];
        Y[i + 1] += alpha * X[i + 1];
        Y[i + 2] += alpha * X[i + 2];
        Y[i + 3] += alpha * X[i + 3];
      }
    } else {
      int ix = 0;  //OFFSET(N, incX);
      int iy = 0;  //OFFSET(N, incY);

      for (i = 0; i < N; i++) {
        Y[iy] += alpha * X[ix];
        ix += incX;
        iy += incY;
      }
    }
    return 0;
  } else {
    cerr << "invalid length in apl_blas_daxpy()" << std::endl;
    exit(0);
  }
}
// ******************************************************
void MVops::cblas_dgemv(const enum APL_ORDER order, const enum APL_TRANSPOSE TransA,
                        const int M, const int N, const double alpha, const double *A,
                        const int lda, const double *X, const int incX,
                        const double beta, double *Y, const int incY) {
  {
    int i, j;
    int lenX, lenY;

    const int Trans = (TransA != aplConjTrans) ? TransA : aplTrans;

    if (M == 0 || N == 0)
      return;

    if ((_iszero(alpha)) && (_isequal(beta, 1.0)))
      return;

    if (Trans == aplNoTrans) {
      lenX = N;
      lenY = M;
    } else {
      lenX = M;
      lenY = N;
    }

    // form  y := beta*y //
    if (_iszero(beta)) {
      int iy = 0;
      for (i = 0; i < lenY; i++) {
        Y[iy] = 0.0;
        iy += incY;
      }
    } else if (beta != 1.0) {
      int iy = 0;
      for (i = 0; i < lenY; i++) {
        Y[iy] *= beta;
        iy += incY;
      }
    }

    if (_iszero(alpha))
      return;

    if ((order == aplRowMajor && Trans == aplNoTrans) || (order == aplColMajor && Trans == aplTrans)) {
      // form  y := alpha*A*x + y //
      int iy = 0;
      for (i = 0; i < lenY; i++) {
        double temp = 0.0;
        int ix = 0;
        for (j = 0; j < lenX; j++) {
          temp += X[ix] * A[lda * i + j];
          ix += incX;
        }
        Y[iy] += alpha * temp;
        iy += incY;
      }
    } else if ((order == aplRowMajor && Trans == aplTrans) || (order == aplColMajor && Trans == aplNoTrans)) {
      // form  y := alpha*A'*x + y //
      int ix = 0;
      for (j = 0; j < lenX; j++) {
        const double temp = alpha * X[ix];
        if (temp != 0.0) {
          int iy = 0;
          for (i = 0; i < lenY; i++) {
            Y[iy] += temp * A[lda * j + i];
            iy += incY;
          }
        }
        ix += incX;
      }
    } else {
      cerr << "unrecognized operation" << std::endl;
      exit(0);
    }
  }
}
// ******************************************************
template <class T>
int MVops::apl_blas_dgemv(APL_TRANSPOSE_t TransA,
                          double alpha,
                          const apl_matrix<T> *A,
                          const apl_vector<T> *X,
                          double beta,
                          apl_vector<T> *Y) {
  const size_t M = A->size1;
  const size_t N = A->size2;

  if ((TransA == aplNoTrans && N == X->size && M == Y->size) || (TransA == aplTrans && M == X->size && N == Y->size)) {
    cblas_dgemv(aplRowMajor, TransA, int(M), int(N), alpha, A->data,
                int(A->tda), X->data, int(X->stride), beta, Y->data,
                int(Y->stride));
    return 0;
  } else {
    cerr << "invalid length" << std::endl;
    exit(0);
  }
}
// ******************************************************
template <class T>
int MVops::apl_blas_ddot(const apl_vector<T> *X, const apl_vector<T> *Y, T *result) {
  if (X->size == Y->size) {
    *result =
        cblas_ddot(int(X->size), X->data, int(X->stride), Y->data,
                   int(Y->stride));
    return 0;
  } else {
    cerr << "invalid length" << std::endl;
    exit(0);
  }
}
// ******************************************************
double
MVops::cblas_ddot(const int N, const double *X, const int incX, const double *Y,
                  const int incY) {
  double r = 0.0;
  int i;
  int ix = 0;
  int iy = 0;

  for (i = 0; i < N; i++) {
    r += X[ix] * Y[iy];
    ix += incX;
    iy += incY;
  }

  return r;
}
// ******************************************************
template <class T>
int MVops::apl_linalg_SV_decomp_mod(apl_matrix<T> *A,
                                    apl_matrix<T> *X,
                                    apl_matrix<T> *V, apl_vector<T> *S, apl_vector<T> *work) {
  size_t i, j;

  const size_t M = A->size1;
  const size_t N = A->size2;

  if (M < N) {
    cerr << "svd of MxN matrix, M<N, is not implemented" << std::endl;
    exit(0);
  } else if (V->size1 != N) {
    cerr << "square matrix V must match second dimension of matrix A" << std::endl;
    exit(0);
  } else if (V->size1 != V->size2) {
    cerr << "matrix V must be square" << std::endl;
    exit(0);
  } else if (X->size1 != N) {
    cerr << "square matrix X must match second dimension of matrix A" << std::endl;
    exit(0);
  } else if (X->size1 != X->size2) {
    cerr << "matrix X must be square" << std::endl;
    exit(0);
  } else if (S->size != N) {
    cerr << "length of vector S must match second dimension of matrix A" << std::endl;
    exit(0);
  } else if (work->size != N) {
    cerr << "length of workspace must match second dimension of matrix A" << std::endl;
    exit(0);
  }
  if (N == 1) {
    apl_vector_view<T> column = apl_matrix_column<T>(A, 0);
    double norm = apl_blas_dnrm2<T>(&column.vector);

    apl_vector_set<T>(S, 0, norm);
    apl_matrix_set<T>(V, 0, 0, 1.0);

    if (norm != 0.0) {
      apl_blas_dscal<T>(1.0 / norm, &column.vector);
    }

    return 0;
  }
  // Convert A into an upper triangular matrix R //

  for (i = 0; i < N; i++) {
    apl_vector_view<T> c = apl_matrix_column<T>(A, i);
    apl_vector_view<T> v = apl_vector_subvector<T>(&c.vector, i, M - i);
    double tau_i = apl_linalg_householder_transform<T>(&v.vector);

    // Apply the transformation to the remaining columns //

    if (i + 1 < N) {
      apl_matrix_view<T> m =
          apl_matrix_submatrix<T>(A, i, i + 1, M - i, N - (i + 1));
      apl_linalg_householder_hm<T>(tau_i, &v.vector, &m.matrix);
    }

    apl_vector_set<T>(S, i, tau_i);
  }

  // Copy the upper triangular part of A into X

  for (i = 0; i < N; i++) {
    for (j = 0; j < i; j++) {
      apl_matrix_set<T>(X, i, j, 0.0);
    }

    {
      double Aii = apl_matrix_get<T>(A, i, i);
      apl_matrix_set<T>(X, i, i, Aii);
    }

    for (j = i + 1; j < N; j++) {
      double Aij = apl_matrix_get<T>(A, i, j);
      apl_matrix_set<T>(X, i, j, Aij);
    }
  }

  // Convert A into an orthogonal matrix L //

  for (j = N; j-- > 0;) {
    // Householder column transformation to accumulate L //
    double tj = apl_vector_get<T>(S, j);
    apl_matrix_view<T> m = apl_matrix_submatrix<T>(A, j, j, M - j, N - j);
    apl_linalg_householder_hm1<T>(tj, &m.matrix);
  }

  // unpack R into X V S //
  apl_linalg_SV_decomp<T>(X, V, S, work);
  // Multiply L by X, to obtain U = L X, stored in U //

  {
    apl_vector_view<T> sum = apl_vector_subvector<T>(work, 0, N);

    for (i = 0; i < M; i++) {
      apl_vector_view<T> L_i = apl_matrix_row<T>(A, i);
      apl_vector_set_all<T>(&sum.vector, 0.0);

      for (j = 0; j < N; j++) {
        double Lij = apl_vector_get<T>(&L_i.vector, j);
        apl_vector_view<T> X_j = apl_matrix_row<T>(X, j);
        apl_blas_daxpy<T>(Lij, &X_j.vector, &sum.vector);
      }

      apl_vector_memcpy<T>(&L_i.vector, &sum.vector);
    }
  }
  return 0;
}
// ******************************************************
template <class T>
int MVops::multifit_wlinear_svd(const apl_matrix<T> *X,
                                const apl_vector<T> *w,
                                const apl_vector<T> *y,
                                double tol,
                                int balance,
                                size_t *rank,
                                apl_vector<T> *c,
                                apl_matrix<T> *cov,
                                double *chisq, apl_multifit_linear_workspace<T> *work) {
  if (X->size1 != y->size) {
    cerr << "number of observations in y does not match rows of matrix X" << std::endl;
    exit(0);
  } else if (X->size2 != c->size) {
    cerr << "number of parameters c does not match columns of matrix X" << std::endl;
    exit(0);
  } else if (w->size != y->size) {
    cerr << "number of weights does not match number of observations" << std::endl;
    exit(0);
  } else if (cov->size1 != cov->size2) {
    cerr << "covariance matrix is not square" << std::endl;
    exit(0);
  } else if (c->size != cov->size1) {
    cerr << "number of parameters does not match size of covariance matrix" << std::endl;
    exit(0);
  } else if (X->size1 != work->n || X->size2 != work->p) {
    cerr << "size of workspace does not match size of observation matrix" << std::endl;
    exit(0);
  } else {
    MVops fit;
    const size_t n = X->size1;
    const size_t p = X->size2;

    size_t i, j, p_eff;

    apl_matrix<T> *A = work->A;
    apl_matrix<T> *Q = work->Q;
    apl_matrix<T> *QSI = work->QSI;
    apl_vector<T> *S = work->S;
    apl_vector<T> *t = work->t;
    apl_vector<T> *xt = work->xt;
    apl_vector<T> *D = work->D;

    // Scale X,  A = sqrt(w) X //
    fit.apl_matrix_memcpy<T>(A, X);

    for (i = 0; i < n; i++) {
      double wi = fit.apl_vector_get<T>(w, i);

      if (wi < 0)
        wi = 0;

      {
        apl_vector_view<T> row = fit.apl_matrix_row<T>(A, i);
        fit.apl_vector_scale<T>(&row.vector, sqrt(wi));
      }
    }
    // Balance the columns of the matrix A if requested //

    if (balance) {
      fit.apl_linalg_balance_columns<T>(A, D);
    } else {
      fit.apl_vector_set_all<T>(D, 1.0);
    }

    // Decompose A into U S Q^T //

    fit.apl_linalg_SV_decomp_mod<T>(A, QSI, Q, S, xt);

    // Solve sqrt(w) y = A c for c, by first computing t = sqrt(w) y //
    for (i = 0; i < n; i++) {
      double wi = fit.apl_vector_get<double>(w, i);
      double yi = fit.apl_vector_get<double>(y, i);
      if (wi < 0)
        wi = 0;
      fit.apl_vector_set<double>(t, i, sqrt(wi) * yi);
    }
    fit.apl_blas_dgemv<T>(aplTrans, 1.0, A, t, 0.0, xt);

    // Scale the matrix Q,  Q' = Q S^-1 //

    fit.apl_matrix_memcpy<T>(QSI, Q);

    {
      double alpha0 = fit.apl_vector_get<T>(S, 0);
      p_eff = 0;

      for (j = 0; j < p; j++) {
        apl_vector_view<T> column = fit.apl_matrix_column<T>(QSI, j);
        double alpha = fit.apl_vector_get<T>(S, j);

        if (alpha <= tol * alpha0) {
          alpha = 0.0;
        } else {
          alpha = 1.0 / alpha;
          p_eff++;
        }

        fit.apl_vector_scale<T>(&column.vector, alpha);
      }

      *rank = p_eff;
    }

    fit.apl_vector_set_all<T>(c, 0.0);

    // Solution //
    fit.apl_blas_dgemv<T>(aplNoTrans, 1.0, QSI, xt, 0.0, c);

    // Unscale the balancing factors //

    fit.apl_vector_div<T>(c, D);

    // Compute chisq, from residual r = y - X c //

    {
      double r2 = 0;

      for (i = 0; i < n; i++) {
        double yi = fit.apl_vector_get<T>(y, i);
        double wi = fit.apl_vector_get<T>(w, i);
        const apl_vector_view<T> row = fit.apl_matrix_const_row<T>(X, i);

        double y_est, ri;
        fit.apl_blas_ddot<T>(&row.vector, c, &y_est);
        ri = yi - y_est;
        r2 += wi * ri * ri;
      }

      *chisq = r2;
    }

    // Form covariance matrix cov = (X^T W X)^-1 = (Q S^-1) (Q S^-1)^T //

    for (i = 0; i < p; i++) {
      apl_vector_view<T> row_i = fit.apl_matrix_row<T>(QSI, i);
      double d_i = fit.apl_vector_get<T>(D, i);

      for (j = i; j < p; j++) {
        apl_vector_view<T> row_j = fit.apl_matrix_row<T>(QSI, j);
        double d_j = fit.apl_vector_get<T>(D, j);
        double s;

        fit.apl_blas_ddot<T>(&row_i.vector, &row_j.vector, &s);

        fit.apl_matrix_set<T>(cov, i, j, s / (d_i * d_j));
        fit.apl_matrix_set<T>(cov, j, i, s / (d_i * d_j));
      }
    }
    return 0;
  }
}
// ******************************************************
// NL fit
template <class utype>
apl_multifit_fdfsolver<utype> *
MVops::apl_multifit_fdfsolver_alloc(const apl_multifit_fdfsolver_type<utype> *T,
                                    size_t n, size_t p) {
  int status;
  apl_multifit_fdfsolver<utype> *s;

  if (n < p) {
    cerr << "insufficient data points, n < p" << std::endl;
    exit(0);
  }

  s = (apl_multifit_fdfsolver<utype> *)malloc(sizeof(apl_multifit_fdfsolver<utype>));
  if (s == 0) {
    cerr << "failed to allocate space for multifit solver struct" << std::endl;
    exit(0);
  }

  s->x = apl_vector_calloc<utype>(p);

  if (s->x == 0) {
    free(s);
    cerr << "failed to allocate space for x" << std::endl;
    exit(0);
  }

  s->f = apl_vector_calloc<utype>(n);

  if (s->f == 0) {
    apl_vector_free<utype>(s->x);
    free(s);
    cerr << "failed to allocate space for f" << std::endl;
    exit(0);
  }

  s->J = apl_matrix_calloc<utype>(n, p);

  if (s->J == 0) {
    apl_vector_free<utype>(s->x);
    apl_vector_free<utype>(s->f);
    free(s);
    cerr << "failed to allocate space for g" << std::endl;
    exit(0);
  }

  s->dx = apl_vector_calloc<utype>(p);

  if (s->dx == 0) {
    apl_matrix_free<utype>(s->J);
    apl_vector_free<utype>(s->x);
    apl_vector_free<utype>(s->f);
    free(s);
    cerr << "failed to allocate space for dx" << std::endl;
    exit(0);
  }

  s->state = malloc(T->size);

  if (s->state == 0) {
    apl_vector_free<utype>(s->dx);
    apl_vector_free<utype>(s->x);
    apl_vector_free<utype>(s->f);
    apl_matrix_free<utype>(s->J);
    free(s);  // exception in constructor, avoid memory leak //

    cerr << "failed to allocate space for multifit solver state" << std::endl;
    exit(0);
  }
  s->type = T;

  status = (s->type->alloc)(s->state, n, p);

  if (status != 0) {
    free(s->state);
    apl_vector_free<utype>(s->dx);
    apl_vector_free<utype>(s->x);
    apl_vector_free<utype>(s->f);
    apl_matrix_free<utype>(s->J);
    free(s);  // exception in constructor, avoid memory leak //

    cerr << "failed to set solver" << std::endl;
    exit(0);
  }

  s->fdf = NULL;
  return s;
}
// ******************************************************
template <class utype>
int MVops::apl_multifit_fdfsolver_set(apl_multifit_fdfsolver<utype> *s,
                                      apl_multifit_function_fdf *f,
                                      const apl_vector<utype> *x) {
  if (s->f->size != f->n) {
    cerr << "function size does not match solver" << std::endl;
    exit(0);
  }

  if (s->x->size != x->size) {
    cerr << "vector length does not match solver" << std::endl;
    exit(0);
  }

  s->fdf = f;
  apl_vector_memcpy(s->x, x);

  return (s->type->set)(s->state, s->fdf, s->x, s->f, s->J, s->dx);
}
// ******************************************************
template <class utype>
int MVops::apl_multifit_fdfsolver_iterate(apl_multifit_fdfsolver<utype> *s) {
  return (s->type->iterate)(s->state, s->fdf, s->x, s->f, s->J, s->dx);
}
// ******************************************************
template <class utype>
int MVops::apl_multifit_test_delta(const apl_vector<utype> *dx, const apl_vector<utype> *x,
                                   double epsabs, double epsrel) {
  size_t i;
  int ok = 1;
  const size_t n = x->size;

  if (epsrel < 0.0) {
    cerr << "relative tolerance is negative" << std::endl;
    exit(0);
  }

  for (i = 0; i < n; i++) {
    double xi = apl_vector_get<utype>(x, i);
    double dxi = apl_vector_get<utype>(dx, i);
    double tolerance = epsabs + epsrel * fabs(xi);

    if (fabs(dxi) < tolerance) {
      ok = 1;
    } else {
      ok = 0;
      break;
    }
  }

  if (ok)
    return 0;

  return -2;
}
// ******************************************************
apl_permutation *
MVops::apl_permutation_alloc(const size_t n) {
  apl_permutation *p;

  if (n == 0) {
    cerr << "permutation length n must be positive integer" << std::endl;
    exit(0);
  }

  p = (apl_permutation *)malloc(sizeof(apl_permutation));

  if (p == 0) {
    cerr << "failed to allocate space for permutation struct" << std::endl;
    exit(0);
  }

  p->data = (size_t *)malloc(n * sizeof(size_t));

  if (p->data == 0) {
    free(p);  // exception in constructor, avoid memory leak //

    cerr << "failed to allocate space for permutation data" << std::endl;
    exit(0);
  }

  p->size = n;

  return p;
}
// ******************************************************
apl_permutation *
MVops::apl_permutation_calloc(const size_t n) {
  size_t i;

  apl_permutation *p = apl_permutation_alloc(n);

  if (p == 0)
    return 0;

  // initialize permutation to identity //

  for (i = 0; i < n; i++) {
    p->data[i] = i;
  }

  return p;
}
// ******************************************************
inline size_t
MVops::apl_permutation_get(const apl_permutation *p, const size_t i) {
#if APL_RANGE_CHECK
  if (i >= p->size) {
    cerr << "index out of range" << std::endl;
    exit(0);
  }
#endif
  return p->data[i];
}
// ******************************************************
void MVops::apl_permutation_free(apl_permutation *p) {
  RETURN_IF_NULL(p);
  free(p->data);
  free(p);
}
// ******************************************************
void MVops::apl_permutation_init(apl_permutation *p) {
  const size_t n = p->size;
  size_t i;

  // initialize permutation to identity //

  for (i = 0; i < n; i++) {
    p->data[i] = i;
  }
}
// ******************************************************
int MVops::apl_permutation_swap(apl_permutation *p, const size_t i, const size_t j) {
  const size_t size = p->size;

  if (i >= size) {
    cerr << "first index is out of range" << std::endl;
    exit(0);
  }

  if (j >= size) {
    cerr << "second index is out of range" << std::endl;
    exit(0);
  }

  if (i != j) {
    size_t tmp = p->data[i];
    p->data[i] = p->data[j];
    p->data[j] = tmp;
  }

  return 0;
}
// ******************************************************
template <class utype>
int MVops::apl_linalg_QRPT_decomp(apl_matrix<utype> *A, apl_vector<utype> *tau, apl_permutation *p, int *signum, apl_vector<utype> *norm) {
  const size_t M = A->size1;
  const size_t N = A->size2;

  if (tau->size != std::min(M, N)) {
    cerr << "size of tau must be MIN(M,N)" << std::endl;
    exit(0);
  } else if (p->size != N) {
    cerr << "permutation size must be N" << std::endl;
    exit(0);
  } else if (norm->size != N) {
    cerr << "norm size must be N" << std::endl;
    exit(0);
  } else {
    size_t i;

    *signum = 1;

    apl_permutation_init(p);  // set to identity //

    // Compute column norms and store in workspace //

    for (i = 0; i < N; i++) {
      apl_vector_view<utype> c = apl_matrix_column<utype>(A, i);
      double x = apl_blas_dnrm2<utype>(&c.vector);
      apl_vector_set<utype>(norm, i, x);
    }

    for (i = 0; i < std::min(M, N); i++) {
      // Bring the column of largest norm into the pivot position //

      double max_norm = apl_vector_get<utype>(norm, i);
      size_t j, kmax = i;

      for (j = i + 1; j < N; j++) {
        double x = apl_vector_get<utype>(norm, j);

        if (x > max_norm) {
          max_norm = x;
          kmax = j;
        }
      }

      if (kmax != i) {
        apl_matrix_swap_columns<utype>(A, i, kmax);
        apl_permutation_swap(p, i, kmax);
        apl_vector_swap_elements<utype>(norm, i, kmax);

        (*signum) = -(*signum);
      }

      // Compute the Householder transformation to reduce the j-th
      //   column of the matrix to a multiple of the j-th unit vector //

      {
        apl_vector_view<utype> c_full = apl_matrix_column<utype>(A, i);
        apl_vector_view<utype> c = apl_vector_subvector<utype>(&c_full.vector,
                                                               i, M - i);
        double tau_i = apl_linalg_householder_transform<utype>(&c.vector);

        apl_vector_set<utype>(tau, i, tau_i);

        // Apply the transformation to the remaining columns //

        if (i + 1 < N) {
          apl_matrix_view<utype> m = apl_matrix_submatrix<utype>(A, i, i + 1, M - i, N - (i + 1));

          apl_linalg_householder_hm<utype>(tau_i, &c.vector, &m.matrix);
        }
      }

      // Update the norms of the remaining columns too //

      if (i + 1 < M) {
        for (j = i + 1; j < N; j++) {
          double x = apl_vector_get<utype>(norm, j);

          if (x > 0.0) {
            double y = 0;
            double temp = apl_matrix_get<utype>(A, i, j) / x;

            if (fabs(temp) >= 1)
              y = 0.0;
            else
              y = x * sqrt(1 - temp * temp);

            // recompute norm to prevent loss of accuracy //

            if (fabs(y / x) < sqrt(20.0) * APL_SQRT_DBL_EPSILON) {
              apl_vector_view<utype> c_full = apl_matrix_column<utype>(A, j);
              apl_vector_view<utype> c =
                  apl_vector_subvector<utype>(&c_full.vector,
                                              i + 1, M - (i + 1));
              y = apl_blas_dnrm2<utype>(&c.vector);
            }

            apl_vector_set<utype>(norm, j, y);
          }
        }
      }
    }

    return 0;
  }
}
// ******************************************************
template <class utype>
int MVops::apl_multifit_covar(const apl_matrix<utype> *J, double epsrel, apl_matrix<utype> *covar) {
  double tolr;

  size_t i, j, k;
  size_t kmax = 0;

  apl_matrix<utype> *r;
  apl_vector<utype> *tau;
  apl_vector<utype> *norm;
  apl_permutation *perm;

  size_t m = J->size1, n = J->size2;

  if (m < n) {
    cerr << "Jacobian be rectangular M x N with M >= N" << std::endl;
    exit(0);
  }

  if (covar->size1 != covar->size2 || covar->size1 != n) {
    cerr << "covariance matrix must be square and match second dimension of jacobian" << std::endl;
    exit(0);
  }

  r = apl_matrix_alloc<utype>(m, n);
  tau = apl_vector_alloc<utype>(n);
  perm = apl_permutation_alloc(n);
  norm = apl_vector_alloc<utype>(n);

  {
    int signum = 0;
    apl_matrix_memcpy(r, J);
    apl_linalg_QRPT_decomp(r, tau, perm, &signum, norm);
  }

  // Form the inverse of R in the full upper triangle of R //

  tolr = epsrel * fabs(apl_matrix_get<utype>(r, 0, 0));

  for (k = 0; k < n; k++) {
    double rkk = apl_matrix_get<utype>(r, k, k);

    if (fabs(rkk) <= tolr) {
      break;
    }

    apl_matrix_set<utype>(r, k, k, 1.0 / rkk);

    for (j = 0; j < k; j++) {
      double t = apl_matrix_get<utype>(r, j, k) / rkk;
      apl_matrix_set<utype>(r, j, k, 0.0);

      for (i = 0; i <= j; i++) {
        double rik = apl_matrix_get<utype>(r, i, k);
        double rij = apl_matrix_get<utype>(r, i, j);

        apl_matrix_set<utype>(r, i, k, rik - t * rij);
      }
    }
    kmax = k;
  }

  // Form the full upper triangle of the inverse of R^T R in the full
  //   upper triangle of R

  for (k = 0; k <= kmax; k++) {
    for (j = 0; j < k; j++) {
      double rjk = apl_matrix_get<utype>(r, j, k);

      for (i = 0; i <= j; i++) {
        double rij = apl_matrix_get<utype>(r, i, j);
        double rik = apl_matrix_get<utype>(r, i, k);

        apl_matrix_set<utype>(r, i, j, rij + rjk * rik);
      }
    }

    {
      double t = apl_matrix_get<utype>(r, k, k);

      for (i = 0; i <= k; i++) {
        double rik = apl_matrix_get<utype>(r, i, k);

        apl_matrix_set<utype>(r, i, k, t * rik);
      };
    }
  }

  // Form the full lower triangle of the covariance matrix in the
  //   strict lower triangle of R and in w

  for (j = 0; j < n; j++) {
    size_t pj = apl_permutation_get(perm, j);

    for (i = 0; i <= j; i++) {
      size_t perm_i = apl_permutation_get(perm, i);

      double rij;

      if (j > kmax) {
        apl_matrix_set<utype>(r, i, j, 0.0);
        rij = 0.0;
      } else {
        rij = apl_matrix_get<utype>(r, i, j);
      }

      if (perm_i > pj) {
        apl_matrix_set<utype>(r, perm_i, pj, rij);
      } else if (perm_i < pj) {
        apl_matrix_set<utype>(r, pj, perm_i, rij);
      }
    }

    {
      double rjj = apl_matrix_get<utype>(r, j, j);
      apl_matrix_set<utype>(covar, pj, pj, rjj);
    }
  }

  // symmetrize the covariance matrix //

  for (j = 0; j < n; j++) {
    for (i = 0; i < j; i++) {
      double rji = apl_matrix_get<utype>(r, j, i);

      apl_matrix_set<utype>(covar, j, i, rji);
      apl_matrix_set<utype>(covar, i, j, rji);
    }
  }

  apl_matrix_free<utype>(r);
  apl_permutation_free(perm);
  apl_vector_free<utype>(tau);
  apl_vector_free<utype>(norm);
  return 0;
}
// ******************************************************
template <class utype>
int MVops::apl_linalg_householder_hv(double tau, const apl_vector<utype> *v, apl_vector<utype> *w) {
  // applies a householder transformation v to vector w //
  const size_t N = v->size;

  if (_iszero(tau))
    return 0;

  {
    // compute d = v'w //

    double d0 = apl_vector_get<utype>(w, 0);
    double d1, d;

    const apl_vector_view<utype> v1 = apl_vector_const_subvector<utype>(v, 1, N - 1);
    apl_vector_view<utype> w1 = apl_vector_subvector<utype>(w, 1, N - 1);

    apl_blas_ddot<utype>(&v1.vector, &w1.vector, &d1);

    d = d0 + d1;

    // compute w = w - tau (v) (v'w) //

    {
      double w0 = apl_vector_get<utype>(w, 0);
      apl_vector_set<utype>(w, 0, w0 - tau * d);
    }

    apl_blas_daxpy<utype>(-tau * d, &v1.vector, &w1.vector);
  }

  return 0;
}
// ******************************************************
template <class utype>
int MVops::apl_linalg_QR_QTvec(const apl_matrix<utype> *QR, const apl_vector<utype> *tau, apl_vector<utype> *v) {
  const size_t M = QR->size1;
  const size_t N = QR->size2;

  if (tau->size != std::min(M, N)) {
    cerr << "size of tau must be MIN(M,N)" << std::endl;
    exit(0);
  } else if (v->size != M) {
    cerr << "vector size must be M" << std::endl;
    exit(0);
  } else {
    size_t i;

    // compute Q^T v //

    for (i = 0; i < std::min(M, N); i++) {
      const apl_vector_view<utype> c = apl_matrix_const_column<utype>(QR, i);
      const apl_vector_view<utype> h = apl_vector_const_subvector<utype>(&(c.vector), i, M - i);
      apl_vector_view<utype> w = apl_vector_subvector<utype>(v, i, M - i);
      double ti = apl_vector_get<utype>(tau, i);
      apl_linalg_householder_hv<utype>(ti, &(h.vector), &(w.vector));
    }
    return 0;
  }
}
// ******************************************************
template <class utype>
void MVops::compute_gradient_direction(const apl_matrix<utype> *r, const apl_permutation *p,
                                       const apl_vector<utype> *qtf, const apl_vector<utype> *diag,
                                       apl_vector<utype> *g) {
  MVops fit;
  const size_t n = r->size2;

  size_t i, j;

  for (j = 0; j < n; j++) {
    double sum = 0;

    for (i = 0; i <= j; i++) {
      sum += fit.apl_matrix_get<utype>(r, i, j) * fit.apl_vector_get<utype>(qtf, i);
    }

    {
      size_t pj = fit.apl_permutation_get(p, j);
      double dpj = fit.apl_vector_get<utype>(diag, pj);

      fit.apl_vector_set<utype>(g, j, sum / dpj);
    }
  }
}
// ******************************************************
template <class utype>
int MVops::apl_permute_vector_inverse(const apl_permutation *p, apl_vector<utype> *v) {
  if (v->size != p->size) {
    cerr << "vector and permutation must be the same length" << std::endl;
    exit(0);
  }

  apl_permute_inverse(p->data, v->data, v->stride, v->size);

  return 0;
}
// ******************************************************
template <class utype>
void MVops::compute_newton_direction(const apl_matrix<utype> *r, const apl_permutation *perm,
                                     const apl_vector<utype> *qtf, apl_vector<utype> *x) {
  // Compute and store in x the Gauss-Newton direction. If the
  //   Jacobian is rank-deficient then obtain a least squares
  //   solution.
  MVops fit;
  const size_t n = r->size2;
  size_t i, j, nsing;

  for (i = 0; i < n; i++) {
    double qtfi = fit.apl_vector_get<utype>(qtf, i);
    fit.apl_vector_set<utype>(x, i, qtfi);
  }

  nsing = fit.count_nsing(r);

  for (i = nsing; i < n; i++) {
    fit.apl_vector_set<utype>(x, i, 0.0);
  }

  if (nsing > 0) {
    for (j = nsing; j > 0 && j--;) {
      double rjj = fit.apl_matrix_get<utype>(r, j, j);
      double temp = fit.apl_vector_get<utype>(x, j) / rjj;

      fit.apl_vector_set<utype>(x, j, temp);

      for (i = 0; i < j; i++) {
        double rij = fit.apl_matrix_get<utype>(r, i, j);
        double xi = fit.apl_vector_get<utype>(x, i);
        fit.apl_vector_set<utype>(x, i, xi - rij * temp);
      }
    }
  }

  fit.apl_permute_vector_inverse<utype>(perm, x);
}
// ******************************************************
template <class utype>
size_t
MVops::count_nsing(const apl_matrix<utype> *r) {
  // Count the number of nonsingular entries. Returns the index of the
  //   first entry which is singular.
  MVops fit;
  size_t n = r->size2;
  size_t i;

  for (i = 0; i < n; i++) {
    double rii = fit.apl_matrix_get<utype>(r, i, i);

    if (_iszero(rii)) {
      break;
    }
  }

  return i;
}
// ******************************************************
double
MVops::scaled_enorm(const apl_vector<double> *d, const apl_vector<double> *f) {
  MVops fit;
  double e2 = 0;
  size_t i, n = f->size;
  for (i = 0; i < n; i++) {
    double fi = fit.apl_vector_get<double>(f, i);
    double di = fit.apl_vector_get<double>(d, i);
    double u = di * fi;
    e2 += u * u;
  }
  return sqrt(e2);
}
// ******************************************************
template <class utype>
void MVops::compute_newton_bound(const apl_matrix<utype> *r, const apl_vector<utype> *x,
                                 double dxnorm, const apl_permutation *perm,
                                 const apl_vector<utype> *diag, apl_vector<utype> *w) {
  // If the jacobian is not rank-deficient then the Newton step
  //   provides a lower bound for the zero of the function. Otherwise
  //   set this bound to zero.
  MVops fit;

  size_t n = r->size2;

  size_t i, j;

  size_t nsing = fit.count_nsing(r);

  if (nsing < n) {
    fit.apl_vector_set_all<utype>(w, 0.0);
    return;
  }

  for (i = 0; i < n; i++) {
    size_t perm_i = fit.apl_permutation_get(perm, i);

    double dpi = fit.apl_vector_get<utype>(diag, perm_i);
    double xpi = fit.apl_vector_get<utype>(x, perm_i);

    fit.apl_vector_set<utype>(w, i, dpi * (dpi * xpi / dxnorm));
  }

  for (j = 0; j < n; j++) {
    double sum = 0;

    for (i = 0; i < j; i++) {
      sum += fit.apl_matrix_get<utype>(r, i, j) * fit.apl_vector_get<utype>(w, i);
    }

    {
      double rjj = fit.apl_matrix_get<utype>(r, j, j);
      double wj = fit.apl_vector_get<utype>(w, j);

      fit.apl_vector_set<utype>(w, j, (wj - sum) / rjj);
    }
  }
}
// ******************************************************
template <class utype>
double
MVops::enorm(const apl_vector<utype> *f) {
  MVops fit;
  double e2 = 0;
  size_t i, n = f->size;
  for (i = 0; i < n; i++) {
    double fi = fit.apl_vector_get<utype>(f, i);
    e2 += fi * fi;
  }
  return sqrt(e2);
}
// ******************************************************
template <class utype>
int MVops::qrsolv(apl_matrix<utype> *r, const apl_permutation *p, const double lambda,
                  const apl_vector<utype> *diag, const apl_vector<utype> *qtb,
                  apl_vector<utype> *x, apl_vector<utype> *sdiag, apl_vector<utype> *wa) {
  MVops fit;
  size_t n = r->size2;

  size_t i, j, k, nsing;

  // Copy r and qtb to preserve input and initialise s. In particular,
  //   save the diagonal elements of r in x

  for (j = 0; j < n; j++) {
    double rjj = fit.apl_matrix_get<utype>(r, j, j);
    double qtbj = fit.apl_vector_get<utype>(qtb, j);

    for (i = j + 1; i < n; i++) {
      double rji = fit.apl_matrix_get<utype>(r, j, i);
      fit.apl_matrix_set<utype>(r, i, j, rji);
    }

    fit.apl_vector_set<utype>(x, j, rjj);
    fit.apl_vector_set<utype>(wa, j, qtbj);
  }

  // Eliminate the diagonal matrix d using a Givens rotation //

  for (j = 0; j < n; j++) {
    double qtbpj;

    size_t pj = fit.apl_permutation_get(p, j);

    double diagpj = lambda * fit.apl_vector_get<utype>(diag, pj);

    if (_iszero(diagpj)) {
      continue;
    }

    fit.apl_vector_set<utype>(sdiag, j, diagpj);

    for (k = j + 1; k < n; k++) {
      fit.apl_vector_set<utype>(sdiag, k, 0.0);
    }

    // The transformations to eliminate the row of d modify only a
    //   single element of qtb beyond the first n, which is initially
    //   zero

    qtbpj = 0;

    for (k = j; k < n; k++) {
      // Determine a Givens rotation which eliminates the
      //   appropriate element in the current row of d

      double sine, cosine;

      double wak = fit.apl_vector_get<utype>(wa, k);
      double rkk = fit.apl_matrix_get<utype>(r, k, k);
      double sdiagk = fit.apl_vector_get<utype>(sdiag, k);

      if (_iszero(sdiagk)) {
        continue;
      }

      if (fabs(rkk) < fabs(sdiagk)) {
        double cotangent = rkk / sdiagk;
        sine = 0.5 / sqrt(0.25 + 0.25 * cotangent * cotangent);
        cosine = sine * cotangent;
      } else {
        double tangent = sdiagk / rkk;
        cosine = 0.5 / sqrt(0.25 + 0.25 * tangent * tangent);
        sine = cosine * tangent;
      }

      // Compute the modified diagonal element of r and the
      //   modified element of [qtb,0]

      {
        double new_rkk = cosine * rkk + sine * sdiagk;
        double new_wak = cosine * wak + sine * qtbpj;

        qtbpj = -sine * wak + cosine * qtbpj;

        fit.apl_matrix_set<utype>(r, k, k, new_rkk);
        fit.apl_vector_set<utype>(wa, k, new_wak);
      }

      // Accumulate the transformation in the row of s //

      for (i = k + 1; i < n; i++) {
        double rik = fit.apl_matrix_get<utype>(r, i, k);
        double sdiagi = fit.apl_vector_get<utype>(sdiag, i);

        double new_rik = cosine * rik + sine * sdiagi;
        double new_sdiagi = -sine * rik + cosine * sdiagi;

        fit.apl_matrix_set<utype>(r, i, k, new_rik);
        fit.apl_vector_set<utype>(sdiag, i, new_sdiagi);
      }
    }

    // Store the corresponding diagonal element of s and restore the
    //   corresponding diagonal element of r //

    {
      double rjj = fit.apl_matrix_get<utype>(r, j, j);
      double xj = fit.apl_vector_get<utype>(x, j);

      fit.apl_vector_set<utype>(sdiag, j, rjj);
      fit.apl_matrix_set<utype>(r, j, j, xj);
    }
  }

  // Solve the triangular system for z. If the system is singular then
  //  obtain a least squares solution

  nsing = n;

  for (j = 0; j < n; j++) {
    double sdiagj = fit.apl_vector_get<utype>(sdiag, j);

    if (_iszero(sdiagj)) {
      nsing = j;
      break;
    }
  }

  for (j = nsing; j < n; j++) {
    fit.apl_vector_set<utype>(wa, j, 0.0);
  }

  for (k = 0; k < nsing; k++) {
    double sum = 0;

    j = (nsing - 1) - k;

    for (i = j + 1; i < nsing; i++) {
      sum += fit.apl_matrix_get<utype>(r, i, j) * fit.apl_vector_get<utype>(wa, i);
    }

    {
      double waj = fit.apl_vector_get<utype>(wa, j);
      double sdiagj = fit.apl_vector_get<utype>(sdiag, j);

      fit.apl_vector_set<utype>(wa, j, (waj - sum) / sdiagj);
    }
  }

  // Permute the components of z back to the components of x //

  for (j = 0; j < n; j++) {
    size_t pj = fit.apl_permutation_get(p, j);
    double waj = fit.apl_vector_get<utype>(wa, j);

    fit.apl_vector_set<utype>(x, pj, waj);
  }

  return 0;
}
// ******************************************************
template <class utype>
void MVops::compute_newton_correction(const apl_matrix<utype> *r, const apl_vector<utype> *sdiag,
                                      const apl_permutation *p, apl_vector<utype> *x,
                                      double dxnorm,
                                      const apl_vector<utype> *diag, apl_vector<utype> *w) {
  MVops fit;
  size_t n = r->size2;
  size_t i, j;

  for (i = 0; i < n; i++) {
    size_t perm_i = fit.apl_permutation_get(p, i);

    double dpi = fit.apl_vector_get<utype>(diag, perm_i);
    double xpi = fit.apl_vector_get<utype>(x, perm_i);

    fit.apl_vector_set<utype>(w, i, dpi * (dpi * xpi) / dxnorm);
  }

  for (j = 0; j < n; j++) {
    double sj = fit.apl_vector_get<utype>(sdiag, j);
    double wj = fit.apl_vector_get<utype>(w, j);

    double tj = wj / sj;

    fit.apl_vector_set<utype>(w, j, tj);

    for (i = j + 1; i < n; i++) {
      double rij = fit.apl_matrix_get<utype>(r, i, j);
      double wi = fit.apl_vector_get<utype>(w, i);

      fit.apl_vector_set<utype>(w, i, wi - rij * tj);
    }
  }
}
// ******************************************************
template <class utype>
int MVops::lmpar(apl_matrix<utype> *r, const apl_permutation *perm, const apl_vector<utype> *qtf,
                 const apl_vector<utype> *diag, double delta, double *par_inout,
                 apl_vector<utype> *newton, apl_vector<utype> *gradient, apl_vector<utype> *sdiag,
                 apl_vector<utype> *x, apl_vector<utype> *w) {
  MVops fit;
  double dxnorm, gnorm, fp, fp_old, par_lower, par_upper, par_c;

  double par = *par_inout;

  size_t iter = 0;

  fit.compute_newton_direction(r, perm, qtf, newton);

  // Evaluate the function at the origin and test for acceptance of
  //   the Gauss-Newton direction.
  dxnorm = fit.scaled_enorm(diag, newton);

  fp = dxnorm - delta;

  if (fp <= 0.1 * delta) {
    fit.apl_vector_memcpy(x, newton);
    *par_inout = 0;
    return 0;
  }

  fit.compute_newton_bound(r, newton, dxnorm, perm, diag, w);

  {
    double wnorm = enorm(w);
    double phider = wnorm * wnorm;

    // w == zero if r rank-deficient,
    //   then set lower bound to zero form MINPACK, lmder.f
    //   Hans E. Plesser 2002-02-25 (hans.plesser@itf.nlh.no)

    if (wnorm > 0)
      par_lower = fp / (delta * phider);
    else
      par_lower = 0.0;
  }

  fit.compute_gradient_direction(r, perm, qtf, diag, gradient);

  gnorm = fit.enorm<utype>(gradient);

  par_upper = gnorm / delta;

  if (par_upper == 0) {
    par_upper = APL_DBL_MIN / std::min(delta, 0.1);
  }

  if (par > par_upper) {
    par = par_upper;
  } else if (par < par_lower) {
    par = par_lower;
  }

  if (par == 0) {
    par = gnorm / dxnorm;
  }

// Beginning of iteration //

iteration:

  iter++;

  // Evaluate the function at the current value of par //

  if (par == 0) {
    par = std::max(0.001 * par_upper, APL_DBL_MIN);
  }

  // Compute the least squares solution of [ R P x - Q^T f, sqrt(par) D x]
  //   for A = Q R P^T //

  {
    double sqrt_par = sqrt(par);

    fit.qrsolv(r, perm, sqrt_par, diag, qtf, x, sdiag, w);
  }

  dxnorm = fit.scaled_enorm(diag, x);

  fp_old = fp;

  fp = dxnorm - delta;

  // If the function is small enough, accept the current value of par //

  if (fabs(fp) <= 0.1 * delta)
    goto line220;

  if (par_lower == 0 && fp <= fp_old && fp_old < 0)
    goto line220;

  // Check for maximum number of iterations //

  if (iter == 10)
    goto line220;

  // Compute the Newton correction //

  fit.compute_newton_correction<utype>(r, sdiag, perm, x, dxnorm, diag, w);

  {
    double wnorm = fit.enorm(w);
    par_c = fp / (delta * wnorm * wnorm);
  }

  // Depending on the sign of the function, update par_lower or par_upper //

  if (fp > 0) {
    if (par > par_lower) {
      par_lower = par;
    }
  } else if (fp < 0) {
    if (par < par_upper) {
      par_upper = par;
    }
  }

  // Compute an improved estimate for par //

  par = std::max(par_lower, par + par_c);

  goto iteration;

line220:

  *par_inout = par;

  return 0;
}
// ******************************************************
template <class utype>
void MVops::compute_trial_step(apl_vector<utype> *x, apl_vector<utype> *dx, apl_vector<utype> *x_trial) {
  size_t i, N = x->size;
  MVops fit;
  for (i = 0; i < N; i++) {
    double perm_i = fit.apl_vector_get<utype>(dx, i);
    double xi = fit.apl_vector_get<utype>(x, i);
    fit.apl_vector_set<utype>(x_trial, i, xi + perm_i);
  }
}
// ******************************************************
double
MVops::compute_actual_reduction(double fnorm, double fnorm1) {
  double actred;

  if (0.1 * fnorm1 < fnorm) {
    double u = fnorm1 / fnorm;
    actred = 1 - u * u;
  } else {
    actred = -1;
  }

  return actred;
}
// ******************************************************
template <class utype>
void MVops::compute_rptdx(const apl_matrix<utype> *r, const apl_permutation *p,
                          const apl_vector<utype> *dx, apl_vector<utype> *rptdx) {
  MVops fit;
  size_t i, j, N = dx->size;

  for (i = 0; i < N; i++) {
    double sum = 0;

    for (j = i; j < N; j++) {
      size_t pj = fit.apl_permutation_get(p, j);

      sum += fit.apl_matrix_get<utype>(r, i, j) * fit.apl_vector_get<utype>(dx, pj);
    }

    fit.apl_vector_set<utype>(rptdx, i, sum);
  }
}
// ******************************************************
template <class utype>
void MVops::update_diag(const apl_matrix<utype> *J, apl_vector<utype> *diag) {
  MVops fit;
  size_t i, j, n = diag->size;

  for (j = 0; j < n; j++) {
    double cnorm, diagj, sum = 0;
    for (i = 0; i < n; i++) {
      double Jij = fit.apl_matrix_get<utype>(J, i, j);
      sum += Jij * Jij;
    }
    if (_iszero(sum))
      sum = 1.0;

    cnorm = sqrt(sum);
    diagj = fit.apl_vector_get<utype>(diag, j);

    if (cnorm > diagj)
      fit.apl_vector_set<utype>(diag, j, cnorm);
  }
}
// ******************************************************
template <class utype>
int MVops::apl_multifit_fdfsolver_dif_df(const apl_vector<utype> *x, apl_multifit_function_fdf *fdf,
                                         const apl_vector<utype> *f, apl_matrix<utype> *J) {
  return fdjac<utype>(x, fdf, f, J);
}
// ******************************************************
template <class utype>
int MVops::iterate(void *vstate, apl_multifit_function_fdf *fdf, apl_vector<utype> *x, apl_vector<utype> *f, apl_matrix<utype> *J, apl_vector<utype> *dx, int scale) {
  MVops fit;

  lmder_state_t<utype> *state = (lmder_state_t<utype> *)vstate;

  apl_matrix<utype> *r = state->r;
  apl_vector<utype> *tau = state->tau;
  apl_vector<utype> *diag = state->diag;
  apl_vector<utype> *qtf = state->qtf;
  apl_vector<utype> *x_trial = state->x_trial;
  apl_vector<utype> *f_trial = state->f_trial;
  apl_vector<utype> *rptdx = state->rptdx;
  apl_vector<utype> *newton = state->newton;
  apl_vector<utype> *gradient = state->gradient;
  apl_vector<utype> *sdiag = state->sdiag;
  apl_vector<utype> *w = state->w;
  apl_vector<utype> *work1 = state->work1;
  apl_permutation *perm = state->perm;

  double prered, actred;
  double pnorm, fnorm1, fnorm1p, gnorm;
  double ratio;
  double dirder;

  int iter = 0;

  double p1 = 0.1, p25 = 0.25, p5 = 0.5, p75 = 0.75, p0001 = 0.0001;

  if (_iszero(state->fnorm)) {
    return 0;
  }

  // Compute qtf = Q^T f //

  fit.apl_vector_memcpy<utype>(qtf, f);
  fit.apl_linalg_QR_QTvec<utype>(r, tau, qtf);

  // Compute norm of scaled gradient //

  fit.compute_gradient_direction<utype>(r, perm, qtf, diag, gradient);

  {
    size_t iamax = fit.apl_blas_idamax<utype>(gradient);

    gnorm = fabs(fit.apl_vector_get<utype>(gradient, iamax) / state->fnorm);
  }

// Determine the Levenberg-Marquardt parameter //

lm_iteration:

  iter++;

  {
    int status = fit.lmpar<utype>(r, perm, qtf, diag, state->delta, &(state->par), newton, gradient, sdiag, dx, w);
    if (status)
      return status;
  }

  // Take a trial step //

  fit.apl_vector_scale<utype>(dx, -1.0);  // reverse the step to go downhill //

  fit.compute_trial_step<utype>(x, dx, state->x_trial);

  pnorm = fit.scaled_enorm(diag, dx);

  if (state->iter == 1) {
    if (pnorm < state->delta) {
      state->delta = pnorm;
    }
  }

  // Evaluate function at x + p //
  // return immediately if evaluation raised error //
  {
    int status = APL_MULTIFIT_FN_EVAL_F(fdf, x_trial, f_trial);
    if (status)
      return status;
  }

  fnorm1 = fit.enorm(f_trial);

  // Compute the scaled actual reduction

  actred = fit.compute_actual_reduction(state->fnorm, fnorm1);

  // Compute rptdx = R P^T dx, noting that |J dx| = |R P^T dx| //

  fit.compute_rptdx(r, perm, dx, rptdx);

  fnorm1p = fit.enorm(rptdx);

  // Compute the scaled predicted reduction = |J dx|^2 + 2 par |D dx|^2 //

  {
    double t1 = fnorm1p / state->fnorm;
    double t2 = (sqrt(state->par) * pnorm) / state->fnorm;

    prered = t1 * t1 + t2 * t2 / p5;
    dirder = -(t1 * t1 + t2 * t2);
  }

  // compute the ratio of the actual to predicted reduction //

  if (prered > 0) {
    ratio = actred / prered;
  } else {
    ratio = 0;
  }

  // update the step bound //

  if (ratio > p25) {
    if (state->par == 0 || ratio >= p75) {
      state->delta = pnorm / p5;
      state->par *= p5;
    }
  } else {
    if (_iszero(actred)) actred = 0.00;
    double temp = (actred >= 0) ? p5 : p5 * dirder / (dirder + p5 * actred);

    if (p1 * fnorm1 >= state->fnorm || temp < p1) {
      temp = p1;
    }

    state->delta = temp * std::min(state->delta, pnorm / p1);

    state->par /= temp;
  }

  // test for successful iteration, termination and stringent tolerances //
  if (ratio >= p0001) {
    fit.apl_vector_memcpy(x, x_trial);
    fit.apl_vector_memcpy(f, f_trial);

    // return immediately if evaluation raised error //
    {
      int status;

      if (fdf->df)
        status = APL_MULTIFIT_FN_EVAL_DF(fdf, x_trial, J);
      else
        status = fit.apl_multifit_fdfsolver_dif_df(x_trial, fdf, f_trial, J);

      if (status)
        return status;
    }
    // wa2_j  = diag_j * x_j //
    state->xnorm = fit.scaled_enorm(diag, x);
    state->fnorm = fnorm1;
    state->iter++;

    // Rescale if necessary //

    if (scale) {
      fit.update_diag(J, diag);
    }

    {
      int signum;
      fit.apl_matrix_memcpy<utype>(r, J);
      fit.apl_linalg_QRPT_decomp<utype>(r, tau, perm, &signum, work1);
    }
    return 0;
  } else if (fabs(actred) <= APL_DBL_EPSILON && prered <= APL_DBL_EPSILON && p5 * ratio <= 1.0) {
    return 29;  // cannot reach the specified tolerance in F //
  } else if (state->delta <= APL_DBL_EPSILON * state->xnorm) {
    return 30;  // cannot reach the specified tolerance in X //
  } else if (gnorm <= APL_DBL_EPSILON) {
    return 31;  // cannot reach the specified tolerance in gradient //
  } else if (iter < 1000) {
    // Repeat inner loop if unsuccessful //
    goto lm_iteration;
  }

  return 37;  //iteration is not making progress towards solution
}
// ******************************************************
template <class utype>
void MVops::apl_multifit_fdfsolver_free(apl_multifit_fdfsolver<utype> *s) {
  RETURN_IF_NULL(s);
  (s->type->free)(s->state);
  free(s->state);
  apl_vector_free<utype>(s->dx);
  apl_vector_free<utype>(s->x);
  apl_vector_free<utype>(s->f);
  apl_matrix_free<utype>(s->J);
  free(s);
}
// ******************************************************
template <class utype>
int MVops::lmder_alloc(void *vstate, size_t n, size_t p) {
  lmder_state_t<utype> *state = (lmder_state_t<utype> *)vstate;
  apl_matrix<utype> *r;
  apl_vector<utype> *tau, *diag, *qtf, *newton, *gradient, *x_trial, *f_trial,
      *df, *sdiag, *rptdx, *w, *work1;
  apl_permutation *perm;

  MVops fit;

  r = fit.apl_matrix_calloc<utype>(n, p);

  if (r == 0) {
    cerr << "failed to allocate space for r" << std::endl;
    exit(0);
  }

  state->r = r;

  tau = fit.apl_vector_calloc<utype>(std::min(n, p));

  if (tau == 0) {
    fit.apl_matrix_free<utype>(r);

    cerr << "failed to allocate space for tau" << std::endl;
    exit(0);
  }

  state->tau = tau;

  diag = fit.apl_vector_calloc<utype>(p);

  if (diag == 0) {
    fit.apl_matrix_free<utype>(r);
    fit.apl_vector_free<utype>(tau);

    cerr << "failed to allocate space for diag" << std::endl;
    exit(0);
  }

  state->diag = diag;

  qtf = fit.apl_vector_calloc<utype>(n);

  if (qtf == 0) {
    fit.apl_matrix_free<utype>(r);
    fit.apl_vector_free<utype>(tau);
    fit.apl_vector_free<utype>(diag);

    cerr << "failed to allocate space for qtf" << std::endl;
    exit(0);
  }

  state->qtf = qtf;

  newton = fit.apl_vector_calloc<utype>(p);

  if (newton == 0) {
    fit.apl_matrix_free<utype>(r);
    fit.apl_vector_free<utype>(tau);
    fit.apl_vector_free<utype>(diag);
    fit.apl_vector_free<utype>(qtf);

    cerr << "failed to allocate space for newton" << std::endl;
    exit(0);
  }

  state->newton = newton;

  gradient = fit.apl_vector_calloc<utype>(p);

  if (gradient == 0) {
    fit.apl_matrix_free<utype>(r);
    fit.apl_vector_free<utype>(tau);
    fit.apl_vector_free<utype>(diag);
    fit.apl_vector_free<utype>(qtf);
    fit.apl_vector_free<utype>(newton);

    cerr << "failed to allocate space for gradient" << std::endl;
    exit(0);
  }

  state->gradient = gradient;

  x_trial = fit.apl_vector_calloc<utype>(p);

  if (x_trial == 0) {
    fit.apl_matrix_free<utype>(r);
    fit.apl_vector_free<utype>(tau);
    fit.apl_vector_free<utype>(diag);
    fit.apl_vector_free<utype>(qtf);
    fit.apl_vector_free<utype>(newton);
    fit.apl_vector_free<utype>(gradient);

    cerr << "failed to allocate space for x_trial" << std::endl;
    exit(0);
  }

  state->x_trial = x_trial;

  f_trial = fit.apl_vector_calloc<utype>(n);

  if (f_trial == 0) {
    fit.apl_matrix_free<utype>(r);
    fit.apl_vector_free<utype>(tau);
    fit.apl_vector_free<utype>(diag);
    fit.apl_vector_free<utype>(qtf);
    fit.apl_vector_free<utype>(newton);
    fit.apl_vector_free<utype>(gradient);
    fit.apl_vector_free<utype>(x_trial);

    cerr << "failed to allocate space for f_trial" << std::endl;
    exit(0);
  }

  state->f_trial = f_trial;

  df = fit.apl_vector_calloc<utype>(n);

  if (df == 0) {
    fit.apl_matrix_free<utype>(r);
    fit.apl_vector_free<utype>(tau);
    fit.apl_vector_free<utype>(diag);
    fit.apl_vector_free<utype>(qtf);
    fit.apl_vector_free<utype>(newton);
    fit.apl_vector_free<utype>(gradient);
    fit.apl_vector_free<utype>(x_trial);
    fit.apl_vector_free<utype>(f_trial);

    cerr << "failed to allocate space for df" << std::endl;
    exit(0);
  }

  state->df = df;

  sdiag = fit.apl_vector_calloc<utype>(p);

  if (sdiag == 0) {
    fit.apl_matrix_free<utype>(r);
    fit.apl_vector_free<utype>(tau);
    fit.apl_vector_free<utype>(diag);
    fit.apl_vector_free<utype>(qtf);
    fit.apl_vector_free<utype>(newton);
    fit.apl_vector_free<utype>(gradient);
    fit.apl_vector_free<utype>(x_trial);
    fit.apl_vector_free<utype>(f_trial);
    fit.apl_vector_free<utype>(df);

    cerr << "failed to allocate space for sdiag" << std::endl;
    exit(0);
  }

  state->sdiag = sdiag;

  rptdx = fit.apl_vector_calloc<utype>(n);

  if (rptdx == 0) {
    fit.apl_matrix_free<utype>(r);
    fit.apl_vector_free<utype>(tau);
    fit.apl_vector_free<utype>(diag);
    fit.apl_vector_free<utype>(qtf);
    fit.apl_vector_free<utype>(newton);
    fit.apl_vector_free<utype>(gradient);
    fit.apl_vector_free<utype>(x_trial);
    fit.apl_vector_free<utype>(f_trial);
    fit.apl_vector_free<utype>(df);
    fit.apl_vector_free<utype>(sdiag);

    cerr << "failed to allocate space for rptdx" << std::endl;
    exit(0);
  }

  state->rptdx = rptdx;

  w = fit.apl_vector_calloc<utype>(n);

  if (w == 0) {
    fit.apl_matrix_free<utype>(r);
    fit.apl_vector_free<utype>(tau);
    fit.apl_vector_free<utype>(diag);
    fit.apl_vector_free<utype>(qtf);
    fit.apl_vector_free<utype>(newton);
    fit.apl_vector_free<utype>(gradient);
    fit.apl_vector_free<utype>(x_trial);
    fit.apl_vector_free<utype>(f_trial);
    fit.apl_vector_free<utype>(df);
    fit.apl_vector_free<utype>(sdiag);
    fit.apl_vector_free<utype>(rptdx);

    cerr << "failed to allocate space for w" << std::endl;
    exit(0);
  }

  state->w = w;

  work1 = fit.apl_vector_calloc<utype>(p);

  if (work1 == 0) {
    fit.apl_matrix_free<utype>(r);
    fit.apl_vector_free<utype>(tau);
    fit.apl_vector_free<utype>(diag);
    fit.apl_vector_free<utype>(qtf);
    fit.apl_vector_free<utype>(newton);
    fit.apl_vector_free<utype>(gradient);
    fit.apl_vector_free<utype>(x_trial);
    fit.apl_vector_free<utype>(f_trial);
    fit.apl_vector_free<utype>(df);
    fit.apl_vector_free<utype>(sdiag);
    fit.apl_vector_free<utype>(rptdx);
    fit.apl_vector_free<utype>(w);

    cerr << "failed to allocate space for work1" << std::endl;
    exit(0);
  }

  state->work1 = work1;

  perm = fit.apl_permutation_calloc(p);

  if (perm == 0) {
    fit.apl_matrix_free<utype>(r);
    fit.apl_vector_free<utype>(tau);
    fit.apl_vector_free<utype>(diag);
    fit.apl_vector_free<utype>(qtf);
    fit.apl_vector_free<utype>(newton);
    fit.apl_vector_free<utype>(gradient);
    fit.apl_vector_free<utype>(x_trial);
    fit.apl_vector_free<utype>(f_trial);
    fit.apl_vector_free<utype>(df);
    fit.apl_vector_free<utype>(sdiag);
    fit.apl_vector_free<utype>(rptdx);
    fit.apl_vector_free<utype>(w);
    fit.apl_vector_free<utype>(work1);

    cerr << "failed to allocate space for perm" << std::endl;
    exit(0);
  }

  state->perm = perm;

  return 0;
}
// ******************************************************
template <class utype>
int MVops::fdjac(const apl_vector<utype> *x, apl_multifit_function_fdf *fdf,
                 const apl_vector<utype> *f, apl_matrix<utype> *J) {
  MVops fit;
  int status = 0;
  size_t i, j;
  double h;
  const double epsfcn = 0.0;
  double eps = sqrt(std::max(epsfcn, APL_DBL_EPSILON));

  for (j = 0; j < fdf->p; ++j) {
    double xj = fit.apl_vector_get<utype>(x, j);

    // use column j of J as temporary storage for f(x + dx) //
    apl_vector_view<utype> v = fit.apl_matrix_column<utype>(J, j);

    h = eps * fabs(xj);
    if (_iszero(h))
      h = eps;

    // perturb x_j to compute forward difference //
    fit.apl_vector_set<utype>((apl_vector<utype> *)x, j, xj + h);

    status += APL_MULTIFIT_FN_EVAL_F(fdf, x, &v.vector);
    if (status)
      return status;

    // restore x_j //
    fit.apl_vector_set<utype>((apl_vector<utype> *)x, j, xj);

    h = 1.0 / h;
    for (i = 0; i < fdf->n; ++i) {
      double fnext = fit.apl_vector_get<utype>(&v.vector, i);
      double fi = fit.apl_vector_get<utype>(f, i);

      fit.apl_matrix_set<utype>(J, i, j, (fnext - fi) * h);
    }
  }

  return status;
}
// ******************************************************
template <class utype>
int MVops::apl_multifit_fdfsolver_dif_fdf(const apl_vector<utype> *x, apl_multifit_function_fdf *fdf,
                                          apl_vector<utype> *f, apl_matrix<utype> *J) {
  int status = 0;

  status = APL_MULTIFIT_FN_EVAL_F(fdf, x, f);
  if (status)
    return status;

  status = fdjac<utype>(x, fdf, f, J);
  if (status)
    return status;

  return status;
}
// ******************************************************
template <class utype>
void MVops::compute_diag(const apl_matrix<utype> *J, apl_vector<utype> *diag) {
  MVops fit;
  size_t i, j, n = J->size1, p = J->size2;

  for (j = 0; j < p; j++) {
    double sum = 0;

    for (i = 0; i < n; i++) {
      double Jij = fit.apl_matrix_get<utype>(J, i, j);
      sum += Jij * Jij;
    }
    if (_iszero(sum))
      sum = 1.0;

    fit.apl_vector_set<utype>(diag, j, sqrt(sum));
  }
}
// ******************************************************
template <class utype>
double
MVops::compute_delta(apl_vector<utype> *diag, apl_vector<utype> *x) {
  MVops fit;
  double Dx = fit.scaled_enorm(diag, x);
  double factor = 100;  // generally recommended value from MINPACK //

  return (Dx > 0) ? factor * Dx : factor;
}
// ******************************************************
template <class utype>
int MVops::set(void *vstate, apl_multifit_function_fdf *fdf, apl_vector<utype> *x, apl_vector<utype> *f, apl_matrix<utype> *J, apl_vector<utype> *dx, int scale) {
  lmder_state_t<utype> *state = (lmder_state_t<utype> *)vstate;

  apl_matrix<utype> *r = state->r;
  apl_vector<utype> *tau = state->tau;
  apl_vector<utype> *diag = state->diag;
  apl_vector<utype> *work1 = state->work1;
  apl_permutation *perm = state->perm;

  int signum;

  MVops fit;
  // Evaluate function at x //
  // return immediately if evaluation raised error //
  {
    int status;
    if (fdf->fdf)
      status = APL_MULTIFIT_FN_EVAL_F_DF(fdf, x, f, J);
    else  // finite difference approximation //
      status = fit.apl_multifit_fdfsolver_dif_fdf(x, fdf, f, J);
    if (status)
      return status;
  }
  state->par = 0;
  state->iter = 1;
  state->fnorm = fit.enorm(f);

  fit.apl_vector_set_all<utype>(dx, 0.0);

  // store column norms in diag //

  if (scale) {
    fit.compute_diag<utype>(J, diag);
  } else {
    fit.apl_vector_set_all<utype>(diag, 1.0);
  }

  // set delta to 100 |D x| or to 100 if |D x| is zero //

  state->xnorm = fit.scaled_enorm(diag, x);
  state->delta = fit.compute_delta<utype>(diag, x);

  // Factorize J into QR decomposition //

  fit.apl_matrix_memcpy<utype>(r, J);
  fit.apl_linalg_QRPT_decomp<utype>(r, tau, perm, &signum, work1);

  fit.apl_vector_set_zero<utype>(state->rptdx);
  fit.apl_vector_set_zero<utype>(state->w);

  // Zero the trial vector, as in the alloc function //

  fit.apl_vector_set_zero<utype>(state->f_trial);

  return 0;
}
// ******************************************************
template <class utype>
int MVops::lmsder_set(void *vstate, apl_multifit_function_fdf *fdf, apl_vector<utype> *x, apl_vector<utype> *f, apl_matrix<utype> *J, apl_vector<utype> *dx) {
  MVops fit;
  int status = fit.set<double>(vstate, fdf, x, f, J, dx, 1);
  return status;
}
// ******************************************************
template <class utype>
int MVops::lmder_set(void *vstate, apl_multifit_function_fdf *fdf, apl_vector<utype> *x, apl_vector<utype> *f, apl_matrix<utype> *J, apl_vector<utype> *dx) {
  MVops fit;
  int status = fit.set<double>(vstate, fdf, x, f, J, dx, 0);
  return status;
}
// ******************************************************
template <class utype>
int MVops::lmsder_iterate(void *vstate, apl_multifit_function_fdf *fdf, apl_vector<utype> *x, apl_vector<utype> *f, apl_matrix<utype> *J, apl_vector<utype> *dx) {
  MVops fit;
  int status = fit.iterate(vstate, fdf, x, f, J, dx, 1);
  return status;
}
// ******************************************************
template <class utype>
int MVops::lmder_iterate(void *vstate, apl_multifit_function_fdf *fdf, apl_vector<utype> *x, apl_vector<utype> *f, apl_matrix<utype> *J, apl_vector<utype> *dx) {
  MVops fit;
  int status = fit.iterate(vstate, fdf, x, f, J, dx, 0);
  return status;
}
// ******************************************************
template <class utype>
void MVops::lmder_free(void *vstate) {
  lmder_state_t<utype> *state = (lmder_state_t<utype> *)vstate;

  MVops fit;

  fit.apl_permutation_free(state->perm);
  fit.apl_vector_free<utype>(state->work1);
  fit.apl_vector_free<utype>(state->w);
  fit.apl_vector_free<utype>(state->rptdx);
  fit.apl_vector_free<utype>(state->sdiag);
  fit.apl_vector_free<utype>(state->df);
  fit.apl_vector_free<utype>(state->f_trial);
  fit.apl_vector_free<utype>(state->x_trial);
  fit.apl_vector_free<utype>(state->gradient);
  fit.apl_vector_free<utype>(state->newton);
  fit.apl_vector_free<utype>(state->qtf);
  fit.apl_vector_free<utype>(state->diag);
  fit.apl_vector_free<utype>(state->tau);
  fit.apl_matrix_free<utype>(state->r);
}
// ******************************************************
apl_eigen_hermv_workspace *
MVops::apl_eigen_hermv_alloc(const size_t n) {
  apl_eigen_hermv_workspace *w;

  if (n == 0) {
    cerr << "matrix dimension must be positive integer";
    exit(0);
  }

  w = (apl_eigen_hermv_workspace *)malloc(sizeof(apl_eigen_hermv_workspace));

  if (w == 0) {
    cerr << "failed to allocate space for workspace";
    exit(0);
  }

  w->d = (double *)malloc(n * sizeof(double));

  if (w->d == 0) {
    free(w);
    cerr << "failed to allocate space for diagonal";
    exit(0);
  }

  w->sd = (double *)malloc(n * sizeof(double));

  if (w->sd == 0) {
    free(w->d);
    free(w);
    cerr << "failed to allocate space for subdiagonal";
    exit(0);
  }

  w->tau = (double *)malloc(2 * n * sizeof(double));

  if (w->tau == 0) {
    free(w->sd);
    free(w->d);
    free(w);
    cerr << "failed to allocate space for tau";
    exit(0);
  }

  w->gc = (double *)malloc(n * sizeof(double));

  if (w->gc == 0) {
    free(w->tau);
    free(w->sd);
    free(w->d);
    free(w);
    cerr << "failed to allocate space for cosines";
    exit(0);
  }

  w->gs = (double *)malloc(n * sizeof(double));

  if (w->gs == 0) {
    free(w->gc);
    free(w->tau);
    free(w->sd);
    free(w->d);
    free(w);
    cerr << "failed to allocate space for sines";
    exit(0);
  }

  w->size = n;

  return w;
}
// ******************************************************
void MVops::apl_eigen_hermv_free(apl_eigen_hermv_workspace *w) {
  RETURN_IF_NULL(w);
  free(w->gs);
  free(w->gc);
  free(w->tau);
  free(w->sd);
  free(w->d);
  free(w);
}
// ******************************************************
template <class T>
double MVops::apl_blas_dznrm2(const apl_vector_complex<T> *X) {
  return cblas_dznrm2(int(X->size), X->data, int(X->stride));
}
// ******************************************************
double
MVops::cblas_dznrm2(const int N, const void *X, const int incX) {
  double scale = 0.0;
  double ssq = 1.0;
  int i;
  int ix = 0;

  if (N == 0 || incX < 1) {
    return 0;
  }

  for (i = 0; i < N; i++) {
    const double x = CONST_REAL(X, ix);
    const double y = CONST_IMAG(X, ix);

    if (x != 0.0) {
      const double ax = std::fabs(x);

      if (scale < ax) {
        ssq = 1.0 + ssq * (scale / ax) * (scale / ax);
        scale = ax;
      } else {
        ssq += (ax / scale) * (ax / scale);
      }
    }

    if (y != 0.0) {
      const double ay = std::fabs(y);

      if (scale < ay) {
        ssq = 1.0 + ssq * (scale / ay) * (scale / ay);
        scale = ay;
      } else {
        ssq += (ay / scale) * (ay / scale);
      }
    }

    ix += incX;
  }

  return scale * sqrt(ssq);
}
// ******************************************************
template <class T>
void MVops::apl_blas_zscal(const apl_complex<T> alpha, apl_vector_complex<T> *X) {
  cblas_zscal(int(X->size), APL_MV_COMPLEX_P(&alpha), X->data,
              int(X->stride));
}
// ******************************************************
void MVops::cblas_zscal(const int N, const void *alpha, void *X, const int incX) {
  int i;
  int ix = 0;
  const double alpha_real = CONST_REAL0(alpha);
  const double alpha_imag = CONST_IMAG0(alpha);

  if (incX <= 0) {
    return;
  }

  for (i = 0; i < N; i++) {
    const double x_real = REAL(X, ix);
    const double x_imag = IMAG(X, ix);
    REAL(X, ix) = x_real * alpha_real - x_imag * alpha_imag;
    IMAG(X, ix) = x_real * alpha_imag + x_imag * alpha_real;
    ix += incX;
  }
}
// ******************************************************
template <class T>
apl_complex<T>
MVops::apl_linalg_complex_householder_transform(apl_vector_complex<T> *v) {
  // replace v[0:n-1] with a householder vector (v[0:n-1]) and
  //   coefficient tau that annihilate v[1:n-1] //

  const size_t n = v->size;

  if (n == 1) {
    apl_complex<T> alpha = apl_vector_complex_get<T>(v, 0);
    double absa = apl_complex_abs<T>(alpha);
    if (_iszero(APL_MV_REAL(alpha))) APL_MV_REAL(alpha) = 0.00;
    double beta_r = -(APL_MV_REAL(alpha) >= 0 ? +1 : -1) * absa;

    apl_complex<T> tau;

    if (_iszero(beta_r)) {
      APL_MV_REAL(tau) = 0.0;
      APL_MV_IMAG(tau) = 0.0;
    } else {
      APL_MV_REAL(tau) = (beta_r - APL_MV_REAL(alpha)) / beta_r;
      APL_MV_IMAG(tau) = -APL_MV_IMAG(alpha) / beta_r;

      {
        apl_complex<T> beta = apl_complex_rect(beta_r, 0.0);
        apl_vector_complex_set<T>(v, 0, beta);
      }
    }

    return tau;
  } else {
    apl_complex<T> tau;
    double beta_r;

    apl_vector_complex_view<T> x = apl_vector_complex_subvector<T>(v, 1, n - 1);
    apl_complex<T> alpha = apl_vector_complex_get<T>(v, 0);
    double absa = apl_complex_abs<T>(alpha);
    double xnorm = apl_blas_dznrm2<T>(&x.vector);
    if ((_iszero(xnorm)) && (_iszero(APL_MV_IMAG(alpha)))) {
      apl_complex<T> zero = apl_complex_rect(0.0, 0.0);
      return zero;  // tau = 0 //
    }

    beta_r = -(APL_MV_REAL(alpha) >= 0 ? +1 : -1) * hypot(absa, xnorm);

    APL_MV_REAL(tau) = (beta_r - APL_MV_REAL(alpha)) / beta_r;
    APL_MV_IMAG(tau) = -APL_MV_IMAG(alpha) / beta_r;

    {
      apl_complex<T> amb = apl_complex_sub_real<T>(alpha, beta_r);
      apl_complex<T> s = apl_complex_inverse(amb);
      apl_blas_zscal(s, &x.vector);
    }

    {
      apl_complex<T> beta = apl_complex_rect(beta_r, 0.0);
      apl_vector_complex_set<T>(v, 0, beta);
    }
    return tau;
  }
}
// ******************************************************
template <class T>
int MVops::apl_blas_zhemv(APL_UPLO Uplo, const apl_complex<T> alpha,
                          const apl_matrix_complex<T> *A, const apl_vector_complex<T> *X,
                          const apl_complex<T> beta, apl_vector_complex<T> *Y) {
  const size_t M = A->size1;
  const size_t N = A->size2;

  if (M != N) {
    cerr << "\nmatrix must be square\n";
    exit(0);
  } else if (N != X->size || N != Y->size) {
    cerr << "\ninvalid length\n";
    exit(0);
  }

  cblas_zhemv(aplRowMajor, Uplo, int(N), APL_MV_COMPLEX_P(&alpha), A->data,
              int(A->tda), X->data, int(X->stride), APL_MV_COMPLEX_P(&beta),
              Y->data, int(Y->stride));
  return 0;
}
// ******************************************************
void MVops::cblas_zhemv(const enum APL_ORDER order, const enum APL_UPLO Uplo,
                        const int N, const void *alpha, const void *A, const int lda,
                        const void *X, const int incX, const void *beta, void *Y,
                        const int incY) {
  const int conj = (order == aplColMajor) ? -1 : 1;
  int i, j;

  //CHECK_ARGS11(CZ_HEMV,order,Uplo,N,alpha,A,lda,X,incX,beta,Y,incY);

  {
    const double alpha_real = CONST_REAL0(alpha);
    const double alpha_imag = CONST_IMAG0(alpha);

    const double beta_real = CONST_REAL0(beta);
    const double beta_imag = CONST_IMAG0(beta);

    if (((_iszero(alpha_real)) && (_iszero(alpha_imag))) && ((_isequal(beta_real, 1.0)) && (_iszero(beta_imag))))
      return;

    // form  y := beta*y //
    if ((_iszero(beta_real)) && (_iszero(beta_imag))) {
      int iy = OFFSET(N, incY);
      for (i = 0; i < N; i++) {
        REAL(Y, iy) = 0.0;
        IMAG(Y, iy) = 0.0;
        iy += incY;
      }
    } else if (!((_isequal(beta_real, 1.0)) && (_iszero(beta_imag)))) {
      int iy = OFFSET(N, incY);
      for (i = 0; i < N; i++) {
        const double y_real = REAL(Y, iy);
        const double y_imag = IMAG(Y, iy);
        const double tmpR = y_real * beta_real - y_imag * beta_imag;
        const double tmpI = y_real * beta_imag + y_imag * beta_real;
        REAL(Y, iy) = tmpR;
        IMAG(Y, iy) = tmpI;
        iy += incY;
      }
    }
    if ((_iszero(alpha_real)) && (_iszero(alpha_imag)))
      return;

    // form  y := alpha*A*x + y //

    if ((order == aplRowMajor && Uplo == aplUpper) || (order == aplColMajor && Uplo == aplLower)) {
      int ix = OFFSET(N, incX);
      int iy = OFFSET(N, incY);
      for (i = 0; i < N; i++) {
        double x_real = CONST_REAL(X, ix);
        double x_imag = CONST_IMAG(X, ix);
        double temp1_real = alpha_real * x_real - alpha_imag * x_imag;
        double temp1_imag = alpha_real * x_imag + alpha_imag * x_real;
        double temp2_real = 0.0;
        double temp2_imag = 0.0;
        const int j_min = i + 1;
        const int j_max = N;
        int jx = OFFSET(N, incX) + j_min * incX;
        int jy = OFFSET(N, incY) + j_min * incY;
        double Aii_real = CONST_REAL(A, lda * i + i);
        // Aii_imag is zero //
        REAL(Y, iy) += temp1_real * Aii_real;
        IMAG(Y, iy) += temp1_imag * Aii_real;
        for (j = j_min; j < j_max; j++) {
          double Aij_real = CONST_REAL(A, lda * i + j);
          double Aij_imag = conj * CONST_IMAG(A, lda * i + j);
          REAL(Y, jy) += temp1_real * Aij_real - temp1_imag * (-Aij_imag);
          IMAG(Y, jy) += temp1_real * (-Aij_imag) + temp1_imag * Aij_real;
          x_real = CONST_REAL(X, jx);
          x_imag = CONST_IMAG(X, jx);
          temp2_real += x_real * Aij_real - x_imag * Aij_imag;
          temp2_imag += x_real * Aij_imag + x_imag * Aij_real;
          jx += incX;
          jy += incY;
        }
        REAL(Y, iy) += alpha_real * temp2_real - alpha_imag * temp2_imag;
        IMAG(Y, iy) += alpha_real * temp2_imag + alpha_imag * temp2_real;
        ix += incX;
        iy += incY;
      }
    } else if ((order == aplRowMajor && Uplo == aplLower) || (order == aplColMajor && Uplo == aplUpper)) {
      int ix = OFFSET(N, incX) + (N - 1) * incX;
      int iy = OFFSET(N, incY) + (N - 1) * incY;
      for (i = N; i > 0 && i--;) {
        double x_real = CONST_REAL(X, ix);
        double x_imag = CONST_IMAG(X, ix);
        double temp1_real = alpha_real * x_real - alpha_imag * x_imag;
        double temp1_imag = alpha_real * x_imag + alpha_imag * x_real;
        double temp2_real = 0.0;
        double temp2_imag = 0.0;
        const int j_min = 0;
        const int j_max = i;
        int jx = OFFSET(N, incX) + j_min * incX;
        int jy = OFFSET(N, incY) + j_min * incY;
        double Aii_real = CONST_REAL(A, lda * i + i);
        // Aii_imag is zero //
        REAL(Y, iy) += temp1_real * Aii_real;
        IMAG(Y, iy) += temp1_imag * Aii_real;

        for (j = j_min; j < j_max; j++) {
          double Aij_real = CONST_REAL(A, lda * i + j);
          double Aij_imag = conj * CONST_IMAG(A, lda * i + j);
          REAL(Y, jy) += temp1_real * Aij_real - temp1_imag * (-Aij_imag);
          IMAG(Y, jy) += temp1_real * (-Aij_imag) + temp1_imag * Aij_real;
          x_real = CONST_REAL(X, jx);
          x_imag = CONST_IMAG(X, jx);
          temp2_real += x_real * Aij_real - x_imag * Aij_imag;
          temp2_imag += x_real * Aij_imag + x_imag * Aij_real;
          jx += incX;
          jy += incY;
        }
        REAL(Y, iy) += alpha_real * temp2_real - alpha_imag * temp2_imag;
        IMAG(Y, iy) += alpha_real * temp2_imag + alpha_imag * temp2_real;
        ix -= incX;
        iy -= incY;
      }
    } else {
      cerr << "unrecognized operation";
      exit(0);
    }
  }
}
// ******************************************************
template <class T>
int MVops::apl_blas_zdotc(const apl_vector_complex<T> *X, const apl_vector_complex<T> *Y,
                          apl_complex<T> *dotc) {
  if (X->size == Y->size) {
    cblas_zdotc_sub(int(X->size), X->data, int(X->stride), Y->data,
                    int(Y->stride), APL_MV_COMPLEX_P(dotc));
    return 0;
  } else {
    cerr << "\n invalid length \n ";
    exit(0);
  }
}
// ******************************************************
void MVops::cblas_zdotc_sub(const int N, const void *X, const int incX, const void *Y,
                            const int incY, void *result) {
#define CONJ_SIGN (-1.0)
  double r_real = 0.0;
  double r_imag = 0.0;
  int i;
  int ix = OFFSET(N, incX);
  int iy = OFFSET(N, incY);
  for (i = 0; i < N; i++) {
    const double x_real = CONST_REAL(X, ix);
    const double x_imag = CONST_IMAG(X, ix);
    const double y_real = CONST_REAL(Y, iy);
    const double y_imag = CONST_IMAG(Y, iy);
    r_real += x_real * y_real - CONJ_SIGN * x_imag * y_imag;
    r_imag += x_real * y_imag + CONJ_SIGN * x_imag * y_real;
    ix += incX;
    iy += incY;
  }
  REAL0(result) = r_real;
  IMAG0(result) = r_imag;
#undef CONJ_SIGN
}
// ******************************************************
template <class T>
int MVops::apl_blas_zaxpy(const apl_complex<T> alpha, const apl_vector_complex<T> *X,
                          apl_vector_complex<T> *Y) {
  if (X->size == Y->size) {
    cblas_zaxpy(int(X->size), APL_MV_COMPLEX_P(&alpha), X->data,
                int(X->stride), Y->data, int(Y->stride));
    return 0;
  } else {
    cerr << "\n invalid length \n ";
    exit(0);
  }
}
// ******************************************************
void MVops::cblas_zaxpy(const int N, const void *alpha, const void *X, const int incX,
                        void *Y, const int incY) {
  int i;
  int ix = OFFSET(N, incX);
  int iy = OFFSET(N, incY);

  const double alpha_real = CONST_REAL0(alpha);
  const double alpha_imag = CONST_IMAG0(alpha);

  if ((_iszero(std::fabs(alpha_real))) && (_iszero(std::fabs(alpha_imag)))) {
    return;
  }

  for (i = 0; i < N; i++) {
    const double x_real = CONST_REAL(X, ix);
    const double x_imag = CONST_IMAG(X, ix);
    REAL(Y, iy) += (alpha_real * x_real - alpha_imag * x_imag);
    IMAG(Y, iy) += (alpha_real * x_imag + alpha_imag * x_real);
    ix += incX;
    iy += incY;
  }
}
// ******************************************************
template <class T>
int MVops::apl_blas_zher2(APL_UPLO_t Uplo, const apl_complex<T> alpha,
                          const apl_vector_complex<T> *X, const apl_vector_complex<T> *Y,
                          apl_matrix_complex<T> *A) {
  const size_t M = A->size1;
  const size_t N = A->size2;

  if (M != N) {
    cerr << "\n matrix must be square \n ";
    exit(0);
  } else if (X->size != N || Y->size != N) {
    cerr << "\n invalid length \n ";
    exit(0);
  }

  cblas_zher2(aplRowMajor, Uplo, int(N), APL_MV_COMPLEX_P(&alpha), X->data,
              int(X->stride), Y->data, int(Y->stride), A->data,
              int(A->tda));
  return 0;
}
// ******************************************************
void MVops::cblas_zher2(const enum APL_ORDER order, const enum APL_UPLO Uplo,
                        const int N, const void *alpha, const void *X, const int incX,
                        const void *Y, const int incY, void *A, const int lda) {
  int i, j;
  const int conj = (order == aplColMajor) ? -1 : 1;

  //CHECK_ARGS10(CZ_HER2,order,Uplo,N,alpha,X,incX,Y,incY,A,lda);

  {
    const double alpha_real = CONST_REAL0(alpha);
    const double alpha_imag = CONST_IMAG0(alpha);

    if ((_iszero(alpha_real)) && (_iszero(alpha_imag)))
      return;

    if ((order == aplRowMajor && Uplo == aplUpper) || (order == aplColMajor && Uplo == aplLower)) {
      int ix = OFFSET(N, incX);
      int iy = OFFSET(N, incY);
      for (i = 0; i < N; i++) {
        const double Xi_real = CONST_REAL(X, ix);
        const double Xi_imag = CONST_IMAG(X, ix);
        /* tmp1 = alpha Xi */
        const double tmp1_real = alpha_real * Xi_real - alpha_imag * Xi_imag;
        const double tmp1_imag = alpha_imag * Xi_real + alpha_real * Xi_imag;

        const double Yi_real = CONST_REAL(Y, iy);
        const double Yi_imag = CONST_IMAG(Y, iy);
        /* tmp2 = conj(alpha) Yi */
        const double tmp2_real = alpha_real * Yi_real + alpha_imag * Yi_imag;
        const double tmp2_imag = -alpha_imag * Yi_real + alpha_real * Yi_imag;

        int jx = ix + incX;
        int jy = iy + incY;

        /* Aij = alpha*Xi*conj(Yj) + conj(alpha)*Yi*conj(Xj) */

        REAL(A, lda * i + i) += 2 * (tmp1_real * Yi_real + tmp1_imag * Yi_imag);
        IMAG(A, lda * i + i) = 0;

        for (j = i + 1; j < N; j++) {
          const double Xj_real = CONST_REAL(X, jx);
          const double Xj_imag = CONST_IMAG(X, jx);
          const double Yj_real = CONST_REAL(Y, jy);
          const double Yj_imag = CONST_IMAG(Y, jy);
          REAL(A, lda * i + j) += ((tmp1_real * Yj_real + tmp1_imag * Yj_imag) + (tmp2_real * Xj_real + tmp2_imag * Xj_imag));
          IMAG(A, lda * i + j) +=
              conj * ((tmp1_imag * Yj_real - tmp1_real * Yj_imag) +
                      (tmp2_imag * Xj_real - tmp2_real * Xj_imag));
          jx += incX;
          jy += incY;
        }
        ix += incX;
        iy += incY;
      }
    } else if ((order == aplRowMajor && Uplo == aplLower) || (order == aplColMajor && Uplo == aplUpper)) {
      int ix = OFFSET(N, incX);
      int iy = OFFSET(N, incY);
      for (i = 0; i < N; i++) {
        const double Xi_real = CONST_REAL(X, ix);
        const double Xi_imag = CONST_IMAG(X, ix);
        const double tmp1_real = alpha_real * Xi_real - alpha_imag * Xi_imag;
        const double tmp1_imag = alpha_imag * Xi_real + alpha_real * Xi_imag;

        const double Yi_real = CONST_REAL(Y, iy);
        const double Yi_imag = CONST_IMAG(Y, iy);
        const double tmp2_real = alpha_real * Yi_real + alpha_imag * Yi_imag;
        const double tmp2_imag = -alpha_imag * Yi_real + alpha_real * Yi_imag;

        int jx = OFFSET(N, incX);
        int jy = OFFSET(N, incY);

        /* Aij = alpha*Xi*conj(Yj) + conj(alpha)*Yi*conj(Xj) */

        for (j = 0; j < i; j++) {
          const double Xj_real = CONST_REAL(X, jx);
          const double Xj_imag = CONST_IMAG(X, jx);
          const double Yj_real = CONST_REAL(Y, jy);
          const double Yj_imag = CONST_IMAG(Y, jy);
          REAL(A, lda * i + j) += ((tmp1_real * Yj_real + tmp1_imag * Yj_imag) + (tmp2_real * Xj_real + tmp2_imag * Xj_imag));
          IMAG(A, lda * i + j) +=
              conj * ((tmp1_imag * Yj_real - tmp1_real * Yj_imag) +
                      (tmp2_imag * Xj_real - tmp2_real * Xj_imag));
          jx += incX;
          jy += incY;
        }

        REAL(A, lda * i + i) += 2 * (tmp1_real * Yi_real + tmp1_imag * Yi_imag);
        IMAG(A, lda * i + i) = 0;

        ix += incX;
        iy += incY;
      }
    } else {
      cerr << "\n unrecognized operation \n";
      exit(0);
    }
  }
}
// ******************************************************
template <class T>
int MVops::apl_linalg_hermtd_decomp(apl_matrix_complex<T> *A, apl_vector_complex<T> *tau) {
  if (A->size1 != A->size2) {
    cerr << "\n hermitian tridiagonal decomposition requires square matrix \n";
    exit(0);
  } else if (tau->size + 1 != A->size1) {
    cerr << "\n size of tau must be (matrix size - 1) \n ";
    exit(0);
  } else {
    const size_t N = A->size1;
    size_t i;

    const apl_complex<T> zero = apl_complex_rect(0.0, 0.0);
    const apl_complex<T> one = apl_complex_rect(1.0, 0.0);
    const apl_complex<T> neg_one = apl_complex_rect(-1.0, 0.0);

    for (i = 0; i < N - 1; i++) {
      apl_vector_complex_view<T> c = apl_matrix_complex_column<T>(A, i);
      apl_vector_complex_view<T> v = apl_vector_complex_subvector<T>(&c.vector, i + 1, N - (i + 1));
      apl_complex<T> tau_i = apl_linalg_complex_householder_transform<T>(&v.vector);
      // Apply the transformation H^T A H to the remaining columns //
      if ((i + 1) < (N - 1) && !((_iszero(APL_MV_REAL(tau_i))) && (_iszero(APL_MV_IMAG(tau_i))))) {
        apl_matrix_complex_view<T> m =
            apl_matrix_complex_submatrix<T>(A, i + 1, i + 1,
                                            N - (i + 1), N - (i + 1));
        apl_complex<T> ei = apl_vector_complex_get<T>(&v.vector, 0);
        apl_vector_complex_view<T> x = apl_vector_complex_subvector<T>(tau, i, N - (i + 1));
        apl_vector_complex_set<T>(&v.vector, 0, one);

        // x = tau * A * v //
        apl_blas_zhemv<T>(aplLower, tau_i, &m.matrix, &v.vector, zero, &x.vector);

        // w = x - (1/2) tau * (x' * v) * v  //
        {
          apl_complex<T> xv, txv, alpha;
          apl_blas_zdotc(&x.vector, &v.vector, &xv);
          txv = apl_complex_mul<T>(tau_i, xv);
          alpha = apl_complex_mul_real<T>(txv, -0.5);
          apl_blas_zaxpy<T>(alpha, &v.vector, &x.vector);
        }

        // apply the transformation A = A - v w' - w v' //
        apl_blas_zher2<T>(aplLower, neg_one, &v.vector, &x.vector, &m.matrix);

        apl_vector_complex_set<T>(&v.vector, 0, ei);
      }

      apl_vector_complex_set<T>(tau, i, tau_i);
    }
    return 0;
  }
}
// ******************************************************
template <class T>
int MVops::apl_linalg_complex_householder_hm(apl_complex<T> tau, const apl_vector_complex<T> *v, apl_matrix_complex<T> *A) {
  // applies a householder transformation v,tau to matrix m //

  size_t i, j;

  if ((_iszero(APL_MV_REAL(tau))) && (_iszero(APL_MV_IMAG(tau)))) {
    return 0;
  }

  // w = (v' A)^T //

  for (j = 0; j < A->size2; j++) {
    apl_complex<T> tauwj;
    apl_complex<T> wj = apl_matrix_complex_get<T>(A, 0, j);

    for (i = 1; i < A->size1; i++)  // note, computed for v(0) = 1 above //
    {
      apl_complex<T> Aij = apl_matrix_complex_get<T>(A, i, j);
      apl_complex<T> vi = apl_vector_complex_get<T>(v, i);
      apl_complex<T> Av = apl_complex_mul<T>(Aij, apl_complex_conjugate<T>(vi));
      wj = apl_complex_add<T>(wj, Av);
    }

    tauwj = apl_complex_mul<T>(tau, wj);

    // A = A - v w^T //

    {
      apl_complex<T> A0j = apl_matrix_complex_get<T>(A, 0, j);
      apl_complex<T> Atw = apl_complex_sub<T>(A0j, tauwj);
      // store A0j - tau  * wj //
      apl_matrix_complex_set<T>(A, 0, j, Atw);
    }

    for (i = 1; i < A->size1; i++) {
      apl_complex<T> vi = apl_vector_complex_get<T>(v, i);
      apl_complex<T> tauvw = apl_complex_mul<T>(vi, tauwj);
      apl_complex<T> Aij = apl_matrix_complex_get<T>(A, i, j);
      apl_complex<T> Atwv = apl_complex_sub<T>(Aij, tauvw);
      // store Aij - tau * vi * wj //
      apl_matrix_complex_set<T>(A, i, j, Atwv);
    }
  }

  return 0;
}
// ******************************************************
template <class T>
int MVops::apl_linalg_hermtd_unpack(const apl_matrix_complex<T> *A,
                                    const apl_vector_complex<T> *tau,
                                    apl_matrix_complex<T> *U,
                                    apl_vector<T> *diag,
                                    apl_vector<T> *sdiag) {
  if (A->size1 != A->size2) {
    cerr << "\n  matrix A must be sqaure \n";
    exit(0);
  } else if (tau->size + 1 != A->size1) {
    cerr << "\n size of tau must be (matrix size - 1) \n";
    exit(0);
  } else if (U->size1 != A->size1 || U->size2 != A->size1) {
    cerr << "\n size of U must match size of A \n";
    exit(0);
  } else if (diag->size != A->size1) {
    cerr << "\n size of diagonal must match size of A \n";
    exit(0);
  } else if (sdiag->size + 1 != A->size1) {
    cerr << "\n  size of subdiagonal must be (matrix size - 1) \n ";
    exit(0);
  } else {
    const size_t N = A->size1;

    size_t i;

    // Initialize U to the identity //

    apl_matrix_complex_set_identity<T>(U);

    for (i = N - 1; i-- > 0;) {
      apl_complex<T> ti = apl_vector_complex_get<T>(tau, i);

      const apl_vector_complex_view<T> c = apl_matrix_complex_const_column<T>(A, i);

      const apl_vector_complex_view<T> h =
          apl_vector_complex_const_subvector<T>(&c.vector, i + 1, N - (i + 1));

      apl_matrix_complex_view<T> m =
          apl_matrix_complex_submatrix<T>(U, i + 1, i + 1, N - (i + 1), N - (i + 1));

      apl_linalg_complex_householder_hm(ti, &h.vector, &m.matrix);
    }
    // Copy diagonal into diag //

    for (i = 0; i < N; i++) {
      apl_complex<T> Aii = apl_matrix_complex_get<T>(A, i, i);
      apl_vector_set<T>(diag, i, APL_MV_REAL(Aii));
    }

    // Copy subdiagonal into sdiag //

    for (i = 0; i < N - 1; i++) {
      apl_complex<T> Aji = apl_matrix_complex_get<T>(A, i + 1, i);
      apl_vector_set<T>(sdiag, i, APL_MV_REAL(Aji));
    }
    return 0;
  }
}
// ******************************************************
void MVops::chop_small_elements1(const size_t N, const double d[], double sd[]) {
  double d_i = d[0];

  size_t i;

  for (i = 0; i < N - 1; i++) {
    double sd_i = sd[i];
    double d_ip1 = d[i + 1];

    if (std::fabs(sd_i) < APL_DBL_EPSILON * (std::fabs(d_i) + std::fabs(d_ip1))) {
      sd[i] = 0.0;
    }
    d_i = d_ip1;
  }
}
// ******************************************************
double
MVops::trailing_eigenvalue1(const size_t n, const double d[], const double sd[]) {
  double ta = d[n - 2];
  double tb = d[n - 1];
  double tab = sd[n - 2];

  double dt = (ta - tb) / 2.0;

  double mu;

  if (dt > 0) {
    mu = tb - tab * (tab / (dt + hypot(dt, tab)));
  } else if (_iszero(dt)) {
    mu = tb - fabs(tab);
  } else {
    mu = tb + tab * (tab / ((-dt) + hypot(dt, tab)));
  }

  return mu;
}
// ******************************************************
void MVops::qrstep1(const size_t n, double d[], double sd[], double gc[], double gs[]) {
  double x, z;
  double ak, bk, zk, ap, bp, aq, bq;
  size_t k;

  double mu = trailing_eigenvalue1(n, d, sd);

  // If mu is large relative to d_0 and sd_0 then the Givens rotation
  //   will have no effect, leading to an infinite loop.

  //   We set mu to zero in this case, which at least diagonalises the
  //   submatrix [d_0, sd_0 ; sd_0, d_0] and allows further progress.

  if (APL_DBL_EPSILON * std::fabs(mu) > (std::fabs(d[0]) + std::fabs(sd[0]))) {
    mu = 0;
  }

  x = d[0] - mu;
  z = sd[0];

  ak = 0;
  bk = 0;
  zk = 0;

  ap = d[0];
  bp = sd[0];

  aq = d[1];

  if (n == 2) {
    double c, s;
    create_givens(x, z, &c, &s);

    if (gc != NULL)
      gc[0] = c;
    if (gs != NULL)
      gs[0] = s;

    {
      double ap1 = c * (c * ap - s * bp) + s * (s * aq - c * bp);
      double bp1 = c * (s * ap + c * bp) - s * (s * bp + c * aq);

      double aq1 = s * (s * ap + c * bp) + c * (s * bp + c * aq);

      ak = ap1;
      bk = bp1;

      ap = aq1;
    }

    d[0] = ak;
    sd[0] = bk;
    d[1] = ap;

    return;
  }

  bq = sd[1];

  for (k = 0; k < n - 1; k++) {
    double c, s;
    create_givens(x, z, &c, &s);

    // store Givens rotation //
    if (gc != NULL)
      gc[k] = c;
    if (gs != NULL)
      gs[k] = s;

    // compute G' T G //

    {
      double bk1 = c * bk - s * zk;

      double ap1 = c * (c * ap - s * bp) + s * (s * aq - c * bp);
      double bp1 = c * (s * ap + c * bp) - s * (s * bp + c * aq);
      double zp1 = -s * bq;

      double aq1 = s * (s * ap + c * bp) + c * (s * bp + c * aq);
      double bq1 = c * bq;

      ak = ap1;
      bk = bp1;
      zk = zp1;

      ap = aq1;
      bp = bq1;

      if (k < n - 2)
        aq = d[k + 2];
      if (k < n - 3)
        bq = sd[k + 2];

      d[k] = ak;

      if (k > 0)
        sd[k - 1] = bk1;

      if (k < n - 2)
        sd[k + 1] = bp;

      x = bk;
      z = zk;
    }
  }
  // k = n - 1 //
  d[k] = ap;
  sd[k - 1] = bk;
}
// ******************************************************
template <class T>
int MVops::apl_eigen_hermv(apl_matrix_complex<T> *A, apl_vector<T> *eval,
                           apl_matrix_complex<T> *evec,
                           apl_eigen_hermv_workspace *w) {
  if (A->size1 != A->size2) {
    cerr << "matrix must be square to compute eigenvalues";
    exit(0);
  } else if (eval->size != A->size1) {
    cerr << "eigenvalue vector must match matrix size";
    exit(0);
  } else if (evec->size1 != A->size1 || evec->size2 != A->size1) {
    cerr << "eigenvector matrix must match matrix size";
    exit(0);
  } else {
    const size_t N = A->size1;
    double *const d = w->d;
    double *const sd = w->sd;

    size_t a, b;

    // handle special case //

    if (N == 1) {
      apl_complex<T> A00 = apl_matrix_complex_get<T>(A, 0, 0);
      apl_vector_set<T>(eval, 0, APL_MV_REAL(A00));
      apl_matrix_complex_set<T>(evec, 0, 0, APL_MV_COMPLEX_ONE);
      return 0;
    }

    // Transform the matrix into a symmetric tridiagonal form
    {
      apl_vector_view<double> d_vec = apl_vector_view_array<double>(d, N);
      apl_vector_view<double> sd_vec = apl_vector_view_array<double>(sd, N - 1);
      apl_vector_complex_view<T> tau_vec = apl_vector_complex_view_array<T>(w->tau, N - 1);
      apl_linalg_hermtd_decomp<T>(A, &tau_vec.vector);
      apl_linalg_hermtd_unpack<T>(A, &tau_vec.vector, evec, &d_vec.vector, &sd_vec.vector);
    }

    // Make an initial pass through the tridiagonal decomposition
    //   to remove off-diagonal elements which are effectively zero

    chop_small_elements1(N, d, sd);
    // Progressively reduce the matrix until it is diagonal //

    b = N - 1;

    while (b > 0) {
      if ((_iszero(sd[b - 1])) || (std::isnan(sd[b - 1]))) {
        b--;
        continue;
      }

      // Find the largest unreduced block (a,b) starting from b
      //    and working backwards

      a = b - 1;

      while (a > 0) {
        if (_iszero(sd[a - 1])) {
          break;
        }
        a--;
      }

      {
        size_t i;
        const size_t n_block = b - a + 1;
        double *d_block = d + a;
        double *sd_block = sd + a;
        double *const gc = w->gc;
        double *const gs = w->gs;

        // apply QR reduction with implicit deflation to the
        //   unreduced block //

        qrstep1(n_block, d_block, sd_block, gc, gs);

        // Apply  Givens rotation Gij(c,s) to matrix Q,  Q <- Q G //

        for (i = 0; i < n_block - 1; i++) {
          const double c = gc[i], s = gs[i];
          size_t k;

          for (k = 0; k < N; k++) {
            apl_complex<T> qki = apl_matrix_complex_get<T>(evec, k, a + i);
            apl_complex<T> qkj = apl_matrix_complex_get<T>(evec, k, a + i + 1);
            // qki <= qki * c - qkj * s //
            // qkj <= qki * s + qkj * c //
            apl_complex<T> x1 = apl_complex_mul_real<T>(qki, c);
            apl_complex<T> y1 = apl_complex_mul_real<T>(qkj, -s);

            apl_complex<T> x2 = apl_complex_mul_real<T>(qki, s);
            apl_complex<T> y2 = apl_complex_mul_real<T>(qkj, c);

            apl_complex<T> qqki = apl_complex_add<T>(x1, y1);
            apl_complex<T> qqkj = apl_complex_add<T>(x2, y2);

            apl_matrix_complex_set<T>(evec, k, a + i, qqki);
            apl_matrix_complex_set<T>(evec, k, a + i + 1, qqkj);
          }
        }

        // remove any small off-diagonal elements //

        chop_small_elements1(n_block, d_block, sd_block);
      }
    }

    {
      apl_vector_view<T> d_vec = apl_vector_view_array<T>(d, N);
      apl_vector_memcpy<T>(eval, &d_vec.vector);
    }

    return 0;
  }
}
// ******************************************************
bool MVops::isorthogonal(apl_matrix_complex<double> *m) {
  size_t nBranches = m->size1;
  for (size_t i = 0; i != nBranches; i++) {
    apl_vector_complex_view<double> v_i = apl_matrix_complex_column<double>(m, i);
    for (size_t j = 0; j != nBranches; j++) {
      apl_vector_complex_view<double> v_j = apl_matrix_complex_column<double>(m, j);
      apl_complex<double> vdot;
      apl_blas_zdotc<double>(&v_i.vector, &v_j.vector, &vdot);  //conj(xT)*y

      double re = APL_MV_REAL(vdot);
      double im = APL_MV_IMAG(vdot);
      if (i == j) {
        if (!((_isequal(re, 1.0)) & (_iszero(im)))) return false;
      } else {
        if (!((_iszero(re)) & (_iszero(im)))) return false;
      }
    }
  }
  return true;
}
// ******************************************************
int MVops::apl_eigen_hermv_sort(apl_vector<double> *eval, apl_matrix_complex<double> *evec,
                                apl_eigen_sort_t sort_type) {
  if (evec->size1 != evec->size2) {
    cerr << "eigenvector matrix must be square \n";
    exit(0);
  } else if (eval->size != evec->size1) {
    cerr << "eigenvalues must match eigenvector matrix \n";
    exit(0);
  } else {
    const size_t N = eval->size;
    size_t i;

    for (i = 0; i < N - 1; i++) {
      size_t j;
      size_t k = i;

      double ek = apl_vector_get<double>(eval, i);

      // search for something to swap //
      for (j = i + 1; j < N; j++) {
        int test;
        const double ej = apl_vector_get<double>(eval, j);

        switch (sort_type) {
          case APL_MV_EIGEN_SORT_VAL_ASC:
            test = (ej < ek);
            break;
          case APL_MV_EIGEN_SORT_VAL_DESC:
            test = (ej > ek);
            break;
          case APL_MV_EIGEN_SORT_ABS_ASC:
            test = (fabs(ej) < fabs(ek));
            break;
          case APL_MV_EIGEN_SORT_ABS_DESC:
            test = (std::fabs(ej) > std::fabs(ek));
            break;
          default:
            cerr << "unrecognized sort type \n ";
            exit(0);
        }

        if (test) {
          k = j;
          ek = ej;
        }
      }

      if (k != i) {
        // swap eigenvalues //
        apl_vector_swap_elements<double>(eval, i, k);

        // swap eigenvectors //
        apl_matrix_complex_swap_columns<double>(evec, i, k);
      }
    }

    return 0;
  }
}
// ******************************************************
void aplEigensystems::eigen_calculation(const aurostd::xmatrix<xcomplex<double> > &M,
                                        aurostd::xvector<double> &apleval,
                                        aurostd::xmatrix<xcomplex<double> > &aplevec) {
  size_t nBranches = M.rows;
  apl_matrix_complex<double> *m = apl_matrix_complex_alloc<double>(nBranches, nBranches);
  apl_matrix_complex<double> *evec = apl_matrix_complex_alloc<double>(nBranches, nBranches);
  apl_vector<double> *eval = apl_vector_alloc<double>(nBranches);

  apl_complex<double> z;
  for (size_t i = 0; i < nBranches; i++) {
    for (size_t j = 0; j < nBranches; j++) {
      APL_MV_SET_COMPLEX(&z, M[i + 1][j + 1].re, M[i + 1][j + 1].im);
      apl_matrix_complex_set<double>(m, i, j, z);
    }
  }
  apl_eigen_hermv_workspace *w = apl_eigen_hermv_alloc(nBranches);
  if (apl_eigen_hermv<double>(m, eval, evec, w) != 0) {
    cout << "\n apl eigenvale is not working" << std::endl;
    exit(0);
  };
  //apl_eigen_hermv_sort (eval, evec, APL_MV_EIGEN_SORT_ABS_ASC);
  apl_eigen_hermv_free(w);

  //if(!isorthogonal(evec))cerr<<"eigenvectors are not orthogonal "<<"\n";
  //PRINTM(evec);

  //copy apl eigenvalues and eigenvectors to CPP variables
  for (uint i = 0; i != nBranches; i++) {
    apleval[i + 1] = apl_vector_get<double>(eval, i);
    for (uint j = 0; j != nBranches; j++) {
      z = apl_matrix_complex_get<double>(evec, i, j);
      aplevec[i + 1][j + 1] = xcomplex<double>(APL_MV_REAL(z), APL_MV_IMAG(z));
    }
  }

  apl_matrix_complex_free<double>(m);
  apl_matrix_complex_free<double>(evec);
  apl_vector_free<double>(eval);
}
// ******************************************************
void aplEigensystems::eigen_calculation(const aurostd::xmatrix<xcomplex<double> > &M,
                                        aurostd::xvector<double> &apleval,
                                        aurostd::xmatrix<xcomplex<double> > &aplevec,
                                        apl_eigen_sort_t t) {
  size_t nBranches = M.rows;
  apl_matrix_complex<double> *m = apl_matrix_complex_alloc<double>(nBranches, nBranches);
  apl_matrix_complex<double> *evec = apl_matrix_complex_alloc<double>(nBranches, nBranches);
  apl_vector<double> *eval = apl_vector_alloc<double>(nBranches);

  apl_complex<double> z;
  for (size_t i = 0; i < nBranches; i++) {
    for (size_t j = 0; j < nBranches; j++) {
      APL_MV_SET_COMPLEX(&z, M[i + 1][j + 1].re, M[i + 1][j + 1].im);
      apl_matrix_complex_set<double>(m, i, j, z);
    }
  }
  apl_eigen_hermv_workspace *w = apl_eigen_hermv_alloc(nBranches);
  if (apl_eigen_hermv<double>(m, eval, evec, w) != 0) {
    cout << "\n apl eigenvale is not working" << std::endl;
    exit(0);
  };

  apl_eigen_hermv_sort(eval, evec, t);
  apl_eigen_hermv_free(w);

  //if(!isorthogonal(evec))cerr<<"eigenvectors are not orthogonal "<<"\n";
  //PRINTM(evec);

  //copy apl eigenvalues and eigenvectors to CPP variables
  for (uint i = 0; i != nBranches; i++) {
    apleval[i + 1] = apl_vector_get<double>(eval, i);
    for (uint j = 0; j != nBranches; j++) {
      z = apl_matrix_complex_get<double>(evec, i, j);
      aplevec[i + 1][j + 1] = xcomplex<double>(APL_MV_REAL(z), APL_MV_IMAG(z));
    }
  }

  apl_matrix_complex_free<double>(m);
  apl_matrix_complex_free<double>(evec);
  apl_vector_free<double>(eval);
}
// ******************************************************
int MVops::apl_blas_zgemv(APL_TRANSPOSE_t TransA, const apl_complex<double> alpha,
                          const apl_matrix_complex<double> *A, const apl_vector_complex<double> *X,
                          const apl_complex<double> beta, apl_vector_complex<double> *Y) {
  const size_t M = A->size1;
  const size_t N = A->size2;

  if ((TransA == aplNoTrans && N == X->size && M == Y->size) || (TransA == aplTrans && M == X->size && N == Y->size) || (TransA == aplConjTrans && M == X->size && N == Y->size)) {
    cblas_zgemv(aplRowMajor, TransA, int(M), int(N),
                APL_MV_COMPLEX_P(&alpha), A->data, int(A->tda), X->data,
                int(X->stride), APL_MV_COMPLEX_P(&beta), Y->data,
                int(Y->stride));
    return 0;
  } else {
    cerr << "invalid length"
         << "\n";
    exit(0);
  }
}
// ******************************************************
void MVops::cblas_zgemv(const enum APL_ORDER order, const enum APL_TRANSPOSE TransA,
                        const int M, const int N, const void *alpha, const void *A,
                        const int lda, const void *X, const int incX, const void *beta,
                        void *Y, const int incY) {
  {
    int i, j;
    int lenX, lenY;

    const double alpha_real = CONST_REAL0(alpha);
    const double alpha_imag = CONST_IMAG0(alpha);

    const double beta_real = CONST_REAL0(beta);
    const double beta_imag = CONST_IMAG0(beta);

    //CHECK_ARGS12(GEMV,order,TransA,M,N,alpha,A,lda,X,incX,beta,Y,incY);

    if (M == 0 || N == 0)
      return;

    if (((_iszero(alpha_real)) && (_iszero(alpha_imag))) && ((_isequal(beta_real, 1.0)) && (_iszero(beta_imag))))
      return;

    if (TransA == aplNoTrans) {
      lenX = N;
      lenY = M;
    } else {
      lenX = M;
      lenY = N;
    }

    // form  y := beta*y

    if ((_iszero(beta_real)) && (_iszero(beta_imag))) {
      int iy = OFFSET(lenY, incY);
      for (i = 0; i < lenY; i++) {
        REAL(Y, iy) = 0.0;
        IMAG(Y, iy) = 0.0;
        iy += incY;
      }
    } else if (!((_isequal(beta_real, 1.0)) && (_iszero(beta_imag)))) {
      int iy = OFFSET(lenY, incY);
      for (i = 0; i < lenY; i++) {
        const double y_real = REAL(Y, iy);
        const double y_imag = IMAG(Y, iy);
        const double tmpR = y_real * beta_real - y_imag * beta_imag;
        const double tmpI = y_real * beta_imag + y_imag * beta_real;
        REAL(Y, iy) = tmpR;
        IMAG(Y, iy) = tmpI;
        iy += incY;
      }
    }
    if ((_iszero(alpha_real)) && (_iszero(alpha_imag)))
      return;

    if ((order == aplRowMajor && TransA == aplNoTrans) || (order == aplColMajor && TransA == aplTrans)) {
      // form  y := alpha*A*x + y //
      int iy = OFFSET(lenY, incY);
      for (i = 0; i < lenY; i++) {
        double dotR = 0.0;
        double dotI = 0.0;
        int ix = OFFSET(lenX, incX);
        for (j = 0; j < lenX; j++) {
          const double x_real = CONST_REAL(X, ix);
          const double x_imag = CONST_IMAG(X, ix);
          const double A_real = CONST_REAL(A, lda * i + j);
          const double A_imag = CONST_IMAG(A, lda * i + j);

          dotR += A_real * x_real - A_imag * x_imag;
          dotI += A_real * x_imag + A_imag * x_real;
          ix += incX;
        }

        REAL(Y, iy) += alpha_real * dotR - alpha_imag * dotI;
        IMAG(Y, iy) += alpha_real * dotI + alpha_imag * dotR;
        iy += incY;
      }
    } else if ((order == aplRowMajor && TransA == aplTrans) || (order == aplColMajor && TransA == aplNoTrans)) {
      // form  y := alpha*A'*x + y //
      int ix = OFFSET(lenX, incX);
      for (j = 0; j < lenX; j++) {
        double x_real = CONST_REAL(X, ix);
        double x_imag = CONST_IMAG(X, ix);
        double tmpR = alpha_real * x_real - alpha_imag * x_imag;
        double tmpI = alpha_real * x_imag + alpha_imag * x_real;

        int iy = OFFSET(lenY, incY);
        for (i = 0; i < lenY; i++) {
          const double A_real = CONST_REAL(A, lda * j + i);
          const double A_imag = CONST_IMAG(A, lda * j + i);
          REAL(Y, iy) += A_real * tmpR - A_imag * tmpI;
          IMAG(Y, iy) += A_real * tmpI + A_imag * tmpR;
          iy += incY;
        }
        ix += incX;
      }
    } else if (order == aplRowMajor && TransA == aplConjTrans) {
      // form  y := alpha*A^H*x + y //
      int ix = OFFSET(lenX, incX);
      for (j = 0; j < lenX; j++) {
        double x_real = CONST_REAL(X, ix);
        double x_imag = CONST_IMAG(X, ix);
        double tmpR = alpha_real * x_real - alpha_imag * x_imag;
        double tmpI = alpha_real * x_imag + alpha_imag * x_real;

        int iy = OFFSET(lenY, incY);
        for (i = 0; i < lenY; i++) {
          const double A_real = CONST_REAL(A, lda * j + i);
          const double A_imag = CONST_IMAG(A, lda * j + i);
          REAL(Y, iy) += A_real * tmpR - (-A_imag) * tmpI;
          IMAG(Y, iy) += A_real * tmpI + (-A_imag) * tmpR;
          iy += incY;
        }
        ix += incX;
      }
    } else if (order == aplColMajor && TransA == aplConjTrans) {
      // form  y := alpha*A^H*x + y //
      int iy = OFFSET(lenY, incY);
      for (i = 0; i < lenY; i++) {
        double dotR = 0.0;
        double dotI = 0.0;
        int ix = OFFSET(lenX, incX);
        for (j = 0; j < lenX; j++) {
          const double x_real = CONST_REAL(X, ix);
          const double x_imag = CONST_IMAG(X, ix);
          const double A_real = CONST_REAL(A, lda * i + j);
          const double A_imag = CONST_IMAG(A, lda * i + j);

          dotR += A_real * x_real - (-A_imag) * x_imag;
          dotI += A_real * x_imag + (-A_imag) * x_real;
          ix += incX;
        }

        REAL(Y, iy) += alpha_real * dotR - alpha_imag * dotI;
        IMAG(Y, iy) += alpha_real * dotI + alpha_imag * dotR;
        iy += incY;
      }
    } else {
      cerr << "unrecognized operation"
           << "\n";
      exit(0);
    }
  }
}
// ******************************************************
xvector<double>
MVops::calculate_gruneisen(xmatrix<xcomplex<double> > &m0, xmatrix<xcomplex<double> > &mP, xmatrix<xcomplex<double> > &mM, double delV, double V_0) {
  //gruneisen
  size_t nBranches = (size_t)m0.rows;
  apl_matrix_complex<double> *aplDM0 = NULL;
  apl_matrix_complex<double> *aplDM1 = NULL;
  apl_matrix_complex<double> *aplDM2 = NULL;
  apl_matrix_complex<double> *ddm = NULL;
  aplDM0 = apl_matrix_complex_alloc<double>(nBranches, nBranches);
  aplDM1 = apl_matrix_complex_alloc<double>(nBranches, nBranches);
  aplDM2 = apl_matrix_complex_alloc<double>(nBranches, nBranches);
  ddm = apl_matrix_complex_alloc<double>(nBranches, nBranches);
  //apl Eigenvectors
  apl_matrix_complex<double> *aplevec0 = NULL;
  aplevec0 = apl_matrix_complex_alloc<double>(nBranches, nBranches);

  //apl Eigenvalues
  apl_vector<double> *apleval0 = NULL;
  apleval0 = apl_vector_alloc<double>(nBranches);

  //apl mat set
  for (size_t i = 0; i < nBranches; i++) {
    for (size_t j = 0; j < nBranches; j++) {
      apl_matrix_complex_set<double>(aplDM0, i, j, apl_complex_rect(m0[i + 1][j + 1].re, m0[i + 1][j + 1].im));
      apl_matrix_complex_set<double>(aplDM1, i, j, apl_complex_rect(mP[i + 1][j + 1].re, mP[i + 1][j + 1].im));
      apl_matrix_complex_set<double>(aplDM2, i, j, apl_complex_rect(mM[i + 1][j + 1].re, mM[i + 1][j + 1].im));
    }
  }

  apl_eigen_hermv_workspace *w = apl_eigen_hermv_alloc(nBranches);
  if (apl_eigen_hermv<double>(aplDM0, apleval0, aplevec0, w) != 0) {
    cout << "\n apl eigenvale is not working" << std::endl;
    exit(0);
  };
  apl_eigen_hermv_sort(apleval0, aplevec0, APL_MV_EIGEN_SORT_ABS_ASC);
  apl_eigen_hermv_free(w);

  //DM derivative
  for (size_t i = 0; i < nBranches; i++) {
    for (size_t j = 0; j < nBranches; j++) {
      apl_complex<double> mij = apl_matrix_complex_get<double>(aplDM1, i, j);
      double a1 = APL_MV_REAL(mij);
      double b1 = APL_MV_IMAG(mij);

      apl_complex<double> nij = apl_matrix_complex_get<double>(aplDM2, i, j);
      double a2 = APL_MV_REAL(nij);
      double b2 = APL_MV_IMAG(nij);

      double da = (a1 - a2) / (2. * delV);
      double db = (b1 - b2) / (2. * delV);
      apl_matrix_complex_set<double>(ddm, i, j, apl_complex_rect(da, db));
    }
  }

  //matrix-vector operations
  apl_complex<double> ONE, ZERO;
  APL_MV_SET_COMPLEX(&ONE, 1., 0.);
  APL_MV_SET_COMPLEX(&ZERO, 0., 0.);

  xvector<double> GP(nBranches, 1);
  vector<double> GPv(nBranches);
  for (size_t i = 0; i < nBranches; i++) {
    apl_vector_complex_view<double> vi = apl_matrix_complex_column<double>(aplevec0, i);
    apl_vector_complex<double> *out = apl_vector_complex_alloc<double>(nBranches);
    apl_blas_zgemv(aplNoTrans, ONE, ddm, &vi.vector, ZERO, out);

    apl_complex<double> vdot;
    apl_blas_zdotc<double>(&vi.vector, out, &vdot);

    double gpr = APL_MV_REAL(vdot);
    double gpi = APL_MV_IMAG(vdot);
    double eval_i = apl_vector_get<double>(apleval0, i);

    //Gruneisen calculations
    gpr = -(V_0 * gpr) / (2.0 * eval_i);
    gpi = -(V_0 * gpi) / (2.0 * eval_i);

    if (std::isnan(gpr)) gpr = 0.0;
    if (std::isinf(gpr)) gpr = 0.00;

    //Big numbers in Gruneisen probably due to numerical erors
    //avoid small number arithmatics
    if (abs(gpr) > 10) {
      if (_iszero(eval_i)) {
        eval_i = 0.0;
        gpr = 0.0;
      }
      if (_iszero(APL_MV_REAL(vdot))) gpr = 0.0;
      //if(abs(APL_MV_REAL(vdot))<1e-4 && (abs(eval_i)<1e-4) ) gpr=0.0;
      //else if( (abs(gpr)>10) && (abs(eval_i)<1e-3) ) gpr=0.0;
    }
    //ignore small negative eigenvalues, big negative eigenvalues assuming is already screened before this point
    if (eval_i < 0.00) {
      eval_i = 0.0;
      gpr = 0.0;
    }
    GPv[i] = gpr;
    apl_vector_complex_free<double>(out);
  }

  //if calculation are in high-symmetric qpoints then sort Gruneisen
  //if(ktype=="path")std::sort(GPv.begin(), GPv.end());
  //store to xvector
  for (uint i = 0; i != GPv.size(); i++) GP[i + 1] = GPv[i];

  apl_matrix_complex_free<double>(ddm);
  apl_matrix_complex_free<double>(aplDM0);
  apl_matrix_complex_free<double>(aplDM1);
  apl_matrix_complex_free<double>(aplDM2);
  apl_matrix_complex_free<double>(aplevec0);
  apl_vector_free<double>(apleval0);
  return GP;
}
// ******************************************************
xvector<int>
MVops::trace_acoustic_mode(const xmatrix<xcomplex<double> > &mk1, const xmatrix<xcomplex<double> > &mk2) {
  size_t nBranches = (size_t)mk1.rows;
  xvector<int> max_overlap_index(nBranches, 1);

  apl_matrix_complex<double> *m1 = apl_matrix_complex_alloc<double>(nBranches, nBranches);
  apl_matrix_complex<double> *m2 = apl_matrix_complex_alloc<double>(nBranches, nBranches);

  //apl mat set
  for (size_t i = 0; i < nBranches; i++) {
    for (size_t j = 0; j < nBranches; j++) {
      apl_matrix_complex_set<double>(m1, i, j, apl_complex_rect(mk1[i + 1][j + 1].re, mk1[i + 1][j + 1].im));
      apl_matrix_complex_set<double>(m2, i, j, apl_complex_rect(mk2[i + 1][j + 1].re, mk2[i + 1][j + 1].im));
    }
  }

  //eigenvectors overlap calculations
  for (size_t i = 0; i < nBranches; i++) {
    apl_vector_complex_view<double> vk1 = apl_matrix_complex_column<double>(m1, i);

    vector<double> abs_values(nBranches, 0.0);
    for (uint j = 0; j < nBranches; j++) {
      apl_vector_complex_view<double> vk2 = apl_matrix_complex_column<double>(m2, j);

      apl_complex<double> vdot;
      apl_blas_zdotc<double>(&vk1.vector, &vk2.vector, &vdot);  //conj(xT)*y
      double vdot_r = APL_MV_REAL(vdot);
      double vdot_i = APL_MV_IMAG(vdot);
      std::complex<double> eigen_pro = std::complex<double>(vdot_r, vdot_i);
      abs_values[j] = std::abs(eigen_pro);
    }
    //store maximum maximum operlap
    double max = abs_values[0];
    uint index = 0;
    for (uint ii = 0; ii < nBranches; ii++) {
      if (max < abs_values[ii]) {
        max = abs_values[ii];
        index = ii;
      }
    }
    max_overlap_index[i + 1] = (int)index;
  }
  apl_matrix_complex_free<double>(m1);
  apl_matrix_complex_free<double>(m2);

  return max_overlap_index;
}

// ******************************************************
}
