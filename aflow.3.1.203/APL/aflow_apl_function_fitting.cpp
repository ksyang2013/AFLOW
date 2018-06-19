#include "aflow_apl.h"

using namespace std;
namespace apl {
bool aflowFITTING::quadraticfit(xvector<double> energy, xvector<double> volume, xvector<double> &params) {
  uint NPT = energy.rows;
  uint NTERM = 3;
  xvector<int> ia(NTERM, 1);
  xvector<double> a(NTERM, 1);
  xvector<double> sig(NPT, 1);
  xmatrix<double> covar(NTERM, NTERM, 1, 1);
  for (uint i = 1; i <= NPT; i++) {
    funcs(volume[i], a);
    sig[i] = FIT_LIMIT;
  }
  for (uint i = 1; i <= NTERM; i++) ia[i] = 1;
  double chisq = 0.;
  void (aflowFITTING::*quadratic)(double, xvector<double> &) = &aflowFITTING::funcs;

  if (lfit(volume, energy, sig, a, ia, covar, chisq, quadratic)) {
    params = a;
    chisq_quadratic = chisq;
    return true;
  } else {
    return false;
  }
}
// ***************************************************
bool aflowFITTING::birch_murnaghan_fitting(xvector<double> energy,
                                           xvector<double> volume,
                                           xvector<double> gues, xvector<double> &params) {
  uint MA = 4;  // number of undetermined coefficients
  Uncertainties_birch_murnaghan.clear();
  Uncertainties_birch_murnaghan.resize(4, 0.0);
  uint NPT = energy.rows;
  int iter, itst, mfit = MA;
  uint Iteration = 0;
  double alamda, chisq_bm, ochisq;
  static xvector<double> a_bm(MA, 1);

  xvector<int> ia_bm(MA, 1);
  xvector<double> sig_bm(NPT, 1);
  xmatrix<double> covar_bm(MA, MA, 1, 1);
  xmatrix<double> alpha(MA, MA, 1, 1);
  for (uint i = 1; i <= NPT; i++) sig_bm[i] = FIT_LIMIT;

  void (aflowFITTING::*myf)(double, const xvector<double>, double *, xvector<double> &) = &aflowFITTING::birch_murnaghan_function;

  for (int i = 1; i <= mfit; i++) ia_bm[i] = 1;
  for (uint i = 1; i <= MA; i++) a_bm[i] = gues[i];
  for (iter = 1; iter <= 2; iter++) {
    alamda = -1;
    if (!mrqmin(volume, energy, sig_bm, a_bm, ia_bm, covar_bm, alpha, chisq_bm, myf, &alamda)) return false;
    itst = 0;
    Iteration = 1;
    for (;;) {
      Iteration++;
      ochisq = chisq_bm;
      if (!mrqmin(volume, energy, sig_bm, a_bm, ia_bm, covar_bm, alpha, chisq_bm, myf, &alamda)) return false;
      if (chisq_bm > ochisq)
        itst = 0;
      else if (fabs(ochisq - chisq_bm) < 0.1)
        itst++;
      if (itst < 4) continue;
      alamda = 0.0;
      mrqmin(volume, energy, sig_bm, a_bm, ia_bm, covar_bm, alpha, chisq_bm, myf, &alamda);
      break;
    }
  }
  Iteration_birch_murnaghan = Iteration;
  alamda_birch_murnaghan = alamda;
  params = a_bm;
  chisq_birch_murnaghan = chisq_bm;
  for (uint i = 1; i <= MA; i++) Uncertainties_birch_murnaghan[i - 1] = sqrt(covar_bm[i][i]);
  return true;
}  //function end
// ***************************************************
template <class utype>
void aflowFITTING::covsrt(xmatrix<utype> &covar,
                          xvector<int> &ia, int &mfit) {
  // Expand in storage the covariance matrix covar[1,ma][1,ma], so as to take into account parameters
  // that are being fixed. (for the latter, return zero covariance.)

  int i, j, k;
  utype temp;

  int ma = covar.rows;
  if (covar.cols != covar.rows) {
    cerr << "covsrt: covar.cols!=covar.rows" << std::endl;
    return;
  }

  for (i = mfit + 1; i <= ma; i++)
    for (j = 1; j <= i; j++) covar[i][j] = covar[j][i] = 0.0;
  k = mfit;
  for (j = ma; j >= 1; j--) {
    if (ia[j]) {
      for (i = 1; i <= ma; i++) { SWAP(covar[i][k], covar[i][j]); }
      for (i = 1; i <= ma; i++) { SWAP(covar[k][i], covar[j][i]); }
      k--;
    }
  }
}
// ***************************************************
template <class utype>
bool aflowFITTING::gaussj(xmatrix<utype> &a,
                          int &n, xmatrix<utype> &b, int m) {
  // linear equation solution by gauss-jordan elimination, a[1,n][1,n] is the input matrix.
  // b[1,n][1,m] is input containing the m right-hand side vectors. On the output a is replaced
  // by its matrix inverse, and b is replaced by the corresponding set of solution vectors.

  int i, icol, irow, j, k, l, ll;
  utype big, dum, pivinv, temp;

  // int n=a.rows,m=b.cols;
  // if(a.cols!=a.rows) {cerr << _AUROSTD_XLIBS_ERROR_ << "gaussj: a.cols!=a.rows" << endl; exit(0);}
  // if(b.rows!=a.rows) {cerr << _AUROSTD_XLIBS_ERROR_ << "gaussj: b.rows!=a.rows" << endl; exit(0);}

  if (n > a.rows) {
    cerr << "gaussj: n>a.rows" << std::endl;
    return false;
  }
  if (n > b.rows) {
    cerr << "gaussj: n>b.rows" << std::endl;
    return false;
  }
  if (m > b.cols) {
    cerr << "gaussj: m>b.cols" << std::endl;
    return false;
  }

  xvector<int> indxc(n, 1);
  xvector<int> indxr(n, 1);
  xvector<int> ipiv(n, 1);

  for (j = 1; j <= n; j++) ipiv[j] = 0;
  for (i = 1; i <= n; i++) {
    big = 0.0;
    for (j = 1; j <= n; j++)
      if (ipiv[j] != 1)
        for (k = 1; k <= n; k++) {
          if (ipiv[k] == 0) {
            if (fabs(a[j][k]) >= big) {
              big = fabs(a[j][k]);
              irow = j;
              icol = k;
            }
          }
        }
    ++(ipiv[icol]);
    if (irow != icol) {
      for (l = 1; l <= n; l++) SWAP(a[irow][l], a[icol][l])
      for (l = 1; l <= m; l++) SWAP(b[irow][l], b[icol][l])
    }
    indxr[i] = irow;
    indxc[i] = icol;
    if (a[icol][icol] == 0.0) {
      cerr << "gaussj: Singular Matrix-1" << std::endl;
      return false;
    }
    pivinv = 1.0 / a[icol][icol];
    a[icol][icol] = 1.0;
    for (l = 1; l <= n; l++) a[icol][l] *= pivinv;
    for (l = 1; l <= m; l++) b[icol][l] *= pivinv;
    for (ll = 1; ll <= n; ll++)
      if (ll != icol) {
        dum = a[ll][icol];
        a[ll][icol] = 0.0;
        for (l = 1; l <= n; l++) a[ll][l] -= a[icol][l] * dum;
        for (l = 1; l <= m; l++) b[ll][l] -= b[icol][l] * dum;
      }
  }
  for (l = n; l >= 1; l--) {
    if (indxr[l] != indxc[l])
      for (k = 1; k <= n; k++)
        SWAP(a[k][indxr[l]], a[k][indxc[l]]);
  }
  return true;
}
// ***************************************************
template <class utype>
bool aflowFITTING::lfit(xvector<utype> &x,
                        xvector<utype> &y,
                        xvector<utype> &sig,
                        xvector<utype> &a,
                        xvector<int> &ia,
                        xmatrix<utype> &covar,
                        utype &chisq, void (aflowFITTING::*funcs)(utype, xvector<utype> &)) {
  // Given a set of data points x[1,ndat],y[1,ndat] with individual standar deviation sig[1,ndat], use chisq minimization to fit for some or all the coefficients
  // a[1,ma] of a function that depends linearly on a, y=sum_i a_i*afunc_i(x). The input array ia[1,ma] indicates by nonzero entries those componends of a
  // that should be fitted for, and by zero entries those components that should be held fiuxed at their input values.
  // The prgram returns value for a[1,ma], chisq,  and the covariance atrix covar[1,ma][1,ma]. (Parameters held fixed will return zero covariances.).
  // The user supplies a routine funcs(x,xvector<afunc>) that returns the ma basis funcions evaluated at x=X in the array afunc[1,ma]
  int ndat = x.rows;

  if (y.rows != x.rows) {
    cerr << "\n lfit: y.rows!=x.rows" << std::endl;
    return false;
  }
  if (sig.rows != x.rows) {
    cerr << "\n lfit: sig.rows!=x.rows" << std::endl;
    return false;
  }

  int ma = a.rows;
  if (ia.rows != a.rows) {
    cerr << _AUROSTD_XLIBS_ERROR_ << "lfit: ia.rows!=a.rows" << std::endl;
    return false;
  }

  int i, j, k, l, m, mfit = 0;
  utype ym, wt, sum, sig2i;

  xmatrix<utype> beta(ma, 1, 1, 1);
  xvector<utype> afunc(ma, 1);

  for (j = 1; j <= ma; j++)
    if (ia[j]) mfit++;
  if (mfit == 0) {
    cerr << "\n lfit: no parameters to be fitted" << std::endl;
    return false;
  }
  for (j = 1; j <= mfit; j++) {
    for (k = 1; k <= mfit; k++) covar[j][k] = 0.0;
    beta[j][1] = 0.0;
  }
  for (i = 1; i <= ndat; i++) {
    (*this.*funcs)(x[i], afunc);  //(this->funcs)(x[i],afunc);
    ym = y[i];
    if (mfit < ma) {
      for (j = 1; j <= ma; j++)
        if (!ia[j]) ym -= a[j] * afunc[j];
    }
    sig2i = 1.0 / sqrt(sig[i]);
    for (j = 0, l = 1; l <= ma; l++) {
      if (ia[l]) {
        wt = afunc[l] * sig2i;
        for (j++, k = 0, m = 1; m <= l; m++)
          if (ia[m]) covar[j][++k] += wt * afunc[m];
        beta[j][1] += ym * wt;
      }
    }
  }
  for (j = 2; j <= mfit; j++)
    for (k = 1; k < j; k++)
      covar[k][j] = covar[j][k];
  if (!gaussj(covar, mfit, beta, 1)) return false;
  for (j = 0, l = 1; l <= ma; l++)
    if (ia[l]) a[l] = beta[++j][1];

  chisq = 0.0;
  for (i = 1; i <= ndat; i++) {
    (this->funcs)(x[i], afunc);
    for (sum = 0.0, j = 1; j <= ma; j++) sum += a[j] * afunc[j];
    chisq += ((y[i] - sum) / sig[i]) * ((y[i] - sum) / sig[i]);
  }
  covsrt(covar, ia, mfit);
  return true;
}
// ***************************************************
void aflowFITTING::funcs(double x, xvector<double> &afunc) {
  afunc[1] = 1.0;
  afunc[2] = x;
  afunc[3] = x * x;
}
// ***************************************************
void aflowFITTING::birch_murnaghan_function(double x, const xvector<double> a, double *y, xvector<double> &dyda) {
  double V, Eo, Bo, Vo, Bp;
  *y = 0.0;
  Eo = a[1];
  Bo = a[2];
  Vo = a[3];
  Bp = a[4];
  V = x;

  double myf = Eo + (Bo * V) / Bp * (pow((Vo / V), Bp) / (Bp - 1.) + 1.) - Vo * Bo / (Bp - 1.);
  *y = myf;

  double D_myf_Eo = 1.0;
  double D_myf_Bo = -(Vo / (-1. + Bp)) + (V * (1. + pow((Vo / V), Bp) / (-1. + Bp))) / Bp;
  double D_myf_Vo = -(Bo / (-1. + Bp)) + (Bo * pow((Vo / V), (-1.0 + Bp))) / (-1.0 + Bp);
  double D_myf_Bp = (Bo * Vo) / ((-1. + Bp) * (-1. + Bp)) - (Bo * V * (1. + pow((Vo / V), Bp) / (-1. + Bp))) / (Bp * Bp) +
                    (Bo * V * (-(pow((Vo / V), Bp) / (-1. + Bp) * (-1. + Bp)) + (pow((Vo / V), Bp) * log(Vo / V)) / (-1.0 + Bp))) / Bp;

  dyda[1] = D_myf_Eo;
  dyda[2] = D_myf_Bo;
  dyda[3] = D_myf_Vo;
  dyda[4] = D_myf_Bp;
}
// ***************************************************
void aflowFITTING::mrqcof(xvector<double> x,
                          xvector<double> y,
                          xvector<double> &sig,
                          xvector<double> &a,
                          xvector<int> &ia,
                          xmatrix<double> &alpha,
                          xvector<double> &beta,
                          double &chisq,
                          void (aflowFITTING::*birch_murnaghan_function)(double, xvector<double>,
                                                                         double *, xvector<double> &)) {
  uint ndata = x.rows;

  if (y.rows != x.rows) {
    cerr << " mrqcof: y.rows!=x.rows" << std::endl;
    return;
  }
  if (sig.rows != x.rows) {
    cerr << "mrqcof: sig.rows!=x.rows" << std::endl;
    return;
  }

  uint ma = a.rows;
  if (ia.rows != a.rows) {
    cerr << "mrqcof: ia.rows!=a.rows" << std::endl;
    return;
  }

  uint i, j, k, l, m, mfit = 0;
  double ymod, wt, sig2i, dy;

  xvector<double> dyda(ma, 1);

  for (j = 1; j <= ma; j++)
    if (ia[j]) mfit++;
  for (j = 1; j <= mfit; j++) {
    for (k = 1; k <= j; k++) alpha[j][k] = 0.0;
    beta[j] = 0.0;
  }
  chisq = 0.0;
  for (i = 1; i <= ndata; i++) {
    //(this->birch_murnaghan_function)(x[i],a,&ymod,dyda);
    (*this.*birch_murnaghan_function)(x[i], a, &ymod, dyda);
    sig2i = 1.0 / (sig[i] * sig[i]);
    dy = y[i] - ymod;
    for (j = 0, l = 1; l <= ma; l++) {
      if (ia[l]) {
        wt = dyda[l] * sig2i;
        for (j++, k = 0, m = 1; m <= l; m++)
          if (ia[m]) alpha[j][++k] += wt * dyda[m];
        beta[j] += dy * wt;
      }
    }
    chisq += dy * dy * sig2i;
  }
  for (j = 2; j <= mfit; j++)
    for (k = 1; k < j; k++) alpha[k][j] = alpha[j][k];
}
// ***************************************************
bool aflowFITTING::mrqmin(xvector<double> x,
                          xvector<double> y,
                          xvector<double> &sig,
                          xvector<double> &a,
                          xvector<int> &ia,
                          xmatrix<double> &covar,
                          xmatrix<double> &alpha,
                          double &chisq,
                          void (aflowFITTING::*birch_murnaghan_function)(double,
                                                                         xvector<double>,
                                                                         double *,
                                                                         xvector<double> &),
                          double *alamda)

{
  //int ndata=x.rows;

  if (y.rows != x.rows) {
    cerr << " mrqmin: y.rows!=x.rows" << std::endl;
    return false;
  }
  if (sig.rows != x.rows) {
    cerr << "mrqmin: sig.rows!=x.rows" << std::endl;
    return false;
  }

  int ma = a.rows;
  if (ia.rows != a.rows) {
    cerr << "mrqmin: ia.rows!=a.rows" << std::endl;
    return false;
  }

  int j, k, l;
  static int mfit;
  //double ochisq=0.0;
  static double ochisq;
  static vector<double> *atry;
  static vector<double> *beta;
  static vector<double> *da;
  static vector<vector<double> > *oneda;

  static xmatrix<double> APL_oneda(ma, 1, 1, 1);
  static xvector<double> APL_atry(ma, 1);
  static xvector<double> APL_beta(ma, 1);
  static xvector<double> APL_da(ma, 1);

  if (*alamda < 0.0) {
    atry = new vector<double>(ma + 1, 0.0);
    beta = new vector<double>(ma + 1, 0.0);
    da = new vector<double>(ma + 1, 0.0);
    for (mfit = 0, j = 1; j <= ma; j++)
      if (ia[j]) mfit++;
    oneda = new vector<vector<double> >((mfit + 1), vector<double>(2));
    *alamda = 0.001;
    mrqcof(x, y, sig, a, ia, alpha, APL_beta, chisq, birch_murnaghan_function);
    ochisq = chisq;
    for (j = 1; j <= ma; j++) {
      (*atry)[j] = a[j];
      (*beta)[j] = APL_beta[j];
      APL_atry[j] = (*atry)[j];
    }
  }
  for (j = 1; j <= mfit; j++) {
    for (k = 1; k <= mfit; k++) covar[j][k] = alpha[j][k];
    covar[j][j] = alpha[j][j] * (1.0 + (*alamda));
    (*oneda)[j][1] = (*beta)[j];
    APL_oneda[j][1] = (*oneda)[j][1];
  }
  if (!gaussj(covar, mfit, APL_oneda, 1)) return false;
  for (j = 1; j <= mfit; j++) {
    (*oneda)[j][1] = APL_oneda[j][1];
    (*da)[j] = (*oneda)[j][1];
    APL_da[j] = (*da)[j];
  }
  if (*alamda == 0.0) {
    covsrt(covar, ia, mfit);
    covsrt(alpha, ia, mfit);
    (*atry).clear();
    (*beta).clear();
    (*da).clear();
    for (int u = 0; u <= mfit; u++) (*oneda)[u].clear();
    (*oneda).clear();
    return false;
  }
  for (j = 0, l = 1; l <= ma; l++)
    if (ia[l]) {
      (*atry)[l] = a[l] + (*da)[++j];
      APL_atry[l] = (*atry)[l];
    }
  mrqcof(x, y, sig, APL_atry, ia, covar, APL_da, chisq, birch_murnaghan_function);
  for (int u = 1; u <= ma; u++) {
    (*atry)[u] = APL_atry[u];
    (*da)[u] = APL_da[u];
  }
  if (chisq < ochisq) {
    *alamda *= 0.1;
    ochisq = (chisq);
    for (j = 1; j <= mfit; j++) {
      for (k = 1; k <= mfit; k++) alpha[j][k] = covar[j][k];
      (*beta)[j] = (*da)[j];
    }
    for (l = 1; l <= ma; l++) a[l] = (*atry)[l];
  } else {
    *alamda *= 10.0;
    chisq = ochisq;
  }
  return true;
}
// ***************************************************
}  //apl end
