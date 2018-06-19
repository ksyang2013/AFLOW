// [OBSOLETE] #include <iostream>
// [OBSOLETE] #include <sstream>
// [OBSOLETE] #include <string>
// [OBSOLETE] #include <limits>
// [OBSOLETE] // [OBSOLETE] // [OBSOLETE] #include <cmath>

#include "aflow_apl.h"

#ifdef USE_MKL
#include "mkl_lapack.h"
#endif

using namespace std;

namespace apl {

// ///////////////////////////////////////////////////////////////////////////

// Based on http://www.feynarts.de/diag/
// http://arxiv.org/abs/physics/0607103

void zheevByJacobiRotation(xmatrix<xcomplex<double> >& A,
                           xvector<double>& d,
                           xmatrix<xcomplex<double> >& U) {
  uint MAX_SWEEPS = 50;
  uint TRESHOLD_SWEEP = 4;
  double epsilon = std::numeric_limits<double>::epsilon();

  // This is only for the square matrices...
  if (!A.issquare) {
    throw APLRuntimeError("apl::zheevByJacobiRotation(); The input matrix is not square.");
  }

  // Get dimension
  uint n = A.rows;
  double REDUCTION = 0.04 / (n * n * n * n);

  // Resize the users variables if they are wrong
  if (d.rows != A.rows) {
    xvector<double> new_d(n, 1);
    d = new_d;
  }

  if (!U.issquare || U.rows != A.rows) {
    xmatrix<xcomplex<double> > new_U(n, n, 1, 1);
    U = new_U;
  }

  // Prepare working arrays and variables
  xmatrix<double> ev(2, n, 1, 1);
  double sum, treshold, t, delta, invc, s;
  xcomplex<double> x, y, Apq;

  // Prepare our outputs, U is unit matrix, and d is current A diagonal elements
  U.clear();
  for (uint p = 1; p <= n; p++) {
    ev(1, p) = 0.0;
    ev(2, p) = A(p, p).re;
    d(p) = ev(2, p);
    U(p, p) = 1.0;
  }

  // Loop over sweeps...

  uint nSweep;
  for (nSweep = 1; nSweep <= MAX_SWEEPS; nSweep++) {
    // Convergence criterion, sum of the squares of the off-diagonal elements
    sum = 0.0;
    for (uint q = 2; q <= n; q++) {
      for (uint p = 1; p <= q - 1; p++) {
        sum = sum + real(A(p, q) * conj(A(p, q)));
      }
    }

    if (sum < 0.5 * epsilon) {
      break;
    }

    // Update the threshold for the applying the rotation for matrix element

    if (nSweep < TRESHOLD_SWEEP) {
      treshold = REDUCTION * sum;
    } else {
      treshold = 0.0;
    }

    // Rotate...
    for (uint q = 2; q <= n; q++) {
      for (uint p = 1; p <= q - 1; p++) {
        sum = real(A(p, q) * conj(A(p, q)));

        if ((nSweep > TRESHOLD_SWEEP) &&
            (sum < 0.5 * epsilon * std::max<double>(ev(2, p) * ev(2, p), ev(2, q) * ev(2, q)))) {
          A(p, q) = 0.0;
        } else {
          if (sum > treshold) {
            t = 0.5 * (ev(2, p) - ev(2, q));
            t = 1.0 / (t + copysign(sqrt(t * t + sum), t));

            delta = t * sum;
            ev(1, p) = ev(1, p) + delta;
            ev(2, p) = d(p) + ev(1, p);
            ev(1, q) = ev(1, q) - delta;
            ev(2, q) = d(q) + ev(1, q);

            invc = sqrt(delta * t + 1.0);
            s = t / invc;
            t = delta / (invc + 1.0);

            Apq = A(p, q);

            for (uint j = 1; j <= p - 1; j++) {
              x = A(j, p);
              y = A(j, q);
              A(j, p) = x + s * (conj(Apq) * y - t * x);
              A(j, q) = y - s * (Apq * x + t * y);
            }

            for (uint j = p + 1; j <= q - 1; j++) {
              x = A(p, j);
              y = A(j, q);
              A(p, j) = x + s * (Apq * conj(y) - t * x);
              A(j, q) = y - s * (Apq * conj(x) + t * y);
            }

            for (uint j = q + 1; j <= n; j++) {
              x = A(p, j);
              y = A(q, j);
              A(p, j) = x + s * (Apq * y - t * x);
              A(q, j) = y - s * (conj(Apq) * x + t * y);
            }

            A(p, q) = 0.0;

            for (uint j = 1; j <= n; j++) {
              x = U(p, j);
              y = U(q, j);
              U(p, j) = x + s * (Apq * y - t * x);
              U(q, j) = y - s * (conj(Apq) * x + t * y);
            }
          }
        }
      }
    }

    for (uint p = 1; p <= n; p++) {
      ev(1, p) = 0.0;
      d(p) = ev(2, p);
    }
  }

  if (nSweep > MAX_SWEEPS) {
    throw APLRuntimeError("apl::zheevByJacobiRotation(); Bad convergence.");
  }

  // Order
  for (uint i = 1; i <= n - 1; i++) {
    for (uint j = i + 1; j <= n; j++) {
      if (d(j) < d(i)) {
        double temp = d(j);
        d(j) = d(i);
        d(i) = temp;
        for (uint k = 1; k <= n; k++) {
          xcomplex<double> xtemp = U(j, k);
          U(j, k) = U(i, k);
          U(i, k) = xtemp;
        }
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

#ifdef USE_MKL

// based on http://software.intel.com/sites/products/documentation/hpc/mkl/lapack/mkl_lapack_examples/zheev_ex.c.htm

void zheevMKL(xmatrix<xcomplex<double> >& A, xvector<double>& d, xmatrix<xcomplex<double> >& U) {
  // This is only for the square matrices...
  if (!A.issquare) {
    throw("IPhononCalculator::zheevMKL(); The input matrix is not square.");
  }

  // Get dimension
  int N = A.rows;
  int LDA = A.rows;

  // Resize the users variables if the are wrong
  if (d.rows != A.rows) {
    xvector<double> new_d(N, 1);
    d = new_d;
  }

  if (!U.issquare || U.rows != A.rows) {
    xmatrix<xcomplex<double> > new_U(N, N, 1, 1);
    U = new_U;
  }

  int n = N, lda = LDA, info, lwork;
  MKL_Complex16 wkopt;
  MKL_Complex16* work;
  // Local arrays
  // rwork dimension should be at least max(1,3*n-2)
  double w[N], rwork[3 * N - 2];
  MKL_Complex16 a[LDA * N];

  // Copy matrix from our format to MKL format
  for (int i = A.lrows; i <= A.rows; i++)
    for (int j = A.lcols; j <= A.cols; j++) {
      a[(i - 1) * A.rows + (j - 1)].real = A(i, j).re;
      a[(i - 1) * A.rows + (j - 1)].imag = A(i, j).im;
    }

  // Query and allocate the optimal workspace
  lwork = -1;
  zheev((char*)("Vectors"), (char*)("Lower"), &n, a, &lda, w, &wkopt, &lwork, rwork, &info);
  lwork = (int)wkopt.real;
  work = (MKL_Complex16*)malloc(lwork * sizeof(MKL_Complex16));

  // Solve eigenproblem
  zheev((char*)("Vectors"), (char*)("Lower"), &n, a, &lda, w, work, &lwork, rwork, &info);
  if (info > 0) {
    throw APLRuntimeError("zheevMKL(); The algorithm failed to compute eigenvalues.");
  }

  // Free workspace
  free((void*)work);

  // Copy result to xmatrix and xvector
  for (int i = d.lrows; i <= d.rows; i++) {
    d(i) = w[i - 1];
  }

  for (int i = A.lrows; i <= A.rows; i++)
    for (int j = A.lcols; j <= A.cols; j++) {
      U(i, j).re = a[(i - 1) * U.rows + (j - 1)].real;
      U(i, j).im = a[(i - 1) * U.rows + (j - 1)].imag;
    }
}

#endif

//////////////////////////////////////////////////////////////////////////////

// Will triadiagonalize the input matrix, not optimized, just the first shoot based on the literature equations...

// http://en.wikipedia.org/wiki/QR_decomposition
// http://en.wikipedia.org/wiki/Householder_transformation
// http://books.google.com/books?id=1oDXWLb9qEkC&printsec=frontcover&dq=Introduction+to+Numerical+Analysis&cd=1#v=onepage&q=&f=false (page226)
void tred2(xmatrix<xcomplex<double> >& A) {
  // This is only for the square matrices...
  if (!A.issquare) {
    throw APLRuntimeError("tred2(); The input matrix is not square.");
  }

  // n is the order of the matrix
  int n = (int)A.rows;

  //
  xmatrix<xcomplex<double> > Q(n, n);
  xvector<xcomplex<double> > u(n), x(n), e1(n);

  //cout << "input A" << std::endl;
  //printXMatrix(A);

  for (int i = 1; i <= n - 2; i++) {
    int l = i + 1;

    // Get the column
    for (int j = 1; j < l; j++)
      x[j] = 0;
    for (int j = l; j <= n; j++)
      x[j] = A(i, j);

    // Norm x
    double normx = 0;
    for (int j = 1; j <= n; j++)
      normx += x[j].re * x[j].re + x[j].im * x[j].im;
    normx = sqrt(normx);

    // Is it already zero column?
    if (normx < _AFLOW_APL_EPS_) return;

    //
    xcomplex<double> alpha = aurostd::polar(-normx, aurostd::arg(x[l]));

    // Get e1
    for (int j = 1; j <= n; j++)
      e1[j] = 0;
    e1[l] = 1;

    // Compute vector u = x - alpha * e1
    for (int j = 1; j < l; j++)
      u[j] = 0;
    for (int j = l; j <= n; j++)
      u[j] = x[j] - alpha * e1[j];

    // Normalize u
    double normu = 0;
    for (int j = 1; j <= n; j++)
      normu += u[j].re * u[j].re + u[j].im * u[j].im;
    normu = sqrt(normu);
    for (int j = 1; j <= n; j++)
      u[j] = u[j] / normu;

    // Compute w = x*^T.u / u*^T.x
    xcomplex<double> w = 0;
    xcomplex<double> s1 = 0;
    xcomplex<double> s2 = 0;
    for (int j = 1; j <= n; j++) {
      s1 = s1 + conj(x[j]) * u[j];
      s2 = s2 + conj(u[j]) * x[j];
    }
    w = s1 / s2;

    // Compute Q = I - ( 1 + w ) vv*^T
    for (int j = 1; j <= n; j++)
      for (int k = 1; k <= n; k++) {
        if (j == k)
          Q(j, k) = 1.0;
        else
          Q(j, k) = 0.0;

        Q(j, k) = Q(j, k) - (1.0 + w) * u[j] * conj(u[k]);
      }

    //cout << "Q" << i << "AQ" << i << std::endl;
    //printXMatrix(Q*A*Q);
    A = Q * A * Q;
  }

  //cout << "output A" << std::endl;
  //printXMatrix(A);
}

//////////////////////////////////////////////////////////////////////////////

}  // namespace apl
