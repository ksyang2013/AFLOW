// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *                                                                         *
// ***************************************************************************
// Written by Stefano Curtarolo 1994-2011

#ifndef _AUROSTD_XMATRIX_H_
#define _AUROSTD_XMATRIX_H_

#define _AUROSTD_XMATRIX_DEFAULT_SIZE_ 3
//#define _AUROSTD_XMATRIX_TOLERANCE_IDENTITY_ (1.0e-6)
//#define _AUROSTD_XMATRIX_TOLERANCE_ROUNDOFF_ (1.0e-6)

//double _AUROSTD_XMATRIX_TOLERANCE_IDENTITY_ = 1.0e-6;
//double _AUROSTD_XMATRIX_TOLERANCE_ROUNDOFF_=1.0e-6;

#ifndef _AUROSTD_XMATRIX_TOLERANCE_ROUNDOFF_
#define _AUROSTD_XMATRIX_TOLERANCE_ROUNDOFF_ AUROSTD_ROUNDOFF_TOL //1.0e-6 // DX 171025
#endif

#ifndef _AUROSTD_XMATRIX_TOLERANCE_IDENTITY_
#define _AUROSTD_XMATRIX_TOLERANCE_IDENTITY_ AUROSTD_IDENTITY_TOL //1.0e-6 // DX 171025
#endif

const int _MAX_ITS=100; // maximum iterations in SVDcmp()

// ----------------------------------------------------------------------------
// -------------------------------------------------------------- class xmatrix

namespace aurostd {
  // namespace aurostd
  template<class utype>
    class xmatrix  {
  public:
    //  xmatrix ();                                        // default constructor
    //  xmatrix (int=3);                                   // default constructor
    xmatrix (int=3,int=3,int=1,int=1);                     // default constructor
    //    xmatrix (int,int,int,int);                       // default constructor
    xmatrix (const xmatrix<utype>&);                          // copy constructor
    xmatrix (int,int,utype*);                                 // copy constructor
    xmatrix<utype>& operator=(const xmatrix<utype>&);               // assignment
    ~xmatrix ();                                            // default destructor
    utype* operator[](int) const;                                 // indicize i,j
    xvector<utype> operator()(int) const;                           // indicize i
    utype& operator()(int,int) const;                             // indicize i,j
    utype& operator()(int,int,bool) const;        // indicize boundary conditions
    //      xmatrix operator()(int,int,int,int) const;    // indicize i1,j1,i2,j2
    // math operators
    xmatrix<utype>& operator +=(const xmatrix<utype>&);
    xmatrix<utype>& operator -=(const xmatrix<utype>&);
    xmatrix<utype>& operator *=(const xmatrix<utype>&);
    // xmatrix<utype>& operator /=(const xmatrix<utype>&);
    // std::ostream operator
    //      friend std::ostream& operator<<<utype>(std::ostream&,const xmatrix<utype>&);
    // friend std::ostream& operator<<<utype>(std::ostream&,const xmatrix<utype>&);
    // friend void balanc(const xmatrix<utype>&); // balanced with eigenvalues
    int rows,lrows,urows;                                  // rows dimensions
    int cols,lcols,ucols;                                  // cols dimensions
    bool issquare,isfloat,iscomplex;
    // operations
    void set(const utype&);
    void reset(void);
    void clear(void);
  private:
    //  friend bool operator==(const xmatrix<utype>&,const xmatrix<utype>&);
    utype** corpus;
    //  bool isfloat,iscomplex;
    char size;
    long int msize;
  };
}

// ----------------------------------------------------------------------------
// ------------------------- primitives for template<class utype> xmatrix<utype>

// ------------------------------------------------- unary and binary operators
namespace aurostd {
  template<class utype> xmatrix<utype>
    operator+(const xmatrix<utype>&,const xmatrix<utype>&) __xprototype;
  
  template<class utype> xmatrix<utype>
    operator-(const xmatrix<utype>&,const xmatrix<utype>&) __xprototype;
  
  template<class utype> xmatrix<utype>                                   // matrix*matrix
    operator*(const xmatrix<utype>&,const xmatrix<utype>&) __xprototype; // matrix*matrix
  
  template<class utype> xvector<utype>                                   // matrix*vector
    operator*(const xmatrix<utype>&,const xvector<utype>&) __xprototype; // matrix*vector
  
  template<class utype> xvector<utype>                                   // vector_row*matrix
    operator*(const xvector<utype>&,const xmatrix<utype>&) __xprototype; // vector_row*matrix

  template<class utype> xmatrix<utype>                                   // scalar*matrix
    operator*(const utype,const xmatrix<utype>&) __xprototype;           // scalar*matrix
  
  template<class utype> xmatrix<utype>                                   // matrix*scalar
    operator*(const xmatrix<utype>&,const utype) __xprototype;           // matrix*scalar
  
  template<class utype> xmatrix<utype>
    operator/(const xmatrix<utype>&,const utype) __xprototype;
  
  template<class utype>
    std::ostream& operator<<(std::ostream&,const xmatrix<utype>&) __xprototype;
  
  template<class utype> xmatrix<utype>
    operator+(const xmatrix<utype>&) __xprototype;
  
  template<class utype> xmatrix<utype>
    operator-(const xmatrix<utype>&) __xprototype;
  
  template<class utype> bool
    identical(const xmatrix<utype>&,const xmatrix<utype>&,const utype&,const char& _mode_) __xprototype;
  
  template<class utype> bool
    identical(const xmatrix<utype>&,const xmatrix<utype>&,const utype&) __xprototype;
  
  template<class utype> bool
    rel_identical(const xmatrix<utype>&,const xmatrix<utype>&,const utype&) __xprototype;
  
  template<class utype> bool
    abs_identical(const xmatrix<utype>&,const xmatrix<utype>&,const utype&) __xprototype;
  
  template<class utype> bool
    identical(const xmatrix<utype>&,const xmatrix<utype>&) __xprototype;
  
  template<class utype> bool
    operator==(const xmatrix<utype>&,const xmatrix<utype>&) __xprototype;
  
  template<class utype> bool
    isdifferent(const xmatrix<utype>&,const xmatrix<utype>&,const utype&) __xprototype;
  
  template<class utype> bool
    isdifferent(const xmatrix<utype>&,const xmatrix<utype>&) __xprototype;
  
  template<class utype> bool
    isequal(const xmatrix<utype>&,const xmatrix<utype>&,const utype&) __xprototype;
  
  template<class utype> bool
    isequal(const xmatrix<utype>&,const xmatrix<utype>&) __xprototype;
  
  template<class utype> bool
    operator!=(const xmatrix<utype>&,const xmatrix<utype>&) __xprototype;

  //CO - START
  template<class utype> bool
    isidentity(const xmatrix<utype>&) __xprototype;
  //CO - END
  
  template<class utype> bool
    isdiagonal(const xmatrix<utype>&,const utype& _eps_=(utype)_AUROSTD_XMATRIX_TOLERANCE_IDENTITY_) __xprototype; // DX 171025

  template<class utype> bool
    isinteger(const xmatrix<utype>&,const utype& tol=(utype)0.01) __xprototype;

  template<class utype> bool
    issymmetric(const xmatrix<utype>&) __xprototype;

  template<class utype> bool
    isantisymmetric(const xmatrix<utype>&) __xprototype;

  template<class utype> bool
    ishermitian(const xmatrix<utype>&) __xprototype;

  template<class utype> bool
    isantihermitian(const xmatrix<utype>&) __xprototype;

}

// ----------------------------------------------------------- xmatrix constrtuction
namespace aurostd {
  template<class utype> xmatrix<utype>
    reshape(const xvector<utype>&) __xprototype;
  
  template<class utype> xmatrix<utype>
    reshape(const xvector<utype>&,const xvector<utype>&) __xprototype;
  
  template<class utype> xmatrix<utype>
    reshape(const xvector<utype>&,const xvector<utype>&,const xvector<utype>&) __xprototype;
  
  template<class utype> xmatrix<utype>
    reshape(const xvector<utype>&,const xvector<utype>&,const xvector<utype>&,const xvector<utype>&) __xprototype;

  template<class utype> xmatrix<utype>
    reshape(const xvector<utype>&,const xvector<utype>&,const xvector<utype>&,const xvector<utype>&,const xvector<utype>&) __xprototype;

  template<class utype> xmatrix<utype>
    reshape(const xvector<utype>&,const xvector<utype>&,const xvector<utype>&,const xvector<utype>&,const xvector<utype>&,const xvector<utype>&) __xprototype;

  template<class utype> xmatrix<utype>
    reshape_cols(const xvector<utype>&) __xprototype;
  
  template<class utype> xmatrix<utype>
    reshape_cols(const xvector<utype>&,const xvector<utype>&) __xprototype;
  
  template<class utype> xmatrix<utype>
    reshape_cols(const xvector<utype>&,const xvector<utype>&,const xvector<utype>&) __xprototype;
  
  template<class utype> xmatrix<utype>
    reshape_cols(const xvector<utype>&,const xvector<utype>&,const xvector<utype>&,const xvector<utype>&) __xprototype;

  template<class utype> xmatrix<utype>
    reshape_cols(const xvector<utype>&,const xvector<utype>&,const xvector<utype>&,const xvector<utype>&,const xvector<utype>&) __xprototype;

  template<class utype> xmatrix<utype>
    reshape_cols(const xvector<utype>&,const xvector<utype>&,const xvector<utype>&,const xvector<utype>&,const xvector<utype>&,const xvector<utype>&) __xprototype;

  template<class utype> xmatrix<utype>
    reshape_rows(const xvector<utype>&) __xprototype;
  
  template<class utype> xmatrix<utype>
    reshape_rows(const xvector<utype>&,const xvector<utype>&) __xprototype;
  
  template<class utype> xmatrix<utype>
    reshape_rows(const xvector<utype>&,const xvector<utype>&,const xvector<utype>&) __xprototype;
  
  template<class utype> xmatrix<utype>
    reshape_rows(const xvector<utype>&,const xvector<utype>&,const xvector<utype>&,const xvector<utype>&) __xprototype;

  template<class utype> xmatrix<utype>
    reshape_rows(const xvector<utype>&,const xvector<utype>&,const xvector<utype>&,const xvector<utype>&,const xvector<utype>&) __xprototype;

  template<class utype> xmatrix<utype>
    reshape_rows(const xvector<utype>&,const xvector<utype>&,const xvector<utype>&,const xvector<utype>&,const xvector<utype>&,const xvector<utype>&) __xprototype;
}

// ----------------------------------------------------------- xmatrix example types
namespace aurostd {
  //171008 - CO
  template<class utype> xmatrix<utype>
    eyes(int=3,int=3,int=1,int=1) __xprototype;
  template<class utype> xmatrix<utype>
    ones(int=3,int=3,int=1,int=1) __xprototype;
}

// ----------------------------------------------------------- xmatrix cast
namespace aurostd {
  template<class utype> xmatrix<long double>
    xlongdouble(const xmatrix<utype>&) __xprototype;
  
  template<class utype> xmatrix<double>
    xdouble(const xmatrix<utype>&) __xprototype;
  
  template<class utype> xmatrix<float>
    xfloat(const xmatrix<utype>&) __xprototype;
  
  template<class utype> xmatrix<long int>
    xlongint(const xmatrix<utype>&) __xprototype;

  template<class utype> xmatrix<int>
    xint(const xmatrix<utype>&) __xprototype;

  template<class utype> xmatrix<char>
    xchar(const xmatrix<utype>&) __xprototype;

  template<class utype> vector<vector<utype> >
    xmatrix2vectorvector(const xmatrix<utype>& xmat) __xprototype;

  template<class utype> xmatrix<utype>
    vectorvector2xmatrix(const vector<vector<utype> >& xmat) __xprototype;
}  

// ----------------------------------------------------------- xmatrix functions
namespace aurostd {
  template<class utype> void
    reset(xmatrix<utype>&) __xprototype;
  
  template<class utype> void
    clear(xmatrix<utype>&) __xprototype;
  
  template<class utype> void
    set(xmatrix<utype>&,const utype&) __xprototype;
  
  //  template<class utype>                        
  //  xvector<utype> vector(const xmatrix<utype>&) __xprototype;
  
  template<class utype> utype
    det(const xmatrix<utype>&) __xprototype;
  
  template<class utype> utype
    determinant(const xmatrix<utype>&) __xprototype;
  
  template<class utype> xmatrix<utype>
    submatrix(const xmatrix<utype>&,const int&,const int&) __xprototype;
  
  template<class utype> utype
    minordet(const xmatrix<utype>&,const int&,const int&) __xprototype;
  
  template<class utype> utype
    minordeterminant(const xmatrix<utype>&,const int&,const int&) __xprototype;
  
  template<class utype> xmatrix<utype>
    inverse(const xmatrix<utype>&) __xprototype;

  template<class utype> bool
    inverse(const xmatrix<utype>&, xmatrix<utype>&) __xprototype; // RETURN ERROR if non invertible

  template<class utype> xmatrix<utype>                  // clear values too small
    roundoff(const xmatrix<utype>&,utype) __xprototype; // claar values too small

  template<class utype> xmatrix<utype>                  // clear values too small
    roundoff(const xmatrix<utype>&) __xprototype;       // clear values too small

  template<class utype> utype
    sum(const xmatrix<utype>&) __xprototype;
  
  template<class utype> xvector<utype>
    sum_column(const xmatrix<utype>&) __xprototype;
  
  template<class utype> xvector<utype>
    mean_column(const xmatrix<utype>&) __xprototype;
  
  template<class utype>
    xvector<utype> sum_row(const xmatrix<utype>&) __xprototype;
  
  template<class utype>
    xvector<utype> mean_row(const xmatrix<utype>&) __xprototype;
  
  template<class utype> utype
    min(const xmatrix<utype>&) __xprototype;
  
  template<class utype> utype
    min(const xmatrix<utype>&,int&,int&) __xprototype;
  
  template<class utype> utype
    max(const xmatrix<utype>&) __xprototype;
  
  template<class utype> utype
    max(const xmatrix<utype>&,int&,int&) __xprototype;
  
  // DX 1/15/18 [OBSOLETE] template<class utype> double
  // DX 1/15/18 [OBSOLETE]  trace(const xmatrix<utype>&) __xprototype;
  template<class utype> utype                     // DX 1/15/18 - double to utype
    trace(const xmatrix<utype>&) __xprototype;    // DX 1/15/18 - double to utype
  
  template<class utype> xmatrix<utype>
    shift_left(const xmatrix<utype>&) __xprototype;
  
  template<class utype> xmatrix<utype>
    shift_right(const xmatrix<utype>&) __xprototype;
  
  template<class utype> xmatrix<utype>
    shift_up(const xmatrix<utype>&) __xprototype;
  
  template<class utype> xmatrix<utype>
    shift_down(const xmatrix<utype>&) __xprototype;
  
  template<class utype> void                                        // swap_cols
    swap_cols(xmatrix<utype>&,const int&,const int&) __xprototype;  // swap_cols

  template<class utype> void                                           // swap_columns
    swap_columns(xmatrix<utype>&,const int&,const int&) __xprototype;  // swap_columns

  template<class utype> void                                           // swap_rows
    swap_rows(xmatrix<utype>&,const int&,const int&) __xprototype;    // swap_rows

  template<class utype> xmatrix<utype>
    identity(const xmatrix<utype>&) __xprototype;
    
  // template<class utype> xmatrix<utype>
  //  identity(const int&,const int&,utype) __xprototype;

  template<class utype> xmatrix<utype>
    identity(const utype&,const int&,const int&) __xprototype;

  template<class utype> xmatrix<utype>  // traspose
    trasp(const xmatrix<utype>&) __xprototype; // 25 january 2000
  
  template<class utype> xmatrix<utype>  // traspose
    trasp(const xvector<utype>&) __xprototype; // 5 febrary 2000

  // Mathematical operations
  
  template<class utype> xmatrix<utype>       // 10 Aug 2007
    abs(const xmatrix<utype>&) __xprototype;

  template<class utype> xmatrix<utype>       // 10 Oct 2007
    mabs(const xmatrix<utype>&) __xprototype;

  template<class utype> xmatrix<utype>       // Sept 2008
    sign(const xmatrix<utype>&) __xprototype;

  template<class utype> xmatrix<utype>       // Sept 2008
    nint(const xmatrix<utype>&) __xprototype;

  template<class utype> xmatrix<utype>       // Oct 2014
    floor(const xmatrix<utype>&) __xprototype;

  template<class utype> xmatrix<utype>       // Oct 2014
    trunc(const xmatrix<utype>&) __xprototype;

  template<class utype> xmatrix<utype>       // Oct 2014
    round(const xmatrix<utype>&) __xprototype;

  template<class utype> xmatrix<utype>       // Oct 2014
    ceil(const xmatrix<utype>&) __xprototype;

  template<class utype> xmatrix<utype>
    exp(const xmatrix<utype>&) __xprototype; // 31 july 1999, fix 10 Aug 2007
  
  template<class utype> xmatrix<utype>
    log(const xmatrix<utype>&) __xprototype;// MISSING, probably never
  
  template<class utype> xmatrix<utype>
    log10(const xmatrix<utype>&) __xprototype;   // MISSING
  
  template<class utype> xmatrix<utype>
    sin(const xmatrix<utype>&) __xprototype; // 10 august 1999
  
  template<class utype> xmatrix<utype>
    asin(const xmatrix<utype>&) __xprototype;   // MISSING
  
  template<class utype> xmatrix<utype>
    cos(const xmatrix<utype>&) __xprototype; // 10 august 1999
  
  template<class utype> xmatrix<utype>
    sec(const xmatrix<utype>&) __xprototype;   // MISSING
  
  template<class utype> xmatrix<utype>
    cosec(const xmatrix<utype>&) __xprototype;   // MISSING
  
  template<class utype> xmatrix<utype>
    asec(const xmatrix<utype>&) __xprototype;   // MISSING
  
  template<class utype> xmatrix<utype>
    acosec(const xmatrix<utype>&) __xprototype;   // MISSING
  
  template<class utype> xmatrix<utype>
    acos(const xmatrix<utype>&) __xprototype;   // MISSING
  
  template<class utype> xmatrix<utype>
    tan(const xmatrix<utype>&) __xprototype;   // MISSING
  
  template<class utype> xmatrix<utype>
    atan(const xmatrix<utype>&) __xprototype;   // MISSING
  
  template<class utype> xmatrix<utype>
    sinh(const xmatrix<utype>&) __xprototype;  // 12 august 1999
  
  template<class utype> xmatrix<utype>
    asinh(const xmatrix<utype>&) __xprototype;   // MISSING
  
  template<class utype> xmatrix<utype>
    cosh(const xmatrix<utype>&) __xprototype;   // 12 august 1999
  
  template<class utype> xmatrix<utype>
    acosh(const xmatrix<utype>&) __xprototype;   // MISSING
  
  template<class utype> xmatrix<utype>
    tanh(const xmatrix<utype>&) __xprototype;   // MISSING
  
  template<class utype> xmatrix<utype>
    atanh(const xmatrix<utype>&) __xprototype;   // MISSING
  
  template<class utype>                          // 25 april 2007
    void GaussJordan(xmatrix<utype>&, xmatrix<utype>&) __xprototype;

  // least square stuff aurostd adaptation of nrecipes
  template<class utype>                          // 1 August 2014
    void gaussj(xmatrix<utype>& a, int n, xmatrix<utype>& b, int m);

  template<class utype>                          // 1 August 2014
    void covsrt(xmatrix<utype>& covar, xvector<int> ia, int mfit);

  template<class utype>                          // 1 August 2014
    void lfit(xvector<utype> x, xvector<utype> y, xvector<utype> sig, 
              xvector<utype>& a, xvector<int> ia, 
              xmatrix<utype>& covar, utype& chisq, 
              void (*funcs)(utype, xvector<utype>&));
}

namespace aurostd {
  // ORTHOGONALITY
  template<class utype> utype
    orthogonality_defect(const xmatrix<utype>& basis) __xprototype;   // imported May 2018
  template<class utype> bool
    gaussian_reduce_two_vectors(xvector<utype>& B, xvector<utype>& C,utype eps) __xprototype;   // imported May 2018
  template<class utype> void
    reduce_A_in_ABC(xvector<utype>& A, xvector<utype>& B, xvector<utype>& C,utype eps) __xprototype;   // imported May 2018
  template<class utype> utype
    reduce_to_shortest_basis(const xmatrix<utype>& IN,xmatrix<utype>& OUT,utype eps,bool VERBOSE) __xprototype;   // imported May 2018
  template<class utype> xmatrix<utype>
    reduce_to_shortest_basis(const xmatrix<utype>& IN,utype eps,bool VERBOSE) __xprototype;   // imported May 2018
  template<class utype> xmatrix<utype>
    reduce_to_shortest_basis(const xmatrix<utype>& IN) __xprototype;   // imported May 2018
}

namespace aurostd {
  // EIGENPROBLEMS
  template<class utype>
    int jacobi(const xmatrix<utype> &ain,xvector<utype> &d,xmatrix<utype> &v) __xprototype;
  // Computes all eigenvalues and eigenvectors of a real symmetric xmatrix a[1..n][1..n].
  // On output, elements of a above the diagonal are destroyed. d[1..n] returns the eigenvalues of a.
  // v[1..n][1..n] is a matrix whose columns contain, on output, the normalized eigenvectors of
  // a. The function returns the number of Jacobi rotations that were required.
  template<class utype>
    void eigsrt(xvector<utype> &d,xmatrix<utype> &v) __xprototype;
  // Given the eigenvalues d[1..n] and eigenvectors v[1..n][1..n] as output from jacobi
  // or tqli,this routine sorts the eigenvalues into descending order, and rearranges
  // the columns of v correspondingly. The method is straight insertion and is N2 rather than NlogN;
  // but since you have just done an N3 procedure to get the eigenvalues, you can afford yourself
  // this little indulgence.
  template<class utype> 
    xmatrix<utype> generalHouseHolderQRDecomposition(xmatrix<utype>& mat,const utype& tol=_AUROSTD_XMATRIX_TOLERANCE_IDENTITY_);
  // general Householder, A is mxn, m>=n
  // output:  Q
  // mat will change to R
  // See Numerical Linear Algebra, Trefethen and Bau, pg. 73
  template<class utype>
    void tred2(const xmatrix<utype> &a,xvector<utype> &d,xvector<utype> &e) __xprototype;
  // Householder reduction of a real, symmetric matrix a[1..n][1..n].
  // On output, a is replaced by the orthogonal matrix Q eﬀecting the
  // transformation. d[1..n] returns the diagonal elments of
  // the tridiagonal matrix, and e[1..n] the oﬀ-diagonal elements, with e[1]=0.
  template<class utype>
    void tqli(xvector<utype> &d,xvector<utype> &e,xmatrix<utype> &z) __xprototype;
  // QL algorithm with implicit shifts, to determine the eigenvalues
  // and eigenvectors of a real, symmetric, tridiagonal matrix, or of a real,
  // symmetric matrix previously reduced by tred2
  // On input, d[1..n] contains the diagonal elements of the tridiagonal
  // matrix. On output, it returns the eigenvalues. The vectore[1..n]
  // inputs the subdiagonal elements of the tridiagonal matrix, with e[1] arbitrary.
  // On output e is destroyed. When finding only the eigenvalues, several lines
  // maybe omitted, as noted in the comments. If the eigenvectors of a tridiagonal
  // matrix are desired, the matrix z[1..n][1..n] is input as the identity
  // matrix. If the eigenvectors of a matrix that has been reduced by tred2
  // are required, then z is input as the matrix output by tred2.
  // In either case, the kth column of z returns the normalized eigenvector
  // corresponding to d[k].
  template<class utype>
    void balanc(xmatrix<utype> &a) __xprototype;
  // Given a matrix a[1..n][1..n], this routine replaces it by
  // a balanced matrix with i dentical eigenvalues. A symmetric matrix
  // is already balanced and is unaﬀected by this procedure. The
  // parameter RADIX should be the machine’s ﬂoating-point radix.
  template<class utype>
    void elmhes(xmatrix<utype> &a) __xprototype;
  // Reduction to Hessenberg form by the elimination method.
  // The real, nonsymmetric matrix a[1..n][1..n] is replaced by an upper
  // Hessenberg matrix with identical eigenvalues.
  // Recommended, but not required, is that this routine be preceded
  // by balanc. On output, the Hessenberg matrix is in elements a[i][j] with i<=j+1.
  // Elements with i>j+1 are to be thought of as zero,
  // but are returned with random values.
  template<class utype>
    void hqr(xmatrix<utype> &a,xvector<utype> &wr,xvector<utype> &wi) __xprototype;
  // Finds all eigenvalues of an upper Hessenberg matrix a[1..n][1..n].
  // On input a can be exactly as output from elmhes; on output it is destroyed.
  // The real and imaginary parts of the eigenvalues are returned in
  //wr[1..n] and wi[1..n], respectively.
  template<class utype>
    void eigen(const xmatrix<utype> &ain,xvector<utype> &wr,xvector<utype> &wi) __xprototype;
  // Finds all eigenvalues of matrix a[1..n][1..n]. The real and imaginary parts
  // of the eigenvalues are returned in wr[1..n] and wi[1..n], respectively.
}


// -------------------------------------------------------------- class cematrix
namespace aurostd {
  // namespace aurostd
  class cematrix {
    // add the least square and singular decomposition methods to xmatrix class
  private:
    int nrow, ncol; // number of rows and columes of matrix A
    xmatrix<double> M; // original matrix
    xvector<double> W; // A = U Z V^T
    xmatrix<double > U, V;
    xvector<double> a_vec; // a vector
    vector<double> a_nvec; // a vector (same as a_vec but in c++ vector form)
    double chisq; // value of the norm square
    xmatrix<double> Cov; // covariance matrix    
  public:
    cematrix();
    cematrix(const xmatrix<double> & A_in);
    ~cematrix();
    void LeastSquare(xvector<double> & y_vec, xvector<double> & y_sigma);
    //void SVDcmp_NR(); // Singular Value Decomposition
    bool SVDcmp_NR(); // Singular Value Decomposition
    void SVDsolve(xvector<double>& b_vec); // Solve vector a
    void SVDFit(xvector<double>& x_sigma, xvector<double>& y); // Least Square fitting using SVD
    void SVDvar(); // get covariance matrix
    xvector<double> GetFitVector() { return a_vec;}; // output a_vec
    xmatrix<double> InverseMatrix(); // inverse of a matrix, not necessary square
    xvector<double> EigenValues(); // calculate eigenvalues of a matrix, not necessary square
    double Pythag2(double a, double b); // (a^2+b^2)^(1/2) avoiding overflow and underflow
    double _sign(double a, double b) { return b > -1? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);}
    vector<double> AVec() { return a_nvec;};
    double ChiSQ() { return chisq;};
  };
}

// -------------------------------------------------------------- CMdet
namespace aurostd {
  template<class utype> utype
    CMdet(const xmatrix<utype>& B);
}

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
#endif

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *                                                                         *
// ***************************************************************************

