// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *                                                                         *
// ***************************************************************************
// Written by Stefano Curtarolo 1994-2011

#ifndef _AUROSTD_XVECTOR_H_
#define _AUROSTD_XVECTOR_H_
#define _AUROSTD_XVECTOR_DEFAULT_SIZE_ 3


#define BOUNDARY_CONDITIONS_NONE 0
#define BOUNDARY_CONDITIONS_PERIODIC 1    

#define _AUROSTD_XVECTOR_TOLERANCE_IDENTITY_ AUROSTD_IDENTITY_TOL //1.0e-6
#define _AUROSTD_XVECTOR_TOLERANCE_ROUNDOFF_ AUROSTD_ROUNDOFF_TOL //1.0e-6

// -------------------------------------------------------------- class xvector

namespace aurostd {
  // namespace aurostd
  template<class utype>
    class xvector {
  public:
    //   xvector();                                        // default constructor
    // xvector(int);                                       // default constructor
    xvector(int=3,int=1);                                  // default constructor
    xvector(const xvector<utype>&);                        // copy constructor
    // xvector (const xmatrix<utype>&);                    // make a vector of a xmatrix
    xvector<utype>& operator=(const xvector<utype>&);	   // assignment
    //     operator xvector<utype>() { return *this;};      // IBM_CPP
    ~xvector();                                            // default destructor
    utype& operator[](int) const;		                    // indicize
    utype& operator()(int) const;		                    // indicize
    utype& operator()(int,bool) const;        // indicize boundary conditions
    xvector<utype>& operator +=(const xvector<utype>&);
    xvector<utype>& operator +=(utype); //CO 180409
    xvector<utype>& operator -=(const xvector<utype>&);
    xvector<utype>& operator -=(utype); //CO 180409
    xvector<utype>& operator *=(utype r); //(const utype& r) //CO 171130 - v*=v[1] doesn't work //CO 180409
    xvector<utype>& operator /=(utype r); //(const utype& r) //CO 171130 - v/=v[1] doesn't work //CO 180409
    // xvector<utype>& operator *=(const xvector<utype>&);
    // xvector<utype>& operator /=(const xvector<utype>&);
    //    friend std::ostream& operator<<<utype>(std::ostream&,const xvector<utype>&);
    //    friend std::ostream& operator< <utype>(std::ostream&,const xvector<utype>&);
    int rows,lrows,urows;  
    bool isfloat,iscomplex;
    // operations
    void set(const utype&);
    void reset(void);
    void clear(void);
  private:
    utype *corpus;
    // bool isfloat,iscomplex;
    char size;
    long int vsize;
  };
}

// ----------------------------------------------------------------------------
// ------------------------- primitives for template<class utype> xvector<utype>

// ----------------------------------------------------------------------------
// ------------------------------------------------- unary and binary operators

namespace aurostd {
  // namespace aurostd
  // template<class utype> std::ostream& operator<<(std::ostream&,const xvector<utype>&);
  
  template<class utype>
    std::ostream& operator<<(std::ostream&,const xvector<utype>&) __xprototype;
  
  template<class utype> xvector<utype>
    operator+(const xvector<utype>&) __xprototype;

  template<class utype> xvector<utype>
    operator-(const xvector<utype>&) __xprototype;

  template<class utype> xvector<utype>
    operator+(const xvector<utype>&,const xvector<utype>&) __xprototype;

  template<class utype> xvector<utype>
    operator-(const xvector<utype>&,const xvector<utype>&) __xprototype;

  template<class utype> xvector<utype>
    operator%(const xvector<utype>&,const xvector<utype>&) __xprototype;

  template<class utype> xvector<utype>
    operator+(const utype,const xvector<utype>&) __xprototype;

  template<class utype> xvector<utype>
    operator+(const xvector<utype>&,const utype) __xprototype;

  template<class utype,class stype> xvector<utype>
    operator+(const stype,const xvector<utype>&) __xprototype;

  template<class utype,class stype> xvector<utype>
    operator+(const xvector<utype>&,const stype) __xprototype;

  template<class utype> xvector<utype>
    operator-(const utype,const xvector<utype>&) __xprototype;

  template<class utype> xvector<utype>
    operator-(const xvector<utype>&,const utype) __xprototype;

  template<class utype,class stype> xvector<utype>
    operator-(const stype,const xvector<utype>&) __xprototype;

  template<class utype,class stype> xvector<utype>
    operator-(const xvector<utype>&,const stype) __xprototype;

  template<class utype> xvector<utype>
    operator*(const utype,const xvector<utype>&) __xprototype;

  template<class utype> xvector<utype>
    operator*(const xvector<utype>&,const utype) __xprototype;

  template<class utype,class stype> xvector<utype>
    operator*(const stype,const xvector<utype>&) __xprototype;

  template<class utype,class stype> xvector<utype>
    operator*(const xvector<utype>&,const stype) __xprototype;

  template<class utype> xvector<utype>
    operator/(const utype,const xvector<utype>&) __xprototype;

  template<class utype> xvector<utype>
    operator/(const xvector<utype>&,const utype) __xprototype;

  template<class utype,class stype> xvector<utype>
    operator/(const stype,const xvector<utype>&) __xprototype;

  template<class utype,class stype> xvector<utype>
    operator/(const xvector<utype>&,const stype) __xprototype;

  template<class utype> xvector<utype>
    operator<<(const xvector<utype>&,const xvector<utype>&) __xprototype;

  template<class utype> xvector<utype>
    operator<<(const xvector<utype>&,const utype) __xprototype;

  template<class utype> xvector<utype>
    operator<<(const utype,const xvector<utype>&) __xprototype;

  template<class utype> utype                        
    operator*(const xvector<utype>&,const xvector<utype>&) __xprototype;

  template<class utype> utype                        
    scalar_product(const xvector<utype>&,const xvector<utype>&) __xprototype;

  template<class utype> xvector<utype>                
    vector_product(const xvector<utype>&,const xvector<utype>&) __xprototype;

  template<class utype> xvector<char>                    // is xvector > scalar ?
    operator>(const xvector<utype>&,const utype&) __xprototype;

  template<class utype> xvector<char>                    // is xvector < scalar ?
    operator<(const xvector<utype>&,const utype&) __xprototype;

  template<class utype> xvector<char>                    // is xvector == scalar ?
    operator==(const xvector<utype>&,const utype&) __xprototype;

  template<class utype> xvector<char>                    // is xvector > xvector ?
    operator>(const xvector<utype>&,const xvector<utype>&) __xprototype;

  template<class utype> xvector<char>                    // is xvector < xvector ?
    operator<(const xvector<utype>&,const xvector<utype>&) __xprototype;

  template<class utype> bool
    identical(const xvector<utype>&,const xvector<utype>&,const utype&) __xprototype;

  template<class utype> bool
    identical(const xvector<utype>&,const xvector<utype>&) __xprototype;

  template<class utype> bool
    operator==(const xvector<utype>&,const xvector<utype>&) __xprototype;

  template<class utype> bool
    isdifferent(const xvector<utype>&,const xvector<utype>&,const utype&) __xprototype;

  template<class utype> bool
    isdifferent(const xvector<utype>&,const xvector<utype>&) __xprototype;

  template<class utype> bool
    isequal(const xvector<utype>&,const xvector<utype>&,const utype&) __xprototype;

  template<class utype> bool
    isequal(const xvector<utype>&,const xvector<utype>&) __xprototype;

  template<class utype> bool
    operator!=(const xvector<utype>&,const xvector<utype>&) __xprototype;
  
  template<class utype> bool
    isinteger(const xvector<utype>&,const utype& tol=(utype)0.01) __xprototype; //CO 180409
  
  template<class utype> bool
    iszero(const xvector<utype>&, const double& tol=1e-7) __xprototype;  // ME180702
  // CONSTRUCTIONS OF VECTORS FROM SCALARS
  
  template<class utype> xvector<utype>
    reshape(const utype&) __xprototype;
  
  template<class utype> xvector<utype>
    reshape(const utype&,const utype&) __xprototype;
  
  template<class utype> xvector<utype>
    reshape(const utype&,const utype&,const utype&) __xprototype;
  
  template<class utype> xvector<utype>
    reshape(const utype&,const utype&,const utype&,const utype&) __xprototype;

  template<class utype> xvector<utype>
    reshape(const utype&,const utype&,const utype&,const utype&,const utype&) __xprototype;

  template<class utype> xvector<utype>
    reshape(const utype&,const utype&,const utype&,const utype&,const utype&,const utype&) __xprototype;


  // CAST ON XVECTORS

  template<class utype> xvector<long double>
    xlongdouble(const xvector<utype>&) __xprototype;
  
  template<class utype> xvector<double>
    xdouble(const xvector<utype>&) __xprototype;
  
  template<class utype> xvector<double>
    floor(const xvector<utype>&) __xprototype;
  
  template<class utype> xvector<double>
    ceil(const xvector<utype>&) __xprototype;
  
  template<class utype> xvector<double>
    round(const xvector<utype>&) __xprototype;
  
  template<class utype> xvector<double>
    trunc(const xvector<utype>&) __xprototype;
  
  template<class utype> xvector<float>
    xfloat(const xvector<utype>&) __xprototype;
  
  template<class utype> xvector<long int>
    xlongint(const xvector<utype>&) __xprototype;

  template<class utype> xvector<int>
    xint(const xvector<utype>&) __xprototype;

  template<class utype> xvector<char>
    xchar(const xvector<utype>&) __xprototype;

  template<class utype> vector<utype>
    xvector2vector(const xvector<utype>&) __xprototype;

  template<class utype> xvector<utype>
    vector2xvector(const vector<utype>&,int lrows=1) __xprototype; //CO 180409

  // OPERATIONS ON XVECTORS

  template<class utype> void
    set(xvector<utype>&, const utype&) __xprototype;

  template<class utype> void
    reset(xvector<utype>&) __xprototype;

  template<class utype> void
    clear(xvector<utype>&) __xprototype;

  template<class utype> xvector<utype>
    vabs(const xvector<utype>&) __xprototype;

  template<class utype> xvector<utype>
    abs(const xvector<utype>&) __xprototype;

  template<class utype> xvector<utype>
    sign(const xvector<utype>&) __xprototype;

  template<class utype> xvector<utype>
    nint(const xvector<utype>&) __xprototype;

  template<class utype> utype
    modulus(const xvector<utype>&) __xprototype;

  template<class utype> utype
    modulussquare(const xvector<utype>&) __xprototype;

  template<class utype> utype
    modulus2(const xvector<utype>&) __xprototype;

  template<class utype> utype
    sum(const xvector<utype>&) __xprototype;

  template<class utype> utype
    min(const xvector<utype>&) __xprototype;

  template<class utype> utype
    min(const xvector<utype>&,int&) __xprototype;

  template<class utype> int
    mini(const xvector<utype>&) __xprototype;

  template<class utype> utype
    max(const xvector<utype>&) __xprototype;

  template<class utype> utype
    max(const xvector<utype>&,int&) __xprototype;

  template<class utype> int
    maxi(const xvector<utype>&) __xprototype;

  template<class utype> xvector<utype>                  // clear values too small
    roundoff(const xvector<utype>&,utype tol=(utype)_AUROSTD_XVECTOR_TOLERANCE_ROUNDOFF_) __xprototype; // claar values too small //CO 180409

  int gcd(const xvector<int>&);                         // get GCD of vector //CO 180409

  template<class utype> xvector<utype>                                       // simply divide by GCD (useful for compounds) //CO 180409
    reduceByGCD(const xvector<utype>& in_V,                                  // simply divide by GCD (useful for compounds) //CO 180409
        const utype& tol=_AUROSTD_XVECTOR_TOLERANCE_IDENTITY_) __xprototype; // simply divide by GCD (useful for compounds) //CO 180409

  template<class utype> xvector<utype> 
    normalizeSumToOne(const xvector<utype>& in_V,
        const utype& tol=(utype)_AUROSTD_XVECTOR_TOLERANCE_ROUNDOFF_) __xprototype; //CO 180801
  
  template<class utype> void                                   // swap
    swap(xvector<utype>&,const int&,const int&) __xprototype;  // swap
  
  template<class utype> void                              //shift lrows so first index is i //CO 180409
    shiftlrows(xvector<utype>&,const int&) __xprototype;  //shift lrows so first index is i //CO 180409

  // Complex operations
  template<class utype> xvector<utype>
    conj(const xvector<utype>&) __xprototype;

  // EXPONENTIAL OPERATIONS

  template<class utype> xvector<utype>
    exp(const xvector<utype>&) __xprototype;

  template<class utype> xvector<utype>
    log(const xvector<utype>&) __xprototype;

  template<class utype> xvector<utype>
    log10(const xvector<utype>&) __xprototype;

  // TRIDIMENSIONAL OPERATIONS
  
  template<class utype> utype
    distance(const xvector<utype>&,const xvector<utype>&) __xprototype;
  
  // TRIGONOMETRIC OPERATIONS
  
  template<class utype> xvector<utype>
    sin(const xvector<utype>&) __xprototype;

  template<class utype> xvector<utype>
    cos(const xvector<utype>&) __xprototype;

  template<class utype> xvector<utype>
    asin(const xvector<utype>&) __xprototype;

  template<class utype> xvector<utype>
    acos(const xvector<utype>&) __xprototype;

  template<class utype> xvector<utype>
    tan(const xvector<utype>&) __xprototype;

  template<class utype> xvector<utype>
    atan(const xvector<utype>&) __xprototype;

  template<class utype> xvector<utype>
    sec(const xvector<utype>&) __xprototype;

  template<class utype> xvector<utype>
    cosec(const xvector<utype>&) __xprototype;

  // HYPERBOLIC TRIGONOMETRIC OPERATIONS

  template<class utype> xvector<utype>
    sinh(const xvector<utype>&) __xprototype;

  template<class utype> xvector<utype>
    asinh(const xvector<utype>&) __xprototype;

  template<class utype> xvector<utype>
    cosh(const xvector<utype>&) __xprototype;

  template<class utype> xvector<utype>
    acosh(const xvector<utype>&) __xprototype;

  template<class utype> xvector<utype>
    tanh(const xvector<utype>&) __xprototype;

  template<class utype> xvector<utype>
    atanh(const xvector<utype>&) __xprototype;

  template<class utype> xvector<utype>
    cotanh(const xvector<utype>&) __xprototype;

  template<class utype> xvector<utype>
    acotanh(const xvector<utype>&) __xprototype;

  template<class utype> xvector<utype>
    sech(const xvector<utype>&) __xprototype;

  template<class utype> xvector<utype>
    cosech(const xvector<utype>&) __xprototype;

  // TRIGONOMETRIC OPERATIONS BETWEEN TWO VECTORS

  template<class utype> double   // cos of angle between two vectors
    cos(const xvector<utype>&,const xvector<utype>&) __xprototype;
  template<class utype> double   // cos of angle between two vectors
    getcos(const xvector<utype>&,const xvector<utype>&) __xprototype;
  
  template<class utype> double   // sin of angle between two vectors
    sin(const xvector<utype>&,const xvector<utype>&) __xprototype;
  template<class utype> double   // sin of angle between two vectors
    getsin(const xvector<utype>&,const xvector<utype>&) __xprototype;

  template<class utype> double   // angle between two vectors in radiants
    angle(const xvector<utype>&,const xvector<utype>&) __xprototype;
  template<class utype> double   // angle between two vectors in radiants
    getangle(const xvector<utype>&,const xvector<utype>&) __xprototype;

  template<class utype> double   // angle between three vectors in radiants
    angle(const xvector<utype>&,const xvector<utype>&,const xvector<utype>&) __xprototype;
  template<class utype> double   // angle between three  vectors in radiants
    getangle(const xvector<utype>&,const xvector<utype>&,const xvector<utype>&) __xprototype;

  template<class utype> bool //determine if two vectors are collinear //CO 180409
    isCollinear(const xvector<utype>& v0,const xvector<utype>& v1,const utype& tol=_AUROSTD_XVECTOR_TOLERANCE_IDENTITY_) __xprototype;

  template<class utype> xvector<utype> //get centroid of data points //CO 180409
    getCentroid(const vector<xvector<utype> >& points) __xprototype;

  // TRIGONOMETRIC OPERATIONS BETWEEN TWO GENERAL VECTORS //CO 180409
  template<class utype> xvector<double> 
    getGeneralAngles(const xvector<utype>& vec,const utype& tol=_AUROSTD_XVECTOR_TOLERANCE_IDENTITY_); //CO 180409

  template<class utype> double
    getGeneralAngle(const xvector<utype>& _vec,int _i,const utype& tol=_AUROSTD_XVECTOR_TOLERANCE_IDENTITY_); //CO 180409

  template<class utype> xvector<utype>
    getGeneralNormal(const vector<xvector<utype> >& _directive_vectors); //CO 180409

  // SIMPLE SORT ROUTINES
  template<class utype> xvector<utype>  // WRAP TO SHELL SHORT
    sort(const xvector<utype>& a) __xprototype;

  template<class utype> xvector<utype> // SHELLSORT
    shellsort(const xvector<utype>& a) __xprototype;

  template<class utype> xvector<utype> // HEAPSORT
    heapsort(const xvector<utype>& a) __xprototype;
  
  template<class utype> xvector<utype>  // QUICKSORT
    quicksort(const xvector<utype>&) __xprototype;
  
  template<class utype1, class utype2> void  // QUICKSORT
    quicksort2(unsigned long n, xvector<utype1>&, xvector<utype2>&) __xprototype;
  
  template<class utype1, class utype2, class utype3> void   // QUICKSORT
    quicksort3(unsigned long n, xvector<utype1>&, xvector<utype2>&, xvector<utype3>&)
    __xprototype;
  
  template<class utype1, class utype2, class utype3> void   // QUICKSORT
    quicksort4(unsigned long n, xvector<utype1>&, xvector<utype2>&, xvector<utype3>&)
    __xprototype;
  
  template<class utype> xvector<utype>             // SSORT shortcup for SHELLSORT
    ssort(const xvector<utype>& a) { return shellsort(a);}
  
  template<class utype> xvector<utype>             // HPSORT shortcup for HEAPSORT  
    hpsort(const xvector<utype>& a) { return heapsort(a);}
  
  //  template<class utype> xvector<utype>             // QSORT shortcup for QUICKSORT  
  //  sort(const xvector<utype>& a) { return heapsort(a);}
  
  template<class utype> xvector<utype>                 // SORT shortcup for HPSORT  
    qsort(const xvector<utype>& a) { return quicksort(a);}

  template<class utype1, class utype2> void  // QUICKSORT
    sort2(unsigned long n, xvector<utype1>&a, xvector<utype2>&b) {
    quicksort2(n,a,b);}

  template<class utype1, class utype2, class utype3> void  // QUICKSORT
    sort3(unsigned long n, xvector<utype1>&a, xvector<utype2>&b, xvector<utype3>&c) {
    quicksort3(n,a,b,c);}

  template<class utype1, class utype2, class utype3, class utype4> void  // QUICKSORT
    sort4(unsigned long n, xvector<utype1>&a, xvector<utype2>&b, xvector<utype3>&c, xvector<utype4>&d) {
    quicksort4(n,a,b,c,d);}

  //SOME STATS STUFF
  template<class utype> void getQuartiles(const xvector<utype>& _a,utype& q1,utype& q2,utype& q3);  //CO 171202
  template<class utype> utype getMAD(const xvector<utype>& _a,utype median=(utype)AUROSTD_NAN);   //CO 171202, absolute deviation around the median (MAD)
}

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

#endif  // _AUROSTD_XVECTOR_H_

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *                                                                         *
// ***************************************************************************

