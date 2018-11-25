// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *                                                                         *
// ***************************************************************************
// Written by Stefano Curtarolo 1994-2011

#ifndef _AUROSTD_XSCALAR_H_
#define _AUROSTD_XSCALAR_H_
#define _AUROSTD_XSCALAR_DEFAULT_SIZE_ 3

#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>

// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
#ifndef __xprototype
#define __xprototype __attribute__((const))
#endif

#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif

// ----------------------------------------------------------------------------
// _AUROSTD_XSCALAR_PROTOTYPES_

// ----------------------------------------------------------------------------
// ---------------------------------------------------------- XXXXXXXXXXXXXXXXX

// ----------------------------------------------------------------------------
// ---------------------------------------------------------- physics constants

#define rad2deg (180.0/3.14159265358979323846)
#define deg2rad (3.14159265358979323846/180.0)
#define angstrom2bohr (1/0.529177249)
#define bohr2angstrom (0.529177249)

//ME181020 - check out https://physics.nist.gov/cuu/Constants/index.html
//#define PI                    3.14159265359
#define PI                      3.14159265358979323846
#define TWOPI                   6.28318530717958647692
#define RTPI                    (sqrt(PI))
#define C_VACUUM                2.99792458E+8                   // m/s
#define EPS_VACUUM              8.854187817E-12                 // C/(N-m2)
#define MU_VACUUM               (4.0*PI*1.0E-7)                 // N/A^2 (T^2m^3/J)
#define AMU2KILOGRAM            1.66054E-27
#define KILOGRAM2AMU            6.022137E+26
#define E_ELECTRON              1.60217662E-19                  // C
#define PLANCKSCONSTANT_h       6.62607E-34                     // J*s
#define PLANCKSCONSTANT_hbar    1.0545718E-34                   // J*s
#define PLANCKSCONSTANTEV_h     (PLANCKSCONSTANT_h/E_ELECTRON)  // eV*s
#define PLANCKSCONSTANTEV_hbar  (PLANCKSCONSTANTEV_h/TWOPI)     // eV*s
#define KBOLTZ                  1.3806504E-23                   // J/K
#define KBOLTZEV                (KBOLTZ/E_ELECTRON)             // eV/K //8.617343E-5
#define eV2K                    (11604.505)                     // 1eV=11604.505 K
#define meV2K                   (11.604505)                     // 1meV=11.604505 K
#define atom2mol                6.0221408E23                    //CO 180329
#define meVatom2kJmol           (E_ELECTRON*atom2mol/1.0e6)     //CO 180329

// ----------------------------------------------------------------------------
// ------------------------------------------------------------------ constants

#define _mm_epsilon      1.0e-10
#define _mm_e            2.71828182845904523536028747135266249
#define _mm_log2e        1.4426950408889634074
#define _mm_log10e       0.43429448190325182765
#define _mm_ln2          0.69314718055994530942
#define _mm_ln10         2.30258509299404568402
#define _mm_pi           3.14159265358979323846
#define pi               3.14159265358979323846
#define PI               3.14159265358979323846
#define _mm_pi_2         1.57079632679489661923
#define _mm_pi_4         0.78539816339744830962
#define _mm_1_pi         0.31830988618379067154
#define _mm_2_pi         0.63661977236758134308
#define _mm_2_sqrtpi     1.12837916709551257390
#define _mm_sqrt2        1.41421356237309504880
#define _mm_sqrt1_2      0.70710678118654752440

#define _AUROSTD_XSCALAR_TOLERANCE_IDENTITY_ 1.0e-6
#define _AUROSTD_XSCALAR_TOLERANCE_ROUNDOFF_ 1.0e-6

// ----------------------------------------------------------------------------
// ------------------------- primitives for template<class utype> xscalar<utype>

// ----------------------------------------------------------------------------
// ------------------------------------------------- unary and binary operators

//--------------------------------------------------------- extra on data types

/* namespace aurostd { */
/*   double atof(string str); */
/*   int atoi(string str); */
/*   long atol(string str); */
/*   // long long atoll(string str); */
/*   // long long atoq(string str); */
/* } */

namespace aurostd {
  // namespace aurostd
  //[CO 180729 OBSOLETE]template<class utype> bool _isfloat(utype) __xprototype;
  template<class utype> bool _isodd(utype) __xprototype;
  template<class utype> bool _iseven(utype) __xprototype;
  template<class utype> bool _isreal(utype) __xprototype;
  template<class utype> bool _iscomplex(utype) __xprototype;
  template<class utype> int _size(utype) __xprototype;
  template<class utype> utype abs(utype) __xprototype;
  template<class utype> utype sqrt(utype) __xprototype;
  template<class utype> utype sign(utype) __xprototype;
  template<class utype> utype mod(utype,utype) __xprototype;
  template<class utype> utype nint(utype) __xprototype;
  int gcd(int a,int b);
  template<class utype> bool _isinteger(utype,utype=(utype)0.01) __xprototype;

  //bool isinteger(bool x,bool tolerance=(bool)0.01);
  //bool isinteger(char x,char tolerance=(char)0.01);
  bool isinteger(int x,int tolerance=(int)0.01);
  bool isinteger(long x,long tolerance=(long)0.01);
  bool isinteger(uint x,uint tolerance=(uint)0.01);
  bool isinteger(float x,float tolerance=(float)0.01);
  bool isinteger(double x,double tolerance=(double)0.01);
  //bool isinteger(long int x,long int tolerance=(long int)0.01);
  bool isinteger(long long int x,long long int tolerance=(long long int)0.01);
  bool isinteger(unsigned long long int x,unsigned long long int tolerance=(unsigned long long int)0.01);
  bool isinteger(long double x,long double tolerance=(long double)0.01);

  template<class utype> utype fact(utype) __xprototype;
  template<class utype> utype factorial(utype) __xprototype;
  template<class utype> utype angle(utype,utype,utype,utype,utype,utype) __xprototype;
  template<class utype> utype angle(utype,utype,utype,utype) __xprototype;
  template<class utype> utype modulus(utype,utype,utype) __xprototype;
  template<class utype> utype modulus(utype,utype) __xprototype;
// [OBSOLETE]  template<class utype> utype round(utype) __xprototype;
// [OBSOLETE]  template<class utype> utype floor(utype) __xprototype;
// [OBSOLETE]  template<class utype> utype ceil(utype) __xprototype;
// [OBSOLETE]  template<class utype> utype trunc(utype) __xprototype;

  double ln(double);
  float lnf(float);
  //  long double lnl(long double);
  float ln(float);
  //  long double ln(long double);
  double log(double);
  float logf(float);
  float log(float);
  //  long double logl(long double);
  //  long double log(long double);
  double log10(double);
  float log10f(float);
  float log10(float);
  // long double log10l(long double);
  // long double log10(long double);
}

// ----------------------------------------------------------------------------
// _isfloat  _isfloat  _isfloat  _isfloat  _isfloat
namespace aurostd {
  // namespace aurostd
  bool _isfloat(bool) __xprototype;
  bool _isfloat(char) __xprototype;
  bool _isfloat(int) __xprototype;
  bool _isfloat(uint) __xprototype;
  bool _isfloat(float) __xprototype;
  bool _isfloat(double) __xprototype;
  bool _isfloat(long int) __xprototype;
  bool _isfloat(long long int) __xprototype;
  bool _isfloat(unsigned long long int) __xprototype;
  bool _isfloat(long double) __xprototype;
#ifdef _AUROSTD_XCOMPLEX_
  bool _isfloat(xcomplex<float>) __xprototype;
  bool _isfloat(xcomplex<double>) __xprototype;
  bool _isfloat(xcomplex<long double>) __xprototype;
#endif
  bool isfloat(const string& in); //CO 180729
}

// ----------------------------------------------------------------------------
// _iscomplex  _iscomplex  _iscomplex  _iscomplex  _iscomplex
namespace aurostd {
  // namespace aurostd
  bool _iscomplex(bool) __xprototype;
  bool _iscomplex(char) __xprototype;
  bool _iscomplex(int) __xprototype;
  bool _iscomplex(uint) __xprototype;
  bool _iscomplex(float) __xprototype;
  bool _iscomplex(double) __xprototype;
  bool _iscomplex(long int) __xprototype;
  bool _iscomplex(long long int) __xprototype;
  bool _iscomplex(unsigned long long int) __xprototype;
  bool _iscomplex(long double) __xprototype;
#ifdef _AUROSTD_XCOMPLEX_
  bool _iscomplex(xcomplex<float>) __xprototype;
  bool _iscomplex(xcomplex<double>) __xprototype;
  bool _iscomplex(xcomplex<long double>) __xprototype;
#endif
}

// ----------------------------------------------------------------------------
// _isreal  _isreal  _isreal  _isreal  _isreal
namespace aurostd {
  // namespace aurostd
  bool _isreal(bool) __xprototype;
  bool _isreal(char) __xprototype;
  bool _isreal(int) __xprototype;
  bool _isreal(uint) __xprototype;
  bool _isreal(float) __xprototype;
  bool _isreal(double) __xprototype;
  bool _isreal(long int) __xprototype;
  bool _isreal(long long int) __xprototype;
  bool _isreal(unsigned long long int) __xprototype;
  bool _isreal(long double) __xprototype;
#ifdef _AUROSTD_XCOMPLEX_
  bool _isreal(xcomplex<float>) __xprototype;
  bool _isreal(xcomplex<double>) __xprototype;
  bool _isreal(xcomplex<long double>) __xprototype;
#endif
}

// ----------------------------------------------------------------------------
// _iseven  _iseven  _iseven  _iseven  _iseven
namespace aurostd {
  // namespace aurostd
  bool _iseven(bool) __xprototype;
  bool _iseven(char) __xprototype;
  bool _iseven(int) __xprototype;
  bool _iseven(uint) __xprototype;
  bool _iseven(float) __xprototype;
  bool _iseven(double) __xprototype;
  bool _iseven(long int) __xprototype;
  bool _iseven(long long int) __xprototype;
  bool _iseven(unsigned long long int) __xprototype;
  bool _iseven(long double) __xprototype;
#ifdef _AUROSTD_XCOMPLEX_
  bool _iseven(xcomplex<float>) __xprototype;
  bool _iseven(xcomplex<double>) __xprototype;
  bool _iseven(xcomplex<long double>) __xprototype;
#endif
}

// ----------------------------------------------------------------------------
// _isodd  _isodd  _isodd  _isodd  _isodd
namespace aurostd {
  // namespace aurostd
  bool _isodd(bool) __xprototype;
  bool _isodd(char) __xprototype;
  bool _isodd(int) __xprototype;
  bool _isodd(uint) __xprototype;
  bool _isodd(float) __xprototype;
  bool _isodd(double) __xprototype;
  bool _isodd(long int) __xprototype;
  bool _isodd(long long int) __xprototype;
  bool _isodd(unsigned long long int) __xprototype;
  bool _isodd(long double) __xprototype;
#ifdef _AUROSTD_XCOMPLEX_
  bool _isodd(xcomplex<float>) __xprototype;
  bool _isodd(xcomplex<double>) __xprototype;
  bool _isodd(xcomplex<long double>) __xprototype;
#endif
}

// ----------------------------------------------------------------------------
// _real _real _real _real _real _real _real _real
namespace aurostd {
  // namespace aurostd
  bool _real(bool) __xprototype;
  char _real(char) __xprototype;
  int _real(int) __xprototype;
  uint _real(uint) __xprototype;
  float _real(float) __xprototype;
  double _real(double) __xprototype;
  long int _real(long int) __xprototype;
  long long int _real(long long int) __xprototype;
  unsigned long long int _real(unsigned long long int) __xprototype;
  long double _real(long double) __xprototype;
#ifdef _AUROSTD_XCOMPLEX_
  float _real(xcomplex<float>) __xprototype;
  double _real(xcomplex<double>) __xprototype;
  long double _real(xcomplex<long double>) __xprototype;
#endif
}

// ----------------------------------------------------------------------------
// sign sign sign sign sign sign sign sign
namespace aurostd {
  // namespace aurostd
  bool sign(bool) __xprototype;
  char sign(char) __xprototype;
  int sign(int) __xprototype;
  uint sign(uint) __xprototype;
  float sign(float) __xprototype;
  double sign(double) __xprototype;
  long int sign(long int) __xprototype;
  long long int sign(long long int) __xprototype;
  long double sign(long double) __xprototype;
#ifdef _AUROSTD_XCOMPLEX_
  xcomplex<float> sign(xcomplex<float>) __xprototype;
  xcomplex<double> sign(xcomplex<double>) __xprototype;
  xcomplex<long double> sign(xcomplex<long double>) __xprototype;
#endif
}

// ----------------------------------------------------------------------------
// signnozero signnozero signnozero signnozero signnozero signnozero signnozero signnozero
namespace aurostd {
  // namespace aurostd
  bool signnozero(bool) __xprototype;
  char signnozero(char) __xprototype;
  int signnozero(int) __xprototype;
  uint signnozero(uint) __xprototype;
  float signnozero(float) __xprototype;
  double signnozero(double) __xprototype;
  long int signnozero(long int) __xprototype;
  long long int signnozero(long long int) __xprototype;
  long double signnozero(long double) __xprototype;
#ifdef _AUROSTD_XCOMPLEX_
  xcomplex<float> signnozero(xcomplex<float>) __xprototype;
  xcomplex<double> signnozero(xcomplex<double>) __xprototype;
  xcomplex<long double> signnozero(xcomplex<long double>) __xprototype;
#endif
}

// ----------------------------------------------------------------------------
// nint nint nint nint nint nint nint nint
namespace aurostd {
  // namespace aurostd
  bool nint(bool) __xprototype;
  char nint(char) __xprototype;
  int nint(int) __xprototype;
  uint nint(uint) __xprototype;
  uint nint(uint) __xprototype;
  float nint(float) __xprototype;
  double nint(double) __xprototype;
  long int nint(long int) __xprototype;
  long long int nint(long long int) __xprototype;
  long long unsigned int nint(long long unsigned int) __xprototype;
  long double nint(long double) __xprototype;
#ifdef _AUROSTD_XCOMPLEX_
  xcomplex<float> nint(xcomplex<float>) __xprototype;
  xcomplex<double> nint(xcomplex<double>) __xprototype;
  xcomplex<long double> nint(xcomplex<long double>) __xprototype;
#endif
}

// ----------------------------------------------------------------------------
// abs  abs  abs  abs  abs
namespace aurostd {
  // namespace aurostd
  // ABS(X)
  // int abs(int x) __xprototype;
  char abs(char x) __xprototype;
  int abs(int x) __xprototype;
  uint abs(uint x) __xprototype;
  float abs(float x) __xprototype;
  double abs(double x) __xprototype;
  long int abs(long int x) __xprototype;
  long long int abs(long long int x) __xprototype;
  unsigned long int abs(unsigned long int x) __xprototype;
  unsigned long long int abs(unsigned long long int x) __xprototype;
  long double abs(long double x) __xprototype;
#ifdef _AUROSTD_XCOMPLEX_
  float abs(xcomplex<float> x) __xprototype;
  double abs(xcomplex<double> x) __xprototype;
  long double abs(xcomplex<long double> x) __xprototype;
#endif
}

// ----------------------------------------------------------------------------
// sqrt  sqrt  sqrt  sqrt  sqrt
namespace aurostd {
  // namespace aurostd
  // SQRT(X)
  // int sqrt(int x) __xprototype;
  char sqrt(char x) __xprototype;
  int sqrt(int x) __xprototype;
  uint sqrt(uint x) __xprototype;
  float sqrt(float x) __xprototype;
  double sqrt(double x) __xprototype;
  long int sqrt(long int x) __xprototype;
  long long int sqrt(long long int x) __xprototype;
  unsigned long long int sqrt(unsigned long long int x) __xprototype;
  long double sqrt(long double x) __xprototype;
#ifdef _AUROSTD_XCOMPLEX_
  //  float sqrt(xcomplex<float> x) __xprototype;
  //  double sqrt(xcomplex<double> x) __xprototype;
  long double sqrt(xcomplex<long double> x) __xprototype;
#endif
}

//--------------------------------------------------------------- isequal
namespace aurostd {
  // with const utype&
  template<class utype> bool identical(const utype&,const utype&,const utype&) __xprototype;
  template<class utype> bool identical(const utype&,const utype&) __xprototype;
  template<class utype> bool isdifferent(const utype&,const utype&,const utype&) __xprototype;
  template<class utype> bool isdifferent(const utype&,const utype&) __xprototype;
  template<class utype> bool isequal(const utype&,const utype&,const utype&) __xprototype;
  template<class utype> bool isequal(const utype&,const utype&) __xprototype;
  /* with utype
     template<class utype> bool identical(utype,utype,utype) __xprototype;
     template<class utype> bool identical(utype,utype) __xprototype;
     template<class utype> bool isdifferent(utype,utype,utype) __xprototype;
     template<class utype> bool isdifferent(utype,utype) __xprototype;
     template<class utype> bool isequal(utype,utype,utype) __xprototype;
     template<class utype> bool isequal(utype,utype) __xprototype;
  */
}

//--------------------------------------------------------------- extra min/max // __XEXTRA_MINMAX_CPP
namespace aurostd {
 // namespace aurostd
  template<class utype> utype min(utype,utype) __xprototype;  // __XEXTRA_MINMAX_CPP
  template<class utype> utype min(utype,utype,utype) __xprototype;  // __XEXTRA_MINMAX_CPP
  template<class utype> utype min(utype,utype,utype,utype) __xprototype;  // __XEXTRA_MINMAX_CPP
  template<class utype> utype min(utype,utype,utype,utype,utype) __xprototype;  // __XEXTRA_MINMAX_CPP
  template<class utype> utype min(utype,utype,utype,utype,utype,utype) __xprototype;  // __XEXTRA_MINMAX_CPP
  template<class utype> utype min(utype,utype,utype,utype,utype,utype,utype) __xprototype;  // __XEXTRA_MINMAX_CPP
  template<class utype> utype min(utype,utype,utype,utype,utype,utype,utype,utype) __xprototype;  // __XEXTRA_MINMAX_CPP
  template<class utype> utype min(utype,utype,utype,utype,utype,utype,utype,utype,utype) __xprototype;  // __XEXTRA_MINMAX_CPP
  template<class utype> utype min(utype,utype,utype,utype,utype,utype,utype,utype,utype,utype) __xprototype;  // __XEXTRA_MINMAX_CPP
  template<class utype> utype min(utype,utype,utype,utype,utype,utype,utype,utype,utype,utype,utype) __xprototype;  // __XEXTRA_MINMAX_CPP
  template<class utype> utype min(utype,utype,utype,utype,utype,utype,utype,utype,utype,utype,utype,utype) __xprototype;  // __XEXTRA_MINMAX_CPP
  template<class utype> utype max(utype,utype) __xprototype;  // __XEXTRA_MINMAX_CPP
  template<class utype> utype max(utype,utype,utype) __xprototype;  // __XEXTRA_MINMAX_CPP
  template<class utype> utype max(utype,utype,utype,utype) __xprototype;  // __XEXTRA_MINMAX_CPP
  template<class utype> utype max(utype,utype,utype,utype,utype) __xprototype;  // __XEXTRA_MINMAX_CPP
  template<class utype> utype max(utype,utype,utype,utype,utype,utype) __xprototype;  // __XEXTRA_MINMAX_CPP
  template<class utype> utype max(utype,utype,utype,utype,utype,utype,utype) __xprototype;  // __XEXTRA_MINMAX_CPP
  template<class utype> utype max(utype,utype,utype,utype,utype,utype,utype,utype) __xprototype;  // __XEXTRA_MINMAX_CPP
  template<class utype> utype max(utype,utype,utype,utype,utype,utype,utype,utype,utype) __xprototype;  // __XEXTRA_MINMAX_CPP
  template<class utype> utype max(utype,utype,utype,utype,utype,utype,utype,utype,utype,utype) __xprototype;  // __XEXTRA_MINMAX_CPP
  template<class utype> utype max(utype,utype,utype,utype,utype,utype,utype,utype,utype,utype,utype) __xprototype;  // __XEXTRA_MINMAX_CPP
  template<class utype> utype max(utype,utype,utype,utype,utype,utype,utype,utype,utype,utype,utype,utype) __xprototype;  // __XEXTRA_MINMAX_CPP
}

namespace aurostd {
  template<class utype> utype _roundoff(const utype& x,utype tolerance);
  int roundoff(int x,int tolerance);
  long roundoff(long x,long tolerance);
  uint roundoff(uint x,uint tolerance);
  float roundoff(float x,float tolerance);
  double roundoff(double x,double tolerance);
  long long int roundoff(long long int x,long long int tolerance);
  unsigned long int roundoff(unsigned long int x,unsigned long int tolerance);
  unsigned long long int roundoff(unsigned long long int x,unsigned long long int tolerance);
  long double roundoff(long double x,long double tolerance);
}

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

#endif  // _AUROSTD_XSCALAR_H_

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *                                                                         *
// ***************************************************************************

