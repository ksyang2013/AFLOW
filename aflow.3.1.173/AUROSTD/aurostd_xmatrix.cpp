// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2015           *
// *                                                                         *
// ***************************************************************************
// Written by Stefano Curtarolo 1994-2014
// fixed for g++ 4.5 on Mar11

#ifndef _AUROSTD_XMATRIX_CPP_
#define _AUROSTD_XMATRIX_CPP_

//#define _AUROSTD_XMATRIX_DEFAULT_SIZE_ 3
//#define _AUROSTD_XMATRIX_TOLERANCE_IDENTITY_ double(1.0e-6)
//#define _AUROSTD_XMATRIX_TOLERANCE_ROUNDOFF_ double(1.0e-6)

//#ifndef _AUROSTD_XMATRIX_TOLERANCE_ROUNDOFF_
//#define _AUROSTD_XMATRIX_TOLERANCE_ROUNDOFF_ 1.0e-6 //DX 171025
//#endif

#ifndef XXEND
#define XXEND 1
#endif

#ifndef _xmatrix_epsilon
#define _xmatrix_epsilon 1.0e-15
#endif

#define _exponential_convergence 1.0e-18

#ifndef _AUROSTD_XCOMPLEX_H_
#include "aurostd_xcomplex.h"
#endif
#ifndef _AUROSTD_XSCALAR_H_
#include "aurostd_xscalar.h"
#endif
#ifndef _AUROSTD_XVECTOR_H_
#include "aurostd_xvector.h"
#endif
#ifndef _AUROSTD_XMATRIX_H_
#include "aurostd_xmatrix.h"
#endif
#ifndef _AUROSTD_XTENSOR6_H_
#include "aurostd_xtensor.h"
#endif

// ----------------------------------------------------------------------------
// --------------------------------------------------------------- constructors
namespace aurostd {  // namespace aurostd
  template<class utype>                                    // constructor
  xmatrix<utype>::xmatrix(int nrh,int nch,int nrl,int ncl) {
    int i,j;
    lrows=std::min(nrl,nrh); //if(!nrh||!nch) lrows=0; this messes up convasp
    urows=std::max(nrl,nrh); //if(!nrh||!nch) urows=0; this messes up convasp
    lcols=std::min(ncl,nch); //if(!nrh||!nch) lcols=0; this messes up convasp
    ucols=std::max(ncl,nch); //if(!nrh||!nch) ucols=0; this messes up convasp
    rows=urows-lrows+1;      //if(!nrh||!nch) rows=0; this messes up convasp
    cols=ucols-lcols+1;      //if(!nrh||!nch) cols=0; this messes up convasp
    issquare=bool(rows == cols);
    isfloat=_isfloat((utype) 0);
    iscomplex=_iscomplex((utype) 0);
    size=(char) (sizeof(utype));
    msize=(long int) size*rows*cols;
#ifdef _XMATH_DEBUG_CONSTRUCTORS
    cerr << "M -> constructor:"
	 << "  lrows=" << lrows << ", urows=" << urows
	 << ", lcols=" << lcols << ", ucols=" << ucols
	 << ", rows="  << rows  << ", cols="  << cols << endl;
#endif
    if(msize>0) {
      //Pass **passes = new ( Pass (*[ 10 ]) );
      //Pass** passes = new Pass*[10];
      corpus=new utype *[rows+XXEND]; // TRY
      if(!corpus)
	{cerr << _AUROSTD_XLIBS_ERROR_ << "allocation failure 1 in xmatrix constructor (int,int,int,int)" << endl;exit(0);}
      corpus+= -lrows+ XXEND;
      corpus[lrows]= new utype[rows*cols+XXEND];
      if(!corpus[lrows])
	{cerr << _AUROSTD_XLIBS_ERROR_ << "allocation failure 2 in xmatrix constructor (int,int,int,int)" << endl;exit(0);}
      corpus[lrows] += -lcols +XXEND;
      for(i=lrows+1;i<=urows;i++)
	corpus[i]=corpus[i-1]+cols;
      for(i=lrows;i<=urows;i++)                      // clear
	for(j=lcols;j<=ucols;j++)                    // clear
	  corpus[i][j]=(utype) 0.0;                    // clear
    }
#ifdef _XMATH_DEBUG_CONSTRUCTORS
    cerr << "issquare=" << issquare << ", isfloat=" << isfloat << ", iscomplex=" << iscomplex
	 << ", sizeof=" << size << ", msize=" << msize << endl;
#endif
  }
}
  
namespace aurostd {  // namespace aurostd
  template<class utype>                                       // copy constructor
  xmatrix<utype>::xmatrix(const xmatrix<utype>& a) {
    int i,j;
    lrows=a.lrows;urows=a.urows;rows=urows-lrows+1;
    lcols=a.lcols;ucols=a.ucols;cols=ucols-lcols+1;
    issquare=bool(rows == cols);
    isfloat=_isfloat((utype) 0);  
    iscomplex=_iscomplex((utype) 0);
    size=(char) (sizeof(utype));
    msize=(long int) size*rows*cols;
#ifdef _XMATH_DEBUG_CONSTRUCTORS
    cerr << "M -> copy constructor:"
	 << "  lrows=" << lrows << ", urows=" << urows
	 << ", lcols=" << lcols << ", ucols=" << ucols
	 << ", rows="  << rows  << ", cols="  << cols << endl;
#endif
    if(msize>0) {
      corpus=new utype *[rows+XXEND];
      if(!corpus)
	{cerr << _AUROSTD_XLIBS_ERROR_ << "allocation failure 1 in xmatrix constructor (xmatrix)" << endl;exit(0);}
      corpus+= -lrows + XXEND;
      corpus[lrows]= new utype[rows*cols+XXEND];
      if(!corpus[lrows])
	{cerr << _AUROSTD_XLIBS_ERROR_ << "allocation failure 2 in xmatrix constructor (xmatrix)" << endl;exit(0);}
      corpus[lrows] += -lcols +XXEND;
      for(i=lrows+1;i<=urows;i++)
	corpus[i]=corpus[i-1]+cols;
      for(i=lrows;i<=urows;i++)          
	for(j=lcols;j<=ucols;j++)        
	  corpus[i][j]=(utype) a.corpus[i][j];
    }
#ifdef _XMATH_DEBUG_CONSTRUCTORS
    cerr << "issquare=" << issquare << ", isfloat=" << isfloat << ", iscomplex=" << iscomplex
	 << ", sizeof=" << size << ", msize=" << msize << endl;
#endif
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype>                                       // copy constructor
  xmatrix<utype>::xmatrix(int vrows,int vcols, utype* a) {  // a starts from 0..
    int i,j;
    lrows=1;urows=vrows;rows=urows-lrows+1;
    lcols=1;ucols=vcols;cols=ucols-lcols+1;
    issquare=bool(rows == cols);
    isfloat=_isfloat((utype) 0);  
    iscomplex=_iscomplex((utype) 0);
    size=(char) (sizeof(utype));
    msize=(long int) size*rows*cols;
#ifdef _XMATH_DEBUG_CONSTRUCTORS
    cerr << "M -> copy constructor:"
	 << "  lrows=" << lrows << ", urows=" << urows
	 << ", lcols=" << lcols << ", ucols=" << ucols
	 << ", rows="  << rows  << ", cols="  << cols << endl;
#endif
    if(msize>0) {
      corpus=new utype *[rows+XXEND];
      if(!corpus)
	{cerr << _AUROSTD_XLIBS_ERROR_ << "allocation failure 1 in xmatrix constructor (xmatrix)" << endl;exit(0);}
      corpus+= -lrows + XXEND;
      corpus[lrows]= new utype[rows*cols+XXEND];
      if(!corpus[lrows])
	{cerr << _AUROSTD_XLIBS_ERROR_ << "allocation failure 2 in xmatrix constructor (xmatrix)" << endl;exit(0);}
      corpus[lrows] += -lcols +XXEND;
      for(i=lrows+1;i<=urows;i++)
	corpus[i]=corpus[i-1]+cols;
      for(i=lrows;i<=urows;i++)          
	for(j=lcols;j<=ucols;j++)        
	  corpus[i][j]=(utype) a[(i-1)*ucols+(j-1)]; // a.corpus[i][j];  // LIKE FORTRAN
      //      delete [] a;
    }
#ifdef _XMATH_DEBUG_CONSTRUCTORS
    cerr << "issquare=" << issquare << ", isfloat=" << isfloat << ", iscomplex=" << iscomplex
	 << ", sizeof=" << size << ", msize=" << msize << endl;
#endif
  }
}

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------- destructor

namespace aurostd {  // namespace aurostd
  template<class utype>                                     // default destructor
  xmatrix<utype>::~xmatrix() {
    // cerr << "problem destructor xmatrix [1]" << endl;
    // free a xmatrix allocated with xmatrix()
#ifdef _XMATH_DEBUG_DESTRUCTORS
    cerr << "M -> default destructor:"
	 << "  lrows=" << lrows << ", urows=" << urows
	 << ", lcols=" << lcols << ", ucols=" << ucols
	 << ", rows="  << rows  << ", cols="  << cols << endl;
#endif
    if(msize>0) {
      delete [] (corpus[lrows]+lcols-XXEND);
      delete [] (corpus+lrows-XXEND);
    }
    // cerr << "problem destructor xmatrix [2]" << endl;
  }
}
// ----------------------------------------------------------------------------
// -------------------------------------------------------- assigment operators
  
namespace aurostd {  // namespace aurostd
  template<class utype>                                             // operator =
  xmatrix<utype>& xmatrix<utype>::operator=(const xmatrix<utype>& a) {
    int i,j;
    if(corpus!=a.corpus||rows!=a.rows||cols!=a.cols||lrows!=a.lrows||urows!=a.urows||lcols!=a.lcols||ucols!=a.ucols) {                // check  for a=a
      //CO 170803 - ODD CORNER CASE, same corpus and rows, but different lrows and urows
      //if(rows!=a.rows||cols!=a.cols) {
      if(rows!=a.rows||cols!=a.cols||lrows!=a.lrows||urows!=a.urows||lcols!=a.lcols||ucols!=a.ucols) {
	// if dims(this)!=dims(a) => build a new xmatrix !!!
	if(msize>0) {
	  delete [] (corpus[lrows]+lcols-XXEND);
	  delete [] (corpus+lrows-XXEND);
	}
	lrows=a.lrows;urows=a.urows;rows=a.rows;
	lcols=a.lcols;ucols=a.ucols;cols=a.cols;
	issquare=bool(rows == cols);
	isfloat=_isfloat((utype) 0);
	iscomplex=_iscomplex((utype) 0);
	size=(char) sizeof(utype);
	msize=(long int) size*rows*cols;
#ifdef _XMATH_DEBUG_OPERATORS
	cerr << "M -> operator =::"
	     << "  lrows=" << lrows << ", urows=" << urows
	     << ", lcols=" << lcols << ", ucols=" << ucols
	     << ", rows="  << rows  << ", cols="  << cols << endl;
#endif
	if(msize>0) {
	  corpus=new utype *[rows+XXEND];
	  if(!corpus)
	    {cerr << _AUROSTD_XLIBS_ERROR_ << "allocation failure 1 in xmatrix assignment" << endl;exit(0);}
	  corpus+= -lrows+ XXEND;
	  corpus[lrows]= new utype[rows*cols+XXEND];
	  if(!corpus[lrows])
	    {cerr << _AUROSTD_XLIBS_ERROR_ << "allocation failure 2 in xmatrix assignment" << endl;exit(0);}
	  corpus[lrows]+= -lcols+XXEND;
	  for(i=lrows+1;i<=urows;i++)
	    corpus[i]=corpus[i-1]+cols;
	}
#ifdef _XMATH_DEBUG_CONSTRUCTORS
	cerr << "issquare=" << issquare << ", isfloat=" << isfloat << ", iscomplex=" << iscomplex
	     << ", sizeof=" << size << ", msize=" << msize << endl;
#endif
      }
      for(j=0;j<cols;j++)
	for(i=0;i<rows;i++)
	  this->corpus[i+lrows][j+lcols] =
	    a.corpus[i+a.lrows][j+a.lcols];
    }
    return *this;
  }
}

// ----------------------------------------------------------------------------
// ------------------------------------------------------------ index operators
// ---------------------------------------------------------------- operator []

namespace aurostd {  // namespace aurostd
  template<class utype>
  // removed inline
  utype* xmatrix<utype>::operator[] (int ir) const {
#ifndef __XOPTIMIZE
    if(ir>urows)  {cerr << _AUROSTD_XLIBS_ERROR_ << "_xmatrix<utype>_rows_high ir=" << ir << ", lrows=" << lrows << ", hrows=" << urows << endl;exit(0);}
    if(ir<lrows) {cerr << _AUROSTD_XLIBS_ERROR_ << "_xmatrix<utype>_rows_low ir=" << ir << ", lrows=" << lrows << ", hrows=" << urows << endl;exit(0);}
#endif
    return corpus[ir];
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype>                                         // operator (i,j)
  // removed inline
  utype& xmatrix<utype>::operator()(int i,int j) const {
    //#ifndef XMATRIX_PERIODIC_BOUNDARY_CONDITIONS
#ifndef __XMATRIX_IGNORE_BOUNDARIES
    if(i>urows) {cerr << _AUROSTD_XLIBS_ERROR_ << "M -> i=" << i << " > urows=" << urows << endl;exit(0);}
    if(i<lrows) {cerr << _AUROSTD_XLIBS_ERROR_ << "M -> i=" << i << " < lrows=" << lrows << endl;exit(0);}
    if(j>ucols) {cerr << _AUROSTD_XLIBS_ERROR_ << "M -> j=" << j << " > ucols=" << ucols << endl;exit(0);}
    if(j<lcols) {cerr << _AUROSTD_XLIBS_ERROR_ << "M -> j=" << j << " < lcols=" << lcols << endl;exit(0);}
#endif // __XMATRIX_IGNORE_BOUNDARIES
    return corpus[i][j];
  }
  /*
  //#else
  # ifdef XMATH_WARNING
  # warning "XMATRIX_PERIODIC_BOUNDARY_CONDITIONS"
  # endif
  int ii=i,jj=j;
  // if(ii>urows) ii=lrows+mod(i-lrows,urows-lrows+1);
  // if(ii<lrows) ii=urows-mod(urows-i,urows-lrows+1);
  // if(jj>ucols) jj=lcols+mod(j-lcols,ucols-lcols+1);
  // if(jj<lcols) jj=ucols-mod(ucols-j,ucols-lcols+1);
  if(ii>urows) ii-=rows;
  if(ii<lrows) ii+=rows;
  if(jj>ucols) jj-=cols;
  if(jj<lcols) jj+=cols;
  return corpus[ii][jj];
  #endif
  }
  */
}
  
namespace aurostd {  // namespace aurostd
  template<class utype>                                         // operator (i)
  xvector<utype> xmatrix<utype>::operator()(int i) const {
    xvector<utype> out(lcols,ucols);
    for(int j=lcols;j<=ucols;j++)
      out(j)=corpus[i][j];
    return out;
  }
}

// ----------------------------------------------------------------------------
// ----------------------------------- index operators with boundary conditions

namespace aurostd {  // namespace aurostd
  template<class utype>                        // operator () boundary conditions
  // removed inline
  utype& xmatrix<utype>::operator()(int i,int j,bool bc) const {
    if(bc==BOUNDARY_CONDITIONS_NONE) {
#ifndef __XOPTIMIZE
      if(i>urows) {cerr << _AUROSTD_XLIBS_ERROR_ << "M -> i=" << i << " > urows=" << urows << endl;exit(0);}
      if(i<lrows) {cerr << _AUROSTD_XLIBS_ERROR_ << "M -> i=" << i << " < lrows=" << lrows << endl;exit(0);}
      if(j>ucols) {cerr << _AUROSTD_XLIBS_ERROR_ << "M -> j=" << j << " > ucols=" << ucols << endl;exit(0);}
      if(j<lcols) {cerr << _AUROSTD_XLIBS_ERROR_ << "M -> j=" << j << " < lcols=" << lcols << endl;exit(0);}
#endif
      return corpus[i][j];
    }
    if(bc==BOUNDARY_CONDITIONS_PERIODIC) {
      int ii=i,jj=j;
      if(ii==urows+1) ii=lrows; // fast switching
      if(ii==lrows-1) ii=urows; // fast switching
      if(ii>urows) ii=lrows+mod(i-lrows,urows-lrows+1);
      if(ii<lrows) ii=urows-mod(urows-i,urows-lrows+1);
      if(jj==ucols+1) jj=lcols; // fast switching
      if(jj==lcols-1) jj=ucols; // fast switching
      if(jj>ucols) jj=lcols+mod(j-lcols,ucols-lcols+1);
      if(jj<lcols) jj=ucols-mod(ucols-j,ucols-lcols+1);
#ifndef __XOPTIMIZE
      if(ii>urows) {cerr << _AUROSTD_XLIBS_ERROR_ << "V -> ii=" << ii << " > urows" << urows << " <<  BC=" << bc << endl;exit(0);}
      if(ii<lrows) {cerr << _AUROSTD_XLIBS_ERROR_ << "V -> ii=" << ii << " < lrows" << lrows << " <<  BC=" << bc << endl;exit(0);}
      if(jj>ucols) {cerr << _AUROSTD_XLIBS_ERROR_ << "V -> jj=" << jj << " > ucols" << ucols << " <<  BC=" << bc << endl;exit(0);}
      if(jj<lcols) {cerr << _AUROSTD_XLIBS_ERROR_ << "V -> jj=" << jj << " < lcols" << lcols << " <<  BC=" << bc << endl;exit(0);}
#endif
      return corpus[ii][jj];
    }
  }
}

// ----------------------------------------------------------------------------
// ------------------------------------------------------- math unary operators

// -------------------------------------------------- operator xmatrix += xmatrix
namespace aurostd {  // namespace aurostd
  template<class utype> xmatrix<utype>&
  // removed inline
  xmatrix<utype>::operator +=(const xmatrix<utype>& r)
  {
#ifdef _XMATH_DEBUG_OPERATORS
    printf("M -> operator +=: ");
    printf("this->lrows=%i, this->urows=%i, ",this->lrows,this->urows);
    printf("this->lcols=%i, this->ucols=%i\n",this->lcols,this->ucols);
    printf("                 ");
    printf("r.lrows=%i, r.urows=%i, ",r.lrows,r.urows);
    printf("r.lcols=%i, r.ucols=%i\n",r.lcols,r.ucols);
#endif
    if(this->rows!=r.rows||this->cols!=r.cols)
      {cerr << _AUROSTD_XLIBS_ERROR_ << "failure in xmatrix operator+=" << endl
	    << " (this->rows!=r.rows||this->cols!=r.cols)=FALSE" << endl;exit(0);}
    for(int i=0;i<rows;i++)
      for(int j=0;j<cols;j++)
	corpus[i+lrows][j+lcols]+=r[i+r.lrows][j+r.lcols];
    return *this;
  }
}

// -------------------------------------------------- operator xmatrix -= xmatrix
namespace aurostd {  // namespace aurostd
  template<class utype> xmatrix<utype>&                  
  // removed inline
  xmatrix<utype>::operator -=(const xmatrix<utype>& r)
  {
#ifdef _XMATH_DEBUG_OPERATORS
    printf("M -> operator -=: ");
    printf("this->lrows=%i, this->urows=%i, ",this->lrows,this->urows);
    printf("this->lcols=%i, this->ucols=%i\n",this->lcols,this->ucols);
    printf("                 ");
    printf("r.lrows=%i, r.urows=%i, ",r.lrows,r.urows);
    printf("r.lcols=%i, r.ucols=%i\n",r.lcols,r.ucols);
#endif
    if(this->rows!=r.rows||this->cols!=r.cols)
      {cerr << _AUROSTD_XLIBS_ERROR_ << "failure in xmatrix operator-=               " << endl
	    << " (this->rows!=r.rows||this->cols!=r.cols)=FALSE" << endl;exit(0);}
    for(int i=0;i<rows;i++)
      for(int j=0;j<cols;j++)
	corpus[i+lrows][j+lcols]-=r[i+r.lrows][j+r.lcols];
    return *this;
  }
}

// -------------------------------------------------- operator xmatrix *= xmatrix
namespace aurostd {  // namespace aurostd
  template<class utype> xmatrix<utype>&                  
  // removed inline
  xmatrix<utype>::operator *=(const xmatrix<utype>& b)
  {
#ifdef _XMATH_DEBUG_OPERATORS
    printf("M -> operator *=: ");
    printf("this->lrows=%i, this->urows=%i, ",this->lrows,this->urows);
    printf("this->lcols=%i, this->ucols=%i\n",this->lcols,this->ucols);
    printf("                 ");
    printf("b.lrows=%i, b.urows=%i, ",b.lrows,b.urows);
    printf("b.lcols=%i, b.ucols=%i\n",b.lcols,b.ucols);
#endif
    if(!this->issquare||!b.issquare||this->rows!=b.rows)
      {cerr << _AUROSTD_XLIBS_ERROR_ << "failure in xmatrix operator*=: "
	    <<  "defined only for square xmatrixes with equal dimensions \n" << endl;exit(0);}
      
    xmatrix<utype> a(this->urows,this->ucols,this->lrows,this->lcols);
    int i,j,k,ii,jj,kk;
    utype *bk,*ai,aik,*thisi;
      
    for(i=this->lrows;i<=this->urows;i++)
      for(j=this->lcols;j<=this->ucols;j++) {
	a.corpus[i][j]=this->corpus[i][j];
	this->corpus[i][j]=(utype) 0;
      }
    for(i=this->lrows,ii=a.lrows;i<=this->urows;i++,ii++) {
      thisi=this->corpus[i];
      ai=a[ii];
      for(k=a.lcols,kk=b.lcols;k<=a.ucols;k++,kk++) {
	bk=b[kk];
	aik=ai[k];
	for(j=this->lrows,jj=b.lcols;j<=this->urows;j++,jj++)
	  thisi[j]+=aik*bk[jj];
      }
    }
    return *this;
  }
}

// ----------------------------------------------------------- operator +xmatrix
namespace aurostd {  // namespace aurostd
  template<class utype>                          
  xmatrix<utype> operator+(const xmatrix<utype>& a) {
    return a;
  }
}

// ----------------------------------------------------------- operator -xmatrix
namespace aurostd {  // namespace aurostd
  template<class utype>                          
  xmatrix<utype> operator-(const xmatrix<utype>& a) {
    xmatrix<utype> c(a.urows,a.ucols,a.lrows,a.lcols);
    for (int i=a.lrows;i<=a.urows;i++)
      for (int j=a.lcols;j<=a.ucols;j++)
	c[i][j]=-a[i][j];
    return c;
  }
}

// ----------------------------------------------------------------------------
// ------------------------------------------------------ math binary operators

// ----------------------------------------------------------------------------
// --------------------------------------------------- operator xmatrix + xmatrix
namespace aurostd {  // namespace aurostd
  template<class utype>                          
  xmatrix<utype> operator+(const xmatrix<utype>& a,const xmatrix<utype>& b) {
    
#ifdef _XMATH_DEBUG_OPERATORS
    printf("M -> operator +: a.lrows=%i, a.urows=%i, a.lcols=%i, a.ucols=%i\n",a.lrows,a.urows,a.lcols,a.ucols);
    printf("M -> operator +: b.lrows=%i, b.urows=%i, b.lcols=%i, b.ucols=%i\n",b.lrows,b.urows,b.lcols,b.ucols);
#endif
    if(a.rows!=b.rows||a.cols!=b.cols)
      {cerr << _AUROSTD_XLIBS_ERROR_ << "failure in xmatrix operator+:"
	    << " (a.rows!=b.rows||a.cols!=b.cols)=FALSE" << endl;exit(0);}
    xmatrix<utype> c(a.rows,a.cols);
    int i,j;
    utype *bi,*ci,*ai;
    for(i=0;i<a.rows;i++) {
      ai=a[i+a.lrows];bi=b[i+b.lrows];ci=c[i+c.lrows];
      for(j=0;j<a.cols;j++)
	ci[j+c.lcols]=ai[j+a.lcols]+bi[j+b.lcols];}
    return c;
  }
}

// ----------------------------------------------------------------------------
// --------------------------------------------------- operator xmatrix - xmatrix
namespace aurostd {  // namespace aurostd
  template<class utype>
  xmatrix<utype> operator-(const xmatrix<utype>& a,const xmatrix<utype>& b) {
    
#ifdef _XMATH_DEBUG_OPERATORS
    printf("M -> operator +: a.lrows=%i, a.urows=%i, a.lcols=%i, a.ucols=%i\n",a.lrows,a.urows,a.lcols,a.ucols);
    printf("M -> operator +: b.lrows=%i, b.urows=%i, b.lcols=%i, b.ucols=%i\n",b.lrows,b.urows,b.lcols,b.ucols);
#endif
    if(a.rows!=b.rows||a.cols!=b.cols)
      {cerr << _AUROSTD_XLIBS_ERROR_ << "failure in xmatrix operator-:"
	    << " (a.rows!=b.rows||a.cols!=b.cols)=FALSE" << endl;exit(0);}
    xmatrix<utype> c(a.rows,a.cols);
    int i,j;
    utype *bi,*ci,*ai;
    for(i=0;i<a.rows;i++) {
      ai=a[i+a.lrows];bi=b[i+b.lrows];ci=c[i+c.lrows];
      for(j=0;j<a.cols;j++)
	ci[j+c.lcols]=ai[j+a.lcols]-bi[j+b.lcols];};
    return c;
  }
}

// ----------------------------------------------------------------------------
// --------------------------------------------------- operator xmatrix * xmatrix
namespace aurostd {  // namespace aurostd
  template<class utype>                  
  xmatrix<utype> operator*(const xmatrix<utype>& a,const xmatrix<utype>& b) {    
#ifdef _XMATH_DEBUG_OPERATORS
    printf("M -> operator *: a.lrows=%i, a.urows=%i, a.lcols=%i, a.ucols=%i\n",a.lrows,a.urows,a.lcols,a.ucols);
    printf("M -> operator *: b.lrows=%i, b.urows=%i, b.lcols=%i, b.ucols=%i\n",b.lrows,b.urows,b.lcols,b.ucols);
#endif    
    if(a.cols!=b.rows)
      {cerr << _AUROSTD_XLIBS_ERROR_ << "failure in xmatrix operator*: (a.cols!=b.rows)=FALSE" << endl;exit(0);}    
    xmatrix<utype> c(a.rows,b.cols);
    int i,j,k,ii,jj,kk;
    // register
    utype *bk,*ci,*ai,aik;    
    for(i=c.lrows,ii=a.lrows;i<=c.urows;i++,ii++) {
      ci=c[i];
      ai=a[ii];
      //for(k=a.lcols,kk=b.lcols;k<=a.ucols;k++,kk++) {
      for(k=a.lcols,kk=b.lrows;k<=a.ucols;k++,kk++) {
	bk=b[kk];
	aik=ai[k];
	for(j=c.lcols,jj=b.lcols;j<=c.ucols;j++,jj++)
	  ci[j]+=aik*bk[jj];
      }
    }    
    /*  
	for(i=c.lrows,ii=a.lrows;i<=c.urows;i++,ii++)          // 48% slower than the
	for(k=a.lcols,kk=b.lrows;k<=a.ucols;k++,kk++)        // previous optimized
	for(j=c.lcols,jj=b.lcols;j<=c.ucols;j++,jj++)      // routine
	c[i][j]+=a[ii][k]*b[kk][jj];	
	for(i=c.lrows,ii=a.lrows;i<=c.urows;i++,ii++)          // 66% slower than the
	for(k=a.lcols,kk=b.lrows;k<=a.ucols;k++,kk++)        // previous optimized
	for(j=c.lcols,jj=b.lcols;j<=c.ucols;j++,jj++)      // routine
	c(i,j)+=a(ii,k)*b(kk,jj);
    */
    return  c;
  }
}

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
namespace aurostd {  // namespace aurostd
  template<class utype>                               // operator xmatrix * xvector
  xvector<utype> operator*(const xmatrix<utype>& a,const xvector<utype>& b) {    
#ifdef _XMATH_DEBUG_OPERATORS
    printf("M -> operator *: a.lrows=%i, a.urows=%i, a.lcols=%i, a.ucols=%i\n",a.lrows,a.urows,a.lcols,a.ucols);
    printf("M -> operator *: b.lrows=%i, b.urows=%i \n",b.lrows,b.urows);
#endif    
    if(a.cols!=b.rows)
      {cerr << _AUROSTD_XLIBS_ERROR_ << "failure in xmatrix operator* xmatrix * xvector: (a.cols!=b.rows)=FALSE " << a.cols << " " << b.rows << endl;exit(0);}    
    xvector<utype> c(a.lrows,a.urows);    
    for(int i=a.lrows;i<=a.urows;i++)
      for(int j=a.lcols;j<=a.ucols;j++)
        //      c[i]+=a[i][j]*b[j-b.lrows+1];
        c(i)+=a(i,j)*b(j-b.lrows+1);   // check... the 1 might be wrong
    return  c;
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype>                               // operator xvector * xmatrix
  xvector<utype> operator*(const xvector<utype>& a,const xmatrix<utype>& b) {    
#ifdef _XMATH_DEBUG_OPERATORS
    printf("M -> operator *: a.lrows=%i, a.urows=%i \n",a.lrows,a.urows);
    printf("M -> operator *: b.lrows=%i, b.urows=%i, b.lcols=%i, b.ucols=%i\n",a.lrows,a.urows,a.lcols,a.ucols);
#endif    
    if(a.rows!=b.rows)
      {cerr << _AUROSTD_XLIBS_ERROR_ << "failure in xmatrix operator* xvector * xmatrix: (a.rows!=b.rows)=FALSE " << a.rows << " " << b.rows << endl;exit(0);}  
    xvector<utype> c(b.lcols,b.ucols);    
    for(int i=b.lcols;i<=b.ucols;i++)
      for(int j=a.lrows;j<=a.urows;j++)
        //      c[i]+=a[j]*b[j-a.lrows+b.lrows][i];
        c(i)+=a(j)*b(j-a.lrows+b.lrows,i);
    return  c;
  }
}

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
namespace aurostd {  // namespace aurostd
  template<class utype> xmatrix<utype>                 // operator xmatrix * scalar
  operator*(const utype s,const xmatrix<utype>& a) {
    xmatrix<utype> c(a.urows,a.ucols,a.lrows,a.lcols);
    for(int i=c.lrows;i<=c.urows;i++)
      for(int j=c.lcols;j<=c.ucols;j++)
	c[i][j]=(utype) a[i][j]*(utype) s;
    return c;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {  // namespace aurostd
  template<class utype> xmatrix<utype>                //  operator scalar * xmatrix
  operator*(const xmatrix<utype>& a,const utype s) {
    return s*a;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {  // namespace aurostd
  template<class utype> xmatrix<utype>                 // operator xmatrix / scalar
  operator/(const xmatrix<utype>& a,const utype s) {
    return (utype) ((utype)1/s)*a;                     //DX 1/15/17 - add utype to 1/s to account for xcomplex
  }
}

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// CONDITIONALS

namespace aurostd {  // namespace aurostd
  template<class utype> bool                             // is xmatrix == xmatrix ?
  __identical(const xmatrix<utype>& a,const xmatrix<utype>& b,const utype& _tol_,const char& _mode_) {
    if(a.rows!=b.rows) {cerr << _AUROSTD_XLIBS_ERROR_ << "failure in xmatrix function identical (xmatrix == xmatrix)[1]" << endl;exit(0);}
    if(a.cols!=b.cols) {cerr << _AUROSTD_XLIBS_ERROR_ << "failure in xmatrix function identical (xmatrix == xmatrix)[2]" << endl;exit(0);}
    bool output=TRUE;
    if(a.isfloat || a.iscomplex) {
      if(_mode_==1) { // relative tolerance
	for(int i=a.lrows,ii=b.lrows;i<=a.urows;i++,ii++)
	  for(int j=a.lcols,jj=b.lcols;j<=a.ucols;j++,jj++) {
	    output=output*(((abs(a[i][j]-b[ii][jj]))/(abs(a[i][j])+abs(b[ii][jj])+_tol_))<=_tol_);
	    if(output==FALSE) return (bool) output;
	  }
      }
      if(_mode_==0) { // absolute tolerance (faster)  DEFAULT
	for(int i=a.lrows,ii=b.lrows;i<=a.urows;i++,ii++)
	  for(int j=a.lcols,jj=b.lcols;j<=a.ucols;j++,jj++) {
	    output=output*(abs(a[i][j]-b[ii][jj])<=_tol_);
	    if(output==FALSE) return (bool) output;
	  }
      }
    } else {
      for(int i=a.lrows,ii=b.lrows;i<=a.urows;i++,ii++)
	for(int j=a.lcols,jj=b.lcols;j<=a.ucols;j++,jj++) {
	  output=output*(a[i][j]==b[ii][jj]);
	  if(output==FALSE) return (bool) output;
	}
    }
    return (bool) output;
  }
  
  template<class utype> bool                             // is xmatrix == xmatrix ?
  identical(const xmatrix<utype>& a,const xmatrix<utype>& b,const utype& _tol_,const char& _mode_) {
    if(a.rows!=b.rows) {cerr << _AUROSTD_XLIBS_ERROR_ << "failure in xmatrix function identical (xmatrix == xmatrix)[1]" << endl;exit(0);}
    if(a.cols!=b.cols) {cerr << _AUROSTD_XLIBS_ERROR_ << "failure in xmatrix function identical (xmatrix == xmatrix)[2]" << endl;exit(0);}
    if(a.isfloat || a.iscomplex) {
      if(_mode_==1) { // relative tolerance
	for(int i=a.lrows,ii=b.lrows;i<=a.urows;i++,ii++)
	  for(int j=a.lcols,jj=b.lcols;j<=a.ucols;j++,jj++)
	    if((abs(a[i][j]-b[ii][jj])/(abs(a[i][j])/2.0+abs(b[ii][jj])/2.0+_tol_))>=_tol_) return FALSE;
      }
      if(_mode_==0) { // absolute tolerance (faster)  DEFAULT
	for(int i=a.lrows,ii=b.lrows;i<=a.urows;i++,ii++)
	  for(int j=a.lcols,jj=b.lcols;j<=a.ucols;j++,jj++)
	    if(abs(a[i][j]-b[ii][jj])>=_tol_) return FALSE;
      }
    } else {
      for(int i=a.lrows,ii=b.lrows;i<=a.urows;i++,ii++)
	for(int j=a.lcols,jj=b.lcols;j<=a.ucols;j++,jj++)
	  if((a[i][j]!=b[ii][jj])) return FALSE;
    }
    return TRUE; // if FALSE has never found....
  }

  // namespace aurostd
  template<class utype> bool                             // is xmatrix == xmatrix ?
  identical(const xmatrix<utype>& a,const xmatrix<utype>& b,const utype& _tol_) {
    return (bool) identical(a,b,_tol_,(char) 0);  // relative
  }
  
  // namespace aurostd
  template<class utype> bool                             // is xmatrix == xmatrix ?
  rel_identical(const xmatrix<utype>& a,const xmatrix<utype>& b,const utype& _tol_) {
    return (bool) identical(a,b,_tol_,(char) 1);  // relative
  }
  
  // namespace aurostd
  template<class utype> bool                             // is xmatrix == xmatrix ?
  abs_identical(const xmatrix<utype>& a,const xmatrix<utype>& b,const utype& _tol_) {
    return (bool) identical(a,b,_tol_,(char) 0);  // absolute
  }
  
  // namespace aurostd
  template<class utype> bool                             // is xmatrix == xmatrix ?
  identical(const xmatrix<utype>& a,const xmatrix<utype>& b) {
    return (bool) identical(a,b,(utype) _AUROSTD_XMATRIX_TOLERANCE_IDENTITY_,(char) 0);
  }
  
  // namespace aurostd
  template<class utype> bool                             // is xmatrix == xmatrix ?
  operator==(const xmatrix<utype>& a,const xmatrix<utype>& b) {
    return (bool) identical(a,b,(utype) _AUROSTD_XMATRIX_TOLERANCE_IDENTITY_,(char) 0);
  }
  
  // namespace aurostd
  template<class utype> bool                             // is xmatrix != xmatrix ?
  isdifferent(const xmatrix<utype>& a,const xmatrix<utype>& b,const utype& _tol_) {
    return (bool) !identical(a,b,_tol_,(char) 0);
  }
  
  // namespace aurostd
  template<class utype> bool                             // is xmatrix != xmatrix ?
  isdifferent(const xmatrix<utype>& a,const xmatrix<utype>& b) {
    return (bool) !identical(a,b,(utype) _AUROSTD_XMATRIX_TOLERANCE_IDENTITY_,(char) 0);
  }
  
  // namespace aurostd
  template<class utype> bool                             // is xmatrix == xmatrix ?
  isequal(const xmatrix<utype>& a,const xmatrix<utype>& b,const utype& _tol_) {
    return (bool) identical(a,b,_tol_,(char) 0);
  }
  
  // namespace aurostd
  template<class utype> bool                             // is xmatrix == xmatrix ?
  isequal(const xmatrix<utype>& a,const xmatrix<utype>& b) {
    return (bool) identical(a,b,(utype) _AUROSTD_XMATRIX_TOLERANCE_IDENTITY_,(char) 0);
  }
  
  // namespace aurostd
  template<class utype> bool                             // is xmatrix != xmatrix ?
  operator!=(const xmatrix<utype>& a,const xmatrix<utype>& b) {
    return (bool) !identical(a,b,(utype) _AUROSTD_XMATRIX_TOLERANCE_IDENTITY_,(char) 0);
  }
  
  // namespace aurostd
  template<class utype> bool
  isinteger(const xmatrix<utype>& a,const utype& tol) {
    if(a.isfloat || a.iscomplex) {
      for(int i=a.lrows;i<=a.urows;i++)
	for(int j=a.lcols;j<=a.ucols;j++)
	  if(isinteger(a[i][j],tol)==FALSE) return FALSE;
    }
    return TRUE;
  }    
  
  // namespace aurostd
  //CO - START
  template<class utype> bool
  isidentity(const xmatrix<utype>& a) {
    if(a.rows!=a.cols) {cerr << _AUROSTD_XLIBS_ERROR_ << "failure in xmatrix issymmetric (must be square)" << endl;exit(0);}
    for(int i=a.lrows;i<=a.urows;i++)
      for(int j=a.lcols;j<=a.ucols;j++)
	if(i-a.lrows+1!=j-a.lcols+1){  // i != j
	  if(aurostd::abs(a[i][j]) >_AUROSTD_XMATRIX_TOLERANCE_IDENTITY_) return FALSE;
        }else{
	  if(aurostd::abs(1.0-a[i][j]) >_AUROSTD_XMATRIX_TOLERANCE_IDENTITY_) return FALSE;
        }
    return TRUE;
  }
  //CO - START

  // namespace aurostd
  template<class utype> bool
  isdiagonal(const xmatrix<utype>& a,const utype& _eps_) { //DX 171025
    for(int i=a.lrows;i<=a.urows;i++)
      for(int j=a.lcols;j<=a.ucols;j++)
	if(i-a.lrows+1!=j-a.lcols+1)  // i != j
	  if(aurostd::abs(a[i][j]) >_eps_) return FALSE;
    return TRUE;
  }
  
  // namespace aurostd
  template<class utype> bool
  issymmetric(const xmatrix<utype>& a) {
    if(a.rows!=a.cols) {cerr << _AUROSTD_XLIBS_ERROR_ << "failure in xmatrix issymmetric (must be square)" << endl;exit(0);}
    for(int i=a.lrows;i<=a.urows;i++)
      for(int j=a.lcols;j<=a.ucols;j++)
	if(aurostd::abs(a[i][j]-a[j][i])> _AUROSTD_XMATRIX_TOLERANCE_IDENTITY_) return FALSE;
    return TRUE;
  }
  
  // namespace aurostd
  template<class utype> bool
  isantisymmetric(const xmatrix<utype>& a) {
    if(a.rows!=a.cols) {cerr << _AUROSTD_XLIBS_ERROR_ << "failure in xmatrix isantisymmetric (must be square)" << endl;exit(0);}
    for(int i=a.lrows;i<=a.urows;i++)
      for(int j=a.lcols;j<=a.ucols;j++)
	if(aurostd::abs(a[i][j]-(-a[j][i]))> _AUROSTD_XMATRIX_TOLERANCE_IDENTITY_) return FALSE;
    return TRUE;
  }

  // namespace aurostd
  template<class utype> bool
  ishermitian(const xmatrix<utype>& a) {
    if(a.rows!=a.cols) {cerr << _AUROSTD_XLIBS_ERROR_ << "failure in xmatrix ishermitian (must be square)" << endl;exit(0);}
    for(int i=a.lrows;i<=a.urows;i++)
      for(int j=a.lcols;j<=a.ucols;j++) {
       	if(aurostd::abs(a[i][j]-aurostd::conj(a[j][i])) > _AUROSTD_XMATRIX_TOLERANCE_IDENTITY_) return FALSE;
      }
    return TRUE;
  }

  // namespace aurostd
  template<class utype> bool
  isantihermitian(const xmatrix<utype>& a) {
    if(a.rows!=a.cols) {cerr << _AUROSTD_XLIBS_ERROR_ << "failure in xmatrix isantihermitian (must be square)" << endl;exit(0);}
    for(int i=a.lrows;i<=a.urows;i++)
      for(int j=a.lcols;j<=a.ucols;j++) {
      	if(aurostd::abs(a[i][j]-(-aurostd::conj(a[j][i]))) > _AUROSTD_XMATRIX_TOLERANCE_IDENTITY_) return FALSE;
      }  
    return TRUE;
  }

}

// ****************************************************************************
// ------------------------------------------------------ xmatrix constrtuction
namespace aurostd {
  // reshape by colums
  template<class utype>
  xmatrix<utype> reshape(const xvector<utype>& v1) {
    xmatrix<utype> c(v1.rows,1);
    for (int i=c.lrows;i<=c.urows;i++)
      c(i,1)=v1(v1.lrows+(i-c.lrows+1));
    return c;
  }
  
  template<class utype>
  xmatrix<utype> reshape(const xvector<utype>& v1,const xvector<utype>& v2) {
    if(v1.rows!=v2.rows) {
      cerr << _AUROSTD_XLIBS_ERROR_ << "reshape(v1,v2): vectors must have same dimension "
	   << v1.rows << " " << v2.rows << " " << endl;
      exit(0);
    }
    xmatrix<utype> c(v1.rows,2);
    for (int i=c.lrows;i<=c.urows;i++) {
      c(i,1)=v1(v1.lrows+(i-c.lrows+1));
      c(i,2)=v2(v2.lrows+(i-c.lrows+1));
    }
    return c;
  }
  
  template<class utype> xmatrix<utype>
  reshape(const xvector<utype>& v1,const xvector<utype>& v2,const xvector<utype>& v3) {
    if(v1.rows!=v2.rows || v2.rows!=v3.rows) {
      cerr << _AUROSTD_XLIBS_ERROR_ << "reshape(v1,v2,v3): vectors must have same dimension "
	   << v1.rows << " " << v2.rows << " " << v3.rows << " " << endl;
      exit(0);
    }
    xmatrix<utype> c(v1.rows,3);
    for (int i=c.lrows;i<=c.urows;i++) {
      c(i,1)=v1(v1.lrows+(i-c.lrows+1));
      c(i,2)=v2(v2.lrows+(i-c.lrows+1));
      c(i,3)=v3(v3.lrows+(i-c.lrows+1));
    }
    return c;
  }
  
  template<class utype> xmatrix<utype>
  reshape(const xvector<utype>& v1,const xvector<utype>& v2,const xvector<utype>& v3,const xvector<utype>& v4) {
    if(v1.rows!=v2.rows || v2.rows!=v3.rows || v3.rows!=v4.rows) {
      cerr << _AUROSTD_XLIBS_ERROR_ << "reshape(v1,v2,v3,v4): vectors must have same dimension "
	   << v1.rows << " " << v2.rows << " " << v3.rows << " "
	   << v4.rows << " " << endl;
      exit(0);
    }
    xmatrix<utype> c(v1.rows,4);
    for (int i=c.lrows;i<=c.urows;i++) {
      c(i,1)=v1(v1.lrows+(i-c.lrows+1));
      c(i,2)=v2(v2.lrows+(i-c.lrows+1));
      c(i,3)=v3(v3.lrows+(i-c.lrows+1));
      c(i,4)=v4(v4.lrows+(i-c.lrows+1));
    }
    return c;
  }

  template<class utype> xmatrix<utype>
  reshape(const xvector<utype>& v1,const xvector<utype>& v2,const xvector<utype>& v3,const xvector<utype>& v4,const xvector<utype>& v5) {
    if(v1.rows!=v2.rows || v2.rows!=v3.rows || v3.rows!=v4.rows || v4.rows!=v5.rows) {
      cerr << _AUROSTD_XLIBS_ERROR_ << "reshape(v1,v2,v3,v4,v5): vectors must have same dimension "
	   << v1.rows << " " << v2.rows << " " << v3.rows << " "
	   << v4.rows << " " << v5.rows << " " << endl;
      exit(0);
    }
    xmatrix<utype> c(v1.rows,5);
    for (int i=c.lrows;i<=c.urows;i++) {
      c(i,1)=v1(v1.lrows+(i-c.lrows+1));
      c(i,2)=v2(v2.lrows+(i-c.lrows+1));
      c(i,3)=v3(v3.lrows+(i-c.lrows+1));
      c(i,4)=v4(v4.lrows+(i-c.lrows+1));
      c(i,5)=v5(v5.lrows+(i-c.lrows+1));
    }
    return c;
  }

  template<class utype> xmatrix<utype>
  reshape(const xvector<utype>& v1,const xvector<utype>& v2,const xvector<utype>& v3,const xvector<utype>& v4,const xvector<utype>& v5,const xvector<utype>& v6) {
    if(v1.rows!=v2.rows || v2.rows!=v3.rows || v3.rows!=v4.rows || v4.rows!=v5.rows || v5.rows!=v6.rows) {
      cerr << _AUROSTD_XLIBS_ERROR_ << "reshape(v1,v2,v3,v4,v5,v6): vectors must have same dimension "
	   << v1.rows << " " << v2.rows << " " << v3.rows << " "
	   << v4.rows << " " << v5.rows << " " << v6.rows << " " << endl;
      exit(0);
    }
    xmatrix<utype> c(v1.rows,6);
    for (int i=c.lrows;i<=c.urows;i++) {
      c(i,1)=v1(v1.lrows+(i-c.lrows+1));
      c(i,2)=v2(v2.lrows+(i-c.lrows+1));
      c(i,3)=v3(v3.lrows+(i-c.lrows+1));
      c(i,4)=v4(v4.lrows+(i-c.lrows+1));
      c(i,5)=v5(v5.lrows+(i-c.lrows+1));
      c(i,6)=v6(v6.lrows+(i-c.lrows+1));
    }
    return c;
  }

  // reshape by colums
  template<class utype>
  xmatrix<utype> reshape_cols(const xvector<utype>& v1) {
    return reshape(v1);
  }
  template<class utype>
  xmatrix<utype> reshape_cols(const xvector<utype>& v1,const xvector<utype>& v2) {
    return reshape(v1,v2);
  }
  template<class utype>
  xmatrix<utype> reshape_cols(const xvector<utype>& v1,const xvector<utype>& v2,const xvector<utype>& v3) {
    return reshape(v1,v2,v3);
  }
  template<class utype>
  xmatrix<utype> reshape_cols(const xvector<utype>& v1,const xvector<utype>& v2,const xvector<utype>& v3,const xvector<utype>& v4) {
    return reshape(v1,v2,v3,v4);
  }
  template<class utype>
  xmatrix<utype> reshape_cols(const xvector<utype>& v1,const xvector<utype>& v2,const xvector<utype>& v3,const xvector<utype>& v4,const xvector<utype>& v5) {
    return reshape(v1,v2,v3,v4,v5);
  }
  template<class utype>
  xmatrix<utype> reshape_cols(const xvector<utype>& v1,const xvector<utype>& v2,const xvector<utype>& v3,const xvector<utype>& v4,const xvector<utype>& v5,const xvector<utype>& v6) {
    return reshape(v1,v2,v3,v4,v5,v6);
  }

  // reshape by rows
  template<class utype>
  xmatrix<utype> reshape_rows(const xvector<utype>& v1) {
    xmatrix<utype> c(1,v1.rows);
    for (int i=c.lcols;i<=c.urows;i++)
      c(i,1)=v1(v1.lrows+(i-c.lcols+1));
    return c;
  }
  
  template<class utype>
  xmatrix<utype> reshape_rows(const xvector<utype>& v1,const xvector<utype>& v2) {
    if(v1.rows!=v2.rows) {
      cerr << _AUROSTD_XLIBS_ERROR_ << "reshape_rows(v1,v2): vectors must have same dimension "
	   << v1.rows << " " << v2.rows << " " << endl;
      exit(0);
    }
    xmatrix<utype> c(2,v1.rows);
    for (int i=c.lcols;i<=c.urows;i++) {
      c(i,1)=v1(v1.lrows+(i-c.lcols+1));
      c(i,2)=v2(v2.lrows+(i-c.lcols+1));
    }
    return c;
  }
  
  template<class utype> xmatrix<utype>
  reshape_rows(const xvector<utype>& v1,const xvector<utype>& v2,const xvector<utype>& v3) {
    if(v1.rows!=v2.rows || v2.rows!=v3.rows) {
      cerr << _AUROSTD_XLIBS_ERROR_ << "reshape_rows(v1,v2,v3): vectors must have same dimension "
	   << v1.rows << " " << v2.rows << " " << v3.rows << " " << endl;
      exit(0);
    }
    xmatrix<utype> c(3,v1.rows);
    for (int i=c.lcols;i<=c.urows;i++) {
      c(i,1)=v1(v1.lrows+(i-c.lcols+1));
      c(i,2)=v2(v2.lrows+(i-c.lcols+1));
      c(i,3)=v3(v3.lrows+(i-c.lcols+1));
    }
    return c;
  }
  
  template<class utype> xmatrix<utype>
  reshape_rows(const xvector<utype>& v1,const xvector<utype>& v2,const xvector<utype>& v3,const xvector<utype>& v4) {
    if(v1.rows!=v2.rows || v2.rows!=v3.rows || v3.rows!=v4.rows) {
      cerr << _AUROSTD_XLIBS_ERROR_ << "reshape_rows(v1,v2,v3,v4): vectors must have same dimension "
	   << v1.rows << " " << v2.rows << " " << v3.rows << " "
	   << v4.rows << " " << endl;
      exit(0);
    }
    xmatrix<utype> c(4,v1.rows);
    for (int i=c.lcols;i<=c.urows;i++) {
      c(i,1)=v1(v1.lrows+(i-c.lcols+1));
      c(i,2)=v2(v2.lrows+(i-c.lcols+1));
      c(i,3)=v3(v3.lrows+(i-c.lcols+1));
      c(i,4)=v4(v4.lrows+(i-c.lcols+1));
    }
    return c;
  }

  template<class utype> xmatrix<utype>
  reshape_rows(const xvector<utype>& v1,const xvector<utype>& v2,const xvector<utype>& v3,const xvector<utype>& v4,const xvector<utype>& v5) {
    if(v1.rows!=v2.rows || v2.rows!=v3.rows || v3.rows!=v4.rows || v4.rows!=v5.rows ) {
      cerr << _AUROSTD_XLIBS_ERROR_ << "reshape_rows(v1,v2,v3,v4,v5): vectors must have same dimension "
	   << v1.rows << " " << v2.rows << " " << v3.rows << " "
	   << v4.rows << " " << v5.rows << " " << endl;
      exit(0);
    }
    xmatrix<utype> c(5,v1.rows);
    for (int i=c.lcols;i<=c.urows;i++) {
      c(i,1)=v1(v1.lrows+(i-c.lcols+1));
      c(i,2)=v2(v2.lrows+(i-c.lcols+1));
      c(i,3)=v3(v3.lrows+(i-c.lcols+1));
      c(i,4)=v4(v4.lrows+(i-c.lcols+1));
      c(i,5)=v5(v5.lrows+(i-c.lcols+1));
    }
    return c;
  }

  template<class utype> xmatrix<utype>
  reshape_rows(const xvector<utype>& v1,const xvector<utype>& v2,const xvector<utype>& v3,const xvector<utype>& v4,const xvector<utype>& v5,const xvector<utype>& v6) {
    if(v1.rows!=v2.rows || v2.rows!=v3.rows || v3.rows!=v4.rows || v4.rows!=v5.rows ) {
      cerr << _AUROSTD_XLIBS_ERROR_ << "reshape_rows(v1,v2,v3,v4,v5,v6): vectors must have same dimension "
	   << v1.rows << " " << v2.rows << " " << v3.rows << " "
	   << v4.rows << " " << v5.rows << " " << v6.rows << " " << endl;
      exit(0);
    }
    xmatrix<utype> c(6,v1.rows);
    for (int i=c.lcols;i<=c.urows;i++) {
      c(i,1)=v1(v1.lrows+(i-c.lcols+1));
      c(i,2)=v2(v2.lrows+(i-c.lcols+1));
      c(i,3)=v3(v3.lrows+(i-c.lcols+1));
      c(i,4)=v4(v4.lrows+(i-c.lcols+1));
      c(i,5)=v5(v5.lrows+(i-c.lcols+1));
      c(i,6)=v6(v6.lrows+(i-c.lcols+1));
    }
    return c;
  }

}

// ****************************************************************************
// -------------------------------------------------------------- xmatrix example types

namespace aurostd {  // namespace aurostd
  //171008 - CO
  template<class utype> xmatrix<utype>
  eyes(int nrh,int nch,int nrl,int ncl) __xprototype {
    xmatrix<utype> a(nrh,nch,nrl,ncl);
    for (int i=a.lrows;i<=a.urows;i++){
      for (int j=a.lcols;j<=a.ucols;j++){
        if(i==j){a[i][j]=(utype)1;}
      }
    }
    return a;
  }
}

// ****************************************************************************
// -------------------------------------------------------------- xmatrix casts

namespace aurostd {  // namespace aurostd
  template<class utype>                                 // conversion to long double
  xmatrix<long double> xlongdouble(const xmatrix<utype> &a) {
    xmatrix<long double> c(a.urows,a.ucols,a.lrows,a.lcols);
    for (int i=a.lrows;i<=a.urows;i++)
      for (int j=a.lcols;j<=a.ucols;j++)
	c[i][j]=(long double) a[i][j];
    return c;
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype>                                 // conversion to double
  xmatrix<double> xdouble(const xmatrix<utype> &a) {
    xmatrix<double> c(a.urows,a.ucols,a.lrows,a.lcols);
    for (int i=a.lrows;i<=a.urows;i++)
      for (int j=a.lcols;j<=a.ucols;j++)
	c[i][j]=(double) a[i][j];
    return c;
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype>                                 // conversion to float
  xmatrix<float> xfloat(const xmatrix<utype> &a) {
    xmatrix<float> c(a.urows,a.ucols,a.lrows,a.lcols);
    for (int i=a.lrows;i<=a.urows;i++)
      for (int j=a.lcols;j<=a.ucols;j++)
	c[i][j]=(float) a[i][j];
    return c;
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype>                                 // conversion to long int
  xmatrix<long int> xlongint(const xmatrix<utype> &a) {
    xmatrix<long int> c(a.urows,a.ucols,a.lrows,a.lcols);
    for (int i=a.lrows;i<=a.urows;i++)
      for (int j=a.lcols;j<=a.ucols;j++)
	c[i][j]=(long int) a[i][j];
    return c;
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype>                                 // conversion to int
  xmatrix<int> xint(const xmatrix<utype> &a) {
    xmatrix<int> c(a.urows,a.ucols,a.lrows,a.lcols);
    for (int i=a.lrows;i<=a.urows;i++)
      for (int j=a.lcols;j<=a.ucols;j++)
	c[i][j]=(int) a[i][j];
    return c;
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype>                                 // conversion to char
  xmatrix<char> xchar(const xmatrix<utype> &a) {
    xmatrix<char> c(a.urows,a.ucols,a.lrows,a.lcols);
    for (int i=a.lrows;i<=a.urows;i++)
      for (int j=a.lcols;j<=a.ucols;j++)
	c[i][j]=(char) a[i][j];
    return c;
  }
}

namespace aurostd {                   // conversion to vector<vector<utype> >
  template<class utype> vector<vector<utype> >
  xmatrix2vectorvector(const xmatrix<utype>& xmat) {
    int isize=xmat.rows,jsize=xmat.cols;
    // vector<vector<utype> > mat; vector<utype> v; mat=vector<vector<utype> > (m,v); // by hand
    // vector<vector<utype> > vectorvector(isize,jsize);              // WORKS WITH gcc/g++ 4.2 and 4.1
    vector<vector<utype> > vectorvector(isize,vector<utype>(jsize));  // WORKS WITH gcc/g++ 4.3

    for(int i=0;i<isize;i++)
      for(int j=0;j<jsize;j++)
	vectorvector[i][j]=xmat(i+xmat.lrows,j+xmat.lcols);
    return vectorvector;
  }
}

namespace aurostd {                   // conversion to xmatrix<utype>
  template<class utype> xmatrix<utype>
  vectorvector2xmatrix(const vector<vector<utype> >& mat) {
    int isize=mat.size(),jsize=mat.at(0).size();
    xmatrix<utype> xmat(isize,jsize);
    for(int i=1;i<=isize;i++)
      for(int j=1;j<=jsize;j++)
	xmat(i,j)=mat.at(i-1).at(j-1);
    return xmat;
  }
}

// ****************************************************************************
// ----------------------------------------------------------------------------
namespace aurostd {  // namespace aurostd
  template<class utype>                               // function reset xmatrix<>
  void xmatrix<utype>::reset(void) {
#ifdef _XMATH_DEBUG_FUNCTIONS
    cout<<"M -> function reset: "
	<<" lrows="<<lrows<<" urows="<<urows<<" lcols="<<lcols<<" ucols="<<ucols<<endl;
#endif
    for(int i=lrows;i<=urows;i++)
      for(int j=lcols;j<=ucols;j++)
	corpus[i][j]=(utype) 0.0;
  }
  template<class utype>                               // function reset xmatrix<>
  void reset(xmatrix<utype>& a) {
#ifdef _XMATH_DEBUG_FUNCTIONS
    cout<<"M -> function reset: "
	<<" a.lrows="<<a.lrows<<" a.urows="<<a.urows<<" a.lcols="<<a.lcols<<" a.ucols="<<a.ucols<<endl;
#endif
    for(int i=a.lrows;i<=a.urows;i++)
      for(int j=a.lcols;j<=a.ucols;j++)
	a[i][j]=(utype) 0.0;
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype>                               // function clear xmatrix<>
  void xmatrix<utype>::clear(void) {
#ifdef _XMATH_DEBUG_FUNCTIONS
    cout<<"M -> function clear: "
	<<" lrows="<<lrows<<" urows="<<urows<<" lcols="<<lcols<<" ucols="<<ucols<<endl;
#endif
    for(int i=lrows;i<=urows;i++)
      for(int j=lcols;j<=ucols;j++)
	corpus[i][j]=(utype) 0.0;
  }
  template<class utype>                               // function clear xmatrix<>
  void clear(xmatrix<utype>& a) {
#ifdef _XMATH_DEBUG_FUNCTIONS
    cout<<"M -> function clear: "
	<<" a.lrows="<<a.lrows<<" a.urows="<<a.urows<<" a.lcols="<<a.lcols<<" a.ucols="<<a.ucols<<endl;
#endif
    for(int i=a.lrows;i<=a.urows;i++)
      for(int j=a.lcols;j<=a.ucols;j++)
	a[i][j]=(utype) 0.0;
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype>                                 // function set xmatrix<>
  void xmatrix<utype>::set(const utype& s) {
#ifdef _XMATH_DEBUG_FUNCTIONS
    cout<<"M -> function set: "
	<<" lrows="<<lrows<<" urows="<<urows<<" lcols="<<lcols<<" ucols="<<ucols<<endl;
#endif
    for(int i=lrows;i<=urows;i++)
      for(int j=lcols;j<=ucols;j++)
	corpus[i][j]=(utype) s;
  }
  template<class utype>                                 // function set xmatrix<>
  void set(xmatrix<utype>& a,const utype& s) {
#ifdef _XMATH_DEBUG_FUNCTIONS
    cout<<"M -> function set: "
	<<" a.lrows="<<a.lrows<<" a.urows="<<a.urows<<" a.lcols="<<a.lcols<<" a.ucols="<<a.ucols<<endl;
#endif
    for(int i=a.lrows;i<=a.urows;i++)
      for(int j=a.lcols;j<=a.ucols;j++)
	a[i][j]=(utype) s;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {  // namespace aurostd
  template<class utype>                           // function vector<xmatrix<>>
  xvector<utype> vector(const xmatrix<utype>& a) {
    int n=(a.rows*a.cols);
    xvector<utype> c(1,n);
    for(int i=0;i<n;i++) {
      c[i+1]=(utype) a(int(i/a.cols)+a.lrows,mod(i,a.cols)+a.lcols);
    }
    return c;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {  // namespace aurostd
  template<class utype>                                 // function det xmatrix<>
  utype det(const xmatrix<utype>& a) {
    /* returns the determinant **/
    if(a.rows != a.cols) {cerr << _AUROSTD_XLIBS_ERROR_ << "Matrix determinant: must have a.rows = a.cols" << endl;exit(0);}
    if(a.lrows!=1)       {cerr << _AUROSTD_XLIBS_ERROR_ << "Matrix determinant: must have a.lrows = 1" << endl;exit(0);}
    if(a.lcols!=1)       {cerr << _AUROSTD_XLIBS_ERROR_ << "Matrix determinant: must have a.lcols = 1" << endl;exit(0);}
    int size=a.rows;
    //  cerr << "DET CALL size="<<size<< endl;
    if(size==1) return (utype) a[1][1];
    if(size==2) return (utype) a[1][1]*a[2][2]-a[1][2]*a[2][1];
    if(size==3)
      return (utype) (a[1][1]*a[2][2]*a[3][3]+
		      a[1][2]*a[2][3]*a[3][1]+
		      a[1][3]*a[2][1]*a[3][2]-
		      a[1][3]*a[2][2]*a[3][1]-
		      a[1][2]*a[2][1]*a[3][3]-
		      a[1][1]*a[2][3]*a[3][2]);
    if(size==4)
      return (utype) (a[1][4]*a[2][3]*a[3][2]*a[4][1]-a[1][3]*a[2][4]*a[3][2]*a[4][1]-a[1][4]*a[2][2]*a[3][3]*a[4][1]+a[1][2]*a[2][4]*a[3][3]*a[4][1]+
		      a[1][3]*a[2][2]*a[3][4]*a[4][1]-a[1][2]*a[2][3]*a[3][4]*a[4][1]-a[1][4]*a[2][3]*a[3][1]*a[4][2]+a[1][3]*a[2][4]*a[3][1]*a[4][2]+
		      a[1][4]*a[2][1]*a[3][3]*a[4][2]-a[1][1]*a[2][4]*a[3][3]*a[4][2]-a[1][3]*a[2][1]*a[3][4]*a[4][2]+a[1][1]*a[2][3]*a[3][4]*a[4][2]+
		      a[1][4]*a[2][2]*a[3][1]*a[4][3]-a[1][2]*a[2][4]*a[3][1]*a[4][3]-a[1][4]*a[2][1]*a[3][2]*a[4][3]+a[1][1]*a[2][4]*a[3][2]*a[4][3]+
		      a[1][2]*a[2][1]*a[3][4]*a[4][3]-a[1][1]*a[2][2]*a[3][4]*a[4][3]-a[1][3]*a[2][2]*a[3][1]*a[4][4]+a[1][2]*a[2][3]*a[3][1]*a[4][4]+
		      a[1][3]*a[2][1]*a[3][2]*a[4][4]-a[1][1]*a[2][3]*a[3][2]*a[4][4]-a[1][2]*a[2][1]*a[3][3]*a[4][4]+a[1][1]*a[2][2]*a[3][3]*a[4][4]);
    utype out=(utype) 0;
    if(size>=5) {
      xmatrix<utype> b(size-1,size-1);
      for(int j=1;j<=size;j++) {
	for(int ib=1;ib<=size-1;ib++)                                         // make sub
	  for(int jb=1;jb<=size-1;jb++)                                       // make sub
	    if(jb<j)  b[ib][jb]=a[ib+1][jb]; else b[ib][jb]=a[ib+1][jb+1];   // make sub --- FASTER
	if(_isodd(j)) out+=a[1][j]*det(b); else out-=a[1][j]*det(b);          // get part --- FASTER
	//if(jb<j)  b(ib,jb)=a(ib+1,jb); else  b(ib,jb)=a(ib+1,jb+1);         // make sub --- SLOWER
	//if(_isodd(j)) out+=a(1,j)*det(b); else out-=a(1,j)*det(b);          // get part --- SLOWER
      }
      return (utype) out;
    }
    return (utype) out;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {  // namespace aurostd
  template<class utype>                       // function determinant xmatrix<>
  utype determinant(const xmatrix<utype>& a) {
    return (utype) det(a);
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {  // namespace aurostd
  template<class utype> xmatrix<utype>                       // minor submatrix
  submatrix(const xmatrix<utype>& a,const int& irow,const int& jcol) {
    //     if(irow<a.lrows || irow>a.urows || jcol<a.lcols || jcol>a.ucols) {
    //       cerr << _AUROSTD_XLIBS_ERROR_ << "submatrix, error in irows,icols: ";
    //       cerr << " a.urows= " << a.urows << " a.lrows= " << a.lrows;
    //       cerr << " a.ucols= " << a.ucols << " a.lcols= " << a.lcols;
    //       cerr << " irow=" << irow << " jcol=" << jcol << endl;
    //       exit(0);
    //     }
    if(irow>=a.lrows && irow<=a.urows && jcol>=a.lcols && jcol<=a.ucols) {
      xmatrix<utype> b(a.urows-1,a.ucols-1,a.lrows,a.lcols);
      for(int i=a.lrows;i<=a.urows;i++)
	for(int j=a.lcols;j<=a.ucols;j++) {
	  if(i<irow && j<jcol) b[i][j]=a[i][j];
	  if(i>irow && j<jcol) b[i-1][j]=a[i][j];
	  if(i<irow && j>jcol) b[i][j-1]=a[i][j];
	  if(i>irow && j>jcol) b[i-1][j-1]=a[i][j];
	}
      return b;
    }
    if((irow<a.lrows || irow>a.urows) && jcol>=a.lcols && jcol<=a.ucols) {
      xmatrix<utype> b(a.urows,a.ucols-1,a.lrows,a.lcols);
      for(int i=a.lrows;i<=a.urows;i++)
	for(int j=a.lcols;j<=a.ucols;j++) {
	  if(j<jcol) b[i][j]=a[i][j];
	  if(j>jcol) b[i][j-1]=a[i][j];
	}
      return b;
    }
    if(irow>=a.lrows && irow<=a.urows && (jcol<a.lcols || jcol>a.ucols)) {
      xmatrix<utype> b(a.urows-1,a.ucols,a.lrows,a.lcols);
      for(int i=a.lrows;i<=a.urows;i++)
	for(int j=a.lcols;j<=a.ucols;j++) {
	  if(i<irow) b[i][j]=a[i][j];
	  if(i>irow) b[i-1][j]=a[i][j];
	}
      return b;
    }
    return a;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {  // namespace aurostd
  template<class utype> utype                       // minor submatrix
  minordet(const xmatrix<utype>& a,const int& irow,const int& jcol) {
    return det(submatrix(a,irow,jcol));
  }
  template<class utype> utype                       // minor submatrix
  minordeterminant(const xmatrix<utype>& a,const int& irow,const int& jcol) {
    return determinant(submatrix(a,irow,jcol));
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {  // namespace aurostd
  template<class utype>                                 // function inverse xmatrix<>
  xmatrix<utype> inverse(const xmatrix<utype>& a) {
    /* returns the inverse **/
    if(a.rows != a.cols) {cerr << _AUROSTD_XLIBS_ERROR_ << "Matrix inverse: must have a.rows = a.cols" << endl;exit(0);}
    if(a.lrows!=1)       {cerr << _AUROSTD_XLIBS_ERROR_ << "Matrix inverse: must have a.lrows = 1" << endl;exit(0);}
    if(a.lcols!=1)       {cerr << _AUROSTD_XLIBS_ERROR_ << "Matrix inverse: must have a.lcols = 1" << endl;exit(0);}
    int size=a.rows;
    xmatrix<utype> b(a.rows,a.cols);
    //  cerr << "DET CALL size="<<size<< endl;
    utype adet=det(a);
    if(adet==(utype) 0)  {cerr << _AUROSTD_XLIBS_ERROR_ << "Matrix inverse: singular matrix" << endl;exit(0);}
    if(size==1) {b[1][1]=1/a[1][1]; return b;}
    if(size==2) {{cerr << _AUROSTD_XLIBS_ERROR_ << "Matrix inverse: 2x2 not written yet" << endl;exit(0);} return b;};
    if(size==3) {
      adet=det(a);
      /*
	b[1][1]=(+a[2][2]*a[3][3]-a[2][3]*a[3][2])/adet;   // with fast index []
	b[1][2]=(-a[1][2]*a[3][3]+a[1][3]*a[3][2])/adet;   // with fast index []
	b[1][3]=(+a[1][2]*a[2][3]-a[1][3]*a[2][2])/adet;   // with fast index []
	b[2][1]=(-a[2][1]*a[3][3]+a[2][3]*a[3][1])/adet;   // with fast index []
	b[2][2]=(+a[1][1]*a[3][3]-a[1][3]*a[3][1])/adet;   // with fast index []
	b[2][3]=(-a[1][1]*a[2][3]+a[1][3]*a[2][1])/adet;   // with fast index []
	b[3][1]=(+a[2][1]*a[3][2]-a[2][2]*a[3][1])/adet;   // with fast index []
	b[3][2]=(-a[1][1]*a[3][2]+a[1][2]*a[3][1])/adet;   // with fast index []
	b[3][3]=(+a[1][1]*a[2][2]-a[1][2]*a[2][1])/adet;   // with fast index []
      */
      b(1,1)=(+a(2,2)*a(3,3)-a(2,3)*a(3,2))/adet;   // with slow index ()
      b(1,2)=(-a(1,2)*a(3,3)+a(1,3)*a(3,2))/adet;   // with slow index ()
      b(1,3)=(+a(1,2)*a(2,3)-a(1,3)*a(2,2))/adet;   // with slow index ()
      b(2,1)=(-a(2,1)*a(3,3)+a(2,3)*a(3,1))/adet;   // with slow index ()
      b(2,2)=(+a(1,1)*a(3,3)-a(1,3)*a(3,1))/adet;   // with slow index ()
      b(2,3)=(-a(1,1)*a(2,3)+a(1,3)*a(2,1))/adet;   // with slow index ()
      b(3,1)=(+a(2,1)*a(3,2)-a(2,2)*a(3,1))/adet;   // with slow index ()
      b(3,2)=(-a(1,1)*a(3,2)+a(1,2)*a(3,1))/adet;   // with slow index ()
      b(3,3)=(+a(1,1)*a(2,2)-a(1,2)*a(2,1))/adet;   // with slow index ()
      return b;
    };
    if(size>=4) {cerr << _AUROSTD_XLIBS_ERROR_ << "Matrix inverse: NxN not written yet" << endl;exit(0);}
    return b;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {  // namespace aurostd
  template<class utype> xmatrix<utype>  // function roundoff clear small elements
  roundoff(const xmatrix<utype>& a,utype _tol_) {
    xmatrix<utype> c(a.urows,a.ucols,a.lrows,a.lcols);
    for(int i=c.lrows;i<=c.urows;i++)
      for(int j=c.lcols;j<=c.ucols;j++) {
	if(abs(a[i][j])<(utype) _tol_) c[i][j]=a[i][j]=(utype) 0.0; else c[i][j]=a[i][j];
	//	c[i][j]=nint(a[i][j]/_tol_)*_tol_;
      }	  
    return c;
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype> xmatrix<utype>  // function roundoff clear small elements
  roundoff(const xmatrix<utype>& a) {
    return roundoff(a,(utype) _AUROSTD_XMATRIX_TOLERANCE_ROUNDOFF_);
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {  // namespace aurostd
  template<class utype>                                 // function sum xmatrix<>
  utype sum(const xmatrix<utype>& a) {
#ifdef _XMATH_DEBUG_FUNCTIONS
    printf("M -> function sum: ");
    printf("a.lrows=%i, a.urows=%i, ",a.lrows,a.urows);
    printf("a.lcols=%i, a.ucols=%i\n",a.lcols,a.ucols);
#endif
    utype c=utype(0.0);
    for(int i=a.lrows;i<=a.urows;i++)
      for(int j=a.lcols;j<=a.ucols;j++)
	c+=a[i][j];
    return c;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {  // namespace aurostd
  template<class utype>                           // function sum_colum xmatrix<>
  xvector<utype> sum_column(const xmatrix<utype>& a) {
    xvector<utype> c(a.lcols,a.ucols);
    for(int j=a.lcols;j<=a.ucols;j++) {
      c[j]=(utype) 0.0;
      for(int i=a.lrows;i<=a.urows;i++)
	c[j]+=a[i][j];
    }
    return c;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {  // namespace aurostd
  template<class utype>                          // function mean_colum xmatrix<>
  xvector<utype> mean_column(const xmatrix<utype>& a) {
    xvector<utype> c(a.lcols,a.ucols);
    for(int j=a.lcols;j<=a.ucols;j++) {
      c[j]=(utype) 0.0;
      for(int i=a.lrows;i<=a.urows;i++)
	c[j]+=a[i][j];
      c[j]/=(a.urows-a.lrows+1);
    }
    return c;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {  // namespace aurostd
  template<class utype>                             // function sum_row xmatrix<>
  xvector<utype> sum_row(const xmatrix<utype>& a) {
    xvector<utype> c(a.lrows,a.urows);
    for(int j=a.lrows;j<=a.urows;j++) {
      c[j]=(utype) 0.0;
      for(int i=a.lcols;i<=a.ucols;i++)
	c[j]+=a[j][i];
    }
    return c;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {  // namespace aurostd
  template<class utype>                            // function mean_row xmatrix<>
  xvector<utype> mean_row(const xmatrix<utype>& a) {
    xvector<utype> c(a.lrows,a.urows);
    for(int j=a.lrows;j<=a.urows;j++) {
      c[j]=(utype) 0.0;
      for(int i=a.lcols;i<=a.ucols;i++)
	c[j]+=a[j][i];
      c[j]/=(a.ucols-a.lcols+1);
    }
    return c;
  }
}

// -------------------------------------------------------- functions of xmatrix
namespace aurostd {  // namespace aurostd
  template<class utype>                                 // function min xmatrix<>
  utype min(const xmatrix<utype>& a) {
#ifdef _XMATH_DEBUG_FUNCTIONS
    printf("M -> function min: ");
    printf("a.lrows=%i, a.urows=%i, ",a.lrows,a.urows);
    printf("a.lcols=%i, a.ucols=%i\n",a.lcols,a.ucols);
#endif
    utype c=a[a.lrows][a.lcols];
    for(int i=a.lrows;i<=a.urows;i++)
      for(int j=a.lcols;j<=a.ucols;j++)
	c = c < a[i][j] ? c:a[i][j];
    return c;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {  // namespace aurostd
  template<class utype>                                 // function min xmatrix<>
  utype min(const xmatrix<utype>& a,int& index_i,int& index_j) {
#ifdef _XMATH_DEBUG_FUNCTIONS
    printf("M -> function min: ");
    printf("a.lrows=%i, a.urows=%i, ",a.lrows,a.urows);
    printf("a.lcols=%i, a.ucols=%i\n",a.lcols,a.ucols);
#endif
    utype c=a[a.lrows][a.lcols];
    index_i=a.lrows,index_j=a.lcols;
    for(int i=a.lrows;i<=a.urows;i++)
      for(int j=a.lcols;j<=a.ucols;j++)
	if(a[i][j] < c) {
	  c = a[i][j];
	  index_i=i;
	  index_j=j;
	}
#ifdef _XMATH_DEBUG_FUNCTIONS
    printf("M -> function max: ");
    printf("index_i=%i, index_j=%i \n",index_i,index_j);
#endif
    return c;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {  // namespace aurostd
  template<class utype>                                 // function max xmatrix<>
  utype max(const xmatrix<utype>& a) {
#ifdef _XMATH_DEBUG_FUNCTIONS
    printf("M -> function max: ");
    printf("a.lrows=%i, a.urows=%i, ",a.lrows,a.urows);
    printf("a.lcols=%i, a.ucols=%i\n",a.lcols,a.ucols);
#endif
    utype c=a[a.lrows][a.lcols];
    for(int i=a.lrows;i<=a.urows;i++)
      for(int j=a.lcols;j<=a.ucols;j++)
	c = c > a[i][j] ? c:a[i][j];
    return c;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {  // namespace aurostd
  template<class utype>                                 // function max xmatrix<>
  utype max(const xmatrix<utype>& a,int& index_i,int& index_j) {
#ifdef _XMATH_DEBUG_FUNCTIONS
    printf("M -> function max: ");
    printf("a.lrows=%i, a.urows=%i, ",a.lrows,a.urows);
    printf("a.lcols=%i, a.ucols=%i\n",a.lcols,a.ucols);
#endif
    utype c=a[a.lrows][a.lcols];
    index_i=a.lrows,index_j=a.lcols;
    for(int i=a.lrows;i<=a.urows;i++)
      for(int j=a.lcols;j<=a.ucols;j++)
	if(a[i][j] > c) {
	  c = a[i][j];
	  index_i=i;
	  index_j=j;
	}
#ifdef _XMATH_DEBUG_FUNCTIONS
    printf("M -> function max: ");
    printf("index_i=%i, index_j=%i \n",index_i,index_j);
#endif
    return c;
  }
}

// ----------------------------------------------------------------------------
/*namespace aurostd {  
// namespace aurostd
template<class utype> double                               // spectral radius
spectral_radius(const xmatrix<utype>& a)
{
#ifdef _XMATH_DEBUG_FUNCTIONS
printf("M -> function spectral radius: ");
printf("a.lrows=%i, a.urows=%i, ",a.lrows,a.urows);
printf("a.lcols=%i, a.ucols=%i\n",a.lcols,a.ucols);
#endif
if(!a.issquare)
{cerr << _AUROSTD_XLIBS_ERROR_ << "failure in xmatrix spectral radius defined for square xmatrixes\n" << endl;exit(0);}
double out=0.0;
for(int i=a.lrows;i<=a.urows;i++)
if(abs(a[i][i])>out)
out=(double) abs(a[i][i]);
return out;
}
}
*/

// ----------------------------------------------------------------------------
namespace aurostd {  // namespace aurostd
  template<class utype> utype  //DX 1/15/17 - double to utype (needed for xcomplex)                                          // trace
  trace(const xmatrix<utype>& a) {
#ifdef _XMATH_DEBUG_FUNCTIONS
    printf("M -> function trace: ");
    printf("a.lrows=%i, a.urows=%i, ",a.lrows,a.urows);
    printf("a.lcols=%i, a.ucols=%i\n",a.lcols,a.ucols);
#endif
    if(!a.issquare)
      {cerr << _AUROSTD_XLIBS_ERROR_ << "failure in xmatrix trace defined for square xmatrixes\n" << endl;exit(0);}
    utype out=0.0; //DX 1/15/17 - double to utype (needed for xcomplex)
    for(int i=a.lrows;i<=a.urows;i++)
      out+=a[i][i];
    return out;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {  // namespace aurostd
  template<class utype> xmatrix<utype>                          // identity xmatrix
  identity(const xmatrix<utype>& a) {
#ifdef _XMATH_DEBUG_FUNCTIONS
    printf("M -> function identity xmatrix: ");
    printf("a.lrows=%i, a.urows=%i, ",a.lrows,a.urows);
    printf("a.lcols=%i, a.ucols=%i\n",a.lcols,a.ucols);
#endif
    if(!a.issquare)
      {cerr << _AUROSTD_XLIBS_ERROR_ << "failure in xmatrix identity defined for square xmatrixes [1]\n" << endl;exit(0);}
    for(int i=a.lrows;i<=a.urows;i++)
      for(int j=a.lcols;j<=a.ucols;j++)
	a[i][j]=utype(0.0);
    for(int i=a.lrows;i<=a.urows;i++)
      a[i][i]=utype(1.0);
    return a;
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype>                                 // identity xmatrix
  xmatrix<utype> identity(const utype& _type,const int& n,const int& m) {
    xmatrix<utype> a(n,m);
    if(!a.issquare)
      {cerr << _AUROSTD_XLIBS_ERROR_ << "failure in xmatrix identity defined for square xmatrixes [2]\n" << endl;exit(0);}
    for(int i=a.lrows;i<=a.urows;i++)
      for(int j=a.lcols;j<=a.ucols;j++)
	if(i==j) a[i][j]=(utype) 1.0; else a[i][j]=(utype) 0.0;
    return a;
    if(_type) {;}  // something phony to keep _type busy !
  }
}

/*
  namespace aurostd {  // namespace aurostd
  xmatrix<double>                          // identity_double xmatrix
  identity_double(const int& n,const int& m) {
  xmatrix<double> a(n,m);
  if(!a.issquare) {cerr << _AUROSTD_XLIBS_ERROR_ << "failure in xmatrix identity_double defined for square xmatrixes [3]\n" << endl;exit(0);}
  for(int i=a.lrows;i<=a.urows;i++)
  for(int j=a.lcols;j<=a.ucols;j++)
  if(i==j) a[i][j]=double(1.0); else a[i][j]=double(0.0);
  return a;
  }
  // namespace aurostd
  xmatrix<double>                          // identity_double xmatrix
  identity_double(const int& n) {
  xmatrix<double> a(n,n);
  for(int i=a.lrows;i<=a.urows;i++)
  for(int j=a.lcols;j<=a.ucols;j++)
  if(i==j) a[i][j]=double(1.0); else a[i][j]=double(0.0);
  return a;
  }
  // namespace aurostd
  xmatrix<float>                          // identity_float xmatrix
  identity_float(const int& n,const int& m) {
  xmatrix<float> a(n,m);
  if(!a.issquare) {cerr << _AUROSTD_XLIBS_ERROR_ << "failure in xmatrix identity_float defined for square xmatrixes [3]\n" << endl;exit(0);}
  for(int i=a.lrows;i<=a.urows;i++)
  for(int j=a.lcols;j<=a.ucols;j++)
  if(i==j) a[i][j]=float(1.0); else a[i][j]=float(0.0);
  return a;
  }
  // namespace aurostd
  xmatrix<float>                          // identity_float xmatrix
  identity_float(const int& n) {
  xmatrix<float> a(n,n);
  for(int i=a.lrows;i<=a.urows;i++)
  for(int j=a.lcols;j<=a.ucols;j++)
  if(i==j) a[i][j]=float(1.0); else a[i][j]=float(0.0);
  return a;
  }
  // namespace aurostd
  xmatrix<int>                          // identity_int xmatrix
  identity_int(const int& n,const int& m) {
  xmatrix<int> a(n,m);
  if(!a.issquare) {cerr << _AUROSTD_XLIBS_ERROR_ << "failure in xmatrix identity_int defined for square xmatrixes [3]\n" << endl;exit(0);}
  for(int i=a.lrows;i<=a.urows;i++)
  for(int j=a.lcols;j<=a.ucols;j++)
  if(i==j) a[i][j]=1; else a[i][j]=0;
  return a;
  }
  // namespace aurostd
  xmatrix<int>                          // identity_int xmatrix
  identity_int(const int& n) {
  xmatrix<int> a(n,n);
  for(int i=a.lrows;i<=a.urows;i++)
  for(int j=a.lcols;j<=a.ucols;j++)
  if(i==j) a[i][j]=1; else a[i][j]=0;
  return a;
  }
  // namespace aurostd
  xmatrix<char>                          // identity_char xmatrix
  identity_char(const int& n,const int& m) {
  xmatrix<char> a(n,m);
  if(!a.issquare) {cerr << _AUROSTD_XLIBS_ERROR_ << "failure in xmatrix identity_char defined for square xmatrixes [3]\n" << endl;exit(0);}
  for(int i=a.lrows;i<=a.urows;i++)
  for(int j=a.lcols;j<=a.ucols;j++)
  if(i==j) a[i][j]=1; else a[i][j]=0;
  return a;
  }
  // namespace aurostd
  xmatrix<char>                          // identity_char xmatrix
  identity_char(const int& n) {
  xmatrix<char> a(n,n);
  for(int i=a.lrows;i<=a.urows;i++)
  for(int j=a.lcols;j<=a.ucols;j++)
  if(i==j) a[i][j]=1; else a[i][j]=0;
  return a;
  }
  // namespace aurostd
  xmatrix<bool>                          // identity_bool xmatrix
  identity_bool(const int& n,const int& m) {
  xmatrix<bool> a(n,m);
  if(!a.issquare) {cerr << _AUROSTD_XLIBS_ERROR_ << "failure in xmatrix identity_bool defined for square xmatrixes [3]\n" << endl;exit(0);}
  for(int i=a.lrows;i<=a.urows;i++)
  for(int j=a.lcols;j<=a.ucols;j++)
  if(i==j) a[i][j]=1; else a[i][j]=0;
  return a;
  }
  // namespace aurostd
  xmatrix<bool>                          // identity_bool xmatrix
  identity_bool(const int& n) {
  xmatrix<bool> a(n,n);
  for(int i=a.lrows;i<=a.urows;i++)
  for(int j=a.lcols;j<=a.ucols;j++)
  if(i==j) a[i][j]=1; else a[i][j]=0;
  return a;
  }
  }
*/

// ----------------------------------------------------------------------------
namespace aurostd {  // namespace aurostd
  template<class utype> xmatrix<utype>                             // trasp xmatrix
  trasp(const xmatrix<utype>& a)
  {
#ifdef _XMATH_DEBUG_FUNCTIONS
    printf("M -> function traspose xmatrix: ");
    printf("a.lrows=%i, a.urows=%i, ",a.lrows,a.urows);
    printf("a.lcols=%i, a.ucols=%i\n",a.lcols,a.ucols);
#endif
    xmatrix<utype> out(a.ucols,a.urows,a.lcols,a.lrows);
    for(int i=a.lrows;i<=a.urows;i++)
      for(int j=a.lcols;j<=a.ucols;j++)
	out[j][i]=a[i][j];
    return out;
  } // BUG aggiungere complesso coniugato se is.float
}

// ----------------------------------------------------------------------------
namespace aurostd {  // namespace aurostd
  template<class utype> xmatrix<utype>                             // trasp xvector
  trasp(const xvector<utype>& a)
  {
#ifdef _XMATH_DEBUG_FUNCTIONS
    printf("M -> function traspose xvector: ");
    printf("a.lrows=%i, a.urows=%i, ",a.lrows,a.urows);
#endif
    xmatrix<utype> out(a.urows,1,a.lrows,1);
    for(int i=a.lrows;i<=a.urows;i++)
      out[i][1]=a[i];
    return out;
  } // BUG aggiungere complesso coniugato se is.float
}

// ****************************************************************************
// ----------------------------------------------------------------------------
namespace aurostd {  // namespace aurostd
  template<class utype> xmatrix<utype>               // function shift_up xmatrix<>
  shift_up(const xmatrix<utype>& a) {
    utype aus;
    for(int j=a.lcols;j<=a.ucols;j++) {
      aus=a[a.lrows][j];
      for(int i=a.lrows;i<=a.urows-1;i++)
	a[i][j]=a[i+1][j];
      a[a.urows][j]=aus;
    }
    return a;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {  // namespace aurostd
  template<class utype> xmatrix<utype>             // function shift_down xmatrix<>
  shift_down(const xmatrix<utype>& a) {
    utype aus;
    for(int j=a.lcols;j<=a.ucols;j++) {
      aus=a[a.urows][j];
      for(int i=a.urows;i>=a.lrows+1;i--)
	a[i][j]=a[i-1][j];
      a[a.lrows][j]=aus;
    }
    return a;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {  // namespace aurostd
  template<class utype> xmatrix<utype>             // function shift_left xmatrix<>
  shift_left(const xmatrix<utype>& a) {
    utype aus;
    for(int i=a.lrows;i<=a.urows;i++) {
      aus=a[i][a.lcols];
      for(int j=a.lcols;j<=a.ucols-1;j++)
	a[i][j]=a[i][j+1];
      a[i][a.ucols]=aus;
    }
    return a;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {  // namespace aurostd
  template<class utype> xmatrix<utype>             // function shift_right xmatrix<>
  shift_right(const xmatrix<utype>& a) {
    utype aus;
    for(int i=a.lrows;i<=a.urows;i++) {
      aus=a[i][a.ucols];
      for(int j=a.ucols;j>=a.lcols+1;j++)
	a[i][j]=a[i][j+1];
      a[i][a.lcols]=aus;
    }
    return a;
  }
}

// ****************************************************************************
// ----------------------------------------------------------------------------
// SWAP THINGS
namespace aurostd {  // namespace aurostd
  template<class utype> void                                // swap_columns
  swap_cols(xmatrix<utype>& a,const int& i,const int& j) {  // swap_columns
    if(i<a.lcols || i>a.ucols) return; // nothing to do, out of boundaries
    if(j<a.lcols || j>a.ucols) return; // nothing to do, out of boundaries
    if(i==j) return; // nothing to do, no swap    
    utype temp;
    for(int k=a.lrows;k<=a.urows;k++) {
      temp=a[k][i];a[k][i]=a[k][j];a[k][j]=temp;
    }
  }
  
  template<class utype> void                                  // swap_columns
  swap_columns(xmatrix<utype>& a,const int& i,const int& j) {  // swap_columns
    swap_cols(a,i,j);
    //     if(i<a.lcols || i>a.ucols) return; // nothing to do, out of boundaries
    //     if(j<a.lcols || j>a.ucols) return; // nothing to do, out of boundaries
    //     if(i==j) return; // nothing to do, no swap    
    //     utype temp;
    //     for(int k=a.lrows;k<=a.urows;k++) {
    //       temp=a[k][i];a[k][i]=a[k][j];a[k][j]=temp;
    //     }
  }
  
  template<class utype> void                                // swap_rows
  swap_rows(xmatrix<utype>& a,const int& i,const int& j) {  // swap_rows
    if(i<a.lrows || i>a.urows) return; // nothing to do, out of boundaries
    if(j<a.lrows || j>a.urows) return; // nothing to do, out of boundaries
    if(i==j) return; // nothing to do, no swap    
    utype temp;
    for(int k=a.lcols;k<=a.ucols;k++) {
      temp=a[i][k];a[i][k]=a[j][k];a[j][k]=temp;
    }
  }
}

// ****************************************************************************
// ----------------------------------------------------------------------------
namespace aurostd {  // namespace aurostd
  template<class utype> xmatrix<utype>                          // sign xmatrix
  sign(const xmatrix<utype>& a) {
    xmatrix<utype> c(a.urows,a.ucols,a.lrows,a.lcols);
    for(int i=c.lrows;i<=c.urows;i++)
      for(int j=c.lcols;j<=c.ucols;j++)
	c[i][j]=aurostd::sign(a[i][j]);
    return c;
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype> xmatrix<utype>                          // nint xmatrix
  nint(const xmatrix<utype>& a) {
    xmatrix<utype> c(a.urows,a.ucols,a.lrows,a.lcols);
    for(int i=c.lrows;i<=c.urows;i++)
      for(int j=c.lcols;j<=c.ucols;j++)
	c[i][j]=aurostd::nint(a[i][j]);
    return c;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {  // namespace aurostd
  template<class utype> xmatrix<utype>                          // floor xmatrix
  floor(const xmatrix<utype>& a) {
    xmatrix<utype> c(a.urows,a.ucols,a.lrows,a.lcols);
    for(int i=c.lrows;i<=c.urows;i++)
      for(int j=c.lcols;j<=c.ucols;j++)
	c[i][j]=std::floor(a[i][j]);
    return c;
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype> xmatrix<utype>                          // trunc xmatrix
  trunc(const xmatrix<utype>& a) {
    xmatrix<utype> c(a.urows,a.ucols,a.lrows,a.lcols);
    for(int i=c.lrows;i<=c.urows;i++)
      for(int j=c.lcols;j<=c.ucols;j++)
	//	c[i][j]=std::trunc(a[i][j]);
	c[i][j]=trunc(a[i][j]);
    return c;
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype> xmatrix<utype>                          // round xmatrix
  round(const xmatrix<utype>& a) {
    xmatrix<utype> c(a.urows,a.ucols,a.lrows,a.lcols);
    for(int i=c.lrows;i<=c.urows;i++)
      for(int j=c.lcols;j<=c.ucols;j++)
	//	c[i][j]=std::round(a[i][j]);
 	c[i][j]=round(a[i][j]);
    return c;
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype> xmatrix<utype>                          // ceil xmatrix
  ceil(const xmatrix<utype>& a) {
    xmatrix<utype> c(a.urows,a.ucols,a.lrows,a.lcols);
    for(int i=c.lrows;i<=c.urows;i++)
      for(int j=c.lcols;j<=c.ucols;j++)
	c[i][j]=std::ceil(a[i][j]);
    return c;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {  // namespace aurostd
  template<class utype> xmatrix<utype>                           // mabs xmatrix
  mabs(const xmatrix<utype>& a) {
    xmatrix<utype> c(a.urows,a.ucols,a.lrows,a.lcols);
    for(int i=c.lrows;i<=c.urows;i++)
      for(int j=c.lcols;j<=c.ucols;j++)
	c[i][j]=aurostd::abs(a[i][j]);
    return c;
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype> xmatrix<utype>                           // abs xmatrix
  abs(const xmatrix<utype>& a) {
    xmatrix<utype> c(a.urows,a.ucols,a.lrows,a.lcols);
    for(int i=c.lrows;i<=c.urows;i++)
      for(int j=c.lcols;j<=c.ucols;j++)
	c[i][j]=aurostd::abs(a[i][j]);
    return c;
  }
}

// ****************************************************************************
// ----------------------------------------------------------------------------
namespace aurostd {  // namespace aurostd
  template<class utype> xmatrix<utype>                               // exp xmatrix
  exp_old(const xmatrix<utype>& a) {
    if(!a.issquare) {cerr << _AUROSTD_XLIBS_ERROR_ << "failure in xmatrix exp defined for square xmatrixes\n" << endl;exit(0);}
    xmatrix<utype> out(a.urows,a.ucols,a.lrows,a.lcols),an(a.urows,a.ucols,a.lrows,a.lcols);
    // UNUSED   bool convergence=FALSE;
    for(int n=0;n<=1000;n++) {
      if(n==0) identity(an);
      else an=(an*a)/utype(n);
      out=out+an;
      //    if(abs(trace(an)/trace(out))<_exponential_convergence)
    }
    return out;
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype> xmatrix<utype>                               // exp xmatrix
  exp(const xmatrix<utype>& a) {
#ifdef _XMATH_DEBUG_FUNCTIONS
    printf("M -> function exponential xmatrix: ");
    printf("a.lrows=%i, a.urows=%i, ",a.lrows,a.urows);
    printf("a.lcols=%i, a.ucols=%i\n",a.lcols,a.ucols);
#endif
    if(!a.issquare) {cerr << _AUROSTD_XLIBS_ERROR_ << "failure in xmatrix exp defined for square xmatrixes\n" << endl;exit(0);}
    xmatrix<utype> out(a.urows,a.ucols,a.lrows,a.lcols),an(a.urows,a.ucols,a.lrows,a.lcols);
    bool convergence=FALSE;
    for(int n=0;!convergence;n++) {
      if(n==0) identity(an);
      else an=(an*a)/utype(n);
      out=out+an;
      // cerr << n << endl;
      // if(abs(trace(an)/trace(out))<_exponential_convergence)
      if(n>30) if(abs(trace(an))<_exponential_convergence) convergence=TRUE;
      if(n>100) convergence=TRUE;
    }
    return out;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {  // namespace aurostd
  template<class utype> xmatrix<utype>                               // sin xmatrix
  sin(const xmatrix<utype>& a) {
#ifdef _XMATH_DEBUG_FUNCTIONS
    printf("M -> function sin xmatrix: ");
    printf("a.lrows=%i, a.urows=%i, ",a.lrows,a.urows);
    printf("a.lcols=%i, a.ucols=%i\n",a.lcols,a.ucols);
#endif
    if(!a.issquare) {cerr << _AUROSTD_XLIBS_ERROR_ << "failure in xmatrix sin defined for square xmatrixes\n" << endl;exit(0);}
    xmatrix<utype> out(a.urows,a.ucols,a.lrows,a.lcols), an(a.urows,a.ucols,a.lrows,a.lcols);
    bool convergence=FALSE;
    for(int n=0;!convergence;n++) {
      if(n==0) an=a;
      else an=(a*a*an)/utype(-1.0*(2.0*n+1.0)*(2.0*n));
      out=out+an;
      //    if(abs(trace(an)/trace(out))<_exponential_convergence)
      if(n>30) if(abs(trace(an))<_exponential_convergence) convergence=TRUE;
      if(n>100) convergence=TRUE;
    }
    return out;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {  // namespace aurostd
  template<class utype> xmatrix<utype>                               // cos xmatrix
  cos(const xmatrix<utype>& a) {
#ifdef _XMATH_DEBUG_FUNCTIONS
    printf("M -> function cos xmatrix: ");
    printf("a.lrows=%i, a.urows=%i, ",a.lrows,a.urows);
    printf("a.lcols=%i, a.ucols=%i\n",a.lcols,a.ucols);
#endif
    if(!a.issquare) {cerr << _AUROSTD_XLIBS_ERROR_ << "failure in xmatrix cos defined for square xmatrixes\n" << endl;exit(0);}
    xmatrix<utype> out(a.urows,a.ucols,a.lrows,a.lcols),an(a.urows,a.ucols,a.lrows,a.lcols);
    bool convergence=FALSE;
    for(int n=0;!convergence;n++) {
      if(n==0) identity(an);
      else an=(a*a*an)/utype(-1.0*(2.0*n)*(2.0*n-1.0));
      out=out+an;
      if(n>30) if(abs(trace(an))<_exponential_convergence) convergence=TRUE;
      if(n>100) convergence=TRUE;
    }
    return out;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {  // namespace aurostd
  template<class utype> xmatrix<utype>                              // sinh xmatrix
  sinh(const xmatrix<utype>& a)
  {
#ifdef _XMATH_DEBUG_FUNCTIONS
    printf("M -> function sinh xmatrix: ");
    printf("a.lrows=%i, a.urows=%i, ",a.lrows,a.urows);
    printf("a.lcols=%i, a.ucols=%i\n",a.lcols,a.ucols);
#endif
    if(!a.issquare) {cerr << _AUROSTD_XLIBS_ERROR_ << "failure in xmatrix sinh defined for square xmatrixes\n" << endl;exit(0);}
    xmatrix<utype> out(a.urows,a.ucols,a.lrows,a.lcols),an(a.urows,a.ucols,a.lrows,a.lcols);
    bool convergence=FALSE;
    for(int n=0;!convergence;n++) {
      if(n==0) an=a;
      else an=(a*a*an)/utype((2.0*n+1.0)*(2.0*n));
      out=out+an;
      // if(abs(trace(an)/trace(out))<_exponential_convergence)
      if(n>30) if(abs(trace(an))<_exponential_convergence) convergence=TRUE;
      if(n>100) convergence=TRUE;
    }
    return out;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {  // namespace aurostd
  template<class utype> xmatrix<utype>                              // cosh xmatrix
  cosh(const xmatrix<utype>& a)
  {
#ifdef _XMATH_DEBUG_FUNCTIONS
    printf("M -> function cosh xmatrix: ");
    printf("a.lrows=%i, a.urows=%i, ",a.lrows,a.urows);
    printf("a.lcols=%i, a.ucols=%i\n",a.lcols,a.ucols);
#endif
    if(!a.issquare) {cerr << _AUROSTD_XLIBS_ERROR_ << "failure in xmatrix cosh defined for square xmatrixes\n" << endl;exit(0);}
    xmatrix<utype> out(a.urows,a.ucols,a.lrows,a.lcols),an(a.urows,a.ucols,a.lrows,a.lcols);
    bool convergence=FALSE;
    for(int n=0;!convergence;n++) {
      if(n==0) identity(an);
      else an=(a*a*an)/utype((2.0*n)*(2.0*n-1.0));
      out=out+an;
      // if(abs(trace(an)/trace(out))<_exponential_convergence)
      if(n>30) if(abs(trace(an))<_exponential_convergence) convergence=TRUE;
      if(n>100) convergence=TRUE;
    }
    return out;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {  // namespace aurostd
  template<class utype>                                    // GaussJordan xmatrix
  void GaussJordan(xmatrix<utype>& A, xmatrix<utype>& B) {
    /// This function uses Gaussian Jordan elimination to solve A*x=b.  It returns the solution x and the inverse of A.
    if(A.lrows!=1) {cerr << _AUROSTD_XLIBS_ERROR_ << "ERROR GaussJordan [1] A.lrows!=1 <<  A.lrows=" << A.lrows << endl;exit(0);}  
    if(A.lcols!=1) {cerr << _AUROSTD_XLIBS_ERROR_ << "ERROR GaussJordan [2] A.lcols!=1 <<  A.lcols=" << A.lcols << endl;exit(0);}  
    if(B.lrows!=1) {cerr << _AUROSTD_XLIBS_ERROR_ << "ERROR GaussJordan [3] B.lrows!=1 <<  B.lrows=" << B.lrows << endl;exit(0);}  
    if(B.lcols!=1) {cerr << _AUROSTD_XLIBS_ERROR_ << "ERROR GaussJordan [4] B.lcols!=1 <<  B.lcols=" << B.lcols << endl;exit(0);}  
    if(A.urows!=A.ucols) {cerr << _AUROSTD_XLIBS_ERROR_ << "ERROR GaussJordan [5] A.urows!=A.ucols <<  A.urows=" << A.urows << " A.ucols=" << A.ucols << endl;exit(0);}  
    if(A.ucols!=B.urows) {cerr << _AUROSTD_XLIBS_ERROR_ << "ERROR GaussJordan [6] A.ucols!=B.urows <<  A.ucols=" << A.ucols << " B.urows=" << B.urows << endl;exit(0);}  
    int n=A.urows;
    int m=B.ucols;

    // cerr << "GaussJordan" << A.urows << " " << A.ucols << endl;
  
    int i,icol=1,irow=1,j,k,l,ll;
    utype big,dum,pivinv,temp;
  
    xvector<int> indxc(n),indxr(n),ipiv(n);

    for(j=1;j<=n;j++) ipiv[j]=0;
    for(i=1;i<=n;i++) {
      big=0.0;
      for(j=1;j<=n;j++)
	if(ipiv[j]!=1)
	  for(k=1;k<=n;k++) {
	    if(ipiv[k] == 0) {
	      if(aurostd::abs(A[j][k])>=big) {
		big=aurostd::abs(A[j][k]);
		irow=j;
		icol=k;
	      }
	    } else if(ipiv[k]>1) {cerr << _AUROSTD_XLIBS_ERROR_ << "ERROR GaussJordan [7]: Singular Matrix-1" << endl;exit(0);}
	  }
      ++(ipiv[icol]);
      if(irow!=icol) {
	for(l=1;l<=n;l++) SWAP(A[irow][l],A[icol][l]);
	for(l=1;l<=m;l++) SWAP(B[irow][l],B[icol][l]);
      }
      indxr[i]=irow;
      indxc[i]=icol;
      if(A[icol][icol]==(double) 0.0) {cerr << _AUROSTD_XLIBS_ERROR_ << "ERROR GaussJordan [8]: Singular Matrix-2" << endl;exit(0);}
      pivinv=1.0/A[icol][icol];
      A[icol][icol]=1.0;
      for(l=1;l<=n;l++) A[icol][l]*=pivinv;
      for(l=1;l<=m;l++) B[icol][l]*=pivinv;
      for(ll=1;ll<=n;ll++)
	if(ll!=icol) {
	  dum=A[ll][icol];
	  A[ll][icol]=0.0;
	  for(l=1;l<=n;l++) A[ll][l]-=A[icol][l]*dum;
	  for(l=1;l<=m;l++) B[ll][l]-=B[icol][l]*dum;
	}
    }
    for(l=n;l>=1;l--) {
      if(indxr[l]!=indxc[l])
	for(k=1;k<=n;k++)
	  SWAP(A[k][indxr[l]],A[k][indxc[l]]);
    }
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {   // least square stuff aurostd adaptation of nrecipes    // 1 August 2014
  // namespace aurostd   
  template<class utype> void gaussj(xmatrix<utype>& a, int n, xmatrix<utype>& b, int m) {  // with indices
    // linear equation solution by gauss-jordan elimination, a[1,n][1,n] is the input matrix.
    // b[1,n][1,m] is input containing the m right-hand side vectors. On the output a is replaced
    // by its matrix inverse, and b is replaced by the corresponding set of solution vectors.
    int i,icol,irow,j,k,l,ll;
    utype big,dum,pivinv,temp;
    
    // int n=a.rows,m=b.cols;
    // if(a.cols!=a.rows) {cerr << _AUROSTD_XLIBS_ERROR_ << "gaussj: a.cols!=a.rows" << endl; exit(0);}
    // if(b.rows!=a.rows) {cerr << _AUROSTD_XLIBS_ERROR_ << "gaussj: b.rows!=a.rows" << endl; exit(0);}
    
    if(n>a.rows) {cerr << _AUROSTD_XLIBS_ERROR_ << "gaussj: n>a.rows" << endl; exit(0);}
    if(n>b.rows) {cerr << _AUROSTD_XLIBS_ERROR_ << "gaussj: n>b.rows" << endl; exit(0);}
    if(m>b.cols) {cerr << _AUROSTD_XLIBS_ERROR_ << "gaussj: m>b.cols" << endl; exit(0);}
    
    xvector<int> indxc(1,n);
    xvector<int> indxr(1,n);
    xvector<int> ipiv(1,n);
    for (j=1;j<=(int) n;j++) ipiv[j]=0;
    for (i=1;i<=(int) n;i++) {
      big=0.0;
      for (j=1;j<=n;j++)
        if(ipiv[j] != 1)
          for (k=1;k<=n;k++) {
            if(ipiv[k] == 0) {
              if(aurostd::abs(a[j][k]) >= big) {
                big=aurostd::abs(a[j][k]);
                irow=j;
                icol=k;
              }
            } else if(ipiv[k] > 1) { cerr << _AUROSTD_XLIBS_ERROR_ << "gaussj: Singular Matrix-1" << endl; exit(0);}
          }
      ++(ipiv[icol]);
      if(irow != icol) {
        for (l=1;l<=n;l++) {SWAP(a[irow][l],a[icol][l]);}
        for (l=1;l<=m;l++) {SWAP(b[irow][l],b[icol][l]);}
      }
      indxr[i]=irow;
      indxc[i]=icol;
      if(a[icol][icol] == 0.0) { cerr << _AUROSTD_XLIBS_ERROR_ << "gaussj: Singular Matrix-2" << endl; exit(0);}
      pivinv=1.0/a[icol][icol];
      a[icol][icol]=1.0;
      for (l=1;l<=n;l++) a[icol][l] *= pivinv;
      for (l=1;l<=m;l++) b[icol][l] *= pivinv;
      for (ll=1;ll<=n;ll++)
        if(ll != icol) {
          dum=a[ll][icol];
          a[ll][icol]=0.0;
          for (l=1;l<=n;l++) a[ll][l] -= a[icol][l]*dum;
          for (l=1;l<=m;l++) b[ll][l] -= b[icol][l]*dum;
        }
    }
    for (l=n;l>=1;l--) {
      if(indxr[l] != indxc[l])
        for (k=1;k<=n;k++) {
          SWAP(a[k][indxr[l]],a[k][indxc[l]]);
        }
    }
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {   // least square stuff aurostd adaptation of nrecipes    // 1 August 2014
  template<class utype> 
  void lfit(xvector<utype> x, xvector<utype> y, xvector<utype> sig, 
            xvector<utype>& a, xvector<int> ia, 
            xmatrix<utype>& covar, utype& chisq,
            void (*funcs)(utype, xvector<utype>&)) {
    // Given a set of data points x[1,ndat],y[1,ndat] with individual standar deviation sig[1,ndat], use chisq minimization to fit for some or all the coefficients
    // a[1,ma] of a function that depends linearly on a, y=sum_i a_i*afunc_i(x). The input array ia[1,ma] indicates by nonzero entries those componends of a
    // that should be fitted for, and by zero entries those components that should be held fiuxed at their input values. 
    // The prgram returns value for a[1,ma], chisq,  and the covariance atrix covar[1,ma][1,ma]. (Parameters held fixed will return zero covariances.). 
    // The user supplies a routine funcs(x,xvector<afunc>) that returns the ma basis funcions evaluated at x=X in the array afunc[1,ma]
    
    int ndat=x.rows;
    if(y.rows!=x.rows) {cerr << _AUROSTD_XLIBS_ERROR_ << "lfit: y.rows!=x.rows" << endl; exit(0);}
    if(sig.rows!=x.rows) {cerr << _AUROSTD_XLIBS_ERROR_ << "lfit: sig.rows!=x.rows" << endl; exit(0);}
    
    int ma=a.rows;
    if(ia.rows!=a.rows) {cerr << _AUROSTD_XLIBS_ERROR_ << "lfit: ia.rows!=a.rows" << endl; exit(0);}
    
    int i,j,k,l,m,mfit=0;
    utype ym,wt,sum,sig2i;

    aurostd::xmatrix<utype> beta(1,ma,1,1);
    aurostd::xvector<utype> afunc(1,ma);
    for (j=1;j<=(int) ma;j++)
      if(ia[j]) mfit++;
    if(mfit == 0) { cerr << _AUROSTD_XLIBS_ERROR_ << "lfit: no parameters to be fitted" << endl; exit(0);}
    for (j=1;j<=mfit;j++) {
      for (k=1;k<=mfit;k++) covar[j][k]=0.0;
      beta[j][1]=0.0;
    }
    for (i=1;i<=(int) ndat;i++) {
      (*funcs)(x[i],afunc);
      ym=y[i];
      if(mfit < (int) ma) {
        for (j=1;j<=(int) ma;j++)
          if(!ia[j]) ym -= a[j]*afunc[j];
      }
      sig2i=1.0/(sig[i]*sig[i]);
      for (j=0,l=1;l<=(int) ma;l++) {
        if(ia[l]) {
          wt=afunc[l]*sig2i;
          for (j++,k=0,m=1;m<=l;m++)
            if(ia[m]) covar[j][++k] += wt*afunc[m];
          beta[j][1] += ym*wt;
        }
      }
    }
    for (j=2;j<=mfit;j++)
      for (k=1;k<j;k++)
        covar[k][j]=covar[j][k];
    gaussj(covar,mfit,beta,1); // operate up to mfit
    for (j=0,l=1;l<=(int) ma;l++)
      if(ia[l]) a[l]=beta[++j][1];
    chisq=0.0;
    for (i=1;i<=(int) ndat;i++) {
      (*funcs)(x[i],afunc);
      for (sum=0.0,j=1;j<=(int) ma;j++) sum += a[j]*afunc[j];
      chisq += (((y[i]-sum)/sig[i])*((y[i]-sum)/sig[i]));
    }
    covsrt(covar,ia,mfit);
  }
  
}
// ----------------------------------------------------------------------------
namespace aurostd {   // least square stuff aurostd adaptation of nrecipes    // 1 August 2014
  template<class utype> void covsrt(xmatrix<utype>&covar, xvector<int> ia, int mfit) {
    // Expand in storage the covariance matrix covar[1,ma][1,ma], so as to take into account parameters
    // that are being fixed. (for the latter, return zero covariance.)
    int i,j,k;
    utype temp;
    
    int ma=covar.rows;
    if(covar.cols!=covar.rows) {cerr << _AUROSTD_XLIBS_ERROR_ << "covsrt: covar.cols!=covar.rows" << endl; exit(0);}
    
    for (i=mfit+1;i<=ma;i++) {
      for (j=1;j<=i;j++) {
        covar[i][j]=covar[j][i]=0.0;
      }
    }
    k=mfit;
    for (j=ma;j>=1;j--) {
      if(ia[j]) {
        for (i=1;i<=ma;i++) {SWAP(covar[i][k],covar[i][j]);}
        for (i=1;i<=ma;i++) {SWAP(covar[k][i],covar[j][i]);}
        k--;
      }
    }
  }
}

// ----------------------------------------------------------------------------

namespace aurostd {  // namespace aurostd
  template<class utype>                                      //  std::cout operator <<
  std::ostream& operator<< (std::ostream& buf,const xmatrix<utype>& x) {
    char buf2[80];
    string iobuf1,iobuf2,iobuf3,iobuf4,iobuf5;         // buffers
    utype xij=0;int i,j;                                             // buffer x[i]
    bool done=FALSE;
    if(_isfloat(xij)) {                                    // floating point mode
      if(!_iscomplex(xij)) {                                     // real numbers
	if(_size(xij)==sizeof(long double)) {                     // long double
	  iobuf1="%13.7lle";                                     // long double
	  iobuf2="long double  ";                                  // long double
	  iobuf3="    0.0      ";                                // long double
	  iobuf4="[%2i]   ";                                       // long double
	  iobuf5="     ";                                        // long double
	  done=TRUE;
	}	
	if(_size(xij)==sizeof(double)) {                               // double
	  iobuf1="%11.4le";                                           // double
	  iobuf2="double     ";                                         // double
	  iobuf3="   0.0     ";                                       // double
	  iobuf4="[%2i]  ";                                             // double
	  iobuf5="    ";                                              // double
	  done=TRUE;
	}	
	if(_size(xij)==sizeof(float)) {                                 // float
	  iobuf1=" % 2.4lf";                                             // float
	  iobuf2="float    ";                                            // float
	  iobuf3="  0.0000";                                             // float
	  iobuf4="[%2i] ";                                               // float
	  iobuf5="    ";
	  done=TRUE;
	}	
      } else {                                                // xcomplex numbers
	if(_size(xij)==sizeof(xcomplex<long double>)) {  // xcomplex<long double>
          // 	  iobuf1=" (% 13.7lle,% 13.7lle)";                // xcomplex<long double>
          // 	  iobuf2="xcomplex<long double>";                 // xcomplex<long double>
          // 	  iobuf3=" (     0.0      ,     0.0      )";      // xcomplex<long double>
          // 	  iobuf4="[%2i]   ";                              // xcomplex<long double>
          // 	  iobuf5="       ";                               // xcomplex<long double>
	  done=TRUE;
	}	
	if(_size(xij)==sizeof(xcomplex<double>)) {            // xcomplex<double>
          // 	  iobuf1=" (% 11.5le,% 11.5le)";                       // xcomplex<double>
          // 	  iobuf2="xcomplex<double>";                           // xcomplex<double>
          // 	  iobuf3=" (   0.0      ,   0.0      )";               // xcomplex<double>
          // 	  iobuf4="[%2i]  ";                                    // xcomplex<double>
          // 	  iobuf5="      ";                                     // xcomplex<double>
	  done=TRUE;
	}	
	if(_size(xij)==sizeof(xcomplex<float>)) {              // xcomplex<float>
          // 	  iobuf1=" (% 10.4le,% 10.4le)";                        // xcomplex<float>
          // 	  iobuf2="xcomplex<float>";                             // xcomplex<float>
          // 	  iobuf3=" (   0.0     ,   0.0     )";                  // xcomplex<float>
          // 	  iobuf4="[%2i] ";                                      // xcomplex<float>
          // 	  iobuf5="      ";                                      // xcomplex<float>
	  done=TRUE;
	}	
      }
    } else {                                                      // integer mode
      if(_size(xij)==sizeof(long int))  {                            // long int
	iobuf1="%11i";                                                // long int
	iobuf2="long int     ";                                       // long int
	iobuf3="        0 ";                                          // long int
	iobuf4="[%2i] ";                                              // long int
	iobuf5="     ";                                               // long int
	done=TRUE;
      }	
      if(_size(xij)==sizeof(int))  {                                      // int
	iobuf1="%11i";                                                     // int
	iobuf2="int          ";                                            // int
	iobuf3="        0 ";                                               // int
	iobuf4="[%2i] ";                                                   // int
	iobuf5="     ";                                                    // int
	done=TRUE;
      }	
      if(_size(xij)==sizeof(char))  {                                    // char
	iobuf1="%3d";                                                     // char
	iobuf2="char ";                                                   // char
	iobuf3="  0 ";                                                    // char
	iobuf4="[%2i] ";                                                  // char
	iobuf5="";                                                        // char
	done=TRUE;
      }	
    }

    //cerr << iobuf2 << endl;
    if(done==FALSE)
      {cerr << _AUROSTD_XLIBS_ERROR_ << "failure in operator <<: no data type available for user type\n" << endl;exit(0);}

#ifdef _XMATH_DEBUG_OUTPUT
    buf << iobuf2;
    for(j=x.lcols;j<=x.ucols;j++) {
      sprintf(buf1,iobuf4.c_str(),j);                                            // above
      buf << buf1 << iobuf5;
    }
    buf << endl ;
#endif
    for(i=x.lrows;i<=x.urows;i++) {
#ifdef _XMATH_DEBUG_OUTPUT
      sprintf(buf1,iobuf4.c_str(),i);                                             // near
      buf << buf1 ;
#endif
#ifdef _XMATH_LATGEN_AL_GULP
      buf << "Al core" ;
#endif  
      for(j=x.lcols;j<=x.ucols;j++) {
	xij=x[i][j];
	if(!_iscomplex(xij)) {
	  if(_isfloat(xij)) {
	    if(abs(xij)> (double) _xmatrix_epsilon) {
	      sprintf(buf2,iobuf1.c_str(),aurostd::_real(xij));
	    } else {
	      //	      sprintf(buf2,iobuf3);
	      sprintf(buf2,iobuf1.c_str(),aurostd::_real(xij));
	    }
	  } else {
	    if(aurostd::_real(xij)!=0) {
	      sprintf(buf2,iobuf1.c_str(),aurostd::_real(xij));
	    } else {
	      //  sprintf(buf2,iobuf3.c_str());
	      sprintf(buf2,iobuf1.c_str(),aurostd::_real(xij));
	    }
	  }
	  buf << string(buf2) << " ";
	} else {
	  //	  if(abs(xij)>  (float) _xmatrix_epsilon)
	  //    sprintf(buf2,iobuf1.c_str(),real(xij),imag(xij));
	  //  else
	  //    sprintf(buf2,iobuf3.c_str());
	  buf << xij << " ";  // problem of printing in xcomplex
	}
	//	if(j<x.ucols) buf;
      }
      if(i<x.urows) buf << endl;
    }
    // cerr << "[3]" << endl;
    return buf;
  }
}

//*****************************************************************************
// EIGENVECTORS EIGENVALUES STUFF

// ****************************************************
namespace aurostd {

#define NRANSI
#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);a[k][l]=h+s*(g-h*tau);
  template<class utype>
  int jacobi(const xmatrix<utype> &ain,xvector<utype> &d,xmatrix<utype> &v) {
    // Computes all eigenvalues and eigenvectors of a real symmetric xmatrix a[1..n][1..n].
    // On output, elements of a above the diagonal are destroyed. d[1..n] returns the eigenvalues of a.
    // v[1..n][1..n] is a matrix whose columns contain, on output, the normalized eigenvectors of
    // a. The function returns the number of Jacobi rotations that were required.
    
    int j,iq,ip,i,n,nrot=0;
    utype tresh,theta,tau,t,sm,s,h,g,c;
    xmatrix<utype> a(ain);
    n=a.rows;
    if(a.rows!=a.cols) {cerr << _AUROSTD_XLIBS_ERROR_ << "JACOBI: 'a' matrix not square  a.rows" << a.rows << " a.cols=" << a.cols << endl;exit(0);}
    if(v.rows!=v.cols) {cerr << _AUROSTD_XLIBS_ERROR_ << "JACOBI: 'v' matrix not square  v.rows" << v.rows << " v.cols=" << v.cols << endl;exit(0);}
    if(a.rows!=v.rows) {cerr << _AUROSTD_XLIBS_ERROR_ << "JACOBI: 'a' and 'v' matrices must have same size  a.rows" << a.rows << " v.rows=" << v.rows << endl;exit(0);}
    if(a.rows!=d.rows) {cerr << _AUROSTD_XLIBS_ERROR_ << "JACOBI: 'a' and 'd' objects must have same size  a.rows" << a.rows << " d.rows=" << d.rows << endl;exit(0);}
  
    xvector<utype> b(1,n);
    xvector<utype> z(1,n);
    for (ip=1;ip<=n;ip++) {
      for (iq=1;iq<=n;iq++) v[ip][iq]=0.0;
      v[ip][ip]=1.0;
    }
    for (ip=1;ip<=n;ip++) {
      b[ip]=d[ip]=a[ip][ip];
      z[ip]=0.0;
    }
    nrot=0;
    for (i=1;i<=50;i++) {
      sm=0.0;
      for (ip=1;ip<=n-1;ip++) {
	for (iq=ip+1;iq<=n;iq++)
	  sm += aurostd::abs(a[ip][iq]);
      }
      if(sm == 0.0) {
	//~z;~b;
	return nrot;
      }
      if(i < 4)
	tresh=0.2*sm/(n*n);
      else
	tresh=0.0;
      for (ip=1;ip<=n-1;ip++) {
	for (iq=ip+1;iq<=n;iq++) {
	  g=100.0*aurostd::abs(a[ip][iq]);
	  if(i > 4 && (utype)(aurostd::abs(d[ip])+g) == (utype)aurostd::abs(d[ip])
             && (utype)(aurostd::abs(d[iq])+g) == (utype)aurostd::abs(d[iq]))
	    a[ip][iq]=0.0;
	  else if(aurostd::abs(a[ip][iq]) > tresh) {
	    h=d[iq]-d[ip];
	    if((utype)(aurostd::abs(h)+g) == (utype)aurostd::abs(h)) {
	      t=(a[ip][iq])/h;
	    } else {
	      theta=0.5*h/(a[ip][iq]);
	      t=1.0/(aurostd::abs(theta)+aurostd::sqrt(1.0+theta*theta));
	      if(theta < 0.0) t = -t;
	    }
	    c=1.0/sqrt(1+t*t);
	    s=t*c;
	    tau=s/(1.0+c);
	    h=t*a[ip][iq];
	    z[ip] -= h;
	    z[iq] += h;
	    d[ip] -= h;
	    d[iq] += h;
	    a[ip][iq]=0.0;
	    for (j=1;j<=ip-1;j++) {ROTATE(a,j,ip,j,iq);}
	    for (j=ip+1;j<=iq-1;j++) {ROTATE(a,ip,j,j,iq);}
	    for (j=iq+1;j<=n;j++) {ROTATE(a,ip,j,iq,j);}
	    for (j=1;j<=n;j++) {ROTATE(v,j,ip,j,iq);}
	    ++(nrot);
	  }
	}
      }
      for (ip=1;ip<=n;ip++) {
	b[ip] += z[ip];
	d[ip]=b[ip];
	z[ip]=0.0;
      }
    }
    cerr << _AUROSTD_XLIBS_ERROR_ << "JACOBI: Too many iterations in routine jacobi" << endl;
    exit(0);
  }
#undef ROTATE  
}

// ****************************************************
namespace aurostd {
  template<class utype>
  void eigsrt(xvector<utype> &d,xmatrix<utype> &v) {
    // Given the eigenvalues d[1..n]and eigenvectors v[1..n][1..n] as output fromjacobi
    // or tqli,this routine sorts the eigenvalues into descending order, and rearranges
    // the columns of v correspondingly. The method is straight insertion and is N2 rather than NlogN;
    // but since you have just done an N3 procedure to get the eigenvalues, you can afford yourself
    // this little indulgence.
    int k,j,i,n;
    utype p;
    
    n=v.rows;
    if(v.rows!=v.cols) {cerr << _AUROSTD_XLIBS_ERROR_ << "EIGSRT: 'v' matrix not square  v.rows" << v.rows << " v.cols=" << v.cols << endl;exit(0);}
    if(v.rows!=d.rows) {cerr << _AUROSTD_XLIBS_ERROR_ << "EIGSRT: 'v' and 'd' objects must have same size  v.rows" << v.rows << " d.rows=" << d.rows << endl;exit(0);}
    
    for (i=1;i<n;i++) {
      p=d[k=i];
      for (j=i+1;j<=n;j++)
	if(d[j] >= p) p=d[k=j];
      if(k != i) {
	d[k]=d[i];
	d[i]=p;
	for (j=1;j<=n;j++) {
	  p=v[j][i];
	  v[j][i]=v[j][k];
	  v[j][k]=p;
	}
      }
    }
  }
}

// ****************************************************
// CO 171129
namespace aurostd {
  template<class utype>
  xmatrix<utype> generalHouseHolderQRDecomposition(xmatrix<utype>& mat,const utype& tol) {
  // mat is mxn, m>=n
  // output:  Q
  // mat will change to R
  // See Numerical Linear Algebra, Trefethen and Bau, pg. 73

  string soliloquy="aurostd::generalHouseHolderQRDecomposition():";
  if(mat.rows<mat.cols){cerr << soliloquy << " ERROR! m<n, please flip the matrix." << endl;exit(1);}

  utype vModulus, xModulus;
  std::vector<xmatrix<utype> > V;

  for (uint i = 1; i < (uint)mat.cols + 1; i++) {
    // create x
    xmatrix<utype> x(mat.rows - i + 1, 1);
    for (uint j = i; j < (uint)mat.rows + 1; j++) {
      x(j - i + 1, 1) = mat(j, i);
    }
    // create identity vector
    xmatrix<utype> e1(x.rows, 1);  // automatically sets all to 0
    e1(1, 1) = (utype)1.0;
    xModulus = sqrt(sum(trasp(x) * x));
    xmatrix<utype> v(x.rows, 1);
    if(xModulus >= tol) {
      v = xModulus * e1;
      if(std::signbit(x(1, 1))) {
        // negative
        v = v * (utype)-1.0;
      }
      v += x;
      // modulus of v
      vModulus = sqrt(sum(trasp(v) * v));
      v = v / vModulus;
      xmatrix<utype> A(mat.rows - i + 1, mat.cols - i + 1);
      for (uint j = i; j < (uint)mat.rows + 1; j++) {
        for (uint k = i; k < (uint)mat.cols + 1; k++) {
          A(j - i + 1, k - i + 1) = mat(j, k);
        }
      }
      A = A - (utype)2.0 * v * trasp(v) * A;
      A = A;
      for (uint j = i; j < (uint)mat.rows + 1; j++) {
        for (uint k = i; k < (uint)mat.cols + 1; k++) {
          mat(j, k) = A(j - i + 1, k - i + 1);
        }
      }
    }
    V.push_back(v);
  }

  xmatrix<utype> Q(mat.rows, mat.rows);
  xmatrix<utype> ek(mat.rows, 1);
  for (uint k = 1; k < (uint)mat.rows + 1; k++) {
    // create identity vector
    for (uint i = 1; i < (uint)mat.rows + 1; i++) {
      if(i == k) {
        ek(i, 1) = (utype)1;
      } else {
        ek(i, 1) = (utype)0;
      }
    }
    for (uint i = (uint)mat.cols; i > 0; i--) {
      xmatrix<utype> x(mat.rows - i + 1, 1);
      for (uint j = i; j < (uint)mat.rows + 1; j++) {
        x(j - i + 1, 1) = ek(j, 1);
      }
      x = x - (utype)2.0 * V.at(i - 1) * trasp(V.at(i - 1)) * x;
      for (uint j = i; j < (uint)mat.rows + 1; j++) {
        ek(j, 1) = x(j - i + 1, 1);
      }
    }
    for (uint i = 1; i < (uint)mat.rows + 1; i++) {
      Q(i, k) = ek(i, 1);
    }
  }
  return Q;
  }
}

// ****************************************************
namespace aurostd {
  template<class utype>
  void tred2(const xmatrix<utype> &a,xvector<utype> &d,xvector<utype> &e) {
    // Householder reduction of a real, symmetric matrix a[1..n][1..n].
    // On output, a is replaced by the orthogonal matrix Q eﬀecting the
    // transformation. d[1..n] returns the diagonal elments of
    // the tridiagonal matrix, and e[1..n] the oﬀ-diagonal elements, with e[1]=0.

    int l,k,j,i,n;
    utype scale,hh,h,g,f;
    
    n=a.rows;
    if(a.rows!=a.cols) {cerr << _AUROSTD_XLIBS_ERROR_ << "TRED2: 'a' matrix not square  a.rows" << a.rows << " a.cols=" << a.cols << endl;exit(0);}
    if(a.rows!=d.rows) {cerr << _AUROSTD_XLIBS_ERROR_ << "TRED2: 'a' and 'd' objects must have same size  a.rows" << a.rows << " d.rows=" << d.rows << endl;exit(0);}
    if(a.rows!=e.rows) {cerr << _AUROSTD_XLIBS_ERROR_ << "TRED2: 'a' and 'e' objects must have same size  a.rows" << a.rows << " e.rows=" << e.rows << endl;exit(0);}

    for (i=n;i>=2;i--) {
      l=i-1;
      h=scale=0.0;
      if(l > 1) {
	for (k=1;k<=l;k++)
	  scale += aurostd::abs(a[i][k]);
	if(scale == 0.0)
	  e[i]=a[i][l];
	else {
	  for (k=1;k<=l;k++) {
	    a[i][k] /= scale;
	    h += a[i][k]*a[i][k];
	  }
	  f=a[i][l];
	  g=(f >= 0.0 ? -sqrt(h) : sqrt(h));
	  e[i]=scale*g;
	  h -= f*g;
	  a[i][l]=f-g;
	  f=0.0;
	  for (j=1;j<=l;j++) {
	    a[j][i]=a[i][j]/h;
	    g=0.0;
	    for (k=1;k<=j;k++)
	      g += a[j][k]*a[i][k];
	    for (k=j+1;k<=l;k++)
	      g += a[k][j]*a[i][k];
	    e[j]=g/h;
	    f += e[j]*a[i][j];
	  }
	  hh=f/(h+h);
	  for (j=1;j<=l;j++) {
	    f=a[i][j];
	    e[j]=g=e[j]-hh*f;
	    for (k=1;k<=j;k++)
	      a[j][k] -= (f*e[k]+g*a[i][k]);
	  }
	}
      } else
	e[i]=a[i][l];
      d[i]=h;
    }
    d[1]=0.0;
    e[1]=0.0;
    // Contents of this loop can be omitted if eigenvectors not
    // wanted except for statement d[i]=a[i][i];
    for (i=1;i<=n;i++) {
      l=i-1;
      if(d[i]) {
	for (j=1;j<=l;j++) {
	  g=0.0;
	  for (k=1;k<=l;k++)
	    g += a[i][k]*a[k][j];
	  for (k=1;k<=l;k++)
	    a[k][j] -= g*a[k][i];
	}
      }
      d[i]=a[i][i];
      a[i][i]=1.0;
      for (j=1;j<=l;j++) a[j][i]=a[i][j]=0.0;
    }
  }

}

// ****************************************************
namespace aurostd {

  template<class utype>
  utype NR_SQR(utype a) {
    if(a==(utype) 0.0) return 0.0; else return a*a;
  }
  
  template<class utype>
  utype pythag(utype a, utype b) {
    utype absa,absb;
    absa=aurostd::abs(a);
    absb=aurostd::abs(b);
    if(absa > absb) return absa*sqrt(1.0+NR_SQR(absb/absa));
    else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+NR_SQR(absa/absb)));
  }
  
  
#define NR_SIGN(a,b) ((b) >= 0.0 ? aurostd::abs(a) : -aurostd::abs(a))
  template<class utype>
  void tqli(xvector<utype> &d,xvector<utype> &e,xmatrix<utype> &z) {
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
    int m,l,iter,i,k,n;
    utype s,r,p,g,f,dd,c,b;
    
    n=z.rows;
    if(z.rows!=z.cols) {cerr << _AUROSTD_XLIBS_ERROR_ << "TQLI: 'z' matrix not square  z.rows" << z.rows << " z.cols=" << z.cols << endl;exit(0);}
    if(z.rows!=d.rows) {cerr << _AUROSTD_XLIBS_ERROR_ << "TQLI: 'z' and 'd' objects must have same size  z.rows" << z.rows << " d.rows=" << d.rows << endl;exit(0);}
    if(z.rows!=e.rows) {cerr << _AUROSTD_XLIBS_ERROR_ << "TQLI: 'z' and 'e' objects must have same size  z.rows" << z.rows << " e.rows=" << e.rows << endl;exit(0);}
  
    for (i=2;i<=n;i++) e[i-1]=e[i];
    e[n]=0.0;
    for (l=1;l<=n;l++) {
      iter=0;
      do {
	for (m=l;m<=n-1;m++) {
	  dd=aurostd::abs(d[m])+aurostd::abs(d[m+1]);
	  if((utype)(aurostd::abs(e[m])+dd) == dd) break;
	}
	if(m != l) {
	  if(iter++ == 30) {cerr << _AUROSTD_XLIBS_ERROR_ << "Too many iterations in tqli" << endl;exit(0);}
	  g=(d[l+1]-d[l])/(2.0*e[l]);
	  r=pythag(g,(utype) 1.0);
	  g=d[m]-d[l]+e[l]/(g+NR_SIGN(r,g));
	  s=c=1.0;
	  p=0.0;
	  for (i=m-1;i>=l;i--) {
	    f=s*e[i];
	    b=c*e[i];
	    e[i+1]=(r=pythag(f,g));
	    if(r == 0.0) {
	      d[i+1] -= p;
	      e[m]=0.0;
	      break;
	    }
	    s=f/r;
	    c=g/r;
	    g=d[i+1]-p;
	    r=(d[i]-g)*s+2.0*c*b;
	    d[i+1]=g+(p=s*r);
	    g=c*r-b;
	    for (k=1;k<=n;k++) {
	      f=z[k][i+1];
	      z[k][i+1]=s*z[k][i]+c*f;
	      z[k][i]=c*z[k][i]-s*f;
	    }
	  }
	  if(r == 0.0 && i >= l) continue;
	  d[l] -= p;
	  e[l]=g;
	  e[m]=0.0;
	}
      } while (m != l);
    }
  }
#undef NR_SIGN
}

// ****************************************************
namespace aurostd {
#define RADIX 2.0
  template<class utype>
  void balanc(xmatrix<utype> &a) {
    // Given a matrix a[1..n][1..n], this routine replaces it by
    // a balanced matrix with i dentical eigenvalues. A symmetric matrix
    // is already balanced and is unaﬀected by this procedure. The
    // parameter RADIX should be the machine’s ﬂoating-point radix.
    int last,j,i,n;
    utype s,r,g,f,c,sqrdx;
    
    n=a.rows;
    if(a.rows!=a.cols) {cerr << _AUROSTD_XLIBS_ERROR_ << "BALANCE: 'a' matrix not square  a.rows" << a.rows << " a.cols=" << a.cols << endl;exit(0);}
    
    sqrdx=RADIX*RADIX;
    last=0;
    while (last == 0) {
      last=1;
      for (i=1;i<=n;i++) {
	r=c=0.0;
	for (j=1;j<=n;j++)
	  if(j != i) {
	    c += aurostd::abs(a[j][i]);
	    r += aurostd::abs(a[i][j]);
	  }
	if(c && r) {
	  g=r/RADIX;
	  f=1.0;
	  s=c+r;
	  while (c<g) {
	    f *= RADIX;
	    c *= sqrdx;
	  }
	  g=r*RADIX;
	  while (c>g) {
	    f /= RADIX;
	    c /= sqrdx;
	  }
	  if((c+r)/f < 0.95*s) {
	    last=0;
	    g=1.0/f;
	    for (j=1;j<=n;j++) a[i][j] *= g;
	    for (j=1;j<=n;j++) a[j][i] *= f;
	  }
	}
      }
    }
  }
}
#undef RADIX

// ****************************************************
namespace aurostd {
#define ELMHES_SWAP(g,h) {y=(g);(g)=(h);(h)=y;}
  template<class utype>
  // Reduction to Hessenberg form by the elimination method.
  // The real, nonsymmetric matrix a[1..n][1..n] is replaced by an upper
  // Hessenberg matrix with identical eigenvalues.
  // Recommended, but not required, is that this routine be preceded
  // by balanc. On output, the Hessenberg matrix is in elements a[i][j] with i<=j+1.
  // Elements with i>j+1 are to be thought of as zero,
  // but are returned with random values.
  void elmhes(xmatrix<utype> &a) {
    int m,j,i,n;
    utype y,x;
    
    n=a.rows;
    if(a.rows!=a.cols) {cerr << _AUROSTD_XLIBS_ERROR_ << "ELMHES: 'a' matrix not square  a.rows" << a.rows << " a.cols=" << a.cols << endl;exit(0);}
    
    for (m=2;m<n;m++) {
      x=0.0;
      i=m;
      for (j=m;j<=n;j++) {
	if(aurostd::abs(a[j][m-1]) > aurostd::abs(x)) {
	  x=a[j][m-1];
	  i=j;
	}
      }
      if(i != m) {
	for (j=m-1;j<=n;j++) ELMHES_SWAP(a[i][j],a[m][j]);
	for (j=1;j<=n;j++) ELMHES_SWAP(a[j][i],a[j][m]);
      }
      if(x) {
	for (i=m+1;i<=n;i++) {
	  if((y=a[i][m-1]) != 0.0) {
	    y /= x;
	    a[i][m-1]=y;
	    for (j=m;j<=n;j++)
	      a[i][j] -= y*a[m][j];
	    for (j=1;j<=n;j++)
	      a[j][m] += y*a[j][i];
	  }
	}
      }
    }
  }
#undef ELMHES_SWAP
}

// ****************************************************
namespace aurostd {
  // Finds all eigenvalues of an upper Hessenberg matrix a[1..n][1..n].
  // On input a can be exactly as output from elmhes; on output it is destroyed.
  // The real and imaginary parts of the eigenvalues are returned in
  //wr[1..n] and wi[1..n], respectively.

#define NR_SIGN(a,b) ((b) >= 0.0 ? aurostd::abs(a) : -aurostd::abs(a))

  template<class utype>
  void hqr(xmatrix<utype> &a,xvector<utype> &wr,xvector<utype> &wi) {
    int nn,m,l,k,j,its,i,mmin,n;
    utype z,y,x,w,v,u,t,s,r=0,q=0,p=0,anorm;

    n=a.rows;
    if(a.rows!=a.cols) {cerr << _AUROSTD_XLIBS_ERROR_ << "HQR: 'a' matrix not square  a.rows" << a.rows << " a.cols=" << a.cols << endl;exit(0);}

    anorm=aurostd::abs(a[1][1]);
    for (i=2;i<=n;i++)
      for (j=(i-1);j<=n;j++)
	anorm += aurostd::abs(a[i][j]);
    nn=n;
    t=0.0;
    while (nn >= 1) {
      its=0;
      do {
	for (l=nn;l>=2;l--) {
	  s=aurostd::abs(a[l-1][l-1])+aurostd::abs(a[l][l]);
	  if(s == 0.0) s=anorm;
	  if((utype)(aurostd::abs(a[l][l-1]) + s) == s) break;
	}
	x=a[nn][nn];
	if(l == nn) {
	  wr[nn]=x+t;
	  wi[nn--]=0.0;
	} else {
	  y=a[nn-1][nn-1];
	  w=a[nn][nn-1]*a[nn-1][nn];
	  if(l == (nn-1)) {
	    p=0.5*(y-x);
	    q=p*p+w;
	    z=sqrt(aurostd::abs(q));
	    x += t;
	    if(q >= 0.0) {
	      z=p+NR_SIGN(z,p);
	      wr[nn-1]=wr[nn]=x+z;
	      if(z) wr[nn]=x-w/z;
	      wi[nn-1]=wi[nn]=0.0;
	    } else {
	      wr[nn-1]=wr[nn]=x+p;
	      wi[nn-1]= -(wi[nn]=z);
	    }
	    nn -= 2;
	  } else {
	    if(its == 30) {cerr << _AUROSTD_XLIBS_ERROR_ << "HQR: Too many iterations in hqr" << endl;exit(0);}
	    if(its == 10 || its == 20) {
	      t += x;
	      for (i=1;i<=nn;i++) a[i][i] -= x;
	      s=aurostd::abs(a[nn][nn-1])+aurostd::abs(a[nn-1][nn-2]);
	      y=x=0.75*s;
	      w = -0.4375*s*s;
	    }
	    ++its;
	    for (m=(nn-2);m>=l;m--) {
	      z=a[m][m];
	      r=x-z;
	      s=y-z;
	      p=(r*s-w)/a[m+1][m]+a[m][m+1];
	      q=a[m+1][m+1]-z-r-s;
	      r=a[m+2][m+1];
	      s=aurostd::abs(p)+aurostd::abs(q)+aurostd::abs(r);
	      p /= s;
	      q /= s;
	      r /= s;
	      if(m == l) break;
	      u=aurostd::abs(a[m][m-1])*(aurostd::abs(q)+aurostd::abs(r));
	      v=aurostd::abs(p)*(aurostd::abs(a[m-1][m-1])
				 +aurostd::abs(z)+aurostd::abs(a[m+1][m+1]));
	      if((utype)(u+v) == v) break;
	    }
	    for (i=m+2;i<=nn;i++) {
	      a[i][i-2]=0.0;
	      if(i != (m+2)) a[i][i-3]=0.0;
	    }
	    for (k=m;k<=nn-1;k++) {
	      if(k != m) {
		p=a[k][k-1];
		q=a[k+1][k-1];
		r=0.0;
		if(k != (nn-1)) r=a[k+2][k-1];
		if((x=aurostd::abs(p)+aurostd::abs(q)+aurostd::abs(r)) != 0.0) {
		  p /= x;
		  q /= x;
		  r /= x;
		}
	      }
	      if((s=NR_SIGN(sqrt(p*p+q*q+r*r),p)) != 0.0) {
		if(k == m) {
		  if(l != m)
		    a[k][k-1] = -a[k][k-1];
		} else
		  a[k][k-1] = -s*x;
		p += s;
		x=p/s;
		y=q/s;
		z=r/s;
		q /= p;
		r /= p;
		for (j=k;j<=nn;j++) {
		  p=a[k][j]+q*a[k+1][j];
		  if(k != (nn-1)) {
		    p += r*a[k+2][j];
		    a[k+2][j] -= p*z;
		  }
		  a[k+1][j] -= p*y;
		  a[k][j] -= p*x;
		}
		mmin = nn<k+3 ? nn : k+3;
		for (i=l;i<=mmin;i++) {
		  p=x*a[i][k]+y*a[i][k+1];
		  if(k != (nn-1)) {
		    p += z*a[i][k+2];
		    a[i][k+2] -= p*r;
		  }
		  a[i][k+1] -= p*q;
		  a[i][k] -= p;
		}
	      }
	    }
	  }
	}
      } while (l < nn-1);
    }
  }
#undef NR_SIGN
}

// ****************************************************
namespace aurostd {
  // Finds all eigenvalues of matrix a[1..n][1..n]. The real and imaginary parts
  // of the eigenvalues are returned in wr[1..n] and wi[1..n], respectively.

  template<class utype>
  void eigen(const xmatrix<utype> &ain,xvector<utype> &wr,xvector<utype> &wi) {
    xmatrix<utype> a(ain);
    balanc(a);
    elmhes(a);
    hqr(a,wr,wi);
  }
}

//*****************************************************************************
// ---------------------------------------------------------- aurostd::cematrix
namespace aurostd { // namespace aurostd
#define cematrix_EXIT_RANK_NOT_MATCH 56
#define cematrix_EQUAL_DOUBLE 1.0e-9 // two doubles are equal if difference is smaller than it

  cematrix::cematrix() { // default constructor
    nrow=1;
    ncol=1;
    M=xmatrix<double>(1,1,1,1);
    W=xvector<double>(1,1);
    U=xmatrix<double>(1,1,1,1);
    V=xmatrix<double>(1,1,1,1);
    a_vec=xvector<double>(1,1);
    a_nvec.clear();
    chisq=0.0;
    Cov=xmatrix<double>(1,1,1,1);
  }
 
  cematrix::cematrix(const xmatrix<double> & A_in) { // copy constructor
    nrow=A_in.rows;
    ncol=A_in.cols;
    M=xmatrix<double>(1,1,nrow,ncol);
    for(int i=1;i<=nrow;i++) 
      for(int j=1;j<=ncol;j++)
	M[i][j]=A_in[i][j];
    W=xvector<double>(1,ncol);
    V=xmatrix<double>(1,1,ncol,ncol);
    U=xmatrix<double>(1,1,nrow,ncol);
    a_vec=xvector<double>(1,1);
    a_nvec.clear();
    chisq=0.0;
    Cov=xmatrix<double>(1,1,ncol,ncol);
  }
 
  cematrix::~cematrix() { // default deconstructor
    a_nvec.clear();
  }

  void cematrix::LeastSquare(xvector<double>& y_vec, xvector<double>& y_sigma) { // function
    if(nrow !=y_vec.rows ) {
      cerr << "ERROR - cematrix::LeastSquare: no match of ranks of b and A " << endl;
      cerr << "ERROR - cematrix::LeastSquare: LeastSquares: input two matrices A (m x n) and b (m x 1) " << endl;
      exit(cematrix_EXIT_RANK_NOT_MATCH);
    }
    //SVDcmp(A);
    SVDFit(y_vec, y_sigma);
  }

  double cematrix::Pythag2(double a, double b) { // calculate (a^2+ b^2)^(1/2)
    // from dlapy2.f in Lapack
    double aabs=abs(a),babs=abs(b);
    double val_min,val_max;
    val_min=min(aabs,babs);
    val_max=max(aabs,babs);
    if(val_min ==0) {
      return val_max;
    } else {
      return val_max*sqrt(1.0+(val_min/val_max)*(val_min/val_max));
    }
  }

  void cematrix::SVDsolve(xvector<double>& b_vec) {
    // solve the least squares problem after SVDcmp()
    // it will use xmatrix operators later
    // the results are stored in a_vec
    xvector<double> temp(1,ncol),a_tmp(1,ncol);
    double wj,s;
    const double _NONZERO=1.0e-8;

    for(int j=1;j<=ncol;j++) {
      //wj=W.at(j-1);
      wj=W[j];
      //bj=b_vec[j];
      s=0.0;
      //if(wj !=0.0 ) {
      if(wj> _NONZERO) {
	// only if W[j] !=0
	for(int i=1;i <=nrow;i++) s+=U[i][j]*b_vec[i];
	s /=wj;
	temp[j]=s;
      }
    }

    for(int j=1;j<=ncol;j++) {  //multiply V
      s=0.0;
      for(int i=1;i<=ncol;i++)
	s+=V[j][i]*temp[i];
      a_tmp[j]=s;
    }

    a_vec=a_tmp;
    a_nvec.clear();
    for(int i=1;i<=ncol;i++) {
      a_nvec.push_back(a_vec[i]);
    }
  }

  //void cematrix::SVDcmp_NR()
  bool cematrix::SVDcmp_NR() {
    // SVD decompose matrix A=U Z V^T
    // decomposed matrices U Z V are stored
    xvector<double> rv1(1,ncol);
    double g, scale,anorm;
    int l,nm,jj,j,k,i;
    double f,c,h,x,y,z,s;
    bool flag;
    bool flag_convergence=true;
    // cerr << "ncol " << ncol << " nrow " << nrow << endl;
    // cerr << M << endl;
  
    l=0;
    g=0.0;
    scale=0.0;
    anorm=0.0;
    xvector<double> W_tmp(1,ncol);
    xmatrix<double> V_tmp(1,1,ncol,ncol);
    xmatrix<double> A(1,1,nrow,ncol);
    for(i=1;i<=nrow;i++) 
      for(j=1;j<=ncol;j++) 
	if(aurostd::abs(M(i,j))<cematrix_EQUAL_DOUBLE)
	  M(i,j)=0;
    A=M; // not destroy input matrix A
    
    // Householder reduction to bidiagonal form
    for(i=1;i<=ncol;i++) {
      l=i+1;
      rv1[i]=scale*g;
      g=0.0;
      s=0.0;
      scale=0.0;
      if(i <=nrow) {
	for(k=i;k <=nrow;k++ )
	  scale+=abs(A[k][i]);
	if(abs(scale)> cematrix_EQUAL_DOUBLE ) { // scale !=0
	  for(k=i;k<=nrow;k++) {
	    A[k][i] /=scale;
	    s+=A[k][i]*A[k][i];
	  }
	  f=A[i][i];
	  g=(-_sign(sqrt(s),f));
	  h=f*g - s;
	  A[i][i]=f - g;
	  for(j=l;j <=ncol;j++) {
	    for(s=0.0,k=i;k <=nrow;k++)
	      s+=A[k][i]*A[k][j];
	    f=s/h;
	    for(k=i;k<=nrow;k++)
	      A[k][j]+=f*A[k][i];
	  }
	  for(k=i;k<=nrow ;k++)
	    A[k][i] *=scale;
	}
      }
      W_tmp[i]=scale*g;
      g=0.0;
      s=0.0;
      scale=0.0;
      if(i <=nrow && i !=ncol ) {
	for(k=l;k<=ncol;k++) 
	  scale+=abs(A[i][k]);
	if(scale !=0.0 ) {
	  for(k=l;k<=ncol;k++) {
	    A[i][k] /=scale;
	    s+=A[i][k]*A[i][k];
	  }
	  f=A[i][l];
	  g=(-_sign(sqrt(s),f));
	  h=f*g - s;
	  A[i][l]=f - g;
	  for(k=l;k<=ncol;k++) 
	    rv1[k]=A[i][k]/h;
	  for(j=l;j<=nrow;j++) {
	    for(s=0.0,k=l;k <=ncol;k++)
	      s+=A[j][k]*A[i][k];
	    for(k=l;k<=ncol;k++)
	      A[j][k]+=s*rv1[k];
	  }
	  for(k=l;k<=ncol;k++)
	    A[i][k] *=scale;
	}
      }
      anorm=max(anorm,(abs(W_tmp[i])+ abs(rv1[i])) );
    } // i
    for(i=ncol;i>=1;i--) { //Accumulation of right-hand trnasformations
      if(i < ncol) {
	if(g !=0.0 ) {
	  for(j=l;j <=ncol;j++ ) 
	    V_tmp[j][i]=(A[i][j]/A[i][l])/g;
	  
	  for(j=l;j <=ncol;j++) {
	    for(s=0.0,k=l;k<=ncol;k++) 
	      s+=A[i][k]*V_tmp[k][j];
	    
	    for(k=l;k<=ncol;k++) 
	      V_tmp[k][j]+=s*V_tmp[k][i];
	    
	  }
	}
	for(j=l;j<=ncol;j++) {
	  V_tmp[i][j]=0.0;
	  V_tmp[j][i]=0.0;
	}
      }
      V_tmp[i][i]=1.0;
      g=rv1[i];
      l=i;
    }
    for(i=min(nrow,ncol);i>=1;i--) { // Accumulation of left-hand transformaions
      l=i+1;
      g=W_tmp[i];
      for(j=l;j<=ncol;j++) 
	A[i][j]=0.0;
      
      if(g !=0.0) {
	g=1.0/g;
	for(j=l;j<=ncol;j++) {
	  for(s=0.0,k=l;k<=nrow;k++) 
	    s+=A[k][i]*A[k][j];
	  f=(s/A[i][i])*g;
	  for(k=i;k<=nrow;k++) 
	    A[k][j]+=f*A[k][i];
	}
	for(j=i;j<=nrow;j++) 
	  A[j][i] *=g;
      } else {
	for(j=i;j<=nrow;j++) 
	  A[j][i]=0.0;
      }
      ++A[i][i];
    }
    for(k=ncol;k>=1;k--) {
      // Diagonalization of the bidiagonal form
      // Loop over singular values and over alowed iterations
      for(int its=1;its<=_MAX_ITS;its++) {
	flag=true;
	//for(l=k;l>=1;l-- ) { // test forsplitting
	// to avoid out of range as nm cannot be 0
	// should check!!!!!
	nm=1;
	for(l=k;l>=2;l-- ) { // test forsplitting
	  nm=l-1;//rv1[1] is always zero
	  if((abs(rv1[l])+ anorm)==anorm) {
	    flag=false;
	    break;
	  }
	  //if(nm !=0 ) {
	  // to avoid abort due to nm ==0
	  // xvector start with 1
	  if((abs(W_tmp[nm])+ anorm)==anorm ) break;
	  //}
	}
	if(flag ) {
	  // cancellation of rv1[l],if l> 1
	  c=0.0;
	  s=1.0;
	  for(i=l;i <=k;i++) {
	    f=s*rv1[i];
	    rv1[i]=c*rv1[i];
	    if((abs(f)+ anorm) ==anorm ) break;
	    g=W_tmp[i];
	    h=Pythag2(f,g);
	    W_tmp[i]=h;
	    h=1.0/h;
	    c=g*h;
	    s=-f*h;
	    for(j=1;j<=nrow;j++) {
	      y=A[j][nm];
	      z=A[j][i];
	      A[j][nm]=y*c+z*s;
	      A[j][i]=z*c-y*s;
	    }
	  }
	}
	z=W_tmp[k];
  
	//   cerr << "SVD its=" << its << " l=" << l << " k=" << endl;
	if(l ==k) { // convergence
	  //   cerr << "SVD its=" << its << " l=" << l << " k=" << k << endl;
	  if(z < 0.0 ) { // singular value is made nonnegative
	    W_tmp[k]=-z;
	    for(j=1;j<=ncol;j++) 
	      V_tmp[j][k]=-V_tmp[j][k];
	  }
	  break;
	}

	//   cerr << "SVD its=" << its << endl;
	if(its ==_MAX_ITS ) {
	  cerr << "ERROR - cematrix::SVDcmp_NR: Not converged in " << _MAX_ITS << " SVDcmp iterations !" << endl;;
	  cerr << "ERROR - cematrix::SVDcmp_NR: ncol " << ncol << " nrow " << nrow << endl;
	  cerr << M << endl;
	  //exit(_EXIT_NO_CONVERGENCE);
	  flag_convergence=false;
	}
	x=W_tmp[l];// shift from bottom 2x2 minor
	nm=k-1;
	//if(nm !=0 ) {
	y=W_tmp[nm];
	g=rv1[nm];
	//} else {
	//  y=1.0;
	//}
	h=rv1[k];
	f=((y-z)*(y+z)+ (g-h)*(g+h))/(2.0*h*y);
	g=Pythag2(f,1.0);
	f=((x-z)*(x+z)+h*((y/(f+_sign(g,f)))-h))/x;
	c=1.0;s=1.0;// next QR transformation
	for(j=l;j<=nm;j++) {
	  i=j+1;g=rv1[i];y=W_tmp[i];h=s*g;g=c*g;
	  z=Pythag2(f,h);rv1[j]=z;c=f/z;s=h/z;f=x*c+g*s;
	  g=g*c-x*s;h=y*s;y*=c;
	  for(jj=1;jj<=ncol;jj++) {
	    x=V_tmp[jj][j];
	    z=V_tmp[jj][i];
	    V_tmp[jj][j]=x*c+z*s;
	    V_tmp[jj][i]=z*c-x*s;
	  }
	  z=Pythag2(f,h);
	  W_tmp[j]=z;
	  if(z!=0.0) { // rotation can be arbitrary if z=0
	    z=1.0/z;c=f*z;s=h*z;
	  }
	  f=c*g+ s*y;
	  x=c*y-s*g;
	  for(int jj=1;jj<=nrow;jj++) {
	    y=A[jj][j];
	    z=A[jj][i];
	    A[jj][j]=y*c+z*s;
	    A[jj][i]=z*c-y*s;
	  }
	}
	rv1[l]=0.0;
	rv1[k]=f;
	W_tmp[k]=x;
      } // its
    } // k

    W=W_tmp;
    V=V_tmp;
    U=A;
    //// output U V W
    //
    //cerr.setf(ios_base::fixed);
    //cerr.precision(6);
    //cerr.width(12);
    ////cout.setw(12);
    //cerr << "*** Decomposition Matrice *** " << endl;;
    //cerr << "*** U matrix *** " << endl;;
    //for(i=1;i <=nrow;i++) { // U
    //  for(j=1;j <=ncol;j++) { cerr << setw(12) << U[i][j] << " "; } cerr << endl;
    //}
    //cerr << "*** W matrix diagonal elements*** " << endl;;
    //for(i=1;i <=ncol;i++) {cerr << setw(12) << W[i] << " ";}
    //cerr << endl;
    //cerr << "*** V matrix *** " << endl;;
    //for(i=1;i <=ncol;i++) { // U
    //  for(j=1;j <=ncol;j++) { cerr << setw(12) << V[i][j] << " ";} cerr << endl;
    //}
    //cerr << "Check the produce against the original matrix " << endl;;
    //cerr << "Original matrix " << endl;;
    //for(i=1;i <=nrow;i++) { // U
    //  for(j=1;j <=ncol;j++) { cerr << setw(12) << M[i][j] << " "; } cerr << endl;
    //}
    //cerr << "Product U W transpose(V) " << endl;;
    for(i=1;i <=nrow;i++) { // U
      for(j=1;j<=ncol;j++) {
	A[i][j]=0.0;
	for(k=1;k <=ncol;k++) 
	  A[i][j]+=U[i][k]*W[k]*V[j][k];
      }
    }
    //for(i=1;i <=nrow;i++) { // U
    //  for(j=1;j <=ncol;j++) {
    //    cerr << setw(12) << A[i][j] << " ";
    //  }
    //  cerr << endl;
    //}
    //cerr << "A and M is equal " << isequal(A,M) << endl;
    if(isequal(A,M) ==false ) {
      cerr << "ERROR - cematrix::SVDcmp_NR: SVD fails, check the code of cematrix::SVDcmp!" << endl;;
      //exit(_EXIT_SVD_FAIL);
      flag_convergence=false;
    }
    return flag_convergence;

  }

  xmatrix<double> cematrix::InverseMatrix() {
    // inverse a general matrix by using SVD
    int i,j,k;
  
    // SVD to get U,V,W
    SVDcmp_NR();//Two functions are essentially the same
    // inverse of matrix
    xmatrix<double> A_inv(1,1,ncol,nrow);
    double DetW=1.0;// determination of diagonal matrix W
    for(i=1;i <=ncol;i++)
      DetW *=W[i];
    if(DetW < cematrix_EQUAL_DOUBLE) {
      cerr << "ERROR - cematrix::InverseMatrix: Singular Matrix. No Inversion" << endl;;
      exit(1);
    }
    for(i=1;i <=ncol;i++) { //row of the inverse Matrix
      for(j=1;j <=nrow;j++) { // colume of the inverse Matrix
	A_inv[i][j]=0.0;
	for(k=1;k <=ncol;k++) // matrix multiplication
	  A_inv[i][j]+=V[i][k]/W[k]*U[j][k];
      }
    }
    cerr << "cematrix::InverseMatrix: Inverse of matrix A" << endl;
    for(i=1;i <=ncol;i++) { // U
      for(j=1;j <=nrow;j++) 
	cerr << setw(12) << A_inv[i][j] << " ";
      cerr << endl;
    }
    xmatrix<double> Iden_tmp(1,1,nrow,nrow);
    xmatrix<double> Iden(1,1,nrow,nrow);
    for(i=1;i<=nrow;i++) {
      for(j=1;j<=nrow;j++) {
	Iden_tmp[i][j]=0.0;
	if(i ==j ) {
	  Iden[i][j]=1.0;
	} else {
	  Iden[i][j]=0.0;
	}
	for(k=1;k<=ncol;k++) {
	  Iden_tmp[i][j]+=M[i][k]*A_inv[k][j];
	}
      }
    }
    cerr << "cematrix::InverseMatrix: A_inv is the inverse of matrix A? ";
    if(isequal(Iden,Iden_tmp) ) {
      cerr << "Yes!" << endl;;
    } else {
      cerr << "No!" << endl;;
    }
    return A_inv;
  }

  void cematrix::SVDFit(xvector<double>& y,xvector<double>& y_sigma) {
    // Least Square fit by using SVD
    // Here x() is afunc() in numerical recipes in C
    int i,j;
    double wmax,tmp,thresh,sum;
    xmatrix<double> A(1,1,nrow,ncol);
    xmatrix<double> M_orig(1,1,nrow,ncol);
    xvector<double> W_orig(1,ncol);
    xmatrix<double> V_orig(1,1,ncol,ncol);
    xmatrix<double> U_orig(1,1,nrow,ncol);
    //const double TOL=1.0e-13;
    const double TOL=1.0e-8;
    xvector<double> y_cal(1,nrow);
    A=M;
    if(A.rows !=y.rows ) {
      cerr << "ERROR - cematrix::SVDFit: Ranks of x vector and y vector does not match!" << endl;;
      exit(cematrix_EXIT_RANK_NOT_MATCH);
    }
    xvector<double> b(1,nrow);
    for(i=1;i <=nrow;i++) { // Accumulate coeeficiens of the fitting matrix
      tmp=1.0/y_sigma[i];
      for(j=1;j<=ncol;j++)
	A[i][j]=M[i][j]*tmp;
      b[i]=y[i]*tmp;
    }
    M_orig=M;
    W_orig=W;
    V_orig=V;
    U_orig=U;
    M=A;
    //SVDcmp_NR();// singular value decomposition
    // singular value decomposition
    if(SVDcmp_NR() ) {
      wmax=0.0;
      for(j=1;j<=ncol;j++) 
	if(W[j]>wmax) wmax=W[j];
      thresh=TOL*wmax;
      for(j=1;j<ncol;j++)
	if(W[j] < thresh ) W[j]=0.0;
      SVDsolve(b);
      chisq=0.0;
      for(i=1;i<=nrow;i++) {
	for(sum=0.0,j=1;j<=ncol;j++)
	  sum+=a_vec[j]*M_orig[i][j];
	tmp=(y[i] - sum)/y_sigma[i];
	chisq+=tmp*tmp;
      }
      chisq=chisq/(nrow - 2.0);// definition in Mathematica
    } else {
      // if svd fails,set the score to a large number
      // to discard it
      chisq=1.0e4;
      for(int i=0;i<y.rows;i++)
	a_nvec.push_back(0.0e0);
    }
    //// output the fitted parameters and square of chi
    //cout << "a_vec " << endl;
    //for(i=1;i<=ncol;i++) {
    //  cout << a_vec[i] << " ";
    //}
    //cout << endl;
    //cout << "chisq " << chisq << endl;
    //cerr << "a_vec " << endl;
    //for(i=1;i<=ncol;i++) {
    //  cerr << a_vec[i] << " ";
    //}
    //cerr << endl;
    //cerr << "chisq " << chisq << endl;
    //cerr << "original y vector" << endl;
    //for(i=1;i<=nrow;i++) {
    //  cerr << y[i] << " ";
    //}
    //  cerr << endl;
    //// fitted value of y vector
    //cerr << "fitted y vector" << endl;
    //for(i=1;i<=nrow;i++) {
    //  for(j=1;j<=ncol;j++) {
    //    y_cal[i]+=A[i][j]*a_vec[j];
    //  }
    //  cerr << y_cal[i] << " ";
    //}
    //cerr << endl;
    //// get the covariance matrix
    //SVDvar();
    //ios_base::fmtflags old_stat=cerr.setf(ios_base::fixed,ios_base::floatfield);
    //cerr.setf(ios_base::scientific);
    //cerr.precision(6);
    //cerr << "Covariance matrix " << endl;;
    //for(i=1;i<=ncol;i++) {
    //  for(j=1;j<=ncol;j++) {
    //    cerr << setw(12) << Cov[i][j] << " ";
    //  }
    //  cerr << endl;
    //}

    // set everything back to original values
    M=M_orig;
    U=U_orig;
    V=V_orig;
    W=W_orig;
    //cerr.setf(old_stat,ios_base::floatfield);
    //cerr.unsetf(ios_base::scientific);

  }

  void cematrix::SVDvar() {
    // get the covariance matrix Cov
    xmatrix<double> Cov_tmp(1,1,ncol,ncol);
    int k,j,i;
    double sum;
    xvector<double> wti(1,ncol);
    for(i=1;i<= ncol;i++) {
      wti[i]=0.0;
      if(W[i] !=0.0 )
	wti[i]=1.0/(W[i]*W[i]);
    }
    for(i=1;i<=ncol;i++) { // sum contributions to covariance matrix
      for(j=1;j<=i;j++) {
	for(sum=0.0,k=1;k<=ncol;k++)
	  sum+=V[i][k]*V[j][k]*wti[k];
	Cov_tmp[j][i]=sum;
	Cov_tmp[i][j]=sum;
      }
    }
    Cov=Cov_tmp*chisq;// definition in Mathematica
  }

  xvector<double> cematrix::EigenValues() {
    // inverse a general matrix by using SVD
    // SVD to get U,V,W
    SVDcmp_NR(); //Two functions are essentially the same
    return W;
  }
}

// **************************************************************************

#endif
// **************************************************************************
// *                                                                        *
// *             STEFANO CURTAROLO - Duke University 2003-2015              *
// *                                                                        *
// **************************************************************************

