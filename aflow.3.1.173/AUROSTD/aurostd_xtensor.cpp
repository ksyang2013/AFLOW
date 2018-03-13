// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2015           *
// *                                                                         *
// ***************************************************************************
// Written by Stefano Curtarolo Nov07
// fixed for g++ 4.5 on Mar11

// ---------------------------------------------------------------------------
// ------------------ implementation for template<class utype> xtensor6<utype>

#ifndef _AUROSTD_XTENSOR_CPP_
#define _AUROSTD_XTENSOR_CPP_

#ifndef XXEND
#define XXEND 1
#endif

#ifndef _AUROSTD_XSCALAR_H_
#include "aurostd_xscalar.h"
#endif
#ifndef _AUROSTD_XCOMPLEX_H_
#include "aurostd_xcomplex.h"
#endif
#ifndef _AUROSTD_XVECTOR_H_
#include "aurostd_xvector.h"
#endif
#ifndef _AUROSTD_XMATRIX_H_
#include "aurostd_xmatrix.h"
#endif
#ifndef _AUROSTD_XTENSOR_H_
#include "aurostd_xtensor.h"
#endif

using std::cout;
using std::cerr;

#define _AUROSTD_XTENSOR3_INDEX_   3
#define _AUROSTD_XTENSOR4_INDEX_   4
#define _AUROSTD_XTENSOR5_INDEX_   5
#define _AUROSTD_XTENSOR6_INDEX_   6
#define _AUROSTD_XTENSOR7_INDEX_   7
#define _AUROSTD_XTENSOR8_INDEX_   8

// ***************************************************************************
// ---------------------------------------------------------------------------

namespace aurostd {  // namespace aurostd
  int eijk(int i,int j,int k) {
    int ii=(i-1)%3+1;
    int jj=(j-1)%3+1;
    int kk=(k-1)%3+1;
    if(ii==1 and jj==2 and kk==3) return  1;
    if(ii==2 and jj==3 and kk==1) return  1;
    if(ii==3 and jj==1 and kk==2) return  1;
    if(ii==1 and jj==3 and kk==2) return -1;
    if(ii==3 and jj==2 and kk==1) return -1;
    if(ii==2 and jj==1 and kk==3) return -1;
    return 0;
  }
  // namespace aurostd
  int eijk(xvector<int> ijk) {
    return eijk(ijk[1],ijk[2],ijk[3]);
  }
  // namespace aurostd
  int estarijk(int i,int j,int k) {
    int ii=(i-1)%3+1;
    int jj=(j-1)%3+1;
    int kk=(k-1)%3+1;
    if(ii==1 and jj==2 and kk==3) return  1;
    if(ii==2 and jj==3 and kk==1) return  1;
    if(ii==3 and jj==1 and kk==2) return  1;
    return 0;
  }
  // namespace aurostd
  int estarijk(xvector<int> ijk) {
    return estarijk(ijk[1],ijk[2],ijk[3]);
  }
}

// ***************************************************************************
// --------------------------------------------------------------------- debug

namespace aurostd {  // namespace aurostd
  template<class utype>
  void xtensor3debug(const xtensor3<utype>& t,const string& str) {
#ifdef _XMATH_DEBUG_CONSTRUCTORS
    cout << "XTENSOR3 -> " << string << ":";
    for(int i=1;i<=_AUROSTD_XTENSOR3_INDEX_;i++)
      cout << "  lindex[i]=" << t.lindex[i] << ", uindex[i]=" << t.uindex[i];
    cout << endl;
#else
    if(t.index[1]) {;} // phony to keep t busy
    if(str.length()) {;} // phony to keep str busy
#endif
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype>
  void xtensor4debug(const xtensor4<utype>& t,const string& str) {
#ifdef _XMATH_DEBUG_CONSTRUCTORS
    cout << "XTENSOR4 -> " << string << ":";
    for(int i=1;i<=_AUROSTD_XTENSOR4_INDEX_;i++)
      cout << "  lindex[i]=" << t.lindex[i] << ", uindex[i]=" << t.uindex[i];
    cout << endl;
#else
    if(t.index[1]) {;} // phony to keep t busy
    if(str.length()) {;} // phony to keep str busy
#endif
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype>
  void xtensor5debug(const xtensor5<utype>& t,const string& str) {
#ifdef _XMATH_DEBUG_CONSTRUCTORS
    cout << "XTENSOR5 -> " << string << ":";
    for(int i=1;i<=_AUROSTD_XTENSOR5_INDEX_;i++)
      cout << "  lindex[i]=" << t.lindex[i] << ", uindex[i]=" << t.uindex[i];
    cout << endl;
#else
    if(t.index[1]) {;} // phony to keep t busy
    if(str.length()) {;} // phony to keep str busy
#endif
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype>
  void xtensor6debug(const xtensor6<utype>& t,const string& str) {
#ifdef _XMATH_DEBUG_CONSTRUCTORS
    cout << "XTENSOR6 -> " << string << ":";
    for(int i=1;i<=_AUROSTD_XTENSOR6_INDEX_;i++)
      cout << "  lindex[i]=" << t.lindex[i] << ", uindex[i]=" << t.uindex[i];
    cout << endl;
#else
    if(t.index[1]) {;} // phony to keep t busy
    if(str.length()) {;} // phony to keep str busy
#endif
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype>
  void xtensor7debug(const xtensor7<utype>& t,const string& str) {
#ifdef _XMATH_DEBUG_CONSTRUCTORS
    cout << "XTENSOR7 -> " << string << ":";
    for(int i=1;i<=_AUROSTD_XTENSOR7_INDEX_;i++)
      cout << "  lindex[i]=" << t.lindex[i] << ", uindex[i]=" << t.uindex[i];
    cout << endl;
#else
    if(t.index[1]) {;} // phony to keep t busy
    if(str.length()) {;} // phony to keep str busy
#endif
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype>
  void xtensor8debug(const xtensor8<utype>& t,const string& str) {
#ifdef _XMATH_DEBUG_CONSTRUCTORS
    cout << "XTENSOR8 -> " << string << ":";
    for(int i=1;i<=_AUROSTD_XTENSOR8_INDEX_;i++)
      cout << "  lindex[i]=" << t.lindex[i] << ", uindex[i]=" << t.uindex[i];
    cout << endl;
#else
    if(t.index[1]) {;} // phony to keep t busy
    if(str.length()) {;} // phony to keep str busy
#endif
  }
}

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

// ***************************************************************************
// --------------------------------------------- allocate_xtensor3corpus xtensor3
namespace aurostd {  // namespace aurostd
  template<class utype>
  void allocate_xtensor3corpus(utype*** &corpus, int* lindex,int* uindex,int* index) {
    int i,j;
    // MEMORY -------------------
    corpus=new utype **[index[1]+XXEND];
    if(!corpus) {cerr << "XTENSOR: allocation failure 1 in xtensor4 constructor" << endl;exit(0);}
    corpus+=-lindex[1]+XXEND;
    for(i=lindex[1];i<=uindex[1];i++) {
      corpus[i]=new utype *[index[2]+XXEND];
      if(!corpus[i]) {cerr << "XTENSOR: aallocation failure 2. " << i << " in xtensor3 constructor" << endl;exit(0);}
      corpus[i]+=-lindex[2]+XXEND;
      corpus[i][lindex[2]]= new utype[index[2]*index[3]+XXEND];
      if(!corpus[i][lindex[2]]) {cerr << "XTENSOR: allocation failure 3." << i << " in xtensor3 constructor" << endl;exit(0);}
      corpus[i][lindex[2]]+=-lindex[3]+XXEND;
      for(j=lindex[2]+1;j<=uindex[2];j++)
	corpus[i][j]=corpus[i][j-1]+index[3];
    }
    // DONE ------------------
  }
}

// --------------------------------------------- allocate_xtensor4corpus xtensor4
namespace aurostd {  // namespace aurostd
  template<class utype>
  void allocate_xtensor4corpus(utype**** &corpus, int* lindex,int* uindex,int* index) {
    int i,j,k;
    // MEMORY -------------------
    corpus=new utype ***[index[1]+XXEND];
    if(!corpus) {cerr << "XTENSOR: allocation failure 1 in xtensor4 constructor" << endl;exit(0);}
    corpus+=-lindex[1]+XXEND;
    for(i=lindex[1];i<=uindex[1];i++) {
      corpus[i]=new utype **[index[2]+XXEND];
      if(!corpus[i]) {cerr << "XTENSOR: aallocation failure 2. " << i << " in xtensor4 constructor" << endl;exit(0);}
      corpus[i]+=-lindex[2]+XXEND;
      for(j=lindex[2];j<=uindex[2];j++) {
	corpus[i][j]=new utype *[index[3]+XXEND];
	if(!corpus[i][j]) {cerr << "XTENSOR: allocation failure 3." << i << "." << j << " in xtensor4 constructor" << endl;exit(0);}
	corpus[i][j]+=-lindex[3]+XXEND;
	corpus[i][j][lindex[3]]= new utype[index[3]*index[4]+XXEND];
	if(!corpus[i][j][lindex[3]]) {cerr << "XTENSOR: allocation failure 4." << i << "." << j  << " in xtensor4 constructor" << endl;exit(0);}
	corpus[i][j][lindex[3]]+=-lindex[4]+XXEND;
	for(k=lindex[3]+1;k<=uindex[3];k++)
	  corpus[i][j][k]=corpus[i][j][k-1]+index[4];
      }
    }
    // DONE ------------------
  }
}

// --------------------------------------------- allocate_xtensor5corpus xtensor5
namespace aurostd {  // namespace aurostd
  template<class utype>
  void allocate_xtensor5corpus(utype***** &corpus, int* lindex,int* uindex,int* index) {
    int i,j,k,l;
    // MEMORY -------------------
    corpus=new utype ****[index[1]+XXEND];
    if(!corpus) {cerr << "XTENSOR: allocation failure 1 in xtensor5 constructor" << endl;exit(0);}
    corpus+=-lindex[1]+XXEND;
    for(i=lindex[1];i<=uindex[1];i++) {
      corpus[i]=new utype ***[index[2]+XXEND];
      if(!corpus[i]) {cerr << "XTENSOR: aallocation failure 2. " << i << " in xtensor5 constructor" << endl;exit(0);}
      corpus[i]+=-lindex[2]+XXEND;
      for(j=lindex[2];j<=uindex[2];j++) {
	corpus[i][j]=new utype **[index[3]+XXEND];
	if(!corpus[i][j]) {cerr << "XTENSOR: allocation failure 3." << i << "." << j << " in xtensor5 constructor" << endl;exit(0);}
	corpus[i][j]+=-lindex[3]+XXEND;
	for(k=lindex[3];k<=uindex[3];k++) {
	  corpus[i][j][k] = new utype *[index[4]+XXEND];
	  if(!corpus[i][j][k]) { cerr << "XTENSOR: allocation failure 4." << i << "." << j << "." << k << " in xtensor5 constructor" << endl;exit(0);}
	  corpus[i][j][k]+=-lindex[4]+XXEND;
	  corpus[i][j][k][lindex[4]]= new utype[index[4]*index[5]+XXEND];
	  if(!corpus[i][j][k][lindex[4]]) {cerr << "XTENSOR: allocation failure 5." << i << "." << j << "."  << k << " in xtensor4 constructor" << endl;exit(0);}
	  corpus[i][j][k][lindex[4]]+=-lindex[5]+XXEND;
	  for(l=lindex[4]+1;l<=uindex[4];l++)
	    corpus[i][j][k][l]=corpus[i][j][k][l-1]+index[5];
	}
      }
    }
    // DONE ------------------
  }
}

// --------------------------------------------- allocate_xtensor6corpus xtensor6
namespace aurostd {  // namespace aurostd
  template<class utype>
  void allocate_xtensor6corpus(utype****** &corpus, int* lindex,int* uindex,int* index) {
    int i,j,k,l,m;
    // MEMORY -------------------
    corpus=new utype *****[index[1]+XXEND];
    if(!corpus) {cerr << "XTENSOR: allocation failure 1 in xtensor6 constructor" << endl;exit(0);}
    corpus+=-lindex[1]+XXEND;
    for(i=lindex[1];i<=uindex[1];i++) {
      corpus[i]=new utype ****[index[2]+XXEND];
      if(!corpus[i]) {cerr << "XTENSOR: aallocation failure 2. " << i << " in xtensor6 constructor" << endl;exit(0);}
      corpus[i]+=-lindex[2]+XXEND;
      for(j=lindex[2];j<=uindex[2];j++) {
	corpus[i][j]=new utype ***[index[3]+XXEND];
	if(!corpus[i][j]) {cerr << "XTENSOR: allocation failure 3." << i << "." << j << " in xtensor6 constructor" << endl;exit(0);}
	corpus[i][j]+=-lindex[3]+XXEND;
	for(k=lindex[3];k<=uindex[3];k++) {
	  corpus[i][j][k] = new utype **[index[4]+XXEND];
	  if(!corpus[i][j][k]) {cerr << "XTENSOR: allocation failure 4." << i << "." << j << "." << k << " in xtensor6 constructor" << endl;exit(0);}
	  corpus[i][j][k]+=-lindex[4]+XXEND;
	  for(l=lindex[4];l<=uindex[4];l++) {
	    corpus[i][j][k][l] = new utype *[index[5]+XXEND];
	    if(!corpus[i][j][k][l]) {cerr << "XTENSOR: allocation failure 5." << i << "." << j << "." << k << "." << l << " in xtensor6 constructor" << endl;exit(0);}
	    corpus[i][j][k][l]+=-lindex[5]+XXEND;	
	    corpus[i][j][k][l][lindex[5]]= new utype[index[5]*index[6]+XXEND];
	    if(!corpus[i][j][k][l][lindex[5]]) {cerr << "XTENSOR: allocation failure 6." << i << "." << j << "." << k << "." << l << " in xtensor6 constructor" << endl;exit(0);}
	    corpus[i][j][k][l][lindex[5]]+=-lindex[6]+XXEND;
	    for(m=lindex[5]+1;m<=uindex[5];m++)
	      corpus[i][j][k][l][m]=corpus[i][j][k][l][m-1]+index[6];
	  }
	}
      }
    }
  }
}

// --------------------------------------------- allocate_xtensor7corpus xtensor7
namespace aurostd {  // namespace aurostd
  template<class utype>
  void allocate_xtensor7corpus(utype******* &corpus, int* lindex,int* uindex,int* index) {
    int i,j,k,l,m,n;
    // MEMORY -------------------
    corpus=new utype ******[index[1]+XXEND];
    if(!corpus) {cerr << "XTENSOR: allocation failure 1 in xtensor7 constructor" << endl;exit(0);}
    corpus+=-lindex[1]+XXEND;
    for(i=lindex[1];i<=uindex[1];i++) {
      corpus[i]=new utype *****[index[2]+XXEND];
      if(!corpus[i]) {cerr << "XTENSOR: aallocation failure 2. " << i << " in xtensor7 constructor" << endl;exit(0);}
      corpus[i]+=-lindex[2]+XXEND;
      for(j=lindex[2];j<=uindex[2];j++) {
	corpus[i][j]=new utype ****[index[3]+XXEND];
	if(!corpus[i][j]) {cerr << "XTENSOR: allocation failure 3." << i << "." << j << " in xtensor7 constructor" << endl;exit(0);}
	corpus[i][j]+=-lindex[3]+XXEND;
	for(k=lindex[3];k<=uindex[3];k++) {
	  corpus[i][j][k] = new utype ***[index[4]+XXEND];
	  if(!corpus[i][j][k]) {cerr << "XTENSOR: allocation failure 4." << i << "." << j << "." << k << " in xtensor7 constructor" << endl;exit(0);}
	  corpus[i][j][k]+=-lindex[4]+XXEND;
	  for(l=lindex[4];l<=uindex[4];l++) {
	    corpus[i][j][k][l] = new utype **[index[5]+XXEND];
	    if(!corpus[i][j][k][l]) {cerr << "XTENSOR: allocation failure 5." << i << "." << j << "." << k << "." << l << " in xtensor7 constructor" << endl;exit(0);}
	    corpus[i][j][k][l]+=-lindex[5]+XXEND;
	    for(m=lindex[5];m<=uindex[5];m++) {
	      corpus[i][j][k][l][m] = new utype *[index[6]+XXEND];
	      if(!corpus[i][j][k][l][m]) {cerr << "XTENSOR: allocation failure 6." << i << "." << j << "." << k << "." << l << "." << m << " in xtensor7 constructor" << endl;exit(0);}
	      corpus[i][j][k][l][m]+=-lindex[6]+XXEND;	
	      corpus[i][j][k][l][m][lindex[6]]= new utype[index[6]*index[7]+XXEND];
	      if(!corpus[i][j][k][l][m][lindex[6]]) {cerr << "XTENSOR: allocation failure 7." << i << "." << j << "." << k << "." << l << "." << m << " in xtensor7 constructor" << endl;exit(0);}
	      corpus[i][j][k][l][m][lindex[6]]+=-lindex[7]+XXEND;
	      for(n=lindex[6]+1;n<=uindex[6];n++)
		corpus[i][j][k][l][m][n]=corpus[i][j][k][l][m][n-1]+index[7];
	    }
	  }
	}
      }
    }
  }
}

// --------------------------------------------- allocate_xtensor8corpus xtensor8
namespace aurostd {  // namespace aurostd
  template<class utype>
  void allocate_xtensor8corpus(utype******** &corpus, int* lindex,int* uindex,int* index) {
    int i,j,k,l,m,n,o;
    // MEMORY -------------------
    corpus=new utype *******[index[1]+XXEND];
    if(!corpus) {cerr << "XTENSOR: allocation failure 1 in xtensor8 constructor" << endl;exit(0);}
    corpus+=-lindex[1]+XXEND;
    for(i=lindex[1];i<=uindex[1];i++) {
      corpus[i]=new utype ******[index[2]+XXEND];
      if(!corpus[i]) {cerr << "XTENSOR: aallocation failure 2. " << i << " in xtensor8 constructor" << endl;exit(0);}
      corpus[i]+=-lindex[2]+XXEND;
      for(j=lindex[2];j<=uindex[2];j++) {
	corpus[i][j]=new utype *****[index[3]+XXEND];
	if(!corpus[i][j]) {cerr << "XTENSOR: allocation failure 3." << i << "." << j << " in xtensor8 constructor" << endl;exit(0);}
	corpus[i][j]+=-lindex[3]+XXEND;
	for(k=lindex[3];k<=uindex[3];k++) {
	  corpus[i][j][k] = new utype ****[index[4]+XXEND];
	  if(!corpus[i][j][k]) {cerr << "XTENSOR: allocation failure 4." << i << "." << j << "." << k << " in xtensor8 constructor" << endl;exit(0);}
	  corpus[i][j][k]+=-lindex[4]+XXEND;
	  for(l=lindex[4];l<=uindex[4];l++) {
	    corpus[i][j][k][l] = new utype ***[index[5]+XXEND];
	    if(!corpus[i][j][k][l]) {cerr << "XTENSOR: allocation failure 5." << i << "." << j << "." << k << "." << l << " in xtensor8 constructor" << endl;exit(0);}
	    corpus[i][j][k][l]+=-lindex[5]+XXEND;
	    for(m=lindex[5];m<=uindex[5];m++) {
	      corpus[i][j][k][l][m] = new utype **[index[6]+XXEND];
	      if(!corpus[i][j][k][l][m]) {cerr << "XTENSOR: allocation failure 6." << i << "." << j << "." << k << "." << l << "." << m << " in xtensor8 constructor" << endl;exit(0);}
	      corpus[i][j][k][l][m]+=-lindex[6]+XXEND;	
	      for(n=lindex[6];n<=uindex[6];n++) {
		corpus[i][j][k][l][m][n] = new utype *[index[7]+XXEND];
		if(!corpus[i][j][k][l][m][n]) {cerr << "XTENSOR: allocation failure 7." << i << "." << j << "." << k << "." << l << "." << m << "." << n << " in xtensor8 constructor" << endl;exit(0);}
		corpus[i][j][k][l][m][n]+=-lindex[7]+XXEND;	
		corpus[i][j][k][l][m][n][lindex[7]]= new utype[index[7]*index[8]+XXEND];
		if(!corpus[i][j][k][l][m][n][lindex[7]]) {cerr << "XTENSOR: allocation failure 8." << i << "." << j << "." << k << "." << l << "." << m << "." << n << " in xtensor8 constructor" << endl;exit(0);}
		corpus[i][j][k][l][m][n][lindex[7]]+=-lindex[8]+XXEND;
		for(o=lindex[7]+1;o<=uindex[7];o++)
		  corpus[i][j][k][l][m][n][o]=corpus[i][j][k][l][m][n][o-1]+index[8];
	      }
	    }
	  }
	}
      }
    }
  }
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

// ***************************************************************************
// ----------------------------------------------------- constructors xtensor3
namespace aurostd {  // namespace aurostd
  template<class utype>                                 // default constructor
  xtensor3<utype>::xtensor3(int index1h,int index2h,int index3h,
			    int index1l,int index2l,int index3l) {
    int i,j,k;
    lindex[1]=min(index1l,index1h);uindex[1]=max(index1l,index1h);index[1]=uindex[1]-lindex[1]+1;
    lindex[2]=min(index2l,index2h);uindex[2]=max(index2l,index2h);index[2]=uindex[2]-lindex[2]+1;
    lindex[3]=min(index3l,index3h);uindex[3]=max(index3l,index3h);index[3]=uindex[3]-lindex[3]+1;
    iscubic=TRUE;
    for(i=1;i<_AUROSTD_XTENSOR3_INDEX_;i++) {iscubic=iscubic && bool(index[i]==index[i+1]);};
    aurostd::xtensor3debug(*this,"default constructor");
    isfloat=_isfloat((utype) 0);
    iscomplex=_iscomplex((utype) 0);
    size=(char) (sizeof(utype));
    tsize=(long int) size;for(i=1;i<=_AUROSTD_XTENSOR3_INDEX_;i++) tsize*=index[i];
    //cout << "XTENSOR3 size=" << tsize << endl;
    // MEMORY -------------------
    aurostd::allocate_xtensor3corpus(corpus,lindex,uindex,index);
    // CLEAR ------------------
    for(i=lindex[1];i<=uindex[1];i++)
      for(j=lindex[2];j<=uindex[2];j++)
	for(k=lindex[3];k<=uindex[3];k++)
	  corpus[i][j][k]=(utype) 0; //i,j,k
    // DONE ------------------
#ifdef _XMATH_DEBUG_CONSTRUCTORS
    cout << "iscubic=" << iscubic << ", isfloat=" << isfloat << ", iscomplex=" << iscomplex << ", sizeof=" << size << ", tsize=" << tsize << endl;
#endif
  }
}

// ----------------------------------------------------- constructors xtensor4
namespace aurostd {  // namespace aurostd
  template<class utype>                                 // default constructor
  xtensor4<utype>::xtensor4(int index1h,int index2h,int index3h,int index4h,
			    int index1l,int index2l,int index3l,int index4l) {
    int i,j,k,l;
    lindex[1]=min(index1l,index1h);uindex[1]=max(index1l,index1h);index[1]=uindex[1]-lindex[1]+1;
    lindex[2]=min(index2l,index2h);uindex[2]=max(index2l,index2h);index[2]=uindex[2]-lindex[2]+1;
    lindex[3]=min(index3l,index3h);uindex[3]=max(index3l,index3h);index[3]=uindex[3]-lindex[3]+1;
    lindex[4]=min(index4l,index4h);uindex[4]=max(index4l,index4h);index[4]=uindex[4]-lindex[4]+1;
    iscubic=TRUE;
    for(i=1;i<_AUROSTD_XTENSOR4_INDEX_;i++) {iscubic=iscubic && bool(index[i]==index[i+1]);};
    aurostd::xtensor4debug(*this,"default constructor");
    isfloat=_isfloat((utype) 0);
    iscomplex=_iscomplex((utype) 0);
    size=(char) (sizeof(utype));
    tsize=(long int) size;for(i=1;i<=_AUROSTD_XTENSOR4_INDEX_;i++) tsize*=index[i];
    //cout << "XTENSOR4 size=" << tsize << endl;
    // MEMORY -------------------
    aurostd::allocate_xtensor4corpus(corpus,lindex,uindex,index);
    // CLEAR ------------------
    for(i=lindex[1];i<=uindex[1];i++)
      for(j=lindex[2];j<=uindex[2];j++)
	for(k=lindex[3];k<=uindex[3];k++)
	  for(l=lindex[4];l<=uindex[4];l++)
	    corpus[i][j][k][l]=(utype) 0; //i,j,k,l
    // DONE ------------------
#ifdef _XMATH_DEBUG_CONSTRUCTORS
    cout << "iscubic=" << iscubic << ", isfloat=" << isfloat << ", iscomplex=" << iscomplex << ", sizeof=" << size << ", tsize=" << tsize << endl;
#endif
  }
}

// ----------------------------------------------------- constructors xtensor5
namespace aurostd {  // namespace aurostd
  template<class utype>                                  // default constructor
  xtensor5<utype>::xtensor5(int index1h,int index2h,int index3h,int index4h,int index5h,
			    int index1l,int index2l,int index3l,int index4l,int index5l) {
    int i,j,k,l,m;
    lindex[1]=min(index1l,index1h);uindex[1]=max(index1l,index1h);index[1]=uindex[1]-lindex[1]+1;
    lindex[2]=min(index2l,index2h);uindex[2]=max(index2l,index2h);index[2]=uindex[2]-lindex[2]+1;
    lindex[3]=min(index3l,index3h);uindex[3]=max(index3l,index3h);index[3]=uindex[3]-lindex[3]+1;
    lindex[4]=min(index4l,index4h);uindex[4]=max(index4l,index4h);index[4]=uindex[4]-lindex[4]+1;
    lindex[5]=min(index5l,index5h);uindex[5]=max(index5l,index5h);index[5]=uindex[5]-lindex[5]+1;
    iscubic=TRUE;
    for(i=1;i<_AUROSTD_XTENSOR5_INDEX_;i++) {iscubic=iscubic && bool(index[i]==index[i+1]);};
    aurostd::xtensor5debug(*this,"default constructor");
    isfloat=_isfloat((utype) 0);
    iscomplex=_iscomplex((utype) 0);
    size=(char) (sizeof(utype));
    tsize=(long int) size;for(i=1;i<=_AUROSTD_XTENSOR5_INDEX_;i++) tsize*=index[i];
    //cout << "XTENSOR5 size=" << tsize << endl;
    // MEMORY -------------------
    aurostd::allocate_xtensor5corpus(corpus,lindex,uindex,index);
    // CLEAR ------------------
    for(i=lindex[1];i<=uindex[1];i++)
      for(j=lindex[2];j<=uindex[2];j++)
	for(k=lindex[3];k<=uindex[3];k++)
	  for(l=lindex[4];l<=uindex[4];l++)
	    for(m=lindex[5];m<=uindex[5];m++)
	      corpus[i][j][k][l][m]=(utype) 0; //i,j,k,l,m
    // DONE ------------------
#ifdef _XMATH_DEBUG_CONSTRUCTORS
    cout << "iscubic=" << iscubic << ", isfloat=" << isfloat << ", iscomplex=" << iscomplex << ", sizeof=" << size << ", tsize=" << tsize << endl;
#endif
  }
}

// ----------------------------------------------------- constructors xtensor6
namespace aurostd {  // namespace aurostd
  template<class utype>                                  // default constructor
  xtensor6<utype>::xtensor6(int index1h,int index2h,int index3h,int index4h,
			    int index5h,int index6h,int index1l,int index2l,
			    int index3l,int index4l,int index5l,int index6l) {
    int i,j,k,l,m,n;
    lindex[1]=min(index1l,index1h);uindex[1]=max(index1l,index1h);index[1]=uindex[1]-lindex[1]+1;
    lindex[2]=min(index2l,index2h);uindex[2]=max(index2l,index2h);index[2]=uindex[2]-lindex[2]+1;
    lindex[3]=min(index3l,index3h);uindex[3]=max(index3l,index3h);index[3]=uindex[3]-lindex[3]+1;
    lindex[4]=min(index4l,index4h);uindex[4]=max(index4l,index4h);index[4]=uindex[4]-lindex[4]+1;
    lindex[5]=min(index5l,index5h);uindex[5]=max(index5l,index5h);index[5]=uindex[5]-lindex[5]+1;
    lindex[6]=min(index6l,index6h);uindex[6]=max(index6l,index6h);index[6]=uindex[6]-lindex[6]+1;
    iscubic=TRUE;
    for(i=1;i<_AUROSTD_XTENSOR6_INDEX_;i++) {iscubic=iscubic && bool(index[i]==index[i+1]);};
    aurostd::xtensor6debug(*this,"default constructor");
    isfloat=_isfloat((utype) 0);
    iscomplex=_iscomplex((utype) 0);
    size=(char) (sizeof(utype));
    tsize=(long int) size;for(i=1;i<=_AUROSTD_XTENSOR6_INDEX_;i++) tsize*=index[i];
    //cout << "XTENSOR6 size=" << tsize << endl;
    // MEMORY -------------------
    aurostd::allocate_xtensor6corpus(corpus,lindex,uindex,index);
    // CLEAR ------------------
    for(i=lindex[1];i<=uindex[1];i++)
      for(j=lindex[2];j<=uindex[2];j++)
	for(k=lindex[3];k<=uindex[3];k++)
	  for(l=lindex[4];l<=uindex[4];l++)
	    for(m=lindex[5];m<=uindex[5];m++)
	      for(n=lindex[6];n<=uindex[6];n++)
		corpus[i][j][k][l][m][n]=(utype) 0; //i,j,k,l,m,n
    // DONE ------------------
#ifdef _XMATH_DEBUG_CONSTRUCTORS
    cout << "iscubic=" << iscubic << ", isfloat=" << isfloat << ", iscomplex=" << iscomplex << ", sizeof=" << size << ", tsize=" << tsize << endl;
#endif
  }
}

// ----------------------------------------------------- constructors xtensor7
namespace aurostd {  // namespace aurostd
  template<class utype>                                  // default constructor
  xtensor7<utype>::xtensor7(int index1h,int index2h,int index3h,int index4h,
			    int index5h,int index6h,int index7h,int index1l,int index2l,
			    int index3l,int index4l,int index5l,int index6l,int index7l) {
    int i,j,k,l,m,n,o;
    lindex[1]=min(index1l,index1h);uindex[1]=max(index1l,index1h);index[1]=uindex[1]-lindex[1]+1;
    lindex[2]=min(index2l,index2h);uindex[2]=max(index2l,index2h);index[2]=uindex[2]-lindex[2]+1;
    lindex[3]=min(index3l,index3h);uindex[3]=max(index3l,index3h);index[3]=uindex[3]-lindex[3]+1;
    lindex[4]=min(index4l,index4h);uindex[4]=max(index4l,index4h);index[4]=uindex[4]-lindex[4]+1;
    lindex[5]=min(index5l,index5h);uindex[5]=max(index5l,index5h);index[5]=uindex[5]-lindex[5]+1;
    lindex[6]=min(index6l,index6h);uindex[6]=max(index6l,index6h);index[6]=uindex[6]-lindex[6]+1;
    lindex[7]=min(index7l,index7h);uindex[7]=max(index7l,index7h);index[7]=uindex[7]-lindex[7]+1;
    iscubic=TRUE;
    for(i=1;i<_AUROSTD_XTENSOR7_INDEX_;i++) {iscubic=iscubic && bool(index[i]==index[i+1]);};
    aurostd::xtensor7debug(*this,"default constructor");
    isfloat=_isfloat((utype) 0);
    iscomplex=_iscomplex((utype) 0);
    size=(char) (sizeof(utype));
    tsize=(long int) size;for(i=1;i<=_AUROSTD_XTENSOR7_INDEX_;i++) tsize*=index[i];
    //cout << "XTENSOR7 size=" << tsize << endl;
    // MEMORY -------------------
    aurostd::allocate_xtensor7corpus(corpus,lindex,uindex,index);
    // CLEAR ------------------
    for(i=lindex[1];i<=uindex[1];i++)
      for(j=lindex[2];j<=uindex[2];j++)
	for(k=lindex[3];k<=uindex[3];k++)
	  for(l=lindex[4];l<=uindex[4];l++)
	    for(m=lindex[5];m<=uindex[5];m++)
	      for(n=lindex[6];n<=uindex[6];n++)
		for(o=lindex[7];o<=uindex[7];o++)
		  corpus[i][j][k][l][m][n][o]=(utype) 0; //i,j,k,l,m,n
    // DONE ------------------
#ifdef _XMATH_DEBUG_CONSTRUCTORS
    cout << "iscubic=" << iscubic << ", isfloat=" << isfloat << ", iscomplex=" << iscomplex << ", sizeof=" << size << ", tsize=" << tsize << endl;
#endif
  }
}

// ----------------------------------------------------- constructors xtensor8
namespace aurostd {  // namespace aurostd
  template<class utype>                                  // default constructor
  xtensor8<utype>::xtensor8(int index1h,int index2h,int index3h,int index4h,
			    int index5h,int index6h,int index7h,int index8h,
			    int index1l,int index2l,int index3l,int index4l,
			    int index5l,int index6l,int index7l,int index8l) {
    int i,j,k,l,m,n,o,p;
    lindex[1]=min(index1l,index1h);uindex[1]=max(index1l,index1h);index[1]=uindex[1]-lindex[1]+1;
    lindex[2]=min(index2l,index2h);uindex[2]=max(index2l,index2h);index[2]=uindex[2]-lindex[2]+1;
    lindex[3]=min(index3l,index3h);uindex[3]=max(index3l,index3h);index[3]=uindex[3]-lindex[3]+1;
    lindex[4]=min(index4l,index4h);uindex[4]=max(index4l,index4h);index[4]=uindex[4]-lindex[4]+1;
    lindex[5]=min(index5l,index5h);uindex[5]=max(index5l,index5h);index[5]=uindex[5]-lindex[5]+1;
    lindex[6]=min(index6l,index6h);uindex[6]=max(index6l,index6h);index[6]=uindex[6]-lindex[6]+1;
    lindex[7]=min(index7l,index7h);uindex[7]=max(index7l,index7h);index[7]=uindex[7]-lindex[7]+1;
    lindex[8]=min(index8l,index8h);uindex[8]=max(index8l,index8h);index[8]=uindex[8]-lindex[8]+1;
    iscubic=TRUE;
    for(i=1;i<_AUROSTD_XTENSOR8_INDEX_;i++) {iscubic=iscubic && bool(index[i]==index[i+1]);};
    aurostd::xtensor8debug(*this,"default constructor");
    isfloat=_isfloat((utype) 0);
    iscomplex=_iscomplex((utype) 0);
    size=(char) (sizeof(utype));
    tsize=(long int) size;for(i=1;i<=_AUROSTD_XTENSOR8_INDEX_;i++) tsize*=index[i];
    //cout << "XTENSOR8 size=" << tsize << endl;
    // MEMORY -------------------
    aurostd::allocate_xtensor8corpus(corpus,lindex,uindex,index);
    // CLEAR ------------------
    for(i=lindex[1];i<=uindex[1];i++)
      for(j=lindex[2];j<=uindex[2];j++)
	for(k=lindex[3];k<=uindex[3];k++)
	  for(l=lindex[4];l<=uindex[4];l++)
	    for(m=lindex[5];m<=uindex[5];m++)
	      for(n=lindex[6];n<=uindex[6];n++)
		for(o=lindex[7];o<=uindex[7];o++)
		  for(p=lindex[8];p<=uindex[8];p++)
		    corpus[i][j][k][l][m][n][o][p]=(utype) 0; //i,j,k,l,m,n,o,p
    // DONE ------------------
#ifdef _XMATH_DEBUG_CONSTRUCTORS
    cout << "iscubic=" << iscubic << ", isfloat=" << isfloat << ", iscomplex=" << iscomplex << ", sizeof=" << size << ", tsize=" << tsize << endl;
#endif
  }
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

// ***************************************************************************
// ------------------------------------------------ copy constructors xtensor3
namespace aurostd {  // namespace aurostd
  template<class utype>                                    // copy constructor
  xtensor3<utype>::xtensor3(const xtensor3<utype>& a) {
    for(int i=1;i<=_AUROSTD_XTENSOR3_INDEX_;i++) {lindex[i]=a.lindex[i];uindex[i]=a.uindex[i];index[i]=uindex[i]-lindex[i]+1;}
    iscubic=TRUE;
    for(int i=1;i<_AUROSTD_XTENSOR3_INDEX_;i++) {iscubic=iscubic && bool(index[i]==index[i+1]);};
    aurostd::xtensor3debug(*this,"copy constructor");
    isfloat=_isfloat((utype) 0);
    iscomplex=_iscomplex((utype) 0);
    size=(char) (sizeof(utype));
    tsize=(long int) size;for(int i=1;i<=_AUROSTD_XTENSOR3_INDEX_;i++) tsize*=index[i];
    //   cout << "size=" << tsize << endl;
    // MEMORY -------------------
    aurostd::allocate_xtensor3corpus(corpus,lindex,uindex,index);
    // CLEAR ------------------
    for(int i=lindex[1];i<=uindex[1];i++)
      for(int j=lindex[2];j<=uindex[2];j++)
	for(int k=lindex[3];k<=uindex[3];k++)
	  corpus[i][j][k]=(utype) a.corpus[i][j][k]; //i,j,k
    // DONE ------------------
#ifdef _XMATH_DEBUG_CONSTRUCTORS
    cout << "iscubic=" << iscubic << ", isfloat=" << isfloat << ", iscomplex=" << iscomplex << ", sizeof=" << size << ", tsize=" << tsize << endl;
#endif
  }
}

// ----------------------------------------------- copy constructors xtensor4
namespace aurostd {  // namespace aurostd
  template<class utype>                                   // copy constructor
  xtensor4<utype>::xtensor4(const xtensor4<utype>& a) {
    for(int i=1;i<=_AUROSTD_XTENSOR4_INDEX_;i++) {lindex[i]=a.lindex[i];uindex[i]=a.uindex[i];index[i]=uindex[i]-lindex[i]+1;}
    iscubic=TRUE;
    for(int i=1;i<_AUROSTD_XTENSOR4_INDEX_;i++) {iscubic=iscubic && bool(index[i]==index[i+1]);};
    aurostd::xtensor4debug(*this,"copy constructor");
    isfloat=_isfloat((utype) 0);
    iscomplex=_iscomplex((utype) 0);
    size=(char) (sizeof(utype));
    tsize=(long int) size;for(int i=1;i<=_AUROSTD_XTENSOR4_INDEX_;i++) tsize*=index[i];
    //   cout << "size=" << tsize << endl;
    // MEMORY -------------------
    aurostd::allocate_xtensor4corpus(corpus,lindex,uindex,index);
    // CLEAR ------------------
    for(int i=lindex[1];i<=uindex[1];i++)
      for(int j=lindex[2];j<=uindex[2];j++)
	for(int k=lindex[3];k<=uindex[3];k++)
	  for(int l=lindex[4];l<=uindex[4];l++)
	    corpus[i][j][k][l]=(utype) a.corpus[i][j][k][l]; //i,j,k,l
    // DONE ------------------
#ifdef _XMATH_DEBUG_CONSTRUCTORS
    cout << "iscubic=" << iscubic << ", isfloat=" << isfloat << ", iscomplex=" << iscomplex << ", sizeof=" << size << ", tsize=" << tsize << endl;
#endif
  }
}

// ------------------------------------------------ copy constructors xtensor5
namespace aurostd {  // namespace aurostd
  template<class utype>                                    // copy constructor
  xtensor5<utype>::xtensor5(const xtensor5<utype>& a) {
    for(int i=1;i<=_AUROSTD_XTENSOR5_INDEX_;i++) {lindex[i]=a.lindex[i];uindex[i]=a.uindex[i];index[i]=uindex[i]-lindex[i]+1;}
    iscubic=TRUE;
    for(int i=1;i<_AUROSTD_XTENSOR5_INDEX_;i++) {iscubic=iscubic && bool(index[i]==index[i+1]);};
    aurostd::xtensor5debug(*this,"copy constructor");
    isfloat=_isfloat((utype) 0);
    iscomplex=_iscomplex((utype) 0);
    size=(char) (sizeof(utype));
    tsize=(long int) size;for(int i=1;i<=_AUROSTD_XTENSOR5_INDEX_;i++) tsize*=index[i];
    //   cout << "size=" << tsize << endl;
    // MEMORY -------------------
    aurostd::allocate_xtensor5corpus(corpus,lindex,uindex,index);
    // CLEAR ------------------
    for(int i=lindex[1];i<=uindex[1];i++)
      for(int j=lindex[2];j<=uindex[2];j++)
	for(int k=lindex[3];k<=uindex[3];k++)
	  for(int l=lindex[4];l<=uindex[4];l++)
	    for(int m=lindex[5];m<=uindex[5];m++)
	      corpus[i][j][k][l][m]=(utype) a.corpus[i][j][k][l][m]; //i,j,k,l,m
    // DONE ------------------
#ifdef _XMATH_DEBUG_CONSTRUCTORS
    cout << "iscubic=" << iscubic << ", isfloat=" << isfloat << ", iscomplex=" << iscomplex << ", sizeof=" << size << ", tsize=" << tsize << endl;
#endif
  }
}

// ------------------------------------------------ copy constructors xtensor6
namespace aurostd {  // namespace aurostd
  template<class utype>                                    // copy constructor
  xtensor6<utype>::xtensor6(const xtensor6<utype>& a) {
    for(int i=1;i<=_AUROSTD_XTENSOR6_INDEX_;i++) {lindex[i]=a.lindex[i];uindex[i]=a.uindex[i];index[i]=uindex[i]-lindex[i]+1;}
    iscubic=TRUE;
    for(int i=1;i<_AUROSTD_XTENSOR6_INDEX_;i++) {iscubic=iscubic && bool(index[i]==index[i+1]);};
    aurostd::xtensor6debug(*this,"copy constructor");
    isfloat=_isfloat((utype) 0);
    iscomplex=_iscomplex((utype) 0);
    size=(char) (sizeof(utype));
    tsize=(long int) size;for(int i=1;i<=_AUROSTD_XTENSOR6_INDEX_;i++) tsize*=index[i];
    //   cout << "size=" << tsize << endl;
    // MEMORY -------------------
    aurostd::allocate_xtensor6corpus(corpus,lindex,uindex,index);
    // CLEAR ------------------
    for(int i=lindex[1];i<=uindex[1];i++)
      for(int j=lindex[2];j<=uindex[2];j++)
	for(int k=lindex[3];k<=uindex[3];k++)
	  for(int l=lindex[4];l<=uindex[4];l++)
	    for(int m=lindex[5];m<=uindex[5];m++)
	      for(int n=lindex[6];n<=uindex[6];n++)
		corpus[i][j][k][l][m][n]=(utype) a.corpus[i][j][k][l][m][n]; //i,j,k,l,m,n
    // DONE ------------------
#ifdef _XMATH_DEBUG_CONSTRUCTORS
    cout << "iscubic=" << iscubic << ", isfloat=" << isfloat << ", iscomplex=" << iscomplex << ", sizeof=" << size << ", tsize=" << tsize << endl;
#endif
  }
}

// ------------------------------------------------ copy constructors xtensor7
namespace aurostd {  // namespace aurostd
  template<class utype>                                    // copy constructor
  xtensor7<utype>::xtensor7(const xtensor7<utype>& a) {
    for(int i=1;i<=_AUROSTD_XTENSOR7_INDEX_;i++) {lindex[i]=a.lindex[i];uindex[i]=a.uindex[i];index[i]=uindex[i]-lindex[i]+1;}
    iscubic=TRUE;
    for(int i=1;i<_AUROSTD_XTENSOR7_INDEX_;i++) {iscubic=iscubic && bool(index[i]==index[i+1]);};
    aurostd::xtensor7debug(*this,"copy constructor");
    isfloat=_isfloat((utype) 0);
    iscomplex=_iscomplex((utype) 0);
    size=(char) (sizeof(utype));
    tsize=(long int) size;for(int i=1;i<=_AUROSTD_XTENSOR7_INDEX_;i++) tsize*=index[i];
    //   cout << "size=" << tsize << endl;
    // MEMORY -------------------
    aurostd::allocate_xtensor7corpus(corpus,lindex,uindex,index);
    // CLEAR ------------------
    for(int i=lindex[1];i<=uindex[1];i++)
      for(int j=lindex[2];j<=uindex[2];j++)
	for(int k=lindex[3];k<=uindex[3];k++)
	  for(int l=lindex[4];l<=uindex[4];l++)
	    for(int m=lindex[5];m<=uindex[5];m++)
	      for(int n=lindex[6];n<=uindex[6];n++)
		for(int o=lindex[7];o<=uindex[7];o++)
		  corpus[i][j][k][l][m][n][o]=(utype) a.corpus[i][j][k][l][m][n][o]; //i,j,k,l,m,n
    // DONE ------------------
#ifdef _XMATH_DEBUG_CONSTRUCTORS
    cout << "iscubic=" << iscubic << ", isfloat=" << isfloat << ", iscomplex=" << iscomplex << ", sizeof=" << size << ", tsize=" << tsize << endl;
#endif
  }
}

// ------------------------------------------------ copy constructors xtensor8
namespace aurostd {  // namespace aurostd
  template<class utype>                                    // copy constructor
  xtensor8<utype>::xtensor8(const xtensor8<utype>& a) {
    for(int i=1;i<=_AUROSTD_XTENSOR8_INDEX_;i++) {lindex[i]=a.lindex[i];uindex[i]=a.uindex[i];index[i]=uindex[i]-lindex[i]+1;}
    iscubic=TRUE;
    for(int i=1;i<_AUROSTD_XTENSOR8_INDEX_;i++) {iscubic=iscubic && bool(index[i]==index[i+1]);};
    aurostd::xtensor8debug(*this,"copy constructor");
    isfloat=_isfloat((utype) 0);
    iscomplex=_iscomplex((utype) 0);
    size=(char) (sizeof(utype));
    tsize=(long int) size;for(int i=1;i<=_AUROSTD_XTENSOR8_INDEX_;i++) tsize*=index[i];
    //   cout << "size=" << tsize << endl;
    // MEMORY -------------------
    aurostd::allocate_xtensor8corpus(corpus,lindex,uindex,index);
    // CLEAR ------------------
    for(int i=lindex[1];i<=uindex[1];i++)
      for(int j=lindex[2];j<=uindex[2];j++)
	for(int k=lindex[3];k<=uindex[3];k++)
	  for(int l=lindex[4];l<=uindex[4];l++)
	    for(int m=lindex[5];m<=uindex[5];m++)
	      for(int n=lindex[6];n<=uindex[6];n++)
		for(int o=lindex[7];o<=uindex[7];o++)
		  for(int p=lindex[8];p<=uindex[8];p++)
		    corpus[i][j][k][l][m][n][o][p]=(utype) a.corpus[i][j][k][l][m][n][o][p]; //i,j,k,l,m,n,o,p
    // DONE ------------------
#ifdef _XMATH_DEBUG_CONSTRUCTORS
    cout << "iscubic=" << iscubic << ", isfloat=" << isfloat << ", iscomplex=" << iscomplex << ", sizeof=" << size << ", tsize=" << tsize << endl;
#endif
  }
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

// ***************************************************************************
// ------------------------------------------------------- destructor xtensor3
namespace aurostd {  // namespace aurostd
  template<class utype>                                  // default destructor
  xtensor3<utype>::~xtensor3(void) {
    // should work, but be careful
    // free a xtensor3 allocated with xtensor3()
    aurostd::xtensor3debug(*this,"default destructor");
    for(int i=lindex[1];i<=uindex[1];i++) {
      delete [] (corpus[i][lindex[2]]+lindex[3]-XXEND);
      delete [] (corpus[i]+lindex[2]-XXEND);
    }
    delete [] (corpus+lindex[1]-XXEND);
  }
}

// ------------------------------------------------------- destructor xtensor4
namespace aurostd {  // namespace aurostd
  template<class utype>                                  // default destructor
  xtensor4<utype>::~xtensor4(void) {
    // should work, but be careful
    // free a xtensor4 allocated with xtensor4()
    aurostd::xtensor4debug(*this,"default destructor");
    for(int i=lindex[1];i<=uindex[1];i++) {
      for(int j=lindex[2];j<=uindex[2];j++) {
	delete [] (corpus[i][j][lindex[3]]+lindex[4]-XXEND);
	delete [] (corpus[i][j]+lindex[3]-XXEND);
      }
      delete [] (corpus[i]+lindex[2]-XXEND);
    }
    delete [] (corpus+lindex[1]-XXEND);
  }
}

// ------------------------------------------------------- destructor xtensor5
namespace aurostd {  // namespace aurostd
  template<class utype>                                  // default destructor
  xtensor5<utype>::~xtensor5(void) {
    // should work, but be careful
    // free a xtensor5 allocated with xtensor5()
    aurostd::xtensor5debug(*this,"default destructor");
    for(int i=lindex[1];i<=uindex[1];i++) {
      for(int j=lindex[2];j<=uindex[2];j++) {
	for(int k=lindex[3];k<=uindex[3];k++) {
	  delete [] (corpus[i][j][k][lindex[4]]+lindex[5]-XXEND);
	  delete [] (corpus[i][j][k]+lindex[4]-XXEND);
	}
	delete [] (corpus[i][j]+lindex[3]-XXEND);
      }
      delete [] (corpus[i]+lindex[2]-XXEND);
    }
    delete [] (corpus+lindex[1]-XXEND);
  }
}

// ------------------------------------------------------- destructor xtensor6
namespace aurostd {  // namespace aurostd
  template<class utype>                                  // default destructor
  xtensor6<utype>::~xtensor6(void) {
    // should work, but be careful
    // free a xtensor6 allocated with xtensor6()
    aurostd::xtensor6debug(*this,"default destructor");
    for(int i=lindex[1];i<=uindex[1];i++) {
      for(int j=lindex[2];j<=uindex[2];j++) {
	for(int k=lindex[3];k<=uindex[3];k++) {
	  for(int l=lindex[4];l<=uindex[4];l++) {
	    delete [] (corpus[i][j][k][l][lindex[5]]+lindex[6]-XXEND);
	    delete [] (corpus[i][j][k][l]+lindex[5]-XXEND);
	  }
	  delete [] (corpus[i][j][k]+lindex[4]-XXEND);
	}
	delete [] (corpus[i][j]+lindex[3]-XXEND);
      }
      delete [] (corpus[i]+lindex[2]-XXEND);
    }
    delete [] (corpus+lindex[1]-XXEND);
  }
}

// ------------------------------------------------------- destructor xtensor7
namespace aurostd {  // namespace aurostd
  template<class utype>                                  // default destructor
  xtensor7<utype>::~xtensor7(void) {
    // should work, but be careful
    // free a xtensor7 allocated with xtensor7()
    aurostd::xtensor7debug(*this,"default destructor");
    for(int i=lindex[1];i<=uindex[1];i++) {
      for(int j=lindex[2];j<=uindex[2];j++) {
	for(int k=lindex[3];k<=uindex[3];k++) {
	  for(int l=lindex[4];l<=uindex[4];l++) {
	    for(int m=lindex[5];m<=uindex[5];m++) {
	      delete [] (corpus[i][j][k][l][m][lindex[6]]+lindex[6]-XXEND);
	      delete [] (corpus[i][j][k][l][m]+lindex[6]-XXEND);
	    }
	    delete [] (corpus[i][j][k][l]+lindex[5]-XXEND);
	  }
	  delete [] (corpus[i][j][k]+lindex[4]-XXEND);
	}
	delete [] (corpus[i][j]+lindex[3]-XXEND);
      }
      delete [] (corpus[i]+lindex[2]-XXEND);
    }
    delete [] (corpus+lindex[1]-XXEND);
  }
}

// ------------------------------------------------------- destructor xtensor8
namespace aurostd {  // namespace aurostd
  template<class utype>                                  // default destructor
  xtensor8<utype>::~xtensor8(void) {
    // should work, but be careful
    // free a xtensor8 allocated with xtensor8()
    aurostd::xtensor8debug(*this,"default destructor");
    for(int i=lindex[1];i<=uindex[1];i++) {
      for(int j=lindex[2];j<=uindex[2];j++) {
	for(int k=lindex[3];k<=uindex[3];k++) {
	  for(int l=lindex[4];l<=uindex[4];l++) {
	    for(int m=lindex[5];m<=uindex[5];m++) {
	      for(int n=lindex[6];n<=uindex[6];n++) {
		delete [] (corpus[i][j][k][l][m][n][lindex[7]]+lindex[7]-XXEND);
		delete [] (corpus[i][j][k][l][m][n]+lindex[7]-XXEND);
	      }
	      delete [] (corpus[i][j][k][l][m]+lindex[6]-XXEND);
	    }
	    delete [] (corpus[i][j][k][l]+lindex[5]-XXEND);
	  }
	  delete [] (corpus[i][j][k]+lindex[4]-XXEND);
	}
	delete [] (corpus[i][j]+lindex[3]-XXEND);
      }
      delete [] (corpus[i]+lindex[2]-XXEND);
    }
    delete [] (corpus+lindex[1]-XXEND);
  }
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ------------------------------------------------------- assigment operators

// ***************************************************************************
// -------------------------------------------------------- operator= xtensor3
namespace aurostd {  // namespace aurostd
  template<class utype>                                         // operator =
  xtensor3<utype>& xtensor3<utype>::operator=(const xtensor3<utype>& a) {
    if(corpus!=a.corpus) {                // check  for a=a
      if(index[1]!=a.index[1] || index[2]!=a.index[2] || index[3]!=a.index[3]) {
	// if dims(this)!=dims(a) => build a new xtensor3 !!!
	// destroy
	for(int i=lindex[1];i<=uindex[1];i++) {
	  delete [] (corpus[i][lindex[2]]+lindex[3]-XXEND);
	  delete [] (corpus[i]+lindex[2]-XXEND);
	}
	delete [] (corpus+lindex[1]-XXEND);
	// build
	for(int i=1;i<=_AUROSTD_XTENSOR3_INDEX_;i++) {index[i]=a.index[i];uindex[i]=a.uindex[i];lindex[i]=a.lindex[i];}
	iscubic=TRUE;
	for(int i=1;i<_AUROSTD_XTENSOR3_INDEX_;i++) {iscubic=iscubic && bool(index[i]==index[i+1]);};
	aurostd::xtensor3debug(*this,"operator =:");
	// MEMORY -------------------
	aurostd::allocate_xtensor3corpus(corpus,lindex,uindex,index);
	// DONE ------------------
	isfloat=_isfloat((utype) 0);
	iscomplex=_iscomplex((utype) 0);
	size=(char) sizeof(utype);
	tsize=(long int) size;for(int i=1;i<=_AUROSTD_XTENSOR3_INDEX_;i++) tsize*=index[i];
#ifdef _XMATH_DEBUG_CONSTRUCTORS
	cout << "iscubic=" << iscubic << ", isfloat=" << isfloat << ", iscomplex=" << iscomplex << ", sizeof=" << size << ", tsize=" << tsize << endl;
#endif
      }
      for(int k=0;k<index[3];k++)
	for(int j=0;j<index[2];j++)
	  for(int i=0;i<index[1];i++)
	    this->corpus[i+lindex[1]][j+lindex[2]][k+lindex[3]]=
	      a.corpus[i+a.lindex[1]][j+a.lindex[2]][k+a.lindex[3]];
    }
    return *this;
  }
}

// -------------------------------------------------------- operator= xtensor4
namespace aurostd {  // namespace aurostd
  template<class utype>                                         // operator =
  xtensor4<utype>& xtensor4<utype>::operator=(const xtensor4<utype>& a) {
    if(corpus!=a.corpus) {                // check  for a=a
      if(index[1]!=a.index[1] || index[2]!=a.index[2] || index[3]!=a.index[3] ||
	 index[4]!=a.index[4]) {
	// if dims(this)!=dims(a) => build a new xtensor4 !!!
	// destroy
	for(int i=lindex[1];i<=uindex[1];i++) {
	  for(int j=lindex[2];j<=uindex[2];j++) {
	    delete [] (corpus[i][j][lindex[3]]+lindex[4]-XXEND);
	    delete [] (corpus[i][j]+lindex[3]-XXEND);
	  }
	  delete [] (corpus[i]+lindex[2]-XXEND);
	}
	delete [] (corpus+lindex[1]-XXEND);
	// build
	for(int i=1;i<=_AUROSTD_XTENSOR4_INDEX_;i++) {index[i]=a.index[i];uindex[i]=a.uindex[i];lindex[i]=a.lindex[i];}
	iscubic=TRUE;
	for(int i=1;i<_AUROSTD_XTENSOR4_INDEX_;i++) {iscubic=iscubic && bool(index[i]==index[i+1]);};
	aurostd::xtensor4debug(*this,"operator =:");

	// MEMORY -------------------
	aurostd::allocate_xtensor4corpus(corpus,lindex,uindex,index);
	// DONE ------------------
	isfloat=_isfloat((utype) 0);
	iscomplex=_iscomplex((utype) 0);
	size=(char) sizeof(utype);
	tsize=(long int) size;for(int i=1;i<=_AUROSTD_XTENSOR4_INDEX_;i++) tsize*=index[i];
#ifdef _XMATH_DEBUG_CONSTRUCTORS
	cout << "iscubic=" << iscubic << ", isfloat=" << isfloat << ", iscomplex=" << iscomplex << ", sizeof=" << size << ", tsize=" << tsize << endl;
#endif
      }
      for(int l=0;l<index[4];l++)
	for(int k=0;k<index[3];k++)
	  for(int j=0;j<index[2];j++)
	    for(int i=0;i<index[1];i++)
	      this->corpus[i+lindex[1]][j+lindex[2]][k+lindex[3]][l+lindex[4]]=
		a.corpus[i+a.lindex[1]][j+a.lindex[2]][k+a.lindex[3]][l+a.lindex[4]];
    }
    return *this;
  }
}

// -------------------------------------------------------- operator= xtensor5
namespace aurostd {  // namespace aurostd
  template<class utype>                                          // operator =
  xtensor5<utype>& xtensor5<utype>::operator=(const xtensor5<utype>& a) {
    if(corpus!=a.corpus) {                // check  for a=a
      if(index[1]!=a.index[1] || index[2]!=a.index[2] || index[3]!=a.index[3] ||
	 index[4]!=a.index[4] || index[5]!=a.index[5]) {
	// if dims(this)!=dims(a) => build a new xtensor5 !!!
	// destroy
	for(int i=lindex[1];i<=uindex[1];i++) {
	  for(int j=lindex[2];j<=uindex[2];j++) {
	    for(int k=lindex[3];k<=uindex[3];k++) {
	      delete [] (corpus[i][j][k][lindex[4]]+lindex[5]-XXEND);
	      delete [] (corpus[i][j][k]+lindex[4]-XXEND);
	    }
	    delete [] (corpus[i][j]+lindex[3]-XXEND);
	  }
	  delete [] (corpus[i]+lindex[2]-XXEND);
	}
	delete [] (corpus+lindex[1]-XXEND);
	// build
	for(int i=1;i<=_AUROSTD_XTENSOR5_INDEX_;i++) {index[i]=a.index[i];uindex[i]=a.uindex[i];lindex[i]=a.lindex[i];}
	iscubic=TRUE;
	for(int i=1;i<_AUROSTD_XTENSOR5_INDEX_;i++) {iscubic=iscubic && bool(index[i]==index[i+1]);};
	aurostd::xtensor5debug(*this,"operator =:");

	// MEMORY -------------------
	aurostd::allocate_xtensor5corpus(corpus,lindex,uindex,index);
	// DONE ------------------
	
	isfloat=_isfloat((utype) 0);
	iscomplex=_iscomplex((utype) 0);
	size=(char) sizeof(utype);
	tsize=(long int) size;for(int i=1;i<=_AUROSTD_XTENSOR5_INDEX_;i++) tsize*=index[i];
#ifdef _XMATH_DEBUG_CONSTRUCTORS
	cout << "iscubic=" << iscubic << ", isfloat=" << isfloat << ", iscomplex=" << iscomplex << ", sizeof=" << size << ", tsize=" << tsize << endl;
#endif
      }
      for(int m=0;m<index[5];m++)
	for(int l=0;l<index[4];l++)
	  for(int k=0;k<index[3];k++)
	    for(int j=0;j<index[2];j++)
	      for(int i=0;i<index[1];i++)
		this->corpus[i+lindex[1]][j+lindex[2]][k+lindex[3]][l+lindex[4]][m+lindex[5]]=
		  a.corpus[i+a.lindex[1]][j+a.lindex[2]][k+a.lindex[3]][l+a.lindex[4]][m+a.lindex[5]];
    }
    return *this;
  }
}

// -------------------------------------------------------- operator= xtensor6
namespace aurostd {  // namespace aurostd
  template<class utype>                                          // operator =
  xtensor6<utype>& xtensor6<utype>::operator=(const xtensor6<utype>& a) {
    if(corpus!=a.corpus) {                // check  for a=a
      if(index[1]!=a.index[1] || index[2]!=a.index[2] || index[3]!=a.index[3] ||
	 index[4]!=a.index[4] || index[5]!=a.index[5] || index[6]!=a.index[6]) {
	// if dims(this)!=dims(a) => build a new xtensor6 !!!
	// destroy
	for(int i=lindex[1];i<=uindex[1];i++) {
	  for(int j=lindex[2];j<=uindex[2];j++) {
	    for(int k=lindex[3];k<=uindex[3];k++) {
	      for(int l=lindex[4];l<=uindex[4];l++) {
		delete [] (corpus[i][j][k][l][lindex[5]]+lindex[6]-XXEND);
		delete [] (corpus[i][j][k][l]+lindex[5]-XXEND);
	      }
	      delete [] (corpus[i][j][k]+lindex[4]-XXEND);
	    }
	    delete [] (corpus[i][j]+lindex[3]-XXEND);
	  }
	  delete [] (corpus[i]+lindex[2]-XXEND);
	}
	delete [] (corpus+lindex[1]-XXEND);
	// build
	for(int i=1;i<=_AUROSTD_XTENSOR6_INDEX_;i++) {index[i]=a.index[i];uindex[i]=a.uindex[i];lindex[i]=a.lindex[i];}
	iscubic=TRUE;
	for(int i=1;i<_AUROSTD_XTENSOR6_INDEX_;i++) {iscubic=iscubic && bool(index[i]==index[i+1]);};
	aurostd::xtensor6debug(*this,"operator =:");

	// MEMORY -------------------
	aurostd::allocate_xtensor6corpus(corpus,lindex,uindex,index);
	// DONE ------------------
	
	isfloat=_isfloat((utype) 0);
	iscomplex=_iscomplex((utype) 0);
	size=(char) sizeof(utype);
	tsize=(long int) size;for(int i=1;i<=_AUROSTD_XTENSOR6_INDEX_;i++) tsize*=index[i];
#ifdef _XMATH_DEBUG_CONSTRUCTORS
	cout << "iscubic=" << iscubic << ", isfloat=" << isfloat << ", iscomplex=" << iscomplex << ", sizeof=" << size << ", tsize=" << tsize << endl;
#endif
      }
      for(int n=0;n<index[6];n++)
	for(int m=0;m<index[5];m++)
	  for(int l=0;l<index[4];l++)
	    for(int k=0;k<index[3];k++)
	      for(int j=0;j<index[2];j++)
		for(int i=0;i<index[1];i++)
		  this->corpus[i+lindex[1]][j+lindex[2]][k+lindex[3]][l+lindex[4]][m+lindex[5]][n+lindex[6]]=
		    a.corpus[i+a.lindex[1]][j+a.lindex[2]][k+a.lindex[3]][l+a.lindex[4]][m+a.lindex[5]][n+a.lindex[6]];
    }
    return *this;
  }
}

// -------------------------------------------------------- operator= xtensor7
namespace aurostd {  // namespace aurostd
  template<class utype>                                          // operator =
  xtensor7<utype>& xtensor7<utype>::operator=(const xtensor7<utype>& a) {
    if(corpus!=a.corpus) {                // check  for a=a
      if(index[1]!=a.index[1] || index[2]!=a.index[2] || index[3]!=a.index[3] ||
	 index[4]!=a.index[4] || index[5]!=a.index[5] || index[6]!=a.index[6] || index[7]!=a.index[7]) {
	// if dims(this)!=dims(a) => build a new xtensor7 !!!
	// destroy
	for(int i=lindex[1];i<=uindex[1];i++) {
	  for(int j=lindex[2];j<=uindex[2];j++) {
	    for(int k=lindex[3];k<=uindex[3];k++) {
	      for(int l=lindex[4];l<=uindex[4];l++) {
		for(int m=lindex[5];m<=uindex[5];m++) {
		  delete [] (corpus[i][j][k][l][m][lindex[6]]+lindex[6]-XXEND);
		  delete [] (corpus[i][j][k][l][m]+lindex[6]-XXEND);
		}
		delete [] (corpus[i][j][k][l]+lindex[5]-XXEND);
	      }
	      delete [] (corpus[i][j][k]+lindex[4]-XXEND);
	    }
	    delete [] (corpus[i][j]+lindex[3]-XXEND);
	  }
	  delete [] (corpus[i]+lindex[2]-XXEND);
	}
	delete [] (corpus+lindex[1]-XXEND);
      	// build
	for(int i=1;i<=_AUROSTD_XTENSOR7_INDEX_;i++) {index[i]=a.index[i];uindex[i]=a.uindex[i];lindex[i]=a.lindex[i];}
	iscubic=TRUE;
	for(int i=1;i<_AUROSTD_XTENSOR7_INDEX_;i++) {iscubic=iscubic && bool(index[i]==index[i+1]);};
	aurostd::xtensor7debug(*this,"operator =:");
	
	// MEMORY -------------------
	aurostd::allocate_xtensor7corpus(corpus,lindex,uindex,index);
	// DONE ------------------
	
	isfloat=_isfloat((utype) 0);
	iscomplex=_iscomplex((utype) 0);
	size=(char) sizeof(utype);
	tsize=(long int) size;for(int i=1;i<=_AUROSTD_XTENSOR7_INDEX_;i++) tsize*=index[i];
#ifdef _XMATH_DEBUG_CONSTRUCTORS
	cout << "iscubic=" << iscubic << ", isfloat=" << isfloat << ", iscomplex=" << iscomplex << ", sizeof=" << size << ", tsize=" << tsize << endl;
#endif
      }
      for(int o=0;o<index[7];o++)
	for(int n=0;n<index[6];n++)
	  for(int m=0;m<index[5];m++)
	    for(int l=0;l<index[4];l++)
	      for(int k=0;k<index[3];k++)
		for(int j=0;j<index[2];j++)
		  for(int i=0;i<index[1];i++)
		    this->corpus[i+lindex[1]][j+lindex[2]][k+lindex[3]][l+lindex[4]][m+lindex[5]][n+lindex[6]][o+lindex[7]]=
		      a.corpus[i+a.lindex[1]][j+a.lindex[2]][k+a.lindex[3]][l+a.lindex[4]][m+a.lindex[5]][n+a.lindex[6]][o+a.lindex[7]];
    }
    return *this;
  }
}

// -------------------------------------------------------- operator= xtensor8
namespace aurostd {  // namespace aurostd
  template<class utype>                                          // operator =
  xtensor8<utype>& xtensor8<utype>::operator=(const xtensor8<utype>& a) {
    if(corpus!=a.corpus) {                // check  for a=a
      if(index[1]!=a.index[1] || index[2]!=a.index[2] || index[3]!=a.index[3] || index[4]!=a.index[4] ||
	 index[5]!=a.index[5] || index[6]!=a.index[6] || index[7]!=a.index[7] || index[8]!=a.index[8]) {
	// if dims(this)!=dims(a) => build a new xtensor8 !!!
	// destroy
	for(int i=lindex[1];i<=uindex[1];i++) {
	  for(int j=lindex[2];j<=uindex[2];j++) {
	    for(int k=lindex[3];k<=uindex[3];k++) {
	      for(int l=lindex[4];l<=uindex[4];l++) {
		for(int m=lindex[5];m<=uindex[5];m++) {
		  for(int n=lindex[6];n<=uindex[6];n++) {
		    delete [] (corpus[i][j][k][l][m][n][lindex[7]]+lindex[7]-XXEND);
		    delete [] (corpus[i][j][k][l][m][n]+lindex[7]-XXEND);
		  }
		  delete [] (corpus[i][j][k][l][m]+lindex[6]-XXEND);
		}
		delete [] (corpus[i][j][k][l]+lindex[5]-XXEND);
	      }
	      delete [] (corpus[i][j][k]+lindex[4]-XXEND);
	    }
	    delete [] (corpus[i][j]+lindex[3]-XXEND);
	  }
	  delete [] (corpus[i]+lindex[2]-XXEND);
	}
	delete [] (corpus+lindex[1]-XXEND);
	// build
	for(int i=1;i<=_AUROSTD_XTENSOR8_INDEX_;i++) {index[i]=a.index[i];uindex[i]=a.uindex[i];lindex[i]=a.lindex[i];}
	iscubic=TRUE;
	for(int i=1;i<_AUROSTD_XTENSOR8_INDEX_;i++) {iscubic=iscubic && bool(index[i]==index[i+1]);};
	aurostd::xtensor8debug(*this,"operator =:");
	
	// MEMORY -------------------
	aurostd::allocate_xtensor8corpus(corpus,lindex,uindex,index);
	// DONE ------------------
	
	isfloat=_isfloat((utype) 0);
	iscomplex=_iscomplex((utype) 0);
	size=(char) sizeof(utype);
	tsize=(long int) size;for(int i=1;i<=_AUROSTD_XTENSOR8_INDEX_;i++) tsize*=index[i];
#ifdef _XMATH_DEBUG_CONSTRUCTORS
	cout << "iscubic=" << iscubic << ", isfloat=" << isfloat << ", iscomplex=" << iscomplex << ", sizeof=" << size << ", tsize=" << tsize << endl;
#endif
      }
      for(int p=0;p<index[8];p++)
	for(int o=0;o<index[7];o++)
	  for(int n=0;n<index[6];n++)
	    for(int m=0;m<index[5];m++)
	      for(int l=0;l<index[4];l++)
		for(int k=0;k<index[3];k++)
		  for(int j=0;j<index[2];j++)
		    for(int i=0;i<index[1];i++)
		      this->corpus[i+lindex[1]][j+lindex[2]][k+lindex[3]][l+lindex[4]][m+lindex[5]][n+lindex[6]][o+lindex[7]][p+lindex[8]]=
			a.corpus[i+a.lindex[1]][j+a.lindex[2]][k+a.lindex[3]][l+a.lindex[4]][m+a.lindex[5]][n+a.lindex[6]][o+a.lindex[7]][p+a.lindex[8]];
    }
    return *this;
  }
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ----------------------------------------------------------- index operators

// ***************************************************************************
// ------------------------------------------------------ operator [] xtensor3
namespace aurostd {  // namespace aurostd
  template<class utype>
  inline utype** xtensor3<utype>::operator[] (int i) const {
#ifndef __OPTIMIZE
    if(i>uindex[1]) {
      cerr << "_xtensor3<utype>_index[1]_high i=" << i << ", lindex[1]=" << lindex[1] << ", hindex[1]=" << uindex[1] << endl;
      exit(0);
    }
    if(i<lindex[1]) {
      cerr << "_xtensor3<utype>_index[1]_low  i=" << i << ", lindex[1]=" << lindex[1] << ", hindex[1]=" << uindex[1] << endl;
      exit(0);
    }
#endif
    return corpus[i];
  }
}

// ------------------------------------------------------ operator [] xtensor4
namespace aurostd {  // namespace aurostd
  template<class utype>
  inline utype*** xtensor4<utype>::operator[] (int i) const {
#ifndef __OPTIMIZE
    if(i>uindex[1]) {
      cerr << "_xtensor4<utype>_index[1]_high i=" << i << ", lindex[1]=" << lindex[1] << ", hindex[1]=" << uindex[1] << endl;
      exit(0);
    }
    if(i<lindex[1]) {
      cerr << "_xtensor4<utype>_index[1]_low  i=" << i << ", lindex[1]=" << lindex[1] << ", hindex[1]=" << uindex[1] << endl;
      exit(0);
    }
#endif
    return corpus[i];
  }
}

// ------------------------------------------------------ operator [] xtensor5
namespace aurostd {  // namespace aurostd
  template<class utype>
  inline utype**** xtensor5<utype>::operator[] (int i) const {
#ifndef __OPTIMIZE
    if(i>uindex[1]) {
      cerr << "_xtensor5<utype>_index[1]_high i=" << i << ", lindex[1]=" << lindex[1] << ", hindex[1]=" << uindex[1] << endl;
      exit(0);
    }
    if(i<lindex[1]) {
      cerr << "_xtensor5<utype>_index[1]_low  i=" << i << ", lindex[1]=" << lindex[1] << ", hindex[1]=" << uindex[1] << endl;
      exit(0);
    }
#endif
    return corpus[i];
  }
}

// ------------------------------------------------------ operator [] xtensor6
namespace aurostd {  // namespace aurostd
  template<class utype>
  inline utype***** xtensor6<utype>::operator[] (int i) const {
#ifndef __OPTIMIZE
    if(i>uindex[1]) {
      cerr << "_xtensor6<utype>_index[1]_high i=" << i << ", lindex[1]=" << lindex[1] << ", hindex[1]=" << uindex[1] << endl;
      exit(0);
    }
    if(i<lindex[1]) {
      cerr << "_xtensor6<utype>_index[1]_low  i=" << i << ", lindex[1]=" << lindex[1] << ", hindex[1]=" << uindex[1] << endl;
      exit(0);
    }
#endif
    return corpus[i];
  }
}

// ------------------------------------------------------ operator [] xtensor7
namespace aurostd {  // namespace aurostd
  template<class utype>
  inline utype****** xtensor7<utype>::operator[] (int i) const {
#ifndef __OPTIMIZE
    if(i>uindex[1]) {
      cerr << "_xtensor7<utype>_index[1]_high i=" << i << ", lindex[1]=" << lindex[1] << ", hindex[1]=" << uindex[1] << endl;
      exit(0);
    }
    if(i<lindex[1]) {
      cerr << "_xtensor7<utype>_index[1]_low  i=" << i << ", lindex[1]=" << lindex[1] << ", hindex[1]=" << uindex[1] << endl;
      exit(0);
    }
#endif
    return corpus[i];
  }
}

// ------------------------------------------------------ operator [] xtensor8
namespace aurostd {  // namespace aurostd
  template<class utype>
  inline utype******* xtensor8<utype>::operator[] (int i) const {
#ifndef __OPTIMIZE
    if(i>uindex[1]) {
      cerr << "_xtensor8<utype>_index[1]_high i=" << i << ", lindex[1]=" << lindex[1] << ", hindex[1]=" << uindex[1] << endl;
      exit(0);
    }
    if(i<lindex[1]) {
      cerr << "_xtensor8<utype>_index[1]_low  i=" << i << ", lindex[1]=" << lindex[1] << ", hindex[1]=" << uindex[1] << endl;
      exit(0);
    }
#endif
    return corpus[i];
  }
}

// ***************************************************************************
// ------------------------------------------------------ operator () xtensor3
namespace aurostd {  // namespace aurostd
  template<class utype>                                         // operator ()
  inline utype& xtensor3<utype>::operator()(int i,int j,int k) const {
#ifndef __OPTIMIZE
    if(i>uindex[1]) {cerr << "XTENSOR3 -> i=" << i << " > uindex[1]=" << uindex[1] << endl;exit(0);}
    if(i<lindex[1]) {cerr << "XTENSOR3 -> i=" << i << " < lindex[1]=" << lindex[1] << endl;exit(0);}
    if(j>uindex[2]) {cerr << "XTENSOR3 -> j=" << j << " > uindex[2]=" << uindex[2] << endl;exit(0);}
    if(j<lindex[2]) {cerr << "XTENSOR3 -> j=" << j << " < lindex[2]=" << lindex[2] << endl;exit(0);}
    if(k>uindex[3]) {cerr << "XTENSOR3 -> k=" << k << " > uindex[3]=" << uindex[3] << endl;exit(0);}
    if(k<lindex[3]) {cerr << "XTENSOR3 -> k=" << k << " < lindex[3]=" << lindex[3] << endl;exit(0);}
#endif
    return corpus[i][j][k];
  }
}

// ------------------------------------------------------ operator () xtensor4
namespace aurostd {  // namespace aurostd
  template<class utype>                                         // operator ()
  inline utype& xtensor4<utype>::operator()(int i,int j,int k,int l) const {
#ifndef __OPTIMIZE
    if(i>uindex[1]) {cerr << "XTENSOR4 -> i=" << i << " > uindex[1]=" << uindex[1] << endl;exit(0);}
    if(i<lindex[1]) {cerr << "XTENSOR4 -> i=" << i << " < lindex[1]=" << lindex[1] << endl;exit(0);}
    if(j>uindex[2]) {cerr << "XTENSOR4 -> j=" << j << " > uindex[2]=" << uindex[2] << endl;exit(0);}
    if(j<lindex[2]) {cerr << "XTENSOR4 -> j=" << j << " < lindex[2]=" << lindex[2] << endl;exit(0);}
    if(k>uindex[3]) {cerr << "XTENSOR4 -> k=" << k << " > uindex[3]=" << uindex[3] << endl;exit(0);}
    if(k<lindex[3]) {cerr << "XTENSOR4 -> k=" << k << " < lindex[3]=" << lindex[3] << endl;exit(0);}
    if(l>uindex[4]) {cerr << "XTENSOR4 -> l=" << l << " > uindex[4]=" << uindex[4] << endl;exit(0);}
    if(l<lindex[4]) {cerr << "XTENSOR4 -> l=" << l << " < lindex[4]=" << lindex[4] << endl;exit(0);}
#endif
    return corpus[i][j][k][l];
  }
}

// ------------------------------------------------------ operator () xtensor5
namespace aurostd {  // namespace aurostd
  template<class utype>                                         // operator ()
  inline utype& xtensor5<utype>::operator()(int i,int j,int k,int l,int m) const {
#ifndef __OPTIMIZE
    if(i>uindex[1]) {cerr << "XTENSOR5 -> i=" << i << " > uindex[1]=" << uindex[1] << endl;exit(0);}
    if(i<lindex[1]) {cerr << "XTENSOR5 -> i=" << i << " < lindex[1]=" << lindex[1] << endl;exit(0);}
    if(j>uindex[2]) {cerr << "XTENSOR5 -> j=" << j << " > uindex[2]=" << uindex[2] << endl;exit(0);}
    if(j<lindex[2]) {cerr << "XTENSOR5 -> j=" << j << " < lindex[2]=" << lindex[2] << endl;exit(0);}
    if(k>uindex[3]) {cerr << "XTENSOR5 -> k=" << k << " > uindex[3]=" << uindex[3] << endl;exit(0);}
    if(k<lindex[3]) {cerr << "XTENSOR5 -> k=" << k << " < lindex[3]=" << lindex[3] << endl;exit(0);}
    if(l>uindex[4]) {cerr << "XTENSOR5 -> l=" << l << " > uindex[4]=" << uindex[4] << endl;exit(0);}
    if(l<lindex[4]) {cerr << "XTENSOR5 -> l=" << l << " < lindex[4]=" << lindex[4] << endl;exit(0);}
    if(m>uindex[5]) {cerr << "XTENSOR5 -> m=" << m << " > uindex[5]=" << uindex[5] << endl;exit(0);}
    if(m<lindex[5]) {cerr << "XTENSOR5 -> m=" << m << " < lindex[5]=" << lindex[5] << endl;exit(0);}
#endif
    return corpus[i][j][k][l][m];
  }
}

// ------------------------------------------------------ operator () xtensor6
namespace aurostd {  // namespace aurostd
  template<class utype>                                         // operator ()
  inline utype& xtensor6<utype>::operator()(int i,int j,int k,int l,int m,int n) const {
#ifndef __OPTIMIZE
    if(i>uindex[1]) {cerr << "XTENSOR6 -> i=" << i << " > uindex[1]=" << uindex[1] << endl;exit(0);}
    if(i<lindex[1]) {cerr << "XTENSOR6 -> i=" << i << " < lindex[1]=" << lindex[1] << endl;exit(0);}
    if(j>uindex[2]) {cerr << "XTENSOR6 -> j=" << j << " > uindex[2]=" << uindex[2] << endl;exit(0);}
    if(j<lindex[2]) {cerr << "XTENSOR6 -> j=" << j << " < lindex[2]=" << lindex[2] << endl;exit(0);}
    if(k>uindex[3]) {cerr << "XTENSOR6 -> k=" << k << " > uindex[3]=" << uindex[3] << endl;exit(0);}
    if(k<lindex[3]) {cerr << "XTENSOR6 -> k=" << k << " < lindex[3]=" << lindex[3] << endl;exit(0);}
    if(l>uindex[4]) {cerr << "XTENSOR6 -> l=" << l << " > uindex[4]=" << uindex[4] << endl;exit(0);}
    if(l<lindex[4]) {cerr << "XTENSOR6 -> l=" << l << " < lindex[4]=" << lindex[4] << endl;exit(0);}
    if(m>uindex[5]) {cerr << "XTENSOR6 -> m=" << m << " > uindex[5]=" << uindex[5] << endl;exit(0);}
    if(m<lindex[5]) {cerr << "XTENSOR6 -> m=" << m << " < lindex[5]=" << lindex[5] << endl;exit(0);}
    if(n>uindex[6]) {cerr << "XTENSOR6 -> n=" << n << " > uindex[6]=" << uindex[6] << endl;exit(0);}
    if(n<lindex[6]) {cerr << "XTENSOR6 -> n=" << n << " < lindex[6]=" << lindex[6] << endl;exit(0);}
#endif
    return corpus[i][j][k][l][m][n];
  }
}

// ------------------------------------------------------ operator () xtensor7
namespace aurostd {  // namespace aurostd
  template<class utype>                                         // operator ()
  inline utype& xtensor7<utype>::operator()(int i,int j,int k,int l,int m,int n,int o) const {
#ifndef __OPTIMIZE
    if(i>uindex[1]) {cerr << "XTENSOR7 -> i=" << i << " > uindex[1]=" << uindex[1] << endl;exit(0);}
    if(i<lindex[1]) {cerr << "XTENSOR7 -> i=" << i << " < lindex[1]=" << lindex[1] << endl;exit(0);}
    if(j>uindex[2]) {cerr << "XTENSOR7 -> j=" << j << " > uindex[2]=" << uindex[2] << endl;exit(0);}
    if(j<lindex[2]) {cerr << "XTENSOR7 -> j=" << j << " < lindex[2]=" << lindex[2] << endl;exit(0);}
    if(k>uindex[3]) {cerr << "XTENSOR7 -> k=" << k << " > uindex[3]=" << uindex[3] << endl;exit(0);}
    if(k<lindex[3]) {cerr << "XTENSOR7 -> k=" << k << " < lindex[3]=" << lindex[3] << endl;exit(0);}
    if(l>uindex[4]) {cerr << "XTENSOR7 -> l=" << l << " > uindex[4]=" << uindex[4] << endl;exit(0);}
    if(l<lindex[4]) {cerr << "XTENSOR7 -> l=" << l << " < lindex[4]=" << lindex[4] << endl;exit(0);}
    if(m>uindex[5]) {cerr << "XTENSOR7 -> m=" << m << " > uindex[5]=" << uindex[5] << endl;exit(0);}
    if(m<lindex[5]) {cerr << "XTENSOR7 -> m=" << m << " < lindex[5]=" << lindex[5] << endl;exit(0);}
    if(n>uindex[6]) {cerr << "XTENSOR7 -> n=" << n << " > uindex[6]=" << uindex[6] << endl;exit(0);}
    if(n<lindex[6]) {cerr << "XTENSOR7 -> n=" << n << " < lindex[6]=" << lindex[6] << endl;exit(0);}
    if(o>uindex[7]) {cerr << "XTENSOR7 -> o=" << o << " > uindex[7]=" << uindex[7] << endl;exit(0);}
    if(o<lindex[7]) {cerr << "XTENSOR7 -> o=" << o << " < lindex[7]=" << lindex[7] << endl;exit(0);}
#endif
    return corpus[i][j][k][l][m][n][o];
  }
}

// ------------------------------------------------------ operator () xtensor8
namespace aurostd {  // namespace aurostd
  template<class utype>                                         // operator ()
  inline utype& xtensor8<utype>::operator()(int i,int j,int k,int l,int m,int n,int o,int p) const {
#ifndef __OPTIMIZE
    if(i>uindex[1]) {cerr << "XTENSOR8 -> i=" << i << " > uindex[1]=" << uindex[1] << endl;exit(0);}
    if(i<lindex[1]) {cerr << "XTENSOR8 -> i=" << i << " < lindex[1]=" << lindex[1] << endl;exit(0);}
    if(j>uindex[2]) {cerr << "XTENSOR8 -> j=" << j << " > uindex[2]=" << uindex[2] << endl;exit(0);}
    if(j<lindex[2]) {cerr << "XTENSOR8 -> j=" << j << " < lindex[2]=" << lindex[2] << endl;exit(0);}
    if(k>uindex[3]) {cerr << "XTENSOR8 -> k=" << k << " > uindex[3]=" << uindex[3] << endl;exit(0);}
    if(k<lindex[3]) {cerr << "XTENSOR8 -> k=" << k << " < lindex[3]=" << lindex[3] << endl;exit(0);}
    if(l>uindex[4]) {cerr << "XTENSOR8 -> l=" << l << " > uindex[4]=" << uindex[4] << endl;exit(0);}
    if(l<lindex[4]) {cerr << "XTENSOR8 -> l=" << l << " < lindex[4]=" << lindex[4] << endl;exit(0);}
    if(m>uindex[5]) {cerr << "XTENSOR8 -> m=" << m << " > uindex[5]=" << uindex[5] << endl;exit(0);}
    if(m<lindex[5]) {cerr << "XTENSOR8 -> m=" << m << " < lindex[5]=" << lindex[5] << endl;exit(0);}
    if(n>uindex[6]) {cerr << "XTENSOR8 -> n=" << n << " > uindex[6]=" << uindex[6] << endl;exit(0);}
    if(n<lindex[6]) {cerr << "XTENSOR8 -> n=" << n << " < lindex[6]=" << lindex[6] << endl;exit(0);}
    if(o>uindex[7]) {cerr << "XTENSOR8 -> o=" << o << " > uindex[7]=" << uindex[7] << endl;exit(0);}
    if(o<lindex[7]) {cerr << "XTENSOR8 -> o=" << o << " < lindex[7]=" << lindex[7] << endl;exit(0);}
    if(p>uindex[8]) {cerr << "XTENSOR8 -> p=" << p << " > uindex[8]=" << uindex[8] << endl;exit(0);}
    if(p<lindex[8]) {cerr << "XTENSOR8 -> p=" << p << " < lindex[8]=" << lindex[8] << endl;exit(0);}
#endif
    return corpus[i][j][k][l][m][n][o][p];
  }
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ------------------------------------------------------ math unary operators

// ***************************************************************************
// --------------------------------------------- operator xtensor3 += xtensor3
namespace aurostd {  // namespace aurostd
  template<class utype> inline xtensor3<utype>&
  xtensor3<utype>::operator +=(const xtensor3<utype>& t) {
    aurostd::xtensor3debug(*this,"operator += :");
    if(this->index[1]!=t.index[1] || this->index[2]!=t.index[2] || this->index[3]!=t.index[3]) {
      cerr << "failure in xtensor3 operator+= " << endl;
      exit(0);
    }
    for(int i=0;i<index[1];i++)
      for(int j=0;j<index[2];j++)
	for(int k=0;k<index[3];k++)
	  corpus[i+lindex[1]][j+lindex[2]][k+lindex[3]]+=
	    t[i+t.lindex[1]][j+t.lindex[2]][k+t.lindex[3]];
    return *this;
  }
}

// --------------------------------------------- operator xtensor4 += xtensor4
namespace aurostd {  // namespace aurostd
  template<class utype> inline xtensor4<utype>&
  xtensor4<utype>::operator +=(const xtensor4<utype>& t) {
    aurostd::xtensor4debug(*this,"operator += :");
    if(this->index[1]!=t.index[1] || this->index[2]!=t.index[2] || this->index[3]!=t.index[3] ||
       this->index[4]!=t.index[4]) {
      cerr << "failure in xtensor4 operator+= " << endl;
      exit(0);
    }
    for(int i=0;i<index[1];i++)
      for(int j=0;j<index[2];j++)
	for(int k=0;k<index[3];k++)
	  for(int l=0;l<index[4];l++)
	    corpus[i+lindex[1]][j+lindex[2]][k+lindex[3]][l+lindex[4]]+=
	      t[i+t.lindex[1]][j+t.lindex[2]][k+t.lindex[3]][l+t.lindex[4]];
    return *this;
  }
}

// --------------------------------------------- operator xtensor5 += xtensor5
namespace aurostd {  // namespace aurostd
  template<class utype> inline xtensor5<utype>&
  xtensor5<utype>::operator +=(const xtensor5<utype>& t) {
    aurostd::xtensor5debug(*this,"operator += :");
    if(this->index[1]!=t.index[1] || this->index[2]!=t.index[2] || this->index[3]!=t.index[3] ||
       this->index[4]!=t.index[4] || this->index[5]!=t.index[5]) {
      cerr << "failure in xtensor5 operator+= " << endl;
      exit(0);
    }
    for(int i=0;i<index[1];i++)
      for(int j=0;j<index[2];j++)
	for(int k=0;k<index[3];k++)
	  for(int l=0;l<index[4];l++)
	    for(int m=0;m<index[5];m++)
	      corpus[i+lindex[1]][j+lindex[2]][k+lindex[3]][l+lindex[4]][m+lindex[5]]+=
		t[i+t.lindex[1]][j+t.lindex[2]][k+t.lindex[3]][l+t.lindex[4]][m+t.lindex[5]];
    return *this;
  }
}

// --------------------------------------------- operator xtensor6 += xtensor6
namespace aurostd {  // namespace aurostd
  template<class utype> inline xtensor6<utype>&
  xtensor6<utype>::operator +=(const xtensor6<utype>& t) {
    aurostd::xtensor6debug(*this,"operator += :");
    if(this->index[1]!=t.index[1] || this->index[2]!=t.index[2] || this->index[3]!=t.index[3] ||
       this->index[4]!=t.index[4] || this->index[5]!=t.index[5] || this->index[6]!=t.index[6]) {
      cerr << "failure in xtensor6 operator+= " << endl;
      exit(0);
    }
    for(int i=0;i<index[1];i++)
      for(int j=0;j<index[2];j++)
	for(int k=0;k<index[3];k++)
	  for(int l=0;l<index[4];l++)
	    for(int m=0;m<index[5];m++)
	      for(int n=0;n<index[6];n++)
		corpus[i+lindex[1]][j+lindex[2]][k+lindex[3]][l+lindex[4]][m+lindex[5]][n+lindex[6]]+=
		  t[i+t.lindex[1]][j+t.lindex[2]][k+t.lindex[3]][l+t.lindex[4]][m+t.lindex[5]][n+t.lindex[6]];
    return *this;
  }
}

// --------------------------------------------- operator xtensor7 += xtensor7
namespace aurostd {  // namespace aurostd
  template<class utype> inline xtensor7<utype>&
  xtensor7<utype>::operator +=(const xtensor7<utype>& t) {
    aurostd::xtensor7debug(*this,"operator += :");
    if(this->index[1]!=t.index[1] || this->index[2]!=t.index[2] || this->index[3]!=t.index[3] ||
       this->index[4]!=t.index[4] || this->index[5]!=t.index[5] || this->index[6]!=t.index[6] || this->index[7]!=t.index[7]) {
      cerr << "failure in xtensor7 operator+= " << endl;
      exit(0);
    }
    for(int i=0;i<index[1];i++)
      for(int j=0;j<index[2];j++)
	for(int k=0;k<index[3];k++)
	  for(int l=0;l<index[4];l++)
	    for(int m=0;m<index[5];m++)
	      for(int n=0;n<index[6];n++)
		for(int o=0;o<index[7];o++)
		  corpus[i+lindex[1]][j+lindex[2]][k+lindex[3]][l+lindex[4]][m+lindex[5]][n+lindex[6]][o+lindex[7]]+=
		    t[i+t.lindex[1]][j+t.lindex[2]][k+t.lindex[3]][l+t.lindex[4]][m+t.lindex[5]][n+t.lindex[6]][o+t.lindex[7]];
    return *this;
  }
}

// --------------------------------------------- operator xtensor8 += xtensor8
namespace aurostd {  // namespace aurostd
  template<class utype> inline xtensor8<utype>&
  xtensor8<utype>::operator +=(const xtensor8<utype>& t) {
    aurostd::xtensor8debug(*this,"operator += :");
    if(this->index[1]!=t.index[1] || this->index[2]!=t.index[2] || this->index[3]!=t.index[3] || this->index[4]!=t.index[4] ||
       this->index[5]!=t.index[5] || this->index[6]!=t.index[6] || this->index[7]!=t.index[7] || this->index[8]!=t.index[8]) {
      cerr << "failure in xtensor8 operator+= " << endl;
      exit(0);
    }
    for(int i=0;i<index[1];i++)
      for(int j=0;j<index[2];j++)
	for(int k=0;k<index[3];k++)
	  for(int l=0;l<index[4];l++)
	    for(int m=0;m<index[5];m++)
	      for(int n=0;n<index[6];n++)
		for(int o=0;o<index[7];o++)
		  for(int p=0;p<index[8];p++)
		    corpus[i+lindex[1]][j+lindex[2]][k+lindex[3]][l+lindex[4]][m+lindex[5]][n+lindex[6]][o+lindex[7]][p+lindex[8]]+=
		      t[i+t.lindex[1]][j+t.lindex[2]][k+t.lindex[3]][l+t.lindex[4]][m+t.lindex[5]][n+t.lindex[6]][o+t.lindex[7]][p+t.lindex[8]];
    return *this;
  }
}

// ***************************************************************************
// --------------------------------------------- operator xtensor3 -= xtensor3
namespace aurostd {  // namespace aurostd
  template<class utype> inline xtensor3<utype>&
  xtensor3<utype>::operator -=(const xtensor3<utype>& t) {
    aurostd::xtensor3debug(*this,"operator -= :");
    if(this->index[1]!=t.index[1] || this->index[2]!=t.index[2] || this->index[3]!=t.index[3]) {
      cerr << "failure in xtensor3 operator-= " << endl;
      exit(0);
    }
    for(int i=0;i<index[1];i++)
      for(int j=0;j<index[2];j++)
	for(int k=0;k<index[3];k++)
	  corpus[i+lindex[1]][j+lindex[2]][k+lindex[3]]-=
	    t[i+t.lindex[1]][j+t.lindex[2]][k+t.lindex[3]];
    return *this;
  }
}

// --------------------------------------------- operator xtensor4 -= xtensor4
namespace aurostd {  // namespace aurostd
  template<class utype> inline xtensor4<utype>&
  xtensor4<utype>::operator -=(const xtensor4<utype>& t) {
    aurostd::xtensor4debug(*this,"operator -= :");
    if(this->index[1]!=t.index[1] || this->index[2]!=t.index[2] || this->index[3]!=t.index[3] ||
       this->index[4]!=t.index[4]) {
      cerr << "failure in xtensor4 operator-= " << endl;
      exit(0);
    }
    for(int i=0;i<index[1];i++)
      for(int j=0;j<index[2];j++)
	for(int k=0;k<index[3];k++)
	  for(int l=0;l<index[4];l++)
	    corpus[i+lindex[1]][j+lindex[2]][k+lindex[3]][l+lindex[4]]-=
	      t[i+t.lindex[1]][j+t.lindex[2]][k+t.lindex[3]][l+t.lindex[4]];
    return *this;
  }
}

// --------------------------------------------- operator xtensor5 -= xtensor5
namespace aurostd {  // namespace aurostd
  template<class utype> inline xtensor5<utype>&
  xtensor5<utype>::operator -=(const xtensor5<utype>& t) {
    aurostd::xtensor5debug(*this,"operator -= :");
    if(this->index[1]!=t.index[1] || this->index[2]!=t.index[2] || this->index[3]!=t.index[3] ||
       this->index[4]!=t.index[4] || this->index[5]!=t.index[5]) {
      cerr << "failure in xtensor5 operator-= " << endl;
      exit(0);
    }
    for(int i=0;i<index[1];i++)
      for(int j=0;j<index[2];j++)
	for(int k=0;k<index[3];k++)
	  for(int l=0;l<index[4];l++)
	    for(int m=0;m<index[5];m++)
	      corpus[i+lindex[1]][j+lindex[2]][k+lindex[3]][l+lindex[4]][m+lindex[5]]-=
		t[i+t.lindex[1]][j+t.lindex[2]][k+t.lindex[3]][l+t.lindex[4]][m+t.lindex[5]];
    return *this;
  }
}

// --------------------------------------------- operator xtensor6 -= xtensor6
namespace aurostd {  // namespace aurostd
  template<class utype> inline xtensor6<utype>&
  xtensor6<utype>::operator -=(const xtensor6<utype>& t) {
    aurostd::xtensor6debug(*this,"operator -= :");
    if(this->index[1]!=t.index[1] || this->index[2]!=t.index[2] || this->index[3]!=t.index[3] ||
       this->index[4]!=t.index[4] || this->index[5]!=t.index[5] || this->index[6]!=t.index[6]) {
      cerr << "failure in xtensor6 operator-= " << endl;
      exit(0);
    }
    for(int i=0;i<index[1];i++)
      for(int j=0;j<index[2];j++)
	for(int k=0;k<index[3];k++)
	  for(int l=0;l<index[4];l++)
	    for(int m=0;m<index[5];m++)
	      for(int n=0;n<index[6];n++)
		corpus[i+lindex[1]][j+lindex[2]][k+lindex[3]][l+lindex[4]][m+lindex[5]][n+lindex[6]]-=
		  t[i+t.lindex[1]][j+t.lindex[2]][k+t.lindex[3]][l+t.lindex[4]][m+t.lindex[5]][n+t.lindex[6]];
    return *this;
  }
}

// --------------------------------------------- operator xtensor7 -= xtensor7
namespace aurostd {  // namespace aurostd
  template<class utype> inline xtensor7<utype>&
  xtensor7<utype>::operator -=(const xtensor7<utype>& t) {
    aurostd::xtensor7debug(*this,"operator -= :");
    if(this->index[1]!=t.index[1] || this->index[2]!=t.index[2] || this->index[3]!=t.index[3] ||
       this->index[4]!=t.index[4] || this->index[5]!=t.index[5] || this->index[6]!=t.index[6] || this->index[7]!=t.index[7]) {
      cerr << "failure in xtensor7 operator-= " << endl;
      exit(0);
    }
    for(int i=0;i<index[1];i++)
      for(int j=0;j<index[2];j++)
	for(int k=0;k<index[3];k++)
	  for(int l=0;l<index[4];l++)
	    for(int m=0;m<index[5];m++)
	      for(int n=0;n<index[6];n++)
		for(int o=0;o<index[7];o++)
		  corpus[i+lindex[1]][j+lindex[2]][k+lindex[3]][l+lindex[4]][m+lindex[5]][n+lindex[6]][o+lindex[7]]-=
		    t[i+t.lindex[1]][j+t.lindex[2]][k+t.lindex[3]][l+t.lindex[4]][m+t.lindex[5]][n+t.lindex[6]][o+t.lindex[7]];
    return *this;
  }
}

// --------------------------------------------- operator xtensor8 -= xtensor8
namespace aurostd {  // namespace aurostd
  template<class utype> inline xtensor8<utype>&
  xtensor8<utype>::operator -=(const xtensor8<utype>& t) {
    aurostd::xtensor8debug(*this,"operator -= :");
    if(this->index[1]!=t.index[1] || this->index[2]!=t.index[2] || this->index[3]!=t.index[3] || this->index[4]!=t.index[4] ||
       this->index[5]!=t.index[5] || this->index[6]!=t.index[6] || this->index[7]!=t.index[7] || this->index[8]!=t.index[8]) {
      cerr << "failure in xtensor8 operator-= " << endl;
      exit(0);
    }
    for(int i=0;i<index[1];i++)
      for(int j=0;j<index[2];j++)
	for(int k=0;k<index[3];k++)
	  for(int l=0;l<index[4];l++)
	    for(int m=0;m<index[5];m++)
	      for(int n=0;n<index[6];n++)
		for(int o=0;o<index[7];o++)
		  for(int p=0;p<index[8];p++)
		    corpus[i+lindex[1]][j+lindex[2]][k+lindex[3]][l+lindex[4]][m+lindex[5]][n+lindex[6]][o+lindex[7]][p+lindex[8]]-=
		      t[i+t.lindex[1]][j+t.lindex[2]][k+t.lindex[3]][l+t.lindex[4]][m+t.lindex[5]][n+t.lindex[6]][o+t.lindex[7]][p+t.lindex[8]];
    return *this;
  }
}

// ***************************************************************************
// -------------------------------------------------------- operator +xtensor3
namespace aurostd {  // namespace aurostd
  template<class utype>  
  xtensor3<utype> operator+(const xtensor3<utype>& a) {
    return a;
  }
}

// -------------------------------------------------------- operator +xtensor4
namespace aurostd {  // namespace aurostd
  template<class utype>  
  xtensor4<utype> operator+(const xtensor4<utype>& a) {
    return a;
  }
}

// -------------------------------------------------------- operator +xtensor5
namespace aurostd {  // namespace aurostd
  template<class utype>  
  xtensor5<utype> operator+(const xtensor5<utype>& a) {
    return a;
  }
}

// -------------------------------------------------------- operator +xtensor6
namespace aurostd {  // namespace aurostd
  template<class utype>  
  xtensor6<utype> operator+(const xtensor6<utype>& a) {
    return a;
  }
}

// -------------------------------------------------------- operator +xtensor7
namespace aurostd {  // namespace aurostd
  template<class utype>  
  xtensor7<utype> operator+(const xtensor7<utype>& a) {
    return a;
  }
}

// -------------------------------------------------------- operator +xtensor8
namespace aurostd {  // namespace aurostd
  template<class utype>  
  xtensor8<utype> operator+(const xtensor8<utype>& a) {
    return a;
  }
}

// ***************************************************************************
// -------------------------------------------------------- operator -xtensor3
namespace aurostd {  // namespace aurostd
  template<class utype>  
  xtensor3<utype> operator-(const xtensor3<utype>& a) {
    xtensor3<utype>c(a.uindex[1],a.uindex[2],a.uindex[3],
		     a.lindex[1],a.lindex[2],a.lindex[3]);
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
	for(int k=a.lindex[3];k<=a.uindex[3];k++)
	  c[i][j][k]=-a[i][j][k];
    return c;
  }
}

// -------------------------------------------------------- operator -xtensor4
namespace aurostd {  // namespace aurostd
  template<class utype>  
  xtensor4<utype> operator-(const xtensor4<utype>& a) {
    xtensor4<utype>c(a.uindex[1],a.uindex[2],a.uindex[3],a.uindex[4],
		     a.lindex[1],a.lindex[2],a.lindex[3],a.lindex[4]);
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
	for(int k=a.lindex[3];k<=a.uindex[3];k++)
	  for(int l=a.lindex[4];l<=a.uindex[4];l++)
	    c[i][j][k][l]=-a[i][j][k][l];
    return c;
  }
}

// -------------------------------------------------------- operator -xtensor5
namespace aurostd {  // namespace aurostd
  template<class utype>  
  xtensor5<utype> operator-(const xtensor5<utype>& a) {
    xtensor5<utype>c(a.uindex[1],a.uindex[2],a.uindex[3],a.uindex[4],a.uindex[5],
		     a.lindex[1],a.lindex[2],a.lindex[3],a.lindex[4],a.lindex[5]);
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
	for(int k=a.lindex[3];k<=a.uindex[3];k++)
	  for(int l=a.lindex[4];l<=a.uindex[4];l++)
	    for(int m=a.lindex[5];m<=a.uindex[5];m++)
	      c[i][j][k][l][m]=-a[i][j][k][l][m];
    return c;
  }
}

// -------------------------------------------------------- operator -xtensor6
namespace aurostd {  // namespace aurostd
  template<class utype>  
  xtensor6<utype> operator-(const xtensor6<utype>& a) {
    xtensor6<utype>c(a.uindex[1],a.uindex[2],a.uindex[3],a.uindex[4],a.uindex[5],a.uindex[6],
		     a.lindex[1],a.lindex[2],a.lindex[3],a.lindex[4],a.lindex[5],a.lindex[6]);
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
	for(int k=a.lindex[3];k<=a.uindex[3];k++)
	  for(int l=a.lindex[4];l<=a.uindex[4];l++)
	    for(int m=a.lindex[5];m<=a.uindex[5];m++)
	      for(int n=a.lindex[6];n<=a.uindex[6];n++)
		c[i][j][k][l][m][n]=-a[i][j][k][l][m][n];
    return c;
  }
}

// -------------------------------------------------------- operator -xtensor7
namespace aurostd {  // namespace aurostd
  template<class utype>  
  xtensor7<utype> operator-(const xtensor7<utype>& a) {
    xtensor7<utype>c(a.uindex[1],a.uindex[2],a.uindex[3],a.uindex[4],a.uindex[5],a.uindex[6],a.uindex[7],
		     a.lindex[1],a.lindex[2],a.lindex[3],a.lindex[4],a.lindex[5],a.lindex[6],a.lindex[7]);
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
	for(int k=a.lindex[3];k<=a.uindex[3];k++)
	  for(int l=a.lindex[4];l<=a.uindex[4];l++)
	    for(int m=a.lindex[5];m<=a.uindex[5];m++)
	      for(int n=a.lindex[6];n<=a.uindex[6];n++)
		for(int o=a.lindex[7];o<=a.uindex[7];o++)
		  c[i][j][k][l][m][n][o]=-a[i][j][k][l][m][n][o];
    return c;
  }
}

// -------------------------------------------------------- operator -xtensor8
namespace aurostd {  // namespace aurostd
  template<class utype>  
  xtensor8<utype> operator-(const xtensor8<utype>& a) {
    xtensor8<utype>c(a.uindex[1],a.uindex[2],a.uindex[3],a.uindex[4],a.uindex[5],a.uindex[6],a.uindex[7],a.uindex[8],
		     a.lindex[1],a.lindex[2],a.lindex[3],a.lindex[4],a.lindex[5],a.lindex[6],a.lindex[7],a.lindex[8]);
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
	for(int k=a.lindex[3];k<=a.uindex[3];k++)
	  for(int l=a.lindex[4];l<=a.uindex[4];l++)
	    for(int m=a.lindex[5];m<=a.uindex[5];m++)
	      for(int n=a.lindex[6];n<=a.uindex[6];n++)
		for(int o=a.lindex[7];o<=a.uindex[7];o++)
		  for(int p=a.lindex[8];p<=a.uindex[8];p++)
		    c[i][j][k][l][m][n][o][p]=-a[i][j][k][l][m][n][o][p];
    return c;
  }
}

// ****************************************************************************
// ****************************************************************************
// ------------------------------------------------------ math binary operators

// ****************************************************************************
// ----------------------------------------------- operator xtensor3 + xtensor3
namespace aurostd {
  template<class utype>
  xtensor3<utype> operator+(const xtensor3<utype>& a,const xtensor3<utype>& b) {
    aurostd::xtensor3debug(a,"operator +");
    aurostd::xtensor3debug(b,"operator +");
    if(a.index[1]!=b.index[1] || a.index[2]!=b.index[2] || a.index[3]!=b.index[3]) {
      cerr << "failure in operator xtensor3+xtensor3 = " << endl;
      exit(0);
    }
    xtensor3<utype>c(a.uindex[1],a.uindex[2],a.uindex[3],
		     a.lindex[1],a.lindex[2],a.lindex[3]);
    for(int i=0;i<a.index[1];i++)
      for(int j=0;j<a.index[2];j++)
	for(int k=0;k<a.index[3];k++)
	  c[i+c.lindex[1]][j+c.lindex[2]][k+c.lindex[3]]=
	    a[i+a.lindex[1]][j+a.lindex[2]][k+a.lindex[3]]+
	    b[i+b.lindex[1]][j+b.lindex[2]][k+b.lindex[3]];
    return c;
  }
}

// ----------------------------------------------- operator xtensor4 + xtensor4
namespace aurostd {
  template<class utype>
  xtensor4<utype> operator+(const xtensor4<utype>& a,const xtensor4<utype>& b) {
    aurostd::xtensor4debug(a,"operator +");
    aurostd::xtensor4debug(b,"operator +");
    if(a.index[1]!=b.index[1] || a.index[2]!=b.index[2] || a.index[3]!=b.index[3] ||
       a.index[4]!=b.index[4] ) {
      cerr << "failure in operator xtensor4+xtensor4 = " << endl;
      exit(0);
    }
    xtensor4<utype>c(a.uindex[1],a.uindex[2],a.uindex[3],a.uindex[4],
		     a.lindex[1],a.lindex[2],a.lindex[3],a.lindex[4]);
    for(int i=0;i<a.index[1];i++)
      for(int j=0;j<a.index[2];j++)
	for(int k=0;k<a.index[3];k++)
	  for(int l=0;l<a.index[4];l++)
	    c[i+c.lindex[1]][j+c.lindex[2]][k+c.lindex[3]][l+c.lindex[4]]=
	      a[i+a.lindex[1]][j+a.lindex[2]][k+a.lindex[3]][l+a.lindex[4]]+
	      b[i+b.lindex[1]][j+b.lindex[2]][k+b.lindex[3]][l+b.lindex[4]];
    return c;
  }
}

// ----------------------------------------------- operator xtensor5 + xtensor5
namespace aurostd {
  template<class utype>
  xtensor5<utype> operator+(const xtensor5<utype>& a,const xtensor5<utype>& b) {
    aurostd::xtensor5debug(a,"operator +");
    aurostd::xtensor5debug(b,"operator +");
    if(a.index[1]!=b.index[1] || a.index[2]!=b.index[2] || a.index[3]!=b.index[3] ||
       a.index[4]!=b.index[4] || a.index[5]!=b.index[5] ) {
      cerr << "failure in operator xtensor5+xtensor5 = " << endl;
      exit(0);
    }
    xtensor5<utype>c(a.uindex[1],a.uindex[2],a.uindex[3],a.uindex[4],a.uindex[5],
		     a.lindex[1],a.lindex[2],a.lindex[3],a.lindex[4],a.lindex[5]);
    for(int i=0;i<a.index[1];i++)
      for(int j=0;j<a.index[2];j++)
	for(int k=0;k<a.index[3];k++)
	  for(int l=0;l<a.index[4];l++)
	    for(int m=0;m<a.index[5];m++)
	      c[i+c.lindex[1]][j+c.lindex[2]][k+c.lindex[3]][l+c.lindex[4]][m+c.lindex[5]]=
		a[i+a.lindex[1]][j+a.lindex[2]][k+a.lindex[3]][l+a.lindex[4]][m+a.lindex[5]]+
		b[i+b.lindex[1]][j+b.lindex[2]][k+b.lindex[3]][l+b.lindex[4]][m+b.lindex[5]];
    return c;
  }
}

// ----------------------------------------------- operator xtensor6 + xtensor6
namespace aurostd {
  template<class utype>
  xtensor6<utype> operator+(const xtensor6<utype>& a,const xtensor6<utype>& b) {
    aurostd::xtensor6debug(a,"operator +");
    aurostd::xtensor6debug(b,"operator +");
    if(a.index[1]!=b.index[1] || a.index[2]!=b.index[2] || a.index[3]!=b.index[3] ||
       a.index[4]!=b.index[4] || a.index[5]!=b.index[5] || a.index[6]!=b.index[6]) {
      cerr << "failure in operator xtensor6+xtensor6 = " << endl;
      exit(0);
    }
    xtensor6<utype>c(a.uindex[1],a.uindex[2],a.uindex[3],a.uindex[4],a.uindex[5],a.uindex[6],
		     a.lindex[1],a.lindex[2],a.lindex[3],a.lindex[4],a.lindex[5],a.lindex[6]);
    for(int i=0;i<a.index[1];i++)
      for(int j=0;j<a.index[2];j++)
	for(int k=0;k<a.index[3];k++)
	  for(int l=0;l<a.index[4];l++)
	    for(int m=0;m<a.index[5];m++)
	      for(int n=0;n<a.index[6];n++)
		c[i+c.lindex[1]][j+c.lindex[2]][k+c.lindex[3]][l+c.lindex[4]][m+c.lindex[5]][n+c.lindex[6]]=
		  a[i+a.lindex[1]][j+a.lindex[2]][k+a.lindex[3]][l+a.lindex[4]][m+a.lindex[5]][n+a.lindex[6]]+
		  b[i+b.lindex[1]][j+b.lindex[2]][k+b.lindex[3]][l+b.lindex[4]][m+b.lindex[5]][n+b.lindex[6]];
    return c;
  }
}

// ----------------------------------------------- operator xtensor7 + xtensor7
namespace aurostd {
  template<class utype>
  xtensor7<utype> operator+(const xtensor7<utype>& a,const xtensor7<utype>& b) {
    aurostd::xtensor7debug(a,"operator +");
    aurostd::xtensor7debug(b,"operator +");
    if(a.index[1]!=b.index[1] || a.index[2]!=b.index[2] || a.index[3]!=b.index[3] ||
       a.index[4]!=b.index[4] || a.index[5]!=b.index[5] || a.index[6]!=b.index[6] || a.index[7]!=b.index[7]) {
      cerr << "failure in operator xtensor7+xtensor7 = " << endl;
      exit(0);
    }
    xtensor7<utype>c(a.uindex[1],a.uindex[2],a.uindex[3],a.uindex[4],a.uindex[5],a.uindex[6],a.uindex[7],
		     a.lindex[1],a.lindex[2],a.lindex[3],a.lindex[4],a.lindex[5],a.lindex[6],a.lindex[7]);
    for(int i=0;i<a.index[1];i++)
      for(int j=0;j<a.index[2];j++)
	for(int k=0;k<a.index[3];k++)
	  for(int l=0;l<a.index[4];l++)
	    for(int m=0;m<a.index[5];m++)
	      for(int n=0;n<a.index[6];n++)
		for(int o=0;o<a.index[7];o++)
		  c[i+c.lindex[1]][j+c.lindex[2]][k+c.lindex[3]][l+c.lindex[4]][m+c.lindex[5]][n+c.lindex[6]][o+c.lindex[7]]=
		    a[i+a.lindex[1]][j+a.lindex[2]][k+a.lindex[3]][l+a.lindex[4]][m+a.lindex[5]][n+a.lindex[6]][o+a.lindex[7]]+
		    b[i+b.lindex[1]][j+b.lindex[2]][k+b.lindex[3]][l+b.lindex[4]][m+b.lindex[5]][n+b.lindex[6]][o+b.lindex[7]];
    return c;
  }
}

// ----------------------------------------------- operator xtensor8 + xtensor8
namespace aurostd {
  template<class utype>
  xtensor8<utype> operator+(const xtensor8<utype>& a,const xtensor8<utype>& b) {
    aurostd::xtensor8debug(a,"operator +");
    aurostd::xtensor8debug(b,"operator +");
    if(a.index[1]!=b.index[1] || a.index[2]!=b.index[2] || a.index[3]!=b.index[3] || a.index[4]!=b.index[4] ||
       a.index[5]!=b.index[5] || a.index[6]!=b.index[6] || a.index[7]!=b.index[7] || a.index[8]!=b.index[8]) {
      cerr << "failure in operator xtensor8+xtensor8 = " << endl;
      exit(0);
    }
    xtensor8<utype>c(a.uindex[1],a.uindex[2],a.uindex[3],a.uindex[4],a.uindex[5],a.uindex[6],a.uindex[7],a.uindex[8],
		     a.lindex[1],a.lindex[2],a.lindex[3],a.lindex[4],a.lindex[5],a.lindex[6],a.lindex[7],a.lindex[8]);
    for(int i=0;i<a.index[1];i++)
      for(int j=0;j<a.index[2];j++)
	for(int k=0;k<a.index[3];k++)
	  for(int l=0;l<a.index[4];l++)
	    for(int m=0;m<a.index[5];m++)
	      for(int n=0;n<a.index[6];n++)
		for(int o=0;o<a.index[7];o++)
		  for(int p=0;p<a.index[8];p++)
		    c[i+c.lindex[1]][j+c.lindex[2]][k+c.lindex[3]][l+c.lindex[4]][m+c.lindex[5]][n+c.lindex[6]][o+c.lindex[7]][p+c.lindex[8]]=
		      a[i+a.lindex[1]][j+a.lindex[2]][k+a.lindex[3]][l+a.lindex[4]][m+a.lindex[5]][n+a.lindex[6]][o+a.lindex[7]][p+a.lindex[8]]+
		      b[i+b.lindex[1]][j+b.lindex[2]][k+b.lindex[3]][l+b.lindex[4]][m+b.lindex[5]][n+b.lindex[6]][o+b.lindex[7]][p+b.lindex[8]];
    return c;
  }
}

// ****************************************************************************
// ----------------------------------------------- operator xtensor3 - xtensor3
namespace aurostd {
  template<class utype>
  xtensor3<utype> operator-(const xtensor3<utype>& a,const xtensor3<utype>& b) {
    aurostd::xtensor3debug(a,"operator -");
    aurostd::xtensor3debug(b,"operator -");
    if(a.index[1]!=b.index[1] || a.index[2]!=b.index[2] || a.index[3]!=b.index[3]) {
      cerr << "failure in operator xtensor3-xtensor3 = " << endl;
      exit(0);
    }
    xtensor3<utype>c(a.uindex[1],a.uindex[2],a.uindex[3],
		     a.lindex[1],a.lindex[2],a.lindex[3]);
    for(int i=0;i<a.index[1];i++)
      for(int j=0;j<a.index[2];j++)
	for(int k=0;k<a.index[3];k++)
	  c[i+c.lindex[1]][j+c.lindex[2]][k+c.lindex[3]]=
	    a[i+a.lindex[1]][j+a.lindex[2]][k+a.lindex[3]]-
	    b[i+b.lindex[1]][j+b.lindex[2]][k+b.lindex[3]];
    return c;
  }
}

// ----------------------------------------------- operator xtensor4 - xtensor4
namespace aurostd {
  template<class utype>
  xtensor4<utype> operator-(const xtensor4<utype>& a,const xtensor4<utype>& b) {
    aurostd::xtensor4debug(a,"operator -");
    aurostd::xtensor4debug(b,"operator -");
    if(a.index[1]!=b.index[1] || a.index[2]!=b.index[2] || a.index[3]!=b.index[3] ||
       a.index[4]!=b.index[4] ) {
      cerr << "failure in operator xtensor4-xtensor4 = " << endl;
      exit(0);
    }
    xtensor4<utype>c(a.uindex[1],a.uindex[2],a.uindex[3],a.uindex[4],
		     a.lindex[1],a.lindex[2],a.lindex[3],a.lindex[4]);
    for(int i=0;i<a.index[1];i++)
      for(int j=0;j<a.index[2];j++)
	for(int k=0;k<a.index[3];k++)
	  for(int l=0;l<a.index[4];l++)
	    c[i+c.lindex[1]][j+c.lindex[2]][k+c.lindex[3]][l+c.lindex[4]]=
	      a[i+a.lindex[1]][j+a.lindex[2]][k+a.lindex[3]][l+a.lindex[4]]-
	      b[i+b.lindex[1]][j+b.lindex[2]][k+b.lindex[3]][l+b.lindex[4]];
    return c;
  }
}

// ----------------------------------------------- operator xtensor5 - xtensor5
namespace aurostd {
  template<class utype>
  xtensor5<utype> operator-(const xtensor5<utype>& a,const xtensor5<utype>& b) {
    aurostd::xtensor5debug(a,"operator -");
    aurostd::xtensor5debug(b,"operator -");
    if(a.index[1]!=b.index[1] || a.index[2]!=b.index[2] || a.index[3]!=b.index[3] ||
       a.index[4]!=b.index[4] || a.index[5]!=b.index[5] ) {
      cerr << "failure in operator xtensor5-xtensor5 = " << endl;
      exit(0);
    }
    xtensor5<utype>c(a.uindex[1],a.uindex[2],a.uindex[3],a.uindex[4],a.uindex[5],
		     a.lindex[1],a.lindex[2],a.lindex[3],a.lindex[4],a.lindex[5]);
    for(int i=0;i<a.index[1];i++)
      for(int j=0;j<a.index[2];j++)
	for(int k=0;k<a.index[3];k++)
	  for(int l=0;l<a.index[4];l++)
	    for(int m=0;m<a.index[5];m++)
	      c[i+c.lindex[1]][j+c.lindex[2]][k+c.lindex[3]][l+c.lindex[4]][m+c.lindex[5]]=
		a[i+a.lindex[1]][j+a.lindex[2]][k+a.lindex[3]][l+a.lindex[4]][m+a.lindex[5]]-
		b[i+b.lindex[1]][j+b.lindex[2]][k+b.lindex[3]][l+b.lindex[4]][m+b.lindex[5]];
    return c;
  }
}

// ----------------------------------------------- operator xtensor6 - xtensor6
namespace aurostd {
  template<class utype>
  xtensor6<utype> operator-(const xtensor6<utype>& a,const xtensor6<utype>& b) {
    aurostd::xtensor6debug(a,"operator -");
    aurostd::xtensor6debug(b,"operator -");
    if(a.index[1]!=b.index[1] || a.index[2]!=b.index[2] || a.index[3]!=b.index[3] ||
       a.index[4]!=b.index[4] || a.index[5]!=b.index[5] || a.index[6]!=b.index[6]) {
      cerr << "failure in operator xtensor6-xtensor6 = " << endl;
      exit(0);
    }
    xtensor6<utype>c(a.uindex[1],a.uindex[2],a.uindex[3],a.uindex[4],a.uindex[5],a.uindex[6],
		     a.lindex[1],a.lindex[2],a.lindex[3],a.lindex[4],a.lindex[5],a.lindex[6]);
    for(int i=0;i<a.index[1];i++)
      for(int j=0;j<a.index[2];j++)
	for(int k=0;k<a.index[3];k++)
	  for(int l=0;l<a.index[4];l++)
	    for(int m=0;m<a.index[5];m++)
	      for(int n=0;n<a.index[6];n++)
		c[i+c.lindex[1]][j+c.lindex[2]][k+c.lindex[3]][l+c.lindex[4]][m+c.lindex[5]][n+c.lindex[6]]=
		  a[i+a.lindex[1]][j+a.lindex[2]][k+a.lindex[3]][l+a.lindex[4]][m+a.lindex[5]][n+a.lindex[6]]-
		  b[i+b.lindex[1]][j+b.lindex[2]][k+b.lindex[3]][l+b.lindex[4]][m+b.lindex[5]][n+b.lindex[6]];
    return c;
  }
}

// ----------------------------------------------- operator xtensor7 - xtensor7
namespace aurostd {
  template<class utype>
  xtensor7<utype> operator-(const xtensor7<utype>& a,const xtensor7<utype>& b) {
    aurostd::xtensor7debug(a,"operator -");
    aurostd::xtensor7debug(b,"operator -");
    if(a.index[1]!=b.index[1] || a.index[2]!=b.index[2] || a.index[3]!=b.index[3] ||
       a.index[4]!=b.index[4] || a.index[5]!=b.index[5] || a.index[6]!=b.index[6] || a.index[7]!=b.index[7]) {
      cerr << "failure in operator xtensor7-xtensor7 = " << endl;
      exit(0);
    }
    xtensor7<utype>c(a.uindex[1],a.uindex[2],a.uindex[3],a.uindex[4],a.uindex[5],a.uindex[6],a.uindex[7],
		     a.lindex[1],a.lindex[2],a.lindex[3],a.lindex[4],a.lindex[5],a.lindex[6],a.lindex[7]);
    for(int i=0;i<a.index[1];i++)
      for(int j=0;j<a.index[2];j++)
	for(int k=0;k<a.index[3];k++)
	  for(int l=0;l<a.index[4];l++)
	    for(int m=0;m<a.index[5];m++)
	      for(int n=0;n<a.index[6];n++)
		for(int o=0;o<a.index[7];o++)
		  c[i+c.lindex[1]][j+c.lindex[2]][k+c.lindex[3]][l+c.lindex[4]][m+c.lindex[5]][n+c.lindex[6]][o+c.lindex[7]]=
		    a[i+a.lindex[1]][j+a.lindex[2]][k+a.lindex[3]][l+a.lindex[4]][m+a.lindex[5]][n+a.lindex[6]][o+a.lindex[7]]-
		    b[i+b.lindex[1]][j+b.lindex[2]][k+b.lindex[3]][l+b.lindex[4]][m+b.lindex[5]][n+b.lindex[6]][o+b.lindex[7]];
    return c;
  }
}

// ----------------------------------------------- operator xtensor8 - xtensor8
namespace aurostd {
  template<class utype>
  xtensor8<utype> operator-(const xtensor8<utype>& a,const xtensor8<utype>& b) {
    aurostd::xtensor8debug(a,"operator -");
    aurostd::xtensor8debug(b,"operator -");
    if(a.index[1]!=b.index[1] || a.index[2]!=b.index[2] || a.index[3]!=b.index[3] || a.index[4]!=b.index[4] ||
       a.index[5]!=b.index[5] || a.index[6]!=b.index[6] || a.index[7]!=b.index[7] || a.index[8]!=b.index[8]) {
      cerr << "failure in operator xtensor8+xtensor8 = " << endl;
      exit(0);
    }
    xtensor8<utype>c(a.uindex[1],a.uindex[2],a.uindex[3],a.uindex[4],a.uindex[5],a.uindex[6],a.uindex[7],a.uindex[8],
		     a.lindex[1],a.lindex[2],a.lindex[3],a.lindex[4],a.lindex[5],a.lindex[6],a.lindex[7],a.lindex[8]);
    for(int i=0;i<a.index[1];i++)
      for(int j=0;j<a.index[2];j++)
	for(int k=0;k<a.index[3];k++)
	  for(int l=0;l<a.index[4];l++)
	    for(int m=0;m<a.index[5];m++)
	      for(int n=0;n<a.index[6];n++)
		for(int o=0;o<a.index[7];o++)
		  for(int p=0;p<a.index[8];p++)
		    c[i+c.lindex[1]][j+c.lindex[2]][k+c.lindex[3]][l+c.lindex[4]][m+c.lindex[5]][n+c.lindex[6]][o+c.lindex[7]][p+c.lindex[8]]=
		      a[i+a.lindex[1]][j+a.lindex[2]][k+a.lindex[3]][l+a.lindex[4]][m+a.lindex[5]][n+a.lindex[6]][o+a.lindex[7]][p+a.lindex[8]]-
		      b[i+b.lindex[1]][j+b.lindex[2]][k+b.lindex[3]][l+b.lindex[4]][m+b.lindex[5]][n+b.lindex[6]][o+b.lindex[7]][p+b.lindex[8]];
    return c;
  }
}

// ****************************************************************************
// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype> xtensor3<utype>             // operator xtensor3 * scalar
  operator*(const utype s,const xtensor3<utype>& a) {
    xtensor3<utype>c(a.uindex[1],a.uindex[2],a.uindex[3],
		     a.lindex[1],a.lindex[2],a.lindex[3]);
    for(int i=c.lindex[1];i<=c.uindex[1];i++)
      for(int j=c.lindex[2];j<=c.uindex[2];j++)
	for(int k=c.lindex[3];k<=c.uindex[3];k++)
	  c[i][j][k]=(utype) a[i][j][k]*(utype) s;
    return c;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype> xtensor4<utype>             // operator xtensor4 * scalar
  operator*(const utype s,const xtensor4<utype>& a) {
    xtensor4<utype>c(a.uindex[1],a.uindex[2],a.uindex[3],a.uindex[4],
		     a.lindex[1],a.lindex[2],a.lindex[3],a.lindex[4]);
    for(int i=c.lindex[1];i<=c.uindex[1];i++)
      for(int j=c.lindex[2];j<=c.uindex[2];j++)
	for(int k=c.lindex[3];k<=c.uindex[3];k++)
	  for(int l=c.lindex[4];l<=c.uindex[4];l++)
	    c[i][j][k][l]=(utype) a[i][j][k][l]*(utype) s;
    return c;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype> xtensor5<utype>             // operator xtensor5 * scalar
  operator*(const utype s,const xtensor5<utype>& a) {
    xtensor5<utype>c(a.uindex[1],a.uindex[2],a.uindex[3],a.uindex[4],a.uindex[5],
		     a.lindex[1],a.lindex[2],a.lindex[3],a.lindex[4],a.lindex[5]);
    for(int i=c.lindex[1];i<=c.uindex[1];i++)
      for(int j=c.lindex[2];j<=c.uindex[2];j++)
	for(int k=c.lindex[3];k<=c.uindex[3];k++)
	  for(int l=c.lindex[4];l<=c.uindex[4];l++)
	    for(int m=c.lindex[5];m<=c.uindex[5];m++)
	      c[i][j][k][l][m]=(utype) a[i][j][k][l][m]*(utype) s;
    return c;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype> xtensor6<utype>             // operator xtensor6 * scalar
  operator*(const utype s,const xtensor6<utype>& a) {
    xtensor6<utype>c(a.uindex[1],a.uindex[2],a.uindex[3],a.uindex[4],a.uindex[5],a.uindex[6],
		     a.lindex[1],a.lindex[2],a.lindex[3],a.lindex[4],a.lindex[5],a.lindex[6]);
    for(int i=c.lindex[1];i<=c.uindex[1];i++)
      for(int j=c.lindex[2];j<=c.uindex[2];j++)
	for(int k=c.lindex[3];k<=c.uindex[3];k++)
	  for(int l=c.lindex[4];l<=c.uindex[4];l++)
	    for(int m=c.lindex[5];m<=c.uindex[5];m++)
	      for(int n=a.lindex[6];n<=a.uindex[6];n++)
		c[i][j][k][l][m][n]=(utype) a[i][j][k][l][m][n]*(utype) s;
    return c;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype> xtensor7<utype>             // operator xtensor7 * scalar
  operator*(const utype s,const xtensor7<utype>& a) {
    xtensor7<utype>c(a.uindex[1],a.uindex[2],a.uindex[3],a.uindex[4],a.uindex[5],a.uindex[6],a.uindex[7],
		     a.lindex[1],a.lindex[2],a.lindex[3],a.lindex[4],a.lindex[5],a.lindex[6],a.lindex[7]);
    for(int i=c.lindex[1];i<=c.uindex[1];i++)
      for(int j=c.lindex[2];j<=c.uindex[2];j++)
	for(int k=c.lindex[3];k<=c.uindex[3];k++)
	  for(int l=c.lindex[4];l<=c.uindex[4];l++)
	    for(int m=c.lindex[5];m<=c.uindex[5];m++)
	      for(int n=a.lindex[6];n<=a.uindex[6];n++)
		for(int o=a.lindex[7];o<=a.uindex[7];o++)
		  c[i][j][k][l][m][n][o]=(utype) a[i][j][k][l][m][n][o]*(utype) s;
    return c;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype> xtensor8<utype>             // operator xtensor8 * scalar
  operator*(const utype s,const xtensor8<utype>& a) {
    xtensor8<utype>c(a.uindex[1],a.uindex[2],a.uindex[3],a.uindex[4],a.uindex[5],a.uindex[6],a.uindex[7],a.uindex[8],
		     a.lindex[1],a.lindex[2],a.lindex[3],a.lindex[4],a.lindex[5],a.lindex[6],a.lindex[7],a.lindex[8]);
    for(int i=c.lindex[1];i<=c.uindex[1];i++)
      for(int j=c.lindex[2];j<=c.uindex[2];j++)
	for(int k=c.lindex[3];k<=c.uindex[3];k++)
	  for(int l=c.lindex[4];l<=c.uindex[4];l++)
	    for(int m=c.lindex[5];m<=c.uindex[5];m++)
	      for(int n=a.lindex[6];n<=a.uindex[6];n++)
		for(int o=a.lindex[7];o<=a.uindex[7];o++)
		  for(int p=a.lindex[8];p<=a.uindex[8];p++)
		    c[i][j][k][l][m][n][o][p]=(utype) a[i][j][k][l][m][n][o][p]*(utype) s;
    return c;
  }
}

// ****************************************************************************
// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype> xtensor3<utype>          //  operator scalar * xtensor3
  operator*(const xtensor3<utype>& a,const utype s) {
    return s*a;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype> xtensor4<utype>          //  operator scalar * xtensor4
  operator*(const xtensor4<utype>& a,const utype s) {
    return s*a;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype> xtensor5<utype>          //  operator scalar * xtensor5
  operator*(const xtensor5<utype>& a,const utype s) {
    return s*a;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype> xtensor6<utype>          //  operator scalar * xtensor6
  operator*(const xtensor6<utype>& a,const utype s) {
    return s*a;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype> xtensor7<utype>          //  operator scalar * xtensor7
  operator*(const xtensor7<utype>& a,const utype s) {
    return s*a;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype> xtensor8<utype>          //  operator scalar * xtensor8
  operator*(const xtensor8<utype>& a,const utype s) {
    return s*a;
  }
}

// ****************************************************************************
// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype> xtensor3<utype>           // operator xtensor3 / scalar
  operator/(const xtensor3<utype>& a,const utype s) {
    return (utype) (1/s)*a;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype> xtensor4<utype>           // operator xtensor4 / scalar
  operator/(const xtensor4<utype>& a,const utype s) {
    return (utype) (1/s)*a;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype> xtensor5<utype>           // operator xtensor5 / scalar
  operator/(const xtensor5<utype>& a,const utype s) {
    return (utype) (1/s)*a;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype> xtensor6<utype>           // operator xtensor6 / scalar
  operator/(const xtensor6<utype>& a,const utype s) {
    return (utype) (1/s)*a;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype> xtensor7<utype>           // operator xtensor7 / scalar
  operator/(const xtensor7<utype>& a,const utype s) {
    return (utype) (1/s)*a;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype> xtensor8<utype>           // operator xtensor8 / scalar
  operator/(const xtensor8<utype>& a,const utype s) {
    return (utype) (1/s)*a;
  }
}

// ****************************************************************************
// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype>                            // function sign xtensor3<>
  xtensor3<utype> sign(const xtensor3<utype>& a) {
    aurostd::xtensor3debug(a,"function sign");
    xtensor3<utype> c(a);
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
	for(int k=a.lindex[3];k<=a.uindex[3];k++)
	  c[i][j][k]=(utype) aurostd::sign(a[i][j][k]);
    return (xtensor3<utype>) c;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype>                            // function sign xtensor4<>
  xtensor4<utype> sign(const xtensor4<utype>& a) {
    aurostd::xtensor4debug(a,"function sign");
    xtensor4<utype> c(a);
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
	for(int k=a.lindex[3];k<=a.uindex[3];k++)
	  for(int l=a.lindex[4];l<=a.uindex[4];l++)
	    c[i][j][k][l]=(utype) aurostd::sign(a[i][j][k][l]);
    return (xtensor4<utype>) c;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype>                            // function sign xtensor5<>
  xtensor5<utype> sign(const xtensor5<utype>& a) {
    aurostd::xtensor5debug(a,"function sign");
    xtensor5<utype> c(a);
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
	for(int k=a.lindex[3];k<=a.uindex[3];k++)
	  for(int l=a.lindex[4];l<=a.uindex[4];l++)
	    for(int m=a.lindex[5];m<=a.uindex[5];m++)
	      c[i][j][k][l][m]=(utype) aurostd::sign(a[i][j][k][l][m]);
    return (xtensor5<utype>) c;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype>                            // function sign xtensor6<>
  xtensor6<utype> sign(const xtensor6<utype>& a) {
    aurostd::xtensor6debug(a,"function sign");
    xtensor6<utype> c(a);
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
	for(int k=a.lindex[3];k<=a.uindex[3];k++)
	  for(int l=a.lindex[4];l<=a.uindex[4];l++)
	    for(int m=a.lindex[5];m<=a.uindex[5];m++)
	      for(int n=a.lindex[6];n<=a.uindex[6];n++)
		c[i][j][k][l][m][n]=(utype) aurostd::sign(a[i][j][k][l][m][n]);
    return (xtensor6<utype>) c;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype>                            // function sign xtensor7<>
  xtensor7<utype> sign(const xtensor7<utype>& a) {
    aurostd::xtensor7debug(a,"function sign");
    xtensor7<utype> c(a);
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
	for(int k=a.lindex[3];k<=a.uindex[3];k++)
	  for(int l=a.lindex[4];l<=a.uindex[4];l++)
	    for(int m=a.lindex[5];m<=a.uindex[5];m++)
	      for(int n=a.lindex[6];n<=a.uindex[6];n++)
		for(int o=a.lindex[7];o<=a.uindex[7];o++)
		  c[i][j][k][l][m][n][o]=(utype) aurostd::sign(a[i][j][k][l][m][n][o]);
    return (xtensor7<utype>) c;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype>                            // function sign xtensor8<>
  xtensor8<utype> sign(const xtensor8<utype>& a) {
    aurostd::xtensor8debug(a,"function sign");
    xtensor8<utype> c(a);
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
	for(int k=a.lindex[3];k<=a.uindex[3];k++)
	  for(int l=a.lindex[4];l<=a.uindex[4];l++)
	    for(int m=a.lindex[5];m<=a.uindex[5];m++)
	      for(int n=a.lindex[6];n<=a.uindex[6];n++)
		for(int o=a.lindex[7];o<=a.uindex[7];o++)
		  for(int p=a.lindex[8];p<=a.uindex[8];p++)
		    c[i][j][k][l][m][n][o][p]=(utype) aurostd::sign(a[i][j][k][l][m][n][o][p]);
    return (xtensor8<utype>) c;
  }
}

// ****************************************************************************
// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype>                            // function nint xtensor3<>
  xtensor3<utype> nint(const xtensor3<utype>& a) {
    aurostd::xtensor3debug(a,"function nint");
    xtensor3<utype> c(a);
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
        for(int k=a.lindex[3];k<=a.uindex[3];k++)
          c[i][j][k]=(utype) aurostd::nint(a[i][j][k]);
    return (xtensor3<utype>) c;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype>                            // function nint xtensor4<>
  xtensor4<utype> nint(const xtensor4<utype>& a) {
    aurostd::xtensor4debug(a,"function nint");
    xtensor4<utype> c(a);
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
        for(int k=a.lindex[3];k<=a.uindex[3];k++)
          for(int l=a.lindex[4];l<=a.uindex[4];l++)
            c[i][j][k][l]=(utype) aurostd::nint(a[i][j][k][l]);
    return (xtensor4<utype>) c;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype>                            // function nint xtensor5<>
  xtensor5<utype> nint(const xtensor5<utype>& a) {
    aurostd::xtensor5debug(a,"function nint");
    xtensor5<utype> c(a);
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
        for(int k=a.lindex[3];k<=a.uindex[3];k++)
          for(int l=a.lindex[4];l<=a.uindex[4];l++)
            for(int m=a.lindex[5];m<=a.uindex[5];m++)
              c[i][j][k][l][m]=(utype) aurostd::nint(a[i][j][k][l][m]);
    return (xtensor5<utype>) c;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype>                            // function nint xtensor6<>
  xtensor6<utype> nint(const xtensor6<utype>& a) {
    aurostd::xtensor6debug(a,"function nint");
    xtensor6<utype> c(a);
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
        for(int k=a.lindex[3];k<=a.uindex[3];k++)
          for(int l=a.lindex[4];l<=a.uindex[4];l++)
            for(int m=a.lindex[5];m<=a.uindex[5];m++)
              for(int n=a.lindex[6];n<=a.uindex[6];n++)
                c[i][j][k][l][m][n]=(utype) aurostd::nint(a[i][j][k][l][m][n]);
    return (xtensor6<utype>) c;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype>                            // function nint xtensor7<>
  xtensor7<utype> nint(const xtensor7<utype>& a) {
    aurostd::xtensor7debug(a,"function nint");
    xtensor7<utype> c(a);
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
        for(int k=a.lindex[3];k<=a.uindex[3];k++)
          for(int l=a.lindex[4];l<=a.uindex[4];l++)
            for(int m=a.lindex[5];m<=a.uindex[5];m++)
              for(int n=a.lindex[6];n<=a.uindex[6];n++)
                for(int o=a.lindex[7];o<=a.uindex[7];o++)
                  c[i][j][k][l][m][n][o]=(utype) aurostd::nint(a[i][j][k][l][m][n][o]);
    return (xtensor7<utype>) c;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype>                            // function nint xtensor8<>
  xtensor8<utype> nint(const xtensor8<utype>& a) {
    aurostd::xtensor8debug(a,"function nint");
    xtensor8<utype> c(a);
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
        for(int k=a.lindex[3];k<=a.uindex[3];k++)
          for(int l=a.lindex[4];l<=a.uindex[4];l++)
            for(int m=a.lindex[5];m<=a.uindex[5];m++)
              for(int n=a.lindex[6];n<=a.uindex[6];n++)
                for(int o=a.lindex[7];o<=a.uindex[7];o++)
                  for(int p=a.lindex[8];p<=a.uindex[8];p++)
                    c[i][j][k][l][m][n][o][p]=(utype) aurostd::nint(a[i][j][k][l][m][n][o][p]);
    return (xtensor8<utype>) c;
  }
}

// ****************************************************************************
// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype>                            // function abs xtensor3<>
  xtensor3<utype> abs(const xtensor3<utype>& a) {
    aurostd::xtensor3debug(a,"function abs");
    xtensor3<utype> c(a);
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
        for(int k=a.lindex[3];k<=a.uindex[3];k++)
          c[i][j][k]=(utype) aurostd::abs(a[i][j][k]);
    return (xtensor3<utype>) c;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype>                            // function abs xtensor4<>
  xtensor4<utype> abs(const xtensor4<utype>& a) {
    aurostd::xtensor4debug(a,"function abs");
    xtensor4<utype> c(a);
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
        for(int k=a.lindex[3];k<=a.uindex[3];k++)
          for(int l=a.lindex[4];l<=a.uindex[4];l++)
            c[i][j][k][l]=(utype) aurostd::abs(a[i][j][k][l]);
    return (xtensor4<utype>) c;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype>                            // function abs xtensor5<>
  xtensor5<utype> abs(const xtensor5<utype>& a) {
    aurostd::xtensor5debug(a,"function abs");
    xtensor5<utype> c(a);
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
        for(int k=a.lindex[3];k<=a.uindex[3];k++)
          for(int l=a.lindex[4];l<=a.uindex[4];l++)
            for(int m=a.lindex[5];m<=a.uindex[5];m++)
              c[i][j][k][l][m]=(utype) aurostd::abs(a[i][j][k][l][m]);
    return (xtensor5<utype>) c;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype>                            // function abs xtensor6<>
  xtensor6<utype> abs(const xtensor6<utype>& a) {
    aurostd::xtensor6debug(a,"function abs");
    xtensor6<utype> c(a);
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
        for(int k=a.lindex[3];k<=a.uindex[3];k++)
          for(int l=a.lindex[4];l<=a.uindex[4];l++)
            for(int m=a.lindex[5];m<=a.uindex[5];m++)
              for(int n=a.lindex[6];n<=a.uindex[6];n++)
                c[i][j][k][l][m][n]=(utype) aurostd::abs(a[i][j][k][l][m][n]);
    return (xtensor6<utype>) c;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype>                            // function abs xtensor7<>
  xtensor7<utype> abs(const xtensor7<utype>& a) {
    aurostd::xtensor7debug(a,"function abs");
    xtensor7<utype> c(a);
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
        for(int k=a.lindex[3];k<=a.uindex[3];k++)
          for(int l=a.lindex[4];l<=a.uindex[4];l++)
            for(int m=a.lindex[5];m<=a.uindex[5];m++)
              for(int n=a.lindex[6];n<=a.uindex[6];n++)
                for(int o=a.lindex[7];o<=a.uindex[7];o++)
                  c[i][j][k][l][m][n][o]=(utype) aurostd::abs(a[i][j][k][l][m][n][o]);
    return (xtensor7<utype>) c;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype>                            // function abs xtensor8<>
  xtensor8<utype> abs(const xtensor8<utype>& a) {
    aurostd::xtensor8debug(a,"function abs");
    xtensor8<utype> c(a);
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
        for(int k=a.lindex[3];k<=a.uindex[3];k++)
          for(int l=a.lindex[4];l<=a.uindex[4];l++)
            for(int m=a.lindex[5];m<=a.uindex[5];m++)
              for(int n=a.lindex[6];n<=a.uindex[6];n++)
                for(int o=a.lindex[7];o<=a.uindex[7];o++)
                  for(int p=a.lindex[8];p<=a.uindex[8];p++)
                    c[i][j][k][l][m][n][o][p]=(utype) aurostd::abs(a[i][j][k][l][m][n][o][p]);
    return (xtensor8<utype>) c;
  }
}


// ****************************************************************************
// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype>                            // function floor xtensor3<>
  xtensor3<utype> floor(const xtensor3<utype>& a) {
    aurostd::xtensor3debug(a,"function floor");
    xtensor3<utype> c(a);
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
        for(int k=a.lindex[3];k<=a.uindex[3];k++)
          c[i][j][k]=(utype) std::floor(a[i][j][k]);
    return (xtensor3<utype>) c;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype>                            // function floor xtensor4<>
  xtensor4<utype> floor(const xtensor4<utype>& a) {
    aurostd::xtensor4debug(a,"function floor");
    xtensor4<utype> c(a);
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
        for(int k=a.lindex[3];k<=a.uindex[3];k++)
          for(int l=a.lindex[4];l<=a.uindex[4];l++)
            c[i][j][k][l]=(utype) std::floor(a[i][j][k][l]);
    return (xtensor4<utype>) c;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype>                            // function floor xtensor5<>
  xtensor5<utype> floor(const xtensor5<utype>& a) {
    aurostd::xtensor5debug(a,"function floor");
    xtensor5<utype> c(a);
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
        for(int k=a.lindex[3];k<=a.uindex[3];k++)
          for(int l=a.lindex[4];l<=a.uindex[4];l++)
            for(int m=a.lindex[5];m<=a.uindex[5];m++)
              c[i][j][k][l][m]=(utype) std::floor(a[i][j][k][l][m]);
    return (xtensor5<utype>) c;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype>                            // function floor xtensor6<>
  xtensor6<utype> floor(const xtensor6<utype>& a) {
    aurostd::xtensor6debug(a,"function floor");
    xtensor6<utype> c(a);
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
        for(int k=a.lindex[3];k<=a.uindex[3];k++)
          for(int l=a.lindex[4];l<=a.uindex[4];l++)
            for(int m=a.lindex[5];m<=a.uindex[5];m++)
              for(int n=a.lindex[6];n<=a.uindex[6];n++)
                c[i][j][k][l][m][n]=(utype) std::floor(a[i][j][k][l][m][n]);
    return (xtensor6<utype>) c;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype>                            // function floor xtensor7<>
  xtensor7<utype> floor(const xtensor7<utype>& a) {
    aurostd::xtensor7debug(a,"function floor");
    xtensor7<utype> c(a);
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
        for(int k=a.lindex[3];k<=a.uindex[3];k++)
          for(int l=a.lindex[4];l<=a.uindex[4];l++)
            for(int m=a.lindex[5];m<=a.uindex[5];m++)
              for(int n=a.lindex[6];n<=a.uindex[6];n++)
                for(int o=a.lindex[7];o<=a.uindex[7];o++)
                  c[i][j][k][l][m][n][o]=(utype) std::floor(a[i][j][k][l][m][n][o]);
    return (xtensor7<utype>) c;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype>                            // function floor xtensor8<>
  xtensor8<utype> floor(const xtensor8<utype>& a) {
    aurostd::xtensor8debug(a,"function floor");
    xtensor8<utype> c(a);
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
        for(int k=a.lindex[3];k<=a.uindex[3];k++)
          for(int l=a.lindex[4];l<=a.uindex[4];l++)
            for(int m=a.lindex[5];m<=a.uindex[5];m++)
              for(int n=a.lindex[6];n<=a.uindex[6];n++)
                for(int o=a.lindex[7];o<=a.uindex[7];o++)
                  for(int p=a.lindex[8];p<=a.uindex[8];p++)
                    c[i][j][k][l][m][n][o][p]=(utype) std::floor(a[i][j][k][l][m][n][o][p]);
    return (xtensor8<utype>) c;
  }
}

// ****************************************************************************
// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype>                            // function round xtensor3<>
  xtensor3<utype> round(const xtensor3<utype>& a) {
    aurostd::xtensor3debug(a,"function round");
    xtensor3<utype> c(a);
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
        for(int k=a.lindex[3];k<=a.uindex[3];k++)
          c[i][j][k]=(utype) round(a[i][j][k]);
    return (xtensor3<utype>) c;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype>                            // function round xtensor4<>
  xtensor4<utype> round(const xtensor4<utype>& a) {
    aurostd::xtensor4debug(a,"function round");
    xtensor4<utype> c(a);
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
        for(int k=a.lindex[3];k<=a.uindex[3];k++)
          for(int l=a.lindex[4];l<=a.uindex[4];l++)
            c[i][j][k][l]=(utype) round(a[i][j][k][l]);
    return (xtensor4<utype>) c;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype>                            // function round xtensor5<>
  xtensor5<utype> round(const xtensor5<utype>& a) {
    aurostd::xtensor5debug(a,"function round");
    xtensor5<utype> c(a);
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
        for(int k=a.lindex[3];k<=a.uindex[3];k++)
          for(int l=a.lindex[4];l<=a.uindex[4];l++)
            for(int m=a.lindex[5];m<=a.uindex[5];m++)
              c[i][j][k][l][m]=(utype) round(a[i][j][k][l][m]);
    return (xtensor5<utype>) c;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype>                            // function round xtensor6<>
  xtensor6<utype> round(const xtensor6<utype>& a) {
    aurostd::xtensor6debug(a,"function round");
    xtensor6<utype> c(a);
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
        for(int k=a.lindex[3];k<=a.uindex[3];k++)
          for(int l=a.lindex[4];l<=a.uindex[4];l++)
            for(int m=a.lindex[5];m<=a.uindex[5];m++)
              for(int n=a.lindex[6];n<=a.uindex[6];n++)
                c[i][j][k][l][m][n]=(utype) round(a[i][j][k][l][m][n]);
    return (xtensor6<utype>) c;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype>                            // function round xtensor7<>
  xtensor7<utype> round(const xtensor7<utype>& a) {
    aurostd::xtensor7debug(a,"function round");
    xtensor7<utype> c(a);
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
        for(int k=a.lindex[3];k<=a.uindex[3];k++)
          for(int l=a.lindex[4];l<=a.uindex[4];l++)
            for(int m=a.lindex[5];m<=a.uindex[5];m++)
              for(int n=a.lindex[6];n<=a.uindex[6];n++)
                for(int o=a.lindex[7];o<=a.uindex[7];o++)
                  c[i][j][k][l][m][n][o]=(utype) round(a[i][j][k][l][m][n][o]);
    return (xtensor7<utype>) c;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype>                            // function round xtensor8<>
  xtensor8<utype> round(const xtensor8<utype>& a) {
    aurostd::xtensor8debug(a,"function round");
    xtensor8<utype> c(a);
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
        for(int k=a.lindex[3];k<=a.uindex[3];k++)
          for(int l=a.lindex[4];l<=a.uindex[4];l++)
            for(int m=a.lindex[5];m<=a.uindex[5];m++)
              for(int n=a.lindex[6];n<=a.uindex[6];n++)
                for(int o=a.lindex[7];o<=a.uindex[7];o++)
                  for(int p=a.lindex[8];p<=a.uindex[8];p++)
                    c[i][j][k][l][m][n][o][p]=(utype) round(a[i][j][k][l][m][n][o][p]);
    return (xtensor8<utype>) c;
  }
}


// ****************************************************************************
// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype>                            // function ceil xtensor3<>
  xtensor3<utype> ceil(const xtensor3<utype>& a) {
    aurostd::xtensor3debug(a,"function ceil");
    xtensor3<utype> c(a);
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
        for(int k=a.lindex[3];k<=a.uindex[3];k++)
          c[i][j][k]=(utype) std::ceil(a[i][j][k]);
    return (xtensor3<utype>) c;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype>                            // function ceil xtensor4<>
  xtensor4<utype> ceil(const xtensor4<utype>& a) {
    aurostd::xtensor4debug(a,"function ceil");
    xtensor4<utype> c(a);
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
        for(int k=a.lindex[3];k<=a.uindex[3];k++)
          for(int l=a.lindex[4];l<=a.uindex[4];l++)
            c[i][j][k][l]=(utype) std::ceil(a[i][j][k][l]);
    return (xtensor4<utype>) c;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype>                            // function ceil xtensor5<>
  xtensor5<utype> ceil(const xtensor5<utype>& a) {
    aurostd::xtensor5debug(a,"function ceil");
    xtensor5<utype> c(a);
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
        for(int k=a.lindex[3];k<=a.uindex[3];k++)
          for(int l=a.lindex[4];l<=a.uindex[4];l++)
            for(int m=a.lindex[5];m<=a.uindex[5];m++)
              c[i][j][k][l][m]=(utype) std::ceil(a[i][j][k][l][m]);
    return (xtensor5<utype>) c;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype>                            // function ceil xtensor6<>
  xtensor6<utype> ceil(const xtensor6<utype>& a) {
    aurostd::xtensor6debug(a,"function ceil");
    xtensor6<utype> c(a);
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
        for(int k=a.lindex[3];k<=a.uindex[3];k++)
          for(int l=a.lindex[4];l<=a.uindex[4];l++)
            for(int m=a.lindex[5];m<=a.uindex[5];m++)
              for(int n=a.lindex[6];n<=a.uindex[6];n++)
                c[i][j][k][l][m][n]=(utype) std::ceil(a[i][j][k][l][m][n]);
    return (xtensor6<utype>) c;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype>                            // function ceil xtensor7<>
  xtensor7<utype> ceil(const xtensor7<utype>& a) {
    aurostd::xtensor7debug(a,"function ceil");
    xtensor7<utype> c(a);
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
        for(int k=a.lindex[3];k<=a.uindex[3];k++)
          for(int l=a.lindex[4];l<=a.uindex[4];l++)
            for(int m=a.lindex[5];m<=a.uindex[5];m++)
              for(int n=a.lindex[6];n<=a.uindex[6];n++)
                for(int o=a.lindex[7];o<=a.uindex[7];o++)
                  c[i][j][k][l][m][n][o]=(utype) std::ceil(a[i][j][k][l][m][n][o]);
    return (xtensor7<utype>) c;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype>                            // function ceil xtensor8<>
  xtensor8<utype> ceil(const xtensor8<utype>& a) {
    aurostd::xtensor8debug(a,"function ceil");
    xtensor8<utype> c(a);
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
        for(int k=a.lindex[3];k<=a.uindex[3];k++)
          for(int l=a.lindex[4];l<=a.uindex[4];l++)
            for(int m=a.lindex[5];m<=a.uindex[5];m++)
              for(int n=a.lindex[6];n<=a.uindex[6];n++)
                for(int o=a.lindex[7];o<=a.uindex[7];o++)
                  for(int p=a.lindex[8];p<=a.uindex[8];p++)
                    c[i][j][k][l][m][n][o][p]=(utype) std::ceil(a[i][j][k][l][m][n][o][p]);
    return (xtensor8<utype>) c;
  }
}



// ****************************************************************************
// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype>                            // function trunc xtensor3<>
  xtensor3<utype> trunc(const xtensor3<utype>& a) {
    aurostd::xtensor3debug(a,"function trunc");
    xtensor3<utype> c(a);
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
        for(int k=a.lindex[3];k<=a.uindex[3];k++)
          c[i][j][k]=(utype) trunc(a[i][j][k]);
    return (xtensor3<utype>) c;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype>                            // function trunc xtensor4<>
  xtensor4<utype> trunc(const xtensor4<utype>& a) {
    aurostd::xtensor4debug(a,"function trunc");
    xtensor4<utype> c(a);
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
        for(int k=a.lindex[3];k<=a.uindex[3];k++)
          for(int l=a.lindex[4];l<=a.uindex[4];l++)
            c[i][j][k][l]=(utype) trunc(a[i][j][k][l]);
    return (xtensor4<utype>) c;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype>                            // function trunc xtensor5<>
  xtensor5<utype> trunc(const xtensor5<utype>& a) {
    aurostd::xtensor5debug(a,"function trunc");
    xtensor5<utype> c(a);
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
        for(int k=a.lindex[3];k<=a.uindex[3];k++)
          for(int l=a.lindex[4];l<=a.uindex[4];l++)
            for(int m=a.lindex[5];m<=a.uindex[5];m++)
              c[i][j][k][l][m]=(utype) trunc(a[i][j][k][l][m]);
    return (xtensor5<utype>) c;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype>                            // function trunc xtensor6<>
  xtensor6<utype> trunc(const xtensor6<utype>& a) {
    aurostd::xtensor6debug(a,"function trunc");
    xtensor6<utype> c(a);
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
        for(int k=a.lindex[3];k<=a.uindex[3];k++)
          for(int l=a.lindex[4];l<=a.uindex[4];l++)
            for(int m=a.lindex[5];m<=a.uindex[5];m++)
              for(int n=a.lindex[6];n<=a.uindex[6];n++)
                c[i][j][k][l][m][n]=(utype) trunc(a[i][j][k][l][m][n]);
    return (xtensor6<utype>) c;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype>                            // function trunc xtensor7<>
  xtensor7<utype> trunc(const xtensor7<utype>& a) {
    aurostd::xtensor7debug(a,"function trunc");
    xtensor7<utype> c(a);
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
        for(int k=a.lindex[3];k<=a.uindex[3];k++)
          for(int l=a.lindex[4];l<=a.uindex[4];l++)
            for(int m=a.lindex[5];m<=a.uindex[5];m++)
              for(int n=a.lindex[6];n<=a.uindex[6];n++)
                for(int o=a.lindex[7];o<=a.uindex[7];o++)
                  c[i][j][k][l][m][n][o]=(utype) trunc(a[i][j][k][l][m][n][o]);
    return (xtensor7<utype>) c;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype>                            // function trunc xtensor8<>
  xtensor8<utype> trunc(const xtensor8<utype>& a) {
    aurostd::xtensor8debug(a,"function trunc");
    xtensor8<utype> c(a);
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
        for(int k=a.lindex[3];k<=a.uindex[3];k++)
          for(int l=a.lindex[4];l<=a.uindex[4];l++)
            for(int m=a.lindex[5];m<=a.uindex[5];m++)
              for(int n=a.lindex[6];n<=a.uindex[6];n++)
                for(int o=a.lindex[7];o<=a.uindex[7];o++)
                  for(int p=a.lindex[8];p<=a.uindex[8];p++)
                    c[i][j][k][l][m][n][o][p]=(utype) trunc(a[i][j][k][l][m][n][o][p]);
    return (xtensor8<utype>) c;
  }
}


// ****************************************************************************
// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype>                            // function reset xtensor3<>
  void xtensor3<utype>::reset(void) {
    aurostd::xtensor3debug(*this,"function reset");
    for(int i=lindex[1];i<=uindex[1];i++)
      for(int j=lindex[2];j<=uindex[2];j++)
	for(int k=lindex[3];k<=uindex[3];k++)
	  corpus[i][j][k]=(utype) 0.0;
  }
  template<class utype>                            // function reset xtensor3<>
  void reset(const xtensor3<utype>& a) {
    aurostd::xtensor3debug(a,"function reset");
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
	for(int k=a.lindex[3];k<=a.uindex[3];k++)
	  a[i][j][k]=(utype) 0.0;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype>                            // function reset xtensor4<>
  void xtensor4<utype>::reset(void) {
    aurostd::xtensor4debug(*this,"function reset");
    for(int i=lindex[1];i<=uindex[1];i++)
      for(int j=lindex[2];j<=uindex[2];j++)
	for(int k=lindex[3];k<=uindex[3];k++)
	  for(int l=lindex[4];l<=uindex[4];l++)
	    corpus[i][j][k][l]=(utype) 0.0;
  }
  template<class utype>                            // function reset xtensor4<>
  void reset(const xtensor4<utype>& a) {
    aurostd::xtensor4debug(a,"function reset");
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
	for(int k=a.lindex[3];k<=a.uindex[3];k++)
	  for(int l=a.lindex[4];l<=a.uindex[4];l++)
	    a[i][j][k][l]=(utype) 0.0;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype>                            // function reset xtensor5<>
  void xtensor5<utype>::reset(void) {
    aurostd::xtensor5debug(*this,"function reset");
    for(int i=lindex[1];i<=uindex[1];i++)
      for(int j=lindex[2];j<=uindex[2];j++)
	for(int k=lindex[3];k<=uindex[3];k++)
	  for(int l=lindex[4];l<=uindex[4];l++)
	    for(int m=lindex[5];m<=uindex[5];m++)
	      corpus[i][j][k][l][m]=(utype) 0.0;
  }
  template<class utype>                            // function reset xtensor5<>
  void reset(const xtensor5<utype>& a) {
    aurostd::xtensor5debug(a,"function reset");
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
	for(int k=a.lindex[3];k<=a.uindex[3];k++)
	  for(int l=a.lindex[4];l<=a.uindex[4];l++)
	    for(int m=a.lindex[5];m<=a.uindex[5];m++)
	      a[i][j][k][l][m]=(utype) 0.0;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype>                            // function reset xtensor6<>
  void xtensor6<utype>::reset(void) {
    aurostd::xtensor6debug(*this,"function reset");
    for(int i=lindex[1];i<=uindex[1];i++)
      for(int j=lindex[2];j<=uindex[2];j++)
	for(int k=lindex[3];k<=uindex[3];k++)
	  for(int l=lindex[4];l<=uindex[4];l++)
	    for(int m=lindex[5];m<=uindex[5];m++)
	      for(int n=lindex[6];n<=uindex[6];n++)
		corpus[i][j][k][l][m][n]=(utype) 0.0;
  }
  template<class utype>                            // function reset xtensor6<>
  void reset(const xtensor6<utype>& a) {
    aurostd::xtensor6debug(a,"function reset");
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
	for(int k=a.lindex[3];k<=a.uindex[3];k++)
	  for(int l=a.lindex[4];l<=a.uindex[4];l++)
	    for(int m=a.lindex[5];m<=a.uindex[5];m++)
	      for(int n=a.lindex[6];n<=a.uindex[6];n++)
		a[i][j][k][l][m][n]=(utype) 0.0;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype>                            // function reset xtensor7<>
  void xtensor7<utype>::reset(void) {
    aurostd::xtensor7debug(*this,"function reset");
    for(int i=lindex[1];i<=uindex[1];i++)
      for(int j=lindex[2];j<=uindex[2];j++)
	for(int k=lindex[3];k<=uindex[3];k++)
	  for(int l=lindex[4];l<=uindex[4];l++)
	    for(int m=lindex[5];m<=uindex[5];m++)
	      for(int n=lindex[6];n<=uindex[6];n++)
		for(int o=lindex[7];o<=uindex[7];o++)
		  corpus[i][j][k][l][m][n][o]=(utype) 0.0;
  }
  template<class utype>                            // function reset xtensor7<>
  void reset(const xtensor7<utype>& a) {
    aurostd::xtensor7debug(a,"function reset");
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
	for(int k=a.lindex[3];k<=a.uindex[3];k++)
	  for(int l=a.lindex[4];l<=a.uindex[4];l++)
	    for(int m=a.lindex[5];m<=a.uindex[5];m++)
	      for(int n=a.lindex[6];n<=a.uindex[6];n++)
		for(int o=a.lindex[7];o<=a.uindex[7];o++)
		  a[i][j][k][l][m][n][o]=(utype) 0.0;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype>                            // function reset xtensor8<>
  void xtensor8<utype>::reset(void) {
    aurostd::xtensor8debug(*this,"function reset");
    for(int i=lindex[1];i<=uindex[1];i++)
      for(int j=lindex[2];j<=uindex[2];j++)
	for(int k=lindex[3];k<=uindex[3];k++)
	  for(int l=lindex[4];l<=uindex[4];l++)
	    for(int m=lindex[5];m<=uindex[5];m++)
	      for(int n=lindex[6];n<=uindex[6];n++)
		for(int o=lindex[7];o<=uindex[7];o++)
		  for(int p=lindex[8];p<=uindex[8];p++)
		    corpus[i][j][k][l][m][n][o][p]=(utype) 0.0;
  }
  template<class utype>                            // function reset xtensor8<>
  void reset(const xtensor8<utype>& a) {
    aurostd::xtensor8debug(a,"function reset");
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
	for(int k=a.lindex[3];k<=a.uindex[3];k++)
	  for(int l=a.lindex[4];l<=a.uindex[4];l++)
	    for(int m=a.lindex[5];m<=a.uindex[5];m++)
	      for(int n=a.lindex[6];n<=a.uindex[6];n++)
		for(int o=a.lindex[7];o<=a.uindex[7];o++)
		  for(int p=a.lindex[8];p<=a.uindex[8];p++)
		    a[i][j][k][l][m][n][o][p]=(utype) 0.0;
  }
}

// ****************************************************************************
// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype>                            // function clear xtensor3<>
  void xtensor3<utype>::clear(void) {
    aurostd::xtensor3debug(*this,"function clear");
    for(int i=lindex[1];i<=uindex[1];i++)
      for(int j=lindex[2];j<=uindex[2];j++)
	for(int k=lindex[3];k<=uindex[3];k++)
	  corpus[i][j][k]=(utype) 0.0;
  }
  template<class utype>                            // function clear xtensor3<>
  void clear(const xtensor3<utype>& a) {
    aurostd::xtensor3debug(a,"function clear");
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
	for(int k=a.lindex[3];k<=a.uindex[3];k++)
	  a[i][j][k]=(utype) 0.0;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype>                            // function clear xtensor4<>
  void xtensor4<utype>::clear(void) {
    aurostd::xtensor4debug(*this,"function clear");
    for(int i=lindex[1];i<=uindex[1];i++)
      for(int j=lindex[2];j<=uindex[2];j++)
	for(int k=lindex[3];k<=uindex[3];k++)
	  for(int l=lindex[4];l<=uindex[4];l++)
	    corpus[i][j][k][l]=(utype) 0.0;
  }
  template<class utype>                            // function clear xtensor4<>
  void clear(const xtensor4<utype>& a) {
    aurostd::xtensor4debug(a,"function clear");
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
	for(int k=a.lindex[3];k<=a.uindex[3];k++)
	  for(int l=a.lindex[4];l<=a.uindex[4];l++)
	    a[i][j][k][l]=(utype) 0.0;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype>                            // function clear xtensor5<>
  void xtensor5<utype>::clear(void) {
    aurostd::xtensor5debug(*this,"function clear");
    for(int i=lindex[1];i<=uindex[1];i++)
      for(int j=lindex[2];j<=uindex[2];j++)
	for(int k=lindex[3];k<=uindex[3];k++)
	  for(int l=lindex[4];l<=uindex[4];l++)
	    for(int m=lindex[5];m<=uindex[5];m++)
	      corpus[i][j][k][l][m]=(utype) 0.0;
  }
  template<class utype>                            // function clear xtensor5<>
  void clear(const xtensor5<utype>& a) {
    aurostd::xtensor5debug(a,"function clear");
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
	for(int k=a.lindex[3];k<=a.uindex[3];k++)
	  for(int l=a.lindex[4];l<=a.uindex[4];l++)
	    for(int m=a.lindex[5];m<=a.uindex[5];m++)
	      a[i][j][k][l][m]=(utype) 0.0;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype>                            // function clear xtensor6<>
  void xtensor6<utype>::clear(void) {
    aurostd::xtensor6debug(*this,"function clear");
    for(int i=lindex[1];i<=uindex[1];i++)
      for(int j=lindex[2];j<=uindex[2];j++)
	for(int k=lindex[3];k<=uindex[3];k++)
	  for(int l=lindex[4];l<=uindex[4];l++)
	    for(int m=lindex[5];m<=uindex[5];m++)
	      for(int n=lindex[6];n<=uindex[6];n++)
		corpus[i][j][k][l][m][n]=(utype) 0.0;
  }
  template<class utype>                            // function clear xtensor6<>
  void clear(const xtensor6<utype>& a) {
    aurostd::xtensor6debug(a,"function clear");
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
	for(int k=a.lindex[3];k<=a.uindex[3];k++)
	  for(int l=a.lindex[4];l<=a.uindex[4];l++)
	    for(int m=a.lindex[5];m<=a.uindex[5];m++)
	      for(int n=a.lindex[6];n<=a.uindex[6];n++)
		a[i][j][k][l][m][n]=(utype) 0.0;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype>                            // function clear xtensor7<>
  void xtensor7<utype>::clear(void) {
    aurostd::xtensor7debug(*this,"function clear");
    for(int i=lindex[1];i<=uindex[1];i++)
      for(int j=lindex[2];j<=uindex[2];j++)
	for(int k=lindex[3];k<=uindex[3];k++)
	  for(int l=lindex[4];l<=uindex[4];l++)
	    for(int m=lindex[5];m<=uindex[5];m++)
	      for(int n=lindex[6];n<=uindex[6];n++)
		for(int o=lindex[7];o<=uindex[7];o++)
		  corpus[i][j][k][l][m][n][o]=(utype) 0.0;
  }
  template<class utype>                            // function clear xtensor7<>
  void clear(const xtensor7<utype>& a) {
    aurostd::xtensor7debug(a,"function clear");
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
	for(int k=a.lindex[3];k<=a.uindex[3];k++)
	  for(int l=a.lindex[4];l<=a.uindex[4];l++)
	    for(int m=a.lindex[5];m<=a.uindex[5];m++)
	      for(int n=a.lindex[6];n<=a.uindex[6];n++)
		for(int o=a.lindex[7];o<=a.uindex[7];o++)
		  a[i][j][k][l][m][n][o]=(utype) 0.0;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype>                            // function clear xtensor8<>
  void xtensor8<utype>::clear(void) {
    aurostd::xtensor8debug(*this,"function clear");
    for(int i=lindex[1];i<=uindex[1];i++)
      for(int j=lindex[2];j<=uindex[2];j++)
	for(int k=lindex[3];k<=uindex[3];k++)
	  for(int l=lindex[4];l<=uindex[4];l++)
	    for(int m=lindex[5];m<=uindex[5];m++)
	      for(int n=lindex[6];n<=uindex[6];n++)
		for(int o=lindex[7];o<=uindex[7];o++)
		  for(int p=lindex[8];p<=uindex[8];p++)
		    corpus[i][j][k][l][m][n][o][p]=(utype) 0.0;
  }
  template<class utype>                            // function clear xtensor8<>
  void clear(const xtensor8<utype>& a) {
    aurostd::xtensor8debug(a,"function clear");
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
	for(int k=a.lindex[3];k<=a.uindex[3];k++)
	  for(int l=a.lindex[4];l<=a.uindex[4];l++)
	    for(int m=a.lindex[5];m<=a.uindex[5];m++)
	      for(int n=a.lindex[6];n<=a.uindex[6];n++)
		for(int o=a.lindex[7];o<=a.uindex[7];o++)
		  for(int p=a.lindex[8];p<=a.uindex[8];p++)
		    a[i][j][k][l][m][n][o][p]=(utype) 0.0;
  }
}

// ****************************************************************************
// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype>                            // function set xtensor3<>
  void xtensor3<utype>::set(const utype& s) {
    aurostd::xtensor3debug(*this,"function set");
    for(int i=lindex[1];i<=uindex[1];i++)
      for(int j=lindex[2];j<=uindex[2];j++)
	for(int k=lindex[3];k<=uindex[3];k++)
	  corpus[i][j][k]=(utype) s;
  }
  template<class utype>                              // function set xtensor3<>
  void set(const xtensor3<utype>& a,const utype& s) {
    aurostd::xtensor3debug(a,"function set");
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
	for(int k=a.lindex[3];k<=a.uindex[3];k++)
	  a[i][j][k]=(utype) s;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype>                            // function set xtensor4<>
  void xtensor4<utype>::set(const utype& s) {
    aurostd::xtensor4debug(*this,"function set");
    for(int i=lindex[1];i<=uindex[1];i++)
      for(int j=lindex[2];j<=uindex[2];j++)
	for(int k=lindex[3];k<=uindex[3];k++)
	  for(int l=lindex[4];l<=uindex[4];l++)
	    corpus[i][j][k][l]=(utype) s;
  }
  template<class utype>                              // function set xtensor4<>
  void set(const xtensor4<utype>& a,const utype& s) {
    aurostd::xtensor4debug(a,"function set");
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
	for(int k=a.lindex[3];k<=a.uindex[3];k++)
	  for(int l=a.lindex[4];l<=a.uindex[4];l++)
	    a[i][j][k][l]=(utype) s;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype>                            // function set xtensor5<>
  void xtensor5<utype>::set(const utype& s) {
    aurostd::xtensor5debug(*this,"function set");
    for(int i=lindex[1];i<=uindex[1];i++)
      for(int j=lindex[2];j<=uindex[2];j++)
	for(int k=lindex[3];k<=uindex[3];k++)
	  for(int l=lindex[4];l<=uindex[4];l++)
	    for(int m=lindex[5];m<=uindex[5];m++)
	      corpus[i][j][k][l][m]=(utype) s;
  }
  template<class utype>                              // function set xtensor5<>
  void set(const xtensor5<utype>& a,const utype& s) {
    aurostd::xtensor5debug(a,"function set");
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
	for(int k=a.lindex[3];k<=a.uindex[3];k++)
	  for(int l=a.lindex[4];l<=a.uindex[4];l++)
	    for(int m=a.lindex[5];m<=a.uindex[5];m++)
	      a[i][j][k][l][m]=(utype) s;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype>                            // function set xtensor6<>
  void xtensor6<utype>::set(const utype& s) {
    aurostd::xtensor6debug(*this,"function set");
    for(int i=lindex[1];i<=uindex[1];i++)
      for(int j=lindex[2];j<=uindex[2];j++)
	for(int k=lindex[3];k<=uindex[3];k++)
	  for(int l=lindex[4];l<=uindex[4];l++)
	    for(int m=lindex[5];m<=uindex[5];m++)
	      for(int n=lindex[6];n<=uindex[6];n++)
		corpus[i][j][k][l][m][n]=(utype) s;
  }
  template<class utype>                              // function set xtensor6<>
  void set(const xtensor6<utype>& a,const utype& s) {
    aurostd::xtensor6debug(a,"function set");
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
	for(int k=a.lindex[3];k<=a.uindex[3];k++)
	  for(int l=a.lindex[4];l<=a.uindex[4];l++)
	    for(int m=a.lindex[5];m<=a.uindex[5];m++)
	      for(int n=a.lindex[6];n<=a.uindex[6];n++)
		a[i][j][k][l][m][n]=(utype) s;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype>                            // function set xtensor7<>
  void xtensor7<utype>::set(const utype& s) {
    aurostd::xtensor7debug(*this,"function set");
    for(int i=lindex[1];i<=uindex[1];i++)
      for(int j=lindex[2];j<=uindex[2];j++)
	for(int k=lindex[3];k<=uindex[3];k++)
	  for(int l=lindex[4];l<=uindex[4];l++)
	    for(int m=lindex[5];m<=uindex[5];m++)
	      for(int n=lindex[6];n<=uindex[6];n++)
		for(int o=lindex[7];o<=uindex[7];o++)
		  corpus[i][j][k][l][m][n][o]=(utype) s;
  }
  template<class utype>                              // function set xtensor7<>
  void set(const xtensor7<utype>& a,const utype& s) {
    aurostd::xtensor7debug(a,"function set");
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
	for(int k=a.lindex[3];k<=a.uindex[3];k++)
	  for(int l=a.lindex[4];l<=a.uindex[4];l++)
	    for(int m=a.lindex[5];m<=a.uindex[5];m++)
	      for(int n=a.lindex[6];n<=a.uindex[6];n++)
		for(int o=a.lindex[7];o<=a.uindex[7];o++)
		  a[i][j][k][l][m][n][o]=(utype) s;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype>                            // function set xtensor8<>
  void xtensor8<utype>::set(const utype& s) {
    aurostd::xtensor8debug(*this,"function set");
    for(int i=lindex[1];i<=uindex[1];i++)
      for(int j=lindex[2];j<=uindex[2];j++)
	for(int k=lindex[3];k<=uindex[3];k++)
	  for(int l=lindex[4];l<=uindex[4];l++)
	    for(int m=lindex[5];m<=uindex[5];m++)
	      for(int n=lindex[6];n<=uindex[6];n++)
		for(int o=lindex[7];o<=uindex[7];o++)
		  for(int p=lindex[8];p<=uindex[8];p++)
		    corpus[i][j][k][l][m][n][o][p]=(utype) s;
  }
  template<class utype>                              // function set xtensor8<>
  void set(const xtensor8<utype>& a,const utype& s) {
    aurostd::xtensor8debug(a,"function set");
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
	for(int k=a.lindex[3];k<=a.uindex[3];k++)
	  for(int l=a.lindex[4];l<=a.uindex[4];l++)
	    for(int m=a.lindex[5];m<=a.uindex[5];m++)
	      for(int n=a.lindex[6];n<=a.uindex[6];n++)
		for(int o=a.lindex[7];o<=a.uindex[7];o++)
		  for(int p=a.lindex[8];p<=a.uindex[8];p++)
		    a[i][j][k][l][m][n][o][p]=(utype) s;
  }
}

// ****************************************************************************
// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype>                              // function sum xtensor3<>
  utype sum(const xtensor3<utype>& a) {
    aurostd::xtensor3debug(a,"function sum");
    utype c=utype(0.0);
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
	for(int k=a.lindex[3];k<=a.uindex[3];k++)
	  c+=a[i][j][k];
    return c;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype>                              // function sum xtensor4<>
  utype sum(const xtensor4<utype>& a) {
    aurostd::xtensor4debug(a,"function sum");
    utype c=utype(0.0);
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
	for(int k=a.lindex[3];k<=a.uindex[3];k++)
	  for(int l=a.lindex[4];l<=a.uindex[4];l++)
	    c+=a[i][j][k][l];
    return c;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype>                              // function sum xtensor5<>
  utype sum(const xtensor5<utype>& a) {
    aurostd::xtensor5debug(a,"function sum");
    utype c=utype(0.0);
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
	for(int k=a.lindex[3];k<=a.uindex[3];k++)
	  for(int l=a.lindex[4];l<=a.uindex[4];l++)
	    for(int m=a.lindex[5];m<=a.uindex[5];m++)
	      c+=a[i][j][k][l][m];
    return c;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype>                              // function sum xtensor6<>
  utype sum(const xtensor6<utype>& a) {
    aurostd::xtensor6debug(a,"function sum");
    utype c=utype(0.0);
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
	for(int k=a.lindex[3];k<=a.uindex[3];k++)
	  for(int l=a.lindex[4];l<=a.uindex[4];l++)
	    for(int m=a.lindex[5];m<=a.uindex[5];m++)
	      for(int n=a.lindex[6];n<=a.uindex[6];n++)
		c+=a[i][j][k][l][m][n];
    return c;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype>                              // function sum xtensor7<>
  utype sum(const xtensor7<utype>& a) {
    aurostd::xtensor7debug(a,"function sum");
    utype c=utype(0.0);
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
	for(int k=a.lindex[3];k<=a.uindex[3];k++)
	  for(int l=a.lindex[4];l<=a.uindex[4];l++)
	    for(int m=a.lindex[5];m<=a.uindex[5];m++)
	      for(int n=a.lindex[6];n<=a.uindex[6];n++)
		for(int o=a.lindex[7];o<=a.uindex[7];o++)
		  c+=a[i][j][k][l][m][n][o];
    return c;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype>                              // function sum xtensor8<>
  utype sum(const xtensor8<utype>& a) {
    aurostd::xtensor8debug(a,"function sum");
    utype c=utype(0.0);
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
	for(int k=a.lindex[3];k<=a.uindex[3];k++)
	  for(int l=a.lindex[4];l<=a.uindex[4];l++)
	    for(int m=a.lindex[5];m<=a.uindex[5];m++)
	      for(int n=a.lindex[6];n<=a.uindex[6];n++)
		for(int o=a.lindex[7];o<=a.uindex[7];o++)
		  for(int p=a.lindex[8];p<=a.uindex[8];p++)
		    c+=a[i][j][k][l][m][n][o][p];
    return c;
  }
}

// ****************************************************************************
// ------------------------------------------------------ functions of xtensor3
namespace aurostd {
  template<class utype>                              // function min xtensor3<>
  utype min(const xtensor3<utype>& a) {
    aurostd::xtensor3debug(a,"function min");
    utype c=a[a.lindex[1]][a.lindex[2]][a.lindex[3]];
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
	for(int k=a.lindex[3];k<=a.uindex[3];k++)
	  if(a[i][j][k]<c) c=a[i][j][k];
    return c;
  }
}

// ------------------------------------------------------ functions of xtensor4
namespace aurostd {
  template<class utype>                              // function min xtensor4<>
  utype min(const xtensor4<utype>& a) {
    aurostd::xtensor4debug(a,"function min");
    utype c=a[a.lindex[1]][a.lindex[2]][a.lindex[3]][a.lindex[4]];
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
	for(int k=a.lindex[3];k<=a.uindex[3];k++)
	  for(int l=a.lindex[4];l<=a.uindex[4];l++)
	    if(a[i][j][k][l]<c) c=a[i][j][k][l];
    return c;
  }
}

// ------------------------------------------------------ functions of xtensor5
namespace aurostd {
  template<class utype>                              // function min xtensor5<>
  utype min(const xtensor5<utype>& a) {
    aurostd::xtensor5debug(a,"function min");
    utype c=a[a.lindex[1]][a.lindex[2]][a.lindex[3]][a.lindex[4]][a.lindex[5]];
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
	for(int k=a.lindex[3];k<=a.uindex[3];k++)
	  for(int l=a.lindex[4];l<=a.uindex[4];l++)
	    for(int m=a.lindex[5];m<=a.uindex[5];m++)
	      if(a[i][j][k][l][m]<c) c=a[i][j][k][l][m];
    return c;
  }
}

// ---------------------------------------------------- functions of xtensor6
namespace aurostd {
  template<class utype>                            // function min xtensor6<>
  utype min(const xtensor6<utype>& a) {
    aurostd::xtensor6debug(a,"function min");
    utype c=a[a.lindex[1]][a.lindex[2]][a.lindex[3]][a.lindex[4]][a.lindex[5]][a.lindex[6]];
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
	for(int k=a.lindex[3];k<=a.uindex[3];k++)
	  for(int l=a.lindex[4];l<=a.uindex[4];l++)
	    for(int m=a.lindex[5];m<=a.uindex[5];m++)
	      for(int n=a.lindex[6];n<=a.uindex[6];n++)
		if(a[i][j][k][l][m][n]<c) c=a[i][j][k][l][m][n];
    return c;
  }
}

// ---------------------------------------------------- functions of xtensor7
namespace aurostd {
  template<class utype>                            // function min xtensor7<>
  utype min(const xtensor7<utype>& a) {
    aurostd::xtensor7debug(a,"function min");
    utype c=a[a.lindex[1]][a.lindex[2]][a.lindex[3]][a.lindex[4]][a.lindex[5]][a.lindex[6]][a.lindex[7]];
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
	for(int k=a.lindex[3];k<=a.uindex[3];k++)
	  for(int l=a.lindex[4];l<=a.uindex[4];l++)
	    for(int m=a.lindex[5];m<=a.uindex[5];m++)
	      for(int n=a.lindex[6];n<=a.uindex[6];n++)
		for(int o=a.lindex[7];o<=a.uindex[7];o++)
		  if(a[i][j][k][l][m][n][o]<c) c=a[i][j][k][l][m][n][o];
    return c;
  }
}

// ---------------------------------------------------- functions of xtensor8
namespace aurostd {
  template<class utype>                            // function min xtensor8<>
  utype min(const xtensor8<utype>& a) {
    aurostd::xtensor8debug(a,"function min");
    utype c=a[a.lindex[1]][a.lindex[2]][a.lindex[3]][a.lindex[4]][a.lindex[5]][a.lindex[6]][a.lindex[7]][a.lindex[8]];
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
	for(int k=a.lindex[3];k<=a.uindex[3];k++)
	  for(int l=a.lindex[4];l<=a.uindex[4];l++)
	    for(int m=a.lindex[5];m<=a.uindex[5];m++)
	      for(int n=a.lindex[6];n<=a.uindex[6];n++)
		for(int o=a.lindex[7];o<=a.uindex[7];o++)
		  for(int p=a.lindex[8];p<=a.uindex[8];p++)
		    if(a[i][j][k][l][m][n][o][p]<c) c=a[i][j][k][l][m][n][o][p];
    return c;
  }
}

// **************************************************************************
// --------------------------------------------------------------------------
namespace aurostd {
  template<class utype>                            // function min xtensor3<>
  utype min(const xtensor3<utype>& a,int& ii,int& jj,int& kk) {
    aurostd::xtensor3debug(a,"function min");
    utype c=a[a.lindex[1]][a.lindex[2]][a.lindex[3]];
    ii=a.lindex[1],jj=a.lindex[2],kk=a.lindex[3];
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
	for(int k=a.lindex[3];k<=a.uindex[3];k++)
	  if(a[i][j][k]<c) {
	    c=a[i][j][k];
	    ii=i;jj=j;kk=k;
	  }
    return c;
  }
}

// --------------------------------------------------------------------------
namespace aurostd {
  template<class utype>                            // function min xtensor4<>
  utype min(const xtensor4<utype>& a,int& ii,int& jj,int& kk,int& ll) {
    aurostd::xtensor4debug(a,"function min");
    utype c=a[a.lindex[1]][a.lindex[2]][a.lindex[3]][a.lindex[4]];
    ii=a.lindex[1],jj=a.lindex[2],kk=a.lindex[3],ll=a.lindex[4];
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
	for(int k=a.lindex[3];k<=a.uindex[3];k++)
	  for(int l=a.lindex[4];l<=a.uindex[4];l++)
	    if(a[i][j][k][l]<c) {
	      c=a[i][j][k][l];
	      ii=i;jj=j;kk=k;ll=l;
	    }
    return c;
  }
}

// --------------------------------------------------------------------------
namespace aurostd {
  template<class utype>                            // function min xtensor5<>
  utype min(const xtensor5<utype>& a,int& ii,int& jj,int& kk,int& ll,int& mm) {
    aurostd::xtensor5debug(a,"function min");
    utype c=a[a.lindex[1]][a.lindex[2]][a.lindex[3]][a.lindex[4]][a.lindex[5]];
    ii=a.lindex[1],jj=a.lindex[2],kk=a.lindex[3],ll=a.lindex[4],mm=a.lindex[5];
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
	for(int k=a.lindex[3];k<=a.uindex[3];k++)
	  for(int l=a.lindex[4];l<=a.uindex[4];l++)
	    for(int m=a.lindex[5];m<=a.uindex[5];m++)
	      if(a[i][j][k][l][m]<c) {
		c=a[i][j][k][l][m];
		ii=i;jj=j;kk=k; ll=l;mm=m;
	      }
    return c;
  }
}

// --------------------------------------------------------------------------
namespace aurostd {
  template<class utype>                            // function min xtensor7<>
  utype min(const xtensor7<utype>& a,int& ii,int& jj,int& kk,int& ll,int& mm,int& nn,int& oo) {
    aurostd::xtensor7debug(a,"function min");
    utype c=a[a.lindex[1]][a.lindex[2]][a.lindex[3]][a.lindex[4]][a.lindex[5]][a.lindex[6]][a.lindex[7]];
    ii=a.lindex[1],jj=a.lindex[2],kk=a.lindex[3],ll=a.lindex[4],mm=a.lindex[5],nn=a.lindex[6],oo=a.lindex[7];
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
	for(int k=a.lindex[3];k<=a.uindex[3];k++)
	  for(int l=a.lindex[4];l<=a.uindex[4];l++)
	    for(int m=a.lindex[5];m<=a.uindex[5];m++)
	      for(int n=a.lindex[6];n<=a.uindex[6];n++)
		for(int o=a.lindex[7];o<=a.uindex[7];o++)
		  if(a[i][j][k][l][m][n][o]<c) {
		    c=a[i][j][k][l][m][n][o];
		    ii=i;jj=j;kk=k;ll=l;mm=m;nn=n;oo=o;
		  }
    return c;
  }
}

// --------------------------------------------------------------------------
namespace aurostd {
  template<class utype>                            // function min xtensor8<>
  utype min(const xtensor8<utype>& a,int& ii,int& jj,int& kk,int& ll,int& mm,int& nn,int& oo,int& pp) {
    aurostd::xtensor8debug(a,"function min");
    utype c=a[a.lindex[1]][a.lindex[2]][a.lindex[3]][a.lindex[4]][a.lindex[5]][a.lindex[6]][a.lindex[7]][a.lindex[8]];
    ii=a.lindex[1],jj=a.lindex[2],kk=a.lindex[3],ll=a.lindex[4],mm=a.lindex[5],nn=a.lindex[6],oo=a.lindex[7],pp=a.lindex[8];
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
	for(int k=a.lindex[3];k<=a.uindex[3];k++)
	  for(int l=a.lindex[4];l<=a.uindex[4];l++)
	    for(int m=a.lindex[5];m<=a.uindex[5];m++)
	      for(int n=a.lindex[6];n<=a.uindex[6];n++)
		for(int o=a.lindex[7];o<=a.uindex[7];o++)
		  for(int p=a.lindex[8];p<=a.uindex[8];p++)
		    if(a[i][j][k][l][m][n][o][p]<c) {
		      c=a[i][j][k][l][m][n][o][p];
		      ii=i;jj=j;kk=k;ll=l;mm=m;nn=n;oo=o;pp=p;
		    }
    return c;
  }
}

// **************************************************************************
// ---------------------------------------------------- functions of xtensor3
namespace aurostd {
  template<class utype>                            // function max xtensor3<>
  utype max(const xtensor3<utype>& a) {
    aurostd::xtensor3debug(a,"function max");
    utype c=a[a.lindex[1]][a.lindex[2]][a.lindex[3]];
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
	for(int k=a.lindex[3];k<=a.uindex[3];k++)
	  if(a[i][j][k]>c) c=a[i][j][k];
    return c;
  }
}

// ---------------------------------------------------- functions of xtensor4
namespace aurostd {
  template<class utype>                            // function max xtensor4<>
  utype max(const xtensor4<utype>& a) {
    aurostd::xtensor4debug(a,"function max");
    utype c=a[a.lindex[1]][a.lindex[2]][a.lindex[3]][a.lindex[4]];
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
	for(int k=a.lindex[3];k<=a.uindex[3];k++)
	  for(int l=a.lindex[4];l<=a.uindex[4];l++)
	    if(a[i][j][k][l]>c) c=a[i][j][k][l];
    return c;
  }
}

// ---------------------------------------------------- functions of xtensor5
namespace aurostd {
  template<class utype>                            // function max xtensor5<>
  utype max(const xtensor5<utype>& a) {
    aurostd::xtensor5debug(a,"function max");
    utype c=a[a.lindex[1]][a.lindex[2]][a.lindex[3]][a.lindex[4]][a.lindex[5]];
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
	for(int k=a.lindex[3];k<=a.uindex[3];k++)
	  for(int l=a.lindex[4];l<=a.uindex[4];l++)
	    for(int m=a.lindex[5];m<=a.uindex[5];m++)
	      if(a[i][j][k][l][m]>c) c=a[i][j][k][l][m];
    return c;
  }
}

// ---------------------------------------------------- functions of xtensor6
namespace aurostd {
  template<class utype>                            // function max xtensor6<>
  utype max(const xtensor6<utype>& a) {
    aurostd::xtensor6debug(a,"function max");
    utype c=a[a.lindex[1]][a.lindex[2]][a.lindex[3]][a.lindex[4]][a.lindex[5]][a.lindex[6]];
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
	for(int k=a.lindex[3];k<=a.uindex[3];k++)
	  for(int l=a.lindex[4];l<=a.uindex[4];l++)
	    for(int m=a.lindex[5];m<=a.uindex[5];m++)
	      for(int n=a.lindex[6];n<=a.uindex[6];n++)
		if(a[i][j][k][l][m][n]>c) c=a[i][j][k][l][m][n];
    return c;
  }
}

// ---------------------------------------------------- functions of xtensor7
namespace aurostd {
  template<class utype>                            // function max xtensor7<>
  utype max(const xtensor7<utype>& a) {
    aurostd::xtensor7debug(a,"function max");
    utype c=a[a.lindex[1]][a.lindex[2]][a.lindex[3]][a.lindex[4]][a.lindex[5]][a.lindex[6]][a.lindex[7]];
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
	for(int k=a.lindex[3];k<=a.uindex[3];k++)
	  for(int l=a.lindex[4];l<=a.uindex[4];l++)
	    for(int m=a.lindex[5];m<=a.uindex[5];m++)
	      for(int n=a.lindex[6];n<=a.uindex[6];n++)
		for(int o=a.lindex[7];o<=a.uindex[7];o++)
		  if(a[i][j][k][l][m][n][o]>c) c=a[i][j][k][l][m][n][o];
    return c;
  }
}

// ---------------------------------------------------- functions of xtensor8
namespace aurostd {
  template<class utype>                            // function max xtensor8<>
  utype max(const xtensor8<utype>& a) {
    aurostd::xtensor8debug(a,"function max");
    utype c=a[a.lindex[1]][a.lindex[2]][a.lindex[3]][a.lindex[4]][a.lindex[5]][a.lindex[6]][a.lindex[7]][a.lindex[8]];
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
	for(int k=a.lindex[3];k<=a.uindex[3];k++)
	  for(int l=a.lindex[4];l<=a.uindex[4];l++)
	    for(int m=a.lindex[5];m<=a.uindex[5];m++)
	      for(int n=a.lindex[6];n<=a.uindex[6];n++)
		for(int o=a.lindex[7];o<=a.uindex[7];o++)
		  for(int p=a.lindex[8];p<=a.uindex[8];p++)
		    if(a[i][j][k][l][m][n][o][p]>c) c=a[i][j][k][l][m][n][o][p];
    return c;
  }
}

// **************************************************************************
// --------------------------------------------------------------------------
namespace aurostd {
  template<class utype>                            // function max xtensor3<>
  utype max(const xtensor3<utype>& a,int& ii,int& jj,int& kk) {
    aurostd::xtensor3debug(a,"function max");
    utype c=a[a.lindex[1]][a.lindex[2]][a.lindex[3]];
    ii=a.lindex[1],jj=a.lindex[2],kk=a.lindex[3];
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
	for(int k=a.lindex[3];k<=a.uindex[3];k++)
	  if(a[i][j][k]>c) {
	    c=a[i][j][k];
	    ii=i;jj=j;kk=k;
	  }
    return c;
  }
}

// --------------------------------------------------------------------------
namespace aurostd {
  template<class utype>                            // function max xtensor4<>
  utype max(const xtensor4<utype>& a,int& ii,int& jj,int& kk,int& ll) {
    aurostd::xtensor4debug(a,"function max");
    utype c=a[a.lindex[1]][a.lindex[2]][a.lindex[3]][a.lindex[4]];
    ii=a.lindex[1],jj=a.lindex[2],kk=a.lindex[3],ll=a.lindex[4];
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
	for(int k=a.lindex[3];k<=a.uindex[3];k++)
	  for(int l=a.lindex[4];l<=a.uindex[4];l++)
	    if(a[i][j][k][l]>c) {
	      c=a[i][j][k][l];
	      ii=i;jj=j;kk=k;ll=l;
	    }
    return c;
  }
}

// --------------------------------------------------------------------------
namespace aurostd {
  template<class utype>                            // function max xtensor5<>
  utype max(const xtensor5<utype>& a,int& ii,int& jj,int& kk,int& ll,int& mm) {
    aurostd::xtensor5debug(a,"function max");
    utype c=a[a.lindex[1]][a.lindex[2]][a.lindex[3]][a.lindex[4]][a.lindex[5]];
    ii=a.lindex[1],jj=a.lindex[2],kk=a.lindex[3],ll=a.lindex[4],mm=a.lindex[5];
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
	for(int k=a.lindex[3];k<=a.uindex[3];k++)
	  for(int l=a.lindex[4];l<=a.uindex[4];l++)
	    for(int m=a.lindex[5];m<=a.uindex[5];m++)
	      if(a[i][j][k][l][m]>c) {
		c=a[i][j][k][l][m];
		ii=i;jj=j;kk=k; ll=l;mm=m;
	      }
    return c;
  }
}

// --------------------------------------------------------------------------
namespace aurostd {
  template<class utype>                            // function max xtensor6<>
  utype max(const xtensor6<utype>& a,int& ii,int& jj,int& kk,int& ll,int& mm,int& nn) {
    aurostd::xtensor6debug(a,"function max");
    utype c=a[a.lindex[1]][a.lindex[2]][a.lindex[3]][a.lindex[4]][a.lindex[5]][a.lindex[6]];
    ii=a.lindex[1],jj=a.lindex[2],kk=a.lindex[3],ll=a.lindex[4],mm=a.lindex[5],nn=a.lindex[6];
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
	for(int k=a.lindex[3];k<=a.uindex[3];k++)
	  for(int l=a.lindex[4];l<=a.uindex[4];l++)
	    for(int m=a.lindex[5];m<=a.uindex[5];m++)
	      for(int n=a.lindex[6];n<=a.uindex[6];n++)
		if(a[i][j][k][l][m][n]>c) {
		  c=a[i][j][k][l][m][n];
		  ii=i;jj=j;kk=k;ll=l;mm=m;nn=n;
		}
    return c;
  }
}

// --------------------------------------------------------------------------
namespace aurostd {
  template<class utype>                            // function max xtensor7<>
  utype max(const xtensor7<utype>& a,int& ii,int& jj,int& kk,int& ll,int& mm,int& nn,int& oo) {
    aurostd::xtensor7debug(a,"function max");
    utype c=a[a.lindex[1]][a.lindex[2]][a.lindex[3]][a.lindex[4]][a.lindex[5]][a.lindex[6]][a.lindex[7]];
    ii=a.lindex[1],jj=a.lindex[2],kk=a.lindex[3],ll=a.lindex[4],mm=a.lindex[5],nn=a.lindex[6],oo=a.lindex[7];
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
	for(int k=a.lindex[3];k<=a.uindex[3];k++)
	  for(int l=a.lindex[4];l<=a.uindex[4];l++)
	    for(int m=a.lindex[5];m<=a.uindex[5];m++)
	      for(int n=a.lindex[6];n<=a.uindex[6];n++)
		for(int o=a.lindex[7];o<=a.uindex[7];o++)
		  if(a[i][j][k][l][m][n][o]>c) {
		    c=a[i][j][k][l][m][n][o];
		    ii=i;jj=j;kk=k;ll=l;mm=m;nn=n;oo=o;
		  }
    return c;
  }
}

// --------------------------------------------------------------------------
namespace aurostd {
  template<class utype>                            // function max xtensor8<>
  utype max(const xtensor8<utype>& a,int& ii,int& jj,int& kk,int& ll,int& mm,int& nn,int& oo,int& pp) {
    aurostd::xtensor8debug(a,"function max");
    utype c=a[a.lindex[1]][a.lindex[2]][a.lindex[3]][a.lindex[4]][a.lindex[5]][a.lindex[6]][a.lindex[7]][a.lindex[8]];
    ii=a.lindex[1],jj=a.lindex[2],kk=a.lindex[3],ll=a.lindex[4],mm=a.lindex[5],nn=a.lindex[6],oo=a.lindex[7],pp=a.lindex[8];
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
	for(int k=a.lindex[3];k<=a.uindex[3];k++)
	  for(int l=a.lindex[4];l<=a.uindex[4];l++)
	    for(int m=a.lindex[5];m<=a.uindex[5];m++)
	      for(int n=a.lindex[6];n<=a.uindex[6];n++)
		for(int o=a.lindex[7];o<=a.uindex[7];o++)
		  for(int p=a.lindex[8];p<=a.uindex[8];p++)
		    if(a[i][j][k][l][m][n][o][p]>c) {
		      c=a[i][j][k][l][m][n][o][p];
		      ii=i;jj=j;kk=k;ll=l;mm=m;nn=n;oo=o;pp=p;
		    }
    return c;
  }
}

// **************************************************************************
// --------------------------------------------------------------------------
namespace aurostd {
  template<class utype> utype                              // trace xtensor3
  trace(const xtensor3<utype>& a) {
    aurostd::xtensor3debug(a,"function trace");
    if(!a.iscubic) {cerr << "failure in xtensor3 trace defined for cubic xtensors" << endl;exit(0);}
    utype out=0;
    for(int i=0;i<a.index[1];i++)
      out+=a[i+a.lindex[1]][i+a.lindex[2]][i+a.lindex[3]];
    return out;
  }
}

// --------------------------------------------------------------------------
namespace aurostd {
  template<class utype> utype                              // trace xtensor4
  trace(const xtensor4<utype>& a) {
    aurostd::xtensor4debug(a,"function trace");
    if(!a.iscubic) {cerr << "failure in xtensor4 trace defined for cubic xtensors" << endl;exit(0);}
    utype out=0;
    for(int i=0;i<a.index[1];i++)
      out+=a[i+a.lindex[1]][i+a.lindex[2]][i+a.lindex[3]][i+a.lindex[4]];
    return out;
  }
}

// --------------------------------------------------------------------------
namespace aurostd {
  template<class utype> utype                              // trace xtensor5
  trace(const xtensor5<utype>& a) {
    aurostd::xtensor5debug(a,"function trace");
    if(!a.iscubic) {cerr << "failure in xtensor5 trace defined for cubic xtensors" << endl;exit(0);}
    utype out=0;
    for(int i=0;i<a.index[1];i++)
      out+=a[i+a.lindex[1]][i+a.lindex[2]][i+a.lindex[3]][i+a.lindex[4]][i+a.lindex[5]];
    return out;
  }
}

// --------------------------------------------------------------------------
namespace aurostd {
  template<class utype> utype                              // trace xtensor6
  trace(const xtensor6<utype>& a) {
    aurostd::xtensor6debug(a,"function trace");
    if(!a.iscubic) {cerr << "failure in xtensor6 trace defined for cubic xtensors" << endl;exit(0);}
    utype out=0;
    for(int i=0;i<a.index[1];i++)
      out+=a[i+a.lindex[1]][i+a.lindex[2]][i+a.lindex[3]][i+a.lindex[4]][i+a.lindex[5]][i+a.lindex[6]];
    return out;
  }
}

// --------------------------------------------------------------------------
namespace aurostd {
  template<class utype> utype                              // trace xtensor7
  trace(const xtensor7<utype>& a) {
    aurostd::xtensor7debug(a,"function trace");
    if(!a.iscubic) {cerr << "failure in xtensor7 trace defined for cubic xtensors" << endl;exit(0);}
    utype out=0;
    for(int i=0;i<a.index[1];i++)
      out+=a[i+a.lindex[1]][i+a.lindex[2]][i+a.lindex[3]][i+a.lindex[4]][i+a.lindex[5]][i+a.lindex[6]][i+a.lindex[7]];
    return out;
  }
}

// --------------------------------------------------------------------------
namespace aurostd {
  template<class utype> utype                              // trace xtensor8
  trace(const xtensor8<utype>& a) {
    aurostd::xtensor8debug(a,"function trace");
    if(!a.iscubic) {cerr << "failure in xtensor8 trace defined for cubic xtensors" << endl;exit(0);}
    utype out=0;
    for(int i=0;i<a.index[1];i++)
      out+=a[i+a.lindex[1]][i+a.lindex[2]][i+a.lindex[3]][i+a.lindex[4]][i+a.lindex[5]][i+a.lindex[6]][i+a.lindex[7]][i+a.lindex[8]];
    return out;
  }
}

// **************************************************************************
// --------------------------------------------------------------------------
namespace aurostd {
  template<class utype> xtensor3<utype>                  // identity xtensor3
  identity(const xtensor3<utype>& a) {
    aurostd::xtensor3debug(a,"function identity xtensor3");
    if(!a.iscubic) {cerr << "failure in xtensor3 identity defined for cubic xtensors" << endl;exit(0);}
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
	for(int k=a.lindex[3];k<=a.uindex[3];k++)
	  a[i][j][k]=(utype) 0;
    for(int i=0;i<a.index[1];i++)
      a[i+a.lindex[1]][i+a.lindex[2]][i+a.lindex[3]]=(utype) 1;
    return a;
  }
}

// --------------------------------------------------------------------------
namespace aurostd {
  template<class utype> xtensor4<utype>                  // identity xtensor4
  identity(const xtensor4<utype>& a) {
    aurostd::xtensor4debug(a,"function identity xtensor4");
    if(!a.iscubic) {cerr << "failure in xtensor4 identity defined for cubic xtensors" << endl;exit(0);}
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
	for(int k=a.lindex[3];k<=a.uindex[3];k++)
	  for(int l=a.lindex[4];l<=a.uindex[4];l++)
	    a[i][j][k][l]=(utype) 0;
    for(int i=0;i<a.index[1];i++)
      a[i+a.lindex[1]][i+a.lindex[2]][i+a.lindex[3]][i+a.lindex[4]]=(utype) 1;
    return a;
  }
}

// --------------------------------------------------------------------------
namespace aurostd {
  template<class utype> xtensor5<utype>                  // identity xtensor5
  identity(const xtensor5<utype>& a) {
    aurostd::xtensor5debug(a,"function identity xtensor5");
    if(!a.iscubic) {cerr << "failure in xtensor5 identity defined for cubic xtensors" << endl;exit(0);}
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
	for(int k=a.lindex[3];k<=a.uindex[3];k++)
	  for(int l=a.lindex[4];l<=a.uindex[4];l++)
	    for(int m=a.lindex[5];m<=a.uindex[5];m++)
	      a[i][j][k][l][m]=(utype) 0;
    for(int i=0;i<a.index[1];i++)
      a[i+a.lindex[1]][i+a.lindex[2]][i+a.lindex[3]][i+a.lindex[4]][i+a.lindex[5]]=(utype) 1;
    return a;
  }
}

// **************************************************************************
// --------------------------------------------------------------------------
namespace aurostd {
  template<class utype> xtensor6<utype>                  // identity xtensor6
  identity(const xtensor6<utype>& a) {
    aurostd::xtensor6debug(a,"function identity xtensor6");
    if(!a.iscubic) {cerr << "failure in xtensor6 identity defined for cubic xtensors" << endl;exit(0);}
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
	for(int k=a.lindex[3];k<=a.uindex[3];k++)
	  for(int l=a.lindex[4];l<=a.uindex[4];l++)
	    for(int m=a.lindex[5];m<=a.uindex[5];m++)
	      for(int n=a.lindex[6];n<=a.uindex[6];n++)
		a[i][j][k][l][m][n]=(utype) 0;
    for(int i=0;i<a.index[1];i++)
      a[i+a.lindex[1]][i+a.lindex[2]][i+a.lindex[3]][i+a.lindex[4]][i+a.lindex[5]][i+a.lindex[6]]=(utype) 1;
    return a;
  }
}

// --------------------------------------------------------------------------
namespace aurostd {
  template<class utype> xtensor7<utype>                  // identity xtensor7
  identity(const xtensor7<utype>& a) {
    aurostd::xtensor7debug(a,"function identity xtensor7");
    if(!a.iscubic) {cerr << "failure in xtensor7 identity defined for cubic xtensors" << endl;exit(0);}
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
	for(int k=a.lindex[3];k<=a.uindex[3];k++)
	  for(int l=a.lindex[4];l<=a.uindex[4];l++)
	    for(int m=a.lindex[5];m<=a.uindex[5];m++)
	      for(int n=a.lindex[6];n<=a.uindex[6];n++)
		for(int o=a.lindex[7];o<=a.uindex[7];o++)
		  a[i][j][k][l][m][n][o]=(utype) 0;
    for(int i=0;i<a.index[1];i++)
      a[i+a.lindex[1]][i+a.lindex[2]][i+a.lindex[3]][i+a.lindex[4]][i+a.lindex[5]][i+a.lindex[6]][i+a.lindex[7]]=(utype) 1;
    return a;
  }
}

// --------------------------------------------------------------------------
namespace aurostd {
  template<class utype> xtensor8<utype>                  // identity xtensor8
  identity(const xtensor8<utype>& a) {
    aurostd::xtensor8debug(a,"function identity xtensor8");
    if(!a.iscubic) {cerr << "failure in xtensor8 identity defined for cubic xtensors" << endl;exit(0);}
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
	for(int k=a.lindex[3];k<=a.uindex[3];k++)
	  for(int l=a.lindex[4];l<=a.uindex[4];l++)
	    for(int m=a.lindex[5];m<=a.uindex[5];m++)
	      for(int n=a.lindex[6];n<=a.uindex[6];n++)
		for(int o=a.lindex[7];o<=a.uindex[7];o++)
		  for(int p=a.lindex[8];p<=a.uindex[8];p++)
		    a[i][j][k][l][m][n][o][p]=(utype) 0;
    for(int i=0;i<a.index[1];i++)
      a[i+a.lindex[1]][i+a.lindex[2]][i+a.lindex[3]][i+a.lindex[4]][i+a.lindex[5]][i+a.lindex[6]][i+a.lindex[7]][i+a.lindex[8]]=(utype) 1;
    return a;
  }
}

#endif

// **************************************************************************
// *                                                                        *
// *             STEFANO CURTAROLO - Duke University 2003-2015              *
// *                                                                        *
// **************************************************************************

