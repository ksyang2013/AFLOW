// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *                                                                         *
// ***************************************************************************
// Stefano Curtarolo

#ifndef _AUROSTD_BOOT_CPP_
#define _AUROSTD_BOOT_CPP_
#include "aurostd.h"

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// **************************************************************************
// dont touch the functions below... their purpose is just to be compiled so
// the templates are constructed. SC

using aurostd::nint;
using aurostd::fact;
using aurostd::factorial;
using aurostd::sign;
using aurostd::_isodd;
using aurostd::_iseven;
using aurostd::_isfloat;
using aurostd::_iscomplex;
using aurostd::_isinteger;
using aurostd::_roundoff;
using aurostd::mod;
using aurostd::min;
using aurostd::max;
using aurostd::xcomplex;
using aurostd::xmatrix;
using aurostd::modulus;
using aurostd::angle;
using aurostd::xvector;
using aurostd::identity;
using aurostd::xmatrix2vectorvector;
using aurostd::vectorvector2xmatrix;
using aurostd::xvector2vector;
using aurostd::vector2xvector;
using aurostd::sign;
using aurostd::set;
using aurostd::reset;
using aurostd::clear;
using aurostd::reshape;
using aurostd::isequal;
using aurostd::isdifferent;
using aurostd::reduceByGCD;
using aurostd::normalizeSumToOne; //CO 180817
//using aurostd::swap;
using aurostd::swap_cols;
using aurostd::swap_columns;
using aurostd::swap_rows;
using aurostd::polar;
using aurostd::_real;
using aurostd::jacobi;
using aurostd::eigsrt;
using aurostd::tred2;
using aurostd::generalHouseHolderQRDecomposition;
using aurostd::tqli;
using aurostd::balanc;
using aurostd::elmhes;
using aurostd::hqr;
using aurostd::eigen;
using aurostd::trasp;
using aurostd::uniform;
using aurostd::gaussian;
using aurostd::expdev;
using aurostd::laplacedev;
using aurostd::utype2string;
using aurostd::GaussJordan;
using aurostd::gaussj;
using aurostd::xoption;
using aurostd::floor;
using aurostd::ceil;

template<class utype> bool initialize_scalar(utype d) {
  string s="";
  utype u=0;
  u+=aurostd::string2utype<utype>(aurostd::utype2string<utype>(utype())+aurostd::utype2string<utype>(utype(),int()));
  u+=aurostd::substring2utype<utype>(s,s,1)+aurostd::substring2utype<utype>(s,s,s,1);
  u+=aurostd::substring2utype<utype>(s,s)+aurostd::substring2utype<utype>(s,s,s);
  double o=0;
  o+=_isfloat(d)+_iscomplex(d);//abs(d);
  o+=max(d,d);o+=max(d,d,d);o+=max(d,d,d,d);o+=max(d,d,d,d,d);o+=min(d,d);o+=min(d,d,d);o+=min(d,d,d,d);o+=min(d,d,d,d,d);o+=_real(d);
  //  char* c=0;
  vector<string> vs;vector<utype> vu;
  deque<string> ds;deque<utype> du;
  o+=aurostd::args2utype(vs,"|",d);
  o+=aurostd::args2flag(vs,"");
  o+=aurostd::args2flag(vs,vs,"");
  o+=aurostd::get_itemized_vector_string_from_input(vs,"",vs,":");
  o+=aurostd::get_itemized_vector_string_from_input(vs,"","",vs,":");
  o+=aurostd::get_itemized_vector_string_from_input(vs,"","","",vs,":");
  o+=aurostd::get_itemized_vector_string_from_input(vs,"","","","",vs,":");
  o+=aurostd::args2attachedutype<utype>(vs,"",0);
  // [OBSOLETE] o+=aurostd::args2attachedutype<utype>(vs,"","",0);
  // [OBSOLETE] o+=aurostd::args2attachedutype<utype>(vs,"","","",0);
  // [OBSOLETE] o+=aurostd::args2attachedutype<utype>(vs,"","","","",0);

  
  //  utype uu; _aflowlib::web2utype("location",uu);
  // [OBSOLETE]  _aflowlib::web2utype("location",s);
  aurostd::url2tokens("location",vu); aurostd::url2tokens("location",vu," ");
  aurostd::url2tokens("location",vs); aurostd::url2tokens("location",vs," "); 
  aurostd::url2tokens("location",du); aurostd::url2tokens("location",du," ");
  aurostd::url2tokens("location",ds); aurostd::url2tokens("location",ds," "); 

  u+=aurostd::execute2utype<int>("");u+=aurostd::execute2utype<int>(string(""));
  cout << aurostd::PaddedPRE(d,1) << aurostd::PaddedPRE("",1) << aurostd::PaddedPRE(s,1);
  cout << aurostd::PaddedPOST(d,1) << aurostd::PaddedPOST("",1) << aurostd::PaddedPOST(s,1);

  xoption opt;
  opt.getattachedutype<utype>("");
  opt.args2addattachedscheme<utype>(vs,vs,"","",u);
  opt.args2addattachedscheme<utype>(vs,"","",u);
  return (o<0);
}

template<class utype> bool initialize_xcomplex(utype d) {
  xcomplex<utype> x(0,0);utype r=d;
  x+x;x+r;r+x;x+=x;x+=r; // plus
  x-x;x-r;r-x;x-=x;x-=r; // minus
  x*x;x*r;r*x;x*=x;x*=r; // multiplication
  x/x;x/r;r/x;x/=x;x/=r; // division
  x=x;x=r;               // equal
  x=abs(x)+arg(x)+polar(r,r)+conj(x)+norm(x)+cos(x)+cosh(x)+exp(x)+log(x)+pow(x,x)+pow(x,int(1))+pow(x,r)+pow(r,x)+sin(x)+sinh(x)+sqrt(x); // functions
  r = magsqr(x);  // ME180907

  aurostd::xvector<xcomplex<utype> > vxfirst,vx(2),vbb=vxfirst,vcc(vxfirst),vxvx,vxvxvxvx(1,2),vb=vxfirst,vc(vxfirst); // DX 1/15/18 - equal operator was missing
  aurostd::xvector<utype> vr;
  vx=vx;vx+vx;vx+=vx;vx-vx;vx-=vx;vx*vx; // DX 1/17/18 - added vx=vx
  //xvector<xcomplex<utype> > va,vb=va,vc(va); // DX 1/15/18 - equal operator was missing
  sin(vx);sinh(vx);cos(vx);cosh(vx);exp(vx);log(vx);//sqrt(vx);
  vx.clear(); // DX 1/15/18 - clear was missing
  cout << vx << endl; // DX 1/15/18 - ostream was missing
  conj(vx);  // ME180904

  aurostd::xmatrix<utype > mx(2),mxmx,mxmxmx(2,2),mxmxmxmxmx(1,2,3,4);
  mx+mx;mx+=mx;mx-mx;mx-=mx;mx*mx;vx(1)=vx(1);vx[1]=vx[1];
  sin(mx);sinh(mx);cos(mx);cosh(mx);exp(mx);
  aurostd::ones<utype>();aurostd::ones<utype>(3);aurostd::ones<utype>(3,3);aurostd::ones<utype>(1,2,3,4);
  aurostd::eyes<utype>();aurostd::eyes<utype>(3);aurostd::eyes<utype>(3,3);aurostd::eyes<utype>(1,2,3,4);
  aurostd::CMdet<utype>(mx);  //CO 180515
  cout << mx << endl;  // DX 1/15/18 - ostream was missing

  aurostd::xmatrix<xcomplex<utype> > a,b(1,1),c(1,2,3,4),m(2),e=a,f(a);
  x+=issymmetric(m)+isantisymmetric(m)+ishermitian(m)+isantihermitian(m);
  trace(m); // DX 1/15/17 - initialize trace for xcomplex
  m(1,1)=m(1,1);m[1][1]=m[1][1];m=m;m=m+m;m=m*m;m=m-m;m.clear();//m=m*r;m=r*m;
  m(1,1)+=m(1,1);m[1][1]-=m[1][1];m[1][1]*=m[1][1];m[1][1]/=m[1][1]; // DX 1/15/18 - operator and equal operator initialzied for xcomplex
  exp(m); // DX 1/15/18 - add exponential or complex matrices
  m=x*m;m=m/x; // DX 1/17/18 - allow for xcomplex * xmatrix<xcomplex>
  cout << m << endl; // DX 1/15/18 - ostream
  vx=m.getcol(1);m=conj(m);trasp(m);trasp(vx);vx=m*vx;  // ME 180904

  //  jacobi(m,vx,m);
 
  return TRUE;
}

template<class utype> bool initialize_eigenproblems(utype x) {
  aurostd::xvector<utype> v(3);v[1]=x;
  aurostd::xmatrix<utype> m(3,3);m[1][1]=x;
  int i=0;i+=jacobi(m,v,m);
  eigsrt(v,m);generalHouseHolderQRDecomposition(m);tred2(m,v,v);tqli(v,v,m);balanc(m);elmhes(m);hqr(m,v,v);eigen(m,v,v); // EIGENVECTORS
  return TRUE;
}

template<class utype> bool initialize_xscalar_xvector_xmatrix_xtensor(utype x) {
  char* c=0;
  int i=0;double d=0,o=0;string s;
  aurostd::xvector<int> xvi(2);
  aurostd::xvector<double> xvd(2);
  
  // initialize vector sort
  vector<utype> vutype;vector<int> vint;vector<uint> vuint;vector<double> vdouble;vector<string> vstring;
  aurostd::sort(vutype);aurostd::sort(vutype,vint);aurostd::sort(vutype,vdouble);
  aurostd::sort(vutype,vutype,vint);aurostd::sort(vutype,vutype,vutype);aurostd::sort(vutype,vint,vint);
  aurostd::sort(vutype,vutype,vdouble);aurostd::sort(vutype,vdouble,vutype);aurostd::sort(vutype,vdouble,vdouble);
  aurostd::sort(vutype,vdouble,vint);aurostd::sort(vutype,vint,vdouble);aurostd::sort(vutype,vint,vint);
  aurostd::sort(vstring);aurostd::sort(vstring,vint);aurostd::sort(vstring,vdouble);aurostd::sort(vstring,vstring);
  aurostd::sort(vdouble);aurostd::sort(vdouble,vint);aurostd::sort(vdouble,vdouble);aurostd::sort(vdouble,vstring);
  aurostd::sort(vstring,vint,vstring);aurostd::sort(vstring,vdouble,vstring);aurostd::sort(vstring,vstring,vstring);
  aurostd::sort(vstring,vstring,vdouble,vstring);aurostd::sort(vstring,vstring,vdouble,vdouble,vstring);
  aurostd::sort_remove_duplicates(vstring);
  // cout << vutype << vint << vdouble << vstring;
  deque<utype> dutype;deque<int> dint;deque<uint> duint;deque<double> ddouble;deque<string> dstring;
  // aurostd::sort(dstring);aurostd::sort(dstring,dint);aurostd::sort(dstring,ddouble);aurostd::sort(dstring,dstring);
  // aurostd::sort(ddouble);aurostd::sort(ddouble,dint);aurostd::sort(ddouble,ddouble);aurostd::sort(ddouble,dstring);
  aurostd::sort(dstring,dint,dstring);aurostd::sort(dstring,ddouble,dstring);aurostd::sort(dstring,dstring,dstring);
  aurostd::sort(dstring,dstring,ddouble,dstring);aurostd::sort(dstring,dstring,ddouble,ddouble,dstring);

  // initialize vector/deque 
  aurostd::string2tokens(s,vstring,"");aurostd::string2tokens(s,dstring,"");
  aurostd::string2tokensAdd(s,vstring,"");aurostd::string2tokensAdd(s,dstring,"");
  aurostd::string2tokens(s,vuint,"");aurostd::string2tokens(s,duint,"");
  aurostd::string2tokensAdd(s,vuint,"");aurostd::string2tokensAdd(s,duint,"");
  aurostd::string2tokens(s,vutype,"");aurostd::string2tokens(s,dutype,"");
  aurostd::string2tokensAdd(s,vutype,"");aurostd::string2tokensAdd(s,dutype,"");
  
  // initialize xvector sort
  xvector<utype> vxu;xvector<int> ixv;xvector<double> dxv;
  aurostd::quicksort2((uint) 0,ixv,dxv);aurostd::quicksort2((uint) 0,dxv,ixv);
  
  // initialize scalars
  o+=nint(x)+factorial(x)+fact(x)+isequal(x,x)+isequal(x,x,(utype) 0)+isequal(x,x,x)+isdifferent(x,x)+isdifferent(x,x,x);
  o+=max(x,x);o+=max(x,x,x);o+=max(x,x,x,x);o+=max(x,x,x,x,x);o+=max(x,x,x,x,x,x);o+=sign(x);
  o+=min(x,x);o+=min(x,x,x);o+=min(x,x,x,x);o+=min(x,x,x,x,x);o+=min(x,x,x,x,x,x);o+=sign(x);
  o+=_isinteger(x,x)+uniform(x)+uniform(x,x)+gaussian(x)+gaussian(x,x)+expdev(x)+laplacedev(x)+laplacedev(x,x);
  o+=_roundoff(x,x);
  o+=combinations(x,x)+Cnk(x,x);

  // initialize xvectors
  aurostd::xvector<utype> v(2),vv,vvvv(1,2);
  std::vector<aurostd::xvector<utype> > vxv;
  v(1);v[1];o+=sum(v);sort(shellsort(heapsort(v+v)));o+=modulus(v);o+=modulus2(v);v=x*v*x/x;v=+v;v=-v;v+=v;v-=v;v+=x;v-=x;v*=x;v/=x;
  v=v*i;v=nint(v);v=sign(v);v=i*v;v=v*d;v=d*v;v=d*v*d/d;o+=min(v);o+=max(v);o+=mini(v);o+=maxi(v);v=v+v;v=v-v;v=-v+v;
  cout<<v<<endl;o+=scalar_product(v,v);o+=cos(v,v);o+=sin(v,v);o+=angle(v,v);v=vector_product(v,v);
  trasp(v);
  //  sin(v);sinh(v);cos(v);cosh(v);exp(v);

  o=+(v==v);o=+(v!=v);o+=identical(v,v);o+=identical(v,v,x);o+=identical(v,v,(utype&) x);o+=identical(v,v,(const utype&) x);o+=isdifferent(v,v);
  o+=isdifferent(v,v,x);v=-v;o+=max(v);v=abs(v);roundoff(v);roundoff(v,x);v+=reduceByGCD(v,x);v+=normalizeSumToOne(v,x);clear(v);floor(v);ceil(v);
  o+=getcos(v,v);v.clear();v.set(x);v.reset();clear(v);reset(v);set(v,x);v=abs(v);v=vabs(v);v=sign(v);
  o+=angle(v,v,v);o+=getangle(v,v,v);isCollinear(v,v,x);v=getCentroid(vxv);o+=distance(v,v);
  getGeneralAngles(v,x);getGeneralAngle(v,0,x);
  getGeneralNormal(vxv);
  v=aurostd::args2xvectorutype<utype>(vstring,c,v);v=aurostd::args2xvectorutype<utype>(vstring,c,i);
  aurostd::args2vectorutype<utype>(vstring,c);aurostd::args2dequeutype<utype>(dstring,c);
  v.clear();clear(v);reshape(x);reshape(x,x);reshape(x,x,x);reshape(x,x,x,x);reshape(x,x,x,x,x);reshape(x,x,x,x,x,x);
  //  sort2((unsigned long)1,v,v);sort2((unsigned long)1,xvd,xvi);sort2((unsigned long)1,xvi,xvd);
  o+=isequal(v,v)+isequal(v,v,(utype) 0)+isequal(v,v,x)+isdifferent(v,v)+isdifferent(v,v,x)+isinteger(v,x);swap(v,1,1);shiftlrows(v,1);
  getQuartiles(v,x,x,x);
  o+=getMAD(v,x);

  // initialize matrices
  utype* mstar;mstar=NULL;
  aurostd::xmatrix<utype> m(2),mm,mmm(2,3),mmmmm(1,2,3,4),m5(1,2,mstar),mkron;
  xdouble(m);xint(m);m=+m;m=-m;o+=m(1)[1];o+=m(1,1);o+=m[1][1];m=identity(m);m=identity(x,1,1);
  vv=m.getcol(1);m=m;m=m+m;m=m-m;m=m*m;m=inverse(m);inverse(m,m);m=reduce_to_shortest_basis(m);
  m=x*m*x/x;o+=(m==m);o+=(m!=m);o+=trace(m);m=-m;m=trasp(m);clear(m);mkron=aurostd::KroneckerProduct(mm,mmm);
  o+=sum(m);m=nint(m);m=sign(m);o+=identical(m,m);o+=identical(m,m,x);o+=isdifferent(m,m);o+=isdifferent(m,m,x);
  o+=isequal(m,m);o+=isequal(m,m,x);cout<<m<<endl;
  roundoff(m);roundoff(m,x);m=exp(m);aurostd::abs(m);o+=max(m);xdouble(v);xint(v);
  m.clear();m.set(x);m.reset();m.reset();clear(m);reset(m);clear(m);set(m,x);abs(m);mabs(m);
  o+=det(m);o+=determinant(m);m=submatrix(m,1,1);o+=minordet(m,1,1);o+=minordeterminant(m,1,1);v*m;m*v;
  vector<utype> vvv=xvector2vector(v);vector2xvector(vvv,1);
  vector<vector<utype> > mvv=xmatrix2vectorvector(m);vectorvector2xmatrix(mvv);
  if(det(m)==0 || sum(m)==0) return FALSE;
  reshape(v);reshape(v,v);reshape(v,v,v);reshape(v,v,v,v);reshape(v,v,v,v,v);reshape(v,v,v,v,v,v);
  reshape_rows(v);reshape_rows(v,v);reshape_rows(v,v,v);reshape_rows(v,v,v,v);reshape_rows(v,v,v,v,v);reshape_rows(v,v,v,v,v,v);
  reshape_cols(v);reshape_cols(v,v);reshape_cols(v,v,v);reshape_cols(v,v,v,v);reshape_cols(v,v,v,v,v);reshape_cols(v,v,v,v,v,v);
  utype tol; // DX 171025
  o+=isequal(m,m)+isequal(m,m,(utype) 0)+isequal(m,m,x)+isdifferent(m,m)+isdifferent(m,m,x)+isinteger(m,x)+isdiagonal(m)+isdiagonal(m,tol)+issymmetric(m)+isantisymmetric(m);
  o+=isidentity(m); //CO
  swap_cols(m,1,1);swap_columns(m,1,1);swap_rows(m,1,1);
  sin(m);sinh(m);cos(m);cosh(m);exp(m);
  aurostd::floor(m);aurostd::ceil(m);
  // aurostd::trunc(m);aurostd::round(m);

  //[ME 180627 START]
  std::vector<int> stdv(3, 3),stdv2(2, 3),vind(3, 1),vind2(2),stdv0(3, 0); std::vector<utype> vut(3);
  aurostd::xvector<int> xv = aurostd::vector2xvector(stdv), xv0(3), xvind2 = aurostd::vector2xvector(vind2);
  aurostd::xvector<utype> xvut(3);aurostd::xvector<int> xvind = aurostd::vector2xvector(vind);aurostd::xmatrix<utype> xmut(3, 3);
  aurostd::xtensor<utype> tdef,t(3, 0),t1(3),t2(stdv2),tv(stdv,stdv0),tv1(stdv),tx(xv,xv0),tx1(xv);
  x=t(vind);aurostd::xtensor<utype> ts(t[1]);ts=t[1];ts=-t[1];ts=+t[1];
  x=t(xvind);x=t[1][1][1];x=t[0](vind2);x=t[0](xvind2);t*=x;t/=x;t+=x;t-=x;t=t/x;t=x*t;t=t*x;t=t;
  t=-t;t=+t;t=t+t;t=t-t;tdef=t;t+=t;t-=t;t[0][0]+=x;t[0][0]-=x;t[0][0]*=x;t[0][0]/=1;t[0]+=t[0];t[0]-=t[0];
  t2+=t[0];t2-=t[0];t[0]+=t2;t[0]-=t2;t2=t2+t[0];t2=t[0]+t2;t2=t[0]+t[0];t2=t2-t[0];t2=t[0]-t2;t2=t[0]-t[0];
  t[0][0]=xvut;t[0][0]=vut;t[0]=xmut;t[0][0][0]=x;t[0]=t[0];t[0]=t2;x=t.get(0);t.set(x, 0);t[0].set(x);
  t(vind)=x;t(xvind)=x;t(vind)=t(vind);t(xvind)=t(xvind);xvut=aurostd::xtensor2xvector(t[0][0]);
  vut=aurostd::xtensor2vector(t[0][0]);xmut=aurostd::xtensor2xmatrix(t[9]);t=aurostd::vector2xtensor(vut);
  t=aurostd::xvector2xtensor(xvut);t=aurostd::xmatrix2xtensor(xmut);t.clear();t.reset();t.set(x);t=aurostd::nint(t);
  t=aurostd::sign(t);t=aurostd::abs(t);t=aurostd::floor(t);t=aurostd::ceil(t);t=aurostd::round(t);x=aurostd::max(t);
  x=aurostd::min(t);x=aurostd::sum(t);x=aurostd::trace(t);t=aurostd::identity_tensor(x,3);t[0]=aurostd::nint(t[0]);
  t[0]=aurostd::sign(t[0]);t[0]=aurostd::abs(t[0]);t[0]=aurostd::floor(t[0]);t[0]=aurostd::ceil(t[0]);
  t[0]=aurostd::round(t[0]);x=aurostd::max(t[0]);x=aurostd::min(t[0]);x=aurostd::sum(t[0]);x=aurostd::trace(t[0]);
  //[ME 180627 END]

  //[OBSOLETE ME180705]// initialize tensor3
  //[OBSOLETE ME180705]aurostd::xtensor3<utype> t3(1),t3a(1),t3b(1),t3c(1,1),t3d(1,1,1);
  //[OBSOLETE ME180705]x=t3[3][3][3];x=t3(3,3,3);t3+=t3;t3-=t3;t3=t3;t3=t3+t3;t3=t3-t3;t3=t3/x;t3=t3*x;
  //[OBSOLETE ME180705]t3=x*t3;t3=-t3;t3=+t3;reset(t3);clear(t3);set(t3,x);x=min(t3);x=max(t3);x=trace(t3);t3=identity(t3);
  //[OBSOLETE ME180705]t3.clear();t3.set(x);t3.reset();t3=nint(t3);t3=sign(t3);t3=abs(t3);t3=floor(t3);t3=ceil(t3);//t3=round(t3);t3=trunc(t3);
  //[OBSOLETE ME180705]int i3[3+1];aurostd::xtensor3debug(t3,"");utype *** b3;aurostd::allocate_xtensor3corpus(b3,i3,i3,i3);
  //[OBSOLETE ME180705]// t3=uniform(t3,(utype) 1,(utype) 2);t3=uniform(t3,(utype) 2);t3=uniform(t3);
  //[OBSOLETE ME180705]// t3=uniform(t3,(utype) 1,(utype) 1,(utype) 1,(utype) 1);t3=uniform(t3,(utype) 1,(utype) 2);t3=uniform(t3);
  //[OBSOLETE ME180705]// t3=gaussian(t3,m,sigma);t3=gaussian(t3,sigma);t3=gaussian(t3);
  //[OBSOLETE ME180705]// t3=gaussian(t3,m1,sigma1,m2,sigma2);t3=gaussian(t3,sigma1,sigma2);
  //[OBSOLETE ME180705]// t3=gaussian(t3,sigma);t3=gaussian(t3);

  //[OBSOLETE ME180705]// initialize tensor4
  //[OBSOLETE ME180705]aurostd::xtensor4<utype> t4(1),t4a(1),t4b(1),t4c(1,1),t4d(1,1,1),t4e(1,1,1,1);
  //[OBSOLETE ME180705]x=t4[4][4][4][4];x=t4(4,4,4,4);t4+=t4;t4-=t4;t4=t4;t4=t4+t4;t4=t4-t4;t4=t4/x;t4=t4*x;
  //[OBSOLETE ME180705]t4=x*t4;t4=-t4;t4=+t4;reset(t4);clear(t4);set(t4,x);x=min(t4);x=max(t4);x=trace(t4);t4=identity(t4);
  //[OBSOLETE ME180705]t4.clear();t4.set(x);t4.reset();t4=nint(t4);t4=sign(t4);t4=abs(t4);t4=floor(t4);t4=ceil(t4);//t4=round(t4);t4=trunc(t4);
  //[OBSOLETE ME180705]int i4[4+1];aurostd::xtensor4debug(t4,"");utype **** b4;aurostd::allocate_xtensor4corpus(b4,i4,i4,i4);

  //[OBSOLETE ME180705]// initialize tensor5
  //[OBSOLETE ME180705]aurostd::xtensor5<utype> t5(1),t5a(1),t5b(1),t5c(1,1),t5d(1,1,1),t5e(1,1,1,1),t5f(1,1,1,1,1);
  //[OBSOLETE ME180705]x=t5[5][5][5][5][5];x=t5(5,5,5,5,5);t5+=t5;t5-=t5;t5=t5;t5=t5+t5;t5=t5-t5;t5=t5/x;t5=t5*x;
  //[OBSOLETE ME180705]t5=x*t5;t5=-t5;t5=+t5;reset(t5);clear(t5);set(t5,x);x=min(t5);x=max(t5);x=trace(t5);t5=identity(t5);
  //[OBSOLETE ME180705]t5.clear();t5.set(x);t5.reset();t5=nint(t5);t5=sign(t5);t5=abs(t5);t5=floor(t5);t5=ceil(t5);//t5=round(t5);t5=trunc(t5);
  //[OBSOLETE ME180705]int i5[5+1];aurostd::xtensor5debug(t5,"");utype ***** b5;aurostd::allocate_xtensor5corpus(b5,i5,i5,i5);

  //[OBSOLETE ME180705]// initialize tensor6
  //[OBSOLETE ME180705]aurostd::xtensor6<utype> t6(1),t6a(1),t6b(1),t6c(1,1),t6d(1,1,1),t6e(1,1,1,1),t6f(1,1,1,1,1),t6g(1,1,1,1,1,1);
  //[OBSOLETE ME180705]x=t6[6][6][6][6][6][6];x=t6(6,6,6,6,6,6);t6+=t6;t6-=t6;t6=t6;t6=t6+t6;t6=t6-t6;t6=t6/x;t6=t6*x;
  //[OBSOLETE ME180705]t6=x*t6;t6=-t6;t6=+t6;reset(t6);clear(t6);set(t6,x);x=min(t6);x=max(t6);x=trace(t6);t6=identity(t6);
  //[OBSOLETE ME180705]t6.clear();t6.set(x);t6.reset();t6=nint(t6);t6=sign(t6);t6=abs(t6);t6=floor(t6);t6=ceil(t6);//t6=round(t6);t6=trunc(t6);
  //[OBSOLETE ME180705]int i6[6+1];aurostd::xtensor6debug(t6,"");utype ****** b6;aurostd::allocate_xtensor6corpus(b6,i6,i6,i6);

  //[OBSOLETE ME180705]// initialize tensor7
  //[OBSOLETE ME180705]aurostd::xtensor7<utype> t7(1),t7a(1),t7b(1),t7c(1,1),t7d(1,1,1),t7e(1,1,1,1),t7f(1,1,1,1,1),t7g(1,1,1,1,1,1,1);
  //[OBSOLETE ME180705]x=t7[7][7][7][7][7][7][7];x=t7(7,7,7,7,7,7,7);t7+=t7;t7-=t7;t7=t7;t7=t7+t7;t7=t7-t7;t7=t7/x;t7=t7*x;
  //[OBSOLETE ME180705]t7=x*t7;t7=-t7;t7=+t7;reset(t7);clear(t7);set(t7,x);x=min(t7);x=max(t7);x=trace(t7);t7=identity(t7);
  //[OBSOLETE ME180705]t7.clear();t7.set(x);t7.reset();t7=nint(t7);t7=sign(t7);t7=abs(t7);t7=floor(t7);t7=ceil(t7);//t7=round(t7);t7=trunc(t7);
  //[OBSOLETE ME180705]int i7[7+1];aurostd::xtensor7debug(t7,"");utype ******* b7;aurostd::allocate_xtensor7corpus(b7,i7,i7,i7);

  //[OBSOLETE ME180705]// initialize tensor8
  //[OBSOLETE ME180705]aurostd::xtensor8<utype> t8(1),t8a(1),t8b(1),t8c(1,1),t8d(1,1,1),t8e(1,1,1,1),t8f(1,1,1,1,1),t8g(1,1,1,1,1,1,1),t8h(1,1,1,1,1,1,1,1);
  //[OBSOLETE ME180705]x=t8[8][8][8][8][8][8][8][8];x=t8(8,8,8,8,8,8,8,8);t8+=t8;t8-=t8;t8=t8;t8=t8+t8;t8=t8-t8;t8=t8/x;t8=t8*x;
  //[OBSOLETE ME180705]t8=x*t8;t8=-t8;t8=+t8;reset(t8);clear(t8);set(t8,x);x=min(t8);x=max(t8);x=trace(t8);t8=identity(t8);
  //[OBSOLETE ME180705]t8.clear();t8.set(x);t8.reset();t8=nint(t8);t8=sign(t8);t8=abs(t8);t8=floor(t8);t8=ceil(t8);//t8=round(t8);t8=trunc(t8);
  //[OBSOLETE ME180705]int i8[8+1];aurostd::xtensor8debug(t8,"");utype ******** b8;aurostd::allocate_xtensor8corpus(b8,i8,i8,i8);
  
  return TRUE;
}

template<class utype> void utype_funcs(utype x,aurostd::xvector<utype>& afunc) {afunc[0]+=x;};

bool initialize_templates_never_call_this_procedure(bool flag) {
  // xoption
  xoption opt1,opt2;
  opt1.clear();opt1.push(string("something"));opt1.pop(string("something"));
  opt2=opt1;

  double o=0;
  o=aurostd::sqrt((double) 0.0)+aurostd::sqrt((int) 0.0)+aurostd::sqrt((float) 0.0);
  // return TRUE;
  // NEVER CALL THIS FUNCTION ... JUST TO INITIALIZE TEMPLATES
  // this function is called with aflow.cpp and it is used to create the various templates used in the whole code.
  // do not call it but just leave it here. (stefano).
  // if(flag) return TRUE; else return TRUE;
  if(flag) {
    //char* c=0;vector<string> vs=0;
    int i=2;double d=2.0;float f=2.0;
    i=max((uint) i,(uint) i);i=min((uint) i,(uint) i);
    o+=modulus(d,d);o+=modulus(d,d,d);o+=angle(d,d,d,d);o+=angle(d,d,d,d,d,d);o+=sign(d);o+=sign(i);o+=sign(f);o+=sign(d);
    o+=_isodd((int) 1)+_iseven((int) 1)+_isodd((uint) 1)+_iseven((uint) 1);
    o+=_isfloat((bool) 1)+_isfloat((char) 1)+_isfloat((int) 1)+_isfloat((float) 1)+_isfloat((double) 1);
    o+=_iscomplex((bool) 1)+_iscomplex((char) 1)+_iscomplex((int) 1)+_iscomplex((float) 1)+_iscomplex((double) 1);
    o+=_isodd((int) i)+_iseven((int) i)+_isodd((uint) i)+_iseven((uint) i);
    o+=mod((int) 20, (int) 2);
    string s;vector<string> vs;
    utype2string("tostring");
    aurostd::string2utype<string>("tostring");

    o+=aurostd::string2utype<bool>("0");
    o+=aurostd::string2utype<char>("0");
    o+=aurostd::string2utype<unsigned int>("0");
    // o+=aurostd::string2utype<uint>("0");
    o+=aurostd::string2utype<int>("0");
    o+=aurostd::string2utype<long int>("0");
    o+=aurostd::string2utype<long long int>("0");
    o+=aurostd::string2utype<float>("0");
    o+=aurostd::string2utype<double>("0");
    o+=aurostd::string2utype<long double>("0");
    
#define AUROSTD_INITIALIZE_BOOL
    //#define AUROSTD_INITIALIZE_CHAR
#define AUROSTD_INITIALIZE_STRING
#define AUROSTD_INITIALIZE_INT
#define AUROSTD_INITIALIZE_UINT
#define AUROSTD_INITIALIZE_FLOAT
#define AUROSTD_INITIALIZE_DOUBLE
#define AUROSTD_INITIALIZE_LONG_DOUBLE
#define AUROSTD_INITIALIZE_LONG_INT
#define AUROSTD_INITIALIZE_UNSIGNED_LONG_INT
#define AUROSTD_INITIALIZE_LONG_LONG_INT
#define AUROSTD_INITIALIZE_UNSIGNED_LONG_LONG_INT
#define AUROSTD_INITIALIZE_COMPLEX_DOUBLE
      

#ifdef AUROSTD_INITIALIZE_BOOL
    o+=initialize_scalar(bool(FALSE));
#endif
#ifdef AUROSTD_INITIALIZE_STRING
    //     string s="";
    // s+=aurostd::args2attachedutype<string>(vs,string(""),string(""));
    // s+=aurostd::args2attachedutype<string>(vs,string(""),string(""),string(""));
    // s+=aurostd::args2attachedutype<string>(vs,string(""),string(""),string(""),string(""));
    // s+=aurostd::args2attachedutype<string>(vs,string(""),string(""),string(""),string(""),string(""));
    aurostd::xcombos xc(1,1); //CO 180627
    vs.clear();vs.push_back("e.g."); //CO 180627
    while(xc.increment()){xc.applyCombo(vs);} //CO 180627
    vector<vector<string> > vvs;
    vvs.push_back(vector<string>(0));vvs[0].push_back("A");
    vvs.push_back(vector<string>(0));vvs[0].push_back("B");vvs[0].push_back("C");
    vector<int> bits;bits.push_back(vvs[0].size());bits.push_back(vvs[1].size());
    xc.reset(bits,'E');
    while(xc.increment()){xc.applyCombo(vvs);}
#endif
#ifdef AUROSTD_INITIALIZE_CHAR
    o+=initialize_scalar(char(1));
    o+=initialize_xscalar_xvector_xmatrix_xtensor(char(1));
#endif
#ifdef AUROSTD_INITIALIZE_INT
    o+=initialize_scalar(int(1));
    o+=initialize_xscalar_xvector_xmatrix_xtensor(int(1));
    aurostd::ones<int>();aurostd::ones<int>(3);aurostd::ones<int>(3,3);aurostd::ones<int>(1,2,3,4); //CO 180515
    aurostd::eyes<int>();aurostd::eyes<int>(3);aurostd::eyes<int>(3,3);aurostd::eyes<int>(1,2,3,4); //CO 180515
    xmatrix<int> mx=aurostd::eyes<int>();aurostd::CMdet<int>(mx);  //CO 180515
#endif
#ifdef AUROSTD_INITIALIZE_UINT
    o+=aurostd::string2utype<uint>(aurostd::utype2string<uint>(uint())+aurostd::utype2string<uint>(uint(),int()));
    o+=initialize_scalar(uint(1));
    //  o+=initialize_xscalar_xvector_xmatrix_xtensor(uint(1));
#endif
#ifdef AUROSTD_INITIALIZE_FLOAT
    if(1) { // AUROSTD_INITIALIZE_FLOAT
      o+=initialize_scalar(float(1));
      o+=initialize_xscalar_xvector_xmatrix_xtensor(float(1));
      o+=initialize_xcomplex(float(1));
      o+=initialize_eigenproblems(float(1));
      xvector<float> x(1);xvector<int> ia;float chi=1,d=1;
      xmatrix<float> z(1,1);gaussj(z,z.rows,z,z.urows);GaussJordan(z,z);
      aurostd::lfit(x,x,x,x,ia,z,chi,utype_funcs);
      z=reduce_to_shortest_basis(z);
      o+=aurostd::isequal<float>(float(d),float(d));
      o+=aurostd::isequal<float>(float(d),float(d),float(d));
      o+=aurostd::isdifferent<float>(float(d),float(d));
      o+=aurostd::isdifferent<float>(float(d),float(d),float(d));
     }
#endif
#ifdef AUROSTD_INITIALIZE_DOUBLE
    if(1) { // AUROSTD_INITIALIZE_DOUBLE
      o+=initialize_scalar(double(1));
      o+=initialize_xscalar_xvector_xmatrix_xtensor(double(1));
      o+=initialize_xcomplex(double(1));
      o+=initialize_eigenproblems(double(1));
      xvector<double> x(1);xvector<int> ia;double chi=1,d=1;
      xmatrix<double> z(1,1);gaussj(z,z.rows,z,z.urows);GaussJordan(z,z);
      aurostd::lfit(x,x,x,x,ia,z,chi,utype_funcs);
      z=reduce_to_shortest_basis(z);
      o+=aurostd::isequal<double>(double(d),double(d));
      o+=aurostd::isequal<double>(double(d),double(d),double(d));
      o+=aurostd::isequal<double>(double(d),double(d),0.0);
      o+=aurostd::isequal<double>(double(d),0.0,0.0);
      o+=aurostd::isequal<double>(0.0,0.0,0.0);
      o+=aurostd::isdifferent<double>(double(d),double(d));
      o+=aurostd::isdifferent<double>(double(d),double(d),double(d));
      o+=aurostd::isdifferent<double>(double(d),double(d),0.0);
      const double dd=0;
      o+=aurostd::identical(dd,dd,dd);
    }
#endif   
#ifdef AUROSTD_INITIALIZE_LONG_DOUBLE
    if(1) { // AUROSTD_INITIALIZE_LONG_DOUBLE
      o+=initialize_scalar((long double)(1));
      o+=initialize_xscalar_xvector_xmatrix_xtensor((long double)(1));
      xvector<long double> x(1);xvector<int> ia;long double chi=1;
      xmatrix<long double> z(1,1);gaussj(z,z.rows,z,z.urows);GaussJordan(z,z);
      z=reduce_to_shortest_basis(z);
      aurostd::lfit(x,x,x,x,ia,z,chi,utype_funcs);
      // o+=initialize_xcomplex((long double)(1));
      // xmatrix<(long double)> mm(1,1);GaussJordan(m,m);
    }
#endif
#ifdef AUROSTD_INITIALIZE_LONG_INT
    o+=initialize_scalar((long int)(1));
    o+=initialize_xscalar_xvector_xmatrix_xtensor((long int)(1));
    // o+=initialize_xcomplex((long int)(1));
    // xmatrix<(long int)> m(1,1);GaussJordan(m,m);
#endif
#ifdef AUROSTD_INITIALIZE_UNSIGNED_LONG_INT
    o+=aurostd::string2utype<unsigned long int>(aurostd::utype2string<unsigned long int>((unsigned long int)(1))+aurostd::utype2string<unsigned long int>((unsigned long int)(1),int()));
    // o+=initialize_scalar(((unsigned long int))(1));
    // o+=initialize_xscalar_xvector_xmatrix_xtensor(((unsigned long int))(1));
    // o+=initialize_xcomplex(((unsigned long int))(1));
    // xmatrix<(long int)> m(1,1);GaussJordan(m,m);
#endif
#ifdef AUROSTD_INITIALIZE_LONG_LONG_INT
    o+=initialize_scalar((long long int)(1));
    // o+=initialize_xscalar_xvector_xmatrix_xtensor((long long int)(1));
    // o+=initialize_xcomplex((long long int)(1));
    // xmatrix<(long long int)> m(1,1);GaussJordan(m,m);
#endif
#ifdef AUROSTD_INITIALIZE_UNSIGNED_LONG_LONG_INT
    o+=aurostd::string2utype<unsigned long long int>(aurostd::utype2string<unsigned long long int>((unsigned long long int)(1))+aurostd::utype2string<unsigned long long int>((unsigned long long int)(1),int()));
    // o+=initialize_scalar((unsigned long long int)(1));
    // o+=initialize_xscalar_xvector_xmatrix_xtensor((unsigned long long int)(1));
    // o+=initialize_xcomplex((unsigned long long int)(1));
    // xmatrix<(long long int)> m(1,1);GaussJordan(m,m);
#endif
  }
  return TRUE;
}


#endif  // _AURO_IMPLEMENTATIONS_

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *                                                                         *
// ***************************************************************************

