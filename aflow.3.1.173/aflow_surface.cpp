// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *                                                                         *
// ***************************************************************************
// Written by Stefano Curtarolo - Sept/Oct/Nov 2007, fixes May 2013

#include "aflow.h"

#define cdebug cerr
#define _EPS_ 0.001
#define _EPS_roundoff_ 1.0e-8
using aurostd::sign;

// prototypes

// implementations
double _sign(const double& x) {return (double) (x<0? -1:1);};
int _sign(const int& x) {return (int) (x<0? -1:1);};

#define _BBFRAC_ 1.3
#define _RRFRAC_ 1.3
#define _HKLDEF_ 4
#define _eps_    0.005
#define _oss_short_precision_aflow_surface_ 6

namespace surface {
  double PointInTriangleContribution(const xvector<double>& _point,const xvector<double>& v1,const xvector<double>& v2,const xvector<double>& v3) {
    double eps=1.1*_eps_; // relax a little bit
    // xvector<double> point(_point);
    xvector<double> point(3);
    //  if(surface::PlaneDistance(_point,v1,v2,v3)>eps)  {cerr << "too far point = " << surface::PlaneDistance(_point,v1,v2,v3) << endl;exit(0);}
    point=surface::PlaneGetProjection(_point,v1,v2,v3);
    if(distance(point,_point)>eps) {cerr << "too far point = " << distance(point,_point) << endl;exit(0);}
    // if(modulus(point-_point)>1e-8)  cout << 1e8*modulus(point-_point) << endl;  // DEBUG OK
    if(distance(v1,v2)<eps) { cerr << "ERROR - surface::PointInTriangleContribution: v1-v2" << endl;exit(0);}
    if(distance(v2,v3)<eps) { cerr << "ERROR - surface::PointInTriangleContribution: v2-v3" << endl;exit(0);}
    if(distance(v3,v1)<eps) { cerr << "ERROR - surface::PointInTriangleContribution: v3-v1" << endl;exit(0);}
    // is in vertices ?
    if(distance(point,v1)<eps) return angle(v1,v2,v3);
    if(distance(point,v2)<eps) return angle(v2,v3,v1);
    if(distance(point,v3)<eps) return angle(v3,v1,v2);
    // is in edges ?
    if(aurostd::abs(sin(v1-point,v2-point))<eps) return angle(point,v1,v2)>eps?pi:0.0;  // return 0 if not in between
    if(aurostd::abs(sin(v2-point,v3-point))<eps) return angle(point,v2,v3)>eps?pi:0.0;  // return 0 if not in between
    if(aurostd::abs(sin(v3-point,v1-point))<eps) return angle(point,v3,v1)>eps?pi:0.0;  // return 0 if not in between
    // is inside, then the angles are 360... otherwise less !
    if(aurostd::abs(angle(point,v1,v2)+angle(point,v2,v3)+angle(point,v3,v1)-2*pi)<eps) return 2.0*pi; // inside, return 2pi
    return 0.0; // not inside
  }
} // namespace surface

namespace surface {
  double PointInRhombusContribution(const xvector<double>& _point,const xvector<double>& v1,const xvector<double>& v2,const xvector<double>& v3,const xvector<double>& v4) {
    double out=0.0;
    out+=surface::PointInTriangleContribution(_point,v2,v3,v4); // double count because pick twice area in rhombi
    out+=surface::PointInTriangleContribution(_point,v3,v4,v1); // double count because pick twice area in rhombi
    out+=surface::PointInTriangleContribution(_point,v4,v1,v2); // double count because pick twice area in rhombi
    out+=surface::PointInTriangleContribution(_point,v1,v2,v3); // double count because pick twice area in rhombi
    return out/2.0;
  }
} // namespace surface

namespace surface {
  double TriangleArea(const xvector<double>& v1,const xvector<double>& v2,const xvector<double>& v3) {
    xvector<double> v12(3),v23(3),v31(3);
    v12=v2-v1;v23=v3-v2;v31=v1-v3;
    return 0.5*sqrt(scalar_product(v12,v12)*scalar_product(v31,v31)-scalar_product(v12,v31)*scalar_product(v12,v31));
  }
} // namespace surface

namespace surface {
  bool PlaneGetABCD(double& a,double& b,double& c,double& d,const xvector<double>& v1,const xvector<double>& v2,const xvector<double>& v3) {
    // *****************************************
    // PLANE
    // http://en.wikipedia.org/wiki/Plane_%28mathematics%29
    // ax+by+cz+d=0;
    a=0.0;b=0.0;c=0.0;d=0.0;
    a= (v2(2)-v1(2))*(v3(3)-v1(3))-(v2(3)-v1(3))*(v3(2)-v1(2));  d+=-v1(1)*a;
    b=-(v2(1)-v1(1))*(v3(3)-v1(3))+(v2(3)-v1(3))*(v3(1)-v1(1));  d+=-v1(2)*b;
    c= (v2(1)-v1(1))*(v3(2)-v1(2))-(v2(2)-v1(2))*(v3(1)-v1(1));  d+=-v1(3)*c;
    // daus=max(a,b,c); a/=daus;b/=daus;c/=daus;d/=daus;
    // cerr << a << " " << b << " " << c << " " << d << " "  << endl;
    // a= det(v1(2),v1(3),1,v2(2),v2(3),1,v3(2),v3(3),1.0);
    // b=-det(v1(1),v1(3),1,v2(1),v2(3),1,v3(1),v3(3),1.0);
    // c= det(v1(1),v1(2),1,v2(1),v2(2),1,v3(1),v3(2),1.0);
    // d=-det(v1(1),v1(2),v1(3),v2(1),v2(2),v2(3),v3(1),v3(2),v3(3));
    // daus=max(aurostd::abs(a),aurostd::abs(b),aurostd::abs(c));a/=daus;b/=daus;c/=daus;d/=daus;
    // cout << endl << a << " " << b << " " << c << " " << d << " "  << endl;
    return TRUE;
  }
} // namespace surface

namespace surface {
  double PlaneDistance(const xvector<double>& r,const double& a,const double& b,const double& c,const double& d) {
    return aurostd::abs(a*r(1)+b*r(2)+c*r(3)+d)/sqrt(a*a+b*b+c*c);
  }
} // namespace surface

namespace surface {
  double PlaneDistance(const xvector<double>& r,const xvector<double>& v1,const xvector<double>& v2,const xvector<double>& v3) {
    double a,b,c,d;
    surface::PlaneGetABCD(a,b,c,d,v1,v2,v3);
    return surface::PlaneDistance(r,a,b,c,d);
  }
} // namespace surface

namespace surface {
  xvector<double> PlaneGetProjection(const xvector<double>& r,const double& a,const double& b,const double& c,const double& d) {
    xvector<double> rproj(3),rorth(3);
    double dist=surface::PlaneDistance(r,a,b,c,d);
    rorth(1)=a;rorth[2]=b;rorth[3]=c;
    rproj=r-dist*rorth/sqrt(a*a+b*b+c*c);
    return rproj;
  }
} // namespace surface

namespace surface {
  xvector<double> PlaneGetProjection(const xvector<double>& r,const xvector<double>& v1,const xvector<double>& v2,const xvector<double>& v3) {
    double a,b,c,d;
    surface::PlaneGetABCD(a,b,c,d,v1,v2,v3);
    return surface::PlaneGetProjection(r,a,b,c,d);
  }
} // namespace surface

namespace surface {
  xvector<double> PlaneGetHKL(const xvector<double>& v1,const xvector<double>& v2,const xvector<double>& v3,
			      const xvector<double>& a1,const xvector<double>& a2,const xvector<double>& a3) {
    xvector<double> hkl(3);
    double eps=_eps_;
    double a,b,c,d;
    surface::PlaneGetABCD(a,b,c,d,v1,v2,v3);
    hkl(1)=-(a*a1(1)+b*a1(2)+c*a1(3))/d;
    hkl(2)=-(a*a2(1)+b*a2(2)+c*a2(3))/d;
    hkl(3)=-(a*a3(1)+b*a3(2)+c*a3(3))/d;
    if(abs(hkl(1))<eps) hkl(1)=0.0;
    if(abs(hkl(2))<eps) hkl(2)=0.0;
    if(abs(hkl(3))<eps) hkl(3)=0.0;
    return hkl;
  }
} // namespace surface

namespace surface {
  bool PlaneGetVVV(const xvector<double>& hkl,double& area,
		   xvector<double>& v1,xvector<double>& v2,xvector<double>& v3,xvector<double>& v4,
		   const xvector<double>& a1,const xvector<double>& a2,const xvector<double>& a3) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    bool isrhombus=TRUE;
    double h=hkl(1),k=hkl(2),l=hkl(3);
    double eps=_eps_;
    v1.clear();v2.clear();v3.clear();v4.clear();

    // *****************************************
    // get the 4 vertices of rhombus
    if(abs(h)>eps && abs(k)>eps && abs(l)>eps) {  // XYZ axis DEFINED
      if(LDEBUG) cerr << "XYZ axis DEFINED" << endl;
      xvector<double> vv1(3),vv2(3),vv3(3);
      v1=a1*(1.0/h);v2=a2*(1.0/k);v3=a3*(1.0/l);
      vv1=-v1+v2+v3;vv2=+v1-v2+v3;vv3=+v1+v2-v3;v1=vv1;v2=vv2;v3=vv3;
      isrhombus=FALSE;area=surface::TriangleArea(v1,v2,v3);
    }
    if(abs(h)<eps && abs(k)>eps && abs(l)>eps) {   // X axis INFINITE
      if(LDEBUG) cerr << "X axis INFINITE - 0kl" << endl;
      v1=a2*(1.0/k);v2=a3*(1.0/l);
      v3=v1+a1;v4=v2+a1;
      if(hkl!=surface::PlaneGetHKL(v1,v2,v3,a1,a2,a3)) {cerr << "ERROR - surface::PlaneGetHKL: GetPLANE: hkl problem in \"[X] axis INFINITE\"" << endl;exit(0);}
      isrhombus=TRUE;area=2.0*surface::TriangleArea(v1,v2,v3);
    }
    if(abs(h)>eps && abs(k)<eps && abs(l)>eps) {   // Y axis INFINITE
      if(LDEBUG) cerr << "Y axis INFINITE - h0l" << endl;
      v1=a1*(1.0/h);v2=a3*(1.0/l);
      v3=v1+a2;v4=v2+a2;  
      if(hkl!=surface::PlaneGetHKL(v1,v2,v3,a1,a2,a3)) {cerr << "ERROR - surface::PlaneGetHKL: GetPLANE: hkl problem in \"[Y] axis INFINITE\"" << endl;exit(0);}
      isrhombus=TRUE;area=2.0*surface::TriangleArea(v1,v2,v3);
    }
    if(abs(h)>eps && abs(k)>eps && abs(l)<eps) {   // Z axis INFINITE
      if(LDEBUG) cerr << "Z axis INFINITE - hk0" << endl;
      v1=a1*(1.0/h);v2=a2*(1.0/k);
      v3=v1+a3;v4=v2+a3;
      if(hkl!=surface::PlaneGetHKL(v1,v2,v3,a1,a2,a3)) {cerr << "ERROR - surface::PlaneGetHKL: GetPLANE: hkl problem in \"[Z] axis INFINITE\"" << endl;exit(0);}
      isrhombus=TRUE;area=2.0*surface::TriangleArea(v1,v2,v3);
    }
    if(abs(h)>eps && abs(k)<eps && abs(l)<eps) {   // YZ axis INFINITE
      if(LDEBUG) cerr << "YZ axis INFINITE - h00" << endl;
      v1=a1*(1.0/h);v2=v1+a2;v3=v1+a3;v4=v1+a2+a3;
      isrhombus=TRUE;area=2.0*surface::TriangleArea(v1,v2,v3);
    }
    if(abs(h)<eps && abs(k)>eps && abs(l)<eps) {   // XZ axis INFINITE
      if(LDEBUG) cerr << "XZ axis INFINITE - 0k0" << endl;
      //v1=a2*(1.0/k);v2=v1+a1;v3=v1+a3;v4=v1+a1+a3;
      v2=a2*(1.0/k);v1=v2+a1;v3=v2+a3;v4=v2+a1+a3;
      isrhombus=TRUE;area=2.0*surface::TriangleArea(v1,v2,v3);
    }
    if(abs(h)<eps && abs(k)<eps && abs(l)>eps) {   // XY axis INFINITE
      if(LDEBUG) cerr << "XY axis INFINITE - 00l" << endl;
      //v1=a3*(1.0/l);v2=v1+a1;v3=v1+a2;v4=v1+a1+a2;
      v3=a3*(1.0/l);v1=v3+a1;v2=v3+a2;v4=v3+a1+a2;
      isrhombus=TRUE;area=2.0*surface::TriangleArea(v1,v2,v3);
    }
    if(abs(h)<eps && abs(k)<eps && abs(l)<eps) {   // XYZ axis INFINITE
      cerr << "ERROR - surface::PlaneGetHKL: aflow_xsurface: h,k,l cannot be 0 0 0" << endl;
      exit(0);
    }
    for(int i=1;i<=3;i++) {
      if(abs(v1[i])<eps) v1[i]=0.0;
      if(abs(v2[i])<eps) v2[i]=0.0;
      if(abs(v3[i])<eps) v3[i]=0.0;
      if(abs(v4[i])<eps) v4[i]=0.0;
    }
    return isrhombus;
  }
} // namespace surface

namespace surface {
  double GetPlaneDensityAtoms(const xstructure& _str,const xvector<double>& hkl,const double& roughness,const int& type_at) {
    xmatrix<double> lattice(3,3);
    lattice=(_str.lattice);
    xvector<double> a1(3),a2(3),a3(3);                    // lattice vectors
    a1=lattice(1);a2=lattice(2);a3=lattice(3);            // a1,a2,a3 are the rows of the lattice matrix
    xvector<double> v1(3),v2(3),v3(3),v4(3);
    double area;
    xstructure str(_str);
    str=BringInCell(str);

    bool isrhombus;

    vector<double> afound(_str.num_each_type.size());
  
    isrhombus=PlaneGetVVV(hkl,area,v1,v2,v3,v4,a1,a2,a3);

    if(0) {
      cout << v1 << " " << v2 << " " << v3 << " " << v4 << " " << endl;
    }
    // *****************************************
    // PLANE ax+by+cz+d=0;
    double a,b,c,d;
    PlaneGetABCD(a,b,c,d,v1,v2,v3);

    // *****************************************
    // create search over all atoms.
    // distance http://en.wikipedia.org/wiki/Plane_(mathematics)
    xvector<int> dims(3);
    double radius,dist;
    xvector<double> rrr(3);
    radius=1.1*max(modulus(v1),modulus(v2),modulus(v3),modulus(v4))+max(modulus(a1),modulus(a2),modulus(a3));
    dims=LatticeDimensionSphere(lattice,radius);

    for(uint itype=0;itype<afound.size();itype++)
      afound.at(itype)=0.0;
  
    //  for(int ii=-dims(1);ii<=dims(1);ii++) {
    //   for(int jj=-dims(2);jj<=dims(2);jj++) {
    //    for(int kk=-dims(3);kk<=dims(3);kk++) {
    for(int ii=dims(1);ii>=-dims(1);ii--) {
      for(int jj=dims(2);jj>=-dims(2);jj--) {
	for(int kk=dims(3);kk>=-dims(3);kk--) {
	  for(uint iat=0;iat<str.atoms.size();iat++) {
	    rrr=((double)ii)*a1+((double)jj)*a2+((double)kk)*a3+str.atoms.at(iat).cpos;
	    dist=aurostd::abs(a*rrr(1)+b*rrr(2)+c*rrr(3)+d)/sqrt(a*a+b*b+c*c);
	    if(dist<roughness) {
	      if(!isrhombus) afound.at(str.atoms.at(iat).type)+=surface::PointInTriangleContribution(rrr,v1,v2,v3);   // 1st triangle
	      if(isrhombus)  afound.at(str.atoms.at(iat).type)+=surface::PointInRhombusContribution(rrr,v1,v2,v3,v4); // rhombus
	    }
	  }
	}
      }
    }
    if(type_at<0) {
      double out=0;
      for(uint itype=0;itype<str.num_each_type.size();itype++)
	out+=afound.at(itype);
      return out/(2*pi*area*str.scale*str.scale);
    } else {
      return afound.at(type_at)/(2*pi*area*str.scale*str.scale);
    }
  }
} // namespace surface

namespace surface {
  double GetPlaneDensityAtoms(const xstructure& _str,const xvector<double>& hkl,const double& roughness) {
    return GetPlaneDensityAtoms(_str,hkl,roughness,-1);
  }
} // namespace surface

namespace surface {
  double GetPlaneDensityBBonds(const xstructure& _str,const xvector<double>& hkl,const double& roughness,const double& bbdistance,const int& type_at1,const int& type_at2) {
    xmatrix<double> lattice(3,3);
    lattice=(_str.lattice);
    xvector<double> a1(3),a2(3),a3(3);                    // lattice vectors
    a1=lattice(1);a2=lattice(2);a3=lattice(3);            // a1,a2,a3 are the rows of the lattice matrix
    xvector<double> v1(3),v2(3),v3(3),v4(3);
    double area;
    xstructure str(_str);
    str=BringInCell(str);
    bool isrhombus;
    double a,b,c,d;
    double afound=0.0;
    double eps=_eps_;
  
    if(roughness) {;}  // phony, just to use

    vector<xvector<double>*> grid_far_atoms_cpos,grid_close_atoms_cpos;
    vector<int> grid_far_atoms_type,grid_close_atoms_type;
    xvector<double> *grid_atoms_cpos_ptr;

    isrhombus=PlaneGetVVV(hkl,area,v1,v2,v3,v4,a1,a2,a3);
    PlaneGetABCD(a,b,c,d,v1,v2,v3);
    //  cerr << "a=" << a << " b=" << b << " c=" << c << " d=" << d << endl;

    // *****************************************
    // create search over all atoms.
    // distance http://en.wikipedia.org/wiki/Plane_(mathematics)
    xvector<int> dims(3);
    double radius,dist,num,den,u;
    xvector<double> rrr(3),rrr1(3),rrr2(3);
    radius=1.1*max(modulus(v1),modulus(v2),modulus(v3),modulus(v4))+max(modulus(a1),modulus(a2),modulus(a3));
    dims=LatticeDimensionSphere(lattice,radius);

    for(int ii=dims(1);ii>=-dims(1);ii--) {
      for(int jj=dims(2);jj>=-dims(2);jj--) {
	for(int kk=dims(3);kk>=-dims(3);kk--) {
	  for(uint iat=0;iat<str.atoms.size();iat++) {
	    rrr=((double)ii)*a1+((double)jj)*a2+((double)kk)*a3+str.atoms.at(iat).cpos;
	    dist=aurostd::abs(a*rrr(1)+b*rrr(2)+c*rrr(3)+d)/sqrt(a*a+b*b+c*c);
	    if(dist<=bbdistance) {
	      if(dist>0.01) {   // FAR
		grid_atoms_cpos_ptr = new xvector<double>(3);
		*grid_atoms_cpos_ptr=rrr;                    
		grid_far_atoms_cpos.push_back(grid_atoms_cpos_ptr);  
		grid_far_atoms_type.push_back(str.atoms.at(iat).type);
	      } else { // CLOSE
		grid_atoms_cpos_ptr = new xvector<double>(3);
		*grid_atoms_cpos_ptr=rrr;                    
		grid_close_atoms_cpos.push_back(grid_atoms_cpos_ptr);  
		grid_close_atoms_type.push_back(str.atoms.at(iat).type);
	      }
	    }
	  }
	}
      }
    }
    //  cout << grid_close_atoms_cpos.size() << " " << grid_far_atoms_cpos.size() << endl;
    int matches=0;
    afound=0.0;
    for(uint iat_close=0;iat_close<grid_close_atoms_cpos.size();iat_close++) {
      rrr1=*grid_close_atoms_cpos.at(iat_close);
      for(uint iat_far=0;iat_far<grid_far_atoms_cpos.size();iat_far++) {
	rrr2=*grid_far_atoms_cpos.at(iat_far);
	if((type_at1<0 && type_at2<0) ||
	   (grid_close_atoms_type.at(iat_close)==type_at1 && grid_far_atoms_type.at(iat_far)==type_at2) ||
	   (grid_close_atoms_type.at(iat_close)==type_at2 && grid_far_atoms_type.at(iat_far)==type_at1)) {
	  dist=distance(rrr2,rrr1);
	  if(dist<=bbdistance && dist>0.2) { // they are close... but not the same
	    // intersection http://local.wasp.uwa.edu.au/~pbourke/geometry/planeline/  
	    // P = P1 + u (P2 - P1)    rrr = rrr1 + u (rrr2 - rrr1)    
	    num=a*rrr1(1)+b*rrr1[2]+c*rrr1[3]+d;
	    den=a*(rrr1(1)-rrr2(1))+b*(rrr1[2]-rrr2[2])+c*(rrr1[3]-rrr2[3]);
	    // if(aurostd::abs(den)<eps/10.0 && aurostd::abs(num)>eps/10.0) {cerr << num << " " << den << endl;exit(0);}
	    // if(aurostd::abs(den)<eps/10.0 && aurostd::abs(num)<=aurostd::abs(den)) {num=0.0;den=1.0;}
	    if(aurostd::abs(den)>eps/10.0) // to avoid rrr1,rrr2,rrr coplanar with the plane
	      { // found
		u=num/den;
		if(u<-2*eps || u>1.0+2*eps) cerr << "u=" << u << " num=" << num << " den=" << den << endl;
		// if(u<0.0 && u>=-eps) u=0.0;
		// if(u>1.0 && u<=1.0+eps) u=1.0;
		// if(u<-eps || u>1+eps) {cerr << "u=" << u << endl; exit(0);}
		if(u<0.0) u=0.0;
		if(u>1.0) u=1.0;
		rrr=rrr1+u*(rrr2-rrr1);
		//  point=PlaneGetProjection(_point,v1,v2,v3);
		if(PlaneDistance(rrr,a,b,c,d)<eps) {
		  if(!isrhombus) afound+=surface::PointInTriangleContribution(rrr,v1,v2,v3);   // 1st triangle
		  if(isrhombus)  afound+=surface::PointInRhombusContribution(rrr,v1,v2,v3,v4); // rhombus
		  matches++;
		}
	      }
	  }
	}
      }
    }
    // EXIT
    // cerr << matches << " " << afound/(2*pi*area*str.scale*str.scale) << endl;
    for(uint i=0;i<grid_far_atoms_cpos.size();i++)
      delete grid_far_atoms_cpos[i]; 
    grid_far_atoms_cpos.clear();
    for(uint i=0;i<grid_close_atoms_cpos.size();i++)
      delete grid_close_atoms_cpos[i]; 
    grid_close_atoms_cpos.clear();
    return afound/(2*pi*area*str.scale*str.scale);
  }
} // namespace surface

namespace surface {
  double GetPlaneDensityBBonds(const xstructure& _str,const xvector<double>& hkl,const double& roughness,const double& bbdistance) {
    return surface::GetPlaneDensityBBonds(_str,hkl,roughness,bbdistance,-1,-1);
  }
} // namespace surface

namespace surface {
  double GetNNeighbours(const xstructure& _str,const int& type_at1,const int& type_at2) {
    xstructure str(_str);
    str=ReScale(BringInCell(_str),1.0);
    xvector<double> a1(3),a2(3),a3(3),rrr1(3),rrr2(3);                // lattice vectors and vectors
    a1=str.lattice(1);a2=str.lattice(2);a3=str.lattice(3);            // a1,a2,a3 are the rows of the lattice matrix
    double r,nndistance,radius=1.5*max(modulus(a1),modulus(a2),modulus(a3));
    xvector<int> dims(3);
    dims=LatticeDimensionSphere(str.lattice,radius);
    nndistance=10.0*radius;
    // for(int i1=-dims(1);i1<=dims(1);i1++)
    // for(int j1=-dims(2);j1<=dims(2);j1++)
    // for(int k1=-dims(3);k1<=dims(3);k1++)
    int i1=0,j1=0,k1=0;
    for(uint iat1=0;iat1<str.atoms.size();iat1++) {
      rrr1=((double)i1)*a1+((double)j1)*a2+((double)k1)*a3+str.atoms.at(iat1).cpos;
      for(int i2=-dims(1);i2<=dims(1);i2++)
	for(int j2=-dims(2);j2<=dims(2);j2++)
	  for(int k2=-dims(3);k2<=dims(3);k2++)
	    for(uint iat2=0;iat2<str.atoms.size();iat2++) {
	      if((type_at1<0 && type_at2<0) ||
		 (str.atoms.at(iat1).type==type_at1 && str.atoms.at(iat2).type==type_at2) ||
		 (str.atoms.at(iat1).type==type_at2 && str.atoms.at(iat2).type==type_at1)) {
		rrr2=((double)i2)*a1+((double)j2)*a2+((double)k2)*a3+str.atoms.at(iat2).cpos;
		r=modulus(rrr1-rrr2);
		if(r<=nndistance && r>_eps_) nndistance=r;
	      }
	    }
    }
    return nndistance;
  }
} // namespace surface

namespace surface {
  double GetNNeighbours(const xstructure& _str) {
    return surface::GetNNeighbours(_str,-1,-1);
  }
} // namespace surface

namespace surface {
  string PrintHKLSigma(int num_types,int num_types_combinations) {
    ostringstream aus;
    if(num_types==1) {
      aus << "     h           k           l        sigma(#/AA)  Nb(#/AA)    " << endl;
    } else {
      // 1st line
      aus << "     h           k           l        ";
      for(int j=0;j<=num_types;j++) aus << "sigma(#/AA) ";
      for(int j=0;j<=num_types_combinations;j++) aus << " Nb(#/AA)   ";
      aus << endl;
      // 2nd line
      aus << "                                      ";
      aus << "   T=(*)    ";      
      for(int j=0;j<num_types;j++) aus << "   T=(" << j << ")    ";
      aus << " TT=(*-*)   ";
      for(int it1=0;it1<num_types;it1++)
	for(int it2=it1;it2<num_types;it2++)
	  aus << " TT=(" << it1 << "-" << it2 << ")   ";
      aus << endl;
      // 3rd line
    }
    return aus.str();
  }
} // namespace surface

namespace surface {
  string PrintHKLSigmaBB(int num_types,int num_types_combinations,const double& bbfrac,const double& bbdistance,const xmatrix<double>& bbdistances) {
    ostringstream aus;
    aus.setf(std::ios::fixed,std::ios::floatfield);
    aus.precision(4);
    aus << surface::PrintHKLSigma(num_types,num_types_combinations);
    if(num_types==1) {
    } else {
      // 3rd line
      aus << "                                       ";      
      for(int j=0;j<num_types+1;j++) aus << "            ";
      aus << "b=" << bbdistance << "    ";
      for(int it1=0;it1<num_types;it1++)
	for(int it2=it1;it2<num_types;it2++)
	  aus << "b=" << bbdistances(it1,it2) << "    ";
      aus << endl;
    }
    if(bbfrac) {;}  // phony, to keep bbfrac busy
    return aus.str();
  }
} // namespace surface

namespace surface {
  string PrintNNdists(int num_types,int num_types_combinations,const double& rrdist,const double& nndist,const xmatrix<double>& nndists) {
    ostringstream aus;
    aus.setf(std::ios::fixed,std::ios::floatfield);
    aus.precision(5);
    aus << "nndist(*-*)=" << nndist << endl;
    if(num_types>1) {
      for(int it1=0;it1<num_types;it1++)
	for(int it2=it1;it2<num_types;it2++)
	  aus << "nndists(" << it1 << "-" << it2 << ")=" << nndists(it1,it2) << " ";
      aus << endl;
    }
    if(num_types_combinations) {;}  // phony, to keep num_types_combinations busy
    if(rrdist) {;}                  // phony, to keep rrdist busy
    return aus.str();
  }
} // namespace surface

namespace surface {
  bool AddPlaneHKL(const xstructure& str,const double& h,const double& k,const double&l,const double& roughness,const double& bbdistance,const xmatrix<double>& bbdistances,const double& eps,const int& jossmax,vector<double>& plane,vector<vector<double> >& planes,const bool& osswrite1,ostream& oss1,const bool& osswrite2,ostream& oss2) {
    double density,bbonds;
    xvector<double> hkl(3);
    bool plane_added=FALSE,hkl_found=FALSE;
    int i,num_types=str.num_each_type.size();
    hkl(1)=h;hkl(2)=k;hkl(3)=l;
    if(abs(hkl(1))<eps && abs(hkl(2))<eps && abs(hkl(3))<eps) return FALSE;
    if((abs(hkl(1))>=1.0 || abs(hkl(1))<eps) &&
       (abs(hkl(2))>=1.0 || abs(hkl(2))<eps) &&
       (abs(hkl(3))>=1.0 || abs(hkl(3))<eps)) {
      for(uint ii=0;ii<planes.size()&&!hkl_found;ii++)
	hkl_found=(abs(h-planes[ii][0])<eps && abs(k-planes[ii][1])<eps && abs(l-planes[ii][2])<eps);
      if(hkl_found==FALSE) {                                 // new operation, generate and save it
	// ---------------------- hkl ----------------------
	plane.at(0)=h;plane.at(1)=k;plane.at(2)=l;
	// ---------------------- density ------------------
	density=surface::GetPlaneDensityAtoms(str,hkl,roughness);
	i=3;plane.at(i++)=density;
	if(num_types>1)
	  for(int it1=0;it1<num_types;it1++)
	    plane.at(i++)=surface::GetPlaneDensityAtoms(str,hkl,roughness,it1);
	// ---------------------- bbonds -------------------
	bbonds=surface::GetPlaneDensityBBonds(str,hkl,roughness,bbdistance);
	plane.at(i++)=bbonds;
	if(num_types>1)
	  for(int it1=0;it1<num_types;it1++)
	    for(int it2=it1;it2<num_types;it2++)
	      plane.at(i++)=surface::GetPlaneDensityBBonds(str,hkl,roughness,bbdistances(it1,it2),it1,it2);
	if(density>0.0) {
	  planes.push_back(plane); // save them
	  plane_added=TRUE;
	  //    oss1 << "*";oss1.flush();
	}
      }
      if(plane_added) {
	for(int j=0;j<=jossmax;j++) {
	  if(osswrite1) oss1 << (plane[j]>=0?"  ":" ") << (abs(plane[j])<10.0?" ":"")<< plane[j] << " ";
	  if(osswrite2) oss2 << (plane[j]>=0?"  ":" ") << (abs(plane[j])<10.0?" ":"")<< plane[j] << " ";
	}
	if(osswrite1) oss1 << planes.size() << " ";
	if(osswrite2) oss2 << planes.size() << " ";
	if(osswrite1) oss1 << endl;
	if(osswrite1) oss1.flush();
	if(osswrite2) oss2 << endl;
	if(osswrite2) oss2.flush();
      }    
    }
    return plane_added;
  }
} // namespace surface

namespace surface {
  bool AddPlaneHKL(const xstructure& str,const xvector<double>& hkl,const double& roughness,const double& bbdistance,const xmatrix<double>& bbdistances,
		   const double& eps,const int& jossmax,vector<double>& plane,vector<vector<double> >& planes,const bool& osswrite1,ostream& oss1,const bool& osswrite2,ostream& oss2) {
    return surface::AddPlaneHKL(str,hkl(1),hkl(2),hkl(3),roughness,bbdistance,bbdistances,eps,jossmax,plane,planes,osswrite1,oss1,osswrite2,oss2);
  }
} // namespace surface

namespace surface {
  bool ReducePrintSurfaces(const xstructure& str,const double& eps,
			   vector<vector<double> >& planesreducible,vector<vector<double> >& planesirreducible,vector<vector<uint> >& planesirreducible_images,
			   const bool& osswrite1,ostream& oss1,const bool& osswrite2,ostream& oss2) {
    uint num_types=str.num_each_type.size(),num_types_combinations=num_types*(num_types+1)/2;
    uint jossmax=4;
    if(num_types_combinations>1) jossmax+=num_types+num_types_combinations;
    vector<double> plane(3+(1+num_types)+(1+num_types_combinations)); // 3 for hkl, 1 for all types and all combinations
    double rdensity,rbbonds,idensity,ibbonds;
    bool hkl_found;

    xvector<double> rhkl(3),ihkl(3);
  
    sort(planesreducible.begin(),planesreducible.end(),aurostd::_sort_double_value012());
    sort(planesreducible.begin(),planesreducible.end(),aurostd::_isort_double_value3());

    //  if(osswrite1) oss1 << "PREFORMING REDUCTION wait (" << planesreducible.size() << " to reduce)" << endl;
    if(osswrite1) oss1 << surface::PrintHKLSigma(num_types,num_types_combinations);
    if(osswrite2) oss2 << surface::PrintHKLSigma(num_types,num_types_combinations);
    // THIS ROUTINE SHOULD BE CHANCED IN MULTITHREADS TO SPEED UP THE REDUCTION
    for(uint i=0;i<planesreducible.size();i++) {
      rhkl(1)=planesreducible[i][0];rhkl(2)=planesreducible[i][1];rhkl(3)=planesreducible[i][2];
      rdensity=planesreducible[i][3];rbbonds =planesreducible[i][4];
      if(rdensity>0.0) { //must make sense
	hkl_found=FALSE;
	for(uint ii=0;ii<planesirreducible.size()&&!hkl_found;ii++) { // check if triplet of vectors are equivalent !
	  ihkl(1)=planesirreducible[ii][0];ihkl(2)=planesirreducible[ii][1];ihkl(3)=planesirreducible[ii][2];
	  idensity=planesirreducible[ii][3];ibbonds =planesirreducible[ii][4];
	  if(abs(rdensity-idensity)<eps && abs(rbbonds-ibbonds)<eps) // check that the density and bbonds makes sense
	    for(uint sg=0;sg<str.fgroup.size()&&!hkl_found;sg++) {
	      if(aurostd::modulus(str.fgroup[sg].ftau)<eps) // only for non shift, otherwise no meaning
		// if(aurostd::modulus(SYM_ApplyFpos(ihkl,str.fgroup[sg],str,FALSE)-rhkl)<eps) hkl_found=TRUE;
		if(aurostd::modulus(str.fgroup[sg].Uf*ihkl-rhkl)<eps) {
		  hkl_found=TRUE; // same but faster than the other one
		}
	    }
	}
	if(hkl_found==FALSE) {                                     // new irreducible operation, generate and save it
	  planesirreducible.push_back(planesreducible[i]);
	  planesirreducible_images.push_back(vector<uint>(0));
	}
	// if(!planesirreducible_images.size()) planesirreducible_images.push_back(vector<uint>(0)); // safety, do we need it ?
	planesirreducible_images.back().push_back(i);
      }
    }
  
    for(uint i=0;i<planesirreducible.size();i++) { // new irreducible operation, generate and save it
      //        if(osswrite1) oss1 << "*";if(osswrite1) oss1.flush();
      for(uint j=0;j<=jossmax;j++) 
	if(osswrite1) oss1 << (planesirreducible.at(i).at(j)>=0?"  ":" ") << (abs(planesirreducible.at(i).at(j))<10.0?" ":"") << planesirreducible.at(i).at(j) << " ";
      for(uint j=0;j<=4;j++)
	if(osswrite2) oss2 << (planesirreducible.at(i).at(j)>=0?"  ":" ") << (abs(planesirreducible.at(i).at(j))<10.0?" ":"") << planesirreducible.at(i).at(j) << " ";
      if(osswrite1) oss1 << planesirreducible.size() << " ";
      if(osswrite2) oss2 << planesirreducible.size() << " ";
      if(osswrite1) oss1 << endl;
      if(osswrite2) oss2 << endl;
      if(osswrite1) {
	oss1 << "           equivalent family " << endl;
	for(uint k=0;k<planesirreducible_images.at(i).size();k++) {
	  for(uint j=0;j<3;j++) {
	    oss1 << "" 
		 << (planesreducible.at(planesirreducible_images.at(i).at(k)).at(j)>=0?"  ":" ") 
		 << (abs(planesreducible.at(planesirreducible_images.at(i).at(k)).at(j))<10.0?" ":"")
		 << planesreducible.at(planesirreducible_images.at(i).at(k)).at(j) << " ";
	  }
	  oss1 << endl;
	}
      }
      //   if(osswrite1) oss1 << endl;
    }
    // if(osswrite1) oss1 << endl;
    // if(osswrite2) oss2 << endl;
    // routine to make the planes as POSITIVE as POSSIBLE
    // ------------------------------------------------------
    // CODE FOR NUM_TYPES AND NUM_TYPES_COMBINATIONS
    return TRUE;
  }
} // namespace surface

// ------------------------------------------------------------------------------------------------------------------
// only one HKL
namespace surface {
  bool GetSurfaceHKL(const xstructure& _str,_aflags& aflags,const xvector<double>& iparams,
		     vector<vector<double> >& planesreducible,vector<vector<double> >& planesirreducible,
		     ostream& oss) {
    bool Krun=TRUE;
    xvector<double> hkl(3);
    int num_types=_str.num_each_type.size(),num_types_combinations=num_types*(num_types+1)/2;
    int jossmax=4,mode=iparams.rows;
    if(num_types_combinations>1) jossmax+=num_types+num_types_combinations;
    vector<double> plane(3+(1+num_types)+(1+num_types_combinations)); // 3 for hkl, 1 for all types and all combinations
    xstructure str(_str);
    // str=BringInCell(_str);
    str=ReScale(BringInCell(_str),1.0);
    // NO ORIGIN if(1) str.ShifOriginToAtom(0);
    double bbfrac=_BBFRAC_;
    // double rrfrac=_RRFRAC_;
    double roughness=_eps_/2.0;
    double bbdistance;
    double density,bbonds,area; //
    xmatrix<double> lattice(3,3);
    lattice=(str.lattice);
    xvector<double> a1(3),a2(3),a3(3);                    // lattice vectors
    a1=lattice(1);a2=lattice(2);a3=lattice(3);            // a1,a2,a3 are the rows of the lattice matrix
    xvector<double> v1(3),v2(3),v3(3),v4(3);              // vectors for symmetry search

    oss.setf(std::ios::fixed,std::ios::floatfield);
    oss.precision(_oss_short_precision_aflow_surface_);

    double nndist=surface::GetNNeighbours(str);
    xmatrix<double> nndists(num_types-1,num_types-1,0,0),bbdistances(num_types-1,num_types-1,0,0);
    for(int it1=0;it1<num_types;it1++)
      for(int it2=it1;it2<num_types;it2++)
	nndists(it1,it2)=surface::GetNNeighbours(str,it1,it2);


    //  if(mode==3 || mode==4) { // all three are given
    if(mode!=3 && mode!=4) {cerr << "only mode 3 and 4 are defined" << endl;exit(0);}
    if(mode==3 || mode==4) {
      if(mode==3) bbfrac=_BBFRAC_;
      if(mode==4) bbfrac=iparams(4);
      bbdistance=bbfrac*nndist;
      bbdistances=bbfrac*nndists;
    
      oss << "HKL CALCULATION" << endl;
      oss << surface::PrintNNdists(num_types,num_types_combinations,bbfrac,nndist,nndists);
      oss << surface::PrintHKLSigmaBB(num_types,num_types_combinations,bbfrac,bbdistance,bbdistances);
      int i=0;
      //    double h,k,l;
      hkl=iparams;
      // hkl
      plane.at(0)=hkl(1);plane.at(1)=hkl(2);plane.at(2)=hkl(3);
      // density
      density=surface::GetPlaneDensityAtoms(str,hkl,roughness);
      i=3;plane.at(i++)=density;
      if(num_types>1)
	for(int it1=0;it1<num_types;it1++)
	  plane.at(i++)=surface::GetPlaneDensityAtoms(str,hkl,roughness,it1);
      // bbonds
      bbonds=surface::GetPlaneDensityBBonds(str,hkl,roughness,bbdistance);
      plane.at(i++)=bbonds;
      if(num_types>1)
	for(int it1=0;it1<num_types;it1++)
	  for(int it2=it1;it2<num_types;it2++)
	    plane.at(i++)=surface::GetPlaneDensityBBonds(str,hkl,roughness,bbdistances(it1,it2),it1,it2);
      for(int j=0;j<=jossmax;j++)
	oss << (plane[j]>=0?"  ":" ") << (abs(plane[j])<10.0?" ":"")<< plane[j] << " ";
      PlaneGetVVV(hkl,area,v1,v2,v3,v4,a1,a2,a3); // oss << v1 << " " << v2 << " " << v3 << " " << v4 << " ";
      //  oss << (max(modulus(v1),modulus(v2),modulus(v3),modulus(v4))-bbdistance)/modulus(a1+a2+a3) << " ";
      //  double radius=max(modulus(v1),modulus(v2),modulus(v3),modulus(v4));
      //  oss << radius << " " << LatticeDimensionSphere(lattice,radius) << " ";
      oss << endl;
      planesreducible.push_back(plane); // save them
      planesirreducible.push_back(plane); // save them
      return Krun;
    }
    if(aflags.QUIET) {;} // phony just to keep aflags busy
    return FALSE;
  }
} // namespace surface

// ------------------------------------------------------------------------------------------------------------------
// Search HKL the trivial/simple/complete
namespace surface {
  bool GetSurfaceHKLSearch(const xstructure& _str,_aflags& aflags,const xvector<double>& iparams,
			   vector<vector<double> >& planesreducible,vector<vector<double> >& planesirreducible,vector<vector<uint> >& planesirreducible_images,
			   ostream& oss,const string& smode) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    bool search_trivial=FALSE,search_simple=FALSE,search_complete=FALSE;
  
    if(smode!="HKL_SEARCH_TRIVIAL" && smode!="HKL_SEARCH_SIMPLE" && smode!="HKL_SEARCH_COMPLETE") {cerr << "Error: surface::GetSurfaceHKLSearch [1]" << endl; exit(0);}
    if(smode=="HKL_SEARCH_TRIVIAL")  {search_trivial=TRUE;search_simple=FALSE;search_complete=FALSE;};
    if(smode=="HKL_SEARCH_SIMPLE" )  {search_trivial=FALSE;search_simple=TRUE;search_complete=FALSE;};
    if(smode=="HKL_SEARCH_COMPLETE") {search_trivial=FALSE;search_simple=FALSE;search_complete=TRUE;};
    if(LDEBUG) cerr << "surface::GetSurfaceHKLSearch: SURFACE begin " <<   endl;
    bool Krun=TRUE;
    double eps=_eps_;
    double roughness=_eps_/2.0;
    double radius,bbdistance,step=0,hklmax=0;
    xvector<int> dims(3);
    xvector<double> dims_hklmax(3),hkl(3);
    int jossmax=4,mode=iparams.rows;
    if(LDEBUG) cerr << "surface::GetSurfaceHKLSearch: [1] mode=" << mode <<  endl;
    if(LDEBUG) cerr << "surface::GetSurfaceHKLSearch: [1] iparams=" << iparams <<  endl;

    int num_types=_str.num_each_type.size(),num_types_combinations=num_types*(num_types+1)/2;
    if(num_types_combinations>1) jossmax+=num_types+num_types_combinations;
    vector<double> plane(3+(1+num_types)+(1+num_types_combinations)); // 3 for hkl, 1 for all types and all combinations

    if(LDEBUG) cerr << "surface::GetSurfaceHKLSearch: [1] " <<   endl;

    string banner="--------------------------------------------------------------";
    for(int i=1;i<=num_types+num_types_combinations&&num_types>1;i++) banner=banner+"------------";
    //  oss << banner << endl; // "------------------------------------------------------------------------------------" << endl;
    // oss << "HKL CALCULATION - Stefano Curtarolo" << endl;
    // if(search_trivial) oss << "TRIVIAL SEARCH" << endl;
    // if(search_simple)  oss << "SIMPLE SEARCH" << endl;
    // if(search_complete) oss << "COMPLETE SEARCH" << endl;
    //  oss << "SURFACE begin " <<   endl;

    xstructure str(_str);
    // str=BringInCell(_str);
    str=ReScale(BringInCell(_str),1.0);
    // str.MinkowskiBasisReduction();
    //  if(1) str.ShifOriginToAtom(0);
    double bbfrac=_BBFRAC_;
    double rrfrac=_RRFRAC_;
    double area,determinant,a,b,c,d;

    xmatrix<double> lattice(3,3);lattice=(str.lattice);
    xvector<double> a1(3),a2(3),a3(3);a1=lattice(1);a2=lattice(2);a3=lattice(3);   // a1,a2,a3 are the rows of the lattice matrix
    xvector<double> v1(3),v2(3),v3(3),v4(3);              // vectors for symmetry search
    xvector<double> rrr(3),rrr1(3),rrr2(3),rrr3(3);
    _atom aaa,aaa1,aaa2,aaa3;
    vector<xvector<double>*> grid_atoms_cpos;
    // vector<xvector<double>*> grid_atoms_fpos;  //USELESS
    vector<int> grid_atoms_number;
    xvector<double> *grid_atoms_cpos_ptr;
    bool PFSWRITE=TRUE;ofstream FileDevNull("/dev/null");
    bool FFFflag=TRUE;
    bool hkl_found;
    int number1,number2,number3;
    uint grid_atoms_size;
    int imin=0,imax=0,jmin=0,jmax=0,kmin=0,kmax=0;
  
    aflags.QUIET=TRUE;
    bool OSSWRITE=FALSE;//TRUE;
  
    oss.setf(std::ios::fixed,std::ios::floatfield);
    oss.precision(_oss_short_precision_aflow_surface_);
    ofstream FFF;
    string FileNameSURFACE;
    if(FFFflag) FileNameSURFACE=_AFLOW_SURFACE_FILE_;
    else FileNameSURFACE="/dev/null";
    FFF.open(FileNameSURFACE.c_str(),std::ios::out);
    FFF.setf(std::ios::fixed,std::ios::floatfield);
    FFF.precision(_oss_short_precision_aflow_surface_);
  
    double nndist=surface::GetNNeighbours(str);
    xmatrix<double> nndists(num_types-1,num_types-1,0,0),bbdistances(num_types-1,num_types-1,0,0);
    for(int it1=0;it1<num_types;it1++)
      for(int it2=it1;it2<num_types;it2++)
	nndists(it1,it2)=surface::GetNNeighbours(str,it1,it2);

    if(LDEBUG) cerr << "surface::GetSurfaceHKLSearch: [2] " <<   endl;

    if(search_trivial) {
      if(mode==0 || mode==1 || mode==2 || mode==3) { // seek
	//  oss << "SURFACE begin " <<   endl;
	if(mode==0) {hklmax=_HKLDEF_;bbfrac=_BBFRAC_;step=1.0;}
	if(mode==1) {hklmax=abs(iparams(1));bbfrac=_BBFRAC_;step=1.0;}
	if(mode==2) {hklmax=abs(iparams(1));bbfrac=iparams(2);step=1.0;}
	if(mode==3) {hklmax=abs(iparams(1));bbfrac=iparams(2);step=abs(iparams(3));}
	dims_hklmax(1)=hklmax;dims_hklmax(2)=hklmax;dims_hklmax(3)=(double) hklmax;
	bbdistance=bbfrac*nndist;
	bbdistances=bbfrac*nndists;

	oss << banner << endl; // ----------------------------------------------------------------
	oss << aflow::Banner("BANNER_TINY") << endl;
	oss << banner << endl; // ----------------------------------------------------------------
	oss << "HKL CALCULATION" << endl;
	if(search_simple)   oss << "SIMPLE SEARCH" << endl;
	if(search_trivial)  oss << "TRIVIAL SEARCH" << endl;
	if(search_complete) oss << "COMPLETE SEARCH" << endl;
	str.LatticeReduction_avoid=TRUE;  // DOES NOT DO LATTICE REDUCTION // NIGGLI and MINK
	str.sgroup_radius=1.05*RadiusSphereLattice(lattice);                                      //CO 171024 - new sym framework
        _kflags kflags; pflow::defaultKFlags4SymCalc(kflags,true);                                //CO 171024 - new sym framework
        pflow::defaultKFlags4SymWrite(kflags,PFSWRITE); kflags.KBIN_SYMMETRY_SGROUP_WRITE=false;  //CO 171024 - new sym framework
        pflow::PerformFullSymmetry(str,FileDevNull,aflags,kflags,OSSWRITE,oss);                   //CO 171024 - new sym framework
	//SYM::CalculatePointGroup(FileDevNull,str,aflags,PFSWRITE,OSSWRITE,oss);                 //CO 171024 - new sym framework
	//SYM::CalculateSitePointGroup(FileDevNull,str,aflags,PFSWRITE,OSSWRITE,oss);             //CO 171024 - new sym framework
	//SYM::CalculateFactorGroup(FileDevNull,str,aflags,PFSWRITE,OSSWRITE,oss);                //CO 171024 - new sym framework
	//str.sgroup_radius=1.05*RadiusSphereLattice(lattice);                                    //CO 171024 - new sym framework
	//SYM::CalculateSpaceGroup(FileDevNull,str,aflags,FALSE,OSSWRITE,oss);                    //CO 171024 - new sym framework
	//   oss << banner << endl; // ----------------------------------------------------------------
	oss << surface::PrintNNdists(num_types,num_types_combinations,bbfrac,nndist,nndists);
	oss << "hklmax=" << hklmax << endl;
	oss << "bbdistance=" << bbdistance/nndist << "  surface::GetNNeighbours=" << surface::GetNNeighbours(str) << " real_bbdistance=" << bbdistance << endl;
	oss << "step=" << step << endl;
	oss << "SCANNING TRIVIAL PLANES" << endl;
	oss << surface::PrintHKLSigmaBB(num_types,num_types_combinations,bbfrac,bbdistance,bbdistances);
	// FFF
	FFF << banner << endl; // ----------------------------------------------------------------
	FFF << aflow::Banner("BANNER_TINY") << endl;
	FFF << banner << endl; // ----------------------------------------------------------------
	FFF << "HKL CALCULATION" << endl;
	if(search_trivial)  FFF << "TRIVIAL SEARCH" << endl;
	if(search_simple)   FFF << "SIMPLE SEARCH" << endl;
	if(search_complete) FFF << "COMPLETE SEARCH" << endl;
	FFF << surface::PrintNNdists(num_types,num_types_combinations,bbfrac,nndist,nndists);
	FFF << "hklmax=" << hklmax << endl;
	FFF << "bbdistance=" << bbdistance/nndist << "  surface::GetNNeighbours=" << surface::GetNNeighbours(str) << " real_bbdistance=" << bbdistance << endl;
	FFF << "step=" << step << endl;
	FFF << "SCANNING TRIVIAL PLANES" << endl;
	FFF << surface::PrintHKLSigmaBB(num_types,num_types_combinations,bbfrac,bbdistance,bbdistances);
	//
	// THIS ROUTINE IS VERY SLOW AND SHOULD BE MADE MULTI-THREADS
	vector<xvector<double> > vhkl; 
	//      for(double h=-dims_hklmax[1];h<=dims_hklmax[1];h+=step)
	//        for(double k=-dims_hklmax[2];k<=dims_hklmax[2];k+=step)
	//          for(double l=-dims_hklmax[3];l<=dims_hklmax[3];l+=step) 

	for(double h=dims_hklmax[1];h>=-dims_hklmax[1];h-=step)
	  for(double k=dims_hklmax[2];k>=-dims_hklmax[2];k-=step)
	    for(double l=dims_hklmax[3];l>=-dims_hklmax[3];l-=step)
	      {
		if(abs(h)<eps) h=0.0;
		if(abs(k)<eps) k=0.0;
		if(abs(l)<eps) l=0.0;
		// cerr << h << " " << k << " " << l << " " << endl;
		// hkl_found=surface::AddPlaneHKL(str, h, k, l,roughness,bbdistance,bbdistances,eps,jossmax,plane,planesreducible,TRUE,oss,FFFflag,FFF);
		vhkl.push_back(xvector<double>(3));
		vhkl.back()(1)=h;vhkl.back()(2)=k;vhkl.back()(3)=l;
		// cerr << vhkl.size() << endl;
	      }
	for(uint i=0;i<vhkl.size();i++) 
	  hkl_found=surface::AddPlaneHKL(str,vhkl.at(i),roughness,bbdistance,bbdistances,eps,jossmax,plane,planesreducible,TRUE,oss,FFFflag,FFF);
      }
    }
  
    if(LDEBUG) cerr << "surface::GetSurfaceHKLSearch: [3] " <<   endl;

    if(search_simple || search_complete) {
      if(mode==0 || mode==1 || mode==2 || mode==3 || mode==4) { // seek
	if(LDEBUG) cerr << "surface::GetSurfaceHKLSearch: [3] mode=" << mode <<  endl;
	if(search_simple)   {step=1.0;};
	if(search_complete) {step=1/3.0/4.0;};
	if(mode==0) {rrfrac=_RRFRAC_;bbfrac=_BBFRAC_;hklmax=-1;}//step=step;}
	if(mode==1) {rrfrac=abs(iparams(1));bbfrac=_BBFRAC_;hklmax=-1;}//step=step;}
	if(mode==2) {rrfrac=abs(iparams(1));bbfrac=abs(iparams(2));hklmax=-1;}//step=step;}
	if(mode==3) {rrfrac=abs(iparams(1));bbfrac=abs(iparams(2));hklmax=abs(iparams(3));}//step=step;}
	if(mode==4) {rrfrac=abs(iparams(1));bbfrac=abs(iparams(2));hklmax=abs(iparams(3));step=abs(iparams(4));}
	radius=rrfrac;
	bbdistance=bbfrac*nndist;
	bbdistances=bbfrac*nndists;
	radius=radius*modulus(a1+a2+a3)+bbdistance;
	dims=LatticeDimensionSphere(lattice,radius);
	if(search_simple)  {imin=0;imax=dims(1);jmin=0;jmax=dims(2);kmin=0;kmax=dims(3);};
	if(search_complete) {imin=-dims(1);imax=dims(1);jmin=-dims(2);jmax=dims(2);kmin=-dims(3);kmax=dims(3);};
	if(hklmax<0) {
	  dims_hklmax(1)=(double) dims(1);dims_hklmax(2)=(double) dims(2);dims_hklmax(3)=(double) dims(3);
	} else {
	  dims_hklmax(1)=hklmax;dims_hklmax(2)=hklmax;dims_hklmax(3)=hklmax;
	}
      
	oss << banner << endl; // ----------------------------------------------------------------
	oss << aflow::Banner("BANNER_TINY") << endl;
	oss << banner << endl; // ----------------------------------------------------------------
	if(search_simple)   oss << "SIMPLE SEARCH" << endl;
	if(search_trivial)  oss << "TRIVIAL SEARCH" << endl;
	if(search_complete) oss << "COMPLETE SEARCH" << endl;
	str.LatticeReduction_avoid=TRUE;
	str.sgroup_radius=1.05*RadiusSphereLattice(lattice);
        _kflags kflags; pflow::defaultKFlags4SymCalc(kflags,true);                                //CO 171024 - new sym framework
        pflow::defaultKFlags4SymWrite(kflags,PFSWRITE); kflags.KBIN_SYMMETRY_SGROUP_WRITE=false;  //CO 171024 - new sym framework
        pflow::PerformFullSymmetry(str,FileDevNull,aflags,kflags,OSSWRITE,oss);                   //CO 171024 - new sym framework
	//SYM::CalculatePointGroup(FileDevNull,str,aflags,PFSWRITE,OSSWRITE,oss);                 //CO 171024 - new sym framework
	//SYM::CalculateSitePointGroup(FileDevNull,str,aflags,PFSWRITE,OSSWRITE,oss);             //CO 171024 - new sym framework
	//SYM::CalculateFactorGroup(FileDevNull,str,aflags,PFSWRITE,OSSWRITE,oss);                //CO 171024 - new sym framework
	//str.sgroup_radius=1.05*RadiusSphereLattice(lattice);                                    //CO 171024 - new sym framework
	//SYM::CalculateSpaceGroup(FileDevNull,str,aflags,FALSE,OSSWRITE,oss);                    //CO 171024 - new sym framework
	// oss
	oss << banner << endl; // ----------------------------------------------------------------
	oss << "HKL CALCULATION" << endl;
	if(search_trivial)  oss << "TRIVIAL SEARCH" << endl;
	if(search_simple)   oss << "SIMPLE SEARCH" << endl;
	if(search_complete) oss << "COMPLETE SEARCH" << endl;
	oss << surface::PrintNNdists(num_types,num_types_combinations,bbfrac,nndist,nndists);
	oss << "CUTOFF rrfrac=" << rrfrac << endl;
	oss << "BOND   bbfrac=" << bbfrac << endl;
	oss << "HKLMAX hklmax=" << hklmax << endl;
	oss << "STEP   step=  " << step << endl;
	oss << " radius=" << radius << endl;
	oss << " dims=(" << dims(1) << "," << dims(2) << "," << dims(3) << ")" << endl;
	oss << " dims_hklmax=(" << dims_hklmax[1] << "," << dims_hklmax[2] << "," << dims_hklmax[3] << ")" << endl;
	oss << "SCANNING " << endl;
	// FFF
	FFF << banner << endl; // ----------------------------------------------------------------
	FFF << aflow::Banner("BANNER_TINY") << endl;
	FFF << banner << endl; // ----------------------------------------------------------------
	FFF << "HKL CALCULATION" << endl;
	if(search_trivial)  FFF << "TRIVIAL SEARCH" << endl;
	if(search_simple)   FFF << "SIMPLE SEARCH" << endl;
	if(search_complete) FFF << "COMPLETE SEARCH" << endl;
	FFF << surface::PrintNNdists(num_types,num_types_combinations,bbfrac,nndist,nndists);
	FFF << "CUTOFF rrfrac=" << rrfrac << endl;
	FFF << "BOND   bbfrac=" << bbfrac << endl;
	FFF << "HKLMAX hklmax=" << hklmax << endl;
	FFF << "STEP   step=  " << step << endl;
	FFF << " radius=" << radius << endl;
	FFF << " dims=(" << dims(1) << "," << dims(2) << "," << dims(3) << ")" << endl;
	FFF << " dims_hklmax=(" << dims_hklmax[1] << "," << dims_hklmax[2] << "," << dims_hklmax[3] << ")" << endl;
	FFF << "SCANNING " << endl;
	// ----
      
	for(int i=imin;i<=imax;i++) {
	  for(int j=jmin;j<=jmax;j++) {
	    for(int k=kmin;k<=kmax;k++) {
	      for(uint iat=0;iat<str.atoms.size();iat++) {
		// if(search_complete || (search_simple && iat==0))
		{
		  rrr=((double)i)*a1+((double)j)*a2+((double)k)*a3+str.atoms.at(iat).cpos;
		  // ddd(1)=(double) i;ddd(2)=(double) j;ddd(3)=(double) k;ddd=ddd+str.atoms.at(iat).fpos;  //USELESS
		  if(modulus(rrr)<=radius && modulus(rrr)>eps) {
		    grid_atoms_cpos_ptr = new xvector<double>(3);
		    *grid_atoms_cpos_ptr=rrr;                    
		    grid_atoms_cpos.push_back(grid_atoms_cpos_ptr);  
		    // grid_atoms_fpos_ptr = new xvector<double>(3);  / /USELESS
		    // *grid_atoms_fpos_ptr=ddd;                       //USELESS
		    // grid_atoms_fpos.push_back(grid_atoms_fpos_ptr); //USELESS
		    grid_atoms_number.push_back(str.atoms.at(iat).number);
		  }
		}
	      }
	    }
	  }
	}
	grid_atoms_size=grid_atoms_cpos.size();
	oss << "grid_atoms_size=" << grid_atoms_size << endl;
	FFF << "grid_atoms_size=" << grid_atoms_size << endl;
	oss << banner << endl; // ----------------------------------------------------------------
	FFF << banner << endl; // ----------------------------------------------------------------
	oss << surface::PrintHKLSigmaBB(num_types,num_types_combinations,bbfrac,bbdistance,bbdistances);
	FFF << surface::PrintHKLSigmaBB(num_types,num_types_combinations,bbfrac,bbdistance,bbdistances);
	//
	// extra juice !
	// if(0) if(search_complete)
	{ // only if complete
	  oss << "SCANNING TRIVIAL PLANES" << endl;
	  FFF << "SCANNING TRIVIAL PLANES" << endl;
	  for(double h=-dims_hklmax[1];h<=dims_hklmax[1];h+=step)
	    for(double k=-dims_hklmax[2];k<=dims_hklmax[2];k+=step)
	      for(double l=-dims_hklmax[3];l<=dims_hklmax[3];l+=step) {
		// for(int i=-5;i<=5;i++) {if(abs(h-i)<eps) h=(double) i;if(abs(k-i)<eps) k=(double) i;if(abs(l-i)<eps) l=(double) i;}
		if(abs(h)<eps) h=0.0;
		if(abs(k)<eps) k=0.0;
		if(abs(l)<eps) l=0.0;
		hkl_found=surface::AddPlaneHKL(str, h, k, l,roughness,bbdistance,bbdistances,eps,jossmax,plane,planesreducible,TRUE,oss,FFFflag,FFF);
	      }
	}
	oss << "SCANNING COMPLICATE PLANES" << endl;
	FFF << "SCANNING COMPLICATE PLANES" << endl;
	for(uint iat1=0;iat1<grid_atoms_size;iat1++) {
	  number1=grid_atoms_number[iat1];
	  rrr1=*grid_atoms_cpos[iat1];
	  // ddd1=*grid_atoms_fpos[iat1];
	  for(uint iat2=iat1+1;iat2<grid_atoms_size;iat2++) {
	    if(iat2!=iat1) {
	      number2=grid_atoms_number[iat2];
	      if(number2==number1) {
		rrr2=*grid_atoms_cpos[iat2];
		// ddd2=*grid_atoms_fpos[iat2];
		for(uint iat3=iat2+1;iat3<grid_atoms_size;iat3++)
		  if(iat3!=iat1 && iat3!=iat2) {
		    number3=grid_atoms_number[iat3];
		    if(number3==number2 && number2==number1) {
		      // all threee numbers identical
		      rrr3=*grid_atoms_cpos[iat3];
		      // ddd3=*grid_atoms_fpos[iat3];
		      determinant=det(rrr1,rrr2,rrr3);   // is zero if origin goes through the triangle (better avoid)
		      if(abs(determinant)>eps) {
			area=surface::TriangleArea(rrr1,rrr2,rrr3);
			if(area>eps) {
			  PlaneGetABCD(a,b,c,d,rrr1,rrr2,rrr3);//  abs(d/sqrt(a*a+b*b+c*c));
			  if(abs(d/sqrt(a*a+b*b+c*c))>eps) {
			    hkl=PlaneGetHKL(rrr1,rrr2,rrr3,a1,a2,a3);
			    double h=hkl(1),k=hkl(2),l=hkl(3);
			    if((abs(h)>=1.0 || abs(h)<eps) &&
			       (abs(k)>=1.0 || abs(k)<eps) &&
			       (abs(l)>=1.0 || abs(l)<eps)) {
			      // for(int i=-5;i<=5;i++) {if(abs(h-i)<eps) h=(double) i;if(abs(k-i)<eps) k=(double) i;if(abs(l-i)<eps) l=(double) i;}
			      if(abs(h)<eps) h=0.0;
			      if(abs(k)<eps) k=0.0;
			      if(abs(l)<eps) l=0.0;
			      // cerr << h << " " << k << " " << l << endl;
			      if(search_simple)   hkl_found=surface::AddPlaneHKL(str, h, k, l,roughness,bbdistance,bbdistances,eps,jossmax,plane,planesreducible,TRUE,oss,FFFflag,FFF);
			      if(search_complete) hkl_found=surface::AddPlaneHKL(str, h, k, l,roughness,bbdistance,bbdistances,eps,jossmax,plane,planesreducible,TRUE,oss,FFFflag,FFF);
			      if(search_complete) hkl_found=surface::AddPlaneHKL(str,-h, k, l,roughness,bbdistance,bbdistances,eps,jossmax,plane,planesreducible,TRUE,oss,FFFflag,FFF);
			      if(search_complete) hkl_found=surface::AddPlaneHKL(str, h,-k, l,roughness,bbdistance,bbdistances,eps,jossmax,plane,planesreducible,TRUE,oss,FFFflag,FFF);
			      if(search_complete) hkl_found=surface::AddPlaneHKL(str, h, k,-l,roughness,bbdistance,bbdistances,eps,jossmax,plane,planesreducible,TRUE,oss,FFFflag,FFF);
			      if(search_complete) hkl_found=surface::AddPlaneHKL(str, h,-k,-l,roughness,bbdistance,bbdistances,eps,jossmax,plane,planesreducible,TRUE,oss,FFFflag,FFF);
			      if(search_complete) hkl_found=surface::AddPlaneHKL(str,-h, k,-l,roughness,bbdistance,bbdistances,eps,jossmax,plane,planesreducible,TRUE,oss,FFFflag,FFF);
			      if(search_complete) hkl_found=surface::AddPlaneHKL(str,-h,-k, l,roughness,bbdistance,bbdistances,eps,jossmax,plane,planesreducible,TRUE,oss,FFFflag,FFF);
			    }
			  }
			}
		      }
		    }
		  }
	      }
	    }
	  }
	}
	//   oss << endl;
      }
    }
    if(LDEBUG) cerr << "surface::GetSurfaceHKLSearch: [4] " <<   endl;

    if(hkl_found) {;} // dummy load

    oss << banner << endl; // ----------------------------------------------------------------
    FFF << banner << endl; // ----------------------------------------------------------------
    oss << "REDUCIBLE: " << planesreducible.size() << endl;
    FFF << "REDUCIBLE: " << planesreducible.size() << endl;
    oss << banner << endl; // ----------------------------------------------------------------
    FFF << banner << endl; // ----------------------------------------------------------------
    // reduce --------------------------------------------------------------------
    surface::ReducePrintSurfaces(str,eps,planesreducible,planesirreducible,planesirreducible_images,TRUE,oss,FFFflag,FFF);
    // print end  ----------------------------------------------------------------
    oss << banner << endl; // ----------------------------------------------------------------
    FFF << banner << endl; // ----------------------------------------------------------------
    oss << "IRREDUCIBLE= " << planesirreducible.size() << "  REDUCIBLE= " << planesreducible.size() << endl;
    FFF << "IRREDUCIBLE= " << planesirreducible.size() << "  REDUCIBLE= " << planesreducible.size() << endl;
    oss << banner << endl; // ----------------------------------------------------------------
    FFF << banner << endl; // ----------------------------------------------------------------
    oss << aflow::Banner("BANNER_TINY") << endl;
    FFF << aflow::Banner("BANNER_TINY") << endl;
    oss << banner << endl; // ----------------------------------------------------------------
    FFF << banner << endl; // ----------------------------------------------------------------
    // END CLEAN EVERYTHING ------------------------------------------------------------------
    FFF.close();
    for(uint i=0;i<grid_atoms_cpos.size();i++) 
      delete grid_atoms_cpos[i];
    grid_atoms_cpos.clear();

    if(LDEBUG) cerr << "surface::GetSurfaceHKLSearch: END " <<   endl;

    return Krun;
  }
} // namespace surface

// ********************************************************************************************************/
// ********************************************************************************************************/
// RC-09: Depending on the initial input POSCAR, the obtained 2 basis vectors in the slab layer may be not as minimum as possible.
//        It should be corrected afterwards, by applying the general methods of unit cell minimization to the obtained here final slab POSCAR.
// RC-09: The slab uc basis vectors coordinates are presented wrt initial POSCAR Cartesian system. So it is easy to find the position of slab uc basis vectors wrt initial POSCAR sites.
// RC-09: If (a) parent latice is undistorted cubic-like and (b) Cartesian axises (wrt which the uc basis vectors of initial POSCAR are determined) are directed along cubic edges
//        then the Cartesian coordinates of third (normal to layers) basis bector of slab uc (S3) determine the slab Miller indices wrt parent CUBIC unit cell.
//        E.g. for L10 with 4at/uc and 2at/us POSCAR definitions: (111)of4/uc=(101)of2/uc, (010)of4/uc=(110)of2/uc, (001)of4/uc=(001)of2/uc, (110)of4/uc=(100)of2/uc

#define _slab_file_OUT        string("POSCAR_slab_hkl")
#define _slab_file2_OUT       string("Plane.dat")
#define _slab_epsilon_        double(1e-4)        // Small nonzero number

// ****************************************************************************************************************/
namespace slab {
  double VectorAbsValue(int Layer, int NinLayer /*IN*/,const xmatrix<double>& UnitCellVector,const vector<vector<vector<int> > >& LayerSitesDirCoords) {  
    int i1, j1;
    double AbsValue=0;
    for(j1=1;j1<=3;j1++) {  
      double x=0;
      for(i1=1;i1<=3;i1++) {
	x+=LayerSitesDirCoords[Layer][NinLayer][i1] * UnitCellVector[i1][j1];}
      AbsValue+=x*x;
    }
    AbsValue=sqrt(AbsValue);
    return AbsValue;
  }
} // namespace slab

// ****************************************************************************************************************/
namespace slab {
  double VectorScalarMult(int Layer1, int NinLayer1, int Layer2, int NinLayer2 /*IN*/,const xmatrix<double>& UnitCellVector,const vector<vector<vector<int> > >& LayerSitesDirCoords) {
    double Help1,Help2;
    int i1,j1;
    double ScalarMult=0;
    for(j1=1;j1<=3;j1++) {  
      Help1=0; 
      for(i1=1;i1<=3;i1++) {
	Help1+=LayerSitesDirCoords[Layer1][NinLayer1][i1] * UnitCellVector[i1][j1];}
      Help2=0;  
      for(i1=1;i1<=3;i1++) {
	Help2+=LayerSitesDirCoords[Layer2][NinLayer2][i1] * UnitCellVector[i1][j1];}
      ScalarMult+=Help1*Help2;
    }
    return ScalarMult;
  }
} // namespace slab
  
// ****************************************************************************************************************/
namespace slab {
  double CosAngle(int Layer1, int NinLayer1, int Layer2, int NinLayer2 /*IN*/,const xmatrix<double>& UnitCellVector,const vector<vector<vector<int> > >& LayerSitesDirCoords) {
    double CosAngleOUT,ScalarMult, AbsValue1, AbsValue2;
    ScalarMult=slab::VectorScalarMult(Layer1,NinLayer1,Layer2,NinLayer2,UnitCellVector,LayerSitesDirCoords);
    AbsValue1=slab::VectorAbsValue(Layer1,NinLayer1,UnitCellVector,LayerSitesDirCoords);
    AbsValue2=slab::VectorAbsValue(Layer2,NinLayer2,UnitCellVector,LayerSitesDirCoords);
    CosAngleOUT=ScalarMult/(AbsValue1*AbsValue2);
    return CosAngleOUT;
  }
} // namespace slab

// ****************************************************************************************************************/
namespace slab {
  double hkl_CartCoord_Length(xvector<double>& hkl_CartCoord,const xmatrix<double>& UnitCellVector,const xvector<double>& hkl) {
    // ReciprUnitCellVector,hkl_CartCoord_Length are in 2PI/LattPar[0] (VolumeReciprUC in LattPar[0]^3) units
    int i,j; double VolumeReciprUC;
    xmatrix<double> ReciprUnitCellVector(3,3);
    double hkl_Length;
    VolumeReciprUC=UnitCellVector(1,1)*(UnitCellVector(2,2)*UnitCellVector(3,3)-UnitCellVector(3,2)*UnitCellVector(2,3))-
      UnitCellVector(1,2)*(UnitCellVector(2,1)*UnitCellVector(3,3)-UnitCellVector(3,1)*UnitCellVector(2,3))+
      UnitCellVector(1,3)*(UnitCellVector(2,1)*UnitCellVector(3,2)-UnitCellVector(3,1)*UnitCellVector(2,2));
    ReciprUnitCellVector(1,1)= UnitCellVector(2,2)*UnitCellVector(3,3)-UnitCellVector(3,2)*UnitCellVector(2,3);
    ReciprUnitCellVector(1,2)=-(UnitCellVector(2,1)*UnitCellVector(3,3)-UnitCellVector(3,1)*UnitCellVector(2,3));
    ReciprUnitCellVector(1,3)= UnitCellVector(2,1)*UnitCellVector(3,2)-UnitCellVector(3,1)*UnitCellVector(2,2);
    
    ReciprUnitCellVector(2,1)= UnitCellVector(3,2)*UnitCellVector(1,3)-UnitCellVector(1,2)*UnitCellVector(3,3);
    ReciprUnitCellVector(2,2)=-(UnitCellVector(3,1)*UnitCellVector(1,3)-UnitCellVector(1,1)*UnitCellVector(3,3));
    ReciprUnitCellVector(2,3)= UnitCellVector(3,1)*UnitCellVector(1,2)-UnitCellVector(1,1)*UnitCellVector(3,2);
    
    ReciprUnitCellVector(3,1)= UnitCellVector(1,2)*UnitCellVector(2,3)-UnitCellVector(2,2)*UnitCellVector(1,3);
    ReciprUnitCellVector(3,2)=-(UnitCellVector(1,1)*UnitCellVector(2,3)-UnitCellVector(2,1)*UnitCellVector(1,3));
    ReciprUnitCellVector(3,3)= UnitCellVector(1,1)*UnitCellVector(2,2)-UnitCellVector(2,1)*UnitCellVector(1,2);
    for(i=1;i<=3;i++) {
      for(j=1;j<=3;j++) {
	ReciprUnitCellVector(i,j)/=VolumeReciprUC;
      }
    }
    // cerr << "Reciprocal Basis: "; for(i=1;i<=3;i++) { cerr << "("; for(j=1;j<=3;j++) {cerr << ReciprUnitCellVector(i,j) << ","; } cerr << ")"; }; cerr << endl << endl;
    for(j=1;j<=3;j++) {
      hkl_CartCoord[j]=0;
      for(i=1;i<=3;i++) {
	hkl_CartCoord[j]+=hkl(i)*ReciprUnitCellVector(i,j);
      } 
    }
    hkl_Length=0;
    for(j=1;j<=3;j++) {
      hkl_Length+=hkl_CartCoord[j]*hkl_CartCoord[j];
    }
    hkl_Length=sqrt(hkl_Length);
    return hkl_Length;
  } // END Reciprocal_Lattice_Calc //
} // namespace slab
  
// ********************************************************************************************************/
namespace slab {
  xstructure MAKE_SLAB(string options,istream& cin) {
    xstructure str_in(cin,IOAFLOW_AUTO);
    xstructure str_out("");
    str_out=MAKE_SLAB(options,str_in);
    //MAKE_SLAB(options,str_in);
    //  cout << str_out;
    return str_out;
  }  
} // namespace slab

namespace slab {
  xstructure MAKE_SLAB(string options,xstructure& str_in) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "slab::MAKE_SLAB: BEGIN" << endl;
    vector<string> tokens;
    aurostd::string2tokens(options,tokens,",");
    if(tokens.size()<3 || tokens.size()>5) {
      init::ErrorOption(cout,options,"slab::MAKE_SLAB","aflow --slab=h,k,l[,#filled_layers[,#vacuum layers]] < POSCAR");
      exit(0);
    }
    int i=0,j=0,k=0;
    //  options[0-2]="h" "k" "l" options[3-4]="NumFilledLayers" "NumEmptyLayers"
 

    xvector<double> hkl(3);
    hkl(1)=1;hkl(2)=1;hkl(3)=1; 
    if(tokens.size()>=1) hkl(1)=aurostd::string2utype<double>(tokens.at(0));
    if(tokens.size()>=2) hkl(2)=aurostd::string2utype<double>(tokens.at(1));
    if(tokens.size()>=3) hkl(3)=aurostd::string2utype<double>(tokens.at(2));

    string file_OUT=_slab_file_OUT;
    stringstream temp;
    temp << "POSCAR_slab_"<< hkl(1)<<hkl(2)<<hkl(3);
    temp >> file_OUT;
    
    int  NumFilledLayers=2, NumEmptyLayers=3;  // number of filled/empty layers in slab    
    int const SearchMax=20+ NumFilledLayers + NumEmptyLayers;  
    if(tokens.size()>=4) NumFilledLayers=aurostd::string2utype<int>(tokens.at(3));
    if(tokens.size()>=5) NumEmptyLayers=aurostd::string2utype<int>(tokens.at(4));

    int const In_Plane_Multiplication[3]={0,1,1};    // In_Plane_Multiplication^2 = number of plane unit cells in slab unit cell
    vector<int> NumAtomsForElementUC(2,0);
    vector< vector< vector<double> > > AtomDirCoords;
    vector< vector< vector<int> > > LayerSitesDirCoords;
    

    //   int const SearchMax=20+ NumFilledLayers + NumEmptyLayers;  
    // [-SearchMax to SearchMax]^3 number of tested sites in lattice=12 + NumFilledLayers + NumEmptyLayers
    
    int NumSitesInPlane[NumFilledLayers+1], LayerBasis[3], Layer0Basis[3], t1=0, t2=0, VECTOR[3], NumLayers, NumBasisSitesInSlab,Layer, l[4], Sum_Int, Element, ElementSite;
    double DET=0.0, s[3], fractpart[3], intpart, Layer0BasisCosAngle=1, Layer0BasisCosAngleCurrent, LastCandidateCoord3, CurrCandidateCoord3, hkl_Length;
    xvector<double> SiteCartCoord(3), RS(3), AbsValueSlab_Basis(3), SiteDirectCoordWRTslab(3), hkl_CartCoord(3), SiteDirectCoord(3);
    xmatrix<double> MATRIX(3,3),INVERSE_MATRIX(3,3);
    xmatrix<double> Slab_Basis_CartCoord(3,3);

    bool Continue, BasisFound;

    int NumberElements=str_in.num_each_type.size();
    string Title=str_in.title;
    double LattPar=str_in.scale;
    
    xmatrix<double> UnitCellVector(3,3);
    
    for(i=1;i<4;i++) {
      for(j=1;j<4;j++) {
	UnitCellVector(i,j)=str_in.lattice(i,j);  
      }
    }

    for(i=1;i<NumberElements+1;i++) {
      NumAtomsForElementUC[i]=str_in.num_each_type[i-1];
    }

    //i=1;NumberElements=0;
    //while (iss >> NumAtomsForElementUC[i]) {i++; NumberElements++; NumAtomsForElementUC.push_back(0);};
    //here
  
    int iatom=-1;
    AtomDirCoords.resize(NumberElements+1);
    for(k=1;k<=NumberElements;k++) {
      AtomDirCoords[k].resize(NumAtomsForElementUC[k]+1);
    }
    for(k=1;k<=NumberElements;k++) {
      for(i=1;i<=NumAtomsForElementUC[k];i++) {
	AtomDirCoords[k][i].resize(4); }
    }
    for(k=1;k<=NumberElements;k++) {
      for(i=1;i<=NumAtomsForElementUC[k];i++) {
	iatom++; 
	for(j=1;j<=3;j++) {
	  AtomDirCoords[k][i][j]=str_in.atoms.at(iatom).fpos[j];
	}
      }
    }

    //  int itmp=0;
    //  AtomDirCoords.resize(NumberElements+1);
    //  for(i=1;i<=NumberElements+1;i++) {
    // itmp=NumAtomsForElementUC[i];
    // AtomDirCoords[i].resize(itmp);
    // for(j=1;j<=itmp;j++)
    //  AtomDirCoords[i][j].resize(4);
    //  }
    //  int iatom=0;
    //  for(i=1;i<=NumberElements;i++) {
    // for(j=1;j<=NumAtomsForElementUC[i];j++) {
    //  iatom++;
    //  for(k=1;k<=3;k++) {
    // AtomDirCoords[i][j][k]=str_in.atoms.at(iatom).fpos[k];
    //  }} }

    vector<int> NumSites(NumberElements+1,0);

    //---------------- Sorting (by layer number) of all sites in [-SearchMax to SearchMax]^3 box -----------------------------------------------//  
    //  from equation Sum_i (l_i*n_i)=Layer number; {l_i}=sites coordinates wrt direct-lattice basis; {n_i}=hkl

    LayerSitesDirCoords.resize(NumFilledLayers+1);  

    for(Layer=0; Layer<=NumFilledLayers; Layer++) {
      NumSitesInPlane[Layer]=0;}
    
    for(l[1]=SearchMax;l[1]>=-SearchMax;l[1]--) {
      for(l[2]=SearchMax;l[2]>=-SearchMax;l[2]--) {
	for(l[3]=SearchMax;l[3]>=-SearchMax;l[3]--) {
	  Sum_Int=l[1]*hkl(1)+l[2]*hkl(2)+l[3]*hkl(3);
	  if((Sum_Int > -1) && (Sum_Int < NumFilledLayers+1)) {
	    NumSitesInPlane[Sum_Int]++;  
	    LayerSitesDirCoords[Sum_Int].resize(NumSitesInPlane[Sum_Int]+1);
	    LayerSitesDirCoords[Sum_Int][NumSitesInPlane[Sum_Int]].resize(4);
	    for(i=1;i<=3;i++) {
	      LayerSitesDirCoords[Sum_Int][NumSitesInPlane[Sum_Int]][i]=l[i];
	    }  
	  }
	}
      }
    }
    if(LDEBUG) {
      for(Layer=0; Layer<=0 ; Layer++) {
	cerr << "Layer= " << Layer << endl;
	for(j=1;j<=NumSitesInPlane[Layer];j++) {
	  cerr << j << " " 
	       << LayerSitesDirCoords [Layer][j][1] << " "
	       << LayerSitesDirCoords [Layer][j][2] << " "
	       << LayerSitesDirCoords[Layer][j][3] << endl;
	}
      }
    }

    stringstream f2out; 
    double x;
    for(k=1;k<=NumSitesInPlane[0];k++) {
      for(j=1;j<=3;j++) {
	x=0.0;  
	for(i=1;i<=3;i++) {
	  x+=LayerSitesDirCoords[0][k][i]*UnitCellVector(i,j);
	}
	f2out << x << " ";
      }
      f2out << k << endl;
    }
    aurostd::stringstream2file(f2out,_slab_file2_OUT);  
    //return 1;
    //cin >> Help_string;
    //----------------  Searching for basis vectors in m=0 layer ------------------------------------------------------//
    //  Condition: any layer-site-vector (k-numbered) is a integer-linear superposition of two layer basis vectors (LayerBasis[1,2]-numbered in a search)
    //  The angle between basis vectors in layer is seek to be closest to 90 degree

    //cerr << "Searching for basis vectors in m=0 layer with closest to 90 degree inter-angle" << endl;

    BasisFound=false;
    LayerBasis[1]=1; 
    while (LayerBasis[1]<=NumSitesInPlane[0]) {//cerr << LayerBasis[1] << " / " << NumSitesInPlane[0] << endl;
      LayerBasis[2]=LayerBasis[1]+1;
      while (LayerBasis[2]<=NumSitesInPlane[0]) {
	for(i=1;i<=2;i++) {
	  for(j=1;j<=3;j++) {
	    MATRIX[j][i]=LayerSitesDirCoords[0][LayerBasis[i]][j];
	  }
	}
	Continue=true;
	i=1; 
	while (i<=2 && Continue==true) {
	  j=i+1;
	  while (j<=3 && Continue==true) {  
	    DET=MATRIX[i][1]*MATRIX[j][2]-MATRIX[i][2]*MATRIX[j][1];  
	    if(fabs(DET)>_slab_epsilon_) {
	      t1=i; 
	      t2=j;
	      Continue=false;
	    }
	    j++;
	  }
	  i++;
	}
	/*
	  if(LayerBasis[1]==48 && LayerBasis[2]==49)
	  {  cerr << " DET= " << DET << endl;
	  for(i=1;i<=2;i++) {for(j=1;j<=3;j++) { cerr << LayerSitesDirCoords[0][LayerBasis[i]][j] << " ";} cerr << endl;}
	  cerr << MATRIX[1][1] << " " << MATRIX[1][2] << endl << MATRIX[2][1] << " " << MATRIX[2][2] << endl;
	  }*/
  
	if(fabs(DET)>_slab_epsilon_) {
	  INVERSE_MATRIX[1][1]= MATRIX[t2][2]/DET;  INVERSE_MATRIX[1][2]=-MATRIX[t1][2]/DET;
	  INVERSE_MATRIX[2][1]=-MATRIX[t2][1]/DET;  INVERSE_MATRIX[2][2]= MATRIX[t1][1]/DET;
	  VECTOR[1]=t1; VECTOR[2]=t2;
	  Continue=true;
	  k=1; 
	  while (k<=NumSitesInPlane[0] && Continue==true) {
	    for(i=1;i<=2;i++) {
	      s[i]=0;
	      for(j=1;j<=2;j++) {
		s[i]+=INVERSE_MATRIX(i,j)*LayerSitesDirCoords[0][k][VECTOR[j]];
	      }
	    }
	    for(i=1;i<=2;i++) {
	      fractpart[i]=modf (s[i] , &intpart);
	    }
  
	    //if(fabs(fractpart[1]-1.0)<_slab_epsilon_ || fabs(fractpart[2]-1.0)<_slab_epsilon_) {for(i=1;i<=2;i++) { cerr << s[i] << " (" << fractpart[i] << ") " << endl; };}
	    //cerr << s[1] << " (" << fractpart[1] << "); "  << s[2] << " (" << fractpart[2] << ") " << endl; cin >> Help_string;
	    //cerr << fabs(fractpart[1]) << "absFrac-_slab_epsilon_" << _slab_epsilon_ << endl;
	    //  if((fabs(fractpart[1])<_slab_epsilon_) &&  (fabs(fractpart[2])<_slab_epsilon_)) { cerr << LayerBasis[1] << " =LayerBasis= " << LayerBasis[2] << endl; cerr << s[1] << " " << s[2] << endl; /*cin >> Help_double;*/}
	    if((fabs(fractpart[1])>_slab_epsilon_ && fabs(fabs(fractpart[1])-1.0)>_slab_epsilon_) ||  (fabs(fractpart[2])>_slab_epsilon_ && fabs(fabs(fractpart[2])-1.0)>_slab_epsilon_)) {
	      Continue=false;
	      /*if(LayerBasis[1]==48 && LayerBasis[2]==49) {
		cerr << LayerBasis[1] << " =LayerBasis= " << LayerBasis[2] << endl;
		cerr << "WRONG-Basis:" << endl; for(i=1;i<=2;i++) {for(j=1;j<=3;j++) { cerr << LayerSitesDirCoords[0][LayerBasis[i]][j] << " ";} cerr << endl;}
		cerr << "WRONG-Plane-point:" << endl; for(j=1;j<=3;j++) {cerr << LayerSitesDirCoords[0][k][j] << " ";} ; cerr << endl;
		cerr << s[1] << " (" << fractpart[1] << "); "  << s[2] << " (" << fractpart[2] << ") " << endl;
		cin >> Help_string;}*/
	    }
	    k++;
	  }
	} else {
	  Continue=false;
	}
	if(Continue==true) { 
	  //cerr << "Basis:" << endl; for(i=1;i<=2;i++) {for(j=1;j<=3;j++) { cerr << LayerSitesDirCoords[0][LayerBasis[i]][j] << " ";} cerr << endl;}
	  Layer0BasisCosAngleCurrent=CosAngle(0,LayerBasis[1],0,LayerBasis[2], UnitCellVector, LayerSitesDirCoords);
	  //cerr << "CosAngle=" << Layer0BasisCosAngleCurrent << endl;
	  if(Layer0BasisCosAngleCurrent>=0 && Layer0BasisCosAngleCurrent<Layer0BasisCosAngle && fabs(Layer0BasisCosAngleCurrent-Layer0BasisCosAngle)>_slab_epsilon_) {
	    for(i=1;i<=2;i++) {
	      Layer0Basis[i]=LayerBasis[i]; 
	    }
	    Layer0BasisCosAngle=Layer0BasisCosAngleCurrent;  BasisFound=true;
	    //cerr << "Cos0Angle=" << Layer0BasisCosAngle << endl; //cin >> Help_double;
	  }
	}
	LayerBasis[2]++;
      }
      LayerBasis[1]++;
    }
    if(BasisFound==false || 180.0/PI*acos(Layer0BasisCosAngle)<10.0) {
      cerr << "slab::MAKE_SLAB: Basis in plane was not found" << endl; }//return 1;}
    for(i=1;i<=2;i++) {
      AbsValueSlab_Basis[i]=slab::VectorAbsValue(0,Layer0Basis[i],UnitCellVector,LayerSitesDirCoords);
    }
    //cerr << "Primitive Basis in plane:" << endl; for(i=1;i<=2;i++) {for(j=1;j<=3;j++) { cerr << LayerSitesDirCoords[0][Layer0Basis[i]][j] << " ";} cerr << endl;} 
    //cerr << "------------------------------" << endl;
    //----------------  Build "vertical" basis vector perpendicular to slab  ------------------------------------------------------//  
    hkl_Length=slab::hkl_CartCoord_Length(hkl_CartCoord,UnitCellVector,hkl);
    //cerr << "hkl_CartCoord="; for(i=1;i<=3;i++) { cerr << hkl_CartCoord[i] << " ";}; cerr << endl;
    //cerr << "hkl_Length=" << hkl_Length << endl;
    double Help_double=(NumFilledLayers + NumEmptyLayers)/(hkl_Length*hkl_Length);
    for(j=1;j<=3;j++) {
      Slab_Basis_CartCoord[3][j]=Help_double*hkl_CartCoord[j];
    }
    AbsValueSlab_Basis[3]=0;
    for(j=1;j<=3;j++) {
      AbsValueSlab_Basis[3]+=Slab_Basis_CartCoord[3][j]*Slab_Basis_CartCoord[3][j];
    };
    AbsValueSlab_Basis[3]=sqrt(AbsValueSlab_Basis[3]);

    //cerr << "AbsValueSlab_Basis[3]=" << AbsValueSlab_Basis[3] << endl;

    //----------------  Build two "horisontal" slab basis vectors  ------------------------------------------------------//  
    for(k=1;k<=2;k++) {
      for(j=1;j<=3;j++) {
	Slab_Basis_CartCoord(k,j)=0;  
	for(i=1;i<=3;i++) {
	  Slab_Basis_CartCoord(k,j)+=LayerSitesDirCoords[0][Layer0Basis[k]][i]*UnitCellVector(i,j);
	}
      }
    }
    for(k=1;k<=2;k++) {
      for(j=1;j<=3;j++) {
	Slab_Basis_CartCoord(k,j)*=In_Plane_Multiplication[k];
      }
    }
    //----------------  Making Slab_Basis to be right-handed: if 1[23]<0 then 1<->2 ------------------------------------------------------//  
    Help_double= Slab_Basis_CartCoord(1,1)*(Slab_Basis_CartCoord(2,2)*Slab_Basis_CartCoord(3,3)-Slab_Basis_CartCoord(3,2)*Slab_Basis_CartCoord(2,3))-
      Slab_Basis_CartCoord(1,2)*(Slab_Basis_CartCoord(2,1)*Slab_Basis_CartCoord(3,3)-Slab_Basis_CartCoord(3,1)*Slab_Basis_CartCoord(2,3))+
      Slab_Basis_CartCoord(1,3)*(Slab_Basis_CartCoord(2,1)*Slab_Basis_CartCoord(3,2)-Slab_Basis_CartCoord(3,1)*Slab_Basis_CartCoord(2,2));  
    if(Help_double<0) {
      for(j=1;j<=3;j++) {
	Help_double=Slab_Basis_CartCoord(1,j);
	Slab_Basis_CartCoord(1,j)=Slab_Basis_CartCoord(2,j);
	Slab_Basis_CartCoord(2,j)=Help_double;
      }
    }
    //for(i=1;i<=3;i++) {for(j=1;j<=3;j++) {cerr << Slab_Basis_CartCoord(i,j) << " ";} cerr << endl;} 
    //----------------  Searching for sites (in filled layers) that are within the slab unit cell  ------------------------------------------------------//  
    MATRIX(1,1)=0;
    for(j=1;j<=3;j++) {
      MATRIX(1,1)+=Slab_Basis_CartCoord(1,j)*Slab_Basis_CartCoord(1,j);
    }
    MATRIX(2,2)=0;
    for(j=1;j<=3;j++) {
      MATRIX(2,2)+=Slab_Basis_CartCoord(2,j)*Slab_Basis_CartCoord(2,j);
    }
    MATRIX(1,2)=0; 
    for(j=1;j<=3;j++) {
      MATRIX(1,2)+=Slab_Basis_CartCoord(1,j)*Slab_Basis_CartCoord(2,j);
    }
    MATRIX(2,1)=MATRIX(1,2);
    DET=MATRIX(1,1)*MATRIX(2,2)-MATRIX(1,2)*MATRIX(2,1);

    //cerr << MATRIX(1,1) << " " << MATRIX(1,2) << endl << MATRIX(2,1) << " " << MATRIX(2,2) << endl << " DET= " << DET << endl;
    INVERSE_MATRIX(1,1)= MATRIX(2,2)/DET;  INVERSE_MATRIX(1,2)=-MATRIX(1,2)/DET;
    INVERSE_MATRIX(2,1)=-MATRIX(2,1)/DET;  INVERSE_MATRIX(2,2)= MATRIX(1,1)/DET;
    vector< vector<double> >  BasisSitesInSlabDirectCoord;
    NumBasisSitesInSlab=0;
    for(Layer=0; Layer<=(NumFilledLayers-1); Layer++) {  
      for(k=1;k<=NumSitesInPlane[Layer];k++) {
	for(j=1;j<=3;j++) {
	  SiteCartCoord[j]=0;
	  for(i=1;i<=3;i++) {
	    SiteCartCoord[j]+=LayerSitesDirCoords[Layer][k][i]*UnitCellVector(i,j);
	  }
	}
	for(i=1;i<=2;i++) {
	  RS[i]=0;
	  for(j=1;j<=3;j++) {
	    RS[i]+=SiteCartCoord[j]*Slab_Basis_CartCoord(i,j);
	  }
	}
	for(i=1;i<=2;i++) {
	  SiteDirectCoord[i]=0;
	  for(j=1;j<=2;j++) {
	    SiteDirectCoord[i]+=INVERSE_MATRIX(i,j)*RS[j];
	  }
	}
	//cerr << SiteDirectCoord[1] << " " << SiteDirectCoord[2] << endl; // cin >> Help_double;
	if(SiteDirectCoord[1]>=0 && SiteDirectCoord[1]<1 && fabs(SiteDirectCoord[1]-1)>_slab_epsilon_ && SiteDirectCoord[2]>=0 && SiteDirectCoord[2]<1 && fabs(SiteDirectCoord[2]-1)>_slab_epsilon_) {
	  //SiteDirectCoord[3]=Layer/hkl_Length/AbsValueSlab_Basis[3];
	  SiteDirectCoord[3]=1.0*Layer/(NumFilledLayers + NumEmptyLayers);
	  // cerr << SiteDirectCoord[1] << " " << SiteDirectCoord[2] << " " << SiteDirectCoord[3] << endl;
	  NumBasisSitesInSlab=NumBasisSitesInSlab+1;
	  BasisSitesInSlabDirectCoord.resize(NumBasisSitesInSlab+1);  
	  BasisSitesInSlabDirectCoord[NumBasisSitesInSlab].resize(4);  
	  for(j=1;j<=3;j++) {
	    BasisSitesInSlabDirectCoord[NumBasisSitesInSlab][j]=SiteDirectCoord[j];
	  }  
	}  
      }
    }
    //---  Coordinates of initial slab unit cell sites (maybe with actual number of layers to be larger than necessary because only layers of one Bravais lattice were layer-enumerated) ------------//  
    vector< vector< vector<double> > > ListSiteDirectCoordWRTslab, ListSiteCartCoord;
    ListSiteDirectCoordWRTslab.resize(NumberElements+1);  
    ListSiteCartCoord.resize(NumberElements+1);
    for(Element=1; Element<=NumberElements; Element++) {
      for(ElementSite=1; ElementSite<=NumAtomsForElementUC[Element]; ElementSite++) {
	for(l[1]=SearchMax;l[1]>=-SearchMax;l[1]--) {
	  for(l[2]=SearchMax;l[2]>=-SearchMax;l[2]--) {
	    for(l[3]=SearchMax;l[3]>=-SearchMax;l[3]--) {
	      for(j=1;j<=3;j++) {
		SiteCartCoord[j]=0;
		for(i=1;i<=3;i++) {
		  SiteCartCoord[j]+=(AtomDirCoords[Element][ElementSite][i]+l[i])*UnitCellVector(i,j);
		}
	      }
	      for(i=1;i<=3;i++) {
		RS[i]=0; for(j=1;j<=3;j++) {
		  RS[i]+=SiteCartCoord[j]*Slab_Basis_CartCoord(i,j);
		} 
	      }
	      //cerr << Element << endl;
	      //for(i=1;i<=3;i++) { cerr << l[i] << " ";}; cerr << " l" << endl;
	      //for(i=1;i<=3;i++) { cerr << AtomDirCoords[Element][ElementSite][i]+l[i] << " ";}; cerr << " AtomDirCoords" << endl;
	      //for(i=1;i<=3;i++) { cerr << SiteCartCoord[i] << " ";}; cerr << " SiteCartCoord" << endl;
  
	      /*
		for(i=1;i<=3;i++) { cerr << RS[i] << " ";}; cerr << " RS" << endl;
		cerr << MATRIX[1][1] << " " << MATRIX[1][2] << endl;
		cerr << MATRIX[2][1] << " " << MATRIX[2][2] << endl;
		cerr  << " DET= " << DET << endl;
		cerr << INVERSE_MATRIX[1][1] << " " << INVERSE_MATRIX[1][2] << endl;
		cerr << INVERSE_MATRIX[2][1] << " " << INVERSE_MATRIX[2][2] << endl;
	      */

	      for(i=1;i<=2;i++) {
		SiteDirectCoordWRTslab[i]=0;
		for(j=1;j<=2;j++) {
		  SiteDirectCoordWRTslab[i]+=INVERSE_MATRIX(i,j)*RS[j]; 
		}
	      }
	      SiteDirectCoordWRTslab[3]=RS[3]/(AbsValueSlab_Basis[3]*AbsValueSlab_Basis[3]);
	      for(i=1;i<=3;i++) {
		if(fabs(SiteDirectCoordWRTslab[i]-0.0)<_slab_epsilon_) {
		  SiteDirectCoordWRTslab[i]=0.0;
		} 
	      }
	      for(i=1;i<=3;i++) {
		if(fabs(SiteDirectCoordWRTslab[i]-1.0)<_slab_epsilon_) {
		  SiteDirectCoordWRTslab[i]=1.0;
		} 
	      }  
	      //for(i=1;i<=3;i++) { cerr << SiteDirectCoordWRTslab[i] << " ";}; cerr << " SiteDirectCoordWRTslab" << endl;  
	      //if(SiteCartCoord[1]==0.5 && Element==2) {cin >> Help_double;}

	      if(SiteDirectCoordWRTslab[1]>=0.0 && SiteDirectCoordWRTslab[1]<1.0 &&
		 SiteDirectCoordWRTslab[2]>=0.0 && SiteDirectCoordWRTslab[2]<1.0 &&
		 SiteDirectCoordWRTslab[3]>=0.0 && SiteDirectCoordWRTslab[3]<1.0) {
		NumSites[Element]++;
		ListSiteDirectCoordWRTslab[Element].resize(NumSites[Element]+1);  
		ListSiteCartCoord[Element].resize(NumSites[Element]+1);  
		ListSiteDirectCoordWRTslab[Element][NumSites[Element]].resize(4);
		ListSiteCartCoord[Element][NumSites[Element]].resize(4);
		for(i=1;i<=3;i++) {
		  ListSiteDirectCoordWRTslab[Element][NumSites[Element]][i]=SiteDirectCoordWRTslab[i];
		}
		for(i=1;i<=3;i++) {
		  ListSiteCartCoord[Element][NumSites[Element]][i]=SiteCartCoord[i];
		}
		//for(i=1;i<=3;i++) { cerr << SiteDirectCoordWRTslab[i] << " ";} cerr << "  Direct" << Element << endl;
		//for(i=1;i<=3;i++) { cerr << SiteCartCoord[i] << " ";} cerr << "  Cart" << Element << endl;
	      }
	    }
	  }
	}
      }  
    }
    //-----------  SORTING atoms inside slab u.c. by Layers -------------------------------------//  
    vector<int>  NumInLayer(1,0);
    vector<double> LayerDirectCoordWRTslab3(1,0);
    vector< vector< vector<int> > > AtomInLayer;
    NumLayers=0;
    for(k=1;k<=NumberElements;k++) {
      for(i=1;i<=NumSites[k];i++) {
	Layer=1;Continue=true;
	while (Layer<=NumLayers && Continue==true) {
	  if(fabs(ListSiteDirectCoordWRTslab[k][i][3]-LayerDirectCoordWRTslab3[Layer])<_slab_epsilon_) {
	    NumInLayer[Layer]++;
	    AtomInLayer[Layer].resize(NumInLayer[Layer]+1);  
	    AtomInLayer[Layer][NumInLayer[Layer]].resize(3);
	    AtomInLayer[Layer][NumInLayer[Layer]][1]=k;
	    AtomInLayer[Layer][NumInLayer[Layer]][2]=i;
	    Continue=false;
	  }
	  Layer=Layer+1;
	}
	if(Continue==true) {
	  NumLayers++;
	  LayerDirectCoordWRTslab3.push_back(ListSiteDirectCoordWRTslab[k][i][3]);
	  NumInLayer.push_back(1);
	  AtomInLayer.resize(NumLayers+1);  
	  AtomInLayer[NumLayers].resize(2);  
	  AtomInLayer[NumLayers][1].resize(3);
	  AtomInLayer[NumLayers][1][1]=k;
	  AtomInLayer[NumLayers][1][2]=i;
	}
      }
    }
    /*
      for(Layer=1; Layer<=NumLayers; Layer++) {
      cerr << "------ " << LayerDirectCoordWRTslab3[Layer] << endl;
      for(k=1;k<=NumInLayer[Layer];k++) {
      cerr << AtomInLayer[Layer][k][1] << "  " << AtomInLayer[Layer][k][2] << endl;  
      for(i=1;i<=3;i++) {
      cerr << ListSiteDirectCoordWRTslab[AtomInLayer[Layer][k][1]][AtomInLayer[Layer][k][2]][i] << " ";
      }      cerr  << endl;      }  }
      cerr << "2--------------------" << endl;
    */
    //-----------  Order Layers by third coordinate wrt S3 ------------------------------------//  
    vector<int> OrderLayers(NumLayers+1);
    vector<double> LayerCoord3(NumLayers+1);
    LayerCoord3[0]=-100;
    for(Layer=1; Layer<=NumLayers; Layer++) {
      LastCandidateCoord3=1000;  
      for(k=1;k<=NumLayers;k++) {  
	CurrCandidateCoord3=LayerDirectCoordWRTslab3[k];
	if(CurrCandidateCoord3>LayerCoord3[Layer-1] && CurrCandidateCoord3<LastCandidateCoord3) {
	  LastCandidateCoord3=CurrCandidateCoord3;OrderLayers[Layer]=k;
	}
      }
      LayerCoord3[Layer]=LayerDirectCoordWRTslab3[OrderLayers[Layer]];
    }


    //for(Layer=1; Layer<=NumLayers; Layer++) {cerr << OrderLayers[Layer]<< "  S3=" << LayerDirectCoordWRTslab3[OrderLayers[Layer]]  << endl;}
    //cerr << "--------------------" << endl;


    //-----------  Shorten S3 (third, "normal to layers" basis vector of slab uc) and project sites coordinates wrt new S3 ------------------------------------//  

    if((NumFilledLayers+NumEmptyLayers) < NumLayers) {
      for(j=1;j<=3;j++) {
	Slab_Basis_CartCoord[3][j]*=LayerDirectCoordWRTslab3[OrderLayers[NumFilledLayers+NumEmptyLayers+1]];
      }

      for(k=1;k<=NumberElements;k++) {
	for(i=1;i<=NumSites[k];i++) {
	  ListSiteDirectCoordWRTslab[k][i][3]/=LayerDirectCoordWRTslab3[OrderLayers[NumFilledLayers+NumEmptyLayers+1]];
	} 
      }
    }

    //Adjust the number of atoms in u.c. for true number of layers
    //for(k=1;k<=NumberElements;k++) {NumSites[k]=0;}
    //for(k=1;k<=NumberElements;k++) {
    //  for(Layer=NumFilledLayers; Layer>=1; Layer--) {
    //  for(i=1;i<=NumInLayer[OrderLayers[Layer]];i++) {
    //  if(AtomInLayer[OrderLayers[Layer]][i][1]==k)
    //  {NumSites[k]++;}
    //} } }
    for(k=1;k<=NumberElements;k++) {
      NumSites[k]=0;
    } //here
    for(k=1;k<=NumberElements;k++) {
      for(Layer=NumFilledLayers; Layer>=1; Layer--) {
	for(i=1;i<=NumInLayer[OrderLayers[Layer]];i++) {
	  if(AtomInLayer[OrderLayers[Layer]][i][1]==k) {
	    NumSites[k]++;
	  }
	}
      }
    }
  
    //-----------  SCREEN OUTPUT ------------------------------------//  
    if(0) {
      cerr << "--------------------" << endl;
      cerr << Title << ", Slab " << hkl(1) <<hkl(2)<<hkl(3)<< ": ("<<NumFilledLayers<<" full + "<<NumEmptyLayers<<" empty)*("<<In_Plane_Multiplication[1] << "x" << In_Plane_Multiplication[2] <<" InPlaneMultipl)"<<endl;
      cerr << LattPar << endl;
      for(i=1;i<=3;i++) {
	for(j=1;j<=3;j++) {
	  cerr << Slab_Basis_CartCoord(i,j) << " ";
	} 
	cerr << endl;
      }
      for(i=1;i<=NumberElements;i++) {
	cerr << NumSites[i] << " ";
      } 
      cerr << endl << "Direct" << endl;
      for(k=1;k<=NumberElements;k++) {
	for(Layer=NumFilledLayers; Layer>=1; Layer--) {
	  
	  for(i=1;i<=NumInLayer[OrderLayers[Layer]];i++) {
	    if(AtomInLayer[OrderLayers[Layer]][i][1]==k) {
	      for(j=1;j<=3;j++) {
		cerr << ListSiteDirectCoordWRTslab[k][AtomInLayer[OrderLayers[Layer]][i][2]][j] << " ";
	      } 
	      cerr << " type=" << k << endl;
	    }
	  } 
	} 
      }
      cerr << "------- With Cartesian Coordinates -------" << endl;
      cerr << Title << ", Slab " << hkl(1) <<hkl(2)<<hkl(3)<< ": ("<<NumFilledLayers<<" full + "<<NumEmptyLayers<<" empty)*("<<In_Plane_Multiplication[1] << "x" << In_Plane_Multiplication[2]<<" InPlaneMultipl)"<<endl;
      cerr << LattPar << endl;
      for(i=1;i<=3;i++) {
	for(j=1;j<=3;j++) {
	  cerr << Slab_Basis_CartCoord(i,j) << " ";} cerr << endl;
      }
      for(i=1;i<=NumberElements;i++) {
	cerr << NumSites[i] << " ";} cerr << endl;
      cerr << "Cartesian" << endl;
      for(k=1;k<=NumberElements;k++) {
	for(Layer=NumFilledLayers; Layer>=1; Layer--) {
          for(i=1;i<=NumInLayer[OrderLayers[Layer]];i++) {
	    if(AtomInLayer[OrderLayers[Layer]][i][1]==k) {
	      for(j=1;j<=3;j++) {
		cerr << ListSiteCartCoord[k][AtomInLayer[OrderLayers[Layer]][i][2]][j] << " ";
	      } 
	      cerr << " type=" << k << endl;
	    }
	  } 
	} 
      }
      cerr << "------- Initial Non-slab POSCAR -------" << endl;
      cerr << Title << endl;
      cerr << LattPar << endl;
      for(i=1;i<=3;i++) {
	for(j=1;j<=3;j++) {
	  cerr << UnitCellVector(i,j) << " ";} cerr << endl;
      }
      for(i=1;i<=NumberElements;i++) {
	cerr << NumAtomsForElementUC[i] << " ";
      } 
      cerr << endl;
      cerr << "Direct  (Must be D for this Code)" << endl;
      
      for(k=1;k<=NumberElements;k++) {
	for(i=1;i<=NumAtomsForElementUC[k];i++) {
	  for(j=1;j<=3;j++) {
	    cerr << AtomDirCoords[k][i][j] << " ";
	  } 
	  cerr << k << endl;
	} 
      }
      cerr << endl;
      for(i=1;i<=3;i++) {
	cerr << LattPar*AbsValueSlab_Basis[i]<< " ";
      } 
      cerr << 180.0/PI*acos(Layer0BasisCosAngle) << endl;
    }
    //-----------  FILE OUTPUT ------------------------------------//  
    if(0) {
      stringstream fout;
      fout << Title << ", Slab " << hkl(1) <<hkl(2)<<hkl(3)<< ": ("<<NumFilledLayers<<" full + "<<NumEmptyLayers<<" empty)*("<<In_Plane_Multiplication[1] << "x" << In_Plane_Multiplication[2]<<" InPlaneMultipl)"<<endl;
      fout << LattPar << endl;
      for(i=1;i<=3;i++) {
	for(j=1;j<=3;j++) {
	  fout << Slab_Basis_CartCoord(i,j) << " ";
	}
	fout << endl;
      }
      for(i=1;i<=NumberElements;i++) {
	fout << NumSites[i] << " ";
      }
      fout << endl;
      fout << "Cartesian" << endl;
      for(k=1;k<=NumberElements;k++) {
	for(Layer=NumFilledLayers; Layer>=1; Layer--) {
	  for(i=1;i<=NumInLayer[OrderLayers[Layer]];i++) {
	    if(AtomInLayer[OrderLayers[Layer]][i][1]==k) {
	      for(j=1;j<=3;j++) {
		fout << ListSiteCartCoord[k][AtomInLayer[OrderLayers[Layer]][i][2]][j] << " ";
	      } 
	      fout << " type=" << k << endl;
	    }
	  } 
	} 
      }
      fout << "------- With Direct Coordinates -------" << endl;  
      fout << Title << ", Slab " << hkl(1) <<hkl(2)<<hkl(3)<< ": ("<<NumFilledLayers<<" full + "<<NumEmptyLayers<<" empty)*("<<In_Plane_Multiplication[1] << "x" << In_Plane_Multiplication[2]<<" InPlaneMultipl)"<<endl;
      fout << LattPar << endl;
      for(i=1;i<=3;i++) {
	for(j=1;j<=3;j++) {
	  fout << Slab_Basis_CartCoord(i,j) << " ";
	}
	fout << endl;
      }
      for(i=1;i<=NumberElements;i++) {
	fout << NumSites[i] << " ";
      }
      fout << endl;
      fout << "Direct" << endl;
      
      for(k=1;k<=NumberElements;k++) {
	for(Layer=NumFilledLayers; Layer>=1; Layer--) {
	  for(i=1;i<=NumInLayer[OrderLayers[Layer]];i++) {
	    if(AtomInLayer[OrderLayers[Layer]][i][1]==k) {
	      for(j=1;j<=3;j++) {
		fout << ListSiteDirectCoordWRTslab[k][AtomInLayer[OrderLayers[Layer]][i][2]][j] << " ";
	      }
	      fout << " type=" << k << endl;
	    }
	  } 
	}
      }
      
      fout << "------- Initial Non-slab POSCAR -------" << endl;
      fout << Title << endl;
      fout << LattPar << endl;
      for(i=1;i<=3;i++) {
	for(j=1;j<=3;j++) {
	  fout << UnitCellVector(i,j) << " ";
	} fout << endl;
      }
      for(i=1;i<=NumberElements;i++) {
	fout << NumAtomsForElementUC[i] << " ";
      }
      fout << endl;
      fout << "Direct  (Must be D for this Code)" << endl;
      for(k=1;k<=NumberElements;k++) {
	for(i=1;i<=NumAtomsForElementUC[k];i++) {
	  for(j=1;j<=3;j++) {
	    fout << AtomDirCoords[k][i][j] << " ";
	  } fout << k << endl;
	}
      }
      fout << "Lengths of basis vectors (A):" << endl;  
      for(i=1;i<=3;i++) { 
	fout << LattPar*AbsValueSlab_Basis[i]<< " ";
      } 
      fout << endl;    
      fout << "Angle between basis vectors in plane:" << endl;  
      fout << 180.0/PI*acos(Layer0BasisCosAngle) << endl;  
      aurostd::stringstream2file(fout,file_OUT); 
    }
    
    //Output xstructure for slab
    xstructure str_out("");

    stringstream sout;
    
    sout << str_in.title << ", Slab " << hkl(1) <<hkl(2)<<hkl(3)<< ": ("<<NumFilledLayers<<" full + "<<NumEmptyLayers<<" empty)*("<<In_Plane_Multiplication[1] << "x" << In_Plane_Multiplication[2]<<" InPlaneMultipl)";
    str_out.title=sout.str();
    str_out.scale= LattPar;
    
    for(i=1;i<=3;i++) {
      for(j=1;j<=3;j++) {
	str_out.lattice(i,j)=Slab_Basis_CartCoord(i,j);
      }
    }
    str_out.FixLattices();  //CO 180202
    
    //str_out.num_each_type.clear();
    //str_out.comp_each_type.clear();
    
    _atom newatom;
    
    iatom=0;
    for(k=1;k<=NumberElements;k++) {
      for(Layer=NumFilledLayers; Layer>=1; Layer--) {
	for(i=1;i<=NumInLayer[OrderLayers[Layer]];i++) {
	  if(AtomInLayer[OrderLayers[Layer]][i][1]==k) {
	    iatom++;
	    for(j=1;j<=3;j++) {
	      newatom.fpos[j]=ListSiteDirectCoordWRTslab[k][AtomInLayer[OrderLayers[Layer]][i][2]][j];
	    }
            newatom.cpos=F2C(str_out.scale,str_out.lattice,newatom.fpos); //CO 180202
	    str_out.AddAtom(newatom);
	  }
	}
      }
    }

    //CO 180202
    if(LDEBUG){
      cerr << "PRINTING OUT STRUCTURE ATTRIBUTES" << endl;
      cerr << "str_out.atoms.size()=" << str_out.atoms.size() << endl;
      cerr << "str_out.num_each_type.size()=" << str_out.num_each_type.size() << endl;
      cerr << "str_out.comp_each_type.size()=" << str_out.comp_each_type.size() << endl;
      for(uint i=0;i<str_out.num_each_type.size();i++){
        cerr << "str_out.num_each_type[i]=" << str_out.num_each_type[i] << endl;
      }
      cerr << str_out << endl;
    }
  
    //CO 180202 - this is obsolete, it is done INSIDE AddAtom()
    //str_out.num_each_type.clear();
    //str_out.comp_each_type.clear();
    //for(i=0;i<NumberElements;i++) {
    //  str_out.num_each_type.push_back(NumSites[i+1]);
    //  str_out.comp_each_type.push_back((double) NumSites[i+1]);
    //}

    for(i=0;i<int(str_out.atoms.size());i++) {
      //str_out.atoms.at(i).name=str_in.SpeciesLabel(i);
      //  cout << str_out.atoms.at(i)<<endl;
      str_out.atoms.at(i).name_is_given=TRUE;
    }

    for(j=0;j<NumberElements;j++) {
      for(i=NumSites[j];i<NumSites[j+1]+NumSites[j];i++) {
        //CO 180202 - added safety
        if(i>(int)str_out.atoms.size()-1){
          cerr << "pflow::MAKE_SLAB: ERROR - not as many atoms were created as cxpected (likely a problem with AddAtom())" << endl;
          cerr << "Exiting!" << endl;
          exit(1);
        }
	str_out.atoms.at(i).name=str_in.SpeciesLabel(j);
      }
    }

    //for(i=0;i<NumberElements;i++)
    //str_out.num_each_type[i]=(NumSites[i+1]);

    if(LDEBUG) cerr << "slab::MAKE_SLAB: END" << endl;
    return str_out;

  }
} // namespace slab

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *                                                                         *
// ***************************************************************************
