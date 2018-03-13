// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *                                                                         *
// ***************************************************************************
// Stefano Curtarolo and Dane Morgan

#ifndef _AFLOW_PFLOW_CPP_
#define _AFLOW_PFLOW_CPP_

#include "aflow_pflow.h"

// ***************************************************************************
// Simple interfaces
// ***************************************************************************
int SignNoZero(const double& x) {return (int) aurostd::signnozero(x);}
int Nint(const double& x) {return (int) aurostd::nint(x);}
int Sign(const double& x) {return (int) aurostd::sign(x);}

// ***************************************************************************
// Matrix classes in Dane Morgan Style
// ***************************************************************************
namespace pflow {
  // constructor
  
  //   template<class utype>
  //   matrix<utype>::matrix(void) { // default
  //   }
  
  template<class utype>
  matrix<utype>::matrix(const int m) { // specifying rows
    vector<utype> v;
    mat=vector<vector<utype> > (m,v);
  }
  
  template<class utype>
  matrix<utype>::matrix(const int m, const int n) { // specifying rows and columns
    vector<utype> v(n);
    mat=vector<vector<utype> > (m,v);
  }
  
  template<class utype>
  matrix<utype>::matrix(const int m, const vector<utype>& inutypevec) { // specifying rows as vectors.
    mat=vector<vector<utype> > (m,inutypevec);
  }
  
  template<class utype>
  matrix<utype>::matrix(const int m, const int n, const utype& inutype) { // specifying rows and columns and initial values
    vector<utype> v(n);
    mat=vector<vector<utype> > (m,v);
    for(uint i=0;i<mat.size();i++) {
      for(uint j=0;j<mat[i].size();j++) {
	mat[i][j]=inutype;
      }
    }
  }
  
  // accessors
  template<class utype>
  void matrix<utype>::print(void) {
    cout.setf(std::ios::fixed,std::ios::floatfield);
    cout.precision(4);
    for(uint i=0;i<mat.size();i++) {
      cout << "  ";
      for(uint j=0;j<mat[i].size();j++) {
	cout << " " << mat[i][j];
      }
      cout << endl;
    }
  }
  
  //   template<class utype>
  //   inline int matrix<utype>::size(void) const{
  //     return (int) mat.size();
  //   }
  

  template<class utype>
  matrix<utype> matrix<utype>::transpose() const{
    matrix<utype> tmat;
    uint m=mat.size();
    if(m==0) return tmat;
    uint n=mat[0].size();
    tmat = matrix<utype> (n,m);
    for(uint i=0;i<m;i++) {
      for(uint j=0;j<n;j++) {
	tmat[j][i]=mat[i][j];
      }
    }
    return tmat;
  }
  
  // template<class utype>
  // vector<vector<utype> >::iterator matrix<utype>::begin() {
  // return mat.begin();
  // }
  //
  // template<class utype>
  // vector<vector<utype> >::iterator matrix<utype>::end() {
  // return mat.end();
  // }
  
  // operator
  //   template<class utype>
  //   vector<utype>& matrix<utype>::operator[] (const int index) {
  //     assert(index>=0 && index<=mat.size());
  //     return mat[index];
  //   }
  //   template<class utype>
  //   const vector<utype>& matrix<utype>::operator[] (const int index) const {
  //     assert(index>=0 && index<=mat.size());
  //     return mat[index];
  
  //   }
  
  template<class utype>
  const matrix<utype>& matrix<utype>::operator=(const matrix<utype> &b) {
    if(this != &b) {
      uint m=b.mat.size();
      uint n=0;
      mat=vector<vector<utype> > (m);
      for(uint i=0;i<m;i++) {
	n=b.mat[i].size();
	mat[i]=vector<utype> (n);
	for(uint j=0;j<n;j++) {
	  mat[i][j]=b.mat[i][j];
	}
      }
    }
    return *this;
  }

  // mutators
  //  template<class utype>
  // void matrix<utype>::push_back(const vector<utype>& inutypevec) {
  //   mat.push_back(inutypevec);
  //  }

  //  template<class utype>
  // void matrix<utype>::pop_back() {
  //   mat.pop_back();
  //  }

  template<class utype>
  void matrix<utype>::vecvec2mat(const vector<vector<utype> >& inVV) {
    mat=vector<vector<utype> > (inVV.size());
    for(uint i=0;i<mat.size();i++) {
      mat[i]=vector<utype> (inVV[i].size());
      for(uint j=0;j<mat[i].size();j++) {
	mat[i][j]=inVV[i][j];
      }
    }
  }

  // template<class utype>
  // void matrix<utype>::clear() {
  //   mat.clear();
  // }

  template<class utype>
  void matrix<utype>::vec2mat(const vector<utype>& inV) {
    mat=vector<vector<utype> > (1);
    mat[0]=vector<utype> (inV.size());
    for(uint j=0;j<mat[0].size();j++) {
      mat[0][j]=inV[j];
    }
  }

  // template<class utype>
  // void matrix<utype>::insert(const int& id, const vector<utype>& inV) {
  //  mat.insert(mat.begin()+id,inV);
  // }

  // template<class utype>
  // void matrix<utype>::erase(const int id) {
  // vector<utype>::iterator p=mat.begin()+id;
  // mat.erase(p);
  // }
  
  // template<class utype>
  // void matrix<utype>::erase_col(const int id) {
  //   for(int i=0;i<mat.size();i++) {
  //           vector<utype>::iterator p=mat[i].begin()+id;
  //     if(id<mat[i].size()) mat[i].erase(p);
  //   }
  // }

  // template<class utype>
  // void matrix<utype>::erase(const int id) {
  //   vector<vector<utype> >::iterator p=mat.begin()+id;
  //   mat.erase(p);
  // }

  template <class utype> matrix<utype>
  xmatrix2matrix(const xmatrix<utype>& _xmatrix) {
    int isize=_xmatrix.rows,jsize=_xmatrix.cols;
    matrix<utype> _matrix(isize,jsize);
    for(register int i=0;i<isize;i++)
      for(register int j=0;j<jsize;j++)
	_matrix[i][j]=_xmatrix(i+_xmatrix.lrows,j+_xmatrix.lcols);
    return _matrix;
  }

//   matrix<double> xmatrix2matrix(const xmatrix<double>& _xmatrix) {
//     int isize=_xmatrix.rows,jsize=_xmatrix.cols;
//     matrix<double> _matrix(isize,jsize);
//     for(register int i=0;i<isize;i++)
//       for(register int j=0;j<jsize;j++)
// 	_matrix[i][j]=_xmatrix(i+_xmatrix.lrows,j+_xmatrix.lcols);
//     return _matrix;
//   }
  
  template <class utype> xmatrix<utype>
  matrix2xmatrix(const matrix<utype>& _matrix) {
    int isize=_matrix.size(),jsize=_matrix[0].size();
    xmatrix<utype> _xmatrix(isize,jsize);
    for(register int i=1;i<=isize;i++)
      for(register int j=1;j<=jsize;j++)
	_xmatrix(i,j)=_matrix[i-1][j-1];
    return _xmatrix;
  }

  

}

// ***************************************************************************
// Function PutInCompact
// ***************************************************************************
// Make a structure where all the atoms are all the
// atoms are mapped through the unit and neighbours cells
// to minimixe the shortest possible bond with an adjacent atom. (SC 6 Aug 04)
xstructure PutInCompact(const xstructure& a) {
  xstructure sstr=a;
  sstr.BringInCompact();
  return sstr;
}

// ***************************************************************************
// Function PutInCell - same as BringInCell
// ***************************************************************************
// Make a structure with all atoms mapped to their images within the
// unit cell.
xstructure PutInCell(const xstructure& a) {
  xstructure sstr=a;
  sstr.BringInCell();
  return sstr;
}

// ***************************************************************************
// Function WignerSeitz
// ***************************************************************************
// Make a structure where all the atoms are
// mapped to their images in the Wigner-Seitz cell.(SC 10Jan04)
xstructure WignerSeitz(const xstructure& a) {
  xstructure sstr=a;
  // zoology of vectors
  xvector<double> rat(3),xoo(3),yoo(3),zoo(3),xyo(3),xzo(3),yzo(3),xyz(3);
  double axoo,ayoo,azoo,axyo,axzo,ayzo,axyz;
  xoo=sstr.lattice(1);axoo=modulus(xoo);
  yoo=sstr.lattice(2);ayoo=modulus(yoo);
  zoo=sstr.lattice(3);azoo=modulus(zoo);
  xyo=sstr.lattice(1)+sstr.lattice(2);axyo=modulus(xyo);
  xzo=sstr.lattice(1)+sstr.lattice(3);axzo=modulus(xzo);
  yzo=sstr.lattice(2)+sstr.lattice(3);ayzo=modulus(yzo);
  xyz=sstr.lattice(1)+sstr.lattice(2)+sstr.lattice(3);axyz=modulus(xyz);
  double projxoo,projyoo,projzoo,projxyo,projxzo,projyzo,projxyz;

  for(uint iat=0;iat<sstr.atoms.size();iat++) {
    for(int i=-1;i<=1;i++)
      for(int j=-1;j<=1;j++)
	for(int k=-1;k<=1;k++) {
	  rat=sstr.atoms.at(iat).cpos+((double)i)*sstr.lattice(1)+((double)j)*sstr.lattice(2)+((double)k)*sstr.lattice(3);
	  projxoo=scalar_product(rat,xoo)/axoo/axoo;
	  projyoo=scalar_product(rat,yoo)/ayoo/ayoo;
	  projzoo=scalar_product(rat,zoo)/azoo/azoo;
	  projxyo=scalar_product(rat,xyo)/axyo/axyo;
	  projxzo=scalar_product(rat,xzo)/axzo/axzo;
	  projyzo=scalar_product(rat,yzo)/ayzo/ayzo;
	  projxyz=scalar_product(rat,xyz)/axyz/axyz;
	  if((projxoo>-0.5 && projxoo<=0.5) &&
	     (projyoo>-0.5 && projyoo<=0.5) &&
	     (projzoo>-0.5 && projzoo<=0.5) &&
	     (projxyo>-0.5 && projxyo<=0.5) &&
	     (projxzo>-0.5 && projxzo<=0.5) &&
	     (projyzo>-0.5 && projyzo<=0.5) &&
	     (projxyz>-0.5 && projxyz<=0.5)) {
	    sstr.atoms.at(iat).cpos(1)=rat(1);
	    sstr.atoms.at(iat).cpos(2)=rat(2);
	    sstr.atoms.at(iat).cpos(3)=rat(3);
	    i=10;j=10;k=10;
	  }
	}
    sstr.atoms.at(iat).fpos=C2F(sstr.lattice,sstr.atoms.at(iat).cpos);
  }
  return sstr;
}

// **************************************************************************
// GetMom1
// **************************************************************************
// this function returns moment_1 position of the atoms in cartesians
xvector<double> GetMom1(const xstructure& a) {    
  // Get's first moment in cartesian coordinates.
  xvector<double> mom1(3);
  for(uint iat=0;iat<a.atoms.size();iat++)
    mom1=mom1+a.atoms.at(iat).cpos;
  if(a.atoms.size()>0) mom1=((double) a.scale/((double) a.atoms.size()))*mom1;
  else clear(mom1);
  return mom1;
}

// **************************************************************************
// SetMom1
// **************************************************************************
// this function sets moment_1 position of atoms
xstructure SetMom1(const xstructure& a, const xvector<double>& mom1_new) {
  xstructure sstr=a;
  // Get old first moment in cart. coords.
  xvector<double> mom1_old(3);mom1_old=GetMom1(sstr);
  // Get change in first moment in cart. coords.
  xvector<double> dmom1(3);dmom1=mom1_new-mom1_old;
  // Shift scaled cart. positions by dmom1.
  for(uint iat=0;iat<sstr.atoms.size();iat++) {
    sstr.atoms.at(iat).cpos=sstr.scale*sstr.atoms.at(iat).cpos;
    sstr.atoms.at(iat).cpos=sstr.atoms.at(iat).cpos+dmom1;
    sstr.atoms.at(iat).cpos=sstr.atoms.at(iat).cpos/sstr.scale;
    sstr.atoms.at(iat).fpos=C2F(sstr.lattice,sstr.atoms.at(iat).cpos);
  }
  return sstr;
}

// **************************************************************************
// Function Function GetCDispFromOrigin
// **************************************************************************
// This function returns displacement from origin
xvector<double> GetCDispFromOrigin(const _atom& atom) {
  xvector<double> diff(3);
  diff=atom.cpos-atom.corigin;
  return diff;
}

// **************************************************************************
// Function GetDistFromOrigin
// **************************************************************************
// This function returns distance from origin
double GetDistFromOrigin(const _atom& atom) {
  return modulus(GetCDispFromOrigin(atom));
}

// [OBSOLETE] // **************************************************************************
// [OBSOLETE] // Function GetUnitCellRep
// [OBSOLETE] // **************************************************************************
// [OBSOLETE] // Dane Morgan /  Stefano Curtarolo
// [OBSOLETE] // This function calculated the representation of a point
// [OBSOLETE] // in terms of a position in cell 0 and its true unit cell.
// [OBSOLETE] // Output position in cell 0 is given with same coordinate
// [OBSOLETE] // type (Cart or Direct) as input position.
// [OBSOLETE] // Note that there is always an ambiguity of being in a given
// [OBSOLETE] // cell with position 0, or one cell over with position 1.  This
// [OBSOLETE] //is broken by forcing all cell positions to be at 0 if they
// [OBSOLETE] //are within TOL of 1.  This gives consistent cell image locations.
// [OBSOLETE] void GetUnitCellRep(const xvector<double>& ppos,xvector<double>& p_cell0,
// [OBSOLETE] 		    xvector<int>& ijk,const xmatrix<double>& lattice,
// [OBSOLETE] 		    const bool coord_flag) {
// [OBSOLETE]   double TOL=1e-11;
// [OBSOLETE]   xvector<double> cpos(3),fpos(3);
// [OBSOLETE]   // p_cell0=xvector<double> (3);
// [OBSOLETE]   // ijk=xvector<int> (3);
// [OBSOLETE]
// [OBSOLETE]   clear(ijk);clear(p_cell0);
// [OBSOLETE]
// [OBSOLETE]   if(coord_flag==TRUE) // ppos is in cartesian coords
// [OBSOLETE]     fpos=C2F(lattice,ppos);
// [OBSOLETE]
// [OBSOLETE]   for(int ic=1;ic<=3;ic++) {
// [OBSOLETE]     int s;
// [OBSOLETE]     s=SignNoZero(fpos(ic));
// [OBSOLETE]     if(s>0) ijk(ic)=int(fpos(ic));
// [OBSOLETE]     else    ijk(ic)=int(fpos(ic))-1;
// [OBSOLETE]     p_cell0(ic)=fpos(ic)-ijk(ic);
// [OBSOLETE]     // Boundary Correction.  p_cell0 is >=0 and <1.  
// [OBSOLETE]     // If p_cell0 is within TOL of 1
// [OBSOLETE]     // then set it to 0 and shift ijk up by 1.
// [OBSOLETE]     if(abs(1-p_cell0(ic))<TOL) {
// [OBSOLETE]       p_cell0(ic)=0;
// [OBSOLETE]       ijk(ic)=ijk(ic)+1;
// [OBSOLETE]     }
// [OBSOLETE]   }
// [OBSOLETE]   if(coord_flag==TRUE) // ppos is in cartesian coords
// [OBSOLETE]     p_cell0=F2C(lattice,p_cell0);
// [OBSOLETE] }

// **************************************************************************
// Function ConvertAtomToLat
// **************************************************************************
// This function makes sure that an atoms fractional coords
// and unit cell parameters are consistent with an input lattice.
_atom ConvertAtomToLat(const _atom& in_at, const xmatrix<double>& lattice) {
  xvector<double> cpos(3);
  xvector<double> p_cell0(3);
  xvector<int> ijk;
  _atom out_at=in_at;
  cpos=in_at.cpos;
  GetUnitCellRep(cpos,p_cell0,ijk,lattice,TRUE);
  out_at.fpos=C2F(lattice,cpos);
  out_at.ijk=ijk;
  return out_at;
}

// [OBSOLETE] // **************************************************************************
// [OBSOLETE] // Function COMPARE_GetNeighData
// [OBSOLETE] // **************************************************************************
// [OBSOLETE] // For sort algorithm in GetNeighData
// [OBSOLETE] class compare_GetNeighData {
// [OBSOLETE] public:  
// [OBSOLETE]   int operator()(const _atom& a, const _atom& b) {
// [OBSOLETE]     bool LDEBUG=(FALSE || XHOST.DEBUG);
// [OBSOLETE]     if(LDEBUG) cerr << "operator()(const _atom& a, const _atom& b)" << endl;
// [OBSOLETE]     double tol=1e-15;
// [OBSOLETE]     // Sort by distance
// [OBSOLETE]     if(aurostd::isequal(GetDistFromOrigin(a),GetDistFromOrigin(b),tol)) {
// [OBSOLETE]       if(LDEBUG) cerr << "operator()(const _atom& a, const _atom& b)  isequal(GetDistFromOrigin(a),GetDistFromOrigin(b),tol)" << endl;
// [OBSOLETE]       // Sort by unit cell values
// [OBSOLETE]       if(a.name==b.name) {
// [OBSOLETE] 	if(LDEBUG) cerr << "operator()(const _atom& a, const _atom& b)  a.name==b.name" << endl;
// [OBSOLETE]         xvector<int> ijka=a.ijk;
// [OBSOLETE]         xvector<int> ijkb=b.ijk;
// [OBSOLETE]         int va=100*ijka(1)+10*ijka(2)+1*ijka(3);
// [OBSOLETE]         int vb=100*ijkb(1)+10*ijkb(2)+1*ijkb(3);
// [OBSOLETE]         return va<vb;
// [OBSOLETE]       } else {
// [OBSOLETE] 	if(LDEBUG) cerr << "operator()(const _atom& a, const _atom& b)  a.name!=b.name" << endl;
// [OBSOLETE] 	return a.name<b.name;
// [OBSOLETE]       }
// [OBSOLETE]     } else{
// [OBSOLETE]       if(LDEBUG) cerr << "operator()(const _atom& a, const _atom& b)  !isequal(GetDistFromOrigin(a),GetDistFromOrigin(b),tol)" << endl;
// [OBSOLETE]       return GetDistFromOrigin(a)<GetDistFromOrigin(b);
// [OBSOLETE]     }
// [OBSOLETE]     // Sort by name
// [OBSOLETE]     // if(a.name==b.name) {
// [OBSOLETE]     //  return GetDistFromOrigin(a)<GetDistFromOrigin(b);
// [OBSOLETE]     // }
// [OBSOLETE]     // else{
// [OBSOLETE]     //  return a.name<b.name;
// [OBSOLETE]     // }
// [OBSOLETE]   }
// [OBSOLETE] };

// [OBSOLETE] // **************************************************************************
// [OBSOLETE] // Function GetNeighData
// [OBSOLETE] // **************************************************************************
// [OBSOLETE] // This function collects all the neighbor data between
// [OBSOLETE] // rmin and rmax and stores it for each atom in a vector of atom
// [OBSOLETE] // objects in order of increasing distance.  
// [OBSOLETE] namespace pflow {
// [OBSOLETE]   void GetNeighData(const deque<_atom>& in_atom_vec, const xstructure& in_str,
// [OBSOLETE] 		    const double& rmin, const double& rmax,
// [OBSOLETE] 		    deque<deque<_atom> >& neigh_mat) {
// [OBSOLETE]     double epsilon=1.0e-3;  // sometimes you have wrong images due to roundoff
// [OBSOLETE]     deque<_atom> neigh_vec;
// [OBSOLETE]     // Get data from str.
// [OBSOLETE]     xstructure sstr=in_str;
// [OBSOLETE]     // Set scale to 1 so you don't need to rescale coordinates.
// [OBSOLETE]     sstr=ReScale(sstr,1.0);
// [OBSOLETE]     xmatrix<double> original_lattice(3,3);
// [OBSOLETE]     original_lattice=sstr.lattice;
// [OBSOLETE]     // Use Niggli cell to avoid missing neighbors in angular cells.
// [OBSOLETE]     // This involves converting your atoms into the 0th Niggli cell,
// [OBSOLETE]     // keeping track of the unit cell you came from, and then shifting back
// [OBSOLETE]     // at the end.  This shifting must be done for both atom_vec (at the
// [OBSOLETE]     // beginning) and neigh_mat (at the end).
// [OBSOLETE]     //sstr.NiggliUnitCellForm();  
// [OBSOLETE]     sstr.MinkowskiBasisReduction();
// [OBSOLETE]     //    cerr << sstr << endl;
// [OBSOLETE]
// [OBSOLETE]     xmatrix<double> lat(3,3);
// [OBSOLETE]     lat=sstr.lattice;
// [OBSOLETE]  
// [OBSOLETE]      //  Convert input atom_vec to Niggli lattice
// [OBSOLETE]     deque<_atom> atom_vec=in_atom_vec;
// [OBSOLETE]     xmatrix<int>  atom_shifts((int) atom_vec.size()-1,3,0,1);
// [OBSOLETE]     for(uint ia=0;ia<atom_vec.size();ia++) {
// [OBSOLETE]       _atom a=atom_vec.at(ia);
// [OBSOLETE]       xvector<double> p_cell0(3);
// [OBSOLETE]       xvector<int>    atom_shifts_temp(3);
// [OBSOLETE]       GetUnitCellRep(a.cpos,p_cell0,atom_shifts_temp,lat,1);
// [OBSOLETE]       atom_vec.at(ia).cpos=p_cell0;
// [OBSOLETE]       atom_vec.at(ia)=ConvertAtomToLat(atom_vec.at(ia),lat);
// [OBSOLETE]       atom_shifts(ia,1)=atom_shifts_temp(1);
// [OBSOLETE]       atom_shifts(ia,2)=atom_shifts_temp(2);
// [OBSOLETE]       atom_shifts(ia,3)=atom_shifts_temp(3);
// [OBSOLETE]     }    
// [OBSOLETE]    
// [OBSOLETE]     // Create vector of all atoms in all periodic images needed so that
// [OBSOLETE]     // for atoms in the unit cell every possible neighbor within rmax
// [OBSOLETE]     // is included.  This can be done by looping over all atoms in
// [OBSOLETE]     // all unit cells such that if the unit cells are given by
// [OBSOLETE]     // n0*v0,n1*v1,n2*v2, then ni_max*vi>rmax for all i=0->2.  This is not
// [OBSOLETE]     // rigorous but seems to work quite well.  It failed for very
// [OBSOLETE]     // angular cells, but by using the Niggli reduced cell this pitfall
// [OBSOLETE]     // seems to be avoided (but no proof, so be careful).
// [OBSOLETE]
// [OBSOLETE]     int imax,jmax,kmax;
// [OBSOLETE]     // Find imax
// [OBSOLETE]     // (algorithm 1) approximate
// [OBSOLETE]     imax=(int)max(rmax/modulus(lat(1)),rmax/modulus(lat(2)),rmax/modulus(lat(3)))+2;
// [OBSOLETE]     jmax=imax;
// [OBSOLETE]     kmax=imax;
// [OBSOLETE]     // (algorithm 2) exact, the requirement is that atoms must in incell.
// [OBSOLETE]     // we should find nlat1,nlat2,nlat3
// [OBSOLETE]     // where nlat1 is the lat1 component that is perpendicular to the plane made by lat2 & lat3.
// [OBSOLETE]     // then imax=ceil(rmax/nlat1)
// [OBSOLETE]     // This algorithm is implemented in function LatticeDimensionSphere
// [OBSOLETE]     xvector<int> dims(3);
// [OBSOLETE]     dims=LatticeDimensionSphere(lat,rmax);
// [OBSOLETE]     imax=dims(1);
// [OBSOLETE]     jmax=dims(2);
// [OBSOLETE]     kmax=dims(3);
// [OBSOLETE]     xvector<int> ijk(3);
// [OBSOLETE]     deque<_atom> all_atom_vec;
// [OBSOLETE]     // latt maybe a rotated version of POSCAR, so need to incell-ized the fpos and cpos
// [OBSOLETE]     // sstr.BringInCell(); // with roundoff
// [OBSOLETE]     sstr.BringInCell(-1.0); // no roundoff
// [OBSOLETE]     for(ijk(1)=-imax;ijk(1)<=imax;ijk(1)++) {
// [OBSOLETE]       for(ijk(2)=-jmax;ijk(2)<=jmax;ijk(2)++) {
// [OBSOLETE] 	for(ijk(3)=-kmax;ijk(3)<=kmax;ijk(3)++) {
// [OBSOLETE] 	  xvector<double> ctpos(3);
// [OBSOLETE] 	  xvector<double> ftpos(3);
// [OBSOLETE] 	  for(uint iat=0;iat<sstr.atoms.size();iat++) {
// [OBSOLETE] 	    _atom a;
// [OBSOLETE] 	    a=sstr.atoms.at(iat);
// [OBSOLETE] 	    a.name=sstr.atoms.at(iat).name;
// [OBSOLETE] 	    a.number=iat;
// [OBSOLETE] 	    a.ijk=ijk;
// [OBSOLETE] 	    for(int ic=1;ic<=3;ic++)
// [OBSOLETE] 	      ctpos(ic)=sstr.atoms.at(iat).cpos(ic)+ijk(1)*lat(1,ic)+ijk(2)*lat(2,ic)+ijk(3)*lat(3,ic);
// [OBSOLETE] 	    ftpos=C2F(lat,ctpos);
// [OBSOLETE] 	    a.cpos=ctpos;
// [OBSOLETE] 	    a.fpos=ftpos;
// [OBSOLETE] 	    a.type=sstr.atoms.at(iat).type;
// [OBSOLETE] 	    all_atom_vec.push_back(a);
// [OBSOLETE] 	  } // iat
// [OBSOLETE] 	} // ijk
// [OBSOLETE]       } // ijk
// [OBSOLETE]     } // ijk
// [OBSOLETE]  
// [OBSOLETE]     // Now build neighbors list for each atom on atom list.  Each
// [OBSOLETE]     // neighbor list will be row in the neigh_mat matrix.
// [OBSOLETE]     for(uint ia1=0;ia1<atom_vec.size();ia1++) {
// [OBSOLETE]       xvector<double> corigin(3);
// [OBSOLETE]       _atom at=atom_vec.at(ia1);
// [OBSOLETE]       corigin = at.cpos;
// [OBSOLETE]       at.corigin=corigin;
// [OBSOLETE]       deque<_atom> neigh_vec;
// [OBSOLETE]       neigh_vec.push_back(at);
// [OBSOLETE]       for(uint ia2=0;ia2<all_atom_vec.size();ia2++) {
// [OBSOLETE] 	double dist = AtomDist(at,all_atom_vec.at(ia2));
// [OBSOLETE] 	if(dist<=rmax && dist>=rmin && dist>=epsilon) {
// [OBSOLETE] 	  all_atom_vec.at(ia2).corigin=corigin;
// [OBSOLETE] 	  neigh_vec.push_back(all_atom_vec.at(ia2));
// [OBSOLETE] 	} // if
// [OBSOLETE]       } // ia2
// [OBSOLETE]       sort(neigh_vec.begin()+1,neigh_vec.end(),compare_GetNeighData());
// [OBSOLETE]       neigh_mat.push_back(neigh_vec);
// [OBSOLETE]     } // ia1
// [OBSOLETE]
// [OBSOLETE]     //  Convert neigh_mat from Niggli cell to original lattice
// [OBSOLETE]     for(uint ia=0;ia<neigh_mat.size();ia++) {
// [OBSOLETE]       for(uint ja=0;ja<neigh_mat.at(ia).size();ja++) {
// [OBSOLETE] 	xvector<double> fpos(3);
// [OBSOLETE] 	xvector<double> cpos(3);
// [OBSOLETE] 	_atom a;
// [OBSOLETE] 	a=neigh_mat.at(ia).at(ja);
// [OBSOLETE] 	fpos=a.fpos;
// [OBSOLETE] 	for(int ic=1;ic<=3;ic++) {
// [OBSOLETE] 	  fpos(ic)=fpos(ic)+atom_shifts(ia,ic);
// [OBSOLETE] 	}
// [OBSOLETE] 	cpos=F2C(lat,fpos);
// [OBSOLETE] 	a.cpos=cpos;
// [OBSOLETE] 	neigh_mat.at(ia).at(ja)=ConvertAtomToLat(a,original_lattice);
// [OBSOLETE]       }
// [OBSOLETE]     }
// [OBSOLETE]   }
// [OBSOLETE] }

// **************************************************************************
// Function GetStrNeighData
// **************************************************************************
// This function collects all the neighbor data out to some
//  cutoff and stores it for each atom in the structure.
namespace pflow {
  void GetStrNeighData(const xstructure& str, const double cutoff,
		       deque<deque<_atom> >& neigh_mat) {
    deque<_atom> atom_vec;
    neigh_mat.clear();
    // Get data from str.
    // Set scale to 1 so you don't need to rescale coordinates.
    xstructure sstr=str;
    sstr=ReScale(sstr,1.0);
    
    // Create atom objects for each atom in structure.
    xvector<int> ijk(3);ijk.clear();
    for(uint iat=0;iat<sstr.atoms.size();iat++) {
      _atom a=sstr.atoms.at(iat);
      a.name=sstr.atoms.at(iat).name;
      a.number=iat;
      a.ijk=sstr.atoms.at(iat).ijk;
      a.cpos=sstr.atoms.at(iat).cpos;
      a.fpos=sstr.atoms.at(iat).fpos;//cerr << sstr.atoms.at(iat).fpos << endl;
      a.type=sstr.atoms.at(iat).type;
      atom_vec.push_back(a);
    }
    double rmin=1e-6;
    // [OBSOLETE]    GetNeighData(atom_vec,sstr,rmin,cutoff,neigh_mat);
    sstr.GetNeighData(atom_vec,rmin,cutoff,neigh_mat);
  }
}

// **************************************************************************
// Function GetXrayScattFactor
// **************************************************************************
double GetXrayScattFactor(const string& name, const double& lambda) {
  if(lambda) {;} // phony just to keep lambda busy
  // Does not use lambda for now.
  double scatt_fact=0.0;
  for(uint iat=0;iat<atom_name_vec.size();iat++)
    if(name==atom_name_vec[iat] || name==atom_symbol_vec[iat]) scatt_fact=xray_scatt_vec[iat];
  return scatt_fact;
}

// ***************************************************************************
// Function GetVol
// ***************************************************************************
// This function returns the volume of the cell from the lattice parameters
namespace pflow {
  double GetVol(const xmatrix<double>& lat) {
    return (double) fabs(pflow::GetSignedVol(lat));
  }
  double GetVol(const matrix<double>& lat) {
    return (double) fabs(pflow::GetSignedVol(lat));
  }
}

// ***************************************************************************
// Function GetSignedVol
// ***************************************************************************
// This function returns the volume of the cell from the lattice parameters
namespace pflow {
  double GetSignedVol(const xmatrix<double>& lat) {
    double vol;
    xvector<double> u(3);
    u=aurostd::vector_product(lat(2),lat(3));
    vol=aurostd::scalar_product(lat(1),u);
    return vol;
  }
  double GetSignedVol(const matrix<double>& lat) {
    return (double) pflow::GetSignedVol(pflow::matrix2xmatrix(lat));
  }
}

// ***************************************************************************
// Function RecipLat xmatrix<double>
// ***************************************************************************
namespace pflow {
  xmatrix<double> RecipLat(const xmatrix<double>& lat) {
    xmatrix<double> rlat(3,3);
    xvector<double> rlat1(3),rlat2(3),rlat3(3);
    double vol = pflow::GetSignedVol(lat);
    rlat1=(2.0*PI/vol)*aurostd::vector_product(lat(2),lat(3));
    rlat2=(2.0*PI/vol)*aurostd::vector_product(lat(3),lat(1));
    rlat3=(2.0*PI/vol)*aurostd::vector_product(lat(1),lat(2));
    for(int i=1;i<=3;i++) rlat(1,i)=rlat1(i);
    for(int i=1;i<=3;i++) rlat(2,i)=rlat2(i);
    for(int i=1;i<=3;i++) rlat(3,i)=rlat3(i);
    return rlat;
  }
}

// ***************************************************************************
// Function RecipLat matrix<double>
// ***************************************************************************
namespace pflow {
  matrix<double> RecipLat(const matrix<double>& lat) {
    matrix<double> rlat (3,3);
    double vol = pflow::GetSignedVol(lat);
    rlat[0]=pflow::SVprod(2.0*PI/vol,pflow::VVcross(lat[1],lat[2]));
    rlat[1]=pflow::SVprod(2.0*PI/vol,pflow::VVcross(lat[2],lat[0]));
    rlat[2]=pflow::SVprod(2.0*PI/vol,pflow::VVcross(lat[0],lat[1]));
    return rlat;
  }
}

// ***************************************************************************
// Function SetCpos
// ***************************************************************************
namespace pflow {
  _atom SetCpos(const _atom& a, const vector<double>& in_cpos) {
    assert (in_cpos.size()==3);
    _atom b;b=a;
    b.cpos(1)=in_cpos.at(0);
    b.cpos(2)=in_cpos.at(1);
    b.cpos(3)=in_cpos.at(2);
    return b;
  }
}

// ***************************************************************************
// Function SetFpos
// ***************************************************************************
namespace pflow {
  _atom SetFpos(const _atom& a, const vector<double>& in_fpos) {
    assert (in_fpos.size()==3);
    _atom b;b=a;
    b.fpos(1)=in_fpos.at(0);
    b.fpos(2)=in_fpos.at(1);
    b.fpos(3)=in_fpos.at(2);
    return b;  
  }
}

// ***************************************************************************
// Function vecF2C
// ***************************************************************************
// This function converts a vector in direct to cartesian
namespace pflow {
  vector<double> vecF2C(const matrix<double>& lat, const vector<double>& vf) {
    vector<double> vc(3,0.0);
    for(int ic=0;ic<3;ic++) {
      for(int jc=0;jc<3;jc++)
	vc[ic]=vc[ic]+vf[jc]*lat[jc][ic];
    }
    return vc;
  }
}

// ***************************************************************************
// Function vecC2F
// ***************************************************************************
// This function converts a vector in direct to cartesian
namespace pflow {
  vector<double> vecC2F(const matrix<double>& lat, const vector<double>& vc) {
    vector<double> vf(3,0.0);
    vf=aurostd::xvector2vector(C2F(matrix2xmatrix(lat),aurostd::vector2xvector(vc)));
    return vf;
  }
}

// ***************************************************************************
// Function SetName
// ***************************************************************************
namespace pflow {
  _atom SetName(const _atom& a, const string& in_name) {
    _atom b;b=a;
    b.name=in_name;
    return b;
  }
}

// ***************************************************************************
// Function SetType
// ***************************************************************************
namespace pflow {
  _atom SetType(const _atom& a, const int in_type) {
    _atom b;b=a;
    b.type=in_type;    // CONVASP_MODE
    return b;
  }
}

// ***************************************************************************
// Function SetNum
// ***************************************************************************
namespace pflow {
  _atom SetNum(const _atom& a,const int in_num) {
    _atom b;b=a;
    b.number=in_num;
    return b;
  }
}

// **************************************************************************
// GetTypes
// **************************************************************************
namespace pflow {
  vector<int> GetTypes(const xstructure& a) {
    vector<int> out_type;
    for(uint iat=0;iat<a.atoms.size();iat++)
      out_type.push_back(a.atoms.at(iat).type);
    return out_type;
  }
}

// **************************************************************************
// GetNames
// **************************************************************************
namespace pflow {
  vector<string> GetNames(const xstructure& a) {
    vector<string> out_name;
    for(uint iat=0;iat<a.atoms.size();iat++)
      out_name.push_back(a.atoms.at(iat).name);
    return out_name;
  }
}

// **************************************************************************
// GetCleanNames
// **************************************************************************
namespace pflow {
  vector<string> GetCleanNames(const xstructure& a) {
    vector<string> out_cleanname;
    for(uint iat=0;iat<a.atoms.size();iat++)
      out_cleanname.push_back(a.atoms.at(iat).cleanname);
    return out_cleanname;
  }
}

// **************************************************************************
// GetSpins
// **************************************************************************
namespace pflow {
  vector<double> GetSpins(const xstructure& a) {
    vector<double> out_spin;
    for(uint iat=0;iat<a.atoms.size();iat++)
      out_spin.push_back(a.atoms.at(iat).spin);
  return out_spin;
  }
}

// ***************************************************************************
// Function GetFpos
// ***************************************************************************
namespace pflow {
  matrix<double> GetFpos(const xstructure& str) {
    int num_atoms=str.atoms.size();
    matrix<double> fpos(num_atoms,3);pflow::VVset(fpos,0.0);
    for(int i=0;i<num_atoms;i++)
      for(int j=0;j<3;j++)
	fpos[i][j]=str.atoms.at(i).fpos(j+1);
    return fpos;
  }
}

// ***************************************************************************
// Function GetCpos
// ***************************************************************************
namespace pflow {
  matrix<double> GetCpos(const xstructure& str) {
    int num_atoms=str.atoms.size();
    matrix<double> cpos(num_atoms,3);pflow::VVset(cpos,0.0);
    for(int i=0;i<num_atoms;i++)
      for(int j=0;j<3;j++)
	cpos[i][j]=str.atoms.at(i).cpos(j+1);
    return cpos;
  }
}

// ***************************************************************************
// Function SetLat
// ***************************************************************************
xstructure SetLat(const xstructure& a, const xmatrix<double>& in_lat) {
  xstructure b(a);
  b.lattice=in_lat;
  return b;
}
namespace pflow {
  xstructure SetLat(const xstructure& a, const matrix<double>& in_lat) {
    xstructure b(a);
    b.lattice=pflow::matrix2xmatrix(in_lat);
    return b;
  }
}

// ***************************************************************************
// Function GetLat
// ***************************************************************************
xmatrix<double> GetLat(const xstructure& a) {
  return a.lattice;
}
namespace pflow {
  matrix<double> GetLat(const xstructure& a) {
    return pflow::xmatrix2matrix(a.lattice);
  }
}

// ***************************************************************************
// Function GetScale
// ***************************************************************************
namespace pflow {
  double GetScale(const xstructure& a) {
    return a.scale;
  }
}

// ***************************************************************************
// Function GetScaledLat
// ***************************************************************************
namespace pflow {
  matrix<double> GetScaledLat(const xstructure& a) {
    return pflow::xmatrix2matrix(a.scale*a.lattice);
  }
}

// ***************************************************************************
// Function SetNumEachType
// ***************************************************************************
namespace pflow {
  xstructure SetNumEachType(const xstructure& a, const deque<int>& in_num_each_type) {
    xstructure b(a);
    b.num_each_type.clear();
    b.comp_each_type.clear();
    for(uint i=0;i<in_num_each_type.size();i++) {
      b.num_each_type.push_back(in_num_each_type.at(i));
      b.comp_each_type.push_back(in_num_each_type.at(i));
     }
    return b;
  }
}

// ***************************************************************************
// Function GetNumEachType
// ***************************************************************************
namespace pflow {
  deque<int> GetNumEachType(const xstructure& a) {
    return a.num_each_type;
  }
}

// ***************************************************************************
// Function AddAllAtomPos
// ***************************************************************************
namespace pflow {
  xstructure AddAllAtomPos(const xstructure& a,
			   const matrix<double>& in_pos,
			   const int in_coord_flag)  {
    assert(in_coord_flag==0 || in_coord_flag==1);
    xstructure b(a);
    b.atoms.clear();
    _atom atom;
    if(in_coord_flag==0) {
      for(uint i=0;i<in_pos.size();i++) {
	// atom=a.atoms.at(i);  // start from scratch
	for(int j=1;j<=3;j++) atom.fpos(j)=in_pos[i][j-1];
	atom.cpos=F2C(b.lattice,atom.fpos);
	b.atoms.push_back(atom);
      }
    }
    if(in_coord_flag==1) {
      for(uint i=0;i<in_pos.size();i++) {
	// atom=a.atoms.at(i);  // start from scratch
	for(int j=1;j<=3;j++) atom.cpos(j)=in_pos[i][j-1];
	atom.fpos=C2F(b.lattice,atom.cpos);
	b.atoms.push_back(atom);
      }
    }
    return b;
  }
}

// ***************************************************************************
// Function SetAllAtomPos
// ***************************************************************************
namespace pflow {
  xstructure SetAllAtomPos(const xstructure& a,
			   const matrix<double>& in_pos,
			   const int in_coord_flag)  {
    assert(in_coord_flag==0 || in_coord_flag==1);
    xstructure b(a);
    b.atoms.clear();
    _atom atom;
    if(in_coord_flag==0) {
      for(uint i=0;i<in_pos.size();i++) {
	atom=a.atoms.at(i);
	for(int j=1;j<=3;j++) atom.fpos(j)=in_pos[i][j-1];
	atom.cpos=F2C(b.lattice,atom.fpos);
	b.atoms.push_back(atom);
      }
    }
    if(in_coord_flag==1) {
      for(uint i=0;i<in_pos.size();i++) {
	atom=a.atoms.at(i);
	for(int j=1;j<=3;j++) atom.cpos(j)=in_pos[i][j-1];
	atom.fpos=C2F(b.lattice,atom.cpos);
	b.atoms.push_back(atom);
      }
    }
    return b;
  }
}

// ***************************************************************************
// Function SetAllAtomNames
// ***************************************************************************
namespace pflow {
  xstructure SetAllAtomNames(const xstructure& a, const vector<string>& in) {
    xstructure b(a);
    if(in.size()==a.num_each_type.size()) {
      cerr << "SetAllAtomNames: this routine must be fixed... it does not work here" << endl;
      exit(0);
      
      for(uint iat=0;iat<b.num_each_type.size();iat++) {
	b.atoms.at(iat).name=in.at(b.atoms.at(iat).type);       // CONVASP_MODE
	if(b.atoms.at(iat).name.size()) {
	  b.atoms.at(iat).name_is_given=TRUE;
	  b.atoms.at(iat).CleanName();
	  //DX 9/21/17 - Need to keep spin info b.atoms.at(iat).CleanSpin();
	}
      }
      return b;
    }
    if(in.size()==a.atoms.size()) {
      for(uint iat=0;iat<b.atoms.size();iat++) {
	b.atoms.at(iat).name=in.at(iat);
	b.atoms.at(iat).name_is_given=TRUE;
	b.atoms.at(iat).CleanName();
	//DX 9/21/17 - Need to keep spin info b.atoms.at(iat).CleanSpin();
      }
      return b;
    }
    cerr << "Must specity as many names as types/numbers: in.size()=" << (int) in.size()
	 << "   =a.num_each_type.size()=" << (int) a.num_each_type.size()
	 << "   =a.atoms.size()=" << (int) a.atoms.size() << endl;
    exit(0);
  }
}

// ***************************************************************************
// Function SetNamesWereGiven
// ***************************************************************************
namespace pflow {
  xstructure SetNamesWereGiven(const xstructure& a, const vector<int>& in) {
    xstructure b(a);
    if(in.size()==a.num_each_type.size()) {
      cerr << "ERROR - SetNamesWereGiven: this routine must be fixed... it does not work here" << endl;
      exit(0);
      for(uint iat=0;iat<b.num_each_type.size();iat++)
	b.atoms.at(iat).name_is_given=in.at(b.atoms.at(iat).type);      // CONVASP_MODE
      return b;
    }
    if(in.size()==a.atoms.size()) {
      for(uint iat=0;iat<b.atoms.size();iat++) {
	b.atoms.at(iat).name_is_given=in.at(iat);
      }
      return b;
    }
    cerr << "Must specity as many names as types/numbers: in.size()=" << (uint) in.size()
	 << "   =a.num_each_type.size()=" << (uint) a.num_each_type.size()
	 << "   =a.atoms.size()=" << (uint) a.atoms.size() << endl;
    exit(0);
  }
}

// ***************************************************************************
// Function SetOrigin
// ***************************************************************************
namespace pflow {
  xstructure SetOrigin(const xstructure& a, const vector<double>& in_origin) {
    xstructure b(a);
    for(uint i=1;i<=3;i++) b.origin[i]=in_origin.at(i-1);
    return b;
  }
  xstructure SetOrigin(const xstructure& a, const xvector<double>& in_origin) {
    xstructure b(a);
    for(uint i=1;i<=3;i++) b.origin[i]=in_origin[i];
    return b;
  }
}

// ***************************************************************************
// Function VVequal
// ***************************************************************************
// This function checks if two vectors are equal. double/int
namespace pflow {
  // vector
  bool VVequal(const vector<double>& a, const vector<double>& b) {
    double tol=1e-15;
    int size=std::min((int) a.size(),(int) b.size());
    for(int i=0;i<size;i++) if(fabs(a.at(i)-b.at(i))>tol) return FALSE;
    return TRUE;
  }
  bool VVequal(const vector<int>& a, const vector<int>& b) {
    double tol=1e-15;
    int size=std::min((int) a.size(),(int) b.size());
    for(int i=0;i<size;i++) if(fabs((double)(a.at(i)-b.at(i)))>tol) return FALSE;
    return TRUE;
  }
  // deque
  bool VVequal(const deque<double>& a, const deque<double>& b) {
    double tol=1e-15;
    int size=std::min((int) a.size(),(int) b.size());
    for(int i=0;i<size;i++) if(fabs(a.at(i)-b.at(i))>tol) return FALSE;
    return TRUE;
  }
  bool VVequal(const deque<int>& a, const deque<int>& b) {
    double tol=1e-15;
    int size=std::min((int) a.size(),(int) b.size());
    for(int i=0;i<size;i++) if(fabs((double)(a.at(i)-b.at(i)))>tol) return FALSE;
    return TRUE;
  }
  
}
// ***************************************************************************
// Function SmoothFunc
// ***************************************************************************
// This function smooths an input functions.
// Uses gaussian smoothing for now.
// Assumes function argument is integer given by vector id.
// Sigma is in units of function arguments (i.e., vector id).
// Dane Morgan
namespace pflow {
  vector<double> SmoothFunc(const vector<double>& func, const double& sigma) {
    if(sigma<1e-10) return func;
    int nx=func.size();
    int range=(int) (5.0*sigma);
    // Get Gaussian weights
    vector<double> norm(nx,0.0);
    matrix<double> wt(nx,2*range+1);pflow::VVset(wt,0.0);

    for(int ix=0;ix<nx;ix++) {
      for(int i=-range;i<=range;i++) {
	if((ix+i)>=0 && (ix+i)<nx) {
	  wt[ix][i+range]=Normal((double)(ix+i),double(ix),sigma);
	  norm[ix]=norm[ix]+wt[ix][i+range];
	}
      }
    }
    // Normalize to one
    for(int ix=0;ix<nx;ix++) {
      for(int i=-range;i<=range;i++) {
	wt[ix][i+range]=wt[ix][i+range]/norm[ix];
      }
    }

    vector<double> sfunc(nx,0.0);
    // Average in weighted nearby bins.
    for(int ix=0;ix<nx;ix++) {
      double sf=0;
      for(int i=-range;i<=range;i++) {
	if((ix+i)>0 && (ix+i)<nx) {
	  sf+=wt[ix][i+range]*func[ix+i];
	}
      }
      sfunc[ix]=sf;
    }
    return sfunc;
  }
}
// ***************************************************************************
//  Function Normal
// ***************************************************************************
// This function returns the value of a normal distribution.
double Normal(const double& x, const double& mu, const double& sigma) {
  double tol=1e-12;
  // If sigma=0 return delta function
  if(abs(sigma)<tol) {
    if(abs(x-mu)<tol) {
      return 1;
    }
    else{
      return 0;
    }
  }
  double arg = -(x-mu)*(x-mu)/(2*sigma*sigma);
  return exp(arg)/(sqrt(TWOPI)*sigma);
}

// ***************************************************************************
// Function Set
// ***************************************************************************
// This function sets a vector of vector with the proper values
namespace pflow {
  void VVset(matrix<double> &mat,const double& value) {
    for(uint i=0;i<mat.size();i++)
      for(uint j=0;j<mat[i].size();j++)
	mat[i][j]=(double) value;
  }
  void VVset(vector<vector<int> > &mat,const int& value) {
    for(uint i=0;i<mat.size();i++)
      for(uint j=0;j<mat[i].size();j++)
	mat[i][j]=(int) value;
  }
}

// ***************************************************************************
// Function norm
// ***************************************************************************
// This function finds the norm of a vector
// Dane Morgan style
namespace pflow {
  //  template<class utype>
  double norm(const vector<double>& v) {
    double sum=0;
    for(uint i=0;i<v.size();i++) {
      sum=sum+v[i]*v[i];
    }
    //    sum=sqrt(fabs(sum));
    sum=sqrt(aurostd::abs(sum));
    return sum;
  }
}

// ***************************************************************************
// Function getcos
// ***************************************************************************
// This function finds the cosine of the angle between
// two vectors.
// Dane Morgan style
namespace pflow {
  // template<class utype>
  double getcos(const vector<double>& a, const vector<double>& b) {
    double sum=0.0;
    uint size=std::min((uint) a.size(),(uint) b.size());
    for(uint i=0;i<size;i++) sum=sum+a[i]*b[i];
    double na=pflow::norm(a);
    double nb=pflow::norm(b);
    assert(na>0.0 && nb>0.0);
    sum=sum/(na*nb);
    return sum;
  }
}

// ***************************************************************************
// Function Getabc_angles
// ***************************************************************************
// This function returns a,b,c,alpha,beta,gamma for the cell
// given the lattice vectors.
// Dane Morgan style
namespace pflow {
  // template<class utype>
  vector<double> Getabc_angles(const matrix<double>& lat) {
    //    cerr << lat[0][0] << " " << lat[0][1] << " " << lat[0][2] << endl;
    //    cerr << lat[1][0] << " " << lat[1][1] << " " << lat[1][2] << endl;
    //   cerr << lat[2][0] << " " << lat[2][1] << " " << lat[2][2] << endl;
    vector<double> data;
    data.push_back(pflow::norm(lat[0])); // cerr << data.at(data.size()-1) << endl;
    data.push_back(pflow::norm(lat[1])); // cerr << data.at(data.size()-1) << endl;
    data.push_back(pflow::norm(lat[2])); // cerr << data.at(data.size()-1) << endl;
    data.push_back(acos(pflow::getcos(lat[1],lat[2]))); // cerr << data.at(data.size()-1) << endl;
    data.push_back(acos(pflow::getcos(lat[0],lat[2]))); // cerr << data.at(data.size()-1) << endl;
    data.push_back(acos(pflow::getcos(lat[0],lat[1]))); // cerr << data.at(data.size()-1) << endl;
    return data;
  }
}

// namespace pflow {
//   void dont_run_this(void) {
//     {
//       vector<float> d;pflow::norm(d);pflow::getcos(d,d);
//       matrix<float> D;pflow::Getabc_angles(D);
//     }
//     {
//       vector<double> d;pflow::norm(d);pflow::getcos(d,d);
//       matrix<double> D;pflow::Getabc_angles(D);
//     }
//     {
//       vector<long double> d;pflow::norm(d);pflow::getcos(d,d);
//       matrix<long double> D;pflow::Getabc_angles(D);
//     }
//   }
// }

// ***************************************************************************
// Function Sort_abc_angles
// ***************************************************************************
// This fucntions sorts lattice vectors to a>=b>=c and
// arranges angles equivalently.
// Dane Morgan style
namespace pflow {
  vector<double> Sort_abc_angles(const vector<double>& abc_angles) {
    vector<double> old = abc_angles;
    vector<double> newv = abc_angles;
    // if b>a then swap a/b: a=b,b=a,c=c,alpha=beta,beta=alpha,gamma=gamma.
    if(old[1]>old[0]) {
      newv[0]=old[1];
      newv[1]=old[0];
      newv[3]=old[4];
      newv[4]=old[3];
      old=newv;
    }
    // if c>a then swap a/c: a=c,b=b,c=a,alpha=gamma,beta=beta,gamma=alpha.
    if(old[2]>old[0]) {
      newv[0]=old[2];
      newv[2]=old[0];
      newv[3]=old[5];
      newv[5]=old[3];
      old=newv;
    }
    // if c>b then swap b/c: a=a,b=c,c=b,alpha=alpha,beta=gamma,gamma=beta.
    if(old[2]>old[1]) {
      newv[1]=old[2];
      newv[2]=old[1];
      newv[4]=old[5];
      newv[5]=old[4];
      old=newv;
    }
    return newv;
  }
}

// ***************************************************************************
// Function Vout
// ***************************************************************************
// This function outputs a vector.
// Dane Morgan style
namespace pflow {
  void Vout(const vector<double>& a, ostream& out) {
    for(uint i=0;i<a.size();i++) out << a.at(i) << " ";  
    out << endl;
  }
  void Vout(const vector<int>& a, ostream& out) {
    for(uint i=0;i<a.size();i++) out << a.at(i) << " ";  
    out << endl;
  }
  void Vout(const vector<string>& a, ostream& out) {
    for(uint i=0;i<a.size();i++) out << a.at(i) << " ";  
    out << endl;
  }
}

// ***************************************************************************
//  Function Mout
// ***************************************************************************
// This function outputs a matrix.
// Dane Morgan style
namespace pflow {
  void Mout(const matrix<double>& m, ostream& out) {
    for(uint i=0;i<m.size();i++) {
      for(uint j=0;j<m[i].size();j++) {
	out << m[i][j] << " ";
      }
      out << endl;
    }
  }
  void Mout(const vector<vector<double> >& m, ostream& out) {
    for(uint i=0;i<m.size();i++) {
      for(uint j=0;j<m[i].size();j++) {
	out << m[i][j] << " ";
      }
      out << endl;
    }
  }
}

// **************************************************
//  Function SVprod
// **************************************************
// This function returns the product of a scalar and a vector.
namespace pflow {
  vector<double> SVprod(const double& s, const vector<double>& b) {
    int size=b.size();
    vector<double> prod(size);
    for(int i=0;i<size;i++) {
      prod[i]=s*b[i];
    }
    return prod;
  }
  vector<int> SVprod(const int& s, const vector<int>& b) {
    int size=b.size();
    vector<int> prod(size);
    for(int i=0;i<size;i++) {
      prod[i]=s*b[i];
    }
    return prod;
  }
}

// **************************************************
//  Function VVsum
// **************************************************
// This function returns the vector sum c=a+b.
namespace pflow {
  vector<double> VVsum(const vector<double>& a, const vector<double>& b) {
    int size=std::min(a.size(),b.size());
    vector<double> c(size);
    for(int i=0;i<size;i++) {
      c[i]=a[i]+b[i];
    }
    return c;
  }
  vector<double> VVsum(const vector<double>& a, const vector<int>& b) {
    int size=std::min(a.size(),b.size());
    vector<double> c(size);
    for(int i=0;i<size;i++) {
      c[i]=a[i]+b[i];
    }
    return c;
  }
}

// **************************************************
// Function VVdiff
// **************************************************
// This function returns the vector c=a-b.
namespace pflow {
  vector<double> VVdiff(const vector<double>& a, const vector<double>& b) {
    int size=std::min(a.size(),b.size());
    vector<double> c(size);
    for(int i=0;i<size;i++) {
      c[i]=a[i]-b[i];
    }
    return c;
  }
}

// **************************************************
//  Function VVprod
// **************************************************
// This function returns the scalar c=a*b.
// Dane Morgan style
namespace pflow {
  double VVprod(const vector<double>& a, const vector<double>& b) {
    int size=std::min(a.size(),b.size());
    double c=0;
    for(int i=0;i<size;i++) {
      c=c+a[i]*b[i];
    }
    return c;
  }
  double VVprod(const vector<double>& a, const vector<int>& b)
  {
    int size=std::min(a.size(),b.size());
    double c=0;
    for(int i=0;i<size;i++) {
      c=c+a[i]*b[i];
    }
    return c;
  }
}

// ***************************************************************************
// Function MMmult
// ***************************************************************************
// This function returns the product of two matrices.
// a=MxN, b=NxM.
// Dane Morgan style
namespace pflow {
  matrix<double> MMmult(const matrix<double>& a, const matrix<double>& b) {
    uint M=std::min((uint) a.size(),(uint) b[0].size());
    uint N=std::min((uint) a[0].size(),(uint) b.size());
    vector<double> v(M,0.0);
    matrix<double> c(M,v);
    for(uint i=0;i<M;i++) {
      for(uint j=0;j<M;j++) {
	double sum=0.0;
	for(uint k=0;k<N;k++) {
	  sum=sum+a[i][k]*b[k][j];
	}
	c[i][j]=sum;
      }
    }
    return c;
  }
}

// ***************************************************************************
// Function MVmult
// ***************************************************************************
//  This function returns the product of a matric and a vector A*v.
//  A=MxN, v=Nx1.
// Dane Morgan style
namespace pflow {
  vector<double> MVmult(const matrix<double>& A, const vector<double>& v) {
    uint M=A.size();
    uint N=std::min((uint) A[0].size(),(uint) v.size());
    vector<double> u(M,0.0);
    for(uint i=0;i<M;i++) {
      double sum=0;
      for(uint j=0;j<N;j++) {
	sum=sum+A[i][j]*v[j];
      }
      u[i]=sum;
    }
    return u;
  }
}

// ***************************************************************************
// Function VMmult
// ***************************************************************************
// This function returns the product of a row vector and a matrix v*A.
// A=MxN, v=1xM. Return vector of length N.
// Dane Morgan style
namespace pflow {
  vector<double> VMmult(const vector<double>& v, const matrix<double>& A) {
    uint M=A.size();
    uint N=std::min((uint) A[0].size(),(uint) v.size());
    vector<double> u(M,0.0);
    for(uint i=0;i<N;i++) {
      double sum=0;
      for(uint j=0;j<M;j++)
	sum=sum+v[j]*A[j][i];
      u[i]=sum;
    }
    return u;
  }
  vector<double> VMmult(const vector<int>& v, const matrix<double>& A) {
    uint M=A.size();
    uint N=std::min((uint) A[0].size(),(uint) v.size());
    vector<double> u(M,0.0);
    for(uint i=0;i<N;i++) {
      double sum=0;
      for(uint j=0;j<M;j++)
	sum=sum+v[j]*A[j][i];
      u[i]=sum;
    }
    return u;
  }
}

// ***************************************************************************
// Function VVcross
// ***************************************************************************
// This function returns the cross product of two 3 dimensional vectors. */
// Dane Morgan style
namespace pflow {
  vector<double> VVcross(const vector<double>& a, const vector<double>& b) {
    if(a.size()!=3 || b.size()!=3) return a;
    vector<double> u(3);
    u[0]=a[1]*b[2]-a[2]*b[1];
    u[1]=-(a[0]*b[2]-a[2]*b[0]);
    u[2]=a[0]*b[1]-a[1]*b[0];
    return u;
  }
}

// ***************************************************************************
// Function VVdot
// ***************************************************************************
// This function returns the dot product of two vectors. */
// Dane Morgan style
namespace pflow {
  double VVdot(const vector<double>& a, const vector<double>& b) {
    double sum=0;
    int size=std::min(a.size(),b.size());
    for(int i=0;i<size;i++) {
      sum=sum+a[i]*b[i];
    }
    return sum;
  }
}

// ***************************************************************************
// Function GetNumAtoms
// ***************************************************************************
// This function returns the number of atoms in a structure
namespace pflow {
  int GetNumAtoms(const xstructure& a) {
    return a.atoms.size();
  }
}

// ***************************************************************************
// Function SetSpline
// ***************************************************************************
// This function finds the second derivates so that
// spline interpolations can be evaluated quickly.
// I have shifted all indices to go from 0->(n-1).
namespace pflow {
  void SetSpline(const vector<double>& x, const vector<double>& y,
		 const double& yp1, const double& ypn, vector<double>& y2) {
    int i,k;
    double p,qn,sig,un;
    int n=x.size();

    vector<double> u(n-1);
    y2=vector<double> (n,0.0);

    if(yp1 > 0.99e30)
      y2[0]=u[0]=0.0;
    else {
      y2[0] = -0.5;
      u[0]=(3.0/(x[2]-x[1]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
    }
    for(i=1;i<=n-2;i++) {
      sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
      p=sig*y2[i-1]+2.0;
      y2[i]=(sig-1.0)/p;
      u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
      u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
    }
    if(ypn > 0.99e30)
      qn=un=0.0;
    else {
      qn=0.5;
      un=(3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
    }
    y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0);
    for(k=n-2;k>=0;k--) {
      y2[k]=y2[k]*y2[k+1]+u[k];
    }
  }
}

// ***************************************************************************
// Function GetSplineInt
// ***************************************************************************
// This function finds the value of a function y at a
// point x from a spline interpolation.  You must run
// SetSpline first.
// I have shifted all indices to go from 0->(n-1).
namespace pflow {
  void GetSplineInt(const vector<double>& xa, const vector<double>& ya,
		    vector<double>& y2a, const double& x, double& y) {
    int klo,khi,k;
    double h,b,a;
    int n=xa.size();
    
    klo=0;
    khi=n-1;
    while (khi-klo > 1) {
      k=(khi+klo) >> 1;
      if(xa[k] > x) khi=k;
      else klo=k;
    }
    h=xa[khi]-xa[klo];
    if(h == 0.0) {
      cerr << "ERROR - pflow::GetSplineInt Bad xa input to routine splint" << endl;
      exit(0);
    }
    a=(xa[khi]-x)/h;
    b=(x-xa[klo])/h;
    y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
  }
}

// ***************************************************************************
// Function PrintSpline
// ***************************************************************************
// This function prints out the original data and the
// spline interpolation with keywords for grep.
namespace pflow {
  void PrintSpline(const vector<double>& x, const vector<double>& y,
		   const int& npts, ostream& outf) {
    outf.setf(std::ios::left, std::ios::adjustfield);
    outf.setf(std::ios::fixed, std::ios::floatfield);
    outf << "******************** Initial Data ********************" << endl;
    int ndat=x.size();
    for(int id=0;id<ndat;id++) {
      outf<< x[id] << " " << y[id]<< "   " << "INPUT" << endl;
    }
    outf << "******************** Cubic Spline Interpolated Points ********************" << endl;
    if(npts==1) {
      outf<< x[0] << " " << y[0] << "   " << "INTERPOLATED" << endl;
    }
    else{
      double yp1=0.0;
      double ypn=0.0;
      vector<double> y2;
      SetSpline(x,y,yp1,ypn,y2);
      double xrange=x[ndat-1]-x[0];
      double dx=xrange/(double)(npts-1);
      for(int ip=0;ip<npts;ip++) {
	double xp=x[0]+(double)ip*dx;
	double yp;
	GetSplineInt(x,y,y2,xp,yp);
	outf << xp << " " << yp << "   " << "INTERPOLATED" << endl;
      }
    }
  }
}

// ***************************************************************************
// XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
// ***************************************************************************

bool never_call_this_function(void) {
  { pflow::matrix<double> m0(1,1,(double)0),m1(1,1),m2(1),m3;
    std::vector<double> v(1);
    m1.push_back(v);m1[0];m1[0][0];m2=m1; }
  { pflow::matrix<int> m0(1,1,(int)0),m1(1,1),m2(1),m3;
    std::vector<int> v(1);
    m1.push_back(v);m1[0];m1[0][0];m2=m1; }
  { pflow::matrix<std::complex<double> > m0(1,1,(std::complex<double>)0.0),m1(1,1),m2(1),m3;
    std::vector<std::complex<double> > v(1);
    m1.push_back(v);m1[0];m1[0][0];m2=m1;}
  { pflow::matrix<double> m0(1,1,(double) 0),m1(1,1),m2(1),m3;
    pflow::matrix<pflow::matrix<double> > mm1(1,1,m0),mm2(1,1,m0);
    mm2=mm1;}
  { pflow::matrix<std::complex<double> > m0(1,1,(std::complex<double>) 0),m1(1,1),m2(1),m3;
    pflow::matrix<pflow::matrix<std::complex<double> > > mm1(1,1,m0),mm2(1,1,m0);
    mm2=mm1;}
  return TRUE;
}

// ***************************************************************************

#endif

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *                                                                         *
// ***************************************************************************
