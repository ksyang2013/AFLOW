// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2008           *
// *                                                                         *
// ***************************************************************************
// Written by Stefano Curtarolo 1994-2008
// Components by Gus Hart 2008

#ifndef _XGUS_CPP_
#define _XGUS_CPP_

#include "aflow.h"

// all vector/matrix must be loaded before

// ----------------------------------------------------------------------------
// ------------------------------------------------- MODULE numerical_utilities

using aurostd::xvector;
using aurostd::xmatrix;
using aurostd::xtensor3;

#define _DEFAULT_EPS_BASIS_REDUCE_ 0.001

namespace gusstd {
  // namespace gusstd
  //****************************************************************************************
  // This function takes two real entities and compares them to see if they are equal
  // within some tolerance. This prevents failed comparisons due to numbers that are "equal"
  // but differ due to small differences arising from finite precision.
  
  bool equal_rank1(const xvector<double>& a, const xvector<double>& b, const double& tolerance) {
    if(a.rows!=b.rows) {cerr << _AUROSTD_XLIBS_ERROR_ << "AFLOW_XGUS.CPP: failure equal_rank1 [1]" << endl;exit(0);}
    for(register int i=a.lrows,ii=b.lrows;i<=a.urows;i++,ii++)
      if(abs(a(i)-b(ii)) > tolerance) return static_cast<bool>(FALSE);    
    return static_cast<bool>(TRUE);
  }
  
  bool equal_rank2(const xmatrix<double>& a, const xmatrix<double>& b, const double& tolerance) {
    if(a.rows!=b.rows) {cerr << _AUROSTD_XLIBS_ERROR_ << "AFLOW_XGUS.CPP: failure equal_rank2 [1]" << endl;exit(0);}
    if(a.cols!=b.cols) {cerr << _AUROSTD_XLIBS_ERROR_ << "AFLOW_XGUS.CPP: failure equal_rank2 [2]" << endl;exit(0);}
    for(register int i=a.lrows,ii=b.lrows;i<=a.urows;i++,ii++)
      for(register int j=a.lcols,jj=b.lcols;j<=a.ucols;j++,jj++)
	if(abs(a(i,j)-b(ii,jj)) > tolerance) return static_cast<bool>(FALSE);    
    return static_cast<bool>(TRUE);
  }
  
  bool equal_rank3(const xtensor3<double>& a, const xtensor3<double>& b, const double& tolerance) {
    if(a.index[1]!=b.index[1]) {cerr << _AUROSTD_XLIBS_ERROR_ << "AFLOW_XGUS.CPP: failure equal_rank3 [1]" << endl;exit(0);}
    if(a.index[2]!=b.index[2]) {cerr << _AUROSTD_XLIBS_ERROR_ << "AFLOW_XGUS.CPP: failure equal_rank3 [2]" << endl;exit(0);}
    if(a.index[3]!=b.index[3]) {cerr << _AUROSTD_XLIBS_ERROR_ << "AFLOW_XGUS.CPP: failure equal_rank3 [3]" << endl;exit(0);}
    for(register int i=a.lindex[1],ii=b.lindex[1];i<=a.uindex[1];i++,ii++)
      for(register int j=a.lindex[2],jj=b.lindex[2];j<=a.uindex[2];j++,jj++)
	for(register int k=a.lindex[3],kk=b.lindex[3];k<=a.uindex[3];k++,kk++)
	  if(abs(a(i,j,k)-b(ii,jj,kk)) > tolerance) return static_cast<bool>(FALSE);    
    return static_cast<bool>(TRUE);
  }
  
  bool equal_scalar(const double& a, const double& b, const double& tolerance) {
    if(abs(a-b) > tolerance) return static_cast<bool>(FALSE);    
    return static_cast<bool>(TRUE);
  }

  bool equal_rank0(const double& a, const double& b, const double& tolerance) {
    if(abs(a-b) > tolerance) return static_cast<bool>(FALSE);    
    return static_cast<bool>(TRUE);
  }

  bool equal_rank1_rank0(const xvector<double>& a, const double& b, const double& tolerance) {
    for(register int i=a.lrows;i<=a.urows;i++)
      if(abs(a(i)-b) > tolerance) return static_cast<bool>(FALSE);    
    return static_cast<bool>(TRUE);
  }

  bool equal_rank2_rank0(const xmatrix<double>& a, const double& b, const double& tolerance) {
    for(register int i=a.lrows;i<=a.urows;i++)
      for(register int j=a.lcols;j<=a.ucols;j++)
	if(abs(a(i,j)-b) > tolerance) return static_cast<bool>(FALSE);    
    return static_cast<bool>(TRUE);
  }
  
  bool equal_rank1_real_int(const xvector<double>& a, const xvector<int>& b, const double& tolerance) {
    if(a.rows!=b.rows) {cerr << _AUROSTD_XLIBS_ERROR_ << "AFLOW_XGUS.CPP: failure equal_rank1_real_int [1]" << endl;exit(0);}
    for(register int i=a.lrows,ii=b.lrows;i<=a.urows;i++,ii++)
      if(abs(a(i)-((double) b(ii))) > tolerance) return static_cast<bool>(FALSE);    
    return static_cast<bool>(TRUE);
  }
  bool equal_rank2_real_int(const xmatrix<double>& a, const xmatrix<int>& b, const double& tolerance) {
    if(a.rows!=b.rows) {cerr << _AUROSTD_XLIBS_ERROR_ << "AFLOW_XGUS.CPP: failure equal_rank2_real_int [1]" << endl;exit(0);}
    if(a.cols!=b.cols) {cerr << _AUROSTD_XLIBS_ERROR_ << "AFLOW_XGUS.CPP: failure equal_rank2_real_int [2]" << endl;exit(0);}
    for(register int i=a.lrows,ii=b.lrows;i<=a.urows;i++,ii++)
      for(register int j=a.lcols,jj=b.lcols;j<=a.ucols;j++,jj++)
	if(abs(a(i,j)-((double) b(ii,jj))) > tolerance) return static_cast<bool>(FALSE);    
    return static_cast<bool>(TRUE);
  }
}

// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// --------------------------------------------- MODULE vector_matrix_utilities

namespace gusstd {
  // namespace gusstd
  
  // ***************************************************************************************************
  //  This function calculates the "orthogonality defect" of the given basis of a 3D lattice.

  double orthogonality_defect(const xmatrix<double>& basis) {
    double od=1.0;
    xvector<double> bj(3);
    for(int j=1;j<=3;j++) {
      bj(1)=basis(1,j);
      bj(2)=basis(2,j);
      bj(3)=basis(3,j);
      od=od*aurostd::modulus(bj);
    }
    od=od/abs(determinant(basis));
    return od;
  }

  
  // ***************************************************************************************************
  //  This routine takes a set of basis vectors (that form a lattice) and reduces them so that they form
  //  the shortest possible basis. See Lecture notes in computer science, ISSN 0302-974, ANTS - VI : algorithmic
  //  number theory, 2004, vol. 3076, pp. 338-357 ISBN 3-540-22156-5
  double reduce_to_shortest_basis(const xmatrix<double>& IN,xmatrix<double>& OUT,double eps,bool VERBOSE) {
    xvector<double> A(3),B(3),C(3);
    xmatrix<double> check(3,3);
    bool err;
    double od,odnew;
    int ii=0,iimax=10000;
    // IN has colum-vectors
    for(int i=1;i<=3;i++) {
      A(i)=IN(i,1); // 1st vector
      B(i)=IN(i,2); // 2nd vector
      C(i)=IN(i,3); // 3rd vector
    }
    odnew=orthogonality_defect(IN);
    if(VERBOSE) cout << "AFLOW_XGUS.CPP: (reduce_to_shortest_basis) Before reduction, the orthogonality defect of the basis was " << odnew << endl;
    bool goexit=FALSE;
    while(goexit==FALSE) {
      od=odnew;
      reduce_A_in_ABC(A,B,C,eps);
      reduce_A_in_ABC(B,C,A,eps);
      reduce_A_in_ABC(C,A,B,eps);
      for(int i=1;i<=3;i++) {
	OUT(i,1)=A(i);OUT(i,2)=B(i);OUT(i,3)=C(i);
      }
      odnew=orthogonality_defect(OUT);
      // write(*,'("OD: ",2(f7.3,1x))') odnew, od
      //      cerr << od << " " << odnew << endl;
      if(abs(od-odnew)<eps) goexit=TRUE;
      if(ii++>iimax) goexit=TRUE;
      //     if(ii++>iimax) {OUT=IN;return orthogonality_defect(OUT);};
    }
    matrix_inverse(OUT,check,err);
    if(err==TRUE) { 
      cerr << _AUROSTD_XLIBS_ERROR_ << "AFLOW_XGUS.CPP: (reduce_to_shortest_basis) OUT matrix is singular in reduce_to_shortest_basis" << endl;
      exit(0);
    }
    //  Check that the conversion from old to new lattice vectors is still an integer matrix
    if(sum(abs(check*IN-nint(check*IN)))>eps) { 
      cerr << _AUROSTD_XLIBS_ERROR_ << "AFLOW_XGUS.CPP: (reduce_to_shortest_basis) ERROR: Reduced lattice vectors in reduce_to_shortest_basis changed the original lattice" << endl;
      exit(0);
    }
    if(VERBOSE) cout << "AFLOW_XGUS.CPP: (reduce_to_shortest_basis) After reduction, the orthogonality defect of the basis is " << orthogonality_defect(OUT) << endl;
    //GH if we have a left-handed basis, then exchange two vectors so that the basis is right-handed (I don't care but VASP does...Why?)
    if(determinant(OUT)<eps) {
      double temp;
      for(int i=1;i<=3;i++) {temp=OUT(i,1);OUT(i,1)=OUT(i,2);OUT(i,2)=temp;} // swap 1st with 2nd vector
    }
    // OUT has colum-vectors
    return orthogonality_defect(OUT);
  }

  xmatrix<double> reduce_to_shortest_basis(const xmatrix<double>& IN,double eps,bool VERBOSE) {
    xmatrix<double> newbasis(3,3);
    reduce_to_shortest_basis(IN,newbasis,eps,VERBOSE);
    return newbasis;
  }
  
  xmatrix<double> reduce_to_shortest_basis(const xmatrix<double>& IN,bool VERBOSE) {
    return reduce_to_shortest_basis(IN,_DEFAULT_EPS_BASIS_REDUCE_,VERBOSE);
  }
  
  xmatrix<double> reduce_to_shortest_basis(const xmatrix<double>& IN,double eps) {
    xmatrix<double> newbasis(3,3);
    reduce_to_shortest_basis(IN,newbasis,eps,FALSE);
    return newbasis;
  }
  
  xmatrix<double> reduce_to_shortest_basis(const xmatrix<double>& IN) {
    return reduce_to_shortest_basis(IN,_DEFAULT_EPS_BASIS_REDUCE_,FALSE);
  }
  

  
  // ***************************************************************************************************
  // This routine takes three vectors, A,B,C, defining a lattice, and reduces the first one so that it
  // is a close as possible to the origin while remaining in the affine plane which is defined by B,C but
  // shifted by A. See Lecture notes in computer science, ISSN 0302-974, ANTS - VI : algorithmic
  // number theory, 2004, vol. 3076, pp. 338-357 ISBN 3-540-22156-5
  void reduce_A_in_ABC(xvector<double>& A, xvector<double>& B, xvector<double>& C,double eps) {
    xvector<double> T(3);  // closest point to origin in B,C+A affine plane
    xmatrix<double> ABC(3,3),ABCinv(3,3),oldABC(3,3); // Matrices of ABC basis vectors and inverse
    xvector<double> dist(4); // the distances from T to enclosing lattice points of B,C (4 corners of the ppiped)
    xvector<double> i(3),i1(3),i2(3),i3(3),i4(3); // lattice coordinates of A, in the affine plane, using the B,C basis vectors
    int idx; // index of the smallest distance from T to a lattice point in B,C
    bool err;
    //integer j
    double lambda;
    //print *,"entering reduction routine..."
    //write(*,'("aurostd::modulus(A): ",f7.3,5x," A ",3(f7.3,1x)A)') aurostd::modulus(A), A
    for(int i=1;i<=3;i++) {
      ABC(i,1)=A(i);ABC(i,2)=B(i);ABC(i,3)=C(i);
    }
    oldABC=ABC;
    // Use Gaussian reduction to reduce the B,C 2D basis so that it is itself Minkowski reduced. If this
    // is done then the closest lattice point (in B,C plane) to projection of A (into the B,C plane) is
    // guaranteed to be one of the corners of the unit cell enclosing the projection of A
    gaussian_reduce_two_vectors(B,C,eps);

    //do j=1,3
    //   write(*,'(3(f11.5,1x))') ABC(j,:)
    //enddo
    //
    // First thing to do is find the (real, not lattice) point in the affine plane B,C + A that is
    // nearest the origin. Call this T.
    lambda=-scalar_product(A,vector_product(B,C))/(aurostd::modulus(vector_product(B,C))*aurostd::modulus(vector_product(B,C)));
    T=A + lambda*vector_product(B,C);
  
    //print *,lambda
    //write(*,'("T (vec in B,C affine plane): ",3(f10.3,1x))') T
    
    // Now find the four points of the B,C lattice, in the affine plane, that enclose the point T
    for(int i=1;i<=3;i++) {//GH We need these 3 lines to load matrix ABC again with the vectors A,B,C
      ABC(i,1)=A(i);ABC(i,2)=B(i);ABC(i,3)=C(i);//GH
    }//GH
    gusstd::matrix_inverse(ABC,ABCinv,err);

    if(err==TRUE)
      {cerr << _AUROSTD_XLIBS_ERROR_ << "AFLOW_XGUS.CPP: reduce_A_in_ABC  A,B,C vectors in reduce_A_in_ABC are co-planar" << endl;exit(0);}
    i=aurostd::nint(ABCinv*T);
    
    //print *,"Lattice coordinates of origin enclosing T:", i
    // Compute the distance from T to each of the four points and pick the one that is the closest.
    i1(1)=i(1);i1(2)=i(2);i1(3)=i(3);dist(1)=aurostd::modulus(T-ABC*i1);
    i2(1)=i(1);i2(2)=i(2)+1;i2(3)=i(3);dist(2)=aurostd::modulus(T-ABC*i2);
    i3(1)=i(1);i3(2)=i(2);i3(3)=i(3)+1;dist(3)=aurostd::modulus(T-ABC*i3);
    i4(1)=i(1);i4(2)=i(2)+1;i4(3)=i(3)+1;dist(4)=aurostd::modulus(T-ABC*i4);
    idx=0;idx=aurostd::mini(dist);
    //write(*,'("Dists: ",4(f10.5,1x))') dist
    
    //if(.not. equal(,origdist,eps)) then // Only change A if the new one really
    
    if(idx==1) A=A-ABC*i1;
    if(idx==2) A=A-ABC*i2;
    if(idx==3) A=A-ABC*i3;
    if(idx==4) A=A-ABC*i4;
    if(idx==0) {cerr << _AUROSTD_XLIBS_ERROR_ << "FLOW_XGUS.CPP: reduce_A_in_ABC  Case failed in reduce_A_in_ABC" << endl;exit(0);}
    //endif
    //write(*,'("aurostd::modulus(A): ",f7.3,5x," A ",3(f7.3,1x)A)') aurostd::modulus(A), A
    for(int i=1;i<=3;i++) {//GH We need these 3 lines to load matrix ABC again with the vectors A,B,C
      ABC(i,1)=A(i);ABC(i,2)=B(i);ABC(i,3)=C(i);//GH
    }//GH
	  
    gusstd::matrix_inverse(ABC,ABCinv,err);
    //
    ABC=ABCinv*oldABC-aurostd::nint(ABCinv*oldABC);
    for(int i=1;i<=3;i++) {
      for(int j=1;j<=3;j++) {
	if(abs(ABC(i,j))>eps) {
	  cerr << "eps=" << endl << eps << endl;
	  cerr << "ABC(i,j)=" << endl << ABC(i,j) << endl;
	  cerr << "ABCinv=" << endl << ABCinv << endl;
	  cerr << "oldABC=" << endl << oldABC << endl;
	  cerr << "ABCinv*oldABC=" << endl << ABCinv*oldABC << endl;
	  cerr << "ABCinv*oldABC-aurostd::nint(ABCinv*oldABC)=" << endl << ABCinv*oldABC-aurostd::nint(ABCinv*oldABC) << endl;
	  cerr << "Lattice was not preserved  in reduce_A_in_ABC" << endl;
	  exit(0);
	}    
      }
    }
    //read(*,*)
  }
  
  // ***************************************************************************************************
  // This routine takes two vectors (in three-space) and reduces them to form a shortest set (Minkowski
  // reduced). The idea is to subtract B from C so that the new C is as close to the origin as any
  // lattice point along the line that passes through C in the direction of B. The process is repeated
  // then for C subtracted from B, and so on, until the new vector isn't shorter than the other. It's
  // pretty obvious if you do an example by hand. Also see 3.1 of Lecture notes in computer science,
  // ISSN 0302-974, ANTS - VI : algorithmic number theory, 2004, vol. 3076, pp. 338-357 ISBN
  // 3-540-22156-5
  bool gaussian_reduce_two_vectors(xvector<double>& B, xvector<double>& C,double eps) {
    xvector<double> temp(3);
    for(int it=0;;it++) { // dont touch this
      int SwapCnt=0;//GH Counter for the number of times B and C are swapped
      if(it > 100) { cerr << "gaussian_reduce_two_vectors failed to converge" << endl;return FALSE;exit(0);}
      if(aurostd::modulus(B) > aurostd::modulus(C)) {
	temp=B;  // Keep C as the longest vector
	B=C;     // Keep C as the longest vector
	C=temp;  // Keep C as the longest vector
	SwapCnt++; //GH Keep track of the number of swaps
      }
      //    cerr << aurostd::modulus(C) << " " << scalar_product(B,C)/aurostd::modulus(B) << endl;
      C=C-nint(scalar_product(B,C)/aurostd::modulus(B)/aurostd::modulus(C))*B;
      if(aurostd::modulus(C) > aurostd::modulus(B)-eps) {
	if(aurostd::mod(SwapCnt,2)!=0) {temp=B;B=C;C=temp;} // GH CORRECT
	// BUG BUG BUG if(aurostd::mod(it,2)!=0) {temp=B;B=C;C=temp;} // GH
	//GH Make sure the routine doesn't change the order of B and C on output
	// In other words, switch B and C again if odd number of swaps so far
	return TRUE; // basis cannot be further reduced
      }
    }
  }
  
  // ****************************************************************************************
  // Given the matrix a, finds its inverse b
  void matrix_inverse(const xmatrix<double>& a, xmatrix<double>& b, bool& err) {
    double c;
    c=aurostd::det(a);
    if(abs(c)<10e-14) {
      err=TRUE;				      
      b.clear();
    } else {
      err=FALSE;
      b=aurostd::inverse(a);
    }
  }

  void matrix_inverse(const xmatrix<double>& a, xmatrix<double>& b) {
    double c;
    c=aurostd::det(a);
    if(abs(c)<10e-14) b.clear();
    else b=aurostd::inverse(a);
  }
  	
  // ****************************************************************************************
  // This routine finds the determinant of 2x2 or 3x3 matrices
  double determinant_real(const xmatrix<double>& a) {
    return aurostd::det(a);
  }

  double determinant_integer(const xmatrix<double>& a) {
    return aurostd::det(a);
  }
	
  
  // ****************************************************************************************
  // This function takes three vectors and returns a "signed" volume of the parallelpiped
  // that they form
  double volume(const xvector<double>& a1,const xvector<double>& a2,const xvector<double>& a3) {
    return (double) scalar_product(a1, vector_product(a2,a3));
  }
}



// ----------------------------------------------------------------------------

#endif  // _XGUS_IMPLEMENTATIONS_



// **************************************************************************
// *                                                                        *
// *             STEFANO CURTAROLO - Duke University 2003-2008              *
// *                                                                        *
// **************************************************************************

