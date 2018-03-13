// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *                                                                         *
// ***************************************************************************
// Written by Stefano Curtarolo 1994-2011
// Components by Gus Hart 2008

#ifndef _XGUS_H_
#define _XGUS_H_
// __xprototype;

using aurostd::xmatrix;
using aurostd::xvector;
using aurostd::xtensor3;


namespace gusstd {
  // namespace gusstd
  // -------------------------------------------------------------------------
  // --------------------------------------------- gus function and procedures
  bool equal_rank1(const xvector<double>& a, const xvector<double>& b, const double& tolerance);
  bool equal_rank2(const xmatrix<double>& a, const xmatrix<double>& b, const double& tolerance);
  bool equal_rank3(const xtensor3<double>& a, const xtensor3<double>& b, const double& tolerance);
  bool equal_scalar(const double& a, const double& b, const double& tolerance);
  bool equal_rank0(const double& a, const double& b, const double& tolerance);
  bool equal_rank1_rank0(const xvector<double>& a, const double& b, const double& tolerance);
  bool equal_rank2_rank0(const xmatrix<double>& a, const double& b, const double& tolerance);
  bool equal_rank1_real_int(const xvector<double>& a, const xvector<int>& b, const double& tolerance);
  bool equal_rank2_real_int(const xmatrix<double>& a, const xmatrix<int>& b, const double& tolerance);
}

namespace gusstd {
  // namespace gusstd
  // -------------------------------------------------------------------------
  // --------------------------------------------- gus function and procedures
  double orthogonality_defect(const xmatrix<double>& basis);
  double reduce_to_shortest_basis(const xmatrix<double>& IN,xmatrix<double>& OUT,double eps,bool VERBOSE);
  xmatrix<double> reduce_to_shortest_basis(const xmatrix<double>& IN,double eps,bool VERBOSE);
  xmatrix<double> reduce_to_shortest_basis(const xmatrix<double>& IN,bool VERBOSE);
  xmatrix<double> reduce_to_shortest_basis(const xmatrix<double>& IN,double eps);
  xmatrix<double> reduce_to_shortest_basis(const xmatrix<double>& IN);
  void reduce_A_in_ABC(xvector<double>& A, xvector<double>& B, xvector<double>& C,double eps);
  bool gaussian_reduce_two_vectors(xvector<double>& B, xvector<double>& C,double eps);
  void matrix_inverse(const xmatrix<double>& a, xmatrix<double>& b, bool& err);
  void matrix_inverse(const xmatrix<double>& a, xmatrix<double>& b);
  double determinant_real(const xmatrix<double>& a);
  double determinant_integer(const xmatrix<double>& a);
  double volume(const xvector<double>& a1,const xvector<double>& a2,const xvector<double>& a3);
}

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

#endif  // _XGUS_H_

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *                                                                         *
// ***************************************************************************

