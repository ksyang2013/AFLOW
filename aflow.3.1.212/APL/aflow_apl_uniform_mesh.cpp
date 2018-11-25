// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2015           *
// *                Aflow PINKU NATH - Duke University 2014-2016             *
// *                                                                         *
// ***************************************************************************
// Written by Pinku Nath
// pn49@duke.edu

#include "aflow_apl.h"

#define _isnegative(a) (a<MIN_EIGEN_TRESHOLD) ? true : false

namespace apl
{
  // ***************************************************************************************
  UniformMesh::UniformMesh(Logger& l):_logger(l)
  {
     _k_index=-10;
     clear();
  }
  // ***************************************************************************************
  UniformMesh::~UniformMesh(){this->clear();}
  // ***************************************************************************************
  void UniformMesh::clear()
  {
    _kpoints.clear();
    _weights.clear();
     _weights.clear();
     _sc_size.clear();
    _rlattice.clear();
     _klattice.clear();
  }
  // ***************************************************************************************
  void UniformMesh::create_uniform_mesh(int na, int nb, int nc, const xstructure& xs)
  {
    _logger<<"Preparing uniform qmesh " <<apl::endl;
    _sc_size.resize(3,0);
    _sc_size[0]=na;
    _sc_size[1]=nb;
    _sc_size[2]=nc;
    xstructure xstr(xs);
    // Setup both lattices
    xstr.FixLattices();
   _rlattice = xstr.lattice;
   _klattice = ReciprocalLattice(_rlattice);

   _kpoints.clear();
    xvector<double> kpoint(3);
    int i=0;
    for(int s = 1; s <=nc; s++)
      for(int r = 1; r <= nb; r++)
        for(int p = 1; p <= na; p++)
          {
            kpoint(1) = ( 2.0 * p - na - 1 ) / ( 2.0 * na );
            kpoint(2) = ( 2.0 * r - nb - 1 ) / ( 2.0 * nb );
            kpoint(3) = ( 2.0 * s - nc - 1 ) / ( 2.0 * nc );
            double k_norm_tmp = std::sqrt(std::pow(kpoint[1], 2)+ std::pow(kpoint[2], 2) + std::pow(kpoint[3], 2));
             if(!_iszero(k_norm_tmp))
             {
                _k_index=i; 
             }
            // Transform to cartesian coordinate
            kpoint = trasp(_klattice) * kpoint;
           _kpoints.push_back(kpoint);
           _weights.push_back(1.0);
          i++;
          }
    //_rlattice.clear();
    //_klattice.clear();
  }
  // ***************************************************************************************
  vector<int> UniformMesh::get_sc_size()
  {
    return _sc_size;
  }
  // ***************************************************************************************
  vector<double> UniformMesh::get_weights()
   {
    return _weights;
   }
  // ***************************************************************************************
   vector<aurostd::xvector<double> > UniformMesh::get_kpoints()
   {
    return _kpoints;
   }
  // ***************************************************************************************
   xmatrix<double> UniformMesh::get_rlattice()
   {
    return _rlattice;
   }
  // ***************************************************************************************
    xmatrix<double> UniformMesh::get_klattice()
    {
      return _klattice;
    }
  // ***************************************************************************************
    int UniformMesh::get_k_index()
    {
        return _k_index;
    }
  // ***************************************************************************************
}//apl namespace end
