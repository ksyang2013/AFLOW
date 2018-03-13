// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *                                                                         *
// ***************************************************************************
// functions written by KESONG YANG
// 2010-2013: kesong.yang@gmail.com

#ifndef _AFLOW_POCCUPATION_FORCEFIELD_CPP_
#define _AFLOW_POCCUPATION_FORCEFIELD_CPP_

#include "aflow.h"        
#include "aflow_contrib_kesong.h"        

//const double unit_kcaltoev=4.336443203200000E-002; // 1(kcal/mol) = 4.33644E-2 eV
const double KCAL_TO_EV =4.336443203200000E-002; // 1(kcal/mol) = 4.33644E-2 eV
//const double KJ_TO_EV =1.0357416650425146E-002 ; // 1 (KJ/mol) = 0.010357 eV
//const double KCAL_TO_KJ = 4.1868; // 1 Kilocalories (Kcal) = 4.1868 Kilojoules (Kj)


// ***************************************************************************
// Bond Class function
// ***************************************************************************

namespace pocc {
  Bond::Bond() {}
  Bond::~Bond() {}
  void Bond::Free() {}
  
  void Bond::Set(xstructure xstr, _atom atomBGN, _atom atomEND) {
    xstr.ReScale(1.0); //safety
    bgn = atomBGN;
    end = atomEND;
    length = AtomDist(bgn, end);
  }
  
  const Bond& Bond::operator=(const Bond& other) {
    if(this != &other) {
      this->bgn = other.bgn;
      this->end = other.end;
      this->length = other.length;
    }
    return *this;
  }
  
  bool Bond::operator==(const Bond &other) const {
    bool FLAG = false;
    if((SameAtom(bgn, other.end) && SameAtom(end, other.bgn)) ||
       (SameAtom(bgn, other.bgn) && SameAtom(end, other.end))) {
        FLAG=true;
    }
    return FLAG;
  }
  
  bool Bond::operator!=(const Bond &other) const {
    return !(*this == other);
  }
  
  ostream& operator<<(ostream& os ,const Bond& bond) {
    cout << bond.bgn << endl;
    cout << bond.end << endl;
    cout << bond.length << endl;
    return os;
  }
}

// ***************************************************************************
// void pocc::RemoveSameBond(vector<Bond>& Bonds_orig, vector<Bond>& Bonds_new)
// ***************************************************************************
namespace pocc {
  void RemoveSameBond(vector<Bond>& Bonds_orig, vector<Bond>& Bonds_new) {
    Bonds_new.push_back(Bonds_orig.at(0));
    for (uint i=1; i<Bonds_orig.size();i++) {
      uint j;
      for (j=0; j<Bonds_new.size();j++) {
	if(Bonds_orig.at(i)==Bonds_new.at(j)) {
	  break;
	}
      }
      if(j==Bonds_new.size()) {
	Bonds_new.push_back(Bonds_orig.at(i));
      }
    }
  }
} // namespace pocc

// ***************************************************************************
// void pocc::SetUFFPara(_atom atomi, _atom atomj, double& R0, double& Kij, double& Xij, double& Dij)
// ***************************************************************************
// Set UFF parameters
namespace pocc {
  void SetUFFPara(_atom atomi, _atom atomj, double& R0, double& Kij, double& Xij, double& Dij) {
    //double bondorder = 1.0; //single bond
    string atomi_name = atomi.name;
    string atomj_name = atomj.name;
    pocc::UFFPara uffparai; uffparai.GetUFFParameters(atomi_name);
    pocc::UFFPara uffparaj; uffparaj.GetUFFParameters(atomj_name);
    double ri = uffparai.r1, rj = uffparaj.r1; //bond length
    double Zi = uffparai.Z1, Zj = uffparaj.Z1; //effective charge
    double chiI = uffparai.Xi, chiJ = uffparaj.Xi; //GMP electronegativity
    double Xi = uffparai.x1, Xj = uffparaj.x1; //nonbond distance
    double Di = uffparai.D1, Dj = uffparaj.D1; //nonbond energy
    //From Equation 3
    //double rbo = -0.1332*(ri+rj)*log(bondorder); //zero for single bond
    //From Equation 4
    double ren = ri*rj*(pow((sqrt(chiI) - sqrt(chiJ)), 2.0)) / (chiI*ri +chiJ*rj);
    //From Equation 2
    //NOTE: See http://towhee.sourceforge.net/forcefields/uff.html
    //There is a typo in the published paper
    //R0 = ri + rj + rbo - ren;
    R0 = ri + rj - ren;
    Kij = 664.12*Zi*Zj/(R0*R0*R0); //From Equation 6
    Xij = sqrt(Xi*Xj); //From Equation 21b
    Dij = sqrt(Di*Dj); //From Equation 22
  }   
} // namespace pocc

  // ***************************************************************************
  // double pocc::CalculateBondEnergy(xstructure xstr, _atom atomi, _atom atomj)
  // ***************************************************************************
namespace pocc {
  double CalculateBondEnergy(xstructure xstr, _atom atomi, _atom atomj) {
    double r0, kij, xij_tmp, dij_tmp;
    pocc::SetUFFPara(atomi, atomj, r0, kij, xij_tmp, dij_tmp);

    xstr.ReScale(1.0); //Safety
    double rij = AtomDist(atomi, atomj); //Notice here rij is r, while r0 is rij in paper
    double delta = rij - r0;
    double delta2 = delta*delta;
    //From Equation 1a
    double energy = 0.5 * kij * delta2;  //unit kcal/mol
    energy = energy*KCAL_TO_EV;
    return energy;
  }
} // namespace pocc

  // ***************************************************************************
  // double pocc::CalculateNonBondEnergy(xstructure xstr, _atom atomi, _atom atomj)
  // ***************************************************************************
namespace pocc {
  double CalculateNonBondEnergy(xstructure xstr, _atom atomi, _atom atomj) {
    //van der Waals (Nonbonded interaction)
    double r0_tmp, kij_tmp, Xij, Dij;
    pocc::SetUFFPara(atomi, atomj, r0_tmp, kij_tmp, Xij, Dij);

    xstr.ReScale(1.0); //Safety
    double x = AtomDist(atomi, atomj); //x distance
    double X6 = pow((Xij/x), 6); 
    double X12 = pow((Xij/x), 12);

    //From Equation 20
    double energy = Dij*(X12-2*X6);
    energy = energy*KCAL_TO_EV;
    return energy;
  }
} // namespace pocc

  // ***************************************************************************
  // double pocc::CalculateUFFEnergy(xstructure xstr)
  // ***************************************************************************
namespace pocc {
  double CalculateUFFEnergy(xstructure xstr) {
    vector<Bond>  Bonds;
    vector<Bond>  NonBonds;
    pocc::AnalyzeBonds(xstr, Bonds, NonBonds);
    double totalenergy=0.0;

    //bond interaction
    double bondenergy=0.0;
    for (uint i=0; i<Bonds.size();i++) {
      Bond bondi = Bonds.at(i);
      _atom atomi = bondi.bgn;
      _atom atomj = bondi.end;
      bondenergy += pocc::CalculateBondEnergy(xstr, atomi, atomj);
    }

    //nonbond interaction, van der Waals
    double nonbondenergy=0.0;
    for (uint i=0; i<NonBonds.size();i++) {
      Bond bondi = NonBonds.at(i);
      _atom atomi = bondi.bgn;
      _atom atomj = bondi.end;
      nonbondenergy += pocc::CalculateNonBondEnergy(xstr, atomi, atomj);
    }

    totalenergy = bondenergy + nonbondenergy;
    return (totalenergy);
  }
} // namespace pocc

// ***************************************************************************
// void pocc::ExtractBonds(const xstructure& xstr, deque<deque<_atom> >& neigh_mat_bonded, deque<deque<_atom> >& neigh_mat_nonbonded)
// ***************************************************************************
namespace pocc {
  void ExtractBonds(const xstructure& xstr, deque<deque<_atom> >& neigh_mat_bonded, deque<deque<_atom> >& neigh_mat_nonbonded) {
    //Extract bonded and nonbonded atoms
    const double radius = 10.00;
    deque<deque<_atom> > neigh_mat;
    xstructure xstr_tmp = xstr;
    xstr_tmp.GetStrNeighData(radius, neigh_mat); // radius 12 angstrom

    deque<_atom>  atom_tmp_bonded;
    deque<_atom>  atom_tmp_nonbonded;
    
    for(uint ia=0; ia<neigh_mat.size();ia++) {
      atom_tmp_bonded.clear(); atom_tmp_nonbonded.clear();
      _atom a = neigh_mat.at(ia).at(0);
      _atom a1 = neigh_mat.at(ia).at(1);
      atom_tmp_bonded.push_back(a); atom_tmp_bonded.push_back(a1); //a, a1 must be bonded
      atom_tmp_nonbonded.push_back(a);

      double dr0 = AtomDist(a,a1);
      for(uint in=2; in<neigh_mat.at(ia).size();in++) {
	_atom an = neigh_mat.at(ia).at(in);
	if(abs(AtomDist(a,an)-dr0)<0.5) {  
	  atom_tmp_bonded.push_back(an);
	}
	else{
	  atom_tmp_nonbonded.push_back(an);
	}
      }
      neigh_mat_bonded.push_back(atom_tmp_bonded);
      neigh_mat_nonbonded.push_back(atom_tmp_nonbonded);
    }
  }
} // namespace pocc


// ***************************************************************************
// void pocc::AnalyzeBonds(const xstructure& xstr, vector<Bond>& Bonds, vector<Bond>& NonBonds)
// ***************************************************************************
namespace pocc {
  void AnalyzeBonds(const xstructure& xstr, vector<Bond>& Bonds, vector<Bond>& NonBonds) {
    deque<deque<_atom> > neigh_mat_bonded;
    deque<deque<_atom> > neigh_mat_nonbonded;
    pocc::ExtractBonds(xstr, neigh_mat_bonded, neigh_mat_nonbonded);
   
    //Bonds 
    vector<vector<Bond> > vv_xstrBonds;
    vector<Bond>  v_xstrBonds;
    for (uint i=0; i<neigh_mat_bonded.size();i++) {
      _atom atomI=neigh_mat_bonded.at(i).at(0);
      v_xstrBonds.clear();
      for (uint j=1; j<neigh_mat_bonded.at(i).size();j++) {
	_atom atomJ = neigh_mat_bonded.at(i).at(j);
	Bond bondIJ; bondIJ.Set(xstr, atomI, atomJ);
	v_xstrBonds.push_back(bondIJ);
      }
      vv_xstrBonds.push_back(v_xstrBonds);
    }
    for (uint i=0;i<vv_xstrBonds.size();i++) {
      for (uint j=0; j<vv_xstrBonds.at(i).size();j++) {
	Bonds.push_back(vv_xstrBonds.at(i).at(j));
      }
    }
    
    //NonBonds
    vector<vector<Bond> > vv_xstrNonBonds;
    vector<Bond>  v_xstrNonBonds;
    for (uint i=0; i<neigh_mat_nonbonded.size();i++) {
      _atom atomI=neigh_mat_nonbonded.at(i).at(0);
      v_xstrNonBonds.clear();
      for (uint j=1; j<neigh_mat_nonbonded.at(i).size();j++) {
	_atom atomJ = neigh_mat_nonbonded.at(i).at(j);
	Bond bondIJ; bondIJ.Set(xstr, atomI, atomJ);
	v_xstrNonBonds.push_back(bondIJ);
      }
      vv_xstrNonBonds.push_back(v_xstrNonBonds);
    }
    for (uint i=0;i<vv_xstrNonBonds.size();i++) {
      for (uint j=0; j<vv_xstrNonBonds.at(i).size();j++) {
	NonBonds.push_back(vv_xstrNonBonds.at(i).at(j));
      }
    }
  }
} // namespace pocc

// ***************************************************************************
// pocc::UFFENERGY(istream& input)
// ***************************************************************************
namespace pocc {
  void UFFENERGY(istream& input) {
    xstructure xstr; 
    xstr.Clear();
    xstr=xstructure(input, IOVASP_POSCAR);
    cout << "total energy: " <<  pocc::CalculateUFFEnergy(xstr) << endl;
  }
} // namespace pocc

#endif // _AFLOW_POCCUPATION_FORCEFIELD_CPP_

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *                                                                         *
// ***************************************************************************

