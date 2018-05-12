// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *                                                                         *
// ***************************************************************************
// Written by Stefano Curtarolo - 2007-2014
#ifndef _AFLOW_XATOM_CPP
#define _AFLOW_XATOM_CPP
#include "aflow.h"
#include "aflow_pflow.h"
#include "AUROSTD/aurostd_xscalar.h"

#define _calculate_symmetry_default_sgroup_radius_   2.0
#define PLATON_MIN_VOLUME_PER_ATOM   6.0   // for symmetry calculation
#define _PRIM_MULTITHREAD_MIN_ATOMS_THRESHOLD_ 50

#define PLATON_TOLERANCE_ANGLE 1.0
#define PLATON_TOLERANCE_D1 0.25
#define PLATON_TOLERANCE_D2 0.25
#define PLATON_TOLERANCE_D3 0.25

#define _EPS_ 0.02

#define _pocc_precision_ 5

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// _XATOM_PROTOTYPES_

// ***************************************************************************
// ***************************************************************************
// ***************************************************************************
// _ATOM
// look into aflow.h for the definitions

// constructors
_atom::_atom() {
  fpos.clear();
  cpos.clear();
  corigin.clear();
  coord.clear();
  type=0;
  spin=0.0;
  noncoll_spin.clear(); //DX 12/5/17 - magnetic sym (non-collinear)
  spin_is_given=FALSE; //DX 9/21/17 - magnetic sym
  noncoll_spin_is_given=FALSE; //DX 12/5/17 - magnetic sym (non-collinear)
  mass=0.0;
  name="";
  info=0;  // (RHT)
  cleanname="";
  atomic_number=0;
  number=0;
  name_is_given=FALSE;
  sd="";
  print_RHT=true; //CHANGE THIS BACK TO FALSE WHEN DONE DEBUGGING  (RHT)
  print_cartesian=FALSE;
  basis=0;
  ijk.clear();
  isincell=FALSE;
  reference=0.0;
  ireference=0;
  equivalent=-1;
  is_inequivalent=TRUE;
  num_equivalents=0;
  index_iatoms=0;
  verbose=FALSE;
  order_parameter_value=0;
  order_parameter_atom=FALSE;
  partial_occupation_value=1.0;
  partial_occupation_flag=FALSE;
  shell = 0;
}

// destructor
_atom::~_atom() {
  free();
}

void _atom::free() {
}

const _atom& _atom::operator=(const _atom& b) {       // operator=
  if(this != &b) {
    free();
    fpos=b.fpos;
    cpos=b.cpos;
    corigin=b.corigin;
    coord=b.coord;
    type=b.type;
    spin=b.spin;
    spin_is_given=b.spin_is_given; //DX 9/21/17 - magnetic sym
    noncoll_spin=b.noncoll_spin; //DX 12/5/17 - magnetic sym (non-collinear)
    noncoll_spin_is_given=b.noncoll_spin_is_given; //DX 12/5/17 - magnetic sym (non-collinear)
    mass=b.mass;
    name=b.name;
    info=b.info;  // (RHT)
    cleanname=b.cleanname;
    atomic_number=b.atomic_number;
    number=b.number;
    name_is_given=b.name_is_given;
    sd=b.sd;
    print_RHT=b.print_RHT;  // (RHT)
    print_cartesian=b.print_cartesian;
    basis=b.basis;
    ijk=b.ijk;
    isincell=b.isincell;
    reference=b.reference;
    ireference=b.ireference;
    equivalent=b.equivalent;
    is_inequivalent=b.is_inequivalent;
    num_equivalents=b.num_equivalents;
    index_iatoms=b.index_iatoms;
    verbose=b.verbose;
    order_parameter_value=b.order_parameter_value;
    order_parameter_atom=b.order_parameter_atom;
    partial_occupation_value=b.partial_occupation_value;
    partial_occupation_flag=b.partial_occupation_flag;
    shell = b.shell;
  }
  return *this;
}

void _atom::clear(){
  _atom a; (*this)=a;
}

ostream& operator<<(ostream& oss,const _atom& atom) {
  oss.setf(std::ios::fixed,std::ios::floatfield);
  oss.precision(10);
  if(atom.print_RHT==true) {
    //oss << "ATOM COUT-RHT" << endl;
    oss << atom.coord << " " << atom.name;
  }
  else{
    if(atom.verbose==TRUE) {
      oss << " " << endl;
      oss << "ATOM COUT" << endl;
      oss << "type=" << atom.type << endl;
      oss << "spin=" << atom.spin << endl;
      oss << "spin_is_given=" << atom.spin_is_given << endl; //DX 9/21/17 - magnetic sym
      oss << "noncoll_spin=" << atom.noncoll_spin << endl;                   //DX 12/5/17 - magnetic sym (non-collinear)
      oss << "noncoll_spin_is_given=" << atom.noncoll_spin_is_given << endl; //DX 12/5/17 - magnetic sym (non-collinear)
      oss << "mass=" << atom.mass << endl;
      oss << "name=" << atom.name << endl;
      oss << "info=" << atom.info << endl;
      oss << "cleanname=" << atom.cleanname << endl;
      oss << "atomic_number=" << atom.atomic_number << endl;
      oss << "number=" << atom.number << endl;
      oss << "name_is_given=" << atom.name_is_given << endl;
      oss << "sd=" <<  atom.sd << endl;
      oss << "print_cartesian" << atom.print_cartesian << endl;
      oss << "print_RHT" << atom.print_RHT << endl;
      oss << "corigin" << atom.corigin(1) << " " << atom.corigin(2) << " " << atom.corigin(3) << endl;
      oss << "coord" << atom.coord(1) << " " << atom.coord(2) << " " << atom.coord(3) << endl;
      oss << "isincell=" << atom.isincell << endl;
      oss << "reference=" << atom.reference << endl;
      oss << "ireference=" << atom.ireference << endl;
      oss << "equivalent=" << atom.equivalent << endl;
      oss << "is_inequivalent=" << atom.is_inequivalent << endl;
      oss << "num_equivalents=" << atom.num_equivalents << endl;
      oss << "index_iatoms=" << atom.index_iatoms << endl;
      oss << "verbose=" << atom.verbose << endl;
      oss << "order_parameter_value=" << atom.order_parameter_value << endl;
      oss << "order_parameter_atom=" << atom.order_parameter_atom << endl;
      oss << "partial_occupation_value=" << atom.partial_occupation_value << endl;
      oss << "partial_occupation_flag=" << atom.partial_occupation_flag << endl;
      oss << "nearest_neighbour_shell_num= " << atom.shell << endl;
    }
    if(atom.print_cartesian==TRUE) {
      if(atom.verbose) oss << "cartesian" << endl;
      oss << "C " << atom.cpos(1) << " " << atom.cpos(2) << " " << atom.cpos(3);// << endl;
    } else {
      if(atom.verbose) oss << "fractional" << endl;
      oss << "F " << atom.fpos(1) << " " << atom.fpos(2) << " " << atom.fpos(3);//  << endl;
    }
    oss << " T=" << atom.type;
    oss << " B=" << atom.basis;
    oss << " N=" << atom.number;
    //  oss << setw(1);
    oss << " ijk=[" << atom.ijk(1) << "," << atom.ijk(2) << "," << atom.ijk(3) << "]";
    if(atom.verbose) oss << " " << endl;
  }
  return oss;
}

bool isequalRHT(const _atom& a, const _atom& b,double tol) {
  bool out = false;
  //DX and CO - START
  if(abs(a.fpos(1)-b.fpos(1)) < tol && 
     abs(a.fpos(2)-b.fpos(2)) < tol && 
     abs(a.fpos(3)-b.fpos(3)) < tol &&
     //DX and CO - END
     a.name == b.name) {
    out = true;
  }
  return out;
}


void _atom::CleanName(void) {
  // This function cleanup the name from VASP stuff
  string name1,name2;
  name1=name+" ";
  name2=name+" ";
  name1=name1.substr(0,1);
  if(!((name1[0]>=65 && name1[0]<=90)||(name1[0]>=97 && name1[0]<=122))) name1="";
  if((name1[0]>=97 && name1[0]<=122)) name1[0]-=-97+65;
  name2=name2.substr(1,1);
  if(!((name2[0]>=65 && name2[0]<=90)||(name2[0]>=97 && name2[0]<=122))) name2="";
  if((name2[0]>=65 && name2[0]<=90)) name2[0]+=-97+65;
  cleanname=name1+name2;
  atomic_number=0;
  for(uint j=0;j< atom_symbol_vec.size();j++) if(cleanname==atom_symbol_vec.at(j)) atomic_number=j;
}

void _atom::CleanSpin(void) {
  spin=0.0;
  spin_is_given=FALSE; //DX 9/21/17 - magnetic sym
  noncoll_spin.clear();            //DX 12/5/17 - magnetic sym (non-collinear)
  noncoll_spin_is_given=FALSE; //DX 12/5/17 - magnetic sym (non-collinear)
  if(name.find("+")!=string::npos) {spin=atof(name.substr(name.find("+")).c_str()); spin_is_given=TRUE;} //DX 9/21/17 - magnetic sym
  if(name.find("-")!=string::npos) {spin=atof(name.substr(name.find("-")).c_str()); spin_is_given=TRUE;} //DX 9/21/17 - magnetic sym
}

std::vector<string> atom_symbol_vec(NUM_ELEMENTS);        // store starting from ONE
std::vector<string> atom_name_vec(NUM_ELEMENTS);          // store starting from ONE
std::vector<double> atom_mass_vec(NUM_ELEMENTS);          // store starting from ONE
std::vector<double> atom_volume_vec(NUM_ELEMENTS);        // store starting from ONE
std::vector<int> atom_valence_iupac_vec(NUM_ELEMENTS);    // store starting from ONE http://en.wikipedia.org/wiki/Valence_(chemistry)
std::vector<int> atom_valence_std_vec(NUM_ELEMENTS);      // store starting from ONE http://en.wikipedia.org/wiki/Valence_(chemistry)
std::vector<double> atom_miedema_phi_star(NUM_ELEMENTS);  // store starting from ONE Miedema Rule Table 1a Physica 100B (1980) 1-28  
std::vector<double> atom_miedema_nws(NUM_ELEMENTS);       // store starting from ONE Miedema Rule Table 1a Physica 100B (1980) 1-28
std::vector<double> atom_miedema_Vm(NUM_ELEMENTS);        // store starting from ONE Miedema Rule Table 1a Physica 100B (1980) 1-28
std::vector<double> atom_miedema_gamma_s(NUM_ELEMENTS);   // store starting from ONE Miedema Rule Table 1a Physica 100B (1980) 1-28
std::vector<double> atom_miedema_BVm(NUM_ELEMENTS);       // store starting from ONE Miedema Rule Table 1a Physica 100B (1980) 1-28
// for lanthines from J.A. Alonso and N.H. March. Electrons in Metals and Alloys, Academic Press, London (1989) (except La)
std::vector<double> atom_radius_vec(NUM_ELEMENTS);        // store starting from ONE
std::vector<double> atom_radius_covalent_vec(NUM_ELEMENTS);// store starting from ONE //DX and CO - 9/4/17 
std::vector<double> atom_electronegativity_vec(NUM_ELEMENTS);        // store starting from ONE
std::vector<string> atom_crystal_vec(NUM_ELEMENTS);        // store starting from ONE
std::vector<double> xray_scatt_vec(NUM_ELEMENTS);         // store starting from ONE
std::vector<double> pettifor_scale(NUM_ELEMENTS);         // store starting from ONE Chemical Scale Pettifor Solid State Communications 51 31-34 1984

void atoms_initialize(void) {
  for(int i=0;i<NUM_ELEMENTS;i++) {       // clear
    xray_scatt_vec.at(i)=(double)i-1;      // shift+1
    atom_mass_vec.at(i)=AMU2KILOGRAM;      // masses in kilos
    atom_volume_vec.at(i)=NNN;             // atomic volume in A^3 from the FCC vasp table and/or successive calculations
    atom_valence_iupac_vec.at(i)=NNN;            // IUPAC Maximum number of univalent atoms that may combine with an atom of the element under consideration, or with a fragment, or for which an atom of this element can be substituted.
    atom_valence_std_vec.at(i)=NNN;            // stanmdard: number electrons minus closed shell at leff (noble gas)
    atom_miedema_phi_star.at(i)=NNN;                 // Miedema Rule Table 1a Physica 100B (1980) 1-28   (phi^\star in (V))
    atom_miedema_nws.at(i)=NNN;                      // Miedema Rule Table 1a Physica 100B (1980) 1-28   n_{ws}^{1/3} in (d.u.)^1/3
    atom_miedema_Vm.at(i)=NNN;                       // Miedema Rule Table 1a Physica 100B (1980) 1-28   V_m^{2/3} in (cm^2)
    atom_miedema_gamma_s.at(i)=NNN;                  // Miedema Rule Table 1a Physica 100B (1980) 1-28   \gamma_s^0 in (mJ/m^2)
    atom_miedema_BVm.at(i)=NNN;                      // Miedema Rule Table 1a Physica 100B (1980) 1-28   BV_m (kJ/mole)
    atom_radius_vec.at(i)=NNN;             // Saxena (nm)
    atom_radius_covalent_vec.at(i)=NNN;    // Codero (Angstroms) //DX and CO - 9/4/17
    atom_electronegativity_vec.at(i)=NNN;  // Saxena
    atom_crystal_vec.at(i)="nnn";          // Ashcroft-Mermin
  }

  // Xray_scatt_vector All data collected from the NIST online tables
  // http://physics.nist.gov/PhysRefData/FFast/html/form.html
  // All data are ideally for f1 values for Cu-alpha (wavelength=1.5418A, E=8.0416keV).
  // These are for E=7.9026keV (Cu-alpha is wavelength=1.5418A, E=8.0416keV).

  // All data collected from the online tables:
  // http://www-cxro.lbl.gov/optical_constants/pert_form.html
  // All data are f1 values for Cu-alpha (wavelength=1.5418A, E=8.0416keV].

  int i;
  // ROW 0
  i=0; atom_symbol_vec[i]="XX"; atom_name_vec[i]="UNDEFINED";    atom_mass_vec[i]*=0.00;      atom_volume_vec[i]=NNN;      atom_valence_std_vec[i]=0; atom_valence_iupac_vec[i]=0;   atom_miedema_phi_star[i]=NNN;   atom_miedema_nws[i]=NNN;   atom_miedema_Vm[i]=NNN;  atom_miedema_gamma_s[i]= NNN;   atom_miedema_BVm[i]=NNN;   atom_radius_vec[i]=NNN;    atom_radius_covalent_vec[i]=NNN;   atom_electronegativity_vec[i]=NNN;  atom_crystal_vec[i]="nnn"; pettifor_scale[i]=0; xray_scatt_vec[i]=0;

  // ROW 1
  // s-electron systems
  i++; atom_symbol_vec[i]="H";  atom_name_vec[i]="Hydrogen";     atom_mass_vec[i]*=1.0079;    atom_volume_vec[i]= 0.75110; atom_valence_std_vec[i]=1; atom_valence_iupac_vec[i]=1;  atom_miedema_phi_star[i]=5.2;   atom_miedema_nws[i]=1.5;   atom_miedema_Vm[i]=NNN;  atom_miedema_gamma_s[i]= NNN;   atom_miedema_BVm[i]=NNN;   atom_radius_vec[i]=0.046;  atom_radius_covalent_vec[i]=0.31;   atom_electronegativity_vec[i]=2.10;  atom_crystal_vec[i]="hex";  pettifor_scale[i]=0; xray_scatt_vec[i]=1.000;  // H volume wrong *dimer*   MIEDEMA = PAUL VAN DER PUT book
  i++; atom_symbol_vec[i]="He"; atom_name_vec[i]="Helium";       atom_mass_vec[i]*=4.0026;    atom_volume_vec[i]= -1.000;  atom_valence_std_vec[i]=0; atom_valence_iupac_vec[i]=0;  atom_miedema_phi_star[i]=NNN;   atom_miedema_nws[i]=NNN;   atom_miedema_Vm[i]=NNN;  atom_miedema_gamma_s[i]= NNN;   atom_miedema_BVm[i]=NNN;   atom_radius_vec[i]=NNN;    atom_radius_covalent_vec[i]=0.28;   atom_electronegativity_vec[i]=NNN;   atom_crystal_vec[i]="hcp";  pettifor_scale[i]=0; xray_scatt_vec[i]=2.000;  // He

  // ROW2
  // s-electron systems
  i++; atom_symbol_vec[i]="Li"; atom_name_vec[i]="Lithium";      atom_mass_vec[i]*=6.941;     atom_volume_vec[i]=20.24110; atom_valence_std_vec[i]=1; atom_valence_iupac_vec[i]=1;  atom_miedema_phi_star[i]=2.85;  atom_miedema_nws[i]=0.98;  atom_miedema_Vm[i]=5.5;  atom_miedema_gamma_s[i]= 530;   atom_miedema_BVm[i]=1.5;   atom_radius_vec[i]=0.152;  atom_radius_covalent_vec[i]=1.28;   atom_electronegativity_vec[i]=0.98;  atom_crystal_vec[i]="bcc";  pettifor_scale[i]=0.45; xray_scatt_vec[i]=3.00145;  // Li
  i++; atom_symbol_vec[i]="Be"; atom_name_vec[i]="Beryllium";    atom_mass_vec[i]*=9.0122;    atom_volume_vec[i]= 7.83290; atom_valence_std_vec[i]=2; atom_valence_iupac_vec[i]=2;  atom_miedema_phi_star[i]=4.20;  atom_miedema_nws[i]=1.60;  atom_miedema_Vm[i]=2.9;  atom_miedema_gamma_s[i]=1900;   atom_miedema_BVm[i]=4.9;   atom_radius_vec[i]=0.114;  atom_radius_covalent_vec[i]=0.96;   atom_electronegativity_vec[i]=1.57;  atom_crystal_vec[i]="hcp";  pettifor_scale[i]=1.50;
  // p-electron systems
  i++; atom_symbol_vec[i]="B";  atom_name_vec[i]="Boron";        atom_mass_vec[i]*=10.81;     atom_volume_vec[i]= 5.88420; atom_valence_std_vec[i]=3; atom_valence_iupac_vec[i]=3;  atom_miedema_phi_star[i]=4.75;  atom_miedema_nws[i]=1.55;  atom_miedema_Vm[i]=2.8;  atom_miedema_gamma_s[i]= NNN;   atom_miedema_BVm[i]=NNN;   atom_radius_vec[i]=0.097;  atom_radius_covalent_vec[i]=0.84;   atom_electronegativity_vec[i]=2.04;  atom_crystal_vec[i]="tet";  pettifor_scale[i]=2.00;
  i++; atom_symbol_vec[i]="C";  atom_name_vec[i]="Carbon";       atom_mass_vec[i]*=12.011;    atom_volume_vec[i]= 5.59490; atom_valence_std_vec[i]=4; atom_valence_iupac_vec[i]=4;  atom_miedema_phi_star[i]=6.20;  atom_miedema_nws[i]=1.90;  atom_miedema_Vm[i]=1.8;  atom_miedema_gamma_s[i]= NNN;   atom_miedema_BVm[i]=NNN;   atom_radius_vec[i]=0.077;  atom_radius_covalent_vec[i]=0.76;   atom_electronegativity_vec[i]=2.55;  atom_crystal_vec[i]="dia";  pettifor_scale[i]=2.50; xray_scatt_vec[i]=6.019;   // C  // DX and CO - 9/4/17 atom_radius_covalent_vec uses sp3 hybridization (most common)
  i++; atom_symbol_vec[i]="N";  atom_name_vec[i]="Nitrogen";     atom_mass_vec[i]*=14.0067;   atom_volume_vec[i]= 7.59940; atom_valence_std_vec[i]=5; atom_valence_iupac_vec[i]=5;  atom_miedema_phi_star[i]=7.00;  atom_miedema_nws[i]=1.60;  atom_miedema_Vm[i]=2.2;  atom_miedema_gamma_s[i]= NNN;   atom_miedema_BVm[i]=NNN;   atom_radius_vec[i]=0.071;  atom_radius_covalent_vec[i]=0.71;   atom_electronegativity_vec[i]=3.04;  atom_crystal_vec[i]="hex";  pettifor_scale[i]=3.00; //JUNKAI CHANGED VALENCE
  i++; atom_symbol_vec[i]="O";  atom_name_vec[i]="Oxygen";       atom_mass_vec[i]*=15.9994;   atom_volume_vec[i]= 7.78230; atom_valence_std_vec[i]=6; atom_valence_iupac_vec[i]=2;  atom_miedema_phi_star[i]=6.97;  atom_miedema_nws[i]=1.70;  atom_miedema_Vm[i]=2.656;atom_miedema_gamma_s[i]= NNN;   atom_miedema_BVm[i]=NNN;   atom_radius_vec[i]=0.060;  atom_radius_covalent_vec[i]=0.66;   atom_electronegativity_vec[i]=3.44;  atom_crystal_vec[i]="cub";  pettifor_scale[i]=3.50; xray_scatt_vec[i]=8.052;  // O Table 27 of JUNKAI
  i++; atom_symbol_vec[i]="F";  atom_name_vec[i]="Fluorine";     atom_mass_vec[i]*=18.9984;   atom_volume_vec[i]= 9.99090; atom_valence_std_vec[i]=7; atom_valence_iupac_vec[i]=1;  atom_miedema_phi_star[i]=NNN;   atom_miedema_nws[i]=NNN;   atom_miedema_Vm[i]=NNN;  atom_miedema_gamma_s[i]= NNN;   atom_miedema_BVm[i]=NNN;   atom_radius_vec[i]=NNN;    atom_radius_covalent_vec[i]=0.57;   atom_electronegativity_vec[i]=3.98;  atom_crystal_vec[i]="mcl";  pettifor_scale[i]=4.00;
  i++; atom_symbol_vec[i]="Ne";atom_name_vec[i]="Neon";        atom_mass_vec[i]*=20.179;   atom_volume_vec[i]=19.9052; atom_valence_std_vec[i]=0; atom_valence_iupac_vec[i]=0; atom_miedema_phi_star[i]=NNN;  atom_miedema_nws[i]=NNN;  atom_miedema_Vm[i]=NNN; atom_miedema_gamma_s[i]= NNN;  atom_miedema_BVm[i]=NNN;  atom_radius_vec[i]=0.160; atom_radius_covalent_vec[i]=0.58;   atom_electronegativity_vec[i]=NNN; atom_crystal_vec[i]="fcc";  pettifor_scale[i]=0; // volume calculated with fcc-pawpbe

  // ROW3
  // s-electron systems
  i++; atom_symbol_vec[i]="Na";atom_name_vec[i]="Sodium";      atom_mass_vec[i]*=22.9898;  atom_volume_vec[i]=36.9135; atom_valence_std_vec[i]=1; atom_valence_iupac_vec[i]=1; atom_miedema_phi_star[i]=2.70; atom_miedema_nws[i]=0.82; atom_miedema_Vm[i]=8.3; atom_miedema_gamma_s[i]= 260;  atom_miedema_BVm[i]=1.6;  atom_radius_vec[i]=0.186; atom_radius_covalent_vec[i]=1.66;   atom_electronegativity_vec[i]=0.93; atom_crystal_vec[i]="bcc";  pettifor_scale[i]=0.40;
  i++; atom_symbol_vec[i]="Mg";atom_name_vec[i]="Magnesium";   atom_mass_vec[i]*=24.305;   atom_volume_vec[i]=22.8178; atom_valence_std_vec[i]=2; atom_valence_iupac_vec[i]=2; atom_miedema_phi_star[i]=3.45; atom_miedema_nws[i]=1.17; atom_miedema_Vm[i]=5.8; atom_miedema_gamma_s[i]= 790;  atom_miedema_BVm[i]=5.0;  atom_radius_vec[i]=0.160; atom_radius_covalent_vec[i]=1.41;   atom_electronegativity_vec[i]=1.31; atom_crystal_vec[i]="hcp";  pettifor_scale[i]=1.28;
  // p-electron systems
  i++; atom_symbol_vec[i]="Al";atom_name_vec[i]="Aluminium";   atom_mass_vec[i]*=26.9815;  atom_volume_vec[i]=16.4000; atom_valence_std_vec[i]=3; atom_valence_iupac_vec[i]=3; atom_miedema_phi_star[i]=4.20; atom_miedema_nws[i]=1.39; atom_miedema_Vm[i]=4.6; atom_miedema_gamma_s[i]=1200;  atom_miedema_BVm[i]=7.2;  atom_radius_vec[i]=0.143; atom_radius_covalent_vec[i]=1.21;   atom_electronegativity_vec[i]=1.61; atom_crystal_vec[i]="fcc";  pettifor_scale[i]=1.66;
  i++; atom_symbol_vec[i]="Si";atom_name_vec[i]="Silicon";     atom_mass_vec[i]*=28.0855;  atom_volume_vec[i]=14.3536; atom_valence_std_vec[i]=4; atom_valence_iupac_vec[i]=4; atom_miedema_phi_star[i]=4.70; atom_miedema_nws[i]=1.50; atom_miedema_Vm[i]=4.2; atom_miedema_gamma_s[i]=1290;  atom_miedema_BVm[i]=11.9; atom_radius_vec[i]=0.117; atom_radius_covalent_vec[i]=1.11;   atom_electronegativity_vec[i]=1.90; atom_crystal_vec[i]="dia";  pettifor_scale[i]=1.92; xray_scatt_vec[i]=14.43; // Si ???
  i++; atom_symbol_vec[i]="P"; atom_name_vec[i]="Phosphorus";  atom_mass_vec[i]*=30.9738;  atom_volume_vec[i]=14.1995; atom_valence_std_vec[i]=5; atom_valence_iupac_vec[i]=5; atom_miedema_phi_star[i]=5.5;  atom_miedema_nws[i]=1.65; atom_miedema_Vm[i]=NNN; atom_miedema_gamma_s[i]= NNN;  atom_miedema_BVm[i]=NNN;  atom_radius_vec[i]=0.109; atom_radius_covalent_vec[i]=1.07;   atom_electronegativity_vec[i]=2.19; atom_crystal_vec[i]="cub";  pettifor_scale[i]=2.18; xray_scatt_vec[i]=15.3133; // P   MIEDEMA = PAUL VAN DER PUT book
  i++; atom_symbol_vec[i]="S"; atom_name_vec[i]="Sulphur";     atom_mass_vec[i]*=32.06;    atom_volume_vec[i]=15.7301; atom_valence_std_vec[i]=6; atom_valence_iupac_vec[i]=6; atom_miedema_phi_star[i]=5.6;  atom_miedema_nws[i]=1.46; atom_miedema_Vm[i]=4.376;atom_miedema_gamma_s[i]= NNN;  atom_miedema_BVm[i]=NNN;  atom_radius_vec[i]=0.106; atom_radius_covalent_vec[i]=1.05;   atom_electronegativity_vec[i]=2.58; atom_crystal_vec[i]="orc";  pettifor_scale[i]=2.44;  // S Table 27 of JUNKAI
  i++; atom_symbol_vec[i]="Cl";atom_name_vec[i]="Chlorine";    atom_mass_vec[i]*=35.453;   atom_volume_vec[i]=21.2947; atom_valence_std_vec[i]=7; atom_valence_iupac_vec[i]=7; atom_miedema_phi_star[i]=5.32;  atom_miedema_nws[i]=0.34;  atom_miedema_Vm[i]=6.71; atom_miedema_gamma_s[i]= 1013;  atom_miedema_BVm[i]=NNN;  atom_radius_vec[i]=0.107; atom_radius_covalent_vec[i]=1.02;   atom_electronegativity_vec[i]=3.16; atom_crystal_vec[i]="orc";  pettifor_scale[i]=2.70; // interpolation phi_star, nws, Vm, gamma  JUNKAI CHANGED VALENCE    
  i++; atom_symbol_vec[i]="Ar";atom_name_vec[i]="Argon";       atom_mass_vec[i]*=39.948;   atom_volume_vec[i]=22.000; atom_valence_std_vec[i]=0; atom_valence_iupac_vec[i]=2; atom_miedema_phi_star[i]=NNN;  atom_miedema_nws[i]=NNN;  atom_miedema_Vm[i]=NNN; atom_miedema_gamma_s[i]= NNN;  atom_miedema_BVm[i]=NNN;  atom_radius_vec[i]=0.192; atom_radius_covalent_vec[i]=1.06;   atom_electronegativity_vec[i]=NNN;  atom_crystal_vec[i]="fcc";  pettifor_scale[i]=0;  // guessed volume, must double check from results  JUNKAI CHANGED VALENCE    

  // ROW4
  // s-electron systems
  i++; atom_symbol_vec[i]="K"; atom_name_vec[i]="Potassium";   atom_mass_vec[i]*=39.0983;  atom_volume_vec[i]=73.9091; atom_valence_std_vec[i]=1; atom_valence_iupac_vec[i]=1; atom_miedema_phi_star[i]=2.25; atom_miedema_nws[i]=0.65; atom_miedema_Vm[i]=12.8;atom_miedema_gamma_s[i]= 150;  atom_miedema_BVm[i]=1.5;  atom_radius_vec[i]=0.231; atom_radius_covalent_vec[i]=2.03;   atom_electronegativity_vec[i]=0.82; atom_crystal_vec[i]="fcc";  pettifor_scale[i]=0.35;
  i++; atom_symbol_vec[i]="Ca";atom_name_vec[i]="Calcium";     atom_mass_vec[i]*=40.08;    atom_volume_vec[i]=42.1927; atom_valence_std_vec[i]=2; atom_valence_iupac_vec[i]=2; atom_miedema_phi_star[i]=2.55; atom_miedema_nws[i]=0.91; atom_miedema_Vm[i]=8.8; atom_miedema_gamma_s[i]= 490;  atom_miedema_BVm[i]=4.0;  atom_radius_vec[i]=0.197; atom_radius_covalent_vec[i]=1.76;   atom_electronegativity_vec[i]=1.00; atom_crystal_vec[i]="bcc";  pettifor_scale[i]=0.60;
  // d-electron systems: transition metals
  i++; atom_symbol_vec[i]="Sc";atom_name_vec[i]="Scandium";    atom_mass_vec[i]*=44.9559;  atom_volume_vec[i]=24.6739; atom_valence_std_vec[i]=3; atom_valence_iupac_vec[i]=3; atom_miedema_phi_star[i]=3.25; atom_miedema_nws[i]=1.27; atom_miedema_Vm[i]=6.1; atom_miedema_gamma_s[i]=1200;  atom_miedema_BVm[i]=6.6;  atom_radius_vec[i]=0.160; atom_radius_covalent_vec[i]=1.70;   atom_electronegativity_vec[i]=1.36; atom_crystal_vec[i]="hcp";  pettifor_scale[i]=0.74; xray_scatt_vec[i]=21.34; // Sc
  i++; atom_symbol_vec[i]="Ti";atom_name_vec[i]="Titanium";    atom_mass_vec[i]*=47.9;     atom_volume_vec[i]=17.1035; atom_valence_std_vec[i]=4; atom_valence_iupac_vec[i]=4; atom_miedema_phi_star[i]=3.65; atom_miedema_nws[i]=1.47; atom_miedema_Vm[i]=4.8; atom_miedema_gamma_s[i]=2050;  atom_miedema_BVm[i]=11.0; atom_radius_vec[i]=0.147; atom_radius_covalent_vec[i]=1.60;   atom_electronegativity_vec[i]=1.54; atom_crystal_vec[i]="hcp";  pettifor_scale[i]=0.79; xray_scatt_vec[i]=22.24; // Ti
  i++; atom_symbol_vec[i]="V"; atom_name_vec[i]="Vanadium";    atom_mass_vec[i]*=50.9415;  atom_volume_vec[i]=13.2086; atom_valence_std_vec[i]=5; atom_valence_iupac_vec[i]=5; atom_miedema_phi_star[i]=4.25; atom_miedema_nws[i]=1.64; atom_miedema_Vm[i]=4.1; atom_miedema_gamma_s[i]=2600;  atom_miedema_BVm[i]=14.0; atom_radius_vec[i]=0.132; atom_radius_covalent_vec[i]=1.53;   atom_electronegativity_vec[i]=1.63; atom_crystal_vec[i]="bcc";  pettifor_scale[i]=0.84;   // V
  i++; atom_symbol_vec[i]="Cr";atom_name_vec[i]="Chromium";    atom_mass_vec[i]*=51.996;   atom_volume_vec[i]=11.4136; atom_valence_std_vec[i]=6; atom_valence_iupac_vec[i]=6; atom_miedema_phi_star[i]=4.65; atom_miedema_nws[i]=1.74; atom_miedema_Vm[i]=3.7; atom_miedema_gamma_s[i]=2400;  atom_miedema_BVm[i]=14.0; atom_radius_vec[i]=0.125; atom_radius_covalent_vec[i]=1.39;   atom_electronegativity_vec[i]=1.66; atom_crystal_vec[i]="bcc";  pettifor_scale[i]=0.89; xray_scatt_vec[i]=23.84; // Cr
  i++; atom_symbol_vec[i]="Mn";atom_name_vec[i]="Manganese";   atom_mass_vec[i]*=54.93805; atom_volume_vec[i]=10.6487; atom_valence_std_vec[i]=7; atom_valence_iupac_vec[i]=7; atom_miedema_phi_star[i]=4.45; atom_miedema_nws[i]=1.61; atom_miedema_Vm[i]=3.8; atom_miedema_gamma_s[i]=1600;  atom_miedema_BVm[i]=4.4;  atom_radius_vec[i]=0.112; atom_radius_covalent_vec[i]=1.61;   atom_electronegativity_vec[i]=1.55; atom_crystal_vec[i]="cub";  pettifor_scale[i]=0.94; xray_scatt_vec[i]=24.46; // xray_scatt_vec[i]=24.3589; Mn JUNKAI CHANGED VALENCE //DX and CO - 9/4/17 atom_radius_covalent_vec[i] uses high spin configuration (most frequent)   
  i++; atom_symbol_vec[i]="Fe";atom_name_vec[i]="Iron";        atom_mass_vec[i]*=55.847;   atom_volume_vec[i]=10.2315; atom_valence_std_vec[i]=8; atom_valence_iupac_vec[i]=6; atom_miedema_phi_star[i]=4.93; atom_miedema_nws[i]=1.77; atom_miedema_Vm[i]=3.7; atom_miedema_gamma_s[i]=2550;  atom_miedema_BVm[i]=12.0; atom_radius_vec[i]=0.124; atom_radius_covalent_vec[i]=1.52;   atom_electronegativity_vec[i]=1.83; atom_crystal_vec[i]="bcc";  pettifor_scale[i]=0.99; xray_scatt_vec[i]=24.85; // xray_scatt_vec[i]=24.6830; Fe JUNKAI CHANGED VALENCE //DX and CO - 9/4/17 atom_radius_covalent_vec[i] uses high spin configuration (most frequent) 
  i++; atom_symbol_vec[i]="Co";atom_name_vec[i]="Cobalt";      atom_mass_vec[i]*=58.9332;  atom_volume_vec[i]=10.3205; atom_valence_std_vec[i]=9; atom_valence_iupac_vec[i]=5; atom_miedema_phi_star[i]=5.10; atom_miedema_nws[i]=1.75; atom_miedema_Vm[i]=3.5; atom_miedema_gamma_s[i]=2550;  atom_miedema_BVm[i]=13.0; atom_radius_vec[i]=0.125; atom_radius_covalent_vec[i]=1.26;   atom_electronegativity_vec[i]=1.88; atom_crystal_vec[i]="hcp";  pettifor_scale[i]=1.04; xray_scatt_vec[i]=24.59; // Co JUNKAI CHANGED VALENCE //DX and CO - 9/4/17 atom_radius_covalent_vec[i] uses low spin configuration (most frequent)   
  i++; atom_symbol_vec[i]="Ni";atom_name_vec[i]="Nickel";      atom_mass_vec[i]*=58.69;    atom_volume_vec[i]=10.8664; atom_valence_std_vec[i]=10; atom_valence_iupac_vec[i]=4; atom_miedema_phi_star[i]=5.20; atom_miedema_nws[i]=1.75; atom_miedema_Vm[i]=3.5; atom_miedema_gamma_s[i]=2450;  atom_miedema_BVm[i]=12.0; atom_radius_vec[i]=0.125; atom_radius_covalent_vec[i]=1.24;   atom_electronegativity_vec[i]=1.91; atom_crystal_vec[i]="fcc";  pettifor_scale[i]=1.09; xray_scatt_vec[i]=25.02; // Ni
  i++; atom_symbol_vec[i]="Cu";atom_name_vec[i]="Copper";      atom_mass_vec[i]*=63.546;   atom_volume_vec[i]=12.0159; atom_valence_std_vec[i]=11; atom_valence_iupac_vec[i]=4; atom_miedema_phi_star[i]=4.55; atom_miedema_nws[i]=1.47; atom_miedema_Vm[i]=3.7; atom_miedema_gamma_s[i]=1850;  atom_miedema_BVm[i]=9.3;  atom_radius_vec[i]=0.128; atom_radius_covalent_vec[i]=1.32;   atom_electronegativity_vec[i]=1.90; atom_crystal_vec[i]="fcc";  pettifor_scale[i]=1.20; xray_scatt_vec[i]=27.03; // Cu JUNKAI CHANGED VALENCE    
  i++; atom_symbol_vec[i]="Zn";atom_name_vec[i]="Zinc";        atom_mass_vec[i]*=65.38;    atom_volume_vec[i]=15.0827; atom_valence_std_vec[i]=12; atom_valence_iupac_vec[i]=2; atom_miedema_phi_star[i]=4.10; atom_miedema_nws[i]=1.32; atom_miedema_Vm[i]=4.4; atom_miedema_gamma_s[i]=1020;  atom_miedema_BVm[i]=5.5;  atom_radius_vec[i]=0.133; atom_radius_covalent_vec[i]=1.22;   atom_electronegativity_vec[i]=1.65; atom_crystal_vec[i]="hcp";  pettifor_scale[i]=1.44; xray_scatt_vec[i]=28.44; // Zn
  // p-electron systems
  i++; atom_symbol_vec[i]="Ga";atom_name_vec[i]="Gallium";     atom_mass_vec[i]*=69.737;   atom_volume_vec[i]=18.9039; atom_valence_std_vec[i]=3; atom_valence_iupac_vec[i]=3; atom_miedema_phi_star[i]=4.10; atom_miedema_nws[i]=1.31; atom_miedema_Vm[i]=5.2; atom_miedema_gamma_s[i]= 830;  atom_miedema_BVm[i]=6.7;  atom_radius_vec[i]=0.135; atom_radius_covalent_vec[i]=1.22;   atom_electronegativity_vec[i]=1.81; atom_crystal_vec[i]="orc";  pettifor_scale[i]=1.68;
  i++; atom_symbol_vec[i]="Ge";atom_name_vec[i]="Germanium";   atom_mass_vec[i]*=72.59;    atom_volume_vec[i]=19.2948; atom_valence_std_vec[i]=4; atom_valence_iupac_vec[i]=4; atom_miedema_phi_star[i]=4.55; atom_miedema_nws[i]=1.37; atom_miedema_Vm[i]=4.6; atom_miedema_gamma_s[i]=1030;  atom_miedema_BVm[i]=10.5; atom_radius_vec[i]=0.122; atom_radius_covalent_vec[i]=1.20;   atom_electronegativity_vec[i]=2.01; atom_crystal_vec[i]="dia";  pettifor_scale[i]=1.92;
  i++; atom_symbol_vec[i]="As";atom_name_vec[i]="Arsenic";     atom_mass_vec[i]*=74.9216;  atom_volume_vec[i]=19.0677; atom_valence_std_vec[i]=5; atom_valence_iupac_vec[i]=5; atom_miedema_phi_star[i]=4.80; atom_miedema_nws[i]=1.44; atom_miedema_Vm[i]=5.2; atom_miedema_gamma_s[i]=1000;  atom_miedema_BVm[i]=5.1;  atom_radius_vec[i]=0.125; atom_radius_covalent_vec[i]=1.19;   atom_electronegativity_vec[i]=2.18; atom_crystal_vec[i]="rhl";  pettifor_scale[i]=2.16;
  i++; atom_symbol_vec[i]="Se";atom_name_vec[i]="Selenium";    atom_mass_vec[i]*=78.96;    atom_volume_vec[i]=20.3733; atom_valence_std_vec[i]=6; atom_valence_iupac_vec[i]=6; atom_miedema_phi_star[i]=5.17; atom_miedema_nws[i]=1.40; atom_miedema_Vm[i]=5.172;atom_miedema_gamma_s[i]= NNN;  atom_miedema_BVm[i]=NNN;  atom_radius_vec[i]=0.116; atom_radius_covalent_vec[i]=1.20;   atom_electronegativity_vec[i]=2.55; atom_crystal_vec[i]="hex";  pettifor_scale[i]=2.40;  // Se Table 27 of JUNKAI
  i++; atom_symbol_vec[i]="Br";atom_name_vec[i]="Bromine";     atom_mass_vec[i]*=79.904;   atom_volume_vec[i]=26.3292; atom_valence_std_vec[i]=7; atom_valence_iupac_vec[i]=7; atom_miedema_phi_star[i]=5.20; atom_miedema_nws[i]=1.35; atom_miedema_Vm[i]=7.31; atom_miedema_gamma_s[i]= 943;  atom_miedema_BVm[i]=3.4;  atom_radius_vec[i]=0.119; atom_radius_covalent_vec[i]=1.20;   atom_electronegativity_vec[i]=2.96; atom_crystal_vec[i]="orc";  pettifor_scale[i]=2.64; // interpolation phi_star, nws, Vm, gamma, BVm JUNKAI CHANGED VALENCE    
  i++; atom_symbol_vec[i]="Kr";atom_name_vec[i]="Krypton";     atom_mass_vec[i]*=83.8;     atom_volume_vec[i]=-1.0000; atom_valence_std_vec[i]=0; atom_valence_iupac_vec[i]=2; atom_miedema_phi_star[i]=NNN;  atom_miedema_nws[i]=NNN;  atom_miedema_Vm[i]=NNN; atom_miedema_gamma_s[i]= NNN;  atom_miedema_BVm[i]=NNN;  atom_radius_vec[i]=0.197; atom_radius_covalent_vec[i]=1.16;   atom_electronegativity_vec[i]=3; atom_crystal_vec[i]="fcc";  pettifor_scale[i]=0;

  // ROW5
  // s-electron systems
  i++; atom_symbol_vec[i]="Rb";atom_name_vec[i]="Rubidium";    atom_mass_vec[i]*=85.4678;  atom_volume_vec[i]=91.2738; atom_valence_std_vec[i]=1; atom_valence_iupac_vec[i]=1; atom_miedema_phi_star[i]=2.10; atom_miedema_nws[i]=0.60; atom_miedema_Vm[i]=14.6;atom_miedema_gamma_s[i]= 120;  atom_miedema_BVm[i]=1.8;  atom_radius_vec[i]=0.251; atom_radius_covalent_vec[i]=2.20;   atom_electronegativity_vec[i]=0.82; atom_crystal_vec[i]="bcc";  pettifor_scale[i]=0.30;
  i++; atom_symbol_vec[i]="Sr";atom_name_vec[i]="Strontium";   atom_mass_vec[i]*=87.62;    atom_volume_vec[i]=55.4105; atom_valence_std_vec[i]=2; atom_valence_iupac_vec[i]=2; atom_miedema_phi_star[i]=2.40; atom_miedema_nws[i]=0.84; atom_miedema_Vm[i]=10.2;atom_miedema_gamma_s[i]= 430;  atom_miedema_BVm[i]=3.9;  atom_radius_vec[i]=0.215; atom_radius_covalent_vec[i]=1.95;   atom_electronegativity_vec[i]=0.95; atom_crystal_vec[i]="fcc";  pettifor_scale[i]=0.55;
  // d-electron systems: transition metals
  i++; atom_symbol_vec[i]="Y"; atom_name_vec[i]="Yttrium";     atom_mass_vec[i]*=88.9059;  atom_volume_vec[i]=32.4546; atom_valence_std_vec[i]=3; atom_valence_iupac_vec[i]=3; atom_miedema_phi_star[i]=3.20; atom_miedema_nws[i]=1.21; atom_miedema_Vm[i]=7.3; atom_miedema_gamma_s[i]=1100;  atom_miedema_BVm[i]=7.2;  atom_radius_vec[i]=0.181; atom_radius_covalent_vec[i]=1.90;   atom_electronegativity_vec[i]=1.22; atom_crystal_vec[i]="hcp";  pettifor_scale[i]=0.70;
  i++; atom_symbol_vec[i]="Zr";atom_name_vec[i]="Zirconium";   atom_mass_vec[i]*=91.22;    atom_volume_vec[i]=23.2561; atom_valence_std_vec[i]=4; atom_valence_iupac_vec[i]=4; atom_miedema_phi_star[i]=3.40; atom_miedema_nws[i]=1.39; atom_miedema_Vm[i]=5.8; atom_miedema_gamma_s[i]=1950;  atom_miedema_BVm[i]=12.0; atom_radius_vec[i]=0.158; atom_radius_covalent_vec[i]=1.75;   atom_electronegativity_vec[i]=1.33; atom_crystal_vec[i]="hcp";  pettifor_scale[i]=0.76;
  i++; atom_symbol_vec[i]="Nb";atom_name_vec[i]="Niobium";     atom_mass_vec[i]*=92.9064;  atom_volume_vec[i]=18.3132; atom_valence_std_vec[i]=5; atom_valence_iupac_vec[i]=5; atom_miedema_phi_star[i]=4.00; atom_miedema_nws[i]=1.62; atom_miedema_Vm[i]=4.9; atom_miedema_gamma_s[i]=2700;  atom_miedema_BVm[i]=18.0; atom_radius_vec[i]=0.143; atom_radius_covalent_vec[i]=1.64;   atom_electronegativity_vec[i]=1.60; atom_crystal_vec[i]="bcc";  pettifor_scale[i]=0.82;
  i++; atom_symbol_vec[i]="Mo";atom_name_vec[i]="Molybdenum";  atom_mass_vec[i]*=95.94;    atom_volume_vec[i]=15.6175; atom_valence_std_vec[i]=6; atom_valence_iupac_vec[i]=6; atom_miedema_phi_star[i]=4.65; atom_miedema_nws[i]=1.77; atom_miedema_Vm[i]=4.4; atom_miedema_gamma_s[i]=2950;  atom_miedema_BVm[i]=26.0; atom_radius_vec[i]=0.136; atom_radius_covalent_vec[i]=1.54;   atom_electronegativity_vec[i]=2.16; atom_crystal_vec[i]="bcc";  pettifor_scale[i]=0.88;
  i++; atom_symbol_vec[i]="Tc";atom_name_vec[i]="Technetium";  atom_mass_vec[i]*=98.9062;  atom_volume_vec[i]=14.4670; atom_valence_std_vec[i]=7; atom_valence_iupac_vec[i]=7; atom_miedema_phi_star[i]=5.30; atom_miedema_nws[i]=1.81; atom_miedema_Vm[i]=4.2; atom_miedema_gamma_s[i]=3050;  atom_miedema_BVm[i]=26.0; atom_radius_vec[i]=NNN;   atom_radius_covalent_vec[i]=1.47;   atom_electronegativity_vec[i]=1.90; atom_crystal_vec[i]="hcp";  pettifor_scale[i]=0.94; //JUNKAI CHANGED VALENCE    
  i++; atom_symbol_vec[i]="Ru";atom_name_vec[i]="Ruthenium";   atom_mass_vec[i]*=101.07;   atom_volume_vec[i]=13.8390; atom_valence_std_vec[i]=8; atom_valence_iupac_vec[i]=8; atom_miedema_phi_star[i]=5.40; atom_miedema_nws[i]=1.83; atom_miedema_Vm[i]=4.1; atom_miedema_gamma_s[i]=3050;  atom_miedema_BVm[i]=26.0; atom_radius_vec[i]=0.134; atom_radius_covalent_vec[i]=1.46;   atom_electronegativity_vec[i]=2.20; atom_crystal_vec[i]="hcp";  pettifor_scale[i]=1.00; //JUNKAI CHANGED VALENCE    
  i++; atom_symbol_vec[i]="Rh";atom_name_vec[i]="Rhodium";     atom_mass_vec[i]*=102.9055; atom_volume_vec[i]=14.1731; atom_valence_std_vec[i]=9; atom_valence_iupac_vec[i]=6; atom_miedema_phi_star[i]=5.40; atom_miedema_nws[i]=1.76; atom_miedema_Vm[i]=4.1; atom_miedema_gamma_s[i]=2750;  atom_miedema_BVm[i]=23.0; atom_radius_vec[i]=0.134; atom_radius_covalent_vec[i]=1.42;   atom_electronegativity_vec[i]=2.28; atom_crystal_vec[i]="fcc";  pettifor_scale[i]=1.06;
  i++; atom_symbol_vec[i]="Pd";atom_name_vec[i]="Palladium";   atom_mass_vec[i]*=106.4;    atom_volume_vec[i]=15.4596; atom_valence_std_vec[i]=10; atom_valence_iupac_vec[i]=4; atom_miedema_phi_star[i]=5.45; atom_miedema_nws[i]=1.67; atom_miedema_Vm[i]=4.3; atom_miedema_gamma_s[i]=2100;  atom_miedema_BVm[i]=16.0; atom_radius_vec[i]=0.137; atom_radius_covalent_vec[i]=1.39;   atom_electronegativity_vec[i]=2.20; atom_crystal_vec[i]="fcc";  pettifor_scale[i]=1.12;
  i++; atom_symbol_vec[i]="Ag";atom_name_vec[i]="Silver";      atom_mass_vec[i]*=107.8682; atom_volume_vec[i]=18.0678; atom_valence_std_vec[i]=11; atom_valence_iupac_vec[i]=4; atom_miedema_phi_star[i]=4.45; atom_miedema_nws[i]=1.39; atom_miedema_Vm[i]=4.7; atom_miedema_gamma_s[i]=1250;  atom_miedema_BVm[i]=10.0; atom_radius_vec[i]=0.144; atom_radius_covalent_vec[i]=1.45;   atom_electronegativity_vec[i]=1.93; atom_crystal_vec[i]="fcc";  pettifor_scale[i]=1.18; xray_scatt_vec[i]=47.18; // Ag JUNKAI CHANGED VALENCE    
  i++; atom_symbol_vec[i]="Cd";atom_name_vec[i]="Cadmium";     atom_mass_vec[i]*=112.41;   atom_volume_vec[i]=22.0408; atom_valence_std_vec[i]=12; atom_valence_iupac_vec[i]=2; atom_miedema_phi_star[i]=4.05; atom_miedema_nws[i]=1.24; atom_miedema_Vm[i]=5.5; atom_miedema_gamma_s[i]= 780;  atom_miedema_BVm[i]=6.10; atom_radius_vec[i]=0.150; atom_radius_covalent_vec[i]=1.44;   atom_electronegativity_vec[i]=1.69; atom_crystal_vec[i]="hcp";  pettifor_scale[i]=1.36;
  // p-electron systems
  i++; atom_symbol_vec[i]="In";atom_name_vec[i]="Indium";      atom_mass_vec[i]*=114.82;   atom_volume_vec[i]=27.5233; atom_valence_std_vec[i]=3; atom_valence_iupac_vec[i]=3; atom_miedema_phi_star[i]=3.90; atom_miedema_nws[i]=1.17; atom_miedema_Vm[i]=6.3; atom_miedema_gamma_s[i]= 690;  atom_miedema_BVm[i]=6.4;  atom_radius_vec[i]=0.157; atom_radius_covalent_vec[i]=1.42;   atom_electronegativity_vec[i]=1.78; atom_crystal_vec[i]="fct";  pettifor_scale[i]=1.60;
  i++; atom_symbol_vec[i]="Sn";atom_name_vec[i]="Tin";         atom_mass_vec[i]*=118.69;   atom_volume_vec[i]=27.5555; atom_valence_std_vec[i]=4; atom_valence_iupac_vec[i]=4; atom_miedema_phi_star[i]=4.15; atom_miedema_nws[i]=1.24; atom_miedema_Vm[i]=6.4; atom_miedema_gamma_s[i]= 710;  atom_miedema_BVm[i]=8.8;  atom_radius_vec[i]=0.158; atom_radius_covalent_vec[i]=1.39;   atom_electronegativity_vec[i]=1.96; atom_crystal_vec[i]="bct";  pettifor_scale[i]=1.84;
  i++; atom_symbol_vec[i]="Sb";atom_name_vec[i]="Antimony";    atom_mass_vec[i]*=121.75;   atom_volume_vec[i]=27.1823; atom_valence_std_vec[i]=5; atom_valence_iupac_vec[i]=5; atom_miedema_phi_star[i]=4.40; atom_miedema_nws[i]=1.26; atom_miedema_Vm[i]=6.6; atom_miedema_gamma_s[i]= 680;  atom_miedema_BVm[i]=7.0;  atom_radius_vec[i]=0.161; atom_radius_covalent_vec[i]=1.39;   atom_electronegativity_vec[i]=2.05; atom_crystal_vec[i]="rhl";  pettifor_scale[i]=2.08;
  i++; atom_symbol_vec[i]="Te";atom_name_vec[i]="Tellurium";   atom_mass_vec[i]*=127.6;    atom_volume_vec[i]=28.1993; atom_valence_std_vec[i]=6; atom_valence_iupac_vec[i]=6; atom_miedema_phi_star[i]=4.72; atom_miedema_nws[i]=1.31;atom_miedema_Vm[i]=6.439;atom_miedema_gamma_s[i]= NNN;  atom_miedema_BVm[i]=NNN;  atom_radius_vec[i]=0.143; atom_radius_covalent_vec[i]=1.38;   atom_electronegativity_vec[i]=2.10; atom_crystal_vec[i]="hex";  pettifor_scale[i]=2.32;  // Te Table 27 of JUNKAI
  i++; atom_symbol_vec[i]="I"; atom_name_vec[i]="Iodine";     atom_mass_vec[i]*=126.9045; atom_volume_vec[i]=34.9784;  atom_valence_std_vec[i]=7; atom_valence_iupac_vec[i]=7; atom_miedema_phi_star[i]=5.33; atom_miedema_nws[i]=0.17;  atom_miedema_Vm[i]=8.72;atom_miedema_gamma_s[i]= NNN;  atom_miedema_BVm[i]=NNN;atom_radius_vec[i]=0.136; atom_radius_covalent_vec[i]=1.39;   atom_electronegativity_vec[i]=2.66; atom_crystal_vec[i]="orc";  pettifor_scale[i]=2.56; // interpolation phi_star, nws, Vm,
  i++; atom_symbol_vec[i]="Xe";atom_name_vec[i]="Xenon";       atom_mass_vec[i]*=131.3;    atom_volume_vec[i]=-1.0000; atom_valence_std_vec[i]=0; atom_valence_iupac_vec[i]=8; atom_miedema_phi_star[i]=NNN;  atom_miedema_nws[i]=NNN;  atom_miedema_Vm[i]=NNN; atom_miedema_gamma_s[i]= NNN;  atom_miedema_BVm[i]=NNN;  atom_radius_vec[i]=0.218; atom_radius_covalent_vec[i]=1.40;   atom_electronegativity_vec[i]=2.60; atom_crystal_vec[i]="fcc";  pettifor_scale[i]=0; //JUNKAI CHANGED VALENCE    

  // ROW6
  // s-electron systems
  i++; atom_symbol_vec[i]="Cs";atom_name_vec[i]="Cesium";      atom_mass_vec[i]*=132.9054; atom_volume_vec[i]=117.281; atom_valence_std_vec[i]=1; atom_valence_iupac_vec[i]=1; atom_miedema_phi_star[i]=1.95; atom_miedema_nws[i]=0.55; atom_miedema_Vm[i]=16.8;atom_miedema_gamma_s[i]=  95;  atom_miedema_BVm[i]=1.4;  atom_radius_vec[i]=0.265; atom_radius_covalent_vec[i]=2.44;   atom_electronegativity_vec[i]=0.79; atom_crystal_vec[i]="bcc";  pettifor_scale[i]=0.25;
  i++; atom_symbol_vec[i]="Ba";atom_name_vec[i]="Barium";      atom_mass_vec[i]*=137.33;   atom_volume_vec[i]=62.6649; atom_valence_std_vec[i]=2; atom_valence_iupac_vec[i]=2; atom_miedema_phi_star[i]=2.32; atom_miedema_nws[i]=0.81; atom_miedema_Vm[i]=11.3;atom_miedema_gamma_s[i]= 370;  atom_miedema_BVm[i]=3.9;  atom_radius_vec[i]=0.217; atom_radius_covalent_vec[i]=2.15;   atom_electronegativity_vec[i]=0.89; atom_crystal_vec[i]="bcc";  pettifor_scale[i]=0.50;
  // d-electron systems: transition metals
  i++; atom_symbol_vec[i]="La";atom_name_vec[i]="Lanthanium";  atom_mass_vec[i]*=138.9055; atom_volume_vec[i]=36.8495; atom_valence_std_vec[i]=3; atom_valence_iupac_vec[i]=3; atom_miedema_phi_star[i]=3.05;/*3.17*/ atom_miedema_nws[i]=1.09;/*1.18*/ atom_miedema_Vm[i]=8.0;/*7.98*/ atom_miedema_gamma_s[i]= 900;  atom_miedema_BVm[i]=5.5;  atom_radius_vec[i]=0.187; atom_radius_covalent_vec[i]=2.07;   atom_electronegativity_vec[i]=1.10; atom_crystal_vec[i]="hex";  pettifor_scale[i]=0.7480;
  // lantanidies
  i++; atom_symbol_vec[i]="Ce";atom_name_vec[i]="Cerium";      atom_mass_vec[i]*=140.12;   atom_volume_vec[i]=26.4729; atom_valence_std_vec[i]=4; atom_valence_iupac_vec[i]=4; atom_miedema_phi_star[i]=3.18;  atom_miedema_nws[i]=1.19;  atom_miedema_Vm[i]=7.76; atom_miedema_gamma_s[i]= NNN;  atom_miedema_BVm[i]=NNN;  atom_radius_vec[i]=0.182; atom_radius_covalent_vec[i]=2.04;   atom_electronegativity_vec[i]=1.12; atom_crystal_vec[i]="fcc";  pettifor_scale[i]=0.7460;  // pettifor linear interpolation // miedema from Alonso-March.
  i++; atom_symbol_vec[i]="Pr";atom_name_vec[i]="Praseodymium";atom_mass_vec[i]*=140.9077; atom_volume_vec[i]=36.4987; atom_valence_std_vec[i]=5; atom_valence_iupac_vec[i]=4; atom_miedema_phi_star[i]=3.19;  atom_miedema_nws[i]=1.20;  atom_miedema_Vm[i]=7.56; atom_miedema_gamma_s[i]= NNN;  atom_miedema_BVm[i]=NNN;  atom_radius_vec[i]=0.183; atom_radius_covalent_vec[i]=2.03;   atom_electronegativity_vec[i]=1.13; atom_crystal_vec[i]="hex";  pettifor_scale[i]=0.7440;  // pettifor linear interpolation
  i++; atom_symbol_vec[i]="Nd";atom_name_vec[i]="Neodymium";   atom_mass_vec[i]*=144.24;   atom_volume_vec[i]=29.6719; atom_valence_std_vec[i]=6; atom_valence_iupac_vec[i]=4; atom_miedema_phi_star[i]=3.19;  atom_miedema_nws[i]=1.20;  atom_miedema_Vm[i]=7.51; atom_miedema_gamma_s[i]= NNN;  atom_miedema_BVm[i]=NNN;  atom_radius_vec[i]=0.182; atom_radius_covalent_vec[i]=2.01;   atom_electronegativity_vec[i]=1.14; atom_crystal_vec[i]="hex";  pettifor_scale[i]=0.7420;  // pettifor linear interpolation JUNKAI CHANGED VALENCE    
  i++; atom_symbol_vec[i]="Pm";atom_name_vec[i]="Promethium";  atom_mass_vec[i]*=146.92;   atom_volume_vec[i]=34.6133; atom_valence_std_vec[i]=7; atom_valence_iupac_vec[i]=3; atom_miedema_phi_star[i]=3.19;  atom_miedema_nws[i]=1.21;  atom_miedema_Vm[i]=7.43; atom_miedema_gamma_s[i]= NNN;  atom_miedema_BVm[i]=NNN;  atom_radius_vec[i]=NNN;   atom_radius_covalent_vec[i]=1.99;   atom_electronegativity_vec[i]=1.13; atom_crystal_vec[i]="hex";  pettifor_scale[i]=0.7400;  // pettifor linear interpolation
  i++; atom_symbol_vec[i]="Sm";atom_name_vec[i]="Samarium";    atom_mass_vec[i]*=150.4;    atom_volume_vec[i]=33.9484; atom_valence_std_vec[i]=8; atom_valence_iupac_vec[i]=3; atom_miedema_phi_star[i]=3.20;  atom_miedema_nws[i]=1.21;  atom_miedema_Vm[i]=7.37; atom_miedema_gamma_s[i]= NNN;  atom_miedema_BVm[i]=NNN;  atom_radius_vec[i]=0.181; atom_radius_covalent_vec[i]=1.98;   atom_electronegativity_vec[i]=1.17; atom_crystal_vec[i]="rhl";  pettifor_scale[i]=0.7380;  // pettifor linear interpolation
  i++; atom_symbol_vec[i]="Eu";atom_name_vec[i]="Europium";    atom_mass_vec[i]*=151.96;   atom_volume_vec[i]=43.1719; atom_valence_std_vec[i]=9; atom_valence_iupac_vec[i]=3; atom_miedema_phi_star[i]=3.20;  atom_miedema_nws[i]=1.21;  atom_miedema_Vm[i]=7.36; atom_miedema_gamma_s[i]= NNN;  atom_miedema_BVm[i]=NNN;  atom_radius_vec[i]=0.204; atom_radius_covalent_vec[i]=1.98;   atom_electronegativity_vec[i]=1.20; atom_crystal_vec[i]="bcc";  pettifor_scale[i]=0.7360;  // pettifor linear interpolation
  i++; atom_symbol_vec[i]="Gd";atom_name_vec[i]="Gadolinium";  atom_mass_vec[i]*=157.25;   atom_volume_vec[i]=32.5777; atom_valence_std_vec[i]=10; atom_valence_iupac_vec[i]=3; atom_miedema_phi_star[i]=3.20;  atom_miedema_nws[i]=1.21;  atom_miedema_Vm[i]=7.34; atom_miedema_gamma_s[i]= NNN;  atom_miedema_BVm[i]=NNN;  atom_radius_vec[i]=0.180; atom_radius_covalent_vec[i]=1.96;   atom_electronegativity_vec[i]=1.20; atom_crystal_vec[i]="hcp";  pettifor_scale[i]=0.7340;  // pettifor linear interpolation
  i++; atom_symbol_vec[i]="Tb";atom_name_vec[i]="Terbium";     atom_mass_vec[i]*=158.9254; atom_volume_vec[i]=32.0200; atom_valence_std_vec[i]=11; atom_valence_iupac_vec[i]=4; atom_miedema_phi_star[i]=3.21;  atom_miedema_nws[i]=1.22;  atom_miedema_Vm[i]=7.20; atom_miedema_gamma_s[i]= NNN;  atom_miedema_BVm[i]=NNN;  atom_radius_vec[i]=0.177; atom_radius_covalent_vec[i]=1.94;   atom_electronegativity_vec[i]=1.10; atom_crystal_vec[i]="hcp";  pettifor_scale[i]=0.7320;  // pettifor linear interpolation
  i++; atom_symbol_vec[i]="Dy";atom_name_vec[i]="Dysprosium";  atom_mass_vec[i]*=162.5;    atom_volume_vec[i]=31.5096; atom_valence_std_vec[i]=12; atom_valence_iupac_vec[i]=4; atom_miedema_phi_star[i]=3.21;  atom_miedema_nws[i]=1.22;  atom_miedema_Vm[i]=7.12; atom_miedema_gamma_s[i]= NNN;  atom_miedema_BVm[i]=NNN;  atom_radius_vec[i]=0.177; atom_radius_covalent_vec[i]=1.92;   atom_electronegativity_vec[i]=1.22; atom_crystal_vec[i]="hcp";  pettifor_scale[i]=0.7300;  // pettifor linear interpolation JUNKAI CHANGED VALENCE    
  i++; atom_symbol_vec[i]="Ho";atom_name_vec[i]="Holmium";     atom_mass_vec[i]*=164.9304; atom_volume_vec[i]=31.0155; atom_valence_std_vec[i]=13; atom_valence_iupac_vec[i]=3; atom_miedema_phi_star[i]=3.22;  atom_miedema_nws[i]=1.22;  atom_miedema_Vm[i]=7.06; atom_miedema_gamma_s[i]= NNN;  atom_miedema_BVm[i]=NNN;  atom_radius_vec[i]=0.176; atom_radius_covalent_vec[i]=1.92;   atom_electronegativity_vec[i]=1.23; atom_crystal_vec[i]="hcp";  pettifor_scale[i]=0.7280;  // pettifor linear interpolation
  i++; atom_symbol_vec[i]="Er";atom_name_vec[i]="Erbium";      atom_mass_vec[i]*=167.26;   atom_volume_vec[i]=30.5431; atom_valence_std_vec[i]=14; atom_valence_iupac_vec[i]=3; atom_miedema_phi_star[i]=3.22;  atom_miedema_nws[i]=1.23;  atom_miedema_Vm[i]=6.98; atom_miedema_gamma_s[i]= NNN;  atom_miedema_BVm[i]=NNN;  atom_radius_vec[i]=0.175; atom_radius_covalent_vec[i]=1.89;   atom_electronegativity_vec[i]=1.24; atom_crystal_vec[i]="hcp";  pettifor_scale[i]=0.7260;  // pettifor linear interpolation
  i++; atom_symbol_vec[i]="Tm";atom_name_vec[i]="Thulium";     atom_mass_vec[i]*=168.9342; atom_volume_vec[i]=30.0016; atom_valence_std_vec[i]=15; atom_valence_iupac_vec[i]=4; atom_miedema_phi_star[i]=3.22;  atom_miedema_nws[i]=1.23;  atom_miedema_Vm[i]=6.90; atom_miedema_gamma_s[i]= NNN;  atom_miedema_BVm[i]=NNN;  atom_radius_vec[i]=0.174; atom_radius_covalent_vec[i]=1.90;   atom_electronegativity_vec[i]=1.25; atom_crystal_vec[i]="hcp";  pettifor_scale[i]=0.7240;  // pettifor linear interpolation JUNKAI CHANGED VALENCE    
  i++; atom_symbol_vec[i]="Yb";atom_name_vec[i]="Ytterbium";   atom_mass_vec[i]*=173.04;   atom_volume_vec[i]=39.4395; atom_valence_std_vec[i]=16; atom_valence_iupac_vec[i]=3; atom_miedema_phi_star[i]=3.22;  atom_miedema_nws[i]=1.23;  atom_miedema_Vm[i]=6.86; atom_miedema_gamma_s[i]= NNN;  atom_miedema_BVm[i]=NNN;  atom_radius_vec[i]=0.193; atom_radius_covalent_vec[i]=1.87;   atom_electronegativity_vec[i]=1.10; atom_crystal_vec[i]="fcc";  pettifor_scale[i]=0.7220;  // pettifor linear interpolation
  i++; atom_symbol_vec[i]="Lu";atom_name_vec[i]="Lutetium";    atom_mass_vec[i]*=174.967;  atom_volume_vec[i]=29.3515; atom_valence_std_vec[i]=17; atom_valence_iupac_vec[i]=3; atom_miedema_phi_star[i]=3.22;  atom_miedema_nws[i]=1.24;  atom_miedema_Vm[i]=6.81; atom_miedema_gamma_s[i]= NNN;  atom_miedema_BVm[i]=NNN;  atom_radius_vec[i]=0.173; atom_radius_covalent_vec[i]=1.87;   atom_electronegativity_vec[i]=1.27; atom_crystal_vec[i]="hcp";  pettifor_scale[i]=0.7200;
  // d-electron systems: transition metalsnnn";
  i++; atom_symbol_vec[i]="Hf";atom_name_vec[i]="Hafnium";     atom_mass_vec[i]*=178.49;   atom_volume_vec[i]=22.0408; atom_valence_std_vec[i]=4; atom_valence_iupac_vec[i]=4; atom_miedema_phi_star[i]=3.55; atom_miedema_nws[i]=1.43; atom_miedema_Vm[i]=5.6; atom_miedema_gamma_s[i]=2200;  atom_miedema_BVm[i]=15.0; atom_radius_vec[i]=0.159; atom_radius_covalent_vec[i]=1.75;   atom_electronegativity_vec[i]=1.30; atom_crystal_vec[i]="hcp";  pettifor_scale[i]=0.775;
  i++; atom_symbol_vec[i]="Ta";atom_name_vec[i]="Tantalum";    atom_mass_vec[i]*=180.9479; atom_volume_vec[i]=18.1100; atom_valence_std_vec[i]=5; atom_valence_iupac_vec[i]=5; atom_miedema_phi_star[i]=4.05; atom_miedema_nws[i]=1.63; atom_miedema_Vm[i]=4.9; atom_miedema_gamma_s[i]=3050;  atom_miedema_BVm[i]=22.0; atom_radius_vec[i]=0.147; atom_radius_covalent_vec[i]=1.70;   atom_electronegativity_vec[i]=1.50; atom_crystal_vec[i]="bcc";  pettifor_scale[i]=0.83;
  i++; atom_symbol_vec[i]="W"; atom_name_vec[i]="Tungsten";    atom_mass_vec[i]*=183.85;   atom_volume_vec[i]=15.9387; atom_valence_std_vec[i]=6; atom_valence_iupac_vec[i]=6; atom_miedema_phi_star[i]=4.80; atom_miedema_nws[i]=1.81; atom_miedema_Vm[i]=4.5; atom_miedema_gamma_s[i]=3300;  atom_miedema_BVm[i]=31.0; atom_radius_vec[i]=0.137; atom_radius_covalent_vec[i]=1.62;   atom_electronegativity_vec[i]=2.36; atom_crystal_vec[i]="bcc";  pettifor_scale[i]=0.885;
  i++; atom_symbol_vec[i]="Re";atom_name_vec[i]="Rhenium";     atom_mass_vec[i]*=186.2;    atom_volume_vec[i]=14.8941; atom_valence_std_vec[i]=7; atom_valence_iupac_vec[i]=7; atom_miedema_phi_star[i]=5.40; atom_miedema_nws[i]=1.86; atom_miedema_Vm[i]=4.3; atom_miedema_gamma_s[i]=3650;  atom_miedema_BVm[i]=33.0; atom_radius_vec[i]=0.138; atom_radius_covalent_vec[i]=1.51;   atom_electronegativity_vec[i]=1.90; atom_crystal_vec[i]="hcp";  pettifor_scale[i]=0.94;
  i++; atom_symbol_vec[i]="Os";atom_name_vec[i]="Osmium";      atom_mass_vec[i]*=190.2;    atom_volume_vec[i]=14.2403; atom_valence_std_vec[i]=8; atom_valence_iupac_vec[i]=8; atom_miedema_phi_star[i]=5.40; atom_miedema_nws[i]=1.85; atom_miedema_Vm[i]=4.2; atom_miedema_gamma_s[i]=3500;  atom_miedema_BVm[i]=35.0; atom_radius_vec[i]=0.135; atom_radius_covalent_vec[i]=1.44;   atom_electronegativity_vec[i]=2.20; atom_crystal_vec[i]="hcp";  pettifor_scale[i]=0.995; //JUNKAI CHANGED VALENCE    
  i++; atom_symbol_vec[i]="Ir";atom_name_vec[i]="Iridium";     atom_mass_vec[i]*=192.22;   atom_volume_vec[i]=14.5561; atom_valence_std_vec[i]=9; atom_valence_iupac_vec[i]=8; atom_miedema_phi_star[i]=5.55; atom_miedema_nws[i]=1.83; atom_miedema_Vm[i]=4.2; atom_miedema_gamma_s[i]=3100;  atom_miedema_BVm[i]=25.0; atom_radius_vec[i]=0.135; atom_radius_covalent_vec[i]=1.41;   atom_electronegativity_vec[i]=2.20; atom_crystal_vec[i]="fcc";  pettifor_scale[i]=1.05; //JUNKAI CHANGED VALENCE    
  i++; atom_symbol_vec[i]="Pt";atom_name_vec[i]="Platinum";    atom_mass_vec[i]*=195.09;   atom_volume_vec[i]=15.7298; atom_valence_std_vec[i]=10; atom_valence_iupac_vec[i]=6; atom_miedema_phi_star[i]=5.65; atom_miedema_nws[i]=1.78; atom_miedema_Vm[i]=4.4; atom_miedema_gamma_s[i]=2550;  atom_miedema_BVm[i]=18.0; atom_radius_vec[i]=0.138; atom_radius_covalent_vec[i]=1.36;   atom_electronegativity_vec[i]=2.28; atom_crystal_vec[i]="fcc";  pettifor_scale[i]=1.105;
  i++; atom_symbol_vec[i]="Au";atom_name_vec[i]="Gold";        atom_mass_vec[i]*=196.9665; atom_volume_vec[i]=18.1904; atom_valence_std_vec[i]=11; atom_valence_iupac_vec[i]=5; atom_miedema_phi_star[i]=5.15; atom_miedema_nws[i]=1.57; atom_miedema_Vm[i]=4.7; atom_miedema_gamma_s[i]=1550;  atom_miedema_BVm[i]=18.0; atom_radius_vec[i]=0.144; atom_radius_covalent_vec[i]=1.36;   atom_electronegativity_vec[i]=2.54; atom_crystal_vec[i]="fcc";  pettifor_scale[i]=1.16; xray_scatt_vec[i]=74.99; // Au
  i++; atom_symbol_vec[i]="Hg";atom_name_vec[i]="Mercury";     atom_mass_vec[i]*=200.59;   atom_volume_vec[i]=29.7156; atom_valence_std_vec[i]=12; atom_valence_iupac_vec[i]=4; atom_miedema_phi_star[i]=4.20; atom_miedema_nws[i]=1.24; atom_miedema_Vm[i]=5.8; atom_miedema_gamma_s[i]= 610;  atom_miedema_BVm[i]=4.0;  atom_radius_vec[i]=0.150; atom_radius_covalent_vec[i]=1.32;   atom_electronegativity_vec[i]=2.00; atom_crystal_vec[i]="rhl";  pettifor_scale[i]=1.32;
  // p-electron systems
  i++; atom_symbol_vec[i]="Tl";atom_name_vec[i]="Thallium";    atom_mass_vec[i]*=204.37;   atom_volume_vec[i]=31.0721; atom_valence_std_vec[i]=3; atom_valence_iupac_vec[i]=3; atom_miedema_phi_star[i]=3.90; atom_miedema_nws[i]=1.12; atom_miedema_Vm[i]=6.6; atom_miedema_gamma_s[i]= 610;  atom_miedema_BVm[i]=6.2;  atom_radius_vec[i]=0.171; atom_radius_covalent_vec[i]=1.45;   atom_electronegativity_vec[i]=1.62; atom_crystal_vec[i]="hcp";  pettifor_scale[i]=1.56;  // electronegativity  2.04=>1.62
  i++; atom_symbol_vec[i]="Pb";atom_name_vec[i]="Lead";        atom_mass_vec[i]*=207.2;    atom_volume_vec[i]=31.6649; atom_valence_std_vec[i]=4; atom_valence_iupac_vec[i]=4; atom_miedema_phi_star[i]=4.10; atom_miedema_nws[i]=1.15; atom_miedema_Vm[i]=6.9; atom_miedema_gamma_s[i]= 610;  atom_miedema_BVm[i]=7.9;  atom_radius_vec[i]=0.175; atom_radius_covalent_vec[i]=1.46;   atom_electronegativity_vec[i]=2.33; atom_crystal_vec[i]="fcc";  pettifor_scale[i]=1.80;
  i++; atom_symbol_vec[i]="Bi";atom_name_vec[i]="Bismuth";     atom_mass_vec[i]*=208.9804; atom_volume_vec[i]=31.5691; atom_valence_std_vec[i]=5; atom_valence_iupac_vec[i]=5; atom_miedema_phi_star[i]=4.15; atom_miedema_nws[i]=1.16; atom_miedema_Vm[i]=7.2; atom_miedema_gamma_s[i]= 550;  atom_miedema_BVm[i]=6.7;  atom_radius_vec[i]=0.182; atom_radius_covalent_vec[i]=1.48;   atom_electronegativity_vec[i]=2.02; atom_crystal_vec[i]="rhl";  pettifor_scale[i]=2.04;
  i++; atom_symbol_vec[i]="Po";atom_name_vec[i]="Polonium";    atom_mass_vec[i]*=209.98;   atom_volume_vec[i]=NNN;     atom_valence_std_vec[i]=6; atom_valence_iupac_vec[i]=6; atom_miedema_phi_star[i]=NNN;  atom_miedema_nws[i]=NNN;  atom_miedema_Vm[i]=NNN; atom_miedema_gamma_s[i]= NNN;  atom_miedema_BVm[i]=NNN;  atom_radius_vec[i]=0.140; atom_radius_covalent_vec[i]=1.40;   atom_electronegativity_vec[i]=2.00; atom_crystal_vec[i]="sc";  pettifor_scale[i]=2.28;
  i++; atom_symbol_vec[i]="At";atom_name_vec[i]="Astatine";    atom_mass_vec[i]*=210;      atom_volume_vec[i]=NNN;     atom_valence_std_vec[i]=7; atom_valence_iupac_vec[i]=7; atom_miedema_phi_star[i]=NNN;  atom_miedema_nws[i]=NNN;  atom_miedema_Vm[i]=NNN; atom_miedema_gamma_s[i]= NNN;  atom_miedema_BVm[i]=NNN;  atom_radius_vec[i]=NNN;   atom_radius_covalent_vec[i]=1.50;   atom_electronegativity_vec[i]=2.20; atom_crystal_vec[i]="nnn";  pettifor_scale[i]=2.52;
  i++; atom_symbol_vec[i]="Rn";atom_name_vec[i]="Radon";       atom_mass_vec[i]*=222;      atom_volume_vec[i]=NNN;     atom_valence_std_vec[i]=0; atom_valence_iupac_vec[i]=6; atom_miedema_phi_star[i]=NNN;  atom_miedema_nws[i]=NNN;  atom_miedema_Vm[i]=NNN; atom_miedema_gamma_s[i]= NNN;  atom_miedema_BVm[i]=NNN;  atom_radius_vec[i]=NNN;   atom_radius_covalent_vec[i]=1.50;   atom_electronegativity_vec[i]=2.2; atom_crystal_vec[i]="fcc";  pettifor_scale[i]=0;

  // ROW7
  // s-electron systems
  i++; atom_symbol_vec[i]="Fr";atom_name_vec[i]="Francium";    atom_mass_vec[i]*=223.02;   atom_volume_vec[i]=NNN;     atom_valence_std_vec[i]=1; atom_valence_iupac_vec[i]=1; atom_miedema_phi_star[i]=NNN;  atom_miedema_nws[i]=NNN;  atom_miedema_Vm[i]=NNN; atom_miedema_gamma_s[i]= NNN;  atom_miedema_BVm[i]=NNN;  atom_radius_vec[i]=NNN;   atom_radius_covalent_vec[i]=2.60;   atom_electronegativity_vec[i]=0.70; atom_crystal_vec[i]="bcc";  pettifor_scale[i]=0;
  i++; atom_symbol_vec[i]="Ra";atom_name_vec[i]="Radium";      atom_mass_vec[i]*=226.0254; atom_volume_vec[i]=-1.0000; atom_valence_std_vec[i]=2; atom_valence_iupac_vec[i]=2; atom_miedema_phi_star[i]=NNN;  atom_miedema_nws[i]=NNN;  atom_miedema_Vm[i]=NNN; atom_miedema_gamma_s[i]= NNN;  atom_miedema_BVm[i]=NNN;  atom_radius_vec[i]=NNN;   atom_radius_covalent_vec[i]=2.21;   atom_electronegativity_vec[i]=0.89; atom_crystal_vec[i]="bct";  pettifor_scale[i]=0;
  // d-electron systems: transition metals
  i++; atom_symbol_vec[i]="Ac";atom_name_vec[i]="Actinium";    atom_mass_vec[i]*=227.03;   atom_volume_vec[i]=45.2437; atom_valence_std_vec[i]=3; atom_valence_iupac_vec[i]=3; atom_miedema_phi_star[i]=NNN;  atom_miedema_nws[i]=NNN;  atom_miedema_Vm[i]=NNN; atom_miedema_gamma_s[i]= NNN;  atom_miedema_BVm[i]=NNN;  atom_radius_vec[i]=NNN;   atom_radius_covalent_vec[i]=2.15;   atom_electronegativity_vec[i]=1.10; atom_crystal_vec[i]="fcc";  pettifor_scale[i]=0;
  // actinidies
  i++; atom_symbol_vec[i]="Th";atom_name_vec[i]="Thorium";     atom_mass_vec[i]*=232.0381; atom_volume_vec[i]=31.9586; atom_valence_std_vec[i]=4; atom_valence_iupac_vec[i]=4; atom_miedema_phi_star[i]=3.30; atom_miedema_nws[i]=1.28; atom_miedema_Vm[i]=7.3; atom_miedema_gamma_s[i]= NNN;  atom_miedema_BVm[i]=NNN;  atom_radius_vec[i]=0.180; atom_radius_covalent_vec[i]=2.06;   atom_electronegativity_vec[i]=1.30; atom_crystal_vec[i]="fcc";  pettifor_scale[i]=0; xray_scatt_vec[i]=86.64; // Th
  i++; atom_symbol_vec[i]="Pa";atom_name_vec[i]="Protoactinium";atom_mass_vec[i]*=231.04;  atom_volume_vec[i]=NNN;     atom_valence_std_vec[i]=5; atom_valence_iupac_vec[i]=5; atom_miedema_phi_star[i]=NNN;  atom_miedema_nws[i]=NNN;  atom_miedema_Vm[i]=NNN; atom_miedema_gamma_s[i]= NNN;  atom_miedema_BVm[i]=NNN;  atom_radius_vec[i]=NNN;   atom_radius_covalent_vec[i]=2.00;   atom_electronegativity_vec[i]=1.50; atom_crystal_vec[i]="bct";  pettifor_scale[i]=0;
  i++; atom_symbol_vec[i]="U"; atom_name_vec[i]="Uranium";     atom_mass_vec[i]*=238.03;   atom_volume_vec[i]=NNN;     atom_valence_std_vec[i]=6; atom_valence_iupac_vec[i]=6; atom_miedema_phi_star[i]=NNN;  atom_miedema_nws[i]=NNN;  atom_miedema_Vm[i]=NNN; atom_miedema_gamma_s[i]= NNN;  atom_miedema_BVm[i]=NNN;  atom_radius_vec[i]=0.138; atom_radius_covalent_vec[i]=1.96;   atom_electronegativity_vec[i]=1.38; atom_crystal_vec[i]="orc";  pettifor_scale[i]=0;
  i++; atom_symbol_vec[i]="Np";atom_name_vec[i]="Neptunium";   atom_mass_vec[i]*=237.05;   atom_volume_vec[i]=NNN;     atom_valence_std_vec[i]=7; atom_valence_iupac_vec[i]=7; atom_miedema_phi_star[i]=NNN;  atom_miedema_nws[i]=NNN;  atom_miedema_Vm[i]=NNN; atom_miedema_gamma_s[i]= NNN;  atom_miedema_BVm[i]=NNN;  atom_radius_vec[i]=NNN;   atom_radius_covalent_vec[i]=1.90;   atom_electronegativity_vec[i]=NNN; atom_crystal_vec[i]="nnn";  pettifor_scale[i]=0;
  i++; atom_symbol_vec[i]="Pu"; atom_name_vec[i]="Plutonium";  atom_mass_vec[i]*=244.06;   atom_volume_vec[i]=NNN;     atom_valence_std_vec[i]=8; atom_valence_iupac_vec[i]=7; atom_miedema_phi_star[i]=NNN;  atom_miedema_nws[i]=NNN;  atom_miedema_Vm[i]=NNN; atom_miedema_gamma_s[i]= NNN;  atom_miedema_BVm[i]=NNN;  atom_radius_vec[i]=NNN; atom_radius_covalent_vec[i]=1.87;   atom_electronegativity_vec[i]=NNN; atom_crystal_vec[i]="nnn";  pettifor_scale[i]=0;
  i++; atom_symbol_vec[i]="Am";atom_name_vec[i]="Americium";   atom_mass_vec[i]*=243.06;   atom_volume_vec[i]=NNN;     atom_valence_std_vec[i]=9; atom_valence_iupac_vec[i]=7; atom_miedema_phi_star[i]=NNN;  atom_miedema_nws[i]=NNN;  atom_miedema_Vm[i]=NNN; atom_miedema_gamma_s[i]= NNN;  atom_miedema_BVm[i]=NNN;  atom_radius_vec[i]=NNN;   atom_radius_covalent_vec[i]=1.80;   atom_electronegativity_vec[i]=NNN; atom_crystal_vec[i]="nnn";  pettifor_scale[i]=0;
  i++; atom_symbol_vec[i]="Cm"; atom_name_vec[i]="Curium";     atom_mass_vec[i]*=247.07;   atom_volume_vec[i]=NNN;     atom_valence_std_vec[i]=10; atom_valence_iupac_vec[i]=8; atom_miedema_phi_star[i]=NNN;  atom_miedema_nws[i]=NNN;  atom_miedema_Vm[i]=NNN; atom_miedema_gamma_s[i]= NNN;  atom_miedema_BVm[i]=NNN;  atom_radius_vec[i]=NNN; atom_radius_covalent_vec[i]=1.69;   atom_electronegativity_vec[i]=NNN; atom_crystal_vec[i]="nnn";  pettifor_scale[i]=0;

  i++; atom_symbol_vec[i]="Bk"; atom_name_vec[i]="Berkelium";  atom_mass_vec[i]*=247.07;   atom_volume_vec[i]=NNN;     atom_valence_std_vec[i]=11; atom_valence_iupac_vec[i]=4; atom_miedema_phi_star[i]=NNN;  atom_miedema_nws[i]=NNN;  atom_miedema_Vm[i]=NNN; atom_miedema_gamma_s[i]= NNN;  atom_miedema_BVm[i]=NNN;  atom_radius_vec[i]=NNN; atom_radius_covalent_vec[i]=NNN;   atom_electronegativity_vec[i]=NNN; atom_crystal_vec[i]="nnn";  pettifor_scale[i]=0;
  i++; atom_symbol_vec[i]="Cf"; atom_name_vec[i]="Californium";atom_mass_vec[i]*=251.08;   atom_volume_vec[i]=NNN;     atom_valence_std_vec[i]=12; atom_valence_iupac_vec[i]=4; atom_miedema_phi_star[i]=NNN;  atom_miedema_nws[i]=NNN;  atom_miedema_Vm[i]=NNN; atom_miedema_gamma_s[i]= NNN;  atom_miedema_BVm[i]=NNN;  atom_radius_vec[i]=NNN; atom_radius_covalent_vec[i]=NNN;   atom_electronegativity_vec[i]=NNN; atom_crystal_vec[i]="nnn";  pettifor_scale[i]=0;
  i++; atom_symbol_vec[i]="Es"; atom_name_vec[i]="Einsteinium";atom_mass_vec[i]*=252.08;   atom_volume_vec[i]=NNN;     atom_valence_std_vec[i]=13; atom_valence_iupac_vec[i]=4; atom_miedema_phi_star[i]=NNN;  atom_miedema_nws[i]=NNN;  atom_miedema_Vm[i]=NNN; atom_miedema_gamma_s[i]= NNN;  atom_miedema_BVm[i]=NNN;  atom_radius_vec[i]=NNN; atom_radius_covalent_vec[i]=NNN;   atom_electronegativity_vec[i]=NNN; atom_crystal_vec[i]="nnn";  pettifor_scale[i]=0;
  i++; atom_symbol_vec[i]="Fm"; atom_name_vec[i]="Fermium";    atom_mass_vec[i]*=257.1;    atom_volume_vec[i]=NNN;     atom_valence_std_vec[i]=14; atom_valence_iupac_vec[i]=3; atom_miedema_phi_star[i]=NNN;  atom_miedema_nws[i]=NNN;  atom_miedema_Vm[i]=NNN; atom_miedema_gamma_s[i]= NNN;  atom_miedema_BVm[i]=NNN;  atom_radius_vec[i]=NNN; atom_radius_covalent_vec[i]=NNN;   atom_electronegativity_vec[i]=NNN; atom_crystal_vec[i]="nnn";  pettifor_scale[i]=0;
  i++; atom_symbol_vec[i]="Md"; atom_name_vec[i]="Mendelevium";atom_mass_vec[i]*=258.1;    atom_volume_vec[i]=NNN;     atom_valence_std_vec[i]=15; atom_valence_iupac_vec[i]=3; atom_miedema_phi_star[i]=NNN;  atom_miedema_nws[i]=NNN;  atom_miedema_Vm[i]=NNN; atom_miedema_gamma_s[i]= NNN;  atom_miedema_BVm[i]=NNN;  atom_radius_vec[i]=NNN; atom_radius_covalent_vec[i]=NNN;   atom_electronegativity_vec[i]=NNN; atom_crystal_vec[i]="nnn";  pettifor_scale[i]=0;
  i++; atom_symbol_vec[i]="No"; atom_name_vec[i]="Nobelium";   atom_mass_vec[i]*=259.1;    atom_volume_vec[i]=NNN;     atom_valence_std_vec[i]=16; atom_valence_iupac_vec[i]=3; atom_miedema_phi_star[i]=NNN;  atom_miedema_nws[i]=NNN;  atom_miedema_Vm[i]=NNN; atom_miedema_gamma_s[i]= NNN;  atom_miedema_BVm[i]=NNN;  atom_radius_vec[i]=NNN; atom_radius_covalent_vec[i]=NNN;   atom_electronegativity_vec[i]=NNN; atom_crystal_vec[i]="nnn";  pettifor_scale[i]=0;
  i++; atom_symbol_vec[i]="Lr"; atom_name_vec[i]="Lawrencium"; atom_mass_vec[i]*=262.11;   atom_volume_vec[i]=NNN;     atom_valence_std_vec[i]=17; atom_valence_iupac_vec[i]=3; atom_miedema_phi_star[i]=NNN;  atom_miedema_nws[i]=NNN;  atom_miedema_Vm[i]=NNN; atom_miedema_gamma_s[i]= NNN;  atom_miedema_BVm[i]=NNN;  atom_radius_vec[i]=NNN; atom_radius_covalent_vec[i]=NNN;   atom_electronegativity_vec[i]=NNN; atom_crystal_vec[i]="nnn";  pettifor_scale[i]=0;

  /*
    int valence_wahyu[]={9999,1,0,1,2,3,4,-3,-2,-1,0,1,2,3,4,-3,-2,-1,0,1,2,3,4,3,3,2,3,3,2,1,2,3,4,-3,-2,-1,0,1,2,3,4,2,4,7,4,3,2,1,2,3,4,-3,-2,-1,0,1,2,3,3,3,3,3,3,2,3,3,3,3,3,3,3,3,4,5,3,4,4,4,2,1,2,3,4,3,0,0,0,1,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3};
    for(uint j=1;j<93;j++) {
    int test=atom_valence_std_vec[j];
    if(j<=9 && test>=5) {test=test-8;} // first row
    if(j<=17 && test>=5) {test=test-8;} // second row
    //   if(test!=valence_wahyu[j]) 
    cerr << atom_symbol_vec[j] << " " << atom_valence_std_vec[j] << " " <<  test  << " " << valence_wahyu[j];
    if(test!=valence_wahyu[j]) cerr << "   ****   ";
    cerr << endl;
    }
    exit(0);
  */


  //   for(int i=1;i<90;i++) {
  //     cerr << i << " " <<  atom_symbol_vec[i] << " " << pettifor_scale[i] << endl;
  //   }
  // cerr << i << " " << atom_symbol_vec[i] << endl; exit(0);

  // aconvasp stuff
  // All data collected from the NIST online tables:
  // http://physics.nist.gov/PhysRefData/FFast/html/form.html
  // All data are ideally for f1 values for Cu-alpha (wavelength=1.5418A, E=8.0416keV).
  // These are for E=7.9026keV (Cu-alpha is wavelength=1.5418A, E=8.0416keV).
  xray_scatt_vec[1+2]=3.00145E+00; // Li
  xray_scatt_vec[1+14]=1.53133E+01; // P
  xray_scatt_vec[1+24]=2.43589E+01; // Mn
  xray_scatt_vec[1+25]=2.46830E+01; // Fe

  // All data collected from the online tables:
  // http://www-cxro.lbl.gov/optical_constants/pert_form.html
  // All data are f1 values for Cu-alpha (wavelength=1.5418A, E=8.0416keV).
  xray_scatt_vec[1+0]=1.000; // H
  xray_scatt_vec[1+1]=2.000; // He
  xray_scatt_vec[1+2]=3.001; // Li
  xray_scatt_vec[1+6]=6.019; // C
  xray_scatt_vec[1+7]=8.052; // O
  xray_scatt_vec[1+13]=14.43; // P
  xray_scatt_vec[1+14]=15.30; // P
  xray_scatt_vec[1+20]=21.34; // Sc
  xray_scatt_vec[1+21]=22.24; // Ti
  xray_scatt_vec[1+23]=23.84; // Cr
  xray_scatt_vec[1+24]=24.46; // Mn
  xray_scatt_vec[1+25]=24.85; // Fe
  xray_scatt_vec[1+26]=24.59; // Co
  xray_scatt_vec[1+27]=25.02; // Ni
  xray_scatt_vec[1+28]=27.03; // Cu
  xray_scatt_vec[1+29]=28.44; // Zn
  xray_scatt_vec[1+46]=47.18; // Ag
  xray_scatt_vec[1+78]=74.99; // Au
  xray_scatt_vec[1+89]=86.64; // Th

  // Atomic masses
  // All indices are the atomic number shifted back by one.
  // All masses are in kilograms

  // not useful anymore all the masses are declared
  //  atom_mass_vec=vector<double> (NUM_ELEMENTS,0.0);
  //  for(i=0;i<NUM_ELEMENTS;i++) {
  //    atom_mass_vec[i]=(double)(2*i)*AMU2KILOGRAM;
  //  }
  //  atom_mass_vec[1+0]=1.0079*AMU2KILOGRAM; //H
  //  atom_mass_vec[1+7]=15.9994*AMU2KILOGRAM; //O
  //  atom_mass_vec[1+24]=54.93805*AMU2KILOGRAM; //Mn

  // finish and copy
  

}

// **************************************************************************
// Function GetAtomNumber
// **************************************************************************
uint GetAtomNumber(const string& symbol) {
  for(uint iat=0;iat<NUM_ELEMENTS;iat++)
    if(symbol==atom_symbol_vec.at(iat) || symbol==atom_name_vec.at(iat))
      return (uint) iat;
  cerr << "GetAtomNumber symbol not found symbol=" << symbol << "  atom_name_vec.size()=" << atom_name_vec.size() << endl;
  return 0; // no symbol found
}

// **************************************************************************
// Function GetAtomName
// **************************************************************************
std::string GetAtomName(const string& symbol) {
  for(int iat=0;iat<NUM_ELEMENTS;iat++)
    if(symbol==atom_symbol_vec.at(iat) || symbol==atom_name_vec.at(iat))
      return atom_name_vec.at(iat);
  return symbol;
}
std::string GetAtomName(const uint& atnum) {
  if(atnum>=atom_name_vec.size() || atnum<=0) {
    cerr << "GetAtomName out of boundary  atnum=" << atnum << "  atom_name_vec.size()=" << atom_name_vec.size() << endl;
    return "not found";
  }
  return atom_name_vec.at(atnum);
}

// **************************************************************************
// Function GetAtomSymbol
// **************************************************************************
std::string GetAtomSymbol(const string& symbol) {
  for(int iat=0;iat<NUM_ELEMENTS;iat++)
    if(symbol==atom_symbol_vec.at(iat) || symbol==atom_symbol_vec.at(iat))
      return atom_symbol_vec.at(iat);
  return symbol;
}
std::string GetAtomSymbol(const uint& atnum) {
  if(atnum>=atom_symbol_vec.size() || atnum<=0) {
    cerr << "GetAtomSymbol out of boundary  atnum=" << atnum << "  atom_symbol_vec.size()=" << atom_symbol_vec.size() << endl;
    return "not found";
  }
  return atom_symbol_vec.at(atnum);
}


// **************************************************************************
// Function GetAtomMass
// **************************************************************************
double GetAtomMass(const string& symbol) {
  for(int iat=0;iat<NUM_ELEMENTS;iat++)
    if(symbol==atom_symbol_vec.at(iat) || symbol==atom_name_vec.at(iat))
      return atom_mass_vec.at(iat);
  return (double) NNN;
}
double GetAtomMass(const uint& atnum) {
  if(atnum>=atom_mass_vec.size() || atnum<=0) {
    cerr << "GetAtomMass out of boundary  atnum=" << atnum << "  atom_mass_vec.size()=" << atom_mass_vec.size() << endl;
    return (double) NNN;
  }
  return atom_mass_vec.at(atnum);
}

// **************************************************************************
// Function GetAtomComptonCrossSection
// **************************************************************************
double GetAtomComptonCrossSection(const string& symbol) {
  //  cout << "symbol=" << symbol << endl;
  //  cout << "GetAtomNumber(symbol)=" << GetAtomNumber(symbol) << endl;
  return (double) GetAtomComptonCrossSection(GetAtomNumber(symbol));
}

double GetAtomComptonCrossSection(const uint& atnum) { //sigma_c at 662 KeV, compton cross section in barn (1 barn = 1e-28 m^2) [ref: Ortiz, Comp Mat Sci 44, 1042 (2009)]
  double sigma_c[]={0.29,0.57,0.86,1.15,1.43,1.72,2.00,2.29,2.58,2.86,3.15,3.43,3.72,4.00,4.29,4.57,4.86,5.14,5.43,5.71,5.99,6.28,6.56,6.84,7.13,7.41,7.69,7.98,8.26,8.54,8.82,9.11,9.39,9.67,9.95,10.24,10.52,10.79,11.08,11.36,11.64,11.92,12.20,12.48,12.76,13.04,13.32,13.60,13.88,14.16,14.44,14.72,15.00,15.27,15.55,15.83,16.11,16.38,16.66,16.94,17.22,17.50,17.77,18.05,18.33,18.60,18.88,19.16,19.44,19.71,19.98,20.26,20.54,20.81,21.08,21.36,21.64,21.91,22.19,22.46,22.73,23.01,23.28,23.55,23.82,24.10,24.38,24.64,24.92,-1.000};
  vector<double> vsigma_c;vsigma_c.push_back(NNN);
  for(uint i=0;sigma_c[i]>0;i++) vsigma_c.push_back(sigma_c[i]);
  if(atnum>=vsigma_c.size() || atnum<=0) {
    cerr << "GetAtomComptonCrossSection out of boundary  atnum=" << atnum << "  vsigma_c.size()=" << vsigma_c.size() << endl;
    return (double) NNN;
  }
  return vsigma_c.at(atnum);
}

// **************************************************************************
// Function GetAtomPhotoelectricCrossSection
// **************************************************************************
double GetAtomPhotoelectricCrossSection(const string& symbol) {
  //  cout << "symbol=" << symbol << endl;
  //  cout << "GetAtomNumber(symbol)=" << GetAtomNumber(symbol) << endl;
  return (double) GetAtomPhotoelectricCrossSection(GetAtomNumber(symbol));
}

double GetAtomPhotoelectricCrossSection(const uint& atnum) { //sigma_pe photoelectric cross section in barn (1 barn = 1e-28 m^2) %ref: Ortiz, Comp Mat Sci 44, 1042 (2009)
  double sigma_pe[]={8.79e-09,2.70e-07,2.62e-06,1.33e-05,4.59e-05,1.22e-04,2.75e-04,5.51e-04,9.85e-04,1.69e-03,2.62e-03,3.97e-03,6.10e-03,8.46e-03,1.18e-02,1.69e-02,2.19e-02,2.87e-02,3.68e-02,4.91e-02,5.94e-02,7.74e-02,9.47e-02,1.16e-01,1.40e-01,1.59e-01,1.94e-01,2.21e-01,2.65e-01,3.10e-01,3.56e-01,4.14e-01,4.68e-01,5.49e-01,6.24e-01,7.09e-01,7.98e-01,9.16e-01,1.02,1.10,1.27,1.49,1.57,1.74,1.92,2.15,2.36,2.60,2.87,3.06,3.30,3.66,4.06,4.31,4.68,5.03,5.47,5.86,6.41,6.95,7.43,7.96,8.56,9.02,9.78,1.05e+01,1.11e+01,1.19e+01,1.27e+01,1.37e+01,1.44e+01,1.48e+01,1.57e+01,1.65e+01,1.78e+01,1.94e+01,2.06e+01,2.11e+01,2.23e+01,2.36e+01,2.56e+01,2.65e+01,2.80e+01,3.00e+01,3.15e+01,3.33e+01,3.47e+01,3.65e+01,3.82e+01,-1.000};
  vector<double> vsigma_pe;vsigma_pe.push_back(NNN);
  for(uint i=0;sigma_pe[i]>0;i++) vsigma_pe.push_back(sigma_pe[i]);
  if(atnum>=vsigma_pe.size() || atnum<=0) {
    cerr << "GetAtomPhotoelectricCrossSection out of boundary  atnum=" << atnum << "  vsigma_pe.size()=" << vsigma_pe.size() << endl;
    return (double) NNN;
  }
  return vsigma_pe.at(atnum);
}

// **************************************************************************
// Function GetAtomVolume
// **************************************************************************
double GetAtomVolume(const string& symbol) {
  for(int iat=0;iat<NUM_ELEMENTS;iat++)
    if(symbol==atom_symbol_vec.at(iat) || symbol==atom_name_vec.at(iat))
      return atom_volume_vec.at(iat);
  return (double) NNN;
}
double GetAtomVolume(const uint& atnum) {
  if(atnum>=atom_volume_vec.size() || atnum<=0) {
    cerr << "GetAtomVolume out of boundary  atnum=" << atnum << "  atom_volume_vec.size()=" << atom_volume_vec.size() << endl;
    return (double) NNN;
  }
  return atom_volume_vec.at(atnum);
}

// **************************************************************************
// Function GetAtomValenceIupac
// **************************************************************************
int GetAtomValenceIupac(const string& symbol) {
  for(int iat=0;iat<NUM_ELEMENTS;iat++)
    if(symbol==atom_symbol_vec.at(iat) || symbol==atom_name_vec.at(iat))
      return atom_valence_iupac_vec.at(iat);
  return (int) NNN;
}

int GetAtomValenceIupac(const uint& atnum) {
  if(atnum>=atom_valence_iupac_vec.size() || atnum<=0) {
    cerr << "GetAtomValenceIupac out of boundary  atnum=" << atnum << "  atom_valence_iupac_vec.size()=" << atom_valence_iupac_vec.size() << endl;
    return (int) NNN;
  }
  return atom_valence_iupac_vec.at(atnum);
}

// **************************************************************************
// Function GetAtomValenceStd
// **************************************************************************
int GetAtomValenceStd(const string& symbol) {
  for(int iat=0;iat<NUM_ELEMENTS;iat++)
    if(symbol==atom_symbol_vec.at(iat) || symbol==atom_name_vec.at(iat))
      return atom_valence_std_vec.at(iat);
  return (int) NNN;
}

int GetAtomValenceStd(const uint& atnum) {
  if(atnum>=atom_valence_std_vec.size() || atnum<=0) {
    cerr << "GetAtomValenceStd out of boundary  atnum=" << atnum << "  atom_valence_std_vec.size()=" << atom_valence_std_vec.size() << endl;
    return (int) NNN;
  }
  return atom_valence_std_vec.at(atnum);
}

// **************************************************************************
// Function GetAtomRadius
// **************************************************************************
double GetAtomRadius(const string& symbol) {
  for(int iat=0;iat<NUM_ELEMENTS;iat++)
    if(symbol==atom_symbol_vec.at(iat) || symbol==atom_name_vec.at(iat))
      return atom_radius_vec.at(iat);
  return (double) NNN;
}
double GetAtomRadius(const uint& atnum) {
  if(atnum>=atom_radius_vec.size() || atnum<=0) {
    cerr << "GetAtomRadius out of boundary  atnum=" << atnum << "  atom_radius_vec.size()=" << atom_radius_vec.size() << endl;
    return (double) NNN;
  }
  return atom_radius_vec.at(atnum);
}

//DX and CO - 9/4/17 - START 
// **************************************************************************
// Function GetAtomRadiusCovalent
// **************************************************************************
double GetAtomRadiusCovalent(const string& symbol) {
  for(int iat=0;iat<NUM_ELEMENTS;iat++)
    if(symbol==atom_symbol_vec.at(iat) || symbol==atom_name_vec.at(iat))
      return atom_radius_covalent_vec.at(iat);
  return (double) NNN;
}
double GetAtomRadiusCovalent(const uint& atnum) {
  if(atnum>=atom_radius_covalent_vec.size() || atnum<=0) {
    cerr << "GetAtomRadiusCovalent out of boundary  atnum=" << atnum << "  atom_radius_covalent_vec.size()=" << atom_radius_covalent_vec.size() << endl;
    return (double) NNN;
  }
  return atom_radius_covalent_vec.at(atnum);
}
//DX and CO - 9/4/17 - END 

// **************************************************************************
// Function GetAtomElectronegativity
// **************************************************************************
double GetAtomElectronegativity(const string& symbol) {
  for(int iat=0;iat<NUM_ELEMENTS;iat++)
    if(symbol==atom_symbol_vec.at(iat) || symbol==atom_name_vec.at(iat))
      return atom_electronegativity_vec.at(iat);
  return (double) NNN;
}
double GetAtomElectronegativity(const uint& atnum) {
  if(atnum>=atom_electronegativity_vec.size() || atnum<=0) {
    cerr << "GetAtomElectronegativity out of boundary  atnum=" << atnum << "  atom_electronegativity_vec.size()=" << atom_electronegativity_vec.size() << endl;
    return (double) NNN;
  }
  return atom_electronegativity_vec.at(atnum);
}

// **************************************************************************
// Function GetAtomCrystal
// **************************************************************************
string GetAtomCrystal(const string& symbol) {
  for(int iat=0;iat<NUM_ELEMENTS;iat++)
    if(symbol==atom_symbol_vec.at(iat) || symbol==atom_name_vec.at(iat))
      return atom_crystal_vec.at(iat);
  return (string) "nnn";
}
string GetAtomCrystal(const uint& atnum) {
  if(atnum>=atom_crystal_vec.size() || atnum<=0) {
    cerr << "GetAtomCrystal out of boundary  atnum=" << atnum << "  atom_crystal_vec.size()=" << atom_crystal_vec.size() << endl;
    return (string) "nnn";
  }
  return atom_crystal_vec.at(atnum);
}

// **************************************************************************
// Function GetAtomPettiforScale
// **************************************************************************
double GetAtomPettiforScale(const string& symbol) {
  for(int iat=0;iat<NUM_ELEMENTS;iat++)
    if(symbol==atom_symbol_vec.at(iat) || symbol==atom_name_vec.at(iat))
      return pettifor_scale.at(iat);
  return 0.0;
}

double GetAtomPettiforScale(const uint& atnum) {
  if(atnum>=atom_crystal_vec.size() || atnum<=0) {
    cerr << "GetAtomPettiforScale out of boundary  atnum=" << atnum << "  atom_crystal_vec.size()=" << atom_crystal_vec.size() << endl;
    return 0.0;
  }
  return pettifor_scale.at(atnum);
}

bool GetAtomPettiforScale(const vector<string>& vsymbol,vector<double>& vvalue) {
  vvalue.clear(); // delete
  for(uint i=0;i<vsymbol.size();i++)
    vvalue.push_back(GetAtomPettiforScale(vsymbol.at(i)));
  return TRUE;
}

bool GetAtomPettiforScale(const vector<uint>& vatnum,vector<double>& vvalue) {
  vvalue.clear(); // delete
  for(uint i=0;i<vatnum.size();i++)
    vvalue.push_back(GetAtomPettiforScale(vatnum.at(i)));
  return TRUE;
}

bool GetAtomPettiforScale(const vector<string>& vsymbol,xvector<double>& vvalue) {
  if(vvalue.rows!=(int) vsymbol.size()) return FALSE; // nothing to be ordered
  if(vvalue.lrows!=1) return FALSE; // start from 1
  for(uint i=0;i<vsymbol.size();i++)
    vvalue[i+1]=GetAtomPettiforScale(vsymbol.at(i));
  return TRUE;
}

bool GetAtomPettiforScale(const vector<uint>& vatnum,xvector<double>& vvalue) {
  if(vvalue.rows!=(int) vatnum.size()) return FALSE; // nothing to be ordered
  if(vvalue.lrows!=1) return FALSE; // start from 1
  for(uint i=0;i<vatnum.size();i++)
    vvalue[i+1]=GetAtomPettiforScale(vatnum.at(i));
  return TRUE;
}

bool SortAtomsPettiforScale(vector<string> &vsymbol,xvector<int> &vorder,xvector<double> &vvalue) {
  if(vorder.rows!=(int) vsymbol.size()) return FALSE; // nothing to be ordered
  if(vvalue.rows!=(int) vsymbol.size()) return FALSE; // nothing to be ordered
  if(vorder.lrows!=1) return FALSE; // start from 1 .. and contains order from 1
  if(vvalue.lrows!=1) return FALSE; // start from 1
  // build
  for(uint i=1;i<=vsymbol.size();i++) {
    vvalue[i]=GetAtomPettiforScale(vsymbol.at(i-1));
    vorder[i]=i;
  }
  aurostd::sort2(vsymbol.size(),vvalue,vorder);
  vector<string> vsymbol_tmp(vsymbol);
  vsymbol.clear();
  for(uint i=1;i<=vsymbol_tmp.size();i++)
    vsymbol.push_back(vsymbol_tmp.at(vorder[i]-1));
  return TRUE;
}

bool SortAtomsPettiforScale(vector<string> &vsymbol,vector<int> &vorder,vector<double> &vvalue) {
  //   vorder.clear();vvalue.clear();
  //   xvector<int> xvorder(vsymbol.size());
  //   xvector<double> xvvalue(vsymbol.size());
  //   SortAtomsPettiforScale(vsymbol,xvorder,xvvalue);
  //   for(uint i=0;i<vsymbol.size();i++) {
  //     vvalue.push_back(xvvalue(i+1));
  //     vorder.push_back(xvorder(i+1));
  //   }
  vector<string> vsymbol_tmp(vsymbol);
  vorder.clear();vvalue.clear();
  for(uint i=0;i<vsymbol.size();i++) {
    vvalue.push_back(GetAtomPettiforScale(vsymbol.at(i)));
    vorder.push_back(i);
  }
  aurostd::sort(vvalue,vorder);
  vsymbol.clear();
  for(uint i=0;i<vsymbol_tmp.size();i++)
    vsymbol.push_back(vsymbol_tmp.at(vorder.at(i)));
  return TRUE;
}

bool SortAtomsPettiforScale(vector<string> &vsymbol,vector<int> &vorder) {
  vector<double> vvalue;
  SortAtomsPettiforScale(vsymbol,vorder,vvalue);
  return TRUE;
}

bool SortAtomsPettiforScale(vector<string> &vsymbol,vector<double> &vvalue) {
  vector<int> vorder;
  SortAtomsPettiforScale(vsymbol,vorder,vvalue);
  return TRUE;
}

bool SortAtomsPettiforScale(vector<string> &vsymbol) {
  vector<double> vvalue;
  vector<int> vorder;
  SortAtomsPettiforScale(vsymbol,vorder,vvalue);
  return TRUE;
}

// **************************************************************************
// Function GetAtomXrayScatt
// **************************************************************************
double GetAtomXrayScatt(const string& symbol) {
  for(int iat=0;iat<NUM_ELEMENTS;iat++)
    if(symbol==atom_symbol_vec.at(iat) || symbol==atom_name_vec.at(iat))
      return xray_scatt_vec.at(iat);
  return (double) NNN;
}
double GetAtomXrayScatt(const uint& atnum) {
  if(atnum>=xray_scatt_vec.size() || atnum<=0) {
    cerr << "GetAtomXrayScatt out of boundary  atnum=" << atnum << "  xray_scatt_vec.size()=" << xray_scatt_vec.size() << endl;
    return (double) NNN;
  }
  return xray_scatt_vec.at(atnum);
}

// **************************************************************************
// Function GetCompoundAttenuationLenght
// **************************************************************************
double GetCompoundAttenuationLenght(const vector<string>& species,const vector<double>& composition,const double& density) { // density in g/cm^3, return in cm
  if(species.size()!=composition.size()) {
    cerr << "ERROR - GetCompoundAttenuationLenght: species.size()[" << species.size() << "]!=composition.size()[" << composition.size() << "]" << endl;
    exit(0);}
  // cout << "Density=" << density << "<br>" << endl;
  double numerator=0.0,denominator=0.0;
  for(uint i=0;i<species.size();i++) {
    numerator+=composition.at(i)*GetAtomMass(species.at(i))*1000;// from Kg to grams  
    denominator+=density*(composition.at(i)*(GetAtomComptonCrossSection(species.at(i))+GetAtomPhotoelectricCrossSection(species.at(i)))*1e-24);
  }
  //  cout << "numerator=" << numerator << endl;
  //  cout << "denominator=" << denominator << endl;
  return (double) numerator/denominator; // in cm
}

double GetCompoundAttenuationLenght(const deque<string>& _species,const deque<int>& _composition,const double& density) { // density in g/cm^3, return in cm
  vector<double> composition;
  for(uint i=0;i<_composition.size();i++) 
    composition.push_back(double(_composition.at(i)));
  vector<string> species;
  for(uint i=0;i<_species.size();i++) 
    species.push_back(string(_species.at(i)));
  return GetCompoundAttenuationLenght(species,composition,density);
}

// **************************************************************************
// Function XATOM_AlphabetizationSpecies & XATOM_AlphabetizationCompound
// **************************************************************************
string XATOM_AlphabetizationSpecies(string speciesA,string speciesB) {
  string system;
  if(speciesA<=speciesB) system=speciesA+speciesB; else system=speciesB+speciesA;
  return system;
}

string XATOM_AlphabetizationSpecies(vector<string> vspecies_in) {
  stringstream system;
  vector<string> vspecies(vspecies_in);
  aurostd::sort(vspecies);
  for(uint i=0;i<vspecies.size();i++)
    system << vspecies.at(i);
  return system.str();
}

string XATOM_AlphabetizationSpecies(vector<string> vspecies_in,vector<double> vnumbers_in) {
  stringstream system;
  vector<string> vspecies(vspecies_in);
  vector<double> vnumbers(vnumbers_in);
  aurostd::sort(vspecies,vnumbers);
  for(uint i=0;i<vspecies.size();i++)
    system << vspecies.at(i) << vnumbers.at(i);
  return system.str();
}

void XATOM_AlphabetizationSpecies(string& system, vector<string>& vspecies,vector<double>& vnumbers) {
  vspecies.clear();vnumbers.clear();
  KBIN::VASP_SplitAlloySpecies(KBIN::VASP_PseudoPotential_CleanName(system),vspecies,vnumbers);
  stringstream systemstream;
  aurostd::sort(vspecies,vnumbers);
  for(uint i=0;i<vspecies.size();i++)
    systemstream << vspecies.at(i);
  system=systemstream.str();
}

void XATOM_AlphabetizationCompound(string& system, vector<string>& vspecies,vector<double>& vnumbers) {
  vspecies.clear();vnumbers.clear();
  KBIN::VASP_SplitAlloySpecies(KBIN::VASP_PseudoPotential_CleanName(system),vspecies,vnumbers);
  stringstream systemstream;
  aurostd::sort(vspecies,vnumbers);
  for(uint i=0;i<vspecies.size();i++)
    systemstream << vspecies.at(i) << vnumbers.at(i);
  system=systemstream.str();
}

void XATOM_AlphabetizationSpecies(string& system, vector<string>& vspecies) {
  vector<double> vnumbers(vspecies.size());
  XATOM_AlphabetizationSpecies(system,vspecies,vnumbers);
}

void XATOM_AlphabetizationSpecies(string& system) {
  vector<string> vspecies;
  vector<double> vnumbers;
  XATOM_AlphabetizationSpecies(system,vspecies,vnumbers);
}

void XATOM_AlphabetizationCompound(string& system) {
  vector<string> vspecies;
  vector<double> vnumbers;
  XATOM_AlphabetizationCompound(system,vspecies,vnumbers);
}

// **************************************************************************
// Function XATOM_SplitAlloySpecies
// **************************************************************************
uint XATOM_SplitAlloySpecies(string alloy_in, vector<string> &speciesX) {
  string alloy=alloy_in;
  alloy=KBIN::VASP_PseudoPotential_CleanName(alloy); // always some cleaning is good
  alloy=KBIN::VASP_PseudoPotential_CleanName(alloy); // always some cleaning is good
  alloy=aurostd::RemoveNumbers(alloy);              // remove composition
  speciesX.clear();
  for(uint i=0;i<alloy.length();i++) {
    if(alloy[i]>='A' && alloy[i]<='Z') speciesX.push_back("");
    speciesX.at(speciesX.size()-1)+=alloy[i];
  }
  for(uint i=0;i<speciesX.size();i++)
    speciesX.at(i)=aurostd::CleanStringASCII(speciesX.at(i));
  return speciesX.size();
}


uint new_XATOM_SplitAlloySpecies(string alloy_in, vector<string> &speciesX, vector<double> &natomsX) {
  string alloy=alloy_in,alloyn;
  string letters="QWERTYUIOPASDFGHJKLZXCVBNMqwertyuiopasdfghjklzxcvbnm";
  string numbers="0123456789";
  cerr << alloy << endl;// exit(0);

  alloy=KBIN::VASP_PseudoPotential_CleanName(alloy); // always some cleaning is good
  cerr << alloy << endl;// exit(0);
  alloy=KBIN::VASP_PseudoPotential_CleanName(alloy); // always some cleaning is good
  cerr << alloy << endl;// exit(0);
  alloy=aurostd::CleanStringASCII(alloy);
  cerr << alloy << endl;// exit(0);
  for(uint i=0;i<alloy.length();i++)
    for(uint j=0;j<letters.length();j++)
      if(alloy[i]==letters[j] && alloy[i]!=0) {cerr << alloy[i] << endl; alloy[i]='_';}
  
  cerr << alloy << endl; exit(0);

  speciesX.clear();
  for(uint i=0;i<alloy.length();i++) {
    if(alloy[i]>='A' && alloy[i]<='Z') speciesX.push_back("");
    speciesX.at(speciesX.size()-1)+=alloy[i];
  }
  for(uint i=0;i<speciesX.size();i++)
    speciesX.at(i)=aurostd::CleanStringASCII(speciesX.at(i));   // clean it up so it does not have problems inside only letters_numbers
  // now the atoms
  natomsX.clear();
  for(uint i=0;i<speciesX.size();i++) {
    if(i<speciesX.size()-1)
      natomsX.push_back(aurostd::string2utype<double>(alloyn.substr(alloyn.find(speciesX.at(i))+speciesX.at(i).length(),alloyn.find(speciesX.at(i+1))-alloyn.find(speciesX.at(i))-speciesX.at(i).length())));
    else
      natomsX.push_back(aurostd::string2utype<double>(alloyn.substr(alloyn.find(speciesX.at(i))+speciesX.at(i).length())));
    if(abs(natomsX.at(natomsX.size()-1))<=0.00001) natomsX.at(natomsX.size()-1)=1.0;  // fix the no number = 1
  }
  //  for(uint i=0;i<natomsX.size();i++)
  //  cerr << natomsX.at(i) << endl;
  return speciesX.size();
}


uint XATOM_SplitAlloySpecies(string alloy_in, vector<string> &speciesX, vector<double> &natomsX) {
  string alloy=alloy_in,alloyn;
  alloy=KBIN::VASP_PseudoPotential_CleanName(alloy); // always some cleaning is good
  alloy=KBIN::VASP_PseudoPotential_CleanName(alloy); // always some cleaning is good
  alloyn=aurostd::CleanStringASCII(alloy);
  alloy=aurostd::RemoveNumbers(alloy);              // remove composition
  speciesX.clear();
  for(uint i=0;i<alloy.length();i++) {
    if(alloy[i]>='A' && alloy[i]<='Z') speciesX.push_back("");
    speciesX.at(speciesX.size()-1)+=alloy[i];
  }
  for(uint i=0;i<speciesX.size();i++)
    speciesX.at(i)=aurostd::CleanStringASCII(speciesX.at(i));   // clean it up so it does not have problems inside only letters_numbers
  // now the atoms
  natomsX.clear();
  for(uint i=0;i<speciesX.size();i++) {
    if(i<speciesX.size()-1)
      natomsX.push_back(aurostd::string2utype<double>(alloyn.substr(alloyn.find(speciesX.at(i))+speciesX.at(i).length(),alloyn.find(speciesX.at(i+1))-alloyn.find(speciesX.at(i))-speciesX.at(i).length())));
    else
      natomsX.push_back(aurostd::string2utype<double>(alloyn.substr(alloyn.find(speciesX.at(i))+speciesX.at(i).length())));
    if(abs(natomsX.at(natomsX.size()-1))<=0.00001) natomsX.at(natomsX.size()-1)=1.0;  // fix the no number = 1
  }
  //  for(uint i=0;i<natomsX.size();i++)
  //  cerr << natomsX.at(i) << endl;
  return speciesX.size();
}


// **************************************************************************
// Function XATOM_SplitAlloyPseudoPotentials
// **************************************************************************
uint XATOM_SplitAlloyPseudoPotentials(string alloy_in, vector<string> &species_ppX) {
  string alloy=alloy_in;
  alloy=aurostd::RemoveNumbers(alloy);              // remove composition
  species_ppX.clear();
  for(uint i=0;i<alloy.length();i++) {
    if(alloy[i]>='A' && alloy[i]<='Z') species_ppX.push_back("");
    species_ppX.at(species_ppX.size()-1)+=alloy[i];
  }
  for(uint i=0;i<species_ppX.size();i++)
    species_ppX.at(i)=aurostd::CleanStringASCII(species_ppX.at(i));
  return species_ppX.size();
}

uint XATOM_SplitAlloyPseudoPotentials(string alloy_in, vector<string> &species_ppX, vector<double> &natomsX) {
  string alloy=alloy_in,alloyn;
  alloyn=aurostd::CleanStringASCII(alloy);
  alloy=aurostd::RemoveNumbers(alloy);              // remove composition
  species_ppX.clear();
  for(uint i=0;i<alloy.length();i++) {
    if(alloy[i]>='A' && alloy[i]<='Z') species_ppX.push_back("");
    species_ppX.at(species_ppX.size()-1)+=alloy[i];
  }
  for(uint i=0;i<species_ppX.size();i++)
    species_ppX.at(i)=aurostd::CleanStringASCII(species_ppX.at(i));   // clean it up so it does not have problems inside only letters_numbers
  // now the atoms
  natomsX.clear();
  for(uint i=0;i<species_ppX.size();i++) {
    if(i<species_ppX.size()-1)
      natomsX.push_back(aurostd::string2utype<double>(alloyn.substr(alloyn.find(species_ppX.at(i))+species_ppX.at(i).length(),alloyn.find(species_ppX.at(i+1))-alloyn.find(species_ppX.at(i))-species_ppX.at(i).length())));
    else
      natomsX.push_back(aurostd::string2utype<double>(alloyn.substr(alloyn.find(species_ppX.at(i))+species_ppX.at(i).length())));
    if(abs(natomsX.at(natomsX.size()-1))<=0.00001) natomsX.at(natomsX.size()-1)=1.0;  // fix the no number = 1
  }
  //  for(uint i=0;i<natomsX.size();i++)
  //  cerr << natomsX.at(i) << endl;
  return species_ppX.size();
}


// ***************************************************************************
// ***************************************************************************
// ***************************************************************************
// _SYM_OP
// look into aflow.h for the definitions

// constructors
_sym_op::_sym_op() {
  Uc.clear();Uf.clear();           // clear stuff
  ctau.clear();ftau.clear();       // clear stuff
  ctrasl.clear();ftrasl.clear();   // clear stuff
  generator.clear();               // clear stuff
  generator_coefficients.clear();  // clear stuff       //DX 12/6/17 - generator coefficients
  su2_coefficients.clear();        // clear stuff       //DX 1/15/18 - su(2) coefficients on Pauli matrices
  SU2_matrix.clear();              // clear stuff	//DX 1/15/18 - 2x2 complex SU(2) matrix
  str_type="";                     // clear stuff
  str_Hermann_Mauguin="";          // clear stuff
  str_Schoenflies="";              // clear stuff
  angle=0.0;                       // clear stuff
  axis.clear();                    // clear stuff
  quaternion_vector.clear();       // clear stuff	//GEENA
  quaternion_matrix.clear();       // clear stuff	//GEENA
  flag_inversion=FALSE;            // clear stuff
  is_pgroup=FALSE;                 // clear stuff
  is_pgroup_xtal=FALSE;            // clear stuff
  is_pgroupk=FALSE;                // clear stuff
  is_pgroupk_xtal=FALSE;           // clear stuff       // DX 12/5/17
  is_fgroup=FALSE;                 // clear stuff
  is_sgroup=FALSE;                 // clear stuff
  is_agroup=FALSE;                 // clear stuff
  site=0;                          // clear stuff       // DX 8/3/17
  basis_atoms_map.clear();         // clear stuff
  basis_types_map.clear();         // clear stuff
  //DX and CO - START
  basis_map_calculated=FALSE;      // clear stuff
  //DX and CO - END
}

// destructor
_sym_op::~_sym_op() {
  free();
}

void _sym_op::free() {
}

const _sym_op& _sym_op::operator=(const _sym_op& b) {       // operator=
  if(this != &b) {
    free();
    Uc=b.Uc;
    Uf=b.Uf;
    ctau=b.ctau;
    ftau=b.ftau;
    ctrasl=b.ctrasl;
    ftrasl=b.ftrasl;
    str_type=b.str_type;
    str_Hermann_Mauguin=b.str_Hermann_Mauguin;
    str_Schoenflies=b.str_Schoenflies;
    generator=b.generator;
    generator_coefficients=b.generator_coefficients; //DX 12/6/17 - generator coefficients
    su2_coefficients=b.su2_coefficients;             //DX 1/15/18 - su(2) coefficients on Pauli matrices
    SU2_matrix=b.SU2_matrix;                         //DX 1/15/18 - 2x2 complex SU(2) matrix 
    angle=b.angle;
    axis=b.axis;
    quaternion_vector=b.quaternion_vector;	//GEENA
    quaternion_matrix=b.quaternion_matrix;	//GEENA
    flag_inversion=b.flag_inversion;
    is_pgroup=b.is_pgroup;
    is_pgroup_xtal=b.is_pgroup_xtal;
    is_pgroupk=b.is_pgroupk;
    is_pgroupk_xtal=b.is_pgroupk_xtal;  //DX 12/5/17
    is_fgroup=b.is_fgroup;
    is_sgroup=b.is_sgroup;
    is_agroup=b.is_agroup;
    site=b.site;                        // DX 8/3/17
    basis_atoms_map.clear();
    for(uint i=0;i<b.basis_atoms_map.size();i++)
      basis_atoms_map.push_back(b.basis_atoms_map.at(i));
    basis_types_map.clear();
    for(uint i=0;i<b.basis_types_map.size();i++)
      basis_types_map.push_back(b.basis_types_map.at(i));
    //DX and CO - START
    basis_map_calculated=b.basis_map_calculated;
    //DX and CP - START
  }
  return *this;
}

ostream& operator<<(ostream& oss,const _sym_op& symop) {
  xmatrix<double> Uexp(3,3);
  // oss.setf(std::ios::fixed,std::ios::floatfield);
  // oss.precision(10);
  if(symop.is_pgroup==TRUE) oss << " pgroup" << endl;
  if(symop.is_pgroup_xtal==TRUE) oss << " pgroup_xtal" << endl;
  if(symop.is_fgroup==TRUE) oss << " fgroup" << endl;
  if(symop.is_sgroup==TRUE) oss << " sgroup" << endl;
  if(symop.is_agroup==TRUE) oss << " agroup" << endl;
  if(symop.is_pgroupk==TRUE) oss << " pgroupk" << endl;
  if(symop.is_pgroupk_xtal==TRUE) oss << " pgroupk_xtal" << endl;  // DX 12/5/17
  if(symop.is_agroup==TRUE) oss << " site:" << symop.site << endl; // DX 8/3/17
  oss << " type: "<< symop.str_type << endl;
  oss << " Hermann_Mauguin: "<< symop.str_Hermann_Mauguin << endl;
  oss << " Schoenflies: "<< symop.str_Schoenflies << endl;
  //DX and CO - START
  oss << "" << roundoff(symop.Uc,1e-8) << " Uc "<< endl; //COREY roundoff for printing
  //oss << "" << symop.Uc << " Uc "<< endl;
  //DX and CO -END
  Uexp=exp(symop.generator);
  // if(! a.symop_inversion[k] ) oss << "" << max(aurostd::abs(Uexp-Uc))
  // << " error " << endl; else oss << "" <<  max(aurostd::abs(Uexp+Uc)) << " error " << endl;
  //oss << "" << xint(symop.Uf) << " Uf "<< endl; //CO, WILL ZERO OUT 0.999999999999999999
  //CO - START
  xmatrix<int> Uf_int(symop.Uf.lrows,symop.Uf.lcols,symop.Uf.urows,symop.Uf.ucols);
  for (int i=symop.Uf.lrows;i<=symop.Uf.urows;i++){
    for (int j=symop.Uf.lcols;j<=symop.Uf.ucols;j++){
      Uf_int[i][j]=nint(symop.Uf[i][j]);
    }
  }
  oss << "" << Uf_int << " Uf "<< endl;
  //CO - END
  oss << "" << symop.generator << " A=generator U=+-exp(A) [not Uc and -1 if inversion]" << endl;
  //DX 12/6/17
  oss << "" << symop.generator_coefficients << "  so(3) expansion coefficients on Lx, Ly, and Lz basis" << endl;
  //DX 12/6/17
  //DX 1/15/18 - adding SU(2) and su(2); specific xcomplex printing - START
  char buf11_re[80],buf11_im[80],buf12_re[80],buf12_im[80],buf21_re[80],buf21_im[80],buf22_re[80],buf22_im[80];
  string iobuf="%11.4le";
  sprintf(buf11_re,iobuf.c_str(),symop.SU2_matrix(1,1).re);sprintf(buf11_im,iobuf.c_str(),symop.SU2_matrix(1,1).im);
  sprintf(buf12_re,iobuf.c_str(),symop.SU2_matrix(1,2).re);sprintf(buf12_im,iobuf.c_str(),symop.SU2_matrix(1,2).im);
  sprintf(buf21_re,iobuf.c_str(),symop.SU2_matrix(2,1).re);sprintf(buf21_im,iobuf.c_str(),symop.SU2_matrix(2,1).im);
  sprintf(buf22_re,iobuf.c_str(),symop.SU2_matrix(2,2).re);sprintf(buf22_im,iobuf.c_str(),symop.SU2_matrix(2,2).im);
  oss << " (" << buf11_re << "," << buf11_im << ") ";
  oss << " (" << buf12_re << "," << buf12_im << ")" << endl; 
  oss << " (" << buf21_re << "," << buf21_im << ") "; 
  oss << " (" << buf22_re << "," << buf22_im << ")" << "  SU(2) complex matrix [(real,imaginary)]" << endl; 
  // DX - formatting issues with complesx: oss << " " << symop.SU2_matrix << "  SU(2) complex matrix [(real,imaginary)]" << endl;
  char buf1_re[80],buf1_im[80],buf2_re[80],buf2_im[80],buf3_re[80],buf3_im[80];
  sprintf(buf1_re,iobuf.c_str(),symop.su2_coefficients(1).re);sprintf(buf1_im,iobuf.c_str(),symop.su2_coefficients(1).im);
  sprintf(buf2_re,iobuf.c_str(),symop.su2_coefficients(2).re);sprintf(buf2_im,iobuf.c_str(),symop.su2_coefficients(2).im);
  sprintf(buf3_re,iobuf.c_str(),symop.su2_coefficients(3).re);sprintf(buf3_im,iobuf.c_str(),symop.su2_coefficients(3).im);
  oss << " (" << buf1_re << "," << buf1_im << ") ";
  oss << " (" << buf2_re << "," << buf2_im << ") ";
  oss << " (" << buf3_re << "," << buf3_im << ")" << "  su(2) expansion coefficients on Pauli matrices [(real,imaginary)]" << endl;
  // DX - formatting issues with complex: oss << " " << symop.su2_coefficients << "  su(2) expansion coefficients on Pauli matrices [(real,imaginary)]" << endl;
  //DX 1/15/18 - adding SU(2) and su(2); specific xcomplex printing - END
  oss << " "<< symop.angle << "  angle " << endl;
  oss << "" << symop.axis <<  "  axis " << endl;
  //GEENA START HERE
  // Quaternion Output
  // Roundoff used to print zeros instead of values e-17
  oss << "" << roundoff(symop.quaternion_vector, 1e-8) << " quaternion_vector " << endl;
  oss << "" << roundoff(symop.quaternion_matrix, 1e-8) << " quaternion_matrix " << endl;
  //GEENA STOP HERE
  oss << " "<< symop.flag_inversion << "   inversion " << endl;
  //   oss << " "<< symop.str_type << endl;
  //   oss << " "<< symop.str_Hermann_Mauguin << endl;
  //   oss << " "<< symop.str_Schoenflies << endl;
  if(symop.is_fgroup || symop.is_sgroup) oss << "" << symop.ctau << "    ctau " << endl;
  if(symop.is_fgroup || symop.is_sgroup) oss << "" << symop.ftau << "    ftau " << endl;
  if(symop.is_sgroup) oss << "" << symop.ctrasl << "    ctrasl " << endl;
  if(symop.is_sgroup) oss << "" << symop.ftrasl << "    ftrasl " << endl;

  //DX and CO - START
  //if(symop.is_fgroup==TRUE||symop.is_agroup==TRUE) {
  if(symop.basis_map_calculated) {
    //DX and CO - END
    oss << " - ";
    for(uint n=0;n<symop.basis_atoms_map.size();n++)
      oss << symop.basis_atoms_map[n] << " ";
    oss << "   basis_atoms_map" << " ";
    oss << " - ";
    for(uint n=0;n<symop.basis_types_map.size();n++)
      oss << symop.basis_types_map[n] << " ";
    oss << "   basis_types_map" << " ";
    oss << endl;
  }

  return oss;
}


// ***************************************************************************
// ***************************************************************************
// ***************************************************************************
// wyckoffsite_ITC
// look into aflow.h for the definitions
// constructors

wyckoffsite_ITC::wyckoffsite_ITC() {
  coord.clear();
  type="";
  wyckoffSymbol="";
}

// destructor
wyckoffsite_ITC::~wyckoffsite_ITC() {
  free();
}

// empty free
void wyckoffsite_ITC::free() {
}

// operator=
const wyckoffsite_ITC& wyckoffsite_ITC::operator=(const wyckoffsite_ITC& b) {       // operator=
  if(this != &b) {
    free();
    coord=b.coord;
    type=b.type;
    wyckoffSymbol=b.wyckoffSymbol;
  }
  return *this;
}

// copy 
wyckoffsite_ITC::wyckoffsite_ITC(const wyckoffsite_ITC& b) {
  free();
  coord=b.coord;
  type=b.type;
  wyckoffSymbol=b.wyckoffSymbol;
}

// operator <<
ostream& operator<<(ostream& oss,const wyckoffsite_ITC& site) {
  //DX 12/12/17 [OBSOLETE] oss << "wyckoffsite_ITC operator<< " << endl;
  oss << " coord: "<< site.coord << endl;
  oss << " type: "<< site.type << endl;
  oss << " wyckoffSymbol: "<< site.wyckoffSymbol << endl;
  return oss;
}

// ***************************************************************************
// ***************************************************************************
// ***************************************************************************
// STRUCTURE
// look into aflow.h for the definitions

void xstructure::Free() {
}

// Constructors
xstructure::xstructure(string structure_title) {
  iomode=IOAFLOW_AUTO;        // what else could we do right now...
  title=structure_title;
  directory="";
  prototype="";
  info="";
  if(title.length()==0) title="NO_TITLE_GIVEN";
  // num_atoms=0;
  // num_types=0;
  scale=1.0;
  neg_scale=FALSE;
  coord_type[0]=coord_type[1]=0;
  coord_flag=_COORDS_FRACTIONAL_; // _COORDS_FRACTIONAL_ (0) fractional, _COORDS_CARTESIAN_ (1) cartesian.
  isd=FALSE; // !=0 => Selective dynamics, =0 => no selective dynamics.
  lattice.clear();
  lattice(1,1)=lattice(2,2)=lattice(3,3)=1.0;
  a=b=c=1.0;alpha=beta=gamma=90.0;
  klattice=ReciprocalLattice(lattice,scale);
  f2c=trasp(lattice);
  c2f=inverse(trasp(lattice));
  origin.clear();
  // TOLERANCES ------------------------
  equiv_fpos_epsilon=_EQUIV_FPOS_EPS_; // standard but you can change
  // NUM_EACH_TYPE ---------------------
  num_each_type.clear();
  comp_each_type.clear();
  stoich_each_type.clear(); //CO 171025
  // SPECIES ---------------------------
  species.clear();species_pp.clear();species_pp_type.clear();species_pp_version.clear();species_pp_ZVAL.clear();species_pp_vLDAU.clear();species_volume.clear();species_mass.clear();
  is_vasp4_poscar_format=TRUE;
  is_vasp5_poscar_format=FALSE;
  // ATOMS -----------------------------
  atoms.clear();
  // FLAGS -----------------------------
  Niggli_calculated=FALSE;
  Niggli_avoid=FALSE;
  Minkowski_calculated=FALSE;
  Minkowski_avoid=FALSE;
  LatticeReduction_calculated=FALSE;
  LatticeReduction_avoid=FALSE;
  // LATTICE stuff ---------------------
  Standard_Lattice_calculated=FALSE;
  Standard_Lattice_avoid=FALSE;
  Standard_Lattice_primitive=FALSE;
  Standard_Lattice_conventional=FALSE;
  Standard_Lattice_has_failed=FALSE;
  bravais_lattice_type="";
  bravais_lattice_variation_type="";//wahyu mod
  bravais_lattice_system="";
  bravais_lattice_lattice_type="";
  bravais_lattice_lattice_variation_type="";//wahyu mod
  bravais_lattice_lattice_system="";
  pearson_symbol="";
  reciprocal_lattice_type="";
  reciprocal_lattice_variation_type="";//wahyu mod
  //reciprocal_conventional_lattice_type="";
  bravais_superlattice_type="";
  bravais_superlattice_variation_type="";
  bravais_superlattice_system="";
  pearson_symbol_superlattice="";
  // GENERAL PURPOSE LABEL -------------
  label_uint=0;
  label_int=0;
  label_double=0;
  // ORDER PARAMETER -------------------
  order_parameter_structure=FALSE;
  order_parameter_atoms.clear();
  order_parameter_orbit=1; // always orbit of itself
  order_parameter_sum=0;
  // PARTIAL OCCUPATION -------------------
  partial_occupation_flag=FALSE;
  partial_occupation_tol=DEFAULT_PARTIAL_OCCUPATION_TOLERANCE;   // DEFAULT
  partial_occupation_HNF=0;
  partial_occupation_sublattice.clear();
  // FORCES/POSITIONS ------------------
  qm_calculated=FALSE;
  qm_scale=1.0;
  qm_lattice.clear();
  qm_lattice(1,1)=qm_lattice(2,2)=qm_lattice(3,3)=1.0;
  qm_klattice=ReciprocalLattice(qm_lattice,qm_scale);
  qm_f2c=trasp(qm_lattice);
  qm_c2f=inverse(trasp(qm_lattice));
  qm_origin.clear();
  qm_atoms.clear();
  qm_forces.clear();
  qm_forces_write=FALSE;
  qm_positions.clear();
  qm_positions_write=FALSE;
  qm_E_cell=0.0;qm_dE_cell=0.0;qm_H_cell=0.0;qm_PV_cell=0.0;qm_mag_cell=0.0;qm_P=0.0;
  qm_E_atom=0.0;qm_dE_atom=0.0;qm_H_atom=0.0;qm_PV_atom=0.0;qm_mag_atom=0.0; 
  // KPOINTS ---------------------------
  kpoints_k1=0;kpoints_k2=0;kpoints_k3=0;
  kpoints_s1=0;kpoints_s2=0;kpoints_s3=0;
  kpoints_kmax=0;kpoints_kppra=0;
  kpoints_mode=0;
  kpoints_kscheme="";
  //DX and CO - START
  dist_nn_min=AUROSTD_NAN;    //CO
  // SYMMETRY TOLERANCE ----------------------------
  sym_eps=AUROSTD_NAN; //DX
  sym_eps_calculated=false; //DX, this means that it was calculated and set by the symmetry routines
  sym_eps_change_count=0; //DX 2/22/18 - added tolerance count specific to structure
  //DX and CO - END
  // PGROUP ----------------------------
  pgroup.clear();            // just initialize
  pgroup_calculated=FALSE;
  // PGROUP_XTAL ----------------------------
  pgroup_xtal.clear();        // just initialize
  pgroup_xtal_calculated=FALSE;
  crystal_family="";crystal_system="";point_group_crystal_class="";
  point_group_Shoenflies="";point_group_Hermann_Mauguin="";point_group_orbifold="";
  point_group_type="";point_group_order="";point_group_structure="";
  // PGROUPK ----------------------------
  pgroupk.clear();            // just initialize
  pgroupk_calculated=FALSE;
  // PGROUPK_XTAL ----------------------------
  pgroupk_xtal.clear();            // just initialize //DX 12/5/17 - Added pgroupk_xtal
  pgroupk_xtal_calculated=FALSE;                      //DX 12/5/17 - Added pgroupk_xtal
  // FGROUP ----------------------------
  fgroup.clear();            // just initialize
  fgroup_calculated=FALSE;
  // SGROUP ----------------------------
  sgroup_radius=-_calculate_symmetry_default_sgroup_radius_; // symmetry not calculated
  sgroup_radius_dims.clear();
  sgroup.clear();            // just initialize
  sgroup_calculated=FALSE;
  // SITE POINT GROUP ------------------
  agroup_calculated=FALSE;
  for(uint i=0;i<agroup.size();i++)
    agroup.at(i).clear();
  agroup.clear();
  // INEQUIVALENT ATOMS ----------------
  iatoms_calculated=FALSE;
  for(uint i=0;i<iatoms.size();i++)
    iatoms.at(i).clear();
  iatoms.clear();
  // SPACE GROUP WITH PLATON/FINDSYM -----------
  spacegroup="";
  spacegrouplabel="";
  spacegroupnumber=0;
  spacegroupnumberoption=0;
  spacegroupoption="";
  is_spacegroup_platon=FALSE;
  is_spacegroup_findsym=FALSE;
  is_spacegroup_aflow=FALSE;
  // SPACE GROUP ITC (RHT)
  crystal_system_ITC="";  // RHT
  point_group_ITC=""; // RHT
  //DX and CO - START
  bravais_label_ITC='X';  // RHT
  lattice_label_ITC='X'; // RHT
  space_group_ITC=0; // RHT
  //DX and CO - END
  wyckoff_library_entry_ITC=""; // RHT
  wyccar_ITC.clear(); // RHT
  standard_lattice_ITC.clear(); // RHT
  standard_basis_ITC.clear(); // RHT
  wyckoff_sites_ITC.clear();   // RHT
  wyckoff_symbols_ITC.clear(); // RHT
  setting_ITC = 0; //DX 8/30/17 - SGDATA
  origin_ITC.clear(); //DX 8/30/17 - SGDATA
  general_position_ITC.clear(); //DX 8/30/17 - SGDATA
  // GRID ATOMS ------------------------
  grid_atoms_calculated=FALSE;
  grid_atoms_dimsL.clear();
  grid_atoms_dimsH.clear();
  grid_atoms.clear();        // just initialize
  grid_atoms_number=0;
  grid_atoms_sc2pcMap.clear(); //CO 171025
  grid_atoms_pc2scMap.clear(); //CO 171025
  // LIJK OBEJCTS ----------------------
  lijk_calculated=FALSE;
  lijk_table.clear();
  lijk_cpos.clear();
  lijk_fpos.clear();
  lijk_dims.clear();
  // NEIGHBOURS ------------------------
  neighbours_calculated=FALSE;
  neighbours_radius=0.0;
  neighbours_dradius=0.0;
  for(uint i=0;i<neighbours_atoms_func_r_vs_nn.size();i++)
    neighbours_atoms_func_r_vs_nn.at(i).clear();
  neighbours_atoms_func_r_vs_nn.clear();
  for(uint i=0;i<neighbours_atoms_func_num_vs_nn.size();i++)
    neighbours_atoms_func_num_vs_nn.at(i).clear();
  neighbours_atoms_func_num_vs_nn.clear();
  neighbours_func_r_vs_nn.clear();                   // contains function distance vs neighbours (all atoms)
  neighbours_func_num_vs_nn.clear();                    // contains function number vs neighbours (all atoms)
  // OUTPUT/ERROR ----------------------
  Niggli_has_failed=FALSE;
  Minkowski_has_failed=FALSE;
  LatticeReduction_has_failed=FALSE;
  write_lattice_flag=FALSE;
  write_klattice_flag=FALSE;
  write_inequivalent_flag=FALSE;
  write_DEBUG_flag=FALSE;
  error_flag=FALSE;
  error_string="";
  // -----------------------------------
}

// ifstream/istream
xstructure::xstructure(istream& _input,int _iomode) {
  (*this).iomode=_iomode;
  _input >> (*this); }

xstructure::xstructure(ifstream& _input,int _iomode) {
  (*this).iomode=_iomode;
  _input >> (*this); }

xstructure::xstructure(stringstream& __input,int _iomode) {
  (*this).iomode=_iomode;
  stringstream _input(__input.str());
  _input >> (*this); }

xstructure::xstructure(const string& _input,int _iomode) {
  stringstream strstream;
  aurostd::efile2stringstream(_input,strstream); //CO 171025
  //bool found=FALSE;
  //if(!found && aurostd::substring2bool(_input,".bz2")) found=aurostd::bz2file2stringstream(_input,strstream);
  //if(!found && aurostd::substring2bool(_input,".gz"))  found=aurostd::gzfile2stringstream(_input,strstream);
  //if(!found) found=aurostd::file2stringstream(_input,strstream);
  (*this).iomode=_iomode;
  // stringstream _input(_input.str()); _input >> (*this); }
  // cerr << _input << endl;exit(0);
  // cerr << strstream.str() << endl;exit(0);
  strstream >> (*this);}

xstructure::xstructure(const string& url,const string& file,int _iomode) {
  stringstream strstream;
  aurostd::url2stringstream(url+"/"+file,strstream);
  (*this).iomode=_iomode;
  strstream >> (*this);}

void xstructure::Copy(const xstructure& bstr) {
  // All the other stuff not set in the constructor
  iomode=bstr.iomode;
  title=bstr.title;
  directory=bstr.directory;
  prototype=bstr.prototype;
  info=bstr.info;
  //  num_types=bstr.num_types;
  //  num_atoms=bstr.num_atoms;
  scale=bstr.scale;
  neg_scale=bstr.neg_scale;
  strcpy(coord_type,bstr.coord_type);
  coord_flag=bstr.coord_flag;
  isd=bstr.isd;
  lattice=bstr.lattice;
  a=bstr.a;b=bstr.b;c=bstr.c;
  alpha=bstr.alpha;beta=bstr.beta;gamma=bstr.gamma;
  klattice=bstr.klattice;
  origin=bstr.origin;
  f2c=bstr.f2c;
  c2f=bstr.c2f;
  // TOLERANCES ------------------------
  equiv_fpos_epsilon=bstr.equiv_fpos_epsilon;
  // NUM_EACH_TYPE ---------------------
  num_each_type.clear();
  for(uint i=0;i<bstr.num_each_type.size();i++)
    num_each_type.push_back(bstr.num_each_type.at(i));
  comp_each_type.clear();
  for(uint i=0;i<bstr.comp_each_type.size();i++)
    comp_each_type.push_back(bstr.comp_each_type.at(i));
  stoich_each_type.clear(); //CO 171025
  for(uint i=0;i<bstr.stoich_each_type.size();i++) //CO 171025
    stoich_each_type.push_back(bstr.stoich_each_type.at(i)); //CO 171025
  // SPECIES ---------------------------
  species.clear();
  for(uint i=0;i<bstr.species.size();i++)
    species.push_back(bstr.species.at(i));
  species_pp.clear();
  for(uint i=0;i<bstr.species_pp.size();i++)
    species_pp.push_back(bstr.species_pp.at(i));
  species_pp_type.clear();
  for(uint i=0;i<bstr.species_pp_type.size();i++)
    species_pp_type.push_back(bstr.species_pp_type.at(i));
  species_pp_version.clear();
  for(uint i=0;i<bstr.species_pp_version.size();i++)
    species_pp_version.push_back(bstr.species_pp_version.at(i));
  species_pp_ZVAL.clear();
  for(uint i=0;i<bstr.species_pp_ZVAL.size();i++)
    species_pp_ZVAL.push_back(bstr.species_pp_ZVAL.at(i));
  species_pp_vLDAU.clear();
  for(uint i=0;i<bstr.species_pp_vLDAU.size();i++)
    species_pp_vLDAU.push_back(bstr.species_pp_vLDAU.at(i));
  species_volume.clear();
  for(uint i=0;i<bstr.species_volume.size();i++)
    species_volume.push_back(bstr.species_volume.at(i));
  species_mass.clear();
  for(uint i=0;i<bstr.species_mass.size();i++)
    species_mass.push_back(bstr.species_mass.at(i));
  is_vasp4_poscar_format=bstr.is_vasp4_poscar_format;
  is_vasp5_poscar_format=bstr.is_vasp5_poscar_format;
  // FLAGS -----------------------------
  Niggli_calculated=bstr.Niggli_calculated;
  Niggli_avoid=bstr.Niggli_avoid;
  Minkowski_calculated=bstr.Minkowski_calculated;
  Minkowski_avoid=bstr.Minkowski_avoid;
  LatticeReduction_calculated=bstr.LatticeReduction_calculated;
  LatticeReduction_avoid=bstr.LatticeReduction_avoid;
  // LATTICE stuff ---------------------
  Standard_Lattice_calculated=bstr.Standard_Lattice_calculated;
  Standard_Lattice_avoid=bstr.Standard_Lattice_avoid;
  Standard_Lattice_primitive=bstr.Standard_Lattice_primitive;
  Standard_Lattice_conventional=bstr.Standard_Lattice_conventional;
  Standard_Lattice_has_failed=bstr.Standard_Lattice_has_failed;
  bravais_lattice_type=bstr.bravais_lattice_type;
  bravais_lattice_variation_type=bstr.bravais_lattice_variation_type;
  bravais_lattice_system=bstr.bravais_lattice_system;
  bravais_lattice_lattice_type=bstr.bravais_lattice_lattice_type;
  bravais_lattice_lattice_variation_type=bstr.bravais_lattice_lattice_variation_type;
  bravais_lattice_lattice_system=bstr.bravais_lattice_lattice_system;
  pearson_symbol=bstr.pearson_symbol;
  reciprocal_lattice_type=bstr.reciprocal_lattice_type;
  reciprocal_lattice_variation_type=bstr.reciprocal_lattice_variation_type;
  bravais_superlattice_type=bstr.bravais_superlattice_type;
  bravais_superlattice_variation_type=bstr.bravais_superlattice_variation_type;
  bravais_superlattice_system=bstr.bravais_superlattice_system;
  pearson_symbol_superlattice=bstr.pearson_symbol_superlattice;
  // ATOMS -----------------------------
  atoms.clear();
  for(uint i=0;i<bstr.atoms.size();i++)
    atoms.push_back(bstr.atoms.at(i));
  // GENERAL PURPOSE LABEL -------------
  label_uint=bstr.label_uint;
  label_int=bstr.label_int;
  label_double=bstr.label_double;
  // ORDER PARAMETER -------------------
  order_parameter_structure=bstr.order_parameter_structure;
  order_parameter_atoms.clear();
  for(uint i=0;i<bstr.order_parameter_atoms.size();i++)
    order_parameter_atoms.push_back(bstr.order_parameter_atoms.at(i));
  order_parameter_orbit=bstr.order_parameter_orbit;
  order_parameter_sum=bstr.order_parameter_sum;
  // PARTIAL OCCUPATION -------------------
  partial_occupation_flag=bstr.partial_occupation_flag;
  partial_occupation_tol=bstr.partial_occupation_tol;
  partial_occupation_HNF=bstr.partial_occupation_HNF;
  partial_occupation_sublattice.clear();
  for(uint i=0;i<bstr.partial_occupation_sublattice.size();i++)
    partial_occupation_sublattice.push_back(bstr.partial_occupation_sublattice.at(i));
  // FORCES/POSITIONS ------------------
  qm_calculated=bstr.qm_calculated;
  qm_scale=bstr.qm_scale;
  qm_lattice=bstr.qm_lattice;
  qm_klattice=bstr.qm_klattice;
  qm_f2c=bstr.qm_f2c;
  qm_c2f=bstr.qm_c2f;
  qm_origin=bstr.qm_origin;
  qm_atoms.clear();
  for(uint i=0;i<bstr.qm_atoms.size();i++)
    qm_atoms.push_back(bstr.qm_atoms.at(i));
  qm_forces.clear();
  for(uint i=0;i<bstr.qm_forces.size();i++)
    qm_forces.push_back(bstr.qm_forces.at(i));
  qm_forces_write=bstr.qm_forces_write;
  qm_positions.clear();
  for(uint i=0;i<bstr.qm_positions.size();i++)
    qm_positions.push_back(bstr.qm_positions.at(i));
  qm_positions_write=bstr.qm_positions_write;
  qm_E_cell=bstr.qm_E_cell;qm_dE_cell=bstr.qm_dE_cell;qm_H_cell=bstr.qm_H_cell;qm_PV_cell=bstr.qm_PV_cell;qm_mag_cell=bstr.qm_mag_cell;qm_P=bstr.qm_P;
  qm_E_atom=bstr.qm_E_atom;qm_dE_atom=bstr.qm_dE_atom;qm_H_atom=bstr.qm_H_atom;qm_PV_atom=bstr.qm_PV_atom;qm_mag_atom=bstr.qm_mag_atom;
  // KPOINTS ---------------------------
  kpoints_k1=bstr.kpoints_k1;kpoints_k2=bstr.kpoints_k2;kpoints_k3=bstr.kpoints_k3;
  kpoints_s1=bstr.kpoints_s1;kpoints_s2=bstr.kpoints_s2;kpoints_s3=bstr.kpoints_s3;
  kpoints_kmax=bstr.kpoints_kmax;kpoints_kppra=bstr.kpoints_kppra;
  kpoints_mode=bstr.kpoints_mode;
  kpoints_kscheme=bstr.kpoints_kscheme;
  //DX and CO - START
  dist_nn_min=bstr.dist_nn_min;    //CO
  // SYMMETRY TOLERANCE ----------------------------
  sym_eps=bstr.sym_eps; //DX
  sym_eps_calculated=bstr.sym_eps_calculated; //DX
  sym_eps_change_count=bstr.sym_eps_change_count; //DX 2/22/18 - added tolerance count specific to structure
  //DX and CO - END
  // PGROUP ----------------------------
  pgroup.clear();
  for(uint i=0;i<bstr.pgroup.size();i++)
    pgroup.push_back(bstr.pgroup.at(i));
  pgroup_calculated=bstr.pgroup_calculated;
  // PGROUP_XTAL ----------------------------
  pgroup_xtal.clear();
  for(uint i=0;i<bstr.pgroup_xtal.size();i++)
    pgroup_xtal.push_back(bstr.pgroup_xtal.at(i));
  pgroup_xtal_calculated=bstr.pgroup_xtal_calculated;
  crystal_family=bstr.crystal_family;
  crystal_system=bstr.crystal_system;
  point_group_crystal_class=bstr.point_group_crystal_class;
  point_group_Shoenflies=bstr.point_group_Shoenflies;
  point_group_Hermann_Mauguin=bstr.point_group_Hermann_Mauguin;
  point_group_orbifold=bstr.point_group_orbifold;
  point_group_type=bstr.point_group_type;
  point_group_order=bstr.point_group_order;
  point_group_structure=bstr.point_group_structure;
  // PGROUPK ----------------------------
  pgroupk.clear();
  for(uint i=0;i<bstr.pgroupk.size();i++)
    pgroupk.push_back(bstr.pgroupk.at(i));
  pgroupk_calculated=bstr.pgroupk_calculated;
  // PGROUPK_XTAL ----------------------------
  pgroupk_xtal.clear();                                    //DX 12/5/17 - Added pgroupk_xtal
  for(uint i=0;i<bstr.pgroupk_xtal.size();i++)             //DX 12/5/17 - Added pgroupk_xtal
    pgroupk_xtal.push_back(bstr.pgroupk_xtal.at(i));       //DX 12/5/17 - Added pgroupk_xtal
  pgroupk_xtal_calculated=bstr.pgroupk_xtal_calculated;    //DX 12/5/17 - Added pgroupk_xtal
  // FGROUP ----------------------------
  fgroup.clear();
  for(uint i=0;i<bstr.fgroup.size();i++)
    fgroup.push_back(bstr.fgroup.at(i));
  fgroup_calculated=bstr.fgroup_calculated;
  // SGROUP ----------------------------
  sgroup_radius=bstr.sgroup_radius;
  sgroup_radius_dims=bstr.sgroup_radius_dims;
  sgroup.clear();
  for(uint i=0;i<bstr.sgroup.size();i++)
    sgroup.push_back(bstr.sgroup.at(i));
  sgroup_calculated=bstr.sgroup_calculated;
  // SITE POINT GROUP ------------------
  agroup_calculated=bstr.agroup_calculated;
  for(uint i=0;i<agroup.size();i++) agroup.at(i).clear();
  agroup.clear();
  agroup=std::vector<std::vector<_sym_op> > (bstr.agroup.size());
  for(uint i=0;i<bstr.agroup.size();i++)
    for(uint j=0;j<bstr.agroup.at(i).size();j++)
      agroup.at(i).push_back(bstr.agroup.at(i).at(j));
  // INEQUIVALENT ATOMS ----------------
  iatoms_calculated=bstr.iatoms_calculated;
  for(uint i=0;i<iatoms.size();i++) iatoms.at(i).clear();
  iatoms.clear();
  for(uint i=0;i<bstr.iatoms.size();i++) {
    iatoms.push_back(std::vector<int>(0));
    for(uint j=0;j<bstr.iatoms.at(i).size();j++)
      iatoms.at(i).push_back(bstr.iatoms.at(i).at(j));
  }
  // SPACE GROUP WITH PLATON/FINDSYM/AFLOW -----------
  spacegroup=bstr.spacegroup;
  spacegrouplabel=bstr.spacegrouplabel;
  spacegroupoption=bstr.spacegroupoption;
  spacegroupnumber=bstr.spacegroupnumber;
  spacegroupnumberoption=bstr.spacegroupnumberoption;
  is_spacegroup_platon=bstr.is_spacegroup_platon;
  is_spacegroup_findsym=bstr.is_spacegroup_findsym;
  is_spacegroup_aflow=bstr.is_spacegroup_aflow;
  // SPACE GROUP ITC -----------
  crystal_system_ITC=bstr.crystal_system_ITC;  // RHT
  point_group_ITC=bstr.point_group_ITC; // RHT
  bravais_label_ITC=bstr.bravais_label_ITC; // RHT
  lattice_label_ITC=bstr.lattice_label_ITC; // RHT
  space_group_ITC=bstr.space_group_ITC; // RHT
  wyckoff_library_entry_ITC=bstr.wyckoff_library_entry_ITC; // RHT
  wyccar_ITC.clear(); for(uint i=0;i<bstr.wyccar_ITC.size();i++) wyccar_ITC.push_back(bstr.wyccar_ITC.at(i)); // RHT
  standard_lattice_ITC=bstr.standard_lattice_ITC; // RHT
  standard_basis_ITC.clear(); for(uint i=0;i<bstr.standard_basis_ITC.size();i++) standard_basis_ITC.push_back(bstr.standard_basis_ITC.at(i)); // RHT
  wyckoff_sites_ITC.clear(); for(uint i=0;i<bstr.wyckoff_sites_ITC.size();i++) wyckoff_sites_ITC.push_back(bstr.wyckoff_sites_ITC.at(i)); // RHT
  wyckoff_symbols_ITC.clear(); for(uint i=0;i<bstr.wyckoff_symbols_ITC.size();i++) wyckoff_symbols_ITC.push_back(bstr.wyckoff_symbols_ITC.at(i)); // RHT
  setting_ITC=bstr.setting_ITC; //DX 8/30/17 - SGDATA
  origin_ITC=bstr.origin_ITC; //DX 8/30/17 - SGDATA
  general_position_ITC=bstr.general_position_ITC; //DX 8/30/17 - SGDATA
  // GRID ATOMS ------------------------
  grid_atoms_calculated=bstr.grid_atoms_calculated;
  grid_atoms_dimsL=bstr.grid_atoms_dimsL;
  grid_atoms_dimsH=bstr.grid_atoms_dimsH;
  grid_atoms.clear();
  for(uint i=0;i<bstr.grid_atoms.size();i++)
    grid_atoms.push_back(bstr.grid_atoms.at(i));
  grid_atoms_number=bstr.grid_atoms_number;
  grid_atoms_sc2pcMap.clear(); for(uint i=0;i<bstr.grid_atoms_sc2pcMap.size();i++){grid_atoms_sc2pcMap.push_back(bstr.grid_atoms_sc2pcMap[i]);} //CO 171025
  grid_atoms_pc2scMap.clear(); for(uint i=0;i<bstr.grid_atoms_pc2scMap.size();i++){grid_atoms_pc2scMap.push_back(bstr.grid_atoms_pc2scMap[i]);} //CO 171025
  // LIJK OBEJCTS ----------------------
  lijk_calculated=bstr.lijk_calculated;
  lijk_table.clear();
  lijk_cpos.clear();
  lijk_fpos.clear();
  for(uint i=0;i<bstr.lijk_table.size();i++) {
    lijk_table.push_back(bstr.lijk_table.at(i));
    lijk_cpos.push_back(bstr.lijk_cpos.at(i));
    lijk_fpos.push_back(bstr.lijk_fpos.at(i));
  }
  lijk_dims=bstr.lijk_dims;
  // NEIGHBOURS ------------------------
  neighbours_calculated=bstr.neighbours_calculated;
  neighbours_radius=bstr.neighbours_radius;
  neighbours_dradius=bstr.neighbours_dradius;
  //  for(uint i=0;i<neighbours_atoms_func_r_vs_nn.size();i++)
  //   neighbours_atoms_func_r_vs_nn.at(i).clear();
  neighbours_atoms_func_r_vs_nn.clear();
  for(uint i=0;i<bstr.neighbours_atoms_func_r_vs_nn.size();i++)
    neighbours_atoms_func_r_vs_nn.push_back(bstr.neighbours_atoms_func_r_vs_nn.at(i));
  //  for(uint i=0;i<neighbours_atoms_func_num_vs_nn.size();i++)
  //   neighbours_atoms_func_num_vs_nn.at(i).clear();
  neighbours_atoms_func_num_vs_nn.clear();
  for(uint i=0;i<bstr.neighbours_atoms_func_num_vs_nn.size();i++)
    neighbours_atoms_func_num_vs_nn.push_back(bstr.neighbours_atoms_func_num_vs_nn.at(i));
  neighbours_func_r_vs_nn.clear();
  for(uint i=0;i<bstr.neighbours_func_r_vs_nn.size();i++)
    neighbours_func_r_vs_nn.push_back(bstr.neighbours_func_r_vs_nn.at(i));
  neighbours_func_num_vs_nn.clear();
  for(uint i=0;i<bstr.neighbours_func_num_vs_nn.size();i++)
    neighbours_func_num_vs_nn.push_back(bstr.neighbours_func_num_vs_nn.at(i));
  // OUTPUT/ERROR ----------------------
  Niggli_has_failed=bstr.Niggli_has_failed;
  Minkowski_has_failed=bstr.Minkowski_has_failed;
  LatticeReduction_has_failed=bstr.LatticeReduction_has_failed;
  write_lattice_flag=bstr.write_lattice_flag;
  write_klattice_flag=bstr.write_klattice_flag;
  write_inequivalent_flag=bstr.write_inequivalent_flag;
  write_DEBUG_flag=bstr.write_DEBUG_flag;
  error_flag=bstr.error_flag;
  error_string=bstr.error_string;
  // ----------------------------------
}

// copy
xstructure::xstructure(const xstructure& b) {
  //  Free();
  // *this=b;
  Copy(b);
}

// destructor
xstructure::~xstructure() {
  atoms.clear();
  qm_forces.clear();
  qm_positions.clear();
  pgroup.clear();
  pgroup_xtal.clear();
  pgroupk.clear();
  pgroupk_xtal.clear(); //DX 12/5/17 - Added pgroupk_xtal
  fgroup.clear();
  sgroup.clear();
  for(uint i=0;i<agroup.size();i++)
    agroup.at(i).clear();
  agroup.clear();
  for(uint i=0;i<iatoms.size();i++)
    iatoms.at(i).clear();
  iatoms.clear();
  grid_atoms_dimsL.clear(); //CO 171025
  grid_atoms_dimsH.clear(); //CO 171025
  grid_atoms.clear();
  grid_atoms_sc2pcMap.clear(); //CO 171025
  grid_atoms_pc2scMap.clear(); //CO 171025
  lijk_table.clear();
  lijk_cpos.clear();
  lijk_fpos.clear();
  for(uint i=0;i<neighbours_atoms_func_r_vs_nn.size();i++)
    neighbours_atoms_func_r_vs_nn.at(i).clear();
  neighbours_atoms_func_r_vs_nn.clear();
  for(uint i=0;i<neighbours_atoms_func_num_vs_nn.size();i++)
    neighbours_atoms_func_num_vs_nn.at(i).clear();
  neighbours_atoms_func_num_vs_nn.clear();
  neighbours_func_r_vs_nn.clear();
  neighbours_func_num_vs_nn.clear();
 
  wyccar_ITC.clear();                   //(RHT)
  wyckoff_sites_ITC.clear();            //(RHT)
  wyckoff_symbols_ITC.clear();          //(RHT)
  standard_basis_ITC.clear();           //(RHT)
  standard_lattice_ITC.clear();         // RHT
  origin_ITC.clear();                   //DX 8/30/17 - SGDATA
  general_position_ITC.clear();         //DX 8/30/17 - SGDATA
  Free();
}

// copies xtructures: b=a
const xstructure& xstructure::operator=(const xstructure& b) {  // operator=
  if(this!=&b) {
    Free();
    Copy(b);
  }
  return *this;
}

void xstructure::Clear() {
  xstructure _tmp;
  (*this)=_tmp;
}

// **************************************************************************
// Xstructure operator<< OUTPUT_XSTRUCTURE_OUTPUT 
// **************************************************************************
// print an xstructure in a variety of forms
ostream& operator<<(ostream& oss,const xstructure& a) { // operator<<
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  int a_iomode=a.iomode;
  //  DEBUG=TRUE;
  // ----------------------------------------------------------------------
  // PUT DEFAULT
  if(a_iomode==IOAFLOW_AUTO) a_iomode=IOVASP_AUTO; // some default
  // ----------------------------------------------------------------------
  // VASP OUTPUT
  if(a_iomode==IOVASP_AUTO || a_iomode==IOVASP_POSCAR || a_iomode==IOVASP_ABCCAR || a_iomode==IOVASP_WYCKCAR) { // VASP POSCAR
    oss.setf(std::ios::fixed,std::ios::floatfield);
    uint _precision_=14; //was 16 stefano 10 dane
    oss.precision(_precision_);
    if(a_iomode==IOAFLOW_AUTO)   oss << a.title <<endl; // << " (AUTO) " << endl;
    if(a_iomode==IOVASP_AUTO)    oss << a.title <<endl; // << " (AUTO) " << endl;
    if(a_iomode==IOVASP_POSCAR)  oss << a.title <<endl; // << " (POSCAR) " << endl;
    if(a_iomode==IOVASP_ABCCAR)  oss << a.title <<endl; // << " (ABCCAR) " << endl;
    if(a_iomode==IOVASP_WYCKCAR) oss << a.title <<endl; // << " (WYCKCAR) " << endl;
    if(a.neg_scale==FALSE) {
      oss.precision(6);  // DANE
      oss << a.scale << endl;
    } else {
      // oss << a.scale << endl;
      double s=a.scale;
      double vol=s*s*s*GetVol(a.lattice);
      oss.precision(6);  // DANE
      oss << -1*vol << endl;
    }
    oss.precision(_precision_);  // STEFANO to cut/paste from matlab in format long
    // ----------------------------------------------------------------------
    if(a_iomode==IOVASP_POSCAR || a_iomode==IOVASP_AUTO) {
      for(uint i=1;i<=3;i++) {
	for(uint j=1;j<=3;j++) {
	  oss << " ";
	  if(abs(a.lattice(i,j))<10.0) oss << " ";
	  if(!std::signbit(a.lattice(i,j))) oss << " ";
	  oss << a.lattice(i,j) << "";
	}
	oss << endl;
      }
    }
    // ----------------------------------------------------------------------
    if(a_iomode==IOVASP_ABCCAR || a_iomode==IOVASP_WYCKCAR) {
      oss << " ";
      oss.precision(10);  // STEFANO to cut/paste from matlab in format long
      if(abs(a.a)<10.0) oss << " ";
      if(!std::signbit(a.a)) oss << " ";
      oss << a.a << "";
      if(abs(a.b)<10.0) oss << " ";
      if(!std::signbit(a.b)) oss << " ";
      oss << a.b << "";
      if(abs(a.c)<10.0) oss << " ";
      if(!std::signbit(a.c)) oss << " ";
      oss << a.c << "";
      oss.precision(4);  // STEFANO to cut/paste from matlab in format long
      if(abs(a.alpha)<10.0) oss << " ";
      if(!std::signbit(a.alpha)) oss << " ";
      oss << a.alpha << "";
      if(abs(a.beta)<10.0)  oss << " ";
      if(!std::signbit(a.beta))  oss << " ";
      oss << a.beta << "";
      if(abs(a.gamma)<10.0) oss << " ";
      if(!std::signbit(a.gamma)) oss << " ";
      oss << a.gamma << "";
      if(a_iomode==IOVASP_WYCKCAR) oss << "  " << a.spacegroupnumber << "";
      oss << endl;
      oss.precision(_precision_);  // STEFANO to cut/paste from matlab in format long
    }
    // ----------------------------------------------------------------------
    if(a.is_vasp4_poscar_format==TRUE) {
    } // nothing to do

    if(a.is_vasp5_poscar_format==TRUE) {
      for(uint i=0;i<a.species.size();i++)
        oss << a.species.at(i)<< " ";
      oss << endl;
    }
    // ----------------------------------------------------------------------
    for(uint i=0;i<a.num_each_type.size();i++)
      oss << a.num_each_type.at(i) << " ";
    oss << endl;
    if(a.isd) oss << "Selective Dynamics" << endl; // DONE YOYO BUG
    // oss << a.coord_type << endl;

    if(a.coord_flag==_COORDS_FRACTIONAL_) oss << "Direct(" << a.atoms.size() << ") ";
    if(a.coord_flag==_COORDS_CARTESIAN_)  oss << "Cartesian(" << a.atoms.size() << ") ";
    //  if(a.partial_occupation_flag==TRUE)  oss << "Pocc ";
    if(a.order_parameter_structure==TRUE)  oss << "OrderParameter(" << a.order_parameter_atoms.size() << ") ";
    if(1) { // write [A1B2C3D4]
      //      oss << "[";for(uint i=0;i<a.num_each_type.size();i++) {oss << char('A'+i) << a.num_each_type.at(i);}oss << "] ";
      oss << "["; for(uint i=0,k=0;i<a.num_each_type.size();k+=a.num_each_type.at(i),i++) { oss << char(a.atoms.at(k).type+65) << a.num_each_type.at(i);} oss << "] ";

      if(a.partial_occupation_flag==TRUE)  oss << "Partial ";
      if(a.partial_occupation_flag==TRUE) {
	oss.precision(_pocc_precision_);
	oss << "[";for(uint i=0;i<a.comp_each_type.size();i++) {oss << char('A'+i) << a.comp_each_type.at(i);}oss << "] ";
	oss.precision(_precision_);
      }
    }
    // done
    oss << endl;

    for(uint iat=0;iat<a.atoms.size();iat++) {
      oss << " ";
      for(uint j=1;j<=3;j++) {
	//	oss << " ";
	if(a.coord_flag==_COORDS_FRACTIONAL_) {if(abs(a.atoms.at(iat).fpos(j))<10.0) oss << " ";if(!std::signbit(a.atoms.at(iat).fpos(j))) oss << " "; oss << a.atoms.at(iat).fpos(j) << " ";}
	if(a.coord_flag==_COORDS_CARTESIAN_)  {if(abs(a.atoms.at(iat).cpos(j))<10.0) oss << " ";if(!std::signbit(a.atoms.at(iat).cpos(j))) oss << " "; oss << a.atoms.at(iat).cpos(j) << " ";}
      }
      //  cout << aurostd::modulus(a.atoms.at(iat).cpos) << " ";
      if(a.isd==TRUE)
	oss << " " << a.atoms.at(iat).sd[0] << " " << a.atoms.at(iat).sd[1] << " " << a.atoms.at(iat).sd[2];
      if(a.atoms.at(iat).name_is_given==TRUE) {
	oss << " " << a.atoms.at(iat).name << " ";
	for(uint j=a.atoms.at(iat).name.length();j<5;j++) oss << " ";
      }
      if(a.partial_occupation_flag==TRUE) {
	oss.precision(_pocc_precision_);
	//if(a.atoms.at(iat).partial_occupation_flag==FALSE) oss << "-      ";
	//	if(a.atoms.at(iat).partial_occupation_flag==TRUE) oss << a.atoms.at(iat).partial_occupation_value << "  ";// << " (" << iat << "/" << a.partial_occupation_flags.size() << ")";
	oss << "pocc=" << a.atoms.at(iat).partial_occupation_value << "  ";
	oss.precision(_precision_);
      }
      if(a.order_parameter_structure==TRUE) {
	if(a.atoms.at(iat).order_parameter_atom==FALSE) oss << "- ";
	if(a.atoms.at(iat).order_parameter_atom==TRUE) oss << a.atoms.at(iat).order_parameter_value << " ";// << " (" << iat << "/" << a.order_parameter_atoms.size() << ")";
      }
      if(a.write_inequivalent_flag==TRUE) {
	oss << " ";
	// ?	if(i<10) oss << "0";
	oss << iat << "[";
	if(a.atoms.at(iat).equivalent<10) oss << "0";
	oss << a.atoms.at(iat).equivalent << "]";
	if(a.atoms.at(iat).is_inequivalent) {
	  oss <<"*";
	  oss << "_(" << a.atoms.at(iat).num_equivalents << ") "; //<< "  index=" << a.atoms.at(iat).index_iatoms << "  ";
	  //  " v" << a.iatoms.size() << "   burp ";
	  // for(uint jat=0;jat<a.iatoms.size();jat++)  oss << a.iatoms.at(jat).size() << " ";
	}
      }
      if(a.qm_forces_write) {
	if(a.qm_calculated==TRUE)  oss << "F *(";
	if(a.qm_calculated==FALSE) oss << "F  (";
	for(uint j=1;j<=3;j++) {
	  if(abs(a.qm_forces.at(iat)(j))<10.0) oss << " ";
	  if(!std::signbit(a.qm_forces.at(iat)(j))) oss << " ";
	  oss << a.qm_forces.at(iat)(j) << " ";
	}
	oss << ")_   ";
      }
      if(a.qm_positions_write) {
	if(a.qm_calculated==TRUE)  oss << "P *(";
	if(a.qm_calculated==FALSE) oss << "P  (";
	for(uint j=1;j<=3;j++) {
	  if(abs(a.qm_positions.at(iat)(j))<10.0) oss << " ";
	  if(!std::signbit(a.qm_positions.at(iat)(j))) oss << " ";
	  oss << a.qm_positions.at(iat)(j) << " ";
	}
 	oss << ")_   ";
      }
      if(a.write_DEBUG_flag) {
	oss << " s"<<a.atoms.at(iat).spin;
	oss << " n"<<a.atoms.at(iat).number;
	oss << " b"<<a.atoms.at(iat).basis;
	oss << " N("<<a.atoms.at(iat).cleanname;
	oss << " "<<a.atoms.at(iat).atomic_number<<" "<<" ["<<a.atoms.at(iat).type<<"] ";
	oss << " ijk("<<a.atoms.at(iat).ijk(1)<<","<<a.atoms.at(iat).ijk(2)<<","<<a.atoms.at(iat).ijk(3)<<")";
      }
      oss << endl;oss.flush();
    } // iat
    if(a.write_lattice_flag) {
      oss << "DIRECT LATTICE per raw" << endl;
      for(uint i=1;i<=3;i++) {
	for(uint j=1;j<=3;j++) {
	  oss << " ";
	  if(!std::signbit(a.scale*a.lattice(i,j))) oss << " ";
	  oss << a.scale*a.lattice(i,j) << " ";
	}
	oss << endl;
      }
    }
    if(a.write_klattice_flag) {
      oss << "RECIPROCAL LATTICE per raw" << endl;
      for(uint i=1;i<=3;i++) {
	for(uint j=1;j<=3;j++) {
	  oss << " ";
	  if(!std::signbit(a.klattice(i,j))) oss << " ";
	  oss << a.klattice(i,j) << " ";
	}
	oss << endl;
      }
    }
    if(a.write_lattice_flag && a.write_klattice_flag) {
      oss << "ORTOGONALITY (a*b')/2pi=I" << endl;
      xmatrix<double> orto(3,3);
      orto=((a.scale*a.lattice))*trasp(a.klattice)/(2.0*pi);
      for(uint i=1;i<=3;i++) {
	for(uint j=1;j<=3;j++) {
	  oss << " ";
	  if(!std::signbit(orto(i,j))) oss << " ";
	  oss << orto(i,j) << " ";
	}
	oss << endl;
      }
    }
    if(a.write_DEBUG_flag) {
      oss << "kpoints_k1,k2,k3 = " << a.kpoints_k1 << "," << a.kpoints_k2 << "," << a.kpoints_k3 << endl;
      oss << "kpoints_kmax     = " << a.kpoints_kmax << endl;
      oss << "kpoints_kppra    = " << a.kpoints_kppra << endl;
      oss << "kpoints_mode     = " << a.kpoints_mode << endl;
      oss << "kpoints_kscheme  = " << a.kpoints_kscheme << endl;
    }
    if(0 && a.partial_occupation_flag==TRUE) { // KESONG
      oss << "*******************************" << endl;
      for(int hnf=1;hnf<10;hnf++) {
	oss << hnf << " ";
	oss << endl;
      }
    }
    return oss;
  } // END OF VASP
  // ----------------------------------------------------------------------
  //  QUANTUM ESPRESSO OUTPUT
  if(a_iomode==IOQE_AUTO || a_iomode==IOQE_GEOM) { // VASP POSCAR
    oss << "! AFLOW::QE BEGIN " << endl;
    uint depthQE=27;
    if(a_iomode==IOQE_AUTO) oss << "! " << a.title <<endl;//<< " (AUTO)" <<endl;
    if(a_iomode==IOQE_GEOM) oss << "! " << a.title <<endl;//<< " (GEOM)" <<endl;
    oss << aurostd::PaddedPOST("&system",depthQE," ") << " ! // aflow " << endl;
    oss << aurostd::PaddedPOST(" ibrav=0,",depthQE," ") << " ! // free " << endl;
    oss << aurostd::PaddedPOST(" nat="+aurostd::utype2string(a.atoms.size())+",",depthQE) << " ! // a.atoms.size() " << endl;
    oss << aurostd::PaddedPOST(" ntyp="+aurostd::utype2string(a.num_each_type.size()),depthQE) << " ! // a.num_each_type.size() " << endl;
    //oss << aurostd::PaddedPOST(" ntyp="+aurostd::utype2string(a.num_each_type.size())+",",depthQE) << " ! // a.num_each_type.size() " << endl;  //CO 171010
    //oss << aurostd::PaddedPOST(" ecutwfc=_AFLOW_ECUTWFC_,",depthQE," ") << " ! // fix these " << endl;  //CO 171010
    //oss << aurostd::PaddedPOST(" ecutrho=_AFLOW_ECUTRHO_",depthQE," ") << " ! // fix these " << endl;   //CO 171010
    oss << " /" << endl;
    oss.setf(std::ios::fixed,std::ios::floatfield);
    uint _precision_=14; //was 16 stefano 10 dane
    oss.precision(_precision_);
    if(a.coord_flag==_COORDS_FRACTIONAL_) oss << "ATOMIC_POSITIONS (crystal)" << endl;
    if(a.coord_flag==_COORDS_CARTESIAN_)  oss << "ATOMIC_POSITIONS (angstrom)" << endl;
    for(uint iat=0;iat<a.atoms.size();iat++) {
      oss << " ";
      if(a.atoms.at(iat).name_is_given==TRUE) {
	oss << " " << aurostd::PaddedPOST(KBIN::VASP_PseudoPotential_CleanName(a.atoms.at(iat).name),5," ") << " ";
      } else {
	cerr << "QE needs atoms species names" << endl; exit(0);
      }
      for(uint j=1;j<=3;j++) {
	//  oss << " ";
	if(a.coord_flag==_COORDS_FRACTIONAL_) {if(abs(a.atoms.at(iat).fpos(j))<10.0) oss << " ";if(!std::signbit(a.atoms.at(iat).fpos(j))) oss << " "; oss << a.atoms.at(iat).fpos(j) << " ";}
	if(a.coord_flag==_COORDS_CARTESIAN_)  {if(abs(a.atoms.at(iat).cpos(j))<10.0) oss << " ";if(!std::signbit(a.atoms.at(iat).cpos(j))) oss << " "; oss << a.atoms.at(iat).cpos(j) << " ";}
      }
      oss << " ! // " << a.atoms.at(iat).cleanname << " ";
      if(a.write_inequivalent_flag==TRUE) {
	oss << " ";
	//	if(i<10) oss << "0";
	oss << iat << "[";
	if(a.atoms.at(iat).equivalent<10) oss << "0";
	oss << a.atoms.at(iat).equivalent << "]";
	if(a.atoms.at(iat).is_inequivalent) {
	  oss <<"*";
	  oss << "_(" << a.atoms.at(iat).num_equivalents << ") "; //<< "  index=" << a.atoms.at(iat).index_iatoms << "  ";
	  //  " v" << a.iatoms.size() << "   burp ";
	  // for(uint jat=0;jat<a.iatoms.size();jat++)  oss << a.iatoms.at(jat).size() << " ";
	}
      }
      oss << endl;
    }
    // ----------------------------------------------------------------------
    oss.precision(_precision_);  // STEFANO to cut/paste from matlab in format long
    oss << "CELL_PARAMETERS (angstrom)" << endl ;
    {
      for(uint i=1;i<=3;i++) {
	for(uint j=1;j<=3;j++) {
	  oss << " ";
	  if(abs(a.lattice(i,j))<10.0) oss << " ";
	  if(!std::signbit(a.lattice(i,j))) oss << " ";
	  oss << a.lattice(i,j)*a.scale << ""; //DX 2/15/18 - added scaling factor
	}
	oss << endl;
      }
    }
    oss << "# AFLOW::QE END " << endl;
    return oss;
  }

  // ----------------------------------------------------------------------
  //  ABINIT OUTPUT
  if(a_iomode==IOABINIT_AUTO || a_iomode==IOABINIT_GEOM) { // VASP POSCAR
    oss << "# AFLOW::ABINIT BEGIN " << endl;
    uint _precision_=14; //was 16 stefano 10 dane
    oss.precision(_precision_);
    oss.setf(std::ios::fixed,std::ios::floatfield);
    if(a_iomode==IOABINIT_AUTO) oss << "# " << a.title <<endl;//<< " (AUTO)" << endl;
    if(a_iomode==IOABINIT_GEOM) oss << "# " << a.title <<endl;//<< " (GEOM)" << endl;
    oss << "acell   " << double(1) << "   " << double(1) << "   " << double(1) << "  ANGSTR" << endl; // scaling of the primitive vectors, in Bohr.
    for(uint j=1;j<=3;j++) { // each COLUMN of this array is one primitive translation
      if(j==1) oss << "rprim";
      if(j==2) oss << "     ";
      if(j==3) oss << "     ";
      for(uint i=1;i<=3;i++) {
	oss << " ";
	if(abs(a.lattice(i,j))<10.0) oss << " ";
	if(!std::signbit(a.lattice(i,j))) oss << " ";
	oss << a.lattice(i,j) << "";
      }
      oss << endl;
    }
    oss << "natom " << a.atoms.size() << endl;
    oss << "typat ";
    //   for(uint i=0;i<a.num_each_type.size();i++) oss << a.num_each_type.at(i) << " ";  oss << endl;
    for(uint i=0;i<a.atoms.size();i++)
      oss << a.atoms.at(i).type+1 << " ";
    oss << endl;
    if(a.coord_flag==_COORDS_FRACTIONAL_) oss << "xred " << endl;
    if(a.coord_flag==_COORDS_CARTESIAN_) oss << "xangst " << endl;
    for(uint iat=0;iat<a.atoms.size();iat++) {
      oss << "      ";
      for(uint j=1;j<=3;j++) {
	if(a.coord_flag==_COORDS_FRACTIONAL_) {
	  if(abs(a.atoms.at(iat).fpos(j))<10.0) oss << " ";
	  if(!std::signbit(a.atoms.at(iat).fpos(j))) oss << " ";
	  oss << a.atoms.at(iat).fpos(j) << " ";
	}
	if(a.coord_flag==_COORDS_CARTESIAN_) {
	  if(abs(a.atoms.at(iat).cpos(j))<10.0) oss << " ";
	  if(!std::signbit(a.atoms.at(iat).cpos(j))) oss << " ";
	  oss << a.atoms.at(iat).cpos(j) << " ";
	}
      }
      oss << " # " << a.atoms.at(iat).cleanname << " ";
      if(a.write_inequivalent_flag==TRUE) {
	oss << " ";
	//	if(i<10) oss << "0";
	oss << iat << "[";
	if(a.atoms.at(iat).equivalent<10) oss << "0";
	oss << a.atoms.at(iat).equivalent << "]";
	if(a.atoms.at(iat).is_inequivalent) {
	  oss <<"*";
	  oss << "_(" << a.atoms.at(iat).num_equivalents << ") "; //<< "  index=" << a.atoms.at(iat).index_iatoms << "  ";
	  //  " v" << a.iatoms.size() << "   burp ";
	  // for(uint jat=0;jat<a.iatoms.size();jat++)  oss << a.iatoms.at(jat).size() << " ";
	}
      }
      oss << endl;
    }
    oss << "# AFLOW::ABINIT END " << endl;
    return oss;
  } 

  // ----------------------------------------------------------------------
  //  AIMS OUTPUT
  if(a_iomode==IOAIMS_AUTO || a_iomode==IOAIMS_GEOM) { // VASP POSCAR
    uint _precision_=14; //was 16 stefano 10 dane
    oss.precision(_precision_);
    oss.setf(std::ios::fixed,std::ios::floatfield);
    oss << "# " << a.title <<endl;//<< " (AUTO)" << endl;
    oss << "# AFLOW::AIMS BEGIN " << endl;  //come after title
    for(uint i=1;i<=3;i++) { // each COLUMN of this array is one primitive translation
      oss << "lattice_vector ";
      for(uint j=1;j<=3;j++) {
        if(abs(a.lattice(i,j))<100.0) oss << " ";
	if(abs(a.lattice(i,j))<10.0) oss << " ";
        if(!std::signbit(a.lattice(i,j))) oss << " ";
	oss << a.lattice(i,j) << "";
      }
      oss << endl;
    }
    for(uint iat=0;iat<a.atoms.size();iat++) {
      if(0||LDEBUG){
        cerr << "XSTRUCTURE<<: AIMS a.coord_flag==" << a.coord_flag << endl;
        cerr << "XSTRUCTURE<<: AIMS _COORDS_FRACTIONAL_==" << _COORDS_FRACTIONAL_ << endl;
        cerr << "XSTRUCTURE<<: AIMS _COORDS_CARTESIAN_==" << _COORDS_CARTESIAN_ << endl;
      }
      oss << (a.coord_flag==_COORDS_FRACTIONAL_ ? "atom_frac" : "atom") << " ";
      for(uint j=1;j<=3;j++) {
	if(a.coord_flag==_COORDS_FRACTIONAL_) {
          if(abs(a.atoms.at(iat).fpos(j))<100.0) oss << " ";
	  if(abs(a.atoms.at(iat).fpos(j))<10.0) oss << " ";
          if(!std::signbit(a.atoms.at(iat).fpos(j))) oss << " ";
	  oss << a.atoms.at(iat).fpos(j) << " ";
	}
	if(a.coord_flag==_COORDS_CARTESIAN_)  {
          if(abs(a.atoms.at(iat).cpos(j))<100.0) oss << " ";
	  if(abs(a.atoms.at(iat).cpos(j))<10.0) oss << " ";
          if(!std::signbit(a.atoms.at(iat).cpos(j))) oss << " ";
	  oss << a.atoms.at(iat).cpos(j) << " ";
	}
      }
      oss << " " << a.atoms.at(iat).cleanname << " ";
      if(a.write_inequivalent_flag==TRUE) {
        oss << " # ";
	//	if(i<10) oss << "0";
	oss << iat << "[";
	if(a.atoms.at(iat).equivalent<10) oss << "0";
	oss << a.atoms.at(iat).equivalent << "]";
	if(a.atoms.at(iat).is_inequivalent) {
	  oss <<"*";
	  oss << "_(" << a.atoms.at(iat).num_equivalents << ") "; //<< "  index=" << a.atoms.at(iat).index_iatoms << "  ";
	  //  " v" << a.iatoms.size() << "   burp ";
	  // for(uint jat=0;jat<a.iatoms.size();jat++)  oss << a.iatoms.at(jat).size() << " ";
	}
      }
      oss << endl;
    }
    oss << "# AFLOW::AIMS END " << endl;
    return oss;
  } 


  // ----------------------------------------------------------------------
  oss << "NOT CODED YET" << endl;
  return oss;
}

// **************************************************************************
// PrintUNCLE
// **************************************************************************
string xstructure::PrintUNCLE(void) {   // Print in uncle format
  stringstream oss;
  this->iomode=IOVASP_POSCAR;
  oss << "# Structure name:" << endl;
  oss << *this;
  return oss.str();
}

// // **************************************************************************
// // PrintADO
// // **************************************************************************
// string xstructure::PrintADO(string strin) {   // Print in ado format
//   stringstream oss;
//   uint _precision_=14; //was 16 stefano 10 dane
//   oss.setf(std::ios::fixed,std::ios::floatfield);
//   oss.precision(_precision_);
//   // oss << "# Structure name number of atoms in the cell:" << endl;
//   //  oss << "1=" << species.at(0) << " 2=" << species.at(1) << " 3=" << species.at(2) << " n=" << species.size() << endl;
//   oss << strin << "# atomic positions in cartesians" << endl;
//   for(uint iat=0;iat<atoms.size();iat++) {
//     oss << strin;
//     for(uint j=1;j<=3;j++) {
//       if(abs(atoms.at(iat).cpos(j))<10.0) oss << " ";if(atoms.at(iat).cpos(j)>=0.0) oss << " "; oss << atoms.at(iat).cpos(j) << " ";
//     }
//     if(species.size()>0) if(atoms.at(iat).name==species.at(0)) oss << " -1";
//     if(species.size()>1) if(atoms.at(iat).name==species.at(1)) oss << "  1";
//     oss << "  " << atoms.at(iat).name;// << " " << atoms.at(iat).type;
//     oss << endl;
//   }  
//   oss << strin << "# Unit cell vectors in units of a" << endl;
//   for(uint ilat=1;ilat<=3;ilat++) {
//     oss << strin;
//     for(uint j=1;j<=3;j++) {
//       if(abs(lattice[ilat][j])<10.0) oss << " ";if(lattice[ilat][j]>=0.0) oss << " "; oss << lattice[ilat][j] << " ";
//     }
//     oss << endl;
//   }
//   return oss.str();
// }


// **************************************************************************
// PAULING DETECTOR
// **************************************************************************
bool PAULING_WyckoffDetector(vector<string> &vinput);

string DeStupidizer(string &strin) {
  aurostd::StringSubst(strin,"\t"," ");
  string s=" ";for(char c=1;c<32;c++) {s[0]=c; aurostd::StringSubst(strin,s," ");}  //destupidization
  aurostd::StringSubst(strin,"  "," "); // cleaning
  return strin;
}

bool sortAtomsTypes(const _atom& a1,const _atom& a2) {
  return a1.type<a2.type;
  //if(a1.type!=a2.type){return a1.type<a2.type;}
  //double dist1=aurostd::modulus(a1.fpos);
  //double dist2=aurostd::modulus(a2.fpos);
  //return dist1<dist2;
  //prettier (fpos standard)
  //if(a1.fpos.rows!=3 || a1.fpos.rows!=a2.fpos.rows){
  //  cerr << "XSTRUCTURE::sortAtomsNames:: bad cartesian coordinates" << endl;
  //  exit(1);
  //}
  //for(uint i=1;i<=3;i++){
  //  if(a1.fpos[i]!=a2.fpos[i]){return a1.fpos[i]<a2.fpos[i];}
  //}
  //return false;
}

bool sortAtomsNames(const _atom& a1,const _atom& a2) {
  return a1.name<a2.name;
  //if(a1.name!=a2.name){return a1.name<a2.name;}
  //double dist1=aurostd::modulus(a1.fpos);
  //double dist2=aurostd::modulus(a2.fpos);
  //return dist1<dist2;
  //prettier (fpos standard)
  //if(a1.fpos.rows!=3 || a1.fpos.rows!=a2.fpos.rows){
  //  cerr << "XSTRUCTURE::sortAtomsNames:: bad cartesian coordinates" << endl;
  //  exit(1);
  //}
  //for(uint i=1;i<=3;i++){
  //  if(a1.fpos[i]!=a2.fpos[i]){return a1.fpos[i]<a2.fpos[i];}
  //}
  //return false;
}

// **************************************************************************
// Xstructure operator>>  INPUT_XSTRUCTURE_INPUT
// **************************************************************************
// loads istream into a xstructure: istream >> xstructure
istream& operator>>(istream& cinput, xstructure& a) {
#define oss cout
  // this is also a constructor so everything should look well defined
  bool LDEBUG=(FALSE || XHOST.DEBUG);

  if(LDEBUG) cerr << "XSTRUCTURE>>: BEGIN" << endl;
  if(LDEBUG) if(a.iomode==IOAFLOW_AUTO) cerr << "XSTRUCTURE>>: a.iomode = IOAFLOW_AUTO" << endl;
  if(LDEBUG) if(a.iomode==IOVASP_AUTO) cerr << "XSTRUCTURE>>: a.iomode = IOVASP_AUTO" << endl;
  if(LDEBUG) if(a.iomode==IOVASP_POSCAR) cerr << "XSTRUCTURE>>: a.iomode = IOVASP_POSCAR" << endl;
  if(LDEBUG) if(a.iomode==IOVASP_ABCCAR) cerr << "XSTRUCTURE>>: a.iomode = IOVASP_ABCCAR" << endl;
  if(LDEBUG) if(a.iomode==IOVASP_WYCKCAR) cerr << "XSTRUCTURE>>: a.iomode = IOVASP_WYCKCAR" << endl;
  if(LDEBUG) if(a.iomode==IOQE_AUTO) cerr << "XSTRUCTURE>>: a.iomode = IOQE_AUTO" << endl;
  if(LDEBUG) if(a.iomode==IOQE_GEOM) cerr << "XSTRUCTURE>>: a.iomode = IOQE_GEOM" << endl;
  if(LDEBUG) if(a.iomode==IOAIMS_AUTO) cerr << "XSTRUCTURE>>: a.iomode = IOAIMS_AUTO" << endl;  //CO 171008
  if(LDEBUG) if(a.iomode==IOAIMS_GEOM) cerr << "XSTRUCTURE>>: a.iomode = IOAIMS_GEOM" << endl;  //CO 171008
  
  if(LDEBUG) cerr << "XSTRUCTURE>>: definitions" << endl;
  uint iline=0;
  vector<string> vinput,tokens;
  aurostd::stream2vectorstring(cinput,vinput);

  if(LDEBUG) cerr << "XSTRUCTURE>>: DeStupidizer" << endl;
  // now clean for comments, tabs, double spaces ... ect
  for(uint i=1;i<vinput.size();i++) { // not 1st line
    aurostd::string2tokens(vinput.at(i),tokens,"//");  // stefano 
    if(i>0 && tokens.size()>0) vinput.at(i)=tokens.at(0);
    aurostd::string2tokens(vinput.at(i),tokens,"#");  // stefano
    if(i>0 && tokens.size()>0) vinput.at(i)=tokens.at(0);
    aurostd::string2tokens(vinput.at(i),tokens,"!");  // for QE
    if(i>0 && tokens.size()>0) vinput.at(i)=tokens.at(0);
    DeStupidizer(vinput.at(i));
  }
  
  //  for(uint i=0;i<vinput.size();i++) cerr << "[" << i << "] " <<  vinput.at(i) << " " << "[X]" << endl;   exit(0);
  string sstring,stmp;
  bool IOMODE_found=FALSE;
  vector<double> poccaus; // partial occupation local host
  a.partial_occupation_sublattice.clear(); // partial occupation local host

  // ----------------------------------------------------------------------
  // SEARTH FOR MODES 
  //  LDEBUG=TRUE;
  // PAULING PROTO DETECTOR
  if(!IOMODE_found) {
    if(LDEBUG) cerr << "XSTRUCTURE>>: PAULING PROTO DETECTOR" << endl;
    aurostd::string2tokens(vinput.at(0),tokens);
    if(tokens.at(0)=="PAULING" || tokens.at(0)=="pauling" || tokens.at(0)=="Pauling") {
      if(LDEBUG) cerr << "XSTRUCTURE>>: PAULING PROTO DETECTOR = TRUE" << endl;
      IOMODE_found=TRUE;
      PAULING_WyckoffDetector(vinput);
      tokens.clear();
    }
  }
  // QUANTUM ESPRESSO FINDER
  if(!IOMODE_found) {
    if(LDEBUG) cerr << "XSTRUCTURE>>: QUANTUM ESPRESSO DETECTOR" << endl;
    uint QE=0;
    bool QE_ERROR=FALSE;
    if(LDEBUG) for(uint i=0;i<vinput.size();i++) cerr << vinput.at(i) << endl;
    for(uint i=0;i<vinput.size();i++) QE+=aurostd::substring2bool(vinput.at(i),"&system",true)+aurostd::substring2bool(vinput.at(i),"&SYSTEM",true); //DX 1/23/18 - added true to clean the spaces in string
    if(LDEBUG) cerr << "XSTRUCTURE>>: QUANTUM ESPRESSO DETECTOR QE(&system)=" << QE << endl;
    for(uint i=0;i<vinput.size();i++) QE+=aurostd::substring2bool(vinput.at(i),"ibrav=",true)+aurostd::substring2bool(vinput.at(i),"IBRAV=",true); //DX 1/23/18 - added true to clean the spaces in string
    if(LDEBUG) cerr << "XSTRUCTURE>>: QUANTUM ESPRESSO DETECTOR QE(ibrav)=" << QE << endl;
    for(uint i=0;i<vinput.size();i++) QE+=aurostd::substring2bool(vinput.at(i),"nat=",true)+aurostd::substring2bool(vinput.at(i),"NAT=",true); //DX 1/23/18 - added true to clean the spaces in string
    if(LDEBUG) cerr << "XSTRUCTURE>>: QUANTUM ESPRESSO DETECTOR QE(nat)=" << QE << endl;
    for(uint i=0;i<vinput.size();i++) QE+=aurostd::substring2bool(vinput.at(i),"ntyp=",true)+aurostd::substring2bool(vinput.at(i),"NTYP=",true); //DX 1/23/18 - added true to clean the spaces in string
    if(LDEBUG) cerr << "XSTRUCTURE>>: QUANTUM ESPRESSO DETECTOR QE(ntyp)=" << QE << endl;
    for(uint i=0;i<vinput.size();i++) QE+=aurostd::substring2bool(vinput.at(i),"atomic_positions",true)+aurostd::substring2bool(vinput.at(i),"ATOMIC_POSITIONS",true); //DX 1/23/18 - added true to clean the spaces in string
    if(LDEBUG) cerr << "XSTRUCTURE>>: QUANTUM ESPRESSO DETECTOR QE(ATOMIC_POSITIONS)=" << QE << endl;
    
    for(uint i=0;i<vinput.size()&&QE==5;i++) {
      if(aurostd::substring2bool(vinput.at(i),"ATOMIC_POSITIONS") && 
	 !aurostd::substring2bool(vinput.at(i),"crystal","CRYSTAL") && 
	 !aurostd::substring2bool(vinput.at(i),"bohr","BOHR") && //DX added
	 !aurostd::substring2bool(vinput.at(i),"angstrom","ANGSTROM")) {
	cerr << "XSTRUCTURE>>: QE input(1) not supported vinput.at(" << i << ")= \"" << vinput.at(i) << "\"" << endl;QE_ERROR=TRUE;
      }
      if( aurostd::substring2bool(vinput.at(i),"&system")) {
	//	cerr << "XSTRUCTURE>>: QE input(2) not supported vinput.at(" << i << ")= \"" << vinput.at(i) << "\"" << endl;QE_ERROR=TRUE;
      }
    }
    if(QE==5 && QE_ERROR) {cerr << "XSTRUCTURE>>: QE input errors..." << endl; exit(0);}
    if(QE==5 && !QE_ERROR) {
      a.iomode=IOQE_AUTO; // might need further discipline but for now it is ok.. 2013 May SC
      if(LDEBUG) cerr << "XSTRUCTURE>>: QUANTUM ESPRESSO DETECTOR = TRUE" << endl;
      IOMODE_found=TRUE;
    }
  }
  
  //for AIMS input - unfortunately, it's very generic so leave for last
  if(!IOMODE_found) {
    vector<string> tokens_line;
    //since it's so generic, let's be super strict, only look at the first word in the line
    //needs to match atom, atom_frac, or lattice
    for(uint i=0;i<vinput.size()&&!IOMODE_found;i++){
      aurostd::string2tokens(vinput[i],tokens_line," ");
      if(tokens_line.size()){
        IOMODE_found=(IOMODE_found || 
            tokens_line[0]=="atom" || 
            tokens_line[0]=="atom_frac" || 
            tokens_line[0]=="lattice_vector");
      }
    }
    if(IOMODE_found){
      if(LDEBUG){cerr << "XSTRUCTURE>>: AIMS GEOM DETECTOR = TRUE" << endl;}
      a.iomode=IOAIMS_GEOM;
      //if atom_frac, then lattice MUST be provided
      //if atom, no need for lattice
      bool lat_found,lat_found_anywhere,lat1_found,lat2_found,lat3_found;
      bool atom_found,atom_found_anywhere,frac_found_anywhere;
      lat_found=lat_found_anywhere=lat1_found=lat2_found=lat3_found=atom_found=atom_found_anywhere=frac_found_anywhere=false;
      for(uint i=0;i<vinput.size();i++){
        tokens_line.clear();
        aurostd::string2tokens(vinput[i],tokens_line," ");
        if(!tokens_line.size()){continue;}
        lat_found=(tokens_line[0]=="lattice_vector");
        lat_found_anywhere=(lat_found_anywhere || lat_found);
        frac_found_anywhere=(frac_found_anywhere || tokens_line[0]=="atom_frac");
        atom_found=(aurostd::substring2bool(tokens_line[0],"atom"));
        atom_found_anywhere=(atom_found_anywhere || atom_found || frac_found_anywhere);
        if(LDEBUG){
          cerr << "XSTRUCTURE>>: AIMS line[" << i+1 << "]: " << vinput[i] << endl;
          cerr << "XSTRUCTURE>>: AIMS lat_found=" << lat_found << endl;
          cerr << "XSTRUCTURE>>: AIMS lat_found_anywhere=" << lat_found_anywhere << endl;
          cerr << "XSTRUCTURE>>: AIMS lat1_found=" << lat1_found << endl;
          cerr << "XSTRUCTURE>>: AIMS lat2_found=" << lat2_found << endl;
          cerr << "XSTRUCTURE>>: AIMS lat3_found=" << lat3_found << endl;
          cerr << "XSTRUCTURE>>: AIMS atom_found=" << atom_found << endl;
          cerr << "XSTRUCTURE>>: AIMS atom_found_anywhere=" << atom_found_anywhere << endl;
          cerr << "XSTRUCTURE>>: AIMS frac_found_anywhere=" << frac_found_anywhere << endl;
        }
        if(lat_found || atom_found){
          if(lat_found && tokens_line.size()<4){  //could be more, but not less
            cerr << "XSTRUCTURE>>: AIMS input error, ";
            cerr << "lattice_vector ";
            cerr << "at line[" << i+1 << "] is ill-defined" << endl;
            cerr << "line: " << vinput[i] << endl;
            exit(1);
          }
          if(atom_found && tokens_line.size()<5){  //could be more, but not less (need name/type in last column)
            cerr << "XSTRUCTURE>>: AIMS input error, ";
            cerr << "atom position ";
            cerr << "at line[" << i+1 << "] is ill-defined" << endl;
            cerr << "line: " << vinput[i] << endl;
            exit(1);
          }
        }
      }
      if(!atom_found_anywhere){
        cerr << "XSTRUCTURE>>: AIMS input error, no atoms found..." << endl;
        exit(1);
      }
      a.coord_flag=_COORDS_CARTESIAN_;
      if(lat_found_anywhere || frac_found_anywhere){
        a.coord_flag=_COORDS_FRACTIONAL_;
        uint lat_count=0;
        vector<string> lat;
        for(uint i=0;i<vinput.size();i++){
          aurostd::string2tokens(vinput[i],tokens_line," ");
          if(tokens_line[0]=="lattice_vector"){
            lat_count++;
            if(lat_count==1){lat1_found=true;}
            else if(lat_count==2){lat2_found=true;}
            else if(lat_count==3){lat3_found=true;}
            else{
              cerr << "XSTRUCTURE>>: AIMS input error, too many lattice vectors found" << endl;
              exit(1);
            }
          }
        }
        if(!lat1_found || !lat2_found || !lat3_found){
          cerr << "XSTRUCTURE>>: AIMS input error, incomplete lattice vector specification (needed if atom_frac found)" << endl;
          exit(1);
        }
      }
    }
    else{if(LDEBUG){cerr << "XSTRUCTURE>>: AIMS GEOM DETECTOR = FALSE" << endl;}}
  }
  
  // DESPERATE FINDING => VASP
  if(!IOMODE_found) {
    if(LDEBUG) cerr << "XSTRUCTURE>>: VASP DETECTOR" << endl;
    if(a.iomode==IOAFLOW_AUTO) {
      a.iomode=IOVASP_AUTO; // still not found, need something to eat.
      //  if(a.iomode==IOVASP_AUTO) a.iomode=IOVASP_POSCAR;
      IOMODE_found=TRUE;
      if(LDEBUG) cerr << "XSTRUCTURE>>: VASP DETECTOR = TRUE" << endl;
    }
  }
  
  // ----------------------------------------------------------------------
  // SOME EXTRA VERBOSE
  if(LDEBUG) if(a.iomode==IOAFLOW_AUTO) cerr << "XSTRUCTURE>>: a.iomode = IOAFLOW_AUTO" << endl;
  if(LDEBUG) if(a.iomode==IOVASP_AUTO) cerr << "XSTRUCTURE>>: a.iomode = IOVASP_AUTO" << endl;
  if(LDEBUG) if(a.iomode==IOVASP_POSCAR) cerr << "XSTRUCTURE>>: a.iomode = IOVASP_POSCAR" << endl;
  if(LDEBUG) if(a.iomode==IOVASP_ABCCAR) cerr << "XSTRUCTURE>>: a.iomode = IOVASP_ABCCAR" << endl;
  if(LDEBUG) if(a.iomode==IOVASP_WYCKCAR) cerr << "XSTRUCTURE>>: a.iomode = IOVASP_WYCKCAR" << endl;
  if(LDEBUG) if(a.iomode==IOQE_AUTO) cerr << "XSTRUCTURE>>: a.iomode = IOQE_AUTO" << endl;
  if(LDEBUG) if(a.iomode==IOQE_GEOM) cerr << "XSTRUCTURE>>: a.iomode = IOQE_GEOM" << endl;
  if(LDEBUG) if(a.iomode==IOAIMS_AUTO) cerr << "XSTRUCTURE>>: a.iomode = IOAIMS_AUTO" << endl;  //CO 171008
  if(LDEBUG) if(a.iomode==IOAIMS_GEOM) cerr << "XSTRUCTURE>>: a.iomode = IOAIMS_GEOM" << endl;  //CO 171008
  // ----------------------------------------------------------------------
  // VASP INPUT
  if(a.iomode==IOVASP_AUTO || a.iomode==IOVASP_POSCAR || a.iomode==IOVASP_ABCCAR || a.iomode==IOVASP_WYCKCAR) { // VASP POSCAR
    // for variable number of items
    bool scale_second_flag=FALSE;//,scale_third_flag=FALSE;
    double scale_second_value=0.0;//,scale_third_value=0.0;;
    //
    a.is_vasp4_poscar_format=FALSE;
    a.is_vasp5_poscar_format=FALSE;
    // -------------- POSCAR EXISTS
    int num_atoms=0;
    // -------------- TITLE
    // input.getline(stmp,MAX_TITLE_SIZE);title=stmp;
    a.title=vinput.at(iline++);
    // -------------- SCALE
    //    input >> a.scale;
    stmp=vinput.at(iline++);
    aurostd::StringSubst(stmp,"\t"," ");aurostd::StringSubst(stmp,"  "," ");aurostd::StringSubst(stmp,"  "," ");
    aurostd::string2tokens(stmp,tokens);
    if(tokens.size()==0) {cerr << "ERROR: Xstructure operator>> missing second line in poscar" << endl; exit(0);}
    // cerr << tokens.size() <<  " = " << tokens.at(0) << endl;exit(0);
    a.scale=0.0;
    if(tokens.size()>0) {a.scale=aurostd::string2utype<double>(tokens.at(0));}
    if(tokens.size()>1) {scale_second_flag=TRUE;scale_second_value=aurostd::string2utype<double>(tokens.at(1));}
    // if(tokens.size()>2) {scale_third_flag=TRUE;scale_third_value=aurostd::string2utype<double>(tokens.at(2));}
    //  cerr << a.scale << endl;
    // if(scale_second_flag) cerr << scale_second_value << endl;
    // if(scale_third_flag) cerr << scale_third_value << endl;
    // cerr << sstring << endl;
    // -------------- UNIT CELL
    stmp=vinput.at(iline++);
    aurostd::string2tokens(stmp,tokens);
    a.iomode=IOVASP_POSCAR;
    if(tokens.size()==3) a.iomode=IOVASP_POSCAR;
    if(tokens.size()==6) a.iomode=IOVASP_ABCCAR;
    if(tokens.size()==7 || tokens.size()==8) a.iomode=IOVASP_WYCKCAR;
    //    cerr << " LDEBUG token.size()=" << tokens.size() << "" << endl;
    // ---------------------------------------------------------------
    if(a.iomode==IOVASP_POSCAR) {
      // cerr << "XSTRUCTURE>>: Chosen IOVASP_POSCAR" << endl;
      //    input >> a.lattice(1,1) >> a.lattice(1,2) >> a.lattice(1,3);
      a.lattice(1,1)=aurostd::string2utype<double>(tokens[0]);
      a.lattice(1,2)=aurostd::string2utype<double>(tokens[1]);
      a.lattice(1,3)=aurostd::string2utype<double>(tokens[2]);
      stringstream input_tmp;
      input_tmp.clear();input_tmp.str(vinput.at(iline++));
      input_tmp >> a.lattice(2,1) >> a.lattice(2,2) >> a.lattice(2,3);
      input_tmp.clear();input_tmp.str(vinput.at(iline++));
      input_tmp >> a.lattice(3,1) >> a.lattice(3,2) >> a.lattice(3,3);
      xvector<double> data(6);
      data=Getabc_angles(a.lattice,DEGREES);
      a.a=data[1];a.b=data[2];a.c=data[3];
      a.alpha=data[4];a.beta=data[5];a.gamma=data[6];
    }
    // ---------------------------------------------------------------
    if(a.iomode==IOVASP_ABCCAR) {
      // cerr << "XSTRUCTURE>>: Chosen IOVASP_ABCCAR" << endl;
      a.a=aurostd::string2utype<double>(tokens[0]);
      a.b=aurostd::string2utype<double>(tokens[1]);
      a.c=aurostd::string2utype<double>(tokens[2]);
      a.alpha=aurostd::string2utype<double>(tokens[3]);
      a.beta=aurostd::string2utype<double>(tokens[4]);
      a.gamma=aurostd::string2utype<double>(tokens[5]);
      a.lattice=GetClat(a.a,a.b,a.c,a.alpha,a.beta,a.gamma);
    }
    // ---------------------------------------------------------------
    if(a.iomode==IOVASP_WYCKCAR) {
      //     cerr << "XSTRUCTURE>>: Chosen IOVASP_WYCKCAR" << endl;
      //DX 9/5/17 - GET SYM_EPS
      vector<string> title_tokens;
      aurostd::string2tokens(a.title,title_tokens);
      if(aurostd::substring2bool(title_tokens,"sym_eps")){
        a.sym_eps=aurostd::string2utype<double>(title_tokens[title_tokens.size()-1]);     
      } 
      //DX 9/5/17 - GET SYM_EPS
      a.a=aurostd::string2utype<double>(tokens[0]);
      a.b=aurostd::string2utype<double>(tokens[1]);
      a.c=aurostd::string2utype<double>(tokens[2]);
      a.alpha=aurostd::string2utype<double>(tokens[3]);
      a.beta=aurostd::string2utype<double>(tokens[4]);
      a.gamma=aurostd::string2utype<double>(tokens[5]);
      a.spacegroupnumber=aurostd::string2utype<int>(tokens[6]);
      a.spacegroupnumberoption=0;// 1; no option
      if(tokens.size()>=8) a.spacegroupnumberoption=aurostd::string2utype<int>(tokens[7]);
      a.spacegrouplabel=GetSpaceGroupLabel(a.spacegroupnumber);
      a.spacegroup=GetSpaceGroupName(a.spacegroupnumber);
      a.lattice=GetClat(a.a,a.b,a.c,a.alpha,a.beta,a.gamma);
    }
    // If scale < 0 then it should be treated as the volume.
    a.neg_scale=FALSE;
    if(a.scale<0.0) {
      double nvol=-1.0*(a.scale);
      double ovol=GetVol(a.lattice);
      a.scale=std::pow((double) (nvol/ovol),(double) 1.0/3.0);
      a.neg_scale=TRUE;
    }
    a.FixLattices(); // Reciprocal/f2c/c2f
    a.kpoints_k1=0;a.kpoints_k2=0;a.kpoints_k3=0;
    a.kpoints_kmax=0;a.kpoints_kppra=0;a.kpoints_kscheme="";
    clear(a.origin);
    // ---------------------------------------------------------------
    // -------------- CHECK VASP4 VASP5 and CHECK DIRECT/CARTESIANS STUFF
    if(a.iomode==IOVASP_POSCAR) stmp=vinput.at(6);
    if(a.iomode==IOVASP_ABCCAR) stmp=vinput.at(4);
    if(a.iomode==IOVASP_WYCKCAR) stmp=vinput.at(4);
    aurostd::StringSubst(stmp,"\t"," ");aurostd::StringSubst(stmp,"  "," ");aurostd::StringSubst(stmp,"  "," ");
    aurostd::string2tokens(stmp,tokens);
    if(tokens.size()==0) {cerr << "ERROR: Xstructure operator>> missing D/C/S line" << endl; exit(0);}
  
    if(tokens.at(0)[0]=='S' || tokens.at(0)[0]=='s' || tokens.at(0)[0]=='D' || tokens.at(0)[0]=='d' || tokens.at(0)[0]=='C' || tokens.at(0)[0]=='c') {
      // VASP 4
      a.is_vasp4_poscar_format=TRUE;
      a.is_vasp5_poscar_format=FALSE;
    } else {
      // VASP 5
      a.is_vasp4_poscar_format=FALSE;
      a.is_vasp5_poscar_format=TRUE;
    }
    a.partial_occupation_flag=FALSE;

    for(uint j=1;j<tokens.size();j++) if(tokens.at(j)[0]=='P' || tokens.at(j)[0]=='P') a.partial_occupation_flag=TRUE;
    if(scale_second_value==0.0) a.partial_occupation_tol=DEFAULT_PARTIAL_OCCUPATION_TOLERANCE; // nothing defined
    if(a.partial_occupation_flag && scale_second_flag) {
      if(scale_second_value>0.0)  a.partial_occupation_tol=scale_second_value;
      if(scale_second_value<0.0)  a.partial_occupation_HNF=(int) (-scale_second_value);
      if(scale_second_value==0.0) a.partial_occupation_tol=DEFAULT_PARTIAL_OCCUPATION_TOLERANCE; // nothing defined
      //     if(abs(a.partial_occupation_tol)>1e-12) cerr << "a.partial_occupation_tol=" << a.partial_occupation_tol << endl;
      //   if(abs(a.partial_occupation_HNF)>0.0) cerr << "a.partial_occupation_HNF=" << a.partial_occupation_HNF << endl;
    }

    if(a.is_vasp5_poscar_format) stmp=vinput.at(iline++);  // to skip toward vasp5
    // ---------------------------------------------------------------
    // -------------- ATOMS
    // Number of atoms of each type and total number of atoms.
    stmp=vinput.at(iline++);
    // The following is necessary because if the last lattice parameter has
    // tabs/spaces after it then the getline just grabs those spaces.  This
    // second getline will then be called to actually read in the number
    // of each type of atom.
    string tmpns;
    tmpns=aurostd::RemoveWhiteSpaces(stmp);
    if(string(tmpns).size()==0) stmp=vinput.at(iline++);
    a.num_each_type.clear();num_atoms=0;
    aurostd::string2tokens(stmp,tokens);
    // -------------- check for partial occupation
    vector<string> tokens_i=tokens,tokens_j,tokens_k;
    int number;double dpocc;
    //    a.partial_occupation_flag=FALSE;
    for(uint i=0;i<tokens_i.size()&&!a.partial_occupation_flag;i++)
      if(aurostd::substring2bool(tokens_i.at(i),"*") || aurostd::substring2bool(tokens_i.at(i),"+"))
	a.partial_occupation_flag=TRUE;
    // -------------- no partial occupation
    if(a.partial_occupation_flag==FALSE) {
      if(LDEBUG) cerr << "XSTRUCTURE>>: PARTIAL OCCUPATION = FALSE" << endl;
      for(uint i=0;i<tokens_i.size();i++) {
	number=aurostd::string2utype<int>(tokens_i.at(i));dpocc=1.0;
	a.num_each_type.push_back(number); num_atoms+=number;
	for(uint l=0;l<(uint) number;l++) poccaus.push_back(dpocc);
	for(uint l=0;l<(uint) number;l++) a.partial_occupation_sublattice.push_back(_pocc_no_sublattice_);
      }
    }
    // -------------- yes partial occupation
    if(a.partial_occupation_flag==TRUE) {
      if(LDEBUG) cerr << "XSTRUCTURE>>: PARTIAL OCCUPATION = TRUE" << endl;
      for(uint i=0;i<tokens_i.size();i++) {
	if(!aurostd::substring2bool(tokens_i.at(i),"*") && !aurostd::substring2bool(tokens_i.at(i),"+")) { // NO POCC KEYWORD
	  number=aurostd::string2utype<int>(tokens_i.at(i));dpocc=1.0;
	  a.num_each_type.push_back(number); num_atoms+=number;
	  for(uint l=0;l<(uint) number;l++) poccaus.push_back(dpocc);
	  for(uint l=0;l<(uint) number;l++) a.partial_occupation_sublattice.push_back(_pocc_no_sublattice_);
	} else {
	  aurostd::string2tokens(tokens_i.at(i),tokens_j,"+");if(LDEBUG) cerr << "tokens_j.size()=" << tokens_j.size() << endl;
	  number=0;int nnumber;
	  for(uint j=0;j<tokens_j.size();j++) {
	    aurostd::string2tokens(tokens_j.at(j),tokens_k,"*");
	    if(tokens_k.size()==0) {cerr << "ERROR XSTRUCTURE>>: PARTIAL OCCUPATION error [1] tokens_k.size()==0, no *" << endl;exit(0);}
	    nnumber=aurostd::string2utype<int>(tokens_k.at(0));
	    if(tokens_k.size()==1) dpocc=1.0;
	    if(tokens_k.size()==2) {dpocc=aurostd::string2utype<double>(tokens_k.at(1));} // SOMETHING FOR ROUNDING DPOCC UP TO FRACTIONS....
	    number+=nnumber; num_atoms+=nnumber;
	    for(uint l=0;l<(uint) nnumber;l++) poccaus.push_back(dpocc);
	    for(uint l=0;l<(uint) nnumber;l++) {
	      if(tokens_k.size()==1) a.partial_occupation_sublattice.push_back(_pocc_no_sublattice_); // no specie number
	      if(tokens_k.size()==2) a.partial_occupation_sublattice.push_back(i); // put specie number
	    }
	    if(tokens_k.size()>=3) {cerr << "ERROR XSTRUCTURE>>: PARTIAL OCCUPATION error [1] tokens_k.size()>=3, too many *" << endl;exit(0);}
	  } // loop on +
	  a.num_each_type.push_back(number);
	}
      }
      if(LDEBUG) {cerr << "P(" << poccaus.size()<< ") = ";for(uint j=0;j<poccaus.size();j++) cerr << poccaus.at(j) << " ";cerr << endl;}
    }
    //  cerr << "num_atoms=" << num_atoms << endl; exit(0); // num_atoms is the SUM of the atoms in the numbers
    // -------------- COORDINATE TYPES
    // Type of coordinates (Fractional or Cartesian) - only 1st character matters (D/d or C/c).
    // This line might also be the Selective Dynamics line so you must check for that (S/s).
    a.isd=FALSE;
    string stmp;
    if(LDEBUG) cerr << "XSTRUCTURE>>: DEBUG [1]" << endl;
    stmp=vinput.at(iline++);
    aurostd::StringSubst(stmp,"\t"," ");
    std::vector<string> stmp_tokens;
    aurostd::string2tokens(stmp,stmp_tokens);
    // Note that if there are spaces at the beginning of the line we have to remove them.
    if(stmp_tokens.size()==0) {
      oss << "ERROR:  Found blank line on line 7. This line should give coordinate type or selective dynamics." << endl;
      exit(1);
    } else {
      if(LDEBUG) cerr << "XSTRUCTURE>>: DEBUG [2]" << endl;
      string sstmp=stmp_tokens.at(0);
      a.order_parameter_structure=FALSE;
      if(sstmp[0]=='S' || sstmp[0]=='s') {
	a.isd=TRUE;
	stmp=vinput.at(iline++);
	sstmp=aurostd::RemoveSpaces(stmp);
      }
      if(LDEBUG) cerr << "XSTRUCTURE>>: DEBUG [3]" << endl;
      if(sstmp[0]=='D' || sstmp[0]=='d') {
	//	cerr << "FRAC" << endl;
	a.coord_type[0]=sstmp[0];
	a.coord_flag=_COORDS_FRACTIONAL_;
      } else{
	if(sstmp[0]=='C' || sstmp[0]=='c') {
	  //	  cerr << "CART" << endl;
	  a.coord_type[0]=sstmp[0];
	  a.coord_flag=_COORDS_CARTESIAN_;
	  if(a.iomode==IOVASP_WYCKCAR) {
	    oss << "ERROR: WYCKOFF mode requires FRACTIONAL coordinates (DIRECT)." << endl;exit(0);
	  }
	} else{
	  oss << "ERROR:  Did not find coordinate type D/d or C/c." << endl;exit(1);
	}
      }
      a.coord_type[1]='\0';
      if(stmp_tokens.size()>1) {
	for(uint ipos=1;ipos<stmp_tokens.size();ipos++) {
	  if(stmp_tokens.at(ipos)[0]=='O' || stmp_tokens.at(ipos)[0]=='o') a.order_parameter_structure=TRUE;
	  //	  if(stmp_tokens.at(ipos)[0]=='P' || stmp_tokens.at(ipos)[0]=='p') a.partial_occupation_flag=TRUE;
	}
      }
      if(LDEBUG) cerr << "XSTRUCTURE>>: DEBUG a.order_parameter_structure=" << a.order_parameter_structure << endl;
    }
    if(LDEBUG) cerr << "XSTRUCTURE>>: DEBUG [4]" << endl;

    // --------------  Basis atom positions and names
    //   int cnt;
    xvector<double> v(3);clear(v);
    xvector<double> vzero(3);clear(vzero);
    // Set default values: names, sd.
    string def_name="H";
    string def_sd = "TTT";
    // NOTHING TO CLEAR     clear(*cpos);clear(*fpos);clear(*types);
    //  Read in all the atom data.
    // The number of fields on a line can vary but must mean something, which depends on isd.
    // For isd=0:
    //   If nf=3 then there are only atomic positions.
    //   If nf>=4 then there are atomic positions and a name.
    // For isd!=0:
    //   If nf=3 then there are only atomic positions.
    //   If nf=4 then there are atomic positions and a name.
    //   If nf=6 then there are atomic positions and sd.
    //   If nf>=7 then there are atomic positions and sd and a name.

    // clear up a little
    // NOTHING TO CLEAR     for(int i=0;i<=MAX_ATOM_UCELL;i++) (*name)[i].clear();          // vector(0,MAX_ATOM_UCELL)
    // NOTHING TO CLEAR     for(int i=0;i<=MAX_ATOM_UCELL;i++) (*sd)[i].clear();             // vector(0,MAX_ATOM_UCELL)
    // NOTHING TO CLEAR     for(int i=1;i<=MAX_ATOM_UCELL;i++) (*name_is_given)(i)=FALSE; // xvector(1,MAX_ATOM_UCELL)

    int iat=-1;uint itype=0,j;uint iline_ref=iline;
    for(itype=0;itype<a.num_each_type.size();itype++) {
      //  cerr << a.num_each_type.at(itype) << endl;
      for(j=0;j<(uint) a.num_each_type.at(itype);j++) {
	if(LDEBUG) cerr << "XSTRUCTURE>>: DEBUG [5] itype,j=" << itype << "," << j << endl;
	//	bool plug_atom=TRUE;
        iat++; // it startf from -1 so the first is ZERO
        //  cerr << iat << " type=" << itype << endl;
        //  for(int iat=1;iat<=(num_atoms);iat++) {
        // cerr << iat << " " << (num_atoms) << endl;
        _atom atom;                                        // create new atom
        string stmp;
  	stmp=vinput.at(iline++);
	if(LDEBUG) cerr << "XSTRUCTURE>>: " << iline << " " << vinput.size() << "," << iline-iline_ref << "," << num_atoms << "," << stmp << endl;
	if(iline==vinput.size() && (iline-iline_ref<(uint) num_atoms)) {
	  oss << "ERROR:  Unsufficient number of atom lines." << endl;
	  exit(1);
	}  
	if(LDEBUG) cerr << "XSTRUCTURE>>: DEBUG [6]" << endl;
        stmp=aurostd::RemoveCharacter(stmp,'\t');
	aurostd::StringSubst(stmp,"\t"," ");
	std::vector<string> stmp_tokens;
	aurostd::string2tokens(stmp,stmp_tokens);
	if(LDEBUG) cerr << "XSTRUCTURE>>: DEBUG [6b] stmp_tokens.size()=" << stmp_tokens.size() << endl;
	if(stmp_tokens.size()<3) {
	  oss << "ERROR:  Unsufficient number of atom entries in atom=" << iline-iline_ref << "" << endl;
	  exit(1);
	}  

        int id=0;
        string name;
        // Read in the atom positions.
        v(1)=atof(stmp_tokens.at(id++).c_str()); // when not a matter of speed, vectors
        v(2)=atof(stmp_tokens.at(id++).c_str()); // should be indicized with () instead of [],
        v(3)=atof(stmp_tokens.at(id++).c_str()); // b/c () does boundary checks !
        if(a.coord_flag==_COORDS_FRACTIONAL_) {
          //    for(int ii=1;ii<=3;ii++)
          atom.fpos=v;
          atom.cpos=F2C(a.lattice,atom.fpos);
        }
        if(a.coord_flag==_COORDS_CARTESIAN_) {
          atom.cpos=v;
          atom.fpos=C2F(a.lattice,atom.cpos);
        }
        atom.number=iat;    // reference position for convasp
        atom.basis=iat;     // position in the basis
        atom.ijk(1)=0;atom.ijk(2)=0;atom.ijk(3)=0; // inside the zero cell...
        atom.corigin(1)=0.0;atom.corigin(2)=0.0;atom.corigin(3)=0.0; // inside the zero cell
	atom.coord(1)=0.0;atom.coord(2)=0.0;atom.coord(3)=0.0; // inside the zero cell
        atom.spin=0.0;
        atom.noncoll_spin.clear(); //DX 12/5/17 - non-collinear spin
        atom.type=itype;                // CONVASP_MODE if I want type 0,1,2,3,...
	atom.order_parameter_value=0;
	atom.order_parameter_atom=FALSE;
	atom.partial_occupation_value=1.0;
	atom.partial_occupation_flag=FALSE;
	//	plug_atom=TRUE;   // not used
	// NO ORDER PARAMETER
	if(LDEBUG) cerr << "XSTRUCTURE>>: DEBUG [7]" << endl;
	if(a.order_parameter_structure==FALSE) {
	  // stmp_tokens.size() = 4 (plus possible comments).
	  // Read in the names.
	  if(stmp_tokens.size()==4 || (stmp_tokens.size()>=4 && a.isd==0)) {
	    atom.name=stmp_tokens.at(id++);
	    atom.CleanName();
	    atom.CleanSpin();
	    atom.name_is_given=TRUE;
	  }
	  // stmp_tokens.size()=6
	  if(a.isd!=0 && stmp_tokens.size()==6) {
	    string sdt;
	    for(int ic=1;ic<=3;ic++)
	      sdt=sdt+stmp_tokens.at(id++);
	    atom.sd=sdt;
	  }
	  // stmp_tokens.size()=7
	  if(a.isd==TRUE && stmp_tokens.size()>=7) {
	    string sdt;
	    for(int ic=1;ic<=3;ic++)
	      sdt=sdt+stmp_tokens.at(id++);
	    atom.sd=sdt;
	    atom.name=stmp_tokens.at(id++);
	    atom.CleanName();
	    atom.CleanSpin();
	    atom.name_is_given=TRUE;
	  }
	}
	// ORDER PARAMETER
	if(a.order_parameter_structure==TRUE) {
	  if(stmp_tokens.size()!=5 && stmp_tokens.size()!=6) {
	    if(a.order_parameter_structure==TRUE) {
	      cerr << "ERROR with the order parameter you must specify " << endl;
	      cerr << "  x y z Name OrderParameter "                     << endl;
	      cerr << " where x,y,z are the coordinates " << endl;
	      cerr << " Name is the symbol of the atom  " << endl;
	      cerr << " Order parameter is -=none, *=consider, -1,0,1 (integer) values " << endl;
	      cerr << " good luck" << endl;
	      exit(0);
	    }
	  }
	  if(stmp_tokens.size()==5 || stmp_tokens.size()==6) {
	    atom.name=stmp_tokens.at(id);
	    atom.CleanName();atom.CleanSpin();
	    atom.name_is_given=TRUE;
	    id++;
	    if(a.order_parameter_structure==TRUE) {
	      if(stmp_tokens.at(id)=="-") {
		atom.order_parameter_value=0;
		atom.order_parameter_atom=FALSE;
	      } else {
		if(stmp_tokens.at(id)=="*") {
		  atom.order_parameter_value=0;
		  atom.order_parameter_atom=TRUE;
		} else {
		  atom.order_parameter_value=atoi(stmp_tokens.at(id).c_str());
		  atom.order_parameter_atom=TRUE;
		  a.order_parameter_sum+=atom.order_parameter_value;
		}
	      }
	      id++;
	    }
	  }
	}
	if(LDEBUG) cerr << "XSTRUCTURE>>: DEBUG [8]" << endl;
	// now plug atom into the atomlist
        //      cerr <<  atom.name << endl;
        //      cerr <<  atom.cleanname << endl;
        //      cerr <<  atom.name_is_given << endl;
        a.atoms.push_back(atom);
      } // iat loop
    }
  } //  end VASP POSCAR MODE


  // ----------------------------------------------------------------------
  // QE INPUT
  if(a.iomode==IOQE_AUTO || a.iomode==IOQE_GEOM) {
    // LDEBUG=TRUE;
    if(LDEBUG) cerr << "XSTRUCTURE>>: QUANTUM ESPRESSO IOQE_AUTO/GEOM" << endl;
    // START FROM CELL
    a.scale=1.0; // standard
    a.neg_scale=FALSE; // standard
    int ibrav_value = 0; //DX 1/23/18 - added ibrav
    int nat_value = 0;   //DX 1/23/18 - added nat
    bool bohr_lat = false;   //DX 2/15/18 - added bohr for lattice
    bool bohr_pos = false;   //DX 2/15/18 - added bohr for positions
    bool alat = false;   //DX 2/15/18 - added alat
    uint iline,jline;
    iline=vinput.size();
    //DX 1/23/18 - added nat - START
    for(uint i=0;i<vinput.size();i++){ 
      if(aurostd::substring2bool(aurostd::toupper(vinput.at(i)),"NAT=",true)){
        vector<string> comma_tokens; 
        aurostd::string2tokens(aurostd::RemoveWhiteSpaces(vinput[i]),comma_tokens,","); //DX - it is possible to have multiple fields per line
        for(uint j=0;j<comma_tokens.size();j++){
          if(aurostd::substring2bool(aurostd::toupper(comma_tokens.at(j)),"NAT=",true)){
            vector<string> tokens; 
            aurostd::string2tokens(aurostd::RemoveWhiteSpaces(comma_tokens.at(j)),tokens,"=");
            nat_value = aurostd::string2utype<int>(tokens.at(1));
          }
        }
      }                                       
    }
    //DX 1/23/18 - added nat - END
    for(uint i=0;i<vinput.size();i++){
      //DX 1/23/18 [OBSOLETE]  if(aurostd::substring2bool(vinput.at(i),"CELL_PARAMETERS") && 
      //DX 1/23/18 [OBSOLETE]    aurostd::substring2bool(vinput.at(i),"angstrom","ANGSTROM")) {iline=i+1;jline=i;}
      if(aurostd::substring2bool(aurostd::toupper(vinput.at(i)),"CELL_PARAMETERS")){ //DX 1/24/17 - added aurostd::toupper
        if(aurostd::substring2bool(aurostd::toupper(vinput.at(i)),"ANGSTROM")){
          iline=i+1;jline=i;
        }
        else if(aurostd::substring2bool(aurostd::toupper(vinput.at(i)),"ALAT")){ //DX 2/15/18 - added alat bool
          iline=i+1;jline=i;                                                     //DX 2/15/18 - added alat bool
          alat = true;                                                           //DX 2/15/18 - added alat bool
        }
        //DX 1/23/18 - added bohr - START
        else if(aurostd::substring2bool(aurostd::toupper(vinput.at(i)),"BOHR")){
          iline=i+1;jline=i;
          bohr_lat = true;
        }
        //DX 1/23/18 - added bohr - END
      }
      //DX 1/23/18 - added ibrav - START
      else if(aurostd::substring2bool(aurostd::toupper(vinput.at(i)),"IBRAV=",true)){
        vector<string> comma_tokens; 
        aurostd::string2tokens(aurostd::RemoveWhiteSpaces(vinput[i]),comma_tokens,","); //DX - it is possible to have multiple fields per line
        for(uint j=0;j<comma_tokens.size();j++){
          if(aurostd::substring2bool(aurostd::toupper(comma_tokens.at(j)),"IBRAV=",true)){
            vector<string> tokens; 
            aurostd::string2tokens(aurostd::RemoveWhiteSpaces(comma_tokens.at(j)),tokens,"=");
            ibrav_value = aurostd::string2utype<int>(tokens.at(1));
          }
        }
      }
      //DX 1/23/18 - added ibrav - END
    }
      xvector<double> parameters(6); 
      uint celldm_count = 0;
      bool isabc = false; // distinguish between celldm and a,b,c,cosAB,cosAC,cosBC

      // celldm variant
      for(uint i=0;i<vinput.size();i++){
        if(aurostd::substring2bool(aurostd::toupper(vinput.at(i)),"CELLDM")){
          vector<string> comma_tokens; 
          aurostd::string2tokens(aurostd::RemoveWhiteSpaces(vinput[i]),comma_tokens,","); //DX - it is possible to have multiple fields per line
          for(uint j=0;j<comma_tokens.size();j++){
            if(aurostd::substring2bool(aurostd::toupper(comma_tokens.at(j)),"CELLDM",true)){ //DX 1/24/17 - added aurostd::toupper
              celldm_count++;
              vector<string> tokens; 
              aurostd::string2tokens(aurostd::RemoveWhiteSpaces(comma_tokens.at(j)),tokens,"=");
              // In case the celldms are not in order, or not all are given; need to check number
              vector<string> parentheses_tokens; 
              aurostd::string2tokens(aurostd::RemoveWhiteSpaces(tokens.at(0)),parentheses_tokens,"(");
              uint index = aurostd::string2utype<uint>(aurostd::RemoveCharacter(parentheses_tokens.at(1),')'));  
              parameters[index] = aurostd::string2utype<double>(tokens.at(1));
            }
          }
        } 
      }
      // a,b,c,cosAB,cosAC,cosBC variant
      if(celldm_count == 0){
        isabc = true;
	for(uint i=0;i<vinput.size();i++){
	  if(aurostd::substring2bool(aurostd::toupper(aurostd::RemoveWhiteSpaces(vinput.at(i))),"A=")){
	    vector<string> comma_tokens; 
	    aurostd::string2tokens(aurostd::toupper(aurostd::RemoveWhiteSpaces(vinput[i])),comma_tokens,","); //DX - it is possible to have multiple fields per line
	    for(uint j=0;j<comma_tokens.size();j++){
	      if(aurostd::substring2bool(aurostd::toupper(comma_tokens.at(j)),"A=",true)){ //DX 1/24/17 - added aurostd::toupper
                vector<string> tokens; 
                aurostd::string2tokens(aurostd::toupper(aurostd::RemoveWhiteSpaces(comma_tokens.at(j))),tokens,"A=");
                if(tokens.size()==1){
		  parameters[1] = aurostd::string2utype<double>(tokens.at(0)); //A only
                }
	      }
	    }
	  } 
	  if(aurostd::substring2bool(aurostd::toupper(aurostd::RemoveWhiteSpaces(vinput.at(i))),"B=")){
	    vector<string> comma_tokens; 
	    aurostd::string2tokens(aurostd::toupper(aurostd::RemoveWhiteSpaces(vinput[i])),comma_tokens,","); //DX - it is possible to have multiple fields per line
	    for(uint j=0;j<comma_tokens.size();j++){
	      if(aurostd::substring2bool(aurostd::toupper(comma_tokens.at(j)),"B=",true)){ //DX 1/24/17 - added aurostd::toupper
                vector<string> tokens; 
                aurostd::string2tokens(aurostd::toupper(aurostd::RemoveWhiteSpaces(comma_tokens.at(j))),tokens,"B=");
                if(tokens.size()==1){
		  parameters[2] = aurostd::string2utype<double>(tokens.at(0)); //B only
                }
	      }
	    }
	  } 
	  if(aurostd::substring2bool(aurostd::toupper(aurostd::RemoveWhiteSpaces(vinput.at(i))),"C=")){
	    vector<string> comma_tokens; 
	    aurostd::string2tokens(aurostd::toupper(aurostd::RemoveWhiteSpaces(vinput[i])),comma_tokens,","); //DX - it is possible to have multiple fields per line
	    for(uint j=0;j<comma_tokens.size();j++){
	      if(aurostd::substring2bool(aurostd::toupper(comma_tokens.at(j)),"C=",true)){ //DX 1/24/17 - added aurostd::toupper
                vector<string> tokens; 
                aurostd::string2tokens(aurostd::toupper(aurostd::RemoveWhiteSpaces(comma_tokens.at(j))),tokens,"C=");
                if(tokens.size()==1){
		  parameters[3] = aurostd::string2utype<double>(tokens.at(0)); //C only
                }
	      }
	    }
	  } 
	  if(aurostd::substring2bool(aurostd::toupper(aurostd::RemoveWhiteSpaces(vinput.at(i))),"COSAB=")){
	    vector<string> comma_tokens; 
	    aurostd::string2tokens(aurostd::toupper(aurostd::RemoveWhiteSpaces(vinput[i])),comma_tokens,","); //DX - it is possible to have multiple fields per line
	    for(uint j=0;j<comma_tokens.size();j++){
	      if(aurostd::substring2bool(aurostd::toupper(comma_tokens.at(j)),"COSAB=",true)){ //DX 1/24/17 - added aurostd::toupper
                vector<string> tokens; 
                aurostd::string2tokens(aurostd::toupper(aurostd::RemoveWhiteSpaces(comma_tokens.at(j))),tokens,"COSAB=");
                if(tokens.size()==1){
		  parameters[4] = aurostd::string2utype<double>(tokens.at(0)); //COSAB only
                }
	      }
	    }
	  } 
	  if(aurostd::substring2bool(aurostd::toupper(aurostd::RemoveWhiteSpaces(vinput.at(i))),"COSAC=")){
	    vector<string> comma_tokens; 
	    aurostd::string2tokens(aurostd::toupper(aurostd::RemoveWhiteSpaces(vinput[i])),comma_tokens,","); //DX - it is possible to have multiple fields per line
	    for(uint j=0;j<comma_tokens.size();j++){
	      if(aurostd::substring2bool(aurostd::toupper(comma_tokens.at(j)),"COSAC=",true)){ //DX 1/24/17 - added aurostd::toupper
                vector<string> tokens; 
                aurostd::string2tokens(aurostd::toupper(aurostd::RemoveWhiteSpaces(comma_tokens.at(j))),tokens,"COSAC=");
                if(tokens.size()==1){
		  parameters[5] = aurostd::string2utype<double>(tokens.at(0)); //COSAC only
                }
	      }
	    }
	  } 
	  if(aurostd::substring2bool(aurostd::toupper(aurostd::RemoveWhiteSpaces(vinput.at(i))),"COSBC=")){
	    vector<string> comma_tokens; 
	    aurostd::string2tokens(aurostd::toupper(aurostd::RemoveWhiteSpaces(vinput[i])),comma_tokens,","); //DX - it is possible to have multiple fields per line
	    for(uint j=0;j<comma_tokens.size();j++){
	      if(aurostd::substring2bool(aurostd::toupper(comma_tokens.at(j)),"COSBC=",true)){ //DX 1/24/17 - added aurostd::toupper
                vector<string> tokens; 
                aurostd::string2tokens(aurostd::toupper(aurostd::RemoveWhiteSpaces(comma_tokens.at(j))),tokens,"COSBC=");
                if(tokens.size()==1){
		  parameters[6] = aurostd::string2utype<double>(tokens.at(0)); //COSBC only
                }
	      }
	    }
	  } 
	}
      }
    //DX 1/23/18 - added ibrav/parameters - START
    if(ibrav_value!=0){
      a.lattice = pflow::QE_ibrav2lattice(ibrav_value, parameters, isabc);
    }
    //DX 1/23/18 - added ibrav/parameters - END
    else{
    if(iline<vinput.size()-3) { // IN a1 / a2 / a3 rows version
      if(LDEBUG) cerr << "XSTRUCTURE>>: QUANTUM ESPRESSO FOUND 3 extra lines after, trying a1/a2/a3 on raws" << endl;
      stringstream input_tmp;
      input_tmp.clear();input_tmp.str(vinput.at(iline++));
      input_tmp >> a.lattice(1,1) >> a.lattice(1,2) >> a.lattice(1,3);
      input_tmp.clear();input_tmp.str(vinput.at(iline++));
      input_tmp >> a.lattice(2,1) >> a.lattice(2,2) >> a.lattice(2,3);
      input_tmp.clear();input_tmp.str(vinput.at(iline++));
      input_tmp >> a.lattice(3,1) >> a.lattice(3,2) >> a.lattice(3,3);
      if(bohr_lat){                     //DX 1/23/18 - added bohr 
	  a.lattice = a.lattice*bohr2angstrom; //DX 1/23/18 - added bohr
	}                             //DX 1/23/18 - added bohr
      if(alat){                     //DX 2/15/18 - added alat 
        if(isabc){
          a.lattice = a.lattice*parameters[1]; //DX 1/23/18 - added alat
        }
        else{
          a.lattice = a.lattice*parameters[1]*bohr2angstrom; //DX 1/23/18 - added alat
        }
      }                             //DX 2/15/18 - added alat
      xvector<double> data(6);
      data=Getabc_angles(a.lattice,DEGREES);
      a.a=data[1];a.b=data[2];a.c=data[3];
      a.alpha=data[4];a.beta=data[5];a.gamma=data[6];
    }
    }
    // ---------------------------------------------------------------
    // SOME FOR ABCCAR
    // ...
    // ---------------------------------------------------------------
    // SOME FOR WYCCAR
    // ...
    a.FixLattices(); // Reciprocal/f2c/c2f
    a.kpoints_k1=0;a.kpoints_k2=0;a.kpoints_k3=0;
    a.kpoints_kmax=0;a.kpoints_kppra=0;a.kpoints_kscheme="";
    clear(a.origin);
    
    // cerr << "iline=" << iline << endl;
    // cerr << "vinput.size()=" << vinput.size() << endl;
    // NOW ADD ATOMS
    for(uint i=0;i<vinput.size();i++) {
      //DX 1/24/18 [OBSOLETE] if(aurostd::substring2bool(vinput.at(i),"ATOMIC_POSITIONS") && aurostd::substring2bool(vinput.at(i),"crystal","CRYSTAL")) {
      if(aurostd::substring2bool(aurostd::toupper(vinput.at(i)),"ATOMIC_POSITIONS") && aurostd::substring2bool(aurostd::toupper(vinput.at(i)),"CRYSTAL")) { //DX 1/24/18 - added aurostd::toupper
	iline=i+1;a.coord_flag=_COORDS_FRACTIONAL_;
        jline=iline+nat_value; //DX 1/23/18 - added nat
      }
      //DX 1/24/18 [OBSOLETE] if(aurostd::substring2bool(vinput.at(i),"ATOMIC_POSITIONS") && aurostd::substring2bool(vinput.at(i),"angstrom","ANGSTROM")) {
      if(aurostd::substring2bool(aurostd::toupper(vinput.at(i)),"ATOMIC_POSITIONS") && aurostd::substring2bool(aurostd::toupper(vinput.at(i)),"ANGSTROM")) { //DX 1/24/18 - added aurostd::toupper
	iline=i+1;a.coord_flag=_COORDS_CARTESIAN_;
        jline=iline+nat_value; //DX 1/23/18 - added nat
      }
      //DX 1/24/18 -- added Bohr case - START
      if(aurostd::substring2bool(aurostd::toupper(vinput.at(i)),"ATOMIC_POSITIONS") && aurostd::substring2bool(aurostd::toupper(vinput.at(i)),"BOHR")) { //DX 1/24/18 - added aurostd::toupper
	iline=i+1;a.coord_flag=_COORDS_CARTESIAN_;
        jline=iline+nat_value; //DX 1/23/18 - added nat
        bohr_pos = true;
      }
      //DX 1/24/18 -- added Bohr case - END
    }

    for(uint i=iline;i<vinput.size()&&i<jline;i++) {
      _atom atom;                                        // create new atom
      string stmp;
      xvector<double> v(3);clear(v);
      stmp=vinput.at(i);       stmp=aurostd::RemoveCharacter(vinput.at(i),'\t');aurostd::StringSubst(stmp,"\t"," ");
      std::vector<string> stmp_tokens;
      aurostd::string2tokens(stmp,stmp_tokens);
      if(LDEBUG) cerr << "XSTRUCTURE>>: DEBUG [6b] stmp_tokens.size()=" << stmp_tokens.size() << endl;
      if(stmp_tokens.size()<4) {
	oss << "ERROR:  Unsufficient number of atom entries in atom=" << i << "" << endl;
	exit(0);
      }  
      int id=0;
      string name;
      // Read in the atom positions.
      name=stmp_tokens.at(id++);
      v(1)=aurostd::string2utype<double>(stmp_tokens.at(id++)); // when not a matter of speed, vectors
      v(2)=aurostd::string2utype<double>(stmp_tokens.at(id++)); // should be indicized with () instead of [],
      v(3)=aurostd::string2utype<double>(stmp_tokens.at(id++)); // b/c () does boundary checks !
      if(a.coord_flag==_COORDS_FRACTIONAL_) {
	atom.fpos=v;
	atom.cpos=F2C(a.lattice,atom.fpos);
      }
      if(a.coord_flag==_COORDS_CARTESIAN_) {
        if(bohr_pos){ //DX 1/24/18 - added bohr case
	  atom.cpos=bohr2angstrom*v; //DX 1/24/18 - added bohr case
	  atom.fpos=C2F(a.lattice,atom.cpos);
        }
        else{ //DX 1/24/18 - added bohr case
	atom.cpos=v;
	atom.fpos=C2F(a.lattice,atom.cpos);
      }
      }
      atom.name=name;
      atom.CleanName();
      atom.CleanSpin();
      atom.name_is_given=TRUE;
       
      atom.number=a.atoms.size();    // reference position for convasp
      atom.basis=a.atoms.size();     // position in the basis
      atom.ijk(1)=0;atom.ijk(2)=0;atom.ijk(3)=0; // inside the zero cell...
      atom.corigin(1)=0.0;atom.corigin(2)=0.0;atom.corigin(3)=0.0; // inside the zero cell
      atom.coord(1)=0.0;atom.coord(2)=0.0;atom.coord(3)=0.0; // inside the zero cell
      atom.spin=0.0;
      atom.noncoll_spin.clear(); //DX 12/5/17 - non-collinear spin
      // FIXED BELOW atom.type=itype;                // CONVASP_MODE if I want type 0,1,2,3,...
      atom.order_parameter_value=0;
      atom.order_parameter_atom=FALSE;
      atom.partial_occupation_value=1.0;
      atom.partial_occupation_flag=FALSE;
      // DONE
      a.AddAtom(atom);
      // NO PARTIAL OCCUPATION
      a.partial_occupation_sublattice.push_back(_pocc_no_sublattice_);
    }
    // FIX TITLE
    // FIX ITYPE
    // for(uint i=0;i<a.atoms.size();i++) cerr << "a.atoms.at(i).type=" << a.atoms.at(i).type << endl;
    a.SpeciesPutAlphabetic(); // fight analphabetization
    uint iat=0;
    for(uint itype=0;itype<a.num_each_type.size();itype++) 
      for(uint j=0;j<(uint) a.num_each_type.at(itype);j++) {
	if(j==0) a.title+=a.atoms.at(iat).name+aurostd::utype2string(a.num_each_type.at(itype));
	a.atoms.at(iat++).type=itype;
      }
    // for(uint i=0;i<a.atoms.size();i++) cerr << "a.atoms.at(i).type=" << a.atoms.at(i).type << endl;
    a.partial_occupation_flag=FALSE;
    //    a.xstructure2vasp();
    a.is_vasp4_poscar_format=FALSE;
    a.is_vasp5_poscar_format=FALSE;
    // DONE ?
  } // QE INPUT
  
  if(a.iomode==IOAIMS_AUTO || a.iomode==IOAIMS_GEOM) { // AIMS GEOM
    //CO 171008
    //remember, we already did all the debugging above
    //if we get here, we can assume geometry.in is solid!
    //if more debugging needed, test above where we detect iomode
    if(LDEBUG) cerr << "XSTRUCTURE>>: AIMS IOAIMS_AUTO/IOAIMS_GEOM" << endl;
    // START FROM CELL
    a.scale=1.0; // standard
    a.neg_scale=FALSE; // standard
    a.lattice=aurostd::eyes<double>();
    
    uint lat_count=0;
    //get lattice first, if available (c2f,f2c)
    for(uint i=0;i<vinput.size();i++){
      aurostd::string2tokens(vinput[i],tokens," ");
      if(tokens[0]=="lattice_vector" && tokens.size()>3){
        lat_count++;
        a.lattice(lat_count,1)=aurostd::string2utype<double>(tokens[1]);
        a.lattice(lat_count,2)=aurostd::string2utype<double>(tokens[2]);
        a.lattice(lat_count,3)=aurostd::string2utype<double>(tokens[3]);
      }
      if(lat_count==3){break;}
    }
    
    a.FixLattices();
    //a.f2c=trasp(a.lattice), a.c2f=inverse(a.f2c);
    
    if(LDEBUG){
      cerr << "XSTRUCTURE>>: AIMS lattice" << endl;
      cerr << a.lattice << endl;
      cerr << "XSTRUCTURE>>: AIMS f2c" << endl;
      cerr << a.f2c << endl;
      cerr << "XSTRUCTURE>>: AIMS c2f" << endl;
      cerr << a.c2f << endl;
    }

    deque<_atom> atoms;
    //now get atoms
    bool atom_found;
    _atom atom;
    bool atom_props_search=false; //AIMS stores more atom properties UNDER "atom" and before next "atom"
    for(uint i=0;i<vinput.size();i++){
      atom_found=false;
      aurostd::string2tokens(vinput[i],tokens," ");
      if(!tokens.size()){continue;}
      if(aurostd::substring2bool(tokens[0],"atom") && tokens.size()>4){
        atom_props_search=false;
        if(tokens[0]=="atom"){
          atom_found=true;atom.clear();
          atom.cpos[1]=aurostd::string2utype<double>(tokens[1]);
          atom.cpos[2]=aurostd::string2utype<double>(tokens[2]);
          atom.cpos[3]=aurostd::string2utype<double>(tokens[3]);
          atom.fpos=a.c2f*atom.cpos;
          atom.name=atom.cleanname=tokens[4];
          if(LDEBUG){
            cerr << "XSTRUCTURE>>: AIMS atom[" << atom.name <<"] found (cartesian):" << endl;
            cerr << "    cpos" << atom.cpos << endl;
            cerr << "    fpos" << atom.fpos << endl;
          }
        }else if(tokens[0]=="atom_frac"){
          atom_found=true;atom.clear();
          atom.fpos[1]=aurostd::string2utype<double>(tokens[1]);
          atom.fpos[2]=aurostd::string2utype<double>(tokens[2]);
          atom.fpos[3]=aurostd::string2utype<double>(tokens[3]);
          atom.cpos=a.f2c*atom.fpos;
          atom.name=atom.cleanname=tokens[4];
          if(LDEBUG){
            cerr << "XSTRUCTURE>>: AIMS atom[" << atom.name <<"] found (fractional):" << endl;
            cerr << "    fpos" << atom.fpos << endl;
            cerr << "    cpos" << atom.cpos << endl;
          }
        }//ignore else, garbage
        if(atom_found){
          atom_props_search=true; //so we can look for other props on the next line
          
          //set some defaults! - START
          atom.name_is_given=true;
          //atom.CleanName();
          if(LDEBUG){
            cerr << "XSTRUCTURE>>: AIMS line=" << vinput[i] << endl;
            cerr << "XSTRUCTURE>>: AIMS atom.cleanname=" << atom.cleanname << endl;
          }
          atom.CleanSpin();
          //FIXED BELOW atom.number=atoms.size();    // reference position for convasp
          //FIXED BELOW atom.basis=atoms.size();     // position in the basis
          atom.ijk(1)=0;atom.ijk(2)=0;atom.ijk(3)=0; // inside the zero cell...
          atom.corigin(1)=0.0;atom.corigin(2)=0.0;atom.corigin(3)=0.0; // inside the zero cell
          atom.coord(1)=0.0;atom.coord(2)=0.0;atom.coord(3)=0.0; // inside the zero cell
          atom.spin=0.0;
          atom.noncoll_spin.clear(); //DX 12/5/17 - non-collinear spin
          // FIXED BELOW atom.type=itype;                // CONVASP_MODE if I want type 0,1,2,3,...
          atom.order_parameter_value=0;
          atom.order_parameter_atom=FALSE;
          atom.partial_occupation_value=1.0;
          atom.partial_occupation_flag=FALSE;
          //set some defaults! - STOP
          
          atoms.push_back(atom);
          //wait till later to add to structure, sort first that we don't screw up species
          //a.AddAtom(atom);
          //a.partial_occupation_sublattice.push_back(_pocc_no_sublattice_);  //default
        }
      }else if(atom_props_search && atoms.size()){
        //look and store additional properties here to atoms.back();
      }//ignore else
    }

		//sort first, then assign types
    //MUST BE STABLE SORT, absolutely critical
    //otherwise the order of equivalent names will be mixed EVERY TIME
    //this will screw up forces
		std::stable_sort(atoms.begin(),atoms.end(),sortAtomsNames);

		//assign types
		uint itype=0;
		atoms[0].type=itype;
		for(uint i=1;i<atoms.size();i++){
			if(atoms[i].name!=atoms[i-1].name){itype++;}
			atoms[i].type=itype;
		}

    a.atoms.clear();
    //now add atoms in right order (for species, etc.)
    for(uint i=0;i<atoms.size();i++){
      a.AddAtom(atoms[i]);  //does num_each_type and comp_each_type
      a.partial_occupation_sublattice.push_back(_pocc_no_sublattice_);
    }

    a.MakeBasis();

    if(LDEBUG){
      for(uint i=0;i<a.num_each_type.size();i++){
        cerr << "XSTRUCTURE>>: AIMS num_each_type[" << i <<"]=" <<a.num_each_type[i] << ", ";
        cerr << "comp_each_type[" << i <<"]=" <<a.comp_each_type[i] << endl;
      }
    }
    
		//grab title last
    a.title.clear();
    if(vinput.size()){
      std::size_t pos=vinput[0].find_first_of('#');
      if(pos!=std::string::npos){
        a.title=vinput[0];
        a.title.erase(pos,1);
        a.title=aurostd::RemoveWhiteSpacesFromTheFrontAndBack(a.title);
      }
    }
    if(a.title.empty()){a.buildGenericTitle();}
    
    if(LDEBUG){cerr << "XSTRUCTURE>>: AIMS title=" << a.title << endl;}

    a.partial_occupation_flag=FALSE;
    a.is_vasp4_poscar_format=FALSE;
    a.is_vasp5_poscar_format=FALSE;
    // DONE ?

  } // AIMS INPUT
  // ----------------------------------------------------------------------
  // COMMON CODE (NEED TO BE PATCHED, THOUGH).
  // FIX NORMAL AND PARTIAL OCCUPAITON
  if(LDEBUG) cerr << "XSTRUCTURE>>: COMMON CODE [9]" << endl;
  if(a.partial_occupation_flag==FALSE) {
    // have partial
    a.comp_each_type.clear();
    for(uint i=0;i<a.num_each_type.size();i++)
      a.comp_each_type.push_back((double) a.num_each_type.at(i));
  } else {
    // have partial
    if(poccaus.size()!=a.atoms.size()) {
      cerr << "ERROR XSTRUCTURE   poccaus.size()=" << poccaus.size() << " a.atoms.size()=" << a.atoms.size() << " " << endl;
      exit(0);
    }
    // create list (empty)
    a.comp_each_type.clear();
    for(uint i=0;i<a.num_each_type.size();i++)
      a.comp_each_type.push_back(0.0);
    // screen each one
    for(uint i=0;i<a.atoms.size();i++) {
      if(poccaus.at(i)<1.0-1e-5) {
	a.atoms.at(i).partial_occupation_flag=TRUE;
	a.atoms.at(i).partial_occupation_value=poccaus.at(i);
	a.comp_each_type.at(a.atoms.at(i).type)+=a.atoms.at(i).partial_occupation_value;
      } else {
	a.atoms.at(i).partial_occupation_flag=FALSE;
	a.atoms.at(i).partial_occupation_value=1.0;
	a.comp_each_type.at(a.atoms.at(i).type)+=a.atoms.at(i).partial_occupation_value;
      }
    }
  }
  a.GetStoich();  //170724 CO
  // ---------------------------------------------------------------
  // -------------- CREATE PARTIAL OCCUPATION STUFF
  // cerr << a.atoms.size() << endl;
  // cerr << a.partial_occupation_sublattice.size() << endl;
  // for(uint i=0;i<a.partial_occupation_sublattice.size();i++) cerr << a.partial_occupation_sublattice.at(i) << endl;
  if(LDEBUG) cerr << "XSTRUCTURE>>: PARTIAL OCCUPATION STUFF [10]" << endl;
  for(uint iatom=0;iatom<a.atoms.size();iatom++) {
    vector<uint> partial_occupation_sublattice_iatom;
    if(a.partial_occupation_sublattice.at(iatom)== _pocc_no_sublattice_) {
      partial_occupation_sublattice_iatom.push_back(iatom);
    } else {
      for(uint jatom=0;jatom<a.atoms.size();jatom++) {
	if(a.atoms.at(jatom).type==a.partial_occupation_sublattice.at(iatom))
	  partial_occupation_sublattice_iatom.push_back(iatom);	  
      }
    }
    if(0) {    cerr << a.partial_occupation_sublattice.at(iatom) << ", ";
      for(uint jatom=0;jatom<partial_occupation_sublattice_iatom.size();jatom++)
	cerr << partial_occupation_sublattice_iatom.at(jatom) << " ";
      cerr << endl;
    }
  }
  // ---------------------------------------------------------------
  // REMOVE ORDER PARAMETER
  if(LDEBUG) cerr << "XSTRUCTURE>>: REMOVE ORDER PARAMETER [11]" << endl;
  a.order_parameter_atoms.clear();
  for(uint i=0;i<a.atoms.size();i++)
    if(a.atoms.at(i).order_parameter_atom==TRUE)
      a.order_parameter_atoms.push_back(i);
  a.order_parameter_structure=(a.order_parameter_atoms.size()>0);
  // ---------------------------------------------------------------
  // -------------- SPECIES
  if(LDEBUG) cerr << "XSTRUCTURE>>: SPECIES [12]" << endl;
  if(LDEBUG) cerr << "XSTRUCTURE>>: a.is_vasp4_poscar_format=" << a.is_vasp4_poscar_format << endl;
  if(LDEBUG) cerr << "XSTRUCTURE>>: a.is_vasp5_poscar_format=" << a.is_vasp5_poscar_format << endl;
  if(a.is_vasp4_poscar_format==TRUE) {
    a.species.clear();a.species_pp.clear();a.species_pp_type.clear();a.species_pp_version.clear();a.species_pp_ZVAL.clear();a.species_pp_vLDAU.clear();a.species_volume.clear();a.species_mass.clear();
    // Plug the species
    for(uint i=0,j=0;i<a.num_each_type.size();j+=a.num_each_type.at(i),i++) {
      a.species.push_back(a.atoms.at(j).name);
      a.species_pp.push_back(a.atoms.at(j).name);
      a.species_pp_type.push_back("");
      a.species_pp_version.push_back("");
      a.species_pp_ZVAL.push_back(0.0);
      a.species_pp_vLDAU.push_back(deque<double>());
      a.species_volume.push_back(0.0);
      a.species_mass.push_back(0.0);
    }
  }
  if(a.is_vasp5_poscar_format==TRUE) {
    a.species.clear();a.species_pp.clear();a.species_pp_type.clear();a.species_pp_version.clear();a.species_pp_ZVAL.clear();a.species_pp_vLDAU.clear();a.species_volume.clear();a.species_mass.clear();
    if(a.iomode==IOVASP_POSCAR) aurostd::string2tokens((vinput.at(5)),tokens);
    if(a.iomode==IOVASP_ABCCAR) aurostd::string2tokens((vinput.at(3)),tokens);
    if(a.iomode==IOVASP_WYCKCAR) aurostd::string2tokens((vinput.at(3)),tokens);
    if(a.num_each_type.size()!=tokens.size()) {
      cerr << "ERRPR in aflow_xatom.cpp you need to specify the same number of species and atoms types" << endl;
      cerr << "      a.num_each_type.size()=" << a.num_each_type.size() << endl;
      cerr << "      tokens.size()=" << tokens.size() << endl;
      exit(0);
    }
    for(uint i=0;i<tokens.size();i++) {
      a.species.push_back(tokens.at(i));
      a.species_pp.push_back(tokens.at(i));
      a.species_pp_type.push_back("");
      a.species_pp_version.push_back("");
      a.species_pp_ZVAL.push_back(0.0);
      a.species_pp_vLDAU.push_back(deque<double>());
      a.species_volume.push_back(0.0);
      a.species_mass.push_back(0.0);
    }
    uint k=0;
    for(uint i=0;i<a.num_each_type.size();i++)
      for(uint j=0;j<(uint) a.num_each_type.at(i);j++) {
	if(a.atoms.at(k).name_is_given==FALSE) {
	  a.atoms.at(k).name=a.species.at(i);
	  a.atoms.at(k).name_is_given=TRUE;
	}
	k++;
      }
  }
  // ---------------------------------------------------------------
  // ALL atoms have been added. Now add the wyckoff ones
  if(LDEBUG) cerr << "XSTRUCTURE>>: WYCKCAR [13]" << endl;
  if(a.iomode==IOVASP_WYCKCAR) {
    // cerr << "[" << a.atoms.size() << "]" << endl;
    a=WyckoffPOSITIONS(a.spacegroupnumber,a.spacegroupnumberoption,a);
    a.title=a.title+" (WICKOFF "+a.spacegroup+" "+a.spacegrouplabel+")";
    // cerr << "[" << a.atoms.size() << "]" << endl;
  }

  // ---------------------------------------------------------------
  // Make spaces and links inside the qm part
  if(LDEBUG) cerr << "XSTRUCTURE>>: QM_SPACE [14]" << endl;
  for(uint i=0;i<a.atoms.size();i++) {
    xvector<double> v(3); v.clear();
    _atom atom;
    atom=a.atoms.at(i);
    atom.cpos.clear();atom.fpos.clear();
    a.qm_forces.push_back(v);
    a.qm_positions.push_back(v);
    a.qm_atoms.push_back(atom);
  }

  // ---------------------------------------------------------------
  if(LDEBUG) cerr << "XSTRUCTURE>>: WRAPPING [13]" << endl;
  // TOLERANCES ------------------------
  a.equiv_fpos_epsilon=_EQUIV_FPOS_EPS_; // standard but you can change
  // MAKE BASIS
  a.MakeBasis();
  // FLAGS -----------------------------
  a.Niggli_calculated=FALSE;
  a.Niggli_avoid=FALSE;
  a.Minkowski_calculated=FALSE;
  a.Minkowski_avoid=FALSE;
  a.LatticeReduction_calculated=FALSE;
  a.LatticeReduction_avoid=FALSE;
  // LATTICE FLAGS -----------------------------
  a.Standard_Lattice_calculated=FALSE;
  a.Standard_Lattice_avoid=FALSE;
  a.Standard_Lattice_primitive=FALSE;
  a.Standard_Lattice_conventional=FALSE;
  a.Standard_Lattice_has_failed=FALSE;
  a.bravais_lattice_type="";
  a.bravais_lattice_variation_type="";
  a.bravais_lattice_system="";
  a.bravais_lattice_lattice_type="";
  a.bravais_lattice_lattice_variation_type="";
  a.bravais_lattice_lattice_system="";
  a.pearson_symbol="";
  a.reciprocal_lattice_type="";
  a.reciprocal_lattice_variation_type="";
  a.bravais_superlattice_type="";
  a.bravais_superlattice_variation_type="";
  a.bravais_superlattice_system="";
  a.pearson_symbol_superlattice="";
  // QM CALCULATED STUFF
  a.qm_origin.clear();
  a.qm_scale=0.0;
  a.qm_lattice.clear();
  a.qm_klattice.clear();
  a.qm_f2c.clear();
  a.qm_c2f.clear();
  a.qm_calculated=FALSE;
  a.qm_forces_write=FALSE;
  a.qm_positions_write=FALSE;
  // CALCULATED STUFF
  //DX and CO - START
  a.dist_nn_min=AUROSTD_NAN;    //CO
  a.sym_eps=AUROSTD_NAN; //DX
  a.sym_eps_calculated=FALSE; //DX
  a.sym_eps_change_count=0; //DX 2/22/18 - added tolerance count specific to structure
  //DX and CO - END
  a.iatoms_calculated=FALSE;
  a.pgroup_calculated=FALSE;
  a.pgroup_xtal_calculated=FALSE;
  a.pgroupk_calculated=FALSE;
  a.pgroupk_xtal_calculated=FALSE; //DX 12/5/17 - Added pgroupk_xtal
  a.fgroup_calculated=FALSE;
  a.sgroup_calculated=FALSE;
  a.grid_atoms_calculated=FALSE;
  a.lijk_calculated=FALSE;
  a.neighbours_calculated=FALSE;
  // OUTPUT STUFF
  a.error_flag=FALSE;
  a.error_string="";
  a.write_lattice_flag=FALSE;
  a.write_inequivalent_flag=FALSE;
  a.write_klattice_flag=FALSE;
  a.write_DEBUG_flag=FALSE;
#ifdef XSTR_DEBUG
  a.write_lattice_flag=TRUE;
  a.write_inequivalent_flag=TRUE;
  a.write_klattice_flag=TRUE;
  a.write_DEBUG_flag=TRUE;
#endif
  // CHECKS
  if(a.atoms.size()!=a.qm_atoms.size())     {cerr << "ERROR operator>>xstructure [1] a.atoms.size()!=a.qm_atoms.size() " << endl;exit(0);}
  if(a.atoms.size()!=a.qm_forces.size())    {cerr << "ERROR operator>>xstructure [2] a.atoms.size()!=a.qm_forces.size() " << endl;exit(0);}
  if(a.atoms.size()!=a.qm_positions.size()) {cerr << "ERROR operator>>xstructure [3] a.atoms.size()!=a.qm_positions.size() " << endl;exit(0);}
  
  // ---------------------------------------------------------------
  if(LDEBUG) cerr << "XSTRUCTURE>>: DONE [99]" << endl;
  // DONE
  return cinput;
}

//DX 1/24/18 - added ibrav to lattice - START
// **************************************************************************
// pflow::QE_ibrav2lattice
// **************************************************************************
namespace pflow {
  xmatrix<double> QE_ibrav2lattice(const int& ibrav, const xvector<double>& parameters, const bool& isabc){
    // Lattice based on ibrav and lattice parameters (see http://www.quantum-espresso.org/wp-content/uploads/Doc/INPUT_PW.html#ibrav)
    xmatrix<double> lattice;

    double a,b,c,gamma,beta,alpha=0;

    if(isabc){
      a = parameters(1);               // parameters(1) = a
      b = parameters(2);               // parameters(2) = b
      c = parameters(3);               // parameters(3) = c
    }
    else{
      a = parameters(1)*bohr2angstrom;               // parameters(1) = a   //DX 2/15/18 - added bohr2angstrom (celldm is in Bohr)
      b = parameters(2)*parameters(1)*bohr2angstrom; // parameters(2) = b/a //DX 2/15/18 - added bohr2angstrom (celldm is in Bohr)
      c = parameters(3)*parameters(1)*bohr2angstrom; // parameters(3) = c/a //DX 2/15/18 - added bohr2angstrom (celldm is in Bohr)
    }
    gamma = acos(parameters(4));       // parameters(4) = cos(gamma) NOTE: Different than AFLOW's order convention of a,b,c,alpha,beta,gamma
    beta = acos(parameters(5));        // parameters(5) = cos(beta)  NOTE: Different than AFLOW's order convention of a,b,c,alpha,beta,gamma
    alpha = acos(parameters(6));       // parameters(6) = cos(alpha) NOTE: Different than AFLOW's order convention of a,b,c,alpha,beta,gamma

    xvector<double> xn(3);   xn(1)=1.0;xn(2)=0.0;xn(3)=0.0;
    xvector<double> yn(3);   yn(1)=0.0;yn(2)=1.0;yn(3)=0.0;
    xvector<double> zn(3);   zn(1)=0.0;zn(2)=0.0;zn(3)=1.0;
    xvector<double> a1(3),a2(3),a3(3);

    if(ibrav == 1){ //CUB
      a1 = a*xn; 
      a2 = a*yn; 
      a3 = a*zn; 
    }
    else if(ibrav == 2){ //FCC
      a1=-(1.0/2.0)*a*xn+(1.0/2.0)*a*zn;
      a2=(1.0/2.0)*a*yn+(1.0/2.0)*a*zn;
      a3=-(1.0/2.0)*a*xn+(1.0/2.0)*a*yn;
    }
    else if(ibrav == 3){ //BCC
      a1=(1.0/2.0)*a*xn+(1.0/2.0)*a*yn+(1.0/2.0)*a*zn;
      a2=-(1.0/2.0)*a*xn+(1.0/2.0)*a*yn+(1.0/2.0)*a*zn;
      a3=-(1.0/2.0)*a*xn-(1.0/2.0)*a*yn+(1.0/2.0)*a*zn;
    }
    else if(ibrav == -3){ //BCC (symmetric axis)
      a1=-(1.0/2.0)*a*xn+(1.0/2.0)*a*yn+(1.0/2.0)*a*zn;
      a2=(1.0/2.0)*a*xn-(1.0/2.0)*a*yn+(1.0/2.0)*a*zn;
      a3=(1.0/2.0)*a*xn+(1.0/2.0)*a*yn-(1.0/2.0)*a*zn;
    }
    else if(ibrav == 4){ //Hexagonal/Trigonal
      a1=a*xn;
      a2=-(1.0/2.0)*a*xn+(sqrt(3.0)/2.0)*a*yn;
      a3=c*zn;
    }
    else if(ibrav == 5){ //Trigonal R
      a1=sqrt((1.0-cos(gamma))/2.0)*a*xn-sqrt((1.0-cos(gamma))/6.0)*a*yn+sqrt((1.0+(2.0*cos(gamma)))/3.0)*a*zn;
      a2=2.0*sqrt((1.0-cos(gamma))/6.0)*a*yn+sqrt((1.0+(2.0*cos(gamma)))/3.0)*a*zn;
      a3=-sqrt((1.0-cos(gamma))/2.0)*a*xn-sqrt((1.0-cos(gamma))/6.0)*a*yn+sqrt((1.0+(2.0*cos(gamma)))/3.0)*a*zn;
    }
    else if(ibrav == -5){ //Trigonal R (3-fold axis)
      a1 = (a/sqrt(3.0))*((sqrt((1.0+(2.0*cos(gamma)))/3.0))-2.0*sqrt(2.0)*(sqrt((1.0-cos(gamma))/6.0)))*xn+(a/sqrt(3.0))*((sqrt((1.0+(2.0*cos(gamma)))/3.0))+sqrt(2.0)*(sqrt((1.0-cos(gamma))/6.0)))*yn+(a/sqrt(3.0))*((sqrt((1.0+(2.0*cos(gamma)))/3.0))+sqrt(2.0)*(sqrt((1.0-cos(gamma))/6.0)))*zn; 
      a2 = (a/sqrt(3.0))*((sqrt((1.0+(2.0*cos(gamma)))/3.0))+sqrt(2.0)*(sqrt((1.0-cos(gamma))/6.0)))*xn+(a/sqrt(3.0))*((sqrt((1.0+(2.0*cos(gamma)))/3.0))-2.0*sqrt(2.0)*(sqrt((1.0-cos(gamma))/6.0)))*yn+(a/sqrt(3.0))*((sqrt((1.0+(2.0*cos(gamma)))/3.0))+sqrt(2.0)*(sqrt((1.0-cos(gamma))/6.0)))*zn; 
      a3 = (a/sqrt(3.0))*((sqrt((1.0+(2.0*cos(gamma)))/3.0))+sqrt(2.0)*(sqrt((1.0-cos(gamma))/6.0)))*xn+(a/sqrt(3.0))*((sqrt((1.0+(2.0*cos(gamma)))/3.0))+sqrt(2.0)*(sqrt((1.0-cos(gamma))/6.0)))*yn+(a/sqrt(3.0))*((sqrt((1.0+(2.0*cos(gamma)))/3.0))-2.0*sqrt(2.0)*(sqrt((1.0-cos(gamma))/6.0)))*zn; 
    }
    else if(ibrav == 6){ //TET
      a1=a*xn;
      a2=a*yn;
      a3=c*zn;
    }
    else if(ibrav == 7){ //BCT
      a1=(1.0/2.0)*a*xn-(1.0/2.0)*a*yn+(1.0/2.0)*c*zn;
      a2=(1.0/2.0)*a*xn+(1.0/2.0)*a*yn+(1.0/2.0)*c*zn;
      a3=-(1.0/2.0)*a*xn-(1.0/2.0)*a*yn+(1.0/2.0)*c*zn;
    }
    else if(ibrav == 8){ //ORC
      a1=a*xn;
      a2=b*yn;
      a3=c*zn;
    }
    else if(ibrav == 9){ //ORCC
      a1=(1.0/2.0)*a*xn+(1.0/2.0)*b*yn;
      a2=-(1.0/2.0)*a*xn+(1.0/2.0)*b*yn;
      a3=c*zn;
    }
    else if(ibrav == -9){ //ORCC - 2
      a1=(1.0/2.0)*a*xn-(1.0/2.0)*b*yn;
      a2=(1.0/2.0)*a*xn+(1.0/2.0)*b*yn;
      a3=c*zn;
    }
    else if(ibrav == 10){ //ORCF
      a1=(1.0/2.0)*a*xn+(1.0/2.0)*c*zn;
      a2=(1.0/2.0)*a*xn+(1.0/2.0)*b*yn;
      a3=(1.0/2.0)*b*yn+(1.0/2.0)*c*zn;
    }
    else if(ibrav == 11){ //ORCI
      a1=(1.0/2.0)*a*xn+(1.0/2.0)*b*yn+(1.0/2.0)*c*zn;
      a2=-(1.0/2.0)*a*xn+(1.0/2.0)*b*yn+(1.0/2.0)*c*zn;
      a3=-(1.0/2.0)*a*xn-(1.0/2.0)*b*yn+(1.0/2.0)*c*zn;
    }
    else if(ibrav == 12){ //MCL (unique axis c)
      a1=a*xn;
      a2=cos(gamma)*b*xn+sin(gamma)*b*yn;
      a3=c*zn;
    }
    else if(ibrav == -12){ //MCL (unique axis b)
      a1=a*xn;
      a2=b*yn;
      a3=cos(beta)*c*xn+sin(beta)*c*zn;
    }
    else if(ibrav == 13){ //MCLC
      a1=(1.0/2.0)*a*xn-(1.0/2.0)*c*zn;
      a2=cos(gamma)*b*xn+sin(gamma)*b*yn;
      a3=(1.0/2.0)*a*xn+(1.0/2.0)*c*zn;
    }
    else if(ibrav == 14){ //TRI
      double cx=c*cos(beta);
      double cy=c*(cos(alpha)-cos(beta)*cos(gamma))/sin(gamma);
      double cz=sqrt(pow(c,2.0)-pow(cx,2.0)-pow(cy,2.0));

      a1=a*xn;
      a2=b*cos(gamma)*xn+b*sin(gamma)*yn;
      a3=cx*xn+cy*yn+cz*zn;
    } 

    lattice(1,1)=a1(1);lattice(1,2)=a1(2);lattice(1,3)=a1(3);
    lattice(2,1)=a2(1);lattice(2,2)=a2(2);lattice(2,3)=a2(3);
    lattice(3,1)=a3(1);lattice(3,2)=a3(2);lattice(3,3)=a3(3);

    return lattice;
  }
}
//DX 1/24/18 - added ibrav to lattice - START

// **************************************************************************
// xstructure::GetStoich
// **************************************************************************
// Get stoichiometries
bool xstructure::GetStoich(void) { //CO 171025
  double total_comp=0.0;
  for(uint i=0;i<comp_each_type.size();i++){
    total_comp+=comp_each_type[i];
  }
  stoich_each_type.clear();
  for(uint i=0;i<comp_each_type.size();i++){
    stoich_each_type.push_back( comp_each_type[i]/total_comp );
  }
  return TRUE;
}

// **************************************************************************
// xstructure::FixLattices
// **************************************************************************
// Fix all the lattices (you can do by hand, but here everything is done)
bool xstructure::FixLattices(void) {
  klattice=ReciprocalLattice(lattice,scale);
  f2c=trasp(lattice);
  c2f=inverse(trasp(lattice));
  // lattice is the reference... everything else depends
  //  if(iomode==IOVASP_POSCAR || iomode==IOVASP_AUTO) {
  xvector<double> data(6);
  data=Getabc_angles(lattice,DEGREES);
  a=data[1];b=data[2];c=data[3];
  alpha=data[4];beta=data[5];gamma=data[6];
  //  }
  //  if(iomode==IOVASP_ABCCAR || iomode==IOVASP_WYCKCAR) {
  //  lattice=GetClat(a,b,c,alpha,beta,gamma);
  // }
  return TRUE;
}


// **************************************************************************
// GetStructure
// **************************************************************************
// get a structure from a directory (POSCAR for VASP
xstructure GetStructure(const int& iomode,ifstream& input) {
  xstructure out;
  if(iomode==IOVASP_POSCAR) { // VASP POSCAR
    out.iomode=IOVASP_POSCAR;
    if(!input) {
      cerr << "EEEEE   File not found" << endl;
      input.clear();input.close();
      out.error_flag=TRUE;
      out.error_string="FILE NOT FOUND";
      return out;
    }
    input >> out;
  }
  return out;
}

// **************************************************************************
// GetStructure
// **************************************************************************
xstructure GetStructure(const int& iomode,const string& Directory) {
  xstructure out;
  if(iomode==IOVASP_POSCAR) { // VASP POSCAR
    out.iomode=IOVASP_POSCAR;
    ifstream input;
    string File;
    File=Directory+"/POSCAR";
    input.open(File.c_str(),std::ios::in);
    out=GetStructure(iomode,input);
    input.clear();input.close();
    return out;
  }
  return out;
}

// **************************************************************************
// xstructure::SetCoordinates
// **************************************************************************
// change coordinates type
void xstructure::SetCoordinates(const int& mode)  {
  switch(mode) {
  case _UPDATE_LATTICE_VECTORS_TO_ABCANGLES_ : {
    cerr << "ERROR SetCoordinate [1] mode=" << mode << endl; exit(0);
    break;
  }
  case _UPDATE_LATTICE_ABCANGLES_TO_VECTORS_ : {
    cerr << "ERROR SetCoordinate [2] mode=" << mode << endl; exit(0);
    break;
  }
  case _COORDS_CARTESIAN_ : {
    coord_flag=_COORDS_CARTESIAN_;
    strcpy(coord_type,"C");
    break;
  }
  case _COORDS_FRACTIONAL_ : {
    coord_flag=_COORDS_FRACTIONAL_;
    strcpy(coord_type,"D");
    break;
  }
  default: {
    cerr << "SetCoordinates NOTHING TO DO  mode=" << mode << endl;
  }
  }
}


// **************************************************************************
// xstructure::MakeBasis
// **************************************************************************
// This adds an atom to the structure.
void xstructure::MakeBasis(void) {
  // need to update NUMBER and BASIS
  for(uint iatom=0;iatom<atoms.size();iatom++) {
    atoms.at(iatom).basis=iatom;
    atoms.at(iatom).number=iatom;
  }
}

// **************************************************************************
// xstructure::AddAtom
// **************************************************************************
// This adds an atom to the structure.

void xstructure::AddAtom(const _atom& atom) {
  bool LDEBUG=(FALSE || XHOST.DEBUG); 
  _atom btom;btom=atom;

  // check that this atom is not already present
  xvector<double> a1(3),a2(3),a3(3),aijk(3);
  a1=lattice(1);a2=lattice(2);a3=lattice(3);

  bool FOUND_POSITION=FALSE;
  for(uint iat=0;iat<atoms.size()&&FOUND_POSITION==FALSE;iat++)
    if(atoms.at(iat).type==atom.type && atoms.at(iat).name==atom.name)
      for(int i=-1;i<=1&&FOUND_POSITION==FALSE;i++)
	for(int j=-1;j<=1&&FOUND_POSITION==FALSE;j++)
	  for(int k=-1;k<=1&&FOUND_POSITION==FALSE;k++) {
	    aijk[1]=i;aijk[2]=j;aijk[3]=k;
	    //	if(aurostd::modulus(atoms.at(iat).cpos-(((double)i)*a1+((double)j)*a2+((double)k)*a3+atom.cpos))<0.1) FOUND_POSITION=TRUE;
            //DX and CO -START    
            //DX if(aurostd::modulus(atoms.at(iat).fpos-(aijk+atom.fpos))<0.01) FOUND_POSITION=TRUE;
            if((*this).sym_eps!=AUROSTD_NAN && (*this).sym_eps<AUROSTD_NAN && (*this).sym_eps>1e-10){ //DX 12/1/17 - Added (*this).sym_eps>1e-10 //DX 2/15/18 - added (*this).sym_eps<AUROSTD_NAN (needed) 
              if(aurostd::modulus((*this).f2c*(atoms.at(iat).fpos-(aijk+atom.fpos)))<=(*this).sym_eps) FOUND_POSITION=TRUE; //DX
            }
            else{ 
            //if(aurostd::modulus(atoms.at(iat).fpos-(aijk+atom.fpos))<1e-10) FOUND_POSITION=TRUE; //DX
	    if(aurostd::modulus(atoms.at(iat).cpos-(((double)i)*a1+((double)j)*a2+((double)k)*a3+atom.cpos))<0.1) FOUND_POSITION=TRUE; //DX 12/1/17
            }
            //DX and CO - END
	  }
  if(FOUND_POSITION==TRUE) return; // found no need to add it further

  // now found that it does not exist check type
  //  cerr << "AddAtom new atom" << endl;
  bool FOUND_SPECIES=FALSE;
  uint species_position=0;
  for(uint isp=0;isp<species.size()&&FOUND_SPECIES==FALSE;isp++)
    if(atom.name==species.at(isp)) {FOUND_SPECIES=TRUE;species_position=isp;}
  
  if(FOUND_SPECIES==FALSE) {
    if(LDEBUG) cerr << "AddAtom new_specie=" << atom.name << endl;
    num_each_type.push_back(1);
    comp_each_type.push_back(atom.partial_occupation_value);
    species.push_back(atom.name); // cerr << "AddAtom=" << atom.name << endl;
    species_pp.push_back(atom.name); // cerr << "AddAtom=" << atom.name << endl;
    species_pp_type.push_back(""); // cerr << "AddAtom=" << atom.name << endl;
    species_pp_version.push_back(""); // cerr << "AddAtom=" << atom.name << endl;
    species_pp_ZVAL.push_back(0.0); // cerr << "AddAtom=" << atom.name << endl;
    species_pp_vLDAU.push_back(deque<double>()); // cerr << "AddAtom=" << atom.name << endl;
    species_volume.push_back(GetAtomVolume(atom.name)); // cerr << "AddAtom=" << atom.name << endl;
    species_mass.push_back(GetAtomMass(atom.name)); // cerr << "AddAtom=" << atom.name << endl;
  } else {
    // cerr <<  num_each_type.size() << " " <<  btom.type << endl;
    // cerr <<  comp_each_type.size() << " " <<  btom.type << endl;
    if(LDEBUG) cerr << "AddAtom increasing species_position " << species_position << endl;
    num_each_type.at(species_position)++;
    comp_each_type.at(species_position)+=atom.partial_occupation_value;
  }
  if(btom.name_is_given) {
    btom.CleanName();
    //DX 9/21/17 - Need to keep spin info  btom.CleanSpin();
  }
  GetStoich();  //170724 - CO

  // OLD STYLE
  //  atoms.push_back(btom);  MakeBasis(); return;

  // NEW STYLE
  bool found=FALSE;
  if(0)  for(uint iat=0;iat<atoms.size()&&!found;iat++) {
      if(iat<atoms.size()-1) {
        if(atoms.at(iat).type==btom.type && atoms.at(iat+1).type!=btom.type) {
          //	if(LDEBUG)
          cerr << "HERE1 iat=" << iat << "  atoms.at(iat).type=" << atoms.at(iat).type << "  btom.type=" << btom.type << endl;//" atoms.begin()=" <<  long(atoms.begin()) << endl;
          atoms.insert(iat+atoms.begin()+1,btom); // potential problem  with CAST
          found=TRUE;
        }
      }
    }
  if(1) {
    std::deque<_atom>::iterator it=atoms.begin();
    for(uint iat=0;iat<atoms.size()&&!found;iat++,it++) {
      if(iat<atoms.size()-1) {
	//	cerr << "HERE0 iat=" << iat << "  atoms.at(iat).type=" << atoms.at(iat).type << "  btom.type=" << btom.type << endl;
	if(atoms.at(iat).type==btom.type && atoms.at(iat+1).type!=btom.type) {
	  //	if(LDEBUG)
	  //	  cerr << "HERE1 iat=" << iat << "  atoms.at(iat).type=" << atoms.at(iat).type << "  btom.type=" << btom.type << endl;//" atoms.begin()=" <<  long(atoms.begin()) << endl;
	  atoms.insert(it+1,btom);  // it is iterator, fine for insert.
	  found=TRUE;
	}
      }
    }
  }
  // if never found add at the end
  if(!found) atoms.push_back(btom);
  // 
  //  atoms.push_back(btom);
  // done
  MakeBasis(); // need to update NUMBER and BASIS
}

// **************************************************************************
// xstructure::RemoveAtom
// **************************************************************************
// This removes an atom to the structure.
// 170721 - CO added some safety checks to make sure we weren't deleting an entry
// of a vector/deque that didn't exist (see if size() > itype)
// this is really important because if we build an xstructure on the fly and only occupy
// some of these attributes, the code will seg fault badly, and the error will not
// show up until the xstructure is deconstructed (super confusing)
// make sure to code safely always!
void xstructure::RemoveAtom(const uint& iatom) {
  if(iatom<atoms.size()) {
    uint itype=atoms.at(iatom).type;
    if(num_each_type.size()>itype) num_each_type.at(itype)--;
    if(comp_each_type.size()>itype) comp_each_type.at(itype)-=atoms.at(iatom).partial_occupation_value;
    if(num_each_type.size()>itype && num_each_type.at(itype)==0) {
      if(num_each_type.size()>itype)      num_each_type.erase(num_each_type.begin()+itype);    // erase num_each_type
      if(comp_each_type.size()>itype)     comp_each_type.erase(comp_each_type.begin()+itype);  // erase comp_each_type
      if(species.size()>itype)            species.erase(species.begin()+itype);  // erase species
      if(species_pp.size()>itype)         species_pp.erase(species_pp.begin()+itype);  // erase species_pp
      if(species_pp_type.size()>itype)    species_pp_type.erase(species_pp_type.begin()+itype);  // erase species_pp_type
      if(species_pp_version.size()>itype) species_pp_version.erase(species_pp_version.begin()+itype);  // erase species_pp_version
      if(species_pp_ZVAL.size()>itype)    species_pp_ZVAL.erase(species_pp_ZVAL.begin()+itype);  // erase species_pp_ZVAL
      if(species_pp_vLDAU.size()>itype)   species_pp_vLDAU.erase(species_pp_vLDAU.begin()+itype);  // erase species_pp_vLDAU
      if(species_volume.size()>itype)     species_volume.erase(species_volume.begin()+itype);  // erase species_volume
      if(species_mass.size()>itype)       species_mass.erase(species_mass.begin()+itype);  // erase species_mass
      //CO 170721 - might be better if we wrote function like MakeBasis() for types and did 
      //at the end, but it is not unsafe (seg fault) in the way it is written
      for(uint i=0;i<atoms.size();i++)
        if(i!=iatom && atoms.at(i).type>(int)itype)
	  atoms.at(i).type--;
    }
    //170721 CO - this is obsolete with MakeBasis() below!
    // do atoms
    //for(uint i=0;i<atoms.size();i++) {
    //  if(i!=iatom && atoms.at(i).number>atoms.at(iatom).number)
    //    atoms.at(i).number--;
    //  if(i!=iatom && atoms.at(i).basis>atoms.at(iatom).basis)
    //    atoms.at(i).basis--;
    //}
    atoms.erase(atoms.begin()+iatom);
    //    // do partial_occupation_flags
    //    for(uint i=0;i<partial_occupation_flags.size();i++)
    //      if(iatom==partial_occupation_flags.at(i))
    //	partial_occupation_flags.erase(partial_occupation_flags.begin()+i);
    // do order_parameter_atoms
    // 170721 CO - this is okay as it won't seg fault (doesn't delete anything bigger than .size() )
    for(uint i=0;i<order_parameter_atoms.size();i++)
      if(iatom==order_parameter_atoms.at(i))
	order_parameter_atoms.erase(order_parameter_atoms.begin()+i);
  }
  GetStoich();  //170724 CO
  // done
  MakeBasis(); // need to update NUMBER and BASIS
}

void xstructure::RemoveAtoms(vector<uint>& v_atoms_to_remove) {
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  std::sort(v_atoms_to_remove.rbegin(),v_atoms_to_remove.rend()); //NOTE the r, reverse sort, that way when we remove, it doesn't affect other indices
  for(uint atom=0;atom<v_atoms_to_remove.size();atom++){
    if(LDEBUG){cerr << "Removing Atom " <<  v_atoms_to_remove[atom] << endl;}
    RemoveAtom(v_atoms_to_remove[atom]);
  }
}

// **************************************************************************
// xstructure::RemoveCopies xstructure::RemoveFractionalCopies xstructure::RemoveCartesianCopies
// **************************************************************************
// This removes atoms too close than tol
void xstructure::RemoveCopies(double tol) {
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  bool flag_isprimitive=FALSE,flag_remove=FALSE;
  uint iat1=0,iat2=0,irm=0;
  while(flag_isprimitive==FALSE) {
    if(iat1==atoms.size()) {
      flag_isprimitive=TRUE;
    } else {
      flag_remove=false;
      for(iat2=iat1+1;iat2<atoms.size()&&flag_remove==FALSE;iat2++)
	if((aurostd::modulus(atoms.at(iat1).fpos-atoms.at(iat2).fpos)<tol || aurostd::modulus(atoms.at(iat1).cpos-atoms.at(iat2).cpos)<tol) && flag_remove==FALSE) {
	  flag_remove=TRUE;
	  irm=iat2;
	}
      if(flag_remove==TRUE) RemoveAtom(irm);
      else iat1++;
    }
    if(LDEBUG) cout << "DEBUG (RemoveCopies) iat1=" << iat1 << endl;
  }
}
void xstructure::RemoveFractionalCopies(double tol) {
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  bool flag_isprimitive=FALSE,flag_remove=FALSE;
  uint iat1=0,iat2=0,irm=0;
  while(flag_isprimitive==FALSE) {
    if(iat1==atoms.size()) {
      flag_isprimitive=TRUE;
    } else {
      flag_remove=false;
      for(iat2=iat1+1;iat2<atoms.size()&&flag_remove==FALSE;iat2++)
	if(aurostd::modulus(atoms.at(iat1).fpos-atoms.at(iat2).fpos)<tol && flag_remove==FALSE) {
	  flag_remove=TRUE;
	  irm=iat2;
	}
      if(flag_remove==TRUE) RemoveAtom(irm);
      else iat1++;
    }
    if(LDEBUG) cout << "DEBUG (RemoveFractionalCopies) iat1=" << iat1 << endl;
  }
}
void xstructure::RemoveCartesianCopies(double tol) {
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  bool flag_isprimitive=FALSE,flag_remove=FALSE;
  uint iat1=0,iat2=0,irm=0;
  while(flag_isprimitive==FALSE) {
    if(iat1==atoms.size()) {
      flag_isprimitive=TRUE;
    } else {
      flag_remove=false;
      for(iat2=iat1+1;iat2<atoms.size()&&flag_remove==FALSE;iat2++)
	if(aurostd::modulus(atoms.at(iat1).cpos-atoms.at(iat2).cpos)<tol && flag_remove==FALSE) {
	  flag_remove=TRUE;
	  irm=iat2;
	}
      if(flag_remove==TRUE) RemoveAtom(irm);
      else iat1++;
    }
    if(LDEBUG) cout << "DEBUG (RemoveCartesianCopies) iat1=" << iat1 << endl;
  }
}

// **************************************************************************
// xstructure::AddCorners
// **************************************************************************
void xstructure::AddCorners(void) {
  xstructure str;
  BringInCell();str=*this;
  while(atoms.size()) RemoveAtom(0);   
  for(uint iat=0;iat<str.atoms.size();iat++) {
    for(double i=0;i<=1;i+=0.99) {
      for(double j=0;j<=1;j+=0.99) {
	for(double k=0;k<=1;k+=0.99) {
	  _atom atom=str.atoms.at(iat);
	  atom.fpos[1]+=i;atom.fpos[2]+=j;atom.fpos[3]+=k;
	  atom.cpos=F2C(lattice,atom.fpos);
	  if(atom.fpos[1]<=1.0 && atom.fpos[2]<=1.0 && atom.fpos[3]<=1.0) AddAtom(atom);
	  //	    if(aurostd::isequal(atom.fpos[1],1.0,0.02) && atom.fpos[2]<1.0 && atom.fpos[3]<1.0)   AddAtom(atom);	    //AddAtom(atom);
	}
      }
    }
  }
  title=title+" with_corners";
}

// **************************************************************************
// xstructure::ShifOriginToAtom
// **************************************************************************
// // Shift the origin to atom(iat)
void xstructure::ShifOriginToAtom(const int& iat) {
  //DX and CO - START
  if(iat<0 || iat>=(int)atoms.size()) {
    cerr << "ERROR void xstructure::ShifOriginToAtom(const int& iat),  iat=" << iat << " out of boundaries (0," << atoms.size()-1 << ")" << endl;
    exit(0);
  }
  xvector<double> frigin(3);
  origin=atoms[iat].cpos;
  frigin=atoms[iat].fpos;
  for(uint i=0;i<atoms.size();i++) {
    atoms[i].fpos=atoms[i].fpos-frigin;
    atoms[i].cpos=atoms[i].cpos-origin;
  }
  //CO TOO SLOW
  //origin=atoms.at(iat).cpos;
  //frigin=atoms.at(iat).fpos;
  //for(uint i=0;i<atoms.size();i++) {
  //  atoms.at(i).fpos=atoms.at(i).fpos-frigin;
  //  atoms.at(i).cpos=atoms.at(i).cpos-origin;
  //}
  //DX and CO - END
}

// **************************************************************************
// SetSDNumbers
// **************************************************************************
xstructure SetSDNumbers(const xstructure& a, const vector<string>& in_sd) {
  // Note that in_sd has one name for each type, not one name for each atom.
  xstructure b(a);
  //  int size=in_sd.size();
  for(uint cnt=0;cnt<b.atoms.size();cnt++) {
    b.atoms[cnt].sd=in_sd[cnt];
    if(in_sd[cnt].size()<2) {
      cerr << "WARNING:  Atom=" << cnt << " you must specify SD strings 3 characters long (switching to TTT)" << endl;
      b.atoms.at(cnt).sd = "TTT";
    }
  }
  return b;
}


// **************************************************************************
// SetSDTypes
// **************************************************************************
xstructure SetSDTypes(const xstructure& a, const vector<string>& in_sd) {
  // Note that in_sd has one name for each type, not one name for each atom.
  xstructure b(a);
  int cnt=-1;
  int size=in_sd.size();
  for(int i=0;i<(int)size;i++) {
    if(i<(int)b.num_each_type.size()) {
      for(int j=0;j<(int)b.num_each_type.at(i);j++) {
        cnt++;
        b.atoms.at(cnt).sd=in_sd[i];
	if(in_sd[i].size()<2) {
	  cerr << "WARNING:  Atom=" << cnt << " you must specify SD strings 3 characters long (switching to TTT)" << endl;
	  b.atoms.at(cnt).sd="TTT";
	}
      }
    }
  }
  return b;
}

// **************************************************************************
// GetTypes
// **************************************************************************
vector<int> GetTypes(const xstructure& a) {
  vector<int> out_type;
  for(int i=0;i<(int)a.atoms.size();i++)
    out_type.push_back(a.atoms.at(i).type);
  return out_type;
}

// **************************************************************************
// GetNames
// **************************************************************************
vector<string> GetNames(const xstructure& a) {
  vector<string> out_name;
  for(int i=0;i<(int)a.atoms.size();i++)
    out_name.push_back(a.atoms.at(i).name);
  return out_name;
}

// **************************************************************************
// GetCleanNames
// **************************************************************************
vector<string> GetCleanNames(const xstructure& a) {
  vector<string> out_cleanname;
  for(int i=0;i<(int)a.atoms.size();i++)
    out_cleanname.push_back(a.atoms.at(i).cleanname);
  return out_cleanname;
}

// **************************************************************************
// GetSpins
// **************************************************************************
vector<double> GetSpins(const xstructure& a) {
  vector<double> out_spin;
  for(int i=0;i<(int)a.atoms.size();i++)
    out_spin.push_back(a.atoms.at(i).spin);
  return out_spin;
}

// ***************************************************************************
// GetElementName
// ***************************************************************************
string GetElementName(string stringin) {
  // need to clean up the _pv stuff of VASP
  for(uint i=0;i<NUM_ELEMENTS;i++)
    if(stringin==atom_symbol_vec.at(i))
      return atom_name_vec.at(i);
  return "NotFound";
}

// ***************************************************************************
// GetSpaceGroupName
// ***************************************************************************
string GetSpaceGroupName(int spacegroupnumber) {
  string spacegroup;
  if(spacegroupnumber < 1 || spacegroupnumber > 230) {
    cerr << "aflow_xatom.cpp GetSpaceGroupName() routine: " << endl;
    cerr << "  space group specified invalid (1-230): " << spacegroupnumber << endl;
    exit(0);
  }
  // OK
  if(spacegroupnumber==1) { // ------------------- 1  P1 #1
    spacegroup="P1";}
  if(spacegroupnumber==2) { // ------------------- 2  P-1 #2
    spacegroup="P-1";}
  if(spacegroupnumber==3) { // ------------------- 3  P2 #3
    spacegroup="P2";}
  if(spacegroupnumber==4) { // ------------------- 4  P2_{1} #4
    spacegroup="P2_{1}";}
  if(spacegroupnumber==5) { // ------------------- 5  C2 #5
    spacegroup="C2";}
  if(spacegroupnumber==6) { // ------------------- 6  Pm #6
    spacegroup="Pm";}
  if(spacegroupnumber==7) { // ------------------- 7  Pc #7
    spacegroup="Pc";}
  if(spacegroupnumber==8) { // ------------------- 8  Cm #8
    spacegroup="Cm";}
  if(spacegroupnumber==9) { // ------------------- 9  Cc #9
    spacegroup="Cc";}
  if(spacegroupnumber==10) { // ------------------- 10  P2/m #10
    spacegroup="P2/m";}
  if(spacegroupnumber==11) { // ------------------- 11  P2_{1}/m #11
    spacegroup="P2_{1}/m";}
  if(spacegroupnumber==12) { // ------------------- 12  C2/m #12
    spacegroup="C2/m";}
  if(spacegroupnumber==13) { // ------------------- 13  P2/c #13
    spacegroup="P2/c";}
  if(spacegroupnumber==14) { // ------------------- 14  P2_{1}/c #14
    spacegroup="P2_{1}/c";}
  if(spacegroupnumber==15) { // ------------------- 15  C2/c #15
    spacegroup="C2/c";}
  if(spacegroupnumber==16) { // ------------------- 16  P222 #16
    spacegroup="P222";}
  if(spacegroupnumber==17) { // ------------------- 17  P222_{1} #17
    spacegroup="P222_{1}";}
  if(spacegroupnumber==18) { // ------------------- 18  P2_{1}2_{1}2 #18
    spacegroup="P2_{1}2_{1}2";}
  if(spacegroupnumber==19) { // ------------------- 19  P2_{1}2_{1}2_{1} #19
    spacegroup="P2_{1}2_{1}2_{1}";}
  if(spacegroupnumber==20) { // ------------------- 20  C222_{1} #20
    spacegroup="C222_{1}";}
  if(spacegroupnumber==21) { // ------------------- 21  C222 #21
    spacegroup="C222";}
  if(spacegroupnumber==22) { // ------------------- 22  F222 #22
    spacegroup="F222";}
  if(spacegroupnumber==23) { // ------------------- 23  I222 #23
    spacegroup="I222";}
  if(spacegroupnumber==24) { // ------------------- 24  I2_{1}2_{1}2_{1} #24
    spacegroup="I2_{1}2_{1}2_{1}";}
  if(spacegroupnumber==25) { // ------------------- 25  Pmm2 #25
    spacegroup="Pmm2";}
  if(spacegroupnumber==26) { // ------------------- 26  Pmc2_{1} #26
    spacegroup="Pmc2_{1}";}
  if(spacegroupnumber==27) { // ------------------- 27  Pcc2 #27
    spacegroup="Pcc2";}
  if(spacegroupnumber==28) { // ------------------- 28  Pma2 #28
    spacegroup="Pma2";}
  if(spacegroupnumber==29) { // ------------------- 29  Pca2_{1} #29
    spacegroup="Pca2_{1}";}
  if(spacegroupnumber==30) { // ------------------- 30  Pnc2 #30
    spacegroup="Pnc2";}
  if(spacegroupnumber==31) { // ------------------- 31  Pmn2_{1} #31
    spacegroup="Pmn2_{1}";}
  if(spacegroupnumber==32) { // ------------------- 32  Pba2 #32
    spacegroup="Pba2";}
  if(spacegroupnumber==33) { // ------------------- 33  Pna2_{1} #33
    spacegroup="Pna2_{1}";}
  if(spacegroupnumber==34) { // ------------------- 34  Pnn2 #34
    spacegroup="Pnn2";}
  if(spacegroupnumber==35) { // ------------------- 35  Cmm2 #35
    spacegroup="Cmm2";}
  if(spacegroupnumber==36) { // ------------------- 36  Cmc2_{1} #36
    spacegroup="Cmc2_{1}";}
  if(spacegroupnumber==37) { // ------------------- 37  Ccc2 #37
    spacegroup="Ccc2";}
  if(spacegroupnumber==38) { // ------------------- 38  Amm2 #38
    spacegroup="Amm2";}
  if(spacegroupnumber==39) { // ------------------- 39  Aem2 #39
    spacegroup="Aem2";}
  if(spacegroupnumber==40) { // ------------------- 40  Ama2 #40
    spacegroup="Ama2";}
  if(spacegroupnumber==41) { // ------------------- 41  Aea2 #41
    spacegroup="Aea2";}
  if(spacegroupnumber==42) { // ------------------- 42  Fmm2 #42
    spacegroup="Fmm2";}
  if(spacegroupnumber==43) { // ------------------- 43  Fdd2 #43
    spacegroup="Fdd2";}
  if(spacegroupnumber==44) { // ------------------- 44  Imm2 #44
    spacegroup="Imm2";}
  if(spacegroupnumber==45) { // ------------------- 45  Iba2 #45
    spacegroup="Iba2";}
  if(spacegroupnumber==46) { // ------------------- 46  Ima2 #46
    spacegroup="Ima2";}
  if(spacegroupnumber==47) { // ------------------- 47  Pmmm #47
    spacegroup="Pmmm";}
  if(spacegroupnumber==48) { // ------------------- 48  Pnnn #48
    spacegroup="Pnnn";}
  if(spacegroupnumber==49) { // ------------------- 49  Pccm #49
    spacegroup="Pccm";}
  if(spacegroupnumber==50) { // ------------------- 50  Pban #50
    spacegroup="Pban";}
  if(spacegroupnumber==51) { // ------------------- 51  Pmma #51
    spacegroup="Pmma";}
  if(spacegroupnumber==52) { // ------------------- 52  Pnna #52
    spacegroup="Pnna";}
  if(spacegroupnumber==53) { // ------------------- 53  Pmna #53
    spacegroup="Pmna";}
  if(spacegroupnumber==54) { // ------------------- 54  Pcca #54
    spacegroup="Pcca";}
  if(spacegroupnumber==55) { // ------------------- 55  Pbam #55
    spacegroup="Pbam";}
  if(spacegroupnumber==56) { // ------------------- 56  Pccn #56
    spacegroup="Pccn";}
  if(spacegroupnumber==57) { // ------------------- 57  Pbcm #57
    spacegroup="Pbcm";}
  if(spacegroupnumber==58) { // ------------------- 58  Pnnm #58
    spacegroup="Pnnm";}
  if(spacegroupnumber==59) { // ------------------- 59  Pmmn #59
    spacegroup="Pmmn";}
  if(spacegroupnumber==60) { // ------------------- 60  Pbcn #60
    spacegroup="Pbcn";}
  if(spacegroupnumber==61) { // ------------------- 61  Pbca #61
    spacegroup="Pbca";}
  if(spacegroupnumber==62) { // ------------------- 62  Pnma #62
    spacegroup="Pnma";}
  if(spacegroupnumber==63) { // ------------------- 63  Cmcm #63
    spacegroup="Cmcm";}
  if(spacegroupnumber==64) { // ------------------- 64  Cmce #64
    spacegroup="Cmce";}
  if(spacegroupnumber==65) { // ------------------- 65  Cmmm #65
    spacegroup="Cmmm";}
  if(spacegroupnumber==66) { // ------------------- 66  Cccm #66
    spacegroup="Cccm";}
  if(spacegroupnumber==67) { // ------------------- 67  Cmme #67
    spacegroup="Cmme";}
  if(spacegroupnumber==68) { // ------------------- 68  Ccce #68
    spacegroup="Ccce";}
  if(spacegroupnumber==69) { // ------------------- 69  Fmmm #69
    spacegroup="Fmmm";}
  if(spacegroupnumber==70) { // ------------------- 70  Fddd #70
    spacegroup="Fddd";}
  if(spacegroupnumber==71) { // ------------------- 71  Immm #71
    spacegroup="Immm";}
  if(spacegroupnumber==72) { // ------------------- 72  Ibam #72
    spacegroup="Ibam";}
  if(spacegroupnumber==73) { // ------------------- 73  Ibca #73
    spacegroup="Ibca";}
  if(spacegroupnumber==74) { // ------------------- 74  Imma #74
    spacegroup="Imma";}
  if(spacegroupnumber==75) { // ------------------- 75  P4 #75
    spacegroup="P4";}
  if(spacegroupnumber==76) { // ------------------- 76  P4_{1} #76
    spacegroup="P4_{1}";}
  if(spacegroupnumber==77) { // ------------------- 77  P4_{2} #77
    spacegroup="P4_{2}";}
  if(spacegroupnumber==78) { // ------------------- 78  P4_{3} #78
    spacegroup="P4_{3}";}
  if(spacegroupnumber==79) { // ------------------- 79  I4 #79
    spacegroup="I4";}
  if(spacegroupnumber==80) { // ------------------- 80  I4_{1} #80
    spacegroup="I4_{1}";}
  if(spacegroupnumber==81) { // ------------------- 81  P-4 #81
    spacegroup="P-4";}
  if(spacegroupnumber==82) { // ------------------- 82  I-4 #82
    spacegroup="I-4";}
  if(spacegroupnumber==83) { // ------------------- 83  P4/m #83
    spacegroup="P4/m";}
  if(spacegroupnumber==84) { // ------------------- 84  P4_{2}/m #84
    spacegroup="P4_{2}/m";}
  if(spacegroupnumber==85) { // ------------------- 85  P4/n #85
    spacegroup="P4/n";}
  if(spacegroupnumber==86) { // ------------------- 86  P4_{2}/n #86
    spacegroup="P4_{2}/n";}
  if(spacegroupnumber==87) { // ------------------- 87  I4/m #87
    spacegroup="I4/m";}
  if(spacegroupnumber==88) { // ------------------- 88  I4_{1}/a #88
    spacegroup="I4_{1}/a";}
  if(spacegroupnumber==89) { // ------------------- 89  P422 #89
    spacegroup="P422";}
  if(spacegroupnumber==90) { // ------------------- 90  P42_{1}2 #90
    spacegroup="P42_{1}2";}
  if(spacegroupnumber==91) { // ------------------- 91  P4_{1}22 #91
    spacegroup="P4_{1}22";}
  if(spacegroupnumber==92) { // ------------------- 92  P4_{1}2_{1}2 #92
    spacegroup="P4_{1}2_{1}2";}
  if(spacegroupnumber==93) { // ------------------- 93  P4_{2}22 #93
    spacegroup="P4_{2}22";}
  if(spacegroupnumber==94) { // ------------------- 94  P4_{2}2_{1}2 #94
    spacegroup="P4_{2}2_{1}2";}
  if(spacegroupnumber==95) { // ------------------- 95  P4_{3}22 #95
    spacegroup="P4_{3}22";}
  if(spacegroupnumber==96) { // ------------------- 96  P4_{3}2_{1}2 #96
    spacegroup="P4_{3}2_{1}2";}
  if(spacegroupnumber==97) { // ------------------- 97  I422 #97
    spacegroup="I422";}
  if(spacegroupnumber==98) { // ------------------- 98  I4_{1}22 #98
    spacegroup="I4_{1}22";}
  if(spacegroupnumber==99) { // ------------------- 99  P4mm #99
    spacegroup="P4mm";}
  if(spacegroupnumber==100) { // ------------------- 100  P4bm #100
    spacegroup="P4bm";}
  if(spacegroupnumber==101) { // ------------------- 101  P4_{2}cm #101
    spacegroup="P4_{2}cm";}
  if(spacegroupnumber==102) { // ------------------- 102  P4_{2}nm #102
    spacegroup="P4_{2}nm";}
  if(spacegroupnumber==103) { // ------------------- 103  P4cc #103
    spacegroup="P4cc";}
  if(spacegroupnumber==104) { // ------------------- 104  P4nc #104
    spacegroup="P4nc";}
  if(spacegroupnumber==105) { // ------------------- 105  P4_{2}mc #105
    spacegroup="P4_{2}mc";}
  if(spacegroupnumber==106) { // ------------------- 106  P4_{2}bc #106
    spacegroup="P4_{2}bc";}
  if(spacegroupnumber==107) { // ------------------- 107  I4mm #107
    spacegroup="I4mm";}
  if(spacegroupnumber==108) { // ------------------- 108  I4cm #108
    spacegroup="I4cm";}
  if(spacegroupnumber==109) { // ------------------- 109  I4_{1}md #109
    spacegroup="I4_{1}md";}
  if(spacegroupnumber==110) { // ------------------- 110  I4_{1}cd #110
    spacegroup="I4_{1}cd";}
  if(spacegroupnumber==111) { // ------------------- 111  P-42m #111
    spacegroup="P-42m";}
  if(spacegroupnumber==112) { // ------------------- 112  P-42c #112
    spacegroup="P-42c";}
  if(spacegroupnumber==113) { // ------------------- 113  P-42_{1}m #113
    spacegroup="P-42_{1}m";}
  if(spacegroupnumber==114) { // ------------------- 114  P-42_{1}c #114
    spacegroup="P-42_{1}c";}
  if(spacegroupnumber==115) { // ------------------- 115  P-4m2 #115
    spacegroup="P-4m2";}
  if(spacegroupnumber==116) { // ------------------- 116  P-4c2 #116
    spacegroup="P-4c2";}
  if(spacegroupnumber==117) { // ------------------- 117  P-4b2 #117
    spacegroup="P-4b2";}
  if(spacegroupnumber==118) { // ------------------- 118  P-4n2 #118
    spacegroup="P-4n2";}
  if(spacegroupnumber==119) { // ------------------- 119  I-4m2 #119
    spacegroup="I-4m2";}
  if(spacegroupnumber==120) { // ------------------- 120  I-4c2 #120
    spacegroup="I-4c2";}
  if(spacegroupnumber==121) { // ------------------- 121  I-42m #121
    spacegroup="I-42m";}
  if(spacegroupnumber==122) { // ------------------- 122  I-42d #122
    spacegroup="I-42d";}
  if(spacegroupnumber==123) { // ------------------- 123  P4/mmm #123
    spacegroup="P4/mmm";}
  if(spacegroupnumber==124) { // ------------------- 124  P4/mcc #124
    spacegroup="P4/mcc";}
  if(spacegroupnumber==125) { // ------------------- 125  P4/nbm #125
    spacegroup="P4/nbm";}
  if(spacegroupnumber==126) { // ------------------- 126  P4/nnc #126
    spacegroup="P4/nnc";}
  if(spacegroupnumber==127) { // ------------------- 127  P4/mbm #127
    spacegroup="P4/mbm";}
  if(spacegroupnumber==128) { // ------------------- 128  P4/mnc #128
    spacegroup="P4/mnc";}
  if(spacegroupnumber==129) { // ------------------- 129  P4/nmm #129
    spacegroup="P4/nmm";}
  if(spacegroupnumber==130) { // ------------------- 130  P4/ncc #130
    spacegroup="P4/ncc";}
  if(spacegroupnumber==131) { // ------------------- 131  P4_{2}/mmc #131
    spacegroup="P4_{2}/mmc";}
  if(spacegroupnumber==132) { // ------------------- 132  P4_{2}/mcm #132
    spacegroup="P4_{2}/mcm";}
  if(spacegroupnumber==133) { // ------------------- 133  P4_{2}/nbc #133
    spacegroup="P4_{2}/nbc";}
  if(spacegroupnumber==134) { // ------------------- 134  P4_{2}/nnm #134
    spacegroup="P4_{2}/nnm";}
  if(spacegroupnumber==135) { // ------------------- 135  P4_{2}/mbc #135
    spacegroup="P4_{2}/mbc";}
  if(spacegroupnumber==136) { // ------------------- 136  P4_{2}/mnm #136
    spacegroup="P4_{2}/mnm";}
  if(spacegroupnumber==137) { // ------------------- 137  P4_{2}/nmc #137
    spacegroup="P4_{2}/nmc";}
  if(spacegroupnumber==138) { // ------------------- 138  P4_{2}/ncm #138
    spacegroup="P4_{2}/ncm";}
  if(spacegroupnumber==139) { // ------------------- 139  I4/mmm #139
    spacegroup="I4/mmm";}
  if(spacegroupnumber==140) { // ------------------- 140  I4/mcm #140
    spacegroup="I4/mcm";}
  if(spacegroupnumber==141) { // ------------------- 141  I4_{1}/amd #141
    spacegroup="I4_{1}/amd";}
  if(spacegroupnumber==142) { // ------------------- 142  I4_{1}/acd #142
    spacegroup="I4_{1}/acd";}
  if(spacegroupnumber==143) { // ------------------- 143  P3 #143
    spacegroup="P3";}
  if(spacegroupnumber==144) { // ------------------- 144  P3_{1} #144
    spacegroup="P3_{1}";}
  if(spacegroupnumber==145) { // ------------------- 145  P3_{2} #145
    spacegroup="P3_{2}";}
  if(spacegroupnumber==146) { // ------------------- 146  R3 #146
    spacegroup="R3";}
  if(spacegroupnumber==147) { // ------------------- 147  P-3 #147
    spacegroup="P-3";}
  if(spacegroupnumber==148) { // ------------------- 148  R-3 #148
    spacegroup="R-3";}
  if(spacegroupnumber==149) { // ------------------- 149  P312 #149
    spacegroup="P312";}
  if(spacegroupnumber==150) { // ------------------- 150  P321 #150
    spacegroup="P321";}
  if(spacegroupnumber==151) { // ------------------- 151  P3_{1}12 #151
    spacegroup="P3_{1}12";}
  if(spacegroupnumber==152) { // ------------------- 152  P3_{1}21 #152
    spacegroup="P3_{1}21";}
  if(spacegroupnumber==153) { // ------------------- 153  P3_{2}12 #153
    spacegroup="P3_{2}12";}
  if(spacegroupnumber==154) { // ------------------- 154  P3_{2}21 #154
    spacegroup="P3_{2}21";}
  if(spacegroupnumber==155) { // ------------------- 155  R32 #155
    spacegroup="R32";}
  if(spacegroupnumber==156) { // ------------------- 156  P3m1 #156
    spacegroup="P3m1";}
  if(spacegroupnumber==157) { // ------------------- 157  P31m #157
    spacegroup="P31m";}
  if(spacegroupnumber==158) { // ------------------- 158  P3c1 #158
    spacegroup="P3c1";}
  if(spacegroupnumber==159) { // ------------------- 159  P31c #159
    spacegroup="P31c";}
  if(spacegroupnumber==160) { // ------------------- 160  R3m #160
    spacegroup="R3m";}
  if(spacegroupnumber==161) { // ------------------- 161  R3c #161
    spacegroup="R3c";}
  if(spacegroupnumber==162) { // ------------------- 162  P-31m #162
    spacegroup="P-31m";}
  if(spacegroupnumber==163) { // ------------------- 163  P-31c #163
    spacegroup="P-31c";}
  if(spacegroupnumber==164) { // ------------------- 164  P-3m1 #164
    spacegroup="P-3m1";}
  if(spacegroupnumber==165) { // ------------------- 165  P-3c1 #165
    spacegroup="P-3c1";}
  if(spacegroupnumber==166) { // ------------------- 166  R-3m #166
    spacegroup="R-3m";}
  if(spacegroupnumber==167) { // ------------------- 167  R-3c #167
    spacegroup="R-3c";}
  if(spacegroupnumber==168) { // ------------------- 168  P6 #168
    spacegroup="P6";}
  if(spacegroupnumber==169) { // ------------------- 169  P6_{1} #169
    spacegroup="P6_{1}";}
  if(spacegroupnumber==170) { // ------------------- 170  P6_{5} #170
    spacegroup="P6_{5}";}
  if(spacegroupnumber==171) { // ------------------- 171  P6_{2} #171
    spacegroup="P6_{2}";}
  if(spacegroupnumber==172) { // ------------------- 172  P6_{4} #172
    spacegroup="P6_{4}";}
  if(spacegroupnumber==173) { // ------------------- 173  P6_{3} #173
    spacegroup="P6_{3}";}
  if(spacegroupnumber==174) { // ------------------- 174  P-6 #174
    spacegroup="P-6";}
  if(spacegroupnumber==175) { // ------------------- 175  P6/m #175
    spacegroup="P6/m";}
  if(spacegroupnumber==176) { // ------------------- 176  P6_{3}/m #176
    spacegroup="P6_{3}/m";}
  if(spacegroupnumber==177) { // ------------------- 177  P622 #177
    spacegroup="P622";}
  if(spacegroupnumber==178) { // ------------------- 178  P6_{1}22 #178
    spacegroup="P6_{1}22";}
  if(spacegroupnumber==179) { // ------------------- 179  P6_{5}22 #179
    spacegroup="P6_{5}22";}
  if(spacegroupnumber==180) { // ------------------- 180  P6_{2}22 #180
    spacegroup="P6_{2}22";}
  if(spacegroupnumber==181) { // ------------------- 181  P6_{4}22 #181
    spacegroup="P6_{4}22";}
  if(spacegroupnumber==182) { // ------------------- 182  P6_{3}22 #182
    spacegroup="P6_{3}22";}
  if(spacegroupnumber==183) { // ------------------- 183  P6mm #183
    spacegroup="P6mm";}
  if(spacegroupnumber==184) { // ------------------- 184  P6cc #184
    spacegroup="P6cc";}
  if(spacegroupnumber==185) { // ------------------- 185  P6_{3}cm #185
    spacegroup="P6_{3}cm";}
  if(spacegroupnumber==186) { // ------------------- 186  P6_{3}mc #186
    spacegroup="P6_{3}mc";}
  if(spacegroupnumber==187) { // ------------------- 187  P-6m2 #187
    spacegroup="P-6m2";}
  if(spacegroupnumber==188) { // ------------------- 188  P-6c2 #188
    spacegroup="P-6c2";}
  if(spacegroupnumber==189) { // ------------------- 189  P-62m #189
    spacegroup="P-62m";}
  if(spacegroupnumber==190) { // ------------------- 190  P-62c #190
    spacegroup="P-62c";}
  if(spacegroupnumber==191) { // ------------------- 191  P6/mmm #191
    spacegroup="P6/mmm";}
  if(spacegroupnumber==192) { // ------------------- 192  P6/mcc #192
    spacegroup="P6/mcc";}
  if(spacegroupnumber==193) { // ------------------- 193  P6_{3}/mcm #193
    spacegroup="P6_{3}/mcm";}
  if(spacegroupnumber==194) { // ------------------- 194  P6_{3}/mmc #194
    spacegroup="P6_{3}/mmc";}
  if(spacegroupnumber==195) { // ------------------- 195  P23 #195
    spacegroup="P23";}
  if(spacegroupnumber==196) { // ------------------- 196  F23 #196
    spacegroup="F23";}
  if(spacegroupnumber==197) { // ------------------- 197  I23 #197
    spacegroup="I23";}
  if(spacegroupnumber==198) { // ------------------- 198  P2_{1}3 #198
    spacegroup="P2_{1}3";}
  if(spacegroupnumber==199) { // ------------------- 199  I2_{1}3 #199
    spacegroup="I2_{1}3";}
  if(spacegroupnumber==200) { // ------------------- 200  Pm-3 #200
    spacegroup="Pm-3";}
  if(spacegroupnumber==201) { // ------------------- 201  Pn-3 #201
    spacegroup="Pn-3";}
  if(spacegroupnumber==202) { // ------------------- 202  Fm-3 #202
    spacegroup="Fm-3";}
  if(spacegroupnumber==203) { // ------------------- 203  Fd-3 #203
    spacegroup="Fd-3";}
  if(spacegroupnumber==204) { // ------------------- 204  Im-3 #204
    spacegroup="Im-3";}
  if(spacegroupnumber==205) { // ------------------- 205  Pa-3 #205
    spacegroup="Pa-3";}
  if(spacegroupnumber==206) { // ------------------- 206  Ia-3 #206
    spacegroup="Ia-3";}
  if(spacegroupnumber==207) { // ------------------- 207  P432 #207
    spacegroup="P432";}
  if(spacegroupnumber==208) { // ------------------- 208  P4_{2}32 #208
    spacegroup="P4_{2}32";}
  if(spacegroupnumber==209) { // ------------------- 209  F432 #209
    spacegroup="F432";}
  if(spacegroupnumber==210) { // ------------------- 210  F4_{1}32 #210
    spacegroup="F4_{1}32";}
  if(spacegroupnumber==211) { // ------------------- 211  I432 #211
    spacegroup="I432";}
  if(spacegroupnumber==212) { // ------------------- 212  P4_{3}32 #212
    spacegroup="P4_{3}32";}
  if(spacegroupnumber==213) { // ------------------- 213  P4_{1}32 #213
    spacegroup="P4_{1}32";}
  if(spacegroupnumber==214) { // ------------------- 214  I4_{1}32 #214
    spacegroup="I4_{1}32";}
  if(spacegroupnumber==215) { // ------------------- 215  P-43m #215
    spacegroup="P-43m";}
  if(spacegroupnumber==216) { // ------------------- 216  F-43m #216
    spacegroup="F-43m";}
  if(spacegroupnumber==217) { // ------------------- 217  I-43m #217
    spacegroup="I-43m";}
  if(spacegroupnumber==218) { // ------------------- 218  P-43n #218
    spacegroup="P-43n";}
  if(spacegroupnumber==219) { // ------------------- 219  F-43c #219
    spacegroup="F-43c";}
  if(spacegroupnumber==220) { // ------------------- 220  I-43d #220
    spacegroup="I-43d";}
  if(spacegroupnumber==221) { // ------------------- 221  Pm-3m #221
    spacegroup="Pm-3m";}
  if(spacegroupnumber==222) { // ------------------- 222  Pn-3n #222
    spacegroup="Pn-3n";}
  if(spacegroupnumber==223) { // ------------------- 223  Pm-3n #223
    spacegroup="Pm-3n";}
  if(spacegroupnumber==224) { // ------------------- 224  Pn-3m #224
    spacegroup="Pn-3m";}
  if(spacegroupnumber==225) { // ------------------- 225  Fm-3m #225
    spacegroup="Fm-3m";}
  if(spacegroupnumber==226) { // ------------------- 226  Fm-3c #226
    spacegroup="Fm-3c";}
  if(spacegroupnumber==227) { // ------------------- 227  Fd-3m #227
    spacegroup="Fd-3m";}
  if(spacegroupnumber==228) { // ------------------- 228  Fd-3c #228
    spacegroup="Fd-3c";}
  if(spacegroupnumber==229) { // ------------------- 229  Im-3m #229
    spacegroup="Im-3m";}
  if(spacegroupnumber==230) { // ------------------- 230  Ia-3d #230
    spacegroup="Ia-3d";}
  // done
  return spacegroup;
}

// ***************************************************************************
// GetSpaceGroupSchoenflies
// ***************************************************************************
string GetSpaceGroupSchoenflies(int spacegroupnumber) {
  string spacegroup;
  if(spacegroupnumber < 1 || spacegroupnumber > 230) {
    cerr << "aflow_xatom.cpp GetSpaceGroupSchoenflies() routine: " << endl;
    cerr << "  space group specified invalid (1-230): " << spacegroupnumber << endl;
    exit(0);
  }
  // OK
  if(spacegroupnumber==1) { // ------------------- 1  C_{1}^{1} #1
    spacegroup="C_{1}^{1}";}
  if(spacegroupnumber==2) { // ------------------- 2  C_{i}^{1} #2
    spacegroup="C_{i}^{1}";}
  if(spacegroupnumber==3) { // ------------------- 3  C_{2}^{1} #3
    spacegroup="C_{2}^{1}";}
  if(spacegroupnumber==4) { // ------------------- 4  C_{2}^{2} #4
    spacegroup="C_{2}^{2}";}
  if(spacegroupnumber==5) { // ------------------- 5  C_{2}^{3} #5
    spacegroup="C_{2}^{3}";}
  if(spacegroupnumber==6) { // ------------------- 6  C_{s}^{1} #6
    spacegroup="C_{s}^{1}";}
  if(spacegroupnumber==7) { // ------------------- 7  C_{s}^{2} #7
    spacegroup="C_{s}^{2}";}
  if(spacegroupnumber==8) { // ------------------- 8  C_{s}^{3} #8
    spacegroup="C_{s}^{3}";}
  if(spacegroupnumber==9) { // ------------------- 9  C_{s}^{4} #9
    spacegroup="C_{s}^{4}";}
  if(spacegroupnumber==10) { // ------------------- 10  C_{2h}^{1} #10
    spacegroup="C_{2h}^{1}";}
  if(spacegroupnumber==11) { // ------------------- 11  C_{2h}^{2} #11
    spacegroup="C_{2h}^{2}";}
  if(spacegroupnumber==12) { // ------------------- 12  C_{2h}^{3} #12
    spacegroup="C_{2h}^{3}";}
  if(spacegroupnumber==13) { // ------------------- 13  C_{2h}^{4} #13
    spacegroup="C_{2h}^{4}";}
  if(spacegroupnumber==14) { // ------------------- 14  C_{2h}^{5} #14
    spacegroup="C_{2h}^{5}";}
  if(spacegroupnumber==15) { // ------------------- 15  C_{2h}^{6} #15
    spacegroup="C_{2h}^{6}";}
  if(spacegroupnumber==16) { // ------------------- 16  D_{2}^{1} #16
    spacegroup="D_{2}^{1}";}
  if(spacegroupnumber==17) { // ------------------- 17  D_{2}^{2} #17
    spacegroup="D_{2}^{2}";}
  if(spacegroupnumber==18) { // ------------------- 18  D_{2}^{3} #18
    spacegroup="D_{2}^{3}";}
  if(spacegroupnumber==19) { // ------------------- 19  D_{2}^{4} #19
    spacegroup="D_{2}^{4}";}
  if(spacegroupnumber==20) { // ------------------- 20  D_{2}^{5} #20
    spacegroup="D_{2}^{5}";}
  if(spacegroupnumber==21) { // ------------------- 21  D_{2}^{6} #21
    spacegroup="D_{2}^{6}";}
  if(spacegroupnumber==22) { // ------------------- 22  D_{2}^{7} #22
    spacegroup="D_{2}^{7}";}
  if(spacegroupnumber==23) { // ------------------- 23  D_{2}^{8} #23
    spacegroup="D_{2}^{8}";}
  if(spacegroupnumber==24) { // ------------------- 24  D_{2}^{9} #24
    spacegroup="D_{2}^{9}";}
  if(spacegroupnumber==25) { // ------------------- 25  C_{2v}^{1} #25
    spacegroup="C_{2v}^{1}";}
  if(spacegroupnumber==26) { // ------------------- 26  C_{2v}^{2} #26
    spacegroup="C_{2v}^{2}";}
  if(spacegroupnumber==27) { // ------------------- 27  C_{2v}^{3} #27
    spacegroup="C_{2v}^{3}";}
  if(spacegroupnumber==28) { // ------------------- 28  C_{2v}^{4} #28
    spacegroup="C_{2v}^{4}";}
  if(spacegroupnumber==29) { // ------------------- 29  C_{2v}^{5} #29
    spacegroup="C_{2v}^{5}";}
  if(spacegroupnumber==30) { // ------------------- 30  C_{2v}^{6} #30
    spacegroup="C_{2v}^{6}";}
  if(spacegroupnumber==31) { // ------------------- 31  C_{2v}^{7} #31
    spacegroup="C_{2v}^{7}";}
  if(spacegroupnumber==32) { // ------------------- 32  C_{2v}^{8} #32
    spacegroup="C_{2v}^{8}";}
  if(spacegroupnumber==33) { // ------------------- 33  C_{2v}^{9} #33
    spacegroup="C_{2v}^{9}";}
  if(spacegroupnumber==34) { // ------------------- 34  C_{2v}^{10} #34
    spacegroup="C_{2v}^{10}";}
  if(spacegroupnumber==35) { // ------------------- 35  C_{2v}^{11} #35
    spacegroup="C_{2v}^{11}";}
  if(spacegroupnumber==36) { // ------------------- 36  C_{2v}^{12} #36
    spacegroup="C_{2v}^{12}";}
  if(spacegroupnumber==37) { // ------------------- 37  C_{2v}^{13} #37
    spacegroup="C_{2v}^{13}";}
  if(spacegroupnumber==38) { // ------------------- 38  C_{2v}^{14} #38
    spacegroup="C_{2v}^{14}";}
  if(spacegroupnumber==39) { // ------------------- 39  C_{2v}^{15} #39
    spacegroup="C_{2v}^{15}";}
  if(spacegroupnumber==40) { // ------------------- 40  C_{2v}^{16} #40
    spacegroup="C_{2v}^{16}";}
  if(spacegroupnumber==41) { // ------------------- 41  C_{2v}^{17} #41
    spacegroup="C_{2v}^{17}";}
  if(spacegroupnumber==42) { // ------------------- 42  C_{2v}^{18} #42
    spacegroup="C_{2v}^{18}";}
  if(spacegroupnumber==43) { // ------------------- 43  C_{2v}^{19} #43
    spacegroup="C_{2v}^{19}";}
  if(spacegroupnumber==44) { // ------------------- 44  C_{2v}^{20} #44
    spacegroup="C_{2v}^{20}";}
  if(spacegroupnumber==45) { // ------------------- 45  C_{2v}^{21} #45
    spacegroup="C_{2v}^{21}";}
  if(spacegroupnumber==46) { // ------------------- 46  C_{2v}^{22} #46
    spacegroup="C_{2v}^{22}";}
  if(spacegroupnumber==47) { // ------------------- 47  D_{2h}^{1} #47
    spacegroup="D_{2h}^{1}";}
  if(spacegroupnumber==48) { // ------------------- 48  D_{2h}^{2} #48
    spacegroup="D_{2h}^{2}";}
  if(spacegroupnumber==49) { // ------------------- 49  D_{2h}^{3} #49
    spacegroup="D_{2h}^{3}";}
  if(spacegroupnumber==50) { // ------------------- 50  D_{2h}^{4} #50
    spacegroup="D_{2h}^{4}";}
  if(spacegroupnumber==51) { // ------------------- 51  D_{2h}^{5} #51
    spacegroup="D_{2h}^{5}";}
  if(spacegroupnumber==52) { // ------------------- 52  D_{2h}^{6} #52
    spacegroup="D_{2h}^{6}";}
  if(spacegroupnumber==53) { // ------------------- 53  D_{2h}^{7} #53
    spacegroup="D_{2h}^{7}";}
  if(spacegroupnumber==54) { // ------------------- 54  D_{2h}^{8} #54
    spacegroup="D_{2h}^{8}";}
  if(spacegroupnumber==55) { // ------------------- 55  D_{2h}^{9} #55
    spacegroup="D_{2h}^{9}";}
  if(spacegroupnumber==56) { // ------------------- 56  D_{2h}^{10} #56
    spacegroup="D_{2h}^{10}";}
  if(spacegroupnumber==57) { // ------------------- 57  D_{2h}^{11} #57
    spacegroup="D_{2h}^{11}";}
  if(spacegroupnumber==58) { // ------------------- 58  D_{2h}^{12} #58
    spacegroup="D_{2h}^{12}";}
  if(spacegroupnumber==59) { // ------------------- 59  D_{2h}^{13} #59
    spacegroup="D_{2h}^{13}";}
  if(spacegroupnumber==60) { // ------------------- 60  D_{2h}^{14} #60
    spacegroup="D_{2h}^{14}";}
  if(spacegroupnumber==61) { // ------------------- 61  D_{2h}^{15} #61
    spacegroup="D_{2h}^{15}";}
  if(spacegroupnumber==62) { // ------------------- 62  D_{2h}^{16} #62
    spacegroup="D_{2h}^{16}";}
  if(spacegroupnumber==63) { // ------------------- 63  D_{2h}^{17} #63
    spacegroup="D_{2h}^{17}";}
  if(spacegroupnumber==64) { // ------------------- 64  D_{2h}^{18} #64
    spacegroup="D_{2h}^{18}";}
  if(spacegroupnumber==65) { // ------------------- 65  D_{2h}^{19} #65
    spacegroup="D_{2h}^{19}";}
  if(spacegroupnumber==66) { // ------------------- 66  D_{2h}^{20} #66
    spacegroup="D_{2h}^{20}";}
  if(spacegroupnumber==67) { // ------------------- 67  D_{2h}^{21} #67
    spacegroup="D_{2h}^{21}";}
  if(spacegroupnumber==68) { // ------------------- 68  D_{2h}^{22} #68
    spacegroup="D_{2h}^{22}";}
  if(spacegroupnumber==69) { // ------------------- 69  D_{2h}^{23} #69
    spacegroup="D_{2h}^{23}";}
  if(spacegroupnumber==70) { // ------------------- 70  D_{2h}^{24} #70
    spacegroup="D_{2h}^{24}";}
  if(spacegroupnumber==71) { // ------------------- 71  D_{2h}^{25} #71
    spacegroup="D_{2h}^{25}";}
  if(spacegroupnumber==72) { // ------------------- 72  D_{2h}^{26} #72
    spacegroup="D_{2h}^{26}";}
  if(spacegroupnumber==73) { // ------------------- 73  D_{2h}^{27} #73
    spacegroup="D_{2h}^{27}";}
  if(spacegroupnumber==74) { // ------------------- 74  D_{2h}^{28} #74
    spacegroup="D_{2h}^{28}";}
  if(spacegroupnumber==75) { // ------------------- 75  C_{4}^{1} #75
    spacegroup="C_{4}^{1}";}
  if(spacegroupnumber==76) { // ------------------- 76  C_{4}^{2} #76
    spacegroup="C_{4}^{2}";}
  if(spacegroupnumber==77) { // ------------------- 77  C_{4}^{3} #77
    spacegroup="C_{4}^{3}";}
  if(spacegroupnumber==78) { // ------------------- 78  C_{4}^{4} #78
    spacegroup="C_{4}^{4}";}
  if(spacegroupnumber==79) { // ------------------- 79  C_{4}^{5} #79
    spacegroup="C_{4}^{5}";}
  if(spacegroupnumber==80) { // ------------------- 80  C_{4}^{6} #80
    spacegroup="C_{4}^{6}";}
  if(spacegroupnumber==81) { // ------------------- 81  S_{4}^{1} #81
    spacegroup="S_{4}^{1}";}
  if(spacegroupnumber==82) { // ------------------- 82  S_{4}^{2} #82
    spacegroup="S_{4}^{2}";}
  if(spacegroupnumber==83) { // ------------------- 83  C_{4h}^{1} #83
    spacegroup="C_{4h}^{1}";}
  if(spacegroupnumber==84) { // ------------------- 84  C_{4h}^{2} #84
    spacegroup="C_{4h}^{2}";}
  if(spacegroupnumber==85) { // ------------------- 85  C_{4h}^{3} #85
    spacegroup="C_{4h}^{3}";}
  if(spacegroupnumber==86) { // ------------------- 86  C_{4h}^{4} #86
    spacegroup="C_{4h}^{4}";}
  if(spacegroupnumber==87) { // ------------------- 87  C_{4h}^{5} #87
    spacegroup="C_{4h}^{5}";}
  if(spacegroupnumber==88) { // ------------------- 88  C_{4h}^{6} #88
    spacegroup="C_{4h}^{6}";}
  if(spacegroupnumber==89) { // ------------------- 89  D_{4}^{1} #89
    spacegroup="D_{4}^{1}";}
  if(spacegroupnumber==90) { // ------------------- 90  D_{4}^{2} #90
    spacegroup="D_{4}^{2}";}
  if(spacegroupnumber==91) { // ------------------- 91  D_{4}^{3} #91
    spacegroup="D_{4}^{3}";}
  if(spacegroupnumber==92) { // ------------------- 92  D_{4}^{4} #92
    spacegroup="D_{4}^{4}";}
  if(spacegroupnumber==93) { // ------------------- 93  D_{4}^{5} #93
    spacegroup="D_{4}^{5}";}
  if(spacegroupnumber==94) { // ------------------- 94  D_{4}^{6} #94
    spacegroup="D_{4}^{6}";}
  if(spacegroupnumber==95) { // ------------------- 95  D_{4}^{7} #95
    spacegroup="D_{4}^{7}";}
  if(spacegroupnumber==96) { // ------------------- 96  D_{4}^{8} #96
    spacegroup="D_{4}^{8}";}
  if(spacegroupnumber==97) { // ------------------- 97  D_{4}^{9} #97
    spacegroup="D_{4}^{9}";}
  if(spacegroupnumber==98) { // ------------------- 98  D_{4}^{10} #98
    spacegroup="D_{4}^{10}";}
  if(spacegroupnumber==99) { // ------------------- 99  C_{4v}^{1} #99
    spacegroup="C_{4v}^{1}";}
  if(spacegroupnumber==100) { // ------------------- 100  C_{4v}^{2} #100
    spacegroup="C_{4v}^{2}";}
  if(spacegroupnumber==101) { // ------------------- 101  C_{4v}^{3} #101
    spacegroup="C_{4v}^{3}";}
  if(spacegroupnumber==102) { // ------------------- 102  C_{4v}^{4} #102
    spacegroup="C_{4v}^{4}";}
  if(spacegroupnumber==103) { // ------------------- 103  C_{4v}^{5} #103
    spacegroup="C_{4v}^{5}";}
  if(spacegroupnumber==104) { // ------------------- 104  C_{4v}^{6} #104
    spacegroup="C_{4v}^{6}";}
  if(spacegroupnumber==105) { // ------------------- 105  C_{4v}^{7} #105
    spacegroup="C_{4v}^{7}";}
  if(spacegroupnumber==106) { // ------------------- 106  C_{4v}^{8} #106
    spacegroup="C_{4v}^{8}";}
  if(spacegroupnumber==107) { // ------------------- 107  C_{4v}^{9} #107
    spacegroup="C_{4v}^{9}";}
  if(spacegroupnumber==108) { // ------------------- 108  C_{4v}^{10} #108
    spacegroup="C_{4v}^{10}";}
  if(spacegroupnumber==109) { // ------------------- 109  C_{4v}^{11} #109
    spacegroup="C_{4v}^{11}";}
  if(spacegroupnumber==110) { // ------------------- 110  C_{4v}^{12} #110
    spacegroup="C_{4v}^{12}";}
  if(spacegroupnumber==111) { // ------------------- 111  D_{2d}^{1} #111
    spacegroup="D_{2d}^{1}";}
  if(spacegroupnumber==112) { // ------------------- 112  D_{2d}^{2} #112
    spacegroup="D_{2d}^{2}";}
  if(spacegroupnumber==113) { // ------------------- 113  D_{2d}^{3} #113
    spacegroup="D_{2d}^{3}";}
  if(spacegroupnumber==114) { // ------------------- 114  D_{2d}^{4} #114
    spacegroup="D_{2d}^{4}";}
  if(spacegroupnumber==115) { // ------------------- 115  D_{2d}^{5} #115
    spacegroup="D_{2d}^{5}";}
  if(spacegroupnumber==116) { // ------------------- 116  D_{2d}^{6} #116
    spacegroup="D_{2d}^{6}";}
  if(spacegroupnumber==117) { // ------------------- 117  D_{2d}^{7} #117
    spacegroup="D_{2d}^{7}";}
  if(spacegroupnumber==118) { // ------------------- 118  D_{2d}^{8} #118
    spacegroup="D_{2d}^{8}";}
  if(spacegroupnumber==119) { // ------------------- 119  D_{2d}^{9} #119
    spacegroup="D_{2d}^{9}";}
  if(spacegroupnumber==120) { // ------------------- 120  D_{2d}^{10} #120
    spacegroup="D_{2d}^{10}";}
  if(spacegroupnumber==121) { // ------------------- 121  D_{2d}^{11} #121
    spacegroup="D_{2d}^{11}";}
  if(spacegroupnumber==122) { // ------------------- 122  D_{2d}^{12} #122
    spacegroup="D_{2d}^{12}";}
  if(spacegroupnumber==123) { // ------------------- 123  D_{4h}^{1} #123
    spacegroup="D_{4h}^{1}";}
  if(spacegroupnumber==124) { // ------------------- 124  D_{4h}^{2} #124
    spacegroup="D_{4h}^{2}";}
  if(spacegroupnumber==125) { // ------------------- 125  D_{4h}^{3} #125
    spacegroup="D_{4h}^{3}";}
  if(spacegroupnumber==126) { // ------------------- 126  D_{4h}^{4} #126
    spacegroup="D_{4h}^{4}";}
  if(spacegroupnumber==127) { // ------------------- 127  D_{4h}^{5} #127
    spacegroup="D_{4h}^{5}";}
  if(spacegroupnumber==128) { // ------------------- 128  D_{4h}^{6} #128
    spacegroup="D_{4h}^{6}";}
  if(spacegroupnumber==129) { // ------------------- 129  D_{4h}^{7} #129
    spacegroup="D_{4h}^{7}";}
  if(spacegroupnumber==130) { // ------------------- 130  D_{4h}^{8} #130
    spacegroup="D_{4h}^{8}";}
  if(spacegroupnumber==131) { // ------------------- 131  D_{4h}^{9} #131
    spacegroup="D_{4h}^{9}";}
  if(spacegroupnumber==132) { // ------------------- 132  D_{4h}^{10} #132
    spacegroup="D_{4h}^{10}";}
  if(spacegroupnumber==133) { // ------------------- 133  D_{4h}^{11} #133
    spacegroup="D_{4h}^{11}";}
  if(spacegroupnumber==134) { // ------------------- 134  D_{4h}^{12} #134
    spacegroup="D_{4h}^{12}";}
  if(spacegroupnumber==135) { // ------------------- 135  D_{4h}^{13} #135
    spacegroup="D_{4h}^{13}";}
  if(spacegroupnumber==136) { // ------------------- 136  D_{4h}^{14} #136
    spacegroup="D_{4h}^{14}";}
  if(spacegroupnumber==137) { // ------------------- 137  D_{4h}^{15} #137
    spacegroup="D_{4h}^{15}";}
  if(spacegroupnumber==138) { // ------------------- 138  D_{4h}^{16} #138
    spacegroup="D_{4h}^{16}";}
  if(spacegroupnumber==139) { // ------------------- 139  D_{4h}^{17} #139
    spacegroup="D_{4h}^{17}";}
  if(spacegroupnumber==140) { // ------------------- 140  D_{4h}^{18} #140
    spacegroup="D_{4h}^{18}";}
  if(spacegroupnumber==141) { // ------------------- 141  D_{4h}^{19} #141
    spacegroup="D_{4h}^{19}";}
  if(spacegroupnumber==142) { // ------------------- 142  D_{4h}^{20} #142
    spacegroup="D_{4h}^{20}";}
  if(spacegroupnumber==143) { // ------------------- 143  C_{3}^{1} #143
    spacegroup="C_{3}^{1}";}
  if(spacegroupnumber==144) { // ------------------- 144  C_{3}^{2} #144
    spacegroup="C_{3}^{2}";}
  if(spacegroupnumber==145) { // ------------------- 145  C_{3}^{3} #145
    spacegroup="C_{3}^{3}";}
  if(spacegroupnumber==146) { // ------------------- 146  C_{3}^{4} #146
    spacegroup="C_{3}^{4}";}
  if(spacegroupnumber==147) { // ------------------- 147  C_{3i}^{1} #147
    spacegroup="C_{3i}^{1}";}
  if(spacegroupnumber==148) { // ------------------- 148  C_{3i}^{2} #148
    spacegroup="C_{3i}^{2}";}
  if(spacegroupnumber==149) { // ------------------- 149  D_{3}^{1} #149
    spacegroup="D_{3}^{1}";}
  if(spacegroupnumber==150) { // ------------------- 150  D_{3}^{2} #150
    spacegroup="D_{3}^{2}";}
  if(spacegroupnumber==151) { // ------------------- 151  D_{3}^{3} #151
    spacegroup="D_{3}^{3}";}
  if(spacegroupnumber==152) { // ------------------- 152  D_{3}^{4} #152
    spacegroup="D_{3}^{4}";}
  if(spacegroupnumber==153) { // ------------------- 153  D_{3}^{5} #153
    spacegroup="D_{3}^{5}";}
  if(spacegroupnumber==154) { // ------------------- 154  D_{3}^{6} #154
    spacegroup="D_{3}^{6}";}
  if(spacegroupnumber==155) { // ------------------- 155  D_{3}^{7} #155
    spacegroup="D_{3}^{7}";}
  if(spacegroupnumber==156) { // ------------------- 156  C_{3v}^{1} #156
    spacegroup="C_{3v}^{1}";}
  if(spacegroupnumber==157) { // ------------------- 157  C_{3v}^{2} #157
    spacegroup="C_{3v}^{2}";}
  if(spacegroupnumber==158) { // ------------------- 158  C_{3v}^{3} #158
    spacegroup="C_{3v}^{3}";}
  if(spacegroupnumber==159) { // ------------------- 159  C_{3v}^{4} #159
    spacegroup="C_{3v}^{4}";}
  if(spacegroupnumber==160) { // ------------------- 160  C_{3v}^{5} #160
    spacegroup="C_{3v}^{5}";}
  if(spacegroupnumber==161) { // ------------------- 161  C_{3v}^{6} #161
    spacegroup="C_{3v}^{6}";}
  if(spacegroupnumber==162) { // ------------------- 162  D_{3d}^{1} #162
    spacegroup="D_{3d}^{1}";}
  if(spacegroupnumber==163) { // ------------------- 163  D_{3d}^{2} #163
    spacegroup="D_{3d}^{2}";}
  if(spacegroupnumber==164) { // ------------------- 164  D_{3d}^{3} #164
    spacegroup="D_{3d}^{3}";}
  if(spacegroupnumber==165) { // ------------------- 165  D_{3d}^{4} #165
    spacegroup="D_{3d}^{4}";}
  if(spacegroupnumber==166) { // ------------------- 166  D_{3d}^{5} #166
    spacegroup="D_{3d}^{5}";}
  if(spacegroupnumber==167) { // ------------------- 167  D_{3d}^{6} #167
    spacegroup="D_{3d}^{6}";}
  if(spacegroupnumber==168) { // ------------------- 168  C_{6}^{1} #168
    spacegroup="C_{6}^{1}";}
  if(spacegroupnumber==169) { // ------------------- 169  C_{6}^{2} #169
    spacegroup="C_{6}^{2}";}
  if(spacegroupnumber==170) { // ------------------- 170  C_{6}^{3} #170
    spacegroup="C_{6}^{3}";}
  if(spacegroupnumber==171) { // ------------------- 171  C_{6}^{4} #171
    spacegroup="C_{6}^{4}";}
  if(spacegroupnumber==172) { // ------------------- 172  C_{6}^{5} #172
    spacegroup="C_{6}^{5}";}
  if(spacegroupnumber==173) { // ------------------- 173  C_{6}^{6} #173
    spacegroup="C_{6}^{6}";}
  if(spacegroupnumber==174) { // ------------------- 174  C_{3h}^{1} #174
    spacegroup="C_{3h}^{1}";}
  if(spacegroupnumber==175) { // ------------------- 175  C_{6h}^{1} #175
    spacegroup="C_{6h}^{1}";}
  if(spacegroupnumber==176) { // ------------------- 176  C_{6h}^{2} #176
    spacegroup="C_{6h}^{2}";}
  if(spacegroupnumber==177) { // ------------------- 177  D_{6}^{1} #177
    spacegroup="D_{6}^{1}";}
  if(spacegroupnumber==178) { // ------------------- 178  D_{6}^{2} #178
    spacegroup="D_{6}^{2}";}
  if(spacegroupnumber==179) { // ------------------- 179  D_{6}^{3} #179
    spacegroup="D_{6}^{3}";}
  if(spacegroupnumber==180) { // ------------------- 180  D_{6}^{4} #180
    spacegroup="D_{6}^{4}";}
  if(spacegroupnumber==181) { // ------------------- 181  D_{6}^{5} #181
    spacegroup="D_{6}^{5}";}
  if(spacegroupnumber==182) { // ------------------- 182  D_{6}^{6} #182
    spacegroup="D_{6}^{6}";}
  if(spacegroupnumber==183) { // ------------------- 183  C_{6v}^{1} #183
    spacegroup="C_{6v}^{1}";}
  if(spacegroupnumber==184) { // ------------------- 184  C_{6v}^{2} #184
    spacegroup="C_{6v}^{2}";}
  if(spacegroupnumber==185) { // ------------------- 185  C_{6v}^{3} #185
    spacegroup="C_{6v}^{3}";}
  if(spacegroupnumber==186) { // ------------------- 186  C_{6v}^{4} #186
    spacegroup="C_{6v}^{4}";}
  if(spacegroupnumber==187) { // ------------------- 187  D_{3h}^{1} #187
    spacegroup="D_{3h}^{1}";}
  if(spacegroupnumber==188) { // ------------------- 188  D_{3h}^{2} #188
    spacegroup="D_{3h}^{2}";}
  if(spacegroupnumber==189) { // ------------------- 189  D_{3h}^{3} #189
    spacegroup="D_{3h}^{3}";}
  if(spacegroupnumber==190) { // ------------------- 190  D_{3h}^{4} #190
    spacegroup="D_{3h}^{4}";}
  if(spacegroupnumber==191) { // ------------------- 191  D_{6h}^{1} #191
    spacegroup="D_{6h}^{1}";}
  if(spacegroupnumber==192) { // ------------------- 192  D_{6h}^{2} #192
    spacegroup="D_{6h}^{2}";}
  if(spacegroupnumber==193) { // ------------------- 193  D_{6h}^{3} #193
    spacegroup="D_{6h}^{3}";}
  if(spacegroupnumber==194) { // ------------------- 194  D_{6h}^{4} #194
    spacegroup="D_{6h}^{4}";}
  if(spacegroupnumber==195) { // ------------------- 195  T^{1} #195
    spacegroup="T^{1}";}
  if(spacegroupnumber==196) { // ------------------- 196  T^{2} #196
    spacegroup="T^{2}";}
  if(spacegroupnumber==197) { // ------------------- 197  T^{3} #197
    spacegroup="T^{3}";}
  if(spacegroupnumber==198) { // ------------------- 198  T^{4} #198
    spacegroup="T^{4}";}
  if(spacegroupnumber==199) { // ------------------- 199  T^{5} #199
    spacegroup="T^{5}";}
  if(spacegroupnumber==200) { // ------------------- 200  T_{h}^{1} #200
    spacegroup="T_{h}^{1}";}
  if(spacegroupnumber==201) { // ------------------- 201  T_{h}^{2} #201
    spacegroup="T_{h}^{2}";}
  if(spacegroupnumber==202) { // ------------------- 202  T_{h}^{3} #202
    spacegroup="T_{h}^{3}";}
  if(spacegroupnumber==203) { // ------------------- 203  T_{h}^{4} #203
    spacegroup="T_{h}^{4}";}
  if(spacegroupnumber==204) { // ------------------- 204  T_{h}^{5} #204
    spacegroup="T_{h}^{5}";}
  if(spacegroupnumber==205) { // ------------------- 205  T_{h}^{6} #205
    spacegroup="T_{h}^{6}";}
  if(spacegroupnumber==206) { // ------------------- 206  T_{h}^{7} #206
    spacegroup="T_{h}^{7}";}
  if(spacegroupnumber==207) { // ------------------- 207  O^{1} #207
    spacegroup="O^{1}";}
  if(spacegroupnumber==208) { // ------------------- 208  O^{2} #208
    spacegroup="O^{2}";}
  if(spacegroupnumber==209) { // ------------------- 209  O^{3} #209
    spacegroup="O^{3}";}
  if(spacegroupnumber==210) { // ------------------- 210  O^{4} #210
    spacegroup="O^{4}";}
  if(spacegroupnumber==211) { // ------------------- 211  O^{5} #211
    spacegroup="O^{5}";}
  if(spacegroupnumber==212) { // ------------------- 212  O^{6} #212
    spacegroup="O^{6}";}
  if(spacegroupnumber==213) { // ------------------- 213  O^{7} #213
    spacegroup="O^{7}";}
  if(spacegroupnumber==214) { // ------------------- 214  O^{8} #214
    spacegroup="O^{8}";}
  if(spacegroupnumber==215) { // ------------------- 215  T_{d}^{1} #215
    spacegroup="T_{d}^{1}";}
  if(spacegroupnumber==216) { // ------------------- 216  T_{d}^{2} #216
    spacegroup="T_{d}^{2}";}
  if(spacegroupnumber==217) { // ------------------- 217  T_{d}^{3} #217
    spacegroup="T_{d}^{3}";}
  if(spacegroupnumber==218) { // ------------------- 218  T_{d}^{4} #218
    spacegroup="T_{d}^{4}";}
  if(spacegroupnumber==219) { // ------------------- 219  T_{d}^{5} #219
    spacegroup="T_{d}^{5}";}
  if(spacegroupnumber==220) { // ------------------- 220  T_{d}^{6} #220
    spacegroup="T_{d}^{6}";}
  if(spacegroupnumber==221) { // ------------------- 221  O_{h}^{1} #221
    spacegroup="O_{h}^{1}";}
  if(spacegroupnumber==222) { // ------------------- 222  O_{h}^{2} #222
    spacegroup="O_{h}^{2}";}
  if(spacegroupnumber==223) { // ------------------- 223  O_{h}^{3} #223
    spacegroup="O_{h}^{3}";}
  if(spacegroupnumber==224) { // ------------------- 224  O_{h}^{4} #224
    spacegroup="O_{h}^{4}";}
  if(spacegroupnumber==225) { // ------------------- 225  O_{h}^{5} #225
    spacegroup="O_{h}^{5}";}
  if(spacegroupnumber==226) { // ------------------- 226  O_{h}^{6} #226
    spacegroup="O_{h}^{6}";}
  if(spacegroupnumber==227) { // ------------------- 227  O_{h}^{7} #227
    spacegroup="O_{h}^{7}";}
  if(spacegroupnumber==228) { // ------------------- 228  O_{h}^{8} #228
    spacegroup="O_{h}^{8}";}
  if(spacegroupnumber==229) { // ------------------- 229  O_{h}^{9} #229
    spacegroup="O_{h}^{9}";}
  if(spacegroupnumber==230) { // ------------------- 230  O_{h}^{10} #230
    spacegroup="O_{h}^{10}";}
  // done
  return spacegroup;
}

// ***************************************************************************
// GetSpaceGroupHall
// ***************************************************************************
string GetSpaceGroupHall(int spacegroupnumber) {
  // DX - Hall distinguishes space group setting.  This table assumes the first 
  //      setting that appears in the ITC.
  //      For more settings, they need to be hard-coded here.
  string spacegroup;
  if(spacegroupnumber < 1 || spacegroupnumber > 230) {
    cerr << "aflow_xatom.cpp GetSpaceGroupHall() routine: " << endl;
    cerr << "  space group specified invalid (1-230): " << spacegroupnumber << endl;
    exit(0);
  }
  // OK
  if(spacegroupnumber==1) { // ------------------- 1  P 1 #1
    spacegroup="P 1";}
  if(spacegroupnumber==2) { // ------------------- 2  -P 1 #2
    spacegroup="-P 1";}
  if(spacegroupnumber==3) { // ------------------- 3  P 2y #3
    spacegroup="P 2y";}
  if(spacegroupnumber==4) { // ------------------- 4  P 2yb #4
    spacegroup="P 2yb";}
  if(spacegroupnumber==5) { // ------------------- 5  C 2y #5
    spacegroup="C 2y";}
  if(spacegroupnumber==6) { // ------------------- 6  P -2y #6
    spacegroup="P -2y";}
  if(spacegroupnumber==7) { // ------------------- 7  P -2yc #7
    spacegroup="P -2yc";}
  if(spacegroupnumber==8) { // ------------------- 8  C -2y #8
    spacegroup="C -2y";}
  if(spacegroupnumber==9) { // ------------------- 9  C -2yc #9
    spacegroup="C -2yc";}
  if(spacegroupnumber==10) { // ------------------- 10  -P 2y #10
    spacegroup="-P 2y";}
  if(spacegroupnumber==11) { // ------------------- 11  -P 2yb #11
    spacegroup="-P 2yb";}
  if(spacegroupnumber==12) { // ------------------- 12  -C 2y #12
    spacegroup="-C 2y";}
  if(spacegroupnumber==13) { // ------------------- 13  -P 2yc #13
    spacegroup="-P 2yc";}
  if(spacegroupnumber==14) { // ------------------- 14  -P 2ybc #14
    spacegroup="-P 2ybc";}
  if(spacegroupnumber==15) { // ------------------- 15  -C 2yc #15
    spacegroup="-C 2yc";}
  if(spacegroupnumber==16) { // ------------------- 16  P 2 2 #16
    spacegroup="P 2 2";}
  if(spacegroupnumber==17) { // ------------------- 17  P 2c 2 #17
    spacegroup="P 2c 2";}
  if(spacegroupnumber==18) { // ------------------- 18  P 2 2ab #18
    spacegroup="P 2 2ab";}
  if(spacegroupnumber==19) { // ------------------- 19  P 2ac 2ab #19
    spacegroup="P 2ac 2ab";}
  if(spacegroupnumber==20) { // ------------------- 20  C 2c 2 #20
    spacegroup="C 2c 2";}
  if(spacegroupnumber==21) { // ------------------- 21  C 2 2 #21
    spacegroup="C 2 2";}
  if(spacegroupnumber==22) { // ------------------- 22  F 2 2 #22
    spacegroup="F 2 2";}
  if(spacegroupnumber==23) { // ------------------- 23  I 2 2 #23
    spacegroup="I 2 2";}
  if(spacegroupnumber==24) { // ------------------- 24  I 2b 2c #24
    spacegroup="I 2b 2c";}
  if(spacegroupnumber==25) { // ------------------- 25  P 2 -2 #25
    spacegroup="P 2 -2";}
  if(spacegroupnumber==26) { // ------------------- 26  P 2c -2 #26
    spacegroup="P 2c -2";}
  if(spacegroupnumber==27) { // ------------------- 27  P 2 -2c #27
    spacegroup="P 2 -2c";}
  if(spacegroupnumber==28) { // ------------------- 28  P 2 -2a #28
    spacegroup="P 2 -2a";}
  if(spacegroupnumber==29) { // ------------------- 29  P 2c -2ac #29
    spacegroup="P 2c -2ac";}
  if(spacegroupnumber==30) { // ------------------- 30  P 2 -2bc #30
    spacegroup="P 2 -2bc";}
  if(spacegroupnumber==31) { // ------------------- 31  P 2ac -2 #31
    spacegroup="P 2ac -2";}
  if(spacegroupnumber==32) { // ------------------- 32  P 2 -2ab #32
    spacegroup="P 2 -2ab";}
  if(spacegroupnumber==33) { // ------------------- 33  P 2c -2n #33
    spacegroup="P 2c -2n";}
  if(spacegroupnumber==34) { // ------------------- 34  P 2 -2n #34
    spacegroup="P 2 -2n";}
  if(spacegroupnumber==35) { // ------------------- 35  C 2 -2 #35
    spacegroup="C 2 -2";}
  if(spacegroupnumber==36) { // ------------------- 36  C 2c -2 #36
    spacegroup="C 2c -2";}
  if(spacegroupnumber==37) { // ------------------- 37  C 2 -2c #37
    spacegroup="C 2 -2c";}
  if(spacegroupnumber==38) { // ------------------- 38  A 2 -2 #38
    spacegroup="A 2 -2";}
  if(spacegroupnumber==39) { // ------------------- 39  A 2 -2c #39
    spacegroup="A 2 -2c";}
  if(spacegroupnumber==40) { // ------------------- 40  A 2 -2a #40
    spacegroup="A 2 -2a";}
  if(spacegroupnumber==41) { // ------------------- 41  A 2 -2ac #41
    spacegroup="A 2 -2ac";}
  if(spacegroupnumber==42) { // ------------------- 42  F 2 -2 #42
    spacegroup="F 2 -2";}
  if(spacegroupnumber==43) { // ------------------- 43  F 2 -2d #43
    spacegroup="F 2 -2d";}
  if(spacegroupnumber==44) { // ------------------- 44  I 2 -2 #44
    spacegroup="I 2 -2";}
  if(spacegroupnumber==45) { // ------------------- 45  I 2 -2c #45
    spacegroup="I 2 -2c";}
  if(spacegroupnumber==46) { // ------------------- 46  I 2 -2a #46
    spacegroup="I 2 -2a";}
  if(spacegroupnumber==47) { // ------------------- 47  -P 2 2 #47
    spacegroup="-P 2 2";}
  // setting 2 if(spacegroupnumber==48) { // ------------------- 48  -P 2ab 2bc #48
  // setting 2   spacegroup="-P 2ab 2bc";}
  if(spacegroupnumber==48) { // ------------------- 48  -P 2ab 2bc #48
    spacegroup="P 2 2 -1n";}
  if(spacegroupnumber==49) { // ------------------- 49  -P 2 2c #49
    spacegroup="-P 2 2c";}
  // setting 2 if(spacegroupnumber==50) { // ------------------- 50  -P 2ab 2b #50
  // setting 2   spacegroup="-P 2ab 2b";}
  if(spacegroupnumber==50) { // ------------------- 50  -P 2ab 2b #50
    spacegroup="P 2 2 -1ab";}
  if(spacegroupnumber==51) { // ------------------- 51  -P 2a 2a #51
    spacegroup="-P 2a 2a";}
  if(spacegroupnumber==52) { // ------------------- 52  -P 2a 2bc #52
    spacegroup="-P 2a 2bc";}
  if(spacegroupnumber==53) { // ------------------- 53  -P 2ac 2 #53
    spacegroup="-P 2ac 2";}
  if(spacegroupnumber==54) { // ------------------- 54  -P 2a 2ac #54
    spacegroup="-P 2a 2ac";}
  if(spacegroupnumber==55) { // ------------------- 55  -P 2 2ab #55
    spacegroup="-P 2 2ab";}
  if(spacegroupnumber==56) { // ------------------- 56  -P 2ab 2ac #56
    spacegroup="-P 2ab 2ac";}
  if(spacegroupnumber==57) { // ------------------- 57  -P 2c 2b #57
    spacegroup="-P 2c 2b";}
  if(spacegroupnumber==58) { // ------------------- 58  -P 2 2n #58
    spacegroup="-P 2 2n";}
  // setting 2 if(spacegroupnumber==59) { // ------------------- 59  -P 2ab 2a #59
  // setting 2   spacegroup="-P 2ab 2a";}
  if(spacegroupnumber==59) { // ------------------- 59  -P 2ab 2a #59
    spacegroup="P 2 2ab -1ab";}
  if(spacegroupnumber==60) { // ------------------- 60  -P 2n 2ab #60
    spacegroup="-P 2n 2ab";}
  if(spacegroupnumber==61) { // ------------------- 61  -P 2ac 2ab #61
    spacegroup="-P 2ac 2ab";}
  if(spacegroupnumber==62) { // ------------------- 62  -P 2ac 2n #62
    spacegroup="-P 2ac 2n";}
  if(spacegroupnumber==63) { // ------------------- 63  -C 2c 2 #63
    spacegroup="-C 2c 2";}
  if(spacegroupnumber==64) { // ------------------- 64  -C 2bc 2 #64
    spacegroup="-C 2bc 2";}
  if(spacegroupnumber==65) { // ------------------- 65  -C 2 2 #65
    spacegroup="-C 2 2";}
  if(spacegroupnumber==66) { // ------------------- 66  -C 2 2c #66
    spacegroup="-C 2 2c";}
  if(spacegroupnumber==67) { // ------------------- 67  -C 2b 2 #67
    spacegroup="-C 2b 2";}
  // setting 2 if(spacegroupnumber==68) { // ------------------- 68  -C 2b 2bc #68
  // setting 2   spacegroup="-C 2b 2bc";}
  if(spacegroupnumber==68) { // ------------------- 68  -C 2b 2bc #68
    spacegroup="C 2 2 -1bc";}
  if(spacegroupnumber==69) { // ------------------- 69  -F 2 2 #69
    spacegroup="-F 2 2";}
  // setting 2 if(spacegroupnumber==70) { // ------------------- 70  -F 2uv 2vw #70
  // setting 2   spacegroup="-F 2uv 2vw";}
  if(spacegroupnumber==70) { // ------------------- 70  -F 2uv 2vw #70
    spacegroup="F 2 2 -1d";}
  if(spacegroupnumber==71) { // ------------------- 71  -I 2 2 #71
    spacegroup="-I 2 2";}
  if(spacegroupnumber==72) { // ------------------- 72  -I 2 2c #72
    spacegroup="-I 2 2c";}
  if(spacegroupnumber==73) { // ------------------- 73  -I 2b 2c #73
    spacegroup="-I 2b 2c";}
  if(spacegroupnumber==74) { // ------------------- 74  -I 2b 2 #74
    spacegroup="-I 2b 2";}
  if(spacegroupnumber==75) { // ------------------- 75  P 4 #75
    spacegroup="P 4";}
  if(spacegroupnumber==76) { // ------------------- 76  P 4w #76
    spacegroup="P 4w";}
  if(spacegroupnumber==77) { // ------------------- 77  P 4c #77
    spacegroup="P 4c";}
  if(spacegroupnumber==78) { // ------------------- 78  P 4cw #78
    spacegroup="P 4cw";}
  if(spacegroupnumber==79) { // ------------------- 79  I 4 #79
    spacegroup="I 4";}
  if(spacegroupnumber==80) { // ------------------- 80  I 4bw #80
    spacegroup="I 4bw";}
  if(spacegroupnumber==81) { // ------------------- 81  P -4 #81
    spacegroup="P -4";}
  if(spacegroupnumber==82) { // ------------------- 82  I -4 #82
    spacegroup="I -4";}
  if(spacegroupnumber==83) { // ------------------- 83  -P 4 #83
    spacegroup="-P 4";}
  if(spacegroupnumber==84) { // ------------------- 84  -P 4c #84
    spacegroup="-P 4c";}
  // setting 2 if(spacegroupnumber==85) { // ------------------- 85  -P 4a #85
  // setting 2   spacegroup="-P 4a";}
  if(spacegroupnumber==85) { // ------------------- 85  -P 4a #85
    spacegroup="P 4ab -1ab";}
  // setting 2 if(spacegroupnumber==86) { // ------------------- 86  -P 4bc #86
  // setting 2   spacegroup="-P 4bc";}
  if(spacegroupnumber==86) { // ------------------- 86  -P 4bc #86
    spacegroup="P 4n -1n";}
  if(spacegroupnumber==87) { // ------------------- 87  -I 4 #87
    spacegroup="-I 4";}
  // setting 2 if(spacegroupnumber==88) { // ------------------- 88  -I 4ad #88
  // setting 2   spacegroup="-I 4ad";}
  if(spacegroupnumber==88) { // ------------------- 88  -I 4ad #88
    spacegroup="I 4bw -1bw";}
  if(spacegroupnumber==89) { // ------------------- 89  P 4 2 #89
    spacegroup="P 4 2";}
  if(spacegroupnumber==90) { // ------------------- 90  P 4ab 2ab #90
    spacegroup="P 4ab 2ab";}
  if(spacegroupnumber==91) { // ------------------- 91  P 4w 2c #91
    spacegroup="P 4w 2c";}
  if(spacegroupnumber==92) { // ------------------- 92  P 4abw 2nw #92
    spacegroup="P 4abw 2nw";}
  if(spacegroupnumber==93) { // ------------------- 93  P 4c 2 #93
    spacegroup="P 4c 2";}
  if(spacegroupnumber==94) { // ------------------- 94  P 4n 2n #94
    spacegroup="P 4n 2n";}
  if(spacegroupnumber==95) { // ------------------- 95  P 4cw 2c #95
    spacegroup="P 4cw 2c";}
  if(spacegroupnumber==96) { // ------------------- 96  P 4nw 2abw #96
    spacegroup="P 4nw 2abw";}
  if(spacegroupnumber==97) { // ------------------- 97  I 4 2 #97
    spacegroup="I 4 2";}
  if(spacegroupnumber==98) { // ------------------- 98  I 4bw 2bw #98
    spacegroup="I 4bw 2bw";}
  if(spacegroupnumber==99) { // ------------------- 99  P 4 -2 #99
    spacegroup="P 4 -2";}
  if(spacegroupnumber==100) { // ------------------- 100  P 4 -2ab #100
    spacegroup="P 4 -2ab";}
  if(spacegroupnumber==101) { // ------------------- 101  P 4c -2c #101
    spacegroup="P 4c -2c";}
  if(spacegroupnumber==102) { // ------------------- 102  P 4n -2n #102
    spacegroup="P 4n -2n";}
  if(spacegroupnumber==103) { // ------------------- 103  P 4 -2c #103
    spacegroup="P 4 -2c";}
  if(spacegroupnumber==104) { // ------------------- 104  P 4 -2n #104
    spacegroup="P 4 -2n";}
  if(spacegroupnumber==105) { // ------------------- 105  P 4c -2 #105
    spacegroup="P 4c -2";}
  if(spacegroupnumber==106) { // ------------------- 106  P 4c -2ab #106
    spacegroup="P 4c -2ab";}
  if(spacegroupnumber==107) { // ------------------- 107  I 4 -2 #107
    spacegroup="I 4 -2";}
  if(spacegroupnumber==108) { // ------------------- 108  I 4 -2c #108
    spacegroup="I 4 -2c";}
  if(spacegroupnumber==109) { // ------------------- 109  I 4bw -2 #109
    spacegroup="I 4bw -2";}
  if(spacegroupnumber==110) { // ------------------- 110  I 4bw -2c #110
    spacegroup="I 4bw -2c";}
  if(spacegroupnumber==111) { // ------------------- 111  P -4 2 #111
    spacegroup="P -4 2";}
  if(spacegroupnumber==112) { // ------------------- 112  P -4 2c #112
    spacegroup="P -4 2c";}
  if(spacegroupnumber==113) { // ------------------- 113  P -4 2ab #113
    spacegroup="P -4 2ab";}
  if(spacegroupnumber==114) { // ------------------- 114  P -4 2n #114
    spacegroup="P -4 2n";}
  if(spacegroupnumber==115) { // ------------------- 115  P -4 -2 #115
    spacegroup="P -4 -2";}
  if(spacegroupnumber==116) { // ------------------- 116  P -4 -2c #116
    spacegroup="P -4 -2c";}
  if(spacegroupnumber==117) { // ------------------- 117  P -4 -2ab #117
    spacegroup="P -4 -2ab";}
  if(spacegroupnumber==118) { // ------------------- 118  P -4 -2n #118
    spacegroup="P -4 -2n";}
  if(spacegroupnumber==119) { // ------------------- 119  I -4 -2 #119
    spacegroup="I -4 -2";}
  if(spacegroupnumber==120) { // ------------------- 120  I -4 -2c #120
    spacegroup="I -4 -2c";}
  if(spacegroupnumber==121) { // ------------------- 121  I -4 2 #121
    spacegroup="I -4 2";}
  if(spacegroupnumber==122) { // ------------------- 122  I -4 2bw #122
    spacegroup="I -4 2bw";}
  if(spacegroupnumber==123) { // ------------------- 123  -P 4 2 #123
    spacegroup="-P 4 2";}
  if(spacegroupnumber==124) { // ------------------- 124  -P 4 2c #124
    spacegroup="-P 4 2c";}
  // setting 2 if(spacegroupnumber==125) { // ------------------- 125  -P 4a 2b #125
  // setting 2   spacegroup="-P 4a 2b";}
  if(spacegroupnumber==125) { // ------------------- 125  -P 4a 2b #125
    spacegroup="P 4 2 -1ab";}
  // setting 2 if(spacegroupnumber==126) { // ------------------- 126  -P 4a 2bc #126
  // setting 2   spacegroup="-P 4a 2bc";}
  if(spacegroupnumber==126) { // ------------------- 126  -P 4a 2bc #126
    spacegroup="P 4 2 -1n";}
  if(spacegroupnumber==127) { // ------------------- 127  -P 4 2ab #127
    spacegroup="-P 4 2ab";}
  if(spacegroupnumber==128) { // ------------------- 128  -P 4 2n #128
    spacegroup="-P 4 2n";}
  // setting 2 if(spacegroupnumber==129) { // ------------------- 129  -P 4a 2a #129
  // setting 2   spacegroup="-P 4a 2a";}
  if(spacegroupnumber==129) { // ------------------- 129  -P 4a 2a #129
    spacegroup="P 4ab 2ab -1ab";}
  // setting 2 if(spacegroupnumber==130) { // ------------------- 130  -P 4a 2ac #130
  // setting 2   spacegroup="-P 4a 2ac";}
  if(spacegroupnumber==130) { // ------------------- 130  -P 4a 2ac #130
    spacegroup="P 4ab 2n -1ab";}
  if(spacegroupnumber==131) { // ------------------- 131  -P 4c 2 #131
    spacegroup="-P 4c 2";}
  if(spacegroupnumber==132) { // ------------------- 132  -P 4c 2c #132
    spacegroup="-P 4c 2c";}
  // setting 2 if(spacegroupnumber==133) { // ------------------- 133  -P 4ac 2b #133
  // setting 2  spacegroup="-P 4ac 2b";}
  if(spacegroupnumber==133) { // ------------------- 133  -P 4ac 2b #133
    spacegroup="P 4n 2c -1n";}
  // setting 2 if(spacegroupnumber==134) { // ------------------- 134  -P 4ac 2bc #134
  // setting 2   spacegroup="-P 4ac 2bc";}
  if(spacegroupnumber==134) { // ------------------- 134  -P 4ac 2bc #134
    spacegroup="P 4n 2 -1n";}
  if(spacegroupnumber==135) { // ------------------- 135  -P 4c 2ab #135
    spacegroup="-P 4c 2ab";}
  if(spacegroupnumber==136) { // ------------------- 136  -P 4n 2n #136
    spacegroup="-P 4n 2n";}
  // setting 2 if(spacegroupnumber==137) { // ------------------- 137  -P 4ac 2a #137
  // setting 2   spacegroup="-P 4ac 2a";}
  if(spacegroupnumber==137) { // ------------------- 137  -P 4ac 2a #137
    spacegroup="P 4n 2n -1n";}
  // setting 2 if(spacegroupnumber==138) { // ------------------- 138  -P 4ac 2ac #138
  // setting 2   spacegroup="-P 4ac 2ac";}
  if(spacegroupnumber==138) { // ------------------- 138  -P 4ac 2ac #138
    spacegroup="P 4n 2ab -1n";}
  if(spacegroupnumber==139) { // ------------------- 139  -I 4 2 #139
    spacegroup="-I 4 2";}
  if(spacegroupnumber==140) { // ------------------- 140  -I 4 2c #140
    spacegroup="-I 4 2c";}
  // setting 2 if(spacegroupnumber==141) { // ------------------- 141  -I 4bd 2 #141
  // setting 2   spacegroup="-I 4bd 2";}
  if(spacegroupnumber==141) { // ------------------- 141  -I 4bd 2 #141
    spacegroup="I 4bw 2bw -1bw";}
  // setting 2 if(spacegroupnumber==142) { // ------------------- 142  -I 4bd 2c #142
  // setting 2   spacegroup="-I 4bd 2c";}
  if(spacegroupnumber==142) { // ------------------- 142  -I 4bd 2c #142
    spacegroup="I 4bw 2aw -1bw";}
  if(spacegroupnumber==143) { // ------------------- 143  P 3 #143
    spacegroup="P 3";}
  if(spacegroupnumber==144) { // ------------------- 144  P 31 #144
    spacegroup="P 31";}
  if(spacegroupnumber==145) { // ------------------- 145  P 32 #145
    spacegroup="P 32";}
  // R setting if(spacegroupnumber==146) { // ------------------- 146  P 3* #146
  // R setting   spacegroup="P 3*";}
  if(spacegroupnumber==146) { // ------------------- 146  R 3 #146
    spacegroup="R 3";}
  if(spacegroupnumber==147) { // ------------------- 147  -P 3 #147
    spacegroup="-P 3";}
  // R setting if(spacegroupnumber==148) { // ------------------- 148  -P 3* #148
  // R setting   spacegroup="-P 3*";}
  if(spacegroupnumber==148) { // ------------------- 148  -R 3 #148
    spacegroup="-R 3";}
  if(spacegroupnumber==149) { // ------------------- 149  P 3 2 #149
    spacegroup="P 3 2";}
  if(spacegroupnumber==150) { // ------------------- 150  P 3 2'' #150
    spacegroup="P 3 2''";}
  if(spacegroupnumber==151) { // ------------------- 151  P 31 2c (0 0 1) #151
    spacegroup="P 31 2c (0 0 1)";}
  if(spacegroupnumber==152) { // ------------------- 152  P 31 2'' #152
    spacegroup="P 31 2''";}
  if(spacegroupnumber==153) { // ------------------- 153  P 32 2c (0 0 -1) #153
    spacegroup="P 32 2c (0 0 -1)";}
  if(spacegroupnumber==154) { // ------------------- 154  P 32 2'' #154
    spacegroup="P 32 2''";}
  if(spacegroupnumber==155) { // ------------------- 155  R 3 2'' #155
    spacegroup="R 3 2''";}
  if(spacegroupnumber==156) { // ------------------- 156  P 3 -2'' #156
    spacegroup="P 3 -2''";}
  if(spacegroupnumber==157) { // ------------------- 157  P 3 -2 #157
    spacegroup="P 3 -2";}
  if(spacegroupnumber==158) { // ------------------- 158  P 3 -2''c #158
    spacegroup="P 3 -2''c";}
  if(spacegroupnumber==159) { // ------------------- 159  P 3 -2c #159
    spacegroup="P 3 -2c";}
  if(spacegroupnumber==160) { // ------------------- 160  R 3 -2'' #160
    spacegroup="R 3 -2''";}
  if(spacegroupnumber==161) { // ------------------- 161  R 3 -2''c #161
    spacegroup="R 3 -2''c";}
  if(spacegroupnumber==162) { // ------------------- 162  -P 3 2 #162
    spacegroup="-P 3 2";}
  if(spacegroupnumber==163) { // ------------------- 163  -P 3 2c #163
    spacegroup="-P 3 2c";}
  if(spacegroupnumber==164) { // ------------------- 164  -P 3 2'' #164
    spacegroup="-P 3 2''";}
  if(spacegroupnumber==165) { // ------------------- 165  -P 3 2''c #165
    spacegroup="-P 3 2''c";}
  if(spacegroupnumber==166) { // ------------------- 166  -R 3 2'' #166
    spacegroup="-R 3 2''";}
  if(spacegroupnumber==167) { // ------------------- 167  -R 3 2''c #167
    spacegroup="-R 3 2''c";}
  if(spacegroupnumber==168) { // ------------------- 168  P 6 #168
    spacegroup="P 6";}
  if(spacegroupnumber==169) { // ------------------- 169  P 61 #169
    spacegroup="P 61";}
  if(spacegroupnumber==170) { // ------------------- 170  P 65 #170
    spacegroup="P 65";}
  if(spacegroupnumber==171) { // ------------------- 171  P 62 #171
    spacegroup="P 62";}
  if(spacegroupnumber==172) { // ------------------- 172  P 64 #172
    spacegroup="P 64";}
  if(spacegroupnumber==173) { // ------------------- 173  P 6c #173
    spacegroup="P 6c";}
  if(spacegroupnumber==174) { // ------------------- 174  P -6 #174
    spacegroup="P -6";}
  if(spacegroupnumber==175) { // ------------------- 175  -P 6 #175
    spacegroup="-P 6";}
  if(spacegroupnumber==176) { // ------------------- 176  -P 6c #176
    spacegroup="-P 6c";}
  if(spacegroupnumber==177) { // ------------------- 177  P 6 2 #177
    spacegroup="P 6 2";}
  if(spacegroupnumber==178) { // ------------------- 178  P 61 2 (0 0 -1) #178
    spacegroup="P 61 2 (0 0 -1)";}
  if(spacegroupnumber==179) { // ------------------- 179  P 65 2 (0 0 1) #179
    spacegroup="P 65 2 (0 0 1)";}
  if(spacegroupnumber==180) { // ------------------- 180  P 62 2c (0 0 1) #180
    spacegroup="P 62 2c (0 0 1)";}
  if(spacegroupnumber==181) { // ------------------- 181  P 64 2c (0 0 -1) #181
    spacegroup="P 64 2c (0 0 -1)";}
  if(spacegroupnumber==182) { // ------------------- 182  P 6c 2c #182
    spacegroup="P 6c 2c";}
  if(spacegroupnumber==183) { // ------------------- 183  P 6 -2 #183
    spacegroup="P 6 -2";}
  if(spacegroupnumber==184) { // ------------------- 184  P 6 -2c #184
    spacegroup="P 6 -2c";}
  if(spacegroupnumber==185) { // ------------------- 185  P 6c -2 #185
    spacegroup="P 6c -2";}
  if(spacegroupnumber==186) { // ------------------- 186  P 6c -2c #186
    spacegroup="P 6c -2c";}
  if(spacegroupnumber==187) { // ------------------- 187  P -6 2 #187
    spacegroup="P -6 2";}
  if(spacegroupnumber==188) { // ------------------- 188  P -6c 2 #188
    spacegroup="P -6c 2";}
  if(spacegroupnumber==189) { // ------------------- 189  P -6 -2 #189
    spacegroup="P -6 -2";}
  if(spacegroupnumber==190) { // ------------------- 190  P -6c -2c #190
    spacegroup="P -6c -2c";}
  if(spacegroupnumber==191) { // ------------------- 191  -P 6 2 #191
    spacegroup="-P 6 2";}
  if(spacegroupnumber==192) { // ------------------- 192  -P 6 2c #192
    spacegroup="-P 6 2c";}
  if(spacegroupnumber==193) { // ------------------- 193  -P 6c 2 #193
    spacegroup="-P 6c 2";}
  if(spacegroupnumber==194) { // ------------------- 194  -P 6c 2c #194
    spacegroup="-P 6c 2c";}
  if(spacegroupnumber==195) { // ------------------- 195  P 2 2 3 #195
    spacegroup="P 2 2 3";}
  if(spacegroupnumber==196) { // ------------------- 196  F 2 2 3 #196
    spacegroup="F 2 2 3";}
  if(spacegroupnumber==197) { // ------------------- 197  I 2 2 3 #197
    spacegroup="I 2 2 3";}
  if(spacegroupnumber==198) { // ------------------- 198  P 2ac 2ab 3 #198
    spacegroup="P 2ac 2ab 3";}
  if(spacegroupnumber==199) { // ------------------- 199  I 2b 2c 3 #199
    spacegroup="I 2b 2c 3";}
  if(spacegroupnumber==200) { // ------------------- 200  -P 2 2 3 #200
    spacegroup="-P 2 2 3";}
  // setting 2 if(spacegroupnumber==201) { // ------------------- 201  -P 2ab 2bc 3 #201
  // setting 2  spacegroup="-P 2ab 2bc 3";}
  if(spacegroupnumber==201) { // ------------------- 201  -P 2ab 2bc 3 #201
    spacegroup="P 2 2 3 -1n";}
  if(spacegroupnumber==202) { // ------------------- 202  -F 2 2 3 #202
    spacegroup="-F 2 2 3";}
  // setting 2 if(spacegroupnumber==203) { // ------------------- 203  -F 2uv 2vw 3 #203
  // setting 2   spacegroup="-F 2uv 2vw 3";}
  if(spacegroupnumber==203) { // ------------------- 203  -F 2uv 2vw 3 #203
    spacegroup="F 2 2 3 -1d";}
  if(spacegroupnumber==204) { // ------------------- 204  -I 2 2 3 #204
    spacegroup="-I 2 2 3";}
  if(spacegroupnumber==205) { // ------------------- 205  -P 2ac 2ab 3 #205
    spacegroup="-P 2ac 2ab 3";}
  if(spacegroupnumber==206) { // ------------------- 206  -I 2b 2c 3 #206
    spacegroup="-I 2b 2c 3";}
  if(spacegroupnumber==207) { // ------------------- 207  P 4 2 3 #207
    spacegroup="P 4 2 3";}
  if(spacegroupnumber==208) { // ------------------- 208  P 4n 2 3 #208
    spacegroup="P 4n 2 3";}
  if(spacegroupnumber==209) { // ------------------- 209  F 4 2 3 #209
    spacegroup="F 4 2 3";}
  if(spacegroupnumber==210) { // ------------------- 210  F 4d 2 3 #210
    spacegroup="F 4d 2 3";}
  if(spacegroupnumber==211) { // ------------------- 211  I 4 2 3 #211
    spacegroup="I 4 2 3";}
  if(spacegroupnumber==212) { // ------------------- 212  P 4acd 2ab 3 #212
    spacegroup="P 4acd 2ab 3";}
  if(spacegroupnumber==213) { // ------------------- 213  P 4bd 2ab 3 #213
    spacegroup="P 4bd 2ab 3";}
  if(spacegroupnumber==214) { // ------------------- 214  I 4bd 2c 3 #214
    spacegroup="I 4bd 2c 3";}
  if(spacegroupnumber==215) { // ------------------- 215  P -4 2 3 #215
    spacegroup="P -4 2 3";}
  if(spacegroupnumber==216) { // ------------------- 216  F -4 2 3 #216
    spacegroup="F -4 2 3";}
  if(spacegroupnumber==217) { // ------------------- 217  I -4 2 3 #217
    spacegroup="I -4 2 3";}
  if(spacegroupnumber==218) { // ------------------- 218  P -4n 2 3 #218
    spacegroup="P -4n 2 3";}
  if(spacegroupnumber==219) { // ------------------- 219  F -4c 2 3 #219
    spacegroup="F -4c 2 3";}
  if(spacegroupnumber==220) { // ------------------- 220  I -4bd 2c 3 #220
    spacegroup="I -4bd 2c 3";}
  if(spacegroupnumber==221) { // ------------------- 221  -P 4 2 3 #221
    spacegroup="-P 4 2 3";}
  // setting 2 if(spacegroupnumber==222) { // ------------------- 222  -P 4a 2bc 3 #222
  // setting 2   spacegroup="-P 4a 2bc 3";}
  if(spacegroupnumber==222) { // ------------------- 222  -P 4a 2bc 3 #222
    spacegroup="P 4 2 3 -1n";}
  if(spacegroupnumber==223) { // ------------------- 223  -P 4n 2 3 #223
    spacegroup="-P 4n 2 3";}
  // setting 2 if(spacegroupnumber==224) { // ------------------- 224  -P 4bc 2bc 3 #224
  // setting 2   spacegroup="-P 4bc 2bc 3";}
  if(spacegroupnumber==224) { // ------------------- 224  -P 4bc 2bc 3 #224
    spacegroup="P 4n 2 3 -1n";}
  if(spacegroupnumber==225) { // ------------------- 225  -F 4 2 3 #225
    spacegroup="-F 4 2 3";}
  if(spacegroupnumber==226) { // ------------------- 226  -F 4c 2 3 #226
    spacegroup="-F 4c 2 3";}
  // setting 2 if(spacegroupnumber==227) { // ------------------- 227  -F 4vw 2vw 3 #227
  // setting 2   spacegroup="-F 4vw 2vw 3";}
  if(spacegroupnumber==227) { // ------------------- 227  -F 4vw 2vw 3 #227
    spacegroup="F 4d 2 3 -1d";}
  // setting 2 if(spacegroupnumber==228) { // ------------------- 228  -F 4cvw 2vw 3 #228
  // setting 2   spacegroup="-F 4cvw 2vw 3";}
  if(spacegroupnumber==228) { // ------------------- 228  -F 4cvw 2vw 3 #228
    spacegroup="F 4d 2 3 -1cd";}
  if(spacegroupnumber==229) { // ------------------- 229  -I 4 2 3 #229
    spacegroup="-I 4 2 3";}
  if(spacegroupnumber==230) { // ------------------- 230  -I 4bd 2c 3 #230
    spacegroup="-I 4bd 2c 3";}
  // done
  return spacegroup;
}

// ***************************************************************************
// GetLaueLabel
// ***************************************************************************
string GetLaueLabel(string& point_group) {
  string laue = "";
  // -1 
  if (point_group=="1" || point_group=="-1"){
    laue = "-1";
  }
  // 2/m
  if (point_group=="2" || point_group=="m" || point_group=="2/m"){
    laue = "2/m";
  }
  // mmm
  if (point_group=="222" || point_group=="mm2" || point_group=="mmm"){
    laue = "mmm";
  }
  // 4/m
  if (point_group=="4" || point_group=="-4" || point_group=="4/m"){
    laue = "4/m";
  }
  // 4/mmmm
  if (point_group=="422" || point_group=="4mm" || point_group=="-42m" || point_group=="-4m2" || point_group=="4/mmm"){
    laue = "4/mmm";
  }
  // -3
  if (point_group=="3" || point_group=="-3"){
    laue = "-3";
  }
  // -3m
  if (point_group=="312" || point_group=="321" || point_group=="32" || point_group=="31m" || point_group=="3m1" || point_group=="-31m" || point_group=="-3m1" || point_group=="-3m" || point_group=="3m"){
    laue = "-3m";
  }
  // 6/m
  if (point_group=="6" || point_group=="-6" || point_group=="6/m"){
    laue = "6/m";
  }
  // 6/mmm
  if (point_group=="622" || point_group=="6mm" || point_group=="-6m2" || point_group=="-62m" || point_group=="6/mmm"){
    laue = "6/mmm";
  }
  // m-3
  if (point_group=="23" || point_group=="m-3"){
    laue = "m-3";
  }
  // m-3m
  if (point_group=="432" || point_group=="-43m" || point_group=="m-3m"){
    laue = "m-3m";
  }
  return laue;
}

// ***************************************************************************
// GetSpaceGroupLabel
// ***************************************************************************
string GetSpaceGroupLabel(int spacegroupnumber) {
  string spacegrouplabel;
  spacegrouplabel="#"+aurostd::StringConvert(spacegroupnumber);
  return spacegrouplabel;
}

// **************************************************************************
// Function MetricTensor
// **************************************************************************
// this function returns the metric tensor
// Corey Oses
xmatrix<double> MetricTensor(const xstructure& a) {return MetricTensor(a.lattice,a.scale);}

xmatrix<double> MetricTensor(const xmatrix<double>& lattice,double scale) {
  if(lattice.rows!=lattice.cols){
    cerr << "metricTensor(): Dimension mismatch, should be square lattice matrix" << endl;
    exit(1);
  }
  xmatrix<double> metric_tensor(lattice.rows,lattice.cols);
  for(int i=1;i<lattice.rows+1;i++){
    for(int j=1;j<lattice.cols+1;j++){
      metric_tensor(i,j)=aurostd::scalar_product(lattice(i),lattice(j));
    }
  }
  return scale*metric_tensor;
}

// **************************************************************************
// Function ReciprocalLattice
// **************************************************************************
// this function returns the reciprocal lattice
// from the original one... it takes the scale !
// Stefano Curtarolo
xmatrix<double> ReciprocalLattice(const xstructure& a){return ReciprocalLattice(a.lattice,a.scale);}

xmatrix<double> ReciprocalLattice(const xmatrix<double>& rlattice,double scale) {        // AFLOW_FUNCTION_IMPLEMENTATION
  xvector<double> a1(3),a2(3),a3(3),b1(3),b2(3),b3(3);
  xmatrix<double> klattice(3,3);  // kvectors are RAWS
  double norm;
  a1[1]=rlattice[1][1];a1[2]=rlattice[1][2];a1[3]=rlattice[1][3];
  a2[1]=rlattice[2][1];a2[2]=rlattice[2][2];a2[3]=rlattice[2][3];
  a3[1]=rlattice[3][1];a3[2]=rlattice[3][2];a3[3]=rlattice[3][3];
  norm=2.0*pi/(det(rlattice)*scale);
  b1=norm*vector_product(a2,a3);
  b2=norm*vector_product(a3,a1);
  b3=norm*vector_product(a1,a2);
  klattice[1][1]=b1[1];klattice[1][2]=b1[2];klattice[1][3]=b1[3];
  klattice[2][1]=b2[1];klattice[2][2]=b2[2];klattice[2][3]=b2[3];
  klattice[3][1]=b3[1];klattice[3][2]=b3[2];klattice[3][3]=b3[3];
  return klattice;
}

//xmatrix<double> ReciprocalLattice(const xmatrix<double>& rlattice) {        // AFLOW_FUNCTION_IMPLEMENTATION
//  return ReciprocalLattice(rlattice,1.0);
//}

// **************************************************************************
// Function KPPRA
// **************************************************************************
// This function calculates k1,k2,k3 starting from real lattice and NK total
// the function does not normalize with number of atoms so the calculation
// must be done somewhere else
string KPPRA(int& k1,int& k2,int& k3,const xmatrix<double>& rlattice,const int& NK) {
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  stringstream aus("");
  aus.precision(5);
  xmatrix<double> klattice(3,3);
  klattice=ReciprocalLattice(rlattice,1.0);
  xvector<double> b1(3),b2(3),b3(3),db1(3),db2(3),db3(3);
  double nb1,nb2,nb3;
  k1=1;k2=1;k3=1;
  int kk1,kk2,kk3,kk;
  b1=klattice(1);nb1=aurostd::modulus(b1);
  b2=klattice(2);nb2=aurostd::modulus(b2);
  b3=klattice(3);nb3=aurostd::modulus(b3);
  if(LDEBUG) aus << "KPPRA LDEBUG:  " << endl << rlattice << endl << endl << klattice << endl << b1 << endl << b2 << endl << b3 << endl << nb1 << endl << nb2 << endl << nb3 << endl;
  if(LDEBUG) aus << "KPPRA LDEBUG:  " << nb1 << " " << nb2 << " " << nb3 << " " << endl;
  if(NK>1) {
    bool found=FALSE;
    double dkdelta,dk;
    dkdelta=0.999;
    dk=aurostd::min(nb1,nb2,nb3);
    kk1=0;kk2=0;kk3=0;
    kk=0;
    int iverbose=0;
    while(!found) {
      kk++;
      if(dk<=1e-5) {
	k1=1;k2=1;k3=1;
	aus << "00000  MESSAGE KPOINTS KPPRA minimun not found k=[" << k1 << "," << k2 << "," << k3 << "]=" << k1*k2*k3 << endl;
      }
      kk1=(int) floor((double) nb1/dk);db1=b1/((double) kk1);
      kk2=(int) floor((double) nb2/dk);db2=b2/((double) kk2);
      kk3=(int) floor((double) nb3/dk);db3=b3/((double) kk3);
      if(kk1+kk2+kk3>iverbose) {
	//  if(!mod(kk,50) || kk1*kk2*kk3>=NK)
	aus << "00000  MESSAGE KPOINTS KPPRA minimizing k=[" << kk1 << "," << kk2 << "," << kk3 << "]=" << kk1*kk2*kk3 << " =[" << aurostd::modulus(db1) << "," << aurostd::modulus(db2) << "," << aurostd::modulus(db3) << "]   dk=" << dk << endl;
	iverbose=kk1+kk2+kk3;
      }
      if(kk1*kk2*kk3>=NK) {
	k1=kk1;k2=kk2;k3=kk3;
	found=TRUE;
      }
      dk=dk*dkdelta;
    }
  } else { // force 1 1 1 for Gamma
    k1=1;k2=1;k3=1;
  }
  db1=b1/((double) k1); db2=b2/((double) k2); db3=b3/((double) k3);
  aus << "00000  MESSAGE KPOINTS KPPRA routine [" << k1 << "," << k2 << "," << k3 << "]=" << k1*k2*k3 << "=[" << aurostd::modulus(db1) << "," << aurostd::modulus(db2) << "," << aurostd::modulus(db3) << "]   "  << endl;
  return aus.str();
}

string KPPRA_LAT(int& k1,int& k2,int& k3,const xmatrix<double>& rlattice,const int& NK) {
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  stringstream aus("");
  aus.precision(5);
  xvector<double> kdata(6),kdatagrid(6);
  xmatrix<double> klattice(3,3);
  klattice=ReciprocalLattice(rlattice,1.0);
  klattice=GetStandardPrimitive(klattice);
  kdata=Getabc_angles(klattice,DEGREES);

  xvector<double> b1(3),b2(3),b3(3),db1(3),db2(3),db3(3);
  double nb1,nb2,nb3;
  k1=1;k2=1;k3=1;
  int kk1,kk2,kk3,kk;
  b1=klattice(1);nb1=aurostd::modulus(b1);
  b2=klattice(2);nb2=aurostd::modulus(b2);
  b3=klattice(3);nb3=aurostd::modulus(b3);
  if(LDEBUG) aus << "KPPRA LDEBUG:  " << endl << rlattice << endl << endl << klattice << endl << b1 << endl << b2 << endl << b3 << endl << nb1 << endl << nb2 << endl << nb3 << endl;
  if(LDEBUG) aus << "KPPRA LDEBUG:  " << nb1 << " " << nb2 << " " << nb3 << " " << endl;
  if(NK>1) {
    bool found=FALSE;
    double dkdelta,dk;
    dkdelta=0.999;
    dk=aurostd::min(nb1,nb2,nb3);
    kk1=0;kk2=0;kk3=0;
    kk=0;
    int iverbose=0;
    while(!found) {
      kk++;
      if(dk<=1e-5) {
	k1=1;k2=1;k3=1;
	aus << "00000  MESSAGE KPOINTS KPPRA minimun not found k=[" << k1 << "," << k2 << "," << k3 << "]=" << k1*k2*k3 << endl;
      }
      kk1=(int) floor((double) nb1/dk);db1=b1/((double) kk1);
      kk2=(int) floor((double) nb2/dk);db2=b2/((double) kk2);
      kk3=(int) floor((double) nb3/dk);db3=b3/((double) kk3);
      kdatagrid=kdata;kdatagrid[1]=1.0;kdatagrid[2]=1.0*aurostd::modulus(db2)/aurostd::modulus(db1);kdatagrid[3]=1.0*aurostd::modulus(db3)/aurostd::modulus(db1);

      if(kk1+kk2+kk3>iverbose) {
	//  if(!mod(kk,50) || kk1*kk2*kk3>=NK)
	aus << "00000  MESSAGE KPOINTS KPPRA minimizing k=[" << kk1 << "," << kk2 << "," << kk3 << "," << GetLatticeType(kdatagrid) << "]=" << kk1*kk2*kk3 << " = [" << aurostd::modulus(db1) << "," << aurostd::modulus(db2) << "," << aurostd::modulus(db3) << "]   dk=" << dk << " " << endl;
	iverbose=kk1+kk2+kk3;
      }
      if(kk1*kk2*kk3>=NK) {
	k1=kk1;k2=kk2;k3=kk3;
	found=TRUE;
      }
      dk=dk*dkdelta;
    }
  } else { // force 1 1 1 for Gamma
    k1=1;k2=1;k3=1;
  }
  db1=b1/((double) k1); db2=b2/((double) k2); db3=b3/((double) k3);
  aus << "00000  MESSAGE KPOINTS KPPRA routine [" << k1 << "," << k2 << "," << k3 << "]=" << k1*k2*k3 << "=[" << aurostd::modulus(db1) << "," << aurostd::modulus(db2) << "," << aurostd::modulus(db3) << "]   "  << endl;
  return aus.str();
}


string KPPRA(xstructure& str,const int& _NK) {
  //  cerr << "KPPRA" << endl;
  int NK=1;
  NK= (int) ((double) _NK/str.atoms.size()+0.5);if(NK<1) NK=1;
  int k1=1,k2=1,k3=1;
  xmatrix<double> rlattice=str.lattice;
  rlattice=str.scale*rlattice;
  // string stringKPPRA=KPPRA_LAT(k1,k2,k3,rlattice,NK);
  string stringKPPRA=KPPRA(k1,k2,k3,rlattice,NK);
  str.kpoints_k1=k1;str.kpoints_k2=k2;str.kpoints_k3=k3;
  str.kpoints_kmax=max(str.kpoints_k1,str.kpoints_k2,str.kpoints_k3);
  str.kpoints_kppra=str.kpoints_k1*str.kpoints_k2*str.kpoints_k3*str.atoms.size();
  return stringKPPRA;
}

// **************************************************************************
// Function KPPRA_DELTA
// **************************************************************************
// This function calculates k1,k2,k3 starting from real lattice and DK
string KPPRA_DELTA(int& k1,int& k2,int& k3,const xmatrix<double>& rlattice,const double& DK) {
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  stringstream aus("");
  aus.precision(5);
  xmatrix<double> klattice(3,3);
  klattice=ReciprocalLattice(rlattice,1.0);
  xvector<double> b1(3),b2(3),b3(3),db1(3),db2(3),db3(3);
  double nb1,nb2,nb3;
  k1=1;k2=1;k3=1;
  int kk1,kk2,kk3,kk;
  b1=klattice(1);nb1=aurostd::modulus(b1);
  b2=klattice(2);nb2=aurostd::modulus(b2);
  b3=klattice(3);nb3=aurostd::modulus(b3);
  if(LDEBUG) aus << "KPPRA LDEBUG:  " << endl << rlattice << endl << endl << klattice << endl << b1 << endl << b2 << endl << b3 << endl << nb1 << endl << nb2 << endl << nb3 << endl;
  if(LDEBUG) aus << "KPPRA LDEBUG:  " << nb1 << " " << nb2 << " " << nb3 << " " << endl;
  if(DK>1.0e-6) {
    bool found=FALSE;
    double dkdelta,dk;
    dkdelta=0.9999;
    dk=aurostd::min(nb1,nb2,nb3);
    kk1=0;kk2=0;kk3=0;
    kk=0;
    int iverbose=0;
    while(!found) {
      kk++;
      kk1=(int) floor((double) nb1/dk);db1=b1/((double) kk1);
      kk2=(int) floor((double) nb2/dk);db2=b2/((double) kk2);
      kk3=(int) floor((double) nb3/dk);db3=b3/((double) kk3);
      if(kk1+kk2+kk3>iverbose) {
        //  if(!mod(kk,50) || kk1*kk2*kk3>=DK)
        aus << "00000  MESSAGE KPOINTS KPPRA minimizing k=[" << kk1 << "," << kk2 << "," << kk3 << "]=" << kk1*kk2*kk3 << " =[" << aurostd::modulus(db1) << "," << aurostd::modulus(db2) << "," << aurostd::modulus(db3) << "]   dk=" << dk << endl;
        iverbose=kk1+kk2+kk3;
      }
      if((aurostd::modulus(db1)<DK) && (aurostd::modulus(db2)<DK) && (aurostd::modulus(db3)<DK)) {
	found=TRUE;k1=kk1;k2=kk2;k3=kk3;}
      dk=dk*dkdelta;
    }
  } else { // force 1 1 1 for Gamma
    k1=1;k2=1;k3=1;
  }
  db1=b1/((double) k1); db2=b2/((double) k2); db3=b3/((double) k3);
  aus << "00000  MESSAGE KPOINTS KPPRA routine [" << k1 << "," << k2 << "," << k3 << "]=" << k1*k2*k3 << "=[" << aurostd::modulus(db1) << "," << aurostd::modulus(db2) << "," << aurostd::modulus(db3) << "]   "  << endl;
  return aus.str();
}

string KPPRA_DELTA(xstructure& str,const double& DK) {
  //  cerr << "KPPRA_DELTA" << endl;
  int k1=1,k2=1,k3=1;
  xmatrix<double> rlattice=str.lattice;
  rlattice=str.scale*rlattice;
  string stringKPPRA=KPPRA_DELTA(k1,k2,k3,rlattice,DK);
  str.kpoints_k1=k1;str.kpoints_k2=k2;str.kpoints_k3=k3;
  str.kpoints_kmax=max(str.kpoints_k1,str.kpoints_k2,str.kpoints_k3);
  str.kpoints_kppra=str.kpoints_k1*str.kpoints_k2*str.kpoints_k3*str.atoms.size();
  return stringKPPRA;
}

// **************************************************************************
// Function GetNBAND
// **************************************************************************
// returns extimated version of NBANDS starting from
// electrons, ions, spin and ispin
int GetNBANDS(int electrons,int nions,int spineach,bool ispin) {
  double out;
  out=max(ceil((electrons+4.0)/1.75)+max(nions/1.75,6.0),ceil(0.80*electrons)); // from VASP
  if(ispin) out=out+(nions*spineach+1)/2;
  //  out=out*1.2;  // safety from vasp
  out=out*1.3;      // safety more
  out=out*1.1+5;    // Thu Jun 11 12:08:42 EDT 2009 // METAL PROJECT
  out=out*1.075;    // Tue Oct 13 07:59:43 EDT 2009 // ICSD PROJECT
  out=out+5;        // Sun Nov  1 10:41:20 EDT 2009 // ICSD PROJECT ORC
  out=out*1.03;     // Tue Feb 26 15:15:36 EST 2013 // HELPS dielectric CALS
  // cerr << "GetNBANDS=" << out << endl;
  out=out*std::pow((double) nions,(double) 0.025);  // rescale so for big numbers of ions you get extra bands // Wed Jun 23 12:29:01 EDT 2010
  //  cerr << "GetNBANDS=" << out << endl;
  // exit(0);
  return (int) ceil(out);
}

// **************************************************************************
// Function GetZVAL from *CAR
// ***************************************************************************
double GetZVAL(const stringstream& sss,vector<double>& vZVAL) {
  xPOTCAR potcar;
  potcar.GetProperties(sss);
  vZVAL.clear();
  for(uint i=0;i<potcar.vZVAL.size();i++) 
    vZVAL.push_back(potcar.vZVAL.at(i));
  return potcar.ZVAL_sum;
  // [OUTCAR] vZVAL.clear();
  // [OUTCAR] string _strline;
  // [OUTCAR] double ZVAL=0.0,_zval_;
  // [OUTCAR] string::size_type idx0;
  // [OUTCAR] bool boolout=FALSE;
  // [OUTCAR] vector<string> vcontent;
  // [OUTCAR] aurostd::string2vectorstring(sss.str(),vcontent);
  // [OUTCAR] for(uint i=0;i<vcontent.size();i++) {
  // [OUTCAR]   if(aurostd::substring2bool(vcontent.at(i),"ZVAL")) 
  // [OUTCAR]     if(aurostd::substring2bool(vcontent.at(i),"POMASS")) {
  // [OUTCAR] 	string aus=vcontent.at(i);
  // [OUTCAR] 	aus=aus.substr(aus.find("ZVAL")+1);
  // [OUTCAR] 	aus=aus.substr(aus.find("=")+1,20);
  // [OUTCAR] 	_zval_=aurostd::string2utype<double>(aus);
  // [OUTCAR] 	ZVAL=ZVAL+_zval_;
  // [OUTCAR] 	vZVAL.push_back(_zval_); }
  // [OUTCAR] }
  // [OUTCAR] return ZVAL;
}

double GetZVAL(const _xvasp& xvasp,vector<double>& vZVAL) {
  return GetZVAL(xvasp.POTCAR,vZVAL);
}

double GetZVAL(const string& directory,vector<double>& vZVAL) {
  stringstream sss("");vector<string> vfile,_aus;
  string file2search="POTCAR,OUTCAR,POTCAR.relax1,POTCAR.relax2,POTCAR.static,POTCAR.bands,OUTCAR.relax1,OUTCAR.relax2,OUTCAR.static,OUTCAR.bands";
  aurostd::string2tokens(file2search,vfile,",");
  _aus=vfile;for(uint i=0;i<_aus.size();i++) {vfile.push_back(_aus.at(i)+".bz2");vfile.push_back(_aus.at(i)+".gz");} // takes care of bz2 and gz
  for(uint i=0;i<vfile.size()&&!sss.str().length();i++)
    if(aurostd::FileExist(directory+"/"+vfile.at(i))) {
      if(aurostd::substring2bool(vfile.at(i),".bz2")) aurostd::bz2file2stringstream(directory+"/"+vfile.at(i),sss);
      if(aurostd::substring2bool(vfile.at(i),".gz")) aurostd::gzfile2stringstream(directory+"/"+vfile.at(i),sss);
      if(!aurostd::substring2bool(vfile.at(i),".bz2") && !aurostd::substring2bool(vfile.at(i),".gz"))
	aurostd::file2stringstream(directory+"/"+vfile.at(i),sss);  
    } //  cerr << sss.str() << endl;
  return GetZVAL(sss,vZVAL);
}

double GetCellAtomZVAL(const stringstream& sss,vector<double>& vZVAL,const stringstream& sstr,vector<double>& sZVAL,string mode) {
  vZVAL.clear();sZVAL.clear();
  GetZVAL(sss,vZVAL);
  stringstream aus(sstr.str()); aus << sstr.str();
  xstructure xstr(aus,IOAFLOW_AUTO);
  if(mode=="CELL" || mode=="") {
    double CellZVAL=xstr.GetZVAL(vZVAL);
    for(uint i=0;i<vZVAL.size();i++) 
      sZVAL.push_back(double(vZVAL.at(i)*double(xstr.num_each_type.at(i))));
    return CellZVAL;
  }
  if(mode=="ATOM") {
    double CellZVAL=xstr.GetZVAL(vZVAL)/double(xstr.atoms.size());
    for(uint i=0;i<vZVAL.size();i++) 
      sZVAL.push_back(double(vZVAL.at(i)*double(xstr.num_each_type.at(i)))/double(xstr.atoms.size()));
    return CellZVAL;
  }
  return 0.0;
}  

double GetCellAtomZVAL(const string& directory,vector<double>& vZVAL,vector<double>& sZVAL,string mode) {  // from directory POT/POS returns total ZVAL cell, vZVAL and sZVAL
  vector<string> vfile,_aus;
  // search for data
  stringstream sss("");
  aurostd::string2tokens("POTCAR,OUTCAR,POTCAR.relax1,POTCAR.relax2,POTCAR.static,POTCAR.bands,OUTCAR.relax1,OUTCAR.relax2,OUTCAR.static,OUTCAR.bands",vfile,",");
  _aus=vfile;for(uint i=0;i<_aus.size();i++) {vfile.push_back(_aus.at(i)+".bz2");vfile.push_back(_aus.at(i)+".gz");} // takes care of bz2 and gz
  for(uint i=0;i<vfile.size()&&!sss.str().length();i++)
    if(aurostd::FileExist(directory+"/"+vfile.at(i))) {
      if(aurostd::substring2bool(vfile.at(i),".bz2")) aurostd::bz2file2stringstream(directory+"/"+vfile.at(i),sss);
      if(aurostd::substring2bool(vfile.at(i),".gz")) aurostd::gzfile2stringstream(directory+"/"+vfile.at(i),sss);
      if(!aurostd::substring2bool(vfile.at(i),".bz2") && !aurostd::substring2bool(vfile.at(i),".gz"))
	aurostd::file2stringstream(directory+"/"+vfile.at(i),sss);  
    } //  cerr << sss.str() << endl;
  // search for xstructure
  stringstream sstr("");
  aurostd::string2tokens("POSCAR,CONTCAR,POSCAR.relax1,POSCAR.relax2,POSCAR.static,POSCAR.bands,CONTCAR.relax1,CONTCAR.relax2,CONTCAR.static,CONTCAR.bands",vfile,",");
  _aus=vfile;for(uint i=0;i<_aus.size();i++) {vfile.push_back(_aus.at(i)+".bz2");vfile.push_back(_aus.at(i)+".gz");} // takes care of bz2 and gz
  for(uint i=0;i<vfile.size()&&!sstr.str().length();i++)
    if(aurostd::FileExist(directory+"/"+vfile.at(i)))  {
      if(aurostd::substring2bool(vfile.at(i),".bz2")) aurostd::bz2file2stringstream(directory+"/"+vfile.at(i),sstr);
      if(aurostd::substring2bool(vfile.at(i),".gz")) aurostd::gzfile2stringstream(directory+"/"+vfile.at(i),sstr);
      if(!aurostd::substring2bool(vfile.at(i),".bz2") && !aurostd::substring2bool(vfile.at(i),".gz"))
	aurostd::file2stringstream(directory+"/"+vfile.at(i),sstr);  
    } // cerr << sstr.str() << endl;
  // done
  return GetCellAtomZVAL(sss,vZVAL,sstr,sZVAL,mode);
}

// ***************************************************************************
// Function GetZVAL
// ***************************************************************************
// Given the ZVAL of each species, it returns total ZVAL of cell
double xstructure::GetZVAL(const vector<double>& vZVAL) {
  if(num_each_type.size()!=vZVAL.size()) {
    cerr << "ERROR GetZVAL num_each_type.size()=" << num_each_type.size() << endl;
    cerr << "ERROR GetZVAL vZVAL.size()=" << vZVAL.size() << endl;
    exit(0);
  }
  double CellZVAL=0.0;
  for(uint i=0;i<vZVAL.size();i++)  
    CellZVAL+=double(vZVAL.at(i)*num_each_type.at(i));
  return CellZVAL;
}

// **************************************************************************
// Function GetPOMASS from *CAR
// ***************************************************************************
double GetPOMASS(const stringstream& sss,vector<double>& vPOMASS) {
  xPOTCAR potcar;
  potcar.GetProperties(sss);
  vPOMASS.clear();
  for(uint i=0;i<potcar.vPOMASS.size();i++) 
    vPOMASS.push_back(potcar.vPOMASS.at(i));
  return potcar.POMASS_sum;
}

double GetPOMASS(const _xvasp& xvasp,vector<double>& vPOMASS) {
  return GetPOMASS(xvasp.POTCAR,vPOMASS);
}

double GetPOMASS(const string& directory,vector<double>& vPOMASS) {
  stringstream sss("");vector<string> vfile,_aus;
  string file2search="POTCAR,OUTCAR,POTCAR.relax1,POTCAR.relax2,POTCAR.static,POTCAR.bands,OUTCAR.relax1,OUTCAR.relax2,OUTCAR.static,OUTCAR.bands";
  aurostd::string2tokens(file2search,vfile,",");
  _aus=vfile;for(uint i=0;i<_aus.size();i++) {vfile.push_back(_aus.at(i)+".bz2");vfile.push_back(_aus.at(i)+".gz");} // takes care of bz2 and gz
  for(uint i=0;i<vfile.size()&&!sss.str().length();i++)
    if(aurostd::FileExist(directory+"/"+vfile.at(i))) {
      if(aurostd::substring2bool(vfile.at(i),".bz2")) aurostd::bz2file2stringstream(directory+"/"+vfile.at(i),sss);
      if(aurostd::substring2bool(vfile.at(i),".gz")) aurostd::gzfile2stringstream(directory+"/"+vfile.at(i),sss);
      if(!aurostd::substring2bool(vfile.at(i),".bz2") && !aurostd::substring2bool(vfile.at(i),".gz"))
	aurostd::file2stringstream(directory+"/"+vfile.at(i),sss);  
    } //  cerr << sss.str() << endl;
  return GetPOMASS(sss,vPOMASS);
}

double GetCellAtomPOMASS(const stringstream& sss,vector<double>& vPOMASS,const stringstream& sstr,vector<double>& sPOMASS,string mode) {
  vPOMASS.clear();sPOMASS.clear();
  GetPOMASS(sss,vPOMASS);
  stringstream aus(sstr.str()); aus << sstr.str();
  xstructure xstr(aus,IOAFLOW_AUTO);
  if(mode=="CELL" || mode=="") {
    double CellPOMASS=xstr.GetPOMASS(vPOMASS);
    for(uint i=0;i<vPOMASS.size();i++) 
      sPOMASS.push_back(double(vPOMASS.at(i)*double(xstr.num_each_type.at(i))));
    return CellPOMASS;
  }
  if(mode=="ATOM") {
    double CellPOMASS=xstr.GetPOMASS(vPOMASS)/double(xstr.atoms.size());
    for(uint i=0;i<vPOMASS.size();i++) 
      sPOMASS.push_back(double(vPOMASS.at(i)*double(xstr.num_each_type.at(i)))/double(xstr.atoms.size()));
    return CellPOMASS;
  }
  return 0.0;
}  

double GetCellAtomPOMASS(const string& directory,vector<double>& vPOMASS,vector<double>& sPOMASS,string mode) {  // from directory POT/POS returns total POMASS cell, vPOMASS and sPOMASS
  vector<string> vfile,_aus;
  // search for data
  stringstream sss("");
  aurostd::string2tokens("POTCAR,OUTCAR,POTCAR.relax1,POTCAR.relax2,POTCAR.static,POTCAR.bands,OUTCAR.relax1,OUTCAR.relax2,OUTCAR.static,OUTCAR.bands",vfile,",");
  _aus=vfile;for(uint i=0;i<_aus.size();i++) {vfile.push_back(_aus.at(i)+".bz2");vfile.push_back(_aus.at(i)+".gz");} // takes care of bz2 and gz
  for(uint i=0;i<vfile.size()&&!sss.str().length();i++)
    if(aurostd::FileExist(directory+"/"+vfile.at(i))) {
      if(aurostd::substring2bool(vfile.at(i),".bz2")) aurostd::bz2file2stringstream(directory+"/"+vfile.at(i),sss);
      if(aurostd::substring2bool(vfile.at(i),".gz")) aurostd::gzfile2stringstream(directory+"/"+vfile.at(i),sss);
      if(!aurostd::substring2bool(vfile.at(i),".bz2") && !aurostd::substring2bool(vfile.at(i),".gz"))
	aurostd::file2stringstream(directory+"/"+vfile.at(i),sss);  
    } //  cerr << sss.str() << endl;
  // search for xstructure
  stringstream sstr("");
  aurostd::string2tokens("POSCAR,CONTCAR,POSCAR.relax1,POSCAR.relax2,POSCAR.static,POSCAR.bands,CONTCAR.relax1,CONTCAR.relax2,CONTCAR.static,CONTCAR.bands",vfile,",");
  _aus=vfile;for(uint i=0;i<_aus.size();i++) {vfile.push_back(_aus.at(i)+".bz2");vfile.push_back(_aus.at(i)+".gz");} // takes care of bz2 and gz
  for(uint i=0;i<vfile.size()&&!sstr.str().length();i++)
    if(aurostd::FileExist(directory+"/"+vfile.at(i)))  {
      if(aurostd::substring2bool(vfile.at(i),".bz2")) aurostd::bz2file2stringstream(directory+"/"+vfile.at(i),sstr);
      if(aurostd::substring2bool(vfile.at(i),".gz")) aurostd::gzfile2stringstream(directory+"/"+vfile.at(i),sstr);
      if(!aurostd::substring2bool(vfile.at(i),".bz2") && !aurostd::substring2bool(vfile.at(i),".gz"))
	aurostd::file2stringstream(directory+"/"+vfile.at(i),sstr);  
    } // cerr << sstr.str() << endl;
  // done
  return GetCellAtomPOMASS(sss,vPOMASS,sstr,sPOMASS,mode);
}

// ***************************************************************************
// Function GetPOMASS
// ***************************************************************************
// Given the POMASS of each species, it returns total POMASS of cell
double xstructure::GetPOMASS(const vector<double>& vPOMASS) {
  if(num_each_type.size()!=vPOMASS.size()) {
    cerr << "ERROR GetPOMASS num_each_type.size()=" << num_each_type.size() << endl;
    cerr << "ERROR GetPOMASS vPOMASS.size()=" << vPOMASS.size() << endl;
    exit(0);
  }
  double CellPOMASS=0.0;
  for(uint i=0;i<vPOMASS.size();i++)  
    CellPOMASS+=double(vPOMASS.at(i)*num_each_type.at(i));
  return CellPOMASS;
}


// **************************************************************************
// Function GetVols, det,
// **************************************************************************
// Wrap up for GetVol functions ... useful to bring
// convasp framework to aflow.

double GetVol(const xmatrix<double>& lat) {
  return abs(det(lat));}        // AFLOW_FUNCTION_IMPLEMENTATION
double det(const xvector<double>& v1,const xvector<double>& v2,const xvector<double>& v3) {
  return scalar_product(v1,vector_product(v2,v3));}
double det(const double& a11,const double& a12,const double& a13,
	   const double& a21,const double& a22,const double& a23,
	   const double& a31,const double& a32,const double& a33) {
  return a11*a22*a33+a12*a23*a31+a13*a21*a32-a13*a22*a31-a12*a21*a33-a11*a23*a32;}
double GetVol(const xvector<double>& v1,const xvector<double>& v2,const xvector<double>& v3) {
  return abs(det(v1,v2,v3));}        // AFLOW_FUNCTION_IMPLEMENTATION


// ***************************************************************************
// Function Getabc_angles
// ***************************************************************************
// This function returns a,b,c,alpha,beta,gamma for the
// cell given the lattice xvectors. Dane Morgan, adjusted by SC

//#define _Getabc_angles Getabc_angles
//#define _Getabc_angles __NO_USE_Sortabc_angles

xvector<double> Getabc_angles(const xmatrix<double>& lat,const int& mode) {  // AFLOW_FUNCTION_IMPLEMENTATION
  xvector<double> data(6);
  data(1)=aurostd::modulus(lat(1));
  data(2)=aurostd::modulus(lat(2));
  data(3)=aurostd::modulus(lat(3));
  data(4)=angle(lat(2),lat(3));
  data(5)=angle(lat(3),lat(1));
  data(6)=angle(lat(1),lat(2));
  if(mode==DEGREES) {
    data(4)*=rad2deg;
    data(5)*=rad2deg;
    data(6)*=rad2deg;
  }
  return data;
}

xvector<double> Getabc_angles(const xmatrix<double>& lat,const xvector<int>& permut,const int& mode) {        // AFLOW_FUNCTION_IMPLEMENTATION
  xvector<double> data(6);
  data(1)=aurostd::modulus(lat(1));
  data(2)=aurostd::modulus(lat(2));
  data(3)=aurostd::modulus(lat(3));
  data(4)=angle(lat(2),lat(3));
  data(5)=angle(lat(3),lat(1));
  data(6)=angle(lat(1),lat(2));
  if(mode==DEGREES) {
    data(4)*=rad2deg;
    data(5)*=rad2deg;
    data(6)*=rad2deg;
  }
  // with permutation - from AVDV
  int i,imin,imax,imid;
  double dmin,dmax;
  dmax=dmin=data(1);
  imax=imin=imid=1;
  for(i=1;i<=3;i++) {
    if(dmax<=data(i)) {dmax=data(i);imax=i; }
    if(dmin>data(i)) {dmin=data(i);imin=i; }}
  // set imid
  for(i=1;i<=3;i++)if(i!=imin && i!=imax) imid=i;
  // if all lattice parameters are equal length, numerical noise may cause imin=imax
  if(imin==imax) for(i=1;i<=3;i++)if(i!=imin && i!=imid) imax=i;
  permut[1]=imin;permut[2]=imid;permut[3]=imax;
  // done permutation - from AVDV
  return data;
}

xvector<double> Getabc_angles(const xvector<double>& r1,      // AFLOW_FUNCTION_IMPLEMENTATION
			      const xvector<double>& r2,      // AFLOW_FUNCTION_IMPLEMENTATION
			      const xvector<double>& r3,      // AFLOW_FUNCTION_IMPLEMENTATION
			      const int& mode) {              // AFLOW_FUNCTION_IMPLEMENTATION
  xmatrix<double> lat(3,3);
  lat(1,1)=r1(1);lat(1,2)=r1(2);lat(1,3)=r1(3);
  lat(2,1)=r2(1);lat(2,2)=r2(2);lat(2,3)=r2(3);
  lat(3,1)=r3(1);lat(3,2)=r3(2);lat(3,3)=r3(3);
  return Getabc_angles(lat,mode);
}

xvector<double> Getabc_angles(const xvector<double>& r1,      // AFLOW_FUNCTION_IMPLEMENTATION
			      const xvector<double>& r2,      // AFLOW_FUNCTION_IMPLEMENTATION
			      const xvector<double>& r3,      // AFLOW_FUNCTION_IMPLEMENTATION
			      const xvector<int>& permut,     // AFLOW_FUNCTION_IMPLEMENTATION
			      const int& mode) {              // AFLOW_FUNCTION_IMPLEMENTATION
  xmatrix<double> lat(3,3);
  lat(1,1)=r1(1);lat(1,2)=r1(2);lat(1,3)=r1(3);
  lat(2,1)=r2(1);lat(2,2)=r2(2);lat(2,3)=r2(3);
  lat(3,1)=r3(1);lat(3,2)=r3(2);lat(3,3)=r3(3);
  return Getabc_angles(lat,permut,mode);
}

xvector<double> Sortabc_angles(const xmatrix<double>& lat,const int& mode) {        // AFLOW_FUNCTION_IMPLEMENTATION
  // with permutation - from AVDV
  int i,imin,imax,imid;
  double dmin,dmax;
  dmax=dmin=aurostd::modulus(lat(1));
  imax=imin=imid=1;
  for(i=1;i<=3;i++) {
    if(dmax<=aurostd::modulus(lat(i))) {dmax=aurostd::modulus(lat(i));imax=i; }
    if(dmin>aurostd::modulus(lat(i))) {dmin=aurostd::modulus(lat(i));imin=i; }}
  // set imid
  for(i=1;i<=3;i++) if(i!=imin && i!=imax) imid=i;
  // if all lattice parameters are equal length, numerical noise may cause imin=imax
  if(imin==imax) for(i=1;i<=3;i++) if(i!=imin && i!=imid) imax=i;

  xvector<double> data(6);
  data(1)=aurostd::modulus(lat(imin));
  data(2)=aurostd::modulus(lat(imid));
  data(3)=aurostd::modulus(lat(imax));
  data(4)=angle(lat(imid),lat(imax));
  data(5)=angle(lat(imin),lat(imax));
  data(6)=angle(lat(imin),lat(imid));
  if(mode==DEGREES) {
    data(4)*=rad2deg;
    data(5)*=rad2deg;
    data(6)*=rad2deg;
  }
  return data;
}

// **************************************************************************
// Function GetClat
// **************************************************************************
// This function returns cartesian lattice xvectors for the cell given
// a,b,c,alpha,beta,gamma. Assumes angles are in radians.
// This routine aligns a along +X, b in XY plane with
// angle gamma from a in +Y direction (bx=b.a_unit=b*cos[gamma],
// by=sqrt(b^2-bx^2)=b*sin[gamma]), and c is then given by
// cx=c.a_unit=c*cos[beta],
// cy=c.by_unit=c*(cos[alpha]-cos[gamma]*cos[beta])/sin[gamma],
// cz=sqrt(c^2-cx^2-cy^2)
// Dane Morgan, adjusted by SC

xmatrix<double> GetClat(const xvector<double>& abc_angles) {   // AFLOW_FUNCTION_IMPLEMENTATION
  xmatrix<double> clattice(3,3);
  double a=abc_angles[1];
  double b=abc_angles[2];
  double c=abc_angles[3];
  double bc= abc_angles[4]*deg2rad; // angle from b to c (remove a)
  double ca= abc_angles[5]*deg2rad; // angle from c to a (remove b)
  double ab= abc_angles[6]*deg2rad; // angle from a to b (remove c)
  //  if(abs(bc)>6.3||abs(ca)>6.3||abs(ab)>6.3) { cerr << _AUROSTD_XLIBS_ERROR_ << "GetClat: angles must be in RADIANT " << endl;exit(0);}
  clattice(1,1)=a;
  clattice(2,1)=b*cos(ab);
  clattice(2,2)=b*sin(ab);
  clattice(3,1)=c*cos(ca);
  if(ab<0.00000001) {
    cerr <<"ERROR: The angle gamma from a to b is too small" << endl;
    cerr <<"ERROR: gamma = " << ab << endl;
    cerr <<"ERROR: STOPPING "<< endl;
    cerr << _AUROSTD_XLIBS_ERROR_ << "ERROR: STOPPING " << endl;
    exit(0);
  }
  clattice(3,2)=c*(cos(bc)-cos(ab)*cos(ca))/sin(ab);
  clattice(3,3)=sqrt(abs(c*c-clattice(3,2)*clattice(3,2)-clattice(3,1)*clattice(3,1)));
  return clattice;
}

xmatrix<double> GetClat(const double &a,const double &b,const double &c,const double &alpha,const double &beta,const double &gamma) {
  xmatrix<double> clattice(3,3);
  double bc= alpha*deg2rad; // angle from b to c (remove a)
  double ca= beta*deg2rad; // angle from c to a (remove b)
  double ab= gamma*deg2rad; // angle from a to b (remove c)
  //  if(abs(bc)>6.3||abs(ca)>6.3||abs(ab)>6.3) { cerr << _AUROSTD_XLIBS_ERROR_ << "GetClat: angles must be in RADIANT " << endl;exit(0);}
  clattice(1,1)=a;
  clattice(2,1)=b*cos(ab);
  clattice(2,2)=b*sin(ab);
  clattice(3,1)=c*cos(ca);
  if(ab<0.00000001) {
    cerr <<"ERROR: The angle gamma from a to b is too small" << endl;
    cerr <<"ERROR: gamma = " << ab << endl;
    cerr <<"ERROR: STOPPING "<< endl;
    cerr << _AUROSTD_XLIBS_ERROR_ << "ERROR: STOPPING " << endl;
    exit(0);
  }
  clattice(3,2)=c*(cos(bc)-cos(ab)*cos(ca))/sin(ab);
  clattice(3,3)=sqrt(abs(c*c-clattice(3,2)*clattice(3,2)-clattice(3,1)*clattice(3,1)));
  return clattice;
}

// **************************************************************************
// Function GetIntpolStr
// **************************************************************************
// This function gets the structure that is linearly interpolated
// a fraction f of the way between strA and strB.  Interpolation
// nincludes lattice parameters.  The scale factors of the interpolated
// structures are all set to 1.

xstructure GetIntpolStr(xstructure strA, xstructure strB, const double& f,const string& path_flag) {
  strA=ReScale(strA,1.0);
  strB=ReScale(strB,1.0);
  // Get new lattice params.
  xmatrix<double> lati(3,3);
  xmatrix<double> latA(3,3);latA=strA.lattice;
  xmatrix<double> latB(3,3);latB=strB.lattice;
  xvector<double> dl(3);
  for(int ic=1;ic<=3;ic++) {
    dl=f*(latB(ic)-latA(ic));
    for(int jc=1;jc<=3;jc++) lati(ic,jc)=latA(ic,jc)+dl(jc);
  }
  // Get new cart. coords.
  if(strA.atoms.size()!=strB.atoms.size()) { 
    cerr << _AUROSTD_XLIBS_ERROR_ << " GetIntpolStr(...) number of atoms must be the same in both structures !!" << endl;
    exit(0);
  }
  int size=strA.atoms.size();
  vector<xvector<double> > cposi(size,3);
  for(int iat=0;iat<size;iat++) {
    xvector<double> dp(3);
    dp=strB.atoms.at(iat).cpos-strA.atoms.at(iat).cpos;
    // If path_flag is n/N then the path is taken between
    // nearest images.  Otherwise path is between the atoms given.
    if(path_flag=="n" || path_flag=="N") {
      xvector<double> ddp(3);ddp=C2F(lati,dp);
      for(int ic=1;ic<=3;ic++) {
        ddp(ic)=ddp(ic)-nint(ddp(ic));
      }
      dp=F2C(lati,ddp);
    }
    dp=f*dp;
    cposi.at(iat)=strA.atoms.at(iat).cpos+dp;
  }
  xstructure stri=strA;
  stri.lattice=lati;
  for(int iat=0;iat<size;iat++) {
    stri.atoms.at(iat).cpos=cposi.at(iat);
    stri.atoms.at(iat).fpos=C2F(stri.lattice,stri.atoms.at(iat).cpos);
  }
  return stri;
}

// **************************************************************************
// Function RadiusSphereLattice
// **************************************************************************
// This function returns the radius of the sphere encompassing the whole cell
double RadiusSphereLattice(const xmatrix<double>& lattice,double scale) {
  double radius=0;
  for(int i=-1;i<=1;i++)
    for(int j=-1;j<=1;j++)
      for(int k=-1;k<=1;k++)
	if(aurostd::modulus((double)i*lattice(1)+(double)j*lattice(2)+(double)k*lattice(3))>radius)
	  radius=aurostd::modulus((double)i*lattice(1)+(double)j*lattice(2)+(double)k*lattice(3));
  return scale*radius;
}

// **************************************************************************
// Function LatticeDimensionSphere
// **************************************************************************
// This function, given a lattice with vectors lat[3][3], finds the
// dimensions along the unit cell vectors such that a sphere of given radius
// fits within a uniform grid of 2dim[1]x2dim[2]x2dim[3] lattice points
// centered at the origin.
//
// The algorithm works by getting the normal (e.g. n1) to each pair of lattice
// vectors (e.g. a2, a3), scaling this normal to have length radius and
// then projecting this normal parallel to the a2,a3 plane onto the
// remaining lattice vector a1. This will tell us the number of a1 vectors
// needed to make a grid to encompass the sphere.
//
// Written as lat_dimension by Anton Van Der Ven
// Adjusted by Stefano Curtarolo

xvector<int> LatticeDimensionSphere(const xmatrix<double>& _lattice, double radius,double scale) {
  // Adapted from AVDV routine
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  xmatrix<double> lattice; lattice=scale*_lattice;
  int i,j,k;
  xmatrix<double> invlattice(3,3),normals(3,3),frac_normals(3,3);
  xvector<int> dim(3);
  double length;
  normals.clear();
  //get the normals to pairs of lattice vectors of length radius
  if(0) { // without eijk
    for(i=1;i<=3;i++)
      for(j=1;j<=3;j++)
	//     normals(i,j)=lattice((i+1)%3+1,(j+1)%3+1)*lattice((i+2)%3+1,(j+2)%3+1)-  // original
	//	lattice((i+1)%3+1,(j+2)%3+1)*lattice((i+2)%3+1,(j+1)%3+1);            // original
	normals(i,j)=lattice((i+1-1)%3+1,(j+1-1)%3+1)*lattice((i+2-1)%3+1,(j+2-1)%3+1)-
	  lattice((i+1-1)%3+1,(j+2-1)%3+1)*lattice((i+2-1)%3+1,(j+1-1)%3+1);
    //  for(i=1;i<=3;i++)
    //   for(j=1;j<=3;j++)
    //    normals(i,j)=vector_product(lattice((i+1-1)%3+1),lattice((i+2-1)%3+1))[j];
  }
  if(1) { // with eijk
    for(int m=1;m<=3;m++)
      for(int n=1;n<=3;n++)
	for(int l=1;l<=3;l++) {
	  normals(1,l)+=aurostd::eijk(l,m,n)*lattice(2,m)*lattice(3,n);
	  normals(2,l)+=aurostd::eijk(l,m,n)*lattice(3,m)*lattice(1,n);
	  normals(3,l)+=aurostd::eijk(l,m,n)*lattice(1,m)*lattice(2,n);
	}
  }
  if(LDEBUG){
    for(uint i=1;i<(uint)normals.rows+1;i++){
      cerr << "LatticeDimensionSphere: normals(" << i << ")=" << normals(i) << endl;
    }
  }
  if(0) { // with eijk and estarijk
    for(int m=1;m<=3;m++)
      for(int n=1;n<=3;n++)
	for(int l=1;l<=3;l++)
	  for(int i=1;i<=3;i++)
	    for(int j=1;j<=3;j++)
	      for(int k=1;k<=3;k++)
		normals(i,l)+=aurostd::eijk(l,m,n)*lattice(j,m)*lattice(k,n)*aurostd::estarijk(i,j,k);
  }
  for(i=1;i<=3;i++) {
    length=aurostd::modulus(normals(i));
    for(j=1;j<=3;j++) normals(i,j)=radius*normals(i,j)/length;
  }
  //get the normals in the coordinates system of the lattice vectors
  invlattice=aurostd::inverse(lattice);

  for(i=1;i<=3;i++) {
    for(j=1;j<=3;j++) {
      frac_normals(i,j)=0.0;
      for(k=1;k<=3;k++)
        frac_normals(i,j)=frac_normals(i,j)+normals(i,k)*invlattice(k,j);
    }
  }
  //the diagonals of frac_normal contain the dimensions of the lattice grid that
  //encompasses a sphere of radius = radius
  for(i=1;i<=3;i++){
    if(LDEBUG){cerr << "LatticeDimensionSphere: abs(frac_normals(i,i))=" << abs(frac_normals(i,i)) << endl;}
    dim(i)=(int)ceil(abs(frac_normals(i,i)));
  }
  if(max(dim)==0) { 
    cerr << _AUROSTD_XLIBS_ERROR_ << "ERROR: LatticeDimensionSphere dim=0 !! " << endl;
    exit(0);
  }
  return dim;
}

xvector<int> LatticeDimensionSphere(const xstructure& str, double radius) {
  return LatticeDimensionSphere(str.lattice,radius,str.scale);//str.scale*str.lattice,radius);  //CO 171024
}

// **************************************************************************
// F2C and C2F exchange transformations
// **************************************************************************
// Stefano Curtarolo
// xf=(CF)*xc and xc=(FC)*xf   (FC)=inv(CF)
//  A is a property so it is a 2 indices tensor
//  x'=A*x and xf'=Af*xf and xc'=Ac*xc
//  xf'=(CF)*xc'=Af*xf=Af*(CF)*xc   mult inv(CF) *
//  xc'=inv(CF)*Af*(CF)*xc HENCE
//  Ac=inv(CF)*Af*(CF);
//  Af=inv(FC)*Ac*(FC)
//  in another way
//  Ac_ij=inv(CF)_il Af_lm (CF)_mj
//  Af_ij=inv(FC)_il Ac_lm (FC)_mj
//  same for operations, although here FC=inv(CF)
//  see notes on operation of force constant tensors

// **************************************************************************
// -------------------------------------------------------- for COLUM vectors
xvector<double> F2C(const double& scale,const xmatrix<double>& lattice,const xvector<double>& fpos) {
  xmatrix<double> f2c(3,3);f2c=scale*trasp(lattice);
  return f2c*fpos;                                 // fpos are F components per COLUMS !
}

xvector<double> F2C(const xmatrix<double>& lattice,const xvector<double>& fpos) {
  return F2C(1.0,lattice,fpos); // VASP cartesian coordinates are intended normalized on the scale
}

xvector<double> C2F(const double& scale,const xmatrix<double>& lattice,const xvector<double>& cpos) {
  xmatrix<double> f2c(3,3);f2c=scale*trasp(lattice);
  return inverse(f2c)*cpos;                         // cpos are C components per COLUMS !
}

xvector<double> C2F(const xmatrix<double>& lattice,const xvector<double>& cpos) {
  return C2F(1.0,lattice,cpos); // VASP cartesian coordinates are intended normalized on the scale
}

// **************************************************************************
// ---------------------------------------------- for a set of COLUM matrices
xmatrix<double> F2C(const double& scale,const xmatrix<double>& lattice,const xmatrix<double>& fpos) {
  xmatrix<double> f2c(3,3);f2c=scale*trasp(lattice);
  return f2c*fpos;                                 // fpos are F components per COLUMS !
}

xmatrix<double> F2C(const xmatrix<double>& lattice,const xmatrix<double>& fpos) {
  return F2C(1.0,lattice,fpos); // VASP cartesian coordinates are intended normalized on the scale
}

xmatrix<double> C2F(const double& scale,const xmatrix<double>& lattice,const xmatrix<double>& cpos) {
  xmatrix<double> f2c(3,3);f2c=scale*trasp(lattice);
  return inverse(f2c)*cpos;                         // cpos are C components per COLUMS !
}

xmatrix<double> C2F(const xmatrix<double>& lattice,const xmatrix<double>& cpos) {
  return C2F(1.0,lattice,cpos); // VASP cartesian coordinates are intended normalized on the scale
}

// **************************************************************************
// -------------------------------------------------------- for _atoms with lattice
_atom F2C(const double& scale,const xmatrix<double>& lattice,const _atom& iatom) {
  xmatrix<double> f2c(3,3);f2c=scale*trasp(lattice);
  _atom oatom(iatom); oatom.cpos=f2c*iatom.fpos; return oatom;
}

_atom F2C(const xmatrix<double>& lattice,const _atom& iatom) {
  return F2C(1.0,lattice,iatom); // VASP cartesian coordinates are intended normalized on the scale
}

_atom C2F(const double& scale,const xmatrix<double>& lattice,const _atom& iatom) {
  xmatrix<double> f2c(3,3);f2c=scale*trasp(lattice);
  _atom oatom(iatom); oatom.fpos=inverse(f2c)*iatom.cpos; return oatom;
}

_atom C2F(const xmatrix<double>& lattice,const _atom& iatom) {
  return C2F(1.0,lattice,iatom); // VASP cartesian coordinates are intended normalized on the scale
}

// **************************************************************************
// -------------------------------------------------------- for _atoms with structure
_atom F2C(const double& scale,const xstructure& str,const _atom& iatom) {
  xmatrix<double> f2c(3,3);f2c=scale*trasp(str.lattice);
  _atom oatom(iatom); oatom.cpos=f2c*iatom.fpos; return oatom;
}

_atom F2C(const xstructure& str,const _atom& iatom) {
  return F2C(1.0,str.lattice,iatom); // VASP cartesian coordinates are intended normalized on the scale
}

_atom C2F(const double& scale,const xstructure& str,const _atom& iatom) {
  xmatrix<double> f2c(3,3);f2c=scale*trasp(str.lattice);
  _atom oatom(iatom); oatom.fpos=inverse(f2c)*iatom.cpos; return oatom;
}

_atom C2F(const xstructure& str,const _atom& iatom) {
  return C2F(1.0,str.lattice,iatom); // VASP cartesian coordinates are intended normalized on the scale
}

// **************************************************************************
// **************************************************************************
// -------------------------------------------------------- for OPERATORS (matrices)
xmatrix<double> FF2CC(const double& scale,const xmatrix<double>& lattice,const xmatrix<double>& fmat) {
  xmatrix<double> f2c(3,3);f2c=scale*trasp(lattice);
  return f2c*fmat*inverse(f2c);                    // fmat is an operation in F coordinates
}

xmatrix<double> FF2CC(const xmatrix<double>& lattice,const xmatrix<double>& fmat) {
  return FF2CC(1.0,lattice,fmat); // VASP cartesian coordinates are intended normalized on the scale
}

xmatrix<double> CC2FF(const double& scale,const xmatrix<double>& lattice,const xmatrix<double>& cmat) {
  xmatrix<double> f2c(3,3);f2c=scale*trasp(lattice);
  return inverse(f2c)*cmat*f2c;                    // cmat is an operation in C coordinates
}

xmatrix<double> CC2FF(const xmatrix<double>& lattice,const xmatrix<double>& cmat) {
  return CC2FF(1.0,lattice,cmat);  // VASP cartesian coordinates are intended normalized on the scale
}

// ***************************************************************************
// Function IdenticalAtoms
// ***************************************************************************
// Makes all atoms the same. Stefano Curtarolo Nov 2008
void xstructure::IdenticalAtoms(void) {
  xvector<double> fpos(3),cpos(3);
  // fix atoms
  for(uint i=0;i<atoms.size();i++) {
    fpos=atoms.at(i).fpos;    // save fpos
    cpos=atoms.at(i).cpos;    // save cpos
    atoms.at(i)=atoms.at(0);  // use this so to copy all info
    atoms.at(i).fpos=fpos;    // fix back fpos
    atoms.at(i).cpos=cpos;    // fix back cpos
    atoms.at(i).partial_occupation_flag=FALSE;
    atoms.at(i).partial_occupation_value=1.0;
  }
  // fix numbers
  num_each_type.clear();num_each_type.push_back(atoms.size());
  comp_each_type.clear();comp_each_type.push_back((double) atoms.size());
  // fix species
  for(uint i=1;i<species.size();i++) {
    species.at(i)="";species_pp.at(i)="";species_volume.at(i)=0.0;species_mass.at(i)=0.0;
  }
}

xstructure IdenticalAtoms(const xstructure& str) {
  xstructure sstr(str);
  sstr.IdenticalAtoms();
  return sstr;
}

// ***************************************************************************
// Function SwapCoordinates
// ***************************************************************************
// Permute Coordinates i with j Stefano Curtarolo Oct 2009
// Wahyu Setyawan: keep the same right-handness of the original
void xstructure::SwapCoordinates(const uint& ii,const uint& jj) {            // Permute Coordinates i with j
  // is obvious
  uint kk=0;
  if((ii==1 && jj==2) || (ii==2 && jj==1)) kk=3;   // permutation
  if((ii==2 && jj==3) || (ii==3 && jj==2)) kk=1;   // permutation
  if((ii==3 && jj==1) || (ii==1 && jj==3)) kk=2;   // permutation
  if((ii==jj) || (ii<1) || (ii>3) || (jj<1) || (jj>3)) return; // nothing to do

  // do the lattice
  double tmp;
  for(uint i=1;i<=3;i++) {
    tmp=lattice[ii][i];lattice[ii][i]=lattice[jj][i];lattice[jj][i]=tmp;   // swap RAW_ii with RAW_jj
    lattice[kk][i]=-lattice[kk][i]; // keep the right-handness of the original
  }
  FixLattices(); // fix the lattices ...
  // cartesian atoms are the same (crystal was not really changed, only definitions of abc were
  // so now I can get the new fractional
  for(uint i=0;i<atoms.size();i++)
    atoms.at(i).fpos=C2F(lattice,atoms.at(i).cpos); // DONE, now it has new fpos with the same cpos !
  // done bring all in cell now
  BringInCell();
}

xstructure SwapCoordinates(const xstructure& str,const uint& ii,const uint& jj) {
  xstructure sstr(str);
  sstr.SwapCoordinates(ii,jj);
  return sstr;
}

// ***************************************************************************
// Function GetLatticeType
// ***************************************************************************
void xstructure::GetLatticeType(void) {
  xstructure str_sp,str_sc;
  GetLatticeType(str_sp,str_sc);
}

void xstructure::GetLatticeType(xstructure& str_sp,xstructure& str_sc) {
  //  bool VERBOSE=TRUE;
  bool VERBOSE=FALSE;
  //DX double eps=0.002,epsang=0.02;
  // double eps=0.02,epsang=0.02;  // JUNKAI
  // DIRECT
  xstructure str_in;//,str_sp,str_sc;
  // start
  str_in=*this;
  str_in.title="NO_RECURSION";
  // str_in.GetPrimitive();
  // str_in.MinkowskiBasisReduction();
  // cerr << str_in << endl;
  // str_in.CalculateSymmetryPointGroup(TRUE);
  // str_in.CalculateSymmetryFactorGroup(TRUE);
  // str_in.CalculateSymmetryPointGroupCrystal(TRUE);
  // DX - START
  bool same_eps = false;
  uint count = 0;
  while(same_eps == false && count++ < 100){
    // DX - END
    if(0) {
      LATTICE::Standard_Lattice_StructureDefault(str_in,str_sp,str_sc); // STD tolerance  // ONLY BRAVAIS_CRYSTAL
      // LATTICE::Standard_Lattice_Structure(str_in,str_sp,str_sc,0.0001,0.001);
    }
    if(1) {
      //    cerr << "LATTICE::Bravais_Lattice_StructureDefault IN" << endl;
      LATTICE::Bravais_Lattice_StructureDefault(str_in,str_sp,str_sc); // STD tolerance  // ONLY BRAVAIS_CRYSTAL
      //  cerr << "LATTICE::Bravais_Lattice_StructureDefault OUT" << endl;
    }
    if(VERBOSE) cerr << "xstructure::GetLatticeType: [4]" << endl;
    if(str_sp.pgroup_calculated==FALSE) str_sp.CalculateSymmetryPointGroup(FALSE);// cerr << "POINT GROUP" << endl;
    if(str_sp.fgroup_calculated==FALSE) str_sp.CalculateSymmetryFactorGroup(FALSE); //cerr << "FACTOR GROUP" << endl;
    if(str_sp.pgroup_xtal_calculated==FALSE) str_sp.CalculateSymmetryPointGroupCrystal(FALSE); //cerr << "POINT GROUP XTAL" << endl;
    //  *this=str_sp; // more obvious but will mess up the structures.... we only want to take the properties
    this->bravais_lattice_type=str_sp.bravais_lattice_type;
    this->bravais_lattice_variation_type=str_sp.bravais_lattice_variation_type;
    this->bravais_lattice_system=str_sp.bravais_lattice_system;
    this->bravais_lattice_lattice_type=str_sp.bravais_lattice_lattice_type;
    this->bravais_lattice_lattice_variation_type=str_sp.bravais_lattice_lattice_variation_type;
    this->bravais_lattice_lattice_system=str_sp.bravais_lattice_lattice_system;
    this->pearson_symbol=str_sp.pearson_symbol;
    this->crystal_family=str_sp.crystal_family;
    this->crystal_system=str_sp.crystal_system;
    this->point_group_crystal_class=str_sp.point_group_crystal_class;
    this->point_group_Shoenflies=str_sp.point_group_Shoenflies;
    this->point_group_Hermann_Mauguin=str_sp.point_group_Hermann_Mauguin;
    this->point_group_orbifold=str_sp.point_group_orbifold;
    this->point_group_type=str_sp.point_group_type;
    this->point_group_order=str_sp.point_group_order;
    this->point_group_structure=str_sp.point_group_structure;
    // RECIPROCAL
    xstructure str_reciprocal_in,str_reciprocal_sp,str_reciprocal_sc;
    str_reciprocal_in.lattice=this->klattice;str_reciprocal_in.FixLattices();
    str_reciprocal_in.title="NO_RECURSION";
    //DX and CO - START
    this->sym_eps=str_reciprocal_in.sym_eps=str_reciprocal_sp.sym_eps=str_reciprocal_sc.sym_eps=str_sp.sym_eps; //DX
    this->sym_eps_calculated=str_reciprocal_in.sym_eps_calculated=str_reciprocal_sp.sym_eps_calculated=str_reciprocal_sc.sym_eps_calculated=str_sp.sym_eps_calculated; //DX
    this->sym_eps_change_count=str_reciprocal_in.sym_eps_change_count=str_reciprocal_sp.sym_eps_change_count=str_reciprocal_sc.sym_eps_change_count=str_sp.sym_eps_change_count; //DX 2/22/18 - added sym_eps change count
    //DX and CO - END
    _atom atom;str_reciprocal_in.AddAtom(atom);
    //LATTICE::Standard_Lattice_Structure(str_reciprocal_in,str_reciprocal_sp,str_reciprocal_sc,eps,epsang); //STEFANO OLD VERSION
    if(VERBOSE) cerr << "xstructure::GetLatticeType: [5]" << endl;
    //DX int ss=0; //JUNKAI
    //DX LATTICE::Standard_Lattice_Structure(str_reciprocal_in,str_reciprocal_sp,str_reciprocal_sc,eps,epsang,ss,_EPS_); //JUNKAI
    //DX 8/14/17 - START - Use real pgroup to calculate pgroupk and then set pgroupk from str_sp to the pgroup and pgroup_xtal of str_reciprocal_in
    //DX 8/14/17 The pgroup and pgroup_xtal are the same for the str_reciprocal structure because there is only one atom at the origin
    //DX 8/14/17 (i.e. lattice and crystal symmetry are the same for the reciprocal space crystal)
    //DX 8/29/17 [OBSOLETE] -since performing full symmetry analysis by default - str_sp.CalculateSymmetryPointGroupKlattice(FALSE);
    str_reciprocal_in.pgroup=str_reciprocal_sp.pgroup=str_reciprocal_sc.pgroup=str_sp.pgroupk;
    str_reciprocal_in.pgroup_calculated=str_reciprocal_sp.pgroup_calculated=str_reciprocal_sc.pgroup_calculated=str_sp.pgroupk_calculated;
    str_reciprocal_in.pgroup_xtal=str_reciprocal_sp.pgroup_xtal=str_reciprocal_sc.pgroup_xtal=str_sp.pgroupk;
    str_reciprocal_in.pgroup_xtal_calculated=str_reciprocal_sp.pgroup_xtal_calculated=str_reciprocal_sc.pgroup_xtal_calculated=str_sp.pgroup_calculated;
    str_reciprocal_in.pgroup_xtal_calculated=str_reciprocal_sp.pgroup_xtal_calculated=str_reciprocal_sc.pgroup_xtal_calculated=str_sp.pgroup_calculated;
    //DX 8/14/17 - END
    LATTICE::Standard_Lattice_StructureDefault(str_reciprocal_in,str_reciprocal_sp,str_reciprocal_sc,false); //DX //DX 2/26/18 - do not need to do full sym for recip
    //DX - START
    if(str_sp.sym_eps == str_reciprocal_sp.sym_eps){
      same_eps = true;
    }
    else{
      str_in.sym_eps = str_sp.sym_eps = str_sc.sym_eps = str_reciprocal_sp.sym_eps;
      str_in.sym_eps_change_count = str_sp.sym_eps_change_count = str_sc.sym_eps_change_count = str_reciprocal_sp.sym_eps_change_count; //DX 2/22/18 - added sym_eps change count
    }
    //DX - END
    this->reciprocal_lattice_type=str_reciprocal_sp.bravais_lattice_type;
    this->reciprocal_lattice_variation_type=str_reciprocal_sp.bravais_lattice_variation_type;
    if(VERBOSE) cerr << "xstructure::GetLatticeType: [6]" << endl;
    // SUPERLATTICE
    xstructure str_superlattice_in,str_superlattice_sp,str_superlattice_sc;
    str_superlattice_in=*this;
    str_superlattice_in.ClearSymmetry();  //DX 8/14/17 - It wasn't cleared, so nothing was being calculated
    str_superlattice_in.title="NO_RECURSION";
    //str_superlattice_in.GetPrimitive(0.01);
    if(VERBOSE) cerr << str_superlattice_in << endl;
    if(VERBOSE) cerr << "xstructure::GetLatticeType: [7]" << endl;
    str_superlattice_in.IdenticalAtoms();  // make superlattice
    if(VERBOSE) cerr << "xstructure::GetLatticeType: [8]" << endl;
    if(VERBOSE) cerr << str_superlattice_in << endl;
    str_superlattice_in.GetPrimitive(0.005);
    if(VERBOSE) cerr << str_superlattice_in << endl;
    if(VERBOSE) cerr << "xstructure::GetLatticeType: [9]" << endl;
    str_superlattice_in.Minkowski_calculated=FALSE;
    if(VERBOSE) cerr << "xstructure::GetLatticeType: [10]" << endl;
    str_superlattice_in.MinkowskiBasisReduction();
    if(VERBOSE) cerr << "xstructure::GetLatticeType: [11]" << endl;
    if(VERBOSE) cerr << str_superlattice_in << endl;
    //  LATTICE::Standard_Lattice_Structure(str_superlattice_in,str_superlattice_sp,str_superlattice_sc); //STEFANO OLD VERSION
    //DX ss=0; //JUNKAI
    //DX and CO - START
    str_superlattice_in.sym_eps=str_superlattice_sp.sym_eps=str_superlattice_sc.sym_eps=str_sp.sym_eps; //DX
    str_superlattice_in.sym_eps_calculated=str_superlattice_sp.sym_eps_calculated=str_superlattice_sc.sym_eps_calculated=str_sp.sym_eps_calculated; //DX
    str_superlattice_in.sym_eps_change_count=str_superlattice_sp.sym_eps_change_count=str_superlattice_sc.sym_eps_change_count=str_sp.sym_eps_change_count; //DX 2/22/18 - added sym_eps change count
    //DX and CO - END
    //DX LATTICE::Standard_Lattice_Structure(str_superlattice_in,str_superlattice_sp,str_superlattice_sc,eps,epsang,ss,_EPS_); //JUNKAI
    LATTICE::Standard_Lattice_StructureDefault(str_superlattice_in,str_superlattice_sp,str_superlattice_sc,false); //DX //DX 2/26/18 - do not need to do full sym for superlattice
    //DX - START
    if(str_sp.sym_eps == str_superlattice_sp.sym_eps){
      same_eps = true;
    }
    else{
      str_sp.sym_eps = str_superlattice_sp.sym_eps;
      str_sp.sym_eps_change_count = str_superlattice_sp.sym_eps_change_count; //DX 2/22/18 - added sym_eps change count
    }
    //DX - END
    if(VERBOSE) cerr << "xstructure::GetLatticeType: [12]" << endl;
    this->bravais_superlattice_type=str_superlattice_sp.bravais_lattice_type;
    this->bravais_superlattice_variation_type=str_superlattice_sp.bravais_lattice_variation_type;
    this->bravais_superlattice_system=str_superlattice_sp.bravais_lattice_system;
    this->pearson_symbol_superlattice=str_superlattice_sp.pearson_symbol;
    if(count==100){
      cerr << "ERROR in Bravais_Lattice_StructureDefault(): Unable to find reliable sym_eps." << endl;
      break;
    }
  } //DX while loop
}

string GetLatticeType(xmatrix<double> lattice) {
  //DX double eps=0.00010,epsang=0.01;
  xstructure str_in,str_sp,str_sc;
  str_in.lattice=lattice;str_in.FixLattices();
  str_in.title="NO_RECURSION";
  _atom atom;str_in.AddAtom(atom);
  // LATTICE::Standard_Lattice_Structure(str_in,str_sp,str_sc,eps,epsang); // STEFANO OLD VERSION
  //DX int ss=0; //JUNKAI
  //DX LATTICE::Standard_Lattice_Structure(str_in,str_sp,str_sc,eps,epsang,ss,_EPS_); //JUNKAI
  LATTICE::Standard_Lattice_StructureDefault(str_in,str_sp,str_sc); //DX
  return str_sp.bravais_lattice_type;
  //  return str_sp.bravais_lattice_variation_type;
}

string GetLatticeType(xvector<double> data) {
  xmatrix<double> lattice(GetClat(data));
  return GetLatticeType(lattice);
}

// ***************************************************************************
// Function Standard_Primitive
// ***************************************************************************
// Lattice Reduction to Max Orthogonality (MINK) and then Niggly Form
xstructure Standard_Primitive_UnitCellForm(const xstructure& a) {
  xstructure str_in(a),str_sp,str_sc;
  if(str_in.Standard_Lattice_avoid==TRUE) return str_in;     // Nothing to do
  if(str_in.Standard_Lattice_primitive==TRUE) return str_in; // already Primitive
  LATTICE::Standard_Lattice_StructureDefault(str_in,str_sp,str_sc);
  return str_sp;
}

xstructure GetStandardPrimitive(const xstructure& a) {
  // cerr << "GetStandardPrimitive(const xstructure& a)" << endl;
  return Standard_Primitive_UnitCellForm(a);
}

void xstructure::Standard_Primitive_UnitCellForm(void) {
  xstructure str_sp,str_sc;
  if(Standard_Lattice_avoid==TRUE) return;     // Nothing to do
  if(Standard_Lattice_primitive==TRUE) return; // already Primitive
  // cerr << bravais_lattice_type << endl;
  LATTICE::Standard_Lattice_StructureDefault(*this,str_sp,str_sc);
  *this=str_sp;
}

void xstructure::GetStandardPrimitive(void) {
  Standard_Primitive_UnitCellForm();
}

xmatrix<double> GetStandardPrimitive(xmatrix<double> lattice) {
  xstructure str;
  str.lattice=lattice;str.scale=1.0;str.FixLattices();
  str.title="NO_RECURSION";
  _atom atom;str.AddAtom(atom);
  str=GetStandardPrimitive(str);
  return str.lattice;
}

xvector<double> GetStandardPrimitive(xvector<double> data) {
  xmatrix<double> lattice(GetClat(data));
  return Getabc_angles(GetStandardPrimitive(lattice),DEGREES);
}

// ***************************************************************************
// Function Standard_Conventional
// ***************************************************************************
// Lattice Reduction to Max Orthogonality (MINK) and then Niggly Form
xstructure Standard_Conventional_UnitCellForm(const xstructure& a) {
  xstructure str_in(a),str_sp,str_sc;
  if(str_in.Standard_Lattice_avoid==TRUE) return str_in;        // Nothing to do
  if(str_in.Standard_Lattice_conventional==TRUE) return str_in; // already Conventional
  LATTICE::Standard_Lattice_StructureDefault(str_in,str_sp,str_sc);
  return str_sc;
}

xstructure GetStandardConventional(const xstructure& a) {
  return Standard_Conventional_UnitCellForm(a);
}

void xstructure::Standard_Conventional_UnitCellForm(void) {
  xstructure str_sp,str_sc;
  if(Standard_Lattice_avoid==TRUE) return;        // Nothing to do
  if(Standard_Lattice_conventional==TRUE) return; // already Conventional
  LATTICE::Standard_Lattice_StructureDefault(*this,str_sp,str_sc);
  *this=str_sc;
}

void xstructure::GetStandardConventional(void) {
  Standard_Conventional_UnitCellForm();
}

xmatrix<double> GetStandardConventional(xmatrix<double> lattice) {
  xstructure str;
  str.lattice=lattice;str.scale=1.0;str.FixLattices();
  str.title="NO_RECURSION";
  _atom atom;str.AddAtom(atom);
  str=GetStandardConventional(str);
  return str.lattice;
}

xvector<double> GetStandardConventional(xvector<double> data) {
  xmatrix<double> lattice(GetClat(data));
  return Getabc_angles(GetStandardConventional(lattice),DEGREES);
}

// ***************************************************************************
// Function SpeciesLabel
// ***************************************************************************
// Returns the name of the specie (the name of the 1st atom of the specie
// Stefano Curtarolo - nov 2008
string xstructure::SpeciesLabel(const uint& A) {
  if(A>num_each_type.size()) return "xxx (outside boundaries)"; // outside boundaries
  string name;vector<string> tokens;

  uint i=0,A_start=0;//,A_stop=0;
  for(i=0;i<num_each_type.size();i++) {
    if(i<A)  {A_start+=num_each_type.at(i);}
    // [UNNECESSARY] if(i==A) {A_start=A_start;}
    //    if(i==A) {A_stop=A_start+num_each_type.at(i);}
  }
  if(atoms.at(A_start).name_is_given==FALSE) return "xxx (name not given)"; // name not given
  name=atoms.at(A_start).name;
  aurostd::StringSubst(name,"+"," ");aurostd::StringSubst(name,"-"," ");
  aurostd::string2tokens(name,tokens," ");
  name=tokens[0];
  return name;
}

string SpeciesLabel(const xstructure& str,const uint& A) {
  xstructure sstr(str);
  return sstr.SpeciesLabel(A);
}

// ***************************************************************************
// Function SpeciesSwap
// ***************************************************************************
// Permute Species A with B (safe for species C) Stefano Curtarolo Nov 2008
void xstructure::SpeciesSwap(const uint& specieA,const uint& specieB) {
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  // some useful checks
  if(num_each_type.size()==0) return; // empty structures have nothing to swap
  if(num_each_type.size()==1) return; // pure structures have nothing to swap
  // check done
  deque<_atom> atoms_buf;
  uint i,j,k;
  uint specieAA=0,specieAA_start=0,specieAA_stop=0;
  uint specieBB=0,specieBB_start=0,specieBB_stop=0;
  // PUT ALPHABETIC
  if(specieA<specieB) {specieAA=specieA;specieBB=specieB;}
  if(specieA>specieB) {specieAA=specieB;specieBB=specieA;}
  // CHECK
  if(specieA==specieB) return;                    // NO SWAP
  if(specieA>num_each_type.size()) return;  // NO SWAP nothing to swap
  if(specieB>num_each_type.size()) return;  // NO SWAP nothing to swap
  // FIND BOUNDARIES
  for(i=0;i<num_each_type.size();i++) {
    if(i<specieAA)  {specieAA_start+=num_each_type.at(i);}
    if(i==specieAA) {specieAA_stop=specieAA_start+num_each_type.at(i);}
    if(i<specieBB)  {specieBB_start+=num_each_type.at(i);}
    if(i==specieBB) {specieBB_stop=specieBB_start+num_each_type.at(i);}
  }
  if(LDEBUG) cerr << "DEBUG  atoms.size()=" << atoms.size() << endl;
  if(LDEBUG) cerr << "DEBUG  specieAA=" << specieAA << endl;
  if(LDEBUG) cerr << "DEBUG  specieBB=" << specieBB << endl;
  if(LDEBUG) cerr << "DEBUG  specieAA_start=" << specieAA_start << endl;
  if(LDEBUG) cerr << "DEBUG  specieAA_stop="  << specieAA_stop << endl;
  if(LDEBUG) cerr << "DEBUG  specieBB_start=" << specieBB_start << endl;
  if(LDEBUG) cerr << "DEBUG  specieBB_stop="  << specieBB_stop << endl;
  // ATOMS  -----------------------------
  atoms_buf.clear();
  for(i=0;i<atoms.size();i++) atoms_buf.push_back(atoms.at(i)); // create buffer
  atoms.clear();   // RECONSTRUCTING
  for(i=0;i<specieAA_start;i++) atoms.push_back(atoms_buf.at(i));                // before specieA preserve
  for(i=specieBB_start;i<specieBB_stop;i++) atoms.push_back(atoms_buf.at(i));    // during specieA swap with specieB
  for(i=specieAA_stop;i<specieBB_start;i++) atoms.push_back(atoms_buf.at(i));    // between specieA and specieB preserve
  for(i=specieAA_start;i<specieAA_stop;i++) atoms.push_back(atoms_buf.at(i));    // during specieB swap with specieA
  for(i=specieBB_stop;i<atoms_buf.size();i++) atoms.push_back(atoms_buf.at(i));  // after specieB preserve
  // swap numbers
  uint iaus;
  iaus=num_each_type.at(specieAA);num_each_type.at(specieAA)=num_each_type.at(specieBB);num_each_type.at(specieBB)=iaus;
  // now fix the types
  k=0;
  for(i=0;i<num_each_type.size();i++)
    for(j=0;j<(uint) num_each_type.at(i);j++)
      atoms.at(k++).type=i; // done
  // swap species
  string saus;double daus;
  if(specieAA<species.size() && specieBB<species.size()) {saus=species.at(specieAA);species.at(specieAA)=species.at(specieBB);species.at(specieBB)=saus;}
  if(specieAA<species_pp.size() && specieBB<species_pp.size()) {saus=species_pp.at(specieAA);species_pp.at(specieAA)=species_pp.at(specieBB);species_pp.at(specieBB)=saus;}
  if(specieAA<species_volume.size() && specieBB<species_volume.size()) {daus=species_volume.at(specieAA);species_volume.at(specieAA)=species_volume.at(specieBB);species_volume.at(specieBB)=daus;}
  if(specieAA<species_mass.size() && specieBB<species_mass.size()) {daus=species_mass.at(specieAA);species_mass.at(specieAA)=species_mass.at(specieBB);species_mass.at(specieBB)=daus;}
}

// ***************************************************************************
// Function SpeciesGetAlphabetic
// ***************************************************************************
// Tell if two species are alphabetic!  Stefano Curtarolo Nov 2008
bool xstructure::SpeciesGetAlphabetic(void) {
  if(num_each_type.size()!=species.size()) {
    cerr << "ERROR: SpeciesGetAlphabetic: num_each_type.size()!=species.size()   ("<<num_each_type.size()<<","<<species.size()<<")" << endl;
    cerr << "ERROR: num_each_type.size()="<<num_each_type.size()<<": "; for(uint i=0;i<num_each_type.size();i++) cerr << num_each_type.at(i) << " "; cerr << endl;
    cerr << "ERROR: species.size()="<<species.size()<< ": "; for(uint i=0;i<species.size();i++) cerr << species.at(i) << " "; cerr << endl;
    cerr << "ERROR: species_pp.size()="<<species_pp.size()<< ": "; for(uint i=0;i<species_pp.size();i++) cerr << species_pp.at(i) << " "; cerr << endl;
    cerr << "ERROR: species_volume.size()="<<species_volume.size()<< ": "; for(uint i=0;i<species_volume.size();i++) cerr << species_volume.at(i) << " "; cerr << endl;
    cerr << "ERROR: species_mass.size()="<<species_mass.size()<< ": "; for(uint i=0;i<species_mass.size();i++) cerr << species_mass.at(i) << " "; cerr << endl;
    exit(0);
  }
  // some useful checks
  if(species.size()==0) return TRUE; // empty structures are always alphabetic
  if(species.size()==1) return TRUE; // pure structures are always alphabetic
  // check done
  vector<string> sspecies;
  for(uint isp=0;isp<species.size();isp++)
    sspecies.push_back(aurostd::RemoveNumbers(KBIN::VASP_PseudoPotential_CleanName(species.at(isp))));
  
  for(uint isp=0;isp<sspecies.size()-1;isp++)
    if(sspecies.at(isp)>sspecies.at(isp+1)) return FALSE;
  // otherwise return TRUE;
  return TRUE;
}

// ***************************************************************************
// Function SpeciesPutAlphabetic
// ***************************************************************************
// Tell if two species are alphabetic!  Stefano Curtarolo Nov 2008
bool xstructure::SpeciesPutAlphabetic(void) {
  if(num_each_type.size()!=species.size()) {
    cerr << "ERROR: SpeciesPutAlphabetic: num_each_type.size()!=species.size()   ("<<num_each_type.size()<<","<<species.size()<<")" << endl;
    exit(0);
  }
  // some useful checks
  if(species.size()==0) return TRUE; // empty structures are always alphabetic
  if(species.size()==1) return TRUE; // pure structures are always alphabetic
  // check done
  vector<string> sspecies;

  while(SpeciesGetAlphabetic()==FALSE) {
    sspecies.clear();
    for(uint isp=0;isp<species.size();isp++)
      sspecies.push_back(aurostd::RemoveNumbers(KBIN::VASP_PseudoPotential_CleanName(species.at(isp))));
    for(uint isp=0;isp<sspecies.size()-1;isp++)
      if(sspecies.at(isp)>sspecies.at(isp+1)) SpeciesSwap(isp,isp+1);
  }
  MakeBasis(); // repetita iuvant
  // otherwise return TRUE;
  return TRUE;
}

// ***************************************************************************
// Function SpeciesString
// ***************************************************************************
// Returns a string with the species.  Stefano Curtarolo Nov 2008
string xstructure::SpeciesString(void) {
  stringstream strstream;
  strstream.clear();strstream.str(std::string());
  for(uint i=0;i<num_each_type.size();i++) {
    strstream << SpeciesLabel(i);
    if(i!=num_each_type.size()-1)
      strstream << " ";
  }
  return strstream.str();
}

// ***************************************************************************
// Function SetSpecies
// ***************************************************************************
// Set the species  Stefano Curtarolo Nov 2014
uint xstructure::SetSpecies(deque<string> vspecies) {
  if(vspecies.size()!=species.size() ) {cerr << "ERROR - xstructure::SetSpecies:  vspecies.size()!=species.size()" << endl;exit(0);}
  if(vspecies.size()!=num_each_type.size() ) {cerr << "ERROR - xstructure::SetSpecies:  vspecies.size()!=num_each_type.size()" << endl;exit(0);}
  int iatom=0;
  for(uint itype=0;itype<num_each_type.size();itype++) {
    //      string species=string(argv.at(2+b.atoms.at(iatom).type));
    string itype_species=vspecies.at(itype);
    species.at(itype)=itype_species;
    for(int j=0;j<num_each_type.at(itype);j++) {
      atoms.at(iatom).name=itype_species;    // CONVASP_MODE
      atoms.at(iatom).CleanName();
      //DX 9/21/17 - Need to keep spin info atoms.at(iatom).CleanSpin();
      atoms.at(iatom).name_is_given=TRUE;
      iatom++;
    }
  }
  return vspecies.size();
}

// ***************************************************************************
// NIGGLI NIGGLI NIGGLI NIGGLI NIGGLI NIGGLI NIGGLI NIGGLI NIGGLI NIGGLI NIGGL
// ***************************************************************************
// GetNiggliCell
// ***************************************************************************
// Calculates the reduced Niggli cell.  It is based on
// a program of Eric Wu's.  Here is his + my documentation.
/*
  FUNCTION:
  Subroutine calculates the reduced cell (Niggli cell) for a given
  primitive cell.  It also gives the transformation matrix to go from
  the primitive cell to the Niggli cell.

  REFERENCES:
  Y. Lepage, J. Appl. Cryst 15 255-229 (1982) and
  Y. Lepage, J. Appl. Cryst 20 264-269 (1987)
  International tables for crystallography Ch 5,9
  I Krivy, B. Gruber. Acta Cryst 1976 A32 297 (1975)
  W. Clegg., Acta Cryst A37 913-915 (1981)
  A.D. Mighell J. Appl Cryst 9 491-498 (1976)

  METHODOLOGY:
  The Niggli cell is UNIQUE and reduced. A Niggli cell is a Buerger cell
  (shortest lattice vectors) but  not necessarily the other way around,
  since the Buerger cell is not unique.  This info is useful because with
  the Niggli cell one knows instantly the type of bravais lattice
  and how to go to the conventional cell if you apply the correct
  algorithm.  There are 2 ways to do this.  The first is the straight
  way -  see International Tables Ch. 9 .  The other involves counting
  the # of 2fold axes and going from there, see the Lepage and Clegg
  references.  This second method is more useful if you want to find
  more than just bravais lattice (symmetry elements or space group
  for example).  The subroutine assumes you input a primitive cell
  and it spits out the Niggli cell and the transformation matrix
  to go from the primitive cell to Niggli cell.
  NB! There is a typo in step 3 of Krivy and Gruber that I corrected.
  [ (ksi*eta*ksi) should be (ksi*eta*zeta) in first part of Eq. 3. ]

  INPUT:
  in_lat: a primitive set of lattice vectors,
  in_lat(1,1)=a_x
  in_lat(1,2)=a_y
  in_lat(2,1)=b_x  ...

  OUTPUT:
  niggli_lat - reduced (Niggli) cell, same format
  P - Matrix to go from Primitive to Niggli cell
  Q - P^-1;

  For more info of the definitions here see Int. Tables Cryst. ch.5
  Let nlat be the niggli lattice, where the first column is a,
  second b, third c.  Let lat be the original cell equivalent.
  The P matrix to transform cells is generally defined by
  nlat = lat * P
  We will follow this convention for defining P and Q.
  Note that is the rest of the code, we define lattices with
  rows equal to a, b, and c.  Therefore, under our conventional
  definitions of nlat, call it nlat_c=nlat^T, and lat, call
  it lat_c=lat^T, we have
  nlat_c = P^T * lat_c
  where ^T means the transpose.
*/

void _sdebug_GetNiggliCell(int,int,double,double,double,double,double,double,ostream&);

bool GetNiggliCell(const xmatrix<double>& in_lat,xmatrix<double>& niggli_lat,xmatrix<double>& P,xmatrix<double>& Q) {
  return GetNiggliCell_20180213(in_lat,niggli_lat,P,Q);
}

//DX 2/13/18 - Niggli algorithm, fixed tolerances and added goto loop in step 6 - START
bool GetNiggliCell_20180213(const xmatrix<double>& in_lat,xmatrix<double>& niggli_lat,xmatrix<double>& P,xmatrix<double>& Q) {
  // return FALSE if failed
  // rerurn TRUE if ok

  //  double RENORM=1E+6;  // DANE not used
  int MAXITER=100000;
  double TOL=_ZERO_TOL_; //DX 2/12/18 - from 1e-15 to 1e-10

  // Initialize matrices for tranformations (3x3).
  xmatrix<double> m1(3,3);
  m1(1,2)=1.0; m1(2,1)=1.0; m1(3,3)=-1.0;
  xmatrix<double> m2(3,3);
  m2(1,1)=-1.0; m2(2,3)=1.0; m2(3,2)=1.0;
  xmatrix<double> m3(3,3);
  // Els 1,1 2,2 3,3 are changed later for m3
  xmatrix<double> m4(3,3);
  // Els 1,1 2,2 3,3 are changed later for m4
  xmatrix<double> m5(3,3);
  m5(1,1)=1.0; m5(2,2)=1.0; m5(3,3)=1.0;
  // El 3,2 is changed later for m5
  xmatrix<double> m6(3,3);
  m6(1,1)=1.0; m6(2,2)=1.0; m6(3,3)=1.0;
  // El 3,1 is changed later for m6
  xmatrix<double> m7(3,3);
  m7(1,1)=1.0; m7(2,2)=1.0; m7(3,3)=1.0;
  // El 2,1 is changed later for m7
  xmatrix<double> m8(3,3);
  m8(1,1)=1.0; m8(2,2)=1.0; m8(1,3)=1.0; m8(2,3)=1.0; m8(3,3)=1.0;

  // Initialize a, b, c, ksi, eta, zeta
  xvector<double> indat(6);
  indat=Getabc_angles(in_lat,RADIANTS);
  double a=indat(1)*indat(1);
  double b=indat(2)*indat(2);
  double c=indat(3)*indat(3);
  double ksi=2.0*indat(2)*indat(3)*cos(indat(4));
  double eta=2.0*indat(1)*indat(3)*cos(indat(5));
  double zeta=2.0*indat(1)*indat(2)*cos(indat(6));
  // Round to RENORM decimal place to eliminate numerical errors
  //  double a=Nint(indat(1)*indat(1)*RENORM);
  //  double b=Nint(indat(2)*indat(2)*RENORM);
  //  double c=Nint(indat(3)*indat(3)*RENORM);
  //  double ksi=Nint(2.0*indat(2)*indat(3)*cos(indat(4))*RENORM);
  //  double eta=Nint(2.0*indat(1)*indat(3)*cos(indat(5))*RENORM);
  //  double zeta=Nint(2.0*indat(1)*indat(2)*cos(indat(6))*RENORM);

  // Dummy variables
  double temp,temp1,temp2,temp3;

  // Initialize tranformation matrix
  //  xmatrix<double> P(3,3);
  P(1,1)=1; P(2,2)=1; P(3,3)=1;

  int cnt=0;
  // tpx
  // _sdebug_GetNiggliCell(0,0,a,b,c,ksi,eta,zeta,cout);

  // Start loop

 LoopHead:

  cnt++;
  //  cerr << "DEBUG: Niggli() cnt=" << cnt << endl;
  if(cnt>MAXITER) {
    //     stringstream oss;
    //     // -------
    //     oss << "EEEEE  aflow " << VERSION << endl;
    //     oss << "EEEEE  ERROR: CellReduceFuncs/GetNiggliCell" << endl;
    //     oss << "EEEEE  ERROR: Number of interations greater the MAXITER = " << MAXITER << endl;
    //     oss << "EEEEE  ERROR: This seems like too many - there is probably some problem." << endl;
    //     oss << "EEEEE  ERROR:." << endl;
    //     //  oss << endl;
    //     // -------
    //     cerr << oss.str();
    //     //  cout << oss.str();
    //     // -------
    return FALSE;
  }

  
  // Step 1
  // DX 2/9/18 [OBSOLETE] if(((a-b)>TOL) || ((abs(a-b)<TOL) && (abs(ksi)>abs(eta)) ) ) {
  if(((a-b)>TOL) || ((abs(a-b)<TOL) && ((abs(ksi)-abs(eta))>TOL) ) ) { //DX 2/9/18 - more robust; precision
    temp=a;a=b;b=temp;
    temp=-ksi;ksi=-eta;eta=temp;
    P=P*m1;
    // tpx
    // _sdebug_GetNiggliCell(1,cnt,a,b,c,ksi,eta,zeta,cout);
  }

  // Step 2
  // DX 2/9/18 [OBSOLETE] if(((b-c)>TOL) || ((abs(b-c)<TOL) && (abs(eta)>abs(zeta)) ) ) {
  if(((b-c)>TOL) || ((abs(b-c)<TOL) && ((abs(eta)-abs(zeta))>TOL) ) ) { //DX 2/9/18 - more robust; precision
    temp=b;b=c;c=temp;
    temp=-eta;eta=-zeta;zeta=temp;
    P=P*m2;
    // tpx
    // _sdebug_GetNiggliCell(2,cnt,a,b,c,ksi,eta,zeta,cout);
    goto LoopHead;
  }

  // Step 3
  // DX 2/9/18 [OBSOLETE] if((ksi*eta*zeta)>0.0) {
  if((ksi*eta*zeta)>TOL) { //DX 2/9/18 - more robust; precision
    m3(1,1)=SignNoZero(ksi);
    m3(2,2)=SignNoZero(eta);
    m3(3,3)=SignNoZero(zeta);
    P=P*m3;
    ksi=abs(ksi);
    eta=abs(eta);
    zeta=abs(zeta);
    // tpx
    // _sdebug_GetNiggliCell(3,cnt,a,b,c,ksi,eta,zeta,cout);
  }

  // Step 4
  // Eric's original if((abs(ksi)<TOL)||(abs(eta)<TOL)||(abs(zeta)<TOL)||(ksi*eta*zeta<=0)) {
  if((abs(ksi)<TOL)||(abs(eta)<TOL)||(abs(zeta)<TOL)||(ksi*eta*zeta<=-TOL)) {  //DX 2/12/18 - if any are zero, then, ksi*eta*zeta is zero
  //if(ksi*eta*zeta<=TOL) { //DX 2/9/18 - more robust; precision
    m4(1,1)=-(SignNoZero(ksi));
    m4(2,2)=-(SignNoZero(eta));
    m4(3,3)=-(SignNoZero(zeta));
    if(abs(ksi)<TOL) m4(1,1)=m4(2,2)*m4(3,3);
    if(abs(eta)<TOL) m4(2,2)=m4(1,1)*m4(3,3);
    if(abs(zeta)<TOL) m4(3,3)=m4(1,1)*m4(2,2);
    P=P*m4;
    ksi=-abs(ksi);
    eta=-abs(eta);
    zeta=-abs(zeta);
    // tpx
    // _sdebug_GetNiggliCell(4,cnt,a,b,c,ksi,eta,zeta,cout);
  }

  // Step 5
  if(((abs(ksi)-b)>TOL) ||
     ((abs(ksi-b)<TOL) && (2.0*eta-zeta)<-TOL) || //DX 2/9/18 - more robust; precision
     ((abs(ksi+b)<TOL) && (zeta<-TOL)) ) { //DX 2/9/18 - more robust; precision
  // DX 2/9/18 [OBSOLETE] if(((abs(ksi)-b)>=TOL) ||
     // DX 2/9/18 [OBSOLETE] ((abs(ksi-b)<TOL) && (2.0*eta<zeta)) ||
     // DX 2/9/18 [OBSOLETE] ((abs(ksi+b)<TOL) && (zeta<0)) ) {
    m5(2,3)=-(SignNoZero(ksi));
    P=P*m5;
    temp1=b+c-ksi*SignNoZero(ksi);
    temp2=eta-zeta*SignNoZero(ksi);
    temp3=ksi-2.0*b*SignNoZero(ksi);
    c=temp1;
    eta=temp2;
    ksi=temp3;
    // tpx
    // _sdebug_GetNiggliCell(5,cnt,a,b,c,ksi,eta,zeta,cout);
    goto LoopHead;
  }

  // Step 6
  if(((abs(eta)-a)>TOL) || 
     ((abs(eta-a)<TOL) && (2.0*ksi-zeta)<-TOL) || //DX 2/9/18 - more robust; precision
     ((abs(eta+a)<TOL) && (zeta<-TOL))) { //DX 2/9/18 - more robust; precision
     //DX 2/9/18 [OBSOLETE] ((abs(eta-a)<TOL) && (2.0*ksi<zeta)) ||
     //DX 2/9/18 [OBSOLETE] ((abs(eta+a)<TOL) && (zeta<0))) {
    m6(1,3)=-SignNoZero(eta);
    P=P*m6;
    temp1=a+c-eta*SignNoZero(eta);
    temp2=ksi-zeta*SignNoZero(eta);
    temp3=eta-2.0*a*SignNoZero(eta);
    c=temp1;
    ksi=temp2;
    eta=temp3;
    // tpx
    // _sdebug_GetNiggliCell(6,cnt,a,b,c,ksi,eta,zeta,cout);
    goto LoopHead; //DX 2/12/18 - this was missing
  }

  // Step 7
  if(((abs(zeta)-a)>TOL) || 
     ((abs(zeta-a)<TOL) && (2.0*ksi-eta)<-TOL) || //DX 2/9/18 - more robust; precision
     ((abs(zeta+a)<TOL) && (eta<-TOL))) { //DX 2/9/18 - more robust; precision
     //DX 2/9/18 [OBSOLETE] ((abs(zeta-a)<TOL) && (2.0*ksi<eta)) ||
     //DX 2/9/18 [OBSOLETE] ((abs(zeta+a)<TOL) && (eta<0))) {
    m7(1,2)=-SignNoZero(zeta);
    P=P*m7;
    temp1=a+b-zeta*SignNoZero(zeta);
    temp2=ksi-eta*SignNoZero(zeta);
    temp3=zeta-2.0*a*SignNoZero(zeta);
    b=temp1;
    ksi=temp2;
    zeta=temp3;
    // tpx
    // _sdebug_GetNiggliCell(7,cnt,a,b,c,ksi,eta,zeta,cout);
    goto LoopHead;
  }

  // Step 8
  //DX 2/9/18 [OBSOLETE] if((ksi+eta+zeta+a+b<0) ||
  //DX 2/9/18 [OBOSLETE]   ((ksi+eta+zeta+a+b<0) && (2.0*(a+eta)+zeta>0))) {
  if((ksi+eta+zeta+a+b<-TOL) || //DX 2/9/18 - more robust; precision
     (abs(ksi+eta+zeta+a+b)<TOL && (2.0*(a+eta)+zeta>TOL))) { //DX 2/9/18 - more robust; precision
    P=P*m8;
    temp1=a+b+c+ksi+eta+zeta;
    temp2=2.0*b+ksi+zeta;
    temp3=2.0*a+eta+zeta;
    c=temp1;
    ksi=temp2;
    eta=temp3;
    // tpx
    // _sdebug_GetNiggliCell(8,cnt,a,b,c,ksi,eta,zeta,cout);
    goto LoopHead;
  }
  
  // Renormalize back to regular cell (divide by RENORM)
  xvector<double> outdat(6);
  //      outdat(1)=sqrt(a/RENORM);
  //      outdat(2)=sqrt(b/RENORM);
  //      outdat(3)=sqrt(c/RENORM);
  //      outdat(4)=acos(ksi/RENORM/2.0/outdat(2)/outdat(3));
  //      outdat(5)=acos(eta/RENORM/2.0/outdat(1)/outdat(3));
  //      outdat(6)=acos(zeta/RENORM/2.0/outdat(1)/outdat(2));
  outdat(1)=sqrt(a);
  outdat(2)=sqrt(b);
  outdat(3)=sqrt(c);
  outdat(4)=acos(ksi/2.0/outdat(2)/outdat(3));
  outdat(5)=acos(eta/2.0/outdat(1)/outdat(3));
  outdat(6)=acos(zeta/2.0/outdat(1)/outdat(2));

  // Get Niggli cell
  niggli_lat=trasp(P)*in_lat;

  // Get Q
  //  Q = P;
  // xmatrix<double> tmat(3,3);
  //  GaussJordan(Q,tmat); // Returns Q=P^-1.
  Q=inverse(P);

  // Checks
  // Make sure that a,b,c,alpha,beta,gamma are the same from
  // direct calculation and from using P to get niggli_lat.
  xvector<double> poutdat(6);
  poutdat=Getabc_angles(niggli_lat,RADIANTS);
  int flag=0;
  for(int i=1;i<=6;i++) {
    if(abs(poutdat(i)-outdat(i))>2*TOL) {flag=1;}
  }
  if(flag) {
    stringstream oss;
    oss << "ERROR: CellReduceFuncs/GetNiggliCell" << endl;
    oss << "ERROR: Lattice parameters/angles as calculated" << endl;
    oss << "ERROR: with Niggli algorithm and P do not match." << endl;
    oss << "ERROR: a,b,c,alpha,beta,gamma from direct algorithm: ";
    oss << outdat << endl;
    oss << "ERROR: a,b,c,alpha,beta,gamma from P matrix: ";
    oss << poutdat << endl;
    oss << "ERROR: Exiting." << endl;
    // output
    cerr << oss.str();
    //  cout << oss.str();
    // done
    return FALSE;
  }
  return TRUE; // perfect
}
//DX 2/13/18 - Niggli algorithm, fixed tolerances and added goto loop in step 6 - END

//DX 2/13/18 - Original Niggli routine below
bool GetNiggliCell_20180101(const xmatrix<double>& in_lat,xmatrix<double>& niggli_lat,xmatrix<double>& P,xmatrix<double>& Q) {
  // return FALSE if failed
  // rerurn TRUE if ok

  //  double RENORM=1E+6;  // DANE not used
  int MAXITER=100000;
  double TOL=1e-14;

  // Initialize matrices for tranformations (3x3).
  xmatrix<double> m1(3,3);
  m1(1,2)=1.0; m1(2,1)=1.0; m1(3,3)=-1.0;
  xmatrix<double> m2(3,3);
  m2(1,1)=-1.0; m2(2,3)=1.0; m2(3,2)=1.0;
  xmatrix<double> m3(3,3);
  // Els 1,1 2,2 3,3 are changed later for m3
  xmatrix<double> m4(3,3);
  // Els 1,1 2,2 3,3 are changed later for m4
  xmatrix<double> m5(3,3);
  m5(1,1)=1.0; m5(2,2)=1.0; m5(3,3)=1.0;
  // El 3,2 is changed later for m5
  xmatrix<double> m6(3,3);
  m6(1,1)=1.0; m6(2,2)=1.0; m6(3,3)=1.0;
  // El 3,1 is changed later for m6
  xmatrix<double> m7(3,3);
  m7(1,1)=1.0; m7(2,2)=1.0; m7(3,3)=1.0;
  // El 2,1 is changed later for m7
  xmatrix<double> m8(3,3);
  m8(1,1)=1.0; m8(2,2)=1.0; m8(1,3)=1.0; m8(2,3)=1.0; m8(3,3)=1.0;

  // Initialize a, b, c, ksi, eta, zeta
  xvector<double> indat(6);
  indat=Getabc_angles(in_lat,RADIANTS);
  double a=indat(1)*indat(1);
  double b=indat(2)*indat(2);
  double c=indat(3)*indat(3);
  double ksi=2.0*indat(2)*indat(3)*cos(indat(4));
  double eta=2.0*indat(1)*indat(3)*cos(indat(5));
  double zeta=2.0*indat(1)*indat(2)*cos(indat(6));
  // Round to RENORM decimal place to eliminate numerical errors
  //  double a=Nint(indat(1)*indat(1)*RENORM);
  //  double b=Nint(indat(2)*indat(2)*RENORM);
  //  double c=Nint(indat(3)*indat(3)*RENORM);
  //  double ksi=Nint(2.0*indat(2)*indat(3)*cos(indat(4))*RENORM);
  //  double eta=Nint(2.0*indat(1)*indat(3)*cos(indat(5))*RENORM);
  //  double zeta=Nint(2.0*indat(1)*indat(2)*cos(indat(6))*RENORM);

  // Dummy variables
  double temp,temp1,temp2,temp3;

  // Initialize tranformation matrix
  //  xmatrix<double> P(3,3);
  P(1,1)=1; P(2,2)=1; P(3,3)=1;

  int cnt=0;
  // tpx
  // _sdebug_GetNiggliCell(0,0,a,b,c,ksi,eta,zeta,cout);

  // Start loop

 LoopHead:

  cnt++;
  //  cerr << "DEBUG: Niggli() cnt=" << cnt << endl;
  if(cnt>MAXITER) {
    //     stringstream oss;
    //     // -------
    //     oss << "EEEEE  aflow " << VERSION << endl;
    //     oss << "EEEEE  ERROR: CellReduceFuncs/GetNiggliCell" << endl;
    //     oss << "EEEEE  ERROR: Number of interations greater the MAXITER = " << MAXITER << endl;
    //     oss << "EEEEE  ERROR: This seems like too many - there is probably some problem." << endl;
    //     oss << "EEEEE  ERROR:." << endl;
    //     //  oss << endl;
    //     // -------
    //     cerr << oss.str();
    //     //  cout << oss.str();
    //     // -------
    return FALSE;
  }
  
  // Step 1
  if(((a-b)>TOL) || ((abs(a-b)<TOL) && (abs(ksi)>abs(eta)) ) ) {
    temp=a;a=b;b=temp;
    temp=-ksi;ksi=-eta;eta=temp;
    P=P*m1;
    // tpx
    // _sdebug_GetNiggliCell(1,cnt,a,b,c,ksi,eta,zeta,cout);
  }

  // Step 2
  if(((b-c)>TOL) || ((abs(b-c)<TOL) && (abs(eta)>abs(zeta)) ) ) {
    temp=b;b=c;c=temp;
    temp=-eta;eta=-zeta;zeta=temp;
    P=P*m2;
    // tpx
    // _sdebug_GetNiggliCell(2,cnt,a,b,c,ksi,eta,zeta,cout);
    goto LoopHead;
  }

  // Step 3
  if((ksi*eta*zeta)>0.0) {
    m3(1,1)=SignNoZero(ksi);
    m3(2,2)=SignNoZero(eta);
    m3(3,3)=SignNoZero(zeta);
    P=P*m3;
    ksi=abs(ksi);
    eta=abs(eta);
    zeta=abs(zeta);
    // tpx
    // _sdebug_GetNiggliCell(3,cnt,a,b,c,ksi,eta,zeta,cout);
  }

  // Step 4
  // Eric's original if((abs(ksi)<TOL)||(abs(eta)<TOL)||(abs(zeta)<TOL)||(ksi*eta*zeta<=0)) {
  if(ksi*eta*zeta<=0) {
    m4(1,1)=-(SignNoZero(ksi));
    m4(2,2)=-(SignNoZero(eta));
    m4(3,3)=-(SignNoZero(zeta));
    if(abs(ksi)<TOL) m4(1,1)=m4(2,2)*m4(3,3);
    if(abs(eta)<TOL) m4(2,2)=m4(1,1)*m4(3,3);
    if(abs(zeta)<TOL) m4(3,3)=m4(1,1)*m4(2,2);
    P=P*m4;
    ksi=-abs(ksi);
    eta=-abs(eta);
    zeta=-abs(zeta);
    // tpx
    // _sdebug_GetNiggliCell(4,cnt,a,b,c,ksi,eta,zeta,cout);
  }

  // Step 5
  if(((abs(ksi)-b)>=TOL) ||
     ((abs(ksi-b)<TOL) && (2.0*eta<zeta)) ||
     ((abs(ksi+b)<TOL) && (zeta<0)) ) {
    m5(2,3)=-(SignNoZero(ksi));
    P=P*m5;
    temp1=b+c-ksi*SignNoZero(ksi);
    temp2=eta-zeta*SignNoZero(ksi);
    temp3=ksi-2.0*b*SignNoZero(ksi);
    c=temp1;
    eta=temp2;
    ksi=temp3;
    // tpx
    // _sdebug_GetNiggliCell(5,cnt,a,b,c,ksi,eta,zeta,cout);
    goto LoopHead;
  }

  // Step 6
  if(((abs(eta)-a)>TOL) || ((abs(eta-a)<TOL) && (2.0*ksi<zeta)) ||
     ((abs(eta+a)<TOL) && (zeta<0))) {
    m6(1,3)=-SignNoZero(eta);
    P=P*m6;
    temp1=a+c-eta*SignNoZero(eta);
    temp2=ksi-zeta*SignNoZero(eta);
    temp3=eta-2.0*a*SignNoZero(eta);
    c=temp1;
    ksi=temp2;
    eta=temp3;
    // tpx
    // _sdebug_GetNiggliCell(6,cnt,a,b,c,ksi,eta,zeta,cout);
  }

  // Step 7
  if(((abs(zeta)-a)>TOL) || ((abs(zeta-a)<TOL) && (2.0*ksi<eta)) ||
     ((abs(zeta+a)<TOL) && (eta<0))) {
    m7(1,2)=-SignNoZero(zeta);
    P=P*m7;
    temp1=a+b-zeta*SignNoZero(zeta);
    temp2=ksi-eta*SignNoZero(zeta);
    temp3=zeta-2.0*a*SignNoZero(zeta);
    b=temp1;
    ksi=temp2;
    zeta=temp3;
    // tpx
    // _sdebug_GetNiggliCell(7,cnt,a,b,c,ksi,eta,zeta,cout);
    goto LoopHead;
  }

  // Step 8
  if((ksi+eta+zeta+a+b<0) ||
     ((ksi+eta+zeta+a+b<0) && (2.0*(a+eta)+zeta>0))) {
    P=P*m8;
    temp1=a+b+c+ksi+eta+zeta;
    temp2=2.0*b+ksi+zeta;
    temp3=2.0*a+eta+zeta;
    c=temp1;
    ksi=temp2;
    eta=temp3;
    // tpx
    // _sdebug_GetNiggliCell(8,cnt,a,b,c,ksi,eta,zeta,cout);
    goto LoopHead;
  }

  // Renormalize back to regular cell (divide by RENORM)
  xvector<double> outdat(6);
  //      outdat(1)=sqrt(a/RENORM);
  //      outdat(2)=sqrt(b/RENORM);
  //      outdat(3)=sqrt(c/RENORM);
  //      outdat(4)=acos(ksi/RENORM/2.0/outdat(2)/outdat(3));
  //      outdat(5)=acos(eta/RENORM/2.0/outdat(1)/outdat(3));
  //      outdat(6)=acos(zeta/RENORM/2.0/outdat(1)/outdat(2));
  outdat(1)=sqrt(a);
  outdat(2)=sqrt(b);
  outdat(3)=sqrt(c);
  outdat(4)=acos(ksi/2.0/outdat(2)/outdat(3));
  outdat(5)=acos(eta/2.0/outdat(1)/outdat(3));
  outdat(6)=acos(zeta/2.0/outdat(1)/outdat(2));

  // Get Niggli cell
  niggli_lat=trasp(P)*in_lat;

  // Get Q
  //  Q = P;
  // xmatrix<double> tmat(3,3);
  //  GaussJordan(Q,tmat); // Returns Q=P^-1.
  Q=inverse(P);

  // Checks
  // Make sure that a,b,c,alpha,beta,gamma are the same from
  // direct calculation and from using P to get niggli_lat.
  xvector<double> poutdat(6);
  poutdat=Getabc_angles(niggli_lat,RADIANTS);
  int flag=0;
  for(int i=1;i<=6;i++) {
    if(abs(poutdat(i)-outdat(i))>2*TOL) {flag=1;}
  }
  if(flag) {
    stringstream oss;
    oss << "ERROR: CellReduceFuncs/GetNiggliCell" << endl;
    oss << "ERROR: Lattice parameters/angles as calculated" << endl;
    oss << "ERROR: with Niggli algorithm and P do not match." << endl;
    oss << "ERROR: a,b,c,alpha,beta,gamma from direct algorithm: ";
    oss << outdat << endl;
    oss << "ERROR: a,b,c,alpha,beta,gamma from P matrix: ";
    oss << poutdat << endl;
    oss << "ERROR: Exiting." << endl;
    // output
    cerr << oss.str();
    //  cout << oss.str();
    // done
    return FALSE;
  }
  return TRUE; // perfect
}

xstructure GetNiggliStr(const xstructure& in_str) {
  xmatrix<double> niggli_lat(3,3);
  xmatrix<double> P(3,3);
  xmatrix<double> Q(3,3);
  xstructure sstr=in_str;
  double scale=sstr.scale;
  sstr=ReScale(sstr,1.0);
  bool is_Niggli_ok=TRUE;
  is_Niggli_ok=GetNiggliCell(sstr.lattice,niggli_lat,P,Q);
  if(is_Niggli_ok==FALSE) {
    sstr=in_str;
    sstr.Niggli_has_failed=TRUE;
    return sstr; // got a problem
  } else {
    // Create new str same as before but with Niggli cell params.
    // Transform cell parameters by xyz_niggli = Q * xyz_original
    // (see ITC, ch.5).
    xstructure niggli_str=sstr;
    niggli_str.Niggli_has_failed=FALSE;
    niggli_str.lattice=niggli_lat;
    niggli_str.FixLattices();
    niggli_str.atoms.clear();

    _atom atom;
    for(int ia=0;ia<(int)sstr.atoms.size();ia++) {
      atom=sstr.atoms.at(ia);
      atom.fpos=Q*sstr.atoms.at(ia).fpos;
      atom.cpos=F2C(niggli_str.lattice,atom.fpos);
      niggli_str.atoms.push_back(atom);
    }
    // Set lattice params, atom positions, put in the 0th cell, and reset scale.
    niggli_str=BringInCell(niggli_str);
    for(int ia=0;ia<(int)sstr.atoms.size();ia++)  // NEED TO CLEAR lattice positions of atoms
      clear(sstr.atoms.at(ia).ijk);           // because the niggly reproduces the same structure
    // DONE and fix it back
    niggli_str=ReScale(niggli_str,scale);
    // DONE
    niggli_str.Niggli_calculated=TRUE;
    return niggli_str;
  }
}

xmatrix<double> GetNiggliStr(const xmatrix<double>& lattice) {
  xmatrix<double> niggli_lat(3,3);
  xmatrix<double> P(3,3);
  xmatrix<double> Q(3,3);
  bool is_Niggli_ok=TRUE;
  is_Niggli_ok=GetNiggliCell(lattice,niggli_lat,P,Q);
  if(is_Niggli_ok==FALSE) return lattice;
  return niggli_lat;
}

void _sdebug_GetNiggliCell(int step,int iter,double a,double b,double c,
			   double ksi,double eta,double zeta,ostream& sout) {
  if(0) {
    // tpx
    sout << "Step: " << step << endl;
    cout << "Iter: " << iter << endl;
    cout << "a,b,c,ksi,eta,zeta: " <<a<<" "<<b<<" "<<c <<" "<<ksi<<" "<<eta<<" "<<zeta<<endl;
    cout << "sum a+b+c: " << a+b+c << endl;
  }
}


// ***************************************************************************
// Function NiggliUnitCellForm
// ***************************************************************************
// Converts the unit cell to the standardized Niggli form.
// It is based on the GetNiggliStr() function written by Eric Wu and
// Dane Morgan, adapted by Stefano Curtarolo and present in aflow_pflow_funcs.cpp
xstructure NiggliUnitCellForm(const xstructure& a) {
  xstructure b(a);
  b.NiggliUnitCellForm();
  return b;
}

void xstructure::NiggliUnitCellForm(void) {
  if(Niggli_avoid==TRUE) return;
  if(Niggli_calculated==FALSE) {
    *this=GetNiggliStr(*this);
    Niggli_calculated=TRUE;
    FixLattices();
    Niggli_has_failed=FALSE;
  } else {
    if(!XHOST.QUIET) cout << "00000  MESSAGE NIGGLI Form already Calculated: skipping" << endl;
  }
}

xmatrix<double> NiggliUnitCellForm(const xmatrix<double>& lattice) {
  return GetNiggliStr(lattice);
}

// **************************************************************************
// Function MinkowskiBasisReduction
// **************************************************************************
// This routine takes a set of basis vectors (that form a lattice)
// and reduces them so that they form the shortest possible basis.
// The reduction is performed so that each vector "a_i" is a close
// as possible to the origin while remaining in the affine plane which
// is defined by "a_j", "a_k" but shifted by "a_i", for any choice
// of even permutations of i,j,k in 1,2,3.
// See Lecture notes in computer science, ISSN 0302-974, ANTS - VI :
// algorithmic number theory, 2004, vol. 3076, pp. 338-357
// ISBN 3-540-22156-5
// Written by Gus Hart in F90, recoded by SC in C++ (Sep/08).
xstructure MinkowskiBasisReduction(const xstructure& a) {
  xstructure b(a);
  b.MinkowskiBasisReduction();
  return b;
}

void xstructure::MinkowskiBasisReduction(void) {
  if(LatticeReduction_avoid==TRUE) return;
  if(Minkowski_avoid==TRUE) return;
  if(Minkowski_calculated==FALSE) {
    xmatrix<double> basis(3,3);
    double old_scale=scale;
    ReScale(1.0);
    basis=trasp(lattice);
    basis=gusstd::reduce_to_shortest_basis(basis);
    lattice=trasp(basis);
    FixLattices(); // get f2c c2f and klattice
    for(uint i=0;i<atoms.size();i++)
      atoms[i].fpos=C2F(lattice,atoms[i].cpos);
    // atoms[i]=C2F(lattice,atoms[i]);  // it works, same just curiosity
    // atoms[i]=C2F(atoms[i]);          // it works, same just curiosity
    ReScale(old_scale);
    FixLattices();
    Minkowski_calculated=TRUE;
    Minkowski_has_failed=FALSE;
  } else {
    //    if(!QUIET) cout << "00000  MESSAGE MINKOWSKI Basis Reduction already Calculated: skipping" << endl;
  }
}

xmatrix<double> MinkowskiBasisReduction(const xmatrix<double>& lattice) {
  xmatrix<double> basis(3,3);
  basis=lattice;
  basis=trasp(basis);
  basis=gusstd::reduce_to_shortest_basis(basis);
  basis=trasp(basis);
  return basis;
}

// ***************************************************************************
// Function LatticeReduction
// ***************************************************************************
// Lattice Reduction to Max Orthogonality (MINK) and then Niggly Form
xstructure LatticeReduction(const xstructure& a) {
  xstructure b(a);
  b.LatticeReduction();
  return b;
}

void xstructure::LatticeReduction(void) {
  if(LatticeReduction_avoid==TRUE) return;
  if(LatticeReduction_calculated==FALSE) {
    MinkowskiBasisReduction();            // Minkowski first
    NiggliUnitCellForm();                 // Niggli Second
    FixLattices();                        // fix f2c c2f and K-space
    BringInCell();                        // bring atoms in new basis
    LatticeReduction_calculated=TRUE;
    LatticeReduction_has_failed=FALSE;
  } else {
    //   if(!QUIET) cout << "00000  MESSAGE LATTICE Basis Reduction already Calculated: skipping" << endl;
  }
}

xmatrix<double> LatticeReduction(const xmatrix<double>& lattice) {
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  xmatrix<double> basis(3,3);
  basis=lattice;
  basis=MinkowskiBasisReduction(basis); // Minkowski first
  basis=NiggliUnitCellForm(basis);      // Niggli Second
  if(LDEBUG) cerr << "WARNING: remeber to FixLattices for f2c,c2f and reciprocal and bring atoms in cell " << endl;
  return basis;
}

// **************************************************************************
// Function BringInCell
// **************************************************************************
// these routines take atoms or structures and brings them
// to the unit cell (incell). There is some overloading for
// structures. SC Aug2007
// EDITED BY COREY AND DAVID to include a tolerance that converts from Cartesian
// space to direct space via covariant and contravariant transformations.
#define _EPS_roundoff_ 0.001
#define _incellcutoff_ (1.0-_EPS_roundoff_)
//DX and CO, IF EPSILON IS PROVIDED, THEN WE ASSUME YOU WANT US TO MOVE THE ATOMS, OTHERWISE IT'S A HARD CUT OFF

xvector<double> BringInCell(const xvector<double>& v_in,double epsilon) {
  //DX and CO - START
  return BringInCell_20160101(v_in,epsilon);
}

xvector<double> BringInCell_20161115(const xvector<double>& v_in,double epsilon) {
  return BringInCell_20160101(v_in,epsilon); //SYM::mod_one_xvec(v_in);
}

xvector<double> BringInCell_20160101(const xvector<double>& v_in,double epsilon) {
  double incelleps=1.0-epsilon;
  xvector<double> v_out(v_in.urows,v_in.lrows);
  for(int i=v_out.lrows;i<=v_out.urows;i++) {
    v_out(i)=v_in(i);
    while(v_out(i)> incelleps) v_out(i)-=1.0;
    while(v_out(i)< 0.0)       v_out(i)+=1.0;
    if(abs(v_out(i))<epsilon)  v_out(i)=0.0;
    if(v_out(i)> incelleps)    v_out(i)=0.0;
  }
  //v_out=roundoff(v_out,_EPS_sym_);
  //DX and CO - END
  return v_out;
}

xvector<double> BringInCell(const xvector<double>& v_in) {
  //DX and CO - START
  return SYM::mod_one_xvec(v_in); //hard cutoff
  //_EPS_sym_ is the universal tolerance here
  //return BringInCell(v_in,_EPS_sym_);
  //DX and CO - END
}

//DX anc CO - START
_atom BringInCell(const _atom& atom_in,const xmatrix<double>& lattice,double epsilon) {
 
  return BringInCell_20160101(atom_in,lattice,epsilon);
}

_atom BringInCell_20161115(const _atom& atom_in,const xmatrix<double>& lattice,double epsilon) {
  return BringInCell_20160101(atom_in,lattice,epsilon);//BringInCell(atom_in,lattice);
}
//DX and CO - END

_atom BringInCell_20160101(const _atom& atom_in,const xmatrix<double>& lattice,double epsilon) {
  double _incelleps=1.0;
  if(epsilon>0.0) _incelleps=1.0-epsilon;else _incelleps=1.0;
  _atom atom;
  atom=atom_in;
  atom.ijk=atom_in.ijk;
  for(int i=1;i<=3;i++) {
    while(atom.fpos(i)> _incelleps) {
      atom.fpos(i)-=1.0;
      atom.ijk(i)++;
    }
    while(atom.fpos(i)<0.0) {
      atom.fpos(i)+=1.0;
      atom.ijk(i)--;
    }
    if(epsilon>0.0) { // roundoff only if epsilon>0.0
      if(abs(atom.fpos(i))<epsilon) atom.fpos(i)=0.0;
      if(atom.fpos(i)>_incelleps) atom.fpos(i)=0.0;
    }
  }
  atom.cpos=F2C(lattice,atom.fpos);
  atom.isincell=TRUE;
  return atom;
}

//DX and CO - START
_atom BringInCell(const _atom& atom_in,const xmatrix<double>& lattice) {
  return BringInCell_20161115(atom_in,lattice);
}

_atom BringInCell_20161115(const _atom& atom_in,const xmatrix<double>& lattice) {
  //_atom atom;
  //atom=atom_in;
  //atom.ijk=atom_in.ijk; //just to be sure //it's part of the assignment operator now
  _atom atom=SYM::mod_one_atom(atom_in);
  atom.cpos=F2C(lattice,atom.fpos);
  atom.isincell=TRUE;
  return atom;
}
//DX and CO - END

_atom BringInCell_20160101(const _atom& atom_in,const xmatrix<double>& lattice) {
  return BringInCell(atom_in,lattice,_EPS_roundoff_);
}

//_atom _BringInCell(const _atom& atom_in,const xmatrix<double>& lattice) {
//  return BringInCell(atom_in,lattice);
//}

//DX and CO - START
xstructure BringInCell(const xstructure& a,double epsilon) {
  return BringInCell_20160101(a,epsilon);
}

xstructure BringInCell_20161115(const xstructure& a,double epsilon) {
  xstructure b(a); // copies everything
  for(int i=0;i<(int)a.atoms.size();i++){
    b.atoms[i]=BringInCell(a.atoms[i],a.lattice,epsilon);//BringInCell(a.atoms[i],a.lattice);
  }
  return b;
}

xstructure BringInCell_20160101(const xstructure& a,double epsilon) {
  xstructure b(a); // copies everything
  for(int i=0;i<(int)a.atoms.size();i++) {
    b.atoms[i]=BringInCell(a.atoms[i],a.lattice,epsilon);//BringInCell(a.atoms[i],a.lattice);
  }
  return b;
}



xstructure BringInCell(const xstructure& a) {
  //DX and CO - START
  return BringInCell_20161115(a);
}

//SUPER SLOW, AVOID!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//USE THE ONE THAT ONLY COPIES ATOMS, COREY
//a.sym_eps is STRONGLY PREFERRED
//however, it requires knowledge (and calculating of the full symmetry, i.e. slow factor group)
//therefore, if you call the default without having a.sym_eps calculated, the assumption
//is that a constant tolerance is good enough (speed)
//if you want to do this correctly, calculate the full symmetry of the structure first (CalculateFullSymmetry() in
//aflow_aconvasp_main.cpp)
xstructure BringInCell_20161115(const xstructure& a) {
  xstructure b=a;
  //  double _eps_;
  //  if(a.sym_eps!=AUROSTD_NAN){ //Tolerance came from user or was calculated
  //    _eps_=a.sym_eps;
  //  }
  //  else{
  //    _eps_=_EPS_sym_;
  //    // Calculate point group/space group to find correct tolerance for system
  //    //SYM::CalculatePointGroup(FileMESSAGE,a,aflags,_write_,osswrite,oss);       
  //    //_eps_=a.sym_eps;
  //  }
  //  return BringInCell(a,_eps_);
  for(uint i=0; i<b.atoms.size(); i++){
    b.atoms[i].fpos = SYM::mod_one_xvec(b.atoms[i].fpos);
  } 
  return b;
}

xstructure BringInCell_20160101(const xstructure& a) {
  return BringInCell(a,_EPS_roundoff_);
}
//DX and CO - END

void xstructure::BringInCell(double epsilon) {
  _atom BringInCell(const _atom& atom_in,const xmatrix<double>& lattice,double epsilon);
  for(int i=0;i<(int)atoms.size();i++){
    atoms[i]=BringInCell(atoms[i],lattice,epsilon);
  }
}
void xstructure::BringInCell(void) {
  //DX and CO - START
  _atom BringInCell(const _atom& atom_in,const xmatrix<double>& lattice);
  for(int i=0;i<(int)atoms.size();i++){
    atoms[i]=BringInCell(atoms[i],lattice);
  }
  //TOO SLOW
  //xstructure BringInCell(const xstructure& a);
  //*this=BringInCell(*this);
  //BringInCell(_EPS_sym_);  //CO
  //DX and CO - END
}

// ***************************************************************************
// Function BringInCompact
// ***************************************************************************
// Make a structure where all the atoms are all the
// atoms are mapped through the unit and neighbours cells
// to minimixe the shortest possible bond with an adjacent atom
// This option is very useful if you run big and complicate
// molecules where atoms exit of the unit cell and you have
// problems understanding where they are because visualization
// packages do not show bonds anymore ...
// Anyway, it is easier to test than to describe. (SC 6 Aug 04)
xstructure BringInCompact(const xstructure& a) {
  xstructure sstr=a;
  sstr.BringInCompact();
  return sstr;
}

xstructure _old_BringInCompact(const xstructure& a) {
  xstructure sstr=a;
  // double TOL=1e-15;
  double min_bond,bond;
  // For direct coordinates
  for(uint iat=0;iat<sstr.atoms.size();iat++)
    sstr.atoms.at(iat)=BringInCell(sstr.atoms.at(iat),sstr.lattice);
  xvector<double> adref1pos(3),adtstpos(3),adtrgpos(3);
  xvector<double> acref1pos(3),actstpos(3),actrgpos(3);
  for(uint i=1;i<sstr.atoms.size();i++) {        // scan all the atoms to move except the first
    adtrgpos=sstr.atoms.at(i).fpos;
    actrgpos=sstr.atoms.at(i).cpos;
    min_bond=1.0e6;
    for(uint ii=0;ii<i;ii++) {                      // scan over all the reference atoms
      adref1pos=sstr.atoms.at(ii).fpos;
      for(int ic=1;ic<=3;ic++) {                     // scan over all the reference atoms
	acref1pos(ic)=0.0;                          // scan over all the reference atoms
	for(int jc=1;jc<=3;jc++)
	  acref1pos(ic)=acref1pos(ic)+adref1pos(jc)*sstr.lattice(jc,ic);
	// scan over all the reference atoms
      }                                             // scan over all the reference atoms

      for(int i1=-1;i1<=1;i1++)                     // roll over first neighbour cells
	for(int j1=-1;j1<=1;j1++)                   // roll over first neighbour cells
	  for(int k1=-1;k1<=1;k1++) {               // roll over first neighbour cells
	    adtstpos(1)=sstr.atoms.at(i).fpos(1)+i1;
	    adtstpos(2)=sstr.atoms.at(i).fpos(2)+j1;
	    adtstpos(3)=sstr.atoms.at(i).fpos(3)+k1;      // test the atom
	    for(int ic=1;ic<=3;ic++) {                // test the atom
	      actstpos(ic)=0.0;                      // test the atom
	      for(int jc=1;jc<=3;jc++)               // test the atom
		actstpos(ic)=actstpos(ic)+adtstpos(jc)*sstr.lattice(jc,ic);   // test the atom
	    }                                        // test the atom  
	    bond=aurostd::modulus(actstpos-acref1pos);    // test the bond distance
	    // if((bond<min_bond) {                  // if it is OK then DO IT !
	    if((bond<1.03*min_bond && abs((int) i-(int) ii)<10) || (bond<0.98*min_bond) ) {
	      // if it is OK then DO IT !
	      min_bond=bond;                                      // update
	      adtrgpos=adtstpos;
	      actrgpos=actstpos;
	    }
	  }
    }
    // now the TRG atom has all the information about the bst position of the atom to be moved.
    // than EAT IT !
    for(int j=1;j<=3;j++) {
      sstr.atoms.at(i).fpos(j)=adtrgpos(j);
      sstr.atoms.at(i).cpos(j)=actrgpos(j);
    }
  }
  // For Cartesian coordinates (get from direct coords);
  for(uint iat=0;iat<sstr.atoms.size();iat++)
    sstr.atoms.at(iat).cpos=F2C(sstr.lattice,sstr.atoms.at(iat).fpos);
  return sstr;
}

void xstructure::BringInCompact(void) {
  double min_bond,bond;
  // For direct coordinates
  BringInCell();
  xvector<double> adref1pos(3),adtstpos(3),adtrgpos(3);
  xvector<double> acref1pos(3),actstpos(3),actrgpos(3);
  for(uint i=1;i<atoms.size();i++) {        // scan all the atoms to move except the first
    adtrgpos=atoms.at(i).fpos;
    actrgpos=atoms.at(i).cpos;
    min_bond=1.0e6;
    for(uint ii=0;ii<i;ii++) {                      // scan over all the reference atoms
      adref1pos=atoms.at(ii).fpos;
      for(int ic=1;ic<=3;ic++) {                     // scan over all the reference atoms
	acref1pos(ic)=0.0;                          // scan over all the reference atoms
	for(int jc=1;jc<=3;jc++)
	  acref1pos(ic)=acref1pos(ic)+adref1pos(jc)*lattice(jc,ic);
	// scan over all the reference atoms
      }                                             // scan over all the reference atoms

      for(int i1=-1;i1<=1;i1++)                     // roll over first neighbour cells
	for(int j1=-1;j1<=1;j1++)                   // roll over first neighbour cells
	  for(int k1=-1;k1<=1;k1++) {               // roll over first neighbour cells
	    adtstpos(1)=atoms.at(i).fpos(1)+i1;
	    adtstpos(2)=atoms.at(i).fpos(2)+j1;
	    adtstpos(3)=atoms.at(i).fpos(3)+k1;      // test the atom
	    for(int ic=1;ic<=3;ic++) {                // test the atom
	      actstpos(ic)=0.0;                      // test the atom        
	      for(int jc=1;jc<=3;jc++)               // test the atom
		actstpos(ic)=actstpos(ic)+adtstpos(jc)*lattice(jc,ic);   // test the atom
	    }                                        // test the atom  
	    bond=aurostd::modulus(actstpos-acref1pos);    // test the bond distance
	    // if((bond<min_bond) {                  // if it is OK then DO IT !
	    if((bond<1.03*min_bond && abs((int) i-(int) ii)<10) || (bond<0.98*min_bond) ) {
	      // if it is OK then DO IT !
	      min_bond=bond;                                      // update
	      adtrgpos=adtstpos;
	      actrgpos=actstpos;
	    }
	  }
    }
    // now the TRG atom has all the information about the bst position of the atom to be moved.
    // than EAT IT !
    for(int j=1;j<=3;j++) {
      atoms.at(i).fpos(j)=adtrgpos(j);
      atoms.at(i).cpos(j)=actrgpos(j);
    }
  }
  // For Cartesian coordinates (get from direct coords);
  for(uint iat=0;iat<atoms.size();iat++)
    atoms.at(iat).cpos=F2C(lattice,atoms.at(iat).fpos);
}

// ***************************************************************************
// Function BringInWignerSeitz
// ***************************************************************************
// Make a structure where all the atoms are
// mapped to their images in the Wigner-Seitz cell.(SC 10Jan04)
void xstructure::BringInWignerSeitz(void) {
  xvector<double> rat(3),a1(3),a2(3),a3(3),a12(3),a31(3),a23(3),a123(3);
  double na1,na2,na3,na12,na31,na23,na123;
  a1=lattice(1);na1=aurostd::modulus(a1);
  a2=lattice(2);na2=aurostd::modulus(a2);
  a3=lattice(3);na3=aurostd::modulus(a3);
  a12=lattice(1)+lattice(2);na12=aurostd::modulus(a12);
  a31=lattice(1)+lattice(3);na31=aurostd::modulus(a31);
  a23=lattice(2)+lattice(3);na23=aurostd::modulus(a23);
  a123=lattice(1)+lattice(2)+lattice(3);na123=aurostd::modulus(a123);
  double proj_a1,proj_a2,proj_a3,proj_a12,proj_a31,proj_a23,proj_a123;

  BringInCell();

  for(uint iat=0;iat<atoms.size();iat++) {
    for(int i=-1;i<=1;i++)
      for(int j=-1;j<=1;j++)
	for(int k=-1;k<=1;k++) {
	  rat=atoms.at(iat).cpos+((double)i)*a1+((double)j)*a2+((double)k)*a3;
	  proj_a1=scalar_product(rat,a1)/na1/na1;
	  proj_a2=scalar_product(rat,a2)/na2/na2;
	  proj_a3=scalar_product(rat,a3)/na3/na3;
	  proj_a12=scalar_product(rat,a12)/na12/na12;
	  proj_a31=scalar_product(rat,a31)/na31/na31;
	  proj_a23=scalar_product(rat,a23)/na23/na23;
	  proj_a123=scalar_product(rat,a123)/na123/na123;
	  if((proj_a1>-0.5 && proj_a1<=0.5) &&
	     (proj_a2>-0.5 && proj_a2<=0.5) &&
	     (proj_a3>-0.5 && proj_a3<=0.5) &&
	     (proj_a12>-0.5 && proj_a12<=0.5) &&
	     (proj_a31>-0.5 && proj_a31<=0.5) &&
	     (proj_a23>-0.5 && proj_a23<=0.5) &&
	     (proj_a123>-0.5 && proj_a123<=0.5)) {
	    atoms.at(iat).cpos(1)=rat(1);
	    atoms.at(iat).cpos(2)=rat(2);
	    atoms.at(iat).cpos(3)=rat(3);
	    i=10;j=10;k=10;
	  }
	}
    atoms.at(iat).fpos=C2F(lattice,atoms.at(iat).cpos);
  }
}

xstructure BringInWignerSeitz(const xstructure& a) {
  xstructure str=a;
  str.BringInWignerSeitz();
  return str;
}


// **************************************************************************
// Function GetPrimitive
// Function IsTranslationFVector& IsTranslationCVector
// ***************************************************************************
// This function returns 1 if the given vector is a translation
// vector of the given structure and 0 otherwise. Input translation
// vector is expectd to be in fractional coordinates.
// written by Stefano Curtarolo (superfast edition)
//#define IsTranslationFVector IsTranslationFVectorFAST
//#define IsTranslationFVector IsTranslationFVectorORIGINAL
//#define IsTranslationFVector IsTranslationFVectorFAST_2011
#define IsTranslationFVector IsTranslationFVectorORIGINAL_2011

bool IsTranslationFVectorFAST(const xstructure& a, const xvector<double>& ftvec) {
  if(a.equiv_fpos_epsilon<1.0e-12) { cerr << "ERROR:  Zero tolerance in aflow_xatom.cpp: IsTranslationFVectorFAST " << a.equiv_fpos_epsilon << endl; exit(0);}
  double tolerance=a.equiv_fpos_epsilon;
  if(aurostd::modulus(ftvec)<=tolerance) return TRUE;

  uint CntGoodTrans=0;
  xvector<double> ftpos(3),diff(3);
  bool found_atom=FALSE,found_tvec=TRUE;
  xvector<int> types(0,a.atoms.size());       // mirror to enhance speed !
  for(uint iat=0;iat<a.atoms.size();iat++)    // mirror to enhance speed !
    types(iat)=a.atoms.at(iat).type;          // mirror to enhance speed !
  for(uint iat=0;iat<a.atoms.size()&&found_tvec;iat++) {
    // Get translated position and shift back into unit cell.
    //  ftpos=a.atoms.at(iat).fpos+ftvec;
    // ftpos=BringInCell(ftpos,tolerance);  // no because it will be needed later...
    // Find closest atom in unit cell to translated position (should only be one atom less than tolerance).
    found_atom=FALSE;
    for(uint jat=0;jat<a.atoms.size()&&!found_atom&&found_tvec;jat++) {
      if(types[iat]==types[jat]) {
	//        diff=ftpos-a.atoms.at(jat).fpos;
 	diff=BringInCell(a.atoms.at(iat).fpos-a.atoms.at(jat).fpos+ftvec,tolerance);
        // If the translated atom maps onto another and its type is the same as the atom mapped onto then increment.
        if(aurostd::modulus(diff)<=tolerance) {CntGoodTrans++;found_atom=TRUE;}
      }
    } // For jat
    if(found_atom==FALSE) found_tvec=FALSE; // at least one atom does not have a copy...
  } // For iat
  // If we counted a good translation for every atoms then we have a translation vector.
  // cout << "GOODs=" << CntGoodTrans << endl;
  if(CntGoodTrans==a.atoms.size()) return TRUE;
  return FALSE;
}

bool IsTranslationFVectorORIGINAL(const xstructure& a, const xvector<double>& ftvec) {
  if(a.equiv_fpos_epsilon<1.0e-12) { cerr << "ERROR:  Zero tolerance in aflow_xatom.cpp: IsTranslationFVectorORIGINAL" << endl; exit(0);}
  double tolerance=a.equiv_fpos_epsilon;
  if(aurostd::modulus(ftvec)<=tolerance) return TRUE;
  uint CntGoodTrans=0;
  xvector<double> ftpos(3),diff(3);
  xvector<int> types(0,a.atoms.size());       // mirror to enhance speed !
  for(uint iat=0;iat<a.atoms.size();iat++)    // mirror to enhance speed !
    types(iat)=a.atoms.at(iat).type;          // mirror to enhance speed !
  for(uint iat=0;iat<a.atoms.size();iat++) {
    // Get translated position and shift back into unit cell.
    // ftpos=a.atoms.at(iat).fpos+ftvec;
    // ftpos=BringInCell(ftpos,tolerance); // no because it will be needed later...
    // Find closest atom in unit cell to translated position (should only be one atom less than tolerance).
    for(uint jat=0;jat<a.atoms.size();jat++) {
      if(types[iat]==types[jat]) {
	//        diff=ftpos-a.atoms.at(jat).fpos;
  	diff=BringInCell(a.atoms.at(iat).fpos-a.atoms.at(jat).fpos+ftvec,tolerance);
	// If the translated atom maps onto another and its type is the same as the atom mapped onto then increment.
        if(aurostd::modulus(diff)<=tolerance) CntGoodTrans++;
      }
    } // For jat
  } // For iat
  // If we counted a good translation for every atoms then we have a translation vector.
  //  cout << "GOODs=" << CntGoodTrans << endl;
  if(CntGoodTrans==a.atoms.size()) return TRUE;
  return FALSE;
}

bool IsTranslationFVectorFAST_2011(const xstructure& a, const xvector<double>& ftvec) {
  if(a.equiv_fpos_epsilon<1.0e-12) { cerr << "ERROR:  Zero tolerance in aflow_xatom.cpp: IsTranslationFVectorFAST " << a.equiv_fpos_epsilon << endl; exit(0);}
  double tolerance=a.equiv_fpos_epsilon;
  if(aurostd::modulus(ftvec)<=tolerance) return TRUE;
  uint CntGoodTrans=0;
  xvector<double> ftpos(3),diff(3);
  bool found_atom=FALSE,found_tvec=TRUE;
  xvector<int> types(0,a.atoms.size());       // mirror to enhance speed !
  for(uint iat=0;iat<a.atoms.size();iat++)    // mirror to enhance speed !
    types(iat)=a.atoms.at(iat).type;          // mirror to enhance speed !
  for(uint iat=0;iat<a.atoms.size()&&found_tvec;iat++) {
    // Get translated position and shift back into unit cell.
    ftpos=a.atoms.at(iat).fpos+ftvec;
    ftpos=BringInCell(ftpos);
    // Find closest atom in unit cell to translated position (should only be one atom less than tolerance).
    found_atom=FALSE;
    for(uint jat=0;jat<a.atoms.size()&&!found_atom&&found_tvec;jat++) {
      if(types[iat]==types[jat]) {
        diff=ftpos-a.atoms.at(jat).fpos;
        // If the translated atom maps onto another and its type is the same as the atom mapped onto then increment.
        if(aurostd::modulus(diff)< tolerance ) {CntGoodTrans++;found_atom=TRUE;}
      }
    } // For jat
    if(found_atom==FALSE) found_tvec=FALSE; // at least one atom does not have a copy...
  } // For iat
  // If we counted a good translation for every atoms then we have a translation vector.
  if(CntGoodTrans==a.atoms.size()) return TRUE;
  return FALSE;
}

bool IsTranslationFVectorORIGINAL_2011(const xstructure& a, const xvector<double>& ftvec) {
  if(a.equiv_fpos_epsilon<1.0e-12) { cerr << "ERROR:  Zero tolerance in aflow_xatom.cpp: IsTranslationFVectorORIGINAL" << endl; exit(0);}
  double tolerance=a.equiv_fpos_epsilon;
  if(aurostd::modulus(ftvec)<=tolerance) return TRUE;

  uint CntGoodTrans=0;
  xvector<double> ftpos(3),diff(3);
  xvector<int> types(0,a.atoms.size());       // mirror to enhance speed !
  for(uint iat=0;iat<a.atoms.size();iat++)    // mirror to enhance speed !
    types(iat)=a.atoms.at(iat).type;          // mirror to enhance speed !
  for(uint iat=0;iat<a.atoms.size();iat++) {
    // Get translated position and shift back into unit cell.
    ftpos=a.atoms.at(iat).fpos+ftvec;
    ftpos=BringInCell(ftpos);
    // Find closest atom in unit cell to translated position (should only be one atom less than tolerance).
    for(uint jat=0;jat<a.atoms.size();jat++) {
      if(types[iat]==types[jat]) {
        diff=ftpos-a.atoms.at(jat).fpos;
        // If the translated atom maps onto another and its type is the same as the atom mapped onto then increment.
        if(aurostd::modulus(diff)< tolerance ) CntGoodTrans++;
      }
    } // For jat
  } // For iat
  // If we counted a good translation for every atoms then we have a translation vector.
  if(CntGoodTrans==a.atoms.size()) return TRUE;
  return FALSE;
}


bool IsTranslationCVector(const xstructure& a, const xvector<double>& ctvec) {
  if(a.equiv_fpos_epsilon<1.0e-12) { cerr << "ERROR:  Zero tolerance in aflow_xatom.cpp: IsTranslationCVector" << endl; exit(0);}
  /* Input translation vector is expectd to be in cartesian coordinates. */
  return IsTranslationFVector(a,C2F(a.lattice,ctvec));
}

// **************************************************************************
// Function GetPrimitive
// **************************************************************************
// This funtion returns a structure object where everything
// has been transformed into a primitive lattice.
// Algorithm
// Construct a list of all candidate primitive lattice vectors.
// This list consists of the present lattice vectors and all basis
// vectors where translation by that basis vector leaves all the
// basis atom types unchanged. Loop over all possible triads in this list.
// For each triad find the volume. Store a list of all triads with the
// minimum volume. Then select the triad for which the first triad vector
// has the maximal projection onto the first lattice vector.  If there
// are more than one of these do the same for the second lattice vector,
// and then the third if needed. Then just pick one if there are still
// multiple candidates.

xstructure GetPrimitiveSINGLE(const xstructure& a,double tolerance);
xstructure GetPrimitiveMULTITHREAD(const xstructure& a,double tolerance);
xstructure GetPrimitive(const xstructure& a);

xstructure GetPrimitive(const xstructure& _a,double tolerance) {
  // return GetPrimitiveMULTITHREAD(a);
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  if(LDEBUG) cerr << "GetPrimitive(const xstructure& _a,double tolerance): BEGIN" << endl;
  
  bool bool_backup=AFLOW_PTHREADS::FLAG;
  uint uint_backup=AFLOW_PTHREADS::MAX_PTHREADS;
  xstructure a(_a);
  if(tolerance<=0.0) a.equiv_fpos_epsilon=_EQUIV_FPOS_EPS_; else a.equiv_fpos_epsilon=tolerance;

  if(AFLOW_PTHREADS::RUNNING==0 && !AFLOW_PTHREADS::FLAG) { // force to save time
    AFLOW_PTHREADS::FLAG=TRUE;   // put multithread by default in
    AFLOW_PTHREADS::MAX_PTHREADS=XHOST.CPU_Cores;
  }
  // calc
  xstructure b;
  if(AFLOW_PTHREADS::FLAG==TRUE) {
    //    cerr << "GetPrimitiveMULTITHREAD a.equiv_fpos_epsilon=" << a.equiv_fpos_epsilon << endl;
    b=GetPrimitiveMULTITHREAD(a,tolerance);
  } else {
    //    cerr << "GetPrimitiveSINGLE a.equiv_fpos_epsilon=" << a.equiv_fpos_epsilon << endl;
    b=GetPrimitiveSINGLE(a,tolerance);
  }
  // restore and exit
  AFLOW_PTHREADS::FLAG=bool_backup;
  AFLOW_PTHREADS::MAX_PTHREADS=uint_backup;
  return b;
}

xstructure GetPrimitive(const xstructure& a) {
  return GetPrimitive(a,_EQUIV_FPOS_EPS_);
}

pthread_mutex_t mutex_XATOM=PTHREAD_MUTEX_INITIALIZER;
#define _PTHREAD_FLUSH_TIME_ 1

typedef struct {
  int      ITHREAD;
  int      THREADS_MAX;
  xstructure *pstr;
  int      ispecie_min;
  vector<xvector<double> > *ptvector;
  xmatrix<double> *polattice;
  int      itbusy;
  int      step;
} _threaded_GETTVECTORS_params;

void *_threaded_GetTvectors(void *ptr) {
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  _threaded_GETTVECTORS_params* pparams;
  xmatrix<double> plattice(3,3);
  xvector<double> fdisp(3),cdisp(3);
  pparams = (_threaded_GETTVECTORS_params*) ptr;
  (pparams->itbusy)=TRUE;
  // pthread_mutex_lock(&mutex_XATOM);
  AFLOW_PTHREADS::RUNNING++;
  // cerr << "AFLOW_PTHREADS::RUNNING=" << AFLOW_PTHREADS::RUNNING << endl;
  // pthread_mutex_unlock(&mutex_XATOM);
  // CODE BEGIN
  //  cerr << "debug " << (pparams->ITHREAD) << "/" << (pparams->THREADS_MAX) << endl;

  double sstr_volume=(*pparams->pstr).Volume();

  // CODE STEP1 BEGIN
  if(pparams->step==1) {
    if(LDEBUG) pthread_mutex_lock(&mutex_XATOM);   // LOCK
    if(LDEBUG) cerr << "*_threaded_GetTvectors STEP1: " << AFLOW_PTHREADS::RUNNING << " " << pparams->ITHREAD << endl;
    if(LDEBUG) pthread_mutex_unlock(&mutex_XATOM); // UNLOCK
    for(uint iat1=0+(pparams->ITHREAD);iat1<(*pparams->pstr).atoms.size();iat1+=(pparams->THREADS_MAX)) { // does modulus thread max
      if((*pparams->pstr).atoms.at(iat1).type==(pparams->ispecie_min))
        for(uint iat2=0;iat2<(*pparams->pstr).atoms.size();iat2++)
          if((*pparams->pstr).atoms.at(iat2).type==(pparams->ispecie_min)) {
            fdisp=(*pparams->pstr).atoms.at(iat2).fpos-(*pparams->pstr).atoms.at(iat1).fpos;
            cdisp=(*pparams->pstr).atoms.at(iat2).cpos-(*pparams->pstr).atoms.at(iat1).cpos;
            if(aurostd::modulus(fdisp)>0.01 && aurostd::modulus(cdisp)>0.01)
              if(IsTranslationFVector((*pparams->pstr),fdisp)) {
                pthread_mutex_lock(&mutex_XATOM);
                (*pparams->ptvector).push_back(cdisp);
                pthread_mutex_unlock(&mutex_XATOM);
              }
          }
    }
  } // CODE STEP1 END
  // CODE STEP2 BEGIN
  if(pparams->step==2) {
    if(LDEBUG) pthread_mutex_lock(&mutex_XATOM);   // LOCK
    if(LDEBUG) cerr << "*_threaded_GetTvectors STEP2: " << AFLOW_PTHREADS::RUNNING << " " << pparams->ITHREAD << endl;
    if(LDEBUG) pthread_mutex_unlock(&mutex_XATOM); // UNLOCK
    for(uint iu=0+(pparams->ITHREAD);iu<(*pparams->ptvector).size();iu+=(pparams->THREADS_MAX)) { // does modulus thread max
      for(uint i=1;i<=3;i++) plattice[1][i]=(*pparams->ptvector).at(iu)[i];
      for(uint iv=0;iv<(*pparams->ptvector).size()&& iv!=iu;iv++) {
        for(uint i=1;i<=3;i++) plattice[2][i]=(*pparams->ptvector).at(iv)[i];
        for(uint iw=0;iw<(*pparams->ptvector).size()&& iw!=iv && iw!=iu;iw++) {
          for(uint i=1;i<=3;i++) plattice[3][i]=(*pparams->ptvector).at(iw)[i];
          if(det(plattice)>0.999 && det(plattice)<sstr_volume) {   // no coplanar and contain at least 1 atom and smaller than the original cell
            if(aurostd::isinteger(sstr_volume/det(plattice))) {    // integer ratio of volumes
              if(det(plattice)<det((*pparams->polattice))) {                    // better than before
                if(LDEBUG) cout<<"DEBUG"<<iu<<","<<iv<<","<<iw<<" "<< sstr_volume<<" "<<det(plattice)<<" "<<sstr_volume/det(plattice)<<endl;
                if(isdifferent(plattice,(*pparams->polattice),0.0001)) {
                  plattice=MinkowskiBasisReduction(plattice);      // Minkowski first
                  plattice=NiggliUnitCellForm(plattice);           // Niggli Second
                  if(isdifferent(plattice,(*pparams->polattice),0.0001)) {
                    pthread_mutex_lock(&mutex_XATOM);   // LOCK
                    (*pparams->polattice)=plattice;
                    pthread_mutex_unlock(&mutex_XATOM); // UNLOCK
                  }
                }
              }
            }
          }
        }
      }
    }
  } // CODE STEP2 END

  (pparams->itbusy)=FALSE;
  AFLOW_PTHREADS::RUNNING--;
  aurostd::Sleep(_PTHREAD_FLUSH_TIME_);
  return NULL;
}

xstructure GetPrimitiveMULTITHREAD(const xstructure& _a,double tolerance) {  // APRIL 2009 JUNE 2012 added tolerance
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  if(LDEBUG) cerr << "GetPrimitiveMULTITHREAD: BEGIN " << endl;
  cout.setf(std::ios::fixed,std::ios::floatfield);
  cout.precision(10);
  xstructure a(_a);
  xstructure sstr=a;
  sstr.SetVolume(sstr.atoms.size());
  sstr=ReScale(sstr,1.0);
  sstr=BringInCell(sstr);
  //  sstr.CalculateSymmetry();
  if(tolerance<=0.0) a.equiv_fpos_epsilon=_EQUIV_FPOS_EPS_; else a.equiv_fpos_epsilon=tolerance;
  double sstr_volume=sstr.Volume();

  if(LDEBUG) cerr << "GetPrimitiveMULTITHREAD: [1] " << endl;

  _aflags aflags;
  // identify the minimum set of atoms
  // bool PGROUPWRITE=FALSE,PGROUPKWRITE=FALSE,FGROUPWRITE=FALSE,IATOMSWRITE=FALSE;
  // bool OSSWRITE=TRUE; // to FileMESSAGE, does not matter as it is /dev/null
  // ofstream FileMESSAGE("/dev/stderr");
  // SYM::CalculatePointGroup(FileMESSAGE,sstr,aflags,PGROUPWRITE,OSSWRITE,cout);
  // SYM::CalculatePointGroupKlattice(FileMESSAGE,sstr,aflags,PGROUPKWRITE,OSSWRITE,cout);
  // SYM::CalculateFactorGroup(FileMESSAGE,sstr,aflags,FGROUPWRITE,OSSWRITE,cout);
  // SYM::CalculateInequivalentAtoms(FileMESSAGE,sstr,aflags,IATOMSWRITE,OSSWRITE,cout);

  xmatrix<double> plattice(3,3),olattice(3,3);
  xvector<double> fdisp(3),cdisp(3);
  std::vector<xvector<double> > candidate_lattice_vector;

  if(LDEBUG) cerr << "GetPrimitiveMULTITHREAD: [2] " << endl;

  int specie_min=sstr.num_each_type.at(0),ispecie_min=0,specie_min_threshold=_PRIM_MULTITHREAD_MIN_ATOMS_THRESHOLD_; // seems a good threshold
  for(uint ispecie=0;ispecie<sstr.num_each_type.size();ispecie++)
    if(sstr.num_each_type.at(ispecie)<specie_min) {
      specie_min=sstr.num_each_type.at(ispecie);
      ispecie_min=ispecie;
    }
  // cerr << "DEBUG specie_min=" << specie_min << endl;

  // generate list of vectors
  candidate_lattice_vector.clear();
  candidate_lattice_vector.push_back(sstr.lattice(1));  // lattice is made of good vectors
  candidate_lattice_vector.push_back(sstr.lattice(2));  // lattice is made of good vectors
  candidate_lattice_vector.push_back(sstr.lattice(3));  // lattice is made of good vectors
  // no threads

  if(LDEBUG) cerr << "GetPrimitiveMULTITHREAD: [3] " << endl;
  
  if(!AFLOW_PTHREADS::FLAG || specie_min<=specie_min_threshold) {
    //   cerr << "NO PTHREADS" << endl;
    for(uint iat1=0;iat1<sstr.atoms.size();iat1++) { //      if(sstr.atoms.at(iat1).type==ispecie_min)
      for(uint iat2=0;iat2<sstr.atoms.size();iat2++) { //           if(sstr.atoms.at(iat2).type==ispecie_min)
	if(iat1!=iat2) {
	  fdisp=sstr.atoms.at(iat2).fpos-sstr.atoms.at(iat1).fpos;
	  cdisp=sstr.atoms.at(iat2).cpos-sstr.atoms.at(iat1).cpos;
	  if(aurostd::modulus(fdisp)>0.01 && aurostd::modulus(cdisp)>0.01)
	    if(IsTranslationFVector(sstr,fdisp)) 
	      candidate_lattice_vector.push_back(cdisp);
	}
      }
    }
  } else { // multithread
    cerr << "START THREADS [1] (GetPrimitiveMULTITHREAD=" << AFLOW_PTHREADS::MAX_PTHREADS << ") [" << specie_min << "]" << endl;
    AFLOW_PTHREADS::Clean_Threads();                                              // multithread clean
    _threaded_GETTVECTORS_params params[MAX_ALLOCATABLE_PTHREADS];                     // multithread
    // _threaded_GETTVECTORS_params params[AFLOW_PTHREADS::MAX_PTHREADS];               // multithread
    // vector<_threaded_GETTVECTORS_params> params(AFLOW_PTHREADS::MAX_PTHREADS);          // multithread
    // prepare
    for(int ithread=0;ithread<AFLOW_PTHREADS::MAX_PTHREADS;ithread++) {              // multithread
      // construction of params[i]                                         // multithread
      params[ithread].ITHREAD=ithread;                                     // multithread
      params[ithread].THREADS_MAX=AFLOW_PTHREADS::MAX_PTHREADS;                      // multithread
      params[ithread].pstr=&sstr;                                          // multithread
      params[ithread].ispecie_min=ispecie_min;                             // multithread
      params[ithread].ptvector=&candidate_lattice_vector;                  // multithread
      params[ithread].polattice=&olattice;                                 // multithread
      params[ithread].itbusy=ithread;                                      // multithread
      params[ithread].step=1;  // 1=gettvectors                            // multithread
    }                                                                      // multithread
    // run
    for(int ithread=0;ithread<AFLOW_PTHREADS::MAX_PTHREADS;ithread++) {
      AFLOW_PTHREADS::viret[ithread]=pthread_create(&(AFLOW_PTHREADS::vpthread[ithread]),NULL,_threaded_GetTvectors,(void*)&params[ithread]);
      //  aurostd::Sleep(10);
    }
    // collect
    for(int ithread=0;ithread<AFLOW_PTHREADS::MAX_PTHREADS;ithread++)
      pthread_join((AFLOW_PTHREADS::vpthread[ithread]),NULL);
  }
  if(LDEBUG) cout << "DEBUG: candidate_lattice_vector.size()=" << candidate_lattice_vector.size() << endl;

  // now split
  int cnt=0;
  olattice=sstr.lattice;                      // the lattice is always a good lattice
  if(!AFLOW_PTHREADS::FLAG || specie_min<=specie_min_threshold) {
    // now generate triplets
    for(uint iu=0;iu<candidate_lattice_vector.size();iu++) {
      for(uint i=1;i<=3;i++) plattice[1][i]=candidate_lattice_vector.at(iu)[i];
      for(uint iv=0;iv<candidate_lattice_vector.size()&& iv!=iu;iv++) {
        for(uint i=1;i<=3;i++) plattice[2][i]=candidate_lattice_vector.at(iv)[i];
        for(uint iw=0;iw<candidate_lattice_vector.size()&& iw!=iv && iw!=iu;iw++) {
          for(uint i=1;i<=3;i++) plattice[3][i]=candidate_lattice_vector.at(iw)[i];
          if(det(plattice)>0.999 && det(plattice)<sstr_volume) {   // no coplanar and contain at least 1 atom and smaller than the original cell
            if(aurostd::isinteger(sstr_volume/det(plattice),0.0001)) {    // integer ratio of volumes
	      //  cerr << sstr_volume/det(plattice) << " " << aurostd::isinteger(sstr_volume/det(plattice),0.001) << endl;
              if(det(plattice)<det(olattice)) {                    // better than before
                if(LDEBUG) cout<<"DEBUG"<<iu<<","<<iv<<","<<iw<<" "<< sstr_volume<<" "<<det(plattice)<<" "<<sstr_volume/det(plattice)<<endl;
                if(isdifferent(plattice,olattice,0.001)) {
                  plattice=MinkowskiBasisReduction(plattice);      // Minkowski first
                  plattice=NiggliUnitCellForm(plattice);           // Niggli Second
                  if(isdifferent(plattice,olattice,0.001)) {
                    olattice=plattice;
                    cnt++;
                  }
                }
              }
            }
          }
        }
      }
    }
  } else { // multithread
    //    if(LDEBUG)
    cerr << "START THREADS [2] (GetPrimitiveMULTITHREAD=" << AFLOW_PTHREADS::MAX_PTHREADS << ") [" << specie_min << "]" << endl;
    AFLOW_PTHREADS::Clean_Threads();                                              // multithread clean
    _threaded_GETTVECTORS_params params[MAX_ALLOCATABLE_PTHREADS];                     // multithread
    // _threaded_GETTVECTORS_params params[AFLOW_PTHREADS::MAX_PTHREADS];               // multithread
    // vector<_threaded_GETTVECTORS_params> params(AFLOW_PTHREADS::MAX_PTHREADS);          // multithread
    // prepare
    for(int ithread=0;ithread<AFLOW_PTHREADS::MAX_PTHREADS;ithread++) {              // multithread
      // construction of params[i]                                         // multithread
      params[ithread].ITHREAD=ithread;                                     // multithread
      params[ithread].THREADS_MAX=AFLOW_PTHREADS::MAX_PTHREADS;                      // multithread
      params[ithread].pstr=&sstr;                                          // multithread
      params[ithread].ispecie_min=ispecie_min;                             // multithread
      params[ithread].ptvector=&candidate_lattice_vector;                  // multithread
      params[ithread].polattice=&olattice;                                 // multithread
      params[ithread].itbusy=ithread;                                      // multithread
      params[ithread].step=2;  // 1=gettvectors                            // multithread
    }                                                                      // multithread
    // run
    for(int ithread=0;ithread<AFLOW_PTHREADS::MAX_PTHREADS;ithread++) {
      AFLOW_PTHREADS::viret[ithread]=pthread_create(&(AFLOW_PTHREADS::vpthread[ithread]),NULL,_threaded_GetTvectors,(void*)&params[ithread]);
      //  aurostd::Sleep(10);
    }
    // collect
    for(int ithread=0;ithread<AFLOW_PTHREADS::MAX_PTHREADS;ithread++)
      pthread_join((AFLOW_PTHREADS::vpthread[ithread]),NULL);
    if(LDEBUG) cerr << "GetPrimitiveMULTITHREAD: END THREADS " << endl;
  }

  plattice=olattice;

  // done

  xstructure b=sstr;
  b.lattice=plattice;//b.lattice=roundoff(b.lattice,_EPS_FPOS_EQUAL_);
  b.FixLattices();
  b.write_lattice_flag=FALSE;
  b.write_klattice_flag=FALSE;
  b.write_DEBUG_flag=FALSE;
  // plug them all
  for(uint iat=0;iat<sstr.atoms.size();iat++) {
    b.atoms.at(iat).fpos=BringInCell(C2F(b.lattice,b.atoms.at(iat).cpos));
    b.atoms.at(iat).cpos=F2C(b.lattice,b.atoms.at(iat).fpos);
  }
  // now remove them
  b=BringInCell(b);b.RemoveFractionalCopies(tolerance);
  b=BringInCell(b);b.RemoveCartesianCopies(0.01);
  // rescale back to original scale.
  b.SetVolume(Volume(a)*b.atoms.size()/a.atoms.size());
  b=ReScale(b,a.scale);
  //  // fix it up with the new Minkowsky and Niggli reductions // CANT DO AUTOMATICALLY
  //  b=LatticeReduction(b);
  // Put everything in new primitive cell.
  b=BringInCell(b);


  // check !

  // no fractional ratio of atoms
  double fraction_atoms=(double) a.atoms.size()/b.atoms.size();
  if(!aurostd::isinteger(fraction_atoms,0.01)) return _a;
  
  // nearest too close
  if(NearestNeighbour(b)<5.0*tolerance)  return _a;
  if(NearestNeighbour(b)<0.1)  return _a;

  // no messed up volume
  double fraction=Volume(a)/Volume(b);
  if(abs(b.atoms.size()*fraction-a.atoms.size())>0.1) {
    cerr << "ERROR   xstructure xstructure::GetPrimitive(void)" << endl;
    cerr << "        supercell has the wrong number of atoms" << endl;
    cerr << "        volume original    = " << Volume(a) << endl;
    cerr << "        volume prim        = " << Volume(b) << endl;
    cerr << "        a.scale            = " << a.scale << endl;
    cerr << "        b.scale            = " << b.scale << endl;
    cerr << "        a.atoms.size()     = " << a.atoms.size() << endl;
    cerr << "        b.atoms.size()     = " << b.atoms.size() << endl;
    cerr << "        fraction           = " << fraction << endl;
    cerr << "        supercell atoms    = " << fraction*b.atoms.size() << endl;
    cerr << b << endl;
    exit(0);
  }
  // everything ok
  if(LDEBUG) cerr << "GetPrimitiveMULTITHREAD: END [ok]=" << fraction_atoms << endl;
  return b;

}


xstructure GetPrimitiveSINGLE(const xstructure& _a,double tolerance) {  // APRIL 2009JUNE 2012 added tolerance
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  cout.setf(std::ios::fixed,std::ios::floatfield);
  cout.precision(10);
  xstructure a(_a);
  xstructure sstr=a;
  sstr.SetVolume(sstr.atoms.size());
  sstr=ReScale(sstr,1.0);
  sstr=BringInCell(sstr);
  if(tolerance<=0.0) a.equiv_fpos_epsilon=_EQUIV_FPOS_EPS_; else a.equiv_fpos_epsilon=tolerance;
  double sstr_volume=sstr.Volume();

  _aflags aflags;

  xmatrix<double> plattice(3,3),olattice(3,3);
  xvector<double> fdisp(3),cdisp(3);
  std::vector<xvector<double> > candidate_lattice_vector;

  int specie_min=100000,ispecie_min=0;
  for(uint ispecie=0;ispecie<sstr.num_each_type.size();ispecie++)
    if(sstr.num_each_type.at(ispecie)<specie_min) {
      specie_min=sstr.num_each_type.at(ispecie);
      ispecie_min=ispecie;
    }

  // generate list of vectors
  candidate_lattice_vector.push_back(sstr.lattice(1));  // lattice is made of good vectors
  candidate_lattice_vector.push_back(sstr.lattice(2));  // lattice is made of good vectors
  candidate_lattice_vector.push_back(sstr.lattice(3));  // lattice is made of good vectors
  for(uint iat1=0;iat1<sstr.atoms.size();iat1++)
    if(sstr.atoms.at(iat1).type==ispecie_min)
      for(uint iat2=0;iat2<sstr.atoms.size();iat2++)
	if(sstr.atoms.at(iat2).type==ispecie_min) {
	  fdisp=sstr.atoms.at(iat2).fpos-sstr.atoms.at(iat1).fpos;
	  cdisp=sstr.atoms.at(iat2).cpos-sstr.atoms.at(iat1).cpos;
	  if(aurostd::modulus(fdisp)>0.01 && aurostd::modulus(cdisp)>0.01)
	    if(IsTranslationFVector(sstr,fdisp))
	      candidate_lattice_vector.push_back(cdisp);
	}
  if(LDEBUG) cout << "DEBUG" << candidate_lattice_vector.size() << endl;
  int cnt=0;
  olattice=sstr.lattice;                      // the lattice is always a good lattice
  // now generate triplets
  for(uint iu=0;iu<candidate_lattice_vector.size();iu++) {
    for(uint i=1;i<=3;i++) plattice(1,i)=candidate_lattice_vector.at(iu)(i);
    for(uint iv=0;iv<candidate_lattice_vector.size()&& iv!=iu;iv++) {
      for(uint i=1;i<=3;i++) plattice(2,i)=candidate_lattice_vector.at(iv)(i);
      for(uint iw=0;iw<candidate_lattice_vector.size()&& iw!=iv && iw!=iu;iw++) {
	for(uint i=1;i<=3;i++) plattice(3,i)=candidate_lattice_vector.at(iw)(i);
	if(det(plattice)>0.999 && det(plattice)<sstr_volume) {   // no coplanar and contain at least 1 atom and smaller than the original cell
	  if(aurostd::isinteger(sstr_volume/det(plattice))) {    // integer ratio of volumes
	    if(det(plattice)<det(olattice)) {                    // better than before
	      if(LDEBUG) cout<<"DEBUG"<<iu<<","<<iv<<","<<iw<<" "<< sstr_volume<<" "<<det(plattice)<<" "<<sstr_volume/det(plattice)<<endl;
	      if(isdifferent(plattice,olattice,0.0001)) {
		plattice=MinkowskiBasisReduction(plattice);      // Minkowski first
		plattice=NiggliUnitCellForm(plattice);           // Niggli Second
		if(isdifferent(plattice,olattice,0.0001)) {
		  olattice=plattice;
		  cnt++;
		}
	      }
	    }
	  }
	}
      }
    }
  }
  plattice=olattice;
  // done

  xstructure b=sstr;
  b.lattice=plattice;//b.lattice=roundoff(b.lattice,_EPS_FPOS_EQUAL_);
  b.FixLattices();
  b.write_lattice_flag=FALSE;
  b.write_klattice_flag=FALSE;
  b.write_DEBUG_flag=FALSE;
  // plug them all
  for(uint iat=0;iat<sstr.atoms.size();iat++) {
    b.atoms.at(iat).fpos=BringInCell(C2F(b.lattice,b.atoms.at(iat).cpos));
    b.atoms.at(iat).cpos=F2C(b.lattice,b.atoms.at(iat).fpos);
  }
  // now remove them
  b.RemoveFractionalCopies();
  // rescale back to original scale.
  b.SetVolume(Volume(a)*b.atoms.size()/a.atoms.size());
  b=ReScale(b,a.scale);
  //  // fix it up with the new Minkowsky and Niggli reductions // CANT DO AUTOMATICALLY
  //  b=LatticeReduction(b);
  // Put everything in new primitive cell.
  b=BringInCell(b);
  // check !
  double fraction=Volume(a)/Volume(b);
  if(abs(b.atoms.size()*fraction-a.atoms.size())>0.1) {
    cerr << "ERROR   xstructure xstructure::GetPrimitive(void)" << endl;
    cerr << "        supercell has the wrong number of atoms" << endl;
    cerr << "        volume original    = " << Volume(a) << endl;
    cerr << "        volume prim        = " << Volume(b) << endl;
    cerr << "        a.scale            = " << a.scale << endl;
    cerr << "        b.scale            = " << b.scale << endl;
    cerr << "        a.atoms.size()     = " << a.atoms.size() << endl;
    cerr << "        b.atoms.size()     = " << b.atoms.size() << endl;
    cerr << "        fraction           = " << fraction << endl;
    cerr << "        supercell atoms    = " << fraction*b.atoms.size() << endl;
    cerr << b << endl;
    exit(0);
  }
  // everything ok
  return b;
}


xstructure GetPrimitive1(const xstructure& a) {  // MARCH 2009
  if(a.equiv_fpos_epsilon<1.0e-12) { cerr << "ERROR:  Zero tolerance in aflow_xatom.cpp: GetPrimitive1" << endl; exit(0);}
  double tolerance=a.equiv_fpos_epsilon;

  cout.setf(std::ios::fixed,std::ios::floatfield);
  cout.precision(10);
  xstructure sstr=a;
  sstr.SetVolume(sstr.atoms.size());
  sstr=ReScale(sstr,1.0);
  sstr=BringInCell(sstr);
  double sstr_volume=sstr.Volume();
  _aflags aflags;
  // identify the minimum set of atoms

  xmatrix<double> plattice(3,3),olattice(3,3);
  xvector<double> fdisp(3),cdisp(3);
  std::vector<xvector<double> > candidate_lattice_vector;
  
  int specie_min=100000,ispecie_min=0;
  for(uint ispecie=0;ispecie<sstr.num_each_type.size();ispecie++)
    if(sstr.num_each_type.at(ispecie)<specie_min) {
      specie_min=sstr.num_each_type.at(ispecie);
      ispecie_min=ispecie;
    }

  // generate list of vectors
  candidate_lattice_vector.push_back(sstr.lattice(1));  // lattice is made of good vectors
  candidate_lattice_vector.push_back(sstr.lattice(2));  // lattice is made of good vectors
  candidate_lattice_vector.push_back(sstr.lattice(3));  // lattice is made of good vectors
  olattice=sstr.lattice;                      // the lattice is always a good lattice
  for(uint iat1=0;iat1<sstr.atoms.size();iat1++) {
    for(uint iat2=0;iat2<sstr.atoms.size();iat2++) {
      if(iat2!=iat1 &&
         sstr.atoms.at(iat1).type==ispecie_min &&
         sstr.atoms.at(iat2).type==ispecie_min) {
        fdisp=sstr.atoms.at(iat2).fpos-sstr.atoms.at(iat1).fpos;
        cdisp=sstr.atoms.at(iat2).cpos-sstr.atoms.at(iat1).cpos;
        if(IsTranslationFVector(sstr,fdisp)) candidate_lattice_vector.push_back(cdisp);
      }
    }
  }
  //  cerr << candidate_lattice_vector.size() << endl;
  int cnt=0;
  // now generate triplets
  for(uint iu=0;iu<candidate_lattice_vector.size();iu++) {
    for(uint i=1;i<=3;i++) plattice(1,i)=candidate_lattice_vector.at(iu)(i);
    for(uint iv=0;iv<candidate_lattice_vector.size()&& iv!=iu;iv++) {
      for(uint i=1;i<=3;i++) plattice(2,i)=candidate_lattice_vector.at(iv)(i);
      for(uint iw=0;iw<candidate_lattice_vector.size()&& iw!=iv && iw!=iu;iw++) {
        for(uint i=1;i<=3;i++) plattice(3,i)=candidate_lattice_vector.at(iw)(i);
        if(det(plattice)>tolerance && det(plattice)<sstr_volume && det(plattice)<det(olattice)) { // well defined
          if(isdifferent(plattice,olattice,0.0001)) {
            plattice=MinkowskiBasisReduction(plattice); // Minkowski first
            plattice=NiggliUnitCellForm(plattice);      // Niggli Second
            if(isdifferent(plattice,olattice,0.0001)) {
              olattice=plattice;
              //      cerr << det(olattice) << " " << cnt<< endl;
              cnt++;
            }
          }
        }
      }
    }
  }
  plattice=olattice;
  // done

  _atom atom;
  xstructure b=sstr;
  b.atoms.clear();
  b.lattice=plattice;b.lattice=roundoff(b.lattice,tolerance);
  b.FixLattices();
  b.write_lattice_flag=FALSE;
  b.write_klattice_flag=FALSE;
  b.write_DEBUG_flag=FALSE;
  bool atom_found=FALSE;
  // for scanning around
  //  double radius=aurostd::modulus(sstr.lattice(1))+aurostd::modulus(sstr.lattice(2))+aurostd::modulus(sstr.lattice(3));
  //  int dims=max(LatticeDimensionSphere(sstr.lattice,radius));
  int dims=max(LatticeDimensionSphere(sstr.lattice,1.5*RadiusSphereLattice(sstr.lattice)));
  for(uint iat=0;iat<b.num_each_type.size();iat++) {
    b.num_each_type.at(iat)=0; // create enough space
    b.comp_each_type.at(iat)=0; // create enough space
  }

  for(uint iat=0;iat<sstr.atoms.size();iat++) {
    for(int i=-dims;i<=dims;i++)
      for(int j=-dims;j<=dims;j++)
        for(int k=-dims;k<=dims;k++) {
          //    atom=BringInCell(sstr.atoms.at(i),sstr.lattice);
          atom=sstr.atoms.at(iat);
          atom.cpos=atom.cpos+((double)i)*sstr.lattice(1)+((double)j)*sstr.lattice(2)+((double)k)*sstr.lattice(3);
          atom.fpos=C2F(b.lattice,atom.cpos);
          if(atom.fpos(1)>=-tolerance && atom.fpos(1)<1.0-tolerance &&
             atom.fpos(2)>=-tolerance && atom.fpos(2)<1.0-tolerance &&
             atom.fpos(3)>=-tolerance && atom.fpos(3)<1.0-tolerance) {      // found something inside
            for(uint ii=0;ii<b.atoms.size()&&!atom_found;ii++)
              atom_found=identical(atom.cpos,b.atoms.at(ii).cpos,0.1);       // look in all the list of operations
            // atom_found=FALSE;
            if(!atom_found) {
              atom.fpos=roundoff(atom.fpos,tolerance);
              atom.cpos=roundoff(atom.cpos,tolerance);
              b.atoms.push_back(atom);
              b.num_each_type.at(atom.type)++;         // CONVASP_MODE
              b.comp_each_type.at(atom.type)+=atom.partial_occupation_value;         // CONVASP_MODE
            }
          }
        }
  }
  b.GetStoich();  //170724 - CO
  // rescale back to original scale.
  b.SetVolume(Volume(a)*b.atoms.size()/a.atoms.size());
  b=ReScale(b,a.scale);
  // fix it up with the new Minkowsky and Niggli reductions
  b=LatticeReduction(b);
  // Put everything in new primitive cell.
  b=BringInCell(b);
  // check !
  double fraction=Volume(a)/Volume(b);
  if(abs(b.atoms.size()*fraction-a.atoms.size())>0.1) {
    cerr << "ERROR   xstructure xstructure::GetPrimitive(void)" << endl;
    cerr << "        supercell has the wrong number of atoms" << endl;
    cerr << "        volume original    = " << Volume(a) << endl;
    cerr << "        volume prim        = " << Volume(b) << endl;
    cerr << "        a.scale            = " << a.scale << endl;
    cerr << "        b.scale            = " << b.scale << endl;
    cerr << "        a.atoms.size()     = " << a.atoms.size() << endl;
    cerr << "        b.atoms.size()     = " << b.atoms.size() << endl;
    cerr << "        fraction           = " << fraction << endl;
    cerr << "        supercell atoms    = " << fraction*b.atoms.size() << endl;
    cerr << b << endl;
    exit(0);
  }
  // everything ok
  return b;
}



// second try
xstructure GetPrimitive2(const xstructure& a) {
  if(a.equiv_fpos_epsilon<1.0e-12) { cerr << "ERROR:  Zero tolerance in aflow_xatom.cpp: GetPrimitive2" << endl; exit(0);}
  double tolerance=a.equiv_fpos_epsilon;

  cout.setf(std::ios::fixed,std::ios::floatfield);
  cout.precision(10);

  xstructure sstr=a;
  //  sstr.lattice=lattice;
  // sstr.scale=scale;
  std::vector<xvector<double> > candidate_lattice_vector;
  xmatrix<double> plattice(3,3);
  vector<xmatrix<double> > plattice_list;
  xvector<double> fdisp(3),cdisp(3);

  // Get all the data from the structure for easy use.
  // Set scale to 1 so you don't need to rescale coordinates.
  sstr=ReScale(sstr,1.0);
  // Put everything in the unit cell.
  //  cerr << sstr.scale << endl;
  sstr=BringInCell(sstr);
  // cerr << sstr.scale << endl;
  // exit(0);
  string title=sstr.title;
  int i;
  for(uint iat=0;iat<sstr.num_each_type.size();iat++) {
    sstr.num_each_type.at(iat)=0; // create enough space
    sstr.comp_each_type.at(iat)=0; // create enough space
  }
  // Create list of candidate plvs.

  candidate_lattice_vector.push_back(sstr.lattice(1));
  candidate_lattice_vector.push_back(sstr.lattice(2));
  candidate_lattice_vector.push_back(sstr.lattice(3));

  for(uint ia=1;ia<sstr.atoms.size();ia++) {
    // Here we calculate all displacements from the first atom, i.e., the
    // first atom is considered the origin.
    cdisp=sstr.atoms.at(ia).cpos-sstr.atoms.at(0).cpos;
    if(IsTranslationCVector(sstr,cdisp)) candidate_lattice_vector.push_back(cdisp);
  }
  // Now the candidate_lattice_vector have been found and we must take all possible
  // traids of *distinct* vectors and get the smallest volume.
  int num_cand=candidate_lattice_vector.size();
  double min_vol,new_vol;
  min_vol=GetVol(sstr.lattice);
  for(int iu=0;iu<num_cand;iu++) {
    for(i=1;i<=3;i++) plattice(1,i)=candidate_lattice_vector.at(iu)(i);
    for(int iv=0;iv<num_cand;iv++) {
      for(i=1;i<=3;i++) plattice(2,i)=candidate_lattice_vector.at(iv)(i);
      for(int iw=0;iw<num_cand;iw++) {
	if(iu!=iw && iu!=iv && iv!=iw) {
	  for(i=1;i<=3;i++) plattice(3,i)=candidate_lattice_vector.at(iw)(i);
	  new_vol=GetVol(plattice);
	  // if(new_vol<=min_vol && new_vol>tolerance) min_vol=new_vol;
	  if(new_vol<=min_vol && new_vol>1.0e-5) min_vol=new_vol;
	}
      } //iw loop
    } //iv loop
  } //iu loop
  // Now that we know the min_vol we must go back through all the triads
  // and store the set that have the min_vol.
  for(int iu=0;iu<num_cand;iu++) {
    for(i=1;i<=3;i++) plattice(1,i)=candidate_lattice_vector.at(iu)(i);
    for(int iv=0;iv<num_cand;iv++) {
      for(i=1;i<=3;i++) plattice(2,i)=candidate_lattice_vector.at(iv)(i);
      for(int iw=0;iw<num_cand;iw++) {
        if(iu!=iw && iu!=iv && iv!=iw) {
	  for(i=1;i<=3;i++) plattice(3,i)=candidate_lattice_vector.at(iw)(i);
          new_vol=GetVol(plattice);
          if(fabs(new_vol-min_vol)<tolerance && new_vol>1.0e-5)
	    plattice_list.push_back(plattice); // Add to plattice_list
	  //if new_vol=min_vol
        } //if iu!=iv!=iw
      } //iw loop
    } //iv loop
  } //iu loop
  // Now we have a set of all possible primitive lattices
  // and we want to pick one that looks as much like the
  // original lattice as possible.
  double paw,pawopt=1e6;
  for(uint i=0;i<plattice_list.size();i++) {
    paw=0;
    for(int j=1;j<=3;j++) paw+=angle(plattice_list.at(i)(j),sstr.lattice(j));
    if(paw<pawopt) {pawopt=paw;plattice=plattice_list.at(i);}; // minimization
  }
  // For now just pick the first one.
  // if(plattice_list.size()>0) plattice=plattice_list.at(1);
  // If volume did not reduce then just keep original lattice
  if(fabs(GetVol(sstr.lattice)-GetVol(plattice))<tolerance) {plattice=sstr.lattice;}

  // Now we have the plvs for the new structure.  Now we must
  // construct a new structure object with the plvs and return it.
  // This new structure object will be just like the old one but
  // with new lattice vectors and new basis atoms (which will be a
  // subset of the original ones, chosen to have direct coordinates
  // less than 1).

  // To create new structure we must set all the parameters.
  // To get new basis atoms, loop over all the original basis atoms,
  // get new direct coordinates, make all less than one, then take unique
  // values.

  _atom atom;
  xstructure b=sstr;
  b.atoms.clear();
  b.lattice=plattice;b.lattice=roundoff(b.lattice,tolerance);
  b.FixLattices();
  b.write_lattice_flag=TRUE;
  b.write_klattice_flag=FALSE;
  b.write_DEBUG_flag=FALSE;
  bool atom_found=FALSE;
  // for scanning around
  double radius;
  radius=aurostd::modulus(sstr.lattice(1))+aurostd::modulus(sstr.lattice(2))+aurostd::modulus(sstr.lattice(3));
  int dims=max(LatticeDimensionSphere(sstr.lattice,radius));
  for(uint iat=0;iat<b.num_each_type.size();iat++) {
    b.num_each_type.at(iat)=0; // create enough space
    b.comp_each_type.at(iat)=0; // create enough space
  }

  for(uint iat=0;iat<sstr.atoms.size();iat++) {
    for(int i=-dims;i<=dims;i++)
      for(int j=-dims;j<=dims;j++)
	for(int k=-dims;k<=dims;k++) {
          //    atom=BringInCell(sstr.atoms.at(i),sstr.lattice);
	  atom=sstr.atoms.at(iat);
	  atom.cpos=atom.cpos+((double)i)*sstr.lattice(1)+((double)j)*sstr.lattice(2)+((double)k)*sstr.lattice(3);
	  atom.fpos=C2F(b.lattice,atom.cpos);
	  if(atom.fpos(1)>=-tolerance && atom.fpos(1)<1.0-tolerance &&
	     atom.fpos(2)>=-tolerance && atom.fpos(2)<1.0-tolerance &&
	     atom.fpos(3)>=-tolerance && atom.fpos(3)<1.0-tolerance) {      // found something inside
	    for(uint ii=0;ii<b.atoms.size()&&!atom_found;ii++)
	      atom_found=identical(atom.cpos,b.atoms.at(ii).cpos,0.1);       // look in all the list of operations
	    // atom_found=FALSE;
	    if(!atom_found) {
	      atom.fpos=roundoff(atom.fpos,tolerance);
	      atom.cpos=roundoff(atom.cpos,tolerance);
	      b.atoms.push_back(atom);
	      b.num_each_type.at(atom.type)++;         // CONVASP_MODE
	      b.comp_each_type.at(atom.type)+=atom.partial_occupation_value;         // CONVASP_MODE
	    }
	  }
	}
  }
  b.GetStoich();  //170724 CO
  // some check
  double fraction=((double) GetVol(sstr.lattice)/GetVol(b.lattice));
  if(abs(b.atoms.size()*fraction-a.atoms.size())>0.1) {
    cerr << "ERROR   xstructure xstructure::GetPrimitive(void)" << endl;
    cerr << "        supercell has the wrong number of atoms" << endl;
    cerr << "        volume original    = " << GetVol(sstr.lattice) << endl;
    cerr << "        volume prim        = " << GetVol(b.lattice) << endl;
    //cerr << "        sstr.scale            = " << sstr.scale << endl;
    //cerr << "        b.scale            = " << b.scale << endl;
    cerr << "        sstr.atoms.size()     = " << a.atoms.size() << endl;
    cerr << "        b.atoms.size()     = " << b.atoms.size() << endl;
    cerr << "        fraction           = " << fraction << endl;
    cerr << "        supercell atoms    = " << fraction*b.atoms.size() << endl;
    cerr << GetVol(b.lattice) << endl;
    cerr << GetVol(a.lattice) << endl;
    //   exit(0);
  }
  // rescale back to original scale.
  b=ReScale(b,a.scale);
  // fix it up with the new Minkowsky and Niggli reductions
  b=LatticeReduction(b);
  // Put everything in new primitive cell.
  b=BringInCell(b);
  // everything ok
  return b;
}

// third try
xstructure GetPrimitive3(const xstructure& a) {
  if(a.equiv_fpos_epsilon<1.0e-12) { cerr << "ERROR:  Zero tolerance in aflow_xatom.cpp: GetPrimitive3" << endl; exit(0);}
  double tolerance=a.equiv_fpos_epsilon;

  cout.setf(std::ios::fixed,std::ios::floatfield);
  cout.precision(10);
  xstructure sstr=a;
  sstr.SetVolume(sstr.atoms.size());
  sstr=ReScale(sstr,1.0);
  sstr=BringInCell(sstr);

  double sstr_volume=sstr.Volume();
  _aflags aflags;aflags.Directory="./";
  // identify the minimum set of atoms
  bool PGROUPWRITE=TRUE,FGROUPWRITE=TRUE;//,IATOMSWRITE=TRUE;
  bool OSSWRITE=FALSE; // to FileMESSAGE, does not matter as it is /dev/null
  ofstream FileMESSAGE("/dev/null");
  bool _QUIET_=XHOST.QUIET;XHOST.QUIET=TRUE;
  SYM::CalculatePointGroup(FileMESSAGE,sstr,aflags,PGROUPWRITE,OSSWRITE,cout);
  SYM::CalculateFactorGroup(FileMESSAGE,sstr,aflags,FGROUPWRITE,OSSWRITE,cout);
  //  SYM::CalculateInequivalentAtoms(FileMESSAGE,sstr,aflags,IATOMSWRITE,OSSWRITE,cout);

  xmatrix<double> plattice(3,3),olattice(3,3);
  xvector<double> fdisp(3),cdisp(3);
  std::vector<xvector<double> > candidate_lattice_vector;

  // generate list of vectors
  candidate_lattice_vector.push_back(sstr.lattice(1));  // lattice is made of good vectors
  candidate_lattice_vector.push_back(sstr.lattice(2));  // lattice is made of good vectors
  candidate_lattice_vector.push_back(sstr.lattice(3));  // lattice is made of good vectors
  olattice=sstr.lattice;                      // the lattice is always a good lattice

  bool sym_found;
  for(uint i=0;i<sstr.fgroup.size();i++) {
    if(aurostd::modulus(sstr.fgroup.at(i).ctau)>0.01) {
      //    cerr << i << " " << sstr.fgroup.at(i).ctau << endl;
      sym_found=FALSE;
      for(uint ii=0;ii<candidate_lattice_vector.size()&&!sym_found;ii++)
	sym_found=identical(sstr.fgroup.at(i).ctau,candidate_lattice_vector[ii],tolerance);    // look in all the list of operations
      if(sym_found==FALSE) {                                          // new operation, generate and save it
	candidate_lattice_vector.push_back(sstr.fgroup.at(i).ctau);
	cerr << i << " " << sstr.fgroup.at(i).ctau << endl;
      }
    }
  }

  cerr << candidate_lattice_vector.size() << endl;
  int cnt=0;
  // now generate triplets
  for(uint iu=0;iu<candidate_lattice_vector.size();iu++) {
    for(uint i=1;i<=3;i++) plattice(1,i)=candidate_lattice_vector.at(iu)(i);
    for(uint iv=0;iv<candidate_lattice_vector.size()&& iv!=iu;iv++) {
      for(uint i=1;i<=3;i++) plattice(2,i)=candidate_lattice_vector.at(iv)(i);
      for(uint iw=0;iw<candidate_lattice_vector.size()&& iw!=iv && iw!=iu;iw++) {
	for(uint i=1;i<=3;i++) plattice(3,i)=candidate_lattice_vector.at(iw)(i);
	if(det(plattice)>tolerance && det(plattice)<sstr_volume && det(plattice)<det(olattice)) { // well defined
	  if(isdifferent(plattice,olattice,0.0001)) {
	    plattice=MinkowskiBasisReduction(plattice); // Minkowski first
	    plattice=NiggliUnitCellForm(plattice);      // Niggli Second
	    if(isdifferent(plattice,olattice,0.0001)) {
	      olattice=plattice;
	      //      cerr << det(olattice) << " " << cnt<< endl;
	      cnt++;
	    }
	  }
	}
      }
    }
  }
  plattice=olattice;
  // done

  _atom atom;
  xstructure b=sstr;
  b.atoms.clear();
  b.lattice=plattice;b.lattice=roundoff(b.lattice,tolerance);
  b.FixLattices();
  b.write_lattice_flag=FALSE;
  b.write_klattice_flag=FALSE;
  b.write_DEBUG_flag=FALSE;
  bool atom_found=FALSE;
  // for scanning around
  int dims=max(LatticeDimensionSphere(sstr.lattice,RadiusSphereLattice(sstr.lattice)));
  for(uint iat=0;iat<b.num_each_type.size();iat++) {
    b.num_each_type.at(iat)=0; // create enough space
    b.comp_each_type.at(iat)=0; // create enough space
  }

  for(uint iat=0;iat<sstr.atoms.size();iat++) {
    for(int i=-dims;i<=dims;i++)
      for(int j=-dims;j<=dims;j++)
	for(int k=-dims;k<=dims;k++) {
          //    atom=BringInCell(sstr.atoms.at(i),sstr.lattice);
	  atom=sstr.atoms.at(iat);
	  atom.cpos=atom.cpos+((double)i)*sstr.lattice(1)+((double)j)*sstr.lattice(2)+((double)k)*sstr.lattice(3);
	  atom.fpos=C2F(b.lattice,atom.cpos);
	  if(atom.fpos(1)>=-tolerance && atom.fpos(1)<1.0-tolerance &&
	     atom.fpos(2)>=-tolerance && atom.fpos(2)<1.0-tolerance &&
	     atom.fpos(3)>=-tolerance && atom.fpos(3)<1.0-tolerance) {      // found something inside
	    for(uint ii=0;ii<b.atoms.size()&&!atom_found;ii++)
	      atom_found=identical(atom.cpos,b.atoms.at(ii).cpos,0.1);       // look in all the list of operations
	    // atom_found=FALSE;
	    if(!atom_found) {
	      atom.fpos=roundoff(atom.fpos,tolerance);
	      atom.cpos=roundoff(atom.cpos,tolerance);
	      b.atoms.push_back(atom);
	      b.num_each_type.at(atom.type)++;         // CONVASP_MODE
	      b.comp_each_type.at(atom.type)+=atom.partial_occupation_value;         // CONVASP_MODE
	    }
	  }
	}
  }
  b.GetStoich();  //170724 CO
  // rescale back to original scale.
  b.SetVolume(Volume(a)*b.atoms.size()/a.atoms.size());
  b=ReScale(b,a.scale);
  // fix it up with the new Minkowsky and Niggli reductions
  b=LatticeReduction(b);
  XHOST.QUIET=_QUIET_;

  // Put everything in new primitive cell.
  b=BringInCell(b);
  // check !
  double fraction=Volume(a)/Volume(b);
  if(abs(b.atoms.size()*fraction-a.atoms.size())>0.1) {
    cerr << "ERROR   xstructure xstructure::GetPrimitive(void)" << endl;
    cerr << "        supercell has the wrong number of atoms" << endl;
    cerr << "        volume original    = " << Volume(a) << endl;
    cerr << "        volume prim        = " << Volume(b) << endl;
    cerr << "        a.scale            = " << a.scale << endl;
    cerr << "        b.scale            = " << b.scale << endl;
    cerr << "        a.atoms.size()     = " << a.atoms.size() << endl;
    cerr << "        b.atoms.size()     = " << b.atoms.size() << endl;
    cerr << "        fraction           = " << fraction << endl;
    cerr << "        supercell atoms    = " << fraction*b.atoms.size() << endl;
    cerr << b << endl;
    exit(0);
  }
  // everything ok
  return b;
}




// ***************************************************************************
// Operator GetPrimitive
// ***************************************************************************

void xstructure::GetPrimitive(void) {
  extern xstructure GetPrimitive(const xstructure& a); // so it does not recurse
  xstructure a(*this);
  *this=GetPrimitive(a);
}

void xstructure::GetPrimitive(double tolerance) {
  extern xstructure GetPrimitive(const xstructure& a,double tolerance); // so it does not recurse
  xstructure a(*this);
  *this=GetPrimitive(a,tolerance);
}

void xstructure::GetPrimitive2(void) {
  extern xstructure GetPrimitive2(const xstructure& a); // so it does not recurse
  xstructure a(*this);
  *this=GetPrimitive2(a);
}

void xstructure::GetPrimitive3(void) {
  extern xstructure GetPrimitive3(const xstructure& a); // so it does not recurse
  xstructure a(*this);
  *this=GetPrimitive3(a);
}

// ***************************************************************************
// Function MinDist
// ***************************************************************************
double xstructure::MinDist(void) {
  dist_nn_min=SYM::minimumDistance(*this);
  return dist_nn_min;
}

// ***************************************************************************
// Function ReScale
// ***************************************************************************
xstructure ReScale(const xstructure& a, const double &in_scale) {
  // This resets scale and changes the cell parameters and coordinates
  // appropriately.  Keeps volume fixed.
  if(in_scale==0.0) {cerr << "structure::ReScale in_scale must be non zero" << endl;exit(0);}
  xstructure b(a);
  if(aurostd::identical(b.scale,in_scale,_ZERO_TOL_)){return b;}  //try hard not to introduce precision errors, currently we print scale with precision 6
  b.lattice=b.lattice*b.scale/in_scale;
  b.origin=b.origin*b.scale/in_scale;
  b.f2c=trasp(b.lattice);
  b.c2f=inverse(trasp(b.lattice));
  // klattice already contained the scale so it does not need to be fixed
  // b.klattice=b.klattice*in_scale/b.scale;
  for(int i=0;i<(int)b.atoms.size();i++) {
    b.atoms.at(i).fpos=a.atoms.at(i).fpos;
    b.atoms.at(i).cpos=a.atoms.at(i).cpos*b.scale/in_scale;
  }
  if(b.fgroup_calculated) {
    for(int fg=0;fg<(int)b.fgroup.size();fg++) {
      b.fgroup[fg].ctau=b.fgroup[fg].ctau*b.scale/in_scale;
    }
  }
  if(b.sgroup_calculated) {
    for(int sg=0;sg<(int)b.sgroup.size();sg++) {
      b.sgroup[sg].ctau=b.sgroup[sg].ctau*b.scale/in_scale;
      b.sgroup[sg].ctrasl=b.sgroup[sg].ctrasl*b.scale/in_scale;
    }
  }
  b.scale=in_scale;
  b.FixLattices();  // touched scale, then fix the lattices
  return b;
}

void xstructure::ReScale(const double &in_scale) {
  if(in_scale==0.0) {cerr << "structure::ReScale in_scale must be non zero" << endl;exit(0);}
  if(aurostd::identical(scale,in_scale,_ZERO_TOL_)){return;}  //try hard not to introduce precision errors, currently we print scale with precision 6
  lattice=lattice*scale/in_scale;
  origin=origin*scale/in_scale;
  f2c=trasp(lattice);
  c2f=inverse(trasp(lattice));
  // klattice already contained the scale so it does not need to be fixed
  // klattice=klattice*in_scale/scale;
  for(int i=0;i<(int)atoms.size();i++){
    atoms.at(i).cpos=atoms.at(i).cpos*scale/in_scale;
  }
  if(fgroup_calculated) {
    for(int fg=0;fg<(int)fgroup.size();fg++) {
      fgroup[fg].ctau=fgroup[fg].ctau*scale/in_scale;
    }
  }
  if(sgroup_calculated) {
    for(int sg=0;sg<(int)sgroup.size();sg++) {
      sgroup[sg].ctau=sgroup[sg].ctau*scale/in_scale;
      sgroup[sg].ctrasl=sgroup[sg].ctrasl*scale/in_scale;
    }
  }
  scale=in_scale;
  FixLattices(); // touched scale, then fix the lattices
}

// ***************************************************************************
// Function SetScale
// ***************************************************************************
void xstructure::SetScale(const double &in_scale) {
  if(in_scale==0.0) { cerr << _AUROSTD_XLIBS_ERROR_ << "structure::SetScale in_scale must be non zero" << endl;exit(0);}
  scale=in_scale;
  FixLattices();  // touched scale, then fix the lattices
}

xstructure SetScale(const xstructure& a,const double &in_scale) {
  // This resets scale.  Keeps volume fixed.
  if(in_scale==0.0) { cerr << _AUROSTD_XLIBS_ERROR_ << "structure::SetScale in_scale must be non zero" << endl;exit(0);}
  xstructure b;b=a;
  b.scale=in_scale;
  b.FixLattices();  // touched scale, then fix the lattices
  return b;
}

// ***************************************************************************
// Function SetVolume
// ***************************************************************************
void xstructure::SetVolume(const double &in_volume) {
  if(in_volume==0.0) { cerr << _AUROSTD_XLIBS_ERROR_ << "structure::SetVolume in_scale must be non zero" << endl;exit(0);}
  scale=std::pow((double) in_volume/det(lattice),(double) 1/3);
  FixLattices();  // touched scale, then fix the lattices
}

xstructure SetVolume(const xstructure& a,const double &in_volume) {
  if(in_volume==0.0) { cerr << _AUROSTD_XLIBS_ERROR_ << "structure::SetVolume in_volume must be non zero" << endl;exit(0);}
  xstructure b;b=a;
  b.scale=std::pow((double) in_volume/det(b.lattice),(double) 1/3);
  b.FixLattices(); // touched scale, need to fix the lattices
  return b;
}

// ***************************************************************************
// Function InflateLattice
// ***************************************************************************
void xstructure::InflateLattice(const double &coefficient) {
  if(coefficient==0.0) { cerr << _AUROSTD_XLIBS_ERROR_ << "structure::InflateLattice coefficient must be non zero" << endl;exit(0);}
  //  scale=coefficient*scale;
  lattice=coefficient*lattice;
  FixLattices();  // touched scale/lattice, then fix the lattices
}

xstructure InflateLattice(const xstructure& a,const double &coefficient) {
  // This resets scale.  Keeps volume fixed.
  if(coefficient==0.0) { cerr << _AUROSTD_XLIBS_ERROR_ << "structure::InflateLattice coefficient must be non zero" << endl;exit(0);}
  xstructure b;b=a;
  //  b.scale=coefficient*b.scale;
  b.lattice=coefficient*b.lattice;
  b.FixLattices();  // touched scale/lattice, then fix the lattices
  return b;
}

// ***************************************************************************
// Function InflateVolume
// ***************************************************************************
void xstructure::InflateVolume(const double &coefficient) {
  if(coefficient==0.0) { cerr << _AUROSTD_XLIBS_ERROR_ << "structure::InflateVolume in_scale must be non zero" << endl;exit(0);}
  // scale=std::pow((double) coefficient,(double) 1/3)*scale;
  lattice=std::pow((double) coefficient,(double) 1/3)*lattice;
  FixLattices();  // touched scale/lattice, then fix the lattices
}

xstructure InflateVolume(const xstructure& a,const double &coefficient) {
  if(coefficient==0.0) { cerr << _AUROSTD_XLIBS_ERROR_ << "structure::InflateVolume coefficient must be non zero" << endl;exit(0);}
  xstructure b;b=a;
  //  b.scale=std::pow((double) coefficient,(double) 1/3)*b.scale;
  b.lattice=std::pow((double) coefficient,(double) 1/3)*b.lattice;
  b.FixLattices(); // touched scale/lattice, need to fix the lattices
  return b;
}


// ***************************************************************************
// Function GetVolume
// ***************************************************************************
double xstructure::GetVolume(void) {
  return scale*scale*scale*det(lattice);
}

double GetVolume(const xstructure& a) {
  return a.scale*a.scale*a.scale*det(a.lattice);
}

// ***************************************************************************
// Function Volume
// ***************************************************************************
double xstructure::Volume(void) {
  return scale*scale*scale*det(lattice);
}

double Volume(const xstructure& a) {
  return a.scale*a.scale*a.scale*det(a.lattice);
}

//DX - START
bool inCell(xvector<double>& pos_vec){
  if(pos_vec(1)>=-_ZERO_TOL_ && pos_vec(1)<1.0-_ZERO_TOL_ &&
     pos_vec(2)>=-_ZERO_TOL_ && pos_vec(2)<1.0-_ZERO_TOL_ &&
     pos_vec(3)>=-_ZERO_TOL_ && pos_vec(3)<1.0-_ZERO_TOL_) {      // found something inside
    return TRUE;
  }
  return FALSE;
}
//DX - END

// ***************************************************************************
// Function GetSuperCell
// ***************************************************************************
// This funtion creates SuperCells with 9 or 3 elements...
// the old routine by Dane Morgan does not work, so I rewrote this
// one from scratch.
// The algorithm is simple: make the bigger cell
// b.lattice=supercell*a.lattice;
// and generate a bunch of atoms around the old cell and check
// if they are inside the new cell... There is a checksum with ERROR
// if the numbers do not match. Stefano Curtarolo (aug07).
// COREY (2017) adding maps to/from supercell/primitive structure, get_symmetry and get_full_basis flags
// get_symmetry will propagate symmetry of primitive cell
// WARNING: if the structure is a derivative structure (i.e., non-uniform expansion, see derivative_structure) then
// NOT all symmetry operations will work, as the symmetry is reduced
// YOU need to check which of these don't work, unless you also calculate get_full_basis,
// which checks the validity of the symmetry operation by finding the full basis map

// COREY - START
xstructure GetSuperCell(const xstructure& aa, const xmatrix<double> &supercell,vector<int>& sc2pcMap,vector<int>& pc2scMap,
                        bool get_symmetry, bool get_full_basis) {
  //xstructure GetSuperCell(const xstructure& aa, const xmatrix<double> &supercell) {
  //#define _eps_scell_ 1.0e-5
  //#define _eps_scell_ 0.0
  // check for error
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  double vol_supercell=det(supercell);
  if(abs(vol_supercell)<0.001) {cerr << "ERROR: singular supercell matrix" << endl;exit(0);}
  xstructure a(aa); a.ReScale(1.0); //the nuclear option, the only way not to mess around with scale EVERYWHERE
  //DO NOT MODIFY STRUCTURE IN HERE, WE WANT TO PROPAGATE SYMMETRY FROM PRIMITIVE STRUCTURE!
  //a.BringInCell();
  xstructure b(a); //actual supercell, need to copy!
  //pflow::CalculateFullSymmetry(a);
  _atom atom,atom2;
  _sym_op pSymOp,fSymOp,aSymOp;  //corey
  int i,j,k,dim;
  uint iat;//,at_indx; //corey
  xvector<int> dims(3);
  xvector<double> cshift(3);  //corey
  vector<xvector<double> > cshifts; //corey
  //double zeroTol=1e-10;
  double radius;
  //cout << a.lattice << endl;
  b.lattice=supercell*a.lattice;
  //cout << b.lattice << endl;
  //b.lattice=a.lattice * trasp(supercell);  //Prof. YANG fixes bug on 05-12-2018
  //cout << b.lattice << endl;
  //exit(0);
  // the scale is kept the same....it is always saved as positive
  b.FixLattices();
  for(i=0;i<(int) b.num_each_type.size();i++) b.num_each_type[i]=0;
  for(i=0;i<(int) b.comp_each_type.size();i++) b.comp_each_type[i]=0;
  b.atoms.clear();
  sc2pcMap.clear(); pc2scMap.clear(); //CO 170722 - clear this everytime

  //calculate dist_nn_min
  if(b.dist_nn_min==AUROSTD_NAN){b.MinDist();}
  bool skew = SYM::isLatticeSkewed(b.lattice,b.dist_nn_min,b.sym_eps);

  double nx=supercell(1,1),ny=supercell(2,2),nz=supercell(3,3);
  //can only propagate if nx==ny==nz and diagonal
  //FALSE means we have a derivative structure
  //bool derivative_structure=!(abs(supercell(1,1)-supercell(2,2))<zeroTol && abs(supercell(1,1)-supercell(3,3))<zeroTol &&
  //                            abs(supercell(1,2))<zeroTol && abs(supercell(1,3))<zeroTol &&
  //                            abs(supercell(2,1))<zeroTol && abs(supercell(2,3))<zeroTol &&
  //                            abs(supercell(3,1))<zeroTol && abs(supercell(3,2))<zeroTol);
  bool derivative_structure=!(aurostd::isdiagonal(supercell) && 
                              aurostd::isequal(nx,ny) && 
                              aurostd::isequal(ny,nz));
  //COREY
  //VERY IMPORTANT, as we need to reduce the symmetry of the derivative structure
  //getting the basis serves as a validation for the symmetry operator
  //get_full_basis = get_full_basis || derivative_structure;

  //corey START
  // symmetry stuff
  b.ClearSymmetry();  //clear first
  //corey END

  radius=RadiusSphereLattice(b.lattice);
  dims=LatticeDimensionSphere(a.lattice,radius);dim=max(dims)+1;
  // if(LDEBUG) cerr << "DEBUG  dims=" << dims << " " << " radius=" << radius << endl;  // DEBUG

  bool match;
  bool pcmap;
  //corey START
  if(a.iatoms_calculated){
    //corey START
    //create bins
    for(uint ii=0;ii<a.iatoms.size();ii++){
      b.iatoms.push_back(vector<int>(0));
    }
    for(uint ia=0;ia<a.iatoms.size();ia++){
      for(uint iia=0;iia<a.iatoms[ia].size();iia++){
        pcmap=false;
        for(i=-dim;i<=dim;i++) {
          for(j=-dim;j<=dim;j++) {
            for(k=-dim;k<=dim;k++) {
              atom=a.atoms[a.iatoms[ia][iia]];
              //cerr << "atom " << a.iatoms[ia][iia] << " fpos_UNrot " << atom.fpos << endl;
              cshift=((double)i)*a.lattice(1)+((double)j)*a.lattice(2)+((double)k)*a.lattice(3);
              atom.cpos=atom.cpos+cshift;
              atom.fpos=b.c2f*atom.cpos; //C2F(b.lattice,atom.cpos);               // put in fractional of new basis
              //  atom.fpos=roundoff(atom.fpos);
              //cerr << "atom " << a.iatoms[ia][iia] << " fpos_rot   " << atom.fpos << endl;
              if(inCell(atom.fpos)){  //hard cut off
                //if(atom.fpos(1)>=-_eps_scell_ && atom.fpos(1)<1.0-_eps_scell_ &&
                //  atom.fpos(2)>=-_eps_scell_ && atom.fpos(2)<1.0-_eps_scell_ &&
                //  atom.fpos(3)>=-_eps_scell_ && atom.fpos(3)<1.0-_eps_scell_) {      // found something inside
                // atom=BringInCell(atom,b.lattice);
                b.num_each_type[atom.type]++;                  // CONVASP_MODE
                b.comp_each_type[atom.type]+=atom.partial_occupation_value;                  // CONVASP_MODE

                //we found a new iatom
                if(b.iatoms[ia].empty()){
                  atom.equivalent=b.atoms.size(); //reference self
                  atom.is_inequivalent=TRUE;  //iatom
                }else{
                  //eq atom
                  atom.equivalent=b.iatoms[ia][0];  //reference first iatom
                  atom.is_inequivalent=FALSE; //equivalent atom
                }

                //ijk
                atom.ijk(1)=i;atom.ijk(2)=j;atom.ijk(3)=k;

                b.atoms.push_back(atom);
                b.iatoms[ia].push_back(b.atoms.size()-1);
                //save cshifts for fgroups later...
                match=false;
                for(uint cf=0;cf<cshifts.size()&&!match;cf++){
                  if(aurostd::isequal(cshift,cshifts[cf],_ZERO_TOL_)){
                    match=true;
                  }
                }
                if(!match){
                  cshifts.push_back(cshift);
                }

                //mapping
                sc2pcMap.push_back(a.iatoms[ia][iia]);
                //pc2scMap is sort of irrelevant, we just need to pick ONE equivalent atom (there are many)
                //so just pick the first one
                if(!pcmap){
                  pc2scMap.push_back(b.atoms.size()-1);
                  pcmap=true;
                }
                //matching cpos does not work!
                //matching by index (i,j,k) does not work because of inCell()
                //if(!i&&!j&&!k) pc2scMap.push_back(b.atoms.size()-1);
                /*if(aurostd::identical(a.atoms[a.iatoms[ia][iia]].cpos,atom.cpos,1e-6)){
                  pc2scMap[a.iatoms[ia][iia]]=b.atoms.size()-1;
                  }*/
              }
            }
          }
        }
      }
    }

    // save the number of equivalents
    uint iequivalent=0;
    for(uint iat=0;iat<b.atoms.size();iat++) {
      if(b.atoms[iat].is_inequivalent) {
        b.atoms[iat].num_equivalents=b.iatoms[iequivalent].size();
        b.atoms[iat].index_iatoms=iequivalent;
        iequivalent++;
      }
    }

    b.iatoms_calculated=TRUE;
    //corey END
  }else{
    for(uint ia=0;ia<a.atoms.size();ia++) {
      pcmap=false;
      for(i=-dim;i<=dim;i++){
        for(j=-dim;j<=dim;j++){
          for(k=-dim;k<=dim;k++){
            atom=a.atoms[ia];
            //atom.cpos=atom.cpos+(((double)i)*a.lattice(1)+((double)j)*a.lattice(2)+((double)k)*a.lattice(3));
            cshift=((double)i)*a.lattice(1)+((double)j)*a.lattice(2)+((double)k)*a.lattice(3);
            atom.cpos=atom.cpos+cshift;
            atom.fpos=b.c2f*atom.cpos; //C2F(b.lattice,atom.cpos);               // put in fractional of new basis
            //  atom.fpos=roundoff(atom.fpos);
            if(inCell(atom.fpos)){
              //if(atom.fpos(1)>=-_eps_scell_ && atom.fpos(1)<1.0-_eps_scell_ &&
              //  atom.fpos(2)>=-_eps_scell_ && atom.fpos(2)<1.0-_eps_scell_ &&
              //  atom.fpos(3)>=-_eps_scell_ && atom.fpos(3)<1.0-_eps_scell_) {      // found something inside
              // atom=BringInCell(atom,b.lattice);
              b.num_each_type[atom.type]++;                  // CONVASP_MODE
              b.comp_each_type[atom.type]+=atom.partial_occupation_value;                  // CONVASP_MODE
              //ijk
              atom.ijk(1)=i;atom.ijk(2)=j;atom.ijk(3)=k;
              b.atoms.push_back(atom);
              //save cshifts for fgroups later...
              match=false;
              for(uint cf=0;cf<cshifts.size()&&!match;cf++){
                if(aurostd::isequal(cshift,cshifts[cf],_ZERO_TOL_)){
                  match=true;
                }
              }
              if(!match){
                cshifts.push_back(cshift);
              }

              //mapping
              sc2pcMap.push_back(ia);
              //pc2scMap is sort of irrelevant, we just need to pick ONE equivalent atom (there are many)
              //so just pick the first one
              if(!pcmap){
                pc2scMap.push_back(b.atoms.size()-1);
                pcmap=true;
              }
              //matching cpos does not work!
              //matching by index (i,j,k) does not work because of inCell()
              //if(!i&&!j&&!k) pc2scMap.push_back(b.atoms.size()-1);
              /*if(aurostd::identical(a.atoms[ia].cpos,atom.cpos,1e-6)){
                pc2scMap[a.iatoms[ia][iia]]=b.atoms.size()-1;
                }*/
            }
          }
        }
      }
    }
  }
  //corey END

  b.GetStoich();  //170724 CO
  b.MakeBasis(); // need to update NUMBER and BASIS

  // some check
  double fraction=((double) GetVol(b.lattice)/GetVol(a.lattice));
  double density_a=a.scale*a.scale*a.scale*abs(det(a.lattice))/a.atoms.size();
  double density_b=b.scale*b.scale*b.scale*abs(det(b.lattice))/b.atoms.size();
  if(abs(b.atoms.size()-a.atoms.size()*fraction)>0.1 || abs(density_a-density_b)>0.001) {
    cerr << "ERROR   xstructure xstructure::GetSuperCell     " << endl;
    cerr << "         supercell has the wrong number of atoms" << endl;
    cerr << "         b.atoms.size()     = " << b.atoms.size() << endl;
    cerr << "         a.atoms.size()     = " << a.atoms.size() << endl;
    cerr << "         b.scale            = " << b.scale << endl;
    cerr << "         a.scale            = " << a.scale << endl;
    cerr << "         fraction           = " << fraction << endl;
    cerr << "         supercell atoms    = " << fraction*a.atoms.size() << endl;
    cerr << "         b.density          = " << density_b << endl;
    cerr << "         a.density          = " << density_a << endl;
    cerr << "         b.lattice          = " << endl;
    cerr << b.lattice << endl;
    cerr << "         a.lattice          = " << endl;
    cerr << a.lattice << endl;
    exit(0);
  }

  bool pretend_uniform=false;//true;  //COREY TEST, REMOVE ME
  get_full_basis = !pretend_uniform && get_full_basis;

  //corey START
  //for now, we focus on pgroup, fgroup, iatoms, and site symmetry (agroup), add as you need
  if(get_symmetry){
    ofstream FileMESSAGE;
    _aflags aflags;
    _kflags kflags;
    bool _write_=false;    //COREY no verbose, annoys jose
    bool osswrite=false;   //COREY no verbose, annoys jose
    bool same_pgroups=true;
    bool calculated_pgroups=false;
    bool KRUN=true;
    //ostream& oss=cout;  //defined in macro at top of file
    if(derivative_structure){
      KRUN = KRUN && SYM::CalculatePointGroup(FileMESSAGE,b,aflags,_write_,osswrite,oss);
      if(LDEBUG){
        if(!KRUN){
          cerr << "Symmetry propagation FAILED with derivative structure at point group" << endl;
        }else{
          cerr << "Symmetry propagation PASSED with derivative structure at point group" << endl;
        }
      }
      same_pgroups=(KRUN && SYM::PointGroupsIdentical(a.pgroup,b.pgroup,b.sym_eps,false));  //DX 12/7/17 - added is_same_lattice
      calculated_pgroups=KRUN;
    }
    if(KRUN && !pretend_uniform && !same_pgroups){
      //we calculate earlier to see if there's a mismatch
      //SYM::CalculatePointGroup(FileMESSAGE,b,aflags,_write_,osswrite,oss);
      //do all at same sym_eps as primitive cell
      KRUN = KRUN && SYM::CalculateFactorGroup(FileMESSAGE,b,aflags,_write_,osswrite,oss);
      if(LDEBUG){
        if(!KRUN){
          cerr << "Symmetry propagation FAILED with derivative structure at factor group" << endl;
        }else{
          cerr << "Symmetry propagation PASSED with derivative structure at factor group" << endl;
        }
      }
      KRUN = KRUN && SYM::CalculatePointGroupCrystal(FileMESSAGE,b,aflags,_write_,osswrite,oss);
      if(LDEBUG){
        if(!KRUN){
          cerr << "Symmetry propagation FAILED with derivative structure at point group crystal" << endl;
        }else{
          cerr << "Symmetry propagation PASSED with derivative structure at point group crystal" << endl;
        }
      }
      //if(!a.iatoms_calculated){
      KRUN = KRUN && SYM::CalculateInequivalentAtoms(FileMESSAGE,b,aflags,_write_,osswrite,oss); //100% necessary, new pgroups means different symmetry, different iatoms
      if(LDEBUG){
        if(!KRUN){
          cerr << "Symmetry propagation FAILED with derivative structure at iatoms" << endl;
        }else{
          cerr << "Symmetry propagation PASSED with derivative structure at iatoms" << endl;
        }
      }
      //}
      //AGAIN, many fgroups, but not many pgroups for derivative structures, let's see if this is faster...
      //if(!SYM::CalculateSitePointGroup(FileMESSAGE,b,true,aflags,_write_,osswrite,oss)){  //iatoms only  
      KRUN = KRUN && SYM::CalculateSitePointGroup(FileMESSAGE,b,(get_full_basis ? 0 : 1 ),aflags,_write_,osswrite,oss);  //don't waste time calculating basis_map for eatoms, really never use them anyway
      if(LDEBUG){
        if(!KRUN){
          cerr << "Symmetry propagation FAILED with derivative structure at agroup" << endl;
        }else{
          cerr << "Symmetry propagation PASSED with derivative structure at agroup" << endl;
        }
      }
      //}
      //validate that we have good symmetry here
    }else if(KRUN){
      //////////////////////////////////////////////////////////////////////////
      //PGROUP
      if(KRUN && a.pgroup_calculated && !calculated_pgroups){
        for(uint i=0;i<a.pgroup.size() && KRUN;i++){
          pSymOp=a.pgroup[i];
          pSymOp.Uf=b.c2f * pSymOp.Uc * b.f2c;
          pSymOp.basis_atoms_map.clear();
          pSymOp.basis_types_map.clear();
          pSymOp.basis_map_calculated=FALSE;
          //CO and DX, getFullSymBasis does not make sense for just rotations
          //if(get_full_basis){
          //  KRUN = KRUN && SYM::getFullSymBasis(b.atoms,b.lattice,b.c2f,b.f2c,pSymOp,FALSE,skew,b.sym_eps,pSymOp.basis_atoms_map,pSymOp.basis_types_map);
          //  if(LDEBUG){
          //    if(!KRUN){
          //      cerr << "Symmetry propagation FAILED with uniform supercell structure at point group" << endl;
          //    }else{
          //      cerr << "Symmetry propagation PASSED with uniform supercell structure at point group" << endl;
            //}
          //  }
          //  pSymOp.basis_map_calculated=KRUN;
          //}
          //if(KRUN){b.pgroup.push_back(pSymOp);}
          if(KRUN){SYM::AddSymmetryToStructure(b,pSymOp.Uc,pSymOp.Uf,pSymOp.ctau,pSymOp.ftau,pSymOp.ctrasl,pSymOp.ftrasl,pSymOp.basis_atoms_map,pSymOp.basis_types_map,pSymOp.basis_map_calculated,_PGROUP_,FALSE);}  //CO 170706 - make sure quaternion is updated
        }
        //if(KRUN && b.pgroup.size()){
        //b.pgroup_calculated=true;
        b.pgroup_calculated=(KRUN && b.pgroup.size());
        //}
      }
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      //FGROUP
      if(KRUN && a.fgroup_calculated) { 
        for(uint fg=0;fg<a.fgroup.size() && KRUN;fg++){
          for(uint cs=0;cs<cshifts.size() && KRUN;cs++){
            fSymOp=a.fgroup[fg];
            fSymOp.ctau=fSymOp.ctau+cshifts[cs];
            fSymOp.ftau=b.c2f*fSymOp.ctau;
            if(inCell(fSymOp.ftau)){ //DX CHANGE HERE; NO MORE TOL_ABC_RES
              // We have to correct the Uf for each symop since we have changed the lattice...
              fSymOp.Uf=b.c2f * fSymOp.Uc * b.f2c;
              fSymOp.basis_atoms_map.clear();
              fSymOp.basis_types_map.clear();
              fSymOp.basis_map_calculated=FALSE;
              for(uint iii=0;iii<b.atoms.size();iii++){
                fSymOp.basis_atoms_map.push_back(0);
                fSymOp.basis_types_map.push_back(0);
              }
              fSymOp.basis_map_calculated=false;
              //calculate basis_atoms_map and basis_types_map
              if(get_full_basis) {
                //NOPE, we will calculate
                KRUN = KRUN && SYM::getFullSymBasis(b.atoms,b.lattice,b.c2f,b.f2c,fSymOp,TRUE,skew,b.sym_eps,fSymOp.basis_atoms_map,fSymOp.basis_types_map);
                if(LDEBUG){
                  if(!KRUN){
                    cerr << "Symmetry propagation FAILED with uniform supercell structure at factor group" << endl;
                  }else{
                    cerr << "Symmetry propagation PASSED with uniform supercell structure at factor group" << endl;
                  }
                }
                //if(!SYM::getFullSymBasis(b.atoms,b.lattice,b.c2f,b.f2c,fSymOp,TRUE,skew,b.sym_eps,fSymOp.basis_atoms_map,fSymOp.basis_types_map)) {
                //cerr << "Unable to find atom/types basis for fgroup" << endl;
                //cerr << fSymOp << endl;
                //exit(0);
                //KRUN = FALSE;
                //}
                fSymOp.basis_map_calculated=KRUN;
              }
              //if(KRUN){b.fgroup.push_back(fSymOp);}
              if(KRUN){SYM::AddSymmetryToStructure(b,fSymOp.Uc,fSymOp.Uf,fSymOp.ctau,fSymOp.ftau,fSymOp.ctrasl,fSymOp.ftrasl,fSymOp.basis_atoms_map,fSymOp.basis_types_map,fSymOp.basis_map_calculated,_FGROUP_,FALSE);}  //CO 170706 - make sure quaternion is updated
            }
          }
        }
        //if(b.fgroup.size()){
        //  b.fgroup_calculated=TRUE;   //there are more, but these are the important ones
        b.fgroup_calculated=(KRUN && b.fgroup.size());   //there are more, but these are the important ones
        //}
      }
      //////////////////////////////////////////////////////////////////////////
      
      //////////////////////////////////////////////////////////////////////////
      //PGROUP_XTAL
      if(KRUN && a.pgroup_xtal_calculated){
        for(uint i=0;i<a.pgroup_xtal.size() && KRUN;i++){
          pSymOp=a.pgroup_xtal[i];
          pSymOp.Uf=b.c2f * pSymOp.Uc * b.f2c;
          pSymOp.basis_atoms_map.clear();
          pSymOp.basis_types_map.clear();
          pSymOp.basis_map_calculated=FALSE;
          //CO and DX, getFullSymBasis does not make sense for just rotations
          //if(get_full_basis){
          //  KRUN = KRUN && SYM::getFullSymBasis(b.atoms,b.lattice,b.c2f,b.f2c,pSymOp,FALSE,skew,b.sym_eps,pSymOp.basis_atoms_map,pSymOp.basis_types_map);
          //  if(LDEBUG){
          //    if(!KRUN){
          //      cerr << "Symmetry propagation FAILED with uniform supercell structure at point group crystal" << endl;
          //    }else{
          //      cerr << "Symmetry propagation PASSED with uniform supercell structure at point group crystal" << endl;
          //    }
          //  }
          //  pSymOp.basis_map_calculated=KRUN;
          //}
          //if(KRUN){b.pgroup_xtal.push_back(pSymOp);}
          if(KRUN){SYM::AddSymmetryToStructure(b,pSymOp.Uc,pSymOp.Uf,pSymOp.ctau,pSymOp.ftau,pSymOp.ctrasl,pSymOp.ftrasl,pSymOp.basis_atoms_map,pSymOp.basis_types_map,pSymOp.basis_map_calculated,_PGROUP_XTAL_,FALSE);}  //CO 170706 - make sure quaternion is updated
        }
        //if(KRUN && b.pgroup_xtal.size()){
        //b.pgroup_xtal_calculated=true;
        b.pgroup_xtal_calculated=(KRUN && b.pgroup_xtal.size());
        //}
      }
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      //IF A.IATOMS, then we already propagated, otherwise, just run full routine for iatoms
      if(KRUN && !a.iatoms_calculated){
        KRUN = KRUN && SYM::CalculateInequivalentAtoms(FileMESSAGE,b,aflags,_write_,osswrite,oss);  //do all atoms
        if(LDEBUG){
          if(!KRUN){
            cerr << "Symmetry propagation FAILED with uniform supercell structure at agroup" << endl;
          }else{
            cerr << "Symmetry propagation PASSED with uniform supercell structure at agroup" << endl;
          }
        }
      }
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      //AGROUP
      if(KRUN && a.agroup_calculated){
        if(b.iatoms_calculated && b.fgroup_calculated && get_full_basis){ //the only way we get a speed up is if we can use fgroup.basis_atoms_map
          //xstructure bb=b;
          deque<_atom> b_atoms=b.atoms;
          xvector<double> origin(3),frigin(3);
          //create space for agroups
          for(uint iia=0;iia<b.atoms.size();iia++){
            b.agroup.push_back(vector<_sym_op>(0));
          }
          for(uint ia=0;ia<b.iatoms.size() && KRUN;ia++){
            iat=b.iatoms[ia][0];
            //let's recycle what we have to improve speed
            origin=b.atoms[iat].cpos;
            frigin=b.atoms[iat].fpos;
            for(uint ii=0;ii<b.atoms.size();ii++){
              //go back to original first, then subtract new origin
              b_atoms[ii].cpos=b.atoms[ii].cpos-origin;
              b_atoms[ii].fpos=b.atoms[ii].fpos-frigin;
              //now bring in cell
              b_atoms[ii]=BringInCell(b_atoms[ii],b.lattice);
            }
            //bb.ShifOriginToAtom(iat);
            //bb.BringInCell();

            //IATOMS ONLY
            for(uint iia=0;iia<a.agroup[sc2pcMap[iat]].size() && KRUN;iia++){
              aSymOp=a.agroup[sc2pcMap[iat]][iia];
              // We have to correct the Uf for each symop since we have changed the lattice...
              aSymOp.Uf=b.c2f * aSymOp.Uc * b.f2c;
              //no longer necessary since we force a basis map calculation
              /*aSymOp.basis_atoms_map.clear();
                aSymOp.basis_types_map.clear();
                aSymOp.basis_map_calculated=FALSE;
                for(uint iii=0;iii<b.atoms.size();iii++){
                aSymOp.basis_atoms_map.push_back(0);
                aSymOp.basis_types_map.push_back(0);
                }
                fSymOp.basis_map_calculated=false;*/
              KRUN = KRUN && SYM::getFullSymBasis(b_atoms,b.lattice,b.c2f,b.f2c,aSymOp,TRUE,skew,b.sym_eps,aSymOp.basis_atoms_map,aSymOp.basis_types_map);
              if(LDEBUG){
                if(!KRUN){
                  cerr << "Symmetry propagation FAILED with uniform supercell structure at agroup" << endl;
                }else{
                  cerr << "Symmetry propagation PASSED with uniform supercell structure at agroup" << endl;
                }
              }
              //if(!SYM::getFullSymBasis(b_atoms,b.lattice,b.c2f,b.f2c,aSymOp,TRUE,skew,b.sym_eps,aSymOp.basis_atoms_map,aSymOp.basis_types_map)) {
              //cerr << "Unable to find atom/types basis for agroup" << endl;
              //cerr << aSymOp << endl;
              //exit(0);
              //KRUN = FALSE;
              //}
              aSymOp.basis_map_calculated=KRUN;
              //if(KRUN){b.agroup[iat].push_back(aSymOp);}
              if(KRUN){SYM::AddSymmetryToStructure(b,iat,aSymOp.Uc,aSymOp.Uf,aSymOp.ctau,aSymOp.ftau,aSymOp.ctrasl,aSymOp.ftrasl,aSymOp.basis_atoms_map,aSymOp.basis_types_map,aSymOp.basis_map_calculated,_AGROUP_,FALSE);}  //CO 170706 - make sure quaternion is updated
            }
          }

          //EATOMS FOLLOW
          KRUN = KRUN && SYM::CalculateSitePointGroup_EquivalentSites(b,get_full_basis,b.sym_eps);
          if(LDEBUG){
            if(!KRUN){
              cerr << "Symmetry propagation FAILED with uniform supercell structure at agroup equivalent" << endl;
            }else{
              cerr << "Symmetry propagation PASSED with uniform supercell structure at agroup equivalent" << endl;
            }
          }
          //KRUN = FALSE;
          //cerr << "Unable to propagate site symmetry to equivalent atoms." << endl;
          //cerr << "Unable to propagate site symmetry to equivalent atoms. Exiting." << endl;
          //exit(0);
          //}
        }else{
          //can be faster than procedure above because there are MANY fgroups
          KRUN = KRUN && SYM::CalculateSitePointGroup(FileMESSAGE,b,1,aflags,_write_,osswrite,oss);  //we already know get_full_basis==FALSE, so don't waste time calculating for eatoms
          if(LDEBUG){
            if(!KRUN){
              cerr << "Symmetry propagation FAILED with uniform supercell structure at agroup" << endl;
            }else{
              cerr << "Symmetry propagation PASSED with uniform supercell structure at agroup" << endl;
            }
          }
        }
        //if(b.agroup.size() && b.agroup[0].size()){      //just a fast check to see we have agroups somewhere (we should always get identity)
        b.agroup_calculated=(KRUN && b.agroup[0].size()); //just a fast check to see we have agroups somewhere (we should always get identity)
        //  b.agroup_calculated=TRUE;
        //}
      }
      //////////////////////////////////////////////////////////////////////////
    }
    if(!KRUN){
      oss << (aflags.QUIET?"":"00000  MESSAGE ") << "SUPERCELL Symmetry propagation FAILED" << Message(aflags,"user,host,time") << endl;
      oss << (aflags.QUIET?"":"00000  MESSAGE ") << "SUPERCELL Symmetry retrying with symmetry scan" << Message(aflags,"user,host,time") << endl;
      pflow::PerformFullSymmetry(b,FileMESSAGE,aflags,kflags,osswrite,oss);
      //FOOLPROOF!!!!!!!!!
      //no need for krun here with force_perform: if it fails with full scan, it will calculate at default tolerance and keep going!
    }
  }
  
  //cerr << "prim: " << a.pgroup.size() << ", sup: " << b.pgroup.size() << endl;    //COREY REMOVE

  //corey END

  b.ReScale(aa.scale);  //the nuclear option, the only way not to mess around with scale EVERYWHERE
  return b;
}
// COREY - END

xstructure GetSuperCell(const xstructure& a, const xvector<double>& supercell,vector<int>& sc2pcMap,vector<int>& pc2scMap,
                        bool get_symmetry, bool get_full_basis) {
  //xstructure GetSuperCell(const xstructure& a, const xvector<double>& supercell) {
  xmatrix<double> _supercell(3,3);
  if(supercell.rows==9) {
    _supercell(1,1)=supercell(1);_supercell(1,2)=supercell(2);_supercell(1,3)=supercell(3);
    _supercell(2,1)=supercell(4);_supercell(2,2)=supercell(5);_supercell(2,3)=supercell(6);
    _supercell(3,1)=supercell(7);_supercell(3,2)=supercell(8);_supercell(3,3)=supercell(9);
    return GetSuperCell(a,_supercell,sc2pcMap,pc2scMap,get_symmetry,get_full_basis);
  }
  if(supercell.rows==3) {
    _supercell(1,1)=supercell(1);_supercell(2,2)=supercell(2);_supercell(3,3)=supercell(3);
    return GetSuperCell(a,_supercell,sc2pcMap,pc2scMap,get_symmetry,get_full_basis);
  }
  cerr << "GetSuperCell - vector must have 9 or 3 elements" << endl;
  exit(0);
}

xstructure GetSuperCell(const xstructure& a, const xvector<int>& supercell,vector<int>& sc2pcMap,vector<int>& pc2scMap,
                        bool get_symmetry, bool get_full_basis) {
  //xstructure GetSuperCell(const xstructure& a, const xvector<int>& supercell) {
  xmatrix<double> _supercell(3,3);
  if(supercell.rows==9) {
    _supercell(1,1)=supercell(1);_supercell(1,2)=supercell(2);_supercell(1,3)=supercell(3);
    _supercell(2,1)=supercell(4);_supercell(2,2)=supercell(5);_supercell(2,3)=supercell(6);
    _supercell(3,1)=supercell(7);_supercell(3,2)=supercell(8);_supercell(3,3)=supercell(9);
    return GetSuperCell(a,_supercell,sc2pcMap,pc2scMap,get_symmetry,get_full_basis);
  }
  if(supercell.rows==3) {
    _supercell(1,1)=supercell(1);_supercell(2,2)=supercell(2);_supercell(3,3)=supercell(3);
    return GetSuperCell(a,_supercell,sc2pcMap,pc2scMap,get_symmetry,get_full_basis);
  }
  cerr << "GetSuperCell - vector must have 9 or 3 elements" << endl;
  exit(0);
}

xstructure GetSuperCell(const xstructure& a, const int& sc11,const int& sc12,const int& sc13, const int& sc21,const int& sc22,const int& sc23, const int& sc31,const int& sc32,const int& sc33,vector<int>& sc2pcMap,vector<int>& pc2scMap,
                        bool get_symmetry, bool get_full_basis) {
  //xstructure GetSuperCell(const xstructure& a, const int& sc11,const int& sc12,const int& sc13, const int& sc21,const int& sc22,const int& sc23, const int& sc31,const int& sc32,const int& sc33) {
  xmatrix<double> _supercell(3,3);
  _supercell.clear();
  _supercell(1,1)=(double) sc11;_supercell(1,2)=(double) sc12;_supercell(1,3)=(double) sc13;
  _supercell(2,1)=(double) sc21;_supercell(2,2)=(double) sc22;_supercell(2,3)=(double) sc23;
  _supercell(3,1)=(double) sc31;_supercell(3,2)=(double) sc32;_supercell(3,3)=(double) sc33;
  return GetSuperCell(a,_supercell,sc2pcMap,pc2scMap,get_symmetry,get_full_basis);
}

xstructure GetSuperCell(const xstructure& a, const int& sc1,const int& sc2,const int& sc3,vector<int>& sc2pcMap,vector<int>& pc2scMap,
                        bool get_symmetry, bool get_full_basis) {
  //xstructure GetSuperCell(const xstructure& a, const int& sc1,const int& sc2,const int& sc3) {
  xmatrix<double> _supercell(3,3);
  _supercell.clear();
  _supercell(1,1)=(double) sc1;
  _supercell(2,2)=(double) sc2;
  _supercell(3,3)=(double) sc3;
  return GetSuperCell(a,_supercell,sc2pcMap,pc2scMap,get_symmetry,get_full_basis);
}

//corey START
xstructure GetSuperCell(const xstructure& aa, const xmatrix<double> &supercell) {
  vector<int> sc2pcMap;
  vector<int> pc2scMap;
  bool get_symmetry=false;
  bool get_full_basis=false;
  return GetSuperCell(aa,supercell,sc2pcMap,pc2scMap,get_symmetry,get_full_basis);
}

xstructure GetSuperCell(const xstructure& a, const xvector<double>& supercell) {
  vector<int> sc2pcMap;
  vector<int> pc2scMap;
  bool get_symmetry=false;
  bool get_full_basis=false;
  return GetSuperCell(a,supercell,sc2pcMap,pc2scMap,get_symmetry,get_full_basis);
}

xstructure GetSuperCell(const xstructure& a, const xvector<int>& supercell) {
  vector<int> sc2pcMap;
  vector<int> pc2scMap;
  bool get_symmetry=false;
  bool get_full_basis=false;
  return GetSuperCell(a,supercell,sc2pcMap,pc2scMap,get_symmetry,get_full_basis);
}

xstructure GetSuperCell(const xstructure& a, const int& sc11,const int& sc12,const int& sc13, const int& sc21,const int& sc22,const int& sc23, const int& sc31,const int& sc32,const int& sc33) {
  vector<int> sc2pcMap;
  vector<int> pc2scMap;
  bool get_symmetry=false;
  bool get_full_basis=false;
  return GetSuperCell(a,sc11,sc12,sc13,sc21,sc22,sc23,sc31,sc32,sc33,sc2pcMap,pc2scMap,get_symmetry,get_full_basis);
}

xstructure GetSuperCell(const xstructure& a, const int& sc1,const int& sc2,const int& sc3) {
  vector<int> sc2pcMap;
  vector<int> pc2scMap;
  bool get_symmetry=false;
  bool get_full_basis=false;
  return GetSuperCell(a,sc1,sc2,sc3,sc2pcMap,pc2scMap,get_symmetry,get_full_basis);
}
//corey END


// ***************************************************************************
// Function ClearSymmetry
// ***************************************************************************
void xstructure::ClearSymmetry(void) {
  // PGROUP ----------------------------
  pgroup.clear();            // just initialize
  pgroup_calculated=FALSE;
  // PGROUP_XTAL ----------------------------
  pgroup_xtal.clear();        // just initialize
  pgroup_xtal_calculated=FALSE;
  crystal_family="";crystal_system="";point_group_crystal_class="";
  point_group_Shoenflies="";point_group_Hermann_Mauguin="";point_group_orbifold="";
  point_group_type="";point_group_order="";point_group_structure="";
  // PGROUPK ----------------------------
  pgroupk.clear();            // just initialize
  pgroupk_calculated=FALSE;
  // PGROUPK_XTAL ----------------------------
  pgroupk_xtal.clear();            // just initialize //DX 12/5/17 - Added pgroupk_xtal
  pgroupk_xtal_calculated=FALSE;                      //DX 12/5/17 - Added pgroupk_xtal
  // FGROUP ----------------------------
  fgroup.clear();            // just initialize
  fgroup_calculated=FALSE;
  // SGROUP ----------------------------
  sgroup_radius=-_calculate_symmetry_default_sgroup_radius_; // symmetry not calculated
  sgroup_radius_dims.clear();
  sgroup.clear();            // just initialize
  sgroup_calculated=FALSE;
  // SITE POINT GROUP ------------------
  agroup_calculated=FALSE;
  for(uint i=0;i<agroup.size();i++)
    agroup.at(i).clear();
  agroup.clear();
  // INEQUIVALENT ATOMS ----------------
  iatoms_calculated=FALSE;
  for(uint i=0;i<iatoms.size();i++)
    iatoms.at(i).clear();
  iatoms.clear();
}

// DX - Consider using pflow::CalculateFullSymmetry in aflow_aconvasp_main.cpp.
//      It contains consistency checks for the symmetry analysis.
// ***************************************************************************
// Function CalculateSymmetry
// ***************************************************************************
bool xstructure::CalculateSymmetry(bool ossverbose,double radius) {
  ofstream FileDevNull("/dev/null");
  
  _aflags aflags; aflags.Directory="./"; aflags.QUIET=TRUE;
  _kflags kflags; pflow::defaultKFlags4SymWrite(kflags,false); pflow::defaultKFlags4SymCalc(kflags,true);
  
  (*this).LatticeReduction_avoid=TRUE;
  (*this).sgroup_radius=radius;
  
  return pflow::PerformFullSymmetry(*this,FileDevNull,aflags,kflags,ossverbose,oss);
  
  //SYM::CalculatePointGroup(FileDevNull,*this,aflags,FALSE,ossverbose,cout);
  //SYM::CalculatePointGroupKlattice(FileDevNull,*this,aflags,FALSE,ossverbose,cout);
  //SYM::CalculateSitePointGroup(FileDevNull,*this,aflags,FALSE,ossverbose,cout);
  //SYM::CalculateFactorGroup(FileDevNull,*this,aflags,FALSE,ossverbose,cout);
  //SYM::CalculateInequivalentAtoms(FileDevNull,*this,aflags,FALSE,ossverbose,cout);
  //SYM::CalculateSpaceGroup(FileDevNull,*this,aflags,FALSE,ossverbose,cout);
}

bool xstructure::CalculateSymmetry(void) {
  LatticeReduction_avoid=TRUE;
  return CalculateSymmetry(FALSE,_calculate_symmetry_default_sgroup_radius_*MaxStructureLattice(*this));
}

bool CalculateSymmetry(xstructure& str,bool ossverbose,ostream& oss,bool fffverbose,double radius) {
  ofstream FileDevNull("/dev/null");
  
  _aflags aflags; aflags.Directory="./"; aflags.QUIET=TRUE;
  _kflags kflags; pflow::defaultKFlags4SymWrite(kflags,fffverbose); pflow::defaultKFlags4SymCalc(kflags,true);
  
  str.LatticeReduction_avoid=TRUE;
  str.sgroup_radius=radius;
  
  return pflow::PerformFullSymmetry(str,FileDevNull,aflags,kflags,ossverbose,oss);
  //SYM::CalculatePointGroup(FileDevNull,str,aflags,fffverbose,ossverbose,oss);
  //SYM::CalculatePointGroupKlattice(FileDevNull,str,aflags,fffverbose,ossverbose,oss);
  //SYM::CalculateSitePointGroup(FileDevNull,str,aflags,fffverbose,ossverbose,oss);
  //SYM::CalculateFactorGroup(FileDevNull,str,aflags,fffverbose,ossverbose,oss);
  //SYM::CalculateInequivalentAtoms(FileDevNull,str,aflags,fffverbose,ossverbose,oss);
  //SYM::CalculateSpaceGroup(FileDevNull,str,aflags,TRUE,ossverbose,oss);
}

bool CalculateSymmetry(xstructure& str,bool ossverbose,ostream& oss,bool fffverbose) {
  str.LatticeReduction_avoid=TRUE;
  return CalculateSymmetry(str,ossverbose,oss,fffverbose,_calculate_symmetry_default_sgroup_radius_*MaxStructureLattice(str));
}

bool CalculateSymmetry(xstructure& str,bool ossverbose,ostream& oss,double radius) {
  str.LatticeReduction_avoid=TRUE;
  _aflags aflags;
  aflags.Directory="./";aflags.QUIET=TRUE;
  return CalculateSymmetry(str,ossverbose,oss,FALSE,radius);
}

bool CalculateSymmetry(xstructure& str,bool ossverbose,double radius) {
  str.LatticeReduction_avoid=TRUE;
  _aflags aflags;
  aflags.Directory="./";aflags.QUIET=TRUE;
  return CalculateSymmetry(str,ossverbose,cout,FALSE,radius);
}

bool CalculateSymmetry(xstructure& str,double radius) {
  str.LatticeReduction_avoid=TRUE;
  _aflags aflags;
  aflags.Directory="./";aflags.QUIET=TRUE;
  return CalculateSymmetry(str,FALSE,cout,FALSE,radius);
}

bool CalculateSymmetry(xstructure& str,bool ossverbose) {
  str.LatticeReduction_avoid=TRUE;
  _aflags aflags;
  aflags.Directory="./";aflags.QUIET=TRUE;
  return CalculateSymmetry(str,ossverbose,cout,FALSE,_calculate_symmetry_default_sgroup_radius_*MaxStructureLattice(str));
}

bool CalculateSymmetry(xstructure& str) {
  str.LatticeReduction_avoid=TRUE;
  _aflags aflags;
  aflags.Directory="./";aflags.QUIET=TRUE;
  return CalculateSymmetry(str,FALSE,cout,FALSE,_calculate_symmetry_default_sgroup_radius_*MaxStructureLattice(str));
}

// DX - Consider using pflow::CalculateFullSymmetry in aflow_aconvasp_main.cpp.
//      It contains consistency checks for the symmetry analysis.
// ***************************************************************************
// Function CalculateSymmetryPointGroup
// ***************************************************************************
void xstructure::CalculateSymmetryPointGroup(bool ossverbose) {
  ofstream FileDevNull("/dev/null");
  _aflags aflags;
  aflags.Directory="./";aflags.QUIET=TRUE;
  LatticeReduction_avoid=TRUE;
  SYM::CalculatePointGroup(FileDevNull,*this,aflags,FALSE,ossverbose,cout);
}

void xstructure::CalculateSymmetryPointGroup(void) {
  LatticeReduction_avoid=TRUE;
  CalculateSymmetryPointGroup(FALSE);
}

void CalculateSymmetryPointGroup(xstructure& str,bool ossverbose,ostream& oss,bool fffverbose) {
  ofstream FileDevNull("/dev/null");
  _aflags aflags;
  aflags.Directory="./";aflags.QUIET=TRUE;
  str.LatticeReduction_avoid=TRUE;
  SYM::CalculatePointGroup(FileDevNull,str,aflags,fffverbose,ossverbose,oss);
}

void CalculateSymmetryPointGroup(xstructure& str,bool ossverbose,ostream& oss) {
  str.LatticeReduction_avoid=TRUE;
  _aflags aflags;
  aflags.Directory="./";aflags.QUIET=TRUE;
  CalculateSymmetryPointGroup(str,ossverbose,oss,FALSE);
}

void CalculateSymmetryPointGroup(xstructure& str,bool ossverbose) {
  str.LatticeReduction_avoid=TRUE;
  _aflags aflags;
  aflags.Directory="./";aflags.QUIET=TRUE;
  CalculateSymmetryPointGroup(str,ossverbose,cout,FALSE);
}

void CalculateSymmetryPointGroup(xstructure& str) {
  str.LatticeReduction_avoid=TRUE;
  _aflags aflags;
  aflags.Directory="./";aflags.QUIET=TRUE;
  CalculateSymmetryPointGroup(str,FALSE,cout,FALSE);
}

// DX - Consider using pflow::CalculateFullSymmetry in aflow_aconvasp_main.cpp.
//      It contains consistency checks for the symmetry analysis.
// ***************************************************************************
// Function CalculateSymmetryPointGroupCrystal
// ***************************************************************************
void xstructure::CalculateSymmetryPointGroupCrystal(bool ossverbose) {
  ofstream FileDevNull("/dev/null");
  _aflags aflags;
  aflags.Directory="./";aflags.QUIET=TRUE;
  LatticeReduction_avoid=TRUE;
  if(pgroup_calculated==FALSE) SYM::CalculatePointGroup(FileDevNull,*this,aflags,FALSE,ossverbose,cout);
  if(fgroup_calculated==FALSE) SYM::CalculateFactorGroup(FileDevNull,*this,aflags,FALSE,ossverbose,cout);
  SYM::CalculatePointGroupCrystal(FileDevNull,*this,aflags,FALSE,ossverbose,cout);
}

void xstructure::CalculateSymmetryPointGroupCrystal(void) {
  LatticeReduction_avoid=TRUE;
  if(pgroup_calculated==FALSE) CalculateSymmetryPointGroup(FALSE);
  if(fgroup_calculated==FALSE) CalculateSymmetryFactorGroup(FALSE);
  CalculateSymmetryPointGroupCrystal(FALSE);
}

void CalculateSymmetryPointGroupCrystal(xstructure& str,bool ossverbose,ostream& oss,bool fffverbose) {
  ofstream FileDevNull("/dev/null");
  _aflags aflags;
  aflags.Directory="./";aflags.QUIET=TRUE;
  str.LatticeReduction_avoid=TRUE;
  if(str.pgroup_calculated==FALSE) SYM::CalculatePointGroup(FileDevNull,str,aflags,fffverbose,ossverbose,oss);
  if(str.fgroup_calculated==FALSE) SYM::CalculateFactorGroup(FileDevNull,str,aflags,fffverbose,ossverbose,oss);
  SYM::CalculatePointGroupCrystal(FileDevNull,str,aflags,fffverbose,ossverbose,oss);
}

void CalculateSymmetryPointGroupCrystal(xstructure& str,bool ossverbose,ostream& oss) {
  str.LatticeReduction_avoid=TRUE;
  _aflags aflags;
  aflags.Directory="./";aflags.QUIET=TRUE;
  CalculateSymmetryPointGroupCrystal(str,ossverbose,oss,FALSE);
}

void CalculateSymmetryPointGroupCrystal(xstructure& str,bool ossverbose) {
  str.LatticeReduction_avoid=TRUE;
  _aflags aflags;
  aflags.Directory="./";aflags.QUIET=TRUE;
  CalculateSymmetryPointGroupCrystal(str,ossverbose,cout,FALSE);
}

void CalculateSymmetryPointGroupCrystal(xstructure& str) {
  str.LatticeReduction_avoid=TRUE;
  _aflags aflags;
  aflags.Directory="./";aflags.QUIET=TRUE;
  if(str.pgroup_calculated==FALSE) CalculateSymmetryPointGroup(str,FALSE,cout,FALSE);
  if(str.fgroup_calculated==FALSE) CalculateSymmetryFactorGroup(str,FALSE,cout,FALSE);
  CalculateSymmetryPointGroupCrystal(str,FALSE,cout,FALSE);
}

// DX - Consider using pflow::CalculateFullSymmetry in aflow_aconvasp_main.cpp.
//      It contains consistency checks for the symmetry analysis.
// ***************************************************************************
// Function CalculateSymmetryFactorGroup
// ***************************************************************************
void xstructure::CalculateSymmetryFactorGroup(bool ossverbose) {
  ofstream FileDevNull("/dev/null");
  _aflags aflags;
  aflags.Directory="./";aflags.QUIET=TRUE;
  LatticeReduction_avoid=TRUE;
  SYM::CalculatePointGroup(FileDevNull,*this,aflags,FALSE,ossverbose,cout);
  SYM::CalculateFactorGroup(FileDevNull,*this,aflags,FALSE,ossverbose,cout);
}

void xstructure::CalculateSymmetryFactorGroup(void) {
  LatticeReduction_avoid=TRUE;
  CalculateSymmetryFactorGroup(FALSE);
}

void CalculateSymmetryFactorGroup(xstructure& str,bool ossverbose,ostream& oss,bool fffverbose) {
  ofstream FileDevNull("/dev/null");
  _aflags aflags;
  aflags.Directory="./";aflags.QUIET=TRUE;
  str.LatticeReduction_avoid=TRUE;
  SYM::CalculatePointGroup(FileDevNull,str,aflags,fffverbose,ossverbose,oss);
  SYM::CalculateFactorGroup(FileDevNull,str,aflags,fffverbose,ossverbose,oss);
}

void CalculateSymmetryFactorGroup(xstructure& str,bool ossverbose,ostream& oss) {
  str.LatticeReduction_avoid=TRUE;
  _aflags aflags;
  aflags.Directory="./";aflags.QUIET=TRUE;
  CalculateSymmetryFactorGroup(str,ossverbose,oss,FALSE);
}

void CalculateSymmetryFactorGroup(xstructure& str,bool ossverbose) {
  str.LatticeReduction_avoid=TRUE;
  _aflags aflags;
  aflags.Directory="./";aflags.QUIET=TRUE;
  CalculateSymmetryFactorGroup(str,ossverbose,cout,FALSE);
}

void CalculateSymmetryFactorGroup(xstructure& str) {
  str.LatticeReduction_avoid=TRUE;
  _aflags aflags;
  aflags.Directory="./";aflags.QUIET=TRUE;
  CalculateSymmetryFactorGroup(str,FALSE,cout,FALSE);
}

// DX - Consider using pflow::CalculateFullSymmetry in aflow_aconvasp_main.cpp.
//      It contains consistency checks for the symmetry analysis.
// ***************************************************************************
// Function CalculateSymmetryPointGroupKlattice
// ***************************************************************************
void xstructure::CalculateSymmetryPointGroupKlattice(bool ossverbose) {
  ofstream FileDevNull("/dev/null");
  _aflags aflags;
  aflags.Directory="./";aflags.QUIET=TRUE;
  // LatticeReduction_avoid=TRUE;  // not necssary for klattice
  SYM::CalculatePointGroupKlattice(FileDevNull,*this,aflags,FALSE,ossverbose,cout);
}

void xstructure::CalculateSymmetryPointGroupKlattice(void) {
  // LatticeReduction_avoid=TRUE;  // not necssary for klattice
  CalculateSymmetryPointGroupKlattice(FALSE);
}

void CalculateSymmetryPointGroupKlattice(xstructure& str,bool ossverbose,ostream& oss,bool fffverbose) {
  ofstream FileDevNull("/dev/null");
  _aflags aflags;
  aflags.Directory="./";aflags.QUIET=TRUE;
  // str.LatticeReduction_avoid=TRUE;  // not necssary for klattice
  SYM::CalculatePointGroupKlattice(FileDevNull,str,aflags,fffverbose,ossverbose,oss);
}

void CalculateSymmetryPointGroupKlattice(xstructure& str,bool ossverbose,ostream& oss) {
  // str.LatticeReduction_avoid=TRUE;  // not necssary for klattice
  _aflags aflags;
  aflags.Directory="./";aflags.QUIET=TRUE;
  CalculateSymmetryPointGroupKlattice(str,ossverbose,oss,FALSE);
}

void CalculateSymmetryPointGroupKlattice(xstructure& str,bool ossverbose) {
  //  str.LatticeReduction_avoid=TRUE;  // not necssary for klattice
  _aflags aflags;
  aflags.Directory="./";aflags.QUIET=TRUE;
  CalculateSymmetryPointGroupKlattice(str,ossverbose,cout,FALSE);
}

void CalculateSymmetryPointGroupKlattice(xstructure& str) {
  //  str.LatticeReduction_avoid=TRUE;  // not necssary for klattice
  _aflags aflags;
  aflags.Directory="./";aflags.QUIET=TRUE;
  CalculateSymmetryPointGroupKlattice(str,FALSE,cout,FALSE);
}

// ***************************************************************************
// Function buildGenericTitle()
// ***************************************************************************
void xstructure::buildGenericTitle(bool vasp_input,bool force_fix){
  if(vasp_input){pflow::fixEmptyAtomNames(*this,force_fix);}
  title.clear();
  uint iat=0;
  //if any names missing from atoms, lets use generic names
  bool atom_names=true;
  for(uint i=0;i<atoms.size()&&atom_names;i++){if(atoms[i].name.empty()){atom_names=false;}}
  for(uint itype=0;itype<num_each_type.size();itype++){
    for(uint j=0;j<(uint) num_each_type.at(itype);j++) {
      if(j==0){
        if(atom_names){title+=atoms.at(iat).name;}
        else{title+=char('A'+itype);}
        title+=aurostd::utype2string(num_each_type.at(itype));
      }
      iat++;
    }
  }
}

// ***************************************************************************
// Function xstructure2qe
// ***************************************************************************
void xstructure::xstructure2qe(void) {
  ReScale(1.0);
  neg_scale=FALSE;
  coord_flag=_COORDS_FRACTIONAL_; 
  iomode=IOQE_GEOM;
  return;
}

// ***************************************************************************
// Function xstructure2vasp
// ***************************************************************************
void xstructure::xstructure2vasp(void) {
  scale=1.0;
  neg_scale=FALSE;
  coord_flag=_COORDS_FRACTIONAL_; 
  iomode=IOVASP_AUTO;
  //  cerr << "title=\"" << title << "\"" << endl;
  if(title.empty()) {buildGenericTitle();}  //CO 171008 - pushed all of this to a function
  //  uint iat=0;
  //  for(uint itype=0;itype<num_each_type.size();itype++) 
  //    for(uint j=0;j<(uint) num_each_type.at(itype);j++) {
  //      if(j==0) title+=atoms.at(iat).name+aurostd::utype2string(num_each_type.at(itype));
  //      //      atoms.at(iat++).type=itype;
  //    }
  //}
  return;
}

// ***************************************************************************
// Function xstructure2aims
// ***************************************************************************
void xstructure::xstructure2aims(void) {
  ReScale(1.0);
  neg_scale=FALSE;
  coord_flag=_COORDS_FRACTIONAL_; 
  iomode=IOAIMS_GEOM;
  if(title.empty()) {buildGenericTitle();}  //CO 171008
  return;
}

// ***************************************************************************
// Function xstructure2abinit
// ***************************************************************************
void xstructure::xstructure2abinit(void) {
  ReScale(1.0);
  neg_scale=FALSE;
  coord_flag=_COORDS_FRACTIONAL_; 
  iomode=IOABINIT_GEOM;
  return;
}

// ***************************************************************************
// Function platon2print
// ***************************************************************************
string xstructure::platon2print(bool P_EQUAL,bool P_EXACT,double P_ang,double P_d1,double P_d2,double P_d3) {
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  stringstream oss;
  oss.setf(std::ios::fixed,std::ios::floatfield);
  xstructure str=*this;

  // Deal with too small volume problem
  if(str.GetVolume() < (double) PLATON_MIN_VOLUME_PER_ATOM*str.atoms.size()) {
    str.SetVolume(PLATON_MIN_VOLUME_PER_ATOM*str.atoms.size()) ;
    if(LDEBUG) cerr << "platon2print:  PLATON FIXED VOLUME PER ATOM = " << PLATON_MIN_VOLUME_PER_ATOM << endl;
  }
  if(LDEBUG) cerr << "platon2print: volume=" << str.GetVolume() << endl;
  // Set scale to 1 so you don't need to rescale coordinates.
  str.ReScale(1.0);
  std::vector<string> el_names(7);
  int k;
  el_names[0]="Ag";
  el_names[1]="Zr";
  el_names[2]="Cd";
  el_names[3]="Mo";
  el_names[4]="Fe";
  el_names[5]="W";
  el_names[6]="O";
  // Print out data in Platon format

  oss << "TITL " << str.title << endl;
  oss.precision(8);
  oss << "CELL "
      << str.scale*aurostd::modulus(str.lattice(1,1),str.lattice(1,2),str.lattice(1,3)) << " "
      << str.scale*aurostd::modulus(str.lattice(2,1),str.lattice(2,2),str.lattice(2,3)) << " "
      << str.scale*aurostd::modulus(str.lattice(3,1),str.lattice(3,2),str.lattice(3,3)) << " "
      << aurostd::angle(str.lattice(2,1),str.lattice(2,2),str.lattice(2,3),str.lattice(3,1),str.lattice(3,2),str.lattice(3,3)) << " "
      << aurostd::angle(str.lattice(1,1),str.lattice(1,2),str.lattice(1,3),str.lattice(3,1),str.lattice(3,2),str.lattice(3,3)) << " "
      << aurostd::angle(str.lattice(1,1),str.lattice(1,2),str.lattice(1,3),str.lattice(2,1),str.lattice(2,2),str.lattice(2,3))
      << endl;
  // oss << str.atoms.size() << endl;
  oss.precision(10);
  k=0;
  if(str.num_each_type.size()!=2) {
    for(uint i=0;i<str.num_each_type.size();i++)
      for(int j=1;j<=str.num_each_type.at(i);j++) {
	//      printf("%c%i ",i+66,j);
	if(str.atoms.at(k).name_is_given) {
	  oss << str.atoms.at(k).cleanname << j << " ";
	}
	else{ // Must make up a name
	  if(i>(el_names.size()-1)) { // No more default names so make everything W
	    oss << "W" << j << " ";
	  }
	  else{ // Use default names
	    oss << el_names[i] << j << " ";
	  }
	}
	oss << str.atoms.at(k).fpos(1)<< " " << str.atoms.at(k).fpos(2) << " " << str.atoms.at(k).fpos(3) << endl;
	k++;
      }
  }
  if(str.num_each_type.size()==2) {
    for(uint i=0;i<str.num_each_type.size();i++)
      for(int j=1;j<=str.num_each_type.at(i);j++) {
	if(str.num_each_type.at(0)>str.num_each_type.at(1))
	  oss << el_names[1-i] << j << " " << str.atoms.at(k).fpos(1)<< " " << str.atoms.at(k).fpos(2) << " " << str.atoms.at(k).fpos(3) << endl;
	else
	  oss << el_names[i] << j << " " << str.atoms.at(k).fpos(1)<< " " << str.atoms.at(k).fpos(2) << " " << str.atoms.at(k).fpos(3) << endl;
	k++;
      }
  }
  // oss <<  num_each_type.at(i) << " ";
  //  int str.atoms.size()=cpos.size();
  // for(i=0;i<str.atoms.size();i++)
  //  if(str.coord_flag==0)
  // oss << " " << str.atoms.at(i).fpos(1)<< " " << str.atoms.at(i).fpos(2) << " " << stry.atoms.at(i).fpos(3) << endl;
  //  printf("CALC ADDSYM ");
  oss << "CALC ADDSYM ";// << endl;
  if(P_EQUAL) oss << "EQUAL ";
  if(P_EXACT) oss << "EXACT ";
  if(!aurostd::isequal(P_ang,PLATON_TOLERANCE_ANGLE,1.0e-6) ||
     !aurostd::isequal(P_d1,PLATON_TOLERANCE_D1,1.0e-6) ||
     !aurostd::isequal(P_d1,PLATON_TOLERANCE_D2,1.0e-6) ||
     !aurostd::isequal(P_d1,PLATON_TOLERANCE_D3,1.0e-6))
    //    printf("%14.10f %14.10f %14.10f %14.10f ",P_ang,P_d1,P_d2,P_d3);  printf("\n");
    oss << P_ang << " " << P_d1 << " " << P_d2 << " " << P_d3 << " ";//
  oss << endl;
  return oss.str();
}


// ***************************************************************************
// Function FakeNames
// ***************************************************************************
void xstructure::FakeNames(void) {
  // fix names
  int iatom=0;
  for(uint itype=0;itype<num_each_type.size();itype++)
    for(int j=0;j<num_each_type.at(itype);j++) {
      //    if(atoms.at(iatom).name_is_given==FALSE)
      {
	atoms.at(iatom).name=atom_symbol_vec[2+atoms.at(iatom).type];
	if(atoms.at(iatom).type==0) atoms.at(iatom).name="Ag"; // works....
	//	if(atoms.at(iatom).type==1) atoms.at(iatom).name="Au"; // works....
	atoms.at(iatom).CleanName();
	//DX 9/21/17 - Need to keep spin info atoms.at(iatom).CleanSpin();
	atoms.at(iatom).name_is_given=TRUE;
      }
      iatom++;
    }
}

// ***************************************************************************
// Function SpaceGroup
// ***************************************************************************
string xstructure::platon2sg(bool P_EQUAL,bool P_EXACT,double P_ang,double P_d1,double P_d2,double P_d3) {
  xstructure str=*this;
  stringstream aus;
  // string directory="/tmp/_aflow_platon_"+XHOST.ostrPID.str();  // dont change what works
  string directory=XHOST.Tmpfs+"/_aflow_platon_"+XHOST.ostrPID.str();  // dont change what works
  string file=directory+"/aflow_platon_";
  string file_spf=file+".spf";
  string file_out=file+".out";
  string output;
  vector<string> space_group;
  aurostd::DirectoryMake(directory);
  str.FakeNames();
  aus << str.platon2print(P_EQUAL,P_EXACT,P_ang,P_d1,P_d2,P_d3);
  aurostd::stringstream2file(aus,file_spf);aus.clear();aus.str(std::string());
  aus << "cd " << directory << endl;
  aus << "platon -o -c " << file_spf << " | grep \"Space Group\" | head -1 | sed \"s/Space Group //g\"  | sed \"s/No:/#/g\" | sed \"s/,/#/g\" | sed \"s/ //g\" > " << file_out << endl;
  aurostd::execute(aus);
  aurostd::string2tokens(aurostd::file2string(file_out),space_group,"#");
  aus.clear();aus.str(std::string());
  aus << "rm -rf " << directory << endl;
  aurostd::execute(aus);
  spacegroup=space_group[0]+" #"+space_group[1];
  spacegrouplabel=space_group[0];
  spacegroupnumber=aurostd::string2utype<int>(space_group[1]);
  spacegroupoption="FIX: to extract from platon";
  is_spacegroup_platon=TRUE;
  return spacegroup;
}

string xstructure::findsym2sg(double tolerance) {
  xstructure str=*this;
  // Read in input file.
  stringstream oss;
  string findsym_in=aurostd::TmpFileCreate("findsym.in");
  oss << str.findsym2print(tolerance);
  aurostd::stringstream2file(oss,findsym_in);
  vector<string> vlines,tokens;
  FINDSYM::Write("data_space.txt","./");
  FINDSYM::Write("data_wyckoff.txt","./");
  aurostd::string2vectorstring(aurostd::execute2string(XHOST.command("findsym")+" < "+findsym_in),vlines);
  FROZSL::Delete("data_space.txt","./");
  FROZSL::Delete("data_wyckoff.txt","./");
  aurostd::RemoveFile(findsym_in);
  aurostd::RemoveFile("./findsym.log");
  for(uint i=0;i<vlines.size();i++) {
    if(aurostd::substring2bool(vlines.at(i),"Space Group"))
      aurostd::string2tokens(vlines.at(i),tokens);
  }
  //    for(uint i=0;i<tokens.size();i++) cerr << tokens.at(i) << endl;
  spacegroup=tokens.at(4)+" #"+tokens.at(2);
  spacegrouplabel=tokens.at(4);
  spacegroupnumber=aurostd::string2utype<int>(tokens.at(2));
  spacegroupoption="FIX: to extract from findsym";
  is_spacegroup_findsym=TRUE;
  return spacegroup;
}

// ***************************************************************************
// xstructure::findsym2execute
// ***************************************************************************
// This funtion executes findsys
string xstructure::findsym2execute(double tolerance) {
  xstructure str=*this;
  // Read in input file.
  stringstream oss;
  string findsym_in=aurostd::TmpFileCreate("findsym.in");
  oss << str.findsym2print(tolerance);
  aurostd::stringstream2file(oss,findsym_in);
  vector<string> vlines,tokens;
  FINDSYM::Write("data_space.txt","./");
  FINDSYM::Write("data_wyckoff.txt","./");
  string out=aurostd::execute2string(XHOST.command("findsym")+" < "+findsym_in);
  FROZSL::Delete("data_space.txt","./");
  FROZSL::Delete("data_wyckoff.txt","./");
  aurostd::RemoveFile(findsym_in);
  aurostd::RemoveFile("./findsym.log");
  return out;
}

// ***************************************************************************
// xstructure::findsym2print
// ***************************************************************************
// This funtion prints out structural data in a format for findsys
string xstructure::findsym2print(double tolerance) {
  xstructure sstr=*this;
  stringstream oss;
  oss.setf(std::ios::fixed,std::ios::floatfield);
  oss.precision(10);
  // Set scale to 1 so you don't need to rescale coordinates.
  sstr.ReScale(1.0);
  // Print out data in Findsym format
  oss << sstr.title << endl;
  oss << tolerance << "                               accuracy (0=highest=1e-6) " << endl;
  if(0) { // old
    oss << "1                                          form lattice: (1=3 vectors) (2=lengths,angles) " << endl;
    oss << sstr.lattice(1,1) << " " << sstr.lattice(1,2) << " " << sstr.lattice(1,3) << "    Bravais lattice" << endl;
    oss << sstr.lattice(2,1) << " " << sstr.lattice(2,2) << " " << sstr.lattice(2,3) << "    Bravais lattice" << endl;
    oss << sstr.lattice(3,1) << " " << sstr.lattice(3,2) << " " << sstr.lattice(3,3) << "    Bravais lattice" << endl;
    oss << "1.0 0.0 0.0                   unit cell in function of lattice parameter" << endl;
    oss << "0.0 1.0 0.0                   unit cell in function of lattice parameter" << endl;
    oss << "0.0 0.0 1.0                   unit cell in function of lattice parameter" << endl;
    oss << sstr.atoms.size() << "                             number of atoms in the primitive unit cell" << endl;
    for(uint i=0;i<sstr.num_each_type.size();i++)
      for(int j=1;j<=sstr.num_each_type.at(i);j++)
	oss <<  sstr.num_each_type.at(i) << " ";
  }
  if(1) { //new
    oss << "2 Form of lattice parameters: to be entered as lengths and angles " << endl;
    oss << sstr.a << " " << sstr.b << " " << sstr.c << " " << sstr.alpha << " " << sstr.beta << " " << sstr.gamma << " a, b, c, alpha, beta, gamma " << endl;
    oss << "1 Vectors defining the unit cell" << endl;
    oss << "1.0 0.0 0.0                   unit cell in function of lattice parameter" << endl;
    oss << "0.0 1.0 0.0                   unit cell in function of lattice parameter" << endl;
    oss << "0.0 0.0 1.0                   unit cell in function of lattice parameter" << endl;
    oss << sstr.atoms.size() << "                            number of atoms in the primitive unit cell" << endl;
    for(uint i=0;i<sstr.atoms.size();i++)
      oss <<  sstr.atoms.at(i).type+1 << " ";
  }  
  oss << endl;
  for(uint i=0;i<sstr.atoms.size();i++) {
    //  if(coord_flag==0)
    oss << " " << sstr.atoms.at(i).fpos(1) << " " << sstr.atoms.at(i).fpos(2) << " " << sstr.atoms.at(i).fpos(3) << endl;
    //  if(coord_flag==1) oss << " " << sstr.atoms.at(i).cpos(1) << " " << sstr.atoms.at(i).cpos(2) << " " << sstr.atoms.at(i).cpos(3) << endl;
  }
  //  oss << endl;
  return oss.str();
}
// Stefano Curtarolo MIT 2002 - DUKE 2013

// ***************************************************************************
// Rotate
// ***************************************************************************
// Rotation is done around the origin of the structure.
// Define origin q, a point p, and a rotation of p around
// q (R_q(p)).  Then R_q(p)=R_0(p-q)+q=R_0(p)-R_0(q)+q.
// We can evalutate the R_0 terms by simply mulitplying by
// a matrix.  Then we must add q to all the final cartesian
// positions.
xstructure Rotate(const xstructure&a, const xmatrix<double>& rm) {
  // Get R_0(p) for all cartesian positions.
  xstructure b(a);
  xmatrix<double> nlattice(3,3);
  nlattice=trasp(a.lattice);
  b.lattice=trasp(rm*nlattice);
  for(int ia=0;ia<(int)b.atoms.size();ia++)
    b.atoms.at(ia).cpos=F2C(b.lattice,b.atoms.at(ia).fpos);
  //Get R_0(q)
  xvector<double> r_orig(3);
  r_orig=rm*b.origin;
  // Assign new cartesian positions
  for(int ia=0;ia<(int)b.atoms.size();ia++)
    b.atoms.at(ia).cpos=b.atoms.at(ia).cpos-r_orig+b.origin;
  // Get all the direct coords.
  for(int ia=0;ia<(int)b.atoms.size();ia++)
    b.atoms.at(ia).fpos=C2F(b.lattice,b.atoms.at(ia).cpos);
  return b;
}


// ***************************************************************************
//  Function GetLTCell
// ***************************************************************************
// This function operates a linear transform on a cell.
xstructure GetLTCell(const xmatrix<double>& lt, const xstructure& str) {
  xstructure nstr(str);
  nstr.lattice=str.lattice*trasp(lt);
  for(int iat=0;iat<(int)nstr.atoms.size();iat++) {
    //    nstr.atoms.at(iat).cpos=str.atoms.at(iat).cpos*trasp(lt);
    //  nstr.atoms.at(iat).fpos=C2F(nstr.lattice,nstr.atoms.at(iat).cpos);
    nstr.atoms.at(iat)=F2C(nstr.lattice,nstr.atoms.at(iat));
  }
  return nstr;
}

// ***************************************************************************
//  Function GetLTFVCell
// ***************************************************************************
// This function operates a linear transform on a cell.  Rotates all vectors
// by phi around nvec.  Formula from Goldstein, 1980 (Actually used version on
// Mathematic web site).

xstructure GetLTFVCell(const xvector<double>& nvec, const double phi, const xstructure& str) {
  //  xmatrix<double> fpos=str.fpos;
  xmatrix<double> nlat(3,3);
  xvector<double> nhat(3),v1(3),v2(3),v3(3),r(3);
  nhat=nvec/aurostd::modulus(nvec);
  double tmp;
  double cphi=std::cos(phi);
  double sphi=std::sin(phi);
  for(int i=1;i<=3;i++) {
    r=str.lattice(i);
    double nhatpr=scalar_product(nhat,r);
    v1=cphi*r;
    tmp=nhatpr*(1.0-cphi);
    v2=tmp*nhat;
    v3=vector_product(r,nhat);
    v3=sphi*v3;
    for(int j=1;j<=3;j++) nlat(i,j)=nlat(i,j)+v1(j);
    for(int j=1;j<=3;j++) nlat(i,j)=nlat(i,j)+v2(j);
    for(int j=1;j<=3;j++) nlat(i,j)=nlat(i,j)+v3(j);
  }
  xstructure nstr=str;
  nstr.lattice=nlat;
  for(int iat=0;iat<(int)nstr.atoms.size();iat++) {
    nstr.atoms.at(iat).fpos=str.atoms.at(iat).fpos;// *trasp(lt);
    nstr.atoms.at(iat).cpos=F2C(nstr.lattice,nstr.atoms.at(iat).fpos);
  }
  return nstr;
}

// ***************************************************************************
//  Function ShiftPos ShiftCpos ShiftFpos
// ***************************************************************************
// This function shifts the position of the atoms !
xstructure ShiftPos(const xstructure& a,const xvector<double>& shift, const int& flag) {
  xstructure b(a);
  if(flag==FALSE) { // Cartesian shift
    b.coord_flag=TRUE;
    for(uint ia=0;ia<b.atoms.size();ia++) {
      b.atoms.at(ia).cpos=a.atoms.at(ia).cpos+shift;
      b.atoms.at(ia).fpos=C2F(b.lattice,b.atoms.at(ia).cpos);
    }
  }
  if(flag==TRUE) { // Direct coords shift
    b.coord_flag=FALSE;
    for(uint ia=0;ia<b.atoms.size();ia++) {
      b.atoms.at(ia).fpos=a.atoms.at(ia).fpos+shift;
      b.atoms.at(ia).cpos=F2C(b.lattice,b.atoms.at(ia).fpos);
    }
  }
  return b;
}

xstructure ShiftCPos(const xstructure& a,const xvector<double>& shift) {
  xstructure b(a);
  b.coord_flag=TRUE;
  for(uint ia=0;ia<b.atoms.size();ia++) {
    b.atoms.at(ia).cpos=b.atoms.at(ia).cpos+shift;
    b.atoms.at(ia).fpos=C2F(b.lattice,b.atoms.at(ia).cpos);
  }
  return b;
}

xstructure ShiftFPos(const xstructure& a,const xvector<double>& shift) {
  xstructure b(a);
  b.coord_flag=FALSE;
  for(uint ia=0;ia<b.atoms.size();ia++) {
    b.atoms.at(ia).fpos=a.atoms.at(ia).fpos+shift;
    b.atoms.at(ia).cpos=F2C(b.lattice,b.atoms.at(ia).fpos);
  }
  return b;
}

// **************************************************************************
// MaxStructureLattice and MinStructureLattice
// **************************************************************************
double MaxStructureLattice(const xstructure& str) {
  return max(aurostd::modulus(str.lattice(1)),aurostd::modulus(str.lattice(2)),aurostd::modulus(str.lattice(3)));
}

double MinStructureLattice(const xstructure& str) {
  return min(aurostd::modulus(str.lattice(1)),aurostd::modulus(str.lattice(2)),aurostd::modulus(str.lattice(3)));
}

// **************************************************************************
// AtomCDisp
// **************************************************************************
// This function returns the Cartesian displacement vector at2-at1 between two atoms.
xvector<double> AtomCDisp(const _atom& at1, const _atom& at2) {
  xvector<double> diff(3);
  diff=at2.cpos-at1.cpos;
  return diff;
}

// **************************************************************************
// AtomDist
// **************************************************************************
// This function returns the distance between two atoms.
double AtomDist(const xstructure& str,const _atom& atom1,const _atom& atom2) {  // with structure
  xvector<double> pos1(3),pos2(3);
  // contains shift WITH the ijk lattice and origin of structure
  pos1=str.scale*(atom1.cpos+atom1.ijk(1)*str.lattice(1)+atom1.ijk(2)*str.lattice(2)+atom1.ijk(3)*str.lattice(3))+str.origin;
  pos2=str.scale*(atom2.cpos+atom2.ijk(1)*str.lattice(1)+atom2.ijk(2)*str.lattice(2)+atom2.ijk(3)*str.lattice(3))+str.origin;
  return aurostd::modulus(pos1-pos2);
}
double AtomDist(const _atom& at1, const _atom& at2) {  // without structure
  return aurostd::modulus(AtomCDisp(at1,at2));
}

// **************************************************************************
// SameAtom
// **************************************************************************
#define _atom_eps_tol_       0.01
bool SameAtom(const xstructure& str,const _atom& atom1,const _atom& atom2) {
  if(AtomDist(str,atom1,atom2)<_atom_eps_tol_ && atom1.type==atom2.type) return TRUE;
  return FALSE;
}
bool SameAtom(const _atom& atom1,const _atom& atom2) {
  if(AtomDist(atom1,atom2)<_atom_eps_tol_ && atom1.type==atom2.type) return TRUE;
  return FALSE;
}

// **************************************************************************
// DifferentAtom
// **************************************************************************
bool DifferentAtom(const xstructure& str,const _atom& atom1,const _atom& atom2) {
  if(AtomDist(str,atom1,atom2)>_atom_eps_tol_ || atom1.type!=atom2.type) return TRUE;
  return FALSE;
}

// **************************************************************************
// GetDistMatrix
// **************************************************************************
//CO 171024
xmatrix<double> GetDistMatrix(const xstructure& aa){
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  xstructure a(aa); a.ReScale(1.0);
  xstructure xstr_cluster=a;
  //int dim=aurostd::max(LatticeDimensionSphere(a.lattice,RadiusSphereLattice(a.lattice))); //OVERKILL
  GenerateGridAtoms(xstr_cluster,RadiusSphereLattice(a.lattice));//dim,dim,dim);            //OVERKILL 
  vector<vector<int> > atom_types;
  atom_types.push_back(vector<int>(0)); atom_types.back().push_back(0);
  bool found;
  uint indx_found;
  for(uint i=1;i<a.atoms.size();i++){
    found=false;
    for(uint j=0;j<atom_types.size()&&!found;j++){
      if(!atom_types[j].size()){
        cerr << "GetDistMatrix: atom_types indices are compromised!" << endl;
        xmatrix<double> dummy;  //dummy
        return dummy;
      }  //should never come up
      if(a.atoms[i].type==a.atoms[atom_types[j][0]].type){
        found=true;
        indx_found=j;
      }
    }
    if(!found){
      atom_types.push_back(vector<int>(0));
      indx_found=atom_types.size()-1;
    }
    atom_types[indx_found].push_back(i);
  }
  if(a.num_each_type.size()){
    if(atom_types.size()!=a.num_each_type.size()){
      cerr << "GetDistMatrix: odd count of atom types, it's probably poorly sorted" << endl;
      xmatrix<double> dummy;  //dummy
      return dummy;
    }
  }
  if(LDEBUG){
    for(uint it1=0;it1<atom_types.size();it1++){cerr << "GetDistMatrix: atom_types[" << it1 <<"][0]=" << atom_types[it1][0] << endl;}
  }
  xmatrix<double> distsij(1,1,atom_types.size(),atom_types.size());
  for(uint it1=0;it1<atom_types.size();it1++){
    for(uint it2=0;it2<atom_types.size();it2++){
      distsij(it1+1,it2+1)=distsij(it2+1,it1+1)=AUROSTD_MAX_DOUBLE;
    }
  }
  uint atom1;
  double distij,min_dist;
  for(uint it1=0;it1<atom_types.size();it1++){      //type 1 (must be in this order)
    for(uint it2=it1;it2<atom_types.size();it2++){  //type 2 (must be in this order)
      if(LDEBUG){cerr << "GetDistMatrix: finding min dist between itype=" << it1 << " and itype=" << it2 << endl;}
      min_dist=AUROSTD_MAX_DOUBLE;                                      //reset min_dist
      for(uint ia1=0;ia1<atom_types[it1].size();ia1++){                 //must go through all atoms of same type
        atom1=xstr_cluster.grid_atoms_pc2scMap[atom_types[it1][ia1]];   //get respective index in cluster
        for(uint ia2=0;ia2<atom_types[it2].size();ia2++){               //must go through all atoms of the same type
          for(uint atom2=0;atom2<(uint)xstr_cluster.grid_atoms_number;atom2++){ //go through all atoms of the cluster
            if(atom1!=atom2 && a.atoms[atom_types[it2][ia2]].type==xstr_cluster.grid_atoms[atom2].type){  //cannot be same index (dist=0), and must be the types we want
              distij=AtomDist(xstr_cluster.grid_atoms[atom1],xstr_cluster.grid_atoms[atom2]); //distance
              if(0&&LDEBUG){
                cerr << "GetDistMatrix: trying dist=" << distij << " between atom[" << atom1 << ",type="<<xstr_cluster.grid_atoms[atom1].type<<"] and atom[" << atom2 << ",type="<<xstr_cluster.grid_atoms[atom2].type<<"]" << endl;
              }
              if(distij<min_dist){  //grab min_dist
                distsij(it1+1,it2+1)=distsij(it2+1,it1+1)=distij;
                min_dist=distij;
                if(LDEBUG){
                  cerr << "GetDistMatrix: found new min dist=" << min_dist << " between atom[" << atom1 << ",type="<<xstr_cluster.grid_atoms[atom1].type<<"] and atom[" << atom2 << ",type="<<xstr_cluster.grid_atoms[atom2].type<<"]" << endl;
                }
              }
            }
          }
        }
      }
    }
  }
  return distsij;
}

// **************************************************************************
// GetNBONDXX
// **************************************************************************
//CO 171024
vector<double> GetNBONDXX(const xstructure& a){
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  xmatrix<double> distsij=GetDistMatrix(a);
  vector<double> dists;
  for(uint it1=1;it1<(uint)distsij.rows+1;it1++){
    for(uint it2=it1;it2<(uint)distsij.cols+1;it2++){
      dists.push_back(distsij(it1,it2));
    }
  }
  if(LDEBUG){
    cerr << "GetNBONDXX: DIST_MATRIX" << endl;
    cerr << distsij << endl;
    cerr << "GetNBONDXX: MIN_DIST=" << SYM::minimumDistance(a) << endl;
  }
  return dists;
}

// **************************************************************************
// GenerateGridAtoms
// **************************************************************************
// make grid of atoms !
int xstructure::GenerateGridAtoms(int i1,int i2,int j1,int j2,int k1,int k2) {
  int GenerateGridAtoms(xstructure& str,int i1,int i2,int j1,int j2,int k1,int k2);
  return GenerateGridAtoms(*this,i1,i2,j1,j2,k1,k2);
}

int xstructure::GenerateGridAtoms(int d1,int d2,int d3) {
  int GenerateGridAtoms(xstructure& str,int i1,int i2,int j1,int j2,int k1,int k2);
  return GenerateGridAtoms(*this,-d1,d1,-d2,d2,-d3,d3);
}

int xstructure::GenerateGridAtoms(int d) {
  int GenerateGridAtoms(xstructure& str,int i1,int i2,int j1,int j2,int k1,int k2);
  return GenerateGridAtoms(*this,-d,d,-d,d,-d,d);
}

int GenerateGridAtoms(xstructure& str,int i1,int i2,int j1,int j2,int k1,int k2) {
  // same scale as before
  str.grid_atoms.clear();
  str.grid_atoms_sc2pcMap.clear(); str.grid_atoms_pc2scMap.clear();
  _atom atom;
  str.BringInCell();  // are INCELL.
  xvector<double> a1(3),a2(3),a3(3);                     // a1,a2,a3 are the rows of the lattice matrix
  a1=str.lattice(1);a2=str.lattice(2);a3=str.lattice(3); // a1,a2,a3 are the rows of the lattice matrix
  for(uint iat=0;iat<str.atoms.size();iat++){
    str.grid_atoms.push_back(str.atoms.at(iat));  // put first the unit cell !
    str.grid_atoms_pc2scMap.push_back(str.grid_atoms.size()-1); //CO 171025 
    str.grid_atoms_sc2pcMap.push_back(iat); //CO 171025
  }
  for(int i=i1;i<=i2;i++) {
    for(int j=j1;j<=j2;j++) {
      for(int k=k1;k<=k2;k++) {
	if(i!=0 || j!=0 || k!=0) {
	  for(uint iat=0;iat<str.atoms.size();iat++) {
	    atom=str.atoms.at(iat);
	    atom.isincell=FALSE; // these are OUT OF CELL
	    atom.cpos=((double)i)*a1+((double)j)*a2+((double)k)*a3+str.atoms.at(iat).cpos;
	    atom.fpos[1]=i+str.atoms.at(iat).fpos[1];
	    atom.fpos[2]=j+str.atoms.at(iat).fpos[2];
	    atom.fpos[3]=k+str.atoms.at(iat).fpos[3];
	    str.grid_atoms.push_back(atom);
            str.grid_atoms_sc2pcMap.push_back(iat); //CO 171025
	  }
	}
      }
    }
  }
  str.grid_atoms_calculated=TRUE;
  str.grid_atoms_dimsL[1]=i1;str.grid_atoms_dimsL[2]=j1;str.grid_atoms_dimsL[3]=k1;
  str.grid_atoms_dimsH[1]=i2;str.grid_atoms_dimsH[2]=j2;str.grid_atoms_dimsH[3]=k2;
  str.grid_atoms_number=str.grid_atoms.size();
  //  for(uint i=0;i<str.grid_atoms.size();i++)
  //   cerr << str.grid_atoms.at(i).cpos << endl;
  // cerr << str.grid_atoms.size() << endl;
  return str.grid_atoms.size();
}

int GenerateGridAtoms(xstructure& str,const double& radius) {
  xvector<int> dims(3);
  dims=LatticeDimensionSphere(str,radius);  // radius is not normalized over the scale
  return GenerateGridAtoms(str,-dims(1),dims(1),-dims(2),dims(2),-dims(3),dims(3));
}

int GenerateGridAtoms(xstructure& str,int d1,int d2,int d3) {
  return GenerateGridAtoms(str,-d1,d1,-d2,d2,-d3,d3);
}

int GenerateGridAtoms(xstructure& str,int d) {
  return GenerateGridAtoms(str,-d,d,-d,d,-d,d);
}

int GenerateGridAtoms(xstructure& str,const xvector<int>& dims) {
  return GenerateGridAtoms(str,-dims(1),dims(1),-dims(2),dims(2),-dims(3),dims(3));
}

int GenerateGridAtoms(xstructure& str) {
  return GenerateGridAtoms(str,-1,1,-1,1,-1,1);
}

// **************************************************************************
// GenerateLIJK table stuff
// **************************************************************************
namespace aurostd {   // INT
  class _ssort_int_value0123 {                    // sorting through reference
  public:
    bool operator()(const vector<int>& v1, const vector<int>& v2) const {
      if(v1[0]<v2[0]) return TRUE;
      if(v1[0]==v2[0]) {
	if(v1[1]>v2[1]) return TRUE;
	if(v1[1]==v2[1]) {
	  if(v1[2]>v2[2]) return TRUE;
	  if(v1[2]==v2[2]) {
	    if(v1[3]>v2[3]) return TRUE;
	  }
	}
      }
      return FALSE;
    }
  };
}

int xstructure::GenerateLIJK(double radius) {
  lijk_table.clear();
  xvector<double> a1(3),a2(3),a3(3); //    a1,a2,a3 are the rows of the lattice matrix
  a1=lattice(1);a2=lattice(2);a3=lattice(3);    // a1,a2,a3 are the rows of the lattice matrix
  vector<int> _lijk(4);
  xvector<int> int_ijk(3);
  xvector<double> cpos(3),fpos(3);
  xvector<int> dims(3);
  dims=LatticeDimensionSphere(lattice,radius);  // for(int i=0;i<=3;i++) dims[i]+=1;
  lijk_dims[1]=dims[1];
  lijk_dims[2]=dims[2];
  lijk_dims[3]=dims[3];
  vector<vector<int> > _temp_lijk_table;
  _temp_lijk_table.clear();

  for(int i=-lijk_dims[1];i<=lijk_dims[1];i++) {
    for(int j=-lijk_dims[2];j<=lijk_dims[2];j++) {
      for(int k=-lijk_dims[3];k<=lijk_dims[3];k++) {
	cpos=((double)i)*a1+((double)j)*a2+((double)k)*a3;
	// 	if(aurostd::modulus(cpos)<=radius)  // I put them all since I will not be scanning through.
	{
	  _lijk[0]=(int) (10000000.0*sqrt(aurostd::modulus(cpos))); // for the next sorting
	  _lijk[1]=i;_lijk[2]=j;_lijk[3]=k;
	  _temp_lijk_table.push_back(_lijk);
	}
      }
    }
  }
  sort(_temp_lijk_table.begin(),_temp_lijk_table.end(),aurostd::_ssort_int_value0123());

  for(uint i=0;i<_temp_lijk_table.size();i++) {
    _temp_lijk_table[i][0]=i; // relabel
    int_ijk[1]=_temp_lijk_table.at(i)[1];
    int_ijk[2]=_temp_lijk_table.at(i)[2];
    int_ijk[3]=_temp_lijk_table.at(i)[3];
    lijk_table.push_back(int_ijk);
  }

  for(uint i=0;i<lijk_table.size();i++) {
    cpos=((double)lijk_table.at(i)[1])*a1+((double)lijk_table.at(i)[2])*a2+((double)lijk_table.at(i)[3])*a3;
    fpos[1]=lijk_table.at(i)[1];fpos[2]=lijk_table.at(i)[2];fpos[3]=lijk_table.at(i)[3];
    lijk_cpos.push_back(cpos);
    lijk_fpos.push_back(fpos);
  }

  lijk_calculated=TRUE;
  return lijk_table.size();
}


// **************************************************************************
// LIJK table: L => IJK
// **************************************************************************
void l2ijk(const xstructure& str,const int &l,int &i,int &j,int &k) {
  if(l<0) {cerr << "l2ijk error: l<0 " << endl; exit(0);}
  if(l>(int) str.lijk_table.size()) {cerr << "l2ijk error: >str.lijk_table.size() " << endl; exit(0);}
  i=str.lijk_table[l][1];j=str.lijk_table[l][2];k=str.lijk_table[l][3];
}

void l2ijk(const xstructure& str,const int &l,xvector<int>& ijk) {
  if(l<0) {cerr << "l2ijk error: l<0 " << endl; exit(0);}
  if(l>(int) str.lijk_table.size()) {cerr << "l2ijk error: >str.lijk_table.size() " << endl; exit(0);}
  ijk[1]=str.lijk_table[l][1];ijk[2]=str.lijk_table[l][2];ijk[3]=str.lijk_table[l][3];
}

xvector<int> l2ijk(const xstructure& str,const int &l) {
  xvector<int> ijk(3);
  int i,j,k;
  l2ijk(str,l,i,j,k);
  ijk(1)=i;ijk(2)=j;ijk(3)=k;
  return ijk;
}


// **************************************************************************
// LIJK table: IJK => L
// **************************************************************************
void ijk2l(const xstructure& str,int &l,const int &i,const int &j,const int &k) {
  l=-1;
  if(i<-str.lijk_dims(1) || i>str.lijk_dims(1)) {cerr << "ijk2l error; i out of boundary" << endl;exit(0);}
  if(j<-str.lijk_dims(2) || j>str.lijk_dims(2)) {cerr << "ijk2l error; j out of boundary" << endl;exit(0);}
  if(k<-str.lijk_dims(3) || k>str.lijk_dims(3)) {cerr << "ijk2l error; k out of boundary" << endl;exit(0);}
  for(uint ll=0;ll<str.lijk_table.size();ll++) // start search
    if(str.lijk_table.at(ll)[1]==i)         // faster comparison one at a time
      if(str.lijk_table.at(ll)[2]==j)       // faster comparison one at a time
	if(str.lijk_table.at(ll)[3]==k)     // faster comparison one at a time
	  l=ll;
  if(l<0) {cerr << "ijk2l error: l not found" << endl; exit(0);}
}

void ijk2l(const xstructure& str,int &l,const xvector<int>& ijk) {
  ijk2l(str,l,ijk(1),ijk(2),ijk(3));
}


int ijk2l(const xstructure& str,const int &i,const int &j,const int &k) {
  int l;
  ijk2l(str,l,i,j,k);
  return l;
}

int ijk2l(const xstructure& str,const xvector<int>& ijk) {
  int l;
  ijk2l(str,l,ijk(1),ijk(2),ijk(3));
  return l;
}

// **************************************************************************
// AQEgeom2aims AQEgeom2abinit AQEgeom2qe AQEgeom2vasp
// ***************************************************************************
xstructure AQEgeom2aims(istream& input) {
  xstructure a(input,IOAFLOW_AUTO);
  //  if(a.iomode==IOVASP_AUTO || a.iomode==IOVASP_POSCAR || a.iomode==IOVASP_ABCCAR || a.iomode==IOVASP_WYCKCAR) 
  a.xstructure2aims();
  return a;
}

xstructure AQEgeom2abinit(istream& input) {
  xstructure a(input,IOAFLOW_AUTO);
  //  if(a.iomode==IOVASP_AUTO || a.iomode==IOVASP_POSCAR || a.iomode==IOVASP_ABCCAR || a.iomode==IOVASP_WYCKCAR) 
  a.xstructure2abinit();
  return a;
}

xstructure AQEgeom2qe(istream& input) {
  xstructure a(input,IOAFLOW_AUTO);
  //  if(a.iomode==IOVASP_AUTO || a.iomode==IOVASP_POSCAR || a.iomode==IOVASP_ABCCAR || a.iomode==IOVASP_WYCKCAR) 
  a.xstructure2qe();
  return a;
}

xstructure AQEgeom2vasp(istream& input) {
  xstructure a(input,IOAFLOW_AUTO);
  //  if(a.iomode==IOQE_AUTO || a.iomode==IOQE_GEOM)
  a.xstructure2vasp();
  //  cerr << a.title << endl;
  return a;
}

// **************************************************************************
// r_lattice
// ***************************************************************************
xvector<double> r_lattice(const xstructure& str,const int &l) {
  xvector<double> rrr(3);
  int i,j,k;
  l2ijk(str,l,i,j,k);
  rrr=((double)i)*str.lattice(1)+((double)j)*str.lattice(2)+((double)k)*str.lattice(3);
  return rrr;
}

xvector<double> r_lattice(const xstructure& str,const int &i,const int &j,const int &k) {
  xvector<double> rrr(3);
  rrr=((double)i)*str.lattice(1)+((double)j)*str.lattice(2)+((double)k)*str.lattice(3);
  return rrr;
}

xvector<double> r_lattice(const xstructure& str,const xvector<int>& ijk) {
  xvector<double> rrr(3);
  rrr=((double)ijk(1))*str.lattice(1)+((double)ijk(2))*str.lattice(2)+((double)ijk(3))*str.lattice(3);
  return rrr;
}

// **************************************************************************
// Function GetUnitCellRep
// **************************************************************************
// Dane Morgan / Stefano Curtarolo
// This function calculated the representation of a point in terms of a position in cell 0 and its true unit cell.
// Output position in cell 0 is given with same coordinate type (Cart or Direct) as input position.
// Note that there is always an ambiguity of being in a given cell with position 0, or one cell over with position 1.  This
// is broken by forcing all cell positions to be at 0 if they are within TOL of 1.  This gives consistent cell image locations.
// CO 170717 - Compare this with PBC() function, this function brings in to 0th cell, PBC brings in to cell between -0.5 and 0.5
// PBC() is good for minimizing overall fpos
void GetUnitCellRep(const xvector<double>& ppos,xvector<double>& p_cell0,
		    xvector<int>& ijk,const xmatrix<double>& lattice,
		    const bool coord_flag) {
  double TOL=1e-11;
  xvector<double> cpos(3),fpos(3);
  // p_cell0=xvector<double> (3);
  // ijk=xvector<int> (3);
  
  clear(ijk);clear(p_cell0);

  if(coord_flag==TRUE) // ppos is in cartesian coords
    fpos=C2F(lattice,ppos);
  
  for(int ic=1;ic<=3;ic++) {
    int s;
    s=SignNoZero(fpos(ic));
    if(s>0) ijk(ic)=int(fpos(ic));
    else    ijk(ic)=int(fpos(ic))-1;
    p_cell0(ic)=fpos(ic)-ijk(ic);
    // Boundary Correction.  p_cell0 is >=0 and <1.  
    // If p_cell0 is within TOL of 1
    // then set it to 0 and shift ijk up by 1.
    if(abs(1-p_cell0(ic))<TOL) {
      p_cell0(ic)=0;
      ijk(ic)=ijk(ic)+1;
    }
  }
  if(coord_flag==TRUE) // ppos is in cartesian coords
    p_cell0=F2C(lattice,p_cell0);
}

// **************************************************************************
// Function COMPARE_GetNeighData
// **************************************************************************
// For sort algorithm in GetNeighData
class compare_GetNeighData {
public:  
  int operator()(const _atom& a, const _atom& b) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "operator()(const _atom& a, const _atom& b)" << endl;
    double tol=1e-15;
    // Sort by distance
    if(aurostd::isequal(GetDistFromOrigin(a),GetDistFromOrigin(b),tol)) {
      if(LDEBUG) cerr << "operator()(const _atom& a, const _atom& b)  isequal(GetDistFromOrigin(a),GetDistFromOrigin(b),tol)" << endl;
      // Sort by unit cell values
      if(a.name==b.name) {
	if(LDEBUG) cerr << "operator()(const _atom& a, const _atom& b)  a.name==b.name" << endl;
        xvector<int> ijka=a.ijk;
        xvector<int> ijkb=b.ijk;
        int va=100*ijka(1)+10*ijka(2)+1*ijka(3);
        int vb=100*ijkb(1)+10*ijkb(2)+1*ijkb(3);
        return va<vb;
      } else {
	if(LDEBUG) cerr << "operator()(const _atom& a, const _atom& b)  a.name!=b.name" << endl;
	return a.name<b.name;
      }
    } else{
      if(LDEBUG) cerr << "operator()(const _atom& a, const _atom& b)  !isequal(GetDistFromOrigin(a),GetDistFromOrigin(b),tol)" << endl;
      return GetDistFromOrigin(a)<GetDistFromOrigin(b);
    }
    // Sort by name
    // if(a.name==b.name) {
    //  return GetDistFromOrigin(a)<GetDistFromOrigin(b);
    // }
    // else{
    //  return a.name<b.name;
    // }
  }
};

// **************************************************************************
// Function GetNeighData
// **************************************************************************
// This function collects all the neighbor data between rmin and rmax and stores 
// it for each atom in a vector of atom objects in order of increasing distance.  

void xstructure::GetNeighData(const deque<_atom>& in_atom_vec,		
			      const double& rmin, const double& rmax,
			      deque<deque<_atom> >& neigh_mat) {
  double epsilon=1.0e-3;  // sometimes you have wrong images due to roundoff
  deque<_atom> neigh_vec;
  // Get data from str.
  xstructure sstr(*this);
  // Set scale to 1 so you don't need to rescale coordinates.
  // [OBSOLETE] sstr=ReScale(sstr,1.0);
  sstr.ReScale(1.0);
  xmatrix<double> original_lattice(3,3);
  original_lattice=sstr.lattice;
  // Use Niggli cell to avoid missing neighbors in angular cells.
  // This involves converting your atoms into the 0th Niggli cell,
  // keeping track of the unit cell you came from, and then shifting back
  // at the end.  This shifting must be done for both atom_vec (at the
  // beginning) and neigh_mat (at the end).
  //sstr.NiggliUnitCellForm();  
  sstr.MinkowskiBasisReduction();
  //    cerr << sstr << endl;
  
  xmatrix<double> lat(3,3);
  lat=sstr.lattice;
  
  //  Convert input atom_vec to Niggli lattice
  deque<_atom> atom_vec=in_atom_vec;
  xmatrix<int>  atom_shifts((int) atom_vec.size()-1,3,0,1);
  for(uint ia=0;ia<atom_vec.size();ia++) {
    _atom a=atom_vec.at(ia);
    xvector<double> p_cell0(3);
    xvector<int>    atom_shifts_temp(3);
    GetUnitCellRep(a.cpos,p_cell0,atom_shifts_temp,lat,1);
    atom_vec.at(ia).cpos=p_cell0;
    atom_vec.at(ia)=ConvertAtomToLat(atom_vec.at(ia),lat);
    atom_shifts(ia,1)=atom_shifts_temp(1);
    atom_shifts(ia,2)=atom_shifts_temp(2);
    atom_shifts(ia,3)=atom_shifts_temp(3);
  }    
  
  // Create vector of all atoms in all periodic images needed so that for atoms in the unit cell every possible neighbor within rmax
  // is included.  This can be done by looping over all atoms in all unit cells such that if the unit cells are given by
  // n0*v0,n1*v1,n2*v2, then ni_max*vi>rmax for all i=0->2.  This is not rigorous but seems to work quite well.  It failed for very
  // angular cells, but by using the Niggli reduced cell this pitfall seems to be avoided (but no proof, so be careful).
  
  int imax,jmax,kmax;
  // Find imax
  // (algorithm 1) approximate
  imax=(int)max(rmax/aurostd::modulus(lat(1)),rmax/aurostd::modulus(lat(2)),rmax/aurostd::modulus(lat(3)))+2;
  jmax=imax;
  kmax=imax;
  // (algorithm 2) exact, the requirement is that atoms must in incell.
  // we should find nlat1,nlat2,nlat3
  // where nlat1 is the lat1 component that is perpendicular to the plane made by lat2 & lat3.
  // then imax=ceil(rmax/nlat1)
  // This algorithm is implemented in function LatticeDimensionSphere
  xvector<int> dims(3);
  dims=LatticeDimensionSphere(lat,rmax);
  imax=dims(1);
  jmax=dims(2);
  kmax=dims(3);
  xvector<int> ijk(3);
  deque<_atom> all_atom_vec;
  // latt maybe a rotated version of POSCAR, so need to incell-ized the fpos and cpos
  // sstr.BringInCell(); // with roundoff
  sstr.BringInCell(-1.0); // no roundoff
  for(ijk(1)=-imax;ijk(1)<=imax;ijk(1)++) {
    for(ijk(2)=-jmax;ijk(2)<=jmax;ijk(2)++) {
      for(ijk(3)=-kmax;ijk(3)<=kmax;ijk(3)++) {
	xvector<double> ctpos(3);
	xvector<double> ftpos(3);
	for(uint iat=0;iat<sstr.atoms.size();iat++) {
	  _atom a;
	  a=sstr.atoms.at(iat);
	  a.name=sstr.atoms.at(iat).name;
	  a.number=iat;
	  a.ijk=ijk;
	  for(int ic=1;ic<=3;ic++)
	    ctpos(ic)=sstr.atoms.at(iat).cpos(ic)+ijk(1)*lat(1,ic)+ijk(2)*lat(2,ic)+ijk(3)*lat(3,ic);
	  ftpos=C2F(lat,ctpos);
	  a.cpos=ctpos;
	  a.fpos=ftpos;
	  a.type=sstr.atoms.at(iat).type;
	  all_atom_vec.push_back(a);
	} // iat
      } // ijk
    } // ijk
  } // ijk
  
  // Now build neighbors list for each atom on atom list.  Each
  // neighbor list will be row in the neigh_mat matrix.
  for(uint ia1=0;ia1<atom_vec.size();ia1++) {
    xvector<double> corigin(3);
    _atom at=atom_vec.at(ia1);
    corigin = at.cpos;
    at.corigin=corigin;
    deque<_atom> neigh_vec;
    neigh_vec.push_back(at);
    for(uint ia2=0;ia2<all_atom_vec.size();ia2++) {
      double dist = AtomDist(at,all_atom_vec.at(ia2));
      if(dist<=rmax && dist>=rmin && dist>=epsilon) {
	all_atom_vec.at(ia2).corigin=corigin;
	neigh_vec.push_back(all_atom_vec.at(ia2));
      } // if
    } // ia2
    sort(neigh_vec.begin()+1,neigh_vec.end(),compare_GetNeighData());
    neigh_mat.push_back(neigh_vec);
  } // ia1
  
    //  Convert neigh_mat from Niggli cell to original lattice
  for(uint ia=0;ia<neigh_mat.size();ia++) {
    for(uint ja=0;ja<neigh_mat.at(ia).size();ja++) {
      xvector<double> fpos(3);
      xvector<double> cpos(3);
      _atom a;
      a=neigh_mat.at(ia).at(ja);
      fpos=a.fpos;
      for(int ic=1;ic<=3;ic++) {
	fpos(ic)=fpos(ic)+atom_shifts(ia,ic);
      }
      cpos=F2C(lat,fpos);
      a.cpos=cpos;
      neigh_mat.at(ia).at(ja)=ConvertAtomToLat(a,original_lattice);
    }
  }
}

// **************************************************************************
// Function GetStrNeighData
// **************************************************************************
// This function collects all the neighbor data out to some
// cutoff and stores it for each atom in the structure.
void xstructure::GetStrNeighData(const double cutoff,deque<deque<_atom> >& neigh_mat) {
  deque<_atom> atom_vec;
  neigh_mat.clear();
  // Get data from str.
  // Set scale to 1 so you don't need to rescale coordinates.
  xstructure sstr(*this);
  // [OBSOLETE]  sstr=ReScale(sstr,1.0);
  sstr.ReScale(1.0);
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
  // [OBSOLETE] GetNeighData(atom_vec,sstr,rmin,cutoff,neigh_mat);
  sstr.GetNeighData(atom_vec,rmin,cutoff,neigh_mat);
}

// **************************************************************************
// GenerateXXXXXXXXXXx
// **************************************************************************
void xstructure::qm_clear(void) {
  qm_calculated=FALSE;
  qm_scale=0.0;
  qm_lattice.clear();
  qm_origin.clear();
  qm_klattice.clear();
  qm_f2c.clear();
  qm_c2f.clear();
  qm_origin.clear();
  qm_forces_write=FALSE;
  qm_positions_write=FALSE;
  qm_atoms.clear();
  qm_forces.clear();
  qm_positions.clear();
  for(uint iat=0;iat<atoms.size();iat++) {
    qm_atoms.push_back(atoms.at(iat));          // just plug something, then I`ll clean
    qm_atoms.at(iat).cpos.clear();
    qm_atoms.at(iat).fpos.clear();
    qm_forces.push_back(atoms.at(iat).cpos);    // just plug something, then I`ll clean
    qm_forces.at(iat).clear();
    qm_positions.push_back(atoms.at(iat).cpos); // just plug something, then I`ll clean
    qm_positions.at(iat).clear();
  }
  qm_E_cell=0.0;qm_dE_cell=0.0;qm_H_cell=0.0;qm_PV_cell=0.0;qm_mag_cell=0.0;qm_P=0.0;
  qm_E_atom=0.0;qm_dE_atom=0.0;qm_H_atom=0.0;qm_PV_atom=0.0;qm_mag_atom=0.0; 
 
  if(atoms.size()!=qm_atoms.size())     {cerr << "ERROR qm_clear [1] atoms.size()!=qm_atoms.size()"<<endl;exit(0);}
  if(atoms.size()!=qm_forces.size())    {cerr << "ERROR qm_clear [2] atoms.size()!=qm_forces.size()"<<endl;exit(0);}
  if(atoms.size()!=qm_positions.size()) {cerr << "ERROR qm_clear [3] atoms.size()!=qm_positions.size()"<<endl;exit(0);}

}

void xstructure::qm_recycle(void) {
  if(atoms.size()!=qm_atoms.size())     {cerr << "ERROR qm_load [1] atoms.size()!=qm_atoms.size()"<<endl;exit(0);}
  if(atoms.size()!=qm_forces.size())    {cerr << "ERROR qm_load [2] atoms.size()!=qm_forces.size()"<<endl;exit(0);}
  if(atoms.size()!=qm_positions.size()) {cerr << "ERROR qm_load [3] atoms.size()!=qm_positions.size()"<<endl;exit(0);}
  scale=qm_scale;
  lattice=qm_lattice;
  klattice=qm_klattice;
  f2c=qm_f2c;
  c2f=qm_c2f;
  origin=qm_origin;
  FixLattices();
  for(uint i=0;i<atoms.size();i++) {                 // copy cpos/fpos
    atoms.at(i).cpos=qm_atoms.at(i).cpos;     // get cpos
    atoms.at(i).fpos=qm_atoms.at(i).fpos;     // get fpos
  }
  qm_clear();
}

void xstructure::qm_load(string Directory,string suffix,int iomode) {
  double data_natoms=double(atoms.size());
  if(iomode!=IOVASP_POSCAR) {cerr << "ERROR xtructure::qm_load, only IOVASP_POSCAR is supported" << endl;exit(0);};
  if(iomode==IOVASP_POSCAR) {
    xOUTCAR outcar;
    //if(aurostd::FileEmpty(Directory+"/OUTCAR"))   {cerr << "ERROR xtructure::qm_load: Empty OUTCAR" << endl;exit(0);}
    if(aurostd::FileEmpty(Directory+"/OUTCAR"+suffix))   {cerr << "ERROR xtructure::qm_load: Empty OUTCAR" << endl;exit(0);}  // PINKU AND JOSE FIXED BUG
 
    outcar.GetPropertiesFile(Directory+"/OUTCAR"+suffix,data_natoms,TRUE);

    if(abs(data_natoms-(double) outcar.natoms)>0.1) { 
      cerr << "ERROR void xstructure::qm_load: data_natoms(" << data_natoms << ")!= (int) outcar.natoms(" << outcar.natoms << ") ..." << endl;
      cerr << "      Directory=" << Directory << endl;
      cerr << "      suffix=" << suffix << endl;
      cerr << "      iomode=" << iomode << endl;
      exit(0);}

    //    cerr << atoms.size() << endl;
    qm_clear();
    if(atoms.size()!=qm_atoms.size())     {cerr << "ERROR xtructure::qm_load [1] atoms.size()!=qm_atoms.size()" << endl;exit(0);}
    if(atoms.size()!=qm_forces.size())    {cerr << "ERROR xtructure::qm_load [2] atoms.size()!=qm_forces.size()" << endl;exit(0);}
    if(atoms.size()!=qm_positions.size()) {cerr << "ERROR xtructure::qm_load [3] atoms.size()!=qm_positions.size()" << endl;exit(0);}

    // NEW WITH xOUTCAR
    qm_forces.clear(); for(uint i=0;i<outcar.vforces.size();i++)  qm_forces.push_back(outcar.vforces.at(i)); 
    qm_positions.clear(); for(uint i=0;i<outcar.vpositions_cartesian.size();i++)  qm_positions.push_back(outcar.vpositions_cartesian.at(i)); 
    if(atoms.size()!=qm_forces.size())    {cerr << "ERROR xtructure::qm_load [4] atoms.size()!=qm_forces.size()" << endl;exit(0);}
    if(atoms.size()!=qm_positions.size()) {cerr << "ERROR xtructure::qm_load [5] atoms.size()!=qm_positions.size()" << endl;exit(0);}

    // NEW WITH xVASPRUNXML
    xVASPRUNXML vasprunxml;
    //if(aurostd::FileEmpty(Directory+"/vasprun.xml"))   {cerr << "ERROR xtructure::qm_load: Empty vasprun.xml" << endl;exit(0);}
    if(aurostd::FileEmpty(Directory+"/vasprun.xml"+suffix))   {cerr << "ERROR xtructure::qm_load: Empty vasprun.xml" << endl;exit(0);} //PINKU AND JOSE FIXED BUG
    //vasprunxml.GetPropertiesFile(Directory+"/vasprun.xml");
    vasprunxml.GetPropertiesFile(Directory+"/vasprun.xml"+suffix); //PINKU AND JOSE FIXED BUG
    qm_forces.clear(); for(uint i=0;i<vasprunxml.vforces.size();i++)  qm_forces.push_back(vasprunxml.vforces.at(i)); 
   
    
    // OLD
    stringstream CONTCAR_stringstream;
    CONTCAR_stringstream.str(std::string()); CONTCAR_stringstream << aurostd::FileToString(Directory+"/CONTCAR"+suffix);
     
    qm_E_cell=outcar.energy_cell;qm_E_atom=outcar.energy_atom;
    qm_H_cell=outcar.enthalpy_cell;qm_H_atom=outcar.enthalpy_atom;
    qm_PV_cell=outcar.PV_cell;qm_PV_atom=outcar.PV_atom;
    qm_P=outcar.pressure;
    qm_mag_cell=outcar.mag_cell;qm_mag_atom=outcar.mag_atom;
    
    // CONTCAR OPERATIONS ---------------------------------------------------------------
    // cerr << "DEBUG: xstructure::qm_clear(void)" << endl;
    // cerr << CONTCAR_stringstream.str() << endl;
    // cerr << "DEBUG: xstructure::qm_clear(void)" << endl;
    // string CONTCAR_string="";
    xstructure b(CONTCAR_stringstream);                               // LOAD it in
    qm_scale=b.scale;
    qm_lattice=b.lattice;
    qm_klattice=b.klattice;
    qm_f2c=b.f2c;
    qm_c2f=b.c2f;
    qm_origin=b.origin;
    for(uint i=0;i<atoms.size();i++) {                 // copy cpos/fpos
      qm_atoms.at(i)=atoms.at(i);
      qm_atoms.at(i).cpos=b.atoms.at(i).cpos;     // get cpos
      qm_atoms.at(i).fpos=b.atoms.at(i).fpos;     // get fpos
    }
    qm_calculated=TRUE;
  }
}

// **************************************************************************
// should be in xlibs
// SEARCH SORT WORLD "sort for vector/deque of double_int"

// namespace aurostd {
//   void sort(vector<double>& varg1,vector<int>& varg2) {
//     vector<aurostd::_double_int_> vv(varg1.size());
//     for(uint i=0;i<varg1.size();i++) {vv.at(i).arg1=varg1.at(i);vv.at(i).arg2=varg2.at(i);}
//     sort(vv.begin(),vv.end(),_sort_double_int_());
//     for(uint i=0;i<varg1.size();i++) {varg1.at(i)=vv.at(i).arg1;varg2.at(i)=vv.at(i).arg2;}
//   }
//   void sort(deque<double>& varg1,deque<int>& varg2) {
//     deque<aurostd::_double_int_> vv(varg1.size());
//     for(uint i=0;i<varg1.size();i++) {vv.at(i).arg1=varg1.at(i);vv.at(i).arg2=varg2.at(i);}
//     sort(vv.begin(),vv.end(),_sort_double_int_());
//     for(uint i=0;i<varg1.size();i++) {varg1.at(i)=vv.at(i).arg1;varg2.at(i)=vv.at(i).arg2;}
//   }
// }

// namespace aurostd {
//   void sort(vector<double>& varg1,vector<double>& varg2) {
//     vector<aurostd::_double_double_> vv(varg1.size());
//     for(uint i=0;i<varg1.size();i++) {vv.at(i).arg1=varg1.at(i);vv.at(i).arg2=varg2.at(i);}
//     sort(vv.begin(),vv.end(),_sort_double_double_());
//     for(uint i=0;i<varg1.size();i++) {varg1.at(i)=vv.at(i).arg1;varg2.at(i)=vv.at(i).arg2;}
//   }
//   void sort(deque<double>& varg1,deque<double>& varg2) {
//     deque<aurostd::_double_double_> vv(varg1.size());
//     for(uint i=0;i<varg1.size();i++) {vv.at(i).arg1=varg1.at(i);vv.at(i).arg2=varg2.at(i);}
//     sort(vv.begin(),vv.end(),_sort_double_double_());
//     for(uint i=0;i<varg1.size();i++) {varg1.at(i)=vv.at(i).arg1;varg2.at(i)=vv.at(i).arg2;}
//   }
// }

// **************************************************************************
// *                                                                        *
// *             STEFANO CURTAROLO - Duke University 2003-2018              *
// *                                                                        *
// **************************************************************************


// **************************************************************************
// PAULING DETECTOR
// **************************************************************************
bool PAULING_WyckoffDetector(vector<string> &vinput) {
  vector<string> tokens,elements;
  aurostd::string2tokens(vinput.at(0),tokens);
  for(uint i=1;i<tokens.size();i++)
    elements.push_back(tokens.at(i));
  for(uint i=0;i<elements.size();i++) {
    //  cout << "scanning " << elements.at(i) << endl;
    for(uint j=1;j<vinput.size();j++) {
      aurostd::string2tokens(vinput.at(j),tokens);
      if(tokens.size()==7) {
	if(tokens.at(1)==elements.at(i)) {
	  cout << "   ";
	  for(uint k=0;k<3;k++) cout << aurostd::PaddedPOST(tokens.at(4+k),7," ") << "  ";
	  cout << char(i+65) << "   "; // shift to ascii
	  cout << aurostd::PaddedPRE(tokens.at(0),7," ") << "  ";
	  for(uint k=1;k<4;k++) cout << aurostd::PaddedPOST(tokens.at(k),4," ") << "  ";
	  cout << endl;
	}
      }
    }
  }

  exit(0);
}

//DX 8/31/17 - xstructure2json - START
// **************************************************************************
// xstructure2json
// **************************************************************************
string xstructure2json(xstructure& xstr) {
  string eendl="";
  bool roff=true; //round off
  bool PRINT_NULL=FALSE;
  stringstream sss;
  stringstream sscontent_json;
  vector<string> vcontent_json;

  // TITLE 
  if(xstr.title.size()){
    sscontent_json << "\"title\":\"" << xstr.title << "\"" << eendl;
  } else{
    if(PRINT_NULL){ sscontent_json << "\"title\":null" << eendl;}
  }
  vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

  // SCALE/VOL 
  if(xstr.scale){
    sscontent_json << "\"scale\":\"" << xstr.scale << "\"" << eendl;
  } else{
    if(PRINT_NULL){ sscontent_json << "\"scale\":null" << eendl;}
  }
  vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
  
  // LATTICE
  if(xstr.lattice.rows){
    sscontent_json << "\"lattice\":[" << aurostd::xmatDouble2String(xstr.lattice,5,roff) << "]" << eendl;
  } else{
    if(PRINT_NULL){ sscontent_json << "\"lattice\":null" << eendl;}
  }
  vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
 
  // SPECIES 
  if(xstr.species.size()){
    sscontent_json << "\"species\":[" << aurostd::joinWDelimiter(aurostd::wrapDeqEntries(xstr.species,"\""),",") << "]" << eendl;
  } else{
    if(PRINT_NULL){ sscontent_json << "\"species\":null" << eendl;}
  }
  vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

  // NUMBER OF TYPES
  if(xstr.num_each_type.size()){
    sscontent_json << "\"number_each_type\":[" << aurostd::joinWDelimiter(xstr.num_each_type,",") << "]" << eendl;
  } else{
    if(PRINT_NULL){ sscontent_json << "\"number_each_type\":null" << eendl;}
  }
  vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

  // COORDINATES TYPE
  if(xstr.coord_flag==_COORDS_FRACTIONAL_){
    sscontent_json << "\"coordinates_type\":\"direct\"" << eendl; //CO 171025
  } else if(xstr.coord_flag==_COORDS_CARTESIAN_){
    sscontent_json << "\"coordinates_type\":\"Cartesian\"" << eendl;
  }
  vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

  // ATOMS
  if(xstr.atoms.size()){
    sscontent_json << "\"atoms\":[" << eendl;
    for(uint i=0;i<xstr.atoms.size();i++){
      sscontent_json << atom2json(xstr.atoms[i],xstr.coord_flag,xstr.partial_occupation_flag) << eendl;
      if(i != xstr.atoms.size()-1){
        sscontent_json << ",";
      }
    }
    sscontent_json << "]" << eendl;
  }
  vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

  sss << "{" << aurostd::joinWDelimiter(vcontent_json,",")  << "}" << eendl;
  return sss.str();
}
//DX 8/31/17 - xstructure2json - END

//DX 8/31/17 - atom2json - START
// **************************************************************************
// atom2json
// **************************************************************************
string atom2json(_atom& atom, int coord_flag, int poccupation) {
  string eendl="";
  bool roff=true; //round off
  bool PRINT_NULL=FALSE;
  stringstream sss;
  stringstream sscontent_json;
  vector<string> vcontent_json;

  // NAME 
  if(atom.name.size()){
    sscontent_json << "\"name\":\"" << atom.name << "\"" << eendl;
  } else{
    if(PRINT_NULL){ sscontent_json << "\"name\":null" << eendl;}
  }
  vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

  // POSITION
  if(coord_flag==_COORDS_FRACTIONAL_){
    sscontent_json << "\"position\":[" << aurostd::joinWDelimiter(xvecDouble2vecString(atom.fpos,5,roff),",") << "]" << eendl;
  } else if(coord_flag==_COORDS_CARTESIAN_){
    sscontent_json << "\"position\":[" << aurostd::joinWDelimiter(xvecDouble2vecString(atom.cpos,5,roff),",") << "]" << eendl;
  } else{
    if(PRINT_NULL){ sscontent_json << "\"position\":null" << eendl;}
  }
  vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
  
  // OCCUPANCY
  if(poccupation==TRUE){
    sscontent_json << "\"occupancy\":" << atom.partial_occupation_value << eendl;
  } else if(poccupation==FALSE){
    sscontent_json << "\"occupancy\":1.0"<< eendl;
  } else{
    if(PRINT_NULL){ sscontent_json << "\"occupancy\":null" << eendl;}
  }
  vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

  sss << "{" << aurostd::joinWDelimiter(vcontent_json,",")  << "}";
  return sss.str();
} 
//DX 8/31/17 - atom2json

#endif  // _AFLOW_XATOM_CPP
