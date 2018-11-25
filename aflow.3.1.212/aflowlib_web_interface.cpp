// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *                                                                         *
// ***************************************************************************
// Stefano Curtarolo
// fixed for XZ - SC

#ifndef _AFLOWLIB_WEB_INTERFACE_CPP_
#define _AFLOWLIB_WEB_INTERFACE_CPP_
#include "aflow.h"
#include "aflowlib_webapp_entry.cpp"  // CO 170622 - bob hanson JMOL stuff
#include "aflowlib_webapp_bands.cpp"  // CO 180305 - geena bands stuff

const string _DEVIL_PROTOTYPES_STRING_ = "64,65,549,550,f8269,f9083,f8819";

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

// ***************************************************************************
namespace aflowlib {
  //  class _aflowlib_entry 
  
  _aflowlib_entry::_aflowlib_entry() {  // constructor PUBLIC
    entry.clear();ventry.clear();
    auid.clear();
    aurl.clear();vaurl.clear();
    keywords.clear();vkeywords.clear();
    aflowlib_date.clear();
    aflowlib_version.clear();
    aflowlib_entries.clear();vaflowlib_entries.clear();
    aflowlib_entries_number=0;
    aflow_version.clear();
    catalog.clear();
    data_api="aapi1.2"; // new version of the API
    data_source="aflowlib";
    data_language="";
    error_status.clear();
    author.clear();vauthor.clear();
    calculation_cores=1;
    calculation_memory=AUROSTD_NAN;
    calculation_time=AUROSTD_NAN;
    corresponding.clear();vcorresponding.clear();
    loop.clear();vloop.clear();
    node_CPU_Cores=AUROSTD_NAN;node_CPU_MHz=AUROSTD_NAN;node_CPU_Model.clear();node_RAM_GB=AUROSTD_NAN;
    Bravais_lattice_orig.clear();Bravais_lattice_relax.clear();
    code.clear();
    composition.clear();vcomposition.clear();
    compound.clear();
    density=AUROSTD_NAN;
    dft_type.clear();vdft_type.clear();
    eentropy_cell=AUROSTD_NAN;eentropy_atom=AUROSTD_NAN;
    Egap=AUROSTD_NAN;Egap_fit=AUROSTD_NAN;
    energy_cell=AUROSTD_NAN;energy_atom=AUROSTD_NAN;energy_atom_relax1=AUROSTD_NAN;
    energy_cutoff=AUROSTD_NAN;
    delta_electronic_energy_convergence=AUROSTD_NAN;
    delta_electronic_energy_threshold=AUROSTD_NAN;
    nkpoints=0;
    nkpoints_irreducible=0;
    kppra=0;
    kpoints.clear();
    kpoints_nnn_relax.clear();
    kpoints_nnn_static.clear();
    kpoints_pairs.clear();
    kpoints_bands_path_grid=0;
    enthalpy_cell=AUROSTD_NAN;enthalpy_atom=AUROSTD_NAN;
    enthalpy_formation_cell=AUROSTD_NAN;enthalpy_formation_atom=AUROSTD_NAN;
    entropic_temperature=AUROSTD_NAN;
    files.clear();vfiles.clear();
    files_LIB.clear();vfiles_LIB.clear();
    files_RAW.clear();vfiles_RAW.clear();
    files_WEB.clear();vfiles_WEB.clear();
    forces.clear();vforces.clear();
    Egap_type.clear();
    geometry.clear();vgeometry.clear();
    lattice_system_orig.clear();lattice_variation_orig.clear();lattice_system_relax.clear();lattice_variation_relax.clear();
    ldau_TLUJ.clear();
    natoms=AUROSTD_NAN;
    nbondxx.clear();vnbondxx.clear();
    nspecies=AUROSTD_NAN;
    Pearson_symbol_orig.clear();Pearson_symbol_relax.clear();
    positions_cartesian.clear();vpositions_cartesian.clear();
    positions_fractional.clear();vpositions_fractional.clear();
    pressure=AUROSTD_NAN;
    stress_tensor.clear();vstress_tensor.clear();
    pressure_residual=AUROSTD_NAN;
    Pulay_stress=AUROSTD_NAN;
    prototype.clear();
    PV_cell=AUROSTD_NAN;PV_atom=AUROSTD_NAN;
    scintillation_attenuation_length=AUROSTD_NAN;
    sg.clear();sg2.clear();vsg.clear();vsg2.clear();  // CO 171202
    spacegroup_orig.clear();spacegroup_relax.clear();
    species.clear();vspecies.clear();
    species_pp.clear();vspecies_pp.clear();
    species_pp_version.clear();vspecies_pp_version.clear();
    species_pp_ZVAL.clear();vspecies_pp_ZVAL.clear();
    spin_cell=AUROSTD_NAN;spin_atom=AUROSTD_NAN;
    spinD.clear();vspinD.clear();
    spinD_magmom_orig.clear();vspinD_magmom_orig.clear();
    spinF=AUROSTD_NAN;
    sponsor.clear();vsponsor.clear();
    stoichiometry.clear();vstoichiometry.clear();
    valence_cell_std=AUROSTD_NAN;valence_cell_iupac=AUROSTD_NAN;
    volume_cell=AUROSTD_NAN;volume_atom=AUROSTD_NAN;
    //DX 20180823 - added more symmetry info - START
    // SYMMETRY
    crystal_family="";
    crystal_system="";
    crystal_class="";
    point_group_Hermann_Mauguin="";
    point_group_Schoenflies="";
    point_group_orbifold="";
    point_group_type="";
    point_group_order=AUROSTD_NAN;
    point_group_structure="";
    Bravais_lattice_lattice_type="";
    Bravais_lattice_lattice_variation_type="";
    Bravais_lattice_lattice_system="";
    Bravais_superlattice_lattice_type="";
    Bravais_superlattice_lattice_variation_type="";
    Bravais_superlattice_lattice_system="";
    Pearson_symbol_superlattice="";
    reciprocal_lattice_type="";
    reciprocal_lattice_variation_type="";
    reciprocal_geometry.clear();vreciprocal_geometry.clear();
    reciprocal_volume_cell=AUROSTD_NAN;
    Wyckoff_letters="";
    Wyckoff_multiplicities="";
    Wyckoff_site_symmetries="";
    //DX 20180823 - added more symmetry info - END
    // AGL/AEL
    agl_thermal_conductivity_300K=AUROSTD_NAN;
    agl_debye=AUROSTD_NAN;
    agl_acoustic_debye=AUROSTD_NAN;
    agl_gruneisen=AUROSTD_NAN;
    agl_heat_capacity_Cv_300K=AUROSTD_NAN;
    agl_heat_capacity_Cp_300K=AUROSTD_NAN;
    agl_thermal_expansion_300K=AUROSTD_NAN;
    agl_bulk_modulus_static_300K=AUROSTD_NAN;
    agl_bulk_modulus_isothermal_300K=AUROSTD_NAN;
    ael_poisson_ratio=AUROSTD_NAN;
    ael_bulk_modulus_voigt=AUROSTD_NAN;
    ael_bulk_modulus_reuss=AUROSTD_NAN;
    ael_shear_modulus_voigt=AUROSTD_NAN;
    ael_shear_modulus_reuss=AUROSTD_NAN;
    ael_bulk_modulus_vrh=AUROSTD_NAN;
    ael_shear_modulus_vrh=AUROSTD_NAN;
    ael_elastic_anistropy=AUROSTD_NAN;
    // BADER
    bader_net_charges.clear();vbader_net_charges.clear();
    bader_atomic_volumes.clear();vbader_atomic_volumes.clear();
    // legacy
    server.clear();vserver.clear();vserverdir.clear();
    icsd.clear();
    stoich.clear();vstoich.clear();
    // apennsy
    structure_name.clear();  // apennsy
    structure_description.clear();  // apennsy
    distance_gnd=AUROSTD_NAN;  // apennsy
    distance_tie=AUROSTD_NAN;  // apennsy
    pureA=FALSE;pureB=FALSE;  // apennsy
    fcc=FALSE; bcc=FALSE;hcp=FALSE;  // apennsy
    stoich_a=AUROSTD_NAN;stoich_b=AUROSTD_NAN;  // apennsy
    bond_aa=AUROSTD_NAN;bond_ab=AUROSTD_NAN;bond_bb=AUROSTD_NAN;  // apennsy
    vNsgroup.clear();  // apennsy
    vsgroup.clear();  // apennsy
    vstr.clear();  // apennsy
  }
  _aflowlib_entry::~_aflowlib_entry() { // destructor PUBLIC
    free();
  }
  
  void _aflowlib_entry::copy(const _aflowlib_entry& b) { // copy PRIVATE
    entry=b.entry;ventry.clear();for(uint i=0;i<b.ventry.size();i++) ventry.push_back(b.ventry.at(i));
    auid=b.auid;
    aurl=b.aurl;vaurl.clear();for(uint i=0;i<b.vaurl.size();i++) vaurl.push_back(b.vaurl.at(i));
    vkeywords.clear();for(uint i=0;i<b.vkeywords.size();i++) vkeywords.push_back(b.vkeywords.at(i));
    aflowlib_date=b.aflowlib_date;
    aflowlib_version=b.aflowlib_version;
    aflowlib_entries=b.aflowlib_entries;
    vaflowlib_entries.clear();for(uint i=0;i<b.vaflowlib_entries.size();i++) vaflowlib_entries.push_back(b.vaflowlib_entries.at(i));
    aflowlib_entries_number=b.aflowlib_entries_number;
    aflow_version=b.aflow_version;
    catalog=b.catalog;
    data_api=b.data_api;
    data_source=b.data_source;
    data_language=b.data_language;
    error_status=b.error_status;
    author=b.author;vauthor.clear();for(uint i=0;i<b.vauthor.size();i++) vauthor.push_back(b.vauthor.at(i));
    calculation_cores=b.calculation_cores;calculation_memory=b.calculation_memory;calculation_time=b.calculation_time;
    corresponding=b.corresponding;vcorresponding.clear();for(uint i=0;i<b.vcorresponding.size();i++) vcorresponding.push_back(b.vcorresponding.at(i));
    loop=b.loop;vloop.clear();for(uint i=0;i<b.vloop.size();i++) vloop.push_back(b.vloop.at(i));
    node_CPU_Cores=b.node_CPU_Cores;node_CPU_MHz=b.node_CPU_MHz;node_CPU_Model=b.node_CPU_Model;node_RAM_GB=b.node_RAM_GB;
    Bravais_lattice_orig=b.Bravais_lattice_orig;Bravais_lattice_relax=b.Bravais_lattice_relax;
    code=b.code;
    composition=b.composition;vcomposition.clear();for(uint i=0;i<b.vcomposition.size();i++) vcomposition.push_back(b.vcomposition.at(i));
    compound=b.compound;
    density=b.density;
    dft_type=b.dft_type;vdft_type.clear();for(uint i=0;i<b.vdft_type.size();i++) vdft_type.push_back(b.vdft_type.at(i));
    eentropy_cell=b.eentropy_cell;eentropy_atom=b.eentropy_atom;
    Egap=b.Egap;Egap_fit=b.Egap_fit;
    energy_cell=b.energy_cell;energy_atom=b.energy_atom;energy_atom_relax1=b.energy_atom_relax1;
    energy_cutoff=b.energy_cutoff;
    delta_electronic_energy_convergence=b.delta_electronic_energy_convergence;
    delta_electronic_energy_threshold=b.delta_electronic_energy_threshold;
    nkpoints=b.nkpoints;
    nkpoints_irreducible=b.nkpoints_irreducible;
    kppra=b.kppra;
    kpoints=b.kpoints;
    kpoints_nnn_relax=b.kpoints_nnn_relax;
    kpoints_nnn_static=b.kpoints_nnn_static;
    kpoints_pairs=b.kpoints_pairs;
    kpoints_bands_path_grid=b.kpoints_bands_path_grid;
    enthalpy_cell=b.enthalpy_cell;enthalpy_atom=b.enthalpy_atom;
    enthalpy_formation_cell=b.enthalpy_formation_cell;enthalpy_formation_atom=b.enthalpy_formation_atom;
    entropic_temperature=b.entropic_temperature;
    files=b.files;vfiles.clear();for(uint i=0;i<b.vfiles.size();i++) vfiles.push_back(b.vfiles.at(i));
    files_LIB=b.files_LIB;vfiles_LIB.clear();for(uint i=0;i<b.vfiles_LIB.size();i++) vfiles_LIB.push_back(b.vfiles_LIB.at(i));
    files_RAW=b.files_RAW;vfiles_RAW.clear();for(uint i=0;i<b.vfiles_RAW.size();i++) vfiles_RAW.push_back(b.vfiles_RAW.at(i));
    files_WEB=b.files_WEB;vfiles_WEB.clear();for(uint i=0;i<b.vfiles_WEB.size();i++) vfiles_WEB.push_back(b.vfiles_WEB.at(i));
    forces=b.forces;vforces.clear();for(uint i=0;i<b.vforces.size();i++) vforces.push_back(b.vforces.at(i));
    Egap_type=b.Egap_type;
    geometry=b.geometry;vgeometry.clear();for(uint i=0;i<b.vgeometry.size();i++) vgeometry.push_back(b.vgeometry.at(i));
    lattice_system_orig=b.lattice_system_orig;lattice_variation_orig=b.lattice_variation_orig;
    lattice_system_relax=b.lattice_system_relax;lattice_variation_relax=b.lattice_variation_relax;
    ldau_TLUJ=b.ldau_TLUJ;
    natoms=b.natoms;
    nbondxx=b.nbondxx;vnbondxx.clear();for(uint i=0;i<b.vnbondxx.size();i++) vnbondxx.push_back(b.vnbondxx.at(i));
    nspecies=b.nspecies;
    Pearson_symbol_orig=b.Pearson_symbol_orig;Pearson_symbol_relax=b.Pearson_symbol_relax;
    positions_cartesian=b.positions_cartesian;vpositions_cartesian.clear();for(uint i=0;i<b.vpositions_cartesian.size();i++) vpositions_cartesian.push_back(b.vpositions_cartesian.at(i));
    positions_fractional=b.positions_fractional;vpositions_fractional.clear();for(uint i=0;i<b.vpositions_fractional.size();i++) vpositions_fractional.push_back(b.vpositions_fractional.at(i));
    pressure=b.pressure;
    stress_tensor=b.stress_tensor;vstress_tensor.clear();for(uint i=0;i<b.vstress_tensor.size();i++) vstress_tensor.push_back(b.vstress_tensor.at(i));
    pressure_residual=b.pressure_residual;
    Pulay_stress=b.Pulay_stress;
    prototype=b.prototype;
    PV_cell=b.PV_cell;PV_atom=b.PV_atom;
    scintillation_attenuation_length=b.scintillation_attenuation_length;
    sg=b.sg;sg2=b.sg2;vsg.clear();for(uint i=0;i<b.vsg.size();i++){vsg.push_back(b.vsg[i]);} vsg2.clear();for(uint i=0;i<b.vsg2.size();i++){vsg2.push_back(b.vsg2[i]);}  // CO 171202
    spacegroup_orig=b.spacegroup_orig;spacegroup_relax=b.spacegroup_relax;
    species=b.species;vspecies.clear();for(uint i=0;i<b.vspecies.size();i++) vspecies.push_back(b.vspecies.at(i));
    species_pp=b.species_pp;vspecies_pp.clear();for(uint i=0;i<b.vspecies_pp.size();i++) vspecies_pp.push_back(b.vspecies_pp.at(i));
    species_pp_version=b.species_pp_version;vspecies_pp_version.clear();for(uint i=0;i<b.vspecies_pp_version.size();i++) vspecies_pp_version.push_back(b.vspecies_pp_version.at(i));
    species_pp_ZVAL=b.species_pp_ZVAL;vspecies_pp_ZVAL.clear();for(uint i=0;i<b.vspecies_pp_ZVAL.size();i++) vspecies_pp_ZVAL.push_back(b.vspecies_pp_ZVAL.at(i));
    spin_cell=b.spin_cell;spin_atom=b.spin_atom;
    spinD=b.spinD;vspinD.clear();for(uint i=0;i<b.vspinD.size();i++) vspinD.push_back(b.vspinD.at(i));
    spinD_magmom_orig=b.spinD_magmom_orig;vspinD_magmom_orig.clear();for(uint i=0;i<b.vspinD_magmom_orig.size();i++) vspinD_magmom_orig.push_back(b.vspinD_magmom_orig.at(i));
    spinF=b.spinF;
    sponsor=b.sponsor;vsponsor.clear();for(uint i=0;i<b.vsponsor.size();i++) vsponsor.push_back(b.vsponsor.at(i));
    stoichiometry=b.stoichiometry;vstoichiometry.clear();for(uint i=0;i<b.vstoichiometry.size();i++) vstoichiometry.push_back(b.vstoichiometry.at(i));
    valence_cell_std=b.valence_cell_std;valence_cell_iupac=b.valence_cell_iupac;
    volume_cell=b.volume_cell;volume_atom=b.volume_atom;
    //DX 20180823 - added more symmetry info - START
    // SYMMETRY
    crystal_family=b.crystal_family;
    crystal_system=b.crystal_system;
    crystal_class=b.crystal_class;
    point_group_Hermann_Mauguin=b.point_group_Hermann_Mauguin;
    point_group_Schoenflies=b.point_group_Schoenflies;
    point_group_orbifold=b.point_group_orbifold;
    point_group_type=b.point_group_type;
    point_group_order=b.point_group_order;
    point_group_structure=b.point_group_structure;
    Bravais_lattice_lattice_type=b.Bravais_lattice_lattice_type;
    Bravais_lattice_lattice_variation_type=b.Bravais_lattice_lattice_variation_type;
    Bravais_lattice_lattice_system=b.Bravais_lattice_lattice_system;
    Bravais_superlattice_lattice_type=b.Bravais_superlattice_lattice_type;
    Bravais_superlattice_lattice_variation_type=b.Bravais_superlattice_lattice_variation_type;
    Bravais_superlattice_lattice_system=b.Bravais_superlattice_lattice_system;
    Pearson_symbol_superlattice=b.Pearson_symbol_superlattice;
    reciprocal_geometry=b.reciprocal_geometry;vreciprocal_geometry.clear();for(uint i=0;i<b.vreciprocal_geometry.size();i++) vreciprocal_geometry.push_back(b.vreciprocal_geometry.at(i));
    reciprocal_volume_cell=b.volume_cell;
    reciprocal_lattice_type=b.reciprocal_lattice_type;
    reciprocal_lattice_variation_type=b.reciprocal_lattice_variation_type;
    Wyckoff_letters=b.Wyckoff_letters;
    Wyckoff_multiplicities=b.Wyckoff_multiplicities;
    Wyckoff_site_symmetries=b.Wyckoff_site_symmetries;
    //DX 20180823 - added more symmetry info - END
    // AGL/AEL
    agl_thermal_conductivity_300K=b.agl_thermal_conductivity_300K;
    agl_debye=b.agl_debye;
    agl_acoustic_debye=b.agl_acoustic_debye;
    agl_gruneisen=b.agl_gruneisen;
    agl_heat_capacity_Cv_300K=b.agl_heat_capacity_Cv_300K;
    agl_heat_capacity_Cp_300K=b.agl_heat_capacity_Cp_300K;
    agl_thermal_expansion_300K=b.agl_thermal_expansion_300K;
    agl_bulk_modulus_static_300K=b.agl_bulk_modulus_static_300K;
    agl_bulk_modulus_isothermal_300K=b.agl_bulk_modulus_isothermal_300K;
    ael_poisson_ratio=b.ael_poisson_ratio;
    ael_bulk_modulus_voigt=b.ael_bulk_modulus_voigt;
    ael_bulk_modulus_reuss=b.ael_bulk_modulus_reuss;
    ael_shear_modulus_voigt=b.ael_shear_modulus_voigt;
    ael_shear_modulus_reuss=b.ael_shear_modulus_reuss;
    ael_bulk_modulus_vrh=b.ael_bulk_modulus_vrh;
    ael_shear_modulus_vrh=b.ael_shear_modulus_vrh;
    ael_elastic_anistropy=b.ael_elastic_anistropy;
    // BADER
    bader_net_charges=b.bader_net_charges;vbader_net_charges.clear();for(uint i=0;i<b.vbader_net_charges.size();i++) vbader_net_charges.push_back(b.vbader_net_charges.at(i));
    bader_atomic_volumes=b.bader_atomic_volumes;vbader_atomic_volumes.clear();for(uint i=0;i<b.vbader_atomic_volumes.size();i++) vbader_atomic_volumes.push_back(b.vbader_atomic_volumes.at(i));
    // legacy
    server=b.server;
    vserver.clear();for(uint i=0;i<b.vserver.size();i++) vserver.push_back(b.vserver.at(i));
    vserverdir.clear();for(uint i=0;i<b.vserverdir.size();i++) vserverdir.push_back(b.vserverdir.at(i));
    icsd=b.icsd;
    stoich=b.stoich;vstoich.clear();for(uint i=0;i<b.vstoich.size();i++) vstoich.push_back(b.vstoich.at(i));
    // apennsy
    structure_name=b.structure_name;  // apennsy
    structure_description=b.structure_description;  // apennsy
    distance_gnd=b.distance_gnd;  // apennsy
    distance_tie=b.distance_tie;  // apennsy
    pureA=b.pureA;pureB=b.pureB;  // apennsy
    fcc=b.fcc;bcc=b.bcc;hcp=b.hcp;  // apennsy
    stoich_a=b.stoich_a;stoich_b=b.stoich_b;  // apennsy
    bond_aa=b.bond_aa;bond_ab=b.bond_ab;bond_bb=b.bond_bb;  // apennsy
    vNsgroup.clear();for(uint i=0;i<b.vNsgroup.size();i++) vNsgroup.push_back(b.vNsgroup.at(i));  // apennsy
    vsgroup.clear();for(uint i=0;i<b.vsgroup.size();i++) vsgroup.push_back(b.vsgroup.at(i));  // apennsy
    vstr.clear();for(uint i=0;i<b.vstr.size();i++) vstr.push_back(b.vstr.at(i));  // apennsy
  }
  

  const _aflowlib_entry& _aflowlib_entry::operator=(const _aflowlib_entry& b) {  // operator= PUBLIC
    if(this!=&b) {free(); copy(b);}
    return *this;
  }
  
  _aflowlib_entry::_aflowlib_entry(const _aflowlib_entry& b) { // copy PUBLIC
    //  free();*this=b;
    copy(b);
  }
  
  void _aflowlib_entry::free() { // free PRIVATE
    ventry.clear();
    vaurl.clear();
    vaflowlib_entries.clear();
    vkeywords.clear();
    vauthor.clear();
    vcorresponding.clear();
    vloop.clear();
    vcomposition.clear(); // clear all vectors
    vdft_type.clear(); // clear all vectors
    vfiles.clear(); // clear all vectors
    vfiles_LIB.clear(); // clear all vectors
    vfiles_RAW.clear(); // clear all vectors
    vfiles_WEB.clear(); // clear all vectors
    vforces.clear(); // clear all vectors
    vgeometry.clear(); // clear all vectors
    vstress_tensor.clear(); // clear all vectors
    vnbondxx.clear(); // clear all vectors
    vpositions_cartesian.clear(); // clear all vectors
    vpositions_fractional.clear(); // clear all vectors
    vspecies.clear(); // clear all vectors
    vspecies_pp.clear(); // clear all vectors
    vspecies_pp_version.clear(); // clear all vectors
    vspecies_pp_ZVAL.clear(); // clear all vectors
    vspinD.clear(); // clear all vectors
    vspinD_magmom_orig.clear(); // clear all vectors
    vsponsor.clear();
    vstoichiometry.clear(); // clear all vectors
    vreciprocal_geometry.clear(); // clear all vectors //DX 20180824 - added reciprocal lattice parameters
    // BADER
    vbader_net_charges.clear();
    vbader_atomic_volumes.clear();
    // legacy
    vserver.clear();  // clear all vectors
    for(uint i=0;i<vserverdir.size();i++)
      vserverdir.at(i).clear();
    vserverdir.clear();  // clear all vectors
    vstoich.clear(); // clear all vectors
  } 
  
  void _aflowlib_entry::clear() {  // clear PRIVATE
    _aflowlib_entry _temp;
    copy(_temp);
  }
  
  // file2aflowlib
  uint _aflowlib_entry::file2aflowlib(const string& file,ostream& oss) {
    if(!aurostd::FileExist(file)) {cerr << "ERROR - _aflowlib_entry::file2aflowlib: " << DEFAULT_FILE_AFLOWLIB_ENTRY_OUT << " not found =" << file << endl;return 0;} //exit(0);} // CO 170609, this is a dud
    string entry;
    aurostd::efile2string(file,entry);
    return Load(entry,oss);
  }

  // Load
  uint _aflowlib_entry::Load(const stringstream& stream,ostream& oss) {
    return Load(stream.str(),oss);
  }

  // LoadWeb
  uint _aflowlib_entry::url2aflowlib(const string& _url,ostream& oss,bool verbose) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string soliloquy="_aflowlib_entry::url2aflowlib():";
    string url=_url;
    if(url.empty()) {cerr << "ERROR - _aflowlib_entry::url2aflowlib: url.empty()" << endl;return 0;} //exit(0);} // CO 170609, this is a dud
    string entry;
    if(aurostd::substring2bool(url,"index") || aurostd::substring2bool(url,"format")) {
      aurostd::StringSubst(url,"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT,"");
      if(!aurostd::url2string(url,entry,verbose)){return 0;}   //corey, this is a dud
    } else {
      aurostd::StringSubst(url,"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT,"");
      if(!aurostd::url2string(url+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT,entry,verbose)){return 0;}  //corey, this is a dud
    }
    if(LDEBUG){cerr << soliloquy << " entry=" << entry << endl;} //CO 180528
    return Load(entry,oss);
  }
  
  // Load overload
  uint _aflowlib_entry::Load(const string& _entry,ostream& oss) {
    clear(); // start from clean
    entry=_entry; // start from loading it up !
    if(entry.empty()) {cerr << "ERROR - _aflowlib_entry::Load: entry.empty()" << endl;return 0;} //exit(0);}  // CO 170609, this is a dud 
    vector<string> tokens,stokens;
    string keyword,content,line;
    aurostd::string2tokens(entry,ventry,"|");
    for(uint i=0;i<ventry.size();i++) {
      line=aurostd::RemoveWhiteSpaces(ventry.at(i));
      aurostd::string2tokens(line,tokens,"=");
      if(tokens.size()>0) {
	keyword=tokens.at(0);
	if(tokens.size()>1) {content=tokens.at(1);} else {continue;} //{content="";}  // CO 180319, content="" screws up string2double(), better to leave as AUROSTD_NAN
        if(content.empty()){continue;}  // CO 180319
        if(content=="null"){continue;}  // CO 180319 - aflux integration!
	aurostd::string2tokens(content,stokens,",");
	if(keyword=="auid") {auid=content;}
  // CO 180409 - added the else if's for speed, no need to go through more checks than necessary
  else if(keyword=="aurl") {aurl=content;aurostd::string2tokens(content,stokens,":");for(uint j=0;j<stokens.size();j++) vaurl.push_back(stokens.at(j));}
  else if(keyword=="keywords") {keywords=content;aurostd::string2tokens(content,stokens,",");for(uint j=0;j<stokens.size();j++) vkeywords.push_back(stokens.at(j));}
  else if(keyword=="aflowlib_date") {aflowlib_date=content;}
  else if(keyword=="aflowlib_version") {aflowlib_version=content;}
  else if(keyword=="aflowlib_entries") {aflowlib_entries=content;aurostd::string2tokens(content,stokens,",");for(uint j=0;j<stokens.size();j++) vaflowlib_entries.push_back(stokens.at(j));}
  else if(keyword=="aflowlib_entries_number") {aflowlib_entries_number=aurostd::string2utype<int>(content);}
  else if(keyword=="aflow_version") {aflow_version=content;}
  else if(keyword=="catalog") {catalog=content;}
  else if(keyword=="data_api") {data_api=content;}
	else if(keyword=="data_source") {data_source=content;}
	else if(keyword=="data_language") {data_language=content;}
	else if(keyword=="error_status") {error_status=content;}
	else if(keyword=="author") {author=content;aurostd::string2tokens(content,stokens,",");for(uint j=0;j<stokens.size();j++) vauthor.push_back(stokens.at(j));}
	else if(keyword=="calculation_cores") {calculation_cores=aurostd::string2utype<int>(content);}
	else if(keyword=="calculation_memory") {calculation_memory=aurostd::string2utype<double>(content);}
	else if(keyword=="calculation_time") {calculation_time=aurostd::string2utype<double>(content);}
	else if(keyword=="corresponding") {corresponding=content;aurostd::string2tokens(content,stokens,",");for(uint j=0;j<stokens.size();j++) vcorresponding.push_back(stokens.at(j));}
	else if(keyword=="loop") {vloop.push_back(content);}  // CHECK THIS OUT IN THE FITURE
	else if(keyword=="node_CPU_Cores") {node_CPU_Cores=aurostd::string2utype<int>(content);}
	else if(keyword=="node_CPU_MHz") {node_CPU_MHz=aurostd::string2utype<double>(content);}
	else if(keyword=="node_CPU_Model") {node_CPU_Model=content;}
	else if(keyword=="node_RAM_GB") {node_RAM_GB=aurostd::string2utype<double>(content);}
	else if(keyword=="Bravais_lattice_orig") {Bravais_lattice_orig=content;}
	else if(keyword=="Bravais_lattice_relax") {Bravais_lattice_relax=content;}
	else if(keyword=="code") {code=content;}
	else if(keyword=="composition") {composition=content;for(uint j=0;j<stokens.size();j++) vcomposition.push_back(aurostd::string2utype<double>(stokens.at(j)));}
	else if(keyword=="compound") {compound=content;}
	else if(keyword=="density") {density=aurostd::string2utype<double>(content);}
	else if(keyword=="dft_type") {dft_type=content;aurostd::string2tokens(content,stokens,",");for(uint j=0;j<stokens.size();j++) vdft_type.push_back(stokens.at(j));}
	else if(keyword=="eentropy_cell") {eentropy_cell=aurostd::string2utype<double>(content);}
	else if(keyword=="eentropy_atom") {eentropy_atom=aurostd::string2utype<double>(content);}
	else if(keyword=="Egap") {Egap=aurostd::string2utype<double>(content);}
	else if(keyword=="Egap_fit") {Egap_fit=aurostd::string2utype<double>(content);}
	else if(keyword=="energy_cell") {energy_cell=aurostd::string2utype<double>(content);}
	else if(keyword=="energy_atom") {energy_atom=aurostd::string2utype<double>(content);energy_atom_relax1=aurostd::string2utype<double>(content);}
	else if(keyword=="energy_cutoff") {energy_cutoff=aurostd::string2utype<double>(content);}
	else if(keyword=="delta_electronic_energy_convergence") {delta_electronic_energy_convergence=aurostd::string2utype<double>(content);}
	else if(keyword=="delta_electronic_energy_threshold") {delta_electronic_energy_threshold=aurostd::string2utype<double>(content);}
	else if(keyword=="nkpoints") {nkpoints=aurostd::string2utype<uint>(content);}
	else if(keyword=="nkpoints_irreducible") {nkpoints_irreducible=aurostd::string2utype<uint>(content);}
	else if(keyword=="kppra") {kppra=aurostd::string2utype<uint>(content);}
	else if(keyword=="kpoints") {kpoints=content;}
	else if(keyword=="enthalpy_cell") {enthalpy_cell=aurostd::string2utype<double>(content);}
	else if(keyword=="enthalpy_atom") {enthalpy_atom=aurostd::string2utype<double>(content);}
	else if(keyword=="enthalpy_formation_cell") {enthalpy_formation_cell=aurostd::string2utype<double>(content);}
	else if(keyword=="enthalpy_formation_atom") {enthalpy_formation_atom=aurostd::string2utype<double>(content);}
	else if(keyword=="entropic_temperature") {entropic_temperature=aurostd::string2utype<double>(content);}
	else if(keyword=="files") {files=content;for(uint j=0;j<stokens.size();j++) vfiles.push_back(stokens.at(j));}
	else if(keyword=="files_LIB") {files_LIB=content;for(uint j=0;j<stokens.size();j++) vfiles_LIB.push_back(stokens.at(j));}
	else if(keyword=="files_RAW") {files_RAW=content;for(uint j=0;j<stokens.size();j++) vfiles_RAW.push_back(stokens.at(j));}
	else if(keyword=="files_WEB") {files_WEB=content;for(uint j=0;j<stokens.size();j++) vfiles_WEB.push_back(stokens.at(j));}
	else if(keyword=="forces") {forces=content;for(uint j=0;j<stokens.size();j++) vforces.push_back(aurostd::string2utype<double>(stokens.at(j)));}  // FIX
	else if(keyword=="geometry") {
	  geometry=content; 
	  vgeometry.push_back(0.0);vgeometry.push_back(0.0);vgeometry.push_back(0.0);
	  vgeometry.push_back(0.0);vgeometry.push_back(0.0);vgeometry.push_back(0.0);
	  if(stokens.size()==6) for(uint j=0;j<stokens.size();j++) vgeometry.at(j)=aurostd::string2utype<double>(stokens.at(j));
	}
	else if(keyword=="lattice_system_orig") {lattice_system_orig=content;}
	else if(keyword=="lattice_variation_orig") {lattice_variation_orig=content;}
	else if(keyword=="lattice_system_relax") {lattice_system_relax=content;}
	else if(keyword=="lattice_variation_relax") {lattice_variation_relax=content;}
	else if(keyword=="ldau_TLUJ") {ldau_TLUJ=content;}
	else if(keyword=="natoms") {natoms=aurostd::string2utype<int>(content);}
	else if(keyword=="nbondxx") {nbondxx=content;for(uint j=0;j<stokens.size();j++) vnbondxx.push_back(aurostd::string2utype<double>(stokens.at(j)));}
	else if(keyword=="nspecies") {nspecies=aurostd::string2utype<int>(content);}
	else if(keyword=="Pearson_symbol_orig") {Pearson_symbol_orig=content;}
	else if(keyword=="Pearson_symbol_relax") {Pearson_symbol_relax=content;}
	else if(keyword=="positions_cartesian") {positions_cartesian=content;for(uint j=0;j<stokens.size();j++) vpositions_cartesian.push_back(aurostd::string2utype<double>(stokens.at(j)));}  // FIX
	else if(keyword=="positions_fractional") {positions_fractional=content;for(uint j=0;j<stokens.size();j++) vpositions_fractional.push_back(aurostd::string2utype<double>(stokens.at(j)));}  // FIX
	else if(keyword=="pressure") {pressure=aurostd::string2utype<double>(content);}
	else if(keyword=="stress_tensor") {
	  stress_tensor=content; 
	  vstress_tensor.push_back(0.0);vstress_tensor.push_back(0.0);vstress_tensor.push_back(0.0);
	  vstress_tensor.push_back(0.0);vstress_tensor.push_back(0.0);vstress_tensor.push_back(0.0);
	  vstress_tensor.push_back(0.0);vstress_tensor.push_back(0.0);vstress_tensor.push_back(0.0);
	  if(stokens.size()==6) for(uint j=0;j<stokens.size();j++) vstress_tensor.at(j)=aurostd::string2utype<double>(stokens.at(j));
	}
	else if(keyword=="pressure_residual") {pressure_residual=aurostd::string2utype<double>(content);}
	else if(keyword=="Pulay_stress") {Pulay_stress=aurostd::string2utype<double>(content);}
	else if(keyword=="prototype") {prototype=content;}  // apennsy
	else if(keyword=="PV_cell") {PV_cell=aurostd::string2utype<double>(content);}
	else if(keyword=="PV_atom") {PV_atom=aurostd::string2utype<double>(content);}
	else if(keyword=="scintillation_attenuation_length") {scintillation_attenuation_length=aurostd::string2utype<double>(content);}
	else if(keyword=="sg") {sg=content;for(uint j=0;j<stokens.size();j++) vsg.push_back(stokens.at(j));} // CO 180101
	else if(keyword=="sg2") {sg2=content;for(uint j=0;j<stokens.size();j++) vsg2.push_back(stokens.at(j));} // CO 180101
	else if(keyword=="spacegroup_orig") {spacegroup_orig=content;}
	else if(keyword=="spacegroup_relax") {spacegroup_relax=content;}
	else if(keyword=="species") {species=content;for(uint j=0;j<stokens.size();j++) vspecies.push_back(stokens.at(j));}
	else if(keyword=="species_pp") {species_pp=content;for(uint j=0;j<stokens.size();j++) vspecies_pp.push_back(stokens.at(j));}
	else if(keyword=="species_pp_version") {species_pp_version=content;for(uint j=0;j<stokens.size();j++) vspecies_pp_version.push_back(stokens.at(j));}
	else if(keyword=="species_pp_ZVAL") {species_pp_ZVAL=content;for(uint j=0;j<stokens.size();j++) vspecies_pp_ZVAL.push_back(aurostd::string2utype<double>(stokens.at(j)));}
	else if(keyword=="spin_cell") {spin_cell=aurostd::string2utype<double>(content);}
	else if(keyword=="spin_atom") {spin_atom=aurostd::string2utype<double>(content);}
	else if(keyword=="spinD") {spinD=content;for(uint j=0;j<stokens.size();j++) vspinD.push_back(aurostd::string2utype<double>(stokens.at(j)));}
	else if(keyword=="spinD_magmom_orig") {spinD_magmom_orig=content;for(uint j=0;j<stokens.size();j++) vspinD_magmom_orig.push_back(aurostd::string2utype<double>(stokens.at(j)));}
	else if(keyword=="spinF") {spinF=aurostd::string2utype<double>(content);}
	else if(keyword=="sponsor") {sponsor=content;aurostd::string2tokens(content,stokens,",");for(uint j=0;j<stokens.size();j++) vsponsor.push_back(stokens.at(j));}
	else if(keyword=="stoichiometry") {stoichiometry=content;for(uint j=0;j<stokens.size();j++) vstoichiometry.push_back(aurostd::string2utype<double>(stokens.at(j)));}
	else if(keyword=="Egap_type") {Egap_type=content;}
	else if(keyword=="valence_cell_std") {valence_cell_std=aurostd::string2utype<double>(content);}
	else if(keyword=="valence_cell_iupac") {valence_cell_iupac=aurostd::string2utype<double>(content);}
	else if(keyword=="volume_cell") {volume_cell=aurostd::string2utype<double>(content);}
	else if(keyword=="volume_atom") {volume_atom=aurostd::string2utype<double>(content);}
 	// legacy
  else if(keyword=="server") {vserver.push_back(content);}
  else if(keyword=="stoich") {aurostd::string2tokens(ventry.at(i),tokens,"=");stoich=tokens.at(1);aurostd::string2tokens(stoich,stokens);for(uint j=0;j<stokens.size();j++) vstoich.push_back(aurostd::string2utype<double>(stokens.at(j)));}
        //DX 20180823 - added more symmetry info - START
        // SYMMETRY
        else if(keyword=="crystal_family") {crystal_family=content;}
        else if(keyword=="crystal_system") {crystal_system=content;}
        else if(keyword=="crystal_class") {crystal_class=content;}
        else if(keyword=="point_group_Hermann_Mauguin") {point_group_Hermann_Mauguin=content;}
        else if(keyword=="point_group_Schoenflies") {point_group_Schoenflies=content;}
        else if(keyword=="point_group_orbifold") {point_group_orbifold=content;}
        else if(keyword=="point_group_type") {point_group_type=content;}
        else if(keyword=="point_group_order") {point_group_order=aurostd::string2utype<uint>(content);}
        else if(keyword=="point_group_structure") {point_group_structure=content;}
        else if(keyword=="Bravais_lattice_lattice_type") {Bravais_lattice_lattice_type=content;}
        else if(keyword=="Bravais_lattice_lattice_variation_type") {Bravais_lattice_lattice_variation_type=content;}
        else if(keyword=="Bravais_lattice_lattice_system") {Bravais_lattice_lattice_system=content;}
        else if(keyword=="Bravais_superlattice_lattice_type") {Bravais_superlattice_lattice_type=content;}
        else if(keyword=="Bravais_superlattice_lattice_variation_type") {Bravais_superlattice_lattice_variation_type=content;}
        else if(keyword=="Bravais_superlattice_lattice_system") {Bravais_superlattice_lattice_system=content;}
        else if(keyword=="Pearson_symbol_superlattice") {Pearson_symbol_superlattice=content;}
	else if(keyword=="reciprocal_geometry") {
	  reciprocal_geometry=content; 
	  vreciprocal_geometry.push_back(0.0);vreciprocal_geometry.push_back(0.0);vreciprocal_geometry.push_back(0.0);
	  vreciprocal_geometry.push_back(0.0);vreciprocal_geometry.push_back(0.0);vreciprocal_geometry.push_back(0.0);
	  if(stokens.size()==6) for(uint j=0;j<stokens.size();j++) vreciprocal_geometry.at(j)=aurostd::string2utype<double>(stokens.at(j));
	}
	else if(keyword=="reciprocal_volume_cell") {reciprocal_volume_cell=aurostd::string2utype<double>(content);}
        else if(keyword=="reciprocal_lattice_type") {reciprocal_lattice_type=content;}
        else if(keyword=="reciprocal_lattice_variation_type") {reciprocal_lattice_variation_type=content;}
        else if(keyword=="Wyckoff_letters") {Wyckoff_letters=content;}
        else if(keyword=="Wyckoff_multiplicities") {Wyckoff_multiplicities=content;}
        else if(keyword=="Wyckoff_site_symmetries") {Wyckoff_site_symmetries=content;}
	//DX 20180823 - added more symmetry info - END
	// AGL/AEL
	else if(keyword=="agl_thermal_conductivity_300K") {agl_thermal_conductivity_300K=aurostd::string2utype<double>(content);}
	else if(keyword=="agl_debye") {agl_debye=aurostd::string2utype<double>(content);}
	else if(keyword=="agl_acoustic_debye") {agl_acoustic_debye=aurostd::string2utype<double>(content);}
	else if(keyword=="agl_gruneisen") {agl_gruneisen=aurostd::string2utype<double>(content);}
	else if(keyword=="agl_heat_capacity_Cv_300K") {agl_heat_capacity_Cv_300K=aurostd::string2utype<double>(content);}
	else if(keyword=="agl_heat_capacity_Cp_300K") {agl_heat_capacity_Cp_300K=aurostd::string2utype<double>(content);}
	else if(keyword=="agl_thermal_expansion_300K") {agl_thermal_expansion_300K=aurostd::string2utype<double>(content);}
	else if(keyword=="agl_bulk_modulus_static_300K") {agl_bulk_modulus_static_300K=aurostd::string2utype<double>(content);}
	else if(keyword=="agl_bulk_modulus_isothermal_300K") {agl_bulk_modulus_isothermal_300K=aurostd::string2utype<double>(content);}
	else if(keyword=="ael_poisson_ratio") {ael_poisson_ratio=aurostd::string2utype<double>(content);}
	else if(keyword=="ael_bulk_modulus_voigt") {ael_bulk_modulus_voigt=aurostd::string2utype<double>(content);}
	else if(keyword=="ael_bulk_modulus_reuss") {ael_bulk_modulus_reuss=aurostd::string2utype<double>(content);}
	else if(keyword=="ael_shear_modulus_voigt") {ael_shear_modulus_voigt=aurostd::string2utype<double>(content);}
	else if(keyword=="ael_shear_modulus_reuss") {ael_shear_modulus_reuss=aurostd::string2utype<double>(content);}
	else if(keyword=="ael_bulk_modulus_vrh") {ael_bulk_modulus_vrh=aurostd::string2utype<double>(content);}
	else if(keyword=="ael_shear_modulus_vrh") {ael_shear_modulus_vrh=aurostd::string2utype<double>(content);}
	else if(keyword=="ael_elastic_anistropy") {ael_elastic_anistropy=aurostd::string2utype<double>(content);}
	// BADER
  else if(keyword=="bader_net_charges") {bader_net_charges=content;aurostd::string2tokens<double>(content,vbader_net_charges,",");}
  else if(keyword=="bader_atomic_volumes") {bader_atomic_volumes=content;aurostd::string2tokens<double>(content,vbader_atomic_volumes,",");}
      }
    }
    // FIX LOOP
    loop="";
    vloop.push_back("aflowlib");
    for(uint j=0;j<vloop.size();j++) loop+=vloop.at(j)+(j<vloop.size()-1?", ":"");
    // FIX SERVER
    server="";
    for(uint j=0;j<vserver.size();j++) {
      server+=vserver.at(j)+(j<vserver.size()-1?", ":"");
      vserverdir.push_back(*(new vector<string>(0))); // space
    }
    // FIX ICSD
    if(aurostd::substring2bool(prototype,"_ICSD_")) {
      aurostd::string2tokens(prototype,tokens,"_");
      icsd=tokens.at(tokens.size()-1);
    }
    // FIX APENNSY
    structure_name=prototype;
    structure_description=prototype;
    distance_gnd=999999; // gotta calculate it
    distance_tie=999999; // gotta calculate it
    pureA=FALSE;pureB=FALSE;
    fcc=FALSE; bcc=FALSE;hcp=FALSE;
    stoich_a=999999;stoich_b=999999;
    bond_aa=999999;bond_ab=999999;bond_bb=999999;
    vNsgroup.clear();vsgroup.clear();vstr.clear();  // apennsy
    // DONE
    if(0) {							
      bool html=TRUE;
      oss << "Keywords" << endl;
      oss << "auid=" << auid << (html?"<br>":"") << endl;
      oss << "aurl=" << aurl << (html?"<br>":"") << endl;
      oss << "keywords=" << keywords << (html?"<br>":"") << "  vkeywords= ";for(uint j=0;j<vkeywords.size();j++) oss << vkeywords.at(j) << " "; oss << (html?"<br>":"") << endl;
      oss << "Optional controls keywords (alphabetic order)" << endl;
      oss << "aflowlib_date=" << aflowlib_date << (html?"<br>":"") << endl; 
      oss << "aflowlib_version=" << aflowlib_version << (html?"<br>":"") << endl; 
      oss << "aflowlib_entries=" << aflowlib_entries << (html?"<br>":"") << endl; 
      oss << "aflowlib_entries_number=" << aflowlib_entries_number << (html?"<br>":"") << endl; 
      oss << "aflow_version=" << aflow_version << (html?"<br>":"") << endl; 
      oss << "catalog=" << catalog << (html?"<br>":"") << endl; 
      oss << "data_api=" << data_api << (html?"<br>":"") << endl; 
      oss << "data_source=" << data_source << (html?"<br>":"") << endl; 
      oss << "data_language=" << data_language << (html?"<br>":"") << endl; 
      oss << "error_status=" << error_status << (html?"<br>":"") << endl; 
      oss << "author=" << author << (html?"<br>":"") << "  vauthor= ";for(uint j=0;j<vauthor.size();j++) oss << vauthor.at(j) << " "; oss << (html?"<br>":"") << endl;
      oss << "calculation_cores=" << calculation_cores << (html?"<br>":"") << endl; 
      oss << "calculation_memory=" << calculation_memory << (html?"<br>":"") << endl; 
      oss << "calculation_time=" << calculation_time << (html?"<br>":"") << endl; 
      oss << "corresponding=" << corresponding << (html?"<br>":"") << "  vcorresponding= ";for(uint j=0;j<vcorresponding.size();j++) oss << vcorresponding.at(j) << " "; oss << (html?"<br>":"") << endl;
      oss << "loop=" << loop << (html?"<br>":"") << "  vloop= ";for(uint j=0;j<vloop.size();j++) oss << vloop.at(j) << " "; oss << (html?"<br>":"") << endl;
      oss << "node_CPU_Cores=" << node_CPU_Cores << (html?"<br>":"") << endl; 
      oss << "node_CPU_MHz=" << node_CPU_MHz << (html?"<br>":"") << endl; 
      oss << "node_CPU_Model=" << node_CPU_Model << (html?"<br>":"") << endl; 
      oss << "node_RAM_GB=" << node_RAM_GB << (html?"<br>":"") << endl; 
      oss << "Optional materials keywords (alphabetic order)" << endl;
      oss << "Bravais_lattice_orig" << Bravais_lattice_orig << (html?"<br>":"") << endl;
      oss << "Bravais_lattice_relax" << Bravais_lattice_relax << (html?"<br>":"") << endl;
      oss << "code=" << code << (html?"<br>":"") << endl;
      oss << "composition=" << composition << "  vcomposition= ";for(uint j=0;j<vcomposition.size();j++) oss << vcomposition.at(j) << " "; oss << (html?"<br>":"") << endl;
      oss << "compound=" << compound << (html?"<br>":"") << endl;
      oss << "density=" << density << (html?"<br>":"") << endl;
      oss << "dft_type=" << dft_type << (html?"<br>":"") << "  vdft_type= ";for(uint j=0;j<vdft_type.size();j++) oss << vdft_type.at(j) << " "; oss << (html?"<br>":"") << endl;
      oss << "eentropy_cell=" << eentropy_cell << (html?"<br>":"") << endl; 
      oss << "eentropy_atom=" << eentropy_atom << (html?"<br>":"") << endl; 
      oss << "Egap=" << Egap << (html?"<br>":"") << endl; 
      oss << "Egap_fit=" << Egap_fit << (html?"<br>":"") << endl; 
      oss << "Egap_type=" << Egap_type << (html?"<br>":"") << endl;
      oss << "energy_cell=" << energy_cell << (html?"<br>":"") << endl; 
      oss << "energy_atom=" << energy_atom << (html?"<br>":"") << endl; 
      oss << "energy_cutoff=" << energy_cutoff << (html?"<br>":"") << endl; 
      oss << "delta_electronic_energy_convergence=" << delta_electronic_energy_convergence << (html?"<br>":"") << endl; 
      oss << "delta_electronic_energy_threshold=" << delta_electronic_energy_threshold << (html?"<br>":"") << endl; 
      oss << "nkpoints=" << nkpoints << (html?"<br>":"") << endl; 
      oss << "nkpoints_irreducible=" << nkpoints_irreducible << (html?"<br>":"") << endl; 
      oss << "kppra=" << kppra << (html?"<br>":"") << endl; 
      oss << "kpoints=" << kpoints << (html?"<br>":"") << endl;      
      oss << "enthalpy_cell=" << enthalpy_cell << (html?"<br>":"") << endl; 
      oss << "enthalpy_atom=" << enthalpy_atom << (html?"<br>":"") << endl; 
      oss << "enthalpy_formation_cell=" << enthalpy_formation_cell << (html?"<br>":"") << endl; 
      oss << "enthalpy_formation_atom=" << enthalpy_formation_atom << (html?"<br>":"") << endl; 
      oss << "entropic_temperature=" << entropic_temperature << (html?"<br>":"") << endl; 
      // oss << "files=" << files << "  vfiles= ";for(uint j=0;j<vfiles.size();j++) oss << vfiles.at(j) << " "; oss << (html?"<br>":"") << endl;
      // oss << "files_LIB=" << files_LIB << "  vfiles_LIB= ";for(uint j=0;j<vfiles_LIB.size();j++) oss << vfiles_LIB.at(j) << " "; oss << (html?"<br>":"") << endl;
      // oss << "files_RAW=" << files_RAW << "  vfiles_RAW= ";for(uint j=0;j<vfiles_RAW.size();j++) oss << vfiles_RAW.at(j) << " "; oss << (html?"<br>":"") << endl;
      // oss << "files_WEB=" << files_WEB << "  vfiles_WEB= ";for(uint j=0;j<vfiles_WEB.size();j++) oss << vfiles_WEB.at(j) << " "; oss << (html?"<br>":"") << endl;
      // oss << "forces=" << forces << "  vforces= ";for(uint j=0;j<vforces.size();j++) oss << vforces.at(j) << " "; oss << (html?"<br>":"") << endl;
      oss << "geometry=" << geometry << "  vgeometry= ";for(uint j=0;j<vgeometry.size();j++) oss << vgeometry.at(j) << " "; oss << (html?"<br>":"") << endl;
      oss << "lattice_system_orig" << lattice_system_orig << (html?"<br>":"") << endl;
      oss << "lattice_variation_orig" << lattice_variation_orig << (html?"<br>":"") << endl;
      oss << "lattice_system_relax" << lattice_system_relax << (html?"<br>":"") << endl;
      oss << "lattice_variation_relax" << lattice_variation_relax << (html?"<br>":"") << endl;
      oss << "ldau_TLUJ=" << ldau_TLUJ << (html?"<br>":"") << endl;      
      oss << "natoms=" << natoms << (html?"<br>":"") << endl;
      oss << "nbondxx=" << nbondxx << "  vnbondxx= ";for(uint j=0;j<vnbondxx.size();j++) oss << vnbondxx.at(j) << " "; oss << (html?"<br>":"") << endl;
      oss << "nspecies=" << nspecies << (html?"<br>":"") << endl;
      oss << "Pearson_symbol_orig" << Pearson_symbol_orig << (html?"<br>":"") << endl;
      oss << "Pearson_symbol_relax" << Pearson_symbol_relax << (html?"<br>":"") << endl;
      // oss << "positions_cartesian=" << positions_cartesian << "  vpositions_cartesian= ";for(uint j=0;j<vpositions_cartesian.size();j++) oss << vpositions_cartesian.at(j) << " "; oss << (html?"<br>":"") << endl;
      // oss << "positions_fractional=" << positions_fractional << "  vpositions_fractional= ";for(uint j=0;j<vpositions_fractional.size();j++) oss << vpositions_fractional.at(j) << " "; oss << (html?"<br>":"") << endl;
      oss << "pressure=" << pressure << (html?"<br>":"") << endl; 
      oss << "stress_tensor=" << stress_tensor << "  vstress_tensor= ";for(uint j=0;j<vstress_tensor.size();j++) oss << vstress_tensor.at(j) << " "; oss << (html?"<br>":"") << endl;
      oss << "pressure_residual=" << pressure_residual << (html?"<br>":"") << endl; 
      oss << "Pulay_stress=" << Pulay_stress << (html?"<br>":"") << endl; 
      oss << "prototype=" << prototype << (html?"<br>":"") << endl;
      oss << "PV_cell=" << PV_cell << (html?"<br>":"") << endl; 
      oss << "PV_atom=" << PV_atom << (html?"<br>":"") << endl; 
      oss << "scintillation_attenuation_length=" << scintillation_attenuation_length << (html?"<br>":"") << endl;
      oss << "sg=" << sg << (html?"<br>":"") << endl;
      oss << "sg2=" << sg2 << (html?"<br>":"") << endl;
      oss << "spacegroup_orig=" << spacegroup_orig << (html?"<br>":"") << endl;
      oss << "spacegroup_relax=" << spacegroup_relax << (html?"<br>":"") << endl;
      oss << "species=" << species << "  vspecies= ";for(uint j=0;j<vspecies.size();j++) oss << vspecies.at(j) << " "; oss << (html?"<br>":"") << endl;
      oss << "species_pp=" << species_pp << "  vspecies_pp= ";for(uint j=0;j<vspecies_pp.size();j++) oss << vspecies_pp.at(j) << " "; oss << (html?"<br>":"") << endl;
      oss << "species_pp_version=" << species_pp_version << "  vspecies_pp_version= ";for(uint j=0;j<vspecies_pp_version.size();j++) oss << vspecies_pp_version.at(j) << " "; oss << (html?"<br>":"") << endl;
      oss << "species_pp_ZVAL=" << species_pp_ZVAL << "  vspecies_pp_ZVAL= ";for(uint j=0;j<vspecies_pp_ZVAL.size();j++) oss << vspecies_pp_ZVAL.at(j) << " "; oss << (html?"<br>":"") << endl;
      oss << "spin_cell=" << spin_cell << (html?"<br>":"") << endl; 
      oss << "spin_atom=" << spin_atom << (html?"<br>":"") << endl; 
      oss << "spinD=" << spinD << "  vspinD= "; for(uint j=0;j<vspinD.size();j++) oss << vspinD.at(j) << " "; oss << (html?"<br>":"") << endl;
      oss << "spinD_magmom_orig=" << spinD_magmom_orig << "  vspinD_magmom_orig= "; for(uint j=0;j<vspinD_magmom_orig.size();j++) oss << vspinD_magmom_orig.at(j) << " "; oss << (html?"<br>":"") << endl;
      oss << "spinF=" << spinF << (html?"<br>":"") << endl;
      oss << "sponsor=" << sponsor << (html?"<br>":"") << "  vsponsor= ";for(uint j=0;j<vsponsor.size();j++) oss << vsponsor.at(j) << " "; oss << (html?"<br>":"") << endl;
      oss << "stoichiometry=" << stoichiometry << "  vstoichiometry= ";for(uint j=0;j<vstoichiometry.size();j++) oss << vstoichiometry.at(j) << " "; oss << (html?"<br>":"") << endl;
      oss << "valence_cell_std=" << valence_cell_std << (html?"<br>":"") << endl; 
      oss << "valence_cell_iupac=" << valence_cell_iupac << (html?"<br>":"") << endl;      
      oss << "volume_cell=" << volume_cell << (html?"<br>":"") << endl; 
      oss << "volume_atom=" << volume_atom << (html?"<br>":"") << endl; 
      //DX 20180823 - added more symmetry info - START
      // SYMMETRY
      oss << "crystal_family" << crystal_family << (html?"<br>":"") << endl;
      oss << "crystal_system" << crystal_system << (html?"<br>":"") << endl;
      oss << "crystal_class" << crystal_class << (html?"<br>":"") << endl;
      oss << "point_group_Hermann_Mauguin" << point_group_Hermann_Mauguin << (html?"<br>":"") << endl;
      oss << "point_group_Schoenflies" << point_group_Schoenflies << (html?"<br>":"") << endl;
      oss << "point_group_orbifold" << point_group_orbifold << (html?"<br>":"") << endl;
      oss << "point_group_type" << point_group_type << (html?"<br>":"") << endl;
      oss << "point_group_order" << point_group_order << (html?"<br>":"") << endl;
      oss << "point_group_structure" << point_group_structure << (html?"<br>":"") << endl;
      oss << "Bravais_lattice_lattice_type" << Bravais_lattice_lattice_type << (html?"<br>":"") << endl;
      oss << "Bravais_lattice_lattice_variation_type" << Bravais_lattice_lattice_variation_type << (html?"<br>":"") << endl;
      oss << "Bravais_lattice_lattice_system" << Bravais_lattice_lattice_system << (html?"<br>":"") << endl;
      oss << "Bravais_superlattice_lattice_type" << Bravais_superlattice_lattice_type << (html?"<br>":"") << endl;
      oss << "Bravais_superlattice_lattice_variation_type" << Bravais_superlattice_lattice_variation_type << (html?"<br>":"") << endl;
      oss << "Bravais_superlattice_lattice_system" << Bravais_superlattice_lattice_system << (html?"<br>":"") << endl;
      oss << "Pearson_symbol_superlattice" << Pearson_symbol_superlattice << (html?"<br>":"") << endl;
      oss << "reciprocal_geometry=" << reciprocal_geometry << "  vreciprocal_geometry= ";for(uint j=0;j<vreciprocal_geometry.size();j++) oss << vreciprocal_geometry.at(j) << " "; oss << (html?"<br>":"") << endl;
      oss << "reciprocal_volume_cell=" << reciprocal_volume_cell << (html?"<br>":"") << endl; 
      oss << "reciprocal_lattice_type" << reciprocal_lattice_type << (html?"<br>":"") << endl;
      oss << "reciprocal_lattice_variation_type" << reciprocal_lattice_variation_type << (html?"<br>":"") << endl;
      oss << "Wyckoff_letters" << Wyckoff_letters << (html?"<br>":"") << endl;
      oss << "Wyckoff_multiplicities" << Wyckoff_multiplicities << (html?"<br>":"") << endl;
      oss << "Wyckoff_site_symmetries" << Wyckoff_site_symmetries << (html?"<br>":"") << endl;
      //DX 20180823 - added more symmetry info - END
      // AGL/AEL
      oss << "agl_thermal_conductivity_300K" << agl_thermal_conductivity_300K << (html?"<br>":"") << endl; 
      oss << "agl_debye" << agl_debye << (html?"<br>":"") << endl; 
      oss << "agl_acoustic_debye" << agl_acoustic_debye << (html?"<br>":"") << endl; 
      oss << "agl_gruneisen" << agl_gruneisen << (html?"<br>":"") << endl; 
      oss << "agl_heat_capacity_Cv_300K" << agl_heat_capacity_Cv_300K << (html?"<br>":"") << endl; 
      oss << "agl_heat_capacity_Cp_300K" << agl_heat_capacity_Cp_300K << (html?"<br>":"") << endl; 
      oss << "agl_thermal_expansion_300K" << agl_thermal_expansion_300K << (html?"<br>":"") << endl; 
      oss << "agl_bulk_modulus_static_300K" << agl_bulk_modulus_static_300K << (html?"<br>":"") << endl; 
      oss << "agl_bulk_modulus_isothermal_300K" << agl_bulk_modulus_isothermal_300K << (html?"<br>":"") << endl; 
      oss << "ael_poisson_ratio" << ael_poisson_ratio << (html?"<br>":"") << endl; 
      oss << "ael_bulk_modulus_voigt" << ael_bulk_modulus_voigt << (html?"<br>":"") << endl; 
      oss << "ael_bulk_modulus_reuss" << ael_bulk_modulus_reuss << (html?"<br>":"") << endl; 
      oss << "ael_shear_modulus_voigt" << ael_shear_modulus_voigt << (html?"<br>":"") << endl; 
      oss << "ael_shear_modulus_reuss" << ael_shear_modulus_reuss << (html?"<br>":"") << endl; 
      oss << "ael_bulk_modulus_vrh" << ael_bulk_modulus_vrh << (html?"<br>":"") << endl; 
      oss << "ael_shear_modulus_vrh" << ael_shear_modulus_vrh << (html?"<br>":"") << endl; 
      oss << "ael_elastic_anistropy" << ael_elastic_anistropy << (html?"<br>":"") << endl; 
      // BADER
      oss << "bader_net_charges" << bader_net_charges << "  vbader_net_charges= ";for(uint j=0;j<vbader_net_charges.size();j++) oss << vbader_net_charges.at(j) << " "; oss << (html?"<br>":"") << endl; 
      oss << "bader_atomic_volumes" << bader_atomic_volumes << "  vbader_atomic_volumes= ";for(uint j=0;j<vbader_atomic_volumes.size();j++) oss << vbader_atomic_volumes.at(j) << " "; oss << (html?"<br>":"") << endl; 
      // legacy
      oss << "server=" << server << (html?"<br>":"") << "  vserver= ";for(uint j=0;j<vserver.size();j++) oss << vserver.at(j) << " "; oss << (html?"<br>":"") << endl;
      oss << "icsd=" << icsd << (html?"<br>":"") << endl;
      oss << "stoich=" << stoich << "  vstoich= ";for(uint j=0;j<vstoich.size();j++) oss << vstoich.at(j) << " "; oss << (html?"<br>":"") << endl;
    }
    return ventry.size();
  }

  // aflowlib2string 
  string _aflowlib_entry::aflowlib2string(string mode) {
    stringstream sss("");
    //  string eendl="\n";
    
    // this is the normal aflowlib.out mode
    if(mode=="" || mode=="out" || mode=="OUT") {
      string eendl="";
      
      if(auid.size()) sss << "" << "aurl=" << aurl << eendl;
      if(auid.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "auid=" << auid << eendl;
      if(data_api.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "data_api=" << data_api << eendl;
      if(data_source.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "data_source=" << data_source << eendl;
      if(data_language.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "data_language=" << data_language << eendl;
      if(error_status.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "error_status=" << error_status << eendl;
      // LOOP
      if(vloop.size()) {
	aurostd::sort(vloop);
	sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "loop=";
	for(uint i=0;i<vloop.size();i++) sss << vloop.at(i) << (i<vloop.size()-1?",":"");
	sss << eendl;
      }
      // MATERIALS
      if(code.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "code=" << code << eendl;
      if(compound.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "compound=" << compound << eendl;
      if(prototype.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "prototype=" << prototype << eendl;
      if(nspecies!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "nspecies=" << nspecies << eendl;
      if(natoms!=AUROSTD_NAN)sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "natoms=" << natoms << eendl;
      if(composition.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "composition=" << composition << eendl;
      if(density!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "density=" << density << eendl;
      if(scintillation_attenuation_length!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "scintillation_attenuation_length=" << scintillation_attenuation_length << eendl;
      if(stoichiometry.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "stoichiometry=" << stoichiometry << eendl;
      if(species.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "species=" << species << eendl;
      if(species_pp.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "species_pp=" << species_pp << eendl;
      if(dft_type.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "dft_type=" << dft_type << eendl;
      // if(species_pp_type.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "species_pp_type=" << species_pp_type << eendl;
      if(species_pp_version.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "species_pp_version=" << species_pp_version << eendl;
      if(species_pp_ZVAL.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "species_pp_ZVAL=" << species_pp_ZVAL << eendl;
      if(ldau_TLUJ.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "ldau_TLUJ=" << ldau_TLUJ << eendl;
      if(valence_cell_iupac!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "valence_cell_iupac=" << valence_cell_iupac << eendl;
      if(valence_cell_std!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "valence_cell_std=" << valence_cell_std << eendl;
      if(volume_cell!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "volume_cell=" << volume_cell << eendl;
      if(volume_atom!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "volume_atom=" << volume_atom << eendl;
      if(pressure!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "pressure=" << pressure << eendl;
      if(stress_tensor.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "stress_tensor=" << stress_tensor << eendl;
      if(pressure_residual!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "pressure_residual=" << pressure_residual << eendl;
      if(Pulay_stress!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "Pulay_stress=" << Pulay_stress << eendl;
      if(geometry.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "geometry=" << geometry << eendl;
      if(Egap!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "Egap=" << Egap << eendl;
      if(Egap_fit!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "Egap_fit=" << Egap_fit << eendl;
      if(Egap_type.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "Egap_type=" << Egap_type << eendl;
      if(energy_cell!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "energy_cell=" << energy_cell << eendl;
      if(energy_atom!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "energy_atom=" << energy_atom << eendl;
      if(energy_cutoff!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "energy_cutoff=" << energy_cutoff << eendl;
      if(delta_electronic_energy_convergence!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "delta_electronic_energy_convergence=" << delta_electronic_energy_convergence << eendl;
      if(delta_electronic_energy_threshold!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "delta_electronic_energy_threshold=" << delta_electronic_energy_threshold << eendl;
      // [NOT_PRINTED]     if(nkpoints!=0) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "nkpoints=" << nkpoints << eendl;
      // [NOT_PRINTED]     if(nkpoints_irreducible!=0) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "nkpoints_irreducible=" << nkpoints_irreducible << eendl;
      // [NOT_PRINTED]     if(kppra!=0) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "kppra=" << kppra << eendl;
      if(kpoints.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "kpoints=" << kpoints << eendl;
      if(enthalpy_cell!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "enthalpy_cell=" << enthalpy_cell << eendl;
      if(enthalpy_atom!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "enthalpy_atom=" << enthalpy_atom << eendl;
      if(eentropy_cell!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "eentropy_cell=" << eentropy_cell << eendl;
      if(eentropy_atom!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "eentropy_atom=" << eentropy_atom << eendl;
      if(enthalpy_formation_cell!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "enthalpy_formation_cell=" << enthalpy_formation_cell << eendl;
      if(enthalpy_formation_atom!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "enthalpy_formation_atom=" << enthalpy_formation_atom << eendl;
      if(entropic_temperature!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "entropic_temperature=" << entropic_temperature << eendl;
      if(PV_cell!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "PV_cell=" << PV_cell << eendl;
      if(PV_atom!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "PV_atom=" << PV_atom << eendl;
      if(spin_cell!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "spin_cell=" << spin_cell << eendl;
      if(spin_atom!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "spin_atom=" << spin_atom << eendl;
      if(spinD.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "spinD=" << spinD << eendl;
      if(spinD_magmom_orig.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "spinD_magmom_orig=" << spinD_magmom_orig << eendl;
      if(spinF!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "spinF=" << spinF << eendl;
      if(stoich.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "stoich=" << stoich << eendl;
      if(calculation_time!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "calculation_time=" << calculation_time << eendl;
      if(calculation_memory!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "calculation_memory=" << calculation_memory << eendl;
      if(calculation_cores!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "calculation_cores=" << calculation_cores << eendl;
      if(nbondxx.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "nbondxx=" << nbondxx << eendl;
      if(sg.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "sg=" << sg << eendl;
      if(sg2.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "sg2=" << sg2 << eendl;
      if(spacegroup_orig.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "spacegroup_orig=" << spacegroup_orig << eendl;
      if(spacegroup_relax.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "spacegroup_relax=" << spacegroup_relax << eendl;
      if(forces.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "forces=" << forces << eendl;
      if(positions_cartesian.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "positions_cartesian=" << positions_cartesian << eendl;
      if(positions_fractional.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "positions_fractional=" << positions_fractional << eendl;
      if(Bravais_lattice_orig.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "Bravais_lattice_orig=" << Bravais_lattice_orig << eendl;
      if(lattice_variation_orig.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "lattice_variation_orig=" << lattice_variation_orig << eendl;
      if(lattice_system_orig.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "lattice_system_orig=" << lattice_system_orig << eendl;
      if(Pearson_symbol_orig.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "Pearson_symbol_orig=" << Pearson_symbol_orig << eendl;
      if(Bravais_lattice_relax.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "Bravais_lattice_relax=" << Bravais_lattice_relax << eendl;
      if(lattice_variation_relax.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "lattice_variation_relax=" << lattice_variation_relax << eendl;
      if(lattice_system_relax.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "lattice_system_relax=" << lattice_system_relax << eendl;
      if(Pearson_symbol_relax.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "Pearson_symbol_relax=" << Pearson_symbol_relax << eendl;
      //DX 20180823 - added more symmetry info - START
      // SYMMETRY
      if(crystal_family.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "crystal_family=" << crystal_family << eendl;
      if(crystal_system.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "crystal_system=" << crystal_system << eendl;
      if(crystal_class.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "crystal_class=" << crystal_class << eendl;
      if(point_group_Hermann_Mauguin.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "point_group_Hermann_Mauguin=" << point_group_Hermann_Mauguin << eendl;
      if(point_group_Schoenflies.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "point_group_Schoenflies=" << point_group_Schoenflies << eendl;
      if(point_group_orbifold.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "point_group_orbifold=" << point_group_orbifold << eendl;
      if(point_group_type.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "point_group_type=" << point_group_type << eendl;
      if(point_group_order!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "point_group_order=" << point_group_order << eendl;
      if(point_group_structure.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "point_group_structure=" << point_group_structure << eendl;
      if(Bravais_lattice_lattice_type.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "Bravais_lattice_lattice_type=" << Bravais_lattice_lattice_type << eendl;
      if(Bravais_lattice_lattice_variation_type.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "Bravais_lattice_lattice_variation_type=" << Bravais_lattice_lattice_variation_type << eendl;
      if(Bravais_lattice_lattice_system.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "Bravais_lattice_lattice_system=" << Bravais_lattice_lattice_system << eendl;
      if(Bravais_superlattice_lattice_type.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "Bravais_superlattice_lattice_type=" << Bravais_superlattice_lattice_type << eendl;
      if(Bravais_superlattice_lattice_variation_type.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "Bravais_superlattice_lattice_variation_type=" << Bravais_superlattice_lattice_variation_type << eendl;
      if(Bravais_superlattice_lattice_system.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "Bravais_superlattice_lattice_system=" << Bravais_superlattice_lattice_system << eendl;
      if(Pearson_symbol_superlattice.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "Pearson_symbol_superlattice=" << Pearson_symbol_superlattice << eendl;
      if(reciprocal_geometry.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "reciprocal_geometry=" << reciprocal_geometry << eendl;
      if(reciprocal_volume_cell!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "reciprocal_volume_cell=" << reciprocal_volume_cell << eendl;
      if(reciprocal_lattice_type.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "reciprocal_lattice_type=" << reciprocal_lattice_type << eendl;
      if(reciprocal_lattice_variation_type.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "reciprocal_lattice_variation_type=" << reciprocal_lattice_variation_type << eendl;
      if(Wyckoff_letters.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "Wyckoff_letters=" << Wyckoff_letters << eendl;
      if(Wyckoff_multiplicities.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "Wyckoff_multiplicities=" << Wyckoff_multiplicities << eendl;
      if(Wyckoff_site_symmetries.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "Wyckoff_site_symmetries=" << Wyckoff_site_symmetries << eendl;
      //DX 20180823 - added more symmetry info - END
      // AGL/AEL
      if(agl_thermal_conductivity_300K!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "agl_thermal_conductivity_300K=" << agl_thermal_conductivity_300K << eendl;
      if(agl_debye!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "agl_debye=" << agl_debye << eendl;
      if(agl_acoustic_debye!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "agl_acoustic_debye=" << agl_acoustic_debye << eendl;
      if(agl_gruneisen!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "agl_gruneisen=" << agl_gruneisen << eendl;
      if(agl_heat_capacity_Cv_300K!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "agl_heat_capacity_Cv_300K=" << agl_heat_capacity_Cv_300K << eendl;
      if(agl_heat_capacity_Cp_300K!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "agl_heat_capacity_Cp_300K=" << agl_heat_capacity_Cp_300K << eendl;
      if(agl_thermal_expansion_300K!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "agl_thermal_expansion_300K=" << agl_thermal_expansion_300K << eendl;
      if(agl_bulk_modulus_static_300K!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "agl_bulk_modulus_static_300K=" << agl_bulk_modulus_static_300K << eendl;
      if(agl_bulk_modulus_isothermal_300K!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "agl_bulk_modulus_isothermal_300K=" << agl_bulk_modulus_isothermal_300K << eendl;
      if(ael_poisson_ratio!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "ael_poisson_ratio=" << ael_poisson_ratio << eendl;
      if(ael_bulk_modulus_voigt!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "ael_bulk_modulus_voigt=" << ael_bulk_modulus_voigt << eendl;
      if(ael_bulk_modulus_reuss!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "ael_bulk_modulus_reuss=" << ael_bulk_modulus_reuss << eendl;
      if(ael_shear_modulus_voigt!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "ael_shear_modulus_voigt=" << ael_shear_modulus_voigt << eendl;
      if(ael_shear_modulus_reuss!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "ael_shear_modulus_reuss=" << ael_shear_modulus_reuss << eendl;
      if(ael_bulk_modulus_vrh!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "ael_bulk_modulus_vrh=" << ael_bulk_modulus_vrh << eendl;
      if(ael_shear_modulus_vrh!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "ael_shear_modulus_vrh=" << ael_shear_modulus_vrh << eendl;
      if(ael_elastic_anistropy!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "ael_elastic_anistropy=" << ael_elastic_anistropy << eendl;
      // BADER
      if(bader_net_charges.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "bader_net_charges=" << bader_net_charges << eendl;
      if(bader_atomic_volumes.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "bader_atomic_volumes=" << bader_atomic_volumes << eendl;
      // FILES
      if(vfiles.size()) {
	aurostd::sort(vfiles);
	sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "files=";
	for(uint i=0;i<vfiles.size();i++) sss << vfiles.at(i) << (i<vfiles.size()-1?",":"");
	sss << eendl;
      }
      // CPUS
      if(node_CPU_Model.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "node_CPU_Model=" << node_CPU_Model << eendl;
      if(node_CPU_Cores!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "node_CPU_Cores=" << node_CPU_Cores << eendl;
      if(node_CPU_MHz!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "node_CPU_MHz=" << node_CPU_MHz << eendl;
      if(node_RAM_GB!=INF) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "node_RAM_GB=" << node_RAM_GB << eendl;
      // VERSION/DATE
      if(aflow_version.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "aflow_version=" << aflow_version << eendl;
      if(catalog.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "catalog=" << catalog << eendl;
      sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "aflowlib_version=" << string(AFLOW_VERSION) << eendl;
      sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "aflowlib_date=" << aurostd::get_datetime() << "_GMT-5" << eendl;
      sss << endl;

    } // out

    // this is the aflowlib.json mode
    if(mode=="json" || mode=="JSON") {  // COREY OPERATE HERE ALL THE STRINGS AS BEFORE
      string eendl="";
      bool PRINT_NULL=FALSE;
      stringstream sscontent_json;
      vector<string> vcontent_json;
      vector<string> sg_tokens;
      stringstream ss_helper;
      vector<vector<string> > vvs;
      vector<string> vs;
      bool odd_xvec_count;
      
      //////////////////////////////////////////////////////////////////////////
      if(auid.size()) {
        sscontent_json << "\"aurl\":\"" << aurl << "\"" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"aurl\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(auid.size()) {
        sscontent_json << "\"auid\":\"" << auid << "\"" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"auid\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////
      
      //////////////////////////////////////////////////////////////////////////
      if(data_api.size()) {
        sscontent_json << "\"data_api\":\"" << data_api << "\"" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"data_api\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(data_source.size()) {
        sscontent_json << "\"data_source\":\"" << data_source << "\"" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"data_source\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(data_language.size()) {
        sscontent_json << "\"data_language\":\"" << data_language << "\"" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"data_language\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(error_status.size()) {
        sscontent_json << "\"error_status\":\"" << error_status << "\"" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"error_status\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      // LOOP
      //////////////////////////////////////////////////////////////////////////
      if(vloop.size()) {
        aurostd::sort(vloop);
        sscontent_json << "\"loop\":[" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(vloop,"\""),",") << "]" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"loop\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      // MATERIALS
      //////////////////////////////////////////////////////////////////////////
      if(code.size()) {
        sscontent_json << "\"code\":\"" << code << "\"" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"code\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(compound.size()) {
        sscontent_json << "\"compound\":\"" << compound << "\"" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"compound\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(prototype.size()) {
        sscontent_json << "\"prototype\":\"" << prototype << "\"" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"prototype\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(nspecies!=AUROSTD_NAN) {
        sscontent_json << "\"nspecies\":" << nspecies << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"nspecies\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////
      
      //////////////////////////////////////////////////////////////////////////
      if(natoms!=AUROSTD_NAN) {
        sscontent_json << "\"natoms\":" << natoms << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"natoms\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////
      
      //////////////////////////////////////////////////////////////////////////
      if(vcomposition.size()) {
        //aflowlib_libraries does not specify precision
        sscontent_json << "\"composition\":[" << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(vcomposition,9),",") << "]" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"composition\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(density!=AUROSTD_NAN) {
        sscontent_json << "\"density\":" << density << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"density\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(scintillation_attenuation_length!=AUROSTD_NAN) {
        sscontent_json << "\"scintillation_attenuation_length\":" << scintillation_attenuation_length << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"scintillation_attenuation_length\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(vstoichiometry.size()) {
        //aflowlib_libraries specifies precision of 9
        sscontent_json << "\"stoichiometry\":[" << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(vstoichiometry,9),",") << "]" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"stoichiometry\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(vspecies.size()) {
        sscontent_json << "\"species\":[" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(vspecies,"\""),",") << "]" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"species\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(vspecies_pp.size()) {
        sscontent_json << "\"species_pp\":[" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(vspecies_pp,"\""),",") << "]" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"species_pp\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(dft_type.size()) {
        // DX and CO - START
        //sscontent_json << "\"dft_type\":\"" << dft_type << "\"" << eendl; // CO, this is technically a vector (RESTAPI paper)
        sscontent_json << "\"dft_type\":[" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(vdft_type,"\""),",") << "]" << eendl;
        // DX and CO - END
      } else {
        if(PRINT_NULL) sscontent_json << "\"dft_type\":null" << dft_type << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      // if(species_pp_type.size()) sscontent_json << "species_pp_type=" << species_pp_type << eendl;
      
      //////////////////////////////////////////////////////////////////////////
      if(vspecies_pp_version.size()) {
        sscontent_json << "\"species_pp_version\":[" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(vspecies_pp_version,"\""),",") << "]" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"species_pp_version\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////
      
      //////////////////////////////////////////////////////////////////////////
      if(vspecies_pp_ZVAL.size()) {
        //aflowlib_libraries does not specify precision
        sscontent_json << "\"species_pp_ZVAL\":[" << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(vspecies_pp_ZVAL,9),",") << "]" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"species_pp_ZVAL\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(ldau_TLUJ.size()) {
        ss_helper.str("");
        vs.clear();
        //only string is available, so we have to parse really fast
        //would be nice if we have vldau_TLUJ already
        vector<string> ldau_TLUJ_tokens;
        aurostd::string2tokens(ldau_TLUJ,ldau_TLUJ_tokens,";");
        if(ldau_TLUJ_tokens.size()==4){
          //conversion to double ENSURES that these are numbers
          //non-numbers without "" will break json
          int T=aurostd::string2utype<int>(ldau_TLUJ_tokens.at(0));
          vector<int> L;
          vector<double> U,J;
          vector<string> ldau_TLUJ_tokens2;
          //breaking up not necessary, but a nice check that we don't have hanging commas
          //the extra space at the end will be removed by joinWDelimiter()
          aurostd::string2tokens(ldau_TLUJ_tokens.at(1),L,",");
          aurostd::string2tokens(ldau_TLUJ_tokens.at(2),U,",");
          aurostd::string2tokens(ldau_TLUJ_tokens.at(3),J,",");
          if(L.size()&&U.size()&&J.size()){
            //no precision needed
            vs.push_back(aurostd::utype2string(T));
            vs.push_back("["+aurostd::joinWDelimiter(L,",")+"]");
            vs.push_back("["+aurostd::joinWDelimiter(aurostd::vecDouble2vecString(U,9),",")+"]");
            vs.push_back("["+aurostd::joinWDelimiter(aurostd::vecDouble2vecString(J,9),",")+"]");
            ss_helper << aurostd::joinWDelimiter(vs,",") << eendl;
            vs.clear();
          }
        }
        vs.clear();
      }
      if(!ss_helper.str().empty()){ // CO 180216 - !empty() is better for strings than !size()
        sscontent_json << "\"ldau_TLUJ\":[" << ss_helper.str() << "]" << eendl; ss_helper.str("");
      } else {
        if(PRINT_NULL) sscontent_json << "\"ldau_TLUJ\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(valence_cell_iupac!=AUROSTD_NAN) {
        sscontent_json << "\"valence_cell_iupac\":" << valence_cell_iupac << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"valence_cell_iupac\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(valence_cell_std!=AUROSTD_NAN) {
        sscontent_json << "\"valence_cell_std\":" << valence_cell_std << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"valence_cell_std\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(volume_cell!=AUROSTD_NAN) {
        sscontent_json << "\"volume_cell\":" << volume_cell << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"volume_cell\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(volume_atom!=AUROSTD_NAN) {
        sscontent_json << "\"volume_atom\":" << volume_atom << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"volume_atom\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(pressure!=AUROSTD_NAN) {
        sscontent_json << "\"pressure\":" << pressure << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"pressure\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////
      
      //////////////////////////////////////////////////////////////////////////
      if(vstress_tensor.size()) {
        //aflowlib_libraries specifies precision of 7
        sscontent_json << "\"stress_tensor\":[" << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(vstress_tensor,7),",") << "]" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"stress_tensor\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(pressure_residual!=AUROSTD_NAN) {
        sscontent_json << "\"pressure_residual\":" << pressure_residual << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"pressure_residual\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(Pulay_stress!=AUROSTD_NAN) {
        sscontent_json << "\"Pulay_stress\":" << Pulay_stress << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"Pulay_stress\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(vgeometry.size()) {
        //aflowlib_libraries specifies precision of 7
        sscontent_json << "\"geometry\":[" << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(vgeometry,7),",") << "]" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"geometry\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(Egap!=AUROSTD_NAN) {
        sscontent_json << "\"Egap\":" << Egap << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"Egap\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(Egap_fit!=AUROSTD_NAN) {
        sscontent_json << "\"Egap_fit\":" << Egap_fit << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"Egap_fit\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(Egap_type.size()) {
        sscontent_json << "\"Egap_type\":\"" << Egap_type << "\"" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"Egap_type\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(energy_cell!=AUROSTD_NAN) {
        sscontent_json << "\"energy_cell\":" << energy_cell << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"energy_cell\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(energy_atom!=AUROSTD_NAN) {
        sscontent_json << "\"energy_atom\":" << energy_atom << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"energy_atom\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(energy_cutoff!=AUROSTD_NAN) {
        sscontent_json << "\"energy_cutoff\":" << energy_cutoff << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"energy_cutoff\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(delta_electronic_energy_convergence!=AUROSTD_NAN) {
        sscontent_json << "\"delta_electronic_energy_convergence\":" << delta_electronic_energy_convergence << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"delta_electronic_energy_convergence\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(delta_electronic_energy_threshold!=AUROSTD_NAN) {
        sscontent_json << "\"delta_electronic_energy_threshold\":" << delta_electronic_energy_threshold << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"delta_electronic_energy_threshold\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      // [NOT_PRINTED]     //////////////////////////////////////////////////////////////////////////
      // [NOT_PRINTED]     if(nkpoints!=0) {
      // [NOT_PRINTED]       sscontent_json << "\"nkpoints\":" << nkpoints << eendl;
      // [NOT_PRINTED]     } else {
      // [NOT_PRINTED]       if(PRINT_NULL) sscontent_json << "\"nkpoints\":null" << eendl;
      // [NOT_PRINTED]     }
      // [NOT_PRINTED]     vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      // [NOT_PRINTED]     //////////////////////////////////////////////////////////////////////////

      // [NOT_PRINTED]     //////////////////////////////////////////////////////////////////////////
      // [NOT_PRINTED]     if(nkpoints_irreducible!=0) {
      // [NOT_PRINTED]       sscontent_json << "\"nkpoints_irreducible\":" << nkpoints_irreducible << eendl;
      // [NOT_PRINTED]     } else {
      // [NOT_PRINTED]       if(PRINT_NULL) sscontent_json << "\"nkpoints_irreducible\":null" << eendl;
      // [NOT_PRINTED]     }
      // [NOT_PRINTED]     vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      // [NOT_PRINTED]     //////////////////////////////////////////////////////////////////////////

      // [NOT_PRINTED]     //////////////////////////////////////////////////////////////////////////
      // [NOT_PRINTED]     if(kppra!=0) {
      // [NOT_PRINTED]       sscontent_json << "\"kppra\":" << kppra << eendl;
      // [NOT_PRINTED]     } else {
      // [NOT_PRINTED]       if(PRINT_NULL) sscontent_json << "\"kppra\":null" << eendl;
      // [NOT_PRINTED]     }
      // [NOT_PRINTED]     vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      // [NOT_PRINTED]     //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(kpoints.size()) {
        //this one is a bit complicated, so we will test if the string was created, and then recreate the json array on the spot
        ss_helper.str("");
        vs.clear();
        if(kpoints_nnn_relax.rows==3){
          vs.push_back("["+aurostd::joinWDelimiter(kpoints_nnn_relax,",")+"]");
        }
        if(kpoints_nnn_static.rows==3){
          vs.push_back("["+aurostd::joinWDelimiter(kpoints_nnn_static,",")+"]");
        }
        if(kpoints_pairs.size()){
          //first for escape characters in \Gamma or \Sigma
          vector<string> kpoints_pairs_new;
          char issue_c='\\';
          stringstream issue_ss; issue_ss << issue_c;
          string fix_s="\\\\";
          for(uint i=0;i<kpoints_pairs.size();i++){
            kpoints_pairs_new.push_back(aurostd::StringSubst(kpoints_pairs.at(i),issue_ss.str(),fix_s));
          }
          vs.push_back("["+aurostd::joinWDelimiter(aurostd::wrapVecEntries(kpoints_pairs_new,"\""),",")+"]");
        }
        if(kpoints_bands_path_grid!=0){
          ss_helper << aurostd::joinWDelimiter(vs,",") << "," << aurostd::utype2string(kpoints_bands_path_grid);
        }
      }
      if(!ss_helper.str().empty()){ // CO 180216 - !empty() is better for strings than !size()
        sscontent_json << "\"kpoints\":[" << ss_helper.str() << "]" << eendl; ss_helper.str("");
      } else {
        if(PRINT_NULL) sscontent_json << "\"kpoints\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(enthalpy_cell!=AUROSTD_NAN) {
        sscontent_json << "\"enthalpy_cell\":" << enthalpy_cell << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"enthalpy_cell\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(enthalpy_atom!=AUROSTD_NAN) {
        sscontent_json << "\"enthalpy_atom\":" << enthalpy_atom << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"enthalpy_atom\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(eentropy_cell!=AUROSTD_NAN) {
        sscontent_json << "\"eentropy_cell\":" << eentropy_cell << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"eentropy_cell\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(eentropy_atom!=AUROSTD_NAN) {
        sscontent_json << "\"eentropy_atom\":" << eentropy_atom << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"eentropy_atom\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(enthalpy_formation_cell!=AUROSTD_NAN) {
        sscontent_json << "\"enthalpy_formation_cell\":" << enthalpy_formation_cell << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"enthalpy_formation_cell\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(enthalpy_formation_atom!=AUROSTD_NAN) {
        sscontent_json << "\"enthalpy_formation_atom\":" << enthalpy_formation_atom << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"enthalpy_formation_atom\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(entropic_temperature!=AUROSTD_NAN) {
        sscontent_json << "\"entropic_temperature\":" << entropic_temperature << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"entropic_temperature\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(PV_cell!=AUROSTD_NAN) {
        sscontent_json << "\"PV_cell\":" << PV_cell << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"PV_cell\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(PV_atom!=AUROSTD_NAN) {
        sscontent_json << "\"PV_atom\":" << PV_atom << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"PV_atom\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(spin_cell!=AUROSTD_NAN) {
        sscontent_json << "\"spin_cell\":" << spin_cell << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"spin_cell\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(spin_atom!=AUROSTD_NAN) {
        sscontent_json << "\"spin_atom\":" << spin_atom << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"spin_atom\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(vspinD.size()) {
        //aflowlib_libraries specifies precision of 5
        sscontent_json << "\"spinD\":[" << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(vspinD,5),",") << "]" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"spinD\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(vspinD_magmom_orig.size()) {
        //aflowlib_libraries specifies precision of 5
        sscontent_json << "\"spinD_magmom_orig\":[" << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(vspinD_magmom_orig,5),",") << "]" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"spinD_magmom_orig\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(spinF!=AUROSTD_NAN) {
        sscontent_json << "\"spinF\":" << spinF << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"spinF\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      // [OBSOLETE] 
      // [OBSOLETE] if(stoich.size()) {
      // [OBSOLETE]   //just use the string stefano made
      // [OBSOLETE]   sscontent_json << "\"stoich\":\"" << stoich << "\"" << eendl;
      // [OBSOLETE]   ////aflowlib_libraries does not specify precision
      // [OBSOLETE]   //sscontent_json << "\"stoich\":[" << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(vstoich,9),",") << "]" << eendl;
      // [OBSOLETE] }else{
      // [OBSOLETE]   if(PRINT_NULL) sscontent_json << "\"stoich\":null" << eendl;
      // [OBSOLETE] }
      // [OBSOLETE] vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(calculation_time!=AUROSTD_NAN) {
        sscontent_json << "\"calculation_time\":" << calculation_time << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"calculation_time\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(calculation_memory!=AUROSTD_NAN) {
        sscontent_json << "\"calculation_memory\":" << calculation_memory << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"calculation_memory\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(calculation_cores!=AUROSTD_NAN) {
        sscontent_json << "\"calculation_cores\":" << calculation_cores << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"calculation_cores\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(vnbondxx.size()) {
        //aflowlib_libraries does not specify precision
        sscontent_json << "\"nbondxx\":[" << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(vnbondxx,9),",") << "]" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"nbondxx\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(sg.size()) {
        aurostd::string2tokens(sg,sg_tokens,",");
        sscontent_json << "\"sg\":[" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(sg_tokens,"\""),",") << "]" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"sg\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(sg2.size()) {
        aurostd::string2tokens(sg2,sg_tokens,",");
        sscontent_json << "\"sg2\":[" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(sg_tokens,"\""),",") << "]" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"sg2\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(spacegroup_orig.size()) {
        sscontent_json << "\"spacegroup_orig\":" << spacegroup_orig << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"spacegroup_orig\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(spacegroup_relax.size()) {
        sscontent_json << "\"spacegroup_relax\":" << spacegroup_relax << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"spacegroup_relax\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(vforces.size()) {
        ss_helper.str("");
        vs.clear();
        vvs.clear();
        odd_xvec_count=false;
        for(uint i=0;i<vforces.size();i++){
          if(vforces.at(i).rows==3){
            //aflowlib_libraries specifies precision of 8
            vvs.push_back(xvecDouble2vecString(vforces.at(i),8));
          } else {
            odd_xvec_count=true;
            break;
          }
        }
        if(!odd_xvec_count&&vvs.size()){
          for(uint i=0;i<vvs.size();i++){
            vs.push_back("["+aurostd::joinWDelimiter(vvs.at(i),",")+"]");
          }
          if(vs.size()){
            ss_helper << aurostd::joinWDelimiter(vs,",") << eendl;
          }
        }
        vs.clear();
        vvs.clear();
      }
      if(!ss_helper.str().empty()){ // CO 180216 - !empty() is better for strings than !size()
        sscontent_json << "\"forces\":[" << ss_helper.str() << "]" << eendl; ss_helper.str("");
      } else {
        if(PRINT_NULL) sscontent_json << "\"forces\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(vpositions_cartesian.size()) {
        ss_helper.str("");
        vs.clear();
        vvs.clear();
        odd_xvec_count=false;
        for(uint i=0;i<vpositions_cartesian.size();i++){
          if(vpositions_cartesian.at(i).rows==3){
            //aflowlib_libraries specifies precision of 8
            vvs.push_back(xvecDouble2vecString(vpositions_cartesian.at(i),8));
          } else {
            odd_xvec_count=true;
            break;
          }
        }
        if(!odd_xvec_count&&vvs.size()){
          for(uint i=0;i<vvs.size();i++){
            vs.push_back("["+aurostd::joinWDelimiter(vvs.at(i),",")+"]");
          }
          if(vs.size()){
            ss_helper << aurostd::joinWDelimiter(vs,",") << eendl;
          }
        }
        vs.clear();
        vvs.clear();
      }
      if(!ss_helper.str().empty()){ // CO 180216 - !empty() is better for strings than !size()
        sscontent_json << "\"positions_cartesian\":[" << ss_helper.str() << "]" << eendl; ss_helper.str("");
      } else {
        if(PRINT_NULL) sscontent_json << "\"positions_cartesian\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(vpositions_fractional.size()) {
        ss_helper.str("");
        vs.clear();
        vvs.clear();
        odd_xvec_count=false;
        for(uint i=0;i<vpositions_fractional.size();i++){
          if(vpositions_fractional.at(i).rows==3){
            //aflowlib_libraries specifies precision of 8
            vvs.push_back(xvecDouble2vecString(vpositions_fractional.at(i),8));
          } else {
            odd_xvec_count=true;
            break;
          }
        }
        if(!odd_xvec_count&&vvs.size()){
          for(uint i=0;i<vvs.size();i++){
            vs.push_back("["+aurostd::joinWDelimiter(vvs.at(i),",")+"]");
          }
          if(vs.size()){
            ss_helper << aurostd::joinWDelimiter(vs,",") << eendl;
          }
        }
        vs.clear();
        vvs.clear();
      }
      if(!ss_helper.str().empty()){ // CO 180216 - !empty() is better for strings than !size()
        sscontent_json << "\"positions_fractional\":[" << ss_helper.str() << "]" << eendl; ss_helper.str("");
      } else {
        if(PRINT_NULL) sscontent_json << "\"positions_fractional\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////
      
      //////////////////////////////////////////////////////////////////////////
      if(Bravais_lattice_orig.size()) {
        sscontent_json << "\"Bravais_lattice_orig\":\"" << Bravais_lattice_orig << "\"" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"Bravais_lattice_orig\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(lattice_variation_orig.size()) {
        sscontent_json << "\"lattice_variation_orig\":\"" << lattice_variation_orig << "\"" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"lattice_variation_orig\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(lattice_system_orig.size()) {
        sscontent_json << "\"lattice_system_orig\":\"" << lattice_system_orig << "\"" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"lattice_system_orig\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(Pearson_symbol_orig.size()) {
        sscontent_json << "\"Pearson_symbol_orig\":\"" << Pearson_symbol_orig << "\"" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"Pearson_symbol_orig\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(Bravais_lattice_relax.size()) {
        sscontent_json << "\"Bravais_lattice_relax\":\"" << Bravais_lattice_relax << "\"" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"Bravais_lattice_relax\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(lattice_variation_relax.size()) {
        sscontent_json << "\"lattice_variation_relax\":\"" << lattice_variation_relax << "\"" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"lattice_variation_relax\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(lattice_system_relax.size()) {
        sscontent_json << "\"lattice_system_relax\":\"" << lattice_system_relax << "\"" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"lattice_system_relax\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(Pearson_symbol_relax.size()) {
        sscontent_json << "\"Pearson_symbol_relax\":\"" << Pearson_symbol_relax << "\"" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"Pearson_symbol_relax\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //DX 20180823 - added more symmetry info - START
      // SYMMETRY
      //////////////////////////////////////////////////////////////////////////
      if(crystal_family.size()){
        sscontent_json << "\"crystal_family\":\"" << crystal_family << "\"" << eendl;
      } else{
        if(PRINT_NULL) sscontent_json << "\"crystal_family\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      
      //////////////////////////////////////////////////////////////////////////
      if(crystal_system.size()){
        sscontent_json << "\"crystal_system\":\"" << crystal_system << "\"" << eendl;
      } else{
        if(PRINT_NULL) sscontent_json << "\"crystal_system\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      
      //////////////////////////////////////////////////////////////////////////
      if(crystal_class.size()){
        sscontent_json << "\"crystal_class\":\"" << crystal_class << "\"" << eendl;
      } else{
        if(PRINT_NULL){ sscontent_json << "\"crystal_class\":null" << eendl;}
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
 
      //////////////////////////////////////////////////////////////////////////
      if(point_group_Hermann_Mauguin.size()){
        sscontent_json << "\"point_group_Hermann_Mauguin\":\"" << point_group_Hermann_Mauguin << "\"" << eendl;
      } else{
        if(PRINT_NULL) sscontent_json << "\"point_group_Hermann_Mauguin\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      
      //////////////////////////////////////////////////////////////////////////
      if(point_group_Schoenflies.size()){
        sscontent_json << "\"point_group_Schoenflies\":\"" << point_group_Schoenflies << "\"" << eendl;
      } else{
        if(PRINT_NULL) sscontent_json << "\"point_group_Schoenflies\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

      //////////////////////////////////////////////////////////////////////////
      if(point_group_orbifold.size()){
        sscontent_json << "\"point_group_orbifold\":\"" << point_group_orbifold << "\"" << eendl;
      } else{
        if(PRINT_NULL) sscontent_json << "\"point_group_orbifold\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      
      //////////////////////////////////////////////////////////////////////////
      if(point_group_type.size()){
        sscontent_json << "\"point_group_type\":\"" << point_group_type << "\"" << eendl;
      } else{
        if(PRINT_NULL) sscontent_json << "\"point_group_type\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

      //////////////////////////////////////////////////////////////////////////
      if(point_group_order!=AUROSTD_NAN){
        sscontent_json << "\"point_group_order\":\"" << point_group_order << "\"" << eendl;
      } else{
        if(PRINT_NULL) sscontent_json << "\"point_group_order\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

      //////////////////////////////////////////////////////////////////////////
      if(point_group_structure.size()){
        sscontent_json << "\"point_group_structure\":\"" << point_group_structure << "\"" << eendl;
      } else{
        if(PRINT_NULL) sscontent_json << "\"point_group_structure\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

      //////////////////////////////////////////////////////////////////////////
      if(Bravais_lattice_lattice_type.size()){
        sscontent_json << "\"Bravais_lattice_lattice_type\":\"" << Bravais_lattice_lattice_type << "\"" << eendl;
      } else{
        if(PRINT_NULL) sscontent_json << "\"Bravais_lattice_lattice_type\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

      //////////////////////////////////////////////////////////////////////////
      if(Bravais_lattice_lattice_variation_type.size()){
        sscontent_json << "\"Bravais_lattice_lattice_variation_type\":\"" << Bravais_lattice_lattice_variation_type << "\"" << eendl;
      } else{
        if(PRINT_NULL) sscontent_json << "\"Bravais_lattice_lattice_variation_type\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

      //////////////////////////////////////////////////////////////////////////
      if(Bravais_lattice_lattice_system.size()){
        sscontent_json << "\"Bravais_lattice_lattice_system\":\"" << Bravais_lattice_lattice_system << "\"" << eendl;
      } else{
        if(PRINT_NULL) sscontent_json << "\"Bravais_lattice_lattice_system\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

      //////////////////////////////////////////////////////////////////////////
      if(Bravais_superlattice_lattice_type.size()){
        sscontent_json << "\"Bravais_superlattice_lattice_type\":\"" << Bravais_superlattice_lattice_type << "\"" << eendl;
      } else{
        if(PRINT_NULL) sscontent_json << "\"Bravais_superlattice_lattice_type\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

      //////////////////////////////////////////////////////////////////////////
      if(Bravais_superlattice_lattice_variation_type.size()){
        sscontent_json << "\"Bravais_superlattice_lattice_variation_type\":\"" << Bravais_superlattice_lattice_variation_type << "\"" << eendl;
      } else{
        if(PRINT_NULL) sscontent_json << "\"Bravais_superlattice_lattice_variation_type\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

      //////////////////////////////////////////////////////////////////////////
      if(Bravais_superlattice_lattice_system.size()){
        sscontent_json << "\"Bravais_superlattice_lattice_system\":\"" << Bravais_superlattice_lattice_system << "\"" << eendl;
      } else{
        if(PRINT_NULL) sscontent_json << "\"Bravais_superlattice_lattice_system\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

      //////////////////////////////////////////////////////////////////////////
      if(Pearson_symbol_superlattice.size()){
        sscontent_json << "\"Pearson_symbol_superlattice\":\"" << Pearson_symbol_superlattice << "\"" << eendl;
      } else{
        if(PRINT_NULL) sscontent_json << "\"Pearson_symbol_superlattice\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      
      //////////////////////////////////////////////////////////////////////////
      if(vreciprocal_geometry.size()) {
        //aflowlib_libraries specifies precision of 7
        sscontent_json << "\"reciprocal_geometry\":[" << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(vreciprocal_geometry,7),",") << "]" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"reciprocal_geometry\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////
      
      //////////////////////////////////////////////////////////////////////////
      if(reciprocal_volume_cell!=AUROSTD_NAN) {
        sscontent_json << "\"reciprocal_volume_cell\":" << reciprocal_volume_cell << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"reciprocal_volume_cell\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////
      
      //////////////////////////////////////////////////////////////////////////
      if(reciprocal_lattice_type.size()){
        sscontent_json << "\"reciprocal_lattice_type\":\"" << reciprocal_lattice_type << "\"" << eendl;
      } else{
        if(PRINT_NULL) sscontent_json << "\"reciprocal_lattice_type\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      
      //////////////////////////////////////////////////////////////////////////
      if(reciprocal_lattice_variation_type.size()){
        sscontent_json << "\"reciprocal_lattice_variation_type\":\"" << reciprocal_lattice_variation_type << "\"" << eendl;
      } else{
        if(PRINT_NULL) sscontent_json << "\"reciprocal_lattice_variation_type\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      
      //////////////////////////////////////////////////////////////////////////
      if(Wyckoff_letters.size()){
        sscontent_json << "\"Wyckoff_letters\":\"" << Wyckoff_letters << "\"" << eendl;
      } else{
        if(PRINT_NULL) sscontent_json << "\"Wyckoff_letters\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

      //////////////////////////////////////////////////////////////////////////
      if(Wyckoff_multiplicities.size()){
        sscontent_json << "\"Wyckoff_multiplicities\":\"" << Wyckoff_multiplicities << "\"" << eendl;
      } else{
        if(PRINT_NULL) sscontent_json << "\"Wyckoff_multiplicities\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

      //////////////////////////////////////////////////////////////////////////
      if(Wyckoff_site_symmetries.size()){
        sscontent_json << "\"Wyckoff_site_symmetries\":\"" << Wyckoff_site_symmetries << "\"" << eendl;
      } else{
        if(PRINT_NULL) sscontent_json << "\"Wyckoff_site_symmetries\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");


      // AGL/AEL
      //////////////////////////////////////////////////////////////////////////
      if(agl_thermal_conductivity_300K!=AUROSTD_NAN) {
        sscontent_json << "\"agl_thermal_conductivity_300K\":" << agl_thermal_conductivity_300K << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"agl_thermal_conductivity_300K\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(agl_debye!=AUROSTD_NAN) {
        sscontent_json << "\"agl_debye\":" << agl_debye << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"agl_debye\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(agl_acoustic_debye!=AUROSTD_NAN) {
        sscontent_json << "\"agl_acoustic_debye\":" << agl_acoustic_debye << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"agl_acoustic_debye\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(agl_gruneisen!=AUROSTD_NAN) {
        sscontent_json << "\"agl_gruneisen\":" << agl_gruneisen << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"agl_gruneisen\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(agl_heat_capacity_Cv_300K!=AUROSTD_NAN) {
        sscontent_json << "\"agl_heat_capacity_Cv_300K\":" << agl_heat_capacity_Cv_300K << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"agl_heat_capacity_Cv_300K\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(agl_heat_capacity_Cp_300K!=AUROSTD_NAN) {
        sscontent_json << "\"agl_heat_capacity_Cp_300K\":" << agl_heat_capacity_Cp_300K << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"agl_heat_capacity_Cp_300K\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(agl_thermal_expansion_300K!=AUROSTD_NAN) {
        sscontent_json << "\"agl_thermal_expansion_300K\":" << agl_thermal_expansion_300K << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"agl_thermal_expansion_300K\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(agl_bulk_modulus_static_300K!=AUROSTD_NAN) {
        sscontent_json << "\"agl_bulk_modulus_static_300K\":" << agl_bulk_modulus_static_300K << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"agl_bulk_modulus_static_300K\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(agl_bulk_modulus_isothermal_300K!=AUROSTD_NAN) {
        sscontent_json << "\"agl_bulk_modulus_isothermal_300K\":" << agl_bulk_modulus_isothermal_300K << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"agl_bulk_modulus_isothermal_300K\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(ael_poisson_ratio!=AUROSTD_NAN) {
        sscontent_json << "\"ael_poisson_ratio\":" << ael_poisson_ratio << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"ael_poisson_ratio\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(ael_bulk_modulus_voigt!=AUROSTD_NAN) {
        sscontent_json << "\"ael_bulk_modulus_voigt\":" << ael_bulk_modulus_voigt << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"ael_bulk_modulus_voigt\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(ael_bulk_modulus_reuss!=AUROSTD_NAN) {
        sscontent_json << "\"ael_bulk_modulus_reuss\":" << ael_bulk_modulus_reuss << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"ael_bulk_modulus_reuss\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(ael_shear_modulus_voigt!=AUROSTD_NAN) {
        sscontent_json << "\"ael_shear_modulus_voigt\":" << ael_shear_modulus_voigt << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"ael_shear_modulus_voigt\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(ael_shear_modulus_reuss!=AUROSTD_NAN) {
        sscontent_json << "\"ael_shear_modulus_reuss\":" << ael_shear_modulus_reuss << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"ael_shear_modulus_reuss\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(ael_bulk_modulus_vrh!=AUROSTD_NAN) {
        sscontent_json << "\"ael_bulk_modulus_vrh=\":" << ael_bulk_modulus_vrh << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"ael_bulk_modulus_vrh=\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(ael_shear_modulus_vrh!=AUROSTD_NAN) {
        sscontent_json << "\"ael_shear_modulus_vrh\":" << ael_shear_modulus_vrh << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"ael_shear_modulus_vrh\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(ael_elastic_anistropy!=AUROSTD_NAN) {
        sscontent_json << "\"ael_elastic_anistropy\":" << ael_elastic_anistropy << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"ael_elastic_anistropy\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      // BADER
      //////////////////////////////////////////////////////////////////////////
      if(vbader_net_charges.size()) {
        sscontent_json << "\"bader_net_charges\":[" << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(vbader_net_charges,6),",") << "]" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"bader_net_charges\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(vbader_atomic_volumes.size()) {
        sscontent_json << "\"bader_atomic_volumes\":[" << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(vbader_atomic_volumes,4),",") << "]" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"bader_atomic_volumes\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      // FILES
      //////////////////////////////////////////////////////////////////////////
      if(vfiles.size()) {
        aurostd::sort(vfiles);
        sscontent_json << "\"files\":[" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(vfiles,"\""),",") << "]" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"files\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      // CPUS
      //////////////////////////////////////////////////////////////////////////
      if(node_CPU_Model.size()) {
        sscontent_json << "\"node_CPU_Model\":\"" << node_CPU_Model << "\"" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"node_CPU_Model\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(node_CPU_Cores!=AUROSTD_NAN) {
        sscontent_json << "\"node_CPU_Cores\":" << node_CPU_Cores << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"node_CPU_Cores\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(node_CPU_MHz!=AUROSTD_NAN) {
        sscontent_json << "\"node_CPU_MHz\":" << node_CPU_MHz << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"node_CPU_MHz\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(node_RAM_GB!=INF) {
        sscontent_json << "\"node_RAM_GB\":" << node_RAM_GB << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"node_RAM_GB\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      // VERSION/DATE
      //////////////////////////////////////////////////////////////////////////
      if(aflow_version.size()) {
        sscontent_json << "\"aflow_version\":\"" << aflow_version << "\"" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"aflow_version\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////
      
      //////////////////////////////////////////////////////////////////////////
      if(catalog.size()) {
        sscontent_json << "\"catalog\":\"" << catalog << "\"" << eendl;
      } else {
        if(PRINT_NULL) sscontent_json << "\"catalog\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      sscontent_json << "\"aflowlib_version\":\"" << string(AFLOW_VERSION) << "\"" << eendl;  // CO 170613
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      sscontent_json << "\"aflowlib_date\":\"" << aurostd::get_datetime() << "_GMT-5\"" << eendl;
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////
      
      sss << "{" << aurostd::joinWDelimiter(vcontent_json,",")  << "}";
      vcontent_json.clear();

      sss << endl;
    } // json

    return sss.str();
  }
  
  //  bool aflowlib2file(
  string _aflowlib_entry::aflowlib2file(string file,string mode) {
    string aflowlib_out=aflowlib2string(mode);
    aurostd::string2file(aflowlib_out,file);
    return aflowlib_out;
  }
}

// CO 171202 - apennsy fixes!
namespace aflowlib {
  void _aflowlib_entry::correctBadDatabase(bool verbose,ostream& oss){
    ofstream FileMESSAGE;
    return correctBadDatabase(FileMESSAGE,verbose,oss);
  }
  void _aflowlib_entry::correctBadDatabase(ofstream& FileMESSAGE,bool verbose,ostream& oss){
    //CO 180828 - LIB2 also contains unaries //so far we only know of bad binaries
    //APENNSY neglect - LIB2 only //CO 180828 - LIB2 also contains unaries  //binaries only
    string soliloquy="_aflowlib_entry::correctBadDatabase():";
    stringstream message;
    if(vspecies_pp.size()==1 || vspecies_pp.size()==2) {
      string pseudoA="",pseudoB="";
      pseudoA = vspecies_pp[0];
      if(vspecies_pp.size()==2){pseudoB = vspecies_pp[1];}
      //[OBSOLETE CO 180828]string pseudoA = vspecies_pp[0];
      //[OBSOLETE CO 180828]string pseudoB = vspecies_pp[1];
      // tiny corrections
      //gamma_IrV
      if(pseudoA == "Cd" && pseudoB == "Pt" && prototype == "181") {
        enthalpy_formation_atom -= 0.0013;
        enthalpy_formation_cell = natoms * enthalpy_formation_atom;
        if(verbose){
          message << "Fixing enthalpy_formation of " << pseudoA << pseudoB << ":" << prototype;
          pflow::logger(soliloquy, message, FileMESSAGE, oss, _LOGGER_MESSAGE_);
        }
      }
      //gamma_IrV
      if(pseudoA == "Ir" && pseudoB == "V_sv" && prototype == "291") {
        enthalpy_formation_cell -= 0.001;
        enthalpy_formation_atom -= 0.0005;
        enthalpy_cell -= 0.001;
        enthalpy_atom -= 0.005;
        if(verbose){
          message << "Fixing enthalpy/enthalpy_formation of " << pseudoA << pseudoB << ":" << prototype;
          pflow::logger(soliloquy, message, FileMESSAGE, oss, _LOGGER_MESSAGE_);
        }
      }
      // HfPd
      if(pseudoA == "Hf_pv" && pseudoB == "Pd_pv" && prototype == "192") {
        enthalpy_formation_atom -= 0.003;
        enthalpy_formation_cell = natoms * enthalpy_formation_atom;
        enthalpy_atom = enthalpy_formation_atom;
        enthalpy_cell = natoms * enthalpy_atom;
        if(verbose){
          message << "Fixing enthalpy/enthalpy_formation of " << pseudoA << pseudoB << ":" << prototype;
          pflow::logger(soliloquy, message, FileMESSAGE, oss, _LOGGER_MESSAGE_);
        }
      }
      // sigma
      if(pseudoA == "Ir" && pseudoB == "Nb_sv" && prototype == "600.ABBAB") {
        enthalpy_formation_cell += 0.001;
        enthalpy_formation_atom += 0.0005;
        enthalpy_cell += 0.001;
        enthalpy_atom += 0.005;
        enthalpy_formation_cell += 0.001;
        enthalpy_formation_atom += 0.0005;
        enthalpy_cell += 0.001;
        enthalpy_atom += 0.005;
        if(verbose){
          message << "Fixing enthalpy/enthalpy_formation of " << pseudoA << pseudoB << ":" << prototype;
          pflow::logger(soliloquy, message, FileMESSAGE, oss, _LOGGER_MESSAGE_);
        }
      }
      // sigma
      if(pseudoA == "Os_pv" && pseudoB == "Re_pv" && prototype == "122") {
        enthalpy_formation_atom -= 0.001;
        enthalpy_formation_cell = natoms * enthalpy_formation_atom;
        enthalpy_atom = enthalpy_formation_atom;
        enthalpy_cell = natoms * enthalpy_atom;
        if(verbose){
          message << "Fixing enthalpy/enthalpy_formation of " << pseudoA << pseudoB << ":" << prototype;
          pflow::logger(soliloquy, message, FileMESSAGE, oss, _LOGGER_MESSAGE_);
        }
      }
    }
  }
  bool _aflowlib_entry::ignoreBadDatabase() const{
    string reason;
    return ignoreBadDatabase(reason);
  }
  bool _aflowlib_entry::ignoreBadDatabase(string& reason) const{
    reason="";
    //grab bad database protos
    vector<string> _DEVIL_PROTOTYPES_;
    aurostd::string2tokens(_DEVIL_PROTOTYPES_STRING_,_DEVIL_PROTOTYPES_,",");
    
    //so far we only know of bad binaries
    //we need something more robust than just exact string match, case: 549 and 549.bis vs. 549.tetra
    bool match=false;
    //DEVIL
    for(uint di=0;di<_DEVIL_PROTOTYPES_.size() && !match;di++){if(pflow::prototypeMatch(prototype, _DEVIL_PROTOTYPES_[di])){match=true;}}
    if(match){
      reason=compound+":"+prototype+" is ill-calculated in the database";
      return true;
    }
    //find .old's
    if(1){
      uint prototype_size=prototype.size();
      string search_string=".old";uint search_string_size=search_string.size();
      if(prototype_size>search_string_size && prototype.substr(prototype_size-search_string_size,search_string_size)==search_string){  //look only at the end of the prototype
        reason=compound+":"+prototype+" is ill-calculated in the database";
        return true;
      }
    }
    //APENNSY neglect - LIB2 only //CO 180828 - LIB2 also contains unaries  //binaries only
    if(vspecies_pp.size()==1 || vspecies_pp.size()==2) {
      string pseudoA="",pseudoB="";
      pseudoA = vspecies_pp[0];
      if(vspecies_pp.size()==2){pseudoB = vspecies_pp[1];}
      //[OBSOLETE CO 180828]string pseudoA = vspecies_pp[0];
      //[OBSOLETE CO 180828]string pseudoB = vspecies_pp[1];
      // bad Ag is a wrong relaxation
      if((pseudoA == "Ag" && pflow::prototypeMatch(prototype, "303")) ||
          (pseudoB == "Ag" && pflow::prototypeMatch(prototype, "304"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Ag is a wrong relaxation
      if((pseudoA == "Ag" && pflow::prototypeMatch(prototype, "323")) ||
          (pseudoB == "Ag" && pflow::prototypeMatch(prototype, "324"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Au is a wrong relaxation
      if((pseudoA == "Au" && pflow::prototypeMatch(prototype, "323")) ||
          (pseudoB == "Au" && pflow::prototypeMatch(prototype, "324"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Al_h pseudopotential !
      if((pseudoA == "Al_h" && pflow::prototypeMatch(prototype, "307")) ||
          (pseudoB == "Al_h" && pflow::prototypeMatch(prototype, "308"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Al_h pseudopotential !
      if((pseudoA == "Al_h" && pflow::prototypeMatch(prototype, "A7.A")) ||
          (pseudoB == "Al_h" && pflow::prototypeMatch(prototype, "A7.B"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Al_h pseudopotential !
      if((pseudoA == "Al_h" && pflow::prototypeMatch(prototype, "323")) ||
          (pseudoB == "Al_h" && pflow::prototypeMatch(prototype, "324"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Ca_sv is a wrong relaxation
      if((pseudoA == "Ca_sv" && pflow::prototypeMatch(prototype, "303")) ||
          (pseudoB == "Ca_sv" && pflow::prototypeMatch(prototype, "304"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Ca_sv is a wrong relaxation
      if((pseudoA == "Ca_sv" && pflow::prototypeMatch(prototype, "323")) ||
          (pseudoB == "Ca_sv" && pflow::prototypeMatch(prototype, "324"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Cd is a wrong relaxation
      if((pseudoA == "Cd" && pflow::prototypeMatch(prototype, "323")) ||
          (pseudoB == "Cd" && pflow::prototypeMatch(prototype, "324"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Cu_pv is a wrong relaxation
      if((pseudoA == "Cu_pv" && pflow::prototypeMatch(prototype, "303")) ||
          (pseudoB == "Cu_pv" && pflow::prototypeMatch(prototype, "304"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Cu_pv is a wrong relaxation
      if((pseudoA == "Cu_pv" && pflow::prototypeMatch(prototype, "323")) ||
          (pseudoB == "Cu_pv" && pflow::prototypeMatch(prototype, "324"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Fe_pv is a wrong relaxation
      if((pseudoA == "Fe_pv" && pflow::prototypeMatch(prototype, "307")) ||
          (pseudoB == "Fe_pv" && pflow::prototypeMatch(prototype, "308"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Fe_pv is a wrong relaxation
      if((pseudoA == "Fe_pv" && pflow::prototypeMatch(prototype, "A7.A")) ||
          (pseudoB == "Fe_pv" && pflow::prototypeMatch(prototype, "A7.B"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Ge_h is a wrong relaxation
      if((pseudoA == "Ge_h" && pflow::prototypeMatch(prototype, "305")) ||
          (pseudoB == "Ge_h" && pflow::prototypeMatch(prototype, "306"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad In_d is a wrong relaxation
      if((pseudoA == "In_d" && pflow::prototypeMatch(prototype, "323")) ||
          (pseudoB == "In_d" && pflow::prototypeMatch(prototype, "324"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Ir is a wrong relaxation
      if((pseudoA == "Ir" && pflow::prototypeMatch(prototype, "303")) ||
          (pseudoB == "Ir" && pflow::prototypeMatch(prototype, "304"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad K_sv is a wrong relaxation
      if((pseudoA == "K_sv" && pflow::prototypeMatch(prototype, "307")) ||
          (pseudoB == "K_sv" && pflow::prototypeMatch(prototype, "308"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad K_sv is a wrong relaxation
      if((pseudoA == "K_sv" && pflow::prototypeMatch(prototype, "A7.A")) ||
          (pseudoB == "K_sv" && pflow::prototypeMatch(prototype, "A7.B"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad La is a wrong relaxation
      if((pseudoA == "La" && pflow::prototypeMatch(prototype, "303")) ||
          (pseudoB == "La" && pflow::prototypeMatch(prototype, "304"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad La is a wrong relaxation
      if((pseudoA == "La" && pflow::prototypeMatch(prototype, "323")) ||
          (pseudoB == "La" && pflow::prototypeMatch(prototype, "324"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Li_sv is a wrong relaxation
      if((pseudoA == "Li_sv" && pflow::prototypeMatch(prototype, "307")) ||
          (pseudoB == "Li_sv" && pflow::prototypeMatch(prototype, "308"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Li_sv is a wrong relaxation
      if((pseudoA == "Li_sv" && pflow::prototypeMatch(prototype, "A7.A")) ||
          (pseudoB == "Li_sv" && pflow::prototypeMatch(prototype, "A7.B"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Na_pv is a wrong relaxation
      if((pseudoA == "Na_pv" && pflow::prototypeMatch(prototype, "307")) ||
          (pseudoB == "Na_pv" && pflow::prototypeMatch(prototype, "308"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Na_pv is a wrong relaxation
      if((pseudoA == "Na_pv" && pflow::prototypeMatch(prototype, "A7.A")) ||
          (pseudoB == "Na_pv" && pflow::prototypeMatch(prototype, "A7.B"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Ni_pv is a wrong relaxation
      if((pseudoA == "Ni_pv" && pflow::prototypeMatch(prototype, "303")) ||
          (pseudoB == "Ni_pv" && pflow::prototypeMatch(prototype, "304"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Ni_pv is a wrong relaxation
      if((pseudoA == "Ni_pv" && pflow::prototypeMatch(prototype, "323")) ||
          (pseudoB == "Ni_pv" && pflow::prototypeMatch(prototype, "324"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Pb_d is a wrong relaxation
      if((pseudoA == "Pb_d" && pflow::prototypeMatch(prototype, "303")) ||
          (pseudoB == "Pb_d" && pflow::prototypeMatch(prototype, "304"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Pb_d is a wrong relaxation
      if((pseudoA == "Pb_d" && pflow::prototypeMatch(prototype, "323")) ||
          (pseudoB == "Pb_d" && pflow::prototypeMatch(prototype, "324"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Pd_pv is a wrong relaxation
      if((pseudoA == "Pd_pv" && pflow::prototypeMatch(prototype, "303")) ||
          (pseudoB == "Pd_pv" && pflow::prototypeMatch(prototype, "304"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Pd_pv is a wrong relaxation
      if((pseudoA == "Pd_pv" && pflow::prototypeMatch(prototype, "323")) ||
          (pseudoB == "Pd_pv" && pflow::prototypeMatch(prototype, "324"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Pt is a wrong relaxation
      if((pseudoA == "Pt" && pflow::prototypeMatch(prototype, "303")) ||
          (pseudoB == "Pt" && pflow::prototypeMatch(prototype, "304"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Pt is a wrong relaxation
      if((pseudoA == "Pt" && pflow::prototypeMatch(prototype, "317")) ||
          (pseudoB == "Pt" && pflow::prototypeMatch(prototype, "318"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }

      // bad Rh_pv is a wrong relaxation
      if((pseudoA == "Rh_pv" && pflow::prototypeMatch(prototype, "303")) ||
          (pseudoB == "Rh_pv" && pflow::prototypeMatch(prototype, "304"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Si_h is a wrong relaxation
      if((pseudoA == "Si_h" && pflow::prototypeMatch(prototype, "305")) ||
          (pseudoB == "Si_h" && pflow::prototypeMatch(prototype, "306"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Si_h is a wrong relaxation
      if((pseudoA == "Si_h" && pflow::prototypeMatch(prototype, "307")) ||
          (pseudoB == "Si_h" && pflow::prototypeMatch(prototype, "308"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Si_h is a wrong relaxation
      if((pseudoA == "Si_h" && pflow::prototypeMatch(prototype, "A7.A")) ||
          (pseudoB == "Si_h" && pflow::prototypeMatch(prototype, "A7.B"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Si_h is a wrong relaxation
      if((pseudoA == "Si_h" && pflow::prototypeMatch(prototype, "323")) ||
          (pseudoB == "Si_h" && pflow::prototypeMatch(prototype, "324"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Ta_pv is a wrong relaxation
      if((pseudoA == "Ta_pv" && pflow::prototypeMatch(prototype, "307")) ||
          (pseudoB == "Ta_pv" && pflow::prototypeMatch(prototype, "308"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Ta_pv is a wrong relaxation
      if((pseudoA == "Ta_pv" && pflow::prototypeMatch(prototype, "A7.A")) ||
          (pseudoB == "Ta_pv" && pflow::prototypeMatch(prototype, "A7.B"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad B_h is a wrong relaxation
      if((pseudoA == "B_h" && pflow::prototypeMatch(prototype, "317")) ||
          (pseudoB == "B_h" && pflow::prototypeMatch(prototype, "318"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }

      // sigma
      if(pseudoA == "Os_pv" && pseudoB == "Re_pv" &&
          pflow::prototypeMatch(prototype, "448")) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // wrong channel, bug
      if(pseudoA == "Rh_pv" && pseudoB == "Zr_sv" &&
          pflow::prototypeMatch(prototype, "381")) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
    }
    return false;
  }
} // namespace aflowlib

namespace aflowlib {
  string _aflowlib_entry::getPathAURL(ostream& oss, bool load_from_common){
    ofstream FileMESSAGE;
    return getPathAURL(FileMESSAGE, oss, load_from_common);
  }
  string _aflowlib_entry::getPathAURL(ofstream& FileMESSAGE,ostream& oss, bool load_from_common){
    string soliloquy = "_aflowlib_entry::getPathAURL():";
    stringstream message;
    string path = "";
    if (aurl.empty()) {return path;}
    vector<string> tokens;
    aurostd::string2tokens(aurl, tokens, ":");
    //LIB1 presents problems here (3 colons): aflowlib.duke.edu:AFLOWDATA/LIB1_RAW/Pt:PAW_PBE:05Jan2001/A6
    if(0){
    if (tokens.size() != 2) {
      message << "Odd AURL format for entry " << auid << ": " << aurl;
      pflow::logger(soliloquy, message, FileMESSAGE, oss, _LOGGER_WARNING_);
      return path;
    }
    }

    //instead just erase first item, join others, assume we're okay...
    tokens.erase(tokens.begin());
    path=aurostd::joinWDelimiter(tokens,":");

    if(load_from_common){return "/www/"+path;}//tokens.at(1);}
    else{
      string server;
      if (XHOST.vflag_control.flag("AFLOWLIB_SERVER")) {
        server = XHOST.vflag_control.getattachedscheme("AFLOWLIB_SERVER");
      } else {
        server = "aflowlib.duke.edu";
      }
      return server+"/"+path;//tokens.at(1);
    }
  }
}

// **************************************************************************
// directory2auid
// auid2present
// **************************************************************************
namespace aflowlib {
  string directory2auid(string directory,string aurl) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "directory2auid: BEGIN" << endl;
    string auid="";
    
    bool conflict=TRUE; 
    while (conflict) {
      uint64_t crc=0;
      // DONT TOUCH THE AUID FLOW
      deque<string> vext; aurostd::string2tokens(".bz2,.xz,.gz",vext,",");
      vector<string> vfiles2;
      for(uint iext=0;iext<vext.size();iext++) {
	vfiles2.push_back("OUTCAR.relax1"+vext.at(iext));
	vfiles2.push_back("OUTCAR.relax2"+vext.at(iext));
	vfiles2.push_back("OUTCAR.relax3"+vext.at(iext));
	vfiles2.push_back("OUTCAR.relax4"+vext.at(iext));
	vfiles2.push_back("OUTCAR.static"+vext.at(iext));
	vfiles2.push_back("OUTCAR.bands"+vext.at(iext));
	vfiles2.push_back("OUTCAR"+vext.at(iext));
      }
      for(uint i=0;i<vfiles2.size();i++) 
	if(aurostd::FileExist(directory+"/"+vfiles2.at(i)))
	  crc=aurostd::crc64(crc,aurostd::efile2string(directory+"/"+vfiles2.at(i))); // DONT TOUCH THIS
      auid="aflow:"+aurostd::crc2string(crc);
      if(LDEBUG) cerr << "directory2auid: auid=" << auid << endl;
      conflict=FALSE;
      uint j=auid2present(auid);
      if(j) {
	if(LDEBUG) cerr << "directory2auid: conflict auid=" << auid << endl;	
	cerr << "[WARNING]  directory2auid: CONFLICT POTENTIAL " << " j=" << j << " " << auid << " " << vAURL.at(j) << " " << aurl << endl;
	if(vAURL.at(j)!=aurl) { // avoid conflict with yourself
	  string salt="AUID_salt["+aurostd::utype2string<long double>(aurostd::get_useconds())+"]";
	  cerr << "[WARNING]  directory2auid: CONFLICT TRUE      " << " j=" << j << " " << auid << " " << vAURL.at(j) << " " << aurl << "  " << salt << endl;
	  string file=vfiles2.at(0);

	  // [OBSOLETE] aurostd::StringSubst(file,".EXT","");
	  // [OBSOLETE] stringstream sss;aurostd:EXTfile2stringstream(directory+"/"+file+".EXT",sss); sss << endl << salt << endl;
	  // [OBSOLETE] aurostd::execute("mv "+directory+"/"+file+".EXT"+" "+directory+"/"+file+".conflict_auid.EXT");
	  // [OBSOLETE] aurostd::stringstream2file(sss,directory+"/"+file);
	  // [OBSOLETE] aurostd::BzipFile(directory+"/"+file);

	  for(uint iext=0;iext<vext.size();iext++) {
	    aurostd::StringSubst(file,vext.at(iext),"");
	  }
	  stringstream sss;aurostd::efile2stringstream(directory+"/"+file+DEFAULT_KZIP_EXT,sss); sss << endl << salt << endl;
	  aurostd::execute("mv "+directory+"/"+file+DEFAULT_KZIP_EXT+" "+directory+"/"+file+".conflict_auid"+DEFAULT_KZIP_EXT);
	  aurostd::stringstream2compressfile(DEFAULT_KZIP_BIN,sss,directory+"/"+file);

	  conflict=TRUE; // recheck
	} else {
	  cerr << "[WARNING]  directory2auid: CONFLICT TRIVIAL   " << " j=" << j << " " << auid << " " << vAURL.at(j) << " " << aurl << endl;
	}
      }
    }
    if(LDEBUG) cerr << "directory2auid: END" << endl;
    return auid;
  }
  
  uint auid2present(string auid) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    uint k=0;
    if(LDEBUG) cerr << "auid2present: BEGIN" << endl;
    if(XHOST_AUID.size()==0) XHOST_AUID=init::InitGlobalObject("vAUID","",FALSE);
    if(XHOST_AURL.size()==0) XHOST_AURL=init::InitGlobalObject("vAURL","",FALSE);
    if(auid=="") return 0;
    if(LDEBUG) cerr << "auid2present: [4] vAURL.size()=" << vAURL.size() << endl;
    if(LDEBUG) cerr << "auid2present: [4] vAUID.size()=" << vAUID.size() << endl;
    bool found=FALSE;
    for(uint j=0;j<vAUID.size()&&!found;j++) {
      if(LDEBUG && vAUID.at(j)==auid) cerr << "[" << auid << "] [" << vAUID.at(j) << "]" << " [" << j << "]" << endl;
      if(vAUID.at(j)==auid) {k=j;found=FALSE;}
    }
    if(LDEBUG) cerr << "auid2present: END k=" << k << endl;
    return k;
  }
}

// **************************************************************************
// Operate on CIFs for JMOL for one structure
// **************************************************************************
namespace aflowlib {
  bool cif2data(string file,double& a,double& b,double& c,double& alpha,double& beta,double& gamma) {
    vector<string> vline,tokens;
    aurostd::file2vectorstring(file,vline);
    for(uint i=0;i<vline.size();i++) {
      aurostd::string2tokens(vline.at(i),tokens," ");
      if(aurostd::substring2bool(vline.at(i),"_cell_length_a") && tokens.size()>1) a=aurostd::string2utype<double>(tokens.at(1));
      if(aurostd::substring2bool(vline.at(i),"_cell_length_b") && tokens.size()>1) b=aurostd::string2utype<double>(tokens.at(1));
      if(aurostd::substring2bool(vline.at(i),"_cell_length_c") && tokens.size()>1) c=aurostd::string2utype<double>(tokens.at(1));
      if(aurostd::substring2bool(vline.at(i),"_cell_angle_alpha") && tokens.size()>1) alpha=aurostd::string2utype<double>(tokens.at(1));
      if(aurostd::substring2bool(vline.at(i),"_cell_angle_beta") && tokens.size()>1) beta=aurostd::string2utype<double>(tokens.at(1));
      if(aurostd::substring2bool(vline.at(i),"_cell_angle_gamma") && tokens.size()>1) gamma=aurostd::string2utype<double>(tokens.at(1));
    }
    return TRUE;
  }

  bool cif2oss(string file,string label,string sgnumber,ostream& oss) {
    double aCIF,bCIF,cCIF,alphaCIF,betaCIF,gammaCIF;
    aflowlib::cif2data(file,aCIF,bCIF,cCIF,alphaCIF,betaCIF,gammaCIF);
    oss << " font echo 14;color echo white; ";
    oss << " set echo 3% 97%; echo \\\"AFLOW-JSmol consortium (AFLOW V" << string(AFLOW_VERSION) << ") | entry="<< label << "  |  | ";
    if(aurostd::string2utype<int>(sgnumber)>0) {
      oss <<"  Spacegroup = "<< GetSpaceGroupName(aurostd::string2utype<int>(sgnumber)) << " (#" <<  sgnumber << ")   |";
    }
    oss << " a=" << aCIF <<"\u212B, b=" << bCIF <<"\u212B, c=" << cCIF <<"\u212B    |" << " \u03B1=" << alphaCIF <<"\u00B0, \u03B2=" << betaCIF <<"\u00B0, \u03B3=" << gammaCIF <<"\u00B0   \\\"; ";
    
    return TRUE;
  }
}

//BEGIN JJPR
// **************************************************************************
// GET BADER jvxl file
// **************************************************************************

namespace aflowlib {
  bool iso2oss(string file,string label, string element,string cutoff,int index,ostream& oss) {
    vector<string> kk;
    kk.resize(9);
    kk[0]="red";
    kk[1]="green";
    kk[2]="yellow";
    kk[3]="blue";
    kk[4]="orange";
    kk[5]="white";
    kk[6]="purple";
    kk[7]="brown";
    kk[8]="pink";
    oss <<"  ISOSURFACE  " << element << " \\\""<< file+"/"+label+"_Bader_"+cutoff+"_"+element+".jvxl\\\"  ;isosurface mesh;  color isosurface " << kk[index] << " translucent;";
    oss << "x = load(\\\""<< file+"/"+label+"_abader.out" << "\\\"); charges = x.split(\\\"=\\\")[2].split(\\\"(\\\")[1].split(\\\",\\\"); {*}.label = charges; label \%[label];";
    // FOR LOCAL TEST: oss <<"  ISOSURFACE  " << element << " "<< "Bader_"+cutoff+"_"+element+".jvxl  ;isosurface mesh;  color isosurface " << kk[index] << " translucent;";
    return TRUE;
  }
}
//END JJPR

// **************************************************************************
// GET SPACE GROUP for one structure
// **************************************************************************
namespace aflowlib {
  uint SGtoNSG(string sgroup) {
    string::size_type idx1;
    string strsub("#");
    idx1=sgroup.find(strsub);
    if(idx1!=string::npos)  
      return (int) atoi(sgroup.substr(sgroup.find(strsub)+strsub.length()).c_str());
    else return 0;
  }
  
  void _aflowlib_entry::GetSGROUP(string aflowlibentry) {
    vector<string> vaflowlib_entry;
    aurostd::string2tokens(aflowlibentry,vaflowlib_entry,"|");
    
    bool VERBOSE_LOCAL=(FALSE || XHOST.DEBUG);
    if(vaflowlib_entry.size()==0) {cerr << "ERROR - aflowlib_entry::GetSGROUP(): " << DEFAULT_FILE_AFLOWLIB_ENTRY_OUT << " file not found " << endl;exit(0);}
    vsgroup.clear();vNsgroup.clear();
    
    string data_aurl="";vector<string> tokens;
    // CHECK FOR HOLES
    if(XGNDSTATE_HOLES==0)  // SAFE NO HOLES IN THE XMATRIX
      if(vaflowlib_entry.size()==0) {cerr << "ERROR - aflowlib_entry::GetSGROUP(): " << DEFAULT_FILE_AFLOWLIB_ENTRY_OUT << " file not found " << endl;exit(0);}
    if(XGNDSTATE_HOLES==1)  // ALLOW HOLES WITH FAKE VALUES
      if(vaflowlib_entry.size()==0) {
	//cerr << "FOUND SG HOLE = " << alloy_dir << "/" << params.structures[2].at(str_number).name << "   === " << str_number << " " << structures_number[str_number]<< endl;
	vsgroup.clear();vsgroup.push_back(NOSG);vsgroup.push_back(NOSG);vsgroup.push_back(NOSG);
	vNsgroup.clear();vNsgroup.push_back(0);vNsgroup.push_back(0);vNsgroup.push_back(0);
	return;
      }
    
    for(uint i=0;i<vaflowlib_entry.size();i++) {
      if(aurostd::substring2bool(vaflowlib_entry.at(i),"aurl=")) data_aurl=aurostd::substring2string(vaflowlib_entry.at(i),"arl=");
      if(aurostd::substring2bool(vaflowlib_entry.at(i),"sg=")) {
	aurostd::string2tokens(vaflowlib_entry.at(i),tokens,",");
	if(tokens.size()==0) {cerr << "ERROR - aflowlib_entry::GetSGROUP(): geometry not enough tokens " << endl;exit(0);}
	if(tokens.size()==3) { // ok
	  vsgroup.clear();
	  aurostd::StringSubst(tokens.at(0),"sg=","");aurostd::StringSubst(tokens.at(0)," ","");aurostd::StringSubst(tokens.at(0),"#"," #");vsgroup.push_back(tokens.at(0));
	  aurostd::StringSubst(tokens.at(1)," ","");aurostd::StringSubst(tokens.at(1),"#"," #");vsgroup.push_back(tokens.at(1));
	  aurostd::StringSubst(tokens.at(2)," ","");aurostd::StringSubst(tokens.at(2),"#"," #");vsgroup.push_back(tokens.at(2));
	  if(VERBOSE_LOCAL) cout << "[5] "<<"["<<tokens.at(0) << "," << tokens.at(1) << "," << tokens.at(2) <<"]" << endl;
	} else {
	  vsgroup.clear();vsgroup.push_back(NOSG);vsgroup.push_back(NOSG);vsgroup.push_back(NOSG);
	  vNsgroup.clear();vNsgroup.push_back(0);vNsgroup.push_back(0);vNsgroup.push_back(0);
	}
      }
    }
    // done, now makes the numbers
    for(uint ii=0;ii<vsgroup.size();ii++) vNsgroup.push_back(SGtoNSG(vsgroup.at(ii)));
    // if(NsgroupPRE==0) {sgroupPRE=sgroupMID;NsgroupPRE=SGtoNSG(sgroupPRE);}	
    // if(NsgroupPRE==0) {sgroupPRE=sgroupPOST;NsgroupPRE=SGtoNSG(sgroupPRE);}	
    for(uint i=vsgroup.size()-2;i<vsgroup.size();i++) {
      if(vNsgroup.at(i)==0) {
	vsgroup.at(i)=vsgroup.at(i-1);
	vNsgroup.at(i)=SGtoNSG(vsgroup.at(i));
      }
    }
    
    if(vNsgroup.size()!=3) for(uint ii=vNsgroup.size();ii<3;ii++) vNsgroup.push_back(0);
    if(vsgroup.size()!=3) for(uint ii=vsgroup.size();ii<3;ii++) vsgroup.push_back("NNN #0");
    
    return;
  }
}


// ***************************************************************************
// aflowlib::AflowlibLocator
// ***************************************************************************
namespace aflowlib { // move to web interface
  bool AflowlibLocator(const string& in,string& out,const string& mode) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "aflowlib::AflowlibLocator: BEGIN" << endl;
    
    if(mode!="AFLOWLIB_AUID2AURL" && mode!="AFLOWLIB_AURL2AUID" && mode!="AFLOWLIB_AUID2LOOP" && mode!="AFLOWLIB_AURL2LOOP") {
      cerr << "ERROR - aflowlib::AflowlibLocator: wrong mode=" << mode << endl;
      exit(0);
    }
    if(XHOST_vAUID.size()==0) aurostd::string2vectorstring(init::InitGlobalObject("vAUID"),XHOST_vAUID); // cerr << XHOST_vAUID.size() << endl;
    if(XHOST_vAURL.size()==0) aurostd::string2vectorstring(init::InitGlobalObject("vAURL"),XHOST_vAURL); // cerr << XHOST_vAURL.size() << endl;
    if(XHOST_vLOOP.size()==0) aurostd::string2vectorstring(init::InitGlobalObject("vLOOP"),XHOST_vLOOP); // cerr << XHOST_vLOOP.size() << endl;
    if(XHOST_vAUID.size()!=XHOST_vAURL.size() || XHOST_vAUID.size()!=XHOST_vLOOP.size()) {
      cerr << "ERROR - aflowlib::AflowlibLocator: XHOST_vAUID.size()!=XHOST_vAURL.size() || XHOST_vAUID.size()!=XHOST_vLOOP.size()" << endl;
      cerr << "                                   XHOST_vAUID.size()=" << XHOST_vAUID.size() << endl;
      cerr << "                                   XHOST_vAURL.size()=" << XHOST_vAURL.size() << endl;
      cerr << "                                   XHOST_vLOOP.size()=" << XHOST_vLOOP.size() << endl;
      exit(0);
    }
    out="";
    for(uint i=0;i<XHOST_vAUID.size()&&out.empty();i++) {
      if(mode=="AFLOWLIB_AUID2AURL" && XHOST_vAUID.at(i)==in) out=XHOST_vAURL.at(i);
      if(mode=="AFLOWLIB_AURL2AUID" && XHOST_vAURL.at(i)==in) out=XHOST_vAUID.at(i);
      if(mode=="AFLOWLIB_AUID2LOOP" && XHOST_vAUID.at(i)==in) out=XHOST_vLOOP.at(i);
      if(mode=="AFLOWLIB_AURL2LOOP" && XHOST_vAURL.at(i)==in) out=XHOST_vLOOP.at(i);
      //     cerr << i << endl;
    }
    if(LDEBUG) cerr << "aflowlib::AflowlibLocator: END" << endl;
    return !out.empty();
  }
} // namespace aflowlib

namespace aflowlib {
  string AflowlibLocator(string options, string mode) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "aflowlib::AflowlibLocator: BEGIN" << endl;
    if(mode!="AFLOWLIB_AUID2AURL" && mode!="AFLOWLIB_AURL2AUID" && mode!="AFLOWLIB_AUID2LOOP" && mode!="AFLOWLIB_AURL2LOOP") {
      cerr << "error - aflowlib::AflowlibLocator: wrong mode=" << mode << endl;
      exit(0);
    }
    vector<string> tokens;
    aurostd::string2tokens(options,tokens,",");
    if(tokens.size()==0) {
      if(mode=="AFLOWLIB_AUID2AURL") init::ErrorOption(cout,options,"aflowlib::AflowlibLocator","aflow --aflowlib_auid2aurl=auid1,auid2....");
      if(mode=="AFLOWLIB_AURL2AUID") init::ErrorOption(cout,options,"aflowlib::AflowlibLocator","aflow --aflowlib_aurl2auid=aurl1,aurl2....");
      if(mode=="AFLOWLIB_AUID2LOOP") init::ErrorOption(cout,options,"aflowlib::AflowlibLocator","aflow --aflowlib_auid2loop=auid1,auid2....");
      if(mode=="AFLOWLIB_AURL2LOOP") init::ErrorOption(cout,options,"aflowlib::AflowlibLocator","aflow --aflowlib_aurl2loop=aurl1,aurl2....");
      exit(0);
    } 
    // move on
    stringstream output;
    string locator;
    for(uint i=0;i<tokens.size();i++) {
      if(aflowlib::AflowlibLocator(tokens.at(i),locator,mode)) {
	output << locator << endl;
      } else {
	output << tokens.at(i) << " not found" << endl;
      }
    }
    
    if(LDEBUG) cerr << "aflowlib::AflowlibLocator: END" << endl;
    return output.str();
  }
} // namespace aflowlib

//AFLUX integration
//FR & CO 180329
namespace aflowlib {
  bool APIget::establish(){
    struct hostent * host = gethostbyname( Domain.c_str() );

    PORT=80;  // CO 180401

    if ( (host == NULL) || (host->h_addr == NULL) ) {
        cerr << "Error retrieving DNS information." << endl;
        return false;
        //exit(1);
    }

    bzero(&client, sizeof(client));
    client.sin_family = AF_INET;
    client.sin_port = htons( PORT );
    memcpy(&client.sin_addr, host->h_addr, host->h_length);

    sock = socket(AF_INET, SOCK_STREAM, 0);

    if (sock < 0) {
        cerr << "Error creating socket." << endl;
        return false;
        //exit(1);
    }

    if ( connect(sock, (struct sockaddr *)&client, sizeof(client)) < 0 ) {
        close(sock);
        cerr << "Could not connect" << endl;
        return false;
        //exit(1);
    }

    stringstream ss;
    ss << "GET " << API_Path << Summons << " HTTP/1.0\r\n" ;
//    cerr << "GET " << API_Path << Summons << " HTTP/1.0\r\n" ;
    ss << "HOST: " << Domain << "\r\n";
    ss << "Connection: close\r\n";
    ss << "\r\n";
    string request = ss.str();

    if (send(sock, request.c_str(), request.length(), 0) != (int)request.length()) {
        cerr << "Error sending request." << endl;
        return false;
        //exit(1);
    }
    return true;
  }
  void APIget::reset( string a_Summons, string a_API_Path, string a_Domain ) {
    if( a_Summons == "#" ) {
        Summons = "";
        API_Path = "/search/API/?";
        Domain = "aflowlib.duke.edu";
    } else {
        Summons = a_Summons;
        if( ! a_API_Path.empty() ) API_Path = a_API_Path;
        if( ! a_Domain.empty() ) Domain = a_Domain;
    }
  }
  ostream& operator<<( ostream& output, APIget& a ) { 
    char cur;
    bool responsedata = false;
    bool waslinefeed = false;
    if( a.establish() ) {
        while ( ! responsedata ) { //discard headers
            read(a.sock, &cur, 1);
            //cerr << cur << ":" << (int)cur << endl;
            if( waslinefeed  && cur == '\r') responsedata = true;
            if( cur == '\n' ) waslinefeed = true;
            else waslinefeed = false;
        };
        read(a.sock, &cur, 1); //discard final \n in header \r\n\r\n
        while ( read(a.sock, &cur, 1) > 0 ) output << cur; //cout << cur;
        close(a.sock);
    }
    return output;
  }
}


// ***************************************************************************
namespace aflowlib {
  uint WEB_Aflowlib_Entry_PHP(string options,ostream& oss) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string soliloquy="aflowlib::WEB_Aflowlib_Entry_PHP():";
    if(LDEBUG) cout << "WEB_Aflowlib_Entry_PHP: begin<br>" << endl;
    
    int atomCOUNT=0;
    stringstream num_prec;
    vector<string> tokens;
    aurostd::string2tokens(options,tokens,",");
    if(tokens.size()==0) {
      init::ErrorOption(cout,options,"aflowlib::WEB_Aflowlib_Entry_PHP","aflow --aflowlib=entry");
      exit(0);
    } 
    
    // move on
    
    string option=tokens.at(0); //aurostd::args2attachedstring(argv,"--aflowlib=",(string) "nan");
    if(option.at(option.size()-1)=='/'|| option.at(option.size()-1)=='.') option.erase(option.end()-1,option.end()-0); //  some demoronization
    if(option.at(0)=='/'|| option.at(0)=='.') option.erase(option.begin(),option.begin()+1); //  some demoronization
    string directory="";
    string directory_RAW="";
    string directory_LIB="";
    string directory_WEB="";
    string url_WEB;
    string label="";
    //string line_gif="<br><img border=0 width=60% height=2 src=http://materials.duke.edu/auro/images/line.gif><br><br>";
    string line_rule="<hr width=\"60%\" style=\"background:black; border:0; height:2px; text-align:left;margin-left:0\" /><br>";
    string art058_link=" [<a href=https://doi.org/10.1016/j.commatsci.2010.05.010 target=\"_blank\"><font color=black><i>cite</i></font></a>]";
    string art064_link=" [<a href=https://doi.org/10.1021/co200012w target=\"_blank\"><font color=black><i>cite</i></font></a>]";
    string icsd_link=" [<a href=https://www.fiz-karlsruhe.com/icsd.html target=\"_blank\"><font color=black><i>info</i></font></a>]";
    string aflow_ael_readme=" [<a href=http://materials.duke.edu/AFLOW/README_AFLOW_AEL.TXT target=\"_blank\"><font color=black><i>info</i></font></a>]"; //CO 180817
    string art096_link=" [<a href=https://doi.org/10.1103/PhysRevB.90.174107 target=\"_blank\"><font color=black><i>cite</i></font></a>]";
    string art100_link=" [<a href=https://www.nature.com/articles/sdata20159 target=\"_blank\"><font color=black><i>cite</i></font></a>]";
    string aflow_agl_readme=" [<a href=http://materials.duke.edu/AFLOW/README_AFLOW_AGL.TXT target=\"_blank\"><font color=black><i>info</i></font></a>]"; //CO 180817
    string art115_link=" [<a href=https://doi.org/10.1103/PhysRevMaterials.1.015401 target=\"_blank\"><font color=black><i>cite</i></font></a>]"; //CO 180817
    string aflow_sym_readme=" [<a href=http://materials.duke.edu/AFLOW/README_AFLOW_SYM.TXT target=\"_blank\"><font color=black><i>info</i></font></a>]"; //CO 180817
    string art135_link=" [<a href=https://doi.org/10.1107/S2053273318003066 target=\"_blank\"><font color=black><i>cite</i></font></a>]"; //CO 180817

    //DX 180817
    string bravais_lattice_orig_wiki_link=" [<a href=http://aflowlib.duke.edu/aflowwiki/doku.php?id=documentation:all_keywords&#bravais_lattice_orig target=\"_blank\"><font color=black><i>info</i></font></a>]";
    string bravais_lattice_relax_wiki_link=" [<a href=http://aflowlib.duke.edu/aflowwiki/doku.php?id=documentation:all_keywords&#bravais_lattice_relax target=\"_blank\"><font color=black><i>info</i></font></a>]";
    string lattice_system_orig_wiki_link=" [<a href=http://aflowlib.duke.edu/aflowwiki/doku.php?id=documentation:all_keywords&#lattice_system_orig target=\"_blank\"><font color=black><i>info</i></font></a>]";
    string lattice_variation_orig_wiki_link=" [<a href=http://aflowlib.duke.edu/aflowwiki/doku.php?id=documentation:all_keywords&#lattice_variation_orig target=\"_blank\"><font color=black><i>info</i></font></a>]";
    string lattice_system_relax_wiki_link=" [<a href=http://aflowlib.duke.edu/aflowwiki/doku.php?id=documentation:all_keywords&#lattice_system_relax target=\"_blank\"><font color=black><i>info</i></font></a>]";
    string lattice_variation_relax_wiki_link=" [<a href=http://aflowlib.duke.edu/aflowwiki/doku.php?id=documentation:all_keywords&#lattice_variation_relax target=\"_blank\"><font color=black><i>info</i></font></a>]";
    string Pearson_symbol_orig_wiki_link=" [<a href=http://aflowlib.duke.edu/aflowwiki/doku.php?id=documentation:all_keywords&#pearson_symbol_orig target=\"_blank\"><font color=black><i>info</i></font></a>]";
    string Pearson_symbol_relax_wiki_link=" [<a href=http://aflowlib.duke.edu/aflowwiki/doku.php?id=documentation:all_keywords&#pearson_symbol_relax target=\"_blank\"><font color=black><i>info</i></font></a>]";
    string sg_wiki_link=" [<a href=http://aflowlib.duke.edu/aflowwiki/doku.php?id=documentation:all_keywords&#sg target=\"_blank\"><font color=black><i>info</i></font></a>]";
    string sg2_wiki_link=" [<a href=http://aflowlib.duke.edu/aflowwiki/doku.php?id=documentation:all_keywords&#sg2 target=\"_blank\"><font color=black><i>info</i></font></a>]";
    string spacegroup_orig_wiki_link=" [<a href=http://aflowlib.duke.edu/aflowwiki/doku.php?id=documentation:all_keywords&#spacegroup_orig target=\"_blank\"><font color=black><i>info</i></font></a>]";
    string spacegroup_relax_wiki_link=" [<a href=http://aflowlib.duke.edu/aflowwiki/doku.php?id=documentation:all_keywords&#spacegroup_relax target=\"_blank\"><font color=black><i>info</i></font></a>]";

    aflowlib::_aflowlib_entry aentry;
 
    xoption vflags;
    vflags.flag("FLAG::PREAMBLE",TRUE);
    vflags.flag("FLAG::CALCULATION",TRUE);
    vflags.flag("FLAG::JMOL",TRUE);
    vflags.flag("FLAG::EDATA_ORIG",FALSE);
    vflags.flag("FLAG::EDATA_RELAX",TRUE);
    vflags.flag("FLAG::THERMODYNAMICS",TRUE);
    vflags.flag("FLAG::MAGNETIC",TRUE);
    vflags.flag("FLAG::ELECTRONIC",FALSE);     // will setup later
    vflags.flag("FLAG::SCINTILLATION",TRUE);   // will setup later
    vflags.flag("FLAG::AGL",FALSE);            // will setup later
    vflags.flag("FLAG::AEL",FALSE);            // will setup later
    vflags.flag("FLAG::BADER",FALSE);          // will setup later
   
    // check if ICSD inside (anyway)
    string lattices[]={"BCC/","BCT/","CUB/","FCC/","HEX/","MCL/","MCLC/","ORC/","ORCC/","ORCF/","ORCI/","RHL/","TET/","TRI/"};
    vector<string> vline;
    // [OBSOLETE]   vector<string> tokens;
    
    string html_TAB=" target=\"_blank\"";

    if(LDEBUG) cout << "WEB_Aflowlib_Entry_PHP: [1]<br>" << endl;

    // check for ICSD
    if(!vflags.flag("FLAG::ICSD") && 
       !vflags.flag("FLAG::LIB1") && 
       !vflags.flag("FLAG::LIB2") && 
       !vflags.flag("FLAG::LIB3") && 
       !vflags.flag("FLAG::LIB4") && 
       !vflags.flag("FLAG::LIB5") && 
       !vflags.flag("FLAG::LIB6") && 
       !vflags.flag("FLAG::LIB7") && 
       !vflags.flag("FLAG::LIB8") && 
       !vflags.flag("FLAG::LIB9") && 
       !vflags.flag("FLAG::AUID") && 
       aurostd::substring2bool(option,"_ICSD_")) {
      vflags.flag("FLAG::ICSD",TRUE);
    }
    
    // check for LIB3
    if(!vflags.flag("FLAG::ICSD") && 
       !vflags.flag("FLAG::LIB1") && 
       !vflags.flag("FLAG::LIB2") && 
       !vflags.flag("FLAG::LIB3") && 
       !vflags.flag("FLAG::LIB4") && 
       !vflags.flag("FLAG::LIB5") && 
       !vflags.flag("FLAG::LIB6") && 
       !vflags.flag("FLAG::LIB7") && 
       !vflags.flag("FLAG::LIB8") && 
       !vflags.flag("FLAG::LIB9") && 
       !vflags.flag("FLAG::AUID") && 
       aurostd::substring2bool(option,"T000")) {
      vflags.flag("FLAG::LIB3",TRUE);
    }//aurostd::StringSubst(option,".T000","/T000");}
    
    // check for LIB4
    if(!vflags.flag("FLAG::ICSD") && 
       !vflags.flag("FLAG::LIB1") && 
       !vflags.flag("FLAG::LIB2") && 
       !vflags.flag("FLAG::LIB3") && 
       !vflags.flag("FLAG::LIB4") && 
       !vflags.flag("FLAG::LIB5") && 
       !vflags.flag("FLAG::LIB6") && 
       !vflags.flag("FLAG::LIB7") && 
       !vflags.flag("FLAG::LIB8") && 
       !vflags.flag("FLAG::LIB9") && 
       !vflags.flag("FLAG::AUID") && 
       aurostd::substring2bool(option,"Q000")) {
      vflags.flag("FLAG::LIB4",TRUE);
    }//aurostd::StringSubst(option,".Q000","/Q000");}
    
    // check for LIB5
    if(!vflags.flag("FLAG::ICSD") && 
       !vflags.flag("FLAG::LIB1") && 
       !vflags.flag("FLAG::LIB2") && 
       !vflags.flag("FLAG::LIB3") && 
       !vflags.flag("FLAG::LIB4") && 
       !vflags.flag("FLAG::LIB5") && 
       !vflags.flag("FLAG::LIB6") && 
       !vflags.flag("FLAG::LIB7") && 
       !vflags.flag("FLAG::LIB8") && 
       !vflags.flag("FLAG::LIB9") && 
       !vflags.flag("FLAG::AUID") && 
       aurostd::substring2bool(option,"P000")) {
      vflags.flag("FLAG::LIB5",TRUE);
    }//aurostd::StringSubst(option,".P000","/P000");}
    
    // check for LIB6
    if(!vflags.flag("FLAG::ICSD") && 
       !vflags.flag("FLAG::LIB1") && 
       !vflags.flag("FLAG::LIB2") && 
       !vflags.flag("FLAG::LIB3") && 
       !vflags.flag("FLAG::LIB4") && 
       !vflags.flag("FLAG::LIB5") && 
       !vflags.flag("FLAG::LIB6") && 
       !vflags.flag("FLAG::LIB7") && 
       !vflags.flag("FLAG::LIB8") && 
       !vflags.flag("FLAG::LIB9") && 
       !vflags.flag("FLAG::AUID") && 
       aurostd::substring2bool(option,"H000")) {
      vflags.flag("FLAG::LIB6",TRUE);
    }//aurostd::StringSubst(option,".H000","/H000");}
    
    // check for AUID
    if(!vflags.flag("FLAG::ICSD") && 
       !vflags.flag("FLAG::LIB1") && 
       !vflags.flag("FLAG::LIB2") && 
       !vflags.flag("FLAG::LIB3") && 
       !vflags.flag("FLAG::LIB4") && 
       !vflags.flag("FLAG::LIB5") && 
       !vflags.flag("FLAG::LIB6") && 
       !vflags.flag("FLAG::LIB7") && 
       !vflags.flag("FLAG::LIB8") && 
       !vflags.flag("FLAG::LIB9") && 
       !vflags.flag("FLAG::AUID") && 
       aurostd::substring2bool(option,"aflow:")) {
      vflags.flag("FLAG::AUID",TRUE);}
    
    // fix AUID
    
    // fix LIB2
    
    // find something
    if(!vflags.flag("FLAG::ICSD") && 
       !vflags.flag("FLAG::LIB1") && 
       !vflags.flag("FLAG::LIB2") && 
       !vflags.flag("FLAG::LIB3") && 
       !vflags.flag("FLAG::LIB4") && 
       !vflags.flag("FLAG::LIB5") && 
       !vflags.flag("FLAG::LIB6") && 
       !vflags.flag("FLAG::LIB7") && 
       !vflags.flag("FLAG::LIB8") && 
       !vflags.flag("FLAG::LIB9") && 
       !vflags.flag("FLAG::AUID") ) {
      string _option=option;
      aurostd::StringSubst(_option,"."," ");
      aurostd::StringSubst(_option,"/"," ");	
      aurostd::string2tokens(_option,tokens," ");
      vector<string> speciesX;
      XATOM_SplitAlloySpecies(tokens.at(0),speciesX);
      if(speciesX.size()==2) vflags.flag("FLAG::LIB2",TRUE);
      if(speciesX.size()==3) vflags.flag("FLAG::LIB3",TRUE);
      if(speciesX.size()==4) vflags.flag("FLAG::LIB4",TRUE);
      if(speciesX.size()==5) vflags.flag("FLAG::LIB5",TRUE);
      if(speciesX.size()==6) vflags.flag("FLAG::LIB6",TRUE);
      if(speciesX.size()==7) vflags.flag("FLAG::LIB7",TRUE);
      if(speciesX.size()==8) vflags.flag("FLAG::LIB8",TRUE);
      if(speciesX.size()==9) vflags.flag("FLAG::LIB9",TRUE);
    }

    // check/fix for LIB2
    if(!vflags.flag("FLAG::ICSD") && 
       !vflags.flag("FLAG::LIB1") && 
       !vflags.flag("FLAG::LIB2") && 
       !vflags.flag("FLAG::LIB3") && 
       !vflags.flag("FLAG::LIB4") && 
       !vflags.flag("FLAG::LIB5") && 
       !vflags.flag("FLAG::LIB6") && 
       !vflags.flag("FLAG::LIB7") && 
       !vflags.flag("FLAG::LIB8") && 
       !vflags.flag("FLAG::LIB9") && 
       !vflags.flag("FLAG::AUID") ) {
      aurostd::string2tokens(option,tokens,":");
      if(tokens.size()>0) {
	vector<string> speciesX;
	XATOM_SplitAlloySpecies(tokens.at(0),speciesX);
	if(speciesX.size()==1) vflags.flag("FLAG::LIB1",TRUE);
      }
    }
    
    if(LDEBUG) cout << "WEB_Aflowlib_Entry_PHP: [2]<br>" << endl;
    // nothing found try ICSD
    if(!vflags.flag("FLAG::ICSD") && 
       !vflags.flag("FLAG::LIB1") && 
       !vflags.flag("FLAG::LIB2") && 
       !vflags.flag("FLAG::LIB3") && 
       !vflags.flag("FLAG::LIB4") && 
       !vflags.flag("FLAG::LIB5") && 
       !vflags.flag("FLAG::LIB6") && 
       !vflags.flag("FLAG::LIB7") && 
       !vflags.flag("FLAG::LIB8") && 
       !vflags.flag("FLAG::LIB9") && 
       !vflags.flag("FLAG::AUID")) { // do something for AURO
      init::InitGlobalObject("Library_CALCULATED_ICSD_RAW");
      aurostd::string2vectorstring(XHOST_Library_CALCULATED_ICSD_RAW,vline);
      for(uint iline=0;iline<vline.size();iline++) {
	aurostd::StringSubst(vline.at(iline)," ","");
	aurostd::string2tokens(vline.at(iline),tokens,"_");
	if(tokens.size()>0)  {
	  if(tokens.at(tokens.size()-1)==option) {
	    directory=vline.at(iline);
	    vflags.flag("FLAG::ICSD",TRUE); 
	    aurostd::string2tokens(vline.at(iline),tokens,"/");
	    if(tokens.size()>0) option=tokens.at(tokens.size()-1);
	  }
	}
      }
    }	    
    
    if(LDEBUG) cout << "WEB_Aflowlib_Entry_PHP: [3]<br>" << endl;
    
    vflags.flag("FLAG::FOUND",(vflags.flag("FLAG::ICSD") || 
			       vflags.flag("FLAG::LIB1") || 
			       vflags.flag("FLAG::LIB2") || 
			       vflags.flag("FLAG::LIB3") || 
			       vflags.flag("FLAG::LIB4") || 
			       vflags.flag("FLAG::LIB5") || 
			       vflags.flag("FLAG::LIB6") || 
			       vflags.flag("FLAG::LIB7") || 
			       vflags.flag("FLAG::LIB8") || 
			       vflags.flag("FLAG::LIB9") || 
			       vflags.flag("FLAG::AUID")));
 
    if(vflags.flag("FLAG::FOUND") && vflags.flag("FLAG::ICSD")) init::InitGlobalObject("Library_CALCULATED_ICSD_RAW");
    if(vflags.flag("FLAG::FOUND") && vflags.flag("FLAG::LIB1")) init::InitGlobalObject("Library_CALCULATED_LIB1_RAW");
    if(vflags.flag("FLAG::FOUND") && vflags.flag("FLAG::LIB2")) init::InitGlobalObject("Library_CALCULATED_LIB2_RAW");
    if(vflags.flag("FLAG::FOUND") && vflags.flag("FLAG::LIB3")) init::InitGlobalObject("Library_CALCULATED_LIB3_RAW");
    if(vflags.flag("FLAG::FOUND") && vflags.flag("FLAG::LIB4")) init::InitGlobalObject("Library_CALCULATED_LIB4_RAW");
    if(vflags.flag("FLAG::FOUND") && vflags.flag("FLAG::LIB5")) init::InitGlobalObject("Library_CALCULATED_LIB5_RAW");
    if(vflags.flag("FLAG::FOUND") && vflags.flag("FLAG::LIB6")) init::InitGlobalObject("Library_CALCULATED_LIB6_RAW");
    if(vflags.flag("FLAG::FOUND") && vflags.flag("FLAG::LIB7")) init::InitGlobalObject("Library_CALCULATED_LIB7_RAW");
    if(vflags.flag("FLAG::FOUND") && vflags.flag("FLAG::LIB8")) init::InitGlobalObject("Library_CALCULATED_LIB8_RAW");
    if(vflags.flag("FLAG::FOUND") && vflags.flag("FLAG::LIB9")) init::InitGlobalObject("Library_CALCULATED_LIB9_RAW");
    //   if(vflags.flag("FLAG::FOUND") && vflags.flag("FLAG::AUID")) {init::InitGlobalObject("vAUID");init::InitGlobalObject("vAURL");}
    if(vflags.flag("FLAG::FOUND")) {init::InitGlobalObject("vAUID");init::InitGlobalObject("vAURL");}

    if(LDEBUG) cout << "WEB_Aflowlib_Entry_PHP: vAUID.size()=" << vAUID.size() << "<br>" << endl;
    if(LDEBUG) cout << "WEB_Aflowlib_Entry_PHP: vAURL.size()=" << vAURL.size() << "<br>" << endl;
    
    if(vflags.flag("FLAG::FOUND") && vflags.flag("FLAG::ICSD")) aurostd::string2vectorstring(XHOST_Library_CALCULATED_ICSD_RAW,vline);
    if(vflags.flag("FLAG::FOUND") && vflags.flag("FLAG::LIB1")) aurostd::string2vectorstring(XHOST_Library_CALCULATED_LIB1_RAW,vline);
    if(vflags.flag("FLAG::FOUND") && vflags.flag("FLAG::LIB3")) aurostd::string2vectorstring(XHOST_Library_CALCULATED_LIB3_RAW,vline);
    if(vflags.flag("FLAG::FOUND") && vflags.flag("FLAG::LIB2")) aurostd::string2vectorstring(XHOST_Library_CALCULATED_LIB2_RAW,vline);
    if(vflags.flag("FLAG::FOUND") && vflags.flag("FLAG::LIB4")) aurostd::string2vectorstring(XHOST_Library_CALCULATED_LIB4_RAW,vline);
    if(vflags.flag("FLAG::FOUND") && vflags.flag("FLAG::LIB5")) aurostd::string2vectorstring(XHOST_Library_CALCULATED_LIB5_RAW,vline);
    if(vflags.flag("FLAG::FOUND") && vflags.flag("FLAG::LIB6")) aurostd::string2vectorstring(XHOST_Library_CALCULATED_LIB6_RAW,vline);
    if(vflags.flag("FLAG::FOUND") && vflags.flag("FLAG::LIB7")) aurostd::string2vectorstring(XHOST_Library_CALCULATED_LIB7_RAW,vline);
    if(vflags.flag("FLAG::FOUND") && vflags.flag("FLAG::LIB8")) aurostd::string2vectorstring(XHOST_Library_CALCULATED_LIB8_RAW,vline);
    if(vflags.flag("FLAG::FOUND") && vflags.flag("FLAG::LIB9")) aurostd::string2vectorstring(XHOST_Library_CALCULATED_LIB9_RAW,vline);
   
    if(LDEBUG) cout << "WEB_Aflowlib_Entry_PHP: [4]<br>" << endl;

    if(vflags.flag("FLAG::FOUND")) {
      if(vflags.flag("FLAG::AUID")) {
	for(uint i=0;i<vAUID.size()&&i<vAURL.size();i++)
	  if(vAUID.at(i)==option) directory=vAURL.at(i);
	aurostd::StringSubst(directory,"aflowlib.duke.edu:","");
	aurostd::StringSubst(directory,"materials.duke.edu:","");
      }
      // for all the others
      for(uint iline=0;iline<vline.size()&&directory.empty();iline++) {
	aurostd::StringSubst(vline.at(iline)," ","");
	if(vflags.flag("FLAG::ICSD")&&directory.empty()) {
	  if(aurostd::substring2bool(vline.at(iline),option)) {
	    for(uint ilat=0;ilat<14&&directory.empty();ilat++) 
	      if(vline.at(iline)==lattices[ilat]+option)  
		directory=vline.at(iline);
	  }
	}
	if(vflags.flag("FLAG::LIB3")) {
	  string _option=option;
	  aurostd::StringSubst(_option,".T000","/T000");
	  aurostd::StringSubst( option,"/T000",".T000");	
	  if(aurostd::substring2bool(vline.at(iline),_option)) {
	    if(vline.at(iline)==_option)
	      directory=vline.at(iline);
	  }
	}
	if(vflags.flag("FLAG::LIB4")) {
	  string _option=option;
	  aurostd::StringSubst(_option,".Q000","/Q000");
	  aurostd::StringSubst( option,"/Q000",".Q000");	
	  if(aurostd::substring2bool(vline.at(iline),_option)) {
	    if(vline.at(iline)==_option)
	      directory=vline.at(iline);
	  }
	}
	if(vflags.flag("FLAG::LIB5")) {
	  string _option=option;
	  aurostd::StringSubst(_option,".P000","/P000");
	  aurostd::StringSubst( option,"/P000",".P000");	
	  if(aurostd::substring2bool(vline.at(iline),_option)) {
	    if(vline.at(iline)==_option)
	      directory=vline.at(iline);
	  }
	}
	if(vflags.flag("FLAG::LIB6")) {
	  string _option=option;
	  aurostd::StringSubst(_option,".H000","/H000");
	  aurostd::StringSubst( option,"/H000",".H000");	
	  if(aurostd::substring2bool(vline.at(iline),_option)) {
	    if(vline.at(iline)==_option)
	      directory=vline.at(iline);
	  }
	}
	// if(vflags.flag("FLAG::LIB7")) {
	//   string _option=option;
	//   aurostd::StringSubst(_option,".H000","/H000");
	//   aurostd::StringSubst( option,"/H000",".H000");	
	//   if(aurostd::substring2bool(vline.at(iline),_option)) {
	//     if(vline.at(iline)==_option)
	//       directory=vline.at(iline);
	//   }
	// }
	// if(vflags.flag("FLAG::LIB8")) {
	//   string _option=option;
	//   aurostd::StringSubst(_option,".H000","/H000");
	//   aurostd::StringSubst( option,"/H000",".H000");	
	//   if(aurostd::substring2bool(vline.at(iline),_option)) {
	//     if(vline.at(iline)==_option)
	//       directory=vline.at(iline);
	//   }
	// }
	// if(vflags.flag("FLAG::LIB9")) {
	//   string _option=option;
	//   aurostd::StringSubst(_option,".H000","/H000");
	//   aurostd::StringSubst( option,"/H000",".H000");	
	//   if(aurostd::substring2bool(vline.at(iline),_option)) {
	//     if(vline.at(iline)==_option)
	//       directory=vline.at(iline);
	//   }
	// }
	if(vflags.flag("FLAG::LIB2") || vflags.flag("FLAG::LIB1")) {
	  string _option=option;
	  aurostd::StringSubst(_option,"."," ");
	  aurostd::StringSubst(_option,"/"," ");	
	  aurostd::string2tokens(_option,tokens," ");
	  directory=tokens.at(0)+"/";
	  for(uint i=1;i<tokens.size();i++)
	    directory+=tokens.at(i)+(i<tokens.size()-1?".":"");
	}
      }
      // now operate on the DIRECTORIES
      if(!directory.empty() && vflags.flag("FLAG::AUID")) {
	if(XHOST_LIBRARY_ICSD!=LIBRARY_NOTHING) { // only if ICSD is defined
	  if(aurostd::substring2bool(directory,"/common/ICSD/LIB/") || aurostd::substring2bool(directory,"AFLOWDATA/ICSD_WEB/")) {
	    aurostd::StringSubst(directory,"/common/ICSD/LIB/","");
	    aurostd::StringSubst(directory,"AFLOWDATA/ICSD_WEB/","");
	    directory_LIB=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_ICSD)+"/LIB/"+directory;
	    directory_WEB=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_ICSD)+"/WEB/"+directory;
	    directory_RAW=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_ICSD)+"/RAW/"+directory;
	    url_WEB="/AFLOWDATA/ICSD_WEB/"+directory;
	    aurostd::string2tokens(directory,tokens,"/");
	    label=tokens.at(tokens.size()-1);
	  }
	} // ICSD
	if(XHOST_LIBRARY_LIB1!=LIBRARY_NOTHING) { // only if LIB1 is defined
	  if(aurostd::substring2bool(directory,"AFLOWDATA/LIB1_RAW/")) {
	    aurostd::StringSubst(directory,"AFLOWDATA/LIB1_RAW/","");
	    directory_LIB=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB1)+"/LIB/"+directory;
	    directory_WEB=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB1)+"/WEB/"+directory;
	    directory_RAW=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB1)+"/RAW/"+directory;
	    url_WEB="/AFLOWDATA/LIB1_WEB/"+directory;
	    aurostd::string2tokens(directory,tokens,"/");
	    label=tokens.at(tokens.size()-2)+"."+tokens.at(tokens.size()-1);
	  }
	} // LIB1
	if(XHOST_LIBRARY_LIB2!=LIBRARY_NOTHING) { // only if LIB2 is defined
	  if(aurostd::substring2bool(directory,"AFLOWDATA/LIB2_RAW/")) {
	    aurostd::StringSubst(directory,"AFLOWDATA/LIB2_RAW/","");
	    directory_LIB=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB2)+"/LIB/"+directory;
	    directory_WEB=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB2)+"/RAW/"+directory;  // WEB=RAW in lib2 June 6 2016
	    directory_RAW=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB2)+"/RAW/"+directory;
	    url_WEB="/AFLOWDATA/LIB2_RAW/"+directory;  // WEB=RAW in lib2 May 21 2014
	    aurostd::string2tokens(directory,tokens,"/");
	    label=tokens.at(tokens.size()-2)+"."+tokens.at(tokens.size()-1);
	  }
	} // LIB2
	if(XHOST_LIBRARY_LIB3!=LIBRARY_NOTHING) { // only if LIB3 is defined
	  if(aurostd::substring2bool(directory,"AFLOWDATA/LIB3_RAW/")) {
	    aurostd::StringSubst(directory,"AFLOWDATA/LIB3_RAW/","");
	    directory_LIB=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB3)+"/LIB/"+directory;
	    directory_WEB=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB3)+"/WEB/"+directory;
	    directory_RAW=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB3)+"/RAW/"+directory;
	    url_WEB="/AFLOWDATA/LIB3_WEB/"+directory;
	    aurostd::string2tokens(directory,tokens,"/");
	    label=tokens.at(tokens.size()-2)+"."+tokens.at(tokens.size()-1);
	  }
	} // LIB3
	if(XHOST_LIBRARY_LIB4!=LIBRARY_NOTHING) { // only if LIB4 is defined
	  if(aurostd::substring2bool(directory,"AFLOWDATA/LIB4_RAW/")) {
	    aurostd::StringSubst(directory,"AFLOWDATA/LIB4_RAW/","");
	    directory_LIB=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB4)+"/LIB/"+directory;
	    directory_WEB=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB4)+"/WEB/"+directory;
	    directory_RAW=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB4)+"/RAW/"+directory;
	    url_WEB="/AFLOWDATA/LIB4_WEB/"+directory;
	    aurostd::string2tokens(directory,tokens,"/");
	    label=tokens.at(tokens.size()-2)+"."+tokens.at(tokens.size()-1);
	  }
	} // LIB4
	if(XHOST_LIBRARY_LIB5!=LIBRARY_NOTHING) { // only if LIB5 is defined
	  if(aurostd::substring2bool(directory,"AFLOWDATA/LIB5_RAW/")) {
	    aurostd::StringSubst(directory,"AFLOWDATA/LIB5_RAW/","");
	    directory_LIB=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB5)+"/LIB/"+directory;
	    directory_WEB=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB5)+"/WEB/"+directory;
	    directory_RAW=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB5)+"/RAW/"+directory;
	    url_WEB="/AFLOWDATA/LIB5_WEB/"+directory;
	    aurostd::string2tokens(directory,tokens,"/");
	    label=tokens.at(tokens.size()-2)+"."+tokens.at(tokens.size()-1);
	  }
	} // LIB5
	if(XHOST_LIBRARY_LIB6!=LIBRARY_NOTHING) { // only if LIB6 is defined
	  if(aurostd::substring2bool(directory,"AFLOWDATA/LIB6_RAW/")) {
	    aurostd::StringSubst(directory,"AFLOWDATA/LIB6_RAW/","");
	    directory_LIB=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB6)+"/LIB/"+directory;
	    directory_WEB=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB6)+"/WEB/"+directory;
	    directory_RAW=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB6)+"/RAW/"+directory;
	    url_WEB="/AFLOWDATA/LIB6_WEB/"+directory;
	    aurostd::string2tokens(directory,tokens,"/");
	    label=tokens.at(tokens.size()-2)+"."+tokens.at(tokens.size()-1);
	  }
	} // LIB6
	if(XHOST_LIBRARY_LIB7!=LIBRARY_NOTHING) { // only if LIB7 is defined
	  if(aurostd::substring2bool(directory,"AFLOWDATA/LIB7_RAW/")) {
	    aurostd::StringSubst(directory,"AFLOWDATA/LIB7_RAW/","");
	    directory_LIB=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB7)+"/LIB/"+directory;
	    directory_WEB=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB7)+"/WEB/"+directory;
	    directory_RAW=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB7)+"/RAW/"+directory;
	    url_WEB="/AFLOWDATA/LIB7_WEB/"+directory;
	    aurostd::string2tokens(directory,tokens,"/");
	    label=tokens.at(tokens.size()-2)+"."+tokens.at(tokens.size()-1);
	  }
	} // LIB7
	if(XHOST_LIBRARY_LIB8!=LIBRARY_NOTHING) { // only if LIB8 is defined
	  if(aurostd::substring2bool(directory,"AFLOWDATA/LIB8_RAW/")) {
	    aurostd::StringSubst(directory,"AFLOWDATA/LIB8_RAW/","");
	    directory_LIB=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB8)+"/LIB/"+directory;
	    directory_WEB=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB8)+"/WEB/"+directory;
	    directory_RAW=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB8)+"/RAW/"+directory;
	    url_WEB="/AFLOWDATA/LIB8_WEB/"+directory;
	    aurostd::string2tokens(directory,tokens,"/");
	    label=tokens.at(tokens.size()-2)+"."+tokens.at(tokens.size()-1);
	  }
	} // LIB8
	if(XHOST_LIBRARY_LIB9!=LIBRARY_NOTHING) { // only if LIB9 is defined
	  if(aurostd::substring2bool(directory,"AFLOWDATA/LIB9_RAW/")) {
	    aurostd::StringSubst(directory,"AFLOWDATA/LIB9_RAW/","");
	    directory_LIB=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB9)+"/LIB/"+directory;
	    directory_WEB=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB9)+"/WEB/"+directory;
	    directory_RAW=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB9)+"/RAW/"+directory;
	    url_WEB="/AFLOWDATA/LIB9_WEB/"+directory;
	    aurostd::string2tokens(directory,tokens,"/");
	    label=tokens.at(tokens.size()-2)+"."+tokens.at(tokens.size()-1);
	  }
	} // LIB9
	if(LDEBUG) cout << "WEB_Aflowlib_Entry_PHP: label=" << label << "<br>" << endl;
	if(LDEBUG) cout << "WEB_Aflowlib_Entry_PHP: directory=" << directory << "<br>" << endl;
	if(LDEBUG) cout << "WEB_Aflowlib_Entry_PHP: directory_LIB=" << directory_LIB << "<br>" << endl;
	if(LDEBUG) cout << "WEB_Aflowlib_Entry_PHP: directory_WEB=" << directory_WEB << "<br>" << endl;
	if(LDEBUG) cout << "WEB_Aflowlib_Entry_PHP: directory_RAW=" << directory_RAW << "<br>" << endl;
	//	exit(0);
      }
      if(!directory.empty() && vflags.flag("FLAG::ICSD")) {
	directory_LIB=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_ICSD)+"/LIB/"+directory;
	directory_WEB=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_ICSD)+"/WEB/"+directory;
	directory_RAW=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_ICSD)+"/RAW/"+directory;
	url_WEB="/AFLOWDATA/ICSD_WEB/"+directory;
	label=option;
      }
      if(!directory.empty() && vflags.flag("FLAG::LIB1")) {
	directory_LIB=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB1)+"/LIB/"+directory;
	directory_WEB=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB1)+"/WEB/"+directory;
	directory_RAW=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB1)+"/RAW/"+directory;
	url_WEB="/AFLOWDATA/LIB1_RAW/"+directory;
	label=option;
      }
	if(!directory.empty() && vflags.flag("FLAG::LIB2")) {
	directory_LIB=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB2)+"/LIB/"+directory;
	directory_WEB=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB2)+"/RAW/"+directory; // June 2016
	directory_RAW=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB2)+"/RAW/"+directory;
	url_WEB="/AFLOWDATA/LIB2_RAW/"+directory; // May 2014
	label=option;
      }
      if(!directory.empty() && vflags.flag("FLAG::LIB3")) {
	directory_LIB=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB3)+"/LIB/"+directory;
	directory_WEB=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB3)+"/WEB/"+directory;
	directory_RAW=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB3)+"/RAW/"+directory;
	url_WEB="/AFLOWDATA/LIB3_WEB/"+directory;
 	label=option;
      }
      if(!directory.empty() && vflags.flag("FLAG::LIB4")) {
	directory_LIB=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB4)+"/LIB/"+directory;
	directory_WEB=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB4)+"/WEB/"+directory;
	directory_RAW=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB4)+"/RAW/"+directory;
	url_WEB="/AFLOWDATA/LIB4_WEB/"+directory;
 	label=option;
      }
      if(!directory.empty() && vflags.flag("FLAG::LIB5")) {
	directory_LIB=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB5)+"/LIB/"+directory;
	directory_WEB=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB5)+"/WEB/"+directory;
	directory_RAW=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB5)+"/RAW/"+directory;
	url_WEB="/AFLOWDATA/LIB5_WEB/"+directory;
 	label=option;
      }
      if(!directory.empty() && vflags.flag("FLAG::LIB6")) {
	directory_LIB=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB6)+"/LIB/"+directory;
	directory_WEB=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB6)+"/WEB/"+directory;
	directory_RAW=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB6)+"/RAW/"+directory;
	url_WEB="/AFLOWDATA/LIB6_WEB/"+directory;
 	label=option;
      }
      if(!directory.empty() && vflags.flag("FLAG::LIB7")) {
	directory_LIB=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB7)+"/LIB/"+directory;
	directory_WEB=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB7)+"/WEB/"+directory;
	directory_RAW=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB7)+"/RAW/"+directory;
	url_WEB="/AFLOWDATA/LIB7_WEB/"+directory;
 	label=option;
      }
      if(!directory.empty() && vflags.flag("FLAG::LIB8")) {
	directory_LIB=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB8)+"/LIB/"+directory;
	directory_WEB=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB8)+"/WEB/"+directory;
	directory_RAW=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB8)+"/RAW/"+directory;
	url_WEB="/AFLOWDATA/LIB8_WEB/"+directory;
 	label=option;
      }
      if(!directory.empty() && vflags.flag("FLAG::LIB9")) {
	directory_LIB=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB9)+"/LIB/"+directory;
	directory_WEB=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB9)+"/WEB/"+directory;
	directory_RAW=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB9)+"/RAW/"+directory;
	url_WEB="/AFLOWDATA/LIB9_WEB/"+directory;
 	label=option;
      }
    }
    
    // now start
    // got it  ?
    if(!aurostd::FileExist(directory_RAW+"/"+_AFLOWIN_)) directory_RAW="";
    
    if(!directory.empty()) { // play with aentry.entry
      aentry.file2aflowlib(directory_RAW+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT,oss);  //   oss << aentry.entry << endl;
      aurostd::string2tokens(aentry.sg2,tokens,"#");if(tokens.size()>0) aentry.sg2=tokens.at(tokens.size()-1);
      if(aentry.vfiles_WEB.size()==0) aentry.vfiles_WEB=aentry.vfiles;
    }
    
    // check AGL/AEL
    vflags.flag("FLAG::ELECTRONIC",aurostd::substring2bool(aentry.vloop,"bands"));
    vflags.flag("FLAG::SCINTILLATION",aurostd::substring2bool(aentry.vloop,"bands"));
    vflags.flag("FLAG::AGL",aurostd::substring2bool(aentry.vloop,"agl"));
    vflags.flag("FLAG::AEL",aurostd::substring2bool(aentry.vloop,"ael"));
    vflags.flag("FLAG::BADER",aurostd::substring2bool(aentry.vloop,"bader"));
    
    if(XHOST.hostname=="nietzsche.mems.duke.edu") {
      oss << "<b>DEBUG: only in " << XHOST.hostname << "</b><br>" << endl;
      oss << "XHOST.hostname=" << XHOST.hostname << "<br>" << endl;
      oss << "option=" << option << "<br>" << endl;
      oss << "label=" << label << "<br>" << endl;
      oss << "directory=" << directory << "<br>" << endl;
      oss << "vflags.flag(\"FLAG::ICSD\")=" << vflags.flag("FLAG::ICSD") << "<br>" << endl;
      oss << "vflags.flag(\"FLAG::LIB1\")=" << vflags.flag("FLAG::LIB1") << "<br>" << endl;
      oss << "vflags.flag(\"FLAG::LIB2\")=" << vflags.flag("FLAG::LIB2") << "<br>" << endl;
      oss << "vflags.flag(\"FLAG::LIB3\")=" << vflags.flag("FLAG::LIB3") << "<br>" << endl;
      oss << "vflags.flag(\"FLAG::LIB4\")=" << vflags.flag("FLAG::LIB4") << "<br>" << endl;
      oss << "vflags.flag(\"FLAG::LIB5\")=" << vflags.flag("FLAG::LIB5") << "<br>" << endl;
      oss << "vflags.flag(\"FLAG::LIB6\")=" << vflags.flag("FLAG::LIB6") << "<br>" << endl;
      oss << "vflags.flag(\"FLAG::LIB7\")=" << vflags.flag("FLAG::LIB7") << "<br>" << endl;
      oss << "vflags.flag(\"FLAG::LIB8\")=" << vflags.flag("FLAG::LIB8") << "<br>" << endl;
      oss << "vflags.flag(\"FLAG::LIB9\")=" << vflags.flag("FLAG::LIB9") << "<br>" << endl;
      oss << "vflags.flag(\"FLAG::AUID\")=" << vflags.flag("FLAG::AUID") << "<br>" << endl;
      oss << "vflags.flag(\"FLAG::FOUND\")=" << vflags.flag("FLAG::FOUND") << "<br>" << endl;
      oss << "vAUID.size()=" << vAUID.size() << "<br>" << endl;
      oss << "vAURL.size()=" << vAURL.size() << "<br>" << endl;
      oss << "directory_LIB=" << directory_LIB << "<br>" << endl;
      oss << "directory_WEB=" << directory_WEB << "<br>" << endl;
      oss << "directory_RAW=" << directory_RAW << "<br>" << endl;
      oss << "vflags.flag(\"FLAG::PREAMBLE\")=" << vflags.flag("FLAG::PREAMBLE")  << "<br>" << endl;
      oss << "vflags.flag(\"FLAG::CALCULATION\")=" << vflags.flag("FLAG::CALCULATION") << "<br>" << endl;
      oss << "vflags.flag(\"FLAG::JMOL\")=" << vflags.flag("FLAG::JMOL")  << "<br>" << endl;
      oss << "vflags.flag(\"FLAG::EDATA_ORIG\")=" << vflags.flag("FLAG::EDATA_ORIG") << "<br>" << endl;
      oss << "vflags.flag(\"FLAG::EDATA_RELAX\")=" << vflags.flag("FLAG::EDATA_RELAX") << "<br>" << endl;
      oss << "vflags.flag(\"FLAG::THERMODYNAMICS\")=" << vflags.flag("FLAG::THERMODYNAMICS") << "<br>" << endl;
      oss << "vflags.flag(\"FLAG::MAGNETIC\")=" << vflags.flag("FLAG::MAGNETIC") << "<br>" << endl;
      oss << "vflags.flag(\"FLAG::ELECTRONIC\")=" << vflags.flag("FLAG::ELECTRONIC") << "<br>" << endl;
      oss << "vflags.flag(\"FLAG::SCINTILLATION\")=" << vflags.flag("FLAG::SCINTILLATION") << "<br>" << endl;
      oss << "vflags.flag(\"FLAG::AGL\")=" << vflags.flag("FLAG::AGL") << "<br>" << endl;
      oss << "vflags.flag(\"FLAG::AEL\")=" << vflags.flag("FLAG::AEL") << "<br>" << endl;
      oss << "vflags.flag(\"FLAG::BADER\")=" << vflags.flag("FLAG::BADER") << "<br>" << endl;
      oss << "aentry.loop=" << aentry.loop << "<br>" << endl;
    }
    
    //CO 180523 - fixing for LIB6 missing from /www directory
    if(aurostd::substring2bool(XHOST.hostname, "aflowlib")){
      vflags.flag("FLAG::JMOL",FALSE);
      string web_path="/www"+url_WEB;
      if(LDEBUG){cerr << soliloquy << " web_path=" << web_path << endl;}
      if(aurostd::IsDirectory(web_path)){
        vector<string> web_path_vfiles;
        aurostd::DirectoryLS(web_path,web_path_vfiles);
        if(aurostd::substring2bool(web_path_vfiles,label+".cif")){vflags.flag("FLAG::JMOL",TRUE);}
      }
      if(LDEBUG){cerr << soliloquy << " vflags.flag(\"FLAG::JMOL\")=" << vflags.flag("FLAG::JMOL") << endl;}
    }
    
    // make ORIG vs RELAX
    vflags.flag("FLAG::EDATA_ORIG",FALSE);
    if(aentry.Bravais_lattice_orig!=aentry.Bravais_lattice_relax && aentry.Bravais_lattice_orig!="nan" && aentry.Bravais_lattice_relax!="nan") vflags.flag("FLAG::EDATA_ORIG",TRUE);
    if(aentry.lattice_variation_orig!=aentry.lattice_variation_relax && aentry.lattice_variation_orig!="nan" && aentry.lattice_variation_relax!="nan") vflags.flag("FLAG::EDATA_ORIG",TRUE);
    if(aentry.lattice_system_orig!=aentry.lattice_system_relax && aentry.lattice_system_orig!="nan" && aentry.lattice_system_relax!="nan") vflags.flag("FLAG::EDATA_ORIG",TRUE);
    if(aentry.Pearson_symbol_orig!=aentry.Pearson_symbol_relax && aentry.Pearson_symbol_orig!="nan" && aentry.Pearson_symbol_relax!="nan") vflags.flag("FLAG::EDATA_ORIG",TRUE);

    // ***************************************************************************
    // not found
    // oss << "Thank you for your query, but unfortunately entry " << label << " has not yet been calculated. If you want to report this omission please email Dr. Rose at aflowdev@aflowlib.duke.ed.<br>" << endl;
    // ***************************************************************************
    // PREAMBLE BEGIN
    string title=label;
    if(vflags.flag("FLAG::FOUND") && vflags.flag("FLAG::PREAMBLE") && !directory.empty()) {
      if(vflags.flag("FLAG::ICSD")) {
	aurostd::string2tokens(label,tokens,"_");
	title=tokens.at(0);
	for(uint i=0;i<=9;i++) aurostd::StringSubst(title,aurostd::utype2string<uint>(i),string("<sub>"+aurostd::utype2string<uint>(i)+"</sub>"));
	title+=" (ICSD# "+tokens.at(2)+")";//+directory;
      }
    }
    if(vflags.flag("FLAG::ELECTRONIC")){ // CO 180502
    oss.setf(std::ios::fixed,std::ios::floatfield);
    oss.precision(4);
    oss << "<script type=\"text/javascript\">" << endl;
    aurostd::xoption aaa;
    stringstream bandsdata;
      oss << "var d3_bands_data = "; estructure::BANDSDATA_JSON(aaa, directory_LIB, bandsdata,true); oss << bandsdata.str();  //GEENA
    oss << ";</script>" << endl;
    }
    oss << "<! HARVEY WORK BEFORE HERE> " << endl;
    oss << "<div id=\"content\">" << endl;
    oss << "<div class=\"title\">" << endl;
    oss << "<FONT SIZE=+3> " << title << " </FONT></div>" << endl;
    // ***************************************************************************
    // COMPOUND
    if(vflags.flag("FLAG::FOUND") && !directory.empty()) {
      oss << "<!-- compound: BEGIN -->" << endl;
      oss << "<FONT SIZE=+3> " << aentry.compound << " </FONT><br>" << endl;
       oss << "<!-- compound: END -->" << endl;
    }	
    // ***************************************************************************
    // LICENSE
    if(vflags.flag("FLAG::FOUND") && !directory.empty()) {
      string LICENSE="The data included within the aflow.org repository is free for scientific, academic and non-commercial purposes. Any other use is prohibited.";
      oss << "<!-- license: BEGIN -->" << endl;
      oss << "<div class = \"url_text\">" << endl;
      oss << "<br><b>LICENSE: <span class=\"url_text\">" << LICENSE << "</span></b>" << endl;
      oss << "</div>" << endl;
      oss << "<!-- license: END -->" << endl;
    }	
    // ***************************************************************************
    // PREAMBLE BEGIN
    if(vflags.flag("FLAG::FOUND") && vflags.flag("FLAG::PREAMBLE") && !directory.empty()) {
      //     string URL=string("http://aflow.org/material.php?proto_name=")+label;
      string URL=string("http://aflow.org/material.php?id=")+label;
      oss << "<!-- preamble: BEGIN -->" << endl;
      //      oss << "<FONT SIZE=+3> " << aentry.compound << " </FONT>" << endl;
      oss << "<div class = \"url\">" << endl;
      oss << "Permanent URL: <span class=\"url_text\">" << URL << "</span>" << endl;
      oss << "</div>" << endl;
      oss << "<br><FONT SIZE=+0>AFLOW.ORG web entry generator V" << string(AFLOW_VERSION) << " [built="  << TODAY << "]" << "</font>" << endl;
      oss << "<!-- preamble: END -->" << endl;
    }	
    // ***************************************************************************
    // NEW JSMOL BEGIN
    if(vflags.flag("FLAG::FOUND") && vflags.flag("FLAG::JMOL") && !directory.empty() && aurostd::substring2bool(aentry.vfiles_WEB,label+".cif")) {  //CO 180523 - fixing for LIB6 missing from /www directory
      //[OBSOLETE CO 170628 - per Bob/JMOL]if label+".cif" is available, assume "_sprim" and "_sconv" are too

      //[OBSOLETE CO 170628 - per Bob/JMOL]space group stuff found in bob's file now
	    //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "<br><br><b>Space Group</b>: " << (aurostd::string2utype<int>(aentry.spacegroup_relax)>0 ? GetSpaceGroupName(aurostd::string2utype<int>(aentry.spacegroup_relax)) : "N/A" ) << "  (#" << aentry.spacegroup_relax << ")" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]if(aurostd::string2utype<int>(aentry.spacegroup_relax)>0) {
	    //[OBSOLETE CO 170628 - per Bob/JMOL]  oss << "<br><br><b>Space Group</b>: " <<  GetSpaceGroupName(aurostd::string2utype<int>(aentry.spacegroup_relax)) << "  (#" << aentry.spacegroup_relax << ")" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]} else {
  	  //[OBSOLETE CO 170628 - per Bob/JMOL]  oss << "<br><br><b>Space Group</b>: " <<  "N/A" << "  (#" << aentry.spacegroup_relax << ")" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]}
      oss << "<!-- jmol: BEGIN -->" << endl;
      oss << "<!--div class = \"jmol\"-->" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "<script type=\"text/javascript\" src=\"/search/Lib/JS/JSmol.min.js\"></script>" << endl;  // CO 170622
      //[OBSOLETE CO 170628 - per Bob/JMOL]string JMOL_PATH="http://aflowlib.duke.edu/users/jmolers/test/jsmol";
      string JMOL_PATH="http://aflowlib.duke.edu/search/Lib/JS/JSMol";
      oss << "<script type=\"text/javascript\" src=\"" << JMOL_PATH << "/JSmol.min.js\"></script>" << endl;  // CO 170622
      oss << "<script type=\"text/javascript\">" << endl;
      // CO 170622 - START
      //build our standard AFLOW object, add from aentry as needed
      oss << "AFLOW={};" << endl;
      oss << "AFLOW.version = \"" << string(AFLOW_VERSION) << "\";" << endl;
      oss << "AFLOW.url_WEB = \"http://aflowlib.duke.edu" << url_WEB << "\";" << endl; // CO check with bob, maybe aflowlib.duke.edu?
      string system_name=KBIN::ExtractSystemName(directory_LIB);
      oss << "AFLOW.label = \"" << system_name << "\";" << endl;
      oss << "AFLOW.spaceGroupNo = " << aentry.spacegroup_relax << ";" << endl;
      oss << "AFLOW.spaceGroupName = \"" << GetSpaceGroupName(aurostd::string2utype<int>(aentry.spacegroup_relax)) << "\";" << endl;
      double aCIF,bCIF,cCIF,alphaCIF,betaCIF,gammaCIF;
      cif2data(directory_WEB+"/"+label+"_sconv.cif",aCIF,bCIF,cCIF,alphaCIF,betaCIF,gammaCIF);
      oss << "AFLOW.cif_sconv = [" << aCIF << "," << bCIF << "," << cCIF << "," << alphaCIF << "," << betaCIF << "," << gammaCIF << "];" << endl;
      cif2data(directory_WEB+"/"+label+".cif",aCIF,bCIF,cCIF,alphaCIF,betaCIF,gammaCIF);
      oss << "AFLOW.cif = [" << aCIF << "," << bCIF << "," << cCIF << "," << alphaCIF << "," << betaCIF << "," << gammaCIF << "];" << endl;
      cif2data(directory_WEB+"/"+label+"_sprim.cif",aCIF,bCIF,cCIF,alphaCIF,betaCIF,gammaCIF);
      oss << "AFLOW.cif_sprim = [" << aCIF << "," << bCIF << "," << cCIF << "," << alphaCIF << "," << betaCIF << "," << gammaCIF << "];" << endl;
	    string sym2json; //PC + DX 180723
	    if(aurostd::FileExist(directory_RAW+"/"+"aflow.fgroup.bands.json.xz")) //PC + DX 180723
      { aurostd::xzfile2string(directory_RAW+"/"+"aflow.fgroup.bands.json.xz",sym2json); //PC + DX 180723
      oss << "AFLOW.sym2json ="; //PC + DX 180723
      oss << sym2json;    // PC + DX 180723
      oss << ";" << endl; } //PC + DX 180723
	    else if (aurostd::FileExist(directory_RAW+"/"+"aflow.fgroup.relax.json.xz")) //PC + DX 180723
      { aurostd::xzfile2string(directory_RAW+"/"+"aflow.fgroup.relax.json.xz",sym2json); //PC + DX 180723
      oss << "AFLOW.sym2json ="; //PC + DX 180723
      oss << sym2json;    // PC + DX 180723
      oss << ";" << endl; }//PC + DX 180723
      else {cerr << "error" << endl; //PC + DX 180723
      }; //PC + DX 180723
      //BEGIN BADER ISOSURFACES
      if(vflags.flag("FLAG::BADER")){ //did we calculate bader?
        if(aurostd::substring2bool(aentry.vfiles_WEB,label+"_Bader_20_"+aentry.vspecies.at(0)+".jvxl")) { //quick (not robust) test that bader loop ran fine
          if(aurostd::substring2bool(aentry.vfiles_WEB,"CONTCAR.relax")) {  //check that we have the right structure
	xstructure xstr(directory_WEB+"/CONTCAR.relax",IOAFLOW_AUTO);
	xstr.ReScale(1.0);
            oss << "AFLOW.baderUnitcell = [";
            string sep = "";
	for(uint i=1;i<=3;i++)
              for(uint j=1;j<=3;j++) {
                oss << sep << xstr.lattice(i,j);
                sep = ",";
	  }
            oss << "];" << endl;
            oss << "AFLOW.baderVSpecies = [" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(aentry.vspecies,"\""),",") << "];" << endl;
	  }
	}
      }
      //END BADER ISOSURFACES
      oss << "AFLOW.jsmolDir = \"" << JMOL_PATH << "\";" << endl;
      
      //adding bob's stuff
      //string aflow_entry_js=AFLOW_ENTRY_JS;
      oss << AFLOW_WEBAPP_ENTRY_JS;  // CO 170622 
      
      // CO 170622 - END
      // CO all that follows is obsolete per bob's js file

      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "var jmolApplet0; // set up in HTML table, below" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "var use = \"HTML5\" // JAVA HTML5 WEBGL IMAGE  are all otions" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "var s = document.location.search;" << endl;

      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "jmol_isReady = function(applet) {" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "        Jmol._getElement(applet, \"appletdiv\").style.border=\"1px solid blue\"" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "}               " << endl;

      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "var Info = {" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "        width:  Math.floor(window.innerWidth*0.60)," << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "        height: Math.floor(window.innerWidth*0.60)," << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "        debug: false," << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "        color: \"black\"," << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "        addSelectionOptions: false," << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "        serverURL: \"http://aflowlib.duke.edu/search/Lib/PHP/jsmol.php\"," << endl;
      ////[OBSOLETE CO 170628 - per Bob/JMOL]oss << "        serverURL: \"search/Lib/PHP/jsmol.php\"," << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "        use: \"HTML5\"," << endl;
      ////[OBSOLETE CO 170628 - per Bob/JMOL]oss << "        j2sPath: \"/search/Lib/JS/j2s\", //Path" << endl; // CO keep
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "        j2sPath: \"../test/jsmol/j2s\", //Path" << endl; // CO remove
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "        readyFunction: jmol_isReady," << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "        script: \"set antialiasDisplay; set platformspeed 3; frank off; unitcell BOUNDBOX; set showUnitCellinfo false; load " << url_WEB << "/" << label << "_sconv.cif packed;"; // CO 170621 - speed up bader!
      //[OBSOLETE CO 170628 - per Bob/JMOL]aflowlib::cif2oss(directory_WEB+"/"+label+"_sconv.cif",label,aentry.spacegroup_relax,oss);
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "\"," << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "        //script: script," << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "        //jarPath: \"java\"," << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "        //jarFile: (useSignedApplet ? \"JmolAppletSigned.jar\" : \"JmolApplet.jar\")," << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "        //isSigned: useSignedApplet," << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "        //disableJ2SLoadMonitor: true," << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "        disableInitialConsole: true" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "}" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "// these are conveniences that mimic behavior of Jmol.js" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "function jmolCheckbox(script1, script0,text,ischecked) {Jmol.jmolCheckbox(jmolApplet0,script1, script0, text, ischecked)}" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "function jmolButton(script, text) {Jmol.jmolButton(jmolApplet0, script,text)}" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "function jmolHtml(s) { document.write(s) };" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "function jmolBr() { jmolHtml(\"<br />\") }" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "</script>" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "<table><tr>" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "<td align=center>" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "<script type=\"text/javascript\">" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "jmolApplet0 = Jmol.getApplet(\"jmolApplet0\", Info)" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "</script>" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "</td><td>" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "<form><!-- (FORM tag is important to automatically set checkbox off when page reloads) -->" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "<script type=\"text/javascript\">" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]//     oss << "<b>Space Group</b>" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]// oss << aentry.spacegroup_relax << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "Jmol.setButtonCss(null,\"style='width:110px'\")<!-- (Left buttons width) -->" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "jmolHtml(\"<b>Visualizer options:</b>\")" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "jmolBr()" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "jmolBr()" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "jmolButton(\"spacefill only;spacefill 23%;wireframe 0.15\",\"Ball & Stick\")" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "jmolButton(\"spacefill #alt:SETTING van der Waals Spheres\", \"Spacefill\");" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "jmolBr()" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "jmolButton(\"spin on\",\"Rotation On\")" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "jmolButton(\"spin off\",\"Rotation Off\")" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "jmolBr()" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "jmolButton(\"label  %a \",\"Label On\")" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "jmolButton(\"labels off \",\"Label Off\")" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "jmolBr()" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "jmolBr()" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "Jmol.setButtonCss(null,\"style='width:180px'\")<!-- (Left buttons width) -->" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "jmolHtml(\"<b>Relaxed structure:</b>\")" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "jmolBr()" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "jmolBr()" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "jmolButton(\"unitcell BOUNDBOX; load " << url_WEB << "/" << label <<".cif {1 1 1} packed;";
      //[OBSOLETE CO 170628 - per Bob/JMOL]aflowlib::cif2oss(directory_WEB+"/"+label+".cif",label,aentry.spacegroup_relax,oss);
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "  set zoomlarge false;zoom {*} 0\", \"As calculated\")" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "jmolBr()" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "jmolButton(\"unitcell BOUNDBOX; load " << url_WEB << "/" << label <<"_sconv.cif {1 1 1} packed;";
      //[OBSOLETE CO 170628 - per Bob/JMOL]aflowlib::cif2oss(directory_WEB+"/"+label+"_sconv.cif",label,aentry.spacegroup_relax,oss);
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "  set zoomlarge false;zoom {*} 0\", \"Standard conventional\")" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "jmolHtml(\"<a href=http://www.sciencedirect.com/science/article/pii/S0927025610002697>[info]</a>\")" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "jmolBr()" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "jmolButton(\"unitcell BOUNDBOX; load " << url_WEB << "/" << label <<"_sprim.cif {1 1 1} packed;";
      //[OBSOLETE CO 170628 - per Bob/JMOL]aflowlib::cif2oss(directory_WEB+"/"+label+"_sprim.cif",label,aentry.spacegroup_relax,oss);
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "  set zoomlarge false;zoom {*} 0\", \"Standard primitive\")" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "jmolHtml(\"<a href=http://www.sciencedirect.com/science/article/pii/S0927025610002697>[info]</a>\")" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "jmolBr()" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "jmolBr()" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "Jmol.setButtonCss(null,\"style='width:110px'\")<!-- (Left buttons width) -->" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "jmolHtml(\"<b>Supercell:</b>\")" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "jmolBr()" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "jmolBr()" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "jmolButton(\"unitcell BOUNDBOX; load '' {2 2 2} packed;  set zoomlarge false;zoom {*} 0\", \"2x2x2\")" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "jmolButton(\"unitcell BOUNDBOX; load '' fill 20;  set zoomlarge false;zoom {*} 0\", \"20&#8491; box\")" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]// BEGIN JJPR: USER CHOSES SUPERCELL SIZE
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "jmolHtml(\"<p><input class='dim' id='dim_1' type='number' value='2'> X <input class='dim' id='dim_2' type='number' value='2'> X <input class='dim' id='dim_3' type='number' value='2'></p><input type='button' id='build_button' value='BUILD' style='width:160px'>\") " << std::endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "  $( \"#build_button\" ).click(function() {" << std::endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "var scriptCommand = \"unitcell BOUNDBOX; load '' {\" +  $(\"#dim_1\").val() + \" \" + $(\"#dim_2\").val() + \" \" + $(\"#dim_3\").val() + \"} packed;  set zoomlarge false;zoom {*} 0\"; " << std::endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "Jmol.script(jmolApplet0, scriptCommand); " << std::endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "}); " << std::endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]// END JJPR: USER CHOSES SUPERCELL SIZE
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "jmolButton(\"unitcell BOUNDBOX; load " << url_WEB << "/" << label <<"_sconv.cif {1 1 1} packed;";
      //[OBSOLETE CO 170628 - per Bob/JMOL]aflowlib::cif2oss(directory_WEB+"/"+label+"_sconv.cif",label,aentry.spacegroup_relax,oss);
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "  set zoomlarge false;zoom {*} 0\", \"RESET\")" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "jmolBr()" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "jmolBr()" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "jmolHtml(\"<b>Crystallographic planes::</b>\")" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "jmolBr()" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "jmolBr()" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "jmolHtml(\"<p> h:<input class='dim' id='plane_1' type='number' value='2'>  k:<input class='dim' id='plane_2' type='number' value='2'>  l:<input class='dim' id='plane_3' type='number' value='2'></p><input type='button' id='plane_button' value='Show plane' style='width:160px'>\") " << std::endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "$( \"#plane_button\" ).click(function() { " << std::endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "var scriptCommand = \"unitcell BOUNDBOX; load '' {444 555 -1} packed;  isosurface hkl {\" +  $(\"#plane_1\").val() + \" \" + $(\"#plane_2\").val() + \" \" + $(\"#plane_3\").val() + \"} colorscheme sets translucent 0.5 green; set zoomlarge false;zoom {*} 0\";" << std::endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "Jmol.script(jmolApplet0, scriptCommand); " << std::endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "}); " << std::endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "jmolButton(\"unitcell BOUNDBOX; load " << url_WEB << "/" << label <<"_sconv.cif {1 1 1} packed;";
      //[OBSOLETE CO 170628 - per Bob/JMOL]aflowlib::cif2oss(directory_WEB+"/"+label+"_sconv.cif",label,aentry.spacegroup_relax,oss);
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "  set zoomlarge false;zoom {*} 0\", \"RESET\")" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "jmolBr()" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]oss << "jmolBr()" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]BEGIN BADER ISOSURFACES
      //[OBSOLETE CO 170628 - per Bob/JMOL]if(aurostd::substring2bool(aentry.vfiles_WEB,"CONTCAR.relax")) {
      //[OBSOLETE CO 170628 - per Bob/JMOL]if(LDEBUG) cerr << "BEGIN BADER ISOSURFACES" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]xstructure xstr(directory_WEB+"/CONTCAR.relax",IOAFLOW_AUTO);
      //[OBSOLETE CO 170628 - per Bob/JMOL]xstr.ReScale(1.0);
      //[OBSOLETE CO 170628 - per Bob/JMOL]stringstream unitcell;
      //[OBSOLETE CO 170628 - per Bob/JMOL]unitcell << std::fixed << std::setprecision(6); // CO 170621 - issues with JMOL
      //[OBSOLETE CO 170628 - per Bob/JMOL]string sep="";  // CO 170621 - issues with JMOL
      //[OBSOLETE CO 170628 - per Bob/JMOL]for(uint i=1;i<=3;i++){
      //[OBSOLETE CO 170628 - per Bob/JMOL]  for(uint j=1;j<=3;j++){
      //[OBSOLETE CO 170628 - per Bob/JMOL]    unitcell << sep << xstr.lattice(i,j); // CO 170621 - issues with JMOL
      //[OBSOLETE CO 170628 - per Bob/JMOL]    sep=",";  // CO 170621 - issues with JMOL
      //[OBSOLETE CO 170628 - per Bob/JMOL]  }
      //[OBSOLETE CO 170628 - per Bob/JMOL]}
      //[OBSOLETE CO 170628 - per Bob/JMOL]if(aurostd::substring2bool(aentry.vfiles_WEB,label+"_Bader_20_"+aentry.vspecies.at(0)+".jvxl")) {
      //[OBSOLETE CO 170628 - per Bob/JMOL]  oss << "Jmol.setButtonCss(null,\"style='width:50px'\")<!-- (Left buttons width) -->" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]  oss << "jmolHtml(\"<b>Bader Isosurfaces:</b>\")" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]  oss << "jmolBr()" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]  oss << "jmolHtml(\"Cutoff = 0.20\")" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]  oss << "jmolBr()" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]  for(uint j=0;j<aentry.vspecies.size();j++) {
      //[OBSOLETE CO 170628 - per Bob/JMOL]    oss << "jmolButton(\"unitcell BOUNDBOX; load " << url_WEB << "/" << label <<".cif {1 1 1} packed UNITCELL["<< unitcell.str() << "];";
      //[OBSOLETE CO 170628 - per Bob/JMOL]    aflowlib::cif2oss(directory_WEB+"/"+label+".cif",label,aentry.spacegroup_relax,oss);
      //[OBSOLETE CO 170628 - per Bob/JMOL]    aflowlib::iso2oss(url_WEB,label,aentry.vspecies.at(j),"20",j,oss);
      //[OBSOLETE CO 170628 - per Bob/JMOL]    oss << "  set zoomlarge false;zoom {*} 0\", \""<<aentry.vspecies.at(j)<<"\")" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]  }
      //[OBSOLETE CO 170628 - per Bob/JMOL]  oss << "jmolButton(\"unitcell BOUNDBOX; load " << url_WEB << "/" << label <<".cif {1 1 1} packed UNITCELL["<< unitcell.str()  << "];";
      //[OBSOLETE CO 170628 - per Bob/JMOL]  aflowlib::cif2oss(directory_WEB+"/"+label+".cif",label,aentry.spacegroup_relax,oss);
      //[OBSOLETE CO 170628 - per Bob/JMOL]  for(uint j=0;j<aentry.vspecies.size();j++) {
      //[OBSOLETE CO 170628 - per Bob/JMOL]    aflowlib::iso2oss(url_WEB,label,aentry.vspecies.at(j),"20",j,oss);
      //[OBSOLETE CO 170628 - per Bob/JMOL]  }
      //[OBSOLETE CO 170628 - per Bob/JMOL]  oss << "  set zoomlarge false;zoom {*} 0\", \"All\")" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]  oss << "jmolBr()" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]  oss << "jmolHtml(\"Cutoff = 0.30\")" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]  oss << "jmolBr()" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]  for(uint j=0;j<aentry.vspecies.size();j++) {
      //[OBSOLETE CO 170628 - per Bob/JMOL]    oss << "jmolButton(\"unitcell BOUNDBOX; load " << url_WEB << "/" << label <<".cif {1 1 1} packed UNITCELL["<< unitcell.str() << "];";
      //[OBSOLETE CO 170628 - per Bob/JMOL]    aflowlib::cif2oss(directory_WEB+"/"+label+".cif",label,aentry.spacegroup_relax,oss);
      //[OBSOLETE CO 170628 - per Bob/JMOL]    aflowlib::iso2oss(url_WEB,label,aentry.vspecies.at(j),"30",j,oss);
      //[OBSOLETE CO 170628 - per Bob/JMOL]    oss << "  set zoomlarge false;zoom {*} 0\", \""<<aentry.vspecies.at(j)<<"\")" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]  }
      //[OBSOLETE CO 170628 - per Bob/JMOL]  oss << "jmolButton(\"unitcell BOUNDBOX; load " << url_WEB << "/" << label <<".cif {1 1 1} packed UNITCELL["<< unitcell.str() << "];";
      //[OBSOLETE CO 170628 - per Bob/JMOL]  aflowlib::cif2oss(directory_WEB+"/"+label+".cif",label,aentry.spacegroup_relax,oss);
      //[OBSOLETE CO 170628 - per Bob/JMOL]  for(uint j=0;j<aentry.vspecies.size();j++) {
      //[OBSOLETE CO 170628 - per Bob/JMOL]    aflowlib::iso2oss(url_WEB,label,aentry.vspecies.at(j),"30",j,oss);
      //[OBSOLETE CO 170628 - per Bob/JMOL]  }
      //[OBSOLETE CO 170628 - per Bob/JMOL]  oss << "  set zoomlarge false;zoom {*} 0\", \"All\")" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]  oss << "jmolBr()" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]  oss << "jmolHtml(\"Cutoff = 0.40\")" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]  oss << "jmolBr()" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]  for(uint j=0;j<aentry.vspecies.size();j++) {
      //[OBSOLETE CO 170628 - per Bob/JMOL]    oss << "jmolButton(\"unitcell BOUNDBOX; load " << url_WEB << "/" << label <<".cif {1 1 1} packed UNITCELL["<< unitcell.str() << "];";
      //[OBSOLETE CO 170628 - per Bob/JMOL]    aflowlib::cif2oss(directory_WEB+"/"+label+".cif",label,aentry.spacegroup_relax,oss);
      //[OBSOLETE CO 170628 - per Bob/JMOL]    aflowlib::iso2oss(url_WEB,label,aentry.vspecies.at(j),"40",j,oss);
      //[OBSOLETE CO 170628 - per Bob/JMOL]    oss << "  set zoomlarge false;zoom {*} 0\", \""<<aentry.vspecies.at(j)<<"\")" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]  }
      //[OBSOLETE CO 170628 - per Bob/JMOL]  oss << "jmolButton(\"unitcell BOUNDBOX; load " << url_WEB << "/" << label <<".cif {1 1 1} packed UNITCELL["<< unitcell.str() << "];";
      //[OBSOLETE CO 170628 - per Bob/JMOL]  aflowlib::cif2oss(directory_WEB+"/"+label+".cif",label,aentry.spacegroup_relax,oss);
      //[OBSOLETE CO 170628 - per Bob/JMOL]  for(uint j=0;j<aentry.vspecies.size();j++) {
      //[OBSOLETE CO 170628 - per Bob/JMOL]    aflowlib::iso2oss(url_WEB,label,aentry.vspecies.at(j),"40",j,oss);
      //[OBSOLETE CO 170628 - per Bob/JMOL]  }
      //[OBSOLETE CO 170628 - per Bob/JMOL]  oss << "  set zoomlarge false;zoom {*} 0\", \"All\")" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]  oss << "jmolBr()" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]  oss << "jmolHtml(\"Cutoff = 0.50\")" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]  oss << "jmolBr()" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]  for(uint j=0;j<aentry.vspecies.size();j++) {
      //[OBSOLETE CO 170628 - per Bob/JMOL]    oss << "jmolButton(\"unitcell BOUNDBOX; load " << url_WEB << "/" << label <<".cif {1 1 1} packed UNITCELL["<< unitcell.str() << "];";
      //[OBSOLETE CO 170628 - per Bob/JMOL]    aflowlib::cif2oss(directory_WEB+"/"+label+".cif",label,aentry.spacegroup_relax,oss);
      //[OBSOLETE CO 170628 - per Bob/JMOL]    aflowlib::iso2oss(url_WEB,label,aentry.vspecies.at(j),"50",j,oss);
      //[OBSOLETE CO 170628 - per Bob/JMOL]    oss << "  set zoomlarge false;zoom {*} 0\", \""<<aentry.vspecies.at(j)<<"\")" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]  }
      //[OBSOLETE CO 170628 - per Bob/JMOL]  oss << "jmolButton(\"unitcell BOUNDBOX; load " << url_WEB << "/" << label <<".cif {1 1 1} packed UNITCELL["<< unitcell.str()  << "];";
      //[OBSOLETE CO 170628 - per Bob/JMOL]  aflowlib::cif2oss(directory_WEB+"/"+label+".cif",label,aentry.spacegroup_relax,oss);
      //[OBSOLETE CO 170628 - per Bob/JMOL]  for(uint j=0;j<aentry.vspecies.size();j++) {
      //[OBSOLETE CO 170628 - per Bob/JMOL]    aflowlib::iso2oss(url_WEB,label,aentry.vspecies.at(j),"50",j,oss);
      //[OBSOLETE CO 170628 - per Bob/JMOL]  }
      //[OBSOLETE CO 170628 - per Bob/JMOL]  oss << "  set zoomlarge false;zoom {*} 0\", \"All\")" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]  oss << "jmolBr()" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]  oss << "jmolBr()" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]}
      //[OBSOLETE CO 170628 - per Bob/JMOL]    }
      //[OBSOLETE CO 170628 - per Bob/JMOL]    //END BADER ISOSURFACES
      //[OBSOLETE CO 170628 - per Bob/JMOL]    oss << "Jmol.setButtonCss(null,\"style='width:140px'\")<!-- (Left buttons width) -->" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]    oss << "jmolHtml(\"<b>Save:</b>\")" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]    oss << "jmolBr()" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]    oss << "jmolBr()" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]    oss << "jmolButton(\"write FILE ?\",\"Save CIF FILE\")" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]    oss << "jmolButton(\"write STATE ?.spt\",\"Save STATE\")" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]    oss << "jmolBr()" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]    oss << "jmolButton(\"write IMAGE ?.jpg\",\"Save JPG\")" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]    oss << "jmolButton(\"write IMAGE ?.png\",\"Save PNG\")" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]    oss << "jmolBr()" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]    oss << "jmolButton(\"write ?.jmol\",\"Save Jmol\")" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]    oss << "jmolButton(\"write PNGJ ?.png\",\"Save PNG+Jmol\")" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]    oss << "jmolBr()" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]    oss << "</script>" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]    oss << "</form>" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]    oss << "</td>" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]    oss << "</tr>" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]    oss << "<tr>" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]    oss << "<td align=center>" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]    oss << "<script type=\"text/javascript\">" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]    oss << "jmolBr()" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]    oss << "jmolBr()" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]    oss << "Jmol.setButtonCss(null,\"style='width:100px'\")" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]    //oss << "jmolButton(\"console\")" << endl;
      //[OBSOLETE CO 170628 - per Bob/JMOL]    //oss << "Jmol.jmolCommandInput(jmolApplet0)" << endl;
      oss << "</script>" << endl;
      oss << "</td></tr></table>" << endl;
      oss << "<!--/div-->" << endl;
      oss << "<!-- jmol: END -->" << endl;
    }
    if(vflags.flag("FLAG::FOUND") && vflags.flag("FLAG::ELECTRONIC") && !directory.empty()) {
      oss << "<!-- geena bands: BEGIN -->" << endl;
      oss << "<script type=\"text/javascript\" src=\"./Lib/JS/d3.min.js\"></script>" << endl;  // CO 170622  ///www/search/Lib/JS/d3.min.js
      oss << "<script type=\"text/javascript\">" << endl;
      oss << AFLOW_WEBAPP_BANDS_JS; //PC 180515
      oss << "</script>" << endl;
      oss << "<!-- geena bands: END -->" << endl;
    }
    // NEW JSMOL END
    // ***************************************************************************
    // CALCULATION
    if(vflags.flag("FLAG::FOUND") && vflags.flag("FLAG::CALCULATION") && !directory.empty()) {
     // oss << line_rule << endl; //JPO 180731
      oss << "<!-- Calculation properties: BEGIN -->" << endl;
      oss << "<div class=\"container\">" << endl; //JPO 180731
      oss << "<div class=\"container-title\"><h1 class=\"section-title\">Calculation details</h1></div>" << endl; //JPO 180731
     // oss << line_rule << endl; //JPO 180731
      //      oss << "<div class=\"calculation_details\">" << endl; //JPO 180731
     // oss << "<div>" << endl; //JPO 180731
      //oss << "<ul>" << endl; //JPO 180731
      if(!aentry.code.empty()) {
	oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\">" << DEFAULT_FILE_AFLOWLIB_ENTRY_OUT << ":</h5><span class=\"value\">" << "[" << "<a href=\"" << url_WEB << "/\"" << html_TAB << ">entry</a>|" << "<a href=\"" << url_WEB << "/" << DEFAULT_FILE_AFLOWLIB_ENTRY_OUT << "\"" << html_TAB << ">raw</a>]" << "</span></div></div>" << endl; //JPO 180731
	oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\">" << DEFAULT_FILE_AFLOWLIB_ENTRY_JSON << ":</h5><span class=\"value\">" << "[" << "<a href=\"" << url_WEB << "/?format=json\"" << html_TAB << ">entry</a>|" << "<a href=\"" << url_WEB << "/" << DEFAULT_FILE_AFLOWLIB_ENTRY_JSON << "\"" << html_TAB << ">raw</a>]" << "</span></div></div>" << endl; //JPO 180731
      }
      if(!aentry.auid.empty()) 
	oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\">AFLOW-UID [<a href=\"" << url_WEB << "/?auid\">auid</a>]:</h5><span class=\"value\">" << aentry.auid << "</span></div></div>" << endl; // PC 180515  //JPO 180731
      if(!aentry.aurl.empty()) 
	oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\">AFLOW-URL [<a href=\"" << url_WEB << "/?aurl\">aurl</a>]:</h5><span class=\"value\">" << aentry.aurl << "</span></div></div>" << endl; // PC 180515 //JPO 180731
      if(!aentry.code.empty()) 
	oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"><i>Ab-initio</i> code [<a href=\"" << url_WEB << "/?code\">code</a>]:</h5><span class=\"value\">" << aentry.code << (aentry.calculation_cores>1?" (MPI) ":"") << "</span></div></div>" << endl;  // PC 180515 //JPO 180731
      if(!aentry.compound.empty()) 
	oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\">Composition per cell [<a href=\"" << url_WEB << "/?compound\">compound</a>]:</h5><span class=\"value\">" << aentry.compound << "</span></div></div>" << endl; // PC 180515 //JPO 180731
      if(!aentry.icsd.empty()) 
	oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\">ICSD entry" << icsd_link << ":</h5><span class=\"value\">" << aentry.icsd << "</span></div></div>" << endl;  // PC 180515 //JPO 180731
      if(!aentry.dft_type.empty()) 
	oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\">Pseudopotentials type [<a href=\"" << url_WEB << "/?pp_type\">pp_type</a>]:</h5><span class=\"value\">" << aentry.dft_type << "</span></div></div>" << endl;  // PC 180515 //JPO 180731
      for(uint i=0;i<aentry.vspecies_pp.size();i++)  
	oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\">PP - Species&middot;Version [<a href=\"" << url_WEB << "/?species_pp_version\">species_pp_version</a>]:</h5><span class=\"value\">" << aentry.vspecies_pp_version.at(i) << "</span></div></div>" << endl; // PC 180515 //JPO 180731
      for(uint i=0;i<aentry.vspecies_pp_ZVAL.size();i++)  
	oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\">PP - Species&middot;ZVAL [<a href=\"" << url_WEB << "/?species_pp_ZVAL\">species_pp_ZVAL</a>]:</h5><span class=\"value\">" << aentry.vspecies_pp_ZVAL.at(i) << "</span></div></div>" << endl; // PC 180515 //JPO 180731
  oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\">LDAU [T;{L};{U};{J}] [<a href=\"" << url_WEB << "/?ldau_TLUJ\">ldau_TLUJ</a>]:</h5><span class=\"value\">" << (!aentry.ldau_TLUJ.empty()?aentry.ldau_TLUJ+"  (type;{angular-1};{eV};{eV}) ":"no-LDAU") << "</span></div></div>" << endl; // PC 180515 //JPO 180731
      if(aentry.calculation_time*aentry.calculation_cores>0.0) 
	oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\">Total CPU&middot;hours [<a href=\"" << url_WEB << "/?calculation_time\">calculation_time</a>*<a href=\"" << url_WEB << "/?calculation_cores\">calculation_cores</a>]/3600:</h5><span class=\"value\">" << aentry.calculation_time*aentry.calculation_cores/3600 << " hours</span></div></div>" << endl;  // PC 180515 //JPO 180731
      if(aentry.calculation_time>0.0) 
	oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\">Total Wall-time [<a href=\"" << url_WEB << "/?calculation_time\">calculation_time/3600</a>]:</h5><span class=\"value\">" << aentry.calculation_time/3600 << " hours</span></div></div>" << endl;  // PC 180515 //JPO 180731
      if(aentry.calculation_memory>0.0)
	oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\">Memory Used [<a href=\"" << url_WEB << "/?calculation_memory\">calculation_memory</a>]:</h5><span class=\"value\">" << aentry.calculation_memory << " MB</span></div></div>" << endl;  // PC 180515 //JPO 180731
      if(aentry.calculation_cores>0) 
	oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\">Number of cores [<a href=\"" << url_WEB << "/?calculation_cores\">calculation_cores</a>]:</h5><span class=\"value\">" << aentry.calculation_cores << (aentry.calculation_cores>1?" (MPI) ":"(serial)") << "</span></div></div>" << endl; // PC 180515 //JPO 180731
      if(!aentry.node_CPU_Model.empty()) 
	oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\">Processor - CPU Model [<a href=\"" << url_WEB << "/?node_CPU_Model\">node_CPU_Model</a>]:</h5><span class=\"value\">" << (!aentry.node_CPU_Model.empty()?aentry.node_CPU_Model:"unavailable") << "</span></div></div>" << endl;  // PC 180515 //JPO 180731
      if(aentry.node_CPU_MHz>0.0) 
	oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\">Processor - CPU MHz [<a href=\"" << url_WEB << "/?node_CPU_MHz\">node_CPU_MHz</a>]:</h5><span class=\"value\">" << (aentry.node_CPU_MHz>1.0?aurostd::utype2string(aentry.node_CPU_MHz)+"MHz ":"unavailable") << " </span></div></div>" << endl;  // PC 180515 //JPO 180731
      if(!aentry.aflow_version.empty()) 
	oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\">AFLOW Version [<a href=\"" << url_WEB << "/?aflow_version\">aflow_version</a>]:</h5><span class=\"value\">" << (!aentry.aflow_version.empty()?aentry.aflow_version:"unavailable") << "</span></div></div>" << endl;  // PC 180515 //JPO 180731
      if(!aentry.catalog.empty()) 
	oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\">AFLOW Catalog [<a href=\"" << url_WEB << "/?catalog\">catalog</a>]:</h5><span class=\"value\">" << (!aentry.catalog.empty()?aentry.catalog:"unavailable") << "</span></div></div>" << endl;  // PC 180515 //JPO 180731
      if(!aentry.data_api.empty()) 
	oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\">AFLOW Data API [<a href=\"" << url_WEB << "/?data_api\">data_api</a>]:</h5><span class=\"value\">" << (!aentry.data_api.empty()?aentry.data_api:"unavailable") << "</span></div></div>" << endl;  // PC 180515 //JPO 180731
      if(!aentry.data_source.empty()) 
	oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\">AFLOW Data Source [<a href=\"" << url_WEB << "/?data_source\">data_source</a>]:</h5><span class=\"value\">" << (!aentry.data_source.empty()?aentry.data_source:"unavailable") << "</span></div></div>" << endl;  // PC 180515 //JPO 180731
      if(!aentry.data_language.empty()) 
	oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\">AFLOW Data Language [<a href=\"" << url_WEB << "/?data_language\">data_language</a>]:</h5><span class=\"value\">" << (!aentry.data_language.empty()?aentry.data_language:"unavailable") << "</span></div></div>" << endl;  // PC 180515 //JPO 180731
      if(!aentry.error_status.empty()) 
	oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\">AFLOW Error Status [<a href=\"" << url_WEB << "/?error_status\">error_status</a>]:</h5><span class=\"value\">" << (!aentry.error_status.empty()?aentry.error_status:"unavailable") << "</span></div></div>" << endl;  // PC 180515 //JPO 180731
      if(!aentry.loop.empty()) 
	oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\">AFLOW loops [<a href=\"" << url_WEB << "/?loop\">loop</a>]:</h5><span class=\"value\">" << aentry.loop << "</span></div></div>" << endl; //JPO 180731
      //     oss << "<li><span class=\"description\">AFLOW precision:</span>" << "unavailable" << "</li>" << endl; //JPO 180731
      //    oss << "<li><span class=\"description\">AFLOW KPPRA:</span>" << "unavailable" << "</li>" << endl; //JPO 180731
      if(!aentry.aflowlib_version.empty()) 
	oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\">AFLOWLIB_entry version [<a href=\"" << url_WEB << "/?aflowlib_version\">aflowlib_version</a>]:</h5><span class=\"value\">" << (!aentry.aflowlib_version.empty()?aentry.aflowlib_version:"unavailable") << "</span></div></div>" << endl; // PC 180515 //JPO 180731
      if(!aentry.aflowlib_date.empty()) 
	oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\">AFLOWLIB_entry date [<a href=\"" << url_WEB << "/?aflowlib_date\">aflowlib_date</a>]:</h5><span class=\"value\">" << (!aentry.aflowlib_date.empty()?aentry.aflowlib_date:"unavailable") << "</span></div></div>" << endl;  // PC 180515 //JPO 180731
      if(!aentry.author.empty()) 
	oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\">AFLOW Version [<a href=\"" << url_WEB << "/?author\">author</a>]:</h5><span class=\"value\">" << (!aentry.author.empty()?aentry.author:"unavailable") << "</span></div></div>" << endl;  // PC 180515 //JPO 180731
      if(!aentry.corresponding.empty()) 
	oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\">AFLOW Version [<a href=\"" << url_WEB << "/?corresponding\">corresponding</a>]:</h5><span class=\"value\">" << (!aentry.corresponding.empty()?aentry.corresponding:"unavailable") << "</span></div></div>" << endl;  // PC 180515 //JPO 180731
      if(!aentry.sponsor.empty()) 
	oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\">AFLOW Version [<a href=\"" << url_WEB << "/?sponsor\">sponsor</a>]:</h5><span class=\"value\">" << (!aentry.sponsor.empty()?aentry.sponsor:"unavailable") << "</span></div></div>" << endl;  // PC 180515 //JPO 180731

      // CORMAC
      if(aentry.energy_cutoff!=AUROSTD_NAN) 
	oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\">energy_cutoff [<a href=\"" << url_WEB << "/?energy_cutoff\">energy_cutoff</a>]:</h5><span class=\"value\">" << aentry.energy_cutoff << " eV</span></div></div>" << endl; // PC 180515 //JPO 180731
      //      oss.precision(6);
      if(aentry.delta_electronic_energy_convergence!=AUROSTD_NAN) 
	oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\">delta_electronic_energy_convergence [<a href=\"" << url_WEB << "/?delta_electronic_energy_convergence\">delta_electronic_energy_convergence</a>]:</h5><span class=\"value\">" << 1000*aentry.delta_electronic_energy_convergence << " meV</span></div></div>" << endl; // PC 180515 //JPO 180731
      if(aentry.delta_electronic_energy_threshold!=AUROSTD_NAN) 
	oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\">delta_electronic_energy_threshold_ [<a href=\"" << url_WEB << "/?delta_electronic_energy_threshold\">delta_electronic_energy_threshold</a>]:</h5><span class=\"value\">" << 1000*aentry.delta_electronic_energy_threshold << " meV</span></div></div>" << endl; // PC 180515 //JPO 180731
      if(aentry.nkpoints!=0) 
	oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\">nkpoints [<a href=\"" << url_WEB << "/?nkpoints\">nkpoints</a>]:</h5><span class=\"value\">" << aentry.nkpoints << " </span></div></div>" << endl; //JPO 180731
      if(aentry.nkpoints_irreducible!=0) 
	oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\">nkpoints_irreducible [<a href=\"" << url_WEB << "/?nkpoints_irreducible\">nkpoints_irreducible</a>]:</h5><span class=\"value\">" << aentry.nkpoints_irreducible << " </span></div></div>" << endl; //JPO 180731
      if(aentry.kppra!=0) 
	oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\">kppra [<a href=\"" << url_WEB << "/?kppra\">kppra</a>]:</h5><span class=\"value\">" << aentry.kppra << " </span></div></div>" << endl; //JPO 180731
      if(aentry.kpoints.empty()) 
	oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\">kpoints [<a href=\"" << url_WEB << "/?kpoints\">kpoints</a>]:</h5><span class=\"value\">" << aentry.kpoints << " </span></div></div>" << endl; //JPO 180731
      //  oss.precision(3);

      //oss << "</ul>" << endl; //JPO 180731
      oss << "</div>" << endl;  //JPO 180731
      oss << "<!-- Calculation properties: END -->" << endl;
    }

    // ***************************************************************************
    if(vflags.flag("FLAG::FOUND") && vflags.flag("FLAG::EDATA_ORIG")&& !directory.empty()) {
      oss << line_rule << endl;
      oss << "<!-- Warning: BEGIN -->" << endl;
      oss << "<span class=\"title\"><FONT SIZE=+3 color=red> Warning!</font></span>" << endl;
      oss << "<span class=\"title\"><FONT SIZE=+3 color=red> Original and relaxed structures</font></span>" << endl;
      oss << "<span class=\"title\"><FONT SIZE=+3 color=red> have different symmetries: listing both.</font></span>" << endl;
      oss << "<!-- Warning: END -->" << endl;
      oss << line_rule << endl;
    }

    // ***************************************************************************
    // EDATA ORIG/RELAX
    if((vflags.flag("FLAG::EDATA_ORIG") || vflags.flag("FLAG::EDATA_RELAX")) && !directory.empty()) {
      for(uint i=0;i<=1;i++)  {
	if((vflags.flag("FLAG::EDATA_ORIG") && i==0) || (vflags.flag("FLAG::EDATA_RELAX") && i==1)) {
	 // oss << line_rule << endl; //JPO 180731
	  oss << "<!-- EDATA: BEGIN -->" << endl;
    oss << "<div class=\"container\">" << endl; //JPO 180731
	  if(i==0) oss << "<div class=\"container-title\"><h1 class=\"section-title\"> Original Structure</h1></div>" << endl; //JPO 180731
	  if(i==1) oss << "<div class=\"container-title\"><h1 class=\"section-title\"> Relaxed Structure</h1></div>" << endl; //JPO 180731
	 // oss << line_rule << endl; //JPO 180731
	 // oss << "<ht /> " << endl; //JPO 180731
	 // oss << "<div class = \"real_space\">" << endl; //JPO 180731
	  if(i==0) oss << "<div class=\"container-subtitle\"><h4 class=\"section-subtitle\"> Real Space Lattice</h4></div>" << endl; //JPO 180731
	  if(i==1) oss << "<div class=\"container-subtitle\"><h4 class=\"section-subtitle\"> Real Space Lattice</h4></div>" << endl; //JPO 180731
	  vector<string> vline_edata;
	  if(i==0 && aurostd::FileExist(directory_RAW+"/"+DEFAULT_FILE_EDATA_ORIG_OUT)) aurostd::file2vectorstring(directory_RAW+"/"+DEFAULT_FILE_EDATA_ORIG_OUT,vline_edata);
	  if(i==1 && aurostd::FileExist(directory_RAW+"/"+DEFAULT_FILE_EDATA_RELAX_OUT)) aurostd::file2vectorstring(directory_RAW+"/"+DEFAULT_FILE_EDATA_RELAX_OUT,vline_edata);
	  // 
	  vector<double> abcR(6);
	  double volumeR,density=0.0,coveraR;
	  string Crystal_Real_space_Bravais_Lattice_Primitive="",Crystal_Real_space_Lattice_Variation="",Crystal_Real_space_Lattice_System="";
	  string Crystal_Real_space_Pearson_Symbol="",Crystal_Real_space_Crystal_Family="",Crystal_Real_space_Crystal_System="";
	  string Crystal_Real_space_Crystal_Class="",Crystal_Real_space_Point_Group_Hermann_Mauguin="",Crystal_Real_space_Point_Group_Schoenflies="";
	  string Crystal_Real_space_Point_Group_Orbifold="",Crystal_Real_space_Point_Group_Type="",Crystal_Real_space_Point_Group_Order="";
	  string Crystal_Real_space_Point_Group_Structure="";
	
	  string Lattice_Real_space_Bravais_Lattice_Primitive="",Lattice_Real_space_Lattice_Variation="",Lattice_Real_space_Lattice_System="";
	
	  string Superattice_Real_space_Bravais_Superlattice_Primitive="",Superattice_Real_space_Superlattice_Variation="",Superattice_Real_space_Superlattice_System="";
	  string Superattice_Real_space_Pearson_Symbol_Superlattice="",Reciprocal_lattice_primitive="",Reciprocal_lattice_variation="";
	
	  vector<double> abcK(6);
	  double volumeK;
	  if(vline_edata.size()>0) {
	    for(uint iline=0;iline<vline_edata.size();iline++) {
	      if(aurostd::substring2bool(vline_edata.at(iline),"Real space a b c alpha beta gamma"))
		if(!aurostd::substring2bool(vline_edata.at(iline),"Bohrs/Degs")) {
		  aurostd::string2tokens(vline_edata.at(iline),tokens);
		  for(uint i=0;i<6;i++) abcR.at(i)=aurostd::string2utype<double>(tokens.at(tokens.size()-6+i));}
	      if(aurostd::substring2bool(vline_edata.at(iline),"Real space Volume")) {
		aurostd::string2tokens(vline_edata.at(iline),tokens);
		volumeR=aurostd::string2utype<double>(tokens.at(tokens.size()-1));}
	      if(aurostd::substring2bool(vline_edata.at(iline),"Real space c/a")) {
		aurostd::string2tokens(vline_edata.at(iline),tokens);
		coveraR=aurostd::string2utype<double>(tokens.at(tokens.size()-1));}
	    
	      if(aurostd::substring2bool(vline_edata.at(iline),"BRAVAIS LATTICE OF THE CRYSTAL")) {
		aurostd::string2tokens(vline_edata.at(iline+1),tokens,"="); Crystal_Real_space_Bravais_Lattice_Primitive=tokens.at(tokens.size()-1);
		aurostd::string2tokens(vline_edata.at(iline+2),tokens,"="); Crystal_Real_space_Lattice_Variation=tokens.at(tokens.size()-1);
		aurostd::string2tokens(vline_edata.at(iline+3),tokens,"="); Crystal_Real_space_Lattice_System=tokens.at(tokens.size()-1);
		aurostd::string2tokens(vline_edata.at(iline+4),tokens,"="); Crystal_Real_space_Pearson_Symbol=tokens.at(tokens.size()-1);}
	      if(aurostd::substring2bool(vline_edata.at(iline),"POINT GROUP CRYSTAL")) {
		aurostd::string2tokens(vline_edata.at(iline+1),tokens,"="); Crystal_Real_space_Crystal_Family=tokens.at(tokens.size()-1);
		aurostd::string2tokens(vline_edata.at(iline+2),tokens,"="); Crystal_Real_space_Crystal_System=tokens.at(tokens.size()-1);
		aurostd::string2tokens(vline_edata.at(iline+3),tokens,"="); Crystal_Real_space_Crystal_Class=tokens.at(tokens.size()-1);
		aurostd::string2tokens(vline_edata.at(iline+4),tokens,"="); Crystal_Real_space_Point_Group_Hermann_Mauguin=tokens.at(tokens.size()-1);
		aurostd::string2tokens(vline_edata.at(iline+5),tokens,"="); Crystal_Real_space_Point_Group_Schoenflies=tokens.at(tokens.size()-1);
		aurostd::string2tokens(vline_edata.at(iline+6),tokens,"="); Crystal_Real_space_Point_Group_Orbifold=tokens.at(tokens.size()-1);
		aurostd::string2tokens(vline_edata.at(iline+7),tokens,"="); Crystal_Real_space_Point_Group_Type=tokens.at(tokens.size()-1);
		aurostd::string2tokens(vline_edata.at(iline+8),tokens,"="); Crystal_Real_space_Point_Group_Order=tokens.at(tokens.size()-1);
		aurostd::string2tokens(vline_edata.at(iline+9),tokens,"="); Crystal_Real_space_Point_Group_Structure=tokens.at(tokens.size()-1);}
	      if(aurostd::substring2bool(vline_edata.at(iline),"BRAVAIS LATTICE OF THE LATTICE")) {
		aurostd::string2tokens(vline_edata.at(iline+1),tokens,"="); Lattice_Real_space_Bravais_Lattice_Primitive=tokens.at(tokens.size()-1);
		aurostd::string2tokens(vline_edata.at(iline+2),tokens,"="); Lattice_Real_space_Lattice_Variation=tokens.at(tokens.size()-1);
		aurostd::string2tokens(vline_edata.at(iline+3),tokens,"="); Lattice_Real_space_Lattice_System=tokens.at(tokens.size()-1);}
	      if(aurostd::substring2bool(vline_edata.at(iline),"SUPERLATTICE")) {
		aurostd::string2tokens(vline_edata.at(iline+1),tokens,"="); Superattice_Real_space_Bravais_Superlattice_Primitive=tokens.at(tokens.size()-1);
		aurostd::string2tokens(vline_edata.at(iline+2),tokens,"="); Superattice_Real_space_Superlattice_Variation=tokens.at(tokens.size()-1);
		aurostd::string2tokens(vline_edata.at(iline+3),tokens,"="); Superattice_Real_space_Superlattice_System=tokens.at(tokens.size()-1);
		aurostd::string2tokens(vline_edata.at(iline+4),tokens,"="); Superattice_Real_space_Pearson_Symbol_Superlattice=tokens.at(tokens.size()-1);}
	      if(aurostd::substring2bool(vline_edata.at(iline),"Reciprocal space a b c alpha beta gamma")) {
		aurostd::string2tokens(vline_edata.at(iline),tokens);
		for(uint i=0;i<6;i++) abcK.at(i)=aurostd::string2utype<double>(tokens.at(tokens.size()-6+i));}
	      if(aurostd::substring2bool(vline_edata.at(iline),"Reciprocal space Volume")) {
		aurostd::string2tokens(vline_edata.at(iline),tokens);
		volumeK=aurostd::string2utype<double>(tokens.at(tokens.size()-1));}
	      if(aurostd::substring2bool(vline_edata.at(iline),"RECIPROCAL LATTICE")) {
		aurostd::string2tokens(vline_edata.at(iline+7),tokens,"="); Reciprocal_lattice_primitive=tokens.at(tokens.size()-1);
		aurostd::string2tokens(vline_edata.at(iline+8),tokens,"="); Reciprocal_lattice_variation=tokens.at(tokens.size()-1);}
	    }
	  }
	
	  if(aentry.vcomposition.size()==aentry.vspecies.size()) {
	    for(uint i=0;i<aentry.vspecies.size();i++) 
	      density=density+aentry.vcomposition.at(i)*GetAtomMass(aentry.vspecies.at(i))/volumeR*1000.0*1e8*1e8*1e8; //grams/cm^3
	  }
	
	  // Print out structural data    
	 // oss << "<ul>" << endl; //JPO 180731
	 // oss << "<li>" << endl; //JPO 180731
	  oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Lattice:</h5>" << endl; //JPO 180731
	 // oss << "<div class=\"lattice_table\">" << endl; //JPO 180731
	 // oss << "<table class=\"lattice\">" << endl; //JPO 180731
	 // oss << "<tbody>" << endl; //JPO 180731
	  oss << "<div class=\"value\">" << endl << "a=" << abcR.at(0) << "&Aring;" << "&nbsp;" << endl << " b=" << abcR.at(1) << "&Aring;" << "&nbsp;" <<  endl << " c=" << abcR.at(2) << "&Aring;" << "&nbsp;" << endl << " c/a=" << coveraR <<  endl << "</div>" << endl; //JPO 180731
	  oss << "<div class=\"value\">" << endl << "&alpha;=" << abcR.at(3) << "&deg" << endl << " &beta;=" << abcR.at(4) << "&deg" << endl << " &gamma;=" << abcR.at(5) << "&deg" << endl << "</div>" << endl; //JPO 180731
          //DX 20180824 - will replace lines above when the database is sufficiently populated - START
          //DX 20180824 vector<double> lattice_params(6);
          //DX 20180824 vector<string> rtokens;
          //DX 20180824 aurostd::string2tokens(aentry.geometry,rtokens,";");
          //DX 20180824 for(uint t=0;t<rtokens.size();t++){ lattice_params[t] = aurostd::string2utype<double>(rtokens[t]); }
          //DX 20180824 double covera = lattice_params.at(0)/lattice_params.at(2);
	  //DX 20180824 oss << "<div class=\"value\">" << endl << "a=" << lattice_params.at(0) << "&Aring;" << "&nbsp;" << endl << " b=" << lattice_params.at(1) << "&Aring;" << "&nbsp;" <<  endl << " c=" << lattice_params.at(2) << "&Aring;" << "&nbsp;" << endl << " c/a=" << covera <<  endl << "</div>" << endl; //JPO 180731
	  //DX 20180824 oss << "<div class=\"value\">" << endl << "&alpha;=" << lattice_params.at(3) << "&deg" << endl << " &beta;=" << lattice_params.at(4) << "&deg" << endl << " &gamma;=" << lattice_params.at(5) << "&deg" << endl << "</div>" << endl; //JPO 180731
          //DX 20180824 - will replace lines above when the database is sufficiently populated - END
	 // oss << "</tbody>" << endl; //JPO 180731
	 // oss << "</table>" << endl; //JPO 180731
	 // oss << "</div>" << endl; //JPO 180731
	 // oss << "</li>" << endl; //JPO 180731
    oss << "</div></div>" << endl; //JPO 180731
	  oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Volume:</h5><span class=\"value\">" << volumeR << "&Aring;<sup>3</sup></span></div></div>" << endl; //JPO 180731
          //DX 20180824 - will replace lines above when the database is sufficiently populated - START
	  //DX 20180824 oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Volume:</h5><span class=\"value\">" << aentry.volume_cell << "&Aring;<sup>3</sup></span></div></div>" << endl; //JPO 180731
          //DX 20180824 - will replace lines above when the database is sufficiently populated - END
	  oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Number of Atoms per Cell:</h5><span class=\"value\">" << aentry.natoms << "</span></div></div>" << endl; //JPO 180731
	  //  if(html && aentry.density==0.0) oss << "<li><span class=\"description\"> Density(calc):</span> " << density << " g/cm<sup>3</sup></li>" << endl;
	  // if(html && aentry.density>0.0) oss << "<li><span class=\"description\"> Density(entry):</span> " << aentry.density << " g/cm<sup>3</sup></li>" << endl;
	  if(aentry.density<0.1)  aentry.density=density;
	  oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Density(calc):</h5><span class=\"value\"> " << density << " g/cm<sup>3</sup></span></div></div>" << endl; //JPO 180731
	  oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Density(aflowlib):</h5><span class=\"value\"> " << aentry.density << " g/cm<sup>3</sup></span></div></div>" << endl; //JPO 180731

	  //	  if(aurostd::substring2bool(aentry.vfiles_WEB,"CONTCAR.relax")) {
	  //    oss << "<li><span class=\"description\"> Relaxed position (aflowlib/VASP):</span>" << "[<a href=\"" << url_WEB << "/CONTCAR.relax\">POSCAR</a>]" << "</li>" << endl;
	  // }
	  
	if (i==1){  //GG 170714 - removed relaxed data from orig
	  if(aurostd::substring2bool(aentry.vfiles_WEB,"CONTCAR.relax.vasp")) {
	    oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Relaxed position (aflowlib/VASP):</h5><span class=\"value\">" 
		<< "[<a href=\"" << url_WEB << "/CONTCAR.relax.vasp\"" << html_TAB << ">VASP-POSCAR</a>]" << "</span></div></div>" << endl; //JPO 180731
	  }
	}

	if (i==1){  //GG 170714 - removed relaxed data from orig
	  if(aurostd::substring2bool(aentry.vfiles_WEB,"CONTCAR.relax.qe")) {
	    oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Relaxed position (aflowlib/QE):</h5><span class=\"value\">" 
		<< "[<a href=\"" << url_WEB << "/CONTCAR.relax.qe\"" << html_TAB << ">QE-GEOMETRY</a>]" << "</span></div></div>" << endl; //JPO 180731
	  }
	}

	if (i==1){  //GG 170714 - removed relaxed data from orig
	  if(aurostd::substring2bool(aentry.vfiles_WEB,"CONTCAR.relax.abinit")) {
	    oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Relaxed position (aflowlib/ABINIT):</h5><span class=\"value\">" 
		<< "[<a href=\"" << url_WEB << "/CONTCAR.relax.abinit\"" << html_TAB << ">ABINIT-GEOMETRY</a>]" << "</span></div></div>" << endl; //JPO 180731
	  }
	}

	if (i==1){  //GG 170714 - removed relaxed data from orig
	  if(aurostd::substring2bool(aentry.vfiles_WEB,"CONTCAR.relax.aims")) {
	    oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Relaxed position (aflowlib/AIMS):</h5><span class=\"value\">" 
		<< "[<a href=\"" << url_WEB << "/CONTCAR.relax.aims\"" << html_TAB << ">AIMS-GEOMETRY</a>]" << "</span></div></div>" << endl; //JPO 180731
	  }
	}

	if (i==1){  //GG 170714 - removed relaxed data from orig
	  if(aurostd::substring2bool(aentry.vfiles_WEB,"INCAR.relax")) {
	    oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> INCAR for relax calculation (aflowlib/VASP):</h5><span class=\"value\">" 
		<< "[<a href=\"" << url_WEB << "/INCAR.relax\"" << html_TAB << ">INCAR.relax</a>]" << "</span></div></div>" << endl; //JPO 180731
	  }
	}

	  if(aurostd::substring2bool(aentry.vfiles_WEB,"INCAR.static")) {
	    oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> INCAR for static calculation (aflowlib/VASP):</h5><span class=\"value\">" 
		<< "[<a href=\"" << url_WEB << "/INCAR.static\"" << html_TAB << ">INCAR.static</a>]" << "</span></div></div>" << endl; //JPO 180731
	  }

	if (i==1){  //GG 170714 - removed relaxed data from orig
	  if(aurostd::substring2bool(aentry.vfiles_WEB,"INCAR.bands")) {
	    oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> INCAR for bands calculation (aflowlib/VASP):</h5><span class=\"value\">"
		<< "[<a href=\"" << url_WEB << "/INCAR.bands\"" << html_TAB << ">INCAR.bands</a>]" << "</span></div></div>" << endl; //JPO 180731
	  }
	}

	if (i==1){  //GG 170714 - removed relaxed data from orig
	  if(aurostd::substring2bool(aentry.vfiles_WEB,"KPOINTS.relax")) {
	    oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> KPOINTS for relax calculation (aflowlib/VASP):</h5><span class=\"value\">" 
		<< "[<a href=\"" << url_WEB << "/KPOINTS.relax\"" << html_TAB << ">KPOINTS.relax</a>]" << "</span></div></div>" << endl; //JPO 180731
	  }
	}
	
  if (i==1){  //GG 170714 - removed relaxed data from orig
	  if(aurostd::substring2bool(aentry.vfiles_WEB,"KPOINTS.static")) {
	    oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> KPOINTS for static calculation (aflowlib/VASP):</h5><span class=\"value\">"
		<< "[<a href=\"" << url_WEB << "/KPOINTS.static\"" << html_TAB << ">KPOINTS.static</a>]" << "</span></div></div>" << endl; //JPO 180731
	  }
	}
  
  if (i==1){  //GG 170714 - removed relaxed data from orig
	  if(aurostd::substring2bool(aentry.vfiles_WEB,"KPOINTS.bands")) {
	    oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> KPOINTS for bands calculation (aflowlib/VASP):</h5><span class=\"value\">" 
		<< "[<a href=\"" << url_WEB << "/KPOINTS.bands\"" << html_TAB << ">KPOINTS.bands</a>]" << "</span></div></div>" << endl; //JPO 180731
	  }
	}
	  if(aurostd::substring2bool(aentry.vfiles_WEB,DEFAULT_FILE_EDATA_ORIG_OUT)) {
	    oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Extended crystallographic data for original structure:</h5><span class=\"value\">" 
		<< "[<a href=\"" << url_WEB << "/" << DEFAULT_FILE_EDATA_ORIG_OUT << "\"" << html_TAB << ">" << DEFAULT_FILE_EDATA_ORIG_OUT << "</a>]" << "</span></div></div>" << endl; //JPO 180731
	  }

  if (i==1){  //GG 170714 - removed relaxed data from orig
	  if(aurostd::substring2bool(aentry.vfiles_WEB,DEFAULT_FILE_EDATA_RELAX_OUT)) {
	    oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Extended crystallographic data for relaxed structure:</h5><span class=\"value\">" 
		<< "[<a href=\"" << url_WEB << "/" << DEFAULT_FILE_EDATA_RELAX_OUT << "\"" << html_TAB << ">" << DEFAULT_FILE_EDATA_RELAX_OUT << "</a>]" << "</span></div></div>" << endl; //JPO 180731
	  }
	}

  if (i==1){  //GG 170714 - removed relaxed data from orig
	  if(aurostd::substring2bool(aentry.vfiles_WEB,DEFAULT_FILE_EDATA_BANDS_OUT)) {
	    oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Extended crystallographic data for band-structure:</h5><span class=\"value\">" 
		<< "[<a href=\"" << url_WEB << "/" << DEFAULT_FILE_EDATA_BANDS_OUT << "\"" << html_TAB << ">" << DEFAULT_FILE_EDATA_BANDS_OUT << "</a>]" << "</span></div></div>" << endl; //JPO 180731
	  }
	}

	 // oss << "</ul>" << endl; //JPO 180731
	 // oss << "</div>" << endl; //JPO 180731
	 // oss << "<hr />" << endl; //JPO 180731
	 // oss << "<div class=\"space_group\">" << endl; //JPO 180731
	  oss << "<div class=\"container-subtitle\"><h4 class=\"section-subtitle\"> Bravais Lattice of the Crystal" << aflow_sym_readme << art135_link << "</h4></div>" << endl; // PC 180620 //JPO 180731  //CO 180817
	 // oss << "<ul>" << endl; //JPO 180731
	  oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Space Group Number:</h5><span class=\"value\">" << aentry.sg2 << "</span></div></div>" << endl; //JPO 180731
	  oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Pearson Symbol:</h5><span class=\"value\">" << Crystal_Real_space_Pearson_Symbol << "</span></div></div>" << endl; //JPO 180731
	  oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Bravais Lattice Primitive:</h5><span class=\"value\">" << Crystal_Real_space_Bravais_Lattice_Primitive << "</span></div></div>" << endl; //JPO 180731
	  oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Bravais Lattice Variation:</h5><span class=\"value\">" << Crystal_Real_space_Lattice_Variation << "</span></div></div>" << endl; //JPO 180731
	  oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Bravais Lattice System:</h5><span class=\"value\">" << Crystal_Real_space_Lattice_System << "</span></div></div>" << endl; //JPO 180731
	 // oss << "</ul>" << endl; //JPO 180731
	 // oss << "</div>" << endl; //JPO 180731
	 // oss << "<hr />" << endl; //JPO 180731
	 // oss << "<div class=\"point_group\">" << endl; //JPO 180731
	  oss << "<div class=\"container-subtitle\"><h4 class=\"section-subtitle\"> Point Group of the Crystal" << aflow_sym_readme << art135_link << "</h4></div>" << endl; // PC 180620 //JPO 180731 //CO 180817
	 // oss << "<ul>" << endl; //JPO 180731
	  oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Crystal Family:</h5><span class=\"value\">" << Crystal_Real_space_Crystal_Family << "</span></div></div>" << endl; //JPO 180731
	  oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Crystal System:</h5><span class=\"value\">" << Crystal_Real_space_Crystal_System << "</span></div></div>" << endl; //JPO 180731
	  oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Crystal Class:</h5><span class=\"value\">" << Crystal_Real_space_Crystal_Class << "</span></div></div>" << endl; //JPO 180731
	  //	oss << "<li><span class=\"description\"> Point Group (Hermann Mauguin):</span>" << Crystal_Real_space_Point_Group_Hermann_Mauguin << "</li><!br>" << endl; //JPO 180731
	  oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Point Group (Herm. Maug.):</h5><span class=\"value\">" << Crystal_Real_space_Point_Group_Hermann_Mauguin << "</span></div></div>" << endl; //JPO 180731
	  oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Point Group (Schoenflies):</h5><span class=\"value\">" << Crystal_Real_space_Point_Group_Schoenflies << "</span></div></div>" << endl; //JPO 180731
	  oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Point Group Orbifold:</h5><span class=\"value\">" << Crystal_Real_space_Point_Group_Orbifold << "</span></div></div>" << endl; //JPO 180731
          //DX 20180827 - new point group type output - START
          string point_group_type = Crystal_Real_space_Point_Group_Type;
	  if(aurostd::RemoveWhiteSpaces(point_group_type) == "-" || aurostd::RemoveWhiteSpaces(point_group_type) == "none"){
            point_group_type = "non-centrosymmetric, non-enantiomorphic, non-polar";
          }
	  oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Point Group Type:</h5><span class=\"value\">" << point_group_type << "</span></div></div>" << endl; //JPO 180731
          //DX 20180827 - new point group type output - END
	  //DX 20180827 [OBSOLETE] oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Point Group Type:</h5><span class=\"value\">" << Crystal_Real_space_Point_Group_Type << "</span></div></div>" << endl; //JPO 180731
	  oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Point Group Type:</h5><span class=\"value\">" << Crystal_Real_space_Point_Group_Type << "</span></div></div>" << endl; //JPO 180731
	  oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Point Group Order:</h5><span class=\"value\">" << Crystal_Real_space_Point_Group_Order << "</span></div></div>" << endl; //JPO 180731
	  oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Point Group Structure:</h5><span class=\"value\">" << Crystal_Real_space_Point_Group_Structure << "</span></div></div>" << endl; //JPO 180731
          //DX 20180824 - will replace lines above when the database is sufficiently populated - START
	  //DX 20180824 oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Crystal Family:</h5><span class=\"value\">" << aentry.crystal_family << "</span></div></div>" << endl; //JPO 180731
	  //DX 20180824 oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Crystal System:</h5><span class=\"value\">" << aentry.crystal_system << "</span></div></div>" << endl; //JPO 180731
	  //DX 20180824 oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Crystal Class:</h5><span class=\"value\">" << aentry.crystal_class << "</span></div></div>" << endl; //JPO 180731
	  //DX 20180824 oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Point Group (Herm. Maug.):</h5><span class=\"value\">" << aentry.point_group_Hermann_Mauguin << "</span></div></div>" << endl; //JPO 180731
	  //DX 20180824 oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Point Group (Schoenflies):</h5><span class=\"value\">" << aentry.point_group_Schoenflies << "</span></div></div>" << endl; //JPO 180731
	  //DX 20180824 oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Point Group Orbifold:</h5><span class=\"value\">" << aentry.point_group_orbifold << "</span></div></div>" << endl; //JPO 180731
	  //DX 20180824 oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Point Group Type:</h5><span class=\"value\">" << aentry.point_group_type << "</span></div></div>" << endl; //JPO 180731
	  //DX 20180824 oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Point Group Order:</h5><span class=\"value\">" << aentry.point_group_order << "</span></div></div>" << endl; //JPO 180731
	  //DX 20180824 oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Point Group Structure:</h5><span class=\"value\">" << aentry.point_group_structure << "</span></div></div>" << endl; //JPO 180731
          //DX 20180824 - will replace lines above when the database is sufficiently populated - END
	 // oss << "</ul>" << endl; //JPO 180731
	 // oss << "</div>" << endl; //JPO 180731
	 // oss << "<hr />" << endl; //JPO 180731
	 // oss << "<div class=\"bravais_lattice\">" << endl; //JPO 180731
	  oss << "<div class=\"container-subtitle\"><h4 class=\"section-subtitle\"> Bravais Lattice of the Lattice" << aflow_sym_readme << art135_link << "</h4></div>" << endl; // PC 180620 //JPO 180731 //CO 180817
	 // oss << "<ul>" << endl; //JPO 180731
	  oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Bravais Lattice Primitive</h5><span class=\"value\">" << Lattice_Real_space_Bravais_Lattice_Primitive << "</span></div>" << endl; //JPO 180731
	  oss << "<div class=\"container__card\"><h5 class=\"value-name\"> Bravais Lattice Variation:</h5><span class=\"value\">" << Lattice_Real_space_Lattice_Variation << "</span></div>" << endl; //JPO 180731
	  oss << "<div class=\"container__card\"><h5 class=\"value-name\"> Bravais Lattice System:</h5><span class=\"value\">" << Lattice_Real_space_Lattice_System << "</span></div>" << endl; //JPO 180731
          //DX 20180824 - will replace lines above when the database is sufficiently populated - START
	  //DX 20180824 oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Bravais Lattice Primitive</h5><span class=\"value\">" << aentry.Bravais_lattice_lattice_type << "</span></div>" << endl; //JPO 180731
	  //DX 20180824 oss << "<div class=\"container__card\"><h5 class=\"value-name\"> Bravais Lattice Variation:</h5><span class=\"value\">" << aentry.Bravais_lattice_lattice_variation_type << "</span></div>" << endl; //JPO 180731
	  //DX 20180824 oss << "<div class=\"container__card\"><h5 class=\"value-name\"> Bravais Lattice System:</h5><span class=\"value\">" << aentry.Bravais_lattice_lattice_system << "</span></div>" << endl; //JPO 180731
          //DX 20180824 - will replace lines above when the database is sufficiently populated - END
	 // oss << "</ul>" << endl; //JPO 180731
          oss << "</div>" << endl; //JPO 180731
	  if(aurostd::substring2bool(aentry.vfiles_WEB,label+"_BZ.png")) {
	    oss << "<div class=\"container__cell\"><div class=\"container__card--img\">" << endl; //JPO 180731
	    oss << "<h5 class=\"value-name\"> Brillouin Zone " << art058_link<< "</h5>" << endl; //JPO 180731
	   // oss << "</div>" << endl; //JPO 180731
	  // oss << "<div class=\"picture_BZ\">" << endl; //JPO 180731
	    //oss << "<img class=\"pic_BZ\" src=\"../SCIENCE/images/brillouin/" << aurostd::RemoveWhiteSpaces(Lattice_Real_space_Lattice_Variation) << ".PNG\" alt=\"Brillouin Zone of " << label << "\" />" << endl; // CO 170621 - relative path
	    oss << "<img class=\"BZ-img\" src=\"http://aflowlib.duke.edu/SCIENCE/images/brillouin/" << aurostd::RemoveWhiteSpaces(Lattice_Real_space_Lattice_Variation) << ".PNG\" alt=\"Brillouin Zone of " << label << "\" />" << endl; // CO 170621 - abs path WORKS  //JPO 180731
	    oss << "</div></div>" << endl; //JPO 180731
	  }
	  //[MOVED UP JPO 180731]// oss << "<hr />" << endl; //JPO 180731
	  //[MOVED UP JPO 180731]// oss << "<div class=\"bravais_lattice\">" << endl; //JPO 180731
	  //[MOVED UP JPO 180731]oss << "<div class=\"container-subtitle\"><h4 class=\"section-subtitle\"> Bravais Lattice of the Lattice" << aflow_sym_readme << art135_link << "</h4></div>" << endl; // PC 180620 //JPO 180731 //CO 180817
	  //[MOVED UP JPO 180731]// oss << "<ul>" << endl; //JPO 180731
	  //[MOVED UP JPO 180731]oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Bravais Lattice Primitive</h5><span class=\"value\">" << Lattice_Real_space_Bravais_Lattice_Primitive << "</span></div>" << endl; //JPO 180731
	  //[MOVED UP JPO 180731]oss << "<div class=\"container__card\"><h5 class=\"value-name\"> Bravais Lattice Variation:</h5><span class=\"value\">" << Lattice_Real_space_Lattice_Variation << "</span></div>" << endl; //JPO 180731
	  //[MOVED UP JPO 180731]oss << "<div class=\"container__card\"><h5 class=\"value-name\"> Bravais Lattice System:</h5><span class=\"value\">" << Lattice_Real_space_Lattice_System << "</span></div>" << endl; //JPO 180731
	  //[MOVED UP JPO 180731]// oss << "</ul>" << endl; //JPO 180731
          //[MOVED UP JPO 180731]oss << "</div>" << endl; //JPO 180731
	 // oss << "<hr />" << endl; //JPO 180731
	 // oss << "<div class=\"superlattice\">" << endl; //JPO 180731
	  oss << "<div class=\"container-subtitle\"><h4 class=\"section-subtitle\"> Superlattice" << aflow_sym_readme << art135_link << "</h4></div>" << endl; // PC 180620 //JPO 180731 //CO 180817
	 // oss << "<ul>" << endl; //JPO 180731
	  oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Superlattice Primitive unit cell:</h5><span class=\"value\">" << Superattice_Real_space_Bravais_Superlattice_Primitive << "</span></div></div>" << endl; //JPO 180731
	  oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Superlattice Variation:</h5><span class=\"value\">" << Superattice_Real_space_Superlattice_Variation << "</span></div></div>" << endl; //JPO 180731
	  oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Superlattice Lattice System :</h5><span class=\"value\">" << Superattice_Real_space_Superlattice_System << "</span></div></div>" << endl; //JPO 180731
	  oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Superlattice Pearson Symbol:</h5><span class=\"value\">" << Superattice_Real_space_Pearson_Symbol_Superlattice << "</span></div></div>" << endl; //JPO 180731
          //DX 20180824 - will replace lines above when the database is sufficiently populated - START
          //DX 20180824 oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Superlattice Primitive unit cell:</h5><span class=\"value\">" << aentry.Bravais_superlattice_lattice_type << "</span></div></div>" << endl; //JPO 180731
	  //DX 20180824 oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Superlattice Variation:</h5><span class=\"value\">" << aentry.Bravais_superlattice_lattice_variation_type << "</span></div></div>" << endl; //JPO 180731
	  //DX 20180824 oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Superlattice Lattice System :</h5><span class=\"value\">" << aentry.Bravais_superlattice_lattice_system << "</span></div></div>" << endl; //JPO 180731
	  //DX 20180824 oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Superlattice Pearson Symbol:</h5><span class=\"value\">" << aentry.Pearson_symbol_superlattice << "</span></div></div>" << endl; //JPO 180731
          //DX 20180824 - will replace lines above when the database is sufficiently populated - END
	 // oss << "</ul>" << endl; //JPO 180731
	 // oss << "</div>" << endl; //JPO 180731
	 // oss << "<hr />" << endl; //JPO 180731
	 // oss << "<div class=\"reciprocal\">" << endl; //JPO 180731
	  oss << "<div class=\"container-subtitle\"><h4 class=\"section-subtitle\"> Reciprocal Space Lattice </h4></div>" << endl; //JPO 180731
	 // oss << "<ul>" << endl; //JPO 180731
	  oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Reciprocal Lattices:</h5>" << endl; //JPO 180731
	 // oss << "<div class=\"lattice_table\">" << endl; //JPO 180731
	 // oss << "<table class=\"reciprocal_lattice\">" << endl; //JPO 180731
	 // oss << "<tbody>" << endl; //JPO 180731
	  oss << "<div class=\"value\">" << endl << "a=" << abcK.at(0) << "&Aring;<sup>-1</sup>" << endl << "b=" << abcK.at(1) << "&Aring;<sup>-1</sup>" << endl << "c=" << abcK.at(2) << "&Aring;<sup>-1</sup>" << endl << "</div>" << endl; //JPO 180731
	 // oss << "<tr>" << endl; //JPO 180731
	  oss << "<div class=\"value\">" << "&alpha;=" << abcK.at(3) << "&deg;" << endl << "&beta;=" << abcK.at(4) << "&deg;" << endl << "&gamma;=" << abcK.at(5) << "&deg;" << "</div>" << endl; //JPO 180731
          //DX 20180824 - will replace lines above when the database is sufficiently populated - START
          //DX 20180824 vector<double> reciprocal_lattice_params(6);
          //DX 20180824 vector<string> ktokens;
          //DX 20180824 aurostd::string2tokens(aentry.reciprocal_geometry,ktokens,";");
          //DX 20180824 for(uint t=0;t<ktokens.size();t++){ reciprocal_lattice_params[t] = aurostd::string2utype<double>(ktokens[t]); }
	  //DX 20180824 oss << "<div class=\"value\">" << endl << "a=" << reciprocal_lattice_params.at(0) << "&Aring;<sup>-1</sup>" << endl << "b=" << reciprocal_lattice_params.at(1) << "&Aring;<sup>-1</sup>" << endl << "c=" << reciprocal_lattice_params.at(2) << "&Aring;<sup>-1</sup>" << endl << "</div>" << endl; //JPO 180731
	  //DX 20180824 oss << "<div class=\"value\">" << "&alpha;=" << reciprocal_lattice_params.at(3) << "&deg;" << endl << "&beta;=" << reciprocal_lattice_params.at(4) << "&deg;" << endl << "&gamma;=" << reciprocal_lattice_params.at(5) << "&deg;" << "</div>" << endl; //JPO 180731
          //DX 20180824 - will replace lines above when the database is sufficiently populated - END
	 // oss << "<tr>" << endl; //JPO 180731
	 // oss << "</tbody>" << endl; //JPO 180731
	 // oss << "</table>" << endl; //JPO 180731
	  oss << "</div></div>" << endl; //JPO 180731
	  oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Volume:</h5><span class=\"value\">" << volumeK << " &Aring;<sup>-3</sup></span></div></div>" << endl; //JPO 180731
	  oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Lattice Primitive:</h5><span class=\"value\">" << Reciprocal_lattice_primitive << "</span></div></div>" << endl; //JPO 180731
	  oss << "<div class=\"container__cell\"><div class =\"container__card\"><h5 class=\"value-name\"> Lattice Variation:</h5><span class=\"value\">" << Reciprocal_lattice_variation << "</span></div></div>" << endl; //JPO 180731
          //DX 20180824 - will replace lines above when the database is sufficiently populated - START
	  //DX 20180824 oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Volume:</h5><span class=\"value\">" << aentry.reciprocal_volume_cell << " &Aring;<sup>-3</sup></span></div></div>" << endl; //JPO 180731
	  //DX 20180824 oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Lattice Primitive:</h5><span class=\"value\">" << aentry.reciprocal_lattice_type << "</span></div></div>" << endl; //JPO 180731
	  //DX 20180824 oss << "<div class=\"container__cell\"><div class =\"container__card\"><h5 class=\"value-name\"> Lattice Variation:</h5><span class=\"value\">" << aentry.reciprocal_lattice_variation_type << "</span></div></div>" << endl; //JPO 180731
          //DX 20180824 - will replace lines above when the database is sufficiently populated - END
	 // oss << "</ul>" << endl; //JPO 180731
	  oss << "</div>" << endl;  //JPO 180731
	  oss << "<!-- EDATA: END -->" << endl;
	}
      }
    }
    // ***************************************************************************
    // THERMODYNAMICS 
    if(vflags.flag("FLAG::FOUND") && vflags.flag("FLAG::THERMODYNAMICS") && !directory.empty()) {
      oss << "<!-- Thermodynamics properties: BEGIN -->" << endl;
     // oss << line_rule << endl; //JPO 180731
     // oss << "<div class=\"Thermodynamics\">" << endl; //JPO 180731
      oss << "<div class=\"container\">" << endl; //JPO 180731
      oss << "<div class=\"container-title\"><h1 class=\"section-title\">Thermodynamics Properties</h1></div>" << endl; //JPO 180731
      //   oss << line_rule << endl; //JPO 180731
      //     oss << "<br>" << endl; //JPO 180731
     // oss << "<span class=\"Thermodynamics_table\">" << endl; //JPO 180731
     // oss << "<ul>" << endl; //JPO 180731
      oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\">Formation Enthalpy_cell:</h5><span class=\"value\">" << (aentry.enthalpy_formation_cell!=0.0?aurostd::utype2string(aentry.enthalpy_formation_cell,5)+" (eV) ":"unavailable") << "</span></div></div>" << endl; //JPO 180731
      oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\">Formation Enthalpy_atom:</h5><span class=\"value\">" << (aentry.enthalpy_formation_atom!=0.0?aurostd::utype2string(aentry.enthalpy_formation_atom,5)+" (eV/atom) ":"unavailable") << "</span></div></div>" << endl; //JPO 180731
      oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"><i>Ab-initio</i> energy_cell:</h5><span class=\"value\">" << (aentry.energy_cell!=0.0?aurostd::utype2string(aentry.energy_cell,5)+" (eV) ":"unavailable") << "</span></div></div>" << endl; //JPO 180731
      oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"><i>Ab-initio</i> energy_atom:</h5><span class=\"value\">" << (aentry.energy_cell!=0.0?aurostd::utype2string(aentry.energy_atom,5)+" (eV/atom) ":"unavailable") << "</span></div></div>" << endl; //JPO 180731
     // oss << "</ul>" << endl; //JPO 180731
      oss << "</div>" << endl; //JPO 180731
      oss << "<!-- Thermodynamics properties: END -->" << endl;
    }
    // ***************************************************************************
    // AGL
    if(vflags.flag("FLAG::FOUND") && vflags.flag("FLAG::AGL") && !directory.empty()) {
      oss << "<!-- AGL properties: BEGIN -->" << endl;
     // oss << line_rule << endl; //JPO 180731
     // oss << "<div class=\"AGL\">" << endl; //JPO 180731
      oss << "<div class=\"container\">" << endl; //JPO 180731
      oss << "<div class=\"container-title\"><h1 class=\"section-title\">AGL Properties (Aflow Gibbs Library)" << aflow_agl_readme << art115_link << "</h1></div>" << endl; // PC 180620 //JPO 180731  //CO 180817
      //   oss << line_rule << endl; //JPO 180731
      //     oss << "<br>" << endl; //JPO 180731
     // oss << "<span class=\"AGL_table\">" << endl; //JPO 180731
     // oss << "<ul>" << endl; //JPO 180731
      if(aurostd::substring2bool(aentry.vfiles_WEB,"aflow.agl.out")) 
	oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> AGL Output" << art096_link << ":</h5><span class=\"value\">" << "[<a href=\"" << url_WEB << "/aflow.agl.out\"" << html_TAB << ">aflow.agl.out</a>]" << "</span></div></div>" << endl; //JPO 180731
      if(aurostd::substring2bool(aentry.vfiles_WEB,"AGL.out")) 
	oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> AGL Complete Output" << art096_link << ":</h5><span class=\"value\">" << "[<a href=\"" << url_WEB << "/AGL.out\"" << html_TAB << ">AGL.out</a>]" << "</span></div></div>" << endl; //JPO 180731
      if(aurostd::substring2bool(aentry.vfiles_WEB,"AGL_energies_temperature.out")) 
	oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> AGL Energy versus Temperature" << art096_link << ":</h5><span class=\"value\">" << "[<a href=\"" << url_WEB << "/AGL_energies_temperature.out\"" << html_TAB << ">AGL_energies_temperature.out</a>]" << "</span></div></div>" << endl; //JPO 180731
      if(aurostd::substring2bool(aentry.vfiles_WEB,"AGL_thermal_properties_temperature.out")) 
	oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> AGL Thermal Properties versus Temperature" << art096_link << ":</h5><span class=\"value\">" << "[<a href=\"" << url_WEB << "/AGL_thermal_properties_temperature.out\"" << html_TAB << ">AGL_thermal_properties_temperature.out</a>]" << "</span></div></div>" << endl; //JPO 180731
      if(aentry.agl_thermal_conductivity_300K<AUROSTD_NAN) 
	oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\" id=\"agl_thermal_conductivity_300K_td\">AGL Thermal Conductivity at 300K" << art096_link << ":</h5><span class=\"value\">" << aentry.agl_thermal_conductivity_300K << " (W/m*K)</span></div></div>" << endl; //JPO 180731
      if(aentry.agl_debye<AUROSTD_NAN) 
	oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\" id=\"agl_debye_td\">AGL Debye Temperature" << art096_link << ":</h5><span class=\"value\">" << aentry.agl_debye << " (K)</span></div></div>" << endl; //JPO 180731
      if(aentry.agl_acoustic_debye<AUROSTD_NAN) 
	oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\" id=\"agl_acoustic_debye_td\">AGL Debye Acoustic Temperature" << art096_link << ":</h5><span class=\"value\">" << aentry.agl_acoustic_debye << " (K)</span></div></div>" << endl; //JPO 180731
      if(aentry.agl_gruneisen<AUROSTD_NAN) 
	oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\" id=\"agl_gruneisen_td\">AGL Gruneisen parameter" << art096_link << ":</h5><span class=\"value\">" << aentry.agl_gruneisen << "</span></div></div>" << endl; //JPO 180731
      if(aentry.agl_heat_capacity_Cv_300K<AUROSTD_NAN) 
	oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\" id=\"agl_heat_capacity_Cv_300K_td\">AGL Specific Heat Cv at 300K" << art096_link << ":</h5><span class=\"value\">" << aentry.agl_heat_capacity_Cv_300K << " (kB/cell)</span></div></div>" << endl; //JPO 180731
      if(aentry.agl_heat_capacity_Cp_300K<AUROSTD_NAN) 
	oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\" id=\"agl_heat_capacity_Cp_300K_td\">AGL Specific Heat Cp at 300K" << art096_link << ":</h5><span class=\"value\">" << aentry.agl_heat_capacity_Cp_300K << " (kB/cell)</span></div></div>" << endl; //JPO 180731
      if(aentry.agl_thermal_expansion_300K<AUROSTD_NAN) 
	oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\" id=\"agl_thermal_expansion_300K_td\">AGL Thermal Expansion at 300K" << art096_link << ":</h5><span class=\"value\">" << aentry.agl_thermal_expansion_300K << " (1/K)</span></div></div>" << endl; //JPO 180731
      if(aentry.agl_bulk_modulus_static_300K<AUROSTD_NAN) 
	oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\" id=\"agl_bulk_modulus_static_300K_td\">AGL Bulk Modulus Static at 300K" << art096_link << ":</h5><span class=\"value\">" << aentry.agl_bulk_modulus_static_300K << " (GPa)</span></div></div>" << endl; //JPO 180731
      if(aentry.agl_bulk_modulus_isothermal_300K<AUROSTD_NAN) 
	oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\" id=\"agl_bulk_modulus_isothermal_300K_td\">AGL Bulk Modulus Isothermal at 300K" << art096_link << ":</h5><span class=\"value\">" << aentry.agl_bulk_modulus_isothermal_300K << " (GPa)</span></div></div>" << endl; //JPO 180731
     // oss << "</ul>" << endl; //JPO 180731
     oss << "</div>" << endl; //JPO 180731
      oss << "<!-- AGL properties: END -->" << endl;
    }
    // ***************************************************************************
    // AEL
    if(vflags.flag("FLAG::FOUND") && vflags.flag("FLAG::AEL") && !directory.empty()) {
      oss << "<!-- AEL properties: BEGIN -->" << endl;
      //oss << line_rule << endl; //JPO 180731
     // oss << "<div class=\"AEL\">" << endl; //JPO 180731
      oss << "<div class=\"container\">" << endl; //JPO 180731
      oss << "<div class=\"container-title\"><h1 class=\"section-title\">AEL Properties (Aflow Elastic Library)</font>" << aflow_ael_readme << art096_link << "</h1></div>" << endl; // PC 180620 //JPO 180731  //CO 180817
     // oss << "<span class=\"AEL_table\">" << endl; //JPO 180731
     // oss << "<ul>" << endl; //JPO 180731
      if(aurostd::substring2bool(aentry.vfiles_WEB,"aflow.ael.out")) 
	oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> AEL Output" << art100_link << ":</h5><span class=\"value\">" << "[<a href=\"" << url_WEB << "/aflow.ael.out\"" << html_TAB << ">aflow.ael.out</a>]" << "</span></div></div>" << endl; //JPO 180731
      if(aurostd::substring2bool(aentry.vfiles_WEB,"AEL_Elastic_constants.out")) 
	oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> AEL Elastic Constants (stiffness tensor)" << art100_link << ":</h5><span class=\"value\">" << "[<a href=\"" << url_WEB << "/AEL_Elastic_constants.out\"" << html_TAB << ">AEL_Elastic_constants.out</a>]" << "</span></div></div>" << endl; //JPO 180731
      if(aurostd::substring2bool(aentry.vfiles_WEB,"AEL_Compliance_tensor.out")) 
	oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> AEL Compliance Constants (compliance tensor)" << art100_link << ":</h5><span class=\"value\">" << "[<a href=\"" << url_WEB << "/AEL_Compliance_tensor.out\"" << html_TAB << ">AEL_Compliance_tensor.out</a>]" << "</span></div></div>" << endl; //JPO 180731
      if(aentry.ael_poisson_ratio<AUROSTD_NAN) 
	oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\" id=\"ael_poisson_ratio_td\">AEL Poisson Ratio" << art100_link << ":</h5><span class=\"value\"> " << aentry.ael_poisson_ratio << "</span></div></div>" << endl; //JPO 180731
      if(aentry.ael_bulk_modulus_voigt<AUROSTD_NAN) 
	oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\" id=\"ael_bulk_modulus_voigt_td\">AEL Bulk Modulus Voigt" << art100_link << ":</h5><span class=\"value\">" << aentry.ael_bulk_modulus_voigt << " (GPa)</span></div></div>" << endl; //JPO 180731
      if(aentry.ael_bulk_modulus_reuss<AUROSTD_NAN) 
	oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\" id=\"ael_bulk_modulus_reuss_td\">AEL Bulk Modulus Reuss" << art100_link << ":</h5><span class=\"value\">" << aentry.ael_bulk_modulus_reuss << " (GPa)</span></div></div>" << endl; //JPO 180731
      if(aentry.ael_bulk_modulus_vrh<AUROSTD_NAN) 
	oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\" id=\"ael_bulk_modulus_vrh_td\">AEL Bulk Modulus Vrh" << art100_link << ":</h5><span class=\"value\">" << aentry.ael_bulk_modulus_vrh << " (GPa)</span></div></div>" << endl; //JPO 180731
      if(aentry.ael_shear_modulus_voigt<AUROSTD_NAN) 
	oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\" id=\"ael_shear_modulus_voigt_td\">AEL Shear Modulus Voigt" << art100_link << ":</h5><span class=\"value\">" << aentry.ael_shear_modulus_voigt << " (GPa)</span></div></div>" << endl; //JPO 180731
      if(aentry.ael_shear_modulus_reuss<AUROSTD_NAN) 
	oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\" id=\"ael_shear_modulus_reuss_td\">AEL Shear Modulus Reuss" << art100_link << ":</h5><span class=\"value\">" << aentry.ael_shear_modulus_reuss << " (GPa)</span></div></div>" << endl; //JPO 180731
      if(aentry.ael_shear_modulus_vrh<AUROSTD_NAN) 
	oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\" id=\"ael_shear_modulus_vrh_td\">AEL Shear Modulus Vrh" << art100_link << ":</h5><span class=\"value\">" << aentry.ael_shear_modulus_vrh << " (GPa)</span></div></div>" << endl; //JPO 180731
      if(aentry.ael_elastic_anistropy<AUROSTD_NAN) 
	oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\" id=\"ael_elastic_anistropy_td\">AEL Elastic Anisotropy" << art100_link << ":</h5><span class=\"value\">" << aentry.ael_elastic_anistropy << "</span></div></div>" << endl; //JPO 180731
     // oss << "</ul>" << endl; //JPO 180731
      oss << "</div>" << endl; //JPO 180731
      oss << "<!-- AEL properties: END -->" << endl;
    }
    // BADER
    if(vflags.flag("FLAG::FOUND") && vflags.flag("FLAG::BADER") && !directory.empty()) {
      oss << "<!-- Bader properties: BEGIN -->" << endl;
     // oss << line_rule << endl; //JPO 180731
     // oss << "<div class=\"Bader\">" << endl; //JPO 180731
      oss << "<div class=\"container\">" << endl; //JPO 180731
      oss << "<div class=\"container-title\"><h1 class=\"section-title\">Bader Atoms in Molecules Properties</font>" << "</h1></div>" << endl; //JPO 180731
      //oss << "<span class=\"title\"><FONT SIZE=+3>Bader Atoms in Molecules Properties</font>" << "<a href=http://materials.duke.edu/AFLOW/README_AFLOW_BADER.TXT>[info]</a>" << "</span>" << endl;
     // oss << "<span class=\"Bader_table\">" << endl; //JPO 180731
     // oss << "<ul>" << endl; //JPO 180731
      string abader_out=aentry.prototype+"_abader.out";
      if(aurostd::substring2bool(aentry.vfiles_WEB,abader_out)) {oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\"> Bader Output" << ":</h5><span class=\"value\">" << "[<a href=\"" << url_WEB << "/" << abader_out << "\"" << html_TAB << ">" << abader_out << "</a>]" << "</span></div></div>" << endl;} //JPO 180731
      //oss << "<li><span class=\"description\"> Bader Output" << art100_link << ":</span>" << "[<a href=\"" << url_WEB << "/" << abader_out << "\"" << html_TAB << ">" << abader_out << "</a>]" << "</li><!br>" << endl;
      string bader_net_charges_wiki_link=" [<a href=http://aflowlib.duke.edu/aflowwiki/doku.php?id=documentation:all_keywords&#bader_net_charges target=\"_blank\"><font color=black><i>info</i></font></a>]";
      string bader_atomic_volumes_wiki_link=" [<a href=http://aflowlib.duke.edu/aflowwiki/doku.php?id=documentation:all_keywords&#bader_atomic_volumes target=\"_blank\"><font color=black><i>info</i></font></a>]";
      if(!aentry.vbader_net_charges.empty()) {
        aurostd::string2tokens(aentry.species,tokens,",");  //tokens=species
        oss << "<div class=\"container__cell--full\"><div class=\"container__card\"><h5 class=\"value-name\" id=\"bader_net_charges_td\">Net Charges" << bader_net_charges_wiki_link << ":</h5><span class=\"value\"> "; //JPO 180731
        atomCOUNT=0;
        for(uint i=0;i<tokens.size();i++) {
          oss << tokens.at(i) << "={";
          for(uint j=0;j<(uint)aentry.vcomposition.at(i);j++) {
            num_prec.str("");
            num_prec << std::fixed << setprecision(4) << aentry.vbader_net_charges.at(atomCOUNT++);
            oss << num_prec.str();
            if(j<(uint)aentry.vcomposition.at(i)-1) {oss << ", ";}
          }
          if(i<tokens.size()-1) {oss << "}; ";} else {oss << "}";}
        }
        oss << " (electrons)</span></div></div>" << endl; //JPO 180731
      }
      if(!aentry.bader_atomic_volumes.empty()) {
        aurostd::string2tokens(aentry.species,tokens,",");  //tokens=species
	oss << "<div class=\"container__cell--full\"><div class=\"container__card\"><h5 class=\"value-name\" id=\"bader_atomic_volumes_td\">Atomic Volumes" << bader_atomic_volumes_wiki_link << ":</h5><span class=\"value\">"; //JPO 180731
        atomCOUNT=0;
        for(uint i=0;i<tokens.size();i++) {
          oss << tokens.at(i) << "={";
          for(uint j=0;j<(uint)aentry.vcomposition.at(i);j++) {
            num_prec.str("");
            num_prec << std::fixed << setprecision(4) << aentry.vbader_atomic_volumes.at(atomCOUNT++);
            oss << num_prec.str();
            if(j<(uint)aentry.vcomposition.at(i)-1) {oss << ", ";}
          }
          if(i<tokens.size()-1) {oss << "}; ";} else {oss << "}";}
        }
        oss << " (Angst<sup>3</sup>)</span></div></div>" << endl; //JPO 180731
      }
     // oss << "</ul>" << endl; //JPO 180731
      oss << "</div>" << endl; //JPO 180731
      oss << "<!-- Bader properties: END -->" << endl;
    }
    // ***************************************************************************
    // MAGNETIC
    if(vflags.flag("FLAG::FOUND") && vflags.flag("FLAG::MAGNETIC") && !directory.empty()) {
      oss << "<!-- Magnetic properties: BEGIN -->" << endl;
     // oss << line_rule << endl; //JPO 180731
     // oss << "<div class=\"Magnetic\">" << endl; //JPO 180731
      oss << "<div class=\"container\">" << endl; //JPO 180731
      oss << "<div class=\"container-title\"><h1 class=\"section-title\">Magnetic Properties </h1></div>" << endl; //JPO 180731
     // oss << "<!span class=\"Magnetic_table\">" << endl; //JPO 180731
     // oss << "<ul>" << endl; //JPO 180731
      // oss << "<li><span class=\"description\" id=\"spin_cell_td\">Magnetic Moment:</span> " << aentry.spin_cell << " &mu;<sub>B</sub></li><!br>" << endl; //JPO 180731
      // oss << "<li><span class=\"description\" id=\"spin_atom_td\">Magnetic Moment/atom:</span> " << aentry.spin_atom << " &mu;<sub>B</sub>/atom</li><!br>" << endl; //JPO 180731
      if(aurostd::substring2bool(aentry.loop,"magnetic")) {
	    oss << "<div class=\"container__cell--full\"><div class=\"container__card\"><h5 class=\"value-name\" id=\"spinD_td\">Spin Decomposition per atoms:</h5><span class=\"value\"> {" << aentry.spinD << "} &mu;<sub>B</sub></span></div></div>" << endl; //JPO 180731
	    oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\" id=\"spinF_td\">Spin Polarization (E<sub>F</sub>):</h5><span class=\"value\">" << aentry.spinF << "</span></div></div>" << endl; //JPO 180731
      }
      oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\" id=\"spin_cell_td\">Magnetic Moment:</h5><span class=\"value\">" << aentry.spin_cell << " &mu;<sub>B</sub></span></div></div>" << endl; //JPO 180731
      oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\" id=\"spin_atom_td\">Magnetic Moment/atom:</h5><span class=\"value\">" << aentry.spin_atom << " &mu;<sub>B</sub>/atom</span></div></div>" << endl; //JPO 180731
     // oss << "</ul>" << endl; //JPO 180731
      oss << "</div>" << endl; //JPO 180731
      oss << "<!-- Magnetic properties: END -->" << endl;
    }
    // ***************************************************************************
    // SCINTILLATION
    if(vflags.flag("FLAG::FOUND") && vflags.flag("FLAG::SCINTILLATION") && !directory.empty()) {
      double scintillation_attenuation_length=GetCompoundAttenuationLenght(aentry.vspecies,aentry.vcomposition,aentry.density);
     // oss << line_rule << endl; //JPO 180731
      oss << "<!-- Scintillation properties: BEGIN -->" << endl;
     // oss << "<div class=\"scintillation\">" << endl; //JPO 180731
      //    oss << "<hr />" << endl; //JPO 180731
      oss << "<div class=\"container\">" << endl; //JPO 180731
      oss << "<div class=\"container-title\"><h1 class=\"section-title\">Scintillation Properties</h1></div>" << endl; //JPO 180731
      //   oss << line_rule << endl; //JPO 180731
     // oss << "<br>" << endl; //JPO 180731
     // oss << "<table class=\"scintillation_table\">" << endl; //JPO 180731
     // oss << "<tbody>" << endl; //JPO 180731
      if(aentry.scintillation_attenuation_length>0.0) {
	oss << "<div class=\"container__cell\"><div class=\"container__card\">" << endl << "<h5 class=\"value-name\">Attenuation Length" << art064_link << ":</h5><span class=\"value\">" << aentry.scintillation_attenuation_length << " cm" << endl << "</span></div></div>" << endl; //JPO 180731
      } else {
	oss << "<div class=\"container__cell\"><div class=\"container__card\">" << endl << "<h5 class=\"value-name\">Attenuation Length" << art064_link << ":</h5><span class=\"value\">" << scintillation_attenuation_length << " cm" << endl << "</span></div></div>" << endl; //JPO 180731
      }
     // oss << "</tbody>" << endl; //JPO 180731
     // oss << "</table>" << endl; //JPO 180731
      oss << "</div>" << endl; //JPO 180731
      oss << "<!-- Scintillation properties: END -->" << endl;
    }
    // ***************************************************************************
    // ELECTRONIC BANDS
    if(vflags.flag("FLAG::FOUND") && vflags.flag("FLAG::ELECTRONIC") && !directory.empty()) {
     // oss << line_rule << endl; //JPO 180731
      oss << "<!-- Electronic properties: BEGIN -->" << endl;
      //   oss << "<hr />" << endl; //JPO 180731
     // oss << "<div class=\"electronic\">" << endl; //JPO 180731
      oss << "<div class=\"container\">" << endl; //JPO 180731
      oss << "<div class=\"container-title\"><h1 class=\"section-title\"> Electronic Properties </h1></div>" << endl; //JPO 180731
     // oss << "<table class=\"electronic_table\">" << endl; //JPO 180731
     // oss << "<tbody>" << endl; //JPO 180731
     // oss << "<tr>" << endl; //JPO 180731
      oss << "<div class=\"container__cell--full\"><div class=\"container__card\"><h5 class=\"value-name\">Spin Decomposition per atoms:</h5><span class=\"value\"> {" << aentry.spinD << "} &mu;<sub>B</sub> </span></div></div>" << endl; //JPO 180731
      oss << "<div class=\"container__cell\"><div class=\"container__card\"><h5 class=\"value-name\">Spin Polarization (E<sub>F</sub>):</h5><span class=\"value\"> " << aentry.spinF << " </span></div></div>" << endl; //JPO 180731
      oss << "<div class=\"container__cell\"><div class=\"container__card\" id=\"band_gap_td\"><h5 class=\"value-name\">Band Gap:</h5><span class=\"value\">" << aentry.Egap << " eV "
	  << (aentry.Egap_type=="insulator_direct"?"(insulator)":"") << (aentry.Egap_type=="insulator_indirect"?"(insulator)":"") << (aentry.Egap_type=="metal"?"(metal)":"") << "</span></div></div>" << endl; //JPO 180731
      oss << "<div class=\"container__cell\"><div class=\"container__card\" id=\"fit_band_gap_td\"><h5 class=\"value-name\">Fit Band Gap:</h5><span class=\"value\">" << aentry.Egap_fit << " eV</span></div></div>" << endl; //JPO 180731
     // oss << "</tr>" << endl;  //JPO 180731
     // oss << "<tr>" << endl; //JPO 180731
      oss << "<div class=\"container__cell\"><div class=\"container__card\" id=\"FIX_m_atom_td\"><h5 class=\"value-name\">Magnetic Moment:</h5><span class=\"value\">" << aentry.spin_cell << " &mu;<sub>B</sub></span></div></div>" << endl; //JPO 180731
      oss << "<div class=\"container__cell\"><div class=\"container__card\" id=\"FIX_m_atom_td\"><h5 class=\"value-name\">Magnetic Moment/atom:</h5><span class=\"value\">" << aentry.spin_atom << " &mu;<sub>B</sub>/atom</span></div></div>" << endl; //JPO 180731
     // oss << "</tr>" << endl; //JPO 180731
     // oss << "<tr>" << endl; //JPO 180731
      oss << "<div class=\"container__cell\"><div class=\"container__card\" id=\"FIX_e_mass_td\"><h5 class=\"value-name\">Electron Mass(FIX):</h5><span class=\"value\"> XXX (m<sub>0</sub>)</span></div></div>" << endl; //JPO 180731
      oss << "<div class=\"container__cell\"><div class=\"container__card\" id=\"FIX_hole_mass_td\"><h5 class=\"value-name\">Hole Mass(FIX):</h5><span class=\"value\"> XXX (m<sub>0</sub>)</span></div></div>" << endl; //JPO 180731
     // oss << "</tr>" << endl; //JPO 180731
     // oss << "<tr>" << endl; //JPO 180731
      if(aentry.Egap_type=="insulator_direct")
	oss << "<div class=\"container__cell\"><div class=\"container__card\" id=\"FIX_band_gap_type_td\"><h5 class=\"value-name\">Band Gap Type:</h5><span class=\"value\">  Direct</span></div></div>" << endl; //JPO 180731
      if(aentry.Egap_type=="insulator_indirect")
	oss << "<div class=\"container__cell\"><div class=\"container__card\" id=\"FIX_band_gap_type_td\"><h5 class=\"value-name\">Band Gap Type:</h5><span class=\"value\">  Indirect</span></div></div>" << endl; //JPO 180731
     // oss << "</tr>" << endl; //JPO 180731
     // oss << "<tr>" << endl; //JPO 180731
     // oss << "<td class=\"electronic_table_td\"><span class=\"table_description\">Spin Polarization (E<sub>F</sub>):</span> " << aentry.spinF << " </td>" << endl; // JPO 180731
     // oss << "<td class=\"electronic_table_td\"><span class=\"table_description\">Spin Decomposition per atoms:</span> {" << aentry.spinD << "} &mu;<sub>B</sub> </td>" << endl; // JPO 180731
     // oss << "</tr>" << endl; //JPO 180731
     // oss << "</tbody>" << endl; //JPO 180731
     // oss << "</table>" << endl; //JPO 180731
     oss << "</div>" << endl; //JPO 180731
      oss << "<ul>" << endl;
      oss << "<li>" << endl;
    //oss << "<div class=\"container-subtitle\">" << endl;
      //oss << "<h4 class=\"section-subtitle\">Band Structure</h4></div>" << endl;

    // ****************************************************************************
    // INTERACTIVE BANDS PLOT
    // GEENA
	//oss << "<div class=\"container__cell--full\" >" << endl;  // PC 180515 // JPO 180731
    oss << "<div class=\"flex-container\">" << endl; // JPO 180731
    oss << "<span class=\"pic_description_band\">Band Structure:</span>" << endl; // JPO 180731
	oss << "<div class=\"DosOptions\">" 				<< endl;  // PC 180515
  oss << "<div> Zoom/Pan type:" << endl;  // PC 180515
	oss << "<select id=\"zoomOptions\"> "								<< endl;  // PC 180515
	oss << "<option value=\"both\">Both X & Y</option>"				<< endl;
	oss << "<option value=\"xOnly\">X Only</option>"				<< endl;
	oss << "<option value=\"yOnly\">Y Only</option>"				<< endl;
	oss << "</select></div>" << endl; // PC 180515
	//if ((aentry.spinF!=AUROSTD_NAN) && (!aurostd::isequal(abs(aentry.spinF),0.0,_ZERO_TOL_))){  // PC 180515
	if ((aentry.spinF!=AUROSTD_NAN) && (!aurostd::isequal(abs(aentry.spin_atom),0.0,_ZERO_TOL_))){  // PC 180515
  oss << "<div>Majority/Minority Spin Selection:**" << endl;  // PC 180525
  oss << "<select id =\"spinBandsOptions\">" << endl; // PC 180515
	oss << "<option value=\"bothS\">Both Spins</option>"				<< endl;  // PC 180515
	oss << "<option value=\"majority\">Majority Spin</option>"				<< endl;  // PC 180515
	oss << "<option value=\"minority\">Minority Spin</option>"				<< endl;  // PC 180515
	oss << "</select></div>"								<< endl;  // PC 180515
  }
  oss << "<div class=\"reset\">Reset Zoom</div>" << endl; // PC 180515
  oss << "</div>"				<< endl;  // PC 180515
	if ((aentry.spinF!=AUROSTD_NAN) && (!aurostd::isequal(abs(aentry.spin_atom),0.0,_ZERO_TOL_))){  // PC 180515
  oss << "<div class=\"star\">**For smoother tracing of bands </div>" << endl; //PC 180525
  } //PC 180525
  oss << "<div class=\"plots\">" << endl; // PC 180515
  oss << "<div class=\"Bands_plot\">" << endl;  // PC 180515
	oss << "<svg id=\"bands_wrapper\" width=\"900\" height=\"550\"></svg>"						<< endl;  // PC 180515
	if ((aentry.spinF!=AUROSTD_NAN) && (!aurostd::isequal(abs(aentry.spin_atom),0.0,_ZERO_TOL_))){  // PC 180515
	oss << "<div id=\"bandLegend\" class=\"legend\">"				<< endl;
 	oss << "<text class=\"legendText \">Majority Spin"				<< endl;  // PC 180515
 	oss << "<line class=\"legendLine\" style=\"border-color:black;\"></line>"	        << endl;  // PC 180515
	oss << "</text>" 								<< endl;  // PC 180515
	oss << "<text class=\"legendText \">Minority Spin" 				<< endl;  // PC 180515
	oss << "<line class=\"legendLine min\"></line>" 					<< endl;  // PC 180515
	oss << "</text></div>"							<< endl;  // PC 180515
  } // PC 180515
	oss << "</div>"  << endl; // PC 180515
  oss << "<div class=\"Dos_plot\">" << endl;  // PC 180515
	oss << "<svg id=\"dos_wrapper\" width=\"300px\" height=\"550px\" ></svg>" << endl;  // PC 180515
  oss << " <div id=\"dosLegend\" class=\"legend\" >" << endl; // PC 180515
	oss << "<text id=\"dosText\"></text>"						<< endl;  // PC 180515
	oss << "</div>" 								<< endl;
	oss << "</div></div></div>" << endl;  // PC 180515
	//[OBSOLETE PC 180515]oss << "<div class=\"legendLine min\"></div>" 					<< endl;
	//[OBSOLETE PC 180515]oss << "</div></div></svg>"							<< endl;
	//[OBSOLETE PC 180515]oss << "<svg id=\"dos_wrapper\">"						<< endl;
	//[OBSOLETE PC 180515]oss << "<div id=\"dosText\"></div>"						<< endl;
	//[OBSOLETE PC 180515]oss << "<div id=\"dosLegend\" class=\"legend \"></div></svg></div></div>"	<< endl;
   // ********************************************************************************** 

      if(aentry.vfiles_WEB.size()>0) 
	for(uint i=0;i<aentry.vfiles_WEB.size();i++)
	  if(aentry.vfiles_WEB.at(i)==label+".png")
	    oss << "<div class = \"picture_band\" id=\"band_dos_pic\"><a id=\"imgPopup\"><img class=\"pic\" src=\"" << url_WEB << "/" << aentry.vfiles_WEB.at(i) << "\" alt=\"Band Structure of " << label << "\" style='display:block;background:white; margin: 0 auto;' /></a></div>" << endl;
      oss << "<div id=\"small_figure\" style=\'display:flex;justify-content:center;\'>" << endl;
      if(aentry.vfiles_WEB.size()>0) 
	for(uint i=0;i<aentry.vfiles_WEB.size();i++)
	  if(aurostd::substring2bool(aentry.vfiles_WEB.at(i),"_PEDOS_") && aurostd::substring2bool(aentry.vfiles_WEB.at(i),".png")) 
	    oss << "<img class=\"pic_small\" src=\"" << url_WEB << "/" << aentry.vfiles_WEB.at(i) << "\" alt=\"Band Structure of pdos_" << label << "\" style='background:white; margin: 0 auto;height:0%;' >" << endl;
      oss << "</div>" << endl;
      oss << "</li></ul>" << endl;
      oss << "</div>" << endl;
      oss << "<!-- Electronic properties: END -->" << endl;
    }

    // ***************************************************************************
    // ELECTRONIC POPUP
    if(vflags.flag("FLAG::FOUND") && vflags.flag("FLAG::ELECTRONIC") && !directory.empty()) {
      oss << "<!-- Electronic popup: BEGIN -->" << endl;
      oss << "</div> <!close global>" << endl;
      oss << "<div id=\"popupImg\">  <! make popup>" << endl;
      oss << "<a id=\"popupImgClose\">" << endl;
      oss << "<img id=\"closeButton\" src=\"img/closeButton.png\" alt=\"close\" ></a>" << endl;
      oss << "<a id=\"backwardArrow\">" << endl;
      oss << "<img src=\"img/backwardArrow.png\" alt=\"backward arrow\" ></a>" << endl;
      oss << "<ul class=\"pic_popup_list\">" << endl;
      oss << "<li class=\"pic_popup\"><img class=\"pic_large\" src=\"" << url_WEB << "/" << label << ".png\" alt=\"Band Structure of " << label << "\"  ></li>" << endl;
      oss << "<li class=\"pic_popup\"><img class=\"pic_large\" src=\"" << url_WEB << "/" << label << "_DOS.png\" alt=\"DOS of " << label << "\"  ></li>" << endl;
      if(aentry.vfiles_WEB.size()>0) 
	for(uint iline=0;iline<aentry.vfiles_WEB.size();iline++)
	  if(aurostd::substring2bool(aentry.vfiles_WEB.at(iline),"PEDOS"))
	    if(aurostd::substring2bool(aentry.vfiles_WEB.at(iline),"png"))
	      oss << "<li class=\"pic_popup\"><img class=\"pic_large\" src=\"" << url_WEB << "/" << aentry.vfiles_WEB.at(iline) << "\" alt=\"Band Structure of pdos_" << label << "\"  ></li>" << endl;
      oss << "</ul>" << endl;
      oss << "<a id=\"forwardArrow\"><img src=\"img/forwardArrow.png\" alt=\"forward arrow\" ></a>" << endl;
      oss << "</div>" << endl;
      oss << "<!-- Electronic popup: END -->" << endl;
    }
    oss << "<!-- Downloadable Files: BEGIN -->" << endl; //JPO 180809
    oss << "<div class=\"container\">" << endl; //JPO 180809
    oss << "<div class=\"container-title\"><h1 class=\"section-title\"> Downloadable Files </h1></div>" << endl; //JPO 180809
    oss << "<div class=\"container__cell--full\"><div class=\"container__card\">" << endl; //JPO 180809
    oss << "<ul class=\"file-list\">" << endl; //JPO 180809
    oss << "<li class=\"file-name\"><a href=\"" << url_WEB << "/\"" << html_TAB << " download=\"aflowlib.out\" >aflowlib.out</a></li>" << endl; //JPO 180809
    oss << "<li class=\"file-name\"><a href=\"" << url_WEB << "/?format=json\"" << html_TAB << "download=\"aflowlib.json\">aflowlib.json</a></li>" << endl; //JPO 180809
    oss << "<li class=\"file-name\"><a href=\"" << url_WEB << "/CONTCAR.relax.vasp\"" << html_TAB << "download=\"CONTCAR.relax.vasp\">Relaxed Position (aflowlib/VASP)</a></li>" << endl; //JPO 180809
    oss << "<li class=\"file-name\"><a href=\"" << url_WEB << "/CONTCAR.relax.qe\"" << html_TAB << " download=\"CONTCAR.relax.ge\">Relaxed Position (aflowlib/QE)</a></li>" << endl; //JPO 180809
    oss << "<li class=\"file-name\"><a href=\"" << url_WEB << "/CONTCAR.relax.abinit\"" << html_TAB << "download=\"CONTCAR.relax.abinit\">Relaxed Position (aflowlib/AIMS)</a></li>" << endl; //JPO 180809
    oss << "<li class=\"file-name\"><a href=\"" << url_WEB << "/CONTCAR.relax.aims\"" << html_TAB << "download=\"CONTCAR.relax.aims\">Relaxed Position (aflowlib/AIMS)</a>" << endl; //JPO 180809
    oss << "<li class=\"file-name\"><a href=\"" << url_WEB << "/INCAR.relax\"" << html_TAB << "download=\"INCAR.relax\">INCAR for relaxed calculation</a></li>" << endl; //JPO 180809
    oss << "<li class=\"file-name\"><a href=\"" << url_WEB << "/INCAR.static\"" << html_TAB << "download=\"INCAR.static\">INCAR for static calculation</a></li>" << endl; //JPO 180809
    oss << "<li class=\"file-name\"><a href=\"" << url_WEB << "/INCAR.bands\"" << html_TAB << "download=\"INCAR.bands\">INCAR for bands calculation</a></li>" << endl; //JPO 180809
    oss << "<li class=\"file-name\"><a href=\"" << url_WEB << "/KPOINTS.relax\"" << html_TAB << "download=\"KPOINTS.relax\">KPOINTS for relaxed calculation</a></li>" << endl; //JPO 180809
    oss << "<li class=\"file-name\"><a href=\"" << url_WEB << "/KPOINTS.static\"" << html_TAB << "download=\"KPOINTS.static\">KPOINTS for static calculation</a></li>" << endl; //JPO 180809
    oss << "<li class=\"file-name\"><a href=\"" << url_WEB << "/KPOINTS.bands\"" << html_TAB << "download=\"KPOINTS.bands\">KPOINTS for bands calculation</a></li>" << endl; //JPO 180809
    oss << "<li class=\"file-name\"><a href=\"" << url_WEB << "/" << DEFAULT_FILE_EDATA_ORIG_OUT << "\"" << html_TAB << "download=\"edata.orig.out\">Extended crystallographic data for original structure</a></li>" << endl; //JPO 180809
    oss << "<li class=\"file-name\"><a href=\"" << url_WEB << "/" << DEFAULT_FILE_EDATA_RELAX_OUT << "\"" << html_TAB << "download=\"edata.relax.out\">Extended crystallographic data for relaxed structure</a></li>" << endl; //JPO 180809
    oss << "<li class=\"file-name\"><a href=\"" << url_WEB << "/" << DEFAULT_FILE_EDATA_BANDS_OUT << "\"" << html_TAB << "download=\"edata.bands.out\">Extended crystallographic data for band-structure</a></li>" << endl; //JPO 180809
    if(aurostd::substring2bool(aentry.vfiles_WEB,"aflow.agl.out"))  //JPO 180809
    oss << "<li class=\"file-name\"><a href=\"" << url_WEB << "/aflow.agl.out\"" << html_TAB << "download=\"aflow.agl.out\">AGL Output</a></li>" << endl; //JPO 180809
    if(aurostd::substring2bool(aentry.vfiles_WEB,"AGL.out")) //JPO 180809
    oss << "<li class=\"file-name\"><a href=\"" << url_WEB << "/AGL.out\"" << html_TAB << "download=\"AGL.out\">AGL complete output</a></li>" << endl; //JPO 180809
    if(aurostd::substring2bool(aentry.vfiles_WEB,"AGL_energies_temperature.out")) //JPO 180809
    oss << "<li class=\"file-name\"><a href=\"" << url_WEB << "/AGL_energies_temperature.out\"" << html_TAB << "download=\"AGL_energies_temperature.out\">AGL Energy versus Temperature</a></li>" << endl; //JPO 180809
    if(aurostd::substring2bool(aentry.vfiles_WEB,"AGL_thermal_properties_temperature.out")) //JPO 180809
    oss << "<li class=\"file-name\"><a href=\"" << url_WEB << "/AGL_thermal_properties_temperature.out\"" << html_TAB << "download=\"AGL_thermal_properties_temperature.out\">AGL Thermal Properties versus Temperature</a></li>" << endl; //JPO 180809
    if(aurostd::substring2bool(aentry.vfiles_WEB,"aflow.ael.out")) //JPO 180809
    oss << "<li class=\"file-name\"><a href=\"" << url_WEB << "/aflow.ael.out\"" << html_TAB << "download=\"aflow.ael.out\">AEL Output</a></li>" << endl; //JPO 180809
    if(aurostd::substring2bool(aentry.vfiles_WEB,"AEL_Elastic_constants.out")) //JPO 180809
    oss << "<li class=\"file-name\"><a href=\"" << url_WEB << "/AEL_Elastic_constants.out\"" << html_TAB << "download=\"AEL_Elastic_constants.out\">AEL Elastic constants</a></li>" << endl; //JPO 180809
    if(aurostd::substring2bool(aentry.vfiles_WEB,"AEL_Compliance_tensor.out")) //JPO 180809
    oss << "<li class=\"file-name\"><a href=\"" << url_WEB << "/AEL_Compliance_tensor.out\"" << html_TAB << "download=\"AEL_Compliance_tensor.out\">AEL Compliance constants</a></li>" << endl; //JPO 180809
    string abader_out=aentry.prototype+"_abader.out"; //JPO 180809
    if(aurostd::substring2bool(aentry.vfiles_WEB,abader_out)) { //JPO 180809
    oss << "<li class=\"file-name\"><a href=\"" << url_WEB << "/" << abader_out << "\"" << html_TAB << "download=\"abader.out\">Bader Output</a></li>" << endl; } //JPO 180809
    oss << "</ul></div></div>" << endl; //JPO 180809
    oss << "<!-- Downloadable Files: END -->" << endl; //JPO 180809
    
    // ***************************************************************************
    // OTHER
    if(0) if((XHOST.hostname=="nietzsche.mems.duke.edu" || XHOST.hostname=="materials.duke.edu" )&& !directory.empty()) {
	oss << line_rule << endl;
	oss << "<!-- DEBUG: BEGIN -->" << endl;
	oss << "<b>DEBUG: only in " << XHOST.hostname << "</b><br>" << endl;
	oss << aentry.entry << "<br><br>" << endl;
	oss << line_rule << endl;
	oss << "<!-- DEBUG: END -->" << endl;
      }
    // DONE

    oss << "<! HARVEY WORK AFTER HERE> " << endl;

    return aentry.ventry.size();
  }
}



#endif  // _AURO_IMPLEMENTATIONS_

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *                                                                         *
// ***************************************************************************

