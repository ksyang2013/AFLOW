// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *                Aflow CORMAC TOHER - Duke University 2013-2018           *
// *                                                                         *
// ***************************************************************************
// Written by Cormac Toher
// cormac.toher@duke.edu
#ifndef _AFLOW_AEL_GET_STRESS_CPP
#define _AFLOW_AEL_GET_STRESS_CPP
#include "aflow.h"
#include "aflow_ael_elasticity.h"

// ###############################################################################
//                  AFLOW Automatic Elasticity Library (AEL) (2014-2018)
// ###############################################################################
//
// Uses strain-stress calculations to obtain elastic constants of materials
// Based on original Python program written by M. de Jong et al.
// See Scientific Data 2, 150009 (2015) for details of original program
// See Phys. Rev. Materials 1, 015401 (2017) for details of this implementation
// Please cite these works in addition to the general AFLOW papers if you use results generated using AEL
//

// *******************************************************************************
// The following functions are for generating _AFLOWIN_ files
// *******************************************************************************

// ***************************************************************************
// AEL_functions::aelvaspflags
// ***************************************************************************
namespace AEL_functions {
  //
  // Function to assign values for VASP input flags from aflow.in file to vaspRun _xvasp class
  // Adapted from section of AFLOW APL function DirectMethodPC::runVASPCalculations()
  //
  uint aelvaspflags(_xvasp& vaspRun, _vflags& _vaspFlags, _kflags& _kbinFlags, string& runname, _AEL_data& AEL_data, ofstream& FileMESSAGE) {
    ostringstream aus;
    vector<string> vfile;
    string vfilename;
    bool vfileexist = false;
    if(AEL_data.relax_static || AEL_data.static_only) {
      aurostd::string2tokens(string("OUTCAR.static.bz2,OUTCAR.static.gz,OUTCAR.static"),vfile,",");
      for(uint ij=0;ij<vfile.size();ij++) {
	if(aurostd::FileExist(vaspRun.Directory+"/"+vfile.at(ij))) {
	  vfilename = vfile.at(ij);
	  vfileexist = true;
	}    
      }  
    } else {
      aurostd::string2tokens(string("OUTCAR.relax2.bz2,OUTCAR.relax2.gz,OUTCAR.relax2"),vfile,",");
      for(uint ij=0;ij<vfile.size();ij++) {
	if(aurostd::FileExist(vaspRun.Directory+"/"+vfile.at(ij))) {
	  vfilename = vfile.at(ij);
	  vfileexist = true;
	}    
      }  
    }
    // SOME WARNINGS: check existence of LOCK and OUTCAR.relax2 files
    // OBSOLETE if( !aurostd::FileExist( vaspRun.Directory + string("/LOCK") ) &&
    if( !aurostd::FileExist( vaspRun.Directory + "/" + _AFLOWLOCK_ ) &&
	( vfileexist ) ) {
      //OBSOLETE	aurostd::FileExist( vaspRun.Directory + string("/OUTCAR.relax2") ) ) {
      aurostd::StringstreamClean(aus);
      // OBSOLETE aus << _AELSTR_WARNING_ + "found OUTCAR.static but no LOCK in " <<  vaspRun.Directory << endl;
      aus << _AELSTR_WARNING_ + "found " << vfilename << " but no " << _AFLOWLOCK_ << " in " <<  vaspRun.Directory << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      return 1;
    }

    // OBSOLETE if( aurostd::FileExist( vaspRun.Directory + string("/LOCK") ) &&
    if( aurostd::FileExist( vaspRun.Directory + "/" + _AFLOWLOCK_ ) &&
	!(vfileexist) ) {
      //OBSOLETE	!aurostd::FileExist( vaspRun.Directory + string("/OUTCAR.relax2") ) ) {
      aurostd::StringstreamClean(aus);
      // OBSOLETE aus << _AELSTR_WARNING_ + "found LOCK but no OUTCAR.static in " <<  vaspRun.Directory << endl;
      aus << _AELSTR_WARNING_ + "found " << _AFLOWLOCK_ << " but no OUTCAR in " <<  vaspRun.Directory << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      return 1;
    }
	  	   	    
    // Switch off autotune
    _kbinFlags.KBIN_MPI_AUTOTUNE = true;

    // Common KPOINTS settings and OVERRIDES
    vaspRun.AVASP_KSCHEME = _vaspFlags.KBIN_VASP_KPOINTS_KSCHEME.content_string;
    vaspRun.AVASP_value_KPPRA = _vaspFlags.KBIN_VASP_KPOINTS_KPPRA.content_int;
    vaspRun.AVASP_STATIC_KSCHEME = _vaspFlags.KBIN_VASP_KPOINTS_STATIC_KSCHEME.content_string;
    vaspRun.AVASP_value_KPPRA_STATIC = _vaspFlags.KBIN_VASP_KPOINTS_STATIC_KPPRA.content_int;


    // Clear old INCAR and set it as we want...
    // Want to relax ions only keeping cell size and shape fixed
    // Might want to create new relaxation type instead of creating INCAR by hand
    vaspRun.INCAR.str(std::string());
    string system;
    for(uint j=0; j < vaspRun.str.species.size(); j++)
      system = system + vaspRun.str.species_pp.at(j);
    system = system + "@" + runname;
    vaspRun.INCAR << "SYSTEM=" << system << std::endl;
    vaspRun.INCAR << "# Added by [AFLOW_AEL] begin" << std::endl;
    vaspRun.INCAR << "NELMIN=4         # The forces have to be well converged" << std::endl;
    vaspRun.INCAR << "NELM = 120       # Many electronic steps (SC2013)" << std::endl;
    vaspRun.INCAR << "ADDGRID=.TRUE.   # For finer forces" << std::endl;
    // OBSOLETE vaspRun.INCAR << "ISIF=2           # To calculate stress tensor including ion relaxation only" << std::endl;
    // OBSOLETE vaspRun.INCAR << "IBRION=2         # Ion relaxation using conjugate gradient" << std::endl;
    // OBSOLETE vaspRun.INCAR << "NSW=51           # Relax ions for long" << std::endl;
    vaspRun.INCAR << "# Added by [AFLOW_AEL] end" << std::endl;

    // Change format of POSCAR
    if( ( !_kbinFlags.KBIN_MPI && ( _kbinFlags.KBIN_BIN.find("46") != string::npos ) ) ||
	(  _kbinFlags.KBIN_MPI && ( _kbinFlags.KBIN_MPI_BIN.find("46") != string::npos ) ) ) {
      vaspRun.str.is_vasp5_poscar_format = false; 
    }

    return 0;
  }
} // namespace AEL_functions

// ***************************************************************************
// AEL_functions::createAFLOWIN
// ***************************************************************************
namespace AEL_functions {
  //
  // Create aflow.in file: makes new directory and writes aflow.in for strained structure file inside it 
  // Adapted from that in AFLOW APL function PhononCalculator::createAFLOWIN()
  //
  uint createAFLOWIN(_xvasp& vaspRun, _xvasp& xvasp, _kflags& _kbinFlags, _vflags& _vaspFlags, _AEL_data& AEL_data, ofstream& FileMESSAGE) {
    bool AFLOWIN_QE_FLAG=FALSE;
    bool SPACES=FALSE;
    ostringstream aus;

    if( !aurostd::FileExist( vaspRun.Directory) ) {
      aurostd::DirectoryMake( vaspRun.Directory );
    }
    // CHMOD Directory 777: change directory permissions to read+write+execute for all users
    aurostd::DirectoryChmod("777", vaspRun.Directory);

    // Create file
    // OBSOLETE string filename =  vaspRun.Directory + string("/aflow.in");
    string filename =  vaspRun.Directory + "/" + _AFLOWIN_;

    // Check if aflow.in file exists in the directory     
    // If the aflow.in file does exist and the overwrite has not been enabled, exit createAFLOWIN for this structure
    if( aurostd::FileExist( filename) && (!AEL_data.aflowin_overwrite) ) {
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "Not overwriting existing file " << _AFLOWIN_ << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      return 0;
    }
    ofstream outfile(filename.c_str(),ios_base::out);

    // Check aflow.in file is open
    if( !outfile.is_open() ) {
      aurostd::StringstreamClean(aus);
      // OBSOLETE aus << _AELSTR_WARNING_ + "Cannot create [aflow.in] file" << endl;
      aus << _AELSTR_WARNING_ + "Cannot create [" << _AFLOWIN_ << "] file" << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      return 1;
    } 
    // CHMOD a+rw aflow.in: change permissions on aflow.in file
    aurostd::ChmodFile("a+rw",filename);

    // Write to aflow.in file
    if(SPACES) { outfile << std::endl; }
    outfile << AFLOWIN_SEPARATION_LINE << std::endl;
    outfile << "[AFLOW]  _  ___ _" << std::endl;
    outfile << "[AFLOW] / \\|   | \\ |" << std::endl;
    outfile << "[AFLOW] | o |-- | " << std::endl;
    outfile << "[AFLOW] |_n_|__ |___| automatic generated file" << std::endl;
    outfile << "[AFLOW]" << std::endl;
    outfile << AFLOWIN_SEPARATION_LINE << std::endl;
    if(SPACES) { outfile << std::endl; }
    outfile << "[AFLOW_MODE=VASP]" << std::endl;
    outfile << "[AFLOW_MODE_ZIP=" << _kbinFlags.KZIP_BIN << "]" << std::endl;
    if(SPACES) { outfile << std::endl; }
    
    //CO 180130 - START
    //adding aflow.rc stuff
    outfile << "[AFLOW_MODE_BINARY=";
    if(!_kbinFlags.KBIN_BIN.empty()){outfile << _kbinFlags.KBIN_BIN;}
    else{outfile << DEFAULT_VASP_BIN;}
    outfile << "]" << std::endl;
    outfile << AFLOWIN_SEPARATION_LINE << std::endl;
    outfile << AFLOWIN_SEPARATION_LINE << std::endl;
    if(!(_kbinFlags.KBIN_MPI || XHOST.MPI)){outfile << "#";}
    outfile << "[AFLOW_MODE_MPI]" << std::endl;
    //be super cautious and avoid empty tags here
    string NCPUS_VAL="MAX";
    if(XHOST.vflag_control.flag("XPLUG_NUM_THREADS")){NCPUS_VAL=XHOST.vflag_control.getattachedscheme("XPLUG_NUM_THREADS");}
    outfile << "[AFLOW_MODE_MPI_MODE]NCPUS=" << NCPUS_VAL << " " << std::endl;
    outfile << "[AFLOW_MODE_MPI_MODE]COMMAND =\"" << MPI_COMMAND_DEFAULT << "\" " << std::endl;
    if( _kbinFlags.KBIN_MPI_AUTOTUNE ) {outfile << "[AFLOW_MODE_MPI_MODE]AUTOTUNE " << std::endl;}
    outfile << "[AFLOW_MODE_MPI_MODE]BINARY=\"";
    if(!_kbinFlags.KBIN_MPI_BIN.empty()){outfile << _kbinFlags.KBIN_MPI_BIN;}
    else{outfile << DEFAULT_VASP_MPI_BIN;}
    outfile << "\"" << std::endl;
    outfile << AFLOWIN_SEPARATION_LINE << std::endl;
    //CO 180130 - STOP

    //CO 180130 - making obsolete with lines above
    //[OBSOLETE]if( _kbinFlags.KBIN_MPI ) {
    //[OBSOLETE]  outfile << AFLOWIN_SEPARATION_LINE << std::endl;
    //[OBSOLETE]  outfile << "[AFLOW_MODE_BINARY=" << _kbinFlags.KBIN_BIN << "]" << std::endl;
    //[OBSOLETE]  outfile << AFLOWIN_SEPARATION_LINE << std::endl;
    //[OBSOLETE]  if(SPACES) { outfile << std::endl; }
    //[OBSOLETE]  outfile << "[AFLOW_MODE_MPI]" << std::endl;
    //[OBSOLETE]  if( _kbinFlags.KBIN_MPI_AUTOTUNE ) {
    //[OBSOLETE]outfile << "[AFLOW_MODE_MPI_MODE]AUTOTUNE" << std::endl;
    //[OBSOLETE]  } else {
    //[OBSOLETE]outfile << "[AFLOW_MODE_MPI_MODE]NCPUS=MAX" << std::endl;
    //[OBSOLETE]  }
    //[OBSOLETE]  outfile << "[AFLOW_MODE_MPI_MODE]BINARY=" << _kbinFlags.KBIN_MPI_BIN << std::endl;
    //[OBSOLETE]} else {
    //[OBSOLETE]  outfile << AFLOWIN_SEPARATION_LINE << std::endl;
    //[OBSOLETE]  outfile << "[AFLOW_MODE_BINARY=" << _kbinFlags.KBIN_BIN << "]" << std::endl;
    //[OBSOLETE]  outfile << AFLOWIN_SEPARATION_LINE << std::endl;
    //[OBSOLETE]  outfile << "[AFLOW_MODE_MPI_MODE]BINARY=\"mpi" << _kbinFlags.KBIN_BIN << "\"" << std::endl;
    //[OBSOLETE]  outfile << "[AFLOW_MODE_MPI_MODE]NCPUS=MAX" << std::endl;
    //[OBSOLETE]  outfile << "[AFLOW_MODE_MPI_MODE]COMMAND=\"mpirun -np\" " << std::endl;
    //[OBSOLETE]  outfile << "[AFLOW_MODE_MPI_MODE]AUTOTUNE" << std::endl;
    //[OBSOLETE]}
    if(SPACES) { outfile << std::endl; }
    
    // Write INCAR lines to aflow.in file 
    outfile << AFLOWIN_SEPARATION_LINE << std::endl;
    if(AEL_data.relax_static) {
      outfile << "[VASP_RUN]RELAX_STATIC=2" << std::endl;
    } else if(AEL_data.static_only) {
      outfile << "[VASP_RUN]STATIC" << std::endl;      
    } else if(AEL_data.relax_only) {
      outfile << "[VASP_RUN]RELAX=2" << std::endl;
    } else {
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_ERROR_ + "No run type selected" << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      return 1;
    }

    if(_vaspFlags.KBIN_VASP_FORCE_OPTION_LDAU1.isentry) {
      outfile << "[VASP_FORCE_OPTION]LDAU1=ON"  << std::endl;
      outfile << "[VASP_FORCE_OPTION]LDAU_PARAMETERS= " << _vaspFlags.KBIN_VASP_LDAU_PARAMETERS  << std::endl;
    }
    if(_vaspFlags.KBIN_VASP_FORCE_OPTION_LDAU2.isentry) {
      outfile << "[VASP_FORCE_OPTION]LDAU2=ON " << std::endl;
      outfile << "[VASP_FORCE_OPTION]LDAU_PARAMETERS= " <<  _vaspFlags.KBIN_VASP_LDAU_PARAMETERS  << std::endl;
    }
    if(SPACES) { outfile << std::endl; }
    outfile << "[VASP_FORCE_OPTION]RELAX_IONS" << std::endl;
    outfile << "[VASP_FORCE_OPTION]WAVECAR=OFF" << std::endl;
    outfile << "[VASP_FORCE_OPTION]CHGCAR=OFF" << std::endl;
    if(AEL_data.precaccalgonorm) {      
      outfile << "[VASP_FORCE_OPTION]PREC=ACCURATE" << std::endl;
      // OBSOLETE outfile << "[VASP_FORCE_OPTION]PREC=HIGH" << std::endl;
      outfile << "[VASP_FORCE_OPTION]ALGO=NORMAL" << std::endl;
      // OBSOLETE outfile << "[VASP_FORCE_OPTION]ALGO=FAST" << std::endl;
    } else {
      outfile << "[VASP_FORCE_OPTION]PREC=" << _vaspFlags.KBIN_VASP_FORCE_OPTION_PREC.content_string << std::endl;
      outfile << "[VASP_FORCE_OPTION]ALGO=" << _vaspFlags.KBIN_VASP_FORCE_OPTION_ALGO.content_string << std::endl;
    }
    //outfile << "[VASP_FORCE_OPTION]RELAX_MODE=ENERGY" << std::endl;

    // Switch off VASP symmetry - this can help when applied strains break the symmetry of the primitive cell
    if(!AEL_data.vasp_symmetry) {
      outfile << "[VASP_FORCE_OPTION]SYM=OFF" << std::endl;
    }

    // OBSOLETE if( _vaspFlags.KBIN_VASP_FORCE_OPTION_SYM.option ) outfile << "[VASP_FORCE_OPTION]SYM=ON" << std::endl;

    if( _vaspFlags.KBIN_VASP_FORCE_OPTION_AUTO_PSEUDOPOTENTIALS.isentry ) { outfile << "[VASP_FORCE_OPTION]AUTO_PSEUDOPOTENTIALS=" << _vaspFlags.KBIN_VASP_FORCE_OPTION_AUTO_PSEUDOPOTENTIALS.scheme << std::endl; }
    if( _vaspFlags.KBIN_VASP_FORCE_OPTION_ABMIX.isentry ) { outfile << "[VASP_FORCE_OPTION]ABMIX=" << _vaspFlags.KBIN_VASP_FORCE_OPTION_ABMIX.scheme << std::endl; }
    if( _vaspFlags.KBIN_VASP_FORCE_OPTION_TYPE.isentry ) { outfile << "[VASP_FORCE_OPTION]TYPE=" << _vaspFlags.KBIN_VASP_FORCE_OPTION_TYPE.scheme << std::endl; }
    if( _vaspFlags.KBIN_VASP_FORCE_OPTION_AUTO_MAGMOM.isentry ) { outfile << "[VASP_FORCE_OPTION]AUTO_MAGMOM=" << (_vaspFlags.KBIN_VASP_FORCE_OPTION_AUTO_MAGMOM.option?"ON":"OFF") << std::endl; }
    if( _vaspFlags.KBIN_VASP_FORCE_OPTION_BADER.isentry &&_vaspFlags.KBIN_VASP_FORCE_OPTION_BADER.option) {
      outfile << "[VASP_FORCE_OPTION]BADER=ON" << std::endl; 
    } else { 
      outfile << "[VASP_FORCE_OPTION]BADER=OFF" << std::endl;
    }
    if( _vaspFlags.KBIN_VASP_FORCE_OPTION_SPIN.isentry ) {
      if(_vaspFlags.KBIN_VASP_FORCE_OPTION_SPIN.option) { 
	outfile << "[VASP_FORCE_OPTION]SPIN=ON" << std::endl;
      } else {
	outfile << "[VASP_FORCE_OPTION]SPIN=OFF" << std::endl;
      }
    }
    else { outfile << "[VASP_FORCE_OPTION]IGNORE_AFIX=NPARC" << std::endl; }

    if(SPACES) { outfile << std::endl; }
    outfile << "[VASP_INCAR_MODE_EXPLICIT]START" << std::endl;
    outfile << vaspRun.INCAR.str();
    outfile << "[VASP_INCAR_MODE_EXPLICIT]STOP" << std::endl;
    if(SPACES) { outfile << std::endl; }
    
    // Write KPOINTS related lines to aflow.in file
    outfile << AFLOWIN_SEPARATION_LINE << std::endl;
    if(SPACES) { outfile << std::endl; }
    outfile << "[VASP_KPOINTS_MODE_IMPLICIT] " << std::endl;
    outfile << "[VASP_KPOINTS_FILE]KSCHEME=" << vaspRun.AVASP_KSCHEME << " " << std::endl;
    outfile << "[VASP_KPOINTS_FILE]KPPRA=" << vaspRun.AVASP_value_KPPRA << std::endl;
    outfile << "[VASP_KPOINTS_FILE]STATIC_KSCHEME=" << vaspRun.AVASP_STATIC_KSCHEME << " " << std::endl;
    outfile << "[VASP_KPOINTS_FILE]STATIC_KPPRA=" << vaspRun.AVASP_value_KPPRA_STATIC << std::endl;
    if(SPACES) { outfile << std::endl; }
    
    // Write POTCAR related lines to aflow.in file
    outfile << AFLOWIN_SEPARATION_LINE << std::endl;
    if(SPACES) { outfile << std::endl; }
    outfile << "[VASP_POTCAR_MODE_IMPLICIT] " << std::endl;
    string pp;
    
    for(uint j=0; j < xvasp.str.species_pp.size(); j++) {
      if(!_vaspFlags.KBIN_VASP_FORCE_OPTION_AUTO_PSEUDOPOTENTIALS.isentry) { outfile << "[VASP_POTCAR_FILE]" << xvasp.str.species_pp.at(j) << std::endl; }
      if(_vaspFlags.KBIN_VASP_FORCE_OPTION_AUTO_PSEUDOPOTENTIALS.isentry)  { outfile << "[VASP_POTCAR_FILE]" << xvasp.str.species.at(j) << std::endl; }
    }
    // OBSOLETE   aurostd::StringstreamClean(aus);
    // OBSOLETE   for(uint j=0; j < xvasp.str.species_pp.size(); j++) {
    // OBSOLETE     aus << _AGLSTR_MESSAGE_ + "Species_pp " << j << " = " << xvasp.str.species_pp.at(j) << endl;
    // OBSOLETE     aus << _AGLSTR_MESSAGE_ + "Species " << j << " = " << xvasp.str.species.at(j) << endl;
    // OBSOLETE   }
    // OBSOLETE   aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    // OBSOLETE   for(uint j=0; j < xvasp.str.species_pp.size(); j++) {
    // OBSOLETE  outfile << "[VASP_POTCAR_FILE]" << xvasp.str.species_pp.at(j) << std::endl;
    // OBSOLETE}

    if(SPACES) { outfile << std::endl; }

    // Write POSCAR lines to aflow.in file
    outfile << AFLOWIN_SEPARATION_LINE << std::endl;
    if(SPACES) { outfile << std::endl; }
    outfile << "[VASP_POSCAR_MODE_EXPLICIT]START " << std::endl;
    vaspRun.str.is_vasp4_poscar_format=TRUE;
    vaspRun.str.is_vasp5_poscar_format=FALSE;
    outfile << vaspRun.str;
    outfile << "[VASP_POSCAR_MODE_EXPLICIT]STOP " << std::endl;
    if(SPACES) { outfile << std::endl; }

    //
    // OBSOLETE outfile << AFLOWIN_SEPARATION_LINE << std::endl;
    // OBSOLETE if(SPACES) outfile << std::endl;
    // OBSOLETE outfile << "[VASP_RUN]STATIC" << std::endl;
    // OBSOLETE if(SPACES) outfile << std::endl;
    outfile << AFLOWIN_SEPARATION_LINE << std::endl;

    if(AFLOWIN_QE_FLAG) {
      outfile << AFLOWIN_SEPARATION_LINE << std::endl; 
      outfile << "[QE_GEOM_MODE_EXPLICIT]START " << std::endl;
      // OBSOLETE xstructure qestr(vaspRun.str);qestr.vasp2qe();
      xstructure qestr(vaspRun.str);qestr.xstructure2qe();
      outfile << qestr;
      outfile << "[QE_GEOM_MODE_EXPLICIT]STOP " << std::endl;
    }

    if(SPACES) { outfile << std::endl; }
    outfile.close();
    outfile.clear();

    return 0;
  }
} // namespace AEL_functions

// ************************************************************************************************
// This set of functions extract stress tensor data from VASP runs
// ************************************************************************************************

// ***************************************************************************
// AEL_functions::extractstress
// ***************************************************************************
namespace AEL_functions {
  //
  // extractstress: Extract stress tensors from the completed VASP calculations
  // Adapted from section of AFLOW APL function DirectMethodPC::runVASPCalculations()
  //
  uint extractstress(vector<_xvasp>& vaspRuns, _AEL_data& AEL_data, ofstream& FileMESSAGE) {
    bool LVERBOSE=(FALSE || XHOST.DEBUG);
    ostringstream aus;
    vector<string> vfile, dfile;
    xOUTCAR outcar;
    xVASPRUNXML vasprunxml;
    string vfilename, dfilename, ffilename;
    bool skipdir = false;
    aurostd::xmatrix<double> stress_tensor(3, 3);
    double pressure_val = 0.0, energy_cell_val = 0.0;
    AEL_data.normal_stress.clear();
    AEL_data.shear_stress.clear();
    AEL_data.normal_stress.resize(AEL_data.normal_strain.size());
    AEL_data.shear_stress.resize(AEL_data.shear_strain.size());
    AEL_data.normal_deformations_complete.resize(AEL_data.normal_strain.size());
    AEL_data.shear_deformations_complete.resize(AEL_data.shear_strain.size());
    AEL_data.energycalculated.clear();
    AEL_data.pressurecalculated.clear();
    AEL_data.stresscalculated.clear();
    AEL_data.structurecalculated.clear();    
    if(AEL_data.relax_static || AEL_data.static_only) {
      if(AEL_data.vasprunxmlstress) {
	aurostd::string2tokens(string("vasprun.xml.static.bz2,vasprun.xml.static.gz,vasprun.xml.static"),vfile,",");
      } else {
	aurostd::string2tokens(string("OUTCAR.static.bz2,OUTCAR.static.gz,OUTCAR.static"),vfile,",");
      }
    } else if(AEL_data.relax_only) {
      if(AEL_data.vasprunxmlstress) {
	aurostd::string2tokens(string("vasprun.xml.relax2.bz2,vasprun.xml.relax2.gz,vasprun.xml.relax2"),vfile,",");
      } else {
	aurostd::string2tokens(string("OUTCAR.relax2.bz2,OUTCAR.relax2.gz,OUTCAR.relax2"),vfile,",");
      }
    } else {
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_ERROR_ + "No run type selected" << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      return 1;
    }

    // Loop over normal strains and read stress tensors
    uint idVaspRun = 0;
    for (uint i = 1; i <= AEL_data.normal_strain.size(); i++) {
      for (uint j = 0; j < AEL_data.normal_deformations.size(); j++) {   
	skipdir = false;
	if(idVaspRun > vaspRuns.size()) {
	  aurostd::StringstreamClean(aus);
	  aus <<  _AELSTR_WARNING_ + "idVaspRun = " << idVaspRun << " > vaspRuns.size() = " << vaspRuns.size() << endl;
	  aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
	  return 1;
	}
	double strainfactor = 1.0 + AEL_data.normal_deformations.at(j);
	aurostd::StringstreamClean(aus);
	// OBSOLETE aus << _AELSTR_MESSAGE_ + "System number = " << idVaspRun << ", Directory = " << vaspRuns.at(idVaspRun).Directory.at(idVaspRun) << endl;
	aus << _AELSTR_MESSAGE_ + "System number = " << idVaspRun << ", Directory = " << vaspRuns.at(idVaspRun).Directory << endl;
	aus << _AELSTR_MESSAGE_ + "Normal stress: i = " << i << ", j = " << j << ", strain factor = " << strainfactor << endl;
	aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
	aurostd::string2tokens(vaspRuns.at(idVaspRun).Directory, dfile, "/");
	dfilename = dfile.at(dfile.size()-1);
	aurostd::StringstreamClean(aus);
	// OBSOLETE aus << _AELSTR_MESSAGE_ + "System number = " << idVaspRun << ", Directory = " << vaspRuns.at(idVaspRun).Directory.at(idVaspRun) << endl;
	aus << _AELSTR_MESSAGE_ + "Directory name = " << dfilename << endl;
	aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
	// Check if structure is on list of failed runs to be skipped
	// If so, then skip reading and continue to next structure
	for (uint ij = 0; ij < AEL_data.failed_arun_list.size(); ij++) {
	  ffilename = AEL_data.failed_arun_list.at(ij);
	  if(LVERBOSE) {
	    aurostd::StringstreamClean(aus);
	    aus << _AELSTR_MESSAGE_ + "dfilename = " << dfilename << endl;
	    aus << _AELSTR_MESSAGE_ + "ffilename = " << ffilename << endl;
	    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
	  }
	  // OBSOLETE if(dfilename == AEL_data.failed_arun_list.at(ij)) continue;
	  // OBSOLETE if(strncmp(dfilename.c_str(), ffilename.c_str(), 12) == 0) continue;
	  if(aurostd::substring2bool(dfilename,ffilename,TRUE)) {
	    aurostd::StringstreamClean(aus);
	    aus << _AELSTR_MESSAGE_ + "Found directory in to-skip list: " << dfilename << endl;
	    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
	    skipdir = true;
	  }
	}
	if(skipdir) {
	  aurostd::StringstreamClean(aus);
	  aus << _AELSTR_MESSAGE_ + "Skipping directory: " << dfilename << endl;
	  aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
	  idVaspRun++;
	  skipdir = false;
	  continue;
	}
	// if(skipdir) continue;
	// continue;

	// If tarred and compressed directory exists...
	string tarfilename = vaspRuns.at(idVaspRun).Directory + ".tar.bz2";
	if( aurostd::FileExist(tarfilename) ) {
	  // Extract all...
	  aurostd::execute( string("tar -xf ") + tarfilename );
	}
      
	// If the LOCK file is missing, then it is probably a corrupted run
	// Do not accept it and wait for the new run
	// OBSOLETE if( !aurostd::FileExist( vaspRuns.at(idVaspRun).Directory + string("/LOCK") ) ) {
	if( !aurostd::FileExist( vaspRuns.at(idVaspRun).Directory + "/" + _AFLOWLOCK_ ) ) {
	  aurostd::StringstreamClean(aus);
	  aus <<  _AELSTR_WARNING_ + "The " << _AFLOWLOCK_ << " file in " << vaspRuns.at(idVaspRun).Directory << " directory is missing." << endl;
	  aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
	  if(AEL_data.autoskipfailedaruns) {
	    aurostd::StringstreamClean(aus);
	    aus << _AELSTR_MESSAGE_ + "Skipping directory: " << dfilename << endl;
	    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
	    idVaspRun++;
	    continue;
	  } else {
	    throw AELStageBreak();
	  }
	}
	if(AEL_data.vasprunxmlstress) {
	  for(uint ij=0;ij<vfile.size()&&(vasprunxml.content=="");ij++) {
	    if(aurostd::FileExist(vaspRuns.at(idVaspRun).Directory+"/"+vfile.at(ij))) {
	      if(LVERBOSE) {
		aurostd::StringstreamClean(aus);
		aus << _AELSTR_MESSAGE_ + "vfile = " << vfile.at(ij) << endl;
		aus << _AELSTR_MESSAGE_ + "file = " << vaspRuns.at(idVaspRun).Directory+"/"+vfile.at(ij) << endl;	    
		aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
	      }
	      vasprunxml.GetPropertiesFile(vaspRuns.at(idVaspRun).Directory+"/"+vfile.at(ij));
	      vfilename = vfile.at(ij);
	    }
	  }
	  if(vasprunxml.content=="") {
	    aurostd::StringstreamClean(aus);
	    aus <<  _AELSTR_WARNING_ + "The " << vfilename << " file in " << vaspRuns.at(idVaspRun).Directory << " directory is missing." << endl;
	    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
	    if(AEL_data.autoskipfailedaruns) {
	      aurostd::StringstreamClean(aus);
	      aus << _AELSTR_MESSAGE_ + "Skipping directory: " << dfilename << endl;
	      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
	      idVaspRun++;
	      continue;
	    } else {
	      throw AELStageBreak();
	    }
	  }
	  stress_tensor = -vasprunxml.stress;
	  vasprunxml.clear();
	} else {
	  for(uint ij=0;ij<vfile.size()&&(outcar.content=="");ij++) {
	    if(aurostd::FileExist(vaspRuns.at(idVaspRun).Directory+"/"+vfile.at(ij))) {
	      if(LVERBOSE) {
		aurostd::StringstreamClean(aus);
		aus << _AELSTR_MESSAGE_ + "vfile = " << vfile.at(ij) << endl;
		aus << _AELSTR_MESSAGE_ + "file = " << vaspRuns.at(idVaspRun).Directory+"/"+vfile.at(ij) << endl;	    
		aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
	      }
	      outcar.GetPropertiesFile(vaspRuns.at(idVaspRun).Directory+"/"+vfile.at(ij));
	      vfilename = vfile.at(ij);
	    }
	  }
	  if(outcar.content=="") {
	    aurostd::StringstreamClean(aus);
	    aus <<  _AELSTR_WARNING_ + "The " << vfilename << " file in " << vaspRuns.at(idVaspRun).Directory << " directory is missing." << endl;
	    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
	    if(AEL_data.autoskipfailedaruns) {
	      aurostd::StringstreamClean(aus);
	      aus << _AELSTR_MESSAGE_ + "Skipping directory: " << dfilename << endl;
	      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
	      idVaspRun++;
	      continue;
	    } else {
	      throw AELStageBreak();
	    }
	  }
	  stress_tensor = -outcar.stress;
	  energy_cell_val = outcar.energy_cell;
	  pressure_val = outcar.pressure_residual;
	  outcar.clear();
	}
	AEL_data.normal_stress.at(i-1).push_back(stress_tensor);
	AEL_data.normal_deformations_complete.at(i-1).push_back(AEL_data.normal_deformations.at(j));
	AEL_data.energycalculated.push_back(energy_cell_val);
	AEL_data.pressurecalculated.push_back(pressure_val);
	AEL_data.stresscalculated.push_back(stress_tensor);
	AEL_data.structurecalculated.push_back(idVaspRun);
	idVaspRun++;
	// Print out stress tensor
	aurostd::StringstreamClean(aus);
	aus << _AELSTR_MESSAGE_ + "System number = " << idVaspRun << ", Normal stress tensor = " << AEL_data.normal_stress.at(i-1).at(AEL_data.normal_stress.at(i-1).size()-1) << endl;
	aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      }
    }

    // Loop over shear strains and read stress tensors
    for (uint i = 1; i <= AEL_data.shear_strain.size(); i++) {
      for (uint j = 0; j < AEL_data.shear_deformations.size(); j++) {   
	skipdir = false;
	if(idVaspRun > vaspRuns.size()) {
	  aurostd::StringstreamClean(aus);
	  aus <<  _AELSTR_WARNING_ + "idVaspRun = " << idVaspRun << " > vaspRuns.size() = " << vaspRuns.size() << endl;
	  aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
	  return 1;
	}
	double strainfactor = 1.0 + AEL_data.shear_deformations.at(j);
	aurostd::StringstreamClean(aus);
	// OBSOLETE aus << _AELSTR_MESSAGE_ + "System number = " << idVaspRun << ", Directory = " << vaspRuns.at(idVaspRun).Directory.at(idVaspRun) << endl;
	aus << _AELSTR_MESSAGE_ + "System number = " << idVaspRun << ", Directory = " << vaspRuns.at(idVaspRun).Directory << endl;
	aus << _AELSTR_MESSAGE_ + "Shear stress: i = " << i << ", j = " << j << ", strain factor = " << strainfactor << endl;
	aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
	aurostd::string2tokens(vaspRuns.at(idVaspRun).Directory, dfile, "/");
	dfilename = dfile.at(dfile.size()-1);
	aurostd::StringstreamClean(aus);
	// OBSOLETE aus << _AELSTR_MESSAGE_ + "System number = " << idVaspRun << ", Directory = " << vaspRuns.at(idVaspRun).Directory.at(idVaspRun) << endl;
	aus << _AELSTR_MESSAGE_ + "Directory name = " << dfilename << endl;
	aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
	// Check if structure is on list of failed runs to be skipped
	// If so, then skip reading and continue to next structure
	for (uint ij = 0; ij < AEL_data.failed_arun_list.size(); ij++) {
	  ffilename = AEL_data.failed_arun_list.at(ij);
	  if(LVERBOSE) {
	    aurostd::StringstreamClean(aus);
	    aus << _AELSTR_MESSAGE_ + "dfilename = " << dfilename << endl;
	    aus << _AELSTR_MESSAGE_ + "ffilename = " << ffilename << endl;
	    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
	  }
	  // OBSOLETE if(dfilename == AEL_data.failed_arun_list.at(ij)) continue;
	  // OBSOLETE if(strncmp(dfilename.c_str(), ffilename.c_str(), 12) == 0) continue;
	  if(aurostd::substring2bool(dfilename,ffilename,TRUE)) {
	    aurostd::StringstreamClean(aus);
	    aus << _AELSTR_MESSAGE_ + "Found directory in to-skip list: " << dfilename << endl;
	    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
	    skipdir = true;
	  }
	}
	if(skipdir) {
	  aurostd::StringstreamClean(aus);
	  aus << _AELSTR_MESSAGE_ + "Skipping directory: " << dfilename << endl;
	  aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
	  skipdir = false;
	  idVaspRun++;
	  continue;
	}

	// If tarred and compressed directory exists...
	string tarfilename = vaspRuns.at(idVaspRun).Directory + ".tar.bz2";
	if( aurostd::FileExist(tarfilename) ) {
	  // Extract all...
	  aurostd::execute( string("tar -xf ") + tarfilename );
	}
      
	// If the LOCK file is missing, then it is probably a corrupted run
	// Do not accept it and wait for the new run
	// OBSOLETE if( !aurostd::FileExist( vaspRuns.at(idVaspRun).Directory + string("/LOCK") ) ) {
	if( !aurostd::FileExist( vaspRuns.at(idVaspRun).Directory + "/" + _AFLOWLOCK_ ) ) {
	  aurostd::StringstreamClean(aus);
	  aus <<  _AELSTR_WARNING_ + "The " << _AFLOWLOCK_ << " file in " << vaspRuns.at(idVaspRun).Directory << " directory is missing." << endl;
	  aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
	  if(AEL_data.autoskipfailedaruns) {
	    aurostd::StringstreamClean(aus);
	    aus << _AELSTR_MESSAGE_ + "Skipping directory: " << dfilename << endl;
	    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
	    idVaspRun++;
	    continue;
	  } else {
	    throw AELStageBreak();
	  }
	}
	if(AEL_data.vasprunxmlstress) {
	  for(uint ij=0;ij<vfile.size()&&(vasprunxml.content=="");ij++) {
	    if(aurostd::FileExist(vaspRuns.at(idVaspRun).Directory+"/"+vfile.at(ij))) {
	      if(LVERBOSE) {
		aurostd::StringstreamClean(aus);
		aus << _AELSTR_MESSAGE_ + "vfile = " << vfile.at(ij) << endl;
		aus << _AELSTR_MESSAGE_ + "file = " << vaspRuns.at(idVaspRun).Directory+"/"+vfile.at(ij) << endl;	    
		aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
	      }
	      vasprunxml.GetPropertiesFile(vaspRuns.at(idVaspRun).Directory+"/"+vfile.at(ij));
	      vfilename = vfile.at(ij);
	    }
	  }
	  if(vasprunxml.content=="") {
	    aurostd::StringstreamClean(aus);
	    aus <<  _AELSTR_WARNING_ + "The " << vfilename << " file in " << vaspRuns.at(idVaspRun).Directory << " directory is missing." << endl;
	    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
	    if(AEL_data.autoskipfailedaruns) {
	      aurostd::StringstreamClean(aus);
	      aus << _AELSTR_MESSAGE_ + "Skipping directory: " << dfilename << endl;
	      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
	      idVaspRun++;
	      continue;
	    } else {
	      throw AELStageBreak();
	    }
	  }
	  stress_tensor = -vasprunxml.stress;
	  vasprunxml.clear();
	} else {
	  for(uint ij=0;ij<vfile.size()&&(outcar.content=="");ij++) {
	    if(aurostd::FileExist(vaspRuns.at(idVaspRun).Directory+"/"+vfile.at(ij))) {
	      if(LVERBOSE) {
		aurostd::StringstreamClean(aus);
		aus << _AELSTR_MESSAGE_ + "vfile = " << vfile.at(ij) << endl;
		aus << _AELSTR_MESSAGE_ + "file = " << vaspRuns.at(idVaspRun).Directory+"/"+vfile.at(ij) << endl;	    
		aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
	      }
	      outcar.GetPropertiesFile(vaspRuns.at(idVaspRun).Directory+"/"+vfile.at(ij));
	      vfilename = vfile.at(ij);
	    }
	  }
	  if(outcar.content=="") {
	    aurostd::StringstreamClean(aus);
	    aus <<  _AELSTR_WARNING_ + "The " << vfilename << " file in " << vaspRuns.at(idVaspRun).Directory << " directory is missing." << endl;
	    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
	    if(AEL_data.autoskipfailedaruns) {
	      aurostd::StringstreamClean(aus);
	      aus << _AELSTR_MESSAGE_ + "Skipping directory: " << dfilename << endl;
	      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
	      idVaspRun++;
	      continue;
	    } else {
	      throw AELStageBreak();
	    }
	  }
	  stress_tensor = -outcar.stress;
	  energy_cell_val = outcar.energy_cell;
	  pressure_val = outcar.pressure_residual;
	  outcar.clear();
	}
	AEL_data.shear_stress.at(i-1).push_back(stress_tensor);
	AEL_data.shear_deformations_complete.at(i-1).push_back(AEL_data.shear_deformations.at(j));
	AEL_data.energycalculated.push_back(energy_cell_val);
	AEL_data.pressurecalculated.push_back(pressure_val);
	AEL_data.stresscalculated.push_back(stress_tensor);
	AEL_data.structurecalculated.push_back(idVaspRun);
	// Print out stress tensor
	aurostd::StringstreamClean(aus);
	aus << _AELSTR_MESSAGE_ + "System number = " << idVaspRun << ", Shear stress tensor = " << AEL_data.shear_stress.at(i-1).at(AEL_data.shear_stress.at(i-1).size()-1)  << endl;
	aus << _AELSTR_MESSAGE_ + "System number = " << idVaspRun << ", energy = " << outcar.energy_cell << endl;
	aus << _AELSTR_MESSAGE_ + "System number = " << idVaspRun << ", pressure = " << outcar.pressure_residual << endl;		
	aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
	idVaspRun++;
      }
    }
    if(AEL_data.calcstrainorigin) {
      if(idVaspRun > vaspRuns.size()) {
	aurostd::StringstreamClean(aus);
	aus <<  _AELSTR_WARNING_ + "idVaspRun = " << idVaspRun << " > vaspRuns.size() = " << vaspRuns.size() << endl;
	aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
	return 1;
      }
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "System number = " << idVaspRun << ", Directory = " << vaspRuns.at(idVaspRun).Directory << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);

      // If tarred and compressed directory exists...
      string tarfilename = vaspRuns.at(idVaspRun).Directory + ".tar.bz2";
      if( aurostd::FileExist(tarfilename) ) {
	// Extract all...
	aurostd::execute( string("tar -xf ") + tarfilename );
      }
      
      // If the LOCK file is missing, then it is probably a corrupted run
      // Do not accept it and wait for the new run
      if( !aurostd::FileExist( vaspRuns.at(idVaspRun).Directory + "/" + _AFLOWLOCK_ ) ) {
	aurostd::StringstreamClean(aus);
	aus <<  _AELSTR_WARNING_ + "The " << _AFLOWLOCK_ << " file in " << vaspRuns.at(idVaspRun).Directory << " directory is missing." << endl;
	aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
	throw AELStageBreak();
      }
      if(AEL_data.vasprunxmlstress) {
	for(uint ij=0;ij<vfile.size()&&(vasprunxml.content=="");ij++) {
	  if(aurostd::FileExist(vaspRuns.at(idVaspRun).Directory+"/"+vfile.at(ij))) {
	    if(LVERBOSE) {
	      aurostd::StringstreamClean(aus);
	      aus << _AELSTR_MESSAGE_ + "vfile = " << vfile.at(ij) << endl;
	      aus << _AELSTR_MESSAGE_ + "file = " << vaspRuns.at(idVaspRun).Directory+"/"+vfile.at(ij) << endl;	    
	      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
	    }
	    vasprunxml.GetPropertiesFile(vaspRuns.at(idVaspRun).Directory+"/"+vfile.at(ij));
	    vfilename = vfile.at(ij);
	  }
	}
	if(vasprunxml.content=="") {
	  aurostd::StringstreamClean(aus);
	  aus <<  _AELSTR_WARNING_ + "The " << vfilename << " file in " << vaspRuns.at(idVaspRun).Directory << " directory is missing." << endl;
	  aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
	  throw AELStageBreak();
	}
	stress_tensor = -vasprunxml.stress;
	vasprunxml.clear();
      } else {
	for(uint ij=0;ij<vfile.size()&&(outcar.content=="");ij++) {
	  if(aurostd::FileExist(vaspRuns.at(idVaspRun).Directory+"/"+vfile.at(ij))) {
	    if(LVERBOSE) {
	      aurostd::StringstreamClean(aus);
	      aus << _AELSTR_MESSAGE_ + "vfile = " << vfile.at(ij) << endl;
	      aus << _AELSTR_MESSAGE_ + "file = " << vaspRuns.at(idVaspRun).Directory+"/"+vfile.at(ij) << endl;	    
	      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
	    }
	    outcar.GetPropertiesFile(vaspRuns.at(idVaspRun).Directory+"/"+vfile.at(ij));
	    vfilename = vfile.at(ij);
	  }
	}
	if(outcar.content=="") {
	  aurostd::StringstreamClean(aus);
	  aus <<  _AELSTR_WARNING_ + "The " << vfilename << " file in " << vaspRuns.at(idVaspRun).Directory << " directory is missing." << endl;
	  aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
	  throw AELStageBreak();
	}
	stress_tensor = -outcar.stress;
	energy_cell_val = outcar.energy_cell;
	pressure_val = outcar.pressure_residual;
	outcar.clear();
      }
      AEL_data.origin_stress.push_back(stress_tensor);
      AEL_data.energycalculated.push_back(energy_cell_val);
      AEL_data.pressurecalculated.push_back(pressure_val);
      AEL_data.stresscalculated.push_back(stress_tensor);
      AEL_data.structurecalculated.push_back(idVaspRun); 
      // Print out stress tensor
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "System number = " << idVaspRun << ", Origin stress tensor = " << AEL_data.origin_stress.at(0)  << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      idVaspRun++;
    } else if(AEL_data.fitrelaxedstruct) {
      string relaxedstructdirname = AEL_data.dirpathname + "/";
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "Initial relaxed structure directory = " << relaxedstructdirname << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);

      // If tarred and compressed directory exists...
      string tarfilename = relaxedstructdirname + ".tar.bz2";
      if( aurostd::FileExist(tarfilename) ) {
	// Extract all...
	aurostd::execute( string("tar -xf ") + tarfilename );
      }
      
      if(AEL_data.vasprunxmlstress) {
	for(uint ij=0;ij<vfile.size()&&(vasprunxml.content=="");ij++) {
	  if(aurostd::FileExist(relaxedstructdirname+"/"+vfile.at(ij))) {
	    if(LVERBOSE) {
	      aurostd::StringstreamClean(aus);
	      aus << _AELSTR_MESSAGE_ + "vfile = " << vfile.at(ij) << endl;
	      aus << _AELSTR_MESSAGE_ + "file = " << relaxedstructdirname+"/"+vfile.at(ij) << endl;	    
	      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
	    }
	    vasprunxml.GetPropertiesFile(relaxedstructdirname+"/"+vfile.at(ij));
	    vfilename = vfile.at(ij);
	  }
	}
	if(vasprunxml.content=="") {
	  aurostd::StringstreamClean(aus);
	  aus <<  _AELSTR_WARNING_ + "The " << vfilename << " file in " << relaxedstructdirname << " directory is missing." << endl;
	  aus <<  _AELSTR_ERROR_ + "The flag [AFLOW_AEL]FITRELAXEDSTRUCT=ON is set in the input file." << endl;
	  aus <<  _AELSTR_ERROR_ + "This requires the results of a relaxed calculation to be present in the directory " << relaxedstructdirname << endl;
	  aus <<  _AELSTR_ERROR_ + "Either set the flag to OFF or run AEL in the appropriate directory." << endl;
	  aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
	  throw AELStageBreak();
	}
	stress_tensor = -vasprunxml.stress;
	vasprunxml.clear();
      } else {
	for(uint ij=0;ij<vfile.size()&&(outcar.content=="");ij++) {
	  if(aurostd::FileExist(relaxedstructdirname+"/"+vfile.at(ij))) {
	    if(LVERBOSE) {
	      aurostd::StringstreamClean(aus);
	      aus << _AELSTR_MESSAGE_ + "vfile = " << vfile.at(ij) << endl;
	      aus << _AELSTR_MESSAGE_ + "file = " << relaxedstructdirname+"/"+vfile.at(ij) << endl;	    
	      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
	    }
	    outcar.GetPropertiesFile(relaxedstructdirname+"/"+vfile.at(ij));
	    vfilename = vfile.at(ij);
	  }
	}
	if(outcar.content=="") {
	  aurostd::StringstreamClean(aus);
	  aus <<  _AELSTR_WARNING_ + "The " << vfilename << " file in " << relaxedstructdirname << " directory is missing." << endl;
	  aus <<  _AELSTR_ERROR_ + "The flag [AFLOW_AEL]FITRELAXEDSTRUCT=ON is set in the input file." << endl;
	  aus <<  _AELSTR_ERROR_ + "This requires the results of a relaxed calculation to be present in the directory " << relaxedstructdirname << endl;
	  aus <<  _AELSTR_ERROR_ + "Either set the flag to OFF or run AEL in the appropriate directory." << endl;
	  aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
	  throw AELStageBreak();
	}
	stress_tensor = -outcar.stress;
	energy_cell_val = outcar.energy_cell;
	pressure_val = outcar.pressure_residual;
	outcar.clear();
      }
      AEL_data.origin_stress.push_back(stress_tensor);
      AEL_data.energycalculated.push_back(energy_cell_val);
      AEL_data.pressurecalculated.push_back(pressure_val);
      AEL_data.stresscalculated.push_back(stress_tensor);
      AEL_data.structurecalculated.push_back(idVaspRun);
      // Print out stress tensor
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "System number = " << idVaspRun << ", Origin stress tensor = " << AEL_data.origin_stress.at(0)  << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      idVaspRun++;
    }        
    return 0;
  }
} // namespace AEL_functions

// **************************************************************************
//  End of AFLOW AEL set-up and extract stress-strain data
// **************************************************************************

#endif  // _AFLOW_AEL_GET_STRESS_CPP
