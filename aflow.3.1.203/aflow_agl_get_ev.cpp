// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *                Aflow CORMAC TOHER - Duke University 2013-2018           *
// *                                                                         *
// ***************************************************************************
// Written by Cormac Toher
// cormac.toher@duke.edu
#ifndef _AFLOW_AGL_GET_EV_CPP
#define _AFLOW_AGL_GET_EV_CPP
#include "aflow.h"
#include "aflow_agl_debye.h"

 
// ###############################################################################
//                  AFLOW Automatic GIBBS Library (AGL) (2013-2018)
// ###############################################################################
//
// Uses quasi-harmonic Debye model to obtain thermodynamic properties of materials
// Based on original Fortran program written by M. A. Blanco et al.
// See Computer Physics Communications 158, 57-72 (2004) and Journal of Molecular Structure (Theochem) 368, 245-255 (1996) for details of original GIBBS program
// See C. Toher et al., Phys. Rev. B 90, 174107 (2014), Phys. Rev. 1, 015401 (2017) and references therein for description of this AGL implementation
// Please cite these works in addition to the general AFLOW papers if you use results generated using AGL
//

// *******************************************************************************
// The following functions are for generating _AFLOWIN_ files
// *******************************************************************************

// ***************************************************************************
// AGL_functions::aglvaspflags
// ***************************************************************************
namespace AGL_functions {
  //
  // Function to assign values for VASP input flags from _AFLOWIN_ file to vaspRun _xvasp class
  // Adapted from section of AFLOW APL function DirectMethodPC::runVASPCalculations()
  //
  uint aglvaspflags(_xvasp& vaspRun, _vflags& _vaspFlags, _kflags& _kbinFlags, string& runname, ofstream& FileMESSAGE) {
    ostringstream aus;
    // SOME WARNINGS: check existence of LOCK and OUTCAR.static files
    if( !aurostd::FileExist( vaspRun.Directory + "/"+_AFLOWLOCK_ ) &&
	aurostd::FileExist( vaspRun.Directory + string("/OUTCAR.static") ) ) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_WARNING_ + "found OUTCAR.static but no LOCK in " <<  vaspRun.Directory << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      return 1;
    }

    if( aurostd::FileExist( vaspRun.Directory + "/"+_AFLOWLOCK_ ) &&
	!aurostd::FileExist( vaspRun.Directory + string("/OUTCAR.static") ) ) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_WARNING_ + "found LOCK but no OUTCAR.static in " <<  vaspRun.Directory << endl;
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
    vaspRun.INCAR.str(std::string());
    string system;
    for(uint j=0; j < vaspRun.str.species.size(); j++)
      system = system + vaspRun.str.species_pp.at(j);
    system = system + "@" + runname;
    vaspRun.INCAR << "SYSTEM=" << system << std::endl;
    vaspRun.INCAR << "# Added by [AFLOW_GIBBS] begin" << std::endl;
    vaspRun.INCAR << "NELMIN=4         # The forces have to be well converged" << std::endl;
    vaspRun.INCAR << "NELM = 120       # Many electronic steps (SC2013)" << std::endl;
    vaspRun.INCAR << "ADDGRID=.TRUE.   # For finer forces" << std::endl;
    vaspRun.INCAR << "# Added by [AFLOW_GIBBS] end" << std::endl;

    // Change format of POSCAR
    if( ( !_kbinFlags.KBIN_MPI && ( _kbinFlags.KBIN_BIN.find("46") != string::npos ) ) ||
	(  _kbinFlags.KBIN_MPI && ( _kbinFlags.KBIN_MPI_BIN.find("46") != string::npos ) ) ) {
      vaspRun.str.is_vasp5_poscar_format = false; 
    }
    return 0;
  }
} // namespace AGL_functions

// ***************************************************************************
// AGL_functions::createAFLOWIN
// ***************************************************************************
namespace AGL_functions {
  //
  // Create _AFLOWIN_ file: makes new directory and writes _AFLOWIN_ for strained structure file inside it 
  // Adapted from that in AFLOW APL function PhononCalculator::createAFLOWIN()
  //
  uint createAFLOWIN(_xvasp& vaspRun, _xvasp& xvasp, _kflags& _kbinFlags, _vflags& _vaspFlags, _AGL_data& AGL_data, ofstream& FileMESSAGE) {
    bool AFLOWIN_QE_FLAG=FALSE;
    bool SPACES=FALSE;
    ostringstream aus;

    if( !aurostd::FileExist( vaspRun.Directory) ) {
      aurostd::DirectoryMake( vaspRun.Directory );
    }
    // CHMOD Directory 777: change directory permissions to read+write+execute for all users
    aurostd::DirectoryChmod("777", vaspRun.Directory);

    // Create file
    string filename =  vaspRun.Directory + "/"+_AFLOWIN_;
    ofstream outfile(filename.c_str(),ios_base::out);

    // Check _AFLOWIN_ file is open
    if( !outfile.is_open() ) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_WARNING_ + "Cannot create [" << _AFLOWIN_ << "] file" << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      return 1;
    }
    //CHMOD a+rw _AFLOWIN_: change permissions on _AFLOWIN_ file
    aurostd::ChmodFile("a+rw",filename);

    // Write to _AFLOWIN_ file
    if(SPACES) { outfile << std::endl; }
    outfile << AFLOWIN_SEPARATION_LINE << std::endl;
    outfile << "[AFLOW] _ ___ _" << std::endl;
    outfile << "[AFLOW] / \\|  || \\ |" << std::endl;
    outfile << "[AFLOW] | o | _ | " << std::endl;
    outfile << "[AFLOW] |_n_|__||___| automatic generated file" << std::endl;
    outfile << "[AFLOW]" << std::endl;
    outfile << AFLOWIN_SEPARATION_LINE << std::endl;
    if(SPACES) { outfile << std::endl; }
    outfile << "[AFLOW_MODE=VASP]" << std::endl;
    outfile << "[AFLOW_MODE_ZIP=" << _kbinFlags.KZIP_BIN << "]" << std::endl;
    if(SPACES) { outfile << std::endl; }

    // CO 180130 - START
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
    // CO 180130 - STOP

    // CO 180130 - making obsolete with lines above
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

    // Write INCAR lines to _AFLOWIN_ file 
    //
    // OBSOLETE outfile << AFLOWIN_SEPARATION_LINE << std::endl;
    // OBSOLETE if(SPACES) outfile << std::endl;
    // OBSOLETE outfile << "[VASP_RUN]STATIC" << std::endl;
    // OBSOLETE if(SPACES) outfile << std::endl;
    outfile << AFLOWIN_SEPARATION_LINE << std::endl;
    if(AGL_data.relax_static) {
      outfile << "[VASP_RUN]RELAX_STATIC=2" << std::endl;
    } else if(AGL_data.static_only) {
      outfile << "[VASP_RUN]STATIC" << std::endl;      
    } else {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_ERROR_ + "No run type selected" << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      return 1;
    }

    outfile << AFLOWIN_SEPARATION_LINE << std::endl;
    if(_vaspFlags.KBIN_VASP_FORCE_OPTION_LDAU1.isentry) {
      outfile << "[VASP_FORCE_OPTION]LDAU1= ON"  << std::endl;
      outfile << "[VASP_FORCE_OPTION]LDAU_PARAMETERS= " << _vaspFlags.KBIN_VASP_LDAU_PARAMETERS  << std::endl;
    }
    if(_vaspFlags.KBIN_VASP_FORCE_OPTION_LDAU2.isentry) {
      outfile << "[VASP_FORCE_OPTION]LDAU2=ON " << std::endl;
      outfile << "[VASP_FORCE_OPTION]LDAU_PARAMETERS= " <<  _vaspFlags.KBIN_VASP_LDAU_PARAMETERS  << std::endl;
    }
    if(SPACES) { outfile << std::endl; }
    outfile << "[VASP_FORCE_OPTION]RELAX_IONS" << std::endl;
    outfile << "[VASP_FORCE_OPTION]WAVECAR=OFF" << std::endl;
    if( _vaspFlags.KBIN_VASP_FORCE_OPTION_BADER.isentry &&_vaspFlags.KBIN_VASP_FORCE_OPTION_BADER.option) { 
      outfile << "[VASP_FORCE_OPTION]CHGCAR=ON" << std::endl; 
    } else {
      outfile << "[VASP_FORCE_OPTION]CHGCAR=OFF" << std::endl; 
    }
    // OBSOLETE outfile << "[VASP_FORCE_OPTION]CHGCAR=OFF" << std::endl;
    outfile << "[VASP_FORCE_OPTION]PREC=ACCURATE" << std::endl;
    outfile << "[VASP_FORCE_OPTION]ALGO=NORMAL" << std::endl;

    if( _vaspFlags.KBIN_VASP_FORCE_OPTION_AUTO_PSEUDOPOTENTIALS.isentry ) { outfile << "[VASP_FORCE_OPTION]AUTO_PSEUDOPOTENTIALS=" << _vaspFlags.KBIN_VASP_FORCE_OPTION_AUTO_PSEUDOPOTENTIALS.xscheme << std::endl; }
    if( _vaspFlags.KBIN_VASP_FORCE_OPTION_ABMIX.isentry ) { outfile << "[VASP_FORCE_OPTION]ABMIX=" << _vaspFlags.KBIN_VASP_FORCE_OPTION_ABMIX.xscheme << std::endl; }
    if( _vaspFlags.KBIN_VASP_FORCE_OPTION_TYPE.isentry ) { outfile << "[VASP_FORCE_OPTION]TYPE=" << _vaspFlags.KBIN_VASP_FORCE_OPTION_TYPE.xscheme << std::endl; }
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
    
    // Write KPOINTS related lines to _AFLOWIN_ file
    outfile << AFLOWIN_SEPARATION_LINE << std::endl;
    if(SPACES) { outfile << std::endl; }
    outfile << "[VASP_KPOINTS_MODE_IMPLICIT] " << std::endl;
    outfile << "[VASP_KPOINTS_FILE]KSCHEME=" << vaspRun.AVASP_KSCHEME << " " << std::endl;
    outfile << "[VASP_KPOINTS_FILE]KPPRA=" << vaspRun.AVASP_value_KPPRA << std::endl;
    outfile << "[VASP_KPOINTS_FILE]STATIC_KSCHEME=" << vaspRun.AVASP_STATIC_KSCHEME << " " << std::endl;
    outfile << "[VASP_KPOINTS_FILE]STATIC_KPPRA=" << vaspRun.AVASP_value_KPPRA_STATIC << std::endl;
    if(SPACES) { outfile << std::endl; }
    
    // Write POTCAR related lines to _AFLOWIN_ file
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
    // OBSOLETE     aus << _AGLSTR_MESSAGE_ << "Species_pp " << j << " = " << xvasp.str.species_pp.at(j) << endl;
    // OBSOLETE     aus << _AGLSTR_MESSAGE_ << "Species " << j << " = " << xvasp.str.species.at(j) << endl;
    // OBSOLETE   }
    // OBSOLETE   aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    // OBSOLETE   for(uint j=0; j < xvasp.str.species_pp.size(); j++) {
    // OBSOLETE  outfile << "[VASP_POTCAR_FILE]" << xvasp.str.species_pp.at(j) << std::endl;
    // OBSOLETE}

    if(SPACES) { outfile << std::endl; }

    // Write POSCAR lines to _AFLOWIN_ file
    outfile << AFLOWIN_SEPARATION_LINE << std::endl;
    if(SPACES) { outfile << std::endl; }
    outfile << "[VASP_POSCAR_MODE_EXPLICIT]START " << std::endl;
    vaspRun.str.is_vasp4_poscar_format=TRUE;
    vaspRun.str.is_vasp5_poscar_format=FALSE;
    outfile << vaspRun.str;
    outfile << "[VASP_POSCAR_MODE_EXPLICIT]STOP " << std::endl;
    if(SPACES) { outfile << std::endl; }

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
} // namespace AGL_functions

// ************************************************************************************************
// This set of functions extract, sort and check (E, V) data from VASP runs
// ************************************************************************************************

// ***************************************************************************
// AGL_functions::extractenerg
// ***************************************************************************
namespace AGL_functions {
  //
  // extractenerg: Extract final energies from the completed VASP calculations
  // Adapted from section of AFLOW APL function DirectMethodPC::runVASPCalculations()
  //
  uint extractenerg(vector<_xvasp>& vaspRuns, _AGL_data& AGL_data, ofstream& FileMESSAGE) {
    bool LVERBOSE=(FALSE || XHOST.DEBUG);
    ostringstream aus;
    vector<string> vfile, dfile;
    aurostd::string2tokens(string("OUTCAR.static.bz2,OUTCAR.static.gz,OUTCAR.static.xz,OUTCAR.static"),vfile,",");
    xOUTCAR outcar;
    string vfilename, dfilename, ffilename;
    bool skipdir = false;
    aurostd::xmatrix<double> stress_tensor(3, 3);
    AGL_data.volumeinput.clear();
    AGL_data.energyinput.clear();
    AGL_data.pressurecalculated.clear();
    AGL_data.stresscalculated.clear();
    for(uint idVaspRun = 0; idVaspRun < vaspRuns.size(); idVaspRun++) {
      skipdir = false;
      aurostd::StringstreamClean(aus);
      // Print out total energy
      aus << _AGLSTR_MESSAGE_ << "System number = " << idVaspRun << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      aurostd::string2tokens(vaspRuns.at(idVaspRun).Directory, dfile, "/");
      dfilename = dfile.at(dfile.size()-1);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ + "Directory name = " << dfilename << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      // Check if structure is on list of failed runs to be skipped
      // If so, then skip reading and continue to next structure
      for (uint ij = 0; ij < AGL_data.failed_arun_list.size(); ij++) {
	ffilename = AGL_data.failed_arun_list.at(ij);
	if(LVERBOSE) {
	  aurostd::StringstreamClean(aus);
	  aus << _AGLSTR_MESSAGE_ + "dfilename = " << dfilename << endl;
	  aus << _AGLSTR_MESSAGE_ + "ffilename = " << ffilename << endl;
	  aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
	}
	if(aurostd::substring2bool(dfilename,ffilename,TRUE)) {
	  aurostd::StringstreamClean(aus);
	  aus << _AGLSTR_MESSAGE_ + "Found directory in to-skip list: " << dfilename << endl;
	  aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
	  skipdir = true;
	}
      }
      if(skipdir) {
	aurostd::StringstreamClean(aus);
	aus << _AGLSTR_MESSAGE_ + "Skipping directory: " << dfilename << endl;
	aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
	skipdir = false;
	continue;
      }

      // If tarred and compressed directory exists...
      // [OBSOLETE]      string tarfilename = vaspRuns.at(idVaspRun).Directory + ".tar.bz2";
      // [OBSOLETE]    if( aurostd::FileExist(tarfilename) ) {
      // [OBSOLETE]   	// Extract all...
      // [OBSOLETE]   	aurostd::execute( string("tar -xf ") + tarfilename );
      // [OBSOLETE]     }
      if( aurostd::FileExist(vaspRuns.at(idVaspRun).Directory + ".tar.bz2") ) { aurostd::execute( string("tar -xf ") + vaspRuns.at(idVaspRun).Directory + ".tar.bz2" ); } // Extract all...
      if( aurostd::FileExist(vaspRuns.at(idVaspRun).Directory + ".tar.gz") ) { aurostd::execute( string("tar -xf ") + vaspRuns.at(idVaspRun).Directory + ".tar.gz" ); } // Extract all...
      if( aurostd::FileExist(vaspRuns.at(idVaspRun).Directory + ".tar.xz") ) { aurostd::execute( string("tar -xf ") + vaspRuns.at(idVaspRun).Directory + ".tar.xz" ); } // Extract all...
      
      // If the LOCK file is missing, then it is probably a corrupted run
      // Do not accept it and wait for the new run
      if( !aurostd::FileExist( vaspRuns.at(idVaspRun).Directory + "/"+_AFLOWLOCK_ ) ) {
	aurostd::StringstreamClean(aus);
	aus <<  _AGLSTR_WARNING_ + "The " << _AFLOWLOCK_ << " file in " << vaspRuns.at(idVaspRun).Directory << " directory is missing." << endl;
	aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
	if(AGL_data.autoskipfailedaruns) {
	  aurostd::StringstreamClean(aus);
	  aus << _AGLSTR_MESSAGE_ + "Skipping directory: " << dfilename << endl;
	  aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
	  continue;
	} else {
	  throw AGLStageBreak();
	}
      }

      //for(uint i=0;i<vfile.size()&&(outcar.outcar=="");i++) {
      for(uint i=0;i<vfile.size()&&(outcar.content=="");i++) {
	if(aurostd::FileExist(vaspRuns.at(idVaspRun).Directory+"/"+vfile.at(i))) {
	  outcar.GetPropertiesFile(vaspRuns.at(idVaspRun).Directory+"/"+vfile.at(i));
	}
      }
      //if(outcar.outcar=="") {
      if(outcar.content=="") {
	aurostd::StringstreamClean(aus);
	aus << _AGLSTR_WARNING_ + "The OUTCAR.static file in " << vaspRuns.at(idVaspRun).Directory << " directory is missing." << endl;
	aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
	if(AGL_data.autoskipfailedaruns) {
	  aurostd::StringstreamClean(aus);
	  aus << _AGLSTR_MESSAGE_ + "Skipping directory: " << dfilename << endl;
	  aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
	  continue;
	} else {
	  throw AGLStageBreak();
	}
      }
 
      AGL_data.energyinput.push_back(outcar.energy_cell);
      AGL_data.volumeinput.push_back(vaspRuns.at(idVaspRun).str.Volume());
      AGL_data.pressurecalculated.push_back(outcar.pressure_residual);
      stress_tensor = -outcar.stress;
      AGL_data.stresscalculated.push_back(stress_tensor);
      AGL_data.structurecalculated.push_back(idVaspRun);

      aurostd::StringstreamClean(aus);
      // Print out total energy, volume and calculated residual pressure
      // OBSOLETE aus << _AGLSTR_MESSAGE_ << "System number = " << idVaspRun << ", Energy (eV) = " << AGL_data.energyinput.at(idVaspRun) << endl;
      aus << _AGLSTR_MESSAGE_ << "System number = " << idVaspRun << ", Energy (eV) = " << AGL_data.energyinput.at(AGL_data.energyinput.size()-1) << endl;
      aus << _AGLSTR_MESSAGE_ << "System number = " << idVaspRun << ", Volume (Ang^3) = " << AGL_data.volumeinput.at(AGL_data.volumeinput.size()-1) << endl;
      aus << _AGLSTR_MESSAGE_ << "System number = " << idVaspRun << ", Pressure (kB) = " << AGL_data.pressurecalculated.at(AGL_data.pressurecalculated.size()-1) << endl;      
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      outcar.clear();
    }
    return 0;
  }
} // namespace AGL_functions

// ***************************************************************************
// AGL_functions::checkmin
// ***************************************************************************
namespace AGL_functions {
  //
  // checkmin: checks position of lowest energy in (E, V)
  // If lowest energy corresponds to smallest volume, sets cmerr = 1
  // If lowest energy corresponds to largest volume, sets cmerr = 2
  // Otherwise, sets cmerr = 0
  //
  uint checkmin(_AGL_data& AGL_data, int& cmerr, ofstream& FileMESSAGE) {
    ostringstream aus;
    vector<double> energy(AGL_data.energyinput.size()), volume(AGL_data.volumeinput.size());
    cmerr = 0;
    uint aglerror = 0;
  
    for (uint i = 0; i < AGL_data.energyinput.size(); i++) {
      energy.at(i) = AGL_data.energyinput.at(i);
      volume.at(i) = AGL_data.volumeinput.at(i);
    }
    // Sorts (E, V) data into order of increasing volume
    //AGL_functions::qcksort (vol, idx, 0, is-1, FileMESSAGE);
    aglerror = AGL_functions::qcksortev (volume, energy, FileMESSAGE);
    if(aglerror != 0) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_ERROR_ + "Failed to sort E(V)" << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      return aglerror;
    }
    // Find global minimum of energy data points
    double etref = energy.at(0);
    uint itmin = 0;
    for (uint i = 0; i < AGL_data.energyinput.size(); i++) {
      if(energy.at(i) < etref) {
	etref = energy.at(i);
	itmin = i;
      }
    }
    // Check that the minimum energy does not correspond to the largest or smallest volume
    // If it does, then this suggests that a larger or smaller volume may have a lower energy
    if(itmin == 0) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_WARNING_ + "Minimum energy is for smallest volume" << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      cmerr = 1;
      return aglerror;
    } else if(itmin == (AGL_data.energyinput.size() - 1) ) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_WARNING_ + "Minimum energy is for largest volume" << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      cmerr = 2;
      return aglerror;
    } else {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Energy minimum is contained in (E, V) data" << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      cmerr = 0;
      return aglerror;
    }
    energy.clear();
    volume.clear();
    return aglerror;
  }
} // namespace AGL_functions

// ***************************************************************************
// AGL_functions::checkconcav
// ***************************************************************************
namespace AGL_functions {
  //
  // checkconcav: checks concavity of (E, V) data
  // If data around global minimum of supplied (E, V) data is not concave, sets ccerr = 1
  // Otherwise, sets ccerr = 0
  //
  uint checkconcav(_AGL_data& AGL_data, int& ccerr, ofstream& FileMESSAGE) {
    ostringstream aus;
    uint j = 1;
    vector<double> energy(AGL_data.energyinput.size()), volume(AGL_data.volumeinput.size());
    ccerr = 0;
    uint aglerror = 0;
    for (uint i = 0; i < AGL_data.energyinput.size(); i++) {
      energy.at(i) = AGL_data.energyinput.at(i);
      volume.at(i) = AGL_data.volumeinput.at(i);
    }
    // Sort (E, V) data in order of increasing volume
    aglerror = AGL_functions::qcksortev (volume, energy, FileMESSAGE);
    if(aglerror != 0) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_ERROR_ + "Failed to sort E(V)" << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      return aglerror;
    }
    // Find global minimum of energy data points
    double etref = energy.at(0);
    uint itmin = 0;
    for (uint i = 0; i < AGL_data.energyinput.size(); i++) {
      if(energy.at(i) < etref) {
	etref = energy.at(i);
	itmin = i;
      }
    }
    // Finds first acceptable point (first point where E-V curve is concave)
    while ( ( (energy.at(j) - energy.at(j-1))/(volume.at(j) - volume.at(j-1)) >= (energy.at(j) - energy.at(j+1))/(volume.at(j) - volume.at(j+1)) ) && j < AGL_data.energyinput.size() - 2) {
      j = j + 1;
    }
    aurostd::StringstreamClean(aus);
    aus << _AGLSTR_MESSAGE_ << "Concavity check: j = " << j << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    // If only the last three points are concave, then the (E, V) data cannot be used for GIBBS
    // A warning is given and the signal is given to rerun the VASP calculations with more k-points
    if(j == AGL_data.energyinput.size() - 1) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_WARNING_ + "All points show convex patterns" << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      ccerr = 1;
      return aglerror;
    }
    // If the first accepted point is already past the global minimum, the (E, V) could cause problems for GIBBS
    // Gives a warning and sends the signal to rerun the VASP calculations with more k-points
    if(j >= itmin) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_WARNING_ + "First concave point is already passed the global minimum" << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      ccerr = 1;
      return aglerror;
    }
    j = j + 1;
    uint jtmax = 2;
    
    // Point j marks the last accepted point, i the new trial point
    for (uint i=j+1; i <= AGL_data.energyinput.size()-1; i++) {
      if( (energy.at(j) - energy.at(j-1))/(volume.at(j) - volume.at(j-1)) < (energy.at(j) - energy.at(i))/(volume.at(j) - volume.at(i)) ) {
	j = j + 1;
	jtmax = i;
      }
    }
    
    // If the global minimum lies outside of the range of the accepted points, then there will be problems using the (E, V) data for GIBBS
    // Gives a warning and then gives the signal to rerun the VASP calculations with more k-points
    // Problems with noise in the (E, V) data are often caused by insufficient K-points or basis set
    if(jtmax <= itmin) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_WARNING_ + "Global minimum of (E, V) data lies outside of initial range of points accepted for concavity" << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      ccerr = 1;
      return aglerror;
    }
    // If the data is concave and includes the global minimum, then data should be good for GIBBS
    // Gives success message and sends signal to continue to GIBBS method
    aurostd::StringstreamClean(aus);
    aus << _AGLSTR_MESSAGE_ << "(E, V) data is concave" << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    ccerr = 0;
    return aglerror;
  }
} // namespace AGL_functions

// ************************************************************************************************
//  This set of functions sorts the (E, V) data, fits it by a polynomial, and finds the minimum
// ************************************************************************************************

// ***************************************************************************
// AGL_functions::qcksortev
// ***************************************************************************
namespace AGL_functions {
  //
  // Quick-sort algorithm: sorts the elements of array in ascending order
  // Calls aurostd::quicksort to actually sort data
  // Returns vol and energ sorted in order of increasing vol
  //
  uint qcksortev(vector<double>& vol, vector<double>& energ, ofstream& FileMESSAGE) {
    int icheck = 0;
    ostringstream aus;
    // Check if data is already in correct order
    for (uint i = 0; i < (vol.size()-1); i++) {
      if(vol.at(i+1) < vol.at(i)) {
	icheck = 1;
      }
    }
    // If data is already in correct order, exits function without sorting
    if(icheck == 0) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ + "qcksort: Data is already arranged in increasing order" << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      return 0;
    }
    // If data is not in the correct order, calls quicksort2 from aurostd library to sort it
    // Sorts both vol and energ so that they are in the order of increasing vol
    aurostd::sort (vol, energ);
    return 0;
  }
} // namespace AGL_functions

// ***************************************************************************
// AGL_functions::qcksortevt
// ***************************************************************************
namespace AGL_functions {
  //
  // Quick-sort algorithm: sorts the elements of array in ascending order
  // Calls aurostd::quicksort to actually sort data
  // Returns vol, energ and tdebye sorted in order of increasing vol
  //
  uint qcksortevt(vector<double>& vol, vector<double>& energ, vector<double>& tdebye, ofstream& FileMESSAGE) {
    int icheck = 0;
    ostringstream aus;
    // Check if data is already in correct order
    for (uint i = 0; i < (vol.size()-1); i++) {
      if(vol.at(i+1) < vol.at(i)) {
	icheck = 1;
      }
    }
    // If data is already in correct order, exits function without sorting
    if(icheck == 0) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "qcksort: Data is already arranged in increasing order" << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      return 0;
    }
    // If data is not in the correct order, calls quicksort2 from aurostd library to sort it
    // Sorts both vol and energ so that they are in the order of increasing vol
    aurostd::sort (vol, energ, tdebye);
    return 0;
  }
} // namespace AGL_functions

// **************************************************************************
//  End of AFLOW AGL set-up, extract, sort and check (E, V) data
// **************************************************************************

#endif  // _AFLOW_AGL_GET_EV_CPP
