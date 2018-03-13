// ***************************************************************************
// *                                                                         *
// *              Aflow KEVIN RASCH - Duke University 2013                   *
// *                                                                         *
// ***************************************************************************
// aflow_bader.cpp
// functions written by KEVIN RASCH
// UPDATED by COREY OSES
// 2013: kevin.rasch@duke.edu
// 2015: corey.oses@duke.edu

#include "aflow.h"
#include "aflow_bader.h"
#include "APL/aflow_apl.h"
#include <iostream>
#include <string>

using namespace apl;

//***************************************************************************//
// bader_functions::BaderCalc
//***************************************************************************//
namespace bader_functions {
string BaderCalc(aurostd::xoption vpflow) {  //corey
  ostringstream oss;
  string soliloquy = "bader_functions::BaderCalc():  ";  // so you know who's speaking
  //debug
  bool LDEBUG = (FALSE || XHOST.DEBUG);
  if(LDEBUG) oss << soliloquy << "BEGIN" << endl;
  string usage_usage = "aflow --bader [bader_options] -D DIRECTORY";
  vector<string> usage_options;
  usage_options.push_back(usage_usage);
  usage_options.push_back("");
  usage_options.push_back("options:");
  usage_options.push_back("--usage");
  usage_options.push_back("--critical_points|--cp");
  usage_options.push_back("--calculate=|--calc=bader|voronoi");
  usage_options.push_back("--nocalculate=|nocalc=bader|voronoi");
  usage_options.push_back("--partition=|--part=neargrid|ongrid");
  usage_options.push_back("--refine_edge_method=|--rem=|--r=-1,-2,-3");
  usage_options.push_back("--reference=|--ref=CHGCAR FILE");
  usage_options.push_back("--vacuum=|--vac=off|auto|VACUUM DENSITY");
  usage_options.push_back("--terminate|--term=known|max");
  usage_options.push_back("--print_all=atom|bader");
  usage_options.push_back("--print_index=|--print_ind=atom|bader");
  usage_options.push_back("--print_select_atom=|--print_sel_atom=1,2,3...|1-5 [volume list (separated by commas) or range (separated by hyphen)]");
  usage_options.push_back("--print_select_bader=|--print_sel_bader=1,2,3...|1-5 [volume list (separated by commas) or range (separated by hyphen)]");
  usage_options.push_back("--print_sum_atom=1,2,3...|1-5 [volume list (separated by commas) or range (separated by hyphen)]");
  usage_options.push_back("--print_sum_bader=1,2,3...|1-5 [volume list (separated by commas) or range (separated by hyphen)]");
  usage_options.push_back("--quiet|--q");
  // output usage
  if(LDEBUG) oss << soliloquy << "CHECK USAGE" << endl;
  if(vpflow.flag("BADER::USAGE")) {
    init::ErrorOption(cout, "--usage", "pflow::BADER()", usage_options);
    return oss.str();
  }
  //FIRST STOP FOR CALLS MADE FROM COMMAND LINE
  string bader_options;  //to store bader options
  //directory
  string directory = ".";  //default
  if(XHOST.vflag_control.flag("DIRECTORY")) {
    directory = XHOST.vflag_control.getattachedscheme("DIRECTORY");  //XHOST.vflag_control.flag("DIRECTORY");//aflags.Directory;//vpflow.getattachedscheme("BADER::DIRECTORY");
  }
  FixDirectory(directory, oss);
  if(!aurostd::FileExist(directory)) {
    oss << endl;
    oss << soliloquy << "ERROR: Unable to locate " << directory << "." << endl;
    oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
    oss << soliloquy << "Exiting." << endl;
    oss << endl;
    return oss.str();
  }
  if(!Flags2BaderCommands(vpflow, bader_options, oss)) {
    oss << endl;
    oss << soliloquy << "ERROR: Unable to gather flags for bader command." << endl;
    oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
    oss << soliloquy << "Exiting." << endl;
    oss << endl;
    return oss.str();
  }
  if(!BaderCalc(vpflow, bader_options, directory, oss)) {
    oss << endl;
    oss << soliloquy << "ERROR: Unable to perform bader calculation as specified." << endl;
    oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
    oss << soliloquy << "Exiting." << endl;
    oss << endl;
    return oss.str();
  }
  if(vpflow.flag("BADER::QUIET")) { oss.str(""); }
  return oss.str();
}
}  // namespace bader_functions

//***************************************************************************//
// bader_functions::BaderCalc
//***************************************************************************//
namespace bader_functions {
bool BaderCalc(aurostd::xoption& vpflow,
               const string& bader_options,
               string& directory,
               ostream& oss) {
  //SECOND STOP FOR CALLS MADE FROM COMMAND_LINE
  string soliloquy = "bader_functions::BaderCalc():  ";          // so you know who's speaking
  string soliloquy_empty = "                                 ";  // so you know who's speaking
  //debug
  bool LDEBUG = (FALSE || XHOST.DEBUG);
  if(LDEBUG) oss << soliloquy << "BEGIN" << endl;

  // get species names and valence charges from OUTCAR
  if(LDEBUG) oss << soliloquy << "LOOKING FOR SUITABLE OUTCAR" << endl;
  string outcar_file;
  if(!BaderExtensionFound("OUTCAR", outcar_file, directory)) {  //replace in place
    oss << endl;
    oss << soliloquy << "ERROR: Cannot find suitable OUTCAR (OUTCAR.static, OUTCAR, or compressed variants)." << endl;
    oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
    oss << soliloquy << "Exiting." << endl;
    oss << endl;
    return FALSE;
  }
  oss << soliloquy << "using " << outcar_file << " to extract system name, species, and valence charges." << endl;
  xOUTCAR outcar(outcar_file);  //can handle compressed files, automatically pulls in system, species (strings), and vZVAL (deque<doubles>)
  //[OBSOLETE as it is automatic now]outcar.GetProperties();                       //we need SYSTEM

  // get number of each species and valence charges from POSCAR
  if(LDEBUG) oss << soliloquy << "LOOKING FOR SUITABLE POSCAR" << endl;
  string poscar_file;
  if(!BaderExtensionFound("POSCAR", poscar_file, directory)) {
    oss << endl;
    oss << soliloquy << "ERROR: Cannot find suitable POSCAR (POSCAR.static, POSCAR, or compressed variants)." << endl;
    oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
    oss << soliloquy << "Exiting." << endl;
    oss << endl;
    return FALSE;
  }
  oss << soliloquy << "using " << poscar_file << " to extract atom count for each species." << endl;
  xstructure xstr_bader(poscar_file, IOVASP_POSCAR);  //automatically pulls in num_each_type (ints); can handle compressed files

  // get cutoffs and downsample_ratios from FLAGS
  if(LDEBUG) oss << soliloquy << "CHECK FOR JVXL COMMAND" << endl;
  vector<double> cutoffs;
  vector<int> downsample_ratios;
  vpflow.flag("BADER::JVXL_CYCLIC", FALSE);
  if(vpflow.flag("BADER::JVXL_ALL_SPECIES")) {
    if(LDEBUG) oss << soliloquy << "vpflow.getattachedscheme(\"BADER::JVXL_ALL_SPECIES\")=" << vpflow.flag("BADER::JVXL_ALL_SPECIES") << endl;
    vector<string> sets, tokens;
    string input = vpflow.getattachedscheme("BADER::JVXL_ALL_SPECIES");

    // case 1: cyclic variable
    // 2a:  CUTOFF1,CUTOFF2…[::DOWNSAMPLE]
    // 2b:  CUTOFF[::DOWNSAMPLE1,DOWNSAMPLE2,...]
    // 2c:  some combination of above
    if(LDEBUG) oss << soliloquy << "CHECK INPUT TYPE" << endl;
    if(aurostd::substring2bool(input, "::")) {
      if(LDEBUG) oss << soliloquy << "CYCLIC FOUND" << endl;
      vpflow.flag("BADER::JVXL_CYCLIC", TRUE);
      aurostd::string2tokens(input, tokens, "::");
      if(tokens.size() < 1 || tokens.size() > 2) {
        oss << endl;
        oss << soliloquy << "ERROR: Incorrect format for input - number of tokens (::)." << endl;
        oss << soliloquy << "Needs to be either: " << endl;
        oss << soliloquy_empty << "CUTOFF1,CUTOFF2…[::DOWNSAMPLE], or " << endl;
        oss << soliloquy_empty << "CUTOFF[::DOWNSAMPLE1,DOWNSAMPLE2,...], or" << endl;
        oss << soliloquy_empty << "some combination." << endl;
        oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
        oss << soliloquy << "Exiting." << endl;
        oss << endl;
        return FALSE;
      }
      //make sure there's no mixing
      if(LDEBUG) oss << soliloquy << "CHECK THAT THERE'S NO MIXING BETWEEN CYCLIC AND SETS" << endl;
      for (uint i = 0; i < tokens.size(); i++) {
        if(aurostd::substring2bool(tokens.at(i), ":")) {
          oss << endl;
          oss << soliloquy << "ERROR: Incorrect format for input, cannot specify sets (:) and cyclic (::) parameters." << endl;
          oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
          oss << soliloquy << "Exiting." << endl;
          oss << endl;
          return FALSE;
        }
      }
      if(aurostd::substring2bool(tokens.at(0), ",")) {
        aurostd::string2tokens<double>(tokens.at(0), cutoffs, ",");
      } else {
        cutoffs.push_back(aurostd::string2utype<double>(tokens.at(0)));
      }
      if(tokens.size() == 2) {
        if(aurostd::substring2bool(tokens.at(1), ",")) {
          aurostd::string2tokens<int>(tokens.at(1), downsample_ratios, ",");
        } else {
          downsample_ratios.push_back(aurostd::string2utype<int>(tokens.at(1)));
        }
      }
      // case 2, sets, CUTOFF1[,DOWNSAMPLE1]:CUTOFF2[,DOWNSAMPLE2]:...
    } else if(aurostd::substring2bool(input, ":")) {
      if(LDEBUG) oss << soliloquy << "SETS FOUND" << endl;
      aurostd::string2tokens(input, sets, ":");
      for (uint i = 0; i < sets.size(); i++) {
        aurostd::string2tokens(sets.at(i), tokens, ",");
        if(tokens.size() < 1 || tokens.size() > 2) {
          oss << endl;
          oss << soliloquy << "ERROR: Incorrect format for input " << i + 1 << "." << endl;
          oss << soliloquy << "Needs to be CUTOFF1[,DOWNSAMPLE1]:CUTOFF2[,DOWNSAMPLE2]:..." << endl;
          oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
          oss << soliloquy << "Exiting." << endl;
          oss << endl;
          return FALSE;
        }
        cutoffs.push_back(aurostd::string2utype<double>(tokens.at(0)));
        if(tokens.size() == 2) { downsample_ratios.push_back(aurostd::string2utype<int>(tokens.at(1))); }
      }
      // case 3: single set
    } else {
      if(LDEBUG) oss << soliloquy << "SINGLE SET FOUND" << endl;
      if(aurostd::substring2bool(input, ",")) {
        aurostd::string2tokens(input, tokens, ",");
        if(tokens.size() != 2) {
          oss << endl;
          oss << soliloquy << "ERROR: Incorrect format for input - number of tokens." << endl;
          oss << soliloquy << "Needs to be CUTOFF[,DOWNSAMPLE]." << endl;
          oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
          oss << soliloquy << "Exiting." << endl;
          oss << endl;
          return FALSE;
        }
        cutoffs.push_back(aurostd::string2utype<double>(tokens.at(0)));
        downsample_ratios.push_back(aurostd::string2utype<int>(tokens.at(1)));
      } else {
        cutoffs.push_back(aurostd::string2utype<double>(input));
      }
    }
  }
  //check if extracted successfully
  if(outcar.SYSTEM.empty()) {
    oss << endl;
    oss << soliloquy << "ERROR: Unable to extract system name from " << outcar_file << "." << endl;
    oss << soliloquy << outcar_file << " may be corrupted." << endl;
    oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
    oss << soliloquy << "Exiting." << endl;
    oss << endl;
    return FALSE; //CO 180220
  }
  if(outcar.species.empty()) {
    oss << endl;
    oss << soliloquy << "ERROR: Unable to extract species types from " << outcar_file << "." << endl;
    oss << soliloquy << outcar_file << " may be corrupted." << endl;
    oss << soliloquy << "OUTCAR may be corrupted." << endl;
    oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
    oss << soliloquy << "Exiting." << endl;
    oss << endl;
    return FALSE; //CO 180220
  }
  if(outcar.vZVAL.empty()) {
    oss << endl;
    oss << soliloquy << "ERROR: Unable to extract valence charges from " << outcar_file << "." << endl;
    oss << soliloquy << outcar_file << " may be corrupted." << endl;
    oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
    oss << soliloquy << "Exiting." << endl;
    oss << endl;
    return FALSE; //CO 180220
  }
  if(xstr_bader.num_each_type.empty()) {
    oss << endl;
    oss << soliloquy << "ERROR: Unable to extract atom count for each species from " << outcar_file << "." << endl;
    oss << soliloquy << poscar_file << " may be corrupted." << endl;
    oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
    oss << soliloquy << "Exiting." << endl;
    oss << endl;
    return FALSE; //CO 180220
  }
  //get vector for vZVAL, we need it for VVdiff later
  //it needs to be as big as the number of atoms in cell
  vector<double> vZVAL;
  for (uint i = 0; i < outcar.vZVAL.size(); i++) {
    for (uint j = 0; j < (uint)xstr_bader.num_each_type.at(i); j++) {
      vZVAL.push_back(outcar.vZVAL.at(i));
    }
  }
  string system_name=KBIN::ExtractSystemName(directory);
  return BaderCalc(vpflow, bader_options, /*outcar.SYSTEM,*/ system_name, outcar.species, xstr_bader.num_each_type, vZVAL, cutoffs, downsample_ratios, directory, oss);
}
}  // namespace bader_functions

//***************************************************************************//
// bader_functions::BaderCalc
//***************************************************************************//
namespace bader_functions {
bool BaderCalc(aurostd::xoption& vpflow,
               const string& bader_options,
               const string& prototype,
               const deque<string>& vspecies,
               const deque<int>& num_each_type,
               const vector<double>& vZVAL,
               const vector<double>& cutoffs,
               const vector<int>& downsample_ratios,
               string& directory,
               ostream& oss) {
  //MAIN FUNCTION
  bool LDEBUG = (FALSE || XHOST.DEBUG);
  string soliloquy = "bader_functions::BaderCalc():  ";          // so you know who's speaking
  string soliloquy_empty = "                                 ";  // so you know who's speaking
  
  //CO 180220 moved up from below
  string execution_path = aurostd::execute2string("pwd");
  FixDirectory(directory, oss);
  oss << soliloquy << "working within " << directory << "." << endl;
  
  //debug
  if(LDEBUG) oss << soliloquy << "BEGIN" << endl;

  //COREY fast checks START - 170613
  if(vspecies.size()!=num_each_type.size()){
    oss << endl;
    oss << soliloquy << "ERROR: Input incorrect, vspecies.size()!=num_each_type.size() (" << vspecies.size() << "!=" << num_each_type.size() << ")" << endl;
    oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
    oss << soliloquy << "Exiting." << endl;
    oss << endl;
    return FALSE;
  }
  uint natoms=0;
  for(uint i=0;i<vspecies.size();i++){
    for(uint j=0;j<(uint)num_each_type[i];j++){
      natoms++;
    }
  }
  if(vZVAL.size()!=natoms){
    oss << endl;
    oss << soliloquy << "ERROR: Input incorrect, vZVAL.size()!=natoms (" << vZVAL.size() << "!=" << natoms << ")" << endl;
    oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
    oss << soliloquy << "Exiting." << endl;
    oss << endl;
    return FALSE;
  }
  //COREY fast checks END - 170613

  vector<string> bzip_files, remove_files, move_files;
  string bzip_file;
  //bzip_files=files to bzip later
  //move_files=files to move later
  //remove_files=files to remove, include directory

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //
  //  Start Directory
  //
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // translate "./" into a full path, might be done twice
  //[CO 180220 - moved up]string execution_path = aurostd::execute2string("pwd");
  //[CO 180220 - moved up]FixDirectory(directory, oss);

  //NO LONGER AN ISSUE, corey
  //FORTRAN has a problem with strings/paths longer than 128 characters, check
  //if(directory.length()>=128) {
  //  oss << endl;
  //  oss << soliloquy << "ERROR: Directory string length longer than 128. This will cause issues in Fortran." << endl;
  //  oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
  //  oss << soliloquy << "Exiting." << endl;
  //  oss << endl;
  //  return FALSE;
  //} else if(directory.length()>100) {
  //  oss << soliloquy << "WARNING: Directory string length longer than 100. This may cause issues with Fortran." << endl;
  //}

  // announce the directory to work in
  // unfortunate location for this message if you are coming from command line, but this is central to both command line/library
  // so I need it here
  //[CO 180220 - moved up]oss << soliloquy << "working within " << directory << "." << endl;

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //
  //  End Directory
  //
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // store expected name of the VASP data files required for Bader analysis
  // in a string

  //populate required_files
  if(LDEBUG) oss << soliloquy << "GATHERING REQUIRED FILES" << endl;
  vector<string> required_files;  //set up vector with TEMPORARY_DECOMPRESSED_FILES that we move later, and a files to move container
  oss << soliloquy << "adding CHGCAR to required files." << endl;
  required_files.push_back("CHGCAR");  //might need to decompress
  //DEFAULT, take AECCAR0+AECCAR2=REFERENCE_FILE (aflow.CHGCAR_sum), otherwise get from FLAGS
  //get REFERENCE_FILE if specified, otherwise DEFAULT
  if(LDEBUG) oss << soliloquy << "CHECK REFERENCE FILE" << endl;
  string ref_file, tmp_file;
  if(vpflow.flag("BADER::REFERENCE")) {
    if(LDEBUG) oss << soliloquy << "vpflow.getattachedscheme(\"BADER::REFERENCE\")=" << vpflow.flag("BADER::REFERENCE") << endl;
    ref_file = execution_path + "/" + vpflow.getattachedscheme("BADER::REFERENCE");  //might need to decompress
    oss << soliloquy << "looking to use " << ref_file << " as reference file." << endl;
    if((!aurostd::FileExist(ref_file)) || aurostd::FileEmpty(ref_file)) {
      oss << endl;
      oss << soliloquy << "ERROR: The reference charge file " << ref_file << " isn't present/is empty." << endl;
      oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
      oss << soliloquy << "Exiting." << endl;
      oss << endl;
      aurostd::RemoveFiles(remove_files);  //CLEAN UP BEFORE RETURNING FALSE
      return FALSE;
    }
    if(aurostd::IsCompressed(ref_file)) {
      oss << soliloquy << "decompressing " << ref_file << "." << endl;
      if(!aurostd::efile2tempfile(ref_file, ref_file)) {
        oss << endl;
        oss << soliloquy << "ERROR: Unable to decompress " << ref_file << "." << endl;
        oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
        oss << soliloquy << "Exiting." << endl;
        oss << endl;
        aurostd::RemoveFiles(remove_files);  //CLEAN UP BEFORE RETURNING FALSE
        return FALSE;
      }
      oss << soliloquy << "created temporary file " << ref_file << "." << endl;
      remove_files.push_back(ref_file);
    }
    oss << soliloquy << "referencing " << ref_file << "." << endl;
  } else {
    required_files.push_back("AECCAR0");
    required_files.push_back("AECCAR2");
    oss << soliloquy << "looking to use AECCAR0 and AECCAR2 to generate reference file." << endl;
  }

  if(LDEBUG) oss << soliloquy << "SEARCHING FOR / DECOMPRESSING REQUIRED FILES" << endl;
  for (uint i = 0; i < required_files.size(); i++) {
    if(!BaderExtensionFound(required_files.at(i), required_files.at(i), directory)) {  //replace in place
      oss << endl;
      oss << soliloquy << "ERROR: Cannot find required file " << required_files.at(i) << " (";
      oss << required_files.at(i) << ".static, " << required_files.at(i) << ", or compressed variants)." << endl;
      oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
      oss << soliloquy << "Exiting." << endl;
      oss << endl;
      aurostd::RemoveFiles(remove_files);  //CLEAN UP BEFORE RETURNING FALSE
      return FALSE;
    }
    if(LDEBUG) oss << soliloquy << "FOUND " << required_files.at(i) << endl;
    if(aurostd::IsCompressed(required_files.at(i))) {
      oss << soliloquy << "decompressing " << required_files.at(i) << "." << endl;
      if(!aurostd::efile2tempfile(required_files.at(i), tmp_file)) {
        oss << endl;
        oss << soliloquy << "ERROR: Unable to decompress " << required_files.at(i) << "." << endl;
        oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
        oss << soliloquy << "Exiting." << endl;
        oss << endl;
        aurostd::RemoveFiles(remove_files);  //CLEAN UP BEFORE RETURNING FALSE
        return FALSE;
      }
      oss << soliloquy << "created temporary file " << tmp_file << "." << endl;
      remove_files.push_back(tmp_file);
      if(aurostd::FileEmpty(tmp_file)) {
        oss << endl;
        oss << soliloquy << "ERROR: " << required_files.at(i) << " is empty." << endl;
        oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
        oss << soliloquy << "Exiting." << endl;
        oss << endl;
        aurostd::RemoveFiles(remove_files);  //CLEAN UP BEFORE RETURNING FALSE
        return FALSE;
      }
      required_files.at(i) = tmp_file;
    }
  }
  if(LDEBUG) oss << soliloquy << "FOUND ALL REQUIRED FILES" << endl;

  /*// Make sub-directory to hold output Bader files
      string bader_directory=directory + "/BADER";
      oss << soliloquy << "creating BADER sub-directory." << endl;
      if(!aurostd::DirectoryMake(bader_directory)) {
      oss << endl;
      oss << soliloquy << "ERROR: Could not open sub-directory in " << directory << "." << endl;
      oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
      oss << soliloquy << "Exiting." << endl;
      oss << endl;
      aurostd::RemoveFiles(remove_files);   //CLEAN UP BEFORE RETURNING FALSE
      return FALSE;
      }
      else {
      oss << soliloquy << "Saving Bader analysis data in " << bader_directory << "." << endl;
      }*/

  //first get species_header
  string species_header;
  for (uint i = 0; i < vspecies.size(); i++) {
    species_header.append(vspecies.at(i));
    if(i < vspecies.size() - 1) { species_header.append(" "); }
  }
  //species_header=species_header.substr(0,species_header.size()-1);  //FIX ME

  if(!vpflow.flag("BADER::REFERENCE")) {
    if(LDEBUG) oss << soliloquy << "BUILDING REFERENCE FILE" << endl;
    // build CHGCAR_sum to use as a reference when running Bader
    //ref_file=bader_directory + "/aflow.CHGCAR_sum";  // path and name for the charge density file
    ref_file = directory + "/aflow.CHGCAR_sum";  // path and name for the charge density file
    oss << soliloquy << "summing core and valence charge densities for reference file." << endl;
    if(!pflow::CHGSUM(species_header, required_files.at(1), required_files.at(2), ref_file, oss)) {
      oss << endl;
      oss << soliloquy << "ERROR: Unable to sum core and valence charge densities." << endl;
      oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
      oss << soliloquy << "Exiting." << endl;
      oss << endl;
      aurostd::RemoveFiles(remove_files);  //CLEAN UP BEFORE RETURNING FALSE
      return FALSE;
    }
    // check files exist
    if((!aurostd::FileExist(ref_file)) || aurostd::FileEmpty(ref_file)) {
      oss << endl;
      oss << soliloquy << "ERROR: The reference charge file " << ref_file << " isn't present/is empty." << endl;
      oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
      oss << soliloquy << "Exiting." << endl;
      oss << endl;
      aurostd::RemoveFiles(remove_files);  //CLEAN UP BEFORE RETURNING FALSE
      return FALSE;
    }
    bzip_files.push_back("aflow.CHGCAR_sum");
    oss << soliloquy << "referencing " << ref_file << "." << endl;
  }

  // run the bader code
  if(vpflow.flag("BADER::AFLOWLIB_LIBRARY")) {
    oss << soliloquy << "calling BADER code." << endl;
  } else {
    oss << soliloquy << "calling \"GRID BASED BADER ANALYSIS\" with" << endl;
    oss << soliloquy_empty << "bader " << bader_options << endl;
    oss << soliloquy_empty << required_files.at(0) << endl;  //CHGCAR file (decompressed)
    oss << soliloquy_empty << "-ref " << ref_file << "." << endl;
  }

  if(LDEBUG) oss << soliloquy << "CHECK BADER COMMAND" << endl;
  if(!aurostd::IsCommandAvailable("bader")) {
    oss << endl;
    oss << soliloquy << "ERROR: bader command NOT found." << endl;
    oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
    oss << soliloquy << "Exiting." << endl;
    oss << endl;
    aurostd::RemoveFiles(remove_files);  //CLEAN UP BEFORE RETURNING FALSE
    return FALSE;
  }
  string BADER_RUN_DIRECTORY = aurostd::TmpDirectoryCreate("BADER");
  oss << soliloquy << "changing directories to " << BADER_RUN_DIRECTORY << " to execute bader command." << endl;
  string work_dir=aurostd::execute2string("pwd");
  chdir(BADER_RUN_DIRECTORY.c_str());
  //FORTRAN has a problem with strings/paths longer than 128 characters, so make everything local
  aurostd::CopyFile(required_files.at(0), "CHGCAR");
  aurostd::CopyFile(ref_file, "REFERENCE");
  stringstream BADER_CODE_OUTPUT;
  //do critical points first, so normal ACF, AVF, etc. get overwritten by the usual bader command (with bader CHGCAR -ref SUM)
  if(vpflow.flag("BADER::CRITICAL_POINTS")) {
    BADER_CODE_OUTPUT << endl
                      << aurostd::execute2string(XHOST.command("bader") + " -cp REFERENCE") << endl;  //local variant, must be SUM_CHARGE file (sum of valence+core charges)
    //not sure if we expect a particular type of error
    //I'm guessing since the didn't exit for maxima in edge refinement (original bader command), we're okay
    if(aurostd::substring2bool(BADER_CODE_OUTPUT.str(), "ERROR: should be no new maxima in edge refinement")) {
      oss << endl;
      oss << soliloquy << "ERROR: bader command issued error" << endl;
      oss << soliloquy << "BADER ERROR: should be no new maxima in edge refinement" << endl;
      oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
      oss << soliloquy << "Exiting." << endl;
      oss << endl;
      aurostd::RemoveFiles(remove_files);  //CLEAN UP BEFORE RETURNING FALSE
      oss << soliloquy << "changing directories back to original directory " << directory << "." << endl;
      chdir(work_dir.c_str());//directory.c_str());
      oss << soliloquy << "removing temporary directory " << BADER_RUN_DIRECTORY << "." << endl;
      if(!aurostd::RemoveDirectory(BADER_RUN_DIRECTORY)) {  //CLEAN UP BEFORE RETURNING FALSE
        oss << endl;
        oss << soliloquy << "WARNING: issues deleting directory " << BADER_RUN_DIRECTORY << "." << endl;
        oss << soliloquy << "Please check." << endl;
        oss << endl;
      }
      //BZIP aflow.CHGCAR_sum (only file created), and it's created in the right directory, so no need to move
      for (uint i = 0; i < bzip_files.size(); i++) {
        //if(aurostd::FileExist(bader_directory+"/"+bzip_files.at(i))) {
        bzip_file = directory + "/" + bzip_files.at(i);
        if(aurostd::FileExist(bzip_file)) {
          oss << soliloquy << "compressing " << bzip_files.at(i) << "." << endl;
          aurostd::BzipFile(bzip_file);
        }
      }
      return FALSE;
    }
  }
  // we don't want output in library runs
  if(!vpflow.flag("BADER::AFLOWLIB_LIBRARY")) {
    oss << BADER_CODE_OUTPUT.str();
  }
  BADER_CODE_OUTPUT.str("");

  //BADER_CODE_OUTPUT << endl << aurostd::execute2string(XHOST.command("bader")+" "+bader_options+" "+required_files.at(0)+" -ref "+ref_file) << endl;
  BADER_CODE_OUTPUT << endl
                    << aurostd::execute2string(XHOST.command("bader") + " " + bader_options + " CHGCAR -ref REFERENCE") << endl;  //local variant
  //BADER_CODE_OUTPUT << endl << aurostd::execute2string("/home/corey/bin/bader "+bader_options+" "+required_files.at(0)+" -ref "+ref_file) << endl;
  if(aurostd::substring2bool(BADER_CODE_OUTPUT.str(), "ERROR: should be no new maxima in edge refinement")) {
    oss << endl;
    oss << soliloquy << "ERROR: bader command issued error" << endl;
    oss << soliloquy << "BADER ERROR: should be no new maxima in edge refinement" << endl;
    oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
    oss << soliloquy << "Exiting." << endl;
    oss << endl;
    aurostd::RemoveFiles(remove_files);  //CLEAN UP BEFORE RETURNING FALSE
    oss << soliloquy << "changing directories back to original directory " << directory << "." << endl;
    chdir(work_dir.c_str());//directory.c_str());
    oss << soliloquy << "removing temporary directory " << BADER_RUN_DIRECTORY << "." << endl;
    if(!aurostd::RemoveDirectory(BADER_RUN_DIRECTORY)) {  //CLEAN UP BEFORE RETURNING FALSE
      oss << endl;
      oss << soliloquy << "WARNING: issues deleting directory " << BADER_RUN_DIRECTORY << "." << endl;
      oss << soliloquy << "Please check." << endl;
      oss << endl;
    }
    //BZIP aflow.CHGCAR_sum (only file created), and it's created in the right directory, so no need to move
    for (uint i = 0; i < bzip_files.size(); i++) {
      //if(aurostd::FileExist(bader_directory+"/"+bzip_files.at(i))) {
      bzip_file = directory + "/" + bzip_files.at(i);
      if(aurostd::FileExist(bzip_file)) {
        oss << soliloquy << "compressing " << bzip_files.at(i) << "." << endl;
        aurostd::BzipFile(bzip_file);
      }
    }
    return FALSE;
  }
  // we don't want output in library runs
  if(!vpflow.flag("BADER::AFLOWLIB_LIBRARY")) {
    oss << BADER_CODE_OUTPUT.str();
  }
  BADER_CODE_OUTPUT.str("");

  //REMOVE FILES, round 1
  if(LDEBUG) oss << soliloquy << "REMOVING FILES" << endl;
  for (uint i = 0; i < remove_files.size(); i++) {
    oss << soliloquy << "removing " << remove_files.at(i) << "." << endl;
    aurostd::RemoveFile(remove_files.at(i));
  }
  remove_files.clear();

  //no bader files expected, dumb case but it's there
  if(vpflow.flag("BADER::NOCALCULATE") && vpflow.getattachedscheme("BADER::NOCALCULATE") == "bader") {
    oss << soliloquy << "no bader files expected." << endl;
    aurostd::RemoveFiles(remove_files);  //CLEAN UP BEFORE RETURNING FALSE
    oss << soliloquy << "changing directories back to original directory " << directory << "." << endl;
    chdir(work_dir.c_str());//directory.c_str());
    oss << soliloquy << "removing temporary directory " << BADER_RUN_DIRECTORY << "." << endl;
    if(!aurostd::RemoveDirectory(BADER_RUN_DIRECTORY)) {
      oss << endl;
      oss << soliloquy << "WARNING: issues deleting directory " << BADER_RUN_DIRECTORY << "." << endl;
      oss << soliloquy << "Please check." << endl;
      oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
      oss << endl;
    }
    //BZIP aflow.CHGCAR_sum (only file created), and it's created in the right directory, so no need to move
    for (uint i = 0; i < bzip_files.size(); i++) {
      //if(aurostd::FileExist(bader_directory+"/"+bzip_files.at(i))) {
      bzip_file = directory + "/" + bzip_files.at(i);
      if(aurostd::FileExist(bzip_file)) {
        oss << soliloquy << "compressing " << bzip_files.at(i) << "." << endl;
        aurostd::BzipFile(bzip_file);
      }
    }
    oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
    //oss << soliloquy << "Exiting." << endl;
    return TRUE;
  }

  //GET LIST OF FILES CREATED TO MOVE AND BZIP
  if(LDEBUG) oss << soliloquy << "GET LIST OF FILES CREATED TO MOVE AND BZIP" << endl;
  if(LDEBUG) oss << soliloquy << "CHECK STANDARD BADER FILES" << endl;
  vector<string> standard_bader_files;
  standard_bader_files.push_back("ACF.dat");
  standard_bader_files.push_back("AVF.dat");
  standard_bader_files.push_back("BCF.dat");
  for (uint i = 0; i < standard_bader_files.size(); i++) {
    //created in current directory
    if(aurostd::FileExist(standard_bader_files.at(i)) && aurostd::FileNotEmpty(standard_bader_files.at(i))) {
      bzip_files.push_back(standard_bader_files.at(i));
      move_files.push_back(standard_bader_files.at(i));
      oss << soliloquy << "file created: " << standard_bader_files.at(i) << "." << endl;
    } else {
      //echo warning, but don't exist
      oss << endl;
      oss << soliloquy << "WARNING: " << standard_bader_files.at(i) << " was not created/is empty." << endl;
      oss << soliloquy << "This is VERY unusual and likely an issue with the Bader analysis." << endl;
      oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
      oss << endl;
    }
  }
  //critical points file
  if(vpflow.flag("BADER::CRITICAL_POINTS")) {
    if(LDEBUG) oss << soliloquy << "CHECK CRITICAL POINTS FILE" << endl;
    if(aurostd::FileExist("CPF.dat") && aurostd::FileNotEmpty("CPF.dat")) {
      bzip_files.push_back("CPF.dat");
      move_files.push_back("CPF.dat");
      oss << soliloquy << "file created: CPF.dat." << endl;
    } else {
      //echo warning, but don't exist
      oss << endl;
      oss << soliloquy << "WARNING: CPF.dat was not created/is empty." << endl;
      oss << soliloquy << "This is VERY unusual and likely an issue with the Bader analysis." << endl;
      oss << endl;
    }
  }
  vector<int> range;
  int atomCOUNT = 1;
  string bader_file;
  //we don't want duplicate files, and print_all_atom and print_select_atom will produce duplicates
  //we also don't want to REedit files
  if(vpflow.flag("BADER::PRINT_ALL") || vpflow.flag("BADER::PRINT_SELECT_ATOM") || vpflow.flag("BADER::PRINT_SELECT_BADER")) {
    //BvAtxxxx.dat
    //know this from vpspecies and num_each_type
    if(vpflow.getattachedscheme("BADER::PRINT_ALL") == "atom" || vpflow.getattachedscheme("BADER::PRINT_ALL") == "both") {
      if(LDEBUG) oss << soliloquy << "CHECKING PRINT ALL ATOM FILES" << endl;
      atomCOUNT = 1;
      for (uint i = 0; i < vspecies.size(); i++) {
        for (uint speciesCOUNT = 0; speciesCOUNT < (uint)num_each_type.at(i); speciesCOUNT++) {
          bader_file = "BvAt" + aurostd::PaddedNumString(atomCOUNT, 4) + ".dat";
          if(aurostd::FileExist(bader_file) && aurostd::FileNotEmpty(bader_file)) {
            bzip_files.push_back(bader_file);
            move_files.push_back(bader_file);
            oss << soliloquy << "file created: " << bader_file << "." << endl;
            //if we keep these files, then edit their header
            if(!vpflow.flag("BADER::REMOVE_BADER_ATOMS")) {
              if(!prepare_CHGCAR_4_Jmol(bader_file, species_header, FALSE, oss)) {
                oss << endl;
                oss << soliloquy << "WARNING: Unable to fix header for " << bader_file << "." << endl;
                oss << soliloquy << "Please check." << endl;
                oss << endl;
              }
            }
          } else {
            //echo warning, do not exit
            oss << endl;
            oss << soliloquy << "WARNING: " << bader_file << " was not created/is empty." << endl;
            oss << soliloquy << "This is VERY unusual." << endl;
            oss << endl;
          }
          atomCOUNT++;
        }
      }
      //LIST/RANGE
    } else if(vpflow.flag("BADER::PRINT_SELECT_ATOM")) {
      if(LDEBUG) oss << soliloquy << "CHECKING SELECT ATOM FILES" << endl;
      if(!listORrange2vec(vpflow.getattachedscheme("BADER::PRINT_SELECT_ATOM"), range, oss)) {
        oss << endl;
        oss << soliloquy << "ERROR: Unable to extract list/range command for print select atom." << endl;
        oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
        oss << soliloquy << "Exiting." << endl;
        oss << endl;
        aurostd::RemoveFiles(remove_files);  //CLEAN UP BEFORE RETURNING FALSE
        oss << soliloquy << "changing directories back to original directory " << directory << "." << endl;
        chdir(work_dir.c_str());//directory.c_str());
        oss << soliloquy << "removing temporary directory " << BADER_RUN_DIRECTORY << "." << endl;
        if(!aurostd::RemoveDirectory(BADER_RUN_DIRECTORY)) {  //CLEAN UP BEFORE RETURNING FALSE
          oss << endl;
          oss << soliloquy << "WARNING: issues deleting directory " << BADER_RUN_DIRECTORY << "." << endl;
          oss << soliloquy << "Please check." << endl;
          oss << endl;
        }
        //BZIP aflow.CHGCAR_sum (only file created), and it's created in the right directory, so no need to move
        for (uint i = 0; i < bzip_files.size(); i++) {
          //if(aurostd::FileExist(bader_directory+"/"+bzip_files.at(i))) {
          bzip_file = directory + "/" + bzip_files.at(i);
          if(aurostd::FileExist(bzip_file)) {
            oss << soliloquy << "compressing " << bzip_files.at(i) << "." << endl;
            aurostd::BzipFile(bzip_file);
          }
        }
        return FALSE;
      }
      for (uint i = 0; i < range.size(); i++) {
        bader_file = "BvAt" + aurostd::PaddedNumString(range.at(i), 4) + ".dat";
        if(aurostd::FileExist(bader_file) && aurostd::FileNotEmpty(bader_file)) {
          bzip_files.push_back(bader_file);
          move_files.push_back(bader_file);
          oss << soliloquy << "file created: " << bader_file << "." << endl;
          if(!prepare_CHGCAR_4_Jmol(bader_file, species_header, FALSE, oss)) {
            oss << endl;
            oss << soliloquy << "WARNING: Unable to fix header for " << bader_file << "." << endl;
            oss << soliloquy << "Please check." << endl;
            oss << endl;
          }
        } else {
          //echo warning, do not exit
          oss << endl;
          oss << soliloquy << "WARNING: " << bader_file << " was not created/is empty." << endl;
          oss << soliloquy << "This is VERY unusual and likely an issue with the Bader analysis." << endl;
          oss << endl;
        }
      }
    }
    //Bvolxxxx.dat
    //no way to know unless we check what bader spits out
    if(vpflow.getattachedscheme("BADER::PRINT_ALL") == "bader" || vpflow.getattachedscheme("BADER::PRINT_ALL") == "both") {
      if(LDEBUG) oss << soliloquy << "CHECKING PRINT ALL BADER FILES" << endl;
      atomCOUNT = 1;
      //check that first exists, other throw an error
      bader_file = "Bvol" + aurostd::PaddedNumString(atomCOUNT, 4) + ".dat";
      if((!aurostd::FileExist(bader_file)) || aurostd::FileEmpty(bader_file)) {
        //echo warning, but do not exist
        oss << endl;
        oss << soliloquy << "WARNING: First expected Bvolxxxx.dat file NOT FOUND/IS EMPTY." << endl;
        oss << soliloquy << "This is VERY unusual and likely an issue with the Bader analysis." << endl;
        oss << endl;
      }
      while (aurostd::FileExist(bader_file) && aurostd::FileNotEmpty(bader_file)) {
        bzip_files.push_back(bader_file);
        move_files.push_back(bader_file);
        oss << soliloquy << "file created: " << bader_file << "." << endl;
        if(!prepare_CHGCAR_4_Jmol(bader_file, species_header, FALSE, oss)) {
          oss << endl;
          oss << soliloquy << "WARNING: Unable to fix header for " << bader_file << "." << endl;
          oss << soliloquy << "Please check." << endl;
          oss << endl;
        }
        bader_file = "Bvol" + aurostd::PaddedNumString(atomCOUNT++, 4) + ".dat";
      }
      //B_wexxxx.dat
      //these have different file names than Bvolxxx.dat even though they are the same files (artifact of bader code?)
      //LIST/RANGE
    } else if(vpflow.flag("BADER::PRINT_SELECT_BADER")) {
      if(LDEBUG) oss << soliloquy << "CHECKING PRINT SELECT BADER FILES" << endl;
      if(!listORrange2vec(vpflow.getattachedscheme("BADER::PRINT_SELECT_BADER"), range, oss)) {
        oss << endl;
        oss << soliloquy << "ERROR: Unable to extract list/range command for print select bader." << endl;
        oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
        oss << soliloquy << "Exiting." << endl;
        oss << endl;
        aurostd::RemoveFiles(remove_files);  //CLEAN UP BEFORE RETURNING FALSE
        oss << soliloquy << "changing directories back to original directory " << directory << "." << endl;
        chdir(work_dir.c_str());//directory.c_str());
        oss << soliloquy << "removing temporary directory " << BADER_RUN_DIRECTORY << "." << endl;
        if(!aurostd::RemoveDirectory(BADER_RUN_DIRECTORY)) {  //CLEAN UP BEFORE RETURNING FALSE
          oss << endl;
          oss << soliloquy << "WARNING: issues deleting directory " << BADER_RUN_DIRECTORY << "." << endl;
          oss << soliloquy << "Please check." << endl;
          oss << endl;
        }
        //BZIP aflow.CHGCAR_sum (only file created), and it's created in the right directory, so no need to move
        for (uint i = 0; i < bzip_files.size(); i++) {
          //if(aurostd::FileExist(bader_directory+"/"+bzip_files.at(i))) {
          bzip_file = directory + "/" + bzip_files.at(i);
          if(aurostd::FileExist(bzip_file)) {
            oss << soliloquy << "compressing " << bzip_files.at(i) << "." << endl;
            aurostd::BzipFile(bzip_file);
          }
        }
        return FALSE;
      }
      for (uint i = 0; i < range.size(); i++) {
        bader_file = "B_we" + aurostd::PaddedNumString(range.at(i), 4) + ".dat";
        if(aurostd::FileExist(bader_file) && aurostd::FileNotEmpty(bader_file)) {
          bzip_files.push_back(bader_file);
          move_files.push_back(bader_file);
          oss << soliloquy << "file created: " << bader_file << "." << endl;
          if(!prepare_CHGCAR_4_Jmol(bader_file, species_header, FALSE, oss)) {
            oss << endl;
            oss << soliloquy << "WARNING: Unable to fix header for " << bader_file << "." << endl;
            oss << soliloquy << "Please check." << endl;
            oss << endl;
          }
        } else {
          //echo warning, do not exit
          oss << endl;
          oss << soliloquy << "WARNING: " << bader_file << " was not created/is empty." << endl;
          oss << soliloquy << "This is VERY unusual and likely an issue with the Bader analysis." << endl;
          oss << endl;
        }
      }
    }
  }
  //BvAt_summed.dat
  if(vpflow.flag("BADER::PRINT_SUM_ATOM")) {
    if(LDEBUG) oss << soliloquy << "CHECKING PRINT SUM ATOM FILE" << endl;
    bader_file = "BvAt_summed.dat";
    if(aurostd::FileExist(bader_file) && aurostd::FileNotEmpty(bader_file)) {
      bzip_files.push_back(bader_file);
      move_files.push_back(bader_file);
      oss << soliloquy << "file created: " << bader_file << "." << endl;
      if(!prepare_CHGCAR_4_Jmol(bader_file, species_header, FALSE, oss)) {
        oss << endl;
        oss << soliloquy << "WARNING: Unable to fix header for " << bader_file << "." << endl;
        oss << soliloquy << "Please check." << endl;
        oss << endl;
      }
    } else {
      //echo warning, do not exit
      oss << endl;
      oss << soliloquy << "WARNING: " << bader_file << " was not created/is empty." << endl;
      oss << soliloquy << "This is VERY unusual and likely an issue with the Bader analysis." << endl;
      oss << endl;
    }
  }
  //Bvol_summed.dat
  if(vpflow.flag("BADER::PRINT_SUM_BADER")) {
    if(LDEBUG) oss << soliloquy << "CHECKING PRINT SUM BADER FILE" << endl;
    bader_file = "Bvol_summed.dat";
    if(aurostd::FileExist(bader_file) && aurostd::FileNotEmpty(bader_file)) {
      bzip_files.push_back(bader_file);
      move_files.push_back(bader_file);
      oss << soliloquy << "file created: " << bader_file << "." << endl;
      if(!prepare_CHGCAR_4_Jmol(bader_file, species_header, FALSE, oss)) {
        oss << endl;
        oss << soliloquy << "WARNING: Unable to fix header for " << bader_file << "." << endl;
        oss << soliloquy << "Please check." << endl;
        oss << endl;
      }
    } else {
      //echo warning, do not exit
      oss << endl;
      oss << soliloquy << "WARNING: " << bader_file << " was not created/is empty." << endl;
      oss << soliloquy << "This is VERY unusual and likely an issue with the Bader analysis." << endl;
      oss << endl;
    }
  }
  if(vpflow.flag("BADER::PRINT_INDEX")) {
    //AtIndex.dat
    if(vpflow.getattachedscheme("BADER::PRINT_INDEX") == "atom" || vpflow.getattachedscheme("BADER::PRINT_INDEX") == "both") {
      if(LDEBUG) oss << soliloquy << "CHECKING PRINT ATOM INDEX FILE" << endl;
      bader_file = "AtIndex.dat";
      if(aurostd::FileExist(bader_file) && aurostd::FileNotEmpty(bader_file)) {
        bzip_files.push_back(bader_file);
        move_files.push_back(bader_file);
        oss << soliloquy << "file created: " << bader_file << "." << endl;
        if(!prepare_CHGCAR_4_Jmol(bader_file, species_header, FALSE, oss)) {
          oss << endl;
          oss << soliloquy << "WARNING: Unable to fix header for " << bader_file << "." << endl;
          oss << soliloquy << "Please check." << endl;
          oss << endl;
        }
      } else {
        //echo warning, do not exit
        oss << endl;
        oss << soliloquy << "WARNING: " << bader_file << " was not created/is empty." << endl;
        oss << soliloquy << "This is VERY unusual and likely an issue with the Bader analysis." << endl;
        oss << endl;
      }
    }
    //BvIndex.dat
    if(vpflow.getattachedscheme("BADER::PRINT_INDEX") == "bader" || vpflow.getattachedscheme("BADER::PRINT_INDEX") == "both") {
      if(LDEBUG) oss << soliloquy << "CHECKING PRINT BADER INDEX FILE" << endl;
      bader_file = "BvIndex.dat";
      if(aurostd::FileExist(bader_file) && aurostd::FileNotEmpty(bader_file)) {
        bzip_files.push_back(bader_file);
        move_files.push_back(bader_file);
        oss << soliloquy << "file created: " << bader_file << "." << endl;
        if(!prepare_CHGCAR_4_Jmol(bader_file, species_header, FALSE, oss)) {
          oss << endl;
          oss << soliloquy << "WARNING: Unable to fix header for " << bader_file << "." << endl;
          oss << soliloquy << "Please check." << endl;
          oss << endl;
        }
      } else {
        //echo warning, do not exit
        oss << endl;
        oss << soliloquy << "WARNING: " << bader_file << " was not created/is empty." << endl;
        oss << soliloquy << "This is VERY unusual and likely an issue with the Bader analysis." << endl;
        oss << endl;
      }
    }
  }
  //move created files into correct directory, we can now work exclusively in this directory
  if(LDEBUG) oss << soliloquy << "MOVING FILES CREATED BY BADER TO DIRECTORY" << endl;
  for (uint i = 0; i < move_files.size(); i++) {
    oss << soliloquy << "moving " << move_files.at(i) << " to " << directory << "." << endl;
    if(!aurostd::MoveFile(move_files.at(i), directory)) {
      oss << endl;
      oss << soliloquy << "ERROR: Unable to move " << move_files.at(i) << "." << endl;
      oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
      oss << soliloquy << "Exiting." << endl;
      oss << endl;
      aurostd::RemoveFiles(remove_files);  //CLEAN UP BEFORE RETURNING FALSE
      oss << soliloquy << "changing directories back to original directory " << directory << "." << endl;
      chdir(work_dir.c_str());//directory.c_str());
      oss << soliloquy << "removing temporary directory " << BADER_RUN_DIRECTORY << "." << endl;
      if(!aurostd::RemoveDirectory(BADER_RUN_DIRECTORY)) {  //CLEAN UP BEFORE RETURNING FALSE
        oss << endl;
        oss << soliloquy << "WARNING: issues deleting directory " << BADER_RUN_DIRECTORY << "." << endl;
        oss << soliloquy << "Please check." << endl;
        oss << endl;
      }
      //BZIP aflow.CHGCAR_sum (only file created), and it's created in the right directory, so no need to move
      for (uint i = 0; i < bzip_files.size(); i++) {
        //if(aurostd::FileExist(bader_directory+"/"+bzip_files.at(i))) {
        bzip_file = directory + "/" + bzip_files.at(i);
        if(aurostd::FileExist(bzip_file)) {
          oss << soliloquy << "compressing " << bzip_files.at(i) << "." << endl;
          aurostd::BzipFile(bzip_file);
        }
      }
      return FALSE;
    }
  }
  oss << soliloquy << "changing directories back to original directory " << directory << "." << endl;
  chdir(work_dir.c_str());//directory.c_str());
  oss << soliloquy << "removing temporary directory " << BADER_RUN_DIRECTORY << "." << endl;
  if(!aurostd::RemoveDirectory(BADER_RUN_DIRECTORY)) {
    oss << endl;
    oss << soliloquy << "WARNING: issues deleting directory " << BADER_RUN_DIRECTORY << "." << endl;
    oss << soliloquy << "Please check." << endl;
    oss << endl;
  }

  if(LDEBUG) oss << soliloquy << "READING ACF.DAT FILE FOR BADER VOLUME AND CHARGE" << endl;
  // containers for AFLOW net charge data
  vector<string> lines, tokens;
  vector<double> bader_charge, volume, net_charge;

  // get bader charge and atomic volume results for each atom from ACF.dat
  //if(aurostd::FileExist(bader_directory+"/ACF.dat")) {
  if(aurostd::FileExist(directory + "/ACF.dat") && aurostd::FileNotEmpty(directory + "/ACF.dat")) {
    //aurostd::efile2vectorstring(bader_directory+"/ACF.dat",lines);
    aurostd::efile2vectorstring(directory + "/ACF.dat", lines);
    for (uint i = 0; i < lines.size(); i++) {
      aurostd::string2tokens(lines.at(i), tokens, " ");
      if(tokens.size() == 7 && tokens.at(4) != "CHARGE") {
        bader_charge.push_back(aurostd::string2utype<double>(tokens.at(4)));
        volume.push_back(aurostd::string2utype<double>(tokens.at(6)));
      }
    }
  } else {
    oss << endl;
    oss << soliloquy << "ERROR: Unable to locate ACF.dat (or its empty)." << endl;
    oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
    oss << soliloquy << "Exiting." << endl;
    oss << endl;
    aurostd::RemoveFiles(remove_files);
    //BZIP those that remain
    for (uint i = 0; i < bzip_files.size(); i++) {
      //if(aurostd::FileExist(bader_directory+"/"+bzip_files.at(i))) {
      bzip_file = directory + "/" + bzip_files.at(i);
      if(aurostd::FileExist(bzip_file)) {
        oss << soliloquy << "compressing " << bzip_files.at(i) << "." << endl;
        aurostd::BzipFile(bzip_file);
      }
    }
    return FALSE;
  }

  if(bader_charge.size() == 0 || volume.size() == 0 || bader_charge.size() != volume.size()) {
    oss << endl;
    oss << soliloquy << "ERROR: Unable to extract bader_charge/volume correctly." << endl;
    oss << soliloquy << "Perhaps the file is corrupt?" << endl;
    oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
    oss << soliloquy << "Exiting." << endl;
    oss << endl;
    aurostd::RemoveFiles(remove_files);
    //BZIP those that remain
    for (uint i = 0; i < bzip_files.size(); i++) {
      //if(aurostd::FileExist(bader_directory+"/"+bzip_files.at(i))) {
      bzip_file = directory + "/" + bzip_files.at(i);
      if(aurostd::FileExist(bzip_file)) {
        oss << soliloquy << "compressing " << bzip_files.at(i) << "." << endl;
        aurostd::BzipFile(bzip_file);
      }
    }
    return FALSE;
  }

  if(LDEBUG) oss << soliloquy << "TAKING DIFFERENCE BETWEEN VALENCE CHARGE AND BADER CHARGE" << endl;
  // calculate net charge on each atom
  net_charge = pflow::VVdiff(vZVAL, bader_charge);

  // make sure it's zero! otherwise spit warning
  double net_sum = 0;
  for (uint i = 0; i < net_charge.size(); i++) {
    net_sum += net_charge.at(i);
  }
  if(abs(net_sum) > 1e-3) {
    oss << endl;
    oss << soliloquy << "WARNING: The sum of the partial charges is non-zero (sum=" << net_sum << ",tol=1e-3)." << endl;
    oss << endl;
  }

  // write net charge on each atom to file abader.out
  string abader_out = prototype + "_abader.out";
  if(LDEBUG) oss << soliloquy << "WRITING RESULTS TO " << abader_out << endl;
  string net_charges_string = "bader_net_charges=";
  string volume_string = "bader_atomic_volumes=";
  string abader_spacer = "[AFLOW] **************************************************************************************************************************";
  //abader_out works for both command line users (non-library use) and library use
  stringstream abader_ss;

  oss << soliloquy << "writing net charges to " << abader_out << "." << endl;
  string column1 = "Atom #";
  string column2 = "Atom Type";
  string column3 = "Net Charge (electrons)";
  string column4 = "Volume (Angst^3)";
  string spaces = "    ";  //4 spaces
  if(vpflow.flag("BADER::AFLOWLIB_LIBRARY")) {
    abader_ss << abader_spacer << endl;
    abader_ss << "[BADER_RESULTS]START" << endl;
  } else {
    abader_ss << column1 << spaces << column2 << spaces << column3 << spaces << column4 << endl;  // header of labels
  }

  stringstream num_prec;
  atomCOUNT = 0;
  for (uint i = 0; i < vspecies.size(); i++) {
    for (uint speciesCOUNT = 0; speciesCOUNT < (uint)num_each_type.at(i); speciesCOUNT++) {
      if(vpflow.flag("BADER::AFLOWLIB_LIBRARY")) {
        num_prec.str("");
        num_prec << std::fixed << setprecision(4) << net_charge.at(atomCOUNT);
        net_charges_string.append(num_prec.str());
        num_prec.str("");
        num_prec << std::fixed << setprecision(4) << volume.at(atomCOUNT);
        volume_string.append(num_prec.str());
        if(!((i == vspecies.size() - 1) && (speciesCOUNT == (uint)num_each_type.at(i) - 1))) {
          net_charges_string.append(",");
          volume_string.append(",");
        }
      } else {
        abader_ss << aurostd::PaddedPRE(aurostd::utype2string(atomCOUNT + 1), column1.length(), " ") << spaces;  //right align
        abader_ss << aurostd::PaddedPOST(vspecies.at(i), column2.length(), " ") << spaces;                       //left align
        num_prec.str("");
        num_prec << std::fixed << setprecision(4) << net_charge.at(atomCOUNT);
        abader_ss << aurostd::PaddedPRE(num_prec.str(), column3.length(), " ") << spaces;  //right align
        num_prec.str("");
        num_prec << std::fixed << setprecision(4) << volume.at(atomCOUNT);
        abader_ss << aurostd::PaddedPRE(num_prec.str(), column4.length(), " ") << endl;  //right align
      }
      atomCOUNT++;
    }
  }

  //FOR AFLOWLIB LABELS
  //vpflow.flag("BADER::AFLOWLIB_LIBRARY",TRUE);
  if(vpflow.flag("BADER::AFLOWLIB_LIBRARY")) {
    //units
    net_charges_string.append(" (electrons)");
    volume_string.append(" (Angst^3)");
    abader_ss << net_charges_string << endl;
    abader_ss << volume_string << endl;
    abader_ss << "[BADER_RESULTS]STOP" << endl;
    abader_ss << abader_spacer << endl;
    aurostd::stringstream2file(abader_ss, directory + "/" + abader_out);
  } else {
    //aurostd::stringstream2file(abader_ss,bader_directory+"/"+abader_out);
    aurostd::stringstream2file(abader_ss, directory + "/" + abader_out);
    bzip_files.push_back(abader_out);  //DO NOT ZIP FOR LIBRARY
  }

  //consolidate, if specified by FLAG
  string species_file;
  if(vpflow.flag("BADER::CONSOLIDATE_ATOMS2SPECIES")) {
    //use sum_files
    vector<string> species_files, bader_atom_files, sum_files;
    atomCOUNT = 1;
    oss << soliloquy << "attempting to add individual atom contributions together for each species." << endl;
    //cycle through to CONSOLIDATE_ATOMS3SPECIES
    for (uint i = 0; i < vspecies.size(); i++) {
      species_file = "BvAt_" + vspecies.at(i) + ".dat";
      bzip_files.push_back(species_file);
      species_files.push_back(directory + "/" + species_file);
      //species_files.push_back(bader_directory+"/"+species_file);
      for (uint speciesCOUNT = 0; speciesCOUNT < (uint)num_each_type.at(i); speciesCOUNT++) {
        bader_file = "BvAt" + aurostd::PaddedNumString(atomCOUNT, 4) + ".dat";
        //bader_atom_files.push_back(bader_directory+"/"+bader_file);
        bader_atom_files.push_back(directory + "/" + bader_file);
        sum_files.push_back(bader_atom_files.back());
        if((!aurostd::FileExist(bader_atom_files.back())) || aurostd::FileEmpty(bader_atom_files.back())) {
          oss << endl;
          oss << soliloquy << "ERROR: " << bader_atom_files.back() << " not found/is empty." << endl;
          oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
          oss << soliloquy << "Exiting." << endl;
          oss << endl;
          aurostd::RemoveFiles(remove_files);
          //BZIP those that remain
          for (uint i = 0; i < bzip_files.size(); i++) {
            //if(aurostd::FileExist(bader_directory+"/"+bzip_files.at(i))) {
            bzip_file = directory + "/" + bzip_files.at(i);
            if(aurostd::FileExist(bzip_file)) {
              oss << soliloquy << "compressing " << bzip_files.at(i) << "." << endl;
              aurostd::BzipFile(bzip_file);
            }
          }
          return FALSE;
        }
        atomCOUNT++;
      }
      if(!pflow::CHGSUM(species_header, sum_files, species_files.back(), oss)) {
        oss << endl;
        oss << soliloquy << "ERROR: Unable to sum atom files to create " << species_files.back() << "." << endl;
        oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
        oss << soliloquy << "Exiting." << endl;
        oss << endl;
        aurostd::RemoveFiles(remove_files);
        //BZIP those that remain
        for (uint i = 0; i < bzip_files.size(); i++) {
          //if(aurostd::FileExist(bader_directory+"/"+bzip_files.at(i))) {
          bzip_file = directory + "/" + bzip_files.at(i);
          if(aurostd::FileExist(bzip_file)) {
            oss << soliloquy << "compressing " << bzip_files.at(i) << "." << endl;
            aurostd::BzipFile(bzip_file);
          }
        }
        return FALSE;
      }
      oss << soliloquy << "file created: " << species_files.back() << "." << endl;
      sum_files.clear();
    }
    //if REMOVE_BADER_ATOMS
    if(vpflow.flag("BADER::REMOVE_BADER_ATOMS")) {
      if(LDEBUG) oss << soliloquy << "REMOVE BADER ATOMS FLAG SPECIFIED" << endl;
      remove_files.insert(remove_files.end(), bader_atom_files.begin(), bader_atom_files.end());
    }
    if(vpflow.flag("BADER::JVXL_ALL_SPECIES")) {
      oss << soliloquy << "attempting to create .jvxl files for visualization in JMOL." << endl;
      string output_file, cutoff_string;
      bool cyclic = FALSE;
      if(vpflow.flag("BADER::JVXL_CYCLIC")) { cyclic = TRUE; }
      //if(downsample_ratios.size()==0) cyclic=FALSE; //OBVIOUS
      //check that they are same size if not cyclic
      if(!cyclic) {
        if(LDEBUG) oss << soliloquy << "PERFORMING SETS ROUTINE WITH CUTOFFS AND DOWNSAMPLE RATIOS" << endl;
        //cutoffs and downsample_ratios share same index
        if(downsample_ratios.size() != 0) {
          if(LDEBUG) oss << soliloquy << "NO DOWNSAMPLE RATIOS SPECIFIED" << endl;
          if(cutoffs.size() != downsample_ratios.size()) {
            oss << endl;
            oss << soliloquy << "ERROR: Number of cutoffs != number of downsample ratios." << endl;
            oss << soliloquy << "This is necessary for sets." << endl;
            oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
            oss << soliloquy << "Exiting." << endl;
            oss << endl;
            aurostd::RemoveFiles(remove_files);  //CLEAN UP BEFORE RETURNING FALSE
            //BZIP those that remain
            for (uint i = 0; i < bzip_files.size(); i++) {
              //if(aurostd::FileExist(bader_directory+"/"+bzip_files.at(i))) {
              bzip_file = directory + "/" + bzip_files.at(i);
              if(aurostd::FileExist(bzip_file)) {
                oss << soliloquy << "compressing " << bzip_files.at(i) << "." << endl;
                aurostd::BzipFile(bzip_file);
              }
            }
            return FALSE;
          }
        }
        for (uint i = 0; i < species_files.size(); i++) {
          for (uint j = 0; j < cutoffs.size(); j++) {
            //get cutoff string, 25 instead of 0.25
            num_prec.str("");
            num_prec << std::fixed << setprecision(2) << cutoffs.at(j);  //avoid issues with setprecision(20) in string2utype<double>
            cutoff_string = num_prec.str();
            //cutoff_string=aurostd::utype2string(cutoffs.at(j));
            if(aurostd::substring2bool(cutoff_string, "0.")) { cutoff_string = aurostd::StringSubst(cutoff_string, "0.", ""); }
            if(aurostd::substring2bool(cutoff_string, ".")) { cutoff_string = aurostd::StringSubst(cutoff_string, ".", ""); }
            //output_file=bader_directory+"/"+prototype+"_Bader_"+cutoff_string+"_";
            output_file = directory + "/" + prototype + "_Bader_" + cutoff_string + "_";
            if(downsample_ratios.size() != 0 && !(vpflow.flag("BADER::AFLOWLIB_LIBRARY"))) {
              output_file.append(aurostd::utype2string(downsample_ratios.at(j)) + "_");
            }
            output_file.append(vspecies.at(i) + ".jvxl");
            oss << soliloquy << "output file will be " << output_file << "." << endl;
            if(downsample_ratios.size() != 0) {
              if(!pflow::CHGCAR2JVXL(species_files.at(i), cutoffs.at(j), downsample_ratios.at(j), output_file, oss)) {
                oss << endl;
                oss << soliloquy << "ERROR: Unable to create " << output_file << "." << endl;
                oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
                oss << soliloquy << "Exiting." << endl;
                oss << endl;
                aurostd::RemoveFiles(remove_files);
                //BZIP those that remain
                for (uint i = 0; i < bzip_files.size(); i++) {
                  //if(aurostd::FileExist(bader_directory+"/"+bzip_files.at(i))) {
                  bzip_file = directory + "/" + bzip_files.at(i);
                  if(aurostd::FileExist(bzip_file)) {
                    oss << soliloquy << "compressing " << bzip_files.at(i) << "." << endl;
                    aurostd::BzipFile(bzip_file);
                  }
                }
                return FALSE;
              }
            } else {
              if(!pflow::CHGCAR2JVXL(species_files.at(i), cutoffs.at(j), output_file, oss)) {
                oss << endl;
                oss << soliloquy << "ERROR: Unable to create " << output_file << "." << endl;
                oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
                oss << soliloquy << "Exiting." << endl;
                oss << endl;
                aurostd::RemoveFiles(remove_files);
                //BZIP those that remain
                for (uint i = 0; i < bzip_files.size(); i++) {
                  //if(aurostd::FileExist(bader_directory+"/"+bzip_files.at(i))) {
                  bzip_file = directory + "/" + bzip_files.at(i);
                  if(aurostd::FileExist(bzip_file)) {
                    oss << soliloquy << "compressing " << bzip_files.at(i) << "." << endl;
                    aurostd::BzipFile(bzip_file);
                  }
                }
                return FALSE;
              }
            }
          }
        }
      } else {
        if(LDEBUG) oss << soliloquy << "PERFORMING CYCLIC ROUTINE WITH CUTOFFS AND DOWNSAMPLE RATIOS" << endl;
        //if cyclic specified, it must have downsample_ratios, what else is it going to cycle over?
        for (uint i = 0; i < species_files.size(); i++) {
          for (uint j = 0; j < cutoffs.size(); j++) {
            for (uint k = 0; k < downsample_ratios.size(); k++) {
              //get cutoff string, 25 instead of 0.25
              num_prec.str("");
              num_prec << std::fixed << setprecision(2) << cutoffs.at(j);  //avoid issues with setprecision(20) in string2utype<double>
              cutoff_string = num_prec.str();
              //cutoff_string=aurostd::utype2string(cutoffs.at(j));
              if(aurostd::substring2bool(cutoff_string, "0.")) { cutoff_string = aurostd::StringSubst(cutoff_string, "0.", ""); }
              if(aurostd::substring2bool(cutoff_string, ".")) { cutoff_string = aurostd::StringSubst(cutoff_string, ".", ""); }
              //output_file=bader_directory+"/"+prototype+"_Bader_"+cutoff_string+"_"
              output_file = directory + "/" + prototype + "_Bader_" + cutoff_string + "_";
              if(!(vpflow.flag("BADER::AFLOWLIB_LIBRARY"))) {
                output_file.append(aurostd::utype2string(downsample_ratios.at(k)) + "_");
              }
              output_file.append(vspecies.at(i) + ".jvxl");
              oss << soliloquy << "output file will be " << output_file << "." << endl;
              if(!pflow::CHGCAR2JVXL(species_files.at(i), cutoffs.at(j), downsample_ratios.at(k), output_file, oss)) {
                oss << endl;
                oss << soliloquy << "ERROR: Unable to create " << output_file << "." << endl;
                oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
                oss << soliloquy << "Exiting." << endl;
                oss << endl;
                aurostd::RemoveFiles(remove_files);
                //BZIP those that remain
                for (uint i = 0; i < bzip_files.size(); i++) {
                  //if(aurostd::FileExist(bader_directory+"/"+bzip_files.at(i))) {
                  bzip_file = directory + "/" + bzip_files.at(i);
                  if(aurostd::FileExist(bzip_file)) {
                    oss << soliloquy << "compressing " << bzip_files.at(i) << "." << endl;
                    aurostd::BzipFile(bzip_file);
                  }
                }
                return FALSE;
              }
            }
          }
        }
      }
      //if KEEP::JVXL_ONLY
      if(vpflow.flag("BADER::KEEP::JVXL_ONLY")) {
        if(LDEBUG) oss << soliloquy << "KEEP JVXL FILES ONLY FLAG SPECIFIED" << endl;
        remove_files.insert(remove_files.end(), species_files.begin(), species_files.end());
      }
    }
  }
  //REMOVE FILES, round 2
  for (uint i = 0; i < remove_files.size(); i++) {
    oss << soliloquy << "removing " << remove_files.at(i) << "." << endl;
    aurostd::RemoveFile(remove_files.at(i));
  }

  //BZIP those that remain
  for (uint i = 0; i < bzip_files.size(); i++) {
    //if(aurostd::FileExist(bader_directory+"/"+bzip_files.at(i))) {
    bzip_file = directory + "/" + bzip_files.at(i);
    if(aurostd::FileExist(bzip_file)) {
      oss << soliloquy << "compressing " << bzip_files.at(i) << "." << endl;
      aurostd::BzipFile(bzip_file);
    }
  }
  return TRUE;
}
}  // namespace bader_functions

//***************************************************************************//
// bader_functions::BaderCalc
//***************************************************************************//
namespace bader_functions {
bool BaderCalc(const string& bader_options,
               const string& prototype,
               const deque<string>& vspecies,
               const deque<int>& num_each_type,
               const vector<double>& vZVAL,
               const vector<double>& cutoffs,
               const vector<int>& downsample_ratios,
               string& directory,
               ostream& oss) {
  //NO FLAGS GIVEN
  aurostd::xoption bader_flags;
  return BaderCalc(bader_flags, bader_options, prototype, vspecies, num_each_type, vZVAL, cutoffs, downsample_ratios, directory, oss);
}
}  // namespace bader_functions

//***************************************************************************//
// bader_functions::FixDirectory
//***************************************************************************//
namespace bader_functions {
void FixDirectory(string& directory, ostream& oss) {
  string soliloquy = "bader_functions::FixDirectory():  ";  // so you know who's speaking
  string execution_path = aurostd::execute2string("pwd");
  // translate "./" into a full path
  bool LDEBUG = (FALSE || XHOST.DEBUG);
  if(LDEBUG) oss << soliloquy << "CHECK DIRECTORY" << endl;
  if(directory.at(directory.length() - 1) == '/') { directory = directory.substr(0, directory.size() - 1); }
  if(directory == "./" || directory == ".") {
    directory = execution_path;
  } else if(directory.at(0) != '/') {
    directory = execution_path + "/" + directory;
  }
}
}

//***************************************************************************//
// bader_functions::Flags2BaderCommands
//***************************************************************************//
namespace bader_functions {
bool Flags2BaderCommands(aurostd::xoption& vpflow, string& bader_options, ostream& oss) {  //corey
  // Perform Bader analysis by means of the code from Henkelman Group at UT, Austin
  // Results for the Bader volume and net charge in the volume are stored in
  // the vectors passed into the function, volume and charge.

  string soliloquy = "bader_functions::BaderCalcFLAGS2COMMANDS():  ";  // so you know who's speaking
  // DEBUG
  bool LDEBUG = (FALSE || XHOST.DEBUG);
  if(LDEBUG) oss << soliloquy << "BEGIN" << endl;

  //CO 180220 - directory stuff for logging
  string directory = ".";  //default
  if(XHOST.vflag_control.flag("DIRECTORY")) {
    directory = XHOST.vflag_control.getattachedscheme("DIRECTORY");  //XHOST.vflag_control.flag("DIRECTORY");//aflags.Directory;//vpflow.getattachedscheme("BADER::DIRECTORY");
  }
  FixDirectory(directory, oss);
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //
  //  Start FLAG automation
  //
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //this is for ease of specifying
  if(vpflow.flag("BADER::AFLOWLIB_LIBRARY")) {
    //vpflow.flag("BADER::JVXL_ALL_SPECIES") has to be turned on previously, need input CUTOFF and DOWNSAMPLES
    //vpflow.flag("BADER::JVXL_CYCLIC")      has to be turned on previously, need input CUTOFF and DOWNSAMPLES
    vpflow.flag("BADER::KEEP::JVXL_ONLY", TRUE);  //only matters if JVXL_ALL_SPECIES turned on already
    //vpflow.flag("BADER::CRITICAL_POINTS",TRUE); //not ready yet
    vpflow.flag("BADER::QUIET", TRUE);
  }
  //if we JVXL_ALL_SPECIES, it should automatically CONSOLIDATE_ATOMS2SPECIES
  //if we KEEP::JVXL_ONLY, it should REMOVE_BADER_ATOMS
  if(vpflow.flag("BADER::JVXL_ALL_SPECIES")) {
    if(!vpflow.flag("BADER::CONSOLIDATE_ATOMS2SPECIES")) { vpflow.flag("BADER::CONSOLIDATE_ATOMS2SPECIES", TRUE); }  //set consolidate_atoms2species if not already set
    if(vpflow.flag("BADER::KEEP::JVXL_ONLY")) {
      if(!vpflow.flag("BADER::REMOVE_BADER_ATOMS")) { vpflow.flag("BADER::REMOVE_BADER_ATOMS", TRUE); }  //set consolidate_atoms4species if not already set
    }
  }
  //same here
  //we can JUST specify CONSOLIDATE_ATOMS2SPECIES,and it automatically knows to print_all=atoms
  if(vpflow.flag("BADER::CONSOLIDATE_ATOMS2SPECIES")) {
    if(!vpflow.flag("BADER::PRINT_ALL")) {
      vpflow.flag("BADER::PRINT_ALL", TRUE);  //set print_all if not already set
      vpflow.addattachedscheme("BADER::PRINT_ALL", "atom", TRUE);
    } else {
      if(vpflow.getattachedscheme("BADER::PRINT_ALL") == "bader") {
        vpflow.purgeattachedscheme("BADER::PRINT_ALL");
        vpflow.addattachedscheme("BADER::PRINT_ALL", "both", TRUE);
      }
    }
  }
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //
  //  End FLAG automation
  //
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // Bader Options -- build a string of optional flags accepted by the bader code
  // discard any option not understood
  //
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  string misc_option;

  // calculate/nocalculate option
  if(LDEBUG) oss << soliloquy << "CHECK CALCULATE | NOCALCULATE" << endl;
  if(vpflow.flag("BADER::CALCULATE")) {
    if(LDEBUG) oss << soliloquy << "vpflow.getattachedscheme(\"BADER::CALCULATE\")=" << vpflow.flag("BADER::CALCULATE") << endl;
    misc_option = vpflow.getattachedscheme("BADER::CALCULATE");
    if(misc_option == "bader") {
      oss << soliloquy << "calculating Bader atoms in molecules." << endl;
      bader_options.append("-c bader ");
    } else if(misc_option == "voronoi") {
      oss << soliloquy << "calculating population analysis based on distance." << endl;
      bader_options.append("-c voronoi ");
    } else {
      oss << endl;
      oss << soliloquy << "ERROR: Unknown calculate option." << endl;
      oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
      oss << soliloquy << "Exiting." << endl;
      oss << endl;
      return FALSE;
    }
  } else if(vpflow.flag("BADER::NOCALCULATE")) {
    if(LDEBUG) oss << soliloquy << "vpflow.getattachedscheme(\"BADER::NOCALCULATE\")=" << vpflow.flag("BADER::NOCALCULATE") << endl;
    misc_option = vpflow.getattachedscheme("BADER::NOCALCULATE");
    if(misc_option == "bader") {
      oss << soliloquy << "NOT calculating Bader atoms in molecules." << endl;
      bader_options.append("-n bader ");
    } else if(misc_option == "voronoi") {
      oss << soliloquy << "NOT calculating population analysis based on distance." << endl;
      bader_options.append("-n voronoi ");
    } else {
      oss << endl;
      oss << soliloquy << "ERROR: Unknown nocalculate option." << endl;
      oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
      oss << soliloquy << "Exiting." << endl;
      oss << endl;
      return FALSE;
    }
  } else {
    oss << soliloquy << "calculating default Bader atoms in molecules." << endl;
  }

  // partition option
  if(LDEBUG) oss << soliloquy << "CHECK PARTITION" << endl;
  if(vpflow.flag("BADER::PARTITION")) {
    if(LDEBUG) oss << soliloquy << "vpflow.getattachedscheme(\"BADER::NOCALCULATE\")=" << vpflow.flag("BADER::PARTITION") << endl;
    misc_option = vpflow.getattachedscheme("BADER::PARTITION");
    if(misc_option == "neargrid") {
      oss << soliloquy << "using the default near-grid bader partitioning algorithm." << endl;
      bader_options.append("-b neargrid ");
    } else if(misc_option == "ongrid") {
      oss << soliloquy << "using the on-grid bader partitioning algorithm." << endl;
      bader_options.append("-b ongrid ");
    } else {
      oss << endl;
      oss << soliloquy << "ERROR: Unknown partition option." << endl;
      oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
      oss << soliloquy << "Exiting." << endl;
      oss << endl;
      return FALSE;
    }
  } else {
    oss << soliloquy << "using the default near-grid bader partitioning algorithm." << endl;
  }

  // refine_edge_method option
  if(LDEBUG) oss << soliloquy << "CHECK REFINE_EDGE_METHOD" << endl;
  if(vpflow.flag("BADER::REFINE_EDGE_METHOD")) {
    if(LDEBUG) oss << soliloquy << "vpflow.getattachedscheme(\"BADER::REFINE_EDGE_METHOD\")=" << vpflow.flag("BADER::REFINE_EDGE_METHOD") << endl;
    misc_option = vpflow.getattachedscheme("BADER::REFINE_EDGE_METHOD");
    if(misc_option == "-1" || misc_option == "-2" || misc_option == "-3") {
      oss << soliloquy << "requesting refinement method " << misc_option << "." << endl;
      bader_options.append("-r " + misc_option + " ");
    } else {
      oss << endl;
      oss << soliloquy << "ERROR: Unknown refine_edge_method option " + misc_option + "." << endl;
      oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
      oss << soliloquy << "Exiting." << endl;
      oss << endl;
      return FALSE;
    }
  } else {
    oss << soliloquy << "using default refinement method -1." << endl;
  }

  // reference option handled internally in BaderCalc function

  // vacuum option
  if(LDEBUG) oss << soliloquy << "CHECK VACUUM" << endl;
  if(vpflow.flag("BADER::VACUUM")) {
    if(LDEBUG) oss << soliloquy << "vpflow.getattachedscheme(\"BADER::VACUUM\")=" << vpflow.flag("BADER::VACUUM") << endl;
    misc_option = vpflow.getattachedscheme("BADER::VACUUM");
    if(misc_option == "off") {
      oss << soliloquy << "NOT assigning low density points to vacuum." << endl;
      bader_options.append("-vac off ");
    } else if(misc_option == "auto") {
      oss << soliloquy << "automatically assigning density points below 1E-3 e/Ang^3 to vacuum." << endl;
      bader_options.append("-vac auto ");
    } else if(misc_option == "") {
      oss << endl;
      oss << soliloquy << "ERROR: No density threshold given." << endl;
      oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
      oss << soliloquy << "Exiting." << endl;
      oss << endl;
      return FALSE;
    } else {
      double vac_test = aurostd::string2utype<double>(misc_option);
      if(vac_test == 0) {
        oss << endl;
        oss << soliloquy << "ERROR: Incorrect input for vacuum density. Must be a number (double)." << endl;
        oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
        oss << soliloquy << "Exiting." << endl;
        oss << endl;
        return FALSE;
      }
      oss << soliloquy << "assigning density points lower than " << misc_option << " to vacuum." << endl;
      bader_options.append("-vac " + misc_option + " ");
    }
  } else {
    oss << soliloquy << "automatically assigning density points below 1E-3 e/Ang^3 to vacuum." << endl;
    bader_options.append("-vac auto ");  //not actually default in bader code
  }

  // terminate option
  if(LDEBUG) oss << soliloquy << "CHECK TERMINATE" << endl;
  if(vpflow.flag("BADER::TERMINATE")) {
    if(LDEBUG) oss << soliloquy << "vpflow.getattachedscheme(\"BADER::TERMINATE\")=" << vpflow.flag("BADER::TERMINATE") << endl;
    misc_option = vpflow.getattachedscheme("BADER::TERMINATE");
    if(misc_option == "known") {
      oss << soliloquy << "terminating trajectories when a point is surrounded by known points." << endl;
      bader_options.append("-m known ");
    } else if(misc_option == "max") {
      oss << soliloquy << "terminating trajectories when charge density maximum is reached." << endl;
      bader_options.append("-m max ");
    } else {
      oss << endl;
      oss << soliloquy << "ERROR: Unknown termination option " + misc_option + "." << endl;
      oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
      oss << soliloquy << "Exiting." << endl;
      oss << endl;
      return FALSE;
    }
  }

  // critical points option
  if(LDEBUG) oss << soliloquy << "CHECK CRITICAL POINTS" << endl;
  if(vpflow.flag("BADER::CRITICAL_POINTS")) {
    if(LDEBUG) oss << soliloquy << "vpflow.flag(\"BADER::CRITICAL_POINTS\")" << endl;
    oss << soliloquy << "finding the critical points." << endl;  // not shown anyway (only in library)
    //bader_options.append("-cp ");   //need to run as a separate command
  }

  // print_all option
  if(LDEBUG) oss << soliloquy << "CHECK PRINT_ALL" << endl;
  if(vpflow.flag("BADER::PRINT_ALL")) {
    if(LDEBUG) oss << soliloquy << "vpflow.getattachedscheme(\"BADER::PRINT_ALL\")=" << vpflow.flag("BADER::PRINT_ALL") << endl;
    misc_option = vpflow.getattachedscheme("BADER::PRINT_ALL");
    //BvAtxxxx.dat
    if(misc_option == "atom") {
      oss << soliloquy << "combining all volumes associated with an atom and writing them to BvAtxxxx.dat." << endl;
      bader_options.append("-p all_atom ");
      //Bvolxxxx.dat
    } else if(misc_option == "bader") {
      oss << soliloquy << "writing each Bader volume to Bvolxxxx.dat." << endl;
      bader_options.append("-p all_bader ");
    } else if(misc_option == "both") {
      oss << soliloquy << "combining all volumes associated with an atom and writing them to BvAtxxxx.dat AND write each Bader volume to Bvolxxxx.dat." << endl;
      bader_options.append("-p all_atom -p all_bader ");
    } else {
      oss << endl;
      oss << soliloquy << "ERROR: Unknown print_all option " + misc_option + "." << endl;
      oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
      oss << soliloquy << "Exiting." << endl;
      oss << endl;
      return FALSE;
    }
  }

  // print_index option
  if(LDEBUG) oss << soliloquy << "CHECK PRINT_INDEX" << endl;
  if(vpflow.flag("BADER::PRINT_INDEX")) {
    if(LDEBUG) oss << soliloquy << "vpflow.getattachedscheme(\"BADER::PRINT_INDEX\")=" << vpflow.flag("BADER::PRINT_INDEX") << endl;
    misc_option = vpflow.getattachedscheme("BADER::PRINT_INDEX");
    //AtIndex.dat
    if(misc_option == "atom") {
      oss << soliloquy << "printing atomic volume indices and writing them to AtIndex.dat." << endl;
      bader_options.append("-p atom_index ");
      //BvIndex.dat
    } else if(misc_option == "bader") {
      oss << soliloquy << "printing bader volume indices and writing them to BvIndex.dat." << endl;
      bader_options.append("-p bader_index ");
    } else if(misc_option == "both") {
      oss << soliloquy << "printing both atomic and bader volume indices and writing them to AtIndex.dat and BvIndex.dat, respectively." << endl;
      bader_options.append("-p atom_index -p bader_index ");
    } else {
      oss << endl;
      oss << soliloquy << "ERROR: Unknown print_index option " + misc_option + "." << endl;
      oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
      oss << soliloquy << "Exiting." << endl;
      oss << endl;
      return FALSE;
    }
  }

  // print_select_atom/print_select_bader option
  vector<string> tokens;
  string push_command;
  if(LDEBUG) oss << soliloquy << "CHECK PRINT_SELECT_ATOM | PRINT_SELECT_BADER" << endl;
  //BvAtxxxx.dat
  if(vpflow.flag("BADER::PRINT_SELECT_ATOM")) {
    if(LDEBUG) oss << soliloquy << "vpflow.getattachedscheme(\"BADER::PRINT_SELECT_ATOM\")=" << vpflow.flag("BADER::PRINT_SELECT_ATOM") << endl;
    if(vpflow.flag("BADER::PRINT_ALL") && vpflow.getattachedscheme("BADER::PRINT_ALL") == "atom") {
      oss << soliloquy << "ignoring print select atom command--already specified print ALL atoms command." << endl;
    } else {
      if(!getPushCommand(vpflow.getattachedscheme("BADER::PRINT_SELECT_ATOM"), push_command, oss)) {
        oss << endl;
        oss << soliloquy << "ERROR: Unable to gather input for print select atom command." << endl;
        oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
        oss << soliloquy << "Exiting." << endl;
        oss << endl;
        return FALSE;
      }
      oss << soliloquy << "combining all volumes associated with atom(s) " + push_command + "and writing them to BvAtxxxx.dat." << endl;  //push command comes with POST-space
      bader_options.append("-p sel_atom " + push_command);
    }
  }
  //B_wexxxx.dat
  if(vpflow.flag("BADER::PRINT_SELECT_BADER")) {
    if(LDEBUG) oss << soliloquy << "vpflow.getattachedscheme(\"BADER::PRINT_SELECT_BADER\")=" << vpflow.flag("BADER::PRINT_SELECT_BADER") << endl;
    if(vpflow.flag("BADER::PRINT_ALL") && vpflow.getattachedscheme("BADER::PRINT_ALL") == "bader") {
      oss << soliloquy << "ignoring print select bader volume command--already specified print ALL bader volumes command." << endl;
    } else {
      if(!getPushCommand(vpflow.getattachedscheme("BADER::PRINT_SELECT_BADER"), push_command, oss)) {
        oss << endl;
        oss << soliloquy << "ERROR: Unable to gather input for print select bader command." << endl;
        oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
        oss << soliloquy << "Exiting." << endl;
        oss << endl;
        return FALSE;
      }
      oss << soliloquy << "writing bader volume(s) " + push_command + "to Bvolxxxx.dat." << endl;  //push command comes with POST-space
      bader_options.append("-p sel_bader " + push_command);
    }
  }

  // print_sum_atom/print_sum_bader option
  //BvAt_summed.dat
  if(LDEBUG) oss << soliloquy << "CHECK PRINT_SUM_ATOM | PRINT_SUM_BADER" << endl;
  if(vpflow.flag("BADER::PRINT_SUM_ATOM")) {
    if(LDEBUG) oss << soliloquy << "vpflow.getattachedscheme(\"BADER::PRINT_SUM_ATOM\")=" << vpflow.flag("BADER::PRINT_SUM_ATOM") << endl;
    if(!getPushCommand(vpflow.getattachedscheme("BADER::PRINT_SUM_ATOM"), push_command, oss)) {
      oss << endl;
      oss << soliloquy << "ERROR: Unable to gather input for print sum atom command." << endl;
      oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
      oss << soliloquy << "Exiting." << endl;
      oss << endl;
      return FALSE;
    }
    oss << soliloquy << "combining all volumes associated with atoms " + push_command + "and writing them to BvAt_summed.dat." << endl;  //push command comes with POST-space
    bader_options.append("-p sum_atom " + push_command);
    //Bvol_summed.dat
  }
  if(vpflow.flag("BADER::PRINT_SUM_BADER")) {
    if(LDEBUG) oss << soliloquy << "vpflow.getattachedscheme(\"BADER::PRINT_SUM_BADER\")=" << vpflow.flag("BADER::PRINT_SUM_BADER") << endl;
    if(!getPushCommand(vpflow.getattachedscheme("BADER::PRINT_SUM_BADER"), push_command, oss)) {
      oss << endl;
      oss << soliloquy << "ERROR: Unable to gather input for print sum bader command." << endl;
      oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
      oss << soliloquy << "Exiting." << endl;
      oss << endl;
      return FALSE;
    }
    oss << soliloquy << "combining bader volumes " + push_command + "and writing them to Bvol_summed.dat." << endl;  //push command comes with POST-space
    bader_options.append("-p sum_bader " + push_command);
  }

  // quiet/verbose option, different from AFLOW's silent (which would silent EVERYTHING, including other functions in AFLOW)
  if(LDEBUG) oss << soliloquy << "CHECK QUIET" << endl;
  if(vpflow.flag("BADER::QUIET")) {
    if(LDEBUG) oss << soliloquy << "vpflow.flag(\"BADER::QUIET\")" << endl;
    oss << soliloquy << "silencing output." << endl;  // not shown anyway (only in library)
  } else {
    if(LDEBUG) oss << soliloquy << "vpflow.flag(\"BADER::VERBOSE\")" << endl;
    oss << soliloquy << "displaying output." << endl;
    bader_options.append("-v ");
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // END Bader Options
  //
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  return TRUE;
}
}  // namespace bader_functions

//***************************************************************************//
// bader_functions::getPushCommand
//***************************************************************************//
namespace bader_functions {
bool getPushCommand(const string& misc_option, string& push_command, ostream& oss) {
  //this converts aflow lists/ranges to bader lists/ranges
  string soliloquy = "bader_functions::getPushCommand():  ";  // so you know who's speaking
  //debug
  bool LDEBUG = (FALSE || XHOST.DEBUG);
  if(LDEBUG) oss << soliloquy << "BEGIN" << endl;
  
  //CO 180220 - directory stuff for logging
  string directory = ".";  //default
  if(XHOST.vflag_control.flag("DIRECTORY")) {
    directory = XHOST.vflag_control.getattachedscheme("DIRECTORY");  //XHOST.vflag_control.flag("DIRECTORY");//aflags.Directory;//vpflow.getattachedscheme("BADER::DIRECTORY");
  }
  FixDirectory(directory, oss);

  vector<int> tokens;
  int test_value;
  push_command = "";
  //check that there is no mixing between list and range
  if(LDEBUG) oss << soliloquy << "CHECK THAT THERE'S NO MIXING BETWEEN LIST AND RANGE" << endl;
  if(aurostd::substring2bool(misc_option, ",") && aurostd::substring2bool(misc_option, "-")) {
    oss << endl;
    oss << soliloquy << "ERROR: Do not mix notation for LIST and RANGE." << endl;
    oss << soliloquy << "You may only have '::' OR ':'." << endl;
    oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
    oss << soliloquy << "Exiting." << endl;
    oss << endl;
    return FALSE;
  }
  //list
  if(LDEBUG) oss << soliloquy << "CHECK FOR LIST" << endl;
  if(aurostd::substring2bool(misc_option, ",")) {
    if(LDEBUG) oss << soliloquy << "LIST FOUND" << endl;
    aurostd::string2tokens<int>(misc_option, tokens, ",");
    for (uint i = 0; i < tokens.size(); i++) {
      if(tokens.at(i) == 0) {
        oss << endl;
        oss << soliloquy << "ERROR: Input for list must be numerical and above 0: " << tokens.at(i) << "." << endl;
        oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
        oss << soliloquy << "Exiting." << endl;
        oss << endl;
        return FALSE;
      }
      push_command.append(aurostd::utype2string(tokens.at(i)) + " ");
    }
    return TRUE;
  }
  //range or single item, no need to change syntax
  //do some checks
  if(LDEBUG) oss << soliloquy << "CHECK FOR RANGE" << endl;
  if(aurostd::substring2bool(misc_option, "-")) {
    if(LDEBUG) oss << soliloquy << "RANGE FOUND" << endl;
    //check that range only has two items
    aurostd::string2tokens<int>(misc_option, tokens, "-");
    if(tokens.size() != 2) {
      oss << endl;
      oss << soliloquy << "ERROR: Cannot specify more than one range." << endl;
      oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
      oss << soliloquy << "Exiting." << endl;
      oss << endl;
      return FALSE;
    }
    for (uint i = 0; i < tokens.size(); i++) {
      if(tokens.at(i) == 0) {
        oss << endl;
        oss << soliloquy << "ERROR: Input for range must be numerical and above 0: " << tokens.at(i) << "." << endl;
        oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
        oss << soliloquy << "Exiting." << endl;
        oss << endl;
        return FALSE;
      }
    }
    //check that first entry !> second entry
    if(tokens.at(0) > tokens.at(1)) {
      oss << endl;
      oss << soliloquy << "ERROR: First range element > second range element." << endl;
      oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
      oss << soliloquy << "Exiting." << endl;
      oss << endl;
      return FALSE;
    }
    push_command.append(aurostd::utype2string(tokens.at(0)) + "-" + aurostd::utype2string(tokens.at(1)) + " ");
    return TRUE;
  }
  //must be singular entry
  test_value = aurostd::string2utype<int>(misc_option);
  if(test_value == 0) {
    oss << endl;
    oss << soliloquy << "ERROR: Input for range must be numerical and above 0: " << test_value << "." << endl;
    oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
    oss << soliloquy << "Exiting." << endl;
    oss << endl;
    return FALSE;
  }
  push_command.append(aurostd::utype2string(test_value) + " ");
  return TRUE;
}
}  // namespace bader_functions

//***************************************************************************//
// bader_functions::listORrange2vec
//***************************************************************************//
namespace bader_functions {
bool listORrange2vec(const string& misc_option, vector<int>& vout, ostream& oss) {
  string soliloquy = "bader_functions::BaderCalc():  ";  // so you know who's speaking
  int test_value;
  //debug
  bool LDEBUG = (FALSE || XHOST.DEBUG);
  if(LDEBUG) oss << soliloquy << "BEGIN" << endl;
  
  //CO 180220 - directory stuff for logging
  string directory = ".";  //default
  if(XHOST.vflag_control.flag("DIRECTORY")) {
    directory = XHOST.vflag_control.getattachedscheme("DIRECTORY");  //XHOST.vflag_control.flag("DIRECTORY");//aflags.Directory;//vpflow.getattachedscheme("BADER::DIRECTORY");
  }
  FixDirectory(directory, oss);

  vout.clear();
  //check that there is no mixing between list and range
  if(LDEBUG) oss << soliloquy << "CHECK THAT THERE'S NO MIXING BETWEEN LIST AND RANGE" << endl;
  if(aurostd::substring2bool(misc_option, ",") && aurostd::substring2bool(misc_option, "-")) {
    oss << endl;
    oss << soliloquy << "ERROR: Do not mix notation for LIST and RANGE." << endl;
    oss << soliloquy << "You may only have '::' OR ':'." << endl;
    oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
    oss << soliloquy << "Exiting." << endl;
    oss << endl;
    return FALSE;
  }
  //list
  if(LDEBUG) oss << soliloquy << "CHECK FOR LIST" << endl;
  if(aurostd::substring2bool(misc_option, ",")) {
    if(LDEBUG) oss << soliloquy << "LIST FOUND" << endl;
    aurostd::string2tokens<int>(misc_option, vout, ",");
    for (uint i = 0; i < vout.size(); i++) {
      if(vout.at(i) == 0) {
        oss << endl;
        oss << soliloquy << "ERROR: Input for list must be numerical and above 0: " << vout.at(i) << "." << endl;
        oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
        oss << soliloquy << "Exiting." << endl;
        oss << endl;
        return FALSE;
      }
    }
    return TRUE;
  }
  //range
  if(LDEBUG) oss << soliloquy << "CHECK FOR RANGE" << endl;
  if(aurostd::substring2bool(misc_option, "-")) {
    if(LDEBUG) oss << soliloquy << "RANGE FOUND" << endl;
    //check that range only has two items
    vector<uint> tokens;
    aurostd::string2tokens<uint>(misc_option, tokens, "-");
    if(tokens.size() != 2) {
      oss << endl;
      oss << soliloquy << "ERROR: Cannot specify more than one range." << endl;
      oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
      oss << soliloquy << "Exiting." << endl;
      oss << endl;
      return FALSE;
    }
    for (uint i = 0; i < tokens.size(); i++) {
      if(tokens.at(i) == 0) {
        oss << endl;
        oss << soliloquy << "ERROR: Input for range must be numerical and above 0: " << tokens.at(i) << "." << endl;
        oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
        oss << soliloquy << "Exiting." << endl;
        oss << endl;
        return FALSE;
      }
    }
    //check that first entry !> second entry
    if(tokens.at(0) > tokens.at(1)) {
      oss << endl;
      oss << soliloquy << "ERROR: First range element > second range element." << endl;
      oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
      oss << soliloquy << "Exiting." << endl;
      oss << endl;
      return FALSE;
    }
    for (uint i = tokens.at(0); i < tokens.at(1) + 1; i++) {
      vout.push_back(i);
    }
    return TRUE;
  }
  //or single item
  test_value = aurostd::string2utype<int>(misc_option);
  if(test_value == 0) {
    oss << endl;
    oss << soliloquy << "ERROR: Input must be numerical and above 0: " << misc_option << "." << endl;
    oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
    oss << soliloquy << "Exiting." << endl;
    oss << endl;
    return FALSE;
  }
  vout.push_back(test_value);
  return TRUE;
}
}  // namespace bader_functions

//***************************************************************************//
// bader_functions::BaderExtensionFound
//***************************************************************************//
namespace bader_functions {
bool BaderExtensionFound(string FileNameIN, string& FileNameOUT, const string& directory) {
  //input filename and directory, output full path of compressed/unccompressed file
  //we only care about .static runs (aflow), so we will only look for files with .static extension
  //or no extension
  if(aurostd::EFileExist(directory + "/" + FileNameIN + ".static", FileNameOUT) && aurostd::FileNotEmpty(FileNameOUT)) { return TRUE; }
  if(aurostd::EFileExist(directory + "/" + FileNameIN, FileNameOUT) && aurostd::FileNotEmpty(FileNameOUT)) { return TRUE; }
  return FALSE;
}
}  // namespace bader_functions

//***************************************************************************//
// bader_functions::BaderExtensionFound
//***************************************************************************//
namespace bader_functions {
bool BaderExtensionFound(string FileNameIN, const string& directory) {
  string FileNameOUT;
  return BaderExtensionFound(FileNameIN, FileNameOUT, directory);
}
}  // namespace bader_functions

//***************************************************************************//
// bader_functions::adjust_header
//***************************************************************************//
namespace bader_functions {
void adjust_header(string& new_header, stringstream& FileIN_ss, ostream& oss) {
  string soliloquy = "bader_functions::adjust_header():  ";  // so you know who's speaking
  //debug
  bool LDEBUG = (FALSE || XHOST.DEBUG);
  if(LDEBUG) oss << soliloquy << "BEGIN" << endl;

  //adds new_header at top of stringstream
  string content;
  vector<string> vcontent;
  content = FileIN_ss.str();
  //cerr << "LOOOOK " << endl;
  //cerr << content << endl;
  //cerr << new_header << endl;
  aurostd::string2vectorstring(content, vcontent);
  //get right number of spaces for formatting
  uint white_space_counter1 = 0, white_space_counter2 = 0;
  while (isspace(vcontent.at(0).at(white_space_counter1))) { white_space_counter1++; }
  while (isspace(new_header.at(white_space_counter2))) { white_space_counter2++; }
  if(white_space_counter1 != white_space_counter2) {
    if(LDEBUG) oss << soliloquy << "ADDING " << white_space_counter1 << " TO HEADER" << endl;
    vcontent.at(0) = new_header.insert(0, white_space_counter1, ' ');
  } else {
    vcontent.at(0) = new_header;
  }
  //rebuild
  FileIN_ss.str("");
  if(LDEBUG) oss << soliloquy << "WRITING OUT STRINGSTREAM WITH ADJUSTED HEADER" << endl;
  for (uint i = 0; i < vcontent.size(); i++) { FileIN_ss << vcontent.at(i) << endl; }
}
}  // namespace bader_functions

//***************************************************************************//
// bader_functions::prepare_CHGCAR_4_Jmol
//***************************************************************************//
namespace bader_functions {
string prepare_CHGCAR_4_Jmol(aurostd::xoption vpflow) {
  ostringstream oss;
  //adjust header of CHGCAR for reading in Jmol
  string soliloquy = "bader_functions::prepare_CHGCAR_4_Jmol():  ";  // so you know who's speaking
  bool LDEBUG = (FALSE || XHOST.DEBUG);
  if(LDEBUG) oss << soliloquy << "BEGIN" << endl;
  string usage_usage = "aflow --prep4jmol [options]";
  string usage_options = aurostd::liststring2string("options = --usage",
                                                    "          --outcar=OUTCAR1,OUTCAR2,...",
                                                    "          --zip");
  // output usage
  if(LDEBUG) oss << soliloquy << "CHECK USAGE" << endl;
  if(vpflow.flag("PREPARE_CHGCAR_4_JMOL::USAGE")) {
    init::ErrorOption(cout, vpflow.getattachedscheme("PREPARE_CHGCAR_4_JMOL"), "bader_functions::prepare_CHGCAR_4_Jmol()", aurostd::liststring2string(usage_usage, usage_options));
    return oss.str();
  }
  if(LDEBUG) oss << soliloquy << "GATHER CHGCAR_FILES" << endl;
  vector<string> chgcar_files;
  string directory = aurostd::execute2string("pwd");
  string input = vpflow.getattachedscheme("PREPARE_CHGCAR_4_JMOL");
  if(aurostd::substring2bool(input, ",")) {
    aurostd::string2tokens(input, chgcar_files, ",");
  } else {
    chgcar_files.push_back(input);
  }

  string outcar_file;
  if(vpflow.flag("PREPARE_CHGCAR_4_JMOL::OUTCAR")) {
    outcar_file = vpflow.getattachedscheme("PREPARE_CHGCAR_4_JMOL::OUTCAR");
    if(outcar_file.at(0) != '/') {  //not root dir, so we can append directory
      outcar_file = directory + "/" + outcar_file;
    }
    oss << soliloquy << "OUTCAR file specied: " << outcar_file << "." << endl;
  }

  bool zip_file = false;
  if(vpflow.flag("PREPARE_CHGCAR_4_JMOL::ZIP")) { zip_file = true; }

  string species_header, path, file, chgcar_file;
  vector<string> _path;
  //TO RUN QUICKLY, SPECIFY OUTCAR
  for (uint i = 0; i < chgcar_files.size(); i++) {
    path = "";
    oss << soliloquy << "working on " << chgcar_files.at(i) << "." << endl;
    //remove file from path so we can search path
    chgcar_file = chgcar_files.at(i);
    if(chgcar_file.at(0) != '/') {  //not root dir, so we can append directory
      chgcar_file = directory + "/" + chgcar_file;
    }
    if(LDEBUG) oss << soliloquy << "CHGCAR " << chgcar_file << endl;
    //only look every time if outcar_file not specified, otherwise, only look once
    if(outcar_file.empty() || species_header.empty()) {
      if(LDEBUG) oss << soliloquy << "LOOKING FOR OUTCAR" << endl;
      if(LDEBUG) oss << soliloquy << "GETTING PATH FROM CHGCAR" << endl;
      if(aurostd::substring2bool(chgcar_file, "/")) {
        aurostd::string2tokens(chgcar_file, _path, "/");
        if(chgcar_file.at(0) == '/') {  //if root, which it will be
          path.append("/");
        }
        for (uint i = 0; i < _path.size() - 1; i++) {
          path.append(_path.at(i) + "/");  //last character will be '/', chgcar_file==_path.back();
        }
        path = path.substr(0, path.size() - 1);  //pop off last '/'
        file = _path.back();
      } else {
        path = ".";
      }  //chgcar_file==chgcar_file;}
      if(LDEBUG) oss << soliloquy << "PATH " << path << endl;
      if(!get_species_string(outcar_file, species_header, path, file, oss)) {  //look in same path as CHGCAR
        oss << endl;
        oss << soliloquy << "ERROR: Unable to gather new header for file." << endl;
        oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
        oss << soliloquy << "Exiting." << endl;
        oss << endl;
        return oss.str();
      }
      if(!vpflow.flag("PREPARE_CHGCAR_4_JMOL::OUTCAR")) {
        outcar_file = "";  //reset so we can look again later (once per CHGCAR)
      }
    }
    if(!prepare_CHGCAR_4_Jmol(chgcar_file, species_header, zip_file, oss)) {
      oss << endl;
      oss << soliloquy << "ERROR: Unable to rewrite file with appropriate header." << endl;
      oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
      oss << soliloquy << "Exiting." << endl;
      oss << endl;
      return oss.str();
    }
    oss << soliloquy << "finished rewriting " << chgcar_file << "." << endl;
  }
  return oss.str();
}
}  // namespace bader_functions

//***************************************************************************//
// bader_functions::get_species_string
//***************************************************************************//
namespace bader_functions {
bool get_species_string(string& outcar_file, string& species_string, const string& dir_to_look, const string& file, ostream& oss) {
  string soliloquy = "bader_functions::get_species_string():  ";  // so you know who's speaking
  //debug
  bool LDEBUG = (FALSE || XHOST.DEBUG);
  if(LDEBUG) oss << soliloquy << "BEGIN" << endl;

  //CO 180220 - directory stuff for logging
  string directory = ".";  //default
  if(XHOST.vflag_control.flag("DIRECTORY")) {
    directory = XHOST.vflag_control.getattachedscheme("DIRECTORY");  //XHOST.vflag_control.flag("DIRECTORY");//aflags.Directory;//vpflow.getattachedscheme("BADER::DIRECTORY");
  }
  FixDirectory(directory, oss);

  species_string = "";  //clear
  string here = ".";
  //CHECK IF OUTCAR CANNOT BE FOUND
  if(!outcar_file.empty()) {
    if(LDEBUG) oss << soliloquy << "OUTCAR SPECIFIED, NO NEED TO LOOK" << endl;
    if((!aurostd::FileExist(outcar_file)) || aurostd::FileEmpty(outcar_file)) {
      oss << endl;
      oss << soliloquy << "ERROR: Unable to locate specified OUTCAR (or its empty): " << outcar_file << "." << endl;
      oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
      oss << soliloquy << "Exiting." << endl;
      oss << endl;
      return FALSE;
    }
    //FIRST CHECK PATH, OTHERWISE CHECK CURRENT DIRECTORY
    //C++ SHOULD SHORT CIRCUIT HERE to yield correct path
  } else {
    if(LDEBUG) oss << soliloquy << "OUTCAR NOT SPECIFIED, LOOKING FOR IT" << endl;
    //first check in CHGCAR file path
    //otherwise check in current path
    if(!BaderExtensionFound("OUTCAR", outcar_file, dir_to_look) && !BaderExtensionFound("OUTCAR", outcar_file, ".")) {
      oss << endl;
      oss << soliloquy << "ERROR: Unable to locate OUTCAR/OUTCAR.static (compressed or otherwise) in path(" << file << ") or pwd." << endl;
      oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
      oss << soliloquy << "Exiting." << endl;
      oss << endl;
      return FALSE;
    }
  }
  oss << soliloquy << "using " << outcar_file << " to gather species." << endl;
  xOUTCAR outcar(outcar_file);
  for (uint i = 0; i < outcar.species.size(); i++) {
    species_string.append(outcar.species.at(i));
    if(i < outcar.species.size() - 1) { species_string.append(" "); }
  }
  //species_string=species_string.substr(0,species_string.size()-1);  //FIX ME
  oss << soliloquy << "new header: " << species_string << "." << endl;
  return TRUE;
}
}  // namespace bader_functions

//***************************************************************************//
// bader_functions::prepare_CHGCAR_4_Jmol
//***************************************************************************//
namespace bader_functions {
bool prepare_CHGCAR_4_Jmol(const string& _chgcar_file, string& species_header, bool zip_file, ostream& oss) {
  string soliloquy = "bader_functions::prepare_CHGCAR_4_Jmol():  ";  // so you know who's speaking
  //debug
  bool LDEBUG = (FALSE || XHOST.DEBUG);
  if(LDEBUG) oss << soliloquy << "BEGIN" << endl;
  
  //CO 180220 - directory stuff for logging
  string directory = ".";  //default
  if(XHOST.vflag_control.flag("DIRECTORY")) {
    directory = XHOST.vflag_control.getattachedscheme("DIRECTORY");  //XHOST.vflag_control.flag("DIRECTORY");//aflags.Directory;//vpflow.getattachedscheme("BADER::DIRECTORY");
  }
  FixDirectory(directory, oss);

  //overwrite file so header has species
  string chgcar_file = aurostd::CleanFileName(_chgcar_file), output_file;
  if((!aurostd::EFileExist(chgcar_file)) || aurostd::FileEmpty(chgcar_file)) {
    oss << endl;
    oss << soliloquy << "ERROR: Unable to find " << chgcar_file << " (compressed or otherwise, or its empty)." << endl;
    oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
    oss << soliloquy << "Exiting." << endl;
    oss << endl;
    return FALSE;
  }
  oss << soliloquy << "editing header of " << chgcar_file << " for Jmol visualization." << endl;
  aurostd::IsCompressed(chgcar_file, output_file);  //get uncompressed variant
  if(LDEBUG) oss << soliloquy << "UNCOMPRESSED VARIANT OF CHGCAR_FILE IS " << output_file << endl;
  stringstream chgcar_ss;
  aurostd::efile2stringstream(chgcar_file, chgcar_ss);  //get stringstream
  adjust_header(species_header, chgcar_ss, oss);        //fix it
  aurostd::stringstream2file(chgcar_ss, output_file);   //write back out
  if(zip_file && !aurostd::MatchCompressed(chgcar_file, output_file)) {
    oss << endl;
    oss << soliloquy << "ERROR: Unable to zip " << output_file << "." << endl;
    oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
    oss << soliloquy << "Exiting." << endl;
    oss << endl;
    return FALSE;
  }
  oss << soliloquy << "finish editing " << output_file << "." << endl;
  oss << soliloquy << output_file << " should be ready for Jmol." << endl;
  return TRUE;
}
}  // namespace bader_functions

//***************************************************************************//
// bader_functions::prepare_CHGCAR_4_Jmol
//***************************************************************************//
namespace bader_functions {
bool prepare_CHGCAR_4_Jmol(string& _chgcar_file, string& species_header, ostream& oss) {
  bool zip_file = false;
  return prepare_CHGCAR_4_Jmol(_chgcar_file, species_header, zip_file, oss);
}
}  // namespace bader_functions

//***************************************************************************//
// bader_functions::prepare_CHGCAR_4_Jmol
//***************************************************************************//
namespace bader_functions {
bool prepare_CHGCAR_4_Jmol(string& _chgcar_file, string& species_header) {
  ostringstream oss;
  return prepare_CHGCAR_4_Jmol(_chgcar_file, species_header, oss);
}
}  // namespace bader_functions

//[CO 180220 - moved to aurostd]***************************************************************************//
//[CO 180220 - moved to aurostd] bader_functions::efile2tempfile
//[CO 180220 - moved to aurostd]***************************************************************************//
//[CO 180220 - moved to aurostd]namespace bader_functions {
//[CO 180220 - moved to aurostd]}  // namespace bader_functions

//***************************************************************************//
// pflow::CHGCAR2JVXL
//***************************************************************************//
namespace pflow {
string CHGCAR2JVXL(aurostd::xoption& vpflow) {  //corey
  ostringstream oss;
  // handles flags for CHGCAR2JVXL
  string soliloquy = "pflow::CHGCAR2JVXL():  ";        // so you know who's speaking
  string soliloquy_empty = "                       ";  // so you know who's speaking
  //debug
  bool LDEBUG = (FALSE || XHOST.DEBUG);
  if(LDEBUG) oss << soliloquy << "BEGIN" << endl;

  //CO 180220 - directory stuff for logging
  string directory = ".";  //default
  if(XHOST.vflag_control.flag("DIRECTORY")) {
    directory = XHOST.vflag_control.getattachedscheme("DIRECTORY");  //XHOST.vflag_control.flag("DIRECTORY");//aflags.Directory;//vpflow.getattachedscheme("BADER::DIRECTORY");
  }
  bader_functions::FixDirectory(directory, oss);
  
  vector<string> tokens;
  string input = vpflow.getattachedscheme("CHGCAR2JVXL");

  // case 1: cyclic variable
  // 2a:  FILE1,FILE2,FILE3...::CUTOFF[::DOWNSAMPLE]
  // 2b:  FILE::CUTOFF1,CUTOFF2…[::DOWNSAMPLE]
  // 2c:  FILE::CUTOFF[::DOWNSAMPLE1,DOWNSAMPLE2,...]
  // 2d:  some combination of above
  vector<string> chgcar_files, sets;
  vector<double> cutoffs;
  vector<int> downsample_ratios;
  bool cyclic = FALSE;
  if(aurostd::substring2bool(input, "::")) {
    if(LDEBUG) oss << soliloquy << "CYCLIC ROUTINE DETECTED" << endl;
    cyclic = TRUE;
    aurostd::string2tokens(input, tokens, "::");
    if(tokens.size() < 2 || tokens.size() > 3) {
      oss << endl;
      oss << soliloquy << "ERROR: Incorrect format for input - number of tokens (::)." << endl;
      oss << soliloquy << "Needs to be either: " << endl;
      oss << soliloquy_empty << "FILE1,FILE2,FILE3...::CUTOFF[::DOWNSAMPLE], " << endl;
      oss << soliloquy_empty << "FILE::CUTOFF1,CUTOFF2…[::DOWNSAMPLE], or " << endl;
      oss << soliloquy_empty << "FILE::CUTOFF[::DOWNSAMPLE1,DOWNSAMPLE2,...]." << endl;
      oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
      oss << soliloquy << "Exiting." << endl;
      oss << endl;
      return oss.str();
    }
    //make sure there's no mixing
    if(LDEBUG) oss << soliloquy << "CYCLIC ROUTINE DETECTED" << endl;
    for (uint i = 0; i < tokens.size(); i++) {
      if(aurostd::substring2bool(tokens.at(i), ":")) {
        oss << endl;
        oss << soliloquy << "ERROR: Incorrect format for input, cannot specify sets (:) and cyclic (::) parameters." << endl;
        oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
        oss << soliloquy << "Exiting." << endl;
        oss << endl;
        return oss.str();
      }
    }
    if(aurostd::substring2bool(tokens.at(0), ",")) {
      aurostd::string2tokens(tokens.at(0), chgcar_files, ",");
    } else {
      chgcar_files.push_back(tokens.at(0));
    }
    if(aurostd::substring2bool(tokens.at(1), ",")) {
      aurostd::string2tokens<double>(tokens.at(1), cutoffs, ",");
    } else {
      cutoffs.push_back(aurostd::string2utype<double>(tokens.at(1)));
    }
    if(tokens.size() == 3) {
      if(aurostd::substring2bool(tokens.at(2), ",")) {
        aurostd::string2tokens<int>(tokens.at(2), downsample_ratios, ",");
      } else {
        downsample_ratios.push_back(aurostd::string2utype<int>(tokens.at(2)));
      }
    }
    if(downsample_ratios.empty()) {
      if(chgcar_files.size() == 1 && chgcar_files.size() == cutoffs.size()) {
        if(LDEBUG) oss << soliloquy << "CYCLIC ROUTINE IGNORED, ONLY ONE SET" << endl;
        cyclic = FALSE;
      }
    } else {
      if(chgcar_files.size() == 1 && chgcar_files.size() == cutoffs.size() && chgcar_files.size() == downsample_ratios.size()) {
        if(LDEBUG) oss << soliloquy << "CYCLIC ROUTINE IGNORED, ONLY ONE SET" << endl;
        cyclic = FALSE;
      }
    }
    //case 2: sets, FILE1,CUTOFF1[,DOWNSAMPLE1]:FILE2,CUTOFF2[,DOWNSAMPLE2]:...
  } else if(aurostd::substring2bool(input, ":")) {
    if(LDEBUG) oss << soliloquy << "SETS ROUTINE DETECTED" << endl;
    aurostd::string2tokens(input, sets, ":");
    for (uint i = 0; i < sets.size(); i++) {
      aurostd::string2tokens(sets.at(i), tokens, ",");
      if(tokens.size() < 2 || tokens.size() > 3) {
        oss << endl;
        oss << soliloquy << "ERROR: Incorrect format for input " << i + 1 << "." << endl;
        oss << soliloquy << "Needs to be FILE1,CUTOFF1[,DOWNSAMPLE1]:FILE2,CUTOFF2[,DOWNSAMPLE2]:..." << endl;
        oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
        oss << soliloquy << "Exiting." << endl;
        oss << endl;
        return oss.str();
      }
      chgcar_files.push_back(tokens.at(0));
      cutoffs.push_back(aurostd::string2utype<double>(tokens.at(1)));
      if(tokens.size() == 3) { downsample_ratios.push_back(aurostd::string2utype<int>(tokens.at(2))); }
    }
    // case 3: single set
  } else if(aurostd::substring2bool(input, ",")) {
    if(LDEBUG) oss << soliloquy << "SINGLE SET ROUTINE DETECTED" << endl;
    aurostd::string2tokens(input, tokens, ",");
    if(tokens.size() < 2 || tokens.size() > 3) {
      oss << endl;
      oss << soliloquy << "ERROR: Incorrect format for input - number of tokens (,)." << endl;
      oss << soliloquy << "Needs to be: FILE,CUTOFF[,DOWNSAMPLE]." << endl;
      oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
      oss << soliloquy << "Exiting." << endl;
      oss << endl;
      return oss.str();
    }
    chgcar_files.push_back(tokens.at(0));
    cutoffs.push_back(aurostd::string2utype<double>(tokens.at(1)));
    if(tokens.size() == 3) { downsample_ratios.push_back(aurostd::string2utype<int>(tokens.at(2))); }
  } else {
    oss << endl;
    oss << soliloquy << "ERROR: Incorrect input, need at least FILE1,CUTOFF1." << endl;
    oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
    oss << soliloquy << "Exiting." << endl;
    oss << endl;
    return oss.str();
  }

  if(vpflow.flag("CHGCAR2JVXL::OUTPUT")) {
    if(!cyclic) {
      if(LDEBUG) oss << soliloquy << "OUTPUT FILES SPECIFIED" << endl;
      //they must be sets, deal with this case first
      vector<string> output_files;
      string output = vpflow.getattachedscheme("CHGCAR2JVXL::OUTPUT");
      if(aurostd::substring2bool(output, "::")) {
        oss << endl;
        oss << soliloquy << "ERROR: Incorrect input delimiter for output_files." << endl;
        oss << soliloquy << "Please remove any \"::\"." << endl;
        oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
        oss << soliloquy << "Exiting." << endl;
        oss << endl;
        return oss.str();
      }
      if(aurostd::substring2bool(output, ",")) {
        aurostd::string2tokens(output, output_files, ",");
      } else {
        output_files.push_back(output);
      }
      if(!CHGCAR2JVXL(chgcar_files, cutoffs, downsample_ratios, output_files, oss)) {
        oss << endl;
        oss << soliloquy << "ERROR: Unable to convert CHGCAR files to .jvxl's." << endl;
        oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
        oss << soliloquy << "Exiting." << endl;
        oss << endl;
        return oss.str();
      }
      return oss.str();
    } else {
      oss << endl;
      oss << soliloquy << "WARNING: An output file name can only be specified for sets." << endl;
      oss << soliloquy << "You have specified CYCLIC mode (::)." << endl;
      oss << soliloquy << "Ignoring output file name input." << endl;
      oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
      oss << endl;
    }
  }
  if(!CHGCAR2JVXL(chgcar_files, cutoffs, downsample_ratios, cyclic, oss)) {
    oss << endl;
    oss << soliloquy << "ERROR: Unable to convert CHGCAR files to .jvxl's." << endl;
    oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
    oss << soliloquy << "Exiting." << endl;
    oss << endl;
    return oss.str();
  }
  return oss.str();
}
}  // namespace pflow

//***************************************************************************//
// pflow::CHGCAR2JVXL
//***************************************************************************//
namespace pflow {
bool CHGCAR2JVXL(vector<string>& chgcar_files, const vector<double>& cutoffs, const vector<int>& downsample_ratios, const bool& cyclic, ostream& oss) {  //corey
  //NO OUTPUT_FILES SPECIFIED
  vector<string> output_files;
  return CHGCAR2JVXL(chgcar_files, cutoffs, downsample_ratios, output_files, cyclic, oss);
}
}  // namespace pflow

//***************************************************************************//
// pflow::CHGCAR2JVXL
//***************************************************************************//
namespace pflow {
bool CHGCAR2JVXL(vector<string>& chgcar_files,
                 const vector<double>& cutoffs,
                 const vector<int>& downsample_ratios,
                 vector<string>& output_files,
                 const bool& cyclic,
                 ostream& oss) {                 //corey
  string soliloquy = "pflow::CHGCAR2JVXL():  ";  // so you know who's speaking
  //debug
  bool LDEBUG = (FALSE || XHOST.DEBUG);
  if(LDEBUG) oss << soliloquy << "BEGIN" << endl;
  
  //CO 180220 - directory stuff for logging
  string directory = ".";  //default
  if(XHOST.vflag_control.flag("DIRECTORY")) {
    directory = XHOST.vflag_control.getattachedscheme("DIRECTORY");  //XHOST.vflag_control.flag("DIRECTORY");//aflags.Directory;//vpflow.getattachedscheme("BADER::DIRECTORY");
  }
  bader_functions::FixDirectory(directory, oss);

  //NO OUTPUT_FILES SPECIFIED, so have to specify whether to cycle through variables or not
  if(cyclic) {
    if(LDEBUG) oss << soliloquy << "CYCLIC ROUTINE DETECTED" << endl;

    if(downsample_ratios.size() != 0) {
      for (uint i = 0; i < chgcar_files.size(); i++) {
        for (uint j = 0; j < cutoffs.size(); j++) {
          for (uint k = 0; k < downsample_ratios.size(); k++) {
            if(!CHGCAR2JVXL(chgcar_files.at(i), cutoffs.at(j), downsample_ratios.at(k), oss)) {
              oss << endl;
              oss << soliloquy << "ERROR: Unable to convert CHGCAR files to .jvxl's." << endl;
              oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
              oss << soliloquy << "Exiting." << endl;
              oss << endl;
              return FALSE;
            }
          }
        }
      }
    } else {
      for (uint i = 0; i < chgcar_files.size(); i++) {
        for (uint j = 0; j < cutoffs.size(); j++) {
          if(!CHGCAR2JVXL(chgcar_files.at(i), cutoffs.at(j), oss)) {
            oss << endl;
            oss << soliloquy << "ERROR: Unable to convert CHGCAR files to .jvxl's." << endl;
            oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
            oss << soliloquy << "Exiting." << endl;
            oss << endl;
            return FALSE;
          }
        }
      }
    }
  } else {
    if(LDEBUG) oss << soliloquy << "SETS ROUTINE DETECTED" << endl;

    //first check if we have the right number of inputs
    if(downsample_ratios.size() != 0) {
      if(chgcar_files.size() != cutoffs.size() || chgcar_files.size() != downsample_ratios.size()) {
        oss << endl;
        oss << soliloquy << "ERROR: Unequal set parameters found." << endl;
        oss << soliloquy << "Must have equal number of CHGCAR files, cutoffs, and downsample ratios." << endl;
        oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
        oss << soliloquy << "Exiting." << endl;
        oss << endl;
        return FALSE;
      }
    } else {
      if(chgcar_files.size() != cutoffs.size()) {
        oss << endl;
        oss << soliloquy << "ERROR: Unequal set parameters found." << endl;
        oss << soliloquy << "Must have equal number of CHGCAR files and cutoffs." << endl;
        oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
        oss << soliloquy << "Exiting." << endl;
        oss << endl;
        return FALSE;
      }
    }

    //figure out if we need output_filenames, otherwise check if we have the right size
    if(output_files.empty()) {
      for (uint i = 0; i < chgcar_files.size(); i++) {
        if(downsample_ratios.size() != 0) {
          output_files.push_back(CHGCAR2JVXL_get_output_filename(chgcar_files.at(i), cutoffs.at(i), downsample_ratios.at(i)));
        } else {
          output_files.push_back(CHGCAR2JVXL_get_output_filename(chgcar_files.at(i), cutoffs.at(i)));
        }
      }
    } else {
      if(downsample_ratios.size() != 0) {
        if(chgcar_files.size() != cutoffs.size() || chgcar_files.size() != downsample_ratios.size() || chgcar_files.size() != output_files.size()) {
          oss << endl;
          oss << soliloquy << "ERROR: Unequal set parameters found." << endl;
          oss << soliloquy << "Number of OUTPUT files must match number of CHGCAR files, cutoffs, and downsample ratios." << endl;
          oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
          oss << soliloquy << "Exiting." << endl;
          oss << endl;
          return FALSE;
        }
      } else {
        //downsample_ratios.size()==0
        if(chgcar_files.size() != cutoffs.size() || chgcar_files.size() != output_files.size()) {
          oss << endl;
          oss << soliloquy << "ERROR: Unequal set parameters found." << endl;
          oss << soliloquy << "Number of OUTPUT files must match number of CHGCAR files and cutoffs." << endl;
          oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
          oss << soliloquy << "Exiting." << endl;
          oss << endl;
          return FALSE;
        }
      }
    }

    //call function
    if(downsample_ratios.size() != 0) {
      for (uint i = 0; i < chgcar_files.size(); i++) {
        if(!CHGCAR2JVXL(chgcar_files.at(i), cutoffs.at(i), downsample_ratios.at(i), output_files.at(i), oss)) {
          oss << endl;
          oss << soliloquy << "ERROR: Unable to convert CHGCAR files to .jvxl's." << endl;
          oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
          oss << soliloquy << "Exiting." << endl;
          oss << endl;
          return FALSE;
        }
      }
      return TRUE;
    }
    //downsample_ratios.size()==0
    for (uint i = 0; i < chgcar_files.size(); i++) {
      if(!CHGCAR2JVXL(chgcar_files.at(i), cutoffs.at(i), output_files.at(i), oss)) {
        oss << endl;
        oss << soliloquy << "ERROR: Unable to convert CHGCAR files to .jvxl's." << endl;
        oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
        oss << soliloquy << "Exiting." << endl;
        oss << endl;
        return FALSE;
      }
    }
    return TRUE;
  }
  return TRUE;
}
}  // namespace pflow

//***************************************************************************//
// pflow::CHGCAR2JVXL
//***************************************************************************//
namespace pflow {
bool CHGCAR2JVXL(vector<string>& chgcar_files, const vector<double>& cutoffs, const vector<int>& downsample_ratios, ostream& oss) {  //corey
  //NO CYCLIC VARIABLE, so assume sets
  bool cyclic = FALSE;
  return CHGCAR2JVXL(chgcar_files, cutoffs, downsample_ratios, cyclic, oss);
}
}  // namespace pflow

//***************************************************************************//
// pflow::CHGCAR2JVXL
//***************************************************************************//
namespace pflow {
bool CHGCAR2JVXL(vector<string>& chgcar_files, const vector<double>& cutoffs, const bool& cyclic, ostream& oss) {  //corey
  //NO DOWNSAMPLE
  vector<int> downsample_ratios;
  return CHGCAR2JVXL(chgcar_files, cutoffs, downsample_ratios, cyclic, oss);
}
}  // namespace pflow

//***************************************************************************//
// pflow::CHGCAR2JVXL
//***************************************************************************//
namespace pflow {
bool CHGCAR2JVXL(vector<string>& chgcar_files, const vector<double>& cutoffs, ostream& oss) {  //corey
  //NO DOWNSAMPLE OR CYCLIC
  bool cyclic = FALSE;
  return CHGCAR2JVXL(chgcar_files, cutoffs, cyclic, oss);
}
}  // namespace pflow

//***************************************************************************//
// pflow::CHGCAR2JVXL
//***************************************************************************//
namespace pflow {
bool CHGCAR2JVXL(vector<string>& chgcar_files,
                 const vector<double>& cutoffs,
                 const vector<int>& downsample_ratios,
                 vector<string>& output_files,
                 ostringstream& oss) {  //corey
  //OUTPUT_FILES SPECIFIED, so must be sets
  bool cyclic = FALSE;
  return CHGCAR2JVXL(chgcar_files, cutoffs, downsample_ratios, output_files, cyclic, oss);
}
}  // namespace pflow

//***************************************************************************//
// pflow::CHGCAR2JVXL
//***************************************************************************//
namespace pflow {
bool CHGCAR2JVXL(vector<string>& chgcar_files,
                 const vector<double>& cutoffs,
                 vector<string>& output_files,
                 ostringstream& oss) {  //corey
  //NO DOWNSAMPLE SPECIFIED
  vector<int> downsample_ratios;
  return CHGCAR2JVXL(chgcar_files, cutoffs, downsample_ratios, output_files, oss);
}
}  // namespace pflow

//***************************************************************************//
// pflow::CHGCAR2JVXL
//***************************************************************************//
namespace pflow {
bool CHGCAR2JVXL(string chgcar_file, const double cutoff, const int downsample_ratio, string output_file, ostream& oss) {
  //ACTUALLY HANDLES JMOL INTERACTION
  string soliloquy = "pflow::CHGCAR2JVXL():  ";        // so you know who's speaking
  string soliloquy_empty = "                       ";  // so you know who's speaking
  bool LDEBUG = (FALSE || XHOST.DEBUG);
  if(LDEBUG) oss << soliloquy << "BEGIN" << endl;
  
  //CO 180220 - directory stuff for logging
  string directory = ".";  //default
  if(XHOST.vflag_control.flag("DIRECTORY")) {
    directory = XHOST.vflag_control.getattachedscheme("DIRECTORY");  //XHOST.vflag_control.flag("DIRECTORY");//aflags.Directory;//vpflow.getattachedscheme("BADER::DIRECTORY");
  }
  bader_functions::FixDirectory(directory, oss);
  
  //string command="java -jar Jmol.jar -ionx";
  string command;
  //Look for JmolData.jar first
  if((XHOST.hostname == "nietzsche.mems.duke.edu" || XHOST.hostname == "aflowlib.mems.duke.edu") && aurostd::FileExist("/usr/local/bin/JmolData.jar")) {
    command = "java -jar /usr/local/bin/JmolData.jar";
    //command="/usr/local/bin/jmol";
  } else if(aurostd::IsCommandAvailable("JmolData.jar")) {
    command = "java -jar " + XHOST.command("JmolData.jar");
    //Look for jmol instead
  } else if(aurostd::IsCommandAvailable("jmol")) {
    command = XHOST.command("jmol");
  } else if(aurostd::IsCommandAvailable("jmol.sh")) {
    command = XHOST.command("jmol.sh");
  } else if(aurostd::IsCommandAvailable("Jmol.jar")) {
    command = "java -jar " + XHOST.command("Jmol.jar");
  } else {
    oss << endl;
    oss << soliloquy << "ERROR: Unable to locate either JmolData (preferred), jmol, jmol.sh, or Jmol.jar in path." << endl;
    oss << soliloquy << "Please install Jmol." << endl;
    oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
    oss << soliloquy << "Exiting." << endl;
    oss << endl;
    return FALSE;
  }
  oss << soliloquy << "using " + command + " to run Jmol." << endl;
  command.append(" -ionx");

  //check if files exist/no negatives are present, throw error
  if(chgcar_file.empty() || (!aurostd::FileExist(chgcar_file)) || aurostd::FileEmpty(chgcar_file)) {
    oss << endl;
    oss << soliloquy << "ERROR: Input CHGCAR file does not exist (or its empty): " << chgcar_file << "." << endl;
    oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
    oss << soliloquy << "Exiting." << endl;
    oss << endl;
    return FALSE;
  }
  if(cutoff < 0) {
    oss << endl;
    oss << soliloquy << "ERROR: Cutoff cannot be negative: " << cutoff << "." << endl;
    oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
    oss << soliloquy << "Exiting." << endl;
    oss << endl;
    return FALSE;
  }
  if(cutoff == 0) {
    oss << endl;
    oss << soliloquy << "ERROR: Input for cutoff be numerical and above 0: " << cutoff << "." << endl;
    oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
    oss << soliloquy << "Exiting." << endl;
    oss << endl;
    return FALSE;
  }
  if(downsample_ratio < 0) {
    oss << endl;
    oss << soliloquy << "ERROR: Downsample ratio cannot be negative: " << downsample_ratio << "." << endl;
    oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
    oss << soliloquy << "Exiting." << endl;
    oss << endl;
    return FALSE;
  }
  if(downsample_ratio == 0) {
    oss << endl;
    oss << soliloquy << "ERROR: Input for downsample ratio be numerical and above 0: " << downsample_ratio << "." << endl;
    oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
    oss << soliloquy << "Exiting." << endl;
    oss << endl;
    return FALSE;
  }
  if(aurostd::FileExist(output_file)) {
    oss << endl;
    oss << soliloquy << "WARNING: Output file " << output_file << " EXISTS." << endl;
    oss << soliloquy << "WARNING: Overwriting output file " << output_file << "." << endl;
    oss << endl;
  }
  if(output_file.substr(output_file.length() - 5, 5) != ".jvxl") {
    oss << endl;
    oss << soliloquy << "ERROR: OUTPUT file must have .jvxl extension." << endl;
    oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
    oss << soliloquy << "Exiting." << endl;
    oss << endl;
    return FALSE;
  }
  if(LDEBUG) oss << soliloquy << "INPUTS SEEM TO BE ENTERED CORRECTLY" << endl;

  bool clean_up = false;
  //check if INPUT is compressed, and decompress
  if(LDEBUG) oss << soliloquy << "FIND / DECOMPRESS INPUT" << endl;
  if(aurostd::IsCompressed(chgcar_file)) {
    clean_up = true;
    oss << soliloquy << "decompressing " << chgcar_file << "." << endl;
    if(!aurostd::efile2tempfile(chgcar_file, chgcar_file)) {
      oss << endl;
      oss << soliloquy << "ERROR: Unable to decompress " << chgcar_file << "." << endl;
      oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
      oss << soliloquy << "Exiting." << endl;
      oss << endl;
      return FALSE;
    }
    oss << soliloquy << "created temporary file " << chgcar_file << "." << endl;
  }
  //NO LONGER NEEDED IF WORKING IN NO DISPLAY VERSION
  //if(!jmoldata) {
  //  command.append(" \""+chgcar_file+"\"");
  //}
  command.append(" -j 'isosurface");
  if(downsample_ratio != 1) {  //already took care of 0
    command.append(" downsample " + aurostd::utype2string(downsample_ratio));
  }
  //command.append(" cutoff "+aurostd::utype2string(cutoff)+" \""+chgcar_file+"\"; write ISOSURFACE");
  stringstream num_prec;
  num_prec << cutoff;
  command.append(" cutoff " + num_prec.str() + " \"" + chgcar_file + "\"; write");
  num_prec.str("");
  command.append(" \"" + output_file + "\";'");
  if(LDEBUG) oss << soliloquy << "Jmol command: " << command << endl;
  string jmol_output;
  jmol_output = aurostd::execute2string(command);
  //jmol_output=aurostd::RemoveWhiteSpaces(aurostd::execute2string(command)); //no longer needed
  if(LDEBUG) oss << soliloquy << "Jmol output:  " << endl;
  if(LDEBUG) oss << jmol_output;
  //I WOULD LIKE TO DO THIS, BUT JMOL'S OUTPUT MAY CHANGE WITH DIFFERENT VERSIONS, SO JUST CHECK FOR JVXL FILE
  //check output and that jvxl file exists
  /*if(!(aurostd::substring2bool(jmol_output,"isosurface1created")&&aurostd::substring2bool(jmol_output,"OKXJVXL")&&
              aurostd::FileExist(output_file))) {
          oss << endl;
          oss << soliloquy << "ERROR: Jmol unable to create .jvxl file." << endl;
          oss << soliloquy << "Please correctly set your display variable and check your version of JmolData/Jmol." << endl;
          oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
          oss << soliloquy << "Exiting." << endl;
          oss << endl;
          aurostd::RemoveFile(chgcar_file);
          return FALSE;
        }*/
  //isosurface1 created
  //NO LONGER APPLICABLE
  /*if(jmol_output!="-5") {    //expected command from Jmol, if it doesn't come out, then likely display variable off
          oss << endl;
          oss << soliloquy << "ERROR: Jmol unable to create .jvxl file." << endl;
          oss << soliloquy << "Please correctly set your display variable." << endl;
          oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
          oss << soliloquy << "Exiting." << endl;
          oss << endl;
          aurostd::RemoveFile(chgcar_file);
          return FALSE;
        }*/
  if((!aurostd::FileExist(output_file)) || aurostd::FileEmpty(output_file)) {
    oss << endl;
    oss << soliloquy << "ERROR: Jmol unable to create .jvxl file." << endl;
    oss << soliloquy << "Please correctly set your display variable and check your version of JmolData/Jmol." << endl;
    oss << soliloquy << "[dir=" << directory << "]" << endl; //CO 180220
    oss << soliloquy << "Exiting." << endl;
    oss << endl;
    aurostd::RemoveFile(chgcar_file);
    return FALSE;
  }
  oss << soliloquy << output_file << " successfully created." << endl;
  //REMOVE TEMP CHGCAR_FILE
  if(clean_up) {
    oss << soliloquy << "removing temporary file " << chgcar_file << "." << endl;
    aurostd::RemoveFile(chgcar_file);
  }
  return TRUE;
}
}  // namespace pflow

//***************************************************************************//
// pflow::CHGCAR2JVXL
//***************************************************************************//
namespace pflow {
bool CHGCAR2JVXL(string chgcar_file, const double& cutoff, const int& downsample_ratio, ostream& oss) {
  //MISSING OUTPUT
  string output_file = CHGCAR2JVXL_get_output_filename(chgcar_file, cutoff, downsample_ratio);
  return CHGCAR2JVXL(chgcar_file, cutoff, downsample_ratio, output_file, oss);
}
}  // namespace pflow

//***************************************************************************//
// pflow::CHGCAR2JVXL
//***************************************************************************//
namespace pflow {
bool CHGCAR2JVXL(string chgcar_file, const double& cutoff, string& output_file, ostream& oss) {
  //MISSING DOWNSAMPLE
  int downsample_ratio = 1;
  return CHGCAR2JVXL(chgcar_file, cutoff, downsample_ratio, output_file, oss);
}
}  // namespace pflow

//***************************************************************************//
// pflow::CHGCAR2JVXL
//***************************************************************************//
namespace pflow {
bool CHGCAR2JVXL(string chgcar_file, const double& cutoff, ostream& oss) {
  //MISSING DOWNSAMPLE AND OUTPUT
  int downsample_ratio = 1;
  return CHGCAR2JVXL(chgcar_file, cutoff, downsample_ratio, oss);
}
}  // namespace pflow

//***************************************************************************//
// pflow::CHGCAR2JVXL_get_output_filename
//***************************************************************************//
namespace pflow {
string CHGCAR2JVXL_get_output_filename(string chgcar_file, const double& cutoff, const int& downsample_ratio) {
  //for handling compressed files, e.g. we don't want .bz2 to be part of output name
  string FileNameIn_decompressed;
  aurostd::IsCompressed(chgcar_file, FileNameIn_decompressed);
  stringstream num_prec;
  num_prec << std::fixed << setprecision(2) << cutoff;
  string cutoff_string = num_prec.str();  //avoid issues with setprecision(20) in string2utype<double>
  //string cutoff_string=aurostd::utype2string(cutoff);
  if(aurostd::substring2bool(cutoff_string, "0.")) { cutoff_string = aurostd::StringSubst(cutoff_string, "0.", ""); }
  if(aurostd::substring2bool(cutoff_string, ".")) { cutoff_string = aurostd::StringSubst(cutoff_string, ".", ""); }
  string output_file = FileNameIn_decompressed + "_Bader_" + cutoff_string;
  if(downsample_ratio != 0) {
    output_file.append("_" + aurostd::utype2string(downsample_ratio));
  }
  output_file.append(".jvxl");
  return output_file;
}
}  // namespace pflow

//***************************************************************************//
// pflow::CHGCAR2JVXL_get_output_filename
//***************************************************************************//
namespace pflow {
string CHGCAR2JVXL_get_output_filename(string chgcar_file, const double& cutoff) {
  //NO DOWNSAMPLE
  int downsample_ratio = 0;
  return CHGCAR2JVXL_get_output_filename(chgcar_file, cutoff, downsample_ratio);
}
}  // namespace pflow

//***************************************************************************//
// pflow::PLOT_PHDISP
//***************************************************************************//
namespace pflow {
bool PLOT_PHDISP(vector<string>& argv) {
  // Generate phonon dispersion plot using GNUPLOT.  Usage:
  // >>aflow --plotphonondispersion | --pphdis [--print=eps | --print=pdf | --print=png | --print=jpg | --print=gif ]  [ --hz | --THz | --rcm | --meV ] DIRECTORY

  // directory is expected as last argument
  string directory = argv.at(argv.size() - 1);

  // exchange "./" for an explicit path
  if(directory == "./") {
    directory = aurostd::execute2string("pwd");
  }

  // announce the directory to be searched for phonon data
  cout << "pflow::PLOT_PHDISP():  working on " << directory << endl;

  // Process output unit selection arguments and initialize a bool for each unit
  bool printHZ = aurostd::args2flag(argv, "--hz|-hz");
  bool printTHZ = aurostd::args2flag(argv, "--THz|-THz");
  bool printRCM = aurostd::args2flag(argv, "--rcm|-rcm");
  bool printMEV = aurostd::args2flag(argv, "--meV|-meV");

  // set and declare default units if none are specified
  if((!printHZ && !printTHZ) && (!printRCM && !printMEV)) {
    cout << "pflow::PLOT_PHDISP:  No output units specified; output will be in units of reciprocal centimeters." << endl;
    printRCM = true;
  }

  // Process output image format selection arguments and initialize a bool for each except .eps
  // optional forms will be converted from .eps

  // announce default behavior if no other selection is made
  if((!XHOST.vflag_control.flag("PRINT_MODE::PDF") && !XHOST.vflag_control.flag("PRINT_MODE::PNG")) && (!XHOST.vflag_control.flag("PRINT_MODE::JPG") && !XHOST.vflag_control.flag("PRINT_MODE::GIF"))) {
    cout << "pflow::PLOT_PHDISP:  No output format specified; .eps output will not be converted." << endl;
  }

  // Use output unit bools to prepare string vectors
  vector<IPCFreqFlags> outFlag;  // implemented in APL/ipc.h
  vector<string> yLabelString;
  vector<string> unitName;
  if(printHZ) {
    outFlag.push_back(HERTZ);
    yLabelString.push_back("Hz");
    unitName.push_back("hz");
  }
  if(printTHZ) {
    outFlag.push_back(THZ);
    yLabelString.push_back("THz");
    unitName.push_back("THz");
  }
  if(printRCM) {
    outFlag.push_back(RECIPROCAL_CM);
    yLabelString.push_back("cm^{-1}");
    unitName.push_back("rcm");
  }
  if(printMEV) {
    outFlag.push_back(MEV);
    yLabelString.push_back("meV");
    unitName.push_back("meV");
  }

  // Expected name of the data files
  string phondisdat = directory + "/PDIS";
  string phondosdat = directory + "/PDOS";
  string phonthrmdat = directory + "/THERMO";

  // Record which files are present and if they were compressed
  bool dis_data = false, dos_data = false, thermo_data = false;
  bool pdis_was_decompressed = false, pdos_was_decompressed = false, thermo_was_decompressed = false;

  // Search for the data files, if not found then search for compressed version
  string PDIS_tmp, PDOS_tmp, THERMO_tmp;  // string to hold path to temporary file
  // look for PDIS
  if(!aurostd::FileExist(phondisdat)) {
    if(aurostd::FileExist(phondisdat + ".bz2")) {
      cout << "pflow::PLOT_PHDISP:  Decompressing PDIS.bz2 . . . " << endl;
      PDIS_tmp = aurostd::TmpFileCreate("PDIS_temp");
      aurostd::execute("bzcat " + phondisdat + ".bz2 > " + PDIS_tmp);
      pdis_was_decompressed = true;
      dis_data = true;
      phondisdat = PDIS_tmp;
    } else if(aurostd::FileExist(phondisdat + ".gz")) {
      cout << "pflow::PLOT_PHDISP:  Decompressing PDIS.gz . . . " << endl;
      PDIS_tmp = aurostd::TmpFileCreate("PDIS_temp");
      aurostd::execute("zcat " + phondisdat + ".gz > " + PDIS_tmp);
      pdis_was_decompressed = true;
      dis_data = true;
      phondisdat = PDIS_tmp;
    } else {
      cout << "pflow::PLOT_PHDISP:  AFLOW found no file named PDIS in " << directory << "." << endl;
    }
  } else {
    dis_data = true;
  }
  // look for PDOS
  if(!aurostd::FileExist(phondosdat)) {
    if(aurostd::FileExist(phondosdat + ".bz2")) {
      cout << "pflow::PLOT_PHDISP:  Decompressing PDOS.bz2 . . ." << endl;
      PDOS_tmp = aurostd::TmpFileCreate("PDOS_temp");
      aurostd::execute("bzcat " + phondosdat + ".bz2 > " + PDOS_tmp);
      pdos_was_decompressed = true;
      dos_data = true;
      phondosdat = PDOS_tmp;
    } else if(aurostd::FileExist(phondosdat + ".gz")) {
      cout << "pflow::PLOT_PHDISP:  Decompressing PDOS.gz . . ." << endl;
      PDOS_tmp = aurostd::TmpFileCreate("PDOS_temp");
      aurostd::execute("zcat " + phondosdat + ".gz > " + PDOS_tmp);
      pdos_was_decompressed = true;
      dos_data = true;
      phondosdat = PDOS_tmp;
    } else {
      cout << "pflow::PLOT_PHDISP:  AFLOW found no file named PDOS in " << directory << "." << endl;
    }
  } else {
    dos_data = true;
  }
  // look for THERMO
  if(!aurostd::FileExist(phonthrmdat)) {
    if(aurostd::FileExist(phonthrmdat + ".bz2")) {
      cout << "pflow::PLOT_PHDISP:  Decompressing THERMO.bz2 . . ." << endl;
      THERMO_tmp = aurostd::TmpFileCreate("THERMO_temp");
      aurostd::execute("bzcat " + phonthrmdat + ".bz2 > " + THERMO_tmp);
      thermo_was_decompressed = true;
      thermo_data = true;
      phonthrmdat = THERMO_tmp;
    } else if(aurostd::FileExist(phonthrmdat + ".gz")) {
      cout << "pflow::PLOT_PHDISP:  Decompressing THERMO.gz . . ." << endl;
      THERMO_tmp = aurostd::TmpFileCreate("THERMO_temp");
      aurostd::execute("zcat " + phonthrmdat + ".gz > " + THERMO_tmp);
      thermo_was_decompressed = true;
      thermo_data = true;
      phonthrmdat = THERMO_tmp;
    } else {
      cout << "pflow::PLOT_PHDISP:  AFLOW found no file named THERMO in " << directory << "." << endl;
    }
  } else {
    thermo_data = true;
  }

  // prepare variables for storing calculation details contained in PDIS
  IPCFreqFlags freqUnits;
  int nbranches;
  vector<string> tokens;
  string chemName = "None";
  string graph_title = "None";
  vector<double> specialTickPositions;
  vector<string> specialTickLabels;
  //vector<double> exactQPointPositions;

  // read PDIS and parse the header
  ifstream fin;
  fin.open(phondisdat.c_str());
  if(fin.is_open()) {
    string line;

    // Read the header
    while (true) {
      getline(fin, line);
      if(fin.eof()) break;        // does the file have contents?
      if(line.empty()) continue;  // was the first line blank?
      if(line[0] != '#') break;   // is the first character #?

      // Read the name of the studied system
      if(line.find("<system>") != string::npos)  // was "<system>" in the line?
      {
        string temp;
        tokenize(line, tokens, string("\""));
        temp = tokens[1].c_str();  // temp is all strings in between " 's
        tokens.clear();
        tokenize(temp, tokens, " ");  // divide by whitespace
        chemName = tokens[0].c_str();
        temp = tokens[0] + " " + tokens[1];  // join chemical and crystal strings
        graph_title = temp.c_str();
        tokens.clear();
        continue;
      }
      // Read the units
      if(line.find("<units>") != string::npos)  // was "<units>" in the line?
      {
        tokenize(line, tokens, string(" "));
        freqUnits = (IPCFreqFlags)atoi(tokens[2].c_str());
        tokens.clear();
        continue;
      }
      // Read number of phonon branches
      if(line.find("<nbranches>") != string::npos)  // was "<nbranches>" in the line?
      {
        // Get value
        tokenize(line, tokens, string(" "));
        nbranches = atoi(tokens[2].c_str());
        tokens.clear();
        continue;
      }
      // Read the labels
      if(line.find("<label>") != string::npos)  // was "<label>" in the line?
      {
        tokenize(line, tokens, string(" "));
        specialTickPositions.push_back(atof(tokens[2].c_str()));
        string label = tokens[3];
        if(label == string("G") || label == string("Gamma") ||
            label == string("\\Gamma"))
          label = string("{/Symbol G}");
        specialTickLabels.push_back(label);
        tokens.clear();
        continue;
      }
      /*
            // Read positions of exact qpoints
            if( line.find("<exact>") != string::npos ) // was "<exact>" in the line?
            {
            tokenize(line, tokens, string(" "));
            exactQPointPositions.push_back( atof( tokens[2].c_str() ) );
            tokens.clear();
            continue;
            }
            */
    }  // End: Read the header
  }    // End: if( fin.open() )
  fin.close();

  // create a list of conversion factors now that we know the units of the data
  vector<double> conversionFactor;
  for (unsigned int a = 0; a < outFlag.size(); ++a) {
    // Conversion from eV/A/A/atomic_mass_unit -> something
    if(freqUnits & apl::RAW) {
      //cout << "Input data is in RAW form." << endl;
      if(outFlag.at(a) & apl::RAW) {
        // Transform to eV/A/A/atomic_mass_unit
        conversionFactor.push_back(1.0);
      } else if(outFlag.at(a) & apl::HERTZ) {
        // Transform to s-1; sqrt(EV_TO_JOULE / (ANGSTROM_TO_METER*ANGSTROM_TO_METER) / AMU_TO_KG);
        conversionFactor.push_back(0.98226977255434387350E14);
      } else if(outFlag.at(a) & apl::THZ) {
        // Transform to THz; (in Herts) / 1E12;
        conversionFactor.push_back(98.226977255434387350);
      } else if(outFlag.at(a) & apl::RECIPROCAL_CM) {
        // Transform to cm-1; 1/lambda(m) = freq.(s-1) / light_speed(m/s);
        conversionFactor.push_back(1E-2 * 0.98226977255434387350E14 / 2.99792458E8);
      } else if(outFlag.at(a) & apl::MEV) {
        // Transform to meV; E(eV) = h(eV.s) * freq(s-1); h[(from J.s) -> (eV.s)] = 4.1356673310E-15
        conversionFactor.push_back(0.98226977255434387350E14 * 4.1356673310E-15 / 1E-3);
      }
    }  // End:  RAW
    // Conversion from Hz -> something
    else if(freqUnits & apl::HERTZ) {
      if(outFlag.at(a) & apl::HERTZ) {
        conversionFactor.push_back(1.0);
      } else if(outFlag.at(a) & apl::THZ) {
        conversionFactor.push_back(1E-12);
      } else if(outFlag.at(a) & apl::RECIPROCAL_CM) {
        conversionFactor.push_back(1E-12 * 33.355640952);
      } else if(outFlag.at(a) & apl::MEV) {
        conversionFactor.push_back(1E-12 * 4.1356673310);
      }
    }  // End:  HERTZ
    // Conversion from THz -> something
    else if(freqUnits & apl::THZ) {
      //cout << "Input data is in THZ form." << endl;
      if(outFlag.at(a) & apl::RAW) {
        // Transform to eV/A/A/atomic_mass_unit
        conversionFactor.push_back(1.0 / 98.226977255434387350);
      } else if(outFlag.at(a) & apl::HERTZ) {
        conversionFactor.push_back(1E12);
      } else if(outFlag.at(a) & apl::THZ) {
        conversionFactor.push_back(1.0);
      } else if(outFlag.at(a) & apl::RECIPROCAL_CM) {
        conversionFactor.push_back(33.355640952);  // http://www.lsbu.ac.uk/water/data3.html
      } else if(outFlag.at(a) & apl::MEV) {
        conversionFactor.push_back(4.1356673310);
      }
    }  //  End:  THZ
    // Conversion from cm-1 --> something
    else if(freqUnits & apl::RECIPROCAL_CM) {
      if(outFlag.at(a) & apl::HERTZ) {
        conversionFactor.push_back(1.0 / (1E-12 * 33.355640952));
      } else if(outFlag.at(a) & apl::THZ) {
        conversionFactor.push_back(1.0 / (33.355640952));
      } else if(outFlag.at(a) & apl::RECIPROCAL_CM) {
        conversionFactor.push_back(1.0);
      } else if(outFlag.at(a) & apl::MEV) {
        conversionFactor.push_back(0.12398419055269108);  // http://www.kbfi.ee/THz/?p=energyconverter
      }
    }  // End:  cm-1
    // Conversion from meV --> something
    else if(freqUnits & apl::MEV) {
      if(outFlag.at(a) & apl::HERTZ) {
        conversionFactor.push_back(1.0 / (1E-12 * 4.1356673310));
      } else if(outFlag.at(a) & apl::THZ) {
        conversionFactor.push_back(1.0 / 4.1356673310);
      } else if(outFlag.at(a) & apl::RECIPROCAL_CM) {
        conversionFactor.push_back(1.0 / 0.12398419055269108);
      } else if(outFlag.at(a) & apl::MEV) {
        conversionFactor.push_back(1.0);
      }
    }  // End:  meV
  }    // End:  for loop over requested units outFlags

  //
  //***********************************GENERATING GNUPLOT SCRIPTS******************************************************
  // write PHON, the quintuplet plot of all data
  if(dis_data && dos_data && thermo_data) {
    for (unsigned int b = 0; b < unitName.size(); ++b) {
      string plot_type = "PHON";
      string gnuplotscript_name = directory + "/" + aurostd::ProperFileName("GNUPLOT_" + chemName + "_" + unitName.at(b) + "_PHON.gp");
      ofstream fout;
      fout.open(gnuplotscript_name.c_str());
      if(fout.is_open()) {
        cout << "pflow::PLOT_PHDISP:  Writing GNUPLOT script for phonon dispersion to " << gnuplotscript_name << endl;
      } else {
        cout << "pflow::PLOT_PHDISP:  ERROR!  Couldn't open " << gnuplotscript_name.c_str() << " for writing." << endl;
        return FALSE;
      }
      fout << gnuplotPhononHeader(directory, chemName, unitName.at(b), plot_type, graph_title);
      fout << gnuplotPhononDos(phondosdat, "0.82", "0.40");
      fout << gnuplotPhononDisp(phondisdat, yLabelString.at(b), conversionFactor.at(b), specialTickLabels, specialTickPositions, nbranches, "0.80", "0.40");
      fout << gnuplotPhononThermo(phonthrmdat, "0.34");
      fout.close();

      // Call gnuplot to plot the bands structure
      aurostd::execute(XHOST.command("gnuplot") + " " + gnuplotscript_name);

      // call any conversion
      callOutputConversion(XHOST.vflag_control.flag("PRINT_MODE::PNG"),
                           XHOST.vflag_control.flag("PRINT_MODE::PDF"),
                           XHOST.vflag_control.flag("PRINT_MODE::JPG"),
                           XHOST.vflag_control.flag("PRINT_MODE::GIF"),
                           directory, chemName, unitName.at(b), plot_type);
      // remove gnuplot scripts if they point to temporary files
      cleanGnuplotScript(pdis_was_decompressed, pdos_was_decompressed, thermo_was_decompressed, gnuplotscript_name);
    }  // End:  for loop over requested unit types
  }    // End:  if all data is present write PHON

  // write PDISDOS, the double plot of dispersion and DoS
  if(dis_data && dos_data) {
    for (unsigned int f = 0; f < unitName.size(); ++f) {
      string plot_type = "PDISDOS";
      string gnuplotscript_name = directory + "/" + aurostd::ProperFileName("GNUPLOT_" + chemName + "_" + unitName.at(f) + "_PDISDOS.gp");
      ofstream fout;
      fout.open(gnuplotscript_name.c_str());
      if(fout.is_open()) {
        cout << "pflow::PLOT_PHDISP:  Writing GNUPLOT script for phonon dispersion and DoSC to " << gnuplotscript_name << endl;
      } else {
        cout << "pflow::PLOT_PHDISP:  ERROR!  Couldn't open " << gnuplotscript_name.c_str() << " for writing." << endl;
        return FALSE;
      }
      fout << gnuplotPhononHeader(directory, chemName, unitName.at(f), plot_type, graph_title);
      fout << gnuplotPhononDos(phondosdat, "0.82", "0.15");
      fout << gnuplotPhononDisp(phondisdat, yLabelString.at(f), conversionFactor.at(f), specialTickLabels, specialTickPositions, nbranches, "0.80", "0.15");
      fout << endl;
      fout << "unset multiplot" << endl;
      fout.close();

      // Call gnuplot to plot
      aurostd::execute(XHOST.command("gnuplot") + " " + gnuplotscript_name);

      // call any conversion
      callOutputConversion(XHOST.vflag_control.flag("PRINT_MODE::PNG"),
                           XHOST.vflag_control.flag("PRINT_MODE::PDF"),
                           XHOST.vflag_control.flag("PRINT_MODE::JPG"),
                           XHOST.vflag_control.flag("PRINT_MODE::GIF"),
                           directory, chemName, unitName.at(f), plot_type);

      // remove gnuplot scripts if they point to temporary files
      cleanGnuplotScript(pdis_was_decompressed, pdos_was_decompressed, thermo_was_decompressed, gnuplotscript_name);
    }  // End:  for loop over requested unit types
  }    // End:  if PDIS and PDOS are present

  // write PDIS graph, the phonon dispersion
  if(dis_data) {
    for (unsigned int f = 0; f < unitName.size(); ++f) {
      // Create scriptname from chemical name
      string plot_type = "PDIS";
      string gnuplotscript_name = directory + "/" + aurostd::ProperFileName("GNUPLOT_" + chemName + "_" + unitName.at(f) + "_PDIS.gp");
      ofstream fout;
      fout.open(gnuplotscript_name.c_str());
      if(fout.is_open()) {
        cout << "pflow::PLOT_PHDISP:  Writing GNUPLOT script for phonon dispersion to " << gnuplotscript_name << endl;
      } else {
        cout << "pflow::PLOT_PHDISP:  ERROR!  Couldn't open " << gnuplotscript_name.c_str() << " for writing." << endl;
        return FALSE;
      }
      fout << gnuplotPhononHeader(directory, chemName, unitName.at(f), plot_type, graph_title);
      fout << gnuplotPhononDisp(phondisdat, yLabelString.at(f), conversionFactor.at(f), specialTickLabels, specialTickPositions, nbranches, "0.92", "0.15");
      fout.close();
      // Call gnuplot to plot
      aurostd::execute(XHOST.command("gnuplot") + " " + gnuplotscript_name);

      // call any conversion
      callOutputConversion(XHOST.vflag_control.flag("PRINT_MODE::PNG"),
                           XHOST.vflag_control.flag("PRINT_MODE::PDF"),
                           XHOST.vflag_control.flag("PRINT_MODE::JPG"),
                           XHOST.vflag_control.flag("PRINT_MODE::GIF"),
                           directory, chemName, unitName.at(f), plot_type);

      // remove gnuplot scripts if they point to temporary files
      cleanGnuplotScript(pdis_was_decompressed, pdos_was_decompressed, thermo_was_decompressed, gnuplotscript_name);
    }  // End:  for loop over requested unit types
  }    // End:  if PDIS data present

  // write PDOS graph, the phonon density of states
  if(dos_data) {
    for (unsigned int f = 0; f < unitName.size(); ++f) {
      string plot_type = "PDOS";
      string gnuplotscript_name = directory + "/" + aurostd::ProperFileName("GNUPLOT_" + chemName + "_" + unitName.at(f) + "_PDOS.gp");
      ofstream fout;
      fout.open(gnuplotscript_name.c_str());
      if(fout.is_open()) {
        cout << "pflow::PLOT_PHDISP:  Writing GNUPLOT script for phonon density of states to " << gnuplotscript_name << endl;
      } else {
        cout << "pflow::PLOT_PHDISP:  ERROR!  Couldn't open " << gnuplotscript_name.c_str() << " for writing." << endl;
        return FALSE;
      }
      fout << gnuplotPhononHeader(directory, chemName, unitName.at(f), plot_type, graph_title);
      fout << gnuplotPhononDos(phondosdat, "0.08", "0.15");
      fout.close();

      // Call gnuplot to plot
      aurostd::execute(XHOST.command("gnuplot") + " " + gnuplotscript_name);

      // call any conversion
      callOutputConversion(XHOST.vflag_control.flag("PRINT_MODE::PNG"),
                           XHOST.vflag_control.flag("PRINT_MODE::PDF"),
                           XHOST.vflag_control.flag("PRINT_MODE::JPG"),
                           XHOST.vflag_control.flag("PRINT_MODE::GIF"),
                           directory, chemName, unitName.at(f), plot_type);

      // remove gnuplot scripts if they point to temporary files
      cleanGnuplotScript(pdis_was_decompressed, pdos_was_decompressed, thermo_was_decompressed, gnuplotscript_name);
    }  // End:  for loop over requested unit types
  }    // End:  if PDOS data present

  // write THERMO graph, the triple plot of F, S, and C_v
  if(thermo_data) {
    for (unsigned int b = 0; b < unitName.size(); ++b) {
      string plot_type = "THERMO";
      string gnuplotscript_name = directory + "/" + aurostd::ProperFileName("GNUPLOT_" + chemName + "_" + unitName.at(b) + "_" + plot_type + ".gp");

      ofstream fout;
      fout.open(gnuplotscript_name.c_str());
      if(fout.is_open()) {
        cout << "pflow::PLOT_PHDISP:  Writing GNUPLOT script for F, S, and C to " << gnuplotscript_name << endl;
      } else {
        cout << "pflow::PLOT_PHDISP:  ERROR!  Couldn't open " << gnuplotscript_name.c_str() << " for writing." << endl;
        return FALSE;
      }
      fout << gnuplotPhononHeader(directory, chemName, unitName.at(b), plot_type, graph_title);
      fout << gnuplotPhononThermo(phonthrmdat, "0.85");
      fout.close();

      // Call gnuplot to plot the bands structure
      aurostd::execute(XHOST.command("gnuplot") + " " + gnuplotscript_name);

      // call any conversion
      callOutputConversion(XHOST.vflag_control.flag("PRINT_MODE::PNG"),
                           XHOST.vflag_control.flag("PRINT_MODE::PDF"),
                           XHOST.vflag_control.flag("PRINT_MODE::JPG"),
                           XHOST.vflag_control.flag("PRINT_MODE::GIF"),
                           directory, chemName, unitName.at(b), plot_type);

      // remove gnuplot scripts if they point to temporary files
      cleanGnuplotScript(pdis_was_decompressed, pdos_was_decompressed, thermo_was_decompressed, gnuplotscript_name);
    }  // End:  for loop over units
  }    // End:  if thermo

  cout << "Exiting pflow::PLOT_PHDIS() . . . " << endl;
  return TRUE;
}  // End: function pflow::PLOT_PHDISP()
}

//***************************************************************************//
// gnuplotPhononHeader
//***************************************************************************//
string gnuplotPhononHeader(string directory, string name, string unit, string plot_type, string title) {
  // a helper function for pflow::PLOT_PHDISP() that returns a string of gnuplot instructions
  // that define an .eps plot
  string temp;
  temp.append("#Generated by AFLOW (Kevin Rasch [kevin.rasch@duke.edu], 2013, Duke) \n");
  temp.append("set term postscript eps enhanced color font \"Times-Roman, 40\" size 18, 10.125 \n");
  temp.append("set output \"" + directory + "/" + name + "_" + unit + "_" + plot_type + ".eps" + "\" \n");
  if((plot_type == "PHON") || (plot_type == "THERMO") || (plot_type == "PDISDOS")) {
    temp.append("set multiplot title \"" + title + "\"  \n");
  } else if((plot_type == "PDIS") || (plot_type == "PDOS")) {
    temp.append("set title \"" + title + "\"  \n");
  }
  temp.append("set label '" + AFLOWLIB_CONSORTIUM_STRING + "' at screen 0.75, 0.02 font \"Times-Roman, 32\"  \n");
  temp.append("\n");
  temp.append("\n");

  return temp;
}

//***************************************************************************//
// gnuplotPhononDisp
//***************************************************************************//
string gnuplotPhononDisp(string dis_path, string ylabel, double conversion,
                         vector<string> specialLabels, vector<double> specialPositions, int nbranches,
                         string rightmargin, string botmargin) {
  // a helper function to pflow::PLOT_PHDISP() that returns a string of gnuplot
  // instructions for a phonon dispersion plot
  string temp = "";
  temp.append("############### \n");
  temp.append("#               \n");
  temp.append("# PH DISP PLOT  \n");
  temp.append("#               \n");
  temp.append("############### \n");
  temp.append("set lmargin at screen 0.08  \n");
  temp.append("set rmargin at screen " + rightmargin + " \n");
  temp.append("set tmargin at screen 0.94  \n");
  temp.append("set bmargin at screen " + botmargin + " \n");
  temp.append("\n");
  temp.append("unset y2tics  \n");
  temp.append("set   ytics  \n");
  temp.append("unset y2label  \n");
  temp.append("set   ylabel 'Frequency (" + ylabel + ")' offset graph 0.002  \n");
  temp.append("set   xrange[*:*]  \n");
  temp.append("set   grid xtics  \n");
  temp.append("set   xtics (");
  for (unsigned int g = 0; g < specialLabels.size(); ++g) {
    temp.append("\"{" + specialLabels.at(g) + "}\" ");
    string dts_pos;
    ostringstream convert;
    convert << specialPositions.at(g);
    dts_pos = convert.str();
    temp.append(dts_pos);
    if(g < specialLabels.size() - 1) {
      temp.append(",");
    }
  }
  temp.append(")  \n");
  temp.append("set   style line 1 lt 1 lc rgb \"black\" lw 5  \n");
  temp.append("\n");
  temp.append("plot \\");
  temp.append("\n");
  // Write 1 line per phonon branch in the form:
  // first phonon branch is the 3rd column of the file
  //    temp.append("     \"PDIS\" u 2:3 pt 7 ps 1 lc rgb \"black\" notitle, \\" << endl;
  string dts;
  ostringstream convert;
  convert << conversion;
  dts = convert.str();
  for (int h = 3; h <= 2 + nbranches; ++h) {
    string its;
    ostringstream convert2;
    convert2 << h;
    its = convert2.str();
    temp.append("     \"" + dis_path + "\" u 2:($" + its + "*" + dts +
                ") pt 7 ps 1 lc rgb \"black\
          \" notitle, \\");
    temp.append("\n");
  }

  return temp;
}

//***************************************************************************//
// gnuplotPhononDos
//***************************************************************************//
string gnuplotPhononDos(string pdos_path, string leftmargin, string botmargin) {
  // Helper function to pflow::PLOT_PHDISP() that returns a string of the gnuplot
  // instructions for plotting the phonon density of states
  string temp = "";
  temp.append("############### \n");
  temp.append("#               \n");
  temp.append("# PH DOS PLOT   \n");
  temp.append("#               \n");
  temp.append("############### \n");
  temp.append("set lmargin at screen " + leftmargin + " \n");
  temp.append("set rmargin at screen 0.92 \n");
  temp.append("set tmargin at screen 0.94 \n");
  temp.append("set bmargin at screen " + botmargin + " \n");
  temp.append("\n");
  temp.append("unset title \n");
  temp.append("unset xtics \n");
  temp.append("set   yrange[0:*] \n");
  temp.append("unset ytics \n");
  temp.append("set   y2tics \n");
  temp.append("set   y2label 'Energy (meV)' offset graph 0.000 rotate by 270 \n");
  temp.append("\n");
  temp.append("plot \\");
  temp.append("\n");
  temp.append("     \"" + pdos_path + "\" u 4:3 lt rgb \"black\" with filledcurves x1 notitle \n");

  return temp;
}

//***************************************************************************//
// gnuplotPhononThermo
//***************************************************************************//
string gnuplotPhononThermo(string thermo_path, string topmargin) {
  // Helper function to pflow::PLOT_PHDISP() that returns a string of the gnuplot
  // instructions for plotting F_vib, S_vib, and C_v
  string temp = "";
  temp.append("###############  \n");
  temp.append("#  \n");
  temp.append("# F_vib \n");
  temp.append("#   \n");
  temp.append("############### \n");
  temp.append("set lmargin at screen 0.08 \n");
  temp.append("set rmargin at screen 0.32 \n");
  temp.append("set tmargin at screen " + topmargin + "\n");
  temp.append("set bmargin at screen 0.15 \n");
  temp.append("\n");
  temp.append("set   tics font \"Times-Roman, 20\"  \n");
  temp.append("unset title \n");
  temp.append("set   ylabel 'F_{vib} (meV/cell)' font \"Times-Roman, 30\" offset graph 0.1,0.0  \n");
  temp.append("set   xlabel 'Temperature (K)' font \"Times-Roman, 30\" \n");
  temp.append("unset y2label  \n");
  temp.append("set   xtics auto \n");
  temp.append("set   grid xtics \n");
  temp.append("set   ytics mirror \n");
  temp.append("set   xrange[*:*] \n");
  temp.append("set   yrange[*:*] \n");
  temp.append("plot \\");
  temp.append("\n");
  temp.append("     \"" + thermo_path + "\" u 1:4 lw 6 with lines notitle \n");
  temp.append("\n");
  temp.append("\n");
  temp.append("###############  \n");
  temp.append("#                \n");
  temp.append("# S_vib          \n");
  temp.append("#                \n");
  temp.append("###############  \n");
  temp.append("set lmargin at screen 0.38 \n");
  temp.append("set rmargin at screen 0.62 \n");
  temp.append("set tmargin at screen " + topmargin + "\n");
  temp.append("set bmargin at screen 0.15 \n");
  temp.append("\n");
  temp.append("unset title \n");
  temp.append("set   ylabel 'S_{vib} (k_B/cell)' font \"Times-Roman, 30\" offset graph 0.1, 0.0 \n");
  temp.append("set   xtics \n");
  temp.append("set   ytics mirror \n");
  temp.append("\n");
  temp.append("plot \\");
  temp.append("\n");
  temp.append("     \"" + thermo_path + "\" u 1:5 lw 6 with lines notitle \n");
  temp.append("\n");
  temp.append("\n");
  temp.append("############### \n");
  temp.append("#               \n");
  temp.append("# C_v           \n");
  temp.append("#               \n");
  temp.append("############### \n");
  temp.append("set lmargin at screen 0.68 \n");
  temp.append("set rmargin at screen 0.92 \n");
  temp.append("set tmargin at screen " + topmargin + "\n");
  temp.append("set bmargin at screen 0.15 \n");
  temp.append("\n");
  temp.append("unset title \n");
  temp.append("unset ylabel \n");
  temp.append("set   y2label 'C_v (k_B/cell)' font \"Times-Roman, 30\" offset graph -0.05,0.0 rotate by 270 \n");
  temp.append("set   xtics \n");
  temp.append("unset ytics \n");
  temp.append("set   y2tics mirror \n");
  temp.append("\n");
  temp.append("plot \\");
  temp.append("\n");
  temp.append("     \"" + thermo_path + "\" u 1:6 lw 6 with lines notitle \n");
  temp.append("\n");
  temp.append("\n");
  temp.append("unset multiplot \n");

  return temp;
}

//***************************************************************************//
// callOutputConversion
//***************************************************************************//
void callOutputConversion(bool PNG, bool PDF, bool JPG, bool GIF, string directory, string name, string unit, string plot) {
  // Helper function for pflow::PLOT_PHDISP() that calls `convert' to produce the requested image types from .eps
  if(PNG) {
    cout << "pflow::PLOT_PHDISP:  Converting image [" << directory << "/" << name << ".eps] to .png . . . " << endl;
    aurostd::execute(XHOST.command("convert") + " " + directory + "/" + name + "_" + unit + "_" + plot + ".eps " + directory + "/" + name + "_" + unit + "_" + plot + ".png");
  }
  if(PDF) {
    cout << "pflow::PLOT_PHDISP:  Converting image [" << directory << "/" << name << ".eps] to .pdf . . . " << endl;
    aurostd::execute(XHOST.command("convert") + " " + directory + "/" + name + "_" + unit + "_" + plot + ".eps " + directory + "/" + name + "_" + unit + "_" + plot + ".pdf");
  }
  if(JPG) {
    cout << "pflow::PLOT_PHDISP:  Converting image [" << directory << "/" << name << ".eps] to .jpg . . . " << endl;
    aurostd::execute(XHOST.command("convert") + " " + directory + "/" + name + "_" + unit + "_" + plot + ".eps " + directory + "/" + name + "_" + unit + "_" + plot + ".jpg");
  }
  if(GIF) {
    cout << "pflow::PLOT_PHDISP:  Converting image [" << directory << "/" << name << ".eps] to .gif . . . " << endl;
    aurostd::execute(XHOST.command("convert") + " " + directory + "/" + name + "_" + unit + "_" + plot + ".eps " + directory + "/" + name + "_" + unit + "_" + plot + ".gif");
  }
}

//***************************************************************************//
// cleanGnuplotScript
//***************************************************************************//
void cleanGnuplotScript(bool file1, bool file2, bool file3, string script) {
  // Helper function for pflow::PLOT_PHDISP() that removes the gnuplot scripts
  if((file1 || file2) || file3) {
    cout << "pflow::PLOT_PHDISP:  Cleaning up . . . Gnuplot script refers to temporary files; removing it." << endl;
    aurostd::execute("rm -f " + script);  // script refers to the above temporary files--saving it is useless
  }
}

// ***************************************************************************
// *                                                                         *
// *              AFlow KEVIN RASCH - Duke University 2013                   *
// *                                                                         *
// ***************************************************************************
