// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *                                                                         *
// ***************************************************************************
// Stefano Curtarolo 
// Dane Morgan (up to 2003)
// Wahyu Setyawan (up to 2009)

#ifndef _AFLOW_PFLOW_MAIN_CPP_
#define _AFLOW_PFLOW_MAIN_CPP_

#include<sys/stat.h>
#include "aflow.h"
#include "aflow_pflow.h"
// [OBSOLETE] #include "aflow_contrib_wahyu.h"
// [OBSOLETE] #include "aflow_contrib_roman_funcs.h"
#include "aflow_contrib_shidong_main.h"
#include "aflow_symmetry_spacegroup.h"
// [OBSOLETE] [JUNKAI] #include "aflow_contrib_junkai_basic.h"
#include "aflow_contrib_kesong.h"
#include "aflow_bader.h"
#include "aflow_chull.h"
// [OBSOLETE] #include "aflow_contrib_cormac.h"
#include "aflowlib.h"

extern double NearestNeighbour(const xstructure& a);

// ***************************************************************************
using aurostd::FileExist;

uint PflowARGs(vector<string> &argv,vector<string> &cmds,aurostd::xoption &vpflow) {  // vpflow is really XHOST.vflag_pflow
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  // GENERAL STUFF
   
  vpflow.flag("PFLOW_HELP",aurostd::args2flag(argv,cmds,"--HELP|--help"));
  vpflow.flag("PROTOS",(aurostd::args2flag(argv,cmds,"--protos|--prototypes") || aurostd::args2flag(argv,cmds,"--proto")));// && (argv.size()==2));
  // [OBSOLETE]  vpflow.flag("PROTOS_ICSD",(aurostd::args2flag(argv,cmds,"--protos_icsd|--prototypes_icsd") || aurostd::args2flag(argv,cmds,"--proto_icsd")) && (argv.size()==2));
  
  vpflow.args2addattachedscheme(argv,cmds,"PROTOS_ICSD","--protos_icsd=|--prototypes_icsd=","");
  
  if(!vpflow.flag("PROTOS_ICSD")) vpflow.flag("PROTOS_ICSD",aurostd::args2flag(argv,cmds,"--protos_icsd|--prototypes_icsd"));
  
  //  bool XXX=aurostd::args2flag(argv,cmds,"--xxx") && argv.at(1)=="--xxx";
  //  bool YYY=aurostd::args2flag(argv,cmds,"--yyy");

  vpflow.flag("ABCCAR",aurostd::args2flag(argv,cmds,"--abccar"));
  vpflow.flag("ACE",aurostd::args2flag(argv,cmds,"--ace"));
  //DX 8/18/17 [OBSOLETE] vpflow.flag("AGROUP",aurostd::args2flag(argv,cmds,"--sitepointgroup|--agroup"));
  //DX 8/18/17 - Added tolerance and no_scan options to Xgroups - START
  vpflow.args2addattachedscheme(argv,cmds,"AGROUP","--sitepointgroup=|--agroup=","");
  if(vpflow.flag("AGROUP")){
    vpflow.flag("SYMMETRY::NO_SCAN",aurostd::args2flag(argv,cmds,"--no_scan"));
    if(aurostd::args2attachedflag(argv,"--sitepointgroup=|--agroup=")){ //DX 8/3/17
      vpflow.args2addattachedscheme(argv,cmds,"SYMMETRY::TOLERANCE","--sitepointgroup=|--agroup=","1");
    }
    vpflow.flag("SYMMETRY::SCREEN_ONLY",aurostd::args2flag(argv,cmds,"--screen_only")); //DX 8/3/17
    //DX 9/21/17 - MAGNETIC SYMMETRY - START
    vpflow.args2addattachedscheme(argv,cmds,"SYMMETRY::MAGNETIC","--mag=|--magnetic=|--magmom=","1"); //DX 8/3/17
    //DX 9/21/17 - MAGNETIC SYMMETRY - END
  }
  //DX 8/18/17 - Added tolerance and no_scan options to Xgroups - END
  vpflow.flag("AGROUP2",aurostd::args2flag(argv,cmds,"--sitepointgroup2|--agroup2"));
  vpflow.flag("AGROUP2m",aurostd::args2flag(argv,cmds,"--sitepointgroup2m|--agroup2m"));
  vpflow.flag("AFLOWIN",aurostd::args2flag(argv,cmds,"--aflowin"));
  vpflow.flag("ALPHABETIC",aurostd::args2flag(argv,cmds,"--alpha|--alphabetic"));
  // [OBSOLETE] vpflow.flag("ALPHA_COMPOUND",aurostd::args2flag(argv,cmds,"--alpha_compound|--alpha_compounds"));
  vpflow.args2addattachedscheme(argv,cmds,"ALPHA_COMPOUND","--alpha_compound=|--alpha_compounds=","");
  // [OBSOLETE] vpflow.flag("ALPHA_SPECIES",aurostd::args2flag(argv,cmds,"--alpha_species|--alpha_specie"));
  vpflow.args2addattachedscheme(argv,cmds,"ALPHA_SPECIES","--alpha_species=|--alpha_specie=","");
  // [OBSOLETE] vpflow.flag("ANGLES",aurostd::args2flag(argv,cmds,"--angle|--angles"));
  vpflow.args2addattachedscheme(argv,cmds,"ANGLES","--angle=|--angles=","0.0");

  vpflow.args2addattachedscheme(argv,cmds,"AFLOWLIB","--aflowlib=","");
  vpflow.args2addattachedscheme(argv,cmds,"AFLOWLIB_AUID2AURL","--aflowlib_auid2aurl=|--auid2aurl=","");
  vpflow.args2addattachedscheme(argv,cmds,"AFLOWLIB_AURL2AUID","--aflowlib_aurl2auid=|--aurl2auid=","");
  vpflow.args2addattachedscheme(argv,cmds,"AFLOWLIB_AUID2LOOP","--aflowlib_auid2loop=|--auid2loop=","");
  vpflow.args2addattachedscheme(argv,cmds,"AFLOWLIB_AURL2LOOP","--aflowlib_aurl2loop=|--aurl2loop=","");

  // Commands for serializing bands and DOS data to JSON
  vpflow.args2addattachedscheme(argv,cmds,"DOSDATA2JSON","--dosdata2json=","./"); // Eric G
  if(vpflow.flag("DOSDATA2JSON")){vpflow.args2addattachedscheme(argv,cmds,"DOSDATA2JSON::PARAMS","--dos_parameters=","");}  //CO 180214 removed params confusion
  vpflow.args2addattachedscheme(argv,cmds,"BANDSDATA2JSON","--bandsdata2json=","./"); // Eric G
  // End commands


  // bool APENNSY=aurostd::args2flag(argv,cmds,"--apennsy|-apennsy|--pennsy|-pennsy");
 
  //corey
  vpflow.flag("BADER",aurostd::args2flag(argv,cmds,"--bader"));
  if(vpflow.flag("BADER")) {    // corey
    //manage directory, ASSUME LOCAL IF NOT SPECIFIED
    if(XHOST.vflag_control.flag("DIRECTORY")) {
      vpflow.addattachedscheme("BADER::DIRECTORY",XHOST.vflag_control.getattachedscheme("DIRECTORY"),TRUE);
    } else {
      vpflow.addattachedscheme("BADER::DIRECTORY",".",TRUE);
    }
    XHOST.vflag_control.getattachedscheme("DIRECTORY");
    vpflow.flag("BADER::USAGE",aurostd::args2flag(argv,cmds,"--usage"));  // usage
    vpflow.flag("BADER::CRITICAL_POINTS",aurostd::args2flag(argv,cmds,"--critical_points|--cp"));  // critical points
    vpflow.args2addattachedscheme(argv,cmds,"BADER::CALCULATE","--calculate=|--calc=","");    // -c
    vpflow.args2addattachedscheme(argv,cmds,"BADER::NOCALCULATE","--nocalculate=|--nocalc=","");  // -n
    vpflow.args2addattachedscheme(argv,cmds,"BADER::PARTITION","--partition=|--part=","");  // -b
    vpflow.args2addattachedscheme(argv,cmds,"BADER::REFINE_EDGE_METHOD","--refine_edge_method=|--rem=|--r=","");  // -r
    vpflow.args2addattachedscheme(argv,cmds,"BADER::REFERENCE","--reference=|--ref=","");     // -ref
    vpflow.args2addattachedscheme(argv,cmds,"BADER::VACUUM","--vacuum=|--vac=","");   // -vac
    vpflow.args2addattachedscheme(argv,cmds,"BADER::TERMINATE","--terminate=|--term=","");    // -m
    vpflow.args2addattachedscheme(argv,cmds,"BADER::PRINT_ALL","--print_all=","");    // -p all_atom | -p all_bader
    vpflow.args2addattachedscheme(argv,cmds,"BADER::PRINT_INDEX","--print_index=|--print_idx=","");     // -p sel_atom | -p sel_bader
    vpflow.args2addattachedscheme(argv,cmds,"BADER::PRINT_SELECT_ATOM","--print_select_atom=|--print_sel_atom=","");    // -p sel_atom
    vpflow.args2addattachedscheme(argv,cmds,"BADER::PRINT_SELECT_BADER","--print_select_bader=|--print_sel_bader=",""); // -p sel_bader
    vpflow.args2addattachedscheme(argv,cmds,"BADER::PRINT_SUM_ATOM","--print_sum_atom=","");  // -p sum_atom
    vpflow.args2addattachedscheme(argv,cmds,"BADER::PRINT_SUM_BADER","--print_sum_bader=","");    // -p sum_bader
    vpflow.flag("BADER::QUIET",aurostd::args2flag(argv,cmds,"--quiet|--q"));  // quiets output
    vpflow.flag("BADER::CONSOLIDATE_ATOMS2SPECIES",aurostd::args2flag(argv,cmds,"--consolidate_atoms2species|--a2s"));  // quiets output
    if(vpflow.flag("BADER::CONSOLIDATE_ATOMS2SPECIES")) {
      //we can consolidate and delete summing chgcars
      vpflow.flag("BADER::REMOVE_BADER_ATOMS",aurostd::args2flag(argv,cmds,"--remove_bader_atoms|--rba"));  // keep only consolidated files
    }
    vpflow.args2addattachedscheme(argv,cmds,"BADER::JVXL_ALL_SPECIES","--jvxl_all_species=|--jvxl=",""); // allows jvxl output for atoms2species
    if(vpflow.flag("BADER::JVXL_ALL_SPECIES")) {
      //when jvxl-ing, we can just remove summing chgcars and keep full consolidated chgcar, or delete all and just keep jvxl
      vpflow.flag("BADER::REMOVE_BADER_ATOMS",aurostd::args2flag(argv,cmds,"--remove_bader_atoms|--rba"));  // keep only consolidated files
      vpflow.flag("BADER::KEEP::JVXL_ONLY",aurostd::args2flag(argv,cmds,"--keep=jvxl_only|--keep_jvxl_only|--jvxl_only"));  // removes consolidated chgcar_files, keeps only jvxl
    }
  }
  
  // [OBSOLETE] vpflow.flag("BANDS",aurostd::args2flag(argv,cmds,"--bands") && argv.at(1)=="--bands");
  vpflow.args2addattachedscheme(argv,cmds,"BANDS","--bands=","");
  
  vpflow.flag("BANDSTRUCTURE",aurostd::args2flag(argv,cmds,"--bandstructure|--bandsstructures|--bands_structures|--band_structures|--bands_structure|--band_structure|--bs"));
  vpflow.flag("BZPLOT",aurostd::args2flag(argv,cmds,"--plotbz|--bzplot"));
  // [OBSOLETE] vpflow.flag("BZPLOTUSEKPOINTS",aurostd::args2flag(argv,cmds,"--bzplotuseKPOINTS|--bzplotdusekpoints"));
  vpflow.args2addattachedscheme(argv,cmds,"BZPLOTUSEKPOINTS","--bzplotuseKPOINTS=|--bzplotdusekpoints=","");   
  vpflow.flag("BZPLOTDATA",aurostd::args2flag(argv,cmds,"--bzplotdata"));
  // [OBSOLETE] vpflow.flag("BZPLOTDATAUSEKPOINTS",aurostd::args2flag(argv,cmds,"--bzplotdatauseKPOINTS|--bzplotdatausekpoints"));
  vpflow.args2addattachedscheme(argv,cmds,"BZPLOTDATAUSEKPOINTS","--bzplotdatauseKPOINTS=|--bzplotdatausekpoints=","");

  // [OBSOLETE]  vpflow.flag("FIX_BANDS",aurostd::args2flag(argv,cmds,"--fix_bands"));// && argv.size()==9);
  vpflow.args2addattachedscheme(argv,cmds,"FIX_BANDS","--fix_bands=","");

  // DX AND COREY - START

  vpflow.args2addattachedscheme(argv,cmds,"FULLSYMMETRY","--aflow-sym=|--AFLOW-SYM=|--AFLOWSYM=|--aflowSYM=|--aflowsym=|--full_symmetry=|--full_sym=|--fullsym=",""); //DX 8/3/17 Added other aliases
  if(vpflow.flag("FULLSYMMETRY")){
    vpflow.flag("FULLSYMMETRY::NO_SCAN",aurostd::args2flag(argv,cmds,"--no_scan"));
    if(aurostd::args2attachedflag(argv,"--aflow-sym=|--AFLOW-SYM=|--AFLOWSYM=|--aflowSYM=|--aflowsym=|--full_symmetry=|--full_sym=|--fullsym=")){ //DX 8/3/17
      vpflow.args2addattachedscheme(argv,cmds,"FULLSYMMETRY::TOLERANCE","--aflow-sym=|--AFLOW-SYM=|--AFLOWSYM=|--aflowSYM=|--aflowsym=|--full_symmetry=|--full_sym=|--fullsym=","1"); //DX 8/3/17
    }
    vpflow.flag("FULLSYMMETRY::SCREEN_ONLY",aurostd::args2flag(argv,cmds,"--screen_only")); //DX 8/3/17
    //DX 9/21/17 - MAGNETIC SYMMETRY - START
    vpflow.args2addattachedscheme(argv,cmds,"FULLSYMMETRY::MAGNETIC","--mag=|--magnetic=|--magmom=","1"); //DX 8/3/17
    //DX 9/21/17 - MAGNETIC SYMMETRY - END
  }
  // DX AND COREY - END

  vpflow.flag("BANDGAP_WAHYU",aurostd::args2flag(argv,cmds,"--bandgap_wahyu|--gap_wahyu"));
  // [OBSOLETE] vpflow.flag("BANDGAP",aurostd::args2flag(argv,cmds,"--bandgap" ) ); // CAMILO
  vpflow.args2addattachedscheme(argv,cmds,"BANDGAP","--bandgap=","./");
  
  vpflow.flag("BANDGAPS",aurostd::args2flag(argv,cmds,"--bandgaps|--gaps"));
  vpflow.flag("BANDGAPDOS",aurostd::args2flag(argv,cmds,"--bandgap_from_dos|--bandgap_from_DOS"));
  vpflow.flag("BANDGAPLISTDOS",aurostd::args2flag(argv,cmds,"--bandgaplist_from_dos|--bandgaplist_from_DOS"));
  // [OBSOLETE] vpflow.flag("BZDIRECTION",aurostd::args2flag(argv,cmds,"--bzdirection|--bzdirections|--bzd"));
  vpflow.args2addattachedscheme(argv,cmds,"BZDIRECTION","--bzdirection=|--bzdirections=|--bzd=","");
  vpflow.flag("BZMAX",aurostd::args2flag(argv,cmds,"--BZmax"));

  // [OBSOLETE] vpflow.flag("CAGES",aurostd::args2flag(argv,cmds,"--cages"));
  vpflow.args2addattachedscheme(argv,cmds,"CAGES","--cages=","-1.0");
  // cerr << "vpflow.flag(\"CAGES\")=" << vpflow.flag("CAGES") << endl;
  // cerr << "vpflow.getattachedscheme(\"CAGES\")=" << vpflow.getattachedscheme("CAGES") << endl;
  // exit(0);
  
  vpflow.flag("CALCULATED_ICSD_RANDOM",(aurostd::args2flag(argv,cmds,"--calculated=icsd")) && aurostd::args2flag(argv,cmds,"--random|--rnd"));
  // [OBSOLETE] string calculated=aurostd::args2attachedstring(argv,"--calculated=","");
  // [OBSOLETE] vpflow.flag("CALCULATED_ALL",aurostd::args2flag(argv,cmds,"--calculated|--calcs"));
  // [OBSOLETE] vpflow.flag("CALCULATED_ICSD",aurostd::substring2bool(calculated,"icsd"));
  // [OBSOLETE] vpflow.flag("CALCULATED_LIB1",aurostd::substring2bool(calculated,"lib1"));
  // [OBSOLETE] vpflow.flag("CALCULATED_LIB2",aurostd::substring2bool(calculated,"lib2"));
  // [OBSOLETE] vpflow.flag("CALCULATED_LIB3",aurostd::substring2bool(calculated,"lib3"));
  // [OBSOLETE] vpflow.flag("CALCULATED_LIB4",aurostd::substring2bool(calculated,"lib4"));
  // [OBSOLETE] vpflow.flag("CALCULATED_LIB5",aurostd::substring2bool(calculated,"lib5"));
  // [OBSOLETE] vpflow.flag("CALCULATED_LIB6",aurostd::substring2bool(calculated,"lib6"));
  // [OBSOLETE] vpflow.flag("CALCULATED_LIB7",aurostd::substring2bool(calculated,"lib7"));
  // [OBSOLETE] vpflow.flag("CALCULATED_LIB8",aurostd::substring2bool(calculated,"lib8"));
  // [OBSOLETE] vpflow.flag("CALCULATED_LIB9",aurostd::substring2bool(calculated,"lib9"));
  vpflow.args2addattachedscheme(argv,cmds,"CALCULATED","--calculated=","all");
  // cerr << "vpflow.flag(\"CALCULATED\")=" << vpflow.flag("CALCULATED") << endl;
  // cerr << vpflow.getattachedscheme("CALCULATED") << endl;

  vpflow.flag("CART",aurostd::args2flag(argv,cmds,"--cart|-cart|-c|--cartesian"));
  vpflow.flag("CHECKINTEGRITIY",aurostd::args2flag(argv, cmds,"--check_integrity|--checki"));
  vpflow.flag("CIF",aurostd::args2flag(argv,cmds,"--cif"));

  // [OBSOLETE] vpflow.flag("CHANGESUFFIX",aurostd::args2attachedflag(argv,cmds,"--suffix=")); //KESONG DEC. 22ND, 2013
  // [OBSOLETE] vpflow.addattachedscheme"CHANGESUFFIX",aurostd::args2attachedstring(argv,"--suffix=",""),vpflow.flag("CHANGESUFFIX"));
  vpflow.args2addattachedscheme(argv,cmds,"CHANGESUFFIX","--suffix=","./");

  //corey
  vpflow.args2addattachedscheme(argv,cmds,"CHGCAR2JVXL","--chgcar2jvxl=|--c2j=","");    // create JVXL from CHGCAR, corey
  if(vpflow.flag("CHGCAR2JVXL")) {
    vpflow.flag("CHGCAR2JVXL::USAGE",aurostd::args2flag(argv,cmds,"--usage"));  // usage
    vpflow.args2addattachedscheme(argv,cmds,"CHGCAR2JVXL::OUTPUT","--output=|--o=","");
  }

  //corey
  vpflow.args2addattachedscheme(argv,cmds,"CHGDIFF","--chgdiff=","");
  if(vpflow.flag("CHGDIFF")) {    // corey
    vpflow.flag("CHGDIFF::USAGE",aurostd::args2flag(argv,cmds,"--usage"));  // usage
    vpflow.args2addattachedscheme(argv,cmds,"CHGDIFF:OUTPUT","--output=|--o=","");
  }

  vpflow.flag("CHGINT",aurostd::args2flag(argv,cmds,"--chgint") && argv.at(1)=="--chgint");
  //corey
  vpflow.args2addattachedscheme(argv,cmds,"CHGSUM","--chgsum=","");
  if(vpflow.flag("CHGSUM")) {    // corey
    vpflow.flag("CHGSUM::USAGE",aurostd::args2flag(argv,cmds,"--usage"));  // usage
    vpflow.args2addattachedscheme(argv,cmds,"CHGSUM::OUTPUT","--output=|--o=","");
  }
  
  vpflow.flag("CLEANALL",aurostd::args2flag(argv,cmds,"--cleanall|--clean_all"));

  // [OBSOLETE] vpflow.flag("COMPARE",aurostd::args2flag(argv,cmds,"--compare"));
  vpflow.args2addattachedscheme(argv,cmds,"COMPARE","--compare=","");
  vpflow.flag("CMPSTR",aurostd::args2flag(argv,cmds,"--cmp_str") && argv.at(1)=="--cmp_str");

  vpflow.flag("CHULL::INIT",aurostd::args2flag(argv,cmds,"--convex_hull|--chull"));  // initiate chull calculation
  //vpflow.args2addattachedscheme(argv,cmds,"CHULL::INIT","--convex_hull=|--chull=","");  //initiate chull calculation
  vpflow.args2addattachedscheme(argv,cmds,"PFLOW::ALLOY","--alloy=","");  //define alloy
  vpflow.flag("PFLOW::LOAD_ENTRIES_ENTRY_OUTPUT",aurostd::args2flag(argv,cmds,"--load_entries_entry_output|--loadentriesentryoutput|--leo"));    //verbose loading entries output
  vpflow.flag("PFLOW::LOAD_API",aurostd::args2flag(argv,cmds,"--load_API|--load_api|--loadapi|--lapi|--api"));  //force load api
  if(vpflow.flag("CHULL::INIT")) {
    //do NOT rearrange flags, order is important!
    vpflow.flag("CHULL::USAGE",aurostd::args2flag(argv,cmds,"--usage"));  // usage
    vpflow.flag("CHULL::LATEX_OUTPUT",aurostd::args2flag(argv,cmds,"--latex_output|--latexoutput"));    //verbose latex output (cout and FileMESSAGE)
    vpflow.flag("CHULL::LATEX_INTERACTIVE",aurostd::args2flag(argv,cmds,"--latex_interactive|--latexinteractive"));  //interact with latex (execute vs. execute2string)
    if(vpflow.flag("CHULL::LATEX_INTERACTIVE")&&!vpflow.flag("CHULL::LATEX_OUTPUT")){
      vpflow.flag("CHULL::LATEX_OUTPUT",TRUE);    //keep verbose latex output
    }
    //vpflow.flag("PFLOW::LOAD_ENTRIES_ONLY_ALPHABETICAL",TRUE);  //for chull, this is always true (non-alphabetic are crap entries), set inside
    if(XHOST.QUIET&&vpflow.flag("PFLOW::LOAD_ENTRIES_ENTRY_OUTPUT")){
      vpflow.flag("PFLOW::LOAD_ENTRIES_ENTRY_OUTPUT",FALSE);    //keep quiet
    }
    vpflow.args2addattachedscheme(argv,cmds,"PFLOW::LOAD_LIBRARY","--load_library=|--loadlibrary=|--ll=","");  //get libraries to load
    vpflow.flag("CHULL::SEE_NEGLECT",aurostd::args2flag(argv,cmds,"--see_neglect|--seeneglect|--sn")); //see why compounds get neglected
    //deal with picking libraries later
    vpflow.args2addattachedscheme(argv,cmds,"CHULL::NEGLECT","--neglect=|--ban=",""); //remove by auid from chull calculation
    vpflow.args2addattachedscheme(argv,cmds,"CHULL::REMOVE_EXTREMA","--remove_extreme_points=|--removeextremepoints=|--remove_extrema=|--removeextrema=|--rep=","");   //set threshold, from bottom for enthalpy of formation, from top for entropic temperature, units of meV OR K
    vpflow.flag("CHULL::LOG",aurostd::args2flag(argv,cmds,"--keep=log|--keep_log|--keeplog|--log"));    //ONLY LDAU CALCULATIONS - ADD FUNCTIONALITY
    vpflow.args2addattachedscheme(argv,cmds,"CHULL::STABILITY_CRITERION","--stability_criterion=|--stabilitycriterion=|--stable_criterion=|--scriterion=|--sc=",""); //calculate stable criterion for point
    vpflow.args2addattachedscheme(argv,cmds,"CHULL::DISTANCE_TO_HULL","--distance_to_hull=|--dist2hull=",""); //calculate stable criterion for point
    vpflow.flag("CHULL::ENTROPIC_TEMPERATURE",aurostd::args2flag(argv,cmds,"--entropic_temperature|--entropictemperature|--entroptemp"));    //use entropic temperature instead of enthalpy of formation
    vpflow.flag("CHULL::SKIP_STRUCTURE_COMPARISON",aurostd::args2flag(argv,cmds,"--skip_structure_comparison|--skipstructruecomparison|--skipstructcomp|--ssc"));    //use entropic temperature instead of enthalpy of formation
    vpflow.flag("CHULL::INCLUDE_UNRELIABLE_HULLS",aurostd::args2flag(argv,cmds,"--include_unreliable_hulls|--include_unreliable|--iuh"));    //use entropic temperature instead of enthalpy of formation
    vpflow.flag("CHULL::INCLUDE_OUTLIERS",aurostd::args2flag(argv,cmds,"--include_outliers|--io"));    //use entropic temperature instead of enthalpy of formation
    vpflow.args2addattachedscheme(argv,cmds,"CHULL::OUTPUT","--output=|--o=|--print=|--p=",""); //determine how to get output
    vpflow.flag("CHULL::SCREEN_ONLY",aurostd::args2flag(argv,cmds,"--screen_only"));    //print to screen
    if(vpflow.flag("CHULL::SCREEN_ONLY")){XHOST.QUIET=true;}
    vpflow.flag("CHULL::GNUPLOT_DOC",FALSE);    //depreciated, NO gnuplot options, all LATEX
    //need to determine output type immediately, so we can accept/ignore output specific flags (SAFE)
    if(vpflow.flag("CHULL::OUTPUT")){
      vector<string> out_forms;
      aurostd::string2tokens(vpflow.getattachedscheme("CHULL::OUTPUT"),out_forms,",");
      for(uint i=0;i<out_forms.size();i++){
        if(out_forms.at(i).at(0)=='A'||out_forms.at(i).at(0)=='a'){
          vpflow.flag("CHULL::APOOL_OUT",TRUE);    //turn on
        }else if(out_forms.at(i).at(0)=='F'||out_forms.at(i).at(0)=='f'){
          vpflow.flag("CHULL::TEXT_DOC",TRUE);    //turn on
          vpflow.flag("CHULL::JSON_DOC",TRUE);    //turn on
          vpflow.flag("CHULL::LATEX_DOC",TRUE);   //turn on default
        }else if(out_forms.at(i).at(0)=='T'||out_forms.at(i).at(0)=='t'){
          vpflow.flag("CHULL::TEXT_DOC",TRUE);    //turn on
        }else if(out_forms.at(i).at(0)=='J'||out_forms.at(i).at(0)=='j'){
          vpflow.flag("CHULL::JSON_DOC",TRUE);    //turn on
        }else if(out_forms.at(i).at(0)=='W'||out_forms.at(i).at(0)=='w'){
          vpflow.flag("CHULL::WEB_DOC",TRUE);     //turn on
          vpflow.flag("CHULL::SKIP_STRUCTURE_COMPARISON",TRUE); //the app cannot handle more than one g-state in the visualization
          vpflow.flag("CHULL::INCLUDE_UNRELIABLE",TRUE); //we include colors on the website
        }else if(out_forms.at(i).at(0)=='L'||out_forms.at(i).at(0)=='l'||out_forms.at(i).at(0)=='P'||out_forms.at(i).at(0)=='p'){    //Latex or Pdf
          vpflow.flag("CHULL::LATEX_DOC",TRUE);   //turn on default
        }//else{ //deal with later
      }
    }else{
      vpflow.flag("CHULL::LATEX_DOC",TRUE);   //default
    }
    if(vpflow.flag("CHULL::STABILITY_CRITERION")||vpflow.flag("CHULL::DISTANCE_TO_HULL")){
      //vpflow.flag("CHULL::TEXT_DOC",FALSE);   //turn off  //leave on, as user might request json/text format output
      //vpflow.flag("CHULL::JSON_DOC",FALSE);   //turn off  //leave on, as user might request json/text format output
      vpflow.flag("CHULL::WEB_DOC",FALSE);    //turn off
      vpflow.flag("CHULL::LATEX_DOC",FALSE);  //turn off
    }
    vpflow.args2addattachedscheme(argv,cmds,"CHULL::PATH","--destination=|--path=",""); //determine how to get output
    if(vpflow.flag("CHULL::LATEX_DOC")) {   //latex specific options
      vpflow.flag("CHULL::DOC_ONLY",aurostd::args2flag(argv,cmds,"--document_only|--documentonly|--doc_only|--doconly|--doc"));   //no convex hull picture
      vpflow.flag("CHULL::NO_DOC",aurostd::args2flag(argv,cmds,"--no_document|--nodocument|--no_doc|--nodoc|--full_page_image|--fullpageimage"));   //no convex hull picture
      if(vpflow.flag("CHULL::NO_DOC")&&vpflow.flag("CHULL::DOC_ONLY")){vpflow.flag("CHULL::DOC_ONLY",FALSE);}
      vpflow.flag("CHULL::KEEP_TEX",aurostd::args2flag(argv,cmds,"--keep=tex|--keep_tex|--keeptex|--tex"));    //keep .tex file to edit
      vpflow.flag("CHULL::IMAGE_ONLY",aurostd::args2flag(argv,cmds,"--image_only|--imageonly|--image"));    //image only, no doc or links, special compilation
      if(vpflow.flag("CHULL::IMAGE_ONLY")&&vpflow.flag("CHULL::DOC_ONLY")) {vpflow.flag("CHULL::IMAGE_ONLY",FALSE);}    //doc takes precedence
      vpflow.flag("CHULL::LIGHT_CONTRAST",aurostd::args2flag(argv,cmds,"--light_contrast|--lightcontrast|--lc")); //lighter blue
      vpflow.flag("CHULL::LARGE_FONT",aurostd::args2flag(argv,cmds,"--large_font|--largefont|--large|--lf"));    //keeps hyperlinks
    }
  }
  
  // [OBSOLETE] vpflow.flag("CLAT",aurostd::args2flag(argv,cmds,"--clat"));
  vpflow.args2addattachedscheme(argv,cmds,"CLAT","--clat=","");
  // [OBSOLETE] vpflow.flag("COMPARE",aurostd::args2flag(argv,cmds,"--compare"));
  vpflow.args2addattachedscheme(argv,cmds,"COMPARE","--compare=","");

  //DAVID
  vpflow.flag("COMPARE_MATERIAL",aurostd::args2attachedflag(argv,cmds,"--compare_material="));
  vpflow.flag("COMPARE_STRUCTURE",aurostd::args2attachedflag(argv,cmds,"--compare_structure="));
  if(vpflow.flag("COMPARE_MATERIAL") || vpflow.flag("COMPARE_STRUCTURE")) {
    vector<string> vinput;
    if(vpflow.flag("COMPARE_MATERIAL")) {
      aurostd::getproto_itemized_vector_string_from_input(argv,"--compare_material",vinput,",");
    }
    else if(vpflow.flag("COMPARE_STRUCTURE")) {
      aurostd::getproto_itemized_vector_string_from_input(argv,"--compare_structure",vinput,",");
    }
    for(uint i=0;i<vinput.size();i++) {
      string scheme_name="COMPARE_STRUCTURE::STRUCTURES_"+aurostd::utype2string(i+1);
      vpflow.addattachedscheme(scheme_name,vinput[i],TRUE);
      vpflow.flag(scheme_name,TRUE);
    }
    vpflow.args2addattachedscheme(argv,cmds,"COMPARE_STRUCTURE::NP","--np=|--num_proc=","");
    vpflow.flag("COMPARE_STRUCTURE::PRINT",aurostd::args2flag(argv,cmds,"--print"));
    vpflow.flag("COMPARE_STRUCTURE::FAST",aurostd::args2flag(argv,cmds,"--fast"));
    vpflow.flag("COMPARE_STRUCTURE::NO_SCALE_VOLUME",aurostd::args2flag(argv,cmds,"--no_scale_volume"));
  }
  vpflow.flag("COMPARE_MATERIAL_DIRECTORY",aurostd::args2flag(argv,cmds,"--compare_material_directory|--compare_material_dir"));
  vpflow.flag("COMPARE_STRUCTURE_DIRECTORY",aurostd::args2flag(argv,cmds,"--compare_structure_directory|--compare_structure_dir"));
  if(vpflow.flag("COMPARE_MATERIAL_DIRECTORY") || vpflow.flag("COMPARE_STRUCTURE_DIRECTORY")) {
    if(XHOST.vflag_control.flag("DIRECTORY")) {
      vpflow.addattachedscheme("COMPARE_STRUCTURE::DIRECTORY",XHOST.vflag_control.getattachedscheme("DIRECTORY"),TRUE);
    } else {
      vpflow.addattachedscheme("COMPARE_STRUCTURE::DIRECTORY",".",true);
    }
    XHOST.vflag_control.getattachedscheme("DIRECTORY");
    vpflow.args2addattachedscheme(argv,cmds,"COMPARE_STRUCTURE::NP","--np=|--num_proc=","");
    vpflow.flag("COMPARE_STRUCTURE::FAST",aurostd::args2flag(argv,cmds,"--fast"));
    vpflow.flag("COMPARE_STRUCTURE::NO_SCALE_VOLUME",aurostd::args2flag(argv,cmds,"--no_scale_volume"));
  }
  //DAVID

  vpflow.flag("CORNERS",aurostd::args2flag(argv,cmds,"--corner|--corners"));

  //DX 9/1/17 [OBSOLETE] vpflow.flag("DATA",aurostd::args2flag(argv,cmds,"--data"));
  vpflow.args2addattachedscheme(argv,cmds,"DATA","--data=",""); //DX 9/1/17 - SGDATA + JSON
  if(vpflow.flag("DATA")){ //DX 9/1/17 - SGDATA + JSON
    vpflow.flag("DATA::NO_SCAN",aurostd::args2flag(argv,cmds,"--no_scan")); //DX 9/1/17 - SGDATA + JSON
    if(aurostd::args2attachedflag(argv,"--data=")){ //DX 9/1/17 - SGDATA + JSON
      vpflow.args2addattachedscheme(argv,cmds,"DATA::TOLERANCE","--data=","1"); //DX 9/1/17 - SGDATA + JSON
    } //DX 9/1/17 - SGDATA + JSON
  } //DX 9/1/17 - SGDATA + JSON
  // [OBSOLETE] vpflow.flag("DATA1",aurostd::args2flag(argv,cmds,"--data1") && argv.at(1)=="--data1");
  vpflow.args2addattachedscheme(argv,cmds,"DATA1","--data1=","");
  vpflow.flag("DATA2",aurostd::args2flag(argv,cmds,"--data2"));
  // [OBSOLETE] vpflow.flag("DEBYE",aurostd::args2flag(argv,cmds,"--debye")); 	
  vpflow.args2addattachedscheme(argv,cmds,"DEBYE","--debye=","");
  // [OBSOLETE] vpflow.flag("DIFF",aurostd::args2attachedflag(argv,cmds,"--diff="));
  vpflow.args2addattachedscheme(argv,cmds,"DIFF","--diff=","./");

  // [OBSOLETE] vpflow.flag("DISP",aurostd::args2flag(argv,cmds,"--disp"));
  vpflow.args2addattachedscheme(argv,cmds,"DISP","--disp=","");
  // [OBSOLETE] vpflow.flag("DIST",aurostd::args2flag(argv,cmds,"--dist"));
  vpflow.args2addattachedscheme(argv,cmds,"DIST","--dist=","");
 
  //DX 9/1/17 [OBSOLETE] vpflow.flag("EDATA",aurostd::args2flag(argv,cmds,"--edata"));
  vpflow.args2addattachedscheme(argv,cmds,"EDATA","--edata=",""); //DX 9/1/17 - SGDATA + JSON
  if(vpflow.flag("EDATA")){ //DX 9/1/17 - SGDATA + JSON
    vpflow.flag("DATA::NO_SCAN",aurostd::args2flag(argv,cmds,"--no_scan")); //DX 9/1/17 - SGDATA + JSON
    if(aurostd::args2attachedflag(argv,"--edata=")){ //DX 9/1/17 - SGDATA + JSON
      vpflow.args2addattachedscheme(argv,cmds,"DATA::TOLERANCE","--edata=","1"); //DX 9/1/17 - SGDATA + JSON
    } //DX 9/1/17 - SGDATA + JSON
    //DX 11/28/17 - MAGNETIC SYMMETRY - START
    vpflow.args2addattachedscheme(argv,cmds,"DATA::MAGNETIC","--mag=|--magnetic=|--magmom=","1"); //DX 11/28/17
    //DX 11/28/17 - MAGNETIC SYMMETRY - END
  } //DX 9/1/17 - SGDATA + JSON
  vpflow.flag("EDOS",aurostd::args2flag(argv,cmds,"--edos"));
  vpflow.flag("EFFMASS",aurostd::args2flag(argv,cmds,"--em")) ; // CAMILO
  // [OBSOLETE] vpflow.flag("EWALD",aurostd::args2flag(argv,cmds,"--ewald") && argv.at(1)=="--ewald");
  vpflow.args2addattachedscheme(argv,cmds,"EWALD","--ewald=","");
 
  //DX 8/18/17 [OBSOLETE] vpflow.flag("EQUIVALENT",aurostd::args2flag(argv,cmds,"--equivalent|--equiv|--inequivalent|--inequiv|--iatoms|--eatoms"));
  //DX 8/18/17 - Added tolerance and no_scan options to Xgroups - START
  vpflow.args2addattachedscheme(argv,cmds,"EQUIVALENT","--equivalent=|--equiv=|--inequivalent=|--inequiv=|--iatoms=|--eatoms=","");
  if(vpflow.flag("EQUIVALENT")){
    vpflow.flag("SYMMETRY::NO_SCAN",aurostd::args2flag(argv,cmds,"--no_scan"));
    if(aurostd::args2attachedflag(argv,"--equivalent=|--equiv=|--inequivalent=|--inequiv=|--iatoms=|--eatoms=")){ //DX 8/3/17
      vpflow.args2addattachedscheme(argv,cmds,"SYMMETRY::TOLERANCE","--equivalent=|--equiv=|--inequivalent=|--inequiv=|--iatoms=|--eatoms=","1");
    }
    //DX 9/21/17 - MAGNETIC SYMMETRY - START
    vpflow.args2addattachedscheme(argv,cmds,"SYMMETRY::MAGNETIC","--mag=|--magnetic=|--magmom=","1"); //DX 8/3/17
    //DX 9/21/17 - MAGNETIC SYMMETRY - END
  }
  //DX 8/18/17 - Added tolerance and no_scan options to Xgroups - END
  vpflow.flag("EXTRACT_SYMMETRY",aurostd::args2flag(argv,cmds,"--extract_symmetry|--xsymmetry"));

  // [OBSOLETE] vpflow.flag("EIGCURV",aurostd::args2flag(argv,cmds,"--eigcurv" ) ); // CAMILO
  vpflow.args2addattachedscheme(argv,cmds,"EIGCURV","--eigcurv=","./") ; // CAMILO

  //DX 8/18/17 [OBSOLETE] vpflow.flag("FGROUP",aurostd::args2flag(argv,cmds,"--factorgroup|--fgroup"));
  //DX 8/18/17 - Added tolerance and no_scan options to Xgroups - START
  vpflow.args2addattachedscheme(argv,cmds,"FGROUP","--factorgroup=|--fgroup=","");
  if(vpflow.flag("FGROUP")){
    vpflow.flag("SYMMETRY::NO_SCAN",aurostd::args2flag(argv,cmds,"--no_scan"));
    if(aurostd::args2attachedflag(argv,"--factorgroup=|--fgroup=")){ //DX 8/3/17
      vpflow.args2addattachedscheme(argv,cmds,"SYMMETRY::TOLERANCE","--factorgroup=|--fgroup=","1");
    }
    vpflow.flag("SYMMETRY::SCREEN_ONLY",aurostd::args2flag(argv,cmds,"--screen_only")); //DX 8/3/17
    //DX 9/21/17 - MAGNETIC SYMMETRY - START
    vpflow.args2addattachedscheme(argv,cmds,"SYMMETRY::MAGNETIC","--mag=|--magnetic=|--magmom=","1"); //DX 8/3/17
    //DX 9/21/17 - MAGNETIC SYMMETRY - END
  }
  //DX 8/18/17 - Added tolerance and no_scan options to Xgroups - END
  
  vpflow.flag("FRAC",aurostd::args2flag(argv,cmds,"--frac|-frac|--fractional|-fract|--fract|--direct|-direct|-f|-d"));
  vpflow.flag("FROZSL_VASPSETUP_AFLOW",aurostd::args2flag(argv,cmds,"--frozsl_vaspsetup_aflow|--frozsl_vaspsetup|--frozsl_vasp|--frozsl_setup|--phvaspsetup"));
  vpflow.flag("FROZSL_VASPSETUP_POSCAR",aurostd::args2flag(argv,cmds,"--frozsl_vaspsetup_poscar"));
  vpflow.flag("FROZSL_INPUT",aurostd::args2flag(argv,cmds,"--frozsl_input"));
  vpflow.flag("FROZSL_OUTPUT",aurostd::args2flag(argv,cmds,"--frozsl_output"));
  vpflow.flag("FROZSL_ANALYZE",aurostd::args2flag(argv,cmds,"--frozsl_analyze"));
  vpflow.flag("FROZSL_README",aurostd::args2flag(argv,cmds,"--readme=frozsl"));

  vpflow.flag("GULP",aurostd::args2flag(argv,cmds,"--gulp"));
  vpflow.flag("GETTEMP",aurostd::args2flag(argv,cmds,"--getTEMP|--getTEMPS|--Init::GetTEMPs|--gettemp|--gettemps|--getemp|--getemps"));
 
  // [OBSOLETE] vpflow.flag("KPOINTS",aurostd::args2flag(argv,cmds,"--kpoints|--kppra|-kpoints|-kppra|-k"));
  vpflow.args2addattachedscheme(argv,cmds,"KPOINTS","--kpoints=|--kppra=|-kpoints=|-kppra=|-k=","1");
  // [OBSOLETE] vpflow.flag("FLAG::XVASP_KPOINTS_DELTA",aurostd::args2flag(argv,cmds,"--delta_kpoints|--dkpoints|-dkpoints|-dk"));
  vpflow.args2addattachedscheme(argv,cmds,"FLAG::XVASP_KPOINTS_DELTA","--delta_kpoints=|--dkpoints=|-dkpoints=|-dk=","0.01"); //CO 171025
 
  vpflow.flag("KPATH",aurostd::args2flag(argv,cmds,"--kpath"));

  vpflow.flag("JOINSTRLIST",aurostd::args2flag(argv,cmds,"--join_strlist") && argv.at(1)=="--join_strlist");

  // [OBSOLETE] vpflow.flag("HKL",aurostd::args2flag(argv,cmds,"--hkl"));
  vpflow.args2addattachedscheme(argv,cmds,"HKL","--hkl=","");
  // [OBSOLETE] vpflow.flag("HKL_SEARCH_TRIVIAL",aurostd::args2flag(argv,cmds,"--hkl_search"));
  vpflow.args2addattachedscheme(argv,cmds,"HKL_SEARCH_TRIVIAL","--hkl_search=","");
  // [OBSOLETE] vpflow.flag("HKL_SEARCH_SIMPLE",aurostd::args2flag(argv,cmds,"--hkl_search_simple"));
  vpflow.args2addattachedscheme(argv,cmds,"HKL_SEARCH_SIMPLE","--hkl_search_simple=","");
  // [OBSOLETE] vpflow.flag("HKL_SEARCH_COMPLETE",aurostd::args2flag(argv,cmds,"--hkl_search_complete"));
  vpflow.args2addattachedscheme(argv,cmds,"HKL_SEARCH_COMPLETE","--hkl_search_complete=","");

  vpflow.flag("ICSD",aurostd::args2flag(argv,cmds,"--icsd"));
  vpflow.flag("ICSD_ALLLESSTHAN",aurostd::args2flag(argv,cmds,"--icsd_alllessthan"));
  vpflow.flag("ICSD_ALLMORETHAN",aurostd::args2flag(argv,cmds,"--icsd_allmorethan"));
  vpflow.flag("ICSD_BASISLT",aurostd::args2flag(argv,cmds,"--icsd_basislessthan"));
  vpflow.flag("ICSD_BASISGT",aurostd::args2flag(argv,cmds,"--icsd_basismorethan"));
  vpflow.flag("ICSD_CHECK_RAW",aurostd::args2flag(argv,cmds,"--icsd_check_raw"));
  vpflow.flag("ICSD_CHEM",aurostd::args2flag(argv,cmds,"--icsd_chem"));
  vpflow.flag("ICSD_CUBIC",aurostd::args2flag(argv,cmds,"--icsd_cubic"));
  vpflow.flag("ICSD_DENSLESSTHAN",aurostd::args2flag(argv,cmds,"--icsd_denslessthan"));
  vpflow.flag("ICSD_DENSMORETHAN",aurostd::args2flag(argv,cmds,"--icsd_densmorethan"));
  vpflow.flag("ICSD_HEXAGONAL",aurostd::args2flag(argv,cmds,"--icsd_hexagonal"));
  vpflow.flag("ICSD_ID",aurostd::args2flag(argv,cmds,"--icsd_id|--icsd_ID"));
  vpflow.flag("ICSD_MAKELABEL",aurostd::args2flag(argv,cmds,"--icsd_makelabel"));
  vpflow.flag("ICSD_LESSTHAN",aurostd::args2flag(argv,cmds,"--icsd_lessthan"));
  vpflow.flag("ICSD_LISTMETALS",aurostd::args2flag(argv,cmds,"--icsd_listmetals"));
  vpflow.flag("ICSD_MONOCLINIC",aurostd::args2flag(argv,cmds,"--icsd_monoclinic"));
  vpflow.flag("ICSD_MORETHAN",aurostd::args2flag(argv,cmds,"--icsd_morethan"));
  vpflow.flag("ICSD_N_ARY",aurostd::args2flag(argv,cmds,"--icsd_n_ary"));
  vpflow.flag("ICSD_NOBROKENBASIS",aurostd::args2flag(argv,cmds,"--icsd_nobrokenbasis|--icsd_completebasis"));
  vpflow.flag("ICSD_NOPARTIALOCC",aurostd::args2flag(argv,cmds,"--icsd_nopartialocc"));
  vpflow.flag("ICSD_ORTHORHOMBIC",aurostd::args2flag(argv,cmds,"--icsd_orthorhombic|--icsd_orthorombic|--icsd_ortorhombic|--icsd_ortorombic"));
  vpflow.flag("ICSD_PROTO",aurostd::args2flag(argv,cmds,"--icsd_proto"));
  vpflow.flag("ICSD_RHOMBOHEDRAL",aurostd::args2flag(argv,cmds,"--icsd_rhombohedral|--icsd_rombohedral|--icsd_trigonal"));
  vpflow.flag("ICSD_REMOVE_AND",aurostd::args2flag(argv,cmds,"--icsd_remove_and"));
  vpflow.flag("ICSD_REMOVE_OR",aurostd::args2flag(argv,cmds,"--icsd_remove_or"));
  vpflow.flag("ICSD_REMOVEMETALS",aurostd::args2flag(argv,cmds,"--icsd_removemetals"));
  vpflow.flag("ICSD_SG",aurostd::args2flag(argv,cmds,"--icsd_sg"));
  vpflow.flag("ICSD_SGLESSTHAN",aurostd::args2flag(argv,cmds,"--icsd_sglessthan"));
  vpflow.flag("ICSD_SGMORETHAN",aurostd::args2flag(argv,cmds,"--icsd_sgmorethan"));
  vpflow.flag("ICSD_TETRAGONAL",aurostd::args2flag(argv,cmds,"--icsd_tetragonal"));
  vpflow.flag("ICSD_TRICLINIC",aurostd::args2flag(argv,cmds,"--icsd_triclinic"));
  vpflow.flag("ICSD_TRIGONAL",aurostd::args2flag(argv,cmds,"--icsd_trigonal"));
  vpflow.flag("ICSD_TRI",aurostd::args2flag(argv,cmds,"--icsd_tri|--icsd_TRI"));
  vpflow.flag("ICSD_MCL",aurostd::args2flag(argv,cmds,"--icsd_mcl|--icsd_MCL"));
  vpflow.flag("ICSD_MCLC",aurostd::args2flag(argv,cmds,"--icsd_mclc|--icsd_MCLC"));
  vpflow.flag("ICSD_ORC",aurostd::args2flag(argv,cmds,"--icsd_orc|--icsd_ORC"));
  vpflow.flag("ICSD_ORCC",aurostd::args2flag(argv,cmds,"--icsd_orcc|--icsd_ORCC"));
  vpflow.flag("ICSD_ORCF",aurostd::args2flag(argv,cmds,"--icsd_orcf|--icsd_ORCF"));
  vpflow.flag("ICSD_ORCI",aurostd::args2flag(argv,cmds,"--icsd_orci|--icsd_ORCI"));
  vpflow.flag("ICSD_TET",aurostd::args2flag(argv,cmds,"--icsd_tet|--icsd_TET"));
  vpflow.flag("ICSD_BCT",aurostd::args2flag(argv,cmds,"--icsd_bct|--icsd_BCT"));
  vpflow.flag("ICSD_RHL",aurostd::args2flag(argv,cmds,"--icsd_rhl|--icsd_RHL"));
  vpflow.flag("ICSD_HEX",aurostd::args2flag(argv,cmds,"--icsd_hex|--icsd_HEX"));
  vpflow.flag("ICSD_CUB",aurostd::args2flag(argv,cmds,"--icsd_cub|--icsd_CUB"));
  vpflow.flag("ICSD_FCC",aurostd::args2flag(argv,cmds,"--icsd_fcc|--icsd_FCC"));
  vpflow.flag("ICSD_BCC",aurostd::args2flag(argv,cmds,"--icsd_bcc|--icsd_BCC"));
  vpflow.flag("ICSD_UNIQUE",aurostd::args2flag(argv,cmds,"--icsd_unique"));
  vpflow.flag("ICSD2POSCAR",aurostd::args2flag(argv,cmds,"--icsd2poscar|--icsd2POSCAR"));
  vpflow.flag("ICSD2PROTO",aurostd::args2flag(argv,cmds,"--icsd2proto"));
  vpflow.flag("ICSD2WYCK",aurostd::args2flag(argv,cmds,"--icsd2wyck|--icsdwyck"));

  vpflow.flag("IDENTICAL",aurostd::args2flag(argv,cmds,"--identical"));
  vpflow.flag("INCELL",aurostd::args2flag(argv,cmds,"--incell"));
  vpflow.flag("INCOMPACT",aurostd::args2flag(argv,cmds,"--incompact"));
  vpflow.flag("INSPHERE",aurostd::args2flag(argv,cmds,"--insphere"));
  // [OBSOLETE] vpflow.flag("INTPOL",(aurostd::args2flag(argv,cmds,"--intpol") && argv.at(1)=="--intpol"));
  vpflow.args2addattachedscheme(argv,cmds,"INTPOL","--intpol=","");

  vpflow.flag("INWS",aurostd::args2flag(argv,cmds,"--inwignerseitz|--inws"));

  // [OBSOLETE] vpflow.flag("INFLATE_LATTICE",aurostd::args2flag(argv,cmds,"--inflate_lattice|--ilattice"));
  vpflow.args2addattachedscheme(argv,cmds,"INFLATE_LATTICE","--inflate_lattice=|--ilattice=","");
  // [OBSOLETE] vpflow.flag("INFLATE_VOLUME",aurostd::args2flag(argv,cmds,"--inflate_volume|--ivolume"));
  vpflow.args2addattachedscheme(argv,cmds,"INFLATE_VOLUME","--inflate_volume=|--ivolume=","");
  
  vpflow.flag("KBAND",aurostd::args2flag(argv,cmds,"--kband"));
  // [OBSOLETE] vpflow.flag("KILL",aurostd::args2flag(argv,cmds,"--kill"));
  vpflow.args2addattachedscheme(argv,cmds,"KILL","--kill=","");

  vpflow.flag("HNF",(aurostd::args2flag(argv,cmds,"--hnf|--HNF|--hfn|--HFN")));
  vpflow.flag("HNFTOL",(aurostd::args2flag(argv,cmds,"--hnftol|--HNFTOL|--hfntol|--HFNTOL")));
  vpflow.flag("HNFCELL",aurostd::args2flag(argv,cmds,"--hnfcell|--hfncell"));
  vpflow.flag("MULTIENUMALL",aurostd::args2flag(argv,cmds,"--multienum|--enum"));
  vpflow.flag("MULTIENUMSORT",aurostd::args2flag(argv,cmds,"--multienumsort|--enumsort"));
  vpflow.flag("POSCAR2ENUM",(aurostd::args2flag(argv,cmds,"--poscar2multienum|--poscar2enum")));
  vpflow.flag("POSCAR2GULP",(aurostd::args2flag(argv,cmds,"--poscar2gulp")));
  vpflow.flag("POCC_INPUT",aurostd::args2flag(argv,cmds,"--pocc_input|--enum_input"));

  vpflow.args2addattachedscheme(argv,cmds,"JMOL","--jmol=","");

  vpflow.flag("JMOLGIF",aurostd::args2flag(argv,cmds,"--jgif"));
  vpflow.args2addattachedscheme(argv,cmds,"JUSTAFTER","--justafter=","");
  vpflow.args2addattachedscheme(argv,cmds,"JUSTBEFORE","--justbefore=","");
  vpflow.args2addattachedscheme(argv,cmds,"JUSTBETWEEN","--justbetween=","");
 
  vpflow.flag("LATTICEREDUCTION",aurostd::args2flag(argv,cmds,"--latticereduction|--latreduction"));

  vpflow.args2addattachedscheme(argv,cmds,"LTCELL","--ltcell=","");
  vpflow.args2addattachedscheme(argv,cmds,"LTCELLFV","--ltcellfv=","");

  vpflow.flag("LATTICE_TYPE",aurostd::args2flag(argv,cmds,"--lattice_type|--lattice|--lattice_crystal"));
  vpflow.flag("LATTICE_LATTICE_TYPE",aurostd::args2flag(argv,cmds,"--lattice_lattice_type|--lattice_lattice"));
  vpflow.flag("LATTICE_HISTOGRAM",aurostd::args2flag(argv,cmds,"--latticehistogram"));

  vpflow.flag("MAGNETICPARAMETERS",aurostd::args2flag(argv,cmds,"--magpara") ||aurostd::args2attachedflag(argv,cmds,"--magpara="));

  // [OBSOLETE] vpflow.flag("MAGNETICPARAMETERSHEUSLER",aurostd::args2flag(argv,cmds,"--magpara_h"));
  // [OBSOLETE] vpflow.flag("MAXATOMS",aurostd::args2flag(argv,cmds,"--maxatoms|--max_atoms|--atomsmax|--atoms_max"));
  vpflow.args2addattachedscheme(argv,cmds,"MAXATOMS","--maxatoms=|--max_atoms=|--atomsmax=|--atoms_max=","");
  
  vpflow.flag("MAKESTRLIST",aurostd::args2flag(argv,cmds,"--make_strlist") && argv.at(1)=="--make_strlist");
  
  // [OBSOLETE] vpflow.flag("MILLER",aurostd::args2flag(argv,cmds,"--miller|--MILLER"));
  // [OBSOLETE] vpflow.args2addattachedscheme(argv,cmds,"MILLER","--miller=","");
  
  vpflow.flag("MINKOWSKI_BASIS_REDUCTION",aurostd::args2flag(argv,cmds,"--minkowski_basis_reduction|--minkowski|--mink"));
  vpflow.flag("MISCIBILITY",aurostd::args2flag(argv,cmds,"--MIX|--mix|--MISCIBILITY|--miscibility|--MISCIBILE|--miscibile"));
  vpflow.flag("MOM",aurostd::args2flag(argv,cmds,"--mom"));
  vpflow.flag("MSI",aurostd::args2flag(argv,cmds,"--msi"));
  
  vpflow.flag("MULTIBZIP2",aurostd::args2flag(argv,cmds,"--multibzip2"));
  vpflow.flag("MULTIBUNZIP2",aurostd::args2flag(argv,cmds,"--multibunzip2"));
  vpflow.flag("MULTIGZIP",aurostd::args2flag(argv,cmds,"--multigzip"));
  vpflow.flag("MULTIGUNZIP",aurostd::args2flag(argv,cmds,"--multigunzip"));
  vpflow.flag("MULTISH",aurostd::args2flag(argv,cmds,"--multish"));
  vpflow.flag("MULTIZIP",aurostd::args2flag(argv,cmds,"--multizip"));
  
  vpflow.flag("NATOMS",aurostd::args2flag(argv,cmds,"--natoms|--numatoms"));
  vpflow.flag("NBONDXX",aurostd::args2flag(argv,cmds,"--nbondxx")); //CO 171025
  vpflow.flag("NAMES",(aurostd::args2flag(argv,cmds,"--names|--species") && (argv.at(1)=="--names" || argv.at(1)=="--species")));
  vpflow.flag("NANOPARTICLE",(aurostd::args2flag(argv,cmds,"--nanoparticle") && argv.at(1)=="--nanoparticle"));
  vpflow.flag("NDATA",aurostd::args2flag(argv,cmds,"--ndata"));
  vpflow.flag("NIGGLI",aurostd::args2flag(argv,cmds,"--niggli"));
  // vpflow.flag("NOSD",aurostd::args2flag(argv,cmds,"--nosd"));
  vpflow.flag("NNDIST",aurostd::args2flag(argv,cmds,"--nn|--nearestneighbour|--nearestneighbor"));
  vpflow.flag("NOORDERPARAMETER",aurostd::args2flag(argv,cmds,"--noorderparameter|--norderparameter|--noorder|--norder"));

  vpflow.flag("NUMNAMES",aurostd::args2flag(argv,cmds,"--numnames") && argv.at(1)=="--numnames");
  vpflow.flag("NSPECIES",aurostd::args2flag(argv,cmds,"--nspecies|--numspecies"));

  vpflow.flag("PEARSON_SYMBOL",aurostd::args2flag(argv,cmds,"--pearson_symbol|--pearson"));
  // vpflow.flag("POCCUPATION",aurostd::args2flag(argv,cmds,"--poccupation|--partial_occupation|--partialoccupation|--pocc"));
  // vpflow.flag("OPARAMETER",aurostd::args2flag(argv,cmds,"--oparameter|--order_parameter|--orderparameter|--opar"));

  vpflow.flag("PDB",aurostd::args2flag(argv,cmds,"--pdb"));
  vpflow.flag("PDOS",aurostd::args2flag(argv,cmds,"--pdos") && argv.at(1)=="--pdos");
  //DX 8/18/17 [OBSOLETE] vpflow.flag("PGROUP",aurostd::args2flag(argv,cmds,"--pointgroup|--pgroup"));
  //DX 8/18/17 - Added tolerance and no_scan options to Xgroups - START
  vpflow.args2addattachedscheme(argv,cmds,"PGROUP","--pointgroup=|--pgroup=","");
  if(vpflow.flag("PGROUP")){
    vpflow.flag("SYMMETRY::NO_SCAN",aurostd::args2flag(argv,cmds,"--no_scan"));
    if(aurostd::args2attachedflag(argv,"--pointgroup=|--pgroup=")){ //DX 8/3/17
      vpflow.args2addattachedscheme(argv,cmds,"SYMMETRY::TOLERANCE","--pointgroup=|--pgroup=","1");
    }
    vpflow.flag("SYMMETRY::SCREEN_ONLY",aurostd::args2flag(argv,cmds,"--screen_only")); //DX 8/3/17
  }
  //DX 8/18/17 - Added tolerance and no_scan options to Xgroups - END
  //DX 8/18/17 [OBSOLETE] vpflow.flag("PGROUPX",aurostd::args2flag(argv,cmds,"--pointgroup_crystal|--pgroup_crystal|--pgroup_xtal|--pgroupx|--pgroupX"));
  //DX 8/18/17 - Added tolerance and no_scan options to Xgroups - START
  vpflow.args2addattachedscheme(argv,cmds,"PGROUPX","--pointgroup_crystal=|--pgroup_crystal=|--pgroup_xtal=|--pgroupx=|--pgroupX=","");
  if(vpflow.flag("PGROUPX")){
    vpflow.flag("SYMMETRY::NO_SCAN",aurostd::args2flag(argv,cmds,"--no_scan"));
    if(aurostd::args2attachedflag(argv,"--pointgroup_crystal=|--pgroup_crystal=|--pgroup_xtal=|--pgroupx=|--pgroupX=")){ //DX 8/3/17
      vpflow.args2addattachedscheme(argv,cmds,"SYMMETRY::TOLERANCE","--pointgroup_crystal=|--pgroup_crystal=|--pgroup_xtal=|--pgroupx=|--pgroupX=","1");
    }
    vpflow.flag("SYMMETRY::SCREEN_ONLY",aurostd::args2flag(argv,cmds,"--screen_only")); //DX 8/3/17
    //DX 9/21/17 - MAGNETIC SYMMETRY - START
    vpflow.args2addattachedscheme(argv,cmds,"SYMMETRY::MAGNETIC","--mag=|--magnetic=|--magmom=","1"); //DX 8/3/17
    //DX 9/21/17 - MAGNETIC SYMMETRY - END
  }
  //DX 8/18/17 - Added tolerance and no_scan options to Xgroups - END
  //DX 8/18/17 [OBSOLETE] vpflow.flag("PGROUPK",aurostd::args2flag(argv,cmds,"--pointgroupklattice|--pgroupk"));
  //DX 8/18/17 - Added tolerance and no_scan options to Xgroups - START
  vpflow.args2addattachedscheme(argv,cmds,"PGROUPK","--pointgroupklattice=|--pgroupk=","");
  if(vpflow.flag("PGROUPK")){
    vpflow.flag("SYMMETRY::NO_SCAN",aurostd::args2flag(argv,cmds,"--no_scan"));
    if(aurostd::args2attachedflag(argv,"--pointgroupklattice=|--pgroupk=")){ //DX 8/3/17
      vpflow.args2addattachedscheme(argv,cmds,"SYMMETRY::TOLERANCE","--pointgroupklattice=|--pgroupk=","1");
    }
    vpflow.flag("SYMMETRY::SCREEN_ONLY",aurostd::args2flag(argv,cmds,"--screen_only")); //DX 8/3/17
  }
  //DX 8/18/17 - Added tolerance and no_scan options to Xgroups - END
  //DX 12/5/17 - Added pgroupk_xtal - START
  vpflow.args2addattachedscheme(argv,cmds,"PGROUPK_XTAL","--pointgroupkcrystal=|--pgroupk_xtal=","");
  if(vpflow.flag("PGROUPK_XTAL")){
    vpflow.flag("SYMMETRY::NO_SCAN",aurostd::args2flag(argv,cmds,"--no_scan"));
    if(aurostd::args2attachedflag(argv,"--pointgroupkcrystal=|--pgroupk_xtal=")){ //DX 8/3/17
      vpflow.args2addattachedscheme(argv,cmds,"SYMMETRY::TOLERANCE","--pointgroupkcrystal=|--pgroupk_xtal=","1");
    }
    vpflow.flag("SYMMETRY::SCREEN_ONLY",aurostd::args2flag(argv,cmds,"--screen_only")); //DX 8/3/17
  }
  //DX 12/5/17 - Added pgroupk_xtal - END
  vpflow.flag("PLANEDENS",aurostd::args2flag(argv,cmds,"--planedens") && argv.at(1)=="--planedens");
  // [OBSOLETE]  vpflow.flag("PLATON",aurostd::args2flag(argv,cmds,"--platon") && argv.at(1)=="--platon");
  vpflow.args2addattachedscheme(argv,cmds,"PLATON","--platon=",""); 

  vpflow.args2addattachedscheme(argv,cmds,"PLOT_BAND","--plotband=","./");
  vpflow.args2addattachedscheme(argv,cmds,"PLOT_BANDSPINSPLIT","--plotband_spinsplit=","./");
  vpflow.args2addattachedscheme(argv,cmds,"PLOT_BAND2","--plotband2=","./");
  vpflow.args2addattachedscheme(argv,cmds,"PLOT_BANDDOS","--plotbanddos=","./");
  vpflow.args2addattachedscheme(argv,cmds,"PLOT_DOS","--plotdos=","./");
  vpflow.args2addattachedscheme(argv,cmds,"PLOT_DOSWEB","--plotdosweb=","./");
  vpflow.args2addattachedscheme(argv,cmds,"PLOT_PEDOS","--plotpedos=","./,1");
  vpflow.args2addattachedscheme(argv,cmds,"PLOT_PEDOSALL","--plotpedosall=","./");
  vpflow.args2addattachedscheme(argv,cmds,"PLOT_PEDOSALL_AFLOWLIB","--plotpedos_nonequivalent=","./");
  
  vpflow.flag("PLOTPHDISP",aurostd::args2flag(argv,cmds,"--plotphonondispersion|--pphdis"));

  vpflow.flag("POCC",aurostd::args2flag(argv,cmds,"--pocc") && argv.at(1)=="--pocc");
  vpflow.args2addattachedscheme(argv,cmds,"POCC_DOS","--pocc_dos=","./");
  vpflow.args2addattachedscheme(argv,cmds,"POCC_MAG","--pocc_mag=","./");
  vpflow.args2addattachedscheme(argv,cmds,"POCC_BANDGAP","--pocc_bandgap=","./");

  vpflow.flag("POSCAR",aurostd::args2flag(argv,cmds,"--poscar"));
  vpflow.flag("POSCAR2AFLOWIN",aurostd::args2flag(argv,cmds,"--poscar2aflowin|--poscaraflowin|--poscar2aflow|--poscaraflow"));
  vpflow.flag("POSCAR2WYCKOFF",aurostd::args2flag(argv,cmds,"--poscar2wyckoff"));
  //corey
  vpflow.args2addattachedscheme(argv,cmds,"PREPARE_CHGCAR_4_JMOL","--prepare_chgcar_4_jmol=|--prep4jmol=","");  //can be multiple
  if(vpflow.flag("PREPARE_CHGCAR_4_JMOL")) {
    vpflow.flag("PREPARE_CHGCAR_4_JMOL::USAGE",aurostd::args2flag(argv,cmds,"--usage"));  // usage
    vpflow.args2addattachedscheme(argv,cmds,"PREPARE_CHGCAR_4_JMOL::OUTCAR","--outcar=",""); //singular
    vpflow.flag("PREPARE_CHGCAR_4_JMOL::ZIP",aurostd::args2flag(argv,cmds,"--zip"));
  }
  vpflow.flag("PRIM",aurostd::args2flag(argv,cmds,"--prim|--prim0"));
  vpflow.flag("PRIM1",aurostd::args2flag(argv,cmds,"--prim1"));
  vpflow.flag("PRIM2",aurostd::args2flag(argv,cmds,"--prim2"));
  vpflow.flag("PRIM3",aurostd::args2flag(argv,cmds,"--prim3"));
  // [OBSOLETE] vpflow.flag("PRIMJ",aurostd::args2flag(argv,cmds,"--primj"));
  // [OBSOLETE] vpflow.flag("PROTO", (aurostd::args2attachedflag(argv,cmds,"--proto=")) ||           // --proto=123:A:B:C ..
  // [OBSOLETE]            (aurostd::args2attachedflag(argv,cmds,"--proto_icsd=")));                 // --proto_icsd=Gd1Mn2Si2_ICSD_54947
  // [OBSOLETE] vpflow.flag("PROTO_AFLOW",(aurostd::args2attachedflag(argv,cmds,"--aflow_proto=")) ||   // --aflow_proto=123:A:B:C ..
  // [OBSOLETE]            (aurostd::args2attachedflag(argv,cmds,"--aflow_proto_icsd=")));         // --aflow_proto_icsd=Gd1Mn2Si2_ICSD_54947
  
  vpflow.args2addattachedscheme(argv,cmds,"PROTO","--proto=|--proto_icsd=","");                                                    // --proto=123:A:B:C --proto_icsd=Gd1Mn2Si2_ICSD_54947
  vpflow.args2addattachedscheme(argv,cmds,"PARAMS","--params=|--parameters=","");                                                  // --proto=123:A:B:C --proto_icsd=Gd1Mn2Si2_ICSD_54947
  if(vpflow.flag("PROTO")) {
    vpflow.flag("PROTO::VASP",aurostd::args2flag(argv,cmds,"--vasp") || (!aurostd::args2flag(argv,cmds,"--qe") && !aurostd::args2flag(argv,cmds,"--abinit") && !aurostd::args2flag(argv,cmds,"--aims")));
    vpflow.flag("PROTO::QE",aurostd::args2flag(argv,cmds,"--qe"));
    vpflow.flag("PROTO::ABINIT",aurostd::args2flag(argv,cmds,"--abinit"));
    vpflow.flag("PROTO::AIMS",aurostd::args2flag(argv,cmds,"--aims"));
    vpflow.flag("PROTO::HEX",aurostd::args2flag(argv,cmds,"--hex"));
    vpflow.flag("PROTO::RHL",aurostd::args2flag(argv,cmds,"--rhl"));
  } 
  vpflow.args2addattachedscheme(argv,cmds,"PROTO_AFLOW","--aflow_proto=|--aflow_proto_icsd=","");      // --aflow_proto=123:A:B:C  --aflow_proto_icsd=Gd1Mn2Si2_ICSD_54947 
  if(vpflow.flag("PROTO_AFLOW")) {
    string vlist="";
    vpflow.args2addattachedscheme(argv,cmds,"PROTO_AFLOW::PRESSURE","--pressure=|--PRESSURE=|--pstress=|--PSTRESS=","");
    vpflow.args2addattachedscheme(argv,cmds,"PROTO_AFLOW::POTENTIAL","--potential=|--pp=|--potentials=","");
    if(vpflow.flag("PROTO_AFLOW::POTENTIAL")) vlist+="--potential="+vpflow.getattachedscheme("PROTO_AFLOW::POTENTIAL")+" "; // recursion is GNU's pleasure (SC2014)
 
    //AFLOW modules
    vpflow.args2addattachedscheme(argv,cmds,"PROTO_AFLOW::MODULE","--module=",""); //CO 180214
    if(vpflow.flag("PROTO_AFLOW::MODULE")) vlist+="--module="+vpflow.getattachedscheme("PROTO_AFLOW::MODULE")+" ";             // recursion is GNU's pleasure (SC2014)
    vpflow.args2addattachedscheme(argv,cmds,"PROTO_AFLOW::APL_SUPERCELL","--apl_supercell=",""); //CO 180214
    if(vpflow.flag("PROTO_AFLOW::APL_SUPERCELL")) vlist+="--apl_supercell="+vpflow.getattachedscheme("PROTO_AFLOW::APL_SUPERCELL")+" ";             // recursion is GNU's pleasure (SC2014)

    vpflow.args2addattachedscheme(argv,cmds,"PROTO_AFLOW::POTIM","--potim=|--POTIM=","");
    if(vpflow.flag("PROTO_AFLOW::POTIM")) vlist+="--potim="+vpflow.getattachedscheme("PROTO_AFLOW::POTIM")+" ";             // recursion is GNU's pleasure (SC2014)
    vpflow.args2addattachedscheme(argv,cmds,"PROTO_AFLOW::RELAX_TYPE","--relax_type=|--RELAX_TYPE=","");   //CO 180214
    if(vpflow.flag("PROTO_AFLOW::RELAX_TYPE")) vlist+="--relax_type="+vpflow.getattachedscheme("PROTO_AFLOW::RELAX_TYPE")+" "; // recursion is GNU's pleasure (SC2014)  //CO 180214
    vpflow.args2addattachedscheme(argv,cmds,"PROTO_AFLOW::RELAX_MODE","--relax_mode=|--RELAX_MODE=","");
    if(vpflow.flag("PROTO_AFLOW::RELAX_MODE")) vlist+="--relax_mode="+vpflow.getattachedscheme("PROTO_AFLOW::RELAX_MODE")+" "; // recursion is GNU's pleasure (SC2014)
    vpflow.args2addattachedscheme(argv,cmds,"PROTO_AFLOW::PRECISION","--prec=|--PREC=|--precision=|--PRECISION=","");
    if(vpflow.flag("PROTO_AFLOW::PRECISION")) vlist+="--precision="+vpflow.getattachedscheme("PROTO_AFLOW::PRECISION")+" "; // recursion is GNU's pleasure (SC2014)
    vpflow.args2addattachedscheme(argv,cmds,"PROTO_AFLOW::ALGORITHM","--algo=|--ALGO=|--algorithm=|--ALGORITHM=","");
    if(vpflow.flag("PROTO_AFLOW::ALGORITHM")) vlist+="--algorithm="+vpflow.getattachedscheme("PROTO_AFLOW::ALGORITHM")+" "; // recursion is GNU's pleasure (SC2014)
    vpflow.args2addattachedscheme(argv,cmds,"PROTO_AFLOW::METAGGA","--metagga=|--METAGGA=","");
    if(vpflow.flag("PROTO_AFLOW::METAGGA")) vlist+="--metagga="+vpflow.getattachedscheme("PROTO_AFLOW::METAGGA")+" "; // recursion is GNU's pleasure (SC2014)
    vpflow.args2addattachedscheme(argv,cmds,"PROTO_AFLOW::IVDW","--ivdw=|--IVDW=","");
    if(vpflow.flag("PROTO_AFLOW::IVDW")) vlist+="--ivdw="+vpflow.getattachedscheme("PROTO_AFLOW::IVDW")+" "; // recursion is GNU's pleasure (SC2014)
    vpflow.args2addattachedscheme(argv,cmds,"PROTO_AFLOW::TYPE","--type=|--TYPE=","");
    if(vpflow.flag("PROTO_AFLOW::TYPE")) vlist+="--type="+vpflow.getattachedscheme("PROTO_AFLOW::TYPE")+" ";                // recursion is GNU's pleasure (SC2014)
    vpflow.args2addattachedscheme(argv,cmds,"PROTO_AFLOW::CONVERT_UNIT_CELL","--convert_unit_cell=|--CONVERT_UNIT_CELL=","");
    if(vpflow.flag("PROTO_AFLOW::CONVERT_UNIT_CELL")) vlist+="--convert_unit_cell="+vpflow.getattachedscheme("PROTO_AFLOW::CONVERT_UNIT_CELL")+" ";  // recursion is GNU's pleasure (SC2014)
    vpflow.args2addattachedscheme(argv,cmds,"PROTO_AFLOW::VOLUME_PLUS_EQUAL","--volume_plus_equal=|--VOLUME_PLUS_EQUAL=","");
    if(vpflow.flag("PROTO_AFLOW::VOLUME_PLUS_EQUAL")) vlist+="--volume_plus_equal="+vpflow.getattachedscheme("PROTO_AFLOW::VOLUME_PLUS_EQUAL")+" ";  // recursion is GNU's pleasure (SC2014)
    vpflow.args2addattachedscheme(argv,cmds,"PROTO_AFLOW::VOLUME_MULTIPLY_EQUAL","--volume_multiply_equal=|--VOLUME_MULTIPLY_EQUAL=","");
    if(vpflow.flag("PROTO_AFLOW::VOLUME_MULTIPLY_EQUAL")) vlist+="--volume_multiply_equal="+vpflow.getattachedscheme("PROTO_AFLOW::VOLUME_MULTIPLY_EQUAL")+" ";   // recursion is GNU's pleasure (SC2014)
    vpflow.flag("PROTO_AFLOW::NO_VOLUME_ADJUSTMENT",aurostd::args2flag(argv,cmds,"--no_volume_adjustment")); //CO 180214
    if(vpflow.flag("PARAMS") && !(vpflow.flag("PROTO_AFLOW::VOLUME_PLUS_EQUAL") || vpflow.flag("PROTO_AFLOW::VOLUME_MULTIPLY_EQUAL"))){  //CO 180214
      vpflow.flag("PROTO_AFLOW::NO_VOLUME_ADJUSTMENT",TRUE);
    }
    if(vpflow.flag("PROTO_AFLOW::NO_VOLUME_ADJUSTMENT")) vlist+="--no_volume_adjustment ";  //CO 180214
    vpflow.args2addattachedscheme(argv,cmds,"PROTO_AFLOW::EDIFFG","--ediffg=|--EDIFFG=","");
    if(vpflow.flag("PROTO_AFLOW::EDIFFG")) vlist+="--ediffg="+vpflow.getattachedscheme("PROTO_AFLOW::EDIFFG")+" ";  // recursion is GNU's pleasure (SC2014)
 
    vpflow.args2addattachedscheme(argv,cmds,"PROTO_AFLOW::KPPRA","--kppra=|--KPPRA=","");
    if(vpflow.flag("PROTO_AFLOW::KPPRA")) vlist+="--kppra="+vpflow.getattachedscheme("PROTO_AFLOW::KPPRA")+" ";   // recursion is GNU's pleasure (SC2017)
    vpflow.args2addattachedscheme(argv,cmds,"PROTO_AFLOW::KPPRA_STATIC","--kppra_static=|--KPPRA_STATIC=","");
    if(vpflow.flag("PROTO_AFLOW::KPPRA_STATIC")) vlist+="--kppra_static="+vpflow.getattachedscheme("PROTO_AFLOW::KPPRA_STATIC")+" ";   // recursion is GNU's pleasure (SC2017)
    vpflow.args2addattachedscheme(argv,cmds,"PROTO_AFLOW::BANDS_GRID","--bands_grid=|--BANDS_GRID=","");
    if(vpflow.flag("PROTO_AFLOW::BANDS_GRID")) vlist+="--bands_grid="+vpflow.getattachedscheme("PROTO_AFLOW::BANDS_GRID")+" ";   // recursion is GNU's pleasure (SC2017)

    vpflow.args2addattachedscheme(argv,cmds,"PROTO_AFLOW::ENMAX_MULTIPLY","--enmax_multiply=|--ENMAX_MULTIPLY=","");
    if(vpflow.flag("PROTO_AFLOW::ENMAX_MULTIPLY")) vlist+="--enmax_multiply="+vpflow.getattachedscheme("PROTO_AFLOW::ENMAX_MULTIPLY")+" ";   // recursion is GNU's pleasure (SC2017)

    vpflow.flag("PROTO_AFLOW::USAGE",aurostd::args2flag(argv,cmds,"--usage"));
    vpflow.flag("PROTO_AFLOW::POTENTIAL_COMPLETE",aurostd::args2flag(argv,cmds,"--potential_complete|--potcomplete|--potentials_complete|--potscomplete|--potc"));
    if(vpflow.flag("PROTO_AFLOW::POTENTIAL_COMPLETE")) vlist+="--potential_complete ";                       // recursion is GNU's pleasure (SC2014)
    vpflow.flag("PROTO_AFLOW::MISSING",aurostd::args2flag(argv,cmds,"--missing"));
    if(vpflow.flag("PROTO_AFLOW::MISSING")) vlist+="--missing ";                                             // recursion is GNU's pleasure (SC2014)
    vpflow.flag("PROTO_AFLOW::NOAUTOPP",aurostd::args2flag(argv,cmds,"--noautopp"));
    if(vpflow.flag("PROTO_AFLOW::NOAUTOPP")) vlist+="--noautopp ";                                           // recursion is GNU's pleasure (SC2014)
    vpflow.flag("PROTO_AFLOW::LDAU",aurostd::args2flag(argv,cmds,"--ldau|--ldau2"));
    if(vpflow.flag("PROTO_AFLOW::LDAU")) vlist+="--ldau ";                                                   // recursion is GNU's pleasure (SC2014)
    vpflow.flag("PROTO_AFLOW::NOLDAU",aurostd::args2flag(argv,cmds,"--noldau|--noldau2"));
    if(vpflow.flag("PROTO_AFLOW::NOLDAU")) vlist+="--noldau ";                                               // recursion is GNU's pleasure (SC2014)
    vpflow.flag("PROTO_AFLOW::BANDS_CALCULATION",aurostd::args2flag(argv,cmds,"--bands|--band"));
    if(vpflow.flag("PROTO_AFLOW::BANDS_CALCULATION")) vlist+="--bands ";                                     // recursion is GNU's pleasure (SC2014)
    vpflow.flag("PROTO_AFLOW::NEGLECT_NOMIX",aurostd::args2flag(argv,cmds,"--neglect_nomix|--neglectnomix"));
    if(vpflow.flag("PROTO_AFLOW::NEGLECT_NOMIX")) vlist+="--neglect_nomix ";                                 // recursion is GNU's pleasure (SC2014)
    vpflow.flag("PROTO_AFLOW::STDOUT",aurostd::args2flag(argv,cmds,"--stdout"));
    if(vpflow.flag("PROTO_AFLOW::STDOUT")) vlist+="--stdout ";                                               // recursion is GNU's pleasure (SC2014)
    vpflow.flag("PROTO_AFLOW::QE",aurostd::args2flag(argv,cmds,"--qe"));
    if(vpflow.flag("PROTO_AFLOW::QE")) vlist+="--qe ";                                                       // recursion is GNU's pleasure (SC2014)
    vpflow.flag("PROTO_AFLOW::ABINIT",aurostd::args2flag(argv,cmds,"--abinit"));
    if(vpflow.flag("PROTO_AFLOW::ABINIT")) vlist+="--abinit ";                                               // recursion is GNU's pleasure (SC2014)
    vpflow.flag("PROTO_AFLOW::AIMS",aurostd::args2flag(argv,cmds,"--aims"));
    if(vpflow.flag("PROTO_AFLOW::AIMS")) vlist+="--aims ";                                                   // recursion is GNU's pleasure (SC2016)
    vpflow.flag("PROTO_AFLOW::HEX",aurostd::args2flag(argv,cmds,"--hex"));
    if(vpflow.flag("PROTO_AFLOW::HEX")) vlist+="--hex ";                                                   // recursion is GNU's pleasure (SC2016)
    vpflow.flag("PROTO_AFLOW::RHL",aurostd::args2flag(argv,cmds,"--rhl"));
    if(vpflow.flag("PROTO_AFLOW::RHL")) vlist+="--rhl ";                                                   // recursion is GNU's pleasure (SC2016)
    vpflow.flag("PROTO_AFLOW::VASP",aurostd::args2flag(argv,cmds,"--vasp"));
    if(vpflow.flag("PROTO_AFLOW::VASP")) vlist+="--vasp ";                                                   // recursion is GNU's pleasure (SC2014)
    vpflow.flag("PROTO_AFLOW::HTQC",aurostd::args2flag(argv,"--htqc"));
    if(vpflow.flag("PROTO_AFLOW::HTQC")) vlist+="--htqc ";                                                   // recursion is GNU's pleasure (SC2014)
    vpflow.flag("PROTO_AFLOW::LIST",aurostd::args2flag(argv,cmds,"--list"));
    //     if(vpflow.flag("PROTO_AFLOW::LIST")) vlist+="--listDEBUG ";
    vpflow.addattachedscheme("PROTO_AFLOW::LIST_VCMD",vlist,vpflow.flag("PROTO_AFLOW::LIST"));               // recursion is GNU's pleasure (SC2014)
    vpflow.flag("KPOINTS",FALSE);   // PRIORITIES
    vpflow.flag("QE",FALSE);        // PRIORITIES
    vpflow.flag("ABINIT",FALSE);    // PRIORITIES
    vpflow.flag("AIMS",FALSE);      // PRIORITIES
  }  

  // [OBSOLETE]  vpflow.flag("PROTOCLASSIFY",aurostd::args2flag(argv,cmds,"--protoclassify|--PROTOCLASSIFY"));
  // [OBSOLETE]  vpflow.flag("PROTOCLASSIFYH",aurostd::args2flag(argv,cmds,"--protoclassify_h|--PROTOCLASSIFY_H"));
  // [OBSOLETE]  vpflow.flag("PROTOCLASSIFYM",aurostd::args2flag(argv,cmds,"--protoclassify_m|--PROTOCLASSIFY_M"));
  // [OBSOLETE]  vpflow.flag("PROTOTYPER",aurostd::args2flag(argv,cmds,"--prototyper|--PROTOTYPER"));
  // [OBSOLETE]  vpflow.flag("PROTOTYPERH",aurostd::args2flag(argv,cmds,"--prototyper_h|--PROTOTYPER_H"));
  // [OBSOLETE]  vpflow.flag("PROTOTYPERM",aurostd::args2flag(argv,cmds,"--prototyper_m|--PROTOTYPER_M"));
  // [OBSOLETE]  vpflow.flag("PROTOTYPER_ANGLE",aurostd::args2flag(argv,cmds,"--prototyper_angle|--PROTOTYPER_ANGLE"));
  // [OBSOLETE]  vpflow.flag("PROTOTYPER_MULT",aurostd::args2flag(argv,cmds,"--prototypermult|--PROTOTYPERMULT"));
  vpflow.flag("PROTOINITIAL",aurostd::args2flag(argv,cmds,"--initial"));

  vpflow.flag("QE",aurostd::args2flag(argv,cmds,"--qe") && !vpflow.flag("PROTO_AFLOW") && !vpflow.flag("PROTO"));
  vpflow.flag("ABINIT",aurostd::args2flag(argv,cmds,"--abinit") && !vpflow.flag("PROTO_AFLOW") && !vpflow.flag("PROTO"));
  vpflow.flag("AIMS",aurostd::args2flag(argv,cmds,"--aims") && !vpflow.flag("PROTO_AFLOW") && !vpflow.flag("PROTO"));
  // [OBSOLETE]  vpflow.flag("QSUB",aurostd::args2flag(argv,cmds,"--qsub|--qstart"));
  vpflow.args2addattachedscheme(argv,cmds,"QSUB","--qsub=|--qstart=|--bsub=|--sbatch=","");
  // [OBSOLETE]  vpflow.flag("QDEL",aurostd::args2flag(argv,cmds,"--qdel|--scancel|--bkill"));
  vpflow.args2addattachedscheme(argv,cmds,"QDEL","--qdel=|--scancel=|--bkill=","");
  
  vpflow.flag("QMVASP",aurostd::args2flag(argv,cmds,"--qmvasp"));
  // [OBSOLETE] vpflow.flag("BSUB",aurostd::args2flag(argv,cmds,"--bsub"));
  // [OBSOLETE] vpflow.args2addattachedscheme(argv,cmds,"BSUB","--bsub=","");
  // [OBSOLETE] vpflow.flag("SBATCH",aurostd::args2flag(argv,cmds,"--sbatch"));
  // [OBSOLETE] vpflow.args2addattachedscheme(argv,cmds,"SBATCH","--sbatch=","");

  // [OBSOLETE]  vpflow.flag("RASMOL",aurostd::args2flag(argv,cmds,"--rasmol"));
  vpflow.args2addattachedscheme(argv,cmds,"RASMOL","--rasmol=","");

  vpflow.flag("RAYTRACE",(aurostd::args2flag(argv,cmds,"--raytrace") && argv.at(1)=="--raytrace"));
  vpflow.flag("RBANAL",aurostd::args2flag(argv,cmds,"--rbanal") && argv.at(1)=="--rbanal");
  vpflow.flag("RBDIST",aurostd::args2flag(argv,cmds,"--rbdist") && argv.at(1)=="--rbdist");
  // [OBSOLETE]  vvpflow.flag("RDF",(aurostd::args2flag(argv,cmds,"--rdf") && argv.at(1)=="--rdf"));
  vpflow.args2addattachedscheme(argv,cmds,"RDF","--rdf=","");
  // [OBSOLETE]  vpflow.flag("RDFCMP",(aurostd::args2flag(argv,cmds,"--rdfcmp") && argv.at(1)=="--rdfcmp"));
  vpflow.args2addattachedscheme(argv,cmds,"RDFCMP","--rdfcmp=","");

  vpflow.flag("RMATOM",aurostd::args2flag(argv,cmds,"--rm_atom") && argv.at(1)=="--rm_atom");
  vpflow.flag("RMCOPIES",aurostd::args2flag(argv,cmds,"--rm_copies") && argv.at(1)=="--rm_copies");
  vpflow.flag("RSM",aurostd::args2flag(argv,cmds,"--rsm"));

  // [OBSOLETE]  vpflow.flag("SCALE",aurostd::args2flag(argv,cmds,"--scale"));
  vpflow.args2addattachedscheme(argv,cmds,"SCALE","--scale=","0.0");
  
  vpflow.flag("SD",(aurostd::args2flag(argv,cmds,"--sd") && argv.at(1)=="--sd"));
  vpflow.flag("SETCM",(aurostd::args2flag(argv,cmds,"--setcm") && argv.at(1)=="--setcm"));
  vpflow.flag("SETORIGIN",(aurostd::args2flag(argv,cmds,"--setorigin") && argv.at(1)=="--setorigin"));
  vpflow.flag("SEWALD",aurostd::args2flag(argv,cmds,"--sewald") && argv.at(1)=="--sewald");
  // [OBSOLETE] vpflow.flag("SHELL",(aurostd::args2flag(argv,cmds,"--shell") && argv.at(1)=="--shell"));
  vpflow.args2addattachedscheme(argv,cmds,"SHELL","--shell=","");
  // [OBSOLETE] vpflow.flag("SHIFT",(aurostd::args2flag(argv,cmds,"--shift") && argv.at(1)=="--shift"));
  vpflow.args2addattachedscheme(argv,cmds,"SHIFT","--shift=","");
  vpflow.flag("SG",aurostd::args2flag(argv,cmds,"--sg|-sg"));
  //DX 8/18/17 [OBSOLETE] vpflow.flag("SGROUP",aurostd::args2flag(argv,cmds,"--spacegroup|--sgroup"));
  //DX 8/18/17 - Added tolerance and no_scan options to Xgroups - START
  vpflow.args2addattachedscheme(argv,cmds,"SGROUP","--spacegroup=|--sgroup=","");
  if(vpflow.flag("SGROUP")){
    vpflow.flag("SYMMETRY::NO_SCAN",aurostd::args2flag(argv,cmds,"--no_scan"));
    if(aurostd::args2attachedflag(argv,"--spacegroup=|--sgroup=")){ //DX 8/3/17
      vpflow.args2addattachedscheme(argv,cmds,"SYMMETRY::TOLERANCE","--spacegroup=|--sgroup=","1");
    }
    vpflow.args2addattachedscheme(argv,cmds,"SYMMETRY::SGROUP_RADIUS","--radius=","");  //DX 8/3/17
    vpflow.flag("SYMMETRY::SCREEN_ONLY",aurostd::args2flag(argv,cmds,"--screen_only")); //DX 8/3/17
    //DX 9/21/17 - MAGNETIC SYMMETRY - START
    vpflow.args2addattachedscheme(argv,cmds,"SYMMETRY::MAGNETIC","--mag=|--magnetic=|--magmom=","1"); //DX 8/3/17
    //DX 9/21/17 - MAGNETIC SYMMETRY - END
  }
  //DX 8/18/17 - Added tolerance and no_scan options to Xgroups - END
  // vpflow.flag("SPLINE",aurostd::args2flag(argv,cmds,"--spline") && argv.at(1)=="--spline");

  // [OBSOLETE] vpflow.flag("SG::AFLOW",aurostd::args2flag(argv,cmds,"--aflowSG") && argv.at(1)=="--aflowSG");
  // [OBSOLETE] vpflow.flag("SG::AFLOW_LABEL",aurostd::args2flag(argv,cmds,"--aflowSG_label") && argv.at(1)=="--aflowSG_label");
  // [OBSOLETE] vpflow.flag("SG::AFLOW_NUMBER",aurostd::args2flag(argv,cmds,"--aflowSG_number|--aflowSGn|--aflowSGN") && argv.at(1)=="--aflowSG_number");
  vpflow.args2addattachedscheme(argv,cmds,"SG::AFLOW","--aflowSG=",""); 
  vpflow.args2addattachedscheme(argv,cmds,"SG::AFLOW_LABEL","--aflowSG_label=",""); 
  vpflow.args2addattachedscheme(argv,cmds,"SG::AFLOW_NUMBER","--aflowSG_number=",""); 
  //DX 9/26/17 - Create flags for SG functions - START
  if(vpflow.flag("SG::AFLOW") || vpflow.flag("SG::AFLOW_LABEL") || vpflow.flag("SG::AFLOW_NUMBER")){
    vpflow.flag("SG::NO_SCAN",aurostd::args2flag(argv,cmds,"--no_scan"));
    if(aurostd::args2attachedflag(argv,"--aflowSG=|--aflowSG_label=|--aflowSG_number=")){
      vpflow.args2addattachedscheme(argv,cmds,"SG::TOLERANCE","--aflowSG=|--aflowSG_label=|--aflowSG_number=","1");
    }
    //DX 9/21/17 - MAGNETIC SYMMETRY - START
    vpflow.args2addattachedscheme(argv,cmds,"SG::MAGNETIC","--mag=|--magnetic=|--magmom=","1"); //DX 8/3/17
    //DX 9/21/17 - MAGNETIC SYMMETRY - END
  }
  //DX 9/26/17 - Create flags for SG functions - END
 
  // [OBSOLETE] vpflow.flag("SG::PLATON",aurostd::args2flag(argv,cmds,"--platonSG") && argv.at(1)=="--platonSG");
  // [OBSOLETE] vpflow.flag("SG::PLATON_LABEL",aurostd::args2flag(argv,cmds,"--platonSG_label") && argv.at(1)=="--platonSG_label");
  // [OBSOLETE] vpflow.flag("SG::PLATON_NUMBER",aurostd::args2flag(argv,cmds,"--platonSG_number|--platonSGn|--platonSGN") && argv.at(1)=="--platonSG_number");
  vpflow.args2addattachedscheme(argv,cmds,"SG::PLATON","--platonSG=",""); 
  vpflow.args2addattachedscheme(argv,cmds,"SG::PLATON_LABEL","--platonSG_label=",""); 
  vpflow.args2addattachedscheme(argv,cmds,"SG::PLATON_NUMBER","--platonSG_number=",""); 
  //DX 9/26/17 - Create flags for SG functions - START
  if(vpflow.flag("SG::PLATON") || vpflow.flag("SG::PLATON_LABEL") || vpflow.flag("SG::PLATON_NUMBER")){
    if(aurostd::args2attachedflag(argv,"--platonSG=|--platonSG_label=|--platonSG_number=")){
      vpflow.args2addattachedscheme(argv,cmds,"SG::TOLERANCE","--platonSG=|--platonSG_label=|--platonSG_number=","1");
    }
  } 
  //DX 9/26/17 - Create flags for SG functions - END 
 
  // [OBSOLETE] vpflow.flag("SG::FINDSYM",aurostd::args2flag(argv,cmds,"--findsymSG") && argv.at(1)=="--findsymSG");
  // [OBSOLETE] vpflow.flag("SG::FINDSYM_LABEL",aurostd::args2flag(argv,cmds,"--findsymSG_label") && argv.at(1)=="--findsymSG_label");
  // [OBSOLETE] vpflow.flag("SG::FINDSYM_NUMBER",aurostd::args2flag(argv,cmds,"--findsymSG_number|--findsymSGn|--findsymSGN") && argv.at(1)=="--findsymSG_number");
  vpflow.args2addattachedscheme(argv,cmds,"SG::FINDSYM","--findsymSG=",""); 
  vpflow.args2addattachedscheme(argv,cmds,"SG::FINDSYM_LABEL","--findsymSG_label=",""); 
  vpflow.args2addattachedscheme(argv,cmds,"SG::FINDSYM_NUMBER","--findsymSG_number=",""); 
  vpflow.args2addattachedscheme(argv,cmds,"SG::FINDSYM_PRINT","--findsym_print=","");
  vpflow.args2addattachedscheme(argv,cmds,"SG::FINDSYM_EXEC","--findsym=","");
  //DX 9/26/17 - Create flags for SG functions - START
  if(vpflow.flag("SG::FINDSYM") || vpflow.flag("SG::FINDSYM_LABEL") || vpflow.flag("SG::FINDSYM_NUMBER") ||
     vpflow.flag("SG::FINDSYM_PRINT") || vpflow.flag("SG::FINDSYM_EXEC")){
    if(aurostd::args2attachedflag(argv,"--findsymSG=|--findsymSG_label=|--findsymSG_number=|--findsym_print=|--findsym=")){
      vpflow.args2addattachedscheme(argv,cmds,"SG::TOLERANCE","--findsymSG=|--findsymSG_label=|--findsymSG_number=|--findsym_print=|--findsym=","1");
    }
  }
  //DX 9/26/17 - Create flags for SG functions - END


  vpflow.args2addattachedscheme(argv,cmds,"SGDATA","--sgdata=|--space_group_data=","");
  if(vpflow.flag("SGDATA")){
    vpflow.flag("SGDATA::NO_SCAN",aurostd::args2flag(argv,cmds,"--no_scan")); //DX 9/1/17 - SGDATA + JSON
    if(aurostd::args2attachedflag(argv,"--sgdata=|--space_group_data=")){
      vpflow.args2addattachedscheme(argv,cmds,"SGDATA::TOLERANCE","--sgdata=|--space_group_data=","1");
    }
    //DX 9/21/17 - MAGNETIC SYMMETRY - START
    vpflow.args2addattachedscheme(argv,cmds,"SGDATA::MAGNETIC","--mag=|--magnetic=|--magmom=","1"); //DX 8/3/17
    //DX 9/21/17 - MAGNETIC SYMMETRY - END
  }
  // [OBSOLETE] vpflow.flag("SLAB",aurostd::args2flag(argv,cmds,"--slab|--SLAB"));
  vpflow.args2addattachedscheme(argv,cmds,"SLAB","--slab=","");

  vpflow.flag("SOF",aurostd::args2flag(argv,cmds,"--sof"));
  vpflow.flag("SPECIES",aurostd::args2flag(argv,cmds,"--species"));
  vpflow.flag("STATDIEL",aurostd::args2flag(argv,cmds,"--statdiel") ); // CAMILO
  vpflow.flag("STDCONVCELL",aurostd::args2flag(argv,cmds,"--sc|--standard_conventional|--std_conv|--sconv"));
  vpflow.flag("STDPRIMCELL",aurostd::args2flag(argv,cmds,"--sp|--standard_primitive|--std_prim|--sprim"));
  vpflow.flag("SUMPDOS",(aurostd::args2flag(argv,cmds,"--sumpdos") && argv.at(1)=="--sumpdos"));

  // [OBSOLETE] vpflow.flag("SUPERCELL",aurostd::args2flag(argv,cmds,"--supercell"));
  vpflow.args2addattachedscheme(argv,cmds,"SUPERCELL","--supercell=","");
  // [OBSOLETE] vpflow.flag("SUPERCELLSTRLIST",aurostd::args2flag(argv,cmds,"--supercell_strlist") && argv.at(1)=="--supercell_strlist");
  vpflow.args2addattachedscheme(argv,cmds,"SUPERCELLSTRLIST","--supercell_strlist=","");

  vpflow.flag("SWAP",aurostd::args2flag(argv,cmds,"--swap"));

  vpflow.flag("TERDATA",aurostd::args2flag(argv,cmds,"--terdata") || aurostd::args2attachedflag(argv,cmds,"--terdata="));
  vpflow.flag("TERDATA_EXIST",aurostd::args2flag(argv,cmds,"--terdata_exist"));
  
  vpflow.flag("UFFENERGY",aurostd::args2flag(argv,cmds,"--uffenergy|--ue"));
  
  vpflow.flag("VASP",aurostd::args2flag(argv,cmds,"--vasp"));

  vpflow.args2addattachedscheme(argv,cmds,"VOLUME::EQUAL","--volume=","");
  vpflow.args2addattachedscheme(argv,cmds,"VOLUME::MULTIPLY_EQUAL","--volume*=","");
  vpflow.args2addattachedscheme(argv,cmds,"VOLUME::PLUS_EQUAL","--volume+=","");

  vpflow.flag("WWW",aurostd::args2flag(argv,cmds,"--web|--www|--http"));
  vpflow.flag("WYCKOFF",aurostd::args2flag(argv,cmds,"--wyckoff|--wy"));

  // [OBSOLETE] vpflow.flag("XRAY",aurostd::args2flag(argv,cmds,"--xray"));
  vpflow.args2addattachedscheme(argv,cmds,"XRAY","--xray=","");
  // [OBSOLETE] vpflow.flag("XYZ",aurostd::args2flag(argv,cmds,"--xyz"));
  vpflow.args2addattachedscheme(argv,cmds,"XYZ","--xyz=","");
  vpflow.flag("XYZWS",aurostd::args2flag(argv,cmds,"--xyzwignerseitz|--xyzws"));

  //  vpflow.flag("XXX",aurostd::args2flag(argv,cmds,"--xxx"));
  vpflow.flag("XFIXX",(aurostd::args2flag(argv,cmds,"--xfixX|--xfixx") && argv.size()==4));

  vpflow.flag("XXX",(aurostd::args2flag(argv,cmds,"--xxx")));
  if(vpflow.flag("XXX")) cout << "XXX" << endl;

  // [OBSOLETE] vpflow.flag("LIB2RAW_ALL",aurostd::args2flag(argv,cmds,"--lib2raw_all|--xraw_all"));
  // [OBSOLETE]  if(vpflow.flag("LIB2RAW_ALL")) vpflow.flag("MULTISH",FALSE);
  // [OBSOLETE] 
  // [OBSOLETE] if(LIB2RAW_ALL) MULTISH=FALSE;
  // [OBSOLETE] ivpflow.flag("LIB2RAW_ALL_FORCE",aurostd::args2flag(argv,cmds,"--lib2raw_all_force|--xraw_all_force"));
  // [OBSOLETE] if(LIB2RAW_ALL_FORCE) MULTISH=FALSE;
  // [OBSOLETE] vpflow.flag("LIB2RAW_FORCE",aurostd::args2flag(argv,cmds,"--lib2raw_force|--xraw_force") || aurostd::args2attachedflag(argv,cmds,"--lib2raw_force=|--xraw_force="));
  // [OBSOLETE] vpflow.flag("LIB2RAW",aurostd::args2flag(argv,cmds,"--lib2raw|--xraw") || aurostd::args2attachedflag(argv,cmds,"--lib2raw=|--xraw="));
  vpflow.args2addattachedscheme(argv,cmds,"LIB2RAW","--lib2raw=","");
  if(vpflow.flag("LIB2RAW")){vpflow.flag("LIB2RAW_LOCAL",aurostd::args2flag(argv,cmds,"--local"));}
  vpflow.flag("FORCE",aurostd::args2flag(argv,cmds,"--force"));

  vpflow.flag("XPLUG",aurostd::args2flag(argv,cmds,"--xplug"));

  // [OBSOLETE] vpflow.flag("XRD_DIST",aurostd::args2flag(argv,cmds,"--xrd_dist"));
  vpflow.args2addattachedscheme(argv,cmds,"XRD_DIST","--xrd_dist=|--XRD_DIST=","");

  // [OBSOLETE] vpflow.flag("ZVAL",aurostd::args2flag(argv,cmds,"--zval|--ZVAL|--zval_cell|--ZVAL_CELL|--pomass|--POMASS|--pomass_cell|--POMASS_CELL"));
  vpflow.args2addattachedscheme(argv,cmds,"ZVAL","--zval=|--ZVAL=","");
  vpflow.args2addattachedscheme(argv,cmds,"ZVAL::CELL","--zval_cell=|--ZVAL_CELL=|--zvalcell=|--ZVALCELL=","");
  vpflow.args2addattachedscheme(argv,cmds,"ZVAL::ATOM","--zval_atom=|--ZVAL_ATOM=|--zvalatom=|--ZVALATOM=","");
  vpflow.args2addattachedscheme(argv,cmds,"POMASS","--pomass=|--POMASS=","");
  vpflow.args2addattachedscheme(argv,cmds,"POMASS::CELL","--pomass_cell=|--POMASS_CELL=","");
  vpflow.args2addattachedscheme(argv,cmds,"POMASS::ATOM","--pomass_atom=|--POMASS_ATOM=|--pomassatom=|--POMASSATOM=","");

  //Richard's symmetry functions (RHT)
  //DX - START
  vpflow.flag("ORTHODEFECT_RHT",aurostd::args2flag(argv,cmds,"--OrthoDefect"));  //RHT
  vpflow.flag("REVERSE_SPACEGROUP_RHT",aurostd::args2flag(argv,cmds,"--revsg")); //RHT
  vpflow.flag("PRIMITIVE_LATTICE_RHT",aurostd::args2flag(argv,cmds, "--primr | --fastprimitivecell | --fprim")); //RHT
  vpflow.args2addattachedscheme(argv,cmds,"WYCCAR_RHT","--wyccar=","");
  // [OBSOLETE]  vpflow.flag("AFLOWSG",aurostd::args2flag(argv,cmds, "--aflowSG")); //RHT  // FIX
  // end Richard's symmetry (RHT)
  // WE MIGHT NEED TO PUT THEM AROUND IN ALPHABETIC ORDER, keep the //RHT
  // DX - END

  // *************************************
  // cluster expansion method
  // [OBSOLETE] vpflow.flag("CE::CLUSTEREXPANSION",aurostd::args2flag(argv,cmds,"--cluster-expansion|--ce") && (argv.size() == 7));
  vpflow.args2addattachedscheme(argv,cmds,"CE::CLUSTEREXPANSION","--cluster-expansion=|--ce=","");
  
  // special Quasirandom Structure (SQS)
  // [OBSOLETE] vpflow.flag("CE::SQS",aurostd::args2flag(argv,cmds,"--special-quasirandom-structure|--sqs") && ((argv.size() == 9) || (argv.size() == 5)));
  vpflow.args2addattachedscheme(argv,cmds,"CE::SQS","--special-quasirandom-structure=|--sqs=",""); 
  // get all clusters
  // [OBSOLETE] vpflow.flag("CE::CLUSTERS",aurostd::args2flag(argv,cmds,"--cluster") && (argv.size() == 7));
  vpflow.args2addattachedscheme(argv,cmds,"CE::CLUSTERS","--cluster=|--clusters=","");
  // get all superlattices
  // [OBSOLETE] vpflow.flag("CE::SUPERLATTICE",aurostd::args2flag(argv,cmds,"--superlattice=") && (argv.size() == 5));
  vpflow.args2addattachedscheme(argv,cmds,"CE::SUPERLATTICE","--superlattice=","");
  
  // *************************************
  // effective mass
  vpflow.flag("EFFECTIVEMASS",aurostd::args2flag(argv,cmds,"--effective-mass|--em"));// && (argv.size() == 3));
  
  if(LDEBUG) cout << "PflowARGs: scheme=" << vpflow.scheme << endl;
  if(LDEBUG) cout << "PflowARGs: vscheme.size()=" << vpflow.vscheme.size() << endl;
  if(LDEBUG) cout << "PflowARGs: argv.size()=" << argv.size() << endl;

  return vpflow.vscheme.size();
}

namespace pflow {
  int main(vector<string> &argv,vector<string> &cmds) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "pflow::main: BEGIN" << endl;
    // cerr << "C username=" << XHOST.User << endl;
    // cerr << "C groupname=" << XHOST.Group << endl;
    // cerr << "C XHOST.ostrPID=" << XHOST.ostrPID.str() << endl;
    // cerr.flush();
    // cerr << argv.size() << endl;
    aurostd::xoption vpflow;
    ostringstream aus;
    // GENERAL **************************************************
    string progname=argv.at(0);
    //  std::vector<string> cmds;
    // INFORMATION **************************************************
  
 
    // cerr << "************************************************************" << endl;
    // cerr << "* AFLOW IN mode ACONVASP                                   *" << endl;
    // cerr << "************************************************************" << endl;

    _aflags aflags;
    aflags.Directory="./";
    aurostd::args2flag(argv,cmds,"--np="); // put them in cmds
    aurostd::args2flag(argv,cmds,"--npmax"); // put them in cmds
  
    // [OBSOLETE]  PflowARGs(argv,cmds,XHOST.vflag_pflow);         // inside init::InitMachine
    vpflow=XHOST.vflag_pflow;

    if(vpflow.flag("PFLOW_HELP") && argv.size() == 2) {
      cout << "**************************************************************************************************" << endl;
      cout << "  aflow --readme=pflow | --readme=processor | --readme=aconvasp | --readme_aconvasp" << endl;
      cout << "     Returns the HELP information for the \"processing machinery\""  << endl;
      //cout<<AFLOW_AConvaspHelp();
      cout<<aflow::Banner("BANNER_BIG");exit(1);}
    if(vpflow.flag("PFLOW_HELP") && argv.size() == 3) helpIndividualOption(argv);
    if(vpflow.flag("PROTOS")) {cout<<aflowlib::PrototypesHelp();cout<<aflow::Banner("BANNER_BIG");exit(1);}

    if(vpflow.flag("FIX_BANDS")) {pflow::FIXBANDS(aflags,vpflow.getattachedscheme("FIX_BANDS"));exit(0);} 
  
    string EXTRACT_KPOINTS=aurostd::args2string(argv,cmds,"--extract_kpoints|--xkpoints","nan");
    string EXTRACT_INCAR=aurostd::args2string(argv,cmds,"--extract_incar|--xincar","nan");
    string EXTRACT_POSCAR=aurostd::args2string(argv,cmds,"--extract_poscar|--xposcar","nan");
    string EXTRACT_POTCAR=aurostd::args2string(argv,cmds,"--extract_potcar|--xpotcar","nan");

    // aflow style operations
    // [OBSOLETE] aflags.AFLOW_PERFORM_CLEAN=aurostd::args2flag(argv,cmds,"--CLEAN|--clean");
    aflags.AFLOW_PERFORM_CLEAN=XHOST.vflag_aflow.flag("CLEAN");
    vpflow.flag("CLEAN",aflags.AFLOW_PERFORM_CLEAN);
    // [OBSOLETE] aflags.AFLOW_PERFORM_DIRECTORY=aurostd::args2flag(argv,cmds,"--DIRECTORY|--D|--d|-D|-d");
    // [OBSOLETE] aflags.AFLOW_PERFORM_FILE=aurostd::args2flag(argv,cmds,"--F|--FILE|--f");
    aflags.AFLOW_PERFORM_DIRECTORY=XHOST.vflag_control.flag("DIRECTORY");
    aflags.AFLOW_PERFORM_FILE=XHOST.vflag_control.flag("FILE");
    if(vpflow.flag("CLEAN")) {
      if(!aflags.AFLOW_PERFORM_DIRECTORY) {
	cerr << "AFLOW: to use --clean, you must specify the directory" << endl; exit(0);
      } else {
	KBIN::Clean(aflags);
      }
    }

    //*************************************
    // DEBUG
    vpflow.args2addattachedscheme(argv,cmds,"QHULL","--qhull=","");
    if(vpflow.flag("QHULL")) {aflowlib::ALIBRARIES(vpflow.getattachedscheme("QHULL"));exit(0);}
    vpflow.args2addattachedscheme(argv,cmds,"EF","--ef=","");
    if(vpflow.flag("EF")) {aflowlib::LIBS_EFormation(vpflow.getattachedscheme("EF"));exit(0);}
  
    //if(XXX) pflow::XXX(argv,cin);
    // if(pflow::CheckCommands(argv,cmds)==FALSE) exit(0);
  
    xstructure a;
    // [OBSOLETE] a.iomode=IOVASP_AUTO;
    //  a.iomode=IOAFLOW_AUTO;
    a.iomode=IOVASP_AUTO;
    xvector<double> _emptyv(1);
    xmatrix<double> _emptym(1,1);
  
    // *********************************************************************
    // always
    //  if(APENNSY) Apennsymain(argv,cmds);
  
    //  if(MILLER) pflow::MILLER(argv,cin);
  
    // check whether a function is run or not
    // if not, give an error message
    bool _PROGRAMRUN=false;
  
    // *********************************************************************
    if(argv.size()==1 && !_PROGRAMRUN) {
      cout << aflow::Banner("BANNER_TINY") << endl;
      //    cout << pflow::Intro_pflow("aconvasp");
      cout << pflow::Intro_pflow("aflow");
      cout << aflow::Banner("BANNER_TINY") << endl;
      _PROGRAMRUN=true;
    }
    //corey
    if(argv.size()>=1 && !_PROGRAMRUN) {
      if(vpflow.flag("BADER")) {cout << bader_functions::BaderCalc(vpflow); _PROGRAMRUN=true;}
      if(vpflow.flag("CHGCAR2JVXL")) {cout << pflow::CHGCAR2JVXL(vpflow); _PROGRAMRUN=true;}
      if(vpflow.flag("CHGDIFF")) {cout << pflow::CHGDIFF(vpflow); _PROGRAMRUN=true;}
      if(vpflow.flag("CHGSUM")) {cout << pflow::CHGSUM(vpflow); _PROGRAMRUN=true;}
      if(vpflow.flag("PREPARE_CHGCAR_4_JMOL")) {cout << bader_functions::prepare_CHGCAR_4_Jmol(vpflow); _PROGRAMRUN=true;}
      //DAVID
      //DX AND COREY - START
      if(vpflow.flag("FULLSYMMETRY")) {pflow::CalculateFullSymmetry(cin,vpflow,cout); _PROGRAMRUN=true;}
      //DX AND COREY - END
      if(vpflow.flag("COMPARE_MATERIAL_DIRECTORY")) {cout << pflow::compareStructureDirectory(vpflow); _PROGRAMRUN=true;}
      if(vpflow.flag("COMPARE_STRUCTURE_DIRECTORY")) {cout << pflow::compareStructureDirectory(vpflow); _PROGRAMRUN=true;}
      if(vpflow.flag("COMPARE_MATERIAL")) {cout << pflow::compareStructures(vpflow); _PROGRAMRUN=true;}
      if(vpflow.flag("COMPARE_STRUCTURE")) {cout << pflow::compareStructures(vpflow); _PROGRAMRUN=true;}
      //DAVID
    }
    // *********************************************************************
    if(argv.size()==2 && !_PROGRAMRUN) {
      // put cartesian or fractional
      if(vpflow.flag("BZPLOT")) {LATTICE::BZPLOTDATA("",cin,1); _PROGRAMRUN=true;}
      if(vpflow.flag("BZPLOTDATA")) {LATTICE::BZPLOTDATA("",cin,0); _PROGRAMRUN=true;}
      // [OBSOLETE] if(vpflow.flag("EWALD")) {pflow::EWALD(argv,cin); _PROGRAMRUN=true;}
      if(vpflow.flag("HNFTOL")) {pflow::HNFTOL(argv,cin,cout); _PROGRAMRUN=true;}
      if(vpflow.flag("ICSD_MAKELABEL")) {pflow::ICSD(argv,cin); _PROGRAMRUN=true;}
      if(vpflow.flag("JMOLGIF")) {pflow::JMOLAnimation(cin,argv); _PROGRAMRUN=true;}
      if(vpflow.flag("KPATH")) {pflow::KPATH(cin,aurostd::args2attachedutype<double>(argv,"--grid=",16.0),vpflow.flag("WWW")); _PROGRAMRUN=true;}
      if(vpflow.flag("NANOPARTICLE")) {cout << pflow::NANOPARTICLE(cin,xvector<double>(0)); _PROGRAMRUN=true;}
      if(vpflow.flag("QMVASP")) {pflow::QMVASP(argv); _PROGRAMRUN=true;}
      // [OBSOLETE] if(vpflow.flag("SG::FINDSYM_PRINT")) {pflow::FINDSYM(vpflow.getattachedscheme("SG::FINDSYM_PRINT"),0,cin); _PROGRAMRUN=true;}
      // [OBSOLETE] if(vpflow.flag("SG::FINDSYM_EXEC")) {pflow::FINDSYM(vpflow.getattachedscheme("SG::FINDSYM_EXEC"),1,cin); _PROGRAMRUN=true;}
      // if(vpflow.flag("PROTO_GUS_CPP")) {pflow::PROTO_GUS_CPP(argv); _PROGRAMRUN=true;}
    }
    // *********************************************************************
    if(argv.size()>=2 && !_PROGRAMRUN) {
      // A
      if(vpflow.flag("ABINIT")) {cout << AQEgeom2abinit(cin); _PROGRAMRUN=true;}
      if(vpflow.flag("AIMS")) {cout << AQEgeom2aims(cin); _PROGRAMRUN=true;}
      if(vpflow.flag("ABCCAR")) {cout << pflow::ABCCAR(cin); _PROGRAMRUN=true;}
      if(vpflow.flag("ACE")) {pflow::ACE(cin); _PROGRAMRUN=true;}
      if(vpflow.flag("AFLOWIN")) {cout << pflow::AFLOWIN(cin); _PROGRAMRUN=true;}
      //DX 8/18/17 [OBSOLETE] if(vpflow.flag("AGROUP")) {pflow::AGROUP(aflags,cin); _PROGRAMRUN=true;}
      if(vpflow.flag("AGROUP")) {pflow::SYMMETRY_GROUPS(aflags,cin,vpflow,cout); _PROGRAMRUN=true;} //DX 8/18/17
      if(vpflow.flag("AGROUP2")) {pflow::AGROUP2(cin); _PROGRAMRUN=true;}
      if(vpflow.flag("AGROUP2m")) {pflow::AGROUP2m(cin); _PROGRAMRUN=true;}
      if(vpflow.flag("ALPHABETIC")) {cout << pflow::ALPHABETIC(cin); _PROGRAMRUN=true;}
      if(vpflow.flag("ANGLES")) {pflow::ANGLES(vpflow.getattachedscheme("ANGLES"),cin); _PROGRAMRUN=true;}
      if(vpflow.flag("ALPHA_COMPOUND")) {cout << pflow::ALPHACompound(vpflow.getattachedscheme("ALPHA_COMPOUND")); _PROGRAMRUN=true;}
      if(vpflow.flag("ALPHA_SPECIES")) {cout << pflow::ALPHASpecies(vpflow.getattachedscheme("ALPHA_SPECIES")); _PROGRAMRUN=true;}
      if(vpflow.flag("AFLOWLIB")) {aflowlib::WEB_Aflowlib_Entry_PHP(vpflow.getattachedscheme("AFLOWLIB"),cout); _PROGRAMRUN=true;}
      if(vpflow.flag("AFLOWLIB_AUID2AURL")) {cout << aflowlib::AflowlibLocator(vpflow.getattachedscheme("AFLOWLIB_AUID2AURL"),"AFLOWLIB_AUID2AURL"); _PROGRAMRUN=true;}
      if(vpflow.flag("AFLOWLIB_AURL2AUID")) {cout << aflowlib::AflowlibLocator(vpflow.getattachedscheme("AFLOWLIB_AURL2AUID"),"AFLOWLIB_AURL2AUID"); _PROGRAMRUN=true;}
      if(vpflow.flag("AFLOWLIB_AUID2LOOP")) {cout << aflowlib::AflowlibLocator(vpflow.getattachedscheme("AFLOWLIB_AUID2LOOP"),"AFLOWLIB_AUID2LOOP"); _PROGRAMRUN=true;}
      if(vpflow.flag("AFLOWLIB_AURL2LOOP")) {cout << aflowlib::AflowlibLocator(vpflow.getattachedscheme("AFLOWLIB_AURL2LOOP"),"AFLOWLIB_AURL2LOOP"); _PROGRAMRUN=true;}
      // B
      if(vpflow.flag("BANDGAP_WAHYU")) {AConvaspBandgap(argv); _PROGRAMRUN=true;}
      if(vpflow.flag("BANDGAP"))       { pflow::BANDGAP(vpflow, cout) ; _PROGRAMRUN=true ; } // CAMILO  //CO 171006
      if(vpflow.flag("BZPLOTDATAUSEKPOINTS")) {LATTICE::BZPLOTDATA(vpflow.getattachedscheme("BZPLOTDATAUSEKPOINTS"),cin,10); _PROGRAMRUN=true;}
      if(vpflow.flag("BZPLOTUSEKPOINTS")) {LATTICE::BZPLOTDATA(vpflow.getattachedscheme("BZPLOTUSEKPOINTS"),cin,11); _PROGRAMRUN=true;}
      if(vpflow.flag("BANDSTRUCTURE")) {pflow::BANDSTRUCTURE(aflags); _PROGRAMRUN=true;}
      if(vpflow.flag("BANDGAPS")) {AConvaspBandgaps(cin,cout); _PROGRAMRUN=true;}
      if(vpflow.flag("BANDGAPDOS")) {AConvaspBandgapFromDOS(cin); _PROGRAMRUN=true;}
      if(vpflow.flag("BANDGAPLISTDOS")) {AConvaspBandgapListFromDOS(cin); _PROGRAMRUN=true;}
      if(vpflow.flag("BZDIRECTION")) {
	cerr << "[" << vpflow.getattachedscheme("BZDIRECTION") << "]" << endl;
	if(vpflow.getattachedscheme("BZDIRECTION").empty()) {
	  cout << pflow::BZDirectionsSTRUCTURE(cin);
	} else {
	  cout << pflow::BZDirectionsLATTICE(vpflow.getattachedscheme("BZDIRECTION"));
	}
	_PROGRAMRUN=true;
      }
      if(vpflow.flag("BZMAX")) {pflow::BZMAX(cin); _PROGRAMRUN=true;}
      if(vpflow.flag("BANDS")) {pflow::BANDS(vpflow.getattachedscheme("BANDS"),cin); _PROGRAMRUN=true;}
      // C
      if(vpflow.flag("CAGES") && !AFLOW_PTHREADS::FLAG) {pflow::CAGES(aflags,vpflow.getattachedscheme("CAGES"),cin); _PROGRAMRUN=true;}
      if(vpflow.flag("CAGES") &&  AFLOW_PTHREADS::FLAG) {pflow::CAGES(aflags,vpflow.getattachedscheme("CAGES"),cin); _PROGRAMRUN=true;}
      if(vpflow.flag("CART")) {cout << pflow::CART(cin); _PROGRAMRUN=true;}
      if(vpflow.flag("CHECKINTEGRITIY")) {pflow::CheckIntegritiy(); _PROGRAMRUN=true;}
      if(vpflow.flag("CHANGESUFFIX")) {pflow::ChangeSuffix(vpflow.getattachedscheme("CHANGESUFFIX")); _PROGRAMRUN=true;} //KESONG Dec.22, 2013
      if(vpflow.flag("CIF")) {pflow::CIF(cin); _PROGRAMRUN=true;}
      if(vpflow.flag("CLEANALL")) {pflow::CLEANALL(cin); _PROGRAMRUN=true;}
      if(vpflow.flag("CORNERS")) {cout << pflow::CORNERS(cin); _PROGRAMRUN=true;}
      if(vpflow.flag("CALCULATED_ICSD_RANDOM")) {cout << aflowlib::CALCULATED_ICSD_RANDOM(); _PROGRAMRUN=true; exit(0);}
      if(vpflow.flag("CALCULATED")) {cout << aflowlib::CALCULATED(vpflow.getattachedscheme("CALCULATED")); _PROGRAMRUN=true;}
      if(vpflow.flag("CLAT")) {pflow::CLAT(vpflow.getattachedscheme("CLAT")); _PROGRAMRUN=true;}
      if(vpflow.flag("CHULL::INIT")) {chull::convexHull(vpflow); _PROGRAMRUN=true;}
      if(vpflow.flag("COMPARE")) {pflow::COMPARE(vpflow.getattachedscheme("COMPARE")); _PROGRAMRUN=true;}
      if(vpflow.flag("CE::CLUSTERS")) {pflow::Cluster(vpflow.getattachedscheme("CE::CLUSTERS")); _PROGRAMRUN=true;}
      if(vpflow.flag("CE::CLUSTEREXPANSION")) {pflow::AClusterExpansionMethodMain(vpflow.getattachedscheme("CE::CLUSTEREXPANSION")); _PROGRAMRUN=true;}
      if(vpflow.flag("CE::SUPERLATTICE")) {pflow::Superlattice(vpflow.getattachedscheme("CE::SUPERLATTICE")); _PROGRAMRUN=true;}
      if(vpflow.flag("CE::SQS")) {pflow::SQS(vpflow.getattachedscheme("CE::SQS")); _PROGRAMRUN=true;}
      // D
      //DX 9/1/17 [OBSOLETE] if(vpflow.flag("DATA")) {pflow::DATA("DATA",cin); _PROGRAMRUN=true;}
      if(vpflow.flag("DATA")) {pflow::DATA("DATA",cin,vpflow,cout); _PROGRAMRUN=true;}
      if(vpflow.flag("DATA1")) {pflow::DATA1(vpflow.getattachedscheme("DATA1"),cin); _PROGRAMRUN=true;}
      if(vpflow.flag("DATA2")) {pflow::DATA2(cin); _PROGRAMRUN=true;}
      if(vpflow.flag("DEBYE")) {pflow::DEBYE(vpflow.getattachedscheme("DEBYE")); _PROGRAMRUN=true;}
      if(vpflow.flag("DIFF")) {pocc::DIFF(vpflow.getattachedscheme("DIFF")); _PROGRAMRUN=true;}
      if(vpflow.flag("DISP")) {pflow::DISP(vpflow.getattachedscheme("DISP"),cin); _PROGRAMRUN=true;}
      if(vpflow.flag("DIST")) {pflow::DIST(vpflow.getattachedscheme("DIST"),cin); _PROGRAMRUN=true;}
      //if(DYNADIEL) {pflow::DYNADIEL(argv) ; _PROGRAMRUN=true ;} // CAMILO
      // E
      //DX 9/1/17 [OBSOLETE] if(vpflow.flag("EDATA")) {pflow::DATA("EDATA",cin); _PROGRAMRUN=true;}
      if(vpflow.flag("EDATA")) {pflow::DATA("EDATA",cin,vpflow,cout); _PROGRAMRUN=true;}
      if(vpflow.flag("EDOS")) {pflow::EDOS(argv); _PROGRAMRUN=true;}
      if(vpflow.flag("EFFMASS")) { pflow::EFFMASS(argv, cout) ; _PROGRAMRUN=true ; } // CAMILO
      //  if(vpflow.flag("EFFECTIVEMASS")) {pflow::EffectiveMass(argv,aurostd::args2string(argv,"--em","./"),cout); _PROGRAMRUN=true;}
      if(vpflow.flag("EIGCURV")) {pflow::EIGCURV(vpflow.getattachedscheme("EIGCURV"),cout) ; _PROGRAMRUN=true ;} // CAMILO
      //DX 8/18/17 [OBSOLETE] if(vpflow.flag("EQUIVALENT")) {cout << pflow::EQUIVALENT(aflags,cin); _PROGRAMRUN=true;}
      if(vpflow.flag("EQUIVALENT")) {cout << pflow::EQUIVALENT(aflags,cin,vpflow); _PROGRAMRUN=true;}
      if(vpflow.flag("EWALD")) {pflow::EWALD(vpflow.getattachedscheme("EWALD"),cin); _PROGRAMRUN=true;}
      // F
      if(vpflow.flag("FROZSL_VASPSETUP_AFLOW")) {cout << pflow::FROZSL_VASPSETUP(argv,0); _PROGRAMRUN=true;}
      if(vpflow.flag("FROZSL_VASPSETUP_POSCAR")) {cout << pflow::FROZSL_VASPSETUP(argv,1); _PROGRAMRUN=true;}
      if(vpflow.flag("FROZSL_ANALYZE")) {cout << pflow::FROZSL_ANALYZE(cin); _PROGRAMRUN=true;}
      if(vpflow.flag("FROZSL_README")) {cout << init::InitGlobalObject("README_AFLOW_FROZSL_TXT") << endl; _PROGRAMRUN=true;}
      if(vpflow.flag("FROZSL_INPUT")) {cout << pflow::FROZSL_INPUT(); _PROGRAMRUN=true;}
      if(vpflow.flag("FROZSL_OUTPUT")) {cout << pflow::FROZSL_OUTPUT(); _PROGRAMRUN=true;}
      //DX 8/18/17 [OBSOLETE] if(vpflow.flag("FGROUP")) {pflow::FGROUP(aflags,cin); _PROGRAMRUN=true;}
      if(vpflow.flag("FGROUP")) {pflow::SYMMETRY_GROUPS(aflags,cin,vpflow,cout); _PROGRAMRUN=true;} //DX 8/18/17
      if(vpflow.flag("FRAC")) {cout << pflow::FRAC(cin); _PROGRAMRUN=true;}
      // G
      if(vpflow.flag("GETTEMP")) {AFLOW_getTEMP(argv); _PROGRAMRUN=true;} //cout << "TEMPERATURE " << Message("host") << endl;exit(0); _PROGRAMRUN=true;}
      if(vpflow.flag("GULP")) {pflow::GULP(cin); _PROGRAMRUN=true;}
      // H
      if(vpflow.flag("HKL")) {pflow::HKL(vpflow.getattachedscheme("HKL"),aflags,cin); _PROGRAMRUN=true;}
      if(vpflow.flag("HKL_SEARCH_TRIVIAL")) {pflow::HKLSearch(vpflow.getattachedscheme("HKL_SEARCH_TRIVIAL"),aflags,cin,"HKL_SEARCH_TRIVIAL"); _PROGRAMRUN=true;}
      if(vpflow.flag("HKL_SEARCH_SIMPLE")) {pflow::HKLSearch(vpflow.getattachedscheme("HKL_SEARCH_SIMPLE"),aflags,cin,"HKL_SEARCH_SIMPLE"); _PROGRAMRUN=true;}
      if(vpflow.flag("HKL_SEARCH_COMPLETE")) {pflow::HKLSearch(vpflow.getattachedscheme("HKL_SEARCH_COMPLETE"),aflags,cin,"HKL_SEARCH_COMPLETE"); _PROGRAMRUN=true;}
      if(vpflow.flag("HNFCELL")) {pocc::HNFCELL(cin); _PROGRAMRUN=true;}
      // K
      if(vpflow.flag("KBAND")) {pflow::KBAND(argv); _PROGRAMRUN=true;}
      if(vpflow.flag("KPOINTS")) {pflow::KPOINTS(vpflow.getattachedscheme("KPOINTS"),cin,cout); _PROGRAMRUN=true;}
      if(vpflow.flag("FLAG::XVASP_KPOINTS_DELTA")) {pflow::KPOINTS_DELTA(vpflow,cin,cout); _PROGRAMRUN=true;}
      if(vpflow.flag("KILL")) {sflow::KILL(vpflow.getattachedscheme("KILL")); _PROGRAMRUN=true;} 
      // J
      if(vpflow.flag("JUSTAFTER")) {sflow::JUST(vpflow.getattachedscheme("JUSTAFTER"),cin,"JUSTAFTER"); _PROGRAMRUN=true;}
      if(vpflow.flag("JUSTBEFORE")) {sflow::JUST(vpflow.getattachedscheme("JUSTBEFORE"),cin,"JUSTBEFORE"); _PROGRAMRUN=true;}
      if(vpflow.flag("JUSTBETWEEN")) {sflow::JUST(vpflow.getattachedscheme("JUSTBETWEEN"),cin,"JUSTBETWEEN"); _PROGRAMRUN=true;}
      if(vpflow.flag("JMOL")) {pflow::JMOL(vpflow.getattachedscheme("JMOL"),cin); _PROGRAMRUN=true;}
      // I
      if(vpflow.flag("ICSD") || vpflow.flag("ICSD_CHEM") || vpflow.flag("ICSD_PROTO") || vpflow.flag("ICSD_ID") || vpflow.flag("ICSD_LESSTHAN") || vpflow.flag("ICSD_MORETHAN") ||
	 vpflow.flag("ICSD_DENSLESSTHAN") || vpflow.flag("ICSD_DENSMORETHAN") || vpflow.flag("ICSD_SG") || vpflow.flag("ICSD_SGLESSTHAN") ||
	 vpflow.flag("ICSD_SGMORETHAN") || vpflow.flag("ICSD_TRICLINIC") || vpflow.flag("ICSD_MONOCLINIC") || vpflow.flag("ICSD_ORTHORHOMBIC") ||
	 vpflow.flag("ICSD_TETRAGONAL") || vpflow.flag("ICSD_RHOMBOHEDRAL") || vpflow.flag("ICSD_TRIGONAL") || vpflow.flag("ICSD_HEXAGONAL") || vpflow.flag("ICSD_CUBIC") ||
	 vpflow.flag("ICSD_UNIQUE") || vpflow.flag("ICSD_BASISLT") || vpflow.flag("ICSD_BASISGT") || vpflow.flag("ICSD_NOBROKENBASIS") ||
	 vpflow.flag("ICSD_NOPARTIALOCC") || vpflow.flag("ICSD_N_ARY") ||
	 vpflow.flag("ICSD_REMOVE_AND") || vpflow.flag("ICSD_REMOVE_OR") || vpflow.flag("ICSD_REMOVEMETALS") ||
	 vpflow.flag("ICSD_ALLLESSTHAN") || vpflow.flag("ICSD_ALLMORETHAN") ||
	 vpflow.flag("ICSD_TRI") || vpflow.flag("ICSD_MCL") || vpflow.flag("ICSD_MCLC") ||
	 vpflow.flag("ICSD_ORC") || vpflow.flag("ICSD_ORCC") || vpflow.flag("ICSD_ORCF") || vpflow.flag("ICSD_ORCI") ||
	 vpflow.flag("ICSD_TET") || vpflow.flag("ICSD_BCT") || vpflow.flag("ICSD_RHL") || vpflow.flag("ICSD_HEX") ||
	 vpflow.flag("ICSD_CUB") || vpflow.flag("ICSD_FCC") || vpflow.flag("ICSD_BCC")) {pflow::ICSD(argv,cin); _PROGRAMRUN=true;}
      if(vpflow.flag("ICSD_CHECK_RAW")) {pflow::ICSD_CheckRaw(argv); _PROGRAMRUN=true;}
      if(vpflow.flag("ICSD_LISTMETALS")) {pflow::ICSD_ListMetals(); _PROGRAMRUN=true;}
      if(vpflow.flag("ICSD2POSCAR")) {pflow::ICSD_2POSCAR(cin); _PROGRAMRUN=true;}
      if(vpflow.flag("ICSD2PROTO")) {pflow::ICSD_2PROTO(cin); _PROGRAMRUN=true;}
      if(vpflow.flag("ICSD2WYCK")) {pflow::ICSD_2WYCK(cin,vpflow.flag("SOF")); _PROGRAMRUN=true;}
      if(vpflow.flag("IDENTICAL")) {cout << pflow::IDENTICAL(cin); _PROGRAMRUN=true;}
      if(vpflow.flag("INCELL")) {cout << pflow::INCELL(cin); _PROGRAMRUN=true;}
      if(vpflow.flag("INCOMPACT")) {cout << pflow::INCOMPACT(cin); _PROGRAMRUN=true;}
      if(vpflow.flag("INTPOL")) {pflow::INTPOL(vpflow.getattachedscheme("INTPOL")); _PROGRAMRUN=true;}
      if(vpflow.flag("INWS")) {cout << pflow::INWS(cin); _PROGRAMRUN=true;}
      // [OBSOLETE] if(vpflow.flag("INFLATE_LATTICE")) {cout << pflow::INFLATE_LATTICE(cin,aurostd::args2utype(argv,"--inflate_lattice|--ilattice",1.0)); _PROGRAMRUN=true;}
      // [OBSOLETE] if(vpflow.flag("INFLATE_VOLUME")) {cout << pflow::INFLATE_VOLUME(cin,aurostd::args2utype(argv,"--inflate_volume|--ivolume",1.0)); _PROGRAMRUN=true;}
      if(vpflow.flag("INFLATE_LATTICE")) {cout << pflow::INFLATE_LATTICE(vpflow.getattachedscheme("INFLATE_LATTICE"),cin); _PROGRAMRUN=true;}
      if(vpflow.flag("INFLATE_VOLUME")) {cout << pflow::INFLATE_VOLUME(vpflow.getattachedscheme("INFLATE_VOLUME"),cin); _PROGRAMRUN=true;}
      // L
      if(vpflow.flag("LATTICEREDUCTION")) {cout << pflow::LATTICEREDUCTION(cin); _PROGRAMRUN=true;}
      if(vpflow.flag("LATTICE_TYPE")) {cout << pflow::LATTICE_TYPE(cin); _PROGRAMRUN=true;}
      if(vpflow.flag("LATTICE_LATTICE_TYPE")) {cout << pflow::LATTICE_LATTICE_TYPE(cin); _PROGRAMRUN=true;}
      if(vpflow.flag("LATTICE_HISTOGRAM")) {CheckLatticeHistogram(); _PROGRAMRUN=true;}
      if(vpflow.flag("LIB2RAW")) {XHOST.sensors_allowed=FALSE;aflowlib::LIB2RAW(vpflow.getattachedscheme("LIB2RAW"),vpflow.flag("FORCE"),vpflow.flag("LIB2RAW_LOCAL"));XHOST.sensors_allowed=TRUE; _PROGRAMRUN=true;}
      if(vpflow.flag("LTCELL")) {cout << pflow::LTCELL(vpflow.getattachedscheme("LTCELL"),cin); _PROGRAMRUN=true;}
      if(vpflow.flag("LTCELLFV")) {cout << pflow::LTCELL(vpflow.getattachedscheme("LTCELLFV"),cin); _PROGRAMRUN=true;}
      // M
      // [OBSOLETE] if(vpflow.flag("MILLER")) {cout << pflow::MILLER(vpflow.getattachedscheme("MILLER"),cin); _PROGRAMRUN=true;}
      if(vpflow.flag("MULTIBZIP2")) {AFLOW_PTHREADS::MULTI_bzip2(argv);_PROGRAMRUN=true;}
      if(vpflow.flag("MULTIBUNZIP2")) {AFLOW_PTHREADS::MULTI_bunzip2(argv);_PROGRAMRUN=true;}
      if(vpflow.flag("MULTIGZIP")) {AFLOW_PTHREADS::MULTI_gzip(argv);_PROGRAMRUN=true;}
      if(vpflow.flag("MULTIGUNZIP")) {AFLOW_PTHREADS::MULTI_gunzip(argv);_PROGRAMRUN=true;}
      if(vpflow.flag("MULTISH")) {AFLOW_PTHREADS::MULTI_sh(argv);_PROGRAMRUN=true;}
      if(vpflow.flag("MULTIZIP")) {AFLOW_PTHREADS::MULTI_zip(argv);_PROGRAMRUN=true;}
      if(vpflow.flag("MAGNETICPARAMETERS")) {pflow::MagneticParameters(aurostd::args2attachedstring(argv,"--magpara=","./"),cout); _PROGRAMRUN=true;}
      if(vpflow.flag("MINKOWSKI_BASIS_REDUCTION")) {cout << pflow::MINKOWSKIBASISREDUCTION(cin); _PROGRAMRUN=true;}
      if(vpflow.flag("MSI")) {pflow::MSI(cin); _PROGRAMRUN=true;}
      if(vpflow.flag("MOM")) {pflow::MOM(cin); _PROGRAMRUN=true;}
      if(vpflow.flag("MULTIENUMALL")) {pocc::MultienumPrintAllXstr(cin); _PROGRAMRUN=true;}
      if(vpflow.flag("MULTIENUMSORT")) {pocc::MultienumPrintSortedXstr(cin); _PROGRAMRUN=true;}
      if(vpflow.flag("MAXATOMS")) {cout << pflow::ATOMSMAX(vpflow.getattachedscheme("MAXATOMS"),cin); _PROGRAMRUN=true;}
      // N   
      if(vpflow.flag("NAMES") && argv.size()>=2) {cout << pflow::NAMES(argv,cin); _PROGRAMRUN=true;}
      if(vpflow.flag("NUMNAMES") && argv.size()>=2) {cout << pflow::NUMNAMES(argv,cin); _PROGRAMRUN=true;}
      if(vpflow.flag("NATOMS")) {cout << pflow::NATOMS(cin) << endl; _PROGRAMRUN=true;}
      if(vpflow.flag("NBONDXX")) {cout << pflow::NBONDXX(cin); _PROGRAMRUN=true;} //CO 171025
      if(vpflow.flag("NDATA")) {pflow::NDATA(cin); _PROGRAMRUN=true;}
      if(vpflow.flag("NIGGLI")) {cout << pflow::NIGGLI(cin); _PROGRAMRUN=true;}
      if(vpflow.flag("NNDIST")) {cout << pflow::NNDIST(cin) << endl; _PROGRAMRUN=true;}
      if(vpflow.flag("NOORDERPARAMETER")) {cout << pflow::NOORDERPARAMETER(cin); _PROGRAMRUN=true;}
      if(vpflow.flag("NSPECIES")) {cout << pflow::NSPECIES(cin) << endl; _PROGRAMRUN=true;}
      // O
      // if(vpflow.flag("OPARAMETER")) {pflow::OPARAMETER(argv,cin); _PROGRAMRUN=true;}
      // P
      
      // Serializers for DOS and bands
      if(vpflow.flag("DOSDATA2JSON")) {estructure::DOSDATA_JSON(vpflow,cout); _PROGRAMRUN=true;} //Eric G
      if(vpflow.flag("BANDSDATA2JSON")) {estructure::BANDSDATA_JSON(vpflow,cout); _PROGRAMRUN=true;} //Eric G
      // End serializers

      if(vpflow.flag("PLOT_BAND")) {estructure::PLOT_BAND(vpflow.getattachedscheme("PLOT_BAND")); _PROGRAMRUN=true;} 
      if(vpflow.flag("PLOT_BANDSPINSPLIT")) {estructure::PLOT_BAND_SPINSPLIT(vpflow.getattachedscheme("PLOT_BANDSPINSPLIT")); _PROGRAMRUN=true;} 
      if(vpflow.flag("PLOT_BAND2")) {estructure::PLOT_BAND2(vpflow.getattachedscheme("PLOT_BAND2")); _PROGRAMRUN=true;} 
      if(vpflow.flag("PLOT_BANDDOS")) {estructure::PLOT_BANDDOS(vpflow.getattachedscheme("PLOT_BANDDOS")); _PROGRAMRUN=true;} 
      if(vpflow.flag("PLOT_DOS")) {estructure::PLOT_DOS(vpflow.getattachedscheme("PLOT_DOS")); _PROGRAMRUN=true;} 
      if(vpflow.flag("PLOT_DOSWEB")) {estructure::PLOT_DOSWEB(vpflow.getattachedscheme("PLOT_DOSWEB")); _PROGRAMRUN=true;} 
      if(vpflow.flag("PLOT_PEDOS")) {estructure::PLOT_PEDOS(vpflow.getattachedscheme("PLOT_PEDOS")); _PROGRAMRUN=true;} 
      if(vpflow.flag("PLOT_PEDOSALL")) {estructure::PLOT_PEDOSALL(vpflow.getattachedscheme("PLOT_PEDOSALL")); _PROGRAMRUN=true;} 
      if(vpflow.flag("PLOT_PEDOSALL_AFLOWLIB")) {estructure::PLOT_PEDOSALL_AFLOWLIB(vpflow.getattachedscheme("PLOT_PEDOSALL_AFLOWLIB"), aflags); _PROGRAMRUN=true;} 
      if(vpflow.flag("PLOTPHDISP")) {pflow::PLOT_PHDISP(argv); _PROGRAMRUN=true;}
      if(vpflow.flag("PROTOS_ICSD")) {cout<<aflowlib::PrototypesIcsdHelp(vpflow.getattachedscheme("PROTOS_ICSD"));cout<<aflow::Banner("BANNER_BIG");exit(1);}
      // if(POCCUPATION) {pflow::POCCUPATION(argv,cin); _PROGRAMRUN=true;}
      if(vpflow.flag("POCC_DOS")) {pocc::POCC_DOS(cout,vpflow.getattachedscheme("POCC_DOS")); _PROGRAMRUN=true;} 
      if(vpflow.flag("POCC_MAG")) {pocc::POCC_MAG(vpflow.getattachedscheme("POCC_MAG")); _PROGRAMRUN=true;} 
      if(vpflow.flag("POCC_BANDGAP")) {pocc::POCC_BANDGAP(vpflow.getattachedscheme("POCC_BANDGAP")); _PROGRAMRUN=true;} 
      // 
      if(vpflow.flag("PROTO")) {cout << pflow::PROTO_LIBRARIES(vpflow); _PROGRAMRUN=true;}
      if(vpflow.flag("PROTO_AFLOW")) {pflow::PROTO_AFLOW(vpflow,FALSE); _PROGRAMRUN=true;} // non reversed
      if(vpflow.flag("PLATON") && argv.size()>=2) {cout << pflow::PLATON(vpflow.getattachedscheme("PLATON"),cin); _PROGRAMRUN=true;}// << endl;
      if(vpflow.flag("POMASS")) {pflow::ZVAL("POMASS,"+vpflow.getattachedscheme("POMASS")); _PROGRAMRUN=true;}
      if(vpflow.flag("POMASS::CELL")) {pflow::ZVAL("POMASS::CELL,"+vpflow.getattachedscheme("POMASS::CELL")); _PROGRAMRUN=true;}
      if(vpflow.flag("POMASS::ATOM")) {pflow::ZVAL("POMASS::ATOM,"+vpflow.getattachedscheme("POMASS::ATOM")); _PROGRAMRUN=true;}
      if(vpflow.flag("PEARSON_SYMBOL")) {cout << pflow::PEARSON_SYMBOL(cin); _PROGRAMRUN=true;}
      if(vpflow.flag("PDB")) {pflow::PDB(cin); _PROGRAMRUN=true;}
      //DX 8/18/17 [OBSOLETE] if(vpflow.flag("PGROUP")) {pflow::PGROUP(aflags,cin); _PROGRAMRUN=true;}
      if(vpflow.flag("PGROUP")) {pflow::SYMMETRY_GROUPS(aflags,cin,vpflow,cout); _PROGRAMRUN=true;} //DX 8/18/17
      //DX 8/18/17 [OBSOLETE] if(vpflow.flag("PGROUPX")) {pflow::PGROUPXTAL(aflags,cin); _PROGRAMRUN=true;}
      if(vpflow.flag("PGROUPX")) {pflow::SYMMETRY_GROUPS(aflags,cin,vpflow,cout); _PROGRAMRUN=true;} //DX 8/18/17
      //DX 8/18/17 [OBSOLETE] if(vpflow.flag("PGROUPK")) {pflow::PGROUPK(aflags,cin); _PROGRAMRUN=true;}
      if(vpflow.flag("PGROUPK")) {pflow::SYMMETRY_GROUPS(aflags,cin,vpflow,cout); _PROGRAMRUN=true;} //DX 8/18/17
      if(vpflow.flag("PGROUPK_XTAL")) {pflow::SYMMETRY_GROUPS(aflags,cin,vpflow,cout); _PROGRAMRUN=true;} //DX 12/5/17
      if(vpflow.flag("POSCAR")) {cout << pflow::POSCAR(cin); _PROGRAMRUN=true;}
      if(vpflow.flag("POSCAR2AFLOWIN")) {cout << pflow::POSCAR2AFLOWIN(cin); _PROGRAMRUN=true;}
      if(vpflow.flag("POSCAR2ENUM")) {pocc::POSCAR2ENUM(cin); _PROGRAMRUN=true;}
      if(vpflow.flag("POSCAR2GULP")) {pocc::POSCAR2GULP(cin); _PROGRAMRUN=true;}
      if(vpflow.flag("POCC_INPUT")) {pflow::POCC_INPUT(); _PROGRAMRUN=true;}
      if(vpflow.flag("POSCAR2WYCKOFF")) {pflow::POSCAR2WYCKOFF(cin); _PROGRAMRUN=true;}
      if(vpflow.flag("PRIM")) {cout << pflow::PRIM(cin,0); _PROGRAMRUN=true;}
      if(vpflow.flag("PRIM1")) {cout << pflow::PRIM(cin,1); _PROGRAMRUN=true;}
      if(vpflow.flag("PRIM2")) {cout << pflow::PRIM(cin,2); _PROGRAMRUN=true;}
      if(vpflow.flag("PRIM3")) {cout << pflow::PRIM(cin,3); _PROGRAMRUN=true;}
      // Q
      if(vpflow.flag("QE")) {cout << AQEgeom2qe(cin); _PROGRAMRUN=true;}
      if(vpflow.flag("QDEL")) {sflow::QDEL(vpflow.getattachedscheme("QDEL")); _PROGRAMRUN=true;} // NEW
      if(vpflow.flag("QSUB")) {sflow::QSUB(vpflow.getattachedscheme("QSUB")); _PROGRAMRUN=true;} // NEW
      // R
      if(vpflow.flag("RDF")) {pflow::RDF(vpflow.getattachedscheme("RDF"),cin); _PROGRAMRUN=true;}
      if(vpflow.flag("RDFCMP")) {pflow::RDFCMP(vpflow.getattachedscheme("RDFCMP")); _PROGRAMRUN=true;}
      if(vpflow.flag("RMCOPIES")) {cout << pflow::RMCOPIES(cin); _PROGRAMRUN=true;}
      if(vpflow.flag("RSM")) {pflow::RSM(argv,cin); _PROGRAMRUN=true;}
      if(vpflow.flag("RASMOL")) {pflow::RASMOL(vpflow.getattachedscheme("RASMOL"),cin); _PROGRAMRUN=true;}
      if(vpflow.flag("RMATOM")) {cout << pflow::RMATOM(cin,aurostd::args2utype(argv,"--rm_atom",(int) (0))); _PROGRAMRUN=true;}
      // S
      if(vpflow.flag("SHELL")) {pflow::SHELL(vpflow.getattachedscheme("SHELL"),cin); _PROGRAMRUN=true;}
      if(vpflow.flag("SHIFT")) {cout << pflow::SHIFT(vpflow.getattachedscheme("SHIFT"),cin); _PROGRAMRUN=true;}
      if(vpflow.flag("SLAB")) {cout << slab::MAKE_SLAB(vpflow.getattachedscheme("SLAB"),cin); _PROGRAMRUN=true;}
      if(vpflow.flag("STATDIEL")) {pflow::STATDIEL(argv) ; _PROGRAMRUN=true ;} // CAMILO
      //DX 9/21/17 [OBSOLETE] if(vpflow.flag("SG::AFLOW") && argv.size()>=2) {cout << pflow::SG(vpflow.getattachedscheme("SG::AFLOW"),cin,"AFLOW","ALL") << endl; _PROGRAMRUN=true;}
      //DX 9/21/17 [OBSOLETE] if(vpflow.flag("SG::AFLOW_LABEL") && argv.size()>=2) {cout << pflow::SG(vpflow.getattachedscheme("SG::AFLOW_LABEL"),cin,"AFLOW","LABEL") << endl; _PROGRAMRUN=true;}
      //DX 9/21/17 [OBSOLETE] if(vpflow.flag("SG::AFLOW_NUMBER") && argv.size()>=2) {cout << pflow::SG(vpflow.getattachedscheme("SG::AFLOW_NUMBER"),cin,"AFLOW","NUMBER") << endl; _PROGRAMRUN=true;}
      //DX 9/21/17 [OBSOLETE] if(vpflow.flag("SG::PLATON") && argv.size()>=2) {cout << pflow::SG(vpflow.getattachedscheme("SG::PLATON"),cin,"PLATON","ALL") << endl; _PROGRAMRUN=true;}
      //DX 9/21/17 [OBSOLETE] if(vpflow.flag("SG::PLATON_LABEL") && argv.size()>=2) {cout << pflow::SG(vpflow.getattachedscheme("SG::PLATON_LABEL"),cin,"PLATON","LABEL") << endl; _PROGRAMRUN=true;}
      //DX 9/21/17 [OBSOLETE] if(vpflow.flag("SG::PLATON_NUMBER") && argv.size()>=2) {cout << pflow::SG(vpflow.getattachedscheme("SG::PLATON_NUMBER"),cin,"PLATON","NUMBER") << endl; _PROGRAMRUN=true;}
      //DX 9/21/17 [OBSOLETE] if(vpflow.flag("SG::FINDSYM") && argv.size()>=2) {cout << pflow::SG(vpflow.getattachedscheme("SG::FINDSYM"),cin,"FINDSYM","ALL") << endl; _PROGRAMRUN=true;}
      //DX 9/21/17 [OBSOLETE] if(vpflow.flag("SG::FINDSYM_LABEL") && argv.size()>=2) {cout << pflow::SG(vpflow.getattachedscheme("SG::FINDSYM_LABEL"),cin,"FINDSYM","LABEL") << endl; _PROGRAMRUN=true;}
      //DX 9/21/17 [OBSOLETE] if(vpflow.flag("SG::FINDSYM_NUMBER") && argv.size()>=2) {cout << pflow::SG(vpflow.getattachedscheme("SG::FINDSYM_NUMBER"),cin,"FINDSYM","NUMBER") << endl; _PROGRAMRUN=true;}

      //DX 9/21/17 [OBSOLETE] if(vpflow.flag("SG::FINDSYM_PRINT")) {pflow::FINDSYM(vpflow.getattachedscheme("SG::FINDSYM_PRINT"),0,cin); _PROGRAMRUN=true;}
      //DX 9/21/17 [OBSOLETE] if(vpflow.flag("SG::FINDSYM_EXEC")) {pflow::FINDSYM(vpflow.getattachedscheme("SG::FINDSYM_EXEC"),1,cin); _PROGRAMRUN=true;}
      
      if(vpflow.flag("SG::AFLOW") && argv.size()>=2) {cout << pflow::SG(vpflow,cin,"AFLOW","ALL") << endl; _PROGRAMRUN=true;} //DX 9/26/17
      if(vpflow.flag("SG::AFLOW_LABEL") && argv.size()>=2) {cout << pflow::SG(vpflow,cin,"AFLOW","LABEL") << endl; _PROGRAMRUN=true;} //DX 9/26/17
      if(vpflow.flag("SG::AFLOW_NUMBER") && argv.size()>=2) {cout << pflow::SG(vpflow,cin,"AFLOW","NUMBER") << endl; _PROGRAMRUN=true;} //DX 9/26/17
      if(vpflow.flag("SG::PLATON") && argv.size()>=2) {cout << pflow::SG(vpflow,cin,"PLATON","ALL") << endl; _PROGRAMRUN=true;} //DX 9/26/17
      if(vpflow.flag("SG::PLATON_LABEL") && argv.size()>=2) {cout << pflow::SG(vpflow,cin,"PLATON","LABEL") << endl; _PROGRAMRUN=true;} //DX 9/26/17
      if(vpflow.flag("SG::PLATON_NUMBER") && argv.size()>=2) {cout << pflow::SG(vpflow,cin,"PLATON","NUMBER") << endl; _PROGRAMRUN=true;} //DX 9/26/17
      if(vpflow.flag("SG::FINDSYM") && argv.size()>=2) {cout << pflow::SG(vpflow,cin,"FINDSYM","ALL") << endl; _PROGRAMRUN=true;} //DX 9/26/17
      if(vpflow.flag("SG::FINDSYM_LABEL") && argv.size()>=2) {cout << pflow::SG(vpflow,cin,"FINDSYM","LABEL") << endl; _PROGRAMRUN=true;} //DX 9/26/17
      if(vpflow.flag("SG::FINDSYM_NUMBER") && argv.size()>=2) {cout << pflow::SG(vpflow,cin,"FINDSYM","NUMBER") << endl; _PROGRAMRUN=true;} //DX 9/26/17

      if(vpflow.flag("SG::FINDSYM_PRINT")) {pflow::FINDSYM(vpflow,0,cin); _PROGRAMRUN=true;} //DX 9/26/17
      if(vpflow.flag("SG::FINDSYM_EXEC")) {pflow::FINDSYM(vpflow,1,cin); _PROGRAMRUN=true;} //DX 9/26/17
    
 
      if(vpflow.flag("SUPERCELL")) {cout << pflow::SUPERCELL(vpflow.getattachedscheme("SUPERCELL"),cin); _PROGRAMRUN=true;}
      if(vpflow.flag("SUPERCELLSTRLIST")) {pflow::SUPERCELLSTRLIST(vpflow.getattachedscheme("SUPERCELLSTRLIST")); _PROGRAMRUN=true;}
      if(vpflow.flag("SD") && argv.size()>=2) {cout << pflow::SD(argv,cin); _PROGRAMRUN=true;}
      if(vpflow.flag("SG")) {pflow::SG(cin); _PROGRAMRUN=true;}
      if(vpflow.flag("SGDATA")) {pflow::SGDATA(cin,vpflow,cout); _PROGRAMRUN=true;} //DX 8/18/17
      //DX 8/18/17 [OBSOLETE] if(vpflow.flag("SGROUP")) {pflow::SGROUP(aflags,cin,KBIN_SYMMETRY_SGROUP_RADIUS_DEFAULT); _PROGRAMRUN=true;}
      if(vpflow.flag("SGROUP")) {pflow::SYMMETRY_GROUPS(aflags,cin,vpflow,cout); _PROGRAMRUN=true;} //DX 8/18/17
      if(vpflow.flag("SPECIES")) {cout << pflow::SPECIES(cin); _PROGRAMRUN=true;}
      if(vpflow.flag("STDCONVCELL")) {cout << GetStandardConventional(xstructure(cin,IOAFLOW_AUTO)); _PROGRAMRUN=true;}
      if(vpflow.flag("STDPRIMCELL")) {cout << GetStandardPrimitive(xstructure(cin,IOAFLOW_AUTO)); _PROGRAMRUN=true;}
      if(vpflow.flag("SCALE")) {cout << pflow::SCALE(vpflow.getattachedscheme("SCALE"),cin); _PROGRAMRUN=true;}

      // T
      // [OBSOLETE] [JUNKAI] if(vpflow.flag("TERDATA")) {INPUTDATAFORTERPHASE(argv); _PROGRAMRUN=true;}
      // [OBSOLETE] [JUNKAI] if(vpflow.flag("TERDATA_EXIST")) {GENERATESTABLELIST(argv); _PROGRAMRUN=true;}
      // U
      if(vpflow.flag("UFFENERGY")) {pocc::UFFENERGY(cin); _PROGRAMRUN=true;}
      // V
      if(vpflow.flag("VASP")) {cout << AQEgeom2vasp(cin); _PROGRAMRUN=true;}
      if(vpflow.flag("VOLUME::EQUAL")) {cout << pflow::VOLUME("VOLUME::EQUAL,"+vpflow.getattachedscheme("VOLUME::EQUAL"),cin); _PROGRAMRUN=true;} 
      if(vpflow.flag("VOLUME::MULTIPLY_EQUAL")) {cout << pflow::VOLUME("VOLUME::MULTIPLY_EQUAL,"+vpflow.getattachedscheme("VOLUME::MULTIPLY_EQUAL"),cin); _PROGRAMRUN=true;} 
      if(vpflow.flag("VOLUME::PLUS_EQUAL")) {cout << pflow::VOLUME("VOLUME::PLUS_EQUAL,"+vpflow.getattachedscheme("VOLUME::PLUS_EQUAL"),cin); _PROGRAMRUN=true;} 
      // X
      if(vpflow.flag("XYZ")) {pflow::XYZ(vpflow.getattachedscheme("XYZ"),cin); _PROGRAMRUN=true;}
      if(vpflow.flag("XYZWS")) {pflow::XYZWS(cin); _PROGRAMRUN=true;}
      if(vpflow.flag("XRAY")) {pflow::XRAY(vpflow.getattachedscheme("XRAY"),cin); _PROGRAMRUN=true;}
      if(vpflow.flag("XRD_DIST")) {pflow::GetAtomicPlaneDist(vpflow.getattachedscheme("XRD_DIST"),cin); _PROGRAMRUN=true;}
      // Y
      // Z
      if(vpflow.flag("ZVAL")) {pflow::ZVAL("ZVAL,"+vpflow.getattachedscheme("ZVAL")); _PROGRAMRUN=true;}
      if(vpflow.flag("ZVAL::CELL")) {pflow::ZVAL("ZVAL::CELL,"+vpflow.getattachedscheme("ZVAL::CELL")); _PROGRAMRUN=true;}
      if(vpflow.flag("ZVAL::ATOM")) {pflow::ZVAL("ZVAL::ATOM,"+vpflow.getattachedscheme("ZVAL::ATOM")); _PROGRAMRUN=true;}
    
      //Richard's Functions:
      if(vpflow.flag("REVERSE_SPACEGROUP_RHT")) {cout << SYM::ReverseSpaceGroup(argv) << endl; _PROGRAMRUN=true;}
      if(vpflow.flag("ORTHODEFECT_RHT")) {SYM::OrthoDefect(cin); _PROGRAMRUN=true;}
      if(vpflow.flag("PRIMITIVE_LATTICE_RHT")) {
        //[OBSOLETE] SetTolerance(argv);
        xstructure str(cin);
        str.GetPrimitiveCell(); cout << str << endl;_PROGRAMRUN=true;}
      //SPACEGROUP FUNCTIONS (AUTO AND MANUAL) (SEE README)
      //[OBSOLETE] (DX) if(vpflow.flag("SPACEGROUP_RHT")) {
      //[OBSOLETE] (DX)  //[OBSOLETE] initsymmats();initglides();initsymops();SetTolerance(argv); 
      //[OBSOLETE] (DX)  xstructure str(cin); double tolerance=SYM::defaultTolerance(str); if(argv.size()==2) {tolerance=SYM::defaultTolerance(str);}; if(argv.size()>=3) {tolerance=atof(argv.at(2).c_str());}; uint sgroup=str.SpaceGroup_ITC(tolerance); cout << GetSpaceGroupName(sgroup) << " " << aurostd::utype2string(sgroup) << endl; _PROGRAMRUN=true;}
      //[OBSOLETE] (DX)if(vpflow.flag("SPACEGROUP_MANUAL_RHT")) {
      //[OBSOLETE] (DX)  //[OBSOLETE] initsymmats();initglides();initsymops();SetTolerance(argv);
      //[OBSOLETE] (DX)  xstructure str(cin); double tolerance=SYM::defaultTolerance(str); int man_int=-1; if(argv.size()==2) {tolerance=SYM::defaultTolerance(str);}; if(argv.size()>=3) {tolerance=atof(argv.at(2).c_str()); man_int=atof(argv.at(3).c_str());}; uint sgroup=str.SpaceGroup_ITC(tolerance,man_int); cout << GetSpaceGroupName(sgroup) << " " << aurostd::utype2string(sgroup) << endl; _PROGRAMRUN=true;}
      //RETURNS SPACE GROUP
      // [OBSOLETE]   if(vpflow.flag("SG::AFLOW") || vpflow.flag("SG::AFLOW_LABEL") || vpflow.flag("SG::AFLOW_NUMBER")) {
      // [OBSOLETE]   initsymmats();initglides();initsymops();SetTolerance(argv);
      // [OBSOLETE]   xstructure str(cin);
      // [OBSOLETE]   if(vpflow.flag("SG::AFLOW")) cout << GetSpaceGroupName(str.SpaceGroup_ITC(false,argv)) << " #" << str.SpaceGroup_ITC(false,argv) << endl;
      // [OBSOLETE]   if(vpflow.flag("SG::AFLOW_LABEL")) cout << GetSpaceGroupName(str.SpaceGroup_ITC(false,argv)) << endl;
      // [OBSOLETE]   if(vpflow.flag("SG::AFLOW_NUMBER")) cout << str.SpaceGroup_ITC(false,argv) << endl;
      // [OBSOLETE]   _PROGRAMRUN=true;}
      //TO CHANGE THE FORMAT OF WYCCAR GO TO SPACEGROUP_ITC FUNCTION AND CHANGE WSS
      aflags.QUIET=TRUE;
      // [OBSOLETE]       bool WRITE=TRUE;
      ofstream File("/dev/null");
      bool verbose=TRUE;
      // DX - START 
      //WYCCAR FUNCTIONS (AUTO and MANUAL)
      // [OBSOLETE]      if(vpflow.flag("WYCCAR_RHT")) {initsymmats();initglides();initsymops();SetTolerance(argv);
      // [OBSOLETE]   xstructure str(cin); str.SpaceGroup_ITC(false,argv);printWyccar(File,str,aflags,WRITE,verbose,cout); _PROGRAMRUN=true;}
      // [OBSOLETE]       if(vpflow.flag("WYCCAR_MANUAL_RHT")) {initsymmats();initglides();initsymops();SetTolerance(argv);
      // [OBSOLETE]   xstructure str(cin); str.SpaceGroup_ITC(true,argv);printWyccar(File,str,aflags,WRITE,verbose,cout); _PROGRAMRUN=true;}

      if(vpflow.flag("WYCCAR_RHT")) {// [OBSOLETE] initsymmats();initglides();initsymops();SetTolerance(argv);
        xstructure str(cin); 
        str.ReScale(1.0);  //170804 DX - rescale from input
        string options = vpflow.getattachedscheme("WYCCAR_RHT"); 
        vector<string> tokens;
        aurostd::string2tokens(options,tokens,",");
        
        if(tokens.size()==1) {
          if(tokens.at(0)=="usage" || tokens.at(0)=="USAGE") {
            init::ErrorOption(cout,options,"WYCCAR_RHT",
                              aurostd::liststring2string("aflow --wyccar[=tolerance| =tight| =loose] < POSCAR  default: (minimum_interatomic_distance)/100.0"));
            exit(0);
          }
        }
        double default_tolerance=SYM::defaultTolerance(str);
        double tolerance = default_tolerance;
        if(tokens.size()==0) {
          tolerance=default_tolerance;
        }
        if(tokens.size()>=1 && tokens.at(0) != "--debug") {
          if(tokens.at(0).at(0) == 't' || tokens.at(0).at(0) == 'T'){ //Tight
            tolerance=default_tolerance;
          }
          else if(tokens.at(0).at(0) == 'l' || tokens.at(0).at(0) == 'L'){ //Loose
            tolerance=default_tolerance*10.0;
          }
          else{ 
            tolerance=aurostd::string2utype<double>(tokens.at(0));
          }
        }
        str.SpaceGroup_ITC(tolerance);
        SYM::printWyccar(File,str,verbose,cout); 
        _PROGRAMRUN=true;
      }
      //[OBSOLETE] (DX)if(vpflow.flag("WYCCAR_MANUAL_RHT")) {// [OBSOLETE] initsymmats();initglides();initsymops();SetTolerance(argv);
      //[OBSOLETE] (DX)  xstructure str(cin); double tolerance=0.001; int man_int=-1; if(argv.size()==2) {tolerance=0.001;}; if(argv.size()>=3) {tolerance=atof(argv.at(2).c_str()); man_int=atof(argv.at(3).c_str());}; str.SpaceGroup_ITC(tolerance,man_int); printWyccar(File,str,verbose,cout); _PROGRAMRUN=true;}
      //End Richard's Functions
      // DX - END
    
    }
    // *********************************************************************
    if(argv.size()==3 && !_PROGRAMRUN) {
      // [OBSOLETE] if(vpflow.flag("ANGLES")) {pflow::ANGLES(cin,aurostd::args2utype(argv,"--angle",0.0)); _PROGRAMRUN=true;}
      // [OBSOLETE] if(vpflow.flag("BANDS")) {pflow::BANDS(argv,cin); _PROGRAMRUN=true;}
      // [OBSOLETE] if(vpflow.flag("BZDIRECTION")) {cout << pflow::BZDirectionsLATTICE(argv); _PROGRAMRUN=true;}
      // [OBSOLETE] if(vpflow.flag("CAGES") && !AFLOW_PTHREADS::FLAG) {pflow::CAGES(aflags,cin,aurostd::args2utype(argv,"--cages",-1.0)); _PROGRAMRUN=true;}
      // [OBSOLETE] if(vpflow.flag("CAGES") && AFLOW_PTHREADS::FLAG) { pflow::CAGES(aflags,cin,aurostd::args2utype(argv,"--cages",-1.0)); _PROGRAMRUN=true;}  // fix around Thu Apr 18 23:29:43 EDT 2013
      if(vpflow.flag("CHGINT")) {pflow::CHGINT(argv); _PROGRAMRUN=true;}
      // [OBSOLETE] if(vpflow.flag("DATA1")) {pflow::DATA1(argv,cin); _PROGRAMRUN=true;}
      // [OBSOLETE] if(vpflow.flag("DISP")) {pflow::DISP(cin,aurostd::args2utype(argv,"--disp",0.0)); _PROGRAMRUN=true;}
      // [OBSOLETE] if(vpflow.flag("DIST")) {pflow::DIST(cin,aurostd::args2utype(argv,"--dist",0.0)); _PROGRAMRUN=true;}
      // [OBSOLETE] if(vpflow.flag("EWALD")) {pflow::EWALD(argv,cin); _PROGRAMRUN=true;}
      if(EXTRACT_KPOINTS!="nan") {cout << pflow::EXTRACT_xcar(aflags,argv,"KPOINTS",EXTRACT_KPOINTS); _PROGRAMRUN=true;}
      if(EXTRACT_INCAR!="nan") {cout << pflow::EXTRACT_xcar(aflags,argv,"INCAR",EXTRACT_INCAR); _PROGRAMRUN=true;}
      if(EXTRACT_POSCAR!="nan") {cout << pflow::EXTRACT_xcar(aflags,argv,"POSCAR",EXTRACT_POSCAR); _PROGRAMRUN=true;}
      if(EXTRACT_POTCAR!="nan") {cout << pflow::EXTRACT_xcar(aflags,argv,"POTCAR",EXTRACT_POTCAR); _PROGRAMRUN=true;}
      if(vpflow.flag("EXTRACT_SYMMETRY")) {cout << pflow::EXTRACT_Symmetry(aflags,argv); _PROGRAMRUN=true;}
      // [OBSOLETE] if(vpflow.flag("SG::FINDSYM_PRINT")) {pflow::FINDSYM(argv,0,cin); _PROGRAMRUN=true;}
      // [OBSOLETE] if(vpflow.flag("SG::FINDSYM_EXEC")) {pflow::FINDSYM(argv,1,cin); _PROGRAMRUN=true;}
      if(vpflow.flag("HNF")) {pflow::HNF(argv,cin,cout); _PROGRAMRUN=true;}
      if(vpflow.flag("HNFTOL")) {pflow::HNFTOL(argv,cin,cout); _PROGRAMRUN=true;}
      if(vpflow.flag("KPATH")) {pflow::KPATH(cin,aurostd::args2attachedutype<double>(argv,"--grid=",16.0),vpflow.flag("WWW")); _PROGRAMRUN=true;}
      // [OBSOLETE] if(vpflow.flag("INFLATE_LATTICE")) {cout << pflow::INFLATE_LATTICE(cin,aurostd::args2utype(argv,"--inflate_lattice|--ilattice",1.0)); _PROGRAMRUN=true;}
      // [OBSOLETE] if(vpflow.flag("INFLATE_VOLUME")) {cout << pflow::INFLATE_VOLUME(cin,aurostd::args2utype(argv,"--inflate_volume|--ivolume",1.0)); _PROGRAMRUN=true;}
      if(vpflow.flag("JMOLGIF")) {pflow::JMOLAnimation(cin,argv); _PROGRAMRUN=true;}
      if(vpflow.flag("INSPHERE")) {pflow::XYZINSPHERE(cin,aurostd::args2utype(argv,"--insphere",0.0)); _PROGRAMRUN=true;}
      // [OBSOLETE] if(vpflow.flag("MAXATOMS")) {cout << pflow::ATOMSMAX(argv,cin); _PROGRAMRUN=true;}
      if(vpflow.flag("NANOPARTICLE")) {cout << pflow::NANOPARTICLE(cin,aurostd::args2xvectorutype<double>(argv,"--nanoparticle",argv.size()-2)); _PROGRAMRUN=true;}
      if(vpflow.flag("POCC")) {pflow::POCC(argv); _PROGRAMRUN=true;}
      // [OBSOLETE] if(vpflow.flag("PROTO")) {cout << pflow::PROTO_LIBRARIES(argv); _PROGRAMRUN=true;} // some options...
      if(vpflow.flag("QMVASP")) {pflow::QMVASP(argv); _PROGRAMRUN=true;}
    
      if(vpflow.flag("RAYTRACE")) {pflow::RAYTRACE(argv); _PROGRAMRUN=true;}
      // [OBSOLETE] if(vpflow.flag("SCALE")) {cout << pflow::SCALE(vpflow.getattachedscheme("SCALE"),cin); _PROGRAMRUN=true;}
      if(vpflow.flag("SETORIGIN")) {cout << pflow::SETORIGIN(cin,(uint) aurostd::string2utype<int>(argv.at(2))); _PROGRAMRUN=true;}
      if(vpflow.flag("SEWALD")) {pflow::SEWALD(argv,cin); _PROGRAMRUN=true;}
      if(vpflow.flag("SGROUP")) { pflow::SGROUP(aflags,cin,aurostd::args2utype(argv,"--sgroup|--spacegroup",KBIN_SYMMETRY_SGROUP_RADIUS_DEFAULT)); _PROGRAMRUN=true;}
 
      if(vpflow.flag("WYCKOFF")) {cout << pflow::WYCKOFF(argv,cin); _PROGRAMRUN=true;}
  
    }
    // *********************************************************************
    if(argv.size()>=3 && !_PROGRAMRUN) {
      // [OBSOLETE]  if(vpflow.flag("ALPHA_COMPOUND")) {cout << pflow::ALPHACompound(argv); _PROGRAMRUN=true;}
      // [OBSOLETE] if(vpflow.flag("ALPHA_SPECIES")) {cout << pflow::ALPHASpecies(argv); _PROGRAMRUN=true;}
      if(vpflow.flag("MISCIBILITY")) {cout << pflow::MISCIBILITY(argv); _PROGRAMRUN=true;}
      // [OBSOLETE] if(vpflow.flag("CE::SQS")) {pflow::SQS(argv); _PROGRAMRUN=true;}
      // calculated stuff
      if(vpflow.flag("XPLUG")) {aflowlib::XPLUG(argv); _PROGRAMRUN=true;}
    }
    // *********************************************************************
    if(argv.size()==4 && !_PROGRAMRUN) {
      // [OBSOLETE] if(vpflow.flag("CAGES") && AFLOW_PTHREADS::FLAG) {pflow::CAGES(aflags,cin,-1.0); _PROGRAMRUN=true;}
      if(vpflow.flag("JOINSTRLIST")) {pflow::JOINSTRLIST(argv); _PROGRAMRUN=true;}
      if(vpflow.flag("MAKESTRLIST")) {pflow::MAKESTRLIST(argv); _PROGRAMRUN=true;}
      if(vpflow.flag("NANOPARTICLE")) {cout << pflow::NANOPARTICLE(cin,aurostd::args2xvectorutype<double>(argv,"--nanoparticle",argv.size()-2)); _PROGRAMRUN=true;}
      if(vpflow.flag("PDOS")) {pflow::PDOS(argv); _PROGRAMRUN=true;}
      if(vpflow.flag("PLANEDENS")) {pflow::PLANEDENS(argv); _PROGRAMRUN=true;}
      if(vpflow.flag("PRIM") && AFLOW_PTHREADS::FLAG) {cout << pflow::PRIM(cin,0); _PROGRAMRUN=true;}
      if(vpflow.flag("QMVASP")) {pflow::QMVASP(argv); _PROGRAMRUN=true;}
      if(vpflow.flag("RBANAL")) {pflow::RBANAL(argv); _PROGRAMRUN=true;}
      if(vpflow.flag("SUMPDOS")) {pflow::SUMPDOS(argv); _PROGRAMRUN=true;}
      if(vpflow.flag("SWAP")) {cout << pflow::xstrSWAP(argv,cin); _PROGRAMRUN=true;}
      if(vpflow.flag("XFIXX")) {aflowlib::XFIX_LIBRARY_ALL("LIB2",argv); _PROGRAMRUN=true;}
    }
    // *********************************************************************
    if(argv.size()==5 && !_PROGRAMRUN) {
      // [OBSOLETE] if(vpflow.flag("CAGES") && AFLOW_PTHREADS::FLAG) {pflow::CAGES(aflags,cin,aurostd::args2utype(argv,"--cages",-1.0)); _PROGRAMRUN=true;}
      // [OBSOLETE] if(vpflow.flag("RDF")) {pflow::RDF(cin,aurostd::args2xvectorutype<double>(argv,"--rdf",3)); _PROGRAMRUN=true;}
      // [OBSOLETE] if(vpflow.flag("MILLER")) {cout << pflow::MILLER(argv,cin) << endl; _PROGRAMRUN=true;}
      if(vpflow.flag("QMVASP")) {pflow::QMVASP(argv); _PROGRAMRUN=true;}
      if(vpflow.flag("RBDIST")) {pflow::RBDIST(argv); _PROGRAMRUN=true;}
      if(vpflow.flag("SETCM")) {cout << pflow::SETCM(cin,aurostd::args2xvectorutype<double>(argv,"--setcm",3)); _PROGRAMRUN=true;}
      if(vpflow.flag("SETORIGIN")) {cout << pflow::SETORIGIN(cin,aurostd::args2xvectorutype<double>(argv,"--setorigin",3)); _PROGRAMRUN=true;}
      if(vpflow.flag("CMPSTR")) {pflow::CMPSTR(argv); _PROGRAMRUN=true;}
      // superlattice
      // [OBSOLETE] if(vpflow.flag("CE::SUPERLATTICE")) {pflow::Superlattice(argv); _PROGRAMRUN=true;}
      // [OBSOLETE] if(vpflow.flag("XRD_DIST")) {pflow::GetAtomicPlaneDist(cin,argv); _PROGRAMRUN=true;}
    }
    // [OBSOLETE] // *********************************************************************
    // [OBSOLETE]  if(argv.size()>=5 && !_PROGRAMRUN) {
    // [OBSOLETE] if(vpflow.flag("SHIFT")) {cout << pflow::SHIFT(argv,cin); _PROGRAMRUN=true;}
    // [OBSOLETE] if(vpflow.flag("JMOL")) {pflow::JMOL( cin,aurostd::args2xvectorutype<int>(argv,"--jmol",3),argv); _PROGRAMRUN=true;}
    // [OBSOLETE] }
    // [OBSOLETE] // *********************************************************************
    // [OBSOLETE] if(argv.size()==6 && !_PROGRAMRUN) {
    // [OBSOLETE] if(vpflow.flag("INTPOL")) {pflow::INTPOL(argv); _PROGRAMRUN=true;}
    // [OBSOLETE]  if(vpflow.flag("MILLER")) {cout << pflow::MILLER(argv,cin) << endl; _PROGRAMRUN=true;}
    // [OBSOLETE] }
    // *********************************************************************
    // [OBSOLETE] if(argv.size()==7 && !_PROGRAMRUN) {
    // [OBSOLETE] if(vpflow.flag("MILLER")) {cout << pflow::MILLER(argv,cin) << endl; _PROGRAMRUN=true;}
    // [OBSOLETE] if(vpflow.flag("SHELL")) {pflow::SHELL(argv,cin); _PROGRAMRUN=true;}
    // [OBSOLETE] if(vpflow.flag("CE::CLUSTERS")) {pflow::Cluster(argv); _PROGRAMRUN=true;}
    // [OBSOLETE] if(vpflow.flag("CE::CLUSTEREXPANSION")) {AClusterExpansionMethodMain(argv); _PROGRAMRUN=true;}
    // [OBSOLETE] }
    // [OBSOLETE] // *********************************************************************
    // [OBSOLETE] if(argv.size()==8 && !_PROGRAMRUN) {
    // [OBSOLETE]     if(vpflow.flag("RDFCMP")) {pflow::RDFCMP(aurostd::args2xvectorutype<double>(argv,"--rdfcmp",4),string(argv.at(6)),string(argv.at(7))); _PROGRAMRUN=true;}
    // [OBSOLETE] if(vpflow.flag("CLAT")) {pflow::CLAT(argv); _PROGRAMRUN=true;}
    // [OBSOLETE] }
    // ----------------------------------------------- ERRORS

    if(!_PROGRAMRUN) {
      cerr << "aflow: the number of arguments is not correct" << endl;
      cerr << "Try \'aflow --help\' for more information" << endl;
      for(uint i=0;i<argv.size();i++)
	// cerr << "argv.at(" << i << ")=" << argv.at(i) << endl; // exit(0);
	cerr << argv.at(i) << " ";
      cerr << endl;
    }
    // *********************************************************************
    // exit(1);
    return 0; //CO
  }
} // namespace pflow

// ***************************************************************************
// pflow::XXX
// ***************************************************************************
/*
  extern "C" {
  int __vector_matrix_utilities_MOD_test_one(int *a);
  int __vector_matrix_utilities_MOD_test_two(int *a, int *b);
  double __vector_matrix_utilities_MOD_test_three(double *a, double *b);
  double __vector_matrix_utilities_MOD_test_four(int *dim, double *a);
  double __vector_matrix_utilities_MOD_test_five(int *dimi, int *dimj, double *a);
  void __cwrapper_MOD_aflow_reduce_to_shortest_basis(double &Ain,double &Bout,const double &eps);
  void __cwrapper_MOD_aflow_reduce_to_shortest_basis2(int &Ai,int &Aj, double &Ain,
  int &Bi,int &Bj, double **Bout,
  const double &eps);
  }


  // 0000000000e23934 T __vector_matrix_utilities_MOD_test_two
  // http://www.neurophys.wisc.edu/comp/docs/notes/not017.html
  //#define    gFortran
  //#include <cfortran.h> //  -DgFortran
  //#include "fortran.h"

  bool XXX(vector<string> argv,istream& input) {
  xstructure str(input,IOAFLOW_AUTO);
  cerr << "XXX TEST GUS" << endl;
  xmatrix<double> lattice(str.lattice),l(str.lattice);
  double *forpus;
  str.MinkowskiBasisReduction();
  cerr << "o" << endl;
  cerr << lattice << endl;
  cerr << "s" << endl;
  cerr << str.lattice << endl;
  double eps=1e-3;
  //__cwrapper_MOD_aflow_reduce_to_shortest_basis(&l[1][1],&lattice[1][1],&eps);
  __cwrapper_MOD_aflow_reduce_to_shortest_basis(l[1][1],lattice[1][1],1e-3);
  int Ai=3,Aj=3,Bi=0,Bj=0; __cwrapper_MOD_aflow_reduce_to_shortest_basis2(Ai,Aj,l[1][1],Bi,Bj,&forpus,eps);
  xmatrix<double> xorpus(Bi,Bj,forpus); free(forpus); // delete [] forpus;
  //  cout << Ai << Aj << Bi << Bj << endl;
  cout << "g" << endl;
  cout << lattice << endl;
  cout << "x" << endl;
  cout << xorpus << endl;
  exit(0);
  return TRUE;
  }

  bool _pflow::XXX(vector<string> argv,istream& input) {
  //  xstructure str(input,IOAFLOW_AUTO);
  cerr << "XXX TEST GUS" << endl;
  int imin=aurostd::string2utype<int>(argv.at(2)),imax=aurostd::string2utype<int>(argv.at(3));
  cerr << imin << endl;
  cerr << imax << endl;
  xvector<double> c(imin,imax);
  for(int i=imin;i<=imax;i++) c(i)=(double) i;
  xmatrix<double> d(imin,imax);
  for(int i=d.lrows;i<=d.urows;i++)                      // clear
  for(int j=d.lcols;j<=d.ucols;j++)                    // clear
  d[i][j]=i+j;
  cerr << d << endl;
  cerr << "GUS test=" << __vector_matrix_utilities_MOD_test_five(&d.rows,&d.cols,&d[d.lrows][d.lcols]) << endl;

  //  int c=__vector_matrix_utilities_MOD_test_one(&a);
  //  int c=__vector_matrix_utilities_MOD_test_three(&a,&b);
  //  cerr << LatticeDimensionSphere(str.lattice,string2utype<int>(argv.at(2))) << endl;
  //  cerr <<  c << endl;   int dim=c.rows; cerr << "c[imin]=" << c[imin] << endl;   cerr << "GUS test=" << __vector_matrix_utilities_MOD_test_four(&dim,&c[imin]) << endl;
  exit(0);
  return TRUE;
  }
*/

// ***************************************************************************
// pflow::CheckCommands
// ***************************************************************************
namespace pflow {
  bool CheckCommands(vector<string> argv,const vector<string> &cmds) {
    string _cmd;
    vector<string> tokens;
    bool found=FALSE;
    // check identities
    for(int i=argv.size()-1;i>=1&&!found;i--) {
      _cmd=aurostd::RemoveWhiteSpaces(string(argv.at(i)));
      for(uint j=0;j<cmds.size()&&!found;j++) {//cerr << _cmd << " " << cmds[j] << endl;
	if(_cmd==cmds.at(j)) found=TRUE;}
    }
    // check with =
    for(int i=argv.size()-1;i>=1&&!found;i--) {
      aurostd::RemoveWhiteSpaces(string(argv.at(i)));
      aurostd::string2tokens(aurostd::RemoveWhiteSpaces(string(argv.at(i))),tokens,"=");
      _cmd=tokens.at(0)+"=";
      for(uint j=0;j<cmds.size()&&!found;j++) {// cerr << _cmd << " " << cmds[j] << endl;
	if(_cmd==cmds.at(j)) found=TRUE;}
    }
    // not found
    if(!found) {
      cerr << aflow::Banner("BANNER_TINY") << endl;
      cerr << "ERROR - pflow::CheckCommands: command not found: " << _cmd << endl;
      return FALSE;
    }
    return TRUE;
  }
} // namespace pflow

// ***************************************************************************
// pflow::Intro_pflow
// ***************************************************************************
namespace pflow {
  string Intro_pflow(string x) {
    stringstream strstream;
    // intro(strstream);
    strstream << "\
******* BEGIN POSTPROCESSING MODE ****************************************************************** \n\
  "<< x<<" --help [-h] option_name \n\
  \n\
  "<< x<<" --abccar < POSCAR | WYCCAR \n\
  "<< x<<" --abinit < POSCAR \n\
  "<< x<<" --aims < POSCAR \n\
  "<< x<<" --ace < POSCAR \n\
  "<< x<<" --use_aflow.in=XXX \n\
  "<< x<<" --aflowin < POSCAR \n\
  "<< x<<" --aflowSG[_label,_number][=tolerance| =tight| =loose] < POSCAR \n\
  "<< x<<" --aflow-sym|--AFLOW-SYM|--AFLOWSYM|--aflowSYM|--aflowsym|--full_symmetry|--full_sym|--fullsym[=tolerance| =tight| =loose] [--no_scan] [--print=txt| =json] [--screen_only] [--mag|--magnetic|--magmom=[m1,m2,...|INCAR|OUTCAR]] < POSCAR \n\
  "<< x<<" --poscar2aflowin < POSCAR \n\
  "<< x<<" --angle=cutoff < POSCAR \n\
  "<< x<<" --agroup | --sitepointgroup[=tolerance| =tight| =loose] [--no_scan] [--print=txt| =json] [--screen_only] [--mag|--magnetic|--magmom=[m1,m2,...|INCAR|OUTCAR]] < POSCAR \n\
  "<< x<<" --agroup2 < POSCAR \n\
  "<< x<<" --agroup2m < POSCAR \n\
  "<< x<<" --alphabetic < POSCAR \n\
  "<< x<<" --alpha_compound=string1,string2.... \n\
  "<< x<<" --alpha_species=string1,string2... \n\
  "<< x<<" --aflowlib=entry \n\
  "<< x<<" --aflowlib_auid2aurl=auid1,auid2.... | --auid2aurl=... \n\
  "<< x<<" --aflowlib_aurl2auid=aurl1,aurl2.... [ --aurl2auid=... \n\
  "<< x<<" --aflowlib_auid2loop=auid1,auid2.... | --auid2loop=... \n\
  "<< x<<" --aflowlib_aurl2loop=aurl1,aurl2.... [ --aurl2loop=... \n\
  "<< x<<" [options] --bader -D DIRECTORY \n\
   options are:  --usage \n\
                 --critical_points|--cp \n\
                 --calculate=|--calc=bader|voronoi \n\
                 --nocalculate=|--nocalc=bader|voronoi \n\
                 --partition=|--part=neargrid|ongrid \n\
                 --refine_edge_method=|--rem=-1|-2|-3 \n\
                 --reference=|--ref=REF_FILE \n\
                 --vacuum=|--vac=off|auto|DENSITY_THRESHOLD \n\
                 --terminate=|--term=known|max \n\
                 --print_all=atom|bader|both \n\
                 --print_index=|--print_idx=atom|bader|both \n\
                 --print_select_atom=|--print_sel_atom=[LIST OR RANGE] \n\
                 --print_select_bader=|--print_sel_bader=[LIST OR RANGE] \n\
                 --print_sum_atom=[LIST OR RANGE] \n\
                 --print_sum_bader=[LIST OR RANGE] \n\
                 --quiet|--q \n\
                 --consolidate_atoms2species|--a2s \n\
                 --remove_bader_atoms|--rba \n\
                 --jvxl_all_species=|--jvxl=CUTOFF1,CUTOFF2...[::DOWNSAMPLE1,DOWNSAMPLE2,..]|CUTOFF1[,DOWNSAMPLE1:CUTOFF2,DOWNSAMPLE2:...]|CUTOFF[,DOWNSAMPLE] \n\
                 --keep=jvxl_only|--jvxl_only \n\
  "<< x<<" --bands=PROOUT < POSCAR [AFLOW: NEED VERIFICATION] \n\
  "<< x<<" --bandgap[=bands_directory1[,bands_directory2,...]] \n\
  "<< x<<" --bandgaps < file_containing_bands_directories \n\
  "<< x<<" --bandgapdos < DOSCAR \n\
  "<< x<<" --bandgaplistdos < DOSCAR_list \n\
  "<< x<<" --bandstructure | --bs \n\
  "<< x<<" --bzdirections | --bzd < POSCAR \n\
  "<< x<<" --bzdirections= | --bzd=LATTICE \n\
  "<< x<<" --BZmax < POSCAR \n\
  "<< x<<" --bzplot | --plotbz < POSCAR \n\
  "<< x<<" --bzplotuseKPOINTS=KPOINTS < POSCAR \n\
  "<< x<<" --bzplotdata < POSCAR \n\
  "<< x<<" --bzplotdatauseKPOINTS=KPOINTS < POSCAR \n\
  "<< x<<" --cart [-c] < POSCAR \n\
  "<< x<<" [options] --chgcar2jvxl=|--c2j=CHGCAR11[,CHGCAR2,...]::CUTOFF1,CUTOFF2...[::DOWNSAMPLE1,DOWNSAMPLE2,...]|CHGCAR1,CUTOFF1[,DOWNSAMPLE1:CHGCAR2,CUTOFF2[,DOWNSAMPLE2:...]]|CHGCAR,CUTOFF[,DOWNSAMPLE] \n\
   options are:  --usage \n\
                 --output=|--o=OUTPUT_FILE \n\
  "<< x<<" [options]--chgdiff=CHGCAR1,CHGCAR2 \n\
   options are:  --usage \n\
                 --output=|--o=CHGCAR_OUT \n\
  "<< x<<" --chgsum=CHGCAR1,CHGCAR2,... \n\
   options are:  --usage \n\
                --output=|--o=CHGCAR_OUT \n\
  "<< x<<" --check_integrity | --checki \n\
  "<< x<<" --cif < POSCAR \n\
  "<< x<<" --clean -D DIRECTORY \n\
  "<< x<<" --clean_all < LIST_DIRECTORIES \n\
  "<< x<<" --convex_hull=|--chull --alloy=MnPdPt[,AlCuZn,...] [chull_options] [--destination=[DIRECTORY]] \n\
    chull_options: \n\
                 --usage \n\
                 --output=|--o=|--print=|--p=latex|pdf|json|text \n\
                 --image_only|--imageonly|--image \n\
                 --no_document|--nodocument|--no_doc|--nodoc|--full_page_image|--fullpageimage \n\
                 --document_only|--documentonly|--doc_only|--doconly|--doc \n\
                 --keep=tex|--keep_tex|--keeptex|--tex \n\
                 --keep=log|--keep_log|--keeplog|--log \n\
                 \n\
                 LOADING OPTIONS: \n\
                 --load_library=|--loadlibrary=|--ll=icsd|lib1|lib2|lib3 \n\
                 --load_API|--load_api|--loadapi|--lapi|--api \n\
                 --load_entries_entry_output|--loadentriesentryoutput|--leo \n\
                 --neglect=|--ban=aflow:bb0d45ab555bc208);aflow:fb9eaa58604ce774 \n\
                 --see_neglect|--seeneglect|--sn \n\
                 --remove_extreme_points=|--removeextremepoints=|--remove_extrema=|--removeextrema=|--rep=-1000 \n\
                 --entropic_temperature|--entropictemperature|--entroptemp \n\
                 \n\
                 ANALYSIS OPTIONS: \n\
                 --stability_criterion=|--stabilitycriterion=|--stable_criterion=|--scriterion=|--sc=aflow:bb0d45ab555bc208,aflow:fb9eaa58604ce774 \n\
                 --distance_to_hull=|--dist2hull=0.25,0.25 \n\
                 --skip_structure_comparison|--skipstructruecomparison|--skipstructcomp|--ssc \n\
                 --include_unreliable_hulls|--include_unreliable|--iuh \n\
                 --include_outliers|--io \n\
                 --force \n\
                 \n\
                 GENERAL LATEX OPTIONS: \n\
                 --latex_output|--latexoutput \n\
                 --latex_interactive|--latexinteractive \n\
                 --light_contrast|--lightcontrast|--lc \n\
                 --large_font|--largefont|--large|--lf \n\
  "<< x<<" --corners | --corner < POSCAR \n\
  "<< x<<" --data[=tolerance| =tight| =loose] [--no_scan] [--print=txt| =json] < POSCAR \n\
  "<< x<<" --data1=rcut < POSCAR \n\
  "<< x<<" --data2 < POSCAR \n\
  "<< x<<" --debye=THERMO[.bz2] \n\
  "<< x<<" --diff=POSCAR1,POSCAR2 \n\
  "<< x<<" --disp=cutoff < POSCAR \n\
  "<< x<<" --dist=cutoff < POSCAR \n\
  "<< x<<" --delta_kpoints=number [or --dkpoints=number | -dkpoints=number | -dk=number] < POSCAR \n\
  "<< x<<" --edata[=tolerance| =tight| =loose] [--no_scan] [--print=txt| =json] < POSCAR \n\
  "<< x<<" --effective-mass | --em DIRECTORY \n\
  "<< x<<" --eigcurv=DIRECTORY(with bands) \n\
  "<< x<<" --enum | --multienum < POPSCAR \n\
  "<< x<<" --enumsort | --multienumsort < POPSCAR \n\
  "<< x<<" --equivalent | --equiv | --inequivalent | --inequiv | --iatoms | --eatoms[=tolerance| =tight| =loose] [--no_scan] [--print=txt| =json] [--mag|--magnetic|--magmom=[m1,m2,...|INCAR|OUTCAR]] < POSCAR \n\
  "<< x<<" --ewald[=eta] < POSCAR \n\
  "<< x<<" --findsym[=tolerance_relative: default " << DEFAULT_FINDSYM_TOL << "] < POSCAR \n\
  "<< x<<" --findsym_print[=tolerance_relative: default " << DEFAULT_FINDSYM_TOL << "] < POSCAR \n\
  "<< x<<" --findsymSG[_label,_number][=tolerance_relative: default " << DEFAULT_FINDSYM_TOL << "] < POSCAR \n\
  "<< x<<" --frac [-f,-d,--fract,--fractional,--direct] < POSCAR \n\
  "<< x<<" --getTEMP [--runstat | --runbar | --refresh=X | --warning_beep=T | --warning_halt=T ] \n\
  "<< x<<" --gulp < POSCAR \n\
  "<< x<<" --identical < POSCAR \n\
  "<< x<<" --incell < POSCAR \n\
  "<< x<<" --incompact < POSCAR \n\
  "<< x<<" --insphere radius < POSCAR \n\
  "<< x<<" --inwignerseitz [--inws] < POSCAR \n\
  "<< x<<" --inflate_lattice=coefficient | --ilattice coefficient=coefficient < POSCAR \n\
  "<< x<<" --inflate_volume=coefficient | --ivolume=coefficient < POSCAR \n\
  "<< x<<" --jgif < POSCAR \n\
  "<< x<<" --jmol[=n1[,n2[,n3[,color[,true/false]]]]] < POSCAR \n\
  "<< x<<" --hnf n < POSCAR \n\
  "<< x<<" --hnftol epsilon < PARTCAR \n\
  "<< x<<" --hnfcell < POSCAR \n\
  "<< x<<" --kpath [--grid=XX] < POSCAR \n\
  "<< x<<" --kpoints=KDENS [or --kppra=KDENS,-k=KDENS] < POSCAR \n\
  "<< x<<" --maxatoms=N | --max_atoms=N | --atoms_max=N | --atomsmax=N < POSCAR \n\
  "<< x<<" --msi < POSCAR \n\
  "<< x<<" --latticehistogram < POSCAR \n\
  "<< x<<" --latticereduction < POSCAR \n\
  "<< x<<" --lattice_type | --lattice | --lattice_crystal < POSCAR \n\
  "<< x<<" --lattice_lattice_type | --lattice_lattice < POSCAR \n\
  "<< x<<" --latticehistogram < POSCAR \n\
  "<< x<<" --use_LOCK=XXX \n\
  "<< x<<" --ltcell=a11,a12,a13,a21,a22,a23,a31,a32,a33 < POSCAR \n\
  "<< x<<" --ltcell=a11,a22,a33 < POSCAR \n\
  "<< x<<" --ltcell=file < POSCAR \n\
  "<< x<<" --ltcellfv=v1,v2,v3,phi < POSCAR \n\
  "<< x<<" --magpara=directory | --magpara \n\
  "<< x<<" --minkowski_basis_reduction | --minkowski| --mink < POSCAR \n\
  "<< x<<" --miscibility | --mix string \n\
  "<< x<<" --mom < POSCAR \n\
  "<< x<<" --natoms | --numatoms < POSCAR \n\
  "<< x<<" --nbondxx < POSCAR \n\
  "<< x<<" --names | --species A1 A2 ... < POSCAR \n\
  "<< x<<" --nanoparticle radius distance < POSCAR \n\
  "<< x<<" --ndata < POSCAR \n\
  "<< x<<" --niggli < POSCAR \n\
  "<< x<<" --nn < POSCAR \n\
  "<< x<<" --noorderparameter < POSCAR \n\
  "<< x<<" --nosd < POSCAR\n\
  "<< x<<" --numnames A1 A2 ... < POSCAR \n\
  "<< x<<" --nspecies | --numspecies< POSCAR \n\
  "<< x<<" --OrthoDefect < POSCAR \n\
  "<< x<<" --pdb < POSCAR \n\
  "<< x<<" --rsm < POSCAR \n\
  "<< x<<" --rsm --z atom_num1 atom_num2 atom_num3 < POSCAR \n\
  "<< x<<" --rsm --z atom_symbol1 atom_symbol2 atom_symbol3 < POSCAR \n\
  "<< x<<" --pearson_symbol | --pearson < POSCAR \n\
  "<< x<<" --platon[=EQUAL | EXACT][,ang,d1,d2,d3] < POSCAR | platonSG \n\
  "<< x<<" --platonSG[_label,_number][=EQUAL | EXACT][,ang,d1,d2,d3] < POSCAR \n\
  "<< x<<" --plotband[=directory[,DOS_Emin[,DOS_Emax[,DOSSCALE]]]]] \n\
  "<< x<<" --plotband_spinsplit[=directory[,DOS_Emin[,DOS_Emax[,DOSSCALE]]]]] \n\
  "<< x<<" --plotband[=directory[,DOS_Emin[,DOS_Emax[,DOSSCALE]]]]] \n\
  "<< x<<" --plotbanddos[=directory[,DOS_Emin[,DOS_Emax[,DOSSCALE]]]]] \n\
  "<< x<<" --plotdos[=directory[,DOS_Emin[,DOS_Emax[,DOSSCALE]]]]] \n\
  "<< x<<" --plotdosweb[=directory[,DOS_Emin[,DOS_Emax[,DOSSCALE]]]] \n\
  "<< x<<" --plotpedos[=directory[,number_atom[,DOS_Emin[,DOS_Emax[,DOSSCALE]]]]] \n\
  "<< x<<" --plotpedosall[=directory[,DOS_Emin[,DOS_Emax[,DOSSCALE]]]] \n\
  "<< x<<" --plotpedosall_nonquivalent[=directory[,DOS_Emin[,DOS_Emax[,DOSSCALE]]]] \n\
  "<< x<<" --plotphonondispersion | --pphdis ( --rcm | --meV | --THz |  --hz ) [--print=eps | --print=pdf | --print=png | --print=jpg | --print=gif ] DIRECTORY\n\
  "<< x<<" --pomass[=directory] \n\
  "<< x<<" --pomass_atom[=directory] \n\
  "<< x<<" --pomass_cell[=directory] \n\
  "<< x<<" --pocc_input \n\
  "<< x<<" --pocc_dos[=directory[,T[,DOS_Emin[,DOS_Emax[,DOSSCALE]]]]] \n\
  "<< x<<" --pocc_mag[=directory[,T[,DOS_Emin[,DOS_Emax[,DOSSCALE]]]]] \n\
  "<< x<<" --pocc_bandgap[=directory[,T[,DOS_Emin[,DOS_Emax[,DOSSCALE]]]]] \n\
  "<< x<<" --poscar < ABCCAR | WYCCAR \n\
  "<< x<<" --poscar2aflowin < POSCAR \n\
  "<< x<<" --poscar2wyckoff < POSCAR \n\
  "<< x<<" --poscar2gulp < POSCAR \n\
  "<< x<<" [options] --prepare_chgcar_4_jmol=|--prep4jmol=CHGCAR1[,CHGACAR2,...] \n\
   options are:  --usage \n\
                 --outcar=OUTCAR \n\
                 --zip \n\
  "<< x<<" --prim < POSCAR \n\
  "<< x<<" --prim2 < POSCAR \n\
  "<< x<<" --primr | --fastprimitivecell | --fprim < POSCAR \n\
  "<< x<<" --qe < POSCAR \n\
  "<< x<<" --qmvasp [--static] [-D directory] \n\
  "<< x<<" --rasmol[=n1[,n2[,n3]]] < POSCAR \n\
  "<< x<<" --revsg [#] [n] [l] [m] \n\
  "<< x<<" --rm_atom iatom < POSCAR \n\
  "<< x<<" --rm_copies < POSCAR \n\
  "<< x<<" --rdf[=rmax[,nbins[,sigma]]] < POSCAR \n\
  "<< x<<" --scale=s < POSCAR \n\
  "<< x<<" --sd A1 A2 ... < POSCAR \n\
  "<< x<<" --setcm cm1 cm2 cm3 < POSCAR \n\
  "<< x<<" --setorigin r1 r2 r3 | atom# < POSCAR \n\
  "<< x<<" --sewald eta < POSCAR \n\
  "<< x<<" --sgdata|--space_group_data[=tolerance| =tight| =loose] [--no_scan] [--print=txt| =json] [--mag|--magnetic|--magmom=[m1,m2,...|INCAR|OUTCAR]] < POSCAR \n\
  "<< x<<" --shell=ns,r1,r2,name,dens < POSCAR \n\
  "<< x<<" --shift=Sx,Sy,Sz[,cCdD] < POSCAR \n\
  "<< x<<" --sitepointgroup | --agroup < POSCAR \n\
  "<< x<<" --slab=h,k,l[,#filled_layers[,#vacuum layers]] < POSCAR \n\
  "<< x<<" --spacegroup radius < POSCAR \n\
  "<< x<<" --species < POSCAR \n\
  "<< x<<" --statdiel OUTCAR* \n\
  "<< x<<" --std_conv | --standard_conventional | --sconv | --sc < POSCAR \n\
  "<< x<<" --std_prim | --standard_primitive | --sprim | --sp < POSCAR \n\
  "<< x<<" --suffix=[directory,]\"from2to\" where from/to=[n=none;r1=relax1;r2=relax2;r3=relax3;s=static;b=bands] or without abbreviations \n\
  "<< x<<" --supercell=a11,a12,a13,a21,a22,a23,a31,a32,a33 < POSCAR \n\
  "<< x<<" --supercell=a11,a22,a33 < POSCAR \n\
  "<< x<<" --supercell=file < POSCAR \n\
  "<< x<<" --supercell_strlist=a11,a12,a13,a21,a22,a23,a31,a32,a33,strlist \n\
  "<< x<<" --supercell_strlist=a11,a22,a33,strlist \n\
  "<< x<<" --supercell_strlist=file,strlist \n\
  "<< x<<" --swap specie0 specie1 < POSCAR \n\
  "<< x<<" --uffenergy | --ue  < POSCAR \n\
  "<< x<<" --vasp < GEOM \n\
  "<< x<<" --volume=x | --volume*=x | --volume+=x < POSCAR \n\
  "<< x<<" --wyccar [TOL] < POSCAR \n\
  "<< x<<" --wyccman | --WyccarManual | --wm [TOL] [ITERATION] < POSCAR \n\
  "<< x<<" --xray=lambda < POSCAR\n\
  "<< x<<" --xrd_dist=h,k,l < POSCAR\n\
  "<< x<<" --xyz[=n1[,n2[,n3]]] < POSCAR \n\
  "<< x<<" --xyzwignerseitz [--xyzws] < POSCAR \n\
  "<< x<<" --zval[=directory] \n\
  "<< x<<" --zval_atom[=directory] \n\
  "<< x<<" --zval_cell[=directory] \n\
  \n\
  "<< x<<" --pointgroup | --pgroup[=tolerance| =tight| =loose] [--no_scan] [--print=txt| =json] [--screen_only] < POSCAR \n\
  "<< x<<" --pointgroup_crystal | --pgroup_crystal | --pgroup_xtal | --pgroupx | --pgroupX[=tolerance| =tight| =loose] [--no_scan] [--print=txt| =json] [--screen_only] [--mag|--magnetic|--magmom=[m1,m2,...|INCAR|OUTCAR]] < POSCAR \n\
  "<< x<<" --pgl < POSCAR \n\
  "<< x<<" --pgx < POSCAR \n\
  "<< x<<" --factorgroup | --fgroup[=tolerance| =tight| =loose] [--no_scan] [--print=txt|=json] [--screen_only] [--mag|--magnetic|--magmom=[m1,m2,...|INCAR|OUTCAR]] < POSCAR \n\
  "<< x<<" --sitepointgroup | --agroup < POSCAR \n\
  "<< x<<" --spacegroup | --sgroup[=tolerance| =tight| =loose] [--no_scan] [--print=txt|=json] [--screen_only] [--radius] [--mag|--magnetic|--magmom=[m1,m2,...|INCAR|OUTCAR]] < POSCAR \n\
  "<< x<<" --pointgroupklattice | --pgroupk[=tolerance| =tight| =loose] [--no_scan] [--print=txt| =json] [--screen_only] < POSCAR \n\
  "<< x<<" --pointgroupkcrystal | --pgroupk_xtal[=tolerance| =tight| =loose] [--no_scan] [--print=txt| =json] [--screen_only] [--mag|--magnetic|--magmom=[m1,m2,...|INCAR|OUTCAR]] < POSCAR \n\
  \n\
 Prototypes HTQC \n\
  "<< x<<" --prototypes | --protos \n\
  "<< x<<" [options] --proto=label*[:speciesA*[:speciesB*]..[:volumeA*[:volumeB*].. | :volume]] [--params=... [--hex]] \n \
  \n\
 Prototypes ICSD (from the local machine or aflowlib servers) \n\
  "<< x<<" [--server=nnnnnnnnnn] --prototypes_icsd[=N] | --protos_icsd[=N] \n\
  "<< x<<" [--server=nnnnnnnnnn] [--vasp | --qe | --abinit | --aims] --proto_icsd=label \n\
  \n\
 Order Parameter [EXPERIMENTA] \n\
  "<< x<<" --order_parameter [--order_parameter_sum XXX] [--Li8Mn16O32] < POSCAR [EXPERIMENTA] \n\
  \n\
 Partial Occupation [EXPERIMENTA] \n\
  "<< x<<" --partial_occupation | --partialoccupation | --pocc < POSCAR [EXPERIMENTA] \n\
  \n\
 AFLOW Operations \n\
  "<< x<<" [options] --aflow_proto=label*:speciesA*[:speciesB*][:volumeA*[:volumeB*] | :volume] \n\
  "<< x<<" [options] --aflow_proto_icsd=label* potential_type (ICSD capable) \n\
   options are: --potential=pot_LDA | pot_GGA | potpaw_LDA | potpaw_GGA | potpaw_PBE \n\
                --potential_complete \n\
                --module=[APL | QHA | AAPL] \n\
                --apl_supercell=NxNxN \n\
                --usage \n\
                --missing \n\
                --noautopp \n\
                --kppra=NNNN (default: DEFAULT_KPPRA in .aflow.rc) --kppra_static=NNNN (default: DEFAULT_KPPRA_STATIC in .aflow.rc) --bands_grid=NNNN  (default: DEFAULT_BANDS_GRID in .aflow.rc ) \n\
                --enmax_multiply=NNNN (default: VASP_PREC_ENMAX_XXXX in .aflow.rc) \n\
                --pressure=0,1,2 (kB) (default:0.0) \n\
                --potim=XXX (default 0.05) (VASP) \n\
                --relax_type=[ALL | IONS | CELL_SHAPE | CELL_VOLUME | IONS_CELL_VOLUME] \n\
                --relax_mode=[ENERGY | FORCES | ENERGY_FORCES | FORCES_ENERGY] (default: DEFAULT_VASP_FORCE_OPTION_RELAX_MODE_SCHEME in .aflow.rc) (VASP) \n\
                --precision=[(LOW | MEDIUM | NORMAL | HIGH | ACCURATE), PRESERVED] (default: DEFAULT_VASP_FORCE_OPTION_PREC_SCHEME in .aflow.rc) (VASP) \n\
                --algorithm=[(NORMAL | VERYFAST | FAST | ALL | DAMPED), PRESERVED] (default: DEFAULT_VASP_FORCE_OPTION_ALGO_SCHEME in .aflow.rc) (VASP) \n\
                --metagga=[TPSS | RTPSS | M06L | MBJL | SCAN | MS0 | MS1 | MS2 | NONE] (defaul: DEFAULT_VASP_FORCE_OPTION_METAGGA_SCHEME in .aflow.rc) (VASP) \n\
                --ivdw=[number_for_VASP_see_manual_for_IVDW | 0] (default: DEFAULT_VASP_FORCE_OPTION_IVDW_SCHEME in .aflow.rc) (VASP) \n\
                --type=[METAL | INSULATOR | SEMICONDUCTOR | DEFAULT] (default: DEFAULT_VASP_FORCE_OPTION_TYPE_SCHEME in .aflow.rc (VASP) \n\
                --convert_unit_cell=[SPRIM, SCONV, NIGGLI, MINK, INCELL, COMPACT, WS, CART, FRAC, PRES] \n\
                --volume_plus_equal=XXX \n\
                --volume_multiply_equal=XXX \n\
                --no_volume_adjustment \n\
                --ediffg=XXX  (default: DEFAULT_VASP_PREC_EDIFFG in .aflow.rc) \n\
                --ldau2 \n\
                --noldau2 \n\
                --bands \n\
                --neglect_nomix \n\
                --stdout \n\
                --qe \n\
                --abinit \n\
                --aims \n\
                --params=... { check aflow --readme=anrl } \n\
                --hex        { check aflow --readme=anrl } \n\
                --list \n\
  "<< x<<" --extract_kpoints | --xkpoints " << _AFLOWIN_ << " \n\
  "<< x<<" --extract_incar | --xincar " << _AFLOWIN_ << " \n\
  "<< x<<" --extract_poscar | --xposcar " << _AFLOWIN_ << " \n\
  "<< x<<" --extract_potcar | --xpotcar " << _AFLOWIN_ << " \n\
  \n\
 CAGES SEARCH \n\
  "<< x<<" --cages[=roughness] < POSCAR \n\
  \n\
 HKL and HKL searches for surfaces \n\
  "<< x<<" --miller=h,k,l[,nlayer[,elayer]] < POSCAR \n\
  "<< x<<" --hkl=h,k,l[,bond] < POSCAR \n\
  "<< x<<" --hkl_search[=khlmax[,bond[,step]]] < POSCAR \n\
  "<< x<<" --hkl_search_simple[=cutoff[,bond[,khlmax[,step]]]] < POSCAR	\n\
  "<< x<<" --hkl_search_complete[=cutoff[,bond[,khlmax[,step]]]] < POSCAR \n\
  \n\
  "<< x<<" --chgint CHGCAR \n\
  "<< x<<" --clat=a,b,c,alpha,beta,gamma \n\
  "<< x<<" --cmp_str POSCAR1 POSCAR2 rcut \n\
  "<< x<<" --compare=a,b,c,d,e,f,g,h,k,j,i,l \n\
  "<< x<<" --intpol=file1,file2,nimages,nearest_image_flag \n\
  "<< x<<" --join_strlist strlist1 strlist2 \n\
  "<< x<<" --make_strlist OUTCAR XDATCAR \n\
  "<< x<<" --pdos pdos.in PROOUT [AFLOW: NEED VERIFICATION]\n\
  "<< x<<" --planedens dens2d.in CHGCAR \n\
  "<< x<<" --pocc PROOUT  [AFLOW: NEED VERIFICATION]\n\
  "<< x<<" --raytrace rtfile \n\
  "<< x<<" --rbanal nim nearest_image_flag \n\
  "<< x<<" --rbdist POSCAR1 POSCAR2 n|N|e|E \n\
  "<< x<<" --rdfcmp=rmax,nbins,sigma,nshmax,POSCAR1,POSCAR2 \n\
  "<< x<<" --spline npt < file \n\
  "<< x<<" --sumpdos pdos.in PROOUT [AFLOW: NEED VERIFICATION] \n\
  \n\
 ICSD MODE \n\
  ICSD input can be ternary.icsd, binary.icsd or other collections \n\
  of databases like sub-databases generated with and/or as follow. \n\
  "<< x<<" --icsd symbol/Z symbol/Z symbol/Z < input.icsd \n\
  "<< x<<" --icsd_alllessthan symbol/Z \n\
  "<< x<<" --icsd_allmorethan symbol/Z \n\
  "<< x<<" --icsd_basislessthan #basis < input.icsd \n\
  "<< x<<" --icsd_basismorethan #basis < input.icsd \n\
  "<< x<<" --icsd_check_raw \n\
  "<< x<<" --icsd_chem MgB4 < input.icsd \n\
  "<< x<<" --icsd_cubic < input.icsd \n\
  "<< x<<" --icsd_triclinic < input.icsd \n\
  "<< x<<" --icsd_monoclinic < input.icsd \n\
  "<< x<<" --icsd_orthorhombic < input.icsd \n\
  "<< x<<" --icsd_tetragonal < input.icsd \n\
  "<< x<<" --icsd_rhombohedral < input.icsd \n\
  "<< x<<" --icsd_trigonal < input.icsd \n\
  "<< x<<" --icsd_hexagonal < input.icsd \n\
  "<< x<<" --icsd_cubic --icsd_orthorhombic < input.icsd \n\
  "<< x<<" --icsd_tri < input.icsd \n\
  "<< x<<" --icsd_mcl < input.icsd \n\
  "<< x<<" --icsd_mclc < input.icsd \n\
  "<< x<<" --icsd_orc < input.icsd \n\
  "<< x<<" --icsd_orcc < input.icsd \n\
  "<< x<<" --icsd_orcf < input.icsd \n\
  "<< x<<" --icsd_orci < input.icsd \n\
  "<< x<<" --icsd_tet < input.icsd \n\
  "<< x<<" --icsd_bct < input.icsd \n\
  "<< x<<" --icsd_rhl < input.icsd \n\
  "<< x<<" --icsd_hex < input.icsd \n\
  "<< x<<" --icsd_cub < input.icsd \n\
  "<< x<<" --icsd_fcc < input.icsd \n\
  "<< x<<" --icsd_bcc < input.icsd \n\
  "<< x<<" --icsd_denslessthan X.X < input.icsd \n\
  "<< x<<" --icsd_densmorethan Y.Y < input.icsd \n\
  "<< x<<" --icsd_denslessthan X.X --icsd_densmorethan Y.Y < input.icsd \n\
  "<< x<<" --icsd_id #ICSD_ID < input.icsd \n\
  "<< x<<" --icsd_makelabel < input \n\
  "<< x<<" --icsd_lessthan symbol/Z < input.icsd \n\
  "<< x<<" --icsd_morethan symbol/Z < input.icsd \n\
  "<< x<<" --icsd_lessthan symbol/Z --icsd_morethan symbol/Z < input.icsd \n\
  "<< x<<" --icsd_listmetals \n\
  "<< x<<" --icsd_nobrokenbasis < input.icsd \n\
  "<< x<<" --icsd_nopartialocc < input.icsd \n\
  "<< x<<" --icsd_n_ary #species < input.icsd \n\
  "<< x<<" --icsd_proto #Nspecies1 #Nspecies2 #Nspecies3 ... < input.icsd \n\
  "<< x<<" --icsd_remove_all symbol/Z symbol/Z \n\
  "<< x<<" --icsd_remove_or symbol/Z symbol/Z \n\
  "<< x<<" --icsd_removemetals \n\
  "<< x<<" --icsd_sg #SG < input.icsd \n\
  "<< x<<" --icsd_sglessthan #SG < input.icsd \n\
  "<< x<<" --icsd_sgmorethan #SG < input.icsd \n\
  "<< x<<" --icsd_sgmorethan #SG --icsd_sglessthan #SG < input.icsd \n\
  "<< x<<" --icsd_unique < input.icsd \n\
  "<< x<<" --icsd2aflowin < input.icsd \n\
  "<< x<<" --icsd2poscar < input.icsd \n\
  "<< x<<" --icsd2proto < input.icsd \n\
  "<< x<<" --icsd2wyck < input.icsd \n\
  "<< x<<" --icsd2wyck --sof < input.icsd \n\
  "<< x<<" --icsdproto2aflowin < input.proto \n\
  \n\
 FROZSL \n\
  scripting for using the FROZSL phonon framework \n\
  "<< x<<" --frozsl_vaspsetup_aflow | --frozsl_vaspsetup < FROZSL.output \n\
  "<< x<<" --frozsl_vaspsetup_aflow --file \n\
  "<< x<<" --frozsl_vaspsetup_poscar < FROZSL.output \n\
  "<< x<<" --frozsl_vaspsetup_poscar --file \n\
  "<< x<<" --frozsl_analyze < aflow.frozsl.out \n\
  "<< x<<" --frozsl_input \n\
  "<< x<<" --frozsl_output \n\
  "<< x<<" --readme=frozsl \n\
  \n\
 ICSD/LIBN MODE (only for duke.edu computers) \n\
  scripting for the ICSD/LIBN libraries \n\
  "<< x<<" [--force] --lib2raw=FCC/La1Se1_ICSD_27104 \n\
  "<< x<<" [--force] --lib2raw=AgCdZr/T0001.A2BC \n\
  "<< x<<" [--force] --lib2raw=all[,directory] \n\
  \n\
 APENNSY MODE (only for duke.edu computers) \n\
  "<< x<<" --xfixX system structure \n\
  \n\
 CALCULATED DATABASE \n\
  "<< x<<" --calculated \n\
  "<< x<<" --calculated=icsd --random \n\
  "<< x<<" --calculated[[=]all|icsd|lib1|lib2|lib3|auro] \n\
   \n\
 CLUSTER EXPANSION (very primitive, only fcc and bcc) \n\
  "<< x<<" --cluster-expansion=... | --ce=structure_type,A,B,EA,EB \n\
  "<< x<<" --superlattice=structure_type,n_min,n_max < POSCAR \n\
  "<< x<<" --superlattice=VASP,structure_type,A,B < superlattice_name (CHECK) ??? \n\
  "<< x<<" --cluster=structure_type,n_min,n_max,m_min,m_max \n\
  "<< x<<" --special-quasirandom-structure=.. | --sqs=structure_type,atom_num,neighbour_num,sl_num_min,sl_num_max,A,B \n\
  "<< x<<" --special-quasirandom-structure=.. | --sqs=structure_type n1 n2 < POSCAR (CHECK) \n\
  \n\
 TERNARY CONVEX HULL (only for duke.edu computers) \n\
  "<< x<<" --terdata=A:B:C  [--fonts=XX | --keep=eps | --print=jpg | --print=gif | --print=png] \n\
  "<< x<<" --terdata_exist list \n\
  \n\
 COMPARE STRUCTURES \n\
  "<< x<<" Refer to aflow --readme=compare\n\
  "<< x<<" --compare_material=POSCAR1,POSCAR2 [--np=xx (default 8)] [--print]\n\
  "<< x<<" --compare_structure=POSCAR1,POSCAR2 [--np=xx (default 8)] [--print]\n\
  "<< x<<" --compare_material_directory|--compare_material_dir [-D \"PATH\"] [--np=xx (default 8)]\n\
  "<< x<<" --compare_structure_directory|--compare_structure_dir [-D \"PATH\"] [--np=xx (default 8)]\n\
  \n\
******* END POSTPROCESSING MODE ******************************************************************** \n\
  " << endl;
    // strstream << "Note: all the routines, except **, are tested to conform to convasp output" << endl;
    return strstream.str();
  }
}

/*
// [OBSOLETE] [JUNKAI] JUNKAI STUFF
// [OBSOLETE] [JUNKAI]   "<< x<<" --magpara label \n\
// [OBSOLETE] [JUNKAI]   "<< x<<" --magpara_h label \n\
// [OBSOLETE] [JUNKAI]   PROTOTYPE CHECKING (we are writing new codes so these instructions will disappear) \n\
// [OBSOLETE] [JUNKAI]   "<< x<<" --prototyper | --PROTOTYPER <POSCAR \n\
// [OBSOLETE] [JUNKAI]   "<< x<<" --prototyper_h | --PROTOTYPER_H < POSCAR \n\
// [OBSOLETE] [JUNKAI]   "<< x<<" --prototyper_m | --PROTOTYPER_M < POSCAR \n\
// [OBSOLETE] [JUNKAI]   "<< x<<" --prototypermult | --PROTOTYPERMULT < POSCAR \n\
// [OBSOLETE] [JUNKAI]   "<< x<<" --protoclassify < list \n\
// [OBSOLETE] [JUNKAI]   "<< x<<" --protoclassify name1,name2,name3 \n\
// [OBSOLETE] [JUNKAI]   "<< x<<" --protoclassify_h < list \n\
// [OBSOLETE] [JUNKAI]   "<< x<<" --protoclassify_m < list \n\
// [OBSOLETE] [JUNKAI]   \n\
*/

// ***************************************************************************
// pflow::ABCCAR
// ***************************************************************************
namespace pflow {
  xstructure ABCCAR(istream& input) {
    xstructure str(input,IOAFLOW_AUTO);
    str.iomode=IOVASP_ABCCAR;
    return str;
  }
} // namespace pflow

// ***************************************************************************
// pflow::ACE
// ***************************************************************************
namespace pflow {
  void ACE(istream& input) {
    //  xstructure a(input,IOAFLOW_AUTO);
    // pflow::PrintACE(a,cout);
    // [OBSOLETE]  pflow::PrintACE(xstructure(input,IOVASP_POSCAR),cout);
    pflow::PrintACE(xstructure(input,IOAFLOW_AUTO),cout);
  }
} // namespace pflow

//DX 9/27/17 - add spin info to xstructure - START
// ***************************************************************************
// pflow::AddSpinToXstructure (collinear version)
// ***************************************************************************
namespace pflow {
  bool AddSpinToXstructure(xstructure& a, vector<double>& vmag){
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(vmag.size()!=a.atoms.size()){
      cerr << "pflow::AddSpinToXstructure (collinear): ERROR: Number of magnetic moments (" << vmag.size() << ") does not match the number of atoms (" << a.atoms.size() << ")." << endl;
      return false;
    }
    // Only collinear for now 9/27/17
    for(uint i=0;i<a.atoms.size();i++){
      a.atoms[i].spin = vmag[i];
      a.atoms[i].spin_is_given = TRUE;
      if(LDEBUG){cerr << "pflow::AddSpinToXstructure (collinear): atom " << i << " magnetic moment: " << a.atoms[i].spin << endl;}
    }
    return true;
  }
}
//DX 9/27/17 - add spin info to xstructure - START

//DX 12/5/17 - add spin info to xstructure - START
// ***************************************************************************
// pflow::AddSpinToXstructure (non-collinear version)
// ***************************************************************************
namespace pflow {
  bool AddSpinToXstructure(xstructure& a, vector<xvector<double> >& vmag_noncoll){
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(vmag_noncoll.size()!=a.atoms.size()){
      cerr << "pflow::AddSpinToXstructure (non-collinear): ERROR: Number of magnetic moments (" << vmag_noncoll.size() << ") does not match the number of atoms (" << a.atoms.size() << ")." << endl;
      return false;
    }
    // Only collinear for now 9/27/17
    for(uint i=0;i<a.atoms.size();i++){
      a.atoms[i].noncoll_spin = vmag_noncoll[i];
      a.atoms[i].noncoll_spin_is_given = TRUE;
      if(LDEBUG){cerr << "pflow::AddSpinToXstructure (non-collinear): atom " << i << " magnetic moment: " << a.atoms[i].noncoll_spin << endl;}
    }
    return true;
  }
}
//DX 12/5/17 - add spin info to xstructure - START

// ***************************************************************************
// pflow::AFLOWIN
// ***************************************************************************
namespace pflow {
  string AFLOWIN(istream& input) {
    stringstream oss;
    xstructure str(input,IOAFLOW_AUTO);
    oss << AFLOWIN_SEPARATION_LINE << endl;
    oss << "[VASP_POSCAR_MODE_EXPLICIT]START" << endl;
    oss << str << "";
    oss << "[VASP_POSCAR_MODE_EXPLICIT]STOP" << endl;
    oss << AFLOWIN_SEPARATION_LINE << endl;
    return oss.str();
  }
} // namespace pflow

// ***************************************************************************
// pflow::POSCAR2AFLOWIN
// ***************************************************************************
namespace pflow {
  string POSCAR2AFLOWIN(istream& input) {
    stringstream oss;
    xstructure str(input,IOAFLOW_AUTO);
    _xvasp xvasp;
    AVASP_DefaultValuesBinary_AFLOWIN(xvasp);
    xvasp.AVASP_prototype_mode=LIBRARY_MODE_PROTOTYPE;
    xvasp.AVASP_flag_PRECISION_scheme="H";
    xvasp.str=str;
    stringstream aflowin;
    AVASP_MakeSingleAFLOWIN(xvasp,aflowin,(bool) FALSE,-1);
    oss << aflowin.str();
    return oss.str();
  }
} // namespace pflow

namespace pflow {
  bool SYMMETRY_GROUPS(_aflags &aflags,istream& input, aurostd::xoption& vpflow, ostream& oss){
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    _kflags kflags;                                   //DX 8/15/17 - Add in consistency checks
    xstructure _a(input,IOAFLOW_AUTO);
    bool osswrite=TRUE;
  
    char mode ='\0'; 
    string aliases = "";
    string sym_specific_options = "";
    string options = "";
    // AGROUP
    if(vpflow.flag("AGROUP")){
      mode = _AGROUP_;
      aliases = "--sitepointgroup|--agroup";
      sym_specific_options = "[--mag|--magnetic|--magmom=[m1,m2,...|INCAR|OUTCAR]]";
      options = vpflow.getattachedscheme("AGROUP");
      kflags.KBIN_SYMMETRY_CALCULATE_PGROUP=TRUE;       //DX 8/15/17 - Add in consistency checks
      kflags.KBIN_SYMMETRY_CALCULATE_PGROUPK=FALSE;     //DX 8/15/17 - Add in consistency checks
      kflags.KBIN_SYMMETRY_CALCULATE_FGROUP=TRUE;       //DX 8/15/17 - Add in consistency checks
      kflags.KBIN_SYMMETRY_CALCULATE_PGROUP_XTAL=TRUE;  //DX 8/15/17 - Add in consistency checks
      kflags.KBIN_SYMMETRY_CALCULATE_PGROUPK_XTAL=FALSE;//DX 8/15/17 - Add in consistency checks //DX 12/5/17 - Added pgroupk_xtal
      kflags.KBIN_SYMMETRY_CALCULATE_SGROUP=FALSE;      //DX 8/15/17 - Add in consistency checks
      kflags.KBIN_SYMMETRY_CALCULATE_IATOMS=TRUE;       //DX 8/15/17 - Add in consistency checks
      kflags.KBIN_SYMMETRY_CALCULATE_AGROUP=TRUE;       //DX 8/15/17 - Add in consistency checks
      kflags.KBIN_SYMMETRY_PGROUP_WRITE=TRUE;
      kflags.KBIN_SYMMETRY_PGROUPK_WRITE=FALSE;
      kflags.KBIN_SYMMETRY_FGROUP_WRITE=TRUE;
      kflags.KBIN_SYMMETRY_PGROUP_XTAL_WRITE=TRUE;
      kflags.KBIN_SYMMETRY_PGROUPK_XTAL_WRITE=FALSE;
      kflags.KBIN_SYMMETRY_SGROUP_WRITE=FALSE;
      kflags.KBIN_SYMMETRY_IATOMS_WRITE=TRUE;
      kflags.KBIN_SYMMETRY_AGROUP_WRITE=TRUE;
    }
    // FGROUP
    else if(vpflow.flag("FGROUP")){
      mode = _FGROUP_;
      aliases = "--factorgroup|--fgroup";
      sym_specific_options = "[--mag|--magnetic|--magmom=[m1,m2,...|INCAR|OUTCAR]]";
      options = vpflow.getattachedscheme("FGROUP");
      kflags.KBIN_SYMMETRY_CALCULATE_PGROUP=TRUE;       //DX 8/15/17 - Add in consistency checks
      kflags.KBIN_SYMMETRY_CALCULATE_PGROUPK=FALSE;     //DX 8/15/17 - Add in consistency checks
      kflags.KBIN_SYMMETRY_CALCULATE_FGROUP=TRUE;       //DX 8/15/17 - Add in consistency checks
      kflags.KBIN_SYMMETRY_CALCULATE_PGROUP_XTAL=FALSE;  //DX 8/15/17 - Add in consistency checks
      kflags.KBIN_SYMMETRY_CALCULATE_PGROUPK_XTAL=FALSE;//DX 8/15/17 - Add in consistency checks //DX 12/5/17 - Added pgroupk_xtal
      kflags.KBIN_SYMMETRY_CALCULATE_SGROUP=FALSE;      //DX 8/15/17 - Add in consistency checks
      kflags.KBIN_SYMMETRY_CALCULATE_IATOMS=FALSE;       //DX 8/15/17 - Add in consistency checks
      kflags.KBIN_SYMMETRY_CALCULATE_AGROUP=FALSE;       //DX 8/15/17 - Add in consistency checks
      kflags.KBIN_SYMMETRY_PGROUP_WRITE=TRUE;
      kflags.KBIN_SYMMETRY_PGROUPK_WRITE=FALSE;
      kflags.KBIN_SYMMETRY_FGROUP_WRITE=TRUE;
      kflags.KBIN_SYMMETRY_PGROUP_XTAL_WRITE=FALSE;
      kflags.KBIN_SYMMETRY_PGROUPK_XTAL_WRITE=FALSE;
      kflags.KBIN_SYMMETRY_SGROUP_WRITE=FALSE;
      kflags.KBIN_SYMMETRY_IATOMS_WRITE=FALSE;
      kflags.KBIN_SYMMETRY_AGROUP_WRITE=FALSE;
    }
    // PGROUP
    else if(vpflow.flag("PGROUP")){
      mode = _PGROUP_;
      aliases = "--pointgroup|--pgroup";
      options = vpflow.getattachedscheme("PGROUP");
      kflags.KBIN_SYMMETRY_CALCULATE_PGROUP=TRUE;       //DX 8/15/17 - Add in consistency checks
      kflags.KBIN_SYMMETRY_CALCULATE_PGROUPK=FALSE;     //DX 8/15/17 - Add in consistency checks
      kflags.KBIN_SYMMETRY_CALCULATE_FGROUP=FALSE;       //DX 8/15/17 - Add in consistency checks
      kflags.KBIN_SYMMETRY_CALCULATE_PGROUP_XTAL=FALSE;  //DX 8/15/17 - Add in consistency checks
      kflags.KBIN_SYMMETRY_CALCULATE_PGROUPK_XTAL=FALSE;//DX 8/15/17 - Add in consistency checks //DX 12/5/17 - Added pgroupk_xtal
      kflags.KBIN_SYMMETRY_CALCULATE_SGROUP=FALSE;      //DX 8/15/17 - Add in consistency checks
      kflags.KBIN_SYMMETRY_CALCULATE_IATOMS=FALSE;       //DX 8/15/17 - Add in consistency checks
      kflags.KBIN_SYMMETRY_CALCULATE_AGROUP=FALSE;       //DX 8/15/17 - Add in consistency checks
      kflags.KBIN_SYMMETRY_PGROUP_WRITE=TRUE;
      kflags.KBIN_SYMMETRY_PGROUPK_WRITE=FALSE;
      kflags.KBIN_SYMMETRY_FGROUP_WRITE=FALSE;
      kflags.KBIN_SYMMETRY_PGROUP_XTAL_WRITE=FALSE;
      kflags.KBIN_SYMMETRY_PGROUPK_XTAL_WRITE=FALSE;
      kflags.KBIN_SYMMETRY_SGROUP_WRITE=FALSE;
      kflags.KBIN_SYMMETRY_IATOMS_WRITE=FALSE;
      kflags.KBIN_SYMMETRY_AGROUP_WRITE=FALSE;
    }
    // PGROUPK
    else if(vpflow.flag("PGROUPK")){
      mode = _PGROUPK_;
      aliases = "--pointgroupklattice|--pgroupk";
      options = vpflow.getattachedscheme("PGROUPK");
      kflags.KBIN_SYMMETRY_CALCULATE_PGROUP=TRUE;       //DX 8/15/17 - Add in consistency checks
      kflags.KBIN_SYMMETRY_CALCULATE_PGROUPK=TRUE;     //DX 8/15/17 - Add in consistency checks
      kflags.KBIN_SYMMETRY_CALCULATE_FGROUP=FALSE;       //DX 8/15/17 - Add in consistency checks
      kflags.KBIN_SYMMETRY_CALCULATE_PGROUP_XTAL=FALSE;  //DX 8/15/17 - Add in consistency checks
      kflags.KBIN_SYMMETRY_CALCULATE_PGROUPK_XTAL=FALSE;//DX 8/15/17 - Add in consistency checks //DX 12/5/17 - Added pgroupk_xtal
      kflags.KBIN_SYMMETRY_CALCULATE_SGROUP=FALSE;      //DX 8/15/17 - Add in consistency checks
      kflags.KBIN_SYMMETRY_CALCULATE_IATOMS=FALSE;       //DX 8/15/17 - Add in consistency checks
      kflags.KBIN_SYMMETRY_CALCULATE_AGROUP=FALSE;       //DX 8/15/17 - Add in consistency checks
      kflags.KBIN_SYMMETRY_PGROUP_WRITE=TRUE;
      kflags.KBIN_SYMMETRY_PGROUPK_WRITE=TRUE;
      kflags.KBIN_SYMMETRY_FGROUP_WRITE=FALSE;
      kflags.KBIN_SYMMETRY_PGROUP_XTAL_WRITE=FALSE;
      kflags.KBIN_SYMMETRY_PGROUPK_XTAL_WRITE=FALSE;
      kflags.KBIN_SYMMETRY_SGROUP_WRITE=FALSE;
      kflags.KBIN_SYMMETRY_IATOMS_WRITE=FALSE;
      kflags.KBIN_SYMMETRY_AGROUP_WRITE=FALSE;
    }
    // PGROUPX
    else if(vpflow.flag("PGROUPX")){
      mode = _PGROUP_XTAL_;
      aliases = "--pointgroup_crystal|--pgroup_crystal|--pgroup_xtal|--pgroupx|--pgroupX";
      sym_specific_options = "[--mag|--magnetic|--magmom=[m1,m2,...|INCAR|OUTCAR]]";
      options = vpflow.getattachedscheme("PGROUPX");
      kflags.KBIN_SYMMETRY_CALCULATE_PGROUP=TRUE;       //DX 8/15/17 - Add in consistency checks
      kflags.KBIN_SYMMETRY_CALCULATE_PGROUPK=FALSE;     //DX 8/15/17 - Add in consistency checks
      kflags.KBIN_SYMMETRY_CALCULATE_FGROUP=TRUE;       //DX 8/15/17 - Add in consistency checks
      kflags.KBIN_SYMMETRY_CALCULATE_PGROUP_XTAL=TRUE;  //DX 8/15/17 - Add in consistency checks
      kflags.KBIN_SYMMETRY_CALCULATE_PGROUPK_XTAL=FALSE;//DX 8/15/17 - Add in consistency checks //DX 12/5/17 - Added pgroupk_xtal
      kflags.KBIN_SYMMETRY_CALCULATE_SGROUP=FALSE;      //DX 8/15/17 - Add in consistency checks
      kflags.KBIN_SYMMETRY_CALCULATE_IATOMS=FALSE;       //DX 8/15/17 - Add in consistency checks
      kflags.KBIN_SYMMETRY_CALCULATE_AGROUP=FALSE;       //DX 8/15/17 - Add in consistency checks
      kflags.KBIN_SYMMETRY_PGROUP_WRITE=TRUE;
      kflags.KBIN_SYMMETRY_PGROUPK_WRITE=FALSE;
      kflags.KBIN_SYMMETRY_FGROUP_WRITE=TRUE;
      kflags.KBIN_SYMMETRY_PGROUP_XTAL_WRITE=TRUE;
      kflags.KBIN_SYMMETRY_PGROUPK_XTAL_WRITE=FALSE;
      kflags.KBIN_SYMMETRY_SGROUP_WRITE=FALSE;
      kflags.KBIN_SYMMETRY_IATOMS_WRITE=FALSE;
      kflags.KBIN_SYMMETRY_AGROUP_WRITE=FALSE;
    }
    // PGROUPK_XTAL
    else if(vpflow.flag("PGROUPK_XTAL")){
      mode = _PGROUPK_XTAL_;
      aliases = "--pointgroupkcrystal|--pgroupk_xtal";
      sym_specific_options = "[--mag|--magnetic|--magmom=[m1,m2,...|INCAR|OUTCAR]]";
      options = vpflow.getattachedscheme("PGROUPK_XTAL");
      kflags.KBIN_SYMMETRY_CALCULATE_PGROUP=TRUE;       //DX 8/15/17 - Add in consistency checks
      kflags.KBIN_SYMMETRY_CALCULATE_PGROUPK=FALSE;     //DX 8/15/17 - Add in consistency checks
      kflags.KBIN_SYMMETRY_CALCULATE_FGROUP=TRUE;       //DX 8/15/17 - Add in consistency checks
      kflags.KBIN_SYMMETRY_CALCULATE_PGROUP_XTAL=TRUE;  //DX 8/15/17 - Add in consistency checks
      kflags.KBIN_SYMMETRY_CALCULATE_PGROUPK_XTAL=TRUE;//DX 8/15/17 - Add in consistency checks //DX 12/5/17 - Added pgroupk_xtal
      kflags.KBIN_SYMMETRY_CALCULATE_SGROUP=FALSE;      //DX 8/15/17 - Add in consistency checks
      kflags.KBIN_SYMMETRY_CALCULATE_IATOMS=FALSE;       //DX 8/15/17 - Add in consistency checks
      kflags.KBIN_SYMMETRY_CALCULATE_AGROUP=FALSE;       //DX 8/15/17 - Add in consistency checks
      kflags.KBIN_SYMMETRY_PGROUP_WRITE=TRUE;
      kflags.KBIN_SYMMETRY_PGROUPK_WRITE=FALSE;
      kflags.KBIN_SYMMETRY_FGROUP_WRITE=TRUE;
      kflags.KBIN_SYMMETRY_PGROUP_XTAL_WRITE=TRUE;
      kflags.KBIN_SYMMETRY_PGROUPK_XTAL_WRITE=TRUE;
      kflags.KBIN_SYMMETRY_SGROUP_WRITE=FALSE;
      kflags.KBIN_SYMMETRY_IATOMS_WRITE=FALSE;
      kflags.KBIN_SYMMETRY_AGROUP_WRITE=FALSE;
    }
    // SGROUP
    else if(vpflow.flag("SGROUP")){
      mode = _SGROUP_;
      aliases = "--spacegroup|--sgroup";
      sym_specific_options = "[--radius] [--mag|--magnetic|--magmom=[m1,m2,...|INCAR|OUTCAR]]";
      options = vpflow.getattachedscheme("SGROUP");
      kflags.KBIN_SYMMETRY_SGROUP_RADIUS = KBIN_SYMMETRY_SGROUP_RADIUS_DEFAULT;
      if(vpflow.flag("SYMMETRY::SGROUP_RADIUS")){
        kflags.KBIN_SYMMETRY_SGROUP_RADIUS = aurostd::string2utype<double>(vpflow.getattachedscheme("SYMMETRY::SGROUP_RADIUS"));
      }
      kflags.KBIN_SYMMETRY_CALCULATE_PGROUP=TRUE;       //DX 8/15/17 - Add in consistency checks
      kflags.KBIN_SYMMETRY_CALCULATE_PGROUPK=FALSE;     //DX 8/15/17 - Add in consistency checks
      kflags.KBIN_SYMMETRY_CALCULATE_FGROUP=TRUE;       //DX 8/15/17 - Add in consistency checks
      kflags.KBIN_SYMMETRY_CALCULATE_PGROUP_XTAL=FALSE;  //DX 8/15/17 - Add in consistency checks
      kflags.KBIN_SYMMETRY_CALCULATE_PGROUPK_XTAL=FALSE;//DX 8/15/17 - Add in consistency checks //DX 12/5/17 - Added pgroupk_xtal
      kflags.KBIN_SYMMETRY_CALCULATE_SGROUP=TRUE;      //DX 8/15/17 - Add in consistency checks
      kflags.KBIN_SYMMETRY_CALCULATE_IATOMS=FALSE;       //DX 8/15/17 - Add in consistency checks
      kflags.KBIN_SYMMETRY_CALCULATE_AGROUP=FALSE;       //DX 8/15/17 - Add in consistency checks
      kflags.KBIN_SYMMETRY_PGROUP_WRITE=TRUE;
      kflags.KBIN_SYMMETRY_PGROUPK_WRITE=FALSE;
      kflags.KBIN_SYMMETRY_FGROUP_WRITE=TRUE;
      kflags.KBIN_SYMMETRY_PGROUP_XTAL_WRITE=FALSE;
      kflags.KBIN_SYMMETRY_PGROUPK_XTAL_WRITE=FALSE;
      kflags.KBIN_SYMMETRY_SGROUP_WRITE=TRUE;
      kflags.KBIN_SYMMETRY_IATOMS_WRITE=FALSE;
      kflags.KBIN_SYMMETRY_AGROUP_WRITE=FALSE;
    }
    // EQUIVALENT / IATOMS (performed in "EQUIVALENT" function)

    vector<string> tokens;
    aurostd::string2tokens(options,tokens,",");
    if(tokens.size()==1) {
      if(tokens.at(0)=="usage" || tokens.at(0)=="USAGE") {
        return init::ErrorOption(cout,options,"pflow::SYMMETRY",
                          aurostd::liststring2string("aflow "+aliases+"[=tolerance| =tight| =loose] [--no_scan] [--print=txt| =json] [--screen_only] "+sym_specific_options+" < POSCAR  default: tolerance=(minimum_interatomic_distance)/100.0, print=txt"));
      }
    }
    //DX 170804 - need to rescale, so we make a fast copy and calculate
    xstructure a(_a);
    a.ReScale(1.0);
    //DX 2/21/18 [OBSOLETE] string directory=".";
    string directory=aurostd::execute2string("pwd"); //DX 2/21/18 - use pwd
    if(XHOST.vflag_control.flag("DIRECTORY")) {
      directory=XHOST.vflag_control.getattachedscheme("DIRECTORY");
    }
    aflags.Directory=directory;
    string print_directory = " [dir=" + a.directory + "]";
    if(!aurostd::FileExist(directory)) {
      oss << "ERROR: Unable to locate " << directory << "." << endl;
      oss << "Exiting." << endl;
      oss << endl;
      //return oss.str();
      return FALSE;
    }

    //DX 9/21/17 - MAGNETIC SYMMETRY - START
    if(vpflow.flag("SYMMETRY::MAGNETIC") && (mode == _AGROUP_ || mode == _FGROUP_ || mode == _PGROUP_XTAL_ || mode == _PGROUPK_XTAL_ || mode == _SGROUP_)){
	string magmom_info = vpflow.getattachedscheme("SYMMETRY::MAGNETIC");
        int num_atoms=a.atoms.size();
        bool is_noncoll=false; 
        vector<xvector<double> > vmag_noncoll;                                          //DX 12/5/17 - added non-collinear
        if(GetNonCollinearMagneticInfo(num_atoms,magmom_info,vmag_noncoll)){       //DX 12/5/17 - added non-collinear
          if(LDEBUG){cerr << "pflow::SYMMETRY: Non-collinear spin system detected." << endl;} //DX 12/5/17 - added non-collinear
          is_noncoll = true;
          if(!AddSpinToXstructure(a,vmag_noncoll)){
            exit(0);
          } 
        }                                                                               //DX 12/5/17 - added non-collinear
        bool is_coll=false; 
        vector<double> vmag;
        if(!is_noncoll){
          if(GetCollinearMagneticInfo(num_atoms,magmom_info,vmag)){
            if(LDEBUG){cerr << "pflow::SYMMETRY: Collinear spin system detected." << endl;} //DX 12/5/17 - added non-collinear
            is_coll = true;
        if(!AddSpinToXstructure(a,vmag)){
              exit(0);
        } 
          } 
        }
        if(!is_noncoll && !is_coll){                                                                                //DX 12/5/17 - added non-collinear
          cerr << "pflow::SYMMETRY: ERROR: Could not detect collinear or non-collinear spin(s). Check spin input." << endl; //DX 12/5/17 - added non-collinear
          exit(0);                                                                                                  //DX 12/5/17 - added non-collinear
        }
      } 
    //DX 9/21/17 - MAGNETIC SYMMETRY - END

    double default_tolerance=SYM::defaultTolerance(a);
    double tolerance = AUROSTD_NAN;
    if(vpflow.flag("SYMMETRY::TOLERANCE")){
      string tolerance_string = vpflow.getattachedscheme("SYMMETRY::TOLERANCE");
      if(tolerance_string[0] == 't' || tolerance_string[0] == 'T'){ //Tight
        tolerance=default_tolerance;
      }
      else if(tolerance_string[0] == 'l' || tolerance_string[0] == 'L'){ //Loose
        tolerance=default_tolerance*10.0;
      }
      else{
        tolerance=aurostd::string2utype<double>(vpflow.getattachedscheme("SYMMETRY::TOLERANCE"));
      }
    }
    else{
      tolerance = default_tolerance;
    }
    if(tolerance < 1e-10){
      cerr << "pflow::SYMMETRY ERROR: Tolerance cannot be zero (i.e. less than 1e-10) " << print_directory << "." << endl;
      return 0;
    }
    // DX 8/3/17 - Add format flag - START
    string format = "txt";
    if(XHOST.vflag_control.flag("PRINT_MODE::TXT")){
        format = "txt";
      }
    else if(XHOST.vflag_control.flag("PRINT_MODE::JSON")){
        format = "json";
    }
    else{ //default is txt
      format = "txt";
    }   
    bool print = false;
    if(vpflow.flag("SYMMETRY::SCREEN_ONLY")) {
      print=true;
      kflags.KBIN_SYMMETRY_PGROUP_WRITE=FALSE;
      kflags.KBIN_SYMMETRY_PGROUPK_WRITE=FALSE;
      kflags.KBIN_SYMMETRY_FGROUP_WRITE=FALSE;
      kflags.KBIN_SYMMETRY_PGROUP_XTAL_WRITE=FALSE;
      kflags.KBIN_SYMMETRY_PGROUPK_XTAL_WRITE=FALSE;
      kflags.KBIN_SYMMETRY_SGROUP_WRITE=FALSE;
      kflags.KBIN_SYMMETRY_IATOMS_WRITE=FALSE;
      kflags.KBIN_SYMMETRY_AGROUP_WRITE=FALSE;
      osswrite=FALSE;
    }  
    // Perform full scan
    bool no_scan = false;
    bool force_perform = true; //if no_scan fails, still return true at default tolerance (even though it cannot be validated)
    if(vpflow.flag("SYMMETRY::NO_SCAN")) {
      no_scan=true;
    }

    bool tobzip2 = TRUE;
    stringstream command;
    ofstream FileMESSAGE("/dev/null");
    //DX 2/21/18 [OBSOLETE] string directory=".";
    //DX 2/21/18 [OBSOLETE] if(XHOST.vflag_control.flag("DIRECTORY")) {
    //DX 2/21/18 [OBSOLETE]   directory=XHOST.vflag_control.getattachedscheme("DIRECTORY");
    //DX 2/21/18 [OBSOLETE] }
    //DX 2/21/18 [OBSOLETE] aflags.Directory=directory;
    //DX 2/21/18 [OBSOLETE] if(!aurostd::FileExist(directory)) {
    //DX 2/21/18 [OBSOLETE]  oss << "ERROR: Unable to locate " << directory << "." << endl;
    //DX 2/21/18 [OBSOLETE]  oss << "Exiting." << endl;
    //DX 2/21/18 [OBSOLETE]  oss << endl;
    //DX 2/21/18 [OBSOLETE]  //return oss.str();
    //DX 2/21/18 [OBSOLETE]  return FALSE;
    //DX 2/21/18 [OBSOLETE]}

    //while(symmetry_commensurate==FALSE){
    if(print == false){  //DX 8/3/17 - PRINT
      if(aurostd::FileExist(directory+string("/aflow.pgroup.out.bz2"))==TRUE ||
	 aurostd::FileExist(directory+string("/aflow.pgroup_xtal.out.bz2"))==TRUE ||
	 aurostd::FileExist(directory+string("/aflow.pgroupk.out.bz2"))==TRUE ||
	 aurostd::FileExist(directory+string("/aflow.pgroupk_xtal.out.bz2"))==TRUE || //DX 1/18/18 - added pgroupk_xtal
	 aurostd::FileExist(directory+string("/aflow.fgroup.out.bz2"))==TRUE ||
	 aurostd::FileExist(directory+string("/aflow.sgroup.out.bz2"))==TRUE ||
	 aurostd::FileExist(directory+string("/aflow.agroup.out.bz2"))==TRUE ||
	 aurostd::FileExist(directory+string("/aflow.iatoms.out.bz2"))==TRUE) {tobzip2=TRUE;}
      if(aurostd::FileExist(directory+string("/aflow.pgroup.out.bz2"))==TRUE || aurostd::FileExist(directory+string("/aflow.pgroup.out"))==TRUE){   command << "rm -f " << directory << "/aflow.pgroup.*" << endl;}
      if(aurostd::FileExist(directory+string("/aflow.pgroup_xtal.out.bz2"))==TRUE || aurostd::FileExist(directory+string("/aflow.pgroup_xtal.out"))==TRUE){   command << "rm -f " << directory << "/aflow.pgroup_xtal.*" << endl;}
      if(aurostd::FileExist(directory+string("/aflow.pgroupk.out.bz2"))==TRUE || aurostd::FileExist(directory+string("/aflow.pgroupk.out"))==TRUE){ command << "rm -f " << directory << "/aflow.pgroupk.*" << endl;}
      if(aurostd::FileExist(directory+string("/aflow.pgroupk_xtal.out.bz2"))==TRUE || aurostd::FileExist(directory+string("/aflow.pgroupk_xtal.out"))==TRUE){ command << "rm -f " << directory << "/aflow.pgroupk_xtal.*" << endl;} //DX 1/18/18 - added pgroupk_xtal
      if(aurostd::FileExist(directory+string("/aflow.fgroup.out.bz2"))==TRUE || aurostd::FileExist(directory+string("/aflow.fgroup.out"))==TRUE){ command << "rm -f " << directory << "/aflow.fgroup.*" << endl;}
      if(aurostd::FileExist(directory+string("/aflow.agroup.out.bz2"))==TRUE || aurostd::FileExist(directory+string("/aflow.agroup.out"))==TRUE){ command << "rm -f " << directory << "/aflow.agroup.*" << endl;}
      if(aurostd::FileExist(directory+string("/aflow.sgroup.out.bz2"))==TRUE || aurostd::FileExist(directory+string("/aflow.sgroup.out"))==TRUE){ command << "rm -f " << directory << "/aflow.sgroup.*" << endl;}
      if(aurostd::FileExist(directory+string("/aflow.iatoms.out.bz2"))==TRUE || aurostd::FileExist(directory+string("/aflow.iatoms.out"))==TRUE){ command << "rm -f " << directory << "/aflow.iatoms.*" << endl;}
	if(aurostd::FileExist(directory+string("/aflow.pgroup.json.bz2"))==TRUE ||
	   aurostd::FileExist(directory+string("/aflow.pgroup_xtal.json.bz2"))==TRUE ||
	   aurostd::FileExist(directory+string("/aflow.pgroupk.json.bz2"))==TRUE ||
	   aurostd::FileExist(directory+string("/aflow.pgroupk_xtal.json.bz2"))==TRUE || //DX 1/18/18 - added pgroupk_xtal
	   aurostd::FileExist(directory+string("/aflow.fgroup.json.bz2"))==TRUE ||
	   aurostd::FileExist(directory+string("/aflow.sgroup.json.bz2"))==TRUE ||
	   aurostd::FileExist(directory+string("/aflow.agroup.json.bz2"))==TRUE ||
	   aurostd::FileExist(directory+string("/aflow.iatoms.json.bz2"))==TRUE) {tobzip2=TRUE;}
	if(aurostd::FileExist(directory+string("/aflow.pgroup.json.bz2"))==TRUE || aurostd::FileExist(directory+string("/aflow.pgroup.json"))==TRUE){   command << "rm -f " << directory << "/aflow.pgroup.*" << endl;}
	if(aurostd::FileExist(directory+string("/aflow.pgroup_xtal.json.bz2"))==TRUE || aurostd::FileExist(directory+string("/aflow.pgroup_xtal.json"))==TRUE){   command << "rm -f " << directory << "/aflow.pgroup_xtal.*" << endl;}
	if(aurostd::FileExist(directory+string("/aflow.pgroupk.json.bz2"))==TRUE || aurostd::FileExist(directory+string("/aflow.pgroupk.json"))==TRUE){ command << "rm -f " << directory << "/aflow.pgroupk.*" << endl;}
	if(aurostd::FileExist(directory+string("/aflow.pgroupk_xtal.json.bz2"))==TRUE || aurostd::FileExist(directory+string("/aflow.pgroupk_xtal.json"))==TRUE){ command << "rm -f " << directory << "/aflow.pgroupk_xtal.*" << endl;} //DX 1/18/18 - added pgroupk_xtal
	if(aurostd::FileExist(directory+string("/aflow.fgroup.json.bz2"))==TRUE || aurostd::FileExist(directory+string("/aflow.fgroup.json"))==TRUE){ command << "rm -f " << directory << "/aflow.fgroup.*" << endl;}
	if(aurostd::FileExist(directory+string("/aflow.agroup.json.bz2"))==TRUE || aurostd::FileExist(directory+string("/aflow.agroup.json"))==TRUE){ command << "rm -f " << directory << "/aflow.agroup.*" << endl;}
	if(aurostd::FileExist(directory+string("/aflow.sgroup.json.bz2"))==TRUE || aurostd::FileExist(directory+string("/aflow.sgroup.json"))==TRUE){ command << "rm -f " << directory << "/aflow.sgroup.*" << endl;}
	if(aurostd::FileExist(directory+string("/aflow.iatoms.json.bz2"))==TRUE || aurostd::FileExist(directory+string("/aflow.iatoms.json"))==TRUE){ command << "rm -f " << directory << "/aflow.iatoms.*" << endl;}
      if(!command.str().empty()){
	aurostd::execute(command);
      }
    }
    if(!pflow::PerformFullSymmetry(a,tolerance,no_scan,force_perform,FileMESSAGE,aflags,kflags,osswrite,oss,format)){
      return FALSE;
    }
    // DX 8/3/17 - Print to symmetry operators to screen - START
    if(print==true){
      KBIN_SymmetryToScreen(a,format,oss,mode);
    }
    // DX 8/3/17 - Print to symmetry operators to screen - END
    // BZIP if necessary
    if(tobzip2) {
      command.str(std::string());
      command << "cd " << directory << endl;
      if(format == "txt"){ //DX 8/3/17 - FORMAT
      if(aurostd::FileExist(directory+string("/aflow.pgroup.out"))==TRUE) command << XHOST.command("bzip2") << " -9 aflow.pgroup.out" << endl;
      if(aurostd::FileExist(directory+string("/aflow.pgroup_xtal.out"))==TRUE) command << XHOST.command("bzip2") << " -9 aflow.pgroup_xtal.out" << endl;
      if(aurostd::FileExist(directory+string("/aflow.pgroupk.out"))==TRUE) command << XHOST.command("bzip2") << " -9 aflow.pgroupk.out" << endl;
      if(aurostd::FileExist(directory+string("/aflow.pgroupk_xtal.out"))==TRUE) command << XHOST.command("bzip2") << " -9 aflow.pgroupk_xtal.out" << endl; //DX 1/18/18 - added pgroupk_xtal
      if(aurostd::FileExist(directory+string("/aflow.fgroup.out"))==TRUE) command << XHOST.command("bzip2") << " -9 aflow.fgroup.out" << endl;
      if(aurostd::FileExist(directory+string("/aflow.sgroup.out"))==TRUE) command << XHOST.command("bzip2") << " -9 aflow.sgroup.out" << endl;
      if(aurostd::FileExist(directory+string("/aflow.agroup.out"))==TRUE) command << XHOST.command("bzip2") << " -9 aflow.agroup.out" << endl;
      if(aurostd::FileExist(directory+string("/aflow.iatoms.out"))==TRUE) command << XHOST.command("bzip2") << " -9 aflow.iatoms.out" << endl;
      }
      else if(format == "json"){ //DX 8/3/17 - FORMAT
        if(aurostd::FileExist(directory+string("/aflow.pgroup.json"))==TRUE) command << XHOST.command("bzip2") << " -9 aflow.pgroup.json" << endl;
        if(aurostd::FileExist(directory+string("/aflow.pgroup_xtal.json"))==TRUE) command << XHOST.command("bzip2") << " -9 aflow.pgroup_xtal.json" << endl;
        if(aurostd::FileExist(directory+string("/aflow.pgroupk.json"))==TRUE) command << XHOST.command("bzip2") << " -9 aflow.pgroupk.json" << endl;
        if(aurostd::FileExist(directory+string("/aflow.pgroupk_xtal.json"))==TRUE) command << XHOST.command("bzip2") << " -9 aflow.pgroupk_xtal.json" << endl; //DX 1/18/18 - added pgroupk_xtal
        if(aurostd::FileExist(directory+string("/aflow.fgroup.json"))==TRUE) command << XHOST.command("bzip2") << " -9 aflow.fgroup.json" << endl;
        if(aurostd::FileExist(directory+string("/aflow.sgroup.json"))==TRUE) command << XHOST.command("bzip2") << " -9 aflow.sgroup.json" << endl;
        if(aurostd::FileExist(directory+string("/aflow.agroup.json"))==TRUE) command << XHOST.command("bzip2") << " -9 aflow.agroup.json" << endl;
        if(aurostd::FileExist(directory+string("/aflow.iatoms.json"))==TRUE) command << XHOST.command("bzip2") << " -9 aflow.iatoms.json" << endl;
      }
      aurostd::execute(command);
    }
    //return oss.str(); //string("AFLOW "+string(AFLOW_VERSION)+" Symmetry Fixed in "+directory+"  (need aflow>=2948)\n");
    return true; //string("AFLOW "+string(AFLOW_VERSION)+" Symmetry Fixed in "+directory+"  (need aflow>=2948)\n");
  }
}

// ***************************************************************************
// pflow::AGROUP
// ***************************************************************************
namespace pflow {
  void AGROUP(_aflags &aflags,istream& input) {
    //  cout << aflow::Banner("BANNER_TINY") << endl;
    aflags.QUIET=TRUE;
    xstructure a(input,IOAFLOW_AUTO);
    bool WRITE=TRUE;
    ofstream File("/dev/null");
    // DX 8/15/17 - Add in consistency checks bool verbose=TRUE;
    // DX 8/15/17 - Add in consistency checks SYM::CalculatePointGroup(File,a,aflags,WRITE,verbose,cout);
    // DX 8/15/17 - Add in consistency checks SYM::CalculateFactorGroup(File,a,aflags,WRITE,verbose,cout);
    // DX 8/15/17 - Add in consistency checks SYM::CalculateSpaceGroup(File,a,aflags,FALSE,verbose,cout);
    // DX 8/15/17 - Add in consistency checks SYM::CalculateInequivalentAtoms(File,a,aflags,WRITE,verbose,cout);
    // DX 8/15/17 - Add in consistency checks SYM::CalculateSitePointGroup(File,a,aflags,WRITE,TRUE,cout);
    _kflags kflags;                                   //DX 8/15/17 - Add in consistency checks
    kflags.KBIN_SYMMETRY_CALCULATE_PGROUP=TRUE;       //DX 8/15/17 - Add in consistency checks
    kflags.KBIN_SYMMETRY_CALCULATE_PGROUPK=FALSE;     //DX 8/15/17 - Add in consistency checks
    kflags.KBIN_SYMMETRY_CALCULATE_FGROUP=TRUE;       //DX 8/15/17 - Add in consistency checks
    kflags.KBIN_SYMMETRY_CALCULATE_PGROUP_XTAL=TRUE;  //DX 8/15/17 - Add in consistency checks
    kflags.KBIN_SYMMETRY_CALCULATE_PGROUPK_XTAL=FALSE;//DX 8/15/17 - Add in consistency checks //DX 12/5/17 - Added pgroupk_xtal
    kflags.KBIN_SYMMETRY_CALCULATE_SGROUP=FALSE;      //DX 8/15/17 - Add in consistency checks
    kflags.KBIN_SYMMETRY_CALCULATE_IATOMS=TRUE;       //DX 8/15/17 - Add in consistency checks
    kflags.KBIN_SYMMETRY_CALCULATE_AGROUP=TRUE;       //DX 8/15/17 - Add in consistency checks
    pflow::PerformFullSymmetry(a,File,aflags,kflags,WRITE,cout); //DX 8/15/17 - Add in consistency checks
  }
} // namespace pflow

// ***************************************************************************
// pflow::AGROUP2
// ***************************************************************************
namespace pflow {
  void AGROUP2(istream& input) {
    xstructure a(input,IOAFLOW_AUTO);
    int Nbasis=a.atoms.size();
    vector<vector<uint > > Bsitesym(Nbasis);
    bool ComMidss=false; //calculate common site symmetry of A,B,and middle_AB?
    SYM::CalculateSitePointGroup2(a,ComMidss);
  }
} // namespace pflow

// ***************************************************************************
// pflow::AGROUP2m
// ***************************************************************************
namespace pflow {
  void AGROUP2m(istream& input) {
    xstructure a(input,IOAFLOW_AUTO);
    int Nbasis=a.atoms.size();
    vector<vector<uint > > Bsitesym(Nbasis);
    bool ComMidss=true; //calculate common site symmetry of A,B,and middle_AB?
    SYM::CalculateSitePointGroup2(a,ComMidss);
  }
} // namespace pflow

// ***************************************************************************
// pflow::ALPHABETIC
// ***************************************************************************
namespace pflow {
  xstructure ALPHABETIC(istream& input) {
    xstructure a(input,IOAFLOW_AUTO);
    //  cerr << "ALPHA=" << a.SpeciesGetAlphabetic() << endl;
    a.SpeciesPutAlphabetic();
    return a;
  }
} // namespace pflow

// ***************************************************************************
// pflow::ALPHACompound
// ***************************************************************************
namespace pflow {
  string ALPHACompound(string options) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "pflow::ALPHACompound: BEGIN" << endl;
    vector<string> tokens;
    aurostd::string2tokens(options,tokens,",");
    if(tokens.size()==0) {
      init::ErrorOption(cout,options,"pflow::ALPHACompound","aflow --alpha_compound=string1,string2,....");
      exit(0);
    } 
    // move on
    stringstream output;
    for(uint i=0;i<tokens.size();i++) {
      string system=string(tokens.at(i));
      vector<string> vsystem;vector<double> vnumber;
      XATOM_AlphabetizationCompound(system,vsystem,vnumber);
      output << system << endl;
    }
    if(LDEBUG) cerr << "pflow::ALPHACompound: END" << endl;
    return output.str();
  }
} // namespace pflow

// ***************************************************************************
// pflow::ALPHASpecies
// ***************************************************************************
namespace pflow {
  string ALPHASpecies(string options) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "pflow::ALPHASpecies: BEGIN" << endl;
    vector<string> tokens;
    aurostd::string2tokens(options,tokens,",");
    if(tokens.size()==0) {
      init::ErrorOption(cout,options,"pflow::ALPHASpecies","aflow --alpha_species=string1,string2,....");
      exit(0);
    } 
    // move on
    stringstream output;
    for(uint i=0;i<tokens.size();i++) {
      string system=string(tokens.at(i));
      vector<string> vsystem;vector<double> vnumber;
      XATOM_AlphabetizationSpecies(system,vsystem,vnumber);
      output << system << endl;
    }
    if(LDEBUG) cerr << "pflow::ALPHASpecies: END" << endl;
    return output.str();
  }
} // namespace pflow

// ***************************************************************************
// pflow::ANGLES
// ***************************************************************************
namespace pflow {
  void ANGLES(string options,istream& input) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "pflow::ANGLES: BEGIN" << endl;
    vector<string> tokens;
    aurostd::string2tokens(options,tokens,",");
    if(tokens.size()!=1) {
      init::ErrorOption(cout,options,"pflow::ANGLES","aflow --angle=cutoff < POSCAR");
      exit(0);
    } 
    // move on
    double cutoff=0.0; // some defaults
    if(tokens.size()>=1) cutoff=aurostd::string2utype<double>(tokens.at(0));
    
    xstructure a(input,IOAFLOW_AUTO);
    cout << aflow::Banner("BANNER_TINY") << endl;
    pflow::PrintAngles(a,cutoff,cout);
    if(LDEBUG) cerr << "pflow::ANGLES: END" << endl;
  }
} // namespace pflow

// ***************************************************************************
// pflow::ATOMSMAX
// ***************************************************************************
namespace pflow {
  string ATOMSMAX(string options,istream& input) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "pflow::ATOMSMAX: BEGIN" << endl;
    vector<string> tokens;
    aurostd::string2tokens(options,tokens,",");
    if(tokens.size()!=1) {
      init::ErrorOption(cout,options,"pflow::ATOMSMAX","aflow --maxatoms=N | --max_atoms=N | --atoms_max=N | --atomsmax=N < POSCAR");
      exit(0);
    } 
    // move on
    uint N=0; // some defaults
    if(tokens.size()>=1) N=aurostd::string2utype<uint>(tokens.at(0));

    xstructure a(input,IOAFLOW_AUTO);
    stringstream oss;
    if(a.atoms.size()>N) {
      oss << "MAX ATOMS SIZE = " << N << endl;
      oss << "your structure has " << a.atoms.size() << " atoms" << endl;
    } else {
      oss << a;
    }
    if(LDEBUG) cerr << "pflow::ATOMSMAX: END" << endl;
    return oss.str();
  }
} // namespace pflow

// ***************************************************************************
// pflow::BANDS
// ***************************************************************************
namespace pflow {
  void BANDS(string options,istream& input) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "pflow::BANDS: BEGIN" << endl;
    vector<string> tokens;
    aurostd::string2tokens(options,tokens,",");
    if(tokens.size()!=1) {
      init::ErrorOption(cout,options,"pflow::BANDS","aflow --bands=PROOUT < POSCAR");
      exit(0);
    } 
    // move on
    string filename=""; // some defaults
    if(tokens.size()>=1) filename=(tokens.at(0));

    //  cout << aflow::Banner("BANNER_TINY") << endl;
    xstructure str(input,IOAFLOW_AUTO);
    pflow::projdata prd;
    prd.PROOUTinfile=filename;
    pflow::ReadInProj(prd);
    prd.lat=pflow::GetScaledLat(str);
    pflow::PrintBands(prd);
    if(LDEBUG) cerr << "pflow::BANDS: END" << endl;
  }
} // namespace pflow

// ***************************************************************************
// pflow::BANDGAP // CAMILO
// ***************************************************************************
namespace pflow {
  void BANDGAP(aurostd::xoption& vpflow,ostream& oss) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "pflow::BANDGAP: BEGIN" << endl;
    string input=aurostd::RemoveWhiteSpaces(vpflow.getattachedscheme("BANDGAP"));
    if(input.empty()){
      init::ErrorOption(cout,vpflow.getattachedscheme("BANDGAP"),"pflow::BANDGAP","aflow --bandgap[=bands_directory1[,bands_directory2,...]]");
      exit(0);
    } 
    vector<string> dirs;
    aurostd::string2tokens(input,dirs,",");
    for(uint i=0;i<dirs.size();i++){
      if(!PrintBandGap(dirs[i], oss)){oss << "pflow::BANDGAP: "+dirs[i]+" failed" << endl;}
    }
    if(LDEBUG) cerr << "pflow::BANDGAP: END" << endl;
  }
} // namespace pflow

// ***************************************************************************
// pflow::BANDSTRUCTURE
// ***************************************************************************
namespace pflow {
  void BANDSTRUCTURE(_aflags &aflags) {
    //  cout << aflow::Banner("BANNER_TINY") << endl;
    aflags.QUIET=TRUE;
    string directory_LIB,directory_RAW;
    directory_LIB=aflags.Directory;
    directory_RAW=directory_LIB+"./BANDS";
    aurostd::execute("rm -rf "+directory_RAW);
    aurostd::DirectoryMake(directory_RAW);
    stringstream aflowlib_out; // will not be used
    vector<string> vspecies,vfiles;
    aflowlib::_aflowlib_entry data;
    aflowlib::GetSpeciesDirectory(directory_LIB,vspecies);
    for(uint i=0;i<vspecies.size();i++) {
      data.species+=vspecies.at(i);if(i<vspecies.size()-1) data.species+=",";
    }
    // [OBSOLETE]  aflowlib::LIB2RAW_Loop_Bands(directory_LIB,directory_RAW,vfiles,aflowlib_out,data,"pflow::BANDSTRUCTURE");
    aflowlib::LIB2RAW_Loop_Bands(directory_LIB,directory_RAW,vfiles,data,"pflow::BANDSTRUCTURE");
  }
} // namespace pflow

// ***************************************************************************
// pflow::BZDirectionsLATTICE
// ***************************************************************************
namespace pflow {
  string BZDirectionsLATTICE(string options) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "pflow::BZDIRECTIONSLATTICE: BEGIN" << endl;
    vector<string> tokens;
    aurostd::string2tokens(options,tokens,",");
    if(tokens.size()!=1) {
      init::ErrorOption(cout,options,"pflow::BZDIRECTIONSLATTICE","aflow --bzdirections= | --bzd=LATTICE");
      exit(0);
    }  else {     // move on
      if(tokens.at(0)!="TRI1a" && tokens.at(0)!="TRI1b" && tokens.at(0)!="TRI2a" && tokens.at(0)!="TRI2b" &&
	 tokens.at(0)!="MCL" &&
	 tokens.at(0)!="MCLC1" && tokens.at(0)!="MCLC2" && tokens.at(0)!="MCLC3" && tokens.at(0)!="MCLC4" && tokens.at(0)!="MCLC5" &&
	 tokens.at(0)!="ORC" &&
	 tokens.at(0)!="ORCC" &&
	 tokens.at(0)!="ORCF1" && tokens.at(0)!="ORCF2" && tokens.at(0)!="ORCF3" &&
	 tokens.at(0)!="ORCI" &&
	 tokens.at(0)!="TET" &&
	 tokens.at(0)!="BCT1" && tokens.at(0)!="BCT2" &&
	 tokens.at(0)!="RHL1" && tokens.at(0)!="RHL2" &&
	 tokens.at(0)!="HEX" &&
	 tokens.at(0)!="CUB" &&
	 tokens.at(0)!="FCC" &&
	 tokens.at(0)!="BCC"	 ) {
	init::ErrorOption(cout,options,"pflow::BZDIRECTIONSLATTICE","aflow --bzdirections= | --bzd=LATTICE");
	cout << "LATTICES = (the ones with \"\") " << endl;
	cout << "1. TRI order: kalpha,kbeta,kgamma  > 90 (kgamma<kalpha, kgamma<kbeta) " << endl;
	cout << "   or kalpha,kbeta,kgamma  < 90 (kgamma>kalpha, kgamma>kbeta) " << endl;
	cout << "   special case when kgamma=90 " << endl;
	cout << "   \"TRI1a\" kalpha>90 kbeta>90 kgamma>90 " << endl;
	cout << "   \"TRI1b\" kalpha<90 kbeta<90 kgamma<90 " << endl;
	cout << "   \"TRI2a\" kalpha>90 kbeta>90 kgamma=90 " << endl;
	cout << "   \"TRI2b\" kalpha<90 kbeta<90 kgamma=90 " << endl;
	cout << "2. \"MCL\" unique (order b<=c) " << endl;
	cout << "3. MCLC (order alpha<90) " << endl;
	cout << "   \"MCLC1\"  kgamma>90 " << endl;
	cout << "   \"MCLC2\"  kgamma=90 " << endl;
	cout << "   \"MCLC3\"  kgamma<90, b*cos(alpha)/c + (b*sin(alpha)/a)^2 < 1 " << endl;
	cout << "   \"MCLC4\"  kgamma<90, b*cos(alpha)/c + (b*sin(alpha)/a)^2 = 1 " << endl;
	cout << "   \"MCLC5\"  kgamma<90, b*cos(alpha)/c + (b*sin(alpha)/a)^2 > 1 " << endl;
	cout << "4. \"ORC\" unique (order a<b<c) " << endl;
	cout << "5. \"ORCC\" unique (order a<b) " << endl;
	cout << "6. ORCF (order a<b<c) " << endl;
	cout << "   \"ORCF1\" \"ORCF_invb2+invc2<inva2\"  for 1/a^2 > 1/b^2 + 1/c^2 " << endl;
	cout << "   \"ORCF2\" \"ORCF_inva2<invb2+invc2\"  for 1/a^2 < 1/b^2 + 1/c^2 " << endl;
	cout << "   \"ORCF3\"                           for 1/a^2 = 1/b^2 + 1/c^2 " << endl;
	cout << "7. \"ORCI\" unique (order a<b<c) " << endl;
	cout << "8. \"TET\" unique (order a a c) " << endl;
	cout << "9. BCT (order a a c) " << endl;
	cout << "   \"BCT1\" \"BCT_c<a\" for c<a " << endl;
	cout << "   \"BCT2\" \"BCT_c>a\" for c>a " << endl;
	cout << "10. \"RHL1\" alpha<90 " << endl;
	cout << "    \"RHL2\" alpha>90 " << endl;
	cout << "11. \"HEX\" unique (order 60 90 90) " << endl;
	cout << "12. \"CUB\" unique " << endl;
	cout << "13. \"FCC\" unique (order 60 60 60) " << endl;
	cout << "14. \"BCC\" unique " << endl;
	exit(0);
      }
    }
    
    double grid=20.0;
    xmatrix<double> rlattice(3,3);
    bool foundBZ;
    rlattice[1][1]=1.0;rlattice[2][2]=1.0;rlattice[3][3]=1.0;
    if(LDEBUG) cerr << "pflow::BZDIRECTIONSLATTICE: END" << endl;
    return LATTICE::KPOINTS_Directions(tokens.at(0),rlattice,grid,IOVASP_AUTO,foundBZ);
  }
} // namespace pflow

// ***************************************************************************
// pflow::BZDirectionsSTRUCTURE
// ***************************************************************************
namespace pflow {
  string BZDirectionsSTRUCTURE(istream& input) {
    xstructure a(input,IOAFLOW_AUTO);
    bool foundBZ;
    double grid=20.0;
    a.GetLatticeType();
    cerr << "LATTICE INFORMATION" << endl;
    cerr << " Real space lattice primitive           = " << a.bravais_lattice_type << endl;
    cerr << " Real space lattice variation           = " << a.bravais_lattice_variation_type << endl; //wahyu mod
    //  cerr << " Real space conventional lattice        = " << a.bravais_conventional_lattice_type << endl;
    cerr << " Real space Pearson symbol              = " << a.pearson_symbol << endl;
    cerr << " Reciprocal lattice primitive           = " << a.reciprocal_lattice_type << endl;
    cerr << " Reciprocal lattice variation           = " << a.reciprocal_lattice_variation_type << endl;//wahyu mod
    //  cerr << " Reciprocal conventional lattice        = " << a.reciprocal_conventional_lattice_type << endl;
    cerr << "KPOINTS FILE" << endl;
    return LATTICE::KPOINTS_Directions(a.bravais_lattice_variation_type,a.lattice,grid,a.iomode,foundBZ);
  }
} // namespace pflow

// ***************************************************************************
// pflow::CAGES
// ***************************************************************************
namespace pflow {
  void CAGES(_aflags &aflags,string options,istream& input) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "pflow::CAGES: BEGIN" << endl;
    vector<string> tokens;
    aurostd::string2tokens(options,tokens,",");
    if(tokens.size()!=0 && tokens.size()!=1) {
      init::ErrorOption(cout,options,"pflow::CAGES","aflow [--np=NP] --cages[=roughness] < POSCAR");
      exit(0);
    }
    double roughness=-1.0;
    if(tokens.size()>=1) roughness=aurostd::string2utype<double>(tokens.at(0));
    xstructure a(input,IOAFLOW_AUTO);
    vector<acage> cagesirreducible,cagesreducible,cages4,cages3,cages2;
    GetCages(a,aflags,cagesirreducible,cagesreducible,cages4,cages3,cages2,roughness,TRUE,cout);
    //  cout << "REDUCIBLE_SIZE " << cagesreducible.size() << endl;
    // cout << "IRREDUCIBLE_SIZE " << cagesirreducible.size() << endl;
    if(LDEBUG) cerr << "pflow::CAGES: END" << endl;
  }
} // namespace pflow

// ***************************************************************************
// pflow::CART
// ***************************************************************************
namespace pflow {
  xstructure CART(istream& input) {
    xstructure a(input,IOAFLOW_AUTO);
    xstructure b(a);
    b.SetCoordinates(_COORDS_CARTESIAN_);
    return b;
  }
} // namespace pflow

// ***************************************************************************
// pflow::ChangeSuffix(string options)
// ***************************************************************************
namespace pflow {
  void ChangeSuffix(string _operation) {
    // aflow --suffix=[directory,]"from2to"
    //  Change the suffixes of VASP files. Easy conversion between AFLOW format and VASP format. (KESONG Dec. 22nd, 2013)
    //	Mnemonic: from2to with from/to =[n=none;r1=relax1;r2=relax2;r3=relax3;s=static;b=bands] or without abbreviations.
    vector<string> vfile,tokens;
    aurostd::string2tokens("AECCAR0,AECCAR1,AECCAR2,INCAR,POSCAR,POTCAR,KPOINTS,CONTCAR,OUTCAR,EIGENVAL,DOSCAR,IBZKPT,OSZICAR,PCDAT,XDATCAR,CHG,CHGCAR,PROCAR,ELFCAR,WAVECAR,LOCPOT,vasprun.xml,vasp.out",vfile,",");
    aurostd::string2tokens(_operation,tokens,",");
    string directory,operation;
    if(tokens.size()==0) {directory="";operation=_operation;}
    if(tokens.size()==1) {directory="";operation=tokens.at(0);}
    if(tokens.size()==2) {directory=tokens.at(0)+"/";operation=tokens.at(1);}
    if(operation!="") {
      for(uint i=0; i<vfile.size();i++) {
	string from="",to="",f=aurostd::CleanFileName(directory+vfile.at(i));
	if(aurostd::FileExist(f+".bz2")) {aurostd::execute(XHOST.command("bunzip2")+" "+f+".bz2");}
	if(aurostd::FileExist(f+".relax1.bz2")) {aurostd::execute(XHOST.command("bunzip2")+" "+f+".relax1.bz2");}
	if(aurostd::FileExist(f+".relax2.bz2")) {aurostd::execute(XHOST.command("bunzip2")+" "+f+".relax2.bz2");}
	if(aurostd::FileExist(f+".relax3.bz2")) {aurostd::execute(XHOST.command("bunzip2")+" "+f+".relax3.bz2");}
	if(aurostd::FileExist(f+".static.bz2")) {aurostd::execute(XHOST.command("bunzip2")+" "+f+".static.bz2");}
	if(aurostd::FileExist(f+".bands.bz2")) {aurostd::execute(XHOST.command("bunzip2")+" "+f+".bands.bz2");}
	if(aurostd::substring2bool(operation,"n2") || aurostd::substring2bool(operation,"none2") || aurostd::substring2bool(operation,"normal2")) from="";
	if(aurostd::substring2bool(operation,"r12") || aurostd::substring2bool(operation,"relax12")) from=".relax1";
	if(aurostd::substring2bool(operation,"r22") || aurostd::substring2bool(operation,"relax22")) from=".relax2";
	if(aurostd::substring2bool(operation,"r32") || aurostd::substring2bool(operation,"relax32")) from=".relax3";
	if(aurostd::substring2bool(operation,"s2") || aurostd::substring2bool(operation,"static2")) from=".static";
	if(aurostd::substring2bool(operation,"b2") || aurostd::substring2bool(operation,"bands2")) from=".bands";
	if(aurostd::substring2bool(operation,"2n") || aurostd::substring2bool(operation,"2none") || aurostd::substring2bool(operation,"2normal")) to="";
	if(aurostd::substring2bool(operation,"2r1") || aurostd::substring2bool(operation,"2relax1")) to=".relax1";
	if(aurostd::substring2bool(operation,"2r2") || aurostd::substring2bool(operation,"2relax2")) to=".relax2";
	if(aurostd::substring2bool(operation,"2r3") || aurostd::substring2bool(operation,"2relax3")) to=".relax3";
	if(aurostd::substring2bool(operation,"2s") || aurostd::substring2bool(operation,"2static")) to=".static";
	if(aurostd::substring2bool(operation,"2b") || aurostd::substring2bool(operation,"2bands")) to=".bands";
	
	if(from!=""||to!="") {
	  if(aurostd::FileExist(f+from)) {
	    cout << "pflow::ChangeSuffix: mv " << f << from << " " << f << to << endl;
	    aurostd::execute(XHOST.command("mv")+" "+f+from+" "+f+to);
	  }
	}
      }
    }
  }
}  // namespace pflow

// ***************************************************************************
// pflow::CheckIntegritiy
// ***************************************************************************
// Shidong Wang 2011
namespace pflow {
  void CheckIntegritiy(void) {
    // Check whether function `isequal', `isdifferent', `identical' etc
    // give correct answer or not
    // issue of g++ optimization (-O -O1 -O2 -O3).

    // xvector operators

    const int test_xv_size = 5;
    const int _max_int = 100;
    const int _tol_int = 0;
    const float _tol_float = 1.0e-6;
    const double _tol_double = 1.0e-6;

    xvector<int> test_int_xv_a(test_xv_size), test_int_xv_b(test_xv_size);
    xvector<float> test_float_xv_a(test_xv_size), test_float_xv_b(test_xv_size);
    xvector<double> test_double_xv_a(test_xv_size), test_double_xv_b(test_xv_size);

    // initialze random seed
    srand(time(NULL));

    // initializae test xvectors
    for(int i=1; i<test_xv_size+1; i++) {
      test_int_xv_a[i] =  rand() % _max_int;
      test_int_xv_b[i] =  rand() % _max_int;
      test_float_xv_a[i] = PI * float((rand()%_max_int));
      test_float_xv_b[i] = PI * float((rand()%_max_int));
      test_double_xv_a[i] = _mm_e * float((rand()%_max_int));
      test_double_xv_b[i] = _mm_e * float((rand()%_max_int));
    }

    const int line_width = 60;
    for(int i=0; i<line_width; i++) { cout << "*"; } cout << endl;
    cout << "Shidong Wang - 2011 " << endl;
    for(int i=0; i<line_width; i++) { cout << "*";} cout << endl;
    cout << "xvector\n";

    // template<class utype> bool
    // identical(const xvector<utype>&,const xvector<utype>&,const utype&) __xprototype;
    //cout << "Indentical Integer with given tolerance " << _tol_int << " : \n"

    if(identical(test_int_xv_a, test_int_xv_a, _tol_int) &&
       (! identical(test_int_xv_a, test_int_xv_b, _tol_int))) {
      cout << ">>> Testing Integer INDENTICAL with given tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Integer INDENTICAL with given tolerance:  ***FAIL***!\n"
	   << " 1) \n"
	   << test_int_xv_a << endl
	   << " and \n"
	   << test_int_xv_a  << endl
	   << " are identical? "
	   << std::boolalpha << identical(test_int_xv_a, test_int_xv_a, _tol_int) << endl;
      cout << "2) \n"
	   << test_int_xv_a << endl
	   << " and \n"
	   << test_int_xv_b  << endl
	   << " are not identical? "
	   << std::boolalpha << (! identical(test_int_xv_a, test_int_xv_b, _tol_int)) << endl;
    }

    if(identical(test_float_xv_a, test_float_xv_a, _tol_float) &&
       (! identical(test_float_xv_a, test_float_xv_b, _tol_float))) {
      cout << ">>> Testing Float INDENTICAL with given tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Float INDENTICAL with given tolerance:  ***FAIL***!\n"
	   << " 1) \n"
	   << test_float_xv_a << endl
	   << " and \n"
	   << test_float_xv_a  << endl
	   << " are identical? "
	   << std::boolalpha << identical(test_float_xv_a, test_float_xv_a, _tol_float) << endl;
      cout << "2) \n"
	   << test_float_xv_a << endl
	   << " and \n"
	   << test_float_xv_b  << endl
	   << " are not identical? "
	   << std::boolalpha << (! identical(test_float_xv_a, test_float_xv_b, _tol_float)) << endl;
    }

    if(identical(test_double_xv_a, test_double_xv_a, _tol_double) &&
       (! identical(test_double_xv_a, test_double_xv_b, _tol_double))) {
      cout << ">>> Testing Double INDENTICAL with default tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Double INDENTICAL with default tolerance:  ***FAIL***!\n"
	   << " 1) \n"
	   << test_double_xv_a << endl
	   << " and \n"
	   << test_double_xv_a  << endl
	   << " are identical? "
	   << std::boolalpha << identical(test_double_xv_a, test_double_xv_a, _tol_double) << endl;
      cout << "2) \n"
	   << test_double_xv_a << endl
	   << " and \n"
	   << test_double_xv_b  << endl
	   << " are not identical? "
	   << std::boolalpha << (! identical(test_double_xv_a, test_double_xv_b, _tol_double)) << endl;
    }

    // template<class utype> bool
    // identical(const xvector<utype>&,const xvector<utype>&) __xprototype;
    if(identical(test_int_xv_a, test_int_xv_a) &&
       (! identical(test_int_xv_a, test_int_xv_b))) {
      cout << ">>> Testing Integer INDENTICAL with default tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Integer INDENTICAL with default tolerance:  ***FAIL***!\n"
	   << " 1) \n"
	   << test_int_xv_a << endl
	   << " and \n"
	   << test_int_xv_a  << endl
	   << " are identical? "
	   << std::boolalpha << identical(test_int_xv_a, test_int_xv_a) << endl;
      cout << "2) \n"
	   << test_int_xv_a << endl
	   << " and \n"
	   << test_int_xv_b  << endl
	   << " are not identical? "
	   << std::boolalpha << (! identical(test_int_xv_a, test_int_xv_b)) << endl;
    }

    if(identical(test_float_xv_a, test_float_xv_a) &&
       (! identical(test_float_xv_a, test_float_xv_b))) {
      cout << ">>> Testing Float INDENTICAL with default tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Float INDENTICAL with default tolerance:  ***FAIL***!\n"
	   << " 1) \n"
	   << test_float_xv_a << endl
	   << " and \n"
	   << test_float_xv_a  << endl
	   << " are identical? "
	   << std::boolalpha << identical(test_float_xv_a, test_float_xv_a) << endl;
      cout << "2) \n"
	   << test_float_xv_a << endl
	   << " and \n"
	   << test_float_xv_b  << endl
	   << " are not identical? "
	   << std::boolalpha << (! identical(test_float_xv_a, test_float_xv_b)) << endl;
    }

    if(identical(test_double_xv_a, test_double_xv_a) &&
       (! identical(test_double_xv_a, test_double_xv_b))) {
      cout << ">>> Testing Double INDENTICAL with default tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Double INDENTICAL with default tolerance:  ***FAIL***!\n"
	   << " 1) \n"
	   << test_double_xv_a << endl
	   << " and \n"
	   << test_double_xv_a  << endl
	   << " are identical? "
	   << std::boolalpha << identical(test_double_xv_a, test_double_xv_a) << endl;
      cout << "2) \n"
	   << test_double_xv_a << endl
	   << " and \n"
	   << test_double_xv_b  << endl
	   << " are not identical? "
	   << std::boolalpha << (! identical(test_double_xv_a, test_double_xv_b)) << endl;
    }

    // template<class utype> bool
    // operator==(const xvector<utype>&,const xvector<utype>&) __xprototype;
    if((test_int_xv_a == test_int_xv_a) &&
       ! (test_int_xv_a == test_int_xv_b)) {
      cout << ">>> Testing Integer == with default tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Integer == with default tolerance:  ***FAIL***!\n"
	   << " 1) \n"
	   << test_int_xv_a << endl
	   << " and \n"
	   << test_int_xv_a  << endl
	   << " are identical? "
	   << std::boolalpha << identical(test_int_xv_a, test_int_xv_a) << endl;
      cout << "2) \n"
	   << test_int_xv_a << endl
	   << " and \n"
	   << test_int_xv_b  << endl
	   << " are not identical? "
	   << std::boolalpha << (! identical(test_int_xv_a, test_int_xv_b)) << endl;
    }

    if((test_float_xv_a == test_float_xv_a) &&
       (! (test_float_xv_a == test_float_xv_b))) {
      cout << ">>> Testing Float == with default tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Float == with default tolerance:  ***FAIL***!\n"
	   << " 1) \n"
	   << test_float_xv_a << endl
	   << " and \n"
	   << test_float_xv_a  << endl
	   << " are identical? "
	   << std::boolalpha << (test_float_xv_a == test_float_xv_a) << endl;
      cout << "2) \n"
	   << test_float_xv_a << endl
	   << " and \n"
	   << test_float_xv_b  << endl
	   << " are not identical? "
	   << std::boolalpha << (! (test_float_xv_a == test_float_xv_b)) << endl;
    }

    if(identical(test_double_xv_a, test_double_xv_a) &&
       (! identical(test_double_xv_a, test_double_xv_b))) {
      cout << ">>> Testing Double == with default tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Double == with default tolerance:  ***FAIL***!\n"
	   << " 1) \n"
	   << test_double_xv_a << endl
	   << " and \n"
	   << test_double_xv_a  << endl
	   << " are identical? "
	   << std::boolalpha << (test_double_xv_a == test_double_xv_a) << endl;
      cout << "2) \n"
	   << test_double_xv_a << endl
	   << " and \n"
	   << test_double_xv_b  << endl
	   << " are not identical? "
	   << std::boolalpha << (! (test_double_xv_a == test_double_xv_b)) << endl;
    }

    // template<class utype> bool
    // isdifferent(const xvector<utype>&,const xvector<utype>&,const utype&) __xprototype;
    if(! isdifferent(test_int_xv_a, test_int_xv_a, _tol_int) &&
       isdifferent(test_int_xv_a,  test_int_xv_b, _tol_int)) {
      cout << ">>> Testing Integer ISDIFFERENT with given tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Integer ISDIFFERENT with default tolerance:  ***FAIL***!\n"
	   << " 1) \n"
	   << test_float_xv_a << endl
	   << " and \n"
	   << test_float_xv_a  << endl
	   << " are different? "
	   << std::boolalpha << isdifferent(test_int_xv_a, test_int_xv_a, _tol_int) << endl;
      cout << "2) \n"
	   << test_float_xv_a << endl
	   << " and \n"
	   << test_float_xv_b  << endl
	   << " are not different? "
	   << std::boolalpha << (! isdifferent(test_int_xv_a, test_int_xv_b, _tol_int)) << endl;
    }

    if(! isdifferent(test_float_xv_a, test_float_xv_a, _tol_float) &&
       (isdifferent(test_float_xv_a, test_float_xv_b, _tol_float))) {
      cout << ">>> Testing Float ISDIFFERENT with given tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Float ISDIFFERENT with given tolerance:  ***FAIL***!\n"
	   << " 1) \n"
	   << test_float_xv_a << endl
	   << " and \n"
	   << test_float_xv_a  << endl
	   << " are different? "
	   << std::boolalpha << isdifferent(test_float_xv_a, test_float_xv_a, _tol_float) << endl;
      cout << "2) \n"
	   << test_float_xv_a << endl
	   << " and \n"
	   << test_float_xv_b  << endl
	   << " are not different? "
	   << std::boolalpha <<(isdifferent(test_float_xv_a, test_float_xv_b, _tol_float)) << endl;
    }

    if(! isdifferent(test_double_xv_a, test_double_xv_a, _tol_double) &&
       (isdifferent(test_double_xv_a, test_double_xv_b, _tol_double))) {
      cout << ">>> Testing Double ISDIFFERENT with given tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Double ISDIFFERENT with given tolerance:  ***FAIL***!\n"
	   << " 1) \n"
	   << test_double_xv_a << endl
	   << " and \n"
	   << test_double_xv_a  << endl
	   << " are different? "
	   << std::boolalpha << isdifferent(test_double_xv_a, test_double_xv_a, _tol_double) << endl;
      cout << "2) \n"
	   << test_double_xv_a << endl
	   << " and \n"
	   << test_double_xv_b  << endl
	   << " are not different? "
	   << std::boolalpha << (! isdifferent(test_double_xv_a, test_double_xv_b, _tol_double)) << endl;
    }

    // template<class utype> bool
    // isdifferent(const xvector<utype>&,const xvector<utype>&) __xprototype;
    if(! isdifferent(test_int_xv_a, test_int_xv_a) &&
       isdifferent(test_int_xv_a,  test_int_xv_b)) {
      cout << ">>> Testing Integer ISDIFFERENT xvector with default tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Integer ISDIFFERENT xvector with default tolerance:  ***FAIL***!\n"
	   << " 1) \n"
	   << test_int_xv_a << endl
	   << " and \n"
	   << test_int_xv_a  << endl
	   << " are different? "
	   << std::boolalpha << isdifferent(test_int_xv_a, test_int_xv_a) << endl;
      cout << "2) \n"
	   << test_int_xv_a << endl
	   << " and \n"
	   << test_int_xv_b  << endl
	   << " are not different? "
	   << std::boolalpha << (! isdifferent(test_int_xv_a, test_int_xv_b)) << endl;
    }

    if(! isdifferent(test_float_xv_a, test_float_xv_a) &&
       (isdifferent(test_float_xv_a, test_float_xv_b))) {
      cout << ">>> Testing Float ISDIFFERENT with default tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Float ISDIFFERENT with default tolerance:  ***FAIL***!\n"
	   << " 1) \n"
	   << test_float_xv_a << endl
	   << " and \n"
	   << test_float_xv_a  << endl
	   << " are different? "
	   << std::boolalpha << isdifferent(test_float_xv_a, test_float_xv_a) << endl;
      cout << "2) \n"
	   << test_int_xv_a << endl
	   << " and \n"
	   << test_int_xv_b  << endl
	   << " are not different? "
	   << std::boolalpha <<(isdifferent(test_float_xv_a, test_float_xv_b)) << endl;
    }

    if(! isdifferent(test_double_xv_a, test_double_xv_a) &&
       (isdifferent(test_double_xv_a, test_double_xv_b))) {
      cout << ">>> Testing Double ISDIFFERENT with default tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Double ISDIFFERENT with default tolerance:  ***FAIL***!\n"
	   << " 1) \n"
	   << test_double_xv_a << endl
	   << " and \n"
	   << test_double_xv_a  << endl
	   << " are different? "
	   << std::boolalpha << isdifferent(test_double_xv_a, test_double_xv_a) << endl;
      cout << "2) \n"
	   << test_int_xv_a << endl
	   << " and \n"
	   << test_int_xv_b  << endl
	   << " are not different? "
	   << std::boolalpha << (! isdifferent(test_double_xv_a, test_double_xv_b)) << endl;
    }

    // template<class utype> bool
    // aurostd::isequal(const xvector<utype>&,const xvector<utype>&,const utype&) __xprototype;
    if(aurostd::isequal(test_int_xv_a, test_int_xv_a, _tol_int) &&
       (!aurostd::isequal(test_int_xv_a, test_int_xv_b, _tol_int))) {
      cout << ">>> Testing Integer ISEQUAL with given tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Integer ISEQUAL with given tolerance:  ***FAIL***!\n"
	   << " 1) \n"
	   << test_int_xv_a << endl
	   << " and \n"
	   << test_int_xv_a  << endl
	   << " are identical? "
	   << std::boolalpha << aurostd::isequal(test_int_xv_a, test_int_xv_a, _tol_int) << endl;
      cout << "2) \n"
	   << test_int_xv_a << endl
	   << " and \n"
	   << test_int_xv_b  << endl
	   << " are not identical? "
	   << std::boolalpha << (!aurostd::isequal(test_int_xv_a, test_int_xv_b, _tol_int)) << endl;
    }

    if(aurostd::isequal(test_float_xv_a, test_float_xv_a, _tol_float) &&
       (!aurostd::isequal(test_float_xv_a, test_float_xv_b, _tol_float))) {
      cout << ">>> Testing Float ISEQUAL with given tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Float ISEQUAL with given tolerance:  ***FAIL***!\n"
	   << " 1) \n"
	   << test_float_xv_a << endl
	   << " and \n"
	   << test_float_xv_a  << endl
	   << " are identical? "
	   << std::boolalpha << aurostd::isequal(test_float_xv_a, test_float_xv_a, _tol_float) << endl;
      cout << "2) \n"
	   << test_float_xv_a << endl
	   << " and \n"
	   << test_float_xv_b  << endl
	   << " are not identical? "
	   << std::boolalpha << (!aurostd::isequal(test_float_xv_a, test_float_xv_b, _tol_float)) << endl;
    }

    if(aurostd::isequal(test_double_xv_a, test_double_xv_a, _tol_double) &&
       (!aurostd::isequal(test_double_xv_a, test_double_xv_b, _tol_double))) {
      cout << ">>> Testing Double ISEQUAL with default tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Double ISEQUAL with default tolerance:  ***FAIL***!\n"
	   << " 1) \n"
	   << test_double_xv_a << endl
	   << " and \n"
	   << test_double_xv_a  << endl
	   << " are identical? "
	   << std::boolalpha << aurostd::isequal(test_double_xv_a, test_double_xv_a, _tol_double) << endl;
      cout << "2) \n"
	   << test_double_xv_a << endl
	   << " and \n"
	   << test_double_xv_b  << endl
	   << " are not identical? "
	   << std::boolalpha << (!aurostd::isequal(test_double_xv_a, test_double_xv_b, _tol_double)) << endl;
    }

    // template<class utype> bool
    // isequal(const xvector<utype>&,const xvector<utype>&) __xprototype;
    if(aurostd::isequal(test_int_xv_a, test_int_xv_a) &&
       (!aurostd::isequal(test_int_xv_a, test_int_xv_b))) {
      cout << ">>> Testing Integer ISEQUAL with default tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Integer ISEQUAL with default tolerance:  ***FAIL***!\n"
	   << " 1) \n"
	   << test_int_xv_a << endl
	   << " and \n"
	   << test_int_xv_a  << endl
	   << " are identical? "
	   << std::boolalpha << aurostd::isequal(test_int_xv_a, test_int_xv_a) << endl;
      cout << "2) \n"
	   << test_int_xv_a << endl
	   << " and \n"
	   << test_int_xv_b  << endl
	   << " are not identical? "
	   << std::boolalpha << (!aurostd::isequal(test_int_xv_a, test_int_xv_b)) << endl;
    }

    if(aurostd::isequal(test_float_xv_a, test_float_xv_a) &&
       (!aurostd::isequal(test_float_xv_a, test_float_xv_b))) {
      cout << ">>> Testing Float ISEQUAL with default tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Float ISEQUAL with default tolerance:  ***FAIL***!\n"
	   << " 1) \n"
	   << test_float_xv_a << endl
	   << " and \n"
	   << test_float_xv_a  << endl
	   << " are identical? "
	   << std::boolalpha << aurostd::isequal(test_float_xv_a, test_float_xv_a) << endl;
      cout << "2) \n"
	   << test_float_xv_a << endl
	   << " and \n"
	   << test_float_xv_b  << endl
	   << " are not identical? "
	   << std::boolalpha << (!aurostd::isequal(test_float_xv_a, test_float_xv_b)) << endl;
    }

    if(aurostd::isequal(test_double_xv_a, test_double_xv_a) &&
       (!aurostd::isequal(test_double_xv_a, test_double_xv_b))) {
      cout << ">>> Testing Double ISEQUAL with default tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Double ISEQUAL with default tolerance:  ***FAIL***!\n"
	   << " 1) \n"
	   << test_double_xv_a << endl
	   << " and \n"
	   << test_double_xv_a  << endl
	   << " are identical? "
	   << std::boolalpha << aurostd::isequal(test_double_xv_a, test_double_xv_a) << endl;
      cout << "2) \n"
	   << test_double_xv_a << endl
	   << " and \n"
	   << test_double_xv_b  << endl
	   << " are not identical? "
	   << std::boolalpha << (!aurostd::isequal(test_double_xv_a, test_double_xv_b)) << endl;
    }

    // template<class utype> bool
    // operator!=(const xvector<utype>&,const xvector<utype>&) __xprototype;
    if(! (test_int_xv_a != test_int_xv_a) &&
       ((test_int_xv_a != test_int_xv_b))) {
      cout << ">>> Testing Integer != with default tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Integer != with default tolerance:  ***FAIL***!\n"
	   << " 1) \n"
	   << test_int_xv_a << endl
	   << " and \n"
	   << test_int_xv_a  << endl
	   << " are different? "
	   << std::boolalpha << (test_int_xv_a != test_int_xv_a) << endl;
      cout << "2) \n"
	   << test_int_xv_a << endl
	   << " and \n"
	   << test_int_xv_b  << endl
	   << " are not different? "
	   << std::boolalpha <<((test_int_xv_a != test_int_xv_b)) << endl;
    }

    if(! (test_float_xv_a != test_float_xv_a) &&
       ((test_float_xv_a != test_float_xv_b))) {
      cout << ">>> Testing Float != with default tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Float != with default tolerance:  ***FAIL***!\n"
	   << " 1) \n"
	   << test_int_xv_a << endl
	   << " and \n"
	   << test_int_xv_a  << endl
	   << " are different? "
	   << std::boolalpha << (test_float_xv_a != test_float_xv_a) << endl;
      cout << "2) \n"
	   << test_int_xv_a << endl
	   << " and \n"
	   << test_int_xv_b  << endl
	   << " are not different? "
	   << std::boolalpha << ((test_float_xv_a != test_float_xv_b)) << endl;
    }

    if(! (test_double_xv_a != test_double_xv_a) &&
       ((test_double_xv_a != test_double_xv_b))) {
      cout << ">>> Testing Double != with default tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Double != with default tolerance:  ***FAIL***!\n"
	   << " 1) \n"
	   << test_int_xv_a << endl
	   << " and \n"
	   << test_int_xv_a  << endl
	   << " are different? "
	   << std::boolalpha << (test_double_xv_a != test_double_xv_a) << endl;
      cout << "2) \n"
	   << test_int_xv_a << endl
	   << " and \n"
	   << test_int_xv_b  << endl
	   << " are not different? "
	   << std::boolalpha << ((test_double_xv_a != test_double_xv_b)) << endl;
    }

    const int test_xm_size = 3;
    xmatrix<int> test_int_xm_a(test_xm_size, test_xm_size), test_int_xm_b(test_xm_size, test_xm_size);
    xmatrix<float> test_float_xm_a(test_xm_size, test_xm_size), test_float_xm_b(test_xm_size, test_xm_size);
    xmatrix<double> test_double_xm_a(test_xm_size, test_xm_size), test_double_xm_b(test_xm_size, test_xm_size);

    // initialze random seed
    srand(time(NULL));

    // initializae test xmatrices
    for(int i=1; i<test_xm_size+1; i++) {
      for(int j=1; j<test_xm_size+1; j++) {
	test_int_xm_a[i][j] =  rand() % _max_int;
	test_int_xm_b[i][j] =  rand() % _max_int;
	test_float_xm_a[i][j] = PI * float((rand()%_max_int));
	test_float_xm_b[i][j] = PI * float((rand()%_max_int));
	test_double_xm_a[i][j] = _mm_e * double((rand()%_max_int));
	test_double_xm_b[i][j] = _mm_e * double((rand()%_max_int));
      }
    }

    for(int i=0; i<line_width; i++) {
      cout << "*";
    }
    cout << endl;

    cout << "xmatrix\n";

    // template<class utype> bool
    // identical(const xmatrix<utype>&,const xmatrix<utype>&,const utype&) __xprototype;
    //cout << "Indentical Integer with given tolerance " << _tol_int << " : \n"

    if(identical(test_int_xm_a, test_int_xm_a, _tol_int) &&
       (! identical(test_int_xm_a, test_int_xm_b, _tol_int))) {
      cout << ">>> Testing Integer INDENTICAL with given tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Integer INDENTICAL with given tolerance:  ***FAIL***!\n"
	   << " 1) \n"
	   << test_int_xm_a << endl
	   << " and \n"
	   << test_int_xm_a  << endl
	   << " are identical? "
	   << std::boolalpha << identical(test_int_xm_a, test_int_xm_a, _tol_int) << endl;
      cout << "2) \n"
	   << test_int_xm_a << endl
	   << " and \n"
	   << test_int_xm_b  << endl
	   << " are not identical? "
	   << std::boolalpha << (! identical(test_int_xm_a, test_int_xm_b, _tol_int)) << endl;
    }

    if(identical(test_float_xm_a, test_float_xm_a, _tol_float) &&
       (! identical(test_float_xm_a, test_float_xm_b, _tol_float))) {
      cout << ">>> Testing Float INDENTICAL with given tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Float INDENTICAL with given tolerance:  ***FAIL***!\n"
	   << " 1) \n"
	   << test_float_xm_a << endl
	   << " and \n"
	   << test_float_xm_a  << endl
	   << " are identical? "
	   << std::boolalpha << identical(test_float_xm_a, test_float_xm_a, _tol_float) << endl;
      cout << "2) \n"
	   << test_float_xm_a << endl
	   << " and \n"
	   << test_float_xm_b  << endl
	   << " are not identical? "
	   << std::boolalpha << (! identical(test_float_xm_a, test_float_xm_b, _tol_float)) << endl;
    }

    if(identical(test_double_xm_a, test_double_xm_a, _tol_double) &&
       (! identical(test_double_xm_a, test_double_xm_b, _tol_double))) {
      cout << ">>> Testing Double INDENTICAL with default tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Double INDENTICAL with default tolerance:  ***FAIL***!\n"
	   << " 1) \n"
	   << test_double_xm_a << endl
	   << " and \n"
	   << test_double_xm_a  << endl
	   << " are identical? "
	   << std::boolalpha << identical(test_double_xm_a, test_double_xm_a, _tol_double) << endl;
      cout << "2) \n"
	   << test_double_xm_a << endl
	   << " and \n"
	   << test_double_xm_b  << endl
	   << " are not identical? "
	   << std::boolalpha << (! identical(test_double_xm_a, test_double_xm_b, _tol_double)) << endl;
    }

    // template<class utype> bool
    // identical(const xmatrix<utype>&,const xmatrix<utype>&) __xprototype;
    if(identical(test_int_xm_a, test_int_xm_a) &&
       (! identical(test_int_xm_a, test_int_xm_b))) {
      cout << ">>> Testing Integer INDENTICAL with default tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Integer INDENTICAL with default tolerance:  ***FAIL***!\n"
	   << " 1) \n"
	   << test_int_xm_a << endl
	   << " and \n"
	   << test_int_xm_a  << endl
	   << " are identical? "
	   << std::boolalpha << identical(test_int_xm_a, test_int_xm_a) << endl;
      cout << "2) \n"
	   << test_int_xm_a << endl
	   << " and \n"
	   << test_int_xm_b  << endl
	   << " are not identical? "
	   << std::boolalpha <<(! identical(test_int_xm_a, test_int_xm_b)) << endl;
    }

    bool pass, pass1;
    pass= identical(test_float_xm_a, test_float_xm_a); if(0) cout << pass;
    pass1= identical(test_float_xm_a, test_float_xm_b); if(0) cout << pass1;

    if(identical(test_float_xm_a, test_float_xm_a) &&
       (! identical(test_float_xm_a, test_float_xm_b))) {
      cout << ">>> Testing Float INDENTICAL with default tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Float INDENTICAL with default tolerance:  ***FAIL***!\n"
	   << " 1) \n"
	   << test_float_xm_a << endl
	   << " and \n"
	   << test_float_xm_a  << endl
	   << " are identical? "
	   << std::boolalpha << identical(test_float_xm_a, test_float_xm_a) << endl;
      cout << "2) \n"
	   << test_float_xm_a << endl
	   << " and \n"
	   << test_float_xm_b  << endl
	   << " are not identical? "
	   << std::boolalpha << (! identical(test_float_xm_a, test_float_xm_b)) << endl;
    }

    if(identical(test_double_xm_a, test_double_xm_a) &&
       (! identical(test_double_xm_a, test_double_xm_b))) {
      cout << ">>> Testing Double INDENTICAL with default tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Double INDENTICAL with default tolerance:  ***FAIL***!\n"
	   << " 1) \n"
	   << test_double_xm_a << endl
	   << " and \n"
	   << test_double_xm_a  << endl
	   << " are identical? "
	   << std::boolalpha << identical(test_double_xm_a, test_double_xm_a) << endl;
      cout << "2) \n"
	   << test_double_xm_a << endl
	   << " and \n"
	   << test_double_xm_b  << endl
	   << " are not identical? "
	   << std::boolalpha << (! identical(test_double_xm_a, test_double_xm_b)) << endl;
    }

    // template<class utype> bool
    // operator==(const xvector<utype>&,const xvector<utype>&) __xprototype;
    if((test_int_xm_a == test_int_xm_a) &&
       ! (test_int_xm_a == test_int_xm_b)) {
      cout << ">>> Testing Integer == with default tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Integer == with default tolerance:  ***FAIL***!\n"
	   << " 1) \n"
	   << test_int_xm_a << endl
	   << " and \n"
	   << test_int_xm_a  << endl
	   << " are identical? "
	   << std::boolalpha << identical(test_int_xm_a, test_int_xm_a) << endl;
      cout << "2) \n"
	   << test_int_xm_a << endl
	   << " and \n"
	   << test_int_xm_b  << endl
	   << " are not identical? "
	   << std::boolalpha << (! identical(test_int_xm_a, test_int_xm_b)) << endl;
    }

    if((test_float_xm_a == test_float_xm_a) &&
       (! (test_float_xm_a == test_float_xm_b))) {
      cout << ">>> Testing Float == with default tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Float == with default tolerance:  ***FAIL***!\n"
	   << " 1) \n"
	   << test_float_xm_a << endl
	   << " and \n"
	   << test_float_xm_a  << endl
	   << " are identical? "
	   << std::boolalpha << (test_float_xm_a == test_float_xm_a) << endl;
      cout << "2) \n"
	   << test_float_xm_a << endl
	   << " and \n"
	   << test_float_xm_b  << endl
	   << " are not identical? "
	   << std::boolalpha << (! (test_float_xm_a == test_float_xm_b)) << endl;
    }

    if(identical(test_double_xm_a, test_double_xm_a) &&
       (! identical(test_double_xm_a, test_double_xm_b))) {
      cout << ">>> Testing Double == with default tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Double == with default tolerance:  ***FAIL***!\n"
	   << " 1) \n"
	   << test_double_xm_a << endl
	   << " and \n"
	   << test_double_xm_a  << endl
	   << " are identical? "
	   << std::boolalpha << (test_double_xm_a == test_double_xm_a) << endl;
      cout << "2) \n"
	   << test_double_xm_a << endl
	   << " and \n"
	   << test_double_xm_b  << endl
	   << " are not identical? "
	   << std::boolalpha << (! (test_double_xm_a == test_double_xm_b)) << endl;
    }

    // template<class utype> bool
    // isdifferent(const xmatrix<utype>&,const xmatrix<utype>&,const utype&) __xprototype;
    if(! isdifferent(test_int_xm_a, test_int_xm_a, _tol_int) &&
       isdifferent(test_int_xm_a,  test_int_xm_b, _tol_int)) {
      cout << ">>> Testing Integer ISDIFFERENT with given tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Integer ISDIFFERENT with default tolerance:  ***FAIL***!\n"
	   << " 1) \n"
	   << test_int_xm_a << endl
	   << " and \n"
	   << test_int_xm_a  << endl
	   << " are different? "
	   << std::boolalpha << isdifferent(test_int_xm_a, test_int_xm_a, _tol_int) << endl;
      cout << "2) \n"
	   << test_int_xm_a << endl
	   << " and \n"
	   << test_int_xm_b  << endl
	   << " are not different? "
	   << std::boolalpha << (! isdifferent(test_int_xm_a, test_int_xm_b, _tol_int)) << endl;
    }

    if(! isdifferent(test_float_xm_a, test_float_xm_a, _tol_float) &&
       (isdifferent(test_float_xm_a, test_float_xm_b, _tol_float))) {
      cout << ">>> Testing Float ISDIFFERENT with given tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Float ISDIFFERENT with given tolerance:  ***FAIL***!\n"
	   << " 1) \n"
	   << test_float_xm_a << endl
	   << " and \n"
	   << test_float_xm_a  << endl
	   << " are different? "
	   << std::boolalpha << isdifferent(test_float_xm_a, test_float_xm_a, _tol_float) << endl;
      cout << "2) \n"
	   << test_float_xm_a << endl
	   << " and \n"
	   << test_float_xm_b  << endl
	   << " are not different? "
	   << std::boolalpha <<(isdifferent(test_float_xm_a, test_float_xm_b, _tol_float)) << endl;
    }

    if(! isdifferent(test_double_xm_a, test_double_xm_a, _tol_double) &&
       (isdifferent(test_double_xm_a, test_double_xm_b, _tol_double))) {
      cout << ">>> Testing Double ISDIFFERENT with given tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Double ISDIFFERENT with given tolerance:  ***FAIL***!\n"
	   << " 1) \n"
	   << test_double_xm_a << endl
	   << " and \n"
	   << test_double_xm_a  << endl
	   << " are different? "
	   << std::boolalpha << isdifferent(test_double_xm_a, test_double_xm_a, _tol_double) << endl;
      cout << "2) \n"
	   << test_double_xm_a << endl
	   << " and \n"
	   << test_double_xm_b  << endl
	   << " are not different? "
	   << std::boolalpha << (! isdifferent(test_double_xm_a, test_double_xm_b, _tol_double)) << endl;
    }

    // template<class utype> bool
    // isdifferent(const xmatrix<utype>&,const xmatrix<utype>&) __xprototype;
    if(! isdifferent(test_int_xm_a, test_int_xm_a) &&
       isdifferent(test_int_xm_a,  test_int_xm_b)) {
      cout << ">>> Testing Integer ISDIFFERENT xvector with default tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Integer ISDIFFERENT xvector with default tolerance:  ***FAIL***!\n"
	   << " 1) \n"
	   << test_int_xm_a << endl
	   << " and \n"
	   << test_int_xm_a  << endl
	   << " are different? "
	   << std::boolalpha << isdifferent(test_int_xm_a, test_int_xm_a) << endl;
      cout << "2) \n"
	   << test_int_xm_a << endl
	   << " and \n"
	   << test_int_xm_b  << endl
	   << " are not different? "
	   << std::boolalpha << (! isdifferent(test_int_xm_a, test_int_xm_b)) << endl;
    }

    if(! isdifferent(test_float_xm_a, test_float_xm_a) &&
       (isdifferent(test_float_xm_a, test_float_xm_b))) {
      cout << ">>> Testing Float ISDIFFERENT with default tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Float ISDIFFERENT with default tolerance:  ***FAIL***!\n"
	   << " 1) \n"
	   << test_float_xm_a << endl
	   << " and \n"
	   << test_float_xm_a  << endl
	   << " are different? "
	   << std::boolalpha << isdifferent(test_float_xm_a, test_float_xm_a) << endl;
      cout << "2) \n"
	   << test_float_xm_a << endl
	   << " and \n"
	   << test_float_xm_b  << endl
	   << " are not different? "
	   << std::boolalpha <<(isdifferent(test_float_xm_a, test_float_xm_b)) << endl;
    }

    if(! isdifferent(test_double_xm_a, test_double_xm_a) &&
       (isdifferent(test_double_xm_a, test_double_xm_b))) {
      cout << ">>> Testing Double ISDIFFERENT with default tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Double ISDIFFERENT with default tolerance:  ***FAIL***!\n"
	   << " 1) \n"
	   << test_double_xm_a << endl
	   << " and \n"
	   << test_double_xm_a  << endl
	   << " are different? "
	   << std::boolalpha << isdifferent(test_double_xm_a, test_double_xm_a) << endl;
      cout << "2) \n"
	   << test_double_xm_a << endl
	   << " and \n"
	   << test_double_xm_b  << endl
	   << " are not different? "
	   << std::boolalpha << (! isdifferent(test_double_xm_a, test_double_xm_b)) << endl;
    }

    // template<class utype> bool
    // isequal(const xmatrix<utype>&,const xmatrix<utype>&,const utype&) __xprototype;
    if(aurostd::isequal(test_int_xm_a, test_int_xm_a, _tol_int) &&
       (!aurostd::isequal(test_int_xm_a, test_int_xm_b, _tol_int))) {
      cout << ">>> Testing Integer ISEQUAL with given tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Integer ISEQUAL with given tolerance:  ***FAIL***!\n"
	   << " 1) \n"
	   << test_int_xm_a << endl
	   << " and \n"
	   << test_int_xm_a  << endl
	   << " are identical? "
	   << std::boolalpha << aurostd::isequal(test_int_xm_a, test_int_xm_a, _tol_int) << endl;
      cout << "2) \n"
	   << test_int_xm_a << endl
	   << " and \n"
	   << test_int_xm_b  << endl
	   << " are not identical? "
	   << std::boolalpha << (!aurostd::isequal(test_int_xm_a, test_int_xm_b, _tol_int)) << endl;
    }

    if(aurostd::isequal(test_float_xm_a, test_float_xm_a, _tol_float) &&
       (!aurostd::isequal(test_float_xm_a, test_float_xm_b, _tol_float))) {
      cout << ">>> Testing Float ISEQUAL with given tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Float ISEQUAL with given tolerance:  ***FAIL***!\n"
	   << " 1) \n"
	   << test_float_xm_a << endl
	   << " and \n"
	   << test_float_xm_a  << endl
	   << " are identical? "
	   << std::boolalpha << aurostd::isequal(test_float_xm_a, test_float_xm_a, _tol_float) << endl;
      cout << "2) \n"
	   << test_float_xm_a << endl
	   << " and \n"
	   << test_float_xm_b  << endl
	   << " are not identical? "
	   << std::boolalpha << (!aurostd::isequal(test_float_xm_a, test_float_xm_b, _tol_float)) << endl;
    }

    if(aurostd::isequal(test_double_xm_a, test_double_xm_a, _tol_double) &&
       (!aurostd::isequal(test_double_xm_a, test_double_xm_b, _tol_double))) {
      cout << ">>> Testing Double ISEQUAL with default tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Double ISEQUAL with default tolerance:  ***FAIL***!\n"
	   << " 1) \n"
	   << test_double_xm_a << endl
	   << " and \n"
	   << test_double_xm_a  << endl
	   << " are identical? "
	   << std::boolalpha << aurostd::isequal(test_double_xm_a, test_double_xm_a, _tol_double) << endl;
      cout << "2) \n"
	   << test_double_xm_a << endl
	   << " and \n"
	   << test_double_xm_b  << endl
	   << " are not identical? "
	   << std::boolalpha << (!aurostd::isequal(test_double_xm_a, test_double_xm_b, _tol_double)) << endl;
    }

    // template<class utype> bool
    // isequal(const xmatrix<utype>&,const xmatrix<utype>&) __xprototype;
    if(aurostd::isequal(test_int_xm_a, test_int_xm_a) &&
       (!aurostd::isequal(test_int_xm_a, test_int_xm_b))) {
      cout << ">>> Testing Integer ISEQUAL with default tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Integer ISEQUAL with default tolerance:  ***FAIL***!\n"
	   << " 1) \n"
	   << test_int_xm_a << endl
	   << " and \n"
	   << test_int_xm_a  << endl
	   << " are identical? "
	   << std::boolalpha << aurostd::isequal(test_int_xm_a, test_int_xm_a) << endl;
      cout << "2) \n"
	   << test_int_xm_a << endl
	   << " and \n"
	   << test_int_xm_b  << endl
	   << " are not identical? "
	   << std::boolalpha << (!aurostd::isequal(test_int_xm_a, test_int_xm_b)) << endl;
    }

    if(aurostd::isequal(test_float_xm_a, test_float_xm_a) &&
       (!aurostd::isequal(test_float_xm_a, test_float_xm_b))) {
      cout << ">>> Testing Float ISEQUAL with default tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Float ISEQUAL with default tolerance:  ***FAIL***!\n"
	   << " 1) \n"
	   << test_float_xm_a << endl
	   << " and \n"
	   << test_float_xm_a  << endl
	   << " are identical? "
	   << std::boolalpha << aurostd::isequal(test_float_xm_a, test_float_xm_a) << endl;
      cout << "2) \n"
	   << test_float_xm_a << endl
	   << " and \n"
	   << test_float_xm_b  << endl
	   << " are not identical? "
	   << std::boolalpha << (!aurostd::isequal(test_float_xm_a, test_float_xm_b)) << endl;
    }

    if(aurostd::isequal(test_double_xm_a, test_double_xm_a) &&
       (!aurostd::isequal(test_double_xm_a, test_double_xm_b))) {
      cout << ">>> Testing Double ISEQUAL with default tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Double ISEQUAL with default tolerance:  ***FAIL***!\n"
	   << " 1) \n"
	   << test_double_xm_a << endl
	   << " and \n"
	   << test_double_xm_a  << endl
	   << " are identical? "
	   << std::boolalpha << aurostd::isequal(test_double_xm_a, test_double_xm_a) << endl;
      cout << "2) \n"
	   << test_double_xm_a << endl
	   << " and \n"
	   << test_double_xm_b  << endl
	   << " are not identical? "
	   << std::boolalpha << (!aurostd::isequal(test_double_xm_a, test_double_xm_b)) << endl;
    }

    // template<class utype> bool
    // operator!=(const xmatrix<utype>&,const xmatrix<utype>&) __xprototype;
    if(! (test_int_xm_a != test_int_xm_a) &&
       ((test_int_xm_a != test_int_xm_b))) {
      cout << ">>> Testing Integer != with default tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Integer != with default tolerance:  ***FAIL***!\n"
	   << " 1) \n"
	   << test_int_xm_a << endl
	   << " and \n"
	   << test_int_xm_a  << endl
	   << " are different? "
	   << std::boolalpha << (test_int_xm_a != test_int_xm_a) << endl;
      cout << "2) \n"
	   << test_int_xm_a << endl
	   << " and \n"
	   << test_int_xm_b  << endl
	   << " are not different? "
	   << std::boolalpha <<((test_int_xm_a != test_int_xm_b)) << endl;
    }

    if(! (test_float_xm_a != test_float_xm_a) &&
       ((test_float_xm_a != test_float_xm_b))) {
      cout << ">>> Testing Float != with default tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Float != with default tolerance:  ***FAIL***!\n"
	   << " 1) \n"
	   << test_float_xm_a << endl
	   << " and \n"
	   << test_float_xm_a  << endl
	   << " are different? "
	   << std::boolalpha << (test_float_xm_a != test_float_xm_a) << endl;
      cout << "2) \n"
	   << test_float_xm_a << endl
	   << " and \n"
	   << test_float_xm_b  << endl
	   << " are not different? "
	   << std::boolalpha << ((test_float_xm_a != test_float_xm_b)) << endl;
    }

    if(! (test_double_xm_a != test_double_xm_a) &&
       ((test_double_xm_a != test_double_xm_b))) {
      cout << ">>> Testing Double != with default tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Double != with default tolerance:  ***FAIL***!\n"
	   << " 1) \n"
	   << test_double_xm_a << endl
	   << " and \n"
	   << test_double_xm_a  << endl
	   << " are different? "
	   << std::boolalpha << (test_double_xm_a != test_double_xm_a) << endl;
      cout << "2) \n"
	   << test_double_xm_a << endl
	   << " and \n"
	   << test_double_xm_b  << endl
	   << " are not different? "
	   << std::boolalpha << ((test_double_xm_a != test_double_xm_b)) << endl;
    }
  }
} // namespace pflow
// ***************************************************************************
// pflow::CHGDIFF
// ***************************************************************************
namespace pflow {
  string CHGDIFF(aurostd::xoption vpflow) {
    // handles flags for CHGDIFF
    string soliloquy="pflow::CHGDIFF():  ";     // so you know who's talking
    string chgcar1_file,chgcar2_file,output_file;
    ostringstream oss;

    // DEBUG
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) oss << soliloquy << "BEGIN FLAGS" << endl; 

    string usage_usage="aflow --chgdiff=CHGCAR1,CHGCAR2 [chgdiff_options]";
    string usage_options=aurostd::liststring2string("options = --usage",
                                                    "          --output=|-o=CHGCAR_OUT");

    // output usage
    if(LDEBUG) oss << soliloquy << "CHECK USAGE" << endl; 
    if(vpflow.flag("CHGDIFF::USAGE")) {
      init::ErrorOption(cout,vpflow.getattachedscheme("CHGDIFF"),"pflow::CHGDIFF()",aurostd::liststring2string(usage_usage,usage_options));
      return oss.str();
    }

    // no input
    vector<string> input_tokens;
    if(LDEBUG) oss << soliloquy << "CHECK INPUT" << endl; 
    if(!vpflow.flag("CHGDIFF")) {
      oss << endl;
      oss << soliloquy << "ERROR: No input given." << endl;
      oss << soliloquy << "Exiting." << endl;
      oss << endl;
      init::ErrorOption(cout,vpflow.getattachedscheme("CHGDIFF"),"pflow::CHGDIFF",aurostd::liststring2string(usage_usage,usage_options));
      return oss.str();
    }

    // check inputs
    string misc_option;
    if(LDEBUG) oss << soliloquy << "vpflow.getattachedscheme(\"CHGDIFF\")=" << vpflow.getattachedscheme("CHGDIFF") << endl;
    misc_option=vpflow.getattachedscheme("CHGDIFF");
    aurostd::string2tokens(misc_option,input_tokens,",");
    if(input_tokens.size()!=2) {
      oss << endl;
      oss << soliloquy << "ERROR: Incorrect input arguments. List two CHGCARs separated by commas." << endl;
      oss << soliloquy << "Exiting." << endl;
      oss << endl;
      init::ErrorOption(cout,vpflow.getattachedscheme("CHGDIFF"),"pflow::CHGDIFF",aurostd::liststring2string(usage_usage,usage_options));
      return oss.str();
    }

    chgcar1_file=input_tokens.at(0);chgcar2_file=input_tokens.at(1);

    // get output, if specified (or standardize)
    if(!vpflow.flag("CHGDIFF::OUTPUT")) {
      output_file="aflow_CHGDIFF.out";
    } else {
      output_file=vpflow.getattachedscheme("CHGDIFF::OUTPUT");
    }
    if(LDEBUG) oss << soliloquy << "CHGCAR_OUT=" << output_file << endl;

    CHGDIFF(chgcar1_file,chgcar2_file,output_file,oss);
    return oss.str();
  }
} // namespace pflow

// ***************************************************************************
// pflow::CHGDIFF
// ***************************************************************************
namespace pflow {
  bool CHGDIFF(const string& chgcar1_file,const string& chgcar2_file,const string& output_file, ostream& oss) {
    // RETURNS CHGCAR_OUT=CHGCAR_INPUT_1-CHGCAR_INPUT_2 
    // Read in CHGCAR or AECCAR files
    // format_dim is as follows: numcolumns chg_tot, npts and numcolumns for 
    // augmentation occupancies, numcolumns chg_diff, npts and numcolumns for
    // augmentation occupancies
    // read chgcars

    stringstream chgcar1_ss,chgcar2_ss,chgcar1_header,chgcar2_header;
    string soliloquy="pflow::CHGDIFF():  ";     // so you know who's talking
    double TOL=1e-5;
    xstructure structure1,structure2;
    vector<int> ngrid1(3),ngrid2(3),format_dim1,format_dim2;
    vector<double> chg_tot1,chg_tot2;
    vector<double> chg_diff1,chg_diff2;

    // DEBUG
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) oss << soliloquy << "BEGIN FUNCTION" << endl; 

    // get input 1
    if(!aurostd::FileExist(chgcar1_file)) {
      oss << endl;
      oss << soliloquy << "ERROR: " << chgcar1_file << " does not exist." << endl;
      oss << soliloquy << "Exiting." << endl;
      oss << endl;
      return FALSE;
    }
    if(LDEBUG) oss << soliloquy << "CHGCAR1=" << chgcar1_file << endl;
    aurostd::efile2stringstream(chgcar1_file,chgcar1_ss);

    // get input 2
    if(!aurostd::FileExist(chgcar2_file)) {
      oss << endl;
      oss << soliloquy << "ERROR: " << chgcar2_file << " does not exist." << endl;
      oss << soliloquy << "Exiting." << endl;
      oss << endl;
      return FALSE;
    } 
    if(LDEBUG) oss << soliloquy << "CHGCAR2=" << chgcar2_file << endl;
    aurostd::efile2stringstream(chgcar2_file,chgcar2_ss);

    if(LDEBUG) oss << soliloquy << "CHECK FORMAT OF " << chgcar1_file << endl;
    if(!pflow::ReadCHGCAR(structure1,chgcar1_header,ngrid1,format_dim1,chg_tot1,chg_diff1,chgcar1_ss,oss)) {
      oss << endl;
      oss << soliloquy << "ERROR: Input " << chgcar1_file << " format not recognized." << endl;
      oss << endl;
      return FALSE;       // CHGCAR format not recognized
    }
    if(LDEBUG) oss << soliloquy << "SUCCESSFULLY READ " << chgcar1_file << endl;
    if(LDEBUG) oss << soliloquy << "CHECK FORMAT OF " << chgcar2_file << endl;
    if(!pflow::ReadCHGCAR(structure2,chgcar2_header,ngrid2,format_dim2,chg_tot2,chg_diff2,chgcar2_ss,oss)) {
      oss << endl;
      oss << soliloquy << "ERROR: Input " << chgcar2_file << " format not recognized." << endl;
      oss << endl;
      return FALSE;       // CHGCAR format not recognized
    }
    if(LDEBUG) oss << soliloquy << "SUCCESSFULLY READ " << chgcar2_file << endl;

    // check formats
    pflow::matrix<double> lat1=pflow::GetLat(structure1);
    pflow::matrix<double> lat2=pflow::GetLat(structure2);
    if(LDEBUG) oss << soliloquy << "CHECK IF FORMATS OF CHGCARS MATCH" << endl;
    if(format_dim1!=format_dim2) {
      oss << endl;
      oss << soliloquy << "WARNING: Format for input " << chgcar1_file << " does not match that of input " << chgcar2_file << "." << endl;
      oss << soliloquy << "WARNING: Using format from input " << chgcar1_file << "." << endl;
      oss << endl;
    }
    if(LDEBUG) oss << soliloquy << "CHECK IF GRIDS OF CHGCARS MATCH" << endl;
    if(ngrid1!=ngrid2) {
      oss << endl;
      oss << soliloquy << "ERROR: Grids for two CHGCAR's do not match. " << endl;
      oss << soliloquy << "ERROR: ngrid of " << chgcar1_file << ": " << ngrid1.at(0) << " " <<  ngrid1.at(1) << " " <<  ngrid1.at(2) << endl;
      oss << soliloquy << "ERROR: ngrid of " << chgcar2_file << ": " << ngrid2.at(0) << " " <<  ngrid2.at(1) << " " <<  ngrid2.at(2) << endl;
      oss << soliloquy << "ERROR: This will give nonsense. " << endl;
      oss << endl;
      return FALSE;
    }
    if(LDEBUG) oss << soliloquy << "CHECK IF LATTICE PARAMETERS OF CHGCARS MATCH" << endl;
    if(pflow::norm(pflow::VVdiff(lat1[0],lat2[0]))>TOL
       || pflow::norm(pflow::VVdiff(lat1[1],lat2[1]))>TOL
       || pflow::norm(pflow::VVdiff(lat1[2],lat2[2]))>TOL) {
      oss << endl;
      oss << soliloquy << "WARNING: Lattice parameters for two CHGCARs do not match. " << endl;
      oss << soliloquy << "WARNING: lattice of " << chgcar1_file << ": " << endl;
      oss << soliloquy << matrix2xmatrix(lat1) << endl;
      oss << soliloquy << "WARNING: lattice of " << chgcar2_file << ": " << endl;
      oss << matrix2xmatrix(lat2) << endl;
      oss << soliloquy << "WARNING: This could give nonsense if there is much difference. " << endl;
      oss << soliloquy << "WARNING: Output will use lattice of " << chgcar1_file << "." << endl;
      oss << endl;
    }
    if(LDEBUG) oss << soliloquy << "OVERALL FORMAT OF CHGCARS LOOKS OK" << endl;

    // Get difference
    vector<double> chg_tot_2m1=pflow::VVdiff(chg_tot1,chg_tot2);
    vector<double> chg_diff_2m1=pflow::VVdiff(chg_diff1,chg_diff2);
    if(LDEBUG) oss << soliloquy << "PRINTING CHGDIFF TO " << output_file << endl;
    pflow::PrintCHGCAR(structure1,chgcar1_header,ngrid1,format_dim1,chg_tot_2m1,chg_diff_2m1,output_file,oss);
    if(LDEBUG) oss << soliloquy << "DONE" << endl;
    oss << soliloquy << output_file << " generated." << endl;
    return TRUE;
  }
} // namespace pflow

//DX and CO - START
// ***************************************************************************
// pflow::CHGINT
// ***************************************************************************
namespace pflow {
  void CHGINT(vector<string> argv) {
    ifstream chgfile(argv.at(2).c_str());
    xstructure str;
    vector<int> ngrid(3);
    // Read in charge
    vector<double> chg_tot;
    vector<double> chg_diff;
    pflow::ReadChg(str,ngrid,chg_tot,chg_diff,chgfile);
    // Integrate charge
    vector<pflow::matrix<double> > rad_chg_int;
    pflow::matrix<double> vor_chg_int;
    pflow::GetChgInt(rad_chg_int,vor_chg_int,str,ngrid,chg_tot,chg_diff);
    // Print results
    pflow::PrintChgInt(rad_chg_int,vor_chg_int,cout);
  }
} // namespace pflow
//DX and CO - END

// ***************************************************************************
// pflow::CHGSUM
// ***************************************************************************
namespace pflow {
  string CHGSUM(aurostd::xoption vpflow) {
    // handles flags for CHGSUM

    string soliloquy="pflow::CHGSUM():  ";     // so you know who's talking
    vector<string> chgcar_files;
    ostringstream oss;
    string output_file;

    // DEBUG
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) oss << soliloquy << "BEGIN FLAGS" << endl; 

    string usage_usage="aflow --chgsum=CHGCAR1,CHGCAR2,... [chgsum_options]";
    string usage_options=aurostd::liststring2string("options = --usage",
                                                    "          --output=|-o=CHGCAR_OUT");

    // output usage
    if(LDEBUG) oss << soliloquy << "CHECK USAGE" << endl; 
    if(vpflow.flag("CHGSUM::USAGE")) {
      init::ErrorOption(cout,vpflow.getattachedscheme("CHGSUM"),"pflow::CHGSUM",aurostd::liststring2string(usage_usage,usage_options));
      return oss.str();
    }

    // no input
    vector<string> input_tokens;
    if(LDEBUG) oss << soliloquy << "CHECK INPUT" << endl; 
    if(!vpflow.flag("CHGSUM")) {
      oss << endl;
      oss << soliloquy << "ERROR: No input given." << endl;
      oss << soliloquy << "Exiting." << endl;
      oss << endl;
      init::ErrorOption(cout,vpflow.getattachedscheme("CHGSUM"),"pflow::CHGSUM",aurostd::liststring2string(usage_usage,usage_options));
      oss << endl;
      return oss.str();
    }

    // check inputs
    string misc_option;
    if(LDEBUG) oss << soliloquy << "vpflow.getattachedscheme(\"CHGSUM\")=" << vpflow.getattachedscheme("CHGSUM") << endl;
    misc_option=vpflow.getattachedscheme("CHGSUM");
    aurostd::string2tokens(misc_option,input_tokens,",");
    if(input_tokens.size()<2) {
      oss << endl;
      oss << soliloquy << "ERROR: Incorrect input arguments. List at least two CHGCARs separated by commas." << endl;
      oss << soliloquy << "Exiting." << endl;
      oss << endl;
      init::ErrorOption(cout,vpflow.getattachedscheme("CHGSUM"),"pflow::CHGSUM",aurostd::liststring2string(usage_usage,usage_options));
      return oss.str();
    } else {
      // get inputs
      for(uint i=0;i<input_tokens.size();i++) {
        chgcar_files.push_back(input_tokens.at(i));
      }
    }

    // get output, if specified (or standardize)
    if(!vpflow.flag("CHGSUM::OUTPUT")) {
      output_file="aflow_CHGSUM.out";
    } else {
      output_file=vpflow.getattachedscheme("CHGSUM::OUTPUT");
    }
    if(LDEBUG) oss << soliloquy << "CHGCAR_OUT=" << output_file << endl;
    CHGSUM(chgcar_files,output_file,oss);
    return oss.str();
  }
} // namespace pflow

// ***************************************************************************
// pflow::CHGSUM
// ***************************************************************************
namespace pflow {
  bool CHGSUM(const string& chgcar_in1,const string& chgcar_in2,ostream& oss) {
    //2 INPUTS, NO OUTPUT
    string output_file="aflow_chgsum.out";
    return CHGSUM(chgcar_in1,chgcar_in2,output_file,oss);
  }
} // namespace pflow

// ***************************************************************************
// pflow::CHGSUM
// ***************************************************************************
namespace pflow {
  bool CHGSUM(string& species_header,const string& chgcar_in1,const string& chgcar_in2,const string& output_file,ostream& oss) {
    //2 INPUTS WITH SPECIES_HEADER
    vector<string> chgcar_files;
    chgcar_files.push_back(chgcar_in1);
    chgcar_files.push_back(chgcar_in2);
    return CHGSUM(species_header,chgcar_files,output_file,oss);
  }
} // namespace pflow

// ***************************************************************************
// pflow::CHGSUM
// ***************************************************************************
namespace pflow {
  bool CHGSUM(const string& chgcar_in1,const string& chgcar_in2,const string& output_file,ostream& oss) {
    //2 INPUTS, NO SPECIES_HEADER
    string species_header;
    return CHGSUM(species_header,chgcar_in1,chgcar_in2,output_file,oss);
  }
} // namespace pflow

// ***************************************************************************
// pflow::CHGSUM
// ***************************************************************************
namespace pflow {
  bool CHGSUM(const vector<string>& chgcar_files,ostream& oss) {
    //VECTOR INPUT, NO SPECIES_HEADER OR OUTPUT
    string species_header;
    return CHGSUM(species_header,chgcar_files,oss);
  }
} // namespace pflow

// ***************************************************************************
// pflow::CHGSUM
// ***************************************************************************
namespace pflow {
  bool CHGSUM(const vector<string>& chgcar_files,const string& output_file,ostream& oss) {
    //VECTOR INPUT, NO SPECIES_HEADER
    string species_header;
    return CHGSUM(species_header,chgcar_files,output_file,oss);
  }
} // namespace pflow

// ***************************************************************************
// pflow::CHGSUM
// ***************************************************************************
namespace pflow {
  bool CHGSUM(string& species_header,const vector<string>& chgcar_files,ostream& oss) {
    //VECTOR INPUT, NO OUTPUT
    string output_file="aflow_chgsum.out";
    return CHGSUM(species_header,chgcar_files,output_file,oss);
  }
} // namespace pflow

// ***************************************************************************
// pflow::CHGSUM
// ***************************************************************************
namespace pflow {
  bool CHGSUM(string& species_header,const vector<string>& chgcar_files,const string& output_file,ostream& oss) {
    // RETURNS CHGCAR_OUT=\sum CHGCAR_INPUTs 
    // Read in CHGCAR or AECCAR files
    // format_dim is as follows: numcolumns chg_tot, npts and numcolumns for 
    // augmentation occupancies, numcolumns chg_diff, npts and numcolumns for
    // augmentation occupancies

    double TOL=1e-5;
    xstructure structure1,structure2;
    stringstream chgcar_ss,chgcar1_header,chgcar2_header;
    string soliloquy="pflow::CHGSUM():  ";     // so you know who's talking

    // DEBUG
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) oss << soliloquy << "BEGIN FUNCTION" << endl; 

    for(uint i=0;i<chgcar_files.size();i++) {
      if(!aurostd::FileExist(chgcar_files.at(i))) {
        oss << endl;
        oss << soliloquy << "ERROR: "<< chgcar_files.at(i) << " does not exist." << endl;
        oss << soliloquy << "Exiting." << endl;
        oss << endl;
        return FALSE;
      }
      if(LDEBUG) oss << soliloquy << "CHGCAR_IN_"<< i+1 << "=" << chgcar_files.at(i) << endl;
    }

    // read first chgcar
    vector<int> ngrid1(3),ngrid2(3),format_dim1,format_dim2;
    vector<double> chg_tot1,chg_tot2;
    vector<double> chg_diff1,chg_diff2;
    pflow::matrix<double> lat1,lat2;

    if(LDEBUG) oss << soliloquy << "CHECK " << chgcar_files.at(0) << endl;
    aurostd::efile2stringstream(chgcar_files.at(0),chgcar_ss);
    if(!pflow::ReadCHGCAR(structure1,chgcar1_header,ngrid1,format_dim1,chg_tot1,chg_diff1,chgcar_ss,oss)) {
      oss << endl;
      oss << soliloquy << "ERROR: Input " << chgcar_files.at(0) << " format not recognized." << endl;
      oss << soliloquy << "Exiting." << endl;
      oss << endl;
      return FALSE;       // CHGCAR format not recognized
    }
    if(LDEBUG) oss << soliloquy << "SUCCESSFULLY READ " << chgcar_files.at(0) << endl;

    // for later checks
    lat1=pflow::GetLat(structure1);

    // scroll through other structures
    for(uint i=1;i<chgcar_files.size();i++) {
      chgcar_ss.str("");
      //            structure2.~xstructure();
      chgcar2_header.str("");
      ngrid2.clear();
      format_dim2.clear();
      chg_tot2.clear();
      chg_diff2.clear();

      if(LDEBUG) oss << soliloquy << "CHECK " << chgcar_files.at(i) << endl;
      aurostd::efile2stringstream(chgcar_files.at(i),chgcar_ss);
      if(!pflow::ReadCHGCAR(structure2,chgcar2_header,ngrid2,format_dim2,chg_tot2,chg_diff2,chgcar_ss,oss)) {
        oss << endl;
        oss << soliloquy << "ERROR: Input " << chgcar_files.at(i) << " format not recognized." << endl;
        oss << soliloquy << "Exiting." << endl;
        oss << endl;
        return FALSE;       // CHGCAR format not recognized
      }
      if(LDEBUG) oss << soliloquy << "SUCCESSFULLY READ CHGCAR_INPUT_" << i << endl;
      lat2=pflow::GetLat(structure2);

      // check formats
      if(LDEBUG) oss << soliloquy << "CHECK IF FORMATS OF " << chgcar_files.at(0) << " and " << chgcar_files.at(i) << " MATCH" << endl;
      if(format_dim1!=format_dim2) {
        oss << endl;
        oss << soliloquy << "WARNING: Format for " << chgcar_files.at(0) << " does not match that of " << chgcar_files.at(i) << "." << endl;
        oss << soliloquy << "WARNING: Using format from " << chgcar_files.at(0) << "." << endl;
        oss << endl;
      }
      if(LDEBUG) oss << soliloquy << "CHECK IF GRIDS OF " << chgcar_files.at(0) << " and " << chgcar_files.at(i) << " MATCH" << endl;
      if(ngrid1!=ngrid2) {
        oss << endl;
        oss << soliloquy << "ERROR: Grid for " << chgcar_files.at(0) << " does not match that of " << chgcar_files.at(i) <<  "." << endl;
        oss << soliloquy << "ERROR: ngrid of " << chgcar_files.at(0) << ": " << ngrid1.at(0) << " " <<  ngrid1.at(1) << " " <<  ngrid1.at(2) << endl;
        oss << soliloquy << "ERROR: ngrid of " << chgcar_files.at(i) << ": " << ngrid2.at(0) << " " <<  ngrid2.at(1) << " " <<  ngrid2.at(2) << endl;
        oss << soliloquy << "ERROR: This will give nonsense." << endl;
        oss << soliloquy << "Exiting." << endl;
        oss << endl;
        return FALSE;
      }
      if(LDEBUG) oss << soliloquy << "CHECK IF LATTICE PARAMETERS OF " << chgcar_files.at(0) << " and " << chgcar_files.at(i) << " MATCH" << endl;
      if(pflow::norm(pflow::VVdiff(lat1[0],lat2[0]))>TOL
         || pflow::norm(pflow::VVdiff(lat1[1],lat2[1]))>TOL
         || pflow::norm(pflow::VVdiff(lat1[2],lat2[2]))>TOL) {
        oss << endl;
        oss << soliloquy << "WARNING: Lattice parameters for " << chgcar_files.at(0) << " and " << chgcar_files.at(i) << " do not match. " << endl;
        oss << soliloquy << "WARNING: lattice of " << chgcar_files.at(0) << ": " << endl;
        oss << soliloquy << matrix2xmatrix(lat1) << endl;
        oss << soliloquy << "WARNING: lattice " << chgcar_files.at(i) << ": " << endl;
        oss << matrix2xmatrix(lat2) << endl;
        oss << soliloquy << "WARNING: This could give nonsense if there is much difference." << endl;
        oss << soliloquy << "WARNING: Output will use lattice of " << chgcar_files.at(0) << "." << endl;
        oss << endl;
      }
      if(LDEBUG) oss << soliloquy << "OVERALL FORMAT OF " << chgcar_files.at(i) << " LOOKS OK" << endl;

      // Get sum
      chg_tot1=pflow::VVsum(chg_tot1,chg_tot2);
      chg_diff1=pflow::VVsum(chg_diff1,chg_diff2);
    }

    if(LDEBUG) oss << soliloquy << "PRINTING CHGSUM TO " << output_file << endl;

    //EDIT CHGCAR_HEADER1 FORMATTING FOR BADER
    if(!species_header.empty()) {bader_functions::adjust_header(species_header,chgcar1_header,oss);}

    //print chgcar
    pflow::PrintCHGCAR(structure1,chgcar1_header,ngrid1,format_dim1,chg_tot1,chg_diff1,output_file,oss);
    if(LDEBUG) oss << soliloquy << "DONE" << endl;
    oss << soliloquy << output_file << " generated." << endl;
    return TRUE;
  }
} // namespace pflow

// ***************************************************************************
// pflow::CIF
// ***************************************************************************
namespace pflow {
  void CIF(istream& input) {
    xstructure a(input,IOAFLOW_AUTO);
    pflow::PrintCIF(cout,a);
  }
} // namespace pflow

// ***************************************************************************
// pflow::CLAT
// ***************************************************************************
namespace pflow {
  void CLAT(string options) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "pflow::CLAT: BEGIN" << endl;
    vector<string> tokens;
    aurostd::string2tokens(options,tokens,",");
    if(tokens.size()!=6) {
      init::ErrorOption(cout,options,"pflow::CLAT","aflow --clat=a,b,c,alpha,beta,gamma");
      exit(0);
    }
    xvector<double> data(6);
    if(tokens.size()>=1) data[1]=aurostd::string2utype<double>(tokens.at(0));
    if(tokens.size()>=2) data[2]=aurostd::string2utype<double>(tokens.at(1));
    if(tokens.size()>=3) data[3]=aurostd::string2utype<double>(tokens.at(2));
    if(tokens.size()>=4) data[1]=aurostd::string2utype<double>(tokens.at(3));
    if(tokens.size()>=5) data[2]=aurostd::string2utype<double>(tokens.at(4));
    if(tokens.size()>=6) data[3]=aurostd::string2utype<double>(tokens.at(5));
    cout << aflow::Banner("BANNER_TINY") << endl;
    pflow::PrintClat(data,cout);
    if(LDEBUG) cerr << "pflow::CLAT: END" << endl;
  }
} // namespace pflow

// ***************************************************************************
// pflow::CLEANALL
// ***************************************************************************
namespace pflow {
  void CLEANALL(istream& input) {
    vector<string> vinput;
    aurostd::stream2vectorstring(input,vinput);
    for(uint iinput=0;iinput<vinput.size();iinput++)
      cout << vinput.at(iinput) << endl;
  }
} // namespace pflow

// ***************************************************************************
// pflow::CMPSTR
// ***************************************************************************
namespace pflow {
  void CMPSTR(vector<string> argv) {
    cout << aflow::Banner("BANNER_TINY") << endl;
    // Read in input file.
    ifstream infile1(argv.at(2).c_str());
    ifstream infile2(argv.at(3).c_str());
    double rcut=atof(argv.at(4).c_str());
    xstructure str1;
    xstructure str2;
    infile1 >> str1;
    infile2 >> str2;
    pflow::PrintCmpStr(str1,str2,rcut,cout);
  }
} // namespace pflow

// ***************************************************************************
// pflow::COMPARE
// ***************************************************************************
namespace pflow {
  void COMPARE(string options) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "pflow::COMPARE: BEGIN" << endl;  
    vector<string> tokens;
    aurostd::string2tokens(options,tokens,",");
    if(tokens.size()!=12) {
      init::ErrorOption(cout,options,"pflow::COMPARE","aflow --compare=a,b,c,d,e,f,g,h,k,j,i,l");
      exit(0);
    }
    xvector<double> aa(12);
    for(int i=1;i<=12;i++) 
      aa(i)=aurostd::string2utype<double>(tokens.at(i-1)); 
    double bb1=0.0,bb2=0.0,bb3=0.0,bb4=0.0,bb5=0.0,bb6=0.0;
    cout << "PERCENTAGES" << "  ";
    cout << (bb1=100*abs(aa(1)-aa(7))/((aa(1)+aa(7))/2.0)) << "  ";
    cout << (bb2=100*abs(aa(2)-aa(8))/((aa(2)+aa(8))/2.0)) << "  ";
    cout << (bb3=100*abs(aa(3)-aa(9))/((aa(3)+aa(9))/2.0)) << "  ";
    cout << (bb4=100*abs(aa(4)-aa(10))/((aa(4)+aa(10))/2.0)) << "  ";
    cout << (bb5=100*abs(aa(5)-aa(11))/((aa(5)+aa(11))/2.0)) << "  ";
    cout << (bb6=100*abs(aa(6)-aa(12))/((aa(6)+aa(12))/2.0)) << "  ";
    cout << endl;
    cout << "MIN " << min(bb1,min(bb2,min(bb3,min(bb4,min(bb5,bb6))))) << endl;
    cout << "MAX " << max(bb1,max(bb2,max(bb3,max(bb4,max(bb5,bb6))))) << endl;
    // exit(0);
  }
} // namespace pflow

// ***************************************************************************
// pflow::CORNERS
// ***************************************************************************
// Stefano Curtarolo (Dec-2009)
namespace pflow {
  xstructure CORNERS(istream& input) {
    xstructure a(input,IOAFLOW_AUTO);
    a.AddCorners();
    return a;

    if(0) {
      xstructure a(input,IOAFLOW_AUTO),b;
      a.BringInCell();b=a;
      while(a.atoms.size()) a.RemoveAtom(0); 
      for(uint iat=0;iat<b.atoms.size();iat++) {
	for(double i=0;i<=1;i+=0.99) {
	  for(double j=0;j<=1;j+=0.99) {
	    for(double k=0;k<=1;k+=0.99) {
	      _atom atom=b.atoms.at(iat);
	      atom.fpos[1]+=i;atom.fpos[2]+=j;atom.fpos[3]+=k;
	      atom.cpos=F2C(a.lattice,atom.fpos);
	      if(atom.fpos[1]<=1.0 && atom.fpos[2]<=1.0 && atom.fpos[3]<=1.0) a.AddAtom(atom);
	      //	    if(aurostd::isequal(atom.fpos[1],1.0,0.02) && atom.fpos[2]<1.0 && atom.fpos[3]<1.0)   a.AddAtom(atom);	    //a.AddAtom(atom);
	    }
	  }
	}
      }
      a.title=a.title+" with_corners";
      return a;
    }
    return a;
  }
} // namespace pflow

// ***************************************************************************
// pflow::CPRIM
// ***************************************************************************
// Stefano Curtarolo (Dec-2009)
namespace pflow {
  xstructure CPRIM(istream& input) {
    cerr << "pflow::CPRIM: THIS IS A DEBUG FUNCTION FOR CODING PURPOSES" << endl;
    xstructure str_in(input,IOAFLOW_AUTO);
    xstructure str_sp,str_sc;
    // str_in.SetCoordinates(_COORDS_CARTESIAN_);
    //DX 8/29/17 [OBSOLETE] LATTICE::Standard_Lattice_StructureDefault(str_in,str_sp,str_sc);
    bool full_sym=false; //DX 8/29/17 - Speed increase
    LATTICE::Standard_Lattice_Structure(str_in,str_sp,str_sc,full_sym); //DX 8/29/17 - Speed increase
    cerr << "ORIGINAL" << endl;
    cerr << Getabc_angles(str_in.lattice,DEGREES) << endl;
    cerr << str_in;
    cerr << "STANDARD_PRIMITIVE" << endl;
    cerr << Getabc_angles(str_sp.lattice,DEGREES) << endl;
    cerr << str_sp;
    cerr << "STANDARD_CONVENTIONAL" << endl;
    cerr << Getabc_angles(str_sc.lattice,DEGREES) << endl;
    cerr << str_sc;
    exit(0);
    return str_sp;
  }
} // namespace pflow

// ***************************************************************************
// pflow::DATA
// ***************************************************************************
namespace pflow {
  //void DATA(string smode,istream& input) {
  bool DATA(string smode, istream& input, aurostd::xoption& vpflow, ostream& oss) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string aliases = "";
    string options = "";
    if(smode == "EDATA"){
      aliases = "--edata";
      options = vpflow.getattachedscheme("EDATA");
    }
    else if(smode == "DATA"){
      aliases = "--data";
      options = vpflow.getattachedscheme("DATA");
    }
    vector<string> tokens;
    aurostd::string2tokens(options,tokens,",");
    if(tokens.size()==1) {
      if(tokens.at(0)=="usage" || tokens.at(0)=="USAGE") {
        return init::ErrorOption(cout,options,"pflow::DATA",
                          aurostd::liststring2string("aflow "+aliases+"[=tolerance| =tight| =loose] [--no_scan] [--print=txt| =json] < POSCAR  default: tolerance=(minimum_interatomic_distance)/100.0, print=txt"));
      }
    }
    xstructure _a(input,IOAFLOW_AUTO);
    xstructure a(_a);
    a.ReScale(1.0);
    string print_directory = " [dir=" + a.directory + "]";
    //DX 11/28/17 - MAGNETIC SYMMETRY - START
    if(vpflow.flag("DATA::MAGNETIC")){
      string magmom_info = vpflow.getattachedscheme("DATA::MAGNETIC");
      int num_atoms=a.atoms.size();
      bool is_noncoll=false; 
      vector<xvector<double> > vmag_noncoll;                                          //DX 12/5/17 - added non-collinear
      if(GetNonCollinearMagneticInfo(num_atoms,magmom_info,vmag_noncoll)){            //DX 12/5/17 - added non-collinear
        if(LDEBUG){cerr << "pflow::DATA: Non-collinear spin system detected." << endl;} //DX 12/5/17 - added non-collinear
        is_noncoll = true;
        if(!AddSpinToXstructure(a,vmag_noncoll)){
          return false;
        } 
      }                                                                               //DX 12/5/17 - added non-collinear
      bool is_coll=false; 
      vector<double> vmag;
      if(!is_noncoll){
        if(GetCollinearMagneticInfo(num_atoms,magmom_info,vmag)){
          if(LDEBUG){cerr << "pflow::DATA: Collinear spin system detected." << endl;} //DX 12/5/17 - added non-collinear
          is_coll = true;
          if(!AddSpinToXstructure(a,vmag)){
            return false;
          } 
        } 
      }
      if(!is_noncoll && !is_coll){                                                                                 //DX 12/5/17 - added non-collinear
        cerr << "pflow::DATA: ERROR: Could not detect collinear or non-collinear spin(s). Check spin input." << endl;//DX 12/5/17 - added non-collinear
        return false;                                                                                              //DX 12/5/17 - added non-collinear
      }
    }
    //DX 11/28/17 - MAGNETIC SYMMETRY - END
    double default_tolerance=SYM::defaultTolerance(a);
    double tolerance = AUROSTD_NAN;
    if(vpflow.flag("DATA::TOLERANCE")){
      string tolerance_string = vpflow.getattachedscheme("DATA::TOLERANCE");
      if(tolerance_string[0] == 't' || tolerance_string[0] == 'T'){ //Tight
        tolerance=default_tolerance;
      }
      else if(tolerance_string[0] == 'l' || tolerance_string[0] == 'L'){ //Loose
        tolerance=default_tolerance*10.0;
      }
      else{
        tolerance=aurostd::string2utype<double>(vpflow.getattachedscheme("DATA::TOLERANCE"));
      }
    }
    else{
      tolerance = default_tolerance;
    }
    if(tolerance < 1e-10){
      cerr << "ERROR: Tolerance cannot be zero (i.e. less than 1e-10) " << print_directory << "." << endl;
      return 0;
    }
    // DX 8/3/17 - Add format flag - START
    string format = "txt";
    if(XHOST.vflag_control.flag("PRINT_MODE::TXT")){
        format = "txt";
      }
    else if(XHOST.vflag_control.flag("PRINT_MODE::JSON")){
        format = "json";
    }
    else{ //default is txt
      format = "txt";
    }   
    bool no_scan = false;
    if(vpflow.flag("DATA::NO_SCAN")){
      no_scan = true;
    }
    if(format=="txt"){
      cout << aflow::Banner("BANNER_TINY") << endl;
    }
    // pflow::PrintData(a,cerr,smode);
    pflow::PrintData(a,oss,smode,tolerance,no_scan,format); //DX cout to oss
    return true;
  }
} // namespace pflow

// ***************************************************************************
// pflow::DATA1
// ***************************************************************************
namespace pflow {
  void DATA1(string options,istream& input) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "pflow::DATA1: BEGIN" << endl;
    vector<string> tokens;
    aurostd::string2tokens(options,tokens,",");
    if(tokens.size()!=1) {
      init::ErrorOption(cout,options,"pflow::DATA1","aflow --data1=rcut < POSCAR");
      exit(0);
    }
    double rcut=0.0;
    if(tokens.size()>=1) rcut=aurostd::string2utype<double>(tokens.at(0));
    
    cout << aflow::Banner("BANNER_TINY") << endl;
    xstructure a(input,IOAFLOW_AUTO);
    // Read in input file.
    // cerr << rcut << endl;exit(0);
    pflow::PrintData1(a,rcut,cout);
    if(LDEBUG) cerr << "pflow::DATA1: END" << endl;
  }
} // namespace pflow

// ***************************************************************************
// pflow::DATA2
// ***************************************************************************
namespace pflow {
  void DATA2(istream& input) {
    xstructure a(input,IOAFLOW_AUTO);
    cout << aflow::Banner("BANNER_TINY") << endl;
    pflow::PrintData2(a,cout);
  }
} // namespace pflow

// ***************************************************************************
// pflow::DISP
// ***************************************************************************
namespace pflow {
  void DISP(string options,istream& input) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "pflow::DISP: BEGIN" << endl;
    vector<string> tokens;
    aurostd::string2tokens(options,tokens,",");
    if(tokens.size()!=1) {
      init::ErrorOption(cout,options,"pflow::DISP","aflow --disp=cutoff < POSCAR");
      exit(0);
    }
    double cutoff=0.0;
    if(tokens.size()>=1) cutoff=aurostd::string2utype<double>(tokens.at(0)); 
    xstructure a(input,IOAFLOW_AUTO);
    cout << aflow::Banner("BANNER_TINY") << endl;
    pflow::PrintDisplacements(a,cutoff,cout);
    if(LDEBUG) cerr << "pflow::DISP: END" << endl;
  }
} // namespace pflow

// ***************************************************************************
// pflow::DIST
// ***************************************************************************
namespace pflow {
  void DIST(string options,istream& input) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "pflow::DIST: BEGIN" << endl;
    vector<string> tokens;
    aurostd::string2tokens(options,tokens,",");
    if(tokens.size()!=1) {
      init::ErrorOption(cout,options,"pflow::DIST","aflow --dist=cutoff < POSCAR");
      exit(0);
    }
    double cutoff=0.0;
    if(tokens.size()>=1) cutoff=aurostd::string2utype<double>(tokens.at(0));
    xstructure a(input,IOAFLOW_AUTO);
    cout << aflow::Banner("BANNER_TINY") << endl;
    pflow::PrintDistances(a,cutoff,cout);
    if(LDEBUG) cerr << "pflow::DIST: END" << endl;
  }
} // namespace pflow

// // ***************************************************************************
// // pflow::DYNADIEL // CAMILO
// // ***************************************************************************
// namespace pflow {
//   void DYNADIEL(vector<string>& argv) {  // loop pflow::DYNADIEL
//     // modify this for the spectrum
//     string outcar ;
//     xvector<double> real(3), imag(3) ;
// 
//     if(argv.size() != 3)
//       {  // user control lines - aflow specific.
//    init::ErrorOption(cout,"","pflow::DYNADIEL","aflow --dynadiel OUTCAR*");
// 	exit(0) ;
//       }
//     outcar = argv.at(2) ;
//     KBIN::GetDynaDiel(outcar, real, imag) ;
//   }  // loop pflow::DYNADIEL
// }

// ***************************************************************************
// pflow::EDOS
// ***************************************************************************
namespace pflow {
  void EDOS(vector<string> argv) {
    // 2008 wahyu setyawan    ws26@duke.edu
    // 2008 roman chepulskyy  rc74@duke.edu
#define _MAXSPECIES_ 10
    /*
      cout<<" udos -h            : help\n"
      <<" udos 1 s p         : nonspin, includes s and p orbitals\n"
      <<" udos 2 s p d f     : spin-polarized, includes s,p,d,and f orbitals.\n\n"
      <<" OUTPUT format:\n"
      <<" udos 1\n"
      <<"   energy DOS\n"
      <<" udos 2\n"
      <<"   energy DOSup DOSdown\n"
      <<" udos 2 f\n"
      <<"   energy fDOSup_elemnt1 fDOSdown_elemnt1 ... fDOSup_elmntN fDOSdown_elmntN\n"
      <<" udos 2 s d\n"
      <<"   energy sDOSup_elmnt1 sDOSdown_elmnt1 ... sDOSup_elmntN sDOSdown_elmntN dDOSup_elmnt1 dDOSdown_elmnt1 .. dDOSup_elmntN dD\
      OSdown_elmntN"
      <<" udos 2 d s\n"
      <<"   energy sDOSup_elmnt1 sDOSdown_elmnt1 ... sDOSup_elmntN sDOSdown_elmntN dDOSup_elmnt1 dDOSdown_elmnt1 .. dDOSup_elmntN dD\
      OSdown_elmntN"
      <<" udos 2 s p d f\n"
      <<"   energy DOSup DOSdown sDOSup_elmnt1 sDOSdown_elmnt1 ...sDOSup_elmntN sDOSdown_elmntN ... fDOSup_elmntN fDOSdown_elmntN\n\\
      n"
      <<" "
      <<" note: DOS for spin down is given in (negative) sign.\n"
      <<"       Splitting of the elements(or species) is according to POSCAR.\n";
    */

    int argc=argv.size();
    if(argc==2) return;
    if(!(atoi(&argv.at(2)[0])==1 or atoi(&argv.at(2)[0])==2)) return;

    bool DO_S,DO_P,DO_D,DO_F;
    int i,j,itmp,ispin,natom,nE,ncol,k,ioffset,iorb,maxOrbital=0;
    int nspec,species[_MAXSPECIES_+1];
    float minE,maxE,Efermi,ftmp;
    string tmpstr;
    ifstream pin,din;

    ispin = atoi(&argv.at(2)[0]);
    if( !(ispin==1 or ispin==2)) return;

    DO_S=false;
    DO_P=false;
    DO_D=false;
    DO_F=false;

    if(argc>3) {
      for(i=3;i<argc;i++) {
	if(argv.at(i)[0]=='s' or argv.at(i)[0]=='S') DO_S=true;
	if(argv.at(i)[0]=='p' or argv.at(i)[0]=='P') DO_P=true;
	if(argv.at(i)[0]=='d' or argv.at(i)[0]=='D') DO_D=true;
	if(argv.at(i)[0]=='f' or argv.at(i)[0]=='F') DO_F=true;
      }
    }
    //getting number of each species from POSCAR to accumulate s,p,d,f atomic DOS
    pin.open("POSCAR");
    for(i=1;i<=5;i++) getline(pin,tmpstr);
    i=0;j=0;
    pin>>tmpstr;
    do {
      species[++i]=atoi(&tmpstr[0]);
      pin>>tmpstr;
    } while(atoi(&tmpstr[0])>0);
    nspec=i;
    pin.close();
    //processing DOSCAR
    din.open("DOSCAR");
    din>>natom;
    itmp=0;
    for(i=1;i<=nspec;i++) itmp+=species[i];
    if(itmp!=natom) {cerr<<"DOSCAR is INcompatible with POSCAR\n";return;}
    getline(din,tmpstr);
    getline(din,tmpstr);
    getline(din,tmpstr);
    getline(din,tmpstr);
    getline(din,tmpstr);
    din>>maxE>>minE>>nE>>Efermi;   getline(din,tmpstr);
    //energy loop for DOS
    //only DOS will be considered and written out because the Integrated DOS
    //can be calculated from DOS and energy bins
    ncol=1;
    if(DO_S) ncol=1+1*nspec;
    if(DO_P) ncol=1+2*nspec;
    if(DO_D) ncol=1+3*nspec;
    if(DO_F) ncol=1+4*nspec;
    if(ispin) ncol=ncol*2;
    ncol++;
    xmatrix<float> EDOS(nE,ncol);
    for(i=1;i<=nE;i++) {
      din>>EDOS[i][1]>>EDOS[i][2]; if(ispin==2) din>>EDOS[i][3];
      getline(din,tmpstr);
    }
    //energy loop for DOS for s,p,d,f
    //sum over all atoms for the same SPECIES!
    //note that we read up to the highest between s,p,d,f
    //even though not all of them will be outputed
    //We will output only the orbitals that are requested
    //at the prompt input.
    if(DO_S) maxOrbital=1;
    if(DO_P) maxOrbital=2;
    if(DO_D) maxOrbital=3;
    if(DO_F) maxOrbital=4;
    if(maxOrbital>0) {
      ioffset=0;
      for(i=1;i<=nspec;i++) {
	ioffset=(i-1)*maxOrbital*ispin;
	for(j=1;j<=species[i];j++) {
	  getline(din,tmpstr);
	  for(k=1;k<=nE;k++) {
	    din>>ftmp;//discard energy
	    if(ispin==1) {
	      for(iorb=1;iorb<=maxOrbital;iorb++) {
		din>>ftmp; EDOS[k][ioffset+iorb+2]+=ftmp;//s
	      }
	      getline(din,tmpstr);
	    }
	    if(ispin==2) {
	      for(iorb=1;iorb<=maxOrbital;iorb++) {
		din>>ftmp; EDOS[k][ioffset+(iorb-1)*2+4]+=ftmp;//up
		din>>ftmp; EDOS[k][ioffset+(iorb-1)*2+5]-=ftmp;//down
	      }
	      getline(din,tmpstr);
	    }
	  }
	}
      }
    }//if maxOrbital>0
    //output
    for(i=1;i<=nE;i++) {
      cout<<"  "<<EDOS[i][1]; //energy
      for(j=1;j<=ispin;j++)
	cout<<"  "<<EDOS[i][j+1]; //DOS
      if(DO_S) {
	//      cerr<<"do S\n";
	for(j=1;j<=nspec;j++) {
	  ioffset=(j-1)*maxOrbital*ispin;
	  for(k=1;k<=ispin;k++)
	    cout<<"  "<<EDOS[i][ioffset+k+1+ispin];
	}
      }
      if(DO_P) {
	//cerr<<"do P\n";
	for(j=1;j<=nspec;j++) {
	  ioffset=(j-1)*maxOrbital*ispin+ispin;
	  for(k=1;k<=ispin;k++)
	    cout<<"  "<<EDOS[i][ioffset+k+1+ispin];
	}
      }
      if(DO_D) {
	//cerr<<"do D\n";
	for(j=1;j<=nspec;j++) {
	  ioffset=(j-1)*maxOrbital*ispin+2*ispin;
	  for(k=1;k<=ispin;k++)
	    cout<<"  "<<EDOS[i][ioffset+k+1+ispin];
	}
      }
      if(DO_F) {
	//cerr<<"do F\n";
	for(j=1;j<=nspec;j++) {
	  ioffset=(j-1)*maxOrbital*ispin+3*ispin;
	  for(k=1;k<=ispin;k++)
	    cout<<"  "<<EDOS[i][ioffset+k+1+ispin];
	}
      }
      cout<<endl;
    }
  }
} // namespace pflow

// [OBSOLETE] // ***************************************************************************
// [OBSOLETE] // pflow::EffectiveMass
// [OBSOLETE] // ***************************************************************************
// [OBSOLETE] namespace pflow {
// [OBSOLETE]   bool EffectiveMass(vector<string> &argv,string directory_name,ostream& oss) {
// [OBSOLETE]     // aflow --effective_mass directory_name
// [OBSOLETE]     // OR
// [OBSOLETE]     // aflow --em directory_name
// [OBSOLETE]     if(argv.size()) {;} // phony just to keep argv busy no complaining about unused
// [OBSOLETE]   
// [OBSOLETE]     // the data should come from the .static step of an AFlow calculation
// [OBSOLETE]     string suffix=".static";
// [OBSOLETE]   
// [OBSOLETE]     // instantiate and load xOUTCAR class
// [OBSOLETE]     xOUTCAR outcar;
// [OBSOLETE]     string outcar_path=directory_name + "OUTCAR" + suffix + ".bz2";
// [OBSOLETE]     stringstream ss_outcar;
// [OBSOLETE]     if(!aurostd::bz2file2stringstream(outcar_path, ss_outcar)) {
// [OBSOLETE]       cerr << "" << endl;
// [OBSOLETE]       return FALSE;
// [OBSOLETE]     }
// [OBSOLETE]     outcar.GetProperties( ss_outcar );
// [OBSOLETE]   
// [OBSOLETE]     // instantiate and load xEIGENVAL class
// [OBSOLETE]     xEIGENVAL eigenval;
// [OBSOLETE]     // assume the file is compressed with bzip2
// [OBSOLETE]     string eigenval_path=directory_name + "EIGENVAL" + suffix + ".bz2";
// [OBSOLETE]     stringstream ss_eigenval;
// [OBSOLETE]     if(!aurostd::bz2file2stringstream(eigenval_path, ss_eigenval)) {
// [OBSOLETE]       cerr << "With no EIGENVAL file, there is no effective mass." << endl;
// [OBSOLETE]       return FALSE;
// [OBSOLETE]     }
// [OBSOLETE]     eigenval.GetProperties( ss_eigenval );
// [OBSOLETE]   
// [OBSOLETE]     // instantiate and load xDOSCAR class
// [OBSOLETE]     xDOSCAR doscar;
// [OBSOLETE]     // assume the file is compressed with bzip2
// [OBSOLETE]     string doscar_path=directory_name + "DOSCAR" + suffix + ".bz2";
// [OBSOLETE]     stringstream ss_doscar;
// [OBSOLETE]     if(!aurostd::bz2file2stringstream(doscar_path, ss_doscar)) {
// [OBSOLETE]       cerr << "With no DOSCAR file, there is no effective mass." << endl;
// [OBSOLETE]       return FALSE;
// [OBSOLETE]     }
// [OBSOLETE]     doscar.GetProperties( ss_doscar );
// [OBSOLETE]   
// [OBSOLETE]     // instantiate and load xstructure class
// [OBSOLETE]     // assume the file is compressed with bzip2
// [OBSOLETE]     string poscar_path=directory_name + "POSCAR" + suffix + ".bz2";
// [OBSOLETE]     stringstream ss_poscar;
// [OBSOLETE]     if(!aurostd::bz2file2stringstream(poscar_path, ss_poscar)) {
// [OBSOLETE]       cerr << "With no POSCAR file, there is no effective mass." << endl;
// [OBSOLETE]       return FALSE;
// [OBSOLETE]     }
// [OBSOLETE]     xstructure xstr(ss_poscar, IOVASP_POSCAR);
// [OBSOLETE]   
// [OBSOLETE]     // WAITING FOR CAMILO    return GetEffectiveMass( outcar, doscar, eigenval, xstr, oss, TRUE);
// [OBSOLETE]   }
// [OBSOLETE] } // namespace pflow

// ***************************************************************************
// pflow::EFFMASS // CAMILO
// ***************************************************************************
namespace pflow {
  void EFFMASS(vector<string>& argv, ostream& oss) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "pflow::EFFMASS: BEGIN" << endl;
    // aflow --effective_mass directory_name
    // aflow --em             directory_name
    string WorkDir = argv.at(2) ;
    PrintEffectiveMass(WorkDir, oss) ;
    if(LDEBUG) cerr << "pflow::EFFMASS: END" << endl;
  }
}

// ***************************************************************************
// pflow::EQUIVALENT
// ***************************************************************************
namespace pflow {
  //DX 8/18/17 [OBSOLETE] xstructure EQUIVALENT(_aflags &aflags,istream& input) {
  string EQUIVALENT(_aflags &aflags,istream& input, aurostd::xoption& vpflow) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "pflow::EQUIVALENT: BEGIN" << endl;
    xstructure _a(input,IOAFLOW_AUTO);
    bool PRINT_SCREEN=FALSE;
    aflags.QUIET=TRUE;
    //  DEBUG=TRUE;
    ofstream FileMESSAGE("/dev/null");
    //DX and CO - START
    _kflags kflags;
    // DX 8/15/17 - Add in consistency checks pflow::PerformFullSymmetry(a,File,aflags,kflags,PRINT_SCREEN,cout);
    kflags.KBIN_SYMMETRY_CALCULATE_PGROUP=TRUE;       //DX 8/15/17 - Add in consistency checks
    kflags.KBIN_SYMMETRY_CALCULATE_PGROUPK=FALSE;     //DX 8/15/17 - Add in consistency checks
    kflags.KBIN_SYMMETRY_CALCULATE_FGROUP=TRUE;       //DX 8/15/17 - Add in consistency checks
    kflags.KBIN_SYMMETRY_CALCULATE_PGROUP_XTAL=TRUE;  //DX 8/15/17 - Add in consistency checks
    kflags.KBIN_SYMMETRY_CALCULATE_PGROUPK_XTAL=FALSE;//DX 8/15/17 - Add in consistency checks //DX 12/5/17 - Added pgroupk_xtal
    kflags.KBIN_SYMMETRY_CALCULATE_SGROUP=FALSE;      //DX 8/15/17 - Add in consistency checks
    kflags.KBIN_SYMMETRY_CALCULATE_IATOMS=TRUE;       //DX 8/15/17 - Add in consistency checks
    kflags.KBIN_SYMMETRY_CALCULATE_AGROUP=FALSE;      //DX 8/15/17 - Add in consistency checks
    string options = vpflow.getattachedscheme("EQUIVALENT");
    vector<string> tokens;
    aurostd::string2tokens(options,tokens,",");
    if(tokens.size()==1) {
      if(tokens.at(0)=="usage" || tokens.at(0)=="USAGE") {
        init::ErrorOption(cout,options,"pflow::EQUIVALENT",
                          aurostd::liststring2string("aflow --equivalent|--equiv|--inequivalent|--inequiv|--iatoms|--eatoms[=tolerance| =tight| =loose] [--no_scan] [--print=txt| =json] [--mag|--magnetic|--magmom=[m1,m2,...|INCAR|OUTCAR]] < POSCAR  default: tolerance=(minimum_interatomic_distance)/100.0"));
        exit(0);
      }
    }
    //DX 170804 - need to rescale, so we make a fast copy and calculate
    xstructure a(_a);
    //DX 2/21/18 - use pwd - START
    if(a.directory == ""){
      if(aflags.Directory != "./"){
        a.directory = aflags.Directory;
      }
      else{
        a.directory = aurostd::execute2string("pwd");
        aflags.Directory = a.directory;
      }
    }
    string print_directory = " [dir=" + a.directory + "]";
    //DX 2/21/18 - use pwd - END
    a.ReScale(1.0);
    //DX 9/21/17 - MAGNETIC SYMMETRY - START
    if(vpflow.flag("SYMMETRY::MAGNETIC")){
      string magmom_info = vpflow.getattachedscheme("SYMMETRY::MAGNETIC");
      int num_atoms=a.atoms.size();
      bool is_noncoll=false; 
      vector<xvector<double> > vmag_noncoll;                                          //DX 12/5/17 - added non-collinear
      if(GetNonCollinearMagneticInfo(num_atoms,magmom_info,vmag_noncoll)){            //DX 12/5/17 - added non-collinear
        if(LDEBUG){cerr << "pflow::EQUIVALENT: Non-collinear spin system detected." << endl;} //DX 12/5/17 - added non-collinear
        is_noncoll = true;
        if(!AddSpinToXstructure(a,vmag_noncoll)){
          exit(0);
        } 
      }                                                                               //DX 12/5/17 - added non-collinear
      bool is_coll=false; 
      vector<double> vmag;
      if(!is_noncoll){
        if(GetCollinearMagneticInfo(num_atoms,magmom_info,vmag)){
          if(LDEBUG){cerr << "pflow::EQUIVALENT: Collinear spin system detected." << endl;}   //DX 12/5/17 - added non-collinear
          is_coll = true;
          if(!AddSpinToXstructure(a,vmag)){
            exit(0);
          } 
        } 
      }
      if(!is_noncoll && !is_coll){                                                                                 //DX 12/5/17 - added non-collinear
        cerr << "pflow::EQUIVALENT: ERROR: Could not detect collinear or non-collinear spin(s). Check spin input." << endl;//DX 12/5/17 - added non-collinear
        exit(0);                                                                                                   //DX 12/5/17 - added non-collinear
      }
    } 
    //DX 9/21/17 - MAGNETIC SYMMETRY - END
    double default_tolerance=SYM::defaultTolerance(a);
    double tolerance = AUROSTD_NAN;
    if(vpflow.flag("SYMMETRY::TOLERANCE")){
      string tolerance_string = vpflow.getattachedscheme("SYMMETRY::TOLERANCE");
      if(tolerance_string[0] == 't' || tolerance_string[0] == 'T'){ //Tight
        tolerance=default_tolerance;
      }
      else if(tolerance_string[0] == 'l' || tolerance_string[0] == 'L'){ //Loose
        tolerance=default_tolerance*10.0;
      }
      else{
        tolerance=aurostd::string2utype<double>(vpflow.getattachedscheme("SYMMETRY::TOLERANCE"));
      }
    }
    else{
      tolerance = default_tolerance;
    }
    if(tolerance < 1e-10){
      cerr << "pflow::EQUIVALENT ERROR: Tolerance cannot be zero (i.e. less than 1e-10) " << print_directory << "." << endl;
      exit(1);
    }
    //DX NOT REALLY NEEDED FOR THIS FUNCTION
    // DX 8/3/17 - Add format flag - START
    string format = "txt";
    if(XHOST.vflag_control.flag("PRINT_MODE::TXT")){
        format = "txt";
      }
    else if(XHOST.vflag_control.flag("PRINT_MODE::JSON")){
        format = "json";
    }
    else{ //default is txt
      format = "txt";
    }
    // Perform full scan
    bool no_scan = false;
    bool force_perform = true; //if no_scan fails, still return true at default tolerance (even though it cannot be validated)
    if(vpflow.flag("SYMMETRY::NO_SCAN")) {
      no_scan=true;
    }

    if(!pflow::PerformFullSymmetry(a,tolerance,no_scan,force_perform,FileMESSAGE,aflags,kflags,PRINT_SCREEN,cout)){
      cerr << "pflow::EQUIVALENT ERROR: Could not find commensurate symmetry at tolerance = " << tolerance << " " <<  print_directory << "." << endl; 
      exit(1);
    }
    //pflow::PerformFullSymmetry(a,File,aflags,kflags,PRINT_SCREEN,cout); //DX 8/15/17 - Add in consistency checks
    //SYM::CalculatePointGroup(File,a,aflags,TRUE,PRINT_SCREEN,cout);
    //SYM::CalculateFactorGroup(File,a,aflags,TRUE,PRINT_SCREEN,cout);
    //SYM::CalculateSpaceGroup(File,a,aflags,FALSE,PRINT_SCREEN,cout);
    //SYM::CalculateInequivalentAtoms(File,a,aflags,TRUE,PRINT_SCREEN,cout);
    //DX and CO - END
    if(LDEBUG) cerr << "pflow::EQUIVALENT: END" << endl;
    stringstream oss;
    if(format == "txt"){
    a.write_inequivalent_flag=TRUE;
      stringstream oss;
      oss << a << endl;
      return oss.str();
    }
    else if(format == "json"){
      char mode = _IATOMS_;
      KBIN_SymmetryToScreen(a,format,oss,mode);
      return oss.str();
    }
    return oss.str();
  }
} // namespace pflow

// ***************************************************************************
// pflow::EWALD
// ***************************************************************************
namespace pflow {
  void EWALD(string options,istream& input) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "pflow::EWALD: BEGIN" << endl;
    vector<string> tokens;
    aurostd::string2tokens(options,tokens,",");
    if(tokens.size()>=2) {
      init::ErrorOption(cout,options,"pflow::EWALD","aflow --ewald[=eta] < POSCAR");
      exit(0);
    } 
    // move on
    double eta=-1.0;
    if(tokens.size()>=1) eta=aurostd::string2utype<double>(tokens.at(0));
    // cout << aflow::Banner("BANNER_TINY") << endl;
    // [OBSOLETE] double eta=-1.0;
    double SUMTOL=1.0e-16;
    // [OBSOLETE] if(argv.size()==3) eta=atof(argv.at(2).c_str());
    xstructure str(input,IOAFLOW_AUTO);
    str = GetNiggliStr(str);
    double epoint,ereal,erecip,eewald;
    pflow::Ewald(str,epoint,ereal,erecip,eewald,eta,SUMTOL);
    pflow::PrintEwald(str,epoint,ereal,erecip,eewald,eta,SUMTOL,cout);

    if(LDEBUG) cerr << "pflow::EWALD: END" << endl;
  }
} // namespace pflow

// ***************************************************************************
// pflow::EXTRACT_xcar
// ***************************************************************************
namespace pflow {
  string EXTRACT_xcar(_aflags &aflags,vector<string> argv,string mode,string file) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "pflow::EXTRACT_xcar: mode=" << mode << endl;
    if(LDEBUG) cerr << "pflow::EXTRACT_xcar: file=" << file << endl;
    if(argv.size()) {;} // phony just to keep argv busy no complaining about unused
    ofstream FileMESSAGE("/dev/null");
    _kflags kflags;kflags.AFLOW_MODE_VASP=TRUE;
    _vflags vflags;_xvasp xvasp;xvasp.clear();
    if(!aurostd::FileExist(file)) {cerr << "pflow::EXTRACT_xcar: mode=" << mode << "  file=" << file << " not found .." << endl;exit(0);};
    string AflowIn; aurostd::file2string(file,AflowIn);
    aflags.QUIET=TRUE;XHOST.QUIET=TRUE;
    if(mode=="POSCAR") {
      vflags=KBIN::VASP_Get_Vflags_from_AflowIN(AflowIn,aflags,kflags);
      //  OBSOLETE]  KBIN::VASP_Produce_POSCAR(xvasp,AflowIn,FileAFLOWIN,FileMESSAGE,aflags,kflags,vflags);
      KBIN::VASP_Produce_POSCAR(xvasp,AflowIn,FileMESSAGE,aflags,kflags,vflags);
      return xvasp.POSCAR.str();
    }
    if(mode=="INCAR") {
      vflags=KBIN::VASP_Get_Vflags_from_AflowIN(AflowIn,aflags,kflags);
      // [OBSOLETE]  KBIN::VASP_Produce_INCAR(xvasp,AflowIn,FileAFLOWIN,FileMESSAGE,aflags,kflags,vflags);
      KBIN::VASP_Produce_INCAR(xvasp,AflowIn,FileMESSAGE,aflags,kflags,vflags);
      KBIN::VASP_Modify_INCAR(xvasp,FileMESSAGE,aflags,kflags,vflags);
      return xvasp.INCAR.str();
    }
    if(mode=="KPOINTS") {
      vflags=KBIN::VASP_Get_Vflags_from_AflowIN(AflowIn,aflags,kflags);
      KBIN::VASP_Produce_POSCAR(xvasp,AflowIn,FileMESSAGE,aflags,kflags,vflags);
      // OBSOLETE]   KBIN::VASP_Produce_KPOINTS(xvasp,AflowIn,FileAFLOWIN,FileMESSAGE,aflags,kflags,vflags);
      KBIN::VASP_Produce_KPOINTS(xvasp,AflowIn,FileMESSAGE,aflags,kflags,vflags);
      KBIN::VASP_Modify_KPOINTS(xvasp,FileMESSAGE,aflags,vflags);
      return xvasp.KPOINTS.str();
    }
    if(mode=="POTCAR") {
      vflags=KBIN::VASP_Get_Vflags_from_AflowIN(AflowIn,aflags,kflags);
      // [OBSOLETE]    KBIN::VASP_Produce_POTCAR(xvasp,AflowIn,FileAFLOWIN,FileMESSAGE,aflags,kflags,vflags);
      KBIN::VASP_Produce_POTCAR(xvasp,AflowIn,FileMESSAGE,aflags,kflags,vflags);
      return xvasp.POTCAR.str();
    }
    return mode;// something must go out
  }
} // namespace pflow


// ***************************************************************************
// pflow::EIGCURV // CAMILO
// ***************************************************************************
namespace pflow {
  void EIGCURV(string options, ostream& oss) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "pflow::EIGCURV: BEGIN" << endl;
    vector<string> tokens;
    aurostd::string2tokens(options,tokens,",");
    if(tokens.size()!=1) {
      init::ErrorOption(cout,options,"pflow::EIGCURV","aflow --eigcurv=DIRECTORY(with bands)");
      exit(0);
    }
    string filename="";
    if(tokens.size()>=1) filename=(tokens.at(0));
    string WorkDir = filename ;
    PrintEigCurv(WorkDir, oss) ;
    if(LDEBUG) cerr << "pflow::EIGCURV: END" << endl;
  }
} // namespace pflow

// DX AND COREY - START
// ***************************************************************************
// pflow::PerformFullSymmetry
// ***************************************************************************
namespace pflow {
  bool PerformFullSymmetry(xstructure& a){
    ofstream FileMESSAGE("/dev/null");
    _aflags aflags; aflags.Directory="."; aflags.QUIET=true;
    _kflags kflags; defaultKFlags4SymWrite(kflags,false); defaultKFlags4SymCalc(kflags,true);
    bool osswrite=false;
    return PerformFullSymmetry(a,FileMESSAGE,aflags,kflags,osswrite,cout);
  }

  bool PerformFullSymmetry(xstructure& a, ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags,const bool& osswrite,ostream& oss, string format){
    double tolerance = a.sym_eps;
    bool no_scan = false;
    bool force_perform = true;  //if no_scan fails, still return true at default tolerance (even though it cannot be validated)
    return PerformFullSymmetry(a,tolerance,no_scan,force_perform,FileMESSAGE,aflags,kflags,osswrite,oss,format);
  }

  bool PerformFullSymmetry(xstructure& a, double& tolerance, bool no_scan, bool force_perform, ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags,const bool& osswrite,ostream& oss, string format){
    xstructure b(a);    //save for later
    a.ReScale(1.0);     //the nuclear option, only way to fix all of the issues with f2c/c2f/ctau/ctrasl/etc.
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    bool symmetry_commensurate = FALSE;
    ostringstream aus;

    //a.directory = aflags.Directory;
    //DX 2/21/18 - use pwd - START
    if(a.directory == ""){
      if(aflags.Directory != "./"){
      a.directory = aflags.Directory;
    }
      else{
        a.directory = aurostd::execute2string("pwd");
        aflags.Directory = a.directory;
      }
    }
    string print_directory = " [dir=" + a.directory + "]";
    //DX 2/21/18 - use pwd - END

    if(LDEBUG) {cerr << "pflow::PerformFullSymmetry: STRUCTURE" << endl;cerr << a << endl;}

    // MOVED DOWN A BIT if(!aflags.QUIET) aus << (aflags.QUIET?"":"00000  MESSAGE ") << "Symmetry: starting tolerance " << _EPS_sym_ << " " << Message(aflags,"user,host,time") << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
    if(1 || a.dist_nn_min == AUROSTD_NAN){  //CO 171024 - always recalculate min_dist (SAFE)
      if(LDEBUG) cerr << "pflow::PerformFullSymmetry: CALCULATING MIN DISTANCE" << print_directory << endl;
      a.MinDist();
      if(LDEBUG) cerr << "pflow::PerformFullSymmetry: MIN DISTANCE DONE" << print_directory << endl;
    }
    double min_dist = a.dist_nn_min;
    //DX 9/5/17 [OBSOLETE] int change_sym_count=1;    
    //if(a.sym_eps!=AUROSTD_NAN){ //Tolerance came from user or was calculated
    //  a.sym_eps;
    //}
    //CO, I changed a bit, if tolerance is specified, it should override
    if(tolerance != AUROSTD_NAN) {
      a.sym_eps=tolerance;
    } else if(!a.sym_eps_calculated || a.sym_eps==AUROSTD_NAN) {
      a.sym_eps = SYM::defaultTolerance(a);
    }
    //if(a.sym_eps == AUROSTD_NAN && tolerance == AUROSTD_NAN){
    //  a.sym_eps = SYM::defaultTolerance(a);
    //}
    //else if(!a.sym_eps_calculated && tolerance != AUROSTD_NAN){
    //  a.sym_eps = tolerance;
    //}
    //a.sym_eps = SYM::defaultTolerance(a);
    double orig_tolerance = a.sym_eps;
    if(!aflags.QUIET) aus << (aflags.QUIET?"":"00000  MESSAGE ") << "Symmetry: starting tolerance " << a.sym_eps << " " << Message(aflags,"user,host,time") << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
    
    while(symmetry_commensurate==FALSE){
      //[DX OBSOLETE] a.SpaceGroup_ITC(a.sym_eps,orig_tolerance,-1,change_sym_count,no_scan); //rescales to 1.0 internally, but doesn't affect a
      //[DX OBSOLETE] if(a.space_group_ITC == 0){
      //[DX OBSOLETE]   //DXreturn FALSE;
      //[DX OBSOLETE]   if(!no_scan){
      //[DX OBSOLETE]     cerr << "pflow::PerformFullSymmetry: ERROR: Space group routine could not find space group at given tolerance." << endl;
      //[DX OBSOLETE]     return FALSE;
      //[DX OBSOLETE]   }else{
      //[DX OBSOLETE]   cerr << "pflow::PerformFullSymmetry: WARNING: Space group routine could not find space group at given tolerance." << endl;
      //[DX OBSOLETE]     //keep going, calculate rest of properties (though they are probably bad)
      //[DX OBSOLETE]   }
      //[DX OBSOLETE] }
      a.sym_eps_calculated=true;
      // Calculate Lattice Point Group 
      if(kflags.KBIN_SYMMETRY_CALCULATE_PGROUP){ //DX 8/14/17
      if(!SYM::CalculatePointGroup(FileMESSAGE,a,aflags,kflags.KBIN_SYMMETRY_PGROUP_WRITE,osswrite,oss,format)){
        if(!no_scan){
          a.ClearSymmetry();
          //DX 9/5/17 [OBSOLETE] if(!SYM::change_tolerance(a,a.sym_eps,orig_tolerance,change_sym_count,min_dist,no_scan)){
          if(!SYM::change_tolerance(a,a.sym_eps,min_dist,no_scan)){
            a=b;  //pretty printing, unmodified structure
            if(force_perform){
              cerr << "pflow::PerformFullSymmetry: Scan failed [0]. Reverting back to original tolerance and recalculating as is (with aforementioned inconsistencies)." << print_directory << endl;
              PerformFullSymmetry(a,orig_tolerance,true,false,FileMESSAGE,aflags,kflags,osswrite,oss,format);
            }else{
              return FALSE;
            }
          }
          if(!aflags.QUIET) aus << (aflags.QUIET?"":"00000  MESSAGE ") << "PGROUP Symmetry: changing tolerance to " << a.sym_eps << " " << Message(aflags,"user,host,time") << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
          continue;
        }
      }
      } //DX 8/14/17
      if(kflags.KBIN_SYMMETRY_CALCULATE_PGROUPK){ //DX 8/14/17
      if(!SYM::CalculatePointGroupKlattice(FileMESSAGE,a,aflags,kflags.KBIN_SYMMETRY_PGROUPK_WRITE,osswrite,oss,format)){
        //cerr << "COREY TESTING POINTGROUP KLATTICE FAILED!!! exiting for now" << endl;
        //exit(1);
        if(!no_scan){
          a.ClearSymmetry();
          //DX 9/5/17 [OBSOLETE] if(!SYM::change_tolerance(a,a.sym_eps,orig_tolerance,change_sym_count,min_dist,no_scan)){
          if(!SYM::change_tolerance(a,a.sym_eps,min_dist,no_scan)){
            a=b;  //pretty printing, unmodified structure
            if(force_perform){
              cerr << "pflow::PerformFullSymmetry: Scan failed [1]. Reverting back to original tolerance and recalculating as is (with aforementioned inconsistencies)." << print_directory << endl;
              PerformFullSymmetry(a,orig_tolerance,true,false,FileMESSAGE,aflags,kflags,osswrite,oss,format);
            }else{
              return FALSE;
            }
          }
          if(!aflags.QUIET) aus << (aflags.QUIET?"":"00000  MESSAGE ") << "PGROUPK Symmetry: changing tolerance to " << a.sym_eps << " " << Message(aflags,"user,host,time") << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
          continue;
        }
      }
      } //DX 8/14/17
      // Check for identity element
      if(kflags.KBIN_SYMMETRY_CALCULATE_PGROUP && SYM::CheckForIdentity(a) == FALSE){
        if(LDEBUG){ 
          cerr << "pflow::PerformFullSymmetry: WARNING: Point group does not contain the identity element (impossible for a crystal)." << print_directory << endl;
        }
        if(!no_scan){
          a.ClearSymmetry();
          //DX 9/5/17 [OBSOLETE] if(!SYM::change_tolerance(a,a.sym_eps,orig_tolerance,change_sym_count,min_dist,no_scan)){
          if(!SYM::change_tolerance(a,a.sym_eps,min_dist,no_scan)){
            a=b;  //pretty printing, unmodified structure
            if(force_perform){
              cerr << "pflow::PerformFullSymmetry: Scan failed [2]. Reverting back to original tolerance and recalculating as is (with aforementioned inconsistencies)." << print_directory << endl;
              PerformFullSymmetry(a,orig_tolerance,true,false,FileMESSAGE,aflags,kflags,osswrite,oss,format);
            }else{
              return FALSE;
            }
          }
          if(!aflags.QUIET) aus << (aflags.QUIET?"":"00000  MESSAGE ") << "PGROUP Symmetry: changing tolerance to " << a.sym_eps << " " << Message(aflags,"user,host,time") << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
          continue;
        }
      }
      // Calculate Factor Group 
      if(kflags.KBIN_SYMMETRY_CALCULATE_FGROUP){ //DX 8/14/17
      if(!SYM::CalculateFactorGroup(FileMESSAGE,a,aflags,kflags.KBIN_SYMMETRY_FGROUP_WRITE,osswrite,oss,format)){
        if(!no_scan){
          a.ClearSymmetry();
          //DX 9/5/17 [OBSOLETE] if(!SYM::change_tolerance(a,a.sym_eps,orig_tolerance,change_sym_count,min_dist,no_scan)){
          if(!SYM::change_tolerance(a,a.sym_eps,min_dist,no_scan)){
            a=b;  //pretty printing, unmodified structure
            if(force_perform){
              cerr << "pflow::PerformFullSymmetry: Scan failed [3]. Reverting back to original tolerance and recalculating as is (with aforementioned inconsistencies)." << print_directory << endl;
              PerformFullSymmetry(a,orig_tolerance,true,false,FileMESSAGE,aflags,kflags,osswrite,oss,format);
            }else{
              return FALSE;
            }
          }
          if(!aflags.QUIET) aus << (aflags.QUIET?"":"00000  MESSAGE ") << "FGROUP Symmetry: changing tolerance to " << a.sym_eps << " " << Message(aflags,"user,host,time") << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
          continue;
        }
      }
      } //DX 8/14/17
      // Calculate Crystallographic Point Group 
      if(kflags.KBIN_SYMMETRY_CALCULATE_PGROUP_XTAL){ //DX 8/14/17
      if(!SYM::CalculatePointGroupCrystal(FileMESSAGE,a,aflags,kflags.KBIN_SYMMETRY_PGROUP_XTAL_WRITE,osswrite,oss,format)){
        if(!no_scan){
	  a.ClearSymmetry();
          //DX 9/5/17 [OBSOLETE] if(!SYM::change_tolerance(a,a.sym_eps,orig_tolerance,change_sym_count,min_dist,no_scan)){
          if(!SYM::change_tolerance(a,a.sym_eps,min_dist,no_scan)){
	    a=b;  //pretty printing, unmodified structure
            if(force_perform){
              cerr << "pflow::PerformFullSymmetry: Scan failed [4]. Reverting back to original tolerance and recalculating as is (with aforementioned inconsistencies)." << print_directory << endl;
              PerformFullSymmetry(a,orig_tolerance,true,false,FileMESSAGE,aflags,kflags,osswrite,oss,format);
            }else{
              return FALSE;
            }
	  }
	  if(!aflags.QUIET) aus << (aflags.QUIET?"":"00000  MESSAGE ") << "PGROUP_XTAL Symmetry: changing tolerance to " << a.sym_eps << " " << Message(aflags,"user,host,time") << endl;
	  aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
	  continue;
        }
      }
      } //DX 8/14/17
      // Check if a point group map was found; if not, change tolerance
      if(kflags.KBIN_SYMMETRY_CALCULATE_PGROUP_XTAL && a.point_group_Hermann_Mauguin.empty() == TRUE){ //DX 8/14/17
        if(LDEBUG){ 
          cerr << "pflow::PerformFullSymmetry: WARNING: Point group crystal operations did not match with any Hermann-Mauguin symbols. (i.e. The set of symmetry elements found are not allowed possible for a crystal.) " << print_directory << endl;;
        }
        if(!no_scan){
	  a.ClearSymmetry();
          //DX 9/5/17 [OBSOLETE] if(!SYM::change_tolerance(a,a.sym_eps,orig_tolerance,change_sym_count,min_dist,no_scan)){
          if(!SYM::change_tolerance(a,a.sym_eps,min_dist,no_scan)){
	    a=b;  //pretty printing, unmodified structure
            if(force_perform){
              cerr << "pflow::PerformFullSymmetry: Scan failed [5]. Reverting back to original tolerance and recalculating as is (with aforementioned inconsistencies)." << print_directory << endl;
              PerformFullSymmetry(a,orig_tolerance,true,false,FileMESSAGE,aflags,kflags,osswrite,oss,format);
            }else{
              return FALSE;
            }
	  }
	  if(!aflags.QUIET) aus << (aflags.QUIET?"":"00000  MESSAGE ") << "PGROUP_XTAL Symmetry: changing tolerance to " << a.sym_eps << " " << Message(aflags,"user,host,time") << endl;
	  aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
	  continue;
        }
      }
      // Check if factor group is integer multiple of the crystallographic point group; if not, change tolerance
      int multiplicity_of_primitive = 0; //DX 8/15/17 - Added consistency checks, need to initialize and calculate if we have those groups
      if(kflags.KBIN_SYMMETRY_CALCULATE_FGROUP && kflags.KBIN_SYMMETRY_CALCULATE_PGROUP_XTAL){ //DX 8/14/17
        multiplicity_of_primitive = a.fgroup.size()/a.pgroup_xtal.size();
        if(a.fgroup.size()%a.pgroup_xtal.size() != 0){ //DX 8/14/17
        if(LDEBUG){ 
          cerr << "pflow::PerformFullSymmetry: WARNING: Number of factor groups is not an integer multiple of the point group crystal." << print_directory << endl;
        }
        if(!no_scan){
	  a.ClearSymmetry();
          //DX 9/5/17 [OBSOLETE] if(!SYM::change_tolerance(a,a.sym_eps,orig_tolerance,change_sym_count,min_dist,no_scan)){
          if(!SYM::change_tolerance(a,a.sym_eps,min_dist,no_scan)){
	    a=b;  //pretty printing, unmodified structure
            if(force_perform){
              cerr << "pflow::PerformFullSymmetry: Scan failed [6]. Reverting back to original tolerance and recalculating as is (with aforementioned inconsistencies)." << print_directory << endl;
              PerformFullSymmetry(a,orig_tolerance,true,false,FileMESSAGE,aflags,kflags,osswrite,oss,format);
            }else{
              return FALSE;
            }
	  }
	  if(!aflags.QUIET) aus << (aflags.QUIET?"":"00000  MESSAGE ") << "PGROUP_XTAL Symmetry: changing tolerance to " << a.sym_eps << " " << Message(aflags,"user,host,time") << endl;
	  aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
	  continue;
        }
      }
      } //DX 8/14/17
      if(kflags.KBIN_SYMMETRY_CALCULATE_PGROUPK_XTAL){ //DX 12/5/17 - Added pgroupk_xtal
        if(!SYM::CalculatePointGroupKCrystal(FileMESSAGE,a,aflags,kflags.KBIN_SYMMETRY_PGROUPK_XTAL_WRITE,osswrite,oss,format)){ //DX 1/18/18 - PGROUPK_XTAL not PGROUPK
          if(!no_scan){
            a.ClearSymmetry();
            //DX 9/5/17 [OBSOLETE] if(!SYM::change_tolerance(a,a.sym_eps,orig_tolerance,change_sym_count,min_dist,no_scan)){
            if(!SYM::change_tolerance(a,a.sym_eps,min_dist,no_scan)){
              a=b;  //pretty printing, unmodified structure
              if(force_perform){
                cerr << "pflow::PerformFullSymmetry: Scan failed [7]. Reverting back to original tolerance and recalculating as is (with aforementioned inconsistencies)." << print_directory << endl;
                PerformFullSymmetry(a,orig_tolerance,true,false,FileMESSAGE,aflags,kflags,osswrite,oss,format);
              }else{
                return FALSE;
              }
            }
            if(!aflags.QUIET) aus << (aflags.QUIET?"":"00000  MESSAGE ") << "PGROUPK_XTAL Symmetry: changing tolerance to " << a.sym_eps << " " << Message(aflags,"user,host,time") << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
            continue;
          }
        }
      } //DX 8/14/17
      //[DX OBSOLETE]// Check if point group and space group are consistent       
      //[DX OBSOLETE]bool derivative_structure = true;
      //[DX OBSOLETE]bool space_and_point_group_match = SYM::ComparePointGroupAndSpaceGroupString(a,multiplicity_of_primitive,derivative_structure);
      //[DX OBSOLETE]if(!space_and_point_group_match && !derivative_structure){
      //[DX OBSOLETE]  if(LDEBUG){ 
      //[DX OBSOLETE]    cerr << "WARNING: Point group crystal and space group are not commensurate. " << endl;
      //[DX OBSOLETE]    cerr << "Point Group Crystal: " << a.point_group_Hermann_Mauguin;
      //[DX OBSOLETE]    if(a.space_group_ITC != 0){ 
      //[DX OBSOLETE]      cerr << " | Space Group: " << GetSpaceGroupName(a.space_group_ITC) << endl; 
      //[DX OBSOLETE]    }
      //[DX OBSOLETE]    else{
      //[DX OBSOLETE]      cerr << " | Space Group: --" << endl; 
      //[DX OBSOLETE]    }
      //[DX OBSOLETE]  }
      //[DX OBSOLETE]  if(!no_scan){
      //[DX OBSOLETE]    a.ClearSymmetry();
      //[DX OBSOLETE]    //DX 9/5/17 [OBSOLETE] if(!SYM::change_tolerance(a,a.sym_eps,orig_tolerance,change_sym_count,min_dist,no_scan)){
      //[DX OBSOLETE]    if(!SYM::change_tolerance(a,a.sym_eps,min_dist,no_scan)){
      //[DX OBSOLETE]	    a=b;  //pretty printing, unmodified structure
      //[DX OBSOLETE]      if(force_perform){
      //[DX OBSOLETE]        cerr << "Scan failed. Reverting back to original tolerance and recalculating as is (with aforementioned inconsistencies)." << endl;
      //[DX OBSOLETE]        PerformFullSymmetry(a,orig_tolerance,true,false,FileMESSAGE,aflags,kflags,osswrite,oss,format);
      //[DX OBSOLETE]      }else{
      //[DX OBSOLETE]        return FALSE;
      //[DX OBSOLETE]      }
      //[DX OBSOLETE]    }
      //[DX OBSOLETE]	 if(!aflags.QUIET) aus << (aflags.QUIET?"":"00000  MESSAGE ") << "PGROUP_XTAL Symmetry: changing tolerance to " << a.sym_eps << " " << Message(aflags,"user,host,time") << endl;
      //[DX OBSOLETE]  aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
      //[DX OBSOLETE]	  continue;
      //[DX OBSOLETE]  }
      //[DX OBSOLETE]}
      //[DX OBSOLETE]else if(!space_and_point_group_match && derivative_structure){
      //[DX OBSOLETE]  if(!aflags.QUIET) aus << (aflags.QUIET?"":"00000  MESSAGE ") << "PGROUP_XTAL WARNING: Point Group Crystal of Original Cell: " << a.point_group_Hermann_Mauguin << " | Space Group of Primitive Cell: " << GetSpaceGroupName(a.space_group_ITC) << " -- This is a derivative structure." << Message(aflags,"user,host,time") << endl;
      //[DX OBSOLETE]  aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
      //[DX OBSOLETE]}
      if(kflags.KBIN_SYMMETRY_CALCULATE_SGROUP){ //DX 8/14/17
      if(kflags.KBIN_SYMMETRY_SGROUP_RADIUS>0.0) {
        if(!aflags.QUIET) aus << (aflags.QUIET?"":"00000  MESSAGE ") << "POSCAR SGROUP: found RADIUS="<<kflags.KBIN_SYMMETRY_SGROUP_RADIUS<<" " << Message(aflags,"user,host,time") << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
      } else {
        kflags.KBIN_SYMMETRY_SGROUP_RADIUS=KBIN_SYMMETRY_SGROUP_RADIUS_DEFAULT;
        if(!aflags.QUIET) aus << (aflags.QUIET?"":"00000  MESSAGE ") << "POSCAR SGROUP: Default RADIUS="<<kflags.KBIN_SYMMETRY_SGROUP_RADIUS<<" " << Message(aflags,"user,host,time") << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
      }
      a.sgroup_radius=kflags.KBIN_SYMMETRY_SGROUP_RADIUS;
      // Calculate Space Group
      if(!SYM::CalculateSpaceGroup(FileMESSAGE,a,aflags,kflags.KBIN_SYMMETRY_SGROUP_WRITE,osswrite,oss,format)){
        if(!no_scan){
	  a.ClearSymmetry();
          //DX 9/5/17 [OBSOLETE] if(!SYM::change_tolerance(a,a.sym_eps,orig_tolerance,change_sym_count,min_dist,no_scan)){
          if(!SYM::change_tolerance(a,a.sym_eps,min_dist,no_scan)){
	    a=b;  //pretty printing, unmodified structure
            if(force_perform){
              cerr << "pflow::PerformFullSymmetry: Scan failed [8]. Reverting back to original tolerance and recalculating as is (with aforementioned inconsistencies)." << print_directory << endl;
              PerformFullSymmetry(a,orig_tolerance,true,false,FileMESSAGE,aflags,kflags,osswrite,oss,format);
            }else{
              return FALSE;
            }
	  }
	  if(!aflags.QUIET) aus << (aflags.QUIET?"":"00000  MESSAGE ") << "SGROUP Symmetry: changing tolerance to " << a.sym_eps << " " << Message(aflags,"user,host,time") << endl;
	  aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
	  continue;
        }
      }
      } //DX 8/14/17
      // Calculate inequivalent atoms
      if(kflags.KBIN_SYMMETRY_CALCULATE_IATOMS){ //DX 8/14/17
      if(!SYM::CalculateInequivalentAtoms(FileMESSAGE,a,aflags,kflags.KBIN_SYMMETRY_IATOMS_WRITE,osswrite,oss,format)){
        if(!no_scan){
	  a.ClearSymmetry();
          //DX 9/5/17 [OBSOLETE] if(!SYM::change_tolerance(a,a.sym_eps,orig_tolerance,change_sym_count,min_dist,no_scan)){
          if(!SYM::change_tolerance(a,a.sym_eps,min_dist,no_scan)){
	    a=b;  //pretty printing, unmodified structure
            if(force_perform){
              cerr << "pflow::PerformFullSymmetry: Scan failed [9]. Reverting back to original tolerance and recalculating as is (with aforementioned inconsistencies)." << print_directory << endl;
              PerformFullSymmetry(a,orig_tolerance,true,false,FileMESSAGE,aflags,kflags,osswrite,oss,format);
            }else{
              return FALSE;
            }
	  }
	  if(!aflags.QUIET) aus << (aflags.QUIET?"":"00000  MESSAGE ") << "IATOMS ATOMS: changing tolerance to " << a.sym_eps << " " << Message(aflags,"user,host,time") << endl;
	  aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
	  continue;
        }
      }
      // Check if number of equivalent atoms is consistent with cell ; if not, change tolerance
      bool iatoms_commensurate = TRUE;
      for (uint i=0; i<a.iatoms.size();i++){
        if(a.iatoms[i].size()%multiplicity_of_primitive != 0){
          iatoms_commensurate = FALSE;
          break; 
        }
      }
      if(iatoms_commensurate == FALSE){
        if(LDEBUG){ 
          cerr << "pflow::PerformFullSymmetry: WARNING: Number of equivalent atoms is not an integer multiple of the number factor groups." << print_directory << endl;
        }
        if(!no_scan){
	  a.ClearSymmetry();
          //DX 9/5/17 [OBSOLETE] if(!SYM::change_tolerance(a,a.sym_eps,orig_tolerance,change_sym_count,min_dist,no_scan)){
          if(!SYM::change_tolerance(a,a.sym_eps,min_dist,no_scan)){
	    a=b;  //pretty printing, unmodified structure
            if(force_perform){
              cerr << "pflow::PerformFullSymmetry: Scan failed [10]. Reverting back to original tolerance and recalculating as is (with aforementioned inconsistencies)." << print_directory << endl;
              PerformFullSymmetry(a,orig_tolerance,true,false,FileMESSAGE,aflags,kflags,osswrite,oss,format);
            }else{
              return FALSE;
            }
	  }
	  if(!aflags.QUIET) aus << (aflags.QUIET?"":"00000  MESSAGE ") << "IATOMS ATOMS: changing tolerance to " << a.sym_eps << " " << Message(aflags,"user,host,time") << endl;
	  aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
	  continue;
        }
      }
      } //DX 8/14/17
      // Calculate site point group
      if(kflags.KBIN_SYMMETRY_CALCULATE_AGROUP){ //DX 8/14/17
      if(!SYM::CalculateSitePointGroup(FileMESSAGE,a,aflags,kflags.KBIN_SYMMETRY_AGROUP_WRITE,osswrite,oss,format)){
        if(!no_scan){
	  a.ClearSymmetry();
          //DX 9/5/17 [OBSOLETE] if(!SYM::change_tolerance(a,a.sym_eps,orig_tolerance,change_sym_count,min_dist,no_scan)){
          if(!SYM::change_tolerance(a,a.sym_eps,min_dist,no_scan)){
	    a=b;  //pretty printing, unmodified structure
            if(force_perform){
              cerr << "pflow::PerformFullSymmetry: Scan failed [11]. Reverting back to original tolerance and recalculating as is (with aforementioned inconsistencies)." << print_directory << endl;
              PerformFullSymmetry(a,orig_tolerance,true,false,FileMESSAGE,aflags,kflags,osswrite,oss,format);
            }else{
              return FALSE;
            }
	  }
	  if(!aflags.QUIET) aus << (aflags.QUIET?"":"00000  MESSAGE ") << "AGROUP Symmetry: changing tolerance to " << a.sym_eps << " " << Message(aflags,"user,host,time") << endl;
	  aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,osswrite,oss);
	  continue;
        }
      }
      } //DX 8/14/17
      symmetry_commensurate = TRUE;  // NOTE: This may not be entirely true if no_scan=TRUE
    }
    a.ReScale(b.scale);                   //the nuclear option, only way to fix all of the issues with f2c/c2f/ctau/ctrasl/etc.
    a.lattice=b.lattice;a.scale=b.scale;  //perfect printing, we should also fix cpos of all atoms, but not really worried since we usually print fpos
    return symmetry_commensurate;
  }
}

namespace pflow {
  void defaultKFlags4SymWrite(_kflags& kflags,bool write) {
    kflags.KBIN_SYMMETRY_PGROUP_WRITE=write;
    kflags.KBIN_SYMMETRY_PGROUPK_WRITE=write;
    kflags.KBIN_SYMMETRY_FGROUP_WRITE=write;
    kflags.KBIN_SYMMETRY_PGROUP_XTAL_WRITE=write;
    kflags.KBIN_SYMMETRY_PGROUPK_XTAL_WRITE=write; //DX 12/5/17 - Added pgroupk_xtal
    kflags.KBIN_SYMMETRY_IATOMS_WRITE=write;
    kflags.KBIN_SYMMETRY_AGROUP_WRITE=write;
    kflags.KBIN_SYMMETRY_SGROUP_WRITE=write;
  }
  void defaultKFlags4SymCalc(_kflags& kflags,bool calc) {
    kflags.KBIN_SYMMETRY_CALCULATE_PGROUP=calc;
    kflags.KBIN_SYMMETRY_CALCULATE_PGROUPK=calc;
    kflags.KBIN_SYMMETRY_CALCULATE_FGROUP=calc;
    kflags.KBIN_SYMMETRY_CALCULATE_PGROUP_XTAL=calc;
    kflags.KBIN_SYMMETRY_CALCULATE_PGROUPK_XTAL=calc; //DX 12/5/17 - Added pgroupk_xtal
    kflags.KBIN_SYMMETRY_CALCULATE_IATOMS=calc;
    kflags.KBIN_SYMMETRY_CALCULATE_AGROUP=calc;
    kflags.KBIN_SYMMETRY_CALCULATE_SGROUP=calc;
  }
}

namespace pflow {
  //COMMAND LINE SYMMETRY CALCULATION, calls main function PerformFullSymmetry()!!!!!!!!!!!
  bool CalculateFullSymmetry(istream& input, aurostd::xoption& vpflow, ostream& oss){ //overload
    xstructure a(input,IOAFLOW_AUTO);
    
    //default aflags for command line
    _aflags aflags;
    aflags.Directory=".";
    aflags.QUIET=false;

    //default kflags
    _kflags kflags; defaultKFlags4SymWrite(kflags,true); defaultKFlags4SymCalc(kflags,true);

    //write output to screen
    bool osswrite=TRUE;

    return CalculateFullSymmetry(aflags,kflags,a,vpflow,osswrite,oss);
  }
  
  //COMMAND LINE SYMMETRY CALCULATION, calls main function PerformFullSymmetry()!!!!!!!!!!!
  bool CalculateFullSymmetry(_aflags &aflags, _kflags& kflags, xstructure& _a, aurostd::xoption& vpflow, bool osswrite,ostream& oss){ //main function
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string options = vpflow.getattachedscheme("FULLSYMMETRY");
    vector<string> tokens;
    aurostd::string2tokens(options,tokens,",");
    if(tokens.size()==1) {
      if(tokens.at(0)=="usage" || tokens.at(0)=="USAGE") {
        return init::ErrorOption(cout,options,"pflow::CalculateFullSymmetry",
                          aurostd::liststring2string("aflow --aflow-sym|--AFLOW-SYM|--AFLOWSYM|--aflowSYM|--aflowsym|--full_symmetry|--full_sym|--fullsym[=tolerance| =tight| =loose] [--no_scan] [--print=txt| =json] [--screen_only] [--mag|--magnetic|--magmom=[m1,m2,...|INCAR|OUTCAR]] < POSCAR  default: tolerance=(minimum_interatomic_distance)/100.0, print=txt"));
      }
    }
    //DX 170804 - need to rescale, so we make a fast copy and calculate
    xstructure a(_a);
    a.ReScale(1.0);

    //DX 9/21/17 - MAGNETIC SYMMETRY - START
    if(vpflow.flag("FULLSYMMETRY::MAGNETIC")){
      string magmom_info = vpflow.getattachedscheme("FULLSYMMETRY::MAGNETIC");
      int num_atoms=a.atoms.size();
      bool is_noncoll=false; 
      vector<xvector<double> > vmag_noncoll;                                          //DX 12/5/17 - added non-collinear
      if(GetNonCollinearMagneticInfo(num_atoms,magmom_info,vmag_noncoll)){            //DX 12/5/17 - added non-collinear
        if(LDEBUG){cerr << "pflow::CalculateFullSymmetry: Non-collinear spin system detected." << endl;} //DX 12/5/17 - added non-collinear
        is_noncoll = true;
        if(!AddSpinToXstructure(a,vmag_noncoll)){
          return false;
        } 
      }                                                                               //DX 12/5/17 - added non-collinear
      bool is_coll=false; 
      vector<double> vmag;
      if(!is_noncoll){
        if(GetCollinearMagneticInfo(num_atoms,magmom_info,vmag)){
          if(LDEBUG){cerr << "pflow::CalculateFullSymmetry: Collinear spin system detected." << endl;}   //DX 12/5/17 - added non-collinear
          is_coll = true;
          if(!AddSpinToXstructure(a,vmag)){
            return false;
          } 
        } 
      }
      if(!is_noncoll && !is_coll){                                                                                 //DX 12/5/17 - added non-collinear
        cerr << "pflow::CalculateFullSymmetry: ERROR: Could not detect collinear or non-collinear spin(s). Check spin input." << endl;//DX 12/5/17 - added non-collinear
        return false;                                                                                              //DX 12/5/17 - added non-collinear
      }
    } 
    //DX 9/21/17 - MAGNETIC SYMMETRY - END

    double default_tolerance=SYM::defaultTolerance(a);
    double tolerance = AUROSTD_NAN;
    if(vpflow.flag("FULLSYMMETRY::TOLERANCE")){
      string tolerance_string = vpflow.getattachedscheme("FULLSYMMETRY::TOLERANCE");
      if(tolerance_string[0] == 't' || tolerance_string[0] == 'T'){ //Tight
        tolerance=default_tolerance;
      }
      else if(tolerance_string[0] == 'l' || tolerance_string[0] == 'L'){ //Loose
        tolerance=default_tolerance*10.0;
      }
      else{
        tolerance=aurostd::string2utype<double>(vpflow.getattachedscheme("FULLSYMMETRY::TOLERANCE"));
      }
    }
    else{
      tolerance = default_tolerance;
    }
    if(tolerance < 1e-10){
      cerr << "pflow::CalculateFullSymmetry: ERROR: Tolerance cannot be zero (i.e. less than 1e-10)" << endl;
      return 0;
    }
    // DX 8/3/17 - Add format flag - START
    string format = "txt";
    if(XHOST.vflag_control.flag("PRINT_MODE::TXT")){
        format = "txt";
      }
    else if(XHOST.vflag_control.flag("PRINT_MODE::JSON")){
        format = "json";
    }
    else{ //default is txt
      format = "txt";
    }
    // DX 8/3/17 - Add format flag - END
    // DX 8/3/17 - Print ouptut to screen - START
    bool print = false;
    if(vpflow.flag("FULLSYMMETRY::SCREEN_ONLY")) {
      print=true;
      defaultKFlags4SymWrite(kflags,false);
      //kflags.KBIN_SYMMETRY_PGROUP_WRITE=FALSE;
      //kflags.KBIN_SYMMETRY_PGROUPK_WRITE=FALSE;
      //kflags.KBIN_SYMMETRY_FGROUP_WRITE=FALSE;
      //kflags.KBIN_SYMMETRY_PGROUP_XTAL_WRITE=FALSE;
      //kflags.KBIN_SYMMETRY_SGROUP_WRITE=FALSE;
      //kflags.KBIN_SYMMETRY_IATOMS_WRITE=FALSE;
      //kflags.KBIN_SYMMETRY_AGROUP_WRITE=FALSE;
      osswrite=FALSE;
    }
    // DX 8/3/17 - Print ouptut to screen - END
    //DX
    //    if(tokens.size()==0) {tolerance=default_tolerance;}
    //    if(tokens.size()>=1 && tokens.at(0) != "--debug" && tokens.at(0) != "--debug" ) {
    //      if(tokens.at(0).at(0) == 't' || tokens.at(0).at(0) == 'T'){ //Tight
    //        tolerance=default_tolerance;
    //      }
    //      else if(tokens.at(0).at(0) == 'l' || tokens.at(0).at(0) == 'L'){ //Loose
    //        tolerance=default_tolerance*10.0;
    //      }
    //      else{ 
    //        tolerance=aurostd::string2utype<double>(tokens.at(0));
    //      }
    //    }
    //DX
    // Perform full scan
    bool no_scan = false;
    if(vpflow.flag("FULLSYMMETRY::NO_SCAN")) {
      no_scan=true;
    }

    bool tobzip2 = TRUE;
    stringstream command;
    ofstream FileMESSAGE("/dev/null");
    //DX 2/21/18 [OBSOLETE] string directory=".";
    string directory=aurostd::execute2string("pwd"); //DX 2/21/18 - use pwd
    if(XHOST.vflag_control.flag("DIRECTORY")) { 
      directory=XHOST.vflag_control.getattachedscheme("DIRECTORY");
    }
    aflags.Directory=directory;
    if(!aurostd::FileExist(directory)) {
      oss << "pflow::CalculateFullSymmetry: ERROR: Unable to locate " << directory << ". Exiting." << endl;
      //return oss.str();
      return FALSE;
    } 

    //while(symmetry_commensurate==FALSE){
    if(print == false){  //DX 8/3/17 - PRINT
    if(aurostd::FileExist(directory+string("/aflow.pgroup.out.bz2"))==TRUE ||
       aurostd::FileExist(directory+string("/aflow.pgroup_xtal.out.bz2"))==TRUE ||
       aurostd::FileExist(directory+string("/aflow.pgroupk.out.bz2"))==TRUE ||
       aurostd::FileExist(directory+string("/aflow.pgroupk_xtal.out.bz2"))==TRUE || //DX 12/5/17 - Added pgroupk_xtal
       aurostd::FileExist(directory+string("/aflow.fgroup.out.bz2"))==TRUE ||
       aurostd::FileExist(directory+string("/aflow.sgroup.out.bz2"))==TRUE ||
       aurostd::FileExist(directory+string("/aflow.agroup.out.bz2"))==TRUE ||
       aurostd::FileExist(directory+string("/aflow.iatoms.out.bz2"))==TRUE) {tobzip2=TRUE;}
    if(aurostd::FileExist(directory+string("/aflow.pgroup.out.bz2"))==TRUE || aurostd::FileExist(directory+string("/aflow.pgroup.out"))==TRUE){   command << "rm -f " << directory << "/aflow.pgroup.*" << endl;}
    if(aurostd::FileExist(directory+string("/aflow.pgroup_xtal.out.bz2"))==TRUE || aurostd::FileExist(directory+string("/aflow.pgroup_xtal.out"))==TRUE){   command << "rm -f " << directory << "/aflow.pgroup_xtal.*" << endl;}
    if(aurostd::FileExist(directory+string("/aflow.pgroupk.out.bz2"))==TRUE || aurostd::FileExist(directory+string("/aflow.pgroupk.out"))==TRUE){ command << "rm -f " << directory << "/aflow.pgroupk.*" << endl;}
    if(aurostd::FileExist(directory+string("/aflow.pgroupk_xtal.out.bz2"))==TRUE || aurostd::FileExist(directory+string("/aflow.pgroupk_xtal.out"))==TRUE){ command << "rm -f " << directory << "/aflow.pgroupk_xtal.*" << endl;} //DX 12/5/17 - Added pgroupk_xtal
    if(aurostd::FileExist(directory+string("/aflow.fgroup.out.bz2"))==TRUE || aurostd::FileExist(directory+string("/aflow.fgroup.out"))==TRUE){ command << "rm -f " << directory << "/aflow.fgroup.*" << endl;}
    if(aurostd::FileExist(directory+string("/aflow.agroup.out.bz2"))==TRUE || aurostd::FileExist(directory+string("/aflow.agroup.out"))==TRUE){ command << "rm -f " << directory << "/aflow.agroup.*" << endl;}
    if(aurostd::FileExist(directory+string("/aflow.sgroup.out.bz2"))==TRUE || aurostd::FileExist(directory+string("/aflow.sgroup.out"))==TRUE){ command << "rm -f " << directory << "/aflow.sgroup.*" << endl;}
    if(aurostd::FileExist(directory+string("/aflow.iatoms.out.bz2"))==TRUE || aurostd::FileExist(directory+string("/aflow.iatoms.out"))==TRUE){ command << "rm -f " << directory << "/aflow.iatoms.*" << endl;}
      if(aurostd::FileExist(directory+string("/aflow.pgroup.json.bz2"))==TRUE ||
	 aurostd::FileExist(directory+string("/aflow.pgroup_xtal.json.bz2"))==TRUE ||
	 aurostd::FileExist(directory+string("/aflow.pgroupk.json.bz2"))==TRUE ||
	 aurostd::FileExist(directory+string("/aflow.pgroupk_xtal.json.bz2"))==TRUE || //DX 12/5/17 - Added pgroupk_xtal
	 aurostd::FileExist(directory+string("/aflow.fgroup.json.bz2"))==TRUE ||
	 aurostd::FileExist(directory+string("/aflow.sgroup.json.bz2"))==TRUE ||
	 aurostd::FileExist(directory+string("/aflow.agroup.json.bz2"))==TRUE ||
	 aurostd::FileExist(directory+string("/aflow.iatoms.json.bz2"))==TRUE) {tobzip2=TRUE;}
      if(aurostd::FileExist(directory+string("/aflow.pgroup.json.bz2"))==TRUE || aurostd::FileExist(directory+string("/aflow.pgroup.json"))==TRUE){   command << "rm -f " << directory << "/aflow.pgroup.*" << endl;}
      if(aurostd::FileExist(directory+string("/aflow.pgroup_xtal.json.bz2"))==TRUE || aurostd::FileExist(directory+string("/aflow.pgroup_xtal.json"))==TRUE){   command << "rm -f " << directory << "/aflow.pgroup_xtal.*" << endl;}
      if(aurostd::FileExist(directory+string("/aflow.pgroupk.json.bz2"))==TRUE || aurostd::FileExist(directory+string("/aflow.pgroupk.json"))==TRUE){ command << "rm -f " << directory << "/aflow.pgroupk.*" << endl;}
      if(aurostd::FileExist(directory+string("/aflow.fgroup.json.bz2"))==TRUE || aurostd::FileExist(directory+string("/aflow.fgroup.json"))==TRUE){ command << "rm -f " << directory << "/aflow.fgroup.*" << endl;}
      if(aurostd::FileExist(directory+string("/aflow.agroup.json.bz2"))==TRUE || aurostd::FileExist(directory+string("/aflow.agroup.json"))==TRUE){ command << "rm -f " << directory << "/aflow.agroup.*" << endl;}
      if(aurostd::FileExist(directory+string("/aflow.sgroup.json.bz2"))==TRUE || aurostd::FileExist(directory+string("/aflow.sgroup.json"))==TRUE){ command << "rm -f " << directory << "/aflow.sgroup.*" << endl;}
      if(aurostd::FileExist(directory+string("/aflow.iatoms.json.bz2"))==TRUE || aurostd::FileExist(directory+string("/aflow.iatoms.json"))==TRUE){ command << "rm -f " << directory << "/aflow.iatoms.*" << endl;}
    if(!command.str().empty()){
      aurostd::execute(command);
    }	
    }	
    
    if(!pflow::PerformFullSymmetry(a,tolerance,no_scan,true,FileMESSAGE,aflags,kflags,osswrite,oss,format)){
      return FALSE;
    }

    // DX 8/3/17 - Print to symmetry operators to screen - START
    if(print==true){
      KBIN_SymmetryToScreen(a,format,oss);
    }
    // DX 8/3/17 - Print to symmetry operators to screen - END

    // BZIP if necessary
    if(tobzip2) {
      command.str(std::string());
      command << "cd " << directory << endl;
      if(format == "txt"){ //DX 8/3/17 - FORMAT
      if(aurostd::FileExist(directory+string("/aflow.pgroup.out"))==TRUE) command << XHOST.command("bzip2") << " -9 aflow.pgroup.out" << endl;
      if(aurostd::FileExist(directory+string("/aflow.pgroup_xtal.out"))==TRUE) command << XHOST.command("bzip2") << " -9 aflow.pgroup_xtal.out" << endl;
      if(aurostd::FileExist(directory+string("/aflow.pgroupk.out"))==TRUE) command << XHOST.command("bzip2") << " -9 aflow.pgroupk.out" << endl;
      if(aurostd::FileExist(directory+string("/aflow.pgroupk_xtal.out"))==TRUE) command << XHOST.command("bzip2") << " -9 aflow.pgroupk_xtal.out" << endl; //DX 12/5/17 - Added pgroupk_xtal
      if(aurostd::FileExist(directory+string("/aflow.fgroup.out"))==TRUE) command << XHOST.command("bzip2") << " -9 aflow.fgroup.out" << endl;
      if(aurostd::FileExist(directory+string("/aflow.sgroup.out"))==TRUE) command << XHOST.command("bzip2") << " -9 aflow.sgroup.out" << endl;
      if(aurostd::FileExist(directory+string("/aflow.agroup.out"))==TRUE) command << XHOST.command("bzip2") << " -9 aflow.agroup.out" << endl;
      if(aurostd::FileExist(directory+string("/aflow.iatoms.out"))==TRUE) command << XHOST.command("bzip2") << " -9 aflow.iatoms.out" << endl;
      }
      else if(format == "json"){ //DX 8/3/17 - FORMAT
	if(aurostd::FileExist(directory+string("/aflow.pgroup.json"))==TRUE) command << XHOST.command("bzip2") << " -9 aflow.pgroup.json" << endl;
	if(aurostd::FileExist(directory+string("/aflow.pgroup_xtal.json"))==TRUE) command << XHOST.command("bzip2") << " -9 aflow.pgroup_xtal.json" << endl;
	if(aurostd::FileExist(directory+string("/aflow.pgroupk.json"))==TRUE) command << XHOST.command("bzip2") << " -9 aflow.pgroupk.json" << endl;
	if(aurostd::FileExist(directory+string("/aflow.pgroupk_xtal.json"))==TRUE) command << XHOST.command("bzip2") << " -9 aflow.pgroupk_xtal.json" << endl; //DX 12/5/17 - Added pgroupk_xtal
	if(aurostd::FileExist(directory+string("/aflow.fgroup.json"))==TRUE) command << XHOST.command("bzip2") << " -9 aflow.fgroup.json" << endl;
	if(aurostd::FileExist(directory+string("/aflow.sgroup.json"))==TRUE) command << XHOST.command("bzip2") << " -9 aflow.sgroup.json" << endl;
	if(aurostd::FileExist(directory+string("/aflow.agroup.json"))==TRUE) command << XHOST.command("bzip2") << " -9 aflow.agroup.json" << endl;
	if(aurostd::FileExist(directory+string("/aflow.iatoms.json"))==TRUE) command << XHOST.command("bzip2") << " -9 aflow.iatoms.json" << endl;
      }
      aurostd::execute(command);
    }
    //return oss.str(); //string("AFLOW "+string(AFLOW_VERSION)+" Symmetry Fixed in "+directory+"  (need aflow>=2948)\n");
    return TRUE; //string("AFLOW "+string(AFLOW_VERSION)+" Symmetry Fixed in "+directory+"  (need aflow>=2948)\n");
  }
} // namespace pflow

// ***************************************************************************
// pflow::fixEmptyAtomNames
// ***************************************************************************
//fixes POSCARs lacking column of atom names
//grabs from POTCAR information of aflow.in
//force_fix: override if a column is present, this may be crucial to fix issues with write_inequivalent_flag
namespace pflow {
  bool fixEmptyAtomNames(xstructure& xstr,bool force_fix) {
    for(uint itype=0;itype<xstr.species.size();itype++) {
      if(xstr.species.size()==xstr.species_pp.size()) {
        if((force_fix || xstr.species.at(itype)=="") && xstr.species_pp.at(itype)!="") 
          xstr.species.at(itype)=KBIN::VASP_PseudoPotential_CleanName(xstr.species_pp.at(itype));
      }
    }  // cormac I`ll write a short pflow for this stuff
    int iatom=0;
    for(uint itype=0;itype<xstr.num_each_type.size();itype++) {
      string species=string(xstr.species.at(itype));
      xstr.species.at(itype)=species;
      for(int j=0;j<xstr.num_each_type.at(itype);j++) {
        xstr.atoms.at(iatom).name=species;    // CONVASP_MODE
        xstr.atoms.at(iatom).CleanName();
        xstr.atoms.at(iatom).CleanSpin();
        xstr.atoms.at(iatom).name_is_given=TRUE;
        iatom++;
      }
    }
    return TRUE;
  }
} // namespace pflow
// DX AND COREY - END

// ***************************************************************************
// pflow::EXTRACT_Symmetry
// ***************************************************************************
namespace pflow {
  string EXTRACT_Symmetry(_aflags &aflags,vector<string> argv) {
    stringstream command;
    string directory=argv.at(2);
    bool tobzip2=FALSE;
    directory=aurostd::RemoveSubStringFirst(directory,_AFLOWIN_);
    aflags.Directory=directory;
    if(aurostd::FileExist(directory+string("/"+_AFLOWIN_))==FALSE) {
      cerr << "pflow::EXTRACT_Symmetry:  File/Directory not found, nothing to do" << endl;
      exit(0);
    }
    command.str(std::string());
    command << "cd " << directory << endl;
    if(aurostd::FileExist(directory+string("/aflow.pgroup.out.bz2"))==TRUE ||
       aurostd::FileExist(directory+string("/aflow.pgroupk.out.bz2"))==TRUE ||
       aurostd::FileExist(directory+string("/aflow.fgroup.out.bz2"))==TRUE ||
       aurostd::FileExist(directory+string("/aflow.sgroup.out.bz2"))==TRUE ||
       aurostd::FileExist(directory+string("/aflow.agroup.out.bz2"))==TRUE ||
       aurostd::FileExist(directory+string("/aflow.iatoms.out.bz2"))==TRUE) tobzip2=TRUE;
    if(aurostd::FileExist(directory+string("/aflow.pgroup.out.bz2"))==TRUE || aurostd::FileExist(directory+string("/aflow.pgroup.out"))==TRUE)   command << "rm -f aflow.pgroup.*" << endl;
    if(aurostd::FileExist(directory+string("/aflow.pgroupk.out.bz2"))==TRUE || aurostd::FileExist(directory+string("/aflow.pgroupk.out"))==TRUE) command << "rm -f aflow.pgroupk.*" << endl;
    if(aurostd::FileExist(directory+string("/aflow.fgroup.out.bz2"))==TRUE || aurostd::FileExist(directory+string("/aflow.fgroup.out"))==TRUE) command << "rm -f aflow.fgroup.*" << endl;
    if(aurostd::FileExist(directory+string("/aflow.agroup.out.bz2"))==TRUE || aurostd::FileExist(directory+string("/aflow.agroup.out"))==TRUE) command << "rm -f aflow.agroup.*" << endl;
    if(aurostd::FileExist(directory+string("/aflow.sgroup.out.bz2"))==TRUE || aurostd::FileExist(directory+string("/aflow.sgroup.out"))==TRUE) command << "rm -f sflow.agroup.*" << endl;
    if(aurostd::FileExist(directory+string("/aflow.iatoms.out.bz2"))==TRUE || aurostd::FileExist(directory+string("/aflow.iatoms.out"))==TRUE) command << "rm -f aflow.iatoms.*" << endl;
    aurostd::execute(command);
    // LOAD _AFLOWIN_
    ofstream FileMESSAGE("/dev/null");
    _kflags kflags;kflags.AFLOW_MODE_VASP=TRUE;
    _vflags vflags;_xvasp xvasp;xvasp.clear();
    ifstream FileAFLOWIN(string(directory+string("/"+_AFLOWIN_)).c_str());
    aurostd::InFileExistCheck("pflow::EXTRACT_Symmetry",string(directory+string("/"+_AFLOWIN_)).c_str(),FileAFLOWIN,cerr);
    string AflowIn;
    AflowIn="";char c; while (FileAFLOWIN.get(c)) AflowIn+=c; // READ _AFLOWIN_ and put into AflowIn
    aflags.QUIET=TRUE;XHOST.QUIET=TRUE;
    vflags=KBIN::VASP_Get_Vflags_from_AflowIN(AflowIn,aflags,kflags);
    KBIN::VASP_Produce_POSCAR(xvasp,AflowIn,FileMESSAGE,aflags,kflags,vflags);
    FileAFLOWIN.close();
    // CREATE PGROUP/FGROUP/SGROUP
    bool PGROUPWRITE=TRUE,PGROUPKWRITE=TRUE,FGROUPWRITE=TRUE,SGROUPWRITE=FALSE,IATOMSWRITE=TRUE,AGROUPWRITE=TRUE;
    bool OSSWRITE=TRUE; // to FileMESSAGE, does not matter as it is /dev/null
    SYM::CalculatePointGroup(FileMESSAGE,xvasp.str,aflags,PGROUPWRITE,OSSWRITE,cout);
    SYM::CalculatePointGroupKlattice(FileMESSAGE,xvasp.str,aflags,PGROUPKWRITE,OSSWRITE,cout);
    SYM::CalculateFactorGroup(FileMESSAGE,xvasp.str,aflags,FGROUPWRITE,OSSWRITE,cout);
    SYM::CalculateSpaceGroup(FileMESSAGE,xvasp.str,aflags,SGROUPWRITE,OSSWRITE,cout);
    SYM::CalculateInequivalentAtoms(FileMESSAGE,xvasp.str,aflags,IATOMSWRITE,OSSWRITE,cout);
    SYM::CalculateSitePointGroup(FileMESSAGE,xvasp.str,aflags,AGROUPWRITE,OSSWRITE,cout);
    // BZIP if necessary
    if(tobzip2) {
      command.str(std::string());
      command << "cd " << directory << endl;
      if(aurostd::FileExist(directory+string("/aflow.pgroup.out"))==TRUE) command << XHOST.command("bzip2") << " -9 aflow.pgroup.out" << endl;
      if(aurostd::FileExist(directory+string("/aflow.pgroupk.out"))==TRUE) command << XHOST.command("bzip2") << " -9 aflow.pgroupk.out" << endl;
      if(aurostd::FileExist(directory+string("/aflow.fgroup.out"))==TRUE) command << XHOST.command("bzip2") << " -9 aflow.fgroup.out" << endl;
      if(aurostd::FileExist(directory+string("/aflow.sgroup.out"))==TRUE) command << XHOST.command("bzip2") << " -9 sflow.agroup.out" << endl;
      if(aurostd::FileExist(directory+string("/aflow.agroup.out"))==TRUE) command << XHOST.command("bzip2") << " -9 aflow.agroup.out" << endl;
      if(aurostd::FileExist(directory+string("/aflow.iatoms.out"))==TRUE) command << XHOST.command("bzip2") << " -9 aflow.iatoms.out" << endl;
      aurostd::execute(command);
    }
    return string("AFLOW "+string(AFLOW_VERSION)+" Symmetry Fixed in "+directory+"  (need aflow>=2948)\n");
  }
} // namespace pflow

// ***************************************************************************
// pflow::FGROUP
// ***************************************************************************
namespace pflow {
  void FGROUP(_aflags &aflags,istream& input) {
    cout << aflow::Banner("BANNER_TINY") << endl;
    aflags.QUIET=TRUE;
    xstructure a(input,IOAFLOW_AUTO);
    bool WRITE=TRUE;
    ofstream FileMESSAGE("/dev/null");
    // DX 8/15/17 - Add in consistency checks SYM::CalculatePointGroup(FileMESSAGE,a,aflags,WRITE,TRUE,cout);
    // DX 8/15/17 - Add in consistency checks SYM::CalculateFactorGroup(FileMESSAGE,a,aflags,WRITE,TRUE,cout);
    _kflags kflags;                                   //DX 8/15/17 - Add in consistency checks
    kflags.KBIN_SYMMETRY_CALCULATE_PGROUP=TRUE;       //DX 8/15/17 - Add in consistency checks
    kflags.KBIN_SYMMETRY_CALCULATE_PGROUPK=FALSE;     //DX 8/15/17 - Add in consistency checks
    kflags.KBIN_SYMMETRY_CALCULATE_FGROUP=TRUE;       //DX 8/15/17 - Add in consistency checks
    kflags.KBIN_SYMMETRY_CALCULATE_PGROUP_XTAL=FALSE; //DX 8/15/17 - Add in consistency checks
    kflags.KBIN_SYMMETRY_CALCULATE_PGROUPK_XTAL=FALSE;//DX 8/15/17 - Add in consistency checks //DX 12/5/17 - Added pgroupk_xtal
    kflags.KBIN_SYMMETRY_CALCULATE_SGROUP=FALSE;      //DX 8/15/17 - Add in consistency checks
    kflags.KBIN_SYMMETRY_CALCULATE_IATOMS=FALSE;      //DX 8/15/17 - Add in consistency checks
    kflags.KBIN_SYMMETRY_CALCULATE_AGROUP=FALSE;      //DX 8/15/17 - Add in consistency checks
    pflow::PerformFullSymmetry(a,FileMESSAGE,aflags,kflags,WRITE,cout); //DX 8/15/17 - Add in consistency checks
  }
} // namespace pflow

// ***************************************************************************
// pflow::FINDSYM
// ***************************************************************************
namespace pflow {
  //DX 9/21/17 [OBSOLETE] void FINDSYM(string options,uint mode,istream& input) {
  void FINDSYM(aurostd::xoption& vpflow,uint mode,istream& input) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "pflow::FINDSYM: BEGIN" << endl;
    string flag_name = "";
    if(mode==0){
      flag_name = "FINDSYM_PRINT";
    }
    if(mode==1){
      flag_name = "FINDSYM_EXEC";
    }
    string options = vpflow.getattachedscheme(flag_name);
    vector<string> tokens;
    aurostd::string2tokens(options,tokens,",");
    if(tokens.size()>=2) {
      init::ErrorOption(cout,options,"pflow::FINDSYM","aflow --findsym[_print][=tolerance_relative] < POSCAR");
      exit(0);
    } 
    // move on

    xstructure a(input,IOAFLOW_AUTO);
    double tolerance=DEFAULT_FINDSYM_TOL;
    if(vpflow.flag("SG::TOLERANCE")){
      string tolerance_string = vpflow.getattachedscheme("SG::TOLERANCE");
      vector<string> tol_tokens;
      aurostd::string2tokens(tolerance_string,tol_tokens,",");
      if(tol_tokens.size()==0) tolerance=DEFAULT_FINDSYM_TOL;
      if(tol_tokens.size()==1) tolerance=aurostd::string2utype<double>(tol_tokens.at(0));
    }
    //DX 9/26/17 [OBSOLETE] if(tol_tokens.size()>=1) tolerance=aurostd::string2utype<double>(tol_tokens.at(0));
    // Read in input file.
    if(mode==0) {cout << a.findsym2print(tolerance) << endl;}
    if(mode==1) {cout << a.findsym2execute(tolerance) << endl;}

    if(LDEBUG) cerr << "pflow::FINDSYM: END" << endl;
  }
} // namespace pflow

// ***************************************************************************
// pflow::FIXBANDS
// ***************************************************************************
namespace pflow {
  bool FIXBANDS(_aflags &aflags,string opts) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    stringstream command;
    stringstream aus;
    vector<string> kpoints_old,eigenval_old,tokens,kpoints_new_string;
    stringstream kpoints_new,eigenval_new;
    string directory=aflags.Directory;
    cout << "pflow::FIXBANDS: aflow --fix_bands=POSCAR,KPOINTS.bands.old,EIGENVAL.bands.old,KPOINTS.bands.new,EIGENVAL.bands.new" << endl;
    double grid=0;
    uint points_bz=0,eigenval_size=0,nbands=0;
    bool foundBZ;
    double pos4=0.0;

    vector<string> vopt;
    aurostd::string2tokens(opts,vopt,",");
    uint vopt_counter=0;
    
    if(vopt.size()<5) return FALSE;
    // FILE_POSCAR_IN ---------------------------------------
    string FILE_POSCAR_IN=directory+"/"+vopt.at(vopt_counter++);
    if(aurostd::FileExist(FILE_POSCAR_IN)==FALSE) {
      cout << "pflow::FIXBANDS: File POSCAR: " << FILE_POSCAR_IN << " not found, nothing to do" << endl;
      return FALSE;
    } else {
      cout << "pflow::FIXBANDS: Load File POSCAR: " << FILE_POSCAR_IN << endl;
    }
    ifstream File_POSCAR_IN(FILE_POSCAR_IN.c_str());
    // [OBSOLETE] xstructure str(File_POSCAR_IN,IOVASP_POSCAR);
    xstructure str(File_POSCAR_IN,IOAFLOW_AUTO);
    File_POSCAR_IN.close();
    // FILE_KPOINTS_BANDS_OLD_IN ---------------------------------------
    string FILE_KPOINTS_BANDS_OLD_IN=directory+"/"+vopt.at(vopt_counter++);
    if(aurostd::FileExist(FILE_KPOINTS_BANDS_OLD_IN)==FALSE) {
      cout << "pflow::FIXBANDS: File KPOINTS_BANDS_OLD: " << FILE_KPOINTS_BANDS_OLD_IN << " not found, nothing to do" << endl;
      return FALSE;
    } else {
      cout << "pflow::FIXBANDS: Load File KPOINTS_BANDS_OLD: " << FILE_KPOINTS_BANDS_OLD_IN << endl;
    }
    aurostd::string2tokens(aurostd::file2string(FILE_KPOINTS_BANDS_OLD_IN),kpoints_old,"\n");
    aus << kpoints_old.at(1);aus >> grid;
    // LATTICE ---------------------------------------
    aurostd::string2tokens(kpoints_old.at(0),tokens," ");
    string LATTICE=tokens.at(0);
    if(LATTICE=="KPOINTS:" && tokens.size()>1) LATTICE=tokens.at(1);
    // FILE_EIGENVAL_BANDS_OLD_IN ---------------------------------------
    string FILE_EIGENVAL_BANDS_OLD_IN=directory+"/"+vopt.at(vopt_counter++);
    if(aurostd::FileExist(FILE_EIGENVAL_BANDS_OLD_IN)==FALSE) {
      cout << "pflow::FIXBANDS: File EIGENVAL_BANDS_OLD: " << FILE_EIGENVAL_BANDS_OLD_IN << " not found, nothing to do" << endl;
      return FALSE;
    } else {
      cout << "pflow::FIXBANDS: Load File EIGENVAL_BANDS_OLD: " << FILE_EIGENVAL_BANDS_OLD_IN << endl;
    }
    string tmp_eigenval=aurostd::file2string(FILE_EIGENVAL_BANDS_OLD_IN);   // to avoid empty lines
    aurostd::StringSubst(tmp_eigenval,"\n"," \n");                       // to avoid empty lines
    aurostd::string2tokens(tmp_eigenval,eigenval_old,"\n");            // to avoid empty lines
    // loaded eigenval now operate
    aurostd::string2tokens(eigenval_old.at(8),tokens," ");eigenval_size=tokens.size()-1;eigenval_size++;
    aurostd::string2tokens(eigenval_old.at(5),tokens," ");points_bz=aurostd::string2utype<uint>(tokens.at(1));
    aurostd::string2tokens(eigenval_old.at(5),tokens," ");nbands=aurostd::string2utype<uint>(tokens.at(2));
    aurostd::string2tokens(eigenval_old.at(7),tokens," ");pos4=aurostd::string2utype<double>(tokens.at(3));

    uint step=7;
    vector<xvector<double> > vbzpoint;
    vector<xmatrix<double> > vbzband;

    for(uint ivbz=0;ivbz<points_bz;ivbz++) {
      xvector<double> bzpoint(3);
      aus.clear();aus.str(eigenval_old.at(step++));
      aus.precision(20);
      aus >> bzpoint[1] >> bzpoint[2] >> bzpoint[3] >> pos4;
      vbzpoint.push_back(bzpoint);
      xmatrix<double> bzband(nbands,eigenval_size);
      for(uint ivbands=0;ivbands<nbands;ivbands++) {
	aus.clear();aus.str(eigenval_old.at(step++));
	for(uint ieig=0;ieig<eigenval_size;ieig++)
	  aus >> bzband[ivbands+1][ieig+1];
      }
      vbzband.push_back(bzband);
      step++;
      //    cerr << step << endl;
    }
    for(uint ivbz=0;ivbz<points_bz;ivbz++) {
      //     cerr << vbzpoint.at(ivbz) << endl;
    }

    ;//   cerr << vbzband.at(0) << endl;
    if(LDEBUG) cerr << "DEBUG eigenval_size=" << eigenval_size << endl;
    if(LDEBUG) cerr << "DEBUG points_bz=" << points_bz << endl;
    if(LDEBUG) cerr << "DEBUG nbands=" << nbands << endl;
    if(LDEBUG) cerr << "DEBUG pos4=" << pos4 << endl;

    cout << "pflow::FIXBANDS: LATTICE=" << LATTICE << endl;
    string tmp_kpoints=LATTICE::KPOINTS_Directions(LATTICE,str.lattice,grid,str.iomode,foundBZ);   // to avoid empty lines
    aurostd::StringSubst(tmp_kpoints,"\n"," \n");                            // to avoid empty lines
    aurostd::string2tokens(tmp_kpoints,kpoints_new_string,"\n");           // to avoid empty lines
    // loaded kpoints now operate
    cout << "pflow::FIXBANDS: PATH= " << kpoints_new_string.at(0) << endl;
    kpoints_new.clear();kpoints_new.str(std::string());
    for(uint i=0;i<=3;i++) kpoints_new << kpoints_new_string.at(i) << endl;
    eigenval_new.clear();eigenval_new.str(std::string());
    for(uint i=0;i<=4;i++) eigenval_new << eigenval_old.at(i) << endl;

    //Patch for EIGENVAL.new (Kesong adds it)
    int kpoints_NEW;
    kpoints_NEW = (kpoints_new_string.size()-3)*grid/3;
    eigenval_new << "   " << int(grid) << "  " << kpoints_NEW << "  " << nbands << endl;
    //eigenval_new << " " << endl; //Comment this, move the blank line before the bands data,
    //Make the total lines of EIGENVAL.aflow exactly same with EIGENVAL.vasp, Kesong

    for(uint ikpz=4;ikpz<kpoints_new_string.size();) {
      aurostd::string2tokens(kpoints_new_string.at(ikpz),tokens);
      if(tokens.size()<2) {
	kpoints_new << kpoints_new_string.at(ikpz++) << endl;
      } else {
	// FROM
	xvector<double> kpoint_from(3);
	string string_from;
	aus.clear();aus.str(kpoints_new_string.at(ikpz));
	string_from=kpoints_new_string.at(ikpz++);
	aus >> kpoint_from[1] >> kpoint_from[2] >> kpoint_from[3];
	// TO
	xvector<double> kpoint_to(3);
	string string_to;
	aus.clear();aus.str(kpoints_new_string.at(ikpz));
	string_to=kpoints_new_string.at(ikpz++);
	aus >> kpoint_to[1] >> kpoint_to[2] >> kpoint_to[3];
	// DELTA
	xvector<double> kpoint_delta(3);
	kpoint_delta=(kpoint_to-kpoint_from)/(grid-1.0);
	double kpoins_delta_prec=1.0*modulus(kpoint_delta);
	for(uint jkpz=0;jkpz<(uint) grid;jkpz++) {
	  // generated and scan
	  uint index=0;
	  double err=1.0e6;
	  for(uint ivbz=0;ivbz<vbzpoint.size();ivbz++) {
	    if(modulus(vbzpoint.at(ivbz)-kpoint_from)<err) {
	      index=ivbz;err=modulus(vbzpoint.at(ivbz)-kpoint_from);
	    }	
	  }
	  if(err>kpoins_delta_prec) {
	    cerr << "************************************************" << endl;
	    cerr << "kpoint_from: not found (" << kpoint_from << ")" << endl;
	    cerr << "closest is vbzpoint.at(index)=" << vbzpoint.at(index) << "   with err=" << err << endl;
	    cerr << "************************************************" << endl;
	    return FALSE;
	    exit(0);
	  }
	  // in index there is the good point

	  eigenval_new << endl;  //Put the blank line before the bands data, Kesong
	  eigenval_new.precision(7);
	  eigenval_new.setf(std::ios::scientific,std::ios::floatfield);
	  eigenval_new << "  " << vbzpoint.at(index)[1] << "  " << vbzpoint.at(index)[2] << "  " << vbzpoint.at(index)[3] << "  " << pos4 << endl;
	  for(uint ivbands=1;ivbands<=nbands;ivbands++) {
	    eigenval_new << "   " << (int) vbzband.at(index)[ivbands][1];
	    eigenval_new.precision(5);
	    eigenval_new.setf(std::ios::fixed,std::ios::floatfield);
	    for(uint ieig=2;ieig<=eigenval_size;ieig++)
	      eigenval_new << "   " << vbzband.at(index)[ivbands][ieig];
	    eigenval_new << endl;
	  }
	  //eigenval_new << endl;  //Remove the additional blank line, Put the blank line before the bands data, Kesong
	  // new one
	  kpoint_from=kpoint_from+kpoint_delta;
	}
	// if OK
	kpoints_new << string_from << endl;
	kpoints_new << string_to << endl;
      }
    }
    kpoints_new << endl;
    // FILE_KPOINTS_BANDS_NEW_OUT ---------------------------------------
    string FILE_KPOINTS_BANDS_NEW_OUT=directory+"/"+vopt.at(vopt_counter++);
    aurostd::stringstream2file(kpoints_new,FILE_KPOINTS_BANDS_NEW_OUT);
    cout << "pflow::FIXBANDS: Save File KPOINTS_BANDS_NEW: " << FILE_KPOINTS_BANDS_NEW_OUT << endl;
    // FILE_EIGENVAL_BANDS_NEW_OUT ---------------------------------------
    string FILE_EIGENVAL_BANDS_NEW_OUT=directory+"/"+vopt.at(vopt_counter++);
    aurostd::stringstream2file(eigenval_new,FILE_EIGENVAL_BANDS_NEW_OUT);
    cout << "pflow::FIXBANDS: Save File EIGENVAL_BANDS_NEW: " << FILE_EIGENVAL_BANDS_NEW_OUT << endl;

    return TRUE;

  }
} // namespace pflow

// ***************************************************************************
// pflow::FRAC
// ***************************************************************************
namespace pflow {
  xstructure FRAC(istream& input) {
    xstructure a(input,IOAFLOW_AUTO);
    xstructure b(a);
    b.SetCoordinates(_COORDS_FRACTIONAL_);
    return b;
  }
} // namespace pflow

// ***************************************************************************
// pflow::FROZSL_ANALYZE
// ***************************************************************************
namespace pflow {
  string FROZSL_ANALYZE(istream& input) {
    ostringstream oss;
    oss.clear();
    oss.setf(std::ios::fixed,std::ios::floatfield);
    uint _precision_=10; //was 16 stefano 10 dane
    oss.precision(_precision_);

    vector<string> vinput,tokens,Mrefs;
    vector<double> Erefs;
    aurostd::stream2vectorstring(input,tokens);
    // take the good ones
    for(uint i=0;i<tokens.size();i++)
      if(aurostd::substring2bool(tokens.at(i),":"))
	vinput.push_back(tokens.at(i));
    // find m0a0
    for(uint i=0;i<vinput.size();i++) {
      if(aurostd::substring2bool(vinput.at(i),"m0a0")) {
	aurostd::string2tokens(vinput.at(i),tokens,"m0a0");
	Mrefs.push_back(tokens.at(0)+"m");
	aurostd::string2tokens(vinput.at(i),tokens,"=");
	Erefs.push_back(aurostd::string2utype<double>(tokens.at(tokens.size()-1)));
      }
    }
    //  for(uint i=0;i<Mrefs.size();i++) cout << Mrefs.at(i) << " " << Erefs.at(i) << endl;
    // now print
    double frozsl_eps=1.0e-6;
    double ev2hartree=1.0/27.211383;
    for(uint i=0;i<vinput.size();i++) {
      for(uint j=0;j<Mrefs.size();j++) {
	if(aurostd::substring2bool(vinput.at(i),Mrefs.at(j))) {
	  //cout << Mrefs.at(j) << " " << Erefs.at(j) << endl;
	  aurostd::string2tokens(vinput.at(i),tokens,"=");
	  if(abs(aurostd::string2utype<double>(tokens.at(tokens.size()-1))-Erefs.at(j))>frozsl_eps)
	    {
	      oss <<  aurostd::PaddedPOST(aurostd::utype2string(1000*ev2hartree*(aurostd::string2utype<double>(tokens.at(tokens.size()-1))-Erefs.at(j))),23," ") << " ! ";
	      aurostd::string2tokens(vinput.at(i),tokens);
	      oss << tokens.at(0) << " - " << Mrefs.at(j) << "0a0" << endl;
	    }
	}
      }
    }
    return oss.str();
  }
} // namespace pflow

namespace pflow {
  string FROZSL_ANALYZE_old(istream& input) {
    // this one was looking for the minimum
    ostringstream oss;
    oss.clear();
    oss.setf(std::ios::fixed,std::ios::floatfield);
    uint _precision_=10; //was 16 stefano 10 dane
    oss.precision(_precision_);

    vector<string> vinput,tokens;
    aurostd::stream2vectorstring(input,tokens);
    // take the good ones
    for(uint i=0;i<tokens.size();i++)
      if(aurostd::substring2bool(tokens.at(i),":"))
	vinput.push_back(tokens.at(i));
    // find minimun
    double minE=1.0e6;
    for(uint i=0;i<vinput.size();i++) {
      aurostd::string2tokens(vinput.at(i),tokens,"=");
      minE=min(minE,aurostd::string2utype<double>(tokens.at(tokens.size()-1)));
    }
    double frozsl_eps=1.0e-6;
    double ev2hartree=1.0/27.211383;
    for(uint i=0;i<vinput.size();i++) {
      aurostd::string2tokens(vinput.at(i),tokens,"=");
      if(abs(aurostd::string2utype<double>(tokens.at(tokens.size()-1))-minE)>frozsl_eps)
	{
	  oss << 1000*ev2hartree*(aurostd::string2utype<double>(tokens.at(tokens.size()-1))-minE) << "   ! ";
	  aurostd::string2tokens(vinput.at(i),tokens);
	  oss << tokens.at(0) << endl;
	}
    }
    return oss.str();
  }
} // namespace pflow

// ***************************************************************************
// pflow::FROZSL_INPUT
// ***************************************************************************
namespace pflow {
  string FROZSL_INPUT(void) {
    _aflags aflags;aflags.Directory="./";
    _kflags kflags;
    ofstream oss;
    string aflowin,MESSAGE="pflow::FROZSL_INPUT ERROR";
    aflowin=string(aflags.Directory+_AFLOWIN_);
    if(!aurostd::FileExist(aflowin)) {cerr << MESSAGE << ": file not found " << aflowin << endl;exit(0);}
  
    string AflowIn;aurostd::file2string(aflowin,AflowIn);
    FROZSL::WGET_INPUT(oss,AflowIn,aflags,kflags);
    // cout << oss << endl;
    return "";
  }
} // namespace pflow

// ***************************************************************************
// pflow::FROZSL_OUTPUT
// ***************************************************************************
namespace pflow {
  string FROZSL_OUTPUT(void) {
    _aflags aflags;aflags.Directory="./";
    _kflags kflags;
    ofstream oss;
    string aflowin,MESSAGE="pflow::FROZSL_OUTPUT ERROR";
    aflowin=string(aflags.Directory+_AFLOWIN_);
    if(!aurostd::FileExist(aflowin)) {cerr << MESSAGE << ": file not found " << aflowin << endl;exit(0);}
    FROZSL::WGET_OUTPUT(oss,aflags,kflags);
    stringstream sss;
    aurostd::file2stringstream(aflags.Directory+"/"+_AFLOW_FROZSL_EIGEN_FILE_,sss);
    //  cerr << sss.str() << endl;
    return "";
  }
} // namespace pflow

// ***************************************************************************
// pflow::FROZSL_VASPSETUP
// ***************************************************************************
namespace pflow {
  string FROZSL_VASPSETUP(vector<string> argv,int mode) {
    string strout="";
    string strfile=_FROZSL_VASPSETUP_FILE_;
    // cerr << "HERE argv.size()=" << argv.size() << endl;
    if(mode==0) strout=init::InitGlobalObject("FROZSL_phvaspsetup_AFLOW");
    if(mode==1) strout=init::InitGlobalObject("FROZSL_phvaspsetup_POSCAR");
    if(argv.size()==2) {
      aurostd::string2file(strout,strfile);
      // cerr << "HERE argv.size()=" << argv.size() << endl;
      aurostd::ChmodFile("755",strfile);
      strout=aurostd::execute2string(strfile);
      aurostd::RemoveFile(strfile);
      return strout;
    }
    return strout;
  }
} // namespace pflow

//DX 9/27/17 - get collinear magnetic info - START
// ***************************************************************************
// pflow::GetCollinearMagneticInfo
// ***************************************************************************
namespace pflow {
  bool GetCollinearMagneticInfo(int& num_atoms, string& magmom_info, vector<double>& vmag){
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(aurostd::substring2bool(magmom_info,"OUTCAR")){
      if(aurostd::FileExist(magmom_info)) {
        xOUTCAR outcar;
        outcar.GetPropertiesFile(magmom_info);
        for(uint i=0;i<outcar.vmag.size();i++){
          vmag.push_back(outcar.vmag[i]);
        }
      }
      else{
        cerr << "pflow::GetCollinearMagneticInfo: ERROR: OUTCAR file does not exist." << endl;
        return false;
      }
    }
    else if(aurostd::substring2bool(magmom_info,"INCAR")){
      if(aurostd::FileExist(magmom_info)) {
        stringstream sss;
        aurostd::efile2stringstream(magmom_info,sss);
        vector<string> vcontent;
        aurostd::string2vectorstring(sss.str(),vcontent);
        bool magmom_found = false;
        for(uint i=0;i<vcontent.size();i++){
          if(vcontent[i].find("MAGMOM=") != std::string::npos){
            if(vcontent[i].find("#") != std::string::npos){
              cerr << "pflow::GetNonCollinearMagneticInfo: ERROR: MAGMOM line in INCAR contains a \"#\" prefix." << endl;
              return false;
            }
            else{
              string magmom_values = aurostd::RemoveSubString(vcontent[i],"MAGMOM=");
              vector<string> mag_tokens;
              aurostd::string2tokens(magmom_values,mag_tokens);
              for(uint m=0;m<mag_tokens.size();m++){
                // INCAR allows multiplication of elements to describe magnetic moment (i.e. 2*2.0)
                if(mag_tokens[m].find("*") != std::string::npos){
                  vector<string> tokens;
                  aurostd::string2tokens(mag_tokens[m],tokens,"*");
                  uint num_magmoms = aurostd::string2utype<int>(tokens[0]);
                  for(uint n=0;n<num_magmoms;n++){
                    vmag.push_back(aurostd::string2utype<double>(tokens[1]));
                  }
                }
                else{
                vmag.push_back(aurostd::string2utype<double>(mag_tokens[m]));
              }
              }
              magmom_found=true;
            }
          }
        }
        if(!magmom_found){
          cerr << "pflow::GetCollinearMagneticInfo: ERROR: MAGMOM tag was not found in the INCAR." << endl;
          return false;
        }
      }
            else {
        cerr << "pflow::GetCollinearMagneticInfo: ERROR: INCAR file does not exist." << endl;
              return false;
      }
    }
    else{
      vector<string> mag_tokens;
      aurostd::string2tokens(magmom_info,mag_tokens,",");
      for(uint m=0;m<mag_tokens.size();m++){
        // Allows multiplication of elements to describe magnetic moment (i.e. 2*2.0)
        if(mag_tokens[m].find("*") != std::string::npos){
          vector<string> tokens;
          aurostd::string2tokens(mag_tokens[m],tokens,"*");
          uint num_magmoms = aurostd::string2utype<int>(tokens[0]);
          for(uint n=0;n<num_magmoms;n++){
            vmag.push_back(aurostd::string2utype<double>(tokens[1]));
          }
        }
        else{
          vmag.push_back(aurostd::string2utype<double>(mag_tokens[m]));
        }
      }
    }
    if((int)vmag.size()!=num_atoms){
      if(LDEBUG){
        cerr << "pflow::GetCollinearMagneticInfo: WARNING: Number of magnetic moments is not equivalent to the number of atoms." << endl;
      }
      return false;
    }
    return true;
  }
}
//DX 9/27/17 - get collinear magnetic info - END

//DX 12/5/17 - get non-collinear magnetic info - START
// ***************************************************************************
// pflow::GetNonCollinearMagneticInfo
// ***************************************************************************
namespace pflow {
  bool GetNonCollinearMagneticInfo(int& num_atoms, string& magmom_info, vector<xvector<double> >& vmag_noncoll){
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(aurostd::substring2bool(magmom_info,"OUTCAR")){
      if(aurostd::FileExist(magmom_info)) {
        xOUTCAR outcar;
        outcar.GetPropertiesFile(magmom_info);
        for(uint i=0;i<outcar.vmag_noncoll.size();i++){
          vmag_noncoll.push_back(outcar.vmag_noncoll[i]);
        }
      }
      else{
        cerr << "pflow::GetNonCollinearMagneticInfo: ERROR: OUTCAR file does not exist." << endl;
        return false;
      }
    }
    else if(aurostd::substring2bool(magmom_info,"INCAR")){
      if(aurostd::FileExist(magmom_info)) {
        stringstream sss;
        aurostd::efile2stringstream(magmom_info,sss);
        vector<string> vcontent;
        aurostd::string2vectorstring(sss.str(),vcontent);
        bool magmom_found = false;
        for(uint i=0;i<vcontent.size();i++){
          if(vcontent[i].find("MAGMOM=") != std::string::npos){
            if(vcontent[i].find("#") != std::string::npos){
              cerr << "pflow::GetNonCollinearMagneticInfo: ERROR: MAGMOM line in INCAR contains a \"#\" prefix." << endl;
              return false;
            }
            else{
              string magmom_values = aurostd::RemoveSubString(vcontent[i],"MAGMOM=");
              vector<string> mag_tokens;
              aurostd::string2tokens(magmom_values,mag_tokens);
              // INCAR allows multiplication of elements to describe magnetic moment (i.e. 2*2.0)
              vector<double> all_magmom_tokens;
              for(uint m=0;m<mag_tokens.size();m++){
                if(mag_tokens[m].find("*") != std::string::npos){
                  vector<string> tokens;
                  aurostd::string2tokens(mag_tokens[m],tokens,"*");
                  uint num_magmoms = aurostd::string2utype<int>(tokens[0]);
                  for(uint n=0;n<num_magmoms;n++){
                    all_magmom_tokens.push_back(aurostd::string2utype<double>(tokens[1]));
                  }
                }
                else{
                  all_magmom_tokens.push_back(aurostd::string2utype<double>(mag_tokens[m]));
                }
              }
              // non-collinear check (should be divisible by 3)
              if((int)all_magmom_tokens.size()!=3*num_atoms){
                if(LDEBUG){
                  cerr << "pflow::GetNonCollinearMagneticInfo: WARNING: From INCAR. Number of magnetic moments not divisible by 3; not non-collinear system." << endl;
                }
                return false;
              }
              xvector<double> mag_xyz; 
              uint index = 1;
              for(uint m=0;m<all_magmom_tokens.size();m++){
                mag_xyz(index)=all_magmom_tokens[m];
                index++;
                if(index==4){
                  vmag_noncoll.push_back(mag_xyz);
                  index=1;
                }
              }
              magmom_found=true;
            }
          }
        }
        if(!magmom_found){
          cerr << "pflow::GetNonCollinearMagneticInfo: ERROR: MAGMOM tag was not found in the INCAR." << endl;
          return false;
        }
      }
      else{
        cerr << "pflow::GetNonCollinearMagneticInfo: ERROR: INCAR file does not exist." << endl;
        return false;
      }
    }
    else{
      vector<string> mag_tokens;
      aurostd::string2tokens(magmom_info,mag_tokens,",");
      // Allows multiplication of elements to describe magnetic moment (i.e. 2*2.0)
      vector<double> all_magmom_tokens;
      for(uint m=0;m<mag_tokens.size();m++){
        if(mag_tokens[m].find("*") != std::string::npos){
          vector<string> tokens;
          aurostd::string2tokens(mag_tokens[m],tokens,"*");
          uint num_magmoms = aurostd::string2utype<int>(tokens[0]);
          for(uint n=0;n<num_magmoms;n++){
            all_magmom_tokens.push_back(aurostd::string2utype<double>(tokens[1]));
          }
        }
        else{
          all_magmom_tokens.push_back(aurostd::string2utype<double>(mag_tokens[m]));
        }
      }
      // non-collinear check (should be divisible by 3)
      if((int)all_magmom_tokens.size()!=3*num_atoms){
        if(LDEBUG){
          cerr << "pflow::GetNonCollinearMagneticInfo: WARNING: From manual input. Number of magnetic moments is not three times the number of atoms; not non-collinear system." << endl;
        }
        return false;
      }
      xvector<double> mag_xyz; 
      uint index = 1;
      for(uint m=0;m<all_magmom_tokens.size();m++){
        mag_xyz(index)=all_magmom_tokens[m];
        index++;
        if(index==4){
          vmag_noncoll.push_back(mag_xyz);
          index=1;
        }
      }
    }
  return true;
  }
}
//DX 9/27/17 - get collinear magnetic info - END

// ***************************************************************************
// pflow::GULP
// ***************************************************************************
namespace pflow {
  void GULP(istream& input) {
    xstructure a(input,IOAFLOW_AUTO);
    pflow::PrintGulp(a,cout);
  }
} // namespace pflow

// ***************************************************************************
// pflow::HKL
// ***************************************************************************
namespace pflow {
  void HKL(string options,_aflags &aflags,istream& input) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "pflow::HKL: BEGIN" << endl;  
    vector<string> tokens;
    aurostd::string2tokens(options,tokens,",");
    if(tokens.size()!=3 && tokens.size()!=4) {
      init::ErrorOption(cout,options,"pflow::HKL","aflow --hkl=h,k,l[,bond] < POSCAR");
      exit(0);
    }
    xstructure a(input,IOAFLOW_AUTO);
    if(LDEBUG) cerr << "pflow::HKL: a=" << a << endl;
    vector<vector<double> > planesreducible,planesirreducible;
    
    if(LDEBUG) cerr << "pflow::HKL: mode=" << tokens.size() << endl;
    xvector<double> iparams(tokens.size(),1);
    if(tokens.size()>0) iparams(1)=aurostd::string2utype<double>(tokens.at(0)); 
    if(tokens.size()>1) iparams(2)=aurostd::string2utype<double>(tokens.at(1)); 
    if(tokens.size()>2) iparams(3)=aurostd::string2utype<double>(tokens.at(2)); 
    if(tokens.size()>3) iparams(4)=aurostd::string2utype<double>(tokens.at(3)); 
    if(LDEBUG) cerr << "pflow::HKL: iparams=" << iparams << endl;
  
    surface::GetSurfaceHKL(a,aflags,iparams,planesreducible,planesirreducible,cout);
    if(LDEBUG) cerr << "pflow::HKL: END" << endl;
  }
 
} // namespace pflow

// ***************************************************************************
// pflow::HKLSearchSearch Trivial/Simple/Complete
// ***************************************************************************
namespace pflow {
  void HKLSearch(string options,_aflags &aflags,istream& input,const string& smode) {
    bool LDEBUG=1;//(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "pflow::HKLSearch: BEGIN" << endl;
    vector<string> tokens;
    aurostd::string2tokens(options,tokens,",");
    
    if(LDEBUG) cerr << "pflow::HKLSearch: smode=" << smode << endl;
    if(LDEBUG) cerr << "pflow::HKLSearch: tokens.size()=" << tokens.size() << endl;
    if(LDEBUG) cerr << "pflow::HKLSearch: options=" << options << endl;
    
    if(smode=="HKL_SEARCH_TRIVIAL" && tokens.size()>3) {
      init::ErrorOption(cout,options,"pflow::HKLSearch: (mode=\"HKL_SEARCH_TRIVIAL\")","aflow --hkl_search[=khlmax[,bond[,step]]] < POSCAR");
      exit(0);
    }
    if(smode=="HKL_SEARCH_SIMPLE" && tokens.size()>4) {
      init::ErrorOption(cout,options,"pflow::HKLSearch: (mode=\"HKL_SEARCH_SIMPLE\")","aflow --hkl_search_simple[=cutoff[,bond[,khlmax[,step]]]] < POSCAR");
      exit(0);
    }
    if(smode=="HKL_SEARCH_COMPLETE" && tokens.size()>4) {
      init::ErrorOption(cout,options,"pflow::HKLSearch: (mode=\"HKL_SEARCH_COMPLETE\")","aflow --hkl_search_complete[=cutoff[,bond[,khlmax[,step]]]] < POSCAR");
      exit(0);
    }
    xstructure a(input,IOAFLOW_AUTO);
    //  if(LDEBUG) cerr << "pflow::HKL: a=" << a << endl;
    vector<vector<double> > planesreducible,planesirreducible;
    vector<vector<uint> > planesirreducible_images;

    if(LDEBUG) cerr << "pflow::HKLSearch: tokens.size()=" << tokens.size() << endl;
    xvector<double> iparams(tokens.size(),(tokens.size()!=0));
    //    cerr << iparams.lrows << endl;
    // cerr << iparams.urows << endl;
    if(LDEBUG) cerr << "pflow::HKLSearch: iparams.rows=" << iparams.rows << endl;
    
    if(tokens.size()>0) iparams(1)=aurostd::string2utype<double>(tokens.at(0)); 
    if(tokens.size()>1) iparams(2)=aurostd::string2utype<double>(tokens.at(1)); 
    if(tokens.size()>2) iparams(3)=aurostd::string2utype<double>(tokens.at(2)); 
    if(tokens.size()>3) iparams(4)=aurostd::string2utype<double>(tokens.at(3)); 
    if(LDEBUG) cerr << "pflow::HKLSearch: iparams=" << iparams << endl;
    //   if(LDEBUG) exit(0);
    
    surface::GetSurfaceHKLSearch(a,aflags,iparams,planesreducible,planesirreducible,planesirreducible_images,cout,smode);
    cout << "REDUCIBLE_SIZE " << planesreducible.size() << endl;
    cout << "IRREDUCIBLE_SIZE " << planesirreducible.size() << endl;
    if(LDEBUG) cerr << "pflow::HKLSearch: END" << endl;
  }
} // namespace pflow

// ***************************************************************************
// Function pflow::HNF
// ***************************************************************************
namespace pflow {
  string HNF(vector<string> argv,istream& input,ostream& oss) {
    oss << AFLOWIN_SEPARATION_LINE << endl;    // --------------------------------
    oss << "[AFLOW] " << aflow::Banner("BANNER_TINY") << endl;
    //  oss << "argv.size()=" << argv.size() << endl;
    xstructure str(input,IOAFLOW_AUTO);
    str.BringInCell();
    //  str.partial_occupation_sum=0;
    // GET SYMMETRY
    str.CalculateSymmetryPointGroup();
    //  str.CalculateSymmetryFactorGroup();
    oss << "[AFLOW] CALCULATION OF HNF xstructures" << endl;
    oss << "[AFLOW] pgroup operations = " << str.pgroup.size() << endl;
    oss.flush();
    // oss << "fgroup operations  = " << str.fgroup.size() << endl;
    // oss << "pgroupk operations = " << str.pgroupk.size() << endl;

    int choice=2;
    if(argv.size()>=3) choice=aurostd::args2utype(argv,"--hnf",(int) 2);
    vector<int> vn;
    if(choice>0) vn.push_back(choice);
    if(choice<0)  for(int i=2;i<=-choice;i++) vn.push_back(i);
    if(choice==0) {oss << "[AFLOW] pflow::HNF n=0" << endl << AFLOWIN_SEPARATION_LINE << endl;}

    vector<xmatrix<double> > sHNF; // contains only the supercells HNF that are uniques
    for(uint in=0;in<vn.size();in++) {
      int n=vn.at(in);
      oss << AFLOWIN_SEPARATION_LINE << endl;    // --------------------------------
      oss << "[AFLOW] n = " << n << endl;
      sHNF=CalculateHNF(str,n);
      // PHYSICAL REVIEW B 77, 224115 2008
      // Algorithm for generating derivative structures
      // Gus L. W. Hart and Rodney W. Forcade
      //xmatrix<double> HNF(3,3);
      //vector<xmatrix<double> > vHNF;
      //xmatrix<double> A(3,3),B(3,3),Bi(3,3),Bj(3,3),R(3,3),H(3,3);A=trasp(str.lattice);
      //vector<xmatrix<double> > vB;   // contains lattices of the potential xstructures
      //vector<xmatrix<double> > sHNF; // contains only the supercells HNF that are uniques
      //long long int i=0;
      //bool found=FALSE;
      //for(int a=1;a<=n;a++)
      //  for(int c=1;c<=n/a;c++)
      //    for(int f=1;f<=n/a/c;f++)
      //      if(a*c*f==n)
      //        for(int b=0;b<c;b++)
      //          for(int d=0;d<f;d++)
      //            for(int e=0;e<f;e++) {
      //              i++;
      //              HNF(1,1)=a;HNF(1,2)=0;HNF(1,3)=0;
      //              HNF(2,1)=b;HNF(2,2)=c;HNF(2,3)=0;
      //              HNF(3,1)=d;HNF(3,2)=e;HNF(3,3)=f;
      //              vHNF.push_back(HNF);
      //            }
      //oss << "[AFLOW] HNF = " << vHNF.size() << endl;
      //oss.flush();
      //for(i=0;i<(int) vHNF.size();i++) {
      //  Bi=A*vHNF.at(i);
      //  found=FALSE;
      //  for(uint istr=0;istr<vB.size()&&found==FALSE;istr++)  {                     // cycle through the unique vBs
      //    Bj=vB.at(istr);
      //    for(uint pgroup=0;pgroup<str.pgroup.size()&&found==FALSE;pgroup++) {      // cycle through the pgroup of str
      //      R=trasp(str.pgroup.at(pgroup).Uc);
      //      H=aurostd::inverse(Bj)*aurostd::inverse(R)*Bi;
      //      if(aurostd::isinteger(H)) found=TRUE;
      //    }
      //  }
      //  if(found==FALSE) { // not found, then plug
      //    vB.push_back(Bi);
      //    sHNF.push_back(vHNF.at(i));
      //    //      cerr << "+";
      //  }
      //}
      //  oss << endl;
      //oss << "[AFLOW] supercells = " << vB.size() << endl;
      oss << "[AFLOW] supercells = " << sHNF.size() << endl;
      //  exit(0);
      oss << AFLOWIN_SEPARATION_LINE << endl;  // --------------------------------
      oss.flush();
      for(uint j=0;j<sHNF.size();j++) {
	xstructure strj=GetSuperCell(str,trasp(sHNF.at(j)));
	// strj.Standard_Conventional_UnitCellForm();
	// strj.iomode=IOVASP_ABCCAR;
	strj.title+="  [HNF("+aurostd::utype2string(n)+")="+aurostd::utype2string(j+1)+"/"+aurostd::utype2string(sHNF.size())+"=";
	for(uint ii=1;ii<=3;ii++) for(uint jj=1;jj<=3;jj++) strj.title+=aurostd::utype2string(sHNF.at(j)(ii,jj))+(ii*jj<9?",":"]");
	oss << "[VASP_POSCAR_MODE_EXPLICIT]START.HNF_" << n << "_" << j+1 << "_" << sHNF.size() << endl;
	oss << strj;// << endl;
	// oss << strj.title << endl;
	oss << "[VASP_POSCAR_MODE_EXPLICIT]STOP.HNF_" << n << "_" << j+1 << "_" << sHNF.size() << endl;
	oss << AFLOWIN_SEPARATION_LINE << endl;  // --------------------------------
	oss.flush();
      }
      // if(0) {    cerr << vB.at(1) << endl << endl;    swap_rows(vB.at(1),1,3);    cerr << vB.at(1) << endl << endl; }
      //  for(uint istr=0;istr<vB.size();istr++)  {   cerr << vB.at(istr) << endl << endl; }
    }
    cerr << aurostd::ostream2string(oss) << endl;
    return aurostd::ostream2string(oss);
  }
} // namespace pflow

// ***************************************************************************
// Function pflow::HNFTOL
// ***************************************************************************
namespace pflow {
  string HNFTOL(vector<string> argv,istream& input,ostream& oss) {
    oss << AFLOWIN_SEPARATION_LINE << endl;    // --------------------------------
    oss << "[AFLOW] " << aflow::Banner("BANNER_TINY") << endl;
    //  oss << "argv.size()=" << argv.size() << endl;
    xstructure str(input,IOAFLOW_AUTO);
    str.BringInCell();
    //  str.partial_occupation_sum=0;
    // GET SYMMETRY
    //str.CalculateSymmetryPointGroup(); //CO, cannot calculate symmetry of a structure with overlapping atoms
    //  str.CalculateSymmetryFactorGroup();
    oss << "[AFLOW] CALCULATION OF HNFTOL xstructures" << endl;
    //oss << "[AFLOW] pgroup operations = " << str.pgroup.size() << endl;
    oss.flush();
    // oss << "fgroup operations  = " << str.fgroup.size() << endl;
    // oss << "pgroupk operations = " << str.pgroupk.size() << endl;
    uint digits=6,digits1=4,digits2=2*digits+11;

    double tolerance=DEFAULT_PARTIAL_OCCUPATION_TOLERANCE;
    if(argv.size()>=3) tolerance=aurostd::args2utype(argv,"--hnftol",(double)  DEFAULT_PARTIAL_OCCUPATION_TOLERANCE);
    if(argv.size()==2) {tolerance=str.partial_occupation_tol;}
    oss << "[AFLOW] hnf_tolerance=" << tolerance << endl;
    oss << AFLOWIN_SEPARATION_LINE << endl;    // --------------------------------
    oss << "[VASP_POSCAR_MODE_EXPLICIT]START" << endl;
    oss << str;// << endl;
    oss << "[VASP_POSCAR_MODE_EXPLICIT]STOP" << endl;
    oss << AFLOWIN_SEPARATION_LINE << endl;    // --------------------------------
    oss.precision(digits);
    double error=1.0,eps;
    int HNFi=0;
    vector<double> error_iatom;  vector<double> effective_pocc_iatom;  vector<uint> M_iatom;
    for(uint iatom=0;iatom<str.atoms.size();iatom++) error_iatom.push_back(1.0);
    for(uint iatom=0;iatom<str.atoms.size();iatom++) effective_pocc_iatom.push_back(0.0);
    for(uint iatom=0;iatom<str.atoms.size();iatom++) M_iatom.push_back(0);
  
    oss << aurostd::PaddedPRE(aurostd::utype2string(0),digits1) << "  " << "--------" << " ";
    for(uint iatom=0;iatom<str.atoms.size();iatom++)
      if(str.atoms.at(iatom).partial_occupation_value<1.0)
	oss << "| "+aurostd::PaddedPOST("iatom="+aurostd::utype2string(iatom)+"/"+aurostd::utype2string(str.atoms.size()),digits2) << " " ;
    oss << " | " << "error" << endl;
    while(error>tolerance) {
      HNFi++;
      error=1.0/HNFi;
      oss << aurostd::PaddedPRE(aurostd::utype2string(HNFi),digits1) << "  " << error << " ";//endl;
      for(uint iatom=0;iatom<str.atoms.size();iatom++) {
	error_iatom.at(iatom)=1.0;
	M_iatom.at(iatom)=0;
	for(int j=0;j<=HNFi;j++) {
	  eps=abs((double) j/HNFi-str.atoms.at(iatom).partial_occupation_value);
	  if(eps<error_iatom.at(iatom)) {
	    error_iatom.at(iatom)=eps;
	    M_iatom.at(iatom)=(uint) j;
	    effective_pocc_iatom.at(iatom)=(double) M_iatom.at(iatom)/HNFi;
	  }
	}
      }
      for(uint iatom=0;iatom<str.atoms.size();iatom++)
	if(str.atoms.at(iatom).partial_occupation_value<1.0) {
	  stringstream aus;aus.precision(digits);aus.setf(std::ios::fixed,std::ios::floatfield);
	  aus << effective_pocc_iatom.at(iatom) << "," << error_iatom.at(iatom) << "," << M_iatom.at(iatom) << "/" << HNFi;
	  oss << "| " << aurostd::PaddedPOST(aus.str(),digits2) << " " ;
	}
      error=aurostd::max(error_iatom);
      oss << " | " << error << endl;
    }
    oss << AFLOWIN_SEPARATION_LINE << endl;    // --------------------------------
    // DONE
    return aurostd::ostream2string(oss);
  }
} // namespace pflow

// ***************************************************************************
// pflow::IDENTICAL
// ***************************************************************************
namespace pflow {
  xstructure IDENTICAL(istream& input) {
    xstructure a(input,IOAFLOW_AUTO);
    xstructure b(a);
    b=IdenticalAtoms(b);
    return b;
  }
} // namespace pflow

// ***************************************************************************
// pflow::INCELL
// ***************************************************************************
namespace pflow {
  xstructure INCELL(istream& input) {
    xstructure a(input,IOAFLOW_AUTO);
    xstructure b(a);
    b=BringInCell(b);
    return b;
  }
} // namespace pflow

// ***************************************************************************
// pflow::INCOMPACR
// ***************************************************************************
namespace pflow {
  xstructure INCOMPACT(istream& input) {
    xstructure a(input,IOAFLOW_AUTO);
    xstructure b;
    b=PutInCompact(a);
    return b;
  }
} // namespace pflow

// ***************************************************************************
// pflow::INTPOL
// ***************************************************************************
namespace pflow {
  void INTPOL(string options) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "pflow::INTPOL: BEGIN" << endl;
    vector<string> tokens;
    aurostd::string2tokens(options,tokens,",");
    if(tokens.size()!=4) {
      init::ErrorOption(cout,options,"pflow::INTPOL","aflow --intpol=file1,file2,nimages,nearest_image_flag");
      exit(0);
    }
    if(!aurostd::FileExist(tokens.at(0))) {
      cerr << "ERROR - pflow::INTPOL: file1=" << tokens.at(0) << " does not exist..." << endl;
      exit(0);
    }
    if(!aurostd::FileExist(tokens.at(1))) {
      cerr << "ERROR - pflow::INTPOL: file2=" << tokens.at(1) << " does not exist..." << endl;
      exit(0);
    }
    xstructure strA(tokens.at(0),IOAFLOW_AUTO);
    xstructure strB(tokens.at(1),IOAFLOW_AUTO);
    
    cout << aflow::Banner("BANNER_TINY") << endl;
    int nimages=aurostd::string2utype<int>(tokens.at(2));
    string nearest_image_flag=tokens.at(3);
    PrintImages(strA,strB,nimages,nearest_image_flag);
    if(LDEBUG) cerr << "pflow::INTPOL: END" << endl;
  }
} // namespace pflow

// ***************************************************************************
// pflow::INWS
// ***************************************************************************
namespace pflow {
  xstructure INWS(istream& input) {
    xstructure a(input,IOAFLOW_AUTO);
    xstructure b;
    b=BringInWignerSeitz(a);
    return b;
  }
} // namespace pflow

// ***************************************************************************
// pflow::JMOL
// ***************************************************************************
// Stefano Curtarolo (Dec-2009) //modification Richard Taylor // modification SC Oct2014
namespace pflow {
  void JMOL(string options,istream& input) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "pflow::JMOL: BEGIN" << endl;
    string COLOR="white";          // default
    bool WRITE=false;              // default
    xvector<int> ijk(3);           // default 
    ijk[1]=1;ijk[2]=1;ijk[3]=1;    // default
    
    vector<string> tokens;
    aurostd::string2tokens(options,tokens,",");
    
    if(tokens.size()>5) {
      init::ErrorOption(cout,options,"pflow::JMOL","aflow --jmol[=n1[,n2[,n3[,color[,true/false]]]]] < POSCAR");
      exit(0);
    }
    
    if(tokens.size()>=1) ijk[1]=aurostd::string2utype<int>(tokens.at(0)); 
    if(tokens.size()>=2) ijk[2]=aurostd::string2utype<int>(tokens.at(1)); 
    if(tokens.size()>=3) ijk[3]=aurostd::string2utype<int>(tokens.at(2)); 
    if(tokens.size()>=4) COLOR=tokens.at(3); 
    if(tokens.size()>=5) {
      if(tokens.at(4)=="true" || tokens.at(4)=="1") {
	WRITE=true;
	cerr << "pflow::JMOL: saving output..." << endl;
      }
    }
    //ofstream script; //file to determine color of background, bond width etc
    ofstream ras;
    string rasFILE="ras."+strPID();
    ras.open(rasFILE.c_str(),std::ios::out);
    ras << "background " << COLOR << endl
	<< "wireframe " << 0.05 << endl; //radius of bond in angstroms
    //anymore? (http://jmol.sourceforge.net/docs/JmolUserGuide/ch04.html)
    ras.close();
    
    ofstream FileOUTPUT;
    string FileOUTPUTName="aflow.jmol.xyz."+strPID();
    FileOUTPUT.open(FileOUTPUTName.c_str(),std::ios::out);
    if(FileOUTPUT) {
      xvector<int> _ijk(vabs(ijk));
      if(max(_ijk)==0) _ijk.set(1);
      ostringstream aus_exec;
      xstructure a(input,IOAFLOW_AUTO);
      if(_ijk(1)!=1 || _ijk(2)!=1 || _ijk(3)!=1) a=GetSuperCell(a,_ijk);
      if(0) {
	xstructure b;
	// a.BringInCell();b=a;
	b=a;
	_atom atom;
	for(uint iat=0;iat<a.atoms.size();iat++) {
	  for(double i=0;i<=1;i+=0.99)
	    for(double j=0;j<=1;j+=0.99)
	      for(double k=0;k<=1;k+=0.99) {
		atom=a.atoms.at(iat);
		atom.fpos[1]+=i;atom.fpos[2]+=j;atom.fpos[3]+=k;
		if(aurostd::isequal(atom.fpos[1],1.0,0.02) || aurostd::isequal(atom.fpos[2],1.0,0.02)|| aurostd::isequal(atom.fpos[3],1.0,0.02))
		  b.AddAtom(F2C(b,atom));
	      }
	}
	a=b;
      }
      pflow::PrintCIF(FileOUTPUT,a);
      FileOUTPUT.close();
      // cerr << _ijk << endl;
      aus_exec << XHOST.command("jmol") << " -s " << rasFILE << " " << FileOUTPUTName << endl;// " && rm -f " << FileOUTPUTName << endl;
      aurostd::execute(aus_exec);
      // aurostd::StringstreamClean(aus_exec);
      if(WRITE == false) {
	aus_exec << "rm -f " <<  FileOUTPUTName << endl;
	aurostd::execute(aus_exec);
      }
      aus_exec << "rm -f " <<  rasFILE << endl;
      aurostd::execute(aus_exec);
    } else FileOUTPUT.close();
  }
} // namespace pflow

// ***************************************************************************
// pflow::KBAND
// ***************************************************************************
// Wahyu Setyawan (Nov-2009)
namespace pflow {
  void KBAND(vector<string> argv) {
    // 2008 wahyu setyawan    ws26@duke.edu
    // 2008 roman chepulskyy  rc74@duke.edu
    /*
      cout<<" uband -h     : help\n"
      <<" uband 1      : nonspin\n"
      <<" uband 2      : spin-polarized.\n";
    */
#define _NSEGMAX_ 50

    int argc=argv.size();
    ifstream kin,ein,din;

    if(argc!=3) return;
    if(!(atoi(&argv.at(2)[0])==1 or atoi(&argv.at(2)[0])==2)) return;

    int i,ind,ispin,iseg,itmp,j,ngrid,nseg,nband;
    float ftmp,kseg[_NSEGMAX_][4];
    float kx1,ky1,kz1,kx2,ky2,kz2,dkx,dky,dkz,dk,koffset;
    string tmpstr;

    ispin=atoi(&argv.at(2)[0]);
    kin.open("KPOINTS");
    getline(kin,tmpstr);
    kin>>ngrid;    getline(kin,tmpstr);
    getline(kin,tmpstr);
    getline(kin,tmpstr);
    i=0;
    while(!kin.eof()) {
      i++;
      kin>>kseg[i][1]>>kseg[i][2]>>kseg[i][3];   getline(kin,tmpstr);
    }
    nseg=i/2;
    kin.close();
    //constructing kline from the norm of kpoints cascaded for continuous plot
    xvector<float> kline(nseg*ngrid);
    koffset=0; ind=1;
    for(iseg=1;iseg<=nseg;iseg++) {
      i = 2*iseg-1;
      kx1=kseg[i][1];    ky1=kseg[i][2];    kz1=kseg[i][3];
      kx2=kseg[i+1][1];  ky2=kseg[i+1][2];  kz2=kseg[i+1][3];
      //    knorm1=sqrt(kx1*kx1 + ky1*ky1 + kz1*kz1); //kstart
      //knorm2=koffset+sqrt(kx2*kx2 + ky2*ky2 + kz2*kz2); //kstop
      dkx = kx2-kx1;
      dky = ky2-ky1;
      dkz = kz2-kz1;
      dk = (sqrt(dkx*dkx + dky*dky + dkz*dkz))/(ngrid-1);
      kline[ind++]=koffset;
      for(j=2;j<=ngrid;j++) {
	kline[ind++] = koffset+(j-1)*dk;
      }
      koffset = koffset+dk*(ngrid-1);
    }
    //creating overall kband data
    //format:
    //             spinup                 spindown
    //kline[1] band1 band2 band3 ...  band1 band2 band3 ...
    //kline[2] band1 band2 band3 ...  band1 band2 band3 ...
    //...
    //this format is nice for plotting with kpoints as the x-axis
    ein.open("EIGENVAL");
    getline(ein,tmpstr);
    getline(ein,tmpstr);
    getline(ein,tmpstr);
    getline(ein,tmpstr);
    getline(ein,tmpstr);
    ein>>itmp>>itmp;
    if(nseg!=itmp/ngrid) {cerr<<"error nseg inconsistent between KPOINTS and EIGENVAL\n"; return;}
    ein>>nband;    getline(ein,tmpstr);
    xmatrix<float> kband(nseg*ngrid,ispin*nband+1);
    for(i=1;i<=nseg*ngrid;i++) {
      ein>>ftmp>>ftmp>>ftmp>>ftmp;
      kband[i][1]=kline[i];
      for(j=1;j<=nband;j++) {
	ein>>ftmp>>kband[i][j+1];
	if(ispin==2) ein>>kband[i][nband+j+1];
      }
    }
    ein.close();
    //output
    for(i=1;i<=nseg*ngrid;i++) {
      for(j=1;j<=nband*ispin+1;j++) {
	cout<<"  "<<kband[i][j];
      }
      cout<<endl;
    }
  }
} // namespace pflow

// ***************************************************************************
// pflow::KPATH
// ***************************************************************************
namespace pflow {
  void KPATH(istream& input,double grid,bool WWW) {
    xstructure str_in(input,IOAFLOW_AUTO);
    xstructure str_sp,str_sc;
    //DX 8/29/17 [OBSOLETE] LATTICE::Standard_Lattice_StructureDefault(str_in,str_sp,str_sc);
    bool full_sym=false; //DX 8/29/17 - Speed increase
    LATTICE::Standard_Lattice_Structure(str_in,str_sp,str_sc,full_sym); //DX 8/29/17 - Speed increase
    string lattice_type;
    stringstream oss;string sss;
    bool foundBZ;
    //    if(grid<=0.0) grid=16.0; // no more default
    lattice_type=str_sp.bravais_lattice_variation_type;
    //  if(WWW) cout << "<pre>" << endl;
    cout << "// STRUCTURE TO RUN ***********************************************************************"<< endl;
    oss.clear();oss << str_sp; sss=oss.str(); 
    if(WWW) {aurostd::StringSubst(sss,"<","&#60;");aurostd::StringSubst(sss,">","&#62;");}
    cout << sss;
    cout << "// KPOINTS TO RUN *************************************************************************"<< endl;
    sss=LATTICE::KPOINTS_Directions(lattice_type,str_sp.lattice,grid,str_sp.iomode,foundBZ); 
    if(WWW) {aurostd::StringSubst(sss,"<","&#60;");aurostd::StringSubst(sss,">","&#62;");}
    cout << sss;
    // cout << LATTICE::KPOINTS_Directions(lattice_type,str_sp.lattice,grid,foundBZ);
    //  if(WWW) cout << "</pre>" << endl;
    if(WWW) {
      cout << "// PATH ***********************************************************************************"<< endl;
      //   cout << lattice_type << endl;
      //   cout << str_sp.lattice << endl;
      cout << "</pre>" << endl;
      cout << "<img height=600 src=http://" << XHOST.AFLOW_MATERIALS_SERVER << "/SCIENCE/images/brillouin/" << lattice_type << ".PNG><br>" << endl;
      cout << "<br> [ ";
      cout << "<a href=http://" << XHOST.AFLOW_MATERIALS_SERVER << "/SCIENCE/images/brillouin/" << lattice_type << ".PNG>png</a>";
      cout << " | ";
      cout << "<a href=http://" << XHOST.AFLOW_MATERIALS_SERVER << "/SCIENCE/images/brillouin/" << lattice_type << ".EPS>eps</a>";
      cout << " ] ";
      cout << "<pre>" << endl;
    }
    cout << "// END ************************************************************************************"<< endl;
  }
} // namespace pflow

// ***************************************************************************
// pflow::KPOINTS
// ***************************************************************************
namespace pflow {
  xstructure KPOINTS(string options, istream& input, ostream& oss) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "pflow::KPOINTS: BEGIN" << endl;
    vector<string> tokens;
    aurostd::string2tokens(options,tokens,",");
    if(tokens.size()!=1 ) {
      init::ErrorOption(cout,options,"pflow::KPOINTS","aflow --kpoints=KDENS [or --kppra=KDENS,-k=KDENS] < POSCAR");
      exit(0);
    }
    int NK=aurostd::string2utype<int>(tokens.at(0));

    oss << aflow::Banner("BANNER_TINY") << endl;
    xstructure str(input,IOAFLOW_AUTO);
    if(NK==0) {
      cerr << "ERROR - pflow::KPOINTS: --kpoints=NUMBER: NUMBER must be bigger than zero" << endl;exit(0);
    } else {
      ofstream FileDUMMY;
      if(NK>10000000) NK=10000000;
      oss << kintro << "KPPRA     = " << NK << " (requested) " << endl;
      double NK_tmp=(int) ((double) NK/str.atoms.size()+0.5);if(NK<1) NK=1;  //CO 180226 - this is done inside KPPRA(), so don't overwrite NK
      oss << kintro << "KPPRA/#AT = "<<NK_tmp<<endl;  //CO 180226
      KPPRA(str,NK);
      //	oss << kintro << "KPOINTS   = [" << str.kpoints_k1 << "," << str.kpoints_k2 << "," << str.kpoints_k3 << "]" << endl;
      oss << kintro << "DKGRID    = ["
	  << modulus(str.klattice(1))/str.kpoints_k1 << ","
	  << modulus(str.klattice(2))/str.kpoints_k2 << ","
	  << modulus(str.klattice(3))/str.kpoints_k3 << "]" << endl;
      oss << kintro << "KPPRA     = " << str.kpoints_kppra << " (found) " << endl;
      // oss << kintro << "next line for automatic scripting (with cat POSCAR | aflow --kpoints | grep -i AUTO | sed \"s/AUTO//g\")" << endl;
      oss << kintro << "KPOINTS   = " << str.kpoints_k1 << " " << str.kpoints_k2 << " " << str.kpoints_k3 <<endl;
      if(LDEBUG) cerr << "pflow::KPOINTS: END" << endl;
      return str;
    }
  }
} // namespace pflow

// ***************************************************************************
// pflow::KPOINTS_DELTA
// ***************************************************************************
namespace pflow {
  xstructure KPOINTS_DELTA(aurostd::xoption& vpflow, istream& input, ostream& oss) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    //[OBSOLETE CO 171010]if(LDEBUG) cerr << "pflow::KPOINTS_DELTA: BEGIN" << endl;
    //[OBSOLETE CO 171010]vector<string> tokens;
    //[OBSOLETE CO 171010]aurostd::string2tokens(options,tokens,",");
    //[OBSOLETE CO 171010]if(tokens.size()!=1 ) {
    //[OBSOLETE CO 171010]  init::ErrorOption(cout,options,"pflow::KPOINTS_DELTA","aflow --delta_kpoints=number [or --dkpoints=number | -dkpoints=number | -dk=number] < POSCAR");
    //[OBSOLETE CO 171010]  exit(0);
    //[OBSOLETE CO 171010]}
    double DK=aurostd::string2utype<double>(vpflow.getattachedscheme("FLAG::XVASP_KPOINTS_DELTA")); //tokens.at(0));  //CO 171010
    
    oss << aflow::Banner("BANNER_TINY") << endl;
    xstructure str(input,IOAFLOW_AUTO);
    if(DK<=1.0E-6) {
      cerr << "ERROR - pflow::KPOINTS_DELTA: --delta_kpoints=NUMBER: NUMBER must be bigger than zero (dk=" << DK << ")" << endl;exit(0);
    } else {
      ofstream FileDUMMY;
      oss << kintro << "DK     = " << DK << " (requested) " << endl;
      KPPRA_DELTA(str,DK);
      //	oss << kintro << "KPOINTS   = [" << str.kpoints_k1 << "," << str.kpoints_k2 << "," << str.kpoints_k3 << "]" << endl;
      oss << kintro << "DKGRID    = ["
	  << modulus(str.klattice(1))/str.kpoints_k1 << ","
	  << modulus(str.klattice(2))/str.kpoints_k2 << ","
	  << modulus(str.klattice(3))/str.kpoints_k3 << "]" << endl;
      oss << kintro << "KPPRA     = " << str.kpoints_kppra << " (found) " << endl;
      // oss << kintro << "next line for automatic scripting (with cat POSCAR | aflow --kpoints | grep -i AUTO | sed \"s/AUTO//g\")" << endl;
      oss << kintro << "KPOINTS   = " << str.kpoints_k1 << " " << str.kpoints_k2 << " " << str.kpoints_k3 <<endl;
      if(LDEBUG) cerr << "pflow::KPOINTS_DELTA: END" << endl;
      return str;
    }
  }
} // namespace pflow

// ***************************************************************************
// pflow::JOINSTRLIST
// ***************************************************************************
namespace pflow {
  void JOINSTRLIST(vector<string> argv) {
    ifstream list1_inf(argv.at(2).c_str());
    aurostd::InFileExistCheck("aflow",argv.at(2).c_str(),list1_inf,cerr);
    ifstream list2_inf(argv.at(3).c_str());
    aurostd::InFileExistCheck("aflow",argv.at(3).c_str(),list2_inf,cerr);
    std::vector<xstructure> str_vec_1(0),str_vec_2(0),str_vec_tot(0);
    pflow::ReadInStrVec(str_vec_1,list1_inf);
    pflow::ReadInStrVec(str_vec_2,list2_inf);
    pflow::JoinStrVec(str_vec_1,str_vec_2,str_vec_tot);
    pflow::PrintStrVec(str_vec_tot,cout);
  }
} // namespace pflow

// ***************************************************************************
// pflow::LATTICEREDUCTION
// ***************************************************************************
namespace pflow {
  xstructure LATTICEREDUCTION(istream& input) {
    xstructure a(input,IOAFLOW_AUTO);
    return LatticeReduction(a);
  }
} // namespace pflow

// ***************************************************************************
// pflow::LATTICE_TYPE
// ***************************************************************************
namespace pflow {
  string LATTICE_TYPE(istream& input) {
    stringstream sss;
    xstructure a(input,IOAFLOW_AUTO);
    //DX 8/24/17 [OBSOLETE] a.GetLatticeType();
    xstructure str_sp,str_sc; //DX 8/24/17 - Speed increase
    bool full_sym=false; //DX 8/29/17 - Speed increase
    LATTICE::Standard_Lattice_Structure(a,str_sp,str_sc,full_sym); //DX 8/29/17 - Speed increase
    // sss << " Real space lattice primitive           = " << a.bravais_lattice_type << endl;
    // sss << " Real space lattice variation           = " << a.bravais_lattice_variation_type << endl;//wahyu mod
    // sss << " Real space conventional lattice        = " << a.bravais_conventional_lattice_type << endl;
    // sss << " Real space Pearson symbol              = " << a.pearson_symbol << endl;
    sss << str_sp.bravais_lattice_type << "," << str_sp.bravais_lattice_variation_type << endl; //wahyu mod //DX 8/24/17 - a to str_sp
    //  sss << a.bravais_lattice_type << "," << a.bravais_conventional_lattice_type << endl;
    return sss.str();
  }
} // namespace pflow

// ***************************************************************************
// pflow::LATTICE_LATTICE_TYPE
// ***************************************************************************
namespace pflow {
  string LATTICE_LATTICE_TYPE(istream& input) {
    stringstream sss;
    xstructure a(input,IOAFLOW_AUTO);
    //DX 8/24/17 [OBSOLETE] a.GetLatticeType();
    xstructure str_sp,str_sc; //DX 8/24/17 - Speed increase
    bool full_sym=false; //DX 8/29/17 - Speed increase
    LATTICE::Bravais_Lattice_StructureDefault(a,str_sp,str_sc,full_sym); //DX 8/29/17 - Speed increase
    // sss << " Real space lattice primitive           = " << a.bravais_lattice_lattice_type << endl;
    // sss << " Real space lattice variation           = " << a.bravais_lattice_lattice_variation_type << endl;//wahyu mod
    // sss << " Real space conventional lattice        = " << a.bravais_conventional_lattice_lattice_type << endl;
    // sss << " Real space Pearson symbol              = " << a.pearson_symbol << endl;
    sss << str_sp.bravais_lattice_lattice_type << "," << str_sp.bravais_lattice_lattice_variation_type << endl; //wahyu mod /DX 8/24/17 - a to str_sp
    //  sss << a.bravais_lattice_lattice_type << "," << a.bravais_conventional_lattice_lattice_type << endl;
    return sss.str();
  }
} // namespace pflow

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//START - all relevent functions for loading entries here
//Added by Corey Oses - May 2017
//load entries is heavily overloaded, mostly to accommodate entries separated as
//vector<vector<vector<> > > entries (unaries vs. binaries, then species-specific, good for convex hull),
//vector<vector<> > entries (unaries vs. binaries), OR
//vector<> entries (all together)
// ***************************************************************************
// pflow::loadEntries(aurostd::xoption& vpflow,vector<string>& velements,string server,
// ofstream& FileMESSAGE, ostream& oss)
// ***************************************************************************
namespace pflow {
bool loadEntries(
    vector<string>& velements, vector<vector<vector<aflowlib::_aflowlib_entry> > >& entries, ostream& oss) {  // overload
  ofstream FileMESSAGE;
  return loadEntries(velements, entries, FileMESSAGE, oss);
}
bool loadEntries(
    vector<string>& velements, vector<vector<vector<aflowlib::_aflowlib_entry> > >& entries, ofstream& FileMESSAGE, ostream& oss) {  // overload
  string server;
  if(XHOST.vflag_control.flag("AFLOWLIB_SERVER")) {
    server = XHOST.vflag_control.getattachedscheme("AFLOWLIB_SERVER");
  } else {
    server = "aflowlib.duke.edu";
  }
  return loadEntries(velements, server, entries, FileMESSAGE, oss);
}
bool loadEntries(
    vector<string>& velements, string server, vector<vector<vector<aflowlib::_aflowlib_entry> > >& entries, ostream& oss) {  // overload
  ofstream FileMESSAGE;
  return loadEntries(velements, server, entries, FileMESSAGE, oss);
}
bool loadEntries(
    vector<string>& velements, string server,
    vector<vector<vector<aflowlib::_aflowlib_entry> > >& entries,
    ofstream& FileMESSAGE, ostream& oss) {  // overload
  string soliloquy = "pflow::loadEntries():";
  aurostd::xoption vpflow;
  pflow::defaultLoadEntriesFlags(vpflow, FileMESSAGE, oss, std::string("A"));
  return loadEntries(vpflow, velements, server, entries, FileMESSAGE, oss);
}
bool loadEntries(
    aurostd::xoption& vpflow, vector<string>& velements,
    vector<vector<vector<aflowlib::_aflowlib_entry> > >& entries,
    ostream& oss) {  // overload
  ofstream FileMESSAGE;
  return loadEntries(vpflow, velements, entries, FileMESSAGE, oss);
}
bool loadEntries(
    aurostd::xoption& vpflow, vector<string>& velements,
    vector<vector<vector<aflowlib::_aflowlib_entry> > >& entries,
    ofstream& FileMESSAGE, ostream& oss) {  // overload
  string server;
  if(XHOST.vflag_control.flag("AFLOWLIB_SERVER")) {
    server = XHOST.vflag_control.getattachedscheme("AFLOWLIB_SERVER");
  } else {
    server = "aflowlib.duke.edu";
  }
  return loadEntries(vpflow, velements, server, entries, FileMESSAGE, oss);
}
bool loadEntries(
    aurostd::xoption& vpflow, vector<string>& velements, string server,
    vector<vector<vector<aflowlib::_aflowlib_entry> > >& entries,
    ostream& oss) {  // overload
  ofstream FileMESSAGE;
  return loadEntries(vpflow, velements, server, entries, FileMESSAGE, oss);
}
bool loadEntries(
    aurostd::xoption& vpflow, vector<string>& velements, string server,
    vector<vector<vector<aflowlib::_aflowlib_entry> > >& entries,
    ofstream& FileMESSAGE, ostream& oss) {  // main function

  string soliloquy = "pflow::loadEntries():";
  stringstream message;
  vpflow.flag("PFLOW::LOAD_ENTRIES_COMING_FROM_LOADENTRIESX",
              true);  // silent some output

  message << "Loading entries for: " << aurostd::joinWDelimiter(velements, ",");
  pflow::logger(soliloquy, message, FileMESSAGE, oss, _LOGGER_MESSAGE_);

  if(vpflow.flag("PFLOW::LOAD_ENTRIES_ONLY_ALPHABETICAL")) {
    pflow::logger(soliloquy, "Loading alphabetized entries ONLY", FileMESSAGE, oss, _LOGGER_OPTION_);
  }

  bool load_from_common = pflow::loadFromCommon(vpflow);

  if( load_from_common ) {
      pflow::logger(soliloquy, "Loading entries from COMMON", FileMESSAGE, oss, _LOGGER_OPTION_);
  } else {
    if(server == "materials.duke.edu") {
      pflow::logger(soliloquy, "Using materials.duke.edu as server", FileMESSAGE, oss, _LOGGER_OPTION_);
    } else if(server == "aflowlib.duke.edu") {
      pflow::logger(soliloquy, "Using aflowlib.duke.edu as server", FileMESSAGE, oss, _LOGGER_OPTION_);
    } else {
      pflow::logger(soliloquy, "Server must be either materials.duke.edu or aflowlib.duke.edu", FileMESSAGE, oss, _LOGGER_ERROR_);
      return false;  //entries;
    }
  }

  string lib_name;
  string lib_count_string;
  vector<vector<string> > combinations;
  vector<vector<vector<aflowlib::_aflowlib_entry> > > _entries; //unsorted

  //////////////////////////////////////////////////////////////////////////////
  // START loadLIBX
  //////////////////////////////////////////////////////////////////////////////

  for (uint lib = 0; lib < velements.size() && lib <= _AFLOW_LIB_MAX_; lib++) {
    lib_count_string = aurostd::utype2string(lib + 1);
    lib_name = std::string("LIB") + lib_count_string;
    if(vpflow.flag("PFLOW::LOAD_ENTRIES_LOAD_" + lib_name)) {
      pflow::logger(soliloquy, "Loading " + lib_name, FileMESSAGE, oss, _LOGGER_MESSAGE_);
      combinations = pflow::elementalCombinations(velements, lib + 1);
      for (uint i = 0; i < combinations.size(); i++) {
        if(!loadAndMergeLIBX(vpflow, combinations[i], lib_count_string, server, _entries, FileMESSAGE, oss)) {
          pflow::logger(soliloquy, "Merging entries for " + lib_name + " failed", FileMESSAGE, oss, _LOGGER_ERROR_);
          return false;
        }
      }
    }
  }

  //////////////////////////////////////////////////////////////////////////////
  // END loadLIBX
  //////////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////////
  // START loadICSD
  //////////////////////////////////////////////////////////////////////////////

  lib_name = std::string("ICSD");
  if(vpflow.flag("PFLOW::LOAD_ENTRIES_LOAD_" + lib_name)) {
    pflow::logger(soliloquy, "Loading " + lib_name, FileMESSAGE, oss, _LOGGER_MESSAGE_);
    if(!loadAndMergeLIBX(vpflow, velements, lib_name, server, _entries, FileMESSAGE, oss)) {
      pflow::logger(soliloquy, "Merging entries for " + lib_name + " failed", FileMESSAGE, oss, _LOGGER_ERROR_);
      return false;
    }
  }

  //////////////////////////////////////////////////////////////////////////////
  // END loadICSD
  //////////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////////
  // START Sort
  //////////////////////////////////////////////////////////////////////////////

  bool found;
  vector<string> vspecies;
  //make space and sort correctly
  for(uint i=0;i<velements.size();i++){
    entries.push_back(vector<vector<aflowlib::_aflowlib_entry> >(0));
    combinations = pflow::elementalCombinations(velements, i+1);
    for(uint j=0;j<combinations.size();j++){
      entries[i].push_back(vector<aflowlib::_aflowlib_entry>(0));
    }
    for(uint j=0;j<combinations.size();j++){
      found=false;
      if(_entries.size()>i){
        for(uint k=0;k<_entries[i].size() && !found;k++){
          if(_entries[i][k].size()) {
            vspecies=_entries[i][k][0].vspecies;
            if(!vpflow.flag("PFLOW::LOAD_ENTRIES_ONLY_ALPHABETICAL")){sort(vspecies.begin(),vspecies.end());}
            if(vspecies == combinations[j]){
              entries[i][j]=_entries[i][k];
              found=true;
            }
          }
        }
      }
    }
  }

  //////////////////////////////////////////////////////////////////////////////
  // END Sort
  //////////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////////
  // START Print loaded entries summary
  //////////////////////////////////////////////////////////////////////////////

  vector<uint> sizes;
  uint totalNum = 0;

  for (uint i = 0; i < entries.size(); i++) {
    sizes.push_back(0);
    for (uint j = 0; j < entries.at(i).size(); j++) {
      for (uint k = 0; k < entries.at(i).at(j).size(); k++) {
        sizes.at(i)++;
      }
    }
    totalNum += sizes.at(i);
  }

  for (uint i = 0; i < entries.size(); i++) {
    message << "Loaded ";
    message << entries.at(i).size() << " ";
    if(i == 0) {
      message << "unaries: ";
    } else if(i == 1) {
      message << "binaries: ";
    } else if(i == 2) {
      message << "ternaries: ";
    } else {
      message << (i + 1) << "-naries: ";
    }
    message << sizes.at(i) << " entries";
    pflow::logger(soliloquy, message, FileMESSAGE, oss, _LOGGER_MESSAGE_);
  }
  message << "Loaded " << totalNum << " entries total";
  pflow::logger(soliloquy, message, FileMESSAGE, oss, _LOGGER_MESSAGE_);

  //////////////////////////////////////////////////////////////////////////////
  // END Print loaded entries summary
  //////////////////////////////////////////////////////////////////////////////

  return true;  //entries;
}
}  // namespace pflow

namespace pflow {
// ***************************************************************************
// pflow::loadEntries(aurostd::xoption& vpflow,vector<string>& velements,string
// server,ostream& oss,ofstream& FileMESSAGE)
// ***************************************************************************
bool loadEntries(
    vector<string>& velements, vector<vector<aflowlib::_aflowlib_entry> >& entries, ostream& oss) {  // overload
  ofstream FileMESSAGE;
  return loadEntries(velements, entries, FileMESSAGE, oss);
}
bool loadEntries(
    vector<string>& velements,
    vector<vector<aflowlib::_aflowlib_entry> >& entries,
    ofstream& FileMESSAGE, ostream& oss) {  // overload
  string server;
  if(XHOST.vflag_control.flag("AFLOWLIB_SERVER")) {
    server = XHOST.vflag_control.getattachedscheme("AFLOWLIB_SERVER");
  } else {
    server = "aflowlib.duke.edu";
  }
  return loadEntries(velements, server, entries, FileMESSAGE, oss);
}
bool loadEntries(
    vector<string>& velements, string server,
    vector<vector<aflowlib::_aflowlib_entry> >& entries,
    ostream& oss) {  // overload
  ofstream FileMESSAGE;
  return loadEntries(velements, server, entries, FileMESSAGE, oss);
}
bool loadEntries(
    vector<string>& velements, string server,
    vector<vector<aflowlib::_aflowlib_entry> >& entries,
    ofstream& FileMESSAGE, ostream& oss) {  // overload
  string soliloquy = "pflow::loadEntries():";
  aurostd::xoption vpflow;
  pflow::defaultLoadEntriesFlags(vpflow, FileMESSAGE, oss, std::string("A"));
  return loadEntries(vpflow, velements, server, entries, FileMESSAGE, oss);
}
bool loadEntries(
    aurostd::xoption& vpflow, vector<string>& velements,
    vector<vector<aflowlib::_aflowlib_entry> >& entries,
    ostream& oss) {  // overload
  ofstream FileMESSAGE;
  return loadEntries(vpflow, velements, entries, FileMESSAGE, oss);
}
bool loadEntries(
    aurostd::xoption& vpflow, vector<string>& velements,
    vector<vector<aflowlib::_aflowlib_entry> >& entries,
    ofstream& FileMESSAGE, ostream& oss) {  // overload
  string server;
  if(XHOST.vflag_control.flag("AFLOWLIB_SERVER")) {
    server = XHOST.vflag_control.getattachedscheme("AFLOWLIB_SERVER");
  } else {
    server = "aflowlib.duke.edu";
  }
  return loadEntries(vpflow, velements, server, entries, FileMESSAGE, oss);
}
bool loadEntries(
    aurostd::xoption& vpflow, vector<string>& velements, string server,
    vector<vector<aflowlib::_aflowlib_entry> >& entries,
    ostream& oss) {  // overload
  ofstream FileMESSAGE;
  return loadEntries(vpflow, velements, server, entries, FileMESSAGE, oss);
}
bool loadEntries(
    aurostd::xoption& vpflow, vector<string>& velements, string server,
    vector<vector<aflowlib::_aflowlib_entry> >& entries,
    ofstream& FileMESSAGE, ostream& oss) {  // main function

  vector<vector<vector<aflowlib::_aflowlib_entry> > > naries;
  if(!loadEntries(vpflow, velements, server, naries, FileMESSAGE, oss)) { return false; }
  if(!mergeEntries(entries, naries, true)) { return false; }
  return true;
}
}  // namespace pflow

namespace pflow {
// ***************************************************************************
// pflow::loadEntries(aurostd::xoption& vpflow,vector<string>
// velements,string server,ostream& oss,ofstream& FileMESSAGE)
// ***************************************************************************
bool loadEntries(
    vector<string>& velements, vector<aflowlib::_aflowlib_entry>& entries, ostream& oss) {  // overload
  ofstream FileMESSAGE;
  return loadEntries(velements, entries, FileMESSAGE, oss);
}
bool loadEntries(
    vector<string>& velements,
    vector<aflowlib::_aflowlib_entry>& entries,
    ofstream& FileMESSAGE, ostream& oss) {  // overload
  string server;
  if(XHOST.vflag_control.flag("AFLOWLIB_SERVER")) {
    server = XHOST.vflag_control.getattachedscheme("AFLOWLIB_SERVER");
  } else {
    server = "aflowlib.duke.edu";
  }
  return loadEntries(velements, server, entries, FileMESSAGE, oss);
}
bool loadEntries(
    vector<string>& velements, string server, vector<aflowlib::_aflowlib_entry>& entries, ostream& oss) {  // overload
  ofstream FileMESSAGE;
  return loadEntries(velements, server, entries, FileMESSAGE, oss);
}
bool loadEntries(
    vector<string>& velements, string server,
    vector<aflowlib::_aflowlib_entry>& entries,
    ofstream& FileMESSAGE, ostream& oss) {  // overload
  string soliloquy = "pflow::loadEntries():";
  aurostd::xoption vpflow;
  pflow::defaultLoadEntriesFlags(vpflow, FileMESSAGE, oss, std::string("A"));
  return loadEntries(vpflow, velements, server, entries, FileMESSAGE, oss);
}
bool loadEntries(
    aurostd::xoption& vpflow, vector<string>& velements,
    vector<aflowlib::_aflowlib_entry>& entries,
    ostream& oss) {  // overload
  ofstream FileMESSAGE;
  return loadEntries(vpflow, velements, entries, FileMESSAGE, oss);
}
bool loadEntries(
    aurostd::xoption& vpflow, vector<string>& velements,
    vector<aflowlib::_aflowlib_entry>& entries,
    ofstream& FileMESSAGE, ostream& oss) {  // overload
  string server;
  if(XHOST.vflag_control.flag("AFLOWLIB_SERVER")) {
    server = XHOST.vflag_control.getattachedscheme("AFLOWLIB_SERVER");
  } else {
    server = "aflowlib.duke.edu";
  }
  return loadEntries(vpflow, velements, server, entries, FileMESSAGE, oss);
}
bool loadEntries(
    aurostd::xoption& vpflow, vector<string>& velements, string server,
    vector<aflowlib::_aflowlib_entry>& entries,
    ostream& oss) {  // overload
  ofstream FileMESSAGE;
  return loadEntries(vpflow, velements, server, entries, FileMESSAGE, oss);
}
bool loadEntries(
    aurostd::xoption& vpflow, vector<string>& velements, string server,
    vector<aflowlib::_aflowlib_entry>& entries,
    ofstream& FileMESSAGE, ostream& oss) {  // main function

  vector<vector<vector<aflowlib::_aflowlib_entry> > > naries;
  if(!loadEntries(vpflow, velements, server, naries, FileMESSAGE, oss)) { return false; }
  if(!mergeEntries(entries, naries)) { return false; }
  return true;
}
}  // namespace pflow

namespace pflow {
// ***************************************************************************
// pflow::loadFromCommon(aurostd::xoption& vpflow)
// simple function for determining if we can load from common
// it's a function because we call it a few times
// ***************************************************************************
bool loadFromCommon(aurostd::xoption& vpflow){
  bool load_from_common = (!vpflow.flag("PFLOW::LOAD_API") &&
                           (aurostd::substring2bool(XHOST.hostname, "nietzsche") ||
                            aurostd::substring2bool(XHOST.hostname, "aflowlib") ||
                            aurostd::substring2bool(XHOST.hostname, "habana")));
  return load_from_common;
}
}  // namespace pflow

namespace pflow {
// ***************************************************************************
// pflow::loadAndMergeLIBX(aurostd::xoption& vpflow, vector<string>& combination, 
// string LIB, string server,vector<vector<vector<aflowlib::_aflowlib_entry> > >& naries, 
// ofstream& FileMESSAGE, ostream& oss
// ***************************************************************************
//helper function for loadEntries, loads select library, LIB can be "LIB2" or "2" or "ICSD"
bool loadAndMergeLIBX(string combination, string LIB, string server,
                      vector<vector<vector<aflowlib::_aflowlib_entry> > >& naries, ostream& oss) {
  ofstream FileMESSAGE;
  return loadAndMergeLIBX(combination, LIB, server, naries, FileMESSAGE, oss);
}
bool loadAndMergeLIBX(string _combination, string LIB, string server,
                      vector<vector<vector<aflowlib::_aflowlib_entry> > >& naries, ofstream& FileMESSAGE, ostream& oss) {
  vector<string> combination;
  combination.push_back(_combination);
  return loadAndMergeLIBX(combination, LIB, server, naries, FileMESSAGE, oss);
}
bool loadAndMergeLIBX(vector<string>& combination, string LIB, string server,
                      vector<vector<vector<aflowlib::_aflowlib_entry> > >& naries, ostream& oss) {
  ofstream FileMESSAGE;
  return loadAndMergeLIBX(combination, LIB, server, naries, FileMESSAGE, oss);
}
bool loadAndMergeLIBX(vector<string>& combination, string LIB, string server,
                      vector<vector<vector<aflowlib::_aflowlib_entry> > >& naries, ofstream& FileMESSAGE, ostream& oss) {
  string soliloquy = "pflow::loadAndMergeLIBX():";
  aurostd::xoption vpflow;
  pflow::defaultLoadEntriesFlags(vpflow, FileMESSAGE, oss, LIB);
  return loadAndMergeLIBX(vpflow, combination, LIB, server, naries, oss);
}
bool loadAndMergeLIBX(aurostd::xoption& vpflow, string combination, string LIB, string server,
                      vector<vector<vector<aflowlib::_aflowlib_entry> > >& naries, ostream& oss) {
  ofstream FileMESSAGE;
  return loadAndMergeLIBX(vpflow, combination, LIB, server, naries, FileMESSAGE, oss);
}
bool loadAndMergeLIBX(aurostd::xoption& vpflow, string _combination, string LIB, string server,
                      vector<vector<vector<aflowlib::_aflowlib_entry> > >& naries, ofstream& FileMESSAGE, ostream& oss) {
  vector<string> combination;
  combination.push_back(_combination);
  return loadAndMergeLIBX(vpflow, combination, LIB, server, naries, FileMESSAGE, oss);
}
bool loadAndMergeLIBX(aurostd::xoption& vpflow, vector<string>& combination, string LIB, string server,
                      vector<vector<vector<aflowlib::_aflowlib_entry> > >& naries, ostream& oss) {
  ofstream FileMESSAGE;
  return loadAndMergeLIBX(vpflow, combination, LIB, server, naries, FileMESSAGE, oss);
}
bool loadAndMergeLIBX(aurostd::xoption& vpflow, vector<string>& combination, string LIB, string server,
                      vector<vector<vector<aflowlib::_aflowlib_entry> > >& naries, ofstream& FileMESSAGE, ostream& oss) {
  string soliloquy = "pflow::loadAndMergeLIBX()";
  LIB = aurostd::toupper(LIB);  //removes ALL ambiguity with case
  vector<vector<aflowlib::_aflowlib_entry> > v_temp;
  if(!loadLIBX(vpflow, LIB, combination, server, v_temp, FileMESSAGE, oss)) { return false; }
  if(vpflow.flag("PFLOW::LOAD_ENTRIES_NARIES_MINUS_ONE")) {
    for (uint i = 0; i < v_temp.size() - 1; i++) {
      if(!mergeEntries(naries, v_temp[i], false)) {  //e.g. for ternary MnPdPt, there are 3 binary combinations
        if(LIB == "ICSD") { // or LIB == "icsd") {
          pflow::logger(soliloquy, "Merging entries for ICSD (" + aurostd::utype2string(i + 1) + "-naries) failed", FileMESSAGE, oss, _LOGGER_ERROR_);
        } else {
          pflow::logger(soliloquy, "Merging entries for LIB" + LIB + " (" + aurostd::utype2string(i + 1) + "-naries) failed", FileMESSAGE, oss, _LOGGER_ERROR_);
        }
        return false;
      }
    }
  }
  if(!mergeEntries(naries, v_temp[v_temp.size() - 1], true)) {  //e.g. for ternary MnPdPt, there is only ONE ternary combination
    if(LIB == "ICSD") { // or LIB == "icsd") {
      pflow::logger(soliloquy, "Merging entries for ICSD (" + aurostd::utype2string(v_temp.size()) + "-naries) failed", FileMESSAGE, oss, _LOGGER_ERROR_);
    } else {
      pflow::logger(soliloquy, "Merging entries for LIB" + LIB + " (" + aurostd::utype2string(v_temp.size()) + "-naries) failed", FileMESSAGE, oss, _LOGGER_ERROR_);
    }
    return false;
  }
  return true;
}
}  // namespace pflow

namespace pflow {
// ***************************************************************************
// pflow::loadLIBX(aurostd::xoption& vpflow,string elements,string server,
// ofstream& FileMESSAGE, ostream& oss)
// ***************************************************************************
// there are MANY loadLIBX overloads, basically variations of string/vector elements
// and whether entries are loaded into vector<> or vector<vector<> > (unaries vs. binaries, etc.)
// loadLIBX string elements
bool loadLIBX(string LIB, string elements,
              vector<aflowlib::_aflowlib_entry>& entries,
              ostream& oss) {  // overload
  ofstream FileMESSAGE;
  return loadLIBX(LIB, elements, entries, FileMESSAGE, oss);
}
bool loadLIBX(string LIB, string elements,
              vector<aflowlib::_aflowlib_entry>& entries,
              ofstream& FileMESSAGE, ostream& oss) {  // overload
  string server;
  if(XHOST.vflag_control.flag("AFLOWLIB_SERVER")) {
    server = XHOST.vflag_control.getattachedscheme("AFLOWLIB_SERVER");
  } else {
    server = "aflowlib.duke.edu";
  }
  return loadLIBX(LIB, elements, server, entries, FileMESSAGE, oss);
}
bool loadLIBX(string LIB, string elements, string server,
              vector<aflowlib::_aflowlib_entry>& entries,
              ostream& oss) {  // overload
  ofstream FileMESSAGE;
  return loadLIBX(LIB, elements, server, entries, FileMESSAGE, oss);
}
bool loadLIBX(string LIB, string elements, string server,
              vector<aflowlib::_aflowlib_entry>& entries,
              ofstream& FileMESSAGE, ostream& oss) {  // overload
  string soliloquy = "pflow::loadLIBX():";
  aurostd::xoption vpflow;
  pflow::defaultLoadEntriesFlags(vpflow, FileMESSAGE, oss, LIB);
  return loadLIBX(vpflow, LIB, elements, server, entries, FileMESSAGE, oss);
}
bool loadLIBX(aurostd::xoption& vpflow,
              string LIB,
              string elements,
              vector<aflowlib::_aflowlib_entry>& entries,
              ostream& oss) {  // overload
  ofstream FileMESSAGE;
  return loadLIBX(vpflow, LIB, elements, entries, FileMESSAGE, oss);
}
bool loadLIBX(aurostd::xoption& vpflow,
              string LIB,
              string elements,
              vector<aflowlib::_aflowlib_entry>& entries,
              ofstream& FileMESSAGE, ostream& oss) {  // overload
  string server;
  if(XHOST.vflag_control.flag("AFLOWLIB_SERVER")) {
    server = XHOST.vflag_control.getattachedscheme("AFLOWLIB_SERVER");
  } else {
    server = "aflowlib.duke.edu";
  }
  return loadLIBX(vpflow, LIB, elements, server, entries, FileMESSAGE, oss);
}
bool loadLIBX(aurostd::xoption& vpflow,
              string LIB,
              string elements, string server,
              vector<aflowlib::_aflowlib_entry>& entries,
              ostream& oss) {  // overload
  ofstream FileMESSAGE;
  return loadLIBX(vpflow, LIB, elements, server, entries, FileMESSAGE, oss);
}
bool loadLIBX(aurostd::xoption& vpflow,
              string LIB,
              string elements, string server,
              vector<aflowlib::_aflowlib_entry>& entries,
              ofstream& FileMESSAGE, ostream& oss) {  // overload

  string soliloquy = "pflow::loadLIBX():";
  vector<string> velements =
      pflow::makeAlphabeticVector(elements, FileMESSAGE, oss);
  return loadLIBX(vpflow, LIB, velements, server, entries, FileMESSAGE, oss);
}
// ***************************************************************************
// pflow::loadLIBX(aurostd::xoption& vpflow,string LIB,vector<string>& velements,string server,
// ofstream& FileMESSAGE, ostream& oss)
// ***************************************************************************
// loadLIBX vector elements
bool loadLIBX(string LIB, vector<string>& velements,
              vector<aflowlib::_aflowlib_entry>& entries,
              ostream& oss) {  // overload
  ofstream FileMESSAGE;
  return loadLIBX(LIB, velements, entries, FileMESSAGE, oss);
}
bool loadLIBX(string LIB, vector<string>& velements,
              vector<aflowlib::_aflowlib_entry>& entries,
              ofstream& FileMESSAGE, ostream& oss) {  // overload
  string server;
  if(XHOST.vflag_control.flag("AFLOWLIB_SERVER")) {
    server = XHOST.vflag_control.getattachedscheme("AFLOWLIB_SERVER");
  } else {
    server = "aflowlib.duke.edu";
  }
  return loadLIBX(LIB, velements, server, entries, FileMESSAGE, oss);
}
bool loadLIBX(string LIB, vector<string>& velements,
              string server,
              vector<aflowlib::_aflowlib_entry>& entries,
              ostream& oss) {  // overload
  ofstream FileMESSAGE;
  return loadLIBX(LIB, velements, server, entries, FileMESSAGE, oss);
}
bool loadLIBX(string LIB, vector<string>& velements,
              string server,
              vector<aflowlib::_aflowlib_entry>& entries,
              ofstream& FileMESSAGE, ostream& oss) {  // overload
  string soliloquy = "pflow::loadLIBX():";
  aurostd::xoption vpflow;
  pflow::defaultLoadEntriesFlags(vpflow, FileMESSAGE, oss, LIB);
  return loadLIBX(vpflow, LIB, velements, server, entries, FileMESSAGE, oss);
}
bool loadLIBX(aurostd::xoption& vpflow,
              string LIB,
              vector<string>& velements,
              vector<aflowlib::_aflowlib_entry>& entries,
              ostream& oss) {  // overload
  ofstream FileMESSAGE;
  return loadLIBX(vpflow, LIB, velements, entries, FileMESSAGE, oss);
}
bool loadLIBX(aurostd::xoption& vpflow,
              string LIB,
              vector<string>& velements,
              vector<aflowlib::_aflowlib_entry>& entries,
              ofstream& FileMESSAGE, ostream& oss) {  // overload
  string server;
  if(XHOST.vflag_control.flag("AFLOWLIB_SERVER")) {
    server = XHOST.vflag_control.getattachedscheme("AFLOWLIB_SERVER");
  } else {
    server = "aflowlib.duke.edu";
  }
  return loadLIBX(vpflow, LIB, velements, server, entries, FileMESSAGE, oss);
}
bool loadLIBX(aurostd::xoption& vpflow,
              string LIB,
              vector<string>& velements,
              string server,
              vector<aflowlib::_aflowlib_entry>& entries,
              ostream& oss) {  // overload
  ofstream FileMESSAGE;
  return loadLIBX(vpflow, LIB, velements, server, entries, FileMESSAGE, oss);
}
bool loadLIBX(
    aurostd::xoption& vpflow, string LIB, vector<string>& velements, string server,
    vector<aflowlib::_aflowlib_entry>& entries,
    ofstream& FileMESSAGE, ostream& oss) {  // main function

  string soliloquy = "pflow::loadLIBX():";
  stringstream message;
  vector<vector<aflowlib::_aflowlib_entry> > naries;  //most intuitive structure from LIBs construction (unary, binaries, etc.), always start here and merge to get other variants

  if(!loadLIBX(vpflow, LIB, velements, server, naries, FileMESSAGE, oss)) { return false; }
  if(!mergeEntries(entries, naries)) { return false; }
  return true;
}
// ***************************************************************************
// pflow::loadLIBX(aurostd::xoption& vpflow,string elements,string server,
// ofstream& FileMESSAGE, ostream& oss)
// ***************************************************************************
// loadLIBX string elements
bool loadLIBX(string LIB, string elements,
              vector<vector<aflowlib::_aflowlib_entry> >& entries,
              ostream& oss) {  // overload
  ofstream FileMESSAGE;
  return loadLIBX(LIB, elements, entries, FileMESSAGE, oss);
}
bool loadLIBX(string LIB, string elements,
              vector<vector<aflowlib::_aflowlib_entry> >& entries,
              ofstream& FileMESSAGE, ostream& oss) {  // overload
  string server;
  if(XHOST.vflag_control.flag("AFLOWLIB_SERVER")) {
    server = XHOST.vflag_control.getattachedscheme("AFLOWLIB_SERVER");
  } else {
    server = "aflowlib.duke.edu";
  }
  return loadLIBX(LIB, elements, server, entries, FileMESSAGE, oss);
}
bool loadLIBX(string LIB, string elements, string server,
              vector<vector<aflowlib::_aflowlib_entry> >& entries,
              ostream& oss) {  // overload
  ofstream FileMESSAGE;
  return loadLIBX(LIB, elements, server, entries, FileMESSAGE, oss);
}
bool loadLIBX(string LIB, string elements, string server,
              vector<vector<aflowlib::_aflowlib_entry> >& entries,
              ofstream& FileMESSAGE, ostream& oss) {  // overload
  string soliloquy = "pflow::loadLIBX():";
  aurostd::xoption vpflow;
  pflow::defaultLoadEntriesFlags(vpflow, FileMESSAGE, oss, LIB);
  return loadLIBX(vpflow, LIB, elements, server, entries, FileMESSAGE, oss);
}
bool loadLIBX(aurostd::xoption& vpflow,
              string LIB,
              string elements,
              vector<vector<aflowlib::_aflowlib_entry> >& entries,
              ostream& oss) {  // overload
  ofstream FileMESSAGE;
  return loadLIBX(vpflow, LIB, elements, entries, FileMESSAGE, oss);
}
bool loadLIBX(aurostd::xoption& vpflow,
              string LIB,
              string elements,
              vector<vector<aflowlib::_aflowlib_entry> >& entries,
              ofstream& FileMESSAGE, ostream& oss) {  // overload
  string server;
  if(XHOST.vflag_control.flag("AFLOWLIB_SERVER")) {
    server = XHOST.vflag_control.getattachedscheme("AFLOWLIB_SERVER");
  } else {
    server = "aflowlib.duke.edu";
  }
  return loadLIBX(vpflow, LIB, elements, server, entries, FileMESSAGE, oss);
}
bool loadLIBX(aurostd::xoption& vpflow,
              string LIB,
              string elements, string server,
              vector<vector<aflowlib::_aflowlib_entry> >& entries,
              ostream& oss) {  // overload
  ofstream FileMESSAGE;
  return loadLIBX(vpflow, LIB, elements, server, entries, FileMESSAGE, oss);
}
bool loadLIBX(aurostd::xoption& vpflow,
              string LIB,
              string elements, string server,
              vector<vector<aflowlib::_aflowlib_entry> >& entries,
              ofstream& FileMESSAGE, ostream& oss) {  // overload

  string soliloquy = "pflow::loadLIBX():";
  vector<string> velements =
      pflow::makeAlphabeticVector(elements, FileMESSAGE, oss);
  return loadLIBX(vpflow, LIB, velements, server, entries, FileMESSAGE, oss);
}
// ***************************************************************************
// pflow::loadLIBX(aurostd::xoption& vpflow,string LIB,vector<string>& velements,string server,
// ofstream& FileMESSAGE, ostream& oss)
// ***************************************************************************
// loadLIBX vector elements
bool loadLIBX(string LIB, vector<string>& velements,
              vector<vector<aflowlib::_aflowlib_entry> >& entries,
              ostream& oss) {  // overload
  ofstream FileMESSAGE;
  return loadLIBX(LIB, velements, entries, FileMESSAGE, oss);
}
bool loadLIBX(string LIB, vector<string>& velements,
              vector<vector<aflowlib::_aflowlib_entry> >& entries,
              ofstream& FileMESSAGE, ostream& oss) {  // overload
  string server;
  if(XHOST.vflag_control.flag("AFLOWLIB_SERVER")) {
    server = XHOST.vflag_control.getattachedscheme("AFLOWLIB_SERVER");
  } else {
    server = "aflowlib.duke.edu";
  }
  return loadLIBX(LIB, velements, server, entries, FileMESSAGE, oss);
}
bool loadLIBX(string LIB, vector<string>& velements,
              string server,
              vector<vector<aflowlib::_aflowlib_entry> >& entries,
              ostream& oss) {  // overload
  ofstream FileMESSAGE;
  return loadLIBX(LIB, velements, server, entries, FileMESSAGE, oss);
}
bool loadLIBX(string LIB, vector<string>& velements,
              string server,
              vector<vector<aflowlib::_aflowlib_entry> >& entries,
              ofstream& FileMESSAGE, ostream& oss) {  // overload
  string soliloquy = "pflow::loadLIBX():";
  aurostd::xoption vpflow;
  pflow::defaultLoadEntriesFlags(vpflow, FileMESSAGE, oss, LIB);
  return loadLIBX(vpflow, LIB, velements, server, entries, FileMESSAGE, oss);
}
bool loadLIBX(aurostd::xoption& vpflow,
              string LIB,
              vector<string>& velements,
              vector<vector<aflowlib::_aflowlib_entry> >& entries,
              ostream& oss) {  // overload
  ofstream FileMESSAGE;
  return loadLIBX(vpflow, LIB, velements, entries, FileMESSAGE, oss);
}
bool loadLIBX(aurostd::xoption& vpflow,
              string LIB,
              vector<string>& velements,
              vector<vector<aflowlib::_aflowlib_entry> >& entries,
              ofstream& FileMESSAGE, ostream& oss) {  // overload
  string server;
  if(XHOST.vflag_control.flag("AFLOWLIB_SERVER")) {
    server = XHOST.vflag_control.getattachedscheme("AFLOWLIB_SERVER");
  } else {
    server = "aflowlib.duke.edu";
  }
  return loadLIBX(vpflow, LIB, velements, server, entries, FileMESSAGE, oss);
}
bool loadLIBX(aurostd::xoption& vpflow,
              string LIB,
              vector<string>& velements,
              string server,
              vector<vector<aflowlib::_aflowlib_entry> >& entries,
              ostream& oss) {  // overload
  ofstream FileMESSAGE;
  return loadLIBX(vpflow, LIB, velements, server, entries, FileMESSAGE, oss);
}
bool loadLIBX(
    aurostd::xoption& vpflow, string LIB, vector<string>& velements, string server,
    vector<vector<aflowlib::_aflowlib_entry> >& entries,
    ofstream& FileMESSAGE, ostream& oss) {  // main function

  string soliloquy = "pflow::loadLIBX():";
  stringstream message;
  //vector<vector<aflowlib::_aflowlib_entry> > entries;
  for (uint i = 0; i < velements.size(); i++) {
    entries.push_back(vector<aflowlib::_aflowlib_entry>(0));
  }

  //categorize LIB
  string lib_name;
  vector<string> symmetries;
  LIB = aurostd::toupper(LIB);  //removes ALL ambiguity with case
  bool isICSD = bool(LIB == "ICSD"); //(LIB == "ICSD" || LIB == "icsd");
  if(isICSD) {
    lib_name = LIB;
    string SYMMETRIES = "BCC BCT CUB FCC HEX MCL MCLC ORC ORCC ORCF ORCI RHL TET TRI";
    aurostd::string2tokens(SYMMETRIES, symmetries, " ");
    //add here other special LIBS, i.e., not LIB2, LIB3, etc.
  } else {
    string lib_count_string;
    //make robust so input can be "LIB2" or "2"
    if(aurostd::substring2bool(LIB, "LIB")) {
      lib_name = LIB;
      lib_count_string = aurostd::RemoveSubString(LIB, "LIB");
    } else {
      lib_name = "LIB" + LIB;
      lib_count_string = LIB;
    }
    //check validity of input
    for (uint i = 0; i < lib_count_string.size(); i++) {
      if(!isdigit(lib_count_string[i])) {
        message << "Unknown LIB specification (" << LIB << "), should be \"LIB1\", \"LIB2\", or \"1\", \"2\", etc.";
        pflow::logger(soliloquy, message, FileMESSAGE, oss, _LOGGER_ERROR_);
        return false;  //entries;
      }
    }
    uint lib_count_uint = aurostd::string2utype<int>(lib_count_string);
    if(velements.size() != lib_count_uint) {
      message << "LIB" << lib_count_uint << " loads " << lib_count_uint << "-naries ONLY";
      pflow::logger(soliloquy, message, FileMESSAGE, oss, _LOGGER_ERROR_);
      return false;  //entries;
    }
  }

  message << "Loading " << lib_name << " entries for: " << aurostd::joinWDelimiter(velements, "");
  pflow::logger(soliloquy, message, FileMESSAGE, oss, _LOGGER_MESSAGE_);

  bool load_from_common = pflow::loadFromCommon(vpflow);
  bool override_load_from_common = false;

  string LIB_path;  //will be actual directory path or URL, depending on machine
  if(load_from_common) {
    if(!vpflow.flag("PFLOW::LOAD_ENTRIES_COMING_FROM_LOADENTRIESX")) {
      pflow::logger(soliloquy, "Loading entries from COMMON", FileMESSAGE, oss, _LOGGER_OPTION_);
    }
    LIB_path = "/common/" + lib_name + "/RAW/";
    if(!aurostd::IsDirectory(LIB_path)) {
      load_from_common = false;
      message << LIB_path << " does not exist! Cannot load from COMMON, switching to API";
      pflow::logger(soliloquy, message, FileMESSAGE, oss, _LOGGER_WARNING_);
      override_load_from_common = true;
    }
  }
  if(!load_from_common) {
    if(server == "materials.duke.edu") {
      if(override_load_from_common || (!vpflow.flag("PFLOW::LOAD_ENTRIES_COMING_FROM_LOADENTRIESX"))) {
        pflow::logger(soliloquy, "Using materials.duke.edu as server", FileMESSAGE, oss, _LOGGER_OPTION_);
      }
    } else if(server == "aflowlib.duke.edu") {
      if(override_load_from_common || (!vpflow.flag("PFLOW::LOAD_ENTRIES_COMING_FROM_LOADENTRIESX"))) {
        pflow::logger(soliloquy, "Using aflowlib.duke.edu as server", FileMESSAGE, oss, _LOGGER_OPTION_);
      }
    } else {
      pflow::logger(soliloquy, "Server must be either materials.duke.edu or aflowlib.duke.edu", FileMESSAGE, oss, _LOGGER_ERROR_);
      return false;  //entries;
    }
    if(isICSD) {
      LIB_path = server + "/AFLOWDATA/ICSD_WEB/";
    } else {
      LIB_path = server + "/AFLOWDATA/" + lib_name + "_RAW/";
    }
  }

  //////////////////////////////////////////////////////////////////////////
  // START Finding matches
  //////////////////////////////////////////////////////////////////////////
  aflowlib::_aflowlib_entry _aflowlib_tmp;
  uint nary;
  uint total_count = 0;
  if(isICSD) {
    bool double_check_icsd = false;  //NOT NECESSARY for ICSD since we load from the calculation layer
    string symmetry_path, clean_icsd;
    vector<string> icsds, tokens, elements;
    for (uint i = 0; i < symmetries.size(); i++) {
      symmetry_path = LIB_path + symmetries[i];
      if(load_from_common && !aurostd::IsDirectory(symmetry_path)) { continue; }
      //not many symmetries, so this boolean placement won't impact much
      if(load_from_common) {
        aurostd::DirectoryLS(symmetry_path, icsds);
      } else {
        aurostd::url2tokens(symmetry_path + "/?aflowlib_entries", icsds, ",");
      }
      for (uint j = 0; j < icsds.size(); j++) {
        aurostd::string2tokens(icsds[j], tokens, "_");   //assumes no '_' before '_ICSD_'
        if(tokens.size() > 2 && tokens[1] == "ICSD") {  //just a quick check
          clean_icsd = tokens[0];
          if(compoundsBelong(velements, clean_icsd,
                                   FileMESSAGE, oss, 
                                   !vpflow.flag("PFLOW::LOAD_ENTRIES_ONLY_ALPHABETICAL"))) {  // NOT recommended, these are BS entries
            // url2aflowlib has bad printing to screen options, so we mimic here

            if(vpflow.flag("PFLOW::LOAD_ENTRIES_ENTRY_OUTPUT")) {
              message << "Loading " << (load_from_common ? "path" : "url") << "=" << symmetry_path << "/" << icsds[j];
              pflow::logger(std::string("aurostd::") + std::string((load_from_common ? "file" : "url")) + std::string("2string():"), message, FileMESSAGE, oss, _LOGGER_MESSAGE_);
            }
            //complicated ternary operator, returns bool, but necessary to avoid double code
            if(
                (load_from_common)
                    ?
                    //if load_from_common, check this bool
                    (
                        //(aurostd::FileExist(symmetry_path + "/" + icsds[j] + "/aflowlib.out")) &&
                        (aurostd::FileNotEmpty(symmetry_path + "/" + icsds[j] + "/aflowlib.out")) &&
                        (_aflowlib_tmp.file2aflowlib(symmetry_path + "/" + icsds[j] + "/aflowlib.out", message) > 0) &&
                        (double_check_icsd ? compoundsBelong(velements, _aflowlib_tmp.vspecies) : true)  //sometimes we find odd entries in the wrong LIBS, better to be safe, NOT NECESSARY for ICSD since we load from the calculation layer
                        )
                    :
                    //if load_from_api, check this bool
                    (
                        (_aflowlib_tmp.url2aflowlib(symmetry_path + "/" + icsds[j], message, false) > 0) &&
                        (double_check_icsd ? compoundsBelong(velements, _aflowlib_tmp.vspecies) : true)  //sometimes we find odd entries in the wrong LIBS, better to be safe, NOT NECESSARY for ICSD since we load from the calculation layer
                        )) {
              nary = _aflowlib_tmp.vspecies.size();
              if(entries.size() < nary) { continue; }  //this entry is garbage (wrong directory)
              entries[nary - 1].push_back(_aflowlib_tmp);
              total_count++;
              if(vpflow.flag("PFLOW::LOAD_ENTRIES_LOAD_XSTRUCTURES")) {
                if(!loadXstructures(entries[nary - 1].back(),FileMESSAGE,oss,vpflow.flag("PFLOW::LOAD_ENTRIES_LOAD_XSTRUCTURES_RELAXED_ONLY"))) {  //CO 171202 - let it decided whether to load from common or not
                  entries[nary - 1].pop_back();
                  total_count--;
                }
              }
            }
            message.str("");
          }
        }
      }
    }
  } else {
    bool double_check_lib = true;  //while NOT NECESSARY for ICSD, LIBs are screwed up sometimes
    string system_path, calculation_path, clean_calculation;
    vector<string> systems, calculations;
    if(load_from_common) {
      aurostd::DirectoryLS(LIB_path, systems);
    } else {
      aurostd::url2tokens(LIB_path + "/?aflowlib_entries", systems, ",");
    }
    for (uint i = 0; i < systems.size(); i++) {
      calculation_path = LIB_path + systems[i];
      if(load_from_common && !aurostd::IsDirectory(calculation_path)) { continue; }
      if(compoundsBelong(velements, systems[i],
                               FileMESSAGE, oss,
                               !vpflow.flag("PFLOW::LOAD_ENTRIES_ONLY_ALPHABETICAL"))) {  // NOT recommended, these are BS entries
        if(load_from_common) {
          aurostd::DirectoryLS(calculation_path, calculations);
        } else {
          aurostd::url2tokens(calculation_path + "/?aflowlib_entries", calculations, ",");
        }
        for (uint j = 0; j < calculations.size(); j++) {
          if(load_from_common && !aurostd::IsDirectory(calculation_path + "/" + calculations[j])) { continue; }
          if(vpflow.flag("PFLOW::LOAD_ENTRIES_ENTRY_OUTPUT")) {
            message << "Loading " << (load_from_common ? "path" : "url") << "=" << calculation_path << "/" << calculations[j];
            pflow::logger(std::string("aurostd::") + std::string((load_from_common ? "file" : "url")) + std::string("2string():"), message, FileMESSAGE, oss, _LOGGER_MESSAGE_);
          }
          //complicated ternary operator, returns bool, but necessary to avoid double code
          if(
              (load_from_common)
                  ?
                  //if load_from_common, check this bool
                  (
                      //(aurostd::FileExist(calculation_path + "/" + calculations[j] + "/aflowlib.out")) &&
                      (aurostd::FileNotEmpty(calculation_path + "/" + calculations[j] + "/aflowlib.out")) &&
                      (_aflowlib_tmp.file2aflowlib(calculation_path + "/" + calculations[j] + "/aflowlib.out", message) > 0) &&
                      (double_check_lib ? compoundsBelong(velements, _aflowlib_tmp.vspecies) : true)  //sometimes we find odd entries in the wrong LIBS, better to be safe
                      )
                  :
                  //if load_from_api, check this bool
                  (
                      (_aflowlib_tmp.url2aflowlib(calculation_path + "/" + calculations[j], message, false) > 0) &&
                      (double_check_lib ? compoundsBelong(velements, _aflowlib_tmp.vspecies) : true)  //sometimes we find odd entries in the wrong LIBS, better to be safe
                      )) {
            nary = _aflowlib_tmp.vspecies.size();
            if(entries.size() < nary) { continue; }  //this entry is garbage (wrong directory)
            entries[nary - 1].push_back(_aflowlib_tmp);
            total_count++;
            if(vpflow.flag("PFLOW::LOAD_ENTRIES_LOAD_XSTRUCTURES")) {
              if(!loadXstructures(entries[nary - 1].back(),FileMESSAGE,oss,vpflow.flag("PFLOW::LOAD_ENTRIES_LOAD_XSTRUCTURES_RELAXED_ONLY"))) {  //CO 171202 - let it decided whether to load from common or not
                entries[nary - 1].pop_back();
                total_count--;
              }
            }
          }
          message.str("");
        }
      }
    }
  }
  //////////////////////////////////////////////////////////////////////////
  // END Finding matches
  //////////////////////////////////////////////////////////////////////////

  // Raises PureA and PureB flags
  if(velements.size() == 2) {
    if(entries.size() && entries[0].size()) {
      for (uint j = 0; j < entries[0].size(); j++) {
        if(entries[0][j].vstoichiometry.size() == 1) {
          if(entries[0][j].species == velements.at(0)) {
            entries[0][j].pureA = true;
          } else if(entries[0][j].species == velements.at(1)) {
            entries[0][j].pureB = true;
          }
        }
      }
    }
  }

  //////////////////////////////////////////////////////////////////////////////
  // START Print loaded entries summary
  //////////////////////////////////////////////////////////////////////////////

  message << "Loaded " << total_count << " entries";
  pflow::logger(soliloquy, message, FileMESSAGE, oss, _LOGGER_MESSAGE_);

  //////////////////////////////////////////////////////////////////////////////
  // END Print loaded entries summary
  //////////////////////////////////////////////////////////////////////////////

  return true;  //entries;
}
}  // namespace pflow

namespace pflow {
// ***************************************************************************
// mergeEntries()
// ***************************************************************************
//simple helper function for loading multiple libraries together, will take
//combinations of nested vectors and convert them to other nested vectors
//naries is output, new_entries is input
//these assume vector<aflowlib::_aflowlib_entry> new_entries is all of the same
//type, e.g., binaries of MnPd (e.g., MnPd2, Mn2Pd, etc., similar structure to LIBS)
//also assumes ordered vector<vector<vector<aflowlib::_aflowlib_entry> > >& naries
//outer most - unary, binary, etc.
//next outer - species match, Mn, Pd, MnPd, etc.
//inner most - entries
bool mergeEntries(vector<vector<vector<aflowlib::_aflowlib_entry> > >& naries,
                  vector<vector<vector<aflowlib::_aflowlib_entry> > >& new_entries) {
  for (uint i = 0; i < new_entries.size(); i++) {
    if(!mergeEntries(naries, new_entries[i], true)) { return false; }  //structured data
  }
  return true;
}

bool mergeEntries(vector<vector<vector<aflowlib::_aflowlib_entry> > >& naries,
                  vector<vector<aflowlib::_aflowlib_entry> >& new_entries, bool assume_same_type) {
  for (uint i = 0; i < new_entries.size(); i++) {
    if(naries.size() <= i + 1) {
      naries.push_back(vector<vector<aflowlib::_aflowlib_entry> >(0));
    }
    if(!mergeEntries(naries[i], new_entries[i], assume_same_type, true)) {  //triple vector<> naries implies this structure
      return false;
    }
  }
  return true;
}

bool mergeEntries(vector<vector<vector<aflowlib::_aflowlib_entry> > >& naries,
                  vector<aflowlib::_aflowlib_entry>& new_entries, bool assume_same_type) {
  if(!new_entries.size()) { return true; }
  int match1, match2;
  if(assume_same_type) {
    if(!mergeEntries(naries, new_entries[0], match1, match2)) { return false; }
    for (uint i = 1; i < new_entries.size(); i++) { naries[match1][match2].push_back(new_entries[i]); }
  } else {
    for (uint i = 0; i < new_entries.size(); i++) {
      if(!new_entries[i].vspecies.size()) { return false; }  //what the heck is this?
      if(!mergeEntries(naries, new_entries[i], match1, match2)) { return false; }
    }
  }
  return true;
}

bool mergeEntries(vector<vector<vector<aflowlib::_aflowlib_entry> > >& naries,
                  aflowlib::_aflowlib_entry& new_entry) {
  int match1 = -1;
  int match2 = -1;
  return mergeEntries(naries, new_entry, match1, match2);
}
bool mergeEntries(vector<vector<vector<aflowlib::_aflowlib_entry> > >& naries,
                  aflowlib::_aflowlib_entry& new_entry, int& match1, int& match2) {
  if(!new_entry.vspecies.size()) { return false; }    //what the heck is this?
  while (naries.size() < new_entry.vspecies.size()) {  //assumes new_entries all have the same vspecies.size()
    naries.push_back(vector<vector<aflowlib::_aflowlib_entry> >(0));
  }
  match1 = new_entry.vspecies.size() - 1;
  return mergeEntries(naries[match1], new_entry, match2, true);  //triple vector<> naries implies this structure
}

//naries takes on two forms depending on sort_by_species
//if sort_by_species==true, then naries is truly a single nary (unary, binary, etc.) with
//the second layer consisting of different species
//otherwise, naries is the total entries (similar to naries above), where second layer
//is unaries, binary, etc. (no layer with different species)
//if sort_by_species and we are coming from LIBs (assume_same_type==true), then we don't need to check every entry, we already
//know they have the same type (binary of same species)
bool mergeEntries(vector<vector<aflowlib::_aflowlib_entry> >& naries,
                  vector<vector<aflowlib::_aflowlib_entry> >& new_entries, bool assume_same_type,
                  bool sort_by_species) {
  for (uint i = 0; i < new_entries.size(); i++) {
    if(!mergeEntries(naries, new_entries[i], assume_same_type, sort_by_species)) { return false; }
  }
  return true;
}

bool mergeEntries(vector<vector<aflowlib::_aflowlib_entry> >& naries,
                  vector<aflowlib::_aflowlib_entry>& new_entries, bool assume_same_type, bool sort_by_species) {
  if(!new_entries.size()) { return true; }
  int match;
  if(assume_same_type) {
    if(!mergeEntries(naries, new_entries[0], match, sort_by_species)) { return false; }
    for (uint i = 1; i < new_entries.size(); i++) { naries[match].push_back(new_entries[i]); }
  } else {
    for (uint i = 0; i < new_entries.size(); i++) {
      if(!mergeEntries(naries, new_entries[i], match, sort_by_species)) { return false; }
    }
  }
  return true;
}

//naries takes on two forms depending on sort_by_species
//if sort_by_species==true, then naries is truly a single nary (unary, binary, etc.) with
//the second layer consisting of different species
//otherwise, naries is the total entries (similar to naries above), where second layer
//is unaries, binary, etc. (no layer with different species)
bool mergeEntries(vector<vector<aflowlib::_aflowlib_entry> >& naries,
                  aflowlib::_aflowlib_entry& new_entry, bool sort_by_species) {
  int match = -1;
  return mergeEntries(naries, new_entry, match, sort_by_species);
}

bool mergeEntries(vector<vector<aflowlib::_aflowlib_entry> >& naries,
                  aflowlib::_aflowlib_entry& new_entry, int& match, bool sort_by_species) {
  if(!new_entry.vspecies.size()) { return false; }
  match = -1;
  if(sort_by_species) {
    //all of naries is unary, binary, etc., now just need to match species
    for (uint i = 0; i < naries.size() && match == -1; i++) {
      //test of stupidity
      if(naries[i][0].vspecies.size() != new_entry.vspecies.size()) { return false; }
      if(naries[i][0].vspecies == new_entry.vspecies) { match = i; }
    }
    if(match == -1) {
      naries.push_back(vector<aflowlib::_aflowlib_entry>(0));
      match = naries.size() - 1;
    }
  } else {
    //just need to create space for unary, binary, etc.
    while (naries.size() < new_entry.vspecies.size()) {
      naries.push_back(vector<aflowlib::_aflowlib_entry>(0));
    }
    match = new_entry.vspecies.size() - 1;
  }
  naries[match].push_back(new_entry);
  return true;
}

//start combining
bool mergeEntries(vector<vector<aflowlib::_aflowlib_entry> >& naries,
                  vector<vector<vector<aflowlib::_aflowlib_entry> > >& new_entries, bool sort_by_species) {
  for (uint i = 0; i < new_entries.size(); i++) {
    if(!mergeEntries(naries, new_entries[i], true, sort_by_species)) { return false; }
  }
  return true;
}

bool mergeEntries(vector<aflowlib::_aflowlib_entry>& naries,
                  vector<vector<vector<aflowlib::_aflowlib_entry> > >& new_entries) {
  for (uint i = 0; i < new_entries.size(); i++) {
    for (uint j = 0; j < new_entries[i].size(); j++) {
      if(!mergeEntries(naries, new_entries[i][j])) { return false; }
    }
  }
  return true;
}

bool mergeEntries(vector<aflowlib::_aflowlib_entry>& naries,
                  vector<vector<aflowlib::_aflowlib_entry> >& new_entries) {
  for (uint i = 0; i < new_entries.size(); i++) {
    if(!mergeEntries(naries, new_entries[i])) { return false; }
  }
  return true;
}

//trivial
bool mergeEntries(vector<aflowlib::_aflowlib_entry>& naries,
                  vector<aflowlib::_aflowlib_entry>& new_entries) {
  for (uint i = 0; i < new_entries.size(); i++) {
    if(!mergeEntries(naries, new_entries[i])) { return false; }
  }
  return true;
}

//trivial
bool mergeEntries(vector<aflowlib::_aflowlib_entry>& naries,
                  aflowlib::_aflowlib_entry& new_entry) {
  naries.push_back(new_entry);
  return true;
}
}  // namespace pflow

namespace pflow {
// ***************************************************************************
// pflow::getCombination(vector<string>& velements,uint nary)
// ***************************************************************************
// for given set of elements, will return nary combinations
// binary combinations of MnPdPt: MnPd, MnPt, PdPt
void getCombination(vector<string>& velements, vector<string>& combination,
                    vector<vector<string> >& combinations, uint offset, uint nary) {
  if(!nary) {
    combinations.push_back(combination);
    return;
  }
  for (uint i = offset; i <= velements.size() - nary; i++) {
    combination.push_back(velements[i]);
    getCombination(velements, combination, combinations, i + 1, nary - 1);
    combination.pop_back();
  }
}

vector<vector<string> > elementalCombinations(vector<string>& velements,
                                              uint nary) {
  vector<string> combination;
  vector<vector<string> > combinations;
  pflow::getCombination(velements, combination, combinations, 0, nary);
  return combinations;
}
}  // namespace pflow

namespace pflow {
// ***************************************************************************
// pflow::compoundsBelong(vector<string>& velements, vector<string>& elements)
// ***************************************************************************
// let's you know if input (or elements) belongs on hull of velements
// if sort_elements==True, will sort(elements) first before matching with velements, 
// sorting is generally NOT preferred as it would match unsorted compounds from database (NOT good)
// NOTE, this will NOT sort velements, it is assumed for speed
bool compoundsBelong(vector<string>& velements, string input, ostream& oss, bool sort_input) {
  ofstream FileMESSAGE;
  return compoundsBelong(velements, input, FileMESSAGE, oss, sort_input);
}
bool compoundsBelong(vector<string>& velements, string input,
                          ofstream& FileMESSAGE, ostream& oss, bool sort_input) {
  vector<string> elements;
  elements = stringElements2VectorElements(input, FileMESSAGE, oss);
  return compoundsBelong(velements, elements, sort_input);
}
bool compoundsBelong(vector<string>& velements, vector<string> elements, bool sort_elements) {
  //simply check if all elements are in velements (in order, sort if necessary)
  if(sort_elements){sort(elements.begin(),elements.end());}
  bool found;
  uint starting_index = 0;  //ensures we search in order
  for (uint i = 0; i < elements.size(); i++) {
    found = false;
    for (uint j = starting_index; j < velements.size() && !found; j++) {
      if(elements[i] == velements[j]) { found = true; starting_index = j+1; }  //start search at the next velement
    }
    if(!found) { return false; }
  }
  return true;
}
}  // namespace pflow

// functions for loading entries
namespace pflow {
// ***************************************************************************
// pflow::loadXstructures(aflowlib::_aflowlib_entry& entry,string path,bool
  // relaxed_only,ostream& oss,ofstream& FileMESSAGE)
// ***************************************************************************
// loads xstructures
bool loadXstructures(aflowlib::_aflowlib_entry& entry, ostream& oss, bool relaxed_only, string path, bool is_url_path) {
    ofstream FileMESSAGE;
  return loadXstructures(entry,FileMESSAGE,oss,relaxed_only,path,is_url_path);
  }
bool loadXstructures(aflowlib::_aflowlib_entry& entry, ofstream& FileMESSAGE, ostream& oss, bool relaxed_only, string path, bool is_url_path) {
  bool LDEBUG=(FALSE || XHOST.DEBUG);
    string soliloquy = "pflow::loadXstructures():";
    stringstream message;

  //get path if not provided
    if (path.empty()) {
    path = entry.getPathAURL(FileMESSAGE,oss,true); is_url_path=false;
    if(!(!path.empty()&&aurostd::IsDirectory(path))){path = "";}  //reset
    }
  if(path.empty()){path = entry.getPathAURL(FileMESSAGE,oss,false);is_url_path=true;}
    if (path.empty()) {
      message << "Path cannot be loaded from entry, skipping loadXstructure()";
    pflow::logger(soliloquy, message, FileMESSAGE, oss, _LOGGER_WARNING_);
      return false;
    }

  xstructure xstrAux;
    stringstream ss;
  vector<string> files;
  if(LDEBUG){cerr << soliloquy << " path=" << path << endl;}
    if(is_url_path){aurostd::url2tokens(path + "/?files", files, ",");}
    else{aurostd::DirectoryLS(path,files);}
  if((!aurostd::substring2bool(files, "POSCAR.relax1") &&
          !aurostd::substring2bool(files, "POSCAR.orig") && !relaxed_only) ||
      (!aurostd::substring2bool(files, "POSCAR.relax2") &&
         !aurostd::substring2bool(files, "CONTCAR.relax1") && !relaxed_only) ||
      (!aurostd::substring2bool(files, "CONTCAR.relax") &&
       !aurostd::substring2bool(files, "CONTCAR.relax2"))) {
      message << "path=" << path << " missing structure files. Ignoring entry";
    pflow::logger(soliloquy, message, FileMESSAGE, oss, _LOGGER_WARNING_);
    return false;
  } else {
      if (!relaxed_only) {
      //////////////////////////////////////////////////////////////////////////
      // START Loading original structures
      //////////////////////////////////////////////////////////////////////////

      if(!xstrAux.atoms.size() && aurostd::substring2bool(files, "POSCAR.orig")) {
          ss.str("");
        if ( (is_url_path ? 
              aurostd::url2stringstream(path + "/POSCAR.orig",ss,false) : 
                aurostd::file2stringstream(path + "/POSCAR.orig",ss)) ) {
        xstrAux = xstructure(ss, IOVASP_AUTO);
      }
      }
      if(!xstrAux.atoms.size() && aurostd::substring2bool(files, "POSCAR.relax1")) {
          ss.str("");
        if ( (is_url_path ? 
              aurostd::url2stringstream(path + "/POSCAR.relax1",ss,false) :
                aurostd::file2stringstream(path + "/POSCAR.relax1",ss)) ) {
        xstrAux = xstructure(ss, IOVASP_AUTO);
      }
        }
        if (!xstrAux.atoms.size()) {
        pflow::logger(soliloquy, "Cannot load original structure", FileMESSAGE, oss, _LOGGER_WARNING_);
          return false;
      }
      entry.vstr.push_back(xstrAux);
      xstrAux.Clear();
      if(LDEBUG){cerr << soliloquy << " loaded ORIGINAL structure" << endl;}

      //////////////////////////////////////////////////////////////////////////
      // END Loading original structures
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      // START Loading singularly-relaxed structures
      //////////////////////////////////////////////////////////////////////////

      if(!xstrAux.atoms.size() && aurostd::substring2bool(files, "POSCAR.relax2")) {
          ss.str("");
        if ( (is_url_path ? 
              aurostd::url2stringstream(path + "/POSCAR.relax2",ss,false) :
                aurostd::file2stringstream(path + "/POSCAR.relax2",ss)) ) {
        xstrAux = xstructure(ss, IOVASP_AUTO);
      }
      }
      if(!xstrAux.atoms.size() && aurostd::substring2bool(files, "CONTCAR.relax1")) {
          ss.str("");
        if ( (is_url_path ? 
              aurostd::url2stringstream(path + "/CONTCAR.relax1",ss,false) :
                aurostd::file2stringstream(path + "/CONTCAR.relax1",ss)) ) {
        xstrAux = xstructure(ss, IOVASP_AUTO);
      }
        }
        if (!xstrAux.atoms.size()) {
        pflow::logger(soliloquy, "Cannot load mid-relaxed structure", FileMESSAGE, oss, _LOGGER_WARNING_);
          return false;
      }
      entry.vstr.push_back(xstrAux);
      xstrAux.Clear();
      if(LDEBUG){cerr << soliloquy << " loaded RELAX1 structure" << endl;}

      //////////////////////////////////////////////////////////////////////////
      // END Loading singularly-relaxed structures
      //////////////////////////////////////////////////////////////////////////
    }

    ////////////////////////////////////////////////////////////////////////////
    // START Loading fully-relaxed structures
    ////////////////////////////////////////////////////////////////////////////

    if(!xstrAux.atoms.size() && aurostd::substring2bool(files, "CONTCAR.relax")) {
        ss.str("");
      if ( (is_url_path ? 
            aurostd::url2stringstream(path + "/CONTCAR.relax",ss,false): 
              aurostd::file2stringstream(path + "/CONTCAR.relax",ss)) ) {
      xstrAux = xstructure(ss, IOVASP_AUTO);
    }
    }
    if(!xstrAux.atoms.size() && aurostd::substring2bool(files, "CONTCAR.relax2")) {
        ss.str("");
      if ( (is_url_path ? 
            aurostd::url2stringstream(path + "/CONTCAR.relax2",ss,false) :
              aurostd::file2stringstream(path + "/CONTCAR.relax2",ss)) ) {
      xstrAux = xstructure(ss, IOVASP_AUTO);
    }
      }
      if (!xstrAux.atoms.size()) {
      pflow::logger(soliloquy, "Cannot load fully-relaxed structure", FileMESSAGE, oss, _LOGGER_WARNING_);
        return false;
    }
    entry.vstr.push_back(xstrAux);
    xstrAux.Clear();
    if(LDEBUG){cerr << soliloquy << " loaded FULLY-RELAXED structure" << endl;}

    ////////////////////////////////////////////////////////////////////////////
    // END Loading fully-relaxed structures
    ////////////////////////////////////////////////////////////////////////////

    return true;
  }
}
}  // namespace pflow

// functions for making input alphabetic
namespace pflow {
// ***************************************************************************
// pflow::stringElements2VectorElements(string input,ostream&
// oss,ofstream& FileMESSAGE)
// ***************************************************************************
// returns UNSORTED vector<string> from string
vector<string> stringElements2VectorElements(string input,
                                             ostream& oss, bool clean) {  // overload
  ofstream FileMESSAGE;
  return stringElements2VectorElements(input, FileMESSAGE, oss, clean);
}
vector<string> stringElements2VectorElements(
    string input, ofstream& FileMESSAGE, ostream& oss, bool clean) {  // main function
  string soliloquy = "pflow::stringElements2VectorElements():";
  vector<string> velements;

  //////////////////////////////////////////////////////////////////////////////
  // START Checks for correct input by counting number of uppercase letters
  //////////////////////////////////////////////////////////////////////////////

  if(input.empty()) {
    pflow::logger(soliloquy, "Empty input", FileMESSAGE, oss, _LOGGER_ERROR_);
    return velements;
  }
  if(islower(input[0])) {
    pflow::logger(soliloquy, "Elements must be properly capitalized", FileMESSAGE, oss, _LOGGER_ERROR_);
    return velements;
  }
  uint numberOfElements = 0;
  for (uint i = 0; i < input.size(); i++) {
    if(isupper(input[i])) {
      numberOfElements++;
    }
  }
  if(numberOfElements == 0) {
    pflow::logger(soliloquy, "Elements must be properly capitalized", FileMESSAGE, oss, _LOGGER_ERROR_);
    return velements;
  }

  //////////////////////////////////////////////////////////////////////////////
  // END Checks for correct input by counting number of uppercase letters
  //////////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////////
  // START Parsing input
  //////////////////////////////////////////////////////////////////////////////

  for (uint i = 0; i < numberOfElements; i++) {
    velements.push_back("");
  }
  uint counter = 0;
  for (uint i = 0; i < input.size(); i++) {
    if(isupper(input[i])) {
      char aux[8];
      int j = 1;
      aux[0] = input[i];
      //while (!isupper(input[i + j]) && (i + j) < input.size()) {  //needs to be more robust to accommodate for PP names, _ICSD_, etc.
      while ((((i + j) < input.size()) && !isupper(input[i + j])) && !(clean && (input[i + j] == '_' || input[i + j] == ':' || isdigit(input[i + j])))) {
        aux[j] = input[i + j];
        j++;
      }
      aux[j] = '\0';
      string auxstr(aux);
      if(clean) { auxstr = KBIN::VASP_PseudoPotential_CleanName(auxstr); }
      velements[counter] = auxstr;
      i += j - 1;  //CO add j to i so we don't check these characters again
      counter++;
    }
  }

  //////////////////////////////////////////////////////////////////////////////
  // END Parsing input
  //////////////////////////////////////////////////////////////////////////////

  return velements;
}
}  // namespace pflow

namespace pflow {
// ***************************************************************************
// pflow::makeAlphabeticVector(string input,ostream& oss,ofstream&
// FileMESSAGE)
// ***************************************************************************
// PdMn -> MnPd, does it by CAPITAL letters
vector<string> makeAlphabeticVector(string input,
                                    ostream& oss) {  // overload
  ofstream FileMESSAGE;
  return makeAlphabeticVector(input, FileMESSAGE, oss);
}
vector<string> makeAlphabeticVector(string input,
                                    ofstream& FileMESSAGE, ostream& oss) {  // main function
  string soliloquy = "pflow::makeAlphabeticVector():";
  vector<string> velements =
      stringElements2VectorElements(input, FileMESSAGE, oss);
  sort(velements.begin(),
       velements.end());  // quicksort is much faster than insertion sort
  return velements;
}
// ***************************************************************************
// pflow::makeAlphabeticString(string input,ostream& oss,ofstream&
// FileMESSAGE)
// ***************************************************************************
// PdMn -> MnPd, does it by CAPITAL letters
string makeAlphabeticString(string input, ostream& oss) {  // overload
  ofstream FileMESSAGE;
  return makeAlphabeticString(input, FileMESSAGE, oss);
}
string makeAlphabeticString(string input,
                            ofstream& FileMESSAGE, ostream& oss) {  // main function
  string soliloquy = "pflow::makeAlphabeticString():";
  string corrected = "";
  vector<string> velements =
      stringElements2VectorElements(input, FileMESSAGE, oss);
  if(!velements.size()) {
    return corrected;
  }
  sort(velements.begin(),
       velements.end());  // quicksort is much faster than insertion sort
  for (uint i = 0; i < velements.size(); i++) {
    corrected += velements.at(i);
  }
  return corrected;
}
}  // namespace pflow

namespace pflow {
//***************************************************************************//
// pflow::defaultLoadEntriesFlags(const string& input)
//***************************************************************************//
// sets some default flags
void defaultLoadEntriesFlags(aurostd::xoption& vpflow, ostream& oss, string input, bool entry_output, bool silent) {  // overload
  ofstream FileMESSAGE;
  return defaultLoadEntriesFlags(vpflow, FileMESSAGE, oss, input, entry_output, silent);
}
void defaultLoadEntriesFlags(aurostd::xoption& vpflow,ofstream& FileMESSAGE, ostream& oss, string input, bool entry_output, bool silent) {  // main function
  string soliloquy = "pflow::defaultLoadEntriesFlags()";
  stringstream message;

  // common
  if(entry_output) {
    vpflow.flag("PFLOW::LOAD_ENTRIES_ENTRY_OUTPUT", true);
    if(!silent) {
      pflow::logger(soliloquy, "PFLOW::LOAD_ENTRIES_ENTRY_OUTPUT set to TRUE", FileMESSAGE, oss, _LOGGER_OPTION_);
    }
  }
  vpflow.flag("PFLOW::LOAD_ENTRIES_ONLY_ALPHABETICAL",
              true);  // un-alphabetized entries are crap
  if(!silent) {
    pflow::logger(soliloquy, "PFLOW::LOAD_ENTRIES_ONLY_ALPHABETICAL set to TRUE", FileMESSAGE, oss, _LOGGER_OPTION_);
  }
  vpflow.flag("PFLOW::LOAD_ENTRIES_NARIES_MINUS_ONE", true);  //if loading ternary, also load relevant binaries and unaries
  if(!silent) {
    pflow::logger(soliloquy, "PFLOW::PFLOW::LOAD_ENTRIES_NARIES_MINUS_ONE set to TRUE", FileMESSAGE, oss, _LOGGER_OPTION_);
  }
  input = aurostd::toupper(input);  //removes ALL ambiguity with case
  //must be specific LIB1, LIB2, etc.
  if(!(input == "A" || input == "ICSD")) { //|| input == "icsd")) {  //put other special libs here
    string lib_name;
    //make robust so input can be "LIB2" or "2"
    if(aurostd::substring2bool(input, "LIB")) {
      lib_name = input;
    } else {
      lib_name = "LIB" + input;
    }
    vpflow.flag("PFLOW::LOAD_ENTRIES_LOAD_" + std::string(lib_name), true);
    if(!silent) {
      pflow::logger(soliloquy, "PFLOW::LOAD_ENTRIES_LOAD_" + std::string(lib_name) + " set to TRUE", FileMESSAGE, oss, _LOGGER_OPTION_);
    }
  } else {
    if(input == "A") {
      //get appropriate size, slightly inefficient (as we got this before), but it's cheap
      //string system = vpflow.getattachedscheme("PFLOW::ALLOY");  //CO 170908 - don't want to have to set this everytime
      //vector<string> velements = pflow::stringElements2VectorElements(system, oss, FileMESSAGE);  //un-sorted, okay
      string lib_count_string;
      //for (uint i = 0; i < velements.size() && i <= _AFLOW_LIB_MAX_; i++) { //CO 170908 - simply load all, LoadEntries() limits appropriately by velements.size()
      for (uint i = 0; i <= _AFLOW_LIB_MAX_; i++) {
        if(i == 0) { continue; }  //skip LIB1 by default
        lib_count_string = aurostd::utype2string(i + 1);
        vpflow.flag("PFLOW::LOAD_ENTRIES_LOAD_LIB" + lib_count_string, true);
        if(!silent) {
          pflow::logger(soliloquy, "PFLOW::LOAD_ENTRIES_LOAD_LIB" + lib_count_string + " set to TRUE", FileMESSAGE, oss, _LOGGER_OPTION_);
        }
      }
    }
    if(input == "A" || input == "ICSD") { //|| input == "icsd") {
      vpflow.flag("PFLOW::LOAD_ENTRIES_LOAD_ICSD", true);
      if(!silent) {
        pflow::logger(soliloquy, "PFLOW::LOAD_ENTRIES_LOAD_ICSD set to TRUE", FileMESSAGE, oss, _LOGGER_OPTION_);
      }
    }
  }
}
}

// ***************************************************************************
// pflow::prototypeMatch(string p1,string p2)
// ***************************************************************************
// better than an exact match of proto_database and proto_search
// we match for 549 and 549.bis
// proto_database is "549.bis" and proto_search is "549"
// only knows "." and ":" for now
namespace pflow {
bool prototypeMatch(string proto_database, string proto_search){
  if(proto_database==proto_search){return TRUE;}
  if(proto_database.size()<(proto_search.size()+1)){return FALSE;}  //we only match "549" + something
  if(proto_database.substr(0,proto_search.size())!=proto_search){return FALSE;} //not something + "549"
  if(proto_database[proto_search.size()]=='.' || proto_database[proto_search.size()]==':'){return TRUE;}
  //cerr << "What is this: " << proto_database << endl;
  return FALSE;
}
}  // namespace pflow
//Added by Corey Oses - May 2017
//END - all relevent functions for loading entries here
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace pflow {
//***************************************************************************//
// pflow::logger(const char& type,const string& function_name,string
// _message,bool silent,ostream& oss,ofstream& FileMESSAGE)
//***************************************************************************//
// added by Corey Oses - May 2017
// effectively logs EVERYTHING, deals with colors, oss, and log files (very robust)
// types include error, warning, options, and more as needed (see below)
// raw is special, prints raw message
// message can be string or stringstream (if stringstream, it gets cleared out)
// oss = cout, cerr (as you prefer)
// FileMESSAGE - logfile
void logger(const string& function_name, stringstream& message, ostream& oss, const char& type, bool silent, const string& message_metadata) {  // overload
  string _message = message.str();
  logger(function_name, _message, oss, type, silent, message_metadata);
  message.str("");
}

void logger(const string& function_name, stringstream& message, ofstream& FileMESSAGE, ostream& oss, const char& type, bool silent, const string& message_metadata) {  // overload
  string _message = message.str();
  logger(function_name, _message, FileMESSAGE, oss, type, silent, message_metadata);
  message.str("");
}

void logger(const string& function_name, stringstream& message, const _aflags& aflags, ostream& oss, const char& type, bool silent, const string& message_metadata) {  // overload
  string _message = message.str();
  logger(function_name, _message, aflags, oss, type, silent, message_metadata);
  message.str("");
}

void logger(const string& function_name, stringstream& message, const _aflags& aflags, ofstream& FileMESSAGE, ostream& oss, const char& type, bool silent, const string& message_metadata) {  // overload
  string _message = message.str();
  logger(function_name, _message, aflags, FileMESSAGE, oss, type, silent, message_metadata);
  message.str("");
}

void logger(const string& function_name, const string& _message, ostream& oss, const char& type, bool silent, const string& message_metadata) {  // overload
  ofstream FileMESSAGE;
  logger(function_name, _message, FileMESSAGE, oss, type, silent, message_metadata);
}

void logger(const string& function_name, const string& _message, ofstream& FileMESSAGE, ostream& oss, const char& type, bool silent, const string& message_metadata) {  // overload
  _aflags aflags; aflags.Directory=aurostd::execute2string(XHOST.command("pwd"))+"/"; //".";  //CO 180220
  logger(function_name, _message, aflags, FileMESSAGE, oss, type, silent, message_metadata);
}

void logger(const string& function_name, const string& _message, const _aflags& aflags, ostream& oss, const char& type, bool silent, const string& message_metadata) {  // overload
  ofstream FileMESSAGE;
  logger(function_name, _message, aflags, FileMESSAGE, oss, type, silent, message_metadata);
}

void logger(const string& function_name, const string& _message, const _aflags& aflags, ofstream& FileMESSAGE, ostream& oss, const char& type, bool silent, const string& message_metadata) {  // main function
  // five types:  M - Message, W - Warning, E - Error, O - Option, C - Complete, R - RAW
  // treat E, W, O, C, and R as special, otherwise treat as a message
  // no need for function name for R, you can put ""

  string message = aurostd::RemoveWhiteSpacesFromTheBack(_message);
  if(message.empty()) {  // && oss.str().empty()) {
    return;
  }

  string soliloquy = aurostd::RemoveWhiteSpaces(function_name);
  string ErrorBarString =
      "EEEEE  "
      "------------------------------------------------------------------------"
      "---------------------------------------------------- ";
  string WarningBarString =
      "WWWWW  "
      "------------------------------------------------------------------------"
      "---------------------------------------------------- ";

  if (type == _LOGGER_ERROR_) {
    ////////////////////////////////////////////////////////////////////////////
    // START Error logger
    ////////////////////////////////////////////////////////////////////////////

    // write to screen if not quiet
    if(!XHOST.QUIET && !silent) {
      // borrowed from APL/aflow_apl.h
      printf("\033[31m");     // red
      oss << ErrorBarString;  // make it clear in log file that an error
      // occurred
      oss << endl;
      if(!message.empty()) {
        oss << "EEEEE";
        oss << "  ";

        printf("\033[0m");   // turn off all cursor attributes
        printf("\033[31m");  // red
        printf("\033[5m");   // bold/blink
        oss << "ERROR";
        printf("\033[0m");   // turn off all cursor attributes
        printf("\033[31m");  // red

        oss << " ";
        oss << soliloquy;
        oss << " ";
        oss << message;
        oss << Message(aflags, message_metadata);
        oss << endl;
      }
      oss << ErrorBarString;  // make it clear in log file that an error
      // occurred
      oss << endl;
      printf("\033[0m");  // turn off all cursor attributes
      oss.flush();
    }
    // write to log
    FileMESSAGE << ErrorBarString;  // make it clear in log file that an error occurred
    FileMESSAGE << endl;
    if(!message.empty()) {
      FileMESSAGE << "EEEEE";
      FileMESSAGE << "  ";
      FileMESSAGE << "ERROR";
      FileMESSAGE << " ";
      FileMESSAGE << soliloquy;
      FileMESSAGE << " ";
      FileMESSAGE << message;
      FileMESSAGE << Message(aflags, message_metadata);
      FileMESSAGE << endl;
    }
    FileMESSAGE << ErrorBarString;  // make it clear in log file that an error occurred
    FileMESSAGE << endl;
    FileMESSAGE.flush();

    ////////////////////////////////////////////////////////////////////////////
    // END Error logger
    ////////////////////////////////////////////////////////////////////////////

  } else if (type == _LOGGER_WARNING_) {
    ////////////////////////////////////////////////////////////////////////////
    // START Warning logger
    ////////////////////////////////////////////////////////////////////////////

    if(!XHOST.QUIET && !silent) {
      // borrowed from APL/aflow_apl.h
      printf("\033[33m\033[1m");  // yellow
      oss << WarningBarString;    // make it clear in log file that an warning
      // occurred
      oss << endl;
      if(!message.empty()) {
        oss << "WWWWW";
        oss << "  ";

        printf("\033[0m");          // turn off all cursor attributes
        printf("\033[33m\033[1m");  // yellow
        printf("\033[5m");          // bold/blink
        oss << "WARNING";
        printf("\033[0m");          // turn off all cursor attributes
        printf("\033[33m\033[1m");  // yellow

        oss << " ";
        oss << soliloquy;
        oss << " ";
        oss << message;
        oss << Message(aflags, message_metadata);
        oss << endl;
      }
      oss << WarningBarString;  // make it clear in log file that an warning
      // occurred
      oss << endl;
      printf("\033[0m");  // turn off all cursor attributes
      oss.flush();
    }
    // write to log
    FileMESSAGE << WarningBarString;  // make it clear in log file that an
    // warning occurred
    FileMESSAGE << endl;
    if(!message.empty()) {
      FileMESSAGE << "WWWWW";
      FileMESSAGE << "  ";
      FileMESSAGE << "WARNING";
      FileMESSAGE << " ";
      FileMESSAGE << soliloquy;
      FileMESSAGE << " ";
      FileMESSAGE << message;
      FileMESSAGE << Message(aflags, message_metadata);
      FileMESSAGE << endl;
    }
    FileMESSAGE << WarningBarString;  // make it clear in log file that a
    // warning occurred
    FileMESSAGE << endl;
    FileMESSAGE.flush();

    ////////////////////////////////////////////////////////////////////////////
    // END Warning logger
    ////////////////////////////////////////////////////////////////////////////

  } else if (type == _LOGGER_COMPLETE_) {
    ////////////////////////////////////////////////////////////////////////////
    // START Complete logger
    ////////////////////////////////////////////////////////////////////////////

    if(!XHOST.QUIET && !silent) {
      // borrowed from APL/aflow_apl.h
      printf("\033[32m");  // green
      if(!message.empty()) {
        oss << "CCCCC";
        oss << "  ";

        printf("\033[0m");   // turn off all cursor attributes
        printf("\033[32m");  // green
        printf("\033[5m");   // bold/blink
        oss << "COMPLETE";
        printf("\033[0m");   // turn off all cursor attributes
        printf("\033[32m");  // green

        oss << " ";
        oss << soliloquy;
        oss << " ";
        oss << message;
        oss << Message(aflags, message_metadata);
        oss << endl;
      }
      printf("\033[0m");  // turn off all cursor attributes
      oss.flush();
    }
    // write to log
    if(!message.empty()) {
      FileMESSAGE << "CCCCC";
      FileMESSAGE << "  ";
      FileMESSAGE << "COMPLETE";
      FileMESSAGE << " ";
      FileMESSAGE << soliloquy;
      FileMESSAGE << " ";
      FileMESSAGE << message;
      FileMESSAGE << Message(aflags, message_metadata);
      FileMESSAGE << endl;
    }
    FileMESSAGE.flush();

    ////////////////////////////////////////////////////////////////////////////
    // END Complete logger
    ////////////////////////////////////////////////////////////////////////////

  } else if (type == _LOGGER_OPTION_) {
    ////////////////////////////////////////////////////////////////////////////
    // START Option logger
    ////////////////////////////////////////////////////////////////////////////

    if(!XHOST.QUIET && !silent) {
      // borrowed from APL/aflow_apl.h
      if(!message.empty()) {
        oss << "-OPT-";
        oss << "  ";
        oss << "MESSAGE-OPTION";
        oss << " ";
        oss << soliloquy;
        oss << " ";
        oss << message;
        oss << Message(aflags, message_metadata);
        oss << endl;
      }
      oss.flush();
    }
    // write to log
    if(!message.empty()) {
      FileMESSAGE << "-OPT-";
      FileMESSAGE << "  ";
      FileMESSAGE << "MESSAGE-OPTION";
      FileMESSAGE << " ";
      FileMESSAGE << soliloquy;
      FileMESSAGE << " ";
      FileMESSAGE << message;
      FileMESSAGE << Message(aflags, message_metadata);
      FileMESSAGE << endl;
    }
    FileMESSAGE.flush();

    ////////////////////////////////////////////////////////////////////////////
    // END Option logger
    ////////////////////////////////////////////////////////////////////////////

  } else if (type == _LOGGER_RAW_) {
    ////////////////////////////////////////////////////////////////////////////
    // START Option logger
    ////////////////////////////////////////////////////////////////////////////

    if (!XHOST.QUIET && !silent) {
      // borrowed from APL/aflow_apl.h
      if (!message.empty()) {
        oss << message;
      }
      oss.flush();
    }
    // write to log
    if (!message.empty()) {
      FileMESSAGE << message;
    }
    FileMESSAGE.flush();

    ////////////////////////////////////////////////////////////////////////////
    // END Raw logger
    ////////////////////////////////////////////////////////////////////////////

  } else {
    ////////////////////////////////////////////////////////////////////////////
    // START Message logger
    ////////////////////////////////////////////////////////////////////////////

    if(!XHOST.QUIET && !silent) {
      // borrowed from APL/aflow_apl.h
      if(!message.empty()) {
        oss << "00000";
        oss << "  ";
        oss << "MESSAGE";
        oss << " ";
        oss << soliloquy;
        oss << " ";
        oss << message;
        oss << Message(aflags, message_metadata);
        oss << endl;
      }
      oss.flush();
    }
    // write to log
    if(!message.empty()) {
      FileMESSAGE << "00000";
      FileMESSAGE << "  ";
      FileMESSAGE << "MESSAGE";
      FileMESSAGE << " ";
      FileMESSAGE << soliloquy;
      FileMESSAGE << " ";
      FileMESSAGE << message;
      FileMESSAGE << Message(aflags, message_metadata);
      FileMESSAGE << endl;
    }
    FileMESSAGE.flush();

    ////////////////////////////////////////////////////////////////////////////
    // END Message logger
    ////////////////////////////////////////////////////////////////////////////
  }
  }
} // namespace pflow

// ***************************************************************************
// pflow::LTCELL
// ***************************************************************************
namespace pflow {
  xstructure LTCELL(string options,istream& input) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "pflow::LTCELL: BEGIN" << endl;
    xstructure str(input,IOAFLOW_AUTO);
    vector<string> tokens;
    aurostd::string2tokens(options,tokens,",");
    xmatrix<double> mlt(3,3);

    if(tokens.size()!=9 && tokens.size()!=3 && tokens.size()!=1 && tokens.size()!=4) {
      init::ErrorOption(cout,options,"pflow::LTCELL",
			aurostd::liststring2string("aflow --ltcell=a11,a12,a13,a21,a22,a23,a31,a32,a33 < POSCAR",
                                                   "aflow --ltcell=a11,a22,a33 < POSCAR",
                                                   "aflow --ltcell=file < POSCAR",
                                                   "aflow --ltcellfv=v1,v2,v3,phi < POSCAR"));
      exit(0);
    }

    if(tokens.size()==9) {
      if(LDEBUG) cerr << "pflow::LTCELL: 9 entries" << endl;
      mlt(1,1)=aurostd::string2utype<double>(tokens.at(0)); 
      mlt(1,2)=aurostd::string2utype<double>(tokens.at(1)); 
      mlt(1,3)=aurostd::string2utype<double>(tokens.at(2)); 
      mlt(2,1)=aurostd::string2utype<double>(tokens.at(3)); 
      mlt(2,2)=aurostd::string2utype<double>(tokens.at(4)); 
      mlt(2,3)=aurostd::string2utype<double>(tokens.at(5)); 
      mlt(3,1)=aurostd::string2utype<double>(tokens.at(6)); 
      mlt(3,2)=aurostd::string2utype<double>(tokens.at(7)); 
      mlt(3,3)=aurostd::string2utype<double>(tokens.at(8)); 
      if(abs(det(mlt))<0.01) {cerr << "ERROR - pflow::LTCELL: singular ltcell matrix" << endl;exit(0);}
      if(LDEBUG) cerr << "pflow::LTCELL: END" << endl;
      return GetLTCell(mlt,str);
    }

    if(tokens.size()==3) {
      if(LDEBUG) cerr << "pflow::LTCELL: 3 entries" << endl;
      mlt(1,1)=aurostd::string2utype<double>(tokens.at(0)); 
      mlt(2,2)=aurostd::string2utype<double>(tokens.at(1)); 
      mlt(3,3)=aurostd::string2utype<double>(tokens.at(2)); 
      if(abs(det(mlt))<0.01) {cerr << "ERROR - pflow::LTCELL: singular ltcell matrix" << endl;exit(0);}
      if(LDEBUG) cerr << "pflow::LTCELL: END" << endl;
      return GetLTCell(mlt,str);
    }

    if(tokens.size()==1) {
      if(LDEBUG) cerr << "pflow::LTCELL: 1 entries" << endl;
      ifstream infile(tokens.at(0).c_str());
      aurostd::InFileExistCheck("pflow::LTCELL",tokens.at(0).c_str(),infile,cerr);
      for(int i=1;i<=3;i++)
	for(int j=1;j<=3;j++)
	  infile >> mlt(i,j);
      if(abs(det(mlt))<0.01) {cerr << "ERROR - pflow::LTCELL: singular ltcell matrix" << endl;exit(0);}
      if(LDEBUG) cerr << "pflow::LTCELL: END" << endl;
      return GetLTCell(mlt,str);
    }
    
    if(tokens.size()==4) {
      if(LDEBUG) cerr << "pflow::LTCELL: 4 entries => LTCELLFV mode" << endl;
      xvector<double> nvec(3);
      nvec(1)=aurostd::string2utype<double>(tokens.at(0));
      nvec(2)=aurostd::string2utype<double>(tokens.at(1));
      nvec(3)=aurostd::string2utype<double>(tokens.at(2));
      double angle=aurostd::string2utype<double>(tokens.at(3))/rad2deg;
      if(LDEBUG) cerr << "pflow::LTCELL: END" << endl;
      return GetLTFVCell(nvec,angle,str);
    }
    return str;
  }
} // namespace pflow

// ***************************************************************************
// pflow::LTCELLFV
// ***************************************************************************
// [OBSOLETE] namespace pflow {
// [OBSOLETE]   xstructure LTCELLFV(string options,istream& input) {
// [OBSOLETE]     //  if(argv.size()!=a.num_each_type.size()+2) {
// [OBSOLETE]     //  cerr << "ERROR - pflow::LTCELLFV: you need to specify as many names as atom types" << endl;
// [OBSOLETE]     //   exit(0);
// [OBSOLETE]     // }
// [OBSOLETE]     // Read in input file.
// [OBSOLETE]     bool LDEBUG=1;//(FALSE || XHOST.DEBUG);
// [OBSOLETE]     if(LDEBUG) cerr << "pflow::LTCELLFV: BEGIN" << endl;
// [OBSOLETE]     xstructure str(input,IOAFLOW_AUTO);
// [OBSOLETE]     vector<string> tokens;
// [OBSOLETE]     aurostd::string2tokens(options,tokens,",");
// [OBSOLETE]  
// [OBSOLETE]    if(tokens.size()!=4) {
// [OBSOLETE]       init::ErrorOption(cout,options,"pflow::LTCELLFV","aflow --ltcellfv=v1,v2,v3,phi < POSCAR");
// [OBSOLETE]       exit(0);
// [OBSOLETE]     }
// [OBSOLETE] 
// [OBSOLETE]     xvector<double> nvec(3);
// [OBSOLETE]     nvec(1)=atof(tokens.at(0).c_str());
// [OBSOLETE]     nvec(2)=atof(tokens.at(1).c_str());
// [OBSOLETE]     nvec(3)=atof(tokens.at(2).c_str());
// [OBSOLETE]     double angle=atof(tokens.at(3).c_str())/rad2deg;
// [OBSOLETE]     return GetLTFVCell(nvec,angle,str);
// [OBSOLETE]   }
// [OBSOLETE] } // namespace pflow

// ***************************************************************************
// pflow::MagneticParameters
// ***************************************************************************
namespace pflow {
  void MagneticParameters(string _directory, ostream& oss) {
    string directory(_directory);
    if(directory=="") directory="./";
    directory=aurostd::CleanFileName(directory);
 
    bool found=FALSE;
    vector<string> vfiles;
    vfiles.clear();
    vfiles.push_back(aurostd::CleanFileName(directory+"/OUTCAR.static.bz2"));
    vfiles.push_back(aurostd::CleanFileName(directory+"/OUTCAR.static"));
    vfiles.push_back(aurostd::CleanFileName(directory+"/OUTCAR"));
    vfiles.push_back(aurostd::CleanFileName(directory+"/OUTCAR.relax2.bz2"));
    vfiles.push_back(aurostd::CleanFileName(directory+"/OUTCAR.relax2"));
    vfiles.push_back(aurostd::CleanFileName(directory+"/OUTCAR.relax1.bz2"));
    vfiles.push_back(aurostd::CleanFileName(directory+"/OUTCAR.relax1"));

    xOUTCAR outcar;
    for(uint i=0;i<vfiles.size()&&!found;i++) {
      if(aurostd::FileExist(vfiles.at(i))) {
	found=TRUE;
	outcar.GetPropertiesFile(vfiles.at(i));
	oss << "pflow::MagneticParameters  : using file=" << vfiles.at(i) << endl;
	
	oss << "MAGNETIC MOMENTUM CELL : " << outcar.mag_cell << endl;
	oss << "MAGNETIC MOMENTUM ATOM : " << outcar.mag_atom << endl;
	oss << "VOLUME CELL            : " << outcar.volume_cell << endl;
	oss << "VOLUME ATOM            : " << outcar.volume_atom << endl;
	if(outcar.vmag.size()>1) {
	  oss << "SPIN DECOMPOSITION     : " << outcar.vmag.at(0);
	  for(uint i=1;i<outcar.vmag.size();i++) oss << "," << outcar.vmag.at(i);
	  oss << endl;
	}
      }
    }

    found=FALSE;
    vfiles.clear();
    vfiles.push_back(aurostd::CleanFileName(directory+"/DOSCAR.static.bz2"));
    vfiles.push_back(aurostd::CleanFileName(directory+"/DOSCAR.static"));
    vfiles.push_back(aurostd::CleanFileName(directory+"/DOSCAR"));

    xDOSCAR doscar;
    for(uint i=0;i<vfiles.size()&&!found;i++) {
      if(aurostd::FileExist(vfiles.at(i))) {
	found=TRUE;
	doscar.GetPropertiesFile(vfiles.at(i));
	oss << "pflow::MagneticParameters  : using file=" << vfiles.at(i) << endl;
	oss << "POLARIZATION FERMI     : " << doscar.spinF << endl;
      }
    }
  }
} // namespace pflow

// ***************************************************************************
// pflow::MAKESTRLIST
// ***************************************************************************
namespace pflow {
  void MAKESTRLIST(vector<string> argv) {
    ifstream outcar_inf(argv.at(2).c_str());
    aurostd::InFileExistCheck("convasp",argv.at(2).c_str(),outcar_inf,cerr);
    ifstream xdatcar_inf(argv.at(3).c_str());
    aurostd::InFileExistCheck("convasp",argv.at(3).c_str(),xdatcar_inf,cerr);
    vector<xstructure> str_vec=pflow::GetStrVecFromOUTCAR_XDATCAR(outcar_inf,xdatcar_inf);
    pflow::PrintStrVec(str_vec,cout);
  }
} // namespace pflow

// [OBSOLETE] // ***************************************************************************
// [OBSOLETE] // pflow::MILLER
// [OBSOLETE] // ***************************************************************************
// [OBSOLETE] namespace pflow {
// [OBSOLETE]   xstructure MILLER(string options,istream& input) {
// [OBSOLETE]     bool LDEBUG=(FALSE || XHOST.DEBUG);
// [OBSOLETE]     if(LDEBUG) cerr << "pflow::MILLER: BEGIN" << endl;  
// [OBSOLETE]     vector<string> tokens;
// [OBSOLETE]     aurostd::string2tokens(options,tokens,",");
// [OBSOLETE]     if(tokens.size()!=3 && tokens.size()!=4 && tokens.size()!=5) {
// [OBSOLETE]       init::ErrorOption(cout,options,"pflow::MILLER","aflow -miller=h,k,l[,nlayer[,elayer] < POSCAR");
// [OBSOLETE]       exit(0);
// [OBSOLETE]     }
// [OBSOLETE] 
// [OBSOLETE]     xstructure a(input,IOAFLOW_AUTO); // LOAD fixes all lattices
// [OBSOLETE]     double minradius=0.0,epsilon=0.01;
// [OBSOLETE]     a.ReScale(1.0);
// [OBSOLETE]     a.FixLattices(); // redundant
// [OBSOLETE]     xvector<int> dims(3);
// [OBSOLETE]     xvector<double> hkl(3),n(3),r(3),r2(3),s1(3),s2(3),s3(3),spoint(3);
// [OBSOLETE]     int nlayers=1,elayers=0;
// [OBSOLETE]     hkl(1)=0.0;hkl(2)=0.0;hkl(3)=0.0;
// [OBSOLETE]     if(tokens.size()>=1) hkl(1)=aurostd::string2utype<double>(tokens.at(0));
// [OBSOLETE]     if(tokens.size()>=2) hkl(2)=aurostd::string2utype<double>(tokens.at(1));
// [OBSOLETE]     if(tokens.size()>=3) hkl(3)=aurostd::string2utype<double>(tokens.at(2));
// [OBSOLETE]     if(tokens.size()>=4) nlayers=aurostd::string2utype<int>(tokens.at(3));
// [OBSOLETE]     if(tokens.size()>=5) elayers=aurostd::string2utype<int>(tokens.at(4));
// [OBSOLETE] 
// [OBSOLETE]     if(LDEBUG) cerr << "DEBUG nlayers=" << nlayers << endl;
// [OBSOLETE]     if(LDEBUG) cerr << "DEBUG elayers=" << elayers << endl;
// [OBSOLETE]     xmatrix<double> plattice(3,3),slattice(3,3);
// [OBSOLETE]     n=hkl*trasp(a.klattice)/(2*pi); // divided by 2pi
// [OBSOLETE]     // double h=hkl[1],k=hkl[2],l=hkl[3];
// [OBSOLETE]     // n=n/modulus(n); // normal
// [OBSOLETE]     dims=LatticeDimensionSphere(a.lattice,modulus(n)*1.1);
// [OBSOLETE]     vector<xvector<double> > nvectors,pvectors;
// [OBSOLETE]     for(int i=dims[1];i>=-dims[1];i--)
// [OBSOLETE]       for(int j=dims[2];j>=-dims[2];j--)
// [OBSOLETE] 	for(int k=dims[3];k>=-dims[3];k--) {
// [OBSOLETE] 	  r=((double)i)*a.lattice(1)+((double)j)*a.lattice(2)+((double)k)*a.lattice(3);
// [OBSOLETE] 	  if(modulus(r)>epsilon && abs(scalar_product(r,n))<epsilon) nvectors.push_back(r);
// [OBSOLETE] 	  if(modulus(r)>epsilon && abs(angle(r,n))<epsilon) pvectors.push_back(r);
// [OBSOLETE] 	}
// [OBSOLETE]     if(LDEBUG) cerr << "DEBUG nvectors.size()=" << nvectors.size() << endl;
// [OBSOLETE]     if(LDEBUG) cerr << "DEBUG pvectors.size()=" << pvectors.size() << endl;
// [OBSOLETE] 
// [OBSOLETE]     // choose shortest pvectors
// [OBSOLETE]     s3=pvectors.at(0);
// [OBSOLETE]     for(uint i=0;i<pvectors.size();i++) if(modulus(pvectors.at(i))<modulus(s3)) s3=pvectors.at(i);
// [OBSOLETE]     if(LDEBUG) cerr << "DEBUG s3=" << s3 << endl;
// [OBSOLETE]     if(LDEBUG) cerr << "DEBUG n=" << n << endl;
// [OBSOLETE]     // test couples in plane to get the in-plane basis, and optimize with respect to minkowski
// [OBSOLETE]     minradius=modulus(nvectors.at(0))+modulus(nvectors.at(1))+modulus(n); // some default to break
// [OBSOLETE]     for(int ii=1;ii<=3;ii++) plattice[ii][3]=s3[ii];
// [OBSOLETE]     for(uint i=0;i<nvectors.size();i++) {
// [OBSOLETE]       for(int ii=1;ii<=3;ii++) plattice[ii][1]=nvectors.at(i)[ii];
// [OBSOLETE]       for(uint j=i+1;j<nvectors.size();j++) {
// [OBSOLETE] 	for(int ii=1;ii<=3;ii++) plattice[ii][2]=nvectors.at(j)[ii];
// [OBSOLETE] 	if(det(plattice)>epsilon) { // check only for well defined lattices
// [OBSOLETE] 	  bool basisflag=TRUE;
// [OBSOLETE] 	  for(uint k=0;k<nvectors.size() && basisflag==TRUE;k++)
// [OBSOLETE] 	    if(k!=i && k!=j && i!=j) {
// [OBSOLETE] 	      spoint=inverse(plattice)*nvectors.at(k);  // get the spoint in fractional coordinates
// [OBSOLETE] 	      if(!aurostd::isinteger(spoint[1]) || !aurostd::isinteger(spoint[2]) || !aurostd::isinteger(spoint[3])) basisflag=FALSE;
// [OBSOLETE] 	    }
// [OBSOLETE] 	  if(basisflag==TRUE && RadiusSphereLattice(trasp(plattice))<minradius) {
// [OBSOLETE] 	    minradius=RadiusSphereLattice(trasp(plattice));
// [OBSOLETE] 	    slattice=trasp(plattice);
// [OBSOLETE] 	  }
// [OBSOLETE] 	}
// [OBSOLETE]       }
// [OBSOLETE]     }
// [OBSOLETE]     // this is a good lattice now s1,s2 are planar s3 is in the orthogonal direction
// [OBSOLETE]     s1=slattice(1);s2=slattice(2);s3=slattice(3);
// [OBSOLETE]     // DEBUG=TRUE;
// [OBSOLETE]     if(LDEBUG) cerr << "DEBUG s1=  " << s1 << endl;
// [OBSOLETE]     if(LDEBUG) cerr << "DEBUG s2=  " << s2 << endl;
// [OBSOLETE]     if(LDEBUG) cerr << "DEBUG s3=  " << s3 << endl;
// [OBSOLETE]     if(LDEBUG) cerr << "DEBUG  n=  " << n << endl;
// [OBSOLETE]     a.lattice=slattice;
// [OBSOLETE]     cerr << det(a.lattice) << " " << det(slattice) << endl;
// [OBSOLETE]     if(LDEBUG) cerr << "pflow::MILLER: END" << endl;  
// [OBSOLETE]     return a;
// [OBSOLETE]   }
// [OBSOLETE] } // namespace pflow

// ***************************************************************************
// pflow::MINKOWSKIBASISREDUCTION
// ***************************************************************************
namespace pflow {
  xstructure MINKOWSKIBASISREDUCTION(istream& input) {
    xstructure a(input,IOAFLOW_AUTO);
    xstructure b(a);
    b.MinkowskiBasisReduction();
    //b.BringInCell(); // not necessary
    return b;
  }
} // namespace pflow

// ***************************************************************************
// pflow::MISCIBILITY
// ***************************************************************************
namespace pflow {
  string MISCIBILITY(vector<string> argv) {
    stringstream output;
    for(uint i=2;i<argv.size();i++) {
      string system=string(argv.at(i));
      XATOM_AlphabetizationSpecies(system);
      //  XATOM_AlphabetizationCompound(system);
      int mix=MiscibilityCheck(system);
      if(mix==MISCIBILITY_SYSTEM_NOMIX) output << system << " " << "MISCIBILITY_SYSTEM_NOMIX" << endl;
      if(mix==MISCIBILITY_SYSTEM_MISCIBLE) output << system << " " << "MISCIBILITY_SYSTEM_MISCIBLE" << endl;
      if(mix==MISCIBILITY_SYSTEM_UNKNOWN) output << system << " " << "MISCIBILITY_SYSTEM_UNKNOWN" << endl;
    }
    return output.str();
  };
}

// ***************************************************************************
// pflow::MOM
// ***************************************************************************
namespace pflow {
  void MOM(istream& input) {
    xstructure str(input,IOAFLOW_AUTO);
    xvector<double> m(3);m=GetMom1(str);
    cout.setf(std::ios::left,std::ios::adjustfield);
    cout.setf(std::ios::fixed,std::ios::floatfield);
    cout << " Moment 1 : "<<setw(15)<<setprecision(10)<<m(1)<<setw(15)<<setprecision(10)<<m(2)<<setw(15)<<setprecision(10)<<m(3)<<endl;
    m=m/str.scale;
    cout << " (unscaled) "<<setw(15)<<setprecision(10)<<m(1)<<setw(15)<<setprecision(10)<<m(2)<<setw(15)<<setprecision(10)<<m(3)<<endl;
  }
} // namespace pflow

// ***************************************************************************
// pflow::MSI
// ***************************************************************************
namespace pflow {
  void MSI(istream& input) {
    xstructure str(input,IOAFLOW_AUTO);
    PrintMSI(str,cout);
  }
} // namespace pflow

// ***************************************************************************
// pflow::NATOMS
// ***************************************************************************
namespace pflow {
  uint NATOMS(istream& input) {
    xstructure a(input,IOAFLOW_AUTO);
    return a.atoms.size();
  }
} // namespace pflow

// ***************************************************************************
// pflow::NBONDXX
// ***************************************************************************
namespace pflow {
  //CO 171025
  string NBONDXX(istream& input,bool aflowlib_legacy_format) {
    xstructure a(input,IOAFLOW_AUTO);
    return NBONDXX(a,aflowlib_legacy_format);
  }

  string NBONDXX(const xstructure& a,bool aflowlib_legacy_format) {
    vector<double> nbondxx=GetNBONDXX(a);
    
    //return aurostd::joinWDelimiter(aurostd::vecDouble2vecString(nbondxx,9),',');
    //print nicely
    
    //are there names in the atoms?
    bool atom_names=true;
    for(uint i=0;i<a.atoms.size()&&atom_names;i++){if(a.atoms[i].name.empty()){atom_names=false;}}

    //get names
    uint iat=0;
    vector<int> first_itypes;
    first_itypes.push_back(0);
    bool found;
    for(uint i=1;i<a.atoms.size();i++){
      found=false;
      for(uint j=0;j<first_itypes.size()&&!found;j++){
        if(a.atoms[i].type==a.atoms[first_itypes[j]].type){found=true;}
      }
      if(!found){
        first_itypes.push_back(i);
      }
    }

    if(aflowlib_legacy_format){atom_names=false;}
    
    vector<string> names;
    stringstream ss;
    if(atom_names){
      for(uint i=0;i<first_itypes.size();i++){
        names.push_back(a.atoms[first_itypes[i]].name);
      }
    }else{
      for(uint i=0;i<first_itypes.size();i++){
        ss.str("");
        ss << char('A'+iat++);
        names.push_back(ss.str());
      }
    }

    iat=0;
    stringstream output; output.str("");
    if(aflowlib_legacy_format){
      for(uint itype=0;itype<first_itypes.size();itype++){
        for(uint jtype=itype;jtype<first_itypes.size();jtype++){
          output << aurostd::PaddedPOST("BOND_"+names[itype]+names[jtype],11) << " " << aurostd::utype2string(nbondxx[iat++],6) << " [Angst]" << endl;
        }
      }
    }else{
      for(uint itype=0;itype<first_itypes.size();itype++){
        for(uint jtype=itype;jtype<first_itypes.size();jtype++){
          output << aurostd::PaddedPOST(names[itype]+"-"+names[jtype]+":",8) << " " << aurostd::utype2string(nbondxx[iat++],6) << " [Angst]" << endl;
        }
      }
    }
    
    return output.str();
  }
} // namespace pflow

// ***************************************************************************
// pflow::NNDIST
// ***************************************************************************
namespace pflow {
  double NNDIST(istream& input) {
    xstructure a(input,IOAFLOW_AUTO);
    return(NearestNeighbour(a));
  }
} // namespace pflow

// ***************************************************************************
// pflow::NSPECIES
// ***************************************************************************
namespace pflow {
  uint NSPECIES(istream& input) {
    xstructure a(input,IOAFLOW_AUTO);
    return a.num_each_type.size();
  }
} // namespace pflow

// ***************************************************************************
// pflow::NAMES
// ***************************************************************************
namespace pflow {
  xstructure NAMES(vector<string> argv,istream& input) {
    xstructure a(input,IOAFLOW_AUTO);
    if(argv.size()!= a.num_each_type.size()+2) {
      cerr << "ERROR - pflow::NAMES: you need to specify as many names as atom types" << endl;
      exit(0);
    }
    xstructure b=a;
    int iatom=0;
    for(uint itype=0;itype<a.num_each_type.size();itype++) {
      string species=string(argv.at(2+b.atoms.at(iatom).type));
      b.species.at(itype)=species;
      for(int j=0;j<a.num_each_type.at(itype);j++) {
	b.atoms.at(iatom).name=species;    // CONVASP_MODE
	b.atoms.at(iatom).CleanName();
	b.atoms.at(iatom).CleanSpin();
	b.atoms.at(iatom).name_is_given=TRUE;
	iatom++;
      }
    }
    return b;
  }
} // namespace pflow

// ***************************************************************************
// pflow::NANOPARTICLE
// ***************************************************************************
namespace pflow {
  xstructure NANOPARTICLE(istream& input,const xvector<double>& iparams) {
    bool LDEBUG=(FALSE || XHOST.DEBUG); //CO 180226
    string soliloquy="pflow::NANOPARTICLE)";
    //  cout << aflow::Banner("BANNER_TINY") << endl;
    double radius=NANOPARTICLE_RADIUS_DEFAULT;
    double distance=NANOPARTICLE_DISTANCE_DEFAULT;
    if(iparams.rows>=1) radius=iparams[1];
    if(iparams.rows>=2) distance=iparams[2];
    if(LDEBUG){ //CO 180226
      cerr << soliloquy << " radius=" << radius << endl;
      cerr << soliloquy << " distance=" << distance << endl;
    }
    xstructure a(input,IOAFLOW_AUTO),b;
    xvector<int> dims(3);
    xvector<double> fshift(3),cshift(3),ucell(3);
    a.ReScale(1.0);

    if(0) { // old origin
      a.BringInCell();
      //  a.SetCoordinates(_COORDS_CARTESIAN_);
      a.BringInCell();
      for(uint iat=0;iat<a.atoms.size();iat++) {                       // CENTER MASS
	cshift=cshift+a.atoms.at(iat).cpos/((double) a.atoms.size());  // CENTER MASS
	fshift=fshift+a.atoms.at(iat).fpos/((double) a.atoms.size());  // CENTER MASS
      }
      for(uint iat=0;iat<a.atoms.size();iat++) {             // SHIFT
	a.atoms.at(iat).fpos=a.atoms.at(iat).fpos-fshift;    // SHIFT
	a.atoms.at(iat).cpos=a.atoms.at(iat).cpos-cshift;    // SHIFT
      }
    }

    b=a;  // COPY
    while(b.atoms.size()>0) {b.RemoveAtom(0);} // EMPTY
    dims=LatticeDimensionSphere(a.lattice,2*radius+distance);
    if(LDEBUG){cerr << soliloquy << " dims=" << dims << endl;}  //CO 180226
    for(int i=1;i<=3;i++) {
      ucell(i)=ceil((2.0*radius+distance)/modulus(b.lattice(i)));
      if(LDEBUG){cerr << soliloquy << " ucell(i=" << i << ")=" << ucell(i) << endl;}  //CO 180226
      for(int j=1;j<=3;j++)
	b.lattice(i,j)=b.lattice(i,j)*ucell(i);
    }
    b.FixLattices();
    for(uint iat=0;iat<a.atoms.size();iat++) {             // SHIFT
      for(int i=-dims(1);i<=dims(1);i++) {
	for(int j=-dims(2);j<=dims(2);j++) {
	  for(int k=-dims(3);k<=dims(3);k++) {
	    _atom atom=a.atoms.at(iat);
	    atom.cpos=((double)i)*a.lattice(1)+((double)j)*a.lattice(2)+((double)k)*a.lattice(3)+a.atoms.at(iat).cpos;
	    atom.fpos=C2F(b.lattice,atom.cpos);               // put in fractional of new basis
	    if(modulus(atom.cpos)<=radius)  b.AddAtom(atom);
	  }
	}
      }
    }
    cerr << "atoms=" << b.atoms.size() << endl;
    //  b.SetCoordinates(_COORDS_FRACTIONAL_);
    b.SetCoordinates(_COORDS_CARTESIAN_);
    return b;
  }
} // namespace pflow

// ***************************************************************************
// pflow::NDATA
// ***************************************************************************
namespace pflow {
  void NDATA(istream& input) {
    xstructure a(input,IOAFLOW_AUTO);
    PrintNdata(a,cout);
  }
} // namespace pflow

// ***************************************************************************
// pflow::NIGGLI
// ***************************************************************************
namespace pflow {
  xstructure NIGGLI(istream& input) {
    xstructure a(input,IOAFLOW_AUTO);
    return GetNiggliStr(a);
  }
} // namespace pflow

// ***************************************************************************
// pflow::NOORDERPARAMETER
// ***************************************************************************
namespace pflow {
  xstructure NOORDERPARAMETER(istream& input) {
    xstructure a(input,IOAFLOW_AUTO);
    a.order_parameter_structure=FALSE;
    a.order_parameter_atoms.clear();
    for(uint i=0;i<a.atoms.size();i++) {
      a.atoms.at(i).order_parameter_atom=FALSE;
      a.atoms.at(i).order_parameter_value=0;
    }
    return a;
  }
} // namespace pflow

// ***************************************************************************
// pflow::NOSD
// ***************************************************************************
namespace pflow {
  xstructure NOSD(istream& input) {
    xstructure a(input,IOAFLOW_AUTO);
    xstructure b=a;
    // Read in input file.
    b.isd=FALSE;
    return b;
  }
} // namespace pflow

// ***************************************************************************
// pflow::NUMNAMES
// ***************************************************************************
namespace pflow {
  xstructure NUMNAMES(vector<string> argv,istream& input) {
    xstructure a(input,IOAFLOW_AUTO);
    if(argv.size()!=a.num_each_type.size()+2) {
      cerr << "ERROR - pflow::NUMNAMES: you need to specify as many names as atom types" << endl;
      exit(0);
    }
    xstructure b=a;
    int iatom=0;
    for(uint itype=0;itype<a.num_each_type.size();itype++) {
      string species=string(argv.at(2+b.atoms.at(iatom).type));
      b.species.at(itype)=species;
      for(int j=0;j<a.num_each_type.at(itype);j++) {
	ostringstream aus;
	aus << species << j+1;      // CONVASP_MODE
	b.atoms.at(iatom).name=aus.str();
	b.atoms.at(iatom).CleanName();
	b.atoms.at(iatom).CleanSpin();
	b.atoms.at(iatom).name_is_given=TRUE;
	iatom++;
      }
    }
    return b;
  }
} // namespace pflow

// ***************************************************************************
// pflow::PDB
// ***************************************************************************
namespace pflow {
  void PDB(istream& input) {
    xstructure a(input,IOAFLOW_AUTO);
    PrintPDB(a,cout);
  }
} // namespace pflow

// ***************************************************************************
// pflow::PDOS
// ***************************************************************************
namespace pflow {
  void PDOS(vector<string> argv) {
    cerr << "# WARNING: THIS REQUIRES AN ALTERED VERSION OF VASP - SEE aflow --help" << endl;
    int only_occ=0; // Use unoccupied states in PDOS.
    pflow::projdata prd;
    pflow::pdosdata pdd;
    prd.PROOUTinfile=argv.at(3);
    pdd.PDOSinfile=argv.at(2);
    pflow::ReadInProj(prd);
    pflow::CalcNeatProj(prd,only_occ);
    pflow::ReadInPDOSData(prd,pdd);
    pflow::CalcPDOS(prd,pdd);
    pdd.PrintPDOS(cout,prd.sp);
    if(pdd.print_params) pdd.PrintParams(cout,prd.LLMnames);
  }
} // namespace pflow

// ***************************************************************************
// pflow::PEARSON_SYMBOL
// ***************************************************************************
namespace pflow {
  string PEARSON_SYMBOL(istream& input) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "pflow::PEARSON_SYMBOL: BEGIN" << endl;
    stringstream sss;
    xstructure a(input,IOAFLOW_AUTO);
    //  cerr << a << endl;
    //  cerr << "here" << endl;
    if(LDEBUG) cerr << "pflow::PEARSON_SYMBOL: X1" << endl;
    //DX 8/24/17 [OBSOLETE] a.GetLatticeType();
    xstructure str_sp,str_sc;
    bool full_sym=false; //DX 8/29/17 - Speed increase
    LATTICE::Standard_Lattice_Structure(a,str_sp,str_sc,full_sym); //DX 8/29/17 - Speed increase
    if(LDEBUG) cerr << "pflow::PEARSON_SYMBOL: X2" << endl;
    if(LDEBUG) cerr << " Real space lattice primitive           = " << str_sp.bravais_lattice_type << endl; //DX 8/24/17 - a to str_sp
    if(LDEBUG) cerr << " Real space lattice variation           = " << str_sp.bravais_lattice_variation_type << endl;//wahyu mod //DX 8/24/17 - a to str_sp
    //  if(LDEBUG) cerr << " Real space conventional lattice        = " << a.bravais_conventional_lattice_type << endl; //DX 8/24/17 - a to str_sp
    if(LDEBUG) cerr << " Real space Pearson symbol              = " << str_sp.pearson_symbol << endl; //DX 8/24/17 - a to str_sp
    sss << str_sp.pearson_symbol << endl; //DX 8/24/17 - a to str_sp
    if(LDEBUG) cerr << "pflow::PEARSON_SYMBOL: END" << endl;
    return sss.str();
  }
} // namespace pflow

// ***************************************************************************
// pflow::PLANEDENS
// ***************************************************************************
namespace pflow {
  void PLANEDENS(vector<string> argv) {
    // Read in charge
    xstructure str;
    ifstream chgfile(argv.at(3).c_str());
    vector<int> ngrid(3);
    vector<double> chg_tot;
    vector<double> chg_diff;
    pflow::ReadChg(str,ngrid,chg_tot,chg_diff,chgfile);
    // Read in planar density parameters.
    pflow::pd_params pdp;
    ifstream pdfile(argv.at(2).c_str());
    pflow::ReadPlaneDensParams(str,pdp,pdfile);
    // Get planar charge density
    vector<double> dens2d_tot;
    vector<double> dens2d_diff;
    pflow::GetPlaneDens(pdp,dens2d_tot,dens2d_diff,str,ngrid,chg_tot,chg_diff);
    pflow::PrintPlaneDens(pdp,dens2d_tot,dens2d_diff,str);
  }
} // namespace pflow

// ***************************************************************************
// pflow::PLATON
// ***************************************************************************
namespace pflow {
  string PLATON(string options,istream& input) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "pflow::PLATON: BEGIN" << endl;
    vector<string> tokens;
    aurostd::string2tokens(options,tokens,",");
    if(tokens.size()==1) {
      if(tokens.at(0)=="usage" || tokens.at(0)=="USAGE") {
	init::ErrorOption(cout,options,"pflow::PLATON",
			  aurostd::liststring2string("aflow --platonSG[_label,_number][=EQUAL| EXACT][,ang,d1,d2,d3] < POSCAR  default:"+
						     string("EQUAL=")+aurostd::utype2string<int>(DEFAULT_PLATON_P_EQUAL)+","+
						     string("EXACT=")+aurostd::utype2string<int>(DEFAULT_PLATON_P_EXACT)+","+
						     aurostd::utype2string(DEFAULT_PLATON_P_ANG,5)+","+
						     aurostd::utype2string(DEFAULT_PLATON_P_D1,5)+","+
						     aurostd::utype2string(DEFAULT_PLATON_P_D2,5)+","+
						     aurostd::utype2string(DEFAULT_PLATON_P_D3,5)+""));
	exit(0);
      } 
    }
    // move on
    
    if(LDEBUG) cerr << "pflow::PLATON: tokens.size()=" << tokens.size() << endl;
    
    xstructure a(input,IOAFLOW_AUTO);
    bool Platon_EQUAL=DEFAULT_PLATON_P_EQUAL;
    bool Platon_EXACT=DEFAULT_PLATON_P_EXACT;
    double Platon_ang=DEFAULT_PLATON_P_ANG;
    double Platon_d1=DEFAULT_PLATON_P_D1;
    double Platon_d2=DEFAULT_PLATON_P_D2;
    double Platon_d3=DEFAULT_PLATON_P_D3;
    //     if(argv.size()==2) Platon_ang=1.0e-2;
    // if(argv.size()==3) Platon_ang=atof(tokens.at(0));
    // Read in input file.
    if((tokens.size()>=1 && tokens.at(0)=="EQUAL") || (tokens.size()>=2  && tokens.at(1)=="EQUAL")) Platon_EQUAL=TRUE;
    if((tokens.size()>=1 && tokens.at(0)=="EXACT") || (tokens.size()>=2 && tokens.at(1)=="EXACT")) Platon_EXACT=TRUE;
    if(tokens.size()>=3) {
      Platon_ang=aurostd::string2utype<double>(tokens.at(tokens.size()-4));
      Platon_d1=aurostd::string2utype<double>(tokens.at(tokens.size()-3));
      Platon_d2=aurostd::string2utype<double>(tokens.at(tokens.size()-2));
      Platon_d3=aurostd::string2utype<double>(tokens.at(tokens.size()-1));
    }
    a.FakeNames();
    return a.platon2print(Platon_EQUAL,Platon_EXACT,Platon_ang,Platon_d1,Platon_d2,Platon_d3);
    // aflow --platon [EQUAL] [EXACT] [ang d1 d2 d3]
    // CALC ADDSYM (EQUAL) (EXACT) (ang d1 d2 d3)
    // EQUAL - Search with all atom type treated as equivalent.
    // EXACT - All atoms should fit for given criteria.
    // ang - Angle criterium in search for metrical symmetry of the
    //       lattice (default 1.0 degree).
    // d1 - Distance criterium for coinciding atoms for non-inversion
    //      (pseudo)symmetry elements (default 0.25 Angstrom).
    // d2 - Distance criterium for coinciding atoms for (pseudo)
    //      inversion symmetry (default 0.45, 0.25 Angstrom).
    // d3 - Distance criterium for coinciding atoms for (pseudo)
    //      translation symmetry (default 0.45, 0.25 Angstrom).
    // SEE: http://www.cryst.chem.uu.nl/platon/pl000401.html
  }
} // namespace pflow

// ***************************************************************************
// pflow::POCC
// ***************************************************************************
namespace pflow {
  void POCC(vector<string> argv) {
    cerr << "# WARNING: THIS REQUIRES AN ALTERED VERSION OF VASP - SEE aflow --help" << endl;
    int only_occ=1; // Use occupied states only in projections.
    pflow::projdata proj_dat;
    proj_dat.PROOUTinfile=argv.at(2);
    pflow::ReadInProj(proj_dat);
    pflow::CalcNeatProj(proj_dat,only_occ);
    PrintNeatProj(proj_dat,cout);
  }
} // namespace pflow

// ***************************************************************************
// pflow::POSCAR
// ***************************************************************************
namespace pflow {
  xstructure POSCAR(istream& input) {
    xstructure str(input,IOAFLOW_AUTO);
    str.iomode=IOVASP_POSCAR;
    return str;
  }
} // namespace pflow

// ***************************************************************************
// pflow::POSCAR2WYCKOFF
// ***************************************************************************
namespace pflow {
  void POSCAR2WYCKOFF(istream & input) {
    // Call findsym to output the wyckoff position from POSCAR
    // Note that findsym must be installed properly
    string findsym_in=aurostd::TmpFileCreate("findsym.in");
  
    stringstream oss;
    // [OBSOLETE] xstructure str(input,IOVASP_POSCAR);
    xstructure str(input,IOAFLOW_AUTO);
    oss << str.findsym2print();
    aurostd::stringstream2file(oss,findsym_in);  
    vector<string> stroutput;
    FINDSYM::Write("data_space.txt","./");
    FINDSYM::Write("data_wyckoff.txt","./");
    aurostd::string2vectorstring(aurostd::execute2string(XHOST.command("findsym")+" < " + findsym_in),stroutput);
    FROZSL::Delete("data_space.txt","./");
    FROZSL::Delete("data_wyckoff.txt","./");

    bool flag=false;
    for(uint i=0;i<stroutput.size();i++) {
      if(aurostd::substring2bool(stroutput.at(i),"----")) flag=!flag;
      if(flag) if(!aurostd::substring2bool(stroutput.at(i),"----")) cout <<stroutput.at(i) << endl;
    }
    aurostd::execute("rm -f findsym.log " + findsym_in);
  }
} // namespace pflow

// ***************************************************************************
// pflow::PGROUP
// ***************************************************************************
namespace pflow {
  void PGROUP(_aflags &aflags,istream& input) {
    cout << aflow::Banner("BANNER_TINY") << endl;
    aflags.QUIET=TRUE;
    xstructure a(input,IOAFLOW_AUTO);
    bool WRITE=TRUE;
    ofstream File("/dev/null");
    _kflags kflags;                                   //DX 8/15/17 - Add in consistency checks
    kflags.KBIN_SYMMETRY_CALCULATE_PGROUP=TRUE;       //DX 8/15/17 - Add in consistency checks
    kflags.KBIN_SYMMETRY_CALCULATE_PGROUPK=FALSE;     //DX 8/15/17 - Add in consistency checks
    kflags.KBIN_SYMMETRY_CALCULATE_FGROUP=FALSE;      //DX 8/15/17 - Add in consistency checks
    kflags.KBIN_SYMMETRY_CALCULATE_PGROUP_XTAL=FALSE; //DX 8/15/17 - Add in consistency checks
    kflags.KBIN_SYMMETRY_CALCULATE_PGROUPK_XTAL=FALSE;//DX 8/15/17 - Add in consistency checks //DX 12/5/17 - Added pgroupk_xtal
    kflags.KBIN_SYMMETRY_CALCULATE_SGROUP=FALSE;      //DX 8/15/17 - Add in consistency checks
    kflags.KBIN_SYMMETRY_CALCULATE_IATOMS=FALSE;      //DX 8/15/17 - Add in consistency checks
    kflags.KBIN_SYMMETRY_CALCULATE_AGROUP=FALSE;      //DX 8/15/17 - Add in consistency checks
    pflow::PerformFullSymmetry(a,File,aflags,kflags,WRITE,cout); //DX 8/15/17 - Add in consistency checks
    // DX 8/15/17 - Add in consistency checks - SYM::CalculatePointGroup(File,a,aflags,WRITE,TRUE,cout);
  }
} // namespace pflow

// ***************************************************************************
// pflow::PGROUPK
// ***************************************************************************
namespace pflow {
  void PGROUPK(_aflags &aflags,istream& input) {
    cout << aflow::Banner("BANNER_TINY") << endl;
    aflags.QUIET=TRUE;
    xstructure a(input,IOAFLOW_AUTO);
    bool WRITE=TRUE;
    ofstream File("/dev/null");
    //DX 8/15/17 SYM::CalculatePointGroupKlattice(File,a,aflags,WRITE,TRUE,cout);
    _kflags kflags;                                   //DX 8/15/17 - Add in consistency checks
    kflags.KBIN_SYMMETRY_CALCULATE_PGROUP=TRUE;       //DX 8/15/17 - Add in consistency checks
    kflags.KBIN_SYMMETRY_CALCULATE_PGROUPK=TRUE;      //DX 8/15/17 - Add in consistency checks
    kflags.KBIN_SYMMETRY_CALCULATE_FGROUP=FALSE;      //DX 8/15/17 - Add in consistency checks
    kflags.KBIN_SYMMETRY_CALCULATE_PGROUP_XTAL=FALSE; //DX 8/15/17 - Add in consistency checks
    kflags.KBIN_SYMMETRY_CALCULATE_PGROUPK_XTAL=FALSE;//DX 8/15/17 - Add in consistency checks //DX 12/5/17 - Added pgroupk_xtal
    kflags.KBIN_SYMMETRY_CALCULATE_SGROUP=FALSE;      //DX 8/15/17 - Add in consistency checks
    kflags.KBIN_SYMMETRY_CALCULATE_IATOMS=FALSE;      //DX 8/15/17 - Add in consistency checks
    kflags.KBIN_SYMMETRY_CALCULATE_AGROUP=FALSE;      //DX 8/15/17 - Add in consistency checks
    pflow::PerformFullSymmetry(a,File,aflags,kflags,WRITE,cout); //DX 8/15/17 - Add in consistency checks
  }
} // namespace pflow

// ***************************************************************************
// pflow::PGROUPXTAL
// ***************************************************************************
namespace pflow {
  void PGROUPXTAL(_aflags &aflags,istream& input) {
    //  cout << aflow::Banner("BANNER_TINY") << endl;
    aflags.QUIET=TRUE;
    xstructure a(input,IOAFLOW_AUTO);
    bool WRITE=TRUE;
    ofstream File("/dev/null");
    // DX 8/15/17 - Add in consistency checks bool verbose=TRUE;
    // DX 8/15/17 - Add in consistency checks SYM::CalculatePointGroup(File,a,aflags,WRITE,verbose,cout);
    //DX SYM::CalculateFactorGroup(File,a,aflags,WRITE,verbose,cout);
    //  SYM::CalculateSpaceGroup(File,a,aflags,FALSE,verbose,cout);
    // SYM::CalculateInequivalentAtoms(File,a,aflags,WRITE,verbose,cout);
    // SYM::CalculateSitePointGroup(File,a,aflags,WRITE,TRUE,cout);
    // DX 8/15/17 - Add in consistency checks SYM::CalculatePointGroupCrystal(File,a,aflags,WRITE,verbose,cout);
    _kflags kflags;                                   //DX 8/15/17 - Add in consistency checks
    kflags.KBIN_SYMMETRY_CALCULATE_PGROUP=TRUE;       //DX 8/15/17 - Add in consistency checks
    kflags.KBIN_SYMMETRY_CALCULATE_PGROUPK=FALSE;     //DX 8/15/17 - Add in consistency checks
    kflags.KBIN_SYMMETRY_CALCULATE_FGROUP=TRUE;       //DX 8/15/17 - Add in consistency checks
    kflags.KBIN_SYMMETRY_CALCULATE_PGROUP_XTAL=TRUE;  //DX 8/15/17 - Add in consistency checks
    kflags.KBIN_SYMMETRY_CALCULATE_PGROUPK_XTAL=FALSE;//DX 8/15/17 - Add in consistency checks //DX 12/5/17 - Added pgroupk_xtal
    kflags.KBIN_SYMMETRY_CALCULATE_SGROUP=FALSE;      //DX 8/15/17 - Add in consistency checks
    kflags.KBIN_SYMMETRY_CALCULATE_IATOMS=FALSE;      //DX 8/15/17 - Add in consistency checks
    kflags.KBIN_SYMMETRY_CALCULATE_AGROUP=FALSE;      //DX 8/15/17 - Add in consistency checks
    pflow::PerformFullSymmetry(a,File,aflags,kflags,WRITE,cout); //DX 8/15/17 - Add in consistency checks
  }
} // namespace pflow

// ***************************************************************************
// pflow::PRIM
// ***************************************************************************
namespace pflow {
  xstructure PRIM(istream& input,uint mode) {
    xstructure a(input,IOAFLOW_AUTO);
    if(mode==0) return GetPrimitive(a,0.005);
    if(mode==1) return GetPrimitive1(a);
    if(mode==2) return GetPrimitive2(a);
    if(mode==3) return GetPrimitive3(a);
    return a;
  }
} // namespace pflow

// ***************************************************************************
bool RequestedAlphabeticLabeling(string &label) {
  if(aurostd::substring2bool(label,".alphabetic")) {aurostd::StringSubst(label,".alphabetic","");return TRUE;}
  if(aurostd::substring2bool(label,".alpha")) {aurostd::StringSubst(label,".alpha","");return TRUE;}
  return FALSE;
}

bool AlphabetizePrototypeLabelSpecies(deque<string> &species,deque<string> &species_pp,deque<double> &vvolume,deque<double> &vmass,string &label) {
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  // DEBUG=TRUE;
  uint nspeciesHTQC=species.size();
  if(LDEBUG) cerr << "AlphabetizePrototypeLabelSpecies: species.size()=" << species.size() << endl;
  if(LDEBUG) cerr << "AlphabetizePrototypeLabelSpecies: species_pp.size()=" << species_pp.size() << endl;
  if(LDEBUG) cerr << "AlphabetizePrototypeLabelSpecies: vvolume.size()=" << vvolume.size() << endl;
  if(LDEBUG) cerr << "AlphabetizePrototypeLabelSpecies: vmass.size()=" << vmass.size() << endl;
  if(LDEBUG) cerr << "AlphabetizePrototypeLabelSpecies: nspeciesHTQC=" << nspeciesHTQC << endl;
  // if(LDEBUG) cerr << "AlphabetizePrototypeLabelSpecies: alphabetic=" << alphabetic << endl;
  if(LDEBUG) cerr << "AlphabetizePrototypeLabelSpecies: label=" << label << endl;
  aurostd::StringSubst(label,".alphabetic","");
  aurostd::StringSubst(label,".alpha","");
  deque<string> rnd_species,rnd_species_pp,rnd_label1,rnd_label2;
  deque<double> rnd_vvolume;
  deque<double> rnd_vmass;
  if(nspeciesHTQC!=species.size()) { cerr << "AlphabetizePrototypeLabelSpecies: ERROR nspeciesHTQC!=species.size" << endl;exit(0);}
  if(nspeciesHTQC!=species_pp.size()) { cerr << "AlphabetizePrototypeLabelSpecies: ERROR nspeciesHTQC!=species_pp.size" << endl;exit(0);}
  if(nspeciesHTQC!=vvolume.size()) { cerr << "AlphabetizePrototypeLabelSpecies: ERROR nspeciesHTQC!=vvolume.size" << endl;exit(0);}
  if(nspeciesHTQC!=vmass.size()) { cerr << "AlphabetizePrototypeLabelSpecies: ERROR nspeciesHTQC!=vmass.size" << endl;exit(0);}
  for(uint i=0;i<nspeciesHTQC;i++) {
    rnd_species.push_back(species.at(i));   // LOAD FROM SPECIES
    rnd_species_pp.push_back(species_pp.at(i));   // LOAD FROM SPECIES
    rnd_vvolume.push_back(vvolume.at(i));  // LOAD FROM SPECIES
    rnd_vmass.push_back(vmass.at(i));  // LOAD FROM SPECIES
    string A="A";A[0]+=i;rnd_label1.push_back(A);rnd_label2.push_back(A); // megatrick thanks C++ to have C inside (SC2011)
  }
  if(0) {
    aurostd::sort(rnd_species,rnd_species_pp,rnd_vvolume,rnd_vmass,rnd_label1);
    label=label+".";for(uint i=0;i<rnd_label1.size();i++) label=label+rnd_label1.at(i);  // FIXING LABEL
  }
  if(1) { // this works
    aurostd::sort(rnd_species,rnd_species_pp,rnd_vvolume,rnd_vmass,rnd_label1);   // how to fix rnd_atomxX by shufflying rnd_label1
    aurostd::sort(rnd_label1,rnd_label2);             // how to go back from rnd_labe1 to rnd_label2
    label=label+".";for(uint i=0;i<rnd_label2.size();i++) label=label+rnd_label2.at(i);  // FIXING LABEL
  }  
  if(LDEBUG) cerr << "AlphabetizePrototypeLabelSpecies: label=" << label << endl;
  for(uint i=0;i<nspeciesHTQC;i++)
    if(LDEBUG) cerr << "AlphabetizePrototypeLabelSpecies: BEFORE species.at(" << i << ")=" 
		    << species.at(i) << " species_pp.at(" << i << ")=" << species_pp.at(i) 
		    <<  " vvolume.at(" << i << ")=" << vvolume.at(i) 
		    <<  " vmass.at(" << i << ")=" << vmass.at(i) << endl;
  for(uint i=0;i<nspeciesHTQC;i++) species.at(i)=rnd_species.at(i);
  for(uint i=0;i<nspeciesHTQC;i++) species_pp.at(i)=rnd_species_pp.at(i);
  for(uint i=0;i<nspeciesHTQC;i++) vvolume.at(i)=rnd_vvolume.at(i);
  for(uint i=0;i<nspeciesHTQC;i++) vmass.at(i)=rnd_vmass.at(i);
  for(uint i=0;i<nspeciesHTQC;i++)
    if(LDEBUG) cerr << "AlphabetizePrototypeLabelSpecies: AFTER species.at(" << i << ")=" 
		    << species.at(i) << " species_pp.at(" << i << ")=" 
		    << species_pp.at(i) <<  " vvolume.at(" << i << ")="
		    << vvolume.at(i) <<  " vmass.at(" << i << ")=" << vmass.at(i) << endl;
  if(LDEBUG) {cerr << "AlphabetizePrototypeLabelSpecies: EXIT" << endl;  exit(0);}
  return TRUE;
}

bool AlphabetizePrototypeLabelSpecies(deque<string> &species,deque<string> &species_pp,string &label) {
  deque<double> vvolume;for(uint i=0;i<species.size();i++) vvolume.push_back(double(0));
  deque<double> vmass;for(uint i=0;i<species.size();i++) vmass.push_back(double(0));
  return AlphabetizePrototypeLabelSpecies(species,species_pp,vvolume,vmass,label);
}

bool AlphabetizePrototypeLabelSpecies(deque<string> &species,deque<double> &vvolume,string &label) {
  deque<string> species_pp;for(uint i=0;i<species.size();i++) species_pp.push_back(species.at(i));
  deque<double> vmass;for(uint i=0;i<species.size();i++) vmass.push_back(double(0));
  return AlphabetizePrototypeLabelSpecies(species,species_pp,vvolume,vmass,label);
}

bool AlphabetizePrototypeLabelSpecies(deque<string> &species,string &label) {
  deque<string> species_pp;for(uint i=0;i<species.size();i++) species_pp.push_back(species.at(i));
  deque<double> vvolume;for(uint i=0;i<species.size();i++) vvolume.push_back(double(0));
  deque<double> vmass;for(uint i=0;i<species.size();i++) vmass.push_back(double(0));
  return AlphabetizePrototypeLabelSpecies(species,species_pp,vvolume,vmass,label);
}

string AlphabetizePrototypeLabelSpeciesArgv(vector<string> &argv) {
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  //  LDEBUG=TRUE;
  string label=argv.at(2);
  uint nspeciesHTQC=aflowlib::PrototypeLibrariesSpeciesNumber(label);
  if(LDEBUG) cerr << "AlphabetizePrototypeLabelSpeciesArgv: nspeciesHTQC=" << nspeciesHTQC << endl;
  if(LDEBUG) cerr << "AlphabetizePrototypeLabelSpeciesArgv: label=" << label << endl;
  deque<string> species;
  for(uint i=0;i<nspeciesHTQC;i++) species.push_back(argv.at(3+i));  // LOAD FROM ARGV
  AlphabetizePrototypeLabelSpecies(species,label);
  for(uint i=0;i<nspeciesHTQC;i++) argv.at(3+i)=species.at(i);
  argv.at(2)=label;
  return label;
}

string AlphabetizePrototypeLabelSpeciesTokens(vector<string> &tokens) {
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  //  LDEBUG=TRUE;
  string label=tokens.at(0);
  uint nspeciesHTQC=aflowlib::PrototypeLibrariesSpeciesNumber(label);
  if(LDEBUG) cerr << "AlphabetizePrototypeLabelSpeciesTokens: nspeciesHTQC=" << nspeciesHTQC << endl;
  if(LDEBUG) cerr << "AlphabetizePrototypeLabelSpeciesTokens: label=" << label << endl;
  deque<string> species;
  for(uint i=0;i<nspeciesHTQC;i++) species.push_back(tokens.at(1+i));  // LOAD FROM TOKENS
  AlphabetizePrototypeLabelSpecies(species,label);
  for(uint i=0;i<nspeciesHTQC;i++) tokens.at(1+i)=species.at(i);
  tokens.at(0)=label;
  return label;
}

// ***************************************************************************
// pflow::PROTO_LIBRARIES
// ***************************************************************************
namespace pflow {
  xstructure PROTO_LIBRARIES(aurostd::xoption vpflow) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "pflow::PROTO_LIBRARIES: BEGIN" << endl;  

    vector<string> params;
    if(LDEBUG) cerr << "pflow::PROTO: vpflow.getattachedscheme(\"PROTO\")=" << vpflow.getattachedscheme("PROTO") << endl;
    if(LDEBUG) cerr << "pflow::PROTO: vpflow.getattachedscheme(\"PROTO_ICSD_AFLOW\")=" << vpflow.getattachedscheme("PROTO_ICSD_AFLOW") << endl;
    if(LDEBUG) cerr << "pflow::PROTO: vpflow.flag(\"PROTO\")=" << vpflow.flag("PROTO") << endl;
    if(LDEBUG) cerr << "pflow::PROTO: vpflow.flag(\"PARAMS\")=" << vpflow.flag("PARAMS") << endl;
    if(LDEBUG) cerr << "pflow::PROTO: vpflow.flag(\"PROTO_ICSD_AFLOW\")=" << vpflow.flag("PROTO_ICSD_AFLOW") << endl;
    aurostd::string2tokens(vpflow.getattachedscheme("PROTO"),params,":");
    
    if(params.size()==0) {
      init::ErrorOption(cout,vpflow.getattachedscheme("PROTO"),
			"pflow::PROTO",
			aurostd::liststring2string("aflow [options] --proto=label*[:speciesA*[:speciesB*]..[:volumeA*[:volumeB*].. | :volume]] [--params=..... [--hex]]",
                                                   "                --proto[_icsd]=ICSD_number.{ABC}[:speciesA*[:speciesB*]..[:volumeA*[:volumeB*].. | :volume]]",
                                                   "                --proto_icsd=label_ICSD_number",
                                                   "options = [--server=xxxxxxx]  [--vasp | --qe | --abinit | --aims] [--params=... | --hex]]",
                                                   "To get the list of prototypes --protos or --protos_icsd ."));
      exit(0);
    }
    
    string label=params.at(0);
    string parameters=vpflow.getattachedscheme("PARAMS");
    
    //  cerr << params.at(0) << " " << params.size() << endl;  for(uint i=0;i<params.size();i++) cerr << params.at(i) << " ";  cerr << endl;

    bool modeLIBRARY=TRUE;//FALSE;
    bool modeFILE=FALSE;
    bool modeCIN=FALSE;
    
    xstructure str;
    deque<string> atomX;
    deque<double> volumeX;
    uint nspecies=0;
    bool alphabetic=TRUE; // DEFAULT

    
    // ***************************************************************************
    // MODE LIBRARY FROM HTQC OR ICSD      
    if(modeLIBRARY) {
      nspecies=aflowlib::PrototypeLibrariesSpeciesNumber(label);
      alphabetic=RequestedAlphabeticLabeling(label);
      
      int mode=LIBRARY_MODE_HTQC;
      if(XHOST.vflag_control.flag("AFLOWLIB_SERVER")) {
	if(aurostd::substring2bool(label,"_ICSD_")) mode=LIBRARY_MODE_ICSD_AFLOWLIB;
	if(aurostd::substring2bool(label,"ICSD_") && !aurostd::substring2bool(label,"_ICSD_")) mode=LIBRARY_MODE_HTQC_ICSD_AFLOWLIB;
      } else {
	if(aurostd::substring2bool(label,"_ICSD_")) mode=LIBRARY_MODE_ICSD;
	if(aurostd::substring2bool(label,"ICSD_") && !aurostd::substring2bool(label,"_ICSD_")) mode=LIBRARY_MODE_HTQC_ICSD;
      } 
	
      // DEBUG=TRUE;
      if(LDEBUG) cerr << "pflow::PROTO_LIBRARIES: params.size()=" << params.size() << endl;
      if(LDEBUG) cerr << "pflow::PROTO_LIBRARIES: nspecies=" << nspecies << endl;
      if(LDEBUG) cerr << "pflow::PROTO_LIBRARIES: alphabetic=" << alphabetic << endl;
      if(LDEBUG) cerr << "pflow::PROTO_LIBRARIES: label=" << label << endl;
      if(LDEBUG) cerr << "pflow::PROTO_LIBRARIES: parameters=" << parameters << endl;
      if(LDEBUG) if(mode==LIBRARY_MODE_HTQC) cerr << "pflow::PROTO_LIBRARIES: mode=LIBRARY_MODE_HTQC" << endl;
      if(LDEBUG) if(mode==LIBRARY_MODE_ICSD) cerr << "pflow::PROTO_LIBRARIES: mode=LIBRARY_MODE_ICSD" << endl;
      if(LDEBUG) if(mode==LIBRARY_MODE_HTQC_ICSD) cerr << "pflow::PROTO_LIBRARIES: mode=LIBRARY_MODE_HTQC_ICSD" << endl;
      if(LDEBUG) if(mode==LIBRARY_MODE_ICSD_AFLOWLIB) cerr << "pflow::PROTO_LIBRARIES: mode=LIBRARY_MODE_ICSD_AFLOWLIB" << endl;
      if(LDEBUG) if(mode==LIBRARY_MODE_HTQC_ICSD_AFLOWLIB) cerr << "pflow::PROTO_LIBRARIES: mode=LIBRARY_MODE_HTQC_ICSD_AFLOWLIB" << endl;
      
      bool found=FALSE;
      
      if(!found && params.size()==1) {  // label: 	// no matter the species, print A
	if(LDEBUG) cerr << "pflow::PROTO_LIBRARIES: label" << endl;
 	str=aflowlib::PrototypeLibraries(cerr,label,parameters,mode);found=TRUE; // good for binary, ternary etc etc...
      }  
      
      if(!found && params.size()==1+nspecies) { // label:species:species...
	if(LDEBUG) cerr << "pflow::PROTO_LIBRARIES: label:species:species..." << endl;
 	atomX.clear();for(uint i=0;i<nspecies;i++) atomX.push_back(params.at(i+1));
	if(LDEBUG) for(uint i=0;i<atomX.size();i++) cerr << "pflow::PROTO_LIBRARIES: atomX.at(" << i << ")=" << atomX.at(i) << endl;
	str=aflowlib::PrototypeLibraries(cerr,label,parameters,atomX,mode);found=TRUE;
      }
      
      if(!found && params.size()==1+nspecies+nspecies) { // label:species:species..:volume:volume...
	if(LDEBUG) cerr << "pflow::PROTO_LIBRARIES: label:species:species...:volume:volume..." << endl;
	atomX.clear();for(uint i=0;i<nspecies;i++) atomX.push_back(params.at(i+1));
	if(LDEBUG) for(uint i=0;i<atomX.size();i++) cerr << "pflow::PROTO_LIBRARIES: atomX.at(" << i << ")=" << atomX.at(i) << endl;
	volumeX.clear();for(uint i=0;i<nspecies;i++) volumeX.push_back(aurostd::string2utype<double>(params.at(i+1+nspecies)));
	if(LDEBUG) for(uint i=0;i<volumeX.size();i++) cerr << "pflow::PROTO_LIBRARIES: volumeX.at(" << i << ")=" << volumeX.at(i) << endl;	
	str=aflowlib::PrototypeLibraries(cerr,label,parameters,atomX,volumeX,-1.0,mode);found=TRUE;
      }
      
      if(!found && params.size()==2 && mode!=LIBRARY_MODE_HTQC_ICSD) {
	atomX.clear();for(uint i=1;i<params.size();i++) atomX.push_back(params.at(i));
	if(LDEBUG) for(uint i=0;i<atomX.size();i++) cerr << "pflow::PROTO_LIBRARIES: atomX.at(" << i << ")=" << atomX.at(i) << endl;
	if(nspecies==1 && mode!=LIBRARY_MODE_HTQC_ICSD) {
	  str=aflowlib::PrototypeLibraries(cerr,label,parameters,atomX,mode);found=TRUE;
	}
	if(mode==LIBRARY_MODE_HTQC_ICSD) {
	  atomX.clear();atomX.push_back(params.at(1));
	  if(LDEBUG) for(uint i=0;i<atomX.size();i++) cerr << "pflow::PROTO_LIBRARIES: atomX.at(" << i << ")=" << atomX.at(i) << endl;
	  str=aflowlib::PrototypeLibraries(cerr,label,parameters,atomX,mode);found=TRUE;
	}
      }
  
      if(!found && params.size()==3) {
	atomX.clear();for(uint i=1;i<params.size();i++) atomX.push_back(params.at(i));
	if(LDEBUG) for(uint i=0;i<atomX.size();i++) cerr << "pflow::PROTO_LIBRARIES: atomX.at(" << i << ")=" << atomX.at(i) << endl;
	if(nspecies==2 && mode!=LIBRARY_MODE_HTQC_ICSD) {
	  str=aflowlib::PrototypeLibraries(cerr,label,parameters,atomX,mode);found=TRUE;
	}
	if(mode==LIBRARY_MODE_HTQC_ICSD) {
	  str=aflowlib::PrototypeLibraries(cerr,label,parameters,atomX,mode);found=TRUE;
	}
      }
  
      if(!found && params.size()==4) {
	atomX.clear();for(uint i=1;i<params.size();i++) atomX.push_back(params.at(i));
	if(nspecies==2 && mode!=LIBRARY_MODE_HTQC_ICSD) {
	  atomX.clear();atomX.push_back(params.at(1));atomX.push_back(params.at(2));
	  if(LDEBUG) for(uint i=0;i<atomX.size();i++) cerr << "pflow::PROTO_LIBRARIES: atomX.at(" << i << ")=" << atomX.at(i) << endl;
	  volumeX.clear();volumeX.push_back(0.0);volumeX.push_back(0.0);
	  if(LDEBUG) for(uint i=0;i<volumeX.size();i++) cerr << "pflow::PROTO_LIBRARIES: volumeX.at(" << i << ")=" << volumeX.at(i) << endl;	
	  str=aflowlib::PrototypeLibraries(cerr,label,parameters,atomX,volumeX,aurostd::string2utype<double>(params.at(3)),mode);found=TRUE;
	}
	if(nspecies==3 && mode!=LIBRARY_MODE_HTQC_ICSD) { // ternaries
	  atomX.clear();for(uint i=0;i<nspecies;i++) atomX.push_back(params.at(1+i));
	  if(LDEBUG) for(uint i=0;i<atomX.size();i++) cerr << "pflow::PROTO_LIBRARIES: atomX.at(" << i << ")=" << atomX.at(i) << endl;
	  str=aflowlib::PrototypeLibraries(cerr,label,parameters,atomX,mode);found=TRUE;
	}
	if(mode==LIBRARY_MODE_HTQC_ICSD) {
	  str=aflowlib::PrototypeLibraries(cerr,label,parameters,atomX,mode);found=TRUE;
	}
      }
  
      if(!found && params.size()==5) {
	if(nspecies==2) {
	  atomX.clear();atomX.push_back(params.at(1));atomX.push_back(params.at(2));
	  if(LDEBUG) for(uint i=0;i<atomX.size();i++) cerr << "pflow::PROTO_LIBRARIES: atomX.at(" << i << ")=" << atomX.at(i) << endl;
	  volumeX.clear();volumeX.push_back(aurostd::string2utype<double>(params.at(3)));volumeX.push_back(aurostd::string2utype<double>(params.at(4)));
	  if(LDEBUG) for(uint i=0;i<volumeX.size();i++) cerr << "pflow::PROTO_LIBRARIES: volumeX.at(" << i << ")=" << volumeX.at(i) << endl;	
	  str=aflowlib::PrototypeLibraries(cerr,label,parameters,atomX,volumeX,-1.0,mode);found=TRUE;
	}
	if(nspecies==4) {  // quaternaries
	  if(alphabetic==TRUE) label=AlphabetizePrototypeLabelSpeciesTokens(params);
	  atomX.clear();for(uint i=0;i<nspecies;i++) atomX.push_back(params.at(1+i));
	  if(LDEBUG) for(uint i=0;i<atomX.size();i++) cerr << "pflow::PROTO_LIBRARIES: atomX.at(" << i << ")=" << atomX.at(i) << endl;
	  str=aflowlib::PrototypeLibraries(cerr,label,parameters,atomX,mode);found=TRUE;
	}
      }
    }
    
    // ***************************************************************************
    // MODE FILE      
    if(modeFILE) {
      cerr << "pflow::PROTO_LIBRARIES: modeFILE=" << modeFILE << endl;
    }
    
    // ***************************************************************************
    // MODE CIN      
    if(modeCIN) {
      cerr << "pflow::PROTO_LIBRARIES: modeCIN=" << modeCIN << endl;
    }
    
    if(LDEBUG) for(uint i=0;i<str.species.size();i++) if(str.species.at(i)!="") cerr << "DEBUG specie(" << i << ")=" << str.species.at(i) << endl;
    //  if(str.AlphabeticSpecie(0,1)==FALSE) str.SwapSpecie(0,1);
    
    // now checking QUANTUM ESPRESSO
    if(vpflow.flag("PROTO::QE")) str.xstructure2qe();

    // now checking ABINIT
    if(vpflow.flag("PROTO::ABINIT")) str.xstructure2abinit();

    // now checking AIMS
    if(vpflow.flag("PROTO::AIMS")) str.xstructure2aims();

    if(LDEBUG) cerr << "pflow::PROTO_LIBRARIES: END" << endl;  
    return str;
  }
} // namespace pflow

// ***************************************************************************
// pflow::PROTO_AFLOW
// ***************************************************************************
// struct _AVASP_PROTO {
//   vector<string> ucell;
//   vector<int> vkppra;
//   vector<double> vpressure;
//   aurostd::xoptions vparams;
// };


namespace pflow {
  bool PROTO_AFLOW(aurostd::xoption vpflow,bool flag_REVERSE) { // too many options
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: BEGIN" << endl; 
    
    // check prototypes
    _AVASP_PROTO PARAMS;
    PARAMS.ucell.clear();
 
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: vpflow.getattachedscheme(\"PROTO_AFLOW\")=" << vpflow.getattachedscheme("PROTO_AFLOW") << endl;
    aurostd::string2tokens(vpflow.getattachedscheme("PROTO_AFLOW"),PARAMS.ucell,":");
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.ucell.size()=" << PARAMS.ucell.size() << endl;

    // check usage
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: CHECK USAGE" << endl; 
    if(vpflow.flag("PROTO_AFLOW::USAGE") || PARAMS.ucell.size()==0) {
      init::ErrorOption(cout,vpflow.getattachedscheme("PROTO_AFLOW"),"pflow::PROTO_AFLOW",
			aurostd::liststring2string("aflow [options] --aflow_proto[|_icsd]=label*:speciesA*[:speciesB*][:volumeA*[:volumeB*]|:volume] [--params=... [--hex]]",
						   "       options:",
						   "                --usage",
						   "                --potential=pot_LDA | pot_GGA | potpaw_LDA | potpaw_GGA | potpaw_PBE",
						   "                --potential_complete",
						   "                --module=[APL | QHA | AAPL]",
						   "                --apl_superce=NxNxN",
						   "                --usage",
						   "                --missing",
						   "                --noautopp",
						   "                --kppra=NNNN (default: DEFAULT_KPPRA in .aflow.rc) --kppra_static=NNNN (default: DEFAULT_KPPRA_STATIC in .aflow.rc) --bands_grid=NNNN (default: DEFAULT_BANDS_GRID in .aflow.rc)",
                                                   "                --enmax_multiply=NNNN (default:VASP_PREC_ENMAX_XXXX in .aflow.rc)", 
						   "                --pressure=0,1,2 (kB) (default:0.0)",
						   "                --potim=XXX (default 0.05)  (VASP) ",
						   "                --relax_type=[ALL | IONS | CELL_SHAPE | CELL_VOLUME | IONS_CELL_VOLUME] ",
						   "                --relax_mode=[ENERGY | FORCES | ENERGY_FORCES | FORCES_ENERGY] (default: DEFAULT_VASP_FORCE_OPTION_RELAX_MODE_SCHEME in .aflow.rc) (VASP) ",
						   "                --precision=[(LOW | MEDIUM | NORMAL | HIGH | ACCURATE), PRESERVED] (default: DEFAULT_VASP_FORCE_OPTION_PREC_SCHEME in .aflow.rc) (VASP) ",
						   "                --algorithm=[(NORMAL | VERYFAST | FAST | ALL | DAMPED), PRESERVED] (default: DEFAULT_VASP_FORCE_OPTION_ALGO_SCHEME in .aflow.rc) (VASP) ",
						   "                --metagga=[TPSS | RTPSS | M06L | MBJL | SCAN | MS0 | MS1 | MS2 | NONE] (default: DEFAULT_VASP_FORCE_OPTION_METAGGA_SCHEME in .aflow.rc) (VASP) ",
						   "                --ivdw=[number_for_VASP_see_manual_for_IVDW | 0] (default: DEFAULT_VASP_FORCE_OPTION_IVDW_SCHEME in .aflow.rc) (VASP) ",
						   "                --type=[METAL | INSULATOR | SEMICONDUCTOR | DEFAULT] (default: DEFAULT_VASP_FORCE_OPTION_TYPE_SCHEME in .aflow.rc) (VASP) ",
						   "                --convert_unit_cell= [SPRIM, SCONV, NIGGLI, MINK, INCELL, COMPACT, WS, CART, FRAC, PRES] ",
						   "                --volume_plus_equal=XXX ",
						   "                --volume_multiply_equal=XXX ",
						   "                --no_volume_adjustment ",
						   "                --ediffg=XXX  (default: DEFAULT_VASP_PREC_EDIFFG in .aflow.rc) (VASP) ",
						   "                --ldau2",
						   "                --noldau2",
						   "                --bands",
						   "                --neglect_nomix",
						   "                --stdout",
						   "                --qe",
						   "                --abinit",
						   "                --aims",
						   "                --list",
						   "                --params=....  { check aflow --readme=anrl }",
						   "                --hex          { check aflow --readme=anrl }"));
      exit(0);
    }
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: vpflow.getattachedscheme(\"PROTO_AFLOW::USAGE\")=" << vpflow.flag("PROTO_AFLOW::USAGE") << endl;

    // label
    string string_LABEL=PARAMS.ucell.at(0);
  
    //    if(LDEBUG) {    cerr << "pflow::PROTO_AFLOW: " << string_LABEL << " " << PARAMS.ucell.size() << endl;     for(uint i=0;i<PARAMS.ucell.size();i++) cerr << PARAMS.ucell.at(i) << " "; cerr << endl;}

    // check potential
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: CHECK POTENTIAL" << endl;
  
    string string_POTENTIAL=_AVASP_PSEUDOPOTENTIAL_AUTO_;
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: vpflow.getattachedscheme(\"PROTO_AFLOW::POTENTIAL\")=" << vpflow.getattachedscheme("PROTO_AFLOW::POTENTIAL") << endl;
    string params_potential=vpflow.getattachedscheme("PROTO_AFLOW::POTENTIAL");
    if(params_potential==_AVASP_PSEUDOPOTENTIAL_AUTO_) string_POTENTIAL=_AVASP_PSEUDOPOTENTIAL_AUTO_;
    aurostd::StringSubst(params_potential,"pot","");aurostd::StringSubst(params_potential,"POT","");
    aurostd::StringSubst(params_potential,"_","");aurostd::StringSubst(params_potential,"paw","PAW");
    if(params_potential=="LDA") string_POTENTIAL=_AVASP_PSEUDOPOTENTIAL_POT_LDA_;
    if(params_potential=="GGA") string_POTENTIAL=_AVASP_PSEUDOPOTENTIAL_POT_GGA_;
    if(params_potential=="PBE") string_POTENTIAL=_AVASP_PSEUDOPOTENTIAL_POT_PBE_;
    if(params_potential=="PAWLDA") string_POTENTIAL=_AVASP_PSEUDOPOTENTIAL_POTPAW_LDA_;
    if(params_potential=="PAWGGA") string_POTENTIAL=_AVASP_PSEUDOPOTENTIAL_POTPAW_GGA_;
    if(params_potential=="PAWPBE") string_POTENTIAL=_AVASP_PSEUDOPOTENTIAL_POTPAW_PBE_;
    if(LDEBUG) {cerr << "pflow::PROTO_AFLOW: string_POTENTIAL=" << string_POTENTIAL << endl;}
    
    //   PARAMS.vparams.flag("AFLOWING_STRING_POTENTIAL",vpflow.flag("PROTO_AFLOW::POTENTIAL"));
    // if(PARAMS.vparams.flag("AFLOWING_STRING_POTENTIAL")) PARAMS.vparams.addattachedscheme("AFLOWING_STRING_POTENTIAL",vpflow.getattachedscheme("PROTO_AFLOW::POTENTIAL"),PARAMS.vparams.flag("AFLOWING_STRING_POTENTIAL"));

    // check potential_complete
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: CHECK POTENTIAL_COMPLETE" << endl; 
    if(vpflow.flag("PROTO_AFLOW::POTENTIAL_COMPLETE")) string_POTENTIAL+=_AVASP_PSEUDOPOTENTIAL_DELIMITER_+_AVASP_PSEUDOPOTENTIAL_POTENTIAL_COMPLETE_;
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: vpflow.getattachedscheme(\"PROTO_AFLOW::POTENTIAL_COMPLETE\")=" << vpflow.flag("PROTO_AFLOW::POTENTIAL_COMPLETE") << endl;
 
    PARAMS.vparams.addattachedscheme("AFLOWIN_STRING::POTENTIAL",string_POTENTIAL,TRUE);

    // reverse
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: CHECK REVERSE" << endl; 
    PARAMS.vparams.flag("AFLOWIN_FLAG::REVERSE",flag_REVERSE);
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.flag(\"AFLOWIN_FLAG::REVERSE\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::REVERSE") << endl;
  
    // check missing
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: CHECK MISSING" << endl; 
    PARAMS.vparams.flag("AFLOWIN_FLAG::MISSING",vpflow.flag("PROTO_AFLOW::MISSING"));
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.flag(\"AFLOWIN_FLAG::MISSING\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::MISSING") << endl;
  
    // check noautopp
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: CHECK NOAUTOPP" << endl; 
    PARAMS.vparams.flag("AFLOWIN_FLAG::NOAUTOPP",vpflow.flag("PROTO_AFLOW::NOAUTOPP"));
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.flag(\"AFLOWIN_FLAG::NOAUTOPP\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::NOAUTOPP") << endl;
    
    // check ldau
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: CHECK LDAU" << endl; 
    PARAMS.vparams.flag("AFLOWIN_FLAG::AUTOLDAU",vpflow.flag("PROTO_AFLOW::LDAU"));
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.flag(\"AFLOWIN_FLAG::AUTOLDAU\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::AUTOLDAU") << endl;
    
    // check noldau
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: CHECK NOLDAU" << endl; 
    PARAMS.vparams.flag("AFLOWIN_FLAG::AUTONOLDAU",vpflow.flag("PROTO_AFLOW::NOLDAU"));
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.flag(\"AFLOWIN_FLAG::AUTONOLDAU\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::AUTONOLDAU") << endl;
    
    // check bands
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: CHECK BANDS" << endl; 
    PARAMS.vparams.flag("AFLOWIN_FLAG::BANDS",vpflow.flag("PROTO_AFLOW::BANDS_CALCULATION"));
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.flag(\"AFLOWIN_FLAG::BANDS\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::BANDS") << endl;
 
    // check neglect_nomix
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: CHECK NEGLECT_NOMIX" << endl; 
    PARAMS.vparams.flag("AFLOWIN_FLAG::NEGLECT_NOMIX",vpflow.flag("PROTO_AFLOW::NEGLECT_NOMIX"));
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.flag(\"AFLOWIN_FLAG::NEGLECT_NOMIX\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::NEGLECT_NOMIX") << endl;
 
    // check stdout
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: CHECK STDOUT" << endl; 
    PARAMS.vparams.flag("AFLOWIN_FLAG::STDOUT",vpflow.flag("PROTO_AFLOW::STDOUT"));
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.flag(\"AFLOWIN_FLAG::STDOUT\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::STDOUT") << endl;
  
    // check module
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: CHECK MODULE" << endl;
    PARAMS.vparams.flag("AFLOWIN_FLAG::MODULE",vpflow.flag("PROTO_AFLOW::MODULE"));
    PARAMS.vparams.addattachedscheme("AFLOWIN_FLAG::MODULE",vpflow.getattachedscheme("PROTO_AFLOW::MODULE"),vpflow.flag("PROTO_AFLOW::MODULE"));
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.flag(\"AFLOWIN_FLAG::MODULE\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::MODULE") << endl;

    // check apl_supercell
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: CHECK APL_SUPERCELL" << endl;
    PARAMS.vparams.flag("AFLOWIN_FLAG::APL_SUPERCELL",vpflow.flag("PROTO_AFLOW::APL_SUPERCELL"));
    PARAMS.vparams.addattachedscheme("AFLOWIN_FLAG::APL_SUPERCELL",vpflow.getattachedscheme("PROTO_AFLOW::APL_SUPERCELL"),vpflow.flag("PROTO_AFLOW::APL_SUPERCELL"));
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.flag(\"AFLOWIN_FLAG::APL_SUPERCELL\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::APL_SUPERCELL") << endl;
    
    // check potim
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: CHECK POTIM" << endl;
    PARAMS.vparams.flag("AFLOWIN_FLAG::POTIM",vpflow.flag("PROTO_AFLOW::POTIM"));
    PARAMS.vparams.addattachedscheme("AFLOWIN_FLAG::POTIM",vpflow.getattachedscheme("PROTO_AFLOW::POTIM"),vpflow.flag("PROTO_AFLOW::POTIM"));
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.flag(\"AFLOWIN_FLAG::POTIM\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::POTIM") << endl;
  
    // check precision
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: CHECK PRECISION" << endl;
    PARAMS.vparams.flag("AFLOWIN_FLAG::PRECISION",vpflow.flag("PROTO_AFLOW::PRECISION"));
    PARAMS.vparams.addattachedscheme("AFLOWIN_FLAG::PRECISION",vpflow.getattachedscheme("PROTO_AFLOW::PRECISION"),vpflow.flag("PROTO_AFLOW::PRECISION"));
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.flag(\"AFLOWIN_FLAG::PRECISION\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::PRECISION") << endl;
  
    // check algorithm
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: CHECK ALGORITHM" << endl;
    PARAMS.vparams.flag("AFLOWIN_FLAG::ALGORITHM",vpflow.flag("PROTO_AFLOW::ALGORITHM"));
    PARAMS.vparams.addattachedscheme("AFLOWIN_FLAG::ALGORITHM",vpflow.getattachedscheme("PROTO_AFLOW::ALGORITHM"),vpflow.flag("PROTO_AFLOW::ALGORITHM"));
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.flag(\"AFLOWIN_FLAG::ALGORITHM\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::ALGORITHM") << endl;
  
    // check metagga
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: CHECK METAGGA" << endl;
    PARAMS.vparams.flag("AFLOWIN_FLAG::METAGGA",vpflow.flag("PROTO_AFLOW::METAGGA"));
    PARAMS.vparams.addattachedscheme("AFLOWIN_FLAG::METAGGA",vpflow.getattachedscheme("PROTO_AFLOW::METAGGA"),vpflow.flag("PROTO_AFLOW::METAGGA"));
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.flag(\"AFLOWIN_FLAG::METAGGA\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::METAGGA") << endl;

    // check ivdw
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: CHECK IVDW" << endl;
    PARAMS.vparams.flag("AFLOWIN_FLAG::IVDW",vpflow.flag("PROTO_AFLOW::IVDW"));
    PARAMS.vparams.addattachedscheme("AFLOWIN_FLAG::IVDW",vpflow.getattachedscheme("PROTO_AFLOW::IVDW"),vpflow.flag("PROTO_AFLOW::IVDW"));
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.flag(\"AFLOWIN_FLAG::IVDW\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::IVDW") << endl;

    // check relax_type
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: CHECK RELAX_TYPE" << endl;
    PARAMS.vparams.flag("AFLOWIN_FLAG::RELAX_TYPE",vpflow.flag("PROTO_AFLOW::RELAX_TYPE"));
    PARAMS.vparams.addattachedscheme("AFLOWIN_FLAG::RELAX_TYPE",vpflow.getattachedscheme("PROTO_AFLOW::RELAX_TYPE"),vpflow.flag("PROTO_AFLOW::RELAX_TYPE"));
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.flag(\"AFLOWIN_FLAG::RELAX_TYPE\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::RELAX_TYPE") << endl;
    
    // check relax_mode
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: CHECK RELAX_MODE" << endl;
    PARAMS.vparams.flag("AFLOWIN_FLAG::RELAX_MODE",vpflow.flag("PROTO_AFLOW::RELAX_MODE"));
    PARAMS.vparams.addattachedscheme("AFLOWIN_FLAG::RELAX_MODE",vpflow.getattachedscheme("PROTO_AFLOW::RELAX_MODE"),vpflow.flag("PROTO_AFLOW::RELAX_MODE"));
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.flag(\"AFLOWIN_FLAG::RELAX_MODE\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::RELAX_MODE") << endl;
  
    // check type
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: CHECK TYPE" << endl;
    PARAMS.vparams.flag("AFLOWIN_FLAG::TYPE",vpflow.flag("PROTO_AFLOW::TYPE"));
    PARAMS.vparams.addattachedscheme("AFLOWIN_FLAG::TYPE",vpflow.getattachedscheme("PROTO_AFLOW::TYPE"),vpflow.flag("PROTO_AFLOW::TYPE"));
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.flag(\"AFLOWIN_FLAG::TYPE\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::TYPE") << endl;
  
    // check CONVERT_UNIT_CELL
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: CHECK CONVERT_UNIT_CELL" << endl;
    PARAMS.vparams.flag("AFLOWIN_FLAG::CONVERT_UNIT_CELL",vpflow.flag("PROTO_AFLOW::CONVERT_UNIT_CELL"));
    PARAMS.vparams.addattachedscheme("AFLOWIN_FLAG::CONVERT_UNIT_CELL",vpflow.getattachedscheme("PROTO_AFLOW::CONVERT_UNIT_CELL"),vpflow.flag("PROTO_AFLOW::CONVERT_UNIT_CELL"));
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.flag(\"AFLOWIN_FLAG::CONVERT_UNIT_CELL\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::CONVERT_UNIT_CELL") << endl;
  
    // check volume_plus_equal
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: CHECK VOLUME_PLUS_EQUAL" << endl;
    PARAMS.vparams.flag("AFLOWIN_FLAG::VOLUME_PLUS_EQUAL",vpflow.flag("PROTO_AFLOW::VOLUME_PLUS_EQUAL"));
    PARAMS.vparams.addattachedscheme("AFLOWIN_FLAG::VOLUME_PLUS_EQUAL",vpflow.getattachedscheme("PROTO_AFLOW::VOLUME_PLUS_EQUAL"),vpflow.flag("PROTO_AFLOW::VOLUME_PLUS_EQUAL"));
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.flag(\"AFLOWIN_FLAG::VOLUME_PLUS_EQUAL\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::VOLUME_PLUS_EQUAL") << endl;
  
    // check volume_multiply_equal
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: CHECK VOLUME_MULTIPLY_EQUAL" << endl;
    PARAMS.vparams.flag("AFLOWIN_FLAG::VOLUME_MULTIPLY_EQUAL",vpflow.flag("PROTO_AFLOW::VOLUME_MULTIPLY_EQUAL"));
    PARAMS.vparams.addattachedscheme("AFLOWIN_FLAG::VOLUME_MULTIPLY_EQUAL",vpflow.getattachedscheme("PROTO_AFLOW::VOLUME_MULTIPLY_EQUAL"),vpflow.flag("PROTO_AFLOW::VOLUME_MULTIPLY_EQUAL"));
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.flag(\"AFLOWIN_FLAG::VOLUME_MULTIPLY_EQUAL\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::VOLUME_MULTIPLY_EQUAL") << endl;
    
    // check no_volume_adjustment
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: CHECK NO_VOLUME_ADJUSTMENT" << endl;
    PARAMS.vparams.flag("AFLOWIN_FLAG::NO_VOLUME_ADJUSTMENT",vpflow.flag("PROTO_AFLOW::NO_VOLUME_ADJUSTMENT"));
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.flag(\"AFLOWIN_FLAG::NO_VOLUME_ADJUSTMENT\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::NO_VOLUME_ADJUSTMENT") << endl;
  
    // check ediffg
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: CHECK EDIFFG" << endl;
    PARAMS.vparams.flag("AFLOWIN_FLAG::EDIFFG",vpflow.flag("PROTO_AFLOW::EDIFFG"));
    PARAMS.vparams.addattachedscheme("AFLOWIN_FLAG::EDIFFG",vpflow.getattachedscheme("PROTO_AFLOW::EDIFFG"),vpflow.flag("PROTO_AFLOW::EDIFFG"));
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.flag(\"AFLOWIN_FLAG::EDIFFG\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::EDIFFG") << endl;
  
    // check kppra
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: CHECK KPPRA" << endl;
    PARAMS.vparams.flag("AFLOWIN_FLAG::KPPRA",vpflow.flag("PROTO_AFLOW::KPPRA"));
    PARAMS.vparams.addattachedscheme("AFLOWIN_FLAG::KPPRA",vpflow.getattachedscheme("PROTO_AFLOW::KPPRA"),vpflow.flag("PROTO_AFLOW::KPPRA"));
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.flag(\"AFLOWIN_FLAG::KPPRA\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::KPPRA") << endl;
  
    // check kppra_static
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: CHECK KPPRA_STATIC" << endl;
    PARAMS.vparams.flag("AFLOWIN_FLAG::KPPRA_STATIC",vpflow.flag("PROTO_AFLOW::KPPRA_STATIC"));
    PARAMS.vparams.addattachedscheme("AFLOWIN_FLAG::KPPRA_STATIC",vpflow.getattachedscheme("PROTO_AFLOW::KPPRA_STATIC"),vpflow.flag("PROTO_AFLOW::KPPRA_STATIC"));
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.flag(\"AFLOWIN_FLAG::KPPRA_STATIC\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::KPPRA_STATIC") << endl;
 
    // check bands_grid
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: CHECK BANDS_GRID" << endl;
    PARAMS.vparams.flag("AFLOWIN_FLAG::BANDS_GRID",vpflow.flag("PROTO_AFLOW::BANDS_GRID"));
    PARAMS.vparams.addattachedscheme("AFLOWIN_FLAG::BANDS_GRID",vpflow.getattachedscheme("PROTO_AFLOW::BANDS_GRID"),vpflow.flag("PROTO_AFLOW::BANDS_GRID"));
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.flag(\"AFLOWIN_FLAG::BANDS_GRID\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::BANDS_GRID") << endl;
  
    // check enmax_multiply
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: CHECK ENMAX_MULTIPLY" << endl;
    PARAMS.vparams.flag("AFLOWIN_FLAG::ENMAX_MULTIPLY",vpflow.flag("PROTO_AFLOW::ENMAX_MULTIPLY"));
    PARAMS.vparams.addattachedscheme("AFLOWIN_FLAG::ENMAX_MULTIPLY",vpflow.getattachedscheme("PROTO_AFLOW::ENMAX_MULTIPLY"),vpflow.flag("PROTO_AFLOW::ENMAX_MULTIPLY"));
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.flag(\"AFLOWIN_FLAG::ENMAX_MULTIPLY\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::ENMAX_MULTIPLY") << endl;

    // check ABINIT QE VASP AIMS
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: CHECK VASP/QE/ABINIT/AIMS" << endl; 
    PARAMS.vparams.flag("AFLOWIN_FLAG::VASP",TRUE); // default
    PARAMS.vparams.flag("AFLOWIN_FLAG::QE",vpflow.flag("PROTO_AFLOW::QE"));
    PARAMS.vparams.flag("AFLOWIN_FLAG::ABINIT",vpflow.flag("PROTO_AFLOW::ABINIT"));
    PARAMS.vparams.flag("AFLOWIN_FLAG::AIMS",vpflow.flag("PROTO_AFLOW::AIMS"));
    if(PARAMS.vparams.flag("AFLOWIN_FLAG::AIMS")) {
      PARAMS.vparams.flag("AFLOWIN_FLAG::QE",FALSE);
      PARAMS.vparams.flag("AFLOWIN_FLAG::VASP",FALSE);
      PARAMS.vparams.flag("AFLOWIN_FLAG::ABINIT",FALSE);
    }
    if(PARAMS.vparams.flag("AFLOWIN_FLAG::ABINIT")) {
      PARAMS.vparams.flag("AFLOWIN_FLAG::QE",FALSE);
      PARAMS.vparams.flag("AFLOWIN_FLAG::VASP",FALSE);
      PARAMS.vparams.flag("AFLOWIN_FLAG::AIMS",FALSE);
    }
    if(PARAMS.vparams.flag("AFLOWIN_FLAG::QE"))     {
      PARAMS.vparams.flag("AFLOWIN_FLAG::VASP",FALSE);
      PARAMS.vparams.flag("AFLOWIN_FLAG::ABINIT",FALSE);
      PARAMS.vparams.flag("AFLOWIN_FLAG::AIMS",FALSE);
    }
    if(PARAMS.vparams.flag("AFLOWIN_FLAG::VASP"))   {
      PARAMS.vparams.flag("AFLOWIN_FLAG::ABINIT",FALSE);
      PARAMS.vparams.flag("AFLOWIN_FLAG::QE",FALSE);
      PARAMS.vparams.flag("AFLOWIN_FLAG::AIMS",FALSE);
    }
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.flag(\"AFLOWIN_FLAG::AIMS\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::AIMS") << endl;
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.flag(\"AFLOWIN_FLAG::ABINIT\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::ABINIT") << endl;
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.flag(\"AFLOWIN_FLAG::QE\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::QE") << endl;
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.flag(\"AFLOWIN_FLAG::VASP\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::VASP") << endl;
    
    // check stdout
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: CHECK STDOUT" << endl; 
    PARAMS.vparams.flag("AFLOWIN_FLAG::STDOUT",vpflow.flag("PROTO_AFLOW::STDOUT"));
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.flag(\"AFLOWIN_FLAG::STDOUT\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::STDOUT") << endl;
  
    // check list
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: CHECK LIST" << endl; 
    PARAMS.vparams.flag("AFLOWIN_FLAG::LIST",vpflow.flag("PROTO_AFLOW::LIST"));
    if(PARAMS.vparams.flag("AFLOWIN_FLAG::LIST")) PARAMS.vparams.addattachedscheme("AFLOWIN_FLAG::LIST_VCMD",vpflow.getattachedscheme("PROTO_AFLOW::LIST_VCMD"),PARAMS.vparams.flag("AFLOWIN_FLAG::LIST"));
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.flag(\"AFLOWIN_FLAG::LIST\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::LIST") << endl;
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.getattachedscheme((\"AFLOWIN_FLAG::LIST_VCMD\")=" << PARAMS.vparams.getattachedscheme("AFLOWIN_FLAG::LIST_VCMD") << endl;
    
    //DX 1/18/18 - Add ANRL functionality - START
    // check params
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: CHECK PARAMS" << endl;
    PARAMS.vparams.flag("AFLOWIN_FLAG::PARAMS",vpflow.flag("PARAMS"));
    PARAMS.vparams.addattachedscheme("AFLOWIN_FLAG::PARAMS",vpflow.getattachedscheme("PARAMS"),vpflow.flag("PARAMS"));
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.flag(\"AFLOWIN_FLAG::PARAMS\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::PARAMS") << endl;
    //DX 1/18/18 - Add ANRL functionality - END
    
    // check htqc
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: CHECK HTQC_ICSD" << endl; 
    PARAMS.vparams.flag("AFLOWIN_FLAG::HTQC_ICSD",vpflow.flag("PROTO_AFLOW::HTQC") || 
			(!aurostd::substring2bool(string_LABEL,"_ICSD_") &&
			 aurostd::substring2bool(string_LABEL,"ICSD_")));
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.flag(\"AFLOWIN_FLAG::HTQC_ICSD\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::HTQC_ICSD") << endl;
 
    // check pressure
    vector<string> params_pressure;
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: vpflow.getattachedscheme(\"PROTO_AFLOW::PRESSURE\")=" << vpflow.getattachedscheme("PROTO_AFLOW::PRESSURE") << endl;
    aurostd::string2tokens(vpflow.getattachedscheme("PROTO_AFLOW::PRESSURE"),params_pressure,",");
    PARAMS.vpressure.clear();
    if(params_pressure.size()!=0) for(uint i=0;i<params_pressure.size();i++) PARAMS.vpressure.push_back(aurostd::string2utype<double>(params_pressure.at(i)));
  
  
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: --- " << endl;
 
    uint nspeciesHTQC=aflowlib::PrototypeLibrariesSpeciesNumber(string_LABEL);
    bool alphabetic=RequestedAlphabeticLabeling(string_LABEL);
    PARAMS.vparams.addattachedscheme("AFLOWIN_STRING::LABEL",string_LABEL,TRUE);

    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: nspeciesHTQC=" << nspeciesHTQC << endl;
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: alphabetic=" << alphabetic << endl;
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vpressure.size()=" << PARAMS.vpressure.size() << endl;
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vkppra.size()=" << PARAMS.vkppra.size() << " " << endl;
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.getattachedscheme((\"AFLOWIN_STRING::LABEL\")=" << PARAMS.vparams.getattachedscheme("AFLOWIN_STRING::LABEL") << endl;  
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.getattachedscheme((\"AFLOWIN_STRING::POTENTIAL\")=" << PARAMS.vparams.getattachedscheme("AFLOWIN_STRING::POTENTIAL") << endl;  
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.flag(\"AFLOWIN_FLAG::MODULE\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::MODULE") << endl; //CO 180214
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.getattachedscheme(\"AFLOWIN_FLAG::MODULE\")=" << PARAMS.vparams.getattachedscheme("AFLOWIN_FLAG::MODULE") << endl; //CO 180214
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.flag(\"AFLOWIN_FLAG::APL_SUPERCELL\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::APL_SUPERCELL") << endl; //CO 180214
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.getattachedscheme(\"AFLOWIN_FLAG::APL_SUPERCELL\")=" << PARAMS.vparams.getattachedscheme("AFLOWIN_FLAG::APL_SUPERCELL") << endl; //CO 180214
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.flag(\"AFLOWIN_FLAG::POTIM\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::POTIM") << endl;
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.getattachedscheme(\"AFLOWIN_FLAG::POTIM\")=" << PARAMS.vparams.getattachedscheme("AFLOWIN_FLAG::POTIM") << endl;
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.flag(\"AFLOWIN_FLAG::PRECISION\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::PRECISION") << endl;
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.getattachedscheme(\"AFLOWIN_FLAG::PRECISION\")=" << PARAMS.vparams.getattachedscheme("AFLOWIN_FLAG::PRECISION") << endl;
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.flag(\"AFLOWIN_FLAG::ALGORITHM\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::ALGORITHM") << endl;
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.getattachedscheme(\"AFLOWIN_FLAG::ALGORITHM\")=" << PARAMS.vparams.getattachedscheme("AFLOWIN_FLAG::ALGORITHM") << endl;
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.flag(\"AFLOWIN_FLAG::METAGGA\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::METAGGA") << endl;
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.getattachedscheme(\"AFLOWIN_FLAG::METAGGA\")=" << PARAMS.vparams.getattachedscheme("AFLOWIN_FLAG::METAGGA") << endl;
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.flag(\"AFLOWIN_FLAG::IVDW\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::IVDW") << endl;
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.getattachedscheme(\"AFLOWIN_FLAG::IVDW\")=" << PARAMS.vparams.getattachedscheme("AFLOWIN_FLAG::IVDW") << endl;
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.flag(\"AFLOWIN_FLAG::RELAX_TYPE\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::RELAX_TYPE") << endl; //CO 180214
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.getattachedscheme(\"AFLOWIN_FLAG::RELAX_TYPE\")=" << PARAMS.vparams.getattachedscheme("AFLOWIN_FLAG::RELAX_TYPE") << endl; //CO 180214
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.flag(\"AFLOWIN_FLAG::RELAX_MODE\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::RELAX_MODE") << endl;
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.getattachedscheme(\"AFLOWIN_FLAG::RELAX_MODE\")=" << PARAMS.vparams.getattachedscheme("AFLOWIN_FLAG::RELAX_MODE") << endl;
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.flag(\"AFLOWIN_FLAG::TYPE\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::TYPE") << endl;
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.getattachedscheme(\"AFLOWIN_FLAG::TYPE\")=" << PARAMS.vparams.getattachedscheme("AFLOWIN_FLAG::TYPE") << endl;
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.flag(\"AFLOWIN_FLAG::CONVERT_UNIT_CELL\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::CONVERT_UNIT_CELL") << endl;
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.getattachedscheme(\"AFLOWIN_FLAG::CONVERT_UNIT_CELL\")=" << PARAMS.vparams.getattachedscheme("AFLOWIN_FLAG::CONVERT_UNIT_CELL") << endl;
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.flag(\"AFLOWIN_FLAG::VOLUME_PLUS_EQUAL\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::VOLUME_PLUS_EQUAL") << endl;
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.getattachedscheme(\"AFLOWIN_FLAG::VOLUME_PLUS_EQUAL\")=" << PARAMS.vparams.getattachedscheme("AFLOWIN_FLAG::VOLUME_PLUS_EQUAL") << endl;
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.flag(\"AFLOWIN_FLAG::VOLUME_MULTIPLY_EQUAL\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::VOLUME_MULTIPLY_EQUAL") << endl;
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.getattachedscheme(\"AFLOWIN_FLAG::VOLUME_MULTIPLY_EQUAL\")=" << PARAMS.vparams.getattachedscheme("AFLOWIN_FLAG::VOLUME_MULTIPLY_EQUAL") << endl;
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.flag(\"AFLOWIN_FLAG::NO_VOLUME_ADJUSTMENT\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::NO_VOLUME_ADJUSTMENT") << endl; //CO 180214
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.flag(\"AFLOWIN_FLAG::EDIFFG\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::EDIFFG") << endl;
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.getattachedscheme(\"AFLOWIN_FLAG::EDIFFG\")=" << PARAMS.vparams.getattachedscheme("AFLOWIN_FLAG::EDIFFG") << endl;
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.flag(\"AFLOWIN_FLAG::KPPRA\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::KPPRA") << endl;
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.getattachedscheme(\"AFLOWIN_FLAG::KPPRA\")=" << PARAMS.vparams.getattachedscheme("AFLOWIN_FLAG::KPPRA") << endl;
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.flag(\"AFLOWIN_FLAG::KPPRA_STATIC\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::KPPRA_STATIC") << endl;
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.getattachedscheme(\"AFLOWIN_FLAG::KPPRA_STATIC\")=" << PARAMS.vparams.getattachedscheme("AFLOWIN_FLAG::KPPRA_STATIC") << endl;
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.flag(\"AFLOWIN_FLAG::BANDS_GRID\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::BANDS_GRID") << endl;
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.getattachedscheme(\"AFLOWIN_FLAG::BANDS_GRID\")=" << PARAMS.vparams.getattachedscheme("AFLOWIN_FLAG::BANDS_GRID") << endl;
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.flag(\"AFLOWIN_FLAG::ENMAX_MULTIPLY\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::ENMAX_MULTIPLY") << endl;
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.getattachedscheme(\"AFLOWIN_FLAG::ENMAX_MULTIPLY\")=" << PARAMS.vparams.getattachedscheme("AFLOWIN_FLAG::ENMAX_MULTIPLY") << endl;
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.flag(\"AFLOWIN_FLAG::USAGE\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::USAGE") << endl;
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.flag(\"AFLOWIN_FLAG::REVERSE\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::REVERSE") << endl;
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.flag(\"AFLOWIN_FLAG::MISSING\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::MISSING") << endl;
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.flag(\"AFLOWIN_FLAG::NOAUTOPP\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::NOAUTOPP") << endl;
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.flag(\"AFLOWIN_FLAG::BANDS\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::BANDS") << endl;
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.flag(\"AFLOWIN_FLAG::AUTOLDAU\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::AUTOLDAU") << endl;
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.flag(\"AFLOWIN_FLAG::AUTONOLDAU\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::AUTONOLDAU") << endl;
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.flag(\"AFLOWIN_FLAG::STDOUT\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::STDOUT") << endl;
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.flag(\"AFLOWIN_FLAG::LIST\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::LIST") << endl;
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.flag(\"AFLOWIN_FLAG::PARAMS\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::PARAMS") << endl; //DX 1/18/18 - Added ANRL functionality
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.getattachedscheme((\"AFLOWIN_FLAG::LIST_VCMD\")=" << PARAMS.vparams.getattachedscheme("AFLOWIN_FLAG::LIST_VCMD") << endl;
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.flag(\"AFLOWIN_FLAG::HTQC_ICSD\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::HTQC_ICSD") << endl;
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.flag(\"AFLOWIN_FLAG::NEGLECT_NOMIX\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::NEGLECT_NOMIX") << endl;
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.flag(\"AFLOWIN_FLAG::VASP\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::VASP") << endl;
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.flag(\"AFLOWIN_FLAG::ABINIT\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::ABINIT") << endl;
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.flag(\"AFLOWIN_FLAG::AIMS\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::AIMS") << endl;
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: PARAMS.vparams.flag(\"AFLOWIN_FLAG::QE\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::QE") << endl;

    // if(LDEBUG) exit(0);
  
    if(PARAMS.ucell.size()==1 && aurostd::substring2bool(PARAMS.vparams.getattachedscheme("AFLOWIN_STRING::LABEL"),"_ICSD_") && !PARAMS.vparams.flag("AFLOWIN_FLAG::HTQC_ICSD")) {
      if(LDEBUG) cerr << "pflow::PROTO_AFLOW: running AVASP_MakePrototypeICSD_AFLOWIN" << endl;
      //  return AVASP_MakePrototypeICSD_AFLOWIN(PARAMS.ucell,FALSE);
      return AVASP_MakePrototypeICSD_AFLOWIN((_AVASP_PROTO*)&PARAMS,FALSE);
    } else {
      // if(alphabetic==TRUE && (nspeciesHTQC==3 || nspeciesHTQC==4)) AlphabetizePrototypeLabelSpeciesArgv(argv); // it is going to fail if you want a big list of protos
      if(LDEBUG) cerr << "pflow::PROTO_AFLOW: running AVASP_MakePrototype_AFLOWIN" << endl;
      return AVASP_MakePrototype_AFLOWIN((_AVASP_PROTO*)&PARAMS);
    }
    if(LDEBUG) cerr << "pflow::PROTO_AFLOW: END" << endl; 
    return FALSE;
  }
} // namespace pflow

// ***************************************************************************
// pflow::STATDIEL // CAMILO
// ***************************************************************************
namespace pflow {
  void STATDIEL(vector<string>& argv) { // loop pflow::STATDIEL
    string outcar ;
    xvector<double> real(3), imag(3) ;
    if(argv.size() != 3) { // user control lines - aflow specific.
      init::ErrorOption(cout,"","pflow::STATDIEL","aflow --statdiel OUTCAR*");
      exit(0) ;
    }
    outcar = argv.at(2) ;
    KBIN::GetStatDiel(outcar, real, imag) ;
    printf ("%11.6f %11.6f %11.6f "  ,real(1),real(2),real(3)) ;
    printf ("%11.6f %11.6f %11.6f \n",imag(1),imag(2),imag(3)) ;
  }  // loop pflow::STATDIEL
}

// ***************************************************************************
// pflow::QMVASP
// ***************************************************************************
namespace pflow {
  bool QMVASP(vector<string> argv) {
    string directory="./";
    vector<string> vCARS,vCARS_relax;
    bool found=FALSE;
    vCARS.push_back("CONTCAR");vCARS.push_back("OSZICAR");vCARS.push_back("OUTCAR");
    if(!found&&argv.size()<2) found=TRUE; // nothing to do
    if(!found&&argv.size()==2) {  // DEFAULT
      vCARS_relax.push_back("CONTCAR.relax2");vCARS_relax.push_back("OSZICAR.relax2");vCARS_relax.push_back("OUTCAR.relax2");
      found=TRUE;
    }
    if(!found&&argv.size()==3 && aurostd::args2flag(argv,"--static")) {  // DEFAULT
      vCARS_relax.push_back("CONTCAR.static");vCARS_relax.push_back("OSZICAR.static");vCARS_relax.push_back("OUTCAR.static");
      found=TRUE;
    }
    // [OBSOLETE] if(!found&&argv.size()==4 && aurostd::args2flag(argv,"--DIRECTORY|--D|--d|./")) {  // DEFAULT
    if(!found&&argv.size()==4 && XHOST.vflag_control.flag("DIRECTORY")) {  // DEFAULT
      // [OBSOLETE]     directory=aurostd::args2string(argv,"--DIRECTORY|--D|--d","./")+"/";    //   cerr << directory << endl;
      directory=XHOST.vflag_control.getattachedscheme("DIRECTORY");
      vCARS_relax.push_back("CONTCAR.relax2");vCARS_relax.push_back("OSZICAR.relax2");vCARS_relax.push_back("OUTCAR.relax2");
      found=TRUE;
    }
    // [OBSOLETE] if(!found&&argv.size()==5 && aurostd::args2flag(argv,"--static") && aurostd::args2flag(argv,"--DIRECTORY|--D|--d|./")) {  // DEFAULT
    if(!found&&argv.size()==5 && aurostd::args2flag(argv,"--static") && XHOST.vflag_control.flag("DIRECTORY")) {  // DEFAULT
      // [OBSOLETE]    directory=aurostd::args2string(argv,"--DIRECTORY|--D|--d","./")+"/";    //   cerr << directory << endl;
      directory=XHOST.vflag_control.getattachedscheme("DIRECTORY");
      vCARS_relax.push_back("CONTCAR.static");vCARS_relax.push_back("OSZICAR.static");vCARS_relax.push_back("OUTCAR.static");
      found=TRUE;
    }
    if(!found) return FALSE; // error
    for(uint i=0;i<vCARS_relax.size();i++) {
      //      cerr <<  vCARS_relax.at(i) << endl;
      if(aurostd::FileExist(directory+"/"+vCARS_relax.at(i)+".bz2")) aurostd::execute(XHOST.command("bzcat")+" "+directory+"/"+vCARS_relax.at(i)+".bz2 > "+directory+"/"+vCARS.at(i));
      if(aurostd::FileExist(directory+"/"+vCARS_relax.at(i)+".gz")) aurostd::execute("zcat "+directory+"/"+vCARS_relax.at(i)+".gz > "+directory+"/"+vCARS.at(i));
      if(aurostd::FileExist(directory+"/"+vCARS_relax.at(i))) aurostd::execute("cat "+directory+"/"+vCARS_relax.at(i)+" > "+directory+"/"+vCARS.at(i));
    }
    cout << "pflow::QMVASP Performing: " << directory << endl;
    _xvasp xvasp;
    xvasp.Directory=directory;
    // [OBSOLETE]  xvasp.str=xstructure(directory+"/CONTCAR",IOVASP_POSCAR);
    xvasp.str=xstructure(directory+"/CONTCAR",IOAFLOW_AUTO);
    //  cout <<
    KBIN::VASP_Analyze(xvasp,TRUE);
    for(uint i=0;i<vCARS.size();i++) aurostd::execute("rm -f "+directory+"/"+vCARS.at(i));
    if(aurostd::FileExist(directory+"/"+vCARS_relax.at(0)+".bz2")) {
      aurostd::execute("rm -f "+xvasp.Directory+"/"+_AFLOW_QMVASP_FILE_+".bz2");    
      aurostd::execute(XHOST.command("bzip2")+" "+xvasp.Directory+"/"+_AFLOW_QMVASP_FILE_);
    }
    if(aurostd::FileExist(directory+"/"+vCARS_relax.at(0)+".gz")) {
      aurostd::execute("rm -f"+xvasp.Directory+"/"+_AFLOW_QMVASP_FILE_+".gz");
      aurostd::execute(XHOST.command("gzip")+" "+xvasp.Directory+"/"+_AFLOW_QMVASP_FILE_);
    }
    //  cout << a << endl;
    return TRUE;
  }
} // namespace pflow

// ***************************************************************************
// pflow::RAYTRACE
// ***************************************************************************
namespace pflow {
  void RAYTRACE(vector<string> argv) {
    // cout << aflow::Banner("BANNER_TINY") << endl;
    pflow::RayTraceManager(argv); // Manage ray tracing related activities.
  }
} // namespace pflow

// ***************************************************************************
// pflow::RASMOL
// ***************************************************************************
namespace pflow {
  void RASMOL(string options,istream& input) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "pflow::RASMOL: BEGIN" << endl;
    xvector<int> ijk(3);           // default 
    ijk[1]=1;ijk[2]=1;ijk[3]=1;    // default
    vector<string> tokens;
    aurostd::string2tokens(options,tokens,",");

    if(tokens.size()>3) {
      init::ErrorOption(cout,options,"pflow::RASMOL","aflow --rasmol[=n1[,n2[,n3]]] < POSCAR");
      exit(0);
    }
    
    if(tokens.size()>=1) ijk[1]=aurostd::string2utype<int>(tokens.at(0)); 
    if(tokens.size()>=2) ijk[2]=aurostd::string2utype<int>(tokens.at(1)); 
    if(tokens.size()>=3) ijk[3]=aurostd::string2utype<int>(tokens.at(2)); 
    
    ofstream FileOUTPUT;
    string FileOUTPUTName="aflow.rasmol.xyz."+strPID();
    FileOUTPUT.open(FileOUTPUTName.c_str(),std::ios::out);
    if(FileOUTPUT) {
      xvector<int> _ijk(vabs(ijk));
      if(max(_ijk)==0) _ijk.set(1);
      ostringstream aus_exec;
      xstructure a(input,IOAFLOW_AUTO);
      PrintXYZ(a,_ijk,FileOUTPUT);
      FileOUTPUT.close();
      //   cerr << _ijk << endl;
      aus_exec << XHOST.command("rasmol") << " -xyz " << FileOUTPUTName << endl;// " && rm -f " << FileOUTPUTName << endl;
      aurostd::execute(aus_exec);
      // aurostd::StringstreamClean(aus_exec);
      aus_exec << "rm -f " <<  FileOUTPUTName << endl;
      aurostd::execute(aus_exec);
    } else FileOUTPUT.close();
  }
} // namespace pflow

// ***************************************************************************
// pflow::RSM
// ***************************************************************************
namespace pflow {
  void RSM(vector<string> argv, istream& input) {
  
    int i,j,Ntypes=0;
    string strtmp="",AtomSymbol="";
    bool Z_flag=false;
    int Z[20]; //to store the z values following --z option

    //FINDING -Z
    int iargv=0;
    Ntypes=0;
    for(i=1;i<(int) argv.size();i++) {
      strtmp=argv.at(i);
      if(strtmp=="--Z" or strtmp=="--z") {
	Z_flag=TRUE; iargv=i+1;
	break;
      }
    }
    Ntypes=(int)argv.size()-iargv;
    cerr<<"Ntypes = "<<Ntypes<<endl;
    for(j=0;j<Ntypes;j++) {
      AtomSymbol=argv.at(j+iargv);
      cerr<<AtomSymbol<<"(";
      if(AtomSymbol[0]<'0' || AtomSymbol[0]>'9')   Z[j]=(int)GetAtomNumber(AtomSymbol);
      else Z[j]=aurostd::string2utype<int>(AtomSymbol);//index in Z starts from ZERO
      cerr<<Z[j]<<")  ";
    }
    cerr<<endl;

    xstructure a(input,IOAFLOW_AUTO);

    if(Ntypes<1) {Z_flag=false; Ntypes=a.num_each_type.size();}
    if((int) a.num_each_type.size()!=Ntypes) {
      cerr<<"STOP: inconsistency in number of types in POSCAR and --z"<<endl
	  <<"a.num_each_type.size() = "<<a.num_each_type.size()<<endl
	  <<"Ntypes = "<<Ntypes<<endl;
      abort();
    }
    if(!Z_flag)
      for(i=0;i<Ntypes;i++) Z[i]=i+1;

    //fixing atoms.type
    int k,ptr=0;
    for(i=0;i<Ntypes;i++) {
      for(j=0;j<a.num_each_type.at(i);j++) {
	ptr=0;
	for(k=0;k<i+1;k++)
	  ptr=ptr+a.num_each_type.at(k);
	ptr=ptr-a.num_each_type.at(i)+j;//        ptr=sum(species(1:i))-species(i)+j;
	a.atoms.at(ptr).type=Z[i];
      }
    }

    int num_atoms=0;
    for(i=0;i<(int) a.num_each_type.size();i++)
      num_atoms+=a.num_each_type.at(i);

    PrintRSM(a,cout);
  }
} // namespace pflow

// ***************************************************************************
// pflow::RBANAL
// ***************************************************************************
namespace pflow {
  void RBANAL(vector<string> argv) {
    // cout << aflow::Banner("BANNER_TINY") << endl;
    // Read in input files.
    int nim=atoi(argv.at(2).c_str());
    string path_flag(argv.at(3));
    PrintRBAnal(nim,path_flag,cout);
  }
} // namespace pflow

// ***************************************************************************
// pflow::RBDIST
// ***************************************************************************
namespace pflow {
  void RBDIST(vector<string> argv) {
    xstructure strA,strB;
    ifstream infileA(argv.at(2).c_str());
    ifstream infileB(argv.at(3).c_str());
    infileA >> strA;
    infileB >> strB;
    string path_flag(argv.at(4).c_str());
    double totdist;
    pflow::matrix<double> cm(2,3,999999);
    xstructure diffstr;
    pflow::RBPoscarDisp(strA,strB,diffstr,totdist,cm,path_flag);
    PrintRBPoscarDisp(diffstr,totdist,cm,path_flag,cout);
  }
} // namespace pflow

// ***************************************************************************
// pflow::RDF
// ***************************************************************************
namespace pflow {
  void RDF(string options,istream& input) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "pflow::RDF: BEGIN" << endl;
    vector<string> tokens;
    aurostd::string2tokens(options,tokens,",");
    if(tokens.size()>3) {
      init::ErrorOption(cout,options,"pflow::RDF","aflow --rdf[=rmax[,nbins[,sigma]]] < POSCAR");
      exit(0);
    } 

    xstructure a(input,IOAFLOW_AUTO);
    double rmax=(double) 5.0;
    int nbins=(int) 25;
    int smooth_width=(int) 0;
    if(LDEBUG) cerr << "pflow::RDF: tokens.size()=" << tokens.size() << endl;
    if(tokens.size()>=1) rmax=aurostd::string2utype<double>(tokens.at(0));
    if(tokens.size()>=1) cerr << "pflow::RDF: tokens.at(0)=" << tokens.at(0) << endl;
    if(tokens.size()>=2) nbins=aurostd::string2utype<int>(tokens.at(1));
    if(tokens.size()>=3) smooth_width=aurostd::string2utype<int>(tokens.at(2));
    
    pflow::matrix<double> rdf_all;
    cerr << "[1]" << endl;
    pflow::GetRDF(a,rmax,nbins,rdf_all);
    cerr << "[1]" << endl;
    // for(int i=0;i<rdf_all.size();i++) {for(int j=0;j<rdf_all[i].size();j++) cerr << rdf_all[i][j] << " "; cerr << endl;};exit(0);
    pflow::matrix<double> rdf_all_sm;
    rdf_all_sm=pflow::GetSmoothRDF(rdf_all,smooth_width);
    pflow::matrix<double> rdfsh_all;
    pflow::matrix<double> rdfsh_loc; // Radial location of rdf shells.
    pflow::GetRDFShells(a,rmax,nbins,smooth_width,rdf_all_sm,rdfsh_all,rdfsh_loc);
    PrintRDF(a,rmax,nbins,smooth_width,rdf_all_sm,rdfsh_all,rdfsh_loc,cout);
    if(LDEBUG) cerr << "pflow::RDF: END" << endl;
  }
} // namespace pflow

// ***************************************************************************
// pflow::RDFCMP
// ***************************************************************************
namespace pflow {
  void RDFCMP(string options) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "pflow::RDFCMP: BEGIN" << endl;
    vector<string> tokens;
    aurostd::string2tokens(options,tokens,",");
    if(tokens.size()!=6) {
      init::ErrorOption(cout,options,"pflow::RDFCMP","aflow --rdfcmp=rmax,nbins,sigma,nshmax,POSCAR1,POSCAR2");
      exit(0);
    } 
    double rmax=(double) aurostd::string2utype<double>(tokens.at(0));
    int nbins=(int) aurostd::string2utype<int>(tokens.at(1));
    int smooth_width=(int) aurostd::string2utype<int>(tokens.at(2));
    int nsh=(int) aurostd::string2utype<int>(tokens.at(3));
    
    if(!aurostd::FileExist(tokens.at(4))) {
      cerr << "ERROR - pflow::RDFCMP: file1=" << tokens.at(4) << " does not exist..." << endl;
      exit(0);
    }
    if(!aurostd::FileExist(tokens.at(5))) {
      cerr << "ERROR - pflow::RDFCMP: file2=" << tokens.at(5) << " does not exist..." << endl;
      exit(0);
    }
    xstructure strA(tokens.at(4),IOAFLOW_AUTO);
    xstructure strB(tokens.at(5),IOAFLOW_AUTO);
    // Get rdfs
    pflow::matrix<double> rdf_all_A;
    pflow::matrix<double> rdf_all_B;
    pflow::GetRDF(strA,rmax,nbins,rdf_all_A);
    pflow::GetRDF(strB,rmax,nbins,rdf_all_B);
    pflow::matrix<double> rdf_all_A_sm=pflow::GetSmoothRDF(rdf_all_A,smooth_width);
    pflow::matrix<double> rdf_all_B_sm=pflow::GetSmoothRDF(rdf_all_B,smooth_width);
    // Get shells
    pflow::matrix<double> rdfsh_all_A;
    pflow::matrix<double> rdfsh_loc_A; // Radial location of rdf shells.
    pflow::matrix<double> rdfsh_all_B;
    pflow::matrix<double> rdfsh_loc_B; // Radial location of rdf shells.
    pflow::GetRDFShells(strA,rmax,nbins,smooth_width,rdf_all_A_sm,rdfsh_all_A,rdfsh_loc_A);
    pflow::GetRDFShells(strB,rmax,nbins,smooth_width,rdf_all_B_sm,rdfsh_all_B,rdfsh_loc_B);
    vector<int> best_match;
    pflow::matrix<double> rms_mat;
    pflow::CmpRDFShells(strA,strB,rdfsh_all_A,rdfsh_all_B,nsh,best_match,rms_mat);
    PrintRDFCmp(strA,strB,rmax,nbins,smooth_width,nsh,rdfsh_all_A,rdfsh_all_B,best_match,rms_mat,cout);
    if(LDEBUG) cerr << "pflow::RDFCMP: END" << endl;
    // exit(1);
  }
} // namespace pflow

// ***************************************************************************
// pflow::RMATOM
// ***************************************************************************
namespace pflow {
  xstructure RMATOM(istream& input,const int& iatom) {
    xstructure a(input,IOAFLOW_AUTO);
    a.RemoveAtom(iatom);
    return a;
  }
} // namespace pflow

// ***************************************************************************
// pflow::RMCOPIES
// ***************************************************************************
namespace pflow {
  xstructure RMCOPIES(istream& input) {
    xstructure a(input,IOAFLOW_AUTO);
    a.RemoveCopies();
    return a;
  }
} // namespace pflow

// ***************************************************************************
// pflow::SCALE
// ***************************************************************************
namespace pflow {
  xstructure SCALE(string options,istream& input) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "pflow::SCALE: BEGIN" << endl;
    vector<string> tokens;
    aurostd::string2tokens(options,tokens,",");
    if(tokens.size()!=1) {
      init::ErrorOption(cout,options,"pflow::SCALE","aflow --scale=s < POSCAR");
      exit(0);
    } 
    // move on
    xstructure a(input,IOAFLOW_AUTO);
    xstructure b;b=a;
    double scale=0.0; // some defaults
    if(tokens.size()>=1) scale=aurostd::string2utype<double>(tokens.at(0));
    
    if(LDEBUG) cerr << "DEBUG  b.scale=" << b.scale << endl;
    b=ReScale(b,scale);
    if(LDEBUG) cerr << "DEBUG  b.scale=" << b.scale << endl;
    b.neg_scale=FALSE;
    if(LDEBUG) cerr << "pflow::SCALE: END" << endl;
    return b;
  }
} // namespace pflow

// ***************************************************************************
// pflow::INFLATE_LATTICE
// ***************************************************************************
namespace pflow {
  xstructure INFLATE_LATTICE(string options,istream& input) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    vector<string> tokens;
    aurostd::string2tokens(options,tokens,",");
    if(tokens.size()!=1) {
      init::ErrorOption(cout,options,"pflow::INFLATE_LATTICE","aflow --inflate_lattice=coefficient | --ilattice=coefficient < POSCAR");
      exit(0);
    } 
    // move on
    xstructure a(input,IOAFLOW_AUTO);
    xstructure b;b=a;
    double coefficient=0.0; // some defaults
    if(tokens.size()>=1) coefficient=aurostd::string2utype<double>(tokens.at(0));
    if(LDEBUG) cerr << "DEBUG  b.scale=" << b.scale << endl;
    b=InflateLattice(b,coefficient);
    if(LDEBUG) cerr << "DEBUG  b.scale=" << b.scale << endl;
    b.neg_scale=FALSE;
    if(LDEBUG) cerr << "pflow::INFLATE_LATTICE: END" << endl;
    return b;
  }
} // namespace pflow

// ***************************************************************************
// pflow::INFLATE_VOLUME
// ***************************************************************************
namespace pflow {
  xstructure INFLATE_VOLUME(string options,istream& input) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    vector<string> tokens;
    aurostd::string2tokens(options,tokens,",");
    if(tokens.size()!=1) {
      init::ErrorOption(cout,options,"pflow::INFLATE_VOLUME","aflow --inflate_volume=coefficient | --ivolume=coefficient < POSCAR");
      exit(0);
    } 
    // move on
    xstructure a(input,IOAFLOW_AUTO);
    xstructure b;b=a;
    double coefficient=0.0; // some defaults
    if(tokens.size()>=1) coefficient=aurostd::string2utype<double>(tokens.at(0));
    if(LDEBUG) cerr << "DEBUG  b.scale=" << b.scale << endl;
    b=InflateVolume(b,coefficient);
    if(LDEBUG) cerr << "DEBUG  b.scale=" << b.scale << endl;
    b.neg_scale=FALSE;
    if(LDEBUG) cerr << "pflow::INFLATE_VOLUME: END" << endl;
    return b;
  }
} // namespace pflow

// ***************************************************************************
// pflow::SD
// ***************************************************************************
namespace pflow {
  xstructure SD(vector<string> argv,istream& input) {
    xstructure a(input,IOAFLOW_AUTO);
    xstructure b;b=a;
    // Read in input file.
    // Set a sd vector
    int num_types=a.num_each_type.size();
    vector<string> sd(num_types);
    for(int i=0;i<num_types;i++) sd[i]="TTT";
    //  vector<string> sd(argv.size()-2);
    for(int i=0;i<min(((int) argv.size()-2),num_types);i++) {
      sd[i]=string(argv.at(i+2));
    }
    b.isd=TRUE;
    b=SetSDTypes(b,sd);
    return b;
  }
} // namespace pflow

// ***************************************************************************
// pflow::SETCM
// ***************************************************************************
namespace pflow {
  xstructure SETCM(istream& input,const xvector<double>& newcm) {
    xstructure str(input,IOAFLOW_AUTO);
    xvector<double> oldcm(GetMom1(str));
    str=ShiftCPos(str,(newcm-oldcm)/str.scale);
    return str;
  }
} // namespace pflow

// ***************************************************************************
// pflow::SETORIGIN
// ***************************************************************************
namespace pflow {
  xstructure SETORIGIN(istream& input,const xvector<double>& neworigin) {
    xstructure str(input,IOAFLOW_AUTO);
    // xvector<double> oldorigin(GetMom1(str));
    xvector<double> oldorigin(3);oldorigin.clear();
    //  str=ShiftCPos(str,(neworigin-oldorigin)/str.scale);
    if(str.coord_flag==_COORDS_CARTESIAN_)  str=ShiftCPos(str,-(neworigin-oldorigin));
    if(str.coord_flag==_COORDS_FRACTIONAL_) str=ShiftFPos(str,-(neworigin-oldorigin));
    //    str=pflow::SetOrigin(str,neworigin);
    return str;
  }
} // namespace pflow

namespace pflow {
  xstructure SETORIGIN(istream& input,const int& natom) {
    xstructure str(input,IOAFLOW_AUTO);
    if(natom<0) return str;
    if(natom<(int)str.atoms.size()) {
      if(str.coord_flag==_COORDS_CARTESIAN_)  str=ShiftCPos(str,-str.atoms.at((uint) natom).cpos);
      if(str.coord_flag==_COORDS_FRACTIONAL_) str=ShiftFPos(str,-str.atoms.at((uint) natom).fpos);
    }
    return str;
  }
} // namespace pflow

// ***************************************************************************
// pflow::SEWALD
// ***************************************************************************
namespace pflow {
  void SEWALD(vector<string> argv,istream& input) {
    // cout << aflow::Banner("BANNER_TINY") << endl;
    double Ks=atof(argv.at(2).c_str());
    double SUMTOL=1.0e-16;
    xstructure str(input,IOAFLOW_AUTO);
    str = GetNiggliStr(str);
    double epoint=0.0,ereal=0.0,erecip=0.0,eewald=0.0;
    ereal=pflow::ScreenedESEner(str,Ks,SUMTOL);
    eewald=ereal;
    pflow::PrintEwald(str,epoint,ereal,erecip,eewald,Ks,SUMTOL,cout);
  }
} // namespace pflow

// ***************************************************************************
// pflow::SG
// ***************************************************************************
namespace pflow {
  //DX 9/21/17 [OBSOLETE] string SG(string options,istream& input,string mode,string print) {
  string SG(aurostd::xoption& vpflow,istream& input,string mode,string print) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string flag_name = "SG::"+mode; //DX 9/26/17
    if(print == "LABEL" || print == "NUMBER"){
      flag_name += "_" + print;
    }
    string options = vpflow.getattachedscheme(flag_name); //DX 9/26/17
    if(LDEBUG) cerr << "pflow::SG: mode=" << mode << endl;
 
    vector<string> tokens;
    aurostd::string2tokens(options,tokens,",");
    if(tokens.size()==1) {
      if(tokens.at(0)=="usage" || tokens.at(0)=="USAGE") {
	init::ErrorOption(cout,options,"pflow::SG",
                          aurostd::liststring2string("aflow --aflowSG[=tolerance| =tight| =loose] [--mag|--magnetic|--magmom=[m1,m2,...|INCAR|OUTCAR]] < POSCAR  default: (minimum_interatomic_distance)/100.0",
						     "aflow --platonSG[_label,_number][=EQUAL| EXACT][,ang,d1,d2,d3] < POSCAR  default:"+
						     string("EQUAL=")+aurostd::utype2string<int>(DEFAULT_PLATON_P_EQUAL)+","+
						     string("EXACT=")+aurostd::utype2string<int>(DEFAULT_PLATON_P_EXACT)+","+
						     aurostd::utype2string(DEFAULT_PLATON_P_ANG,5)+","+
						     aurostd::utype2string(DEFAULT_PLATON_P_D1,5)+","+
						     aurostd::utype2string(DEFAULT_PLATON_P_D2,5)+","+
						     aurostd::utype2string(DEFAULT_PLATON_P_D3,5)+"",
						     "aflow --findsymSG[_label,_number][=tolerance] < POSCAR   default:"+aurostd::utype2string(DEFAULT_FINDSYM_TOL,5)));
	exit(0);
      } 
    }
    // move on
    
    if(LDEBUG) cerr << "pflow::SG: tokens.size()=" << tokens.size() << endl;
    
    // check usage
    //   if(LDEBUG) cerr << "pflow::SG: vpflow.getattachedscheme(\"SG::USAGE\")=" << vpflow.flag("SG::USAGE") << endl;
    
    // [OBSOLETE] xstructure a(input,IOVASP_POSCAR);
    if(input.peek() == EOF) {
      cerr << "File is empty. Check POSCAR." << endl;
      exit(0);
    }
    xstructure a(input,IOAFLOW_AUTO);
    //DX 1/24/18 [OBSOLETE] a.is_vasp4_poscar_format=TRUE; a.is_vasp5_poscar_format=FALSE;
    //   cerr << a << endl; exit(0);
    // AFLOW ENGINE RHT
    if(mode=="AFLOW" || mode=="aflow") { // RHT
      if(LDEBUG) cerr << "pflow::SG: aflow" << endl;
      // DX - START
      a.ReScale(1.0);
      //DX 9/21/17 - MAGNETIC SYMMETRY - START
      if(vpflow.flag("SG::MAGNETIC")){
	string magmom_info = vpflow.getattachedscheme("SG::MAGNETIC");
        int num_atoms=a.atoms.size();
        bool is_noncoll=false; 
        vector<xvector<double> > vmag_noncoll;                                          //DX 12/5/17 - added non-collinear
        if(GetNonCollinearMagneticInfo(num_atoms,magmom_info,vmag_noncoll)){            //DX 12/5/17 - added non-collinear
          if(LDEBUG){cerr << "pflow::SG: Non-collinear spin system detected." << endl;} //DX 12/5/17 - added non-collinear
          is_noncoll = true;
          if(!AddSpinToXstructure(a,vmag_noncoll)){
            exit(0);
          } 
        }                                                                               //DX 12/5/17 - added non-collinear
        bool is_coll=false; 
        vector<double> vmag;
        if(!is_noncoll){
          if(GetCollinearMagneticInfo(num_atoms,magmom_info,vmag)){
            if(LDEBUG){cerr << "pflow::SG: Collinear spin system detected." << endl;}   //DX 12/5/17 - added non-collinear
            is_coll = true;
            if(!AddSpinToXstructure(a,vmag)){
               exit(0);
            } 
          } 
        }
        if(!is_noncoll && !is_coll){                                                                                 //DX 12/5/17 - added non-collinear
          cerr << "pflow::SG: ERROR: Could not detect collinear or non-collinear spin(s). Check spin input." << endl;//DX 12/5/17 - added non-collinear
          exit(0);                                                                                                   //DX 12/5/17 - added non-collinear
        }
        //DX [OBSOLETE] 12/5/17 - if(!AddSpinToXstructure(a,vmag)){
        //DX [OBSOLETE] 12/5/17 -   exit(0);
        //DX [OBSOLETE] 12/5/17 - } 
      }
      //DX 9/21/17 - MAGNETIC SYMMETRY - END
      //DX [OBSOLETE] 9/21/17 - double default_tolerance=SYM::defaultTolerance(a);
      //DX [OBSOLETE] 9/21/17 - double tolerance = default_tolerance;
      //DX [OBSOLETE] 9/21/17 - if(tokens.size()==0) {tolerance=default_tolerance;}
      //DX [OBSOLETE] 9/21/17 -if(tokens.size()>=1 && tokens.at(0) != "--debug") {
      //DX [OBSOLETE] 9/21/17 -  if(tokens.at(0).at(0) == 't' || tokens.at(0).at(0) == 'T'){ //Tight
      //DX [OBSOLETE] 9/21/17 -    tolerance=default_tolerance;
      //DX [OBSOLETE] 9/21/17 -  }
      //DX [OBSOLETE] 9/21/17 -  else if(tokens.at(0).at(0) == 'l' || tokens.at(0).at(0) == 'L'){ //Loose
      //DX [OBSOLETE] 9/21/17 -    tolerance=default_tolerance*10.0;
      //DX [OBSOLETE] 9/21/17 -  }
      //DX [OBSOLETE] 9/21/17 -  else{
      //DX [OBSOLETE] 9/21/17 -    tolerance=aurostd::string2utype<double>(tokens.at(0));
      //DX [OBSOLETE] 9/21/17 -  }
      //DX [OBSOLETE] 9/21/17 -}
      // DX - END
      double default_tolerance=SYM::defaultTolerance(a);
      double tolerance = AUROSTD_NAN;
      if(vpflow.flag("SG::TOLERANCE")){
	string tolerance_string = vpflow.getattachedscheme("SG::TOLERANCE");
	if(tolerance_string[0] == 't' || tolerance_string[0] == 'T'){ //Tight
          tolerance=default_tolerance;
        }
	else if(tolerance_string[0] == 'l' || tolerance_string[0] == 'L'){ //Loose
          tolerance=default_tolerance*10.0;
        }
        else{
	  tolerance=aurostd::string2utype<double>(vpflow.getattachedscheme("SG::TOLERANCE"));
        }
      }
      else{
	tolerance = default_tolerance;
      }
      if(tolerance < 1e-10){
	cerr << "pflow::SG::ERROR: Tolerance cannot be zero (i.e. less than 1e-10)." << endl;
	return 0;
      }
      //DX 9/26/17 - NO SCAN - START
      bool no_scan = false;
      if(vpflow.flag("SG::NO_SCAN")){
        no_scan = true;
      }
      //DX 9/26/17 - NO SCAN - END
      uint sgroup=a.SpaceGroup_ITC(tolerance,no_scan);
      // [OBSOLETE] uint sgroup=a.SpaceGroup_ITC(false,argv);    
      a.spacegroup=GetSpaceGroupName(sgroup)+" #"+aurostd::utype2string(sgroup);
      // [OBSOLETE] a.spacegroup=GetSpaceGroupName(a.SpaceGroup_ITC(false,argv))+" #"+aurostd::utype2string(a.SpaceGroup_ITC(false,argv)); // RHT
      //  return a.spacegroup; // RHT
    }
    
    if(mode=="PLATON" || mode=="platon") {
      if(LDEBUG) cerr << "pflow::SG: platon" << endl;
      // SIMPLE CALCULATION
      // output << a.platon2sg() << endl;
      // PERFECT CALCULATION
      bool Platon_EQUAL=DEFAULT_PLATON_P_EQUAL;
      bool Platon_EXACT=DEFAULT_PLATON_P_EXACT;
      double Platon_ang=DEFAULT_PLATON_P_ANG;
      double Platon_d1=DEFAULT_PLATON_P_D1;
      double Platon_d2=DEFAULT_PLATON_P_D2;
      double Platon_d3=DEFAULT_PLATON_P_D3;
      //DX 9/26/17 - FLAGS - START
      if(vpflow.flag("SG::TOLERANCE")){
	string tolerance_string = vpflow.getattachedscheme("SG::TOLERANCE");
	vector<string> tol_tokens;
	aurostd::string2tokens(tolerance_string,tol_tokens,",");
      //     if(tokens.size()==0) Platon_ang=1.0e-2;
      // if(tokens.size()==1) Platon_ang=aurostd::string2utype<double>(tokens.at(0));
      // Read in input file.
	if((tol_tokens.size()>=1 && tol_tokens.at(0)=="EQUAL") || (tol_tokens.size()>=2 && tol_tokens.at(1)=="EQUAL")) Platon_EQUAL=TRUE;
	if((tol_tokens.size()>=1 && tol_tokens.at(0)=="EXACT") || (tol_tokens.size()>=2 && tol_tokens.at(1)=="EXACT")) Platon_EXACT=TRUE;
	if(tol_tokens.size()>=3) {
	  Platon_ang=aurostd::string2utype<double>(tol_tokens.at(tol_tokens.size()-4));
	  Platon_d1=aurostd::string2utype<double>(tol_tokens.at(tol_tokens.size()-3));
	  Platon_d2=aurostd::string2utype<double>(tol_tokens.at(tol_tokens.size()-2));
	  Platon_d3=aurostd::string2utype<double>(tol_tokens.at(tol_tokens.size()-1));
      }
      }
      //DX 9/26/17 - FLAGS - END
      a.platon2sg(Platon_EQUAL,Platon_EXACT,Platon_ang,Platon_d1,Platon_d2,Platon_d3);
      //  return a.spacegroup;
    }
    
    if(mode=="FINDSYM" || mode=="findsym") {
      if(LDEBUG) cerr << "pflow::SG: findsym" << endl;
      // SIMPLE CALCULATION
      string out;
      double tolerance=DEFAULT_FINDSYM_TOL;
      if(vpflow.flag("SG::TOLERANCE")){
        string tolerance_string = vpflow.getattachedscheme("SG::TOLERANCE");
        vector<string> tol_tokens;
        aurostd::string2tokens(tolerance_string,tol_tokens,",");
        if(tol_tokens.size()==0) tolerance=DEFAULT_FINDSYM_TOL;
        if(tol_tokens.size()==1) tolerance=aurostd::string2utype<double>(tol_tokens.at(0));
      }
      a.findsym2sg(tolerance);
      //   return a.spacegroup;
    }
    // what to print 
    
    //   vector<string> tokens;
    if(print=="LABEL" || print=="label") {aurostd::string2tokens(a.spacegroup,tokens," ");return tokens[0];}
    if(print=="NUMBER" || print=="number") {aurostd::string2tokens(a.spacegroup,tokens," ");return tokens[1];}
    return a.spacegroup;
    
    //  return NOSG;
  }
} // namespace pflow

// ***************************************************************************
// pflow::SHELL
// ***************************************************************************
namespace pflow {
  void SHELL(string options,istream& input) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "pflow::SHELL: BEGIN" << endl;
    vector<string> tokens;
    aurostd::string2tokens(options,tokens,",");
    if(tokens.size()!=5) {
      init::ErrorOption(cout,options,"pflow::SHELL","aflow --shell=ns,r1,r2,name,dens < POSCAR");
      exit(0);
    } 
    
    int ns=4;  // some defaults
    double r1=1.8,r2=2.2; // some defaults
    string name=""; // some defaults
    int dens=20; // some defaults
    if(tokens.size()>=1) ns=aurostd::string2utype<int>(tokens.at(0));
    if(tokens.size()>=2) r1=aurostd::string2utype<double>(tokens.at(1));
    if(tokens.size()>=3) r2=aurostd::string2utype<double>(tokens.at(2));
    if(tokens.size()>=4) name=tokens.at(3);
    if(tokens.size()>=5) ns=aurostd::string2utype<int>(tokens.at(4));
    
    cout << aflow::Banner("BANNER_TINY") << endl;
    xstructure a(input,IOAFLOW_AUTO);
    if(0) if(tokens.size()!=a.num_each_type.size()) {
	cerr << "ERROR - pflow::SHELL: you need to specify as many names as atom types" << endl;
	exit(0);
      }

    PrintShell(a,ns,r1,r2,name,dens,cout);
  
    if(LDEBUG) cerr << "pflow::SHELL: END" << endl;
  }
} // namespace pflow

// ***************************************************************************
// pflow::SHIFT
// ***************************************************************************
namespace pflow {
  xstructure SHIFT(string options,istream& input) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "pflow::SHIFT: BEGIN" << endl;
    vector<string> tokens;
    aurostd::string2tokens(options,tokens,",");
    if(tokens.size()!=3 && tokens.size()!=4) {
      init::ErrorOption(cout,options,"pflow::SHIFT","aflow --shift=Sx,Sy,Sz[,cCdD] < POSCAR");
      exit(0);
    }      
    // cout << aflow::Banner("BANNER_TINY") << endl;
    // Get shift vector
    xvector<double> shift(3);    
    for(int ic=1;ic<=3;ic++) {shift[ic]=aurostd::string2utype<double>(tokens.at(ic-1));}
    // argv(4)=c cart,d for direct shift
    bool flag=FALSE;
    if(tokens.size()==3) {
      cerr << "WARNING - pflow::SHIFT: 4th argument can be [cCfFdD], defaulting to [cC]=Cartesian shift" << endl;
    } else {
      if(tokens.at(3)=="c" || tokens.at(3)=="C") {
	flag=FALSE;
      } else {
	if(tokens.at(3)=="d" || tokens.at(3)=="D" || tokens.at(3)=="f" || tokens.at(3)=="F") {
	  flag=TRUE;
	} else{
	  cerr << "WARNING - pflow::SHIFT: 4th argument can be [cCfFdD], Defaulting to [cC]=Cartesian shift" << endl;
	  flag=FALSE;
	} //else
      } //else
    }//else
    // Read in input file.
    xstructure str(input,IOAFLOW_AUTO);
    str=ShiftPos(str,shift,flag);
    if(LDEBUG) cerr << "pflow::SHIFT: END" << endl;
    return str;
  }
} // namespace pflow

// ***************************************************************************
// pflow::SG
// ***************************************************************************
namespace pflow {
  void SG(istream& input) {
    xstructure a(input,IOAFLOW_AUTO);
    a.CalculateSymmetryFactorGroup(TRUE);
    spacegroup::SpaceGroupInitialize();
    cerr << a.pgroup.size() << "," << a.fgroup.size() << endl;
    cerr.flush();
    xmatrix<double> I(3,3),U(3,3),A(3,3);
    xvector<double> tau(3);
    I=identity(I);
    for(uint i=0;i<a.fgroup.size();i++) {
      U=a.fgroup.at(i).Uf;
      tau=a.fgroup.at(i).ftau;
      A=I-U;
      if(abs(det(A))>0.01) {
	cerr << inverse(A)*tau << endl;
      }
    }
    //  spacegroup::SpaceGroupNumberStructure(a);
  }
} // namespace pflow

// ***************************************************************************
// pflow::SGDATA
// ***************************************************************************
namespace pflow {
  bool SGDATA(istream& input, aurostd::xoption& vpflow, ostream& oss) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string options = vpflow.getattachedscheme("SGDATA");
    vector<string> tokens;
    aurostd::string2tokens(options,tokens,",");
    if(tokens.size()==1) {
      if(tokens.at(0)=="usage" || tokens.at(0)=="USAGE") {
        return init::ErrorOption(cout,options,"pflow::SGDATA",
                          aurostd::liststring2string("aflow --sgdata|--space_group_data[=tolerance| =tight| =loose] [--no_scan] [--print=txt| =json] [--mag|--magnetic|--magmom=[m1,m2,...|INCAR|OUTCAR]] < POSCAR  default: tolerance=(minimum_interatomic_distance)/100.0, print=txt, non-magnetic"));
      }
    }
    xstructure _a(input,IOAFLOW_AUTO);
    xstructure a(_a);
    a.ReScale(1.0);
    //DX 2/21/18 - use pwd - START
    if(a.directory == ""){
      a.directory = aurostd::execute2string("pwd");
    }
    string print_directory = " [dir=" + a.directory + "]";
    //DX 2/21/18 - use pwd - END
    //DX 9/21/17 - MAGNETIC SYMMETRY - START
    if(vpflow.flag("SGDATA::MAGNETIC")){
      string magmom_info = vpflow.getattachedscheme("SGDATA::MAGNETIC");
      int num_atoms=a.atoms.size();
      bool is_noncoll=false; 
      vector<xvector<double> > vmag_noncoll;                                          //DX 12/5/17 - added non-collinear
      if(GetNonCollinearMagneticInfo(num_atoms,magmom_info,vmag_noncoll)){            //DX 12/5/17 - added non-collinear
        if(LDEBUG){cerr << "pflow::SG: Non-collinear spin system detected." << endl;} //DX 12/5/17 - added non-collinear
        is_noncoll = true;
        if(!AddSpinToXstructure(a,vmag_noncoll)){
          return false;
        } 
      }                                                                               //DX 12/5/17 - added non-collinear
      bool is_coll=false; 
      vector<double> vmag;
      if(!is_noncoll){
        if(GetCollinearMagneticInfo(num_atoms,magmom_info,vmag)){
          if(LDEBUG){cerr << "pflow::SG: Collinear spin system detected." << endl;}   //DX 12/5/17 - added non-collinear
          is_coll = true;
          if(!AddSpinToXstructure(a,vmag)){
            return false;
          } 
        } 
      }
      if(!is_noncoll && !is_coll){                                                                                 //DX 12/5/17 - added non-collinear
        cerr << "pflow::SG: ERROR: Could not detect collinear or non-collinear spin(s). Check spin input." << endl;//DX 12/5/17 - added non-collinear
        return false;                                                                                              //DX 12/5/17 - added non-collinear
      } 
    } 
    //DX 9/21/17 - MAGNETIC SYMMETRY - END

    double default_tolerance=SYM::defaultTolerance(a);
    double tolerance = AUROSTD_NAN;
    if(vpflow.flag("SGDATA::TOLERANCE")){
      string tolerance_string = vpflow.getattachedscheme("SGDATA::TOLERANCE");
      if(tolerance_string[0] == 't' || tolerance_string[0] == 'T'){ //Tight
        tolerance=default_tolerance;
      }
      else if(tolerance_string[0] == 'l' || tolerance_string[0] == 'L'){ //Loose
        tolerance=default_tolerance*10.0;
      }
      else{
        tolerance=aurostd::string2utype<double>(vpflow.getattachedscheme("SGDATA::TOLERANCE"));
      }
    }
    else{
      tolerance = default_tolerance;
    }
    if(tolerance < 1e-10){
      cerr << "pflow::SGDATA::ERROR: Tolerance cannot be zero (i.e. less than 1e-10) " << print_directory << "." << endl;
      return 0;
    }
    // DX 8/3/17 - Add format flag - START
    string format = "txt";
    if(XHOST.vflag_control.flag("PRINT_MODE::TXT")){
        format = "txt";
      }
    else if(XHOST.vflag_control.flag("PRINT_MODE::JSON")){
        format = "json";
    }
    else{ //default is txt
      format = "txt";
    }   
    bool no_scan = false;
    if(vpflow.flag("SGDATA::NO_SCAN")){
      no_scan = true;
    }
    if(format=="txt"){
      cout << aflow::Banner("BANNER_TINY") << endl;
    }
    // pflow::PrintData(a,cerr,smode);
    pflow::PrintSGData(a,tolerance,oss,no_scan,true,format,false); //CO 171027
    return true;
  }
} // namespace pflow

// ***************************************************************************
// pflow::SGROUP
// ***************************************************************************
namespace pflow {
  void SGROUP(_aflags &aflags,istream& input,double radius) {
    cout << aflow::Banner("BANNER_TINY") << endl;
    aflags.QUIET=TRUE;
    xstructure a(input,IOAFLOW_AUTO);
    bool WRITE=TRUE;
    ofstream File("/dev/null");
    SYM::CalculatePointGroup(File,a,aflags,WRITE,TRUE,cout);
    SYM::CalculatePointGroupKlattice(File,a,aflags,WRITE,TRUE,cout);
    SYM::CalculateFactorGroup(File,a,aflags,WRITE,TRUE,cout);
    SYM::CalculatePointGroupCrystal(File,a,aflags,WRITE,TRUE,cout);
    a.sgroup_radius=radius;
    SYM::CalculateSpaceGroup(File,a,aflags,WRITE,TRUE,cout);
  }
} // namespace pflow

// ***************************************************************************
// pflow::SPECIES
// ***************************************************************************
namespace pflow {
  string SPECIES(istream& input) {
    xstructure a(input,IOAFLOW_AUTO);
    string strout=a.SpeciesString()+"\n";
    return strout;
  }
} // namespace pflow

// ***************************************************************************
// pflow::SPLINE
// ***************************************************************************
namespace pflow {
  void SPLINE(vector<string> argv) {
    // Read in input data.
    vector<double> x,y;
    double p,q;
    char dum[500];
    while(cin >> p >> q) {
      x.push_back(p);
      y.push_back(q);
      cin.getline(dum,500); // Get anything else on line.
    }
    int npts=atoi(argv.at(2).c_str());
    pflow::PrintSpline(x,y,npts,cout);
  }
} // namespace pflow

// ***************************************************************************
// pflow::SUMPDOS
// ***************************************************************************
namespace pflow {
  void SUMPDOS(vector<string> argv) {
    cerr << "# WARNING: THIS REQUIRES THAT YOU HAVE PDOS IN DOSCAR FILE - THIS IS OBTAINED BY RUNNING WITH LORBIT=2 - SEE aflow --h" << endl;
    ifstream SumPDOSParams_infile(argv.at(2).c_str());
    aurostd::InFileExistCheck("convasp.cc",argv.at(2),SumPDOSParams_infile,cerr);
    ifstream PDOS_infile(argv.at(3).c_str());
    aurostd::InFileExistCheck("convasp.cc",argv.at(3),PDOS_infile,cerr);
    pflow::pdosdata pdd;
    pflow::matrix<pflow::matrix<double> > allpdos;
    pflow::ReadSumDOSParams(SumPDOSParams_infile,pdd);
    //      projdata prd;
    //      pdd.PrintParams(cout,prd.LMnames);
    pflow::ReadInPDOSData(allpdos,pdd,PDOS_infile);
    SumPDOS(allpdos,pdd);
  }
} // namespace pflow

// ***************************************************************************
// pflow::SUPERCELL
// ***************************************************************************
namespace pflow {
  xstructure SUPERCELL(string options,istream& input) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "pflow::SUPERCELL: BEGIN" << endl;
    xstructure str(input,IOAFLOW_AUTO);
    vector<string> tokens;
    aurostd::string2tokens(options,tokens,",");
    xmatrix<double> msc(3,3);

    if(tokens.size()!=9 && tokens.size()!=3 && tokens.size()!=1) {
      init::ErrorOption(cout,options,"pflow::SUPERCELL",
			aurostd::liststring2string("aflow --supercell=a11,a12,a13,a21,a22,a23,a31,a32,a33 < POSCAR",
                                                   "aflow --supercell=a11,a22,a33 < POSCAR",
                                                   "aflow --supercell=file < POSCAR"));
      exit(0);
    }

    if(tokens.size()==9) {
      if(LDEBUG) cerr << "pflow::SUPERCELL: 9 entries" << endl;
      msc(1,1)=aurostd::string2utype<double>(tokens.at(0)); 
      msc(1,2)=aurostd::string2utype<double>(tokens.at(1)); 
      msc(1,3)=aurostd::string2utype<double>(tokens.at(2)); 
      msc(2,1)=aurostd::string2utype<double>(tokens.at(3)); 
      msc(2,2)=aurostd::string2utype<double>(tokens.at(4)); 
      msc(2,3)=aurostd::string2utype<double>(tokens.at(5)); 
      msc(3,1)=aurostd::string2utype<double>(tokens.at(6)); 
      msc(3,2)=aurostd::string2utype<double>(tokens.at(7)); 
      msc(3,3)=aurostd::string2utype<double>(tokens.at(8)); 
      if(abs(det(msc))<0.01) {cerr << "ERROR - pflow::SUPERCELL: singular supercell matrix" << endl;exit(0);}
      if(LDEBUG) cerr << "pflow::SUPERCELL: END" << endl;
      return GetSuperCell(str,msc);
    }
    
    if(tokens.size()==3) {
      if(LDEBUG) cerr << "pflow::SUPERCELL: 3 entries" << endl;
      msc(1,1)=aurostd::string2utype<double>(tokens.at(0)); 
      msc(2,2)=aurostd::string2utype<double>(tokens.at(1)); 
      msc(3,3)=aurostd::string2utype<double>(tokens.at(2)); 
      if(abs(det(msc))<0.01) {cerr << "ERROR - pflow::SUPERCELL: singular supercell matrix" << endl;exit(0);}
      if(LDEBUG) cerr << "pflow::SUPERCELL: END" << endl;
      return GetSuperCell(str,msc);
    }
    
    if(tokens.size()==1) {
      if(LDEBUG) cerr << "pflow::SUPERCELL: 1 entries" << endl;
      ifstream infile(tokens.at(0).c_str());
      aurostd::InFileExistCheck("pflow::SUPERCELL",tokens.at(0).c_str(),infile,cerr);
      for(int i=1;i<=3;i++)
        for(int j=1;j<=3;j++)
          infile >> msc(i,j);
      if(abs(det(msc))<0.01) {cerr << "ERROR - pflow::SUPERCELL: singular supercell matrix" << endl;exit(0);}
      if(LDEBUG) cerr << "pflow::SUPERCELL: END" << endl;
      return GetSuperCell(str,msc);
    }
    return str;
  }
} // namespace pflow

// ***************************************************************************
// pflow::SUPERCELLSTRLIST
// ***************************************************************************
namespace pflow {
  void SUPERCELLSTRLIST(string options) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "pflow::SUPERCELLSTRLIST: BEGIN" << endl;
    vector<string> tokens;
    aurostd::string2tokens(options,tokens,",");
    xmatrix<double> msc(3,3);

    string infile_name="";

    if(tokens.size()!=10 && tokens.size()!=4 && tokens.size()!=2) {
      init::ErrorOption(cout,options,"pflow::SUPERCELLSTRLIST",
			aurostd::liststring2string("aflow --supercell_strlist=a11,a12,a13,a21,a22,a23,a31,a32,a33,strlist",
                                                   "aflow --supercell_strlist=a11,a22,a33,strlist",
                                                   "aflow --supercell_strlist=file,strlist"));
      exit(0);
    }

    if(tokens.size()==10) {
      if(LDEBUG) cerr << "pflow::SUPERCELLSTRLIST: 10 entries" << endl;
      msc(1,1)=aurostd::string2utype<double>(tokens.at(0)); 
      msc(1,2)=aurostd::string2utype<double>(tokens.at(1)); 
      msc(1,3)=aurostd::string2utype<double>(tokens.at(2)); 
      msc(2,1)=aurostd::string2utype<double>(tokens.at(3)); 
      msc(2,2)=aurostd::string2utype<double>(tokens.at(4)); 
      msc(2,3)=aurostd::string2utype<double>(tokens.at(5)); 
      msc(3,1)=aurostd::string2utype<double>(tokens.at(6)); 
      msc(3,2)=aurostd::string2utype<double>(tokens.at(7)); 
      msc(3,3)=aurostd::string2utype<double>(tokens.at(8)); 
      infile_name=tokens.at(9);
      if(abs(det(msc))<0.01) {cerr << "ERROR - pflow::SUPERCELLSTRLIST: singular supercell matrix" << endl;exit(0);}
    }
    if(tokens.size()==4) {
      if(LDEBUG) cerr << "pflow::SUPERCELLSTRLIST: 4 entries" << endl;
      msc(1,1)=aurostd::string2utype<double>(tokens.at(0)); 
      msc(2,2)=aurostd::string2utype<double>(tokens.at(1)); 
      msc(3,3)=aurostd::string2utype<double>(tokens.at(2)); 
      infile_name=tokens.at(3);
      if(abs(det(msc))<0.01) {cerr << "ERROR - pflow::SUPERCELLSTRLIST: singular supercell matrix" << endl;exit(0);}
    }    
    if(tokens.size()==2) {
      if(LDEBUG) cerr << "pflow::SUPERCELLSTRLIST: 2 entries" << endl;
      ifstream infile(tokens.at(0).c_str());
      aurostd::InFileExistCheck("pflow::SUPERCELLSTRLIST",tokens.at(0).c_str(),infile,cerr);
      for(int i=1;i<=3;i++)
        for(int j=1;j<=3;j++)
          infile >> msc(i,j);
      infile_name=tokens.at(1);
      if(abs(det(msc))<0.01) {cerr << "ERROR - pflow::SUPERCELLSTRLIST: singular supercell matrix" << endl;exit(0);}
    }
    //    ifstream infile(infile_name.c_str());
    // aurostd::InFileExistCheck("aflow",infile_name,infile,cerr);
    pflow::matrix<double> mmsc(3,3);
    // if(LDEBUG) cerr << "pflow::SUPERCELLSTRLIST: msc=" << msc << endl;
    mmsc=pflow::xmatrix2matrix(msc);
    //  pflow::Mout(mmsc,cout);
    ifstream list_inf(infile_name.c_str());
    aurostd::InFileExistCheck("pflow::SUPERCELLSTRLIST",infile_name,list_inf,cerr);
    vector<xstructure> vstr;
    vector<xstructure> vstr_sc;
    pflow::ReadInStrVec(vstr,list_inf);
    for(uint i=0;i<vstr.size();i++) vstr_sc.push_back(vstr.at(i));
    pflow::SuperCellStrVec(vstr_sc,mmsc);
    pflow::PrintStrVec(vstr_sc,cout);
    cerr << "pflow::SUPERCELLSTRLIST: vstr_sc.size()=" << vstr_sc.size() << endl;
  } 
} // namespace pflow

// ***************************************************************************
// pflow::xstrSWAP
// ***************************************************************************
namespace pflow {
  xstructure xstrSWAP(vector<string> argv,istream& input) {
    xstructure a(input,IOAFLOW_AUTO);
    int speciesA=aurostd::string2utype<int>(argv.at(2));
    int speciesB=aurostd::string2utype<int>(argv.at(3));
    if(speciesA<0) {cerr << "pflow::xstrSWAP: Error, speciesA<0 (speciesA=" << speciesA << ")" << endl;exit(0);}
    if(speciesA>=(int) a.num_each_type.size()) {cerr << "pflow::xstrSWAP: Error, speciesA>=num_each_type.size() (speciesA=" << speciesA << ")" << endl;exit(0);}
    if(speciesB<0) {cerr << "pflow::xstrSWAP: Error, speciesB<0 (speciesB=" << speciesB << ")" << endl;exit(0);}
    if(speciesB>=(int) a.num_each_type.size()) {cerr << "pflow::xstrSWAP: Error, speciesB>=num_each_type.size() (speciesB=" << speciesB << ")" << endl;exit(0);}
    a.SpeciesSwap(speciesA,speciesB);
    return a;
  }
} // namespace pflow

// ***************************************************************************
// pflow::VOLUME
// ***************************************************************************
namespace pflow {
  xstructure VOLUME(string options, istream& input) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "pflow::VOLUME: BEGIN" << endl;  
    vector<string> tokens;
    aurostd::string2tokens(options,tokens,",");
    if(tokens.size()!=2) {
      init::ErrorOption(cout,options,"pflow::VOLUME","aflow --volume[|*|+]=x < POSCAR");
      exit(0);
    }
    xstructure a(input,IOAFLOW_AUTO);
    if(tokens.at(0)=="VOLUME::EQUAL") a=SetVolume(a,aurostd::string2utype<double>(tokens.at(1)));
    if(tokens.at(0)=="VOLUME::MULTIPLY_EQUAL") a=SetVolume(a,a.Volume()*aurostd::string2utype<double>(tokens.at(1)));
    if(tokens.at(0)=="VOLUME::PLUS_EQUAL") a=SetVolume(a,a.Volume()+aurostd::string2utype<double>(tokens.at(1)));
    if(LDEBUG) cerr << "pflow::VOLUME: END" << endl;  
    return a;
  }
} // namespace pflow

// ***************************************************************************
// pflow::WYCKOFF
// ***************************************************************************
namespace pflow {
  xstructure WYCKOFF(vector<string> argv,istream& input) {
    xstructure a(input,IOAFLOW_AUTO);
    int sg=aurostd::string2utype<int>(argv.at(2));
    a=WyckoffPOSITIONS(sg,a);
    cerr << a.spacegroup << endl;
    cerr << a.spacegroupnumber << endl;
    cerr << a.spacegroupoption << endl;
    return a;
  }
} // namespace pflow

// ***************************************************************************
// pflow::XRAY
// ***************************************************************************
namespace pflow {
  void XRAY(string options,istream& input) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "pflow::XRAY: BEGIN" << endl;  
    vector<string> tokens;
    aurostd::string2tokens(options,tokens,",");
    if(tokens.size()!=1) {
      init::ErrorOption(cout,options,"pflow::XRAY","aflow --xray=lambda < POSCAR");
      exit(0);
    }
    double l=0.0;
    if(tokens.size()>=1) l=aurostd::string2utype<double>(tokens.at(0)); 
 
    xstructure a(input,IOAFLOW_AUTO);
    cout << aflow::Banner("BANNER_TINY") << endl;
    PrintXray(a,l,cout);
    if(LDEBUG) cerr << "pflow::XRAY: END" << endl;  
  }
} // namespace pflow


// ***************************************************************************
// pflow::XYZ
// ***************************************************************************
namespace pflow {
  void XYZ(string options,istream& input) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "pflow::XYZ: BEGIN" << endl;
    xvector<int> ijk(3);           // default 
    ijk[1]=1;ijk[2]=1;ijk[3]=1;    // default
    vector<string> tokens;
    aurostd::string2tokens(options,tokens,",");

    if(tokens.size()>3) {
      init::ErrorOption(cout,options,"pflow::XYZ","aflow --xyz[=n1[,n2[,n3]]] < POSCAR");
      exit(0);
    }
    
    if(tokens.size()>=1) ijk[1]=aurostd::string2utype<int>(tokens.at(0)); 
    if(tokens.size()>=2) ijk[2]=aurostd::string2utype<int>(tokens.at(1)); 
    if(tokens.size()>=3) ijk[3]=aurostd::string2utype<int>(tokens.at(2)); 
    
    xstructure a(input,IOAFLOW_AUTO);
    PrintXYZ(a,ijk,cout);
    if(LDEBUG) cerr << "pflow::XYZ: END" << endl;
  }
} // namespace pflow

// ***************************************************************************
// pflow::XYZINSPHERE
// ***************************************************************************
namespace pflow {
  void XYZINSPHERE(istream& input,double radius) {
    cout << aflow::Banner("BANNER_TINY") << endl;
    xstructure a(input,IOAFLOW_AUTO);
    PrintXYZInSphere(a,radius,cout);
  }
} // namespace pflow

// ***************************************************************************
// pflow::XYZWS
// ***************************************************************************
namespace pflow {
  void XYZWS(istream& input) {
    cout << aflow::Banner("BANNER_TINY") << endl;
    xstructure a(input,IOAFLOW_AUTO);
    PrintXYZws(a,cout);
  }
} // namespace pflow

// ***************************************************************************
// pflow::ZVAL
// ***************************************************************************
namespace pflow {
  void ZVAL(string options) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "pflow::ZVAL: BEGIN" << endl;  
    vector<string> tokens;
    aurostd::string2tokens(options,tokens,",");
    if(tokens.size()==0 || tokens.size()>2) {
      init::ErrorOption(cout,options,"pflow::ZVAL","aflow "+tokens.at(0)+"[=directory]");
      exit(0);
    }
    string Directory;
    if(tokens.size()==1) Directory="./";
    if(tokens.size()==2) Directory=tokens.at(1);
   
    vector<double> vZVAL,sZVAL,vPOMASS,sPOMASS;
    if(tokens.at(0)=="ZVAL")        cout << "Total ZVAL (from PP) = " << GetZVAL(Directory,vZVAL) << endl;   
    if(tokens.at(0)=="ZVAL::CELL")   cout << "Total ZVAL_CELL (from PP and xstructure) = " << GetCellAtomZVAL(Directory,vZVAL,sZVAL,"CELL") << endl;
    if(tokens.at(0)=="ZVAL::ATOM")   cout << "Total ZVAL_ATOM (from PP and xstructure) = " << GetCellAtomZVAL(Directory,vZVAL,sZVAL,"ATOM") << endl;
    if(tokens.at(0)=="POMASS")      cout << "Total POMASS (from PP) = " << GetPOMASS(Directory,vPOMASS) << endl;   
    if(tokens.at(0)=="POMASS::CELL") cout << "Total POMASS_CELL (from PP and xstructure) = " << GetCellAtomPOMASS(Directory,vPOMASS,sPOMASS,"CELL") << endl;
    if(tokens.at(0)=="POMASS::ATOM") cout << "Total POMASS_ATOM (from PP and xstructure) = " << GetCellAtomPOMASS(Directory,vPOMASS,sPOMASS,"ATOM") << endl;
    if(LDEBUG) cerr << "pflow::ZVAL: END" << endl;  
  }
} // namespace pflow


// ***************************************************************************
// pflow::OPARAMETER
// ***************************************************************************

bool AreAtomsEquivalent(const deque<_atom> &atoms1,const deque<_atom> &atoms2,bool oparameter_check,bool poccupation_check) {
  double epsilon=0.0001,epsilon_pocc=0.05;
  if(atoms1.size()!=atoms2.size()) return FALSE; // cant be equivalent if different numbers
  for(uint iat1=0;iat1<atoms1.size();iat1++) {
    bool found=FALSE;
    for(uint iat2=0;iat2<atoms2.size() && found==FALSE;iat2++) {
      if(atoms1.at(iat1).type==atoms2.at(iat2).type &&
	 aurostd::isequal(atoms1.at(iat1).fpos,atoms2.at(iat2).fpos,epsilon)) {
	if(oparameter_check==FALSE && poccupation_check==FALSE) found=TRUE;
	if(oparameter_check==TRUE && poccupation_check==FALSE)
	  if(atoms1.at(iat1).order_parameter_value==atoms2.at(iat2).order_parameter_value) found=TRUE;
	//  atoms1.at(iat1).order_parameter_atom==TRUE && atoms2.at(iat2).order_parameter_atom==TRUE &&
	if(oparameter_check==FALSE && poccupation_check==TRUE)
	  if(aurostd::isequal(atoms1.at(iat1).partial_occupation_value,atoms2.at(iat2).partial_occupation_value,epsilon_pocc)) found=TRUE;
	//  atoms1.at(iat1).partial_occupation_atom==TRUE && atoms2.at(iat2).partial_occupation_atom==TRUE &&
	if(oparameter_check==TRUE && poccupation_check==TRUE)
	  if(atoms1.at(iat1).order_parameter_value==atoms2.at(iat2).order_parameter_value &&
	     aurostd::isequal(atoms1.at(iat1).partial_occupation_value,atoms2.at(iat2).partial_occupation_value,epsilon_pocc)) found=TRUE;
	//  atoms1.at(iat1).order_parameter_atom==TRUE && atoms2.at(iat2).order_parameter_atom==TRUE &&
	//  atoms1.at(iat1).partial_occupation_atom==TRUE && atoms2.at(iat2).partial_occupation_atom==TRUE &&
      }
    }
    if(found==FALSE) return FALSE;
  }
  return TRUE; // survived everything
}

bool AreAtomsEquivalent(const deque<_atom> &atoms1,const deque<_atom> &atoms2) {
  return AreAtomsEquivalent(atoms1,atoms2,FALSE,FALSE); // no order check
}

bool AreAtomsOrderEquivalent(const deque<_atom> &atoms1,const deque<_atom> &atoms2) {
  return AreAtomsEquivalent(atoms1,atoms2,TRUE,FALSE); // no order check
}

bool AreAtomsOccupationEquivalent(const deque<_atom> &atoms1,const deque<_atom> &atoms2) {
  return AreAtomsEquivalent(atoms1,atoms2,FALSE,TRUE); // no order check
}

bool AreAtomsOccupationOrderEquivalent(const deque<_atom> &atoms1,const deque<_atom> &atoms2) {
  return AreAtomsEquivalent(atoms1,atoms2,TRUE,TRUE); // no order check
}

class _sort_xstructure_order_parameter_orbit {                    // sorting through reference
public:
  bool operator()(const xstructure& str1, const xstructure& str2) const {
    return (bool) (str1.order_parameter_orbit<str2.order_parameter_orbit);}
};

class _rsort_xstructure_order_parameter_orbit {                    // sorting through reference
public:
  bool operator()(const xstructure& str1, const xstructure& str2) const {
    return (bool) (str1.order_parameter_orbit>str2.order_parameter_orbit);}
};

int order_parameter_sum(const xstructure& str) {
  int _sum=0;
  for(uint iat=0;iat<str.order_parameter_atoms.size();iat++)
    _sum+=(int) str.atoms.at(str.order_parameter_atoms.at(iat)).order_parameter_value;
  return _sum;
}

/*
  #define VSTR_WRITE_OPAR      1
  #define VSTR_WRITE_POCC      2
  #define VSTR_WRITE_POCC_OPAR 3

  bool pflow::VSTR_Write(string directory,vector<xstructure> &vstr,bool FLAG_POSCAR,int mode) {
  uint index=0;
  for(uint istr=0;istr<vstr.size();istr++) {                       // cycle through the unique vstrs
  vstr.at(istr).title+=" [orbits="+int2string(vstr.at(istr).order_parameter_orbit)+"]";
  index+=vstr.at(istr).order_parameter_orbit;
  //  cout << istr << "  " << vstr.at(istr).order_parameter_orbit << " " << vstr.at(istr).order_parameter_sum << endl;
  //    if(vstr.at(istr).order_parameter_orbit==vstr.at(0).order_parameter_orbit)
  _xvasp xvasp;
  AVASP_DefaultValuesTernaryLSDAU_AFLOWIN(xvasp);
  xvasp.AVASP_dirbase=directory;
  if(mode==VSTR_WRITE_OPAR) xvasp.AVASP_label=string("RUN_os")+PaddedPRE(vstr.at(istr).order_parameter_sum,2,"0")+string("_od")+PaddedPRE(int2string(vstr.at(istr).order_parameter_orbit),3,"0")+"_t"+PaddedPRE(int2string(istr),4,"0");
  if(mode==VSTR_WRITE_POCC) xvasp.AVASP_label=string("RUN_po")+PaddedPRE(vstr.at(istr).partial_occupation_sum,2,"0")+string("_od")+PaddedPRE(int2string(vstr.at(istr).order_parameter_orbit),3,"0")+"_t"+PaddedPRE(int2string(istr),4,"0");
  if(mode==VSTR_WRITE_POCC_OPAR) xvasp.AVASP_label=string("RUN_os")+PaddedPRE(vstr.at(istr).order_parameter_sum,2,"0")+string("_po")+PaddedPRE(vstr.at(istr).partial_occupation_sum,2,"0")+string("_od")+PaddedPRE(int2string(vstr.at(istr).order_parameter_orbit),3,"0")+"_t"+PaddedPRE(int2string(istr),4,"0");
  xvasp.str=vstr.at(istr);
  xvasp.str.species.clear();xvasp.str.species_pp.clear();xvasp.str.species_pp_type.clear();xvasp.str.species_pp_version.clear();xvasp.str.species_pp_ZVAL.clear();xvasp.str.species_pp_vLDAU.clear();
  xvasp.str.species.push_back("Li_sv");xvasp.str.species.push_back("Mn");xvasp.str.species.push_back("O_s");
  xvasp.str.species_pp=xvasp.str.species;
  xvasp.str.species_pp_type.push_back("");xvasp.str.species_pp_type.push_back("");xvasp.str.species_pp_type.push_back("");
  xvasp.str.species_pp_version.push_back("");xvasp.str.species_pp_version.push_back("");xvasp.str.species_pp_version.push_back("");
  xvasp.str.species_pp_ZVAL.push_back(0.0);xvasp.str.species_pp_ZVAL.push_back(0.0);xvasp.str.species_pp_ZVAL.push_back(0.0);
  xvasp.str.species_pp_vLDAU.push_back(deque<double>());xvasp.str.species_pp_vLDAU.push_back(deque<double>());xvasp.str.species_pp_vLDAU.push_back(deque<double>());
  xvasp.AVASP_volume_in=558.7229;
  xvasp.AVASP_flag_AUTO_MAGMOM=TRUE;
  xvasp.AVASP_value_KPPRA=7000;
  xvasp.AVASP_flag_NBANDS=FALSE;
  xvasp.AVASP_flag_LSDAU=TRUE;
  xvasp.AVASP_Get_LDAU2_ParametersU_MODE=2;
  xvasp.AVASP_flag_PRECISION_flag=FALSE;
  xvasp.AVASP_flag_METAGGA_flag=FALSE;
  xvasp.AVASP_flag_IVDW_flag=FALSE;
  xvasp.AVASP_flag_TYPE_METAL=FALSE;
  xvasp.AVASP_flag_TYPE_INSULATOR=FALSE;
  xvasp.AVASP_flag_TYPE_SEMICONDUCTOR=FALSE;
  xvasp.AVASP_flag_EXTRA_INCAR=TRUE;
  xvasp.AVASP_value_NSW=101;
  stringstream _incar("");
  _incar << "#Electronic convergence" << endl;
  _incar << "PREC =  med" << endl;
  _incar << "ENMAX = 400" << endl;
  _incar << "#Initialization" << endl;
  _incar << "ISTART = 0" << endl;
  _incar << "ISMEAR = 1" << endl;
  _incar << "SIGMA = 0.2" << endl;
  _incar << "#Miscellaneous" << endl;
  _incar << "NPAR = 2" << endl;
  _incar << "#Output" << endl;
  _incar << "LORBIT = 10" << endl;
  xvasp.AVASP_EXTRA_INCAR << _incar.str();
  // cerr << xvasp.AVASP_label << endl;
  // cout << xvasp.str << endl;
  if(FLAG_POSCAR) AVASP_MakeSinglePOSCAR(xvasp);
  else AVASP_MakeSingleAFLOWIN(xvasp);  // obsolete stuff
  }
  return TRUE;
  }

  bool OPARAMETER(vector<string> argv,istream& input) {
  cout << "--------------------------------------------------------- "  << endl;
  cout << aflow::Banner("BANNER_TINY") << endl;
  // FLAGS
  bool FLAG_LIMNO=FALSE,FLAG_POSCAR=FALSE;
  if(argv.size()>=3 && aurostd::args2flag(argv,"--Li8Mn16O32|--LiMnO")) {
  FLAG_LIMNO=TRUE;
  }
  if(argv.size()>=3 && aurostd::args2flag(argv,"--poscar")) {
  FLAG_POSCAR=TRUE;
  }
  bool FLAG_OPARSUM=FALSE;int oparsum=0;
  if(argv.size()>=4 && aurostd::args2flag(argv,"--opar_sum|--order_parameter_sum|--oparameter_sum")) {
  oparsum=aurostd::args2utype(argv,"--opar_sum|--order_parameter_sum|--oparameter_sum",(int) 0);
  FLAG_OPARSUM=TRUE;
  cout << "order_parameter_sum=" << oparsum << endl;
  }
  // priorities
  if(FLAG_LIMNO) FLAG_OPARSUM=FALSE;
  // LOAD
  xstructure str(input,IOAFLOW_AUTO);
  str.BringInCell();
  str.order_parameter_sum=0;
  // GET SYMMETRY
  _aflags aflags;aflags.Directory="./";
  aflags.QUIET=TRUE;XHOST.QUIET=TRUE;
  ofstream FileMESSAGE("/dev/null");
  //  bool PGROUPWRITE=FALSE,FGROUPWRITE=FALSE,SGROUPWRITE=FALSE,IATOMSWRITE=TRUE,AGROUPWRITE=TRUE;
  bool PGROUPWRITE=FALSE,FGROUPWRITE=FALSE;//,SGROUPWRITE=FALSE,IATOMSWRITE=TRUE,AGROUPWRITE=TRUE;
  bool OSSWRITE=FALSE; // to FileMESSAGE, does not matter as it is /dev/null
  //  SGROUPWRITE=FALSE,IATOMSWRITE=TRUE,AGROUPWRITE=TRUE;
  SYM::CalculatePointGroup(FileMESSAGE,str,aflags,PGROUPWRITE,OSSWRITE,cout);
  cout << "pgroup operations  =" << str.pgroup.size() << endl;
  SYM::CalculateFactorGroup(FileMESSAGE,str,aflags,FGROUPWRITE,OSSWRITE,cout);
  cout << "fgroup operations  =" << str.fgroup.size() << endl;
  //  SYM::CalculatePointGroupKlattice(FileMESSAGE,str,aflags,PGROUPWRITE,OSSWRITE,cout);
  // cout << "pgroupk operations =" << str.pgroup.size() << endl;
  // SYM::CalculateSpaceGroup(FileMESSAGE,str,aflags,SGROUPWRITE,OSSWRITE,cout);
  //  cout << "sgroup operations=" << str.sgroup.size() << endl;
  // SYM::CalculateInequivalentAtoms(FileMESSAGE,str,aflags,IATOMSWRITE,OSSWRITE,cout);
  // SYM::CalculateSitePointGroup(FileMESSAGE,str,aflags,AGROUPWRITE,OSSWRITE,cout);
  //  cout << str.order_parameter_atoms.size() << endl;
  //for(uint i=0;i<str.order_parameter_atoms.size();i++)   cout << str.order_parameter_atoms.at(i) << endl;

  // fix the order parameter atoms
  for(uint iat=0;iat<str.order_parameter_atoms.size();iat++)
  str.atoms.at(str.order_parameter_atoms.at(iat)).order_parameter_atom=TRUE;

  // for(uint op=56;op<=64;op++) // hack to make them all
  {
  //    oparsum=op;  FLAG_OPARSUM=TRUE; // hack to make them all
  cout << "--------------------------------------------------------- "  << endl;
  cout << "order_parameter_sum=" << oparsum << endl;

  vector<xstructure> vstr;
  xstructure str_ord,str_aus;
  uint maxindex_Mn=(uint) pow(2,str.order_parameter_atoms.size());

  cout << "max_index_Mn=" << maxindex_Mn << endl;
  str_ord=str; // str_ord contains some order stuff
  bool found=FALSE;
  for(uint index_Mn=0;index_Mn<maxindex_Mn;index_Mn++) {
  // generate structure
  //   str_ord=str; // str_ord contains some order stuff
  // str_ord.title+=" ["+int2string(index_Mn)+"]";
  // str_ord.BringInCell(); // not really necessary
  str_ord.order_parameter_sum=0;
  for(uint iat=0,uindex_Mn=index_Mn;iat<str.order_parameter_atoms.size();iat++,uindex_Mn/=2) {
  //  str_ord.atoms.at(str.order_parameter_atoms.at(iat)).order_parameter_value=(int) 1-2*mod((int) uindex_Mn,2); // flip +3,+4 through the Mn !
  str_ord.atoms.at(str.order_parameter_atoms.at(iat)).order_parameter_value=(int) 3+mod((int) uindex_Mn,2); // flip +3,+4 through the Mn !
  str_ord.atoms.at(str.order_parameter_atoms.at(iat)).order_parameter_atom=TRUE;
  str_ord.order_parameter_sum+=str_ord.atoms.at(str.order_parameter_atoms.at(iat)).order_parameter_value;
  }
  //   str_ord.order_parameter_sum=order_parameter_sum(str_ord);
  // structure created

  if((FLAG_LIMNO && (str_ord.order_parameter_sum==56 || str_ord.order_parameter_sum==60)) || // 8Mn4+ 8Mn3+ fully-Li gs; 12Mn4+ 4Mn3+ half-Li state, metastabless
  (FLAG_OPARSUM && str_ord.order_parameter_sum==oparsum) ||                               // choice of oparsum
  (!FLAG_LIMNO && !FLAG_OPARSUM && str_ord.order_parameter_sum>=56 && str_ord.order_parameter_sum<=64)) {                                                       // do everything
  found=FALSE;
  for(uint ifgroup=0;ifgroup<str.fgroup.size()&&found==FALSE;ifgroup++) {                  // cycle through the fgroup of str
  str_aus=SYM_ApplyXstructure(str.fgroup.at(ifgroup),str_ord,TRUE);
  for(uint istr=0;istr<vstr.size()&&found==FALSE;istr++) {                               // cycle through the unique vstrs
  if(str_ord.order_parameter_sum==vstr.at(istr).order_parameter_sum)
  if(AreAtomsEquivalent(str_aus.atoms,vstr.at(istr).atoms,TRUE,FALSE)) {
  found=TRUE;
  vstr.at(istr).order_parameter_orbit++;
  }
  }
  }
  if(found==FALSE) { // not found, then plug
  if(1) { // chose this way, I`m afraid of loosing the memory
  xstructure *str_ptr;
  str_ord.order_parameter_orbit=1;
  str_ptr = new xstructure(str_ord);
  vstr.push_back(*str_ptr);
  }
  if(0) {
  str_ord.order_parameter_orbit=1;
  vstr.push_back(str_ord);
  }
  cout << vstr.size() << "(" << index_Mn << ") ";
  cout.flush();
  }
  }
  if(mod((int) index_Mn,1000)==0) {cout << "*";cout.flush();}// << str_ord.order_parameter_sum;
  }
  cout << endl << "unique orbits = " << vstr.size() << endl;
  sort(vstr.begin(),vstr.end(),_rsort_xstructure_order_parameter_orbit());
  pflow::VSTR_Write("./",vstr,FLAG_POSCAR,VSTR_WRITE_OPAR);
  }
  return TRUE;
  }
  // ******************** UNSAFE

  bool SHIRLEY(vector<string> argv,istream& input) {
  cout << "--------------------------------------------------------- "  << endl;
  cout << aflow::Banner("BANNER_TINY") << endl;
  cout << "argv.size()=" << argv.size() << endl;
  // FLAGS
  bool FLAG_LIMNO=FALSE,FLAG_POSCAR=FALSE;
  if(argv.size()>=3 && aurostd::args2flag(argv,"--Li8Mn16O32|--LiMnO")) {
  FLAG_LIMNO=TRUE;
  }
  if(argv.size()>=3 && aurostd::args2flag(argv,"--poscar|-poscar")) {
  FLAG_POSCAR=TRUE;
  }
  bool FLAG_POCCSUM=FALSE;
  uint poccsum=0;
  if(argv.size()>=4 && aurostd::args2flag(argv,"--partial_occupation_sum|--pocc_sum|--Li|-Li")) {
  poccsum=(uint) aurostd::args2utype(argv,"--partial_occupation_sum|--pocc_sum|--Li|-Li",(int) 0);
  FLAG_POCCSUM=TRUE;
  cout << "partial_occupation_sum=" << poccsum << endl;
  if((int)poccsum < 0 || poccsum > 8) {
  cerr << "opar_sum in [0,8]" << endl;
  exit(0);
  }
  }
  // priorities
  if(FLAG_LIMNO) FLAG_POCCSUM=FALSE;
  // LOAD
  xstructure str(input,IOAFLOW_AUTO);
  str.BringInCell();
  str.order_parameter_sum=0;
  str.partial_occupation_sum=0.0;
  double eps_partial_occupation=0.1;
  // GET SYMMETRY
  _aflags aflags;aflags.Directory="./";
  aflags.QUIET=TRUE;XHOST.QUIET=TRUE;
  ofstream FileMESSAGE("/dev/null");
  bool PGROUPWRITE=FALSE,FGROUPWRITE=FALSE;
  bool OSSWRITE=TRUE; // to FileMESSAGE, does not matter as it is /dev/null
  SYM::CalculatePointGroup(FileMESSAGE,str,aflags,PGROUPWRITE,OSSWRITE,cout);
  cout << "pgroup operations=" << str.pgroup.size() << endl;
  SYM::CalculateFactorGroup(FileMESSAGE,str,aflags,FGROUPWRITE,OSSWRITE,cout);
  cout << "fgroup operations=" << str.fgroup.size() << endl;

  // fix the order parameter atoms
  for(uint iat=0;iat<str.order_parameter_atoms.size();iat++)
  str.atoms.at(str.order_parameter_atoms.at(iat)).order_parameter_atom=TRUE;

  uint Mn_max=str.order_parameter_atoms.size();
  uint Li_max=str.partial_occupation_atoms.size();
  vector<vector<vector<uint> > > vMn_configurations(Li_max+1);vector<uint> Mn_data(Mn_max+1);
  vector<vector<vector<uint> > > vLi_configurations(Li_max+1);vector<uint> Li_data(Li_max+1);
  vector<vector<uint> > vLiMn_configurations;vector<uint> LiMn_data(Li_max+Mn_max+1+1+1); // sum,orbit,equivalence

  //   // (Li_max+1)*(appropriate size)
  //   for(uint i=0;i<=Li_max;i++) {
  //     uint Li=i,Mn4=Mn_max*4-Mn_max*3-Li,opsum=Mn_max*4-Li;  // disbalance between 3 and 4
  //     double confsMn=((double) fact((double) Mn_max)/fact((double) Mn4)/fact((double) Mn_max-Mn4));
  //     double confsLi=((double) fact((double) Li_max)/fact((double) Li)/fact((double) Li_max-Li));
  //     cerr << Li << " " << Mn_max << " " << Mn4 << " " << opsum << " " << confsMn << endl;
  //     for(uint j=0;j<(uint) confsMn;j++) vMn_configurations.at(Li).push_back(*(new vector<uint>(Mn_max+1)));
  //     for(uint j=0;j<(uint) confsLi;j++) vLi_configurations.at(Li).push_back(*(new vector<uint>(Li_max+1)));
  //   }

  uint maxindex_Mn=(uint) pow(2,Mn_max);
  uint maxindex_Li=(uint) pow(2,Li_max);
  cout << "max_index_Mn=" << maxindex_Mn << "   Mn_max=" << Mn_max << endl;
  cout << "max_index_Li=" << maxindex_Li << "   Li_max=" << Li_max << endl;

  for(uint index_Mn=0;index_Mn<maxindex_Mn;index_Mn++) {
  int Li;
  aurostd::reset(Mn_data);
  for(uint iat=0,uindex_Mn=index_Mn;iat<Mn_max;iat++,uindex_Mn/=2)
  //  Mn_data.at(iat)=(int) 1-2*mod((int) uindex_Mn,2); // flip +3,+4 through the Mn !
  Mn_data.at(iat)=(int) 3+mod((int) uindex_Mn,2); // flip +3,+4 through the Mn !
  Mn_data.at(Mn_max)=aurostd::sum(Mn_data);
  Li=Mn_max*4-Mn_data.at(Mn_max);
  if(Li>=0 && Li<=(int) Li_max) vMn_configurations.at((uint) Li).push_back(Mn_data);
  }
  for(uint index_Li=0;index_Li<maxindex_Li;index_Li++) {
  int Li;
  aurostd::reset(Li_data);
  for(uint iat=0,uindex_Li=index_Li;iat<Li_max;iat++,uindex_Li/=2)
  Li_data.at(iat)=(int) ((double) mod((int) uindex_Li,2)); // flip 0,1
  Li_data.at(Li_max)=aurostd::sum(Li_data);
  Li=Li_data.at(Li_max);
  if(Li>=0 && Li<=(int) Li_max) vLi_configurations.at((uint) Li).push_back(Li_data);
  }
  //  for(uint i=0;i<=Li_max;i++)
  //  cerr << vMn_configurations.at(i).size() << " " << vLi_configurations.at(i).size() << endl;

  uint valueuint;
  for(uint index_Li=0;index_Li<=Li_max;index_Li++)
  if((FLAG_LIMNO && (index_Li==8 || index_Li==4)) || // 8Mn4+ 8Mn3+ fully-Li gs; 12Mn4+ 4Mn3+ half-Li state, metastabless
  (FLAG_POCCSUM && index_Li==poccsum) ||                               // choice of poccsum
  (!FLAG_LIMNO && !FLAG_POCCSUM && (int) index_Li>=0 && index_Li<=8)) {                         // do everything
  cout << "--------------------------------------------------------- "  << endl;
  cout << "spanning Mn configurations = " << vMn_configurations.at(index_Li).size() << endl;
  cout << "spanning Li configurations = " << vLi_configurations.at(index_Li).size() << endl;
  vector<xstructure> vstr;
  xstructure str_ord,str_ord_pocc,str_aus;
  str_ord=str; // str_ord contains some order stuff
  int imax=0;
  for(uint i=0;i<vMn_configurations.at((uint) index_Li).size();i++) {
  // generate structure Mn
  str_ord.order_parameter_sum=0;
  aurostd::reset(LiMn_data);
  for(uint iat=0;iat<Mn_max;iat++) {
  valueuint=vMn_configurations.at((uint) index_Li).at(i).at(iat); // flip +3,+4 through the Mn !
  str_ord.atoms.at(str.order_parameter_atoms.at(iat)).order_parameter_value=valueuint;
  str_ord.atoms.at(str.order_parameter_atoms.at(iat)).order_parameter_atom=TRUE;
  str_ord.order_parameter_sum+=valueuint;
  LiMn_data.at(iat+Li_max)=valueuint;
  }
  // structure created
  // generated structure Li
  for(uint j=0;j<vLi_configurations.at((uint) index_Li).size();j++) {
  str_ord_pocc=str_ord;
  str_ord_pocc.partial_occupation_sum=0.0;	
  for(uint iat=0;iat<Li_max;iat++) {
  valueuint=vLi_configurations.at((uint) index_Li).at(j).at(iat); // flip 0.0, or 1.0
  str_ord_pocc.atoms.at(str.partial_occupation_atoms.at(iat)).partial_occupation_value=valueuint;
  str_ord_pocc.atoms.at(str.partial_occupation_atoms.at(iat)).partial_occupation_atom=TRUE;
  str_ord_pocc.partial_occupation_sum+=valueuint;
  LiMn_data.at(iat)=valueuint;
  }
  LiMn_data.at(Li_max+Mn_max)=0;
  for(uint iat=0;iat<Li_max+0*Mn_max;iat++)
  LiMn_data.at(Li_max+Mn_max)+=LiMn_data.at(iat); // sum
  if(mod((int) ++imax,10)==0) {cout << "*";cout.flush();}// << str_ord.order_parameter_sum;
  // remove atoms
  str.partial_occupation_structure_flag=TRUE;
  if(0) for(int iat=str.partial_occupation_atoms.at(Li_max-1);iat>=0;iat--) {
  str.partial_occupation_structure_flag=FALSE;
  if(aurostd::isequal(str_ord_pocc.atoms.at(str.partial_occupation_atoms.at(iat)).partial_occupation_value,0.0,eps_partial_occupation)) {
  str_ord_pocc.partial_occupation_sum-=str_ord_pocc.atoms.at(str.partial_occupation_atoms.at(iat)).partial_occupation_value;
  LiMn_data.at(Li_max+Mn_max)-=(uint) str_ord_pocc.atoms.at(str.partial_occupation_atoms.at(iat)).partial_occupation_value;
  str_ord_pocc.RemoveAtom(str.partial_occupation_atoms.at(iat));
  }
  }
  // search
  if(LDEBUG) {cerr << "DEBUG  " << str_ord_pocc << endl; exit(0);}
  // check if already in the set
  bool found=FALSE;
  for(uint ifgroup=0;ifgroup<str.fgroup.size()&&found==FALSE;ifgroup++) {                  // cycle through the fgroup of str
  str_aus=SYM_ApplyXstructure(str.fgroup.at(ifgroup),str_ord_pocc,TRUE);                 // slow, do the least number !
  for(uint istr=0;istr<vstr.size()&&found==FALSE;istr++) {                               // cycle through the unique vstrs
  // if(str_ord_pocc.order_parameter_sum==vstr.at(istr).order_parameter_sum && aurostd::isequal(str_ord_pocc.partial_occupation_sum,vstr.at(istr).partial_occupation_sum,eps_partial_occupation))
  {
  if(AreAtomsEquivalent(str_aus.atoms,vstr.at(istr).atoms,TRUE,TRUE)) {
  found=TRUE;
  vstr.at(istr).order_parameter_orbit++;
  LiMn_data.at(Li_max+Mn_max+1)=0; // orbit  FIND
  LiMn_data.at(Li_max+Mn_max+2)=vstr.at(istr).label_uint; // equivalence FIND
  vLiMn_configurations.push_back(LiMn_data);
  vLiMn_configurations.at(vstr.at(istr).label_uint).at(Li_max+Mn_max+1)+=1;
  }
  }
  }
  }
  if(found==FALSE) { // not found, then plug
  str_ord_pocc.order_parameter_orbit=1;
  str_ord_pocc.label_uint=vLiMn_configurations.size(); // osition in the list
  // xstructure *str_ptr;str_ptr = new xstructure(str_ord_pocc);vstr.push_back(*str_ptr);
  vstr.push_back(str_ord_pocc);
  LiMn_data.at(Li_max+Mn_max+1)=1; // orbit
  LiMn_data.at(Li_max+Mn_max+2)=vLiMn_configurations.size(); // equivalence
  vLiMn_configurations.push_back(LiMn_data);
  cout << vstr.size() << "(" << i << "/" << vMn_configurations.at((uint) index_Li).size() << ") ";cout.flush();
  }
  // search done
  //	cout << "[" <<  index_Li << " " << j << " " << i << "] ";
  } // order parameter cycle
  } // partial occupation cycle
    // FINISHED, Got them all, now order them
    cout << endl;
    cout << "unique orbits = " << vstr.size() << endl;
    if(1) {
    ofstream FileOUT(string("configurations.Li"+int2string((int) index_Li)+"Mn16O32.data").c_str());
    bool found_unique;
    for(uint i=0;i<vLiMn_configurations.size();i++) {
    FileOUT << aurostd::PaddedPRE(aurostd::utype2string(i),5) << " ";
    for(uint j=0;j<Li_max+Mn_max;j++)
    FileOUT << vLiMn_configurations.at(i).at(j) << " ";
    FileOUT << vLiMn_configurations.at(i).at(Li_max+Mn_max) << " ";
    FileOUT << aurostd::PaddedPRE(aurostd::utype2string(vLiMn_configurations.at(i).at(Li_max+Mn_max+1)),4) << " ";
    FileOUT << aurostd::PaddedPRE(aurostd::utype2string(vLiMn_configurations.at(i).at(Li_max+Mn_max+2)),4) << " ";
    found_unique=FALSE;
    for(uint istr=0;istr<vstr.size();istr++)
    if(vstr.at(istr).label_uint==i) {
    FileOUT << "  t" << aurostd::PaddedPRE(aurostd::utype2string(istr),4,"0") << "   E    ";
    FileOUT <<string("  RUN_os")+PaddedPRE(vstr.at(istr).order_parameter_sum,2,"0")+string("_po")+PaddedPRE(vstr.at(istr).partial_occupation_sum,2,"0")+string("_od")+PaddedPRE(int2string(vstr.at(istr).order_parameter_orbit),3,"0")+"_t"+PaddedPRE(int2string(istr),4,"0");
    found_unique=TRUE;
    FileOUT << endl;
    }
    if(found_unique==FALSE) FileOUT << "          -    " << endl;
    }
    FileOUT.flush();FileOUT.close();
    }
    //    sort(vstr.begin(),vstr.end(),_rsort_xstructure_order_parameter_orbit());
    pflow::VSTR_Write("./",vstr,FLAG_POSCAR,VSTR_WRITE_POCC_OPAR);
    } // Li cycle
    return TRUE;
    }

    // ******************** SAFE
    bool SHIRLEY2(vector<string> argv,istream& input) {
    cout << "--------------------------------------------------------- "  << endl;
    cout << aflow::Banner("BANNER_TINY") << endl;
    cout << "argv.size()=" << argv.size() << endl;
    // FLAGS
    bool FLAG_LIMNO=FALSE,FLAG_POSCAR=FALSE;
    int partial_occupation_total_atoms=0;
    if(argv.size()>=3 && aurostd::args2flag(argv,"--Li8Mn16O32|--LiMnO")) {
    FLAG_LIMNO=TRUE;
    }
    if(argv.size()>=3 && aurostd::args2flag(argv,"--poscar")) {
    FLAG_POSCAR=TRUE;
    }
    bool FLAG_OPARSUM=FALSE;int oparsum=0;
    if(argv.size()>=4 && aurostd::args2flag(argv,"--opar_sum|--order_parameter_sum|--oparameter_sum")) {
    oparsum=aurostd::args2utype(argv,"--opar_sum|--order_parameter_sum|--oparameter_sum",(int) 0);
    FLAG_OPARSUM=TRUE;
    cout << "order_parameter_sum=" << oparsum << endl;
    partial_occupation_total_atoms=64-oparsum;
    cout << "partial_occupation_total_atoms=" << partial_occupation_total_atoms << endl;
    if(oparsum < 56 || oparsum > 64) {
    cerr << "opar_sum in [56,64]" << endl;
    exit(0);
    }
    }
    // priorities
    if(FLAG_LIMNO) FLAG_OPARSUM=FALSE;
    // LOAD
    xstructure str(input,IOAFLOW_AUTO);
    str.BringInCell();
    str.order_parameter_sum=0;
    str.partial_occupation_sum=0.0;
    double eps_partial_occupation=0.1;
    // GET SYMMETRY
    _aflags aflags;aflags.Directory="./";
    aflags.QUIET=TRUE;XHOST.QUIET=TRUE;
    ofstream FileMESSAGE("/dev/null");
    bool PGROUPWRITE=FALSE,FGROUPWRITE=FALSE;
    bool OSSWRITE=FALSE; // to FileMESSAGE, does not matter as it is /dev/null
    SYM::CalculatePointGroup(FileMESSAGE,str,aflags,PGROUPWRITE,OSSWRITE,cout);
    cout << "pgroup operations=" << str.pgroup.size() << endl;
    SYM::CalculateFactorGroup(FileMESSAGE,str,aflags,FGROUPWRITE,OSSWRITE,cout);
    cout << "fgroup operations=" << str.fgroup.size() << endl;

    // fix the order parameter atoms
    for(uint iat=0;iat<str.order_parameter_atoms.size();iat++)
    str.atoms.at(str.order_parameter_atoms.at(iat)).order_parameter_atom=TRUE;

    // for(uint op=56;op<=64;op++) // hack to make them all
    {
    //    oparsum=op;  FLAG_OPARSUM=TRUE; // hack to make them all
    cout << "--------------------------------------------------------- "  << endl;

    vector<xstructure> vstr;
    xstructure str_ord,str_ord_pocc,str_aus;
    uint maxindex_Mn=(uint) pow(2,str.order_parameter_atoms.size());
    uint maxindex_Li=(uint) pow(2,str.partial_occupation_atoms.size());
    cout << "max_index_Mn=" << maxindex_Mn << "   str.order_parameter_atoms.size())=" << str.order_parameter_atoms.size() << endl;
    cout << "max_index_Li=" << maxindex_Li << "   str.partial_occupation_atoms.size())=" << str.partial_occupation_atoms.size() << endl;
    str_ord=str; // str_ord contains some order stuff
    bool found=FALSE;

    for(uint index_Mn=0;index_Mn<maxindex_Mn;index_Mn++) {
    // generate structure
    str_ord.order_parameter_sum=0;
    for(uint iat=0,uindex_Mn=index_Mn;iat<str.order_parameter_atoms.size();iat++,uindex_Mn/=2) {
    //  str_ord.atoms.at(str.order_parameter_atoms.at(iat)).order_parameter_value=(int) 1-2*mod((int) uindex_Mn,2); // flip +3,+4 through the Mn !
    str_ord.atoms.at(str.order_parameter_atoms.at(iat)).order_parameter_value=(int) 3+mod((int) uindex_Mn,2); // flip +3,+4 through the Mn !
    str_ord.atoms.at(str.order_parameter_atoms.at(iat)).order_parameter_atom=TRUE;
    str_ord.order_parameter_sum+=str_ord.atoms.at(str.order_parameter_atoms.at(iat)).order_parameter_value;
    partial_occupation_total_atoms=64-str_ord.order_parameter_sum; // 60=>4, 56=>8  (dum from 56 (full Li) to 64 (fully de-Li))
    }
    // structure created
    if((FLAG_LIMNO && (str_ord.order_parameter_sum==56 || str_ord.order_parameter_sum==60)) || // 8Mn4+ 8Mn3+ fully-Li gs; 12Mn4+ 4Mn3+ half-Li state, metastabless
    (FLAG_OPARSUM && str_ord.order_parameter_sum==oparsum) ||                               // choice of oparsum
    (!FLAG_LIMNO && !FLAG_OPARSUM && str_ord.order_parameter_sum>=56 && str_ord.order_parameter_sum<=64)) {                         // do everything
    //	cout << "[" << index_Mn << "] "; cout.flush();
    // now I need to bust out half Li
    for(uint index_Li=0;index_Li<maxindex_Li;index_Li++) {
    // create structure from index_Li
    str_ord_pocc=str_ord;
    str_ord_pocc.partial_occupation_sum=0.0;
    for(uint iat=0,uindex_Li=index_Li;iat<str.partial_occupation_atoms.size();iat++,uindex_Li/=2) {
    str_ord_pocc.atoms.at(str.partial_occupation_atoms.at(iat)).partial_occupation_value=(double) mod((int) uindex_Li,2); // flip 0.0, or 1.0
    str_ord_pocc.atoms.at(str.partial_occupation_atoms.at(iat)).partial_occupation_atom=TRUE;
    str_ord_pocc.partial_occupation_sum+=str_ord_pocc.atoms.at(str.partial_occupation_atoms.at(iat)).partial_occupation_value;
    }
    if(aurostd::isequal(str_ord_pocc.partial_occupation_sum, (double) partial_occupation_total_atoms,eps_partial_occupation)) {
    // cout << "str_ord_pocc.order_parameter_sum=" << str_ord_pocc.order_parameter_sum << endl;
    // cout << "str_ord_pocc.partial_occupation_sum=" << str_ord_pocc.partial_occupation_sum << endl;
    // remove atoms
    for(int iat=str.partial_occupation_atoms.at(str.partial_occupation_atoms.size()-1);iat>=0;iat--) {
    if(aurostd::isequal(str_ord_pocc.atoms.at(str.partial_occupation_atoms.at(iat)).partial_occupation_value,0.0,eps_partial_occupation)) {
    str_ord_pocc.partial_occupation_sum-=str_ord_pocc.atoms.at(str.partial_occupation_atoms.at(iat)).partial_occupation_value;
    str_ord_pocc.RemoveAtom(str.partial_occupation_atoms.at(iat));
    }
    }
    if(LDEBUG) {cerr << "DEBUG  " << str_ord_pocc << endl; exit(0);}
    // check if already in the set
    found=FALSE;
    for(uint ifgroup=0;ifgroup<str.fgroup.size()&&found==FALSE;ifgroup++) {                  // cycle through the fgroup of str
    str_aus=SYM_ApplyXstructure(str.fgroup.at(ifgroup),str_ord_pocc,TRUE);                 // slow, do the least number !
    for(uint istr=0;istr<vstr.size()&&found==FALSE;istr++) {                               // cycle through the unique vstrs
    if(str_ord_pocc.order_parameter_sum==vstr.at(istr).order_parameter_sum &&
    aurostd::isequal(str_ord_pocc.partial_occupation_sum,vstr.at(istr).partial_occupation_sum,eps_partial_occupation)) {
    if(AreAtomsEquivalent(str_aus.atoms,vstr.at(istr).atoms,TRUE,TRUE)) {
    found=TRUE;
    vstr.at(istr).order_parameter_orbit++;
    }
    }
    }
    }
    if(found==FALSE) { // not found, then plug
    str_ord_pocc.order_parameter_orbit=1;
    // xstructure *str_ptr;str_ptr = new xstructure(str_ord_pocc);vstr.push_back(*str_ptr);
    vstr.push_back(str_ord_pocc);
    // done
    cout << vstr.size() << "(" << index_Mn << "," << index_Li << ") ";cout.flush();
    }
    }
    }
    }
    if(mod((int) index_Mn,1000)==0) {cout << "*";cout.flush();}// << str_ord.order_parameter_sum;
    }
    // FINISHED, Got them all, now order them
    cout << endl << "unique orbits = " << vstr.size() << endl;
    sort(vstr.begin(),vstr.end(),_rsort_xstructure_order_parameter_orbit());
    pflow::VSTR_Write("./",vstr,FLAG_POSCAR,VSTR_WRITE_POCC_OPAR);
    }
    return TRUE;
    }

    // ***************************************************************************
    // pflow::POCCUPATION
    // ***************************************************************************
    class _sort_xstructure_partial_occupation_orbit {                    // sorting through reference
    public:
    bool operator()(const xstructure& str1, const xstructure& str2) const {
    return (bool) (str1.partial_occupation_orbit<str2.partial_occupation_orbit);}
    };

    class _rsort_xstructure_partial_occupation_orbit {                    // sorting through reference
    public:
    bool operator()(const xstructure& str1, const xstructure& str2) const {
    return (bool) (str1.partial_occupation_orbit>str2.partial_occupation_orbit);}
    };

    bool POCCUPATION(vector<string> argv,istream& input) {
    cout << "--------------------------------------------------------- "  << endl;
    cout << aflow::Banner("BANNER_TINY") << endl;
    cout << "argv.size()=" << argv.size() << endl;
    // LOAD
    xstructure str(input,IOAFLOW_AUTO);
    str.BringInCell();
    str.partial_occupation_sum=0;
    // FLAGS
    bool FLAG_POSCAR=FALSE;
    if(argv.size()>=3 && aurostd::args2flag(argv,"--poscar")) {
    FLAG_POSCAR=TRUE;
    if(FLAG_POSCAR) cout << "dummy: FLAG_POSCAR" << endl;
    }
    bool FLAG_POCCSUM=FALSE;int poccsum=0;
    if(argv.size()>=4 && aurostd::args2flag(argv,"--pocc_sum|--partial_occupation_sum|--poccameter_sum")) {
    poccsum=aurostd::args2utype(argv,"--pocc_sum|--partial_occupation_sum|--poccameter_sum",(int) 0);
    FLAG_POCCSUM=TRUE;
    cout << "partial_occupation_sum=" << poccsum << endl;
    if(poccsum < 0 || poccsum > (int) str.partial_occupation_atoms.size()) {
    cerr << "pocc_sum in [0," << str.partial_occupation_atoms.size() << "]" << endl;
    exit(0);
    }
    }
    // GET SYMMETRY
    _aflags aflags;aflags.Directory="./";
    aflags.QUIET=TRUE;XHOST.QUIET=TRUE;
    ofstream FileMESSAGE("/dev/null");
    bool PGROUPWRITE=FALSE,FGROUPWRITE=FALSE;//,SGROUPWRITE=FALSE,IATOMSWRITE=TRUE,AGROUPWRITE=TRUE;
    bool OSSWRITE=FALSE; // to FileMESSAGE, does not matter as it is /dev/null
    //  SGROUPWRITE=FALSE,IATOMSWRITE=TRUE,AGROUPWRITE=TRUE;
    SYM::CalculatePointGroup(FileMESSAGE,str,aflags,PGROUPWRITE,OSSWRITE,cout);
    cout << "pgroup operations=" << str.pgroup.size() << endl;
    SYM::CalculateFactorGroup(FileMESSAGE,str,aflags,FGROUPWRITE,OSSWRITE,cout);
    cout << "fgroup operations=" << str.fgroup.size() << endl;
    // SYM::CalculateSpaceGroup(FileMESSAGE,str,aflags,SGROUPWRITE,OSSWRITE,cout);
    //  cout << "sgroup operations=" << str.sgroup.size() << endl;
    // SYM::CalculateInequivalentAtoms(FileMESSAGE,str,aflags,IATOMSWRITE,OSSWRITE,cout);
    // SYM::CalculateSitePointGroup(FileMESSAGE,str,aflags,AGROUPWRITE,OSSWRITE,cout);
    //  cout << str.partial_occupation_atoms.size() << endl;
    //for(uint i=0;i<str.partial_occupation_atoms.size();i++)   cout << str.partial_occupation_atoms.at(i) << endl;

    // fix the partial occupation atoms
    for(uint iat=0;iat<str.partial_occupation_atoms.size();iat++)
    str.atoms.at(str.partial_occupation_atoms.at(iat)).partial_occupation_atom=TRUE;

    vector<xstructure> vstr;
    xstructure str_ord,str_aus;
    uint maxindex=(uint) pow(2,str.partial_occupation_atoms.size());

    cout << "max_index=" << maxindex << endl;
    str_ord=str; // str_ord contains some partial stuff
    for(uint index=0;index<maxindex;index++) {
    // generate structure
    //   str_ord=str; // str_ord contains some partial stuff
    // str_ord.title+=" ["+int2string(index)+"]";
    // str_ord.BringInCell(); // not really necessary
    str_ord.partial_occupation_sum=0;
    for(uint iat=0,uindex=index;iat<str.partial_occupation_atoms.size();iat++,uindex/=2) {
    //  str_ord.atoms.at(str.partial_occupation_atoms.at(iat)).partial_occupation_value=(int) 1-2*mod((int) uindex,2); // flip +3,+4 through the Mn !
    //      str_ord.atoms.at(str.partial_occupation_atoms.at(iat)).partial_occupation_value=(int) 3+mod((int) uindex,2); // flip +3,+4 through the Mn !
    str_ord.atoms.at(str.partial_occupation_atoms.at(iat)).partial_occupation_value=(int) 0+mod((int) uindex,2); // flip +3,+4 through the Mn !
    str_ord.atoms.at(str.partial_occupation_atoms.at(iat)).partial_occupation_atom=TRUE;
    str_ord.partial_occupation_sum+=str_ord.atoms.at(str.partial_occupation_atoms.at(iat)).partial_occupation_value;
    }
    if(FLAG_POCCSUM==FALSE || (FLAG_POCCSUM==TRUE && str_ord.partial_occupation_sum==poccsum)) {
    //   str_ord.partial_occupation_sum=partial_occupation_sum(str_ord);
    // structure created
    bool found=FALSE;
    for(uint ifgroup=0;ifgroup<str.fgroup.size()&&found==FALSE;ifgroup++) {                  // cycle through the fgroup of str
    str_aus=SYM_ApplyXstructure(str.fgroup.at(ifgroup),str_ord,TRUE);
    for(uint istr=0;istr<vstr.size()&&found==FALSE;istr++) {                               // cycle through the unique vstrs
    if(str_ord.partial_occupation_sum==vstr.at(istr).partial_occupation_sum)
    if(AreAtomsEquivalent(str_aus.atoms,vstr.at(istr).atoms,FALSE,TRUE)) {
    found=TRUE;
    vstr.at(istr).partial_occupation_orbit++;
    }
    }
    }
    if(found==FALSE) { // not found, then plug
    str_ord.partial_occupation_orbit=1;
    // xstructure *str_ptr;str_ptr = new xstructure(str_ord);vstr.push_back(*str_ptr);
    vstr.push_back(str_ord);
    cout << vstr.size() << "(" << index << ") ";
    cout.flush();
    }
    if(mod((int) index,1000)==0) {cout << "*";cout.flush();}// << str_ord.partial_occupation_sum;
    } // 3
    }
    cout << endl << "unique orbits = " << vstr.size() << endl;
    sort(vstr.begin(),vstr.end(),_rsort_xstructure_partial_occupation_orbit());

    uint index=0;
    for(uint istr=0;istr<vstr.size();istr++) {                       // cycle through the unique vstrs
    vstr.at(istr).title+=" [orbits="+int2string(vstr.at(istr).partial_occupation_orbit)+"]";
    index+=vstr.at(istr).partial_occupation_orbit;
    cout << "************************************************************************************" << endl;
    cout << "info=" << istr << "  " << vstr.at(istr).partial_occupation_orbit << " " << vstr.at(istr).partial_occupation_sum << endl;
    cout << vstr.at(istr) << endl;
    cout << "************************************************************************************" << endl;
    //    if(vstr.at(istr).partial_occupation_orbit==vstr.at(0).partial_occupation_orbit)
    if(0) {
    _xvasp xvasp;
    xvasp.AVASP_dirbase="./";
    xvasp.AVASP_label=string("RUN_os")+PaddedPRE(vstr.at(istr).partial_occupation_sum,2,"0")+string("_od")+PaddedPRE(int2string(vstr.at(istr).partial_occupation_orbit),3,"0")+"_t"+PaddedPRE(int2string(istr),3,"0");
    xvasp.str=vstr.at(istr);

    xvasp.str.species.clear();xvasp.str.species_pp.clear();xvasp.str.species_pp_type.clear();xvasp.str.species_pp_version.clear();xvasp.str.species_pp_ZVAL.clear();xvasp.str.species_pp_vLDAU.clear();
    xvasp.str.species.push_back("Li_sv");xvasp.str.species.push_back("Mn");xvasp.str.species.push_back("O_s");
    xvasp.str.species_pp=xvasp.str.species;
    xvasp.str.species_pp_type.push_back("");xvasp.str.species_pp_type.push_back("");xvasp.str.species_pp_type.push_back("");
    xvasp.str.species_pp_version.push_back("");xvasp.str.species_pp_version.push_back("");xvasp.str.species_pp_version.push_back("");
    xvasp.str.species_pp_ZVAL.push_back(0.0);xvasp.str.species_pp_ZVAL.push_back(0.0);xvasp.str.species_pp_ZVAL.push_back(0.0);
    xvasp.str.species_pp_vLDAU.push_back(deque<double>());xvasp.str.species_pp_vLDAU.push_back(deque<double>());xvasp.str.species_pp_vLDAU.push_back(deque<double>());

    xvasp.AVASP_volume_in=558.7229;
    xvasp.str.SetVolume(xvasp.AVASP_volume_in);
    xvasp.AVASP_prototype_from_library_=FALSE;
    //cerr << xvasp.AVASP_label << endl;
    //	AVASP_MakeSinglePOSCAR(xvasp);
    }
    }
    return TRUE;
    }
*/

// ***************************************************************************
// helpIndividualOption
// ***************************************************************************
void helpIndividualOption(vector<string> & argv) {
  // "aflow --help option"
  // to print out the help information of given option
  // use the information stored in README_AFLOW_PFLOW.txt

  stringstream help_file;
  help_file << init::InitGlobalObject("README_AFLOW_PFLOW_TXT");

  bool help_found=false;
  string option = argv.at(2);

  const string seperation_line_pattern = "******";

  string line;

  // to get standard form of option --option+" "
  option.append(" ");
  option.insert(0, "--");

  bool output_flag=false;
  bool begin_section_flag=false;
  while(getline(help_file, line)) {

    // get the begining of the help section of an option
    // it begins with "aflow --option_name"
    // it may have one or several space between the two words

    size_t pos = line.find(option);
    string program_name = "aflow";
    string option_prefix = "--";
    size_t pos2 = line.find(program_name);
    size_t pos3 = line.find(option_prefix);
    string::iterator str_itr;

    if((pos2 < pos3) &&(pos3 != string::npos)) {
      for(str_itr = line.begin() + pos2 + program_name.size();
	  str_itr < line.begin() + pos3; str_itr++) {
	if(*str_itr != ' ') {
	  begin_section_flag=false;
	  break;
	}
      }
      begin_section_flag=true;
    } else {
      begin_section_flag=false;
    }

    if(begin_section_flag && (pos != string::npos) && (! output_flag)) {
      // find the line with --option
      // output this line and lines below
      // until next line begining with "aflow --" pattern

      int last_space_pos = line.find_first_not_of(" ");
      line.erase(line.begin(), line.begin()+last_space_pos);
      cout << line << endl;
      output_flag=true;
      help_found=true;
      continue;
    } else if(begin_section_flag && output_flag) {
      // reset the flag and output
      // the other part of help information if found

      pos = line.find(option);
      if((pos != string::npos)) {
	int last_space_pos = line.find_first_not_of(" ");
	line.erase(line.begin(), line.begin()+last_space_pos);
	cout << line << endl;
	output_flag=true;
      } else {
	output_flag=false;
      }
    } else if(!begin_section_flag && output_flag) {
      // output the part of help section not
      // beginning with "aflow --option" pattern

      if(line.find(seperation_line_pattern) == string::npos) {
	// do not output seperation line
	int last_space_pos = line.find_first_not_of(" ");
	line.erase(line.begin(), line.begin()+last_space_pos);
	cout << "    " << line << endl;
      }
    }
  }

  if(!help_found) {
    // not supported option
    size_t pos = option.find_last_of(" ");
    string option_raw = option.substr(2, pos-2);
    cerr << "No option \"" << option_raw << "\" is found." << endl;
    cerr << "Try \"aflow --help\" to see all supported options." << endl;
    exit(1);
  } else {
    exit(0);
  }
}

// ***************************************************************************

// ----------------------------------------------------------------------------
// get_itemized_vector_string_from_input
// Stefano Curtarolo
// namespace aurostd {
//   bool get_itemized_vector_string_from_input(vector<string> &argv,const string& s0,vector<string>& tokens,const string& delimiter) {// =":") {
//     string s0neq=s0,s0equ;aurostd::StringSubst(s0neq,"=","");s0equ=s0neq+"=";
//     if(aurostd::args2attachedflag(argv,s0equ)) {aurostd::string2tokens(aurostd::args2attachedstring(argv,s0equ,EMPTY_WORDING),tokens,delimiter);}
//     if(aurostd::args2flag(argv,s0neq)) {aurostd::args2string(argv,s0neq,EMPTY_WORDING);tokens=aurostd::args2vectorstring(argv,s0neq,EMPTY_WORDING);}
//     if(tokens.size()==1 && aurostd::substring2bool(tokens.at(0),delimiter)) {s0equ=tokens.at(0);aurostd::string2tokens(s0equ,tokens,delimiter);}
//     if(tokens.size()==0) return FALSE;
//     return TRUE;
//   }
//   bool get_itemized_vector_string_from_input(vector<string> &argv,const string& s0,const string& s1,vector<string>& tokens,const string& delimiter) {// =":") {
//     string s0neq=s0,s0equ;aurostd::StringSubst(s0neq,"=","");s0equ=s0neq+"=";
//     string s1neq=s1,s1equ;aurostd::StringSubst(s1neq,"=","");s1equ=s1neq+"=";
//     if(aurostd::args2attachedflag(argv,s0equ)) {aurostd::string2tokens(aurostd::args2attachedstring(argv,s0equ,EMPTY_WORDING),tokens,delimiter);}
//     if(aurostd::args2attachedflag(argv,s1equ)) {aurostd::string2tokens(aurostd::args2attachedstring(argv,s1equ,EMPTY_WORDING),tokens,delimiter);}
//     if(aurostd::args2flag(argv,s0neq)) {aurostd::args2string(argv,s0neq,EMPTY_WORDING);tokens=aurostd::args2vectorstring(argv,s0neq,EMPTY_WORDING);}
//     if(aurostd::args2flag(argv,s1neq)) {aurostd::args2string(argv,s1neq,EMPTY_WORDING);tokens=aurostd::args2vectorstring(argv,s1neq,EMPTY_WORDING);}
//     if(tokens.size()==1 && aurostd::substring2bool(tokens.at(0),delimiter)) {s0equ=tokens.at(0);aurostd::string2tokens(s0equ,tokens,delimiter);}
//     if(tokens.size()==0) return FALSE;
//     return TRUE;
//   }
//   bool get_itemized_vector_string_from_input(vector<string> &argv,const string& s0,const string& s1,const string& s2,vector<string>& tokens,const string& delimiter) {// =":") {
//     string s0neq=s0,s0equ;aurostd::StringSubst(s0neq,"=","");s0equ=s0neq+"=";
//     string s1neq=s1,s1equ;aurostd::StringSubst(s1neq,"=","");s1equ=s1neq+"=";
//     string s2neq=s2,s2equ;aurostd::StringSubst(s2neq,"=","");s2equ=s2neq+"=";
//     if(aurostd::args2attachedflag(argv,s0equ)) {aurostd::string2tokens(aurostd::args2attachedstring(argv,s0equ,EMPTY_WORDING),tokens,delimiter);}
//     if(aurostd::args2attachedflag(argv,s1equ)) {aurostd::string2tokens(aurostd::args2attachedstring(argv,s1equ,EMPTY_WORDING),tokens,delimiter);}
//     if(aurostd::args2attachedflag(argv,s2equ)) {aurostd::string2tokens(aurostd::args2attachedstring(argv,s2equ,EMPTY_WORDING),tokens,delimiter);}
//     if(aurostd::args2flag(argv,s0neq)) {aurostd::args2string(argv,s0neq,EMPTY_WORDING);tokens=aurostd::args2vectorstring(argv,s0neq,EMPTY_WORDING);}
//     if(aurostd::args2flag(argv,s1neq)) {aurostd::args2string(argv,s1neq,EMPTY_WORDING);tokens=aurostd::args2vectorstring(argv,s1neq,EMPTY_WORDING);}
//     if(aurostd::args2flag(argv,s2neq)) {aurostd::args2string(argv,s2neq,EMPTY_WORDING);tokens=aurostd::args2vectorstring(argv,s2neq,EMPTY_WORDING);}
//     if(tokens.size()==1 && aurostd::substring2bool(tokens.at(0),delimiter)) {s0equ=tokens.at(0);aurostd::string2tokens(s0equ,tokens,delimiter);}
//     if(tokens.size()==0) return FALSE;
//     return TRUE;
//   }
//   bool get_itemized_vector_string_from_input(vector<string> &argv,const string& s0,const string& s1,const string& s2,const string& s3,vector<string>& tokens,const string& delimiter) {// =":") {
//     string s0neq=s0,s0equ;aurostd::StringSubst(s0neq,"=","");s0equ=s0neq+"=";
//     string s1neq=s1,s1equ;aurostd::StringSubst(s1neq,"=","");s1equ=s1neq+"=";
//     string s2neq=s2,s2equ;aurostd::StringSubst(s2neq,"=","");s2equ=s2neq+"=";
//     string s3neq=s3,s3equ;aurostd::StringSubst(s3neq,"=","");s3equ=s3neq+"=";
//     if(aurostd::args2attachedflag(argv,s0equ)) {aurostd::string2tokens(aurostd::args2attachedstring(argv,s0equ,EMPTY_WORDING),tokens,delimiter);}
//     if(aurostd::args2attachedflag(argv,s1equ)) {aurostd::string2tokens(aurostd::args2attachedstring(argv,s1equ,EMPTY_WORDING),tokens,delimiter);}
//     if(aurostd::args2attachedflag(argv,s2equ)) {aurostd::string2tokens(aurostd::args2attachedstring(argv,s2equ,EMPTY_WORDING),tokens,delimiter);}
//     if(aurostd::args2attachedflag(argv,s3equ)) {aurostd::string2tokens(aurostd::args2attachedstring(argv,s3equ,EMPTY_WORDING),tokens,delimiter);}
//     if(aurostd::args2flag(argv,s0neq)) {aurostd::args2string(argv,s0neq,EMPTY_WORDING);tokens=aurostd::args2vectorstring(argv,s0neq,EMPTY_WORDING);}
//     if(aurostd::args2flag(argv,s1neq)) {aurostd::args2string(argv,s1neq,EMPTY_WORDING);tokens=aurostd::args2vectorstring(argv,s1neq,EMPTY_WORDING);}
//     if(aurostd::args2flag(argv,s2neq)) {aurostd::args2string(argv,s2neq,EMPTY_WORDING);tokens=aurostd::args2vectorstring(argv,s2neq,EMPTY_WORDING);}
//     if(aurostd::args2flag(argv,s3neq)) {aurostd::args2string(argv,s3neq,EMPTY_WORDING);tokens=aurostd::args2vectorstring(argv,s3neq,EMPTY_WORDING);}
//     if(tokens.size()==1 && aurostd::substring2bool(tokens.at(0),delimiter)) {s0equ=tokens.at(0);aurostd::string2tokens(s0equ,tokens,delimiter);}
//     if(tokens.size()==0) return FALSE;
//     return TRUE;
//   }
// }

// // ----------------------------------------------------------------------------
// // getproto_itemized_vector_string_from_input
// // Stefano Curtarolo
// namespace aurostd {
//   bool getproto_itemized_vector_string_from_input(vector<string> &argv,const string& s0,vector<string>& tokens,const string& delimiter) {// =":") {
//     string ss;tokens.clear();vector<string> stokens;
//     string s0neq=s0,s0equ;aurostd::StringSubst(s0neq,"=","");s0equ=s0neq+"=";
//     if(aurostd::args2attachedflag(argv,s0equ)) ss=aurostd::args2attachedstring(argv,s0equ,EMPTY_WORDING);
//     if(aurostd::args2flag(argv,s0neq)) ss=aurostd::args2string(argv,s0neq,EMPTY_WORDING);
//     if(aurostd::substring2bool(ss,"./")) aurostd::StringSubst(ss,"./","");
//     if(ss!="") {
//       if(!aurostd::substring2bool(ss,"/")) { return get_itemized_vector_string_from_input(argv,s0,tokens,delimiter);
//       } else {
// 	aurostd::string2tokens(ss,stokens,"/");
// 	tokens.push_back(stokens.at(1));
// 	KBIN::VASP_SplitAlloyPseudoPotentials(stokens.at(0),stokens);
// 	for(uint i=0;i<stokens.size();i++) tokens.push_back(stokens.at(i));
//       }
//     }
//     if(tokens.size()==1 && aurostd::substring2bool(tokens.at(0),delimiter)) {s0equ=tokens.at(0);aurostd::string2tokens(s0equ,tokens,delimiter);}
//     if(tokens.size()==0) return FALSE;
//     return TRUE;
//   }
//   bool getproto_itemized_vector_string_from_input(vector<string> &argv,const string& s0,const string& s1,vector<string>& tokens,const string& delimiter) {// =":") {
//     string ss;tokens.clear();vector<string> stokens;
//     string s0neq=s0,s0equ;aurostd::StringSubst(s0neq,"=","");s0equ=s0neq+"=";
//     string s1neq=s1,s1equ;aurostd::StringSubst(s1neq,"=","");s1equ=s1neq+"=";
//     if(aurostd::args2attachedflag(argv,s0equ)) ss=aurostd::args2attachedstring(argv,s0equ,EMPTY_WORDING);
//     if(aurostd::args2flag(argv,s0neq)) ss=aurostd::args2string(argv,s0neq,EMPTY_WORDING);
//     if(aurostd::args2attachedflag(argv,s1equ)) ss=aurostd::args2attachedstring(argv,s1equ,EMPTY_WORDING);
//     if(aurostd::args2flag(argv,s1neq)) ss=aurostd::args2string(argv,s1neq,EMPTY_WORDING);
//     if(aurostd::substring2bool(ss,"./")) aurostd::StringSubst(ss,"./","");
//     if(ss!="") {
//       if(!aurostd::substring2bool(ss,"/")) { return get_itemized_vector_string_from_input(argv,s0,s1,tokens,delimiter);
//       } else {
// 	aurostd::string2tokens(ss,stokens,"/");
// 	tokens.push_back(stokens.at(1));
// 	KBIN::VASP_SplitAlloyPseudoPotentials(stokens.at(0),stokens);
// 	for(uint i=0;i<stokens.size();i++) tokens.push_back(stokens.at(i));
//       }
//     }
//     if(tokens.size()==1 && aurostd::substring2bool(tokens.at(0),delimiter)) {s0equ=tokens.at(0);aurostd::string2tokens(s0equ,tokens,delimiter);}
//     if(tokens.size()==0) return FALSE;
//     return TRUE;
//   }
//   bool getproto_itemized_vector_string_from_input(vector<string> &argv,const string& s0,const string& s1,const string& s2,vector<string>& tokens,const string& delimiter) {// =":") {
//     string ss;tokens.clear();vector<string> stokens;
//     string s0neq=s0,s0equ;aurostd::StringSubst(s0neq,"=","");s0equ=s0neq+"=";
//     string s1neq=s1,s1equ;aurostd::StringSubst(s1neq,"=","");s1equ=s1neq+"=";
//     string s2neq=s2,s2equ;aurostd::StringSubst(s2neq,"=","");s2equ=s2neq+"=";
//     if(aurostd::args2attachedflag(argv,s0equ)) ss=aurostd::args2attachedstring(argv,s0equ,EMPTY_WORDING);
//     if(aurostd::args2flag(argv,s0neq)) ss=aurostd::args2string(argv,s0neq,EMPTY_WORDING);
//     if(aurostd::args2attachedflag(argv,s1equ)) ss=aurostd::args2attachedstring(argv,s1equ,EMPTY_WORDING);
//     if(aurostd::args2flag(argv,s1neq)) ss=aurostd::args2string(argv,s1neq,EMPTY_WORDING);
//     if(aurostd::args2attachedflag(argv,s2equ)) ss=aurostd::args2attachedstring(argv,s2equ,EMPTY_WORDING);
//     if(aurostd::args2flag(argv,s2neq)) ss=aurostd::args2string(argv,s2neq,EMPTY_WORDING);
//     if(aurostd::substring2bool(ss,"./")) aurostd::StringSubst(ss,"./","");
//     if(ss!="") {
//       if(!aurostd::substring2bool(ss,"/")) { return get_itemized_vector_string_from_input(argv,s0,s1,s2,tokens,delimiter);
//       } else {
// 	aurostd::string2tokens(ss,stokens,"/");
// 	tokens.push_back(stokens.at(1));
// 	KBIN::VASP_SplitAlloyPseudoPotentials(stokens.at(0),stokens);
// 	for(uint i=0;i<stokens.size();i++) tokens.push_back(stokens.at(i));
//       }
//     }
//     if(tokens.size()==1 && aurostd::substring2bool(tokens.at(0),delimiter)) {s0equ=tokens.at(0);aurostd::string2tokens(s0equ,tokens,delimiter);}
//     if(tokens.size()==0) return FALSE;
//     return TRUE;
//   }
//   bool getproto_itemized_vector_string_from_input(vector<string> &argv,const string& s0,const string& s1,const string& s2,const string& s3,vector<string>& tokens,const string& delimiter) {// =":") {
//     string ss;tokens.clear();vector<string> stokens;
//     string s0neq=s0,s0equ;aurostd::StringSubst(s0neq,"=","");s0equ=s0neq+"=";
//     string s1neq=s1,s1equ;aurostd::StringSubst(s1neq,"=","");s1equ=s1neq+"=";
//     string s2neq=s2,s2equ;aurostd::StringSubst(s2neq,"=","");s2equ=s2neq+"=";
//     string s3neq=s3,s3equ;aurostd::StringSubst(s3neq,"=","");s3equ=s3neq+"=";
//     if(aurostd::args2attachedflag(argv,s0equ)) ss=aurostd::args2attachedstring(argv,s0equ,EMPTY_WORDING);
//     if(aurostd::args2flag(argv,s0neq)) ss=aurostd::args2string(argv,s0neq,EMPTY_WORDING);
//     if(aurostd::args2attachedflag(argv,s1equ)) ss=aurostd::args2attachedstring(argv,s1equ,EMPTY_WORDING);
//     if(aurostd::args2flag(argv,s1neq)) ss=aurostd::args2string(argv,s1neq,EMPTY_WORDING);
//     if(aurostd::args2attachedflag(argv,s2equ)) ss=aurostd::args2attachedstring(argv,s2equ,EMPTY_WORDING);
//     if(aurostd::args2flag(argv,s2neq)) ss=aurostd::args2string(argv,s2neq,EMPTY_WORDING);
//     if(aurostd::args2attachedflag(argv,s3equ)) ss=aurostd::args2attachedstring(argv,s3equ,EMPTY_WORDING);
//     if(aurostd::args2flag(argv,s3neq)) ss=aurostd::args2string(argv,s3neq,EMPTY_WORDING);
//     if(aurostd::substring2bool(ss,"./")) aurostd::StringSubst(ss,"./","");
//     if(ss!="") {
//       if(!aurostd::substring2bool(ss,"/")) { return get_itemized_vector_string_from_input(argv,s0,s1,s2,s3,tokens,delimiter);
//       } else {
// 	aurostd::string2tokens(ss,stokens,"/");
// 	tokens.push_back(stokens.at(1));
// 	KBIN::VASP_SplitAlloyPseudoPotentials(stokens.at(0),stokens);
// 	for(uint i=0;i<stokens.size();i++) tokens.push_back(stokens.at(i));
//       }
//     }
//     if(tokens.size()==1 && aurostd::substring2bool(tokens.at(0),delimiter)) {s0equ=tokens.at(0);aurostd::string2tokens(s0equ,tokens,delimiter);}
//     if(tokens.size()==0) return FALSE;
//     return TRUE;
//   }
// }


// ***************************************************************************
// RICHARD stuff on XRD
namespace pflow {
  double GetAtomicPlaneDist(string options,istream & cin) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "pflow::GetAtomicPlaneDist: BEGIN" << endl;
    vector<string> tokens;
    aurostd::string2tokens(options,tokens,",");
    if(tokens.size()!=3) {
      init::ErrorOption(cout,options,"pflow::GetAtomicPlaneDist","aflow --xrd_dist=h,k,l < POSCAR");
      exit(0);
    } 
    // move on
    double h=0,k=0,l=0;
    if(tokens.size()>=1) h=aurostd::string2utype<double>(tokens.at(0));
    if(tokens.size()>=2) k=aurostd::string2utype<double>(tokens.at(1));
    if(tokens.size()>=3) l=aurostd::string2utype<double>(tokens.at(2));
    

    _aflags aflags;
    aflags.QUIET=TRUE;
    xstructure a(cin,IOVASP_AUTO);
    double tol = 1e-6;  
    double dist;
    //  a = GetStandardConventional(a);
    
    xvector<double> origin;origin(1)=0;origin(2)=0;origin(3)=0;
      
    xvector<double> A(3),B(3),C(3);
  
    //NO ZEROS
    if(h > tol && k > tol && l > tol) {
      A = (1/h)*a.lattice(1);
      B = (1/k)*a.lattice(2);
      C = (1/l)*a.lattice(3);
    }
    //ONE ZERO
    if(h < tol && k > tol && l > tol) { // 0 k l
      B = (1/k)*a.lattice(2);
      A = a.lattice(1)+B;
      C = (1/l)*a.lattice(3);
    }
    if(h > tol && k < tol && l > tol) { // h 0 l
      A = (1/h)*a.lattice(1);
      B =  a.lattice(2)+A;
      C = (1/l)*a.lattice(3);
    }
    if(h > tol && k > tol && l < tol) { // h k 0
      A = (1/h)*a.lattice(1);
      B = (1/k)*a.lattice(2);
      C = a.lattice(3)+B;
    }
    //TWO ZEROS
    if(h < tol && k < tol && l > tol) { // 0 0 l
      C = (1/l)*a.lattice(3);
      A = a.lattice(1)+C;
      B = a.lattice(2)+C;
    }
    if(h < tol && k > tol && l < tol) { // 0 k 0
      B = (1/k)*a.lattice(2);
      A = a.lattice(1)+B;
      C = a.lattice(3)+B;
    }
    if(h > tol && k < tol && l < tol) { // h 0 0
      A = (1/h)*a.lattice(1);
      B = a.lattice(2)+A;
      C = a.lattice(3)+A;
    }
    //GET PLANE NORMAL VECTOR & NORMAL ALONG A
    xvector<double> n =  1/aurostd::modulus(aurostd::vector_product(A-B,A-C))*aurostd::vector_product(A-B,A-C);
    xvector<double> nA = 1/aurostd::modulus(A)*A;
    dist = aurostd::modulus(A)*aurostd::scalar_product(nA,n);
    cout << setprecision(6) << dist << endl;
    if(LDEBUG) cerr << "pflow::GetAtomicPlaneDist: END" << endl;
    return dist;
  }
} // namespace pflow

namespace pflow {
  int whereischar(string str, char c) {
    int out=-1;
    for(uint i=0;i<str.size();i++) {
      if(str[i] == c)
	out=i;
    }
    return out;
  }
} // namespace pflow

namespace pflow {
  bool havechar(string str_in, char c) {
    bool contains=false;
    for(uint i=0;i<str_in.length();i++) {
      if(str_in[i]==c)
	contains=true;
    }
    return contains;
  }
} // namespace pflow
 
namespace pflow {
  void cleanupstring(string & str) {
    uint i = str.size()-1;
    while(str[i]==' ') {
      i--;
    }
    str.erase(i+1,str.size()-1);
    i=0;
    while(str[i]==' ') {
      i++;
    }
    str.erase(0,i);
  }
} // namespace pflow


namespace pflow {
  double frac2dbl(string str) {
    double out;
    bool neg = false;
    cleanupstring(str);
    if(str[0]=='-') {
      neg = true;
      str.erase(str.begin(),str.begin()+1);
    }
    if(str[0]=='+') {
      str.erase(str.begin(),str.begin()+1);
    }
    if(havechar(str,'/')) {
      double numerator;
      double denominator;
      uint slash = whereischar(str, '/');
      char num [256];
      char den [256];
      for(uint i=0;i<slash;i++)
	num[i] = str[i];
      for(uint i=0;i<(str.size()-slash);i++)
	den[i] = str[i+slash+1];
      numerator = atof(num);
      denominator = atof(den);
      if(neg ==true) {
	out = -numerator/denominator;
      }
      else{
	out =numerator/denominator;
      }
    }
    else{
      out = atof(str.c_str());
    }
    return out;
  }
} // namespace pflow

#endif  // 

// **************************************************************************
// *                                                                        *
// *             STEFANO CURTAROLO - Duke University 2003-2018              *
// *                                                                        *
// **************************************************************************
