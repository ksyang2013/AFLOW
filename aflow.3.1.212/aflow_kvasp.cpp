// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *                                                                         *
// ***************************************************************************
// this file contains the routines to run VASP in KBIN
// Stefano Curtarolo - 2007-2012 Duke
/*
   GENERATE, STATIC, KPOINTS, RELAX, RELAX_STATIC, RELAX_STATIC_BANDS, STATIC_BANDS, DIELECTRIC_STATIC, DIELECTRIC_DYNAMIC, DSCF
   */

#ifndef _AFLOW_KVASP_CPP
#define _AFLOW_KVASP_CPP
#define _incarpad_ 26

#include "aflow.h"

#define _KVASP_VASP_SLEEP_   2
#define _KVASP_WAIT_SLEEP_   10
#define _KVASP_CHECK_SLEEP_  60
//#define _KVASP_CHECK_SLEEP_  10
#define KBIN_WRONG_ENTRY_STRING string("WRONG_ENTRY")
#define KBIN_WRONG_ENTRY_NUMBER -123

using aurostd::RemoveWhiteSpaces;
using aurostd::RemoveWhiteSpacesFromTheBack;
using aurostd::FileExist;

pthread_mutex_t mutex_KVASP=PTHREAD_MUTEX_INITIALIZER;

#define _VASP_CONTCAR_SAVE_  TRUE

#define _STROPT_ string("[VASP_FORCE_OPTION]")
#define LDAU_ADIABATIC_RELAX_DEFAULT 6
#define DUKE_MATERIALS_VASP5_CORES_DIELECTRIC 16
#define AFLOWLIB_VASP5_CORES_DIELECTRIC 16
#define DUKE_BETA_VASP5_CORES_DIELECTRIC 16


// ******************************************************************************************************************************************************
// ******************************************************************************************************************************************************

namespace KBIN {
    _vflags VASP_Get_Vflags_from_AflowIN(const string &AflowIn,_aflags &aflags,_kflags &kflags) {
        ofstream FileMESSAGE("/dev/null");
        return KBIN::VASP_Get_Vflags_from_AflowIN(AflowIn,FileMESSAGE,aflags,kflags);
    }
} // namespace KBIN


namespace KBIN {
    _vflags VASP_Get_Vflags_from_AflowIN(const string &_AflowIn,ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags) {
        bool LDEBUG=(FALSE || XHOST.DEBUG);
        _vflags vflags;
        string AflowIn=aurostd::RemoveComments(_AflowIn); // for safety // CO 180502
        vector<string> vAflowIn;aurostd::string2vectorstring(AflowIn,vAflowIn);
        string BflowIn=AflowIn;

        if(LDEBUG) cerr << "DEBUG: KBIN::VASP_Get_Vflags_from_AflowIN (START)" << endl;
        // HOW TO RUN
        vflags.KBIN_VASP_RUN_NRELAX=0;
        // [OBSOLETE]  vflags.KBIN_VASP_RUN_GENERATE           =(aurostd::substring2bool(AflowIn,"[VASP_RUN_GENERATE]") || aurostd::substring2bool(AflowIn,"[VASP_RUN]GENERATE")) || aflags.KBIN_GEN_VASP_FROM_AFLOWIN;
        // [OBSOLETE] vflags.KBIN_VASP_RUN_STATIC              =(aurostd::substring2bool(AflowIn,"[VASP_RUN_STATIC]") || aurostd::substring2bool(AflowIn,"[VASP_RUN]STATIC"));
        // [OBSOLETE] vflags.KBIN_VASP_RUN_KPOINTS             =(aurostd::substring2bool(AflowIn,"[VASP_RUN_KPOINTS]") || aurostd::substring2bool(AflowIn,"[VASP_RUN]KPOINTS"));

        vflags.KBIN_VASP_RUN.clear();
        if((aurostd::substring2bool(AflowIn,"[VASP_RUN_GENERATE]") || aurostd::substring2bool(AflowIn,"[VASP_RUN]GENERATE")) || aflags.KBIN_GEN_VASP_FROM_AFLOWIN) 
            vflags.KBIN_VASP_RUN.push("GENERATE");
        if((aurostd::substring2bool(AflowIn,"[VASP_RUN_STATIC]") || aurostd::substring2bool(AflowIn,"[VASP_RUN]STATIC")))
            vflags.KBIN_VASP_RUN.push("STATIC");
        if((aurostd::substring2bool(AflowIn,"[VASP_RUN_KPOINTS]") || aurostd::substring2bool(AflowIn,"[VASP_RUN]KPOINTS"))) 
            vflags.KBIN_VASP_RUN.push("KPOINTS");

        for(uint i=0;i<vAflowIn.size();i++) {
            if(aurostd::substring2bool(vAflowIn.at(i),"VASP_RUN")) {
                string vasp_run_string=vAflowIn.at(i);
                if(vasp_run_string.find("#")!=string::npos) vasp_run_string=vasp_run_string.substr(0,vasp_run_string.find("#"));
                if(vasp_run_string.find("//")!=string::npos) vasp_run_string=vasp_run_string.substr(0,vasp_run_string.find("//"));
                if(vasp_run_string.find("!")!=string::npos) vasp_run_string=vasp_run_string.substr(0,vasp_run_string.find("!"));

                //      cout << vasp_run_string << endl;
                vector<string> aflowin_tokens;
                aurostd::string2tokens(vasp_run_string,aflowin_tokens,",");
                //      cerr << aflowin_tokens.size() << endl;exit(0);
                if(aflowin_tokens.size()>0) {
                    // [OBSOLETE] vflags.KBIN_VASP_RUN_RELAX               =(aurostd::substring2bool(aflowin_tokens.at(0),"[VASP_RUN_RELAX=") || aurostd::substring2bool(aflowin_tokens.at(0),"[VASP_RUN]RELAX="));
                    if((aurostd::substring2bool(aflowin_tokens.at(0),"[VASP_RUN_RELAX=") || aurostd::substring2bool(aflowin_tokens.at(0),"[VASP_RUN]RELAX="))) vflags.KBIN_VASP_RUN.push("RELAX");
                    if(aurostd::substring2bool(aflowin_tokens.at(0),"[VASP_RUN_RELAX=")) vflags.KBIN_VASP_RUN_NRELAX=aurostd::substring2utype<int>(aflowin_tokens.at(0),"[VASP_RUN_RELAX=");
                    if(aurostd::substring2bool(aflowin_tokens.at(0),"[VASP_RUN]RELAX=")) vflags.KBIN_VASP_RUN_NRELAX=aurostd::substring2utype<int>(aflowin_tokens.at(0),"[VASP_RUN]RELAX=");
                    // [OBSOLETE] vflags.KBIN_VASP_RUN_RELAX_STATIC        =(aurostd::substring2bool(aflowin_tokens.at(0),"[VASP_RUN_RELAX_STATIC=") || aurostd::substring2bool(aflowin_tokens.at(0),"[VASP_RUN]RELAX_STATIC="));
                    if((aurostd::substring2bool(aflowin_tokens.at(0),"[VASP_RUN_RELAX_STATIC=") || aurostd::substring2bool(aflowin_tokens.at(0),"[VASP_RUN]RELAX_STATIC="))) vflags.KBIN_VASP_RUN.push("RELAX_STATIC");
                    if(aurostd::substring2bool(aflowin_tokens.at(0),"[VASP_RUN_RELAX_STATIC=")) vflags.KBIN_VASP_RUN_NRELAX=aurostd::substring2utype<int>(aflowin_tokens.at(0),"[VASP_RUN_RELAX_STATIC=");
                    if(aurostd::substring2bool(aflowin_tokens.at(0),"[VASP_RUN]RELAX_STATIC=")) vflags.KBIN_VASP_RUN_NRELAX=aurostd::substring2utype<int>(aflowin_tokens.at(0),"[VASP_RUN]RELAX_STATIC=");
                    // [OBSOLETE] vflags.KBIN_VASP_RUN_RELAX_STATIC_BANDS  =(aurostd::substring2bool(aflowin_tokens.at(0),"[VASP_RUN_RELAX_STATIC_BANDS=") || aurostd::substring2bool(aflowin_tokens.at(0),"[VASP_RUN]RELAX_STATIC_BANDS="));
                    if((aurostd::substring2bool(aflowin_tokens.at(0),"[VASP_RUN_RELAX_STATIC_BANDS=") || aurostd::substring2bool(aflowin_tokens.at(0),"[VASP_RUN]RELAX_STATIC_BANDS="))) vflags.KBIN_VASP_RUN.push("RELAX_STATIC_BANDS");
                    if(aurostd::substring2bool(aflowin_tokens.at(0),"[VASP_RUN_RELAX_STATIC_BANDS=")) vflags.KBIN_VASP_RUN_NRELAX=aurostd::substring2utype<int>(aflowin_tokens.at(0),"[VASP_RUN_RELAX_STATIC_BANDS=");
                    if(aurostd::substring2bool(aflowin_tokens.at(0),"[VASP_RUN]RELAX_STATIC_BANDS=")) vflags.KBIN_VASP_RUN_NRELAX=aurostd::substring2utype<int>(aflowin_tokens.at(0),"[VASP_RUN]RELAX_STATIC_BANDS=");
                    // [OBSOLETE] vflags.KBIN_VASP_RUN_STATIC_BANDS        =(aurostd::substring2bool(aflowin_tokens.at(0),"[VASP_RUN_STATIC_BANDS]") || aurostd::substring2bool(aflowin_tokens.at(0),"[VASP_RUN]STATIC_BANDS"));
                    if((aurostd::substring2bool(aflowin_tokens.at(0),"[VASP_RUN_STATIC_BANDS]") || aurostd::substring2bool(aflowin_tokens.at(0),"[VASP_RUN]STATIC_BANDS"))) vflags.KBIN_VASP_RUN.push("STATIC_BANDS");
                }
                for(uint j=1;j<aflowin_tokens.size();j++) {
                    // [OBSOLETE] vflags.KBIN_VASP_RUN_DIELECTRIC_STATIC   =(vflags.KBIN_VASP_RUN_DIELECTRIC_STATIC || aurostd::substring2bool(aflowin_tokens.at(j),"DS") || aurostd::substring2bool(aflowin_tokens.at(j),"DIELECTRIC_STATIC"));
                    if((vflags.KBIN_VASP_RUN.flag("DIELECTRIC_STATIC") || aurostd::substring2bool(aflowin_tokens.at(j),"DS") || aurostd::substring2bool(aflowin_tokens.at(j),"DIELECTRIC_STATIC"))) vflags.KBIN_VASP_RUN.push("DIELECTRIC_STATIC");
                    // [OBSOLETE] vflags.KBIN_VASP_RUN_DIELECTRIC_DYNAMIC  =(vflags.KBIN_VASP_RUN_DIELECTRIC_DYNAMIC || aurostd::substring2bool(aflowin_tokens.at(j),"DD") || aurostd::substring2bool(aflowin_tokens.at(j),"DIELECTRIC_DYNAMIC"));
                    if((vflags.KBIN_VASP_RUN.flag("DIELECTRIC_DYNAMIC") || aurostd::substring2bool(aflowin_tokens.at(j),"DD") || aurostd::substring2bool(aflowin_tokens.at(j),"DIELECTRIC_DYNAMIC"))) vflags.KBIN_VASP_RUN.push("DIELECTRIC_DYNAMIC");
                    // [OBSOLETE] vflags.KBIN_VASP_RUN_DSCF                =(vflags.KBIN_VASP_RUN_DSCF || aurostd::substring2bool(aflowin_tokens.at(j),"DSCF"));
                    if((vflags.KBIN_VASP_RUN.flag("DSCF") || aurostd::substring2bool(aflowin_tokens.at(j),"DSCF"))) vflags.KBIN_VASP_RUN.push("DSCF");
                }
                if(vflags.KBIN_VASP_RUN.flag("DSCF")) vflags.KBIN_VASP_RUN.push("DIELECTRIC_DYNAMIC");
                if(vflags.KBIN_VASP_RUN.flag("DIELECTRIC_DYNAMIC")) vflags.KBIN_VASP_RUN.push("DIELECTRIC_STATIC");
            }
        }
        if(vflags.KBIN_VASP_RUN.xscheme!="") vflags.KBIN_VASP_RUN.isentry=TRUE;

        // if(vflags.KBIN_VASP_RUN.flag("DIELECTRIC_STATIC")) cout << "vflags.KBIN_VASP_RUN.flag(\"DIELECTRIC_STATIC\")" << endl;
        // if(vflags.KBIN_VASP_RUN.flag("DIELECTRIC_DYNAMIC")) cout << "vflags.KBIN_VASP_RUN.flag(\"DIELECTRIC_DYNAMIC\")" << endl;
        // if(vflags.KBIN_VASP_RUN.flag("DSCF")) cout << "vflags.KBIN_VASP_RUN.flag(\"DSCF\")" << endl;
        // exit(0);

        // [OBSOLETE] vflags.KBIN_VASP_REPEAT_BANDS        = aurostd::FileExist(aflags.Directory+string("/REPEAT_BANDS")) || aurostd::substring2bool(AflowIn,"[VASP_RUN_REPEAT_BANDS]") || aurostd::substring2bool(AflowIn,"[VASP_RUN]REPEAT_BANDS]");
        // [OBSOLETE] vflags.KBIN_VASP_REPEAT_STATIC_BANDS = aurostd::FileExist(aflags.Directory+string("/REPEAT_STATIC_BANDS")) || aurostd::substring2bool(AflowIn,"[VASP_RUN_REPEAT_STATIC_BANDS]") || aurostd::substring2bool(AflowIn,"[VASP_RUN]REPEAT_STATIC_BANDS]");
        // [OBSOLETE] vflags.KBIN_VASP_REPEAT_DELSOL       = aurostd::FileExist(aflags.Directory+string("/REPEAT_DELSOL")) || aurostd::substring2bool(AflowIn,"[VASP_RUN_REPEAT_DELSOL]") || aurostd::substring2bool(AflowIn,"[VASP_RUN]REPEAT_DELSOL]");

        vflags.KBIN_VASP_REPEAT.clear();
        if(aurostd::FileExist(aflags.Directory+string("/REPEAT_BANDS")) || aurostd::substring2bool(AflowIn,"[VASP_RUN_REPEAT_BANDS]") || aurostd::substring2bool(AflowIn,"[VASP_RUN]REPEAT_BANDS]")) vflags.KBIN_VASP_REPEAT.push("REPEAT_BANDS");
        if(aurostd::FileExist(aflags.Directory+string("/REPEAT_STATIC_BANDS")) || aurostd::substring2bool(AflowIn,"[VASP_RUN_REPEAT_STATIC_BANDS]") || aurostd::substring2bool(AflowIn,"[VASP_RUN]REPEAT_STATIC_BANDS]")) vflags.KBIN_VASP_REPEAT.push("REPEAT_STATIC_BANDS");
        if(aurostd::FileExist(aflags.Directory+string("/REPEAT_DELSOL")) || aurostd::substring2bool(AflowIn,"[VASP_RUN_REPEAT_DELSOL]") || aurostd::substring2bool(AflowIn,"[VASP_RUN]REPEAT_DELSOL]")) vflags.KBIN_VASP_REPEAT.push("REPEAT_DELSOL");


        if(vflags.KBIN_VASP_REPEAT.flag("REPEAT_STATIC_BANDS")) cout << "vflags.KBIN_VASP_REPEAT.flag(\"REPEAT_STATIC_BANDS\")" << endl;

        // priorities about RUN
        if(vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_BANDS")) {  // RELAX_STATIC_BANDS
            //  cerr << "[DEBUG] vflags.KBIN_VASP_RUN.flag(\"RELAX_STATIC_BANDS\")==TRUE" << endl;
            vflags.KBIN_VASP_RUN.flag("STATIC",FALSE);
            vflags.KBIN_VASP_RUN.flag("STATIC_BANDS",FALSE);
            vflags.KBIN_VASP_RUN.flag("KPOINTS",FALSE);
            vflags.KBIN_VASP_RUN.flag("RELAX",FALSE);
            vflags.KBIN_VASP_RUN.flag("RELAX_STATIC",FALSE);
            vflags.KBIN_VASP_RUN.push("RELAX_STATIC_BANDS");
        } else {
            if(vflags.KBIN_VASP_RUN.flag("RELAX_STATIC")) {  // RELAX_STATIC
                //  cerr << "[DEBUG] vflags.KBIN_VASP_RUN.flag(\"RELAX_STATIC\")==TRUE" << endl;
                vflags.KBIN_VASP_RUN.flag("STATIC",FALSE);
                vflags.KBIN_VASP_RUN.flag("STATIC_BANDS",FALSE);
                vflags.KBIN_VASP_RUN.flag("KPOINTS",FALSE);
                vflags.KBIN_VASP_RUN.flag("RELAX",FALSE);
                vflags.KBIN_VASP_RUN.push("RELAX_STATIC");
                vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_BANDS",FALSE);
            } else {
                if(vflags.KBIN_VASP_RUN.flag("STATIC_BANDS")) {  // STATIC_BANDS
                    //  cerr << "[DEBUG] vflags.KBIN_VASP_RUN.flag(\"STATIC_BANDS\")==TRUE" << endl;
                    vflags.KBIN_VASP_RUN.flag("STATIC",FALSE);
                    vflags.KBIN_VASP_RUN.push("STATIC_BANDS");
                    vflags.KBIN_VASP_RUN.flag("KPOINTS",FALSE);
                    vflags.KBIN_VASP_RUN.flag("RELAX",FALSE);
                    vflags.KBIN_VASP_RUN.flag("RELAX_STATIC",FALSE);
                    vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_BANDS",FALSE);
                } else {                                  
                    if(vflags.KBIN_VASP_RUN.flag("STATIC")) {  // STATIC
                        //	  cerr << "[DEBUG] vflags.KBIN_VASP_RUN.flag(\"STATIC\")==TRUE" << endl;
                        vflags.KBIN_VASP_RUN.flag("KPOINTS",FALSE);
                        vflags.KBIN_VASP_RUN.flag("RELAX",FALSE);
                        vflags.KBIN_VASP_RUN.push("STATIC");
                        vflags.KBIN_VASP_RUN.flag("STATIC_BANDS",FALSE);
                        vflags.KBIN_VASP_RUN.flag("RELAX_STATIC",FALSE);
                        vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_BANDS",FALSE);
                    } else {                            
                        if(vflags.KBIN_VASP_RUN.flag("KPOINTS")) {  // KPOINTS
                            //  cerr << "[DEBUG] vflags.KBIN_VASP_RUN.flag(\"KPOINTS\")==TRUE" << endl;
                            vflags.KBIN_VASP_RUN.flag("STATIC",FALSE);
                            vflags.KBIN_VASP_RUN.flag("STATIC_BANDS",FALSE);
                            vflags.KBIN_VASP_RUN.flag("RELAX",FALSE);
                            vflags.KBIN_VASP_RUN.push("KPOINTS");
                            vflags.KBIN_VASP_RUN.flag("RELAX_STATIC",FALSE);
                            vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_BANDS",FALSE);
                        } else {
                            //   cerr << "[DEBUG] DEFAULT" << endl;
                            vflags.KBIN_VASP_RUN.push("RELAX");
                            vflags.KBIN_VASP_RUN.flag("STATIC",FALSE);
                            vflags.KBIN_VASP_RUN.flag("STATIC_BANDS",FALSE);
                            vflags.KBIN_VASP_RUN.flag("KPOINTS",FALSE);
                            vflags.KBIN_VASP_RUN.flag("RELAX_STATIC",FALSE);
                            vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_BANDS",FALSE);
                        }
                    }
                }
            }
        }
        // priorities about REPEAT
        if(vflags.KBIN_VASP_REPEAT.flag("REPEAT_STATIC_BANDS")) {
            vflags.KBIN_VASP_RUN.flag("STATIC",FALSE);
            vflags.KBIN_VASP_RUN.flag("STATIC_BANDS",FALSE);
            vflags.KBIN_VASP_RUN.flag("KPOINTS",FALSE);
            vflags.KBIN_VASP_RUN.flag("RELAX",FALSE);
            vflags.KBIN_VASP_RUN.flag("RELAX_STATIC",FALSE);
            vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_BANDS",FALSE);
            vflags.KBIN_VASP_REPEAT.flag("REPEAT_BANDS",FALSE);
        }
        if(vflags.KBIN_VASP_REPEAT.flag("REPEAT_BANDS")) {
            vflags.KBIN_VASP_RUN.flag("STATIC",FALSE);
            vflags.KBIN_VASP_RUN.flag("STATIC_BANDS",FALSE);
            vflags.KBIN_VASP_RUN.flag("KPOINTS",FALSE);
            vflags.KBIN_VASP_RUN.flag("RELAX",FALSE);
            vflags.KBIN_VASP_RUN.flag("RELAX_STATIC",FALSE);
            vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_BANDS",FALSE);
        }
        if(vflags.KBIN_VASP_REPEAT.flag("REPEAT_DELSOL")) {
            vflags.KBIN_VASP_RUN.flag("STATIC",FALSE);
            vflags.KBIN_VASP_RUN.flag("STATIC_BANDS",FALSE);
            vflags.KBIN_VASP_RUN.flag("KPOINTS",FALSE);
            vflags.KBIN_VASP_RUN.flag("RELAX",FALSE);
            vflags.KBIN_VASP_RUN.flag("RELAX_STATIC",FALSE);
            vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_BANDS",FALSE);
            vflags.KBIN_VASP_REPEAT.flag("REPEAT_BANDS",FALSE);
        }

        if(kflags.KBIN_PHONONS_CALCULATION_APL || kflags.KBIN_PHONONS_CALCULATION_QHA || kflags.KBIN_PHONONS_CALCULATION_AAPL || kflags.KBIN_PHONONS_CALCULATION_FROZSL || kflags.KBIN_PHONONS_CALCULATION_AGL || kflags.KBIN_PHONONS_CALCULATION_AEL) {  // CO 170601
            vflags.KBIN_VASP_REPEAT.flag("REPEAT_DELSOL",FALSE);
            vflags.KBIN_VASP_RUN.flag("STATIC",FALSE);
            vflags.KBIN_VASP_RUN.flag("STATIC_BANDS",FALSE);
            vflags.KBIN_VASP_RUN.flag("KPOINTS",FALSE);
            vflags.KBIN_VASP_RUN.flag("RELAX",FALSE);
            vflags.KBIN_VASP_RUN.flag("RELAX_STATIC",FALSE);
            vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_BANDS",FALSE);
            vflags.KBIN_VASP_REPEAT.flag("REPEAT_BANDS",FALSE); 
            kflags.KBIN_SYMMETRY_CALCULATION=FALSE;
        }

        // RELAX_MODE AND PRIORITIES  // ENERGY | FORCES | ENERGY_FORCES | FORCES_ENERGY (default ENERGY) "
        vflags.KBIN_VASP_FORCE_OPTION_RELAX_MODE.options2entry(AflowIn,_STROPT_+"RELAX_MODE=",FALSE,DEFAULT_VASP_FORCE_OPTION_RELAX_MODE_SCHEME);
        vflags.KBIN_VASP_FORCE_OPTION_RELAX_MODE.xscheme=KBIN_WRONG_ENTRY_STRING;
        vflags.KBIN_VASP_FORCE_OPTION_RELAX_MODE.scheme2scheme("ENERGY","ENERGY");
        vflags.KBIN_VASP_FORCE_OPTION_RELAX_MODE.scheme2scheme("FORCES","FORCES");vflags.KBIN_VASP_FORCE_OPTION_RELAX_MODE.scheme2scheme("FORCE","FORCES");
        vflags.KBIN_VASP_FORCE_OPTION_RELAX_MODE.scheme2scheme("ENERGY_FORCES","ENERGY_FORCES");vflags.KBIN_VASP_FORCE_OPTION_RELAX_MODE.scheme2scheme("ENERGY_FORCE","ENERGY_FORCES");
        vflags.KBIN_VASP_FORCE_OPTION_RELAX_MODE.scheme2scheme("FORCES_ENERGY","FORCES_ENERGY");vflags.KBIN_VASP_FORCE_OPTION_RELAX_MODE.scheme2scheme("FORCE_ENERGY","FORCES_ENERGY");
        if(vflags.KBIN_VASP_FORCE_OPTION_RELAX_MODE.isentry && vflags.KBIN_VASP_FORCE_OPTION_RELAX_MODE.xscheme==KBIN_WRONG_ENTRY_STRING) {
            cerr << "ERROR: KBIN::VASP_Get_Vflags_from_AflowIN:  vflags.KBIN_VASP_FORCE_OPTION_RELAX_MODE.content_string="  <<  vflags.KBIN_VASP_FORCE_OPTION_RELAX_MODE.content_string << endl;
            exit(0);
        }
        //  cerr << "vflags.KBIN_VASP_FORCE_OPTION_RELAX_MODE.xscheme=" << vflags.KBIN_VASP_FORCE_OPTION_RELAX_MODE.xscheme << endl; exit(0);

        // FORCE OPTIONS
        vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.options2entry(AflowIn,_STROPT_+"NOTUNE");

        // FORCE OPTIONS SYSTEM_AUTO
        vflags.KBIN_VASP_FORCE_OPTION_SYSTEM_AUTO.options2entry(AflowIn,_STROPT_+"SYSTEM_AUTO");

        // FORCE OPTIONS STATIC RELAX_ALL RELAX_IONS RELAX CELL_VOLUME 
        // cerr << "vflags.KBIN_VASP_FORCE_OPTION_STATIC= " << vflags.KBIN_VASP_FORCE_OPTION_STATIC << endl;
        vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.clear();
        if(aurostd::substring2bool(AflowIn,_STROPT_+"STATIC",TRUE)) vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.push("STATIC");
        if(aurostd::substring2bool(AflowIn,_STROPT_+"RELAX_ALL",TRUE) || aurostd::substring2bool(AflowIn,_STROPT_+"RELAX",TRUE))  vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.push("ALL");
        if(aurostd::substring2bool(AflowIn,_STROPT_+"RELAX_IONS",TRUE)) vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.push("IONS");
        if(aurostd::substring2bool(AflowIn,_STROPT_+"RELAX_CELL_SHAPE",TRUE) || aurostd::substring2bool(AflowIn,_STROPT_+"RELAX_SHAPE",TRUE)) vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.push("CELL_SHAPE");
        if(aurostd::substring2bool(AflowIn,_STROPT_+"RELAX_CELL_VOLUME",TRUE) || aurostd::substring2bool(AflowIn,_STROPT_+"RELAX_VOLUME",TRUE)) vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.push("CELL_VOLUME");
        if(aurostd::substring2bool(AflowIn,_STROPT_+"RELAX_IONS_CELL_VOLUME",TRUE) || aurostd::substring2bool(AflowIn,_STROPT_+"RELAX_IONS_VOLUME",TRUE)) vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.push("IONS_CELL_VOLUME");
        if(vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.xscheme!="") vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.isentry=TRUE;

        // [OBSOLETE] vflags.KBIN_VASP_FORCE_OPTION_STATIC=aurostd::substring2bool(AflowIn,_STROPT_+"STATIC",TRUE);
        // [OBSOLETE] cerr<< "vflags.KBIN_VASP_FORCE_OPTION_STATIC= " << vflags.KBIN_VASP_FORCE_OPTION_STATIC << endl;
        // [OBSOLETE] vflags.KBIN_VASP_FORCE_OPTION_RELAX_ALL = aurostd::substring2bool(AflowIn,_STROPT_+"RELAX",TRUE) || aurostd::substring2bool(AflowIn,_STROPT_+"RELAX_ALL",TRUE);
        // [OBSOLETE] vflags.KBIN_VASP_FORCE_OPTION_RELAX_IONS = aurostd::substring2bool(AflowIn,_STROPT_+"RELAX_IONS",TRUE);
        // [OBSOLETE] vvflags.KBIN_VASP_FORCE_OPTION_RELAX_CELL_SHAPE = aurostd::substring2bool(AflowIn,_STROPT_+"RELAX_CELL_SHAPE",TRUE) || aurostd::substring2bool(AflowIn,_STROPT_+"RELAX_SHAPE",TRUE);
        // [OBSOLETE] vvflags.KBIN_VASP_FORCE_OPTION_RELAX_CELL_VOLUME = aurostd::substring2bool(AflowIn,_STROPT_+"RELAX_CELL_VOLUME",TRUE) || aurostd::substring2bool(AflowIn,_STROPT_+"RELAX_VOLUME",TRUE);
        // [OBSOLETE] vvflags.KBIN_VASP_FORCE_OPTION_RELAX_IONS_CELL_VOLUME = aurostd::substring2bool(AflowIn,_STROPT_+"RELAX_IONS_CELL_VOLUME",TRUE) || aurostd::substring2bool(AflowIn,_STROPT_+"RELAX_IONS_VOLUME",TRUE);
        // [OBSOLETE] if(vflags.KBIN_VASP_FORCE_OPTION_RELAX_IONS_CELL_VOLUME) vflags.KBIN_VASP_FORCE_OPTION_RELAX_IONS=FALSE;
        // [OBSOLETE] if(vflags.KBIN_VASP_FORCE_OPTION_RELAX_ALL && (vflags.KBIN_VASP_FORCE_OPTION_RELAX_IONS || vflags.KBIN_VASP_FORCE_OPTION_RELAX_CELL_SHAPE ||
        // [OBSOLETE]  vflags.KBIN_VASP_FORCE_OPTION_RELAX_CELL_VOLUME || vflags.KBIN_VASP_FORCE_OPTION_RELAX_IONS_CELL_VOLUME)) vflags.KBIN_VASP_FORCE_OPTION_RELAX_ALL=FALSE;

        if(vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("IONS_CELL_VOLUME"))
            vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("IONS",FALSE);
        if(vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("ALL") && (vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("STATIC") || vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("IONS") || 
                    vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("CELL_SHAPE") || vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("CELL_VOLUME") ||
                    vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("IONS_CELL_VOLUME")))    
            vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("ALL",FALSE);
        if(vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("STATIC")) {
            vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.clear();
            vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.push("STATIC");
            vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.isentry=TRUE;}

        // cerr << "KBIN::VASP_Get_Vflags_from_AflowIN:  vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.xscheme="  <<  vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.xscheme << endl; exit(0);

        // PRECISION AND PRIORITIES // (LOW | MEDIUM | NORMAL | HIGH | ACCURATE), PRESERVED
        vflags.KBIN_VASP_FORCE_OPTION_PREC.options2entry(AflowIn,_STROPT_+"PREC=",FALSE,DEFAULT_VASP_FORCE_OPTION_PREC_SCHEME);
        vflags.KBIN_VASP_FORCE_OPTION_PREC.xscheme=KBIN_WRONG_ENTRY_STRING;
        vflags.KBIN_VASP_FORCE_OPTION_PREC.scheme2scheme('L',"LOW");
        vflags.KBIN_VASP_FORCE_OPTION_PREC.scheme2scheme('M',"MEDIUM");
        vflags.KBIN_VASP_FORCE_OPTION_PREC.scheme2scheme('N',"NORMAL");
        vflags.KBIN_VASP_FORCE_OPTION_PREC.scheme2scheme('H',"HIGH");
        vflags.KBIN_VASP_FORCE_OPTION_PREC.scheme2scheme('A',"ACCURATE");
        vflags.KBIN_VASP_FORCE_OPTION_PREC.scheme2scheme('P',"PHONONS"); // JJPR Modification
        if(vflags.KBIN_VASP_FORCE_OPTION_PREC.isentry && vflags.KBIN_VASP_FORCE_OPTION_PREC.xscheme==KBIN_WRONG_ENTRY_STRING) {
            cerr << "ERROR: KBIN::VASP_Get_Vflags_from_AflowIN:  vflags.KBIN_VASP_FORCE_OPTION_PREC.content_string="  <<  vflags.KBIN_VASP_FORCE_OPTION_PREC.content_string << endl;
            exit(0);
        }

        // ALGO AND PRIORITIES // (NORMAL | VERYFAST | FAST | ALL | DAMPED), PRESERVED
        vflags.KBIN_VASP_FORCE_OPTION_ALGO.options2entry(AflowIn,_STROPT_+"ALGO=",FALSE,DEFAULT_VASP_FORCE_OPTION_ALGO_SCHEME);
        vflags.KBIN_VASP_FORCE_OPTION_ALGO.xscheme=KBIN_WRONG_ENTRY_STRING;
        vflags.KBIN_VASP_FORCE_OPTION_ALGO.scheme2scheme('N',"NORMAL");
        vflags.KBIN_VASP_FORCE_OPTION_ALGO.scheme2scheme('V',"VERYFAST");
        vflags.KBIN_VASP_FORCE_OPTION_ALGO.scheme2scheme('F',"FAST");
        vflags.KBIN_VASP_FORCE_OPTION_ALGO.scheme2scheme('A',"ALL");
        vflags.KBIN_VASP_FORCE_OPTION_ALGO.scheme2scheme('D',"DAMPED");
        if(vflags.KBIN_VASP_FORCE_OPTION_ALGO.isentry && vflags.KBIN_VASP_FORCE_OPTION_ALGO.xscheme==KBIN_WRONG_ENTRY_STRING) {
            cerr << "ERROR: KBIN::VASP_Get_Vflags_from_AflowIN:  vflags.KBIN_VASP_FORCE_OPTION_ALGO.content_string="  <<  vflags.KBIN_VASP_FORCE_OPTION_ALGO.content_string << endl;
            exit(0);
        }
        vflags.KBIN_VASP_FORCE_OPTION_ALGO.preserved= vflags.KBIN_VASP_FORCE_OPTION_ALGO.preserved || aurostd::substring2bool(AflowIn,_STROPT_+"ALGO_PRESERVED",TRUE); // FIX ALGO_PRESERVED

        // ABMIX AND PRIORITIES // empty | [AUTO | US | PAW | #AMIX,#BMIX[,#AMIX_MAG,#BMIX_MAG]]
        vflags.KBIN_VASP_FORCE_OPTION_ABMIX.options2entry(AflowIn,_STROPT_+"ABMIX=",FALSE,DEFAULT_VASP_FORCE_OPTION_ABMIX_SCHEME);
        vflags.KBIN_VASP_FORCE_OPTION_ABMIX.xscheme=KBIN_WRONG_ENTRY_STRING;
        vflags.KBIN_VASP_FORCE_OPTION_ABMIX.scheme2scheme('A',"AUTO");
        vflags.KBIN_VASP_FORCE_OPTION_ABMIX.scheme2scheme('U',"US");
        vflags.KBIN_VASP_FORCE_OPTION_ABMIX.scheme2scheme('L',"US"); // LDA
        vflags.KBIN_VASP_FORCE_OPTION_ABMIX.scheme2scheme('G',"US"); // GGA
        vflags.KBIN_VASP_FORCE_OPTION_ABMIX.scheme2scheme('P',"PAW");  // something with PAW..
        if(vflags.KBIN_VASP_FORCE_OPTION_ABMIX.isentry && vflags.KBIN_VASP_FORCE_OPTION_ABMIX.xscheme==KBIN_WRONG_ENTRY_STRING) {
            cerr << "ERROR: KBIN::VASP_Get_Vflags_from_AflowIN:  vflags.KBIN_VASP_FORCE_OPTION_ABMIX.content_string="  <<  vflags.KBIN_VASP_FORCE_OPTION_ABMIX.content_string << endl;
            exit(0);
        }

        // METAGGA AND PRIORITIES // TPSS | RTPSS | M06L | MBJL | SCAN | MS0 | MS1 | MS2 | NONE
        if(LDEBUG) cerr << "KBIN::VASP_Get_Vflags_from_AflowIN: METAGGA" << endl;
        vflags.KBIN_VASP_FORCE_OPTION_METAGGA.options2entry(AflowIn,_STROPT_+"METAGGA=",FALSE,DEFAULT_VASP_FORCE_OPTION_METAGGA_SCHEME);
        // vflags.KBIN_VASP_FORCE_OPTION_METAGGA.xscheme=KBIN_WRONG_ENTRY_STRING;
        // vflags.KBIN_VASP_FORCE_OPTION_METAGGA.scheme2scheme('T',"TPSS");
        // vflags.KBIN_VASP_FORCE_OPTION_METAGGA.scheme2scheme('R',"RTPSS");
        // vflags.KBIN_VASP_FORCE_OPTION_METAGGA.scheme2scheme('S',"SCAN");
        // vflags.KBIN_VASP_FORCE_OPTION_METAGGA.scheme2scheme('N',"NONE");
        if(vflags.KBIN_VASP_FORCE_OPTION_METAGGA.isentry && vflags.KBIN_VASP_FORCE_OPTION_METAGGA.xscheme==KBIN_WRONG_ENTRY_STRING) {
            cerr << "ERROR: KBIN::VASP_Get_Vflags_from_AflowIN:  vflags.KBIN_VASP_FORCE_OPTION_METAGGA.content_string="  <<  vflags.KBIN_VASP_FORCE_OPTION_METAGGA.content_string << endl;
            exit(0);
        } 
        if(LDEBUG) cerr << "KBIN::VASP_Get_Vflags_from_AflowIN: METAGGA vflags.KBIN_VASP_FORCE_OPTION_METAGGA.isentry=" << vflags.KBIN_VASP_FORCE_OPTION_METAGGA.isentry << endl;
        if(LDEBUG) cerr << "KBIN::VASP_Get_Vflags_from_AflowIN: METAGGA vflags.KBIN_VASP_FORCE_OPTION_METAGGA.xscheme=" << vflags.KBIN_VASP_FORCE_OPTION_METAGGA.xscheme << endl;

        // IVDW AND PRIORITIES // [number_for_VASP_see_manual_for_IVDW | 0] 
        if(LDEBUG) cerr << "KBIN::VASP_Get_Vflags_from_AflowIN: IVDW" << endl;
        vflags.KBIN_VASP_FORCE_OPTION_IVDW.options2entry(AflowIn,_STROPT_+"IVDW=",FALSE,DEFAULT_VASP_FORCE_OPTION_IVDW_SCHEME);
        if(vflags.KBIN_VASP_FORCE_OPTION_IVDW.isentry && vflags.KBIN_VASP_FORCE_OPTION_IVDW.xscheme==KBIN_WRONG_ENTRY_STRING) {
            cerr << "ERROR: KBIN::VASP_Get_Vflags_from_AflowIN:  vflags.KBIN_VASP_FORCE_OPTION_IVDW.content_string="  <<  vflags.KBIN_VASP_FORCE_OPTION_IVDW.content_string << endl;
            exit(0);
        } 
        if(LDEBUG) cerr << "KBIN::VASP_Get_Vflags_from_AflowIN: IVDW vflags.KBIN_VASP_FORCE_OPTION_IVDW.isentry=" << vflags.KBIN_VASP_FORCE_OPTION_IVDW.isentry << endl;
        if(LDEBUG) cerr << "KBIN::VASP_Get_Vflags_from_AflowIN: IVDW vflags.KBIN_VASP_FORCE_OPTION_IVDW.xscheme=" << vflags.KBIN_VASP_FORCE_OPTION_IVDW.xscheme << endl;

        // NEGLECT_NOMIX
        vflags.KBIN_VASP_FORCE_OPTION_SKIP_NOMIX.options2entry(AflowIn,string(_STROPT_+"NEGLECT_IMMISCIBLE"+"|"+_STROPT_+"NEGLECT_NOMIX"+"|"+_STROPT_+"SKIP_NOMIX"));

        // AUTO_PSEUDOPOTENTIALS and AUTO_PSEUDOPOTENTIALS_TYPE
        vflags.KBIN_VASP_FORCE_OPTION_AUTO_PSEUDOPOTENTIALS.options2entry(AflowIn,_STROPT_+"AUTO_PSEUDOPOTENTIALS=",FALSE,DEFAULT_VASP_PSEUDOPOTENTIAL_TYPE);

        // POTIM
        // cerr << "POTIM" << endl;
        vflags.KBIN_VASP_FORCE_OPTION_POTIM_EQUAL.options2entry(AflowIn,_STROPT_+"POTIM=",TRUE,vflags.KBIN_VASP_FORCE_OPTION_POTIM_EQUAL.xscheme); // scheme already loaded in aflow_xclasses.cpp is "DEFAULT_VASP_PREC_POTIM"

        // PSTRESS
        // cerr << "PSTRESS" << endl;
        vflags.KBIN_VASP_FORCE_OPTION_PSTRESS_EQUAL.options2entry(AflowIn,_STROPT_+"PSTRESS=",TRUE,vflags.KBIN_VASP_FORCE_OPTION_PSTRESS_EQUAL.xscheme); // scheme already loaded in aflow_xclasses.cpp is "0.0"  

        // EDIFFG
        // cerr << "EDIFFG" << endl;
        vflags.KBIN_VASP_FORCE_OPTION_EDIFFG_EQUAL.options2entry(AflowIn,_STROPT_+"EDIFFG=",TRUE,vflags.KBIN_VASP_FORCE_OPTION_EDIFFG_EQUAL.xscheme); // scheme already loaded in aflow_xclasses.cpp is "DEFAULT_VASP_PREC_EDIFFG"  

        // NBANDS and/or NBANDS=
        //  cerr << "NBANDS_AUTO" << endl;
        vflags.KBIN_VASP_FORCE_OPTION_NBANDS_AUTO_isentry = aurostd::substring2bool(AflowIn,_STROPT_+"NBANDS",TRUE);
        // cerr << "vflags.KBIN_VASP_FORCE_OPTION_NBANDS_AUTO_isentry=" << vflags.KBIN_VASP_FORCE_OPTION_NBANDS_AUTO_isentry << endl;
        // cerr << "NBANDS_EQUAL" << endl;
        // cerr << "vflags.KBIN_VASP_FORCE_OPTION_NBANDS_EQUAL.xscheme=" << vflags.KBIN_VASP_FORCE_OPTION_NBANDS_EQUAL.xscheme << endl;
        vflags.KBIN_VASP_FORCE_OPTION_NBANDS_EQUAL.options2entry(AflowIn,_STROPT_+"NBANDS=",TRUE,vflags.KBIN_VASP_FORCE_OPTION_NBANDS_EQUAL.xscheme); // scheme already loaded in aflow_xclasses.cpp is "0"  
        // [OBSOLETE]  vflags.KBIN_VASP_FORCE_OPTION_NBANDS_EQUAL_isentry  = aurostd::substring2bool(AflowIn,_STROPT_+"NBANDS=",TRUE);
        if(vflags.KBIN_VASP_FORCE_OPTION_NBANDS_EQUAL.isentry) vflags.KBIN_VASP_FORCE_OPTION_NBANDS_AUTO_isentry=FALSE;
        // cerr << "vflags.KBIN_VASP_FORCE_OPTION_NBANDS_AUTO_isentry=" << vflags.KBIN_VASP_FORCE_OPTION_NBANDS_AUTO_isentry << endl;
        // cerr << "vflags.KBIN_VASP_FORCE_OPTION_NBANDS_EQUAL.isentry=" << vflags.KBIN_VASP_FORCE_OPTION_NBANDS_EQUAL.isentry << endl;
        // cerr << "vflags.KBIN_VASP_FORCE_OPTION_NBANDS_EQUAL.xscheme=" << vflags.KBIN_VASP_FORCE_OPTION_NBANDS_EQUAL.xscheme << endl;
        // exit(0);

        // cerr << "ENMAX_MULTIPLY" << endl;
        vflags.KBIN_VASP_FORCE_OPTION_ENMAX_MULTIPLY_EQUAL.options2entry(AflowIn,_STROPT_+"ENMAX_MULTIPLY=",TRUE,vflags.KBIN_VASP_FORCE_OPTION_ENMAX_MULTIPLY_EQUAL.xscheme); // scheme already loaded in aflow_xclasses.cpp is "0.0"  

        // RWIGS_STATIC
        // cerr << "RWIGS" << endl;
        vflags.KBIN_VASP_FORCE_OPTION_RWIGS_STATIC   =
            aurostd::substring2bool(AflowIn,_STROPT_+"RWIGS_STATIC",TRUE);

        // SPIN AND PRIORITIES // ON | OFF
        vflags.KBIN_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_1=DEFAULT_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_1;   // DEFAULT
        vflags.KBIN_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_2=DEFAULT_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_2;   // DEFAULT
        vflags.KBIN_VASP_FORCE_OPTION_SPIN.options2entry(AflowIn,_STROPT_+"SPIN=",DEFAULT_VASP_FORCE_OPTION_SPIN);
        if(aurostd::substring2bool(vflags.KBIN_VASP_FORCE_OPTION_SPIN.content_string,"REMOVE_RELAX_1")) vflags.KBIN_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_1=TRUE;
        if(aurostd::substring2bool(vflags.KBIN_VASP_FORCE_OPTION_SPIN.content_string,"REMOVE_RELAX_2")) vflags.KBIN_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_2=TRUE;
        if(!vflags.KBIN_VASP_FORCE_OPTION_SPIN.option) vflags.KBIN_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_1=FALSE; // nothing to remove
        if(!vflags.KBIN_VASP_FORCE_OPTION_SPIN.option) vflags.KBIN_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_2=FALSE; // nothing to remove

        // BADER AND PRIORITIES // ON | OFF
        vflags.KBIN_VASP_FORCE_OPTION_BADER.options2entry(AflowIn,_STROPT_+"BADER=",DEFAULT_VASP_FORCE_OPTION_BADER);
        if(vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_BANDS")) vflags.KBIN_VASP_FORCE_OPTION_BADER.isentry=TRUE; // DEFAULT 
        if(vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_BANDS")) vflags.KBIN_VASP_FORCE_OPTION_BADER.option=TRUE; // DEFAULT 
        // if(vflags.KBIN_VASP_RUN.flag("RELAX_STATIC")) vflags.KBIN_VASP_FORCE_OPTION_BADER.option=TRUE; // DEFAULT 
        // if(vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_BANDS")) vflags.KBIN_VASP_FORCE_OPTION_BADER.option=TRUE; // DEFAULT 
        // if(vflags.KBIN_VASP_RUN.flag("STATIC")) vflags.KBIN_VASP_FORCE_OPTION_BADER.option=TRUE; // DEFAULT 
        // if(vflags.KBIN_VASP_DIELECTRIC_STATIC) vflags.KBIN_VASP_FORCE_OPTION_BADER.option=TRUE; // DEFAULT 
        // if(vflags.KBIN_VASP_DIELECTRIC_DYNAMIC) vflags.KBIN_VASP_FORCE_OPTION_BADER.option=TRUE; // DEFAULT 
        // if(vflags.KBIN_VASP_DSCF) vflags.KBIN_VASP_FORCE_OPTION_BADER.option=TRUE; // DEFAULT 

        // ELF AND PRIORITIES // ON | OFF
        vflags.KBIN_VASP_FORCE_OPTION_ELF.options2entry(AflowIn,_STROPT_+"ELF=",DEFAULT_VASP_FORCE_OPTION_ELF);
        //  if(vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_BANDS")) vflags.KBIN_VASP_FORCE_OPTION_ELF.isentry=TRUE; // DEFAULT 
        // if(vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_BANDS")) vflags.KBIN_VASP_FORCE_OPTION_ELF.option=TRUE; // DEFAULT 

        // AUTO_MAGMOM AND PRIORITIES  // ON | OFF
        vflags.KBIN_VASP_FORCE_OPTION_AUTO_MAGMOM.options2entry(AflowIn,_STROPT_+"AUTO_MAGMOM=",DEFAULT_VASP_FORCE_OPTION_AUTO_MAGMOM);
        // LSCOUPLING AND PRIORITIES  // ON | OFF
        vflags.KBIN_VASP_FORCE_OPTION_LSCOUPLING.options2entry(AflowIn,_STROPT_+"LSCOUPLING=",DEFAULT_VASP_FORCE_OPTION_LSCOUPLING);
        if(vflags.KBIN_VASP_FORCE_OPTION_LSCOUPLING.option) {
            if(!aurostd::substring2bool(kflags.KBIN_BIN,"LS") && !aurostd::substring2bool(kflags.KBIN_BIN,"ls")) kflags.KBIN_BIN+="LS";
            if(!aurostd::substring2bool(kflags.KBIN_MPI_BIN,"LS") && !aurostd::substring2bool(kflags.KBIN_MPI_BIN,"ls")) kflags.KBIN_MPI_BIN+="LS";
            kflags.KBIN_BIN=aurostd::RemoveCharacter(kflags.KBIN_BIN,' ');            // if there is junk
            kflags.KBIN_MPI_BIN=aurostd::RemoveCharacter(kflags.KBIN_MPI_BIN,' ');    // if there is junk
        }
        // SYM AND PRIORITIES
        vflags.KBIN_VASP_FORCE_OPTION_SYM.options2entry(AflowIn,_STROPT_+"SYM=",DEFAULT_VASP_FORCE_OPTION_SYM);
        // WAVECAR AND PRIORITIES
        vflags.KBIN_VASP_FORCE_OPTION_WAVECAR.options2entry(AflowIn,_STROPT_+"WAVECAR=",DEFAULT_VASP_FORCE_OPTION_WAVECAR);
        // CHGCAR AND PRIORITIES
        vflags.KBIN_VASP_FORCE_OPTION_CHGCAR.options2entry(AflowIn,_STROPT_+"CHGCAR=",DEFAULT_VASP_FORCE_OPTION_CHGCAR);

        // LDAU2 AND PRIORITIES
        vflags.KBIN_VASP_LDAU_SPECIES="";
        vflags.KBIN_VASP_LDAU_PARAMETERS="";
        vflags.KBIN_VASP_LDAU_AFLOW_AUTO_flag=TRUE;

        BflowIn=AflowIn;aurostd::StringSubst(BflowIn,"LDAU1=","LDAU=");aurostd::StringSubst(BflowIn,"LDAU2=","LDAU=");
        vflags.KBIN_VASP_FORCE_OPTION_LDAU0.options2entry(BflowIn,string(_STROPT_+"LDAU=OFF"+"|"+_STROPT_+"LDAU=0"+"|"+_STROPT_+"LDAU=N"+"|"+_STROPT_+"LDAU=FALSE"));
        vflags.KBIN_VASP_FORCE_OPTION_LDAU1.options2entry(AflowIn,string(_STROPT_+"LDAU1=ON"+"|"+_STROPT_+"LDAU1=1"+"|"+"LDAU1=Y"+"|"+_STROPT_+"LDAU1=TRUE"+"|"+_STROPT_+"LDAU1=ADIABATIC"+"|"+_STROPT_+"LDAU1=CUTOFF"));
        vflags.KBIN_VASP_FORCE_OPTION_LDAU2.options2entry(AflowIn,string(_STROPT_+"LDAU2=ON"+"|"+_STROPT_+"LDAU2=1"+"|"+"LDAU2=Y"+"|"+_STROPT_+"LDAU2=TRUE"+"|"+_STROPT_+"LDAU2=ADIABATIC"+"|"+_STROPT_+"LDAU2=CUTOFF"));
        if(vflags.KBIN_VASP_FORCE_OPTION_LDAU1.isentry || vflags.KBIN_VASP_FORCE_OPTION_LDAU2.isentry)  vflags.KBIN_VASP_FORCE_OPTION_LDAU0.isentry=FALSE;
        if(vflags.KBIN_VASP_FORCE_OPTION_LDAU1.isentry || vflags.KBIN_VASP_FORCE_OPTION_LDAU2.isentry) {
            if(aurostd::substring2bool(AflowIn,_STROPT_+"LDAU_SPECIES=",TRUE))
                vflags.KBIN_VASP_LDAU_SPECIES=aurostd::substring2string(AflowIn,_STROPT_+"LDAU_SPECIES=",FALSE);
            if(aurostd::substring2bool(AflowIn,_STROPT_+"LDAU1_SPECIES=",TRUE))
                vflags.KBIN_VASP_LDAU_SPECIES=aurostd::substring2string(AflowIn,_STROPT_+"LDAU1_SPECIES=",FALSE);
            if(aurostd::substring2bool(AflowIn,_STROPT_+"LDAU2_SPECIES=",TRUE))
                vflags.KBIN_VASP_LDAU_SPECIES=aurostd::substring2string(AflowIn,_STROPT_+"LDAU2_SPECIES=",FALSE);
            if(aurostd::substring2bool(AflowIn,_STROPT_+"LDAU_PARAMETERS=",TRUE)) 
                vflags.KBIN_VASP_LDAU_PARAMETERS=RemoveWhiteSpaces(aurostd::substring2string(AflowIn,_STROPT_+"LDAU_PARAMETERS=",FALSE));
            if(vflags.KBIN_VASP_LDAU_SPECIES!="") 
                vflags.KBIN_VASP_LDAU_AFLOW_AUTO_flag=TRUE;
            if(vflags.KBIN_VASP_LDAU_PARAMETERS!="") 
                vflags.KBIN_VASP_LDAU_AFLOW_AUTO_flag=FALSE;
        }
        // ADIABATIC
        vflags.KBIN_VASP_FORCE_OPTION_LDAU_ADIABATIC.options2entry(AflowIn,string(_STROPT_+"LDAU1=ADIABATIC"+"|"+_STROPT_+"LDAU2=ADIABATIC"+"|"+_STROPT_+"LDAU=ADIABATIC"));
        if(vflags.KBIN_VASP_FORCE_OPTION_LDAU_ADIABATIC.isentry) {
            if(vflags.KBIN_VASP_RUN_NRELAX<LDAU_ADIABATIC_RELAX_DEFAULT)
                vflags.KBIN_VASP_RUN_NRELAX=LDAU_ADIABATIC_RELAX_DEFAULT;
            vflags.KBIN_VASP_FORCE_OPTION_LDAU_ADIABATIC.content_int=vflags.KBIN_VASP_RUN_NRELAX;
        }
        // CUTOFF
        vflags.KBIN_VASP_FORCE_OPTION_LDAU_CUTOFF.options2entry(AflowIn,string(_STROPT_+"LDAU1=CUTOFF"+"|"+_STROPT_+"LDAU2=CUTOFF"+"|"+_STROPT_+"LDAU=CUTOFF"));
        if(vflags.KBIN_VASP_FORCE_OPTION_LDAU_CUTOFF.isentry) {
            vflags.KBIN_VASP_RUN_NRELAX++;
        }
        // KPOINTS
        BflowIn=AflowIn;aurostd::StringSubst(BflowIn,"=","_");aurostd::StringSubst(BflowIn,"KPOINTS_","KPOINTS="); // bypass for getting all "_"
        vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.options2entry(BflowIn,string(_STROPT_+"KPOINTS="),aurostd_xoptionMULTI,""); // stack them all
        if(0) {
            cerr << "vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.content_string=" << vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.content_string << endl;
            cerr << "vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.flag(\"KEEPK\")=" << vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.flag("KEEPK") << endl;
            cerr << "vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.flag(\"EVEN\")=" << vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.flag("EVEN") << endl;
            cerr << "vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.flag(\"ODD\")=" << vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.flag("ODD") << endl;
            cerr << "vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.flag(\"KSHIFT_GAMMA_EVEN\")=" << vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.flag("KSHIFT_GAMMA_EVEN") << endl;
            cerr << "vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.flag(\"KSHIFT_GAMMA_ODD\")=" << vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.flag("KSHIFT_GAMMA_ODD") << endl;
            cerr << "vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.flag(\"GAMMA\")=" << vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.flag("GAMMA") << endl;
            cerr << "vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.flag(\"KSCHEME_MONKHORST_PACK\")=" << vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.flag("KSCHEME_MONKHORST_PACK") << endl;
            cerr << "vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.flag(\"KSCHEME_GAMMA\")=" << vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.flag("KSCHEME_GAMMA") << endl;
            cerr << "vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.flag(\"KSCHEME_AUTO\")=" << vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.flag("KSCHEME_AUTO") << endl;
            cerr << "vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.flag(\"IBZKPT\")=" << vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.flag("IBZKPT") << endl;
        }

        // TYPE AND PRIORITIES // METAL | INSULATOR | SEMICONDUCTOR | DEFAULT
        vflags.KBIN_VASP_FORCE_OPTION_TYPE.options2entry(AflowIn,_STROPT_+"TYPE=",FALSE,DEFAULT_VASP_FORCE_OPTION_TYPE_SCHEME);
        vflags.KBIN_VASP_FORCE_OPTION_TYPE.xscheme=KBIN_WRONG_ENTRY_STRING;
        vflags.KBIN_VASP_FORCE_OPTION_TYPE.scheme2scheme('M',"METAL");
        vflags.KBIN_VASP_FORCE_OPTION_TYPE.scheme2scheme('I',"INSULATOR");
        vflags.KBIN_VASP_FORCE_OPTION_TYPE.scheme2scheme('S',"SEMICONDUCTOR");
        vflags.KBIN_VASP_FORCE_OPTION_TYPE.scheme2scheme('D',"DEFAULT");
        if(vflags.KBIN_VASP_FORCE_OPTION_TYPE.isentry && vflags.KBIN_VASP_FORCE_OPTION_TYPE.xscheme==KBIN_WRONG_ENTRY_STRING) {
            cerr << "ERROR: KBIN::VASP_Get_Vflags_from_AflowIN:  vflags.KBIN_VASP_FORCE_OPTION_TYPE.content_string="  <<  vflags.KBIN_VASP_FORCE_OPTION_TYPE.content_string << endl;
            exit(0);
        }

        // PARAMETERS FOR INCAR
        // NSW=
        vflags.KBIN_VASP_FORCE_OPTION_NSW_EQUAL =
            aurostd::substring2bool(AflowIn,_STROPT_+"NSW=",TRUE);
        if(vflags.KBIN_VASP_FORCE_OPTION_NSW_EQUAL)
            vflags.KBIN_VASP_FORCE_OPTION_NSW_EQUAL_VALUE=aurostd::substring2utype<int>(AflowIn,_STROPT_+"NSW=",-1);
        else
            vflags.KBIN_VASP_FORCE_OPTION_NSW_EQUAL_VALUE=0;

        // IGNORE_AFIX stuff
        BflowIn=AflowIn;aurostd::StringSubst(BflowIn,"=","_");aurostd::StringSubst(BflowIn,"IGNORE_AFIX_","IGNORE_AFIX="); // bypass for getting all "_"
        vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.options2entry(BflowIn,string(_STROPT_+"IGNORE_AFIX="),aurostd_xoptionMULTI,""); // stack them all
        if(0) {
            cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"BRMIX\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("BRMIX") << endl;
            cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"CSLOSHING\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("CSLOSHING") << endl;
            cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"DAV\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("DAV") << endl;
            cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"DENTET\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("DENTET") << endl;
            cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"EDDDAV\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("EDDDAV") << endl;
            cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"EDDRMM\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("EDDRMM") << endl;
            cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"EXCCOR\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("EXCCOR") << endl;
            cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"GAMMA_SHIFT\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("GAMMA_SHIFT") << endl;
            cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"IBZKPT\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("IBZKPT") << endl;
            cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"INVGRP\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("INVGRP") << endl;
            cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"LREAL\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("LREAL") << endl;
            cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"LRF_COMMUTATOR\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("LRF_COMMUTATOR") << endl;
            cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"MEMORY\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("MEMORY") << endl;
            cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"MPICH11\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("MPICH11") << endl;
            cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"MPICH139\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("MPICH139") << endl;
            cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"NATOMS\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("NATOMS") << endl;
            cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"NBANDS\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("NBANDS") << endl;
            cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"NKXYZ_IKPTD\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("NKXYZ_IKPTD") << endl;
            cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"NONE\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("NONE") << endl;
            cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"NPAR\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("NPAR") << endl;
            cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"NPARC\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("NPARC") << endl;
            cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"NPARN\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("NPARN") << endl;
            cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"NPAR_REMOVE\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("NPAR_REMOVE") << endl;
            cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"PSMAXN\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("PSMAXN") << endl;
            cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"READ_KPOINTS_RD_SYM\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("READ_KPOINTS_RD_SYM") << endl;
            cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"ROTMAT\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("ROTMAT") << endl;
            cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"SGRCON\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("SGRCON") << endl;
            cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"SYMPREC\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("SYMPREC") << endl;
            cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"ZPOTRF\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("ZPOTRF") << endl;
        }

        // INPUT FILES
        // [OBSOLETE] vflags.KBIN_VASP_INCAR_FILE_KEYWORD                =    aurostd::substring2bool(AflowIn,"[VASP_INCAR_FILE]");      
        // [OBSOLETE] vflags.KBIN_VASP_INCAR_FILE_SYSTEM_AUTO            =    aurostd::substring2bool(AflowIn,"[VASP_INCAR_FILE]SYSTEM_AUTO",TRUE);
        // [OBSOLETE] vflags.KBIN_VASP_INCAR_FILE_FILE                   =    aurostd::substring2bool(AflowIn,"[VASP_INCAR_FILE]FILE=",TRUE);
        // [OBSOLETE] vflags.KBIN_VASP_INCAR_FILE_COMMAND                =    aurostd::substring2bool(AflowIn,"[VASP_INCAR_FILE]COMMAND=",TRUE);

        if(aurostd::substring2bool(AflowIn,"[VASP_INCAR_FILE]")) vflags.KBIN_VASP_INCAR_FILE.push("KEYWORD");
        if(aurostd::substring2bool(AflowIn,"[VASP_INCAR_FILE]SYSTEM_AUTO",TRUE)) vflags.KBIN_VASP_INCAR_FILE.push("SYSTEM_AUTO");
        if(aurostd::substring2bool(AflowIn,"[VASP_INCAR_FILE]FILE=",TRUE)) vflags.KBIN_VASP_INCAR_FILE.push("FILE");
        if(aurostd::substring2bool(AflowIn,"[VASP_INCAR_FILE]COMMAND=",TRUE)) vflags.KBIN_VASP_INCAR_FILE.push("COMMAND");

        // [OBSOLETE] vflags.KBIN_VASP_INCAR_MODE_EXPLICIT             =    aurostd::substring2bool(AflowIn,"[VASP_INCAR_MODE_EXPLICIT]");
        // [OBSOLETE] vflags.KBIN_VASP_INCAR_MODE_EXPLICIT_START_STOP  =    aurostd::substring2bool(AflowIn,"[VASP_INCAR_MODE_EXPLICIT]START") && aurostd::substring2bool(AflowIn,"[VASP_INCAR_MODE_EXPLICIT]STOP");
        // [OBSOLETE] vflags.KBIN_VASP_INCAR_MODE_IMPLICIT             =    aurostd::substring2bool(AflowIn,"[VASP_INCAR_MODE_IMPLICIT]");
        // [OBSOLETE] vflags.KBIN_VASP_INCAR_MODE_EXTERNAL             =    aurostd::substring2bool(AflowIn,"[VASP_INCAR_MODE_EXTERNAL]");
        if(aurostd::substring2bool(AflowIn,"[VASP_INCAR_MODE_EXPLICIT]")) vflags.KBIN_VASP_INCAR_MODE.push("EXPLICIT");
        if(aurostd::substring2bool(AflowIn,"[VASP_INCAR_MODE_EXPLICIT]START") && aurostd::substring2bool(AflowIn,"[VASP_INCAR_MODE_EXPLICIT]STOP")) vflags.KBIN_VASP_INCAR_MODE.push("EXPLICIT_START_STOP");
        if(aurostd::substring2bool(AflowIn,"[VASP_INCAR_MODE_IMPLICIT]")) vflags.KBIN_VASP_INCAR_MODE.push("IMPLICIT");
        if(aurostd::substring2bool(AflowIn,"[VASP_INCAR_MODE_EXTERNAL]")) vflags.KBIN_VASP_INCAR_MODE.push("EXTERNAL");

        // [OBSOLETE] vflags.KBIN_VASP_KPOINTS_FILE_KEYWORD                      =    aurostd::substring2bool(AflowIn,"[VASP_KPOINTS_FILE]");      
        if(aurostd::substring2bool(AflowIn,"[VASP_KPOINTS_FILE]")) vflags.KBIN_VASP_KPOINTS_FILE.push("KEYWORD");
        //  vflags.KBIN_VASP_KPOINTS_FILE_SYSTEM_AUTO          =  aurostd::substring2bool(AflowIn,"[VASP_KPOINTS_FILE]SYSTEM_AUTO",TRUE);

        // [OBSOLETE] vflags.KBIN_VASP_KPOINTS_MODE_EXPLICIT             =    aurostd::substring2bool(AflowIn,"[VASP_KPOINTS_MODE_EXPLICIT]");
        // [OBSOLETE] vflags.KBIN_VASP_KPOINTS_MODE_EXPLICIT_START_STOP  =    aurostd::substring2bool(AflowIn,"[VASP_KPOINTS_MODE_EXPLICIT]START") && aurostd::substring2bool(AflowIn,"[VASP_KPOINTS_MODE_EXPLICIT]STOP");
        // [OBSOLETE] vflags.KBIN_VASP_KPOINTS_MODE_IMPLICIT             =    aurostd::substring2bool(AflowIn,"[VASP_KPOINTS_MODE_IMPLICIT]");
        // [OBSOLETE] vflags.KBIN_VASP_KPOINTS_MODE_EXTERNAL             =    aurostd::substring2bool(AflowIn,"[VASP_KPOINTS_MODE_EXTERNAL]");
        if(aurostd::substring2bool(AflowIn,"[VASP_KPOINTS_MODE_EXPLICIT]")) vflags.KBIN_VASP_KPOINTS_MODE.push("EXPLICIT");
        if(aurostd::substring2bool(AflowIn,"[VASP_KPOINTS_MODE_EXPLICIT]START") && aurostd::substring2bool(AflowIn,"[VASP_KPOINTS_MODE_EXPLICIT]STOP")) vflags.KBIN_VASP_KPOINTS_MODE.push("EXPLICIT_START_STOP");
        if(aurostd::substring2bool(AflowIn,"[VASP_KPOINTS_MODE_IMPLICIT]")) vflags.KBIN_VASP_KPOINTS_MODE.push("IMPLICIT");
        if(aurostd::substring2bool(AflowIn,"[VASP_KPOINTS_MODE_EXTERNAL]")) vflags.KBIN_VASP_KPOINTS_MODE.push("EXTERNAL");
        // [OBSOLETE] vflags.KBIN_VASP_KPOINTS_FILE_FILE                 =    aurostd::substring2bool(AflowIn,"[VASP_KPOINTS_FILE]FILE=",TRUE);
        // [OBSOLETE] vflags.KBIN_VASP_KPOINTS_FILE_COMMAND              =    aurostd::substring2bool(AflowIn,"[VASP_KPOINTS_FILE]COMMAND=",TRUE);

        if(aurostd::substring2bool(AflowIn,"[VASP_KPOINTS_FILE]FILE=",TRUE)) vflags.KBIN_VASP_KPOINTS_FILE.push("FILE");
        if(aurostd::substring2bool(AflowIn,"[VASP_KPOINTS_FILE]COMMAND=",TRUE)) vflags.KBIN_VASP_KPOINTS_FILE.push("COMMAND");

        // KPOINTS FOR RELAX
        vflags.KBIN_VASP_KPOINTS_KMODE.options2entry(AflowIn,"[VASP_KPOINTS_FILE]KMODE=",FALSE,vflags.KBIN_VASP_KPOINTS_KMODE.xscheme); // scheme already loaded in aflow_xclasses.cpp is "0"  
        // cerr << "vflags.KBIN_VASP_KPOINTS_KMODE.isentry=" << vflags.KBIN_VASP_KPOINTS_KMODE.isentry << endl << "vflags.KBIN_VASP_KPOINTS_KMODE.xscheme=" << vflags.KBIN_VASP_KPOINTS_KMODE.xscheme << endl;
        vflags.KBIN_VASP_KPOINTS_KPPRA.options2entry(AflowIn,"[VASP_KPOINTS_FILE]KPPRA=",FALSE,vflags.KBIN_VASP_KPOINTS_KPPRA.xscheme); // scheme already loaded in aflow_xclasses.cpp is "1"
        if(vflags.KBIN_VASP_KPOINTS_KPPRA.isentry==FALSE) {vflags.KBIN_VASP_KPOINTS_KPPRA.clear();vflags.KBIN_VASP_KPOINTS_KPPRA.isentry=TRUE;vflags.KBIN_VASP_KPOINTS_KPPRA.push("100");}
        // cerr << "vflags.KBIN_VASP_KPOINTS_KPPRA.isentry=" << vflags.KBIN_VASP_KPOINTS_KPPRA.isentry << endl << "vflags.KBIN_VASP_KPOINTS_KPPRA.xscheme=" << vflags.KBIN_VASP_KPOINTS_KPPRA.xscheme << endl;
        vflags.KBIN_VASP_KPOINTS_KSCHEME.options2entry(AflowIn,"[VASP_KPOINTS_FILE]KSCHEME=",FALSE,vflags.KBIN_VASP_KPOINTS_KSCHEME.xscheme); // scheme already loaded in aflow_xclasses.cpp is "1"
        if(vflags.KBIN_VASP_KPOINTS_KSCHEME.isentry==FALSE) {vflags.KBIN_VASP_KPOINTS_KSCHEME.clear();vflags.KBIN_VASP_KPOINTS_KSCHEME.isentry=TRUE;vflags.KBIN_VASP_KPOINTS_KSCHEME.push(DEFAULT_KSCHEME);}
        // cerr << "vflags.KBIN_VASP_KPOINTS_KSCHEME.isentry=" << vflags.KBIN_VASP_KPOINTS_KSCHEME.isentry << endl << "vflags.KBIN_VASP_KPOINTS_KSCHEME.xscheme=" << vflags.KBIN_VASP_KPOINTS_KSCHEME.xscheme << endl;
        vflags.KBIN_VASP_KPOINTS_KSHIFT.options2entry(AflowIn,"[VASP_KPOINTS_FILE]KSHIFT=",FALSE,vflags.KBIN_VASP_KPOINTS_KSHIFT.xscheme); // scheme already loaded in aflow_xclasses.cpp is "0 0 0"
        // cerr << "vflags.KBIN_VASP_KPOINTS_KSHIFT.isentry=" << vflags.KBIN_VASP_KPOINTS_KSHIFT.isentry << endl << "vflags.KBIN_VASP_KPOINTS_KSHIFT.xscheme=" << vflags.KBIN_VASP_KPOINTS_KSHIFT.xscheme << endl;

        // KPOINTS FOR STATIC
        vflags.KBIN_VASP_KPOINTS_STATIC_KMODE.options2entry(AflowIn,"[VASP_KPOINTS_FILE]STATIC_KMODE=",FALSE,vflags.KBIN_VASP_KPOINTS_STATIC_KMODE.xscheme); // scheme already loaded in aflow_xclasses.cpp is "0"
        // cerr << "vflags.KBIN_VASP_KPOINTS_STATIC_KMODE.isentry=" << vflags.KBIN_VASP_KPOINTS_STATIC_KMODE.isentry << endl << "vflags.KBIN_VASP_KPOINTS_STATIC_KMODE.xscheme=" << vflags.KBIN_VASP_KPOINTS_STATIC_KMODE.xscheme << endl;
        vflags.KBIN_VASP_KPOINTS_STATIC_KPPRA.options2entry(AflowIn,"[VASP_KPOINTS_FILE]STATIC_KPPRA=",FALSE,vflags.KBIN_VASP_KPOINTS_STATIC_KPPRA.xscheme); // scheme already loaded in aflow_xclasses.cpp is "1"
        // cerr << "vflags.KBIN_VASP_KPOINTS_STATIC_KPPRA.isentry=" << vflags.KBIN_VASP_KPOINTS_STATIC_KPPRA.isentry << endl << "vflags.KBIN_VASP_KPOINTS_STATIC_KPPRA.xscheme=" << vflags.KBIN_VASP_KPOINTS_STATIC_KPPRA.xscheme << endl;
        vflags.KBIN_VASP_KPOINTS_STATIC_KSCHEME.options2entry(AflowIn,"[VASP_KPOINTS_FILE]STATIC_KSCHEME=",FALSE,vflags.KBIN_VASP_KPOINTS_STATIC_KSCHEME.xscheme); // scheme already loaded in aflow_xclasses.cpp is "Monkhorst-Pack"
        // cerr << "vflags.KBIN_VASP_KPOINTS_STATIC_KSCHEME.isentry=" << vflags.KBIN_VASP_KPOINTS_STATIC_KSCHEME.isentry << endl << "vflags.KBIN_VASP_KPOINTS_STATIC_KSCHEME.xscheme=" << vflags.KBIN_VASP_KPOINTS_STATIC_KSCHEME.xscheme << endl;
        vflags.KBIN_VASP_KPOINTS_STATIC_KSHIFT.options2entry(AflowIn,"[VASP_KPOINTS_FILE]STATIC_KSHIFT=",FALSE,vflags.KBIN_VASP_KPOINTS_STATIC_KSHIFT.xscheme); // scheme already loaded in aflow_xclasses.cpp is "0 0 0"
        //cerr << "vflags.KBIN_VASP_KPOINTS_STATIC_KSHIFT.isentry=" << vflags.KBIN_VASP_KPOINTS_STATIC_KSHIFT.isentry << endl << "vflags.KBIN_VASP_KPOINTS_STATIC_KSHIFT.xscheme=" << vflags.KBIN_VASP_KPOINTS_STATIC_KSHIFT.xscheme << endl;
        //  exit(0);

        // [OBSOLETE] vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE_FLAG        =    aurostd::substring2bool(AflowIn,"[VASP_KPOINTS_FILE]BANDS_LATTICE=",TRUE);
        // [OBSOLETE] if(vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE_FLAG) vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE_VALUE=aurostd::RemoveWhiteSpaces(aurostd::substring2string(AflowIn,"[VASP_KPOINTS_FILE]BANDS_LATTICE=",TRUE));
        // [OBSOLETE]  else {vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE_VALUE=DEFAULT_BANDS_LATTICE;vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE_AUTO_FLAG=TRUE;} // DEFAULT FIX

        vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE.options2entry(AflowIn,"[VASP_KPOINTS_FILE]BANDS_LATTICE=",FALSE,vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE.xscheme); // scheme already loaded in aflow_xclasses.cpp is ""
        if(!vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE.isentry &&
                !vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE.flag(DEFAULT_BANDS_LATTICE))
            vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE.push(DEFAULT_BANDS_LATTICE);

        if((vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_BANDS") || vflags.KBIN_VASP_RUN.flag("STATIC_BANDS")) && !vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE.isentry) {
            cerr << "WARNING: if you use vflags.KBIN_VASP_RUN.flag(\"RELAX_STATIC_BANDS\") or vflags.KBIN_VASP_RUN.flag(\"STATIC_BANDS\"), you must specify KBIN_VASP_KPOINTS_BANDS_LATTICE" << endl;
            cerr << "         Taking defauls vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE.content_string=" << vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE.content_string << endl;
        }
        vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE_AUTO_FLAG    =    aurostd::substring2bool(AflowIn,"[VASP_KPOINTS_FILE]BANDS_LATTICE=AUTO",TRUE) ||    aurostd::substring2bool(AflowIn,"[VASP_KPOINTS_FILE]BANDS_LATTICE=auto",TRUE);
        // [OBSOLETE]  if(vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE_AUTO_FLAG) vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE_FLAG=FALSE;
        if(vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE.content_string=="AUTO") vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE_AUTO_FLAG=TRUE;
        if(vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE_AUTO_FLAG) vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE.clear();

        vflags.KBIN_VASP_KPOINTS_BANDS_GRID_FLAG    =     aurostd::substring2bool(AflowIn,"[VASP_KPOINTS_FILE]BANDS_GRID=",TRUE);
        if(vflags.KBIN_VASP_KPOINTS_BANDS_GRID_FLAG)
            vflags.KBIN_VASP_KPOINTS_BANDS_GRID_VALUE=aurostd::substring2utype<int>(AflowIn,"[VASP_KPOINTS_FILE]BANDS_GRID=",TRUE);
        else {vflags.KBIN_VASP_KPOINTS_BANDS_GRID_VALUE=16;vflags.KBIN_VASP_KPOINTS_BANDS_GRID_FLAG=TRUE;
            //    cerr << "WARNING: if you use VASP_RUN_RELAX_STATIC_BANDS or vflags.KBIN_VASP_RUN.flag(\"STATIC_BANDS\"), you must specify KBIN_VASP_KPOINTS_BANDS_GRID_FLAG" << endl;
            //   cerr << "         Taking defauls vflags.KBIN_VASP_KPOINTS_BANDS_GRID_VALUE=" << vflags.KBIN_VASP_KPOINTS_BANDS_GRID_VALUE << endl;
        } // DEFAULT FIX
        if((vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_BANDS") || vflags.KBIN_VASP_RUN.flag("STATIC_BANDS")) && !vflags.KBIN_VASP_KPOINTS_BANDS_GRID_FLAG) {
            cerr << "WARNING: if you use VASP_RUN_RELAX_STATIC_BANDS or vflags.KBIN_VASP_RUN.flag(\"STATIC_BANDS\"), you must specify KBIN_VASP_KPOINTS_BANDS_GRID_FLAG" << endl;
            cerr << "         Taking defauls vflags.KBIN_VASP_KPOINTS_BANDS_GRID_VALUE=" << vflags.KBIN_VASP_KPOINTS_BANDS_GRID_VALUE << endl;
        }

        // [OBSOLETE] vflags.KBIN_VASP_POSCAR_FILE_KEYWORD                 =    aurostd::substring2bool(AflowIn,"[VASP_POSCAR_FILE]"); 
        if(aurostd::substring2bool(AflowIn,"[VASP_POSCAR_FILE]"))  vflags.KBIN_VASP_POSCAR_FILE.push("KEYWORD");

        // [OBSOLETE] vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT                  =    aurostd::substring2bool(AflowIn,"[VASP_POSCAR_MODE_EXPLICIT]");
        if(aurostd::substring2bool(AflowIn,"[VASP_POSCAR_MODE_EXPLICIT]")) vflags.KBIN_VASP_POSCAR_MODE.push("EXPLICIT");

        // [OBSOLETE] vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_START_STOP       =    aurostd::substring2bool(AflowIn,"[VASP_POSCAR_MODE_EXPLICIT]START") &&  aurostd::substring2bool(AflowIn,"[VASP_POSCAR_MODE_EXPLICIT]STOP");
        if(aurostd::substring2bool(AflowIn,"[VASP_POSCAR_MODE_EXPLICIT]START") &&  aurostd::substring2bool(AflowIn,"[VASP_POSCAR_MODE_EXPLICIT]STOP")) vflags.KBIN_VASP_POSCAR_MODE.push("EXPLICIT_START_STOP");

        // [OBSOLETE] vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_START_STOP_POINT =    (aurostd::substring2bool(AflowIn,"[VASP_POSCAR_MODE_EXPLICIT]START.") && aurostd::substring2bool(AflowIn,"[VASP_POSCAR_MODE_EXPLICIT]STOP."));
        if((aurostd::substring2bool(AflowIn,"[VASP_POSCAR_MODE_EXPLICIT]START.") && aurostd::substring2bool(AflowIn,"[VASP_POSCAR_MODE_EXPLICIT]STOP."))) vflags.KBIN_VASP_POSCAR_MODE.push("EXPLICIT_START_STOP_POINT");

        if(vflags.KBIN_VASP_POSCAR_MODE.flag("EXPLICIT_START_STOP_POINT") && !kflags.KBIN_FROZSL) {  // NO FROZSL
            if(LDEBUG) cerr << "DEBUG: vflags.KBIN_VASP_POSCAR_MODE.flag(\"EXPLICIT_START_STOP_POINT\")=" << vflags.KBIN_VASP_POSCAR_MODE.flag("EXPLICIT_START_STOP_POINT") << endl;
            if(LDEBUG) cerr << "DEBUG: kflags.KBIN_PHONONS_CALCULATION_FROZSL=" << kflags.KBIN_PHONONS_CALCULATION_FROZSL << endl;
            if(LDEBUG) cerr << "DEBUG: kflags.KBIN_FROZSL_DOWNLOAD=" << kflags.KBIN_FROZSL_DOWNLOAD << endl;
            if(LDEBUG) cerr << "DEBUG: kflags.KBIN_FROZSL_FILE=" << kflags.KBIN_FROZSL_FILE << endl;
            stringstream input_file;
            input_file.clear();
            // loading
            if(vflags.KBIN_VASP_POSCAR_MODE.flag("EXPLICIT_START_STOP_POINT")) input_file.str(AflowIn);
            if(kflags.KBIN_PHONONS_CALCULATION_FROZSL) {
                FROZSL::Setup_frozsl_init_input(AflowIn,FileMESSAGE,input_file,aflags,kflags);
                FROZSL::Extract_INPUT(AflowIn,FileMESSAGE,input_file,aflags,kflags);
            }
            if(kflags.KBIN_FROZSL_DOWNLOAD)    FROZSL::Setup_frozsl_init_input(AflowIn,FileMESSAGE,input_file,aflags,kflags);
            if(kflags.KBIN_FROZSL_FILE)        FROZSL::File_INPUT(AflowIn,FileMESSAGE,input_file,aflags,kflags);

            vflags.KBIN_VASP_POSCAR_MODE.push("EXPLICIT_START_STOP_POINT");
            // done loading now load structures up
            vflags.KBIN_VASP_POSCAR_MODE.flag("EXPLICIT_START_STOP",FALSE); // some default
            aurostd::substring2strings(input_file.str(),vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING,"[VASP_POSCAR_MODE_EXPLICIT]START.");
            // some verbose
            for(uint i=0;i<vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING.size();i++)
                if(LDEBUG) cerr << "DEBUG= " << vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING.at(i) << endl;
            // load up the structures
            for(uint i=0;i<vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING.size();i++) {
                string START="[VASP_POSCAR_MODE_EXPLICIT]START";
                string STOP="[VASP_POSCAR_MODE_EXPLICIT]STOP";
                START="[VASP_POSCAR_MODE_EXPLICIT]START."+vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING.at(i);
                STOP="[VASP_POSCAR_MODE_EXPLICIT]STOP."+vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING.at(i);
                stringstream POSCAR;POSCAR.clear();POSCAR.str(std::string());
                if(aurostd::substring2bool(input_file.str(),START) && aurostd::substring2bool(input_file.str(),STOP))
                    aurostd::ExtractToStringstreamEXPLICIT(input_file.str(),POSCAR,START,STOP);
                vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRUCTURE.push_back(xstructure(POSCAR,IOVASP_AUTO));
            }
            if(LDEBUG) cerr << "DEBUG " << vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING.size() << endl;
            if(LDEBUG) cerr << "DEBUG " << vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRUCTURE.size() << endl;
            if(vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING.size() != vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRUCTURE.size()) {
                cerr << "ERROR (KBIN::VASP_Get_Vflags_from_AflowIN) IN " << _AFLOWIN_ << " in Directory=" << aflags.Directory << endl;
                cerr << "vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING.size()=" << vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING.size() << endl;
                cerr << "vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRUCTURE.size()=" << vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRUCTURE.size() << endl;
                exit(0);
            }
            for(uint i=0;i<vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING.size();i++)
                if(LDEBUG) cerr << "DEBUG= " << vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRUCTURE.at(i) << endl;
        } else {
            vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING.clear();
            vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRUCTURE.clear();
        }
        // the rest for POSCAR
        // [OBSOLETE]  vflags.KBIN_VASP_POSCAR_MODE_IMPLICIT             =    aurostd::substring2bool(AflowIn,"[VASP_POSCAR_MODE_IMPLICIT]");
        if(aurostd::substring2bool(AflowIn,"[VASP_POSCAR_MODE_IMPLICIT]")) vflags.KBIN_VASP_POSCAR_MODE.push("IMPLICIT");

        // vflags.KBIN_VASP_POSCAR_FILE_SYSTEM_AUTO                =   aurostd::substring2bool(AflowIn,"[VASP_POSCAR_FILE]SYSTEM_AUTO",TRUE);
        // [OBSOLETE] vflags.KBIN_VASP_POSCAR_FILE_PROTOTYPE                     =    aurostd::substring2bool(AflowIn,"[VASP_POSCAR_FILE]PROTOTYPE=",TRUE);
        if(aurostd::substring2bool(AflowIn,"[VASP_POSCAR_FILE]PROTOTYPE=",TRUE)) vflags.KBIN_VASP_POSCAR_FILE.push("PROTOTYPE");

        // [OBSOLETE]  vflags.KBIN_VASP_POSCAR_MODE_EXTERNAL                      =    aurostd::substring2bool(AflowIn,"[VASP_POSCAR_MODE_EXTERNAL]");
        if(aurostd::substring2bool(AflowIn,"[VASP_POSCAR_MODE_EXTERNAL]")) vflags.KBIN_VASP_POSCAR_MODE.push("EXTERNAL");

        // [OBSOLETE]  vflags.KBIN_VASP_POSCAR_FILE_FILE                          =    aurostd::substring2bool(AflowIn,"[VASP_POSCAR_FILE]FILE=",TRUE);
        if(aurostd::substring2bool(AflowIn,"[VASP_POSCAR_FILE]FILE=",TRUE)) vflags.KBIN_VASP_POSCAR_FILE.push("FILE");
        // [OBSOLETE]  vflags.KBIN_VASP_POSCAR_FILE_COMMAND                       =    aurostd::substring2bool(AflowIn,"[VASP_POSCAR_FILE]COMMAND=",TRUE);
        if(aurostd::substring2bool(AflowIn,"[VASP_POSCAR_FILE]COMMAND=",TRUE)) vflags.KBIN_VASP_POSCAR_FILE.push("COMMAND");

        // VOLUMES
        // [OBSOLETE] vflags.KBIN_VASP_POSCAR_FILE_VOLUME_EQUAL_EQUAL            =    aurostd::substring2bool(AflowIn,"[VASP_POSCAR_FILE]VOLUME=",TRUE);
        // [OBSOLETE] vflags.KBIN_VASP_POSCAR_FILE_VOLUME_PLUS_EQUAL             =    aurostd::substring2bool(AflowIn,"[VASP_POSCAR_FILE]VOLUME+=",TRUE);
        // [OBSOLETE] vflags.KBIN_VASP_POSCAR_FILE_VOLUME_MINUS_EQUAL            =    aurostd::substring2bool(AflowIn,"[VASP_POSCAR_FILE]VOLUME-=",TRUE);
        // [OBSOLETE] vflags.KBIN_VASP_POSCAR_FILE_VOLUME_MULTIPLY_EQUAL         =    aurostd::substring2bool(AflowIn,"[VASP_POSCAR_FILE]VOLUME*=",TRUE);
        // [OBSOLETE] vflags.KBIN_VASP_POSCAR_FILE_VOLUME_DIVIDE_EQUAL           =    aurostd::substring2bool(AflowIn,"[VASP_POSCAR_FILE]VOLUME/=",TRUE);
        vflags.KBIN_VASP_POSCAR_FILE_VOLUME.clear();
        if(aurostd::substring2bool(AflowIn,"[VASP_POSCAR_FILE]VOLUME=",TRUE)) vflags.KBIN_VASP_POSCAR_FILE_VOLUME.push("EQUAL_EQUAL");
        if(aurostd::substring2bool(AflowIn,"[VASP_POSCAR_FILE]VOLUME+=",TRUE)) vflags.KBIN_VASP_POSCAR_FILE_VOLUME.push("PLUS_EQUAL");
        // [OBSOLETE] if(aurostd::substring2bool(AflowIn,"[VASP_POSCAR_FILE]VOLUME-=",TRUE)) vflags.KBIN_VASP_POSCAR_FILE_VOLUME.push("MINUS_EQUAL");
        if(aurostd::substring2bool(AflowIn,"[VASP_POSCAR_FILE]VOLUME*=",TRUE)) vflags.KBIN_VASP_POSCAR_FILE_VOLUME.push("MULTIPLY_EQUAL");
        // [OBSOLETE] if(aurostd::substring2bool(AflowIn,"[VASP_POSCAR_FILE]VOLUME/=",TRUE)) vflags.KBIN_VASP_POSCAR_FILE_VOLUME.push("DIVIDE_EQUAL");
        if(vflags.KBIN_VASP_POSCAR_FILE_VOLUME.xscheme!="") vflags.KBIN_VASP_POSCAR_FILE_VOLUME.isentry=TRUE;

        // CONVERT_UNIT_CELL stuff
        BflowIn=AflowIn;aurostd::StringSubst(BflowIn,"=","_");aurostd::StringSubst(BflowIn,"CONVERT_UNIT_CELL_","CONVERT_UNIT_CELL="); // bypass for getting all "_"
        vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.options2entry(BflowIn,string(_STROPT_+"CONVERT_UNIT_CELL="),aurostd_xoptionMULTI,""); // stack them all
        // // PRIORITIES
        if(vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("STANDARD_PRIMITIVE") || vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("STANDARD_CONVENTIONAL")) {
            vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("MINKOWSKI",FALSE);
            vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("INCELL",FALSE);
            vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("COMPACT",FALSE);
            vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("WIGNERSEITZ",FALSE);
        } // some PRIORITIES
        if(vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("STANDARD_PRIMITIVE") && vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("STANDARD_CONVENTIONAL")) {
            vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("STANDARD_CONVENTIONAL",FALSE);
        }

        // DEBUG
        if(LDEBUG) cerr << "KBIN::VASP_Get_Vflags_from_AflowIN: vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.content_string=" << vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.content_string  << endl;
        if(LDEBUG) cerr << "KBIN::VASP_Get_Vflags_from_AflowIN: vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag(\"STANDARD_PRIMITIVE\")=" <<  vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("STANDARD_PRIMITIVE") << endl;
        if(LDEBUG) cerr << "KBIN::VASP_Get_Vflags_from_AflowIN: vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag(\"STANDARD_CONVENTIONAL\")=" <<  vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("STANDARD_CONVENTIONAL") << endl;
        if(LDEBUG) cerr << "KBIN::VASP_Get_Vflags_from_AflowIN: vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag(\"NIGGLI\")=" <<  vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("NIGGLI") << endl;
        if(LDEBUG) cerr << "KBIN::VASP_Get_Vflags_from_AflowIN: vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag(\"MINKOWSKI\")=" <<  vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("MINKOWSKI") << endl;
        if(LDEBUG) cerr << "KBIN::VASP_Get_Vflags_from_AflowIN: vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag(\"INCELL\")=" <<  vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("INCELL") << endl;
        if(LDEBUG) cerr << "KBIN::VASP_Get_Vflags_from_AflowIN: vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag(\"COMPACT\")=" <<  vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("COMPACT") << endl;
        if(LDEBUG) cerr << "KBIN::VASP_Get_Vflags_from_AflowIN: vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag(\"WIGNERSEITZ\")=" <<  vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("WIGNERSEITZ") << endl;
        if(LDEBUG) cerr << "KBIN::VASP_Get_Vflags_from_AflowIN: vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag(\"CARTESIAN\")=" <<  vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("CARTESIAN") << endl;
        if(LDEBUG) cerr << "KBIN::VASP_Get_Vflags_from_AflowIN: vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag(\"FRACTIONAL\")=" <<  vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("FRACTIONAL") << endl;
        if(LDEBUG) cerr << "KBIN::VASP_Get_Vflags_from_AflowIN: vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag(\"DIRECT\")=" <<  vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("DIRECT") << endl;
        if(LDEBUG) cerr << "KBIN::VASP_Get_Vflags_from_AflowIN: vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag(\"PRESERVE\")=" <<  vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("PRESERVE") << endl;

        // VOLUMES
        // [OBSOLETE] vflags.KBIN_VASP_FORCE_OPTION_VOLUME_EQUAL_EQUAL      =    aurostd::substring2bool(AflowIn,_STROPT_+"VOLUME=",TRUE);
        // [OBSOLETE] vflags.KBIN_VASP_FORCE_OPTION_VOLUME_PLUS_EQUAL       =    aurostd::substring2bool(AflowIn,_STROPT_+"VOLUME+=",TRUE);
        // [OBSOLETE] vflags.KBIN_VASP_FORCE_OPTION_VOLUME_MINUS_EQUAL      =    aurostd::substring2bool(AflowIn,_STROPT_+"VOLUME-=",TRUE);
        // [OBSOLETE] vflags.KBIN_VASP_FORCE_OPTION_VOLUME_MULTIPLY_EQUAL   =    aurostd::substring2bool(AflowIn,_STROPT_+"VOLUME*=",TRUE);
        // [OBSOLETE] vflags.KBIN_VASP_FORCE_OPTION_VOLUME_DIVIDE_EQUAL     =    aurostd::substring2bool(AflowIn,_STROPT_+"VOLUME/=",TRUE);
        vflags.KBIN_VASP_FORCE_OPTION_VOLUME.clear();

        // [OBSOLETE]  if(aurostd::substring2bool(AflowIn,_STROPT_+"VOLUME=",TRUE)) vflags.KBIN_VASP_FORCE_OPTION_VOLUME.push("EQUAL_EQUAL");
        // [OBSOLETE]  if(aurostd::substring2bool(AflowIn,_STROPT_+"VOLUME+=",TRUE)) vflags.KBIN_VASP_FORCE_OPTION_VOLUME.push("PLUS_EQUAL");
        // [OBSOLETE]  if(aurostd::substring2bool(AflowIn,_STROPT_+"VOLUME-=",TRUE)) vflags.KBIN_VASP_FORCE_OPTION_VOLUME.push("MINUS_EQUAL");
        // [OBSOLETE]  if(aurostd::substring2bool(AflowIn,_STROPT_+"VOLUME*=",TRUE)) vflags.KBIN_VASP_FORCE_OPTION_VOLUME.push("MULTIPLY_EQUAL");
        // [OBSOLETE]  if(aurostd::substring2bool(AflowIn,_STROPT_+"VOLUME/=",TRUE)) vflags.KBIN_VASP_FORCE_OPTION_VOLUME.push("DIVIDE_EQUAL");

        vflags.KBIN_VASP_FORCE_OPTION_VOLUME.args2addattachedscheme(vAflowIn,"EQUAL_EQUAL",_STROPT_+"VOLUME=","");
        vflags.KBIN_VASP_FORCE_OPTION_VOLUME.args2addattachedscheme(vAflowIn,"PLUS_EQUAL",_STROPT_+"VOLUME+=","");
        // [OBSOLETE] vflags.KBIN_VASP_FORCE_OPTION_VOLUME.args2addattachedscheme(vAflowIn,"MINUS_EQUAL",_STROPT_+"VOLUME-=","");
        vflags.KBIN_VASP_FORCE_OPTION_VOLUME.args2addattachedscheme(vAflowIn,"MULTIPLY_EQUAL",_STROPT_+"VOLUME*=","");
        // [OBSOLETE] vflags.KBIN_VASP_FORCE_OPTION_VOLUME.args2addattachedscheme(vAflowIn,"DIVIDE_EQUAL",_STROPT_+"VOLUME/=","");
        if(LDEBUG) cerr << "CORMAC STUFF  " << "vflags.KBIN_VASP_FORCE_OPTION_VOLUME.flag(\"EQUAL_EQUAL\")=" << vflags.KBIN_VASP_FORCE_OPTION_VOLUME.flag("EQUAL_EQUAL") << endl;
        if(LDEBUG) cerr << "CORMAC STUFF  " << "vflags.KBIN_VASP_FORCE_OPTION_VOLUME.getattachedscheme(\"EQUAL_EQUAL\")=" << vflags.KBIN_VASP_FORCE_OPTION_VOLUME.getattachedscheme("EQUAL_EQUAL") << endl;
        if(LDEBUG) cerr << "CORMAC STUFF  " << "vflags.KBIN_VASP_FORCE_OPTION_VOLUME.flag(\"PLUS_EQUAL\")=" << vflags.KBIN_VASP_FORCE_OPTION_VOLUME.flag("PLUS_EQUAL") << endl;
        if(LDEBUG) cerr << "CORMAC STUFF  " << "vflags.KBIN_VASP_FORCE_OPTION_VOLUME.getattachedscheme(\"PLUS_EQUAL\")=" << vflags.KBIN_VASP_FORCE_OPTION_VOLUME.getattachedscheme("PLUS_EQUAL") << endl;
        if(LDEBUG) cerr << "CORMAC STUFF  " << "vflags.KBIN_VASP_FORCE_OPTION_VOLUME.flag(\"MULTIPLY_EQUAL\")=" << vflags.KBIN_VASP_FORCE_OPTION_VOLUME.flag("MULTIPLY_EQUAL") << endl;
        if(LDEBUG) cerr << "CORMAC STUFF  " << "vflags.KBIN_VASP_FORCE_OPTION_VOLUME.getattachedscheme(\"MULTIPLY_EQUAL\")=" << vflags.KBIN_VASP_FORCE_OPTION_VOLUME.getattachedscheme("MULTIPLY_EQUAL") << endl;

        if(vflags.KBIN_VASP_FORCE_OPTION_VOLUME.xscheme!="") vflags.KBIN_VASP_FORCE_OPTION_VOLUME.isentry=TRUE;

        // [OBSOLETE] vflags.KBIN_VASP_POTCAR_FILE_KEYWORD              =    aurostd::substring2bool(AflowIn,"[VASP_POTCAR_FILE]");      
        // [OBSOLETE] vflags.KBIN_VASP_POTCAR_FILE_SYSTEM_AUTO          =    aurostd::substring2bool(AflowIn,"[VASP_POTCAR_FILE]SYSTEM_AUTO",TRUE);
        // [OBSOLETE] vflags.KBIN_VASP_POTCAR_FILE_PREFIX               =    aurostd::substring2bool(AflowIn,"[VASP_POTCAR_FILE]PREFIX=",TRUE);
        // [OBSOLETE] vflags.KBIN_VASP_POTCAR_FILE_SUFFIX               =    aurostd::substring2bool(AflowIn,"[VASP_POTCAR_FILE]SUFFIX=",TRUE);
        // [OBSOLETE] vflags.KBIN_VASP_POTCAR_FILE_FILE                 =    aurostd::substring2bool(AflowIn,"[VASP_POTCAR_FILE]FILE=",TRUE);
        // [OBSOLETE] vflags.KBIN_VASP_POTCAR_FILE_COMMAND              =    aurostd::substring2bool(AflowIn,"[VASP_POTCAR_FILE]COMMAND=",TRUE);

        if(aurostd::substring2bool(AflowIn,"[VASP_POTCAR_FILE]")) vflags.KBIN_VASP_POTCAR_FILE.push("KEYWORD");
        if(aurostd::substring2bool(AflowIn,"[VASP_POTCAR_FILE]SYSTEM_AUTO",TRUE)) vflags.KBIN_VASP_POTCAR_FILE.push("SYSTEM_AUTO");
        if(aurostd::substring2bool(AflowIn,"[VASP_POTCAR_FILE]PREFIX=",TRUE)) vflags.KBIN_VASP_POTCAR_FILE.push("PREFIX");
        if(aurostd::substring2bool(AflowIn,"[VASP_POTCAR_FILE]SUFFIX=",TRUE)) vflags.KBIN_VASP_POTCAR_FILE.push("SUFFIX");
        if(aurostd::substring2bool(AflowIn,"[VASP_POTCAR_FILE]FILE=",TRUE)) vflags.KBIN_VASP_POTCAR_FILE.push("FILE");
        if(aurostd::substring2bool(AflowIn,"[VASP_POTCAR_FILE]COMMAND=",TRUE)) vflags.KBIN_VASP_POTCAR_FILE.push("COMMAND");

        // [OBSOLETE] vflags.KBIN_VASP_POTCAR_MODE_EXPLICIT             =    aurostd::substring2bool(AflowIn,"[VASP_POTCAR_MODE_EXPLICIT]");
        // [OBSOLETE] vflags.KBIN_VASP_POTCAR_MODE_IMPLICIT             =    aurostd::substring2bool(AflowIn,"[VASP_POTCAR_MODE_IMPLICIT]");
        // [OBSOLETE] vflags.KBIN_VASP_POTCAR_MODE_EXTERNAL             =    aurostd::substring2bool(AflowIn,"[VASP_POTCAR_MODE_EXTERNAL]");

        if(aurostd::substring2bool(AflowIn,"[VASP_POTCAR_MODE_EXPLICIT]")) vflags.KBIN_VASP_POTCAR_MODE.push("EXPLICIT");
        if(aurostd::substring2bool(AflowIn,"[VASP_POTCAR_MODE_IMPLICIT]")) vflags.KBIN_VASP_POTCAR_MODE.push("IMPLICIT");
        if(aurostd::substring2bool(AflowIn,"[VASP_POTCAR_MODE_EXTERNAL]")) vflags.KBIN_VASP_POTCAR_MODE.push("EXTERNAL");

        // APL ENTRIES
        if(LDEBUG) cerr << "DEBUG: KBIN::VASP_Get_Vflags_from_AflowIN (APL)" << endl;

        // CO - START 170601
        // CO make backwards and forwards compatible with all possible workflows
        vflags.KBIN_VASP_KPOINTS_PHONONS_KPPRA.options2entry(AflowIn,"[AFLOW_APL]KPPRA=|[AFLOW_QHA]KPPRA=|[AFLOW_AAPL]KPPRA=|[AFLOW_PHONONS]KPPRA=",
                FALSE,vflags.KBIN_VASP_KPOINTS_PHONONS_KPPRA.xscheme); // scheme already loaded in aflow_xclasses.cpp is "1"
        // cerr << "vflags.KBIN_VASP_KPOINTS_PHONONS_KPPRA.isentry=" << vflags.KBIN_VASP_KPOINTS_PHONONS_KPPRA.isentry << endl << "vflags.KBIN_VASP_KPOINTS_PHONONS_KPPRA.xscheme=" << vflags.KBIN_VASP_KPOINTS_PHONONS_KPPRA.xscheme << endl;
        vflags.KBIN_VASP_KPOINTS_PHONONS_KSCHEME.options2entry(AflowIn,"[AFLOW_APL]KSCHEME=|[AFLOW_QHA]KSCHEME=|[AFLOW_AAPL]KSCHEME=|[AFLOW_PHONONS]KSCHEME=",
                FALSE,vflags.KBIN_VASP_KPOINTS_PHONONS_KSCHEME.xscheme); // scheme already loaded in aflow_xclasses.cpp is "DEFAULT_SCHEME"
        // cerr << "vflags.KBIN_VASP_KPOINTS_PHONONS_KSCHEME.isentry=" << vflags.KBIN_VASP_KPOINTS_PHONONS_KSCHEME.isentry << endl << "vflags.KBIN_VASP_KPOINTS_PHONONS_KSCHEME.xscheme=" << vflags.KBIN_VASP_KPOINTS_PHONONS_KSCHEME.xscheme << endl;


        vflags.KBIN_VASP_FORCE_OPTION_KPOINTS_PHONONS_PARITY.clear();
        if(aurostd::substring2bool(AflowIn,"[AFLOW_APL]KPOINTS=EVEN",TRUE) || 
                aurostd::substring2bool(AflowIn,"[AFLOW_QHA]KPOINTS=EVEN",TRUE) ||
                aurostd::substring2bool(AflowIn,"[AFLOW_AAPL]KPOINTS=EVEN",TRUE) ||
                aurostd::substring2bool(AflowIn,"[AFLOW_PHONONS]KPOINTS_EVEN",TRUE)) 
        {vflags.KBIN_VASP_FORCE_OPTION_KPOINTS_PHONONS_PARITY.push("EVEN");}
        if(aurostd::substring2bool(AflowIn,"[AFLOW_APL]KPOINTS=ODD",TRUE) || 
                aurostd::substring2bool(AflowIn,"[AFLOW_QHA]KPOINTS=ODD",TRUE) ||
                aurostd::substring2bool(AflowIn,"[AFLOW_AAPL]KPOINTS=ODD",TRUE) ||
                aurostd::substring2bool(AflowIn,"[AFLOW_PHONONS]KPOINTS_ODD",TRUE)) 
        {vflags.KBIN_VASP_FORCE_OPTION_KPOINTS_PHONONS_PARITY.push("ODD");}
        //  vflags.KBIN_VASP_FORCE_OPTION_KPOINTS_PHONONS_PARITY_EVEN=aurostd::substring2bool(AflowIn,"[AFLOW_PHONONS]KPOINTS=EVEN",TRUE)||aurostd::substring2bool(AflowIn,"[AFLOW_PHONONS]KPOINTS_EVEN",TRUE);
        // vflags.KBIN_VASP_FORCE_OPTION_KPOINTS_PHONONS_PARITY_ODD=aurostd::substring2bool(AflowIn,"[AFLOW_PHONONS]KPOINTS=ODD",TRUE)||aurostd::substring2bool(AflowIn,"[AFLOW_PHONONS]KPOINTS_ODD",TRUE);
        // CO - END 170601

        // FROZSL ENTRIES
        if(LDEBUG) cerr << "DEBUG: KBIN::VASP_Get_Vflags_from_AflowIN (FROZSL)" << endl;

        if(LDEBUG) cerr << "DEBUG: KBIN::VASP_Get_Vflags_from_AflowIN (STOP)" << endl;

        return vflags;
    }
} // namespace KBIN

// ******************************************************************************************************************************************************
// ******************************************************************************************************************************************************
namespace KBIN {
    bool VASP_ExtractNGF(string OUTCAR,int &NGXF,int &NGYF,int &NGZF);
} // namespace KBIN

//KESONG where VASP starts
namespace KBIN {
    bool VASP_Directory(ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags) { // AFLOW_FUNCTION_IMPLEMENTATION
        bool LDEBUG=(FALSE || XHOST.DEBUG);
        if(LDEBUG) cerr << "DEBUG: KBIN::VASP_Directory (BEGIN)" << endl;
        //  bool KBIN_MPI_LOCAL;KBIN_MPI_LOCAL=MPI;
        // bool KBIN_VASP_WRITE_KPOINTS;
        // string::size_type sub_size1,sub_size2;
        string subS,subS1,subS2;
        ostringstream aus;
        string::iterator pos;
        bool Krun=TRUE;

        deque<string> vext; aurostd::string2tokens(".bz2,.xz,.gz",vext,",");
        deque<string> vzip; aurostd::string2tokens("bzip2,xz,gzip",vzip,",");

        ifstream FileAFLOWIN;
        string FileNameAFLOWIN;
        string AflowIn;
        FileNameAFLOWIN=aflags.Directory+"/"+_AFLOWIN_;
        FileAFLOWIN.open(FileNameAFLOWIN.c_str(),std::ios::in);
        FileAFLOWIN.clear();FileAFLOWIN.seekg(0);
        AflowIn="";char c; while (FileAFLOWIN.get(c)) AflowIn+=c;  // READ _AFLOWIN_ and put into AflowIn
        FileAFLOWIN.clear();FileAFLOWIN.seekg(0);
        AflowIn=aurostd::RemoveComments(AflowIn); // NOW Clean AFLOWIN
        if(!FileAFLOWIN) {                                                                                      // ******* _AFLOWIN_ does not exist
            aus << "EEEEE  " << _AFLOWIN_ << " ABSENT   = " << Message(aflags,"user,host,time") << endl;
            aurostd::PrintMessageStream(aus,XHOST.QUIET);
            return FALSE;
        }
        aflags.QUIET=FALSE;
        _vflags vflags;
        vflags=KBIN::VASP_Get_Vflags_from_AflowIN(AflowIn,FileMESSAGE,aflags,kflags);

        // *********************************************************************************************************************
        // OPERATIONS related to PARTICULAR MACHINES ***************************************************************************

        if(LDEBUG) cerr << "[DEBUG] aflags.AFLOW_MACHINE_GLOBAL=" << aflags.AFLOW_MACHINE_GLOBAL << endl;
        if(LDEBUG) cerr << "[DEBUG] aflags.AFLOW_MACHINE_LOCAL=" << aflags.AFLOW_MACHINE_LOCAL << endl;

        // ***************************************************************************
        if(LDEBUG) cerr << "[DEBUG] vflags.KBIN_VASP_REPEAT.flag(\"REPEAT_BANDS\")=" << vflags.KBIN_VASP_REPEAT.flag("REPEAT_BANDS") << endl;
        if(LDEBUG) cerr << "[DEBUG] vflags.KBIN_VASP_REPEAT.flag(\"REPEAT_STATIC_BANDS\")=" << vflags.KBIN_VASP_REPEAT.flag("REPEAT_STATIC_BANDS") << endl;
        if(LDEBUG) cerr << "[DEBUG] vflags.KBIN_VASP_REPEAT.flag(\"REPEAT_DELSOL\")=" << vflags.KBIN_VASP_REPEAT.flag("REPEAT_DELSOL") << endl;

        // ***************************************************************************
        // Get the KBIN_BIN name
        aurostd::StringstreamClean(aus);
        aus << "00000  MESSAGE KBIN::VASP_Directory Running KBIN_BIN=\"" << kflags.KBIN_BIN << "\" " << Message(aflags,"user,host,time") << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        // ***************************************************************************
        // Some verbose
        if(kflags.KBIN_POCC) {aus << "00000  MESSAGE KBIN::VASP_Directory Running POCC_CALCULATION" << Message(aflags,"user,host,time") << endl;} // CO 180419 //POCC is special needs to run first because there is NO poscar defined yet
        else if(kflags.KBIN_PHONONS_CALCULATION_APL) aus << "00000  MESSAGE KBIN::VASP_Directory Running PHONONS_CALCULATION_APL" << Message(aflags,"user,host,time") << endl;
        else if(kflags.KBIN_PHONONS_CALCULATION_QHA) aus << "00000  MESSAGE KBIN::VASP_Directory Running PHONONS_CALCULATION_QHA" << Message(aflags,"user,host,time") << endl;   // CO 170601
        else if(kflags.KBIN_PHONONS_CALCULATION_AAPL) aus << "00000  MESSAGE KBIN::VASP_Directory Running PHONONS_CALCULATION_AAPL" << Message(aflags,"user,host,time") << endl; // CO 170601
        else if(kflags.KBIN_PHONONS_CALCULATION_AGL) aus << "00000  MESSAGE KBIN::VASP_Directory Running PHONONS_CALCULATION_AGL (Debye Model)" << Message(aflags,"user,host,time") << endl;
        else if(kflags.KBIN_PHONONS_CALCULATION_AEL) aus << "00000  MESSAGE KBIN::VASP_Directory Running PHONONS_CALCULATION_AEL (Elastic constants)" << Message(aflags,"user,host,time") << endl;
        else if(kflags.KBIN_PHONONS_CALCULATION_FROZSL) aus << "00000  MESSAGE KBIN::VASP_Directory Running PHONONS_CALCULATION_FROZSL" << Message(aflags,"user,host,time") << endl;
        else{
            if(vflags.KBIN_VASP_RUN.flag("GENERATE")) aus << "00000  MESSAGE KBIN::VASP_Directory Running RUN_GENERATE" << Message(aflags,"user,host,time") << endl;
            if(vflags.KBIN_VASP_RUN.flag("STATIC")) aus << "00000  MESSAGE KBIN::VASP_Directory Running RUN_STATIC" << Message(aflags,"user,host,time") << endl;
            if(vflags.KBIN_VASP_RUN.flag("KPOINTS")) aus << "00000  MESSAGE KBIN::VASP_Directory Running RUN_KPOINTS" << Message(aflags,"user,host,time") << endl;
            if(vflags.KBIN_VASP_RUN.flag("RELAX")) aus << "00000  MESSAGE KBIN::VASP_Directory Running RUN_RELAX" << Message(aflags,"user,host,time") << endl;
            if(vflags.KBIN_VASP_RUN.flag("RELAX_STATIC")) aus << "00000  MESSAGE KBIN::VASP_Directory Running RUN_RELAX_STATIC" << Message(aflags,"user,host,time") << endl;
            if(vflags.KBIN_VASP_RUN.flag("STATIC_BANDS")) aus << "00000  MESSAGE KBIN::VASP_Directory Running RUN_STATIC_BANDS" << Message(aflags,"user,host,time") << endl;
            if(vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_BANDS")) aus << "00000  MESSAGE KBIN::VASP_Directory Running RUN_RELAX_STATIC_BANDS" << Message(aflags,"user,host,time") << endl;
            if(vflags.KBIN_VASP_RUN.flag("DIELECTRIC_STATIC")) aus << "00000  MESSAGE KBIN::VASP_Directory Running RUN_DIELECTRIC_STATIC" << Message(aflags,"user,host,time") << endl;
            if(vflags.KBIN_VASP_RUN.flag("DIELECTRIC_DYNAMIC")) aus << "00000  MESSAGE KBIN::VASP_Directory Running RUN_DIELECTRIC_DYNAMIC" << Message(aflags,"user,host,time") << endl;
            if(vflags.KBIN_VASP_RUN.flag("DSCF")) aus << "00000  MESSAGE KBIN::VASP_Directory Running RUN_DSCF" << Message(aflags,"user,host,time") << endl;
        }
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET); //exit(0);
        // ***************************************************************************
        uint ntasks=0;
        ntasks=1; // default
        if(vflags.KBIN_VASP_POSCAR_MODE.flag("EXPLICIT_START_STOP_POINT")) {
            ntasks=vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING.size();
            aus << "00000  MESSAGE Loaded ntasks = " << ntasks << " - " << Message(aflags,"user,host,time") << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            for(uint i=0;i<vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING.size();i++) {
                aus << "00000  MESSAGE task " << i << "/" << ntasks << " in subdirectory " << vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING.at(i) << endl;
                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            }
        }
        // ***************************************************************************
        // start the loop !
        _aflags aflags_backup;aflags_backup=aflags;
        _kflags kflags_backup;kflags_backup=kflags;
        //  _vflags vflags_backup;vflags_backup=vflags;


        for(uint ixvasp=0;ixvasp<ntasks;ixvasp++) {  // LOOP ixvasp
            // declarations
            _xvasp xvasp;xvasp.clear();
            xvasp.POSCAR_index=ixvasp;
            aflags=aflags_backup;kflags=kflags_backup; // load it up
            // some verbose
            if(vflags.KBIN_VASP_POSCAR_MODE.flag("EXPLICIT_START_STOP_POINT")) {
                aus << "00000  MESSAGE START loop " << xvasp.POSCAR_index << "/" << vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING.size() << " - " << Message(aflags,"user,host,time") << endl;
                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            }
            if(LDEBUG) cerr << "KBIN::VASP_Directory: [1]" << xvasp.str << endl; 
            // ------------------------------------------
            // now start for each xvasp
            Krun=TRUE;  // guess everything is intelligent !
            xvasp.Directory=aflags.Directory;
            if(vflags.KBIN_VASP_POSCAR_MODE.flag("EXPLICIT_START_STOP_POINT")) {
                xvasp.Directory=aflags.Directory+"/"+KBIN_SUBDIRECTORIES+vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING.at(xvasp.POSCAR_index);
                aus << "00000  MESSAGE Taking loop directory = " << xvasp.Directory << " - " << Message(aflags,"user,host,time") << endl;
                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            }
            // check for directory KESONG CHEC THIS (if Krun=FALSE, everything stops).
            if(Krun && vflags.KBIN_VASP_POSCAR_MODE.flag("EXPLICIT_START_STOP_POINT")) {
                if(aurostd::FileExist(xvasp.Directory)) {
                    Krun=FALSE; // avoid rerunning
                    aus << "00000  MESSAGE Skipping loop directory = " << xvasp.Directory << " - " << Message(aflags,"user,host,time") << endl;
                    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                } else {
                    // before making it, check it again... NFS problem... check LOCK again
                    if(Krun && aurostd::FileExist(xvasp.Directory+"/"+_AFLOWLOCK_)) Krun=FALSE;    // to fight against NFS cache
                    if(Krun && aurostd::EFileExist(xvasp.Directory+"/"+_AFLOWLOCK_)) Krun=FALSE;     // to fight against NFS cache
                    if(Krun && aurostd::FileExist(xvasp.Directory+"/LLOCK")) Krun=FALSE;     // to fight against NFS cache
                    if(Krun && aurostd::EFileExist(xvasp.Directory+"/LLOCK")) Krun=FALSE;     // to fight against NFS cache
                    if(Krun) {
                        aurostd::DirectoryMake(xvasp.Directory);
                        aus << "00000  MESSAGE Creating loop directory = " << xvasp.Directory << " - " << Message(aflags,"user,host,time") << endl;
                        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                        aurostd::execute("echo \"NNNNN  KBIN LLOCK ASAP for NFS concurrent jobs (aflow"+string(AFLOW_VERSION)+")\" >> "+xvasp.Directory+"/LLOCK");
                    }
                }
            }


            if(Krun) {
                aflags.Directory=xvasp.Directory; // so we are set ! since there are plenty of routines with aflags.Directory inside
                aus << "00000  MESSAGE Performing loop directory = " << xvasp.Directory << " - " << Message(aflags,"user,host,time") << endl;
                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            }
            // ------------------------------------------
            // do the flags
            if(LDEBUG) cerr << "KBIN::VASP_Directory: [2]" << xvasp.str << endl;
            vflags.KBIN_VASP_INCAR_VERBOSE=TRUE; // ALWAYS
            if(vflags.KBIN_VASP_RUN.flag("STATIC")) vflags.KBIN_VASP_INCAR_VERBOSE=FALSE; // TURN OFF VERBOSITY
            if(vflags.KBIN_VASP_RUN.flag("STATIC_BANDS")) vflags.KBIN_VASP_INCAR_VERBOSE=FALSE; // TURN OFF VERBOSITY
            if(vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_BANDS")) vflags.KBIN_VASP_INCAR_VERBOSE=FALSE; // TURN OFF VERBOSITY
            if(vflags.KBIN_VASP_REPEAT.flag("REPEAT_BANDS")) vflags.KBIN_VASP_INCAR_VERBOSE=FALSE; // TURN OFF VERBOSITY
            if(vflags.KBIN_VASP_REPEAT.flag("REPEAT_STATIC_BANDS")) vflags.KBIN_VASP_INCAR_VERBOSE=FALSE; // TURN OFF VERBOSITY
            if(vflags.KBIN_VASP_REPEAT.flag("REPEAT_DELSOL")) vflags.KBIN_VASP_INCAR_VERBOSE=FALSE; // TURN OFF VERBOSITY

            // COREY, come back here at some point
            //think about taking pocc structure and reducing first if possible (and requested)
            //need to first convert to non-pocc structure (already programmed in aflow_pocc.cpp)
            //plug this in as xvasp.str (also create xvasp.str_pocc then)
            //reduce as requested and re-pocc the structure
            //think about if we need a separate flag for reducing pocc vs. reducing derivative structures
            // produce BEFORE NOMIX
            if(!(
                        kflags.KBIN_POCC ||
                        //kflags.KBIN_PHONONS_CALCULATION_APL ||    // CO 180515 - KEEP APL/QHA/AAPL to grab structure
                        //kflags.KBIN_PHONONS_CALCULATION_QHA ||    // CO 180515 - KEEP APL/QHA/AAPL to grab structure
                        //kflags.KBIN_PHONONS_CALCULATION_AAPL ||   // CO 180515 - KEEP APL/QHA/AAPL to grab structure
                        kflags.KBIN_PHONONS_CALCULATION_FROZSL ||   // CO 180503 - KEEP AEL/AGL stuff running per normal
                        //kflags.KBIN_PHONONS_CALCULATION_AGL ||    // CO 180503 - KEEP AEL/AGL stuff running per normal
                        //kflags.KBIN_PHONONS_CALCULATION_AEL ||    // CO 180503 - KEEP AEL/AGL stuff running per normal
                        FALSE)  //identity
              ) { // CO 180419, do NOT produce POSCAR for POCC //KSYANG
                if(Krun) Krun=(Krun && KBIN::VASP_Produce_INPUT(xvasp,AflowIn,FileMESSAGE,aflags,kflags,vflags));  //first produce, then modify, then, submit
                if(Krun) Krun=(Krun && KBIN::VASP_Modify_INPUT(xvasp,FileMESSAGE,aflags,kflags,vflags));
                if(Krun && kflags.KBIN_QSUB) Krun=(Krun && KBIN::QSUB_Extract(xvasp.xqsub,AflowIn,FileAFLOWIN,FileMESSAGE,aflags,kflags));
                if(Krun && kflags.KBIN_QSUB_MODE1) Krun=(Krun && KBIN::QSUB_Extract_Mode1(xvasp.xqsub,FileMESSAGE,aflags,kflags));
                if(Krun && kflags.KBIN_QSUB_MODE2) Krun=(Krun && KBIN::QSUB_Extract_Mode2(xvasp.xqsub,FileMESSAGE,aflags,kflags));
                if(Krun && kflags.KBIN_QSUB_MODE3) Krun=(Krun && KBIN::QSUB_Extract_Mode3(xvasp.xqsub,FileMESSAGE,aflags,kflags));
            }
            if(Krun && vflags.KBIN_VASP_FORCE_OPTION_SKIP_NOMIX.isentry) {
                string potentials=xvasp.POTCAR_POTENTIALS.str();
                if(!aurostd::substring2bool(aurostd::CleanFileName(xvasp.Directory+"/"),"/1/") &&
                        !aurostd::substring2bool(aurostd::CleanFileName(xvasp.Directory+"/"),"/2/") &&
                        !aurostd::substring2bool(aurostd::CleanFileName(xvasp.Directory+"/"),"/3/") &&
                        !aurostd::substring2bool(aurostd::CleanFileName(xvasp.Directory+"/"),"/58/") &&
                        !aurostd::substring2bool(aurostd::CleanFileName(xvasp.Directory+"/"),"/59/") &&
                        !aurostd::substring2bool(aurostd::CleanFileName(xvasp.Directory+"/"),"/60/") &&
                        !aurostd::substring2bool(aurostd::CleanFileName(xvasp.Directory+"/"),"/115/") &&
                        !aurostd::substring2bool(aurostd::CleanFileName(xvasp.Directory+"/"),"/116/") &&
                        !aurostd::substring2bool(aurostd::CleanFileName(xvasp.Directory+"/"),"/117/")
                  ) {
                    aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]SKIP_NOMIX (NEGLECT_NOMIX, NEGLECT_IMMISCIBLE) - " << Message(aflags,"user,host,time") << endl;
                    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                    //	cerr << "potentials=" << potentials << endl;
                    if(MiscibilityCheck(potentials)==MISCIBILITY_SYSTEM_NOMIX) {
                        aus << "00000  MESSAGE Skipping system: " << KBIN::VASP_PseudoPotential_CleanName(potentials) << " is known to be immiscible (aflow_nomix.cpp) - " << Message(aflags,"user,host,time") << endl;
                        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                        stringstream command("");
                        command << "cat " << xvasp.Directory << "/" << _AFLOWLOCK_ << " > " << xvasp.Directory << "/" << DEFAULT_AFLOW_IMMISCIBILITY_OUT << endl;
                        aurostd::execute(command);
                        Krun=FALSE;
                    }
                }
                if(MiscibilityCheck(potentials)==MISCIBILITY_SYSTEM_MISCIBLE) {
                    aus << "00000  MESSAGE Running system: " << KBIN::VASP_PseudoPotential_CleanName(potentials) << " is known to be miscible (aflow_nomix.cpp) - " << Message(aflags,"user,host,time") << endl;
                    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                    Krun=TRUE;
                }
                if(MiscibilityCheck(potentials)==MISCIBILITY_SYSTEM_UNKNOWN) {
                    aus << "00000  MESSAGE Running system: " << KBIN::VASP_PseudoPotential_CleanName(potentials) << " is unknown (aflow_nomix.cpp) - " << Message(aflags,"user,host,time") << endl;
                    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                    Krun=TRUE;
                }
            } 

            // produce AFTER NOMIX
            // if(Krun) Krun=(Krun && KBIN::VASP_Produce_INPUT(xvasp,AflowIn,FileAFLOWIN,FileMESSAGE,aflags,kflags,vflags));
            // if(Krun) Krun=(Krun && KBIN::VASP_Modify_INPUT(xvasp,FileMESSAGE,aflags,kflags,vflags));
            // if(Krun && kflags.KBIN_QSUB) Krun=(Krun && KBIN::QSUB_Extract(xvasp.xqsub,AflowIn,FileAFLOWIN,FileMESSAGE,aflags,kflags));
            // if(Krun && kflags.KBIN_QSUB_MODE1) Krun=(Krun && KBIN::QSUB_Extract_Mode1(xvasp.xqsub,FileMESSAGE,aflags,kflags));
            // if(Krun && kflags.KBIN_QSUB_MODE2) Krun=(Krun && KBIN::QSUB_Extract_Mode2(xvasp.xqsub,FileMESSAGE,aflags,kflags));
            // if(Krun && kflags.KBIN_QSUB_MODE3) Krun=(Krun && KBIN::QSUB_Extract_Mode3(xvasp.xqsub,FileMESSAGE,aflags,kflags));


            // ***************************************************************************
            // READY TO RUN
            if(Krun) {
                if(LDEBUG) cerr << "KBIN::VASP_Directory: [3]" << endl;
                if(LDEBUG) cerr << xvasp.str << endl;
                xvasp.NRELAX=0;
                bool Krun=true;
                ostringstream aus;
                bool PAWGGA2=FALSE;
                // ***************************************************************************
                // directory check
                ifstream DirectoryStream;
                DirectoryStream.open(xvasp.Directory.c_str(),std::ios::in);
                if(!DirectoryStream) {
                    //   aus << "EEEEE  DIRECTORY_NOT_FOUND = " << Message(aflags,"user,host,time") << endl;
                    aus << "XXXXX  MAKING DIRECTORY = " << xvasp.Directory << "  " << Message(aflags,"user,host,time") << endl;
                    aurostd::PrintMessageStream(aus,XHOST.QUIET); // return FALSE;
                    string str="mkdir "+xvasp.Directory;
                    system(str.c_str());
                }
                // some check
                if(!FileAFLOWIN) {                                                                                      // ******* _AFLOWIN_ does not exist
                    //    aus << "EEEEE  " << _AFLOWIN_ << " ABSENT   = " << Message(aflags,"user,host,time") << endl;
                    //    aurostd::PrintMessageStream(aus,XHOST.QUIET);
                    //    return FALSE;
                }
                // ***************************************************************************
                // DO THE SYMMETRY NEIGHBOURS CALCULATION
                //if(!kflags.KBIN_PHONONS_CALCULATION_FROZSL) {
                // DX

                if(!(
                            kflags.KBIN_POCC ||
                            kflags.KBIN_PHONONS_CALCULATION_APL ||
                            kflags.KBIN_PHONONS_CALCULATION_QHA || 
                            kflags.KBIN_PHONONS_CALCULATION_AAPL || 
                            kflags.KBIN_PHONONS_CALCULATION_FROZSL ||     // CO 180503 - KEEP AEL/AGL stuff running per normal
                            //kflags.KBIN_PHONONS_CALCULATION_AGL ||      // CO 180503 - KEEP AEL/AGL stuff running per normal
                            //kflags.KBIN_PHONONS_CALCULATION_AEL ||      // CO 180503 - KEEP AEL/AGL stuff running per normal
                            FALSE) || // CO 170601 //identity
                        aflags.KBIN_GEN_SYMMETRY_OF_AFLOWIN
                  ) {  // CO, do internally
                    // DX
                    if(Krun) Krun=KBIN_StepSymmetryPerform(xvasp.str,AflowIn,FileMESSAGE,aflags,kflags,TRUE,cout); // DO THE SYMMETRY CALCULATION
                    if(Krun) Krun=StepNeighboursPerform(xvasp.str,AflowIn,FileMESSAGE,aflags,kflags); // DO THE NEIGHBOURS CALCULATION
                    // DX
                    //cerr << "KBIN GEN SYMMETRY OF AFLOWIN: " << aflags.KBIN_GEN_SYMMETRY_OF_AFLOWIN << endl;
                    if(aflags.KBIN_GEN_SYMMETRY_OF_AFLOWIN){
                        return Krun;
                    }
                    // DX
                }
                // VASP VASP WRITE
                //   if(Krun) Krun=(Krun && KBIN::VASP_Write_INPUT(xvasp,vflags));
                // ***************************************************************************
                // VASP INPUT FILES ARE DONE, NOW WE CAN USE OR MODYFYING THEM
                if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry) {
                    aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]NOTUNE, no tuning xCARs - ";
                    aus << xvasp.Directory << " - K=[" << xvasp.str.kpoints_k1 << " " << xvasp.str.kpoints_k2 << " " << xvasp.str.kpoints_k3 << " " << xvasp.str.kpoints_kmax << "] - ";
                    aus << XHOST.hostname << " - " << aflow_get_time_string() << endl;
                    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                }	  
                // ***************************************************************************
                // VASP HOW TO RUN ??
                // GENERATE ONLY -------------------------------------------------------------
                if(vflags.KBIN_VASP_RUN.flag("GENERATE")) {
                    KBIN::VASP_Write_INPUT(xvasp,vflags); // VASP VASP WRITE
                    aus << "00000  MESSAGE VASP generation files ONLY " << Message(aflags,"user,host,time") << endl;
                    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                    Krun=FALSE;
                    xvasp.NRELAX=0;
                } else {
                    // CO 180420 - added if/else-if for workflows that need to PRECEDE relax/static/etc.
                    // RUN SOMETHING
                    if(kflags.KBIN_POCC) {  // RUN POCC ------------------------  // CO 180419 //POCC is special, run as priority
                        aus << "00000  MESSAGE PERFORMING POCC_CALCULATION " << Message(aflags,"user,host,time") << endl;
                        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                        xvasp.NRELAX=-3;}
                    else if(kflags.KBIN_PHONONS_CALCULATION_APL) {  // RUN PHONONS APL ------------------------
                        aus << "00000  MESSAGE PERFORMING PHONONS_CALCULATION_APL " << Message(aflags,"user,host,time") << endl;
                        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                        xvasp.NRELAX=-3;}
                    // CO - START 170601
                    else if(kflags.KBIN_PHONONS_CALCULATION_QHA) {  // RUN PHONONS QHA ------------------------
                        aus << "00000  MESSAGE PERFORMING PHONONS_CALCULATION_QHA " << Message(aflags,"user,host,time") << endl;
                        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                        xvasp.NRELAX=-3;}
                    else if(kflags.KBIN_PHONONS_CALCULATION_AAPL) {  // RUN PHONONS AAPL ------------------------
                        aus << "00000  MESSAGE PERFORMING PHONONS_CALCULATION_AAPL " << Message(aflags,"user,host,time") << endl;
                        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                        xvasp.NRELAX=-3;}
                    // CO - END 170601
                    else if(kflags.KBIN_PHONONS_CALCULATION_AGL) {  // RUN PHONONS AGL ------------------------
                        aus << "00000  MESSAGE PERFORMING PHONONS_CALCULATION_AGL (Debye Model) " << Message(aflags,"user,host,time") << endl;
                        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                        xvasp.NRELAX=-3;}
                    else if(kflags.KBIN_PHONONS_CALCULATION_AEL) {  // RUN PHONONS AEL ------------------------
                        aus << "00000  MESSAGE PERFORMING PHONONS_CALCULATION_AEL (Elastic constants) " << Message(aflags,"user,host,time") << endl;
                        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                        xvasp.NRELAX=-3;}
                    else if(kflags.KBIN_PHONONS_CALCULATION_FROZSL) {  // RUN PHONONS FROZSL ------------------------
                        aus << "00000  MESSAGE PERFORMING PHONONS_CALCULATION_FROZSL " << Message(aflags,"user,host,time") << endl;
                        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                        xvasp.NRELAX=-3;}
                    else{
                        if(vflags.KBIN_VASP_RUN.flag("STATIC")) {  // RUN STATIC ------------------------
                            aus << "00000  MESSAGE Performing Static RUN " << Message(aflags,"user,host,time") << endl;
                            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                            xvasp.NRELAX=-1;
                        }
                        if(vflags.KBIN_VASP_RUN.flag("KPOINTS")) {  // RUN KPOINTS ------------------------
                            aus << "00000  MESSAGE Running KPOINTS swap " << Message(aflags,"user,host,time") << endl;
                            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                            xvasp.NRELAX=-2;
                        }
                        if(vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_BANDS")) {  // RUN RELAX_STATIC_BANDS ------------------------
                            xvasp.NRELAX=vflags.KBIN_VASP_RUN_NRELAX; //  aurostd::substring2utype<int>(AflowIn,"[VASP_RUN_RELAX_STATIC_BANDS=");
                            if(xvasp.NRELAX<0)  {
                                aus << "EEEEE  No relaxation to run or nrelax<0 [nrelax=" << xvasp.NRELAX << "]  " << Message(aflags,"user,host,time") << endl;
                                aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);    
                                Krun=FALSE;	   //	  FileINPUT.clear();FileINPUT.close();FileMESSAGE.clear();FileMESSAGE.close();
                                xvasp.NRELAX=0;
                            }
                            {
                                aus << "00000  MESSAGE RELAX_STATIC_BANDS Running [nrelax=" << xvasp.NRELAX << "]  " << Message(aflags,"user,host,time") << endl;
                                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
                            }
                        }
                        if(vflags.KBIN_VASP_RUN.flag("RELAX_STATIC")) {  // RUN RELAX_STATIC ------------------------
                            xvasp.NRELAX=vflags.KBIN_VASP_RUN_NRELAX; //  aurostd::substring2utype<int>(AflowIn,"[VASP_RUN_RELAX_STATIC=");
                            if(xvasp.NRELAX<0)  {
                                aus << "EEEEE  No relaxation to run or nrelax<0 [nrelax=" << xvasp.NRELAX << "]  " << Message(aflags,"user,host,time") << endl;
                                aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);    
                                Krun=FALSE;	   //	  FileINPUT.clear();FileINPUT.close();FileMESSAGE.clear();FileMESSAGE.close();
                                xvasp.NRELAX=0;
                            }
                            aus << "00000  MESSAGE RELAX_STATIC Running [nrelax=" << xvasp.NRELAX << "]  " << Message(aflags,"user,host,time") << endl;
                            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
                        }
                        if(vflags.KBIN_VASP_RUN.flag("STATIC_BANDS")) { // RUN STATIC_BANDS ------------------------
                            xvasp.NRELAX=-1;	
                            aus << "00000  MESSAGE STATIC_BANDS Running  " << Message(aflags,"user,host,time") << endl;
                            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
                        }
                        if(vflags.KBIN_VASP_RUN.flag("RELAX")) { // RUN RELAX ------------------------
                            if(!(aurostd::substring2bool(AflowIn,"[VASP_RUN_RELAX=") || aurostd::substring2bool(AflowIn,"[VASP_RUN]RELAX="))) {
                                xvasp.NRELAX=2;
                                aus << "00000  MESSAGE Running DEFAULT [nrelax=" << xvasp.NRELAX << "]  " << Message(aflags,"user,host,time") << endl;
                                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                            } else  { 	  
                                xvasp.NRELAX=vflags.KBIN_VASP_RUN_NRELAX; //  aurostd::substring2utype<int>(AflowIn,"[VASP_RUN_RELAX=");
                            }
                            if(xvasp.NRELAX==0 || xvasp.NRELAX<0)  {
                                aus << "EEEEE  No relaxation to run or nrelax<0 [nrelax=" << xvasp.NRELAX << "]  " << Message(aflags,"user,host,time") << endl;
                                aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);    
                                Krun=FALSE;	   //	  FileINPUT.clear();FileINPUT.close();FileMESSAGE.clear();FileMESSAGE.close();
                                xvasp.NRELAX=0;
                            }
                            aus << "00000  MESSAGE RELAX Running [nrelax=" << xvasp.NRELAX << "]  " << Message(aflags,"user,host,time") << endl;
                            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
                        }
                        if(vflags.KBIN_VASP_REPEAT.flag("REPEAT_BANDS")) { // RUN REPEAT_BANDS ------------------------
                            xvasp.NRELAX=-1;	
                            aus << "00000  MESSAGE REPEAT_BANDS Running  " << Message(aflags,"user,host,time") << endl;
                            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
                        }
                        if(vflags.KBIN_VASP_REPEAT.flag("REPEAT_STATIC_BANDS")) { // RUN REPEAT_STATIC_BANDS ------------------------
                            xvasp.NRELAX=-1;	
                            aus << "00000  MESSAGE REPEAT_STATIC_BANDS Running  " << Message(aflags,"user,host,time") << endl;
                            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
                        }
                        if(vflags.KBIN_VASP_REPEAT.flag("REPEAT_DELSOL")) { // RUN REPEAT_DELSOL ------------------------
                            xvasp.NRELAX=-1;	
                            aus << "00000  MESSAGE REPEAT_DELSOL Running  " << Message(aflags,"user,host,time") << endl;
                            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
                        }
                    }

                    // ***************************************************************************
                    // READY TO RUN
                    if(Krun) {   // survived all troubles
                        // ***************************************************************************
                        // START
                        if(LDEBUG) cerr << "KBIN::VASP_Directory: [4]" << xvasp.str << endl;
                        // ***************************************************************************
                        // FIX BLANC SPECIES
                        if(xvasp.str.species.size()>0) {
                            if(xvasp.str.species.at(0)=="") {
                                pflow::fixEmptyAtomNames(xvasp.str);  //corey moved to pflow
                            }
                        }
                        // ***************************************************************************
                        // PRESCRIPT
                        if(kflags.AFLOW_MODE_PRESCRIPT_EXPLICIT || kflags.AFLOW_MODE_PRESCRIPT_EXPLICIT_START_STOP)
                            KBIN::RUN_DirectoryScript(aflags,DEFAULT_AFLOW_PRESCRIPT_COMMAND,DEFAULT_AFLOW_PRESCRIPT_OUT);
                        // ***************************************************************************
                        // CO 180419 - POCC always comes first (NO POSCAR), need to convert PARTCAR -> POSCARs
                        //other workflows follow, all of these precede relaxation/static/etc.
                        if(kflags.KBIN_POCC){KBIN::VASP_RunPOCC(xvasp,AflowIn,aflags,kflags,vflags,FileMESSAGE);} // CO 180419
                        else if(kflags.KBIN_PHONONS_CALCULATION_APL || kflags.KBIN_PHONONS_CALCULATION_QHA || kflags.KBIN_PHONONS_CALCULATION_AAPL) {KBIN::VASP_RunPhonons_APL(xvasp,AflowIn,aflags,kflags,vflags,FileMESSAGE);} // PHONONIC PHONONIC PHONONIC // CO 170601
                        else if(kflags.KBIN_PHONONS_CALCULATION_AGL==TRUE) {KBIN::VASP_RunPhonons_AGL(xvasp,AflowIn,aflags,kflags,vflags,FileMESSAGE);} // here CORMAC WILL POOP an IF=>AGL
                        else if(kflags.KBIN_PHONONS_CALCULATION_AEL==TRUE) {KBIN::VASP_RunPhonons_AEL(xvasp,AflowIn,aflags,kflags,vflags,FileMESSAGE);}
                        else if(kflags.KBIN_PHONONS_CALCULATION_FROZSL) {KBIN::VASP_RunPhonons_FROZSL(xvasp,AflowIn,aflags,kflags,vflags,FileMESSAGE);}
                        else{
                            // --------------------------------------------------------------------------------------------------------------------
                            // STATIC STATIC STATIC
                            if(vflags.KBIN_VASP_RUN.flag("STATIC")) {    // xvasp.RELAX=-1
                                xvasp.aopts.flag("FLAG::POSCAR_PRESERVED",TRUE); // in case of errors it is not lost bur recycled
                                KBIN::VASP_Write_INPUT(xvasp,vflags); // VASP VASP WRITE
                                aus << 11111 << "  STATIC - " <<  xvasp.Directory << " - K=[" << xvasp.str.kpoints_k1 << " " << xvasp.str.kpoints_k2 << " " << xvasp.str.kpoints_k3 << "]" << " - " << kflags.KBIN_BIN << " - " << Message("user,host,time") << endl;
                                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                                Krun=KBIN::VASP_Run(xvasp,aflags,kflags,vflags,FileMESSAGE);
                                if(!Krun) {KBIN::VASP_Error(xvasp,FileMESSAGE,"EEEEE  runtime error [STATIC] ");return Krun;}
                                bool qmwrite=TRUE;
                                KBIN::VASP_Backup(xvasp,qmwrite,string("static"));
                            }		
                            // --------------------------------------------------------------------------------------------------------------------
                            // RELAX RELAX RELAX
                            if(vflags.KBIN_VASP_RUN.flag("RELAX")) {    // xvasp.RELAX>0
                                KBIN::VASP_Write_INPUT(xvasp,vflags); // VASP VASP WRITE
                                if(PAWGGA2) {  // WAS A BUG IN PAW MAYBE IT IS FIXED
                                    // STEP 1
                                    aus << "11111  RELAXATION - " <<  xvasp.Directory << " - K=[" << xvasp.str.kpoints_k1 << " " << xvasp.str.kpoints_k2 << " " << xvasp.str.kpoints_k3 << "]" << " - " << kflags.KBIN_BIN << " - " << Message("user,host,time") << endl;
                                    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                                    Krun=KBIN::VASP_Run(xvasp,aflags,kflags,vflags,"relax2paw_gga",TRUE,FileMESSAGE);
                                    if(!Krun) {KBIN::VASP_Error(xvasp,FileMESSAGE,"EEEEE  runtime error [PAWGGA2 REL]");return Krun;}
                                    aus << "22222  END        - " <<  xvasp.Directory << " - K=[" << xvasp.str.kpoints_k1 << " " << xvasp.str.kpoints_k2 << " " << xvasp.str.kpoints_k3 << "]" << " - " << kflags.KBIN_BIN << " - " << Message("user,host,time") << endl;
                                    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
                                } else {
                                    if(xvasp.NRELAX==0) {
                                        cerr << "STATIC RUN FIX INCAR: should not be here" << endl;
                                        return FALSE;
                                    } else { // DYNAMIC RUN
                                        for(xvasp.NRELAXING=1;xvasp.NRELAXING<=xvasp.NRELAX;xvasp.NRELAXING++) {
                                            aus << 11111*xvasp.NRELAXING << "  RELAXATION - " <<  xvasp.Directory << " - K=[" << xvasp.str.kpoints_k1 << " " << xvasp.str.kpoints_k2 << " " << xvasp.str.kpoints_k3 << "]" << " - " << kflags.KBIN_BIN << " - " << Message("user,host,time") << endl;
                                            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
                                            if(vflags.KBIN_VASP_FORCE_OPTION_LDAU_ADIABATIC.content_int>0) KBIN::XVASP_INCAR_LDAU_ADIABATIC(xvasp,xvasp.NRELAXING); // ADIABATIC
                                            if(xvasp.NRELAXING<xvasp.NRELAX)  {
                                                Krun=KBIN::VASP_Run(xvasp,aflags,kflags,vflags,"relax"+aurostd::utype2string(xvasp.NRELAXING),"relax"+aurostd::utype2string(xvasp.NRELAXING),TRUE,FileMESSAGE);
                                                if(!Krun) {KBIN::VASP_Error(xvasp,FileMESSAGE,"EEEEE  runtime error [RELAXATION<]");return Krun;}
                                                KBIN::XVASP_INCAR_SPIN_REMOVE_RELAX(xvasp,aflags,vflags,xvasp.NRELAXING,FileMESSAGE);         // check if it is the case of turning off spin
                                                KBIN::XVASP_KPOINTS_IBZKPT_UPDATE(xvasp,aflags,vflags,xvasp.NRELAXING,FileMESSAGE);           // check if it is the case of updating IBZKPT
                                            }
                                            if(xvasp.NRELAXING==xvasp.NRELAX) {
                                                Krun=KBIN::VASP_Run(xvasp,aflags,kflags,vflags,"relax"+aurostd::utype2string(xvasp.NRELAXING),TRUE,FileMESSAGE);
                                                if(!Krun) {KBIN::VASP_Error(xvasp,FileMESSAGE,"EEEEE  runtime error [RELAXATION=]");return Krun;}
                                            }
                                        }
                                        xvasp.NRELAXING=xvasp.NRELAX;
                                        xvasp.NRELAXING++;
                                        aus << 11111*xvasp.NRELAXING << "  END        - " <<  xvasp.Directory << " - K=[" << xvasp.str.kpoints_k1 << " " << xvasp.str.kpoints_k2 << " " << xvasp.str.kpoints_k3 << "]" << " - " << kflags.KBIN_BIN << " - " << Message("user,host,time") << endl;
                                        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET); 
                                    }
                                }
                            }
                            // --------------------------------------------------------------------------------------------------------------------
                            // --------------------------------------------------------------------------------------------------------------------
                            // --------------------------------------------------------------------------------------------------------------------
                            // RELAX_STATIC_BANDS RELAX_STATIC_BANDS RELAX_STATIC_BANDS REPEAT_STATIC_BANDS REPEAT_BANDS
                            // STATIC_BANDS STATIC_BANDS STATIC_BANDS 	
                            // RELAX_STATIC RELAX_STATIC RELAX_STATIC
                            // REPEAT_STATIC_BANDS REPEAT_STATIC_BANDS
                            // REPEAT_BANDS REPEAT_BANDS REPEAT_BANDS
                            if(vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_BANDS") || vflags.KBIN_VASP_RUN.flag("STATIC_BANDS") || vflags.KBIN_VASP_RUN.flag("RELAX_STATIC") || vflags.KBIN_VASP_REPEAT.flag("REPEAT_STATIC_BANDS") || vflags.KBIN_VASP_REPEAT.flag("REPEAT_BANDS")) {    // xvasp.RELAX>0
                                vector<double> xvasp_spin_evolution;
                                xmatrix<double> rlattice(xvasp.str.lattice);

                                string STRING_TO_SHOW="";
                                if(vflags.KBIN_VASP_RUN.flag("STATIC")) STRING_TO_SHOW="STATIC";
                                if(vflags.KBIN_VASP_RUN.flag("RELAX_STATIC")) STRING_TO_SHOW="RELAX_STATIC";
                                if(vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_BANDS")) STRING_TO_SHOW="RELAX_STATIC_BANDS";
                                if(vflags.KBIN_VASP_RUN.flag("STATIC_BANDS")) STRING_TO_SHOW="STATIC_BANDS";
                                if(vflags.KBIN_VASP_REPEAT.flag("REPEAT_STATIC_BANDS")) STRING_TO_SHOW="REPEAT_STATIC_BANDS";
                                if(vflags.KBIN_VASP_REPEAT.flag("REPEAT_BANDS")) STRING_TO_SHOW="REPEAT_BANDS";
                                aus << "00000  MESSAGE MODE= (" << STRING_TO_SHOW << ") - " << xvasp.Directory << " - " << Message("user,host,time") << endl;
                                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    

                                xvasp.aopts.flag("FLAG::POSCAR_PRESERVED",TRUE); // in case of errors it is not lost bur recycled
                                xvasp.aopts.flag("FLAG::CHGCAR_PRESERVED",TRUE); // in case of errors it is not lost bur recycled

                                if(vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_BANDS") || vflags.KBIN_VASP_RUN.flag("STATIC_BANDS") || vflags.KBIN_VASP_RUN.flag("RELAX_STATIC") || vflags.KBIN_VASP_RUN.flag("STATIC")) {
                                    // DO THE RELAX PART (IF ANY)
                                    KBIN::VASP_Write_INPUT(xvasp,vflags); // VASP VASP WRITE
                                    if(vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_BANDS") || vflags.KBIN_VASP_RUN.flag("RELAX_STATIC")) {
                                        for(xvasp.NRELAXING=1;xvasp.NRELAXING<=xvasp.NRELAX;xvasp.NRELAXING++) {
                                            aus << 11111*xvasp.NRELAXING << "  RELAXATION (" << STRING_TO_SHOW << ") - " <<  xvasp.Directory << " - K=[" << xvasp.str.kpoints_k1 << " " << xvasp.str.kpoints_k2 << " " << xvasp.str.kpoints_k3 << "]" << " - " << kflags.KBIN_BIN << " - " << Message("user,host,time") << endl;
                                            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                                            if(vflags.KBIN_VASP_FORCE_OPTION_LDAU_ADIABATIC.content_int>0) KBIN::XVASP_INCAR_LDAU_ADIABATIC(xvasp,xvasp.NRELAXING); // ADIABATIC
                                            if(xvasp.NRELAXING<xvasp.NRELAX)  {
                                                Krun=KBIN::VASP_Run(xvasp,aflags,kflags,vflags,"relax"+aurostd::utype2string(xvasp.NRELAXING),"relax"+aurostd::utype2string(xvasp.NRELAXING),TRUE,FileMESSAGE);
                                                if(!Krun) {KBIN::VASP_Error(xvasp,FileMESSAGE,"EEEEE  runtime error [RELAX_STATIC_BANDS RELAXATION<]");return Krun;}
                                            }
                                            if(xvasp.NRELAXING==xvasp.NRELAX) {
                                                Krun=KBIN::VASP_Run(xvasp,aflags,kflags,vflags,"relax"+aurostd::utype2string(xvasp.NRELAXING),TRUE,FileMESSAGE);
                                                if(!Krun) {KBIN::VASP_Error(xvasp,FileMESSAGE,"EEEEE  runtime error [RELAX_STATIC_BANDS RELAXATION=]");return Krun;}
                                            }
                                            xvasp_spin_evolution.push_back(xvasp.str.qm_mag_atom); // keep track of spins
                                            aus << "00000  MESSAGE RESULT SPIN=" << xvasp_spin_evolution.at(xvasp_spin_evolution.size()-1) << " - " << Message(aflags,"user,host,time") << endl;
                                            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                                            if(xvasp.NRELAXING<xvasp.NRELAX) KBIN::XVASP_INCAR_SPIN_REMOVE_RELAX(xvasp,aflags,vflags,xvasp.NRELAXING,FileMESSAGE); 	// check if it is the case of turning off spin
                                        }
                                        if(xvasp.NRELAX>0) KBIN::VASP_Recycle(xvasp,"relax"+aurostd::utype2string(xvasp.NRELAX));  // bring back the stuff
                                        if(xvasp.NRELAX==2) KBIN::XVASP_INCAR_SPIN_REMOVE_RELAX(xvasp,aflags,vflags,xvasp.NRELAX,FileMESSAGE); 	// check if it is the case of turning off spin
                                    }
                                    if(vflags.KBIN_VASP_RUN.flag("RELAX_STATIC")) {
                                        aus << "00000  NO RELAXATION IN (" << STRING_TO_SHOW << ") - " << xvasp.Directory << " - " << Message("user,host,time") << endl;
                                        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
                                    }
                                    if(vflags.KBIN_VASP_RUN.flag("STATIC")) {
                                        aus << "00000  NO RELAXATION IN (" << STRING_TO_SHOW << ") - " << xvasp.Directory << " - " << Message("user,host,time") << endl;
                                        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                                        xvasp.NRELAX=0;
                                    }
                                    xvasp.NRELAXING=xvasp.NRELAX;
                                } // vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_BANDS") || vflags.KBIN_VASP_RUN.flag("STATIC_BANDS") || vflags.KBIN_VASP_RUN.flag("RELAX_STATIC")
                                // REPEAT_STATIC_BANDS PART ----------------------------------------------------------------------------
                                if(vflags.KBIN_VASP_REPEAT.flag("REPEAT_STATIC_BANDS")) {
                                    // LOAD FORMER LOCK
                                    if(aurostd::FileExist(xvasp.Directory+string("/REPEAT_STATIC_BANDS"))) {
                                        stringstream lock_recycled;
                                        aurostd::file2stringstream(xvasp.Directory+"/REPEAT_STATIC_BANDS",lock_recycled);
                                        aus << "XXXXX ---------------------------------------------------------------------------------------------- " << endl;
                                        aus << "XXXXX FORMER LOCK BEGIN, recycled (" << STRING_TO_SHOW << ") - " << xvasp.Directory << " - " << Message("user,host,time") << endl;
                                        //	aus << lock_recycled.str();
                                        aus << "XXXXX FORMER LOCK END, recycled (" << STRING_TO_SHOW << ") - " << xvasp.Directory << " - " << Message("user,host,time") << endl;
                                        aus << "XXXXX ---------------------------------------------------------------------------------------------- " << endl;
                                        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                                    }
                                    // UNZIP EVERYTHING
                                    for(uint iext=0;iext<vext.size();iext++) { 
                                        // aurostd::execute("cd "+xvasp.Directory+" && "+"bzip2 -dfq *bz2 "); // ORIGINAL
                                        aus << "DEBUG - KBIN::VASP_Directory: EXT POINT [1] " << endl; aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);		  
                                        aurostd::execute(vzip.at(iext)+" -dfq `find "+xvasp.Directory+" -name \"*"+vext.at(iext)+"\"`");
                                    }		
                                    if(aurostd::FileExist(xvasp.Directory+string("/POSCAR.relax2"))) {
                                        KBIN::VASP_Recycle(xvasp,"relax2");
                                    } else {
                                        if(aurostd::FileExist(xvasp.Directory+string("/POSCAR.relax1"))) {
                                            KBIN::VASP_Recycle(xvasp,"relax1");
                                        } else {
                                            cerr << "REPEAT_STATIC_BANDS: ERROR: RELAX2 or RELAX1 must be present.. error" << endl;
                                            exit(0);
                                        }
                                    }
                                    // clean up worthless stuff
                                    aurostd::execute("cd "+xvasp.Directory+" && rm -f *bands* *static*");
                                } // vflags.KBIN_VASP_REPEAT.flag("REPEAT_STATIC_BANDS")
                                // REPEAT_BANDS PART ----------------------------------------------------------------------------
                                if(vflags.KBIN_VASP_REPEAT.flag("REPEAT_BANDS")) {
                                    // LOAD FORMER LOCK
                                    if(aurostd::FileExist(xvasp.Directory+string("/REPEAT_BANDS"))) {
                                        stringstream lock_recycled;
                                        aurostd::file2stringstream(xvasp.Directory+"/REPEAT_BANDS",lock_recycled);
                                        aus << "XXXXX ---------------------------------------------------------------------------------------------- " << endl;
                                        aus << "XXXXX FORMER LOCK BEGIN, recycled (" << STRING_TO_SHOW << ") - " << xvasp.Directory << " - " << Message("user,host,time") << endl;
                                        //	aus << lock_recycled.str();
                                        aus << "XXXXX FORMER LOCK END, recycled (" << STRING_TO_SHOW << ") - " << xvasp.Directory << " - " << Message("user,host,time") << endl;
                                        aus << "XXXXX ---------------------------------------------------------------------------------------------- " << endl;
                                        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                                    }
                                    // UNZIP EVERYTHING
                                    for(uint iext=0;iext<vext.size();iext++) { 
                                        // aurostd::execute("cd "+xvasp.Directory+" && "+"bzip2 -dfq *bz2 "); // ORIGINAL
                                        aus << "DEBUG - KBIN::VASP_Directory: EXT POINT [2] " << endl; aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);		  
                                        aurostd::execute(vzip.at(iext)+" -dfq `find "+xvasp.Directory+" -name \"*"+vext.at(iext)+"\"`");
                                    }		

                                    if(aurostd::FileExist(xvasp.Directory+string("/POSCAR.static"))) {
                                        KBIN::VASP_Recycle(xvasp,"static");
                                    } else {
                                        cerr << "REPEAT_BANDS: ERROR: STATIC must be present.. error" << endl;
                                        exit(0);
                                    }
                                    // clean up worthless stuff
                                    aurostd::execute("cd "+xvasp.Directory+" && rm -f *bands* ");
                                } // vflags.KBIN_VASP_REPEAT.flag("REPEAT_BANDS")
                                // STATIC PART ----------------------------------------------------------------------------
                                // STATIC PART ----------------------------------------------------------------------------
                                // STATIC PART ----------------------------------------------------------------------------
                                // STATIC PART ----------------------------------------------------------------------------
                                // NOW DO THE STATIC PATCHING POSCAR
                                if(vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_BANDS") || vflags.KBIN_VASP_RUN.flag("STATIC_BANDS") || vflags.KBIN_VASP_RUN.flag("RELAX_STATIC") || vflags.KBIN_VASP_REPEAT.flag("REPEAT_STATIC_BANDS") || vflags.KBIN_VASP_RUN.flag("STATIC")) {
                                    aus << "00000  MESSAGE Patching POSCAR  " << Message(aflags,"user,host,time") << endl;
                                    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                                    vflags.KBIN_VASP_RUN.push("STATIC");        // force to suck them up from STATIC_KPPRA....
                                    // LOAD THE RELAXED STRUCTURE WHICH WILL BE USED FOR THE DIRECTIONS
                                    stringstream straus;
                                    aurostd::file2stringstream(xvasp.Directory+"/POSCAR",straus);
                                    xvasp.str=xstructure(straus,IOVASP_AUTO);
                                    xvasp.str.FixLattices();
                                    rlattice=xvasp.str.lattice; // in rlattice I`ve always the final structure
                                    bool STATIC_DEBUG=FALSE;//TRUE;
                                    // RECREATE CONVENTIONAL OR PRIMITIVE
                                    if(vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_BANDS") || vflags.KBIN_VASP_RUN.flag("STATIC_BANDS") || vflags.KBIN_VASP_REPEAT.flag("REPEAT_STATIC_BANDS")) {
                                        // shall we restandardize ?
                                    }
                                    // WAHYU_LOOK
                                    //    STATIC_DEBUG=TRUE;
                                    if(vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_BANDS") || vflags.KBIN_VASP_RUN.flag("STATIC_BANDS") || vflags.KBIN_VASP_REPEAT.flag("REPEAT_STATIC_BANDS")) {
                                        if(STATIC_DEBUG) {aus << "STATIC_DEBUG: " << endl;}
                                        if(STATIC_DEBUG) {aus << "STATIC_DEBUG: vflags.KBIN_VASP_RUN.flag(\"RELAX_STATIC_BANDS\")=" << vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_BANDS") << endl;}
                                        if(STATIC_DEBUG) {aus << "STATIC_DEBUG: vflags.KBIN_VASP_RUN.flag(\"STATIC_BANDS\")=" << vflags.KBIN_VASP_RUN.flag("STATIC_BANDS") << endl;}
                                        if(STATIC_DEBUG) {aus << "STATIC_DEBUG: vflags.KBIN_VASP_REPEAT.flag(\"REPEAT_STATIC_BANDS\")=" << vflags.KBIN_VASP_REPEAT.flag("REPEAT_STATIC_BANDS") << endl;}
                                        if(STATIC_DEBUG) {aus << "STATIC_DEBUG: vflags.KBIN_VASP_REPEAT.flag(\"REPEAT_BANDS\")=" << vflags.KBIN_VASP_REPEAT.flag("REPEAT_BANDS") << endl;}
                                        if(STATIC_DEBUG) {aus << "STATIC_DEBUG: vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE.isentry=" << vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE.isentry << endl;}
                                        if(STATIC_DEBUG) {aus << "STATIC_DEBUG: vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE.content_string=" << vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE.content_string << endl;}
                                        if(STATIC_DEBUG) {aus << "STATIC_DEBUG: vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE_AUTO_FLAG=" << vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE_AUTO_FLAG << endl;}
                                        if(STATIC_DEBUG) {aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
                                        string stringBZ="";
                                        bool foundBZ=FALSE;
                                        if(STATIC_DEBUG) {aus << "STATIC_DEBUG: " << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
                                        if(STATIC_DEBUG) {aus << "STATIC_DEBUG: foundBZ=" << foundBZ << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
                                        if(STATIC_DEBUG) {aus << "STATIC_DEBUG: " << xvasp.str << endl;}
                                        if(STATIC_DEBUG) {aus << "STATIC_DEBUG: " << rlattice << endl;}
                                        if(vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE_AUTO_FLAG==FALSE) {
                                            stringBZ=LATTICE::KPOINTS_Directions(vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE.content_string,rlattice,vflags.KBIN_VASP_KPOINTS_BANDS_GRID_VALUE,xvasp.str.iomode,foundBZ); // rlattice = updated structure
                                        } else {
                                            foundBZ=FALSE;
                                        }
                                        if(STATIC_DEBUG) {aus << "STATIC_DEBUG: " << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
                                        if(STATIC_DEBUG) {aus << "STATIC_DEBUG: " << stringBZ << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
                                        if(STATIC_DEBUG) {aus << "STATIC_DEBUG: foundBZ=" << foundBZ << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
                                        if(STATIC_DEBUG) {aus << "STATIC_DEBUG: " << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}

                                        // always recalculate standardization
                                        if(vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("PRESERVE")==FALSE) {
                                            aus << "00000  MESSAGE WARNING RECALCULATING STANDARD STRUCTURE" << " - " << Message(aflags,"user,host,time") << endl;
                                            aus << "00000  MESSAGE BEFORE: a,b,c,alpha,beta,gamma " << xvasp.str.a << "," << xvasp.str.b << "," << xvasp.str.c << "," << xvasp.str.alpha << "," << xvasp.str.beta << "," << xvasp.str.gamma << endl;
                                            aus << "00000  MESSAGE BEFORE: lattice: " << vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE.content_string << endl;
                                            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                                            // reshuffle the structure
                                            if(vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("STANDARD_CONVENTIONAL")) {xvasp.str.Standard_Conventional_UnitCellForm();}
                                            else {xvasp.str.Standard_Primitive_UnitCellForm();}
                                            xvasp.POSCAR.str(std::string());xvasp.POSCAR.clear();
                                            xvasp.POSCAR << xvasp.str;
                                            xvasp.aopts.flag("FLAG::XVASP_POSCAR_generated",TRUE);xvasp.aopts.flag("FLAG::XVASP_POSCAR_changed",TRUE);
                                            aurostd::stringstream2file(xvasp.POSCAR,string(xvasp.Directory+"/POSCAR"));
                                            xvasp.str.FixLattices();
                                            rlattice=xvasp.str.lattice; // in rlattice I`ve always the final structure
                                            // [OBSOLETE] vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE_VALUE=xvasp.str.bravais_lattice_variation_type;//wahyu mod
                                            // vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE_VALUE=xvasp.str.bravais_conventional_lattice_type;
                                            vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE.clear();vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE.push(xvasp.str.bravais_lattice_variation_type); //wahyu mod

                                            aus << "00000  MESSAGE AFTER: a,b,c,alpha,beta,gamma " << xvasp.str.a << "," << xvasp.str.b << "," << xvasp.str.c << "," << xvasp.str.alpha << "," << xvasp.str.beta << "," << xvasp.str.gamma << endl;
                                            aus << "00000  MESSAGE AFTER: lattice: " << vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE.content_string << endl;
                                            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                                            // 		  cerr << "vasp.str.bravais_lattice_variation_type=" << xvasp.str.bravais_lattice_variation_type << endl;
                                            // 		  cerr << "vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE.content_string=" << vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE.content_string << endl;		  exit(0);
                                        } else {
                                            // nothing
                                            aus << "00000  MESSAGE PRESERVING ORIGINAL STRUCTURE" << " - " << Message(aflags,"user,host,time") << endl;
                                            aus << "00000  MESSAGE ORIGINAL: a,b,c,alpha,beta,gamma " << xvasp.str.a << "," << xvasp.str.b << "," << xvasp.str.c << "," << xvasp.str.alpha << "," << xvasp.str.beta << "," << xvasp.str.gamma << endl;
                                            aus << "00000  MESSAGE ORIGINAL: lattice: " << vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE.content_string << endl;
                                            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                                            // 		  cerr << "vasp.str.bravais_lattice_variation_type=" << xvasp.str.bravais_lattice_variation_type << endl;
                                            // 		  cerr << "vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE.content_string=" << vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE.content_string << endl;		  exit(0);  
                                        }
                                        stringBZ=LATTICE::KPOINTS_Directions(vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE.content_string,rlattice,vflags.KBIN_VASP_KPOINTS_BANDS_GRID_VALUE,xvasp.str.iomode,foundBZ); // rlattice = updated structure
                                        if(STATIC_DEBUG) {aus << "STATIC_DEBUG: " << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
                                        if(STATIC_DEBUG) {aus << "STATIC_DEBUG: " << stringBZ << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
                                        if(STATIC_DEBUG) {aus << "STATIC_DEBUG: foundBZ=" << foundBZ << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
                                        if(STATIC_DEBUG) {aus << "STATIC_DEBUG: " << xvasp.str << endl;}
                                        if(STATIC_DEBUG) {aus << "STATIC_DEBUG: " << rlattice << endl;}
                                        if(STATIC_DEBUG) {aus << "STATIC_DEBUG: " << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
                                        if(foundBZ==FALSE) {
                                            aus << "EEEEE  ERROR (KBIN::VASP_Directory) Unrecoverable error, lattice not found " << " - " << Message(aflags,"user,host,time") << endl;
                                            aus << xvasp.str << endl;
                                            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                                            exit(0);
                                        }
                                    }
                                    // done with the fixing
                                    xvasp.str.FixLattices();
                                    rlattice=xvasp.str.lattice; // in rlattice I`ve always the final structure
                                    // NOW DO THE STATIC PATCHING KPOINTS
                                    aus << "00000  MESSAGE Patching KPOINTS  " << Message(aflags,"user,host,time") << endl;
                                    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                                    //
                                    // [OBSOLETE]	      KBIN::VASP_Produce_KPOINTS(xvasp,AflowIn,FileAFLOWIN,FileMESSAGE,aflags,kflags,vflags);
                                    KBIN::VASP_Produce_KPOINTS(xvasp,AflowIn,FileMESSAGE,aflags,kflags,vflags);
                                    KBIN::VASP_Modify_KPOINTS(xvasp,FileMESSAGE,aflags,vflags);
                                    aurostd::stringstream2file(xvasp.KPOINTS,string(xvasp.Directory+"/KPOINTS"));
                                    // NOW DO THE STATIC PATCHING INCAR
                                    aus << "00000  MESSAGE [" << STRING_TO_SHOW << "] Patching INCAR (static_patching) " << Message(aflags,"user,host,time") << endl;
                                    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
                                    // KBIN::VASP_Produce_INCAR(xvasp,AflowIn,FileAFLOWIN,FileMESSAGE,aflags,kflags,vflags); // BETTER than produce, SHOULD reread it
                                    KBIN::VASP_Reread_INCAR(xvasp,FileMESSAGE,aflags); // REREAD IT
                                    // KBIN::VASP_Modify_INCAR(xvasp,FileMESSAGE,aflags,kflags,vflags);  // MODIFY ACCORDINGLY
                                    KBIN::XVASP_INCAR_Relax_Static_ON(xvasp,vflags);     // FIX
                                    // do the RWIGS ON
                                    if(vflags.KBIN_VASP_FORCE_OPTION_RWIGS_STATIC) KBIN::XVASP_INCAR_RWIGS_Static(xvasp,vflags,FileMESSAGE,ON);
                                    // done write INCAR
                                    aurostd::stringstream2file(xvasp.INCAR,string(xvasp.Directory+"/INCAR"));
                                    // NOW DO THE STATIC RUN
                                    if(vflags.KBIN_VASP_RUN.flag("STATIC_BANDS")) xvasp.NRELAXING=xvasp.NRELAX; //0;
                                    if(vflags.KBIN_VASP_RUN.flag("STATIC")) xvasp.NRELAXING=xvasp.NRELAX; // 0;
                                    xvasp.NRELAXING++;
                                    aus << aurostd::PaddedPRE(aurostd::utype2string(11111*xvasp.NRELAXING),5,"0") << "  STATIC (" << STRING_TO_SHOW << ") - " <<  xvasp.Directory 
                                        << " - K=[" << xvasp.str.kpoints_k1 << " " << xvasp.str.kpoints_k2 << " " << xvasp.str.kpoints_k3 << "]" << " - " << kflags.KBIN_BIN << " - " << Message("user,host,time") << endl;
                                    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                                    Krun=KBIN::VASP_Run(xvasp,aflags,kflags,vflags,FileMESSAGE);
                                    if(!Krun) {KBIN::VASP_Error(xvasp,FileMESSAGE,"EEEEE  runtime error [RELAX_STATIC_BANDS STATIC]");return Krun;}
                                    //	    if(_VASP_CONTCAR_SAVE_) KBIN::VASP_CONTCAR_Save(xvasp,string("static"));
                                    bool qmwrite=TRUE;
                                    KBIN::VASP_Backup(xvasp,qmwrite,string("static"));
                                    xvasp_spin_evolution.push_back(xvasp.str.qm_mag_atom); // keep track of spins
                                    aus << "00000  MESSAGE RESULT SPIN=" << xvasp_spin_evolution.at(xvasp_spin_evolution.size()-1) << " - " << Message(aflags,"user,host,time") << endl;
                                    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                                } // vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_BANDS") || vflags.KBIN_VASP_RUN.flag("STATIC_BANDS") || vflags.KBIN_VASP_RUN.flag("RELAX_STATIC") || vflags.KBIN_VASP_REPEAT.flag("REPEAT_STATIC_BANDS")
                                // BANDS PART ----------------------------------------------------------------------------
                                // BANDS PART ----------------------------------------------------------------------------
                                // BANDS PART ----------------------------------------------------------------------------
                                // BANDS PART ----------------------------------------------------------------------------
                                if(vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_BANDS") || vflags.KBIN_VASP_RUN.flag("STATIC_BANDS") || vflags.KBIN_VASP_REPEAT.flag("REPEAT_STATIC_BANDS") || vflags.KBIN_VASP_REPEAT.flag("REPEAT_BANDS")) {
                                    // NOW DO THE BANDS PATCHING KPOINTS (if necessary...)
                                    bool foundBZ;
                                    KBIN::VASP_Recycle(xvasp,"static");  // bring back the stuff
                                    KBIN::VASP_RecycleExtraFile(xvasp,"CHGCAR","static");  // bring back the stuff
                                    xvasp.aopts.flag("FLAG::CHGCAR_PRESERVED",TRUE); // in case of errors it is not lost bur recycled
                                    aus << "00000  MESSAGE Patching KPOINTS with BANDS LATTICE = \"" << vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE.content_string << "\" - " << Message(aflags,"user,host,time") << endl;
                                    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                                    // KBIN::VASP_Produce_KPOINTS(xvasp,AflowIn,FileAFLOWIN,FileMESSAGE,aflags,kflags,vflags);
                                    // KBIN::VASP_Modify_KPOINTS(xvasp,FileMESSAGE,aflags,vflags);
                                    // poscar was already conventionalized in the static part
                                    xvasp.KPOINTS.clear();xvasp.KPOINTS.str(std::string());
                                    //	      xvasp.KPOINTS <<
                                    string stringBZ;
                                    stringBZ=LATTICE::KPOINTS_Directions(vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE.content_string,rlattice,vflags.KBIN_VASP_KPOINTS_BANDS_GRID_VALUE,xvasp.str.iomode,foundBZ); // rlattice = updated structure
                                    // removed stuff BELOW
                                    xvasp.KPOINTS << stringBZ;
                                    aurostd::stringstream2file(xvasp.KPOINTS,string(xvasp.Directory+"/KPOINTS"));
                                    xvasp.aopts.flag("FLAG::KPOINTS_PRESERVED",TRUE); // don`t touch kpoints if there are flaws
                                    // NOW DO THE BANDS PATCHING INCAR
                                    aus << "00000  MESSAGE [" << STRING_TO_SHOW << "] Patching INCAR (bands_patching) " << Message(aflags,"user,host,time") << endl;
                                    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
                                    //  KBIN::VASP_Produce_INCAR(xvasp,AflowIn,FileAFLOWIN,FileMESSAGE,aflags,kflags,vflags); // BETTER than produce, SHOULD reread it
                                    KBIN::VASP_Reread_INCAR(xvasp,FileMESSAGE,aflags); // REREAD IT
                                    // KBIN::VASP_Modify_INCAR(xvasp,FileMESSAGE,aflags,kflags,vflags); // MODIFY ACCORDINGLY
                                    KBIN::XVASP_INCAR_Relax_Static_Bands_ON(xvasp,vflags);     // FIX
                                    // do the RWIGS OFF
                                    if(vflags.KBIN_VASP_FORCE_OPTION_RWIGS_STATIC)
                                        KBIN::XVASP_INCAR_RWIGS_Static(xvasp,vflags,FileMESSAGE,OFF);
                                    // done write INCAR
                                    aurostd::stringstream2file(xvasp.INCAR,string(xvasp.Directory+"/INCAR"));
                                    if(0)  {
                                        stringstream command;
                                        command << "cd " <<  xvasp.Directory << endl;
                                        command << "cat INCAR | grep -v NGXF | grep -v NGYF | grep -v NGZF > INCAR.new" << endl;
                                        command << "cat OUTCAR.static | grep NGXF | grep dimension | sed \"s/NG/\nNG/g\" | grep -v dimension | sed \"s/ //g\" >> INCAR.new" << endl;
                                        command << "mv INCAR.new INCAR " << endl;
                                        aurostd::execute(command);
                                    }
                                    // NOW DO THE BANDS RUN
                                    xvasp.NRELAXING++;
                                    aus << 11111*xvasp.NRELAXING << "  BANDS (" << STRING_TO_SHOW << ") - " <<  xvasp.Directory << " - K=[" << vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE.content_string << "," << vflags.KBIN_VASP_KPOINTS_BANDS_GRID_VALUE << "] - " << Message("user,host,time") << endl;
                                    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                                    vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.push("ROTMAT");	// dont mess up KPOINTS in bands
                                    vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.push("IBZKPT");    // dont mess up KPOINTS in bands
                                    vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.push("EDDRMM");	// dont mess up KPOINTS in bands
                                    Krun=KBIN::VASP_Run(xvasp,aflags,kflags,vflags,FileMESSAGE);
                                    if(!Krun) {KBIN::VASP_Error(xvasp,FileMESSAGE,"EEEEE  runtime error [RELAX_STATIC_BANDS BANDS]");return Krun;}
                                    //  if(_VASP_CONTCAR_SAVE_) KBIN::VASP_CONTCAR_Save(xvasp,string("bands"));
                                    bool qmwrite=FALSE;
                                    KBIN::VASP_Backup(xvasp,qmwrite,string("bands"));
                                    xvasp_spin_evolution.push_back(xvasp.str.qm_mag_atom); // keep track of spins
                                    aus << "00000  MESSAGE RESULT SPIN=" << xvasp_spin_evolution.at(xvasp_spin_evolution.size()-1) << " - " << Message(aflags,"user,host,time") << endl;
                                    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                                } // vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_BANDS") || vflags.KBIN_VASP_RUN.flag("STATIC_BANDS") || vflags.KBIN_VASP_REPEAT.flag("REPEAT_STATIC_BANDS") || vflags.KBIN_VASP_REPEAT.flag("REPEAT_BANDS")
                                // DIELECTRIC PART - DSCF ----------------------------------------------------------------
                                // DIELECTRIC PART -----------------------------------------------------------------------
                                // DIELECTRIC PART -----------------------------------------------------------------------
                                // DIELECTRIC PART -----------------------------------------------------------------------
                                xvasp.aopts.flag("FLAG::KPOINTS_PRESERVED",FALSE); // bands are done... I can refix the KPOINTS
                                if(vflags.KBIN_VASP_RUN.flag("RELAX_STATIC") || vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_BANDS") || vflags.KBIN_VASP_RUN.flag("STATIC_BANDS") || vflags.KBIN_VASP_RUN.flag("STATIC")) {
                                    // have static
                                    if(vflags.KBIN_VASP_RUN.flag("DIELECTRIC_STATIC")) {  // check for DIELECTRIC STATIC
                                        // check VASP version
                                        double dversion=0.0;
                                        string sversion=aurostd::execute2string("cat "+xvasp.Directory+"/OUTCAR.static | grep vasp | head -1 | sed \"s/ /\\n/g\" | grep vasp | sed \"s/vasp\\.//g\"");
                                        vector<string> tokens; aurostd::string2tokensAdd(sversion,tokens,".");
                                        if(tokens.size()>0) dversion+=aurostd::string2utype<double>(tokens.at(0));
                                        if(tokens.size()>1) dversion+=aurostd::string2utype<double>(tokens.at(1))/10.0;
                                        aus << "00000  MESSAGE Found VASP version=" << sversion << "  " << Message(aflags,"user,host,time") << endl;
                                        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                                        if(dversion<5.2) { // cant do it
                                            aus << "EEEEE  ERROR: Dielectric calculations need VASP >=5.2 " << Message(aflags,"user,host,time") << endl;
                                            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                                            Krun=FALSE;return Krun;}
                                        // PROCEED
                                        xvasp.NRELAXING++;
                                        aus << 11111*xvasp.NRELAXING << "  RUN_DIELECTRIC_STATIC (" << STRING_TO_SHOW << ") - " <<  xvasp.Directory << " - " << Message("user,host,time") << endl;
                                        kflags.KBIN_MPI_NCPUS_BUFFER=kflags.KBIN_MPI_NCPUS;
                                        if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_MATERIALS") && kflags.KBIN_MPI_NCPUS==24) {
                                            uint ncpus_before=kflags.KBIN_MPI_NCPUS;
                                            if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_MATERIALS")) kflags.KBIN_MPI_NCPUS=DUKE_MATERIALS_VASP5_CORES_DIELECTRIC;  // bug in mpivasp5
                                            aflags.AFLOW_GLOBAL_NCPUS=-kflags.KBIN_MPI_NCPUS;
                                            aus << "00000  MESSAGE Running RUN_DIELECTRIC_STATIC fixing mpivasp5 with " << ncpus_before << "-AMD cores to " << kflags.KBIN_MPI_NCPUS << "-AMD cores " << Message(aflags,"user,host,time") << endl;
                                        }	
                                        if((aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_BETA_MPICH") || aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_BETA_OPENMPI") || aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_AFLOWLIB"))
                                                && (kflags.KBIN_MPI_NCPUS==64 || kflags.KBIN_MPI_NCPUS==48 || kflags.KBIN_MPI_NCPUS==32)) {
                                            uint ncpus_before=kflags.KBIN_MPI_NCPUS;
                                            kflags.KBIN_MPI_NCPUS=DUKE_BETA_VASP5_CORES_DIELECTRIC; // something
                                            if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_BETA_MPICH")) kflags.KBIN_MPI_NCPUS=DUKE_BETA_VASP5_CORES_DIELECTRIC;  // bug in mpivasp5
                                            if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_BETA_OPENMPI")) kflags.KBIN_MPI_NCPUS=DUKE_BETA_VASP5_CORES_DIELECTRIC;  // bug in mpivasp5
                                            if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_AFLOWLIB")) kflags.KBIN_MPI_NCPUS=AFLOWLIB_VASP5_CORES_DIELECTRIC;  // bug in mpivasp5
                                            if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_QRATS_MPICH")) kflags.KBIN_MPI_NCPUS=AFLOWLIB_VASP5_CORES_DIELECTRIC;  // bug in mpivasp5
                                            if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_QFLOW_OPENMPI")) kflags.KBIN_MPI_NCPUS=AFLOWLIB_VASP5_CORES_DIELECTRIC;  // bug in mpivasp5
                                            if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MPCDF_EOS")) kflags.KBIN_MPI_NCPUS=AFLOWLIB_VASP5_CORES_DIELECTRIC;  // bug in mpivasp5
                                            if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MPCDF_DRACO")) kflags.KBIN_MPI_NCPUS=AFLOWLIB_VASP5_CORES_DIELECTRIC;  // bug in mpivasp5
                                            if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MPCDF_COBRA")) kflags.KBIN_MPI_NCPUS=AFLOWLIB_VASP5_CORES_DIELECTRIC;  // bug in mpivasp5
                                            if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MPCDF_HYDRA")) kflags.KBIN_MPI_NCPUS=AFLOWLIB_VASP5_CORES_DIELECTRIC;  // bug in mpivasp5
                                            if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DOD_CONRAD")) kflags.KBIN_MPI_NCPUS=AFLOWLIB_VASP5_CORES_DIELECTRIC;  // bug in mpivasp5 //DX 5/2/18 - DoD CONRAD
                                            aflags.AFLOW_GLOBAL_NCPUS=-kflags.KBIN_MPI_NCPUS;
                                            aus << "00000  MESSAGE Running RUN_DIELECTRIC_STATIC fixing mpivasp5 with " << ncpus_before << "-AMD cores to " << kflags.KBIN_MPI_NCPUS << "-AMD cores " << Message(aflags,"user,host,time") << endl;
                                        }	
                                        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                                        // a. Reuse the INCAR.static  NELM = 0, and remove the NBANDS parameter.
                                        // a. Set the k-point DKGRID < 0.10 (I've been using "aflow -k" for this).
                                        // b. Retain the following static run entries and their values: ALGO, LREAL, NSIM, ISYM, IBRION, NSW, NELM, NELMIN, ENMAX, ISPIN, ISMEAR, SIGMA, and everything LDA+U related.
                                        // c. Set NBANDS to a value that is around 10x the VASP default obtained in STEP 00.
                                        // d. Eliminate PSTRESS, EMIN, EMAX, LORBIT, ISIF, NEDOS.
                                        aus << "00000  MESSAGE Running RUN_DIELECTRIC_STATIC recycling static " << Message(aflags,"user,host,time") << endl;
                                        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                                        KBIN::VASP_Recycle(xvasp,"static");  // bring back the stuff
                                        KBIN::VASP_Reread_INCAR(xvasp,FileMESSAGE,aflags); // REREAD IT
                                        KBIN::VASP_Reread_KPOINTS(xvasp,FileMESSAGE,aflags); // REREAD IT
                                        KBIN::XVASP_INCAR_KPOINTS_Dielectric_SET(xvasp,kflags,vflags,"STATIC");     // FIX
                                        xvasp.aopts.flag("FLAG::XVASP_INCAR_generated",TRUE);xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
                                        aurostd::stringstream2file(xvasp.INCAR,string(xvasp.Directory+"/INCAR"));
                                        xvasp.aopts.flag("FLAG::XVASP_KPOINTS_generated",TRUE);xvasp.aopts.flag("FLAG::XVASP_KPOINTS_changed",TRUE);
                                        aurostd::stringstream2file(xvasp.KPOINTS,string(xvasp.Directory+"/KPOINTS"));
                                        Krun=KBIN::VASP_Run(xvasp,aflags,kflags,vflags,FileMESSAGE);
                                        if(!Krun) {KBIN::VASP_Error(xvasp,FileMESSAGE,"EEEEE  runtime error [RUN_DIELECTRIC_STATIC]");return Krun;}
                                        xvasp.aopts.flag("FLAG::WAVECAR_PRESERVED",TRUE); // WAVECAR.dielectric_static
                                        bool qmwrite=TRUE;
                                        KBIN::VASP_Backup(xvasp,qmwrite,string("dielectric_static"));
                                        //		kflags.KBIN_MPI_NCPUS=kflags.KBIN_MPI_NCPUS_BUFFER;
                                    }
                                    if(vflags.KBIN_VASP_RUN.flag("DIELECTRIC_DYNAMIC") && vflags.KBIN_VASP_RUN.flag("DIELECTRIC_STATIC")) {  // check for DIELECTRIC DYNAMIC
                                        xvasp.NRELAXING++;
                                        aus << 11111*xvasp.NRELAXING << "  RUN_DIELECTRIC_DYNAMIC (" << STRING_TO_SHOW << ") - " <<  xvasp.Directory << " - " << Message("user,host,time") << endl;
                                        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                                        // a. Reuse the STEP 01 WAVECAR + ALGO=EXACT  NELM=1 LOPTICS=.TRUE. CSHIFT=0.15 OMEGAMAX=25 NEDOS=12500  Remove LEPSILON and LRPA
                                        aus << "00000  MESSAGE Running RUN_DIELECTRIC_DYNAMIC recycling dielectric_static " << Message(aflags,"user,host,time") << endl;
                                        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                                        KBIN::VASP_Recycle(xvasp,"dielectric_static"); // bring back the stuff
                                        KBIN::VASP_RecycleExtraFile(xvasp,"WAVECAR","dielectric_static");  // bring back the stuff
                                        KBIN::VASP_Reread_INCAR(xvasp,FileMESSAGE,aflags); // REREAD IT
                                        KBIN::XVASP_INCAR_KPOINTS_Dielectric_SET(xvasp,kflags,vflags,"DYNAMIC");   // FIX
                                        xvasp.aopts.flag("FLAG::XVASP_INCAR_generated",TRUE);xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
                                        aurostd::stringstream2file(xvasp.INCAR,string(xvasp.Directory+"/INCAR"));
                                        Krun=KBIN::VASP_Run(xvasp,aflags,kflags,vflags,FileMESSAGE);
                                        if(!Krun) {KBIN::VASP_Error(xvasp,FileMESSAGE,"EEEEE  runtime error [RUN_DIELECTRIC_DYNAMIC]");return Krun;}
                                        aurostd::execute("rm -f "+xvasp.Directory+"/WAVECAR.dielectric_static");
                                        xvasp.aopts.flag("FLAG::WAVECAR_PRESERVED",FALSE); // all gone
                                        bool qmwrite=TRUE;
                                        KBIN::VASP_Backup(xvasp,qmwrite,string("dielectric_dynamic"));
                                    }
                                    if(vflags.KBIN_VASP_RUN.flag("DSCF") && vflags.KBIN_VASP_RUN.flag("DIELECTRIC_DYNAMIC") && vflags.KBIN_VASP_RUN.flag("DIELECTRIC_STATIC")) {  // check for DIELECTRIC DYNAMIC
                                        xvasp.NRELAXING++;
                                        aus << 11111*xvasp.NRELAXING << "  RUN_DSCF (" << STRING_TO_SHOW << ") - " <<  xvasp.Directory << " - " << Message("user,host,time") << endl;
                                        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                                        exit(0);
                                    }
                                }
                                // FINISHED
                                xvasp.NRELAXING++;
                                aus << 11111*xvasp.NRELAXING << "  END (" << STRING_TO_SHOW << ")        - " <<  xvasp.Directory << " - " << Message("user,host,time") << endl;
                                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                                vflags.KBIN_VASP_RUN.flag("STATIC",FALSE);  // put back the options

                                // CLEAN-UP BY WAHYU
                                ostringstream xaus;
                                xaus << "cd " << xvasp.Directory << endl;
                                // if(vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_BANDS")) xaus << "rm -f CHG.relax* CHGCAR.relax* POTCAR.* CHGCAR.bands CHG.bands" << endl;
                                // if(vflags.KBIN_VASP_RUN.flag("STATIC_BANDS"))  xaus << "rm -f POTCAR.* CHGCAR.bands CHG.bands" << endl;
                                if(vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_BANDS")) xaus << "rm -f CHG.relax* CHGCAR.relax* CHGCAR.bands CHG.bands" << endl;
                                if(vflags.KBIN_VASP_RUN.flag("STATIC_BANDS"))  xaus << "rm -f CHGCAR.bands CHG.bands" << endl;
                                if(vflags.KBIN_VASP_REPEAT.flag("REPEAT_STATIC_BANDS"))  xaus << "rm -f CHGCAR.bands CHG.bands" << endl;
                                if(vflags.KBIN_VASP_REPEAT.flag("REPEAT_STATIC_BANDS")) {;}
                                aurostd::execute(xaus);
                                // done ....
                                // wahyu, you might ask something more here
                            } // vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_BANDS") || vflags.KBIN_VASP_RUN.flag("STATIC_BANDS") || vflags.KBIN_VASP_REPEAT.flag("REPEAT_STATIC_BANDS") || vflags.KBIN_VASP_REPEAT.flag("REPEAT_BANDS")
                            // --------------------------------------------------------------------------------------------------------------------
                            // --------------------------------------------------------------------------------------------------------------------
                            // --------------------------------------------------------------------------------------------------------------------
                            // REPEAT_DELSOL REPEAT_DELSOL REPEAT_DELSOL REPEAT_DELSOL REPEAT_DELSOL REPEAT_DELSOL 	
                            // REPEAT_DELSOL REPEAT_DELSOL REPEAT_DELSOL REPEAT_DELSOL REPEAT_DELSOL REPEAT_DELSOL 	
                            // REPEAT_DELSOL REPEAT_DELSOL REPEAT_DELSOL REPEAT_DELSOL REPEAT_DELSOL REPEAT_DELSOL
                            // PRL 105, 196403 (2010)
                            if(vflags.KBIN_VASP_REPEAT.flag("REPEAT_DELSOL")) {
                                bool delsol_d=FALSE;
                                float NELECT=0,Nr=0;
                                xmatrix<double> rlattice(3,3);
                                string STRING_TO_SHOW="",stmp="";
                                string fnamedelsol=xvasp.Directory+string("/delsol.tmp");
                                stringstream command,strdelsol;
                                ifstream fdelsol;
                                if(vflags.KBIN_VASP_REPEAT.flag("REPEAT_DELSOL")) STRING_TO_SHOW="REPEAT_DELSOL";
                                aus << "00000  MESSAGE MODE= (" << STRING_TO_SHOW << ") - " << xvasp.Directory << " - " << Message("user,host,time") << endl;
                                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                                // LOAD FORMER LOCK
                                if(aurostd::FileExist(xvasp.Directory+string("/REPEAT_DELSOL"))) {
                                    stringstream lock_recycled;
                                    aurostd::file2stringstream(xvasp.Directory+"/REPEAT_DELSOL",lock_recycled);
                                    aus << "XXXXX ---------------------------------------------------------------------------------------------- " << endl;
                                    aus << "XXXXX FORMER LOCK BEGIN, recycled (" << STRING_TO_SHOW << ") - " << xvasp.Directory << " - " << Message("user,host,time") << endl;
                                    aus << lock_recycled.str();
                                    aus << "XXXXX FORMER LOCK END, recycled (" << STRING_TO_SHOW << ") - " << xvasp.Directory << " - " << Message("user,host,time") << endl;
                                    aus << "XXXXX ---------------------------------------------------------------------------------------------- " << endl;
                                    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                                }
                                // UNZIP EVERYTHING	      
                                for(uint iext=0;iext<vext.size();iext++) { 
                                    // command.clear();command.str(std::string());  // ORIGINAL
                                    // command << "cd " <<  xvasp.Directory << endl; command << "bzip2" << " -dfq *bz2 " << endl;  // ORIGINAL
                                    // aurostd::execute(command);  // ORIGINAL
                                    aus << "DEBUG - KBIN::VASP_Directory: EXT POINT [3] " << endl; aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);		  
                                    aurostd::execute(vzip.at(iext)+" -dfq `find "+xvasp.Directory+" -name \"*"+vext.at(iext)+"\"`");
                                }		

                                // cp INCAR, POSCAR, KPOINTS, POTCAR from *.static
                                KBIN::VASP_Recycle(xvasp,"static");  // bring back the stuff
                                //Scanning whether it is sp or spd from the POTCAR
                                command.clear();command.str(std::string());
                                command << "cd " <<  xvasp.Directory << endl;
                                command << "grep VRHFIN POTCAR.static | sed \'s/:/\\n/g\' | grep -v VRHFIN > delsol.tmp" << endl;
                                aurostd::execute(command);	
                                strdelsol.clear();strdelsol.str(std::string());
                                aurostd::file2stringstream(xvasp.Directory+"/delsol.tmp",strdelsol);
                                command.clear();command.str(std::string());
                                command << "rm -f " << xvasp.Directory << "/delsol.tmp" << endl;
                                aurostd::execute(command);
                                delsol_d=FALSE;
                                if((aurostd::substring2bool(strdelsol.str(),"d"))) delsol_d=TRUE;
                                //Scanning NELECT from OUTCAR.static
                                command.clear();command.str(std::string());
                                command << "cd " << xvasp.Directory << endl;
                                command << "grep NELECT OUTCAR.static | sed \'s/=/\\n/g\' | grep -v NELECT > delsol.tmp" << endl;
                                aurostd::execute(command);
                                strdelsol.clear();strdelsol.str(std::string());
                                aurostd::file2stringstream(xvasp.Directory+"/delsol.tmp",strdelsol);
                                command.clear();command.str(std::string());
                                command << "rm -f " << xvasp.Directory << "/delsol.tmp" << endl;
                                aurostd::execute(command);
                                strdelsol>>NELECT;
                                //if(NELECT<1.0) ;//need to add error handling here
                                Nr=NELECT/68.0;
                                if(delsol_d) Nr=NELECT/72.0;
                                aus << "DELSOL: NELECT N0=" << NELECT << endl << "DELSOL: NELECT n=" << Nr << endl;
                                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    

                                // NOW MODIFY THE INCAR
                                aus << "00000  MESSAGE [" << STRING_TO_SHOW << "] modifying INCAR (delsol_patching) " << Message(aflags,"user,host,time") << endl;
                                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
                                KBIN::VASP_Reread_INCAR(xvasp,FileMESSAGE,aflags); // REREAD IT
                                stmp="NELECT="+aurostd::utype2string(NELECT+Nr);
                                xvasp.INCAR << aurostd::PaddedPOST(stmp,_incarpad_) << " # NELECT = N0 + n for DELSOL plus" << endl;
                                // do the RWIGS OFF
                                if(vflags.KBIN_VASP_FORCE_OPTION_RWIGS_STATIC)
                                    KBIN::XVASP_INCAR_RWIGS_Static(xvasp,vflags,FileMESSAGE,OFF);
                                aurostd::stringstream2file(xvasp.INCAR,string(xvasp.Directory+"/INCAR"));

                                /* 		//Reread POSCAR
                                        aus << "00000  MESSAGE [" << STRING_TO_SHOW << "] rereading POSCAR (delsol) " << Message(aflags,"user,host,time") << endl;
                                        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
                                        KBIN::VASP_Reread_POSCAR(xvasp,FileMESSAGE,aflags); // REREAD IT
                                //Reread KPOINTS
                                aus << "00000  MESSAGE [" << STRING_TO_SHOW << "] rereading KPOINTS (delsol) " << Message(aflags,"user,host,time") << endl;
                                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
                                KBIN::VASP_Reread_KPOINTS(xvasp,FileMESSAGE,aflags); // REREAD IT
                                */
                                KBIN::VASP_RecycleExtraFile(xvasp,"POSCAR","static");  // bring back the stuff
                                KBIN::VASP_RecycleExtraFile(xvasp,"KPOINTS","static");  // bring back the stuff

                                // NOW RUN DELSOL plus
                                uint vrelax=7;
                                aus << 11111*vrelax << "  DELSOL plus (" << STRING_TO_SHOW << ") - " <<  xvasp.Directory << " - " << stmp << " "<< Message("user,host,time") << endl;
                                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                                vrelax++;
                                Krun=KBIN::VASP_Run(xvasp,aflags,kflags,vflags,FileMESSAGE);
                                if(!Krun) {KBIN::VASP_Error(xvasp,FileMESSAGE,"EEEEE  runtime error [DELSOL]");return Krun;}
                                //	    if(_VASP_CONTCAR_SAVE_) KBIN::VASP_CONTCAR_Save(xvasp,string("dsolp"));
                                bool qmwrite=FALSE;
                                KBIN::VASP_Backup(xvasp,qmwrite,string("dsolp"));

                                //NOW DO DELSOL minus
                                KBIN::VASP_Recycle(xvasp,"static");
                                KBIN::VASP_Reread_INCAR(xvasp,FileMESSAGE,aflags);
                                stmp="NELECT="+aurostd::utype2string(NELECT-Nr);
                                xvasp.INCAR << aurostd::PaddedPOST(stmp,_incarpad_) << " # NELECT = N0 - n for DELSOL minus" << endl;
                                // do the RWIGS OFF
                                if(vflags.KBIN_VASP_FORCE_OPTION_RWIGS_STATIC)
                                    KBIN::XVASP_INCAR_RWIGS_Static(xvasp,vflags,FileMESSAGE,OFF);
                                aurostd::stringstream2file(xvasp.INCAR,string(xvasp.Directory+"/INCAR"));

                                KBIN::VASP_RecycleExtraFile(xvasp,"POSCAR","static");  // bring back the stuff
                                KBIN::VASP_RecycleExtraFile(xvasp,"KPOINTS","static");  // bring back the stuff
                                // NOW RUN DELSOL minus
                                aus << 11111*vrelax << "  DELSOL minus (" << STRING_TO_SHOW << ") - " <<  xvasp.Directory << " - " << stmp << " "<< Message("user,host,time") << endl;
                                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                                vrelax++;
                                Krun=KBIN::VASP_Run(xvasp,aflags,kflags,vflags,FileMESSAGE);
                                if(!Krun) {KBIN::VASP_Error(xvasp,FileMESSAGE,"EEEEE  runtime error [DELSOL minus]");return Krun;}
                                //	    if(_VASP_CONTCAR_SAVE_) KBIN::VASP_CONTCAR_Save(xvasp,string("dsolm"));
                                qmwrite=FALSE;
                                KBIN::VASP_Backup(xvasp,qmwrite,string("dsolm"));		
                                // FINISHED
                                aus << 11111*vrelax << "  END (" << STRING_TO_SHOW << ")        - " <<  xvasp.Directory << " - " << Message("user,host,time") << endl;
                                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                                vflags.KBIN_VASP_RUN.flag("STATIC",FALSE);  // put back the options
                                // CLEAN-UP BY WAHYU
                                ostringstream xaus;
                                xaus << "cd " << xvasp.Directory << endl;
                                if(vflags.KBIN_VASP_REPEAT.flag("REPEAT_BANDS"))  xaus << "rm -f CHGCAR.dsol* CHG.dsol*" << endl;
                                aurostd::execute(xaus);
                                // done ....
                                // wahyu, you might ask something more here
                            } // (vflags.KBIN_VASP_REPEAT.flag("REPEAT_DELSOL")
                            //************* END OF REPEAT_DELSOL ************************
                            //************* END OF REPEAT_DELSOL ************************
                            //************* END OF REPEAT_DELSOL ************************
                            //************* END OF REPEAT_DELSOL ************************
                            //************* END OF REPEAT_DELSOL ************************
                            //************* END OF REPEAT_DELSOL ************************
                            // --------------------------------------------------------------------------------------------------------------------
                            // KPOINTS KPOINTS KPOINTS
                            if(vflags.KBIN_VASP_RUN.flag("KPOINTS")) {            // NON THREADED
                                KBIN::VASP_Write_INPUT(xvasp,vflags); // VASP VASP WRITE
                                xvasp.aopts.flag("FLAG::XVASP_KPOINTS_changed",TRUE);                                                                // BACKUP KPOINTS
                                aurostd::stringstream2file(xvasp.KPOINTS_orig,string(xvasp.Directory+"/KPOINTS.orig"));   // BACKUP KPOINTS
                                int kbak_k1=xvasp.str.kpoints_k1;
                                int kbak_k2=xvasp.str.kpoints_k2;
                                int kbak_k3=xvasp.str.kpoints_k3;
                                //	    int kbak_kmax=xvasp.str.kpoints_kmax;   kbak_kmax=max(kbak_k1,kbak_k2,kbak_k3);
                                int kk1,kk2,kk3;
                                string relax,relaxfile;
                                // 1st step: 1/2
                                kk1=(kbak_k1+1)/2;kk2=(kbak_k2+1)/2;kk3=(kbak_k3+1)/2;
                                relax="11111a ";relaxfile="relax1";
                                aus << relax << "RELAXATION - " << xvasp.Directory << " - K=[" << kk1 << " " << kk2 << " " << kk3 << "]" << " - " << Message("user,host,time") << endl;
                                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
                                xvasp.str.kpoints_k1=min(kk1,kbak_k1);xvasp.str.kpoints_k2=min(kk2,kbak_k2);xvasp.str.kpoints_k3=min(kk3,kbak_k3);
                                xvasp.str.kpoints_kmax=max(xvasp.str.kpoints_k1,xvasp.str.kpoints_k2,xvasp.str.kpoints_k3);
                                KBIN::XVASP_numbers2string(xvasp);
                                aurostd::stringstream2file(xvasp.KPOINTS,string(xvasp.Directory+"/KPOINTS"));   // BACKUP KPOINTS
                                Krun=KBIN::VASP_Run(xvasp,aflags,kflags,vflags,relaxfile,relaxfile,FALSE,FileMESSAGE);
                                if(!Krun) {KBIN::VASP_Error(xvasp,FileMESSAGE,"EEEEE  runtime error [KPOINTS 1]");return Krun;}
                                // 2nd step: 1/2
                                kk1=3*(kbak_k1+1)/4;kk2=3*(kbak_k2+1)/4;kk3=3*(kbak_k3+1)/4;
                                relax="11111b ";relaxfile="relax1";
                                aus << relax << "RELAXATION - " << xvasp.Directory << " - K=[" << kk1 << " " << kk2 << " " << kk3 << "]" << " - " << Message("user,host,time") << endl;
                                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
                                xvasp.str.kpoints_k1=min(kk1,kbak_k1);xvasp.str.kpoints_k2=min(kk2,kbak_k2);xvasp.str.kpoints_k3=min(kk3,kbak_k3);
                                xvasp.str.kpoints_kmax=max(xvasp.str.kpoints_k1,xvasp.str.kpoints_k2,xvasp.str.kpoints_k3);
                                KBIN::XVASP_numbers2string(xvasp);
                                aurostd::stringstream2file(xvasp.KPOINTS,string(xvasp.Directory+"/KPOINTS"));   // BACKUP KPOINTS
                                Krun=KBIN::VASP_Run(xvasp,aflags,kflags,vflags,relaxfile,relaxfile,FALSE,FileMESSAGE);
                                if(!Krun) {KBIN::VASP_Error(xvasp,FileMESSAGE,"EEEEE  runtime error [KPOINTS 2]");return Krun;}
                                // 3rd step: 1/2
                                kk1=kbak_k1;kk2=kbak_k2;kk3=kbak_k3;
                                relax="11111c ";relaxfile="relax1";
                                aus << relax << "RELAXATION - " << xvasp.Directory << " - K=[" << kk1 << " " << kk2 << " " << kk3 << "]" << " - " << Message("user,host,time") << endl;
                                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
                                xvasp.str.kpoints_k1=min(kk1,kbak_k1);xvasp.str.kpoints_k2=min(kk2,kbak_k2);xvasp.str.kpoints_k3=min(kk3,kbak_k3);
                                xvasp.str.kpoints_kmax=max(xvasp.str.kpoints_k1,xvasp.str.kpoints_k2,xvasp.str.kpoints_k3);
                                KBIN::XVASP_numbers2string(xvasp);
                                aurostd::stringstream2file(xvasp.KPOINTS,string(xvasp.Directory+"/KPOINTS"));   // BACKUP KPOINTS
                                // ----------------------------------------------------------
                                // with PAWGGA2
                                if(PAWGGA2) {
                                    // STEP 1
                                    aus << "11111  RELAXATION - " <<  xvasp.Directory << " - K=[" << xvasp.str.kpoints_k1 << " " << xvasp.str.kpoints_k2 << " " << xvasp.str.kpoints_k3 << "]" << " - " << kflags.KBIN_BIN << " - " << Message("user,host,time") << endl;
                                    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
                                    Krun=KBIN::VASP_Run(xvasp,aflags,kflags,vflags,"relax2paw_gga",FALSE,FileMESSAGE);
                                    if(!Krun) {KBIN::VASP_Error(xvasp,FileMESSAGE,"EEEEE  runtime error [KPOINTS PAWGGA2]");return Krun;}
                                    aus << "22222  END        - " <<  xvasp.Directory << " - K=[" << xvasp.str.kpoints_k1 << " " << xvasp.str.kpoints_k2 << " " << xvasp.str.kpoints_k3 << "]" << " - " << kflags.KBIN_BIN << " - " << Message("user,host,time") << endl;
                                    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
                                }
                                // ----------------------------------------------------------
                                // norma, without PAWGGA2
                                if(!PAWGGA2) {
                                    int vrelax=1;
                                    for(int i=1;i<=xvasp.NRELAX;i++) {
                                        aus << 11111*vrelax << "   RELAXATION - " << xvasp.Directory << " - K=[" << kk1 << " " << kk2 << " " << kk3
                                            << "]" << " - " << kflags.KBIN_BIN << " - " << Message("user,host,time") << endl;
                                        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
                                        if(vflags.KBIN_VASP_FORCE_OPTION_LDAU_ADIABATIC.content_int>0) KBIN::XVASP_INCAR_LDAU_ADIABATIC(xvasp,i); // ADIABATIC
                                        if(i<xvasp.NRELAX)  {
                                            Krun=KBIN::VASP_Run(xvasp,aflags,kflags,vflags,"relax"+aurostd::utype2string(vrelax),"relax"+aurostd::utype2string(vrelax),TRUE,FileMESSAGE);
                                            if(!Krun) {KBIN::VASP_Error(xvasp,FileMESSAGE,"EEEEE  runtime error [KPOINT 4]");return Krun;}
                                        }
                                        if(i==xvasp.NRELAX) {
                                            Krun=KBIN::VASP_Run(xvasp,aflags,kflags,vflags,"relax"+aurostd::utype2string(vrelax),TRUE,FileMESSAGE);
                                            if(!Krun) {KBIN::VASP_Error(xvasp,FileMESSAGE,"EEEEE  runtime error [KPOINTS 5]");return Krun;}
                                        }
                                        vrelax++;
                                    }
                                    aus << 11111*vrelax << "   END        - " << xvasp.Directory << " - " << kflags.KBIN_BIN << " - " << Message("user,host,time") << endl;
                                    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
                                }
                            }
                            // ***************************************************************************
                            // POSTSCRIPT
                            if(!vflags.KBIN_VASP_POSCAR_MODE.flag("EXPLICIT_START_STOP_POINT"))
                                if(kflags.AFLOW_MODE_POSTSCRIPT_EXPLICIT || kflags.AFLOW_MODE_POSTSCRIPT_EXPLICIT_START_STOP)
                                    KBIN::RUN_DirectoryScript(aflags,DEFAULT_AFLOW_POSTSCRIPT_COMMAND,DEFAULT_AFLOW_POSTSCRIPT_OUT);
                            // ***************************************************************************
                        }
                    }
                }
            }
            // **********
            // some verbose
            if(vflags.KBIN_VASP_POSCAR_MODE.flag("EXPLICIT_START_STOP_POINT")) {
                aus << "00000  MESSAGE END loop " << xvasp.POSCAR_index << "/" << vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING.size()
                    << " - " << Message(aflags,"user,host,time") << endl;
                aus << "00000  MESSAGE END loop in directory =" << xvasp.Directory << " - " << Message(aflags,"user,host,time") << endl;
                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                // compress the subdirectories
                if(Krun && kflags.KZIP_COMPRESS) KBIN::CompressDirectory(aflags,kflags);
            }
            aflags=aflags_backup;kflags=kflags_backup; // RESTORE
            } // LOOP ixvasp
            // ***************************************************************************
            aflags=aflags_backup;kflags=kflags_backup; // RESTORE
            // POSTSCRIPT
            if(vflags.KBIN_VASP_POSCAR_MODE.flag("EXPLICIT_START_STOP_POINT"))
                if(kflags.AFLOW_MODE_POSTSCRIPT_EXPLICIT || kflags.AFLOW_MODE_POSTSCRIPT_EXPLICIT_START_STOP)
                    KBIN::RUN_DirectoryScript(aflags,DEFAULT_AFLOW_POSTSCRIPT_COMMAND,DEFAULT_AFLOW_POSTSCRIPT_OUT);
            // ***************************************************************************
            FileAFLOWIN.clear();FileAFLOWIN.close();
            return Krun;
        }
    } // namespace

    // ******************************************************************************************************************************************************
    // ******************************************************************************************************************************************************

    int CheckStringInFile(string FileIn,string str,int PID) {
        int out;
        ostringstream aus_exec;
        aus_exec << "cat " << FileIn << " | grep -c \"" << str << "\" > aflow.out." << PID << endl;
        aurostd::execute(aus_exec);
        ifstream FileCHECK;
        stringstream FineNameAflowTmpPID;
        FineNameAflowTmpPID << "aflow.out." << PID;
        FileCHECK.open(FineNameAflowTmpPID.str().c_str(),std::ios::in);
        FileCHECK >> out;
        FileCHECK.clear();FileCHECK.close();
        aus_exec << "rm -f aflow.out." << PID << endl;
        aurostd::execute(aus_exec);
        return out;
    }

    // ******************************************************************************************************************************************************
    // ******************************************************************************************************************************************************
    // FLAG Class for KBIN_VASP_Run

    namespace KBIN {
        bool VASP_Run(_xvasp &xvasp,_aflags &aflags,_kflags &kflags,_vflags &vflags,ofstream &FileMESSAGE) {        // AFLOW_FUNCTION_IMPLEMENTATION
            bool LDEBUG=(FALSE || XHOST.DEBUG);
            if(LDEBUG) cerr << "KBIN::VASP_Run: " << Message("time") << "  BEGIN" << endl;

            ostringstream aus_exec,aus;
            xoption xwarning,xfixed,xmessage;
            int nbands,counter_ZPOTRF=0;
            double enmax;
            bool vasp_start=TRUE;
            aurostd::StringstreamClean(aus_exec);
            aurostd::StringstreamClean(aus);
            int nrun=0,maxrun=15;
            int fix_NIRMAT=0;
            int kpoints_k1=xvasp.str.kpoints_k1; double kpoints_s1=xvasp.str.kpoints_s1;
            int kpoints_k2=xvasp.str.kpoints_k2; double kpoints_s2=xvasp.str.kpoints_s2;
            int kpoints_k3=xvasp.str.kpoints_k3; double kpoints_s3=xvasp.str.kpoints_s3;

            // get CPUS from PBS/SLURM
            // string ausenv;
            aus << "DDDDD  PBS_NUM_PPN=" << XHOST.PBS_NUM_PPN << "  PBS_NNODES=" << XHOST.PBS_NNODES << " - " << Message(aflags,"user,host,time") << endl;
            aus << "DDDDD  SLURM_CPUS_ON_NODE=" << XHOST.SLURM_CPUS_ON_NODE << "  SLURM_NNODES=" << XHOST.SLURM_NNODES << " - " << Message(aflags,"user,host,time") << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_BETA_MPICH")) { 	//corey
                if(kflags.KBIN_MPI_NCPUS==0) kflags.KBIN_MPI_NCPUS=XHOST.PBS_NUM_PPN;}
            if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_BETA_OPENMPI")) {kflags.KBIN_MPI_NCPUS=XHOST.PBS_NUM_PPN; if(kflags.KBIN_MPI_NCPUS==0) kflags.KBIN_MPI_NCPUS=XHOST.CPU_Cores;}
            if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_QRATS_MPICH")) {	//corey
                if(kflags.KBIN_MPI_NCPUS==0) kflags.KBIN_MPI_NCPUS=XHOST.PBS_NUM_PPN;}
            if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_QFLOW_OPENMPI")) {	//corey
                if(kflags.KBIN_MPI_NCPUS==0) kflags.KBIN_MPI_NCPUS=XHOST.PBS_NUM_PPN;}
            if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MPCDF_EOS")) {	//corey
                if(kflags.KBIN_MPI_NCPUS==0) kflags.KBIN_MPI_NCPUS=XHOST.SLURM_CPUS_ON_NODE;}
            if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MPCDF_DRACO")) {	//corey
                if(kflags.KBIN_MPI_NCPUS==0) kflags.KBIN_MPI_NCPUS=XHOST.SLURM_CPUS_ON_NODE;}
            if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MPCDF_COBRA")) {	//corey
                if(kflags.KBIN_MPI_NCPUS==0) kflags.KBIN_MPI_NCPUS=XHOST.SLURM_CPUS_ON_NODE;}
            if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MPCDF_HYDRA")) {	//corey
                if(kflags.KBIN_MPI_NCPUS==0) kflags.KBIN_MPI_NCPUS=XHOST.SLURM_CPUS_ON_NODE;}
            if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DOD_CONRAD")) {
                if(kflags.KBIN_MPI_NCPUS==0) kflags.KBIN_MPI_NCPUS=XHOST.PBS_NUM_PPN;}
            if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::TERAGRID_RANGER")) {kflags.KBIN_MPI_NCPUS=XHOST.PBS_NNODES; if(kflags.KBIN_MPI_NCPUS==0) kflags.KBIN_MPI_NCPUS=XHOST.CPU_Cores;}
            if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::TERAGRID_KRAKEN")) {kflags.KBIN_MPI_NCPUS=XHOST.PBS_NNODES; if(kflags.KBIN_MPI_NCPUS==0) kflags.KBIN_MPI_NCPUS=XHOST.CPU_Cores;}

            // for reducint CPUs on the fly
            if(aflags.AFLOW_GLOBAL_NCPUS<0) kflags.KBIN_MPI_NCPUS=-aflags.AFLOW_GLOBAL_NCPUS; // this to force things on reducing CPUS

            // for for LS coupling
            if(vflags.KBIN_VASP_FORCE_OPTION_LSCOUPLING.option) {
                if(!aurostd::substring2bool(kflags.KBIN_BIN,VASPLS_BIN_POSTFIX_DEFAULT)) kflags.KBIN_BIN=kflags.KBIN_BIN+VASPLS_BIN_POSTFIX_DEFAULT; // standard LS
                if(!aurostd::substring2bool(kflags.KBIN_MPI_BIN,VASPLS_BIN_POSTFIX_DEFAULT)) kflags.KBIN_MPI_BIN=kflags.KBIN_MPI_BIN+VASPLS_BIN_POSTFIX_DEFAULT; // standard LS
                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);	
                aus << "00000  MESSAGE SPIN-ORBIT TYPE CALCULATIONS , adding " << VASPLS_BIN_POSTFIX_DEFAULT << " to BIN " << " - " << Message(aflags,"user,host,time") << endl;
                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            }
            if(kflags.KBIN_MPI) kflags.KBIN_BIN=kflags.KBIN_MPI_BIN; // forcing, no matter what

            while(vasp_start) {
                // ********* RUN VASP                
                { // ERRORS
                    bool error=FALSE;
                    if(aurostd::FileEmpty(xvasp.Directory+"/INCAR"))   {KBIN::VASP_Error(xvasp,FileMESSAGE,"EEEEE  ERROR KBIN::VASP_Run: "+Message("time")+"  Empty INCAR ");error=TRUE;return FALSE;}
                    if(aurostd::FileEmpty(xvasp.Directory+"/POSCAR"))  {KBIN::VASP_Error(xvasp,FileMESSAGE,"EEEEE  ERROR KBIN::VASP_Run: "+Message("time")+"  Empty POSCAR ");error=TRUE;return FALSE;}
                    if(aurostd::FileEmpty(xvasp.Directory+"/KPOINTS")) {KBIN::VASP_Error(xvasp,FileMESSAGE,"EEEEE  ERROR KBIN::VASP_Run: "+Message("time")+"  Empty KPOINTS ");error=TRUE;return FALSE;}
                    if(aurostd::FileEmpty(xvasp.Directory+"/POTCAR"))  {KBIN::VASP_Error(xvasp,FileMESSAGE,"EEEEE  ERROR KBIN::VASP_Run: "+Message("time")+"  Empty POTCAR ");error=TRUE;return FALSE;}
                    if(error) return FALSE;

                    if(LDEBUG) cerr << "KBIN::VASP_Run: " << Message("time") << "  [1]" << endl;

                    // FIX INCAR if alternating
                    if(vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("IONS_CELL_VOLUME")) {
                        if(aurostd::_isodd(xvasp.NRELAXING))  aus << "00000  MESSAGE Alternating: RELAX_CELL_VOLUME - " << Message(aflags,"user,host,time") << endl;
                        if(aurostd::_iseven(xvasp.NRELAXING)) aus << "00000  MESSAGE Alternating: RELAX_IONS - " << Message(aflags,"user,host,time") << endl;
                        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                        vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("STATIC",FALSE);
                        vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("ALL",FALSE);
                        vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("IONS",FALSE);
                        vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("CELL_SHAPE",FALSE);
                        vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("CELL_VOLUME",FALSE);
                        vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.push("IONS_CELL_VOLUME");
                        KBIN::XVASP_INCAR_Relax_ON(xvasp,vflags,xvasp.NRELAXING);
                    }
                    // RUN VASP NON QUEUE ------------------------------------------------------------------------
                    if(kflags.KBIN_QSUB==FALSE) {
                        nrun++;
                        aus_exec << "cd " << xvasp.Directory << endl;
                        aus_exec << "rm -f vasp.out" << endl;
                        if(kflags.KBIN_MPI==FALSE) {
                            aus << "00000  MESSAGE SERIAL job - [" << xvasp.str.atoms.size() << "atoms] - " << Message(aflags,"user,host,time") << endl;
                            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                            aus_exec << kflags.KBIN_BIN << " > vasp.out " << endl;
                            aus << "00000  MESSAGE Executing: " << kflags.KBIN_BIN << " > vasp.out " << Message(aflags,"user,host,time,memory") << endl;  // CO 170628 - SLOW WITH MEMORY
                            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                            aurostd::execute(aus_exec);
                            aurostd::Sleep(_KVASP_VASP_SLEEP_);
                        } else {
                            aus << "00000  MESSAGE MPI PARALLEL job - [" << xvasp.str.atoms.size() << "atoms] - " << " MPI=" << kflags.KBIN_MPI_NCPUS << "CPUs  " << Message(aflags,"user,host,time") << endl;
                            if(kflags.KBIN_MPI_OPTIONS!="") aus << "00000  MESSAGE MPI OPTIONS=[" << kflags.KBIN_MPI_OPTIONS << "]" << Message(aflags,"user,host,time") << endl;	      
                            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);	
                            if(LDEBUG) cerr << "KBIN::VASP_Run: " << Message("time") << "  aflags.AFLOW_MACHINE_GLOBAL=" << aflags.AFLOW_MACHINE_GLOBAL << endl;
                            if(LDEBUG) cerr << "KBIN::VASP_Run: " << Message("time") << "  aflags.AFLOW_MACHINE_LOCAL=" << aflags.AFLOW_MACHINE_LOCAL << endl;
                            // NO HOST ------------------------------------------------------------------------
                            if(!aflags.AFLOW_MACHINE_LOCAL.flag()) {
                                aus << "00000  MESSAGE Executing: ";
                                if(!kflags.KBIN_MPI_OPTIONS.empty()){
                                    aus_exec << kflags.KBIN_MPI_OPTIONS << endl;
                                    aus << kflags.KBIN_MPI_OPTIONS << "; ";
                                }
                                if(!kflags.KBIN_MPI_START.empty()){
                                    aus_exec << kflags.KBIN_MPI_START << " > vasp.out " << endl;
                                    aus << kflags.KBIN_MPI_START << " > vasp.out; ";
                                }
                                aus_exec << kflags.KBIN_MPI_COMMAND << " " << kflags.KBIN_MPI_NCPUS << " " << kflags.KBIN_MPI_BIN << " >> vasp.out " << endl;
                                aus << kflags.KBIN_MPI_COMMAND << " " << kflags.KBIN_MPI_NCPUS << " " << kflags.KBIN_MPI_BIN << " >> vasp.out; ";
                                if(!kflags.KBIN_MPI_STOP.empty()){
                                    aus_exec << kflags.KBIN_MPI_STOP << " >> vasp.out " << endl;
                                    aus << kflags.KBIN_MPI_STOP << " >> vasp.out ";
                                }
                                aus << Message(aflags,"user,host,time,memory") << endl; // CO 170628 - SLOW WITH MEMORY
                                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);	
                                aurostd::execute(aus_exec);
                            }
                            // HOST DUKE_BETA_MPICH ------------------------------------------------------------------------
                            if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_BETA_MPICH")) {
                                // verbosization
                                aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL << "  MPI PARALLEL job - [" << xvasp.str.atoms.size() << "atoms] - " << " MPI=" << kflags.KBIN_MPI_NCPUS << "CPUs  " << Message(aflags,"user,host,time","aflow_kvasp.cpp") << endl;
                                aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL << "  Executing: " << MPI_COMMAND_DUKE_BETA_MPICH << " " << kflags.KBIN_MPI_NCPUS << " " << MPI_BINARY_DIR_DUKE_BETA_MPICH << kflags.KBIN_MPI_BIN << " >> vasp.out " << Message(aflags,"user,host,time,memory","aflow_kvasp.cpp") << endl;  // CO 170628 - SLOW WITH MEMORY
                                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                                // run
                                aus_exec << kflags.KBIN_MPI_OPTIONS << endl;
                                aus_exec << MPI_OPTIONS_DUKE_BETA_MPICH << endl;
                                aus_exec << MPI_COMMAND_DUKE_BETA_MPICH << " " << kflags.KBIN_MPI_NCPUS << " " << MPI_BINARY_DIR_DUKE_BETA_MPICH << kflags.KBIN_MPI_BIN << " >> vasp.out " << endl;
                                //	    aurostd::PrintMessageStream(FileMESSAGE,aus_exec,XHOST.QUIET);
                                aurostd::execute(aus_exec);
                            }
                            // HOST DUKE_BETA_OPENMPI ------------------------------------------------------------------------
                            if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_BETA_OPENMPI")) {
                                if(!aurostd::substring2bool(kflags.KBIN_MPI_BIN,"_openmpi")) kflags.KBIN_MPI_BIN=kflags.KBIN_MPI_BIN+"_openmpi"; // fix the OPENMPI
                                // verbosization
                                aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL << "  MPI PARALLEL job - [" << xvasp.str.atoms.size() << "atoms] - " << " MPI=" << kflags.KBIN_MPI_NCPUS << "CPUs  " << Message(aflags,"user,host,time","aflow_kvasp.cpp") << endl;
                                aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL << "  Executing: " << MPI_COMMAND_DUKE_BETA_OPENMPI << " " << kflags.KBIN_MPI_NCPUS << " " << MPI_BINARY_DIR_DUKE_BETA_OPENMPI << kflags.KBIN_MPI_BIN << " >> vasp.out " << Message(aflags,"user,host,time,memory","aflow_kvasp.cpp") << endl;  // CO 170628 - SLOW WITH MEMORY
                                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                                // run
                                aus_exec << kflags.KBIN_MPI_OPTIONS << endl;
                                aus_exec << MPI_OPTIONS_DUKE_BETA_OPENMPI << endl;
                                aus_exec << MPI_COMMAND_DUKE_BETA_OPENMPI << " " << kflags.KBIN_MPI_NCPUS << " " << MPI_BINARY_DIR_DUKE_BETA_OPENMPI << kflags.KBIN_MPI_BIN << " >> vasp.out " << endl;
                                //	    aurostd::PrintMessageStream(FileMESSAGE,aus_exec,XHOST.QUIET);
                                aurostd::execute(aus_exec);
                            }
                            // HOST DUKE_QRATS_MPICH ------------------------------------------------------------------------
                            if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_QRATS_MPICH")) {
                                // verbosization
                                aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL << "  MPI PARALLEL job - [" << xvasp.str.atoms.size() << "atoms] - " << " MPI=" << kflags.KBIN_MPI_NCPUS << "CPUs  " << Message(aflags,"user,host,time","aflow_kvasp.cpp") << endl;
                                aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL << "  Executing: " << MPI_COMMAND_DUKE_QRATS_MPICH << " " << kflags.KBIN_MPI_NCPUS << " " << MPI_BINARY_DIR_DUKE_QRATS_MPICH << kflags.KBIN_MPI_BIN << " >> vasp.out " << Message(aflags,"user,host,time,memory","aflow_kvasp.cpp") << endl;  // CO 170628 - SLOW WITH MEMORY
                                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                                // run
                                aus_exec << kflags.KBIN_MPI_OPTIONS << endl;
                                aus_exec << MPI_OPTIONS_DUKE_QRATS_MPICH << endl;
                                aus_exec << MPI_COMMAND_DUKE_QRATS_MPICH << " " << kflags.KBIN_MPI_NCPUS << " " << MPI_BINARY_DIR_DUKE_QRATS_MPICH << kflags.KBIN_MPI_BIN << " >> vasp.out " << endl;
                                //	    aurostd::PrintMessageStream(FileMESSAGE,aus_exec,XHOST.QUIET);
                                aurostd::execute(aus_exec);
                            }
                            // HOST DUKE_QFLOW_OPENMPI ------------------------------------------------------------------------
                            if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_QFLOW_OPENMPI")) {
                                // verbosization
                                aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL << "  MPI PARALLEL job - [" << xvasp.str.atoms.size() << "atoms] - " << " MPI=" << kflags.KBIN_MPI_NCPUS << "CPUs  " << Message(aflags,"user,host,time","aflow_kvasp.cpp") << endl;
                                aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL << "  Executing: " << MPI_COMMAND_DUKE_QFLOW_OPENMPI << " " << kflags.KBIN_MPI_NCPUS << " " << MPI_BINARY_DIR_DUKE_QFLOW_OPENMPI << kflags.KBIN_MPI_BIN << " >> vasp.out " << Message(aflags,"user,host,time,memory","aflow_kvasp.cpp") << endl;  // CO 170628 - SLOW WITH MEMORY
                                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                                // run
                                aus_exec << kflags.KBIN_MPI_OPTIONS << endl;
                                aus_exec << MPI_OPTIONS_DUKE_QFLOW_OPENMPI << endl;
                                aus_exec << MPI_COMMAND_DUKE_QFLOW_OPENMPI << " " << kflags.KBIN_MPI_NCPUS << " " << MPI_BINARY_DIR_DUKE_QFLOW_OPENMPI << kflags.KBIN_MPI_BIN << " >> vasp.out " << endl;
                                //	    aurostd::PrintMessageStream(FileMESSAGE,aus_exec,XHOST.QUIET);
                                aurostd::execute(aus_exec);
                            }
                            // HOST MPCDF_EOS_MPI ------------------------------------------------------------------------
                            if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MPCDF_EOS")) {
                                // verbosization
                                int local_NCPUS=kflags.KBIN_MPI_NCPUS;
                                if(MPI_NCPUS_MPCDF_EOS>0) {
                                    local_NCPUS=MPI_NCPUS_MPCDF_EOS;
                                    aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL << "  Forcing: kflags.KBIN_MPI_NCPUS=MPI_NCPUS_MPCDF_EOS" << Message(aflags,"user,host,time","aflow_kvasp.cpp") << endl;
                                    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                                }
                                if(MPI_HYPERTHREADING_MPCDF_EOS=="FALSE" || MPI_HYPERTHREADING_MPCDF_EOS=="OFF") {
                                    local_NCPUS=local_NCPUS/2;
                                    aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL << "  Forcing: HYPERTHREADING = OFF" << Message(aflags,"user,host,time","aflow_kvasp.cpp") << endl;
                                    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                                }
                                if(MPI_HYPERTHREADING_MPCDF_EOS=="TRUE" || MPI_HYPERTHREADING_MPCDF_EOS=="ON") {
                                    local_NCPUS=local_NCPUS*2;
                                    aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL << "  Forcing: HYPERTHREADING = ON" << Message(aflags,"user,host,time","aflow_kvasp.cpp") << endl;
                                    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                                }
                                aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL << "  MPI PARALLEL job - [" << xvasp.str.atoms.size() << "atoms] - " << " MPI=" << local_NCPUS << "CPUs  " << Message(aflags,"user,host,time","aflow_kvasp.cpp") << endl;
                                aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL << "  Executing: " << MPI_COMMAND_MPCDF_EOS << " " << local_NCPUS << " " << MPI_BINARY_DIR_MPCDF_EOS << kflags.KBIN_MPI_BIN << " >> vasp.out " << Message(aflags,"user,host,time,memory","aflow_kvasp.cpp") << endl;  // CO 170628 - SLOW WITH MEMORY
                                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                                // run
                                aus_exec << kflags.KBIN_MPI_OPTIONS << endl;
                                aus_exec << MPI_OPTIONS_MPCDF_EOS << endl;
                                aus_exec << MPI_COMMAND_MPCDF_EOS << " " << local_NCPUS << " " << MPI_BINARY_DIR_MPCDF_EOS << kflags.KBIN_MPI_BIN << " >> vasp.out " << endl;
                                //	    aurostd::PrintMessageStream(FileMESSAGE,aus_exec,XHOST.QUIET);
                                aurostd::execute(aus_exec);
                            }
                            // HOST MPCDF_DRACO_MPI ------------------------------------------------------------------------
                            if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MPCDF_DRACO")) {
                                // verbosization
                                int local_NCPUS=kflags.KBIN_MPI_NCPUS;
                                if(MPI_NCPUS_MPCDF_DRACO>0) {
                                    local_NCPUS=MPI_NCPUS_MPCDF_DRACO;
                                    aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL << "  Forcing: kflags.KBIN_MPI_NCPUS=MPI_NCPUS_MPCDF_DRACO" << Message(aflags,"user,host,time","aflow_kvasp.cpp") << endl;
                                    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                                }
                                if(MPI_HYPERTHREADING_MPCDF_DRACO=="FALSE" || MPI_HYPERTHREADING_MPCDF_DRACO=="OFF") {
                                    local_NCPUS=local_NCPUS/2;
                                    aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL << "  Forcing: HYPERTHREADING = OFF" << Message(aflags,"user,host,time","aflow_kvasp.cpp") << endl;
                                    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                                }
                                if(MPI_HYPERTHREADING_MPCDF_DRACO=="TRUE" || MPI_HYPERTHREADING_MPCDF_DRACO=="ON") {
                                    local_NCPUS=local_NCPUS*2;
                                    aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL << "  Forcing: HYPERTHREADING = ON" << Message(aflags,"user,host,time","aflow_kvasp.cpp") << endl;
                                    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                                }
                                aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL << "  MPI PARALLEL job - [" << xvasp.str.atoms.size() << "atoms] - " << " MPI=" << local_NCPUS << "CPUs  " << Message(aflags,"user,host,time","aflow_kvasp.cpp") << endl;
                                aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL << "  Executing: " << MPI_COMMAND_MPCDF_DRACO << " " << local_NCPUS << " " << MPI_BINARY_DIR_MPCDF_DRACO << kflags.KBIN_MPI_BIN << " >> vasp.out " << Message(aflags,"user,host,time,memory","aflow_kvasp.cpp") << endl;  // CO 170628 - SLOW WITH MEMORY
                                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                                // run
                                aus_exec << kflags.KBIN_MPI_OPTIONS << endl;
                                aus_exec << MPI_OPTIONS_MPCDF_DRACO << endl;
                                aus_exec << MPI_COMMAND_MPCDF_DRACO << " " << local_NCPUS << " " << MPI_BINARY_DIR_MPCDF_DRACO << kflags.KBIN_MPI_BIN << " >> vasp.out " << endl;
                                //	    aurostd::PrintMessageStream(FileMESSAGE,aus_exec,XHOST.QUIET);
                                aurostd::execute(aus_exec);
                            }
                            // HOST MPCDF_COBRA_MPI ------------------------------------------------------------------------
                            if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MPCDF_COBRA")) {
                                // verbosization 
                                int local_NCPUS=kflags.KBIN_MPI_NCPUS;
                                if(MPI_NCPUS_MPCDF_COBRA>0) {
                                    local_NCPUS=MPI_NCPUS_MPCDF_COBRA;
                                    aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL << "  Forcing: kflags.KBIN_MPI_NCPUS=MPI_NCPUS_MPCDF_COBRA" << Message(aflags,"user,host,time","aflow_kvasp.cpp") << endl;
                                    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                                }
                                if(MPI_HYPERTHREADING_MPCDF_COBRA=="FALSE" || MPI_HYPERTHREADING_MPCDF_COBRA=="OFF") {
                                    local_NCPUS=local_NCPUS/2;
                                    aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL << "  Forcing: HYPERTHREADING = OFF" << Message(aflags,"user,host,time","aflow_kvasp.cpp") << endl;
                                    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                                }
                                if(MPI_HYPERTHREADING_MPCDF_COBRA=="TRUE" || MPI_HYPERTHREADING_MPCDF_COBRA=="ON") {
                                    local_NCPUS=local_NCPUS*2;
                                    aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL << "  Forcing: HYPERTHREADING = ON" << Message(aflags,"user,host,time","aflow_kvasp.cpp") << endl;
                                    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                                }
                                aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL << "  MPI PARALLEL job - [" << xvasp.str.atoms.size() << "atoms] - " << " MPI=" << local_NCPUS << "CPUs  " << Message(aflags,"user,host,time","aflow_kvasp.cpp") << endl;
                                aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL << "  Executing: " << MPI_COMMAND_MPCDF_COBRA << " " << local_NCPUS << " " << MPI_BINARY_DIR_MPCDF_COBRA << kflags.KBIN_MPI_BIN << " >> vasp.out " << Message(aflags,"user,host,time,memory","aflow_kvasp.cpp") << endl;  // CO 170628 - SLOW WITH MEMORY
                                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                                // run
                                aus_exec << kflags.KBIN_MPI_OPTIONS << endl;
                                aus_exec << MPI_OPTIONS_MPCDF_COBRA << endl;
                                aus_exec << MPI_COMMAND_MPCDF_COBRA << " " << local_NCPUS << " " << MPI_BINARY_DIR_MPCDF_COBRA << kflags.KBIN_MPI_BIN << " >> vasp.out " << endl;
                                //	    aurostd::PrintMessageStream(FileMESSAGE,aus_exec,XHOST.QUIET);
                                aurostd::execute(aus_exec);
                            }
                            // HOST MPCDF_HYDRA_MPI ------------------------------------------------------------------------
                            if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MPCDF_HYDRA")) {
                                // verbosization 
                                int local_NCPUS=kflags.KBIN_MPI_NCPUS;
                                if(MPI_NCPUS_MPCDF_HYDRA>0) {
                                    local_NCPUS=MPI_NCPUS_MPCDF_HYDRA;
                                    aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL << "  Forcing: kflags.KBIN_MPI_NCPUS=MPI_NCPUS_MPCDF_HYDRA" << Message(aflags,"user,host,time","aflow_kvasp.cpp") << endl;
                                    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                                }
                                if(MPI_HYPERTHREADING_MPCDF_HYDRA=="FALSE" || MPI_HYPERTHREADING_MPCDF_HYDRA=="OFF") {
                                    local_NCPUS=local_NCPUS/2;
                                    aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL << "  Forcing: HYPERTHREADING = OFF" << Message(aflags,"user,host,time","aflow_kvasp.cpp") << endl;
                                    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                                }
                                if(MPI_HYPERTHREADING_MPCDF_HYDRA=="TRUE" || MPI_HYPERTHREADING_MPCDF_HYDRA=="ON") {
                                    local_NCPUS=local_NCPUS*2;
                                    aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL << "  Forcing: HYPERTHREADING = ON" << Message(aflags,"user,host,time","aflow_kvasp.cpp") << endl;
                                    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                                }
                                aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL << "  MPI PARALLEL job - [" << xvasp.str.atoms.size() << "atoms] - " << " MPI=" << local_NCPUS << "CPUs  " << Message(aflags,"user,host,time","aflow_kvasp.cpp") << endl;
                                //	      aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL << "  Executing: " << MPI_COMMAND_MPCDF_HYDRA << " " << local_NCPUS << " " << MPI_BINARY_DIR_MPCDF_HYDRA << kflags.KBIN_MPI_BIN << " >> vasp.out " << Message(aflags,"user,host,time,memory","aflow_kvasp.cpp") << endl;  
                                aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL << "  Executing: " << MPI_COMMAND_MPCDF_HYDRA << " " << MPI_BINARY_DIR_MPCDF_HYDRA << kflags.KBIN_MPI_BIN << " >> vasp.out " << Message(aflags,"user,host,time,memory","aflow_kvasp.cpp") << endl;   // poe not MPI run
                                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                                // run
                                aus_exec << kflags.KBIN_MPI_OPTIONS << endl;
                                aus_exec << MPI_OPTIONS_MPCDF_HYDRA << endl;
                                //	      aus_exec << MPI_COMMAND_MPCDF_HYDRA << " " << local_NCPUS << " " << MPI_BINARY_DIR_MPCDF_HYDRA << kflags.KBIN_MPI_BIN << " >> vasp.out " << endl;
                                aus_exec << MPI_COMMAND_MPCDF_HYDRA << " " << MPI_BINARY_DIR_MPCDF_HYDRA << kflags.KBIN_MPI_BIN << " >> vasp.out " << endl;  // poe not MPI run
                                //	    aurostd::PrintMessageStream(FileMESSAGE,aus_exec,XHOST.QUIET);
                                aurostd::execute(aus_exec);
                            }
                            // HOST CONRAD_MPICH ------------------------------------------------------------------------
                            if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DOD_CONRAD")) {
                                // verbosization
                                aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL << "  MPI PARALLEL job - [" << xvasp.str.atoms.size() << "atoms] - " << " MPI=" << kflags.KBIN_MPI_NCPUS << "CPUs  " << Message(aflags,"user,host,time","aflow_kvasp.cpp") << endl;
                                aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL << "  Executing: " << MPI_COMMAND_DOD_CONRAD << " " << kflags.KBIN_MPI_NCPUS << " " << MPI_BINARY_DIR_DOD_CONRAD << kflags.KBIN_MPI_BIN << " >> vasp.out " << Message(aflags,"user,host,time,memory","aflow_kvasp.cpp") << endl;  //CO 170628 - SLOW WITH MEMORY
                                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                                // run
                                aus_exec << kflags.KBIN_MPI_OPTIONS << endl;
                                aus_exec << MPI_OPTIONS_DOD_CONRAD << endl;
                                aus_exec << MPI_COMMAND_DOD_CONRAD << " " << kflags.KBIN_MPI_NCPUS << " " << MPI_BINARY_DIR_DOD_CONRAD << kflags.KBIN_MPI_BIN << " >> vasp.out " << endl;
                                aurostd::execute(aus_exec);
                            }
                            // HOST DUKE_MATERIALS ------------------------------------------------------------------------
                            if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_MATERIALS")) {
                                // verbosization
                                aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL << "  MPI PARALLEL job - [" << xvasp.str.atoms.size() << "atoms] - " << " MPI=" << kflags.KBIN_MPI_NCPUS << "CPUs  " << Message(aflags,"user,host,time","aflow_kvasp.cpp") << endl;
                                aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL << "  Executing: " << MPI_COMMAND_DUKE_MATERIALS << " " << kflags.KBIN_MPI_NCPUS << " " << MPI_BINARY_DIR_DUKE_MATERIALS << kflags.KBIN_MPI_BIN << " >> vasp.out " << Message(aflags,"user,host,time,memory","aflow_kvasp.cpp") << endl;  // CO 170628 - SLOW WITH MEMORY
                                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                                // run
                                aus_exec << kflags.KBIN_MPI_OPTIONS << endl;
                                aus_exec << MPI_OPTIONS_DUKE_MATERIALS << endl;
                                aus_exec << MPI_COMMAND_DUKE_MATERIALS << " " << kflags.KBIN_MPI_NCPUS << " " << MPI_BINARY_DIR_DUKE_MATERIALS << kflags.KBIN_MPI_BIN << " >> vasp.out " << endl;
                                aurostd::execute(aus_exec);
                            }
                            // HOST DUKE_AFLOWLIB ------------------------------------------------------------------------
                            if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_AFLOWLIB")) {
                                // verbosization
                                aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL << "  MPI PARALLEL job - [" << xvasp.str.atoms.size() << "atoms] - " << " MPI=" << kflags.KBIN_MPI_NCPUS << "CPUs  " << Message(aflags,"user,host,time","aflow_kvasp.cpp") << endl;
                                aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL << "  Executing: " << MPI_COMMAND_DUKE_AFLOWLIB << " " << kflags.KBIN_MPI_NCPUS << " " << MPI_BINARY_DIR_DUKE_AFLOWLIB << kflags.KBIN_MPI_BIN << " >> vasp.out " << Message(aflags,"user,host,time,memory","aflow_kvasp.cpp") << endl;  // CO 170628 - SLOW WITH MEMORY
                                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                                // run
                                aus_exec << kflags.KBIN_MPI_OPTIONS << endl;
                                aus_exec << MPI_OPTIONS_DUKE_AFLOWLIB << endl;
                                aus_exec << MPI_COMMAND_DUKE_AFLOWLIB << " " << kflags.KBIN_MPI_NCPUS << " " << MPI_BINARY_DIR_DUKE_AFLOWLIB << kflags.KBIN_MPI_BIN << " >> vasp.out " << endl;
                                aurostd::execute(aus_exec);
                            }
                            // HOST DUKE_HABANA ------------------------------------------------------------------------
                            if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_HABANA")) {
                                // verbosization
                                aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL << "  MPI PARALLEL job - [" << xvasp.str.atoms.size() << "atoms] - " << " MPI=" << kflags.KBIN_MPI_NCPUS << "CPUs  " << Message(aflags,"user,host,time","aflow_kvasp.cpp") << endl;
                                aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL << "  Executing: " << MPI_COMMAND_DUKE_HABANA << " " << kflags.KBIN_MPI_NCPUS << " " << MPI_BINARY_DIR_DUKE_HABANA << kflags.KBIN_MPI_BIN << " >> vasp.out " << Message(aflags,"user,host,time,memory","aflow_kvasp.cpp") << endl;  // CO 170628 - SLOW WITH MEMORY
                                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                                // run
                                aus_exec << kflags.KBIN_MPI_OPTIONS << endl;
                                aus_exec << MPI_OPTIONS_DUKE_HABANA << endl;
                                aus_exec << MPI_COMMAND_DUKE_HABANA << " " << kflags.KBIN_MPI_NCPUS << " " << MPI_BINARY_DIR_DUKE_HABANA << kflags.KBIN_MPI_BIN << " >> vasp.out " << endl;
                                aurostd::execute(aus_exec);
                            }
                            // HOST TERAGRID_RANGER ------------------------------------------------------------------------
                            if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::TERAGRID_RANGER")) {
                                // verbosization
                                aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL << "  MPI PARALLEL job - [" << xvasp.str.atoms.size() << "atoms] - " << " MPI=" << kflags.KBIN_MPI_NCPUS << "CPUs  " << Message(aflags,"user,host,time","aflow_kvasp.cpp") << endl;
                                aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL << "  Executing: " << MPI_COMMAND_TERAGRID_RANGER << " " << MPI_BINARY_DIR_TERAGRID_RANGER << kflags.KBIN_MPI_BIN << " >> vasp.out " << Message(aflags,"user,host,time,memory","aflow_kvasp.cpp") << endl;  // CO 170628 - SLOW WITH MEMORY
                                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                                // run
                                aus_exec << kflags.KBIN_MPI_OPTIONS << endl;
                                aus_exec << MPI_COMMAND_TERAGRID_RANGER << " " << MPI_BINARY_DIR_TERAGRID_RANGER << kflags.KBIN_MPI_BIN << " >> vasp.out " << endl;
                                aurostd::execute(aus_exec);
                            }
                            // HOST TERAGRID_KRAKEN ------------------------------------------------------------------------
                            if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::TERAGRID_KRAKEN")) {
                                // verbosization
                                aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL << "  MPI PARALLEL job - [" << xvasp.str.atoms.size() << "atoms] - " << " MPI=" << kflags.KBIN_MPI_NCPUS << "CPUs  " << Message(aflags,"user,host,time","aflow_kvasp.cpp") << endl;
                                aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL << "  Executing: " << MPI_COMMAND_TERAGRID_KRAKEN << " " << kflags.KBIN_MPI_NCPUS << " " << MPI_BINARY_DIR_TERAGRID_KRAKEN << kflags.KBIN_MPI_BIN << " >> vasp.out " << Message(aflags,"user,host,time,memory","aflow_kvasp.cpp") << endl;  // CO 170628 - SLOW WITH MEMORY
                                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                                // run
                                aus_exec << kflags.KBIN_MPI_OPTIONS << endl;
                                aus_exec << MPI_COMMAND_TERAGRID_KRAKEN << " " << kflags.KBIN_MPI_NCPUS << " " << MPI_BINARY_DIR_TERAGRID_KRAKEN << kflags.KBIN_MPI_BIN << " >> vasp.out " << endl;
                                aurostd::execute(aus_exec);
                            }
                            // HOST TRINITY_PARSONS ------------------------------------------------------------------------
                            if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::TRINITY_PARSONS")) {
                                // verbosization
                                aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL << "  MPI PARALLEL job - [" << xvasp.str.atoms.size() << "atoms] - " << " MPI=" << kflags.KBIN_MPI_NCPUS << "CPUs  " << Message(aflags,"user,host,time","aflow_kvasp.cpp") << endl;
                                aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL << "  Executing: " << MPI_COMMAND_TRINITY_PARSONS << " " << kflags.KBIN_MPI_NCPUS << " " << MPI_BINARY_DIR_TRINITY_PARSONS << kflags.KBIN_MPI_BIN << " >> vasp.out " << Message(aflags,"user,host,time,memory","aflow_kvasp.cpp") << endl;  // CO 170628 - SLOW WITH MEMORY
                                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                                // run
                                aus_exec << kflags.KBIN_MPI_OPTIONS << endl;
                                aus_exec << MPI_COMMAND_TRINITY_PARSONS << " " << kflags.KBIN_MPI_NCPUS << " " << MPI_BINARY_DIR_TRINITY_PARSONS << kflags.KBIN_MPI_BIN << " >> vasp.out " << endl;
                                aurostd::execute(aus_exec);
                            }
                            // HOST FULTON_MARYLOU ------------------------------------------------------------------------
                            if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::FULTON_MARYLOU")) {
                                // verbosization
                                aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL << "  MPI PARALLEL job - [" << xvasp.str.atoms.size() << "atoms] - " << " MPI=" << kflags.KBIN_MPI_NCPUS << "CPUs  " << Message(aflags,"user,host,time","aflow_kvasp.cpp") << endl;
                                //	      aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL << "  Executing: " << MPI_COMMAND_FULTON_MARYLOU << " " << kflags.KBIN_MPI_NCPUS << " " << MPI_BINARY_DIR_FULTON_MARYLOU << kflags.KBIN_MPI_BIN << " >> vasp.out " << Message(aflags,"user,host,time,memory","aflow_kvasp.cpp") << endl;  
                                aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL << "  Executing: " << MPI_COMMAND_FULTON_MARYLOU << " " << MPI_BINARY_DIR_FULTON_MARYLOU << kflags.KBIN_MPI_BIN << " >> vasp.out " << Message(aflags,"user,host,time,memory","aflow_kvasp.cpp") << endl;  
                                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                                // run
                                aus_exec << kflags.KBIN_MPI_OPTIONS << endl;
                                aus_exec << MPI_OPTIONS_FULTON_MARYLOU << endl;
                                aus_exec << MPI_COMMAND_FULTON_MARYLOU << " " << MPI_BINARY_DIR_FULTON_MARYLOU << kflags.KBIN_MPI_BIN << " >> vasp.out " << endl;
                                // aus_exec << MPI_COMMAND_FULTON_MARYLOU << " " << kflags.KBIN_MPI_NCPUS << " " << MPI_BINARY_DIR_FULTON_MARYLOU << kflags.KBIN_MPI_BIN << " >> vasp.out " << endl;  // with --np
                                aurostd::execute(aus_exec);
                            }
                            // HOST OHAD ------------------------------------------------------------------------
                            if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::OHAD")) {
                                // verbosization
                                aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL << "  MPI PARALLEL job - [" << xvasp.str.atoms.size() << "atoms] - " << " MPI=" << kflags.KBIN_MPI_NCPUS << "CPUs  " << Message(aflags,"user,host,time","aflow_kvasp.cpp") << endl;
                                aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL << "  Executing: " << MPI_COMMAND_MACHINE2 << " " << MPI_BINARY_DIR_MACHINE2 << kflags.KBIN_MPI_BIN << " >> vasp.out " << Message(aflags,"user,host,time,memory","aflow_kvasp.cpp") << endl;  // CO 170628 - SLOW WITH MEMORY
                                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                                // run
                                aus_exec << kflags.KBIN_MPI_OPTIONS << endl;
                                aus_exec << MPI_COMMAND_MACHINE2 << " " << MPI_BINARY_DIR_MACHINE2 << kflags.KBIN_MPI_BIN << " >> vasp.out " << endl;
                                aurostd::execute(aus_exec);
                            }
                            // HOST HOST1 ------------------------------------------------------------------------
                            if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::HOST1")) {
                                // verbosization
                                aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL << "  MPI PARALLEL job - [" << xvasp.str.atoms.size() << "atoms] - " << " MPI=" << kflags.KBIN_MPI_NCPUS << "CPUs  " << Message(aflags,"user,host,time","aflow_kvasp.cpp") << endl;
                                aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL << "  Executing: " << MPI_COMMAND_MACHINE1 << " " << MPI_BINARY_DIR_MACHINE1 << kflags.KBIN_MPI_BIN << " >> vasp.out " << Message(aflags,"user,host,time,memory","aflow_kvasp.cpp") << endl;  // CO 170628 - SLOW WITH MEMORY
                                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                                // run
                                aus_exec << kflags.KBIN_MPI_OPTIONS << endl;
                                aus_exec << MPI_COMMAND_MACHINE1 << " " << MPI_BINARY_DIR_MACHINE1 << kflags.KBIN_MPI_BIN << " >> vasp.out " << endl;
                                aurostd::execute(aus_exec);
                            }
                            // DONE ------------------------------------------------------------------------
                        }
                        aurostd::Sleep(_KVASP_VASP_SLEEP_);
                        vasp_start=FALSE;
                    }
                    // RUN VASP QUEUED ------------------------------------------------------------------------
                    if(kflags.KBIN_QSUB) {
                        nrun++;
                        aus_exec << "cd " << xvasp.Directory << endl;
                        if(kflags.KBIN_MPI==FALSE) {
                            aus << "00000  MESSAGE QUEUED SERIAL job - [" << xvasp.str.atoms.size() << "atoms] - " << Message(aflags,"user,host,time") << endl;
                            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                            aurostd::RemoveFile(string(xvasp.Directory+"/aflow.qsub.done"));
                            aurostd::stringstream2file(xvasp.xqsub.QSUB,string(xvasp.Directory+"/aflow.qsub.run"));
                            aurostd::ChmodFile("755",string(xvasp.Directory+"/aflow.qsub.run"));
                            //	  exit(0);
                            aus_exec << kflags.KBIN_QSUB_COMMAND << " " << kflags.KBIN_QSUB_PARAMS << " " << "./aflow.qsub.run &" << endl;
                            aurostd::execute(aus_exec);
                            KBIN::QSUB_WaitFinished(aflags,FileMESSAGE,FALSE);
                            aurostd::RemoveFile(string(xvasp.Directory+"/aflow.qsub.done"));
                        } else {
                            aus << "00000  MESSAGE QUEUED MPI PARALLEL job - [" << xvasp.str.atoms.size() << "atoms] - " << " MPI=" << kflags.KBIN_MPI_NCPUS << "CPUs  " << Message(aflags,"user,host,time") << endl;
                            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                            aurostd::RemoveFile(string(xvasp.Directory+"/aflow.qsub.done"));
                            aurostd::stringstream2file(xvasp.xqsub.QSUB,string(xvasp.Directory+"/aflow.qsub.run"));
                            aurostd::ChmodFile("755",string(xvasp.Directory+"/aflow.qsub.run"));
                            aus_exec << kflags.KBIN_QSUB_COMMAND << " " << kflags.KBIN_QSUB_PARAMS << " " << "./aflow.qsub.run &" << endl;
                            aurostd::execute(aus_exec);
                            KBIN::QSUB_WaitFinished(aflags,FileMESSAGE,FALSE);
                            aurostd::RemoveFile(string(xvasp.Directory+"/aflow.qsub.done"));
                        }	
                        aurostd::Sleep(_KVASP_VASP_SLEEP_);
                        vasp_start=FALSE;
                    }
                }
                if(LDEBUG) cerr << "KBIN::VASP_Run: " << Message("time") << "  [2]" << endl;

                if(aurostd::FileEmpty(xvasp.Directory+"/vasp.out"))  {KBIN::VASP_Error(xvasp,FileMESSAGE,"EEEEE  ERROR KBIN::VASP_Run: "+Message("time")+"  Empty vasp.out ");return FALSE;}
                if(aurostd::FileEmpty(xvasp.Directory+"/OUTCAR"))  {KBIN::VASP_Error(xvasp,FileMESSAGE,"EEEEE  ERROR KBIN::VASP_Run: "+Message("time")+"  Empty OUTCAR ");return FALSE;}
                // DONT CHECK CONTCAR it can be empty
                // DONT CHECK OSZICAR it can be empty

                // update kpoints table

                kpoints_k1=xvasp.str.kpoints_k1; kpoints_s1=xvasp.str.kpoints_s1;
                kpoints_k2=xvasp.str.kpoints_k2; kpoints_s2=xvasp.str.kpoints_s2;
                kpoints_k3=xvasp.str.kpoints_k3; kpoints_s3=xvasp.str.kpoints_s3;

                // ***************** CHECK FOR ERRORS *********
                if(LDEBUG) cerr << "KBIN::VASP_Run: " << Message("time") << "  [3a]  nrun=" << nrun<< endl;
                if(LDEBUG) cerr << "KBIN::VASP_Run: " << Message("time") << "  [3b]  maxrun=" << maxrun<< endl;

                // check VASP version
                double DVERSION=0.0;
                xwarning.push_attached("SVERSION",aurostd::execute2string("cat "+xvasp.Directory+"/OUTCAR | grep vasp | head -1 | sed \"s/ /\\n/g\" | grep vasp | sed \"s/vasp\\.//g\""));
                vector<string> vtokens; aurostd::string2tokensAdd(xwarning.getattachedscheme("SVERSION"),vtokens,".");
                if(vtokens.size()>0) DVERSION+=aurostd::string2utype<double>(vtokens.at(0));
                if(vtokens.size()>1) DVERSION+=aurostd::string2utype<double>(vtokens.at(1))/10.0;
                xwarning.push_attached("DVERSION",aurostd::utype2string((double) DVERSION));

                //check OUTCAR converged?
                string stmp_toten = aurostd::execute2string("grep 'TOTEN' "+xvasp.Directory+"/OUTCAR");
                vector<string> tokens, vlines, vdata;
                tokens.clear(), vlines.clear(), vdata.clear();
                aurostd::string2tokens(stmp_toten, vlines, "\n");
                vector<double> vtoten; vtoten.clear();
                string line; line ="";
                for(uint i=0;i<vlines.size();i++){
                    line = vlines.at(i);
                    aurostd::string2tokens(line, vdata);
                    double toten = aurostd::string2utype<double>(vdata.at(4));
                    vtoten.push_back(toten);
                }
                int n = vtoten.size();
                double ediff = abs(vtoten.at(n-1) - vtoten.at(n-2));
                bool isCONVERGED = FALSE;
                if (ediff < 1.01E-4) {
                    isCONVERGED = TRUE;
                }

                if(nrun<maxrun) {  
                    if(LDEBUG) cerr << "KBIN::VASP_Run: " << Message("time") << "  checking warnings" << endl;
                    xmessage.flag("REACHED_ACCURACY",aurostd::substring_present_file_FAST(xvasp.Directory+"/vasp.out","reached required accuracy"));
                    xwarning.flag("KKSYM",aurostd::substring_present_file_FAST(xvasp.Directory+"/vasp.out","Reciprocal lattice and k-lattice belong to different class of lattices"));
                    xwarning.flag("SGRCON",aurostd::substring_present_file_FAST(xvasp.Directory+"/vasp.out","VERY BAD NEWS! internal error in subroutine SGRCON"));
                    xwarning.flag("NIRMAT",aurostd::substring_present_file_FAST(xvasp.Directory+"/vasp.out","Found some non-integer element in rotation matrix"));
                    xwarning.flag("BRMIX",aurostd::substring_present_file_FAST(xvasp.Directory+"/vasp.out","BRMIX: very serious problems"));
                    xwarning.flag("DAV",(!xmessage.flag("REACHED_ACCURACY") &&
                                aurostd::substring_present_file_FAST(xvasp.Directory+"/vasp.out","WARNING: Sub-Space-Matrix is not hermitian in DAV")) );
                    xwarning.flag("EDDDAV",aurostd::substring_present_file_FAST(xvasp.Directory+"/vasp.out","Error EDDDAV: Call to ZHEGV failed. Returncode"));
                    xwarning.flag("ZPOTRF",aurostd::substring_present_file_FAST(xvasp.Directory+"/vasp.out","LAPACK: Routine ZPOTRF failed"));
                    xwarning.flag("NBANDS",(aurostd::substring_present_file_FAST(xvasp.Directory+"/vasp.out","NBANDS") && // GET ALL TIMES
                                !aurostd::substring_present_file_FAST(xvasp.Directory+"/vasp.out","The number of bands has been changed from the values supplied") &&  // AVOID SELF HEALING
                                !aurostd::substring_present_file_FAST(xvasp.Directory+"/vasp.out","now  NBANDS  =")) );  // AVOID SELF HEALING
                    xwarning.flag("LRF_COMMUTATOR",aurostd::substring_present_file_FAST(xvasp.Directory+"/vasp.out","LRF_COMMUTATOR internal error: the vector")); // GET ALL TIMES
                    xwarning.flag("EXCCOR",aurostd::substring_present_file_FAST(xvasp.Directory+"/vasp.out","ERROR FEXCF: supplied exchange-correlation table")); // look for problem at the correlation
                    xwarning.flag("NATOMS",aurostd::substring_present_file_FAST(xvasp.Directory+"/vasp.out","The distance between some ions is very small")); // look for problem for distance
                    xwarning.flag("MEMORY",aurostd::substring_present_file_FAST(xvasp.Directory+"/vasp.out","AFLOW ERROR: AFLOW_MEMORY=")); // look for problem for distance
                    // xwarning.flag("PSMAXN",aurostd::substring_present_file_FAST(xvasp.Directory+"/vasp.out","WARNING: PSMAXN for non-local potential too small")); // look for problem for distance
                    xwarning.flag("PSMAXN",aurostd::substring_present_file_FAST(xvasp.Directory+"/vasp.out","REAL_OPT: internal ERROR"));
                    xwarning.flag("IBZKPT",(!xmessage.flag("REACHED_ACCURACY") &&
                                aurostd::substring_present_file_FAST(xvasp.Directory+"/vasp.out","VERY BAD NEWS! internal error in subroutine IBZKPT")) );
                    xwarning.flag("EDDRMM",(!xmessage.flag("REACHED_ACCURACY") &&
                                aurostd::substring_present_file_FAST(xvasp.Directory+"/vasp.out","WARNING in EDDRMM: call to ZHEGV failed, returncode")) ); // && !xwarning.flag("ZPOTRF");
                    xwarning.flag("REAL_OPTLAY_1",aurostd::substring_present_file_FAST(xvasp.Directory+"/vasp.out","REAL_OPTLAY: internal error (1)"));
                    xwarning.flag("REAL_OPT",aurostd::substring_present_file_FAST(xvasp.Directory+"/vasp.out","REAL_OPT: internal ERROR"));
                    xwarning.flag("SYMPREC",aurostd::substring_present_file_FAST(xvasp.Directory+"/vasp.out","inverse of rotation matrix was not found (increase SYMPREC)"));
                    xwarning.flag("SYMPREC",aurostd::substring_present_file_FAST(xvasp.Directory+"/vasp.out","VERY BAD NEWS! internal error in subroutine PRICEL"));
                    xwarning.flag("SYMPREC2",aurostd::substring_present_file_FAST(xvasp.Directory+"/vasp.out","Sorry, number of cells and number of vectors did not agree"));
                    xwarning.flag("INVGRP",aurostd::substring_present_file_FAST(xvasp.Directory+"/vasp.out","VERY BAD NEWS! internal error in subroutine INVGRP")); //KESONG, 2018-11-26
                    xwarning.flag("NPAR",aurostd::substring_present_file_FAST(xvasp.Directory+"/vasp.out","please rerun with NPAR=")); // not only npar==1
                    xwarning.flag("NPARC",(aurostd::substring_present_file_FAST(xvasp.Directory+"/vasp.out","NPAR = 4") &&
                                aurostd::substring_present_file_FAST(xvasp.Directory+"/vasp.out","NPAR=number of cores")) ); // fix with NPAR=cores in MPI
                    xwarning.flag("NPARN",(aurostd::substring_present_file_FAST(xvasp.Directory+"/vasp.out","NPAR = 4") &&
                                aurostd::substring_present_file_FAST(xvasp.Directory+"/vasp.out","NPAR=number of nodes")) ); // fix with NPAR=nodes in MPI
                    xwarning.flag("NPAR_REMOVE",aurostd::substring_present_file_FAST(xvasp.Directory+"/vasp.out","Please remove the tag NPAR from the INCAR file and restart the"));
                    xwarning.flag("GAMMA_SHIFT",aurostd::substring_present_file_FAST(xvasp.Directory+"/vasp.out","shift your grid to Gamma"));
                    xwarning.flag("CSLOSHING",KBIN::VASP_CheckUnconvergedOSZICAR(xvasp.Directory)); // check from OSZICAR
                    xwarning.flag("NKXYZ_IKPTD",(aurostd::substring_present_file_FAST(xvasp.Directory+"/vasp.out","NKX>IKPTD") ||
                                aurostd::substring_present_file_FAST(xvasp.Directory+"/vasp.out","NKY>IKPTD") ||
                                aurostd::substring_present_file_FAST(xvasp.Directory+"/vasp.out","NKZ>IKPTD")) );
                    xwarning.flag("DENTET",aurostd::substring_present_file_FAST(xvasp.Directory+"/vasp.out","WARNING: DENTET: can't reach specified precision")); // not only npar==1
                    xwarning.flag("EFIELD_PEAD",aurostd::substring_present_file_FAST(xvasp.Directory+"/vasp.out","EFIELD_PEAD is too large")); // EFIELD_PEAD
                    xwarning.flag("READ_KPOINTS_RD_SYM",(aurostd::substring_present_file_FAST(xvasp.Directory+"/vasp.out","ERROR in RE_READ_KPOINTS_RD") &&
                                aurostd::substring_present_file_FAST(xvasp.Directory+"/vasp.out","switch off symmetry")) );
                    xwarning.flag("MPICH11",(aurostd::substring_present_file_FAST(xvasp.Directory+"/vasp.out","BAD TERMINATION OF ONE OF YOUR APPLICATION PROCESSES") &&
                                aurostd::substring_present_file_FAST(xvasp.Directory+"/vasp.out","EXIT CODE: 11")) );
                    xwarning.flag("MPICH139",(aurostd::substring_present_file_FAST(xvasp.Directory+"/vasp.out","BAD TERMINATION OF ONE OF YOUR APPLICATION PROCESSES") &&
                                aurostd::substring_present_file_FAST(xvasp.Directory+"/vasp.out","EXIT CODE: 139")) );
                    if(xwarning.flag("MPICH11")) xwarning.flag("NBANDS",FALSE); // fix MPICH11 first

                    if(xwarning.flag("NPARC") && (aurostd::substring_present_file_FAST(xvasp.Directory+"/INCAR","NPAR=2") || // dont bother for small NPAR
                                aurostd::substring_present_file_FAST(xvasp.Directory+"/INCAR","LRPA") ||
                                aurostd::substring_present_file_FAST(xvasp.Directory+"/INCAR","LEPSILON") ||
                                aurostd::substring_present_file_FAST(xvasp.Directory+"/INCAR","LOPTICS"))) xwarning.flag("NPARC",FALSE);  // dont touch NPARC if LRPA or LEPSILON or LOPTICS necessary	
                    if(xwarning.flag("NPARN") && (aurostd::substring_present_file_FAST(xvasp.Directory+"/INCAR","LRPA") ||
                                aurostd::substring_present_file_FAST(xvasp.Directory+"/INCAR","LEPSILON") ||
                                aurostd::substring_present_file_FAST(xvasp.Directory+"/INCAR","LOPTICS"))) xwarning.flag("NPARN",FALSE);  // dont touch NPARN if LRPA or LEPSILON or LOPTICS necessary

                    int NBANDS_OUTCAR=0;
                    xOUTCAR OUTCAR_NBANDS(xvasp.Directory+"/OUTCAR");
                    NBANDS_OUTCAR=OUTCAR_NBANDS.NBANDS;
                    if(xwarning.flag("NBANDS") && NBANDS_OUTCAR>1000) xwarning.flag("NBANDS",FALSE); // for safety
                    if(xwarning.flag("NBANDS") && aurostd::substring_present_file_FAST(xvasp.Directory+"/INCAR","DIELECTRIC_STATIC") && NBANDS_OUTCAR>1000) xwarning.flag("NBANDS",FALSE); // for safety
                    if(xwarning.flag("NBANDS") && aurostd::substring_present_file_FAST(xvasp.Directory+"/INCAR","DIELECTRIC_DYNAMIC") && NBANDS_OUTCAR>1000) xwarning.flag("NBANDS",FALSE); // for safety
                    if(xwarning.flag("NBANDS") && aurostd::substring_present_file_FAST(xvasp.Directory+"/INCAR","DSCF") && NBANDS_OUTCAR>1000) xwarning.flag("NBANDS",FALSE); // for safety

                    // fix troubles
                    if(xmessage.flag("REACHED_ACCURACY") && xwarning.flag("IBZKPT")) xwarning.flag("IBZKPT",FALSE);  // priority
                    if(xwarning.flag("NKXYZ_IKPTD")) xwarning.flag("IBZKPT",FALSE); // priority
                    if (isCONVERGED && xwarning.flag("IBZKPT")) xwarning.flag("IBZKPT",FALSE);

                    xfixed.flag("ALL",FALSE);
                    xfixed.flag("MPICH11",FALSE); // all the items that must be restarted until they work
                    xfixed.flag("MPICH139",FALSE); // all the items that must be restarted until they work
                    vasp_start=FALSE;

                    // ********* CHECK NBANDS PROBLEMS ******************
                    if(LDEBUG) cerr << "KBIN::VASP_Run: " << Message("time") << "  [CHECK NBANDS PROBLEMS]" << endl;
                    if(!vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("NBANDS") && !xfixed.flag("ALL")) { // check NBANDS
                        if(xwarning.flag("NBANDS")) {
                            KBIN::VASP_Error(xvasp,"WWWWW  ERROR KBIN::VASP_Run: "+Message("time")+"  NBANDS problems ");
                            KBIN::XVASP_Afix_NBANDS(xvasp,nbands,!XHOST.QUIET);  // here it does the nbands_update
                            xfixed.flag("NBANDS",TRUE);xfixed.flag("ALL",TRUE);
                            aus << "WWWWW  FIX NBANDS = [" << nbands << "] - " << Message(aflags,"user,host,time") << endl;
                            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                        }
                    }
                    // ********* CHECK LRF_COMMUTATOR PROBLEMS ******************
                    if(LDEBUG) cerr << "KBIN::VASP_Run: " << Message("time") << "  [CHECK LRF_COMMUTATOR PROBLEMS]" << endl;
                    if(!vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("LRF_COMMUTATOR") && !xfixed.flag("ALL")) { // check LRF_COMMUTATOR
                        if(0 && xwarning.flag("LRF_COMMUTATOR")) {
                            KBIN::VASP_Error(xvasp,"WWWWW  ERROR KBIN::VASP_Run: "+Message("time")+"  LRF_COMMUTATOR problems ");
                            //	  KBIN_XVASP_Afix_LRF_COMMUTATOR(xvasp,nbands);  // here it does the nbands_update
                            xfixed.flag("LRF_COMMUTATOR",TRUE);xfixed.flag("ALL",TRUE);
                            aus << "WWWWW  FIX LRF_COMMUTATOR = [" << nbands << "] - " << Message(aflags,"user,host,time") << endl;
                            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                        }
                    }
                    // // ********* CHECK SGRCON AND NIRMAT PROBLEMS ******************
                    //	if(LDEBUG) cerr << "KBIN::VASP_Run: " << Message("time") << "  [CHECK SGRCON AND NIRMAT PROBLEMS]" << endl;
                    // if(!vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("SGRCON") && !xfixed.flag("ALL")) { // OPTIONS FOR SYMMETRY
                    // 	if(xwarning.flag("SGRCON") && xwarning.flag("NIRMAT")) {
                    // 	  KBIN::VASP_Error(xvasp,"WWWWW  ERROR KBIN::VASP_Run: "+Message("time")+"  SGRCON/NIRMAT problems ");
                    // 	  aus << "WWWWW  FIX SGRCON/NIRMAT - " << Message(aflags,"user,host,time") << endl;
                    // 	  aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                    // 	  KBIN::XVASP_Afix_GENERIC("SGRCON/NIRMAT",xvasp,kflags,vflags);
                    // 	  xfixed.flag("SGRCON",TRUE);xfixed.flag("ALL",TRUE);
                    // 	  xfixed.flag("NIRMAT",TRUE);xfixed.flag("ALL",TRUE);
                    // 	  // if(nrun<maxrun) vasp_start=TRUE;
                    // 	}
                    // }
                    // ********* CHECK SYMMETRY PROBLEMS ******************
                    if(LDEBUG) cerr << "KBIN::VASP_Run: " << Message("time") << "  [CHECK SYMMETRY PROBLEMS]" << endl;
                    if(!vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("ROTMAT") && !xfixed.flag("ALL")) { // OPTIONS FOR SYMMETRY
                        if(xwarning.flag("KKSYM") || xwarning.flag("SGRCON") || xwarning.flag("NIRMAT")) {
                            if(xwarning.flag("NIRMAT")) {
                                fix_NIRMAT++;
                                xvasp.str.kpoints_k1=kpoints_k1;xvasp.str.kpoints_s1=kpoints_s1;
                                xvasp.str.kpoints_k2=kpoints_k2;xvasp.str.kpoints_s2=kpoints_s2;
                                xvasp.str.kpoints_k3=kpoints_k3;xvasp.str.kpoints_s3=kpoints_s3;
                                KBIN::VASP_Error(xvasp,"WWWWW  ERROR KBIN::VASP_Run: "+Message("time")+"  NIRMAT problems ");
                                if(fix_NIRMAT<=6) {
                                    aus << "WWWWW  FIX NIRMAT (" << fix_NIRMAT << ")- " << Message(aflags,"user,host,time") << endl;
                                    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                                    KBIN::XVASP_Afix_ROTMAT(xvasp,fix_NIRMAT,!XHOST.QUIET,aflags,FileMESSAGE);
                                    xfixed.flag("ALL",TRUE);
                                } else {
                                    aus << "WWWWW  FIX NIRMAT (" << 0 << ")- " << Message(aflags,"user,host,time") << endl;
                                    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                                    KBIN::XVASP_Afix_ROTMAT(xvasp,0,!XHOST.QUIET,aflags,FileMESSAGE);
                                    xfixed.flag("KKSYM",TRUE);xfixed.flag("NIRMAT",TRUE);xfixed.flag("ALL",TRUE);
                                    // if(nrun<maxrun) vasp_start=TRUE;
                                }
                            } else { // JAN 2012
                                /*
                                //NOT WORKING
                                KBIN::VASP_Error(xvasp,"WWWWW  ERROR KBIN::VASP_Run: "+Message("time")+"  \"Reciprocal lattice and k-lattice belong to different class of lattices\" problem");
                                aus << "WWWWW  FIX K1=K2=K3=KMAX - " << Message(aflags,"user,host,time") << endl;
                                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                                xvasp.str.kpoints_k1=kpoints_k1;xvasp.str.kpoints_s1=kpoints_s1;
                                xvasp.str.kpoints_k2=kpoints_k2;xvasp.str.kpoints_s2=kpoints_s2;
                                xvasp.str.kpoints_k3=kpoints_k3;xvasp.str.kpoints_s3=kpoints_s3;
                                KBIN::XVASP_Afix_ROTMAT(xvasp,3,!XHOST.QUIET,aflags,FileMESSAGE);
                                xfixed.flag("KKSYM",TRUE);xfixed.flag("ALL",TRUE);
                                */
                            }
                        }
                    }
                    // ********* CHECK SGRCON PROBLEMS ******************
                    if(LDEBUG) cerr << "KBIN::VASP_Run: " << Message("time") << "  [CHECK SGRCON PROBLEMS]" << endl;
                    if(!vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("SGRCON") && !xfixed.flag("ALL")) { // OPTIONS FOR SYMMETRY
                        if(xwarning.flag("SGRCON") && !xfixed.flag("SGRCON")) {
                            KBIN::VASP_Error(xvasp,"WWWWW  ERROR KBIN::VASP_Run: "+Message("time")+"  SGRCON problems ");
                            aus << "WWWWW  FIX SGRCON - " << Message(aflags,"user,host,time") << endl;
                            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                            KBIN::XVASP_Afix_GENERIC("SGRCON",xvasp,kflags,vflags);
                            xfixed.flag("SGRCON",TRUE);xfixed.flag("ALL",TRUE);
                            // if(nrun<maxrun) vasp_start=TRUE;
                        }
                    }
                    // ********* CHECK GAMMA_SHIFT PROBLEMS ******************
                    if(LDEBUG) cerr << "KBIN::VASP_Run: " << Message("time") << "  [CHECK GAMMA_SHIFT PROBLEMS]" << endl;
                    if(!vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("GAMMA_SHIFT") && !xfixed.flag("ALL")) {
                        if(xwarning.flag("GAMMA_SHIFT") && !xfixed.flag("GAMMA_SHIFT")) {
                            KBIN::VASP_Error(xvasp,"WWWWW  ERROR KBIN::VASP_Run: "+Message("time")+"  GAMMA_SHIFT problems ");
                            aus << "WWWWW  FIX GAMMA_SHIFT - " << Message(aflags,"user,host,time") << endl;
                            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                            KBIN::XVASP_Afix_GENERIC("GAMMA_SHIFT",xvasp,kflags,vflags);
                            xfixed.flag("GAMMA_SHIFT",TRUE);xfixed.flag("ALL",TRUE);
                            // if(nrun<maxrun) vasp_start=TRUE;
                        }
                    }
                    // ********* CHECK MPICH11 PROBLEMS ******************
                    if(LDEBUG) cerr << "KBIN::VASP_Run: " << Message("time") << "  [CHECK MPICH11 PROBLEMS]" << endl;
                    if(!vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("MPICH11") && !xfixed.flag("ALL")) {
                        if(xwarning.flag("MPICH11") && !xfixed.flag("MPICH11")) {
                            KBIN::VASP_Error(xvasp,"WWWWW  ERROR KBIN::VASP_Run: "+Message("time")+"  MPICH11 problems ");
                            aus << "WWWWW  FIX MPICH11 - " << Message(aflags,"user,host,time") << endl;
                            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                            KBIN::XVASP_Afix_GENERIC("MPICH11",xvasp,kflags,vflags);
                            xfixed.flag("MPICH11",TRUE);xfixed.flag("ALL",TRUE);
                            // if(nrun<maxrun) vasp_start=TRUE;
                        }
                    }
                    if(!vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("MPICH139") && !xfixed.flag("ALL")) {
                        if(xwarning.flag("MPICH139") && !xfixed.flag("MPICH139")) {
                            KBIN::VASP_Error(xvasp,"WWWWW  ERROR KBIN::VASP_Run: "+Message("time")+"  MPICH139 problems ");
                            aus << "WWWWW  FIX MPICH139 - " << Message(aflags,"user,host,time") << endl;
                            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                            KBIN::XVASP_Afix_GENERIC("MPICH139",xvasp,kflags,vflags);
                            xfixed.flag("MPICH139",TRUE);xfixed.flag("ALL",TRUE);
                            // if(nrun<maxrun) vasp_start=TRUE;
                        }
                    }
                    // ********* CHECK IBZKPT PROBLEMS ******************
                    if(LDEBUG) cerr << "KBIN::VASP_Run: " << Message("time") << "  [CHECK IBZKPT PROBLEMS]" << endl;
                    if(!vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("IBZKPT") && !xfixed.flag("ALL")) { // OPTIONS FOR SYMMETRY
                        if(xwarning.flag("IBZKPT") && !xfixed.flag("IBZKPT")) {
                            KBIN::VASP_Error(xvasp,"WWWWW  ERROR KBIN::VASP_Run: "+Message("time")+"  IBZKPT problems ");
                            aus << "WWWWW  FIX IBZKPT - " << Message(aflags,"user,host,time") << endl;
                            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                            KBIN::XVASP_Afix_GENERIC("IBZKPT",xvasp,kflags,vflags);
                            xfixed.flag("IBZKPT",TRUE);xfixed.flag("ALL",TRUE);
                            // if(nrun<maxrun) vasp_start=TRUE;
                        }
                    }
                    // ********* CHECK NKXYZ_IKPTD PROBLEMS ******************
                    if(LDEBUG) cerr << "KBIN::VASP_Run: " << Message("time") << "  [CHECK NKXYZ_IKPTD PROBLEMS]" << endl;
                    if(!vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("NKXYZ_IKPTD") && !xfixed.flag("ALL")) { // OPTIONS FOR SYMMETRY
                        if(xwarning.flag("NKXYZ_IKPTD") && !xfixed.flag("NKXYZ_IKPTD")) {
                            KBIN::VASP_Error(xvasp,"WWWWW  ERROR KBIN::VASP_Run: "+Message("time")+"  NKXYZ_IKPTD problems ");
                            aus << "WWWWW  FIX NKXYZ_IKPTD - " << Message(aflags,"user,host,time") << endl;
                            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                            KBIN::XVASP_Afix_GENERIC("NKXYZ_IKPTD",xvasp,kflags,vflags);
                            xfixed.flag("NKXYZ_IKPTD",TRUE);xfixed.flag("ALL",TRUE);
                            // if(nrun<maxrun) vasp_start=TRUE;
                        }
                    }
                    // ********* CHECK SYMPREC PROBLEMS ******************
                    if(LDEBUG) cerr << "KBIN::VASP_Run: " << Message("time") << "  [CHECK SYMPREC PROBLEMS]" << endl;
                    if(!vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("SYMPREC") && !xfixed.flag("ALL")) { // OPTIONS FOR SYMMETRY
                        if(xwarning.flag("SYMPREC") && !xfixed.flag("SYMPREC")) {
                            KBIN::VASP_Error(xvasp,"WWWWW  ERROR KBIN::VASP_Run: "+Message("time")+"  SYMPREC problems ");
                            aus << "WWWWW  FIX SYMPREC - " << Message(aflags,"user,host,time") << endl;
                            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                            KBIN::XVASP_Afix_GENERIC("SYMPREC",xvasp,kflags,vflags);
                            xfixed.flag("SYMPREC",TRUE);xfixed.flag("ALL",TRUE);
                            // if(nrun<maxrun) vasp_start=TRUE;
                        }
                    }
                    // ********* CHECK SYMPREC2 PROBLEMS ******************
                    if(LDEBUG) cerr << "KBIN::VASP_Run: " << Message("time") << "  [CHECK SYMPREC PROBLEMS]" << endl;
                    if(!vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("SYMPREC2") && !xfixed.flag("ALL")) { // OPTIONS FOR SYMMETRY
                        if(xwarning.flag("SYMPREC2") && !xfixed.flag("SYMPREC2")) {
                            KBIN::VASP_Error(xvasp,"WWWWW  ERROR KBIN::VASP_Run: "+Message("time")+"  SYMPREC2 problems ");
                            aus << "WWWWW  FIX SYMPREC2 - " << Message(aflags,"user,host,time") << endl;
                            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                            KBIN::XVASP_Afix_GENERIC("SYMPREC2",xvasp,kflags,vflags);
                            xfixed.flag("SYMPREC2",TRUE);xfixed.flag("ALL",TRUE);
                            // if(nrun<maxrun) vasp_start=TRUE;
                        }
                    }
                    // ********* CHECK INVGRP PROBLEMS ******************
                    if(LDEBUG) cerr << "KBIN::VASP_Run: " << Message("time") << "  [CHECK INVGRP PROBLEMS]" << endl;
                    if(!vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("INVGRP") && !xfixed.flag("ALL")) { // OPTIONS FOR SYMMETRY
                        if(xwarning.flag("INVGRP") && !xfixed.flag("INVGRP")) {
                            KBIN::VASP_Error(xvasp,"WWWWW  ERROR KBIN::VASP_Run: "+Message("time")+"  INVGRP problems ");
                            aus << "WWWWW  FIX INVGRP - " << Message(aflags,"user,host,time") << endl;
                            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                            KBIN::XVASP_Afix_GENERIC("INVGRP",xvasp,kflags,vflags);
                            xfixed.flag("INVGRP",TRUE);xfixed.flag("ALL",TRUE);
                            // if(nrun<maxrun) vasp_start=TRUE;
                        }
                    }
                    // ********* CHECK EDDRMM PROBLEMS ******************
                    if(LDEBUG) cerr << "KBIN::VASP_Run: " << Message("time") << "  [CHECK EDDRMM PROBLEMS]" << endl;
                    if(!vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("EDDRMM") && !xfixed.flag("ALL")) { // OPTIONS FOR EDDRMM
                        if(xwarning.flag("EDDRMM") && !xfixed.flag("EDDRMM")) {
                            KBIN::VASP_Error(xvasp,"WWWWW  ERROR KBIN::VASP_Run: "+Message("time")+"  EDDRMM problems ");
                            aus << "WWWWW  FIX EDDRMM - " << Message(aflags,"user,host,time") << endl;
                            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                            KBIN::XVASP_Afix_GENERIC("EDDRMM",xvasp,kflags,vflags);
                            xfixed.flag("EDDRMM",TRUE);xfixed.flag("ALL",TRUE);
                            // if(nrun<maxrun) vasp_start=TRUE;
                        }
                    }
                    // ********* CHECK REAL_OPTLAY_1 REAL_OPT PROBLEMS ******************
                    if(LDEBUG) cerr << "KBIN::VASP_Run: " << Message("time") << "  [CHECK REAL_OPTLAY_1 REAL_OPT PROBLEMS]" << endl;
                    if(!vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("LREAL") && !xfixed.flag("ALL")) { // OPTIONS FOR REAL_OPTLAY_1 OPTLAY
                        if(xwarning.flag("REAL_OPTLAY_1") && !xfixed.flag("REAL_OPTLAY_1")) {
                            KBIN::VASP_Error(xvasp,"WWWWW  ERROR KBIN::VASP_Run: "+Message("time")+"  REAL_OPTLAY_1 problems ");
                            aus << "WWWWW  FIX REAL_OPTLAY_1 - " << Message(aflags,"user,host,time") << endl;
                            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                            KBIN::XVASP_Afix_GENERIC("LREAL",xvasp,kflags,vflags);
                            xfixed.flag("REAL_OPTLAY_1",TRUE);xfixed.flag("ALL",TRUE);
                            // if(nrun<maxrun) vasp_start=TRUE;
                        }
                        if(xwarning.flag("REAL_OPT") && !xfixed.flag("REAL_OPT")) {
                            KBIN::VASP_Error(xvasp,"WWWWW  ERROR KBIN::VASP_Run: "+Message("time")+"  REAL_OPT problems ");
                            aus << "WWWWW  FIX REAL_OPT - " << Message(aflags,"user,host,time") << endl;
                            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                            KBIN::XVASP_Afix_GENERIC("LREAL",xvasp,kflags,vflags);
                            xfixed.flag("REAL_OPT",TRUE);xfixed.flag("ALL",TRUE);
                            // if(nrun<maxrun) vasp_start=TRUE;
                        }
                    }
                    // ********* CHECK BRMIX PROBLEMS ******************
                    if(LDEBUG) cerr << "KBIN::VASP_Run: " << Message("time") << "  [CHECK BRMIX PROBLEMS]" << endl;
                    if(!vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("BRMIX") && !xfixed.flag("ALL")) { // check BRMIX
                        if(xwarning.flag("BRMIX") && !xfixed.flag("BRMIX")) { // Apply only ONCE
                            KBIN::VASP_Error(xvasp,"WWWWW  ERROR KBIN::VASP_Run: "+Message("time")+"  BRMIX problems ");
                            aus << "WWWWW  FIX BRMIX - " << Message(aflags,"user,host,time") << endl;
                            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                            /*	  if(aurostd::substring2bool(xvasp.INCAR,"IALGO=48",TRUE)) {
                                  KBIN::VASP_Error(xvasp,"WWWWW  ERROR KBIN::VASP_Run: "+Message("time")+"  BRMIX problems not performed with IALGO=48 ");
                                  aus << "WWWWW  FIX BRMIX NOT PERFORMED with IALGO=48 - " << Message(aflags,"user,host,time") << endl;
                                  aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                                  xfixed.flag("BRMIX",TRUE);xfixed.flag("ALL",TRUE);
                                  } else {	
                                  */
                            KBIN::XVASP_Afix_GENERIC("BRMIX",xvasp,kflags,vflags);
                            xfixed.flag("BRMIX",TRUE);xfixed.flag("ALL",TRUE);
                        }
                        }
                        // ********* CHECK DAV PROBLEMS ******************
                        if(LDEBUG) cerr << "KBIN::VASP_Run: " << Message("time") << "  [CHECK DAV PROBLEMS]" << endl;
                        if(!vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("DAV") && !xfixed.flag("ALL")) { // check DAV
                            if(xwarning.flag("DAV") && !xfixed.flag("DAV")) { // Apply only ONCE
                                KBIN::VASP_Error(xvasp,"WWWWW  ERROR KBIN::VASP_Run: "+Message("time")+"  DAV problems ");
                                aus << "WWWWW  FIX DAV - " << Message(aflags,"user,host,time") << endl;
                                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                                KBIN::XVASP_Afix_GENERIC("DAV",xvasp,kflags,vflags);
                                xfixed.flag("DAV",TRUE);xfixed.flag("ALL",TRUE);
                            }
                        }
                        // ********* CHECK EDDDAV PROBLEMS ******************
                        if(LDEBUG) cerr << "KBIN::VASP_Run: " << Message("time") << "  [CHECK EDDDAV PROBLEMS]" << endl;
                        if(!vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("EDDDAV") && !xfixed.flag("ALL")) { // check EDDDAV
                            if(xwarning.flag("EDDDAV") && !xfixed.flag("EDDDAV")) { // Apply only ONCE
                                KBIN::VASP_Error(xvasp,"WWWWW  ERROR KBIN::VASP_Run: "+Message("time")+"  EDDDAV problems ");
                                aus << "WWWWW  FIX EDDDAV - " << Message(aflags,"user,host,time") << endl;
                                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                                KBIN::XVASP_Afix_GENERIC("EDDDAV",xvasp,kflags,vflags);
                                xfixed.flag("EDDDAV",TRUE);xfixed.flag("ALL",TRUE);
                            }
                        }
                        // ********* CHECK EFIELD_PEAD PROBLEMS ******************
                        if(LDEBUG) cerr << "KBIN::VASP_Run: " << Message("time") << "  [CHECK EFIELD_PEAD PROBLEMS]" << endl;
                        if(!vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("EFIELD_PEAD") && !xfixed.flag("ALL")) { // check EFIELD_PEAD
                            if(xwarning.flag("EFIELD_PEAD")) { // can be applied many times && !xfixed.flag("EFIELD_PEAD")) { // Apply only ONCE
                                KBIN::VASP_Error(xvasp,"WWWWW  ERROR KBIN::VASP_Run: "+Message("time")+"  EFIELD_PEAD problems ");
                                aus << "WWWWW  FIX EFIELD_PEAD - " << Message(aflags,"user,host,time") << endl;
                                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                                KBIN::XVASP_Afix_GENERIC("EFIELD_PEAD",xvasp,kflags,vflags);
                                xfixed.flag("EFIELD_PEAD",TRUE);xfixed.flag("ALL",TRUE);
                            }
                            }
                            // ********* CHECK ZPOTRF PROBLEMS ******************
                            if(LDEBUG) cerr << "KBIN::VASP_Run: " << Message("time") << "  [CHECK ZPOTRF PROBLEMS]" << endl;
                            if(0) if(!vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("ZPOTRF") && !xfixed.flag("ALL")) { // check ZPOTRF
                                if(xwarning.flag("ZPOTRF") && !xfixed.flag("ZPOTRF")) { // Apply only ONCE
                                    counter_ZPOTRF++;
                                    KBIN::VASP_Error(xvasp,"WWWWW  ERROR KBIN::VASP_Run: "+Message("time")+"  ZPOTRF problems ");
                                    aus << "WWWWW  FIX ZPOTRF - " << Message(aflags,"user,host,time") << endl;
                                    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                                    aus << "00000  MESSAGE WARNING SWAPPING TO CONVENTIONAL STRUCTURE" << " - " << Message(aflags,"user,host,time") << endl;
                                    aus << "00000  MESSAGE BEFORE: a,b,c,alpha,beta,gamma " << xvasp.str.a << "," << xvasp.str.b << "," << xvasp.str.c << "," << xvasp.str.alpha << "," << xvasp.str.beta << "," << xvasp.str.gamma << endl;
                                    aus << "00000  MESSAGE BEFORE: structure: " << endl;
                                    aus << xvasp.str;
                                    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                                    KBIN::XVASP_Afix_GENERIC("ZPOTRF",xvasp,kflags,vflags);
                                    aus << "00000  MESSAGE AFTER: a,b,c,alpha,beta,gamma " << xvasp.str.a << "," << xvasp.str.b << "," << xvasp.str.c << "," << xvasp.str.alpha << "," << xvasp.str.beta << "," << xvasp.str.gamma << endl;
                                    aus << "00000  MESSAGE AFTER: structure: " << endl;
                                    aus << xvasp.str;
                                    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                                    xfixed.flag("ZPOTRF",TRUE);xfixed.flag("ALL",TRUE);
                                }
                            }
                            if(!vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("ZPOTRF") && !xfixed.flag("ALL")) { // check ZPOTRF
                                if(xwarning.flag("ZPOTRF") && !xfixed.flag("ZPOTRF")) { // Apply only ONCE
                                    counter_ZPOTRF++;
                                    KBIN::VASP_Error(xvasp,"WWWWW  ERROR KBIN::VASP_Run: "+Message("time")+"  ZPOTRF problems ");
                                    aus << "WWWWW  FIX ZPOTRF - " << Message(aflags,"user,host,time") << endl;
                                    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);	    
                                    aus << "00000  MESSAGE WARNING CHANGING POTIM" << " - " << Message(aflags,"user,host,time") << endl;
                                    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                                    KBIN::XVASP_Afix_GENERIC("ZPOTRF_POTIM",xvasp,kflags,vflags);
                                    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                                    xfixed.flag("ZPOTRF",TRUE);xfixed.flag("ALL",TRUE);
                                }
                            }

                            // ********* CHECK EXCHANGE_CORRELATION PROBLEMS ******************
                            if(LDEBUG) cerr << "KBIN::VASP_Run: " << Message("time") << "  [CHECK EXCHANGE_CORRELATION PROBLEMS]" << endl;
                            if(!vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("EXCCOR") && !xfixed.flag("ALL")) { // OPTIONS FOR EXCCOR
                                if(xwarning.flag("EXCCOR") && !xfixed.flag("EXCCOR")) {  // Apply only ONCE
                                    KBIN::VASP_Error(xvasp,"WWWWW  ERROR KBIN::VASP_Run: "+Message("time")+"  EXCHANGE-CORRELATION problems ");
                                    aus << "WWWWW  FIX  EXCHANGE-CORRELATION - " << Message(aflags,"user,host,time") << endl;
                                    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                                    KBIN::XVASP_Afix_GENERIC("EXCCOR",xvasp,kflags,vflags);
                                    xfixed.flag("EXCCOR",TRUE);xfixed.flag("ALL",TRUE);
                                }
                            }
                            // ********* NEAREST NEIGHBOUR ATOMS PROBLEMS ******************
                            if(LDEBUG) cerr << "KBIN::VASP_Run: " << Message("time") << "  [CHECK NEAREST NEIGHBOUR ATOMS PROBLEMS]" << endl;
                            if(!vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("NATOMS") && !xfixed.flag("ALL")) { // OPTIONS FOR NATOMS
                                if(xwarning.flag("NATOMS") && !xfixed.flag("NATOMS")) {  // Apply only ONCE
                                    KBIN::VASP_Error(xvasp,"WWWWW  ERROR KBIN::VASP_Run: "+Message("time")+"  NEAREST-NEIGHBOURS ATOMS problems ");
                                    aus << "WWWWW  FIX  NEAREST-NEIGHBOURS ATOMS - " << Message(aflags,"user,host,time") << endl;
                                    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                                    KBIN::XVASP_Afix_GENERIC("NATOMS",xvasp,kflags,vflags);
                                    xfixed.flag("NATOMS",TRUE);xfixed.flag("ALL",TRUE);
                                }
                            }
                            // ********* MEMORY PROBLEMS ******************
                            if(LDEBUG) cerr << "KBIN::VASP_Run: " << Message("time") << "  [CHECK MEMORY PROBLEMS]" << endl;
                            if(!vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("MEMORY") && !xfixed.flag("ALL")) { // OPTIONS FOR MEMORY
                                if(xwarning.flag("MEMORY") && !xfixed.flag("MEMORY")) {  // Apply only ONCE
                                    KBIN::VASP_Error(xvasp,"WWWWW  ERROR KBIN::VASP_Run: "+Message("time")+"  MEMORY problems ");
                                    aus << "WWWWW  FIX  MEMORY - " << Message(aflags,"user,host,time") << endl;
                                    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                                    stringstream command("");
                                    command << "cat " << xvasp.Directory << "/vasp.out | grep AFLOW > " << xvasp.Directory << "/" << DEFAULT_AFLOW_MEMORY_OUT << endl;
                                    command << "cat " << xvasp.Directory << "/vasp.out | grep AFLOW > " << xvasp.Directory << "/SKIP" << endl;	
                                    command << "cat " << xvasp.Directory << "/vasp.out | grep AFLOW >> " << xvasp.Directory << DEFAULT_AFLOW_ERVASP_OUT << endl;	
                                    aurostd::execute(command);
                                    xfixed.flag("MEMORY",TRUE);xfixed.flag("ALL",TRUE);
                                    KBIN::VASP_Error(xvasp,FileMESSAGE,"EEEEE  ERROR MEMORY PROBLEMS "); return FALSE;
                                }
                            }
                            // ********* CHECK PSMAXN PROBLEMS ******************
                            if(LDEBUG) cerr << "KBIN::VASP_Run: " << Message("time") << "  [CHECK PSMAXN PROBLEMS]" << endl;
                            if(!vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("PSMAXN") && !xfixed.flag("ALL")) { // CHECK FOR PSMAXN
                                if(xwarning.flag("PSMAXN") && !xfixed.flag("PSMAXN")) {  // Apply only ONCE
                                    KBIN::VASP_Error(xvasp,"WWWWW  ERROR KBIN::VASP_Run: "+Message("time")+"  PSMAXN problems ");
                                    aus << "WWWWW  FIX PSMAXN - " << Message(aflags,"user,host,time") << endl;
                                    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                                    enmax=KBIN::XVASP_Afix_GENERIC("PSMAXN",xvasp,kflags,vflags,enmax);
                                    xfixed.flag("PSMAXN",TRUE);xfixed.flag("ALL",TRUE);
                                }
                            }
                            // ********* CHECK NPAR PROBLEMS ******************
                            if(LDEBUG) cerr << "KBIN::VASP_Run: " << Message("time") << "  [CHECK NPAR PROBLEM]" << endl;
                            if(!vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("NPAR") && !xfixed.flag("ALL")) { // CHECK FOR NPAR
                                if(xwarning.flag("NPAR") && !xfixed.flag("NPAR")) {  // Apply only ONCE
                                    KBIN::VASP_Error(xvasp,"WWWWW  ERROR KBIN::VASP_Run: "+Message("time")+"  NPAR problems ");
                                    aus << "WWWWW  FIX NPAR - " << Message(aflags,"user,host,time") << endl;
                                    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                                    KBIN::XVASP_Afix_GENERIC("NPAR",xvasp,kflags,vflags);
                                    xfixed.flag("NPAR",TRUE);xfixed.flag("ALL",TRUE);
                                }
                            }
                            // ********* CHECK NPARC PROBLEMS ******************
                            if(LDEBUG) cerr << "KBIN::VASP_Run: " << Message("time") << "  [CHECK NPARC PROBLEMS]" << endl;
                            if(!vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("NPARC") && !xfixed.flag("ALL")) { // CHECK FOR NPARC
                                if(xwarning.flag("NPARC") && !xfixed.flag("NPARC")) {  // Apply only ONCE
                                    KBIN::VASP_Error(xvasp,"WWWWW  ERROR KBIN::VASP_Run: "+Message("time")+"  NPARC problems ");
                                    aus << "WWWWW  FIX NPARC - " << Message(aflags,"user,host,time") << endl;
                                    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                                    KBIN::XVASP_Afix_GENERIC("NPARC",xvasp,kflags,vflags,0.0,kflags.KBIN_MPI_NCPUS);  // this will kill the system as NBANDS_MPI=NPAR*NBANDS_SERIAL
                                    xfixed.flag("NPARC",TRUE);xfixed.flag("ALL",TRUE);
                                }
                            }
                            // ********* CHECK NPARN PROBLEMS ******************
                            if(LDEBUG) cerr << "KBIN::VASP_Run: " << Message("time") << "  [CHECK NPARN PROBLEMS]" << endl;
                            if(!vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("NPARN") && !xfixed.flag("ALL")) { // CHECK FOR NPARN
                                if(xwarning.flag("NPARN") && !xfixed.flag("NPARN")) {  // Apply only ONCE
                                    KBIN::VASP_Error(xvasp,"WWWWW  ERROR KBIN::VASP_Run: "+Message("time")+"  NPARN problems ");
                                    aus << "WWWWW  FIX NPARN - " << Message(aflags,"user,host,time") << endl;
                                    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                                    KBIN::XVASP_Afix_GENERIC("NPARN",xvasp,kflags,vflags,0.0,kflags.KBIN_MPI_NCPUS);  // this will kill the system as NBANDS_MPI=NPAR*NBANDS_SERIAL
                                    xfixed.flag("NPARN",TRUE);xfixed.flag("ALL",TRUE);
                                }
                            }
                            // ********* CHECK NPAR_REMOVE PROBLEMS ******************
                            if(LDEBUG) cerr << "KBIN::VASP_Run: " << Message("time") << "  [CHECK NPAR_REMOVE PROBLEMS]" << endl;
                            if(!vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("NPAR_REMOVE") && !xfixed.flag("ALL")) { // CHECK FOR NPAR_REMOVE
                                if(xwarning.flag("NPAR_REMOVE") && !xfixed.flag("NPAR_REMOVE")) {  // Apply only ONCE
                                    KBIN::VASP_Error(xvasp,"WWWWW  ERROR KBIN::VASP_Run: "+Message("time")+"  NPAR_REMOVE problems ");
                                    aus << "WWWWW  FIX NPAR_REMOVE - " << Message(aflags,"user,host,time") << endl;
                                    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                                    KBIN::XVASP_Afix_GENERIC("NPAR_REMOVE",xvasp,kflags,vflags,0.0,kflags.KBIN_MPI_NCPUS);  // this will be asked when NPAR is needed
                                    xfixed.flag("NPAR_REMOVE",TRUE);xfixed.flag("ALL",TRUE);
                                }
                            }
                            // ********* CHECK CSLOSHING PROBLEMS ******************
                            if(LDEBUG) cerr << "KBIN::VASP_Run: " << Message("time") << "  [CHECK CSLOSHING PROBLEMS]" << endl;
                            if(!vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("CSLOSHING") && !xfixed.flag("ALL")) { // check CSLOSHING
                                if(xwarning.flag("CSLOSHING") && !xfixed.flag("CSLOSHING")) { // Apply only ONCE
                                    KBIN::VASP_Error(xvasp,"WWWWW  ERROR KBIN::VASP_Run: "+Message("time")+"  CSLOSHING problems ");
                                    aus << "WWWWW  FIX CSLOSHING - " << Message(aflags,"user,host,time") << endl;
                                    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                                    KBIN::XVASP_Afix_GENERIC("CSLOSHING",xvasp,kflags,vflags);
                                    //	  KBIN_XVASP_Afix_CSLOSHING(xvasp,kflags,vflags);
                                    xfixed.flag("CSLOSHING",TRUE);xfixed.flag("ALL",TRUE);
                                }
                            }
                            // ********* CHECK DENTET PROBLEMS ******************
                            if(LDEBUG) cerr << "KBIN::VASP_Run: " << Message("time") << "  [CHECK DENTET PROBLEMS]" << endl;
                            if(!vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("DENTET") && !xfixed.flag("ALL")) { // CHECK FOR DENTET
                                if(xwarning.flag("DENTET") && !xfixed.flag("DENTET")) {  // Apply only ONCE
                                    KBIN::VASP_Error(xvasp,"WWWWW  ERROR KBIN::VASP_Run: "+Message("time")+"  DENTET problems ");
                                    aus << "WWWWW  FIX DENTET - " << Message(aflags,"user,host,time") << endl;
                                    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                                    KBIN::XVASP_Afix_GENERIC("DENTET",xvasp,kflags,vflags);
                                    xfixed.flag("DENTET",TRUE);xfixed.flag("ALL",TRUE);
                                }
                            }
                            // ********* VASP TO BE RESTARTED *********
                            if(LDEBUG) cerr << "KBIN::VASP_Run: " << Message("time") << "  [DONE WIHT CHECKS]" << endl;
                            if(xfixed.flag("ALL")) vasp_start=TRUE;
                            if(vasp_start) {
                                if(LDEBUG) cerr << "KBIN::VASP_Run: " << Message("time") << "  [VASP TO BE RESTARTED]" << endl;
                                aus << "00000  RESTART VASP   - " << Message(aflags,"user,host,time") << endl;
                                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                            }
                        }
                    }

                    if(LDEBUG) aus << "MMMMM  MESSAGE tested all the errors" << endl;
                    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);

                    bool Krun=TRUE;
                    if(!aurostd::FileExist(xvasp.Directory+"/vasp.out")) {Krun=FALSE;KBIN::VASP_Error(xvasp,"EEEEE  file does not exist=vasp.out ");}
                    if(aurostd::FileEmpty(xvasp.Directory+"/vasp.out")) {Krun=FALSE;KBIN::VASP_Error(xvasp,"EEEEE  file empty=vasp.out ");}
                    if(!aurostd::FileExist(xvasp.Directory+"/OUTCAR")) {Krun=FALSE;KBIN::VASP_Error(xvasp,"EEEEE  file does not exist=OUTCAR ");}
                    if(aurostd::FileEmpty(xvasp.Directory+"/OUTCAR")) {Krun=FALSE;KBIN::VASP_Error(xvasp,"EEEEE  file empty=OUTCAR ");}
                    if(!aurostd::FileExist(xvasp.Directory+"/CONTCAR")) {Krun=FALSE;KBIN::VASP_Error(xvasp,"EEEEE  file does not exist=CONTCAR ");}
                    if(aurostd::FileEmpty(xvasp.Directory+"/CONTCAR")) {Krun=FALSE;KBIN::VASP_Error(xvasp,"EEEEE  file xsempty=CONTCAR ");}

                    if(LDEBUG) cerr << "KBIN::VASP_Run: " << Message("time") << "  Krun=" << Krun << endl;
                    if(LDEBUG) {
                        aus << "KBIN::VASP_Run: " << Message("time") << "  Krun=" << Krun << endl;
                        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                    }

                    if(LDEBUG) cerr << "KBIN::VASP_Run: " << Message("time") << "  END" << endl;

                    if(LDEBUG) {
                        aus << "KBIN::VASP_Run: " << Message("time") << "  END" << endl;
                        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                    }

                    return Krun;
                    // ********* FINISH
                }
            } // namespace KBIN

            namespace KBIN {
                bool VASP_Run(_xvasp &xvasp,_aflags &aflags,_kflags &kflags,_vflags &vflags,string relax,bool qmwrite,ofstream &FileMESSAGE) {        // AFLOW_FUNCTION_IMPLEMENTATION
                    bool Krun=TRUE;
                    Krun=KBIN::VASP_Run(xvasp,aflags,kflags,vflags,FileMESSAGE);
                    if(!Krun) {KBIN::VASP_Error(xvasp,"EEEEE  Error in  \"KBIN::VASP_Run(_xvasp &xvasp,_aflags &aflags,_kflags &kflags,_vflags &vflags,string relax,bool qmwrite,ofstream &FileMESSAGE)\"");}
                    if(_VASP_CONTCAR_SAVE_) KBIN::VASP_CONTCAR_Save(xvasp,string(relax));
                    KBIN::VASP_Backup(xvasp,qmwrite,relax);
                    return Krun;
                }
            } // namespace KBIN

            namespace KBIN {
                bool VASP_Run(_xvasp &xvasp,_aflags &aflags,_kflags &kflags,_vflags &vflags,string relaxA,string relaxB,bool qmwrite,ofstream &FileMESSAGE) {        // AFLOW_FUNCTION_IMPLEMENTATION
                    bool Krun=TRUE;
                    if(relaxA!=relaxB)  {cerr << "ERROR KBIN::VASP_Run: " << Message("time") << " " << relaxA << "!=" << relaxB << endl;exit(0);}
                    Krun=KBIN::VASP_Run(xvasp,aflags,kflags,vflags,FileMESSAGE);
                    if(!Krun) {KBIN::VASP_Error(xvasp,"EEEEE  Error in  \"KBIN::VASP_Run(_xvasp &xvasp,_aflags &aflags,_kflags &kflags,_vflags &vflags,string relaxA,string relaxB,bool qmwrite,ofstream &FileMESSAGE)\"");}
                    // if(!Krun) return Krun;
                    if(_VASP_CONTCAR_SAVE_) KBIN::VASP_CONTCAR_Save(xvasp,string(relaxA));
                    KBIN::VASP_Backup(xvasp,qmwrite,relaxA);
                    KBIN::VASP_Recycle(xvasp,relaxB);
                    return Krun;
                }
            } // namespace KBIN

            namespace KBIN {
                bool VASP_RunFinished(_xvasp &xvasp,_aflags &aflags,ofstream &FileMESSAGE,bool verbose) {
                    ostringstream aus_exec,aus;
                    aurostd::StringstreamClean(aus_exec);
                    aurostd::StringstreamClean(aus);
                    // if(verbose) aus << "00000  MESSAGE RUN CHECK FINISHED : " << Message(aflags,"user,host,time") << endl;
                    // if(verbose) aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                    // OUTCAR DOES NOT EXIST
                    if(!aurostd::FileExist(xvasp.Directory+"/OUTCAR")) {
                        if(verbose) aus << "00000  MESSAGE RUN NOT FINISHED (OUTCAR does not exist) : " << Message(aflags,"user,host,time") << endl;
                        if(verbose) aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                        return FALSE;
                    }
                    // OTUCAR ESISTS BUT EMPTY
                    if(aurostd::FileEmpty(xvasp.Directory+"/OUTCAR")) {
                        if(verbose) aus << "00000  MESSAGE RUN NOT FINISHED (OUTCAR is empty) : " << Message(aflags,"user,host,time") << endl;
                        if(verbose) aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                        return FALSE;
                    }
                    // OTUCAR EXISTS
                    aus_exec << "cd " << xvasp.Directory << endl;
                    aus_exec << "cat OUTCAR | grep CPU > aflow.check_outcar.tmp " << endl;
                    aurostd::execute(aus_exec);
                    aurostd::StringstreamClean(aus_exec);
                    aus_exec << "cd " << xvasp.Directory << endl;
                    aus_exec << "rm -f aflow.check_outcar.tmp " << endl;
                    if(aurostd::substring_present_file(xvasp.Directory+"/aflow.check_outcar.tmp",aurostd::RemoveWhiteSpaces("Total CPU time used (sec)"),TRUE)) {
                        if(verbose) aus << "00000  MESSAGE RUN FINISHED (OUTCAR is complete) : " << Message(aflags,"user,host,time") << endl;
                        if(verbose) aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                        aurostd::execute(aus_exec);
                        return TRUE;
                    }
                    aurostd::execute(aus_exec);
                    if(verbose) aus << "00000  MESSAGE RUN NOT FINISHED (OUTCAR is incomplete) : " << Message(aflags,"user,host,time") << endl;
                    if(verbose) aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);

                    return FALSE;
                }
            } // namespace KBIN

            namespace KBIN {
                void WaitFinished(_xvasp &xvasp,_aflags &aflags,ofstream &FileMESSAGE,bool verbose) {
                    while(!KBIN::VASP_RunFinished(xvasp,aflags,FileMESSAGE,verbose)) {
                        aurostd::Sleep(_KVASP_CHECK_SLEEP_);
                    }
                }
            } // namespace KBIN

            namespace KBIN {
                void VASP_Error(_xvasp &xvasp,string message1,string message2,string message3) {
                    //  pthread_mutex_lock( &mutex_KVASP );
                    ofstream FileXVASP;
                    string FileNameXVASP;
                    FileNameXVASP=xvasp.Directory+"/"+DEFAULT_AFLOW_ERVASP_OUT;
                    if(1) {
                        FileXVASP.open(FileNameXVASP.c_str(),std::ios::out|std::ios::app);
                        // FileXVASP.open(FileNameXVASP.c_str(),std::ios::app);
                        FileXVASP << message1 << message2 << message3 << "  - " <<endl;
                        //   cerr << message1 << message2 << message3 << "  - " <<endl;
                        FileXVASP.flush(); FileXVASP.clear();
                        FileXVASP.close();
                    }
                    if(0) {
                        ostringstream aus;
                        aus << "echo \"" << message1 << message2 << message3 << " \" >> " << FileNameXVASP << endl;//endl;
                        aurostd::execute(aus);
                        //    cerr << aus.str();
                    }
                    //  pthread_mutex_unlock( &mutex_KVASP);
                }
            } // namespace KBIN

            namespace KBIN {
                void VASP_Error(_xvasp &xvasp,ofstream &FileMESSAGE,string message1,string message2,string message3) {
                    KBIN::VASP_Error(xvasp,message1,message2,message3);
                    ostringstream aus; aus << message1 << message2 << message3;
                    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                }
            } // namespace KBIN

            namespace KBIN {
                string VASP_Analyze(_xvasp &xvasp,bool qmwrite) {       // AFLOW_FUNCTION_IMPLEMENTATION
                    // CHECK ERRORS
                    bool error=FALSE;
                    // if(aurostd::FileEmpty(xvasp.Directory+"/EIGENVAL")) {KBIN::VASP_Error(xvasp,"EEEEE  ERROR KBIN::VASP_Analyze: Empty EIGENVAL ");error=TRUE;}
                    // if(aurostd::FileEmpty(xvasp.Directory+"/CHG"))      {KBIN::VASP_Error(xvasp,"EEEEE  ERROR KBIN::VASP_Analyze: Empty CHG ");error=TRUE;}
                    // if(aurostd::FileEmpty(xvasp.Directory+"/CHGCAR"))   {KBIN::VASP_Error(xvasp,"EEEEE  ERROR KBIN::VASP_Analyze: Empty CHGCAR ");error=TRUE;}
                    // if(aurostd::FileEmpty(xvasp.Directory+"/DOSCAR"))   {KBIN::VASP_Error(xvasp,"EEEEE  ERROR KBIN::VASP_Analyze: Empty DOSCAR ");error=TRUE;}
                    if(aurostd::FileEmpty(xvasp.Directory+"/CONTCAR"))  {KBIN::VASP_Error(xvasp,"EEEEE  ERROR KBIN::VASP_Analyze: Empty CONTCAR ");error=TRUE;}
                    if(aurostd::FileEmpty(xvasp.Directory+"/OUTCAR"))   {KBIN::VASP_Error(xvasp,"EEEEE  ERROR KBIN::VASP_Analyze: Empty OUTCAR  ");error=TRUE;}
                    if(aurostd::FileEmpty(xvasp.Directory+"/INCAR"))    {KBIN::VASP_Error(xvasp,"EEEEE  ERROR KBIN::VASP_Analyze: Empty INCAR   ");error=TRUE;}
                    if(aurostd::FileEmpty(xvasp.Directory+"/vasprun.xml"))    {KBIN::VASP_Error(xvasp,"EEEEE  ERROR KBIN::VASP_Analyze: Empty vasprun.xml   ");error=TRUE;}
                    if(error) return "";

                    // cerr << "# KBIN::VASP_Analyze BEGIN" << endl;
                    xvasp.str.qm_clear();
                    xvasp.str.qm_load(xvasp.Directory);
                    // LOAD OUTPUTS
                    stringstream strstream;
                    strstream.clear();
                    strstream.str(std::string());
                    strstream.setf(std::ios::fixed,std::ios::floatfield);
                    // OUTCAR OPERATIONS ---------------------------------------------------------------
                    strstream << "[AFLOW] **************************************************************************************************************************" << endl;
                    strstream << "[KBIN_ANALYZE]START_" << xvasp.AnalyzeLabel << endl;
                    if(xvasp.AnalyzeLabel!="dielectric_static" && xvasp.AnalyzeLabel!="dielectric_dynamic") {
                        strstream << "[AFLOW] **************************************************************************************************************************" << endl;
                        strstream << "# POSITION                                       TOTAL-FORCE (eV/Angst)               " << endl;
                        strstream << "[AFLOW] **************************************************************************************************************************" << endl;
                        strstream.precision(12);
                        for(uint i=0;i<xvasp.str.atoms.size();i++) {   // clear        (from the previous step)
                            for(uint j=1;j<=3;j++) {if(abs(xvasp.str.qm_positions.at(i)[j])<10.0) strstream << " ";if(xvasp.str.qm_positions.at(i)[j]>=0.0) strstream << " "; strstream << "   " << xvasp.str.qm_positions.at(i)[j] << " ";}
                            for(uint j=1;j<=3;j++) {if(abs(xvasp.str.qm_forces.at(i)[j])<10.0) strstream << " ";if(xvasp.str.qm_forces.at(i)[j]>=0.0) strstream << " "; strstream << "   " << xvasp.str.qm_forces.at(i)[j] << " ";}
                            strstream << endl;
                        }
                        strstream << "[AFLOW] **************************************************************************************************************************" << endl;
                        // OUTCAR OPERATIONS ---------------------------------------------------------------
                        strstream.setf(std::ios::scientific,std::ios::floatfield);
                        strstream.setf(std::ios::left,std::ios::adjustfield);
                        strstream << "E_cell=" << xvasp.str.qm_E_cell << "  (eV/cell)" << endl; 
                        strstream << "E_atom=" << xvasp.str.qm_E_atom << "  (eV/at)" << endl;
                        strstream << "H_cell=" << xvasp.str.qm_H_cell << "  (eV/cell)" << endl; 
                        strstream << "H_atom=" << xvasp.str.qm_H_atom << "  (eV/at)" << endl;
                        strstream << "PV_cell=" << xvasp.str.qm_PV_cell << "  (eV/cell)" << endl; 
                        strstream << "PV_atom=" << xvasp.str.qm_PV_atom << "  (eV/at)" << endl;
                        strstream << "mag_cell=" << xvasp.str.qm_mag_cell << "  (mu/cell)" << endl;
                        strstream << "mag_atom="<< xvasp.str.qm_mag_atom << "  (mu/at)" << endl;
                        xstructure qm_str(xvasp.str);    // suck it in !
                        // qm_str=xvasp.str;
                        qm_str.qm_recycle();
                        strstream << "[AFLOW] **************************************************************************************************************************" << endl;
                        strstream << "[VASP_POSCAR_MODE_EXPLICIT]START " << endl;
                        strstream << qm_str;
                        strstream << "[VASP_POSCAR_MODE_EXPLICIT]STOP " << endl;
                    }

                    //  if(aurostd::substring2bool(aurostd::execute2string("grep LEPSILON "+xvasp.Directory+"/OUTCAR"),"LEPSILON=T",TRUE))
                    if(xvasp.AnalyzeLabel=="dielectric_static") {
                        strstream << "[AFLOW] **************************************************************************************************************************" << endl;
                        strstream << "[KBIN_ANALYZE]START_DIELECTRIC_STATIC" << endl;
                        vector<string> vlines,tokens;
                        aurostd::string2vectorstring(aurostd::execute2string("grep -A 4 \"MACROSCOPIC STATIC DIELECTRIC TENSOR\" "+xvasp.Directory+"/OUTCAR  | tail -n 3"),vlines);
                        xmatrix<double> epsilon(3,3);
                        if(vlines.size()==3) {
                            for(uint i=1;i<=3;i++) {
                                aurostd::string2tokens(vlines.at(i-1),tokens," ");
                                if(tokens.size()==3) {
                                    for(uint j=1;j<=3;j++)
                                        epsilon(i,j)=aurostd::string2utype<double>(tokens.at(j-1));
                                }
                            }
                        }
                        strstream << " epsilon = " << endl;
                        strstream << "  " << epsilon(1,1) << "  " << epsilon(1,2) << "  " << epsilon(1,3) << "  " << endl;
                        strstream << "  " << epsilon(2,1) << "  " << epsilon(2,2) << "  " << epsilon(2,3) << "  " << endl;
                        strstream << "  " << epsilon(3,1) << "  " << epsilon(3,2) << "  " << epsilon(3,3) << "  " << endl;
                        strstream << "[KBIN_ANALYZE]STOP_DIELECTRIC_STATIC" << endl;
                    }
                    //  if(aurostd::substring2bool(aurostd::execute2string("grep LOPTICS "+xvasp.Directory+"/OUTCAR"),"LOPTICS=T",TRUE)) {
                    if(xvasp.AnalyzeLabel=="dielectric_dynamic") {
                        strstream << "[AFLOW] **************************************************************************************************************************" << endl;
                        strstream << " DIELECTRIC DYNAMIC " << endl;  
                        vector<string> vlines,tokens;string line;
                        vector<xvector<double> > vepsilonIMAG,vepsilonREAL;
                        aurostd::file2vectorstring(xvasp.Directory+"/OUTCAR",vlines);
                        uint IMAG_start=0,IMAG_end=0,REAL_start=0,REAL_end=0;
                        for(uint i=0;i<vlines.size();i++) {
                            if(aurostd::substring2bool(vlines.at(i),"IMAGINARY DIELECTRIC FUNCTION")) IMAG_start=i;
                            if(aurostd::substring2bool(vlines.at(i),"REAL DIELECTRIC FUNCTION")) {REAL_start=i;IMAG_end=i-2;REAL_end=(IMAG_end-IMAG_start)+REAL_start;}
                        }
                        // WRITING DIELECTRIC_IMAGINARY
                        strstream << "[KBIN_ANALYZE]START_DIELECTRIC_IMAGINARY" << endl;
                        for(uint i=IMAG_start;i<=IMAG_end;i++) {
                            line=vlines.at(i);
                            aurostd::StringSubst(line,"-0.000000","0");aurostd::StringSubst(line,"0.000000","0");aurostd::StringSubst(line,"\t"," ");aurostd::StringSubst(line,"  "," ");aurostd::StringSubst(line,"  "," ");aurostd::StringSubst(line,"  "," ");
                            strstream << line << endl;
                            if(i>=IMAG_start+3) {
                                xvector<double> epsilonIMAG(7);
                                aurostd::string2tokens(vlines.at(i),tokens," ");
                                if(tokens.size()==7) for(uint j=0;j<tokens.size();j++) epsilonIMAG(j+1)=aurostd::string2utype<double>(tokens.at(j));
                                vepsilonIMAG.push_back(epsilonIMAG);
                            }
                        }
                        strstream << "[KBIN_ANALYZE]STOPT_DIELECTRIC_IMAGINARY" << endl;
                        // WRITING DIELECTRIC_REAL
                        strstream << "[KBIN_ANALYZE]START_DIELECTRIC_REAL" << endl;
                        for(uint i=REAL_start;i<=REAL_end;i++) {
                            line=vlines.at(i);
                            aurostd::StringSubst(line,"-0.000000","0");aurostd::StringSubst(line,"0.000000","0");aurostd::StringSubst(line,"\t"," ");aurostd::StringSubst(line,"  "," ");aurostd::StringSubst(line,"  "," ");aurostd::StringSubst(line,"  "," ");
                            strstream << line << endl;
                            if(i>=REAL_start+3) {
                                xvector<double> epsilonREAL(7);
                                aurostd::string2tokens(vlines.at(i),tokens," ");
                                if(tokens.size()==7) for(uint j=0;j<tokens.size();j++) epsilonREAL(j+1)=aurostd::string2utype<double>(tokens.at(j));
                                vepsilonREAL.push_back(epsilonREAL);
                            }
                        }
                        strstream << "[KBIN_ANALYZE]STOPT_DIELECTRIC_REAL" << endl;
                    }  
                    strstream << "[AFLOW] **************************************************************************************************************************" << endl;
                    strstream << "[KBIN_ANALYZE]STOP_" << xvasp.AnalyzeLabel << endl;
                    strstream << "[AFLOW] **************************************************************************************************************************" << endl;

                    if(qmwrite) {
                        string FileNameXVASP=xvasp.Directory+"/"+DEFAULT_AFLOW_QMVASP_OUT;
                        stringstream FileXVASPout;
                        if(aurostd::FileExist(xvasp.Directory+"/"+DEFAULT_AFLOW_QMVASP_OUT)) { //RECYCLE PREVIOUS STUFF
                            stringstream FileXVASPin;
                            aurostd::file2stringstream(xvasp.Directory+"/"+DEFAULT_AFLOW_QMVASP_OUT,FileXVASPin);
                            FileXVASPout << FileXVASPin.str();
                        }
                        FileXVASPout << strstream.str();
                        aurostd::stringstream2file(FileXVASPout,xvasp.Directory+"/"+DEFAULT_AFLOW_QMVASP_OUT);
                    }
                    xvasp.str.qm_calculated=TRUE;
                    //  cerr << "# KBIN::VASP_Analyze END" << endl;
                    return strstream.str();
                }
                }  // namespace KBIN

                namespace KBIN {
                    void GenerateAflowinFromVASPDirectory(_aflags &aflags) {        // AFLOW_FUNCTION_IMPLEMENTATION
                        ifstream FileSUBDIR;string FileNameSUBDIR;
                        FileNameSUBDIR=aflags.Directory;
                        FileSUBDIR.open(FileNameSUBDIR.c_str(),std::ios::in);
                        FileSUBDIR.clear();FileSUBDIR.close();
                        ostringstream aus;

                        if(aflags.Directory.at(0)!='/' && aflags.Directory.at(0)!='.' && aflags.Directory.at(0)!=' ') aflags.Directory="./"+aflags.Directory;

                        if(!FileSUBDIR) {                                                                                           // ******* Directory is non existent
                            aus << "EEEEE  DIRECTORY_NOT_FOUND = " << Message(aflags,"user,host,time") << endl;
                            aurostd::PrintMessageStream(aus,XHOST.QUIET);
                        } else {                                                                                                    // ******* Directory EXISTS
                            // Check LOCK again
                            // ifstream FileLOCK0;string FileNameLOCK0=aflags.Directory+"/"+_AFLOWLOCK_;    FileLOCK0.open(FileNameLOCK0.c_str(),std::ios::in);FileLOCK0.close();
                            // ifstream FileLOCK1;string FileNameLOCK1=aflags.Directory+"/"+_AFLOWLOCK_+".gz"; FileLOCK1.open(FileNameLOCK1.c_str(),std::ios::in);FileLOCK1.close();
                            // ifstream FileLOCK2;string FileNameLOCK2=aflags.Directory+"/"+_AFLOWLOCK_+".bz2";FileLOCK2.open(FileNameLOCK2.c_str(),std::ios::in);FileLOCK2.close();
                            // ifstream FileSKIP0;string FileNameSKIP0=aflags.Directory+"/SKIP";    FileSKIP0.open(FileNameSKIP0.c_str(),std::ios::in);FileSKIP0.close();
                            // ifstream FileSKIP1;string FileNameSKIP1=aflags.Directory+"/SKIP.gz"; FileSKIP1.open(FileNameSKIP1.c_str(),std::ios::in);FileSKIP1.close();
                            // ifstream FileSKIP2;string FileNameSKIP2=aflags.Directory+"/SKIP.bz2";FileSKIP2.open(FileNameSKIP2.c_str(),std::ios::in);FileSKIP2.close();
                            // // CHECK FOR LOCK
                            // if(FileLOCK0 || FileLOCK1 || FileLOCK2) {                                                                 // ******* Directory is locked
                            // 	// LOCK exist, then RUN already RUN
                            // 	aus << "EEEEE  LOCKED " << Message(aflags,"user,host,time") << endl;
                            // 	aurostd::PrintMessageStream(aus,XHOST.QUIET);
                            // 	exit(1);
                            // }
                            // if(FileSKIP0 || FileSKIP1 || FileSKIP2) {                                                                 // ******* Directory is skipped
                            // 	// SKIP exist, then RUN already RUN
                            // 	aus << "EEEEE  SKIPPED " << Message(aflags,"user,host,time") << endl;
                            // 	aurostd::PrintMessageStream(aus,XHOST.QUIET);
                            // 	exit(1);
                            // }
                            if(aurostd::FileExist(aflags.Directory+"/"+_AFLOWLOCK_) || aurostd::EFileExist(aflags.Directory+"/"+_AFLOWLOCK_))	{ // ******* Directory is locked
                                // LOCK exist, then RUN already RUN
                                aus << "EEEEE  LOCKED " << Message(aflags,"user,host,time") << endl;
                                aurostd::PrintMessageStream(aus,XHOST.QUIET);
                                exit(1);
                            }
                            if(aurostd::FileExist(aflags.Directory+"/SKIP") || aurostd::EFileExist(aflags.Directory+"/SKIP")) {	// ******* Directory is skipped
                                // SKIP exist, then RUN already RUN
                                aus << "EEEEE  SKIPPED " << Message(aflags,"user,host,time") << endl;
                                aurostd::PrintMessageStream(aus,XHOST.QUIET);
                                exit(1);
                            }

                            // ******* Directory is un locked/skipped
                            /// ******************************************************************
                            // RESET LOCK
                            ofstream FileLOCK;
                            string FileNameLOCK=aflags.Directory+"/"+_AFLOWLOCK_;
                            FileLOCK.open(FileNameLOCK.c_str(),std::ios::out);
                            /// ******************************************************************
                            // CHECK FOR INCAR KPOINTS POSCAR POTCAR
                            ifstream FileINCAR;string FileNameINCAR=aflags.Directory+"/INCAR";FileINCAR.open(FileNameINCAR.c_str(),std::ios::in);
                            if(!FileINCAR)  {                                                                                        // ******* INCAR does not exist
                                aus << "EEEEE  INCAR ABSENT  = " << Message(aflags,"user,host,time") << endl;
                                aurostd::PrintErrorStream(FileLOCK,aus,XHOST.QUIET);
                            }
                            ifstream FileKPOINTS;string FileNameKPOINTS=aflags.Directory+"/KPOINTS";FileKPOINTS.open(FileNameKPOINTS.c_str(),std::ios::in);
                            if(!FileKPOINTS)  {                                                                                        // ******* KPOINTS does not exist
                                aus << "EEEEE  KPOINTS ABSENT  = " << Message(aflags,"user,host,time") << endl;
                                aurostd::PrintErrorStream(FileLOCK,aus,XHOST.QUIET);
                            }
                            ifstream FilePOSCAR;string FileNamePOSCAR=aflags.Directory+"/POSCAR";FilePOSCAR.open(FileNamePOSCAR.c_str(),std::ios::in);
                            if(!FilePOSCAR)  {                                                                                        // ******* POSCAR does not exist
                                aus << "EEEEE  POSCAR ABSENT  = " << Message(aflags,"user,host,time") << endl;
                                aurostd::PrintErrorStream(FileLOCK,aus,XHOST.QUIET);
                            }
                            ifstream FilePOTCAR;string FileNamePOTCAR=aflags.Directory+"/POTCAR";FilePOTCAR.open(FileNamePOTCAR.c_str(),std::ios::in);
                            if(!FilePOTCAR)  {                                                                                        // ******* POTCAR does not exist
                                aus << "EEEEE  POTCAR ABSENT  = " << Message(aflags,"user,host,time") << endl;
                                aurostd::PrintErrorStream(FileLOCK,aus,XHOST.QUIET);
                            }
                            // ----------------------------------------------------------------------------------------------------
                            if(FileINCAR && FileKPOINTS && FilePOSCAR && FilePOTCAR) {
                                // VASP INCAR KPOINTS POSCAR POTCAR ARE PRESENT
                                /// ******************************************************************
                                // WRITE LOCK
                                aus << "MMMMM  AFLOW VERSION " << string(AFLOW_VERSION) << " Automatic-Flow - " << Message(aflags,"user,host,time") << endl;
                                aus << "MMMMM  (C) " << XHOST.Copyright_Years << ", Stefano Curtarolo - Duke University   - " << Message(aflags,"user,host,time") << endl;
                                aus << "MMMMM  High-Throughput ab-initio Computing - " << Message(aflags,"user,host,time") << endl;
                                aurostd::PrintMessageStream(FileLOCK,aus,XHOST.QUIET);
                                aus << "00000  MESSAGE GENERATING " << _AFLOWIN_ << " from VASP-xCARs files " << Message(aflags,"user,host,time") << endl;
                                aurostd::PrintMessageStream(FileLOCK,aus,XHOST.QUIET);
                                /// ******************************************************************
                                // RESET AFLOWIN
                                ofstream FileAFLOWIN;
                                string FileNameAFLOWIN=aflags.Directory+"/"+_AFLOWIN_;
                                FileAFLOWIN.open(FileNameAFLOWIN.c_str(),std::ios::out);
                                /// ******************************************************************
                                // WRITE AFLOWIN
                                // WRITE TITLE
                                string str1,str2;
                                getline(FileINCAR,str1);
                                FileINCAR.seekg(0);
                                FileAFLOWIN << AFLOWIN_SEPARATION_LINE << endl;
                                FileAFLOWIN << "[AFLOW] Automatically generated from XCARS by aflow/aflowd " << string(AFLOW_VERSION) << endl;
                                FileAFLOWIN << "[AFLOW] Automatic-Flow - " << Message(aflags,"user,host,time") << endl;
                                FileAFLOWIN << AFLOWIN_SEPARATION_LINE << endl;
                                // WRITE HEADER
                                FileAFLOWIN << AFLOWIN_SEPARATION_LINE << endl;
                                FileAFLOWIN << "[AFLOW] " << str1 << endl;
                                FileAFLOWIN << AFLOWIN_SEPARATION_LINE << endl;
                                FileAFLOWIN << "[AFLOW] input file for aflow " << endl;
                                FileAFLOWIN << "[AFLOW] comments with label " << endl;
                                FileAFLOWIN << "[AFLOW] separating with __ the options makes them ignored " << endl;
                                FileAFLOWIN << "[AFLOW_MODE=VASP] " << endl;
                                FileAFLOWIN << "[VASP] *************************************************** " << endl;
                                for(int i=0;i<(int) XHOST.argv.size()-1;i++) {
                                    str1=XHOST.argv.at(i);
                                    str2=XHOST.argv.at(i+1);
                                    if(str1=="--set" && str2.at(0)=='[') {
                                        aus << "00000  MESSAGE Adding " << str2 << " to " << _AFLOWIN_ << " " << Message(aflags,"user,host,time") << endl;
                                        aurostd::PrintMessageStream(FileLOCK,aus,XHOST.QUIET);
                                        FileAFLOWIN << str2 << endl;
                                    }
                                }
                                // WRITE INCAR
                                FileAFLOWIN << AFLOWIN_SEPARATION_LINE << endl;
                                FileAFLOWIN << "[VASP_INCAR_MODE_EXPLICIT]" << endl;
                                while (getline(FileINCAR,str1)) FileAFLOWIN << "[VASP_INCAR_FILE]" << str1 << endl;
                                // WRITE KPOINTS
                                FileAFLOWIN << AFLOWIN_SEPARATION_LINE << endl;
                                FileAFLOWIN << "[VASP_KPOINTS_MODE_EXPLICIT]" << endl;
                                while (getline(FileKPOINTS,str1)) FileAFLOWIN << "[VASP_KPOINTS_FILE]" << str1 << endl;
                                // WRITE POSCAR
                                FileAFLOWIN << AFLOWIN_SEPARATION_LINE << endl;
                                FileAFLOWIN << "[VASP_POSCAR_MODE_EXPLICIT]" << endl;
                                while (getline(FilePOSCAR,str1)) FileAFLOWIN << "[VASP_POSCAR_FILE]" << str1 << endl;
                                // WRITE POTCAR
                                FileAFLOWIN << AFLOWIN_SEPARATION_LINE << endl;
                                FileAFLOWIN << "[VASP_POTCAR_MODE_EXPLICIT]" << endl;
                                while (getline(FilePOTCAR,str1)) FileAFLOWIN << str1 << endl;
                                // close everything.
                                FileINCAR.clear();FileINCAR.close();
                                FileKPOINTS.clear();FileKPOINTS.close();
                                FilePOSCAR.clear();FilePOSCAR.close();
                                FilePOTCAR.clear();FilePOTCAR.close();
                                FileAFLOWIN.flush();FileAFLOWIN.clear();FileAFLOWIN.close();
                                /// ******************************************************************
                                // everything is done. check if we nned to delete VASP FILES
                                if(aurostd::args2flag(XHOST.argv,"--delete_xcars")) {
                                    aus << "00000  MESSAGE Removing vasp files in " << Message(aflags,"user,host,time") << endl;
                                    aurostd::PrintMessageStream(FileLOCK,aus,XHOST.QUIET);
                                    aus << "cd " << aflags.Directory << endl;
                                    aus << "rm -f `ls | grep -v " << _AFLOWIN_ << " | grep -v LOCK ` " << endl;
                                    aurostd::execute(aus);
                                }
                            }
                            FileLOCK.flush();FileLOCK.clear();FileLOCK.close();
                        }
                    }
                } // namespace KBIN

                namespace KBIN {
                    void VASP_Backup(_xvasp& xvasp,bool qmwrite,string relax) {        // AFLOW_FUNCTION_IMPLEMENTATION
                        xvasp.AnalyzeLabel=relax;
                        KBIN::VASP_Analyze(xvasp,qmwrite);
                        ostringstream aus;
                        deque<string> vext; aurostd::string2tokens(".bz2,.xz,.gz",vext,",");

                        if(aurostd::FileExist(xvasp.Directory+"/core")) aurostd::execute("rm -f "+xvasp.Directory+"/core");
                        for(uint iext=0;iext<vext.size();iext++) { 
                            if(aurostd::FileExist(xvasp.Directory+"/core"+vext.at(iext))) aurostd::execute("rm -f "+xvasp.Directory+"/core"+vext.at(iext));
                        }
                        if(!xvasp.aopts.flag("FLAG::WAVECAR_PRESERVED") && aurostd::FileExist(xvasp.Directory+"/WAVECAR")) aurostd::execute("rm -f  WAVECAR");
                        if(!xvasp.aopts.flag("FLAG::WAVEDER_PRESERVED") && aurostd::FileExist(xvasp.Directory+"/WAVEDER")) aurostd::execute("rm -f  WAVEDER");
                        if(aurostd::FileExist(xvasp.Directory+"/aflow.qsub.run")) aurostd::execute("rm -f  aflow.qsub.run");
                        if(aurostd::FileExist(xvasp.Directory+"/aflow.qsub.out")) aurostd::execute("rm -f  aflow.qsub.out");
                        if(aurostd::FileExist(xvasp.Directory+"/AECCAR0")) aurostd::execute("mv "+xvasp.Directory+"/AECCAR0 "+xvasp.Directory+"/AECCAR0"+"."+relax);  // BADER
                        if(aurostd::FileExist(xvasp.Directory+"/AECCAR1")) aurostd::execute("mv "+xvasp.Directory+"/AECCAR1 "+xvasp.Directory+"/AECCAR1"+"."+relax);  // BADER
                        if(aurostd::FileExist(xvasp.Directory+"/AECCAR2")) aurostd::execute("mv "+xvasp.Directory+"/AECCAR2 "+xvasp.Directory+"/AECCAR2"+"."+relax);  // BADER
                        if(aurostd::FileExist(xvasp.Directory+"/CHG")) aurostd::execute("mv "+xvasp.Directory+"/CHG "+xvasp.Directory+"/CHG"+"."+relax);
                        if(aurostd::FileExist(xvasp.Directory+"/CHGCAR")) aurostd::execute("mv "+xvasp.Directory+"/CHGCAR "+xvasp.Directory+"/CHGCAR"+"."+relax);
                        if(aurostd::FileExist(xvasp.Directory+"/CONTCAR")) aurostd::execute("mv "+xvasp.Directory+"/CONTCAR "+xvasp.Directory+"/CONTCAR"+"."+relax);
                        if(aurostd::FileExist(xvasp.Directory+"/DYNMAT")) aurostd::execute("mv "+xvasp.Directory+"/DYNMAT "+xvasp.Directory+"/DYNMAT"+"."+relax);
                        if(aurostd::FileExist(xvasp.Directory+"/DOSCAR")) aurostd::execute("mv "+xvasp.Directory+"/DOSCAR "+xvasp.Directory+"/DOSCAR"+"."+relax);
                        if(aurostd::FileExist(xvasp.Directory+"/ELFCAR")) aurostd::execute("mv "+xvasp.Directory+"/ELFCAR "+xvasp.Directory+"/ELFCAR"+"."+relax);  // ELF
                        if(aurostd::FileExist(xvasp.Directory+"/EIGENVAL")) aurostd::execute("mv "+xvasp.Directory+"/EIGENVAL "+xvasp.Directory+"/EIGENVAL"+"."+relax);
                        if(aurostd::FileExist(xvasp.Directory+"/IBZKPT")) aurostd::execute("mv "+xvasp.Directory+"/IBZKPT "+xvasp.Directory+"/IBZKPT"+"."+relax);
                        if(aurostd::FileExist(xvasp.Directory+"/INCAR")) aurostd::execute("mv "+xvasp.Directory+"/INCAR "+xvasp.Directory+"/INCAR"+"."+relax);
                        if(aurostd::FileExist(xvasp.Directory+"/KPOINTS")) aurostd::execute("mv "+xvasp.Directory+"/KPOINTS "+xvasp.Directory+"/KPOINTS"+"."+relax);
                        if(aurostd::FileExist(xvasp.Directory+"/OSZICAR")) aurostd::execute("mv "+xvasp.Directory+"/OSZICAR "+xvasp.Directory+"/OSZICAR"+"."+relax);
                        if(aurostd::FileExist(xvasp.Directory+"/OUTCAR")) aurostd::execute("mv "+xvasp.Directory+"/OUTCAR "+xvasp.Directory+"/OUTCAR"+"."+relax);
                        if(aurostd::FileExist(xvasp.Directory+"/PCDAT")) aurostd::execute("mv "+xvasp.Directory+"/PCDAT "+xvasp.Directory+"/PCDAT"+"."+relax);
                        if(aurostd::FileExist(xvasp.Directory+"/POSCAR")) aurostd::execute("mv "+xvasp.Directory+"/POSCAR "+xvasp.Directory+"/POSCAR"+"."+relax);
                        if(aurostd::FileExist(xvasp.Directory+"/POTCAR")) aurostd::execute("mv "+xvasp.Directory+"/POTCAR "+xvasp.Directory+"/POTCAR"+"."+relax);
                        if(aurostd::FileExist(xvasp.Directory+"/PROCAR")) aurostd::execute("mv "+xvasp.Directory+"/PROCAR "+xvasp.Directory+"/PROCAR"+"."+relax);
                        if(aurostd::FileExist(xvasp.Directory+"/XDATCAR")) aurostd::execute("mv "+xvasp.Directory+"/XDATCAR "+xvasp.Directory+"/XDATCAR"+"."+relax);
                        if(xvasp.aopts.flag("FLAG::WAVECAR_PRESERVED") && aurostd::FileExist(xvasp.Directory+"/WAVECAR")) aurostd::execute("mv "+xvasp.Directory+"/WAVECAR "+xvasp.Directory+"/WAVECAR"+"."+relax);
                        if(xvasp.aopts.flag("FLAG::WAVEDER_PRESERVED") && aurostd::FileExist(xvasp.Directory+"/WAVEDER")) aurostd::execute("mv "+xvasp.Directory+"/WAVEDER "+xvasp.Directory+"/WAVEDER"+"."+relax);
                        if(aurostd::FileExist(xvasp.Directory+"/vasp.out")) aurostd::execute("mv "+xvasp.Directory+"/vasp.out "+xvasp.Directory+"/vasp.out"+"."+relax);
                        if(aurostd::FileExist(xvasp.Directory+"/vasprun.xml")) aurostd::execute("mv "+xvasp.Directory+"/vasprun.xml "+xvasp.Directory+"/vasprun.xml"+"."+relax);
                    }
                } // namespace KBIN

                namespace KBIN {
                    void VASP_CONTCAR_Save(_xvasp xvasp,string relax) {        // AFLOW_FUNCTION_IMPLEMENTATION
                        if(aurostd::FileExist(xvasp.Directory+string("/CONTCAR")))
                            if(!aurostd::FileEmpty(xvasp.Directory+string("/CONTCAR"))) {
                                ostringstream aus;
                                aus << "cd " << xvasp.Directory << endl;
                                aus << "echo \"[AFLOW] SELF-MODIFICATION \" >> " << _AFLOWIN_ << " " << endl;
                                aus << "echo \"[AFLOW] Recycling CONTCAR of " << relax << " \" >> " << _AFLOWIN_ << " " << endl;

                                /*//KESONG fixes potential aflow command problem in future, 2018-11-25
                                  string str_tmp = aurostd::file2string(xvasp.Directory+string("/CONTCAR"));   
                                  stringstream ss_tmp; ss_tmp << ""; 
                                  ss_tmp << str_tmp << endl;
                                  string str_tmp2 = pflow::AFLOWIN(ss_tmp) + " >> ";
                                  aurostd::string2file(str_tmp2, "test.txt");
                                  aus <<   str_tmp2  <<  _AFLOWIN_ << " " << endl;
                                  */

                                aus << "cat CONTCAR | aflow --aflowin  >> " << _AFLOWIN_ << " " << endl;
                                aurostd::execute(aus);
                                aus << "cd " << xvasp.Directory << endl;
                                aus << "cat " << _AFLOWIN_ << " | sed \"s/\\[VASP_FORCE_OPTION\\]VOLUME/#\\[VASP_FORCE_OPTION\\]VOLUME/g\" | sed \"s/##\\[/#\\[/g\" > aflow.tmp && mv aflow.tmp " << _AFLOWIN_ << "" << endl; // PRESERVE VOLUME
                                aurostd::execute(aus);
                            }
                    }
                } // namespace KBIN

                namespace KBIN {
                    void VASP_Recycle(_xvasp xvasp,string relax) {        // AFLOW_FUNCTION_IMPLEMENTATION
                        ostringstream aus;
                        aus << "cd " << xvasp.Directory << endl;
                        aus << "cp CONTCAR." << relax << " POSCAR" << endl;
                        aus << "cp INCAR." << relax << " INCAR" << endl;
                        aus << "cp KPOINTS." << relax << " KPOINTS" << endl;
                        aus << "cp POTCAR." << relax << " POTCAR" << endl;
                        aurostd::execute(aus);
                    }
                } // namespace KBIN

                namespace KBIN {
                    void VASP_Recycle(_xvasp xvasp,int relax_number) {        // AFLOW_FUNCTION_IMPLEMENTATION
                        ostringstream aus;

                        deque<string> vext; aurostd::string2tokens(".bz2,.xz,.gz",vext,",");
                        deque<string> vzip; aurostd::string2tokens("bzip2,xz,gzip",vzip,",");

                        aus << "cd " << xvasp.Directory << endl;
                        for(uint iext=0;iext<vext.size();iext++) { 
                            aus << vzip.at(iext) << " -dqf *" << vext.at(iext) << endl;
                        }
                        aus << "cp CONTCAR.relax" << relax_number << " POSCAR" << endl;
                        aus << "cp INCAR.relax" << relax_number << " INCAR" << endl;
                        aus << "cp KPOINTS.relax" << relax_number << " KPOINTS" << endl;
                        aus << "cp POTCAR.relax" << relax_number << " POTCAR" << endl;
                        aurostd::execute(aus);
                    }
                } // namespace KBIN

                namespace KBIN {
                    void VASP_RecycleExtraFile(_xvasp xvasp,string xfile,string relax) {        // AFLOW_FUNCTION_IMPLEMENTATION
                        ostringstream aus;
                        aus << "cd " << xvasp.Directory << endl;
                        aus << "cp " << xfile << "." << relax << " " << xfile << endl;
                        aurostd::execute(aus);
                    }
                } // namespace KBIN

                namespace KBIN {
                    void VASP_RecycleExtraFile(_xvasp xvasp,string xfile,int relax_number) {        // AFLOW_FUNCTION_IMPLEMENTATION
                        ostringstream aus;

                        deque<string> vext; aurostd::string2tokens(".bz2,.xz,.gz",vext,",");
                        deque<string> vzip; aurostd::string2tokens("bzip2,xz,gzip",vzip,",");

                        aus << "cd " << xvasp.Directory << endl;
                        for(uint iext=0;iext<vext.size();iext++) { 
                            aus << vzip.at(iext) << " -dqf " <<  xfile << vext.at(iext) << endl;
                        }
                        aus << "cp " <<  xfile << ".relax" << relax_number << " " << xfile << endl;
                        aurostd::execute(aus);
                    }
                } // namespace KBIN

                namespace KBIN {
                    void VASP_BackupOriginal(_xvasp xvasp) {        // AFLOW_FUNCTION_IMPLEMENTATION
                        ostringstream aus;
                        aus << "cd " << xvasp.Directory << endl;
                        aus << "cp KPOINTS KPOINTS.orig" << endl;
                        aus << "cp INCAR   INCAR.orig" << endl;
                        aus << "cp POSCAR  POSCAR.orig" << endl;
                        aurostd::execute(aus);
                    }
                } // namespace KBIN

                namespace KBIN {
                    bool VASP_CheckUnconvergedOSZICAR(string dir) {
                        uint ielectrons=0,issues=0,cutoff=3;
                        vector<string> vlines,vrelax,tokens;
                        aurostd::file2vectorstring(dir+"/OSZICAR",vlines);
                        for(uint i=0;i<vlines.size();i++)
                            if(aurostd::substring2bool(vlines.at(i),"F="))
                                vrelax.push_back(vlines.at(i-1));
                        if(vrelax.size()<cutoff) return FALSE; // no problem
                        // otherwise check for issues.
                        for(uint i=0;i<vrelax.size()&&i<cutoff;i++) {
                            aurostd::string2tokens(vrelax.at(i),tokens," ");
                            ielectrons=aurostd::string2utype<uint>(tokens.at(1));
                            if(ielectrons==60 || ielectrons==120) issues++;
                        }
                        if(issues==cutoff) return TRUE;
                        return FALSE;
                    }
                } // namespace KBIN

                // ***************************************************************************
                // functions written by CAMILO CALDERON
                // 2013: camilo.calderon@duke.edu

                // todo:
                // Finish the DYNADIEL tag
                // OUTCAR file & type as a separate subroutine
                // Add more options to the statdiel tag (various dielectric tensor types)

                // ***************************************************************************
                namespace KBIN {
                    void GetStatDiel(string& outcar, xvector<double>& eigr, xvector<double>& eigi) { // loop GetStatDiel
                        int PATH_LENGTH_MAX = 1024 ;
                        char work_dir[PATH_LENGTH_MAX] ;
                        string outcarfile, outcarpath ;
                        string outcarpath_tmp = aurostd::TmpFileCreate("OUTCARc1.tmp") ;
                        vector<string> outcarlines, endline, startline, vasptoken ;
                        xmatrix<double> statdiel(3,3), eigenvec(3,3) ;
                        double eps = 1.0E-5 ; // need to define this more rigorously
                        getcwd(work_dir, PATH_LENGTH_MAX) ;

                        if(!aurostd::FileExist(outcar)) {
                            cout << "KBIN::GetStatDiel: check filename || file missing" << endl ;
                            exit(0) ;
                        } else {
                            outcarpath = "/" + outcar ;
                            outcarpath = work_dir + outcarpath ;
                            vector<string> outcardata ;
                            aurostd::string2tokens(outcarpath, outcardata, ".") ;
                            if(outcardata.at(outcardata.size()-1) == "bz2") { // compressed option
                                aurostd::execute("bzcat " + outcarpath + " > " + outcarpath_tmp) ;
                                aurostd::file2vectorstring(outcarpath_tmp, outcarlines) ;
                            } else if(outcardata.at(outcardata.size()-1) == "xz") { // compressed option
                                aurostd::execute("xzcat " + outcarpath + " > " + outcarpath_tmp) ;
                                aurostd::file2vectorstring(outcarpath_tmp, outcarlines) ;
                            } else if(outcardata.at(outcardata.size()-1) == "gz") { // compressed option
                                aurostd::execute("gzcat " + outcarpath + " > " + outcarpath_tmp) ;
                                aurostd::file2vectorstring(outcarpath_tmp, outcarlines) ;
                            } else { // plain text option
                                aurostd::execute("cat " + outcarpath + " > " + outcarpath_tmp) ;
                                aurostd::file2vectorstring(outcarpath_tmp, outcarlines) ;
                            }
                        }
                        // check the loaded OUTCAR
                        aurostd::string2tokens(outcarlines.at(0),startline," ");
                        aurostd::string2tokens(startline.at(0),vasptoken,".");
                        aurostd::string2tokens(outcarlines.at(outcarlines.size()-1),endline," ");
                        if(vasptoken.at(0) != "vasp" || endline.at(0) != "Voluntary") { // first and last line check
                            cerr << "OUTCAR file is probably corrupt" << endl ;
                            exit(0) ;
                        }
                        uint sec_count = 0 ;
                        for (uint ii=outcarlines.size()-12 ; ii<outcarlines.size() ; ii++) { // presence timing information check
                            vector<string> timetoken ;
                            aurostd::string2tokens(outcarlines.at(ii),timetoken," ") ;
                            if(timetoken.size() > 0) {
                                for (uint jj=0 ; jj<timetoken.size() ; jj++)
                                { if(timetoken.at(jj) == "(sec):") sec_count+=1 ; }
                            }
                        }
                        if(sec_count != 4) { // first and last line check
                            cerr << "OUTCAR file is probably corrupt" << endl ;
                            exit(0) ;
                        }
                        // OUTCAR is now in memory, now parse the info
                        vector<string> words_line ;
                        vector<string> vec1, vec2, vec3 ;
                        bool  check_digit = false ;
                        uint  refline ;
                        for (uint ii=outcarlines.size()-1 ; ii > 1 ; ii--) { // line contents
                            for (uint jj=0 ; jj < words_line.size() ; jj++) {
                                string search_term = "MACROSCOPIC" ;
                                string test_word = words_line.at(jj) ;
                                if(test_word == search_term) { // start of dielectric tensor
                                    refline = ii + 2 ;
                                    check_digit = true ;
                                }
                            }
                            if(check_digit) { // put the tensor info into the string vectors
                                aurostd::string2tokens(outcarlines.at(refline+0),vec1," ") ;
                                aurostd::string2tokens(outcarlines.at(refline+1),vec2," ") ;
                                aurostd::string2tokens(outcarlines.at(refline+2),vec3," ") ;
                                for (uint jj=1 ; jj <= 3 ; jj++) { // string to double, 3x3 matrix, be careful with array bounds
                                    statdiel(1,jj) = atof(vec1.at(jj-1).c_str()) ;
                                    statdiel(2,jj) = atof(vec2.at(jj-1).c_str()) ;
                                    statdiel(3,jj) = atof(vec3.at(jj-1).c_str()) ;
                                }
                                break ;
                            }
                        }
                        if(!check_digit) {
                            cerr << outcar << " lacks MACROSCOPIC statement" << endl ;
                            exit(0) ;
                        } // DONE PARSING //
                        bool matcheck = false ;
                        for (uint ii = 1 ; ii <= 3 ; ii++) { // clean up spuriously small values: e.g. "-0.000001"
                            if( abs(statdiel(1,ii)) < eps ) statdiel(1,ii) = 0.00 ;
                            if( abs(statdiel(2,ii)) < eps ) statdiel(2,ii) = 0.00 ;
                            if( abs(statdiel(3,ii)) < eps ) statdiel(3,ii) = 0.00 ;
                        }
                        for (uint ii = 1 ; ii <= 3 ; ii++) { // check if it is asymmetric & if large off-diags exist
                            for (uint jj = 1 ; jj <= 3 ; jj++) {
                                double testdiff = statdiel[ii][jj] - statdiel[jj][ii] ;
                                if(testdiff >= eps) { // eps is a bit arbitrary right now ..
                                    // serious issues with VASP calculation here: 
                                    cerr << "asymmetric dielectric tensor, exiting..." << endl ;
                                    exit(0) ;
                                } else { // only if small
                                    statdiel(ii,jj) = statdiel(jj,ii) ;
                                }
                                if(ii != jj) {
                                    if(abs(statdiel(ii,jj)) > 0 || abs(statdiel(jj,ii)) > 0) {
                                        matcheck = true ;
                                        break ;
                                    }
                                }
                            }
                            if(matcheck) break ;
                        }
                        matcheck = true ;
                        if(matcheck)
                        { // diagonalize the 3x3 matrix
                            aurostd::eigen(statdiel,eigr,eigi) ;
                        }
                    } // loop GetStatDiel
                } // namespace KBIN

                // ***************************************************************************
                namespace KBIN {
                    void GetDynaDiel(string& outcar, xvector<double>& eigr, xvector<double>& eigi) {
                        int PATH_LENGTH_MAX = 1024 ;
                        char work_dir[PATH_LENGTH_MAX] ;
                        string outcarfile, outcarpath ;
                        string outcarpath_tmp = aurostd::TmpFileCreate("OUTCARc2.tmp") ;
                        vector<string> outcarlines, endline, startline, vasptoken ;
                        xmatrix<double> statdiel(3,3), eigenvec(3,3) ;
                        // double eps = 1.0E-5 ; // need to define this more rigorously
                        getcwd(work_dir, PATH_LENGTH_MAX) ;

                        eigr(1)=0.00 ; eigr(2)=0.00 ; eigr(3)=0.00 ;
                        eigi(1)=0.00 ; eigi(2)=0.00 ; eigi(3)=0.00 ;

                        if(!aurostd::FileExist(outcar)) {
                            cout << "KBIN::GetDynaDiel: check filename || file missing" << endl ;
                            exit(0) ;
                        } else {
                            outcarpath = "/" + outcar ;
                            outcarpath = work_dir + outcarpath ;
                            vector<string> outcardata ;
                            aurostd::string2tokens(outcarpath, outcardata, ".") ;
                            if(outcardata.at(outcardata.size()-1) == "bz2") { // compressed option
                                aurostd::execute("bzcat " + outcarpath + " > " + outcarpath_tmp) ;
                                aurostd::file2vectorstring(outcarpath_tmp, outcarlines) ;
                            } else if(outcardata.at(outcardata.size()-1) == "xz") { // compressed option
                                aurostd::execute("xzcat " + outcarpath + " > " + outcarpath_tmp) ;
                                aurostd::file2vectorstring(outcarpath_tmp, outcarlines) ;
                            } else if(outcardata.at(outcardata.size()-1) == "gz") { // compressed option
                                aurostd::execute("gzcat " + outcarpath + " > " + outcarpath_tmp) ;
                                aurostd::file2vectorstring(outcarpath_tmp, outcarlines) ;
                            } else { // plain text option
                                aurostd::execute("cat " + outcarpath + " > " + outcarpath_tmp) ;
                                aurostd::file2vectorstring(outcarpath_tmp, outcarlines) ;
                            }
                        }
                        exit(0) ;
                        // check the loaded OUTCAR
                        aurostd::string2tokens(outcarlines.at(0),startline," ");
                        aurostd::string2tokens(startline.at(0),vasptoken,".");
                        aurostd::string2tokens(outcarlines.at(outcarlines.size()-1),endline," ");
                        if(vasptoken.at(0) != "vasp" || endline.at(0) != "Voluntary") { // first and last line check
                            cerr << "OUTCAR file is probably corrupt" << endl ;
                            exit(0) ;
                        }
                        uint sec_count = 0 ;
                        for (uint ii=outcarlines.size()-12 ; ii<outcarlines.size() ; ii++) { // presence timing information check
                            vector<string> timetoken ;
                            aurostd::string2tokens(outcarlines.at(ii),timetoken," ") ;
                            if(timetoken.size() > 0) {
                                for (uint jj=0 ; jj<timetoken.size() ; jj++)
                                { if(timetoken.at(jj) == "(sec):") sec_count+=1 ; }
                            }
                        }
                        if(sec_count != 4) { // first and last line check
                            cerr << "OUTCAR file is probably corrupt" << endl ;
                            exit(0) ;
                        }
                        // OUTCAR is now in memory, now parse the info
                        cout << "test" << endl ;
                    } // loop GetDynaDiel
                } // namespace KBIN

                // ***************************************************************************

#endif

                // ***************************************************************************
                // *                                                                         *
                // *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
                // *                                                                         *
                // ***************************************************************************
