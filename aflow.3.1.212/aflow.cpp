// ***************************************************************************
// *                                                                         *
// *         aflow - Automatic FLOW for materials discovery project          *
// *             Stefano Curtarolo - Duke University - 2003-2018             *
// *                                                                         *
// ***************************************************************************
//
//  Copyright 2003-2018 - Stefano Curtarolo - AFLOW.ORG consortium
//
//  This file is part of AFLOW software.
//
//  AFLOW is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
// 
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
// 
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// ***************************************************************************

#include "aflow.h"
#include "aflow_pflow.h"
// [OBSOLETE] #include "aflow_aqe.h"

//#define  __XOPTIMIZE
//#include "aflow_array.h"

// ***************************************************************************

// bool AFLOW_PTHREADS::FLAG;
// int  AFLOW_PTHREADS::MAX_PTHREADS;
// int  AFLOW_PTHREADS::RUNNING;
// pthread_t thread[MAX_ALLOCATABLE_PTHREADS];
// int iret[MAX_ALLOCATABLE_PTHREADS];
// bool thread_busy[MAX_ALLOCATABLE_PTHREADS];

#include "aflow_test.cpp"
#include <sys/ioctl.h>
// #include <linux/kd.h>
//   0x00004B2F   KIOCSOUND     int

namespace aflowlib {
  bool aflowlib2stream(const aflowlib::_aflowlib_entry& data,const string& file,stringstream& stream,bool VERBOSE) {
    return aurostd::url2stringstream(data.aurl+"/"+file,stream,VERBOSE);
  }
}

namespace aflowlib {
  bool aflowlib2file(const aflowlib::_aflowlib_entry& data,const string& file,bool VERBOSE) {
    stringstream stream;
    bool out=aurostd::url2stringstream(data.aurl+"/"+file,stream,VERBOSE);
    aurostd::stringstream2file(stream,file);
    return out;
  }
}

int main(int _argc,char **_argv) {
  string soliloquy="main():"; // CO 180419
  ostream& oss=cout;  // CO 180419
  try{
  bool LDEBUG=FALSE; // TRUE;
  if(LDEBUG) cerr << "AFLOW-MAIN [1]" << endl;
  std::vector<string> argv(aurostd::get_arguments_from_input(_argc,_argv));
  if(LDEBUG) cerr << "AFLOW-MAIN [2]" << endl;
  std::vector<string> cmds;  
  // MACHINE
  init::InitMachine(FALSE,argv,cmds,cerr);    

  atoms_initialize();
  // spacegroup::SpaceGroupInitialize(); only if necessary
  // INFORMATION **************************************************
  AFLOW_PTHREADS::FLAG=AFLOW_PTHREADS::Check_Threads(argv,!XHOST.QUIET);
 
  bool Arun=FALSE;
  if(!Arun && aurostd::args2flag(argv,cmds,"--prx|--prx="))  {Arun=TRUE;PERFORM_PRX(cout);}
  //only for test, aflow --test
  if(!Arun && aurostd::args2flag(argv,cmds,"--test_stefano")) {
    uint y=2017,m=11;
    m+=1;
    for(uint i=0;i<200;i++) {
      if(m==0) {
 	cout << "mv \"unknown.pdf\" stefano_" << y << m << ".pdf" << endl;
      } else {
	if(m<10) {
	  cout << "mv \"unknown(" << i << ").pdf\" stefano_" << y << "0" << m << ".pdf" << endl;
	} else {
	  cout << "mv \"unknown(" << i << ").pdf\" stefano_" << y << m << ".pdf" << endl;
	}
      }
      m--;
      if(m==0) {y--;m+=12;} 
    }
      //exit(0);
      return 0; // CO 180419
   }
  /*
  if(!Arun && aurostd::args2flag(argv,cmds,"--test")) {
    deque<string> vext; aurostd::string2tokens(".bz2,.xz,.gz",vext,",");vext.push_front("");
    deque<string> vcat; aurostd::string2tokens("cat,bzcat,xzcat,gzcat",vcat,",");
    if(vext.size()!=vcat.size()) { cerr << "ERROR - aflow.cpp:main: vext.size()!=vcat.size(), aborting." << endl; exit(0); }
	
    for(uint iext=0;iext<vext.size();iext++) { 
      cout << "\"" << vext.at(iext) << "\"" << " " << "\"" << vcat.at(iext) << "\"" << endl;
    }
    
    int dim=7;
    cout << "dimm=" << dim << endl;

    xmatrix<double> m(dim,dim),mi(dim,dim);
    for(int i=1;i<=dim;i++) 
      for(int j=1;j<=dim;j++)
	m(i,j)=aurostd::ran0();
    cout << "m=" << endl << m << endl;
    mi=inverse(m);
    cout << "mi=" << endl << mi << endl;
    cout << "mi*m=" << endl << det(mi*m) << endl; 
     
    string test="2.730747137  -2.730747137-12.397646334";
    vector<string> _tokens;
    aurostd::string2tokens(test,_tokens,"-");
    for(uint i=0;i<_tokens.size();i++){
      cerr << _tokens[i] << endl;
    }
      //exit(0);
      return 0; // CO 180419
    // COREY START 170614 - some SQLITE tests
    //http://zetcode.com/db/sqlitec/ - more tests here
    //this will create test.db file
    sqlite3 *db;
    char *err_msg = 0;
    int rc = sqlite3_open("test.db", &db);
    if(rc != SQLITE_OK) {
        fprintf(stderr, "Cannot open database: %s\n", sqlite3_errmsg(db));
        sqlite3_close(db);
        return 1;
    }
    char const *sql = "DROP TABLE IF EXISTS Cars;" 
                "CREATE TABLE Cars(Id INT, Name TEXT, Price INT);" 
                "INSERT INTO Cars VALUES(1, 'Audi', 52642);" 
                "INSERT INTO Cars VALUES(2, 'Mercedes', 57127);" 
                "INSERT INTO Cars VALUES(3, 'Skoda', 9000);" 
                "INSERT INTO Cars VALUES(4, 'Volvo', 29000);" 
                "INSERT INTO Cars VALUES(5, 'Bentley', 350000);" 
                "INSERT INTO Cars VALUES(6, 'Citroen', 21000);" 
                "INSERT INTO Cars VALUES(7, 'Hummer', 41400);" 
                "INSERT INTO Cars VALUES(8, 'Volkswagen', 21600);";
    rc = sqlite3_exec(db, sql, 0, 0, &err_msg);
    if(rc != SQLITE_OK ) {
        fprintf(stderr, "SQL error: %s\n", err_msg);
        sqlite3_free(err_msg);        
        sqlite3_close(db);
        return 1;
    } 
    sqlite3_close(db);
    //return 0;

    //MORE TESTS
    //printf("%s\n,sqlite3_libversion()");
    //    sqlite3 *db;
    //sqlite3_stmt *res;
    //int rc = sqlite3_open(":memory:", &db);
    //if(rc != SQLITE_OK) {
    //    fprintf(stderr, "Cannot open database: %s\n", sqlite3_errmsg(db));
    //    sqlite3_close(db);
    //    return 1;
    //}
    //rc = sqlite3_prepare_v2(db, "SELECT SQLITE_VERSION()", -1, &res, 0);    
    //if(rc != SQLITE_OK) {
    //    fprintf(stderr, "Failed to fetch data: %s\n", sqlite3_errmsg(db));
    //    sqlite3_close(db);
    //    return 1;
    //}    
    //rc = sqlite3_step(res);
    //if(rc == SQLITE_ROW) {
    //    printf("%s\n", sqlite3_column_text(res, 0));
    //}
    //sqlite3_finalize(res);
    //sqlite3_close(db);
    //return 0;
    
    //quick easy test
    cerr << sqlite3_libversion() << endl;
    // COREY END 170614 - some SQLITE tests
    aurostd::xcomplex<double> x(123.0,456.0);
    cout << x.re << "," << x.im << " - " << x.real() << "," << x.imag() << " - " << x << endl;
    x.re=111;x.im=222;
    cout << x.re << "," << x.im << " - " << x.real() << "," << x.imag() << " - " << x << endl;
    cout << aurostd::PaddedPOST("EMIN= -30.0",10) << endl;;
    stringstream for_corey;
    for_corey << "scatter/use mapped color={draw=black,fill=mapped color,solid}";
    string corey=for_corey.str();
    cout << corey << endl;
    stringstream aus;
    aus << "************************   00000  MESSAGE KPOINTS KSHIFT=[" << 1 << " " << 2 << " " << 3 << "]" << " ************************ " << endl;
    cout << aus.str() << endl;
    //exit(0);
    return 0; // CO 180419
  }
  */
  
  if(!Arun && aurostd::args2attachedflag(argv,"--bin2base64=")) {
    string b64String;
    aurostd::bin2base64(aurostd::args2attachedstring(argv,"--bin2base64=",""),b64String);
    cout << b64String << endl;
    exit(0);
  }

  if(!Arun && aurostd::args2flag(argv,cmds,"--test=POTCAR|--test=POTCAR.relax1"+DEFAULT_KZIP_EXT+"|--test=POTCAR.relax2"+DEFAULT_KZIP_EXT+"|--test=POTCAR.static"+DEFAULT_KZIP_EXT+"|--test=POTCAR.bands"+DEFAULT_KZIP_EXT+"")) {
      XHOST.DEBUG=TRUE;xPOTCAR(aurostd::args2attachedstring(argv,"--test=",""));/*exit(0)*/return 0;} // CO 180419
  if(!Arun && aurostd::args2flag(argv,cmds,"--test=DOSCAR|--test=DOSCAR.relax1"+DEFAULT_KZIP_EXT+"|--test=DOSCAR.relax2"+DEFAULT_KZIP_EXT+"|--test=DOSCAR.static"+DEFAULT_KZIP_EXT+"|--test=DOSCAR.bands"+DEFAULT_KZIP_EXT+"")) {
      XHOST.DEBUG=TRUE;xDOSCAR(aurostd::args2attachedstring(argv,"--test=",""));/*exit(0)*/return 0;} // CO 180419
  if(!Arun && aurostd::args2flag(argv,cmds,"--test=EIGENVAL|--test=EIGENVAL.relax1"+DEFAULT_KZIP_EXT+"|--test=EIGENVAL.relax2"+DEFAULT_KZIP_EXT+"|--test=EIGENVAL.static"+DEFAULT_KZIP_EXT+"|--test=EIGENVAL.bands"+DEFAULT_KZIP_EXT+"")) {
      XHOST.DEBUG=TRUE;xEIGENVAL(aurostd::args2attachedstring(argv,"--test=",""));/*exit(0)*/return 0;} // CO 180419
  if(!Arun && aurostd::args2flag(argv,cmds,"--test=OUTCAR|--test=OUTCAR.relax1"+DEFAULT_KZIP_EXT+"|--test=OUTCAR.relax2"+DEFAULT_KZIP_EXT+"|--test=OUTCAR.static"+DEFAULT_KZIP_EXT+"|--test=OUTCAR.bands"+DEFAULT_KZIP_EXT+"")) {
      XHOST.DEBUG=TRUE;xOUTCAR(aurostd::args2attachedstring(argv,"--test=",""));/*exit(0)*/return 0;} // CO 180419
  if(!Arun && aurostd::args2flag(argv,cmds,"--test=KPOINTS|--test=KPOINTS.relax1"+DEFAULT_KZIP_EXT+"|--test=KPOINTS.relax2"+DEFAULT_KZIP_EXT+"|--test=KPOINTS.static"+DEFAULT_KZIP_EXT+"|--test=KPOINTS.bands"+DEFAULT_KZIP_EXT+"")) {
      XHOST.DEBUG=TRUE;xKPOINTS(aurostd::args2attachedstring(argv,"--test=",""));/*exit(0)*/return 0;}  // CO 180419
  if(!Arun && aurostd::args2flag(argv,cmds,"--test=vasprun|--test=vasprun.xml.relax1"+DEFAULT_KZIP_EXT+"|--test=vasprun.xml.relax2"+DEFAULT_KZIP_EXT+"|--test=vasprun.xml.static"+DEFAULT_KZIP_EXT+"|--test=vasprun.xml.bands"+DEFAULT_KZIP_EXT+"")) {
      XHOST.DEBUG=TRUE;xVASPRUNXML(aurostd::args2attachedstring(argv,"--test=",""));/*exit(0)*/return 0;} // CO 180419
  if(!Arun && aurostd::args2flag(argv,cmds,"--test=IBZKPT|--test=IBZKPT.relax1"+DEFAULT_KZIP_EXT+"|--test=IBZKPT.relax2"+DEFAULT_KZIP_EXT+"|--test=IBZKPT.static"+DEFAULT_KZIP_EXT+"|--test=IBZKPT.bands"+DEFAULT_KZIP_EXT+"")) {
      XHOST.DEBUG=TRUE;xIBZKPT(aurostd::args2attachedstring(argv,"--test=",""));/*exit(0)*/return 0;} // CO 180419
 if(!Arun && aurostd::args2flag(argv,cmds,"--test=CHGCAR|--test=CHGCAR.relax1"+DEFAULT_KZIP_EXT+"|--test=CHGCAR.relax2"+DEFAULT_KZIP_EXT+"|--test=CHGCAR.static"+DEFAULT_KZIP_EXT+"|--test=CHGCAR.bands"+DEFAULT_KZIP_EXT+"")) {
      XHOST.DEBUG=TRUE;xCHGCAR(aurostd::args2attachedstring(argv,"--test=",""));/*exit(0)*/return 0;} // CO 180419

  if(!Arun && aurostd::args2flag(argv,cmds,"--test_proto1")) {
    vector<xstructure> vstr;
    vector<string> vlattice;aurostd::string2tokens("BCC,BCT,CUB,FCC,HEX,MCL,MCLC,ORC,ORCC,ORCF,ORCI,RHL,TET,TRI",vlattice,",");
    aflowlib::_aflowlib_entry data;
    vector<aflowlib::_aflowlib_entry> vdata;
    for(uint i=4;i<vlattice.size();i++) {
      data.clear();
      data.url2aflowlib("materials.duke.edu:AFLOWDATA/ICSD_WEB/"+vlattice.at(i)+"/?format=text",cout,FALSE);vdata.push_back(data);
      cout << "AFLOWLIB " << vlattice.at(i) << "=" << data.vaflowlib_entries.size() << endl;
      for(uint j=0;j<data.vaflowlib_entries.size();j++) {
	aflowlib::_aflowlib_entry dataj;
	dataj.url2aflowlib("materials.duke.edu:AFLOWDATA/ICSD_WEB/"+vlattice.at(i)+"/"+data.vaflowlib_entries.at(j),cout,TRUE);
       	aurostd::StringSubst(dataj.aurl,"aflowlib","materials");
	if(dataj.aurl!="") {
	  xstructure str(dataj.aurl,"CONTCAR.relax.vasp",IOAFLOW_AUTO);
	  xEIGENVAL xEIGENVAL; xEIGENVAL.GetPropertiesUrlFile(dataj.aurl,"EIGENVAL.bands"+DEFAULT_KZIP_EXT+"",FALSE);
	  xOUTCAR xOUTCAR; xOUTCAR.GetPropertiesUrlFile(dataj.aurl,"OUTCAR.static"+DEFAULT_KZIP_EXT+"",FALSE);
	  xDOSCAR xDOSCAR; xDOSCAR.GetPropertiesUrlFile(dataj.aurl,"DOSCAR.static"+DEFAULT_KZIP_EXT+"",FALSE);
	  // if(aurostd::args2flag(argv,cmds,"--vasp")) aurostd::url2stringstream(dataj.aurl+"/CONTCAR.relax.vasp",stream);
	  // if(aurostd::args2flag(argv,cmds,"--qe")) aurostd::url2stringstream(dataj.aurl+"/CONTCAR.relax.qe",stream);
	  // if(aurostd::args2flag(argv,cmds,"--abinit")) aurostd::url2stringstream(dataj.aurl+"/CONTCAR.relax.abinit",stream);
	  // if(aurostd::args2flag(argv,cmds,"--aims")) aurostd::url2stringstream(dataj.aurl+"/CONTCAR.relax.aims",stream);
	  vstr.push_back(str);
	  cerr << "vstr.size()=" << vstr.size() << "  "
	       << "str.atoms.size()=" << str.atoms.size() << "  "
	       << "OUTCAR.static"+DEFAULT_KZIP_EXT+".size()=" << xOUTCAR.vcontent.size() << "  "
	       << "DOSCAR.static"+DEFAULT_KZIP_EXT+".size()=" << xDOSCAR.vcontent.size() << "  "
	       << "EIGENVAL.static"+DEFAULT_KZIP_EXT+".size()=" << xEIGENVAL.vcontent.size() << "  "
	       << endl;
	  //	  cerr << str << endl;
	}
      }
    }
      //exit(0);
      return 0; // CO 180419
  }

  if(!Arun && aurostd::args2flag(argv,cmds,"--testJ")) {Arun=TRUE;PERFORM_TESTJ(cout);}
  if(!Arun && aurostd::args2flag(argv,cmds,"--test1")) {Arun=TRUE;PERFORM_TEST1(cout);}
  if(!Arun && aurostd::args2flag(argv,cmds,"--test3")) {Arun=TRUE;PERFORM_TEST3(cout);}
  if(!Arun && XHOST.vflag_control.flag("MACHINE"))  {Arun=TRUE;init::InitMachine(TRUE,argv,cmds,cout);}
  
  // **************************************************************
  // INTERCEPT AFLOW
  if(!Arun && XHOST.vflag_control.flag("SWITCH_AFLOW")) {Arun=TRUE;AFLOW_main(argv);}
  // DX
  if(!Arun && XHOST.vflag_control.flag("AFLOWIN_SYM")) {Arun=TRUE;AFLOW_main(argv);} 
  // DX
  if(!Arun && (XHOST.vflag_aflow.flag("CLEAN") || XHOST.vflag_aflow.flag("XCLEAN") || XHOST.AFLOW_RUNDIRflag || XHOST.AFLOW_MULTIflag || XHOST.AFLOW_RUNXflag)) {
    Arun=TRUE;AFLOW_main(argv);
  }
  //  // **************************************************************
  // // INTERCEPT AFLOWLIB
  // MOVED INSIDE PFLOW
  // if(!Arun && XHOST.vflag_control.flag("SWITCH_AFLOWLIB") && !XHOST.vflag_pflow.flag("PROTOS") && !XHOST.vflag_pflow.flag("PROTOS_ICSD") && !XHOST.vflag_pflow.flag("PROTO")) {Arun=TRUE;aflowlib::WEB_Aflowlib_Entry_PHP(argv,cout);}
  // **************************************************************
  // INTERCEPT APENNSY
  if(!Arun && XHOST.vflag_control.flag("SWITCH_APENNSY1")) {Arun=TRUE;Apennsymain(argv,cmds);}
  if(!Arun && XHOST.vflag_control.flag("SWITCH_APENNSY2")) {Arun=TRUE;Apennsymain(argv,cmds);}
  
  // **************************************************************
  // INTERCEPT aconvasp/aqe/apennsy by title
  if(!Arun && aurostd::substring2bool(XHOST.Progname,"aconvasp","convasp")) {Arun=TRUE;pflow::main(argv,cmds);}
  if(!Arun && aurostd::substring2bool(XHOST.Progname,"aqe")) {Arun=TRUE;pflow::main(argv,cmds);}
  if(!Arun && aurostd::substring2bool(XHOST.Progname,"apennsy")) {Arun=TRUE;Apennsymain(argv,cmds);}
  
  // **************************************************************
  // intercept commands
    if(!Arun && XHOST.vflag_control.flag("MULTI=SH")) {Arun=TRUE;AFLOW_PTHREADS::MULTI_sh(argv);/*exit(0)*/return 0;}  // CO 180419
    if(!Arun && XHOST.vflag_control.flag("MULTI=BZIP2")) {Arun=TRUE;AFLOW_PTHREADS::MULTI_compress("bzip2",argv);/*exit(0)*/return 0;}  // CO 180419
    if(!Arun && XHOST.vflag_control.flag("MULTI=BUNZIP2")) {Arun=TRUE;AFLOW_PTHREADS::MULTI_compress("bunzip2",argv);/*exit(0)*/return 0;}  // CO 180419
    if(!Arun && XHOST.vflag_control.flag("MULTI=GZIP")) {Arun=TRUE;AFLOW_PTHREADS::MULTI_compress("gzip",argv);/*exit(0)*/return 0;}  // CO 180419
    if(!Arun && XHOST.vflag_control.flag("MULTI=GUNZIP")) {Arun=TRUE;AFLOW_PTHREADS::MULTI_compress("gunzip",argv);/*exit(0)*/return 0;}  // CO 180419
    if(!Arun && XHOST.vflag_control.flag("MULTI=XZIP")) {Arun=TRUE;AFLOW_PTHREADS::MULTI_compress("xz",argv);/*exit(0)*/return 0;}  // CO 180419
    if(!Arun && XHOST.vflag_control.flag("MULTI=XUNZIP")) {Arun=TRUE;AFLOW_PTHREADS::MULTI_compress("xunzip",argv);/*exit(0)*/return 0;}  // CO 180419
    if(!Arun && XHOST.vflag_control.flag("MULTI=BZ2XZ")) {Arun=TRUE;AFLOW_PTHREADS::MULTI_bz2xz(argv);/*exit(0)*/return 0;}  // CO 180419
    if(!Arun && XHOST.vflag_control.flag("MULTI=GZ2XZ")) {Arun=TRUE;AFLOW_PTHREADS::MULTI_gz2xz(argv);/*exit(0)*/return 0;}  // CO 180419
    if(!Arun && XHOST.vflag_control.flag("MULTI=ZIP")) {Arun=TRUE;AFLOW_PTHREADS::MULTI_zip(argv);/*exit(0)*/return 0;}  // CO 180419
    if(!Arun && XHOST.vflag_control.flag("MONITOR")) {Arun=TRUE;AFLOW_monitor(argv);/*exit(0)*/return 0;} // CO 180419
    if(!Arun && XHOST.vflag_control.flag("GETTEMP")) {Arun=TRUE;AFLOW_getTEMP(argv);/*exit(0)*/return 0;} // CO 180419

  // **************************************************************
  // INTERCEPT HELP
  if(XHOST.vflag_control.flag("AFLOW_HELP")) {
    cout << aflow::Banner("BANNER_BIG") << endl << aflow::Intro_HELP("aflow") << aflow::Banner("BANNER_BIG") << endl;
      /*exit(1)*/return 0;} // << endl; // CO 180419
  if(XHOST.vflag_control.flag("AFLOW_EXCEPTIONS")) {
    cout << aflow::Banner("BANNER_BIG") << endl << aflow::Banner("EXCEPTIONS") << endl;
    return 0;
  }  // ME180531
  if(XHOST.vflag_control.flag("README_AFLOW_LICENSE_GPL3"))  {
    cout << aflow::License_Preamble_aflow() << endl;
    cout << " " << endl;
    cout << init::InitGlobalObject("README_AFLOW_LICENSE_GPL3_TXT") << endl;
    cout << " " << endl;
    cout << "*************************************************************************** " << endl;
      /*exit(1)*/return 0;} // << endl; // CO 180419
  if(XHOST.vflag_control.flag("README_AFLOW"))  {
    cout << init::InitGlobalObject("README_AFLOW_TXT") << endl;
      /*exit(1)*/return 0;} // << endl; // CO 180419
  if(XHOST.vflag_control.flag("README_AFLOW_VERSIONS_HISTORY"))  {
    cout << init::InitGlobalObject("README_AFLOW_VERSIONS_HISTORY_TXT") << endl;
      /*exit(1)*/return 0;} // << endl; // CO 180419
  if(XHOST.vflag_control.flag("README_AFLOW_PFLOW"))  {
    cout << init::InitGlobalObject("README_AFLOW_PFLOW_TXT") << endl;
      /*exit(1)*/return 0;} // << endl; // CO 180419
  if(XHOST.vflag_control.flag("README_FROZSL"))  {
    cout << init::InitGlobalObject("README_AFLOW_FROZSL_TXT") << endl;
      /*exit(1)*/return 0;} // << endl; // CO 180419
  if(XHOST.vflag_control.flag("README_APL"))  {
    cout << init::InitGlobalObject("README_AFLOW_APL_TXT") << endl;
      /*exit(1)*/return 0;} // << endl; // CO 180419
  if(XHOST.vflag_control.flag("README_QHA"))  {
    cout << init::InitGlobalObject("README_AFLOW_QHA_SCQHA_QHA3P_TXT") << endl;
      /*exit(1)*/return 0;} // << endl; // CO 180419
  if(XHOST.vflag_control.flag("README_AAPL"))  {
    cout << init::InitGlobalObject("README_AFLOW_APL_TXT") << endl;
      /*exit(1)*/return 0;} // << endl; // CO 180419
  if(XHOST.vflag_control.flag("README_AGL"))  {
    cout << init::InitGlobalObject("README_AFLOW_AGL_TXT") << endl;
      /*exit(1)*/return 0;} // << endl; // CO 180419
  if(XHOST.vflag_control.flag("README_AEL"))  {
    cout << init::InitGlobalObject("README_AFLOW_AEL_TXT") << endl;
      /*exit(1)*/return 0;} // << endl; // CO 180419
  if(XHOST.vflag_control.flag("README_ANRL"))  {
    cout << init::InitGlobalObject("README_AFLOW_ANRL_TXT") << endl;
      /*exit(1)*/return 0;} // << endl; // CO 180419
  if(XHOST.vflag_control.flag("README_COMPARE"))  {
    cout << init::InitGlobalObject("README_AFLOW_COMPARE_TXT") << endl;
      /*exit(1)*/return 0;} // << endl; // CO 180419
  if(XHOST.vflag_control.flag("README_SYMMETRY"))  {
    cout << init::InitGlobalObject("README_AFLOW_SYM_TXT") << endl;
      /*exit(1)*/return 0;} // << endl; // CO 180419
  if(XHOST.vflag_control.flag("README_CHULL"))  {
    cout << init::InitGlobalObject("README_AFLOW_CHULL_TXT") << endl;
      /*exit(1)*/return 0;} // << endl; // CO 180419
  if(XHOST.vflag_control.flag("README_PARTIAL_OCCUPATION")) {
    cout << init::InitGlobalObject("README_AFLOW_POCC_TXT") << endl;
      /*exit(1)*/return 0;} // << endl; // CO 180419
  if(XHOST.vflag_control.flag("README_APENNSY"))  {
    cout << init::InitGlobalObject("README_AFLOW_APENNSY_TXT") << endl;
      /*exit(1)*/return 0;} // << endl; // CO 180419
  if(XHOST.vflag_control.flag("README_SCRIPTING"))  {
    cout << init::InitGlobalObject("README_AFLOW_SCRIPTING_TXT") << endl;
      /*exit(1)*/return 0;} // << endl; // CO 180419
  if(XHOST.vflag_control.flag("README_EXCEPTIONS"))  {
    cout << init::InitGlobalObject("README_AFLOW_EXCEPTIONS_TXT") << endl;
    return 0;}  // ME 180531
  if(XHOST.vflag_control.flag("README_HTRESOURCES"))  {
    cout << init::InitGlobalObject("README_AFLOW_HTRESOURCES_TXT") << endl;
      /*exit(1)*/return 0;} // << endl; // CO 180419
  if(XHOST.vflag_control.flag("README_XAFLOW"))  {
    cout << init::InitGlobalObject("README_AFLOW_XAFLOW_TXT") << endl;
      /*exit(1)*/return 0;} // << endl; // CO 180419
  if(XHOST.vflag_control.flag("README_AFLOWRC"))  {
    cout << init::InitGlobalObject("README_AFLOW_AFLOWRC_TXT") << endl;
      /*exit(1)*/return 0;} // << endl; // CO 180419
 
  // **************************************************************
  // PHP-WEB AND CURRICULUM AND HIGH-THROUGHPUT STUFF
  ProcessPhpLatexCv();
  //  ProcessSecurityOptions(argv,cmds); OLD STUFF AFLOW SECURITY
  
  // **************************************************************
  bool VVERSION=aurostd::args2flag(argv,cmds,"-v|--version");
    if(!Arun && VVERSION)  {Arun=TRUE; cout << aflow::Banner("AFLOW_VERSION");/*exit(0)*/return 0;} // look for version IMMEDIATELY // CO 180419
    if(!Arun && XHOST.TEST) { Arun=TRUE;cerr << "test" << endl; /*exit(0)*/return 0;} // CO 180419
  
  if(!Arun && XHOST.argv.size()==1 && (aurostd::substring2bool(XHOST.Progname,"aflow")  || aurostd::substring2bool(XHOST.Progname,"aflowd"))) {   
    //   Arun=TRUE;AFLOW_main(argv);
    Arun=TRUE;
    //    cout << "******************************************************************************************************" << endl;
    //  cout << aflow::Banner("BANNER_TINY") << endl;
    cout << aflow::Banner("BANNER_BIG") << endl;
    cout << aflow::Intro_aflow("aflow") << endl;
    cout << pflow::Intro_pflow("aflow") << endl;
    cout << aflow::Intro_sflow("aflow") << endl;
    cout << aflow::Intro_HELP("aflow") << endl;
    cout << aflow::Banner("BANNER_BIG") << endl;
    //    cout << "XHOST.argv.size()=" << XHOST.argv.size()<< endl;
    //    cout << "******************************************************************************************************" << endl;
  }

  // **************************************************************
  // LAST RESOURCE PFLOW
  if(!Arun) { 
    Arun=TRUE;pflow::main(argv,cmds);
  }
  // **************************************************************
  // END
  return (int) !Arun; //Arun==TRUE is 1, so flip because return 0 is normal
}
  // CO 180729 - OBSOLETE - use xerror
  //[OBSOLETE]catch(AFLOWRuntimeError& re){
  //[OBSOLETE]  pflow::logger(soliloquy, "AFLOWRuntimeError detected. Report on the AFLOW Forum: aflow.org/forum.", oss, _LOGGER_ERROR_);
  //[OBSOLETE]  pflow::logger(re.where(), re.what(), oss, _LOGGER_ERROR_);
  //[OBSOLETE]  return 1;
  //[OBSOLETE]}
  //[OBSOLETE]catch(AFLOWLogicError& le){
  //[OBSOLETE]  pflow::logger(soliloquy, "AFLOWLogicError detected. Adjust your inputs accordingly.", oss, _LOGGER_ERROR_);
  //[OBSOLETE]  pflow::logger(le.where(), le.what(), oss, _LOGGER_ERROR_);
  //[OBSOLETE]  return 1;
  //[OBSOLETE]}
  catch (aurostd::xerror& excpt) {
    pflow::logger(excpt.where(), excpt.error_message, oss, _LOGGER_ERROR_);
    return excpt.error_code;
  }
}


// ***************************************************************************
// AFLOW_main
// ***************************************************************************
int AFLOW_main(vector<string> &argv) {
  if(!XHOST.QUIET) cout << aflow::Banner("INTRODUCTION");// << endl;
  KBIN::KBIN_Main(argv);
  return 0; 
}

// ***************************************************************************
// aflow::License_aflow
// ***************************************************************************
namespace aflow {
  string License_Preamble_aflow(void) {
    //( C) 2003-2018 Stefano Curtarolo, MIT-Duke University stefano@duke.edu
    stringstream strstream;
    strstream << "\
   \n\
   *************************************************************************** \n\
   *                                                                         * \n\
   *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           * \n\
   *                                                                         * \n\
   *************************************************************************** \n\
    Copyright 2003-2018 - Stefano Curtarolo - AFLOW.ORG consortium \n\
    \n\
    This file is part of AFLOW software. \n\
    \n\
    AFLOW is free software: you can redistribute it and/or modify \n\
    it under the terms of the GNU General Public License as published by \n\
    the Free Software Foundation, either version 3 of the License, or \n\
    (at your option) any later version. \n\
    \n\
    This program is distributed in the hope that it will be useful, \n\
    but WITHOUT ANY WARRANTY; without even the implied warranty of \n\
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the \n\
    GNU General Public License for more details. \n\
    \n\
    You should have received a copy of the GNU General Public License \n\
    along with this program.  If not, see <http://www.gnu.org/licenses/>. \n\
    \n\
   *************************************************************************** \n\
 " << endl;
    return strstream.str();
  }
} // namespace aflow


// ***************************************************************************
// aflow::Intro_aflow
// ***************************************************************************
namespace aflow {
  string Intro_aflow(string x) {
    //( C) 2003-2018 Stefano Curtarolo, MIT-Duke University stefano@duke.edu
    stringstream strstream;
    strstream << "\
   \n\
******* BEGIN INFORMATION MODE ********************************************************************* \n\
  "<< x<<" -h | --help | --readme_aflow   CHECK \n\
  "<< x<<" --version | -v                 CHECK \n\
  "<< x<<" --machine                      CHECK \n\
******* END INFORMATION MODE *********************************************************************** \n\
  \n\
******* BEGIN RUNNING MODE ************************************************************************* \n\
  "<< x<<" --run | --run=multi | --run=N   COMMANDS for running \n\
  "<< x<<" --clean | -c                    COMMAND for cleaning   \n\
  \n\
  MODIFIERS \n\
   --DIRECTORY[=| ]dir | --D[=| ]dir | --d[=| ]dir \n\
   --FILE[=| ]file | --F[=| ]file | --f[=| ]file  \n\
   --quiet | -quiet | -q  \n\
   --loop  \n\
   --sort | -sort \n\
   --reverse | -rsort \n\
   --random | -rnd \n\
   --force | -force \n\
   --mem=XX | --maxmem=XX \n\
   --readme= xaflow | aflow | aconvasp | aflowrc | scripting | apennsy | apl | agl | ael | anrl | compare | symmetry | chull | errors | exceptions | frozsl CHECK !!!! \n\
   --np=NUMBER |   --npmax     \n\
   --generate_aflowin_from_vasp \n\
   --generate_vasp_from_aflowin | --generate \n\
   --use_aflow.in=XXX \n\
   --use_LOCK=XXX \n\
 \n\
  MODIFIERS MPI/SERIAL PARAMETERS \n\
   --mpi | -nompi | --serial  \n\
 \n\
  HOST ORIENTED OPTION \n\
   --machine=beta | beta_openmpi | qrats | qflow | conrad | eos | materials | habana | aflowlib | ranger | kraken \n\
             marylou | parsons | jellium | ohad | host1 \n\
             raptor --np=N | diamond --np=N \n\
******* END RUNNING MODE *************************************************************************** \n\
 " << endl;
    // --readme=htresources
    return strstream.str();
  }
} // namespace aflow

// ***************************************************************************
// aflow::Intro_sflow
// ***************************************************************************
namespace aflow {
  string Intro_sflow(string x) {
    stringstream strstream;
    strstream << "\
******* BEGIN SCRIPTING MODE *********************************************************************** \n\
 AFLOW SCRIPTING COMMANDS \n\
  "<< x<<" --justafter=string \n\
  "<< x<<" --justbefore=string \n\
  "<< x<<" --justbetween=string_from[,string_to] \n\
  "<< x<<" --qsub=N,file \n\
  "<< x<<" --qdel=aaa,nnn:mmm,aaa,bbb,ccc \n\
  "<< x<<" --bsub=N,file \n\
  "<< x<<" --bkill=aaa,nnn:mmm,aaa,bbb,ccc \n\
  "<< x<<" --sbatch=N,file \n\
  "<< x<<" --scancel=aaa,nnn:mmm,aaa,bbb,ccc \n\
  "<< x<<" --kill=aaa,nnn:mmm,aaa,bbb,ccc \n\
  "<< x<<" --multi=sh [--np=NUMBER | npmax | nothing] [--F[ILE]] file \n\
  "<< x<<" --multi=zip [--prefix PREFIX] [--size SSSS] --F[ILE] file | --D[IRECTORY directory1 directory2 .... \n\
  "<< x<<" --multi=bzip2 [--np=NUMBER | npmax | nothing] --F[ILE] file1 file2 file3 ....  \n\
  "<< x<<" --multi=bunzip2 [--np=NUMBER | npmax | nothing] --F[ILE] file1.bz2 file2.bz2 file3.bz2 .... \n\
  "<< x<<" --multi=gzip [--np=NUMBER | npmax | nothing] --F[ILE] file1 file2 file3 ....  \n\
  "<< x<<" --multi=gunzip [--np=NUMBER | npmax | nothing] --F[ILE] file1.gz file2.gz file3.gz .... \n\
  "<< x<<" --multi=xzip [--np=NUMBER | npmax | nothing] --F[ILE] file1 file2 file3 ....  \n\
  "<< x<<" --multi=xunzip [--np=NUMBER | npmax | nothing] --F[ILE] file1.xz file2.xz file3.xz .... \n\
  "<< x<<" --multi=bz2xz [--np=NUMBER | npmax | nothing] --F[ILE] file1.bz2 file2.bz2 file3.bz2 .... \n\
  "<< x<<" --multi=gz2xz [--np=NUMBER | npmax | nothing] --F[ILE] file1.gz file2.gz file3.gz .... \n\
  "<< x<<" --getTEMP [--runstat | --runbar | --refresh=X | --warning_beep=T | --warning_halt=T | --mem=XX ] \n\
  "<< x<<" --monitor [--mem=XX] \n\
******* END SCRIPTING MODE ************************************************************************* \n\
 " << endl;
    return strstream.str();
  }
} // namespace aflow

// ***************************************************************************
// aflow::Intro_HELP
// ***************************************************************************
namespace aflow {
  string Intro_HELP(string x) {
    stringstream strstream;
    strstream << "\
******* BEGIN HELP MODE **************************************************************************** \n\
 AFLOW HELP AVAILABLE HELPS  \n\
  "<< x<<" --license \n\
           License information.  \n\
  "<< x<<" --help \n\
           This help.  \n\
  "<< x<<" --readme \n\
           The list of all the commands available.  \n\
  "<< x<<" --readme=xaflow \n\
           Returns the HELP information for the installation of aflow. \n\
  "<< x<<" --readme=aflow | --readme=run | --readme_aflow  \n\
           Returns the HELP information for the \"running machinery\".  \n\
  "<< x<<" --readme=aflowrc \n\
           Returns the HELP information for the installation of aflow. \n\
  "<< x<<" --readme=pflow | --readme=processor | --readme=aconvasp | --readme_aconvasp \n\
           Returns the HELP information for the \"processing machinery\".  \n\
  "<< x<<" --readme=scripting | --readme_scripting \n\
           Returns the HELP information for the \"scripting\" operations.  \n\
  "<< x<<" --readme=apl | --readme_apl \n\
           Returns the HELP information for the \"aflow-harmonic-phonon-library\".  \n\
  "<< x<<" --readme=qha | --readme_qha | --readme=qha3p | --readme_qha3p | --readme=scqha | --readme_scqha \n\
           Returns the HELP information for the \"aflow-quasi-harmonic-library\".  \n\
  "<< x<<" --readme=aapl | --readme_aapl \n\
           Returns the HELP information for the \"aflow-anharmonic-phonon-library (AFLOW-AAPL)\".  \n\
  "<< x<<" --readme=agl | --readme_agl \n\
           Returns the HELP information for the \"aflow-gibbs-library (AFLOW-AGL)\".  \n\
  "<< x<<" --readme=ael | --readme_ael \n\
           Returns the HELP information for the \"aflow-elastic-library (AFLOW-AEL)\".  \n\
  "<< x<<" --readme=anrl | --readme_anrl \n\
           Returns the HELP information for the \"aflow library of prototypes\".  \n\
  "<< x<<" --readme=compare | --readme_compare \n\
           Returns the HELP information for the \"structure geometry comparison code\".  \n\
  "<< x<<" --readme=symmetry | --readme_symmetry \n\
           Returns the HELP information for the \"symmetry library (AFLOW-SYM)\".  \n\
  "<< x<<" --readme=chull | --readme_chull \n\
           Returns the HELP information for the \"convex hull library (AFLOW-hull)\".  \n\
  "<< x<<" --readme=errors | --readme=exceptions | --readme_errors | --readme_exceptions \n\
           Returns the HELP information for exception handling in AFLOW.  \n\
  "<< x<<" --readme=partial_occupation | --readme=pocc | --readme_pocc \n\
           Returns the HELP information for the \"partial occupation library\".  \n\
  "<< x<<" --readme=apennsy | --readme_apennsy \n\
           Returns the HELP information for the \"apennsy\" binary phase diagram generator.  \n\
  "<< x<<" --readme=frozsl|--readme_frozsl \n\
           Returns the HELP information for the \"frozsl\" add ons. \n\
******* END HELP MODE ****************************************************************************** \n\
 " << endl;
    return strstream.str();
  }
} // namespace aflow


// ***************************************************************************
// aflow::Banner
// ***************************************************************************
namespace aflow {
  string Banner(string type) {
    stringstream oss;
    if(type=="VERSION" || type=="AFLOW_VERSION") {
      oss << "" << endl;
      oss << "AFLOW VERSION " << string(AFLOW_VERSION) << " Automatic-FLOW [(C) "<<XHOST.Copyright_Years<<" aflow.org consortium]" << endl;
      oss << "New versions are available here: <http://" << XHOST.AFLOW_MATERIALS_SERVER << "/AFLOW/>" << endl;
      oss << "" << endl;
      oss << "AFLOW is free software: you can redistribute it and/or modify it under the terms of the" << endl;
      oss << "GNU General Public License as published by the Free Software Foundation, either version 3" << endl;
      oss << "of the License, or (at your option) any later version." << endl;
      oss << "" << endl;
      oss << "This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;" << endl;
      oss << "without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE." << endl;
      oss << "See the GNU General Public License for more details." << endl;
      oss << "" << endl;
      oss << "You should have received a copy of the GNU General Public License along with this program." << endl;
      oss << "If not, see <http://www.gnu.org/licenses/>." << endl;
      oss << "" << endl;
      oss << "AFLOW V" << string(AFLOW_VERSION) << " [" << XHOST.hostname << "] [" << XHOST.MachineType << "] ["<< XHOST.CPU_Cores << "] [" << XHOST.Find_Parameters << "]" << endl;
      //     oss << endl;
      return oss.str();
    }
    if(type=="INTRODUCTION") {
      stringstream oss;
      oss << "MMMMM  AFLOW VERSION " << string(AFLOW_VERSION) << " Automatic-FLOW for materials discovery - [" << TODAY << "] -" << aflow_get_time_string() << endl;
      //    oss << "MMMMM  AFLOW VERSION " << aurostd::PaddedPOST(string(AFLOW_VERSION),5) <<" - BUILT ["<<TODAY<<"] - (C) " << XHOST.Copyright_Years << "    " << endl;
      oss << "MMMMM  AFLOW.org consortium - High-Throughput ab-initio Computing Project - (C) " << XHOST.Copyright_Years << "    " << endl;
      oss << "MMMMM  ";// << endl;
      oss << "MMMMM  AFLOW is free software: you can redistribute it and/or modify it under the terms of the" << endl;
      oss << "MMMMM  GNU General Public License as published by the Free Software Foundation, either version 3" << endl;
      oss << "MMMMM  of the License, or (at your option) any later version." << endl;
      oss << "MMMMM  " << endl;
      oss << "MMMMM  This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;" << endl;
      oss << "MMMMM  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE." << endl;
      oss << "MMMMM  See the GNU General Public License for more details." << endl;
      oss << "MMMMM  " << endl;
      oss << "MMMMM  You should have received a copy of the GNU General Public License along with this program." << endl;
      oss << "MMMMM  If not, see <http://www.gnu.org/licenses/>." << endl;
      oss << "MMMMM " << endl;
      return oss.str();
    }
    if(type=="BANNER_NORMAL") {
      oss << "****************************************************************************************************" << endl;
      oss << "MMMMM  AFLOW V" << string(AFLOW_VERSION) << " Automatic-FLOW [" << TODAY << "] - " << endl; // << aflow_get_time_string() << endl; // CO
      oss << "****************************************************************************************************" << endl;
      return oss.str();
    }
    if(type=="BANNER_BIG") {
      oss << "****************************************************************************************************" << endl;
      oss << "*                                                                                                  *" << endl;
      oss << "*                          aflow - Automatic-FLOW for materials discovery                          *" << endl;
      oss << "*                aflow.org consortium - High-Throughput ab-initio Computing Project                *" << endl;
      oss << "*" << aurostd::PaddedCENTER(string("version "+string(AFLOW_VERSION)+" - g++/gcc "+aurostd::utype2string(__GNUC__)+"."+aurostd::utype2string(__GNUC_MINOR__)+"."+aurostd::utype2string(__GNUC_PATCHLEVEL__)+" - built ["+string(TODAY)+"] - (C) " +XHOST.Copyright_Years),100) << "*" << endl;
      oss << "*                                                                                                  *" << endl;
      oss << "****************************************************************************************************" << endl;
      oss << "*                                                                                                  *" << endl;
      oss << "*    AFLOW is free software: you can redistribute it and/or modify it under the terms of the       *" << endl;
      oss << "*    GNU General Public License as published by the Free Software Foundation, either version 3     *" << endl;
      oss << "*    of the License, or (at your option) any later version.                                        *" << endl;
      oss << "*                                                                                                  *" << endl;
      oss << "*    This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;     *" << endl;
      oss << "*    without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.     *" << endl;
      oss << "*    See the GNU General Public License for more details.                                          *" << endl;
      oss << "*                                                                                                  *" << endl;
      oss << "*    You should have received a copy of the GNU General Public License along with this program.    *" << endl;
      oss << "*    If not, see <http://www.gnu.org/licenses/>.                                                   *" << endl;
      oss << "*                                                                                                  *" << endl;
      oss << "****************************************************************************************************" << endl;
      oss << "*                                                                                                  *" << endl;
      oss << "*     Use of AFLOW software and repositories welcomes references to the following publications     *" << endl;
      oss << "*                                                                                                  *" << endl;
      oss << "*  Oses et al.      https://arxiv.org/abs/1806.06901 submitted                       (AFLOW-CHULL) *" << endl;
      oss << "*  Hicks et al.     https://arxiv.org/abs/1806.07864 submitted                       (ANRL proto2) *" << endl;
      oss << "*  Gossett et al.   Comp. Mat. Sci. (2018)           10.1016/j.commatsci.2018.03.075 (AFLOW-ML)    *" << endl;
      oss << "*  Hicks et al.     Acta Cryst. A74, 184-203 (2018)  10.1107/S2053273318003066       (AFLOW-SYM)   *" << endl;
      oss << "*  MBNardelli et al Comp. Mat. Sci. 143, 462 (2018)  10.1016/j.commatsci.2017.11.034 (PAOFLOW)     *" << endl;
      oss << "*  Rose et al.      Comp. Mat. Sci. 137, 362 (2017)  10.1016/j.commatsci.2017.04.036 (AFLUX lang)  *" << endl;
      oss << "*  Supka et al.     Comp. Mat. Sci. 136, 76 (2017)   10.1016/j.commatsci.2017.03.055 (AFLOWpi)     *" << endl;
      oss << "*  Plata et al.     npj Comput. Mater. 3, 45 (2017)  10.1038/s41524-017-0046-7       (AAPL kappa)  *" << endl;
      oss << "*  Toher et al.     Phys. Rev.Mater.1, 015401 (2017) 10.1103/PhysRevMaterials.1.015401 (AEL elast) *" << endl;
      oss << "*  Mehl et al.      Comp. Mat. Sci. 136, S1 (2017)   10.1016/j.commatsci.2017.01.017 (ANRL proto1) *" << endl;
      oss << "*  Calderon et al.  Comp. Mat. Sci. 108A, 233 (2015) 10.1016/j.commatsci.2015.07.019 (standard)    *" << endl;
      oss << "*  Toher et al.     Phys. Rev. B 90, 174107 (2014)   10.1103/PhysRevB.90.174107      (AGL Gibbs)   *" << endl;
      oss << "*  Taylor et al.    Comp. Mat. Sci. 93, 178 (2014)   10.1016/j.commatsci.2014.05.014 (REST-API)    *" << endl;
      oss << "*  Curtarolo et al. Comp. Mat. Sci. 58, 227 (2012)   10.1016/j.commatsci.2012.02.002 (AFLOW.org)   *" << endl;
      oss << "*  Curtarolo et al. Comp. Mat. Sci. 58, 218 (2012)   10.1016/j.commatsci.2012.02.005 (AFLOW C++)   *" << endl;
      oss << "*                                                                                                  *" << endl;
      oss << "****************************************************************************************************" << endl;
      oss << "*                                                                                                  *" << endl;
      oss << "*" << aurostd::PaddedCENTER(string("aflow/aflow.org - CONTRIBUTORS"),100) << "*" << endl;
      oss << "*  2000-2018 Stefano Curtarolo (aflow); 2002-2004 Dane Morgan (convasp); 2007-2011 Wahyu Setyawan  *" << endl;
      oss << "*  (--rsm --edos --kband --icsd*); 2008-2011 Roman Chepulskyy (--edos --kband  surfaces);          *" << endl;
      oss << "*  2008 Gus Hart (lattice reductions - prototypes); 2009-2011, Ohad Levy (prototypes);             *" << endl;
      oss << "*  2009-2010, Michal Jahnatek (APL); 2010-2013 Shidong Wang (cluster expansion); 2010-2013         *" << endl;
      oss << "*  Richard Taylor (surfaces, apennsy); 2010-2013 Junkai Xue (prototyper); 2010-2013 Kesong Yang    *" << endl;
      oss << "*  (findsym, frozsl, plotband/dos); 2013-2018 Cormac Toher (AGL Debye-Gruneisen, AEL elastic);     *" << endl;
      oss << "*  2013-2018 Frisco Rose (API, Aflux); 2013-2018 Pinku Nath (Quasi-harmonic approximation);        *" << endl;
      oss << "*  2013-2017 Jose J. Plata (AAPL, thermal cond.); 2014-2018 David Hicks (symmetry, structure       *" << endl;
      oss << "*  comparison, prototypes); 2014-2018 Corey Oses (bader, chull, APL, pocc);                        *" << endl;
      oss << "*                                                                                                  *" << endl;
      oss << "****************************************************************************************************" << endl;
      oss << "*" << aurostd::PaddedCENTER(string("version "+string(AFLOW_VERSION)+" - g++/gcc "+aurostd::utype2string(__GNUC__)+"."+aurostd::utype2string(__GNUC_MINOR__)+"."+aurostd::utype2string(__GNUC_PATCHLEVEL__)+" - built ["+string(TODAY)+"] - (C) " +XHOST.Copyright_Years),100) << "*" << endl;
      oss << "****************************************************************************************************" << endl;
      return oss.str();
    }
    if(type=="BANNER_TINY") {
      // called 63 times
      oss << "AFLOW VERSION "<<AFLOW_VERSION<<":  [aflow.org consortium - 2003-2018] ";
      return oss.str();
    }
    if(type == "EXCEPTIONS") {
      oss << "List of AFLOW exceptions with error codes. See README_AFLOW_EXCEPTIONS.TXT for more information" << endl;
      oss << endl;
      oss << "----------------------------------------------------------------------------------------------------" << endl;
      oss << "Error Code      Error Type             Error                              Name of Constant          " << endl;
      oss << "----------------------------------------------------------------------------------------------------" << endl;
      oss << "        1       N/A                    Generic error                      _GENERIC_ERROR_           " << endl;
      oss << "        2                              Illegal error code                 _ILLEGAL_CODE_            " << endl;
      oss << "       10       Input Error            generic                            _INPUT_ERROR_             " << endl;
      oss << "       11                              unknown flag                       _INPUT_UNKNOWN_           " << endl;
      oss << "       12                              missing flag                       _INPUT_MISSING_           " << endl;
      oss << "       13                              input ambiguous                    _INPUT_AMBIGUOUS_         " << endl;
      oss << "       14                              illegal parameter                  _INPUT_ILLEGAL_           " << endl;
      oss << "       15                              number of parameters               _INPUT_NUMBER_            " << endl;
      oss << "       20       File Error             generic                            _FILE_ERROR_              " << endl;
      oss << "       21                              file not found                     _FILE_NOT_FOUND_          " << endl;
      oss << "       22                              wrong format                       _FILE_WRONG_FORMAT_       " << endl;
      oss << "       23                              file corrupt                       _FILE_CORRUPT_            " << endl;
      oss << "       30       Value Error            generic                            _VALUE_ERROR_             " << endl;
      oss << "       31                              illegal value                      _VALUE_ILLEGAL_           " << endl;
      oss << "       32                              out of range                       _VALUE_RANGE_             " << endl;
      oss << "       40       Index Error            generic                            _INDEX_ERROR_             " << endl;
      oss << "       41                              illegal value                      _INDEX_ILLEGAL_           " << endl;
      oss << "       42                              out of bounds                      _INDEX_BOUNDS_            " << endl;
      oss << "       43                              mismatch                           _INDEX_MISMATCH_          " << endl;
      oss << "       50       Runtime Error          generic                            _RUNTIME_ERROR_           " << endl;
      oss << "       51                              not initialized                    _RUNTIME_INIT_            " << endl;
      oss << "       60       Allocation Error       generic                            _ALLOC_ERROR_             " << endl;
      oss << "       61                              could not allocate memory          _ALLOC_ALLOCATE_          " << endl;
      oss << "       62                              insufficient memory                _ALLOC_INSUFFICIENT_      " << endl;
      oss << "----------------------------------------------------------------------------------------------------" << endl;
      return oss.str();
    }
    cerr << "aflow::Banner type=" << type << " not found..." << endl;
    oss << "aflow::Banner type=" << type << " not found..." << endl;
    //exit(0);
    return 0; // CO 180419
    return oss.str();
  }
} // namespace aflow

// ***************************************************************************
// *                                                                         *
// *         aflow - Automatic FLOW for materials discovery project          *
// *                                                                         *
// ***************************************************************************

// Update Mon Mar  7 14:05:40 EST 2011 (by wahyu):
// ICSD prototypes of compounds used in aflow are generated using the following
// command:
// xzcat /common/NIST/$1ary.icsd.xz | pflow --icsd_nobrokenbasis |
// pflow --icsd_nopartialocc | pflow --icsd2proto > README_LIBRARY_ICSD$1.TXT
