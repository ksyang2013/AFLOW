// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *                                                                         *
// ***************************************************************************
// A1,A2,A3 for all potpaw_GGA
/*
   aflow --noldau --potential=potpaw_GGA --potential_complete --aflow_proto=A1,A2,A3:Ac,Au,Bi_d,Cd,Co,Fe_sv,Ge,H.75,I,La,Mn,Na_sv,Ni_pv,O_s,Pd,Pu,Rh,Sc,Sm_3,Tc,Ti_sv,V,Yb_2,Ac_s,B,Br,Ce,Cr,Dy_3,F_h,Ge_d,He,In,La_s,Mn_pv,Nb_pv,Np,Os_pv,Pd_pv,Pu_s,Rh_pv,Sc_sv,Sn,Tc_pv,Tl,V_pv,Y_sv,Ag,Ba_sv,B_s,Ce_3,Cr_pv,Er_3,F_s,Ge_h,Hf,In_d,Li,Mo,Nb_sv,Np_s,P,P_h,Rb_pv,Ru,Se,Sn_d,Te,Tl_d,V_sv,Zn,Al,Be,C,Ce_s,C_s,Eu_2,Ga,H,Hf_pv,Ir,Li_sv,Mo_pv,Nd_3,N_s,Pa,Pm_3,Rb_sv,Ru_pv,S_h,Sr_sv,Th,Tm,W,Zr,Al_h,Be_sv,Ca,C_h,Cs_sv,F,Ga_d,H1.25,Hg,K_pv,Lu_3,N,Ne,O,Pa_s,Re,Ru_sv,Si,Ta,Th_s,Tm_3,W_pv,Zr_sv,Ar,B_h,Ca_pv,Cl,Cu,Fe,Ga_h,H1.5,H_h,Kr,Mg,Na,N_h,O_h,Pb,Pr_3,Si_h,Ta_pv,Ti,U,Xe,As,Bi,Ca_sv,Cl_h,Cu_pv,Fe_pv,Gd_3,H.5,Ho_3,K_sv,Mg_pv,Na_pv,Ni,Os,Pb_d,Pt,Re_pv,Sb,Sm_2,Tb_3,Ti_pv,U_s,Yb


   aflow --noldau --potential=potpaw_LDA --potential_complete --aflow_proto=A1,A2,A3:Ac,As,B_h,Ca_pv,Cl,Cs_sv,Fe_pv,Ge,H.5,H_h,Kr,Mg,Mo_sv,Nb_sv,Np_s,Os_pv,Pd,Pu_s,Rh,Sc_sv,Sn_d,Te,Tl,V_sv,Zr_sv,Ac_s,Au,Bi,Ca_sv,Cl_h,Cu,F_h,Ge_d,H.75,I,K_sv,Mg_pv,N,Ne,N_s,P,Pd_pv,Rb_pv,Rh_pv,Se,Sr_sv,Th,Tl_d,W,Ag,B,Bi_d,Cd,Co,Cu_pv,F_s,Ge_h,He,In,La,Mn,Na,N_h,O,Pa,P_h,Rb_sv,Ru,S_h,Ta,Th_s,U,W_pv,Al,Ba_sv,Br,Ce,Cr,Ga,H,Hf,In_d,La_s,Mn_pv,Na_pv,Ni,O_h,Pa_s,Re,Ru_pv,Si,Ta_pv,Ti,U_s,Xe,Al_h,Be,B_s,Ce_s,Cr_pv,F,Ga_d,H1.25,Hf_pv,Ir,Li,Mo,Na_sv,Ni_pv,Os,Pb,Pt,S,Si_h,Tc,Ti_pv,V,Y_sv,Ar,Be_sv,C,C_h,C_s,Fe,Ga_h,H1.5,Hg,K_pv,Li_sv,Mo_pv,Nb_pv,Np,O_s,Pb_d,Pu,Re_pv,Sb,Sn,Tc_pv,Ti_sv,V_pv,Zn


   aflow --noldau --potential=potpaw_PBE --potential_complete --aflow_proto=A1,A2,A3:Ac,Au,Bi_d,Ce,Cr_pv,Dy_3,Fe_pv,Gd_3,H.5,Ho_3,K_sv,Mg,Na,Ne,O,Pa_s,Pm_3,Rb_pv,Ru,Si,Ta,Th_s,Tm_3,W_pv,Zr_sv,Ac_s,B,Br,Ce_3,C_s,Er_2,F_h,Ge,H.75,I,La,Mg_pv,Na_pv,N_h,O_h,Pb,Rb_sv,Ru_pv,Si_h,Ta_pv,Ti,U,Xe,Ag,Ba_sv,B_s,C_h,Cs_sv,Er_3,F_s,Ge_d,He,In,La_s,Mn,Na_sv,Ni,Os,Pb_d,Pr,Re,S,Sm,Tb_3,Ti_pv,U_s,Yb,Al,Be,C,Cl,Cu,Eu,Ga,Ge_h,Hf,In_d,Li,Mn_pv,Nb_pv,Ni_pv,O_s,Pd,Pr_3,Sb,Sm_3,Tc,Ti_sv,V,Yb_2,Al_h,Be_sv,Ca_pv,Cl_h,Cu_pv,Eu_2,Ga_d,H,Hf_pv,Ir,Li_sv,Mo,Nb_sv,Np,Os_pv,Pd_pv,Pt,Re_pv,Sc_sv,Sn,Tc_pv,Tl,V_pv,Y_sv,Ar,B_h,Ca_sv,Co,F,Ga_h,H1.25,Hg,K_pv,Lu,Mo_pv,Nd,Np_s,P,P_h,Pu,Rh,Se,Sn_d,Te,Tl_d,V_sv,Zn,As,Bi,Cd,Cr,Fe,Gd,H1.5,H_h,Kr,Lu_3,N,Nd_3,N_s,Pa,Pm,Pu_s,Rh_pv,S_h,Sr_sv,Th,Tm,W,Zr


   aflow --noldau --potential=pot_LDA --potential_complete --aflow_proto=A1,A2,A3:Ag,As,Be,B_s,Cd,Cs,Ge,I,K,Li,Mg_h,Mo,Ne,Os,Pb_d,Rb,Rh,Sc,Si_h,Sr_pv,Ti_pv,W,Zn,Al,Au,Bi,C,Cl,C_s,F_s,H1.25,Hf,In,K_pv,Li_h,Mg_pv,Mo_pv,Ni,O_s,Pd,Rb_pv,Ru,Sc_pv,Sn,Tc,Tl_d,Xe,Zr,Al_h,B,Bi_d,Ca,Co,Cs_pv,F,Ga,H_200eV,Hg,In_d,Kr,Li_pv,Mn,N,Nb,N_s,P,Rb_s,S,Se,Sn_d,Te,V,Y,Zr_pv,Ar,Ba_pv,Br,Ca_pv,Cr,Cu,Fe,Ga_d,H.75,H_soft,Ir,K_s,Mg,Mn_pv,Na,Nb_pv,O,Pb,Pt,Re,Sb,Si,Sr,Ti,V_pv,Y_pv
// removed Na_h Na_pv He

aflow --noldau --potential=pot_GGA --potential_complete --aflow_proto=A1,A2,A3:Ag,As,Be,C,Cl,C_s,Fe,Ge,H_soft,Ir,Li,Mg_h,Mo_pv,Na_pv,Ni,O_s,Pd,Rb_pv,S,Se,Sn_d,Tc,Tl_d,Xe,Zr,Al,Au,Bi,Ca,Co,Cs_pv,F_s,H_200eV,I,K,Li_h,Mg_pv,N,Nb,N_s,P,Re,Sb,Si,Sr,Te,V,Y,Zr_pv,Al_h,B,Br,Ca_pv,Cr,Cu,Ga,Hf,In,K_pv,Li_pv,Mn,Na,Nb_pv,O,Pb,Pt,Rh,Sc,Si_h,Sr_pv,Ti,V_pv,Y_pv,Ar,Ba_pv,B_s,Cd,Cs,F,Ga_d,Hg,In_d,Kr,Mg,Mo,Na_h,Ne,Os,Pb_d,Rb,Ru,Sc_pv,Sn,Ta,Ti_pv,W,Zn

// removed Na_h Na_pv


aflow --noldau --potential=potpaw_LDA --potential_complete --aflow_proto=A4:Si,Si_h,Ge,Ge_d,Ge_h,C,C_s,C_h
aflow --noldau --potential=potpaw_GGA --potential_complete --aflow_proto=A4:Si,Si_h,Ge,Ge_d,Ge_h,C,C_s,C_h
aflow --noldau --potential=potpaw_PBE --potential_complete --aflow_proto=A4:Si,Si_h,Ge,Ge_d,Ge_h,C,C_s,C_h
aflow --noldau --potential=pot_LDA --potential_complete --aflow_proto=A4:Si,Si_h,Ge,C,C_s
aflow --noldau --potential=pot_GGA --potential_complete --aflow_proto=A4:Si,Si_h,Ge,C,C_s


OHAD 600 series
aflow --neglect_nomix --aflow_proto=600.XXXXX:Ir:Nb,Ta,Mo &
aflow --neglect_nomix --aflow_proto=600.XXXXX:Mo,Nb:Os &
aflow --neglect_nomix --aflow_proto=600.XXXXX:Nb:Pd,Pt,Rh &
aflow --neglect_nomix --aflow_proto=600.XXXXX:Os,Pd,Pt,Rh:Ta &
aflow --neglect_nomix --aflow_proto=600.XXXXX:Ir,Os,Ru:W &
aflow --neglect_nomix --aflow_proto=600.XXXXX:Cr:Ru &


aflow --potential=potpaw_GGA --aflow_proto=ICSD_64654.AB,ICSD_44587.AB,ICSD_602748.AB,ICSD_615455.AB,ICSD_615459.AB,ICSD_615461.AB,ICSD_615462.AB,ICSD_615466.AB,ICSD_615469.AB,ICSD_615471.AB,ICSD_615473.AB,ICSD_615474.AB,ICSD_615476.AB,ICSD_615478.AB,ICSD_615441.AB,ICSD_615444.AB,ICSD_615446.AB,ICSD_615449.AB,ICSD_615453.AB,ICSD_248240.AB,ICSD_248241.AB,ICSD_248242.AB,ICSD_248243.AB:B_h:Sm_3 &
aflow --potential=potpaw_GGA --aflow_proto=ALL:B_h:Sm_3 &

ICSD_64654.AB,ICSD_44587.AB,ICSD_602748.AB,ICSD_615455.AB,ICSD_615459.AB,ICSD_615461.AB,ICSD_615462.AB,ICSD_615466.AB,ICSD_615469.AB,ICSD_615471.AB,ICSD_615473.AB,ICSD_615474.AB,ICSD_615476.AB,ICSD_615478.AB,ICSD_615441.AB,ICSD_615444.AB,ICSD_615446.AB,ICSD_615449.AB,ICSD_615453.AB,ICSD_248240.AB,ICSD_248241.AB,ICSD_248242.AB,ICSD_248243.AB

*/

// this file contains the routines to prepare VASP input files
// Stefano Curtarolo - 2007-2011 Duke
// XBANDS put the DEBUG later (search BANDS)
// ./aflow --missing --aflow_proto bcc,fcc,hcp Ac,Ar,Be,Br,Cd,Cl_h,Cs_sv,Dy_3,F,Ga,Ge,H1.5,Hf_pv,In,K_sv,Lu,Mn_pv,Na_pv,Nd_3,Np,Os,Pa_s,P_h,Pr_3,Rb_sv,Rh_pv,Sc_sv,Sm,Ta,Te,Ti_sv,U,W,Y_sv,Ac_s,As,Be_sv,B_s,Ce,Co,Cu,Er_2,Fe,Ga_d,Ge_d,H.5,Hg,In_d,La,Lu_3,Mo,Na_sv,Ne,Np_s,O_s,Pb,Pm,Pt,Re,Ru,Se,Sm_3,Ta_pv,Th,Tl,U_s,W_pv,Zn,Ag,Au,B_h,C,Ce_3,Cr,Cu_pv,Er_3,Fe_pv,Ga_h,Ge_h,H.75,H_h,Ir,La_s,Mg,Mo_pv,Nb_pv,N_h,N_s,Os_pv,Pb_d,Pm_3,Pu,Ru_pv,S_h,Sn,Tb_3,Th_s,Tl_d,V,Xe,Zr,Al,B,Bi,Ca_pv,C_h,Cr_pv,Eu,F_h,Gd,H,He,Ho_3,K_pv,Li,Mg_pv,N,Nb_sv,Ni,O,P,Pd,Pu_s,Re_pv,S,Si,Sn_d,Tc,Ti,Tm,V_pv,Yb,Zr_sv,Al_h,Ba_sv,Bi_d,Ca_sv,Cl,C_s,Eu_2,F_s,Gd_3,H1.25,Hf,I,Kr,Li_sv,Mn,Na,Nd,Ni_pv,O_h,Pa,Pd_pv,Pr,Rb_pv,Rh,Sb,Si_h,Sr_sv,Tc_pv,Ti_pv,Tm_3,V_sv,Yb_2
// ./aflow  --aflow_proto A4,A9 C,C_h,C_s,S,S_h,Si,Si_h,Sn,Sn_d

#ifndef _AFLOW_AVASP_CPP
#define _AFLOW_AVASP_CPP

#include "aflow.h"
#include "aflow_pflow.h"

#define _ICSD_DIRBASE_ "ICSD"

#define _aflowinpad_ 60

string lattices[]={"","BCC","FCC","CUB","HEX","RHL","BCT","TET","ORC","ORCC","ORCF","ORCI","MCL","MCLC","TRI","XXX"};

#define _AVASP_DOUBLE2STRING_PRECISION_ 9

// Lances tests show that: (07/13/2010)
// * PREC=med and kppra of 6000 makes enthalpy errors of 2~3 meV or so. (Not good enough for a CE if the average enthalpy is ~10 meV)
// * With PREC=high and *equivalent set* of around kppra=2000 (normally we go 2000-6000 depending on lattice and other issues), the errors are about 10 times smaller.
// * PREC=high and a dense MP set of kppra of ~10000 is about as good as the equivalent set (with kppra ~2000)
// * PREC=high and a dense *equivalent set* of ~12000 kppra has an average error of a few hundredths of a meV/atom (which is overkill so we don't normally do that)
// The piece of information that we are running will be to compare PREC=med for both MP and equivalent sets.

#define CONVENTIONAL_FORCE "b46,b138,b144,b150,b156,b162,b168,b174,b180,b186,b192,b198,b204,b222,XXX"

#define LABEL_RAW_LIBRARY_PROTOTYPES "f1,f2,b1,b2,h1,h3,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,58,59,60,61,62,63,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,115,117,116,118,121,119,120,123,122,124,125,127,126,128,132,129,134,130,137,131,135,133,140,136,138,139,141,145,142,147,143,150,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,201,202,203,204,205,210,211,216,217,218,230,231,232,233,234,235,238,239,242,243,244,245,246,247,248,252,253,254,255,256,258,270,208,209,219,220,221,222,223,226,227,228,229,249,257,259,260,261,262,263,264,269,270,271,272,263,264,267,268,273,274,275,276,277,278,279,280,281,282,283,284,285,286,287,288,289,290,291,292,301,302,303,304,305,306,A7.A,A7.B,309,310,311,312,313,314,315,316,317,318,323,324,359,360,363,364,365,366,371,372,373,374,375,376,381,382,389,390,405,406,407,408,447,448,471,472,473,474,475,476,477,478,479,480,537,538,539,540,541,542,543,545,546,547,548,551,583,584,595,596,411,412,559,560,613,614,654.AB,654.BA,657.AB,657.BA,658.AB,658.BA,662.AB,662.BA,669.AB,670.AB,670.BA,678.AB,678.BA,681.AB,681.BA,689.AB,689.BA,692.AB,692.BA,693.AB,693.BA,694.AB,689.BA,695.AB,695.BA,696.AB,696.BA,701.AB,701.BA,707.AB,707.BA,708.AB,708.BA,715.AB,715.BA,716.AB,716.BA,720.AB,720.BA,721.AB,721.BA,722.AB,722.BA,724.AB,724.BA,726.AB,726.BA,727.AB,727.BA"

// 64,65,549,550 removed.. PAW bug

// mike aflow --aflow_proto f1,f2,b1,b2,h1,h3,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,58,59,60,61,62,63,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,115,117,116,118,121,119,120,123,122,124,125,127,126,128,132,129,134,130,137,131,135,133,140,136,138,139,141,145,142,147,143,150,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,201,202,203,204,205,210,211,216,217,218,230,231,232,233,234,235,238,239,242,243,244,245,246,247,248,252,253,254,255,256,258,270,208,209,219,220,221,222,223,226,227,228,229,249,257,259,260,261,262,263,264,269,270,271,272,263,264,267,268,273,274,275,276,277,278,279,280,281,282,283,284,285,286,287,288,289,290,291,292,301,302,303,304,305,306,A7.A,A7.B,309,310,311,312,313,314,315,316,317,318,323,324,359,360,363,364,365,366,371,372,373,374,375,376,381,382,389,390,405,406,407,408,447,448,471,472,473,474,475,476,477,478,479,480,537,538,539,540,541,542,543,545,546,547,548,551,583,584,595,596,411,412,559,560,613,614,654.AB,654.BA,657.AB,657.BA,658.AB,658.BA,662.AB,662.BA,669.AB,670.AB,670.BA,675.AB,675.BA,678.AB,678.BA,679.AB,680.AB,681.AB,682.AB,683.AB,684.AB,685.AB,686.AB,687.AB,688.AB,689.AB,690.AB,691.AB,692.AB,693.AB,694.AB,695.AB,696.AB,697.AB,698.AB,699.AB,700.AB,701.AB,702.AB,703.AB,704.AB,705.AB,706.AB,707.AB,708.AB,709.AB,710.AB,711.AB,712.AB,713.AB,714.AB,715.AB,716.AB,717.AB,718.AB,719.AB,720.AB,721.AB,722.AB,723.AB,724.AB,725.AB,726.AB,727.AB,728.A,679.BA,680.BA,681.BA,682.BA,683.BA,684.BA,685.BA,686.BA,687.BA,688.BA,689.BA,690.BA,691.BA,692.BA,693.BA,694.BA,695.BA,696.BA,697.BA,698.BA,699.BA,700.BA,701.BA,702.BA,703.BA,704.BA,705.BA,706.BA,707.BA,708.BA,709.BA,710.BA,711.BA,712.BA,713.BA,714.BA,715.BA,716.BA,717.BA,718.BA,719.BA,720.BA,721.BA,722.BA,723.BA,724.BA,725.BA,726.BA,727.BA,728.B  N W


// 325,326

//[CO 180705 - moved to aflow.h for broader access]#define SPECIE_RAW_LIB2U "Ag,Au,Cd,Co,Cr_pv,Cu_pv,Fe_pv,Hf_pv,Hg,Ir,La,Mn_pv,Mo_pv,Nb_sv,Ni_pv,Os_pv,Pd_pv,Pt,Re_pv,Rh_pv,Ru_pv,Sc_sv,Ta_pv,Tc_pv,Ti_sv,V_sv,W_pv,Y_sv,Zn,Zr_sv"

//[CO 180705 - moved to aflow.h for broader access]#define SPECIE_RAW_LIB2 "Ag,Al,As,Au,B_h,Ba_sv,Be_sv,Bi_d,Br,Ca_sv,Cd,Cl,Co,Cr_pv,Cu_pv,Fe_pv,Ga_h,Ge_h,Hf_pv,Hg,In_d,Ir,K_sv,La,Li_sv,Mg_pv,Mn_pv,Mo_pv,Na_pv,Nb_sv,Ni_pv,Os_pv,P,Pb_d,Pd_pv,Pt,Re_pv,Rh_pv,Ru_pv,Sb,Sc_sv,Se,Si,Sn,Sr_sv,Ta_pv,Tc_pv,Te,Ti_sv,Tl_d,V_sv,W_pv,Y_sv,Zn,Zr_sv"

//[CO 180705 - moved to aflow.h for broader access]//#define SPECIE_RAW_LIB3 "Ag,Au,Cd,Co,Cr_pv,Cu_pv,Fe_pv,Hf_pv,Hg,Ir,La,Mn_pv,Mo_pv,Nb_sv,Ni_pv,Os_pv,Pd_pv,Pt,Re_pv,Rh_pv,Ru_pv,Sc_sv,Ta_pv,Tc_pv,Ti_sv,V_sv,W_pv,Y_sv,Zn,Zr_sv"

//[CO 180705 - moved to aflow.h for broader access]#define SPECIE_RAW_LIB3 "Ag,Al,As,Au,B_h,Ba_sv,Be_sv,Bi_d,Br,Ca_sv,Cd,Cl,Co,Cr_pv,Cu_pv,Fe_pv,Ga_h,Ge_h,Hf_pv,Hg,In_d,Ir,K_sv,La,Li_sv,Mg_pv,Mn_pv,Mo_pv,Na_sv,Nb_sv,Ni_pv,Os_pv,P,Pb_d,Pd_pv,Pt,Re_pv,Rh_pv,Ru_pv,Sb,Sc_sv,Se,Si,Sn,Sr_sv,Ta_pv,Tc_pv,Te,Ti_sv,Tl_d,V_sv,W_pv,Y_sv,Zn,Zr_sv"
//[CO 180705 - moved to aflow.h for broader access]//#define SPECIE_RAW_LIB3 "Ag,Al,As,Au,B_h,Bi_d,Cd,Co,Cr_pv,Cu_pv,Fe_pv,Ga_h,Ge_h,Hf_pv,Hg,In_d,Ir,La,Mg_pv,Mn_pv,Mo_pv,Nb_sv,Ni_pv,Os_pv,P,Pb_d,Pd_pv,Pt,Re_pv,Rh_pv,Ru_pv,Sb,Sc_sv,Se,Si,Sn,Ta_pv,Te,Tc_pv,Ti_sv,V_sv,W_pv,Y_sv,Zn,Zr_sv"


// Ga_hNi_pv

string AVASP_Shortcuts_for_Binaries(string &label);
string AVASP_Shortcuts_for_Ternaries(string &label);

double VolumeSpecie(string specie) {
    return GetAtomVolume(KBIN::VASP_PseudoPotential_CleanName(specie));
}
double MassSpecie(string specie) {
    return GetAtomMass(KBIN::VASP_PseudoPotential_CleanName(specie));
}

bool AVASP_CheckHeuslersSanvito(_xvasp& xvasp, string& msg);
bool AVASP_ADD_LDAU(_xvasp &xvasp);
bool AVASP_REMOVE_LDAU(_xvasp &xvasp);

typedef struct {
    deque<_xvasp>  *dxvasp;     // CONTENT
    bool     flag_WRITE;  // CONTENT
    int      itbusy;      // FOR AVASP_AFLOWIN (ALL)
    int      ITHREAD;     // FOR MULTISH_TIMESHARING
    int      THREADS_MAX; // FOR MULTISH_TIMESHARING
} _threaded_AVASP_AFLOWIN_params;

pthread_mutex_t mutex_AVASP=PTHREAD_MUTEX_INITIALIZER;

void *_threaded_AVASP_AFLOWIN_CONCURRENT(void *ptr);
void *_threaded_AVASP_AFLOWIN_SEQUENTIAL(void *ptr);
bool AVASP_MakePrototype_AFLOWIN_LOOP(deque<_xvasp>& dxvasp,bool flag_WRITE);

// ***************************************************************************
void AVASP_Get_LDAU_Parameters(string _species,bool &LDAU,vector<string>& vLDAUspecies,vector<uint>& vLDAUtype,vector<int>& vLDAUL, vector<double>& vLDAUU, vector<double> &vLDAUJ) {
    // vLDAUtype = type: 1 Liechtenstein ; 2 Dudarev
    // vLDAUL = orbital s,p,d,f.. 0,1,2,3 (multiply times 2 to get LMAXMIX)
    int p=1,d=2,f=3;
    //  int s=0;
    // vLDAUU = Ueff
    // vLDAUJ = Jeff
    int DLTYPE=2; // force TYPE1 to be 2...  Dudarev`s
    //  bool is_O=FALSE;
    //  bool is_pure=FALSE;
    string species=KBIN::VASP_PseudoPotential_CleanName(_species);

    // f systems
    if(species=="La") { // REF: PRB 73, 115403 (2006)
        LDAU=TRUE;vLDAUspecies.push_back(species);vLDAUtype.push_back(DLTYPE);vLDAUL.push_back(f);vLDAUU.push_back(8.1);vLDAUJ.push_back(0.6);return;}  // #La-f orbital,U eff, J parameter
    if(species=="Ce") { // REF J. Chem. Phys. 123, 064701 (2005) PHYSICAL REVIEW B 75, 035115 (2007)
        LDAU=TRUE;vLDAUspecies.push_back(species);vLDAUtype.push_back(DLTYPE);vLDAUL.push_back(f);vLDAUU.push_back(7.0);vLDAUJ.push_back(0.7);return;}  // #Ce-f orbital,U eff, J parameter
    if(species=="Pr") { // REF [56] of setyawan ht-bands
        LDAU=TRUE;vLDAUspecies.push_back(species);vLDAUtype.push_back(DLTYPE);vLDAUL.push_back(f);vLDAUU.push_back(6.5);vLDAUJ.push_back(1.0);return;}  // #Pr-f orbital,U eff, J parameter
    if(species=="Gd") { // REF J. Phys.:Condens. Matter 9, 767 (1997) PHYSICAL REVIEW B 73, 094410 (2006)
        LDAU=TRUE;vLDAUspecies.push_back(species);vLDAUtype.push_back(DLTYPE);vLDAUL.push_back(f);vLDAUU.push_back(6.7);vLDAUJ.push_back(0.7);return;}  // #Gd-f orbital,U eff, J parameter
    if(species=="Nd") { // REF: wahyu fitting to XPS-BIS 4f levels: J. Phys. F Met. Phys. 11, 121 (1981)
        LDAU=TRUE;vLDAUspecies.push_back(species);vLDAUtype.push_back(DLTYPE);vLDAUL.push_back(f);vLDAUU.push_back(7.2);vLDAUJ.push_back(1.0);return;}  // #Nd-f orbital
    if(species=="Sm") { // REF: wahyu fitting to XPS-BIS 4f levels: J. Phys. F Met. Phys. 11, 121 (1981)
        LDAU=TRUE;vLDAUspecies.push_back(species);vLDAUtype.push_back(DLTYPE);vLDAUL.push_back(f);vLDAUU.push_back(7.4);vLDAUJ.push_back(1.0);return;}  // #Sm-f orbital
    if(species=="Eu") { // REF: wahyu fitting to XPS-BIS 4f levels: J. Phys. F Met. Phys. 11, 121 (1981)
        LDAU=TRUE;vLDAUspecies.push_back(species);vLDAUtype.push_back(DLTYPE);vLDAUL.push_back(f);vLDAUU.push_back(6.4);vLDAUJ.push_back(1.0);return;}  // #Eu-f orbital
    if(species=="Tm") { // REF: wahyu fitting to XPS-BIS 4f levels: J. Phys. F Met. Phys. 11, 121 (1981)
        LDAU=TRUE;vLDAUspecies.push_back(species);vLDAUtype.push_back(DLTYPE);vLDAUL.push_back(f);vLDAUU.push_back(7.0);vLDAUJ.push_back(1.0);return;}  // #Tm-f orbital
    if(species=="Dy") { // {10.1016/j.physc.2009.06.003},
        LDAU=TRUE;vLDAUspecies.push_back(species);vLDAUtype.push_back(2);vLDAUL.push_back(f);vLDAUU.push_back(5.6);vLDAUJ.push_back(0.0);return;}  // #Dy-f orbital
    if(species=="Yb") { // REF J. Phys.: Condens. Matter 18, 6769-6775 (2006)
        LDAU=TRUE;vLDAUspecies.push_back(species);vLDAUtype.push_back(DLTYPE);vLDAUL.push_back(f);vLDAUU.push_back(7.0);vLDAUJ.push_back(0.67);return;}  // #Yb-f orbital,U eff, J parameter
    if(species=="Yb2+") { // REF Antonov PRB 58 9752 (1998) YbX X=As Sb Bi
        LDAU=TRUE;vLDAUspecies.push_back(species);vLDAUtype.push_back(DLTYPE);vLDAUL.push_back(f);vLDAUU.push_back(5.3);vLDAUJ.push_back(0.0);return;}  // #Yb-f orbital,U eff for Yb2+
    if(species=="Yb3+") { // REF Antonov PRB 58 9752 (1998)
        LDAU=TRUE;vLDAUspecies.push_back(species);vLDAUtype.push_back(DLTYPE);vLDAUL.push_back(f);vLDAUU.push_back(8.8);vLDAUJ.push_back(0.0);return;}  // #Yb-f orbital,U eff for Yb3+
    if(species=="Lu") { // REF: PRB 73, 115403 (2006)
        LDAU=TRUE;vLDAUspecies.push_back(species);vLDAUtype.push_back(DLTYPE);vLDAUL.push_back(f);vLDAUU.push_back(4.8);vLDAUJ.push_back(0.95);return;}  // #Lu-f orbital,U eff, J parameter
    if(species=="U") { // http://cms.mpi.univie.ac.at/vasp-forum/forum_viewtopic.php?3.789
        LDAU=TRUE;vLDAUspecies.push_back(species);vLDAUtype.push_back(2);vLDAUL.push_back(f);vLDAUU.push_back(4.0);vLDAUJ.push_back(0.0);return;}  
    // #U-f orbital,U eff, J parameter
    //   LDAU=TRUE;vLDAUspecies.push_back(species);vLDAUtype.push_back(DLTYPE);vLDAUL.push_back(f);vLDAUU.push_back(4.5);vLDAUJ.push_back(0.51);return;
    //   // #U-f orbital,U eff  10.1103/PhysRevB.84.014116,
    if(species=="Th") { // {10.1103/PhysRevB.80.014108},
        LDAU=TRUE;vLDAUspecies.push_back(species);vLDAUtype.push_back(2);vLDAUL.push_back(f);vLDAUU.push_back(5.0);vLDAUJ.push_back(0.0);return;}  // #Th-f orbital

    // d systems
    if(species=="Sc") { // [14]
        LDAU=TRUE;vLDAUspecies.push_back(species);vLDAUtype.push_back(2);vLDAUL.push_back(d);vLDAUU.push_back(2.9);vLDAUJ.push_back(0.0);return;}  // #Sc-d orbital,U eff
    if(species=="Ti") { // [5] SrTiO3, PRL 98, 115503 (2007)
        LDAU=TRUE;vLDAUspecies.push_back(species);vLDAUtype.push_back(2);vLDAUL.push_back(d);vLDAUU.push_back(4.4);vLDAUJ.push_back(0.0);return;}  // #Ti-d orbital,U eff CEDER SAYS 0.0
    if(species=="V") { // REF Pickett PRB 58, 1201 (1998)  (Anisimov Zaanem Anderson 6.7)
        LDAU=TRUE;vLDAUspecies.push_back(species);vLDAUtype.push_back(2);vLDAUL.push_back(d);vLDAUU.push_back(2.7);vLDAUJ.push_back(0.0);return;}  // #V-d orbital,U eff CEDER SAYS 3.1
    if(species=="Cr") { // [4] VO, MnO, FeO, CoO, NiO, CuO, Cr2O3, PRB 73, 195107 (2006)
        LDAU=TRUE;vLDAUspecies.push_back(species);vLDAUtype.push_back(2);vLDAUL.push_back(d);vLDAUU.push_back(3.5);vLDAUJ.push_back(0.0);return;}  // #Cr-d orbital,U eff
    if(species=="Mn") { // [4] VO, MnO, FeO, CoO, NiO, CuO, Cr2O3, PRB 73, 195107 (2006)
        LDAU=TRUE;vLDAUspecies.push_back(species);vLDAUtype.push_back(2);vLDAUL.push_back(d);vLDAUU.push_back(4.0);vLDAUJ.push_back(0.0);return;}  // #Mn-d orbital,U eff  CEDER SAYS 3.9
    if(species=="Fe") { // [7] FeO, Fe2SiO4, PRB 71, 035105
        LDAU=TRUE;vLDAUspecies.push_back(species);vLDAUtype.push_back(2);vLDAUL.push_back(d);vLDAUU.push_back(4.6);vLDAUJ.push_back(0.0);return;}  // #Fe-d orbital,U eff
    if(species=="Co") { // [8] MnO, FeO, CoO, NiO, PRB 58, 1201 (1998)
        LDAU=TRUE;vLDAUspecies.push_back(species);vLDAUtype.push_back(2);vLDAUL.push_back(d);vLDAUU.push_back(5.0);vLDAUJ.push_back(0.0);return;}  // #Co-d orbital,U eff  CEDER SAYS 5.7
    if(species=="Ni") { // REF Pickett PRB 58, 1201 (1998)  (Anisimov Zaanem Anderson 8.0) (Cococcioni PRB71 035105 4.6)
        LDAU=TRUE;vLDAUspecies.push_back(species);vLDAUtype.push_back(2);vLDAUL.push_back(d);vLDAUU.push_back(5.1);vLDAUJ.push_back(0.0);return;}  // #Ni-d orbital,U eff CEDER SAYS 6.0
    if(species=="Cu") { // [4] VO, MnO, FeO, CoO, NiO, CuO, Cr2O3, PRB 73, 195107 (2006)
        LDAU=TRUE;vLDAUspecies.push_back(species);vLDAUtype.push_back(2);vLDAUL.push_back(d);vLDAUU.push_back(4.0);vLDAUJ.push_back(0.0);return;}  // #Cu-d orbital,U eff
    if(species=="Zn") { // REF Erhart, PRB 73 205203 (2006)
        LDAU=TRUE;vLDAUspecies.push_back(species);vLDAUtype.push_back(2);vLDAUL.push_back(d);vLDAUU.push_back(7.5);vLDAUJ.push_back(0.0);return;}  // #Zn-d orbital,U eff
    //  if(species=="Y") { // NONE
    //   LDAU=TRUE;vLDAUspecies.push_back(species);vLDAUtype.push_back(2);vLDAUL.push_back(d);vLDAUU.push_back(5.1);vLDAUJ.push_back(0.0);return;}  // #Y-d orbital,U eff
    //  if(species=="Zr") { // NONE
    //   LDAU=TRUE;vLDAUspecies.push_back(species);vLDAUtype.push_back(2);vLDAUL.push_back(d);vLDAUU.push_back(5.1);vLDAUJ.push_back(0.0);return;}  // #Zr-d orbital,U eff
    if(species=="Nb") { // [11] d-impurities in Rb, PRB 50, 16861 (1994)
        LDAU=TRUE;vLDAUspecies.push_back(species);vLDAUtype.push_back(2);vLDAUL.push_back(d);vLDAUU.push_back(2.1);vLDAUJ.push_back(0.0);return;}  // #Nb-d orbital,U eff CEDER SAYS 1.5
    if(species=="Mo") { // [11] d-impurities in Rb, PRB 50, 16861 (1994)
        LDAU=TRUE;vLDAUspecies.push_back(species);vLDAUtype.push_back(2);vLDAUL.push_back(d);vLDAUU.push_back(2.4);vLDAUJ.push_back(0.0);return;}  // #Mo-d orbital,U eff CEDER SAYS 3.5
    if(species=="Tc") { // [11] d-impurities in Rb, PRB 50, 16861 (1994)
        LDAU=TRUE;vLDAUspecies.push_back(species);vLDAUtype.push_back(2);vLDAUL.push_back(d);vLDAUU.push_back(2.7);vLDAUJ.push_back(0.0);return;}  // #Tc-d orbital,U eff
    if(species=="Ru") { // [11] d-impurities in Rb, PRB 50, 16861 (1994)
        LDAU=TRUE;vLDAUspecies.push_back(species);vLDAUtype.push_back(2);vLDAUL.push_back(d);vLDAUU.push_back(3.0);vLDAUJ.push_back(0.0);return;}  // #Ru-d orbital,U eff
    if(species=="Rh") { // [11] d-impurities in Rb, PRB 50, 16861 (1994)
        LDAU=TRUE;vLDAUspecies.push_back(species);vLDAUtype.push_back(2);vLDAUL.push_back(d);vLDAUU.push_back(3.3);vLDAUJ.push_back(0.0);return;}  // #Rh-d orbital,U eff
    if(species=="Pd") { // [11] d-impurities in Rb, PRB 50, 16861 (1994)
        LDAU=TRUE;vLDAUspecies.push_back(species);vLDAUtype.push_back(2);vLDAUL.push_back(d);vLDAUU.push_back(3.6);vLDAUJ.push_back(0.0);return;}  // #Pd-d orbital,U eff
    // if(species=="Ag") { // NONE
    // LDAU=TRUE;vLDAUspecies.push_back(species);vLDAUtype.push_back(2);vLDAUL.push_back(d);vLDAUU.push_back(5.1);vLDAUJ.push_back(0.0);return;}  // #Ag-d orbital,U eff
    if(species=="Ag") { // L. Wang, T. Maxisch, G. Ceder PRB 73 2006 19 WRONG
        // LDAU=TRUE;vLDAUspecies.push_back(species);vLDAUtype.push_back(2);vLDAUL.push_back(d);vLDAUU.push_back(4.0);vLDAUJ.push_back(0.0);return;  // #Ag-d orbital,U eff  // guessed 4
        // LDAU=TRUE;vLDAUspecies.push_back(species);vLDAUtype.push_back(2);vLDAUL.push_back(d);vLDAUU.push_back(1.5);vLDAUJ.push_back(0.0);return;  // #Ag-d orbital,U eff // L. Wang, T. Maxisch, G. Ceder PRB 73 2006 19 WRONG
        LDAU=TRUE;vLDAUspecies.push_back(species);vLDAUtype.push_back(2);vLDAUL.push_back(d);vLDAUU.push_back(5.8);vLDAUJ.push_back(0.0);return;}  // #Ag-d orbital,U eff // 10.1103/PhysRevB.83.035202
    if(species=="Cd") { // [2] ZnO, CdO, GaN, InN, PRB 74, 045202 (2006)
        LDAU=TRUE;vLDAUspecies.push_back(species);vLDAUtype.push_back(2);vLDAUL.push_back(d);vLDAUU.push_back(2.1);vLDAUJ.push_back(0.0);return;}  // #Cd-d orbital,U eff
    // if(species=="Hf")  // NONE
    // LDAU=TRUE;vLDAUspecies.push_back(species);vLDAUtype.push_back(2);vLDAUL.push_back(d);vLDAUU.push_back(5.1);vLDAUJ.push_back(0.0);return;  // #Hf-d orbital,U eff
    if(species=="Ta") { // [11] d-impurities in Rb, PRB 50, 16861 (1994)
        LDAU=TRUE;vLDAUspecies.push_back(species);vLDAUtype.push_back(2);vLDAUL.push_back(d);vLDAUU.push_back(2.0);vLDAUJ.push_back(0.0);return;}  // #Ta-d orbital,U eff
    if(species=="W") { // [11] d-impurities in Rb, PRB 50, 16861 (1994)
        LDAU=TRUE;vLDAUspecies.push_back(species);vLDAUtype.push_back(2);vLDAUL.push_back(d);vLDAUU.push_back(2.2);vLDAUJ.push_back(0.0);return;}  // #W-d orbital,U eff
    if(species=="Re") { // [11] d-impurities in Rb, PRB 50, 16861 (1994)
        LDAU=TRUE;vLDAUspecies.push_back(species);vLDAUtype.push_back(2);vLDAUL.push_back(d);vLDAUU.push_back(2.4);vLDAUJ.push_back(0.0);return;}  // #Re-d orbital,U eff
    if(species=="Os") { // [11] d-impurities in Rb, PRB 50, 16861 (1994)
        LDAU=TRUE;vLDAUspecies.push_back(species);vLDAUtype.push_back(2);vLDAUL.push_back(d);vLDAUU.push_back(2.6);vLDAUJ.push_back(0.0);return;}  // #Os-d orbital,U eff
    if(species=="Ir") { // [11] d-impurities in Rb, PRB 50, 16861 (1994)
        LDAU=TRUE;vLDAUspecies.push_back(species);vLDAUtype.push_back(2);vLDAUL.push_back(d);vLDAUU.push_back(2.8);vLDAUJ.push_back(0.0);return;}  // #Ir-d orbital,U eff
    if(species=="Pt") { // [11] d-impurities in Rb, PRB 50, 16861 (1994)
        LDAU=TRUE;vLDAUspecies.push_back(species);vLDAUtype.push_back(2);vLDAUL.push_back(d);vLDAUU.push_back(3.0);vLDAUJ.push_back(0.0);return;}  // #Pt-d orbital,U eff
    // if(species=="Au") { // NONE
    //  LDAU=TRUE;vLDAUspecies.push_back(species);vLDAUtype.push_back(2);vLDAUL.push_back(d);vLDAUU.push_back(5.1);vLDAUJ.push_back(0.0);return;}  // #Au-d orbital,U eff
    if(species=="Au") { // guessed
        LDAU=TRUE;vLDAUspecies.push_back(species);vLDAUtype.push_back(2);vLDAUL.push_back(d);vLDAUU.push_back(4.0);vLDAUJ.push_back(0.0);return;}  // #Au-d orbital,U eff
    // if(species=="Hg") { // NONE
    //   LDAU=TRUE;vLDAUspecies.push_back(species);vLDAUtype.push_back(2);vLDAUL.push_back(d);vLDAUU.push_back(5.1);vLDAUJ.push_back(0.0);return;}  // #Hg-d orbital,U eff
    // if(species=="Ce" && is_pure) { // REF Cococcioni PRB71 035105 Ce elemental 4.5  (4.4 from 4f to 6s  and 6.4 from 4f to 5d)
    //  LDAU=TRUE;vLDAUspecies.push_back(species);vLDAUtype.push_back(2);vLDAUL.push_back(d);vLDAUU.push_back(4.5);vLDAUJ.push_back(0.0);return;}  // #Ce-d orbital,U eff

    if(species=="Ga") { // [2] ZnO, CdO, GaN, InN, PRB 74, 045202 (2006)
        LDAU=TRUE;vLDAUspecies.push_back(species);vLDAUtype.push_back(2);vLDAUL.push_back(d);vLDAUU.push_back(3.9);vLDAUJ.push_back(0.0);return;}  // #Ga-d orbital,U eff
    if(species=="In") { // [2] ZnO, CdO, GaN, InN, PRB 74, 045202 (2006)
        LDAU=TRUE;vLDAUspecies.push_back(species);vLDAUtype.push_back(2);vLDAUL.push_back(d);vLDAUU.push_back(1.9);vLDAUJ.push_back(0.0);return;}  // #In-d orbital,U eff

    if(species=="Bi") { // ceder
        LDAU=TRUE;vLDAUspecies.push_back(species);vLDAUtype.push_back(2);vLDAUL.push_back(p);vLDAUU.push_back(0.0);vLDAUJ.push_back(0.0);return;}  // #p
    if(species=="Pb") { // ceder
        LDAU=TRUE;vLDAUspecies.push_back(species);vLDAUtype.push_back(2);vLDAUL.push_back(p);vLDAUU.push_back(0.0);vLDAUJ.push_back(0.0);return;}  // #p
    if(species=="Sn") { // [Sn d-orbital U = 3.5eV 10.1103/PhysRevLett.101.055502,
        LDAU=TRUE;vLDAUspecies.push_back(species);vLDAUtype.push_back(2);vLDAUL.push_back(d);vLDAUU.push_back(3.5);vLDAUJ.push_back(0.0);return;}  // #p

    // ELSE
    if(species=="Np") {
        cerr << "AVASP_Get_LDAU2_Parameters: LDAU for " << species << " is not implemented yet, exiting (WAHYU): " << species  << endl;
        exit(0);
        return;
    }

    // LDAU=FALSE; // dont modify
    vLDAUspecies.push_back(species);vLDAUtype.push_back(0);vLDAUL.push_back(-1);vLDAUU.push_back(0.0);vLDAUJ.push_back(0.0);  // #nothing NO LDAU
    // the LDAUL=0 will be mapped in -1 as in http://cms.mpi.univie.ac.at/vasp/vasp/node161.html
    return;
}

// ------------------------------------------------------------------------------------------------------
// Sc     |  Ti    | V      | Cr     | Mn     | Fe     | Co     | Ni     | Cu     | Zn    | Ga    |
// 2.9[9] |  4.4[3]| 2.7[12]| 3.5[19]| 4.0[19]| 4.6[2] | 5.0[12]| 5.1[12]| 4.0[19]| 7.5[4]| 3.9[7]|
// ------------------------------------------------------------------------------------------------------
// Y      |  Zr    | Nb     | Mo     | Tc     | Ru     | Rh     | Pd     | Ag     | Cd    | In    | Sn
// -      |  -     | 2.1[17]| 2.4[17]| 2.7[17]| 3.0[17]| 3.3[4] | 3.6[4] | 5.8[18]| 2.1[4]| 1.9[4]| 3.5[16]
// ------------------------------------------------------------------------------------------------------
//        |  Hf    | Ta     | W      | Re     | Os     | Ir     | Pt     | Au     | Hg    | Tl    | Pb
//        |  -     | 2.0[17]| 2.2[17]| 2.4[17]| 2.6[17]| 2.8[17]| 3.0[17]| -      | -     | -     | -
// ------------------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------------------
// La     |  Ce    | Pr     | Nd     | Sm     | Eu     | Gd     | Dy     | Tm     | Yb    | Lu
// 7.5[20]| 6.3[10]| 5.5[5] | 6.2[14]| 6.4[14]| 5.4[14]| 6.0[6] | 5.6[11]| 6.0[1] | 6.3[8]| 3.9[21]
// ------------------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------------------
// Ac     | Th     | Pa     | U
// -      | 5.0[15]| -      | 4.0[13]
// ------------------------------------------------------------------------------------------------------
//
// [1]  PRB 63, 205112 (2001)
// [2]  PRB 71, 035105 (2005)
// [3]  PRL 98, 115503 (2007)
// [4]  PRB 73, 205203 (2006)
// [5]  JMMM 226, 83 (2010)
// [6]  JPCS 56, 1521 (1995)
// [7]  PRB 74, 045202 (2006)
// [8]  JPCM 18, 6769 (2006)
// [9]  PRB 82, 045108 (2011); 10.1103/PhysRevB.82.045108; http://link.aps.org/doi/10.1103/PhysRevB.82.045108
// [10] JCP 123, 064701 (2005)
// [11] Physica C: Superconductivity 469, 1892 (2009)
// [12] PRB 58, 1201 (1998)
// [13] PRB 84, 014116 (2011); 10.1103/PhysRevB.84.014116; http://link.aps.org/doi/10.1103/PhysRevB.84.014116
// [14] ACSCombSci 13, 382 (2011)
// [15] PRB 80, 014108 (2009); 10.1103/PhysRevB.80.014108; http://link.aps.org/doi/10.1103/PhysRevB.80.014108
// [16] PRL 101, 055502 (2008); 10.1103/PhysRevLett.101.055502; http://link.aps.org/doi/10.1103/PhysRevLett.101.055502
// [17] PRB 50, 16861 (1994)
// [18] PRB 83, 035202 (2011); 10.1103/PhysRevB.83.035202; http://link.aps.org/doi/10.1103/PhysRevB.83.035202
// [19] PRB 73, 195107 (2006)
// [20] PRB 73, 115403 (2006)


//actual usage in generating potential
// ***************************************************************************
string AVASP_Get_PseudoPotential_PAW_PBE(string species) {
    string error;
    bool ALLOW_ACTINIDIES=TRUE; //FALSE;
    if(ALLOW_ACTINIDIES==FALSE) {
        if(species=="Ra" || species=="Ac" || species=="Th" || species=="Pa" || species=="U" ||
                species=="Np" || species=="Pu" || species=="Am" || species=="Cm" || species=="Bk" || species=="Cf" ||
                species=="Es" || species=="Fm" || species=="Md" || species=="No" || species=="Lw") {
            error="WARNING AVASP_Get_PseudoPotential_PAW_PBE: not producing Actinides";
            aurostd::execute("echo "+error+" >> "+DEFAULT_AFLOW_ERVASP_OUT);
            cerr << error << endl;
            exit(0);
        }
        if(species=="D" || species=="T") {
            error="WARNING AVASP_Get_PseudoPotential_PAW_PBE: not producing heavy hydrogen";
            aurostd::execute("echo "+error+" >> "+DEFAULT_AFLOW_ERVASP_OUT);
            cerr << error << endl;
            exit(0);
        }
    }

    if(species=="Ac") return "Ac"; // UNTESTED STEFANO of Ac,Ac_s take the hardest one [Wed Nov 23 EST 2011]     [PAW_PBE]
    if(species=="Ag") return "Ag"; // unique choice...     [PAW_PBE]
    if(species=="Al") return "Al"; // STEFANO Al_h has ghost states     [PAW_PBE]
    if(species=="Ar") return "Ar"; // unique choice...     [PAW_PBE]
    if(species=="As") return "As"; // unique choice...     [PAW_PBE]
    if(species=="Au") return "Au"; // unique choice...     [PAW_PBE]
    if(species=="B")  return "B_h"; // STEFANO     [PAW_PBE]
    //if(species=="B")  return "B"; //vasp recommends     [PAW_PBE]
    if(species=="Ba") return "Ba_sv"; // unique choice...     [PAW_PBE]
    if(species=="Be") return "Be_sv"; // STEFANO     [PAW_PBE]
    //if(species=="Bi") return "Bi_d"; //STEFANO,WAHYU     [PAW_PBE]
    if(species=="Bi") return "Bi"; //HY
    if(species=="Br") return "Br"; // unique choice...     [PAW_PBE]
    if(species=="C")  return "C"; // STEFANO, vasp recommends     [PAW_PBE]
    if(species=="Ca") return "Ca_sv"; // STEFANO     [PAW_PBE]
    if(species=="Cd") return "Cd"; // unique choice..     [PAW_PBE]
    if(species=="Ce") return "Ce_3"; //YANG; 2022-11-02
    //if(species=="Ce") return "Ce"; //WAHYU     [PAW_PBE]
    if(species=="Cl") return "Cl"; //WAHYU     [PAW_PBE]
    if(species=="Co") return "Co"; // unique choice..     [PAW_PBE]
    if(species=="Cr") return "Cr_pv"; // STEFANO     [PAW_PBE]
    if(species=="Cs") return "Cs_sv"; // unique choice..     [PAW_PBE]
    if(species=="Cu") return "Cu_pv"; // STEFANO     [PAW_PBE]
    if(species=="Dy") return "Dy_3"; // unique choice...     [PAW_PBE]
    if(species=="Er") return "Er_3"; // WAHYU, teragrid project     [PAW_PBE]
    if(species=="Eu") return "Eu_3"; // WAHYU, frozen f     [PAW_PBE] //YANG easy to get converged, 2024-02-11
    //if(species=="Eu") return "Eu"; // WAHYU, with f-states, teragrid project     [PAW_PBE]
    if(species=="F")  return "F"; // WAHYU     [PAW_PBE]
    if(species=="Fe") return "Fe_pv"; // STEFANO     [PAW_PBE]
    //if(species=="Ga") return "Ga_h"; // STEFANO     [PAW_PBE]
    if(species=="Ga") return "Ga"; // HY
    if(species=="Ge") return "Ge_h"; // STEFANO     [PAW_PBE]
    // if(species=="Ge") return "Ge_d"; // WAHYU     [PAW_PBE]
    //if(species=="Gd") return "Gd"; // WAHYU // LDAU     [PAW_PBE]
    if(species=="Gd") return "Gd_3"; // WAHYU // LDAU     [PAW_PBE]  //YANG 2024-02-11
    if(species=="H")  return "H"; // WAHYU     [PAW_PBE]
    if(species=="He") return "He"; // unique choice...     [PAW_PBE]
    if(species=="Hf") return "Hf_pv"; // STEFANO     [PAW_PBE]
    if(species=="Hg") return "Hg"; // unique choice...     [PAW_PBE]
    if(species=="Ho") return "Ho_3"; // unique choice...     [PAW_PBE]
    if(species=="I")  return "I"; // unique choice...     [PAW_PBE]
    if(species=="In") return "In_d"; // STEFANO     [PAW_PBE]
    if(species=="Ir") return "Ir"; // unique choice...     [PAW_PBE]
    if(species=="K")  return "K_sv"; // STEFANO     [PAW_PBE]
    if(species=="Kr") return "Kr"; // unique choice...     [PAW_PBE]
    if(species=="La") return "La"; // STEFANO,WAHYU // LDAU     [PAW_PBE]
    if(species=="Li") return "Li_sv"; // STEFANO     [PAW_PBE]
    if(species=="Lu") return "Lu_3"; // WAHYU, frozen f     [PAW_PBE]  //Kesong fix convergence issue
    //if(species=="Lu") return "Lu"; //WAHYU, with f-states,  LDAU     [PAW_PBE]
    if(species=="Mg") return "Mg_pv"; // STEFANO     [PAW_PBE]
    if(species=="Mn") return "Mn_pv"; // STEFANO     [PAW_PBE]
    if(species=="Mo") return "Mo_pv"; // STEFANO     [PAW_PBE]
    //if(species=="Na") return "Na_sv"; // was pv "; // STEFANO,WAHYU     [PAW_PBE]
    if(species=="Na") return "Na"; // Kesong
    if(species=="N")  return "N"; // vasp recommends     [PAW_PBE]
    if(species=="Nb") return "Nb_sv"; // STEFANO     [PAW_PBE]
    if(species=="Nd") return "Nd_3"; // WAHYU, frozen f     [PAW_PBE] YANG, 2022-11-02, halide perovskit
    //if(species=="Nd") return "Nd"; // WAHYU, with f-states, teragrid project     [PAW_PBE]
    if(species=="Ne") return "Ne"; // unique choice...     [PAW_PBE]
    if(species=="Ni") return "Ni_pv"; // STEFANO     [PAW_PBE]
    if(species=="Np") return "Np_s"; // STEFANO NEVER USED JUST FOR COMPLETENESS     [PAW_PBE]
    if(species=="O")  return "O"; //WAHYU     [PAW_PBE]
    if(species=="Os") return "Os_pv"; // STEFANO     [PAW_PBE]
    if(species=="P")  return "P"; // WAHYU, teragrid project     [PAW_PBE]
    if(species=="Pa") return "Pa"; // UNTESTED STEFANO of Pa,Pa_s take the hardest one [Wed Nov 23 EST 2011]     [PAW_PBE]
    if(species=="Pb") return "Pb_d"; //STEFANO,WAHYU     [PAW_PBE]
    if(species=="Pd") return "Pd_pv"; // STEFANO     [PAW_PBE]
    //  if(species=="Pm") return "Pm_3"; // WAHYU, frozen f     [PAW_PBE]
    if(species=="Pm") return "Pm"; // WAHYU, with f-states, teragrid project     [PAW_PBE]
    if(species=="Pt") return "Pt"; // STEFANO     [PAW_PBE]
    if(species=="Pr") return "Pr_3"; // WAHYU, frozen f     [PAW_PBE] //YANG 2024-02-11
    //if(species=="Pr") return "Pr";  // WAHYU, with f-states, teragrid project     [PAW_PBE]
    if(species=="Pu") return "Pu_s"; // STEFANO NEVER USED JUST FOR COMPLETENESS     [PAW_PBE]
    if(species=="Rb") return "Rb_sv"; // vasp recommends     [PAW_PBE]
    if(species=="Re") return "Re_pv"; // STEFANO     [PAW_PBE]
    if(species=="Rh") return "Rh_pv"; // STEFANO     [PAW_PBE]
    if(species=="Ru") return "Ru_pv"; // STEFANO     [PAW_PBE]
    if(species=="S")  return "S"; // WAHYU, teragrid project     [PAW_PBE]
    if(species=="Sb") return "Sb"; // unique choice...     [PAW_PBE]
    if(species=="Sc") return "Sc_sv"; // unique choice...     [PAW_PBE]
    if(species=="Se") return "Se"; // unique choice...     [PAW_PBE]
    if(species=="Si") return "Si"; // was _h"; //STEFANO,WAHYU     [PAW_PBE]
    if(species=="Sm") return "Sm_3"; // WAHYU, frozen f     [PAW_PBE] //YANG 2024-02-11
    //if(species=="Sm") return "Sm"; // WAHYU, with f-states, teragrid project     [PAW_PBE]
    if(species=="Sn") return "Sn";//STEFANO     [PAW_PBE]
    if(species=="Sr") return "Sr_sv"; // unique choice...     [PAW_PBE]
    if(species=="Ta") return "Ta_pv"; // STEFANO     [PAW_PBE]
    if(species=="Tb") return "Tb_3"; //  unique choice...     [PAW_PBE]
    if(species=="Tc") return "Tc_pv"; // STEFANO     [PAW_PBE]
    if(species=="Te") return "Te"; // unique choice...     [PAW_PBE]
    if(species=="Th") return "Th_s"; // STEFANO NEVER USED JUST FOR COMPLETENESS     [PAW_PBE]
    if(species=="Ti") return "Ti_pv"; // STEFANO    //KESONG changes Ti_sv to Ti_pv, 
    if(species=="Tl") return "Tl_d"; // STEFANO     [PAW_PBE]
    if(species=="Tm") return "Tm_3"; // WAHYU, frozen f     [PAW_PBE] // YANG 2024-02-11
    //if(species=="Tm") return "Tm";  // WAHYU, with f-states, teragrid project     [PAW_PBE]
    if(species=="U") return "U"; // UNTESTED STEFANO of U,U_s take the hardest one [Wed Nov 23 EST 2011]     [PAW_PBE]
    if(species=="V")  return "V_sv"; // STEFANO     [PAW_PBE]
    if(species=="W")  return "W_pv"; // STEFANO     [PAW_PBE]
    if(species=="Xe") return "Xe"; // unique choice...     [PAW_PBE]
    if(species=="Y")  return "Y_sv"; // unique choice...     [PAW_PBE]
    if(species=="Yb") return "Yb_3"; // WAHYU, frozen f     [PAW_PBE]  //YANG 2024-02-17
    //if(species=="Yb") return "Yb_2"; // WAHYU, frozen f     [PAW_PBE]  //YANG 2024-02-1
    //if(species=="Yb") return "Yb";  // WAHYU, with f-states, teragrid project     [PAW_PBE]
    if(species=="Zn") return "Zn"; // unique choice...     [PAW_PBE]
    if(species=="Zr") return "Zr_sv"; // STEFANO     [PAW_PBE]

    if(species=="Po") {
        error="ERROR AVASP_Get_PseudoPotential_PAW_PBE: No pseudopotential available: "+species;
        aurostd::execute("echo "+error+" >> "+DEFAULT_AFLOW_ERVASP_OUT);
        cerr << error << endl;
        exit(0);
    }
    // If not found then UNKNOWN
    error="ERROR AVASP_Get_PseudoPotential_PAW_PBE: Potential Not found: "+species;
    aurostd::execute("echo "+error+" >> "+DEFAULT_AFLOW_ERVASP_OUT);
    cerr << error << endl;

    exit(0);
    return species; // you shall hope that there is something compatible
}

// ***************************************************************************
string AVASP_Get_PseudoPotential_PAW_GGA(string species) {
    if(species=="Gd") return "Gd_3"; // in potpaw_GGA only the Gd !!

    if(species=="Pr") return "Pr"; //"Pr_3";  // wahyu, teragrid project  // same as potpaw_PBE, compatibility
    if(species=="Lu") return "Lu"; //"Lu_3"; // wahyu // LDAU  // same as potpaw_PBE, compatibility
    if(species=="Sm") return "Sm_3"; //"Sm_3"; // wahyu, teragrid project // same as potpaw_PBE, compatibility
    if(species=="Nd") return "Nd"; //"Nd_3"; // wahyu, teragrid project // same as potpaw_PBE, compatibility
    if(species=="Eu") return "Eu"; //"Eu_2"; // wahyu, teragrid project // same as potpaw_PBE, compatibility
    if(species=="Pm") return "Pm"; //"Pm_3"; // wahyu, teragrid project // same as potpaw_PBE, compatibility
    return AVASP_Get_PseudoPotential_PAW_PBE(species);
}

// ***************************************************************************
string AVASP_Get_PseudoPotential_PAW_LDA(string species) {
    // not defined... let`s hope the PAWGGA is good naming
    return  AVASP_Get_PseudoPotential_PAW_GGA(species);
}

// ***************************************************************************
string AVASP_Get_PseudoPotential_GGA(string species) {

    if(species=="B")  return "B"; // vasp recommends    [GGA]
    if(species=="Ba") return "Ba_pv"; // unique choice...    [GGA]
    if(species=="Be") return "Be"; //  unique choice...    [GGA]
    if(species=="Bi") return "Bi"; //  unique choice...    [GGA]

    if(species=="C")  return "C"; // stefano    [GGA]
    if(species=="Ca") return "Ca_pv"; // stefano    [GGA]
    if(species=="Cl") return "Cl"; // unique choice..    [GGA]
    if(species=="Cr") return "Cr"; // unique choice..    [GGA]
    if(species=="Cs") return "Cs_pv"; // unique choice..    [GGA]
    if(species=="Cu") return "Cu"; // unique choice..    [GGA]

    if(species=="Fe") return "Fe"; // unique choice..    [GGA]
    if(species=="Ga") return "Ga_d"; // stefano    [GGA]
    if(species=="Ge") return "Ge"; // wahyu    [GGA]
    if(species=="H")  return "H_soft"; // wahyu    [GGA]
    if(species=="Hf") return "Hf"; // unique choice...    [GGA]
    if(species=="Hg") return "Hg"; // unique choice...    [GGA]

    if(species=="K")  return "K_pv"; // stefano    [GGA]
    if(species=="Kr") return "Kr"; // unique choice...    [GGA]

    if(species=="Li") return "Li_pv"; // stefano    [GGA]
    if(species=="Mg") return "Mg_pv"; // stefano    [GGA]
    if(species=="Mn") return "Mn"; // unique choice...    [GGA]

    if(species=="N")  return "N"; // vasp recommends    [GGA]
    //if(species=="Na") return "Na_pv"; // there is no _sv in potGGA    [GGA]
    if(species=="Na") return "Na"; // Kesong
    if(species=="Nb") return "Nb_pv"; // stefano    [GGA]
    if(species=="Ne") return "Ne"; // unique choice...    [GGA]
    if(species=="Ni") return "Ni"; // unique choice...    [GGA]
    if(species=="O")  return "O"; //wahyu    [GGA]
    if(species=="Os") return "Os"; // unique choice...    [GGA]

    if(species=="Pd") return "Pd"; // stefano    [GGA]
    if(species=="Pt") return "Pt"; // stefano    [GGA]
    if(species=="Rb") return "Rb_pv"; // vasp recommends    [GGA]
    if(species=="Re") return "Re"; // stefano    [GGA]
    if(species=="Rh") return "Rh"; // stefano    [GGA]
    if(species=="Ru") return "Ru"; // stefano    [GGA]
    if(species=="Sc") return "Sc_pv"; // unique choice...    [GGA]
    if(species=="Sm") return "Sm_3"; //   [GGA]
    if(species=="Sr") return "Sr_pv"; // unique choice...    [GGA]
    if(species=="Ta") return "Ta"; // stefano    [GGA]
    if(species=="Tc") return "Tc"; // stefano    [GGA]
    if(species=="Ti") return "Ti_pv"; // stefano    [GGA]
    if(species=="V")  return "V_pv"; // stefano    [GGA]
    if(species=="W")  return "W"; // stefano    [GGA]
    if(species=="Y")  return "Y_pv"; // unique choice...    [GGA]
    if(species=="Zr") return "Zr_pv"; // stefano    [GGA]

    return  AVASP_Get_PseudoPotential_PAW_PBE(species);
}

// ***************************************************************************
string AVASP_Get_PseudoPotential_LDA(string species) {
    // not defined... let`s hope the GGA is good naming
    return  AVASP_Get_PseudoPotential_GGA(species);
}

// ***************************************************************************
string AVASP_Get_PseudoPotential_PBE(string species) {
    // not defined... let`s hope the GGA is good naming
    return  AVASP_Get_PseudoPotential_GGA(species);
}


//Low priority  (first use above setting, if not found, then use this?) //KESONG 2019-08-03
// ***************************************************************************
string AVASP_Get_PseudoPotential_PAW_PBE_KIN(string species) {
    string error;
    bool ALLOW_ACTINIDIES=TRUE; //FALSE;
    if(ALLOW_ACTINIDIES==FALSE) {
        if(species=="Ra" || species=="Ac" || species=="Th" || species=="Pa" || species=="U" ||
                species=="Np" || species=="Pu" || species=="Am" || species=="Cm" || species=="Bk" || species=="Cf" ||
                species=="Es" || species=="Fm" || species=="Md" || species=="No" || species=="Lw") {
            error="WARNING AVASP_Get_PseudoPotential_PAW_PBE_KIN: not producing Actinides";
            aurostd::execute("echo "+error+" >> "+DEFAULT_AFLOW_ERVASP_OUT);
            cerr << error << endl;
            exit(0);
        }
        if(species=="D" || species=="T") {
            error="WARNING AVASP_Get_PseudoPotential_PAW_PBE_KIN: not producing heavy hydrogen";
            aurostd::execute("echo "+error+" >> "+DEFAULT_AFLOW_ERVASP_OUT);
            cerr << error << endl;
            exit(0);
        }
    }

    if(species=="Ac") return "Ac"; // unique choice    [PAW_PBE_KIN]
    if(species=="Ag") return "Ag_pv"; // injecting pv inside...           [PAW_PBE_KIN]
    if(species=="Al") return "Al"; // seems right           [PAW_PBE_KIN]
    if(species=="Ar") return "Ar"; // unique choice...           [PAW_PBE_KIN]
    if(species=="As") return "As_d"; // need test           [PAW_PBE_KIN]
    if(species=="Au") return "Au"; // unique choice...           [PAW_PBE_KIN]
    if(species=="B")  return "B_h"; //            [PAW_PBE_KIN]
    //if(species=="B")  return "B"; //vasp recommends           [PAW_PBE_KIN]
    if(species=="Ba") return "Ba_sv"; // unique choice...           [PAW_PBE_KIN]
    if(species=="Be") return "Be_sv"; // Rico suggests Be           [PAW_PBE_KIN]
    if(species=="Bi") return "Bi"; //           [PAW_PBE_KIN]
    if(species=="Br") return "Br"; // unique choice...           [PAW_PBE_KIN]
    if(species=="C")  return "C"; // vasp recommends           [PAW_PBE_KIN]
    if(species=="Ca") return "Ca_sv"; // stefano           [PAW_PBE_KIN]
    if(species=="Cd") return "Cd"; // unique choice..           [PAW_PBE_KIN]
    if(species=="Ce") return "Ce_3"; // right on the middle           [PAW_PBE_KIN]
    if(species=="Cl") return "Cl"; //           [PAW_PBE_KIN]
    if(species=="Co") return "Co_sv"; // let`s see           [PAW_PBE_KIN]
    if(species=="Cr") return "Cr_pv"; //            [PAW_PBE_KIN]
    if(species=="Cs") return "Cs_sv"; // unique choice..           [PAW_PBE_KIN]
    if(species=="Cu") return "Cu_pv"; //            [PAW_PBE_KIN]
    if(species=="Dy") return "Dy"; // a harder one           [PAW_PBE_KIN]
    if(species=="Er") return "Er"; // a harder one           [PAW_PBE_KIN]
    //  if(species=="Eu") return "Eu_3"; //  frozen f           [PAW_PBE_KIN]
    if(species=="Eu") return "Eu"; //  with f-states           [PAW_PBE_KIN]
    if(species=="F")  return "F"; //            [PAW_PBE_KIN]
    if(species=="Fe") return "Fe_pv"; //            [PAW_PBE_KIN]
    //if(species=="Ga") return "Ga_h"; // pick the _h instead of _d           [PAW_PBE_KIN]
    if(species=="Ga") return "Ga"; // HY
    if(species=="Ge") return "Ge_h"; // pick the _h instead of _d           [PAW_PBE_KIN]
    // if(species=="Ge") return "Ge_d"; //            [PAW_PBE_KIN]
    if(species=="Gd") return "Gd"; //  // LDAU           [PAW_PBE_KIN]
    if(species=="H")  return "H"; // commented since not ready to be used with SCAN; Rico April 12 2018           [PAW_PBE_KIN]
    if(species=="He") return "He"; // unique choice... // commented since not ready to be used with SCAN; Rico April 12 2018           [PAW_PBE_KIN]
    if(species=="Hf") return "Hf_sv"; //            [PAW_PBE_KIN]
    if(species=="Hg") return "Hg"; // unique choice...           [PAW_PBE_KIN]
    if(species=="Ho") return "Ho"; // finally a harder one           [PAW_PBE_KIN]
    if(species=="I")  return "I"; // unique choice...           [PAW_PBE_KIN]
    if(species=="In") return "In_d"; // harder           [PAW_PBE_KIN]
    if(species=="Ir") return "Ir"; // unique choice...           [PAW_PBE_KIN]
    if(species=="K")  return "K_sv"; //            [PAW_PBE_KIN]
    if(species=="Kr") return "Kr"; // unique choice...           [PAW_PBE_KIN]
    if(species=="La") return "La"; //  // LDAU           [PAW_PBE_KIN]
    if(species=="Li") return "Li"; // Rico suggests Li            [PAW_PBE_KIN]
    //  if(species=="Lu") return "Lu_3"; //  frozen f           [PAW_PBE_KIN]
    if(species=="Lu") return "Lu"; //, with f-states,  LDAU           [PAW_PBE_KIN]
    if(species=="Mg") return "Mg_pv"; //            [PAW_PBE_KIN]
    if(species=="Mn") return "Mn_pv"; // vasp suggested           [PAW_PBE_KIN]
    if(species=="Mo") return "Mo_sv"; // curious about _sv           [PAW_PBE_KIN]
    //if(species=="Na") return "Na_pv"; // Rico recommends _pv           [PAW_PBE_KIN]
    if(species=="Na") return "Na"; //  Kesong
    if(species=="N")  return "N"; // vasp recommends           [PAW_PBE_KIN]
    if(species=="Nb") return "Nb_sv"; //            [PAW_PBE_KIN]
    //  if(species=="Nd") return "Nd_3"; //  frozen f           [PAW_PBE_KIN]
    if(species=="Nd") return "Nd"; //  with f-states           [PAW_PBE_KIN]
    if(species=="Ne") return "Ne"; // unique choice...           [PAW_PBE_KIN]
    if(species=="Ni") return "Ni_pv"; //            [PAW_PBE_KIN]
    if(species=="Np") return "Np"; //  a harder one, NEVER USED JUST FOR COMPLETENESS           [PAW_PBE_KIN]
    if(species=="O")  return "O"; //           [PAW_PBE_KIN]
    if(species=="Os") return "Os_pv"; // stefano           [PAW_PBE_KIN]
    if(species=="P")  return "P"; //            [PAW_PBE_KIN]
    if(species=="Pa") return "Pa"; // UNTESTED  of Pa,Pa_s take the hardest one [Wed Nov 23 EST 2011]           [PAW_PBE_KIN]
    if(species=="Pb") return "Pb_d"; // the harder           [PAW_PBE_KIN]
    if(species=="Pd") return "Pd_pv"; //            [PAW_PBE_KIN]
    //  if(species=="Pm") return "Pm_3"; //  frozen f           [PAW_PBE_KIN]
    if(species=="Pm") return "Pm"; //  with f-states           [PAW_PBE_KIN]
    if(species=="Po") return "Po_d"; // the harder           [PAW_PBE_KIN]
    if(species=="Pt") return "Pt_pv"; // let`s try it           [PAW_PBE_KIN]
    //  if(species=="Pr") return "Pr_3"; //  frozen f           [PAW_PBE_KIN]
    if(species=="Pr") return "Pr";  //  with f-states           [PAW_PBE_KIN]
    if(species=="Pu") return "Pu"; //  NEVER USED JUST FOR COMPLETENESS           [PAW_PBE_KIN]
    if(species=="Rb") return "Rb_sv"; // vasp recommends           [PAW_PBE_KIN]
    if(species=="Re") return "Re_pv"; //            [PAW_PBE_KIN]
    if(species=="Rh") return "Rh_pv"; //            [PAW_PBE_KIN]
    if(species=="Ru") return "Ru_pv"; //            [PAW_PBE_KIN]
    if(species=="S")  return "S"; //            [PAW_PBE_KIN]
    if(species=="Sb") return "Sb"; // unique choice...           [PAW_PBE_KIN]
    if(species=="Sc") return "Sc_sv"; // pick the _sv           [PAW_PBE_KIN]
    if(species=="Se") return "Se"; // unique choice...           [PAW_PBE_KIN]
    if(species=="Si") return "Si"; // unique ?           [PAW_PBE_KIN]
    if(species=="Sm") return "Sm_3"; // WAHYU, frozen f //KESONG changes Sm to Sm_3 
    if(species=="Sn") return "Sn_d"; // the harder           [PAW_PBE_KIN]
    if(species=="Sr") return "Sr_sv"; // unique choice...           [PAW_PBE_KIN]
    if(species=="Ta") return "Ta_pv"; //            [PAW_PBE_KIN]
    if(species=="Tb") return "Tb"; //  a stronger one           [PAW_PBE_KIN]
    if(species=="Tc") return "Tc_pv"; //            [PAW_PBE_KIN]
    if(species=="Te") return "Te"; // unique choice...           [PAW_PBE_KIN]
    if(species=="Th") return "Th"; //  an harder one           [PAW_PBE_KIN]
    if(species=="Ti") return "Ti_pv"; // STEFANO    //KESONG changes Ti_sv to Ti_pv, 
    if(species=="Tl") return "Tl_d"; //            [PAW_PBE_KIN]
    //  if(species=="Tm") return "Tm_3"; //  frozen f           [PAW_PBE_KIN]
    if(species=="Tm") return "Tm";  //  with f-states           [PAW_PBE_KIN]
    if(species=="U") return "U"; // UNTESTED  of U,U_s take the hardest one [Wed Nov 23 EST 2011]           [PAW_PBE_KIN]
    if(species=="V")  return "V_sv"; // vasp recommended           [PAW_PBE_KIN]
    if(species=="W")  return "W_pv"; //            [PAW_PBE_KIN]
    if(species=="Xe") return "Xe"; // unique choice...           [PAW_PBE_KIN]
    if(species=="Y")  return "Y_sv"; // unique choice...           [PAW_PBE_KIN]
    if(species=="Yb") return "Yb";  //  with f-states           [PAW_PBE_KIN]
    if(species=="Zn") return "Zn"; // unique choice...           [PAW_PBE_KIN]
    if(species=="Zr") return "Zr_sv"; // unique           [PAW_PBE_KIN]

    // If not found then UNKNOWN
    error="ERROR AVASP_Get_PseudoPotential_PAW_PBE_KIN: Potential Not found: "+species;
    aurostd::execute("echo "+error+" >> "+DEFAULT_AFLOW_ERVASP_OUT);
    cerr << error << endl;

    exit(0);
    return species; // you shall hope that there is something compatible
}

// ***************************************************************************
string AVASP_Get_PseudoPotential_PAW_LDA_KIN(string species) {
    // not defined... let`s hope the PAW_PBE_KIN is good naming
    return  AVASP_Get_PseudoPotential_PAW_PBE_KIN(species);
}


// ***************************************************************************

bool AVASP_populateXVASP(const _aflags& aflags,const _kflags& kflags,const _vflags& vflags,_xvasp& xvasp){
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string soliloquy="AVASP_populateXVASP():";

    if(aflags.QUIET){}
    if(kflags.AFLOW_MODE_VASP){}
    if(vflags.KBIN_VASP_FORCE_OPTION_NBANDS_AUTO_isentry){}

    if(LDEBUG){cerr << soliloquy << " BEGIN" << endl;}

    xvasp.AVASP_potential=vflags.KBIN_VASP_FORCE_OPTION_AUTO_PSEUDOPOTENTIALS.xscheme;


    if(LDEBUG){cerr << soliloquy << " xvasp.AVASP_EXTRA_INCAR=" << endl << xvasp.AVASP_EXTRA_INCAR.str() << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.AVASP_KSCHEME=" << xvasp.AVASP_KSCHEME << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.AVASP_LDAU_PARAMETERS_STRING=" << xvasp.AVASP_LDAU_PARAMETERS_STRING << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.AVASP_STATIC_KSCHEME=" << xvasp.AVASP_STATIC_KSCHEME << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.AVASP_aflowin_only_if_missing=" << xvasp.AVASP_aflowin_only_if_missing << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.AVASP_alpha_fix=" << xvasp.AVASP_alpha_fix << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.AVASP_dirbase=" << xvasp.AVASP_dirbase << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.AVASP_directory_from_library_=" << xvasp.AVASP_directory_from_library_ << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.AVASP_flag_ABMIX_scheme=" << xvasp.AVASP_flag_ABMIX_scheme << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.AVASP_flag_ALGO_scheme=" << xvasp.AVASP_flag_ALGO_scheme << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.AVASP_flag_GENERATE=" << xvasp.AVASP_flag_GENERATE << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.AVASP_flag_MPI=" << xvasp.AVASP_flag_MPI << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.AVASP_flag_PRECISION_scheme=" << xvasp.AVASP_flag_PRECISION_scheme << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.AVASP_flag_RUN_RELAX=" << xvasp.AVASP_flag_RUN_RELAX << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.AVASP_flag_RUN_RELAX_STATIC=" << xvasp.AVASP_flag_RUN_RELAX_STATIC << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.AVASP_flag_RUN_RELAX_STATIC_BANDS=" << xvasp.AVASP_flag_RUN_RELAX_STATIC_BANDS << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.AVASP_flag_RUN_STATIC=" << xvasp.AVASP_flag_RUN_STATIC << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.AVASP_flag_RUN_STATIC_BANDS=" << xvasp.AVASP_flag_RUN_STATIC_BANDS << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.AVASP_flag_TYPE.xscheme=" << xvasp.AVASP_flag_TYPE.xscheme << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.AVASP_label=" << xvasp.AVASP_label << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.AVASP_libbase=" << xvasp.AVASP_libbase << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.AVASP_parameters=" << xvasp.AVASP_parameters << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.AVASP_path_BANDS=" << xvasp.AVASP_path_BANDS << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.AVASP_potential=" << xvasp.AVASP_potential << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.AVASP_prototype_from_library_=" << xvasp.AVASP_prototype_from_library_ << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.AVASP_prototype_mode=" << xvasp.AVASP_prototype_mode << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.AVASP_value_BANDS_GRID=" << xvasp.AVASP_value_BANDS_GRID << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.AVASP_value_KPPRA=" << xvasp.AVASP_value_KPPRA << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.AVASP_value_KPPRA_STATIC=" << xvasp.AVASP_value_KPPRA_STATIC << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.AVASP_value_NSW=" << xvasp.AVASP_value_NSW << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.AVASP_volume_in=" << xvasp.AVASP_volume_in << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.Directory=" << xvasp.Directory << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.POTCAR_TYPE_DATE_PRINT_flag=" << xvasp.POTCAR_TYPE_DATE_PRINT_flag << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.aopts.flag(\"AFLOWIN_FLAG::ABINIT\")=" << xvasp.aopts.flag("AFLOWIN_FLAG::ABINIT") << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.aopts.flag(\"AFLOWIN_FLAG::AIMS\")=" << xvasp.aopts.flag("AFLOWIN_FLAG::AIMS") << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.aopts.flag(\"AFLOWIN_FLAG::ALGORITHM\")=" << xvasp.aopts.flag("AFLOWIN_FLAG::ALGORITHM") << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.aopts.flag(\"AFLOWIN_FLAG::APL_SUPERCELL\")=" << xvasp.aopts.flag("AFLOWIN_FLAG::APL_SUPERCELL") << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.aopts.flag(\"AFLOWIN_FLAG::BANDS_GRID\")=" << xvasp.aopts.flag("AFLOWIN_FLAG::BANDS_GRID") << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.aopts.flag(\"AFLOWIN_FLAG::CONVERT_UNIT_CELL\")=" << xvasp.aopts.flag("AFLOWIN_FLAG::CONVERT_UNIT_CELL") << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.aopts.flag(\"AFLOWIN_FLAG::EDIFFG\")=" << xvasp.aopts.flag("AFLOWIN_FLAG::EDIFFG") << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.aopts.flag(\"AFLOWIN_FLAG::ENMAX_MULTIPLY\")=" << xvasp.aopts.flag("AFLOWIN_FLAG::ENMAX_MULTIPLY") << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.aopts.flag(\"AFLOWIN_FLAG::ENMAX_MULTYPLY\")=" << xvasp.aopts.flag("AFLOWIN_FLAG::ENMAX_MULTYPLY") << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.aopts.flag(\"AFLOWIN_FLAG::IVDW\")=" << xvasp.aopts.flag("AFLOWIN_FLAG::IVDW") << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.aopts.flag(\"AFLOWIN_FLAG::KPPRA\")=" << xvasp.aopts.flag("AFLOWIN_FLAG::KPPRA") << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.aopts.flag(\"AFLOWIN_FLAG::KPPRA_STATIC\")=" << xvasp.aopts.flag("AFLOWIN_FLAG::KPPRA_STATIC") << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.aopts.flag(\"AFLOWIN_FLAG::METAGGA\")=" << xvasp.aopts.flag("AFLOWIN_FLAG::METAGGA") << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.aopts.flag(\"AFLOWIN_FLAG::MODULE\")=" << xvasp.aopts.flag("AFLOWIN_FLAG::MODULE") << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.aopts.flag(\"AFLOWIN_FLAG::NBANDS\")=" << xvasp.aopts.flag("AFLOWIN_FLAG::NBANDS") << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.aopts.flag(\"AFLOWIN_FLAG::NO_VOLUME_ADJUSTMENT\")=" << xvasp.aopts.flag("AFLOWIN_FLAG::NO_VOLUME_ADJUSTMENT") << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.aopts.flag(\"AFLOWIN_FLAG::PARAMS\")=" << xvasp.aopts.flag("AFLOWIN_FLAG::PARAMS") << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.aopts.flag(\"AFLOWIN_FLAG::POTIM\")=" << xvasp.aopts.flag("AFLOWIN_FLAG::POTIM") << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.aopts.flag(\"AFLOWIN_FLAG::PRECISION\")=" << xvasp.aopts.flag("AFLOWIN_FLAG::PRECISION") << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.aopts.flag(\"AFLOWIN_FLAG::PSTRESS\")=" << xvasp.aopts.flag("AFLOWIN_FLAG::PSTRESS") << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.aopts.flag(\"AFLOWIN_FLAG::QE\")=" << xvasp.aopts.flag("AFLOWIN_FLAG::QE") << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.aopts.flag(\"AFLOWIN_FLAG::RELAX_MODE\")=" << xvasp.aopts.flag("AFLOWIN_FLAG::RELAX_MODE") << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.aopts.flag(\"AFLOWIN_FLAG::RELAX_TYPE\")=" << xvasp.aopts.flag("AFLOWIN_FLAG::RELAX_TYPE") << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.aopts.flag(\"AFLOWIN_FLAG::RELAX_TYPE\")=" << xvasp.aopts.flag("AFLOWIN_FLAG::RELAX_TYPE") << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.aopts.flag(\"AFLOWIN_FLAG::TYPE\")=" << xvasp.aopts.flag("AFLOWIN_FLAG::TYPE") << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.aopts.flag(\"AFLOWIN_FLAG::VASP\")=" << xvasp.aopts.flag("AFLOWIN_FLAG::VASP") << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.aopts.flag(\"AFLOWIN_FLAG::VOLUME_MULTIPLY_EQUAL\")=" << xvasp.aopts.flag("AFLOWIN_FLAG::VOLUME_MULTIPLY_EQUAL") << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.aopts.flag(\"AFLOWIN_FLAG::VOLUME_PLUS_EQUAL\")=" << xvasp.aopts.flag("AFLOWIN_FLAG::VOLUME_PLUS_EQUAL") << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.aopts.flag(\"AVASP_flag_CONVERT_UNIT_CELL_CARTESIAN\")=" << xvasp.aopts.flag("AVASP_flag_CONVERT_UNIT_CELL_CARTESIAN") << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.aopts.flag(\"AVASP_flag_CONVERT_UNIT_CELL_COMPACT\")=" << xvasp.aopts.flag("AVASP_flag_CONVERT_UNIT_CELL_COMPACT") << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.aopts.flag(\"AVASP_flag_CONVERT_UNIT_CELL_FRACTIONAL\")=" << xvasp.aopts.flag("AVASP_flag_CONVERT_UNIT_CELL_FRACTIONAL") << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.aopts.flag(\"AVASP_flag_CONVERT_UNIT_CELL_INCELL\")=" << xvasp.aopts.flag("AVASP_flag_CONVERT_UNIT_CELL_INCELL") << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.aopts.flag(\"AVASP_flag_CONVERT_UNIT_CELL_INCOMPACT\")=" << xvasp.aopts.flag("AVASP_flag_CONVERT_UNIT_CELL_INCOMPACT") << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.aopts.flag(\"AVASP_flag_CONVERT_UNIT_CELL_MINKOWSKI\")=" << xvasp.aopts.flag("AVASP_flag_CONVERT_UNIT_CELL_MINKOWSKI") << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.aopts.flag(\"AVASP_flag_CONVERT_UNIT_CELL_NIGGLI\")=" << xvasp.aopts.flag("AVASP_flag_CONVERT_UNIT_CELL_NIGGLI") << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.aopts.flag(\"AVASP_flag_CONVERT_UNIT_CELL_PRESERVE\")=" << xvasp.aopts.flag("AVASP_flag_CONVERT_UNIT_CELL_PRESERVE") << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.aopts.flag(\"AVASP_flag_CONVERT_UNIT_CELL_STANDARD_CONVENTIONAL\")=" << xvasp.aopts.flag("AVASP_flag_CONVERT_UNIT_CELL_STANDARD_CONVENTIONAL") << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.aopts.flag(\"AVASP_flag_CONVERT_UNIT_CELL_STANDARD_PRIMITIVE\")=" << xvasp.aopts.flag("AVASP_flag_CONVERT_UNIT_CELL_STANDARD_PRIMITIVE") << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.aopts.flag(\"AVASP_flag_CONVERT_UNIT_CELL_WIGNERSEITZ\")=" << xvasp.aopts.flag("AVASP_flag_CONVERT_UNIT_CELL_WIGNERSEITZ") << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.aopts.flag(\"FLAG::ABMIX_SET\")=" << xvasp.aopts.flag("FLAG::ABMIX_SET") << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.aopts.flag(\"FLAG::ALGO_PRESERVED\")=" << xvasp.aopts.flag("FLAG::ALGO_PRESERVED") << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.aopts.flag(\"FLAG::ALGO_SET\")=" << xvasp.aopts.flag("FLAG::ALGO_SET") << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.aopts.flag(\"FLAG::AVASP_AUTO_MAGMOM\")=" << xvasp.aopts.flag("FLAG::AVASP_AUTO_MAGMOM") << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.aopts.flag(\"FLAG::AVASP_AUTO_PSEUDOPOTENTIALS\")=" << xvasp.aopts.flag("FLAG::AVASP_AUTO_PSEUDOPOTENTIALS") << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.aopts.flag(\"FLAG::AVASP_BADER\")=" << xvasp.aopts.flag("FLAG::AVASP_BADER") << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.aopts.flag(\"FLAG::AVASP_CHGCAR\")=" << xvasp.aopts.flag("FLAG::AVASP_CHGCAR") << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.aopts.flag(\"FLAG::AVASP_ELF\")=" << xvasp.aopts.flag("FLAG::AVASP_ELF") << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.aopts.flag(\"FLAG::AVASP_FORCE_NOLDAU\")=" << xvasp.aopts.flag("FLAG::AVASP_FORCE_NOLDAU") << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.aopts.flag(\"FLAG::AVASP_KPPRA\")=" << xvasp.aopts.flag("FLAG::AVASP_KPPRA") << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.aopts.flag(\"FLAG::AVASP_LDAU1\")=" << xvasp.aopts.flag("FLAG::AVASP_LDAU1") << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.aopts.flag(\"FLAG::AVASP_LDAU2\")=" << xvasp.aopts.flag("FLAG::AVASP_LDAU2") << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.aopts.flag(\"FLAG::AVASP_LDAU_ADIABATIC\")=" << xvasp.aopts.flag("FLAG::AVASP_LDAU_ADIABATIC") << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.aopts.flag(\"FLAG::AVASP_LDAU_CUTOFF\")=" << xvasp.aopts.flag("FLAG::AVASP_LDAU_CUTOFF") << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.aopts.flag(\"FLAG::AVASP_LSCOUPLING\")=" << xvasp.aopts.flag("FLAG::AVASP_LSCOUPLING") << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.aopts.flag(\"FLAG::AVASP_RELAX_FORCES\")=" << xvasp.aopts.flag("FLAG::AVASP_RELAX_FORCES") << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.aopts.flag(\"FLAG::AVASP_SKIP_NOMIX\")=" << xvasp.aopts.flag("FLAG::AVASP_SKIP_NOMIX") << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.aopts.flag(\"FLAG::AVASP_SPIN\")=" << xvasp.aopts.flag("FLAG::AVASP_SPIN") << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.aopts.flag(\"FLAG::AVASP_SPIN_REMOVE_RELAX_1\")=" << xvasp.aopts.flag("FLAG::AVASP_SPIN_REMOVE_RELAX_1") << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.aopts.flag(\"FLAG::AVASP_SPIN_REMOVE_RELAX_2\")=" << xvasp.aopts.flag("FLAG::AVASP_SPIN_REMOVE_RELAX_2") << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.aopts.flag(\"FLAG::AVASP_WAVECAR\")=" << xvasp.aopts.flag("FLAG::AVASP_WAVECAR") << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.aopts.flag(\"FLAG::EXTRA_INCAR\")=" << xvasp.aopts.flag("FLAG::EXTRA_INCAR") << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.aopts.flag(\"FLAG::IVDW_SET\")=" << xvasp.aopts.flag("FLAG::IVDW_SET") << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.aopts.flag(\"FLAG::PRECISION_PRESERVED\")=" << xvasp.aopts.flag("FLAG::PRECISION_PRESERVED") << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.aopts.flag(\"FLAG::PRECISION_SET\")=" << xvasp.aopts.flag("FLAG::PRECISION_SET") << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.aopts.flag(\"FLAG::VOLUME_PRESERVED\")=" << xvasp.aopts.flag("FLAG::VOLUME_PRESERVED") << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.aopts.flag(\"FLAGS::AVASP_AAPL=OFF\")=" << xvasp.aopts.flag("FLAGS::AVASP_AAPL=OFF") << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.aopts.flag(\"FLAGS::AVASP_APL=OFF\")=" << xvasp.aopts.flag("FLAGS::AVASP_APL=OFF") << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.aopts.flag(\"FLAGS::AVASP_NEIGHBOURS=OFF\")=" << xvasp.aopts.flag("FLAGS::AVASP_NEIGHBOURS=OFF") << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.aopts.flag(\"FLAGS::AVASP_QHA=OFF\")=" << xvasp.aopts.flag("FLAGS::AVASP_QHA=OFF") << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.aopts.flag(\"FLAGS::AVASP_SYMMMETRY=OFF\")=" << xvasp.aopts.flag("FLAGS::AVASP_SYMMMETRY=OFF") << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.str=" << endl << xvasp.str << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.str.bravais_lattice_type=" << xvasp.str.bravais_lattice_type << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.str.bravais_lattice_variation_type=" << xvasp.str.bravais_lattice_variation_type << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.str.comp_each_type=";for(uint i=0;i<xvasp.str.comp_each_type.size();i++){cerr << xvasp.str.comp_each_type[i];} cerr << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.str.num_each_type=";for(uint i=0;i<xvasp.str.num_each_type.size();i++){cerr << xvasp.str.num_each_type[i];} cerr << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.str.spacegroupnumber=" << xvasp.str.spacegroupnumber << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.str.species=";for(uint i=0;i<xvasp.str.species.size();i++){cerr << xvasp.str.species[i];} cerr << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.str.species_mass=";for(uint i=0;i<xvasp.str.species_mass.size();i++){cerr << xvasp.str.species_mass[i];} cerr << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.str.species_pp=";for(uint i=0;i<xvasp.str.species_pp.size();i++){cerr << xvasp.str.species_pp[i];} cerr << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.str.species_pp_ZVAL=";for(uint i=0;i<xvasp.str.species_pp_ZVAL.size();i++){cerr << xvasp.str.species_pp_ZVAL[i];} cerr << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.str.species_pp_type=";for(uint i=0;i<xvasp.str.species_pp_type.size();i++){cerr << xvasp.str.species_pp_type[i];} cerr << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.str.species_pp_vLDAU=";for(uint i=0;i<xvasp.str.species_pp_vLDAU.size();i++){for(uint j=0;j<xvasp.str.species_pp_vLDAU[i].size();j++){cerr << xvasp.str.species_pp_vLDAU[i][j];} cerr << endl;} cerr << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.str.species_pp_version=";for(uint i=0;i<xvasp.str.species_pp_version.size();i++){cerr << xvasp.str.species_pp_version[i];} cerr << endl;}
    if(LDEBUG){cerr << soliloquy << " xvasp.str.species_volume=";for(uint i=0;i<xvasp.str.species_volume.size();i++){cerr << xvasp.str.species_volume[i];} cerr << endl;}

    if(LDEBUG){cerr << soliloquy << " END" << endl;}

    return true;
}

bool AVASP_MakeSingleAFLOWIN(_xvasp& xvasp_in,stringstream &_aflowin,bool flag_WRITE,int pthread,bool flag_PRINT) {
    //  if(flag_WRITE==FALSE) DEBUG=TRUE;
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string soliloquy="AVASP_MakeSingleAFLOWIN():";
    stringstream message;

    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " " << "[0]" << endl;
    if(LDEBUG && xvasp_in.AVASP_prototype_mode==LIBRARY_MODE_HTQC)       cerr << "DEBUG - " << soliloquy << " xvasp_in.AVASP_prototype_mode==LIBRARY_MODE_HTQC" << endl;
    if(LDEBUG && xvasp_in.AVASP_prototype_mode==LIBRARY_MODE_ICSD)       cerr << "DEBUG - " << soliloquy << " xvasp_in.AVASP_prototype_mode==LIBRARY_MODE_ICSD" << endl;
    if(LDEBUG && xvasp_in.AVASP_prototype_mode==LIBRARY_MODE_HTQC_ICSD)  cerr << "DEBUG - " << soliloquy << " xvasp_in.AVASP_prototype_mode==LIBRARY_MODE_HTQC_ICSD" << endl;
    if(LDEBUG && xvasp_in.AVASP_prototype_mode==LIBRARY_MODE_LIB3)       cerr << "DEBUG - " << soliloquy << " xvasp_in.AVASP_prototype_mode==LIBRARY_MODE_LIB3" << endl;
    if(LDEBUG && xvasp_in.AVASP_prototype_mode==LIBRARY_MODE_PROTOTYPE)  cerr << "DEBUG - " << soliloquy << " xvasp_in.AVASP_prototype_mode==LIBRARY_MODE_PROTOTYPE" << endl;
    if(LDEBUG && xvasp_in.AVASP_prototype_mode==LIBRARY_MODE_XSTRUCTURE) cerr << "DEBUG - " << soliloquy << " xvasp_in.AVASP_prototype_mode==LIBRARY_MODE_XSTRUCTURE" << endl;

    _xvasp xvasp=xvasp_in; // copy

    // cerr << "[6] xvasp.AVASP_potential=" << xvasp.AVASP_potential << endl; exit(0);

    // xvasp.aopts.flag("AFLOWIN_FLAG::PSTRESS",TRUE);
    //  if(xvasp.aopts.flag("AFLOWIN_FLAG::PSTRESS",TRUE)) xvasp.aopts.push_attached("AFLOWIN_FLAG::PSTRESS","7.65");
    // xvasp.AVASP_value_PSTRESS=7.65;

    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " xvasp.aopts.flag(\"AFLOWIN_FLAG::VASP\")=" << xvasp.aopts.flag("AFLOWIN_FLAG::VASP") << endl;
    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " xvasp.aopts.flag(\"FLAG::AVASP_FORCE_NOLDAU\")=" << xvasp.aopts.flag("FLAG::AVASP_FORCE_NOLDAU") << endl;
    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " xvasp.aopts.flag(\"FLAG::AVASP_FORCE_NOLDAU\")=" << xvasp.aopts.flag("FLAG::AVASP_FORCE_NOLDAU") << endl;
    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " xvasp.aopts.flag(\"AFLOWIN_FLAG::PSTRESS\")=" << xvasp.aopts.flag("AFLOWIN_FLAG::PSTRESS") << endl;
    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " xvasp.aopts.getattachedscheme(\"AFLOWIN_FLAG::PSTRESS\")=" << xvasp.aopts.getattachedscheme("AFLOWIN_FLAG::PSTRESS") << endl;
    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " xvasp.aopts.flag(\"AFLOWIN_FLAG::MODULE\")=" << xvasp.aopts.flag("AFLOWIN_FLAG::MODULE") << endl;  // CO 180214
    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " xvasp.aopts.getattachedscheme(\"AFLOWIN_FLAG::MODULE\")=" << xvasp.aopts.getattachedscheme("AFLOWIN_FLAG::MODULE") << endl;  // CO 180214
    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " xvasp.aopts.flag(\"AFLOWIN_FLAG::APL_SUPERCELL\")=" << xvasp.aopts.flag("AFLOWIN_FLAG::APL_SUPERCELL") << endl;  // CO 180214
    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " xvasp.aopts.getattachedscheme(\"AFLOWIN_FLAG::APL_SUPERCELL\")=" << xvasp.aopts.getattachedscheme("AFLOWIN_FLAG::APL_SUPERCELL") << endl;  // CO 180214
    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " xvasp.aopts.flag(\"AFLOWIN_FLAG::EDIFFG\")=" << xvasp.aopts.flag("AFLOWIN_FLAG::EDIFFG") << endl;
    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " xvasp.aopts.getattachedscheme(\"AFLOWIN_FLAG::EDIFFG\")=" << xvasp.aopts.getattachedscheme("AFLOWIN_FLAG::EDIFFG") << endl;
    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " xvasp.aopts.flag(\"AFLOWIN_FLAG::POTIM\")=" << xvasp.aopts.flag("AFLOWIN_FLAG::POTIM") << endl;
    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " xvasp.aopts.getattachedscheme(\"AFLOWIN_FLAG::POTIM\")=" << xvasp.aopts.getattachedscheme("AFLOWIN_FLAG::POTIM") << endl;
    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " xvasp.aopts.flag(\"AFLOWIN_FLAG::PRECISION\")=" << xvasp.aopts.flag("AFLOWIN_FLAG::PRECISION") << endl;
    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " xvasp.aopts.getattachedscheme(\"AFLOWIN_FLAG::PRECISION\")=" << xvasp.aopts.getattachedscheme("AFLOWIN_FLAG::PRECISION") << endl;
    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " xvasp.aopts.flag(\"AFLOWIN_FLAG::ALGORITHM\")=" << xvasp.aopts.flag("AFLOWIN_FLAG::ALGORITHM") << endl;
    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " xvasp.aopts.getattachedscheme(\"AFLOWIN_FLAG::ALGORITHM\")=" << xvasp.aopts.getattachedscheme("AFLOWIN_FLAG::ALGORITHM") << endl;
    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " xvasp.aopts.flag(\"AFLOWIN_FLAG::METAGGA\")=" << xvasp.aopts.flag("AFLOWIN_FLAG::METAGGA") << endl;
    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " xvasp.aopts.getattachedscheme(\"AFLOWIN_FLAG::METAGGA\")=" << xvasp.aopts.getattachedscheme("AFLOWIN_FLAG::METAGGA") << endl;
    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " xvasp.aopts.flag(\"AFLOWIN_FLAG::IVDW\")=" << xvasp.aopts.flag("AFLOWIN_FLAG::IVDW") << endl;
    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " xvasp.aopts.getattachedscheme(\"AFLOWIN_FLAG::IVDW\")=" << xvasp.aopts.getattachedscheme("AFLOWIN_FLAG::IVDW") << endl;
    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " xvasp.aopts.flag(\"AFLOWIN_FLAG::RELAX_TYPE\")=" << xvasp.aopts.flag("AFLOWIN_FLAG::RELAX_TYPE") << endl;  // CO 180214
    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " xvasp.aopts.getattachedscheme(\"AFLOWIN_FLAG::RELAX_TYPE\")=" << xvasp.aopts.getattachedscheme("AFLOWIN_FLAG::RELAX_TYPE") << endl;  // CO 180214
    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " xvasp.aopts.flag(\"AFLOWIN_FLAG::RELAX_MODE\")=" << xvasp.aopts.flag("AFLOWIN_FLAG::RELAX_MODE") << endl;
    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " xvasp.aopts.getattachedscheme(\"AFLOWIN_FLAG::RELAX_MODE\")=" << xvasp.aopts.getattachedscheme("AFLOWIN_FLAG::RELAX_MODE") << endl;
    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " xvasp.aopts.flag(\"AFLOWIN_FLAG::TYPE\")=" << xvasp.aopts.flag("AFLOWIN_FLAG::TYPE") << endl;
    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " xvasp.aopts.getattachedscheme(\"AFLOWIN_FLAG::TYPE\")=" << xvasp.aopts.getattachedscheme("AFLOWIN_FLAG::TYPE") << endl;
    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " xvasp.aopts.flag(\"AFLOWIN_FLAG::CONVERT_UNIT_CELL\")=" << xvasp.aopts.flag("AFLOWIN_FLAG::CONVERT_UNIT_CELL") << endl;
    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " xvasp.aopts.getattachedscheme(\"AFLOWIN_FLAG::CONVERT_UNIT_CELL\")=" << xvasp.aopts.getattachedscheme("AFLOWIN_FLAG::CONVERT_UNIT_CELL") << endl;
    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " xvasp.aopts.flag(\"AFLOWIN_FLAG::VOLUME_PLUS_EQUAL\")=" << xvasp.aopts.flag("AFLOWIN_FLAG::VOLUME_PLUS_EQUAL") << endl;
    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " xvasp.aopts.getattachedscheme(\"AFLOWIN_FLAG::VOLUME_PLUS_EQUAL\")=" << xvasp.aopts.getattachedscheme("AFLOWIN_FLAG::VOLUME_PLUS_EQUAL") << endl;
    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " xvasp.aopts.flag(\"AFLOWIN_FLAG::VOLUME_MULTIPLY_EQUAL\")=" << xvasp.aopts.flag("AFLOWIN_FLAG::VOLUME_MULTIPLY_EQUAL") << endl;
    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " xvasp.aopts.getattachedscheme(\"AFLOWIN_FLAG::VOLUME_MULTIPLY_EQUAL\")=" << xvasp.aopts.getattachedscheme("AFLOWIN_FLAG::VOLUME_MULTIPLY_EQUAL") << endl;
    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " xvasp.aopts.flag(\"AFLOWIN_FLAG::NO_VOLUME_ADJUSTMENT\")=" << xvasp.aopts.flag("AFLOWIN_FLAG::NO_VOLUME_ADJUSTMENT") << endl;  // CO 180214

    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " xvasp.aopts.flag(\"AFLOWIN_FLAG::KPPRA\")=" << xvasp.aopts.flag("AFLOWIN_FLAG::KPPRA") << endl;
    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " xvasp.aopts.getattachedscheme(\"AFLOWIN_FLAG::KPPRA\")=" << xvasp.aopts.getattachedscheme("AFLOWIN_FLAG::KPPRA") << endl;
    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " xvasp.aopts.flag(\"AFLOWIN_FLAG::KPPRA_STATIC\")=" << xvasp.aopts.flag("AFLOWIN_FLAG::KPPRA_STATIC") << endl;
    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " xvasp.aopts.getattachedscheme(\"AFLOWIN_FLAG::KPPRA_STATIC\")=" << xvasp.aopts.getattachedscheme("AFLOWIN_FLAG::KPPRA_STATIC") << endl;
    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " xvasp.aopts.flag(\"AFLOWIN_FLAG::BANDS_GRID\")=" << xvasp.aopts.flag("AFLOWIN_FLAG::BANDS_GRID") << endl;
    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " xvasp.aopts.getattachedscheme(\"AFLOWIN_FLAG::BANDS_GRID\")=" << xvasp.aopts.getattachedscheme("AFLOWIN_FLAG::BANDS_GRID") << endl;

    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " xvasp.aopts.flag(\"AFLOWIN_FLAG::ENMAX_MULTYPLY\")=" << xvasp.aopts.flag("AFLOWIN_FLAG::ENMAX_MULTYPLY") << endl;
    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " xvasp.aopts.getattachedscheme(\"AFLOWIN_FLAG::ENMAX_MULTYPLY\")=" << xvasp.aopts.getattachedscheme("AFLOWIN_FLAG::ENMAX_MULTYPLY") << endl;
    // if(LDEBUG) exit(0);

    //  cerr << aurostd::vectorstring2string(xvasp.str.species) << "/" << xvasp.AVASP_label << endl;
    // checks if everything is different
    ofstream oaus;
    //  xstructure str;
    string directory,system,formula,label_MIX,label_NOMIX,label_SKIP,label_ALREADY;
    bool Krun=TRUE; //,swapspeciesAB=FALSE,swapspeciesBC=FALSE,swapspeciesCA=FALSE;
    stringstream aflowin;aflowin.clear();aflowin.str(std::string());
    stringstream aflowin_qe; aflowin_qe.clear();aflowin_qe.str(std::string());
    stringstream aflowin_abinit; aflowin_abinit.clear();aflowin_abinit.str(std::string());
    stringstream aflowin_aims; aflowin_aims.clear();aflowin_aims.str(std::string());
    // for HTQC
    // xstructure str1[100];_xvasp xvsp[100];

    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " " << "[1]" << endl;
    //clean up name
    xvasp.AVASP_label=aurostd::RemoveSubString(xvasp.AVASP_label,"/");
    xvasp.AVASP_label=aurostd::RemoveSubString(xvasp.AVASP_label," ");
    // loading BUGS
    vector<string> vstr;
    aurostd::string2tokens(string(CONVENTIONAL_FORCE),vstr,",");
    bool fixed_lattice_vasp=FALSE;
    for(uint istr=0;istr<vstr.size()&&!fixed_lattice_vasp;istr++) {
        if(vstr.at(istr)==xvasp.AVASP_label) {
            //   cerr << vstr.at(istr) << endl;
            fixed_lattice_vasp=TRUE;
            if(xvasp.aopts.flag("AVASP_flag_CONVERT_UNIT_CELL_STANDARD_PRIMITIVE")) {
                xvasp.aopts.flag("AVASP_flag_CONVERT_UNIT_CELL_STANDARD_PRIMITIVE",FALSE);
                xvasp.aopts.flag("AVASP_flag_CONVERT_UNIT_CELL_STANDARD_CONVENTIONAL",TRUE);
            }
        }
    }

    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " " << "[2]" << endl;

    if(LDEBUG) {
        cerr << "xvasp.str.species.size()=" << xvasp.str.species.size() << endl;
        cerr << "xvasp.str.species_pp.size()=" << xvasp.str.species_pp.size() << "   ";
        for(uint i=0;i<xvasp.str.species_pp.size();i++) cerr << xvasp.str.species_pp.at(i) << " "; 
        cerr  << endl;
        cerr << "xvasp.str.species_pp_type.size()=" << xvasp.str.species_pp_type.size() << "   ";
        for(uint i=0;i<xvasp.str.species_pp_type.size();i++) cerr << xvasp.str.species_pp_type.at(i) << " "; 
        cerr  << endl;
        cerr << "xvasp.str.species_pp_version.size()=" << xvasp.str.species_pp_version.size() << "   ";
        for(uint i=0;i<xvasp.str.species_pp_version.size();i++) cerr << xvasp.str.species_pp_version.at(i) << " "; 
        cerr  << endl;
        cerr << "xvasp.str.species_pp_ZVAL.size()=" << xvasp.str.species_pp_ZVAL.size() << "   ";
        for(uint i=0;i<xvasp.str.species_pp_ZVAL.size();i++) cerr << xvasp.str.species_pp_ZVAL.at(i) << " "; 
        cerr  << endl;
        cerr << "xvasp.str.species_pp_vLDAU.size()=" << xvasp.str.species_pp_vLDAU.size() << "   ";
        for(uint i=0;i<xvasp.str.species_pp_vLDAU.size();i++) cerr << xvasp.str.species_pp_vLDAU.at(i).size() << " "; 
        cerr  << endl;
        cerr << "xvasp.str.num_each_type.size()=" << xvasp.str.num_each_type.size() << endl;
        cerr << "xvasp.str.comp_each_type.size()=" << xvasp.str.comp_each_type.size() << endl;
    }
    if(xvasp.str.num_each_type.size()!=xvasp.str.species.size()) {
        if(LDEBUG) cerr << soliloquy << " WARNING 1 (fixing)" << endl;
        xvasp.str.num_each_type.clear(); 
        for(uint i=0;i<xvasp.str.species.size();i++) xvasp.str.num_each_type.push_back(0);
        xvasp.str.comp_each_type.clear(); 
        for(uint i=0;i<xvasp.str.species.size();i++) xvasp.str.comp_each_type.push_back(0);
    }
    if(xvasp.str.species_pp.size()!=xvasp.str.species.size()) {
        if(LDEBUG) cerr <<"" << soliloquy << " WARNING 2 (fixing)" << endl;
        xvasp.str.species_pp=xvasp.str.species; // so we have the same dimension !
    }
    if(xvasp.str.species_volume.size()!=xvasp.str.species.size()) {
        if(LDEBUG) cerr <<"" << soliloquy << " WARNING 3 (fixing)" << endl;
        xvasp.str.species_volume.clear(); 
        for(uint i=0;i<xvasp.str.species.size();i++) xvasp.str.species_volume.push_back(0);
    }
    if(xvasp.str.species_mass.size()!=xvasp.str.species.size()) {
        if(LDEBUG) cerr <<"" << soliloquy << " WARNING 4 (fixing)" << endl;
        xvasp.str.species_mass.clear(); 
        for(uint i=0;i<xvasp.str.species.size();i++) xvasp.str.species_mass.push_back(0);
    }
    deque<double> vmass(xvasp.str.species_mass); // BACKUP
    //  for(uint i=0;i<xvasp.str.species.size();i++) cerr << "xvasp.str.species_mass=" << xvasp.str.species_mass.at(i) << endl;

    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " " << "[3]" << endl;

    //  cerr << "[7] xvasp.AVASP_potential=" << xvasp.AVASP_potential << endl; exit(0);

    // DEBUG for XBANDS put 1
    if(0) { // DEBUG XBANDS
        xvasp.aopts.flag("FLAG::AVASP_AUTO_PSEUDOPOTENTIALS",FALSE);
        if(xvasp.aopts.flag("FLAG::AVASP_AUTO_PSEUDOPOTENTIALS")==FALSE) {
            for(uint isp=0;isp<xvasp.str.species.size();isp++) {
                xvasp.str.species_pp.at(isp)=xvasp.str.species.at(isp);
            }
        }
    }

    //  cerr << "xvasp.str.species_pp.size()=" << xvasp.str.species_pp << endl;for(uint i=0;i<xvasp.str.species_pp.size();i++) cerr << xvasp.str.species_pp.at(i) << endl; exit(0);

    // FIX THE SPECIES IF AUTO_PSEUDOPOTENTIALS
    if(xvasp.aopts.flag("FLAG::AVASP_AUTO_PSEUDOPOTENTIALS")) {
        if(1 && xvasp.str.species.size()==2) // AlMg potpaw_GGA
            if((KBIN::VASP_PseudoPotential_CleanName(xvasp.str.species.at(0))=="Ga" && KBIN::VASP_PseudoPotential_CleanName(xvasp.str.species.at(1))=="Mg") ||
                    (KBIN::VASP_PseudoPotential_CleanName(xvasp.str.species.at(0))=="Ge" && KBIN::VASP_PseudoPotential_CleanName(xvasp.str.species.at(1))=="Mg") ||
                    (KBIN::VASP_PseudoPotential_CleanName(xvasp.str.species.at(0))=="Al" && KBIN::VASP_PseudoPotential_CleanName(xvasp.str.species.at(1))=="Mg") ||
                    (KBIN::VASP_PseudoPotential_CleanName(xvasp.str.species.at(0))=="Mg" && KBIN::VASP_PseudoPotential_CleanName(xvasp.str.species.at(1))=="Si")) {
                cout << soliloquy << " Pseudopotential exception: " << xvasp.str.species_pp.at(0) << xvasp.str.species_pp.at(1) << " -> potpaw_GGA " << endl;
                xvasp.AVASP_potential="potpaw_GGA";
            }
        // FIX THE SPECIES IF AUTO_PSEUDOPOTENTIALS && LIBRARY_MODE_HTQC
        if(xvasp.aopts.flag("FLAG::AVASP_AUTO_PSEUDOPOTENTIALS") && (xvasp.AVASP_prototype_mode==LIBRARY_MODE_HTQC || xvasp.AVASP_prototype_mode==LIBRARY_MODE_HTQC_ICSD)) {
            // patches for potpaw_GGA
            if(xvasp.AVASP_potential==DEFAULT_VASP_POTCAR_DIR_POTPAW_GGA) {
                for(uint isp=0;isp<xvasp.str.species.size();isp++)
                    xvasp.str.species_pp.at(isp)=AVASP_Get_PseudoPotential_PAW_GGA(KBIN::VASP_PseudoPotential_CleanName(xvasp.str.species.at(isp)));
                if(xvasp.str.species.size()==2) {
                    // check for EXCEPTION GaMg
                    if(KBIN::VASP_PseudoPotential_CleanName(xvasp.str.species.at(0))=="Ga" && KBIN::VASP_PseudoPotential_CleanName(xvasp.str.species.at(1))=="Mg") {
                        cout << soliloquy << " Pseudopotential exception: " << xvasp.str.species_pp.at(0) << xvasp.str.species_pp.at(1) << " -> ";
                        xvasp.str.species_pp.at(0)="Ga";xvasp.str.species_pp.at(1)="Mg_pv";cout << xvasp.str.species_pp.at(0) << xvasp.str.species_pp.at(1) << endl;
                    }
                    // check for EXCEPTION GeMg
                    if(KBIN::VASP_PseudoPotential_CleanName(xvasp.str.species.at(0))=="Ge" && KBIN::VASP_PseudoPotential_CleanName(xvasp.str.species.at(1))=="Mg") {
                        cout << soliloquy << " Pseudopotential exception: " << xvasp.str.species_pp.at(0) << xvasp.str.species_pp.at(1) << " -> ";
                        xvasp.str.species_pp.at(0)="Ge";xvasp.str.species_pp.at(1)="Mg_pv";cout << xvasp.str.species_pp.at(0) << xvasp.str.species_pp.at(1) << endl;
                    }
                    if(KBIN::VASP_PseudoPotential_CleanName(xvasp.str.species.at(0))=="Mg" && KBIN::VASP_PseudoPotential_CleanName(xvasp.str.species.at(1))=="Si") {
                        cout << soliloquy << " Pseudopotential exception: " << xvasp.str.species_pp.at(0) << xvasp.str.species_pp.at(1) << " -> ";
                        xvasp.str.species_pp.at(0)="Mg_pv";xvasp.str.species_pp.at(1)="Si";cout << xvasp.str.species_pp.at(0) << xvasp.str.species_pp.at(1) << endl;
                    }
                }
            }
            // patches for potpaw_PBE
            if(xvasp.AVASP_potential==DEFAULT_VASP_POTCAR_DIR_POTPAW_PBE) {
                for(uint isp=0;isp<xvasp.str.species.size();isp++){
                    xvasp.str.species_pp.at(isp)=AVASP_Get_PseudoPotential_PAW_PBE(KBIN::VASP_PseudoPotential_CleanName(xvasp.str.species.at(isp)));
                    if(LDEBUG){cerr << soliloquy << " xvasp.str.species_pp.at(" << isp << ")=" << xvasp.str.species_pp.at(isp) << endl;}  //CO 180705
                }
                if(xvasp.str.species.size()==2) {
                    // check for EXCEPTION BHf
                    if(KBIN::VASP_PseudoPotential_CleanName(xvasp.str.species.at(0))=="B" && KBIN::VASP_PseudoPotential_CleanName(xvasp.str.species.at(1))=="Hf") {
                        cout << soliloquy << " Pseudopotential exception: " << xvasp.str.species_pp.at(0) << xvasp.str.species_pp.at(1) << " -> ";
                        xvasp.str.species_pp.at(0)="B_s";xvasp.str.species_pp.at(1)="Hf_pv";cout << xvasp.str.species_pp.at(0) << xvasp.str.species_pp.at(1) << endl;
                    }
                    // check for EXCEPTION BeHf
                    if(KBIN::VASP_PseudoPotential_CleanName(xvasp.str.species.at(0))=="Be" && KBIN::VASP_PseudoPotential_CleanName(xvasp.str.species.at(1))=="Hf") {
                        cout << soliloquy << " Pseudopotential exception: " << xvasp.str.species_pp.at(0) << xvasp.str.species_pp.at(1) << " -> ";
                        xvasp.str.species_pp.at(0)="Be";xvasp.str.species_pp.at(1)="Hf_pv";cout << xvasp.str.species_pp.at(0) << xvasp.str.species_pp.at(1) << endl;
                    }
                    // check for EXCEPTION HfK
                    if(KBIN::VASP_PseudoPotential_CleanName(xvasp.str.species.at(0))=="Hf" && KBIN::VASP_PseudoPotential_CleanName(xvasp.str.species.at(1))=="K") {
                        cout << soliloquy << " Pseudopotential exception: " << xvasp.str.species_pp.at(0) << xvasp.str.species_pp.at(1) << " -> ";
                        xvasp.str.species_pp.at(0)="Hf_pv";xvasp.str.species_pp.at(1)="K_pv";cout << xvasp.str.species_pp.at(0) << xvasp.str.species_pp.at(1) << endl;
                    }
                    // check for EXCEPTION HfSn
                    if(KBIN::VASP_PseudoPotential_CleanName(xvasp.str.species.at(0))=="Hf" && KBIN::VASP_PseudoPotential_CleanName(xvasp.str.species.at(1))=="Sn") {
                        cout << soliloquy << " Pseudopotential exception: " << xvasp.str.species_pp.at(0) << xvasp.str.species_pp.at(1) << " -> ";
                        xvasp.str.species_pp.at(0)="Hf_pv";xvasp.str.species_pp.at(1)="Sn_d";cout << xvasp.str.species_pp.at(0) << xvasp.str.species_pp.at(1) << endl;
                    }
                    // check for EXCEPTION BSm
                    if(KBIN::VASP_PseudoPotential_CleanName(xvasp.str.species.at(0))=="B" && KBIN::VASP_PseudoPotential_CleanName(xvasp.str.species.at(1))=="Sm") {
                        cout << soliloquy << " Pseudopotential exception: " << xvasp.str.species_pp.at(0) << xvasp.str.species_pp.at(1) << " -> ";
                        xvasp.str.species_pp.at(0)="B_h";xvasp.str.species_pp.at(1)="Sm_3";cout << xvasp.str.species_pp.at(0) << xvasp.str.species_pp.at(1) << endl;
                        aurostd::StringSubst(xvasp.AVASP_potential,"potpaw_PBE","potpaw_GGA");
                        cout << soliloquy << " Pseudopotential exception: \"potpaw_PBE\"=>\"potpaw_GGA\"" << endl;
                        xvasp.aopts.flag("FLAG::AVASP_AUTO_PSEUDOPOTENTIALS",FALSE);
                        cout << soliloquy << " turning OFF pseudopotential" << endl;	    
                    }
                }
            }
            // patches for potpaw_PBE_KIN
            if(xvasp.AVASP_potential==DEFAULT_VASP_POTCAR_DIR_POTPAW_PBE_KIN) {
                for(uint isp=0;isp<xvasp.str.species.size();isp++)
                    xvasp.str.species_pp.at(isp)=AVASP_Get_PseudoPotential_PAW_PBE_KIN(KBIN::VASP_PseudoPotential_CleanName(xvasp.str.species.at(isp)));
                // no exceptions yet
            }
            // patches for potpaw_LDA_KIN
            if(xvasp.AVASP_potential==DEFAULT_VASP_POTCAR_DIR_POTPAW_LDA_KIN) {
                for(uint isp=0;isp<xvasp.str.species.size();isp++)
                    xvasp.str.species_pp.at(isp)=AVASP_Get_PseudoPotential_PAW_LDA_KIN(KBIN::VASP_PseudoPotential_CleanName(xvasp.str.species.at(isp)));
                // no exceptions yet
            }
            // patches for potpaw_LDA pot_GGA pot_LDA
            for(uint isp=0;isp<xvasp.str.species.size();isp++)  {
                //	if(aurostd::substring2bool(xvasp.AVASP_potential,"potpaw_GGA"))
                //	  xvasp.str.species_pp.at(isp)=AVASP_Get_PseudoPotential_PAW_GGA(KBIN::VASP_PseudoPotential_CleanName(xvasp.str.species.at(isp)));
                if(xvasp.AVASP_potential==DEFAULT_VASP_POTCAR_DIR_POTPAW_LDA)
                    xvasp.str.species_pp.at(isp)=AVASP_Get_PseudoPotential_PAW_LDA(KBIN::VASP_PseudoPotential_CleanName(xvasp.str.species.at(isp)));
                if(xvasp.AVASP_potential==DEFAULT_VASP_POTCAR_DIR_POT_GGA)
                    xvasp.str.species_pp.at(isp)=AVASP_Get_PseudoPotential_GGA(KBIN::VASP_PseudoPotential_CleanName(xvasp.str.species.at(isp)));
                if(xvasp.AVASP_potential==DEFAULT_VASP_POTCAR_DIR_POT_LDA)
                    xvasp.str.species_pp.at(isp)=AVASP_Get_PseudoPotential_LDA(KBIN::VASP_PseudoPotential_CleanName(xvasp.str.species.at(isp)));
            }
        }
        // FIX THE SPECIES IF   AUTO_PSEUDOPOTENTIALS && (LIBRARY_MODE_ICSD ||LIBRARY_MODE_LIB3)
        if(xvasp.aopts.flag("FLAG::AVASP_AUTO_PSEUDOPOTENTIALS") && (xvasp.AVASP_prototype_mode==LIBRARY_MODE_ICSD || xvasp.AVASP_prototype_mode==LIBRARY_MODE_LIB3)) {
            for(uint isp=0;isp<xvasp.str.species.size();isp++) {
                if(xvasp.AVASP_potential==DEFAULT_VASP_POTCAR_DIR_POTPAW_PBE_KIN)
                    xvasp.str.species_pp.at(isp)=AVASP_Get_PseudoPotential_PAW_PBE_KIN(KBIN::VASP_PseudoPotential_CleanName(xvasp.str.species.at(isp)));
                if(xvasp.AVASP_potential==DEFAULT_VASP_POTCAR_DIR_POTPAW_LDA_KIN)
                    xvasp.str.species_pp.at(isp)=AVASP_Get_PseudoPotential_PAW_LDA_KIN(KBIN::VASP_PseudoPotential_CleanName(xvasp.str.species.at(isp)));
                if(xvasp.AVASP_potential==DEFAULT_VASP_POTCAR_DIR_POTPAW_PBE)
                    xvasp.str.species_pp.at(isp)=AVASP_Get_PseudoPotential_PAW_PBE(KBIN::VASP_PseudoPotential_CleanName(xvasp.str.species.at(isp)));
                if(xvasp.AVASP_potential==DEFAULT_VASP_POTCAR_DIR_POTPAW_GGA)
                    xvasp.str.species_pp.at(isp)=AVASP_Get_PseudoPotential_PAW_GGA(KBIN::VASP_PseudoPotential_CleanName(xvasp.str.species.at(isp)));
                if(xvasp.AVASP_potential==DEFAULT_VASP_POTCAR_DIR_POTPAW_LDA)
                    xvasp.str.species_pp.at(isp)=AVASP_Get_PseudoPotential_PAW_LDA(KBIN::VASP_PseudoPotential_CleanName(xvasp.str.species.at(isp)));
                if(xvasp.AVASP_potential==DEFAULT_VASP_POTCAR_DIR_POT_GGA)
                    xvasp.str.species_pp.at(isp)=AVASP_Get_PseudoPotential_GGA(KBIN::VASP_PseudoPotential_CleanName(xvasp.str.species.at(isp)));
                if(xvasp.AVASP_potential==DEFAULT_VASP_POTCAR_DIR_POT_LDA)
                    xvasp.str.species_pp.at(isp)=AVASP_Get_PseudoPotential_LDA(KBIN::VASP_PseudoPotential_CleanName(xvasp.str.species.at(isp)));
            }
        }

        //    cerr << "[8] xvasp.AVASP_potential=" << xvasp.AVASP_potential << endl; exit(0);

        // type and version fixed
        xvasp.str.species_pp_type.clear();for(uint isp=0;isp<xvasp.str.species.size();isp++) xvasp.str.species_pp_type.push_back(xvasp.AVASP_potential);
        xvasp.str.species_pp_version.clear();for(uint isp=0;isp<xvasp.str.species.size();isp++) xvasp.str.species_pp_version.push_back("");
        xvasp.str.species_pp_ZVAL.clear();for(uint isp=0;isp<xvasp.str.species.size();isp++) xvasp.str.species_pp_ZVAL.push_back(0.0);
        xvasp.str.species_pp_vLDAU.clear();for(uint isp=0;isp<xvasp.str.species.size();isp++) xvasp.str.species_pp_vLDAU.push_back(deque<double>());

    } // NAMES FIXED !

    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " " << "[4]" << endl;

    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " [5] xvasp.str.species.size()=" << xvasp.str.species.size() << endl;

    // cerr << "xvasp.str.species_pp.size()=" << xvasp.str.species_pp.size() << endl; for(uint i=0;i<xvasp.str.species_pp.size();i++) cerr << xvasp.str.species_pp.at(i) << endl; exit(0);

    // START THE HTQC STUFF
    if((xvasp.AVASP_prototype_mode==LIBRARY_MODE_HTQC || xvasp.AVASP_prototype_mode==LIBRARY_MODE_HTQC_ICSD || xvasp.AVASP_prototype_mode==LIBRARY_MODE_LIB3)) {
        // check for EXCEPTION Gd and Ce metals PUT LDAU
        if(aurostd::substring2bool(xvasp.str.species,"Gd")) {xvasp.aopts.flag("FLAG::AVASP_LDAU1",FALSE);xvasp.aopts.flag("FLAG::AVASP_LDAU2",TRUE);}
        if(aurostd::substring2bool(xvasp.str.species,"Ce")) {xvasp.aopts.flag("FLAG::AVASP_LDAU1",FALSE);xvasp.aopts.flag("FLAG::AVASP_LDAU2",TRUE);}
        if(LDEBUG)   cerr << "xvasp.aopts.flag(\"FLAG::AVASP_FORCE_NOLDAU\")=" << xvasp.aopts.flag("FLAG::AVASP_FORCE_NOLDAU") << endl;
        if(xvasp.aopts.flag("FLAG::AVASP_FORCE_NOLDAU"))  {xvasp.aopts.flag("FLAG::AVASP_LDAU1",FALSE);xvasp.aopts.flag("FLAG::AVASP_LDAU2",FALSE);}

        if(xvasp.AVASP_directory_from_library_) {
            for(uint isp1=0;isp1<xvasp.str.species.size();isp1++)
                for(uint isp2=isp1+1;isp2<xvasp.str.species.size();isp2++)
                    if(KBIN::VASP_PseudoPotential_CleanName(xvasp.str.species.at(isp1))==KBIN::VASP_PseudoPotential_CleanName(xvasp.str.species.at(isp2))) {
                        ostringstream aus;
                        if(LDEBUG) aus << "EEEEE  PROTO_Generation_Binary_Aflowin error xvasp.str.species.at(isp1)==xvasp.str.species.at(isp2) " << xvasp.str.species.at(isp1) << endl;
                        if(LDEBUG) aurostd::PrintErrorStream(oaus,aus,XHOST.QUIET);
                        Krun=FALSE;
                        return Krun; // empty
                    }
        }
        //  cerr << xvasp.AVASP_volumeC; // keep it bisy
        //KESONG FIXED THE DIRECTORY, 12/22/2013
        _xvasp xvasp_dir=xvasp;
        if(xvasp.AVASP_alpha_fix==TRUE) {
            xvasp_dir.str.SpeciesPutAlphabetic(); // should put alphabetic also the PSEUDOPOTENTIALS //fixed!!!!!! //KESONG 
        }
        //xvasp.str.SpeciesPutAlphabetic(); // should put alphabetic also the PSEUDOPOTENTIALS //fixed!!!!!! //KESONG 
        //Although this can make directory alphabetic, but it obfucates POSCAR!!!!!!!!


        if(xvasp.AVASP_prototype_from_library_) {
            string label_AFLOWLIB_LIB1,label_AFLOWLIB_LIB2,label_AFLOWLIB_LIB3,label_AFLOWLIB_LIB4,label_AFLOWLIB_LIB5,label_AFLOWLIB_LIB6,label_AFLOWLIB_LIB7,label_AFLOWLIB_LIB8,label_AFLOWLIB_LIB9;
            directory=xvasp.AVASP_dirbase+"/";
            for(uint i=0;i<xvasp.str.species_pp.size();i++) {
                directory+=xvasp_dir.str.species_pp.at(i);   //FIXED DIRECTORY, //KESONG
                if(xvasp.POTCAR_TYPE_DATE_PRINT_flag) { // add potential type and date
                    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " [5a.2]" << endl;
                    string FilePotcar,DataPotcar;
                    if(!KBIN::VASP_Find_FILE_POTCAR(xvasp.AVASP_potential+"/"+xvasp.str.species_pp.at(i),FilePotcar,DataPotcar)) {
                        cerr << "EEEEE  POTCAR [" << xvasp.AVASP_potential+"/"+xvasp.str.species_pp.at(i) << "] not found! " << Message("user,host,time") << endl;
                        return FALSE; // dont die
                        //exit(0);
                    } else {
                        if(LDEBUG) cerr << "DEBUG - " << soliloquy << "  POTCAR [" << xvasp.AVASP_potential+"/"+xvasp.str.species_pp.at(i) << "] = " << FilePotcar << endl;
                        vector<string> tokens;
                        string pottype="",date="",sgrep=aurostd::execute2string("cat "+FilePotcar+" | grep TITEL");
                        aurostd::string2tokens(sgrep,tokens," ");
                        if(tokens.size()!=4 && tokens.size()!=5) {
                            cerr << "EEEEE  POTCAR [" << xvasp.AVASP_potential+"/"+xvasp.str.species_pp.at(i) << "] = " << FilePotcar << "  wrong TITEL: " << sgrep << endl; 
                            return FALSE; // dont die
                            //exit(0);
                        }
                        if(tokens.size()==4) { // no indications of time
                            if(tokens.at(2)=="US" && aurostd::substring2bool(FilePotcar,"LDA","lda")) {pottype="LDA";date=DEFAULT_VASP_POTCAR_DATE_POT_LDA;}
                            if(tokens.at(2)=="US" && aurostd::substring2bool(FilePotcar,"GGA","gga")) {pottype="GGA";date=DEFAULT_VASP_POTCAR_DATE_POT_GGA;}
                        }
                        if(tokens.size()==5) {
                            if(tokens.at(2)=="PAW") {pottype="PAW_LDA";date=tokens.at(4);}
                            if(tokens.at(2)=="PAW_GGA") {pottype="PAW_GGA";date=tokens.at(4);}
                            if(tokens.at(2)=="PAW_RPBE") {pottype="PAW_RPBE";date=tokens.at(4);}  // potpaw_GGA/DEFAULT_VASP_POTCAR_DATE/Ge_h
                            if(tokens.at(2)=="PAW_PBE") {pottype="PAW_PBE";date=tokens.at(4);}
                            if(tokens.at(2)=="PAW_PBE") {pottype="PAW_PBE_KIN";date=tokens.at(4);} // FIX COREY/STEFANIO PBE_KIN CHECK PRESENCE OF "mkinetic energy-density pseudized"
                            if(tokens.at(2)=="PAW_LDA") {pottype="PAW_LDA_KIN";date=tokens.at(4);} // FIX COREY/STEFANIO LDA_KIN CHECK PRESENCE OF "mkinetic energy-density pseudized"
                            // SEE https://cms.mpi.univie.ac.at/wiki/index.php/METAGGA
                        }
                        if(pottype.empty()) {
                            cerr << "EEEEE  POTCAR [" << xvasp.AVASP_potential+"/"+xvasp.str.species_pp.at(i) << "] = " << FilePotcar << "  wrong pottype:" << sgrep << endl; 
                            return FALSE; // dont die
                            //  exit(0);
                        }
                        if(date.empty()) {
                            cerr << "EEEEE  POTCAR [" << xvasp.AVASP_potential+"/"+xvasp.str.species_pp.at(i) << "] = " << FilePotcar << "  wrong date:" << sgrep << endl;
                            return FALSE; // dont die
                            //  exit(0);
                        }
                        directory+=":"+pottype+":"+date;
                        if(i<xvasp.str.species_pp.size()-1) directory+=":";
                    }
                    if(LDEBUG) cerr << xvasp.AVASP_potential << "/" << xvasp.str.species_pp.at(i) << endl;	
                }
            }
            if(LDEBUG) cerr << "DEBUG - " << soliloquy << " [5a.3]" << endl;
            directory+="/"+xvasp.AVASP_label;
            //  cerr << "DEBUG - " << soliloquy << " xvasp.AVASP_prototype_from_library_ directory(1)=" << directory << endl;    
            system=""+aurostd::vectorstring2string(xvasp.str.species_pp)+"."+xvasp.AVASP_label;
            // label_text0 - MIX
            label_MIX=xvasp.AVASP_libbase+"/"+aurostd::vectorstring2string(xvasp.str.species_pp)+"/"+xvasp.AVASP_label+"/"+_AFLOWIN_;
            if(LDEBUG) cerr << "DEBUG - " << soliloquy << " [5a.4]" << endl;
            if(XHOST.hostname==XHOST.AFLOW_MATERIALS_SERVER || XHOST.hostname==XHOST.AFLOW_WEB_SERVER) {
                if(xvasp.AVASP_libbase.size()==0) {
                    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " [5a.5]" << endl;
                    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " [5a.6] XHOST_LIBRARY_LIB2=" << XHOST_LIBRARY_LIB2 << endl;
                    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " [5a.7] vAFLOW_PROJECTS_DIRECTORIES.size()=" << vAFLOW_PROJECTS_DIRECTORIES.size() << endl;
                    label_AFLOWLIB_LIB1=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB1)+"/RAW"+label_MIX;
                    label_AFLOWLIB_LIB2=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB2)+"/RAW"+label_MIX;
                    label_AFLOWLIB_LIB4=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB4)+"/RAW"+label_MIX;
                    label_AFLOWLIB_LIB5=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB5)+"/RAW"+label_MIX;
                    label_AFLOWLIB_LIB6=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB6)+"/RAW"+label_MIX;
                    label_AFLOWLIB_LIB7=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB7)+"/RAW"+label_MIX;
                    label_AFLOWLIB_LIB8=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB8)+"/RAW"+label_MIX;
                    label_AFLOWLIB_LIB9=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB9)+"/RAW"+label_MIX;
                    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " [5a.6] label_AFLOWLIB_LIB2=" << label_AFLOWLIB_LIB2 << endl;
                    if(XHOST_LIBRARY_LIB3!=LIBRARY_NOTHING)  {label_AFLOWLIB_LIB3=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB3)+"/LIB"+label_MIX;} else {label_AFLOWLIB_LIB3="./LIB"+label_MIX;}
                } else {
                    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " [5d.02] else CHECK THESE BY HAND" << endl;
                    label_AFLOWLIB_LIB2=label_MIX;aurostd::StringSubst(label_AFLOWLIB_LIB2,vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB3),vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB2));
                    label_AFLOWLIB_LIB3=label_MIX;aurostd::StringSubst(label_AFLOWLIB_LIB3,vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB2),vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB3));
                    aurostd::StringSubst(label_AFLOWLIB_LIB3,"LIB2U","");
                    aurostd::StringSubst(label_AFLOWLIB_LIB2,"//","/");
                    aurostd::StringSubst(label_AFLOWLIB_LIB3,"//","/");
                    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " [5d.02] label_AFLOWLIB_LIB1=" << label_AFLOWLIB_LIB1 << endl;
                    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " [5d.02] label_AFLOWLIB_LIB2=" << label_AFLOWLIB_LIB2 << endl;
                    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " [5d.02] label_AFLOWLIB_LIB3=" << label_AFLOWLIB_LIB3 << endl;	
                    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " [5d.02] label_AFLOWLIB_LIB4=" << label_AFLOWLIB_LIB4 << endl;
                    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " [5d.02] label_AFLOWLIB_LIB5=" << label_AFLOWLIB_LIB5 << endl;
                    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " [5d.02] label_AFLOWLIB_LIB6=" << label_AFLOWLIB_LIB6 << endl;
                    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " [5d.02] label_AFLOWLIB_LIB7=" << label_AFLOWLIB_LIB7 << endl;
                    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " [5d.02] label_AFLOWLIB_LIB8=" << label_AFLOWLIB_LIB8 << endl;
                    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " [5d.02] label_AFLOWLIB_LIB9=" << label_AFLOWLIB_LIB9 << endl;
                }
                if(LDEBUG) cerr << "DEBUG - " << soliloquy << " [5d.1]" << endl;
                //     if(xvasp.AVASP_libbase.size()==0) label_MIX=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB2)+"/LIB"+label_MIX;
                if(xvasp.AVASP_libbase.size()==0 && xvasp.AVASP_prototype_mode==LIBRARY_MODE_HTQC) label_MIX=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB2)+"/RAW"+label_MIX;
                if(xvasp.AVASP_libbase.size()==0 && xvasp.AVASP_prototype_mode==LIBRARY_MODE_HTQC_ICSD) label_MIX=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB2)+"/RAW"+label_MIX;
                if(xvasp.AVASP_libbase.size()==0 && xvasp.AVASP_prototype_mode==LIBRARY_MODE_LIB3) label_MIX=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB3)+"/LIB"+label_MIX;
                //  cerr << "LIB2=" << label_AFLOWLIB_LIB2 << endl;
                //  cerr << "LIB3=" << label_AFLOWLIB_LIB3 << endl;
            } else {
                if(LDEBUG) cerr << "DEBUG - " << soliloquy << " [5d.03] no check for MIX in non aflowlib servers " << endl;
            }	
            //    exit(0);

            if(LDEBUG) cerr << "DEBUG - AVASP_MakeSingleAFLOWIN [label_MIX] " << label_MIX << endl;
            if(xvasp.AVASP_aflowin_only_if_missing==TRUE) if(aurostd::FileExist(label_MIX)) {if(LDEBUG) cerr << label_MIX << endl; return TRUE;}  // ALREADY WRITTEN
            if(xvasp.AVASP_aflowin_only_if_missing==TRUE) if(aurostd::FileExist(label_AFLOWLIB_LIB1)) {if(LDEBUG) cerr << label_AFLOWLIB_LIB1 << endl; return TRUE;}  // IN AFLOWLIB
            if(xvasp.AVASP_aflowin_only_if_missing==TRUE) if(aurostd::FileExist(label_AFLOWLIB_LIB2)) {if(LDEBUG) cerr << label_AFLOWLIB_LIB2 << endl; return TRUE;}  // IN AFLOWLIB
            if(xvasp.AVASP_aflowin_only_if_missing==TRUE) if(aurostd::FileExist(label_AFLOWLIB_LIB3)) {if(LDEBUG) cerr << label_AFLOWLIB_LIB3 << endl; return TRUE;}  // IN AFLOWLIB
            if(xvasp.AVASP_aflowin_only_if_missing==TRUE) if(aurostd::FileExist(label_AFLOWLIB_LIB4)) {if(LDEBUG) cerr << label_AFLOWLIB_LIB4 << endl; return TRUE;}  // IN AFLOWLIB
            if(xvasp.AVASP_aflowin_only_if_missing==TRUE) if(aurostd::FileExist(label_AFLOWLIB_LIB5)) {if(LDEBUG) cerr << label_AFLOWLIB_LIB5 << endl; return TRUE;}  // IN AFLOWLIB
            if(xvasp.AVASP_aflowin_only_if_missing==TRUE) if(aurostd::FileExist(label_AFLOWLIB_LIB6)) {if(LDEBUG) cerr << label_AFLOWLIB_LIB6 << endl; return TRUE;}  // IN AFLOWLIB
            if(xvasp.AVASP_aflowin_only_if_missing==TRUE) if(aurostd::FileExist(label_AFLOWLIB_LIB7)) {if(LDEBUG) cerr << label_AFLOWLIB_LIB7 << endl; return TRUE;}  // IN AFLOWLIB
            if(xvasp.AVASP_aflowin_only_if_missing==TRUE) if(aurostd::FileExist(label_AFLOWLIB_LIB8)) {if(LDEBUG) cerr << label_AFLOWLIB_LIB8 << endl; return TRUE;}  // IN AFLOWLIB
            if(xvasp.AVASP_aflowin_only_if_missing==TRUE) if(aurostd::FileExist(label_AFLOWLIB_LIB9)) {if(LDEBUG) cerr << label_AFLOWLIB_LIB9 << endl; return TRUE;}  // IN AFLOWLIB

            if(xvasp.AVASP_aflowin_only_if_missing==TRUE) if(aurostd::FileExist("/common/BANDS/LIB/"+aurostd::vectorstring2string(xvasp.str.species_pp)+"/"+xvasp.AVASP_label+"/"+_AFLOWIN_)) {if(LDEBUG) cerr << "BANDS" << endl; return TRUE;}
            if(xvasp.AVASP_aflowin_only_if_missing==TRUE) if(aurostd::FileExist("/common/BANDS/RAW/"+aurostd::vectorstring2string(xvasp.str.species_pp)+"/"+xvasp.AVASP_label+"/"+_AFLOWIN_)) {if(LDEBUG) cerr << "BANDS" << endl; return TRUE;}
            // label_NOMIX - NOMIX
            label_NOMIX=label_MIX;aurostd::StringSubst(label_NOMIX,"LIB2","LIB2/NOMIX");
            if(LDEBUG) cerr << "DEBUG - AVASP_MakeSingleAFLOWIN [label_NOMIX] " << label_NOMIX << endl;
            if(xvasp.AVASP_aflowin_only_if_missing==TRUE) if(aurostd::FileExist(label_NOMIX)) {if(LDEBUG) cerr << label_NOMIX << endl; return TRUE;}
            // label_SKIP - SKIP
            label_SKIP=label_MIX;aurostd::StringSubst(label_SKIP,"LIB2","LIB2/SKIP");
            if(LDEBUG) cerr << "DEBUG - AVASP_MakeSingleAFLOWIN [label_SKIP] " << label_SKIP << endl;
            if(xvasp.AVASP_aflowin_only_if_missing==TRUE) if(aurostd::FileExist(label_SKIP)) {if(LDEBUG) cerr << label_SKIP << endl; return TRUE;}
            // label_ALREADY - ALREADY
            label_ALREADY=directory+"/"+_AFLOWIN_;
            if(LDEBUG) cerr << "DEBUG - AVASP_MakeSingleAFLOWIN [label_ALREADY] " << label_ALREADY << endl;
            if(xvasp.AVASP_aflowin_only_if_missing==TRUE) if(aurostd::FileExist(label_ALREADY)) {if(LDEBUG) cerr << label_ALREADY << endl; return TRUE;}
            // done generate
            //     cerr << "[1]" << xvasp.str.species_pp.at(0) << " " << xvasp.str.species_pp.at(1) << endl;
            //      cerr << xvasp.str.species.size() << endl; exit(0);
            if(LDEBUG) cerr << "DEBUG - " << soliloquy << " [6] xvasp.str.species.size()=" << xvasp.str.species.size() << endl;

            if(xvasp.str.species.size()==1) { // because specie 1 is part of the binary stuff
                if(LDEBUG) cerr << "DEBUG - " << soliloquy << " [6a] xvasp.str.species.size()=" << xvasp.str.species.size() << endl;	
                // DX [OBSOLETE] - Not sure why we add 2 of each; breaks ANRL unaries: deque<string> species_pp;      species_pp.push_back(xvasp.str.species_pp.at(0));species_pp.push_back(xvasp.str.species_pp.at(0));
                // DX [OBSOLETE] - Not sure why we add 2 of each; breaks ANRL unaries: deque<string> species_pp_backup; species_pp_backup.push_back(xvasp.str.species_pp.at(0));species_pp_backup.push_back(xvasp.str.species_pp.at(0));  // to prevent mess up of species
                // DX [OBSOLETE] - Not sure why we add 2 of each; breaks ANRL unaries: deque<double> species_volume;  species_volume.push_back(xvasp.str.species_volume.at(0));species_volume.push_back(xvasp.str.species_volume.at(0));
                // DX [OBSOLETE] - Not sure why we add 2 of each; breaks ANRL unaries: deque<double> species_mass;    species_mass.push_back(xvasp.str.species_mass.at(0));species_mass.push_back(xvasp.str.species_mass.at(0));
                deque<string> species_pp;      species_pp.push_back(xvasp.str.species_pp.at(0)); // DX 1/18/18 - instead of adding 2 (see above)
                deque<string> species_pp_backup; species_pp_backup.push_back(xvasp.str.species_pp.at(0)); // DX 1/18/18 - instead of adding 2 (see above)
                deque<double> species_volume;  species_volume.push_back(xvasp.str.species_volume.at(0)); // DX 1/18/18 - instead of adding 2 (see above)
                deque<double> species_mass;    species_mass.push_back(xvasp.str.species_mass.at(0)); // DX 1/18/18 - instead of adding 2 (see above)
                //	xvasp.str=aflowlib::PrototypeLibraries(oaus,xvasp.AVASP_label,xvasp.AVASP_parameters,species_pp,species_volume,xvasp.AVASP_volume_in); // FOR HTQC
                xvasp.str=aflowlib::PrototypeLibraries(oaus,xvasp.AVASP_label,xvasp.AVASP_parameters,species_pp,species_volume,xvasp.AVASP_volume_in,xvasp_in.AVASP_prototype_mode); 
                for(uint i=0;i<xvasp.str.species_pp.size()&&i<species_pp_backup.size();i++) xvasp.str.species_pp.at(i)=species_pp_backup.at(i); // to restore mess up of species
                // xvasp.str=PrototypeBinary(oaus,xvasp.AVASP_label,xvasp.str.species_pp.at(0),xvasp.str.species_volume.at(0),xvasp.str.species_pp.at(0),xvasp.str.species_volume.at(0),xvasp.AVASP_volume_in); // OLD WAY
                xvasp.str.species_mass=vmass;
            }

            if(xvasp.str.species.size()>=2) {
                if(LDEBUG) cerr << "DEBUG - " << soliloquy << " [7]  xvasp.str.species.size()>=2" << endl;
                for(uint i=0;i<xvasp.str.species_pp.size();i++) if(LDEBUG) cerr << "DEBUG - " << soliloquy << " [7.0a]: xvasp.str.species_pp.at(i) " << xvasp.str.species_pp.at(i) << endl;
                for(uint i=0;i<xvasp.str.species_pp_type.size();i++) if(LDEBUG) cerr << "DEBUG - " << soliloquy << " [7.0b]: xvasp.str.species_pp_type.at(i) " << xvasp.str.species_pp_type.at(i) << endl;
                for(uint i=0;i<xvasp.str.species_pp_version.size();i++) if(LDEBUG) cerr << "DEBUG - " << soliloquy << " [7.0c]: xvasp.str.species_pp_version.at(i) " << xvasp.str.species_pp_version.at(i) << endl;
                for(uint i=0;i<xvasp.str.species_pp_ZVAL.size();i++) if(LDEBUG) cerr << "DEBUG - " << soliloquy << " [7.0c]: xvasp.str.species_pp_ZVAL.at(i) " << xvasp.str.species_pp_ZVAL.at(i) << endl;
                for(uint i=0;i<xvasp.str.species_pp_vLDAU.size();i++) if(LDEBUG) cerr << "DEBUG - " << soliloquy << " [7.0c]: xvasp.str.species_pp_vLDAU.at(i) " << xvasp.str.species_pp_vLDAU.at(i).size() << endl;
                for(uint i=0;i<xvasp.str.species_volume.size();i++) if(LDEBUG) cerr << "DEBUG - " << soliloquy << " [7.0d]: xvasp.str.species_volume.at(i) " << xvasp.str.species_volume.at(i) << endl;
                for(uint i=0;i<xvasp.str.species_mass.size();i++) if(LDEBUG) cerr << "DEBUG - " << soliloquy << " [7.0e]: xvasp.str.species_mass.at(i) " << xvasp.str.species_mass.at(i) << endl;
                if(LDEBUG) cerr << "DEBUG - " << soliloquy << " [7.1a1]: xvasp.AVASP_volume_in=" << xvasp.AVASP_volume_in << endl;
                if(LDEBUG) cerr << "DEBUG - " << soliloquy << " [7.1a2]: xvasp.AVASP_label=" << xvasp.AVASP_label << endl;
                if(LDEBUG) {cerr << "DEBUG - " << soliloquy << " [7.1a3]: xvasp.str.species_pp.size()=" << xvasp.str.species_pp.size() << ": "; for(uint i=0;i<xvasp.str.species_pp.size();i++) cerr << " " << xvasp.str.species_pp.at(i); cerr << endl;}
                if(LDEBUG) {cerr << "DEBUG - " << soliloquy << " [7.1a3]: xvasp.str.species_volume.size()=" << xvasp.str.species_volume.size() << ": "; for(uint i=0;i<xvasp.str.species_volume.size();i++) cerr << " " << xvasp.str.species_volume.at(i); cerr << endl;}
                xstructure a;
                //	xvasp.str=aflowlib::PrototypeLibraries(oaus,xvasp.AVASP_label,xvasp.AVASP_parameters,xvasp.str.species_pp,xvasp.str.species_volume,xvasp.AVASP_volume_in);   // FOR HTQC
                xvasp.str=aflowlib::PrototypeLibraries(oaus,xvasp.AVASP_label,xvasp.AVASP_parameters,xvasp.str.species_pp,xvasp.str.species_volume,xvasp.AVASP_volume_in,xvasp_in.AVASP_prototype_mode); 

                if(LDEBUG) cerr << "DEBUG - " << soliloquy << " [7.1b]  xvasp.str.species.size()>=2" << endl;
                if(LDEBUG) cerr << xvasp.str << endl;
                if(LDEBUG) cerr << "DEBUG - " << soliloquy << " [7.1c]  xvasp.str.species.size()>=2" << endl;
                // xvasp.str=PrototypeBinary(oaus,xvasp.AVASP_label,xvasp.str.species_pp.at(0),xvasp.str.species_volume.at(0),xvasp.str.species_pp.at(1),xvasp.str.species_volume.at(1),xvasp.AVASP_volume_in); // OLD WAY
                xvasp.str.species_mass=vmass;
            }
            if(LDEBUG) cerr << "DEBUG - " << soliloquy << " [8] xvasp.str.species.size()=" << xvasp.str.species.size() << endl;
            //   cerr << "[7]" << xvasp.str.species_pp.at(0) << " " << xvasp.str.species_pp.at(1) << endl;
        } else {
            directory=xvasp.AVASP_dirbase+"/"+aurostd::vectorstring2string(xvasp.str.species_pp)+"/"+xvasp.AVASP_label;
            // cerr << "DEBUG - " << soliloquy << " directory(2)=" << directory << endl;
            system=aurostd::vectorstring2string(xvasp.str.species_pp)+"."+xvasp.AVASP_label;
            label_MIX=xvasp.AVASP_libbase+"/"+aurostd::vectorstring2string(xvasp.str.species_pp)+"/"+xvasp.AVASP_label+"/"+_AFLOWIN_;
        }
        if(xvasp.AVASP_directory_from_library_==FALSE) {
            directory=xvasp.AVASP_dirbase+"/"+xvasp.AVASP_label;
            //   cerr << "DEBUG - " << soliloquy << " directory(3)=" << directory << endl;
            system=xvasp.AVASP_label;
        }
        if(LDEBUG) cerr << "DEBUG - " << soliloquy << " [9]" << endl;

        // MAKE Directory
        //  Krun=aurostd::DirectoryMake(directory); // do later
        // if(Krun==FALSE) return Krun;             // do later
        // UNDERSTAND SPECIES !
        //  cerr << str.num_each_type.size() << endl;
        //  cerr << str.comp_each_type.size() << endl;
        if(xvasp.str.SpeciesGetAlphabetic()==FALSE) xvasp.str.SpeciesPutAlphabetic();
    }

    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " " << "[10]" << endl;

    // cerr << "xvasp.str.species_pp.size()=" << xvasp.str.species_pp.size() << endl; for(uint i=0;i<xvasp.str.species_pp.size();i++) cerr << xvasp.str.species_pp.at(i) << endl; exit(0);

    // for ICSD
    if(xvasp.AVASP_prototype_mode==LIBRARY_MODE_ICSD) {
        if(LDEBUG) cerr << "DEBUG: xvasp.AVASP_prototype_mode==LIBRARY_MODE_ICSD" << endl;
        // xvasp.AVASP_dirbase=_ICSD_DIRBASE_;
        // xvasp.AVASP_libbase=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_ICSD);
        // directory=xvasp.AVASP_dirbase+"/"+xvasp.AVASP_label;
        system=xvasp.AVASP_label;
        //    if(xvasp.str.spacegroupnumber!=0) {
        //       cerr << "DEBUG: calling LATTICE::Lattice_SpaceGroup from AVASP_MakeSingleAFLOWIN" << endl;
        //       directory=xvasp.AVASP_dirbase+"/"+LATTICE::Lattice_SpaceGroup(xvasp.str.spacegroupnumber)+"/"+xvasp.AVASP_label;}
        //     else
        directory=xvasp.AVASP_dirbase+"/"+xvasp.str.bravais_lattice_type+"/"+xvasp.AVASP_label;
        //  cerr << "DEBUG - " << soliloquy << " directory(4)=" << directory << endl;
        //    cerr << directory << endl;
    }

    // for ICSD
    if(xvasp.AVASP_prototype_mode==LIBRARY_MODE_PROTOTYPE || xvasp.AVASP_prototype_mode==LIBRARY_MODE_XSTRUCTURE) {
        system="";
        formula="";
        directory=xvasp.Directory;
        for(uint i=0;i<xvasp.str.species_pp.size();i++) {
            system+=xvasp.str.species_pp.at(i);
            formula+=xvasp.str.species_pp.at(i)+aurostd::utype2string(xvasp.str.num_each_type.at(i));
        }
        //if(xvasp.AVASP_prototype_mode==LIBRARY_MODE_PROTOTYPE) directory=xvasp.Directory;
        if(xvasp.AVASP_prototype_mode==LIBRARY_MODE_XSTRUCTURE) system+="/"+xvasp.Directory;
        aurostd::StringSubst(system,"/./","/");aurostd::StringSubst(system,"//","/");
    }

    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " " << "[11]" << endl;

    bool AFLOWIN_OLD=FALSE;
    bool AFLOWIN_NEW=TRUE;

    // DX 1/23/18 - Add ANRL parameters to label - START
    // check for ANRL
    if(xvasp.aopts.flag("AFLOWIN_FLAG::PARAMS")==TRUE) {
        string anrl_add_on=":ANRL="+xvasp.AVASP_parameters;
        system+=anrl_add_on;directory+=anrl_add_on;xvasp.AVASP_label+=anrl_add_on;
    }
    // DX 1/23/18 - Add ANRL parameters to label - END

    // check for PSTRESS/PRESSURE   //  if(xvasp.aopts.flag("AFLOWIN_FLAG::PSTRESS",TRUE)) xvasp.aopts.push_attached("AFLOWIN_FLAG::PSTRESS","7.65");
    if(xvasp.aopts.flag("AFLOWIN_FLAG::PSTRESS")==TRUE) {
        double PSTRESS=aurostd::string2utype<double>(xvasp.aopts.getattachedscheme("AFLOWIN_FLAG::PSTRESS"));
        if(abs(PSTRESS)>0.00001) {
            if(LDEBUG) cerr << "DEBUG - " << soliloquy << " with pressure " << "[12p]" << endl;
            string pressure_add_on=":P="+aurostd::utype2string(PSTRESS,_AVASP_DOUBLE2STRING_PRECISION_)+"kB";
            if(LDEBUG) cerr << "DEBUG - " << soliloquy << " with pressure " << "[12p2] pressure_add_on="  << pressure_add_on << endl;
            system+=pressure_add_on;directory+=pressure_add_on;xvasp.AVASP_label+=pressure_add_on;
        }
    }

    // check for LDAU1
    if(0 && xvasp.aopts.flag("FLAG::AVASP_LDAU1") && !aurostd::substring2bool(system,"_ICSD")) {
        if(LDEBUG) cerr << "DEBUG - " << soliloquy << " with LDAU1 " << "[12p3]" << endl;
        string ldau1_add_on=":LDAU1";
        system+=ldau1_add_on;directory+=ldau1_add_on;xvasp.AVASP_label+=ldau1_add_on;
    }
    // check for LDAU2
    if((xvasp.aopts.flag("FLAG::AVASP_LDAU1") || xvasp.aopts.flag("FLAG::AVASP_LDAU2") ) && !aurostd::substring2bool(system,"_ICSD")) {
        if(LDEBUG) cerr << "DEBUG - " << soliloquy << " with LDAU2 " << "[12p4]" << endl;
        string ldau2_add_on=":LDAU2";
        system+=ldau2_add_on;directory+=ldau2_add_on;xvasp.AVASP_label+=ldau2_add_on;
    }

    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " " << "[12.2]" << endl;

    // cerr << "[9] xvasp.AVASP_potential=" << xvasp.AVASP_potential << endl; exit(0);

    // skipped if necessary
    // PREAMBLE
    aflowin << AFLOWIN_SEPARATION_LINE << endl; // [AFLOW] **************************************************
    aflowin << "[AFLOW]                                                                                     " << endl;
    aflowin << "[AFLOW]                     .o.        .o88o. oooo                                          " << endl;
    aflowin << "[AFLOW]                    .888.       888 `` `888                                          " << endl;
    aflowin << "[AFLOW]                   .8'888.     o888oo   888   .ooooo.  oooo oooo    ooo              " << endl;
    aflowin << "[AFLOW]                  .8' `888.     888     888  d88' `88b  `88. `88.  .8'               " << endl;
    aflowin << "[AFLOW]                 .88ooo8888.    888     888  888   888   `88..]88..8'                " << endl;
    aflowin << "[AFLOW]                .8'     `888.   888     888  888   888    `888'`888'                 " << endl;
    aflowin << "[AFLOW]               o88o     o8888o o888o   o888o `Y8bod8P'     `8'  `8'  .in             " << endl;
    aflowin << "[AFLOW]                                                                                     " << endl;
    aflowin << AFLOWIN_SEPARATION_LINE << endl; // [AFLOW]  **************************************************
    aflowin << "[AFLOW] * Stefano Curtarolo - (AFLOW V" << string(AFLOW_VERSION) << ") " << endl;
    aflowin << "[AFLOW] * D. Morgan, W. Setyawan, G. Hart, M. Jahnatek, S. Wang, O. Levy, K. Yang, J. Xue,  " << endl;
    aflowin << "[AFLOW] * R. Taylor, C. Calderon, C. Toher, C. Oses, J. J. Plata, D. Hicks, P. Nath, F. Rose,  " << endl;
    aflowin << "[AFLOW] * E. Gossett, E. Perim, R. Chepulskyy, K. Rasch,  M. Fornari, M. Buongiorno Nardelli " << endl;
    aflowin << AFLOWIN_SEPARATION_LINE << endl; // [AFLOW] **************************************************
    aflowin << "[AFLOW] Aflow automatically generated (aflow_avasp.cpp) " << endl;
    aflowin << "[AFLOW] GENERATOR = " << XHOST.User << " " << endl;  //CO 180622
    aflowin << AFLOWIN_SEPARATION_LINE << endl; // [AFLOW] **************************************************
    aflowin << AFLOWIN_SEPARATION_LINE << endl; // [AFLOW] **************************************************
    aflowin << "[AFLOW]SYSTEM=" << system << endl;
    if(xvasp.AVASP_prototype_mode==LIBRARY_MODE_PROTOTYPE || xvasp.AVASP_prototype_mode==LIBRARY_MODE_XSTRUCTURE)
        aflowin << "[AFLOW]FORMULA=" << formula << endl;
    if(xvasp.str.species.size()==1)
        aflowin << "#[AFLOW] single element calculation" << endl;
    aflowin << AFLOWIN_SEPARATION_LINE << endl; // [AFLOW] **************************************************
    aflowin << "[AFLOW] input file for aflow " << endl;
    aflowin << "[AFLOW_MODE=VASP] " << endl;
    aflowin << AFLOWIN_SEPARATION_LINE << endl; // [AFLOW] **************************************************
    aflowin << "[AFLOW_MODE_ZIP="+DEFAULT_KZIP_BIN+"] " << endl;
    aflowin << "[AFLOW_MODE_BINARY=" << DEFAULT_VASP_BIN << "] " << endl;
    aflowin << AFLOWIN_SEPARATION_LINE << endl; // [AFLOW] **************************************************
    aflowin << AFLOWIN_SEPARATION_LINE << endl; // [AFLOW] **************************************************
    if(xvasp.AVASP_flag_MPI || XHOST.MPI) {
        aflowin << "[AFLOW_MODE_MPI]" << endl;
    } else {
        aflowin << "#[AFLOW_MODE_MPI]" << endl;
    }
    string NCPUS_VAL="MAX";
    if(XHOST.vflag_control.flag("XPLUG_NUM_THREADS")) NCPUS_VAL=XHOST.vflag_control.getattachedscheme("XPLUG_NUM_THREADS");
    aflowin << "[AFLOW_MODE_MPI_MODE]NCPUS=" << NCPUS_VAL << " " << endl;
#ifdef MPI_LAM
    aflowin << "[AFLOW_MODE_MPI_MODE]START=\"lamboot\" " << endl;
    aflowin << "[AFLOW_MODE_MPI_MODE]STOP=\"lamhalt\" " << endl;
#endif
    aflowin << "[AFLOW_MODE_MPI_MODE]COMMAND =\"" << MPI_COMMAND_DEFAULT << "\" " << endl;
    aflowin << "[AFLOW_MODE_MPI_MODE]AUTOTUNE " << endl;
    aflowin << "[AFLOW_MODE_MPI_MODE]BINARY=\"" << DEFAULT_VASP_MPI_BIN << "\" " << endl;
    aflowin << AFLOWIN_SEPARATION_LINE << endl; // [AFLOW] **************************************************
    // SYMMETRY WRITE
    if(!xvasp.aopts.flag("FLAGS::AVASP_SYMMMETRY=OFF")) {
        aflowin << "[AFLOW_SYMMETRY]CALC " << endl;
        aflowin << "#[AFLOW_SYMMETRY]SGROUP_WRITE " << endl;
        aflowin << "#[AFLOW_SYMMETRY]SGROUP_RADIUS=7.77 " << endl;
        aflowin << AFLOWIN_SEPARATION_LINE << endl; // [AFLOW] **************************************************
    }
    // NEIGHBOURS WRITE
    if(!xvasp.aopts.flag("FLAGS::AVASP_NEIGHBOURS=OFF")) {
        aflowin << "#[AFLOW_NEIGHBOURS]CALC " << endl;
        aflowin << "[AFLOW_NEIGHBOURS]RADIUS=7.7 " << endl;
        aflowin << "[AFLOW_NEIGHBOURS]DRADIUS=0.1 " << endl;
        aflowin << AFLOWIN_SEPARATION_LINE << endl; // [AFLOW] **************************************************
    }
    string MODULE = aurostd::toupper(xvasp.aopts.getattachedscheme("AFLOWIN_FLAG::MODULE"));
    // APL WRITING
    if(!xvasp.aopts.flag("FLAGS::AVASP_APL=OFF")) { // CO 180214 - I interpret this flag to refer to WRITING APL options, not if they are on
        aflowin << aurostd::PaddedPOST((MODULE=="APL"?string(""):string("#"))+"[AFLOW_APL]CALC",_aflowinpad_) << "// README_AFLOW_APL.TXT" << endl; // CO 180214
        aflowin << aurostd::PaddedPOST("[AFLOW_APL]ENGINE=DM",_aflowinpad_) << "// README_AFLOW_APL.TXT" << endl;
        aflowin << aurostd::PaddedPOST("[AFLOW_APL]DMAG=0.015",_aflowinpad_) << "// README_AFLOW_APL.TXT" << endl;
        // DX and CO - START
        aflowin << aurostd::PaddedPOST((xvasp.aopts.flag("AFLOWIN_FLAG::APL_SUPERCELL")?string("#"):string(""))+"[AFLOW_APL]MINATOMS=100",_aflowinpad_) << "// README_AFLOW_APL.TXT" << endl;
        aflowin << aurostd::PaddedPOST((xvasp.aopts.flag("AFLOWIN_FLAG::APL_SUPERCELL")?string(""):string("#"))+"[AFLOW_APL]SUPERCELL="+
                (xvasp.aopts.flag("AFLOWIN_FLAG::APL_SUPERCELL")?xvasp.aopts.getattachedscheme("AFLOWIN_FLAG::APL_SUPERCELL"):"3x3x3"),_aflowinpad_) << "// README_AFLOW_APL.TXT" << endl;  // CO 180214
        aflowin << aurostd::PaddedPOST("[AFLOW_APL]DC=y",_aflowinpad_) << "// README_AFLOW_APL.TXT" << endl;
        aflowin << aurostd::PaddedPOST("[AFLOW_APL]DPM=y",_aflowinpad_) << "// README_AFLOW_APL.TXT" << endl;
        aflowin << aurostd::PaddedPOST("[AFLOW_APL]ZEROSTATE=y",_aflowinpad_) << "// README_AFLOW_APL.TXT" << endl;
        // DX AND CO - END
        aflowin << aurostd::PaddedPOST("[AFLOW_APL]DOS=y",_aflowinpad_) << "// README_AFLOW_APL.TXT" << endl;
        aflowin << aurostd::PaddedPOST("[AFLOW_APL]TP=y",_aflowinpad_) << "// README_AFLOW_APL.TXT" << endl;
        aflowin << aurostd::PaddedPOST("[AFLOW_APL]TPT=0:2000:10",_aflowinpad_) << "// README_AFLOW_APL.TXT" << endl;
        aflowin << AFLOWIN_SEPARATION_LINE << endl; // [AFLOW] **************************************************
    }
    if(!xvasp.aopts.flag("FLAGS::AVASP_QHA=OFF")) {
        aflowin << aurostd::PaddedPOST((MODULE=="QHA"?string(""):string("#"))+"[AFLOW_QHA]CALC",_aflowinpad_) << "// README_AFLOW_QHA_SCQHA_QHA3P.TXT" << endl;
        aflowin << aurostd::PaddedPOST("[AFLOW_QHA]MODE=QHA3P",_aflowinpad_) << "// README_AFLOW_QHA_SCQHA_QHA3P.TXT" << endl;
        aflowin << aurostd::PaddedPOST("[AFLOW_QHA]EOS=y",_aflowinpad_) << "// README_AFLOW_QHA_SCQHA_QHA3P.TXT" << endl;
        //[OBSOLETE PN180717]aflowin << aurostd::PaddedPOST("[AFLOW_QHA]SCQHA=y",_aflowinpad_) << "// README_AFLOW_QHA_SCQHA_QHA3P.TXT" << endl;
        //[OBSOLETE PN180717]aflowin << aurostd::PaddedPOST("[AFLOW_QHA]GP_DISTORTION=0.03",_aflowinpad_) << "// README_AFLOW_QHA_SCQHA_QHA3P.TXT" << endl;
        //[OBSOLETE PN180717]aflowin << aurostd::PaddedPOST("[AFLOW_QHA]SCQHA_DISTORTION=3",_aflowinpad_) << "// README_AFLOW_QHA_SCQHA_QHA3P.TXT" << endl;
        //[OBSOLETE PN180717]aflowin << aurostd::PaddedPOST("[AFLOW_QHA]EOS_DISTORTION_RANGE=-3:6:1",_aflowinpad_) << "// README_AFLOW_QHA_SCQHA_QHA3P.TXT" << endl;
        //[OBSOLETE PN180717]aflowin << aurostd::PaddedPOST("[AFLOW_QHA]SCQHA_PDIS_T=50,100,200",_aflowinpad_) << "// README_AFLOW_QHA_SCQHA_QHA3P.TXT" << endl;
        //[OBSOLETE CO 180705]aflowin << aurostd::PaddedPOST("[AFLOW_QHA]GP_VOL_DISTORTION_PERCENTAGE=0.03",_aflowinpad_) << "// README_AFLOW_APL.TXT" << endl;
        //[OBSOLETE CO 180705]aflowin << aurostd::PaddedPOST("[AFLOW_QHA]DISPLACEMENTS=y",_aflowinpad_) << "// README_AFLOW_APL.TXT" << endl;
        //[OBSOLETE CO 180705]aflowin << aurostd::PaddedPOST("[AFLOW_QHA]PROJECTION_DIR=1:1:1",_aflowinpad_) << "// README_AFLOW_APL.TXT" << endl;
        //[OBSOLETE CO 180705]aflowin << aurostd::PaddedPOST("[AFLOW_QHA]EOS=n",_aflowinpad_) << "// README_AFLOW_APL.TXT" << endl;
        //[OBSOLETE CO 180705]aflowin << aurostd::PaddedPOST("[AFLOW_QHA]EOS_VOLRANGE_DIST=-2:4:0.5",_aflowinpad_) << "// README_AFLOW_APL.TXT" << endl;
        //[OBSOLETE CO 180705]aflowin << aurostd::PaddedPOST("[AFLOW_QHA]EOS_KPOINTS_MODE=32768:10000:20:100000",_aflowinpad_) << "// README_AFLOW_APL.TXT" << endl;
        aflowin << AFLOWIN_SEPARATION_LINE << endl; // [AFLOW] **************************************************
    }
    if(!xvasp.aopts.flag("FLAGS::AVASP_AAPL=OFF")) {
        aflowin << aurostd::PaddedPOST((MODULE=="AAPL"?string(""):string("#"))+"[AFLOW_AAPL]CALC",_aflowinpad_) << "// README_AFLOW_APL.TXT" << endl;
        aflowin << aurostd::PaddedPOST("[AFLOW_AAPL]TDMAG=0.015",_aflowinpad_) << "// README_AFLOW_APL.TXT" << endl;
        aflowin << aurostd::PaddedPOST("[AFLOW_AAPL]CUT_SHELL=4",_aflowinpad_) << "// README_AFLOW_APL.TXT" << endl;
        aflowin << aurostd::PaddedPOST("[AFLOW_AAPL]CUT_RAD=4.5",_aflowinpad_) << "// README_AFLOW_APL.TXT" << endl;
        aflowin << aurostd::PaddedPOST("[AFLOW_AAPL]SUMRULE=1E-5",_aflowinpad_) << "// README_AFLOW_APL.TXT" << endl;
        aflowin << aurostd::PaddedPOST("[AFLOW_AAPL]BTE=FULL",_aflowinpad_) << "// README_AFLOW_APL.TXT" << endl;
        aflowin << aurostd::PaddedPOST("[AFLOW_AAPL]THERMALGRID=21x21x21",_aflowinpad_) << "// README_AFLOW_APL.TXT" << endl;
        aflowin << aurostd::PaddedPOST("[AFLOW_AAPL]ISOTOPE=y",_aflowinpad_) << "// README_AFLOW_APL.TXT" << endl;
        aflowin << aurostd::PaddedPOST("[AFLOW_AAPL]CUMULATIVEK=y",_aflowinpad_) << "// README_AFLOW_APL.TXT" << endl;
        aflowin << aurostd::PaddedPOST("[AFLOW_AAPL]BOUNDARY=n",_aflowinpad_) << "// README_AFLOW_APL.TXT" << endl;
        aflowin << aurostd::PaddedPOST("[AFLOW_AAPL]NANO_SIZE=100",_aflowinpad_) << "// README_AFLOW_APL.TXT" << endl;
        aflowin << aurostd::PaddedPOST("[AFLOW_AAPL]TCT=200:700:20",_aflowinpad_) << "// README_AFLOW_APL.TXT" << endl;
        aflowin << AFLOWIN_SEPARATION_LINE << endl; // [AFLOW] **************************************************
    }
    if(0) {
        if(xvasp.str.species.size()==1 || xvasp.AVASP_prototype_mode==LIBRARY_MODE_HTQC || xvasp.AVASP_prototype_mode==LIBRARY_MODE_HTQC_ICSD) {
            xvasp.AVASP_flag_RUN_RELAX=FALSE;
            xvasp.AVASP_flag_RUN_RELAX_STATIC_BANDS=TRUE;
        }
    }

    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " " << "[12.3]" << endl;

    // RELAX WRITING
    if(xvasp.AVASP_flag_RUN_RELAX) aflowin << aurostd::PaddedPOST("[VASP_RUN]RELAX=2",_aflowinpad_);
    if(xvasp.AVASP_flag_RUN_RELAX_STATIC) aflowin << aurostd::PaddedPOST("[VASP_RUN]RELAX_STATIC=2",_aflowinpad_);
    if(xvasp.AVASP_flag_RUN_RELAX_STATIC_BANDS) aflowin << aurostd::PaddedPOST("[VASP_RUN]RELAX_STATIC_BANDS=2",_aflowinpad_);
    if(xvasp.AVASP_flag_RUN_STATIC_BANDS) aflowin << aurostd::PaddedPOST("[VASP_RUN]STATIC_BANDS",_aflowinpad_);
    if(xvasp.AVASP_flag_RUN_STATIC) aflowin << aurostd::PaddedPOST("[VASP_RUN]STATIC",_aflowinpad_);
    if(xvasp.AVASP_flag_GENERATE) aflowin << aurostd::PaddedPOST("[VASP_RUN]GENERATE",_aflowinpad_);
    aflowin << "// GENERATE | STATIC | RELAX=N | RELAX_STATIC=N | STATIC_BANDS | RELAX_STATIC_BANDS=N | REPEAT_BANDS [,DS[,DD[,DSCF]]] " << endl;

    if(xvasp.AVASP_flag_RUN_RELAX_STATIC ||
            xvasp.AVASP_flag_RUN_RELAX_STATIC_BANDS ||
            xvasp.AVASP_flag_RUN_STATIC_BANDS ||
            xvasp.AVASP_flag_RUN_STATIC) {
        if(aurostd::substring2bool(xvasp.AVASP_potential,"pot_LDA")) aflowin << "[VASP_FORCE_OPTION]RWIGS_STATIC" << endl; // dont need for paw
        if(aurostd::substring2bool(xvasp.AVASP_potential,"pot_GGA")) aflowin << "[VASP_FORCE_OPTION]RWIGS_STATIC" << endl; // dont need for paw
    }

    // KPOINTS WRITING
    if(xvasp.AVASP_value_NSW>0) {
        aflowin << "[VASP_FORCE_OPTION]NSW=101 " << endl;
    }

    // NOMIX WRITING
    if(xvasp.aopts.flag("FLAG::AVASP_SKIP_NOMIX")) {
        aflowin << "[VASP_FORCE_OPTION]NEGLECT_NOMIX " << endl;
    } else {
        aflowin << "#[VASP_FORCE_OPTION]NEGLECT_NOMIX " << endl;
    }

    // WAVECAR WRITING
    if(xvasp.aopts.flag("FLAG::AVASP_WAVECAR")) {
        aflowin << aurostd::PaddedPOST("[VASP_FORCE_OPTION]WAVECAR=ON",_aflowinpad_) << "// ON | OFF (default OFF)" << endl;
    } else {
        // aflowin << aurostd::PaddedPOST("[VASP_FORCE_OPTION]WAVECAR=OFF",_aflowinpad_) << "// ON | OFF (default OFF)" << endl;
    }

    // CHGCAR WRITING
    if(xvasp.aopts.flag("FLAG::AVASP_CHGCAR")) {
        // aflowin << aurostd::PaddedPOST("[VASP_FORCE_OPTION]CHGCAR=ON",_aflowinpad_) << "// ON | OFF (default ON)" << endl;
    } else {
        aflowin << aurostd::PaddedPOST("[VASP_FORCE_OPTION]CHGCAR=OFF",_aflowinpad_) << "// ON | OFF (default ON)" << endl;
    }

    // KPOINTS WRITING
    aflowin << aurostd::PaddedPOST("#[VASP_FORCE_OPTION]KPOINTS=keyword[,keyword]",_aflowinpad_) << "// EVEN | ODD | KSHIFT_GAMMA_EVEN | KSHIFT_GAMMA_ODD | KSCHEME_MONKHORST_PACK | KSCHEME_GAMMA | GAMMA | KEEPK | IBZKPT" << endl;

    aflowin << aurostd::PaddedPOST("[VASP_FORCE_OPTION]SYM=ON",_aflowinpad_) << "// ON | OFF  (default ON)" << endl;

    // cerr << "xvasp.AVASP_potential=" << xvasp.AVASP_potential << endl; exit(0);

    string aus_PP="";
    if(xvasp.AVASP_potential==DEFAULT_VASP_POTCAR_DIR_POT_LDA) aus_PP="pot_LDA";
    if(xvasp.AVASP_potential==DEFAULT_VASP_POTCAR_DIR_POT_GGA) aus_PP="pot_GGA";
    if(xvasp.AVASP_potential==DEFAULT_VASP_POTCAR_DIR_POT_PBE) aus_PP="pot_PBE";
    if(xvasp.AVASP_potential==DEFAULT_VASP_POTCAR_DIR_POTPAW_LDA) aus_PP="potpaw_LDA";
    if(xvasp.AVASP_potential==DEFAULT_VASP_POTCAR_DIR_POTPAW_GGA) aus_PP="potpaw_GGA";
    if(xvasp.AVASP_potential==DEFAULT_VASP_POTCAR_DIR_POTPAW_PBE) aus_PP="potpaw_PBE";
    if(xvasp.AVASP_potential==DEFAULT_VASP_POTCAR_DIR_POTPAW_LDA_KIN) aus_PP="potpaw_LDA_KIN";
    if(xvasp.AVASP_potential==DEFAULT_VASP_POTCAR_DIR_POTPAW_PBE_KIN) aus_PP="potpaw_PBE_KIN";

    if(xvasp.aopts.flag("FLAG::AVASP_AUTO_PSEUDOPOTENTIALS")) {
        aflowin << aurostd::PaddedPOST("[VASP_FORCE_OPTION]AUTO_PSEUDOPOTENTIALS="+aus_PP,_aflowinpad_) << "// pot_LDA | pot_GGA | potpaw_LDA | potpaw_GGA | potpaw_PBE | potpaw_LDA_KIN | potpaw_PBE_KIN  " << endl;

    } else {
        //   cerr << soliloquy << " AUTO_PSEUDOPOTENTIALS=FALSE" << endl;
        aflowin << aurostd::PaddedPOST("#[VASP_FORCE_OPTION]AUTO_PSEUDOPOTENTIALS="+aus_PP,_aflowinpad_) << "// pot_LDA | pot_GGA | potpaw_LDA | potpaw_GGA | potpaw_PBE | potpaw_LDA_KIN | potpaw_PBE_KIN  " << endl;
    }

    // NBANDS WRITING
    if(xvasp.aopts.flag("AFLOWIN_FLAG::NBANDS")) {
        aflowin << aurostd::PaddedPOST("[VASP_FORCE_OPTION]NBANDS",_aflowinpad_) << "// Estimate Bands (better than VASP)" << endl;
    } else {
        //  aflowin << aurostd::PaddedPOST("#[VASP_FORCE_OPTION]NBANDS",_aflowinpad_) << "// Estimate Bands (better than VASP)" << endl;
    }

    // PSTRESS WRITING
    if(xvasp.aopts.flag("AFLOWIN_FLAG::PSTRESS")) {
        double PSTRESS=aurostd::string2utype<double>(xvasp.aopts.getattachedscheme("AFLOWIN_FLAG::PSTRESS"));
        aflowin << aurostd::PaddedPOST("[VASP_FORCE_OPTION]PSTRESS="+aurostd::utype2string(PSTRESS,_AVASP_DOUBLE2STRING_PRECISION_),_aflowinpad_) << "// Pressure in kBar (1kB=0.1GPa)" << endl;
    } else {
        aflowin << aurostd::PaddedPOST("#[VASP_FORCE_OPTION]PSTRESS=0.0",_aflowinpad_) << "// Pressure in kBar (1kB=0.1GPa)" << endl;
    }

    // EDIFFG WRITING
    if(xvasp.aopts.flag("AFLOWIN_FLAG::EDIFFG")) {
        double EDIFFG=aurostd::string2utype<double>(xvasp.aopts.getattachedscheme("AFLOWIN_FLAG::EDIFFG"));
        aflowin << aurostd::PaddedPOST("[VASP_FORCE_OPTION]EDIFFG="+aurostd::utype2string(EDIFFG,_AVASP_DOUBLE2STRING_PRECISION_),_aflowinpad_) << "// EDIFFG for relaxed forces" << endl;
    } else {
        aflowin << aurostd::PaddedPOST("#[VASP_FORCE_OPTION]EDIFFG="+aurostd::utype2string(DEFAULT_VASP_PREC_EDIFFG,_AVASP_DOUBLE2STRING_PRECISION_),_aflowinpad_) << "// EDIFFG for relaxed forces" << endl;
    }

    // ENMAX_MULTIPLY WRITING
    if(xvasp.aopts.flag("AFLOWIN_FLAG::ENMAX_MULTIPLY")) {
        double ENMAX_MULTIPLY=aurostd::string2utype<double>(xvasp.aopts.getattachedscheme("AFLOWIN_FLAG::ENMAX_MULTIPLY"));
        aflowin << aurostd::PaddedPOST("[VASP_FORCE_OPTION]ENMAX_MULTIPLY="+aurostd::utype2string(ENMAX_MULTIPLY,_AVASP_DOUBLE2STRING_PRECISION_),_aflowinpad_) << "// extra multiplication" << endl;
    } else {
        aflowin << aurostd::PaddedPOST("#[VASP_FORCE_OPTION]ENMAX_MULTIPLY="+aurostd::utype2string(DEFAULT_VASP_PREC_ENMAX_HIGH,4),_aflowinpad_) << "// Multiplication of the max(pseudopotential_cutoffs)" << endl;
    }

    // POTIM WRITING
    if(xvasp.aopts.flag("AFLOWIN_FLAG::POTIM")) {
        double POTIM=aurostd::string2utype<double>(xvasp.aopts.getattachedscheme("AFLOWIN_FLAG::POTIM"));
        aflowin << aurostd::PaddedPOST("[VASP_FORCE_OPTION]POTIM="+aurostd::utype2string(POTIM,_AVASP_DOUBLE2STRING_PRECISION_),_aflowinpad_) << "// for ionic time-step" << endl;
    } else {
        aflowin << aurostd::PaddedPOST("#[VASP_FORCE_OPTION]POTIM="+aurostd::utype2string(DEFAULT_VASP_PREC_POTIM,_AVASP_DOUBLE2STRING_PRECISION_),_aflowinpad_) << "// ionic time-step" << endl;
    }

    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " " << "[12.4]" << endl;

    // SPIN WRITING
    if(AFLOWIN_OLD) { // OLD
        if(xvasp.aopts.flag("FLAG::AVASP_SPIN")) {
            aflowin << aurostd::PaddedPOST("[VASP_FORCE_OPTION]SPIN=ON",_aflowinpad_) << "// ON | OFF  (default NEGLECT)" << endl;
            aflowin << "[VASP_FORCE_OPTION]SPIN_REMOVE_RELAX_1 " << endl;
        } else {
            aflowin << aurostd::PaddedPOST("[VASP_FORCE_OPTION]SPIN=OFF",_aflowinpad_) << "// ON | OFF  (default NEGLECT)" << endl;
        }
    }
    if(AFLOWIN_NEW) { // NEW
        vector<string> vstr;string tstr;
        vstr.clear();tstr="";
        if(xvasp.aopts.flag("FLAG::AVASP_SPIN"))  vstr.push_back("ON");
        if(!xvasp.aopts.flag("FLAG::AVASP_SPIN")) vstr.push_back("OFF");
        if(xvasp.aopts.flag("FLAG::AVASP_SPIN") && xvasp.aopts.flag("FLAG::AVASP_SPIN_REMOVE_RELAX_1")) vstr.push_back("REMOVE_RELAX_1");
        if(xvasp.aopts.flag("FLAG::AVASP_SPIN") && xvasp.aopts.flag("FLAG::AVASP_SPIN_REMOVE_RELAX_2")) vstr.push_back("REMOVE_RELAX_2");
        if(vstr.size()>0) {
            for(uint i=0;i<vstr.size();i++) tstr+=vstr.at(i)+(i!=vstr.size()-1 ? "," : " ");
            aflowin << aurostd::PaddedPOST("[VASP_FORCE_OPTION]SPIN="+tstr,_aflowinpad_) << "// (ON | OFF  (default ON)), REMOVE_RELAX_1 | _2" << endl;
        }
    }

    // BADER WRITING
    if(xvasp.aopts.flag("FLAG::AVASP_BADER")) {
        aflowin << aurostd::PaddedPOST("[VASP_FORCE_OPTION]BADER=ON",_aflowinpad_) << "// ON | OFF (default OFF)" << endl;
    } else {
        //    aflowin << aurostd::PaddedPOST("[VASP_FORCE_OPTION]BADER=OFF",_aflowinpad_) << "// ON | OFF (default OFF)" << endl;
    }
    // ELF WRITING
    if(xvasp.aopts.flag("FLAG::AVASP_ELF")) {
        aflowin << aurostd::PaddedPOST("[VASP_FORCE_OPTION]ELF=ON",_aflowinpad_) << "// ON | OFF (default OFF)" << endl;
    } else {
        //    aflowin << aurostd::PaddedPOST("[VASP_FORCE_OPTION]ELF=OFF",_aflowinpad_) << "// ON | OFF (default OFF)" << endl;
    }

    // LSCOUPLING WRITING
    if(xvasp.aopts.flag("FLAG::AVASP_LSCOUPLING")) {
        aflowin << aurostd::PaddedPOST("[VASP_FORCE_OPTION]LSCOUPLING=ON",_aflowinpad_) << "// ON | OFF (default OFF)" << endl;
    } else {
        // aflowin << aurostd::PaddedPOST("[VASP_FORCE_OPTION]LSCOUPLING=OFF",_aflowinpad_) << "// ON | OFF (default OFF)" << endl;
    }

    // MAGMOM WRITING
    if(xvasp.aopts.flag("FLAG::AVASP_AUTO_MAGMOM")) {
        aflowin << aurostd::PaddedPOST("[VASP_FORCE_OPTION]AUTO_MAGMOM=ON",_aflowinpad_) << "// ON | OFF (default OFF)" << endl;
    } else {
        aflowin << aurostd::PaddedPOST("#[VASP_FORCE_OPTION]AUTO_MAGMOM=ON",_aflowinpad_) << "// ON | OFF (default OFF)" << endl;
    }

    // RELAX_TYPE
    if(xvasp.aopts.flag("AFLOWIN_FLAG::RELAX_TYPE")){  // CO 180214
        aflowin << aurostd::PaddedPOST("[VASP_FORCE_OPTION]RELAX_"+xvasp.aopts.getattachedscheme("AFLOWIN_FLAG::RELAX_TYPE"),_aflowinpad_);
        aflowin << "// ALL | IONS | CELL_SHAPE | CELL_VOLUME | IONS_CELL_VOLUME " << endl;
    }else{
        // RELAX_MODE
        if(xvasp.aopts.flag("AFLOWIN_FLAG::RELAX_MODE")) {
            string RELAX_MODE=xvasp.aopts.getattachedscheme("AFLOWIN_FLAG::RELAX_MODE");
            //   if(RELAX_MODE.at(0)=='E' || RELAX_MODE.at(0)=='e')
            aflowin << aurostd::PaddedPOST("[VASP_FORCE_OPTION]RELAX_MODE="+RELAX_MODE,_aflowinpad_) << "// (ENERGY | FORCES | ENERGY_FORCES | FORCES_ENERGY) (default: DEFAULT_VASP_FORCE_OPTION_RELAX_MODE_SCHEME in .aflow.rc) " << endl;
        } else {
            if(xvasp.aopts.flag("FLAG::AVASP_RELAX_FORCES")) {
                aflowin << aurostd::PaddedPOST("[VASP_FORCE_OPTION]RELAX_MODE=FORCES",_aflowinpad_) << "// (ENERGY | FORCES | ENERGY_FORCES | FORCES_ENERGY) (default: DEFAULT_VASP_FORCE_OPTION_RELAX_MODE_SCHEME in .aflow.rc) " << endl;
            } else {
                aflowin << aurostd::PaddedPOST("[VASP_FORCE_OPTION]RELAX_MODE=ENERGY",_aflowinpad_) << "// (ENERGY | FORCES | ENERGY_FORCES | FORCES_ENERGY) (default: DEFAULT_VASP_FORCE_OPTION_RELAX_MODE_SCHEME in .aflow.rc) " << endl;
            }
        }
    }

    // PRECISION WRITING
    string PRECISION=xvasp.AVASP_flag_PRECISION_scheme;
    if(xvasp.str.species.size()==1) PRECISION="A";
    if(xvasp.aopts.flag("AFLOWIN_FLAG::PRECISION"))
        PRECISION=xvasp.aopts.getattachedscheme("AFLOWIN_FLAG::PRECISION");
    if(xvasp.aopts.flag("FLAG::PRECISION_SET") || xvasp.aopts.flag("AFLOWIN_FLAG::PRECISION")) {
        if(AFLOWIN_OLD) { // OLD
            if(PRECISION.at(0)=='L' || PRECISION.at(0)=='l')
                aflowin << aurostd::PaddedPOST("[VASP_FORCE_OPTION]PREC=LOW",_aflowinpad_) << "// (LOW | MEDIUM | NORMAL | HIGH | ACCURATE), PRESERVED (default: DEFAULT_VASP_FORCE_OPTION_PREC_SCHEME in .aflow.rc" << endl;
            if(PRECISION.at(0)=='M' || PRECISION.at(0)=='m')
                aflowin << aurostd::PaddedPOST("[VASP_FORCE_OPTION]PREC=MEDIUM",_aflowinpad_) << "// (LOW | MEDIUM | NORMAL | HIGH | ACCURATE), PRESERVED (default: DEFAULT_VASP_FORCE_OPTION_PREC_SCHEME in .aflow.rc" << endl;
            if(PRECISION.at(0)=='N' || PRECISION.at(0)=='n')
                aflowin << aurostd::PaddedPOST("[VASP_FORCE_OPTION]PREC=NORMAL",_aflowinpad_) << "// (LOW | MEDIUM | NORMAL | HIGH | ACCURATE), PRESERVED (default: DEFAULT_VASP_FORCE_OPTION_PREC_SCHEME in .aflow.rc" << endl;
            if(PRECISION.at(0)=='H' || PRECISION.at(0)=='h')
                aflowin << aurostd::PaddedPOST("[VASP_FORCE_OPTION]PREC=HIGH",_aflowinpad_) << "// (LOW | MEDIUM | NORMAL | HIGH | ACCURATE), PRESERVED (default: DEFAULT_VASP_FORCE_OPTION_PREC_SCHEME in .aflow.rc" << endl;
            if(PRECISION.at(0)=='A' || PRECISION.at(0)=='a')
                aflowin << aurostd::PaddedPOST("[VASP_FORCE_OPTION]PREC=ACCURATE",_aflowinpad_) << "// (LOW | MEDIUM | NORMAL | HIGH | ACCURATE), PRESERVED (default: DEFAULT_VASP_FORCE_OPTION_PREC_SCHEME in .aflow.rc" << endl;
            aflowin << "#[VASP_FORCE_OPTION]PREC_preserved " << endl;
        }
        if(AFLOWIN_NEW) { // NEW
            vector<string> vstr;string tstr;
            vstr.clear();tstr="";
            if(PRECISION.at(0)=='L' || PRECISION.at(0)=='l')
                vstr.push_back("LOW");
            if(PRECISION.at(0)=='M' || PRECISION.at(0)=='m')
                vstr.push_back("MEDIUM");
            if(PRECISION.at(0)=='N' || PRECISION.at(0)=='n')
                vstr.push_back("NORMAL");
            if(PRECISION.at(0)=='H' || PRECISION.at(0)=='h')
                vstr.push_back("HIGH");
            if(PRECISION.at(0)=='A' || PRECISION.at(0)=='a')
                vstr.push_back("ACCURATE");
            if(xvasp.aopts.flag("FLAG::PRECISION_PRESERVED")) vstr.push_back("PRESERVED");
            if(vstr.size()>0) {
                for(uint i=0;i<vstr.size();i++)
                    tstr+=vstr.at(i)+(i!=vstr.size()-1 ? "," : " ");
                aflowin << aurostd::PaddedPOST("[VASP_FORCE_OPTION]PREC="+tstr,_aflowinpad_) << "// (LOW | MEDIUM | NORMAL | HIGH | ACCURATE), PRESERVED (default: DEFAULT_VASP_FORCE_OPTION_PREC_SCHEME in .aflow.rc)" << endl;
            }
            if(vstr.size()==0) {
                aflowin << aurostd::PaddedPOST("#[VASP_FORCE_OPTION]PREC=something",_aflowinpad_) << "// (LOW | MEDIUM | NORMAL | HIGH | ACCURATE), PRESERVED (default: DEFAULT_VASP_FORCE_OPTION_PREC_SCHEME in .aflow.rc)" << endl;
            }
        }
    }

    // ALGO WRITING
    string ALGORITHM=xvasp.AVASP_flag_ALGO_scheme;
    if(xvasp.str.species.size()==1) ALGORITHM="N";
    if(xvasp.aopts.flag("AFLOWIN_FLAG::ALGORITHM"))
        ALGORITHM=xvasp.aopts.getattachedscheme("AFLOWIN_FLAG::ALGORITHM");
    if(xvasp.aopts.flag("FLAG::ALGO_SET") || xvasp.aopts.flag("AFLOWIN_FLAG::ALGORITHM")) {
        if(AFLOWIN_OLD) { // OLD
            if(ALGORITHM.at(0)=='N' || ALGORITHM.at(0)=='n')
                aflowin << aurostd::PaddedPOST("[VASP_FORCE_OPTION]ALGO=NORMAL",_aflowinpad_) << "// (NORMAL | VERYFAST | FAST | ALL | DAMPED), PRESERVED (default: DEFAULT_VASP_FORCE_OPTION_ALGO_SCHEME in .aflow.rc)" << endl;
            if(ALGORITHM.at(0)=='V' || ALGORITHM.at(0)=='v')
                aflowin << aurostd::PaddedPOST("[VASP_FORCE_OPTION]ALGO=VERYFAST",_aflowinpad_) << "// (NORMAL | VERYFAST | FAST | ALL | DAMPED), PRESERVED (default: DEFAULT_VASP_FORCE_OPTION_ALGO_SCHEME in .aflow.rc)" << endl;
            if(ALGORITHM.at(0)=='F' || ALGORITHM.at(0)=='f')
                aflowin << aurostd::PaddedPOST("[VASP_FORCE_OPTION]ALGO=FAST",_aflowinpad_) << "// (NORMAL | VERYFAST | FAST | ALL | DAMPED), PRESERVED (default: DEFAULT_VASP_FORCE_OPTION_ALGO_SCHEME in .aflow.rc)" << endl;
            if(ALGORITHM.at(0)=='A' || ALGORITHM.at(0)=='a')
                aflowin << aurostd::PaddedPOST("[VASP_FORCE_OPTION]ALGO=ALL",_aflowinpad_) << "// (NORMAL | VERYFAST | FAST | ALL | DAMPED), PRESERVED (default: DEFAULT_VASP_FORCE_OPTION_ALGO_SCHEME in .aflow.rc)" << endl;
            if(ALGORITHM.at(0)=='D' || ALGORITHM.at(0)=='d')
                aflowin << aurostd::PaddedPOST("[VASP_FORCE_OPTION]ALGO=DAMPED",_aflowinpad_) << "// (NORMAL | VERYFAST | FAST | ALL | DAMPED), PRESERVED (default: DEFAULT_VASP_FORCE_OPTION_ALGO_SCHEME in .aflow.rc)" << endl;
            aflowin << "#[VASP_FORCE_OPTION]ALGO_PRESERVED " << endl;
        }
        if(AFLOWIN_NEW) { // NEW
            vector<string> vstr;string tstr;
            vstr.clear();tstr="";
            if(ALGORITHM.at(0)=='N' || ALGORITHM.at(0)=='n')
                vstr.push_back("NORMAL");
            if(ALGORITHM.at(0)=='V' || ALGORITHM.at(0)=='v')
                vstr.push_back("VERYFAST");
            if(ALGORITHM.at(0)=='F' || ALGORITHM.at(0)=='f')
                vstr.push_back("FAST");
            if(ALGORITHM.at(0)=='A' || ALGORITHM.at(0)=='a')
                vstr.push_back("ALL");
            if(ALGORITHM.at(0)=='D' || ALGORITHM.at(0)=='d')
                vstr.push_back("DAMPED");
            if(xvasp.aopts.flag("FLAG::ALGO_PRESERVED")) vstr.push_back("PRESERVED");
            if(vstr.size()>0) {
                for(uint i=0;i<vstr.size();i++) tstr+=vstr.at(i)+(i!=vstr.size()-1 ? "," : " ");
                aflowin << aurostd::PaddedPOST("[VASP_FORCE_OPTION]ALGO="+tstr,_aflowinpad_) << "// (NORMAL | VERYFAST | FAST | ALL | DAMPED), PRESERVED (default: DEFAULT_VASP_FORCE_OPTION_ALGO_SCHEME in .aflow.rc)" << endl;
            }
            if(vstr.size()==0) {
                aflowin << aurostd::PaddedPOST("#[VASP_FORCE_OPTION]ALGO=something",_aflowinpad_) << "// (NORMAL | VERYFAST | FAST | ALL | DAMPED), PRESERVED (default: DEFAULT_VASP_FORCE_OPTION_ALGO_SCHEME in .aflow.rc)" << endl;
            }
        }
    }


    // METAGGA WRITING
    string METAGGA=DEFAULT_VASP_FORCE_OPTION_METAGGA_SCHEME; // default

    if(xvasp.aopts.flag("AFLOWIN_FLAG::METAGGA"))
        METAGGA=xvasp.aopts.getattachedscheme("AFLOWIN_FLAG::METAGGA");

    if(xvasp.aopts.flag("AFLOWIN_FLAG::METAGGA")) {
        vector<string> vstr;string tstr;
        vstr.clear();tstr="";
        if(aurostd::toupper(METAGGA)=="TPSS") { vstr.push_back("TPSS"); }
        if(aurostd::toupper(METAGGA)=="RTPSS") { vstr.push_back("RTPSS"); }
        if(aurostd::toupper(METAGGA)=="M06L") { vstr.push_back("M06L"); }
        if(aurostd::toupper(METAGGA)=="MBJL") { vstr.push_back("MBJL"); }
        if(aurostd::toupper(METAGGA)=="SCAN") { vstr.push_back("SCAN"); }
        if(aurostd::toupper(METAGGA)=="MS0") { vstr.push_back("MS0"); }
        if(aurostd::toupper(METAGGA)=="MS1") { vstr.push_back("MS1"); }
        if(aurostd::toupper(METAGGA)=="MS2") { vstr.push_back("MS2"); }
        if(aurostd::toupper(METAGGA)=="NONE" || aurostd::toupper(METAGGA)=="" || aurostd::toupper(METAGGA)=="OFF" || aurostd::toupper(METAGGA)=="0") { }
        if(vstr.size()>0) {
            for(uint i=0;i<vstr.size();i++)
                tstr+=vstr.at(i)+(i!=vstr.size()-1 ? "," : " ");
            aflowin << aurostd::PaddedPOST("[VASP_FORCE_OPTION]METAGGA="+tstr,_aflowinpad_) << "// (TPSS | RTPSS | M06L | MBJL | SCAN | MS0 | MS1 | MS2 | NONE) (default: DEFAULT_VASP_FORCE_OPTION_METAGGA_SCHEME in .aflow.rc)" << endl;
        }
        if(vstr.size()==0) {
            aflowin << aurostd::PaddedPOST("#[VASP_FORCE_OPTION]METAGGA=NONE",_aflowinpad_) << "// (TPSS | RTPSS | M06L | MBJL | SCAN | MS0 | MS1 | MS2 | NONE) (default: DEFAULT_VASP_FORCE_OPTION_METAGGA_SCHEME in .aflow.rc)" << endl;
        }
    } else {
        aflowin << aurostd::PaddedPOST("#[VASP_FORCE_OPTION]METAGGA=NONE",_aflowinpad_) << "// (TPSS | RTPSS | M06L | MBJL | SCAN | MS0 | MS1 | MS2 | NONE) (default: DEFAULT_VASP_FORCE_OPTION_METAGGA_SCHEME in .aflow.rc)" << endl;
    }

    // IVDW WRITING
    string IVDW=DEFAULT_VASP_FORCE_OPTION_IVDW_SCHEME; // default

    if(xvasp.aopts.flag("AFLOWIN_FLAG::IVDW"))
        IVDW=xvasp.aopts.getattachedscheme("AFLOWIN_FLAG::IVDW");

    if(xvasp.aopts.flag("FLAG::IVDW_SET") || xvasp.aopts.flag("AFLOWIN_FLAG::IVDW")) {
        vector<string> vstr;string tstr;
        vstr.clear();tstr="";
        if(aurostd::toupper(IVDW)=="OFF" || aurostd::toupper(IVDW)=="NONE" || aurostd::toupper(IVDW)=="0" || aurostd::toupper(IVDW)=="") { ; // nothing
        } else {
            vstr.push_back(aurostd::toupper(IVDW));
        }
        if(vstr.size()>0) {
            for(uint i=0;i<vstr.size();i++)
                tstr+=vstr.at(i)+(i!=vstr.size()-1 ? "," : " ");
            aflowin << aurostd::PaddedPOST("[VASP_FORCE_OPTION]IVDW="+tstr,_aflowinpad_) << "// (number_for_VASP_see_manual_for_IVDW | 0) (default: DEFAULT_VASP_FORCE_OPTION_IVDW_SCHEME in .aflow.rc)" << endl;
        }
        if(vstr.size()==0) {
            aflowin << aurostd::PaddedPOST("#[VASP_FORCE_OPTION]IVDW=0",_aflowinpad_) << "// (number_for_VASP_see_manual_for_IVDW | 0) (default: DEFAULT_VASP_FORCE_OPTION_IVDW_SCHEME in .aflow.rc)" << endl;
        }
    } else {
        aflowin << aurostd::PaddedPOST("#[VASP_FORCE_OPTION]IVDW=0",_aflowinpad_) << "// (number_for_VASP_see_manual_for_IVDW | 0) (default: DEFAULT_VASP_FORCE_OPTION_IVDW_SCHEME in .aflow.rc)" << endl;
    }

    // ABMIX WRITING
    if(xvasp.str.species.size()==1) xvasp.AVASP_flag_ABMIX_scheme="N";
    if(xvasp.aopts.flag("FLAG::ABMIX_SET")) { //  unspecified | [AUTO | US | PAW | #AMIX,#BMIX[,#AMIX_MAG,#BMIX_MAG]]
        vector<string> vstr;string tstr;
        vstr.clear();tstr="";
        if(xvasp.AVASP_flag_ABMIX_scheme.at(0)=='A' || xvasp.AVASP_flag_ABMIX_scheme.at(0)=='a')
            vstr.push_back("AUTO");
        if(xvasp.AVASP_flag_ABMIX_scheme.at(0)=='U' || xvasp.AVASP_flag_ABMIX_scheme.at(0)=='u')
            vstr.push_back("US");
        if(xvasp.AVASP_flag_ABMIX_scheme.at(0)=='P' || xvasp.AVASP_flag_ABMIX_scheme.at(0)=='p')
            vstr.push_back("PAW");
        if(vstr.size()>0) {
            for(uint i=0;i<vstr.size();i++)
                tstr+=vstr.at(i)+(i!=vstr.size()-1 ? "," : " ");
            aflowin << aurostd::PaddedPOST("[VASP_FORCE_OPTION]ABMIX="+tstr,_aflowinpad_) << "// unspecified | [AUTO | US | PAW | (exp) #AMIX,#BMIX[,#AMIX_MAG,#BMIX_MAG]] (default=unspecified " << DEFAULT_VASP_FORCE_OPTION_ABMIX_SCHEME<< ")" << endl;
        }
        if(vstr.size()==0) {
            aflowin << aurostd::PaddedPOST("#[VASP_FORCE_OPTION]#ABMIX=...",_aflowinpad_) << "// unspecified | [AUTO | US | PAW | (exp) #AMIX,#BMIX[,#AMIX_MAG,#BMIX_MAG]] (default=unspecified " << DEFAULT_VASP_FORCE_OPTION_ABMIX_SCHEME<< ")" << endl;
        }
    }

    // RELAX WRITING
    if(xvasp.aopts.flag("FLAG::VOLUME_PRESERVED")==TRUE || xvasp.AVASP_flag_RUN_STATIC==TRUE || xvasp.AVASP_flag_RUN_STATIC_BANDS==TRUE) {
    } else {
        if(!xvasp.aopts.flag("AFLOWIN_FLAG::RELAX_TYPE")) aflowin << "[VASP_FORCE_OPTION]RELAX_ALL " << endl; // CO 180214
    }

    // NOTUNE WRITING
    aflowin << "#[VASP_FORCE_OPTION]NOTUNE " << endl;

    // TYPE WRITING
    string TYPE=xvasp.AVASP_flag_TYPE.xscheme;
    if(xvasp.aopts.flag("AFLOWIN_FLAG::TYPE")) 
        TYPE=xvasp.aopts.getattachedscheme("AFLOWIN_FLAG::TYPE");
    if(TYPE.at(0)!='M' && TYPE.at(0)!='I' && TYPE.at(0)!='S' && TYPE.at(0)!='m' && TYPE.at(0)!='i' && TYPE.at(0)!='s')
        aflowin << aurostd::PaddedPOST("[VASP_FORCE_OPTION]TYPE=DEFAULT",_aflowinpad_) << "// (METAL | INSULATOR | SEMICONDUCTOR | DEFAULT) (default " << DEFAULT_VASP_FORCE_OPTION_TYPE_SCHEME<< ") " << endl;
    if(TYPE.at(0)=='M' || TYPE.at(0)=='m')
        aflowin << aurostd::PaddedPOST("[VASP_FORCE_OPTION]TYPE=METAL",_aflowinpad_) << "// (METAL | INSULATOR | SEMICONDUCTOR | DEFAULT) (default " << DEFAULT_VASP_FORCE_OPTION_TYPE_SCHEME<< ") " << endl;
    if(TYPE.at(0)=='I' || TYPE.at(0)=='i')
        aflowin << aurostd::PaddedPOST("[VASP_FORCE_OPTION]TYPE=INSULATOR",_aflowinpad_) << "// (METAL | INSULATOR | SEMICONDUCTOR | DEFAULT) (default " << DEFAULT_VASP_FORCE_OPTION_TYPE_SCHEME<< ") " << endl;
    if(TYPE.at(0)=='S' || TYPE.at(0)=='s')
        aflowin << aurostd::PaddedPOST("[VASP_FORCE_OPTION]TYPE=SEMICONDUCTOR",_aflowinpad_) << "// (METAL | INSULATOR | SEMICONDUCTOR | DEFAULT) (default " << DEFAULT_VASP_FORCE_OPTION_TYPE_SCHEME<< ") " << endl;

    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " " << "[12.5]" << endl;

    if(AFLOWIN_OLD) { // OLD STYLE
        if(xvasp.aopts.flag("FLAG::AVASP_LDAU1")) aflowin << "[VASP_FORCE_OPTION]LDAU1=ON" << endl;
        if(xvasp.aopts.flag("FLAG::AVASP_LDAU2")) aflowin << "[VASP_FORCE_OPTION]LDAU2=ON" << endl;
        if(xvasp.aopts.flag("FLAG::AVASP_LDAU1") || xvasp.aopts.flag("FLAG::AVASP_LDAU2")) {
            stringstream aus;
            aus << "[VASP_FORCE_OPTION]LDAU_SPECIES=";
            for(uint i=0;i<xvasp.str.species.size();i++) aus << KBIN::VASP_PseudoPotential_CleanName(xvasp.str.species.at(i)) << " ";
            aflowin << aurostd::PaddedPOST(aus.str(),_aflowinpad_) << "// LDAU SPECIES separated by spaces" << endl;
        }
    }
    if(AFLOWIN_NEW) { // NEW STYLE
        if(xvasp.aopts.flag("FLAG::AVASP_LDAU_ADIABATIC")==FALSE && xvasp.aopts.flag("FLAG::AVASP_LDAU_CUTOFF")==FALSE) {
            if(xvasp.aopts.flag("FLAG::AVASP_LDAU1") && !xvasp.aopts.flag("FLAG::AVASP_LDAU2")) aflowin << aurostd::PaddedPOST("[VASP_FORCE_OPTION]LDAU1=ON",_aflowinpad_) << "// ON | OFF | ADIABATIC | CUTOFF  " << endl;
            if(xvasp.aopts.flag("FLAG::AVASP_LDAU2") && !xvasp.aopts.flag("FLAG::AVASP_LDAU1")) aflowin << aurostd::PaddedPOST("[VASP_FORCE_OPTION]LDAU2=ON",_aflowinpad_) << "// ON | OFF | ADIABATIC | CUTOFF  " << endl;
            if(xvasp.aopts.flag("FLAG::AVASP_LDAU2") && xvasp.aopts.flag("FLAG::AVASP_LDAU1"))  aflowin << aurostd::PaddedPOST("[VASP_FORCE_OPTION]LDAU2=ON",_aflowinpad_) << "// ON | OFF | ADIABATIC | CUTOFF  " << endl;
        }
        if(xvasp.aopts.flag("FLAG::AVASP_LDAU_ADIABATIC")==TRUE && xvasp.aopts.flag("FLAG::AVASP_LDAU_CUTOFF")==FALSE) {
            if(xvasp.aopts.flag("FLAG::AVASP_LDAU1") && !xvasp.aopts.flag("FLAG::AVASP_LDAU2")) aflowin << aurostd::PaddedPOST("[VASP_FORCE_OPTION]LDAU1=ADIABATIC",_aflowinpad_) << "// ON | OFF | ADIABATIC | CUTOFF  " << endl;
            if(xvasp.aopts.flag("FLAG::AVASP_LDAU2") && !xvasp.aopts.flag("FLAG::AVASP_LDAU1")) aflowin << aurostd::PaddedPOST("[VASP_FORCE_OPTION]LDAU2=ADIABATIC",_aflowinpad_) << "// ON | OFF | ADIABATIC | CUTOFF  " << endl;
            if(xvasp.aopts.flag("FLAG::AVASP_LDAU2") && xvasp.aopts.flag("FLAG::AVASP_LDAU1"))  aflowin << aurostd::PaddedPOST("[VASP_FORCE_OPTION]LDAU2=ADIABATIC",_aflowinpad_) << "// ON | OFF | ADIABATIC | CUTOFF  " << endl;
        }
        if(xvasp.aopts.flag("FLAG::AVASP_LDAU_ADIABATIC")==FALSE && xvasp.aopts.flag("FLAG::AVASP_LDAU_CUTOFF")==TRUE) {
            if(xvasp.aopts.flag("FLAG::AVASP_LDAU1") && !xvasp.aopts.flag("FLAG::AVASP_LDAU2")) aflowin << aurostd::PaddedPOST("[VASP_FORCE_OPTION]LDAU1=CUTOFF",_aflowinpad_) << "// ON | OFF | ADIABATIC | CUTOFF  " << endl;
            if(xvasp.aopts.flag("FLAG::AVASP_LDAU2") && !xvasp.aopts.flag("FLAG::AVASP_LDAU1")) aflowin << aurostd::PaddedPOST("[VASP_FORCE_OPTION]LDAU2=CUTOFF",_aflowinpad_) << "// ON | OFF | ADIABATIC | CUTOFF  " << endl;
            if(xvasp.aopts.flag("FLAG::AVASP_LDAU2") && xvasp.aopts.flag("FLAG::AVASP_LDAU1"))  aflowin << aurostd::PaddedPOST("[VASP_FORCE_OPTION]LDAU2=CUTOFF",_aflowinpad_) << "// ON | OFF | ADIABATIC | CUTOFF  " << endl;
        }
        if(xvasp.aopts.flag("FLAG::AVASP_LDAU1") || xvasp.aopts.flag("FLAG::AVASP_LDAU2")) aflowin << aurostd::PaddedPOST("[VASP_FORCE_OPTION]LDAU_PARAMETERS="+xvasp.AVASP_LDAU_PARAMETERS_STRING,_aflowinpad_) << "// species;Ls;Us;Js " << endl;
    }

    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " " << "[12.6]" << endl;

    // CONVERT_UNIT_CELL
    // gus check
    {
        vector<string> tokens;aurostd::string2tokens(system,tokens,".");
        if(tokens.size()>1) {
            if(tokens.at(1).at(0)=='f') xvasp.aopts.flag("AVASP_flag_CONVERT_UNIT_CELL_STANDARD_PRIMITIVE",FALSE);
            if(tokens.at(1).at(0)=='b') xvasp.aopts.flag("AVASP_flag_CONVERT_UNIT_CELL_STANDARD_PRIMITIVE",FALSE);
            if(tokens.at(1).at(0)=='s') xvasp.aopts.flag("AVASP_flag_CONVERT_UNIT_CELL_STANDARD_PRIMITIVE",FALSE);
        }
    }

    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " " << "[12.7]" << endl;
    if(AFLOWIN_NEW) {
        vector<string> vstr;string tstr;
        vstr.clear();tstr="";

        if(xvasp.aopts.flag("AFLOWIN_FLAG::CONVERT_UNIT_CELL")) {
            vector<string>  vCONVERT_UNIT_CELL;
            aurostd::string2tokens(aurostd::toupper(xvasp.aopts.getattachedscheme("AFLOWIN_FLAG::CONVERT_UNIT_CELL")),vCONVERT_UNIT_CELL,",");
            for(uint i=0;i<vCONVERT_UNIT_CELL.size();i++) {
                if(vCONVERT_UNIT_CELL.at(i)=="SPRIM") vstr.push_back("STANDARD_PRIMITIVE");
                if(vCONVERT_UNIT_CELL.at(i)=="STD_PRIM") vstr.push_back("STANDARD_PRIMITIVE");
                if(vCONVERT_UNIT_CELL.at(i)=="STANDARD_PRIMITIVE") vstr.push_back("STANDARD_PRIMITIVE");
                if(vCONVERT_UNIT_CELL.at(i)=="SCONV") vstr.push_back("STANDARD_CONVENTIONAL");
                if(vCONVERT_UNIT_CELL.at(i)=="STD_CONV") vstr.push_back("STANDARD_CONVENTIONAL");
                if(vCONVERT_UNIT_CELL.at(i)=="STANDARD_CONVENTIONAL") vstr.push_back("STANDARD_CONVENTIONAL");
                if(vCONVERT_UNIT_CELL.at(i)=="NIGGLI") vstr.push_back("NIGGLI");
                if(vCONVERT_UNIT_CELL.at(i)=="MINK") vstr.push_back("MINKOWSKI");
                if(vCONVERT_UNIT_CELL.at(i)=="MINKOWSKI") vstr.push_back("MINKOWSKI");
                if(vCONVERT_UNIT_CELL.at(i)=="INCELL") vstr.push_back("INCELL");
                if(vCONVERT_UNIT_CELL.at(i)=="COMPACT") vstr.push_back("COMPACT");
                if(vCONVERT_UNIT_CELL.at(i)=="INCOMPACT") vstr.push_back("COMPACT");
                if(vCONVERT_UNIT_CELL.at(i)=="INWIGNERSEITZ") vstr.push_back("WIGNERSEITZ");
                if(vCONVERT_UNIT_CELL.at(i)=="WS") vstr.push_back("WIGNERSEITZ");
                if(vCONVERT_UNIT_CELL.at(i)=="WIGNERSEITZ") vstr.push_back("WIGNERSEITZ");
                if(vCONVERT_UNIT_CELL.at(i)=="C") vstr.push_back("CARTESIAN");
                if(vCONVERT_UNIT_CELL.at(i)=="CART") vstr.push_back("CARTESIAN");
                if(vCONVERT_UNIT_CELL.at(i)=="CARTESIAN") vstr.push_back("CARTESIAN");
                if(vCONVERT_UNIT_CELL.at(i)=="F") vstr.push_back("FRACTIONAL");
                if(vCONVERT_UNIT_CELL.at(i)=="FRAC") vstr.push_back("FRACTIONAL");
                if(vCONVERT_UNIT_CELL.at(i)=="FRACTIONAL") vstr.push_back("FRACTIONAL");
                if(vCONVERT_UNIT_CELL.at(i)=="D") vstr.push_back("DIRECT");
                if(vCONVERT_UNIT_CELL.at(i)=="DIR") vstr.push_back("DIRECT");
                if(vCONVERT_UNIT_CELL.at(i)=="DIRECT") vstr.push_back("DIRECT");
                if(vCONVERT_UNIT_CELL.at(i)=="PRE") vstr.push_back("PRESERVE");
                if(vCONVERT_UNIT_CELL.at(i)=="PRES") vstr.push_back("PRESERVE");
                if(vCONVERT_UNIT_CELL.at(i)=="PRESERVE") vstr.push_back("PRESERVE");
            }
        } else {
            if(xvasp.aopts.flag("AVASP_flag_CONVERT_UNIT_CELL_STANDARD_PRIMITIVE")) vstr.push_back("SPRIM");
            //if(xvasp.aopts.flag("AVASP_flag_CONVERT_UNIT_CELL_STANDARD_PRIMITIVE")) vstr.push_back("STANDARD_PRIMITIVE");
            if(xvasp.aopts.flag("AVASP_flag_CONVERT_UNIT_CELL_STANDARD_CONVENTIONAL")) vstr.push_back("SCONV");
            //if(xvasp.aopts.flag("AVASP_flag_CONVERT_UNIT_CELL_STANDARD_CONVENTIONAL")) vstr.push_back("STANDARD_CONVENTIONAL");
            if(xvasp.aopts.flag("AVASP_flag_CONVERT_UNIT_CELL_NIGGLI")) vstr.push_back("NIGGLI");
            if(xvasp.aopts.flag("AVASP_flag_CONVERT_UNIT_CELL_MINKOWSKI")) vstr.push_back("MINK");
            //if(xvasp.aopts.flag("AVASP_flag_CONVERT_UNIT_CELL_MINKOWSKI")) vstr.push_back("MINKOWSKI");
            if(xvasp.aopts.flag("AVASP_flag_CONVERT_UNIT_CELL_INCELL")) vstr.push_back("INCELL");
            if(xvasp.aopts.flag("AVASP_flag_CONVERT_UNIT_CELL_INCOMPACT")) vstr.push_back("COMPACT");
            if(xvasp.aopts.flag("AVASP_flag_CONVERT_UNIT_CELL_COMPACT")) vstr.push_back("COMPACT");
            if(xvasp.aopts.flag("AVASP_flag_CONVERT_UNIT_CELL_CARTESIAN")) vstr.push_back("CART");
            //if(xvasp.aopts.flag("AVASP_flag_CONVERT_UNIT_CELL_CARTESIAN")) vstr.push_back("CARTESIAN");
            if(xvasp.aopts.flag("AVASP_flag_CONVERT_UNIT_CELL_FRACTIONAL")) vstr.push_back("FRAC");
            //if(xvasp.aopts.flag("AVASP_flag_CONVERT_UNIT_CELL_FRACTIONAL")) vstr.push_back("FRACTIONAL");
            if(xvasp.aopts.flag("AVASP_flag_CONVERT_UNIT_CELL_WIGNERSEITZ")) vstr.push_back("WS");
            //if(xvasp.aopts.flag("AVASP_flag_CONVERT_UNIT_CELL_WIGNERSEITZ")) vstr.push_back("WIGNERSEITZ");
            if(xvasp.aopts.flag("AVASP_flag_CONVERT_UNIT_CELL_PRESERVE")) vstr.push_back("PRES");
            //if(xvasp.aopts.flag("AVASP_flag_CONVERT_UNIT_CELL_PRESERVE")) vstr.push_back("PRESERVE");
        }
        if(vstr.size()>0) {
            for(uint i=0;i<vstr.size();i++) tstr+=vstr.at(i)+(i!=vstr.size()-1 ? "," : " ");
            aflowin << aurostd::PaddedPOST("[VASP_FORCE_OPTION]CONVERT_UNIT_CELL="+tstr,_aflowinpad_) << "// (SPRIM, SCONV, NIGGLI, MINK, INCELL, COMPACT, WS, CART, FRAC, PRES) " << endl;
        }
        if(vstr.size()==0) {
            aflowin << aurostd::PaddedPOST("#[VASP_FORCE_OPTION]CONVERT_UNIT_CELL=something",_aflowinpad_) << "// (SPRIM, SCONV, NIGGLI, MINK, INCELL, COMPACT, WS, CART, FRAC, PRES)" << endl;
        }
    }

    // VOLUME_PLUS_EQUAL
    if(xvasp.aopts.flag("AFLOWIN_FLAG::VOLUME_PLUS_EQUAL") && !xvasp.aopts.flag("AFLOWIN_FLAG::NO_VOLUME_ADJUSTMENT")) {
        aflowin << "[VASP_FORCE_OPTION]VOLUME+=" << xvasp.aopts.getattachedscheme("AFLOWIN_FLAG::VOLUME_PLUS_EQUAL") << " " << endl;
        xvasp.aopts.pop_attached("AFLOWIN_FLAG::VOLUME_MULTIPLY_EQUAL");
    } else {
        aflowin << "#[VASP_FORCE_OPTION]VOLUME+=10.0 " << endl;
    }

    // VOLUME_MULTIPLY_EQUAL
    if(xvasp.aopts.flag("AFLOWIN_FLAG::VOLUME_MULTIPLY_EQUAL") && !xvasp.aopts.flag("AFLOWIN_FLAG::NO_VOLUME_ADJUSTMENT")) {
        aflowin << "[VASP_FORCE_OPTION]VOLUME*=" << xvasp.aopts.getattachedscheme("AFLOWIN_FLAG::VOLUME_MULTIPLY_EQUAL") << " " << endl;
        xvasp.aopts.pop_attached("AFLOWIN_FLAG::VOLUME_PLUS_EQUAL");
    } else {
        if(xvasp.aopts.flag("FLAG::VOLUME_PRESERVED")==TRUE || xvasp.AVASP_flag_RUN_STATIC==TRUE || xvasp.AVASP_flag_RUN_STATIC_BANDS==TRUE || xvasp.aopts.flag("AFLOWIN_FLAG::NO_VOLUME_ADJUSTMENT")) {
            aflowin << "#[VASP_FORCE_OPTION]VOLUME*=1.05 " << endl;
        } else {
            if(!xvasp.aopts.flag("AFLOWIN_FLAG::VOLUME_PLUS_EQUAL")) {aflowin << "[VASP_FORCE_OPTION]VOLUME*=1.05 " << endl;} else {aflowin << "#[VASP_FORCE_OPTION]VOLUME*=1.05 " << endl;}
        }
    }

    if(0) { // write the AFIX
        aflowin << AFLOWIN_SEPARATION_LINE << endl; // [AFLOW] **************************************************
        aflowin << "[AFLOW] # Do not uncomment the AFIX unless you know what  you are doing. " << endl;
        aflowin <<  aurostd::PaddedPOST("[VASP_FORCE_OPTION]IGNORE_AFIX=NONE",_aflowinpad_) << "// ROTMAT, SGRCON, IBZKPT, SYMPREC, INVGRP, EDDRMM, LREAL, BRMIX, DAV, EDDDAV, EFIELD_PEAD, ZPOTRF, EXCCOR, NATOMS, NBANDS, MEMORY, PSMAXN, NPAR" << endl;
    }
    aflowin << AFLOWIN_SEPARATION_LINE << endl; // [AFLOW] **************************************************
    aflowin << AFLOWIN_SEPARATION_LINE << endl; // [AFLOW] **************************************************
    aflowin << "[VASP_INCAR_MODE_EXPLICIT]START " << endl;
    aflowin << "SYSTEM=" << system << endl;
    //aflowin << " " << endl;
    if(xvasp.aopts.flag("FLAG::EXTRA_INCAR")==TRUE) {
        aflowin << xvasp.AVASP_EXTRA_INCAR.str(); // << endl;
    }
    aflowin << aurostd::PaddedPOST("#PSTRESS=000     # Pressure in kBar (1kB=0.1GPa) ",_aflowinpad_) << "# for hand modification" << endl;
    aflowin << aurostd::PaddedPOST("#EDIFFG="+aurostd::utype2string(DEFAULT_VASP_PREC_EDIFFG,_AVASP_DOUBLE2STRING_PRECISION_)+"     # For relaxed forces ",_aflowinpad_) << "# for hand modification" << endl;
    aflowin << aurostd::PaddedPOST("#POTIM="+aurostd::utype2string(DEFAULT_VASP_PREC_EDIFFG,_AVASP_DOUBLE2STRING_PRECISION_)+"      # default ",_aflowinpad_) << "# for hand modification" << endl;
    aflowin << aurostd::PaddedPOST("#NBANDS=XX  ",_aflowinpad_) << "# for hand modification" << endl;
    aflowin << aurostd::PaddedPOST("#IALGO=48   ",_aflowinpad_) << "# for hand modification" << endl;
    aflowin << "[VASP_INCAR_MODE_EXPLICIT]STOP " << endl;
    aflowin << AFLOWIN_SEPARATION_LINE << endl; // [AFLOW] **************************************************
    aflowin << "[VASP_KPOINTS_MODE_IMPLICIT] " << endl;
    aflowin << "[VASP_KPOINTS_FILE]KSCHEME=" << xvasp.AVASP_KSCHEME << " " << endl;

    // MODIFIERS
    if(xvasp.aopts.flag("AFLOWIN_FLAG::KPPRA")) xvasp.AVASP_value_KPPRA=aurostd::string2utype<int>(xvasp.aopts.getattachedscheme("AFLOWIN_FLAG::KPPRA"));    // KPPRA
    if(xvasp.aopts.flag("AFLOWIN_FLAG::KPPRA_STATIC")) xvasp.AVASP_value_KPPRA_STATIC=aurostd::string2utype<int>(xvasp.aopts.getattachedscheme("AFLOWIN_FLAG::KPPRA_STATIC")); // KPPRA_STATIC
    if(xvasp.aopts.flag("AFLOWIN_FLAG::BANDS_GRID")) xvasp.AVASP_value_BANDS_GRID=aurostd::string2utype<uint>(xvasp.aopts.getattachedscheme("AFLOWIN_FLAG::BANDS_GRID")); // BANDS_GRID

    if(xvasp.str.species.size()==1) if(!xvasp.aopts.flag("AFLOWIN_FLAG::KPPRA")) xvasp.AVASP_value_KPPRA=DEFAULT_UNARY_KPPRA; // DEFAULT KPPRA
    // if(xvasp.aopts.flag("FLAG::AVASP_KPPRA") && PARAMS->vkppra.size()>0) {
    //   aflowin << "[VASP_KPOINTS_FILE]KPPRA=" << PARAMS->vkppra.front() << endl;
    //   PARAMS->vkppra.pop_front();
    // } else {
    aflowin << "[VASP_KPOINTS_FILE]KPPRA=" << xvasp.AVASP_value_KPPRA << endl;
    //  }
    if(xvasp.AVASP_prototype_mode==LIBRARY_MODE_ICSD) {
        if(xvasp.AVASP_flag_RUN_RELAX_STATIC || xvasp.AVASP_flag_RUN_RELAX_STATIC_BANDS) {
            aflowin << "[VASP_KPOINTS_FILE]STATIC_KSCHEME=" << xvasp.AVASP_STATIC_KSCHEME << " " << endl;
            aflowin << "[VASP_KPOINTS_FILE]STATIC_KPPRA=" << xvasp.AVASP_value_KPPRA_STATIC << endl;
            if(xvasp.AVASP_flag_RUN_RELAX_STATIC_BANDS) {
                aflowin << "[VASP_KPOINTS_FILE]BANDS_LATTICE=" << xvasp.AVASP_path_BANDS << endl;
                aflowin << "[VASP_KPOINTS_FILE]BANDS_GRID=" << xvasp.AVASP_value_BANDS_GRID << endl;
            }
        }
    }
    // HTQC do AUTO
    if(xvasp.AVASP_prototype_mode==LIBRARY_MODE_HTQC || xvasp.AVASP_prototype_mode==LIBRARY_MODE_HTQC_ICSD || xvasp.AVASP_prototype_mode==LIBRARY_MODE_LIB3) {
        if(!xvasp.aopts.flag("AFLOWIN_FLAG::KPPRA_STATIC")) xvasp.AVASP_value_KPPRA_STATIC=DEFAULT_KPPRA_STATIC; 
        xvasp.AVASP_STATIC_KSCHEME=DEFAULT_STATIC_KSCHEME;
        xvasp.AVASP_path_BANDS=DEFAULT_BANDS_LATTICE;
        if(!xvasp.aopts.flag("AFLOWIN_FLAG::BANDS_GRID")) xvasp.AVASP_value_BANDS_GRID=DEFAULT_BANDS_GRID;
        if(xvasp.str.species.size()==1) if(!xvasp.aopts.flag("AFLOWIN_FLAG::KPPRA_STATIC")) xvasp.AVASP_value_KPPRA_STATIC=DEFAULT_UNARY_KPPRA_STATIC;
        if(xvasp.str.species.size()==1) if(!xvasp.aopts.flag("AFLOWIN_FLAG::BANDS_GRID")) xvasp.AVASP_value_BANDS_GRID=DEFAULT_UNARY_BANDS_GRID;
        aflowin << "[VASP_KPOINTS_FILE]STATIC_KSCHEME=" << xvasp.AVASP_STATIC_KSCHEME << " " << endl;
        aflowin << "[VASP_KPOINTS_FILE]STATIC_KPPRA=" << xvasp.AVASP_value_KPPRA_STATIC << endl;
        aflowin << "[VASP_KPOINTS_FILE]BANDS_LATTICE=" << xvasp.AVASP_path_BANDS << endl;
        aflowin << "[VASP_KPOINTS_FILE]BANDS_GRID=" << xvasp.AVASP_value_BANDS_GRID << endl;
    }
    // PROTOTYPE check the stuff out
    if(xvasp.AVASP_prototype_mode==LIBRARY_MODE_PROTOTYPE) {
        cerr << "DEBUG - AVASP calculating parameters for bands" << endl;
        if(!xvasp.aopts.flag("AFLOWIN_FLAG::KPPRA_STATIC")) xvasp.AVASP_value_KPPRA_STATIC=DEFAULT_KPPRA_STATIC;                           
        xvasp.AVASP_STATIC_KSCHEME=DEFAULT_STATIC_KSCHEME;
        // kpoints for the brillouin zone
        xvasp.str.GetLatticeType(); // takes care of everything
        xvasp.AVASP_path_BANDS=xvasp.str.bravais_lattice_variation_type;  // ICSD BASTARDS
        cerr << "DEBUG - AVASP xvasp.AVASP_path_BANDS=" << xvasp.AVASP_path_BANDS << endl;
        if(!xvasp.aopts.flag("AFLOWIN_FLAG::BANDS_GRID")) xvasp.AVASP_value_BANDS_GRID=DEFAULT_BANDS_GRID;                               
        if(xvasp.AVASP_path_BANDS=="HEX" || xvasp.AVASP_path_BANDS=="FCC") {
            xvasp.AVASP_KSCHEME="G";           // HEXAGONAL/FCC SYSTEMS GET GAMMA KPOINTS GRID
            xvasp.AVASP_STATIC_KSCHEME="G";    // HEXAGONAL/FCC SYSTEMS GET GAMMA KPOINTS GRID
        }
        if(xvasp.AVASP_flag_RUN_RELAX_STATIC || xvasp.AVASP_flag_RUN_RELAX_STATIC_BANDS || xvasp.AVASP_prototype_mode==LIBRARY_MODE_PROTOTYPE) {
            aflowin << "[VASP_KPOINTS_FILE]STATIC_KSCHEME=" << xvasp.AVASP_STATIC_KSCHEME << " " << endl;
            aflowin << "[VASP_KPOINTS_FILE]STATIC_KPPRA=" << xvasp.AVASP_value_KPPRA_STATIC << endl;
            if(xvasp.AVASP_flag_RUN_RELAX_STATIC_BANDS || xvasp.AVASP_prototype_mode==LIBRARY_MODE_PROTOTYPE) {
                aflowin << "[VASP_KPOINTS_FILE]BANDS_LATTICE=" << xvasp.AVASP_path_BANDS << endl;
                aflowin << "[VASP_KPOINTS_FILE]BANDS_GRID=" << xvasp.AVASP_value_BANDS_GRID << endl;
            }
        }
    }

    aflowin << AFLOWIN_SEPARATION_LINE << endl; // [AFLOW] **************************************************
    // STRUCTURE VASP
    if(xvasp.aopts.flag("AFLOWIN_FLAG::VASP")) {
        aflowin << AFLOWIN_SEPARATION_LINE << endl; // [AFLOW] **************************************************
        aflowin << "[VASP_POSCAR_MODE_EXPLICIT]START " << endl;
        aflowin << xvasp.str;
        aflowin << "[VASP_POSCAR_MODE_EXPLICIT]STOP " << endl;
    }
    // STRUCTURE QUANTUM ESPRESSO
    //  xvasp.aopts.flag("AFLOWIN_FLAG::QE",FALSE); // APL patch
    if(xvasp.aopts.flag("AFLOWIN_FLAG::QE")) {
        aflowin_qe << AFLOWIN_SEPARATION_LINE << endl; // [AFLOW] **************************************************
        aflowin_qe << "[QE_GEOM_MODE_EXPLICIT]START " << endl;
        xstructure _str(xvasp.str);_str.xstructure2qe();
        aflowin_qe << _str;
        aflowin_qe << "[QE_GEOM_MODE_EXPLICIT]STOP " << endl;
        //   aflowin_qe << AFLOWIN_SEPARATION_LINE << endl; // [AFLOW] **************************************************
    }
    // STRUCTURE ABINIT
    //  xvasp.aopts.flag("AFLOWIN_FLAG::ABINIT",FALSE); // APL patch
    if(xvasp.aopts.flag("AFLOWIN_FLAG::ABINIT")) {
        aflowin_abinit << AFLOWIN_SEPARATION_LINE << endl; // [AFLOW] **************************************************
        aflowin_abinit << "[ABINIT_GEOM_MODE_EXPLICIT]START " << endl;
        xstructure _str(xvasp.str);_str.xstructure2abinit();
        aflowin_abinit << _str;
        aflowin_abinit << "[ABINIT_GEOM_MODE_EXPLICIT]STOP " << endl;
        //   aflowin_abinit << AFLOWIN_SEPARATION_LINE << endl; // [AFLOW] **************************************************
    }
    // STRUCTURE AIMS
    //  xvasp.aopts.flag("AFLOWIN_FLAG::AIMS",FALSE); // APL patch
    if(xvasp.aopts.flag("AFLOWIN_FLAG::AIMS")) {
        aflowin_aims << AFLOWIN_SEPARATION_LINE << endl; // [AFLOW] **************************************************
        aflowin_aims << "[AIMS_GEOM_MODE_EXPLICIT]START " << endl;
        xstructure _str(xvasp.str);_str.xstructure2aims();
        aflowin_aims << _str;
        aflowin_aims << "[AIMS_GEOM_MODE_EXPLICIT]STOP " << endl;
        //   aflowin_aims << AFLOWIN_SEPARATION_LINE << endl; // [AFLOW] **************************************************
    }
    // POTENTIAL
    aflowin << AFLOWIN_SEPARATION_LINE << endl; // [AFLOW] **************************************************
    aflowin << "[VASP_POTCAR_MODE_IMPLICIT] " << endl;
    if(LDEBUG) cerr << "DEBUG - AVASP_MakeSingleAFLOWIN [13a]" << endl;
    for(uint i=0;i<xvasp.str.species.size();i++) {
        if(xvasp.aopts.flag("FLAG::AVASP_AUTO_PSEUDOPOTENTIALS")==FALSE) {
            if(xvasp.str.species.at(i)!="" && xvasp.str.species.at(i)!="X") {// need because some times we have the "" around
                aflowin << "[VASP_POTCAR_FILE]" << xvasp.AVASP_potential << "/" << xvasp.str.species_pp.at(i) << endl;
                if(LDEBUG) cerr << soliloquy << " AUTO_PSEUDOPOTENTIALS=FALSE - " << xvasp.AVASP_potential << "/" << xvasp.str.species_pp.at(i) << endl;
            }
        } else {
            if(xvasp.str.species.at(i)!="" && xvasp.str.species.at(i)!="X") // need because some times we have the "" around
                aflowin << "[VASP_POTCAR_FILE]" << KBIN::VASP_PseudoPotential_CleanName(xvasp.str.species.at(i)) << endl;
        }
    }
    if(LDEBUG) cerr << "DEBUG - AVASP_MakeSingleAFLOWIN [13b]" << endl;
    if(xvasp.aopts.flag("FLAG::AVASP_AUTO_PSEUDOPOTENTIALS")==TRUE) {
        aflowin << "[AFLOW] " << xvasp.AVASP_potential << ": ";
        for(uint i=0;i<xvasp.str.species_pp.size();i++) aflowin << xvasp.str.species_pp.at(i) << " ";
        aflowin << endl;
    }
    if(1) {
        aflowin << "[AFLOW] COMPOSITION_PP=|";  // COMPOSITION_PP
        for(uint i=0;i<xvasp.str.species_pp.size();i++)
            aflowin << xvasp.str.species_pp.at(i) << xvasp.str.num_each_type.at(i) << "|";
        aflowin << endl;
        aflowin << "[AFLOW] COMPOSITION=|"; // COMPOSITION
        for(uint i=0;i<xvasp.str.species_pp.size();i++) 
            aflowin << KBIN::VASP_PseudoPotential_CleanName(xvasp.str.species.at(i)) << xvasp.str.num_each_type.at(i) << "|";
        aflowin << endl;
        aflowin << "[AFLOW] VOLUME(A^3)=|"; // VOLUME
        for(uint i=0;i<xvasp.str.species_pp.size();i++)
            aflowin << aurostd::utype2string<double>(xvasp.str.species_volume.at(i),9) << "|";
        aflowin << endl;
        aflowin << "[AFLOW] MASS(amu)=|"; // MASS
        for(uint i=0;i<xvasp.str.species_pp.size();i++)
            aflowin << aurostd::utype2string<double>(xvasp.str.species_mass.at(i)/AMU2KILOGRAM,9) << "|";
        aflowin << endl;
    }
    if(LDEBUG) cerr << "DEBUG - AVASP_MakeSingleAFLOWIN [13c]" << endl;

    // add quantum espresso LOGICS (if any)
    if(xvasp.aopts.flag("AFLOWIN_FLAG::QE")) {
        aflowin << AFLOWIN_SEPARATION_LINE << endl; // [AFLOW] **************************************************
        aflowin << "[AFLOW] Quantum Espresso information (experimental)" << endl;
        //    aflowin << AFLOWIN_SEPARATION_LINE << endl; // [AFLOW] **************************************************
        aflowin << aflowin_qe.str();
    }
    // add abinit LOGICS (if any)
    if(xvasp.aopts.flag("AFLOWIN_FLAG::ABINIT")) {
        aflowin << AFLOWIN_SEPARATION_LINE << endl; // [AFLOW] **************************************************
        aflowin << "[AFLOW] Abinit information (experimental)" << endl;
        //    aflowin << AFLOWIN_SEPARATION_LINE << endl; // [AFLOW] **************************************************
        aflowin << aflowin_abinit.str();
    }
    // add aims LOGICS (if any)
    if(xvasp.aopts.flag("AFLOWIN_FLAG::AIMS")) {
        aflowin << AFLOWIN_SEPARATION_LINE << endl; // [AFLOW] **************************************************
        aflowin << "[AFLOW] Aims information (experimental)" << endl;
        //    aflowin << AFLOWIN_SEPARATION_LINE << endl; // [AFLOW] **************************************************
        aflowin << aflowin_aims.str();
    }
    // now close
    aflowin << AFLOWIN_SEPARATION_LINE << endl; // [AFLOW] **************************************************
    aflowin << "[AFLOW] Aflow automatically generated (aflow_avasp.cpp) " << endl;
    aflowin << AFLOWIN_SEPARATION_LINE << endl; // [AFLOW] **************************************************
    aflowin << "[AFLOW] AFLOW V(" << string(AFLOW_VERSION) << ") in " << directory << endl;
    aflowin << AFLOWIN_SEPARATION_LINE << endl; // [AFLOW] **************************************************

    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " " << "[14]" << endl;

    string msg;
    if(xvasp.AVASP_prototype_mode==LIBRARY_MODE_LIB3)
        flag_WRITE=AVASP_CheckHeuslersSanvito(xvasp,msg);

    // NO WRITE
    if(flag_WRITE==FALSE && 0) {
        cout << "AFLOW V(" << string(AFLOW_VERSION) << ") skipping " << directory;
        cout << " [" << msg << "] *************** ";
        cout << endl;cout.flush();
    }

    // YES WRITE
    if(flag_WRITE==TRUE) {
        // DONE and write
        if(1) {
            if(directory=="") {
                cerr << "AVASP_MakeSingleAFLOWIN empty directory" << endl;exit(0);
            }
            if(pthread>=0)  pthread_mutex_lock(&mutex_AVASP);
            cout << "AFLOW V(" << string(AFLOW_VERSION) << ") creating " << directory;
            if(xvasp.aopts.flag("FLAG::AVASP_LDAU1")) cout << " [LDAU1]";
            if(xvasp.aopts.flag("FLAG::AVASP_LDAU2")) cout << " [LDAU2]";
            if(fixed_lattice_vasp) cout << " [fixed_lattice_vasp=" << xvasp.AVASP_label << "]";
            if(pthread>=0) cout << "   - pthread=" << pthread;
            cout << endl;cout.flush();
            if(pthread>=0)  pthread_mutex_unlock(&mutex_AVASP);
        }
        // cout << aflowin.str() << endl; exit(0);
        // MAKE Directory and WRITE _AFLOWIN_
        Krun=aurostd::DirectoryMake(directory);
        if(Krun==FALSE) return Krun;
        // CHMOD Directory 777
        Krun=aurostd::DirectoryChmod("777",directory);
        if(Krun==FALSE) return Krun;
        // Write _AFLOWIN_
        Krun=aurostd::stringstream2file(aflowin,string(directory+"/"+_AFLOWIN_));
        if(Krun==FALSE) return Krun;
        // CHMOD a+w _AFLOWIN_
        Krun=aurostd::ChmodFile("a+w",string(directory+"/"+_AFLOWIN_));
        if(Krun==FALSE) return Krun;
        //  cout << str;
        if(LDEBUG) cerr << "DEBUG - AVASP_MakeSingleAFLOWIN [15]" << endl;
    }
    if(flag_WRITE==FALSE&&flag_PRINT==TRUE) {
        cout << aflowin.str() << endl;cout.flush();
    }

    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " " << "[16]" << endl;

    // make _aflowin
    _aflowin.clear();_aflowin.str(std::string());
    _aflowin << aflowin.str();


    return TRUE;
}

bool AVASP_MakeSingleAFLOWIN(_xvasp& xvasp_in,bool flag_WRITE,int pthread,bool flag_PRINT) {
    stringstream _aflowin;
    return AVASP_MakeSingleAFLOWIN(xvasp_in,_aflowin,flag_WRITE,pthread,flag_PRINT);
}

bool AVASP_MakeSingleAFLOWIN(_xvasp& xvasp_in,int pthread,bool flag_PRINT) {
    stringstream _aflowin;
    return AVASP_MakeSingleAFLOWIN(xvasp_in,_aflowin,TRUE,pthread,flag_PRINT);
}

// ***************************************************************************
bool AVASP_MakeSinglePOSCAR(_xvasp& xvasp) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    // checks if everything is different
    ofstream oaus;
    stringstream poscar;
    xstructure str;
    string directory,system;
    bool Krun=TRUE; //,swapspeciesAB=FALSE,swapspeciesBC=FALSE,swapspeciesCA=FALSE;
    poscar.clear();poscar.str(std::string());
    // check if 2 are equal
    for(uint isp1=0;isp1<xvasp.str.species.size();isp1++)
        for(uint isp2=isp1+1;isp2<xvasp.str.species.size();isp2++)
            if(KBIN::VASP_PseudoPotential_CleanName(xvasp.str.species.at(isp1))==KBIN::VASP_PseudoPotential_CleanName(xvasp.str.species.at(isp2))) {
                ostringstream aus;
                if(LDEBUG) aus << "EEEEE PROTO_Generation_Binary_Poscar error xvasp.str.species.at(isp1)==xvasp.str.species.at(isp2) " << xvasp.str.species.at(isp1) << endl;
                if(LDEBUG) aurostd::PrintErrorStream(oaus,aus,XHOST.QUIET);
                Krun=FALSE;
                return Krun; // empty
            }

    //  cerr << xvasp.AVASP_volumeC; // keep it bisy
    if(xvasp.AVASP_alpha_fix==TRUE) {
        xvasp.str.SpeciesPutAlphabetic(); // should put alphabetic also the PSEUDOPOTENTIALS
    }

    if(xvasp.AVASP_prototype_from_library_) {
        directory=xvasp.AVASP_dirbase+"/"+aurostd::vectorstring2string(xvasp.str.species_pp)+"/"+xvasp.AVASP_label;
        //  cerr << "DEBUG - AVASP_MakeSingleAFLOWIN: directory(4)=" << directory << endl;
        system=aurostd::vectorstring2string(xvasp.str.species_pp)+"."+xvasp.AVASP_label;
    }
    // MAKE Directory
    Krun=aurostd::DirectoryMake(directory);
    if(Krun==FALSE) return Krun;
    // STRUCTURE
    poscar << xvasp.str;
    // DONE and write
    cout << "AFLOW V(" << string(AFLOW_VERSION) << ") creating " << directory << endl;
    aurostd::stringstream2file(poscar,string(directory+"/POSCAR"));
    //  cout << str;
    return TRUE;
}

// ***************************************************************************
bool AVASP_DefaultValuesBinary_AFLOWIN(_xvasp &xvasp) {
    xvasp.clear();
    xvasp.AVASP_dirbase="./AFLOWDATA";
    xvasp.AVASP_prototype_mode=LIBRARY_MODE_HTQC;
    xvasp.AVASP_alpha_fix=FALSE;
    xvasp.AVASP_prototype_from_library_=TRUE;
    xvasp.AVASP_directory_from_library_=TRUE;
    xvasp.aopts.flag("AFLOWIN_FLAG::NBANDS",TRUE);
    xvasp.aopts.flag("FLAG::AVASP_SKIP_NOMIX",TRUE);
    xvasp.aopts.flag("FLAG::AVASP_CHGCAR",TRUE);
    xvasp.aopts.flag("FLAG::AVASP_SPIN",TRUE);
    xvasp.aopts.flag("FLAG::AVASP_SPIN_REMOVE_RELAX_1",DEFAULT_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_1);
    xvasp.aopts.flag("FLAG::AVASP_SPIN_REMOVE_RELAX_2",DEFAULT_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_2);
    xvasp.aopts.flag("FLAG::AVASP_BADER",FALSE);
    xvasp.aopts.flag("FLAG::AVASP_ELF",FALSE);
    xvasp.aopts.flag("FLAG::AVASP_AUTO_MAGMOM",FALSE);
    xvasp.AVASP_value_KPPRA=DEFAULT_KPPRA;
    if(xvasp.aopts.flag("AFLOWIN_FLAG::KPPRA")) xvasp.AVASP_value_KPPRA=aurostd::string2utype<uint>(xvasp.aopts.getattachedscheme("AFLOWIN_FLAG::KPPRA"));
    xvasp.AVASP_KSCHEME=DEFAULT_KSCHEME;
    xvasp.AVASP_value_KPPRA_STATIC=DEFAULT_KPPRA_STATIC;
    if(xvasp.aopts.flag("AFLOWIN_FLAG::KPPRA_STATIC")) xvasp.AVASP_value_KPPRA_STATIC=aurostd::string2utype<uint>(xvasp.aopts.getattachedscheme("AFLOWIN_FLAG::KPPRA_STATIC"));
    xvasp.AVASP_STATIC_KSCHEME=DEFAULT_STATIC_KSCHEME;
    xvasp.AVASP_value_NSW=-1;
    xvasp.aopts.flag("FLAG::PRECISION_SET",TRUE);
    xvasp.AVASP_flag_PRECISION_scheme=DEFAULT_VASP_FORCE_OPTION_PREC_SCHEME;
    xvasp.aopts.flag("FLAG::ALGO_SET",TRUE);
    xvasp.AVASP_flag_ALGO_scheme="F";
    // [OBSOLETE] xvasp.aopts.flag("FLAG::METAGGA_SET",FALSE);
    // [OBSOLETE] xvasp.AVASP_flag_METAGGA_scheme=DEFAULT_VASP_FORCE_OPTION_METAGGA_SCHEME;
    // [OBSOLETE] xvasp.aopts.flag("FLAG::IVDW_SET",FALSE);
    // [OBSOLETE] xvasp.AVASP_flag_IVDW_scheme=DEFAULT_VASP_FORCE_OPTION_IVDW_SCHEME;
    xvasp.aopts.flag("AVASP_flag_CONVERT_UNIT_CELL_STANDARD_PRIMITIVE",TRUE);
    xvasp.aopts.flag("AVASP_flag_CONVERT_UNIT_CELL_MINKOWSKI",TRUE); 
    xvasp.aopts.flag("FLAG::VOLUME_PRESERVED",FALSE);
    xvasp.aopts.flag("FLAG::EXTRA_INCAR",FALSE);
    xvasp.AVASP_EXTRA_INCAR.clear();
    xvasp.AVASP_volume_in=-1.0;
    xvasp.AVASP_flag_RUN_RELAX=TRUE;
    xvasp.AVASP_path_BANDS="";         // DEFAULT VALUES
    xvasp.AVASP_value_BANDS_GRID=DEFAULT_BANDS_GRID;  // DEFAULT VALUES
    if(xvasp.aopts.flag("AFLOWIN_FLAG::BANDS_GRID")) xvasp.AVASP_value_BANDS_GRID=aurostd::string2utype<uint>(xvasp.aopts.getattachedscheme("AFLOWIN_FLAG::BANDS_GRID"));
    return TRUE;
}

// ***************************************************************************
bool AVASP_CheckHeuslersSanvito(_xvasp& xvasp, string& msg) { // STEFANO SANVITO WITH HEUSLERS
    return TRUE;
    bool flag_WRITE=TRUE;
    if(xvasp.AVASP_label=="T0001.A2BC" || xvasp.AVASP_label=="T0001.AB2C" || xvasp.AVASP_label=="T0001.ABC2" ||
            xvasp.AVASP_label=="T0002.A2BC" || xvasp.AVASP_label=="T0002.AB2C" || xvasp.AVASP_label=="T0002.ABC2" ||
            xvasp.AVASP_label=="T0003.ABC" || xvasp.AVASP_label=="T0003.BCA" || xvasp.AVASP_label=="T0003.CAB") {
        uint occ=0;
        for(uint i=0;i<xvasp.str.species_pp.size();i++) {
            if(KBIN::VASP_PseudoPotential_CleanName(xvasp.str.species.at(i))=="B")  {occ++;msg+="B";if(xvasp.str.num_each_type.at(i)==2) {msg+="2";flag_WRITE=FALSE;}}
            if(KBIN::VASP_PseudoPotential_CleanName(xvasp.str.species.at(i))=="Al") {occ++;msg+="Al";if(xvasp.str.num_each_type.at(i)==2) {msg+="2";flag_WRITE=FALSE;}}
            if(KBIN::VASP_PseudoPotential_CleanName(xvasp.str.species.at(i))=="Si") {occ++;msg+="Si";if(xvasp.str.num_each_type.at(i)==2) {msg+="2";flag_WRITE=FALSE;}}
            if(KBIN::VASP_PseudoPotential_CleanName(xvasp.str.species.at(i))=="Ga") {occ++;msg+="Ga";if(xvasp.str.num_each_type.at(i)==2) {msg+="2";flag_WRITE=FALSE;}}
            if(KBIN::VASP_PseudoPotential_CleanName(xvasp.str.species.at(i))=="Ge") {occ++;msg+="Ge";if(xvasp.str.num_each_type.at(i)==2) {msg+="2";flag_WRITE=FALSE;}}
            if(KBIN::VASP_PseudoPotential_CleanName(xvasp.str.species.at(i))=="As") {occ++;msg+="As";if(xvasp.str.num_each_type.at(i)==2) {msg+="2";flag_WRITE=FALSE;}}
            if(KBIN::VASP_PseudoPotential_CleanName(xvasp.str.species.at(i))=="In") {occ++;msg+="In";if(xvasp.str.num_each_type.at(i)==2) {msg+="2";flag_WRITE=FALSE;}}
            if(KBIN::VASP_PseudoPotential_CleanName(xvasp.str.species.at(i))=="Sn") {occ++;msg+="Sn";if(xvasp.str.num_each_type.at(i)==2) {msg+="2";flag_WRITE=FALSE;}}
            if(KBIN::VASP_PseudoPotential_CleanName(xvasp.str.species.at(i))=="Sb") {occ++;msg+="Sb";if(xvasp.str.num_each_type.at(i)==2) {msg+="2";flag_WRITE=FALSE;}}
            if(KBIN::VASP_PseudoPotential_CleanName(xvasp.str.species.at(i))=="Pb") {occ++;msg+="Pb";if(xvasp.str.num_each_type.at(i)==2) {msg+="2";flag_WRITE=FALSE;}}
            if(KBIN::VASP_PseudoPotential_CleanName(xvasp.str.species.at(i))=="Bi") {occ++;msg+="Bi";if(xvasp.str.num_each_type.at(i)==2) {msg+="2";flag_WRITE=FALSE;}}
            if(KBIN::VASP_PseudoPotential_CleanName(xvasp.str.species.at(i))=="Mg") {occ++;msg+="Mg";if(xvasp.str.num_each_type.at(i)==2) {msg+="2";flag_WRITE=FALSE;}}
            if(KBIN::VASP_PseudoPotential_CleanName(xvasp.str.species.at(i))=="P")  {occ++;msg+="P";if(xvasp.str.num_each_type.at(i)==2) {msg+="2";flag_WRITE=FALSE;}}
            if(KBIN::VASP_PseudoPotential_CleanName(xvasp.str.species.at(i))=="Se") {occ++;msg+="Se";if(xvasp.str.num_each_type.at(i)==2) {msg+="2";flag_WRITE=FALSE;}}
            if(KBIN::VASP_PseudoPotential_CleanName(xvasp.str.species.at(i))=="Te") {occ++;msg+="Te";if(xvasp.str.num_each_type.at(i)==2) {msg+="2";flag_WRITE=FALSE;}}
        }
        if(occ==2 || occ==3) {flag_WRITE=FALSE;}
    }
    return flag_WRITE;
}


// ***************************************************************************
// typedef struct {
//   deque<_xvasp>  *dxvasp;     // CONTENT
//   bool     flag_WRITE;  // CONTENT
//   int      itbusy;      // FOR AVASP_AFLOWIN (ALL)
//   int      ITHREAD;     // FOR MULTISH_TIMESHARING
//   int      THREADS_MAX; // FOR MULTISH_TIMESHARING
// } _threaded_AVASP_AFLOWIN_params;

// pthread_mutex_t mutex_AVASP=PTHREAD_MUTEX_INITIALIZER;

#define AVASP_AFLOWIN_MAX_JOBS 256

void *_threaded_AVASP_AFLOWIN_CONCURRENT(void *ptr) {
    bool FRONT=TRUE;
    _threaded_AVASP_AFLOWIN_params* pparams=(_threaded_AVASP_AFLOWIN_params*) ptr;
    _xvasp xvasp;
    AFLOW_PTHREADS::vpthread_busy[pparams->itbusy]=TRUE;
    AFLOW_PTHREADS::RUNNING++;
    while((*pparams->dxvasp).size()>0) {
        pthread_mutex_lock(&mutex_AVASP);
        if(FRONT)  {xvasp=(*pparams->dxvasp).at(0);(*pparams->dxvasp).pop_front();}  // from the front
        if(!FRONT) {xvasp=(*pparams->dxvasp).at((*pparams->dxvasp).size()-1);(*pparams->dxvasp).pop_back();}  // from the back
        pthread_mutex_unlock(&mutex_AVASP);
        AVASP_MakeSingleAFLOWIN(xvasp,pparams->flag_WRITE,pparams->itbusy);
    }
    AFLOW_PTHREADS::vpthread_busy[pparams->itbusy]=FALSE;
    AFLOW_PTHREADS::RUNNING--;
    return NULL;
}

void *_threaded_AVASP_AFLOWIN_SEQUENTIAL(void *ptr) {
    _threaded_AVASP_AFLOWIN_params* pparams=(_threaded_AVASP_AFLOWIN_params*) ptr;
    _xvasp xvasp;
    AFLOW_PTHREADS::vpthread_busy[pparams->itbusy]=TRUE;
    AFLOW_PTHREADS::RUNNING++;
    for(uint ithread=(pparams->ITHREAD);ithread<(*pparams->dxvasp).size();ithread+=(pparams->THREADS_MAX)) {
        xvasp=(*pparams->dxvasp).at(ithread);
        AVASP_MakeSingleAFLOWIN(xvasp,pparams->flag_WRITE,pparams->itbusy);  
    }
    AFLOW_PTHREADS::vpthread_busy[pparams->itbusy]=FALSE;
    AFLOW_PTHREADS::RUNNING--;
    return NULL;
}

bool AVASP_MakePrototype_AFLOWIN_LOOP(deque<_xvasp>& dxvasp,bool flag_WRITE) {
    AFLOW_PTHREADS::Clean_Threads();      
    int NUM_THREADS=min(8,XHOST.CPU_Cores);  
    deque<_xvasp> dxvasp_local;
    cerr << dxvasp.size() << endl;
    uint _MOD=AVASP_AFLOWIN_MAX_JOBS;
    _threaded_AVASP_AFLOWIN_params params[MAX_ALLOCATABLE_PTHREADS];                 // prepare
    for(uint i=0;i<dxvasp.size();i+=_MOD) {   //   cerr << i << endl;
        dxvasp_local.clear();
        for(uint j=i;j<i+_MOD&&j<dxvasp.size();j++)
            dxvasp_local.push_back(dxvasp.at(j));
        // now run with dxvasp_local_local
        bool SEQUENTIAL=FALSE;
        if(dxvasp_local.size()>5000) SEQUENTIAL=TRUE;  // should be faster for longer deques but I dont know the threshold (SC)
        if((int) dxvasp_local.size()<=NUM_THREADS) NUM_THREADS=(uint) dxvasp_local.size(); // SAFETY    
        if(NUM_THREADS<=1)                                                     // run singular
            for(uint k=0;k<dxvasp_local.size();k++)                              // run singular
                AVASP_MakeSingleAFLOWIN(dxvasp_local.at(k),flag_WRITE);            // run singular
        if(NUM_THREADS>=2) {                                                   // multithread
            AFLOW_PTHREADS::FLAG=TRUE;AFLOW_PTHREADS::MAX_PTHREADS=NUM_THREADS;             // prepare
            if(AFLOW_PTHREADS::MAX_PTHREADS>MAX_ALLOCATABLE_PTHREADS) AFLOW_PTHREADS::MAX_PTHREADS=MAX_ALLOCATABLE_PTHREADS; // check max
            AFLOW_PTHREADS::Clean_Threads();                                            // clean threads
            //  _threaded_AVASP_AFLOWIN_params params[MAX_ALLOCATABLE_PTHREADS];             // prepare
            for(int ithread=0;ithread<AFLOW_PTHREADS::MAX_PTHREADS;ithread++) {            // prepare loop
                params[ithread].dxvasp=&dxvasp_local;                              // prepare params
                params[ithread].flag_WRITE=flag_WRITE;                             // prepare params
                params[ithread].ITHREAD=ithread;                                   // prepare params
                params[ithread].THREADS_MAX=AFLOW_PTHREADS::MAX_PTHREADS;                    // prepare params
                params[ithread].itbusy=ithread;                                    // prepare params
            }                                                
            for(int ithread=0;ithread<AFLOW_PTHREADS::MAX_PTHREADS;ithread++) {            // run threads
                if( SEQUENTIAL) AFLOW_PTHREADS::viret[ithread]=pthread_create(&AFLOW_PTHREADS::vpthread[ithread],NULL,_threaded_AVASP_AFLOWIN_SEQUENTIAL,(void*)&params[ithread]); // run threads
                if(!SEQUENTIAL) AFLOW_PTHREADS::viret[ithread]=pthread_create(&AFLOW_PTHREADS::vpthread[ithread],NULL,_threaded_AVASP_AFLOWIN_CONCURRENT,(void*)&params[ithread]); // run threads
            }
            for(int ithread=0;ithread<AFLOW_PTHREADS::MAX_PTHREADS;ithread++)              // flush
                pthread_join(AFLOW_PTHREADS::vpthread[ithread],NULL);               // flush
            //  cerr << "FLUSHED" << endl;
        }
    }
    dxvasp.clear();
    return TRUE;
}

// ***************************************************************************
// struct _AVASP_PROTO {
//   vector<string> ucell;
//   string label;
//   deque<int> vkppra;
//   vector<double> vpressure;
//   aurostd::xoption vparams;
// };

bool AVASP_MakePrototype_AFLOWIN(_AVASP_PROTO *PARAMS) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string soliloquy="AVASP_MakePrototype_AFLOWIN():";
    stringstream message;
    if(LDEBUG) cerr << " " << soliloquy << " BEGIN" << endl; 

    _xvasp xvasp;
    AVASP_DefaultValuesBinary_AFLOWIN(xvasp);
    xvasp.AVASP_prototype_mode=LIBRARY_MODE_HTQC;
    xvasp.AVASP_parameters=PARAMS->vparams.getattachedscheme("AFLOWIN_FLAG::PARAMS"); // DX 1/18/18 - added parameters to AVASP

    if(PARAMS->vparams.flag("AFLOWIN_FLAG::AUTOLDAU")) xvasp.aopts.flag("FLAG::AVASP_FORCE_LDAU",TRUE);
    if(PARAMS->vparams.flag("AFLOWIN_FLAG::AUTONOLDAU")) xvasp.aopts.flag("FLAG::AVASP_FORCE_NOLDAU",TRUE);
    if(PARAMS->vparams.flag("AFLOWIN_FLAG::VASP")) xvasp.aopts.flag("AFLOWIN_FLAG::VASP",TRUE);
    if(PARAMS->vparams.flag("AFLOWIN_FLAG::ABINIT")) xvasp.aopts.flag("AFLOWIN_FLAG::ABINIT",TRUE);
    if(PARAMS->vparams.flag("AFLOWIN_FLAG::AIMS")) xvasp.aopts.flag("AFLOWIN_FLAG::AIMS",TRUE);
    if(PARAMS->vparams.flag("AFLOWIN_FLAG::QE")) xvasp.aopts.flag("AFLOWIN_FLAG::QE",TRUE);
    if(PARAMS->vparams.flag("AFLOWIN_FLAG::HTQC_ICSD")) xvasp.AVASP_prototype_mode=LIBRARY_MODE_HTQC_ICSD;
    if(PARAMS->vparams.flag("AFLOWIN_FLAG::MODULE")) {xvasp.aopts.flag("AFLOWIN_FLAG::MODULE",TRUE);xvasp.aopts.push_attached("AFLOWIN_FLAG::MODULE",PARAMS->vparams.getattachedscheme("AFLOWIN_FLAG::MODULE"));}  // CO 180214
    if(PARAMS->vparams.flag("AFLOWIN_FLAG::APL_SUPERCELL")) {xvasp.aopts.flag("AFLOWIN_FLAG::APL_SUPERCELL",TRUE);xvasp.aopts.push_attached("AFLOWIN_FLAG::APL_SUPERCELL",PARAMS->vparams.getattachedscheme("AFLOWIN_FLAG::APL_SUPERCELL"));}  // CO 180214
    if(PARAMS->vparams.flag("AFLOWIN_FLAG::POTIM")) {xvasp.aopts.flag("AFLOWIN_FLAG::POTIM",TRUE);xvasp.aopts.push_attached("AFLOWIN_FLAG::POTIM",PARAMS->vparams.getattachedscheme("AFLOWIN_FLAG::POTIM"));}
    if(PARAMS->vparams.flag("AFLOWIN_FLAG::PRECISION")) {xvasp.aopts.flag("AFLOWIN_FLAG::PRECISION",TRUE);xvasp.aopts.push_attached("AFLOWIN_FLAG::PRECISION",PARAMS->vparams.getattachedscheme("AFLOWIN_FLAG::PRECISION"));}
    if(PARAMS->vparams.flag("AFLOWIN_FLAG::ALGORITHM")) {xvasp.aopts.flag("AFLOWIN_FLAG::ALGORITHM",TRUE);xvasp.aopts.push_attached("AFLOWIN_FLAG::ALGORITHM",PARAMS->vparams.getattachedscheme("AFLOWIN_FLAG::ALGORITHM"));}
    if(PARAMS->vparams.flag("AFLOWIN_FLAG::METAGGA")) {xvasp.aopts.flag("AFLOWIN_FLAG::METAGGA",TRUE);xvasp.aopts.push_attached("AFLOWIN_FLAG::METAGGA",PARAMS->vparams.getattachedscheme("AFLOWIN_FLAG::METAGGA"));}  
    if(PARAMS->vparams.flag("AFLOWIN_FLAG::IVDW")) {xvasp.aopts.flag("AFLOWIN_FLAG::IVDW",TRUE);xvasp.aopts.push_attached("AFLOWIN_FLAG::IVDW",PARAMS->vparams.getattachedscheme("AFLOWIN_FLAG::IVDW"));}  
    if(PARAMS->vparams.flag("AFLOWIN_FLAG::RELAX_TYPE")) {xvasp.aopts.flag("AFLOWIN_FLAG::RELAX_TYPE",TRUE);xvasp.aopts.push_attached("AFLOWIN_FLAG::RELAX_TYPE",PARAMS->vparams.getattachedscheme("AFLOWIN_FLAG::RELAX_TYPE"));}  // CO 180214
    if(PARAMS->vparams.flag("AFLOWIN_FLAG::RELAX_MODE")) {xvasp.aopts.flag("AFLOWIN_FLAG::RELAX_MODE",TRUE);xvasp.aopts.push_attached("AFLOWIN_FLAG::RELAX_MODE",PARAMS->vparams.getattachedscheme("AFLOWIN_FLAG::RELAX_MODE"));}
    if(PARAMS->vparams.flag("AFLOWIN_FLAG::TYPE")) {xvasp.aopts.flag("AFLOWIN_FLAG::TYPE",TRUE);xvasp.aopts.push_attached("AFLOWIN_FLAG::TYPE",PARAMS->vparams.getattachedscheme("AFLOWIN_FLAG::TYPE"));}
    if(PARAMS->vparams.flag("AFLOWIN_FLAG::CONVERT_UNIT_CELL")) {xvasp.aopts.flag("AFLOWIN_FLAG::CONVERT_UNIT_CELL",TRUE);xvasp.aopts.push_attached("AFLOWIN_FLAG::CONVERT_UNIT_CELL",PARAMS->vparams.getattachedscheme("AFLOWIN_FLAG::CONVERT_UNIT_CELL"));}
    if(PARAMS->vparams.flag("AFLOWIN_FLAG::VOLUME_PLUS_EQUAL")) {xvasp.aopts.flag("AFLOWIN_FLAG::VOLUME_PLUS_EQUAL",TRUE);xvasp.aopts.push_attached("AFLOWIN_FLAG::VOLUME_PLUS_EQUAL",PARAMS->vparams.getattachedscheme("AFLOWIN_FLAG::VOLUME_PLUS_EQUAL"));}
    if(PARAMS->vparams.flag("AFLOWIN_FLAG::VOLUME_MULTIPLY_EQUAL")) {xvasp.aopts.flag("AFLOWIN_FLAG::VOLUME_MULTIPLY_EQUAL",TRUE);xvasp.aopts.push_attached("AFLOWIN_FLAG::VOLUME_MULTIPLY_EQUAL",PARAMS->vparams.getattachedscheme("AFLOWIN_FLAG::VOLUME_MULTIPLY_EQUAL"));}
    if(PARAMS->vparams.flag("AFLOWIN_FLAG::NO_VOLUME_ADJUSTMENT")) xvasp.aopts.flag("AFLOWIN_FLAG::NO_VOLUME_ADJUSTMENT",TRUE); // CO 180214
    if(PARAMS->vparams.flag("AFLOWIN_FLAG::EDIFFG")) {xvasp.aopts.flag("AFLOWIN_FLAG::EDIFFG",TRUE);xvasp.aopts.push_attached("AFLOWIN_FLAG::EDIFFG",PARAMS->vparams.getattachedscheme("AFLOWIN_FLAG::EDIFFG"));}
    if(PARAMS->vparams.flag("AFLOWIN_FLAG::KPPRA")) {xvasp.aopts.flag("AFLOWIN_FLAG::KPPRA",TRUE);xvasp.aopts.push_attached("AFLOWIN_FLAG::KPPRA",PARAMS->vparams.getattachedscheme("AFLOWIN_FLAG::KPPRA"));}
    if(PARAMS->vparams.flag("AFLOWIN_FLAG::KPPRA_STATIC")) {xvasp.aopts.flag("AFLOWIN_FLAG::KPPRA_STATIC",TRUE);xvasp.aopts.push_attached("AFLOWIN_FLAG::KPPRA_STATIC",PARAMS->vparams.getattachedscheme("AFLOWIN_FLAG::KPPRA_STATIC"));}
    if(PARAMS->vparams.flag("AFLOWIN_FLAG::BANDS_GRID")) {xvasp.aopts.flag("AFLOWIN_FLAG::BANDS_GRID",TRUE);xvasp.aopts.push_attached("AFLOWIN_FLAG::BANDS_GRID",PARAMS->vparams.getattachedscheme("AFLOWIN_FLAG::BANDS_GRID"));}
    if(PARAMS->vparams.flag("AFLOWIN_FLAG::ENMAX_MULTIPLY")) {xvasp.aopts.flag("AFLOWIN_FLAG::ENMAX_MULTIPLY",TRUE);xvasp.aopts.push_attached("AFLOWIN_FLAG::ENMAX_MULTIPLY",PARAMS->vparams.getattachedscheme("AFLOWIN_FLAG::ENMAX_MULTIPLY"));}
    if(PARAMS->vparams.flag("AFLOWIN_FLAG::PARAMS")) xvasp.aopts.flag("AFLOWIN_FLAG::PARAMS",TRUE); // DX 1/23/18 - added ANRL parameters to PARAMS

    string string_POTENTIAL=PARAMS->vparams.getattachedscheme("AFLOWIN_STRING::POTENTIAL");
    if(LDEBUG)  cerr << "DEBUG - " << soliloquy << " string_POTENTIAL=" << string_POTENTIAL << endl;

    if(PARAMS->vpressure.size()==0) {
        xvasp.aopts.flag("AFLOWIN_FLAG::PSTRESS",FALSE);
        // [OBSOLETE]    xvasp.AVASP_value_PSTRESS=0.0;
        PARAMS->vpressure.push_back(0.0);
    } else {
        xvasp.aopts.flag("AFLOWIN_FLAG::PSTRESS",TRUE);
    }

    if(PARAMS->vkppra.size()==0) {
        xvasp.aopts.flag("FLAG::AVASP_KPPRA",FALSE);
        xvasp.AVASP_value_KPPRA=DEFAULT_KPPRA;
        if(xvasp.aopts.flag("AFLOWIN_FLAG::KPPRA")) xvasp.AVASP_value_KPPRA=aurostd::string2utype<int>(xvasp.aopts.getattachedscheme("AFLOWIN_FLAG::KPPRA"));    // KPPRA
        //   PARAMS->vkppra.push_back(xvasp.AVASP_value_KPPRA=DEFAULT_KPPRA);
    } else {
        xvasp.aopts.flag("FLAG::AVASP_KPPRA",TRUE);
        //   xvasp.AVASP_value_KPPRA=PARAMS->vkppra.front();  // once forever
    }

    // loop on pressures

    // recompile and ahave fun

    //  DEBUG=TRUE;

    if(PARAMS->vparams.flag("AFLOWIN_FLAG::BANDS")) {xvasp.AVASP_flag_RUN_RELAX_STATIC_BANDS=TRUE;xvasp.AVASP_flag_RUN_RELAX=FALSE;}

    bool HTQC=FALSE,flag_LIB2U=FALSE,flag_LIB3=FALSE,flag_LIB2=FALSE;
    if(!PARAMS->vparams.flag("AFLOWIN_FLAG::MISSING")) xvasp.AVASP_aflowin_only_if_missing=FALSE;
    if(PARAMS->vparams.flag("AFLOWIN_FLAG::MISSING")) xvasp.AVASP_aflowin_only_if_missing=TRUE;
    if(PARAMS->vparams.flag("AFLOWIN_FLAG::NEGLECT_NOMIX")==TRUE) xvasp.aopts.flag("FLAG::AVASP_SKIP_NOMIX",FALSE);


    uint nspecies=0;
    uint nspeciesHTQC=aflowlib::PrototypeLibrariesSpeciesNumber(PARAMS->ucell.at(0));

    bool alphabetic=TRUE;
    if(!PARAMS->vparams.flag("AFLOWIN_FLAG::HTQC_ICSD")) alphabetic=RequestedAlphabeticLabeling(PARAMS->ucell.at(0));
    if(PARAMS->vparams.flag("AFLOWIN_FLAG::HTQC_ICSD")) alphabetic=TRUE;

    if(nspeciesHTQC==3) {xvasp.AVASP_flag_RUN_RELAX_STATIC_BANDS=TRUE;xvasp.AVASP_flag_RUN_RELAX=FALSE;} // force bands for ternary
    if(nspeciesHTQC==4) {xvasp.AVASP_flag_RUN_RELAX_STATIC_BANDS=TRUE;xvasp.AVASP_flag_RUN_RELAX=FALSE;} // force bands for quaternary


    // I`m not sure shifts for quaternary is correct... I need to try

    if(PARAMS->ucell.size()==2) nspecies=1; else nspecies=2;

    //  if(PARAMS->ucell.size()==2 || PARAMS->ucell.size()==1+nspeciesHTQC) AUTOVOLS=TRUE;
    if(nspeciesHTQC==3) nspecies=3;
    if(nspeciesHTQC==4) nspecies=4;
    if(PARAMS->vparams.flag("AFLOWIN_FLAG::HTQC_ICSD")) {nspecies=nspeciesHTQC=PARAMS->ucell.size()-1;}

    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " labels=" << PARAMS->ucell.at(0) << endl;
    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " nspecies=" << nspecies << endl;
    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " nspeciesHTQC=" << nspeciesHTQC << endl;
    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " alphabetic=" << alphabetic << endl;
    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " PARAMS->ucell.at(0)=" << PARAMS->ucell.at(0) << endl;
    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " PARAMS->ucell.size()=" << PARAMS->ucell.size() << endl;
    for(uint i=0;i<PARAMS->ucell.size();i++)
        if(LDEBUG) cerr << "DEBUG - " << soliloquy << " PARAMS->ucell.at(" << i << ")=" << PARAMS->ucell.at(i) << endl;
    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " PARAMS->vparams.flag(\"AFLOWIN_FLAG::HTQC_ICSD\")=" << PARAMS->vparams.flag("AFLOWIN_FLAG::HTQC_ICSD") << endl;
    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " 1+nspecies=" << 1+nspecies << endl;
    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " 1+nspecies+1=" << 1+nspecies+1 << endl;
    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " 2+nspecies=" << 2+nspecies << endl;
    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " xvasp.AVASP_potential=" << xvasp.AVASP_potential << endl;
    if(LDEBUG) {cerr << "DEBUG - " << soliloquy << " ";for(uint i=0;i<PARAMS->ucell.size();i++) cerr << PARAMS->ucell.at(i) << " ";cerr << endl;}

    vector<string> label;
    vector<vector<string> > specieX;
    vector<vector<double> > volumeX;
    vector<vector<double> > massX;
    vector<string> specieX_raw;
    string label_raw;
    // some wrap up and definitions
    if(PARAMS->ucell.size()!=1+nspecies && PARAMS->ucell.size()!=1+nspecies+1 && PARAMS->ucell.size()!=1+nspecies*2) {
        cerr << soliloquy << " need to specify more arguments... exiting" << endl;
        cerr << " aflow --aflow_proto[=]label*:specieA*[:specieB*][:volumeA*[:volumeB*] | :volume]" << endl;
        exit(0);
    }

    if(PARAMS->ucell.size()==1) {
        vector<string> austokens;
        aurostd::RemoveSubString(PARAMS->ucell.at(0),"/"+_AFLOWIN_);
        aurostd::RemoveSubString(PARAMS->ucell.at(0),"/"+_AFLOWLOCK_);
        aurostd::string2tokens(PARAMS->ucell.at(0),austokens,"/");
        KBIN::VASP_SplitAlloySpecies(austokens.at(0),specieX_raw);
        nspecies=specieX_raw.size()+1;  // backward
        if(nspeciesHTQC==3) nspecies=3;       // new
        if(nspeciesHTQC==4) nspecies=4;       // new
        specieX.clear();volumeX.clear();massX.clear();
        for(uint ispecies=0;ispecies<nspecies;ispecies++) {
            //  cerr << specieX_raw.at(ispecies) << endl;
            specieX.push_back(*(new vector<string>(0))); specieX.at(specieX.size()-1).clear();
            volumeX.push_back(*(new vector<double>(0))); volumeX.at(volumeX.size()-1).clear();
            massX.push_back(*(new vector<double>(0))); massX.at(massX.size()-1).clear();
        }
    }

    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " [x1] " << endl;

    // some wrap up and definitions
    if(PARAMS->ucell.size()==1+nspecies || PARAMS->ucell.size()==1+nspecies+1 || PARAMS->ucell.size()==1+nspecies*2) {
        for(uint ispecies=0;ispecies<nspecies;ispecies++) {
            specieX.push_back(*(new vector<string>(0))); specieX.at(specieX.size()-1).clear();
            specieX_raw.push_back("");
            volumeX.push_back(*(new vector<double>(0))); volumeX.at(volumeX.size()-1).clear();
            massX.push_back(*(new vector<double>(0))); massX.at(massX.size()-1).clear();
        }

        if(PARAMS->vparams.flag("AFLOWIN_FLAG::REVERSE")==FALSE) {
            label_raw=PARAMS->ucell.at(0);
            for(uint isp=1;isp<=nspecies;isp++)
                if(nspecies>=isp) specieX_raw.at(isp-1)=PARAMS->ucell.at(0+isp);
        }
        if(PARAMS->vparams.flag("AFLOWIN_FLAG::REVERSE")==TRUE) {
            specieX_raw.at(0)=PARAMS->ucell.at(0);
            if(nspecies==1) {label_raw=PARAMS->ucell.at(1);}
            if(nspecies>=2) {specieX_raw.at(1)=PARAMS->ucell.at(1);label_raw=PARAMS->ucell.at(2);}  // for backward compatibility
        }
    }

    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " [x2] " << endl;

    // SEARCH FOR LABELS
    label_raw=AVASP_Shortcuts_for_Binaries(label_raw);
    label_raw=AVASP_Shortcuts_for_Ternaries(label_raw);
    //  cerr << "label_raw=" << label_raw << endl;
    // SEARCH FOR SPECIES
    if(PARAMS->ucell.size()==1+nspecies)
    {
        if(LDEBUG) cerr << "DEBUG - " << soliloquy << " [x2a] " << endl;
        bool present_all,present_here;
        // LIB3
        present_all=TRUE;present_here=FALSE;
        for(uint isp=0;isp<specieX_raw.size();isp++) {
            present_here=aurostd::substring2bool(specieX_raw.at(isp),"LIB3");
            if(present_here) aurostd::StringSubst(specieX_raw.at(isp),"LIB3",SPECIE_RAW_LIB3);
            present_all=present_all && present_here;
        }
        if(present_all) {
            cout << "LIB3 Compatibility mode" << endl;
            xvasp.AVASP_prototype_mode=LIBRARY_MODE_LIB3;
            flag_LIB3=TRUE;xvasp.AVASP_dirbase="./LIB3/LIB";xvasp.AVASP_libbase=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB3)+"/LIB";
        }
        // LIB2U
        present_all=TRUE;present_here=FALSE;
        for(uint isp=0;isp<specieX_raw.size();isp++) {
            present_here=aurostd::substring2bool(specieX_raw.at(isp),"LIB2U");
            if(present_here) aurostd::StringSubst(specieX_raw.at(isp),"LIB2U",SPECIE_RAW_LIB2U);
            present_all=present_all && present_here;
        }
        if(present_all) {
            cout << "LIB2U Compatibility mode" << endl;
            xvasp.AVASP_prototype_mode=LIBRARY_MODE_HTQC;
            flag_LIB2U=TRUE;xvasp.AVASP_dirbase="./LIB2/LIB";xvasp.AVASP_libbase=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB2)+"/LIB";
        }
        // LIB2
        present_all=TRUE;present_here=FALSE;
        for(uint isp=0;isp<specieX_raw.size();isp++) {
            present_here=aurostd::substring2bool(specieX_raw.at(isp),"LIB2");
            if(present_here) aurostd::StringSubst(specieX_raw.at(isp),"LIB2",SPECIE_RAW_LIB2);
            present_all=present_all && present_here;
        }
        if(present_all) {
            cout << "LIB2 Compatibility mode" << endl;
            xvasp.AVASP_prototype_mode=LIBRARY_MODE_HTQC;
            flag_LIB2=TRUE;xvasp.AVASP_dirbase="./LIB2/LIB";xvasp.AVASP_libbase=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB2)+"/LIB";
        }
    }

    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " [x3] " << endl;

    aurostd::string2tokens(label_raw,label,",");
    for(uint ispecies=0;ispecies<nspecies;ispecies++) {
        //  cerr << "specieX_raw.at(ispecies)=" << specieX_raw.at(ispecies) << endl;
        aurostd::string2tokens(specieX_raw.at(ispecies),specieX.at(ispecies),",");
    }
    ostringstream aus;
    if(nspecies>=2) if(specieX.at(0).size()>1 || specieX.at(1).size()>1) HTQC=TRUE;
    if(HTQC) xvasp.AVASP_alpha_fix=TRUE;
    if(flag_LIB3) xvasp.AVASP_alpha_fix=TRUE;   // some fix for historic reasons
    if(flag_LIB2U) xvasp.AVASP_alpha_fix=TRUE;   // some fix for historic reasons
    if(flag_LIB2) xvasp.AVASP_alpha_fix=TRUE;   // some fix for historic reasons
    // get the xvasp.AVASP_potential
    if(PARAMS->ucell.size()==1+nspecies) {  // USE AUTOMATIC VOLUMES
        if(LDEBUG) cerr << "DEBUG - " << soliloquy << " [x3b] - USE AUTOMATIC VOLUMES" << endl;
        for(uint ispecies=0;ispecies<nspecies;ispecies++)
            for(uint i=0;i<specieX.at(ispecies).size();i++) {
                volumeX.at(ispecies).push_back(VolumeSpecie(specieX.at(ispecies).at(i)));
                massX.at(ispecies).push_back(MassSpecie(specieX.at(ispecies).at(i)));
                //	cerr << 1.0e25*MassSpecie(specieX.at(ispecies).at(i)) << endl;
            }
    }
    double vol; //CO 180705
    if(PARAMS->ucell.size()==1+nspecies+1) {  // USE ONE VOLUME FITS ALL
        if(LDEBUG) cerr << "DEBUG - " << soliloquy << " [x3d] - USE ONE VOLUME FITS ALL" << endl;
        for(uint ispecies=0;ispecies<nspecies;ispecies++) {
            vector<string> tokens;tokens.clear();
            aurostd::string2tokens(PARAMS->ucell.at(1+nspecies+1-1),tokens,",");
            for(uint i=0;i<specieX.at(ispecies).size();i++) {
                if(!aurostd::isfloat(tokens.at(tokens.size()-1))){ //CO 180729 - check for string stupidity
                    message << "Invalid volume specification (params[" << tokens.size()-1 << "]=" << tokens.at(tokens.size()-1) << "), must be float input";
                    throw aurostd::xerror(soliloquy,message,_INPUT_ILLEGAL_);
                }
                vol=aurostd::string2utype<double>(tokens.at(tokens.size()-1));  //CO 180705
                if(vol==0.0){ //CO 180705 - check for volume stupidity
                    message << "Invalid volume specification (params[" << 1+nspecies+1-1 << "]=" << tokens.at(tokens.size()-1) << "), must be >0";
                    throw aurostd::xerror(soliloquy,message,_INPUT_ILLEGAL_);
                }
                volumeX.at(ispecies).push_back(vol);
                //volumeX.at(ispecies).push_back(aurostd::string2utype<double>(tokens.at(tokens.size()-1)));
                massX.at(ispecies).push_back(MassSpecie(specieX.at(ispecies).at(i)));
                //	cerr << 1.0e25*MassSpecie(specieX.at(ispecies).at(i)) << endl;
            }
        }
    }
    if(PARAMS->ucell.size()==1+nspecies*2) {
        if(LDEBUG) cerr << "DEBUG - " << soliloquy << " [x3e] - EXTRACT VOLUMES FROM STRING" << endl;
        for(uint ispecies=0;ispecies<nspecies;ispecies++) {
            vector<string> tokens;tokens.clear();
            aurostd::string2tokens(PARAMS->ucell.at(1+nspecies+ispecies),tokens,",");
            for(uint i=0;i<tokens.size();i++) {
                if(!aurostd::isfloat(tokens.at(i))){ //CO 180729 - check for string stupidity
                    message << "Invalid volume specification (params[" << i << "]=" << tokens.at(i) << "), must be float input";
                    throw aurostd::xerror(soliloquy,message,_INPUT_ILLEGAL_);
                }
                vol=aurostd::string2utype<double>(tokens.at(i));
                if(vol==0.0){ //CO 180705 - check for stupidity
                    message << "Invalid volume specification (params[" << i << "]=" << tokens.at(i) << "), must be >0";
                    throw aurostd::xerror(soliloquy,message,_INPUT_ILLEGAL_);
                }
                volumeX.at(ispecies).push_back(vol);
                //volumeX.at(ispecies).push_back(aurostd::string2utype<double>(tokens.at(i)));
                massX.at(ispecies).push_back(MassSpecie(specieX.at(ispecies).at(i)));
                //	cerr << 1.0e25*MassSpecie(specieX.at(ispecies).at(i)) << endl;
            }
        }
        for(uint ispecies=0;ispecies<nspecies;ispecies++) {
            if(specieX.at(ispecies).size() != volumeX.at(ispecies).size()) {
                aus << "EEEEE  PROTO_Generation_Binary_Aflowin error \"specieX.at(" << ispecies << ").size()!=volumeX.at(" << ispecies << ").size()\" " << specieX.at(ispecies).size() << "," << volumeX.at(ispecies).size() << endl;
                aurostd::PrintErrorStream(cerr,aus,XHOST.QUIET);
                return FALSE;
            }
            if(specieX.at(ispecies).size() != massX.at(ispecies).size()) {
                aus << "EEEEE  PROTO_Generation_Binary_Aflowin error \"specieX.at(" << ispecies << ").size()!=massX.at(" << ispecies << ").size()\" " << specieX.at(ispecies).size() << "," << massX.at(ispecies).size() << endl;
                aurostd::PrintErrorStream(cerr,aus,XHOST.QUIET);
                return FALSE;
            }
        }
    }
    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " [0] xvasp.AVASP_potential=" << xvasp.AVASP_potential << endl;
    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " [0] string_POTENTIAL=" << string_POTENTIAL << endl;

    aurostd::StringSubst(string_POTENTIAL,  // TEMPORARY PATCH Sun Apr 14 22:05:01 EDT 2013
            _AVASP_PSEUDOPOTENTIAL_AUTO_+_AVASP_PSEUDOPOTENTIAL_DELIMITER_+_AVASP_PSEUDOPOTENTIAL_POTENTIAL_COMPLETE_, // TEMPORARY PATCH Sun Apr 14 22:05:01 EDT 2013
            DEFAULT_VASP_POTCAR_DIR_POTPAW_PBE+_AVASP_PSEUDOPOTENTIAL_DELIMITER_+_AVASP_PSEUDOPOTENTIAL_POTENTIAL_COMPLETE_); // TEMPORARY PATCH Sun Apr 14 22:05:01 EDT 2013
    xvasp.AVASP_potential=DEFAULT_VASP_POTCAR_DIR_POTPAW_PBE;
    if(aurostd::substring2bool(string_POTENTIAL,_AVASP_PSEUDOPOTENTIAL_AUTO_)) {
        xvasp.AVASP_potential=DEFAULT_VASP_POTCAR_DIR_POTPAW_PBE;
        if(DEFAULT_VASP_PSEUDOPOTENTIAL_TYPE=="pot_LDA") xvasp.AVASP_potential=DEFAULT_VASP_POTCAR_DIR_POT_LDA;
        if(DEFAULT_VASP_PSEUDOPOTENTIAL_TYPE=="pot_GGA") xvasp.AVASP_potential=DEFAULT_VASP_POTCAR_DIR_POT_GGA;
        if(DEFAULT_VASP_PSEUDOPOTENTIAL_TYPE=="pot_PBE") xvasp.AVASP_potential=DEFAULT_VASP_POTCAR_DIR_POT_PBE;
        if(DEFAULT_VASP_PSEUDOPOTENTIAL_TYPE=="potpaw_LDA") xvasp.AVASP_potential=DEFAULT_VASP_POTCAR_DIR_POTPAW_LDA;
        if(DEFAULT_VASP_PSEUDOPOTENTIAL_TYPE=="potpaw_GGA") xvasp.AVASP_potential=DEFAULT_VASP_POTCAR_DIR_POTPAW_GGA;
        if(DEFAULT_VASP_PSEUDOPOTENTIAL_TYPE=="potpaw_PBE") xvasp.AVASP_potential=DEFAULT_VASP_POTCAR_DIR_POTPAW_PBE;
        if(DEFAULT_VASP_PSEUDOPOTENTIAL_TYPE=="potpaw_LDA_KIN") xvasp.AVASP_potential=DEFAULT_VASP_POTCAR_DIR_POTPAW_LDA_KIN;
        if(DEFAULT_VASP_PSEUDOPOTENTIAL_TYPE=="potpaw_PBE_KIN") xvasp.AVASP_potential=DEFAULT_VASP_POTCAR_DIR_POTPAW_PBE_KIN;
        //   cerr << "MAKING AUTO" << endl; exit(0);
    }
    //  if((!aurostd::substring2bool(string_POTENTIAL,DEFAULT_VASP_POTCAR_DIR_POTPAW_PBE) &&
    if((string_POTENTIAL!=DEFAULT_VASP_POTCAR_DIR_POTPAW_PBE &&
                !aurostd::substring2bool(string_POTENTIAL,_AVASP_PSEUDOPOTENTIAL_AUTO_)) ||
            aurostd::substring2bool(string_POTENTIAL,_AVASP_PSEUDOPOTENTIAL_DELIMITER_)) {
        xvasp.AVASP_potential=string_POTENTIAL;
        if(LDEBUG) cerr << "DEBUG - " << soliloquy << " [1] xvasp.AVASP_potential=" << xvasp.AVASP_potential << endl;
        if(LDEBUG) cerr << "DEBUG - " << soliloquy << " [1] string_POTENTIAL=" << string_POTENTIAL << endl;

        if(aurostd::substring2bool(string_POTENTIAL,_AVASP_PSEUDOPOTENTIAL_DELIMITER_)) {
            vector<string> tokens_string_POTENTIAL;
            cerr << "PARAMS->vparams.getattachedscheme(\"AFLOWIN_STRING::POTENTIAL\")=" << string_POTENTIAL << endl;
            aurostd::string2tokens(string_POTENTIAL,tokens_string_POTENTIAL,_AVASP_PSEUDOPOTENTIAL_DELIMITER_);
            if(tokens_string_POTENTIAL.size()>=1)
                xvasp.AVASP_potential=tokens_string_POTENTIAL.at(0);
            xvasp.POTCAR_TYPE_DATE_PRINT_flag=FALSE;
            if(tokens_string_POTENTIAL.size()>=2)
                if(tokens_string_POTENTIAL.at(1)==_AVASP_PSEUDOPOTENTIAL_POTENTIAL_COMPLETE_)
                    xvasp.POTCAR_TYPE_DATE_PRINT_flag=TRUE;
        }
    } 

    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " [X] xvasp.AVASP_potential=" << xvasp.AVASP_potential << endl;
    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " [X] string_POTENTIAL=" << string_POTENTIAL << endl;
    // for(uint i=0;i<PARAMS->ucell.size();i++) cerr << PARAMS->ucell.at(i) << endl;
    // exit(0);


    if(xvasp.AVASP_potential!=DEFAULT_VASP_POTCAR_DIR_POT_LDA &&
            xvasp.AVASP_potential!=DEFAULT_VASP_POTCAR_DIR_POT_GGA &&
            xvasp.AVASP_potential!=DEFAULT_VASP_POTCAR_DIR_POT_PBE &&
            xvasp.AVASP_potential!=DEFAULT_VASP_POTCAR_DIR_POTPAW_LDA &&
            xvasp.AVASP_potential!=DEFAULT_VASP_POTCAR_DIR_POTPAW_GGA &&
            xvasp.AVASP_potential!=DEFAULT_VASP_POTCAR_DIR_POTPAW_PBE &&
            xvasp.AVASP_potential!=DEFAULT_VASP_POTCAR_DIR_POTPAW_LDA_KIN &&
            xvasp.AVASP_potential!=DEFAULT_VASP_POTCAR_DIR_POTPAW_PBE_KIN) {
        // bad potential
        aus << "EEEEE  PROTO_Generation_Binary_Aflowin potential is not VASP style: " << xvasp.AVASP_potential << endl;
        aurostd::PrintErrorStream(cerr,aus,XHOST.QUIET);
        return FALSE;
    }
    // good potential

    //  bool CHECK_LDAU=TRUE;
    //  if(CHECK_LDAU) AVASP_ADD_LDAU(xvasp);
    //  cerr << xvasp.aopts.flag("FLAG::AVASP_LDAU1") << endl;
    //  cerr << xvasp.aopts.flag("FLAG::AVASP_LDAU2") << endl;

    if(LDEBUG) cerr << "DEBUG - " << soliloquy << " [x4] " << endl;

    // ***************************************************************************
    deque<_xvasp> dxvasp;
    bool MULTITHREAD=FALSE; //TRUE;
    // ***************************************************************************
    // 1 SPECIES
    if(nspecies==1) {
        if(LDEBUG) cerr <<"DEBUG - " << soliloquy << " [5] 1species" << endl;
        uint ii[nspecies];
        for(ii[0]=0;ii[0]<specieX.at(0).size();ii[0]++) {
            for(uint k=0;k<label.size();k++) {
                for(uint ip=0;ip<PARAMS->vpressure.size();ip++) {
                    _xvasp xaus(xvasp);
                    xaus.AVASP_label=label.at(k);
                    xaus.str.species.clear();xaus.str.species_pp.clear();xaus.str.species_pp_type.clear();xaus.str.species_pp_version.clear();xaus.str.species_pp_ZVAL.clear();xaus.str.species_pp_vLDAU.clear();xaus.str.species_volume.clear();xaus.str.species_mass.clear();
                    for(uint i=0;i<nspecies;i++) {xaus.str.species.push_back(specieX.at(i).at(ii[i]));xaus.str.species_volume.push_back(volumeX.at(i).at(ii[i]));xaus.str.species_mass.push_back(massX.at(i).at(ii[i]));}
                    xaus.str.species_pp=xaus.str.species;
                    for(uint i=0;i<nspecies;i++) xaus.str.species_pp_type.push_back(xaus.AVASP_potential);
                    for(uint i=0;i<nspecies;i++) xaus.str.species_pp_version.push_back("");
                    for(uint i=0;i<nspecies;i++) xaus.str.species_pp_ZVAL.push_back(0.0);
                    for(uint i=0;i<nspecies;i++) xaus.str.species_pp_vLDAU.push_back(deque<double>());
                    if(!aurostd::substring2bool(string_POTENTIAL,_AVASP_PSEUDOPOTENTIAL_AUTO_)) xaus.aopts.flag("FLAG::AVASP_AUTO_PSEUDOPOTENTIALS",FALSE);
                    if(xaus.aopts.flag("FLAG::AVASP_FORCE_LDAU")) AVASP_ADD_LDAU(xaus);
                    if(xaus.aopts.flag("FLAG::AVASP_FORCE_NOLDAU")) AVASP_REMOVE_LDAU(xaus);	
                    xaus.aopts.push_attached("AFLOWIN_FLAG::PSTRESS",aurostd::utype2string(PARAMS->vpressure.at(ip)));
                    if(PARAMS->vpressure.size()==1 && abs(aurostd::string2utype<double>(xaus.aopts.getattachedscheme("AFLOWIN_FLAG::PSTRESS")))<0.00001) {
                        xaus.aopts.flag("AFLOWIN_FLAG::PSTRESS",FALSE);xaus.aopts.push_attached("AFLOWIN_FLAG::PSTRESS","0.0");};
                    if(!PARAMS->vparams.flag("AFLOWIN_FLAG::LIST")) {
                        if(!MULTITHREAD) AVASP_MakeSingleAFLOWIN(xaus,!PARAMS->vparams.flag("AFLOWIN_FLAG::STDOUT"));   // SINGLE_THREAD
                        if(MULTITHREAD) {dxvasp.push_back(xaus);if(dxvasp.size()==AVASP_AFLOWIN_MAX_JOBS) AVASP_MakePrototype_AFLOWIN_LOOP(dxvasp,!PARAMS->vparams.flag("AFLOWIN_FLAG::STDOUT"));}
                    } else {
                        if(xaus.aopts.flag("AFLOWIN_FLAG::PSTRESS")==FALSE) {
                            cout << "aflow " << PARAMS->vparams.getattachedscheme("AFLOWIN_FLAG::LIST_VCMD") << "--aflow_proto=" << xaus.AVASP_label; 
                            for(uint i=0;i<nspecies;i++) cout << ":" << KBIN::VASP_PseudoPotential_CleanName(xaus.str.species.at(i));	
                            cout << endl;
                        } else {
                            cout << "aflow " << PARAMS->vparams.getattachedscheme("AFLOWIN_FLAG::LIST_VCMD") << "--pressure=" << xaus.aopts.getattachedscheme("AFLOWIN_FLAG::PSTRESS") << " --aflow_proto=" << xaus.AVASP_label; 
                            for(uint i=0;i<nspecies;i++) cout << ":" << KBIN::VASP_PseudoPotential_CleanName(xaus.str.species.at(i)); 
                            cout << endl;
                        }
                    }
                }
            }
        }
        if(MULTITHREAD) AVASP_MakePrototype_AFLOWIN_LOOP(dxvasp,!PARAMS->vparams.flag("AFLOWIN_FLAG::STDOUT")); // FLUSH OUT
    }
    // ***************************************************************************
    // 2 SPECIES
    if(nspecies==2) {
        if(LDEBUG) cerr <<"DEBUG - " << soliloquy << " [5] 2species" << endl;
        uint ii[nspecies];
        for(ii[0]=0;ii[0]<specieX.at(0).size();ii[0]++) {
            if(LDEBUG) cerr <<"DEBUG - " << soliloquy << " [5a] specieX.at(0).size()=" << specieX.at(0).size() << "  -  specieX.at(0).at(ii[0])=" << specieX.at(0).at(ii[0]) << endl;
            for(ii[1]=0;ii[1]<specieX.at(1).size();ii[1]++) {
                if(LDEBUG) cerr <<"DEBUG - " << soliloquy << " [5b] specieX.at(1).size()=" << specieX.at(1).size() << "  -  specieX.at(1).at(ii[1])=" << specieX.at(1).at(ii[1]) << endl;
                if(specieX.at(0).at(ii[0])<specieX.at(1).at(ii[1])) {
                    for(uint k=0;k<label.size();k++) {
                        for(uint ip=0;ip<PARAMS->vpressure.size();ip++) {
                            _xvasp xaus(xvasp);
                            xaus.AVASP_label=label.at(k);
                            xaus.str.species.clear();
                            xaus.str.species_pp.clear();
                            xaus.str.species_pp_type.clear();
                            xaus.str.species_pp_version.clear();
                            xaus.str.species_pp_ZVAL.clear();
                            xaus.str.species_pp_vLDAU.clear();
                            xaus.str.species_volume.clear();
                            xaus.str.species_mass.clear();
                            for(uint i=0;i<nspecies;i++) {
                                xaus.str.species.push_back(specieX.at(i).at(ii[i]));
                                xaus.str.species_volume.push_back(volumeX.at(i).at(ii[i]));
                                xaus.str.species_mass.push_back(massX.at(i).at(ii[i]));
                            }
                            xaus.str.species_pp=xaus.str.species;
                            for(uint i=0;i<nspecies;i++) xaus.str.species_pp_type.push_back(xaus.AVASP_potential);
                            for(uint i=0;i<nspecies;i++) xaus.str.species_pp_version.push_back("");
                            for(uint i=0;i<nspecies;i++) xaus.str.species_pp_ZVAL.push_back(0.0);
                            for(uint i=0;i<nspecies;i++) xaus.str.species_pp_vLDAU.push_back(deque<double>());
                            if(!aurostd::substring2bool(string_POTENTIAL,_AVASP_PSEUDOPOTENTIAL_AUTO_)) xaus.aopts.flag("FLAG::AVASP_AUTO_PSEUDOPOTENTIALS",FALSE);
                            if(xaus.aopts.flag("FLAG::AVASP_FORCE_LDAU")) AVASP_ADD_LDAU(xaus);
                            if(xaus.aopts.flag("FLAG::AVASP_FORCE_NOLDAU")) AVASP_REMOVE_LDAU(xaus);
                            xaus.aopts.push_attached("AFLOWIN_FLAG::PSTRESS",aurostd::utype2string(PARAMS->vpressure.at(ip)));
                            if(PARAMS->vpressure.size()==1 && abs(aurostd::string2utype<double>(xaus.aopts.getattachedscheme("AFLOWIN_FLAG::PSTRESS")))<0.00001) {
                                xaus.aopts.flag("AFLOWIN_FLAG::PSTRESS",FALSE);xaus.aopts.push_attached("AFLOWIN_FLAG::PSTRESS","0.0");};
                            if(!PARAMS->vparams.flag("AFLOWIN_FLAG::LIST")) {
                                if(!MULTITHREAD) AVASP_MakeSingleAFLOWIN(xaus,!PARAMS->vparams.flag("AFLOWIN_FLAG::STDOUT"));   // SINGLE_THREAD
                                if(MULTITHREAD) {dxvasp.push_back(xaus);if(dxvasp.size()==AVASP_AFLOWIN_MAX_JOBS) AVASP_MakePrototype_AFLOWIN_LOOP(dxvasp,!PARAMS->vparams.flag("AFLOWIN_FLAG::STDOUT"));}
                            } else {
                                if(xaus.aopts.flag("AFLOWIN_FLAG::PSTRESS")==FALSE) {
                                    cout << "aflow " << PARAMS->vparams.getattachedscheme("AFLOWIN_FLAG::LIST_VCMD") << "--aflow_proto=" << xaus.AVASP_label;
                                    for(uint i=0;i<nspecies;i++) cout << ":" << KBIN::VASP_PseudoPotential_CleanName(xaus.str.species.at(i)); 
                                    cout << endl;
                                } else {cout << "aflow " << PARAMS->vparams.getattachedscheme("AFLOWIN_FLAG::LIST_VCMD") << "--pressure=" << xaus.aopts.getattachedscheme("AFLOWIN_FLAG::PSTRESS") << " --aflow_proto=" << xaus.AVASP_label; 
                                    for(uint i=0;i<nspecies;i++) cout << ":" << KBIN::VASP_PseudoPotential_CleanName(xaus.str.species.at(i)); 
                                    cout << endl;
                                }
                            }
                        }
                    }
                } else { // no alphabetic
                    if(flag_LIB2U==FALSE && flag_LIB2==FALSE && flag_LIB3==FALSE) {
                        for(uint i=0;i<nspecies-1;i++) {
                            if(specieX.at(i).at(ii[i])>specieX.at(i+1).at(ii[i+1])) {
                                cerr << soliloquy << " system ";
                                for(uint j=0;j<nspecies;j++) {
                                    cerr << specieX.at(j).at(ii[j]);
                                }
                                cerr << " neglected (not alphabetic)" << endl;
                                // exit(0);
                            }
                        }	
                    } // no alphabetic
                }
                }
            }
            if(MULTITHREAD) AVASP_MakePrototype_AFLOWIN_LOOP(dxvasp,!PARAMS->vparams.flag("AFLOWIN_FLAG::STDOUT")); // FLUSH OUT
        }
        // **********************************************************n*****************
        // 3 SPECIES
        if(nspecies==3) {
            //   LDEBUG=TRUE;
            if(LDEBUG) cerr <<"DEBUG - " << soliloquy << " [5] 3species" << endl;
            uint ii[nspecies];
            for(ii[0]=0;ii[0]<specieX.at(0).size();ii[0]++) {
                for(ii[1]=0;ii[1]<specieX.at(1).size();ii[1]++) {
                    for(ii[2]=0;ii[2]<specieX.at(2).size();ii[2]++) {
                        if(LDEBUG) cerr <<"DEBUG - " << soliloquy << " [5.0.3] 3species specieX.at(2).at(ii[2])=" << specieX.at(2).at(ii[2]) <<  endl;
                        if(LDEBUG) { 
                            cerr << "DEBUG - " << soliloquy << " [5.0.3] 3species = " << specieX.at(0).at(ii[0]) << " " << specieX.at(1).at(ii[1]) << " " << specieX.at(2).at(ii[2]);
                            if(specieX.at(0).at(ii[0])<specieX.at(1).at(ii[1])) cerr << " 0<1 ";
                            if(specieX.at(1).at(ii[1])<specieX.at(2).at(ii[2])) cerr << " 1<2 ";
                            cerr << endl;
                        }
                        if(specieX.at(0).at(ii[0])<specieX.at(1).at(ii[1]) && specieX.at(1).at(ii[1])<specieX.at(2).at(ii[2])) { // check if identical is done later
                            if(LDEBUG) cerr <<"DEBUG - " << soliloquy << " [5.0.4] 3species " << specieX.at(0).at(ii[0]) << specieX.at(1).at(ii[1]) << specieX.at(2).at(ii[2]) <<  endl;
                            for(uint k=0;k<label.size();k++) {
                                if(LDEBUG) cerr <<"DEBUG - " << soliloquy << " [5.1] 3species label.at(" << k << ")=" << label.at(k) << endl;
                                for(uint ip=0;ip<PARAMS->vpressure.size();ip++) {
                                    if(LDEBUG) cerr <<"DEBUG - " << soliloquy << " [5.1] 3species PARAMS->vpressure.at(" << ip << ")=" << PARAMS->vpressure.at(ip) << endl;
                                    _xvasp xaus(xvasp);
                                    xaus.AVASP_label=label.at(k);
                                    xaus.str.species.clear();xaus.str.species_pp.clear();xaus.str.species_pp_type.clear();xaus.str.species_pp_version.clear();xaus.str.species_pp_ZVAL.clear();xaus.str.species_pp_vLDAU.clear();xaus.str.species_volume.clear();xaus.str.species_mass.clear();
                                    for(uint i=0;i<nspecies;i++) {xaus.str.species.push_back(specieX.at(i).at(ii[i]));xaus.str.species_volume.push_back(volumeX.at(i).at(ii[i]));xaus.str.species_mass.push_back(massX.at(i).at(ii[i]));}
                                    xaus.str.species_pp=xaus.str.species;
                                    for(uint i=0;i<nspecies;i++) xaus.str.species_pp_type.push_back(xaus.AVASP_potential);
                                    for(uint i=0;i<nspecies;i++) xaus.str.species_pp_version.push_back("");
                                    for(uint i=0;i<nspecies;i++) xaus.str.species_pp_ZVAL.push_back(0.0);
                                    for(uint i=0;i<nspecies;i++) xaus.str.species_pp_vLDAU.push_back(deque<double>());
                                    if(!aurostd::substring2bool(string_POTENTIAL,_AVASP_PSEUDOPOTENTIAL_AUTO_)) xaus.aopts.flag("FLAG::AVASP_AUTO_PSEUDOPOTENTIALS",FALSE);
                                    if(xaus.aopts.flag("FLAG::AVASP_FORCE_LDAU")) AVASP_ADD_LDAU(xaus);
                                    if(xaus.aopts.flag("FLAG::AVASP_FORCE_NOLDAU")) AVASP_REMOVE_LDAU(xaus);
                                    xaus.aopts.push_attached("AFLOWIN_FLAG::PSTRESS",aurostd::utype2string(PARAMS->vpressure.at(ip)));
                                    if(PARAMS->vpressure.size()==1 && abs(aurostd::string2utype<double>(xaus.aopts.getattachedscheme("AFLOWIN_FLAG::PSTRESS")))<0.00001) {
                                        xaus.aopts.flag("AFLOWIN_FLAG::PSTRESS",FALSE);xaus.aopts.push_attached("AFLOWIN_FLAG::PSTRESS","0.0");};
                                    if(!PARAMS->vparams.flag("AFLOWIN_FLAG::LIST")) {
                                        if(!MULTITHREAD) AVASP_MakeSingleAFLOWIN(xaus,!PARAMS->vparams.flag("AFLOWIN_FLAG::STDOUT"));   // SINGLE_THREAD
                                        if(MULTITHREAD) {dxvasp.push_back(xaus);if(dxvasp.size()==AVASP_AFLOWIN_MAX_JOBS) AVASP_MakePrototype_AFLOWIN_LOOP(dxvasp,!PARAMS->vparams.flag("AFLOWIN_FLAG::STDOUT"));}
                                    } else {
                                        if(xaus.aopts.flag("AFLOWIN_FLAG::PSTRESS")==FALSE) {
                                            cout << "aflow " << PARAMS->vparams.getattachedscheme("AFLOWIN_FLAG::LIST_VCMD") << "--aflow_proto=" << xaus.AVASP_label; 
                                            for(uint i=0;i<nspecies;i++) cout << ":" << KBIN::VASP_PseudoPotential_CleanName(xaus.str.species.at(i)); 
                                            cout << endl;
                                        } else {
                                            cout << "aflow " << PARAMS->vparams.getattachedscheme("AFLOWIN_FLAG::LIST_VCMD") << "--pressure=" << xaus.aopts.getattachedscheme("AFLOWIN_FLAG::PSTRESS") << " --aflow_proto=" << xaus.AVASP_label;
                                            for(uint i=0;i<nspecies;i++) cout << ":" << KBIN::VASP_PseudoPotential_CleanName(xaus.str.species.at(i)); 
                                            cout << endl;
                                        }
                                    }
                                }
                            }
                        } else { // no alphabetic
                            if(flag_LIB2U==FALSE && flag_LIB2==FALSE && flag_LIB3==FALSE) {
                                for(uint i=0;i<nspecies-1;i++) {
                                    if(specieX.at(i).at(ii[i])>specieX.at(i+1).at(ii[i+1])) {
                                        cerr << soliloquy << " system ";
                                        for(uint j=0;j<nspecies;j++) {
                                            cerr << specieX.at(j).at(ii[j]);
                                        }
                                        cerr << " neglected (not alphabetic)" << endl;
                                        // exit(0);
                                    }
                                }	
                            }
                        } // no alphabetic
                    }
                }
            }
            if(MULTITHREAD) AVASP_MakePrototype_AFLOWIN_LOOP(dxvasp,!PARAMS->vparams.flag("AFLOWIN_FLAG::STDOUT")); // FLUSH OUT
        }
        // ***************************************************************************
        // 4 SPECIES
        if(nspecies==4) {
            if(LDEBUG) cerr <<"DEBUG - " << soliloquy << " [5] 4species" << endl;
            uint ii[nspecies];
            for(ii[0]=0;ii[0]<specieX.at(0).size();ii[0]++) {
                for(ii[1]=0;ii[1]<specieX.at(1).size();ii[1]++) {
                    for(ii[2]=0;ii[2]<specieX.at(2).size();ii[2]++) {
                        for(ii[3]=0;ii[3]<specieX.at(3).size();ii[3]++) {
                            if(specieX.at(0).at(ii[0])<specieX.at(1).at(ii[1])  && specieX.at(1).at(ii[1])<specieX.at(2).at(ii[2])  && specieX.at(2).at(ii[2])<specieX.at(3).at(ii[3])) { // check if identical is done later
                                for(uint k=0;k<label.size();k++) {
                                    for(uint ip=0;ip<PARAMS->vpressure.size();ip++) {
                                        _xvasp xaus(xvasp);
                                        xaus.AVASP_label=label.at(k);
                                        xaus.str.species.clear();xaus.str.species_pp.clear();xaus.str.species_pp_type.clear();xaus.str.species_pp_version.clear();xaus.str.species_pp_ZVAL.clear();xaus.str.species_pp_vLDAU.clear();xaus.str.species_volume.clear();xaus.str.species_mass.clear();
                                        for(uint i=0;i<nspecies;i++) {xaus.str.species.push_back(specieX.at(i).at(ii[i]));xaus.str.species_volume.push_back(volumeX.at(i).at(ii[i]));xaus.str.species_mass.push_back(massX.at(i).at(ii[i]));}
                                        xaus.str.species_pp=xaus.str.species;
                                        for(uint i=0;i<nspecies;i++) xaus.str.species_pp_type.push_back(xaus.AVASP_potential);
                                        for(uint i=0;i<nspecies;i++) xaus.str.species_pp_version.push_back("");
                                        for(uint i=0;i<nspecies;i++) xaus.str.species_pp_ZVAL.push_back(0.0);
                                        for(uint i=0;i<nspecies;i++) xaus.str.species_pp_vLDAU.push_back(deque<double>());
                                        // check alpha
                                        if(alphabetic==TRUE) AlphabetizePrototypeLabelSpecies(xaus.str.species,xaus.str.species_pp,xaus.str.species_volume,xaus.str.species_mass,xaus.AVASP_label);
                                        if(!aurostd::substring2bool(string_POTENTIAL,_AVASP_PSEUDOPOTENTIAL_AUTO_)) xaus.aopts.flag("FLAG::AVASP_AUTO_PSEUDOPOTENTIALS",FALSE);
                                        if(1) { 
                                            xaus.aopts.flag("FLAG::AVASP_SKIP_NOMIX",FALSE);
                                            xaus.aopts.flag("FLAG::VOLUME_PRESERVED",TRUE);
                                            if(1) AVASP_ADD_LDAU(xaus);
                                            if(xaus.aopts.flag("FLAG::AVASP_FORCE_NOLDAU")) AVASP_REMOVE_LDAU(xaus);
                                        }
                                        xaus.aopts.push_attached("AFLOWIN_FLAG::PSTRESS",aurostd::utype2string(PARAMS->vpressure.at(ip)));
                                        if(PARAMS->vpressure.size()==1 && abs(aurostd::string2utype<double>(xaus.aopts.getattachedscheme("AFLOWIN_FLAG::PSTRESS")))<0.00001) {
                                            xaus.aopts.flag("AFLOWIN_FLAG::PSTRESS",FALSE);xaus.aopts.push_attached("AFLOWIN_FLAG::PSTRESS","0.0");};
                                        if(!PARAMS->vparams.flag("AFLOWIN_FLAG::LIST")) {
                                            if(!MULTITHREAD) AVASP_MakeSingleAFLOWIN(xaus,!PARAMS->vparams.flag("AFLOWIN_FLAG::STDOUT"));   // SINGLE_THREAD
                                            if(MULTITHREAD) {dxvasp.push_back(xaus);if(dxvasp.size()==AVASP_AFLOWIN_MAX_JOBS) AVASP_MakePrototype_AFLOWIN_LOOP(dxvasp,!PARAMS->vparams.flag("AFLOWIN_FLAG::STDOUT"));}
                                        } else {
                                            if(xaus.aopts.flag("AFLOWIN_FLAG::PSTRESS")==FALSE) {
                                                cout << "aflow " << PARAMS->vparams.getattachedscheme("AFLOWIN_FLAG::LIST_VCMD") << "--aflow_proto=" << xaus.AVASP_label; 
                                                for(uint i=0;i<nspecies;i++) cout << ":" << KBIN::VASP_PseudoPotential_CleanName(xaus.str.species.at(i));
                                                cout << endl;
                                            } else {
                                                cout << "aflow " << PARAMS->vparams.getattachedscheme("AFLOWIN_FLAG::LIST_VCMD") << "--pressure=" << xaus.aopts.getattachedscheme("AFLOWIN_FLAG::PSTRESS") << " --aflow_proto=" << xaus.AVASP_label; 
                                                for(uint i=0;i<nspecies;i++) cout << ":" << KBIN::VASP_PseudoPotential_CleanName(xaus.str.species.at(i)); 
                                                cout << endl;
                                            }
                                        }
                                    }
                                }
                            } else { // no alphabetic
                                for(uint i=0;i<nspecies-1;i++) {
                                    if(specieX.at(i).at(ii[i])>specieX.at(i+1).at(ii[i+1])) {
                                        cerr << soliloquy << " system ";
                                        for(uint j=0;j<nspecies;j++) {
                                            cerr << specieX.at(j).at(ii[j]);
                                        }
                                        cerr << " neglected (not alphabetic)" << endl;
                                        // exit(0);
                                    }
                                }
                            }
                        }
                    }
                }
            }
            if(MULTITHREAD) AVASP_MakePrototype_AFLOWIN_LOOP(dxvasp,!PARAMS->vparams.flag("AFLOWIN_FLAG::STDOUT")); // FLUSH OUT
        }
        // ***************************************************************************

        return TRUE;
}

// ***************************************************************************
bool AVASP_DefaultValuesICSD_AFLOWIN(_xvasp &xvasp) {        
    xvasp.clear();
    xvasp.AVASP_dirbase="./ICSD";
    xvasp.AVASP_prototype_mode=LIBRARY_MODE_ICSD;
    xvasp.AVASP_prototype_from_library_=FALSE;
    xvasp.AVASP_directory_from_library_=FALSE;
    xvasp.aopts.flag("AFLOWIN_FLAG::NBANDS",TRUE);
    xvasp.aopts.flag("FLAG::AVASP_CHGCAR",FALSE);
    xvasp.aopts.flag("FLAG::AVASP_SPIN",TRUE);
    xvasp.aopts.flag("FLAG::AVASP_SPIN_REMOVE_RELAX_1",DEFAULT_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_1);
    xvasp.aopts.flag("FLAG::AVASP_SPIN_REMOVE_RELAX_2",DEFAULT_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_2);
    xvasp.aopts.flag("FLAG::AVASP_AUTO_MAGMOM",FALSE);
    xvasp.AVASP_value_KPPRA=DEFAULT_KPPRA;
    if(xvasp.aopts.flag("AFLOWIN_FLAG::KPPRA")) xvasp.AVASP_value_KPPRA=aurostd::string2utype<int>(xvasp.aopts.getattachedscheme("AFLOWIN_FLAG::KPPRA"));    // KPPRA
    xvasp.AVASP_KSCHEME=DEFAULT_KSCHEME;
    xvasp.AVASP_value_KPPRA_STATIC=DEFAULT_KPPRA_STATIC;
    if(xvasp.aopts.flag("AFLOWIN_FLAG::KPPRA_STATIC")) xvasp.AVASP_value_KPPRA_STATIC=aurostd::string2utype<int>(xvasp.aopts.getattachedscheme("AFLOWIN_FLAG::KPPRA_STATIC")); // KPPRA_STATIC
    xvasp.AVASP_STATIC_KSCHEME=DEFAULT_STATIC_KSCHEME;
    xvasp.AVASP_value_NSW=-1;
    xvasp.aopts.flag("FLAG::PRECISION_SET",TRUE);
    xvasp.AVASP_flag_PRECISION_scheme=DEFAULT_VASP_FORCE_OPTION_PREC_SCHEME;
    xvasp.aopts.flag("FLAG::ALGO_SET",TRUE);
    xvasp.AVASP_flag_ALGO_scheme="N";
    // [OBSOLETE] xvasp.aopts.flag("FLAG::METAGGA_SET",FALSE);
    // [OBSOLETE] xvasp.AVASP_flag_METAGGA_scheme=DEFAULT_VASP_FORCE_OPTION_METAGGA_SCHEME;
    // [OBSOLETE] xvasp.aopts.flag("FLAG::IVDW_SET",FALSE);
    // [OBSOLETE] xvasp.AVASP_flag_IVDW_scheme=DEFAULT_VASP_FORCE_OPTION_IVDW_SCHEME;
    //  xvasp.AVASP_flag_TYPE.xscheme.at(0)=='I'=TRUE; wrong
    xvasp.aopts.flag("FLAG::VOLUME_PRESERVED",TRUE);
    xvasp.aopts.flag("FLAG::EXTRA_INCAR",FALSE);
    xvasp.AVASP_EXTRA_INCAR.clear();
    xvasp.AVASP_volume_in=-1.0;
    xvasp.AVASP_flag_RUN_RELAX=TRUE;
    xvasp.AVASP_path_BANDS="fcc";
    xvasp.AVASP_value_BANDS_GRID=DEFAULT_BANDS_GRID;
    if(xvasp.aopts.flag("AFLOWIN_FLAG::BANDS_GRID")) xvasp.AVASP_value_BANDS_GRID=aurostd::string2utype<uint>(xvasp.aopts.getattachedscheme("AFLOWIN_FLAG::BANDS_GRID")); // BANDS_GRID
    return TRUE;
}


bool AVASP_ADD_LDAU(_xvasp &xvasp) {
    xvasp.aopts.flag("FLAG::AVASP_LDAU1",FALSE); // DEFAULT search
    xvasp.aopts.flag("FLAG::AVASP_LDAU2",FALSE); // DEFAULT search
    xvasp.AVASP_LDAU_PARAMETERS_STRING="";
    xvasp.AVASP_LDAU_PARAMETERS_UJSUM=0.0;
    vector<string> vLDAUspecies;vector<uint> vLDAUtype;vector<int> vLDAUL;vector<double> vLDAUU,vLDAUJ; // objects for the search
    bool LDAU=FALSE;vLDAUspecies.clear();vLDAUtype.clear();vLDAUL.clear();vLDAUU.clear();vLDAUJ.clear();
    for(uint i=0;i<xvasp.str.species.size();i++) {
        AVASP_Get_LDAU_Parameters(aurostd::CleanStringASCII(KBIN::VASP_PseudoPotential_CleanName(xvasp.str.species.at(i))),LDAU,vLDAUspecies,vLDAUtype,vLDAUL,vLDAUU,vLDAUJ); // parameters for LDAU2
    }
    // check if LDAU really needed 12/08/27
    for(uint i=0;i<xvasp.str.species.size();i++) xvasp.AVASP_LDAU_PARAMETERS_UJSUM+=aurostd::abs(vLDAUU.at(i))+aurostd::abs(vLDAUJ.at(i));
    if(aurostd::abs(xvasp.AVASP_LDAU_PARAMETERS_UJSUM)<0.01) LDAU=FALSE;
    if(aurostd::abs(xvasp.AVASP_LDAU_PARAMETERS_UJSUM)>0.01) LDAU=TRUE;

    if(LDAU==TRUE) { // got LDAU
        for(uint i=0;i<xvasp.str.species.size();i++) {
            if(vLDAUtype.at(i)==0) {;} // do nothing
            if(vLDAUtype.at(i)==1) xvasp.aopts.flag("FLAG::AVASP_LDAU1",TRUE);
            if(vLDAUtype.at(i)==2) xvasp.aopts.flag("FLAG::AVASP_LDAU2",TRUE);
        }
        LDAU=(xvasp.aopts.flag("FLAG::AVASP_LDAU1") || xvasp.aopts.flag("FLAG::AVASP_LDAU2"));
        if(xvasp.aopts.flag("FLAG::AVASP_LDAU1")==TRUE && xvasp.aopts.flag("FLAG::AVASP_LDAU2")==TRUE) {
            cerr << "AVASP_MakePrototypeICSD_AFLOWIN: you can not be here: xvasp.aopts.flag(\"FLAG::AVASP_LDAU1\")==TRUE && xvasp.aopts.flag(\"FLAG::AVASP_LDAU2\")==TRUE: need to get different LDAU parameterization" << endl;
            exit(0);
        }
        stringstream aus;
        aus.clear();aus.str(std::string());
        for(uint i=0;i<xvasp.str.species.size();i++) aus << vLDAUspecies.at(i) << (i!=xvasp.str.species.size()-1 ? "," : ";");
        for(uint i=0;i<xvasp.str.species.size();i++) aus << vLDAUL.at(i) << (i!=xvasp.str.species.size()-1 ? "," : ";");
        for(uint i=0;i<xvasp.str.species.size();i++) aus << vLDAUU.at(i) << (i!=xvasp.str.species.size()-1 ? "," : ";");
        for(uint i=0;i<xvasp.str.species.size();i++) aus << vLDAUJ.at(i) << (i!=xvasp.str.species.size()-1 ? "," : " ");

        // KBIN::VASP_PseudoPotential_CleanName(xvasp.str.species.at(i)) << (i!=xvasp.str.species.size()-1 ? "," : ";") ;

        xvasp.AVASP_LDAU_PARAMETERS_STRING=aus.str();
    }
    xvasp.aopts.flag("FLAG::AVASP_FORCE_LDAU",TRUE);
    //  cerr << xvasp.AVASP_LDAU_PARAMETERS_UJSUM << endl;
    return TRUE;
}

bool AVASP_REMOVE_LDAU(_xvasp &xvasp) {
    xvasp.aopts.flag("FLAG::AVASP_LDAU1",FALSE); // DEFAULT search
    xvasp.aopts.flag("FLAG::AVASP_LDAU2",FALSE); // DEFAULT search
    xvasp.AVASP_LDAU_PARAMETERS_STRING="";
    xvasp.AVASP_LDAU_PARAMETERS_UJSUM=0.0;
    xvasp.str.species_pp_vLDAU.clear();
    xvasp.aopts.flag("FLAG::AVASP_FORCE_NOLDAU",TRUE);
    //  cerr << "AVASP_REMOVE_LDAU" << endl;
    return TRUE;
}

//bool AVASP_MakePrototypeICSD_AFLOWIN(vector<string> params_UCELL_SPECIES,bool flag_AFLOW_IN_ONLY_IF_MISSING) {
bool AVASP_MakePrototypeICSD_AFLOWIN(_AVASP_PROTO *PARAMS,bool flag_AFLOW_IN_ONLY_IF_MISSING) {

    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string soliloquy="AVASP_MakePrototypeICSD_AFLOWIN():";
    if(LDEBUG) cerr << soliloquy << " 0a" << endl;
    _xvasp xvasp;
    AVASP_DefaultValuesICSD_AFLOWIN(xvasp);
    string parameters="";

    if(PARAMS->vparams.flag("AFLOWIN_FLAG::AUTOLDAU")) xvasp.aopts.flag("FLAG::AVASP_FORCE_LDAU",TRUE);
    if(PARAMS->vparams.flag("AFLOWIN_FLAG::AUTONOLDAU")) xvasp.aopts.flag("FLAG::AVASP_FORCE_NOLDAU",TRUE);
    if(PARAMS->vparams.flag("AFLOWIN_FLAG::VASP")) xvasp.aopts.flag("AFLOWIN_FLAG::VASP",TRUE);
    if(PARAMS->vparams.flag("AFLOWIN_FLAG::ABINIT")) xvasp.aopts.flag("AFLOWIN_FLAG::ABINIT",TRUE);
    if(PARAMS->vparams.flag("AFLOWIN_FLAG::AIMS")) xvasp.aopts.flag("AFLOWIN_FLAG::AIMS",TRUE);
    if(PARAMS->vparams.flag("AFLOWIN_FLAG::QE")) xvasp.aopts.flag("AFLOWIN_FLAG::QE",TRUE);
    if(PARAMS->vparams.flag("AFLOWIN_FLAG::MODULE")) {xvasp.aopts.flag("AFLOWIN_FLAG::MODULE",TRUE);xvasp.aopts.push_attached("AFLOWIN_FLAG::MODULE",PARAMS->vparams.getattachedscheme("AFLOWIN_FLAG::MODULE"));}  // CO 180214
    if(PARAMS->vparams.flag("AFLOWIN_FLAG::APL_SUPERCELL")) {xvasp.aopts.flag("AFLOWIN_FLAG::APL_SUPERCELL",TRUE);xvasp.aopts.push_attached("AFLOWIN_FLAG::APL_SUPERCELL",PARAMS->vparams.getattachedscheme("AFLOWIN_FLAG::APL_SUPERCELL"));}  // CO 180214
    if(PARAMS->vparams.flag("AFLOWIN_FLAG::POTIM")) {xvasp.aopts.flag("AFLOWIN_FLAG::POTIM",TRUE);xvasp.aopts.push_attached("AFLOWIN_FLAG::POTIM",PARAMS->vparams.getattachedscheme("AFLOWIN_FLAG::POTIM"));}
    if(PARAMS->vparams.flag("AFLOWIN_FLAG::PRECISION")) {xvasp.aopts.flag("AFLOWIN_FLAG::PRECISION",TRUE);xvasp.aopts.push_attached("AFLOWIN_FLAG::PRECISION",PARAMS->vparams.getattachedscheme("AFLOWIN_FLAG::PRECISION"));}
    if(PARAMS->vparams.flag("AFLOWIN_FLAG::ALGORITHM")) {xvasp.aopts.flag("AFLOWIN_FLAG::ALGORITHM",TRUE);xvasp.aopts.push_attached("AFLOWIN_FLAG::ALGORITHM",PARAMS->vparams.getattachedscheme("AFLOWIN_FLAG::ALGORITHM"));}
    if(PARAMS->vparams.flag("AFLOWIN_FLAG::METAGGA")) {xvasp.aopts.flag("AFLOWIN_FLAG::METAGGA",TRUE);xvasp.aopts.push_attached("AFLOWIN_FLAG::METAGGA",PARAMS->vparams.getattachedscheme("AFLOWIN_FLAG::METAGGA"));}  
    if(PARAMS->vparams.flag("AFLOWIN_FLAG::IVDW")) {xvasp.aopts.flag("AFLOWIN_FLAG::IVDW",TRUE);xvasp.aopts.push_attached("AFLOWIN_FLAG::IVDW",PARAMS->vparams.getattachedscheme("AFLOWIN_FLAG::IVDW"));}  
    if(PARAMS->vparams.flag("AFLOWIN_FLAG::RELAX_TYPE")) {xvasp.aopts.flag("AFLOWIN_FLAG::RELAX_TYPE",TRUE);xvasp.aopts.push_attached("AFLOWIN_FLAG::RELAX_TYPE",PARAMS->vparams.getattachedscheme("AFLOWIN_FLAG::RELAX_TYPE"));}  // CO 180214
    if(PARAMS->vparams.flag("AFLOWIN_FLAG::RELAX_MODE")) {xvasp.aopts.flag("AFLOWIN_FLAG::RELAX_MODE",TRUE);xvasp.aopts.push_attached("AFLOWIN_FLAG::RELAX_MODE",PARAMS->vparams.getattachedscheme("AFLOWIN_FLAG::RELAX_MODE"));}
    if(PARAMS->vparams.flag("AFLOWIN_FLAG::TYPE")) {xvasp.aopts.flag("AFLOWIN_FLAG::TYPE",TRUE);xvasp.aopts.push_attached("AFLOWIN_FLAG::TYPE",PARAMS->vparams.getattachedscheme("AFLOWIN_FLAG::TYPE"));}
    if(PARAMS->vparams.flag("AFLOWIN_FLAG::CONVERT_UNIT_CELL")) {xvasp.aopts.flag("AFLOWIN_FLAG::CONVERT_UNIT_CELL",TRUE);xvasp.aopts.push_attached("AFLOWIN_FLAG::CONVERT_UNIT_CELL",PARAMS->vparams.getattachedscheme("AFLOWIN_FLAG::CONVERT_UNIT_CELL"));}
    if(PARAMS->vparams.flag("AFLOWIN_FLAG::VOLUME_PLUS_EQUAL")) {xvasp.aopts.flag("AFLOWIN_FLAG::VOLUME_PLUS_EQUAL",TRUE);xvasp.aopts.push_attached("AFLOWIN_FLAG::VOLUME_PLUS_EQUAL",PARAMS->vparams.getattachedscheme("AFLOWIN_FLAG::VOLUME_PLUS_EQUAL"));}
    if(PARAMS->vparams.flag("AFLOWIN_FLAG::VOLUME_MULTIPLY_EQUAL")) {xvasp.aopts.flag("AFLOWIN_FLAG::VOLUME_MULTIPLY_EQUAL",TRUE);xvasp.aopts.push_attached("AFLOWIN_FLAG::VOLUME_MULTIPLY_EQUAL",PARAMS->vparams.getattachedscheme("AFLOWIN_FLAG::VOLUME_MULTIPLY_EQUAL"));}
    if(PARAMS->vparams.flag("AFLOWIN_FLAG::NO_VOLUME_ADJUSTMENT")) xvasp.aopts.flag("AFLOWIN_FLAG::NO_VOLUME_ADJUSTMENT",TRUE); // CO 180214
    if(PARAMS->vparams.flag("AFLOWIN_FLAG::EDIFFG")) {xvasp.aopts.flag("AFLOWIN_FLAG::EDIFFG",TRUE);xvasp.aopts.push_attached("AFLOWIN_FLAG::EDIFFG",PARAMS->vparams.getattachedscheme("AFLOWIN_FLAG::EDIFFG"));}
    if(PARAMS->vparams.flag("AFLOWIN_FLAG::KPPRA")) {xvasp.aopts.flag("AFLOWIN_FLAG::KPPRA",TRUE);xvasp.aopts.push_attached("AFLOWIN_FLAG::KPPRA",PARAMS->vparams.getattachedscheme("AFLOWIN_FLAG::KPPRA"));}
    if(PARAMS->vparams.flag("AFLOWIN_FLAG::KPPRA_STATIC")) {xvasp.aopts.flag("AFLOWIN_FLAG::KPPRA_STATIC",TRUE);xvasp.aopts.push_attached("AFLOWIN_FLAG::KPPRA_STATIC",PARAMS->vparams.getattachedscheme("AFLOWIN_FLAG::KPPRA_STATIC"));}
    if(PARAMS->vparams.flag("AFLOWIN_FLAG::BANDS_GRID")) {xvasp.aopts.flag("AFLOWIN_FLAG::BANDS_GRID",TRUE);xvasp.aopts.push_attached("AFLOWIN_FLAG::BANDS_GRID",PARAMS->vparams.getattachedscheme("AFLOWIN_FLAG::BANDS_GRID"));}
    if(PARAMS->vparams.flag("AFLOWIN_FLAG::ENMAX_MULTIPLY")) {xvasp.aopts.flag("AFLOWIN_FLAG::ENMAX_MULTIPLY",TRUE);xvasp.aopts.push_attached("AFLOWIN_FLAG::ENMAX_MULTIPLY",PARAMS->vparams.getattachedscheme("AFLOWIN_FLAG::ENMAX_MULTIPLY"));}

    // recompile and ahave fun
    // cerr << soliloquy << " 0b" << endl;
    if(flag_AFLOW_IN_ONLY_IF_MISSING==FALSE) xvasp.AVASP_aflowin_only_if_missing=FALSE;
    if(flag_AFLOW_IN_ONLY_IF_MISSING==TRUE) xvasp.AVASP_aflowin_only_if_missing=TRUE;
    string label=PARAMS->vparams.getattachedscheme("AFLOWIN_STRING::LABEL");
    // get the prototype
    xvasp.AVASP_label=label;
    if(LDEBUG) cerr << soliloquy << " 0c" << endl;

    xvasp.AVASP_dirbase=_ICSD_DIRBASE_;
    xvasp.AVASP_libbase="./";
    if(XHOST_LIBRARY_ICSD+1<vAFLOW_PROJECTS_DIRECTORIES.size()) xvasp.AVASP_libbase=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_ICSD);
    xvasp.AVASP_prototype_mode=LIBRARY_MODE_ICSD;

    bool DEBUG_SKIP=TRUE;

    if(LDEBUG) cerr << soliloquy << " 0d" << endl;

    if(xvasp.AVASP_aflowin_only_if_missing) { 
        // this test is redundant, but since the Prototype generation is slow, it is nice to avoid unless you really need it
        // TEST IF CALCULATED
        if(aurostd::substring2bool(init::InitGlobalObject("Library_CALCULATED_ICSD_LIB"),xvasp.AVASP_label+" ")) {if(DEBUG_SKIP) cerr << "SKIP (calculated): " << xvasp.AVASP_label << endl; return TRUE;}
        for(uint ilattice=1;ilattice<=14;ilattice++) {
            if(aurostd::substring2bool(init::InitGlobalObject("Library_CALCULATED_ICSD_LIB"),lattices[ilattice]+"/"+xvasp.AVASP_label+" ")) {
                if(DEBUG_SKIP) cerr << "SKIP (calculated): " << lattices[ilattice] << "/" << xvasp.AVASP_label << endl;
                return TRUE;
            }
        }
        // TEST IF PREPARED
        string label;
        label="./ICSD/"+xvasp.AVASP_label+"/"+_AFLOWIN_;if(aurostd::FileExist(label)) {if(DEBUG_SKIP) cerr << "SKIP (prepared): " << xvasp.AVASP_label << endl; return TRUE;}
        for(uint ilattice=1;ilattice<=14;ilattice++) {
            label="./ICSD/"+lattices[ilattice]+"/"+xvasp.AVASP_label+"/"+_AFLOWIN_;
            if(aurostd::FileExist(label)) {
                if(DEBUG_SKIP) cerr << "SKIP (prepared): " << lattices[ilattice] << "/" << xvasp.AVASP_label << endl;
                return TRUE;
            }
        }
    }

    if(LDEBUG) cerr << soliloquy << " 0d-1" << endl;
    xvasp.str=aflowlib::PrototypeLibraries(cerr,label,parameters,LIBRARY_MODE_ICSD);
    if(LDEBUG) cerr << soliloquy << " 0d-2" << endl;
    if(LDEBUG) cerr << soliloquy << " DEBUG: xvasp.str.bravais_lattice_type=" << xvasp.str.bravais_lattice_type << endl;
    if(LDEBUG) cerr << soliloquy << " DEBUG: xvasp.str.bravais_lattice_variation_type=" << xvasp.str.bravais_lattice_variation_type << endl;
    if(LDEBUG) cerr << soliloquy << " 0d-3" << endl;

    if(xvasp.str.bravais_lattice_type.empty()){xvasp.str.GetLatticeType();} //CO 180622 - if we download the structure from online, bravais_lattice_type and bravais_lattice_variation_type are empty

    if(LATTICE::lattice_is_working(xvasp.str.bravais_lattice_type)==FALSE) exit(0); // DO NOT PRODUCE THE ONES NOT READY TO BRILLOUIN
    if(LDEBUG) cerr << soliloquy << " 0d-4" << endl;
    if(xvasp.str.error_flag==TRUE) {         // flag TRUE is error
        cerr << "ERROR " << soliloquy << " " << endl << xvasp.str.error_string << endl;
        exit(0);
    }
    if(LDEBUG) cerr << soliloquy << " 0d-5" << endl;

    // get the species
    bool potpaw_PBE=TRUE,potpaw_PBE_KIN=FALSE,potpaw_LDA_KIN=FALSE,potpaw_GGA=FALSE,potpaw_LDA=FALSE,pot_GGA=FALSE,pot_PBE=FALSE,pot_LDA=FALSE;
    if(LDEBUG) cerr << soliloquy << " DEBUG: label=" << label << endl;
    //  if(aurostd::substring2bool(label,"Cs")) {potpaw_PBE=FALSE;potpaw_PBE_KIN=FALSE;potpaw_LDA_KIN=FALSE;potpaw_GGA=FALSE;pot_GGA=TRUE;pot_PBE=FALSE;pot_LDA=FALSE;}
    if(aurostd::substring2bool(label,"Cs")) {potpaw_PBE=TRUE;potpaw_PBE_KIN=FALSE;potpaw_LDA_KIN=FALSE;potpaw_GGA=FALSE;pot_GGA=FALSE;pot_PBE=FALSE;pot_LDA=FALSE;}
    //  if(aurostd::substring2bool(label,"Sm")) {potpaw_PBE=FALSE;potpaw_PBE_KIN=FALSE;potpaw_LDA_KIN=FALSE;potpaw_GGA=TRUE;pot_GGA=FALSE;pot_PBE=FALSE;pot_LDA=FALSE;}
    if(potpaw_PBE_KIN==TRUE) xvasp.AVASP_potential="potpaw_PBE_KIN";
    if(potpaw_LDA_KIN==TRUE) xvasp.AVASP_potential="potpaw_LDA_KIN";
    if(potpaw_PBE==TRUE)     xvasp.AVASP_potential="potpaw_PBE";
    if(potpaw_GGA==TRUE)     xvasp.AVASP_potential="potpaw_GGA";
    if(potpaw_LDA==TRUE)     xvasp.AVASP_potential="potpaw_LDA";
    if(pot_GGA==TRUE)        xvasp.AVASP_potential="pot_GGA";
    if(pot_PBE==TRUE)        xvasp.AVASP_potential="pot_PBE";
    if(pot_LDA==TRUE)        xvasp.AVASP_potential="pot_LDA";
    if(LDEBUG) cerr << soliloquy << " DEBUG: potential=" << xvasp.AVASP_potential << endl;

    // LDAU
    if(!PARAMS->vparams.flag("AFLOWIN_FLAG::AUTONOLDAU")) {
        bool CHECK_LDAU=TRUE;
        if(CHECK_LDAU) AVASP_ADD_LDAU(xvasp);
    }

    // POTENTIALS
    for(uint i=0;i<xvasp.str.species.size();i++) {
        if(potpaw_PBE_KIN==TRUE) xvasp.str.species.at(i)=AVASP_Get_PseudoPotential_PAW_PBE_KIN(aurostd::CleanStringASCII(xvasp.str.species.at(i)));
        if(potpaw_LDA_KIN==TRUE) xvasp.str.species.at(i)=AVASP_Get_PseudoPotential_PAW_LDA_KIN(aurostd::CleanStringASCII(xvasp.str.species.at(i)));
        if(potpaw_PBE==TRUE) xvasp.str.species.at(i)=AVASP_Get_PseudoPotential_PAW_PBE(aurostd::CleanStringASCII(xvasp.str.species.at(i)));
        if(potpaw_GGA==TRUE) xvasp.str.species.at(i)=AVASP_Get_PseudoPotential_PAW_GGA(aurostd::CleanStringASCII(xvasp.str.species.at(i)));
        if(potpaw_LDA==TRUE) xvasp.str.species.at(i)=AVASP_Get_PseudoPotential_PAW_LDA(aurostd::CleanStringASCII(xvasp.str.species.at(i)));
        if(pot_GGA==TRUE) xvasp.str.species.at(i)=AVASP_Get_PseudoPotential_GGA(aurostd::CleanStringASCII(xvasp.str.species.at(i)));
        if(pot_PBE==TRUE) xvasp.str.species.at(i)=AVASP_Get_PseudoPotential_PBE(aurostd::CleanStringASCII(xvasp.str.species.at(i)));
        if(pot_LDA==TRUE) xvasp.str.species.at(i)=AVASP_Get_PseudoPotential_LDA(aurostd::CleanStringASCII(xvasp.str.species.at(i)));
    }

    // Parameters
    if(LDEBUG) cerr << soliloquy << " 1a" << endl;
    if(!xvasp.aopts.flag("FLAG::AVASP_KPPRA")) xvasp.AVASP_value_KPPRA=DEFAULT_KPPRA_ICSD;                                     
    if(xvasp.aopts.flag("AFLOWIN_FLAG::KPPRA")) xvasp.AVASP_value_KPPRA=aurostd::string2utype<int>(xvasp.aopts.getattachedscheme("AFLOWIN_FLAG::KPPRA"));    // KPPRA
    xvasp.AVASP_KSCHEME=DEFAULT_KSCHEME;                                                  
    xvasp.AVASP_value_KPPRA_STATIC=DEFAULT_KPPRA_STATIC;                            
    if(xvasp.aopts.flag("AFLOWIN_FLAG::KPPRA_STATIC")) xvasp.AVASP_value_KPPRA_STATIC=aurostd::string2utype<int>(xvasp.aopts.getattachedscheme("AFLOWIN_FLAG::KPPRA_STATIC")); // KPPRA_STATIC
    xvasp.AVASP_STATIC_KSCHEME=DEFAULT_STATIC_KSCHEME;
    // kpoints for the brillouin zone
    xvasp.AVASP_path_BANDS=xvasp.str.bravais_lattice_variation_type;  // ICSD BASTARDS
    xvasp.AVASP_value_BANDS_GRID=DEFAULT_BANDS_GRID;                               
    if(xvasp.aopts.flag("AFLOWIN_FLAG::BANDS_GRID")) xvasp.AVASP_value_BANDS_GRID=aurostd::string2utype<uint>(xvasp.aopts.getattachedscheme("AFLOWIN_FLAG::BANDS_GRID")); // BANDS_GRID
    if(xvasp.AVASP_path_BANDS=="HEX") {
        xvasp.AVASP_KSCHEME="G";           // HEXAGONAL SYSTEMS GET GAMMA KPOINTS GRID
        xvasp.AVASP_STATIC_KSCHEME="G";    // HEXAGONAL SYSTEMS GET GAMMA KPOINTS GRID
    }
    // MINKOWSKI
    xvasp.aopts.flag("AVASP_flag_CONVERT_UNIT_CELL_MINKOWSKI",FALSE);

    // rest
    xvasp.aopts.flag("FLAG::AVASP_SKIP_NOMIX",FALSE);                                        
    // xvasp.AVASP_flag_TYPE.xscheme.at(0)=='M'=FALSE;
    // xvasp.AVASP_flag_TYPE.xscheme.at(0)=='I'=TRUE;
    xvasp.AVASP_flag_MPI=FALSE; // will fix it by itself
    xvasp.AVASP_flag_RUN_RELAX=FALSE;
    xvasp.AVASP_flag_RUN_RELAX_STATIC=FALSE;
    xvasp.AVASP_flag_RUN_RELAX_STATIC_BANDS=TRUE;
    xvasp.AVASP_flag_RUN_STATIC_BANDS=FALSE;

    xvasp.aopts.flag("FLAG::PRECISION_SET",TRUE);
    xvasp.AVASP_flag_PRECISION_scheme=DEFAULT_VASP_FORCE_OPTION_PREC_SCHEME;
    xvasp.aopts.flag("FLAG::ALGO_SET",TRUE);
    xvasp.AVASP_flag_ALGO_scheme="F";
    // [OBSOLETE] xvasp.aopts.flag("FLAG::METAGGA_SET",FALSE);
    // [OBSOLETE] xvasp.AVASP_flag_METAGGA_scheme=DEFAULT_VASP_FORCE_OPTION_METAGGA_SCHEME;
    // [OBSOLETE] xvasp.aopts.flag("FLAG::IVDW_SET",FALSE);
    // [OBSOLETE] xvasp.AVASP_flag_IVDW_scheme=DEFAULT_VASP_FORCE_OPTION_IVDW_SCHEME;

    // extra INCAR
    xvasp.aopts.flag("FLAG::EXTRA_INCAR",TRUE);                               
    xvasp.AVASP_EXTRA_INCAR.clear();                                          
    xvasp.AVASP_EXTRA_INCAR << "NELM = 120" << endl;                           // WAHYU DEFAULT
    xvasp.AVASP_EXTRA_INCAR << "NELMIN=2"  << endl;                            // WAHYU DEFAULT
    xvasp.AVASP_EXTRA_INCAR << "LPLANE=.TRUE." << endl;                        // WAHYU DEFAULT
    if(xvasp.str.atoms.size()<=10)                     // cutoff for LREAL     // WAHYU DEFAULT
        xvasp.AVASP_EXTRA_INCAR << "LREAL=.FALSE." << endl;                      // WAHYU DEFAULT
    else
        xvasp.AVASP_EXTRA_INCAR << "LREAL=Auto" << endl;                         // WAHYU DEFAULT
    xvasp.AVASP_EXTRA_INCAR << "LSCALU=.FALSE." << endl;                       // WAHYU DEFAULT

    // make the _AFLOWIN_

    if(LDEBUG) cerr << soliloquy << " 2a" << endl;
    if(LDEBUG) cerr << xvasp.str << endl;
    if(LDEBUG) cerr << soliloquy << " Running AVASP_MakeSingleAFLOWIN" << endl;
    AVASP_MakeSingleAFLOWIN(xvasp,!PARAMS->vparams.flag("AFLOWIN_FLAG::STDOUT"));
    if(LDEBUG) cerr << soliloquy << " 3a" << endl;

    return TRUE;
}

// ***************************************************************************
string AVASP_Shortcuts_for_Binaries(string &label) {

    // SEARCH FOR LABELS
    if(label=="LIBRARY_PROTOTYPES" || label=="LIBRARY_PROTOS" || label=="PROTOS" || label=="LIB2"|| label=="lib2" || label=="ALL"|| label=="all") label=LABEL_RAW_LIBRARY_PROTOTYPES;

    // SIGMA PHASE STUFF
    if(label=="600.XXXXX" || label=="600.XXXX" || label=="600.XXX" || label=="600.XX" || label=="600.X")
        label="600.AAAAA,600.AAAAB,600.AAABA,600.AAABB,600.AABAA,600.AABAB,600.AABBA,600.AABBB,600.ABAAA,600.ABAAB,600.ABABA,600.ABABB,600.ABBAA,600.ABBAB,600.ABBBA,600.ABBBB,600.BAAAA,600.BAAAB,600.BAABA,600.BAABB,600.BABAA,600.BABAB,600.BABBA,600.BABBB,600.BBAAA,600.BBAAB,600.BBABA,600.BBABB,600.BBBAA,600.BBBAB,600.BBBBA,600.BBBBB"; // all the  sigma.ABCDE CrFe stuff

    // GUS STUFF some FCC/BCC/HCP shortcuts
    stringstream ss;
    // FCC
    if(aurostd::substring2bool(label,"fccnAt1-16")) aurostd::StringSubst(label,"fccnAt1-16","fccnAt1,fccnAt2,fccnAt3,fccnAt4,fccnAt5,fccnAt6,fccnAt7,fccnAt8,fccnAt9,fccnAt10,fccnAt11,fccnAt12,fccnAt13,fccnAt14,fccnAt15,fccnAt16");
    if(aurostd::substring2bool(label,"fccnAt1-15")) aurostd::StringSubst(label,"fccnAt1-15","fccnAt1,fccnAt2,fccnAt3,fccnAt4,fccnAt5,fccnAt6,fccnAt7,fccnAt8,fccnAt9,fccnAt10,fccnAt11,fccnAt12,fccnAt13,fccnAt14,fccnAt15");
    if(aurostd::substring2bool(label,"fccnAt1-14")) aurostd::StringSubst(label,"fccnAt1-14","fccnAt1,fccnAt2,fccnAt3,fccnAt4,fccnAt5,fccnAt6,fccnAt7,fccnAt8,fccnAt9,fccnAt10,fccnAt11,fccnAt12,fccnAt13,fccnAt14");
    if(aurostd::substring2bool(label,"fccnAt1-13")) aurostd::StringSubst(label,"fccnAt1-13","fccnAt1,fccnAt2,fccnAt3,fccnAt4,fccnAt5,fccnAt6,fccnAt7,fccnAt8,fccnAt9,fccnAt10,fccnAt11,fccnAt12,fccnAt13");
    if(aurostd::substring2bool(label,"fccnAt1-12")) aurostd::StringSubst(label,"fccnAt1-12","fccnAt1,fccnAt2,fccnAt3,fccnAt4,fccnAt5,fccnAt6,fccnAt7,fccnAt8,fccnAt9,fccnAt10,fccnAt11,fccnAt12");
    if(aurostd::substring2bool(label,"fccnAt1-11")) aurostd::StringSubst(label,"fccnAt1-11","fccnAt1,fccnAt2,fccnAt3,fccnAt4,fccnAt5,fccnAt6,fccnAt7,fccnAt8,fccnAt9,fccnAt10,fccnAt11");
    if(aurostd::substring2bool(label,"fccnAt1-10")) aurostd::StringSubst(label,"fccnAt1-10","fccnAt1,fccnAt2,fccnAt3,fccnAt4,fccnAt5,fccnAt6,fccnAt7,fccnAt8,fccnAt9,fccnAt10");
    if(aurostd::substring2bool(label,"fccnAt1-9")) aurostd::StringSubst(label,"fccnAt1-9","fccnAt1,fccnAt2,fccnAt3,fccnAt4,fccnAt5,fccnAt6,fccnAt7,fccnAt8,fccnAt9");
    if(aurostd::substring2bool(label,"fccnAt1-8")) aurostd::StringSubst(label,"fccnAt1-8","fccnAt1,fccnAt2,fccnAt3,fccnAt4,fccnAt5,fccnAt6,fccnAt7,fccnAt8");
    if(aurostd::substring2bool(label,"fccnAt1-7")) aurostd::StringSubst(label,"fccnAt1-7","fccnAt1,fccnAt2,fccnAt3,fccnAt4,fccnAt5,fccnAt6,fccnAt7");
    if(aurostd::substring2bool(label,"fccnAt1-6")) aurostd::StringSubst(label,"fccnAt1-6","fccnAt1,fccnAt2,fccnAt3,fccnAt4,fccnAt5,fccnAt6");
    if(aurostd::substring2bool(label,"fccnAt1-5")) aurostd::StringSubst(label,"fccnAt1-5","fccnAt1,fccnAt2,fccnAt3,fccnAt4,fccnAt5");
    if(aurostd::substring2bool(label,"fccnAt1-4")) aurostd::StringSubst(label,"fccnAt1-4","fccnAt1,fccnAt2,fccnAt3,fccnAt4");
    if(aurostd::substring2bool(label,"fccnAt1-3")) aurostd::StringSubst(label,"fccnAt1-3","fccnAt1,fccnAt2,fccnAt3");
    if(aurostd::substring2bool(label,"fccnAt1-2")) aurostd::StringSubst(label,"fccnAt1-2","fccnAt1,fccnAt2");
    if(aurostd::substring2bool(label,"fccnAt1-1")) aurostd::StringSubst(label,"fccnAt1-1","fccnAt1");
    string fcc_str="f";
    uint fcc_start1=1; uint fcc_start2=3; uint fcc_start3=5; uint fcc_start4=11;
    uint fcc_start5=30; uint fcc_start6=58; uint fcc_start7=138; uint fcc_start8=242;
    uint fcc_start9=632; uint fcc_start10=1136; uint fcc_start11=2347; uint fcc_start12=3711;
    uint fcc_start13=10851; uint fcc_start14=16099; uint fcc_start15=34369; uint fcc_start16=67537;
    uint fcc_start17=163373+1;
    if(aurostd::substring2bool(label,"fccnAt16")) {ss.clear();ss.str(std::string());for(uint i=fcc_start16;i<=fcc_start17-1;i++) if(i<fcc_start17-1) ss<<fcc_str<<i<<","; else ss<<fcc_str<<i;aurostd::StringSubst(label,"fccnAt16",ss.str());}
    if(aurostd::substring2bool(label,"fccnAt15")) {ss.clear();ss.str(std::string());for(uint i=fcc_start15;i<=fcc_start16-1;i++) if(i<fcc_start16-1) ss<<fcc_str<<i<<","; else ss<<fcc_str<<i;aurostd::StringSubst(label,"fccnAt15",ss.str());}
    if(aurostd::substring2bool(label,"fccnAt14")) {ss.clear();ss.str(std::string());for(uint i=fcc_start14;i<=fcc_start15-1;i++) if(i<fcc_start15-1) ss<<fcc_str<<i<<","; else ss<<fcc_str<<i;aurostd::StringSubst(label,"fccnAt14",ss.str());}
    if(aurostd::substring2bool(label,"fccnAt13")) {ss.clear();ss.str(std::string());for(uint i=fcc_start13;i<=fcc_start14-1;i++) if(i<fcc_start14-1) ss<<fcc_str<<i<<","; else ss<<fcc_str<<i;aurostd::StringSubst(label,"fccnAt13",ss.str());}
    if(aurostd::substring2bool(label,"fccnAt12")) {ss.clear();ss.str(std::string());for(uint i=fcc_start12;i<=fcc_start13-1;i++) if(i<fcc_start13-1) ss<<fcc_str<<i<<","; else ss<<fcc_str<<i;aurostd::StringSubst(label,"fccnAt12",ss.str());}
    if(aurostd::substring2bool(label,"fccnAt11")) {ss.clear();ss.str(std::string());for(uint i=fcc_start11;i<=fcc_start12-1;i++) if(i<fcc_start12-1) ss<<fcc_str<<i<<","; else ss<<fcc_str<<i;aurostd::StringSubst(label,"fccnAt11",ss.str());}
    if(aurostd::substring2bool(label,"fccnAt10")) {ss.clear();ss.str(std::string());for(uint i=fcc_start10;i<=fcc_start11-1;i++) if(i<fcc_start11-1) ss<<fcc_str<<i<<","; else ss<<fcc_str<<i;aurostd::StringSubst(label,"fccnAt10",ss.str());}
    if(aurostd::substring2bool(label,"fccnAt9")) {ss.clear();ss.str(std::string());for(uint i=fcc_start9;i<=fcc_start10-1;i++) if(i<fcc_start10-1) ss<<fcc_str<<i<<","; else ss<<fcc_str<<i;aurostd::StringSubst(label,"fccnAt9",ss.str());}
    if(aurostd::substring2bool(label,"fccnAt8")) {ss.clear();ss.str(std::string());for(uint i=fcc_start8;i<=fcc_start9-1;i++) if(i<fcc_start9-1) ss<<fcc_str<<i<<","; else ss<<fcc_str<<i;aurostd::StringSubst(label,"fccnAt8",ss.str());}
    if(aurostd::substring2bool(label,"fccnAt7")) {ss.clear();ss.str(std::string());for(uint i=fcc_start7;i<=fcc_start8-1;i++) if(i<fcc_start8-1) ss<<fcc_str<<i<<","; else ss<<fcc_str<<i;aurostd::StringSubst(label,"fccnAt7",ss.str());}
    if(aurostd::substring2bool(label,"fccnAt6")) {ss.clear();ss.str(std::string());for(uint i=fcc_start6;i<=fcc_start7-1;i++) if(i<fcc_start7-1) ss<<fcc_str<<i<<","; else ss<<fcc_str<<i;aurostd::StringSubst(label,"fccnAt6",ss.str());}
    if(aurostd::substring2bool(label,"fccnAt5")) {ss.clear();ss.str(std::string());for(uint i=fcc_start5;i<=fcc_start6-1;i++) if(i<fcc_start6-1) ss<<fcc_str<<i<<","; else ss<<fcc_str<<i;aurostd::StringSubst(label,"fccnAt5",ss.str());}
    if(aurostd::substring2bool(label,"fccnAt4")) {ss.clear();ss.str(std::string());for(uint i=fcc_start4;i<=fcc_start5-1;i++) if(i<fcc_start5-1) ss<<fcc_str<<i<<","; else ss<<fcc_str<<i;aurostd::StringSubst(label,"fccnAt4",ss.str());}
    if(aurostd::substring2bool(label,"fccnAt3")) {ss.clear();ss.str(std::string());for(uint i=fcc_start3;i<=fcc_start4-1;i++) if(i<fcc_start4-1) ss<<fcc_str<<i<<","; else ss<<fcc_str<<i;aurostd::StringSubst(label,"fccnAt3",ss.str());}
    if(aurostd::substring2bool(label,"fccnAt2")) {ss.clear();ss.str(std::string());for(uint i=fcc_start2;i<=fcc_start3-1;i++) if(i<fcc_start3-1) ss<<fcc_str<<i<<","; else ss<<fcc_str<<i;aurostd::StringSubst(label,"fccnAt2",ss.str());}
    if(aurostd::substring2bool(label,"fccnAt1")) {ss.clear();ss.str(std::string());for(uint i=fcc_start1;i<=fcc_start2-1;i++) if(i<fcc_start2-1) ss<<fcc_str<<i<<","; else ss<<fcc_str<<i;aurostd::StringSubst(label,"fccnAt1",ss.str());}
    // BCC
    if(aurostd::substring2bool(label,"bccnAt1-16")) aurostd::StringSubst(label,"bccnAt1-16","bccnAt1,bccnAt2,bccnAt3,bccnAt4,bccnAt5,bccnAt6,bccnAt7,bccnAt8,bccnAt9,bccnAt10,bccnAt11,bccnAt12,bccnAt13,bccnAt14,bccnAt15,bccnAt16");
    if(aurostd::substring2bool(label,"bccnAt1-15")) aurostd::StringSubst(label,"bccnAt1-15","bccnAt1,bccnAt2,bccnAt3,bccnAt4,bccnAt5,bccnAt6,bccnAt7,bccnAt8,bccnAt9,bccnAt10,bccnAt11,bccnAt12,bccnAt13,bccnAt14,bccnAt15");
    if(aurostd::substring2bool(label,"bccnAt1-14")) aurostd::StringSubst(label,"bccnAt1-14","bccnAt1,bccnAt2,bccnAt3,bccnAt4,bccnAt5,bccnAt6,bccnAt7,bccnAt8,bccnAt9,bccnAt10,bccnAt11,bccnAt12,bccnAt13,bccnAt14");
    if(aurostd::substring2bool(label,"bccnAt1-13")) aurostd::StringSubst(label,"bccnAt1-13","bccnAt1,bccnAt2,bccnAt3,bccnAt4,bccnAt5,bccnAt6,bccnAt7,bccnAt8,bccnAt9,bccnAt10,bccnAt11,bccnAt12,bccnAt13");
    if(aurostd::substring2bool(label,"bccnAt1-12")) aurostd::StringSubst(label,"bccnAt1-12","bccnAt1,bccnAt2,bccnAt3,bccnAt4,bccnAt5,bccnAt6,bccnAt7,bccnAt8,bccnAt9,bccnAt10,bccnAt11,bccnAt12");
    if(aurostd::substring2bool(label,"bccnAt1-11")) aurostd::StringSubst(label,"bccnAt1-11","bccnAt1,bccnAt2,bccnAt3,bccnAt4,bccnAt5,bccnAt6,bccnAt7,bccnAt8,bccnAt9,bccnAt10,bccnAt11");
    if(aurostd::substring2bool(label,"bccnAt1-10")) aurostd::StringSubst(label,"bccnAt1-10","bccnAt1,bccnAt2,bccnAt3,bccnAt4,bccnAt5,bccnAt6,bccnAt7,bccnAt8,bccnAt9,bccnAt10");
    if(aurostd::substring2bool(label,"bccnAt1-9")) aurostd::StringSubst(label,"bccnAt1-9","bccnAt1,bccnAt2,bccnAt3,bccnAt4,bccnAt5,bccnAt6,bccnAt7,bccnAt8,bccnAt9");
    if(aurostd::substring2bool(label,"bccnAt1-8")) aurostd::StringSubst(label,"bccnAt1-8","bccnAt1,bccnAt2,bccnAt3,bccnAt4,bccnAt5,bccnAt6,bccnAt7,bccnAt8");
    if(aurostd::substring2bool(label,"bccnAt1-7")) aurostd::StringSubst(label,"bccnAt1-7","bccnAt1,bccnAt2,bccnAt3,bccnAt4,bccnAt5,bccnAt6,bccnAt7");
    if(aurostd::substring2bool(label,"bccnAt1-6")) aurostd::StringSubst(label,"bccnAt1-6","bccnAt1,bccnAt2,bccnAt3,bccnAt4,bccnAt5,bccnAt6");
    if(aurostd::substring2bool(label,"bccnAt1-5")) aurostd::StringSubst(label,"bccnAt1-5","bccnAt1,bccnAt2,bccnAt3,bccnAt4,bccnAt5");
    if(aurostd::substring2bool(label,"bccnAt1-4")) aurostd::StringSubst(label,"bccnAt1-4","bccnAt1,bccnAt2,bccnAt3,bccnAt4");
    if(aurostd::substring2bool(label,"bccnAt1-3")) aurostd::StringSubst(label,"bccnAt1-3","bccnAt1,bccnAt2,bccnAt3");
    if(aurostd::substring2bool(label,"bccnAt1-2")) aurostd::StringSubst(label,"bccnAt1-2","bccnAt1,bccnAt2");
    if(aurostd::substring2bool(label,"bccnAt1-1")) aurostd::StringSubst(label,"bccnAt1-1","bccnAt1");
    string bcc_str="b";
    uint bcc_start1=1; uint bcc_start2=3; uint bcc_start3=5; uint bcc_start4=11;
    uint bcc_start5=30; uint bcc_start6=58; uint bcc_start7=138; uint bcc_start8=242;
    uint bcc_start9=632; uint bcc_start10=1136; uint bcc_start11=2347; uint bcc_start12=3711;
    uint bcc_start13=10851; uint bcc_start14=16099; uint bcc_start15=34369; uint bcc_start16=67537;
    uint bcc_start17=163373+1;
    if(aurostd::substring2bool(label,"bccnAt16")) {ss.clear();ss.str(std::string());for(uint i=bcc_start16;i<=bcc_start17-1;i++) if(i<bcc_start17-1) ss<<bcc_str<<i<<","; else ss<<bcc_str<<i;aurostd::StringSubst(label,"bccnAt16",ss.str());}
    if(aurostd::substring2bool(label,"bccnAt15")) {ss.clear();ss.str(std::string());for(uint i=bcc_start15;i<=bcc_start16-1;i++) if(i<bcc_start16-1) ss<<bcc_str<<i<<","; else ss<<bcc_str<<i;aurostd::StringSubst(label,"bccnAt15",ss.str());}
    if(aurostd::substring2bool(label,"bccnAt14")) {ss.clear();ss.str(std::string());for(uint i=bcc_start14;i<=bcc_start15-1;i++) if(i<bcc_start15-1) ss<<bcc_str<<i<<","; else ss<<bcc_str<<i;aurostd::StringSubst(label,"bccnAt14",ss.str());}
    if(aurostd::substring2bool(label,"bccnAt13")) {ss.clear();ss.str(std::string());for(uint i=bcc_start13;i<=bcc_start14-1;i++) if(i<bcc_start14-1) ss<<bcc_str<<i<<","; else ss<<bcc_str<<i;aurostd::StringSubst(label,"bccnAt13",ss.str());}
    if(aurostd::substring2bool(label,"bccnAt12")) {ss.clear();ss.str(std::string());for(uint i=bcc_start12;i<=bcc_start13-1;i++) if(i<bcc_start13-1) ss<<bcc_str<<i<<","; else ss<<bcc_str<<i;aurostd::StringSubst(label,"bccnAt12",ss.str());}
    if(aurostd::substring2bool(label,"bccnAt11")) {ss.clear();ss.str(std::string());for(uint i=bcc_start11;i<=bcc_start12-1;i++) if(i<bcc_start12-1) ss<<bcc_str<<i<<","; else ss<<bcc_str<<i;aurostd::StringSubst(label,"bccnAt11",ss.str());}
    if(aurostd::substring2bool(label,"bccnAt10")) {ss.clear();ss.str(std::string());for(uint i=bcc_start10;i<=bcc_start11-1;i++) if(i<bcc_start11-1) ss<<bcc_str<<i<<","; else ss<<bcc_str<<i;aurostd::StringSubst(label,"bccnAt10",ss.str());}
    if(aurostd::substring2bool(label,"bccnAt9")) {ss.clear();ss.str(std::string());for(uint i=bcc_start9;i<=bcc_start10-1;i++) if(i<bcc_start10-1) ss<<bcc_str<<i<<","; else ss<<bcc_str<<i;aurostd::StringSubst(label,"bccnAt9",ss.str());}
    if(aurostd::substring2bool(label,"bccnAt8")) {ss.clear();ss.str(std::string());for(uint i=bcc_start8;i<=bcc_start9-1;i++) if(i<bcc_start9-1) ss<<bcc_str<<i<<","; else ss<<bcc_str<<i;aurostd::StringSubst(label,"bccnAt8",ss.str());}
    if(aurostd::substring2bool(label,"bccnAt7")) {ss.clear();ss.str(std::string());for(uint i=bcc_start7;i<=bcc_start8-1;i++) if(i<bcc_start8-1) ss<<bcc_str<<i<<","; else ss<<bcc_str<<i;aurostd::StringSubst(label,"bccnAt7",ss.str());}
    if(aurostd::substring2bool(label,"bccnAt6")) {ss.clear();ss.str(std::string());for(uint i=bcc_start6;i<=bcc_start7-1;i++) if(i<bcc_start7-1) ss<<bcc_str<<i<<","; else ss<<bcc_str<<i;aurostd::StringSubst(label,"bccnAt6",ss.str());}
    if(aurostd::substring2bool(label,"bccnAt5")) {ss.clear();ss.str(std::string());for(uint i=bcc_start5;i<=bcc_start6-1;i++) if(i<bcc_start6-1) ss<<bcc_str<<i<<","; else ss<<bcc_str<<i;aurostd::StringSubst(label,"bccnAt5",ss.str());}
    if(aurostd::substring2bool(label,"bccnAt4")) {ss.clear();ss.str(std::string());for(uint i=bcc_start4;i<=bcc_start5-1;i++) if(i<bcc_start5-1) ss<<bcc_str<<i<<","; else ss<<bcc_str<<i;aurostd::StringSubst(label,"bccnAt4",ss.str());}
    if(aurostd::substring2bool(label,"bccnAt3")) {ss.clear();ss.str(std::string());for(uint i=bcc_start3;i<=bcc_start4-1;i++) if(i<bcc_start4-1) ss<<bcc_str<<i<<","; else ss<<bcc_str<<i;aurostd::StringSubst(label,"bccnAt3",ss.str());}
    if(aurostd::substring2bool(label,"bccnAt2")) {ss.clear();ss.str(std::string());for(uint i=bcc_start2;i<=bcc_start3-1;i++) if(i<bcc_start3-1) ss<<bcc_str<<i<<","; else ss<<bcc_str<<i;aurostd::StringSubst(label,"bccnAt2",ss.str());}
    if(aurostd::substring2bool(label,"bccnAt1")) {ss.clear();ss.str(std::string());for(uint i=bcc_start1;i<=bcc_start2-1;i++) if(i<bcc_start2-1) ss<<bcc_str<<i<<","; else ss<<bcc_str<<i;aurostd::StringSubst(label,"bccnAt1",ss.str());}
    // SC
    if(aurostd::substring2bool(label,"scnAt1-16")) aurostd::StringSubst(label,"scnAt1-16","scnAt1,scnAt2,scnAt3,scnAt4,scnAt5,scnAt6,scnAt7,scnAt8,scnAt9,scnAt10,scnAt11,scnAt12,scnAt13,scnAt14,scnAt15,scnAt16");
    if(aurostd::substring2bool(label,"scnAt1-15")) aurostd::StringSubst(label,"scnAt1-15","scnAt1,scnAt2,scnAt3,scnAt4,scnAt5,scnAt6,scnAt7,scnAt8,scnAt9,scnAt10,scnAt11,scnAt12,scnAt13,scnAt14,scnAt15");
    if(aurostd::substring2bool(label,"scnAt1-14")) aurostd::StringSubst(label,"scnAt1-14","scnAt1,scnAt2,scnAt3,scnAt4,scnAt5,scnAt6,scnAt7,scnAt8,scnAt9,scnAt10,scnAt11,scnAt12,scnAt13,scnAt14");
    if(aurostd::substring2bool(label,"scnAt1-13")) aurostd::StringSubst(label,"scnAt1-13","scnAt1,scnAt2,scnAt3,scnAt4,scnAt5,scnAt6,scnAt7,scnAt8,scnAt9,scnAt10,scnAt11,scnAt12,scnAt13");
    if(aurostd::substring2bool(label,"scnAt1-12")) aurostd::StringSubst(label,"scnAt1-12","scnAt1,scnAt2,scnAt3,scnAt4,scnAt5,scnAt6,scnAt7,scnAt8,scnAt9,scnAt10,scnAt11,scnAt12");
    if(aurostd::substring2bool(label,"scnAt1-11")) aurostd::StringSubst(label,"scnAt1-11","scnAt1,scnAt2,scnAt3,scnAt4,scnAt5,scnAt6,scnAt7,scnAt8,scnAt9,scnAt10,scnAt11");
    if(aurostd::substring2bool(label,"scnAt1-10")) aurostd::StringSubst(label,"scnAt1-10","scnAt1,scnAt2,scnAt3,scnAt4,scnAt5,scnAt6,scnAt7,scnAt8,scnAt9,scnAt10");
    if(aurostd::substring2bool(label,"scnAt1-9")) aurostd::StringSubst(label,"scnAt1-9","scnAt1,scnAt2,scnAt3,scnAt4,scnAt5,scnAt6,scnAt7,scnAt8,scnAt9");
    if(aurostd::substring2bool(label,"scnAt1-8")) aurostd::StringSubst(label,"scnAt1-8","scnAt1,scnAt2,scnAt3,scnAt4,scnAt5,scnAt6,scnAt7,scnAt8");
    if(aurostd::substring2bool(label,"scnAt1-7")) aurostd::StringSubst(label,"scnAt1-7","scnAt1,scnAt2,scnAt3,scnAt4,scnAt5,scnAt6,scnAt7");
    if(aurostd::substring2bool(label,"scnAt1-6")) aurostd::StringSubst(label,"scnAt1-6","scnAt1,scnAt2,scnAt3,scnAt4,scnAt5,scnAt6");
    if(aurostd::substring2bool(label,"scnAt1-5")) aurostd::StringSubst(label,"scnAt1-5","scnAt1,scnAt2,scnAt3,scnAt4,scnAt5");
    if(aurostd::substring2bool(label,"scnAt1-4")) aurostd::StringSubst(label,"scnAt1-4","scnAt1,scnAt2,scnAt3,scnAt4");
    if(aurostd::substring2bool(label,"scnAt1-3")) aurostd::StringSubst(label,"scnAt1-3","scnAt1,scnAt2,scnAt3");
    if(aurostd::substring2bool(label,"scnAt1-2")) aurostd::StringSubst(label,"scnAt1-2","scnAt1,scnAt2");
    if(aurostd::substring2bool(label,"scnAt1-1")) aurostd::StringSubst(label,"scnAt1-1","scnAt1");
    string sc_str="s";
    uint sc_start1=1; uint sc_start2=3; uint sc_start3=6; uint sc_start4=12;
    uint sc_start5=36; uint sc_start6=64; uint sc_start7=168; uint sc_start8=272;
    uint sc_start9=763; uint sc_start10=1267; uint sc_start11=2761; uint sc_start12=4125;
    uint sc_start13=12859; uint sc_start14=18107; uint sc_start15=39717; uint sc_start16=72885;
    uint sc_start17=188729+1;
    if(aurostd::substring2bool(label,"scnAt16")) {ss.clear();ss.str(std::string());for(uint i=sc_start16;i<=sc_start17-1;i++) if(i<sc_start17-1) ss<<sc_str<<i<<","; else ss<<sc_str<<i;aurostd::StringSubst(label,"scnAt16",ss.str());}
    if(aurostd::substring2bool(label,"scnAt15")) {ss.clear();ss.str(std::string());for(uint i=sc_start15;i<=sc_start16-1;i++) if(i<sc_start16-1) ss<<sc_str<<i<<","; else ss<<sc_str<<i;aurostd::StringSubst(label,"scnAt15",ss.str());}
    if(aurostd::substring2bool(label,"scnAt14")) {ss.clear();ss.str(std::string());for(uint i=sc_start14;i<=sc_start15-1;i++) if(i<sc_start15-1) ss<<sc_str<<i<<","; else ss<<sc_str<<i;aurostd::StringSubst(label,"scnAt14",ss.str());}
    if(aurostd::substring2bool(label,"scnAt13")) {ss.clear();ss.str(std::string());for(uint i=sc_start13;i<=sc_start14-1;i++) if(i<sc_start14-1) ss<<sc_str<<i<<","; else ss<<sc_str<<i;aurostd::StringSubst(label,"scnAt13",ss.str());}
    if(aurostd::substring2bool(label,"scnAt12")) {ss.clear();ss.str(std::string());for(uint i=sc_start12;i<=sc_start13-1;i++) if(i<sc_start13-1) ss<<sc_str<<i<<","; else ss<<sc_str<<i;aurostd::StringSubst(label,"scnAt12",ss.str());}
    if(aurostd::substring2bool(label,"scnAt11")) {ss.clear();ss.str(std::string());for(uint i=sc_start11;i<=sc_start12-1;i++) if(i<sc_start12-1) ss<<sc_str<<i<<","; else ss<<sc_str<<i;aurostd::StringSubst(label,"scnAt11",ss.str());}
    if(aurostd::substring2bool(label,"scnAt10")) {ss.clear();ss.str(std::string());for(uint i=sc_start10;i<=sc_start11-1;i++) if(i<sc_start11-1) ss<<sc_str<<i<<","; else ss<<sc_str<<i;aurostd::StringSubst(label,"scnAt10",ss.str());}
    if(aurostd::substring2bool(label,"scnAt9")) {ss.clear();ss.str(std::string());for(uint i=sc_start9;i<=sc_start10-1;i++) if(i<sc_start10-1) ss<<sc_str<<i<<","; else ss<<sc_str<<i;aurostd::StringSubst(label,"scnAt9",ss.str());}
    if(aurostd::substring2bool(label,"scnAt8")) {ss.clear();ss.str(std::string());for(uint i=sc_start8;i<=sc_start9-1;i++) if(i<sc_start9-1) ss<<sc_str<<i<<","; else ss<<sc_str<<i;aurostd::StringSubst(label,"scnAt8",ss.str());}
    if(aurostd::substring2bool(label,"scnAt7")) {ss.clear();ss.str(std::string());for(uint i=sc_start7;i<=sc_start8-1;i++) if(i<sc_start8-1) ss<<sc_str<<i<<","; else ss<<sc_str<<i;aurostd::StringSubst(label,"scnAt7",ss.str());}
    if(aurostd::substring2bool(label,"scnAt6")) {ss.clear();ss.str(std::string());for(uint i=sc_start6;i<=sc_start7-1;i++) if(i<sc_start7-1) ss<<sc_str<<i<<","; else ss<<sc_str<<i;aurostd::StringSubst(label,"scnAt6",ss.str());}
    if(aurostd::substring2bool(label,"scnAt5")) {ss.clear();ss.str(std::string());for(uint i=sc_start5;i<=sc_start6-1;i++) if(i<sc_start6-1) ss<<sc_str<<i<<","; else ss<<sc_str<<i;aurostd::StringSubst(label,"scnAt5",ss.str());}
    if(aurostd::substring2bool(label,"scnAt4")) {ss.clear();ss.str(std::string());for(uint i=sc_start4;i<=sc_start5-1;i++) if(i<sc_start5-1) ss<<sc_str<<i<<","; else ss<<sc_str<<i;aurostd::StringSubst(label,"scnAt4",ss.str());}
    if(aurostd::substring2bool(label,"scnAt3")) {ss.clear();ss.str(std::string());for(uint i=sc_start3;i<=sc_start4-1;i++) if(i<sc_start4-1) ss<<sc_str<<i<<","; else ss<<sc_str<<i;aurostd::StringSubst(label,"scnAt3",ss.str());}
    if(aurostd::substring2bool(label,"scnAt2")) {ss.clear();ss.str(std::string());for(uint i=sc_start2;i<=sc_start3-1;i++) if(i<sc_start3-1) ss<<sc_str<<i<<","; else ss<<sc_str<<i;aurostd::StringSubst(label,"scnAt2",ss.str());}
    if(aurostd::substring2bool(label,"scnAt1")) {ss.clear();ss.str(std::string());for(uint i=sc_start1;i<=sc_start2-1;i++) if(i<sc_start2-1) ss<<sc_str<<i<<","; else ss<<sc_str<<i;aurostd::StringSubst(label,"scnAt1",ss.str());}
    // HCP
    if(aurostd::substring2bool(label,"hcpnAt1-8")) aurostd::StringSubst(label,"hcpnAt1-8","hcpnAt1,hcpnAt2,hcpnAt3,hcpnAt4,hcpnAt5,hcpnAt6,hcpnAt7,hcpnAt8");
    if(aurostd::substring2bool(label,"hcpnAt1-7")) aurostd::StringSubst(label,"hcpnAt1-7","hcpnAt1,hcpnAt2,hcpnAt3,hcpnAt4,hcpnAt5,hcpnAt6,hcpnAt7");
    if(aurostd::substring2bool(label,"hcpnAt1-6")) aurostd::StringSubst(label,"hcpnAt1-6","hcpnAt1,hcpnAt2,hcpnAt3,hcpnAt4,hcpnAt5,hcpnAt6");
    if(aurostd::substring2bool(label,"hcpnAt1-5")) aurostd::StringSubst(label,"hcpnAt1-5","hcpnAt1,hcpnAt2,hcpnAt3,hcpnAt4,hcpnAt5");
    if(aurostd::substring2bool(label,"hcpnAt1-4")) aurostd::StringSubst(label,"hcpnAt1-4","hcpnAt1,hcpnAt2,hcpnAt3,hcpnAt4");
    if(aurostd::substring2bool(label,"hcpnAt1-3")) aurostd::StringSubst(label,"hcpnAt1-3","hcpnAt1,hcpnAt2,hcpnAt3");
    if(aurostd::substring2bool(label,"hcpnAt1-2")) aurostd::StringSubst(label,"hcpnAt1-2","hcpnAt1,hcpnAt2");
    if(aurostd::substring2bool(label,"hcpnAt1-1")) aurostd::StringSubst(label,"hcpnAt1-1","hcpnAt1");
    string hcp_str="h";
    uint hcp_start1=1; uint hcp_start2=4; uint hcp_start3=14; uint hcp_start4=64;
    uint hcp_start5=334; uint hcp_start6=985; uint hcp_start7=5778; uint hcp_start8=15796;
    uint hcp_start9=1643380+1;
    if(aurostd::substring2bool(label,"hcpnAt8")) {ss.clear();ss.str(std::string());for(uint i=hcp_start8;i<=hcp_start9-1;i++) if(i<hcp_start9-1) ss<<hcp_str<<i<<","; else ss<<hcp_str<<i;aurostd::StringSubst(label,"hcpnAt8",ss.str());}
    if(aurostd::substring2bool(label,"hcpnAt7")) {ss.clear();ss.str(std::string());for(uint i=hcp_start7;i<=hcp_start8-1;i++) if(i<hcp_start8-1) ss<<hcp_str<<i<<","; else ss<<hcp_str<<i;aurostd::StringSubst(label,"hcpnAt7",ss.str());}
    if(aurostd::substring2bool(label,"hcpnAt6")) {ss.clear();ss.str(std::string());for(uint i=hcp_start6;i<=hcp_start7-1;i++) if(i<hcp_start7-1) ss<<hcp_str<<i<<","; else ss<<hcp_str<<i;aurostd::StringSubst(label,"hcpnAt6",ss.str());}
    if(aurostd::substring2bool(label,"hcpnAt5")) {ss.clear();ss.str(std::string());for(uint i=hcp_start5;i<=hcp_start6-1;i++) if(i<hcp_start6-1) ss<<hcp_str<<i<<","; else ss<<hcp_str<<i;aurostd::StringSubst(label,"hcpnAt5",ss.str());}
    if(aurostd::substring2bool(label,"hcpnAt4")) {ss.clear();ss.str(std::string());for(uint i=hcp_start4;i<=hcp_start5-1;i++) if(i<hcp_start5-1) ss<<hcp_str<<i<<","; else ss<<hcp_str<<i;aurostd::StringSubst(label,"hcpnAt4",ss.str());}
    if(aurostd::substring2bool(label,"hcpnAt3")) {ss.clear();ss.str(std::string());for(uint i=hcp_start3;i<=hcp_start4-1;i++) if(i<hcp_start4-1) ss<<hcp_str<<i<<","; else ss<<hcp_str<<i;aurostd::StringSubst(label,"hcpnAt3",ss.str());}
    if(aurostd::substring2bool(label,"hcpnAt2")) {ss.clear();ss.str(std::string());for(uint i=hcp_start2;i<=hcp_start3-1;i++) if(i<hcp_start3-1) ss<<hcp_str<<i<<","; else ss<<hcp_str<<i;aurostd::StringSubst(label,"hcpnAt2",ss.str());}
    if(aurostd::substring2bool(label,"hcpnAt1")) {ss.clear();ss.str(std::string());for(uint i=hcp_start1;i<=hcp_start2-1;i++) if(i<hcp_start2-1) ss<<hcp_str<<i<<","; else ss<<hcp_str<<i;aurostd::StringSubst(label,"hcpnAt1",ss.str());}
    // SEARCH FOR SPECIES
    //  cerr << label << endl;
    return label;
}
// cat aflow_library_gus.dat | grep hcp | grep "  1   n " | head -1


string AVASP_Shortcuts_for_Ternaries(string &label) {
    // SEARCH FOR LABELS
    if(label=="T0001") label="T0001.AB2C,T0001.A2BC,T0001.ABC2";  // HEUSLERS AlCu2Mn-Heusler
    if(label=="T0002") label="T0002.AB2C,T0002.A2BC,T0002.ABC2";  // HEUSLERS Cu1Li2Sn1-antiHeusler
    if(label=="T0003") label="T0003.ABC,T0003.BCA,T0003.CAB";     // HALF HEUSLERS
    if(label=="T0004") label="T0004.AB2C,T0004.A2BC,T0004.ABC2";  // HEUSLERS Cu1Li2Sn1-antiHeusler
    if(label=="T0005") label="T0005.ABC,T0005.BAC,T0005.BCA,T0005.CBA,T0005.CAB,T0005.ACB";  // Co1Ge1Mn1_ICSD_623495
    if(label=="T0006") label="T0006.ABC,T0006.BAC,T0006.BCA,T0006.CBA,T0006.CAB,T0006.ACB";  // Co3Ge1Mn2_ICSD_52972
    if(label=="T0007") label="T0007.ABC,T0007.BAC,T0007.BCA,T0007.CBA,T0007.CAB,T0007.ACB";  // Cl1La1Se1_ICSD_425686
    if(label=="T0008") label="T0008.ABC,T0008.BAC,T0008.BCA,T0008.CBA,T0008.CAB,T0008.ACB";  // Cd1Na1P1_ICSD_44858
    if(label=="TS001") label="TS001.ABC,TS001.BCA,TS001.CAB";      // GUS SQS_L12
    if(label=="T0009") label="T0009.ABC,T0009.BAC,T0009.BCA,T0009.CBA,T0009.CAB,T0009.ACB";  // perovskite
    if(label=="T0010") label="T0010.ABC,T0010.BCA,T0010.CAB";     // C11
    if(label=="T0011") label="T0011.ABC,T0011.ACB,T0011.BAC,T0011.BCA,T0011.CAB,T0011.CBA";

    // FCC SYSTEMS
    if(label=="TFCC001") label="TFCC001.ABC";
    if(label=="TFCC002") label="TFCC002.ABC";
    if(label=="TFCC003") label="TFCC003.ABC";
    if(label=="TFCC004") label="TFCC004.ABC,TFCC004.BCA,TFCC004.CAB";
    if(label=="TFCC005") label="TFCC005.ABC,TFCC005.BCA,TFCC005.CAB";
    if(label=="TFCC006") label="TFCC006.ABC,TFCC006.BCA,TFCC006.CAB";
    if(label=="TFCC007") label="TFCC007.ABC,TFCC007.BCA,TFCC007.CAB";
    if(label=="TFCC008") label="TFCC008.ABC,TFCC008.BCA,TFCC008.CAB";
    if(label=="TFCC009") label="TFCC009.ABC,TFCC009.BCA,TFCC009.CAB";
    if(label=="TFCC010") label="TFCC010.ABC,TFCC010.BCA,TFCC010.CAB";
    if(label=="TFCC011") label="TFCC011.ABC,TFCC011.BCA,TFCC011.CAB";
    if(label=="TFCC012") label="TFCC012.ABC,TFCC012.BCA,TFCC012.CAB";
    if(label=="TFCC013") label="TFCC013.ABC,TFCC013.BCA,TFCC013.CAB";
    if(label=="TFCC014") label="TFCC014.ABC,TFCC014.BCA,TFCC014.CAB";
    if(label=="TFCC015") label="TFCC015.ABC,TFCC015.BCA,TFCC015.CAB";
    if(label=="TFCC016") label="TFCC016.ABC,TFCC016.BCA,TFCC016.CAB";

    if(label=="TFCC00X") label="TFCC001.ABC,TFCC002.ABC,TFCC003.ABC,TFCC004.ABC,TFCC004.BCA,TFCC004.CAB,TFCC005.ABC,TFCC005.BCA,TFCC005.CAB,TFCC006.ABC,TFCC006.BCA,TFCC006.CAB,TFCC007.ABC,TFCC007.BCA,TFCC007.CAB,TFCC008.ABC,TFCC008.BCA,TFCC008.CAB,TFCC009.ABC,TFCC009.BCA,TFCC009.CAB,TFCC010.ABC,TFCC010.BCA,TFCC010.CAB,TFCC011.ABC,TFCC011.BCA,TFCC011.CAB,TFCC012.ABC,TFCC012.BCA,TFCC012.CAB,TFCC013.ABC,TFCC013.BCA,TFCC013.CAB,TFCC014.ABC,TFCC014.BCA,TFCC014.CAB,TFCC015.ABC,TFCC015.BCA,TFCC015.CAB,TFCC016.ABC,TFCC016.BCA,TFCC016.CAB";

    // BCC SYSTEMS
    if(label=="TBCC001") label="TBCC001.ABC";
    if(label=="TBCC002") label="TBCC002.ABC";
    if(label=="TBCC003") label="TBCC003.ABC";
    if(label=="TBCC004") label="TBCC004.ABC,TBCC004.BCA,TBCC004.CAB";
    if(label=="TBCC005") label="TBCC005.ABC,TBCC005.BCA,TBCC005.CAB";
    if(label=="TBCC006") label="TBCC006.ABC,TBCC006.BCA,TBCC006.CAB";
    if(label=="TBCC007") label="TBCC007.ABC,TBCC007.BCA,TBCC007.CAB";
    //  if(label=="TBCC008") label="TBCC008.ABC,TBCC008.BCA,TBCC008.CAB";
    if(label=="TBCC009") label="TBCC009.ABC,TBCC009.BCA,TBCC009.CAB";
    if(label=="TBCC010") label="TBCC010.ABC,TBCC010.BCA,TBCC010.CAB";
    if(label=="TBCC011") label="TBCC011.ABC,TBCC011.BCA,TBCC011.CAB";
    if(label=="TBCC012") label="TBCC012.ABC,TBCC012.BCA,TBCC012.CAB";
    //  if(label=="TBCC013") label="TBCC013.ABC,TBCC013.BCA,TBCC013.CAB";
    if(label=="TBCC014") label="TBCC014.ABC,TBCC014.BCA,TBCC014.CAB";
    if(label=="TBCC015") label="TBCC015.ABC,TBCC015.BCA,TBCC015.CAB";
    if(label=="TBCC016") label="TBCC016.ABC,TBCC016.BCA,TBCC016.CAB";

    if(label=="TBCC00X") label="TBCC001.ABC,TBCC002.ABC,TBCC003.ABC,TBCC004.ABC,TBCC004.BCA,TBCC004.CAB,TBCC005.ABC,TBCC005.BCA,TBCC005.CAB,TBCC006.ABC,TBCC006.BCA,TBCC006.CAB,TBCC007.ABC,TBCC007.BCA,TBCC007.CAB,TBCC009.ABC,TBCC009.BCA,TBCC009.CAB,TBCC010.ABC,TBCC010.BCA,TBCC010.CAB,TBCC011.ABC,TBCC011.BCA,TBCC011.CAB,TBCC012.ABC,TBCC012.BCA,TBCC012.CAB,TBCC014.ABC,TBCC014.BCA,TBCC014.CAB,TBCC015.ABC,TBCC015.BCA,TBCC015.CAB,TBCC016.ABC,TBCC016.BCA,TBCC016.CAB";

    if(label=="Q0001all") label="Q0001.ABCD,Q0001.ABDC,Q0001.ACDB,Q0001.ACBD,Q0001.ADBC,Q0001.ADCB,Q0001.BCDA,Q0001.BCAD,Q0001.BDAC,Q0001.BDCA,Q0001.BACD,Q0001.BADC,Q0001.CDAB,Q0001.CDBA,Q0001.CABD,Q0001.CADB,Q0001.CBDA,Q0001.CBAD,Q0001.DABC,Q0001.DACB,Q0001.DBCA,Q0001.DBAC,Q0001.DCAB,Q0001.DCBA";  // Elpasolite

    return label;
}


#endif

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *                                                                         *
// ***************************************************************************

