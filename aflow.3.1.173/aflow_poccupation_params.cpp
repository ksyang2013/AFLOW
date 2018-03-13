// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *                                                                         *
// ***************************************************************************
// functions written by KESONG YANG
// 2010-2012: kesong.yang@gmail.com

#ifndef _AFLOW_POCCUPATION_PARAMS_CPP_
#define _AFLOW_POCCUPATION_PARAMS_CPP_

#include "aflow.h"

//References: 
//"UFF, a full periodic table force field for molecular mechanics and molecular dynamics simulations"
//A. K. Rappe, C. J. Casewit, K. S. Colwell, W. A. Goddard III, W. M. Skiff, JACS, 114, 10024 (1992)

// ***************************************************************************
// string pocc::ReturnAtomSpecies(string& atom)
// ***************************************************************************
// input for GULP program
namespace pocc {
  string ReturnAtomSpecies(string atom) {
    if(atom.compare("H")==0)    {return  "H  core";}
    else if(atom.compare("He")==0)   {return  "He core";}
    else if(atom.compare("Li")==0)   {return  "Li core";}
    else if(atom.compare("Be")==0)   {return  "Be core";}
    else if(atom.compare("B")==0)    {return  "B  core";}
    else if(atom.compare("C")==0)    {return  "C  core";}
    else if(atom.compare("N")==0)    {return  "N  core";}
    else if(atom.compare("O")==0)    {return  "O  core";}
    else if(atom.compare("F")==0)    {return  "F  core";}
    else if(atom.compare("Ne")==0)   {return  "Ne core";}
    else if(atom.compare("Na")==0)   {return  "Na core";}
    else if(atom.compare("Mg")==0)   {return  "Mg core";}
    else if(atom.compare("Al")==0)   {return  "Al core";}
    else if(atom.compare("Si")==0)   {return  "Si core";}
    else if(atom.compare("P")==0)    {return  "P  core";}
    else if(atom.compare("S")==0)    {return  "S  core";}
    else if(atom.compare("Cl")==0)   {return  "Cl core";}
    else if(atom.compare("Ar")==0)   {return  "Ar core";}
    else if(atom.compare("K")==0)    {return  "K  core";}
    else if(atom.compare("Ca")==0)   {return  "Ca core";}
    else if(atom.compare("Sc")==0)   {return  "Sc core";}
    else if(atom.compare("Ti")==0)   {return  "Ti core";}
    else if(atom.compare("V")==0)    {return  "V  core";}
    else if(atom.compare("Cr")==0)   {return  "Cr core";}
    else if(atom.compare("Mn")==0)   {return  "Mn core";}
    else if(atom.compare("Fe")==0)   {return  "Fe core";}
    else if(atom.compare("Co")==0)   {return  "Co core";}
    else if(atom.compare("Ni")==0)   {return  "Ni core";}
    else if(atom.compare("Cu")==0)   {return  "Cu core";}
    else if(atom.compare("Zn")==0)   {return  "Zn core";}
    else if(atom.compare("Ga")==0)   {return  "Ga core";}
    else if(atom.compare("Ge")==0)   {return  "Ge core";}
    else if(atom.compare("As")==0)   {return  "As core";}
    else if(atom.compare("Se")==0)   {return  "Se core";}
    else if(atom.compare("Br")==0)   {return  "Br core";}
    else if(atom.compare("Kr")==0)   {return  "Kr core";}
    else if(atom.compare("Rb")==0)   {return  "Rb core";}
    else if(atom.compare("Sr")==0)   {return  "Sr core";}
    else if(atom.compare("Y")==0)    {return  "Y  core";}
    else if(atom.compare("Zr")==0)   {return  "Zr core";}
    else if(atom.compare("Nb")==0)   {return  "Nb core";}
    else if(atom.compare("Mo")==0)   {return  "Mo core";}
    else if(atom.compare("Tc")==0)   {return  "Tc core";}
    else if(atom.compare("Ru")==0)   {return  "Ru core";}
    else if(atom.compare("Rh")==0)   {return  "Rh core";}
    else if(atom.compare("Pd")==0)   {return  "Pd core";}
    else if(atom.compare("Ag")==0)   {return  "Ag core";}
    else if(atom.compare("Cd")==0)   {return  "Cd core";}
    else if(atom.compare("In")==0)   {return  "In core";}
    else if(atom.compare("Sn")==0)   {return  "Sn core";}
    else if(atom.compare("Sb")==0)   {return  "Sb core";}
    else if(atom.compare("Te")==0)   {return  "Te core";}
    else if(atom.compare("I")==0)    {return  "I  core";}
    else if(atom.compare("Xe")==0)   {return  "Xe core";}
    else if(atom.compare("Cs")==0)   {return  "Cs core";}
    else if(atom.compare("Ba")==0)   {return  "Ba core";}
    else if(atom.compare("La")==0)   {return  "La core";}
    else if(atom.compare("Ce")==0)   {return  "Ce core";}
    else if(atom.compare("Pr")==0)   {return  "Pr core";}
    else if(atom.compare("Nd")==0)   {return  "Nd core";}
    else if(atom.compare("Pm")==0)   {return  "Pm core";}
    else if(atom.compare("Sm")==0)   {return  "Sm core";}
    else if(atom.compare("Eu")==0)   {return  "Eu core";}
    else if(atom.compare("Gd")==0)   {return  "Gd core";}
    else if(atom.compare("Tb")==0)   {return  "Tb core";}
    else if(atom.compare("Dy")==0)   {return  "Dy core";}
    else if(atom.compare("Ho")==0)   {return  "Ho core";}
    else if(atom.compare("Er")==0)   {return  "Er core";}
    else if(atom.compare("Tm")==0)   {return  "Tm core";}
    else if(atom.compare("Yb")==0)   {return  "Yb core";}
    else if(atom.compare("Lu")==0)   {return  "Lu core";}
    else if(atom.compare("Hf")==0)   {return  "Hf core";}
    else if(atom.compare("Ta")==0)   {return  "Ta core";}
    else if(atom.compare("W")==0)    {return  "W  core";}
    else if(atom.compare("Re")==0)   {return  "Re core";}
    else if(atom.compare("Os")==0)   {return  "Os core";}
    else if(atom.compare("Ir")==0)   {return  "Ir core";}
    else if(atom.compare("Pt")==0)   {return  "Pt core";}
    else if(atom.compare("Au")==0)   {return  "Au core";}
    else if(atom.compare("Hg")==0)   {return  "Hg core";}
    else if(atom.compare("Tl")==0)   {return  "Tl core";}
    else if(atom.compare("Pb")==0)   {return  "Pb core";}
    else if(atom.compare("Bi")==0)   {return  "Bi core";}
    else if(atom.compare("Po")==0)   {return  "Po core";}
    else if(atom.compare("At")==0)   {return  "At core";}
    else if(atom.compare("Rn")==0)   {return  "Rn core";}
    else if(atom.compare("Fr")==0)   {return  "Fr core";}
    else if(atom.compare("Ra")==0)   {return  "Ra core";}
    else if(atom.compare("Ac")==0)   {return  "Ac core";}
    else if(atom.compare("Th")==0)   {return  "Th core";}
    else if(atom.compare("Pa")==0)   {return  "Pa core";}
    else if(atom.compare("U")==0)    {return  "U  core";}
    else if(atom.compare("Np")==0)   {return  "Np core";}
    else if(atom.compare("Pu")==0)   {return  "Pu core";}
    else if(atom.compare("Am")==0)   {return  "Am core";}
    else if(atom.compare("Cm")==0)   {return  "Cm core";}
    else if(atom.compare("Bk")==0)   {return  "Bk core";}
    else if(atom.compare("Cf")==0)   {return  "Cf core";}
    else if(atom.compare("Es")==0)   {return  "Es core";}
    else if(atom.compare("Fm")==0)   {return  "Fm core";}
    else if(atom.compare("Md")==0)   {return  "Md core";}
    else if(atom.compare("No")==0)   {return  "No core";}
    else if(atom.compare("Lr")==0)   {return  "Lr core";}
    else {return "Unknown";}
  }
} // namespace pocc 
// ***************************************************************************
// string pocc::ReturnAtomSpecies(string& atom)
// ***************************************************************************
// Potential for GULP, only consider sp1
namespace pocc {
  string ReturnAtomSpeciesPotential(string atom) {
    if(atom.compare("H")==0)        {return  "H    0.354   180.000 2.886   0.044   12.000 0.7120 0 0.000  0.0  0.0000  4.528";}      
    else if(atom.compare("He")==0)  {return  "He   0.849   90.000  2.362   0.056   15.240 0.0972 4 0.000  0.0  0.0000  9.660";}      
    else if(atom.compare("Li")==0)  {return  "Li   1.336   180.000 2.451   0.025   12.000 1.0255 0 0.000  0.0  0.0000  3.006";}      
    else if(atom.compare("Be")==0)  {return  "Be   1.074   109.470 2.745   0.085   12.000 1.5650 3 0.000  0.0  0.0000  4.877";}      
    else if(atom.compare("B")==0)   {return  "B    0.838   109.470 4.083   0.180   12.052 1.7550 3 0.000  0.0  0.0000  4.068";}      
    else if(atom.compare("C")==0)   {return  "C    0.706   180.000 3.851   0.105   12.730 1.9120 3 2.119  0.0  0.0000  5.343";}      
    else if(atom.compare("N")==0)   {return  "N    0.656   180.000 3.660   0.069   13.407 2.5438 3 0.450  0.0 61.2230  6.899";}      
    else if(atom.compare("O")==0)   {return  "O    0.639   180.000 3.500   0.060   14.085 2.2998 3 0.018  0.0  0.0000  8.741";}      
    else if(atom.compare("F")==0)   {return  "F    0.668   180.000 3.364   0.050   14.762 1.7350 0 0.000  0.0  0.0000 10.874";}      
    else if(atom.compare("Ne")==0)  {return  "Ne   0.920   90.000  3.243   0.042   15.440 0.1944 4 0.000  0.0  0.0000 11.040";}      
    else if(atom.compare("Na")==0)  {return  "Na   1.539   180.000 2.983   0.030   12.000 1.0809 0 0.000  0.0  0.0000  2.843";}      
    else if(atom.compare("Mg")==0)  {return  "Mg   1.421   109.470 3.021   0.111   12.000 1.7866 3 0.000  0.0  0.0000  3.951";}      
    else if(atom.compare("Al")==0)  {return  "Al   1.244   109.470 4.499   0.505   11.278 1.7924 3 0.000  0.0  0.0000  3.041";}      
    else if(atom.compare("Si")==0)  {return  "Si   1.117   109.470 4.295   0.402   12.175 2.3232 3 1.225  0.0  0.0000  4.168";}      
    else if(atom.compare("P")==0)   {return  "P    1.101   93.800  4.147   0.305   13.072 2.8627 3 2.400 22.0 84.4339  5.463";}      
    else if(atom.compare("S")==0)   {return  "S    1.064   92.1000 4.035   0.274   13.969 2.7032 3 0.484  0.0  0.0000  6.928";}      
    else if(atom.compare("Cl")==0)  {return  "Cl   1.044   180.000 3.947   0.227   14.866 2.3484 0 0.000  0.0  0.0000  8.564";}      
    else if(atom.compare("Ar")==0)  {return  "Ar   1.032   90.000  3.868   0.185   15.763 0.2994 4 0.000  0.0  0.0000  9.465";}      
    else if(atom.compare("K")==0)   {return  "K    1.953   180.000 3.812   0.035   12.000 1.1645 0 0.000  0.0  0.0000  2.421";}      
    else if(atom.compare("Ca")==0)  {return  "Ca   1.761   90.000  3.399   0.238   12.000 2.1414 6 0.000  0.0  0.0000  3.231";}      
    else if(atom.compare("Sc")==0)  {return  "Sc   1.513   109.470 3.295   0.019   12.000 2.5924 3 0.000  0.0  0.0000  3.395";}      
    else if(atom.compare("Ti")==0)  {return  "Ti   1.412   109.470 3.175   0.017   12.000 2.6595 3 0.000  0.0  0.0000  3.470";}      
    else if(atom.compare("V")==0)   {return  "V    1.402   109.470 3.144   0.016   12.000 2.6789 3 0.000  0.0  0.0000  3.650";}      
    else if(atom.compare("Cr")==0)  {return  "Cr   1.345   90.000  3.023   0.015   12.000 2.4631 6 0.000  0.0  0.0000  3.415";}      
    else if(atom.compare("Mn")==0)  {return  "Mn   1.382   90.000  2.961   0.013   12.000 2.4301 6 0.000  0.0  0.0000  3.325";}      
    else if(atom.compare("Fe")==0)  {return  "Fe   1.270   109.470 2.912   0.013   12.000 2.4301 3 0.000  0.0  0.0000  3.760";}      
    else if(atom.compare("Co")==0)  {return  "Co   1.241   90.000  2.872   0.014   12.000 2.4301 6 0.000  0.0  0.0000  4.105";}      
    else if(atom.compare("Ni")==0)  {return  "Ni   1.164   90.000  2.834   0.015   12.000 2.4301 4 0.000  0.0  0.0000  4.465";}      
    else if(atom.compare("Cu")==0)  {return  "Cu   1.302   109.470 3.495   0.005   12.000 1.7565 3 0.000  0.0  0.0000  3.729";}      
    else if(atom.compare("Zn")==0)  {return  "Zn   1.193   109.470 2.763   0.124   12.000 1.3084 3 0.000  0.0  0.0000  5.106";}      
    else if(atom.compare("Ga")==0)  {return  "Ga   1.260   109.470 4.383   0.415   11.000 1.8206 3 0.000  0.0  0.0000  2.999";}      
    else if(atom.compare("Ge")==0)  {return  "Ge   1.197   109.470 4.280   0.379   12.000 2.7888 3 0.701  0.0  0.0000  4.051";}      
    else if(atom.compare("As")==0)  {return  "As   1.211   92.100  4.230   0.309   13.000 2.8640 3 1.500 22.0 86.9735  5.188";}      
    else if(atom.compare("Se")==0)  {return  "Se   1.190   90.600  4.205   0.291   14.000 2.7645 3 0.335  0.0  0.0000  6.428";}      
    else if(atom.compare("Br")==0)  {return  "Br   1.192   180.000 4.189   0.251   15.000 2.5186 0 0.000  0.0  0.0000  7.790";}      
    else if(atom.compare("Kr")==0)  {return  "Kr   1.147   90.000  4.141   0.220   16.000 0.4520 4 0.000  0.0  0.0000  8.505";}      
    else if(atom.compare("Rb")==0)  {return  "Rb   2.260   180.000 4.114   0.040   12.000 1.5922 0 0.000  0.0  0.0000  2.331";}      
    else if(atom.compare("Sr")==0)  {return  "Sr   2.052   90.000  3.641   0.235   12.000 2.4486 6 0.000  0.0  0.0000  3.024";}      
    else if(atom.compare("Y")==0)   {return  "Y    1.698   109.470 3.345   0.072   12.000 3.2573 3 0.000  0.0  0.0000  3.830";}      
    else if(atom.compare("Zr")==0)  {return  "Zr   1.564   109.470 3.124   0.069   12.000 3.6675 3 0.000  0.0  0.0000  3.400";}      
    else if(atom.compare("Nb")==0)  {return  "Nb   1.473   109.470 3.165   0.059   12.000 3.6179 3 0.000  0.0  0.0000  3.550";}      
    else if(atom.compare("Mo")==0)  {return  "Mo   1.467   90.000  3.052   0.056   12.000 3.4021 6 0.000  0.0  0.0000  3.465";}      
    else if(atom.compare("Tc")==0)  {return  "Mo   1.322   90.000  2.998   0.048   12.000 3.4021 3 0.000  0.0  0.0000  3.465";}      
    else if(atom.compare("Ru")==0)  {return  "Tc   1.478   90.000  2.963   0.056   12.000 3.4021 6 0.000  0.0  0.0000  3.290";}      
    else if(atom.compare("Rh")==0)  {return  "Ru   1.332   90.000  2.929   0.053   12.000 3.4021 6 0.000  0.0  0.0000  3.575";}      
    else if(atom.compare("Pd")==0)  {return  "Pd   1.338   90.000  2.899   0.048   12.000 3.2077 4 0.000  0.0  0.0000  4.320";}      
    else if(atom.compare("Ag")==0)  {return  "Ag   1.386   180.000 3.148   0.036   12.000 1.9557 1 0.200  0.0  0.0000  4.436";}      
    else if(atom.compare("Cd")==0)  {return  "Cd   1.403   109.470 2.848   0.228   12.000 1.6525 3 0.000  0.0  0.0000  5.034";}      
    else if(atom.compare("In")==0)  {return  "In   1.459   109.470 4.463   0.599   11.000 2.0704 3 0.000  0.0  0.0000  2.997";}      
    else if(atom.compare("Sn")==0)  {return  "Sn   1.398   109.470 4.392   0.567   12.000 2.9608 3 0.199  0.0  0.0000  3.987";}      
    else if(atom.compare("Sb")==0)  {return  "Sb   1.407   91.600  4.420   0.449   13.000 2.7042 3 1.100 22.0 87.7047  4.899";}      
    else if(atom.compare("Te")==0)  {return  "Te   1.386   90.250  4.470   0.398   14.000 2.8821 3 0.300  0.0  0.0000  5.816";}     
    else if(atom.compare("I")==0)   {return  "I    1.382   180.000 4.500   0.339   15.000 2.6537 0 0.000  0.0  0.0000  6.822";}     
    else if(atom.compare("Xe")==0)  {return  "Xe   1.267   90.000  4.404   0.332   12.000 0.5560 4 0.000  0.0  0.0000  7.595";}     
    else if(atom.compare("Cs")==0)  {return  "Cs   2.570   180.000 4.517   0.045   12.000 1.5728 0 0.000  0.0  0.0000  2.183";}     
    else if(atom.compare("Ba")==0)  {return  "Ba   2.277   90.000  3.703   0.364   12.000 2.7266 6 0.000  0.0  0.0000  2.814";}     
    else if(atom.compare("La")==0)  {return  "La   1.943   109.470 3.522   0.017   12.000 3.3049 3 0.000  0.0  0.0000  2.8355";}     
    else if(atom.compare("Ce")==0)  {return  "Ce   1.841   90.000  3.556   0.013   12.000 3.3049 6 0.000  0.0  0.0000  2.774";}     
    else if(atom.compare("Pr")==0)  {return  "Pr   1.823   90.000  3.606   0.010   12.000 3.3049 6 0.000  0.0  0.0000  2.858";}     
    else if(atom.compare("Nd")==0)  {return  "Nd   1.816   90.000  3.575   0.010   12.000 3.3049 6 0.000  0.0  0.0000  2.8685";}     
    else if(atom.compare("Pm")==0)  {return  "Pm   1.801   90.000  3.547   0.009   12.000 3.3049 6 0.000  0.0  0.0000  2.881";}     
    else if(atom.compare("Sm")==0)  {return  "Sm   1.780   90.000  3.520   0.008   12.000 3.3049 6 0.000  0.0  0.0000  2.9115";}     
    else if(atom.compare("Eu")==0)  {return  "Eu   1.771   90.000  3.493   0.008   12.000 3.3049 6 0.000  0.0  0.0000  2.8785";}     
    else if(atom.compare("Gd")==0)  {return  "Gd   1.735   90.000  3.368   0.009   12.000 3.3049 6 0.000  0.0  0.0000  3.1665";}     
    else if(atom.compare("Tb")==0)  {return  "Tb   1.732   90.000  3.451   0.007   12.000 3.3049 6 0.000  0.0  0.0000  3.018";}     
    else if(atom.compare("Dy")==0)  {return  "Dy   1.710   90.000  3.428   0.007   12.000 3.3049 6 0.000  0.0  0.0000  3.0555";}     
    else if(atom.compare("Ho")==0)  {return  "Ho   1.696   90.000  3.409   0.007   12.000 3.4157 6 0.000  0.0  0.0000  3.127";}     
    else if(atom.compare("Er")==0)  {return  "Er   1.673   90.000  3.391   0.007   12.000 3.3049 6 0.000  0.0  0.0000  3.1865";}      
    else if(atom.compare("Tm")==0)  {return  "Tm   1.660   90.000  3.374   0.006   12.000 3.3049 6 0.000  0.0  0.0000  3.2514";}      
    else if(atom.compare("Yb")==0)  {return  "Yb   1.637   90.000  3.355   0.228   12.000 2.6177 6 0.000  0.0  0.0000  3.2889";}      
    else if(atom.compare("Lu")==0)  {return  "Lu   1.671   90.000  3.640   0.041   12.000 3.2709 6 0.000  0.0  0.0000  2.9629";}      
    else if(atom.compare("Hf")==0)  {return  "Hf   1.611   109.470 3.141   0.072   12.000 3.9212 3 0.000  0.0  0.0000  3.7000";}      
    else if(atom.compare("Ta")==0)  {return  "Ta   1.511   109.470 3.170   0.081   12.000 4.0748 3 0.000  0.0  0.0000  5.1000";}      
    else if(atom.compare("W")==0)   {return  "W    1.392   90.000  3.069   0.067   12.000 3.6937 6 0.000  0.0  0.0000  4.6300";}      
    else if(atom.compare("Re")==0)  {return  "Re   1.372   90.000  2.954   0.066   12.000 3.6937 6 0.000  0.0  0.0000  3.9600";}      
    else if(atom.compare("Os")==0)  {return  "Os   1.372   90.000  3.120   0.037   12.000 3.6937 6 0.000  0.0  0.0000  5.1400";}      
    else if(atom.compare("Ir")==0)  {return  "Ir   1.371   90.000  2.840   0.073   12.000 3.7307 6 0.000  0.0  0.0000  5.0000";}      
    else if(atom.compare("Pt")==0)  {return  "Pt   1.364   90.000  2.754   0.080   12.000 3.3817 4 0.000  0.0  0.0000  4.7900";}      
    else if(atom.compare("Au")==0)  {return  "Au   1.262   90.000  3.293   0.039   12.000 2.6255 4 0.000  0.0  0.0000  4.8940";}      
    else if(atom.compare("Hg")==0)  {return  "Hg   1.340   180.000 2.705   0.385   12.000 1.7497 1 0.000  0.0  0.0000  6.2700";}      
    else if(atom.compare("Tl")==0)  {return  "Tl   1.518   120.000 4.347   0.680   11.000 2.0685 3 0.000  0.0  0.0000  3.2000";}      
    else if(atom.compare("Pb")==0)  {return  "Pb   1.459   109.470 4.297   0.663   12.000 2.8461 3 0.100  0.0  0.0000  3.9000";}      
    else if(atom.compare("Bi")==0)  {return  "Bi   1.512   90.000  4.370   0.518   13.000 2.4700 3 1.000 22.0 90.0000  4.6900";}      
    else if(atom.compare("Po")==0)  {return  "Po   1.500   90.000  4.709   0.325   14.000 2.3329 3 0.300  0.0  0.0000  4.2100";}      
    else if(atom.compare("At")==0)  {return  "At   1.545   180.000 4.750   0.284   15.000 2.2357 0 0.000  0.0  0.0000  4.7500";}      
    else if(atom.compare("Rn")==0)  {return  "Rn   1.420   90.000  4.765   0.248   16.000 0.5832 4 0.000  0.0  0.0000  5.3700";}      
    else if(atom.compare("Fr")==0)  {return  "Fr   2.880   180.000 4.900   0.050   12.000 1.8469 0 0.000  0.0  0.0000  2.0000";}      
    else if(atom.compare("Ra")==0)  {return  "Ra   2.512   90.000  3.677   0.404   12.000 2.9161 6 0.000  0.0  0.0000  2.8430";}      
    else if(atom.compare("Ac")==0)  {return  "Ac   1.983   90.000  3.478   0.033   12.000 3.8882 6 0.000  0.0  0.0000  2.8350";}      
    else if(atom.compare("Th")==0)  {return  "Th   1.721   90.000  3.396   0.026   12.000 4.2021 6 0.000  0.0  0.0000  3.1750";}      
    else if(atom.compare("Pa")==0)  {return  "Pa   1.711   90.000  3.424   0.022   12.000 3.8882 6 0.000  0.0  0.0000  2.9850";}      
    else if(atom.compare("U")==0)   {return  "U    1.684   90.000  3.395   0.022   12.000 3.8882 6 0.000  0.0  0.0000  3.3410";}      
    else if(atom.compare("Np")==0)  {return  "Np   1.666   90.000  3.424   0.019   12.000 3.8882 6 0.000  0.0  0.0000  3.5490";}      
    else if(atom.compare("Pu")==0)  {return  "Pu   1.657   90.000  3.424   0.016   12.000 3.8882 6 0.000  0.0  0.0000  3.2430";}      
    else if(atom.compare("Am")==0)  {return  "Am   1.660   90.000  3.381   0.014   12.000 3.8882 6 0.000  0.0  0.0000  2.9895";}      
    else if(atom.compare("Cm")==0)  {return  "Cm   1.801   90.000  3.326   0.013   12.000 3.8882 6 0.000  0.0  0.0000  2.8315";}      
    else if(atom.compare("Bk")==0)  {return  "Bk   1.761   90.000  3.339   0.013   12.000 3.8882 6 0.000  0.0  0.0000  3.1935";}      
    else if(atom.compare("Cf")==0)  {return  "Cf   1.750   90.000  3.313   0.013   12.000 3.8882 6 0.000  0.0  0.0000  3.1970";}      
    else if(atom.compare("Es")==0)  {return  "Es   1.724   90.000  3.299   0.012   12.000 3.8882 6 0.000  0.0  0.0000  3.3330";}      
    else if(atom.compare("Fm")==0)  {return  "Fm   1.712   90.000  3.286   0.012   12.000 3.8882 6 0.000  0.0  0.0000  3.4000";}      
    else if(atom.compare("Md")==0)  {return  "Md   1.689   90.000  3.274   0.011   12.000 3.8882 6 0.000  0.0  0.0000  3.4700";}      
    else if(atom.compare("No")==0)  {return  "No   1.679   90.000  3.248   0.011   12.000 3.8882 6 0.000  0.0  0.0000  3.4750";}      
    else if(atom.compare("Lr")==0)  {return  "Lr   1.698   90.000  3.236   0.011   12.000 3.8882 6 0.000  0.0  0.0000  3.5000";}      
    else {return "Unknown";}                                                                                                    
  }                                                                                                                            
} // namespace pocc                                                                                                            
                                                                                                                               
// ***************************************************************************                                                 
// string pocc::ReturnUFFParameters(string atom)                                                                               
// ***************************************************************************                                                 
// UFF parameters                                                                                                              
// Force field parameters for UFF, the Universal Force Field                                                                   
// Atom   r1 theta0  x1  D1  zeta    Z1  Vi  Uj  Xi  Hard    Radius                                                            
namespace pocc {                                                                                                               
  string ReturnUFFParameters(string atom) {                                                                                    
    if(atom=="H")       {return "H   0.354   180.000 2.886   0.044   12.000  0.712   0.000   0.000   4.528   6.945   0.371   ";}   
    else if(atom=="He") {return "He  0.849   90.000  2.362   0.056   15.240  0.098   0.000   0.000   9.660   14.920  1.300   ";}   
    else if(atom=="Li") {return "Li  1.336   180.000 2.451   0.025   12.000  1.026   0.000   2.000   3.006   2.386   1.557   ";}   
    else if(atom=="Be") {return "Be  1.074   109.470 2.745   0.085   12.000  1.565   0.000   2.000   4.877   4.443   1.240   ";}  
    else if(atom=="B")  {return "B   0.838   109.470 4.083   0.180   12.052  1.755   0.000   2.000   5.110   4.750   0.822   ";}  
    else if(atom=="C")  {return "C   0.706   180.000 3.851   0.105   12.730  1.912   0.000   2.000   5.343   5.063   0.759   ";}  
    else if(atom=="N")  {return "N   0.656   180.000 3.660   0.069   13.407  2.544   0.000   2.000   6.899   5.880   0.715   ";}  
    else if(atom=="O")  {return "O   0.639   180.000 3.500   0.060   14.085  2.300   0.000   2.000   8.741   6.682   0.669   ";}  
    else if(atom=="F")  {return "F   0.668   180.000 3.364   0.050   14.762  1.735   0.000   2.000   10.874  7.474   0.706   ";}  
    else if(atom=="Ne") {return "Ne  0.920   90.000  3.243   0.042   15.440  0.194   0.000   2.000   11.040  10.550  1.768   ";}  
    else if(atom=="Na") {return "Na  1.539   180.000 2.983   0.030   12.000  1.081   0.000   1.250   2.843   2.296   2.085   ";}  
    else if(atom=="Mg") {return "Mg  1.421   109.470 3.021   0.111   12.000  1.787   0.000   1.250   3.951   3.693   1.500   ";}  
    else if(atom=="Al") {return "Al  1.244   109.470 4.499   0.505   11.278  1.792   0.000   1.250   4.060   3.590   1.201   ";}  
    else if(atom=="Si") {return "Si  1.117   109.470 4.295   0.402   12.175  2.323   1.225   1.250   4.168   3.487   1.176   ";}  
    else if(atom=="P")  {return "P   1.101   93.800  4.147   0.305   13.072  2.863   2.400   1.250   5.463   4.000   1.102   ";}  
    else if(atom=="S")  {return "S   1.064   92.1000 4.035   0.274   13.969  2.703   0.000   1.250   6.928   4.486   1.047   ";}  
    else if(atom=="Cl") {return "Cl  1.044   180.000 3.947   0.227   14.866  2.348   0.000   1.250   8.564   4.946   0.994   ";}  
    else if(atom=="Ar") {return "Ar  1.032   90.000  3.868   0.185   15.763  0.300   0.000   1.250   9.465   6.355   2.108   ";}  
    else if(atom=="K")  {return "K   1.953   180.000 3.812   0.035   12.000  1.165   0.000   0.700   2.421   1.920   2.586   ";}  
    else if(atom=="Ca") {return "Ca  1.761   90.000  3.399   0.238   12.000  2.141   0.000   0.700   3.231   2.880   2.000   ";}  
    else if(atom=="Sc") {return "Sc  1.513   109.470 3.295   0.019   12.000  2.592   0.000   0.700   3.395   3.080   1.750   ";}  
    else if(atom=="Ti") {return "Ti  1.412   109.470 3.175   0.017   12.000  2.659   0.000   0.700   3.470   3.380   1.607   ";}  
    else if(atom=="V")  {return "V   1.402   109.470 3.144   0.016   12.000  2.679   0.000   0.700   3.650   3.410   1.470   ";}  
    else if(atom=="Cr") {return "Cr  1.345   90.000  3.023   0.015   12.000  2.463   0.000   0.700   3.415   3.865   1.402   ";}  
    else if(atom=="Mn") {return "Mn  1.382   90.000  2.961   0.013   12.000  2.430   0.000   0.700   3.325   4.105   1.533   ";}  
    else if(atom=="Fe") {return "Fe  1.270   109.470 2.912   0.013   12.000  2.430   0.000   0.700   3.760   4.140   1.393   ";}  
    else if(atom=="Co") {return "Co  1.241   90.000  2.872   0.014   12.000  2.430   0.000   0.700   4.105   4.175   1.406   ";}  
    else if(atom=="Ni") {return "Ni  1.164   90.000  2.834   0.015   12.000  2.430   0.000   0.700   4.465   4.205   1.398   ";}  
    else if(atom=="Cu") {return "Cu  1.302   109.470 3.495   0.005   12.000  1.756   0.000   0.700   4.200   4.220   1.434   ";}  
    else if(atom=="Zn") {return "Zn  1.193   109.470 2.763   0.124   12.000  1.308   0.000   0.700   5.106   4.285   1.400   ";}  
    else if(atom=="Ga") {return "Ga  1.260   109.470 4.383   0.415   11.000  1.821   0.000   0.700   3.641   3.160   1.211   ";}  
    else if(atom=="Ge") {return "Ge  1.197   109.470 4.280   0.379   12.000  2.789   0.701   0.700   4.051   3.438   1.189   ";}  
    else if(atom=="As") {return "As  1.211   92.100  4.230   0.309   13.000  2.864   1.500   0.700   5.188   3.809   1.204   ";}  
    else if(atom=="Se") {return "Se  1.190   90.600  4.205   0.291   14.000  2.764   0.335   0.700   6.428   4.131   1.224   ";}  
    else if(atom=="Br") {return "Br  1.192   180.000 4.189   0.251   15.000  2.519   0.000   0.700   7.790   4.425   1.141   ";}  
    else if(atom=="Kr") {return "Kr  1.147   90.000  4.141   0.220   16.000  0.452   0.000   0.700   8.505   5.715   2.270   ";}  
    else if(atom=="Rb") {return "Rb  2.260   180.000 4.114   0.040   12.000  1.592   0.000   0.200   2.331   1.846   2.770   ";}  
    else if(atom=="Sr") {return "Sr  2.052   90.000  3.641   0.235   12.000  2.449   0.000   0.200   3.024   2.440   2.415   ";}  
    else if(atom=="Y")  {return "Y   1.698   109.470 3.345   0.072   12.000  3.257   0.000   0.200   3.830   2.810   1.998   ";}  
    else if(atom=="Zr") {return "Zr  1.564   109.470 3.124   0.069   12.000  3.667   0.000   0.200   3.400   3.550   1.758   ";}  
    else if(atom=="Nb") {return "Nb  1.473   109.470 3.165   0.059   12.000  3.618   0.000   0.200   3.550   3.380   1.603   ";}  
    else if(atom=="Mo") {return "Mo  1.467   90.000  3.052   0.056   12.000  3.400   0.000   0.200   3.465   3.755   1.530   ";}  
    else if(atom=="Tc") {return "Tc  1.322   90.000  2.998   0.048   12.000  3.400   0.000   0.200   3.290   3.990   1.500   ";}  
    else if(atom=="Ru") {return "Ru  1.478   90.000  2.963   0.056   12.000  3.400   0.000   0.200   3.575   4.015   1.500   ";}  
    else if(atom=="Rh") {return "Rh  1.332   90.000  2.929   0.053   12.000  3.500   0.000   0.200   3.975   4.005   1.509   ";}  
    else if(atom=="Pd") {return "Pd  1.338   90.000  2.899   0.048   12.000  3.210   0.000   0.200   4.320   4.000   1.544   ";}  
    else if(atom=="Ag") {return "Ag  1.386   180.000 3.148   0.036   12.000  1.956   0.000   0.200   4.436   3.134   1.622   ";}  
    else if(atom=="Cd") {return "Cd  1.403   109.470 2.848   0.228   12.000  1.650   0.000   0.200   5.034   3.957   1.600   ";}  
    else if(atom=="In") {return "In  1.459   109.470 4.463   0.599   11.000  2.070   0.000   0.200   3.506   2.896   1.404   ";}  
    else if(atom=="Sn") {return "Sn  1.398   109.470 4.392   0.567   12.000  2.961   0.199   0.200   3.987   3.124   1.354   ";}  
    else if(atom=="Sb") {return "Sb  1.407   91.600  4.420   0.449   13.000  2.704   1.100   0.200   4.899   3.342   1.404   ";}  
    else if(atom=="Te") {return "Te  1.386   90.250  4.470   0.398   14.000  2.882   0.300   0.200   5.816   3.526   1.380   ";}  
    else if(atom=="I")  {return "I   1.382   180.000 4.500   0.339   15.000  2.650   0.000   0.200   6.822   3.762   1.333   ";}  
    else if(atom=="Xe") {return "Xe  1.267   90.000  4.404   0.332   12.000  0.556   0.000   0.200   7.595   4.975   2.459   ";}  
    else if(atom=="Cs") {return "Cs  2.570   180.000 4.517   0.045   12.000  1.573   0.000   0.100   2.183   1.711   2.984   ";}  
    else if(atom=="Ba") {return "Ba  2.277   90.000  3.703   0.364   12.000  2.727   0.000   0.100   2.814   2.396   2.442   ";}  
    else if(atom=="La") {return "La  1.943   109.470 3.522   0.017   12.000  3.300   0.000   0.100   2.836   2.742   2.071   ";}  
    else if(atom=="Ce") {return "Ce  1.841   90.000  3.556   0.013   12.000  3.300   0.000   0.100   2.774   2.692   1.925   ";}  
    else if(atom=="Pr") {return "Pr  1.823   90.000  3.606   0.010   12.000  3.300   0.000   0.100   2.858   2.564   2.007   ";}  
    else if(atom=="Nd") {return "Nd  1.816   90.000  3.575   0.010   12.000  3.300   0.000   0.100   2.869   2.621   2.007   ";}  
    else if(atom=="Pm") {return "Pm  1.801   90.000  3.547   0.009   12.000  3.300   0.000   0.100   2.881   2.673   2.000   ";}  
    else if(atom=="Sm") {return "Sm  1.780   90.000  3.520   0.008   12.000  3.300   0.000   0.100   2.912   2.720   1.978   ";}  
    else if(atom=="Eu") {return "Eu  1.771   90.000  3.493   0.008   12.000  3.300   0.000   0.100   2.879   2.788   2.227   ";}  
    else if(atom=="Gd") {return "Gd  1.735   90.000  3.368   0.009   12.000  3.300   0.000   0.100   3.167   2.975   1.968   ";}  
    else if(atom=="Tb") {return "Tb  1.732   90.000  3.451   0.007   12.000  3.300   0.000   0.100   3.018   2.834   1.954   ";}  
    else if(atom=="Dy") {return "Dy  1.710   90.000  3.428   0.007   12.000  3.300   0.000   0.100   3.056   2.872   1.934   ";}  
    else if(atom=="Ho") {return "Ho  1.696   90.000  3.409   0.007   12.000  3.416   0.000   0.100   3.127   2.891   1.925   ";}  
    else if(atom=="Er") {return "Er  1.673   90.000  3.391   0.007   12.000  3.300   0.000   0.100   3.187   2.915   1.915   ";}  
    else if(atom=="Tm") {return "Tm  1.660   90.000  3.374   0.006   12.000  3.300   0.000   0.100   3.251   2.933   2.000   ";}  
    else if(atom=="Yb") {return "Yb  1.637   90.000  3.355   0.228   12.000  2.618   0.000   0.100   3.289   2.965   2.158   ";}  
    else if(atom=="Lu") {return "Lu  1.671   90.000  3.640   0.041   12.000  3.271   0.000   0.100   2.963   2.463   1.896   ";}  
    else if(atom=="Hf") {return "Hf  1.611   109.470 3.141   0.072   12.000  3.921   0.000   0.100   3.700   3.400   1.759   ";}  
    else if(atom=="Ta") {return "Ta  1.511   109.470 3.170   0.081   12.000  4.075   0.000   0.100   5.100   2.850   1.605   ";}  
    else if(atom=="W")  {return "W   1.392   90.000  3.069   0.067   12.000  3.700   0.000   0.100   4.630   3.310   1.538   ";}  
    else if(atom=="Re") {return "Re  1.372   90.000  2.954   0.066   12.000  3.700   0.000   0.100   3.960   3.920   1.600   ";}  
    else if(atom=="Os") {return "Os  1.372   90.000  3.120   0.037   12.000  3.700   0.000   0.100   5.140   3.630   1.700   ";}  
    else if(atom=="Ir") {return "Ir  1.371   90.000  2.840   0.073   12.000  3.731   0.000   0.100   5.000   4.000   1.866   ";}  
    else if(atom=="Pt") {return "Pt  1.364   90.000  2.754   0.080   12.000  3.382   0.000   0.100   4.790   4.430   1.557   ";}  
    else if(atom=="Au") {return "Au  1.262   90.000  3.293   0.039   12.000  2.625   0.000   0.100   4.894   2.586   1.618   ";}  
    else if(atom=="Hg") {return "Hg  1.340   180.000 2.705   0.385   12.000  1.750   0.000   0.100   6.270   4.160   1.600   ";}  
    else if(atom=="Tl") {return "Tl  1.518   120.000 4.347   0.680   11.000  2.068   0.000   0.100   3.200   2.900   1.530   ";}  
    else if(atom=="Pb") {return "Pb  1.459   109.470 4.297   0.663   12.000  2.846   0.100   0.100   3.900   3.530   1.444   ";}  
    else if(atom=="Bi") {return "Bi  1.512   90.000  4.370   0.518   13.000  2.470   1.000   0.100   4.690   3.740   1.514   ";}  
    else if(atom=="Po") {return "Po  1.500   90.000  4.709   0.325   14.000  2.330   0.300   0.100   4.210   4.210   1.480   ";}  
    else if(atom=="At") {return "At  1.545   180.000 4.750   0.284   15.000  2.240   0.000   0.100   4.750   4.750   1.470   ";}  
    else if(atom=="Rn") {return "Rn  1.420   90.000  4.765   0.248   16.000  0.583   0.000   0.100   5.370   5.370   2.200   ";}  
    else if(atom=="Fr") {return "Fr  2.880   180.000 4.900   0.050   12.000  1.847   0.000   0.000   2.000   2.000   2.300   ";}  
    else if(atom=="Ra") {return "Ra  2.512   90.000  3.677   0.404   12.000  2.920   0.000   0.000   2.843   2.434   2.200   ";}  
    else if(atom=="Ac") {return "Ac  1.983   90.000  3.478   0.033   12.000  3.900   0.000   0.000   2.835   2.835   2.108   ";}  
    else if(atom=="Th") {return "Th  1.721   90.000  3.396   0.026   12.000  4.202   0.000   0.000   3.175   2.905   2.018   ";}  
    else if(atom=="Pa") {return "Pa  1.711   90.000  3.424   0.022   12.000  3.900   0.000   0.000   2.985   2.905   1.800   ";}  
    else if(atom=="U")  {return "U   1.684   90.000  3.395   0.022   12.000  3.900   0.000   0.000   3.341   2.853   1.713   ";}  
    else if(atom=="Np") {return "Np  1.666   90.000  3.424   0.019   12.000  3.900   0.000   0.000   3.549   2.717   1.800   ";}  
    else if(atom=="Pu") {return "Pu  1.657   90.000  3.424   0.016   12.000  3.900   0.000   0.000   3.243   2.819   1.840   ";}  
    else if(atom=="Am") {return "Am  1.660   90.000  3.381   0.014   12.000  3.900   0.000   0.000   2.990   3.004   1.942   ";}  
    else if(atom=="Cm") {return "Cm  1.801   90.000  3.326   0.013   12.000  3.900   0.000   0.000   2.832   3.190   1.900   ";}  
    else if(atom=="Bk") {return "Bk  1.761   90.000  3.339   0.013   12.000  3.900   0.000   0.000   3.194   3.036   1.900   ";}  
    else if(atom=="Cf") {return "Cf  1.750   90.000  3.313   0.013   12.000  3.900   0.000   0.000   3.197   3.101   1.900   ";}  
    else if(atom=="Es") {return "Es  1.724   90.000  3.299   0.012   12.000  3.900   0.000   0.000   3.333   3.089   1.900   ";}  
    else if(atom=="Fm") {return "Fm  1.712   90.000  3.286   0.012   12.000  3.900   0.000   0.000   3.400   3.100   1.900   ";}  
    else if(atom=="Md") {return "Md  1.689   90.000  3.274   0.011   12.000  3.900   0.000   0.000   3.470   3.110   1.900   ";}  
    else if(atom=="No") {return "No  1.679   90.000  3.248   0.011   12.000  3.900   0.000   0.000   3.475   3.175   1.900   ";}  
    else if(atom=="Lw") {return "Lw  1.698   90.000  3.236   0.011   12.000  3.900   0.000   0.000   3.500   3.200   1.900   ";}  
    else {return "Unknown";}                                                                                                      
  }                                                                                                                             
} // namespace pocc                                                                                                             
                                                                                                                                
// ***********************************************************************************************************************************************************
// ***********************************************************************************************************************************************************
// constructor for Class UFFPara                                                                                                
namespace pocc {                                                                                                                
  UFFPara::UFFPara() {
    symbol = "";
    r1 = 0;
    theta0 = 0;
    x1 = 0;
    D1 = 0;
    zeta = 0;
    Z1 = 0;
    Vi = 0;
    Uj = 0;
    Xi = 0;
    hard = 0;
    radius = 0;
  }
} // namespace pocc 

// destructor
namespace pocc {
  UFFPara::~UFFPara() {
    Free();
  }} // namespace pocc 


namespace pocc {
  void UFFPara::Free() {}
} // namespace pocc 

namespace pocc {
  void UFFPara::GetUFFParameters(string atom) {
    //string strData=ReturnUFFParameters(atom);
    string strData=ReturnUFFParameters(KBIN::VASP_PseudoPotential_CleanName(atom));
    vector<string> UFFData;
    aurostd::string2tokens(strData, UFFData, " ");
    symbol = UFFData.at(0); //Symbol
    r1 = aurostd::string2utype<double>(UFFData.at(1)); 
    theta0 = aurostd::string2utype<double>(UFFData.at(2)); 
    x1 = aurostd::string2utype<double>(UFFData.at(3));  //nonbond distance
    D1 = aurostd::string2utype<double>(UFFData.at(4));  //nonbond energy
    zeta = aurostd::string2utype<double>(UFFData.at(5));  //scale
    Z1 = aurostd::string2utype<double>(UFFData.at(6)); 
    Vi = aurostd::string2utype<double>(UFFData.at(7)); 
    Uj = aurostd::string2utype<double>(UFFData.at(8)); 
    Xi = aurostd::string2utype<double>(UFFData.at(9)); //electronegativity
    hard = aurostd::string2utype<double>(UFFData.at(10)); 
    radius = aurostd::string2utype<double>(UFFData.at(11)); 
  }
} // namespace pocc 

// ***********************************************************************************************************************************************************
// ***********************************************************************************************************************************************************

// ***************************************************************************
// string pocc::ReturnAtomProperties(string& atom)
// ***************************************************************************
// http://en.wikipedia.org/wiki/Atomic_radii_of_the_elements_(data_page)                
// name; symbol; atomic number; atomic mass; atomic radius; pauling electronegativity
namespace pocc {
  string ReturnAtomProperties(string atom) {
    if(atom.compare("H")==0)        {return  "Hydrogen         H     1   1.007940    0.53    2.20     ";}  
    else if(atom.compare("He")==0)  {return  "Helium           He    2   4.002602    0.31    4.16     ";}  
    else if(atom.compare("Li")==0)  {return  "Lithium          Li    3   6.941000    1.67    0.98     ";}  
    else if(atom.compare("Be")==0)  {return  "Beryllium        Be    4   9.012182    1.12    1.57     ";}  
    else if(atom.compare("B")==0)   {return  "Boron            B     5   10.811000   0.87    2.04     ";}  
    else if(atom.compare("C")==0)   {return  "Carbon           C     6   12.010700   0.67    2.55     ";}  
    else if(atom.compare("N")==0)   {return  "Nitrogen         N     7   14.006700   0.56    3.04     ";}  
    else if(atom.compare("O")==0)   {return  "Oxygen           O     8   15.999400   0.48    3.44     ";}  
    else if(atom.compare("F")==0)   {return  "Fluorine         F     9   18.998403   0.42    3.98     ";}  
    else if(atom.compare("Ne")==0)  {return  "Neon             Ne    10  20.179700   0.38    4.79     ";}  
    else if(atom.compare("Na")==0)  {return  "Sodium           Na    11  22.989770   1.90    0.93     ";}  
    else if(atom.compare("Mg")==0)  {return  "Magnesium        Mg    12  24.305000   1.45    1.31     ";}  
    else if(atom.compare("Al")==0)  {return  "Aluminium        Al    13  26.981538   1.18    1.61     ";}  
    else if(atom.compare("Si")==0)  {return  "Silicon          Si    14  28.085500   1.11    1.90     ";}  
    else if(atom.compare("P")==0)   {return  "Phosphorus       P     15  30.973761   0.98    2.19     ";}  
    else if(atom.compare("S")==0)   {return  "Sulphur          S     16  32.065000   0.88    2.58     ";}  
    else if(atom.compare("Cl")==0)  {return  "Chlorine         Cl    17  35.453000   0.79    3.16     ";}  
    else if(atom.compare("Ar")==0)  {return  "Argon            Ar    18  39.948000   0.71    3.24     ";}  
    else if(atom.compare("K")==0)   {return  "Potassium        K     19  39.098300   2.43    0.82     ";}  
    else if(atom.compare("Ca")==0)  {return  "Calcium          Ca    20  40.078000   1.94    1.00     ";}  
    else if(atom.compare("Sc")==0)  {return  "Scandium         Sc    21  44.955910   1.84    1.36     ";}  
    else if(atom.compare("Ti")==0)  {return  "Titanium         Ti    22  47.867000   1.76    1.54     ";}  
    else if(atom.compare("V")==0)   {return  "Vanadium         V     23  50.941500   1.71    1.63     ";}  
    else if(atom.compare("Cr")==0)  {return  "Chromium         Cr    24  51.996100   1.66    1.66     ";}  
    else if(atom.compare("Mn")==0)  {return  "Manganese        Mn    25  54.938049   1.61    1.55     ";}  
    else if(atom.compare("Fe")==0)  {return  "Iron             Fe    26  55.845000   1.56    1.83     ";}  
    else if(atom.compare("Co")==0)  {return  "Cobalt           Co    27  58.933200   1.52    1.88     ";}  
    else if(atom.compare("Ni")==0)  {return  "Nickel           Ni    28  58.693400   1.49    1.91     ";}  
    else if(atom.compare("Cu")==0)  {return  "Copper           Cu    29  63.546000   1.45    1.90     ";}  
    else if(atom.compare("Zn")==0)  {return  "Zinc             Zn    30  65.390000   1.42    1.65     ";}  
    else if(atom.compare("Ga")==0)  {return  "Gallium          Ga    31  69.723000   1.36    1.81     ";}  
    else if(atom.compare("Ge")==0)  {return  "Germanium        Ge    32  72.640000   1.25    2.01     ";}  
    else if(atom.compare("As")==0)  {return  "Arsenic          As    33  74.921600   1.14    2.18     ";}  
    else if(atom.compare("Se")==0)  {return  "Selenium         Se    34  78.960000   1.03    2.55     ";}  
    else if(atom.compare("Br")==0)  {return  "Bromine          Br    35  79.904000   0.94    2.96     ";}  
    else if(atom.compare("Kr")==0)  {return  "Krypton          Kr    36  83.800000   0.88    3.00     ";}  
    else if(atom.compare("Rb")==0)  {return  "Rubidium         Rb    37  85.467800   2.65    0.82     ";}  
    else if(atom.compare("Sr")==0)  {return  "Strontium        Sr    38  87.620000   2.19    0.95     ";}  
    else if(atom.compare("Y")==0)   {return  "Yttrium          Y     39  88.905850   2.12    1.22     ";}  
    else if(atom.compare("Zr")==0)  {return  "Zirconium        Zr    40  91.224000   2.06    1.33     ";}  
    else if(atom.compare("Nb")==0)  {return  "Niobium          Nb    41  92.906380   1.98    1.60     ";}  
    else if(atom.compare("Mo")==0)  {return  "Molybdenum       Mo    42  95.940000   1.90    2.16     ";}  
    else if(atom.compare("Tc")==0)  {return  "Technetium       Mo    43  98.000000   1.83    1.90     ";}  
    else if(atom.compare("Ru")==0)  {return  "Ruthenium        Tc    44  101.070000  1.78    2.20     ";}  
    else if(atom.compare("Rh")==0)  {return  "Rhodium          Ru    45  102.905500  1.73    2.28     ";}  
    else if(atom.compare("Pd")==0)  {return  "Palladium        Pd    46  106.420000  1.69    2.20     ";}  
    else if(atom.compare("Ag")==0)  {return  "Silver           Ag    47  107.868200  1.65    1.93     ";}  
    else if(atom.compare("Cd")==0)  {return  "Cadmium          Cd    48  112.411000  1.61    1.69     ";}  
    else if(atom.compare("In")==0)  {return  "Indium           In    49  114.818000  1.56    1.78     ";}  
    else if(atom.compare("Sn")==0)  {return  "Tin              Sn    50  118.710000  1.45    1.96     ";}  
    else if(atom.compare("Sb")==0)  {return  "Antimony         Sb    51  121.760000  1.33    2.05     ";}  
    else if(atom.compare("Te")==0)  {return  "Tellurium        Te    52  127.600000  1.23    2.10     ";}  
    else if(atom.compare("I")==0)   {return  "Iodine           I     53  126.904470  1.15    2.66     ";}  
    else if(atom.compare("Xe")==0)  {return  "Xenon            Xe    54  131.293000  1.08    2.60     ";}  
    else if(atom.compare("Cs")==0)  {return  "Cesium           Cs    55  132.905450  2.98    0.79     ";}  
    else if(atom.compare("Ba")==0)  {return  "Barium           Ba    56  137.327000  2.53    0.89     ";}  
    else if(atom.compare("La")==0)  {return  "Lanthanium       La    57  138.905500  1.95    1.10     ";}   
    else if(atom.compare("Ce")==0)  {return  "Cerium           Ce    58  140.116000  1.85    1.12     ";}  
    else if(atom.compare("Pr")==0)  {return  "Praseodymium     Pr    59  140.907650  2.47    1.13     ";}  
    else if(atom.compare("Nd")==0)  {return  "Neodymium        Nd    60  144.240000  2.06    1.14     ";}   
    else if(atom.compare("Pm")==0)  {return  "Promethium       Pm    61  145.000000  2.05    1.07     ";}  
    else if(atom.compare("Sm")==0)  {return  "Samarium         Sm    62  150.360000  2.38    1.17     ";}   
    else if(atom.compare("Eu")==0)  {return  "Europium         Eu    63  151.964000  2.31    1.01     ";}    
    else if(atom.compare("Gd")==0)  {return  "Gadolinium       Gd    64  157.250000  2.33    1.20     ";}   
    else if(atom.compare("Tb")==0)  {return  "Terbium          Tb    65  158.925340  2.25    1.10     ";}  
    else if(atom.compare("Dy")==0)  {return  "Dysprosium       Dy    66  162.500000  2.28    1.22     ";}   
    else if(atom.compare("Ho")==0)  {return  "Holmium          Ho    67  164.930320  2.26    1.23     ";}  
    else if(atom.compare("Er")==0)  {return  "Erbium           Er    68  167.259000  2.26    1.24     ";}   
    else if(atom.compare("Tm")==0)  {return  "Thulium          Tm    69  168.934210  2.22    1.25     ";}   
    else if(atom.compare("Yb")==0)  {return  "Ytterbium        Yb    70  173.040000  2.22    1.06     ";}   
    else if(atom.compare("Lu")==0)  {return  "Lutetium         Lu    71  174.967000  2.17    1.27     ";}   
    else if(atom.compare("Hf")==0)  {return  "Hafnium          Hf    72  178.490000  2.08    1.30     ";}   
    else if(atom.compare("Ta")==0)  {return  "Tantalum         Ta    73  180.947900  2.00    1.50     ";}   
    else if(atom.compare("W")==0)   {return  "Tungsten         W     74  183.840000  1.93    2.36     ";}   
    else if(atom.compare("Re")==0)  {return  "Rhenium          Re    75  186.207000  1.88    1.90     ";}   
    else if(atom.compare("Os")==0)  {return  "Osmium           Os    76  190.230000  1.85    2.20     ";}   
    else if(atom.compare("Ir")==0)  {return  "Iridium          Ir    77  192.217000  1.80    2.20     ";}   
    else if(atom.compare("Pt")==0)  {return  "Platinum         Pt    78  195.078000  1.77    2.28     ";}   
    else if(atom.compare("Au")==0)  {return  "Gold             Au    79  196.966550  1.74    2.54     ";}   
    else if(atom.compare("Hg")==0)  {return  "Mercury          Hg    80  200.590000  1.71    2.00     ";}   
    else if(atom.compare("Tl")==0)  {return  "Thallium         Tl    81  204.383300  1.56    1.62     ";}   
    else if(atom.compare("Pb")==0)  {return  "Lead             Pb    82  207.200000  1.54    2.33     ";} 
    else if(atom.compare("Bi")==0)  {return  "Bismuth          Bi    83  208.980380  1.43    2.02     ";}   
    else if(atom.compare("Po")==0)  {return  "Polonium         Po    84  209.000000  1.35    2.00     ";}   
    else if(atom.compare("At")==0)  {return  "Astatine         At    85  210.000000  1.27    2.20     ";}   
    else if(atom.compare("Rn")==0)  {return  "Radon            Rn    86  222.000000  1.20    2.59     ";}   
    else if(atom.compare("Fr")==0)  {return  "Francium         Fr    87  223.000000  2.18    0.70     ";}   
    else if(atom.compare("Ra")==0)  {return  "Radium           Ra    88  226.000000  2.40    0.90     ";}   
    else if(atom.compare("Ac")==0)  {return  "Actinium         Ac    89  227.000000  2.20    1.10     ";}   
    else if(atom.compare("Th")==0)  {return  "Thorium          Th    90  232.038100  1.65    1.30     ";}   
    else if(atom.compare("Pa")==0)  {return  "Protoactinium    Pa    91  231.035880  1.65    1.50     ";}   
    else if(atom.compare("U")==0)   {return  "Uranium          U     92  238.028910  1.42    1.38     ";}   
    //else {return "Hydrogen         H     1   1.007940   0.53    2.20     ";}                                        
    else {return "Unknown";}                                        
  }                                 
                                    
  // ***********************************************************************************************************************************************************
  // Class for Atom                 
  Atom::Atom() {                     
    name = "";                      
    symbol = "";                    
    number =0;                      
    mass = 0 ;                      
    radius = 0;                     
    Xi = 0;
  }

  Atom::~Atom() {
    Free(); 
  }

  void Atom::Free() {}

  void Atom::GetAtomicProperties(string atom) {
    //string strData=pocc::ReturnAtomProperties(atom);
    string strData=pocc::ReturnAtomProperties(KBIN::VASP_PseudoPotential_CleanName(atom));
    vector<string> atomData;
    aurostd::string2tokens(strData, atomData, " ");
    name = atomData.at(0);
    symbol = atomData.at(1);
    number = aurostd::string2utype<int>(atomData.at(2));
    mass = aurostd::string2utype<double>(atomData.at(3));
    radius = aurostd::string2utype<double>(atomData.at(4));
    Xi = aurostd::string2utype<double>(atomData.at(5)); //pauling electronegativity, Greek chi
  }
} // namespace pocc 

// ***********************************************************************************************************************************************************

#endif  // _AFLOW_POCCUPATION_PARAMS_CPP_

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *                                                                         *
// ***************************************************************************

