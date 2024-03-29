#ifndef _AFLOW_GNUPLOT_FUNCS_CPP_
#define _AFLOW_GNUPLOT_FUNCS_CPP_
// aflow_gnuplot_funcs.cpp automatic generated from GNUPLOT/*
#include <sstream>
// aflow_gnuplot_funcs.cpp automatic generated from GNUPLOT/plotbz.sh
//std::string GNUPLOT_FUNCS_plotbz(std::string DIRECTORY,std::string OPTION1){
std::string GNUPLOT_FUNCS_plotbz(void){
   std::stringstream strstream; strstream << "\
#!/bin/ksh" << std::endl;
   strstream << "# A script to plot brillouin zone and kpath." << std::endl;
   strstream << "# It needs aflow with --bzplotdata option" << std::endl;
   strstream << "# input files: POSCAR.bands" << std::endl;
   strstream << "#" << std::endl;
   strstream << "# written: 2010 wahyu@alumni.duke.edu" << std::endl;
   strstream << "" << std::endl;
   strstream << "cdir=./" << std::endl;
   strstream << "echo plotbz.sh::DIR=  $cdir" << std::endl;
   strstream << "kfile=${cdir}/KPOINTS.bands" << std::endl;
   strstream << "posfile=${cdir}/POSCAR.bands" << std::endl;
   strstream << "datfile=${cdir}/plotbz.dat" << std::endl;
   strstream << "recipfile=${cdir}/recip.dat #reciprocal vectors" << std::endl;
   strstream << "kptsfile=${cdir}/kpts.dat #kpts" << std::endl;
   strstream << "kpathfile=${cdir}/kpath.dat #kpath" << std::endl;
   strstream << "nonkpathfile=${cdir}/nonkpath.dat #irbz non kpath" << std::endl;
   strstream << "fbzfile=${cdir}/fbz.dat #front bz" << std::endl;
   strstream << "bbzfile=${cdir}/bbz.dat #back bz" << std::endl;
   strstream << "" << std::endl;
   strstream << "aflow --bzplotdatauseKPOINTS=$kfile < $posfile > $datfile" << std::endl;
   strstream << "" << std::endl;
   strstream << "#-----reading datfile-----" << std::endl;
   strstream << "#---reciprocal lattice vectors b1 b2 b3---" << std::endl;
   strstream << "exec 0< $datfile" << std::endl;
   strstream << "$IFS read glatt" << std::endl;
   strstream << "$IFS read gpath" << std::endl;
   strstream << "$IFS read xrotview zrotview srest" << std::endl;
   strstream << "read srest" << std::endl;
   strstream << "$IFS read Lb1 Lb2 Lb3 srest" << std::endl;
   strstream << "$IFS read Tb1 Tb2 Tb3 srest" << std::endl;
   strstream << "read srest" << std::endl;
   strstream << "$IFS read b11 b12 b13" << std::endl;
   strstream << "$IFS read b21 b22 b23" << std::endl;
   strstream << "$IFS read b31 b32 b33" << std::endl;
   strstream << "((hb11=b11/2.0)); ((hb12=b12/2.0)); ((hb13=b13/2.0)); " << std::endl;
   strstream << "((hb21=b21/2.0)); ((hb22=b22/2.0)); ((hb23=b23/2.0)); " << std::endl;
   strstream << "((hb31=b31/2.0)); ((hb32=b32/2.0)); ((hb33=b33/2.0)); " << std::endl;
   strstream << "#--length of b1 b2 b3 arrows--" << std::endl;
   strstream << "((db11=b11*Lb1)); ((db12=b12*Lb1)); ((db13=b13*Lb1)); " << std::endl;
   strstream << "((db21=b21*Lb2)); ((db22=b22*Lb2)); ((db23=b23*Lb2)); " << std::endl;
   strstream << "((db31=b31*Lb3)); ((db32=b32*Lb3)); ((db33=b33*Lb3)); " << std::endl;
   strstream << "echo $hb11 $hb12 $hb13 $db11 $db12 $db13 > $recipfile" << std::endl;
   strstream << "echo $hb21 $hb22 $hb23 $db21 $db22 $db23 >> $recipfile" << std::endl;
   strstream << "echo $hb31 $hb32 $hb33 $db31 $db32 $db33 >> $recipfile" << std::endl;
   strstream << "((tb11=hb11+b11*Tb1)); ((tb12=hb12+b12*Tb1)); ((tb13=hb13+b13*Tb1));" << std::endl;
   strstream << "((tb21=hb21+b21*Tb2)); ((tb22=hb22+b22*Tb2)); ((tb23=hb23+b23*Tb2));" << std::endl;
   strstream << "((tb31=hb31+b31*Tb3)); ((tb32=hb32+b32*Tb3)); ((tb33=hb33+b33*Tb3));" << std::endl;
   strstream << "" << std::endl;
   strstream << "#---kpts---" << std::endl;
   strstream << "$IFS read Nkpts srest" << std::endl;
   strstream << "rm -f $kptsfile" << std::endl;
   strstream << "for ((i=0;i<Nkpts;i++))" << std::endl;
   strstream << "do" << std::endl;
   strstream << "read kptsx[$i] kptsy[$i] kptsz[$i] klab[$i]" << std::endl;
   strstream << "echo ${kptsx[$i]} ${kptsy[$i]} ${kptsz[$i]} ${klab[$i]} >> $kptsfile" << std::endl;
   strstream << "done" << std::endl;
   strstream << "" << std::endl;
   strstream << "#---kpath---" << std::endl;
   strstream << "read srest" << std::endl;
   strstream << "read srest" << std::endl;
   strstream << "$IFS read N srest" << std::endl;
   strstream << "rm -f $kpathfile" << std::endl;
   strstream << "for ((i=0;i<N;i++))" << std::endl;
   strstream << "do" << std::endl;
   strstream << "read srest" << std::endl;
   strstream << "echo $srest >> $kpathfile" << std::endl;
   strstream << "done" << std::endl;
   strstream << "#---irbz non kpath---" << std::endl;
   strstream << "$IFS read Nnonkpath srest" << std::endl;
   strstream << "rm -f $nonkpathfile" << std::endl;
   strstream << "for ((i=0;i<Nnonkpath;i++))" << std::endl;
   strstream << "do" << std::endl;
   strstream << "read srest" << std::endl;
   strstream << "echo $srest >> $nonkpathfile" << std::endl;
   strstream << "done" << std::endl;
   strstream << "#---front bz---" << std::endl;
   strstream << "$IFS read N srest" << std::endl;
   strstream << "rm -f $fbzfile" << std::endl;
   strstream << "for ((i=0;i<N;i++))" << std::endl;
   strstream << "do" << std::endl;
   strstream << "read srest" << std::endl;
   strstream << "echo $srest >> $fbzfile" << std::endl;
   strstream << "done" << std::endl;
   strstream << "#---back bz---" << std::endl;
   strstream << "$IFS read N srest" << std::endl;
   strstream << "rm -f $bbzfile" << std::endl;
   strstream << "for ((i=0;i<N;i++))" << std::endl;
   strstream << "do" << std::endl;
   strstream << "read srest" << std::endl;
   strstream << "echo $srest >> $bbzfile" << std::endl;
   strstream << "done" << std::endl;
   strstream << "" << std::endl;
   strstream << "#getting the structure name and icsd number" << std::endl;
   strstream << "name=$(pwd)" << std::endl;
   strstream << "echo $name | sed \"s/\\//\\n/g\" | grep '_ICSD_' > wahyutmp" << std::endl;
   strstream << "exec 0< wahyutmp" << std::endl;
   strstream << "$IFS read name" << std::endl;
   strstream << "rm -f wahyutmp" << std::endl;
   strstream << "" << std::endl;
   strstream << "echo \"set terminal postscript eps color enhanced \\\"Helvetica\\\" 18 \"> plotbz.gnu" << std::endl;
   strstream << "echo \"set output \\\"bz_$name.eps\\\" \">> plotbz.gnu" << std::endl;
   strstream << "echo \"set key off \">> plotbz.gnu" << std::endl;
   strstream << "echo \"unset border \">> plotbz.gnu" << std::endl;
   strstream << "echo \"unset xtics \">> plotbz.gnu" << std::endl;
   strstream << "echo \"unset ytics \">> plotbz.gnu" << std::endl;
   strstream << "echo \"unset ztics \">> plotbz.gnu" << std::endl;
   strstream << "echo \"set view equal xyz \">> plotbz.gnu" << std::endl;
   strstream << "echo \"set view $xrotview,$zrotview \">> plotbz.gnu" << std::endl;
   strstream << "echo \"set title \\\"$glatt  path: $gpath\\\" \">> plotbz.gnu" << std::endl;
   strstream << "for ((i=0;i<Nkpts;i++))" << std::endl;
   strstream << "do" << std::endl;
   strstream << "echo \"set label \\\"${klab[$i]}\\\" at ${kptsx[i]},${kptsy[i]},${kptsz[i]}\" >> plotbz.gnu" << std::endl;
   strstream << "done" << std::endl;
   strstream << "echo \"set label \\\"b{/Symbol _1}\\\" at $tb11,$tb12,$tb13 font \\\"Helvetica,14\\\" \">> plotbz.gnu" << std::endl;
   strstream << "echo \"set label \\\"b{/Symbol _2}\\\" at $tb21,$tb22,$tb23 font \\\"Helvetica,14\\\" \">> plotbz.gnu" << std::endl;
   strstream << "echo \"set label \\\"b{/Symbol _3}\\\" at $tb31,$tb32,$tb33 font \\\"Helvetica,14\\\" \">> plotbz.gnu" << std::endl;
   strstream << "echo \"splot '$recipfile' using 1:2:3:4:5:6 with vectors head filled lt 1 lw 1 lc rgb \\\"#000000\\\", \\\\\" >> plotbz.gnu" << std::endl;
   strstream << "echo \"'$fbzfile' using 1:2:3:4:5:6 with vector nohead lt 1 lw 1 lc rgb \\\"#000000\\\", \\\\\" >> plotbz.gnu" << std::endl;
   strstream << "echo \"'$bbzfile' using 1:2:3:4:5:6 with vector nohead lt 0 lw 1.5 lc rgb \\\"#000000\\\", \\\\\" >> plotbz.gnu" << std::endl;
   strstream << "echo \"'$kpathfile' using 1:2:3:4:5:6 with vector nohead lt 1 lw 1 lc rgb \\\"#FF0000\\\", \\\\\" >> plotbz.gnu" << std::endl;
   strstream << "if ((Nnonkpath>0))" << std::endl;
   strstream << "then" << std::endl;
   strstream << "echo \"'$nonkpathfile' using 1:2:3:4:5:6 with vector nohead lt 2 lw 1 lc rgb \\\"#FF0000\\\", \\\\\" >> plotbz.gnu" << std::endl;
   strstream << "fi" << std::endl;
   strstream << "echo \"'$kptsfile' using 1:2:3 with points pt 7 ps 0.8 lc rgb \\\"#FF0000\\\" \" >> plotbz.gnu" << std::endl;
   strstream << "" << std::endl;
   strstream << "gnuplot plotbz.gnu" << std::endl;
   strstream << "convert -density 200 bz_$name.eps bz_$name.png" << std::endl;
   strstream << "rm -f plotbz.gnu $datfile $recipfile $kptsfile $kpathfile $nonkpathfile $fbzfile $bbzfile" << std::endl;
   strstream << "" << std::endl;
   strstream << "" << std::endl;
   strstream << "" << std::endl;
   strstream << "" << std::endl; return strstream.str();};
#endif // _AFLOW_GNUPLOT_FUNCS_CPP_
