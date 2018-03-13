// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *                                                                         *
// ***************************************************************************
// Stefano Curtarolo
// Dane Morgan

#ifndef _AFLOW_PFLOW_PRINT_CPP_
#define _AFLOW_PFLOW_PRINT_CPP_

#include "aflow.h"
#include "aflow_pflow.h"

// ***************************************************************************
string _CleanElementName(const string& name) {
  // Dane Morgan - Stefano Curtarolo
  string name1,name2,nameout;
  name1=name+" ";
  name2=name+" ";
  name1=name1.substr(0,1);
  if(!((name1[0]>=65 && name1[0]<=90)||(name1[0]>=97 && name1[0]<=122))) name1="";
  if((name1[0]>=97 && name1[0]<=122)) name1[0]-=-97+65;
  name2=name2.substr(1,1);
  if(!((name2[0]>=65 && name2[0]<=90)||(name2[0]>=97 && name2[0]<=122))) name2="";
  if((name2[0]>=65 && name2[0]<=90)) name2[0]+=-97+65;
  nameout=name1+name2;
  return nameout;
}

// ***************************************************************************
// pflow::PrintACE
// ***************************************************************************
// This funtion prints out structural data in a ace format (for CaRIne).
// Dane Morgan - Stefano Curtarolo
namespace pflow {
  void PrintACE(const xstructure& str, ostream& oss) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "pflow::PrintACE: BEGIN" << endl;
    // Print out data in ace format
    oss.setf(std::ios::fixed,std::ios::floatfield);
    oss.precision(10);
    oss << "CaRIne Crystallography 3.0" << endl;
    oss << "Cell Standard ASCII File=*.ACE" << endl;
    oss << "File Version Number=1" << endl;
    oss << endl;
    //    oss << "Colors Refered to: Mendele" << (char) 239 << "ev Table" << endl;
    oss << "Colors Refered to: RGB Components" << endl;
    oss << endl;
    oss << "----------------------" << endl;
    oss << "Cell Definition" << endl;
    oss << "----------------------" << endl;
    oss << endl;
    oss << "Cell Parameters(" << (char) 197 << " and " << (char) 176 << ")"  << endl;
    // Get cell data
    xvector<double> data(6);
    data=Getabc_angles(str.lattice,DEGREES);data(1)*=str.scale;data(2)*=str.scale;data(3)*=str.scale;
    oss << "a=" << data(1) << '\t'
        << "b=" << data(2) << '\t'
        << "c=" << data(3) << endl;
    oss << "alpha=" << data(4) << '\t'
        << "beta=" << data(5) << '\t'
        << "gamma=" << data(6) << endl;
    oss << endl;
    oss << "System=not used" << endl;
    oss << "Space Group Number=not used" << endl;
    oss << "Number of positions in Cell=" << str.atoms.size() << endl;
    oss << endl;
    // If cartesian convert to direct
    oss << "Atom" << '\t' << "Oxi." << '\t' << "X" << '\t' << "Y" << '\t' << "Z" << '\t'
        << "R(A)" << '\t' << "Occ." << endl;
    for(uint i=0;i<str.atoms.size();i++) {
      oss << (str.atoms.at(i)).cleanname << '\t' << "0" << '\t' << str.atoms.at(i).fpos(1) << '\t' << str.atoms.at(i).fpos(2) << '\t' << str.atoms.at(i).fpos(3) << '\t'
          << "1.00" << '\t' << "1.0000" << endl;
    }
  }
} // namespace pflow

// ***************************************************************************
// pflow::PrintAngles
// ***************************************************************************
// This funtion prints out 3 lattice vectors in cartesian coordinates.
// Dane Morgan - Stefano Curtarolo
namespace pflow {
  void PrintAngles(xstructure str, const double& cutoff, ostream& oss) {  
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "pflow::PrintAngles: BEGIN" << endl;
    uint MAX_NUM_ANGLE=21;
    oss.setf(std::ios::left, std::ios::adjustfield);
    oss.setf(std::ios::fixed, std::ios::floatfield);
    //  oss.setf(std::ios::fixed,std::ios::floatfield);
    //  oss.precision(10);
    deque<deque<_atom> > neigh_mat;
    // [OBSOLETE]  GetStrNeighData(str,cutoff,neigh_mat);
    str.GetStrNeighData(cutoff,neigh_mat);   // once PrintAngles goes in xstructure I can remove the copy
    double tol=1e-15;

    // Output header
    oss << "All atoms identified by" << endl;
    oss << "[basis number]  [name]  [n1 n2 n3 (unit cell)]" << endl;

    // Print out angles
    oss << endl;
    oss << "********************" << " ANGLES " << "********************" << endl;
    oss << "Atoms of increasing indentation are A, B, and C." << endl;
    oss << "The given angles are those made by the lines A to B and B to C (atom B is at the vertex)." << endl;
    oss << "The given angles are always taken 0<=angle<=180." << endl;
    oss << "Angles are only calculated for the first " << MAX_NUM_ANGLE-1 << " neighbors within " << cutoff <<  endl;
    xmatrix<int> ijk(3,3);
    xvector<int> pflag(3);
    for(uint ia1=0;ia1<neigh_mat.size();ia1++) {
      pflag(1)=1;
      _atom a1=neigh_mat.at(ia1).at(0);
      xvector<double> a1pos=a1.cpos;
      for(int i=1;i<=3;i++) ijk(1,i)=a1.ijk(i);
      for(uint ia2=0;(ia2<neigh_mat.at(ia1).size() && ia2<MAX_NUM_ANGLE);ia2++) {
        pflag(2)=1;
        _atom a2=neigh_mat.at(ia1).at(ia2);
        xvector<double> a2pos=a2.cpos;
        for(int i=1;i<=3;i++) ijk(2,i)=a2.ijk(i);
        for(uint ia3=0;(ia3<neigh_mat.at(a2.number).size() && ia3<MAX_NUM_ANGLE);ia3++) {
          pflag(3)=1;
          _atom a3=neigh_mat.at(a2.number).at(ia3);      
          xvector<double> a3pos=a3.cpos;
          for(int i=1;i<=3;i++) ijk(3,i)=a3.ijk(i);
          xvector<double> a1toa2(3),a2toa3(3),a1toa3(3);
          for(int ic=1;ic<=3;ic++) {
            a1toa2(ic)=a2pos(ic)-a1pos(ic);
            a2toa3(ic)=a3pos(ic)-a2pos(ic);
            a1toa3(ic)=a3pos(ic)-a1pos(ic);
          }
          if(modulus(a1toa2)>tol && modulus(a2toa3)>tol && modulus(a1toa3)>tol) {
            double angle = 180.0-rad2deg*acos(getcos(a1toa2,a2toa3));
            if(fabs(getcos(a1toa2,a2toa3)-1)<tol) angle=180;
            if(fabs(getcos(a1toa2,a2toa3)+1)<tol) angle=0;
            if(pflag(1)) {
              oss << setw(4) << a1.number+1
                  << " " << setw(4) << a1.name.c_str()
                  << " " << setw(3) << ijk(1,1)
                  << " " << setw(3) << ijk(1,2)
                  << " " << setw(3) << ijk(1,3)
                  << "   ";
              oss << endl;
            } // if plag(1)
            if(pflag(2)) {
              oss << "     " << setw(4) << a2.number+1
                  << " " << setw(4) << a2.name.c_str()
                  << " " << setw(3) << ijk(2,1)
                  << " " << setw(3) << ijk(2,2)
                  << " " << setw(3) << ijk(2,3)
                  << "   ";
              oss << endl;
            } // if plag(2)
            if(pflag(3)) {
              oss << "          " << setw(4) << a3.number+1
                  << " " << setw(4) << a3.name.c_str()
                  << " " << setw(3) << ijk(3,1)
                  << " " << setw(3) << ijk(3,2)
                  << " " << setw(3) << ijk(3,3)
                  << "   ";
              oss << "     "  << std::setprecision(4) << angle << endl;
            } // if plag(3)
            // Set print flags not to print until next looping.
            for(int i=1;i<=3;i++) pflag(i)=0;	  
          } // if norms>tol
        } // a3
      } // a2
    } // a1
  } // end routine
} // namespace pflow

// ***************************************************************************
// pflow::PrintBands
// ***************************************************************************
// Dane Morgan - Stefano Curtarolo
namespace pflow {
  void CalcPathLength(const pflow::projdata& pd, vector<double>& pathl,
                      vector<double>& dpathl) {
    pflow::matrix<double> rlat=RecipLat(pd.lat);
    pathl[0]=0;
    dpathl[0]=0;
    for(int ik=1;ik<pd.nkpts;ik++) {
      vector<double> disp(3);
      disp=pflow::VVdiff(pd.kpts[ik],pd.kpts[ik-1]);
      dpathl[ik]=dpathl[ik-1]+pflow::norm(disp);
      disp=pflow::vecF2C(rlat,disp);
      pathl[ik]=pathl[ik-1]+pflow::norm(disp);
    }
  }// end
} // namespace pflow

namespace pflow {
  void PrintBands(const pflow::projdata& pd) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "pflow::PrintBands: BEGIN" << endl;
    ofstream outf_up("band.up.out");
    outf_up.precision(5);
    outf_up.setf(std::ios::fixed,std::ios::floatfield);
    outf_up.setf(std::ios::left,std::ios::adjustfield);
    vector<double> pathl(pd.nkpts);
    vector<double> dpathl(pd.nkpts);
    CalcPathLength(pd,pathl,dpathl);
    for(int ik=0;ik<pd.nkpts;ik++) {
      outf_up << pathl[ik] << " " << dpathl[ik];
      for(int ib=0;ib<pd.nbands;ib++) {
        outf_up << " " << pd.ener_k_b_u[ik][ib];
      }
      outf_up << "          " << ik;
      for(int ic=0;ic<3;ic++) {
        outf_up << " " << pd.kpts[ik][ic];
      }
      outf_up<<endl;
    }//kpt
    if(pd.sp) {
      ofstream outf_dn("band.dn.out");
      outf_dn.precision(5);
      outf_dn.setf(std::ios::fixed,std::ios::floatfield);
      outf_dn.setf(std::ios::left,std::ios::adjustfield);
      for(int ik=0;ik<pd.nkpts;ik++) {
        outf_dn << pathl[ik] << " " << dpathl[ik];
        for(int ib=0;ib<pd.nbands;ib++) {
          outf_dn << " " << pd.ener_k_b_d[ik][ib];
        }
        outf_dn << "          " << ik;
        for(int ic=0;ic<3;ic++) {
          outf_dn << " " << pd.kpts[ik][ic];
        }
        outf_dn<<endl;
      }//kpt
    }// if sp
  }// end
} // namespace pflow

// ***************************************************************************
// pflow::PrintCHGCAR
// ***************************************************************************
// Dane Morgan - Stefano Curtarolo
// edited by Corey 1/8/14
namespace pflow {
  bool PrintCHGCAR(const xstructure& str,
                   const stringstream& chgcar_header,
                   const vector<int>& ngrid,
                   const vector<int>& format_dim,
                   const vector<double>& chg_tot,
                   const vector<double>& chg_diff,
                   const string& output_name,
                   ostream& oss) {  

    string soliloquy="pflow::PrintCHGCAR():  ";     // so you know who's talking

    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) oss << soliloquy << "BEGIN" << endl;
    stringstream CHGCARout_ss;
    uint npts,natoms=pflow::GetNumAtoms(str),indx=0,numcolumns;

    // print POSCAR
    if(LDEBUG) oss << soliloquy << "PRINTING POSCAR" << endl;
    CHGCARout_ss << chgcar_header.str();
    // empty line
    CHGCARout_ss << " " << endl;    //put space so AFLOW pick this up as a line
    // grid
    if(ngrid.size()!=3) {
      oss << endl;
      oss << soliloquy << "ERROR: ngrid must be length 3." << endl;
      oss << soliloquy << "Exiting." << endl;
      oss << endl;
      return FALSE;
    }
    if(LDEBUG) oss << soliloquy << "PRINTING GRID" << endl;
    for(uint i=0;i<ngrid.size();i++) {
      CHGCARout_ss << "   " << ngrid.at(i);
    }
    CHGCARout_ss << " " << endl;    //put space so AFLOW pick this up as a line
    npts=ngrid.at(0)*ngrid.at(1)*ngrid.at(2);
    // chg_tot
    if(LDEBUG) oss << soliloquy << "PRINTING CHG_TOT" << endl;
    for(uint i=0;i<npts;i+=format_dim.at(indx)) {
      // make sure we get exactly format_dim.at(indx) 0's and not more
      if(npts-i>=(uint)format_dim.at(indx)) {
        numcolumns=format_dim.at(indx);
      } else {
        numcolumns=npts-i;
      }
      for(uint j=0;j<numcolumns;j++) {    
        //CO - START
        //chg_tot CAN be negative too 
        //http://cms.mpi.univie.ac.at/vasp-forum/viewtopic.php?t=4194&p=4194
        CHGCARout_ss << " " << std::setw(18) << std::setprecision(11) << std::scientific << std::uppercase << chg_tot.at(i+j);
        //CO - END
      }
      CHGCARout_ss << " " << endl;    //put space so AFLOW pick this up as a line
    }
    indx++;
    // augmentation occupancies, all set to 0
    if(format_dim.size()>1) {   // we have CHGCAR, not AECCAR
      if(LDEBUG) oss << soliloquy << "THIS IS A CHGCAR FILE, NOT AN AECCAR FILE" << endl;
      if(LDEBUG) oss << soliloquy << "PRINTING AUGMENTATION OCCUPANCIES (ALL SET TO 0)" << endl;
      for(uint i=1;i<(natoms+1);i++) {
        CHGCARout_ss << "augmentation occupancies " << i << " " << format_dim.at(indx) << endl;
        for(uint j=0;j<(uint)format_dim.at(indx);j+=(uint)format_dim.at(indx+1)) {
          // make sure we get exactly format_dim.at(indx) 0's and not more
          if((uint)format_dim.at(indx)-j>=(uint)format_dim.at(indx+1)) {
            numcolumns=format_dim.at(indx+1);
          } else {
            numcolumns=(uint)format_dim.at(indx)-j;
          }
          for(uint k=0;k<numcolumns;k++) {
            CHGCARout_ss << "  " << std::setprecision(7) << std::scientific << std::uppercase << 0.0;
          }
          CHGCARout_ss << " " << endl;    //put space so AFLOW pick this up as a line
        }
        indx+=2;
      }
    }
    if(chg_diff.size()>0) {     // we have spin-polarized
      if(LDEBUG) oss << soliloquy << "THIS IS A SPIN-POLARIZED CALCULATION" << endl;
      // add line of 0's
      for(uint i=0;i<(uint)format_dim.at(indx);i++) { //mimic formatting of chg_diff for 0's line
        // standard shows +1 more 0
        CHGCARout_ss << "  " << std::setprecision(12) << std::scientific << std::uppercase << 0.0;
      }
      CHGCARout_ss << " " << endl;    //put space so AFLOW pick this up as a line
      // grid
      if(LDEBUG) oss << soliloquy << "PRINTING GRID" << endl;
      for(uint i=0;i<ngrid.size();i++) {
        CHGCARout_ss << "   " << ngrid.at(i);
      }
      CHGCARout_ss << " " << endl;    //put space so AFLOW pick this up as a line
      // chg_diff
      if(LDEBUG) oss << soliloquy << "PRINTING CHG_DIFF" << endl;
      for(uint i=0;i<npts;i+=format_dim.at(indx)) {
        // make sure we get exactly format_dim.at(indx) 0's and not more
        if(npts-i>=(uint)format_dim.at(indx)) {
          numcolumns=format_dim.at(indx);
        } else {
          numcolumns=npts-i;
        }
        for(uint j=0;j<numcolumns;j++) {    
          // chg_diff needs more space for minus sign
          CHGCARout_ss << " " << std::setw(18) << std::setprecision(11) << std::scientific << std::uppercase << chg_diff.at(i+j);
        }
        CHGCARout_ss << " " << endl;    //put space so AFLOW pick this up as a line
      }
      indx++;
      // check for augmentation occupancies
      if(format_dim.size()>indx) { // we have more augmentation occupancies to write out, all set to 0
        if(LDEBUG) oss << soliloquy << "PRINTING AUGMENTATION OCCUPANCIES (ALL SET TO 0)" << endl;
        for(uint i=1;i<(natoms+1);i++) {
          CHGCARout_ss << "augmentation occupancies " << i << " " << format_dim.at(indx) << endl;
          for(uint j=0;j<(uint)format_dim.at(indx);j+=(uint)format_dim.at(indx+1)) {
            // make sure we get exactly format_dim.at(indx) 0's and not more
            if((uint)format_dim.at(indx)-j>=(uint)format_dim.at(indx+1)) {
              numcolumns=format_dim.at(indx+1);
            } else {
              numcolumns=(uint)format_dim.at(indx)-j;
            }
            for(uint k=0;k<numcolumns;k++) {
              CHGCARout_ss << "  " << std::setprecision(7) << std::scientific << std::uppercase << 0.0;
            }
            CHGCARout_ss << " " << endl;    //put space so AFLOW pick this up as a line
          }
          indx+=2;
        }

      }
    }

    // check if output file exists, delete it first
    if(aurostd::FileExist(output_name)) {
      if(LDEBUG) oss << soliloquy << "DELETING EXISTING OUTPUT FILE" << endl;
      aurostd::RemoveFile(output_name);
    }

    // write out to file
    if(LDEBUG) oss << soliloquy << "WRITING OUTPUT FILE" << endl;
    aurostd::stringstream2file(CHGCARout_ss,output_name);
    if(LDEBUG) oss << soliloquy << "DONE" << endl;
    return TRUE;
  }
} // namespace pflow

// ***************************************************************************
// pflow::PrintChgInt
// ***************************************************************************
// This funtion prints out 3 lattice vectors in cartesian coordinates.
// Dane Morgan - Stefano Curtarolo
namespace pflow {
  void PrintChgInt(vector<pflow::matrix<double> >& rad_chg_int,
                   pflow::matrix<double>& vor_chg_int, ostream& oss) {  
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "pflow::PrintChgInt: BEGIN" << endl;
    oss.setf(std::ios::left, std::ios::adjustfield);
    oss.setf(std::ios::fixed, std::ios::floatfield);

    // Output Voronoi data
    oss << "******************** Charge in Voronoi volume: Tot Diff Up Dn ********************" << endl;
    for(uint iat=0;iat<vor_chg_int.size();iat++) {
      oss << iat+1;
      for(int i=0;i<4;i++) {
        oss << " " << vor_chg_int[iat][i];
      }
      oss << endl;
    }

    // Output Radial data
    oss << "******************** Charge in sphere: Tot Diff Up Dn ********************" << endl;
    for(uint iat=0;iat<rad_chg_int.size();iat++) {
      oss << "Atom " << iat+1 << endl;
      for(uint ib=0;ib<rad_chg_int[iat].size();ib++) {
        oss << "  " << rad_chg_int[iat][ib][0];
        for(int i=1;i<5;i++) {
          oss << " " << rad_chg_int[iat][ib][i];
        }
        oss << endl;
      }
    }
  } // end routine
} // namespace pflow

// ***************************************************************************
// pflow::PrintCIF
// ***************************************************************************
// This funtion prints out structural data in a cif format.
// Dane Morgan - Stefano Curtarolo
namespace pflow {
  void PrintCIF(ostream& oss,const xstructure& str,int _spacegroupnumber) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "pflow::PrintCIF: BEGIN" << endl;
    // Print out data in cif format
    oss.setf(std::ios::fixed,std::ios::floatfield);
    oss.precision(10);
    int spacegroupnumber=_spacegroupnumber;
    if(spacegroupnumber<1 || spacegroupnumber>230) spacegroupnumber=1;
    //title
    //CO MODS START 170613 - fixes for BOB cif
    string compound="",pd_name="";
    string tmp;
    if(!str.title.empty()){
      vector<string> tokens;
      aurostd::string2tokens(str.title,tokens," ");
      if(tokens.size()&&!tokens[0].empty()){
        if(aurostd::substring2bool(tokens[0],"/")){ //something out of LIBX: Ru_pvTa_pvY_sv/TFCC012.ABC, ICSD looks like this: Ag1F6Sb1
          tmp=tokens[0];
          aurostd::string2tokens(tmp,tokens,"/");
          compound=tokens[0];
          if(tokens.size()>1){pd_name=tokens[1];}  //just safety
        }else{
          compound=pd_name=tokens[0];
        }
      }
    }
    if(compound.empty()){
      if(!pd_name.empty()){compound=pd_name;}
      else{compound="1";}
    }
    //oss << "data_" << str.title << endl;
    oss << "# AFLOW.org Repositories" << endl;  //CO fix 170606
    oss << "# " << str.title << endl;  //CO fix 170606
    oss << "data_" << compound << endl;  //CO fix 170606
    //CO try to get phase name, usually written as first part of title in AFLOW
    if(!pd_name.empty()){oss << "_pd_phase_name " << pd_name << endl;}  //just safety
    //CO MODS END 170613 - fixes for BOB cif
    xvector<double> data(6);
    data=Getabc_angles(str.lattice,DEGREES);data(1)*=str.scale;data(2)*=str.scale;data(3)*=str.scale;
    oss << "_cell_length_a  " << setw(12) << data(1) <<endl;
    oss << "_cell_length_b  " << setw(12) << data(2) <<endl;
    oss << "_cell_length_c  " << setw(12) << data(3) <<endl;
    oss << "_cell_angle_alpha  " << setw(9) << data(4) <<endl;
    oss << "_cell_angle_beta  " << setw(9) << data(5) <<endl;
    oss << "_cell_angle_gamma  " << setw(9) << data(6) <<endl;
    // Output space group - Use P1 for now
    oss << "_symmetry_space_group_name_H-M  '" << GetSpaceGroupName(spacegroupnumber) << "'" << endl;
    oss << "_symmetry_Int_Tables_Number  " << spacegroupnumber << endl;
    oss << "loop_" << endl;
    oss << "_symmetry_equiv_pos_site_id" << endl;
    //oss << "_symmetry_equiv_pos_as_xyz_" << endl;
    oss << "_symmetry_equiv_pos_as_xyz" << endl;  //CO fix 170606
    oss << "  1  x,y,z" << endl;

    //atom posns header
    oss << "loop_" << endl;
    oss << " _atom_site_label" << endl;
    oss << " _atom_site_occupancy" << endl;
    oss << " _atom_site_fract_x" << endl;
    oss << " _atom_site_fract_y" << endl;
    oss << " _atom_site_fract_z" << endl;
    oss << " _atom_site_thermal_displace_type" << endl;
    oss << " _atom_site_B_iso_or_equiv" << endl;
    oss << " _atom_site_type_symbol" << endl;

    // Coordinates in direct

    oss.precision(10);
    for(uint i=0;i<str.atoms.size();i++) {
      oss << str.atoms.at(i).cleanname.c_str() << i+1 << " 1.0 "
          << setw(12) << str.atoms.at(i).fpos(1) <<" "
          << setw(12) << str.atoms.at(i).fpos(2) <<" "
          << setw(12) << str.atoms.at(i).fpos(3) <<" "
          << " Biso 1.0 " << str.atoms.at(i).cleanname.c_str() <<endl;
    }
  }
} // namespace pflow

// ***************************************************************************
// pflow::PrintClat
// ***************************************************************************
// This funtion prints out 3 lattice vectors in cartesian coordinates.
// Dane Morgan - Stefano Curtarolo
namespace pflow {
  void PrintClat(const xvector<double>& data, ostream& oss) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "pflow::PrintClat: BEGIN" << endl;
    if(data.rows!=6) {
      init::ErrorOption(cout,"","pflow::PrintClat",aurostd::liststring2string("data.size()=",aurostd::utype2string(data.rows)));
      exit(0);
    }
    oss.setf(std::ios::fixed,std::ios::floatfield);
    oss.precision(10);
    xvector<double> abc_angles(data);
    xmatrix<double> lattice = GetClat(abc_angles); // angles in radiants
    roundoff(lattice,1.0e-4);
    oss << lattice(1,1) << " " << lattice(1,2) << " " << lattice(1,3) << endl;
    oss << lattice(2,1) << " " << lattice(2,2) << " " << lattice(2,3) << endl;
    oss << lattice(3,1) << " " << lattice(3,2) << " " << lattice(3,3) << endl;
    if(LDEBUG) cerr << "pflow::PrintClat: BEGIN" << endl;
  }
} // namespace pflow

// ***************************************************************************
// pflow::PrintCmpStr
// ***************************************************************************
// Prints out the comparison between two structures.
// Dane Morgan - Stefano Curtarolo
namespace pflow {
  void PrintCmpStr(const xstructure& str1, const xstructure& str2, const double& rcut, ostream& oss) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "pflow::PrintCmpStr: BEGIN" << endl;
    oss.setf(std::ios::left, std::ios::adjustfield);
    oss.setf(std::ios::fixed, std::ios::floatfield);
    oss.precision(4);
    double cutoff=rcut;
    //  oss.width(10);
    // Rescale to avoid scale problems
    xstructure sstr1=str1;
    xstructure sstr2=str2;
    sstr1=ReScale(sstr1,1.0);
    sstr2=ReScale(sstr2,1.0);
    int nat1=sstr1.atoms.size();
    int nat2=sstr2.atoms.size();
    std::deque<int> num_each_type1=sstr1.num_each_type;
    std::deque<int> num_each_type2=sstr2.num_each_type;
    int ntyp_min;
    int ntyp_max;
    int nprs_min;
    int nprs_max;

    oss << "*************************  Structure Comparison *************************" << endl;

    // Nun atoms
    oss << "Number of atoms" << endl;
    oss << "  Num Atom Str1: " << nat1 << endl;
    oss << "  Num Atom Str2: " << nat2 << endl;
    oss << "  Num Atom %Err: " << 100*2*(nat2-nat1)/(nat2+nat1) << endl;

    // Num types
    int ntyp1=sstr1.num_each_type.size();
    int ntyp2=sstr2.num_each_type.size();
    oss << "Number of types" << endl;
    oss << "  Num Type Str1: " << ntyp1 << endl;
    oss << "  Num Type Str2: " << ntyp2 << endl;
    oss << "  Num Type Err : " << ntyp2-ntyp1 << endl;

    // Stoichiometry
    int wid=8;
    ntyp_max=max(ntyp1,ntyp2);
    vector<double> tmp_stoich1(ntyp_max,0);
    vector<double> tmp_stoich2(ntyp_max,0);
    oss << "Stoichiometry" << endl;
    oss << "  Type:        ";
    for(int it=0;it<ntyp_max;it++) {
      oss << setw(wid) << it << " ";
    }
    oss << endl;
    oss << "  Stoich Str1: ";
    for(int it=0;it<ntyp1;it++) {
      tmp_stoich1[it]=(double) num_each_type1[it]/ (double) nat1;
      oss << setw(wid) << tmp_stoich1[it] << " ";
    }
    oss << endl;
    oss << "  Stoich Str2: ";
    for(int it=0;it<ntyp2;it++) {
      tmp_stoich2[it]=(double) num_each_type2[it]/ (double) nat2;
      oss << setw(wid) << tmp_stoich2[it] << " ";
    }
    oss << endl;
    oss << "  Stoich Diff: ";
    for(int it=0;it<ntyp_max;it++) {
      oss << setw(wid) << tmp_stoich2[it]-tmp_stoich1[it] << " ";
    }
    oss << endl;

    // Volume
    double v1=GetVol(sstr1.lattice)/nat1;
    double v2=GetVol(sstr2.lattice)/nat2;
    oss << "Volume per atom" << endl;
    oss << "  Vol Str1: " << v1 << endl;
    oss << "  Vol Str2: " << v2 << endl;
    oss << "  Vol %Err: " << 100*2*(v2-v1)/(v2+v1) << endl;

    // lat par and angles
    // vector<double> ldat1=pflow::Getabc_angles(pflow::xmatrix2matrix(sstr1.lattice));
    vector<double> ldat1=xvector2vector(Getabc_angles(sstr1.lattice,DEGREES));
    // vector<double> ldat2=pflow::Getabc_angles(pflow::xmatrix2matrix(sstr2.lattice));
    vector<double> ldat2=xvector2vector(Getabc_angles(sstr2.lattice,DEGREES));

    ldat1=pflow::Sort_abc_angles(ldat1);  
    ldat2=pflow::Sort_abc_angles(ldat2);  
    vector<double> ldat_err(7,0.0);
    ldat_err[0]=100*2*(v2*nat2-v1*nat1)/(v2*nat2+v1*nat1);
    for(int i=0;i<6;i++) {
      double tmp=(abs(ldat2[i])+abs(ldat1[i]));
      if(tmp>0) {
        ldat_err[i+1]=100*2*(ldat2[i]-ldat1[i])/(ldat2[i]+ldat1[i]);
      }
      else{
        ldat_err[i+1]=0.0;
      }
    }
    oss << "Unit cell: tot_vol a b c alpha beta gamma" << endl;
    oss << "  Cell Str1: ";
    oss << v1*nat1 << " ";
    pflow::Vout(ldat1,oss);
    oss << "  Cell Str2: ";
    oss << v2*nat2 << " ";
    pflow::Vout(ldat2,oss);
    oss << "  Cell %Err: ";
    pflow::Vout(ldat_err,oss);

    // Volume adjusted distances
    double vavg=(v1+v2)/2.0;
    double rescale1=std::pow((double) vavg/v1,(double) 1.0/3.0);
    double rescale2=std::pow((double) vavg/v2,(double) 1.0/3.0);
    xstructure str1_newvol=sstr1;
    xstructure str2_newvol=sstr2;
    str1_newvol=SetScale(str1_newvol,rescale1);
    str2_newvol=SetScale(str2_newvol,rescale2);
    pflow::matrix<double> dist1;
    pflow::matrix<double> dist2;
    pflow::matrix<double> dist_diff;
    pflow::matrix<double> dist_diff_n;
    // If cutoff < 0 use default of 4*(avg rad of an atom).  Note
    // that in a rescaled str where each atom has vavg=1, rad=.6204.
    // In the true str, where each atom has vol vavg, rad=0.6204*vavg^(1/3).
    if(rcut<0) cutoff=4*0.6204*std::pow((double) vavg,(double) 1.0/3.0);
    pflow::CmpStrDist(str1_newvol,str2_newvol,cutoff,dist1,dist2,dist_diff,dist_diff_n);
    ntyp_min=min(ntyp1,ntyp2);
    nprs_min=dist_diff.size();
    oss << "Avg of Distance Magnitude Differences" << endl;
    oss << "  Rescaling vol/atom of Str1 and Str2 to the average vol/atom: " << vavg << endl;
    oss << "  Working with the all the types up to num types: " << ntyp_min << endl;
    oss << "  This corresponds to the following number of distinct pairs: " << nprs_min << endl;
    double tadiff=0;
    double tadiff_n=0;
    int tnn=0;
    for(int it1=0;it1<ntyp_min;it1++) {
      for(int it2=it1;it2<ntyp_min;it2++) {
        int id=it2-it1+ntyp_min*it1-max(0,it1*(it1-1)/2);
        int nn=dist_diff[id].size();      
        oss << "  P" << id << " Pair between types: "
            << it1 << " " << it2 <<endl;
        oss << "    P" << id << " Comparing neighbors within cutoff (number): "
            << cutoff << " (" << nn << ")" <<  endl;
        double adiff=0;
        double adiff_n=0;
        // double perc_err=0;  // DANE not used
        for(int in=0;in<nn;in++) {
          adiff=adiff+abs(dist_diff[id][in]);
          adiff_n=adiff_n+abs(dist_diff_n[id][in]);
        }
        tadiff=tadiff+adiff;
        tadiff_n=tadiff_n+adiff_n;
        if(nn>0) {
          adiff=adiff/nn;
          adiff_n=adiff_n/nn;
        }
        else{
          adiff=0.0;
          adiff_n=0.0;
        }
        tnn=tnn+nn;
        oss << "    P" << id << " Avg Difference: " << adiff << endl;
        oss << "    P" << id << " Avg Diff/vol_avg^1/3 %: "
            << 100*adiff/std::pow((double) vavg,(double) 1.0/3.0) << endl;
        oss << "    P" << id << " Avg Difference/bond length %: " << 100*adiff_n << endl;
      }// it2
    }// it1
    if(tnn>0) {
      tadiff=tadiff/tnn;
      tadiff_n=tadiff_n/tnn;
    }
    else{
      tadiff=0;
      tadiff_n=0;
    }
    oss << "  Totals for all pairs, compared by type, within cutoff (number): "
        << cutoff << " (" << tnn << ")" <<  endl;
    oss << "    Tot Avg Difference: " << tadiff << endl;
    oss << "    Tot Avg Diff/vol_avg^1/3 %: "
        << 100*tadiff/std::pow((double) vavg,(double) 1.0/3.0) << endl;
    oss << "    Tot Avg Diff/bond_length %: "
        << 100*tadiff_n << endl;

    /* THIS NEEDS WORK - need to sum the rdfsh over all atoms of one type,
       then output average numbers.  This is still not perfect, since you could
       get cancelling errors.  Could try to do site to site comparison.
       // Compare radial distribution functions
       ntyp_max=max(ntyp1,ntyp2);
       nprs_max=(ntyp_max+1)*(ntyp_max+2)/2;
       vector<double> tmpsh1(nprs_max,0);
       vector<double> tmpsh2(nprs_max,0);
       double bin_width=0.01;
       int nbins=Nint(cutoff/bin_width);
       int smooth_width=5;
       int nsh_max=2;
       pflow::matrix<double> rdf_all_1;
       pflow::matrix<double> rdf_all_2;
       GetRDF(sstr1,cutoff,nbins,rdf_all_1);
       GetRDF(sstr2,cutoff,nbins,rdf_all_2);
       pflow::matrix<double> rdf_all_1_sm=GetSmoothRDF(rdf_all_1,smooth_width);
       pflow::matrix<double> rdf_all_2_sm=GetSmoothRDF(rdf_all_2,smooth_width);
       // Get shells
       pflow::matrix<double> rdfsh_all_1;
       pflow::matrix<double> rdfsh_loc_1; // Radial location of rdf shells.
       pflow::matrix<double> rdfsh_all_2;
       pflow::matrix<double> rdfsh_loc_2; // Radial location of rdf shells.
       GetRDFShells(sstr1,cutoff,nbins,smooth_width,
       rdf_all_1_sm,rdfsh_all_1,rdfsh_loc_1);
       GetRDFShells(sstr2,cutoff,nbins,smooth_width,
       rdf_all_2_sm,rdfsh_all_2,rdfsh_loc_2);
       // print out position, number for each radial shell up to nsh_max.
       // Print out for each type of pair
       oss << "Radial shell comparison" << endl;
       for(int ish=0;ish<nsh_max;ish++) { // Radial shells loop.
       oss << "  Shell " << ish << endl;
       oss << "                    ";
       for(int it1=0;it1<ntyp_max;it1++) {
       for(int it2=it1;it2<ntyp_max;it2++) {
       oss << it1 << "-" << it2 << "  ";
       }
       }
       oss << endl;
       oss << "    Shell Occ Str1: ";
       for(int it1=0;it1<ntyp1;it1++) {
       for(int it2=0;it2<ntyp1;it2++) {
       int id=it2+(ntyp1+1)*it1;
       int idb=ntyp1*nat1+it2+(ntyp1+1)*it1;
       cout << id<<" "<<idb<< " "<<rdfsh_all_1.size()<<" "<<ntyp1<<" "<< nat1<< endl;
       tmpsh1[id]=rdfsh_all_1[ntyp1*nat1+id][ish];
       oss << tmpsh1[id] << " ";
       }
       }
       oss << endl;
       }    

    */

    // Number of each type up to rcut.
    nprs_max=ntyp_max*(ntyp_max+1)/2;
    vector<int> n1(nprs_max,0);
    vector<int> n2(nprs_max,0);
    oss << "Number neighbors of each type within cutoff for rescaled volume" << endl;
    oss << "  Cutoff: " << cutoff << endl;
    oss << "  Rescaling vol/atom of Str1 and Str2 to the average vol/atom: " << vavg << endl;
    oss << "           ";
    for(int it1=0;it1<ntyp_max;it1++) {
      for(int it2=it1;it2<ntyp_max;it2++) {
        oss << it1 << "-" << it2 << "  ";
      }
    }
    oss << endl;
    oss << "  NN Str1  ";
    for(int it1=0;it1<ntyp_max;it1++) {
      for(int it2=it1;it2<ntyp_max;it2++) {
        if(it1<ntyp1 && it2<ntyp1) {
          int id1=it2-it1+ntyp1*it1-max(0,it1*(it1-1)/2);
          n1[id1]=dist1[id1].size();
          oss << setw(4) << n1[id1] << " ";
        }
      }
    }
    oss << endl;
    oss << "  NN Str2  ";
    for(int it1=0;it1<ntyp_max;it1++) {
      for(int it2=it1;it2<ntyp_max;it2++) {
        if(it1<ntyp2 && it2<ntyp2) {
          int id2=it2-it1+ntyp2*it1-max(0,it1*(it1-1)/2);
          n2[id2]=dist2[id2].size();
          oss << setw(4) << n2[id2] << " ";
        }
      }
    }
    oss << endl;  
    oss << "  NN Err   ";
    for(int it1=0;it1<ntyp_max;it1++) {
      for(int it2=it1;it2<ntyp_max;it2++) {
        int id=it2-it1+ntyp_max*it1-max(0,it1*(it1-1)/2);
        oss << setw(4) << n2[id]-n1[id] << " ";
      }
    }
    oss << endl;  

    // Space group
    bool Platon_EQUAL=FALSE,Platon_EXACT=FALSE;
    double Platon_ang=1.0,Platon_d1=0.25,Platon_d2=0.25,Platon_d3=0.25;
    system("touch pflow_tmp_in"); // Create input file to connect to stream.
    ofstream outfile("pflow_tmp_out");
    ifstream infile("pflow_tmp_in");
    string spcgrp1a,spcgrp2a,spcgrp1b,spcgrp2b;
    platon2print(str1,Platon_EQUAL,Platon_EXACT,Platon_ang,Platon_d1,Platon_d2,Platon_d3,outfile);
    system("cat pflow_tmp_out | platonSG > pflow_tmp_in");
    infile >> spcgrp1a >> spcgrp1b;
    // Puts file position pointers back to beginning.
    outfile.seekp(0);
    infile.seekg(0);
    platon2print(sstr2,Platon_EQUAL,Platon_EXACT,Platon_ang,Platon_d1,Platon_d2,Platon_d3,outfile);
    system("cat pflow_tmp_out | platonSG > pflow_tmp_in");
    infile >> spcgrp2a >> spcgrp2b;
    system("/bin/rm pflow_tmp_out");
    system("/bin/rm pflow_tmp_in");
    oss << "Space groups " << endl;
    oss << "  SG Str1: " << spcgrp1a << " " << spcgrp1b << endl;
    oss << "  SG Str2: " << spcgrp2a << " " << spcgrp2b << endl;

    // Center of mass  
  } // end routine
} // namespace pflow

// ***************************************************************************
// pflow::PrintData
// ***************************************************************************
// This funtion prints out structural data.
// Stefano Curtarolo
namespace pflow {
  void PrintData(const xstructure& str,xstructure& str_sym,xstructure& str_sp,xstructure& str_sc, ostream& oss,string smode, const string& format,bool already_calculated) { //CO171027
    //DX 12/15/18 [OBSOLETE] double tolerance = AUROSTD_NAN;
    double tolerance = SYM::defaultTolerance(str_sym); //DX 2/15/18 - start with default tolerance //DX 2/26/18 - str_sym not str_sp
    if(already_calculated){tolerance=str_sp.sym_eps;} //CO 171025
    bool no_scan = false;
    PrintData(str,str_sym,str_sp,str_sc,oss,smode,tolerance,no_scan,format,already_calculated); //CO171027
  }
}

namespace pflow {
  void PrintData(const xstructure& str,xstructure& str_sym,xstructure& str_sp,xstructure& str_sc, ostream& oss_final,string smode, double tolerance, bool no_scan, const string& format,bool already_calculated) { //CO171027
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "pflow::PrintData: BEGIN" << endl;
    // smode=="DATA" or "EDATA"
    // Print out structural data
    stringstream oss;
    oss.setf(std::ios::fixed,std::ios::floatfield);
    double vol;
    vol=GetVol(str.lattice)*str.scale*str.scale*str.scale;

    bool symmetry_commensurate = false;
    bool force_perform = true;  // Force a result to be found at the end, even if incommensurate 
    double orig_tolerance = tolerance;
    uint sym_eps_change_count = 0; //DX 2/26/18 - added sym eps change count

    xstructure str_aus(str);//,str_sp,str_sc; //CO171027
    while(symmetry_commensurate == false){
    if(smode!="EDATA" || no_scan){
      symmetry_commensurate=true;
    }
    // FORMAT = TXT
    if(format=="txt"){
    oss << "REAL LATTICE" << endl;
    oss << " Real space a b c alpha beta gamma: ";
    oss.precision(10);
    xvector<double> data(6);
    data=Getabc_angles(str.lattice,DEGREES);data(1)*=str.scale;data(2)*=str.scale;data(3)*=str.scale;
    oss.precision(10);oss << data(1) << " " << data(2) << " " << data(3) << " ";
    oss.precision(3); oss << data(4) << " " << data(5) << " " << data(6) << endl;
    oss.precision(4);
    if(smode=="EDATA") {
      oss << " Real space a b c alpha beta gamma: ";
      oss.precision(10);oss << data(1)*angstrom2bohr << " " << data(2)*angstrom2bohr << " " << data(3)*angstrom2bohr << " ";
      oss.precision(3); oss << data(4) << " " << data(5) << " " << data(6) << "   Bohrs/Degs " << endl;
      oss.precision(4);
    }
    oss << " Real space Volume: " << vol << endl;
    oss << " Real space c/a = " << data(3)/data(1) << endl;
        str_aus=str;//,str_sp,str_sc; //CO171027
      //DX 9/1/17 -Add tolerance - START
      str_aus.sym_eps = str_sp.sym_eps = str_sc.sym_eps = tolerance; 
      str_aus.sym_eps_change_count = str_sp.sym_eps_change_count = str_sc.sym_eps_change_count = sym_eps_change_count; //DX 2/26/18 - added sym eps change count
      //DX 9/1/17 -Add tolerance - END
    if(smode=="EDATA") {
      // str_aus.CalculateSymmetryPointGroup(FALSE);
      // str_aus.CalculateSymmetryFactorGroup(FALSE);
      // str_aus.CalculateSymmetryPointGroupCrystal(FALSE);
          if(!already_calculated){  //CO 171025
      str_aus.GetLatticeType(str_sp,str_sc);
      if(orig_tolerance == AUROSTD_NAN){
        orig_tolerance = str_sp.sym_eps;  // Use new tolerance calculated here 
      }
      tolerance = str_sp.sym_eps;
      sym_eps_change_count = str_aus.sym_eps_change_count = str_sp.sym_eps_change_count; //DX 2/26/18 - added sym eps change count
          }
      oss << "BRAVAIS LATTICE OF THE CRYSTAL (pgroup_xtal)" << endl;
      oss << " Real space: Bravais Lattice Primitive        = " << str_aus.bravais_lattice_type << endl;// " " << str.title << endl;
      oss << " Real space: Lattice Variation                = " << str_aus.bravais_lattice_variation_type << endl;//wahyu mod
      oss << " Real space: Lattice System                   = " << str_aus.bravais_lattice_system << endl;
      oss << " Real space: Pearson Symbol                   = " << str_aus.pearson_symbol << endl;
      oss << "POINT GROUP CRYSTAL" << endl;
      oss << " Real space: Crystal Family                   = " << str_aus.crystal_family << endl;
      oss << " Real space: Crystal System                   = " << str_aus.crystal_system << endl;
      oss << " Real space: Crystal Class                    = " << str_aus.point_group_crystal_class << endl;
      oss << " Real space: Point Group (Hermann Mauguin)    = " << str_aus.point_group_Hermann_Mauguin << endl; // "      PGXTAL" << endl;
      oss << " Real space: Point Group (Schoenflies)        = " << str_aus.point_group_Shoenflies << endl;
      oss << " Real space: Point Group Orbifold             = " << str_aus.point_group_orbifold << endl;
      oss << " Real space: Point Group Type                 = " << str_aus.point_group_type << endl;
      oss << " Real space: Point Group Order                = " << str_aus.point_group_order << endl;
      oss << " Real space: Point Group Structure            = " << str_aus.point_group_structure << endl;
          PrintSGData(str_aus,tolerance,oss,no_scan,false,format,already_calculated);  //DX 8/31/17 - SGDATA //CO 171025 //CO171027
      sym_eps_change_count = str_aus.sym_eps_change_count; //DX 2/26/18 - added sym eps change count
        //DX 9/1/17 - Add consistency check for input symmetry method and ITC method - START
        int multiplicity_of_primitive=str_sp.fgroup.size()/str_sp.pgroup_xtal.size();
        bool derivative_structure=false;
        if(SYM::ComparePointGroupAndSpaceGroupString(str_aus,multiplicity_of_primitive,derivative_structure) || no_scan){
          symmetry_commensurate=true;
        }
        else{
          if(LDEBUG){
            cerr << "pflow::PrintData: WARNING: Space group symbol and point group symbol do not match. (sg=" << GetSpaceGroupName(str_aus.space_group_ITC) << ", pg=" << str_aus.point_group_Hermann_Mauguin << ")"  << endl;
          }
          str_aus.ClearSymmetry();
          oss.clear();
          //if(!SYM::change_tolerance(str_aus,tolerance,orig_tolerance,change_sym_count,str_aus.dist_nn_min,no_scan)){
          if(!SYM::change_tolerance(str_aus,tolerance,str_aus.dist_nn_min,no_scan)){
            if(force_perform){
              if(LDEBUG){
                cerr << "pflow::PrintData: WARNING: Scan failed. Reverting back to original tolerance and recalculating as is (with aforementioned inconsistencies)." << str.directory << endl;
              }
              no_scan=true; //Force it to continue 
              tolerance = orig_tolerance;
              str_aus.ClearSymmetry();
              oss.str("");
              continue;
            }
          }
          else{
            str_aus.ClearSymmetry();
            oss.str("");
            sym_eps_change_count = str_aus.sym_eps_change_count; //DX 2/26/18 - added sym eps change count
            continue;
          }
        }
        //DX 9/1/17 - Add consistency check for input symmetry method and ITC method - END
      oss << "BRAVAIS LATTICE OF THE LATTICE (pgroup)" << endl;
      oss << " Real space: Bravais Lattice Primitive        = " << str_aus.bravais_lattice_lattice_type << endl;// " " << str.title << endl;
      oss << " Real space: Lattice Variation                = " << str_aus.bravais_lattice_lattice_variation_type << endl;//wahyu mod
      oss << " Real space: Lattice System                   = " << str_aus.bravais_lattice_lattice_system << endl;
      oss << "SUPERLATTICE (equally decorated)" << endl;
      oss << " Real space: Bravais Superlattice Primitive   = " << str_aus.bravais_superlattice_type << endl;
      oss << " Real space: Superlattice Variation           = " << str_aus.bravais_superlattice_variation_type << endl;
      //DX oss << " Real space: Superlattice System              = " << str_aus.bravais_lattice_system << endl;
      oss << " Real space: Superlattice System              = " << str_aus.bravais_superlattice_system << endl; //DX - fixed mistake in line above
      oss << " Real space: Pearson Symbol Superlattice      = " << str_aus.pearson_symbol_superlattice << endl;
    }
    oss << "RECIPROCAL LATTICE" << endl;
    oss << " Reciprocal space lattice:" << endl;
    oss.precision(10);
    oss <<"   "<<str.klattice(1,1)<<"  "<<str.klattice(1,2) <<" "<<str.klattice(1,3)<<endl;
    oss <<"   "<<str.klattice(2,1)<<"  "<<str.klattice(2,2) <<" "<<str.klattice(2,3)<<endl;
    oss <<"   "<<str.klattice(3,1)<<"  "<<str.klattice(3,2) <<" "<<str.klattice(3,3)<<endl;
    oss << " Reciprocal space a b c alpha beta gamma: ";
    data=Getabc_angles(str.klattice,DEGREES);
    oss.precision(10);oss << data(1) << " " << data(2) << " " << data(3) << " ";
    oss.precision(3); oss << data(4) << " " << data(5) << " " << data(6) << endl;
    oss.precision(4);
    double kvol=GetVol(str.klattice);
    oss << " Reciprocal space Volume: " << kvol << endl;
    if(smode=="EDATA") {
      oss << " Reciprocal lattice primitive            = " << str_aus.reciprocal_lattice_type << endl;
      oss << " Reciprocal lattice variation            = " << str_aus.reciprocal_lattice_variation_type << endl;//wahyu mod
      //oss << " Reciprocal conventional lattice         = " << str_aus.reciprocal_conventional_lattice_type << endl;
      oss << "SPRIM" << endl;
      oss << str_sp << endl;
      oss << "SCONV" << endl;
      oss << str_sc << endl;
    }
  }
    // FORMAT = JSON
    if(format=="json"){
      string eendl="";
      bool roff=true; //round off
      bool PRINT_NULL=FALSE;
      stringstream sscontent_json;
      vector<string> vcontent_json;
      
      // Real lattice parameters
      xvector<double> data(6);
      data=Getabc_angles(str.lattice,DEGREES);data(1)*=str.scale;data(2)*=str.scale;data(3)*=str.scale;
      if(data.rows){
        sscontent_json << "\"lattice_parameters\":[" << aurostd::joinWDelimiter(xvecDouble2vecString(data,5,roff),",") << "]" << eendl;
      } else{
        if(PRINT_NULL){ sscontent_json << "\"lattice_parameters\":null" << eendl;}
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
     
      // Real lattice parameters (Bohr/Deg)    
      if(smode=="EDATA") {
        xvector<double> Bohr_Degs_data(6);
        Bohr_Degs_data(1) =  data(1)*angstrom2bohr;
        Bohr_Degs_data(2) =  data(2)*angstrom2bohr;
        Bohr_Degs_data(3) =  data(3)*angstrom2bohr;
        Bohr_Degs_data(4) =  data(4);
        Bohr_Degs_data(5) =  data(5);
        Bohr_Degs_data(6) =  data(6);
        if(data.rows){
          sscontent_json << "\"lattice_parameters_Bohr_deg\":[" << aurostd::joinWDelimiter(xvecDouble2vecString(Bohr_Degs_data,5,roff),",") << "]" << eendl;
        } else{
          if(PRINT_NULL){ sscontent_json << "\"lattice_parameters_Bohr_deg\":null" << eendl;}
        }
        vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      }
      
      // Real space volume
      sscontent_json << "\"volume\":" << vol << eendl;
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      
      // Real space c/a
      sscontent_json << "\"c_over_a\":" << data(3)/data(1) << eendl;
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
     
          str_aus=str;//,str_sp,str_sc; //CO171027
      str_aus.sym_eps = str_sp.sym_eps = str_sc.sym_eps = tolerance; 
      if(smode=="EDATA") {
          if(!already_calculated){  //CO 171025
	str_aus.GetLatticeType(str_sp,str_sc);
        if(orig_tolerance == AUROSTD_NAN){
          orig_tolerance = str_sp.sym_eps;  // Use new tolerance calculated here 
        }
        tolerance = str_sp.sym_eps;
        sym_eps_change_count = str_aus.sym_eps_change_count = str_sp.sym_eps_change_count; //DX 2/26/18 - added sym eps change count
          }
       
        // Real space: bravais lattice primitive
        if(str_aus.bravais_lattice_type.size()){
          sscontent_json << "\"Bravais_lattice_type\":\"" << str_aus.bravais_lattice_type << "\"" << eendl;
        } else{
          sscontent_json << "\"Bravais_lattice_type\":null" << eendl;
        }
        vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

        // Real space: bravais lattice variation
        if(str_aus.bravais_lattice_variation_type.size()){
          sscontent_json << "\"Bravais_lattice_variation_type\":\"" << str_aus.bravais_lattice_variation_type << "\"" << eendl;
        } else{
          sscontent_json << "\"Bravais_lattice_variation_type\":null" << eendl;
        }
        vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
        
        // Real space: lattice system
        if(str_aus.bravais_lattice_system.size()){
          sscontent_json << "\"Bravais_lattice_system\":\"" << str_aus.bravais_lattice_system << "\"" << eendl;
        } else{
          sscontent_json << "\"Bravais_lattice_system\":null" << eendl;
        }
        vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

        // Real space: Pearson symbol
        if(str_aus.pearson_symbol.size()){
          sscontent_json << "\"Pearson_symbol\":\"" << str_aus.pearson_symbol << "\"" << eendl;
        } else{
          sscontent_json << "\"Pearson_symbol\":null" << eendl;
        }
        vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
        
        // Real space: crystal family
        if(str_aus.crystal_family.size()){
          sscontent_json << "\"crystal_family\":\"" << str_aus.crystal_family << "\"" << eendl;
        } else{
          sscontent_json << "\"crystal_family\":null" << eendl;
        }
        vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
        
        // Real space: crystal system
        if(str_aus.crystal_system.size()){
          sscontent_json << "\"crystal_system\":\"" << str_aus.crystal_system << "\"" << eendl;
        } else{
          sscontent_json << "\"crystal_system\":null" << eendl;
        }
        vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
 
        // Real space: point group crystal class
        if(str_aus.point_group_crystal_class.size()){
          sscontent_json << "\"point_group_crystal_class\":\"" << str_aus.point_group_crystal_class << "\"" << eendl;
        } else{
          sscontent_json << "\"point_group_crystal_class\":null" << eendl;
        }
        vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

        // Real space: point group Hermann Mauguin
        if(str_aus.point_group_Hermann_Mauguin.size()){
          sscontent_json << "\"point_group_Hermann_Mauguin\":\"" << str_aus.point_group_Hermann_Mauguin << "\"" << eendl;
        } else{
          sscontent_json << "\"point_group_Hermann_Mauguin\":null" << eendl;
        }
        vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
        
        // Real space: point group Shoenflies
        if(str_aus.point_group_Shoenflies.size()){
          sscontent_json << "\"point_group_Shoenflies\":\"" << str_aus.point_group_Shoenflies << "\"" << eendl;
        } else{
          sscontent_json << "\"point_group_Shoenflies\":null" << eendl;
        }
        vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

        // Real space: point group orbifold
        if(str_aus.point_group_orbifold.size()){
          sscontent_json << "\"point_group_orbifold\":\"" << str_aus.point_group_orbifold << "\"" << eendl;
        } else{
          sscontent_json << "\"point_group_orbifold\":null" << eendl;
        }
        vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
        
        // Real space: point group type
        if(str_aus.point_group_type.size()){
          sscontent_json << "\"point_group_type\":\"" << str_aus.point_group_type << "\"" << eendl;
        } else{
          sscontent_json << "\"point_group_type\":null" << eendl;
        }
        vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

        // Real space: point group order
        if(str_aus.point_group_order.size()){
          sscontent_json << "\"point_group_order\":\"" << str_aus.point_group_order << "\"" << eendl;
        } else{
          sscontent_json << "\"point_group_order\":null" << eendl;
        }
        vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

        // Real space: point group structure
        if(str_aus.point_group_structure.size()){
          sscontent_json << "\"point_group_structure\":\"" << str_aus.point_group_structure << "\"" << eendl;
        } else{
          sscontent_json << "\"point_group_structure\":null" << eendl;
        }
        vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

        // SGDATA
            PrintSGData(str_aus,tolerance,sscontent_json,no_scan,false,format,already_calculated);  //DX 8/31/17 - SGDATA //CO171027
        sym_eps_change_count = str_aus.sym_eps_change_count; //DX 2/26/18 - added sym eps change count
        vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
        //DX 9/1/17 - Add consistency check for input symmetry method and ITC method - START
        int multiplicity_of_primitive =str_sp.fgroup.size()/str_sp.pgroup_xtal.size();
        bool derivative_structure=false;
        if(SYM::ComparePointGroupAndSpaceGroupString(str_aus,multiplicity_of_primitive,derivative_structure) || no_scan == true){
          symmetry_commensurate=true;
            }else{ //CO171027
          if(LDEBUG){
            cerr << "pflow::PrintData: WARNING: Space group symbol and point group symbol do not match." << endl;
          }
          str_aus.ClearSymmetry();
          vcontent_json.clear();
          //if(!SYM::change_tolerance(str_aus,tolerance,orig_tolerance,change_sym_count,str_aus.dist_nn_min,no_scan)){
          if(!SYM::change_tolerance(str_aus,tolerance,str_aus.dist_nn_min,no_scan)){
            if(force_perform){
              if(LDEBUG){
                cerr << "pflow::PrintData: WARNING: Scan failed. Reverting back to original tolerance and recalculating as is (with aforementioned inconsistencies)." << str.directory << endl;
              }
              no_scan=true; //Force it to continue 
              tolerance = orig_tolerance;
              str_aus.ClearSymmetry();
              sscontent_json.str("");
              vcontent_json.clear();
              continue;
            }
          }
          else{
            str_aus.ClearSymmetry();
            sscontent_json.str("");
            vcontent_json.clear();
            sym_eps_change_count = str_aus.sym_eps_change_count; //DX 2/26/18 - added sym eps change count
            continue;
          }
        }
        //DX 9/1/17 - Add consistency check for input symmetry method and ITC method - END
        
        // Real space: bravais lattice lattice type
        if(str_aus.bravais_lattice_lattice_type.size()){
          sscontent_json << "\"Bravais_lattice_lattice_type\":\"" << str_aus.bravais_lattice_lattice_type << "\"" << eendl;
        } else{
          sscontent_json << "\"Bravais_lattice_lattice_type\":null" << eendl;
        }
        vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

        // Real space: bravais lattice lattice variation type
        if(str_aus.bravais_lattice_lattice_variation_type.size()){
          sscontent_json << "\"Bravais_lattice_lattice_variation_type\":\"" << str_aus.bravais_lattice_lattice_variation_type << "\"" << eendl;
        } else{
          sscontent_json << "\"Bravais_lattice_lattice_variation_type\":null" << eendl;
        }
        vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

        // Real space: bravais lattice lattice system
        if(str_aus.bravais_lattice_lattice_system.size()){
          sscontent_json << "\"Bravais_lattice_lattice_system\":\"" << str_aus.bravais_lattice_lattice_system << "\"" << eendl;
        } else{
          sscontent_json << "\"Bravais_lattice_lattice_system\":null" << eendl;
        }
        vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

        // Real space: bravais superlattice type
        if(str_aus.bravais_superlattice_type.size()){
          sscontent_json << "\"Bravais_superlattice_type\":\"" << str_aus.bravais_superlattice_type << "\"" << eendl;
        } else{
          sscontent_json << "\"Bravais_superlattice_type\":null" << eendl;
        }
        vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

        // Real space: bravais superlattice variation type
        if(str_aus.bravais_superlattice_variation_type.size()){
          sscontent_json << "\"Bravais_superlattice_variation_type\":\"" << str_aus.bravais_superlattice_variation_type << "\"" << eendl;
        } else{
          sscontent_json << "\"Bravais_superlattice_variation_type\":null" << eendl;
        }
        vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

        // Real space: bravais superlattice system
        if(str_aus.bravais_superlattice_system.size()){
          sscontent_json << "\"Bravais_superlattice_system\":\"" << str_aus.bravais_superlattice_system << "\"" << eendl;
        } else{
          sscontent_json << "\"Bravais_superlattice_system\":null" << eendl;
        }
        vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

        // Real space: pearson symbol superlattice
        if(str_aus.pearson_symbol_superlattice.size()){
          sscontent_json << "\"Pearson_symbol_superlattice\":\"" << str_aus.pearson_symbol_superlattice << "\"" << eendl;
        } else{
          sscontent_json << "\"Pearson_symbol_superlattice\":null" << eendl;
        }
        vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      }
      //RECIP
      // Reciprocal space lattice
        if(str_aus.klattice.rows){
          sscontent_json << "\"reciprocal_lattice_vectors\":[" << aurostd::xmatDouble2String(str_aus.klattice,5,roff) << "]" << eendl;
        } else{
          sscontent_json << "\"reciprocal_lattice_vectors\":null" << eendl;
        }
        vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      
      // Reciprocal lattice parameters
      data=Getabc_angles(str.klattice,DEGREES);
      if(data.rows){
        sscontent_json << "\"reciprocal_lattice_parameters\":[" << aurostd::joinWDelimiter(xvecDouble2vecString(data,5,roff),",") << "]" << eendl;
      } else{
        sscontent_json << "\"reciprocal_lattice_parameters\":null" << eendl;
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

      // Reciprocal space volume
      double kvol=GetVol(str.klattice);
      sscontent_json << "\"reciprocal_volume\":" << kvol << eendl;
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      if(smode=="EDATA") {
        // Reciprocal space: reciprocal lattice type
        if(str_aus.reciprocal_lattice_type.size()){
          sscontent_json << "\"reciprocal_lattice_type\":\"" << str_aus.reciprocal_lattice_type << "\"" << eendl;
        } else{
          sscontent_json << "\"reciprocal_lattice_type\":null" << eendl;
        }
        vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
        
        // Reciprocal space: reciprocal lattice variation type
        if(str_aus.reciprocal_lattice_variation_type.size()){
          sscontent_json << "\"reciprocal_lattice_variation_type\":\"" << str_aus.reciprocal_lattice_variation_type << "\"" << eendl;
        } else{
          sscontent_json << "\"reciprocal_lattice_variation_type\":null" << eendl;
        }
        vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

        // SPRIM
        sscontent_json << "\"standard_primitive_structure\":" << xstructure2json(str_sp) << eendl;
        vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
        
        // SCONV
        sscontent_json << "\"standard_conventional_structure\":" << xstructure2json(str_sc) << eendl;
        vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      }
      oss << "{" << aurostd::joinWDelimiter(vcontent_json,",")  << "}" << endl;
    }
    }
    oss_final << oss.str();
        str_sym=str_aus; //CO171027
  }
} // namespace pflow

namespace pflow {
      void PrintData(const xstructure& str,ostream& oss,string smode, const string& format) {
        xstructure str_sym,str_sp,str_sc; //CO171027
    double tolerance = AUROSTD_NAN;
    bool no_scan = false;
        PrintData(str,str_sym,str_sp,str_sc,oss,smode,tolerance,no_scan,format); //CO171027
  }
} // namespace pflow

namespace pflow {
      void PrintData(const xstructure& str,ostream& oss,string smode, double tolerance, bool no_scan, const string& format) {
        xstructure str_sym,str_sp,str_sc; //CO171027
        PrintData(str,str_sym,str_sp,str_sc,oss,smode,tolerance,no_scan,format); //CO171027
  }
} // namespace pflow

// ***************************************************************************
// pflow::PrintData1
// ***************************************************************************
// Prints out the comparison between two xstructures.
// Dane Morgan - Stefano Curtarolo
namespace pflow {
  void PrintData1(const xstructure& str1, const double& rcut, ostream& oss) {
    oss << PrintData1(str1,rcut);
  }
} // namespace pflow

namespace pflow {
  string PrintData1(const xstructure& str1, const double& rcut) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "pflow::PrintData1: BEGIN" << endl;
    // LDEBUG=TRUE;
    stringstream oss;
    oss.setf(std::ios::left, std::ios::adjustfield);
    oss.setf(std::ios::fixed, std::ios::floatfield);
    oss.precision(4);
    double cutoff=rcut;
    //  oss.width(10);
    // Rescale to avoid scale problems
    xstructure sstr1=str1;
    sstr1=ReScale(sstr1,1.0);
    int nat1=sstr1.atoms.size();
    int ntyp1=sstr1.num_each_type.size();
    std::deque<int> num_each_type1=sstr1.num_each_type;

    if(LDEBUG) cerr << "PrintData1 [2]" << endl;
    oss << "*************************  Structure Comparison *************************" << endl;

    // Nun atoms
    oss << "Number of atoms" << endl;
    oss << "  Num Atom Str1: " << nat1 << endl;

    // Num types
    oss << "Number of types" << endl;
    oss << "  Num Type Str1: " << ntyp1 << endl;

    // Stoichiometry
    int wid=8;
    vector<double> tmp_stoich1(ntyp1,0);
    oss << "Stoichiometry" << endl;
    oss << "  Type:        ";
    for(int it=0;it<ntyp1;it++) {
      oss << setw(wid) << it << " ";
    }
    oss << endl;
    oss << "  Stoich Str1: ";
    for(int it=0;it<ntyp1;it++) {
      tmp_stoich1[it]=(double) num_each_type1[it]/ (double) nat1;
      oss << setw(wid) << tmp_stoich1[it] << " ";
    }
    oss << endl;

    if(LDEBUG) cerr << "PrintData1 [3]" << endl;
    // Volume
    double v1=sstr1.scale*sstr1.scale*sstr1.scale*GetVol(sstr1.lattice)/nat1;
    oss << "Volume per atom" << endl;
    oss << "  Vol Str1: " << v1 << endl;

    // lat par and angles
    //  vector<double> ldat1=pflow::Getabc_angles(pflow::xmatrix2matrix(sstr1.lattice));
    vector<double> ldat1=xvector2vector(Getabc_angles(sstr1.lattice,DEGREES));
    ldat1=pflow::Sort_abc_angles(ldat1);  
    oss << "Unit cell: tot_vol a b c alpha beta gamma" << endl;
    oss << "  Cell Str1: ";
    oss << v1*nat1 << " ";
    pflow::Vout(ldat1,oss);

    // Volume adjusted distances

    double rescale1=std::pow((double) 1.0/v1,(double) 1.0/3.0);
    xstructure str1_newvol=sstr1;
    // str1_newvol=SetScale(str1_newvol,rescale1); // Sets vol/atom=1
    str1_newvol.scale=rescale1;
    pflow::matrix<double> dist1;
    pflow::matrix<double> dist2;
    pflow::matrix<double> dist_diff;
    pflow::matrix<double> dist_diff_n;
    // If cutoff < 0 use default of 4*(avg rad of an atom).  Note
    // that in the rescaled str, where each atom has vavg=1, rad=.6204
    // In the true str, where each atom has vol = v1, rad=0.6204*v1^(1/3).
    if(LDEBUG) cerr << "PrintData1 [3]" << endl;
    if(rcut<0) cutoff=4*0.6204*std::pow((double) v1,(double) 1.0/3.0);  // THAT`S WRONG
    //  cerr << cutoff << endl;
    if(rcut<0) cutoff=1.01*RadiusSphereLattice(str1_newvol.lattice); // better way to nail it for sure I`ll hit something. 1% to be sure of roundoff
    //  cerr << cutoff << endl;
    //  if(rcut<0) cutoff=max(modulus(str1_newvol.lattice(1)),modulus(str1_newvol.lattice(2)),modulus(str1_newvol.lattice(3))); // better way to nail it for sure I`ll hit something.
    //  cerr << cutoff << endl;
    if(LDEBUG) {cerr << "PrintData1 [4]  " << 4*0.6204*std::pow((double) v1,(double) 1.0/3.0)<< " " << RadiusSphereLattice(str1_newvol.lattice) << " " << max(modulus(str1_newvol.lattice(1)),modulus(str1_newvol.lattice(2)),modulus(str1_newvol.lattice(3))) << endl;}// exit(0);}
    pflow::CmpStrDist(str1_newvol,str1_newvol,cutoff*rescale1,dist1,dist2,dist_diff,dist_diff_n);
    int nprs1=dist1.size();
    oss << "Avg of Distance Magnitude Differences" << endl;
    oss << "  Rescaling vol/atom of Str1 to the average vol/atom (return to orig scale factor): "
        << str1_newvol.scale*str1_newvol.scale*str1_newvol.scale*GetVol(str1_newvol.lattice)/nat1 << " (" << 1/rescale1 << ")" << endl;
    oss << "  Working with the all the types up to num types (num pair types): " << ntyp1 << " (" << nprs1 << ")" << endl;
    oss << "  Including all pairs within cutoff: " << cutoff << endl;

    // double tadiff=0;   // DANE not used
    // double tadiff_n=0; // DANE not used
    // int tnn=0;  // DANE not used
    for(int it1=0;it1<ntyp1;it1++) {
      for(int it2=it1;it2<ntyp1;it2++) {
        int id=it2-it1+ntyp1*it1-std::max(0,it1*(it1-1)/2);
        int nn=dist1[id].size();      
        oss << "  P" << id << "  " << nn << " Pairs between types: "
            << it1 << " " << it2 << "   ";
        for(int in=0;in<nn;in++) {
          oss << " " << dist1[id][in];
        }
        oss << endl;
      }// it2
    }// it1

    //    // Number of each type up to rcut.
    //    nprs_max=ntyp_max*(ntyp_max+1)/2;
    //    vector<int> n1(nprs_max,0);
    //    vector<int> n2(nprs_max,0);
    //    oss << "Number neighbors of each type within cutoff for rescaled volume" << endl;
    //    oss << "  Cutoff: " << cutoff << endl;
    //    oss << "  Rescaling vol/atom of Str1 and Str2 to the average vol/atom: " << vavg << endl;
    //    oss << "           ";
    //    for(int it1=0;it1<ntyp_max;it1++) {
    //      for(int it2=it1;it2<ntyp_max;it2++) {
    //        oss << it1 << "-" << it2 << "  ";
    //      }
    //    }
    //    oss << endl;
    //    oss << "  NN Str1  ";
    //    for(int it1=0;it1<ntyp_max;it1++) {
    //      for(int it2=it1;it2<ntyp_max;it2++) {
    //        if(it1<ntyp1 && it2<ntyp1) {
    //  	int id1=it2-it1+ntyp1*it1-std::max(0,it1*(it1-1)/2);
    //  	n1[id1]=dist1[id1].size();
    //  	oss << setw(4) << n1[id1] << " ";
    //        }
    //      }
    //    }
    //    oss << endl;
    //    oss << "  NN Str2  ";
    //    for(int it1=0;it1<ntyp_max;it1++) {
    //      for(int it2=it1;it2<ntyp_max;it2++) {
    //        if(it1<ntyp2 && it2<ntyp2) {
    //  	int id2=it2-it1+ntyp2*it1-std::max(0,it1*(it1-1)/2);
    //  	n2[id2]=dist2[id2].size();
    //  	oss << setw(4) << n2[id2] << " ";
    //        }
    //      }
    //    }
    //    oss << endl;  
    //    oss << "  NN Err   ";
    //    for(int it1=0;it1<ntyp_max;it1++) {
    //      for(int it2=it1;it2<ntyp_max;it2++) {
    //        int id=it2-it1+ntyp_max*it1-std::max(0,it1*(it1-1)/2);
    //        oss << setw(4) << n2[id]-n1[id] << " ";
    //      }
    //    }
    //    oss << endl;  

    //    // Space group
    //    bool Platon_EQUAL=FALSE,Platon_EXACT=FALSE;
    //    double Platon_ang=1.0,Platon_d1=0.25,Platon_d2=0.25,Platon_d3=0.25;
    //    system("touch pflow_tmp_in"); // Create input file to connect to stream.
    //    ofstream outfile("pflow_tmp_out");
    //    ifstream infile("pflow_tmp_in");
    //    string spcgrp1a,spcgrp2a,spcgrp1b,spcgrp2b;
    //    platon2print(str1,Platon_EQUAL,Platon_EXACT,Platon_ang,Platon_d1,Platon_d2,Platon_d3,outfile);
    //    system("cat pflow_tmp_out | platonSG > pflow_tmp_in");
    //    infile >> spcgrp1a >> spcgrp1b;
    //    // Puts file position pointers back to beginning.
    //    outfile.seekp(0);
    //    infile.seekg(0);
    //    platon2print(sstr2,Platon_EQUAL,Platon_EXACT,Platon_ang,Platon_d1,Platon_d2,Platon_d3,outfile);
    //    system("cat pflow_tmp_out | platonSG > pflow_tmp_in");
    //    infile >> spcgrp2a >> spcgrp2b;
    //    system("/bin/rm pflow_tmp_out");
    //    system("/bin/rm pflow_tmp_in");
    //    oss << "Space groups " << endl;
    //    oss << "  SG Str1: " << spcgrp1a << " " << spcgrp1b << endl;
    //    oss << "  SG Str2: " << spcgrp2a << " " << spcgrp2b << endl;

    // Center of mass
    return oss.str();
  } // end routine
} // namespace pflow

// ***************************************************************************
// pflow::PrintData2
// ***************************************************************************
// This funtion prints out structural data.
// Stefano Curtarolo
namespace pflow {
  void PrintData2(const xstructure& str, ostream& oss) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "pflow::PrintData2: BEGIN" << endl;
    oss.setf(std::ios::fixed,std::ios::floatfield);
    oss.precision(6);
    xvector<double> data(6);
    oss << endl;
    oss << "DIRECT LATTICE" << endl;
    data=Getabc_angles(str.lattice,DEGREES);
    oss << "volume = " << str.scale*str.scale*str.scale*GetVol(str.lattice) << endl;
    oss << "a,b,c,alpha,beta,gamma = "
        << str.scale*data(1) << " " << str.scale*data(2) << " " << str.scale*data(3) << " "
        << data(4) << " " << data(5) << " " << data(6) << endl;
    oss << "b/a = " << data(2)/data(1) << "  c/a = " << data(3)/data(1) << endl;
    oss << "Lattice =" << endl;
    oss << str.scale*str.lattice << endl;
    oss << endl;      
    oss << "RECIPROCAL LATTICE" << endl;
    data=Getabc_angles(str.klattice,DEGREES);
    oss << "volume = " << GetVol(str.klattice) << endl;
    oss << "a,b,c,alpha,beta,gamma = "
        << str.scale*data(1) << " " << str.scale*data(2) << " " << str.scale*data(3) << " "
        << data(4) << " " << data(5) << " " << data(6) << endl;
    oss << "K-Lattice =" << endl;
    oss << str.klattice << endl;
    oss << endl;
  }
} // namespace pflow

// ***************************************************************************
// pflow::PrintDisplacements
// ***************************************************************************
//  This funtion prints out displacements
// Dane Morgan - Stefano Curtarolo
namespace pflow {
  void PrintDisplacements(xstructure str, const double cutoff, ostream& oss) {  
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "pflow::PrintDisplacements: BEGIN" << endl;
    oss.setf(std::ios::left, std::ios::adjustfield);
    oss.setf(std::ios::fixed, std::ios::floatfield);
    //  oss.setf(std::ios::fixed,std::ios::floatfield);
    //  oss.precision(10);
    deque<deque<_atom> > neigh_mat;
    // [OBSOLETE]  pflow::GetStrNeighData(str,cutoff,neigh_mat);
    str.GetStrNeighData(cutoff,neigh_mat);
    // double tol=1e-15;   // DANE not used
    xmatrix<double> lattice(3);
    lattice=str.lattice;
    int coord_flag=str.coord_flag;

    // Output header
    oss << "All neighbor data given as" << endl;
    oss << "[basis number]  [name]  [n1 n2 n3 (unit cell)] [d1 d2 d3 (displacement)] [distance]" << endl;
    oss << "Distances are only calculated for the neighbors within " << cutoff <<  endl;

    // Print out dislpacements and distances
    oss << endl;
    oss << "********************" << " DISPLACEMENTS " << "********************" << endl;
    for(uint ia=0;ia<neigh_mat.size();ia++) {
      _atom a = neigh_mat.at(ia).at(0);

      // Output reference atom info.
      oss << setw(4) << a.number+1 << " " << setw(4) << a.name.c_str();
      xvector<double> pos(3);
      if(str.coord_flag==FALSE) { // direct
        pos=a.fpos;
      } else{
        pos=a.cpos;
      }
      oss << "   " << setw(9) << setprecision(4) << pos(1)
          << setw(9) << setprecision(4) << pos(2)
          << setw(9) << setprecision(4) << pos(3)
          << endl;

      for(uint in=1;in<neigh_mat.at(ia).size();in++) {
        _atom an = neigh_mat.at(ia).at(in);
        xvector<int> ijk(3);
        ijk=an.ijk;
        oss << "        ";
        oss << setw(4) << an.number+1 << " ";
        oss << setw(4) << an.name.c_str() << "   ";
        oss << setw(3) << ijk(1) << " " << setw(3) << ijk(2) << " " << setw(3) << ijk(3) << "   ";
        xvector<double> disp(3);
        disp=AtomCDisp(a,an);
        if(coord_flag==0) { // direct
          disp=C2F(lattice,disp);
        }
        oss << setw(9) << setprecision(4) << disp(1) << " " << setw(9) << setprecision(4) << disp(2) << " " << setw(9) << setprecision(4) << disp(3);
        oss << "   " << setw(9) << setprecision(4) << AtomDist(a,an);
        oss << endl;
      } // in
    } // ia

  } // end routine
} // namespace pflow

// ***************************************************************************
// pflow::PrintDistances
// ***************************************************************************
// This funtion prints out distances
// Dane Morgan - Stefano Curtarolo
namespace pflow {
  void PrintDistances(xstructure str, const double cutoff, ostream& oss) {  
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "pflow::PrintDistances: BEGIN" << endl;
    oss.setf(std::ios::left, std::ios::adjustfield);
    oss.setf(std::ios::fixed, std::ios::floatfield);
    //  oss.setf(std::ios::fixed,std::ios::floatfield);
    oss.precision(10);
    deque<deque<_atom> > neigh_mat;
    // [OBSOLETE] pflow::GetStrNeighData(str,cutoff,neigh_mat);
    str.GetStrNeighData(cutoff,neigh_mat);
    //  cerr << "DEBUG: neigh_mat.size()=" << neigh_mat.size() << endl;
    //  for(uint i=0;i<neigh_mat.size();i++)
    //   cerr << "DEBUG: neigh_mat.at(i).size()=" << neigh_mat.at(i).size() << endl;

    // double tol=1e-15;  // DANE not used
    uint xtra=0;//3;

    // Output header
    oss << "All atoms identified by" << endl;
    oss << "[basis number]  [name]  [n1 n2 n3 (unit cell)]" << endl;
    oss << "Distances are only calculated for the neighbors within " << cutoff <<  endl;

    // Print out distances
    oss << endl;
    oss << "********************" << " DISTANCES " << "********************" << endl;
    for(uint ia=0;ia<neigh_mat.size();ia++) {
      _atom a = neigh_mat.at(ia).at(0);

      // Output reference atom info.
      oss << setw(4) << a.number+1 << " " << setw(4) << a.name.c_str();
      xvector<double> pos(3);
      if(str.coord_flag==FALSE) { // direct
        pos=a.fpos;
      } else{
        pos=a.cpos;
      }
      oss << "   "
          << setw(9+xtra) << setprecision(4+xtra) << pos(1)
          << setw(9+xtra) << setprecision(4+xtra) << pos(2)
          << setw(9+xtra) << setprecision(4+xtra) << pos(3)
          << endl;

      for(uint in=1;in<neigh_mat.at(ia).size();in++) {
        xvector<int> ijk(3);ijk=neigh_mat.at(ia).at(in).ijk;
        xvector<double> fpos(3);fpos=neigh_mat.at(ia).at(in).fpos;
        xvector<double> cpos(3);cpos=neigh_mat.at(ia).at(in).cpos;
        oss << "      ";
        oss << setw(4) << neigh_mat.at(ia).at(in).number+1 << " ";
        oss << setw(4) << neigh_mat.at(ia).at(in).name.c_str() << "   ";
        oss << setw(3+xtra) <<  ijk(1) << " " << setw(3+xtra) <<  ijk(2) << " " << setw(3+xtra) <<  ijk(3) << "   ";
        //     oss << setw(3+xtra) << "F " << fpos(1) << " " << setw(3+xtra) << fpos(2) << " " << setw(3+xtra) << fpos(3) << "   ";
        //   oss << setw(3+xtra) << "C " << cpos(1) << " " << setw(3+xtra) << cpos(2) << " " << setw(3+xtra) << cpos(3) << "   ";
        oss << setprecision(4+xtra) << AtomDist(a,neigh_mat.at(ia).at(in));
        oss << endl;
      } // in
    } // ia

  } // end routine
} // namespace pflow

// ***************************************************************************
// pflow::PrintEwald
// ***************************************************************************
// Dane Morgan - Stefano Curtarolo
namespace pflow {
  void PrintEwald(const xstructure& in_str, double& epoint,
                  double& ereal, double& erecip, double& eewald,
                  double& eta, const double& SUMTOL, ostream& oss) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "pflow::PrintEwald: BEGIN" << endl;
    oss.setf(std::ios::left, std::ios::adjustfield);
    oss.setf(std::ios::fixed, std::ios::floatfield);
    pflow::matrix<double> lat=pflow::GetScaledLat(in_str);
    int nat = pflow::GetNumAtoms(in_str);
    double vol=pflow::GetVol(lat);
    double d=std::pow((double) vol,(double) 1.0/3.0)*1E-10; // In meters
    double eV2NoUnit=2*TWOPI*EPS_VACUUM*d/(E_ELECTRON*nat);
    int wid;
    int prec;
    wid=15;
    prec=6;
    oss << "*************************** Ewald Results ***********************" << endl;
    oss << "E(NoUnit)=PerAtom(eV)*4*pi*eps_0*vol^1/3/e_electron=Charge/Distance unit independent."<< endl;
    oss << "Using eta= " << setprecision(2*prec) << eta << endl;
    oss << "All real and k-space sums stop when maximum term in shell is less than SUMTOL= " << std::setiosflags(std::ios::scientific) << SUMTOL << endl;
    oss << "SUMTOL is the following percentage of your total Ewald energy: " << std::setiosflags(std::ios::scientific) << 100*SUMTOL/eewald << endl;
    oss.setf(std::ios::fixed, std::ios::floatfield);
    oss << "                         "
        << setw(wid) << "Total(eV)"
        << setw(wid) << "Per atom(eV)"
        << setw(wid) << "Total(NoUnit)" << endl;
    oss << "Point Energy:            "
        << " " << setw(wid) << setprecision(prec) << epoint
        << " " << setw(wid) << setprecision(prec) << epoint/nat
        << " " << setw(wid) << setprecision(prec) << eV2NoUnit*epoint
        << endl;
    oss << "Real Space Energy:       "
        << " " << setw(wid) << setprecision(prec) << ereal
        << " " << setw(wid) << setprecision(prec) << ereal/nat
        << " " << setw(wid) << setprecision(prec) << eV2NoUnit*ereal
        << endl;
    oss << "Recipricol Space Energy: "
        << " " << setw(wid) << setprecision(prec) << erecip
        << " " << setw(wid) << setprecision(prec) << erecip/nat
        << " " << setw(wid) << setprecision(prec) << eV2NoUnit*erecip
        << endl;
    oss << "Total Ewald Energy:      "
        << " " << setw(wid) << setprecision(prec) << eewald
        << " " << setw(wid) << setprecision(prec) << eewald/nat
        << " " << setw(wid) << setprecision(prec) << eV2NoUnit*eewald
        << endl;
  }
} // namespace pflow

// ***************************************************************************
// pflow::PrintGulp
// ***************************************************************************
// This funtion prints out structural data in a format for
// gulp (file will cause gulp to calculate distances between
// atoms when used as input).
// Stefano Curtarolo
namespace pflow {
  void PrintGulp(const xstructure& str, ostream& oss) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "pflow::PrintGulp: BEGIN" << endl;
    oss.setf(std::ios::fixed,std::ios::floatfield);
    oss.precision(10);
    xstructure sstr=str;
    // Set scale to 1 so you don't need to rescale coordinates.
    sstr=ReScale(sstr,1.0);

    // Print out data in Gulp format
    oss << "distance" << endl;
    oss << "cutd 5" << endl;
    oss << "vectors" << endl;
    oss << sstr.lattice(1,1) << " " << sstr.lattice(1,2) << " " << sstr.lattice(1,3) << "    Bravais lattice" << endl;
    oss << sstr.lattice(2,1) << " " << sstr.lattice(2,2) << " " << sstr.lattice(2,3) << "    Bravais lattice" << endl;
    oss << sstr.lattice(3,1) << " " << sstr.lattice(3,2) << " " << sstr.lattice(3,3) << "    Bravais lattice" << endl;

    for(uint i=0;i<sstr.num_each_type.size();i++)
      oss <<  sstr.num_each_type.at(i) << " ";

    oss << endl;
    if(sstr.coord_flag==0) oss << "Fractional" << endl;
    if(sstr.coord_flag==1) oss << "Cartesian" << endl;
    for(uint i=0;i<sstr.atoms.size();i++) {
      if(sstr.coord_flag==0)
        oss << " " << sstr.atoms.at(i).fpos(1) << " " << sstr.atoms.at(i).fpos(2) << " " << sstr.atoms.at(i).fpos(3) << endl;
      if(sstr.coord_flag==1) oss << " " << sstr.atoms.at(i).cpos(1) << " " << sstr.atoms.at(i).cpos(2) << " " << sstr.atoms.at(i).cpos(3) << endl;
    }
  }
} // namespace pflow

// ***************************************************************************
// PrintKmesh
// ***************************************************************************
//  This funtion prints a set of Monkhurst-Pack kpoint mesh values
// (divisions along each reciprical lattice parameter)
// based on the desired kpt density and the lattice parameters.  The
// values assure the most even distribution of kpts along the lattice
// params consistent with the kpt density.  Returns 4 kmesh vectors of length
// 3
//    1. Rounded to Nint
//    2. Exact
//    3. Rounded down
//    4. Rounded up
// Dane Morgan - Stefano Curtarolo
void PrintKmesh(const xmatrix<double>& kmesh, ostream& oss) {
  oss.setf(std::ios::fixed,std::ios::floatfield);
  oss.precision(3);
  vector<string> s(4);
  s[0]="RoundNint:";
  s[1]="Exact:    ";
  s[2]="RoundDown:";
  s[3]="RoundUp:  ";
  for(int i=1;i<=4;i++)
    oss << s[i-1] << "  " << kmesh(i,1) << "  " << kmesh(i,2) << "  " << kmesh(i,3) << endl;
}

// ***************************************************************************
// PrintImages
// ***************************************************************************
//  This function prints out all the image files.
// Dane Morgan - Stefano Curtarolo
void PrintImages(xstructure strA, xstructure strB, const int& ni, const string& path_flag) {  
  // Get cm of strA. Get image at i/(n-1) between
  // the A and B.  Set mom1 of image to be that of strA.  
  // Output that image as image i.
  xvector<double> mom1A=GetMom1(strA);
  for(int im=0;im<ni+2;im++) {
    double f;
    if((ni+2)>1) {
      f=(double)im/double(ni+1);
    } else{
      f=0;
    }
    xstructure imstr=GetIntpolStr(strA,strB,f,path_flag);
    imstr=SetMom1(imstr,mom1A);
    // Set title to include image name.
    ostringstream titlestream;
    if(im<10) {
      titlestream << "POSCAR_0" << im << ":" << imstr.title << ends;
    } else{
      titlestream << "POSCAR_" << im << ":" << imstr.title << ends;
    }
    imstr.title=titlestream.str();
    // Set output file name.
    ostringstream namestream;
    if(im<10) {
      namestream << "POSCAR_0" << im << ends;
    } else{
      namestream << "POSCAR_" << im << ends;
    }
    //  char* namechar = namestream.str().c_str();
    //  ofstream outfile(namechar);
    ofstream outfile(namestream.str().c_str());
    outfile << imstr;
    // Create directory
    ostringstream aus;
    if(im<10) {
      aus << "mkdir 0" << im << ends;
      aurostd::execute(aus);
      aus << "cp " << namestream.str().c_str() << " " << "./" << "0" << im << "/POSCAR" << endl;
      aurostd::execute(aus);
    } else{
      aus << "mkdir " << im << ends;
      aurostd::execute(aus);
      aus << "cp " << namestream.str().c_str() << " " << "./" << im << "/POSCAR" << endl;
      aurostd::execute(aus);
    }
  }
} // end routine

// ***************************************************************************
// PrintMSI
// ***************************************************************************
//  This funtion prints out structural data in a msi format (for cerius).
// Dane Morgan - Stefano Curtarolo
void PrintMSI(const xstructure& str, ostream& oss) {
  oss.setf(std::ios::fixed,std::ios::floatfield);
  oss.precision(10);
  xstructure sstr=str;

  // Set scale to 1 so you don't need to rescale coordinates.
  sstr=ReScale(sstr,1.0);

  // Print out data in msi format

  // Tolerance
  double Tol=1e-10;

  // First we must determine the atomic numbers.
  for(uint i=0;i<sstr.atoms.size();i++) {
    sstr.atoms.at(i).CleanName();
    if(sstr.atoms.at(i).atomic_number<1) {
      cerr << "ERROR - PrintMSI: atomic_number not found: sstr.atoms.at(" << i << ").cleanname=" << sstr.atoms.at(i).cleanname << endl;
      exit(0);
    }
  }

  //   Now we have to orient the lattice vectors for msi.
  //   The desired orientation is obtained by
  //   rotating c around Z into the XZ plane (R1a), then rotating c around
  //   Y to align it with Z (R1b), then rotating b around Z (leaving c fixed)
  //   until it is in the YZ plane (R2).  The total rotation needed is then
  //   Rtot=R2*R1b*R1a.

  xvector<double> v(3);
  xmatrix<double> R1a(3,3);  
  xmatrix<double> R1b(3,3);  
  xmatrix<double> R2(3,3);  
  xmatrix<double> Rtot(3,3);  
  xvector<double> a(3);a=sstr.lattice(1);//double an=modulus(a); // DANE not used
  xvector<double> b(3);b=sstr.lattice(2);//double bn=modulus(b); // DANE not used
  xvector<double> c(3);c=sstr.lattice(3);double cn=modulus(c);

  double cxyn=sqrt(cn*cn-c(3)*c(3));
  // double bxyn=sqrt(bn*bn-b(3)*b(3)); // DANE not used
  // double cosa=cos(a,b); // DANE not used
  // double cosb=cos(a,c); // DANE not used
  // double cosg=cos(b,c); // DANE not used

  // Build R1a.
  double cosp,sinp;
  if(cxyn>Tol && (abs(c(1))>Tol || abs(c(2))>Tol)) {
    cosp=c(1)/cxyn;
    sinp=c(2)/cxyn;
  } else{
    cosp=1.0;
    sinp=0.0;
  }
  clear(R1a);
  R1a(1,1)=cosp;
  R1a(1,2)=sinp;
  R1a(2,1)=-sinp;
  R1a(2,2)=cosp;
  R1a(3,3)=1.0;

  // tpx
  //  Mout(R1a,oss);
  //  oss << "cosp, sinp" << cosp << " " << sinp << endl;
  // Build R1b.
  double costh=c(3)/cn;
  double sinth=sqrt(1-costh*costh); // This is OK since sinth>0.
  R1b(1,1)=costh;
  R1b(1,3)=-sinth;
  R1b(3,1)=sinth;
  R1b(3,3)=costh;
  R1b(2,2)=1.0;
  xmatrix<double> RR(3,3);
  RR=R1b*R1a;
  // Build R2.  
  // Using rotated b vector, called bp, (bp=R1b*R1a*b) to define the angles.
  xvector<double> bp(3);
  bp=RR*b;
  double bpn=modulus(bp);
  double bpxyn=sqrt(bpn*bpn-bp(3)*bp(3));
  double cospp,sinpp;
  if(bpxyn>Tol && (b(1)>Tol || b(2)>Tol)) {
    cospp=bp(2)/bpxyn;
    sinpp=bp(1)/bpxyn;
  } else {
    cospp=1.0;
    sinpp=0.0;
  }
  clear(R2);
  R2(1,1)=cospp;
  R2(1,2)=-sinpp;
  R2(2,1)=sinpp;
  R2(2,2)=cospp;
  R2(3,3)=1.0;
  // Calculate final rotation xmatrix.
  // tpx
  //  for(int i=1;i<=3;i++)
  //  for(int j=1;j<=3;j++)
  //  oss << "R2 " << R2(i,j) << endl;
  Rtot=R2*RR;

  //tpx
  //Mout(Rtot,oss);
  //oss << "HERE " << cxyn << " " << cosa << " " << cosb << " " << cosg << " " << cosp << " " << sinp << endl;
  // Now rotate all the lattice vectors and the positions.
  xmatrix<double> c2lat(3,3);
  std::vector<xvector<double> > c2pos(sstr.atoms.size());

  v=Rtot*a;
  for(int i=1;i<=3;i++)
    c2lat(1,i)=v(i);
  v=Rtot*b;
  for(int i=1;i<=3;i++)
    c2lat(2,i)=v(i);
  v=Rtot*c;
  for(int i=1;i<=3;i++)
    c2lat(3,i)=v(i);

  for(uint i=0;i<sstr.atoms.size();i++) {
    v=Rtot*sstr.atoms.at(i).cpos;
    c2pos.at(i)(1)=v(1);
    c2pos.at(i)(2)=v(2);
    c2pos.at(i)(3)=v(3);
  }

  // Write output.
  // Title.
  oss << "# MSI CERIUS2 DataModel File Version 3 5" << endl;

  // Stuff.
  oss << "(1 Model" << endl;

  // Lattice vectors.
  oss << " (A D A3 (" << c2lat(1,1) << " " << c2lat(1,2) << " " << c2lat(1,3) << " " << "))" << endl;
  oss << " (A D B3 (" << c2lat(2,1) << " " << c2lat(2,2) << " " << c2lat(2,3) << " " << "))" << endl;
  oss << " (A D C3 (" << c2lat(3,1) << " " << c2lat(3,2) << " " << c2lat(3,3) << " " << "))" << endl;

  // Stuff.
  oss << " (A I Id 1)" << endl;
  oss << " (A C Label \"Model1\")" << endl;

  // The atoms.
  for(uint i=0;i<sstr.atoms.size();i++) {
    oss << " (" << i+2 << " Atom" << endl;
    oss << "  (A C ACL \"" << sstr.atoms.at(i).atomic_number << " " << sstr.atoms.at(i).cleanname << "\")" << endl;
    oss << "  (A D XYZ (" << c2pos.at(i)(1) << " " << c2pos.at(i)(2) << " " << c2pos.at(i)(3) << " " << "))" << endl;
    oss << "  (A I Id " << i+1 << ")" << endl;
    oss << " )" << endl;
  }
  // Final parenthesis
  oss << ")" << endl;

}

// ***************************************************************************
// PrintNData
// ***************************************************************************
// This funtion prints out structural data.
// STefano Curtarolo
void PrintNdata(const xstructure& str, ostream& oss) {
  oss.setf(std::ios::fixed,std::ios::floatfield);
  oss.precision(10);
  xstructure sstr=str;

  // Set scale to 1 so you don't need to rescale coordinates.
  sstr=ReScale(sstr,1.0);
  // Print out data in Platon format
  double volume=GetVol(sstr.lattice);
  double normalization=std::pow((double) volume,(double) 1.0/3.0);

  printf(" %6.5f %6.5f %6.5f %6.5f %6.5f %6.5f \n",
         modulus(sstr.lattice(1))/normalization,
         modulus(sstr.lattice(2))/normalization,
         modulus(sstr.lattice(3))/normalization,
         aurostd::angle(sstr.lattice(2,1),sstr.lattice(2,2),sstr.lattice(2,3),sstr.lattice(3,1),sstr.lattice(3,2),sstr.lattice(3,3)),
         aurostd::angle(sstr.lattice(1,1),sstr.lattice(1,2),sstr.lattice(1,3),sstr.lattice(3,1),sstr.lattice(3,2),sstr.lattice(3,3)),
         aurostd::angle(sstr.lattice(1,1),sstr.lattice(1,2),sstr.lattice(1,3),sstr.lattice(2,1),sstr.lattice(2,2),sstr.lattice(2,3)));//,volume);
}

// ***************************************************************************
// PrintNeatProj
// ***************************************************************************
// Dane Morgan - Stefano Curtarolo
void PrintNeatProj(pflow::projdata& pd, ostream& oss) {
  //  namespace pflow
  pd.Print(oss);
}

// ***************************************************************************
// PrintPDB
// ***************************************************************************
// This funtion prints out structural data in a PDB format.
// Note that I don't really know anything about PDB but it seems
// to be an annoyingly column formatted.  Therefore, all the widths
// below must tbe kept as is.  This means that if we have width W
// for variable X, and prec P, then X takes up P+1 for for the decimal
// part and decimal point, and we have only R=W-(P+1)-1=W-P-2 remaining digits
// (the -1 is because if X takes up all of W then you run into
// the previous field).  So we have the constraints
// Cell vectors: W=9,P=3,R=4 => <=10^5
// Cell angles: W=7,P=2,R=3 => <=10^3 (which always works since angles are
// given as >=0 and <=360)
// Cartesian positions: W=12,8,P=3,R=7,3 => <=10^8,10^4 (>0)
// and <=10^7,10^3 (<0 since you need a
// space for the - sign)
// This is all probably fine unless an atom makes it to more negative than
// -999.999.
// Dane Morgan - Stefano Curtarolo
void PrintPDB(const xstructure& str, ostream& oss) {
  oss.setf(std::ios::fixed,std::ios::floatfield);
  oss.precision(10);
  xstructure sstr=str;

  // Set scale to 1 so you don't need to rescale coordinates.
  sstr=ReScale(sstr,1.0);

  // Print out data in pdb format
  oss.setf(std::ios::fixed,std::ios::floatfield);
  oss.precision(3);
  oss << "HEADER    " << sstr.title << endl;
  oss << "REMARK   1" << endl;
  oss << "REMARK   PDB Format: Created by convasp from POSCAR file format." << endl;
  oss << "REMARK   1" << endl;
  xvector<double> data(6);
  data=Getabc_angles(sstr.lattice,DEGREES);//  no need to rescale
  oss << "CRYST1" << setw(9) << data(1)
      << setw(9) << data(2)
      << setw(9) << data(3);
  oss.precision(2);
  oss << setw(7) << data(4)
      << setw(7) << data(5)
      << setw(7) << data(6);
  // Output space group - Use P1 for now
  oss << setw(11) << "P1" << endl;
  // oss << "              P1" << endl;
  // Coordinates in cartesian
  oss.precision(3);
  for(uint i=0;i<sstr.atoms.size();i++) {
    oss << "ATOM  " << setw(5) << i+1 << setw(4) << sstr.atoms.at(i).cleanname.c_str()
        << " " << "MOL" << " H" << setw(4) << 0 << " " << "   "
        << setw(8) << sstr.atoms.at(i).cpos(1)
        << setw(8) << sstr.atoms.at(i).cpos(2)
        << setw(8) << sstr.atoms.at(i).cpos(3)
        << "  1.00" << "  0.00"
        << setw(12) << sstr.atoms.at(i).cleanname.c_str() << "0" << endl;
  }
  oss << "TER" << endl;
}

// ***************************************************************************
// platon2print
// ***************************************************************************
// This funtion prints out structural data in a format for PLATON
// Stefano Curtarolo
void platon2print(xstructure str,bool P_EQUAL,bool P_EXACT,double P_ang,double P_d1,double P_d2,double P_d3, ostream& oss) {
  oss << str.platon2print(P_EQUAL,P_EXACT,P_ang,P_d1,P_d2,P_d3);
}

// **************************************************************************
//  Function PrintRBAnal
// **************************************************************************
// This function prints out analysis of the rubber band run.
// Dane Morgan - Stefano Curtarolo
void PrintRBAnal(const int& nim, const string& path_flag, ostream& oss) {
  int DEFAULT_NPT_FOR_SPLINE=50;
  oss.setf(std::ios::left, std::ios::adjustfield);
  oss.setf(std::ios::fixed, std::ios::floatfield);
  int w=10;
  int p=4;

  // Get data
  vector<string> dirvec=pflow::GetRBDir(nim);
  vector<double> enervec=pflow::GetRBEner(nim);
  vector<xstructure> strvec=pflow::GetRBStruct(nim);
  vector<double> dcum(nim+2);
  vector<double> d_from_init(nim+2);
  vector<double> d_from_final(nim+2);

  // Print data
  oss << "********* Analysis of rubber band run *********" << endl;

  // Exact flag
  string path_flag_temp="E";
  dcum=pflow::GetRBDistCum(strvec,path_flag_temp);
  d_from_init=pflow::GetRBDistFromStrI(strvec,strvec[0],path_flag_temp);
  d_from_final=pflow::GetRBDistFromStrI(strvec,strvec[nim+1],path_flag_temp);
  oss << "** Here we calculate with distance method (eE=exact, nN=nearest): "
      << path_flag_temp << " **" << endl;
  oss << "Dir    Dist(Cum)   Dist(00)    Dist(END)   Ener" <<endl;  
  for(int im=0;im<nim+2;im++) {
    oss << setprecision(p) << setw(5) << dirvec[im]
        << "  " << setw(w) << dcum[im]
        << "  " << setw(w) << d_from_init[im]
        << "  " << setw(w) << d_from_final[im]
        << "  " << setw(w) << enervec[im] << endl;
  }
  // Nearest flag
  path_flag_temp="N";
  dcum=pflow::GetRBDistCum(strvec,path_flag_temp);
  d_from_init=pflow::GetRBDistFromStrI(strvec,strvec[0],path_flag_temp);
  d_from_final=pflow::GetRBDistFromStrI(strvec,strvec[nim+1],path_flag_temp);
  oss << "** Here we calculate with distance method (eE=exact, nN=nearest): "
      << path_flag_temp << " **" << endl;
  oss << "Dir    Dist(Cum)   Dist(00)    Dist(END)   Ener" <<endl;  
  for(int im=0;im<nim+2;im++) {
    oss << setprecision(p) << setw(5) << dirvec[im]
        << "  " << setw(w) << dcum[im]
        << "  " << setw(w) << d_from_init[im]
        << "  " << setw(w) << d_from_final[im]
        << "  " << setw(w) << enervec[im] << endl;
  }

  oss << endl;
  oss << "*******  SPLINE FIT *******" << endl;
  oss << "Distance used in spline (eE=exact, nN=nearest): "
      << path_flag << endl;  
  oss << endl;
  pflow::PrintSpline(dcum,enervec,DEFAULT_NPT_FOR_SPLINE,oss);
}  

// ***************************************************************************
// Function PrintRBPoscarDisp
// ***************************************************************************
// This function prints out the displacement info for 2 POSCAR.
// Dane Morgan - Stefano Curtarolo
void PrintRBPoscarDisp(const xstructure& diffstr, double& totdist,
                       pflow::matrix<double>& cm, const string& path_flag,
                       ostream& oss) {
  oss << "*******  BEGIN:  Structure 2 - Structure 1  *******" << endl;
  oss << diffstr;
  oss << "*******  END:  Structure 2 - Structure 1  *******" << endl;
  oss << endl;
  oss << "Center of mass of POSCARs: "<< endl;
  oss << "  POSCAR1  ";
  pflow::Vout(cm[0],oss);
  oss << "  POSCAR2  ";
  pflow::Vout(cm[1],oss);
  oss << "Total Distance: " << totdist << endl;
  oss << "Method (eE=exact, nN=nearest): " << path_flag << endl;
}

// **************************************************************************
// PrintRDF
// **************************************************************************
// This function prints out the RDF information.
// Dane Morgan - Stefano Curtarolo
void PrintRDF(const xstructure& str, const double& rmax,
              const int& nbins,
              const int& smooth_width,
              const pflow::matrix<double>& rdf_all,
              pflow::matrix<double>& rdfsh_all,
              pflow::matrix<double>& rdfsh_loc,
              ostream& oss) {
  int natoms=str.atoms.size();
  //  vector<int> types=str.GetTypes();
  //  vector<int> num_each_type=str.GetNumEachType();
  int ntypes=str.num_each_type.size();
  oss.setf(std::ios::left, std::ios::adjustfield);
  oss.setf(std::ios::fixed, std::ios::floatfield);
  int w1=5;
  int _rdi=rdfsh_loc.size();
  int _rdj=nbins;
  pflow::matrix<double> rdfsh_bin(_rdi,_rdj);pflow::VVset(rdfsh_bin,0.0);

  // Put rdfsh_all into bins so you can write out in same
  // way as rdf_all
  for(uint i=0;i<rdfsh_loc.size();i++) {
    for(uint j=0;j<rdfsh_loc[i].size();j++) {
      int ib=int((rdfsh_loc[i][j]/rmax)*nbins); // first bin is [0,rmax/nbins).
      rdfsh_bin[i][ib]+=rdfsh_all[i][j];
    }
  }

  // RDF
  oss <<"# Radial Distribution Function (RDF)"<< endl;
  oss << "# rmax= " << rmax << "     nbins = " << nbins << "     sigma= " << smooth_width << endl;
  oss << "# Format: bin  rad  type1:rdf shells  type2:rdf shells ... " << endl;
  double bw=rmax/(double)nbins;
  for(int ia=0;ia<natoms;ia++) {
    oss << "# atom / type: " << ia+1 << " / " << str.atoms.at(ia).type+1 << endl;  // CONVASP_MODE
    oss << "# bin    rad ";
    oss << "       type" << 1;
    for(int it=1;it<ntypes;it++) {
      oss << "          type" << it+1;
    }
    oss <<   "          typeAll";
    oss << endl;
    for(int ib=0;ib<nbins;ib++) {
      oss << "  " << setw(w1) << ib+1;
      oss << "  " << ib*bw;
      oss << "  ";
      for(int it=0;it<ntypes;it++) {
        oss <<" " << setw(w1) << rdf_all[(ntypes+1)*ia+it][ib];
        oss <<" " << setw(w1) << Nint(rdfsh_bin[(ntypes+1)*ia+it][ib]);
      }
      oss <<" " << setw(w1) << rdf_all[(ntypes+1)*ia+ntypes][ib];
      oss <<" " << setw(w1) << Nint(rdfsh_bin[(ntypes+1)*ia+ntypes][ib]);
      oss << endl;
    }

    // RDF SHELLS
    // Header line
    oss << "# atom / type: " << ia+1 << " / " << str.atoms.at(ia).type+1 << endl;  // CONVASP_MODE
    oss << "# NN   ";
    for(int it=0;it<ntypes;it++) {
      oss << "  type" << it+1;
      oss << ":rad/size";
    }
    oss <<   "  typeAll";
    oss << ":rad/size";
    oss << endl;
    // Get nshmx for this atom
    int nshmx=0;
    for(int i=0;i<ntypes+1;i++) {
      int temp=rdfsh_all[(ntypes+1)*ia+i].size();
      if(temp>nshmx) nshmx=temp;
    }
    // Output shell data
    for(int is=0;is<nshmx;is++) {
      oss << "  " << setw(w1) << is+1;
      for(int it=0;it<ntypes+1;it++) {
        int id=(ntypes+1)*ia+it;
        if(is<(int) rdfsh_all[id].size()) {
          oss << "  " << setw(w1+2) << setprecision(4) << rdfsh_loc[id][is];
          oss << "  " << setw(w1) << Nint(rdfsh_all[id][is]);
        }
        else{
          oss << "  " << setw(w1) << "NA";
          oss << "  " << setw(w1) << "NA";
        }
      } // for it
      oss << endl;
    } // for is
  } // for ia

}

// **************************************************************************
// PrintRDFCmp
// **************************************************************************
// This function prints out the RDF shells comparison information.
// Dane Morgan - Stefano Curtarolo
void PrintRDFCmp(const xstructure& str_A, const xstructure& str_B,
                 const double& rmax, const int nbins,
                 const double& smooth_width,const int nsh,
                 const pflow::matrix<double>& rdfsh_all_A,
                 const pflow::matrix<double>& rdfsh_all_B,
                 const vector<int>& best_match,
                 const pflow::matrix<double>& rms_mat, ostream& oss) {
  int nat=str_A.atoms.size();   if(str_B.atoms.size()) {;} // phony to keep str_B busy
  std::deque<int> num_each_type=str_A.num_each_type;
  int nt=num_each_type.size();
  int w1=4;
  oss.setf(std::ios::left, std::ios::adjustfield);
  oss.setf(std::ios::fixed, std::ios::floatfield);

  // Info on RDF calculation parameters.
  oss <<"# Comparison of structures based on radial distribution functions (RDFs)"<< endl;
  oss << "# rmax= " << rmax << "     nbins = " << nbins << endl;
  oss << "# sigma= " << smooth_width << "     nsh = " << nsh << endl;
  oss << endl;
  oss << "# The matrix of rms error between atoms in str A and str B" << endl;

  // The RMS matrix.

  // header line
  oss << setw(w1+1) << " ";
  for(int i=0;i<nat;i++) {
    oss << " " << setw(w1) << i+1;
  }
  oss << endl;

  // The RMS matrix as a row/column matrix.
  for(int i=0;i<nat;i++) {
    oss << setw(w1) << i+1 << " ";
    for(int j=0;j<nat;j++) {
      oss << " " << setprecision(2) << rms_mat[i][j];
    }
    oss << endl;
  }

  // The fit data for each atom's best fit.
  double rms_tot=0;
  oss << endl;

  // Basic atoms and rms.
  for(int iA=0;iA<nat;iA++) {
    int iB=best_match[iA];
    rms_tot+=rms_mat[iA][iB];
    oss << "A atom / Best matching B atom / RMS: "
        << iA+1 << " " << iB+1 << " " << rms_mat[iA][iB] << endl;

    // Detailed data on the A and B shells
    // Header line
    oss << "# CN ";
    oss << "   type" << 1;
    for(int it=1;it<nt;it++) {
      oss << "             type" << it+1;
    }
    oss <<   "             typeAll";
    oss << endl;

    // Get nshmx for this atom
    int nshmx=0;
    int temp;
    for(int i=0;i<nt+1;i++) {
      temp=rdfsh_all_A[(nt+1)*iA+i].size();
      if(temp>nshmx) nshmx=temp;
      temp=rdfsh_all_B[(nt+1)*iB+i].size();
      if(temp>nshmx) nshmx=temp;
    }

    // Output shell data
    for(int is=0;is<nshmx;is++) {
      oss << "  " << setw(w1) << is+1;
      for(int it=0;it<nt+1;it++) {
        int idA=(nt+1)*iA+it;
        int idB=(nt+1)*iB+it;
        if(is<(int) rdfsh_all_A[idA].size()) {
          oss << " " << setw(w1) << Nint(rdfsh_all_A[idA][is]);
        }
        else{
          oss << " " << setw(w1) << "NA";
        }
        if(is<(int) rdfsh_all_B[idB].size()) {
          oss << " " << setw(w1) << Nint(rdfsh_all_B[idB][is]);
        }
        else{
          oss << " " << setw(w1) << "NA";
        }
        if(is<(int) rdfsh_all_A[idA].size() && is<(int) rdfsh_all_B[idB].size()) {
          oss << " " << setw(w1) << Nint(rdfsh_all_A[idA][is]-rdfsh_all_B[idB][is]);
        }
        else{
          oss << " " << setw(w1) << "NA";
        }
        oss << "   ";
      } // for it
      oss << endl;
    } // for is

  } // for iA

  // Total RMS
  oss << "RMS TOTAL: " << rms_tot/nat;
}

// ***************************************************************************
// PrintRSM
// ***************************************************************************
// This funtion prints out structural data in a Rasmol (RSM) format with unit cell wireframe
// Wahyu Setyawan
void PrintRSM(const xstructure& str, ostream& oss) {
  string atName[]={"H    ","He   ","Li   ","Be   ","B    ","C    ","N    ","O    ","F    ","Ne   ",
                   "Na   ","Mg   ","Al   ","Si   ","P    ","S    ","Cl   ","Ar   ","K    ","Ca   ",
                   "Sc   ","Ti   ","V    ","Cr   ","Mn   ","Fe   ","Co   ","Ni   ","Cu   ","Zn   ",
                   "Ga   ","Ge   ","As   ","Se   ","Br   ","Kr   ","Rb   ","Sr   ","Y    ","Zr   ",
                   "Nb   ","Mo   ","Tc   ","Ru   ","Rh   ","Pd   ","Ag   ","Cd   ","In   ","Sn   ",
                   "Sb   ","Te   ","I    ","Xe   ","Cs   ","Ba   ","La   ","Ce   ","Pr   ","Nd   ",
                   "Pm   ","Sm   ","Eu   ","Gd   ","Tb   ","Dy   ","Ho   ","Er   ","Tm   ","Yb   ",
                   "Lu   ","Hf   ","Ta   ","W    ","Re   ","Os   ","Ir   ","Pt   ","Au   ","Hg   ",
                   "Tl   ","Pb   ","Bi   ","Po   ","At   ","Rn   ","Fr   ","Ra   ","Ac   ","Th   ",
                   "Pa   ","U    ","Np   ","Pu   ","Am   ","Cm   ","Bk   ","Cf   ","Es   ","Fm   "};
  oss<<"#Title"<<endl
     <<"#23456#8901##4567#9012#456#890#2345678#0123456#8901234#67890#23456789"<<endl;
  int num_atoms;
  num_atoms=0;
  for(uint i=0;i<str.num_each_type.size();i++)
    num_atoms+=str.num_each_type.at(i);

  oss.setf(std::ios::fixed);
  for(int i=0;i<num_atoms;i++) {
    oss<<"ATOM  "
       <<setw(5)<<i+9
       <<" "
       <<atName[-1+str.atoms.at(i).type]
       <<"UNK  "
       <<setw(4)<<str.atoms.at(i).type
       <<"    "
       <<setprecision(3)<<setw(8)<<str.scale*str.atoms.at(i).cpos(1)
       <<setprecision(3)<<setw(8)<<str.scale*str.atoms.at(i).cpos(2)
       <<setprecision(3)<<setw(8)<<str.scale*str.atoms.at(i).cpos(3)<<endl;
    //        fprintf(fid,'%6s%5d%1s%5s%5s %3d%4s%8.4f%8.4f%8.4f\n', ...^M
    //            'ATOM  ',ptr+8,' ',atName((label(i)-1)*5+1:label(i)*5),'UNK  ',label(i),'    ',x,y,z);^M
  }

  double latt1[3]={str.scale*str.lattice(1,1),str.scale*str.lattice(1,2),str.scale*str.lattice(1,3)};
  double latt2[3]={str.scale*str.lattice(2,1),str.scale*str.lattice(2,2),str.scale*str.lattice(2,3)};
  double latt3[3]={str.scale*str.lattice(3,1),str.scale*str.lattice(3,2),str.scale*str.lattice(3,3)};

  double cellframe[8][3]={{0,0,0},
                          {latt1[0],latt1[1],latt1[2]}, //latt(1,:)
                          {latt2[0],latt2[1],latt2[2]}, //latt(2,:)
                          {latt3[0],latt3[1],latt3[2]}, //latt(3,:)
                          {latt1[0]+latt2[0],latt1[1]+latt2[1],latt1[2]+latt2[2]}, //latt(1,:)+latt(2,:);
                          {latt1[0]+latt3[0],latt1[1]+latt3[1],latt1[2]+latt3[2]}, //latt(1,:)+latt(3,:);
                          {latt3[0]+latt2[0],latt3[1]+latt2[1],latt3[2]+latt2[2]}, //latt(3,:)+latt(2,:);
                          {latt1[0]+latt2[0]+latt3[0],latt1[1]+latt2[1]+latt3[1],latt1[2]+latt2[2]+latt3[2]}};

  for(int i=0;i<8;i++) {
    oss<<"ATOM  "
       <<setw(5)<<i+1
       <<" "
       <<"X    "
       <<"UNK  "
       <<setw(4)<<0
       <<"    "
       <<setprecision(3)<<setw(8)<<cellframe[i][0]
       <<setprecision(3)<<setw(8)<<cellframe[i][1]
       <<setprecision(3)<<setw(8)<<cellframe[i][2]<<endl;
    //        fprintf(fid,'%6s%5d%1s%5s%5s %3d%4s%8.4f%8.4f%8.4f\n', ...^M
    //            'ATOM  ',i,' ','X    ','UNK  ',0,'    ',x,y,z);^M
  }

  int lineframe[12][2]={{1,2},  {1,3},  {1,4},  {8,7},  {8,6},  {8,5},  {2,5},  {2,6},  {3,5},  {3,7},  {4,6},  {4,7}};
  for(int i=0;i<12;i++) {
    oss<<"CONECT "
       <<setw(4)<<lineframe[i][0]
       <<" "
       <<setw(4)<<lineframe[i][1]<<endl;
    //    fprintf(fid,'%6s %4d %4d\n','CONECT',lineframe(i,1),lineframe(i,2));^M
  }

  oss<<"END"<<endl
     <<"#!rasmol -rsm file"<<endl
     <<"load pdb inline"<<endl
     <<"spacefill 150"<<endl
     <<"select 0"<<endl
     <<"spacefill 10"<<endl
     <<"exit"<<endl;

  /*
    %#23456#8901##4567#9012#456#890#2345678#0123456#8901234#67890#23456789^M
    %COLO      1 Al   UNK  0001       0.500   0.500   0.000   0.2^M
    %ATOM      1 Al   UNK  0001       0.000   0.000   0.000^M
    %ATOM      2 O    UNK  0002       2.000   0.000   0.000^M
    %CONECTxxxxxyyyyy^M
  */
}

//DX 8/30/17 - SGDATA
// ***************************************************************************
// pflow::PrintSGData
// ***************************************************************************
// This funtion prints out space group data.
// David Hicks (DX)
namespace pflow {
    bool PrintSGData(xstructure& str_sg, ostream& oss, bool standalone, const string& format,bool already_calculated) {
  double tolerance=SYM::defaultTolerance(str_sg);
      if(already_calculated){tolerance=str_sg.sym_eps;} //CO 171025
  bool no_scan=false;
      return PrintSGData(str_sg,tolerance,oss,no_scan,standalone,format,already_calculated);
  }
}
namespace pflow {
    bool PrintSGData(xstructure& str_sg, double& tolerance, ostream& oss_final, bool no_scan, bool standalone, const string& format,bool already_calculated) { // DX 2/26/18 - added & to tolerance
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "pflow::PrintSGData: BEGIN" << endl;
    stringstream oss;
    // smode=="DATA" or "EDATA"
    // Print out structural data
    oss.setf(std::ios::fixed,std::ios::floatfield);
    //DX TEST str_sg.SpaceGroup_ITC(tolerance);
      if(!already_calculated){str_sg.SpaceGroup_ITC(tolerance, -1, no_scan);} //CO 171025
    // FORMAT = TXT
    if(format=="txt"){
      oss << "SPACE GROUP OF THE CRYSTAL" << endl;
      oss << " Space group number                           = " << str_sg.space_group_ITC << endl;
      oss << " Space group label (Hermann Mauguin)          = " << GetSpaceGroupName(str_sg.space_group_ITC) << endl;
      oss << " Space group label (Hall)                     = " << GetSpaceGroupHall(str_sg.space_group_ITC) << endl;
      oss << " Space group label (Schoenflies)              = " << GetSpaceGroupSchoenflies(str_sg.space_group_ITC) << endl;
      oss << " Laue class                                   = " << GetLaueLabel(str_sg.point_group_ITC) << endl;
      oss << " Crystal class                                = " << str_sg.point_group_ITC << endl;
      oss << "ITC REPRESENTATION OF CRYSTAL" << endl;
      oss << " Setting                                      = " << str_sg.setting_ITC << endl;
      oss << " Origin                                       = " << roundoff(str_sg.origin_ITC,1e-8) << endl;
      oss << " General Wyckoff position" << endl;
      for(uint i=0;i<str_sg.general_position_ITC.size();i++){
	oss << "  " << i+1 << " " << str_sg.general_position_ITC[i] << endl;
      }
      oss << " Representative Wyckoff positions" << endl;
      for(uint i=0;i<str_sg.wyccar_ITC.size();i++){
        if(i>4 && i!=str_sg.wyccar_ITC.size()-1){ //Skip title, scale, lattice parameters, number of atoms, and coordinate type, and last newline
	  oss << str_sg.wyccar_ITC[i] << endl;
        }
      }
      oss << "WYCCAR" << endl;
      // convert vector<string> of WYCCAR to xstructure
      stringstream wss; 
      for(uint i=0;i<str_sg.wyccar_ITC.size();i++){
        wss << str_sg.wyccar_ITC[i] << endl;
      }
      xstructure xstr_wyccar(wss,IOVASP_WYCKCAR);
      oss << xstr_wyccar << endl;
    }
    // FORMAT = JSON
    else if(format=="json"){
      string eendl="";
      bool roff=true; //round off
      bool PRINT_NULL=FALSE;
      stringstream sscontent_json;
      vector<string> vcontent_json;

      // space group number
      if(str_sg.space_group_ITC){
        sscontent_json << "\"space_group_number\":\"" << str_sg.space_group_ITC << "\"" << eendl;
      } else{
        if(PRINT_NULL){ sscontent_json << "\"space_group_number\":null" << eendl;}
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

      // space group label (HM)
      string sg_HM = GetSpaceGroupName(str_sg.space_group_ITC);
      if(sg_HM.size()){
        sscontent_json << "\"space_group_Hermann_Mauguin\":\"" << sg_HM << "\"" << eendl;
      } else{
        if(PRINT_NULL){ sscontent_json << "\"space_group_Hermann_Mauguin\":null" << eendl;}
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

      // space group label (Hall)
      string sg_Hall = GetSpaceGroupHall(str_sg.space_group_ITC);
      if(sg_Hall.size()){
        sscontent_json << "\"space_group_Hall\":\"" << sg_Hall << "\"" << eendl;
      } else{
        if(PRINT_NULL){ sscontent_json << "\"space_group_Hall\":null" << eendl;}
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      
      // space group label (Schoenflies)
      string sg_Schoenflies = GetSpaceGroupSchoenflies(str_sg.space_group_ITC);
      if(sg_Schoenflies.size()){
        sscontent_json << "\"space_group_Schoenflies\":\"" << sg_Schoenflies << "\"" << eendl;
      } else{
        if(PRINT_NULL){ sscontent_json << "\"space_group_Schoenflies\":null" << eendl;}
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      
      // Laue
      string laue = GetLaueLabel(str_sg.point_group_ITC);
      if(laue.size()){
        sscontent_json << "\"Laue\":\"" << laue << "\"" << eendl;
      } else{
        if(PRINT_NULL){ sscontent_json << "\"Laue\":null" << eendl;}
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      
      // crystal class
      if(str_sg.point_group_ITC.size()){
        sscontent_json << "\"crystal_class\":\"" << str_sg.point_group_ITC << "\"" << eendl;
      } else{
        if(PRINT_NULL){ sscontent_json << "\"crystal_class\":null" << eendl;}
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      
      // ITC setting
      if(str_sg.setting_ITC){
        sscontent_json << "\"setting_ITC\":\"" << str_sg.setting_ITC << "\"" << eendl;
      } else{
        if(PRINT_NULL){ sscontent_json << "\"setting_ITC\":null" << eendl;}
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

      // ITC origin
      if(str_sg.origin_ITC.rows){
        sscontent_json << "\"origin_ITC\":[" << aurostd::joinWDelimiter(xvecDouble2vecString(str_sg.origin_ITC,5,roff),",") << "]" << eendl;
      } else{
        if(PRINT_NULL){ sscontent_json << "\"origin_ITC\":null" << eendl;}
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      
      // general Wyckoff position
      if(str_sg.general_position_ITC.size()){
        sscontent_json << "\"general_position_ITC\":[";
        for(uint i=0;i<str_sg.general_position_ITC.size();i++){
          vector<string> tokens;
          aurostd::string2tokens(str_sg.general_position_ITC[i],tokens,",");
          sscontent_json << "[" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(tokens,"\""),",") << "]" << eendl;
          if(i != str_sg.general_position_ITC.size()-1){
            sscontent_json << ",";
          }
        }
        sscontent_json << "]";
      } else{
        if(PRINT_NULL){ sscontent_json << "\"setting_ITC\":null" << eendl;}
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
     
      // representative Wyckoff positions
      if(str_sg.wyccar_ITC.size()){
        sscontent_json << "\"Wyckoff_positions\":["<< eendl;
        vector<string> wyckoff_set;
        for(uint i=0;i<str_sg.wyccar_ITC.size();i++){
          if(i>4 && i!=str_sg.wyccar_ITC.size()-1){ //Skip title, scale, lattice parameters, number of atoms, and coordinate type, and last newline
            vector<string> tokens;
            aurostd::string2tokens(str_sg.wyccar_ITC[i],tokens," ");
            string name = "\"name\":";
            xvector<double> position;
            string multiplicity = "\"multiplicity\":";
            string wyckoff_letter = "\"Wyckoff_letter\":";
            string site_symmetry = "\"site_symmetry\":";
            stringstream sswyckoff;
            vector<string> vwyckoff_json;
            for(uint t=0;t<tokens.size();t++){
              if(t==0){
                position(1) = aurostd::string2utype<double>(tokens[t]);
              }
              if(t==1){
                position(2) = aurostd::string2utype<double>(tokens[t]);
              }
              if(t==2){
                position(3) = aurostd::string2utype<double>(tokens[t]);
              }
              if(t==3){
                name += "\""+tokens[t]+"\"";
              }
              if(t==4){
                multiplicity += tokens[t];
              }
              if(t==5){
                wyckoff_letter += "\""+tokens[t]+"\"";
              }
              if(t==6){
                site_symmetry += "\""+tokens[t]+"\"";
              }
            }
            sswyckoff << "\"position\":[" << aurostd::joinWDelimiter(xvecDouble2vecString(position,5,roff),",") << "]" << eendl;
            vwyckoff_json.push_back(sswyckoff.str()); sswyckoff.str("");
            sswyckoff << name << eendl;
            vwyckoff_json.push_back(sswyckoff.str()); sswyckoff.str("");
            sswyckoff << multiplicity << eendl;
            vwyckoff_json.push_back(sswyckoff.str()); sswyckoff.str("");
            sswyckoff << wyckoff_letter << eendl;
            vwyckoff_json.push_back(sswyckoff.str()); sswyckoff.str("");
            sswyckoff << site_symmetry << eendl;
            vwyckoff_json.push_back(sswyckoff.str()); sswyckoff.str("");
            sswyckoff << "{" << aurostd::joinWDelimiter(vwyckoff_json,",") << "}" << eendl;
            wyckoff_set.push_back(sswyckoff.str()); sswyckoff.str("");
          }
        }
        sscontent_json << aurostd::joinWDelimiter(wyckoff_set,",") << eendl;  	    
        sscontent_json << "]" << eendl;  	    
        vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      }
 
      // WYCCAR
      if(str_sg.wyccar_ITC.size()){
        // convert vector<string> of WYCCAR to xstructure
        stringstream wss; 
        for(uint i=0;i<str_sg.wyccar_ITC.size();i++){
          wss << str_sg.wyccar_ITC[i] << endl;
        }
        xstructure xstr_wyccar(wss,IOVASP_WYCKCAR);
        
        sscontent_json << "\"wyccar\":" << xstructure2json(xstr_wyccar) << eendl;
      } else{
        if(PRINT_NULL){ sscontent_json << "\"wyccar\":null" << eendl;}
      }
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

      if(standalone){
        oss << "{" << aurostd::joinWDelimiter(vcontent_json,",")  << "}" << endl;
      }
      else{
        oss << aurostd::joinWDelimiter(vcontent_json,",");
      }
    }
    oss_final << oss.str();
    return true;
  }
}      
//DX 8/30/17 - SGDATA

// **************************************************************************
// PrintShell PrintShell
// **************************************************************************
// This function prints out all points in unit cell with
//  ns atoms between rmin and rmax.
// --------------------------------------------------
// GetAngleMat
// Given a point q and list of points, {p_i} returns the xmatrix
// M_ij = Angle from pi-q-pj.
// Dane Morgan - Stefano Curtarolo
pflow::matrix<double> GetAngleMat(const xvector<double>& q, const vector<xvector<double> >& p) {
  double TOL=1e-6;
  int np=p.size();
  pflow::matrix<double> M(np,np);pflow::VVset(M,0.0);
  for(int i=0;i<np;i++) {
    for(int j=0;j<np;j++) {
      xvector<double> a(3);a=p.at(i)-q;
      xvector<double> b(3);b=p.at(j)-q;
      //  double na=aurostd::modulus(a); // DANE not used
      //  double nb=aurostd::modulus(b); // DANE not used
      double angle;
      double xcos=getcos(a,b);
      if(xcos>(1-TOL)) {
        angle=0;
      } else{
        if(xcos<(-1+TOL)) {
          angle=180;
        } else{
          angle=rad2deg*acos(xcos);
        }
      }
      M[i][j]=angle;
    }
  }
  return M;
}

// ***************************************************************************
// PrintShell GetGoodShellPoints
// ***************************************************************************
// This function collects the information about atoms
// within a shell around a set of points.  The points are
// input in fractional coordinates.  The points with the
// correct shells (ns atoms between rmin and rmax) are
// stored, with their neighbors, as atoms in neigh_mat.
// Dane Morgan - Stefano Curtarolo
void GetGoodShellPoints(vector<xvector<double> >& points, const xstructure& str,
                        const int& ns, const double& rmin,const double& rmax,
                        const string& sname,deque<deque<_atom> >& neigh_mat) {
  xstructure sstr=str;
  sstr=ReScale(sstr,1.0);
  xmatrix<double> lat(3,3);
  lat=sstr.lattice;

  // define fake atoms for each point
  deque<_atom> atvec;
  for(uint ipt=0;ipt<points.size();ipt++) {
    _atom a;
    // Set positions, other defaults are fine for a.
    a.fpos=points.at(ipt);
    a.cpos=F2C(lat,a.fpos);
    atvec.push_back(a);
  }

  // Send atoms to GetNeighData to get their neigh info.
  //[OBSOLETE]  pflow::GetNeighData(atvec,str,rmin,rmax,neigh_mat);
  sstr.GetNeighData(atvec,rmin,rmax,neigh_mat);

  // Remove every point that does not have ns atoms in shell.
  deque<deque<_atom> > nnm;
  for(uint ip=0;ip<neigh_mat.size();ip++) {
    if((int) neigh_mat.at(ip).size()==(ns+1)) { // size is number of neighbors + the point.
      nnm.push_back(neigh_mat.at(ip));
    }
  }
  // Remove every point that does not have neighbors with the required name.
  // If sname==NONE then keep all the points.
  deque<deque<_atom> > nnnm;
  string none="NONE";
  if(!(sname==none)) {
    for(uint ip=0;ip<nnm.size();ip++) {
      int names_match=1;
      for(uint in=1;in<nnm.at(ip).size();in++) {
        if(!(sname==nnm.at(ip).at(in).name)) {
          names_match=0;
        }
      } // for in
      if(names_match) nnnm.push_back(nnm.at(ip));
    } // for ip
    neigh_mat=nnnm;
  } // if sname
  else{
    neigh_mat=nnm;
  } // if sname
}

// ***************************************************************************
// PrintShell new_point
// ***************************************************************************
// return true for a point with a different average
// neighbor position.
// Dane Morgan - Stefano Curtarolo
bool new_point(const vector<xvector<double> >& avnfpos, int na) {
  for(int ia=0;ia<na;ia++)
    if(identical(avnfpos.at(ia),avnfpos.at(na),1e-15)) return FALSE;
  return TRUE;
}

// --------------------------------------------------
// PrintShell
// Dane Morgan - Stefano Curtarolo
void PrintShell(const xstructure& str, const int& ns,const double& rmin, const double& rmax,const string& sname, const int lin_dens, ostream& oss) {  
  oss.setf(std::ios::left, std::ios::adjustfield);
  oss.setf(std::ios::fixed, std::ios::floatfield);
  xstructure sstr=str;
  // double scale=str.scale; // DANE not used
  sstr=ReScale(sstr,1.0);
  xmatrix<double> lattice=sstr.lattice;
  //  oss.setf(std::ios::fixed,std::ios::floatfield);
  //  oss.precision(10);
  deque<deque<_atom> > neigh_mat;
  vector<xvector<double> > points;
  for(int ix=0;ix<lin_dens;ix++) {
    for(int iy=0;iy<lin_dens;iy++) {
      for(int iz=0;iz<lin_dens;iz++) {
        xvector<double> pt(3);
        pt(1)=ix/(double)lin_dens;
        pt(2)=iy/(double)lin_dens;
        pt(3)=iz/(double)lin_dens;
        points.push_back(pt);
      }
    }
  }
  GetGoodShellPoints(points,str,ns,rmin,rmax,sname,neigh_mat);
  // double tol=1e-15; // DANE not used

  // Get average fpos of neighbors of each point.
  vector<xvector<double> > avnfpos;
  for(uint ia=0;ia<neigh_mat.size();ia++) {
    xvector<double> av(3);
    for(uint in=1;in<neigh_mat.at(ia).size();in++) {
      _atom a=neigh_mat.at(ia).at(in);
      //   xvector<double> pos(3);
      //  pos=a.fpos;
      av=av+a.fpos;
    } // for in
    for(int ic=1;ic<=3;ic++) {
      av(ic)=av(ic)/(neigh_mat.at(ia).size()-1);
    } // for ic
    avnfpos.push_back(av);
  }
  // Output header
  oss << "All atoms identified by:  ";
  oss << "[basis number]  [name]  [n1 n2 n3 (unit cell)]" << endl;
  oss << "Points must have " << ns << " atoms around them within a shell bounded by " << rmin << " and " << rmax << endl;
  oss << "Atoms must be of type: " << sname << ", ";
  oss << "Linear density of searched points is: " << lin_dens << endl;
  oss << "All coordinates are direct (cartesian) if the input POSCAR is direct (cartesian)." << endl;

  // Print out all points with unique environments.
  oss << endl;
  oss << "********** Unique environment points meeting shell constraints **********" << endl;
  int newcnt=0;
  vector<xvector<double> > all_new_pts;
  for(uint ia=0;ia<neigh_mat.size();ia++) {
    vector<xvector<double> > p; // For neighbors of averge point q.
    if(new_point(avnfpos,ia)) {
      newcnt++;
      _atom a;
      xvector<double> pt(3);pt=avnfpos.at(ia);
      all_new_pts.push_back(pt);
      a.fpos=pt;
      a.cpos=F2C(lattice,pt);
      if(sstr.coord_flag==FALSE) { // direct
        oss << setprecision(10) << setw(7) << ia+1 << "(" << newcnt << ")"
            << "  UAvg "
            << " " << setw(10) << pt(1)
            << " " << setw(10) << pt(2)
            << " " << setw(10) << pt(3) << endl;
      }
      if(sstr.coord_flag==TRUE) { // cart
        xvector<double> cv(3);cv=F2C(lattice,pt);
        oss << setprecision(10) << setw(7) << ia+1 << "(" << newcnt << ")"
            << "  UAvg "
            << " " << setw(10) << cv(1)
            << " " << setw(10) << cv(2)
            << " " << setw(10) << cv(3) << endl;
      }

      // Print out avg position of neighbors.
      //      oss << "        ";
      //      oss << setw(4) << "Avg ";
      //      oss << setw(4) << " ";
      //      if(sstr.coord_flag==0) { // direct
      //      oss << setprecision(10)<< setw(10) << avnfpos.at(ia)(1) << " " << setw(10) << avnfpos.at(ia)(2) << " " << setw(10) << avnfpos.at(ia)(3) << endl;
      //      }
      //      if(sstr.coord_flag==1) { // cart
      //      vector<double> cv=F2C(lattice,avnfpos.at(ia));
      //      oss << setprecision(10)<< setw(10) << cv(1) << " " << setw(10) << cv(2) << " " << setw(10) << cv(3) << endl;
      //      }
      // Print out neighbors.
      for(uint in=1;in<neigh_mat.at(ia).size();in++) {
        _atom an=neigh_mat.at(ia)[in];
        xvector<int> ijk(3);
        ijk=an.ijk;
        oss << "        ";
        oss << setw(4) << an.number+1 << " ";
        oss << setw(4) << an.name.c_str() << "   ";
        oss << setw(3) << ijk(1) << " " << setw(3) << ijk(2) << " " << setw(3) << ijk(3) << "   ";
        oss << setprecision(4) << AtomDist(a,an);
        oss << endl;
        p.push_back(an.cpos);
      } // in

      // Print out angles with neighbors.
      xvector<double> q(3);q=F2C(lattice,avnfpos.at(ia));
      pflow::matrix<double> ang_mat=GetAngleMat(q,p);
      int np=p.size();
      oss << "        ";
      oss << "Angles"<<endl;
      oss << "           ";
      for(int i=0;i<np;i++) {
        oss << " " << setw(6) << i+1;
      }
      oss << endl;
      for(int i=0;i<np;i++) {
        oss << "        ";
        oss << setw(2) << i+1 ;
        for(int j=0;j<np;j++) {
          oss << " " << setprecision(2) << setw(6) << ang_mat[i][j];
        }
        oss << endl;
      }
    } // if newpoint
  } //ia

  // Print out all points.
  oss << endl;
  oss << "********** All points meeting shell constraints **********" << endl;
  for(uint ia=0;ia<neigh_mat.size();ia++) {
    _atom a=neigh_mat.at(ia).at(0);
    xvector<double> pt(3);pt=a.fpos;
    if(sstr.coord_flag==FALSE) { // direct
      oss << setprecision(10) << setw(7) << ia+1
          << " " << setw(10) << pt(1)
          << " " << setw(10) << pt(2)
          << " " << setw(10) << pt(3) << endl;
    }
    if(sstr.coord_flag==TRUE) { // cart
      xvector<double> cv(3);cv=F2C(lattice,pt);
      oss << setprecision(10) << setw(7) << ia+1
          << " " << setw(10) << cv(1)
          << " " << setw(10) << cv(2)
          << " " << setw(10) << cv(3) << endl;
    }
    for(uint in=1;in<neigh_mat.at(ia).size();in++) {
      _atom an=neigh_mat.at(ia).at(in);
      xvector<int> ijk(3);
      ijk=an.ijk;
      oss << "        ";
      oss << setw(4) << an.number+1 << " ";
      oss << setw(4) << an.name.c_str() << "   ";
      oss << setw(3) << ijk(1) << " " << setw(3) << ijk(2) << " " << setw(3) << ijk(3) << "   ";
      oss << setprecision(4) << AtomDist(a,an);
      oss << endl;
    } // in
    // Print out avg position of neighbors.
    oss << "        ";
    oss << setw(4) << "Avg ";
    oss << setw(4) << " ";
    if(sstr.coord_flag==FALSE) { // direct
      oss << setprecision(10)<< setw(10) << avnfpos.at(ia)(1) << " " << setw(10) << avnfpos.at(ia)(2) << " " << setw(10) << avnfpos.at(ia)(3) << endl;
    }
    if(sstr.coord_flag==TRUE) { // cart
      xvector<double> cv(3);cv=F2C(lattice,avnfpos.at(ia));
      oss << setprecision(10)<< setw(10) << cv(1) << " " << setw(10) << cv(2) << " " << setw(10) << cv(3) << endl;
    }
  } // ia

  // WORKING HERE ???  Finish this and check consistency of
  //  setting names, pos, nat, and particularly types - xstructure should
  //    keep all these consistent when any one is changed.
  //  Add new points to the POSCAR file and output.
  oss << endl;
  oss << "********** POSCAR file with unique environment points meeting shell constraints **********" << endl;
  // Number of atoms

  uint nnat=all_new_pts.size();
  uint nat=str.atoms.size();
  // uint tnat=nat+nnat;  // DANE not used
  // Number of types
  sstr.num_each_type.push_back(nnat);
  sstr.comp_each_type.push_back((double) nnat);
  for(uint iat=0;iat<nnat;iat++) {
    _atom newatom;
    newatom.fpos=all_new_pts.at(iat);
    newatom.cpos=F2C(sstr.lattice,newatom.fpos);
    newatom.corigin=sstr.origin;
    newatom.type=sstr.num_each_type.size()-1;   // CONVASP_MODE
    ostringstream ostmp;ostmp<<iat+1;
    newatom.name="SHELL"+ostmp.str();
    newatom.name_is_given=TRUE;
    newatom.cleanname="XX";
    newatom.sd="";
    newatom.atomic_number=-1;
    newatom.number=nat+iat;                 // reference position for convasp
    newatom.basis=nat+iat;                  // position in the basis	
    clear(newatom.ijk);                     // position in the unit cell
    sstr.atoms.push_back(newatom);
  }
  //  sstr.write_DEBUG_flag=TRUE;
  oss << sstr;
} // end routine

// void PrintShell(const xstructure& str, const int& ns,const double& rmin,
// const double& rmax,const string& sname, const int lin_dens) {  

// ***************************************************************************
// PrintSumDOS
// ***************************************************************************
// Prints out the summed pdos in user friendly format.
// Dane Morgan - Stefano Curtarolo
void PrintSumPDOS(pflow::pdosdata& pdd, ostream& oss) {
  int spin=pdd.spin-1;
  pdd.PrintPDOS(oss,spin);
  if(pdd.print_params==1) {
    pflow::projdata prd;
    std::vector<string> LMnames=prd.LMnames;
    pdd.PrintParams(oss,LMnames);
  }
}

// ***************************************************************************
// PrintXYZ  PrintXYZws PrintXYZInSphere PrintXYZNanoparticle
// ***************************************************************************
// This funtion prints out structural data in an XYZ format
// Stefano Curtarolo
void PrintXYZ(const xstructure& a, const xvector<int>& n, ostream& oss) {
  oss.setf(std::ios::fixed,std::ios::floatfield);
  oss.precision(10);
  xstructure sstr=a;
  sstr=ReScale(sstr,1.0);
  // Print out data in XYZ format
  oss << sstr.atoms.size()*n(1)*n(2)*n(3) << endl;
  oss << sstr.title << endl;
  for(int i=0;i<n(1);i++) {
    for(int j=0;j<n(2);j++) {
      for(int k=0;k<n(3);k++) {
        for(uint iat=0;iat<sstr.atoms.size();iat++) {
          oss << sstr.atoms.at(iat).cleanname << " "
              << (sstr.atoms.at(iat).cpos(1)+i*sstr.lattice(1,1)+j*sstr.lattice(2,1)+k*sstr.lattice(3,1)) << " "
              << (sstr.atoms.at(iat).cpos(2)+i*sstr.lattice(1,2)+j*sstr.lattice(2,2)+k*sstr.lattice(3,2)) << " "
              << (sstr.atoms.at(iat).cpos(3)+i*sstr.lattice(1,3)+j*sstr.lattice(2,3)+k*sstr.lattice(3,3)) << " " << endl;
        } // iat
      } // k
    } // j
  } // i
}

// ***************************************************************************
// PrintXYZws
// ***************************************************************************
// This funtion moves atoms inside WS cell and prints out structural
// data in an XYZ format
// Stefano Curtarolo
void PrintXYZws(const xstructure& a, ostream& oss) {
  // SC 10Jan2003
  oss.setf(std::ios::fixed,std::ios::floatfield);
  oss.precision(10);
  xstructure sstr=a;
  // Set scale to 1 so you don't need to rescale coordinates.
  sstr=ReScale(sstr,1.0);
  // zoology of vectors
  xvector<double> rat(3),xoo(3),yoo(3),zoo(3),xyo(3),xzo(3),yzo(3),xyz(3);
  double axoo,ayoo,azoo,axyo,axzo,ayzo,axyz;
  xoo=sstr.lattice(1);axoo=aurostd::modulus(xoo);
  yoo=sstr.lattice(2);ayoo=aurostd::modulus(yoo);
  zoo=sstr.lattice(3);azoo=aurostd::modulus(zoo);
  xyo=sstr.lattice(1)+sstr.lattice(2);axyo=aurostd::modulus(xyo);
  xzo=sstr.lattice(1)+sstr.lattice(3);axzo=aurostd::modulus(xzo);
  yzo=sstr.lattice(2)+sstr.lattice(3);ayzo=aurostd::modulus(yzo);
  xyz=sstr.lattice(1)+sstr.lattice(2)+sstr.lattice(3);axyz=aurostd::modulus(xyz);
  double projxoo,projyoo,projzoo,projxyo,projxzo,projyzo,projxyz;

  // Print out data in XYZ format
  oss << sstr.atoms.size() << endl;
  oss << sstr.title << endl;
  for(uint iat=0;iat<sstr.atoms.size();iat++)
    for(int i=-1;i<=1;i++)
      for(int j=-1;j<=1;j++)
        for(int k=-1;k<=1;k++)
          {
            rat=sstr.atoms.at(iat).cpos+((double)i)*sstr.lattice(1)+((double)j)*sstr.lattice(2)+((double)k)*sstr.lattice(3);
            projxoo=scalar_product(rat,xoo)/axoo/axoo;
            projyoo=scalar_product(rat,yoo)/ayoo/ayoo;
            projzoo=scalar_product(rat,zoo)/azoo/azoo;
            projxyo=scalar_product(rat,xyo)/axyo/axyo;
            projxzo=scalar_product(rat,xzo)/axzo/axzo;
            projyzo=scalar_product(rat,yzo)/ayzo/ayzo;
            projxyz=scalar_product(rat,xyz)/axyz/axyz;
            //  if((projxoo>-0.5 && projxoo<=0.5))
            //  oss << projxoo << " " << projyoo << " " << projzoo << endl;
            if((projxoo>-0.5 && projxoo<=0.5) &&
               (projyoo>-0.5 && projyoo<=0.5) &&
               (projzoo>-0.5 && projzoo<=0.5) &&
               (projxyo>-0.5 && projxyo<=0.5) &&
               (projxzo>-0.5 && projxzo<=0.5) &&
               (projyzo>-0.5 && projyzo<=0.5) &&
               (projxyz>-0.5 && projxyz<=0.5))
              oss << sstr.atoms.at(iat).cleanname << " " << rat(1) << " " << rat(2) << " " << rat(3) << " " << endl;
          }
}

// ***************************************************************************
// PrintXYZInSphere
// ***************************************************************************
// This funtion prints out structural data inside a sphere with radius radius
// and 0,0,0 as origin, in an XYZ format
// Stefano Curtarolo
void PrintXYZInSphere(const xstructure& a, const double& radius, ostream& oss) {
  oss.setf(std::ios::fixed,std::ios::floatfield);
  oss.precision(10);
  xstructure sstr=a;
#define _SPHERE_CUT_OFF_ 5

  // Set scale to 1 so you don't need to rescale coordinates.
  sstr=ReScale(sstr,1.0);

  int atoms_inside=0;
  xvector<double> x(3);
  // Print out data in InSphere format
  for(int i=-_SPHERE_CUT_OFF_;i<_SPHERE_CUT_OFF_;i++) {
    for(int j=-_SPHERE_CUT_OFF_;j<_SPHERE_CUT_OFF_;j++) {
      for(int k=-_SPHERE_CUT_OFF_;k<_SPHERE_CUT_OFF_;k++) {
        for(uint iat=0;iat<sstr.atoms.size();iat++) {
          x(1)=(sstr.atoms.at(iat).cpos(1)+i*sstr.lattice(1,1)+j*sstr.lattice(2,1)+k*sstr.lattice(3,1));
          x(2)=(sstr.atoms.at(iat).cpos(2)+i*sstr.lattice(1,2)+j*sstr.lattice(2,2)+k*sstr.lattice(3,2));
          x(3)=(sstr.atoms.at(iat).cpos(3)+i*sstr.lattice(1,3)+j*sstr.lattice(2,3)+k*sstr.lattice(3,3));
          if(aurostd::modulus(x)<=radius) atoms_inside++;
        } // iat
      } // k
    } // j
  } // i

  oss << atoms_inside << endl;
  oss << sstr.title << endl;
  // Print out data in InSphere format
  for(int i=-_SPHERE_CUT_OFF_;i<_SPHERE_CUT_OFF_;i++) {
    for(int j=-_SPHERE_CUT_OFF_;j<_SPHERE_CUT_OFF_;j++) {
      for(int k=-_SPHERE_CUT_OFF_;k<_SPHERE_CUT_OFF_;k++) {
        for(uint iat=0;iat<sstr.atoms.size();iat++) {
          x(1)=(sstr.atoms.at(iat).cpos(1)+i*sstr.lattice(1,1)+j*sstr.lattice(2,1)+k*sstr.lattice(3,1));
          x(2)=(sstr.atoms.at(iat).cpos(2)+i*sstr.lattice(1,2)+j*sstr.lattice(2,2)+k*sstr.lattice(3,2));
          x(3)=(sstr.atoms.at(iat).cpos(3)+i*sstr.lattice(1,3)+j*sstr.lattice(2,3)+k*sstr.lattice(3,3));
          if(aurostd::modulus(x)<=radius)
            oss << sstr.atoms.at(iat).cleanname << " " << x(1) << " " << x(2) << " " << x(3) << " " << endl;
        } // iat
      } // k
    } // j
  } // i
}

// ***************************************************************************
// PrintXRAY ids_cmp
// ***************************************************************************
// This funtion prints out the xray scattering pattern
class ids_cmp{
public:
  int operator()(const vector<double>& a, const vector<double>& b)
  {return a[0]>b[0];} // Sorts in increasing order.
};

class hkl_cmp{
public:
  int operator()(const vector<int>& a, const vector<int>& b)
  {
    int na=a[0]*100+a[1]*10+a[2];
    int nb=b[0]*100+b[1]*10+b[2];
    return na>nb;
  }
};
//      int t=1,f=1;
//      if(a[0]<b[0]) {
//	return false;
//      }
//      else{
//	if(a[1]<b[1]) {
//	  return false;
//	}
//	else{
//	  if(a[2]<b[2]) {
//	    return false;
//	  }
//	  else{
//	    return true;
//	  }
//	}
//      }
//    }

// ***************************************************************************
// PrintXRAY CorrectionFactor
// ***************************************************************************
// This function contains any correction factors that need
// to be multiply the structure factors.  Theta should be in rads.
// Includes - Lorentz-Polarization
// Dane Morgan - Stefano Curtarolo
double CorrectionFactor(const double& th) {
  //  double th=theta*TWOPI/360.0; // rad->degrees
  double correction_factor=1;
  // Lorentz-polarization
  double lp;
  double c2t=cos(2.0*th);
  double st=sin(th);
  if(th>0 || th<0) {
    lp=(1.0+c2t*c2t)/(st*st);
  }
  else{
    lp=1;
  }
  correction_factor=correction_factor*lp;
  return lp;
}  

// ***************************************************************************
// PrintXRAY PrintXray
// ***************************************************************************
// Dane Morgan - Stefano Curtarolo
void PrintXray(const xstructure& str, const double& lambda, ostream& oss) {  
  oss.setf(std::ios::left, std::ios::adjustfield);
  oss.setf(std::ios::fixed, std::ios::floatfield);
  //  oss.setf(std::ios::fixed,std::ios::floatfield);
  int PREC_DEFAULT=4;
  oss.precision(PREC_DEFAULT);

  int num_atoms=str.atoms.size();
  //  vector<string> names=str.GetNames();
  vector<double> dist,sf;
  vector<double> scatt_fact(num_atoms,0.0);
  vector<double> mass(num_atoms,0.0);
  vector<double> twoB_vec(num_atoms,0.0);
  pflow::GetXray(str,dist,sf,lambda,scatt_fact,mass,twoB_vec);

  double tol=1.0E-5;
  int w1=4; // int
  int w2=12; // some doubles
  int w3=20; // Integrated intensities
  int tlen=sf.size();
  int len=Nint(std::pow((double) tlen,(double) 1.0/3.0));
  int kmx=(len-1)/2; // len should be odd.

  // Sort by theta (reverse sort by distance).
  // Define an id pointer to sort.
  vector<double> v(5);
  // pflow::matrix<double> ids(tlen,v);
  vector<vector<double> > ids(tlen,v);
  for(int i0=-kmx;i0<=kmx;i0++) {
    for(int i1=-kmx;i1<=kmx;i1++) {
      for(int i2=-kmx;i2<=kmx;i2++) {
        int ii0=i0+kmx;
        int ii1=i1+kmx;
        int ii2=i2+kmx;
        int id=ii2+ii1*len+ii0*len*len;
        ids[id][0]=dist[id];
        ids[id][1]=id;
        ids[id][2]=i0;
        ids[id][3]=i1;
        ids[id][4]=i2;
      } // i2
    } // i1
  } // i0

  sort(ids.begin(),ids.end(),ids_cmp());

  // Add corrections to all the amplitudes.
  // Get max amplitude for normalizing and percentages.
  double ampmax=1E-8;
  for(int i0=-kmx;i0<=kmx;i0++) {
    for(int i1=-kmx;i1<=kmx;i1++) {
      for(int i2=-kmx;i2<=kmx;i2++) {
        int ii0=i0+kmx;
        int ii1=i1+kmx;
        int ii2=i2+kmx;
        int id1=ii2+ii1*len+ii0*len*len;
        int id=(int) ids[id1][1];
        ii0=(int) ids[id1][2];
        ii1=(int) ids[id1][3];
        ii2=(int) ids[id1][4];
        double theta=0;
        if(dist[id]>0) {
          double term=lambda/(2.0*dist[id]);
          if(term<=1) {
            theta=asin(term);
            if(theta>tol) {
              sf[id]=sf[id]*CorrectionFactor(theta);
              if(sf[id]>ampmax) {ampmax=sf[id];}
            } // if theta>tol
          } // if term<=1
        } // if dist>0
      } // i2
    } // i1
  } // i0

  // Print out all data.
  oss << "Wavelength (Ang) = " << lambda << endl;
  oss << "Atom Name  ScattFact   Mass(amu)   B(Ang)(DW=exp(-B*sin(theta)^2/lambda^2))" << endl;
  for(uint iat=0;iat<(uint) num_atoms;iat++) {
    oss <<setw(4)<<iat+1<<" "<<setw(4)<<str.atoms.at(iat).cleanname << " ";
    if(str.atoms.at(iat).cleanname.length()>1) oss << " ";
    oss <<setw(w2)<<scatt_fact[iat]<<setw(w2)<<KILOGRAM2AMU*mass[iat]<<setw(w2)<<1E+20*twoB_vec[iat]/2.0<<endl;
  }
  oss << "******************** All data ********************" << endl;
  oss << "2*theta      Intensity            h    k    l    dist         keyword " << endl;
  for(int i0=-kmx;i0<=kmx;i0++) {
    for(int i1=-kmx;i1<=kmx;i1++) {
      for(int i2=-kmx;i2<=kmx;i2++) {
        int ii0=i0+kmx;
        int ii1=i1+kmx;
        int ii2=i2+kmx;
        int id1=ii2+ii1*len+ii0*len*len;
        int id=(int) ids[id1][1];
        ii0=(int) ids[id1][2];
        ii1=(int) ids[id1][3];
        ii2=(int) ids[id1][4];
        double theta=0;
        if(dist[id]>0) {
          double term=lambda/(2.0*dist[id]);
          if(term<=1) {
            theta=asin(term);
            theta=theta*360.0/TWOPI; // rad->degrees
            if(theta>tol) oss
                            <<setw(w2)<<2.0*theta<<" " // angle
                            <<setw(w3)<<sf[id]<<" " // sf
                            <<setw(w1)<<ii0<<" " // h
                            <<setw(w1)<<ii1<<" " // k
                            <<setw(w1)<<ii2<<" " // l
                            <<setw(w2)<<dist[id]<<" " // dist
                            << "SINGLE"
                            << endl;
          } // if term<=1
        } // if dist>0
      } // i2
    } // i1
  } // i0

  // Now group everything at the same distance together and only store one entry for each distance.
  // Choose hkl such that h is the largest, k the second, and l the third.
  // Multiply the sf by the number of degenerate points at that distance.
  // Get max integrated intensity for normalizations and percentages.
  double intmax=1E-8;
  double odist=dist[(int)ids[0][1]]; // Initialize odsit to first distance.
  double osf=sf[(int)ids[0][1]]; // Initialize osf to first distance.
  vector<vector<int> > hkl_list;
  vector<int> hkl(3);
  pflow::matrix<double> data;
  for(int i0=-kmx;i0<=kmx;i0++) {
    for(int i1=-kmx;i1<=kmx;i1++) {
      for(int i2=-kmx;i2<=kmx;i2++) {
        int ii0=i0+kmx;
        int ii1=i1+kmx;
        int ii2=i2+kmx;
        int id1=ii2+ii1*len+ii0*len*len;
        int id=(int) ids[id1][1];
        ii0=(int) ids[id1][2];
        ii1=(int) ids[id1][3];
        ii2=(int) ids[id1][4];
        double pdist=dist[id]; // Get present distance.
        // Create vector of all the hkl values with the same distance
        if(fabs(pdist-odist)<tol) { // Add present h,k,l to hkl_list.
          hkl[0]=ii0;		    
          hkl[1]=ii1;		    
          hkl[2]=ii2;		    
          hkl_list.push_back(hkl);
        }
        else{ // Store one hkl, dist, total sf, multiplicity of point in data vector and then reset hkl_list vector to new hkl.
          vector<double> datav(6);
          // Sort hkl
          sort(hkl_list.begin(),hkl_list.end(),hkl_cmp());
          datav[0]=(double) hkl_list[0][0];  
          datav[1]=(double) hkl_list[0][1];  
          datav[2]=(double) hkl_list[0][2];  
          datav[3]=odist;
          datav[4]=osf*hkl_list.size();
          datav[5]=hkl_list.size();
          data.push_back(datav);
          if(datav[4]>intmax) {intmax=datav[4];}
          vector<int> v(0);
          hkl_list= vector<vector<int> > (0,v);
          hkl[0]=ii0;		    
          hkl[1]=ii1;		    
          hkl[2]=ii2;		    
          hkl_list.push_back(hkl);
        }
        odist=pdist;
        osf=sf[id];
      } // i2
    } // i1
  } // i0

  // Output grouped data
  oss << "******************** Grouped data ********************" << endl;
  oss << "2*theta      IntIntensity         %ofMaxInt    h    k    l    dist         mult. correction    keyword " << endl;
  for(uint i=0;i<data.size();i++) {
    double theta=0;
    if(data[i][3]>0) {
      double term=lambda/(2.0*data[i][3]);
      if(term<=1) {
        theta=asin(term);
        theta=theta*360.0/TWOPI; // rad->degrees
        if(theta>tol) oss
                        <<setw(w2)<<2.0*theta<<" " // angle
                        <<setw(w3)<<setprecision(2)<<data[i][4]<<setprecision(PREC_DEFAULT)<<" " // sf
                        <<setw(w2)<<setprecision(2)<<100*data[i][4]/intmax<<setprecision(PREC_DEFAULT)<<" " // % max sf
                        <<setw(w1)<<(int)data[i][0]<<" " // h
                        <<setw(w1)<<(int)data[i][1]<<" " // k
                        <<setw(w1)<<(int)data[i][2]<<" " // l
                        <<setw(w2)<<data[i][3]<<" " // dist
                        <<setw(5)<<(int)data[i][5]<<" " // mult.
                        <<setw(w2)<<CorrectionFactor(theta*TWOPI/360.0)<<" " // correction.
                        << " GROUP"
                        << endl;
      } // if term<=1
    } // if dist>0
  } // for

  // Output data to plot
  oss << "******************** To Plot data ********************" << endl;
  oss << "2*theta      Amplitude    keyword " << endl;
  for(uint i=0;i<data.size();i++) {
    double theta=0;
    if(data[i][3]>0) {
      double term=lambda/(2.0*data[i][3]);
      if(term<=1) {
        theta=asin(term);
        theta=theta*360.0/TWOPI; // rad->degrees
        if(theta>tol) {
          // initial 0.
          oss <<setw(w2)<<2.0*theta<<" " // angle
              <<setw(w2)<<"0"<<" " // sf
              << "TOPLOT"
              << endl;
          // true value of sf/intmax.
          oss <<setw(w2)<<2.0*theta<<" " // angle
              <<setw(w2)<<setprecision(2)<<100*data[i][4]/intmax<<setprecision(PREC_DEFAULT)<<" " // sf
              << "TOPLOT"
              << endl;
          // final 0.
          oss <<setw(w2)<<2.0*theta<<" " // angle
              <<setw(w2)<<"0"<<" " // sf
              << "TOPLOT"
              << endl;
          // tpx
        } // if theta<tol
      } // if term<=1
    } // if dist>0
  } // for
} // end routine

#endif

// **************************************************************************
// *                                                                        *
// *             STEFANO CURTAROLO - Duke University 2003-2018              *
// *                                                                        *
// **************************************************************************

