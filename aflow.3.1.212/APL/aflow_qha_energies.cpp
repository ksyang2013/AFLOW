// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2015           *
// *                Aflow PINKU NATH - Duke University 2014-2018             *
// *                                                                         *
// ***************************************************************************
// Written by Pinku Nath
// pn49@duke.edu

#include "aflow_apl.h"
#include <iterator>
#define _isnegative(a) (a<MIN_EIGEN_TRESHOLD) ? true : false

/*
  This class storing all static energies, phonon dispersion curves and electronic curves from various distorted directories.
*/


namespace apl
{
  // ***************************************************************************************
  QH_ENERGIES::QH_ENERGIES( IPhononCalculator& pc, QHA_AFLOWIN_CREATOR& runeos, Logger& l): _pc(pc),_runeos(runeos),_logger(l)
  {
    _logger << "Calculating electronic and thermal energies to compute the thermodynamic properties"<< apl::endl; 
    //clear all memories before using
    clear();
  }
  // ***************************************************************************************
  QH_ENERGIES::~QH_ENERGIES()
  {
    this->clear();
  }
  // ***************************************************************************************
  void QH_ENERGIES::clear()
  {
    _is_magnetic=false;
    _ph_vols.clear();
    _ph_dirs.clear();
    _ele_vols.clear();
    _ele_dirs.clear();
    _pdos.clear();
    _edos.clear();
    _cedos.clear();
    _eo.clear();
    _tmp_dir="";
    _fermi_energies.clear();
    _pV.clear();
    _mag_cell.clear();
    _eqm_ele_dir_index=-11;
    _index_imag_freq.clear();
    _atomic_species.clear();
  }
  // ***************************************************************************************
  //polupating xtracture
  void QH_ENERGIES::get_xtracture(const xstructure& xs)
  {
    _logger<<"Populationg xtracture in QH_ENERGIES" <<apl::endl;
    xstructure xstr(xs);
    for(uint i=0; i!=xstr.atoms.size(); i++){
      _atomic_species.push_back(xstr.atoms[i].name);
    }
  }
  // ***************************************************************************************
  //get phonon, electronin energies
  bool QH_ENERGIES::get_qha_energies()
  {
   //phonon directory names
    _ph_dirs  =  _runeos.get_ph_dir_names();
   //static energies directory names
    _ele_dirs =  _runeos.get_eos_dir_names();
   //volumes of primotive cessa
    _ele_vols =  _runeos.get_eos_volumes();
   //index of relaxed volume
    _eqm_ele_dir_index=_runeos.get_zero_index();
   //polulate pdos
    if(!get_pdos()) return false;
   //populate static energies
    if(!getE0K())   return false;
    //check compound is magnetic
    if(_is_magnetic) _logger << "Compound is magnetic and skipping electronic dos calculation "<< apl::endl;
    //populate electronic dos
    if(!get_edos()) return false;
    //check sizes of all variables
    if(!size_check()) return false;
    //List of configurations with imaginary frequencies and remove them from list
    get_imaginary_freq_index();
    //correct noise from electronic bands
    //get_electronic_corrected_energies();
    //write energies
    if(!write_energies())return false;
    //if imaginary frequengy exist print message
    print_imaginary_freq_msg();
    //remove imaginary frequency directory from list if exist
    remove_imaginary_freqs();
    //check sizes of all variables again
    if(!size_check()) return false;
    //correct noise from electronic bands
    //get_electronic_corrected_energies();
    return true;
  }
  // ***************************************************************************************
  //get phonon, electronin energies for SCQHA and QHA3P methods
  bool QH_ENERGIES::get_scqha_energies()
  {
    _ele_dirs =  _runeos.get_eos_dir_names();
    _ele_vols =  _runeos.get_eos_volumes();
    _eqm_ele_dir_index=_runeos.get_zero_index();
    if(!getE0K())   return false;
    if(_is_magnetic) _logger << "Compound is magnetic and skipping electronic dos calculation "<< apl::endl;
    if(!get_edos()) return false;
    if(!size_scqha_check()) return false;
    //get_electronic_corrected_energies();
    if(!write_energies())return false;
    return true;
  }
  // ***************************************************************************************
  //check sizes for all vectors
  bool QH_ENERGIES::size_check()
  {
    if(_pdos.size()==0)
      {
	_logger << apl::error << "_pdos.size()==0 "<< apl::endl;
	return false;
      }
    if(_pdos.size()!=_edos.size())
      {
	_logger << apl::error << "_pdos.size()!=_edos.size() "<< apl::endl;
	return false;
      }
    if(_eo.size()==0)
      {
	_logger << apl::error << "_eo.size()==0 "<< apl::endl;
	return false;
      }
    if(_eo.size()!=_fermi_energies.size())
      {
	_logger << apl::error << "_eo.size()!=_fermi_energies.size() "<< apl::endl;
	return false;
      }
    if(_eo.size()!=_pV.size())
      {
	_logger << apl::error << "_eo.size()!=_pV.size() "<< apl::endl;
	return false;
      }
    if(_eo.size()!=_mag_cell.size())
      {
	_logger << apl::error << "_eo.size()!=_mag_cell.size() "<< apl::endl;
	return false;
      }
    if(_eo.size()!=_ele_dirs.size())
      {
	_logger << apl::error << "_eo.size()!=_ele_dirs.size() "<< apl::endl;
	return false;
      }
    if(_eo.size()!=_ele_vols.size())
      {
	_logger << apl::error << "_eo.size()!=_ele_vols.size() "<< apl::endl;
	return false;
      }
    if(_atomic_species.size()==0)
      {
	_logger << apl::error << "_atomic_species()==0 "<< apl::endl;
	return false;
      }
    return true;
  }
  // ***************************************************************************************
  //check sizes for all vectors for SCQHA and QHA3P methods
  bool QH_ENERGIES::size_scqha_check()
  {
    if(_eo.size()==0)
      {
	_logger << apl::error << "_eo.size()==0 "<< apl::endl;
	return false;
      }
    if(_eo.size()!=_fermi_energies.size())
      {
	_logger << apl::error << "_eo.size()!=_fermi_energies.size() "<< apl::endl;
	return false;
      }
    if(_eo.size()!=_pV.size())
      {
	_logger << apl::error << "_eo.size()!=_pV.size() "<< apl::endl;
	return false;
      }
    if(_eo.size()!=_mag_cell.size())
      {
	_logger << apl::error << "_eo.size()!=_mag_cell.size() "<< apl::endl;
	return false;
      }
    if(_eo.size()!=_ele_dirs.size())
      {
	_logger << apl::error << "_eo.size()!=_ele_dirs.size() "<< apl::endl;
	return false;
      }
    if(_eo.size()!=_ele_vols.size())
      {
	_logger << apl::error << "_eo.size()!=_ele_vols.size() "<< apl::endl;
	return false;
      }
    return true;
  }
  // ***************************************************************************************
  //get directory indices if imaginary frequency exist
  void QH_ENERGIES::print_imaginary_freq_msg()
  {
    uint j=0;
    if(_index_imag_freq.size()>0)
      {
	for(uint i=0; i!=_ele_dirs.size(); i++)
	  {
	    if(i==_index_imag_freq[j])
	      {
		_logger << "Imaginary freqs corresponding to the phonon dir = "<<_ele_dirs[i]<<" and excluded." << apl::endl;
		j++;
	      }
	  }
      }
  }
  // ***************************************************************************************
  //Check imaginary frequencies
  void QH_ENERGIES::get_imaginary_freq_index()
  {
    for(uint i=0; i!=_pdos.size(); i++)
      {
	for(uint j=0; j!=_pdos[i].size(); j++)
	  {
	    if(_isnegative(_pdos[i][j][0]))
	      {
		if(_isnegative(_pdos[i][j][0]))
		  {
		    _index_imag_freq.push_back(i);
		    break;
		  }
		_pdos[i][j][0]=0.0;
	      }
	  }
      }
  }
  // ***************************************************************************************
  //remove imaginary frequency directory from the list
  bool QH_ENERGIES::remove_imaginary_freqs()
  {
    if(_index_imag_freq.size()>0)
      {
        vector<vector<vector<double> > > tmp_pdos = _pdos;
        vector<vector<vector<double> > > tmp_edos = _edos;
        vector<double> tmp_eo=_eo;
        vector<double> tmp_fermi_energies=_fermi_energies;
        vector<double> tmp_pV=_pV;
        vector<double> tmp_mag_cell=_mag_cell;
        vector<string> tmp_ele_dirs=_ele_dirs;
        vector<double> tmp_ele_vols=_ele_vols;
        _pdos.clear();
        _edos.clear();
        _eo.clear();
        _fermi_energies.clear();
        _pV.clear();
        _mag_cell.clear();
        _ele_dirs.clear();
        _ele_vols.clear();
        //
        uint j=0;
        for(uint i=0; i!=tmp_pdos.size(); i++)
          {
            if(i!=_index_imag_freq[j])
              {
                _pdos.push_back(tmp_pdos[i]);
                _edos.push_back(tmp_edos[i]);
                _eo.push_back(tmp_eo[i]);
                _fermi_energies.push_back(tmp_fermi_energies[i]);
                _pV.push_back(tmp_pV[i]);
                _mag_cell.push_back(tmp_mag_cell[i]);
                _ele_dirs.push_back(tmp_ele_dirs[i]);
                _ele_vols.push_back(tmp_ele_vols[i]);
              }else j++;
          }

        if(!write_imag_freq_corrected_energies(tmp_eo)) return false;
 
        //clear temporary variables
        tmp_pdos.clear();
        tmp_edos.clear();
        tmp_eo.clear();
        tmp_fermi_energies.clear();
        tmp_pV.clear();
        tmp_mag_cell.clear();
        tmp_ele_dirs.clear();
        tmp_ele_vols.clear();
        //
	//remove possible noise from electronic bands
	//get_electronic_corrected_energies();
      }
    return true;
  }
  // ***************************************************************************************
  //write static energies and volumes
  bool QH_ENERGIES::write_energies()
  {
    _logger << "Writing aflow.qha.static_energies.out file "<< apl::endl;
    stringstream out;
    out << std::setiosflags(std::ios::fixed|std::ios::showpoint|std::ios::right);
    out << setprecision(6);
    std::string STAR = std::string(50, '*');
    out<<"[AFLOW] "<<STAR<<"\n";
    out<<"#"<<"Total phonon dos found = "<<_pdos.size() << '\n';
    out<<"#"<<"Total electronic dos found = "<<_edos.size() << '\n';
    //out<<"#"<<"Electronic Gruneisen = "<< get_electronic_gruneisen() << '\n';
    out<<"[AFLOW] "<<STAR<<"\n";
    out<< "[AFLOW_QHA_STATIC_ENERGIES]START" <<"\n";
    out<<"#"<<setw(30)<<"dir"<<setw(15)<<"E0 (eV/Cell)"<<setw(15)<<"Vol (A^3)"
       <<setw(15)<<"Ef (eV)"<<setw(20)<<"pV (eV/Cell)"<<setw(20)<<"magcell (mu/cell)"<<setw(15)<<"freq_type" <<'\n';

    uint j=0;
    for(uint i=0; i!=_eo.size(); i++)
      {
	out<<setw(30) << _ele_dirs[i];
        out<<setw(15)<<_eo[i]<<setw(15)<<_ele_vols[i]
           <<setw(15)<<_fermi_energies[i]
	  //<<setw(15)<<_edosATfermi[i]
           <<setw(20)<<_pV[i]<<setw(15)<<_mag_cell[i];

	if(_index_imag_freq.size()>0)
	  {
	    if(i==_index_imag_freq[j])
	      {
		out<<setw(15)<<"-";
		j++;
	      }else out<<setw(15)<<"+";
	  }else out<<setw(15)<<"+";

        out<<'\n';
      }
    out<< "[AFLOW_QHA_STATIC_ENERGIES]END" <<"\n";
    out<<"[AFLOW] "<<STAR<<"\n";
    string eos_out =  "aflow.qha.static_energies.out";
    if(!aurostd::stringstream2file(out, eos_out, "WRITE")) {
      throw APLRuntimeError("Cannot write aflow.qha.static_energies.out");
    }
    aurostd::StringstreamClean(out);
    return true;
  }
  // ***************************************************************************************
  //write list of directories with positive frequencies
  bool QH_ENERGIES::write_imag_freq_corrected_energies(const vector<double> &tmp_eo)
  {
    if(_index_imag_freq.size()>0){
      _logger << "Writing aflow.apl.static_corrected.out file "<< apl::endl;

      stringstream out;
      string outfile="aflow.apl.static_corrected.out";     
      out << std::setprecision(6) << std::fixed;
      std::string STAR = std::string(50, '*');
      out<<"[AFLOW] "<<STAR<<"\n";
      out<< "[APL_STATIC_ENERGIES]START" <<"\n";
      out<<"#"<<setw(30)<<"dir"<<setw(15)<<"E0 (eV/Cell)"<<setw(15)<<"Vol (A^3)"
	 <<setw(15)<<"E_f (eV)"<<setw(20)<<"pV (eV/Cell)"<<setw(20)<<"magcell (mu/cell)"<<setw(15)<<"freq_type" <<'\n';

      uint j=0;
      uint k=0;
      for(uint i=0; i!=tmp_eo.size(); i++)
	{
	  if(i!=_index_imag_freq[j])
	    {
	      out<<setw(30) << _ele_dirs[k];
	      out<<setw(15)<<_eo[i]<<setw(15)<<_ele_vols[k]
		 <<setw(15)<<_fermi_energies[k]
		 <<setw(20)<<_pV[i]<<setw(15)<<_mag_cell[k];
	      out<<setw(15)<<"+"<<'\n';
	      k++;
	    }else
	    {
	      j++;
	    }
	}
      out<< "[APL_STATIC_ENERGIES]END" <<"\n";
      out<<"[AFLOW] "<<STAR<<"\n";
      if(!aurostd::stringstream2file(out, outfile, "WRITE")) {
	throw APLRuntimeError("Cannot write aflow.apl.static_corrected.out");
      }
      aurostd::StringstreamClean(out);
    }
    return true;
  }
  // ***************************************************************************************
  //read pdos from various distorted directories
  bool QH_ENERGIES::get_pdos()
  {
    string file="";
    _logger <<"Reading phonon dos files " << apl::endl;
    if(_ele_dirs.size()!=_ph_dirs.size()+1)
      {
	_logger << apl::error <<" _ele_dirs.size()!=_ph_dirs.size()+1 " << apl::endl;
	return false;
      }
    for(uint i=0; i!=_ph_dirs.size(); i++)
      {
	file="";
	//adding equilibrium directory
	if((int)i==_eqm_ele_dir_index){
	  file="PDOS";
	  if(!exists_test0(file) && !aurostd::EFileExist(file)) return false;
	  _pdos.push_back(get_pdos(file));
	}
	file=_tmp_dir+"/PDOS."+_ph_dirs[i];
	if((!exists_test0(file) && !aurostd::EFileExist(file))) return false;
	_pdos.push_back(get_pdos(file));
	
      }
    return true;
  }
  // ***************************************************************************************
  //read edos from various distorted directories
  bool QH_ENERGIES::get_edos()
  {
    _logger<<"Reading electronic dos files " << apl::endl;

    for(uint i=0; i!=_ele_dirs.size(); i++)
      {
	string file=_ele_dirs[i]+string("/")+string("DOSCAR.static");
	if(!exists_test0(file) && !aurostd::EFileExist(file)) return false;
	double fermi=0.0;
	_edos.push_back(get_edos(file, fermi));
	_fermi_energies.push_back(fermi);
      }
    return true;
  }
  // ***************************************************************************************
  //read static energies from various distorted directories
  bool QH_ENERGIES::getE0K()
  {
    _logger <<"Reading static energies " << apl::endl;
    for(uint i=0; i!=_ele_dirs.size(); i++)
      {
	string file=_ele_dirs[i]+"/aflow.qmvasp.out";
	if(!exists_test0(file) && !aurostd::EFileExist(file)) return false;
	double pv=0.;
	double magcell=0.;
	_eo.push_back(getE0K(file, pv, magcell));
	_pV.push_back(pv);
	_mag_cell.push_back(magcell);
	if(std::abs(magcell)>0) _is_magnetic=true;
      }
    return true;
  }
  // ***************************************************************************************
  //populate pdos
  vector<vector<double> > QH_ENERGIES::get_pdos(const string file)
  {
    if(!exists_test0(file) && !aurostd::EFileExist(file)) {
      throw apl::APLRuntimeError("QH_ENERGIES:: Missing file: "+file);
    }
    vector<string> vlines;
    aurostd::efile2vectorstring(file, vlines);
    if (!vlines.size()) {
      throw apl::APLRuntimeError("QH_ENERGIES:: Missing file: "+file);
    }
    uint line_count = 0;
    string line;
    vector<vector<double> > return_vec; return_vec.clear();

    while (line_count < vlines.size()){
      line = vlines[line_count++];
      if(line=="")continue;
      if(line[0]=='#')continue;
      vector<string> vstr=split<string>(line);
      if(vstr.size()!=4)
	{
	  _logger << apl::error << file<<" Wrong format." << apl::endl; exit(0);
	}
      vector<double>  tmp(2, 0.0);
      tmp[0]=atof(vstr[0].c_str());
      tmp[1]=atof(vstr[3].c_str());
      return_vec.push_back(tmp);
    }
    return return_vec;
  }
  // ***************************************************************************************
  void QH_ENERGIES::get_tmp_dir_name(const string dir)
  {
    _tmp_dir=dir;
  }
  // ***************************************************************************************
  //populate edos
  vector<vector<double> > QH_ENERGIES::get_edos(const string file, double &fermi)
  {
    if (!exists_test0(file) && !aurostd::EFileExist(file)) {
      throw apl::APLRuntimeError("QH_ENERGIES:: Missing file: "+file);
    }
    vector<string> vlines;
    aurostd::efile2vectorstring(file, vlines);
    if (!vlines.size()) {
      throw apl::APLRuntimeError("QH_ENERGIES:: Missing file: "+file);
    }

    fermi=0.0;

    string line;
    uint line_count = 0;
    uint LC=0;
    vector<vector<double> > tmp_edos; tmp_edos.clear();
    while (line_count < vlines.size()){
      line = vlines[line_count++];
      if(line[0]=='#')continue;
      if(line=="")continue;
      LC++;
      vector<string> v=split<string>(line);
      if(LC==6)fermi=atof(v[v.size()-2].c_str());
      else if(LC>6)
	{
	  if(v.size()!=3)break;
	  vector<double> tmp_d(2, 0.0);
	  tmp_d[0]=atof(v[0].c_str());
	  tmp_d[1]=atof(v[1].c_str());
	  tmp_edos.push_back(tmp_d);
	}
    }
    return tmp_edos;
  }
  // ***************************************************************************************
  //read 0K energies with magnetic properties
  double QH_ENERGIES::getE0K(const string file, double &pv, double &mag)
  {
    if (!exists_test0(file) && !aurostd::EFileExist(file)) {
      throw apl::APLRuntimeError("QH_ENERGIES:: Missing file: "+file);
    }
    vector<string> vlines;
    aurostd::efile2vectorstring(file, vlines);
    if (!vlines.size()) {
      throw apl::APLRuntimeError("QH_ENERGIES:: Missing file: "+file);
    }

    pv=0.0;
    mag=0.0;
    vector<string> WORD;
    vector<string> PV;
    vector<string> mag_cell;
    string line;
    uint line_count = 0;

    while (line_count < vlines.size()) {
      line = vlines[line_count++];
      if(line=="")continue;
      if((line[0]=='#') && (line[0]=='/'))continue;
      if(line.find("H_cell") != std::string::npos)WORD.push_back(line);
      else if(line.find("PV_cell") != std::string::npos)PV.push_back(line);
      else if(line.find("mag_cell") != std::string::npos)mag_cell.push_back(line);
    }
    //depends on aflow.qmvasp.out format
    string s=WORD[WORD.size()-1];
    vector<string> vec=split<string>(s);
    string str=vec[0];
    str.erase (str.begin(), str.begin()+7);
    vector<double> vec1=split<double>(str);

    string s1=PV[PV.size()-1];
    vector<string> vecPV=split<string>(s1);
    string str1=vecPV[0];
    str1.erase (str1.begin(), str1.begin()+8);
    vector<double> vecPV1=split<double>(str1);
    pv=vecPV1[0];

    string s2=mag_cell[mag_cell.size()-1];
    vector<string> vecMag=split<string>(s2);
    string str2=vecMag[0];
    str2.erase (str2.begin(), str2.begin()+9);
    vector<double> tmpMag=split<double>(str2);
    mag=tmpMag[0];
    return vec1[0];
  }
  // ***************************************************************************************
  vector<uint>  QH_ENERGIES::get_index_imag_freq()
  { 
    return _index_imag_freq;
  }
  // ***************************************************************************************
  vector<double>  QH_ENERGIES::get_mag_cell()
  {
    return _mag_cell;
  }
  // ***************************************************************************************
  vector<double>  QH_ENERGIES::get_pV()
  {
    return _pV;
  }
  // ***************************************************************************************
  vector<double>  QH_ENERGIES::get_fermi_energies()
  {
    return _fermi_energies;
  }
  // ***************************************************************************************
  vector<double>  QH_ENERGIES::get_eo()
  {
    return _eo;
  }
  // ***************************************************************************************
  vector<vector<vector<double> > >  QH_ENERGIES::get_edos_data()
  {
    return _edos;
  }
  // ***************************************************************************************
  vector<vector<vector<double> > >  QH_ENERGIES::get_cedos_data()
  {
    //return _cedos;
    return _edos;
  }
  // ***************************************************************************************
  vector<vector<vector<double> > >  QH_ENERGIES::get_pdos_data()
  {
    return _pdos;
  }
  // ***************************************************************************************
  vector<double>  QH_ENERGIES::get_ele_vols()
  {
    return _ele_vols;
  }
  // ***************************************************************************************
  bool  QH_ENERGIES::get_is_magnetic()
  {
    return _is_magnetic;
  }
  // ***************************************************************************************
  int  QH_ENERGIES::get_eqm_ele_dir_index()
  {
    return _eqm_ele_dir_index;
  }
  // ***************************************************************************************
  vector<string> QH_ENERGIES::get_atomic_species()
  {
    return _atomic_species;
  }
  // ***************************************************************************************
  template<typename T>
  std::vector<T> QH_ENERGIES::split(const std::string& line)
  {
    std::istringstream is(line);
    return std::vector<T>(std::istream_iterator<T>(is), std::istream_iterator<T>());
  }
  // ***************************************************************************************
  //vector sorting
  template<class T>
  vector<uint> QH_ENERGIES::sorted_order (const vector<T> & arr)
  {
    uint n=arr.size();
    vector<uint> idx(n, 0);

    for (uint i=0; i<n; i++)
      {
	idx[i] = i;
      }

    for (uint i=0; i<n; i++)
      {
	for (uint j=i+1; j<n; j++)
	  {
	    if (arr[idx[i]] > arr[idx[j]])
	      {
		std::swap (idx[i], idx[j]);
	      }
	  }
      }

    return idx;
  }
  // ***************************************************************************************
  //remove possible noise from electronic bands
  void  QH_ENERGIES::get_electronic_corrected_energies()
  {
    _cedos=_edos;
    //fermi centering
    for(uint i=0; i!=_edos.size(); i++){
      double dE =abs(_edos[i][1][0]-_edos[i][0][0]);
      //if deltaE is very small then take the average
      if(_iszero(dE)){
	dE=0.0;
	for(uint j=0; j<20; j++)dE+=_edos[i][j+1][0]-_edos[i][j][0];
	dE/=20.0;
      }
      vector<double> dEv;dEv.clear();
      vector<uint> index; index.clear();
      for(uint j=0; j!=_edos[i].size(); j++){
	double delE=(_edos[i][j][0]-_fermi_energies[i]);
	if(std::abs(delE)>dE)continue;
	dEv.push_back(abs(delE));
	index.push_back(j);
      }
      double correction=0.0;
      vector<uint> index_sorted = sorted_order(dEv);
      correction=(_edos[i][index[index_sorted[0]]][0]- _fermi_energies[i]);
      //edos energies are shifting towards fremi energy
      //it is less than 2-3 meV shift to make the edos integration noise free
      for(uint j=0; j!=_edos[i].size(); j++){_cedos[i][j][0]= _edos[i][j][0]-correction;}
    }
  }
  // ***************************************************************************************
  //file existence check
  bool QH_ENERGIES::exists_test0 (const std::string& name)
  {
    ifstream f(name.c_str());
    if (f.good()) {
      f.close();
      return true;
    } else {
      f.close();
      return false;
    }
  }
  // ***************************************************************************************
}
