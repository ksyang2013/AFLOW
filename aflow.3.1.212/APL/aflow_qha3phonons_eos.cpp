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

#if GCC_VERSION >= 40400   // added two zeros
#define AFLOW_APL_MULTITHREADS_ENABLE 1
#include <thread>
#else
#warning "The multithread parts of APL will be not included, since they need gcc 4.4 and higher (C++0x support)."
#endif
/*
  It is QHA3P implementation. Details methodology is explained in QHA3P paper.
  In this method the thermodynamic properties are calculated by using
  minimization of total free energy function, similar to QHA method. 
*/

namespace apl
{
  // ***************************************************************************************
  // constructor
  QHA3POINTS::QHA3POINTS(SCQHA_QHA3P& scqha, QH_ENERGIES& qh_energies, Logger& l):_scqha(scqha), _qh_energies(qh_energies), _logger(l)
  {
    clear();
    _logger<<"QHA3P calculation is started, "<<apl::endl;
  }
  // ***************************************************************************************
  //destructor
  QHA3POINTS::~QHA3POINTS()
  {
    this->clear();
  }
  // ***************************************************************************************
  //clear
  void QHA3POINTS::clear()
  {
    _nBranches=0;
    _pext=0.0;
    _eo.clear();
    _ele_vols.clear();
    _freq0.clear();
    _freqM.clear();
    _freqP.clear();
    _d1fdv1.clear();
    _d2fdv2.clear();
    _scqha_volumes.clear();
    _weights.clear();
    _edos.clear();
    _pV.clear();
    _fermi_energies.clear();
    _ep_freqs.clear();
    _include_ele=false;
    _is_magnetic=false;
    _Eeq=0.0, _Beq=0.0, _Veq=0.0, _Bp=0.0;
    _TF.clear();
    _atomic_species.clear();
  }
  // ***************************************************************************************
  //populating external veriables
  bool QHA3POINTS::import_variables()
  {
    _logger<<"QHA3P Importing variables "<<apl::endl;
    _eo=_qh_energies.get_eo();
    _ele_vols=_qh_energies.get_ele_vols();
    _pV=_qh_energies.get_pV();
    _freq0=_scqha.get_freqs_mesh();
    _freqM=_scqha.get_freqs_meshM();
    _freqP=_scqha.get_freqs_meshP();
    _scqha_volumes=_scqha.get_qha_gpvol();
    _weights=_scqha.get_weights();
    _edos=_qh_energies.get_cedos_data();
    _fermi_energies=_qh_energies.get_fermi_energies();
    _is_magnetic=_qh_energies.get_is_magnetic();
    _atomic_species=_qh_energies.get_atomic_species();
    return check_size();
  }
  // ***************************************************
  //check size of all vectors
  bool QHA3POINTS::check_size()
  {
    if(_scqha_volumes.size()!=3)
      {
	_logger << apl::error << "QHA3POINTS::check_size() _scqha_volumes.size()!=3 "<< apl::endl;
	return false;
      }
    if(_eo.size()==0)
      {
	_logger << apl::error << "QHA3POINTS::check_size() _eo.size()==0 "<< apl::endl;
	return false;
      }
    if(_ele_vols.size()==0)
      {
	_logger << apl::error << "QHA3POINTS::check_size() _ele_vols.size()==0 "<< apl::endl;
	return false;
      }
    if(_freq0.size()==0)
      {
	_logger << apl::error << "QHA3POINTS::check_size() _freq0.size()==0 "<< apl::endl;
	return false;
      }
    if(_freqM.size()==0)
      {
	_logger << apl::error << "QHA3POINTS::check_size() _freqM.size()==0 "<< apl::endl;
	return false;
      }
    if(_freqP.size()==0)
      {
	_logger << apl::error << "QHA3POINTS::check_size() _freqP.size()==0 "<< apl::endl;
	return false;
      }
    if(_weights.size()==0)
      {
	_logger << apl::error << "QHA3POINTS::check_size() _weights.size()==0 "<< apl::endl;
	return false;
      }
    if(_edos.size()==0)
      {
	_logger << apl::error << "QHA3POINTS::check_size() _edos.size()==0 "<< apl::endl;
	return false;
      }
    if(_fermi_energies.size()==0)
      {
	_logger << apl::error <<"QHA3POINTS::check_size() _fermi_energies.size()==0 "<< apl::endl;
	return false;
      }
    if(_atomic_species.size()==0)
      {
	_logger << apl::error << "QHA3POINTS::check_size() _atomic_species()==0 "<< apl::endl;
	return false;
      }
    _nBranches=_freq0[0].rows;
    return true;
  }
  // ***************************************************************************************
  //Energy volume calling fitting function
  void QHA3POINTS::fitting()
  {
    _logger<<"SCQHA Fitting E-V data"<<apl::endl;

    //energy values
    xvector<double> E(_eo.size(), 1);
    //volume values
    xvector<double> V(_eo.size(), 1);

    for(uint i=0; i!=_eo.size(); i++)
      {
	E[i+1]=_eo[i];
	V[i+1]=_ele_vols[i];
      }
    //perform fitting
    md_lsquares_call(V,E);
  }
  // ***************************************************************************************
  //populate energy and volume and perform both linear and non-linear fit
  void QHA3POINTS::md_lsquares_call(const xvector<double> &V, const xvector<double> &E)
  {
    md_lsquares mdfit;
    mdfit.clear();
    for(int i=1; i<=V.rows; i++)
      {
        mdfit.Xdata.push_back(V[i]);
        mdfit.Ydata.push_back(E[i]);
      }
   // it calls both linear and nonlinear fit functions
    mdfit.cubic_polynomial_fit(); 
    _Veq = mdfit.nleqmV0;
    _Eeq = mdfit.nleqmE0;
    _Beq = mdfit.nleqmB0;
    _Bp  = mdfit.nleqmBp;
    //more refinements in fitting
    xvector<double> guess(4,1);
    guess[1]=_Eeq;
    guess[2]=_Beq;
    guess[3]=_Veq;
    guess[4]=_Bp;
    xvector<double> out(4,1);
    more_refinement(E, V, guess, out);
  }
  // ***************************************************************************************
  //doing more refinments after fitting to get better fitting parameters
  bool QHA3POINTS::more_refinement(const xvector<double> &E, const xvector<double> &V,
				   xvector<double> &guess, xvector<double> &out)
  {
    apl::aflowFITTING fit;

    if(!fit.birch_murnaghan_fitting (E, V, guess, out))return false;
    _Eeq  =out[1];
    _Beq  =out[2];
    _Veq  =out[3];
    _Bp   =out[4];
    fit.clear();
    return true;
  }
  // ***************************************************************************************
  //Calculate Taylor coefficients
  void QHA3POINTS::calculate_freq_derivative()
  {
    xvector<double> d(_nBranches,1);
    _d1fdv1.clear(); _d2fdv2.clear();
    _d1fdv1.resize(_freq0.size(),d);
    _d2fdv2.resize(_freq0.size(),d);

#ifdef AFLOW_APL_MULTITHREADS_ENABLE
    // Get the number of CPUS
    int ncpus = sysconf(_SC_NPROCESSORS_ONLN);// AFLOW_MachineNCPUs;
//    int qpointsPerCPU = _freq0.size() / ncpus;  OBSOLETE ME180801
    // Show info 
    if( ncpus == 1 )
      _logger.initProgressBar("Calculating taylor coefficients for QHA3POINTS");
    else
      _logger.initProgressBar("Calculating taylor coefficients for QHA3POINTS  (" + stringify(ncpus) + " threads)");

    // Distribute the calculation
    int startIndex, endIndex;
    std::vector< std::thread* > threads;
    vector<vector<int> > thread_dist = getThreadDistribution((int) _freq0.size(), ncpus);
    for (int icpu = 0; icpu < ncpus; icpu++) {
      startIndex = thread_dist[icpu][0];
      endIndex = thread_dist[icpu][1];
      threads.push_back( new std::thread(&QHA3POINTS::calculate_derivative,this,startIndex,endIndex) );
    }

/* OBSOLETE ME180801
    for(int icpu = 0; icpu < ncpus; icpu++) {
      startIndex = icpu * qpointsPerCPU;
      endIndex = startIndex + qpointsPerCPU;
      if( ( (uint)endIndex > _freq0.size() ) ||
          ( ( icpu == ncpus-1 ) && ( (uint)endIndex < _freq0.size() ) ) )
        endIndex = _freq0.size();
      threads.push_back( new std::thread(&QHA3POINTS::calculate_derivative,this,startIndex,endIndex) );
    }
*/

    // Wait to finish all threads here!
    for(uint i = 0; i < threads.size(); i++) {
      threads[i]->join();
      delete threads[i];
    }
    // Done
    _logger.finishProgressBar();

#else
    _logger<<"Calculating taylor coefficients for QHA3POINTS without threads "<<apl::endl;
    calculate_derivative(0, _freq0.size());
#endif
  }
  // ***************************************************************************************
  //Calculate Taylor expansion in threads
  void QHA3POINTS::calculate_derivative(int startIndex, int endIndex)
  {
    double dv=_scqha_volumes[0]-_scqha_volumes[2];
    double dv2=0.25*dv*dv; 

    for(int i=startIndex; i<endIndex; i++){
      for(uint j=1; j<=_nBranches; j++)
	{
          //1st derivative
	  _d1fdv1[i][j]= (_freqM[i][j] - _freqP[i][j])/dv;     
          //2nd derivative
	  _d2fdv2[i][j]= (_freqM[i][j] + _freqP[i][j] - 2.0* _freq0[i][j])/dv2;  
	}
    }
  }
  // ***************************************************************************************
  //Calculate thermodynamic properties in a given temperature range
  void QHA3POINTS::qha3pts_temperature_loop(double Tmin, double Tmax, double delta_T, ThermalPropertiesCalculator &e)
  {
    _logger<<"Writing aflow.qha3P.thermo.out file"<<apl::endl;
    if(!_include_ele)_logger<<"Writing aflow.qha3P.FVT.out file"<<apl::endl;

    std::string STAR80 = std::string(80, '*');
    std::string STAR130 = std::string(130, '*');
    std::string STAR180 = std::string(180, '*');
    std::string STAR150 = std::string(150, '*');

    //calculate Taylor coefficients
    calculate_freq_derivative();
    //Energy and volume variables
    xvector<double> E(_eo.size(), 1);
    xvector<double> V(_eo.size(), 1);

    //populate energy and volume
    for(uint i=0; i!=_eo.size(); i++)
      {
        E[i+1]=_eo[i];
        V[i+1]=_ele_vols[i];
      }
    //printing variables
    stringstream outfile;
    outfile << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
    outfile << setprecision(4);

    stringstream osfvt;
    bool includePV=false;

    //check expernal pressure is applied or not
    for(uint i=0; i!=_pV.size(); i++){
      if(!_iszero(_pV[i]))includePV=true;
    }

    //variable to store extrapolated frequencies
    _ep_freqs.resize(_ele_vols.size());

    //performing frequency extrapolation 
    for(uint v=0; v!=_ele_vols.size(); v++){
      //initializing
      _ep_freqs[v]=_freq0;

      for(uint i=0; i!=_freq0.size(); i++){
	for(int j=1; j<=_freq0[i].rows; j++){
	  if(_freq0[i][j]>0.0){
	    _ep_freqs[v][i][j] = _freq0[i][j]+ _d1fdv1[i][j]*(_ele_vols[v]-_scqha_volumes[1]) + 0.5*_d2fdv2[i][j]*pow((_ele_vols[v]-_scqha_volumes[1]), 2.0);
	  }else{
	    _ep_freqs[v][i][j]=0.0;
	  }
	}
      }
    }

    //start EOS calculations

    outfile<<"[AFLOW] "<<STAR130<<"\n";
    outfile<<"# U     => vibrational internal energy\n";
    outfile<<"# F     => vibrational free energy\n";
    outfile<<"# S     => vibrational energy \n";
    outfile<<"# Fqh   => QHA-3pts free energy \n";
    outfile<<"# V     => volume \n";
    outfile<<"# B     => Bulk Modulus \n";
    outfile<<"# Bp    => Pressure derivative Bulk Modulus \n";
    outfile<<"# gamma => Gruneisen Parameter \n";
    outfile<<"# alpha => Volume expansion cofficient x 10^6 \n";
    //outfile<<"# Bpp => Pressure derivative of Bp  \n";
    if(_include_ele){
      outfile<<"# cv_e => electronic cv \n";
      outfile<<"# alpha_e => electronic volume expansion cofficient x 10^6 \n";
    }
    outfile<<"[AFLOW] "<<STAR130<<"\n";
    outfile<<"[AFLOW_QHA3P_THERMO]START"<<"\n";

    outfile<<"#"<<setw(15)<<"T(K)";
    outfile<<setw(15)<<"U(eV/cell)";
    outfile<<setw(15)<<"F(eV/cell)";
    outfile<<setw(15)<<"S(kB/cell)";
    outfile<<setw(15)<<"Fqh(eV/cell)";
    outfile<<setw(15)<<"B(GPa)";
    outfile<<setw(15)<<"V(A^3)";
    outfile<<setw(12)<<"Bp";
    outfile<<setw(15)<<"gamma";
    outfile<<setw(20)<<"alpha x10^6(1/K)";
    outfile<<setw(15)<<"Cv(kB/cell)";
    outfile<<setw(15)<<"Cp(kB/cell)";
    if(_include_ele){
      outfile<<aurostd::PaddedPRE("Cv_e(kB/cell)",15," ");
      outfile<<aurostd::PaddedPRE("alpha_e (1/K)",15," ");
    }
    outfile<<'\n';
    //indev for relaxed volume
    int eq_index= _qh_energies.get_eqm_ele_dir_index();
    double ele_cv=0.0;
    double ele_CTE=0.0;

    //temperature loop
    for(double t1=Tmin; t1<=Tmax; t1+=delta_T)
      {
	if(!_include_ele)
	  {
	    osfvt<<"[AFLOW] "<<STAR80<<"\n";
	    osfvt << "[QHA3P_EOS_ENERGIES  T="<<setprecision(2)<<std::fixed<<t1<<" K ]START" <<"\n";
	    osfvt<<"#"<<setw(15)<<aurostd::PaddedPRE("V(A^3)",12," ")
		 <<aurostd::PaddedPRE("Ftot(eV/Cell)",15," ")
		 <<aurostd::PaddedPRE("E0K(eV/Cell)",15," ")
		 <<aurostd::PaddedPRE("Fvib(eV/Cell)",15," ");
	    if(includePV)osfvt<<aurostd::PaddedPRE("pV(eV/Cell)",15," ");
	    osfvt<<"\n";
	  }

        //Internal energy in meV
	double THERMO_Ut=e.getInternalEnergy(t1,apl::meV);
        //Total free energy in meV
	double THERMO_Ft=e.getVibrationalFreeEnergy(t1,apl::meV);
        //Total enetropy in kB/cell
	double THERMO_St=e.getVibrationalEntropy(t1, apl::kB);

	double cv=0.0;
	for(uint v=0; v!=_ele_vols.size(); v++){
	  //vibrational free energy
	  double fv=0.0;
	  if(eq_index==(int)v)cv=0.0;

	  for(uint i=0; i!=_freq0.size(); i++){
	    for(int j=1; j<=_freq0[i].rows; j++){
	      if(_freq0[i][j]>0.0){
		fv+=free_energy(_ep_freqs[v][i][j], t1);
		if(eq_index==(int)v)cv+=heat_capacity(_ep_freqs[v][i][j], t1);
	      }
	    }
	  }
	  //meV to eV
	  fv=(fv/(double)_freq0.size())*0.001;

	  if(eq_index==(int)v)cv=(cv/(double)_freq0.size());

	  //electronic cv and thermal expansion
	  if(_include_ele){
	    ele_cv=Electronic_Cv(t1,eq_index);
	    ele_CTE=(2.0/(3.0*_Beq*_Veq))*ele_cv;
            if(t1<1.0){
	    ele_cv=0.0;
	    ele_CTE=0.0;
            }
	  }


	  E[v+1]=_eo[v]+fv;

          //print 
	  if(!_include_ele){
	    osfvt <<setw(15)<<std::setprecision(8)  <<V[v+1]
		  <<setw(15)<<std::setprecision(8)  <<E[v+1]
		  <<setw(15)<<std::setprecision(8)  <<_eo[v]
		  <<setw(15)<<std::setprecision(8)  <<fv;
	    if(includePV)osfvt<<setw(15)<<std::setprecision(8)  <<_pV[v];
	    osfvt<<"\n";
	  }
	}
	if(!_include_ele){
	  osfvt << "[QHA3P_EOS_ENERGIES  T="<<setprecision(2)<<std::fixed<<t1<<" K ]END" <<"\n";
	  osfvt<<"[AFLOW] "<<STAR80<<"\n";
	}
        //energy volume fiiting for each temperature
	md_lsquares_call(V, E);
        //average Gruneisen parameter
	double avg_gp=_scqha.average_gruneisen_parameter(t1);
        //Coefficient of thermal expansion
	double CTE=(avg_gp*cv*0.001)/(_Beq*_Veq);
        //Specific heat at constant pressure
	double cp=cv*0.001 + CTE* CTE * _Beq* _Veq* t1;

	if(t1<0.1)
	  {
	    cv=0.0;
	    cp=0.0;
	    CTE=0.0;
	  }

	outfile<<setw(15)<<t1;
	outfile<<setw(15)<<THERMO_Ut/1000.00;
	outfile<<setw(15)<<THERMO_Ft/1000.00;
	outfile<<setw(15)<<THERMO_St;  
	outfile<<setw(15)<<_Eeq;
	outfile<<setw(15)<<_Beq*160.2176487;
	outfile<<setw(15)<<_Veq;
	outfile<<setw(15)<< _Bp;
	outfile<<setw(15)<< avg_gp;
	outfile<<setw(15)<< CTE*1e6;
	outfile<<setw(15)<< cv/(8.6173324*1e-2);
	outfile<<setw(15)<< cp/(8.6173324*1e-5);
	if(_include_ele){
	  outfile <<setw(15)<< ele_cv/(8.6173324*1e-5); //converted to kB/cell
	  outfile <<setw(15)<< ele_CTE*1.0E6;
	}
	outfile<<'\n';
       vector<double> tmp(2, 0);
       tmp[0]=t1;
       tmp[1]=_Eeq;
       _TF.push_back(tmp);
      }

    string FVTfile =  "aflow.qha3P.FVT.out";
    if(!_include_ele){
      if(!aurostd::stringstream2file(osfvt, FVTfile, "WRITE")) {
	throw APLRuntimeError("Cannot write aflow.qha3P.FVT.out");
      }
      aurostd::StringstreamClean(osfvt);
    }


    outfile<<"[AFLOW_QHA3P_THERMO]END"<<"\n";
    outfile<<"[AFLOW] "<<STAR130<<"\n";
    FVTfile="aflow.qha3P.thermo.out";
    if(!aurostd::stringstream2file(outfile, FVTfile, "WRITE")) {
      throw APLRuntimeError("Cannot write aflow.qha3P.thermo.out");
    }
    aurostd::StringstreamClean(outfile);
 
    //total enthalpy calculation
    if(!_include_ele){
      //without including temperature dependent electronic effect
      total_enthalpy();
    }else{
      //with temperature dependent electronic effect
       enthalpy_incuding_ele(Tmin, Tmax, delta_T);
    }
  }
  // ***************************************************************************************
  //internal energy in meV
  double QHA3POINTS::internal_energy(const double omeg, const double temp)
  { 

    if(omeg< 0.001){
      _logger<<apl::error <<"Frequency too small (<0.001 THz) for U_vib(T)"<<apl::endl;
      exit(0);
    }

    double betaa=47.9924*omeg/temp;
    double U=4.14129*omeg*(1.0/(exp(betaa)-1.0)+0.5);
    return U;
  }
  // ***************************************************************************************
  //Vibrational free energy in meV
  double QHA3POINTS::free_energy(const double omeg, const double temp)
  {
    if(omeg<0.001)
      {
	_logger<<apl::error<<"Frequency too small (<0.001 THz) for F_vib(T)"<<apl::endl;
	exit(0);
      }

    double betaa=47.9924*omeg/temp;
    double F=2.07065*omeg+0.0861733*temp*log(1.0-exp(-1.0*betaa));
    return F;
  }
  // ***************************************************************************************
  //Vibrational entropy
  double QHA3POINTS::entropy(const double omeg, const double temp)
  {
    if(omeg<0.001){
      _logger<<apl::error<< "Frequency too small (<0.001 THz) for S_vib(T)" <<apl::endl;
      exit(0);
    }
    double  betaa=47.9924*omeg/temp;
    double  S=4.14129*omeg/temp/(exp(betaa)-1.0)-0.0861733* log(1.0-exp(-1.0*betaa));
    return S;
  }
  // ***************************************************************************************
  //heat capacity at constant volume
  double QHA3POINTS::heat_capacity(const double omeg, const double temp)
  {
    if(omeg<0.001)
      {
	_logger<<apl::error<<"Frequency too small (<0.001 THz) for Cv(T)" <<apl::endl;
	exit(0);
      }

    double betaa=47.9924*omeg/temp;
    double dumm=(exp(betaa/2.0)-exp(-1.0*betaa/2.0));
    double C=0.0861733*betaa*betaa*1.0/(dumm*dumm);
    return C;
  }
  // ***************************************************************************************
  //get electronic energies for proper file index
  double QHA3POINTS::ElectronicEnergy(double temperature_in_kelvins, uint dir_index)
  {
    if(temperature_in_kelvins<1.0) return 0.0;
    vector<vector<double> > edos=_edos[dir_index];
    double fermi = _fermi_energies[dir_index];
    double sum1=0.0;
    double sum2=0.0;
    double sum3=0.0;

    for(uint j=0; j<edos.size()-1; j++)
      {
        double dE=edos[j+1][0]-edos[j][0];
        double Ediff=edos[j][0]-fermi;
        double f= fermi_dirac_distribution(Ediff, temperature_in_kelvins);
        //continue if edos is .gt. 1.0e-5
        if(edos[j][1]<1.0e-5)continue;
        if(edos[j][0]<=fermi)sum2+=edos[j][1]*edos[j][0]*dE;
        //continue if f is .gt. 1.0e-6
        if(f<1.0e-6)continue;
        sum1+=edos[j][1]*f*edos[j][0]*dE;

        if(!_isequal(f,1.0))
          {
            double x=f*log(f)+(1.0-f)*log(1.0-f);
            sum3+=edos[j][1]*x*dE;
          }
      }
    sum3*=8.6173303*1.0e-5*temperature_in_kelvins;
    return (sum1-sum2+sum3);
  }
  // ***************************************************************************************
  //Total enthalpy calculation
  void QHA3POINTS::total_enthalpy()
  {
    _logger<<"Writing aflow.qha3P.enthalpy.out file"<<apl::endl;

    std::string STAR50 = std::string(50, '*');
    stringstream outfile;
    outfile << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
    outfile << setprecision(4);

    outfile<<"[AFLOW] "<< STAR50 <<"\n";
    outfile<<"# T     => temperature [K]\n";
    outfile<<"# H     => total enthalpy [eV/cell]\n";
    outfile<<"[AFLOW] "<< STAR50 <<"\n";
    outfile<<"[AFLOW_QHA3P_ENTHALPY]START"<<"\n";
    outfile<<"#"<<setw(15)<<"T(K)"<<setw(15)<<"H(eV/cell)"<<"\n";
    /////////////////////////////////////////////////////////////////
    vector< vector <double> > data = _TF; _TF.clear();
    double e_im=0;
    for(uint i=0; i<data.size()-1; i++){
      double t_i=data[i][0];
      double e_i=data[i][1];
      double t_ip=data[i+1][0];
      double e_ip=data[i+1][1];
      //entropy
      double S=(e_ip-e_im)/(2.0*(t_ip-t_i));
      //enthalpy
      double H= e_i - t_i * S;
      e_im=e_i;
      outfile<<setw(15)<<t_i<<setw(15)<<H<<'\n';
    }

    outfile<<"[AFLOW_QHA3P_ENTHALPY]END"<<"\n";
    outfile<<"[AFLOW] "<<STAR50<<"\n";
    string file="aflow.qha3P.enthalpy.out";
    if(!aurostd::stringstream2file(outfile, file, "WRITE")) {
      throw APLRuntimeError("Cannot write aflow.qha3P.enthalpy.out");
    }
    aurostd::StringstreamClean(outfile);
    data.clear();
  }
  // ***************************************************************************************
  //total enthalpy including electronic effect
  void QHA3POINTS::enthalpy_incuding_ele(double USER_TP_TSTART, double USER_TP_TEND, double USER_TP_TSTEP)
  {
    stringstream osfvt;
    std::string STAR80 = std::string(80, '*');

    xvector<double> Ftotal(_eo.size(), 1);
    xvector<double> volume(_eo.size(), 1);
    xvector<double> E0(_eo.size(), 1);

    for(uint i=0; i!=_eo.size(); i++)
      {
	E0[i+1]=ElectronicEnergy(2.0,i);
      }

    vector<vector<double> > data;
    double e_i=0.0;
    for(double K=(double)USER_TP_TSTART; K<=(double)USER_TP_TEND; K+=USER_TP_TSTEP)
      {
	osfvt<<"[AFLOW] "<<STAR80<<"\n";
	osfvt << "[AFLOW_QHA_ENERGIES  T="<<setprecision(2)<<std::fixed<<K<<" K ]START" <<"\n";
	osfvt<<"#"<<setw(15)<<aurostd::PaddedPRE("V(A^3)",12," ")
	     <<aurostd::PaddedPRE("Ftot(eV/Cell)",15," ")
	     <<aurostd::PaddedPRE("E0K(eV/Cell)",15," ")
	     <<aurostd::PaddedPRE("Fvib(eV/Cell)",15," ")
	     <<aurostd::PaddedPRE("Fele(eV/Cell)",15," ");
	osfvt<<"\n";


	for(uint i=0; i!=_eo.size(); i++)
	  {

	    double vib_i=0.0;
	    //vibrational energy calculation
	    for(uint j=0; j!=_freq0.size(); j++){
	      for(int k=1; k<=_freq0[i].rows; k++){
		if(_freq0[j][k]>0.0){
		  vib_i+=free_energy(_ep_freqs[i][j][k], K);
		}
	      }
	    }

	    //meV to eV
	    vib_i=(vib_i/(double)_freq0.size())*0.001;

	    //electronic energy calculation
	    e_i=ElectronicEnergy(K,i)-E0[i+1];

	    if(K<1.0)e_i=0.0;
	    if(abs(e_i)<0.001)e_i=0.0;

	    Ftotal[i+1]=_eo[i]+vib_i+e_i;
	    volume[i+1]=_ele_vols[i];

	    osfvt <<setw(15)<<std::setprecision(8)  <<volume[i+1]
		  <<setw(15)<<std::setprecision(8)  <<Ftotal[i+1]
		  <<setw(15)<<std::setprecision(8)  <<_eo[i]
		  <<setw(15)<<std::setprecision(8)  <<vib_i
		  <<setw(15)<<std::setprecision(8)  <<e_i;
	    osfvt<<"\n";
	  }

	//fitting
	md_lsquares_call(volume, Ftotal);
	vector<double> tmp(2,0);
	tmp[0]=K;
	tmp[1]=_Eeq;
	data.push_back(tmp);
	osfvt << "[AFLOW_QHA_ENERGIES  T="<<setprecision(2)<<std::fixed<<K<<" K ]END" <<"\n";
	osfvt<<"[AFLOW] "<<STAR80<<"\n";
      }

    string FVTfile =  "aflow.qha3P.FVT.out";
    if(!aurostd::stringstream2file(osfvt, FVTfile, "WRITE")) {
      throw APLRuntimeError("Cannot write aflow.qha3P.FVT.out");
    }
    aurostd::StringstreamClean(osfvt);


    _logger<<"Writing aflow.qha3P.enthalpy.out file"<<apl::endl;

    std::string STAR50 = std::string(50, '*');
    stringstream outfile;
    outfile << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
    outfile << setprecision(4);

    outfile<<"[AFLOW] "<< STAR50 <<"\n";
    outfile<<"# T     => temperature [K]\n";
    outfile<<"# H     => total enthalpy [eV/cell]\n";
    outfile<<"[AFLOW] "<< STAR50 <<"\n";
    outfile<<"[AFLOW_QHA3P_ENTHALPY]START"<<"\n";
    outfile<<"#"<<setw(15)<<"T(K)"<<setw(15)<<"H(eV/cell)"<<"\n";

    double e_im=data[0][1];
    for(uint i=0; i<data.size()-1; i++){
      double t_i=data[i][0];
      double e_i=data[i][1];
      double t_ip=data[i+1][0];
      double e_ip=data[i+1][1];
      double S=(e_ip-e_im)/(2.0*(t_ip-t_i));
      double H= e_i - t_i * S;
      e_im=e_i;
      outfile<<setw(15)<<t_i<<setw(15)<<H<<'\n';
    }

    outfile<<"[AFLOW_QHA3P_ENTHALPY]END"<<"\n";
    outfile<<"[AFLOW] "<<STAR50<<"\n";
    string file="aflow.qha3P.enthalpy.out";
    if(!aurostd::stringstream2file(outfile, file, "WRITE")) {
      throw APLRuntimeError("Cannot write aflow.qha3P.enthalpy.out");
    }
    aurostd::StringstreamClean(outfile);
    data.clear();
  }
  // ***************************************************************************************
  //get electronic specific heat for a proper file index
  double QHA3POINTS::Electronic_Cv(double temperature_in_kelvins, uint dir_index)
  {
    if(temperature_in_kelvins<1.0) return 0.0;
    vector<vector<double> > edos=_edos[dir_index];
    double fermi = _fermi_energies[dir_index];
    double sum_t=0.0;
    double sum_tp=0.0;

    double dE=edos[1][0]-edos[0][0];
    for(uint j=0; j<edos.size(); j++)
      {
        double Ediff=edos[j][0]-fermi;
        double f= fermi_dirac_distribution(Ediff, temperature_in_kelvins);
        //continue if edos is .gt. 1.0e-5
        if(edos[j][1]<1.0e-5)continue;
        //continue if f is .gt. 1.0e-6
        if(f<1.0e-6)continue;

        if(!_isequal(f,1.0))
          {
            double x=f*log(f)+(1.0-f)*log(1.0-f);
            sum_t+=edos[j][1]*x;
          }
      }
    sum_t*=-8.6173303*1.0e-5*dE;
    for(uint j=0; j<edos.size(); j++)
      {
        double Ediff=edos[j][0]-fermi;
        double f= fermi_dirac_distribution(Ediff, temperature_in_kelvins+10.0);
        //continue if edos is .gt. 1.0e-5
        if(edos[j][1]<1.0e-5)continue;
        //continue if f is .gt. 1.0e-6
        if(f<1.0e-6)continue;

        if(!_isequal(f,1.0))
          {
            double x=f*log(f)+(1.0-f)*log(1.0-f);
            sum_tp+=edos[j][1]*x;
          }
      }
    sum_tp*=-8.6173303*1.0e-5*dE;

    double cv=temperature_in_kelvins*(sum_tp-sum_t)/10.0;
    return (cv);
  }
  // ***************************************************************************************
  //Fermi-Dirac distribution function
  double QHA3POINTS::fermi_dirac_distribution(const double delE, const double t)
  {
    double beta=1.0/(KBOLTZEV*t); // (_kB*t); //CO181019
    double x=(delE)*beta;
    x=exp(x)+1.0;
    return 1.0/x;
  }
  // ***************************************************************************************
  void QHA3POINTS::set_include_ele(bool b)
  {
    _include_ele=b;
  }
  // ***************************************************************************************
  //check file existence
  bool QHA3POINTS::exists_test0 (const std::string& name)
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
}//apl namespace end
