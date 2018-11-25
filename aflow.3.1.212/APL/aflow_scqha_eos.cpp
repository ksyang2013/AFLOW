// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2015           *
// *                Aflow PINKU NATH - Duke University 2014-2018             *
// *                                                                         *
// ***************************************************************************
// Written by Pinku Nath
// pn49@duke.edu

#include "aflow_apl.h"

#define _isnegative(a) (a<MIN_EIGEN_TRESHOLD) ? true : false

#if GCC_VERSION >= 40400   // added two zeros
#define AFLOW_APL_MULTITHREADS_ENABLE 1
#include <thread>
#else
#warning "The multithread parts of APL will be not included, since they need gcc 4.4 and higher (C++0x support)."
#endif

/*
  This is the SCQHA implementation [ Ref. Comp. Mat. Sci. 120 (2016) 84–93].
  In this technique the thermodynamic properties are calculated by self-consistently
  minimization of volume.
*/

namespace apl
{
  // ***************************************************************************************
  //constructor
  SCQHAEOS::SCQHAEOS(SCQHA_QHA3P& scqha, QH_ENERGIES& qh_energies, Logger& l):_scqha(scqha), _qh_energies(qh_energies), _logger(l)
  {
    clear();
    _logger<<"SCQHA calculation is started, "<<apl::endl;
  }
  // ***************************************************************************************
  //destructor
  SCQHAEOS::~SCQHAEOS()
  {
    this->clear();
  }
  // ***************************************************************************************
  //clear variables
  void SCQHAEOS::clear()
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
    _TF.clear();
    _Eeq=0.0, _Beq=0.0, _Veq=0.0, _Bp=0.0;
    _dE_dV=0.0, _d2E_dV2=0.0;
    _TV.clear();
    _inpiut_T.clear();
  }
  // ***************************************************************************************
  //initialize all variables with this function
  bool SCQHAEOS::import_variables()
  {
    _logger <<"SCQHA Importing variables, "<<apl::endl;
    _eo=_qh_energies.get_eo();
    _ele_vols=_qh_energies.get_ele_vols();
    _freq0=_scqha.get_freqs_mesh();
    _freqM=_scqha.get_freqs_meshM();
    _freqP=_scqha.get_freqs_meshP();
    _scqha_volumes=_scqha.get_qha_gpvol();
    _weights=_scqha.get_weights();
    return check_size();
  }
  // ***************************************************
  //check sizes of vectors before using them
  bool SCQHAEOS::check_size()
  {
    if(_scqha_volumes.size()!=3)
      {
	_logger << apl::error << "_scqha_volumes.size()!=3 "<< apl::endl;
	return false;
      }
    if(_eo.size()==0)
      {
	_logger << apl::error << "_eo.size()==0 "<< apl::endl;
	return false;
      }
    if(_ele_vols.size()==0)
      {
	_logger << apl::error << "_ele_vols.size()==0 "<< apl::endl;
	return false;
      }
    if(_freq0.size()==0)
      {
	_logger << apl::error << "_freq0.size()==0 "<< apl::endl;
	return false;
      }
    if(_freqM.size()==0)
      {
	_logger << apl::error << "_freqM.size()==0 "<< apl::endl;
	return false;
      }
    if(_freqP.size()==0)
      {
	_logger << apl::error << "_freqP.size()==0 "<< apl::endl;
	return false;
      }
    if(_weights.size()==0)
      {
	_logger << apl::error << "_weights.size()==0 "<< apl::endl;
	return false;
      }
    _nBranches=_freq0[0].rows;
    return true;
  }
  // ***************************************************************************************
  //energy volume fitting function
  void SCQHAEOS::fitting()
  {
    _logger<<"SCQHA Fitting E-V data"<<apl::endl;

    xvector<double> E(_eo.size(), 1);
    xvector<double> V(_eo.size(), 1);

    for(uint i=0; i!=_eo.size(); i++)
      {
	E[i+1]=_eo[i];
	V[i+1]=_ele_vols[i];
      }
    md_lsquares_call(V,E);
  }
  // ***************************************************************************************
  //energy volume nonlinear fitting function
  void SCQHAEOS::md_lsquares_call(const xvector<double> &V, const xvector<double> &E)
  {
    md_lsquares mdfit;
    mdfit.clear();
    for(int i=1; i<=V.rows; i++)
      {
        mdfit.Xdata.push_back(V[i]);
        mdfit.Ydata.push_back(E[i]);
      }
    // it does both linear and nonlinear fit
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
  //doing more refinments after fitting 
  bool SCQHAEOS::more_refinement(const xvector<double> &E, const xvector<double> &V,
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
  //Calculate Taylor cofficients 
  void SCQHAEOS::calculate_freq_derivative()
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
      _logger.initProgressBar("Calculating taylor coefficients for SCQHA");
    else
      _logger.initProgressBar("Calculating taylor coefficients for SCQHA  (" + stringify(ncpus) + " threads)");

    // Distribute the calculation
    int startIndex, endIndex;
    std::vector< std::thread* > threads;
    vector<vector<int> > thread_dist = getThreadDistribution((int) _freq0.size(), ncpus);
    for (int icpu = 0; icpu < ncpus; icpu++) {
      startIndex = thread_dist[icpu][0];
      endIndex = thread_dist[icpu][1];
      threads.push_back( new std::thread(&SCQHAEOS::calculate_derivative,this,startIndex,endIndex) );
    }

/* OBSOLETE ME180801
    for(int icpu = 0; icpu < ncpus; icpu++) {
      startIndex = icpu * qpointsPerCPU;
      endIndex = startIndex + qpointsPerCPU;
      if( ( (uint)endIndex > _freq0.size() ) ||
          ( ( icpu == ncpus-1 ) && ( (uint)endIndex < _freq0.size() ) ) )
        endIndex = _freq0.size();
      threads.push_back( new std::thread(&SCQHAEOS::calculate_derivative,this,startIndex,endIndex) );
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
    _logger<<"Calculating taylor coefficients for SCQHA without threads "<<apl::endl;
    calculate_derivative(0, _freq0.size());
#endif
  }
  // ***************************************************************************************
  //Calculate Taylor cofficients in threads
  void SCQHAEOS::calculate_derivative(int startIndex, int endIndex)
  {
    double dv=_scqha_volumes[0]-_scqha_volumes[2];
    double dv2=0.25*dv*dv; 

    for(int i=startIndex; i<endIndex; i++){
      for(uint j=1; j<=_nBranches; j++)
	{
	  _d1fdv1[i][j]= (_freqM[i][j] - _freqP[i][j])/dv;     
	  _d2fdv2[i][j]= (_freqM[i][j] + _freqP[i][j] - 2.0* _freq0[i][j])/dv2;  
	}
    }
  }
  // ***************************************************************************************
  //SCQHA temperature loop
  void SCQHAEOS::sccycle(double Tmin, double Tmax, double delta_T)
  {
    if(_iszero(Tmin))Tmin=1.0;

    _logger<<"Entered into SCQHA cycle "<<apl::endl;

    //E-V fitting
    fitting();

    double E_min=E_V(_Veq);
    std::string STAR150 = std::string(150, '*');
    std::string STAR10 = std::string(10, '*');
    stringstream scf_thermo, scf_thermo_p, scf_err;

    scf_err << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
    scf_err << setprecision(8);

    scf_err<<std::setprecision(4)<<std::fixed;

    scf_thermo << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
    scf_thermo << setprecision(8);

    scf_thermo_p << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
    scf_thermo_p << setprecision(8);

    _logger <<"Writng aflow.scqha.iter.out file"<<apl::endl; 
    _logger <<"Writng aflow.scqha.thermo.out file"<<apl::endl; 
    _logger <<"Writng aflow.scqha.pressure.out file"<<apl::endl;
 

    scf_thermo_p<<std::setprecision(6)<<std::fixed;
    scf_err <<"[AFLOW] "<<STAR10<<"\n";
    scf_err <<"[AFLOW_SCQHA_VOLUME_ITER]START"<<"\n";
    scf_err <<"#V_T = volume in A^3 uint "<<"\n";
    scf_thermo <<"[AFLOW_SCQHA_THERMO]START" <<"\n";
    scf_thermo<<std::setprecision(3)<<std::fixed;
    scf_thermo <<"[AFLOW] "<<STAR150<<"\n";
    scf_thermo <<"#defination of variables \n";
    scf_thermo<<"#P_ext(meV/A^3)= "<<_pext<<"\n";
    scf_thermo<<"#E_min(eV) = "<<E_min<<"\n";
    scf_thermo<<"#Ee - E_min = Electronic_energy - Emin(meV)"<<"\n";
    scf_thermo<<"#Fvib = Free energy of vibration\n";
    scf_thermo<<"#U = Vibrational internal energy\n";
    scf_thermo<<"#G = U+P_ext*V\n";
    scf_thermo<<"#Svib = entropy of vib\n";
    scf_thermo<<"#Cv = specific heat at const vol\n";
    scf_thermo<<"#Cp = specific heat at const Presure\n";
    scf_thermo<<"#Ce = electronic specific heat\n";
    scf_thermo<<"#alpha_V = volume thermal expansion cofficient *10^(-6)\n";
    scf_thermo<<"#V_T = equilibrium volume at temperature T \n";
    scf_thermo<<"# gamma => Gruneisen Parameter \n";
    scf_thermo <<"[AFLOW] "<<STAR150<<"\n";
    scf_thermo <<"[AFLOW_SCQHA_THERMO]START" <<"\n";
    scf_thermo<<"#"<<setw(15)<<"T[K]"
	      <<setw(15)<<"Ee(meV/Cell)"
	      <<setw(15)<<"Fvib(meV/Cell)"
	      <<setw(15)<<"U(meV/Cell)"
	      <<setw(15)<<"G(meV/Cell)"
	      <<setw(15)<<"Svib(kB/Cell)"
	      <<setw(15)<<"Cv(kB/Cell)"
	      <<setw(15)<<"Cp(kB/Cell)"
	      <<setw(15)<<"Ce(kB/Cell)"
	      <<setw(15)<<"alpha_V(1/K)"
	      <<setw(15)<<"V_T(A^3/Cell)"
	      <<setw(15)<<"gamma";
    scf_thermo<<"\n";
    scf_thermo<<std::setprecision(6)<<std::fixed;
    scf_thermo_p <<"[AFLOW] "<<STAR150<<"\n";
    scf_thermo_p <<"#defination of variables\n";
    scf_thermo_p<<"#B_T = isothermal bulk modulus\n";
    scf_thermo_p<<"#B_e = electronic bulk modulus\n";
    scf_thermo_p<<"#B_gamma = first order bulk modulus due to change in Gruneisen parameter\n";
    scf_thermo_p<<"#B_Dgamma = second order bulk modulus due to change in Gruneisen parameter\n";
    scf_thermo_p<<"#P_e = electronic pressure\n";
    scf_thermo_p<<"#P_gamma = anharmonic phonon pressure\n";
    scf_thermo_p <<"[AFLOW] "<<STAR150<<"\n";
    scf_thermo_p <<"[AFLOW_SCQHA_PRESSURE]START" <<"\n";
    scf_thermo_p<<"#"<<setw(15)<<"T[K]"
		<<setw(15)<<"B_T(GPa)"
		<<setw(15)<<"B_e(GPa)"
		<<setw(15)<<"B_gamma(GPa)"
		<<setw(15)<<"B_Dgamma(GPa)"
		<<setw(15)<<"P_e(GPa)"
		<<setw(15)<<"P_gamma(GPa)";
    scf_thermo_p<<"\n";

    //calculating Taylor coefficients
    calculate_freq_derivative();
    //Print Taylor coefficients
    print_freq_taylor_cofficients();

    //number of temperature points
    int N_T=int((Tmax-Tmin)/delta_T)+1;

    //V(T)
    vector<double> V_T; V_T.clear();
    //alpha(T)
    vector<double> alpha_T; alpha_T.clear();
    //B(T)
    vector<double> B_T; B_T.clear();

    V_T.resize(2,0.0);
    alpha_T.resize(2,0.0);
    B_T.resize(2,0.0);

    //initialization of temperature dependent volume
    V_T[0] = (1.0+0.1)*_Veq;


    
    uint ksize=_freq0.size();
    vector<xvector<double> > dvx;dvx.clear();
    //Gruneisen Parameter
    vector<vector<xvector<double> > > gruneisen;
    //grequency(T)
    vector<vector<xvector<double> > > omega_T;
    gruneisen.clear(); 
    omega_T.clear(); 
    xvector<double> dx(_nBranches,1);
    dvx.resize(ksize,dx);
    gruneisen.resize(2,dvx);
    omega_T.resize(2,dvx);
    //specific heat 
    vector<vector<xvector<double> > >Cv; Cv.clear();
    Cv.resize(2,dvx);

    double temperature = Tmin;
    double ph_pressure=0.0;
    double P_gamma=0.0;
    double P_e=0.0;
    double B_e=0.0;
    double B_gamma=0.0;
    double B_2=0.0;
    double B_2_dw2=0.0;

    //self consistent cycle to minimize the volume
    for(uint i=0; i!=1000; i++)
      {
	vector<xvector<double> > U_ph; U_ph.clear();
	U_ph.resize(ksize,dx);
	for(uint j=0; j!=ksize; j++)
	  {
	    for(uint k=1; k<=_nBranches; k++)
	      {
		if(_freq0[j][k]<0.01) continue;
		if( (abs(V_T[0]-ph_pressure)/V_T[0]) > 1.0e-3)
		  {
                    //initializing frequencies
		    omega_T[0][j][k]= _freq0[j][k];
		  }else{
                  //updating frequencies
		  omega_T[0][j][k]= _freq0[j][k]+ _d1fdv1[j][k]*(V_T[0]-_scqha_volumes[1]) + 0.5*_d2fdv2[j][k]*pow((V_T[0]-_scqha_volumes[1]), 2.0);
		}
		if(omega_T[0][j][k]<0.01) continue;
                //updating Gruneisen parameter
		gruneisen[0][j][k]= (-1.0*V_T[0])/omega_T[0][j][k]* (_d1fdv1[j][k] + _d2fdv2[j][k]*(V_T[0]-_scqha_volumes[1]));
                //updating internal energies
		U_ph[j][k]=internal_energy(omega_T[0][j][k],temperature);
                //total phonon pressure updating
		ph_pressure+=U_ph[j][k]*gruneisen[0][j][k]*_weights[j];
	      }
	  }
        ph_pressure/=(double)(ksize);
        //calculating electronic pressure
        derivatives(V_T[0]);
        //total pressure updating
        ph_pressure=ph_pressure/(_dE_dV+_pext);

        //checking self-consistent volume is converged or not
	if(abs(V_T[0]-ph_pressure)/V_T[0] > 1e-5)
	  {
	    V_T[0]=V_T[0]+(ph_pressure-V_T[0])*0.001;
	    ph_pressure=0.0;
	    scf_err<<"V_T = "<<V_T[0]<<"\n";
	  }else{
          //Thermodynamic properties calculations
	  B_e=V_T[0]*_d2E_dV2;
	  P_e=-1.0*_dE_dV;
	  B_gamma=0.0;
	  B_2=0.0;
	  P_gamma=0.0;

          for(uint j=0; j!=ksize; j++){
	    for(uint k=1; k<=_nBranches; k++){
	      if(_freq0[j][k]<0.1) continue;
	      Cv[0][j][k]=heat_capacity(omega_T[0][j][k],temperature);
	      B_gamma+=(U_ph[j][k]-temperature*Cv[0][j][k])*gruneisen[0][j][k]*gruneisen[0][j][k]*_weights[j];
	      B_2+=U_ph[j][k]*((1.0+gruneisen[0][j][k])*gruneisen[0][j][k]-V_T[0]*V_T[0]/omega_T[0][j][k]*_d2fdv2[j][k])*_weights[j];
	      B_2_dw2+=U_ph[j][k]* V_T[0]*V_T[0]/omega_T[0][j][k]*_d2fdv2[j][k]*_weights[j];
	      P_gamma+=U_ph[j][k]*gruneisen[0][j][k]*_weights[j];
            }
          }
	  B_gamma=B_gamma/(double)(ksize)/V_T[0];
	  B_2=-1.0*B_2/(double)(ksize)/V_T[0];
	  B_2_dw2=B_2_dw2/(double)(ksize)/V_T[0];
	  P_gamma=P_gamma/(double)(ksize)/V_T[0];
           
	  B_T[0]=B_e+B_gamma+B_2+P_gamma;
	  scf_err<<"#Total iteration =  "<<i<<" Steps."<<"\n";
	  break;
	}
      }//iteration loop
 
    //calculating thermodynamic properties

    //temperature independent electronic enegy
    double Ee_T=E_V(V_T[0]);
    //temperature dependent electronic enegy
    double    Ee_tt=0.0;
    //total internal energy
    double    U_tot=0.0;
    //total entropy
    double    S_tot=0.0;
    //total specific heat
    double    Cv_tot=0.0;
    //total specific heat at constant pressure
    double    Cp_tot=0.0;
    //total Helmotz energy
    double    H_tot=0.0;
    //phonon pressure
    P_gamma=0.0;
    //electronic pressure
    P_e=0.0;
    //electronic Bulk modulus
    B_e=0.0;
    //phonon Bulk modulus
    B_gamma=0.0;
    B_2=0.0;
    B_2_dw2=0.0;
    _logger<<"Entering into Temperature loop "<<apl::endl;
    for(int uT=1; uT!=N_T; uT++)
      {
        double avg_gp=0.0;
        uint i=1;
	double    F_tot=0.0;
	double temperature=Tmin+(double)(uT)*delta_T;

        for(uint j=0; j!=ksize; j++){
          for(uint k=1; k<=_nBranches; k++){
	    alpha_T[i]+=Cv[i-1][j][k]*gruneisen[i-1][j][k]*_weights[j];
	  }
	}
        //thermal expansion
        alpha_T[i]=alpha_T[i]/(double)(ksize)/(B_T[i-1]*V_T[i-1]);

        V_T[i]=(1.0+alpha_T[i]*delta_T)*V_T[i-1];
   
        derivatives(V_T[i]);

        B_e=V_T[i]*_d2E_dV2;
        P_e=-1.0*_dE_dV;

        for(uint j=0; j!=ksize; j++){
          for(uint k=1; k<=_nBranches; k++){
          
	    if(_freq0[j][k]<0.01)continue;
            
            //frequency updating
	    omega_T[i][j][k]=_freq0[j][k]+_d1fdv1[j][k]*(V_T[i]-_scqha_volumes[1])+0.5*_d2fdv2[j][k]*(V_T[i]-_scqha_volumes[1])*(V_T[i]-_scqha_volumes[1]);

	    if(omega_T[i][j][k]<0.01)continue;

            //Gruneisen parameter updating
	    gruneisen[i][j][k]=-1.0*V_T[i]/omega_T[i][j][k]*(_d1fdv1[j][k]+
							     _d2fdv2[j][k]*(V_T[i]-_scqha_volumes[1]));

	    if(omega_T[i][j][k]<0.01){
	      _logger<<apl::error<<"Frequency too small at "<<temperature<<" K "<<apl::endl;
	      exit(0);
	    }
	    double F_ph_i=free_energy(omega_T[i][j][k],temperature);
	    double U_ph_i=internal_energy(omega_T[i][j][k],temperature);
	    double S_ph_i=entropy(omega_T[i][j][k],temperature);
            Cv[i][j][k]=heat_capacity(omega_T[i][j][k],temperature);

	    B_gamma+=(U_ph_i-temperature*Cv[i][j][k])*gruneisen[i][j][k]*gruneisen[i][j][k]*_weights[j];

	    B_2+=U_ph_i*((1.0+gruneisen[i][j][k])*gruneisen[i][j][k]-V_T[i]*V_T[i]/omega_T[i][j][k]*_d2fdv2[j][k])*_weights[j];

	    B_2_dw2+=U_ph_i*V_T[i]*V_T[i]/omega_T[i][j][k]*_d2fdv2[j][k]*_weights[j];

	    P_gamma+=U_ph_i*gruneisen[i][j][k]*_weights[j];

	    F_tot+=F_ph_i*_weights[j];
	    U_tot+=U_ph_i*_weights[j];
	    S_tot+=S_ph_i*_weights[j];
	    Cv_tot+=Cv[i][j][k]*_weights[j];
            avg_gp+=gruneisen[i][j][k]*Cv[i][j][k]*_weights[j];
	  }
	}

        B_gamma=B_gamma/(double)(ksize)/V_T[i];
        B_2=-1.0*B_2/(double)(ksize)/V_T[i];
        B_2_dw2=B_2_dw2/(double)(ksize)/V_T[i];
        P_gamma=P_gamma/(double)(ksize)/V_T[i];
        avg_gp=avg_gp/Cv_tot;

        B_T[i]=B_e+B_gamma+B_2+P_gamma;

        F_tot=F_tot/double(ksize);
        U_tot=U_tot/double(ksize);
        H_tot=U_tot+_pext*V_T[i];
        S_tot=S_tot/(double)(ksize);
        Cv_tot=Cv_tot/(double)(ksize);
        Cp_tot=Cv_tot+alpha_T[i]*alpha_T[i]*B_T[i]*V_T[i]*temperature;

        Ee_tt=Ee_T;

        Ee_T=E_V(V_T[i]);
        //writing to files
        scf_thermo<<std::setprecision(2)<<setw(15)<<temperature;
        scf_thermo<<std::setprecision(6);
        scf_thermo<<setw(15)<<Ee_T-E_min;
        scf_thermo<<setw(15)<<F_tot;
        scf_thermo<<setw(15)<<U_tot;
        scf_thermo<<setw(15)<<H_tot;
        //converting to kB/cell
        scf_thermo<<setw(15)<<S_tot/(8.6173324*1e-2);
        //converting to kB/cell
        scf_thermo<<setw(15)<<Cv_tot/(8.6173324*1e-2);
        //converting to kB/cell
        scf_thermo<<setw(15)<<Cp_tot/(8.6173324*1e-2);
        //converting to meV/cell
        scf_thermo<<setw(15)<<((Ee_T-Ee_tt)/delta_T)/(8.6173324*1e-2);
        //
        scf_thermo<<setw(15)<<alpha_T[i]*1e6;
        //
        scf_thermo<<setw(15)<<V_T[i];
        //
        scf_thermo<<setw(15)<<avg_gp;
        scf_thermo<<"\n";

        scf_thermo_p<<std::setprecision(2)<<setw(15)<<temperature;
        scf_thermo_p<<std::setprecision(6);
        //cobvering to GPa
        scf_thermo_p<<setw(15)<<0.16022*B_T[i];
        scf_thermo_p<<setw(15)<<0.16022*B_e;
        scf_thermo_p<<setw(15)<<0.16022*B_gamma;
        scf_thermo_p<<setw(15)<<0.16022*B_2;
        scf_thermo_p<<setw(15)<<0.16022*P_e;
        scf_thermo_p<<setw(15)<<0.16022*P_gamma;
        scf_thermo_p<<"\n";

	//save current itaration 
	alpha_T[0]=alpha_T[1];
	Cv[0]=Cv[1];
	gruneisen[0]=gruneisen[1];
	B_T[0]=B_T[1];
	V_T[0]=V_T[1];
        vector<double> tmp(2, 0);
        tmp[0]=temperature;
        tmp[1]=Ee_T+F_tot;
	_TF.push_back(tmp);
	tmp.clear();
	tmp.resize(2, 0);
	//saving data to calculate temperature dependent PDIS
	for (uint w=0; w<_inpiut_T.size(); w++)
	  {
	    if (abs(temperature-_inpiut_T[w])<0.001){
	      tmp[0]=temperature;
	      tmp[1]=V_T[i];
	      _TV.push_back(tmp);
	    }
	  }  


      }//temperature loop ends



    _logger<<"SCQHA cycle completed"<<apl::endl;
    scf_thermo <<"[AFLOW_SCQHA_THERMO]END" <<"\n";
    scf_thermo <<"[AFLOW] "<<STAR150<<"\n";

    scf_err <<"[AFLOW_SCQHA_VOLUME_ITER]END" <<"\n";
    scf_err <<"[AFLOW] "<<STAR10<<"\n";

    scf_thermo_p <<"[AFLOW_SCQHA_PRESSURE]END" <<"\n";
    scf_thermo_p <<"[AFLOW] "<<STAR150<<"\n";

    string filename = "aflow.scqha.iter.out";
    aurostd::stringstream2file(scf_err, filename);
    if (!aurostd::FileExist(filename)) {
      throw apl::APLRuntimeError("Cannot open output aflow.scqha.iter.out file.");
    }


    filename = "aflow.scqha.thermo.out";
    aurostd::stringstream2file(scf_thermo, filename);
    if (!aurostd::FileExist(filename)) {
      throw apl::APLRuntimeError("Cannot open output aflow.scqha.thermo.out file.");
    }

    filename = "aflow.scqha.pressure.out";
    aurostd::stringstream2file(scf_thermo_p, filename);
    if (!aurostd::FileExist(filename)) {
      throw apl::APLRuntimeError("Cannot open output aflow.scqha.pressure.out file.");
    }

    //calculating total enthalpy
    total_enthalpy();
  }
  // ***************************************************************************************
  //Calculate total enthalpy
  void SCQHAEOS::total_enthalpy()
  {
    _logger<<"Writing aflow.scqha.enthalpy.out file"<<apl::endl;

    std::string STAR50 = std::string(50, '*');
    stringstream outfile;
    outfile << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
    outfile << setprecision(4);

    outfile<<"[AFLOW] "<< STAR50 <<"\n";
    outfile<<"# T     => temperature [K]\n";
    outfile<<"# H     => total enthalpy [eV/cell]\n";
    outfile<<"[AFLOW] "<< STAR50 <<"\n";
    outfile<<"[AFLOW_SCQHA_ENTHALPY]START"<<"\n";
    outfile<<"#"<<setw(15)<<"T(K)"<<setw(15)<<"H(eV/cell)"<<"\n";
    /////////////////////////////////////////////////////////////////
    vector< vector <double> > data = _TF; _TF.clear(); 
    double e_im=0;
    for(uint i=0; i<data.size()-1; i++){
      double t_i=data[i][0];
      double e_i=data[i][1];
      double t_ip=data[i+1][0];
      double e_ip=data[i+1][1];
      double S=(e_ip-e_im)/(2.0*(t_ip-t_i));
      double H= e_i - t_i * S;
      e_im=e_i;
      if(i==0)outfile<<setw(15)<<t_i<<setw(15)<<data[i][1]/1000.00<<'\n';
      else outfile<<setw(15)<<t_i<<setw(15)<<H/1000.00<<'\n';

    }

    outfile<<"[AFLOW_SCQHA_ENTHALPY]END"<<"\n";
    outfile<<"[AFLOW] "<<STAR50<<"\n";
    string file="aflow.scqha.enthalpy.out";
    if(!aurostd::stringstream2file(outfile, file, "WRITE")) {
      throw APLRuntimeError("Cannot write aflow.scqha.enthalpy.out");
    }
    aurostd::StringstreamClean(outfile);
    data.clear();
  }
  // ***************************************************************************************
  //mode internal energy in meV
  double SCQHAEOS::internal_energy(const double omeg, const double temp)
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
  //energy derivative
  void SCQHAEOS::derivatives(double vol)
  {
    xvector<double> coeff(3,1);
    coeff(1)= _Beq*1000.0;
    coeff(2)= _Veq;
    coeff(3)= _Bp;
    _dE_dV = (coeff(1) *(1.0 + pow((coeff(2)/vol),coeff(3))/(-1.0 + coeff(3))))/coeff(3) 
      - (coeff(1)*coeff(2) * pow((coeff(2)/vol), (-1.0 + coeff(3))))/((-1.0 + coeff(3)) * vol);

    _d2E_dV2 = (coeff(1)* pow(coeff(2), 2.0 ) * pow((coeff(2)/vol), (-2.0 + coeff(3))))/ pow(vol, 3.0);
  }
  // ***************************************************************************************
  //mode free energy in meV
  double SCQHAEOS::free_energy(const double omeg, const double temp)
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
  //mode entropy in meV/K
  double SCQHAEOS::entropy(const double omeg, const double temp)
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
  //mode heat capacity
  double SCQHAEOS::heat_capacity(const double omeg, const double temp)
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
  //energy as a function of volume
  double SCQHAEOS::E_V(const double volume)
  {
    double V0  = _Veq;
    double E0  = _Eeq*1000.0;
    double B0  = _Beq*1000.0;
    double B0p = _Bp;

    double V = volume;
    double t1 = pow((V0/V), (1.0/3.0));
    double t2 = t1*t1;
    double t3 = t2-1.0;
    double energy = 9.0/8.0*B0*V0*t3*t3*(B0p*t3/2.0 - 2.0*t2 + 3.0) + E0;
    return energy;
  }
  // ***************************************************************************************
  //write Taylor cofficients
  void SCQHAEOS::print_freq_taylor_cofficients()
  {
    _logger<<"Writing aflow.mesh.taylor_cofficients.out file," <<apl::endl;
    std::string STAR50 = std::string(50, '*');

    stringstream out;
    out << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
    out << std::setprecision(8) << std::fixed;
    out<<"[AFLOW]"<<STAR50<<'\n';
    out<<"#f=frequency [THz]"<<"\n";
    out<<"#df/dv=First Taylor coefficient ([THz]/A^3)"<<"\n";
    out<<"#d^2f/dv^2=Second Taylor coefficient (THz/A^6)"<<"\n";
    out<<"#branches=phonon branch"<<"\n";
    out<<"[AFLOW]"<<STAR50<<'\n';
    out <<"[AFLOW_SCQHA_TAYLOR_COEFFICIENTS]START"<<"\n";
    out<<"#"<<setw(15)<<"f"<<setw(15)<<"df/dv"<<setw(15)<<"d^2f/dv^2"<<setw(15)<<"branches"<<"\n";
    for(uint i=0; i!=_d1fdv1.size(); i++){
      for(int j=1; j<=_d1fdv1[i].rows; j++){
	out<<setw(15)<<_freq0[i][j]<<setw(15)<<_d1fdv1[i][j]<<setw(15)<<_d2fdv2[i][j]<<setw(15)<<j<<'\n';
      }out<<"#\n";}
    out <<"[AFLOW_SCQHA_TAYLOR_COEFFICIENTS]END"<<"\n";
    out<<"[AFLOW]"<<STAR50<<'\n';
    string filename = "aflow.mesh.taylor_cofficients.out";
    aurostd::stringstream2file(out, filename);
    if (!aurostd::FileExist(filename)) {
      throw apl::APLRuntimeError("Cannot open aflow.mesh.taylor_cofficients.out file.");
    }
  }
  // ***************************************************************************************
  vector<vector<double> > SCQHAEOS::get_TV_data()
  {
    return _TV;
  }
  // ***************************************************************************************
  void  SCQHAEOS::set_input_temperature(const vector<double> &a)
  {
    _inpiut_T=a;
  }
  // ***************************************************************************************
}//apl namespace end
