// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2015           *
// *                Aflow PINKU NATH - Duke University 2014-2018             *
// *                                                                         *
// ***************************************************************************
// Written by Pinku Nath
// pn49@duke.edu

/*
   This class computes EOS using QHA method [Ref. Comp. Mat. Sci. 125, 82–91 (2016)]
*/

#include "aflow_apl.h"
using namespace std;
#define _kBeV 8.6173303
#define _heV  4.135667662
namespace apl
{
  QHAEOS::QHAEOS(QHA& qha, QH_ENERGIES& qhen, Logger& l): _qha(qha),_qhen(qhen),_logger(l)
  {
    _logger << "QHAEOS calculation is started "<< apl::endl; 
    //clear all memories before using
    clear();
  }
  // ***************************************************
  QHAEOS::~QHAEOS()
  {
    this->clear();
  }
  // ***************************************************
  void QHAEOS::clear()
  {
    _eqm_ele_dir_index=0;
    _eo.clear();
    _ele_vols.clear();
    _pdos.clear();
    _edos.clear();
    _fermi_energies.clear();
    _pV.clear();
    _zpe.clear();
    _data_read_error=false;
    _nl_success_status=0;
    _nl_err_msg="";
    _cutoff_freq=0.0;
    _luncertanity_V0=0.0;
    _luncertanity_E0=0.0;
    _luncertanity_B0=0.0;
    _lchisq=0.0;
    _leqmV0=0.0;
    _leqmE0=0.0;
    _leqmB0=0.0;
    _TF.clear();
    _uncertanity_V0=0.0;
    _uncertanity_E0=0.0;
    _uncertanity_B0=0.0;
    _uncertanity_Bp=0.0;
    _uncertanity_Bpp=0.0;
    _chisq_dof=0.0; 
    _nleqmV0=0.0;
    _nleqmE0=0.0;
    _nleqmB0=0.0;
    _nleqmBp=0.0;
    _nleqmBpp=0.0;
    _fdfsolver_name="";
    _fitting_type="BM1";
    _Feqm=0.0, _Beqm=0.0, _Veqm=0.0, _Bp=0.0, _Bpp=0.0;
    _atomic_species.clear();
  }
  // ***************************************************
  void QHAEOS::set_fitting_type(string s)
  {
    _fitting_type=s;
  }
  // ***************************************************
  void QHAEOS::set_include_ele(bool b)
  {
    _include_ele=b;
  }
  // ***************************************************
  //initialize all variables with this function
  bool QHAEOS::setvariables()
  {
    _logger<<"QHAEOS importing variables, "<<apl::endl;
    _eqm_ele_dir_index=_qhen.get_eqm_ele_dir_index();
    _eo=_qhen.get_eo();
    _ele_vols=_qhen.get_ele_vols();
    _pdos=_qhen.get_pdos_data();
    _edos=_qhen.get_cedos_data();
    _fermi_energies=_qhen.get_fermi_energies();
    _pV=_qhen.get_pV();
    //_edosATfermi=_qhen.get_edosATfermi();
    _is_magnetic=_qhen.get_is_magnetic();
    _zpe.resize(_eo.size(),0.0);
    calculate_zpe();
    //_eweights=_qhen.get_eweights();
    //_eband_energies=_qhen.get_eband_energies();
    _atomic_species=_qhen.get_atomic_species();
    return check_size();
  }
  // ***************************************************
  bool QHAEOS::check_size()
  {
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
    if(_pdos.size()==0)
      {
	_logger << apl::error << "_pdos.size()==0 "<< apl::endl;
	return false;
      }
    if(_fermi_energies.size()==0)
      {
	_logger << apl::error << "_fermi_energies.size()==0 "<< apl::endl;
	return false;
      }
    if(_pV.size()==0)
      {
	_logger << apl::error << "_pV.size()==0 "<< apl::endl;
	return false;
      }
    if(_zpe.size()==0)
      {
	_logger << apl::error << "_zpe.size()==0 "<< apl::endl;
	return false;
      }
    if(!_is_magnetic)
      {
	if(_edos.size()==0)
	  {
	    _logger << apl::error << "_edos.size()==0 "<< apl::endl;
	    return false;
	  }
      }
    if(_atomic_species.size()==0)
      {
	_logger << apl::error << "_atomic_species()==0 "<< apl::endl;
	return false;
      }
    return true;
  }
  // ***************************************************************************************
  //calculate EOS for a range of temperatures
  void QHAEOS::cal_qheos(double USER_TP_TSTART, double USER_TP_TEND, double USER_TP_TSTEP, ThermalPropertiesCalculator &e)
  {
    std::string STAR50 = std::string(50, '*');
    std::string STAR180 = std::string(180, '*');
    std::string STAR150 = std::string(150, '*');

    _logger << "Writing EOS into file aflow.qha.thermo.out"<<apl::endl;
    _logger << "Writing errors of fitting into file aflow.qha.err_fit.out"<<apl::endl;
    if(!_include_ele)_logger << "Writing EOS energies into file aflow.qha.FVT.out"<<apl::endl;

    //calculating zpe for all configurations

    stringstream osfvt, osfit, oseos;
    bool includePV=false;
    for(uint i=0; i!=_pV.size(); i++)if(!_iszero(_pV[i]))includePV=true;

    osfit<<"# Estimating Errors in Volume-Energy Curve Fitting\n";
    osfit<<"# lchisq    => linear fit chi square\n";
    osfit<<"# lERR(X)   => error associated with linear fit of variable X\n";
    osfit<<"# lmsder    => Levenberg-Marquardt algorithm\n";
    osfit<<"# chisq_dof => nonlinear fit chi square per degrees of freedom\n";
    osfit<<"# nERR(X)   => error associated with nonlinear fit of variable X\n";
    osfit<<"[AFLOW] "<<STAR180<<"\n";
    osfit << "[AFLOW_QHA_FIT_ERROR]START" <<"\n";
    osfit <<setprecision(2) << std::fixed << std::showpoint;
    osfit<<"#"<<aurostd::PaddedPRE("T(K)",14," ")
         <<aurostd::PaddedPRE("status",15," ")
         <<aurostd::PaddedPRE("lchisq",15," ")
         <<aurostd::PaddedPRE("lERR(V)",17," ")
         <<aurostd::PaddedPRE("lERR(E)",14," ")
         <<aurostd::PaddedPRE("lERR(B)",15," ")
         <<aurostd::PaddedPRE("name",12," ")
         <<aurostd::PaddedPRE("chisq_dof",19," ")
         <<aurostd::PaddedPRE("nERR(V)",15," ")
         <<aurostd::PaddedPRE("nERR(E)",15," ")
         <<aurostd::PaddedPRE("nERR(B)",15," ")
         <<aurostd::PaddedPRE("nERR(Bp)",15," ")
         <<"\n";

    oseos<<"[AFLOW] "<<STAR150<<"\n";
    oseos<<"# U     => vibrational internal energy\n";
    oseos<<"# F     => vibrational free energy\n";
    oseos<<"# S     => vibrational energy \n";
    oseos<<"# Fqh   => quasi-harmonic free energy \n";
    oseos<<"# V     => volume \n";
    oseos<<"# B     => Bulk Modulus \n";
    oseos<<"# Bp    => Pressure derivative Bulk Modulus \n";
    oseos<<"# gamma => Gruneisen Parameter \n";
    oseos<<"# alpha => Volume expansion cofficient x 10^6 \n";
    //oseos<<"# Bpp => Pressure derivative of Bp  \n";
    if(_include_ele){
      oseos<<"# cv_e => electronic cv \n";
      oseos<<"# alpha_e => electronic volume expansion cofficient x 10^6 \n";
    }

    oseos<<"[AFLOW] "<<STAR150<<"\n";
    oseos << "[AFLOW_QHA_THERMO]START" <<"\n";
    oseos<<"#"<<aurostd::PaddedPRE("T(K)",6," ")
         <<aurostd::PaddedPRE("U(eV/cell)",19," ")
         <<aurostd::PaddedPRE("F(eV/cell)",15," ")
         <<aurostd::PaddedPRE("S(kB/cell)",14," ")
         <<aurostd::PaddedPRE("Fqh(eV/cell)",16," ")
         <<aurostd::PaddedPRE("V(A^3)",10," ")
         <<aurostd::PaddedPRE("B(GPa)",15," ")
         <<aurostd::PaddedPRE("Bp",11," ")
         <<aurostd::PaddedPRE("gamma",19," ")
         <<aurostd::PaddedPRE("alpha (1/K)",20," ")
         <<aurostd::PaddedPRE("Cv(kB/cell)",15," ")
         <<aurostd::PaddedPRE("Cp(kB/cell)",15," ");
    if(_include_ele){
      oseos<<aurostd::PaddedPRE("Cv_e(kB/cell)",15," ");
      oseos<<aurostd::PaddedPRE("alpha_e (1/K)",15," ");
    }

    if(_fitting_type=="BM3")
      {
        oseos <<aurostd::PaddedPRE("Bpp(1/GPa)",15," ")<<"\n";
      }else{
      oseos<<"\n";
    }

    for(double K=(double)USER_TP_TSTART; K<=(double)USER_TP_TEND; K+=USER_TP_TSTEP) //temperature loop
      {
	if(!_include_ele){osfvt<<"[AFLOW] "<<STAR50<<"\n";
	  osfvt << "[AFLOW_QHA_ENERGIES  T="<<setprecision(2)<<std::fixed<<K<<" K ]START" <<"\n";
        }
	xvector<double> Ftotal(_eo.size(), 1);
	xvector<double> volume(_eo.size(), 1);

        if(!_include_ele){
	  osfvt<<"#"<<setw(15)<<aurostd::PaddedPRE("V(A^3)",12," ")
	       <<aurostd::PaddedPRE("Ftot(eV/Cell)",15," ")
	       <<aurostd::PaddedPRE("E0K(eV/Cell)",15," ")
	       <<aurostd::PaddedPRE("Fvib(eV/Cell)",15," ");
        }
	if(includePV && !_include_ele)osfvt<<aurostd::PaddedPRE("pV(eV/Cell)",15," ")<<"\n";
	else osfvt<<"\n";

	for(uint i=0; i!=_eo.size(); i++)
	  {
            double vib_i=VibrationEnergy(K, i);

	    if(includePV)Ftotal[i+1]=_eo[i]+vib_i+_pV[i];
	    else         Ftotal[i+1]=_eo[i]+vib_i;

	    volume[i+1]=_ele_vols[i];

            if(!_include_ele){
	      osfvt <<setw(15)<<std::setprecision(8)  <<volume[i+1]
		    <<setw(15)<<std::setprecision(8)  <<Ftotal[i+1]
		    <<setw(15)<<std::setprecision(8)  <<_eo[i]
		    <<setw(15)<<std::setprecision(8)  <<vib_i;
	      if(includePV)osfvt<<setw(15)<<std::setprecision(8)  <<_pV[i];
	      osfvt<<"\n";
            }
	  }
        if(!_include_ele){
	  osfvt << "[AFLOW_QHA_ENERGIES  T="<<setprecision(2)<<std::fixed<<K<<" K ]END" <<"\n";
	  osfvt<<"[AFLOW] "<<STAR50<<"\n";
        }
	//fitting
	initialize_output_variables();
	md_lsquares_call(volume, Ftotal);
        double error_sum=std::abs(_uncertanity_V0)+std::abs(_uncertanity_E0)+std::abs(_uncertanity_B0)+std::abs(_uncertanity_Bp); 

	//calculate CTE and Cp
	double avgGP=_qha.average_gruneisen_parameter(K);
	double cv=getIsochoricSpecificHeat(K, _eqm_ele_dir_index);
	double CTE=((avgGP)*cv)/(_Beqm*_Veqm);
	double Cp=cv + CTE* CTE * _Beqm* _Veqm* K;

	double THERMO_Ut=e.getInternalEnergy(K,apl::meV);
	double THERMO_Ft=e.getVibrationalFreeEnergy(K,apl::meV);
	double THERMO_St=e.getVibrationalEntropy(K, apl::kB);
	double THERMO_Cvt=e.getIsochoricSpecificHeat(K, apl::kB);

        //electronic cv 
	double ele_cv=0.0;
	double ele_CTE=0.0;
	if(_include_ele){
	  ele_cv=Electronic_Cv(K,_eqm_ele_dir_index);
	  ele_CTE=(2.0/(3.0*_Beqm*_Veqm))*ele_cv;
	}

        oseos <<setw(8) << setprecision(2)<<std::fixed << std::showpoint << K
              <<setw(15)<< setprecision(6)<<std::fixed << std::showpoint << THERMO_Ut/1000.00
              <<setw(15)<< setprecision(6)<<std::fixed << std::showpoint << THERMO_Ft/1000.00
              <<setw(15)<< setprecision(6)<<std::fixed << std::showpoint << THERMO_St
              <<setw(15)<< setprecision(6)<<std::fixed << std::showpoint << _Feqm
              <<setw(15)<<setprecision(6) <<std::fixed << std::showpoint << _Veqm
              <<setw(15)<<setprecision(6) <<std::fixed << std::showpoint << _Beqm*160.2176487
              <<setw(15)<<setprecision(6) <<std::fixed << std::showpoint << _Bp
              <<setw(15)<<setprecision(6) <<std::fixed << std::showpoint << avgGP
              <<setw(15)<<setprecision(6)<<std::fixed << std::showpoint << CTE*1.0E6
              <<setw(15)<<setprecision(6)<<std::fixed << std::showpoint << THERMO_Cvt 
              <<setw(15)<<setprecision(6)<<std::fixed << std::showpoint << (Cp/0.000086173324); //converted to Kb/cell

	if(_include_ele){
	  oseos <<setw(15)<<setprecision(6)<<std::fixed << std::showpoint << ele_cv/0.000086173324; //converted to Kb/cell
	  oseos <<setw(15)<<setprecision(6)<<std::fixed << std::showpoint << ele_CTE*1.0E6;
	}

        if(_fitting_type=="BM3")
          oseos<<setw(15)<<setprecision(6)<<std::fixed << std::showpoint << _Bpp/160.2176487<<"\n";
        else oseos<<"\n";
	//printing all output
	if(!_data_read_error){
	  if((!std::isnan(error_sum)) || (error_sum <allowed_fit_error)){
	    osfit<<setw(15)<<K<<setw(15)<<"sucess";
	    osfit<<setw(15)<<_lchisq <<setw(15)<<_luncertanity_V0<<setw(15)<<_luncertanity_E0<<setw(15)<<_luncertanity_B0;
	    osfit << setw(15) << _fdfsolver_name <<setw(15)<<_chisq_dof<<setw(15)<<_uncertanity_V0
		  <<setw(15)<<_uncertanity_E0<<setw(15)<<_uncertanity_B0<<setw(15)<<_uncertanity_Bp;
	    osfit<<"\n";
	  }else{
	    osfit<<"WARNING at T= "<< K <<" data is not well behaved "<<"\n";
	    _logger << apl::warning<<" problem with FVT data at T= "<<K << "[K] check  err_fit.out "<<apl::endl;
	  }
	}else{
	  osfit<<"error in data reading"<<"\n"; 
	  apl::APLRuntimeError("problem with volume and energy data");
	}
       vector<double> tmp(2, 0);
       tmp[0]=K;
       tmp[1]=_Feqm;
      _TF.push_back(tmp);
      }//K loop

    osfit << "[AFLOW_QHA_FIT_ERROR]END" <<"\n";
    osfit<<"[AFLOW] "<<STAR180<<"\n";
    oseos << "[AFLOW_QHA_THERMO]END" <<"\n";
    oseos<<"[AFLOW] "<<STAR150<<"\n";

    if(!_include_ele){
      string FVTfile =  "aflow.qha.FVT.out";
      if(!aurostd::stringstream2file(osfvt, FVTfile, "WRITE")) {
	throw APLRuntimeError("Cannot write aflow.qha.FVT.out");
      }
      aurostd::StringstreamClean(osfvt);
    }
    string fit_out =  "aflow.qha.err_fit.out";
    if(!aurostd::stringstream2file(osfit, fit_out, "WRITE")) {
      throw APLRuntimeError("Cannot write aflow.qha.err_fit.out");
    }
    aurostd::StringstreamClean(osfit);

    string eos_out =  "aflow.qha.thermo.out";
    if(!aurostd::stringstream2file(oseos, eos_out, "WRITE")) {
      throw APLRuntimeError("Cannot write aflow.qha.thermo.out");
    }
    aurostd::StringstreamClean(oseos);

    //total enthalpy calculations
    if(!_include_ele)total_enthalpy();
    else{ 
      enthalpy_incuding_ele(USER_TP_TSTART, USER_TP_TEND, USER_TP_TSTEP);
    }
  }//fn end
  // ***************************************************
  //get vibrational energies for proper file indix
  double QHAEOS::VibrationEnergy(double temperature_in_kelvins, uint dir_index)
  {
    if( temperature_in_kelvins < _AFLOW_APL_EPS_ ) return _zpe[dir_index];

    double steps=_pdos[dir_index][1][0]-_pdos[dir_index][0][0];
    double f = 0.0;
    double beta = 1.0 / ( 0.0861734315 * temperature_in_kelvins ); // beta = 1/kBT = 1 / ([meV/K] [K])
    for(uint i = 0; i != _pdos[dir_index].size(); i++)
      {
        double hni = 4.1356673310 * _pdos[dir_index][i][0]; // hplanck in [eV.s] * 10^-15 * freq in [Hz] * 10^12 => hni in [meV].
        f += _pdos[dir_index][i][1] * aurostd::ln( 1.0 - exp( -beta * hni ) ) / beta ;
      }
    return ( _zpe[dir_index] + ( f * steps )/1000.0 );
  }
  // ***************************************************
  void QHAEOS::calculate_zpe()
  {
    for(uint i=0; i!=_pdos.size(); i++)
      {
	double sum=0.0;
	double d_omega=_pdos[i][1][0]-_pdos[i][0][0];
	for(uint j=0; j!=_pdos[i].size(); j++)
	  {
	    sum+=_pdos[i][j][1]*_pdos[i][j][0];
	  }
	sum*=0.5*d_omega*4.135667662*1.0E-3;//1/2 * h * nu
	_zpe[i]=sum;
      }
  }
  // ***************************************************
  //calculate zero point energies at each distorted volumes point
  void QHAEOS::getZeroPointVibrationEnergy()
  {
    for(uint i=0; i!=_pdos.size(); i++)
      {
	double sum=0.0;
	double steps=_pdos[i][1][0]-_pdos[i][0][0];
	for(uint j=0; j!=_pdos[i].size(); j++){
	  sum += _pdos[i][j][0] * _pdos[i][j][1];
	}
	sum *= steps * 0.5 * 4.1356651596736798;//  0.5*h*nu
	_zpe[i]=sum/1000.0;
      }
  }
  // ***************************************************
  //get electronic energies for a proper file index
  double QHAEOS::ElectronicEnergy(double temperature_in_kelvins, uint dir_index)
  {
    if(temperature_in_kelvins<1.0) return 0.0;
    vector<vector<double> > edos=_edos[dir_index];
    double fermi = _fermi_energies[dir_index];
    double sum1=0.0;
    double sum2=0.0;
    double sum3=0.0;
            
    double dE=edos[1][0]-edos[0][0];
    for(uint j=0; j<edos.size(); j++)
      {
        double Ediff=edos[j][0]-fermi;
        double f= fermi_dirac_distribution(Ediff, temperature_in_kelvins);
        //continue if edos is .gt. 1.0e-5
        if(edos[j][1]<1.0e-5)continue;
        if(edos[j][0]<=fermi)sum2+=edos[j][1]*edos[j][0];
        //continue if f is .gt. 1.0e-6
        if(f<1.0e-6)continue;
        sum1+=edos[j][1]*f*edos[j][0];

        if(!_isequal(f,1.0))
          {
            double x=f*log(f)+(1.0-f)*log(1.0-f);
            sum3+=edos[j][1]*x;
          }
      }
    sum1*=dE;
    sum2*=dE;
    sum3*=8.6173303*1.0e-5*temperature_in_kelvins*dE;
    return (sum1-sum2+sum3);
  } 
  // ***************************************************
  //get electronic specific heat for a proper file index
  double QHAEOS::Electronic_Cv(double temperature_in_kelvins, uint dir_index)
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
  // ***************************************************
  //Fermi-Dirac distribution function
  double QHAEOS::fermi_dirac_distribution(const double delE, const double t)
  {
    double beta=1.0/(_kBeV*t*1.0e-5);
    double x=(delE)*beta;
    x=exp(x)+1.0;
    return 1.0/x;
  }
  // ***************************************************
  //initialize all fittiung variables to zero
  void QHAEOS::initialize_output_variables()
  {
    _data_read_error=false;
    _nl_success_status=0;
    _nl_err_msg="";
    _cutoff_freq=0.0;
    _luncertanity_V0=0.0;
    _luncertanity_E0=0.0;
    _luncertanity_B0=0.0;
    _lchisq=0.0;
    _leqmV0=0.0;
    _leqmE0=0.0;
    _leqmB0=0.0;
    _uncertanity_V0=0.0;
    _uncertanity_E0=0.0;
    _uncertanity_B0=0.0;
    _uncertanity_Bp=0.0;
    _uncertanity_Bpp=0.0;
    _chisq_dof=0.0;
    _nleqmV0=0.0;
    _nleqmE0=0.0;
    _nleqmB0=0.0;
    _nleqmBp=0.0;
    _nleqmBpp=0.0;
    _fdfsolver_name="";
    _Feqm=0.0, _Beqm=0.0, _Veqm=0.0, _Bp=0.0, _Bpp=0.0;
  }
  // ***************************************************
  //it will call non-linear fitting fuctions for fitting
  void QHAEOS::md_lsquares_call(const xvector<double> &V, const xvector<double> &E)
  {
    md_lsquares mdfit;
    mdfit.clear();
    for(int i=1; i<=V.rows; i++)
      {
	mdfit.Xdata.push_back(V[i]);
	mdfit.Ydata.push_back(E[i]);
      }
    mdfit.cubic_polynomial_fit(); // it calls both linear and nonlinear fit functions
    _data_read_error=mdfit.data_read_error;
    _nl_success_status=mdfit.nl_success_status;
    _nl_err_msg=mdfit.nl_err_msg;
    _luncertanity_V0=mdfit.luncertanity_V0;
    _luncertanity_E0=mdfit.luncertanity_E0;
    _luncertanity_B0=mdfit.luncertanity_B0;
    _lchisq=mdfit.lchisq;
    _leqmV0=mdfit.leqmV0;
    _leqmE0=mdfit.leqmE0;
    _leqmB0=mdfit.leqmB0;
    _uncertanity_V0=mdfit.uncertanity_V0;
    _uncertanity_E0=mdfit.uncertanity_E0;
    _uncertanity_B0=mdfit.uncertanity_B0;    
    _uncertanity_Bp=mdfit.uncertanity_Bp;
    _chisq_dof=mdfit.chisq_dof; //chi square per degrees of freedom    
    _nleqmV0=mdfit.nleqmV0;
    _nleqmE0=mdfit.nleqmE0;
    _nleqmB0=mdfit.nleqmB0;
    _nleqmBp=mdfit.nleqmBp;
    _fdfsolver_name=mdfit.fdfsolver_name;
    _Feqm=_nleqmE0, _Beqm=_nleqmB0, _Veqm=_nleqmV0, _Bp=_nleqmBp;

    //more refinements in fitting
    xvector<double> guess(4,1);
    guess[1]=_nleqmE0;
    guess[2]=_nleqmB0;
    guess[3]=_nleqmV0;
    guess[4]=_nleqmBp;
    xvector<double> out(4,1);
    more_refinement(E, V, guess, out);
    //initial guess to BM 4th order fit from previous fit
    if(_fitting_type=="BM3"){
      xvector<double> testguess(4,1);
      testguess[1]=_nleqmV0;
      testguess[2]=_nleqmE0;
      testguess[3]=_nleqmB0;
      testguess[4]=_nleqmBp;
      mdfit.birch_murnaghan_4th_order_fit(testguess);
      _nleqmV0  = mdfit.nleqmV0;
      _nleqmE0  = mdfit.nleqmE0;
      _nleqmB0  = mdfit.nleqmB0;
      _nleqmBp  = mdfit.nleqmBp;
      _nleqmBpp = mdfit.nleqmBpp;
      _Feqm=_nleqmE0, _Beqm=_nleqmB0, _Veqm=_nleqmV0, _Bp=_nleqmBp, _Bpp=_nleqmBpp;
      _uncertanity_V0  =mdfit.uncertanity_V0;
      _uncertanity_E0  =mdfit.uncertanity_E0;
      _uncertanity_B0  =mdfit.uncertanity_B0;    
      _uncertanity_Bp  =mdfit.uncertanity_Bp;
      _uncertanity_Bpp =mdfit.uncertanity_Bpp;
      _chisq_dof=mdfit.chisq_dof;
    } 
    if(_fitting_type=="BM2"){
      xvector<double> testguess(4,1);
      testguess[1]=_nleqmV0;
      testguess[2]=_nleqmE0;
      testguess[3]=_nleqmB0;
      testguess[4]=_nleqmBp;
      mdfit.birch_murnaghan_3rd_order_fit(testguess);
      _nleqmV0  = mdfit.nleqmV0;
      _nleqmE0  = mdfit.nleqmE0;
      _nleqmB0  = mdfit.nleqmB0;
      _nleqmBp  = mdfit.nleqmBp;
      _nleqmBpp = mdfit.nleqmBpp;
      _Feqm=_nleqmE0, _Beqm=_nleqmB0, _Veqm=_nleqmV0, _Bp=_nleqmBp;
      _uncertanity_V0  =mdfit.uncertanity_V0;
      _uncertanity_E0  =mdfit.uncertanity_E0;
      _uncertanity_B0  =mdfit.uncertanity_B0;    
      _uncertanity_Bp  =mdfit.uncertanity_Bp;
      _chisq_dof=mdfit.chisq_dof; 
    }
    mdfit.clear();
  }
  // ***************************************************
  //doing more refinments after fitting 
  bool QHAEOS::more_refinement(const xvector<double> &E, const xvector<double> &V, xvector<double> &guess, xvector<double> &out)
  {
    apl::aflowFITTING fit;

    if(!fit.birch_murnaghan_fitting (E, V, guess, out))return false;
    _Feqm  =out[1];
    _Beqm  =out[2];
    _Veqm  =out[3];
    _Bp    =out[4];

    chisq_quadratic=fit.get_chisq_quadratic();
    chisq_birch_murnaghan=fit.get_chisq_birch_murnaghan();
    Iteration_birch_murnaghan=fit.getIteration_birch_murnaghan();
    alamda_birch_murnaghan=fit.getalamda_birch_murnaghan();
    Uncertainties_birch_murnaghan=fit.getUncertainties_birch_murnaghan();
    fit.clear();
    return true;
  }
  // ***************************************************
  //get specific heat for proper file-index
  double QHAEOS::getIsochoricSpecificHeat(double temperature_in_kelvins, uint dir_index)
  {
    double cv = 0.0;
    double beta = 1.0 / ( 0.0861734315 * temperature_in_kelvins ); // beta = 1/kBT = 1 / ([meV/K] [K])
    double steps=_pdos[dir_index][1][0]-_pdos[dir_index][0][0];
    for(uint i = 0; i < _pdos[dir_index].size(); i++)
      {
        double bhni = beta * 4.1356673310 * _pdos[dir_index][i][0];
        double ebhni = exp(bhni);
        cv += _pdos[dir_index][i][1] * 0.0861734315 * bhni * bhni / ( ( 1.0 - 1.0 / ebhni ) * ( ebhni - 1.0 ) );
      }
    if( isnan(cv) ) return 0.0;
    return (cv * steps)*0.001; //[meV/K]->[eV/K]
  }
  // ***************************************************
  void QHAEOS::total_enthalpy()
  {
    _logger<<"Writing aflow.qha.enthalpy.out file"<<apl::endl;

    std::string STAR50 = std::string(50, '*');
    stringstream outfile;
    outfile << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
    outfile << setprecision(4);

    outfile<<"[AFLOW] "<< STAR50 <<"\n";
    outfile<<"# T     => temperature [K]\n";
    outfile<<"# H     => total enthalpy [eV/cell]\n";
    outfile<<"[AFLOW] "<< STAR50 <<"\n";
    outfile<<"[AFLOW_QHA_ENTHALPY]START"<<"\n";
    outfile<<"#"<<setw(15)<<"T(K)"<<setw(15)<<"H(eV/cell)"<<"\n";
    /////////////////////////////////////////////////////////////////
    vector< vector <double> > data = _TF; _TF.clear(); 
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

    outfile<<"[AFLOW_QHA_ENTHALPY]END"<<"\n";
    outfile<<"[AFLOW] "<<STAR50<<"\n";
    string file="aflow.qha.enthalpy.out";
    if(!aurostd::stringstream2file(outfile, file, "WRITE")) {
      throw APLRuntimeError("Cannot write aflow.qha.enthalpy.out");
    }
    aurostd::StringstreamClean(outfile);
    data.clear();
  }
  // ***************************************************
  void QHAEOS::enthalpy_incuding_ele(double USER_TP_TSTART, double USER_TP_TEND, double USER_TP_TSTEP)
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
	    double vib_i=VibrationEnergy(K, i);
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
	initialize_output_variables();
	md_lsquares_call(volume, Ftotal);
        vector<double> tmp(2,0);
        tmp[0]=K;
        tmp[1]=_Feqm;
        data.push_back(tmp);


        osfvt << "[AFLOW_QHA_ENERGIES  T="<<setprecision(2)<<std::fixed<<K<<" K ]END" <<"\n";
        osfvt<<"[AFLOW] "<<STAR80<<"\n";
      }

    string FVTfile =  "aflow.qha.FVT.out";
    if(!aurostd::stringstream2file(osfvt, FVTfile, "WRITE")) {
      throw APLRuntimeError("Cannot write aflow.qha.FVT.out");
    }
    aurostd::StringstreamClean(osfvt);


    _logger<<"Writing aflow.qha.enthalpy.out file"<<apl::endl;

    std::string STAR50 = std::string(50, '*');
    stringstream outfile;
    outfile << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
    outfile << setprecision(4);

    outfile<<"[AFLOW] "<< STAR50 <<"\n";
    outfile<<"# T     => temperature [K]\n";
    outfile<<"# H     => total enthalpy [eV/cell]\n";
    outfile<<"[AFLOW] "<< STAR50 <<"\n";
    outfile<<"[AFLOW_QHA_ENTHALPY]START"<<"\n";
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

    outfile<<"[AFLOW_QHA_ENTHALPY]END"<<"\n";
    outfile<<"[AFLOW] "<<STAR50<<"\n";
    string file="aflow.qha.enthalpy.out";
    if(!aurostd::stringstream2file(outfile, file, "WRITE")) {
      throw APLRuntimeError("Cannot write aflow.qha.enthalpy.out");
    }
    aurostd::StringstreamClean(outfile);
    data.clear();
  }
  // ***************************************************
  bool QHAEOS::exists_test0 (const std::string& name)
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
  // ***************************************************
}//apl end
