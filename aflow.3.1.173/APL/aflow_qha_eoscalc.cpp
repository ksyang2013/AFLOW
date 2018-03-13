// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *                Aflow PINKU NATH - Duke University 2014-2017             *
// *                                                                         *
// ***************************************************************************
// Written by Pinku Nath
// pn49@duke.edu
//this cpp file generate EOS data

#include "aflow_apl.h"
using namespace std;
#define _isnegative(a) (a < MIN_FREQ_TRESHOLD) ? true : false
namespace apl {
// ***************************************************
EOS::EOS(IPhononCalculator &pc, eos_calculate &run, Logger &l) : QuasiHarmonicGruneisen(pc, run, l) {
  phdir.clear();
  eosdir.clear();
  E0K.clear();
  pV.clear();
  zpe.clear();
  fermi.clear();
  eosvol.clear();
  pDOS.clear();
  eDOS.clear();
  NEGATIVE_FREQ.clear();
  _ismagnetic = false;
}
EOS::~EOS() { clear(); }
void EOS::clear() {
  phdir.clear();
  eosdir.clear();
  E0K.clear();
  pV.clear();
  zpe.clear();
  fermi.clear();
  eosvol.clear();
  pDOS.clear();
  eDOS.clear();
  NEGATIVE_FREQ.clear();
}
// ***************************************************
//initialize all variables with this function
bool EOS::setvariables() {
  _logger << "calculating EOS " << apl::endl;
  dirfrefix = QuasiHarmonicGruneisen::dirfrefix;
  phdir = QuasiHarmonicGruneisen::phdir;
  eosdir = QuasiHarmonicGruneisen::eosdir;
  eosvol = QuasiHarmonicGruneisen::eosvol;
  CUTOFF_FREQ = QuasiHarmonicGruneisen::CUTOFF_FREQ;
  if ((phdir.size() != 0) && (eosdir.size() != 0) && (eosvol.size() != 0)) {
    ZERO_STATIC_DIR_INDEX = QuasiHarmonicGruneisen::ZERO_STATIC_DIR_INDEX;
    if (!calculatePDOS()) return false;  //vibrartional energies of each configuration
    calculateE0K();                      //Total energies at 0K of each configurations
    geteDOS();                           //electronic energies of each configurations
    if (_ismagnetic) _logger << apl::warning << " COMPOUND is magnetic, and electronic energy calculations skipped " << apl::endl;
    if (NEGATIVE_FREQ.size() != 0) resize_EOSnPH();
    getZeroPointVibrationEnergy();
  } else {
    _logger << apl::error << " EOS::setvariables() either dir of volume size is zero" << apl::endl;
    return false;
  }
  return true;
}
// ***************************************************
//calculate EOS for a range of temperatures
void EOS::calEOS(double USER_TP_TSTART, double USER_TP_TEND, double USER_TP_TSTEP, ThermalPropertiesCalculator &e) {
  std::string STAR80 = std::string(80, '*');
  std::string STAR130 = std::string(130, '*');
  std::string STAR180 = std::string(180, '*');
  std::string STAR150 = std::string(150, '*');

  _logger << "Writing EOS into file aflow.apl.eos.out" << apl::endl;
  _logger << "Writing errors of fitting into file aflow.apl.err_fit.out" << apl::endl;
  _logger << "Writing EOS energies into file aflow.apl.FVT.out" << apl::endl;

  stringstream osfvt, osfit, oseos;
  bool includePV = false;

  {  //properties at 300 K
    oseos << setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
    oseos << setprecision(6);
    oseos << "[AFLOW] " << STAR150 << "\n";
    oseos << "[APL_EOS_300K]START"
          << "\n";

    xvector<double> Ftotal(eDOS.size(), 1);
    xvector<double> volume(eDOS.size(), 1);
    for (uint i = 0; i != eDOS.size(); i++) {
      if (!_ismagnetic)
        Ftotal[i + 1] = E0K[i] + VibrationEnergy(300.00, i) + ElectronicEnergy(300.00, i) + pV[i];
      else
        Ftotal[i + 1] = E0K[i] + VibrationEnergy(300.00, i) + pV[i];
      volume[i + 1] = eosvol[i];
    }

    for (uint i = 0; i != pV.size(); i++)
      if (!_iszero(pV[i])) includePV = true;
    //performing fitting
    initialize_output_variables();
    md_lsquares_call(volume, Ftotal);

    //calculate EOS at 300K
    double avgGP = QuasiHarmonicGruneisen::averageGP(300.0);
    double cv = getIsochoricSpecificHeat(300.0, ZERO_STATIC_DIR_INDEX);
    double CTE = (avgGP * cv) / (Beqm * Veqm);
    double Cp = cv + CTE * CTE * Beqm * Veqm * 300.0;

    double U0 = e.getZeroPointVibrationEnergy(apl::meV);
    double U = e.getInternalEnergy(300, apl::meV);
    double Fvib = e.getVibrationalFreeEnergy(300, apl::meV);
    double Svib = e.getVibrationalEntropy(300, apl::kB);
    double Cv = e.getIsochoricSpecificHeat(300, apl::kB);

    //oseos<<"#cutoff frequency used = "<<CUTOFF_FREQ<<" amu" <<"\n";
    //os_thermo<<setprecision(8) << std::fixed << std::showpoint;
    oseos << "apl_zero_point_energy = " << U0 / 1000.00 << " (eV/cell)"
          << "\n";
    oseos << "apl_internal_energy_300K = " << U / 1000.00 << " (eV/cell)"
          << "\n";
    oseos << "apl_vibrational_free_energy_300K = " << Fvib / 1000.00 << " (eV/cell)"
          << "\n";
    oseos << "apl_vibrational_entropy_300K = " << Svib << " (kB/cell)"
          << "\n";
    oseos << "apl_heat_capacity_Cv_300K = " << Cv << " (kB/cell)"
          << "\n";
    oseos << "apl_bulk_modulus_isothermal_300K = " << Beqm * 160.2176487 << " (GPa)"
          << "\n";
    oseos << "apl_heat_capacity_Cp_300K  = " << Cp / 0.000086173324 << " (kB/cell)"
          << "\n";  // converting in kB/cell
    oseos << "apl_volume_expansion_cofficient_300K  = " << CTE * 1E5 << " (10^5 * 1/K)"
          << "\n";
    oseos << "[APL_EOS_300K]STOP"
          << "\n";
    oseos << "[AFLOW] " << STAR130 << "\n";
  }
  osfit << "# Estimating Errors in Volume-Energy Curve Fitting\n";
  osfit << "# lchisq    => linear fit chi square\n";
  osfit << "# lERR(X)   => error associated with linear fit of variable X\n";
  osfit << "# lmsder    => Levenberg-Marquardt algorithm\n";
  osfit << "# chisq_dof => nonlinear fit chi square per degrees of freedom\n";
  osfit << "# nERR(X)   => error associated with nonlinear fit of variable X\n";
  osfit << "[AFLOW] " << STAR180 << "\n";
  osfit << "[APL_FIT_ERROR]START"
        << "\n";
  osfit << setprecision(2) << std::fixed << std::showpoint;
  osfit << "#" << aurostd::PaddedPRE("T(K)", 14, " ")
        << aurostd::PaddedPRE("status", 15, " ")
        << aurostd::PaddedPRE("lchisq", 15, " ")
        << aurostd::PaddedPRE("lERR(V)", 17, " ")
        << aurostd::PaddedPRE("lERR(E)", 14, " ")
        << aurostd::PaddedPRE("lERR(B)", 15, " ")
        << aurostd::PaddedPRE("name", 12, " ")
        << aurostd::PaddedPRE("chisq_dof", 19, " ")
        << aurostd::PaddedPRE("nERR(V)", 15, " ")
        << aurostd::PaddedPRE("nERR(E)", 15, " ")
        << aurostd::PaddedPRE("nERR(B)", 15, " ")
        << aurostd::PaddedPRE("nERR(Bp)", 15, " ")
        << "\n";

  oseos << "# U     => vibrational internal energy\n";
  oseos << "# F     => vibrational free energy\n";
  oseos << "# S     => vibrational energy \n";
  oseos << "# Fqh   => quasi-harmonic free energy \n";
  oseos << "# V     => volume \n";
  oseos << "# B     => Bulk Modulus \n";
  oseos << "# Bp    => Pressure derivative Bulk Modulus \n";
  oseos << "# gamma => Gruneisen Parameter \n";
  oseos << "# alpha => Volume expansion cofficient x 10^5 \n";
  oseos << "# Bpp => Pressure derivative of Bp  \n";
  oseos << "[AFLOW] " << STAR130 << "\n";
  if (FITTING_TYPE == "BM3") {
    oseos << "[AFLOW] "
          << "Birch-Murnaghan 4th-order function is used \n";
  } else if (FITTING_TYPE == "BM2") {
    oseos << "[AFLOW] "
          << "Birch-Murnaghan 3rd-order function is used \n";
  } else {
    oseos << "[AFLOW] "
          << "Murnaghan function is used \n";
  }
  oseos << "[APL_EOS]START"
        << "\n";
  oseos << "#" << aurostd::PaddedPRE("T(K)", 6, " ")
        << aurostd::PaddedPRE("U(eV/cell)", 19, " ")
        << aurostd::PaddedPRE("F(eV/cell)", 15, " ")
        << aurostd::PaddedPRE("S(kB/cell)", 14, " ")
        << aurostd::PaddedPRE("Fqh(eV/cell)", 16, " ")
        << aurostd::PaddedPRE("V(A^3)", 10, " ")
        << aurostd::PaddedPRE("B(GPa)", 15, " ")
        << aurostd::PaddedPRE("Bp", 11, " ")
        << aurostd::PaddedPRE("gamma", 19, " ")
        << aurostd::PaddedPRE("alpha (1/K)", 20, " ")
        << aurostd::PaddedPRE("Cv(kB/cell)", 15, " ")
        << aurostd::PaddedPRE("Cp(kB/cell)", 15, " ");
  if (FITTING_TYPE == "BM3") {
    oseos << aurostd::PaddedPRE("Bpp(1/GPa)", 15, " ") << "\n";
  } else {
    oseos << "\n";
  }

  for (double K = (double)USER_TP_TSTART; K <= (double)USER_TP_TEND; K += USER_TP_TSTEP)  //temperature loop
  {
    osfvt << "[AFLOW] " << STAR80 << "\n";
    osfvt << "[APL_EOS_ENERGIES  T=" << setprecision(2) << std::fixed << K << " K ]START"
          << "\n";
    xvector<double> Ftotal(eDOS.size(), 1);
    xvector<double> volume(eDOS.size(), 1);
    osfvt << "#" << setw(15) << aurostd::PaddedPRE("V(A^3)", 12, " ")
          << aurostd::PaddedPRE("Ftot(eV/Cell)", 15, " ")
          << aurostd::PaddedPRE("E0K(eV/Cell)", 15, " ")
          << aurostd::PaddedPRE("Fvib(eV/Cell)", 15, " ");
    if (!_ismagnetic) osfvt << aurostd::PaddedPRE("Fele(eV/Cell)", 15, " ");
    if (includePV)
      osfvt << aurostd::PaddedPRE("pV(eV/Cell)", 15, " ") << "\n";
    else
      osfvt << "\n";

    for (uint i = 0; i != eDOS.size(); i++) {
      Ftotal[i + 1] = E0K[i] + VibrationEnergy(K, i) + ElectronicEnergy(K, i) + pV[i];
      volume[i + 1] = eosvol[i];
      osfvt << setw(15) << std::setprecision(8) << volume[i + 1]
            << setw(15) << std::setprecision(8) << Ftotal[i + 1]
            << setw(15) << std::setprecision(8) << E0K[i]
            << setw(15) << std::setprecision(8) << VibrationEnergy(K, i);
      if (!_ismagnetic)
        osfvt << setw(15) << std::setprecision(8) << ElectronicEnergy(K, i);
      if (includePV)
        osfvt << setw(15) << std::setprecision(8) << pV[i] << "\n";
      else
        osfvt << "\n";
    }
    osfvt << "[APL_EOS_ENERGIES  " << setprecision(2) << K << " K ]STOP"
          << "\n";
    osfvt << "[AFLOW] " << STAR80 << "\n";

    //fitting
    initialize_output_variables();
    md_lsquares_call(volume, Ftotal);
    double error_sum = uncertanity_V0 + uncertanity_E0 + uncertanity_B0 + uncertanity_Bp;

    //calculate CTE and Cp
    double avgGP = QuasiHarmonicGruneisen::averageGP(K);
    double cv = getIsochoricSpecificHeat(K, ZERO_STATIC_DIR_INDEX);
    double CTE = (avgGP * cv) / (Beqm * Veqm);
    double Cp = cv + CTE * CTE * Beqm * Veqm * K;

    double THERMO_Ut = e.getInternalEnergy(K, apl::meV);
    double THERMO_Ft = e.getVibrationalFreeEnergy(K, apl::meV);
    double THERMO_St = e.getVibrationalEntropy(K, apl::kB);
    double THERMO_Cvt = e.getIsochoricSpecificHeat(K, apl::kB);

    oseos << setw(8) << setprecision(2) << std::fixed << std::showpoint << K
          << setw(15) << setprecision(6) << std::fixed << std::showpoint << THERMO_Ut / 1000.00
          << setw(15) << setprecision(6) << std::fixed << std::showpoint << THERMO_Ft / 1000.00
          << setw(15) << setprecision(6) << std::fixed << std::showpoint << THERMO_St
          << setw(15) << setprecision(6) << std::fixed << std::showpoint << Feqm
          << setw(15) << setprecision(6) << std::fixed << std::showpoint << Veqm
          << setw(15) << setprecision(6) << std::fixed << std::showpoint << Beqm * 160.2176487
          << setw(15) << setprecision(6) << std::fixed << std::showpoint << Bp
          << setw(15) << setprecision(6) << std::fixed << std::showpoint << avgGP
          << setw(15) << setprecision(6) << std::fixed << std::showpoint << CTE * 1E5
          << setw(15) << setprecision(6) << std::fixed << std::showpoint << THERMO_Cvt
          << setw(15) << setprecision(6) << std::fixed << std::showpoint << (Cp / 0.000086173324);  //converted to Kb/cell
    if (FITTING_TYPE == "BM3")
      oseos << setw(15) << setprecision(6) << std::fixed << std::showpoint << Bpp / 160.2176487 << "\n";
    else
      oseos << "\n";
    //printing all output
    if (!data_read_error) {
      if ((!std::isnan(error_sum)) || (error_sum < allowed_fit_error)) {
        osfit << setw(15) << K << setw(15) << "sucess";
        osfit << setw(15) << lchisq << setw(15) << luncertanity_V0 << setw(15) << luncertanity_E0 << setw(15) << luncertanity_B0;
        osfit << setw(15) << fdfsolver_name << setw(15) << chisq_dof << setw(15) << uncertanity_V0
              << setw(15) << uncertanity_E0 << setw(15) << uncertanity_B0 << setw(15) << uncertanity_Bp;
        osfit << "\n";
      } else {
        osfit << "WARNING at T= " << K << " data is not well behaved "
              << "\n";
        _logger << apl::warning << " problem with FVT data at T= " << K << "[K] check  err_fit.out " << apl::endl;
      }
    } else {
      osfit << "error in data reading"
            << "\n";
      apl::APLRuntimeError("problem with volume and energy data");
    }

  }  //K loop

  osfit << "[APL_FIT_ERROR]STOP"
        << "\n";
  osfit << "[AFLOW] " << STAR180 << "\n";
  oseos << "[APL_EOS]STOP"
        << "\n";
  oseos << "[AFLOW] " << STAR150 << "\n";

  string FVTfile = "aflow.apl.FVT.out";
  if (!aurostd::stringstream2file(osfvt, FVTfile, "WRITE")) {
    throw APLRuntimeError("Cannot write aflow.apl.FVT.out");
  }
  aurostd::StringstreamClean(osfvt);

  string fit_out = "aflow.apl.err_fit.out";
  if (!aurostd::stringstream2file(osfit, fit_out, "WRITE")) {
    throw APLRuntimeError("Cannot write aflow.apl.err_fit.out");
  }
  aurostd::StringstreamClean(osfit);

  string eos_out = "aflow.apl.eos.out";
  if (!aurostd::stringstream2file(oseos, eos_out, "WRITE")) {
    throw APLRuntimeError("Cannot write aflow.apl.eos.out");
  }
  aurostd::StringstreamClean(oseos);
}  //fn end
// ***************************************************
//initialize all fittiung variables to zero
void EOS::initialize_output_variables() {
  data_read_error = false;  //no error
  nl_success_status = 0;    //success
  nl_err_msg = "";
  luncertanity_V0 = 0.0;
  luncertanity_E0 = 0.0;
  luncertanity_B0 = 0.0;
  lchisq = 0.0;
  leqmV0 = 0.0;
  leqmE0 = 0.0;
  leqmB0 = 0.0;
  uncertanity_V0 = 0.0;
  uncertanity_E0 = 0.0;
  uncertanity_B0 = 0.0;
  uncertanity_Bp = 0.0;
  chisq_dof = 0.0;  //chi square per degrees of freedom
  nleqmV0 = 0.0;
  nleqmE0 = 0.0;
  nleqmB0 = 0.0;
  nleqmBp = 0.0;
  fdfsolver_name = "";
  Feqm = 0., Beqm = 0., Veqm = 0., Bp = 0.;
}
// ***************************************************
//it will call non-linear fitting fuctions for fitting
void EOS::md_lsquares_call(const xvector<double> &V, const xvector<double> &E) {
  md_lsquares mdfit;
  mdfit.clear();
  for (int i = 1; i <= V.rows; i++) {
    mdfit.Xdata.push_back(V[i]);
    mdfit.Ydata.push_back(E[i]);
  }
  mdfit.cubic_polynomial_fit();  // it calls both linear and nonlinear fit functions
  data_read_error = mdfit.data_read_error;
  nl_success_status = mdfit.nl_success_status;
  nl_err_msg = mdfit.nl_err_msg;
  luncertanity_V0 = mdfit.luncertanity_V0;
  luncertanity_E0 = mdfit.luncertanity_E0;
  luncertanity_B0 = mdfit.luncertanity_B0;
  lchisq = mdfit.lchisq;
  leqmV0 = mdfit.leqmV0;
  leqmE0 = mdfit.leqmE0;
  leqmB0 = mdfit.leqmB0;
  uncertanity_V0 = mdfit.uncertanity_V0;
  uncertanity_E0 = mdfit.uncertanity_E0;
  uncertanity_B0 = mdfit.uncertanity_B0;
  uncertanity_Bp = mdfit.uncertanity_Bp;
  chisq_dof = mdfit.chisq_dof;  //chi square per degrees of freedom
  nleqmV0 = mdfit.nleqmV0;
  nleqmE0 = mdfit.nleqmE0;
  nleqmB0 = mdfit.nleqmB0;
  nleqmBp = mdfit.nleqmBp;
  fdfsolver_name = mdfit.fdfsolver_name;
  Feqm = nleqmE0, Beqm = nleqmB0, Veqm = nleqmV0, Bp = nleqmBp;

  //more refinements in fitting
  xvector<double> guess(4, 1);
  guess[1] = nleqmE0;
  guess[2] = nleqmB0;
  guess[3] = nleqmV0;
  guess[4] = nleqmBp;
  xvector<double> out(4, 1);
  more_refinement(E, V, guess, out);
  //initial guess to BM 4th order fit from previous fit
  if (FITTING_TYPE == "BM3") {
    xvector<double> testguess(4, 1);
    testguess[1] = nleqmV0;
    testguess[2] = nleqmE0;
    testguess[3] = nleqmB0;
    testguess[4] = nleqmBp;
    mdfit.birch_murnaghan_4th_order_fit(testguess);
    nleqmV0 = mdfit.nleqmV0;
    nleqmE0 = mdfit.nleqmE0;
    nleqmB0 = mdfit.nleqmB0;
    nleqmBp = mdfit.nleqmBp;
    nleqmBpp = mdfit.nleqmBpp;
    Feqm = nleqmE0, Beqm = nleqmB0, Veqm = nleqmV0, Bp = nleqmBp, Bpp = nleqmBpp;
    uncertanity_V0 = mdfit.uncertanity_V0;
    uncertanity_E0 = mdfit.uncertanity_E0;
    uncertanity_B0 = mdfit.uncertanity_B0;
    uncertanity_Bp = mdfit.uncertanity_Bp;
    uncertanity_Bpp = mdfit.uncertanity_Bpp;
    chisq_dof = mdfit.chisq_dof;
  }
  if (FITTING_TYPE == "BM2") {
    xvector<double> testguess(4, 1);
    testguess[1] = nleqmV0;
    testguess[2] = nleqmE0;
    testguess[3] = nleqmB0;
    testguess[4] = nleqmBp;
    mdfit.birch_murnaghan_3rd_order_fit(testguess);
    nleqmV0 = mdfit.nleqmV0;
    nleqmE0 = mdfit.nleqmE0;
    nleqmB0 = mdfit.nleqmB0;
    nleqmBp = mdfit.nleqmBp;
    nleqmBpp = mdfit.nleqmBpp;
    Feqm = nleqmE0, Beqm = nleqmB0, Veqm = nleqmV0, Bp = nleqmBp;
    uncertanity_V0 = mdfit.uncertanity_V0;
    uncertanity_E0 = mdfit.uncertanity_E0;
    uncertanity_B0 = mdfit.uncertanity_B0;
    uncertanity_Bp = mdfit.uncertanity_Bp;
    chisq_dof = mdfit.chisq_dof;
  }
  mdfit.clear();
}
// ***************************************************
//collecting fitting parameters
bool EOS::EVfit(const xvector<double> &E, const xvector<double> &V) {
  apl::aflowFITTING fit;
  xvector<double> P(3, 1);
  if (!fit.quadraticfit(E, V, P)) return false;
  double V0 = -P[2] / (2. * P[3]);
  double E0 = P[3] * V0 * V0 + P[2] * V0 + P[1];
  double B0 = 2. * P[3] * V0;
  xvector<double> gues(P.rows + 1, 1);
  gues[1] = E0;
  gues[2] = B0;
  gues[3] = V0;
  gues[4] = 4.0;
  xvector<double> P_bm(4, 1);

  if (!fit.birch_murnaghan_fitting(E, V, gues, P_bm)) return false;
  Feqm = P_bm[1];
  Beqm = P_bm[2];
  Veqm = P_bm[3];
  Bp = P_bm[4];

  chisq_quadratic = fit.get_chisq_quadratic();
  chisq_birch_murnaghan = fit.get_chisq_birch_murnaghan();
  Iteration_birch_murnaghan = fit.getIteration_birch_murnaghan();
  alamda_birch_murnaghan = fit.getalamda_birch_murnaghan();
  Uncertainties_birch_murnaghan = fit.getUncertainties_birch_murnaghan();
  fit.clear();
  return true;
}
// ***************************************************
//doing more refinments after fitting
bool EOS::more_refinement(const xvector<double> &E, const xvector<double> &V, xvector<double> &guess, xvector<double> &out) {
  apl::aflowFITTING fit;

  if (!fit.birch_murnaghan_fitting(E, V, guess, out)) return false;
  Feqm = out[1];
  Beqm = out[2];
  Veqm = out[3];
  Bp = out[4];

  chisq_quadratic = fit.get_chisq_quadratic();
  chisq_birch_murnaghan = fit.get_chisq_birch_murnaghan();
  Iteration_birch_murnaghan = fit.getIteration_birch_murnaghan();
  alamda_birch_murnaghan = fit.getalamda_birch_murnaghan();
  Uncertainties_birch_murnaghan = fit.getUncertainties_birch_murnaghan();
  fit.clear();
  return true;
}
// ***************************************************
//get specific heat for proper file-index
double EOS::getIsochoricSpecificHeat(double temperature_in_kelvins, uint dir_index) {
  double cv = 0.0;
  double beta = 1.0 / (0.0861734315 * temperature_in_kelvins);  // beta = 1/kBT = 1 / ([meV/K] [K])
  double steps = pDOS[dir_index][1][0] - pDOS[dir_index][0][0];
  for (uint i = 0; i < pDOS[dir_index].size(); i++) {
    double bhni = beta * 4.1356673310 * pDOS[dir_index][i][0];
    double ebhni = exp(bhni);
    cv += pDOS[dir_index][i][1] * 0.0861734315 * bhni * bhni / ((1.0 - 1.0 / ebhni) * (ebhni - 1.0));
  }
  if (isnan(cv)) return 0.0;
  return (cv * steps) * 0.001;  //[meV/K]->[eV/K]
}
// ***************************************************
//get vibrational energies for proper file indix
double EOS::VibrationEnergy(double temperature_in_kelvins, uint dir_index) {
  if (temperature_in_kelvins < _AFLOW_APL_EPS_) return zpe[dir_index];
  double steps = pDOS[dir_index][1][0] - pDOS[dir_index][0][0];
  double f = 0.0;
  double beta = 1.0 / (0.0861734315 * temperature_in_kelvins);  // beta = 1/kBT = 1 / ([meV/K] [K])
  for (uint i = 0; i != pDOS[dir_index].size(); i++) {
    double hni = 4.1356673310 * pDOS[dir_index][i][0];  // hplanck in [eV.s] * 10^-15 * freq in [Hz] * 10^12 => hni in [meV].
    f += pDOS[dir_index][i][1] * aurostd::ln(1.0 - exp(-beta * hni)) / beta;
  }
  return (zpe[dir_index] + (f * steps) / 1000.0);
}
// ***************************************************
//get electronic energies for proper file index
double EOS::ElectronicEnergy(double temperature_in_kelvins, uint dir_index) {
  if (temperature_in_kelvins < 1) return 0;

  long double epsilon = 1e-5;
  double kb = 8.6173324e-5;  //eV/K
  double ENERGY = 0.0, ENERGY_EF = 0.0, S_ele = 0.0;
  double fermi_level = fermi[dir_index];

  for (uint i = 0; i < eDOS[dir_index].size() - 1; i++)  //integration over DOS
  {
    double Eele_i = eDOS[dir_index][i][0];
    double Eele_ip = eDOS[dir_index][i + 1][0];
    double dos_i = eDOS[dir_index][i][1];
    if (dos_i < epsilon) continue;
    double delE = std::abs(Eele_ip - Eele_i);  //delE is not uniform in DOScar file
    double fi = fermi_dirac_distribution(Eele_i, fermi_level, temperature_in_kelvins);
    ENERGY += dos_i * fi * Eele_i * delE;
    S_ele += -dos_i * (fi * aurostd::ln(fi) + (1.0 - fi) * aurostd::ln(1.0 - fi)) * delE;
    if (Eele_i < fermi_level) ENERGY_EF += dos_i * Eele_i * delE;
  }
  if (std::isnan(S_ele)) S_ele = 0.00;

  return ((ENERGY - ENERGY_EF) + kb * temperature_in_kelvins * kb * S_ele);
}
// ***************************************************
//Fermi-Dirac distribution function
double EOS::fermi_dirac_distribution(double Ei, double fermi_energy, double temperature_in_kelvins) {
  long double kb = 8.6173324e-5;  //eV/K
  double f = exp((Ei - fermi_energy) / (kb * temperature_in_kelvins)) + 1;
  return 1. / f;
}
// ***************************************************
//calculate phonon density of states and exclude negative density of states files
bool EOS::calculatePDOS() {
  pDOS.clear();
  NEGATIVE_FREQ.clear();
  string file = "";
  vector<uint> rc(2, 0);
  rc[0] = 0;
  rc[1] = 3;
  //reading freq and pdos from PDOS file
  for (uint i = 0; i != phdir.size(); i++) {
    if (i == (uint)ZERO_STATIC_DIR_INDEX) {
      file = "PDOS";
      pDOS.push_back(readfile<double>(file, rc));
    }
    file = dirfrefix + "/PDOS." + phdir[i];
    pDOS.push_back(readfile<double>(file, rc));
  }

  //checking -ve frequency in PDOS file
  for (uint i = 0; i != pDOS.size(); i++) {
    for (uint j = 0; j != pDOS[i].size(); j++) {
      if (pDOS[i][j][0] < 0.0) {
        if (_isnegative(pDOS[i][j][0])) {
          NEGATIVE_FREQ.push_back(i);  //identify pDOS files with -ve frequencies
          break;
        } else
          pDOS[i][j][0] = 0.00;
      }
    }
  }

  //indefying -ve frequencies index in PDOS file
  if (ZERO_STATIC_DIR_INDEX < 0) {
    file = "PDOS";
    pDOS.push_back(readfile<double>(file, rc));
    ZERO_STATIC_DIR_INDEX = pDOS.size() - 1;
  }

  if (NEGATIVE_FREQ.size() != 0) {
    for (uint i = 0; i != NEGATIVE_FREQ.size(); i++) {
      uint index = NEGATIVE_FREQ[i];
      if ((ZERO_STATIC_DIR_INDEX > 0) && (NEGATIVE_FREQ[i] > (uint)ZERO_STATIC_DIR_INDEX))
        index = NEGATIVE_FREQ[i] - 1;
      file = dirfrefix + "/PDOS." + phdir[index];
      _logger << "Negative frequencies in " << file << " and excluded" << apl::endl;
    }
  }

  double WARNING = (double)(NEGATIVE_FREQ.size() / phdir.size());

  if (WARNING > 0.6) {
    _logger << "There is no enough data to calculate EOS" << apl::endl;
    return false;
  }
  return true;
}  //fn end
// ***************************************************
//get 0K energies from different directories
void EOS::calculateE0K() {
  E0K.clear();
  pV.clear();
  E0K.resize(eosdir.size(), 0.0);
  pV.resize(eosdir.size(), 0.0);
  MagCell.resize(eosdir.size(), 0.0);

  for (uint i = 0; i != eosdir.size(); i++) {
    string file = eosdir[i] + "/aflow.qmvasp.out";
    double pv = 0.;
    double magcell = 0.;
    E0K[i] = getE0K(file, pv, magcell);
    pV[i] = pv;
    MagCell[i] = magcell;
  }

  if (MagCell[0] > 0.01) {
    _ismagnetic = true;
    _logger << apl::warning << " COMPOUND is magnetic, magcell :: " << MagCell[0] << apl::endl;
    _logger << apl::warning << " check your FVT data carefully " << apl::endl;
  }
}
// ***************************************************
//after excliuding negative phonon density of states resize vector sizes
void EOS::resize_EOSnPH() {
  vector<vector<vector<double> > > tmp_pDOS;
  vector<vector<vector<double> > > tmp_eDOS;
  vector<double> tmp_E0K;
  vector<double> tmp_pV;
  vector<double> tmp_fermi;
  vector<double> tmp_eosvol;
  tmp_pDOS = pDOS;
  tmp_eDOS = eDOS;
  tmp_E0K = E0K;
  tmp_pV = pV;
  tmp_fermi = fermi;
  tmp_eosvol = eosvol;
  E0K.clear();
  pV.clear();
  fermi.clear();
  eosvol.clear();
  uint size = 0;
  if (ZERO_STATIC_DIR_INDEX > 0)
    size = pDOS.size();
  else
    size = eDOS.size();
  pDOS.clear();
  eDOS.clear();

  uint j = 0;
  for (uint i = 0; i != size; i++) {
    if (i != NEGATIVE_FREQ[j]) {
      if ((ZERO_STATIC_DIR_INDEX > 0) && (i == (uint)ZERO_STATIC_DIR_INDEX))
        ZERO_STATIC_DIR_INDEX = pDOS.size();
      E0K.push_back(tmp_E0K[i]);
      pV.push_back(tmp_pV[i]);
      pDOS.push_back(tmp_pDOS[i]);
      eDOS.push_back(tmp_eDOS[i]);
      fermi.push_back(tmp_fermi[i]);
      eosvol.push_back(tmp_eosvol[i]);
    } else
      j++;
  }

  if (ZERO_STATIC_DIR_INDEX < 0) {
    pDOS.push_back(tmp_pDOS[tmp_pDOS.size() - 1]);
    ZERO_STATIC_DIR_INDEX = pDOS.size() - 1;
  }

  tmp_pDOS.clear();
  tmp_eDOS.clear();
  tmp_E0K.clear();
  tmp_pV.clear();
  tmp_fermi.clear();
  tmp_eosvol.clear();
}
// ***************************************************
//calculate zero point energies at each distorted volumes point
void EOS::getZeroPointVibrationEnergy() {
  zpe.clear();
  zpe.resize(pDOS.size(), 0.0);
  for (uint i = 0; i != pDOS.size(); i++) {
    double sum = 0.0;
    double steps = pDOS[i][1][0] - pDOS[i][0][0];
    for (uint j = 0; j != pDOS[i].size(); j++) {
      sum += pDOS[i][j][0] * pDOS[i][j][1];
    }
    sum *= steps * 0.5 * 4.1356651596736798;  //  0.5*h*nu
    zpe[i] = sum / 1000.0;
  }
}
// ***************************************************
//get electronic density of states and fermi energies and store to vectors
void EOS::geteDOS() {
  eDOS.clear();
  fermi.clear();
  for (uint i = 0; i != eosdir.size(); i++) {
    string filename = eosdir[i] + string("/") + string("DOSCAR.static");
    double Fermi = 0.0;
    eDOS.push_back(geteDOS<double>(Fermi, filename));
    fermi.push_back(Fermi);
  }
}
// ***************************************************
//read electronic density of states
template <class T>
vector<vector<T> > EOS::geteDOS(T &fermi, string file) {
  //CO - START
  if (!aurostd::EFileExist(file)) {
    _logger << apl::error << file << " does not exist" << apl::endl;
    exit(0);
  }
  vector<string> vlines;
  aurostd::efile2vectorstring(file, vlines);
  if (!vlines.size()) {
    _logger << apl::error << file << " does not exist" << apl::endl;
    exit(0);
  }
  //string command="";
  //string bz2=string(file)+string(".bz2");
  //bool t1=aurostd::FileExist(file);
  //bool t2=aurostd::FileExist(bz2);
  //if(!t1)
  //  {
  //if(t2){command=string("bzip2 -d ")+string(bz2);aurostd::execute2string(command);}
  //else{
  //  _logger << apl::error << file<<" does't exist" << apl::endl; exit(0);
  //}
  //    }
  vector<vector<T> > DOS;
  string line;
  //uint LC=0;
  uint line_count = 0;
  bool LC6 = false;
  //ifstream in (file.c_str());
  //if(!in.is_open()) {

  //while ( getline (in,line) ){
  while (line_count < vlines.size()) {
    line = vlines[line_count++];
    if (line == "") continue;
    if (line[0] == '#') continue;
    vector<T> vec = split<T>(line);
    if (line_count == 6) LC6 = true;
    if (line_count == 5) fermi = vec[3];
    if (LC6) { DOS.push_back(vec); }
    //LC++;
  }

  //in.clear();
  //in.close();
  //CO - END
  //if((t1+t2)!=2){command=string("bzip2 ")+string(file);aurostd::execute2string(command);}
  return DOS;
}
// ***************************************************
//read files with specific columns numbers
template <class T>
vector<vector<T> >
EOS::readfile(string file, const vector<uint> &readCOL) {
  //CO - START
  if (!aurostd::EFileExist(file)) {
    _logger << apl::error << file << " does not exist" << apl::endl;
    exit(0);
  }
  vector<string> vlines;
  aurostd::efile2vectorstring(file, vlines);
  if (!vlines.size()) {
    _logger << apl::error << file << " does not exist" << apl::endl;
    exit(0);
  }
  //string command="";
  //string bz2=string(file)+string(".bz2");
  //bool t1=aurostd::FileExist(file);
  //bool t2=aurostd::FileExist(bz2);
  //if(!t1){
  //  if(t2){command=string("bzip2 -d ")+string(bz2);aurostd::execute2string(command);}
  //  else{
  //_logger << apl::error << file<<" does't exist" << apl::endl; exit(0);
  //  }
  //}

  string line;
  vector<vector<T> > return_vec;
  return_vec.clear();
  uint line_count = 0;
  //ifstream in (file.c_str());
  //if (!in.is_open()){
  //  _logger << apl::error << file<<" does't exist" << apl::endl; exit(0);
  //}
  //while ( getline (in,line) ){
  while (line_count < vlines.size()) {
    line = vlines[line_count++];
    if (line == "") continue;
    if (line[0] == '#') continue;
    vector<T> vec = split<T>(line);
    vector<T> vec_readCOL;
    for (uint i = 0; i < readCOL.size(); i++) {
      for (uint j = 0; j < vec.size(); j++) {
        if (readCOL[i] == j) {
          vec_readCOL.push_back(vec[readCOL[i]]);
        }
      }
    }
    return_vec.push_back(vec_readCOL);
  }
  //in.clear();
  //in.close();
  //CO - END
  //if((t1+t2)!=2){command=string("bzip2 ")+string(file);aurostd::execute2string(command);}
  return return_vec;
}
// ***************************************************
//read 0K energies with magnetic properties
double EOS::getE0K(string file, double &pv, double &MAGCELL) {
  //CO - START
  if (!aurostd::EFileExist(file)) {
    _logger << apl::error << file << " does not exist" << apl::endl;
    exit(0);
  }
  vector<string> vlines;
  aurostd::efile2vectorstring(file, vlines);
  if (!vlines.size()) {
    _logger << apl::error << file << " does not exist" << apl::endl;
    exit(0);
  }
  pv = 0.0;
  MAGCELL = 0.0;
  //string command="";
  //string bz2=string(file)+string(".bz2");
  //bool t1=aurostd::FileExist(file);
  //bool t2=aurostd::FileExist(bz2);
  //if(!t1){
  //  if(t2){command=string("bzip2 -d ")+string(bz2);aurostd::execute2string(command);}
  //  else{
  //_logger << apl::error << file<<" does't exist" << apl::endl; exit(0);
  //  }
  //}
  //CO - END

  vector<string> WORD;
  vector<string> PV;
  vector<string> mag_cell;
  string line;
  uint line_count = 0;  //CO

  //CO - START
  //ifstream myfile (file.c_str());
  //if (!myfile.is_open()){
  //  _logger << apl::error << file<<" does't exist" << apl::endl; exit(0);
  //}
  //while ( getline (myfile,line) ){
  while (line_count < vlines.size()) {
    line = vlines[line_count++];
    if (line == "") continue;
    if ((line[0] == '#') || (line[0] == '/')) continue;
    if (line.find("H_cell") != std::string::npos)
      WORD.push_back(line);
    else if (line.find("PV_cell") != std::string::npos)
      PV.push_back(line);
    else if (line.find("mag_cell") != std::string::npos)
      mag_cell.push_back(line);
  }
  //myfile.clear();
  //myfile.close();
  //CO - END

  //depends on aflow.qmvasp.out format
  string s = WORD[WORD.size() - 1];
  vector<string> vec = split<string>(s);
  string str = vec[0];
  str.erase(str.begin(), str.begin() + 7);
  vector<double> vec1 = split<double>(str);

  string s1 = PV[PV.size() - 1];
  vector<string> vecPV = split<string>(s1);
  string str1 = vecPV[0];
  str1.erase(str1.begin(), str1.begin() + 8);
  vector<double> vecPV1 = split<double>(str1);
  pv = vecPV1[0];

  string s2 = mag_cell[mag_cell.size() - 1];
  vector<string> vecMag = split<string>(s2);
  string str2 = vecMag[0];
  str2.erase(str2.begin(), str2.begin() + 9);
  vector<double> tmpMag = split<double>(str2);
  MAGCELL = tmpMag[0];
  //if((t1+t2)!=2){command=string("bzip2 ")+string(file);aurostd::execute2string(command);}
  return vec1[0];
}
// ***************************************************
template <typename T>
std::vector<T> EOS::split(const std::string &line) {
  std::istringstream is(line);
  return std::vector<T>(std::istream_iterator<T>(is), std::istream_iterator<T>());
}
// ***************************************************
}  //apl end
