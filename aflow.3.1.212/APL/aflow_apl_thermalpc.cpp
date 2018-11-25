// [OBSOLETE] #include <iostream>

#include "aflow_apl.h"

using namespace std;

namespace apl {

// ///////////////////////////////////////////////////////////////////////////

ThermalPropertiesCalculator::ThermalPropertiesCalculator(IDOSCalculator& dosc, Logger& l)
    : _dosc(dosc), _logger(l) {
  // Copy to us, since it is more quick than call these function again and again
  _bins = _dosc.getBins();  //_bins == omega
  _dos = _dosc.getDOS();

  // Get step
  if (_bins.size() > 2)
    _stepDOS = _bins[1] - _bins[0];
  else
    throw APLRuntimeError("ThermalPropertiesCalculator::ThermalPropertiesCalculator(); Problem to obtain the step of DOS.");

  // Precompute it, since it wil be used many times
  _isCalcZeroPointVibrationEnergy_meV = false;
  _zeroPointVibrationEnergy_meV = getZeroPointVibrationEnergy(apl::meV);
  _isCalcZeroPointVibrationEnergy_meV = true;
}

ThermalPropertiesCalculator::~ThermalPropertiesCalculator() {
  this->clear();
}

// ///////////////////////////////////////////////////////////////////////////

void ThermalPropertiesCalculator::clear() {
  _bins.clear();
  _dos.clear();
  _zeroPointVibrationEnergy_meV = 0.0;
}

// ///////////////////////////////////////////////////////////////////////////

double ThermalPropertiesCalculator::getZeroPointVibrationEnergy(ThermalPropertiesUnits return_unit) {
  // Return precomputed value if the required units are the same
  if (_isCalcZeroPointVibrationEnergy_meV && (return_unit == apl::meV)) {
    return _zeroPointVibrationEnergy_meV;
  }

  // Compute
  double zpe = 0.0;
  for (uint i = 0; i < _bins.size(); i++) {
    zpe += _bins[i] * _dos[i];  //_bins == omega
  }
  zpe *= _stepDOS * 0.5 * 6.6260689633 / 1.60217733;  // 1/2 * hplanck/e-charge (e-charge converts hplanck to units of eV vs. J)
  return getScalingFactor(return_unit) * zpe;
}

// ///////////////////////////////////////////////////////////////////////////

double ThermalPropertiesCalculator::getInternalEnergy(double temperature_in_kelvins, ThermalPropertiesUnits return_unit) {
  if (temperature_in_kelvins < _AFLOW_APL_EPS_) return _zeroPointVibrationEnergy_meV;

  double u = 0.0;
  double beta = 1.0 / (0.0861734315 * temperature_in_kelvins);  // beta = 1/kBT = 1 / ([meV/K] [K])
  for (uint i = 0; i < _bins.size(); i++) {
    double hni = 4.1356673310 * _bins[i];  // hplanck in [eV.s] * 10^-15 * freq in [Hz] * 10^12 => hni in [meV].
    u += _dos[i] * hni / (exp(beta * hni) - 1.0);
  }
  return (_zeroPointVibrationEnergy_meV + getScalingFactor(return_unit) * u * _stepDOS);
}

// ///////////////////////////////////////////////////////////////////////////

double ThermalPropertiesCalculator::getVibrationalFreeEnergy(double temperature_in_kelvins, ThermalPropertiesUnits return_unit) {
  if (temperature_in_kelvins < _AFLOW_APL_EPS_) return _zeroPointVibrationEnergy_meV;

  double f = 0.0;
  double beta = 1.0 / (0.0861734315 * temperature_in_kelvins);  // beta = 1/kBT = 1 / ([meV/K] [K])
  for (uint i = 0; i < _bins.size(); i++) {
    double hni = 4.1356673310 * _bins[i];  // hplanck in [eV.s] * 10^-15 * freq in [Hz] * 10^12 => hni in [meV].
    //f += _dos[i] * ( 0.5 * hni + aurostd::ln( 1.0 - exp( -beta * hni ) ) / beta );
    f += _dos[i] * aurostd::ln(1.0 - exp(-beta * hni)) / beta;
  }
  return (_zeroPointVibrationEnergy_meV + (getScalingFactor(return_unit) * f * _stepDOS));
}

// ///////////////////////////////////////////////////////////////////////////

double ThermalPropertiesCalculator::getVibrationalEntropy(double temperature_in_kelvins, ThermalPropertiesUnits return_unit) {
  if (temperature_in_kelvins < _AFLOW_APL_EPS_) return 0.0;

  double u = getInternalEnergy(temperature_in_kelvins, apl::meV);
  double f = getVibrationalFreeEnergy(temperature_in_kelvins, apl::meV);

  return (getScalingFactor(return_unit) * (u - f) / temperature_in_kelvins);
}

// ///////////////////////////////////////////////////////////////////////////

double ThermalPropertiesCalculator::getIsochoricSpecificHeat(double temperature_in_kelvins, ThermalPropertiesUnits return_unit) {
  double cv = 0.0;
  double beta = 1.0 / (0.0861734315 * temperature_in_kelvins);  // beta = 1/kBT = 1 / ([meV/K] [K])
  for (uint i = 0; i < _bins.size(); i++) {
    double bhni = beta * 4.1356673310 * _bins[i];   // hplanck in [eV.s] * 10^-15 * freq in [Hz] * 10^12 => hni in [meV].
    double ebhni = exp(bhni);
    cv += _dos[i] * 0.0861734315 * bhni * bhni / ((1.0 - 1.0 / ebhni) * (ebhni - 1.0)); // beta = 1/kBT = 1 / ([meV/K] [K])
  }
  if (isnan(cv)) return 0.0;
  return getScalingFactor(return_unit) * cv * _stepDOS;
}

// ///////////////////////////////////////////////////////////////////////////

double ThermalPropertiesCalculator::getScalingFactor(ThermalPropertiesUnits units) {
  switch (units) {
    case eV:
    case eVK:
      return 0.001;
    case meV:
    case meVK:
      return 1.0;
    case ueV:
    case ueVK:
      return 1000.0;
    case kB:
      return 1.0 / 0.0861734315;
  }

  return 1.0;
}

// ///////////////////////////////////////////////////////////////////////////

void ThermalPropertiesCalculator::writeTHERMO(double USER_TP_TSTART, double USER_TP_TEND, double USER_TP_TSTEP) {
  // The output file THERMO
  //CO - START
  //ofstream outfile("THERMO",ios_base::out);
  stringstream outfile;
  //if( !outfile.is_open() )
  //{
  //    throw apl::APLRuntimeError("ThermalPropertiesCalculator::writeTHERMO(); Cannot open output THERMO file.");
  //}
  //CO - END

  _logger << "Writing thermodynamic properties into file " << DEFAULT_APL_THERMO_FILE << "." << apl::endl;

  outfile << "#  T(K)     U0(meV/cell)    U(meV/cell)     F(meV/cell)      S(kB/cell)      Cv(kB/cell)" << std::endl;
  outfile << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
  int n = (int)((USER_TP_TEND - USER_TP_TSTART) / USER_TP_TSTEP);
  for (int i = 0; i <= n; i++) {
    double TEMPERATURE_IN_K = USER_TP_TSTART + i * USER_TP_TSTEP;
    double U0 = getZeroPointVibrationEnergy(apl::meV);
    double U = getInternalEnergy(TEMPERATURE_IN_K, apl::meV);
    double Fvib = getVibrationalFreeEnergy(TEMPERATURE_IN_K, apl::meV);
    double Svib = getVibrationalEntropy(TEMPERATURE_IN_K, apl::kB);
    double Cv = getIsochoricSpecificHeat(TEMPERATURE_IN_K, apl::kB);

    outfile << setw(8) << setprecision(2) << TEMPERATURE_IN_K << setprecision(8)
            << setw(15) << U0 << "\t"
            << setw(15) << U << "\t"
            << setw(15) << Fvib << "\t"
            << setw(15) << Svib << "\t"
            << setw(15) << Cv << std::endl;
  }
  //CO - START
  aurostd::stringstream2file(outfile, DEFAULT_APL_THERMO_FILE);
  if (!aurostd::FileExist(DEFAULT_APL_THERMO_FILE)) {
    string function = "ThermalPropertiesCalculator::writeTHERMO()";
    string message = "Cannot open output file " + DEFAULT_APL_THERMO_FILE + ".";
    throw aurostd::xerror(function, message, _FILE_ERROR_);
//    throw apl::APLRuntimeError("ThermalPropertiesCalculator::writeTHERMO(); Cannot open output THERMO file.");
  }
  //outfile.clear();
  //outfile.close();
  //CO - START
}

// ///////////////////////////////////////////////////////////////////////////

}  // namespace APL
