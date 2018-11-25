// [OBSOLETE] #include <iostream>
// [OBSOLETE] #include <sstream>
// [OBSOLETE] #include <string>
// [OBSOLETE] #include <limits>

#include "aflow_apl.h"

//CO - START
// Some parts are written within the C++0x support in GCC, especially the std::thread,
// which is implemented in gcc 4.4 and higher.... For multithreads with std::thread see:
// http://www.justsoftwaresolutions.co.uk/threading/multithreading-in-c++0x-part-1-starting-threads.html
#if GCC_VERSION >= 40400  // added two zeros
#define AFLOW_APL_MULTITHREADS_ENABLE 1
#include <thread>
#else
#warning "The multithread parts of APL will be not included, since they need gcc 4.4 and higher (C++0x support)."
#endif
//CO - END

#define MIN_FREQ_TRESHOLD -0.1

using namespace std;

namespace apl {

// ///////////////////////////////////////////////////////////////////////////

DOSCalculator::DOSCalculator(IPhononCalculator& pc, IReciprocalPointGrid& rg, Logger& l)
    : _pc(pc), _rg(rg), _logger(l) {
  clear();
  calculateFrequencies();
}

// ///////////////////////////////////////////////////////////////////////////

DOSCalculator::~DOSCalculator() {
  clear();
}

// ///////////////////////////////////////////////////////////////////////////

void DOSCalculator::clear() {
  _qpoints.clear();
  _qweights.clear();
  _freqs.clear();
  _bins.clear();
  _dos.clear();
}

// ///////////////////////////////////////////////////////////////////////////

//CO - START
void DOSCalculator::calculateInOneThread(int startIndex, int endIndex) {
  //cout << "Thread: from " << startIndex << " to " <<  endIndex << std::endl;
  for (int iqp = startIndex; iqp < endIndex; iqp++) {
    _logger.updateProgressBar(iqp, _qpoints.size());
    _freqs[iqp] = _pc.getFrequency(_qpoints[iqp], apl::THZ | apl::ALLOW_NEGATIVE);
    //std::this_thread::yield();
  }
}
//CO - END

//////////////////////////////////////////////////////////////////////////////

void DOSCalculator::calculateFrequencies() {
  // Get q-points for which to calculate the frequencies
  _qpoints = _rg.getPoints();
  _qweights = _rg.getWeights();

// 07-08-2010 We do not need it anymore, since mpmesh has all points in catesian coords now
// Transform points to the form as it is expected by CPC:
// Transform from the reciprocal lattice of the primitive cell to cartesian coordinates
// for(uint t = 0; t < _qpoints.size(); t++)
//    _qpoints[t] = trasp(ReciprocalLattice(_pc.getPrimitiveCellStructure().lattice)) * _qpoints[t];

//CO - START
#ifdef AFLOW_APL_MULTITHREADS_ENABLE

  // Get the number of CPUS
  int ncpus; //= sysconf(_SC_NPROCESSORS_ONLN);  // AFLOW_MachineNCPUs;  //CO 180214
  _pc.get_NCPUS(ncpus);  //CO 180214
  if (ncpus < 1) ncpus = 1;
//  int qpointsPerCPU = _qpoints.size() / ncpus;  OBSOLETE ME 180801

  // Show info
  if (ncpus == 1)
    _logger.initProgressBar("Calculating frequencies for DOS");
  else
    _logger.initProgressBar("Calculating frequencies for DOS (" + stringify(ncpus) + " threads)");

  // Prepare storage
  _freqs.clear();
  xvector<double> zero(_pc.getNumberOfBranches());
  for (uint i = 0; i < _qpoints.size(); i++)
    _freqs.push_back(zero);

  // Distribute the calculation
  int startIndex, endIndex;
  std::vector<std::thread*> threads;
  vector<vector<int> > thread_dist = getThreadDistribution((int) _qpoints.size(), ncpus);
  for (int icpu = 0; icpu < ncpus; icpu++) {
    startIndex = thread_dist[icpu][0];
    endIndex = thread_dist[icpu][1];
    threads.push_back(new std::thread(&DOSCalculator::calculateInOneThread, this, startIndex, endIndex));
  }

/* OBSOLETE ME 180801
  for (int icpu = 0; icpu < ncpus; icpu++) {
    startIndex = icpu * qpointsPerCPU;
    endIndex = startIndex + qpointsPerCPU;
    if (((uint)endIndex > _qpoints.size()) ||
        ((icpu == ncpus - 1) && ((uint)endIndex < _qpoints.size())))
      endIndex = _qpoints.size();
    threads.push_back(new std::thread(&DOSCalculator::calculateInOneThread, this, startIndex, endIndex));
  }
*/

  // Wait to finish all threads here!
  for (uint i = 0; i < threads.size(); i++) {
    threads[i]->join();
    delete threads[i];
  }
  threads.clear();

  // Done
  _logger.finishProgressBar();

#else

  // Calculate frequencies
  _logger.initProgressBar("Calculating frequencies for DOS");
  for (uint iqp = 0; iqp < _qpoints.size(); iqp++) {
    _logger.updateProgressBar(iqp, _qpoints.size());
    _freqs.push_back(_pc.getFrequency(_qpoints[iqp], apl::THZ | apl::ALLOW_NEGATIVE));
  }
  _logger.finishProgressBar();

#endif
  //CO - END

  //if freq > MIN_FREQ_TRESHOLD considerd as +ve freq [PINKU]
  for (uint i = 0; i < _freqs.size(); i++) {
    for (int j = _freqs[i].lrows; j <= _freqs[i].urows; j++) {
      if ((_freqs[i][j] < 0.00) && (_freqs[i][j] > MIN_FREQ_TRESHOLD)) _freqs[i][j] = 0.00;
    }
  }
  //PINKU END

  // Get min and max values
  _maxFreq = -1.0;
  _minFreq = 0.0;
  for (uint i = 0; i < _freqs.size(); i++) {
    for (int j = _freqs[i].lrows; j <= _freqs[i].urows; j++) {
      if (_freqs[i](j) > _maxFreq) _maxFreq = _freqs[i](j);
      if (_freqs[i](j) < _minFreq) _minFreq = _freqs[i](j);
    }
  }
  _maxFreq += 1.0;
  if (_minFreq < MIN_FREQ_TRESHOLD) _minFreq -= 1.0;
}

// ///////////////////////////////////////////////////////////////////////////

void DOSCalculator::smearWithGaussian(vector<double>& dos, double h, double sigma) {
  // Construct table for gaussian function
  int ng = (int)(6.0 * sigma / h + 1.0);
  double fact = 1.0 / (sqrt(2.0 * M_PI) * sigma);
  vector<double> gauss;
  double gnorm = 0.0;
  for (int ig = -ng; ig <= ng; ig++) {
    double eg = ig * h;
    double arg = eg * eg / (sigma * sigma) / 2.0;
    gauss.push_back(fact * exp(-arg));
    gnorm += gauss.back();
  }

  // Norm gauss table to one
  gnorm *= h;
  for (int ig = -ng; ig <= ng; ig++) {
    gauss[ig + ng] /= gnorm;
  }

  // Prepare new dos
  vector<double> newdos;
  for (uint i = 0; i < dos.size(); i++) {
    newdos.push_back(0.0);
  }

  // Convolute...
  for (int ie = 0; ie < (int)dos.size(); ie++) {
    double wt = dos[ie] * h;
    for (int jg = -ng; jg <= ng; jg++) {
      int je = ie + jg;
      if (je < 0) continue;
      if (je >= (int)dos.size()) continue;
      newdos[je] += gauss[jg + ng] * wt;
    }
  }

  //
  dos.clear();
  for (uint i = 0; i < newdos.size(); i++) {
    dos.push_back(newdos[i]);
  }
  newdos.clear();
  gauss.clear();
}

// ///////////////////////////////////////////////////////////////////////////

void DOSCalculator::calc(int USER_DOS_NPOINTS) {
  calc(USER_DOS_NPOINTS, 0.0);
}

// ///////////////////////////////////////////////////////////////////////////

void DOSCalculator::calc(int USER_DOS_NPOINTS, double USER_DOS_SMEAR) {
  // Calculate steps
  _stepDOS = (_maxFreq - _minFreq) / (double)USER_DOS_NPOINTS;
  _halfStepDOS = 0.5 * _stepDOS;

  // Clear old stuff
  _dos.clear();
  _bins.clear();

  // Prepare storagearrays
  for (int k = 0; k < USER_DOS_NPOINTS; k++) {
    _dos.push_back(0);
    _bins.push_back(_minFreq + k * _stepDOS + _halfStepDOS);
  }

  // Perform the raw specific calculation by method
  rawCalc(USER_DOS_NPOINTS);

  // Smooth DOS by gaussians
  if (USER_DOS_SMEAR > 1E-6)
    smearWithGaussian(_dos, _stepDOS, USER_DOS_SMEAR);

  // Normalize to number of branches
  double sum = 0.0;
  for (int k = 0; k < USER_DOS_NPOINTS; k++)
    sum += _dos[k];
  sum /= _pc.getNumberOfBranches();

  for (int k = 0; k < USER_DOS_NPOINTS; k++)
    _dos[k] /= (sum * _stepDOS);
}

// ///////////////////////////////////////////////////////////////////////////

void DOSCalculator::writePDOS() {
  // Write PHDOS file
  //CO - START
  //ofstream outfile("PDOS",ios_base::out);
  stringstream outfile;
  //if( !outfile.is_open() )
  //{
  //    throw apl::APLRuntimeError("DOSCalculator::writePDOS(); Cannot open output PDOS file.");
  //}
  //CO - END

  double factorTHz2Raw = _pc.getFrequencyConversionFactor(apl::THZ, apl::RAW);
  double factorRaw2rcm = _pc.getFrequencyConversionFactor(apl::RAW, apl::RECIPROCAL_CM);
  double factorRaw2meV = _pc.getFrequencyConversionFactor(apl::RAW, apl::MEV);

  _logger << "Writing phonon density of states into file " << DEFAULT_APL_PDOS_FILE << "." << apl::endl;
  //outfile << "############### ############### ############### ###############" << std::endl;
  outfile << "#    f(THz)      1/lambda(cm-1)      E(meV)          pDOS      " << std::endl;
  outfile << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
  outfile << setprecision(8);
  for (uint k = 0; k < _dos.size(); k++) {
    outfile << setw(15) << _bins[k] << " "
            << setw(15) << _bins[k] * factorTHz2Raw * factorRaw2rcm << " "
            << setw(15) << _bins[k] * factorTHz2Raw * factorRaw2meV << " "
            << setw(15) << _dos[k] << std::endl;
  }

  //CO - START
  aurostd::stringstream2file(outfile, DEFAULT_APL_PDOS_FILE);
  if (!aurostd::FileExist(DEFAULT_APL_PDOS_FILE)) {
    string function = "DOSCalculator::writePDOS()";
    string message = "Cannot open output file " + DEFAULT_APL_PDOS_FILE + ".";
    throw aurostd::xerror(function, message, _FILE_ERROR_);
//    throw apl::APLRuntimeError("DOSCalculator::writePDOS(); Cannot open output PDOS file.");
  }
  //outfile.clear();
  //outfile.close();
  //CO - END
}

// ///////////////////////////////////////////////////////////////////////////

vector<double> DOSCalculator::getBins() {
  return _bins;
}

vector<double> DOSCalculator::getDOS() {
  return _dos;
}

bool DOSCalculator::hasNegativeFrequencies() {
  return (_minFreq < MIN_FREQ_TRESHOLD ? true : false);
}

// ///////////////////////////////////////////////////////////////////////////
//PINKU - START
void DOSCalculator::writePDOS(string path, string ex)  //[PINKU]
{
  //CO - START
  // Write PHDOS file
  //ofstream outfile(file.c_str(),ios_base::out);
  stringstream outfile;
  //if( !outfile.is_open() )
  //{
  //    throw apl::APLRuntimeError("DOSCalculator::writePDOS(); Cannot open output PDOS file.");
  //}
  //CO - END

  double factorTHz2Raw = _pc.getFrequencyConversionFactor(apl::THZ, apl::RAW);
  double factorRaw2rcm = _pc.getFrequencyConversionFactor(apl::RAW, apl::RECIPROCAL_CM);
  double factorRaw2meV = _pc.getFrequencyConversionFactor(apl::RAW, apl::MEV);

  _logger << "Writing phonon density of states into file " << DEFAULT_APL_PDOS_FILE << "." << apl::endl;
  //outfile << "############### ############### ############### ###############" << std::endl;
  outfile << "#    f(THz)      1/lambda(cm-1)      E(meV)          pDOS      " << std::endl;
  outfile << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
  outfile << setprecision(8);
  for (uint k = 0; k < _dos.size(); k++) {
    outfile << setw(15) << _bins[k] << " "
            << setw(15) << _bins[k] * factorTHz2Raw * factorRaw2rcm << " "
            << setw(15) << _bins[k] * factorTHz2Raw * factorRaw2meV << " "
            << setw(15) << _dos[k] << std::endl;
  }

  //CO - START
  string file = path + "/" + DEFAULT_APL_PDOS_FILE + "." + ex;
  aurostd::stringstream2file(outfile, file);
  if (!aurostd::FileExist(file)) {
    string function = "DOSCalculator::writePDOS()";
    string message = "Cannot open output file " + DEFAULT_APL_PDOS_FILE + ".";
    throw aurostd::xerror(function, message, _FILE_ERROR_);
//    throw apl::APLRuntimeError("DOSCalculator::writePDOS(); Cannot open output PDOS file.");
  }
  //outfile.clear();
  //outfile.close();
  //CO - END
}
//PINKU - END
// ///////////////////////////////////////////////////////////////////////////

}  // namespace apl
