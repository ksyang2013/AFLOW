

#include "aflow_apl.h"

// Some parts are written within the C++0x support in GCC, expecially the std::thread,
// which is implemented in gcc 4.4 and higher.... For multithreads with std::thread see:
// http://www.justsoftwaresolutions.co.uk/threading/multithreading-in-c++0x-part-1-starting-threads.html
#if GCC_VERSION >= 40400  // added two zeros
#define AFLOW_APL_MULTITHREADS_ENABLE 1
#include <thread>
#else
#warning "The multithread parts of APL will be not included, since they need gcc 4.4 and higher (C++0x support)."
#endif

using namespace std;

namespace apl {

// ///////////////////////////////////////////////////////////////////////////

PhononDispersionCalculator::PhononDispersionCalculator(IPhononCalculator& pc, Logger& l) : _pc(pc), _logger(l) {
}

// ///////////////////////////////////////////////////////////////////////////

PhononDispersionCalculator::~PhononDispersionCalculator() {
  clear();
}

// ///////////////////////////////////////////////////////////////////////////

void PhononDispersionCalculator::clear() {
  _qpoints.clear();
  _freqs.clear();
}

//////////////////////////////////////////////////////////////////////////////

void PhononDispersionCalculator::initPathCoords(  //CO 180406
    const string& USER_DC_INITCOORDS,
    const string& USER_DC_INITLABELS,
    int USER_DC_NPOINTS, 
    bool CARTESIAN_COORDS) {
  if(USER_DC_INITCOORDS.empty() || USER_DC_INITLABELS.empty()){throw APLRuntimeError("apl::PhononDispersionCalculator::initPathCoords; Inputs are empty.");}
  _pb.defineCustomPoints(USER_DC_INITCOORDS,USER_DC_INITLABELS,_pc.getSupercell(),CARTESIAN_COORDS);
  _pb.setDensity(USER_DC_NPOINTS);
  _qpoints = _pb.getPath(); // Get points
}

void PhononDispersionCalculator::initPathLattice(const string& USER_DC_INITLATTICE,int USER_DC_NPOINTS){
  string lattice = USER_DC_INITLATTICE;
  if (lattice.empty()) {
    xstructure a(_pc.getInputCellStructure());
    //CO - START
    if (a.bravais_lattice_variation_type == "") {
      if (a.spacegroup == "") {
        if(a.space_group_ITC<1 || a.space_group_ITC>230){a.space_group_ITC = a.SpaceGroup_ITC();} //if (a.space_group_ITC == 0) { //CO 180214 - if not set then it could be 32767
        a.spacegroup = GetSpaceGroupName(a.space_group_ITC) + " #" + aurostd::utype2string(a.space_group_ITC);  //will break here if spacegroup is bad
      }
      // Use PLATON to get space group number if user did not get it...
      //a.platon2sg(_PLATON_P_EQUAL_DEFAULT,    //corey
      //	      _PLATON_P_EXACT_DEFAULT,
      //	      _PLATON_P_ANG_DEFAULT,
      //	      _PLATON_P_D1_DEFAULT,
      //	      _PLATON_P_D2_DEFAULT,
      //            _PLATON_P_D3_DEFAULT);
      //if( a.spacegroup.empty() )
      //throw apl::APLRuntimeError("apl::PhononDispersionCalculator::initPath(); The PLATON call to get spacegroup number failed. You have to specify it by DCINITSG in "+_AFLOWIN_);

      //vector<string> tokens;
      //aurostd::string2tokens(a.spacegroup,tokens,"#");
      //int spacegroupNumber = aurostd::string2utype<int>(tokens[1]);
      //tokens.clear();

      //lattice = LATTICE_Lattice_Variation_SpaceGroup(spacegroupNumber,_pc.getInputCellStructure());
      //lattice = LATTICE::SpaceGroup2LatticeVariation(spacegroupNumber,_pc.getInputCellStructure());
      lattice = LATTICE::SpaceGroup2LatticeVariation(a.space_group_ITC, a);
    }else{lattice = a.bravais_lattice_variation_type;}
    _logger << "The phonon dispersion curves will be generated for lattice variation " << lattice << "." << apl::endl;
  }
  //CO - END

  // cerr << "LATTICE=" << lattice << std::endl;
  // Suck point definition from the electronic structure part of AFLOW...
  _pb.takeAflowElectronicPath(lattice,_pc.getSupercell());             //CO 180406
                              //_pc.getInputCellStructure(),        //CO 180406
                              //_pc.getSuperCellStructure());       //CO 180406
  
  _pb.setDensity(USER_DC_NPOINTS);
  _qpoints = _pb.getPath(); // Get points
}

//////////////////////////////////////////////////////////////////////////////

void PhononDispersionCalculator::setPath(const string& USER_DC_OWNPATH) {
  // Get user's path...
  if (!USER_DC_OWNPATH.empty()) {
    if (USER_DC_OWNPATH.find('|') != string::npos)
      _qpoints = _pb.getPath(apl::PathBuilder::COUPLE_POINT_MODE, USER_DC_OWNPATH);
    else
      _qpoints = _pb.getPath(apl::PathBuilder::SINGLE_POINT_MODE, USER_DC_OWNPATH);
  }
}

//////////////////////////////////////////////////////////////////////////////

void PhononDispersionCalculator::calculateInOneThread(int startIndex, int endIndex) {
  //cout << "Thread: from " << startIndex << " to " <<  endIndex << std::endl;
  for (int iqp = startIndex; iqp < endIndex; iqp++) {
    _logger.updateProgressBar(iqp, _qpoints.size());
    _freqs[iqp] = _pc.getFrequency(_qpoints[iqp], _frequencyFormat);
    //std::this_thread::yield();
  }
}

//////////////////////////////////////////////////////////////////////////////

void PhononDispersionCalculator::calc(const IPCFreqFlags frequencyFormat) {
  // Save
  _frequencyFormat = frequencyFormat;

  // Maybe there was some error and the list of q-points is empty, hence bye-bye...
  if (_qpoints.empty())
    throw apl::APLRuntimeError("There are no points for calculation.");

// Compute frequencies for each q-point

#ifdef AFLOW_APL_MULTITHREADS_ENABLE

  // Get the number of CPUS
  int ncpus; //= sysconf(_SC_NPROCESSORS_ONLN);  // AFLOW_MachineNCPUs;  //CO 180214
  _pc.get_NCPUS(ncpus);  //CO 180214
  if (ncpus < 1) ncpus = 1;
  int qpointsPerCPU = _qpoints.size() / ncpus;

  // Show info
  if (ncpus == 1)
    _logger.initProgressBar("Calculating frequencies for PDIS");
  else
    _logger.initProgressBar("Calculating frequencies for PDIS (" + stringify(ncpus) + " threads)");

  // Prepare storage
  _freqs.clear();
  xvector<double> zero(_pc.getNumberOfBranches());
  for (uint i = 0; i < _qpoints.size(); i++)
    _freqs.push_back(zero);

  // Distribute the calculation
  int startIndex, endIndex;
  std::vector<std::thread*> threads;
  for (int icpu = 0; icpu < ncpus; icpu++) {
    startIndex = icpu * qpointsPerCPU;
    endIndex = startIndex + qpointsPerCPU;
    if (((uint)endIndex > _qpoints.size()) ||
        ((icpu == ncpus - 1) && ((uint)endIndex < _qpoints.size())))
      endIndex = _qpoints.size();
    threads.push_back(new std::thread(&PhononDispersionCalculator::calculateInOneThread, this, startIndex, endIndex));
  }

  // Wait to finish all threads here!
  for (uint i = 0; i < threads.size(); i++) {
    threads[i]->join();
    delete threads[i];
  }
  threads.clear();

  // Done
  _logger.finishProgressBar();

#else

  _logger.initProgressBar("Calculating frequencies for PDIS");
  for (uint iqp = 0; iqp < _qpoints.size(); iqp++) {
    _logger.updateProgressBar(iqp, _qpoints.size());
    _freqs.push_back(_pc.getFrequency(_qpoints[iqp], _frequencyFormat));
  }
  _logger.finishProgressBar();

#endif
}

//////////////////////////////////////////////////////////////////////////////

void PhononDispersionCalculator::writePDIS() {
  _logger << "Writing dispersion curves into file PDIS." << apl::endl;

  //CO - START
  //ofstream outfile("PDIS",ios_base::out);
  stringstream outfile;
  //if( !outfile.is_open() ) {
  //  throw apl::APLRuntimeError("Cannot open output PDIS file.");
  //}
  //CO - END

  // Write header
  outfile << "# Phonon dispersion curves calculated by Aflow" << std::endl;
  outfile << "#" << std::endl;
  outfile << "# <system>    \"" << _pc.getSuperCellStructure().title << "\"" << std::endl;
  outfile << "#" << std::endl;
  outfile << "# <units>     " << _frequencyFormat << std::endl;
  outfile << "# <nbranches> " << _pc.getNumberOfBranches() << std::endl;
  outfile << "# <npoints>   " << _freqs.size() << std::endl;
  outfile << "# <nsubpathp> " << _pb.getDensity() + 1 << std::endl;
  outfile << "#" << std::endl;

  // Write table of label points
  outfile << setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
  outfile << setprecision(8);
  double x = 0.0;
  double wholePathLength = _pb.getPathLength();
  map<double, string> labelMap;
  for (uint i = 1; i < _pb.getPointSize(); i++) {
    outfile << "# <label>     " << x << " "
            << setw(5) << _pb.getPointLabel(i)
            << std::endl;
    labelMap.insert(std::pair<double, string>(x, _pb.getPointLabel(i)));
    x += _pb.getPathLength(i) / wholePathLength;
  }
  outfile << "# <label>     " << 1.0 << " " << setw(5) << _pb.getPointLabel(_pb.getPointSize()) << std::endl;
  labelMap.insert(std::pair<double, string>(1.0, _pb.getPointLabel(_pb.getPointSize())));
  outfile << "#" << std::endl;

  // Write table of exact _qpoints + use label map to identify its labels
  x = 0.0;
  int subpath = 0;
  double xstep = 0.0;
  int p = 0;
  vector<double> exactPointPositions;
  for (uint i = 0; i < _qpoints.size(); i++) {
    // Check it
    if (isExactQPoint(_qpoints[i], _pc.getSuperCellStructure().lattice)) {
      // Is it new exact points
      uint j = 0;
      for (; j < exactPointPositions.size(); j++)
        if (exactPointPositions[j] == x) break;

      // If yes, add it....
      if (j == exactPointPositions.size()) {
        exactPointPositions.push_back(x);
        string name = "-";
        std::map<double, string>::iterator iter = labelMap.begin();
        for (; iter != labelMap.end(); iter++)
          if (fabs(iter->first - x) < _AFLOW_APL_EPS_) break;
        if (iter != labelMap.end())
          name = iter->second;
        outfile << "# <exact>     " << x << " "
                << setw(5) << name
                << std::endl;
      }
    }

    // Step of x will change in each subpath
    if (i % (_pb.getDensity() + 1) == 0) {
      if (i + 1 != _freqs.size())
        xstep = _pb.getPathLength(++subpath) / wholePathLength / (_pb.getDensity());
    }
    x += xstep;
    if (p != 0 && (p % _pb.getDensity()) == 0) {
      x -= xstep;
      p = -1;
    }
    p++;
  }
  outfile << "#" << std::endl;
  labelMap.clear();
  exactPointPositions.clear();

  // Write frequencies
  path_segment.clear();  //[PINKU]
  path.clear();          //[PINKU]
  outfile << setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
  outfile << setprecision(8);
  x = 0.0;
  subpath = 0;
  xstep = 0.0;
  p = 0;
  for (uint i = 0; i < _freqs.size(); i++) {
    outfile << setw(4) << p << " ";
    path_segment.push_back(p);  //[PINKU]
    outfile << setw(15) << x << " ";
    path.push_back(x);  //[PINKU]
    for (uint j = 1; j <= _pc.getNumberOfBranches(); j++)
      outfile << setw(15) << _freqs[i](j) << " ";
    outfile << std::endl;

    // Step of x will change in each subpath
    if (i % (_pb.getDensity() + 1) == 0) {
      if (i + 1 != _freqs.size())
        xstep = _pb.getPathLength(++subpath) / wholePathLength / (_pb.getDensity());
    }
    x += xstep;
    if (p != 0 && (p % _pb.getDensity()) == 0) {
      x -= xstep;
      p = -1;
    }
    p++;
  }

  //CO - START
  string filename = "PDIS";
  aurostd::stringstream2file(outfile, filename);
  if (!aurostd::FileExist(filename)) {
    throw apl::APLRuntimeError("Cannot open output PDIS file.");
  }
  //
  //outfile.clear();
  //outfile.close();
  //CO - END
}

//////////////////////////////////////////////////////////////////////////////

bool PhononDispersionCalculator::isExactQPoint(const xvector<double>& qpoint,
                                               const xmatrix<double>& lattice) {
  xcomplex<double> iONE(0.0, 1.0);

  bool isExact = false;
  for (_AFLOW_APL_REGISTER_ int i = 1; i <= 1; i++) {
    for (_AFLOW_APL_REGISTER_ int j = 1; j <= 1; j++) {
      for (_AFLOW_APL_REGISTER_ int k = 1; k <= 1; k++) {
        xvector<double> L = (((double)i) * lattice(1) +
                             ((double)j) * lattice(2) +
                             ((double)k) * lattice(3));
        xcomplex<double> p = exp(iONE * scalar_product(qpoint, L));
        if ((fabs(p.imag()) < _AFLOW_APL_EPS_) &&
            (fabs(p.real() - 1.0) < _AFLOW_APL_EPS_)) {
          isExact = true;
          break;
        }
      }
      if (isExact) break;
    }
    if (isExact) break;
  }

  return (isExact);
}

// ///////////////////////////////////////////////////////////////////////////

}  // namespace apl
