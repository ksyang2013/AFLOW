#include "aflow_apl.h"

using namespace std;

namespace apl {

  // ///////////////////////////////////////////////////////////////////////////

  LinearTetrahedronMethod::LinearTetrahedronMethod(IPhononCalculator& pc, IReciprocalPointGrid& rg, Logger& l)
    : DOSCalculator(pc,rg,l) {
    // OK. Everybody knows we need at least 4 points for 1 tetrahedra, if it
    // is sufficient, it is an another question...
    if( _qpoints.size() < 4 ) {
      APLRuntimeError("LinearTetrahedronMethod: The number of qpoints is not sufficient.");
    }
    generateTetrahedras();
  }

  // ///////////////////////////////////////////////////////////////////////////

  LinearTetrahedronMethod::~LinearTetrahedronMethod() {
    clear();
  }

  // ///////////////////////////////////////////////////////////////////////////

  void LinearTetrahedronMethod::clear() {
    _irrTetrahedraList.clear();
    _irrTetrahedraWeightList.clear();
  }

  // ///////////////////////////////////////////////////////////////////////////

  void LinearTetrahedronMethod::generateTetrahedras() {
    // Define 6 tetrahedras
    vector< vector< xvector<double> > > tetrahedra0, tetrahedra;
    vector< xvector<double>  > th;
    xvector<double> p(3);
    // Tetrahedra 1348 or 2437 or 3126 or 4215
    p(1) = 0; p(2) = 0; p(3) = 0; th.push_back(p);
    p(1) = 0; p(2) = 1; p(3) = 0; th.push_back(p);
    p(1) = 1; p(2) = 1; p(3) = 0; th.push_back(p);
    p(1) = 1; p(2) = 1; p(3) = 1; th.push_back(p);
    tetrahedra0.push_back(th); th.clear();
    // Tetrahedra 1248 or 2137 or 3426 or 4315
    p(1) = 0; p(2) = 0; p(3) = 0; th.push_back(p);
    p(1) = 1; p(2) = 0; p(3) = 0; th.push_back(p);
    p(1) = 1; p(2) = 1; p(3) = 0; th.push_back(p);
    p(1) = 1; p(2) = 1; p(3) = 1; th.push_back(p);
    tetrahedra0.push_back(th); th.clear();
    // Tetrahedra 1268 or 2157 or 3486 or 4375
    p(1) = 0; p(2) = 0; p(3) = 0; th.push_back(p);
    p(1) = 1; p(2) = 0; p(3) = 0; th.push_back(p);
    p(1) = 1; p(2) = 0; p(3) = 1; th.push_back(p);
    p(1) = 1; p(2) = 1; p(3) = 1; th.push_back(p);
    tetrahedra0.push_back(th); th.clear();
    // Tetrahedra 1378 or 2487 or 3156 or 4265
    p(1) = 0; p(2) = 0; p(3) = 0; th.push_back(p);
    p(1) = 0; p(2) = 1; p(3) = 0; th.push_back(p);
    p(1) = 0; p(2) = 1; p(3) = 1; th.push_back(p);
    p(1) = 1; p(2) = 1; p(3) = 1; th.push_back(p);
    tetrahedra0.push_back(th); th.clear();
    // Tetrahedra 1578 or 2687 or 3756 or 4865
    p(1) = 0; p(2) = 0; p(3) = 0; th.push_back(p);
    p(1) = 0; p(2) = 0; p(3) = 1; th.push_back(p);
    p(1) = 0; p(2) = 1; p(3) = 1; th.push_back(p);
    p(1) = 1; p(2) = 1; p(3) = 1; th.push_back(p);
    tetrahedra0.push_back(th); th.clear();
    // Tetrahedra 1568 or 2657 or 3786 or 4875
    p(1) = 0; p(2) = 0; p(3) = 0; th.push_back(p);
    p(1) = 0; p(2) = 0; p(3) = 1; th.push_back(p);
    p(1) = 1; p(2) = 0; p(3) = 1; th.push_back(p);
    p(1) = 1; p(2) = 1; p(3) = 1; th.push_back(p);
    tetrahedra0.push_back(th); th.clear();

    // Now each corner of our defined tetrahedrons will be shifted by x, y, and xy and
    // the resulted tetrahedron will be measured. That shift, which will result in the
    // shortest tetrahedron will be used in next calculation. The prefer of choice of
    // shortestegdes for all tetrahedra is due to avoding the long interpolation
    // distances.
    int lxx = 0;
    int lyy = 0;
    double gmax = 1E30;
    double gmin = 0.0;

    for(int lx = 0; lx <= 1; lx++)
      for(int ly = 0; ly <= 1; ly++) {
        for(uint i = 0; i < tetrahedra.size(); i++) tetrahedra[i].clear();
        tetrahedra.clear();

        for(int i = 0; i < 6; i++) {
	  th.clear();
	  for(int j = 0; j < 4; j++) {
	    p = tetrahedra0[i][j];
	    if( lx == 1 ) p(1) = 1 - p(1);
	    if( ly == 1 ) p(2) = 1 - p(2);
	    th.push_back(p);
	  }
	  tetrahedra.push_back(th);
        }

        double lmax = 0.0;
        double lmin = 1E30;
        for(int i = 0; i < 6; i++) {
	  for(int j = 0; j < 4; j++) {
	    tetrahedra[i][j] = F2C(_rg.getReciprocalLattice(),tetrahedra[i][j]);
	  }

	  for(int j = 0; j < 4-1; j++) {
	    for(int k = j+1; k < 4; k++) {
	      double d = aurostd::modulus( tetrahedra[i][j] - tetrahedra[i][k] );
	      if( d > lmax ) lmax = d;
	      if( d < lmin ) lmin = d;
	    }
	  }
        }

        if( lmax < gmax ) {
	  lxx = lx;
	  lyy = ly;
	  gmax = lmax;
	  gmin = lmin;
        }

      }
    if(gmin) {;} // dummy load

    // OK. We know now which shift will produce the most compact tetrahedra, hence
    // setup this configuration
    for(uint i = 0; i < tetrahedra.size(); i++) tetrahedra[i].clear();
    tetrahedra.clear();
    for(int i = 0; i < 6; i++) {
      th.clear();
      for(int j = 0; j < 4; j++) {
	p = tetrahedra0[i][j];
	if( lxx == 1 ) p(1) = 1 - p(1);
	if( lyy == 1 ) p(2) = 1 - p(2);
	th.push_back(p);
      }
      tetrahedra.push_back(th);
    }
    th.clear();
    for(uint i = 0; i < tetrahedra0.size(); i++) tetrahedra0[i].clear();
    tetrahedra0.clear();

    // We need this map because we are loking for irreducible tetrahedras
    aurostd::xtensor3<int> allToIrrPointMap = _rg.getAllToIrrPointMap();

    // OK. Go throught the whole k-point mesh and generate all tetrahedra...
    _irrTetrahedraList.clear();
    _irrTetrahedraWeightList.clear();
    for(int i3 = 1; i3 <= _rg.getN(3); i3++)
      for(int i2 = 1; i2 <= _rg.getN(2); i2++)
	for(int i1 = 1; i1 <= _rg.getN(1); i1++) {

	  // Get the mapping between the microcell's 8 corners and the ID of kpoint
	  // in the irreducible list of points
	  xtensor3<int> tetraCornerToIrrPoint(2,2,2,0,0,0);        
	  for(int k3 = 0; k3 <= 1; k3++) {
            _AFLOW_APL_REGISTER_ int j3 = ( ( i3 + k3 - 1 ) % _rg.getN(3) ) + 1;
            for(int k2 = 0; k2 <= 1; k2++) {
	      _AFLOW_APL_REGISTER_ int j2 = ( ( i2 + k2 - 1 ) % _rg.getN(2) ) + 1;
	      for(int k1 = 0; k1 <= 1; k1++) {
		_AFLOW_APL_REGISTER_ int j1 = ( ( i1 + k1 - 1 ) % _rg.getN(1) ) + 1;
		tetraCornerToIrrPoint(k1,k2,k3) = allToIrrPointMap(j1,j2,j3);
		//cout << k1 << k2 << k3 << " -> " << j1 << " " << j2 << " " << j3 << " -> " << tetraCornerToIrrPoint(k1,k2,k3) << std::endl;
	      }
	    }
	  }

	  // Get all 6 tetrahedras now...
	  vector<int> tetraIrrPointID;
	  for(int i = 0; i < 6; i++) {
            tetraIrrPointID.clear();
            for(int j = 0; j < 4; j++) {
	      tetraIrrPointID.push_back( tetraCornerToIrrPoint(
							       (int)tetrahedra[i][j](1), (int)tetrahedra[i][j](2),
							       (int)tetrahedra[i][j](3) ) );
	    }

            // Order
            for(int j = 0; j < 4-1; j++) {
	      for(int k = j+1; k < 4; k++) {
		if( tetraIrrPointID[k] < tetraIrrPointID[j] ) {
		  int temp = tetraIrrPointID[k];
		  tetraIrrPointID[k] = tetraIrrPointID[j];
		  tetraIrrPointID[j] = temp;
		}
	      }
	    }
            //for(int j = 0; j < 4; j++) cout << setw(5) << tetraIrrPointID[j] << " ";

            // Did we generate it already?
            uint l = 0;
            for( ; l < _irrTetrahedraList.size(); l++) {                
	      uint m = 0;
	      for( ; m < 4; m++)
		if( _irrTetrahedraList[l][m] != tetraIrrPointID[m] ) break;
	      if( m == 4 ) break;
	    }

            if( l == _irrTetrahedraList.size() ) {
	      _irrTetrahedraList.push_back(tetraIrrPointID);
	      _irrTetrahedraWeightList.push_back(1);
	      //cout << "new";
	    } else {
	      _irrTetrahedraWeightList[l]++;
	      //cout << irrTetrahedraWeightList[l] << "th in the list";
	    }
            //cout << std::endl;

	  }
	  tetraIrrPointID.clear();
	}

    _weightVolumeOfEachTetrahedron = 1.0 / (double)( 6 * _rg.getN(1) * _rg.getN(2) * _rg.getN(3) );
    _logger << "Generated " << (int)_irrTetrahedraList.size() << " irreducible tetrahedras." << apl::endl;

    // Clear stuff
    allToIrrPointMap.clear();
    for(uint i = 0; i < tetrahedra.size(); i++) tetrahedra[i].clear();
    tetrahedra.clear();
  }

  // ///////////////////////////////////////////////////////////////////////////

  void LinearTetrahedronMethod::rawCalc(int USER_DOS_NPOINTS) {
    if(USER_DOS_NPOINTS) {;} // dummy load
    // The loop over the all tetrahedras
    for(uint itet = 0; itet < _irrTetrahedraList.size(); itet++) {
      // We are going over the each phonon branch
      for(int ibranch = _freqs[0].lrows; ibranch <= _freqs[0].urows; ibranch++) {
	// Get frequency value at all 4 corners
	vector<double> f(4);
	for(uint icorner = 0; icorner < 4; icorner++)
	  f[icorner] = _freqs[_irrTetrahedraList[itet][icorner]](ibranch);

	// Order frequencies
	for(uint i = 0; i < 3; i++) {
	  for(uint j = i+1; j < 4; j++) {
	    if( f[j] < f[i] ) {
	      double temp = f[j];
	      f[j] = f[i];
	      f[i] = temp;    
	    }
	  }
	}
	double &fmin = f[0];
	double &fmax = f[3];

	// If the max frequency is higher that _maxFreq there is no contribution to DOS
	if( fmax > ( _bins.back() + _halfStepDOS ) ) continue;

	// If the min frequencies is lower than _minFreq there is no contribution to DOS
	if( fmin < ( _bins.front() - _halfStepDOS ) ) continue;

	// OK. We are in, so calculate the contributions to DOS
	// Calculate the indices of bins which will be affected by the currect frequencies
	// at current tetrahedra
	int kstart = (int)( ( fmin - _bins.front() - _halfStepDOS ) / _stepDOS );
	if( kstart < 0 ) kstart = 0;
	if( kstart > (int)_bins.size() ) kstart = _bins.size() - 1;
	int kstop = (int)( ( fmax - _bins.front() - _halfStepDOS ) / _stepDOS );
	if( kstop < 0 ) kstop = 0;
	if( kstop > (int)_bins.size() ) kstop = _bins.size() - 1;

	// Precompute some constants participating in the DOS constribution formulas
	double f21 = f[1] - f[0];
	double f31 = f[2] - f[0];
	double f32 = f[2] - f[1];
	double f41 = f[3] - f[0];
	double f42 = f[3] - f[1];
	double f43 = f[3] - f[2];
	double VTG = _irrTetrahedraWeightList[itet] * _weightVolumeOfEachTetrahedron;
	double cc12 = 3.0 * VTG / ( f21 * f31 * f41 );
	double cc34 = 3.0 * VTG / ( f41 * f42 * f43 );
	double cc23 = VTG / ( f31 * f41 );
	double cc23a = 3.0 * f21 * cc23;
	double cc23b = 6.0 * cc23;
	double cc23c = -3.0 * cc23 * ( f31 + f42 ) / ( f32 * f42 );

	// Calculate DOS constributions from this tetrahedra
	for(int k = kstart; k <= kstop; k++) {
	  double fbin = _minFreq + k*_stepDOS + _halfStepDOS;
	  if( fbin < f[0] ) {
	    _dos[k] += 0.0;
	  } else if( ( f[0] <= fbin ) && ( fbin <= f[1] ) ) {
	    double x = fbin - f[0];
	    _dos[k] += cc12 * x * x;
	  } else if( ( f[1] < fbin ) && ( fbin <= f[2] ) ) {
	    double x = fbin - f[1];
	    _dos[k] += cc23a + cc23b * x + cc23c * x * x;
	  } else if( ( f[2] < fbin ) && ( fbin <= f[3] ) ) {
	    double x = f[3] - fbin;
	    _dos[k] += cc34 * x * x;
	  } else if( fbin > f[3] ) {
	    _dos[k] += 0.0;
	  }
	}

	// Clear
	f.clear();
      }
    }
  }

  // ///////////////////////////////////////////////////////////////////////////

} // namespace apl
