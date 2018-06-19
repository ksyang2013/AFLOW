#include "aflow_apl.h"

using namespace std;

namespace apl {

// ///////////////////////////////////////////////////////////////////////////

MonkhorstPackMesh::MonkhorstPackMesh(int na, int nb, int nc, const xstructure& xs, Logger& l) : _logger(l) {
    // Get out copy of the structure
    xstructure xstr(xs);

    // Setup both lattices
    xstr.FixLattices();
    _rlattice = xstr.lattice;
    _klattice = ReciprocalLattice(_rlattice);

    // Set correct dimension of arrays
    setDensity(na,nb,nc);

    // Generate the whole grid of kpoints
    generateAllGridPoints();

    // Is the symmetry of the reciprocal lattice calculated?
    if( !xstr.pgroupk_xtal_calculated )
    {
        ofstream fileDevNull("/dev/null");
        if( !fileDevNull.is_open() )
        {
            throw APLRuntimeError("MonkhorstPackMesh::MonkhorstPackMesh(): Cannot open output stream /dev/null.");
        }

        _aflags af;
        af.Directory = "./";
//        af.hostname = "";
        af.QUIET = TRUE;
        xstr.LatticeReduction_avoid = TRUE;

        //SYM::CalculatePointGroupKlattice(fileDevNull, xstr, af, FALSE, FALSE, _logger.getOutputStream()); //CO 171215 pgroupk_xtal NOT pgroupk
        SYM::CalculatePointGroupKCrystal(fileDevNull, xstr, af, FALSE, FALSE, _logger.getOutputStream());
        xstr.pgroupk_xtal_calculated = TRUE;

        fileDevNull.clear();
        fileDevNull.close();
    }

    // OK. Do reduction to irreducible representation, or let it be in case of some problems...
    if( xstr.pgroupk_xtal_calculated ) {
      makeIrreducible(xstr.pgroupk_xtal);
      _isIrreducible = true;
    } else {
      _isIrreducible = false;
      _logger << apl::warning  << "The points group of the reciprocal lattice is missing. No reduction to the irreducible representation is done." << apl::endl;
    }
}

// ///////////////////////////////////////////////////////////////////////////

MonkhorstPackMesh::~MonkhorstPackMesh()
{
    this->clear();
}

// ///////////////////////////////////////////////////////////////////////////

void MonkhorstPackMesh::clear()
{
    xvector<double> zero(3);
    zero(1) = 0; zero(2) = 0; zero(3) = 0;
    _shift = zero;

    xvector<int> nzero(3);
    nzero(1) = 0; nzero(2) = 0; nzero(3) = 0;
    _n = nzero;

    _allToIrrPointMap.clear();

    _kpoints.clear();
    _weights.clear();

     xmatrix<double> zeroMatrix(3,3);
    _rlattice = zeroMatrix;
    _klattice = zeroMatrix;

    _isIrreducible = false;
}

// ///////////////////////////////////////////////////////////////////////////

void MonkhorstPackMesh::setDensity(int n)
{
    setDensity(n,n,n);
}

void MonkhorstPackMesh::setDensity(int na, int nb, int nc)
{
    _n(1) = na;
    _n(2) = nb;
    _n(3) = nc;

    // Setup also map
    xtensor3<int> temp(na,nb,nc,1,1,1);
    _allToIrrPointMap = temp;

    // TODO: Be aware of the bug in aflow_xtensor.cpp in operator= function...
    if( _allToIrrPointMap.uindex[1] == 3 ||
        _allToIrrPointMap.uindex[2] == 3 ||
        _allToIrrPointMap.uindex[3] == 3 )
    {
        throw APLRuntimeError("Be aware of the bug in aflow_xtensor.cpp in operator= function...");
    }
}

// ///////////////////////////////////////////////////////////////////////////

void MonkhorstPackMesh::generateAllGridPoints()
{
    // Clear old mesh
    _kpoints.clear();
    _weights.clear();

    /*    
    xvector<double> kpoint(3);
    int id = 0;
    for(int s = 1; s <= _n(3); s++)
    for(int r = 1; r <= _n(2); r++)
    for(int p = 1; p <= _n(1); p++)
    {
        kpoint(1) = -0.5 + (double)(p-1)/(_n(1)-1);
        kpoint(2) = -0.5 + (double)(r-1)/(_n(2)-1);
        kpoint(3) = -0.5 + (double)(s-1)/(_n(3)-1);
        // Transform to cartesian coordinate
        kpoint = trasp(_klattice) * kpoint;
        _kpoints.push_back(kpoint);
        _weights.push_back(1.0);
        _allToIrrPointMap(p,r,s) = id++;
    }
    */

    // Monkhorst-Pack formula
    xvector<double> kpoint(3);
    int id = 0;
    for(int s = 1; s <= _n(3); s++)
    for(int r = 1; r <= _n(2); r++)
    for(int p = 1; p <= _n(1); p++)
    {
        kpoint(1) = ( 2.0 * p - _n(1) - 1 ) / ( 2.0 * _n(1) );
        kpoint(2) = ( 2.0 * r - _n(2) - 1 ) / ( 2.0 * _n(2) );
        kpoint(3) = ( 2.0 * s - _n(3) - 1 ) / ( 2.0 * _n(3) );
        // Transform to cartesian coordinate
        kpoint = trasp(_klattice) * kpoint;
        _kpoints.push_back(kpoint);
        _weights.push_back(1.0);
        _allToIrrPointMap(p,r,s) = id++;
    }
}

// ///////////////////////////////////////////////////////////////////////////

void MonkhorstPackMesh::makeIrreducible(const vector<_sym_op>& pgroupk_xtal)
{
    vector< xvector<double> > ibzkpts;
    vector<double> weights;
    vector<int> map;

    _logger.initProgressBar("Counting the irreducible k-points");

    //
    xvector<double> newkpoint(3);
    for(uint i = 0; i < _kpoints.size(); i++)
    {
        _logger.updateProgressBar(i,_kpoints.size());
        bool isNew = true;
        for(uint symOpID = 0; symOpID < pgroupk_xtal.size(); symOpID++)
        {
            newkpoint = pgroupk_xtal[symOpID].Uc * _kpoints[i];

            uint j = 0;
            for( ; j < ibzkpts.size(); j++)
            {
                if( ( sign( newkpoint(1) ) == sign( ibzkpts[j](1) ) ) && ( fabs( fabs( newkpoint(1) ) - fabs( ibzkpts[j](1) ) ) < _AFLOW_APL_EPS_ ) &&
                    ( sign( newkpoint(2) ) == sign( ibzkpts[j](2) ) ) && ( fabs( fabs( newkpoint(2) ) - fabs( ibzkpts[j](2) ) ) < _AFLOW_APL_EPS_ ) &&
                    ( sign( newkpoint(3) ) == sign( ibzkpts[j](3) ) ) && ( fabs( fabs( newkpoint(3) ) - fabs( ibzkpts[j](3) ) ) < _AFLOW_APL_EPS_ ) )
                {
                    break;
                }
            }

            // OK. Not a new kpoint, next one...
            if( j != ibzkpts.size() )
            {
                weights[j]++;
                isNew = false;
                map.push_back(j);
                break;
            }
        }
        // OK. NEW!
        if( isNew )
        {
            ibzkpts.push_back(newkpoint);
            weights.push_back(1);
            map.push_back(ibzkpts.size()-1);
        }
    }
    _logger.finishProgressBar();

    _logger << "Found " << (int)(ibzkpts.size()) << " irreducible qpoints." << apl::endl;
    //for(uint i = 0; i < ibzkpts.size(); i++)
    //    { cout << "[" << i << "] "; printXVector(ibzkpts[i],false); cout << " weight = " << weights[i] << std::endl; }

    // Copy
    _kpoints.clear();
    _weights.clear();
    for(uint i = 0; i < ibzkpts.size(); i++)
    {
        _kpoints.push_back(ibzkpts[i]);
        _weights.push_back(weights[i]);
    }
    ibzkpts.clear();
    weights.clear();

    // Update map
    int id = 0;
    for(int s = 1; s <= _n(3); s++)
    for(int r = 1; r <= _n(2); r++)
    for(int p = 1; p <= _n(1); p++)
        _allToIrrPointMap(p,r,s) = map[id++];
    map.clear();
    
}

// INTERFACE /////////////////////////////////////////////////////////////////

int MonkhorstPackMesh::getN(int i)
{
    return _n(i);
}

// ///////////////////////////////////////////////////////////////////////////

aurostd::xmatrix<double> MonkhorstPackMesh::getReciprocalLattice()
{
    return _klattice;
}

// ///////////////////////////////////////////////////////////////////////////

std::vector< aurostd::xvector<double> > MonkhorstPackMesh::getPoints()
{
    return _kpoints;
}

// ///////////////////////////////////////////////////////////////////////////

std::vector<double> MonkhorstPackMesh::getWeights()
{
    return _weights;
}

// ///////////////////////////////////////////////////////////////////////////

aurostd::xtensor3<int> MonkhorstPackMesh::getAllToIrrPointMap()
{
    return _allToIrrPointMap;
}

// ///////////////////////////////////////////////////////////////////////////

} // namespace apl
