#include "aflow_apl.h"

using namespace std;

namespace apl {

// ///////////////////////////////////////////////////////////////////////////

GeneralizedSupercellApproach::GeneralizedSupercellApproach(
    Supercell& sc, vector<ClusterSet>& clst, _xinput& xinput,
    _aflags& aflags, _kflags& kflags,
    _xflags& xflags, //_vflags& vflags, 
    string& AflowIn, Logger& l)
    : DirectMethodPC(sc, clst, xinput, aflags, kflags, xflags, AflowIn, l) {
}

// ///////////////////////////////////////////////////////////////////////////

GeneralizedSupercellApproach::~GeneralizedSupercellApproach() {
  clear();
}

// ///////////////////////////////////////////////////////////////////////////

void GeneralizedSupercellApproach::clear() {
}

//////////////////////////////////////////////////////////////////////////////

}  // namespace apl
