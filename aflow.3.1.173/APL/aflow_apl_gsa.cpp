#include "aflow_apl.h"

using namespace std;

namespace apl {

// ///////////////////////////////////////////////////////////////////////////

GeneralizedSupercellApproach::GeneralizedSupercellApproach(
    Supercell& sc, StrPairs& strpair, _xinput& xinput, //_xvasp& xvasp, 
    _aflags& aflags, _kflags& kflags,
    _xflags& xflags, //_vflags& vflags, 
    string& AflowIn,
    Logger& l)
    : DirectMethodPC(sc, strpair, xinput, aflags, kflags, xflags, AflowIn, l) { //xvasp, aflags, kflags, vflags, l) {
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
