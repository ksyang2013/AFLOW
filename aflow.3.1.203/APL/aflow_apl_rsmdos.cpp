#include "aflow_apl.h"

using namespace std;

namespace apl
{

// ///////////////////////////////////////////////////////////////////////////

RootSamplingMethod::RootSamplingMethod(IPhononCalculator& pc, IReciprocalPointGrid& rg, Logger& l)
                                      : DOSCalculator(pc,rg,l)
{
}

// ///////////////////////////////////////////////////////////////////////////

RootSamplingMethod::~RootSamplingMethod()
{
}

// ///////////////////////////////////////////////////////////////////////////

void RootSamplingMethod::rawCalc(int USER_DOS_NPOINTS)
{
    for(int k = 0; k < USER_DOS_NPOINTS; k++)
    {
        for(uint i = 0; i < _freqs.size(); i++)
        {
            for(int j = _freqs[i].lrows; j <= _freqs[i].urows; j++)
            {
                if( ( _freqs[i](j) > ( _bins[k] - _halfStepDOS ) ) &&
                    ( _freqs[i](j) < ( _bins[k] + _halfStepDOS ) ) )
                {
                    _dos[k] += _qweights[i];
                }
            }
        }
    }
}

// ///////////////////////////////////////////////////////////////////////////

} // namespace apl
