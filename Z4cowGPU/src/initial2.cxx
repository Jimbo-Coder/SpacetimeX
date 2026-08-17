#include <cctk.h>
#include <cctk_Arguments.h>
#include <loop_device.hxx>

#include <array>
#include <cmath>

#include <cx_derivsinline.hxx>
#include <cx_powerinline.hxx>

namespace Z4cowGPU {
using namespace Arith;
using namespace Loop;
using namespace CXUtils;

extern "C" void Z4cowGPU_Initial2(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Z4cowGPU_Initial2;
  DECLARE_CCTK_PARAMETERS;

  const std::array<CCTK_REAL, 3> invDxyz{1. / CCTK_DELTA_SPACE(0),
                                    1. / CCTK_DELTA_SPACE(1),
                                    1. / CCTK_DELTA_SPACE(2)};

  const GF3D2layout layout2(cctkGH, {0, 0, 0});

  // Input grid functions
  const std::array<const CCTK_REAL *, 6> gf_gamt{gammatxx, gammatxy, gammatxz,
                                            gammatyy, gammatyz, gammatzz};

  // Output grid functions
  const std::array<CCTK_REAL *, 3> gf_trGt{Gamtx, Gamty, Gamtz};

  // Derivs Lambdas
#include "../wolfram/Z4cowGPU_derivs1st.hxx"

  // Loop
  const Loop::GridDescBaseDevice grid(cctkGH);

#include "../wolfram/Z4cowGPU_initial2.hxx"
}

} // namespace Z4cowGPU
