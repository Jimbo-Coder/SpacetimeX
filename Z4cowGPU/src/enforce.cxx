#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <loop_device.hxx>

#include <array>
#include <cmath>

#include <cx_powerinline.hxx>

namespace Z4cowGPU {
using namespace Arith;
using namespace Loop;
using namespace CXUtils;

extern "C" void Z4cowGPU_Enforce(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Z4cowGPU_Enforce;
  DECLARE_CCTK_PARAMETERS;

  const GF3D2layout layout2(cctkGH, {0, 0, 0});

  // Output grid functions
  const std::array<CCTK_REAL *, 6> gf_gamt{gammatxx, gammatxy, gammatxz,
                                      gammatyy, gammatyz, gammatzz};
  const std::array<CCTK_REAL *, 6> gf_exAt{Atxx, Atxy, Atxz, Atyy, Atyz, Atzz};
  CCTK_REAL *gf_alpha = alphaG;

  // Parameters
  const CCTK_REAL cWfloor = W_floor;
  const CCTK_REAL calphafloor = alphaG_floor;

  // Loop
  const Loop::GridDescBaseDevice grid(cctkGH);

#include "../wolfram/Z4cowGPU_enforce.hxx"
}

} // namespace Z4cowGPU
