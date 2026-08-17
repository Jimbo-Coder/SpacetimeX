#include <cctk.h>
#include <cctk_Arguments.h>
#include <loop_device.hxx>

#include <array>
#include <cmath>

#include <cx_powerinline.hxx>

namespace Z4cowGPU {
using namespace Arith;
using namespace Loop;
using namespace CXUtils;

extern "C" void Z4cowGPU_Initial1(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Z4cowGPU_Initial1;

  const GF3D2layout layout2(cctkGH, {0, 0, 0});

  // Input grid functions
  const std::array<const CCTK_REAL *, 6> gf_ADMgam{gxx, gxy, gxz, gyy, gyz, gzz};
  const std::array<const CCTK_REAL *, 6> gf_ADMK{kxx, kxy, kxz, kyy, kyz, kzz};
  const CCTK_REAL *gf_ADMalpha = alp;
  const std::array<const CCTK_REAL *, 3> gf_ADMbeta{betax, betay, betaz};

  // Output grid functions
  const std::array<CCTK_REAL *, 6> gf_gamt{gammatxx, gammatxy, gammatxz,
                                      gammatyy, gammatyz, gammatzz};
  CCTK_REAL *gf_exKh = Kh;
  const std::array<CCTK_REAL *, 6> gf_exAt{Atxx, Atxy, Atxz, Atyy, Atyz, Atzz};
  CCTK_REAL *gf_alpha = alphaG;
  const std::array<CCTK_REAL *, 3> gf_beta{betaGx, betaGy, betaGz};

  // Loop
  const Loop::GridDescBaseDevice grid(cctkGH);

#include "../wolfram/Z4cowGPU_initial1.hxx"
}

} // namespace Z4cowGPU
