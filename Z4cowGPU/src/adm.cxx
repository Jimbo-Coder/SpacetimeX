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

extern "C" void Z4cowGPU_ADM(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Z4cowGPU_ADM;
  DECLARE_CCTK_PARAMETERS;

  const GF3D2layout layout2(cctkGH, {0, 0, 0});

  // Input grid functions
  const std::array<const CCTK_REAL *, 6> gf_gamt{gammatxx, gammatxy, gammatxz,
                                            gammatyy, gammatyz, gammatzz};
  const CCTK_REAL *gf_exKh = Kh;
  const std::array<const CCTK_REAL *, 6> gf_exAt{Atxx, Atxy, Atxz, Atyy, Atyz, Atzz};
  const CCTK_REAL *gf_alpha = alphaG;
  const std::array<const CCTK_REAL *, 3> gf_beta{betaGx, betaGy, betaGz};

  // Output grid functions
  const std::array<CCTK_REAL *, 6> gf_ADMgam{gxx, gxy, gxz, gyy, gyz, gzz};
  const std::array<CCTK_REAL *, 6> gf_ADMK{kxx, kxy, kxz, kyy, kyz, kzz};
  CCTK_REAL *gf_ADMalpha = alp;
  const std::array<CCTK_REAL *, 3> gf_ADMbeta{betax, betay, betaz};

  // Loop
  const Loop::GridDescBaseDevice grid(cctkGH);

#include "../wolfram/Z4cowGPU_set_ADM.hxx"
}

} // namespace Z4cowGPU
