#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <loop_device.hxx>
#include <cx_utils.hxx>

#include <array>
#include <cmath>

#include <cx_derivsGF3D5.hxx>
#include <cx_derivsinline.hxx>
#include <cx_dissinline.hxx>
#include <cx_powerinline.hxx>

namespace Z4cowGPU {
using namespace Arith;
using namespace Loop;
using namespace CXUtils;

extern "C" void Z4cowGPU_RHS(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Z4cowGPU_RHS;
  DECLARE_CCTK_PARAMETERS;

  for (int d = 0; d < 3; ++d)
    if (cctk_nghostzones[d] < deriv_order / 2)
      CCTK_VERROR("Need at least %d ghost zones", deriv_order / 2);

  const std::array<CCTK_REAL, 3> invDxyz{1. / CCTK_DELTA_SPACE(0),
                                    1. / CCTK_DELTA_SPACE(1),
                                    1. / CCTK_DELTA_SPACE(2)};

  const GF3D2layout layout2(cctkGH, {0, 0, 0});

  // Input grid functions
  const CCTK_REAL *gf_W = W;
  const std::array<const CCTK_REAL *, 6> gf_gamt{gammatxx, gammatxy, gammatxz,
                                            gammatyy, gammatyz, gammatzz};
  const CCTK_REAL *gf_exKh = Kh;
  const std::array<const CCTK_REAL *, 6> gf_exAt{Atxx, Atxy, Atxz, Atyy, Atyz, Atzz};
  const std::array<const CCTK_REAL *, 3> gf_trGt{Gamtx, Gamty, Gamtz};
  const CCTK_REAL *gf_Theta = Theta;
  const CCTK_REAL *gf_alpha = alphaG;
  const std::array<const CCTK_REAL *, 3> gf_beta{betaGx, betaGy, betaGz};

  // More input grid functions
  // const CCTK_REAL *gf_eTtt = eTtt;
  const std::array<const CCTK_REAL *, 3> gf_eTt{eTtx, eTty, eTtz};
  const std::array<const CCTK_REAL *, 6> gf_eT{eTxx, eTxy, eTxz, eTyy, eTyz, eTzz};

  // Output grid functions
  CCTK_REAL *gf_dtW = W_rhs;
  const std::array<CCTK_REAL *, 6> gf_dtgamt{gammatxx_rhs, gammatxy_rhs,
                                        gammatxz_rhs, gammatyy_rhs,
                                        gammatyz_rhs, gammatzz_rhs};
  CCTK_REAL *gf_dtexKh = Kh_rhs;
  const std::array<CCTK_REAL *, 6> gf_dtexAt{Atxx_rhs, Atxy_rhs, Atxz_rhs,
                                        Atyy_rhs, Atyz_rhs, Atzz_rhs};
  const std::array<CCTK_REAL *, 3> gf_dttrGt{Gamtx_rhs, Gamty_rhs, Gamtz_rhs};
  CCTK_REAL *gf_dtTheta = Theta_rhs;
  CCTK_REAL *gf_dtalpha = alphaG_rhs;
  const std::array<CCTK_REAL *, 3> gf_dtbeta{betaGx_rhs, betaGy_rhs, betaGz_rhs};

  // Parameters
  const CCTK_REAL cpi = M_PI;
  const CCTK_REAL ckappa1 = kappa1;
  const CCTK_REAL ckappa2 = kappa2;
  const CCTK_REAL cmuL = f_mu_L;
  const CCTK_REAL cmuS = f_mu_S;
  // const CCTK_REAL ceta = eta;
  const auto calceta =
      [=] CCTK_DEVICE(const CCTK_REAL x, const CCTK_REAL y, const CCTK_REAL z)
          CCTK_ATTRIBUTE_ALWAYS_INLINE {
            const CCTK_REAL r2 = x * x + y * y + z * z;
            const CCTK_REAL r4 = r2 * r2;
            const CCTK_REAL is4 =
                1.0 / (veta_width * veta_width * veta_width * veta_width);
            return (veta_central - veta_outer) * exp(-r4 * is4) + veta_outer;
          };

  // Loop
  const Loop::GridDescBaseDevice grid(cctkGH);

  if (calc_derivs_live) {

    // Derivs Lambdas
#include "../wolfram/Z4cowGPU_derivs1st.hxx"
#include "../wolfram/Z4cowGPU_derivs2nd.hxx"

#include "../wolfram/Z4cowGPU_set_rhs.hxx"

  } else {

    // Tile variables for derivatives and so on
    const int ntmps = 132;
    int itmp = 0;
    const GF3D5layout layout5 = CXUtils::get_GF3D5layout<0, 0, 0>(cctkGH);
    CXUtils::GF3D5Factory<CCTK_REAL> fct(layout5, ntmps, itmp);

    const auto tl_dW = fct.make_vec_gf();
    const auto tl_ddW = fct.make_smat_gf();
    const auto tl_dgamt = fct.make_smat_vec_gf();
    const auto tl_ddgamt = fct.make_smat_smat_gf();
    const auto tl_dexKh = fct.make_vec_gf();
    const auto tl_dexAt = fct.make_smat_vec_gf();
    const auto tl_dtrGt = fct.make_vec_vec_gf();
    const auto tl_dTheta = fct.make_vec_gf();
    const auto tl_dalpha = fct.make_vec_gf();
    const auto tl_ddalpha = fct.make_smat_gf();
    const auto tl_dbeta = fct.make_vec_vec_gf();
    const auto tl_ddbeta = fct.make_vec_smat_gf();

    if (itmp != ntmps)
      CCTK_VERROR("Wrong number of temporary variables: ntmps=%d itmp=%d",
                  ntmps, itmp);
    itmp = -1;

    // Define derivs lambdas
    const auto calcderivs1st = [&](const auto &dgf, const auto &gf_) {
      calc_derivs1st<0, 0, 0>(grid, layout5, dgf, layout2, gf_, invDxyz,
                              deriv_order);
    };
    const auto calcderivs2nd = [&](const auto &dgf, const auto &ddgf,
                                   const auto &gf_) {
      calc_derivs2nd<0, 0, 0>(grid, layout5, dgf, ddgf, layout2, gf_, invDxyz,
                              deriv_order);
    };

    calcderivs2nd(tl_dW, tl_ddW, gf_W);
    calcderivs2nd(tl_dgamt, tl_ddgamt, gf_gamt);
    calcderivs1st(tl_dexKh, gf_exKh);
    calcderivs1st(tl_dexAt, gf_exAt);
    calcderivs1st(tl_dtrGt, gf_trGt);
    calcderivs1st(tl_dTheta, gf_Theta);
    calcderivs2nd(tl_dalpha, tl_ddalpha, gf_alpha);
    calcderivs2nd(tl_dbeta, tl_ddbeta, gf_beta);

#include "../wolfram/Z4cowGPU_set_rhs_GF3D5.hxx"
  }

  // Dissipation
#include "../wolfram/Z4cowGPU_applydiss.hxx"

  applydiss(gf_W, gf_dtW);
  for (int a = 0; a < 6; ++a)
    applydiss(gf_gamt[a], gf_dtgamt[a]);
  applydiss(gf_exKh, gf_dtexKh);
  for (int a = 0; a < 6; ++a)
    applydiss(gf_exAt[a], gf_dtexAt[a]);
  for (int a = 0; a < 3; ++a)
    applydiss(gf_trGt[a], gf_dttrGt[a]);
  if (!set_Theta_zero)
    applydiss(gf_Theta, gf_dtTheta);
  applydiss(gf_alpha, gf_dtalpha);
  for (int a = 0; a < 3; ++a)
    applydiss(gf_beta[a], gf_dtbeta[a]);
}

extern "C" void Z4cowGPU_Sync(CCTK_ARGUMENTS) {
  // do nothing
}

} // namespace Z4cowGPU
