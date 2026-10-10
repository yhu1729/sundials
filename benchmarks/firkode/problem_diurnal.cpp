/* -----------------------------------------------------------------
 * Programmer(s): Yifan Hu @ UMBC
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2025-2026, Lawrence Livermore National Security,
 * University of Maryland Baltimore County, and the SUNDIALS contributors.
 * Copyright (c) 2013-2025, Lawrence Livermore National Security
 * and Southern Methodist University.
 * Copyright (c) 2002-2013, Lawrence Livermore National Security.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * 2-species diurnal kinetics advection-diffusion PDE on an MX x MY
 * grid (port of examples/firkode/C_serial/firk_diurnal_kry.c with a
 * run-time grid). Matrix-free: Jacobian-vector product and a
 * block-diagonal (2x2 per grid point) preconditioner.
 * -----------------------------------------------------------------*/

#include <cmath>
#include <vector>

#include <nvector/nvector_serial.h>
#include <sundials/sundials_dense.h>

#include "problem.hpp"

namespace benchmark {

namespace {

constexpr sunrealtype ZERO     = SUN_RCONST(0.0);
constexpr sunrealtype ONE      = SUN_RCONST(1.0);
constexpr sunrealtype TWO      = SUN_RCONST(2.0);
constexpr int NUM_SPECIES      = 2;
constexpr sunrealtype KH       = SUN_RCONST(4.0e-6);
constexpr sunrealtype VEL      = SUN_RCONST(0.001);
constexpr sunrealtype KV0      = SUN_RCONST(1.0e-8);
constexpr sunrealtype Q1       = SUN_RCONST(1.63e-16);
constexpr sunrealtype Q2       = SUN_RCONST(4.66e-16);
constexpr sunrealtype C3       = SUN_RCONST(3.7e16);
constexpr sunrealtype A3       = SUN_RCONST(22.62);
constexpr sunrealtype A4       = SUN_RCONST(7.601);
constexpr sunrealtype C1_SCALE = SUN_RCONST(1.0e6);
constexpr sunrealtype C2_SCALE = SUN_RCONST(1.0e12);
constexpr sunrealtype HALFDAY  = SUN_RCONST(4.32e4);
constexpr sunrealtype PI       = SUN_RCONST(3.1415926535898);
constexpr sunrealtype XMIN     = ZERO;
constexpr sunrealtype XMAX     = SUN_RCONST(20.0);
constexpr sunrealtype YMIN     = SUN_RCONST(30.0);
constexpr sunrealtype YMAX     = SUN_RCONST(50.0);
constexpr sunrealtype XMID     = SUN_RCONST(10.0);
constexpr sunrealtype YMID     = SUN_RCONST(40.0);

struct DiurnalData
{
  int mx = 10, my = 10;
  sunrealtype q4 = ZERO, om = ZERO, dx = ZERO, dy = ZERO;
  sunrealtype hdco = ZERO, haco = ZERO, vdco = ZERO;
  std::vector<sunrealtype**> P, Jbd;
  std::vector<sunindextype*> pivot;

  sunrealtype& at(sunrealtype* v, int i, int jx, int jy) const
  {
    return v[i - 1 + jx * NUM_SPECIES + jy * NUM_SPECIES * mx];
  }

  int blk(int jx, int jy) const { return jx + jy * mx; }

  void allocate(int mx_, int my_)
  {
    release();
    mx = mx_;
    my = my_;
    P.resize((size_t)mx * my);
    Jbd.resize((size_t)mx * my);
    pivot.resize((size_t)mx * my);
    for (size_t k = 0; k < P.size(); k++)
    {
      P[k]     = SUNDlsMat_newDenseMat(NUM_SPECIES, NUM_SPECIES);
      Jbd[k]   = SUNDlsMat_newDenseMat(NUM_SPECIES, NUM_SPECIES);
      pivot[k] = SUNDlsMat_newIndexArray(NUM_SPECIES);
    }
    om   = PI / HALFDAY;
    dx   = (XMAX - XMIN) / (mx - 1);
    dy   = (YMAX - YMIN) / (my - 1);
    hdco = KH / (dx * dx);
    haco = VEL / (TWO * dx);
    vdco = (ONE / (dy * dy)) * KV0;
  }

  void release()
  {
    for (size_t k = 0; k < P.size(); k++)
    {
      SUNDlsMat_destroyMat(P[k]);
      SUNDlsMat_destroyMat(Jbd[k]);
      SUNDlsMat_destroyArray(pivot[k]);
    }
    P.clear();
    Jbd.clear();
    pivot.clear();
  }
};

int diurnal_rhs(sunrealtype t, N_Vector u, N_Vector udot, void* user_data)
{
  DiurnalData* d      = static_cast<DiurnalData*>(user_data);
  sunrealtype* udata  = N_VGetArrayPointer(u);
  sunrealtype* dudata = N_VGetArrayPointer(udot);

  const sunrealtype s = std::sin(d->om * t);
  sunrealtype q3;
  if (s > ZERO)
  {
    q3    = std::exp(-A3 / s);
    d->q4 = std::exp(-A4 / s);
  }
  else
  {
    q3    = ZERO;
    d->q4 = ZERO;
  }
  const sunrealtype q4coef = d->q4, dely = d->dy, verdco = d->vdco,
                    hordco = d->hdco, horaco = d->haco;

  for (int jy = 0; jy < d->my; jy++)
  {
    const sunrealtype ydn  = YMIN + (jy - SUN_RCONST(0.5)) * dely;
    const sunrealtype yup  = ydn + dely;
    const sunrealtype cydn = verdco * std::exp(SUN_RCONST(0.2) * ydn);
    const sunrealtype cyup = verdco * std::exp(SUN_RCONST(0.2) * yup);
    const int idn          = (jy == 0) ? 1 : -1;
    const int iup          = (jy == d->my - 1) ? -1 : 1;
    for (int jx = 0; jx < d->mx; jx++)
    {
      const sunrealtype c1  = d->at(udata, 1, jx, jy);
      const sunrealtype c2  = d->at(udata, 2, jx, jy);
      const sunrealtype qq1 = Q1 * c1 * C3, qq2 = Q2 * c1 * c2, qq3 = q3 * C3,
                        qq4   = q4coef * c2;
      const sunrealtype rkin1 = -qq1 - qq2 + TWO * qq3 + qq4;
      const sunrealtype rkin2 = qq1 - qq2 - qq4;

      const sunrealtype c1dn   = d->at(udata, 1, jx, jy + idn);
      const sunrealtype c2dn   = d->at(udata, 2, jx, jy + idn);
      const sunrealtype c1up   = d->at(udata, 1, jx, jy + iup);
      const sunrealtype c2up   = d->at(udata, 2, jx, jy + iup);
      const sunrealtype vertd1 = cyup * (c1up - c1) - cydn * (c1 - c1dn);
      const sunrealtype vertd2 = cyup * (c2up - c2) - cydn * (c2 - c2dn);

      const int ileft          = (jx == 0) ? 1 : -1;
      const int iright         = (jx == d->mx - 1) ? -1 : 1;
      const sunrealtype c1lt   = d->at(udata, 1, jx + ileft, jy);
      const sunrealtype c2lt   = d->at(udata, 2, jx + ileft, jy);
      const sunrealtype c1rt   = d->at(udata, 1, jx + iright, jy);
      const sunrealtype c2rt   = d->at(udata, 2, jx + iright, jy);
      const sunrealtype hord1  = hordco * (c1rt - TWO * c1 + c1lt);
      const sunrealtype hord2  = hordco * (c2rt - TWO * c2 + c2lt);
      const sunrealtype horad1 = horaco * (c1rt - c1lt);
      const sunrealtype horad2 = horaco * (c2rt - c2lt);

      d->at(dudata, 1, jx, jy) = vertd1 + hord1 + horad1 + rkin1;
      d->at(dudata, 2, jx, jy) = vertd2 + hord2 + horad2 + rkin2;
    }
  }
  return 0;
}

int diurnal_jtv(N_Vector v, N_Vector Jv, sunrealtype t, N_Vector u, N_Vector,
                void* user_data, N_Vector)
{
  DiurnalData* d      = static_cast<DiurnalData*>(user_data);
  sunrealtype* udata  = N_VGetArrayPointer(u);
  sunrealtype* vdata  = N_VGetArrayPointer(v);
  sunrealtype* Jvdata = N_VGetArrayPointer(Jv);

  const sunrealtype s      = std::sin(d->om * t);
  d->q4                    = (s > ZERO) ? std::exp(-A4 / s) : ZERO;
  const sunrealtype q4coef = d->q4, dely = d->dy, verdco = d->vdco,
                    hordco = d->hdco, horaco = d->haco;

  for (int jy = 0; jy < d->my; jy++)
  {
    const sunrealtype ydn  = YMIN + (jy - SUN_RCONST(0.5)) * dely;
    const sunrealtype yup  = ydn + dely;
    const sunrealtype cydn = verdco * std::exp(SUN_RCONST(0.2) * ydn);
    const sunrealtype cyup = verdco * std::exp(SUN_RCONST(0.2) * yup);
    const int idn          = (jy == 0) ? 1 : -1;
    const int iup          = (jy == d->my - 1) ? -1 : 1;
    for (int jx = 0; jx < d->mx; jx++)
    {
      const sunrealtype c1   = d->at(udata, 1, jx, jy);
      const sunrealtype c2   = d->at(udata, 2, jx, jy);
      const sunrealtype v1   = d->at(vdata, 1, jx, jy);
      const sunrealtype v2   = d->at(vdata, 2, jx, jy);
      const sunrealtype v1dn = d->at(vdata, 1, jx, jy + idn);
      const sunrealtype v2dn = d->at(vdata, 2, jx, jy + idn);
      const sunrealtype v1up = d->at(vdata, 1, jx, jy + iup);
      const sunrealtype v2up = d->at(vdata, 2, jx, jy + iup);
      const int ileft        = (jx == 0) ? 1 : -1;
      const int iright       = (jx == d->mx - 1) ? -1 : 1;
      const sunrealtype v1lt = d->at(vdata, 1, jx + ileft, jy);
      const sunrealtype v2lt = d->at(vdata, 2, jx + ileft, jy);
      const sunrealtype v1rt = d->at(vdata, 1, jx + iright, jy);
      const sunrealtype v2rt = d->at(vdata, 2, jx + iright, jy);

      sunrealtype Jv1 = -(Q1 * C3 + Q2 * c2) * v1 + (q4coef - Q2 * c1) * v2;
      sunrealtype Jv2 = (Q1 * C3 - Q2 * c2) * v1 - (q4coef + Q2 * c1) * v2;
      Jv1 += -(cyup + cydn) * v1 + cyup * v1up + cydn * v1dn;
      Jv2 += -(cyup + cydn) * v2 + cyup * v2up + cydn * v2dn;
      Jv1 += hordco * (v1rt - TWO * v1 + v1lt);
      Jv2 += hordco * (v2rt - TWO * v2 + v2lt);
      Jv1 += horaco * (v1rt - v1lt);
      Jv2 += horaco * (v2rt - v2lt);
      d->at(Jvdata, 1, jx, jy) = Jv1;
      d->at(Jvdata, 2, jx, jy) = Jv2;
    }
  }
  return 0;
}

int diurnal_psetup(sunrealtype, N_Vector u, N_Vector, sunbooleantype jok,
                   sunbooleantype* jcurPtr, sunrealtype gamma, void* user_data)
{
  DiurnalData* d     = static_cast<DiurnalData*>(user_data);
  sunrealtype* udata = N_VGetArrayPointer(u);

  if (jok)
  {
    for (int jy = 0; jy < d->my; jy++)
    {
      for (int jx = 0; jx < d->mx; jx++)
      {
        SUNDlsMat_denseCopy(d->Jbd[(size_t)d->blk(jx, jy)],
                            d->P[(size_t)d->blk(jx, jy)], NUM_SPECIES,
                            NUM_SPECIES);
      }
    }
    *jcurPtr = SUNFALSE;
  }
  else
  {
    const sunrealtype q4coef = d->q4, dely = d->dy, verdco = d->vdco,
                      hordco = d->hdco;
    for (int jy = 0; jy < d->my; jy++)
    {
      const sunrealtype ydn  = YMIN + (jy - SUN_RCONST(0.5)) * dely;
      const sunrealtype yup  = ydn + dely;
      const sunrealtype cydn = verdco * std::exp(SUN_RCONST(0.2) * ydn);
      const sunrealtype cyup = verdco * std::exp(SUN_RCONST(0.2) * yup);
      const sunrealtype diag = -(cydn + cyup + TWO * hordco);
      for (int jx = 0; jx < d->mx; jx++)
      {
        const sunrealtype c1 = d->at(udata, 1, jx, jy);
        const sunrealtype c2 = d->at(udata, 2, jx, jy);
        sunrealtype** j      = d->Jbd[(size_t)d->blk(jx, jy)];
        sunrealtype** a      = d->P[(size_t)d->blk(jx, jy)];
        j[0][0]              = (-Q1 * C3 - Q2 * c2) + diag;
        j[1][0]              = -Q2 * c1 + q4coef;
        j[0][1]              = Q1 * C3 - Q2 * c2;
        j[1][1]              = (-Q2 * c1 - q4coef) + diag;
        SUNDlsMat_denseCopy(j, a, NUM_SPECIES, NUM_SPECIES);
      }
    }
    *jcurPtr = SUNTRUE;
  }

  for (size_t k = 0; k < d->P.size(); k++)
  {
    SUNDlsMat_denseScale(-gamma, d->P[k], NUM_SPECIES, NUM_SPECIES);
    SUNDlsMat_denseAddIdentity(d->P[k], NUM_SPECIES);
    if (SUNDlsMat_denseGETRF(d->P[k], NUM_SPECIES, NUM_SPECIES, d->pivot[k]) != 0)
    {
      return 1;
    }
  }
  return 0;
}

int diurnal_psolve(sunrealtype, N_Vector, N_Vector, N_Vector r, N_Vector z,
                   sunrealtype, sunrealtype, int, void* user_data)
{
  DiurnalData* d     = static_cast<DiurnalData*>(user_data);
  sunrealtype* zdata = N_VGetArrayPointer(z);
  N_VScale(ONE, r, z);
  for (int jx = 0; jx < d->mx; jx++)
  {
    for (int jy = 0; jy < d->my; jy++)
    {
      sunrealtype* v = &(d->at(zdata, 1, jx, jy));
      SUNDlsMat_denseGETRS(d->P[(size_t)d->blk(jx, jy)], NUM_SPECIES,
                           d->pivot[(size_t)d->blk(jx, jy)], v);
    }
  }
  return 0;
}

// IMEX splitting "imex": the kinetics are implicit, the transport
// (vertical diffusion, horizontal diffusion and advection) explicit.
// Separate functions; the kinetics set q4 as diurnal_rhs does.
int diurnal_kinetics(sunrealtype t, N_Vector u, N_Vector udot, void* user_data)
{
  DiurnalData* d      = static_cast<DiurnalData*>(user_data);
  sunrealtype* udata  = N_VGetArrayPointer(u);
  sunrealtype* dudata = N_VGetArrayPointer(udot);
  const sunrealtype s = std::sin(d->om * t);
  sunrealtype q3;
  if (s > ZERO)
  {
    q3    = std::exp(-A3 / s);
    d->q4 = std::exp(-A4 / s);
  }
  else
  {
    q3    = ZERO;
    d->q4 = ZERO;
  }
  const sunrealtype q4coef = d->q4;
  for (int jy = 0; jy < d->my; jy++)
  {
    for (int jx = 0; jx < d->mx; jx++)
    {
      const sunrealtype c1  = d->at(udata, 1, jx, jy);
      const sunrealtype c2  = d->at(udata, 2, jx, jy);
      const sunrealtype qq1 = Q1 * c1 * C3, qq2 = Q2 * c1 * c2, qq3 = q3 * C3,
                        qq4    = q4coef * c2;
      d->at(dudata, 1, jx, jy) = -qq1 - qq2 + TWO * qq3 + qq4;
      d->at(dudata, 2, jx, jy) = qq1 - qq2 - qq4;
    }
  }
  return 0;
}

int diurnal_transport(sunrealtype, N_Vector u, N_Vector udot, void* user_data)
{
  DiurnalData* d         = static_cast<DiurnalData*>(user_data);
  sunrealtype* udata     = N_VGetArrayPointer(u);
  sunrealtype* dudata    = N_VGetArrayPointer(udot);
  const sunrealtype dely = d->dy, verdco = d->vdco, hordco = d->hdco,
                    horaco = d->haco;
  for (int jy = 0; jy < d->my; jy++)
  {
    const sunrealtype ydn  = YMIN + (jy - SUN_RCONST(0.5)) * dely;
    const sunrealtype yup  = ydn + dely;
    const sunrealtype cydn = verdco * std::exp(SUN_RCONST(0.2) * ydn);
    const sunrealtype cyup = verdco * std::exp(SUN_RCONST(0.2) * yup);
    const int idn          = (jy == 0) ? 1 : -1;
    const int iup          = (jy == d->my - 1) ? -1 : 1;
    for (int jx = 0; jx < d->mx; jx++)
    {
      const sunrealtype c1     = d->at(udata, 1, jx, jy);
      const sunrealtype c2     = d->at(udata, 2, jx, jy);
      const sunrealtype c1dn   = d->at(udata, 1, jx, jy + idn);
      const sunrealtype c2dn   = d->at(udata, 2, jx, jy + idn);
      const sunrealtype c1up   = d->at(udata, 1, jx, jy + iup);
      const sunrealtype c2up   = d->at(udata, 2, jx, jy + iup);
      const sunrealtype vertd1 = cyup * (c1up - c1) - cydn * (c1 - c1dn);
      const sunrealtype vertd2 = cyup * (c2up - c2) - cydn * (c2 - c2dn);
      const int ileft          = (jx == 0) ? 1 : -1;
      const int iright         = (jx == d->mx - 1) ? -1 : 1;
      const sunrealtype c1lt   = d->at(udata, 1, jx + ileft, jy);
      const sunrealtype c2lt   = d->at(udata, 2, jx + ileft, jy);
      const sunrealtype c1rt   = d->at(udata, 1, jx + iright, jy);
      const sunrealtype c2rt   = d->at(udata, 2, jx + iright, jy);
      const sunrealtype hord1  = hordco * (c1rt - TWO * c1 + c1lt);
      const sunrealtype hord2  = hordco * (c2rt - TWO * c2 + c2lt);
      const sunrealtype horad1 = horaco * (c1rt - c1lt);
      const sunrealtype horad2 = horaco * (c2rt - c2lt);
      d->at(dudata, 1, jx, jy) = vertd1 + hord1 + horad1;
      d->at(dudata, 2, jx, jy) = vertd2 + hord2 + horad2;
    }
  }
  return 0;
}

int diurnal_jtv_kinetics(N_Vector v, N_Vector Jv, sunrealtype t, N_Vector u,
                         N_Vector, void* user_data, N_Vector)
{
  DiurnalData* d           = static_cast<DiurnalData*>(user_data);
  sunrealtype* udata       = N_VGetArrayPointer(u);
  sunrealtype* vdata       = N_VGetArrayPointer(v);
  sunrealtype* Jvdata      = N_VGetArrayPointer(Jv);
  const sunrealtype s      = std::sin(d->om * t);
  d->q4                    = (s > ZERO) ? std::exp(-A4 / s) : ZERO;
  const sunrealtype q4coef = d->q4;
  for (int jy = 0; jy < d->my; jy++)
  {
    for (int jx = 0; jx < d->mx; jx++)
    {
      const sunrealtype c1     = d->at(udata, 1, jx, jy);
      const sunrealtype c2     = d->at(udata, 2, jx, jy);
      const sunrealtype v1     = d->at(vdata, 1, jx, jy);
      const sunrealtype v2     = d->at(vdata, 2, jx, jy);
      d->at(Jvdata, 1, jx, jy) = -(Q1 * C3 + Q2 * c2) * v1 +
                                 (q4coef - Q2 * c1) * v2;
      d->at(Jvdata, 2, jx, jy) = (Q1 * C3 - Q2 * c2) * v1 -
                                 (q4coef + Q2 * c1) * v2;
    }
  }
  return 0;
}

// P = I - gamma J_kinetics, block by block: the exact inverse of the
// implicit part's Newton matrix (q4 from t, not from the last f call)
int diurnal_psetup_kinetics(sunrealtype t, N_Vector u, N_Vector,
                            sunbooleantype jok, sunbooleantype* jcurPtr,
                            sunrealtype gamma, void* user_data)
{
  DiurnalData* d           = static_cast<DiurnalData*>(user_data);
  sunrealtype* udata       = N_VGetArrayPointer(u);
  const sunrealtype s      = std::sin(d->om * t);
  const sunrealtype q4coef = (s > ZERO) ? std::exp(-A4 / s) : ZERO;
  for (int jy = 0; jy < d->my; jy++)
  {
    for (int jx = 0; jx < d->mx; jx++)
    {
      sunrealtype** j = d->Jbd[(size_t)d->blk(jx, jy)];
      sunrealtype** a = d->P[(size_t)d->blk(jx, jy)];
      if (!jok)
      {
        const sunrealtype c1 = d->at(udata, 1, jx, jy);
        const sunrealtype c2 = d->at(udata, 2, jx, jy);
        j[0][0]              = -Q1 * C3 - Q2 * c2;
        j[1][0]              = -Q2 * c1 + q4coef;
        j[0][1]              = Q1 * C3 - Q2 * c2;
        j[1][1]              = -Q2 * c1 - q4coef;
      }
      SUNDlsMat_denseCopy(j, a, NUM_SPECIES, NUM_SPECIES);
    }
  }
  *jcurPtr = jok ? SUNFALSE : SUNTRUE;
  for (size_t k = 0; k < d->P.size(); k++)
  {
    SUNDlsMat_denseScale(-gamma, d->P[k], NUM_SPECIES, NUM_SPECIES);
    SUNDlsMat_denseAddIdentity(d->P[k], NUM_SPECIES);
    if (SUNDlsMat_denseGETRF(d->P[k], NUM_SPECIES, NUM_SPECIES, d->pivot[k]) != 0)
    {
      return 1;
    }
  }
  return 0;
}

class Diurnal : public Problem
{
public:
  Diurnal()
  {
    spec_.name         = "diurnal";
    spec_.t0           = ZERO;
    spec_.tf           = SUN_RCONST(86400.0);
    spec_.rtol_default = SUN_RCONST(1.0e-5);
    spec_.atol_default = SUN_RCONST(1.0e-3);
    spec_.ls_default   = LinSolKind::SPGMR;
    spec_.has_jtimes   = true;
    spec_.has_prec     = true;
    set_grid(10, 10);
  }

  ~Diurnal() override { data_.release(); }

  void configure(const Options& o) override
  {
    const int mx = o.mx > 0 ? o.mx : 10;
    const int my = o.my > 0 ? o.my : mx;
    set_grid(mx, my);
  }

  void init(N_Vector u) const override
  {
    sunrealtype* udata = N_VGetArrayPointer(u);
    for (int jy = 0; jy < data_.my; jy++)
    {
      const sunrealtype yv = YMIN + jy * data_.dy;
      sunrealtype cy       = (SUN_RCONST(0.1) * (yv - YMID)) *
                       (SUN_RCONST(0.1) * (yv - YMID));
      cy = ONE - cy + SUN_RCONST(0.5) * cy * cy;
      for (int jx = 0; jx < data_.mx; jx++)
      {
        const sunrealtype xv = XMIN + jx * data_.dx;
        sunrealtype cx       = (SUN_RCONST(0.1) * (xv - XMID)) *
                         (SUN_RCONST(0.1) * (xv - XMID));
        cx                         = ONE - cx + SUN_RCONST(0.5) * cx * cx;
        data_.at(udata, 1, jx, jy) = C1_SCALE * cx * cy;
        data_.at(udata, 2, jx, jy) = C2_SCALE * cx * cy;
      }
    }
  }

  void* user_data() override { return &data_; }

  FIRKRhsFn rhs() const override { return diurnal_rhs; }

  FIRKLsJacTimesVecFn jtimes() const override { return diurnal_jtv; }

  FIRKLsPrecSetupFn psetup() const override { return diurnal_psetup; }

  FIRKLsPrecSolveFn psolve() const override { return diurnal_psolve; }

  std::vector<RhsSplit> splits() const override
  {
    RhsSplit imex;
    imex.name     = "imex";
    imex.fe       = diurnal_transport;
    imex.fi       = diurnal_kinetics;
    imex.jtimes_i = diurnal_jtv_kinetics;
    imex.psetup_i = diurnal_psetup_kinetics;
    imex.psolve_i = diurnal_psolve;
    return {imex};
  }

private:
  void set_grid(int mx, int my)
  {
    data_.allocate(mx, my);
    spec_.N = (sunindextype)NUM_SPECIES * mx * my;
  }

  DiurnalData data_;
};

} // namespace

std::unique_ptr<Problem> make_diurnal()
{
  return std::unique_ptr<Problem>(new Diurnal());
}

} // namespace benchmark
