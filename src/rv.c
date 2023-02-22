#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void getbounds(double bounds[], double E_tab[], const double e) {
  // Credit https://github.com/t-brandt/orvara

  const double pi = 3.14159265358979323846264338327950288;
  const double pi_d_12 = 0.26179938779914943653855361527329190667;
  const double pi_d_6 = 0.52359877559829887307710723054658381333;
  const double pi_d_4 = 0.78539816339744830961566084581987572;
  const double pi_d_3 = 1.0471975511965977461542144610931676267;
  const double fivepi_d_12 = 1.3089969389957471826927680763664595333;
  const double pi_d_2 = 1.57079632679489661923132169163975144;
  const double sevenpi_d_12 = 1.8325957145940460557698753069130433467;
  const double twopi_d_3 = 2.0943951023931954923084289221863352533;
  const double threepi_d_4 = 2.35619449019234492884698253745962716;
  const double fivepi_d_6 = 2.6179938779914943653855361527329190667;
  const double elevenpi_d_12 = 2.8797932657906438019240897680062109733;
  // const double pi_d_12 = 3.14159265358979323846264338327950288 / 12;
  // const double pi_d_6 = 3.14159265358979323846264338327950288 / 6;
  // const double pi_d_4 = 3.14159265358979323846264338327950288 / 4;
  // const double pi_d_3 = 3.14159265358979323846264338327950288 / 3;
  // const double fivepi_d_12 = 3.14159265358979323846264338327950288 * 5. / 12;
  // const double pi_d_2 = 3.14159265358979323846264338327950288 / 2;
  // const double sevenpi_d_12 = 3.14159265358979323846264338327950288 * 7. /
  // 12; const double twopi_d_3 = 3.14159265358979323846264338327950288 * 2. /
  // 3; const double threepi_d_4 = 3.14159265358979323846264338327950288 * 3. /
  // 4; const double fivepi_d_6 = 3.14159265358979323846264338327950288 * 5. /
  // 6; const double elevenpi_d_12 = 3.14159265358979323846264338327950288 * 11.
  // / 12;

  const double g2s_e = 0.2588190451025207623489 * e;
  const double g3s_e = 0.5 * e;
  const double g4s_e = 0.7071067811865475244008 * e;
  const double g5s_e = 0.8660254037844386467637 * e;
  const double g6s_e = 0.9659258262890682867497 * e;
  const double g2c_e = g6s_e;
  const double g3c_e = g5s_e;
  const double g4c_e = g4s_e;
  const double g5c_e = g3s_e;
  const double g6c_e = g2s_e;

  bounds[0] = 0;
  bounds[1] = pi_d_12 - g2s_e;
  bounds[2] = pi_d_6 - g3s_e;
  bounds[3] = pi_d_4 - g4s_e;
  bounds[4] = pi_d_3 - g5s_e;
  bounds[5] = fivepi_d_12 - g6s_e;
  bounds[6] = pi_d_2 - e;
  bounds[7] = sevenpi_d_12 - g6s_e;
  bounds[8] = twopi_d_3 - g5s_e;
  bounds[9] = threepi_d_4 - g4s_e;
  bounds[10] = fivepi_d_6 - g3s_e;
  bounds[11] = elevenpi_d_12 - g2s_e;
  bounds[12] = pi;

  double x;

  E_tab[1] = 1. / (1. - e);
  E_tab[2] = 0.;

  x = 1. / (1 - g2c_e);
  E_tab[7] = x;
  E_tab[8] = -0.5 * g2s_e * x * x * x;
  x = 1. / (1 - g3c_e);
  E_tab[13] = x;
  E_tab[14] = -0.5 * g3s_e * x * x * x;
  x = 1. / (1. - g4c_e);
  E_tab[19] = x;
  E_tab[20] = -0.5 * g4s_e * x * x * x;
  x = 1. / (1. - g5c_e);
  E_tab[25] = x;
  E_tab[26] = -0.5 * g5s_e * x * x * x;
  x = 1. / (1. - g6c_e);
  E_tab[31] = x;
  E_tab[32] = -0.5 * g6s_e * x * x * x;

  E_tab[37] = 1.;
  E_tab[38] = -0.5 * e;

  x = 1. / (1. + g6c_e);
  E_tab[43] = x;
  E_tab[44] = -0.5 * g6s_e * x * x * x;
  x = 1. / (1. + g5c_e);
  E_tab[49] = x;
  E_tab[50] = -0.5 * g5s_e * x * x * x;
  x = 1. / (1. + g4c_e);
  E_tab[55] = x;
  E_tab[56] = -0.5 * g4s_e * x * x * x;
  x = 1. / (1. + g3c_e);
  E_tab[61] = x;
  E_tab[62] = -0.5 * g3s_e * x * x * x;
  x = 1. / (1. + g2c_e);
  E_tab[67] = x;
  E_tab[68] = -0.5 * g2s_e * x * x * x;

  E_tab[73] = 1. / (1 + e);
  E_tab[74] = 0;

  double B0, B1, B2, idx, idx2;
  int k;
  for (int i = 0; i < 12; i++) {
    idx = 1. / (bounds[i + 1] - bounds[i]);
    idx2 = idx * idx;
    k = 6. * i;
    E_tab[k] = i * pi_d_12;

    B0 = idx * -E_tab[k + 2] - idx2 * (E_tab[k + 1] - idx * pi_d_12);
    // B0 = idx * (-E_tab[k + 2] - idx * (E_tab[k + 1] - idx * pi_d_12));
    B1 = idx * (-2 * E_tab[k + 2]) - idx2 * (E_tab[k + 1] - E_tab[k + 7]);
    // B1 = idx * (-2 * E_tab[k + 2] - idx * (E_tab[k + 1] - E_tab[k + 7]));
    B2 = idx * (E_tab[k + 8] - E_tab[k + 2]);

    E_tab[k + 3] = B2 - 4. * B1 + 10. * B0;
    E_tab[k + 4] = (-2. * B2 + 7. * B1 - 15. * B0) * idx;
    E_tab[k + 5] = (B2 - 3. * B1 + 6. * B0) * idx2;
  }

  return;
}

inline double shortsin(const double x) {
  // Credit: https://github.com/t-brandt/orvara

  // const double if3 = 1. / 6;
  // const double if5 = 1. / (6. * 20);
  // const double if7 = 1. / (6. * 20 * 42);
  // const double if9 = 1. / (6. * 20 * 42 * 72);
  // const double if11 = 1. / (6. * 20 * 42 * 72 * 110);
  // const double if13 = 1. / (6. * 20 * 42 * 72 * 110 * 156);
  // const double if15 = 1. / (6. * 20 * 42 * 72 * 110 * 156 * 210);
  const double if3 = 0.16666666666666666666666666666666666667;
  const double if5 = 0.0083333333333333333333333333333333333333;
  const double if7 = 1.984126984126984126984126984126984127e-4;
  const double if9 = 2.7557319223985890652557319223985890653e-6;
  const double if11 = 2.5052108385441718775052108385441718775e-8;
  const double if13 = 1.6059043836821614599392377170154947933e-10;
  const double if15 = 7.6471637318198164759011319857880704442e-13;

  const double x2 = x * x;

  return x *
         (1 - x2 * (if3 -
                    x2 * (if5 -
                          x2 * (if7 -
                                x2 * (if9 - x2 * (if11 -
                                                  x2 * (if13 - x2 * if15)))))));
}

inline double Estart(const double M, const double e) {
  // Credit: https://github.com/t-brandt/orvara

  const double ome = 1. - e;
  const double sqrt_ome = sqrt(ome);

  const double chi = M / (sqrt_ome * ome);
  const double Lam = sqrt(8 + 9 * chi * chi);
  const double S = cbrt(Lam + 3 * chi);
  const double sigma = 6 * chi / (2 + S * S + 4. / (S * S));
  const double s2 = sigma * sigma;
  const double denom = s2 + 2;
  const double E =
      sigma *
      (1 + s2 * ome *
               ((s2 + 20) / (60. * denom) +
                s2 * ome * (s2 * s2 * s2 + 25 * s2 * s2 + 340 * s2 + 840) /
                    (1400 * denom * denom * denom)));
  return E * sqrt_ome;
}

/*=================================================================
 * eccanom_orvara   Invert Kepler's time equation for elliptical orbits
 *                  using orvara's method.
 *
 * eccanom_orvara(E, sinE, cosE, M, e, n)
 *   E       nx1     Output array - Eccentric anomalies (rad)
 *   sinE    nx1     Output array - Sine of eccentric anomalies (rad)
 *   cosE    nx1     Output array - Cosine of eccentric anomalies (rad)
 *   M       nx1     Mean anomalies (rad)
 *   e       1       Eccentricity
 *   n       1       Lengths of E, sinE, cosE, M
 *=================================================================*/
void eccanom_orvara(double E[], double sinE[], double cosE[], const double M[],
                    const double e, const int n) {
  double E_tab[6 * 13];
  double bounds[13];
  getbounds(bounds, E_tab, e);

  const double pi = 3.14159265358979323846264338327950288;
  const double pi_d_4 = 0.78539816339744830961566084581987572;
  const double pi_d_2 = 1.57079632679489661923132169163975144;
  const double threepi_d_4 = 2.35619449019234492884698253745962716;
  const double twopi = 6.28318530717958647692528676655900576;
  int j, k;
  double dx;

  const double one_over_ecc = 1.0 / fmax(1e-17, e);
  double _M, _E, _sinE, _cosE, dE, one_minus_dE2_d_2, dEsq_d6;
  int Esign;
  double num, denom;

  const double one_sixth = 1. / 6;
  if (e < 0.78) {
    for (int i = 0; i < n; ++i) {
      _M = M[i];

      // Cut mean anomaly between 0 and pi to use shorter Taylor series
      if (_M > pi) {
        Esign = -1;
        _M = twopi - _M;
      } else {
        Esign = 1;
      }

      // Find the relevant interval, searching backwards
      for (j = 11;; --j) {
        if (_M > bounds[j]) {
          break;
        }
      }
      k = 6 * j;
      dx = _M - bounds[j];

      // Initial guess from lookup table
      _E = E_tab[k] +
           dx * (E_tab[k + 1] +
                 dx * (E_tab[k + 2] +
                       dx * (E_tab[k + 3] +
                             dx * (E_tab[k + 4] + dx * E_tab[k + 5]))));

      // Calculate _sinE and _cosE using the short sin function and sqrt call
      // ALMOST ALL COMPUTATION TIME IS SPENT HERE
      if (!(_E > pi_d_4)) {
        _sinE = shortsin(_E);
        _cosE = sqrt(1. - _sinE * _sinE);
      } else if (_E < threepi_d_4) {
        _cosE = shortsin(pi_d_2 - _E);
        _sinE = sqrt(1. - _cosE * _cosE);
      } else {
        _sinE = shortsin(pi - _E);
        _cosE = -sqrt(1. - _sinE * _sinE);
      }
      num = (_M - _E) * one_over_ecc + _sinE;
      denom = one_over_ecc - _cosE;

      // Get the second order approximation of dE
      dE = num * denom / (denom * denom + 0.5 * _sinE * num);
      one_minus_dE2_d_2 = 1 - (dE * dE / 2);

      // Apply correction to E, sinE, and _cosE with second order
      // approximation
      E[i] = fmod(Esign * (_E + dE) + twopi, twopi);
      sinE[i] = Esign * (_sinE * one_minus_dE2_d_2 + dE * _cosE);
      cosE[i] = _cosE * one_minus_dE2_d_2 - dE * _sinE;
    }
  }
  // For higher eccentricities we need to go to third order
  else {
    for (int i = 0; i < n; ++i) {
      _M = M[i];
      // Cut mean anomaly between 0 and pi to use shorter Taylor series
      if (_M > pi) {
        Esign = -1;
        _M = twopi - _M;
      } else {
        Esign = 1;
      }
      if ((2 * _M + (1 - e)) > 0.2) {
        for (j = 11;; --j) {
          if (_M > bounds[j]) {
            break;
          }
        }
        k = 6 * j;
        dx = _M - bounds[j];

        // Initial guess from lookup table
        _E = E_tab[k] +
             dx * (E_tab[k + 1] +
                   dx * (E_tab[k + 2] +
                         dx * (E_tab[k + 3] +
                               dx * (E_tab[k + 4] + dx * E_tab[k + 5]))));

      } else {
        _E = Estart(_M, e);
      }

      // Calculate _sinE and _cosE using the short sin function and sqrt call
      if (!(_E > pi_d_4)) {
        _sinE = shortsin(_E);
        _cosE = sqrt(1. - _sinE * _sinE);
      } else if (_E < threepi_d_4) {
        _cosE = shortsin(pi_d_2 - _E);
        _sinE = sqrt(1 - _cosE * _cosE);
      } else {
        _sinE = shortsin(pi - _E);
        _cosE = -sqrt(1 - _sinE * _sinE);
      }

      num = (_M - _E) * one_over_ecc + _sinE;
      denom = one_over_ecc - _cosE;

      if (_M > 0.4) {
        dE = num * denom / (denom * denom + 0.5 * _sinE * num);
      } else {
        dE = num * (denom * denom + 0.5 * num * _sinE);
        dE /= denom * denom * denom +
              num * (denom * _sinE + one_sixth * num * _cosE);
      }
      dEsq_d6 = dE * dE * one_sixth;

      E[i] = fmod(Esign * (_E + dE) + twopi, twopi);
      sinE[i] =
          Esign * (_sinE * (1 - 3 * dEsq_d6) + dE * (1 - dEsq_d6) * _cosE);
      cosE[i] = _cosE * (1 - 3 * dEsq_d6) - dE * (1 - dEsq_d6) * _sinE;
    }
  }
}
double timetrans_to_timeperi_c(const double tc, const double per,
                               const double ecc, const double omega) {
  const double PI = 3.141592653589793;
  double f, ee, tp;
  f = PI / 2.0 - omega;
  ee = 2.0 * atan(tan(f / 2.0) * sqrt((1.0 - ecc) / (1.0 + ecc)));
  tp = tc - per / (2.0 * PI) * (ee - ecc * sin(ee));
  return tp;
}

void meananom(double M[], const double t[], const double tp, const double per,
              const int n) {
  const double twopi = 6.283185307179586;
  double phase;
  for (int i = 0; i < n; ++i) {
    phase = (t[i] - tp) / per;
    M[i] = twopi * (phase - floorf(phase));
  }
  return;
}

void RV_from_time(double rv[], const double t[], const double tp[],
                  const double per[], const double e[], const double w[],
                  const double K[], const int ntimes, const int nplan,
                  const int specific_planet) {
  /*Finds radial velocity for a single object at the desired epochs

    Args:
        rv (ndarray):
            Preexisting radial velocities, can also be zeros (rad)
        t (ndarray):
            Times of to calculate RV at (jd)
        tp (float):
            Time of periastron
        per (float):
            Period
        e (float):
            Eccentricity
        w (float):
            Argument of periapsis (rad)
        K (float):
            RV semi-amplitude (m/s)
  */

  const double pi = 3.14159265358979323846264338327950288;
  const double twopi = 6.28318530717958647692528676655900576;
  const double pi_d_2 = 1.57079632679489661923132169163975144;
  const double one_d_24 = 0.041666666666666666666666666666666666667; // 1. / 24;
  const double one_d_240 =
      0.0041666666666666666666666666666666666667; // 1. / 240;

  int p0, pf;
  if (!specific_planet) {
    p0 = 0;
    pf = nplan;
  } else {
    p0 = specific_planet - 1;
    pf = specific_planet;
  }
  double _tp, _per, _e, _w, _K, _mean_motion;
  double M[ntimes], E[ntimes], sinE[ntimes], cosE[ntimes];
  double sqrt1pe, sqrt1me, cosarg, sinarg, ecccosarg, sqrt1pe_div_sqrt1me;
  double TA, ratio, fac, tanEAd2;

  double _E;
  for (int j = p0; j < pf; ++j) {
    _tp = tp[j];
    _per = per[j];
    _e = e[j];
    _w = w[j];
    _K = K[j];
    if (_K == 0) {
      // Trivial case that can be ignored
      continue;
    } else {

      // Calculate mean anomaly
      meananom(M, t, _tp, _per, ntimes);
      cosarg = cos(_w);
      sinarg = sin(_w);

      // Calculating E, sinE, and cosE from M
      eccanom_orvara(E, sinE, cosE, M, _e, ntimes);

      sqrt1pe = sqrt(1.0 + _e);
      sqrt1me = sqrt(1.0 - _e);

      cosarg = cos(_w);
      sinarg = sin(_w);
      ecccosarg = _e * cosarg;
      sqrt1pe_div_sqrt1me = sqrt1pe / sqrt1me;

      // ##################################################################
      // #Trickery with trig identities.The code below is mathematically
      // #identical to the use of the true anomaly.If sin(EA) is small
      // #and cos(EA) is close to - 1, no problem as long as sin(EA) is not
      // #precisely zero(set tan(EA / 2) = 1e100 in this case).If sin(EA)
      // #is small and EA is close to zero, use the fifth - order Taylor
      // #expansion for tangent.This is good to ~1e-15 for EA within
      // #~0.015 of 0. Assume eccentricity is not precisely unity(this
      // #should be forbidden by the priors).Very, very high
      // #eccentricities(significantly above 0.9999) may be problematic.
      // #This routine assumes range reduction of the eccentric anomaly to
      // #(- pi, pi] and will throw an error if this is violated.
      // ##################################################################

      for (int i = 0; i < ntimes; ++i) {
        _E = E[i];
        if (_E > pi) {
          _E = twopi - _E;
        }
        if (fabs(sinE[i]) > 1.5e-2) {
          tanEAd2 = (1.0 - cosE[i]) / sinE[i];
        } else if (fabs(_E) < pi_d_2) {
          tanEAd2 = _E * (0.5 + _E * _E * (one_d_24 + one_d_240 * _E * _E));
        } else if (sinE[i] != 0) {
          tanEAd2 = (1.0 - cosE[i]) / sinE[i];
        } else {
          tanEAd2 = 1e100;
        }
        ratio = sqrt1pe_div_sqrt1me * tanEAd2;
        fac = 2.0 / (1.0 + ratio * ratio);
        rv[i] += _K * (cosarg * (fac - 1.0) - sinarg * ratio * fac + ecccosarg);
      }
    }
  }
  return;
}

void convert_basis(const double params[], const char basis[], const int nplan,
                   double per[], double tp[], double e[], double w[],
                   double k[]) {
  /*
  Working with the radvel basis set directly here. Standardize to "per tp e w k"
  */
  if (strcmp(basis, "per tp e w k") == 0) {
    for (int i = 0; i < nplan; ++i) {
      // param array is (nplan*5+n_derived_params)x4, index w/ one value for
      // speed
      const int per_pos = 5 * i * 4;
      const int tp_pos = per_pos + 4;
      const int e_pos = tp_pos + 4;
      const int w_pos = e_pos + 4;
      const int k_pos = w_pos + 4;
      per[i] = params[per_pos];
      tp[i] = params[tp_pos];
      e[i] = params[e_pos];
      w[i] = params[w_pos];
      k[i] = params[k_pos];
    }
    return;
  } else if (strcmp(basis, "per tc secosw sesinw k") == 0) {
    for (int i = 0; i < nplan; ++i) {
      // param array is (nplan*5+n_derived_params)x4, index w/ one value for
      // speed
      const int per_pos = 5 * i * 4;
      const int tc_pos = per_pos + 4;
      const int secosw_pos = tc_pos + 4;
      const int sesinw_pos = secosw_pos + 4;
      const int k_pos = sesinw_pos + 4;
      const double _tc = params[tc_pos];
      const double _secosw = params[secosw_pos];
      const double _sesinw = params[sesinw_pos];
      per[i] = params[per_pos];
      e[i] = _secosw * _secosw + _sesinw * _sesinw;
      w[i] = atan2(_sesinw, _secosw);
      tp[i] = timetrans_to_timeperi_c(_tc, params[per_pos], e[i], w[i]);
      k[i] = params[k_pos];
    }
    return;
  } else if (strcmp(basis, "logper tc secosw sesinw k") == 0) {
    for (int i = 0; i < nplan; ++i) {
      const int logper_pos = 5 * i * 4;
      const int tc_pos = logper_pos + 4;
      const int secosw_pos = tc_pos + 4;
      const int sesinw_pos = secosw_pos + 4;
      const int k_pos = sesinw_pos + 4;
      const double _tc = params[tc_pos];
      const double _secosw = params[secosw_pos];
      const double _sesinw = params[sesinw_pos];
      per[i] = exp(params[logper_pos]);
      e[i] = _secosw * _secosw + _sesinw * _sesinw;
      w[i] = atan2(_sesinw, _secosw);
      tp[i] = timetrans_to_timeperi_c(_tc, per[i], e[i], w[i]);
      k[i] = params[k_pos];
    }
    return;
  } else {
    printf("Basis %s not found, decent chance it isn't programmed yet.", basis);
    exit(0);
  }
}

void rv_calc(double rv[], const double t[], const double params[],
             const char basis[], const int ntimes, const int nplan,
             const int specific_planet) {
  /*
  Calculating radial velocities for a system of n planets at n times
  t (nx1)
  params - Radvel parameter array to be used
  basis - string with the basis set (e.g. 'per tc se w k') to be converted to
  'per tp e w k'
  */
  double per[nplan], tp[nplan], e[nplan], w[nplan], k[nplan];
  convert_basis(params, basis, nplan, per, tp, e, w, k);

  RV_from_time(rv, t, tp, per, e, w, k, ntimes, nplan, specific_planet);
}

void model_call(double mod[], const double t[], const double params[],
                const char basis[], const double time_base, const int jit_ind,
                const int gamma_ind, const int ntimes, const int nplan,
                const int specific_planet) {
  const double pi = 3.14159265358979323846264338327950288;
  const double twopi = 2. * pi;

  const int dvdt_ind = 5 * nplan * 4;
  const double dvdt = params[dvdt_ind];
  const double curv = params[dvdt_ind + 4];

  // Calculate rv values with rv_calc and add dvdt+curv to match the
  // RVLikelihood __call__ method
  rv_calc(mod, t, params, basis, ntimes, nplan, specific_planet);
  double tdiff;
  for (int i = 0; i < ntimes; ++i) {
    // for (int i = ntimes; i--;) {
    tdiff = t[i] - time_base;
    mod[i] += dvdt * tdiff + curv * tdiff * tdiff;
  }
  return;
}

void residuals(double res[], const double t[], const double rv[],
               const double rv_err[], const double params[], const char basis[],
               const double time_base, const int jit_ind, const int gamma_ind,
               const int ntimes, const int nplan) {
  const double twopi = 6.28318530717958647692528676655900576;

  /***********************
  | RESIDUAL CALCULATION *
  ***********************/
  // Replaces line: mod = self.model(t)
  double mod[ntimes]; // Model's RV values
  for (int i = 0; i < ntimes; ++i) {
    // Memory problems if not initialized this way
    mod[i] = 0.0;
  }

  const double jit = params[jit_ind * 4];
  const double jit2 = jit * jit;
  const int true_gamma_ind = gamma_ind * 4;
  double gamma_val = params[true_gamma_ind];

  model_call(mod, t, params, basis, time_base, jit_ind, gamma_ind, ntimes,
             nplan, 0);
  // const int dvdt_ind = 5 * nplan * 4;
  // const double dvdt = params[dvdt_ind];
  // const double curv = params[dvdt_ind + 4];

  // // Calculate rv values with rv_calc and add dvdt+curv to match the
  // // RVLikelihood __call__ method
  // rv_calc(mod, t, params, basis, ntimes, nplan, 0);
  // double tdiff;
  // for (int i = 0; i < ntimes; ++i) {
  //   tdiff = t[i] - time_base;
  //   mod[i] += dvdt * tdiff + curv * tdiff * tdiff;
  // }

  // Replacing the if statement for the gamma value
  double rv_err2[ntimes];
  double sum_sig_quad[ntimes];
  double rv_mod_diff[ntimes];
  if ((params[true_gamma_ind + 3]) && !(params[true_gamma_ind + 1])) {
    double sum_sq;
    double ztil_num = 0.0;
    double summed_sum_sq = 0.0;
    for (int i = 0; i < ntimes; ++i) {
      rv_err2[i] = rv_err[i] * rv_err[i];
      sum_sig_quad[i] = rv_err2[i] + jit2;
      rv_mod_diff[i] = rv[i] - mod[i];

      sum_sq = 1.0 / (sum_sig_quad[i]);
      summed_sum_sq += sum_sq;
      ztil_num += rv_mod_diff[i] * sum_sq;
    }
    gamma_val = ztil_num / summed_sum_sq;
  } else {
    for (int i = 0; i < ntimes; ++i) {
      rv_mod_diff[i] = rv[i] - mod[i];
    }
  }
  for (int i = 0; i < ntimes; ++i) {
    res[i] = rv_mod_diff[i] - gamma_val;
  }
  return;
}

double logprob(const double t[], const double rv[], const double rv_err[],
               const double params[], const char basis[],
               const double time_base, const int jit_ind, const int gamma_ind,
               const int ntimes, const int nplan) {
  double _logprob = 0.0;
  const double twopi = 6.28318530717958647692528676655900576;

  /***********************
  | RESIDUAL CALCULATION *
  ***********************/
  double res[ntimes]; // Residuals
  residuals(res, t, rv, rv_err, params, basis, time_base, jit_ind, gamma_ind,
            ntimes, nplan);
  // double mod[ntimes]; // Model's RV values

  // for (int i = 0; i < ntimes; ++i) {
  //   // Memory problems if not initialized this way
  //   mod[i] = 0.0;
  // }

  const double jit = params[jit_ind * 4];
  const double jit2 = jit * jit;
  const int true_gamma_ind = gamma_ind * 4;
  double gamma_val = params[true_gamma_ind];

  // const int dvdt_ind = 5 * nplan * 4;
  // const double dvdt = params[dvdt_ind];
  // const double curv = params[dvdt_ind + 4];

  // // Calculate rv values with rv_calc and add dvdt+curv to match the
  // // RVLikelihood __call__ method
  // // Replaces line: mod = self.model(t)
  // rv_calc(mod, t, params, basis, ntimes, nplan, 0);
  // double tdiff;
  // for (int i = 0; i < ntimes; ++i) {
  //   tdiff = t[i] - time_base;
  //   mod[i] += dvdt * tdiff + curv * tdiff * tdiff;
  // }

  // // Replacing the if statement for the gamma value
  // int rv_err2_calculated = 0;
  // // double rv_err2[ntimes];
  // if ((params[true_gamma_ind + 3]) && !(params[true_gamma_ind + 1])) {
  //   int rv_err2_calculated = 1;
  //   double ztil, ztil_num;
  //   double sum_sq;
  //   double summed_sum_sq = 0.0;
  //   for (int i = 0; i < ntimes; ++i) {
  //     rv_err2[i] = rv_err[i] * rv_err[i];
  //     sum_sq = 1.0 / (rv_err2[i] + jit2);
  //     summed_sum_sq += sum_sq;
  //     ztil_num += (rv[i] - mod[i]) * sum_sq;
  //   }
  //   gamma_val = ztil_num / summed_sum_sq;
  // }
  // for (int i = 0; i < ntimes; ++i) {
  //   res[i] = rv[i] - gamma_val - mod[i];
  // }

  // WHERE THE DECORRELATION STUFF WOULD GO

  /**********************
  | LOGLIKE CALCULATION *
  **********************/
  // Call is loglike_jitter(residuals, rv_err, jit)
  // double res2;
  // double rv_err2[ntimes];
  double sum_sig_quad[ntimes];
  double penalty = 0;
  double chi2 = 0;
  for (int i = 0; i < ntimes; ++i) {
    // rv_err2[i] = rv_err[i] * rv_err[i];
    sum_sig_quad[i] = rv_err[i] * rv_err[i] + jit2;
    penalty += log(sqrt(twopi * sum_sig_quad[i]));
    // Chi2
    // res2 = res[i] * res[i];
    chi2 += res[i] * res[i] / sum_sig_quad[i];
  }
  _logprob = -0.5 * chi2 - penalty;

  /*****************
  | END OF LOGPROB *
  *****************/
  double sigz, sigz_denom = 0.0;
  if ((params[true_gamma_ind + 3]) && !(params[true_gamma_ind + 1])) {
    for (int i = 0; i < ntimes; ++i) {
      sigz_denom += 1.0 / sum_sig_quad[i];
    }
    sigz = 1.0 / sigz_denom;
    _logprob += log(sqrt(twopi * sigz));
  }
  return _logprob;
}

// https://github.com/scipy/scipy/blob/v1.10.1/scipy/optimize/_minimize.py#L687
// res = _minimize_powell(fun, x0, args, callback, bounds, **options)
// https://github.com/scipy/scipy/blob/v1.10.1/scipy/optimize/_optimize.py
// def _minimize_powell(func, x0, args=(), callback=None, bounds=None,
//                      xtol=1e-4, ftol=1e-4, maxiter=None, maxfev=None,
//                      disp=False, direc=None, return_all=False,
//                      **unknown_options):
