#include "../microfacet.h"

static const int pMax = 3;

inline Real SafeSqrt(Real x) {
    return std::sqrt(std::max(Real(0), x));
}

inline Real SafeASin(Real x) {
    return std::asin(std::clamp(x, Real(-1), Real(1)));
}

inline Real Sqr(Real v) { return v * v; }

template <int n>
static Real Pow(Real v) {
    Real n2 = Pow<n / 2>(v);
    return n2 * n2 * Pow<n & 1>(v);
}
template <> Real Pow<1>(Real v) { return v; }
template <> Real Pow<0>(Real v) { return 1; }

inline Real I0(Real x) {
    Real val = 0;
    Real x2i = 1;
    int64_t ifact = 1;
    int i4 = 1;
    // I0(x) \approx Sum_i x^(2i) / (4^i (i!)^2)
    for (int i = 0; i < 10; ++i) {
        if (i > 1) ifact *= i;
        val += x2i / (i4 * Sqr(ifact));
        x2i *= x * x;
        i4 *= 4;
    }
    return val;
}

inline Real LogI0(Real x) {
    if (x > Real(12))
        return x + Real(0.5) * (-std::log(2 * c_PI) + std::log(1 / x) + 1 / (8 * x));
    else
        return std::log(I0(x));
}

inline Real Phi(int p, Real gammaO, Real gammaT) {
    return 2 * p * gammaT - 2 * gammaO + p * c_PI;
}

inline Real Logistic(Real x, Real s) {
    x = std::abs(x);
    return std::exp(-x / s) / (s * Sqr(1 + std::exp(-x / s)));
}

inline Real LogisticCDF(Real x, Real s) {
    return 1 / (1 + std::exp(-x / s));
}

inline Real TrimmedLogistic(Real x, Real s, Real a, Real b) {
    return Logistic(x, s) / (LogisticCDF(b, s) - LogisticCDF(a, s));
}

inline Real Np(Real phi, int p, Real s, Real gammaO, Real gammaT) {
    Real dphi = phi - Phi(p, gammaO, gammaT);
    // Remap _dphi_ to $[-\pi,\pi]$
    while (dphi > c_PI) dphi -= 2 * c_PI;
    while (dphi < -c_PI) dphi += 2 * c_PI;
    return TrimmedLogistic(dphi, s, -c_PI, c_PI);
}

static Real SampleTrimmedLogistic(Real u, Real s, Real a, Real b) {
    Real k = LogisticCDF(b, s) - LogisticCDF(a, s);
    Real x = -s * std::log(1 / (u * k + LogisticCDF(a, s)) - 1);
    return std::clamp(x, a, b);
}

// Longitudinal Scattering Functions
static Real Mp(Real cosThetaI, Real cosThetaO, Real sinThetaI,
               Real sinThetaO, Real v) {
    Real a = cosThetaI * cosThetaO / v;
    Real b = sinThetaI * sinThetaO / v;
    Real mp =
        (v <= Real(0.1))
            ? (std::exp(LogI0(a) - b - 1 / v + Real(0.6931) + std::log(1 / (2 * v))))
            : (std::exp(-b) * I0(a)) / (std::sinh(1 / v) * 2 * v);
    return mp;
}

Real FrDielectric(Real cosThetaI, Real etaI, Real etaT) {
    cosThetaI = std::clamp(cosThetaI, Real(-1), Real(1));
    // Potentially swap indices of refraction
    bool entering = cosThetaI > Real(0);
    if (!entering) {
        std::swap(etaI, etaT);
        cosThetaI = std::abs(cosThetaI);
    }

    // Compute _cosThetaT_ using Snell's law
    Real sinThetaI = std::sqrt(std::max(Real(0), 1 - cosThetaI * cosThetaI));
    Real sinThetaT = etaI / etaT * sinThetaI;

    // Handle total internal reflection
    if (sinThetaT >= 1) return 1;
    Real cosThetaT = std::sqrt(std::max(Real(0), 1 - sinThetaT * sinThetaT));
    Real Rparl = ((etaT * cosThetaI) - (etaI * cosThetaT)) /
                  ((etaT * cosThetaI) + (etaI * cosThetaT));
    Real Rperp = ((etaI * cosThetaI) - (etaT * cosThetaT)) /
                  ((etaI * cosThetaI) + (etaT * cosThetaT));
    return (Rparl * Rparl + Rperp * Rperp) / 2;
}

static std::array<Spectrum, pMax + 1> Ap(Real cosThetaO, Real eta, Real h,
                                         const Spectrum &T) {
    std::array<Spectrum, pMax + 1> ap;
    // Compute $p=0$ attenuation at initial cylinder intersection
    Real cosGammaO = SafeSqrt(1 - h * h);
    Real cosTheta = cosThetaO * cosGammaO;
    Real f = FrDielectric(cosTheta, Real(1), eta);
    ap[0] = make_const_spectrum(f);

    // Compute $p=1$ attenuation term
    ap[1] = Sqr(Real(1) - f) * T;

    // Compute attenuation terms up to $p=_pMax_$
    for (int p = 2; p < pMax; ++p) ap[p] = ap[p - 1] * T * f;

    // Compute attenuation term accounting for remaining orders of scattering
    ap[pMax] = ap[pMax - 1] * f * T / (make_const_spectrum(1) - T * f);
    return ap;
}

std::array<Real, pMax + 1> ComputeApPdf(Real cosThetaO, Real eta, Real h, Spectrum sigma_a) {
    // Compute array of $A_p$ values for _cosThetaO_
    Real sinThetaO = SafeSqrt(1 - cosThetaO * cosThetaO);

    // Compute $\cos \thetat$ for refracted ray
    Real sinThetaT = sinThetaO / eta;
    Real cosThetaT = SafeSqrt(1 - Sqr(sinThetaT));

    // Compute $\gammat$ for refracted ray
    Real etap = std::sqrt(eta * eta - Sqr(sinThetaO)) / cosThetaO;
    Real sinGammaT = h / etap;
    Real cosGammaT = SafeSqrt(1 - Sqr(sinGammaT));

    // Compute the transmittance _T_ of a single path through the cylinder
    Spectrum T = exp(-sigma_a * (2 * cosGammaT / cosThetaT));
    std::array<Spectrum, pMax + 1> ap = Ap(cosThetaO, eta, h, T);

    // Compute $A_p$ PDF from individual $A_p$ terms
    std::array<Real, pMax + 1> apPdf;

    Real sumY = Real(0);
    for ( int i = 0; i <= pMax; i++ ) {
        sumY += luminance(ap[i]);
    }
    // Real sumY =
    //     std::accumulate(ap.begin(), ap.end(), Real(0),
    //                     [](Real s, const Spectrum &ap) { return s + ap[1]; });

    for (int i = 0; i <= pMax; ++i) apPdf[i] = luminance(ap[i]) / sumY;
    return apPdf;
}

// https://fgiesen.wordpress.com/2009/12/13/decoding-morton-codes/
static uint32_t Compact1By1(uint32_t x) {
    // TODO: as of Haswell, the PEXT instruction could do all this in a
    // single instruction.
    // x = -f-e -d-c -b-a -9-8 -7-6 -5-4 -3-2 -1-0
    x &= 0x55555555;
    // x = --fe --dc --ba --98 --76 --54 --32 --10
    x = (x ^ (x >> 1)) & 0x33333333;
    // x = ---- fedc ---- ba98 ---- 7654 ---- 3210
    x = (x ^ (x >> 2)) & 0x0f0f0f0f;
    // x = ---- ---- fedc ba98 ---- ---- 7654 3210
    x = (x ^ (x >> 4)) & 0x00ff00ff;
    // x = ---- ---- ---- ---- fedc ba98 7654 3210
    x = (x ^ (x >> 8)) & 0x0000ffff;
    return x;
}

static Vector2 DemuxFloat(float f) {
    uint64_t v = f * (1ull << 32);
    uint32_t bits[2] = {Compact1By1(v), Compact1By1(v >> 1)};
    return {bits[0] / float(1 << 16), bits[1] / float(1 << 16)};
}

Spectrum eval_op::operator()(const HairBSDF &bsdf) const {
    bool reflect = dot(vertex.geometry_normal, dir_in) *
                   dot(vertex.geometry_normal, dir_out) > 0;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometry_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    (void)reflect; // silence unuse warning, remove this when implementing hw


    // Get the parameterized positions 
    Real u = vertex.st.x; // u=[0,1] along the curve
    Real v = Real(0.5); //vertex.st.y; // v=[0,1] along the width of the tube, should be 0 for RTC_GEOMETRY_TYPE_ROUND_BEZIER_CURVE but not
    Real h = Real(-1) + Real(2) * v; // h=[-1, 1] point on the diameter of the intersection circle
    Real gammaO = SafeASin(h);

    // Compute hair coordinate system terms related to dir_out
    Real sinThetaO = dir_in.x;
    Real cosThetaO = SafeSqrt(1 - Sqr(sinThetaO));
    Real phiO = std::atan2(dir_in.z, dir_in.y);

    // Compute hair coordinate system terms related to dir_in
    Real sinThetaI = dir_out.x;
    Real cosThetaI = SafeSqrt(1 - Sqr(sinThetaI));
    Real phiI = std::atan2(dir_out.z, dir_out.y);

    // Precompute the roughness v based on beta_m at each bounce
    Real longitudinal_v[pMax + 1];
    longitudinal_v[0] = Sqr(Real(0.726) * bsdf.beta_m + Real(0.812) * Sqr(bsdf.beta_m) + Real(3.7) * Pow<20>(bsdf.beta_m));
    longitudinal_v[1] = Real(0.25) * longitudinal_v[0];
    longitudinal_v[2] = 4 * longitudinal_v[0];
    for (int p = 3; p <= pMax; ++p) {
        longitudinal_v[p] = longitudinal_v[2];
    }


    // 2. Compute Absorption

    // Compute the direction of the refracted ray using Snell's law
    Real sinThetaT = sinThetaO / bsdf.eta;
    Real cosThetaT = SafeSqrt(1 - Sqr(sinThetaT));

    Real etap = std::sqrt(bsdf.eta * bsdf.eta - Sqr(sinThetaO)) / cosThetaO;
    Real sinGammaT = h / etap;
    Real cosGammaT = SafeSqrt(1 - Sqr(sinGammaT));
    Real gammaT = SafeASin(sinGammaT);

    Spectrum T = exp(-bsdf.sigma_a * (Real(2) * cosGammaT / cosThetaT));


    // 3. Compute Azimuthal Scattering

    Real s = sqrt(c_PI/Real(8)) * (Real(0.265) * bsdf.beta_n + Real(1.194) * Sqr(bsdf.beta_n) + Real(5.372) * Pow<22>(bsdf.beta_n));

    // Compute alpha terms for hair scales
    Real sin2kAlpha[3], cos2kAlpha[3];
    sin2kAlpha[0] = std::sin(radians(bsdf.alpha));
    cos2kAlpha[0] = SafeSqrt(1 - Sqr(sin2kAlpha[0]));
    for (int i = 1; i < 3; ++i) {
        sin2kAlpha[i] = 2 * cos2kAlpha[i - 1] * sin2kAlpha[i - 1];
        cos2kAlpha[i] = Sqr(cos2kAlpha[i - 1]) - Sqr(sin2kAlpha[i - 1]);
    }


    // Evaluate hair BSDF
    Real phi = phiI - phiO;
    std::array<Spectrum, pMax + 1> ap = Ap(cosThetaO, bsdf.eta, h, T);
    Spectrum fsum = make_zero_spectrum();
    for (int p = 0; p < 1; ++p) {
    // for (int p = 0; p < pMax; ++p) {
        // Compute $\sin \thetao$ and $\cos \thetao$ terms accounting for scales
        Real sinThetaOp, cosThetaOp;
        if (p == 0) {
            sinThetaOp = sinThetaO * cos2kAlpha[1] - cosThetaO * sin2kAlpha[1];
            cosThetaOp = cosThetaO * cos2kAlpha[1] + sinThetaO * sin2kAlpha[1];
        }

        // Handle remainder of $p$ values for hair scale tilt
        else if (p == 1) {
            sinThetaOp = sinThetaO * cos2kAlpha[0] + cosThetaO * sin2kAlpha[0];
            cosThetaOp = cosThetaO * cos2kAlpha[0] - sinThetaO * sin2kAlpha[0];
        } else if (p == 2) {
            sinThetaOp = sinThetaO * cos2kAlpha[2] + cosThetaO * sin2kAlpha[2];
            cosThetaOp = cosThetaO * cos2kAlpha[2] - sinThetaO * sin2kAlpha[2];
        } else {
            sinThetaOp = sinThetaO;
            cosThetaOp = cosThetaO;
        }

        // Handle out-of-range $\cos \thetao$ from scale adjustment
        cosThetaOp = std::abs(cosThetaOp);
        fsum += Mp(cosThetaI, cosThetaOp, sinThetaI, sinThetaOp, longitudinal_v[p]) * ap[p] *
                Np(phi, p, s, gammaO, gammaT);
    }

    // Compute contribution of remaining terms after _pMax_
    fsum += Mp(cosThetaI, cosThetaO, sinThetaI, sinThetaO, longitudinal_v[pMax]) * ap[pMax] /
            (Real(2) * c_PI);
    if (std::abs(dir_out.z) > 0) fsum /= std::abs(dir_out.z);

    return fsum;
}

Real pdf_sample_bsdf_op::operator()(const HairBSDF &bsdf) const {
    bool reflect = dot(vertex.geometry_normal, dir_in) *
                   dot(vertex.geometry_normal, dir_out) > 0;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometry_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    (void)reflect; // silence unuse warning, remove this when implementing hw

    // Get the parameterized positions 
    Real v = vertex.st.y; // v=[0,1] along the width of the tube, should be 0 for RTC_GEOMETRY_TYPE_ROUND_BEZIER_CURVE but not
    Real h = Real(-1) + Real(2) * v; // h=[-1, 1] point on the diameter of the intersection circle


    // Compute hair coordinate system terms related to _wo_
    Real sinThetaO = dir_in.x;
    Real cosThetaO = SafeSqrt(1 - Sqr(sinThetaO));
    Real phiO = std::atan2(dir_in.z, dir_in.y);

    // Compute hair coordinate system terms related to _wi_
    Real sinThetaI = dir_out.x;
    Real cosThetaI = SafeSqrt(1 - Sqr(sinThetaI));
    Real phiI = std::atan2(dir_out.z, dir_out.y);

    // Compute $\gammat$ for refracted ray
    Real etap = std::sqrt(bsdf.eta * bsdf.eta - Sqr(sinThetaO)) / cosThetaO;
    Real sinGammaT = h / etap;
    Real gammaT = SafeASin(sinGammaT);

    // Precompute the roughness v based on beta_m at each bounce
    Real longitudinal_v[pMax + 1];
    longitudinal_v[0] = Sqr(Real(0.726) * bsdf.beta_m + Real(0.812) * Sqr(bsdf.beta_m) + Real(3.7) * Pow<20>(bsdf.beta_m));
    longitudinal_v[1] = Real(0.25) * longitudinal_v[0];
    longitudinal_v[2] = 4 * longitudinal_v[0];
    for (int p = 3; p <= pMax; ++p) {
        longitudinal_v[p] = longitudinal_v[2];
    }

    Real s = sqrt(c_PI/Real(8)) * (Real(0.265) * bsdf.beta_n + Real(1.194) * Sqr(bsdf.beta_n) + Real(5.372) * Pow<22>(bsdf.beta_n));
    Real gammaO = SafeASin(h);

    // Compute $\alpha$ terms for hair scales
    Real sin2kAlpha[3];
    Real cos2kAlpha[3];
    sin2kAlpha[0] = std::sin(radians(bsdf.alpha));
    cos2kAlpha[0] = SafeSqrt(1 - Sqr(sin2kAlpha[0]));
    for (int i = 1; i < 3; ++i) {
        sin2kAlpha[i] = 2 * cos2kAlpha[i - 1] * sin2kAlpha[i - 1];
        cos2kAlpha[i] = Sqr(cos2kAlpha[i - 1]) - Sqr(sin2kAlpha[i - 1]);
    }

    std::array<Real, pMax + 1> apPdf = ComputeApPdf(cosThetaO, bsdf.eta, h, bsdf.sigma_a);

    // Compute PDF sum for hair scattering events
    Real phi = phiI - phiO;
    Real pdf = 0;
    for (int p = 0; p < pMax; ++p) {
        // Compute $\sin \thetao$ and $\cos \thetao$ terms accounting for scales
        Real sinThetaOp, cosThetaOp;
        if (p == 0) {
            sinThetaOp = sinThetaO * cos2kAlpha[1] - cosThetaO * sin2kAlpha[1];
            cosThetaOp = cosThetaO * cos2kAlpha[1] + sinThetaO * sin2kAlpha[1];
        }

        // Handle remainder of $p$ values for hair scale tilt
        else if (p == 1) {
            sinThetaOp = sinThetaO * cos2kAlpha[0] + cosThetaO * sin2kAlpha[0];
            cosThetaOp = cosThetaO * cos2kAlpha[0] - sinThetaO * sin2kAlpha[0];
        } else if (p == 2) {
            sinThetaOp = sinThetaO * cos2kAlpha[2] + cosThetaO * sin2kAlpha[2];
            cosThetaOp = cosThetaO * cos2kAlpha[2] - sinThetaO * sin2kAlpha[2];
        } else {
            sinThetaOp = sinThetaO;
            cosThetaOp = cosThetaO;
        }

        // Handle out-of-range $\cos \thetao$ from scale adjustment
        cosThetaOp = std::abs(cosThetaOp);
        pdf += Mp(cosThetaI, cosThetaOp, sinThetaI, sinThetaOp, longitudinal_v[p]) *
               apPdf[p] * Np(phi, p, s, gammaO, gammaT);
    }
    pdf += Mp(cosThetaI, cosThetaO, sinThetaI, sinThetaO, longitudinal_v[pMax]) *
           apPdf[pMax] * (1 / (2 * c_PI));

    // return Real(1);
    return pdf;

}



std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const HairBSDF &bsdf) const {
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometry_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!

    // Get the parameterized positions 
    Real v = vertex.st.y; // v=[0,1] along the width of the tube, should be 0 for RTC_GEOMETRY_TYPE_ROUND_BEZIER_CURVE but not
    Real h = Real(-1) + Real(2) * v; // h=[-1, 1] point on the diameter of the intersection circle


    // Compute hair coordinate system terms related to _wo_
    Real sinThetaO = dir_in.x;
    Real cosThetaO = SafeSqrt(1 - Sqr(sinThetaO));
    Real phiO = std::atan2(dir_in.z, dir_in.y);

    // Derive four random samples from _u2_
    Vector2 u[2] = {DemuxFloat(float(rnd_param_uv[0])), DemuxFloat(float(rnd_param_uv[1]))};

    // Determine which term $p$ to sample for hair scattering
    std::array<Real, pMax + 1> apPdf = ComputeApPdf(cosThetaO, bsdf.eta, h, bsdf.sigma_a);
    int p;
    for (p = 0; p < pMax; ++p) {
        if (u[0][0] < apPdf[p]) break;
        u[0][0] -= apPdf[p];
    }

    // Compute $\alpha$ terms for hair scales
    Real sin2kAlpha[3];
    Real cos2kAlpha[3];
    sin2kAlpha[0] = std::sin(radians(bsdf.alpha));
    cos2kAlpha[0] = SafeSqrt(1 - Sqr(sin2kAlpha[0]));
    for (int i = 1; i < 3; ++i) {
        sin2kAlpha[i] = 2 * cos2kAlpha[i - 1] * sin2kAlpha[i - 1];
        cos2kAlpha[i] = Sqr(cos2kAlpha[i - 1]) - Sqr(sin2kAlpha[i - 1]);
    }

    // Rotate $\sin \thetao$ and $\cos \thetao$ to account for hair scale tilt
    Real sinThetaOp, cosThetaOp;
    if (p == 0) {
        sinThetaOp = sinThetaO * cos2kAlpha[1] - cosThetaO * sin2kAlpha[1];
        cosThetaOp = cosThetaO * cos2kAlpha[1] + sinThetaO * sin2kAlpha[1];
    }
    else if (p == 1) {
        sinThetaOp = sinThetaO * cos2kAlpha[0] + cosThetaO * sin2kAlpha[0];
        cosThetaOp = cosThetaO * cos2kAlpha[0] - sinThetaO * sin2kAlpha[0];
    } else if (p == 2) {
        sinThetaOp = sinThetaO * cos2kAlpha[2] + cosThetaO * sin2kAlpha[2];
        cosThetaOp = cosThetaO * cos2kAlpha[2] - sinThetaO * sin2kAlpha[2];
    } else {
        sinThetaOp = sinThetaO;
        cosThetaOp = cosThetaO;
    }


    // Precompute the roughness v based on beta_m at each bounce
    Real longitudinal_v[pMax + 1];
    longitudinal_v[0] = Sqr(Real(0.726) * bsdf.beta_m + Real(0.812) * Sqr(bsdf.beta_m) + Real(3.7) * Pow<20>(bsdf.beta_m));
    longitudinal_v[1] = Real(0.25) * longitudinal_v[0];
    longitudinal_v[2] = 4 * longitudinal_v[0];
    for (int p = 3; p <= pMax; ++p) {
        longitudinal_v[p] = longitudinal_v[2];
    }

    Real s = sqrt(c_PI/Real(8)) * (Real(0.265) * bsdf.beta_n + Real(1.194) * Sqr(bsdf.beta_n) + Real(5.372) * Pow<22>(bsdf.beta_n));
    Real gammaO = SafeASin(h);


    // Sample $M_p$ to compute $\thetai$
    u[1][0] = std::max(u[1][0], Real(1e-5));
    Real cosTheta =
        1 + longitudinal_v[p] * std::log(u[1][0] + (1 - u[1][0]) * std::exp(-2 / longitudinal_v[p]));
    Real sinTheta = SafeSqrt(1 - Sqr(cosTheta));
    Real cosPhi = std::cos(2 * c_PI * u[1][1]);
    Real sinThetaI = -cosTheta * sinThetaOp + sinTheta * cosPhi * cosThetaOp;
    Real cosThetaI = SafeSqrt(1 - Sqr(sinThetaI));

    // Sample $N_p$ to compute $\Delta\phi$

    // Compute $\gammat$ for refracted ray
    Real etap = std::sqrt(bsdf.eta * bsdf.eta - Sqr(sinThetaO)) / cosThetaO;
    Real sinGammaT = h / etap;
    Real gammaT = SafeASin(sinGammaT);
    Real dphi;
    if (p < pMax)
        dphi =
            Phi(p, gammaO, gammaT) + SampleTrimmedLogistic(u[0][1], s, -c_PI, c_PI);
    else
        dphi = 2 * c_PI * u[0][1];
    

    // Compute _wi_ from sampled hair scattering angles
    Real phiI = phiO + dphi;
    Vector3 dir_out = Vector3(sinThetaI, cosThetaI * std::cos(phiI),
                   cosThetaI * std::sin(phiI));

    // Compute PDF for sampled hair scattering direction _wi_
    Real pdf = 0;
    for (int p = 0; p < pMax; ++p) {
        // Compute $\sin \thetao$ and $\cos \thetao$ terms accounting for scales
        Real sinThetaOp, cosThetaOp;
        if (p == 0) {
            sinThetaOp = sinThetaO * cos2kAlpha[1] - cosThetaO * sin2kAlpha[1];
            cosThetaOp = cosThetaO * cos2kAlpha[1] + sinThetaO * sin2kAlpha[1];
        }

        // Handle remainder of $p$ values for hair scale tilt
        else if (p == 1) {
            sinThetaOp = sinThetaO * cos2kAlpha[0] + cosThetaO * sin2kAlpha[0];
            cosThetaOp = cosThetaO * cos2kAlpha[0] - sinThetaO * sin2kAlpha[0];
        } else if (p == 2) {
            sinThetaOp = sinThetaO * cos2kAlpha[2] + cosThetaO * sin2kAlpha[2];
            cosThetaOp = cosThetaO * cos2kAlpha[2] - sinThetaO * sin2kAlpha[2];
        } else {
            sinThetaOp = sinThetaO;
            cosThetaOp = cosThetaO;
        }

        // Handle out-of-range $\cos \thetao$ from scale adjustment
        cosThetaOp = std::abs(cosThetaOp);
        pdf += Mp(cosThetaI, cosThetaOp, sinThetaI, sinThetaOp, longitudinal_v[p]) *
                apPdf[p] * Np(dphi, p, s, gammaO, gammaT);
    }
    pdf += Mp(cosThetaI, cosThetaO, sinThetaI, sinThetaO, longitudinal_v[pMax]) *
            apPdf[pMax] * (1 / (2 * c_PI));
    
    return BSDFSampleRecord{
        dir_out,
        Real(0) /* eta */, Real(1) /* roughness */};
}

TextureSpectrum get_texture_op::operator()(const HairBSDF &bsdf) const {
    return {};
}
