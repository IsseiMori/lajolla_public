#include "../microfacet.h"

Real schlick_fresnel(Vector3 half_vector, Vector3 dir_out) {
    Real eta = Real(1.5);
    Real R_0 = pow(eta - Real(1), 2) / pow(eta + Real(1), 2);
    Real F_c = R_0 + (Real(1) - R_0) * pow(Real(1) - fabs(dot(half_vector, dir_out)), 5);
    return F_c;
}

Real compute_D(Real clearcoat_gloss, Real hlz2) {
    Real alpha = (Real(1) - clearcoat_gloss) * 
                  Real(0.1) + clearcoat_gloss * Real(0.001);

    Real alpha2 = alpha * alpha;

    return (alpha2 - Real(1)) / (c_PI * log(alpha2) * (Real(1) + (alpha - Real(1) * hlz2)));
}

Spectrum eval_op::operator()(const DisneyClearcoat &bsdf) const {
    if (dot(vertex.geometry_normal, dir_in) < 0 ||
            dot(vertex.geometry_normal, dir_out) < 0) {
        // No light below the surface
        return make_zero_spectrum();
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!

    // Pre-compute useful values
    Vector3 half_vector = normalize(dir_in + dir_out);
    Real n_dot_h = dot(frame.n, half_vector);
    Real n_dot_in = dot(frame.n, dir_in);
    if ( n_dot_h <= 0 ) {
        return make_zero_spectrum();
    }

    Real clearcoat_gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);
    
    Real F = schlick_fresnel(half_vector, dir_out);
    Real D = compute_D(clearcoat_gloss, pow(dot(frame.n, half_vector),2));
    Real G = smith_masking_gtr2(to_local(frame, dir_in), Real(0.05)) *
             smith_masking_gtr2(to_local(frame, dir_out), Real(0.05));

    return make_const_spectrum(F * D * G / (Real(4) * fabs(n_dot_in)));
}

Real pdf_sample_bsdf_op::operator()(const DisneyClearcoat &bsdf) const {
    if (dot(vertex.geometry_normal, dir_in) < 0 ||
            dot(vertex.geometry_normal, dir_out) < 0) {
        // No light below the surface
        return 0;
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!

    // Pre-compute useful values
    Vector3 half_vector = normalize(dir_in + dir_out);
    Real n_dot_h = dot(frame.n, half_vector);
    Real n_dot_out = dot(frame.n, dir_out);

    Real clearcoat_gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);
    
    Real D = compute_D(clearcoat_gloss, pow(dot(frame.n, half_vector),2));
   
    return D * fabs(n_dot_h) / (Real(4) * fabs(n_dot_out));
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneyClearcoat &bsdf) const {
    if (dot(vertex.geometry_normal, dir_in) < 0) {
        // No light below the surface
        return {};
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!

    return BSDFSampleRecord{
        to_world(frame, sample_cos_hemisphere(rnd_param_uv)),
        Real(0) /* eta */, Real(1) /* roughness */};
}

TextureSpectrum get_texture_op::operator()(const DisneyClearcoat &bsdf) const {
    return make_constant_spectrum_texture(make_zero_spectrum());
}
