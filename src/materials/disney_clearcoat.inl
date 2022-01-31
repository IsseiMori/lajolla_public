#include "../microfacet.h"

Real schlick_fresnel(Vector3 half_vector, Vector3 dir_out) {
    Real eta = Real(1.5);
    Real R_0 = pow(eta - Real(1), 2) / pow(eta + Real(1), 2);
    Real F_c = R_0 + (Real(1) - R_0) * pow(Real(1) - fabs(dot(half_vector, dir_out)), 5);
    return F_c;
}

Real compute_Dc(Real clearcoat_gloss, Real hlz2) {
    Real a = (Real(1) - clearcoat_gloss) * 
                  Real(0.1) + clearcoat_gloss * Real(0.001);
    Real a2 = a * a;

    return (a2 - Real(1)) / ( c_PI * log(a2) * (Real(1) + (a2 - Real(1)) * hlz2) );
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
    Real D = compute_Dc(clearcoat_gloss, pow(dot(frame.n, half_vector),2));
    Real G = smith_masking_gtr2(to_local(frame, dir_in), Real(0.5)) *
             smith_masking_gtr2(to_local(frame, dir_out), Real(0.5));

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

    Real clearcoat_gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);
    
    Real D = compute_Dc(clearcoat_gloss, pow(dot(frame.n, half_vector),2));
   
    return D * fabs(n_dot_h) / (Real(4) * fabs(dot(half_vector, dir_out)));
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

    Real clearcoat_gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);
    
    Real a = (Real(1) - clearcoat_gloss) * Real(0.1) + clearcoat_gloss * Real(0.001);
    Real a2 = a * a;

    Real cos_h_elevation = sqrt((Real(1) - pow(a2, Real(1) - rnd_param_uv[0])) / (Real(1) - a2));
    Real h_elevation = acos(cos_h_elevation);
    Real h_azimuth = Real(2) * c_PI * rnd_param_uv[1];
    Real hlx = sin(h_elevation) * cos(h_azimuth);
    Real hly = sin(h_elevation) * sin(h_azimuth);
    Real hlz = cos(h_elevation);
    Vector3 h = normalize(Vector3(hlx, hly, hlz));

    Vector3 half_vector = to_world(frame, h);

    Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, half_vector) * half_vector);

    return BSDFSampleRecord{
        reflected,
        Real(0) /* eta */, Real(1) /* roughness */};
}

TextureSpectrum get_texture_op::operator()(const DisneyClearcoat &bsdf) const {
    return make_constant_spectrum_texture(make_zero_spectrum());
}
