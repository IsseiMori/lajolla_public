#include "../microfacet.h"

Real smithG_GGX_aniso(Real NdotW, Real WdotX, Real WdotY, Real ax, Real ay) {

    Real lambda = Real(0.5) * (sqrt(Real(1) + (pow(WdotX * ax, 2) + pow(WdotY * ay, 2)) / pow(NdotW, 2)) - Real(1));

    return Real(1) / (Real(1) + lambda);
}

Real GTR2_aniso(Real ax, Real ay, Frame frame, Vector3 h) {
    Real ax2 = ax * ax;
    Real ay2 = ay * ay;
    Real hlx2 = pow(dot(frame.x, h),2);
    Real hly2 = pow(dot(frame.y, h),2);
    Real hlz2 = pow(dot(frame.n, h),2);
    Real Dm = Real(1) / (c_PI * ax * ay * pow(hlx2 / ax2 + hly2 / ay2 + hlz2, 2));

    return Dm;
}

Spectrum eval_op::operator()(const DisneyMetal &bsdf) const {
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
    Spectrum base_color = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    roughness = std::clamp(roughness, Real(0.01), Real(1));

    // Pre-compute reusable general terms
    Vector3 h = normalize(dir_in + dir_out); // half vector
    Vector3 n = frame.n;
    Real NdotIn = fabs(dot(n, dir_in));
    Real NdotOut = fabs(dot(n, dir_out));

    // Compute F
    Spectrum Fm = base_color + (Real(1) - base_color)
                * pow(Real(1) - fabs(dot(h, dir_out)), 5);

    // Compute Dm
    Real aspect = sqrt(Real(1) - Real(0.9) * anisotropic);
    Real a_min = Real(0.0001);
    Real ax = fmax(a_min, roughness * roughness / aspect);
    Real ay = fmax(a_min, roughness * roughness * aspect);
    Real Dm = GTR2_aniso(ax, ay, frame, h);

    // Compute G
    Real Gin = smithG_GGX_aniso(NdotIn, dot(dir_in, frame.x), dot(dir_in, frame.y), ax, ay);
    Real Gout = smithG_GGX_aniso(NdotOut, dot(dir_out, frame.x), dot(dir_out, frame.y), ax, ay);

    return Fm * Dm * Gin * Gout / (Real(4) * fabs(dot(n, dir_in)));
}

Real pdf_sample_bsdf_op::operator()(const DisneyMetal &bsdf) const {
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

    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    roughness = std::clamp(roughness, Real(0.01), Real(1));

    // Pre-compute reusable general terms
    Vector3 h = normalize(dir_in + dir_out); // half vector
    Vector3 n = frame.n;
    Real NdotIn = fabs(dot(n, dir_in));

    // Compute Dm
    Real aspect = sqrt(Real(1) - Real(0.9) * anisotropic);
    Real a_min = Real(0.0001);
    Real ax = fmax(a_min, roughness * roughness / aspect);
    Real ay = fmax(a_min, roughness * roughness * aspect);
    Real Dm = GTR2_aniso(ax, ay, frame, h);

    // Compute G
    Real Gin = smithG_GGX_aniso(NdotIn, dot(dir_in, frame.x), dot(dir_in, frame.y), ax, ay);

    return Dm * Gin / (Real(4) * fabs(dot(n, dir_in)));
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneyMetal &bsdf) const {
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

    // Convert the incoming direction to local coordinates
    Vector3 local_dir_in = to_local(frame, dir_in);
    Real roughness = eval(
        bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    
    // Clamp roughness to avoid numerical issues.
    roughness = std::clamp(roughness, Real(0.01), Real(1));
    Real alpha = roughness * roughness;
    Vector3 local_micro_normal =
        sample_visible_normals(local_dir_in, alpha, rnd_param_uv);
    
    // Transform the micro normal to world space
    Vector3 half_vector = to_world(frame, local_micro_normal);
    // Reflect over the world space normal
    Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, half_vector) * half_vector);
    return BSDFSampleRecord{
        reflected,
        Real(0) /* eta */, roughness /* roughness */
    };
}

TextureSpectrum get_texture_op::operator()(const DisneyMetal &bsdf) const {
    return bsdf.base_color;
}
