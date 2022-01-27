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

Vector3 sample_visible_normals_aniso(const Vector3 &local_dir_in, Real ax, Real ay, const Vector2 &rnd_param) {
    // The incoming direction is in the "ellipsodial configuration" in Heitz's paper
    if (local_dir_in.z < 0) {
        // Ensure the input is on top of the surface.
        return -sample_visible_normals_aniso(-local_dir_in, ax, ay, rnd_param);
    }

    // Transform the incoming direction to the "hemisphere configuration".
    Vector3 hemi_dir_in = normalize(
        Vector3{ax * local_dir_in.x, ay * local_dir_in.y, local_dir_in.z});

    // Parameterization of the projected area of a hemisphere.
    // First, sample a disk.
    Real r = sqrt(rnd_param.x);
    Real phi = 2 * c_PI * rnd_param.y;
    Real t1 = r * cos(phi);
    Real t2 = r * sin(phi);
    // Vertically scale the position of a sample to account for the projection.
    Real s = (1 + hemi_dir_in.z) / 2;
    t2 = (1 - s) * sqrt(1 - t1 * t1) + s * t2;
    // Point in the disk space
    Vector3 disk_N{t1, t2, sqrt(max(Real(0), 1 - t1*t1 - t2*t2))};

    // Reprojection onto hemisphere -- we get our sampled normal in hemisphere space.
    Frame hemi_frame(hemi_dir_in);
    Vector3 hemi_N = to_world(hemi_frame, disk_N);

    // Transforming the normal back to the ellipsoid configuration
    return normalize(Vector3{ax * hemi_N.x, ay * hemi_N.y, max(Real(0), hemi_N.z)});
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
    Vector3 half_vector = normalize(dir_in + dir_out);
    Real h_dot_out = dot(half_vector, dir_out);

    // Compute F
    Spectrum Fm = base_color + (Real(1) - base_color)
                * pow(Real(1) - fabs(h_dot_out), 5);

    // Compute Dm
    Real aspect = sqrt(Real(1) - Real(0.9) * anisotropic);
    Real a_min = Real(0.0001);
    Real ax = fmax(a_min, roughness * roughness / aspect);
    Real ay = fmax(a_min, roughness * roughness * aspect);
    Real Dm = GTR2_aniso(ax, ay, frame, half_vector);

    // Compute G
    Real Gin = smithG_GGX_aniso(dot(dir_in, frame.n), dot(dir_in, frame.x), dot(dir_in, frame.y), ax, ay);
    Real Gout = smithG_GGX_aniso(dot(dir_out, frame.n), dot(dir_out, frame.x), dot(dir_out, frame.y), ax, ay);

    return Fm * Dm * Gin * Gout / (Real(4) * fabs(dot(dir_in, frame.n)));
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
    Vector3 half_vector = normalize(dir_in + dir_out);

    // Compute Dm
    Real aspect = sqrt(Real(1) - Real(0.9) * anisotropic);
    Real a_min = Real(0.0001);
    Real ax = fmax(a_min, roughness * roughness / aspect);
    Real ay = fmax(a_min, roughness * roughness * aspect);
    Real Dm = GTR2_aniso(ax, ay, frame, half_vector);

    // Compute G
    Real Gin = smithG_GGX_aniso(dot(dir_in, frame.n), dot(dir_in, frame.x), dot(dir_in, frame.y), ax, ay);

    return Dm * Gin / (Real(4) * fabs(dot(dir_in, frame.n)));
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
    
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    roughness = std::clamp(roughness, Real(0.01), Real(1));

    // Clamp roughness to avoid numerical issues.
    Real aspect = sqrt(Real(1) - Real(0.9) * anisotropic);
    Real a_min = Real(0.0001);
    Real ax = fmax(a_min, roughness * roughness / aspect);
    Real ay = fmax(a_min, roughness * roughness * aspect);

    Vector3 local_micro_normal =
        sample_visible_normals_aniso(local_dir_in, ax, ay, rnd_param_uv);
    
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
