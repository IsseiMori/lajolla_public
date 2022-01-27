Spectrum eval_op::operator()(const DisneyDiffuse &bsdf) const {
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
    Real subsurface = eval(bsdf.subsurface, vertex.uv, vertex.uv_screen_size, texture_pool);

    // Pre-compute reusable general terms
    Vector3 half_vector = normalize(dir_in + dir_out);
    Real h_dot_in = dot(half_vector, dir_in);
    Real h_dot_out = dot(half_vector, dir_out);

    // Compute the base diffuse term
    Real FD90 = Real(0.5) + Real(2) * roughness * h_dot_out * h_dot_out;
    Real FD_in = Real(1) + (FD90 - Real(1)) * (Real(1) - pow(h_dot_in,5));
    Real FD_out = Real(1) + (FD90 - Real(1)) * (Real(1) - pow(h_dot_out,5));
    Spectrum f_d = base_color * FD_in * FD_out * h_dot_out / c_PI;

    // Compute the subsurface term
    Real FSS90 = roughness * h_dot_out * h_dot_out;
    Real FSS_in = Real(1) * (FSS90 - Real(1)) * (Real(1) - pow(h_dot_in, 5));
    Real FSS_out = Real(1) * (FSS90 - Real(1)) * (Real(1) - pow(h_dot_out, 5));
    Spectrum f_ss = Real(1.25) * base_color 
                    * (FSS_in * FSS_out * (Real(1) / (h_dot_in + h_dot_out) - Real(0.5)) + Real(0.5))
                    * h_dot_out / c_PI;

    return (Real(1) - subsurface) * f_d + subsurface * f_ss;
}

Real pdf_sample_bsdf_op::operator()(const DisneyDiffuse &bsdf) const {
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
    return fmax(dot(frame.n, dir_out), Real(0)) / c_PI;
}

std::optional<BSDFSampleRecord> sample_bsdf_op::operator()(const DisneyDiffuse &bsdf) const {
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

TextureSpectrum get_texture_op::operator()(const DisneyDiffuse &bsdf) const {
    return bsdf.base_color;
}
