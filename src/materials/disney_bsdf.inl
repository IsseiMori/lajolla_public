#include "../microfacet.h"

Spectrum eval_op::operator()(const DisneyBSDF &bsdf) const {
    bool reflect = dot(vertex.geometry_normal, dir_in) *
                   dot(vertex.geometry_normal, dir_out) > 0;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometry_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    (void)reflect; // silence unuse warning, remove this when implementing hw

    Spectrum base_color = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real specular_transmission = eval(bsdf.specular_transmission, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real metallic = eval(bsdf.metallic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real subsurface = eval(bsdf.subsurface, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real specular = eval(bsdf.specular, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real specular_tint = eval(bsdf.specular_tint, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real sheen_tint = eval(bsdf.sheen_tint, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat = eval(bsdf.clearcoat, vertex.uv, vertex.uv_screen_size, texture_pool);
    roughness = std::clamp(roughness, Real(0.01), Real(1));

    Spectrum f_diffuse;
    Spectrum f_metal;
    Spectrum f_glass;
    Spectrum f_clearcoat;
    Spectrum f_sheen;

    /*-------------------- Diffuse --------------------*/
    {
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

        f_diffuse = (Real(1) - subsurface) * f_d + subsurface * f_ss;
    }


    /*-------------------- Metal --------------------*/
    {
        // Pre-compute reusable general terms
        Vector3 h = normalize(dir_in + dir_out); // half vector
        Vector3 n = frame.n;
        Real NdotIn = fabs(dot(n, dir_in));
        Real NdotOut = fabs(dot(n, dir_out));

        Real eta = dot(vertex.geometry_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;

        // Compute F
        // Spectrum Fm = base_color + (Real(1) - base_color)
        //             * pow(Real(1) - fabs(dot(h, dir_out)), 5);
        
        // Include an achromatic specular component
        Spectrum C_tint = base_color / luminance(base_color);
        if ( luminance(base_color) <= 0 ) C_tint = make_const_spectrum(Real(1));

        Spectrum Ks = (Real(1) - specular_tint) + specular_tint * C_tint;
        Spectrum C0 = specular * eta * (Real(1) - metallic) * Ks + metallic * base_color;
        Spectrum Fm = C0 + (Real(1) - C0) * pow(Real(1) - fabs(dot(h, dir_out)), 5);

        // Compute Dm
        Real aspect = sqrt(Real(1) - Real(0.9) * anisotropic);
        Real a_min = Real(0.0001);
        Real ax = fmax(a_min, roughness * roughness / aspect);
        Real ay = fmax(a_min, roughness * roughness * aspect);
        Real Dm = GTR2_aniso(ax, ay, frame, h);

        // Compute G
        Real Gin = smithG_GGX_aniso(NdotIn, dot(dir_in, frame.x), dot(dir_in, frame.y), ax, ay);
        Real Gout = smithG_GGX_aniso(NdotOut, dot(dir_out, frame.x), dot(dir_out, frame.y), ax, ay);

        f_metal = Fm * Dm * Gin * Gout / (Real(4) * fabs(dot(n, dir_in)));
    }


    /*-------------------- Glass --------------------*/
    {
        Real eta = dot(vertex.geometry_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;

        Vector3 half_vector;
        if (reflect) {
            half_vector = normalize(dir_in + dir_out);
        } else {
            // "Generalized half-vector" from Walter et al.
            // See "Microfacet Models for Refraction through Rough Surfaces"
            half_vector = normalize(dir_in + dir_out * eta);
        }

        // Flip half-vector if it's below surface
        if (dot(half_vector, frame.n) < 0) {
            half_vector = -half_vector;
        }

        // Clamp roughness to avoid numerical issues.
        roughness = std::clamp(roughness, Real(0.01), Real(1));

        Real h_dot_in = dot(half_vector, dir_in);
        Real F = fresnel_dielectric(h_dot_in, eta);
        Real D = GTR2(dot(frame.n, half_vector), roughness);
        Real G = smith_masking_gtr2(to_local(frame, dir_in), roughness) *
                smith_masking_gtr2(to_local(frame, dir_out), roughness);
        if (reflect) {
            f_glass = base_color * (F * D * G) / (4 * fabs(dot(frame.n, dir_in)));
        } else {
            Real h_dot_out = dot(half_vector, dir_out);
            
            f_glass = sqrt(base_color) * (Real(1) - F) * D * G * fabs(h_dot_out * h_dot_in) / 
                (fabs(dot(frame.n, dir_in)) * pow(h_dot_in + eta * h_dot_out, 2));
        }
    }

    /*-------------------- Clearcoat --------------------*/
    {
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

        f_clearcoat = make_const_spectrum(F * D * G / (Real(4) * fabs(n_dot_in)));
    }

    /*-------------------- Sheen --------------------*/
    {
        // Pre-compute useful values
        Vector3 half_vector = normalize(dir_in + dir_out);
        Real n_dot_out = dot(frame.n, dir_out);


        Spectrum C_tint = base_color / luminance(base_color);
        if ( luminance(base_color) <= 0 ) C_tint = make_const_spectrum(Real(1));

        Spectrum C_sheen = (Real(1) - sheen_tint) + sheen_tint * C_tint;
        f_sheen = C_sheen * pow(Real(1) - fabs(dot(half_vector, dir_out)), 5) * fabs(n_dot_out);
    }


    Spectrum f_bsdf = (Real(1) - specular_transmission) * (Real(1) - metallic) * f_diffuse
        + (Real(1) - metallic) * f_sheen
        + (Real(1) - specular_transmission * (Real(1) - metallic)) * f_metal
        + Real(0.25) * clearcoat * f_clearcoat
        + (Real(1) - metallic) * specular_transmission * f_glass;


    return f_bsdf;
}

Real pdf_sample_bsdf_op::operator()(const DisneyBSDF &bsdf) const {
    bool reflect = dot(vertex.geometry_normal, dir_in) *
                   dot(vertex.geometry_normal, dir_out) > 0;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometry_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    (void)reflect; // silence unuse warning, remove this when implementing hw

    Real specular_transmission = eval(bsdf.specular_transmission, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real metallic = eval(bsdf.metallic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat = eval(bsdf.clearcoat, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real eta = dot(vertex.geometry_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;
    roughness = std::clamp(roughness, Real(0.01), Real(1));

    
    Real diffuse_weight = (Real(1) - metallic) * (Real(1) - specular_transmission);
    Real metal_weight = (Real(1) - specular_transmission * (Real(1) - metallic));
    Real glass_weight = (Real(1) - metallic) * specular_transmission;
    Real clearcoat_weight = Real(0.25) * clearcoat;

    if ( dot(vertex.geometry_normal, dir_in) <= 0 ) {
        diffuse_weight = 0;
        metal_weight = 0;
        clearcoat_weight = 0;
    }

    Real weight_total = diffuse_weight + metal_weight + glass_weight + clearcoat_weight;
    diffuse_weight /= weight_total;
    metal_weight /= weight_total;
    glass_weight /= weight_total;
    clearcoat_weight /= weight_total;


    Real sum_pdf = Real(0);

    /*-------------------- Diffuse --------------------*/
    {
        sum_pdf += diffuse_weight * fmax(dot(frame.n, dir_out), Real(0)) / c_PI;
    }
    


    /*-------------------- Metal --------------------*/
    {
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

        sum_pdf +=  metal_weight * Dm * Gin / (Real(4) * fabs(dot(n, dir_in)));
    }



    /*-------------------- Glass --------------------*/
    {
        Vector3 half_vector;
        if (reflect) {
            half_vector = normalize(dir_in + dir_out);
        } else {
            // "Generalized half-vector" from Walter et al.
            // See "Microfacet Models for Refraction through Rough Surfaces"
            half_vector = normalize(dir_in + dir_out * eta);
        }

        // Flip half-vector if it's below surface
        if (dot(half_vector, frame.n) < 0) {
            half_vector = -half_vector;
        }

        Real roughness = eval(
            bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
        // Clamp roughness to avoid numerical issues.
        roughness = std::clamp(roughness, Real(0.01), Real(1));

        // We sample the visible normals, also we use F to determine
        // whether to sample reflection or refraction
        // so PDF ~ F * D * G_in for reflection, PDF ~ (1 - F) * D * G_in for refraction.
        Real h_dot_in = dot(half_vector, dir_in);
        Real F = fresnel_dielectric(h_dot_in, eta);
        Real D = GTR2(dot(half_vector, frame.n), roughness);
        Real G_in = smith_masking_gtr2(to_local(frame, dir_in), roughness);
        if (reflect) {
            sum_pdf +=  glass_weight *  (F * D * G_in) / (4 * fabs(dot(frame.n, dir_in)));
        } else {
            Real h_dot_out = dot(half_vector, dir_out);
            Real sqrt_denom = h_dot_in + eta * h_dot_out;
            Real dh_dout = eta * eta * h_dot_out / (sqrt_denom * sqrt_denom);
            sum_pdf +=  glass_weight *  (1 - F) * D * G_in * fabs(dh_dout * h_dot_in / dot(frame.n, dir_in));
        }
    }

    /*-------------------- Clearcoat --------------------*/
    {
        // Pre-compute useful values
        Vector3 half_vector = normalize(dir_in + dir_out);
        Real n_dot_h = dot(frame.n, half_vector);
        Real n_dot_out = dot(frame.n, dir_out);

        Real clearcoat_gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);
        
        Real D = compute_D(clearcoat_gloss, pow(dot(frame.n, half_vector),2));
    
        sum_pdf +=  clearcoat_weight * D * fabs(n_dot_h) / (Real(4) * fabs(n_dot_out));
    }

    return sum_pdf;
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneyBSDF &bsdf) const {
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometry_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!

    Real specular_transmission = eval(bsdf.specular_transmission, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real metallic = eval(bsdf.metallic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat = eval(bsdf.clearcoat, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat_gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real eta = dot(vertex.geometry_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;
    roughness = std::clamp(roughness, Real(0.01), Real(1));
    
    Real diffuse_weight = (Real(1) - metallic) * (Real(1) - specular_transmission);
    Real metal_weight = (Real(1) - specular_transmission * (Real(1) - metallic));
    Real glass_weight = (Real(1) - metallic) * specular_transmission;
    Real clearcoat_weight = Real(0.25) * clearcoat;

    if ( dot(vertex.geometry_normal, dir_in) <= 0 ) {
        diffuse_weight = 0;
        metal_weight = 0;
        clearcoat_weight = 0;
    }

    Real weight_total = diffuse_weight + metal_weight + glass_weight + clearcoat_weight;
    diffuse_weight /= weight_total;
    metal_weight /= weight_total;
    glass_weight /= weight_total;
    clearcoat_weight /= weight_total;

    Real rand = rnd_param_uv[0];

    if ( rand < diffuse_weight ) {

        /*-------------------- Diffuse --------------------*/

        return BSDFSampleRecord{
            to_world(frame, sample_cos_hemisphere(rnd_param_uv)),
            Real(0) /* eta */, Real(1) /* roughness */};

    } else if ( rand < diffuse_weight + metal_weight ) {

        /*-------------------- Metal --------------------*/

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

    } else if ( rand < diffuse_weight + metal_weight + glass_weight ) {

        /*-------------------- Glass --------------------*/

        // Sample a micro normal and transform it to world space -- this is our half-vector.
        Real alpha = roughness * roughness;
        Vector3 local_dir_in = to_local(frame, dir_in);
        Vector3 local_micro_normal =
            sample_visible_normals(local_dir_in, alpha, rnd_param_uv);

        Vector3 half_vector = to_world(frame, local_micro_normal);
        // Flip half-vector if it's below surface
        if (dot(half_vector, frame.n) < 0) {
            half_vector = -half_vector;
        }

        // Now we need to decide whether to reflect or refract.
        // We do this using the Fresnel term.
        Real h_dot_in = dot(half_vector, dir_in);
        Real F = fresnel_dielectric(h_dot_in, eta);

        if (rnd_param_w <= F) {
            // Reflection
            Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, half_vector) * half_vector);
            // set eta to 0 since we are not transmitting
            return BSDFSampleRecord{reflected, Real(0) /* eta */, roughness};
        } else {
            // Refraction
            // https://en.wikipedia.org/wiki/Snell%27s_law#Vector_form
            // (note that our eta is eta2 / eta1, and l = -dir_in)
            Real h_dot_out_sq = 1 - (1 - h_dot_in * h_dot_in) / (eta * eta);
            if (h_dot_out_sq <= 0) {
                // Total internal reflection
                // This shouldn't really happen, as F will be 1 in this case.
                return {};
            }
            // flip half_vector if needed
            if (h_dot_in < 0) {
                half_vector = -half_vector;
            }
            Real h_dot_out= sqrt(h_dot_out_sq);
            Vector3 refracted = -dir_in / eta + (fabs(h_dot_in) / eta - h_dot_out) * half_vector;
            return BSDFSampleRecord{refracted, eta, roughness};
        }

    } else if ( rand < diffuse_weight + metal_weight + glass_weight + clearcoat_weight ) {

        /*-------------------- Clearcoat --------------------*/

        Real alpha = (Real(1) - clearcoat_gloss) * Real(0.1) + 
                 clearcoat_gloss * Real(0.001);
        Real alpha2 = alpha * alpha;

        Real cos_h_elevation = sqrt((Real(1) - pow(alpha2, Real(1) - rnd_param_uv[0])) / (Real(1) - alpha2));
        Real h_elevation = acos(cos_h_elevation);
        Real h_azimuth = Real(2) * c_PI * rnd_param_uv[1];
        Real hlx = sin(h_elevation) * cos(h_azimuth);
        Real hly = sin(h_elevation) * sin(h_azimuth);
        Real hlz = cos(h_elevation);
        Vector3 h = Vector3(hlx, hly, hlz);
        Vector3 hemi_N = to_world(h, dir_in);

        return BSDFSampleRecord{
            hemi_N,
            Real(0) /* eta */, Real(1) /* roughness */};
    }

    return BSDFSampleRecord{
            Vector3(0,0,0),
            Real(0) /* eta */, Real(1) /* roughness */};
}

TextureSpectrum get_texture_op::operator()(const DisneyBSDF &bsdf) const {
    return bsdf.base_color;
}
