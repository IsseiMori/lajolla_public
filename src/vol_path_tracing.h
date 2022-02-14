#pragma once

// The simplest volumetric renderer: 
// single absorption only homogeneous volume
// only handle directly visible light sources
Spectrum vol_path_tracing_1(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!

    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                      (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    RayDifferential ray_diff = RayDifferential{Real(0), Real(0)};

    std::optional<PathVertex> vertex_ = intersect(scene, ray, ray_diff);

    if ( vertex_ ) {

        PathVertex vertex = *vertex_;

        if ( vertex.exterior_medium_id != -1 ) {
            const Medium exterior_media = scene.media[vertex.exterior_medium_id];
            Real t_hit = distance(vertex.position, ray.org);
            Spectrum sigma_a = get_sigma_a(exterior_media, vertex.position);

            Spectrum transmittance = exp(-sigma_a * t_hit);
            Spectrum Le = make_zero_spectrum();
            
            if (is_light(scene.shapes[vertex.shape_id])) {
                Le = emission(vertex, -ray.dir, scene);
            }

            return transmittance * Le;
        }
    }

    return make_zero_spectrum();
}

// The second simplest volumetric renderer: 
// single monochromatic homogeneous volume with single scattering,
// no need to handle surface lighting, only directly visible light source
Spectrum vol_path_tracing_2(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!

    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                      (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    RayDifferential ray_diff = RayDifferential{Real(0), Real(0)};

    std::optional<PathVertex> vertex_ = intersect(scene, ray, ray_diff);


    PathVertex vertex = *vertex_;

    int medium_id;
    Real t_hit;
    if ( vertex_ ) {
        medium_id = vertex.exterior_medium_id;
        t_hit = distance(vertex.position, ray.org);
    }
    else {
        medium_id = scene.camera.medium_id;
        t_hit = Real(INFINITY);
    }
    
    const Medium exterior_media = scene.media[medium_id];


    Spectrum sigma_s = get_sigma_s(exterior_media, vertex.position);
    Spectrum sigma_a = get_sigma_a(exterior_media, vertex.position);
    Spectrum sigma_t = sigma_s + sigma_a;
    
    Real u = next_pcg32_real<Real>(rng);
    Real t = -log(1 - u) / sigma_t.x;

    if ( t < t_hit ) {
        Spectrum trans_pdf = exp(-sigma_t * t) * sigma_t;
        Spectrum transmittance = exp(-sigma_t * t);

        Vector3 p = ray.org + t * ray.dir;
        
        // Equation 7
        // L_s1_estimate, L_s1_pdf = L_s1(p, sample_point_on_light(rng))
        // return (transmittance / trans_pdf) * sigma_s * (L_s1_estimate / L_s1_pdf)

        // Sample a point on the light source
        // by picking a light source, then pick a point on it
        Vector2 light_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
        Real light_w = next_pcg32_real<Real>(rng);
        Real shape_w = next_pcg32_real<Real>(rng);
        int light_id = sample_light(scene, light_w);
        const Light &light = scene.lights[light_id];
        PointAndNormal point_on_light =
            sample_point_on_light(light, p, light_uv, shape_w, scene);
        

        PhaseFunction phase_function = get_phase_function(exterior_media);
        Vector3 dir_view = -ray.dir;
        Vector3 dir_light = normalize(point_on_light.position - p);
        Spectrum rho = eval(phase_function, dir_view, dir_light);

        Spectrum Le = emission(light, -dir_light, Real(0), point_on_light, scene);

        Spectrum exp_term = exp(-sigma_t * (distance(p, point_on_light.position)));

        Real visibility = Real(1);
        Ray shadow_ray{p, dir_light, 
                    get_shadow_epsilon(scene),
                    (1 - get_shadow_epsilon(scene)) *
                        distance(point_on_light.position, p)};
        if (occluded(scene, shadow_ray)) {
            visibility = Real(0);
        }

        Real jacobian = fabs(dot(dir_light, point_on_light.normal)) / 
                            distance_squared(p, point_on_light.position) * 
                            visibility;
        
        Spectrum L_s1_estimate = rho * Le * exp_term * jacobian;

        Real L_s1_pdf = light_pmf(scene, light_id) *
            pdf_point_on_light(light, point_on_light, p, scene);

        return (transmittance / trans_pdf) * sigma_s * (L_s1_estimate / L_s1_pdf);

    }
    else {
        Spectrum trans_pdf = exp(-sigma_t * t_hit);
        Spectrum transmittance = exp(-sigma_t * t_hit);
        Spectrum Le = make_const_spectrum(0);

        if (is_light(scene.shapes[vertex.shape_id])) {
            Le = emission(vertex, -ray.dir, scene);
        }

        return transmittance / trans_pdf * Le;
    }
}

// The third volumetric renderer (not so simple anymore): 
// multiple monochromatic homogeneous volumes with multiple scattering
// no need to handle surface lighting, only directly visible light source
Spectrum vol_path_tracing_3(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
    return make_zero_spectrum();
}

// The fourth volumetric renderer: 
// multiple monochromatic homogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// still no surface lighting
Spectrum vol_path_tracing_4(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
    return make_zero_spectrum();
}

// The fifth volumetric renderer: 
// multiple monochromatic homogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// with surface lighting
Spectrum vol_path_tracing_5(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
    return make_zero_spectrum();
}

// The final volumetric renderer: 
// multiple chromatic heterogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// with surface lighting
Spectrum vol_path_tracing(const Scene &scene,
                          int x, int y, /* pixel coordinates */
                          pcg32_state &rng) {
    // Homework 2: implememt this!
    return make_zero_spectrum();
}
