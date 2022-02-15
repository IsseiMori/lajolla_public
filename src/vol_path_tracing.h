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
        // Scatter once before hitting a surface
        // Compute the radiative transfer equation with single recursion

        Spectrum trans_pdf = exp(-sigma_t * t) * sigma_t;
        Spectrum transmittance = exp(-sigma_t * t);

        Vector3 p = ray.org + t * ray.dir;

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

        // hit a surface, account for surface emission

        Spectrum trans_pdf = exp(-sigma_t * t_hit);
        Spectrum transmittance = exp(-sigma_t * t_hit);
        Spectrum Le = make_zero_spectrum();

        if (is_light(scene.shapes[vertex.shape_id])) {
            Le = emission(vertex, -ray.dir, scene);
        }

        return transmittance / trans_pdf * Le;
    }
}

int update_medium(PathVertex vertex, Ray ray) {

    int medium = -1;

    if ( vertex.interior_medium_id != vertex.exterior_medium_id ) {

        if ( dot(ray.dir, vertex.geometry_normal) > 0 ) {
            medium = vertex.exterior_medium_id;
        }
        else {
            medium = vertex.interior_medium_id;
        }
    }

    return medium;

}

// The third volumetric renderer (not so simple anymore): 
// multiple monochromatic homogeneous volumes with multiple scattering
// no need to handle surface lighting, only directly visible light source
Spectrum vol_path_tracing_3(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!

    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                      (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    RayDifferential ray_diff = RayDifferential{Real(0), Real(0)};

    int current_medium = scene.camera.medium_id;
    Spectrum current_path_throughput = make_const_spectrum(1);
    Spectrum radiance = make_zero_spectrum();
    int bounces = 0;


    while ( true ) {
        
        bool scatter = false;

        std::optional<PathVertex> vertex_ = intersect(scene, ray, ray_diff);
        PathVertex vertex = *vertex_;

        Spectrum transmittance = make_const_spectrum(1);
        Spectrum trans_pdf = make_const_spectrum(1);


        // if ( current_medium ) 
        {
            const Medium media = scene.media[current_medium];
            Spectrum sigma_s = get_sigma_s(media, vertex.position);
            Spectrum sigma_a = get_sigma_a(media, vertex.position);
            Spectrum sigma_t = sigma_s + sigma_a;
            
            Real u = next_pcg32_real<Real>(rng);
            Real t = -log(1 - u) / sigma_t.x;

            trans_pdf = exp(-sigma_t * t) * sigma_t;
            transmittance = exp(-sigma_t * t);

            Real t_hit = distance(vertex.position, ray.org);
            if ( !vertex_ ) t_hit = Real(INFINITY);

            if ( t < t_hit ) {
                scatter = true;

                Vector3 p = ray.org + t * ray.dir;

                // Sample a point on the light source
                // by picking a light source, then pick a point on it
                Vector2 light_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
                Real light_w = next_pcg32_real<Real>(rng);
                Real shape_w = next_pcg32_real<Real>(rng);
                int light_id = sample_light(scene, light_w);
                const Light &light = scene.lights[light_id];
                PointAndNormal point_on_light =
                    sample_point_on_light(light, p, light_uv, shape_w, scene);
                

                PhaseFunction phase_function = get_phase_function(media);
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

                current_path_throughput *= sigma_s * (L_s1_estimate / L_s1_pdf);
            }

            ray.org = ray.org + t * ray.dir + get_intersection_epsilon(scene);
        }


        current_path_throughput *= transmittance / trans_pdf;


        if ( !scatter ) {
            // reach a surface, include emission
            
            Spectrum Le = make_zero_spectrum();

            if (vertex_ && is_light(scene.shapes[vertex.shape_id])) {
                Le = emission(vertex, -ray.dir, scene);
            }

            radiance += current_path_throughput * Le;
        }

        if ( bounces == scene.options.max_depth - 1 && scene.options.max_depth != -1 ) {
            break;
        }

        if ( !scatter && vertex_ ) {
            if ( vertex.material_id == -1 ) {
                current_medium = update_medium(vertex, ray);
                bounces++;
                continue;
            }
        }

        if ( scatter ) {
            // return current_path_throughput * Spectrum(0,0,1);
            // sample next direct & update path throughput
            const Medium media = scene.media[current_medium];
            PhaseFunction phase_function = get_phase_function(media);
            const Vector2 phase_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
            std::optional<Vector3> next_dir_ = sample_phase_function(phase_function, -ray.dir, phase_uv);
            Vector3 &next_dir = *next_dir_;

            Spectrum sigma_s = get_sigma_s(media, vertex.position);

            Real phase_pdf = pdf_sample_phase(phase_function, -ray.dir, next_dir);

            current_path_throughput *= eval(phase_function, -ray.dir, next_dir) /
                                       phase_pdf *
                                       sigma_s;

            ray.dir = next_dir;
        }
        else {
            // Hit a surface -- don’t need to deal with this yet
            break;
        }

        Real rr_prob = Real(1);
        if ( bounces >= scene.options.rr_depth ) {
            // TODO: choosing R channel for now
            rr_prob = min(current_path_throughput.x, 0.95);
            Real u = next_pcg32_real<Real>(rng);
            if ( u > rr_prob ) {
                break;
            }
            else {
                current_medium /= rr_prob;
            }
        }

        bounces++;
    }

    return radiance;

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
