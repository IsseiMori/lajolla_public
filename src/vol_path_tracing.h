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
        // intersect a surface

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

int update_medium(PathVertex vertex, Ray ray, int medium) {

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



/*
Compute transmittance to light. Skip through index-matching shapes.
*/
Spectrum next_event_estimation(const Scene &scene, 
                           pcg32_state &rng, 
                           Vector3 p, 
                           int current_medium,
                           int bounces,
                           Vector3 dir_view, // previous vertex to current point
                           PathVertex vertex) {

    // Sample a point on light
    // p' (p_prime) is the point on the light source
    Vector2 light_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
    Real light_w = next_pcg32_real<Real>(rng);
    Real shape_w = next_pcg32_real<Real>(rng);
    int light_id = sample_light(scene, light_w);
    const Light &light = scene.lights[light_id];
    PointAndNormal point_on_light =
        sample_point_on_light(light, p, light_uv, shape_w, scene);
    Vector3 dir_light = normalize(point_on_light.position - p);
    Vector3 p_prime = point_on_light.position;
    Vector3 p_origin = p; // remember where p started



    int shadow_medium = current_medium;
    int shadow_bounces = 0;
    Spectrum transmittance_light = make_const_spectrum(1);
    Spectrum p_trans_dir = make_const_spectrum(1); // for multiple importance sampling

    PathVertex shadow_vertex;
    Real next_t = Real(0);


    while ( true ) {
        
        // Ray intersect towards p'
        Ray shadow_ray = Ray{p, dir_light, 
                             get_shadow_epsilon(scene), 
                             (1 - get_shadow_epsilon(scene)) *
                             distance(p, p_prime)};
        RayDifferential ray_diff = RayDifferential{Real(0), Real(0)};
        std::optional<PathVertex> shadow_vertex_ = intersect(scene, shadow_ray, ray_diff);
        shadow_vertex = *shadow_vertex_;
        

        // Set the distance to the nearest intersection
        next_t = distance(p, p_prime);
        if ( shadow_vertex_ ) {
            next_t = distance(p, shadow_vertex.position);
        }

        // Account for the transmittance to next_t
        if ( shadow_medium != -1 ) {
            const Medium media = scene.media[shadow_medium];
            Spectrum sigma_s = get_sigma_s(media, vertex.position);
            Spectrum sigma_a = get_sigma_a(media, vertex.position);
            Spectrum sigma_t = sigma_s + sigma_a;
            
            transmittance_light *= exp(-sigma_t * next_t);
            p_trans_dir *= exp(-sigma_t * next_t);
        }

        if ( !shadow_vertex_ ) {

            // Nothing is blocking, we’re done
            break;

        }
        else {

            // Something is blocking: is it an opaque surface?
            if ( shadow_vertex.material_id >= 0 ) {
                // we’re blocked
                return make_zero_spectrum();
            }

            // otherwise, it’s an index-matching surface and
            // we want to pass through -- this introduces
            // one extra connection vertex
            shadow_bounces++;

            if ( scene.options.max_depth != -1 && bounces + shadow_bounces >= scene.options.max_depth ) {
                // Reach the max no. of vertices
                return make_zero_spectrum();
            }

            shadow_medium = update_medium(shadow_vertex, shadow_ray, shadow_medium);
            p = p + next_t * dir_light;
        }
    }


    // TODO: checking only R channel 
    if ( max(transmittance_light) > 0 ) {
        
        // Compute T_light * G * f * L & pdf_nee

        const Medium media = scene.media[current_medium];

        PhaseFunction phase_function = get_phase_function(media);
        const Vector2 phase_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
        Spectrum rho = eval(phase_function, dir_view, dir_light);

        Spectrum Le = emission(light, -dir_light, Real(0), point_on_light, scene);
        
        Real jacobian = fabs(dot(dir_light, point_on_light.normal)) / 
                            distance_squared(p_origin, p_prime);

        Real pdf_nee = light_pmf(scene, light_id) *
                        pdf_point_on_light(light, point_on_light, p_origin, scene);

        Spectrum contrib = transmittance_light * rho * Le * jacobian / pdf_nee;


        Spectrum pdf_phase = pdf_sample_phase(phase_function, dir_view, dir_light) * jacobian * p_trans_dir;

        Spectrum w = (pdf_nee * pdf_nee) / (pdf_nee * pdf_nee + pdf_phase * pdf_phase);

        return contrib * w;

    }

    return make_zero_spectrum();

}


/*
Compute transmittance to light. Skip through index-matching shapes.
*/
Spectrum next_event_estimation_final(const Scene &scene, 
                           pcg32_state &rng, 
                           Vector3 p, 
                           int current_medium,
                           int bounces,
                           Vector3 dir_view, // previous vertex to current point
                           bool is_surface,
                           PathVertex vertex
                           ) {

    // Sample a point on light
    // p' (p_prime) is the point on the light source
    Vector2 light_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
    Real light_w = next_pcg32_real<Real>(rng);
    Real shape_w = next_pcg32_real<Real>(rng);
    int light_id = sample_light(scene, light_w);
    const Light &light = scene.lights[light_id];
    PointAndNormal point_on_light =
        sample_point_on_light(light, p, light_uv, shape_w, scene);
    Vector3 dir_light = normalize(point_on_light.position - p);
    Vector3 p_prime = point_on_light.position;
    Vector3 p_origin = p; // remember where p started



    int shadow_medium = current_medium;
    int shadow_bounces = 0;
    Spectrum transmittance_light = make_const_spectrum(1);

    Spectrum p_trans_nee = make_const_spectrum(1);
    Spectrum p_trans_dir = make_const_spectrum(1); // for multiple importance sampling

    PathVertex shadow_vertex;
    Real next_t = Real(0);


    while ( true ) {
        
        // Ray intersect towards p'
        Ray shadow_ray = Ray{p, dir_light, 
                             get_shadow_epsilon(scene), 
                             (1 - get_shadow_epsilon(scene)) *
                             distance(p, p_prime)};
        RayDifferential ray_diff = RayDifferential{Real(0), Real(0)};
        std::optional<PathVertex> shadow_vertex_ = intersect(scene, shadow_ray, ray_diff);
        shadow_vertex = *shadow_vertex_;
        

        // Set the distance to the nearest intersection
        next_t = distance(p, p_prime);
        if ( shadow_vertex_ ) {
            next_t = distance(p, shadow_vertex.position);
        }

        // Account for the transmittance to next_t
        if ( shadow_medium != -1 ) {

            const Medium &media = scene.media[shadow_medium];
            Spectrum majorant = get_majorant(media, shadow_ray);

            Real u = next_pcg32_real<Real>(rng);
            int channel = std::clamp(int(u * 3), 0, 2);
            Real accum_t = Real(0);
            int iteration = 0;

            while ( true ) {

                if ( majorant[channel] <= 0 ) {
                    break;
                }

                if ( iteration >= scene.options.max_null_collisions ) {
                    break;
                }

                Real t = -log(1 - next_pcg32_real<Real>(rng)) / majorant[channel];
                Real dt = next_t - accum_t;

                // Update accumulated distance
                accum_t = std::min(accum_t + t, next_t);

                if ( t < dt ) {

                    Vector3 pos = p + dir_light * accum_t;

                    Spectrum sigma_s = get_sigma_s(media, pos);
                    Spectrum sigma_a = get_sigma_a(media, pos);
                    Spectrum sigma_t = sigma_s + sigma_a;

                    Spectrum sigma_n = majorant * ( 1 - sigma_t / majorant );

                    Spectrum real_prob = sigma_t / majorant;
                    transmittance_light *= exp(-majorant * t) * sigma_n / max(majorant);
                    p_trans_nee *= exp(-majorant * t) * majorant / max(majorant);
                    p_trans_dir *= exp(-majorant * t) * majorant * (1- real_prob) / max(majorant);

                    if ( max(transmittance_light) <= 0 ) {
                        break;
                    }

                }
                else {
                    transmittance_light *= exp(-majorant * dt);
                    p_trans_nee *= exp(-majorant * dt);
                    p_trans_dir *= exp(-majorant * dt);
                    break;
                }

                iteration++;
            }

        }

        if ( !shadow_vertex_ ) {

            // Nothing is blocking, we’re done
            break;

        }
        else {

            // Something is blocking: is it an opaque surface?
            if ( shadow_vertex.material_id >= 0 ) {
                // we’re blocked
                return make_zero_spectrum();
            }

            // otherwise, it’s an index-matching surface and
            // we want to pass through -- this introduces
            // one extra connection vertex
            shadow_bounces++;

            if ( scene.options.max_depth != -1 && bounces + shadow_bounces >= scene.options.max_depth ) {
                // Reach the max no. of vertices
                return make_zero_spectrum();
            }

            shadow_medium = update_medium(shadow_vertex, shadow_ray, shadow_medium);
            p = p + next_t * dir_light;
        }
    }


    // TODO: checking only R channel 
    if ( max(transmittance_light) > 0 ) {
        

        Spectrum pdf_dir = make_zero_spectrum();
        Spectrum evaluated_f = make_zero_spectrum();

        Spectrum Le = emission(light, -dir_light, Real(0), point_on_light, scene);
        Real jacobian = max(-dot(dir_light, point_on_light.normal), Real(0)) / 
                        distance_squared(p_origin, p_prime);

        
        Spectrum pdf_nee = light_pmf(scene, light_id) *
                           pdf_point_on_light(light, point_on_light, p_origin, scene) * 
                           p_trans_nee;

        if ( is_surface ) {

            const Material &mat = scene.materials[vertex.material_id];

            evaluated_f = eval(mat, dir_view, dir_light, vertex, scene.texture_pool);

            Real pdf_bsdf = pdf_sample_bsdf(mat, dir_view, dir_light, vertex, scene.texture_pool);

            if (pdf_bsdf <= 0) {
                // Numerical issue -- we generated some invalid rays.
                return make_zero_spectrum();
            }

            pdf_dir = pdf_bsdf * jacobian * p_trans_dir;

        }
        else {

            const Medium &media = scene.media[current_medium];

            PhaseFunction phase_function = get_phase_function(media);
            const Vector2 phase_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
            evaluated_f = eval(phase_function, dir_view, dir_light);

            Real pdf_phase = pdf_sample_phase(phase_function, dir_view, dir_light);
            pdf_dir = pdf_phase * jacobian * p_trans_dir;

        }

        Spectrum contrib = transmittance_light * evaluated_f * Le * jacobian / average(pdf_nee);
        Spectrum w = (pdf_nee * pdf_nee) / (pdf_nee * pdf_nee + pdf_dir * pdf_dir);

        return contrib * w;
    }

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


    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                      (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    ray.tnear = get_intersection_epsilon(scene);
    RayDifferential ray_diff = RayDifferential{Real(0), Real(0)};

    int current_medium = scene.camera.medium_id;
    Spectrum current_path_throughput = make_const_spectrum(1);
    Spectrum radiance = make_zero_spectrum();
    int bounces = 0;

    Real dir_pdf = 0; // the pdf of the latest phase function sampling
    Vector3 nee_p_cache; // the last position p that can issue a next event estimation

    // flag to record if nee is issued to avoid including nee pdf contribution
    // when no nee has been issued
    bool is_nee_issued = false;
    
    // The product PDF of transmittance sampling going through 
    // several index-matching surfaces from the last phase function sampling
    Spectrum multi_trans_pdf = make_const_spectrum(1); 

    Real eta_scale = Real(1);


    // if ( !(x > 400 && x < 404) ) {
    //     return Spectrum(0,0,1);
    // }

    while ( true ) {
        
        bool scatter = false;

        std::optional<PathVertex> vertex_ = intersect(scene, ray, ray_diff);
        PathVertex vertex = *vertex_;

        Real t_hit = distance(vertex.position, ray.org);
        if ( !vertex_ ) t_hit = Real(INFINITY);

        Spectrum transmittance = make_const_spectrum(1);
        Spectrum trans_dir_pdf = make_const_spectrum(1); // PDF for free-flight sampling
        Spectrum trans_nee_pdf = make_const_spectrum(1); // PDF for next event estimation


        if ( current_medium != -1 ) {
            
            const Medium &media = scene.media[current_medium];
            Spectrum majorant = get_majorant(media, ray);

            
            Real u = next_pcg32_real<Real>(rng);
            int channel = std::clamp(int(u * 3), 0, 2);
            Real accum_t = Real(0);
            int iteration = 0;

            while ( true ) {

                if ( majorant[channel] <= 0 ) {
                    break;
                }

                if ( iteration >= scene.options.max_null_collisions ) {
                    break;
                }

                Real t = -log(1 - next_pcg32_real<Real>(rng)) / majorant[channel];
                Real dt = t_hit - accum_t;

                // Update accumulated distance
                accum_t = std::min(accum_t + t, t_hit);

                if ( t < dt ) {

                    Vector3 p = ray.org + ray.dir * accum_t;

                    Spectrum sigma_s = get_sigma_s(media, p);
                    Spectrum sigma_a = get_sigma_a(media, p);
                    Spectrum sigma_t = sigma_s + sigma_a;


                    Spectrum sigma_n = majorant * ( 1 - sigma_t / majorant );

                    Spectrum real_prob = sigma_t / majorant;
                    if ( next_pcg32_real<Real>(rng) < real_prob[channel] ) {
                        // Hit a "real" particle

                        scatter = true;
                        transmittance *= exp(-majorant * t) / max(majorant);
                        trans_dir_pdf *= exp(-majorant * t) * majorant * real_prob / max(majorant);


                        ray.org = p;

                        break;
                    }
                    else {
                        transmittance *= exp(-majorant * t) * sigma_n / max(majorant);
                        trans_dir_pdf *= exp(-majorant * t) * majorant * (1- real_prob) / max(majorant);
                        trans_nee_pdf *= exp(-majorant * t) * majorant / max(majorant);
                    }

                }
                else {
                    transmittance *= exp(-majorant * dt);
                    trans_dir_pdf *= exp(-majorant * dt);
                    trans_nee_pdf *= exp(-majorant * dt);
                    ray.org = vertex.position;
                    break;
                }

                iteration++;
            }


            // Accumulate the pdf of the trasmittance sampling
            // since the previous phase function sampling
            multi_trans_pdf *= trans_dir_pdf;


        }
        else {
            // For vaccume volume, we simply move the ray origin to the intersection point
            // If no intersection, 
            if ( vertex_ ) {
                ray.org = vertex.position;
            }
            else {
                return make_zero_spectrum();
            }
        }

        // if ( bounces > 0 )
        //     return ray.org;
        
        // std::cout << bounces << std::endl;
        
        // if ( x == 200 && y == 200 ) {
        //     std::cout << bounces << std::endl;
        // }



        current_path_throughput *= transmittance / average(trans_dir_pdf);


        // Hit a light source.
        // Add light contribution.
        if ( !scatter && vertex_ && is_light(scene.shapes[vertex.shape_id]) ) {
            // reach a surface, include emission

            Spectrum Le = emission(vertex, -ray.dir, scene);

            if ( bounces == 0 ) {

                // This is the only way we can see the light source, so
                // we don’t need multiple importance sampling.
                radiance += current_path_throughput * Le;

                return radiance;

            } else {
                // Need to account for next event estimation

        
                // Add light contribution only if the surface is emissive

                // Compute the probability of sampling the current intersected light point
                // from the point the next event estimation was issued 
                PointAndNormal light_point = PointAndNormal{vertex.position, vertex.geometry_normal};
                int light_id = get_area_light_id(scene.shapes[vertex.shape_id]);
                assert(light_id >= 0);
                const Light &light = scene.lights[light_id];
                
                // Compute the pdf of the nee only when at least one nee has been issued
                // TODO: use is_nee_issued condition if error
                Spectrum pdf_nee = light_pmf(scene, light_id) *
                                   pdf_point_on_light(light, light_point, nee_p_cache, scene) * 
                                   trans_nee_pdf;
                

                // Next, compute the PDF for sampling the current intersected light point
                // using the latest phase function sampling + the all trasmittance sampling
                // after the last phase function sampling.


                // The geometry term (=jacobian)
                Real jacobian = max(-dot(-ray.dir, light_point.normal), Real(0)) /
                                distance_squared(nee_p_cache, light_point.position);

                Spectrum pdf_phase = dir_pdf * multi_trans_pdf * jacobian;


                // Compute the multi importance sampling between
                // the next event estimation and the phase function sampling
                Spectrum w = (pdf_phase * pdf_phase) / (pdf_phase * pdf_phase + pdf_nee * pdf_nee);

                // Add the emission weighted by the multi importance sampling
                radiance += current_path_throughput * Le * w;

            }
            
        }



        // Hit a index-matchcing surface
        if ( !scatter && vertex_ && vertex.material_id == -1 ) {
            // If the intersected surface is a index-matching surface
            // update the current medium index and 
            // pass through without scattering
            current_medium = update_medium(vertex, ray, current_medium);

            ray.org = vertex.position;

            bounces++;
            continue;
        }


        // Reached the maximum bounces. Terminate.
        if ( bounces >= scene.options.max_depth - 1 && scene.options.max_depth != -1 ) {
            break;
        }



        if ( scatter && current_medium != -1 ) {
            

            const Medium &media = scene.media[current_medium];
            Spectrum sigma_s = get_sigma_s(media, ray.org);

            Spectrum nee = next_event_estimation_final(scene, 
                                                rng, 
                                                ray.org, 
                                                current_medium, 
                                                bounces, 
                                                -ray.dir,
                                                false,
                                                vertex);

            radiance += current_path_throughput * sigma_s * nee;


            // Record the last position that can issue a next event estimation
            // NEE is 0 and invalid if it is blocked by something
            // or does not reach the surface before the bounce limit
            if ( max(nee) > 0 ) {
                nee_p_cache = ray.org;
                is_nee_issued = true;
            }


            // return current_path_throughput * Spectrum(0,0,1);
            // sample next direct & update path throughput
            PhaseFunction phase_function = get_phase_function(media);
            const Vector2 phase_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
            std::optional<Vector3> next_dir_ = sample_phase_function(phase_function, -ray.dir, phase_uv);
            Vector3 &next_dir = *next_dir_;


            Real phase_pdf = pdf_sample_phase(phase_function, -ray.dir, next_dir);

            current_path_throughput *= eval(phase_function, -ray.dir, next_dir) /
                                       phase_pdf *
                                       sigma_s;

            ray.dir = next_dir;

            // Record the pdf of the latest phase function sampling for importance sampling
            // Also reset the accumulation buffer of the transmittance pdf
            // since the last phase function sampling
            dir_pdf = phase_pdf;
            multi_trans_pdf = make_const_spectrum(1); 
        }
        else if ( vertex_ ) {
            // Hit a surface

            Spectrum nee = next_event_estimation_final(scene, 
                                                rng, 
                                                ray.org, 
                                                current_medium, 
                                                bounces, 
                                                -ray.dir,
                                                true,
                                                vertex);

            // After computing bsdf contribution?
            radiance += current_path_throughput * nee;



            // Record the last position that can issue a next event estimation
            // NEE is 0 and invalid if it is blocked by something
            // or does not reach the surface before the bounce limit
            if ( max(nee) > 0 ) {
                nee_p_cache = ray.org;
                is_nee_issued = true;
            }




            const Material &mat = scene.materials[vertex.material_id];
            
            Vector3 dir_view = -ray.dir;
            Vector2 bsdf_rnd_param_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
            Real bsdf_rnd_param_w = next_pcg32_real<Real>(rng);
            std::optional<BSDFSampleRecord> bsdf_sample_ =
                sample_bsdf(mat,
                            dir_view,
                            vertex,
                            scene.texture_pool,
                            bsdf_rnd_param_uv,
                            bsdf_rnd_param_w);

            if (!bsdf_sample_) {
                // BSDF sampling failed. Abort the loop.
                break;
            }

            const BSDFSampleRecord &bsdf_sample = *bsdf_sample_;
            Vector3 dir_bsdf = bsdf_sample.dir_out;
            ray.dir = dir_bsdf;

            if (bsdf_sample.eta == 0) {
                ray_diff.spread = reflect(ray_diff, vertex.mean_curvature, bsdf_sample.roughness);
            } else {
                ray_diff.spread = refract(ray_diff, vertex.mean_curvature, bsdf_sample.eta, bsdf_sample.roughness);
                eta_scale /= (bsdf_sample.eta * bsdf_sample.eta);
                current_medium = update_medium(vertex, ray, current_medium);
            }

            Spectrum f = eval(mat, dir_view, dir_bsdf, vertex, scene.texture_pool);
            Real pfd_bsdf = pdf_sample_bsdf(mat, dir_view, dir_bsdf, vertex, scene.texture_pool);

            current_path_throughput *= f / pfd_bsdf;

        }


        Real rr_prob = Real(1);
        if ( bounces >= scene.options.rr_depth ) {
            // TODO: choosing R channel for now
            rr_prob = min(max((1 / eta_scale) * current_path_throughput), Real(0.95));
            Real u = next_pcg32_real<Real>(rng);
            if ( u > rr_prob ) {
                break;
            }
            else {
                current_path_throughput /= rr_prob;
            }
        }

        bounces++;

    }

    return radiance;
}


// The third volumetric renderer (not so simple anymore): 
// multiple monochromatic homogeneous volumes with multiple scattering
// no need to handle surface lighting, only directly visible light source
Spectrum vol_path_tracing_3(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!

    return vol_path_tracing(scene, x, y, rng);

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


        if ( current_medium != -1 ) {
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


        if ( !scatter && vertex_ ) {
            if ( vertex.material_id == -1 ) {
                current_medium = update_medium(vertex, ray, current_medium);
                bounces++;
                continue;
            }
        }

        if ( bounces == scene.options.max_depth - 1 && scene.options.max_depth != -1 ) {
            break;
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

    return vol_path_tracing(scene, x, y, rng);


    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                      (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    ray.tnear = get_intersection_epsilon(scene);
    RayDifferential ray_diff = RayDifferential{Real(0), Real(0)};

    int current_medium = scene.camera.medium_id;
    Spectrum current_path_throughput = make_const_spectrum(1);
    Spectrum radiance = make_zero_spectrum();
    int bounces = 0;

    Real dir_pdf = 0; // the pdf of the latest phase function sampling
    Vector3 nee_p_cache; // the last position p that can issue a next event estimation

    // flag to record if nee is issued to avoid including nee pdf contribution
    // when no nee has been issued
    bool is_nee_issued = false;
    
    // The product PDF of transmittance sampling going through 
    // several index-matching surfaces from the last phase function sampling
    Spectrum multi_trans_pdf = make_const_spectrum(1); 

    while ( true ) {
        
        bool scatter = false;

        std::optional<PathVertex> vertex_ = intersect(scene, ray, ray_diff);
        PathVertex vertex = *vertex_;

        Spectrum transmittance = make_const_spectrum(1);
        Spectrum trans_pdf = make_const_spectrum(1);


        if ( current_medium != -1 ) {
            
            const Medium media = scene.media[current_medium];
            Spectrum sigma_s = get_sigma_s(media, vertex.position);
            Spectrum sigma_a = get_sigma_a(media, vertex.position);
            Spectrum sigma_t = sigma_s + sigma_a;
            
            Real u = next_pcg32_real<Real>(rng);
            Real t = -log(1 - u) / sigma_t.x;

            trans_pdf = exp(-sigma_t * t) * sigma_t;
            transmittance = exp(-sigma_t * t);

            // Accumulate the pdf of the trasmittance sampling
            // since the previous phase function sampling
            multi_trans_pdf *= trans_pdf;

            Real t_hit = distance(vertex.position, ray.org);
            if ( !vertex_ ) t_hit = Real(INFINITY);

            if ( t < t_hit ) {
                scatter = true;

                Vector3 p = ray.org + t * ray.dir;

                ray.org = ray.org + t * ray.dir;   
            }
            else {
                ray.org = vertex.position;
            }

        }
        else {
            // For vaccume volume, we simply move the ray origin to the intersection point
            // If no intersection, 
            if ( vertex_ ) {
                ray.org = vertex.position;
            }
            else {
                return make_zero_spectrum();
            }
        }


        current_path_throughput *= transmittance / trans_pdf;


        // Hit a light source.
        // Add light contribution.
        if ( !scatter && vertex_ && is_light(scene.shapes[vertex.shape_id]) ) {
            // reach a surface, include emission

            Spectrum Le = emission(vertex, -ray.dir, scene);

            if ( bounces == 0 ) {

                // This is the only way we can see the light source, so
                // we don’t need multiple importance sampling.
                radiance += current_path_throughput * Le;

            } else {
                // Need to account for next event estimation

        
                // Add light contribution only if the surface is emissive

                // Compute the probability of sampling the current intersected light point
                // from the point the next event estimation was issued 
                PointAndNormal light_point = PointAndNormal{vertex.position, vertex.geometry_normal};
                int light_id = get_area_light_id(scene.shapes[vertex.shape_id]);
                assert(light_id >= 0);
                const Light &light = scene.lights[light_id];
                
                // Compute the pdf of the nee only when at least one nee has been issued
                Real pdf_nee = Real(0);
                if ( is_nee_issued ) {
                    // Need to add light_pmf(scene, vertex.shape_id) if
                    // there is more than 1 light source
                    // Otherwise, adding it causes an assertion error
                    pdf_nee = pdf_point_on_light(light, light_point, nee_p_cache, scene);
                }
                

                // Next, compute the PDF for sampling the current intersected light point
                // using the latest phase function sampling + the all trasmittance sampling
                // after the last phase function sampling.


                // The geometry term (=jacobian)
                Real jacobian = fabs(dot(-ray.dir, light_point.normal)) / 
                                distance_squared(nee_p_cache, light_point.position);

                Spectrum pdf_phase = dir_pdf * multi_trans_pdf * jacobian;


                // Compute the multi importance sampling between
                // the next event estimation and the phase function sampling
                Spectrum w = (pdf_phase * pdf_phase) / (pdf_phase * pdf_phase + pdf_nee * pdf_nee);

                // Add the emission weighted by the multi importance sampling
                radiance += current_path_throughput * Le * w;

        

            }
            
        }

        // Reached the maximum bounces. Terminate.
        if ( bounces == scene.options.max_depth - 1 && scene.options.max_depth != -1 ) {
            break;
        }


        // Hit a index-matchcing surface
        if ( !scatter && vertex_ ) {
            // If the intersected surface is a index-matching surface
            // update the current medium index and 
            // pass through without scattering
            if ( vertex.material_id == -1 ) {
                current_medium = update_medium(vertex, ray, current_medium);

                bounces++;
                continue;
            }
        }


        if ( scatter && current_medium != -1 ) {

            const Medium media = scene.media[current_medium];
            Spectrum sigma_s = get_sigma_s(media, vertex.position);

            Spectrum nee = next_event_estimation(scene, 
                                                rng, 
                                                ray.org, 
                                                current_medium, 
                                                bounces, 
                                                -ray.dir,
                                                vertex);

            // return transmittance * sigma_s * nee / trans_pdf;
            radiance += transmittance * sigma_s * nee / trans_pdf;

            // Record the last position that can issue a next event estimation
            // NEE is 0 and invalid if it is blocked by something
            // or does not reach the surface before the bounce limit
            if ( nee.x > 0 ) {
                nee_p_cache = ray.org;
                is_nee_issued = true;
            }


            // return current_path_throughput * Spectrum(0,0,1);
            // sample next direct & update path throughput
            PhaseFunction phase_function = get_phase_function(media);
            const Vector2 phase_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
            std::optional<Vector3> next_dir_ = sample_phase_function(phase_function, -ray.dir, phase_uv);
            Vector3 &next_dir = *next_dir_;


            Real phase_pdf = pdf_sample_phase(phase_function, -ray.dir, next_dir);

            current_path_throughput *= eval(phase_function, -ray.dir, next_dir) /
                                       phase_pdf *
                                       sigma_s;

            ray.dir = next_dir;

            // Record the pdf of the latest phase function sampling for importance sampling
            // Also reset the accumulation buffer of the transmittance pdf
            // since the last phase function sampling
            dir_pdf = phase_pdf;
            multi_trans_pdf = make_const_spectrum(1); 
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

// The fifth volumetric renderer: 
// multiple monochromatic homogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// with surface lighting
Spectrum vol_path_tracing_5(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!

    return vol_path_tracing(scene, x, y, rng);


    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                      (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    // ray.tnear = get_intersection_epsilon(scene);
    RayDifferential ray_diff = RayDifferential{Real(0), Real(0)};

    int current_medium = scene.camera.medium_id;
    Spectrum current_path_throughput = make_const_spectrum(1);
    Spectrum radiance = make_zero_spectrum();
    int bounces = 0;

    Real dir_pdf = 0; // the pdf of the latest phase function sampling
    Vector3 nee_p_cache; // the last position p that can issue a next event estimation

    // flag to record if nee is issued to avoid including nee pdf contribution
    // when no nee has been issued
    bool is_nee_issued = false;
    
    // The product PDF of transmittance sampling going through 
    // several index-matching surfaces from the last phase function sampling
    Spectrum multi_trans_pdf = make_const_spectrum(1); 

    Real eta_scale = Real(1);

    while ( true ) {
        
        bool scatter = false;

        std::optional<PathVertex> vertex_ = intersect(scene, ray, ray_diff);
        PathVertex vertex = *vertex_;

        Spectrum transmittance = make_const_spectrum(1);
        Spectrum trans_pdf = make_const_spectrum(1);


        if ( current_medium != -1 ) {
            
            const Medium media = scene.media[current_medium];
            Spectrum sigma_s = get_sigma_s(media, vertex.position);
            Spectrum sigma_a = get_sigma_a(media, vertex.position);
            Spectrum sigma_t = sigma_s + sigma_a;
            
            Real u = next_pcg32_real<Real>(rng);
            Real t = -log(1 - u) / sigma_t.x;

            trans_pdf = exp(-sigma_t * t) * sigma_t;
            transmittance = exp(-sigma_t * t);

            // Accumulate the pdf of the trasmittance sampling
            // since the previous phase function sampling
            multi_trans_pdf *= trans_pdf;

            Real t_hit = distance(vertex.position, ray.org);
            if ( !vertex_ ) t_hit = Real(INFINITY);

            if ( t < t_hit ) {
                scatter = true;

                Vector3 p = ray.org + t * ray.dir;

                ray.org = ray.org + t * ray.dir;   
            }
            else {
                ray.org = vertex.position;
            }

        }
        else {
            // For vaccume volume, we simply move the ray origin to the intersection point
            // If no intersection, 
            if ( vertex_ ) {
                ray.org = vertex.position;
            }
            else {
                return make_zero_spectrum();
            }
        }


        current_path_throughput *= transmittance / trans_pdf;


        // Hit a light source.
        // Add light contribution.
        if ( !scatter && vertex_ && is_light(scene.shapes[vertex.shape_id]) ) {
            // reach a surface, include emission

            Spectrum Le = emission(vertex, -ray.dir, scene);

            if ( bounces == 0 ) {

                // This is the only way we can see the light source, so
                // we don’t need multiple importance sampling.
                radiance += current_path_throughput * Le;

            } else {
                // Need to account for next event estimation

        
                // Add light contribution only if the surface is emissive

                // Compute the probability of sampling the current intersected light point
                // from the point the next event estimation was issued 
                PointAndNormal light_point = PointAndNormal{vertex.position, vertex.geometry_normal};
                int light_id = get_area_light_id(scene.shapes[vertex.shape_id]);
                assert(light_id >= 0);
                const Light &light = scene.lights[light_id];
                
                // Compute the pdf of the nee only when at least one nee has been issued
                Real pdf_nee = Real(0);
                if ( is_nee_issued ) {
                    // Need to add light_pmf(scene, vertex.shape_id) if
                    // there is more than 1 light source
                    // Otherwise, adding it causes an assertion error
                    pdf_nee = pdf_point_on_light(light, light_point, nee_p_cache, scene);
                }
                

                // Next, compute the PDF for sampling the current intersected light point
                // using the latest phase function sampling + the all trasmittance sampling
                // after the last phase function sampling.


                // The geometry term (=jacobian)
                Real jacobian = fabs(dot(-ray.dir, light_point.normal)) / 
                                distance_squared(nee_p_cache, light_point.position);

                Spectrum pdf_phase = dir_pdf * multi_trans_pdf * jacobian;


                // Compute the multi importance sampling between
                // the next event estimation and the phase function sampling
                Spectrum w = (pdf_phase * pdf_phase) / (pdf_phase * pdf_phase + pdf_nee * pdf_nee);

                // Add the emission weighted by the multi importance sampling
                radiance += current_path_throughput * Le * w;

        

            }
            
        }

        // Reached the maximum bounces. Terminate.
        if ( bounces == scene.options.max_depth - 1 && scene.options.max_depth != -1 ) {
            break;
        }


        // Hit a index-matchcing surface
        if ( !scatter && vertex_ ) {
            // If the intersected surface is a index-matching surface
            // update the current medium index and 
            // pass through without scattering
            if ( vertex.material_id == -1 ) {
                current_medium = update_medium(vertex, ray, current_medium);

                // return Spectrum(0,1,0);

                bounces++;
                continue;
            }
        }


        if ( scatter && current_medium != -1 ) {

            // if ( vertex_ ) return Spectrum(0,10,0);
            // else return Spectrum(0,0,1);
            

            const Medium media = scene.media[current_medium];
            Spectrum sigma_s = get_sigma_s(media, vertex.position);

            Spectrum nee = next_event_estimation(scene, 
                                                rng, 
                                                ray.org, 
                                                current_medium, 
                                                bounces, 
                                                -ray.dir,
                                                vertex);

            // return transmittance * sigma_s * nee / trans_pdf;
            radiance += current_path_throughput * sigma_s * nee;

            // Record the last position that can issue a next event estimation
            // NEE is 0 and invalid if it is blocked by something
            // or does not reach the surface before the bounce limit
            if ( nee.x > 0 ) {
                nee_p_cache = ray.org;
                is_nee_issued = true;
            }


            // return current_path_throughput * Spectrum(0,0,1);
            // sample next direct & update path throughput
            PhaseFunction phase_function = get_phase_function(media);
            const Vector2 phase_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
            std::optional<Vector3> next_dir_ = sample_phase_function(phase_function, -ray.dir, phase_uv);
            Vector3 &next_dir = *next_dir_;


            Real phase_pdf = pdf_sample_phase(phase_function, -ray.dir, next_dir);

            current_path_throughput *= eval(phase_function, -ray.dir, next_dir) /
                                       phase_pdf *
                                       sigma_s;

            ray.dir = next_dir;

            // Record the pdf of the latest phase function sampling for importance sampling
            // Also reset the accumulation buffer of the transmittance pdf
            // since the last phase function sampling
            dir_pdf = phase_pdf;
            multi_trans_pdf = make_const_spectrum(1); 
        }
        else if ( vertex_ ) {
            // Hit a surface

            Spectrum nee = next_event_estimation(scene, 
                                                rng, 
                                                ray.org, 
                                                current_medium, 
                                                bounces, 
                                                -ray.dir,
                                                vertex);
                                                

            // Record the last position that can issue a next event estimation
            // NEE is 0 and invalid if it is blocked by something
            // or does not reach the surface before the bounce limit
            if ( max(nee) > 0 ) {
                nee_p_cache = ray.org;
                is_nee_issued = true;
            }


            const Material &mat = scene.materials[vertex.material_id];
            

            Vector3 dir_view = -ray.dir;
            Vector2 bsdf_rnd_param_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
            Real bsdf_rnd_param_w = next_pcg32_real<Real>(rng);
            std::optional<BSDFSampleRecord> bsdf_sample_ =
                sample_bsdf(mat,
                            dir_view,
                            vertex,
                            scene.texture_pool,
                            bsdf_rnd_param_uv,
                            bsdf_rnd_param_w);

            if (!bsdf_sample_) {
                // BSDF sampling failed. Abort the loop.
                break;
            }

            const BSDFSampleRecord &bsdf_sample = *bsdf_sample_;
            Vector3 dir_bsdf = bsdf_sample.dir_out;
            ray.dir = dir_bsdf;

            if (bsdf_sample.eta == 0) {
                ray_diff.spread = reflect(ray_diff, vertex.mean_curvature, bsdf_sample.roughness);
            } else {
                ray_diff.spread = refract(ray_diff, vertex.mean_curvature, bsdf_sample.eta, bsdf_sample.roughness);
                eta_scale /= (bsdf_sample.eta * bsdf_sample.eta);
                current_medium = update_medium(vertex, ray, current_medium);
            }

            Spectrum f = eval(mat, dir_view, dir_bsdf, vertex, scene.texture_pool);
            Real pfd_bsdf = pdf_sample_bsdf(mat, dir_view, dir_bsdf, vertex, scene.texture_pool);

            current_path_throughput *= f / pfd_bsdf;

            // After computing bsdf contribution?
            radiance += current_path_throughput * nee;



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
                current_path_throughput /= rr_prob;
            }
        }

        bounces++;
    }

    return radiance;
}