#include "../microfacet.h"

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

    // std::cout << vertex.geometry_normal << std::endl;

    Real u = vertex.st.x; // u=[0,1] along the curve
    Real v = vertex.st.y; // v=[0,1] along the width of the tube, should be 0 for RTC_GEOMETRY_TYPE_ROUND_BEZIER_CURVE but not
    Real h = Real(-1) + Real(2) * v; // h=[-1, 1] point on the diameter of the intersection circle

    static const int pMax = 3;

    std::cout << bsdf.beta_m << std::endl;





    return Spectrum(vertex.geometry_normal[0] / Real(2) + Real(0.5), vertex.geometry_normal[1] / Real(2) + Real(0.5), vertex.geometry_normal[2] / Real(2) + Real(0.5));



    return make_zero_spectrum();
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

    // return 1;

    return 0;
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const HairBSDF &bsdf) const {
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometry_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!

    return {};
}

TextureSpectrum get_texture_op::operator()(const HairBSDF &bsdf) const {
    return {};
}
