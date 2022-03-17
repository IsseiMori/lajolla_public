uint32_t register_embree_op::operator()(const BezierCurve &curve) const {
    // RTC_GEOMETRY_TYPE_FLAT_BEZIER_CURVE
    // RTC_GEOMETRY_TYPE_NORMAL_ORIENTED_BEZIER_CURVE
    // RTC_GEOMETRY_TYPE_ROUND_BEZIER_CURVE
    // https://manpages.debian.org/experimental/libembree-dev/RTC_GEOMETRY_TYPE_CURVE.3.en.html

    RTCGeometry rtc_geom = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_ROUND_BEZIER_CURVE);
    // A geomID is the ID associated with the shape inside Embree.
    uint32_t geomID = rtcAttachGeometry(scene, rtc_geom);
    Vector4f *points = (Vector4f*)rtcSetNewGeometryBuffer(
        rtc_geom, RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT4,
        sizeof(Vector4f), curve.points.size());
    
    uint32_t *indices = (uint32_t*)rtcSetNewGeometryBuffer(
        rtc_geom, RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT,
        sizeof(uint32_t), curve.num);
    
    // vertex
    for (int i = 0; i < curve.points.size(); i++) {
        Vector3 point = curve.points[i];
        points[i] = Vector4f{(float)point[0], (float)point[1], (float)point[2], curve.radius};
    }
    
    // indices pointing to first control vertex
    for (int i = 0; i < curve.num; i++) {
        indices[i] = (uint32_t)(i * 4);
    }

    rtcSetGeometryVertexAttributeCount(rtc_geom, 1);
    rtcCommitGeometry(rtc_geom);
    rtcReleaseGeometry(rtc_geom);
    return geomID;
}

PointAndNormal sample_point_on_shape_op::operator()(const BezierCurve &curve) const {
    
}

Real surface_area_op::operator()(const BezierCurve &curve) const {
    
}

Real pdf_point_on_shape_op::operator()(const BezierCurve &curve) const {
    
}

void init_sampling_dist_op::operator()(BezierCurve &curve) const {
   
}

ShadingInfo compute_shading_info_op::operator()(const BezierCurve &curve) const {
    
}

