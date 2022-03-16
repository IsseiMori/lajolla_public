uint32_t register_embree_op::operator()(const BezierCurve &curve) const {
    // RTC_GEOMETRY_TYPE_FLAT_BEZIER_CURVE
    // RTC_GEOMETRY_TYPE_NORMAL_ORIENTED_BEZIER_CURVE
    // RTC_GEOMETRY_TYPE_ROUND_BEZIER_CURVE
    RTCGeometry rtc_geom = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_ROUND_BEZIER_CURVE);
    // A geomID is the ID associated with the shape inside Embree.
    uint32_t geomID = rtcAttachGeometry(scene, rtc_geom);
    Vector4f *points = (Vector4f*)rtcSetNewGeometryBuffer(
        rtc_geom, RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT3,
        sizeof(Vector4f), curve.points.size());
    
    Vector3i *curves = (Vector3i*)rtcSetNewGeometryBuffer(
        rtc_geom, RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT3,
        sizeof(Vector3i), mesh.indices.size());
    
    // only vertex
    for (int i = 0; i < curve.points.size(); i++) {
        Vector3 point = curve.points[i];
        points[i] = Vector4f{(float)point[0], (float)point[1], (float)point[2], 0.f};
    }

    // dummy indices
    for (int i = 0; i < (int)curve.points.size(); i++) {
        curves[i] = i;
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

