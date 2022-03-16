#pragma once

#include "lajolla.h"
#include "shape.h"

/// Load Mitsuba's serialized file format.
BezierCurve load_hair(const fs::path &filename, float radius);