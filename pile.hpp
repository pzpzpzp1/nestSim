#ifndef PILE_HPP
#define PILE_HPP

// Structures/functions to deal with metrics of a pile
//  All indexed midpoints
//  All indexed orientations
//  All double-indexed contact points
//  All indexed stability counts
//  xyz mass distribution

void printFloatVector(std::vector<float> vec) {
    printf("[");
    for (int i = 0; i < vec.size() - 1; i++) {
        printf("%f, ", vec[i]);
    }
    printf("%f]\n", vec[vec.size() - 1]);
}


float highestMidpoint(void) {
    float height = 0;

    for (int i = 0; i < obj.size(); i++) {
        const dReal* pos = dGeomGetPosition(obj[i].geom);
        if (pos[2] > height) { height = pos[2]; }
    }

    return height;
}


// Global packing fraction
//    Set height = 0 to use highest midpoint as the height
float packingFraction(float height) {
    float total_vol, obj_vol;

    if (nwalls < 5) {
        printf("Can't calculate packing fraction without boundaries.\n");
        return 0;
    }

    if (height == 0) { height = highestMidpoint(); }

    // Calculate the volume of the box/cylinder
    if (nwalls == 5) {
        total_vol = height * pow(bound * 2, 2);
    }
    else {
        total_vol = height * M_PI * pow(bound, 2);
    }

    // Calculate the volume of the rods
    obj_vol = obj.size() * ((M_PI * pow(rad, 2)) * (rad * AR)
                            + (4/3) * M_PI * pow(rad, 3));

    return obj_vol / total_vol;
}


// No cotan function in C++!
float cot(float theta) {
    return cos(theta) / sin(theta);
}


// The indefinite integral of b * sqrt(1 - x^2 / a^2), which defines an
//      ellipse of semi-major length a and semi-minor length b
// x = point to evaluate the integral
// a = semi-major axis of ellipse = rad / cos(theta)
// b = semi-minor axis of ellipse
float ellipseIntegral(float x, float a, float b) {
    assert(std::abs(x) <= a);

    float res = 0.;

    if (a >= 0) {
        res += b * (x * sqrt(1 - pow(x / a, 2)) + a * asin(x / a));
    } else {
        a = -a;
        res -= b * (x * sqrt(1 - pow(x / a, 2)) + a * asin(x / a));
    }

    assert(!isnan(res));

//    printf("\t\t(x, a, b, res): (%f, %f, %f, %f)\n", x, a, b, res);

    return res;
}


// Calculate the mass density for height-parameterized slices dh-thick.
// Based on area-density of a slice; no integration in the z-direction.
// Because of this, make sure that dh < rad.
std::vector<float> massDensityByHeight(float dh) {
    bool log = false;

    int nslices = ceil(highestMidpoint() / dh);
    std::vector<float> density(nslices, 0);

    // Area of a slice
    float slice_area;
    if (nwalls == 5) {
        slice_area = pow(bound * 2, 2);
    }
    else if (nwalls > 5) {
        slice_area = M_PI * pow(bound, 2);
    }
    else {
        slice_area = 1.;
        // throw std::runtime_error("Can't take slices without boundaries.");
    }

    // Iterate through the objects
    for (int i = 0; i < obj.size(); i++) {
        if (log) { printf("Rod %d:\n", i); }

        std::vector<float> pos = position(i);
        std::vector<float> angles = orientationAngles(i);

        if (isnan(angles[0])) {
            throw std::runtime_error("NaN angle encountered.");
            continue;
        }

//        printf("orientation angles: (%f, %f)\n", angles[0], angles[1]);

//        printf("sin(theta): %f\n", std::abs(sin(angles[0])));
//        printf("offset: %f\n", rad * AR * sin(angles[0]));

//        printf("\tzpos: %f\n", pos[2]);

        float h_high = pos[2] + rad + (0.5 * rad * AR * std::abs(cos(angles[0])));
        float h_low = pos[2] - rad - (0.5 * rad * AR * std::abs(cos(angles[0])));

        // Sometimes the rods end up "intersecting with the floor". Push
        // them above the floor.
        if (h_low < 0) {
            h_high = h_high + h_low;
            h_low = 0;
        }

//        printf("nominal length: %f\n", (h_high - h_low - (2 * rad)) / cos(angles[0]));

//        printf("\tlower, higher: [%f, %f]\n", h_low + rad, h_high - rad);
//        printf("\tlowest, highest: [%f, %f]\n", h_low, h_high);

        float mid_of_upper_cap = h_high - rad;
        float top_of_upper_cap = h_high - rad + rad * std::abs(sin(angles[0]));
        float bot_of_upper_cap = h_high - rad - rad * std::abs(sin(angles[0]));

        float mid_of_lower_cap = h_low + rad;
        float top_of_lower_cap = h_low + rad + rad * std::abs(sin(angles[0]));
        float bot_of_lower_cap = h_low + rad - rad * std::abs(sin(angles[0]));

//        printf("upper cap: [%f, %f, %f]\n", bot_of_upper_cap, mid_of_upper_cap, top_of_upper_cap);
//        printf("lower cap: [%f, %f, %f]\n\n", bot_of_lower_cap, mid_of_lower_cap, top_of_lower_cap);

        assert(top_of_upper_cap >= bot_of_upper_cap);
        assert(top_of_lower_cap >= bot_of_lower_cap);

        // Highest slice may be higher than the nominal pile height, currently
        // defined by the highest midpoint.
        int lowest_slice = floor(h_low / dh);
        int highest_slice = floor(h_high / dh);
        if (highest_slice >= nslices) {
            highest_slice = nslices - 1;
        }

//        printf("low/high slice of %d: %d/%d\n", nslices, lowest_slice, highest_slice);

        // The end of the rod is below zero for some reason
        if (lowest_slice < 0) {
            throw std::runtime_error("Slice too low.");
        }

        // Prepare for iterating through planes
        float h0 = dh / 2. + EPSILON;
        float h, area, offset;
        float a, offset_low, offset_high, d_expr, d_low, d_high,
              e_low, e_high, n_low, n_high;
        bool contained;

        // Iterate through slices
        for (int j = lowest_slice; j <= highest_slice; j++) {
            h = h0 + j * dh;
            area = 0.;

            // Some precalculated values
            a             = rad / cos(angles[0]);
            offset_low    = std::abs(h - mid_of_lower_cap);
            offset_high   = std::abs(h - mid_of_upper_cap);
            d_expr        = cot(angles[0]) + 1. / cot(angles[0]);
            d_low         = offset_low * d_expr;
            d_high        = offset_high * d_expr;
            e_low         = std::abs(h - bot_of_lower_cap) * d_expr;
            e_high        = std::abs(top_of_upper_cap - h) * d_expr;
            n_low         = rad * sqrt(1 - pow(offset_low / rad, 2));
            n_high        = rad * sqrt(1 - pow(offset_high / rad, 2));

            // Intersection code
            contained = (h_low <= h && h <= h_high);

            // The rod is essentially horizontal (w/in fraction of a degree)
            if (contained && std::abs(M_PI/2 - angles[0]) < EPSILON) {
                if (log) { printf("\thorizontal          \t"); }
                offset = std::abs(h - pos[2]);
                area = 2 * rad * sqrt(1. - pow(offset / rad, 2))
                         * (rad * AR)
                       + 2 * M_PI * rad * sqrt(1. - pow(offset / rad, 2));
            }

            // Rectangular only
            else if (contained
                     && h >= bot_of_lower_cap
                     && h >= top_of_lower_cap
                     && h <= bot_of_upper_cap
                     && h <= top_of_upper_cap) {
                if (log) { printf("\trectangular only\t"); }
                area = M_PI * rad * a; // elliptical area
            }

            // Bottom cap only
            else if (contained
                     && h <= bot_of_lower_cap
                     && h <= top_of_lower_cap
                     && h <= bot_of_upper_cap
                     && h <= top_of_upper_cap) {
                if (log) { printf("\tbottom cap only    \t"); }
                area = M_PI * pow(n_low, 2);
            }

            // Top cap only
            else if (contained
                     && h >= bot_of_lower_cap
                     && h >= top_of_lower_cap
                     && h >= bot_of_upper_cap
                     && h >= top_of_upper_cap) {
                if (log) { printf("\ttop cap only        \t"); }
                area = M_PI * pow(n_high, 2);
            }

            // Bottom cap + rectangular
            else if (contained
                     && h >= bot_of_lower_cap
                     && h <= top_of_lower_cap
                     && h <= bot_of_upper_cap
                     && h <= top_of_upper_cap) {
                if (log) { printf("\tbot + rectangular"); }
                if (h < mid_of_lower_cap) {
                    if (log) { printf("*"); }
                    area =    ellipseIntegral(a, a, rad)
                            - ellipseIntegral(a - e_low, a, rad)
                            + ellipseIntegral(n_low, n_low, n_low)
                            - ellipseIntegral(-n_low * cos(angles[0]), n_low, n_low);
                } else {
                    area =    ellipseIntegral(a, a, rad)
                            - ellipseIntegral(-d_low, a, rad)
                            + ellipseIntegral(n_low, n_low, n_low)
                            - ellipseIntegral(n_low * cos(angles[0]), n_low, n_low);
                }
                if (log) { printf("\t"); }
            }

            // Both caps
            else if (contained
                     && h >= bot_of_lower_cap
                     && h <= top_of_lower_cap
                     && h >= bot_of_upper_cap
                     && h <= top_of_upper_cap) {
                if (log) { printf("\tboth caps        \t"); }

                // Higher than both midpoints
                if (h > mid_of_upper_cap) {
                    // rectangular region, lower cap, upper cap
                    area =    ellipseIntegral(d_high + rad * AR / sin(angles[0]), a, rad)
                            - ellipseIntegral(d_high, a, rad)
                            + ellipseIntegral(n_low, n_low, n_low)
                            - ellipseIntegral(n_low * cos(angles[0]), n_low, n_low)
                            + ellipseIntegral(n_high, n_high, n_high)
                            - ellipseIntegral(-n_high * cos(angles[0]), n_high, n_high);
                }
                // Higher than the lower midpoint, lower than the higher midpoint
                else if (mid_of_lower_cap <= h && h <= mid_of_upper_cap) {
                    // rectangular region, lower cap, upper cap
                    area =    ellipseIntegral(d_high, a, rad)
                            - ellipseIntegral(-d_low, a, rad)
                            + ellipseIntegral(n_low, n_low, n_low)
                            - ellipseIntegral(n_low * cos(angles[0]), n_low, n_low)
                            + ellipseIntegral(n_high, n_high, n_high)
                            - ellipseIntegral(-n_high * cos(angles[0]), n_high, n_high);
                }
                // Lower than both midpoints
                else {
                    // rectangular region, lower cap, upper cap
                    area =    ellipseIntegral(d_low + rad * AR / sin(angles[0]), a, rad)
                            - ellipseIntegral(d_low, a, rad)
                            + ellipseIntegral(n_low, n_low, n_low)
                            - ellipseIntegral(-n_low * cos(angles[0]), n_low, n_low)
                            + ellipseIntegral(n_high, n_high, n_high)
                            - ellipseIntegral(n_high * cos(angles[0]), n_high, n_high);
                }
            }

            // Rectangular + top cap
            else if (contained
                     && h >= bot_of_lower_cap
                     && h >= top_of_lower_cap
                     && h >= bot_of_upper_cap
                     && h <= top_of_upper_cap) {
                if (log) { printf("\ttop + rectangular"); }
                if (h > mid_of_upper_cap) {
                    if (log) { printf("*"); }
                    area =  ellipseIntegral(a, a, rad)
                            - ellipseIntegral(a - e_high, a, rad)
                            + ellipseIntegral(n_high, n_high, n_high)
                            - ellipseIntegral(-n_high * cos(angles[0]), n_high, n_high);
                } else {
                    area =  ellipseIntegral(a, a, rad)
                            - ellipseIntegral(-d_high, a, rad)
                            + ellipseIntegral(n_high, n_high, n_high)
                            - ellipseIntegral(n_high * cos(angles[0]), n_high, n_high);
                }
                if (log) { printf("\t"); }
            }

            else {
                if (log) { printf("\tno intersection\n"); }
                continue;
            }

            if (isnan(area)) {
                throw std::runtime_error("NaN area encountered.");
                break;
            }

            if (log) { printf("\t%f\n", area); }

            // This is the maximum possible area
            if (area > (2 * rad) * (rad * AR) + M_PI * pow(rad, 2) + 100 * EPSILON) {
                throw std::runtime_error("Area too large. Check your calculations.");
            }

            density[j] += area;
        }
    }

    for (int i = 0; i < nslices; i++) {
        density[i] /= slice_area;
    }

    return density;
}


// Calculate orientation density as a function of height
// Return the mass-averaged orientation (in S^2) for a given slice.
// dh = spacing between slices
// FYI, this is almost exactly copy-pasted from massDensityByHeight.
std::vector<float> orientationDensityByHeight(float dh) {
    bool log = false;

    int nslices = ceil(highestMidpoint() / dh);
    std::vector<float> density(nslices, 0);
    std::vector<float> theta(nslices, 0);

    // Placeholder vars
    float theta_ = 0.;

    // Area of a slice
    float slice_area;
    if (nwalls == 5) {
        slice_area = pow(bound * 2, 2);
    }
    else if (nwalls > 5) {
        slice_area = M_PI * pow(bound, 2);
    }
    else {
        slice_area = 1.;
        // throw std::runtime_error("Can't take slices without boundaries.");
    }

    // Iterate through the objects
    for (int i = 0; i < obj.size(); i++) {
        if (log) { printf("Rod %d:\n", i); }

        std::vector<float> pos = position(i);
        std::vector<float> angles = orientationAngles(i);

        if (isnan(angles[0])) {
            throw std::runtime_error("NaN angle encountered.");
            continue;
        }

        theta_ = angles[0];

        //        printf("orientation angles: (%f, %f)\n", angles[0], angles[1]);

        //        printf("sin(theta): %f\n", std::abs(sin(angles[0])));
        //        printf("offset: %f\n", rad * AR * sin(angles[0]));

        //        printf("\tzpos: %f\n", pos[2]);

        float h_high = pos[2] + rad + (0.5 * rad * AR * std::abs(cos(angles[0])));
        float h_low = pos[2] - rad - (0.5 * rad * AR * std::abs(cos(angles[0])));

        // Sometimes the rods end up "intersecting with the floor". Push
        // them above the floor.
        if (h_low < 0) {
            h_high = h_high + h_low;
            h_low = 0;
        }

        //        printf("nominal length: %f\n", (h_high - h_low - (2 * rad)) / cos(angles[0]));

        //        printf("\tlower, higher: [%f, %f]\n", h_low + rad, h_high - rad);
        //        printf("\tlowest, highest: [%f, %f]\n", h_low, h_high);

        float mid_of_upper_cap = h_high - rad;
        float top_of_upper_cap = h_high - rad + rad * std::abs(sin(angles[0]));
        float bot_of_upper_cap = h_high - rad - rad * std::abs(sin(angles[0]));

        float mid_of_lower_cap = h_low + rad;
        float top_of_lower_cap = h_low + rad + rad * std::abs(sin(angles[0]));
        float bot_of_lower_cap = h_low + rad - rad * std::abs(sin(angles[0]));

        //        printf("upper cap: [%f, %f, %f]\n", bot_of_upper_cap, mid_of_upper_cap, top_of_upper_cap);
        //        printf("lower cap: [%f, %f, %f]\n\n", bot_of_lower_cap, mid_of_lower_cap, top_of_lower_cap);

        assert(top_of_upper_cap >= bot_of_upper_cap);
        assert(top_of_lower_cap >= bot_of_lower_cap);

        // Highest slice may be higher than the nominal pile height, currently
        // defined by the highest midpoint.
        int lowest_slice = floor(h_low / dh);
        int highest_slice = floor(h_high / dh);
        if (highest_slice >= nslices) {
            highest_slice = nslices - 1;
        }

        //        printf("low/high slice of %d: %d/%d\n", nslices, lowest_slice, highest_slice);

        // The end of the rod is below zero for some reason
        if (lowest_slice < 0) {
            throw std::runtime_error("Slice too low.");
        }

        // Prepare for iterating through planes
        float h0 = dh / 2. + EPSILON;
        float h, area, offset;
        float a, offset_low, offset_high, d_expr, d_low, d_high,
        e_low, e_high, n_low, n_high;
        bool contained;

        // Iterate through slices
        for (int j = lowest_slice; j <= highest_slice; j++) {
            h = h0 + j * dh;
            area = 0.;

            // Some precalculated values
            a             = rad / cos(angles[0]);
            offset_low    = std::abs(h - mid_of_lower_cap);
            offset_high   = std::abs(h - mid_of_upper_cap);
            d_expr        = cot(angles[0]) + 1. / cot(angles[0]);
            d_low         = offset_low * d_expr;
            d_high        = offset_high * d_expr;
            e_low         = std::abs(h - bot_of_lower_cap) * d_expr;
            e_high        = std::abs(top_of_upper_cap - h) * d_expr;
            n_low         = rad * sqrt(1 - pow(offset_low / rad, 2));
            n_high        = rad * sqrt(1 - pow(offset_high / rad, 2));

            // Intersection code
            contained = (h_low <= h && h <= h_high);

            // The rod is essentially horizontal (w/in fraction of a degree)
            if (contained && std::abs(M_PI/2 - angles[0]) < EPSILON) {
                if (log) { printf("\thorizontal          \t"); }
                offset = std::abs(h - pos[2]);
                area = 2 * rad * sqrt(1. - pow(offset / rad, 2))
                * (rad * AR)
                + 2 * M_PI * rad * sqrt(1. - pow(offset / rad, 2));
            }

            // Rectangular only
            else if (contained
                     && h >= bot_of_lower_cap
                     && h >= top_of_lower_cap
                     && h <= bot_of_upper_cap
                     && h <= top_of_upper_cap) {
                if (log) { printf("\trectangular only\t"); }
                area = M_PI * rad * a; // elliptical area
            }

            // Bottom cap only
            else if (contained
                     && h <= bot_of_lower_cap
                     && h <= top_of_lower_cap
                     && h <= bot_of_upper_cap
                     && h <= top_of_upper_cap) {
                if (log) { printf("\tbottom cap only    \t"); }
                area = M_PI * pow(n_low, 2);
            }

            // Top cap only
            else if (contained
                     && h >= bot_of_lower_cap
                     && h >= top_of_lower_cap
                     && h >= bot_of_upper_cap
                     && h >= top_of_upper_cap) {
                if (log) { printf("\ttop cap only        \t"); }
                area = M_PI * pow(n_high, 2);
            }

            // Bottom cap + rectangular
            else if (contained
                     && h >= bot_of_lower_cap
                     && h <= top_of_lower_cap
                     && h <= bot_of_upper_cap
                     && h <= top_of_upper_cap) {
                if (log) { printf("\tbot + rectangular"); }
                if (h < mid_of_lower_cap) {
                    if (log) { printf("*"); }
                    area =    ellipseIntegral(a, a, rad)
                    - ellipseIntegral(a - e_low, a, rad)
                    + ellipseIntegral(n_low, n_low, n_low)
                    - ellipseIntegral(-n_low * cos(angles[0]), n_low, n_low);
                } else {
                    area =    ellipseIntegral(a, a, rad)
                    - ellipseIntegral(-d_low, a, rad)
                    + ellipseIntegral(n_low, n_low, n_low)
                    - ellipseIntegral(n_low * cos(angles[0]), n_low, n_low);
                }
                if (log) { printf("\t"); }
            }

            // Both caps
            else if (contained
                     && h >= bot_of_lower_cap
                     && h <= top_of_lower_cap
                     && h >= bot_of_upper_cap
                     && h <= top_of_upper_cap) {
                if (log) { printf("\tboth caps        \t"); }

                // Higher than both midpoints
                if (h > mid_of_upper_cap) {
                    // rectangular region, lower cap, upper cap
                    area =    ellipseIntegral(d_high + rad * AR / sin(angles[0]), a, rad)
                    - ellipseIntegral(d_high, a, rad)
                    + ellipseIntegral(n_low, n_low, n_low)
                    - ellipseIntegral(n_low * cos(angles[0]), n_low, n_low)
                    + ellipseIntegral(n_high, n_high, n_high)
                    - ellipseIntegral(-n_high * cos(angles[0]), n_high, n_high);
                }
                // Higher than the lower midpoint, lower than the higher midpoint
                else if (mid_of_lower_cap <= h && h <= mid_of_upper_cap) {
                    // rectangular region, lower cap, upper cap
                    area =    ellipseIntegral(d_high, a, rad)
                    - ellipseIntegral(-d_low, a, rad)
                    + ellipseIntegral(n_low, n_low, n_low)
                    - ellipseIntegral(n_low * cos(angles[0]), n_low, n_low)
                    + ellipseIntegral(n_high, n_high, n_high)
                    - ellipseIntegral(-n_high * cos(angles[0]), n_high, n_high);
                }
                // Lower than both midpoints
                else {
                    // rectangular region, lower cap, upper cap
                    area =    ellipseIntegral(d_low + rad * AR / sin(angles[0]), a, rad)
                    - ellipseIntegral(d_low, a, rad)
                    + ellipseIntegral(n_low, n_low, n_low)
                    - ellipseIntegral(-n_low * cos(angles[0]), n_low, n_low)
                    + ellipseIntegral(n_high, n_high, n_high)
                    - ellipseIntegral(n_high * cos(angles[0]), n_high, n_high);
                }
            }

            // Rectangular + top cap
            else if (contained
                     && h >= bot_of_lower_cap
                     && h >= top_of_lower_cap
                     && h >= bot_of_upper_cap
                     && h <= top_of_upper_cap) {
                if (log) { printf("\ttop + rectangular"); }
                if (h > mid_of_upper_cap) {
                    if (log) { printf("*"); }
                    area =  ellipseIntegral(a, a, rad)
                    - ellipseIntegral(a - e_high, a, rad)
                    + ellipseIntegral(n_high, n_high, n_high)
                    - ellipseIntegral(-n_high * cos(angles[0]), n_high, n_high);
                } else {
                    area =  ellipseIntegral(a, a, rad)
                    - ellipseIntegral(-d_high, a, rad)
                    + ellipseIntegral(n_high, n_high, n_high)
                    - ellipseIntegral(n_high * cos(angles[0]), n_high, n_high);
                }
                if (log) { printf("\t"); }
            }

            else {
                if (log) { printf("\tno intersection\n"); }
                continue;
            }

            if (isnan(area)) {
                throw std::runtime_error("NaN area encountered.");
                break;
            }

            if (log) { printf("\t%f\n", area); }

            // This is the maximum possible area
            if (area > (2 * rad) * (rad * AR) + M_PI * pow(rad, 2) + 100 * EPSILON) {
                throw std::runtime_error("Area too large. Check your calculations.");
            }
            
            density[j] += area;
            theta[j] += theta_ * area;
        }
    }
    
    for (int i = 0; i < nslices; i++) {
        // density[i] /= slice_area;
        theta[i] /= density[i];
        assert(0 <= theta[i] && theta[i] <= M_PI / 2);
    }
    
    return theta;
}






// Calculate the mass density for radially-parameterized rings
//std::vector<float> massDensityByRadius(float dr) {
//
//}

// Contact density by height

// Contact density by radius























#endif // PILE_HPP