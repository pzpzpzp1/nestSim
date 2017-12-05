#ifndef PILE_HPP
#define PILE_HPP

// Structures to deal with metrics of a pile
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

float cot(float theta) {
    return cos(theta) / sin(theta);
}

// x = boundary value requested
// a = semi-major axis of ellipse = rad / cos(theta)
float ellipseIntegral(float x, float a, float b) {
    return b * (x * sqrt(1 - pow(x / a, 2)) + a * asin(x / a));
}

// Calculate the mass density for height-parameterized slices dh-thick.
// Based on area-density of a slice; no integration in the z-direction.
// Because of this, make sure that dh < rad.
std::vector<float> massDensityByHeight(float dh) {
    int nslices = ceil(highestMidpoint() / dh);
    std::vector<float> density(nslices, 0);

    // Area of a slice
    float area;
    if (nwalls == 5) {
        area = pow(bound * 2, 2);
    }
    else if (nwalls > 5) {
        area = M_PI * pow(bound, 2);
    }
    else {
        throw std::runtime_error("Can't take slices without boundaries.");
    }

    // Iterate through the objects
    for (int i = 0; i < obj.size(); i++) {
        printf("Rod %d:\n", i);

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

        // Sometimes the rods end up "intersecting with the floor"
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

        int lowest_slice = floor(h_low / dh);
        int highest_slice = floor(h_high / dh);

//        printf("low/high slice of %d: %d/%d\n", nslices, lowest_slice, highest_slice);

        // The end of the rod is below zero for some reason
        if (lowest_slice < 0) {
            printf("\tLowest slice: %d\n", lowest_slice);
            throw std::runtime_error("Slice too low.");
            continue;
        }

        // The end of the rod sticks above the nominal height of the pile
        // don't skip here; this is wrong
//        if (highest_slice >= nslices) {
//            printf("\tSlice too high.\n");
//            continue;
//        }

        // the plane we're intersecting against the rod
        float h = lowest_slice * dh + dh / 2.;
        for (int j = lowest_slice; j <= highest_slice; j++) {
            float area = 0;

            // The rod is essentially horizontal (w/in fraction of a degree)
            if (std::abs(M_PI/2 - angles[0]) < 0.01) {
                if (h_low < h && h < h_high) {
                    printf("\thorizontal     \t");
                    float offset = std::abs(h - pos[2]);
                    area = 2 * rad * sqrt(1. - pow(offset / rad, 2))
                             * (rad * AR)
                           + 2 * M_PI * rad * sqrt(1. - pow(offset / rad, 2));
                }
                else {
                    printf("\tno intersection 0\t");
                }
            }
            else {
                if (h_low < h && h < h_high) {
                    float a = rad / cos(angles[0]);

                    assert(a > 0.);

                    float offset = pos[2] - h; // must be signed

                    float offset_low = std::abs(h - mid_of_lower_cap);
                    float n_low = a * sqrt(1 - pow(offset_low / rad, 2));
                    float d_low = offset_low * (cot(angles[0]) + 1. / cot(angles[0]))
                                  + offset / cot(angles[0]);

                    float offset_high = std::abs(h - mid_of_upper_cap);
                    float n_high = a * sqrt(1 - pow(offset_high / rad, 2));
                    float d_high = offset_high * (cot(angles[0]) + 1. / cot(angles[0]))
                                   - offset / cot(angles[0]);

                    // Only intersects the rectangular portion of cylinder
                    if (top_of_lower_cap <= h && h <= bot_of_upper_cap) {
                        printf("\trectangular only\t");
                        // Area is an ellipse
                        area = M_PI * rad * a;
                    }
                    // Intersects low cap
                    else if (h <= top_of_lower_cap) {
                        assert(offset_low <= rad);

                        // Intersects both caps
                        if (bot_of_upper_cap <= h) {
                            assert(offset_high <= rad);
                            assert(d_low <= rad * sqrt(1 + AR / 4));
                            assert(d_high <= rad * sqrt(1 + AR / 4));

                            printf("\tboth caps       \t");
                            area += ellipseIntegral(d_high, a, rad)
                                    - ellipseIntegral(-d_low, a, rad);

                            area += ellipseIntegral(n_low, n_low, n_low)
                                   - ellipseIntegral(n_low * cos(angles[0]),
                                                     n_low, n_low);

                            area += ellipseIntegral(n_high, n_high, n_high)
                                    - ellipseIntegral(n_high * cos(angles[0]),
                                                      n_high, n_high);
                        }
                        // Intersects rectangular portion + low cap
                        else if (h <= bot_of_upper_cap) {
                            assert(d_low <= rad * sqrt(1 + AR / 4));

                            printf("\tlow + rectangular\t");
                            area += ellipseIntegral(a, a, rad)
                                    - ellipseIntegral(-d_low, a, rad);

                            area += ellipseIntegral(n_low, n_low, n_low)
                                   - ellipseIntegral(n_low * cos(angles[0]),
                                                     n_low, n_low);
                        }
                        // Intersects only low cap
                        else if (h <= bot_of_lower_cap) {
                            printf("\tlow cap only    \t");
                            area += M_PI * pow(n_low, 2);
                        }
                    }
                    // Intersects high cap but not low cap
                    else if (h >= bot_of_upper_cap) {
                        assert(offset_high <= rad);
                        // Intersects rectangular portion + high cap
                        if (h <= top_of_upper_cap) {
                            printf("\thigh + rectangular\t");
                            area += ellipseIntegral(-a, a, rad)
                                    - ellipseIntegral(d_high, a, rad);

                            area += ellipseIntegral(n_high, n_high, n_high)
                                    - ellipseIntegral(n_high * cos(angles[0]),
                                                      n_high, n_high);
                        }
                        // Intersects only high cap
                        else if (h >= top_of_upper_cap) {
                            printf("\thigh cap only     \t");
                            area += M_PI * pow(n_high, 2);
                        }
                    }
                }
            }

            if (isnan(area)) {
                printf("\tNaN area encountered.\n\n");
                break;
            }

            printf("\t%f\n", area);

            density[j] += area;
            h += dh;
        }
    }

    for (int i = 0; i < nslices; i++) {
        density[i] /= area;
    }

    return density;
}

// Calculate the mass density for radially-parameterized rings
//std::vector<float> massDensityByRadius(float dr) {
//
//}

// Contact density by height

// Contact density by radius























#endif // PILE_HPP