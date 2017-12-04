#ifndef PILE_HPP
#define PILE_HPP

// Structures to deal with metrics of a pile
//  All indexed midpoints
//  All indexed orientations
//  All double-indexed contact points
//  All indexed stability counts
//  xyz mass distribution

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

// Calculate the mass density for height-parameterized slices
// Based on area-density of a slice; no integration in the z-direction
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

    for (int i = 0; i < obj.size(); i++) {
        std::vector<float> pos = position(i);
        std::vector<float> angles = orientationAngles(i);

        float h_high = pos[2] + rad + rad * AR * cos(angles[0]);
        float h_low = pos[2] - rad - rad * AR * cos(angles[0]);

        float top_of_upper_cap = h_high - rad + rad * std::abs(sin(angles[0]));
        float bot_of_upper_cap = h_high - rad - rad * std::abs(sin(angles[0]));
        float top_of_lower_cap = h_low - rad + rad * std::abs(sin(angles[0]));
        float bot_of_lower_cap = h_low - rad - rad * std::abs(sin(angles[0]));

        assert(top_of_upper_cap > bot_of_upper_cap);
        assert(top_of_lower_cap > bot_of_lower_cap);

        int lowest_slice = floor(h_low / dh);
        int highest_slice = floor(h_high / dh);

        assert(lowest_slice >= 0 && highest_slice < nslices);

        // the plane we're intersecting against the rod
        float h = lowest_slice * dh;
        for (int j = lowest_slice; j <= highest_slice; j++) {
            float area = 0;

            if (h_low < h && h < h_high) {
                // Only intersects the rectangular portion of cylinder
                if (top_of_lower_cap <= h && h <= bot_of_upper_cap) {

                }
                // Intersects low cap
                else if (h <= top_of_lower_cap) {
                    // Intersects both caps
                    if (bot_of_upper_cap <= h) {

                    }
                    // Intersects rectangular portion + low cap
                    else if (h <= bot_of_upper_cap) {

                    }
                    // Intersects only low cap
                    else if (h <= bot_of_lower_cap) {

                    }
                }
                // Intersects high cap but not low cap
                else if (h >= bot_of_upper_cap) {
                    // Intersects rectangular portion + high cap
                    if (h <= top_of_upper_cap) {

                    }
                    // Intersects only high cap
                    else if (h >= top_of_upper_cap) {

                    }
                }
            }

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