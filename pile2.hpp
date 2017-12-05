float a = rad / cos(angles[0]);

float offset_low = std::abs(h - mid_of_lower_cap);
float offset_high = std::abs(h - mid_of_upper_cap);

// h > mid_of_lower_cap --> x = offset_low
// else                 --> x = std::abs(h - bot_of_lower_cap)
float d_low = x * (cot(angles[0] + 1. / cot(angles[0])));
float n_low = rad * sqrt(1 - pow(offset_low / rad, 2));

// h < mid_of_upper_cap --> x = offset_high
// else                 --> x = std::abs(top_of_upper_cap - h)
float d_high = x * (cot(angles[0] + 1. / cot(angles[0])));
float n_high = rad * sqrt(1 - pow(offset_high / rad, 2)));


if (h_low < h && h < h_high) {

    // Only intersects the rectangular portion of cylinder
    if (top_of_lower_cap <= h && h <= bot_of_upper_cap) {
        printf("\trectangular only\t");

        area = M_PI * rad * a; // elliptical area
    }
    // Intersects low cap
    else if (h <= top_of_lower_cap) {
        // Intersects both caps
        if (bot_of_upper_cap <= h) {
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
            printf("\tlow + rectangular\t");

            assert(d_low < a);

            //                            area += ellipseIntegral(a, a, rad)
            //                                    - ellipseIntegral(-d_low, a, rad);

            if (h < mid_of_lower_cap) {
                area += ellipseIntegral(n_low, n_low, n_low)
                - ellipseIntegral(-n_low * cos(angles[0]),
                                  n_low, n_low);
            } else {
                area += ellipseIntegral(n_low, n_low, n_low)
                - ellipseIntegral(n_low * cos(angles[0]),
                                  n_low, n_low);
            }
        }
        // Intersects only low cap
        else if (h <= bot_of_lower_cap) {
            printf("\tlow cap only    \t");

            area += M_PI * pow(n_low, 2);
        }
    }
    // Intersects high cap but not low cap
    else if (h >= bot_of_upper_cap) {
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