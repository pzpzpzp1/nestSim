#ifndef JENGA_HPP
#define JENGA_HPP

// Setup:
// - No gravity
// - No friction
// - Cylindrical boundary conditions
// - Cylindrical piston as "ceiling"
// - Initialize ## rods all at once, oriented completely randomly & placed
//   completely randomly within a radius of (1/2)*rad*AR + rad
// - Drop ceiling with a constant force


// Setup functions

void removeWall(int wall_index) {

}


MyObject createCeiling(void) {
    // Initial height of 3 * bound
    // Give it an intial velocity

    float ceil_radius = 0.98 * bound;
    float ceil_length = 3.;

    size_t i;
    int j,k;
    dMass m;
    int setBody = 0;

    MyObject new_obj;

    new_obj.body = dBodyCreate(world);
    dBodySetLinearVel(new_obj.body, 0., 0., -0.5);
    dBodySetData (new_obj.body,(void*) i);
    dMassSetCylinder (&m,DENSITY,3, ceil_radius, ceil_length);

    new_obj.geom = dCreateCylinder(space, ceil_radius, ceil_length);
    dGeomSetBody(new_obj.geom, new_obj.body);

    dBodySetMass (new_obj.body, &m);

    dMatrix3 R;

    dBodySetPosition (new_obj.body,
                      0., 0., 3. * bound + (0.5) * ceil_length);

//    float xy_angle, h_angle, r_angle;
//
//        xy_angle = dRandReal() * 2. * M_PI;
//        h_angle = dRandReal() * 2. * M_PI;
//        r_angle = dRandReal() * 2. * M_PI;
//    }
//
//    //dRFromAxisAndAngle(R, 1, 0, 0, 0.01);
//    dRFromAxisAndAngle (R, sin(xy_angle)*sin(h_angle), cos(xy_angle)*
//                        sin(h_angle), cos(h_angle),r_angle);
//    
//    dBodySetRotation (new_obj.body, R);

    return new_obj;

}


// Return true if the ceiling cylinder has stopped moving
bool ceilingStopped(MyObject ceiling) {
    return true;
}


float ceilingHeight(MyObject ceiling) {
    return 0.0;
}


// Place N rods, randomly oriented, randomly placed within (0.5)*rad*AR + rad
// distance of the walls & at a max height of 2 * bound
void placeRods(int N) {
    float height_step = 2. * bound / (N + 1.);

    for (int i = 0; i < N; i++) {
        drop(height_step * (i + 1));
    }

    return;
}


// N = # of rods to initialize
void jenga_setup(int N) {
    // Remove all objects
    while (num > 0) {
        removeObject(num - 1);
    }

    // Enforce cylindrical walls -- hopefully this is already done
    // ...

    // Turn off gravity
    dWorldSetGravity(world, 0., 0., 0.);

    // Increase damping
    dWorldSetLinearDamping(world, 0.001);

    // Turn off friction
    MU = 0;
    MU2 = 0;

    // Create pile of N rods
    placeRods(N);

    // Create ceiling and give it an initial velocity downward
    jenga_ceiling = createCeiling();

    return;
}


// Dump data after a jenga run.
void jenga_dump(void) {

}













































#endif // JENGA_HPP