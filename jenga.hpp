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


//MyObject createCeiling(void) {
//    // Initial height of 2.5 * bound
//    // Give it an intial velocity
//}


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

    // Turn off friction
    MU = 0;
    MU2 = 0;

    // Create pile of N rods
    placeRods(N);

    // Create ceiling and give it an initial velocity downward
//    MyObject ceiling = createCeiling();

    return;
}


// Dump data after a jenga run.
void jenga_dump(void) {

}













































#endif // JENGA_HPP