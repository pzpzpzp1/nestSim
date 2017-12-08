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


// Return true if the ceiling cylinder has stopped moving
bool ceilingStopped(MyObject ceiling) {

}


float ceilingHeight(MyObject ceiling) {

}


// N = # of rods to initialize
void jenga_setup(int N) {
    // Remove all objects
    if (num > 0) {
        ;
    }

    // Enforce cylindrical walls
    if (nwalls <= 5) {
        ;
    }

    // Turn off gravity
    dWorldSetGravity(world, 0, 0, 0);

    // Turn off friction: stupid but works
    MU = 0;
    MU2 = 0;

    // Create pile of N rods
    for (int i = 0; i < N; i++) {
        ;
    }

    // Create ceiling and give it an initial velocity downward
    // ...

    return;
}


// Dump data after a jenga run.
void jenga_dump(void) {

}













































#endif // JENGA_HPP