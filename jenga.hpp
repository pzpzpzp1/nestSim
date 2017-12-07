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

void jenga_setup(void) {
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

    // Create ceiling
    
}













































#endif // JENGA_HPP