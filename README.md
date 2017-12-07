# nestSim

## Notes

* Aspect ratio is really 2x true aspect ratio because it multiplies the  
  radius, not diameter
* Mass/orientation density functions are *not* normalized

## To-do

* "Drop protocols"
* Height/radius/orientation/density metrics
* Other papers review
* Calculate nematic order
* Assumptions -- uniform size, BC's, rigid rods
    * ODE equations of motion
    * dropping randomly: kill rods that intersect w walls, or just drop them within 
      1/2 * rad * AR of the wall
    * totally classical, nothing thermal. what stochasticity and why.
* Movies

## Contact Forces
* These are not stable ([1](https://groups.google.com/d/msg/ode-users/kPfQIo-QOlE/I3EwRFI6BwAJ), [2](https://groups.google.com/d/msg/ode-users/OK1V4SXrb_k/ykJB7n7j4HIJ)) and can fluctuate wildly in value (e.g., select an object and hit 'f' multiple times to output the joint feedbacks)
* Must be low-pass filtered

## Pile Metrics

### Volumetric packing fraction
* Maximum packing fraction would be the same as that of circles ~ 0.9
* N=200 piles --> packing fraction of 0.05
* For capsules, the "length" does not include the caps

### Height/radial spread

(Radial spread will require cylindrical walls)

* Highest midpoint position (discard if rod is leaning against a wall)
* Midpoint distribution (as a function of height/radius)
* Mass distribution (as a function of height/radius)
* Contact density (as a function of height/radius)
* Orientation density (as a function of height/radius)

## Stuff to try

### Plots

* Height/density as a function of AR
* Number of rattlers after compressing w/o gravity. Start with isotropic 
  distribution in height.
    * What % can be removed? Disks ~ 1-5%
* No friction, gravity, or rattlers: (# contacts) = 5 * (# DOF)
    * As a function of aspect ratio

### Effect of a single object removal

1. Start with a pile and get metrics
2. Save state
3. Remove a single rod, noting its information
4. Measure the pile metrics again
5. Reload state and repeat

### Drop angles

1. Drop $N$ single rods from a particular angular distribution
2. Get pile metrics

### Drop positions

1. Drop $N$ single rods from a particular spatial distribution
2. Get pile metrics

### Drop heights

1. Drop $N$ single rods from a particular height
2. Get pile metrics

## Parameters

* Friction
* Longitudinal translation
* Branching
* Gravity
* Boundary conditions
* Aspect ratio

## Questions

* How to establish local stability?
* What's the relationship between local and global stability?
* What are properties of "minimally stable" structures?
* How is the stability distributed throughout the nest? E.g., how many of the rods are truly crucial to the nest's structure?
* How do real birds' nests compare to a rigid rod model?
* What are some properties of "regular" structures made by rods? How do these differ from structures of randomly placed rods?
* What is a good, procedural way to build a nest?
* What's a good model/explanation for the "drop into a pile, kick out the middle" building behavior?

## Local stability

Note: all of the counting arguments ignore longitudinal translation and metastable edge cases.

### Reducibility

Question: In a given $n$-contact configuration that is stable, is there always a subset of 4 contacts that is stable? How "redundant" is each contact?

### Counting arguments

For a given cylindrical rod with contacts distributed (uniformly) randomly along the length, it's easy to determine if the rod will be fixed in place by these contacts. One way is to look at "pivot points" between any pairs of contacts along the axis. If there exists an "open hemisphere" on one side of the pivot, and the opposite "open hemisphere" on the other side, then the rod is not stable. After checking all possible pivots, the rod is guaranteed to be stable.

4+ contacts are required for stability, and these contacts must be arranged in a particular way depending on their ordering along the rod's axis.

## Code

### Stability checking

* `bool CheckStable(int rodInd)`
  - Identifies the number and positions of contacts on a given rod
    - Returns false if the number of contacts is <= 3
  - Contact positions are shifted so that the rod is centered at (0, 0, 0), vertically oriented
  - Contact angles are assigned based on their projection onto the x-y plane
  - Contact orderings are assigned based on their z-value, then shifted so that the first contact has an angle of 0
    - Returns false if a pair of adjacent angles has a gap larger than pi
  - Check each pivot:
    - `getFreeAngles`
      - Args: z-sorted angles, # contacts, pivot index, front flag
      - Returns: a vector of pairs corresponding to ranges within 0..2pi where that end of the rod is free to move
    - `Suspend`
      - Takes a vector of range pairs and appends the flipped ranges as well
    - `processRanges`
      - Does some bounds/order checking on range pairs
      - Returns an equivalent vector of range pairs, where all pairs lie within 0..2pi
    - `hasIntersectingFreeAngles`
      - Takes two processed vectors of range pairs and check if there is any intersection

### Notes

* `fmod` does not work with negative numbers; use `fmodulo` instead.

## Miscellaneous

* Basmati rice/spaghetti
* String case

## Other

* Shoestring code
* [Project Chrono](http://api.projectchrono.org/tutorial_demo_bricks.html)
* [Open Dynamics Engine (ODE)](https://www.ode-wiki.org/wiki/)
* [Demura](http://demura.net/english) (ODE demos/tutorials)

* [Corey O'Hern](http://jamming.research.yale.edu/)
* Microtubule networks (formed w/ crosslinks)
* Hair simulations for collision detection
* [Nest documentary](https://www.youtube.com/watch?v=vxC85hSerkU)
* Beaver dams
