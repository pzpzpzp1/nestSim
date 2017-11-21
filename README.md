# nestSim

## TODO

### Code
* Make functions that dump a nice output of contact orientations
* Make friction dynamically adjustable
* 

### Other
* Rewrite Mathematica constraint counting in terms of differences

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

* Corey O'Hern's code
* Maha's shoestring code
* [Project Chrono](http://api.projectchrono.org/tutorial_demo_bricks.html)
* [Open Dynamics Engine (ODE)](https://www.ode-wiki.org/wiki/)
    * [Demura](http://demura.net/english) (ODE demos/tutorials)

## Miscellaneous

* Basmati rice/spaghetti
* String case

## Other

* [Corey O'Hern](http://jamming.research.yale.edu/)
* Microtubule networks (formed w/ crosslinks)
* Hair simulations for collision detection
* [Nest documentary](https://www.youtube.com/watch?v=vxC85hSerkU)
* Beaver dams
