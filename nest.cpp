// ./nest -notex -noshadows

#include <ode/ode.h>
#include <drawstuff/drawstuff.h>
#include "texturepath.h"
#include <assert.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <map>
#include <iostream>
#include <algorithm>
//#include <Eigen/Dense>

#define EPSILON 0.000001

#ifdef _MSC_VER
#pragma warning(disable:4244 4305)  // for VC++, no precision loss complaints
#endif

#include "icosahedron_geom.h"

float fmodulo(float a, float b) {
    float res = a - b * floor(a / b);
    assert(0 <= res && res < b);
    return res;
}

bool orderbyfloat2(float *a, float * b)
{
    return a[0] < b[0];
}

bool orderbyfloat(float a, float b)
{
    return a < b;
}


// Sort vec1 by the contents of vec2.
template <class T>
std::vector<T> sortVectorBy(const std::vector<T> vec1,
                            const std::vector<T> vec2) {
    assert(vec1.size() == vec2.size());

    std::vector<T> result;
    std::vector< std::vector<T> > pairs;

    struct comparator {
        static bool compare(const std::vector<T> &pair1,
                            const std::vector<T> &pair2) {
            return pair1[1] < pair2[1];
        }
    };

    for (int i = 0; i < vec2.size(); i++) {
        std::vector<T> pair;
        pair.push_back(vec1[i]);
        pair.push_back(vec2[i]);
        pairs.push_back(pair);
    }

    std::sort(pairs.begin(), pairs.end(), comparator::compare);

    for (int i = 0; i < vec2.size(); i++) {
        result.push_back(pairs[i][0]);
    }

    return result;
}

// returns true if front has intersecting range with back. assumes ranges have been processed.
bool hasIntersectingFreeAngles(std::vector< std::pair<float, float> > front, std::vector< std::pair<float, float> > back) {
    for (int i = 0; i < front.size(); i++) {
        for (int j = 0; j < back.size(); j++) {
            std::pair<float, float> p1 = front[i];
            std::pair<float, float> p2 = back[j];

            float l = std::max(p1.first, p2.first);
            float r = std::min(p1.second, p2.second);

            if (l < r) { return true; }
        }
    }

    return false;
}

// takes ranges, returns ranges with the flipped versions
std::vector< std::pair<float, float> > Suspend(std::vector< std::pair<float, float> > ranges)
{
    std::vector< std::pair<float, float> > res;
    res.clear();

    for (int i = 0; i < ranges.size(); i++) {
        std::pair<float,float> thets = ranges.at(i);
        assert(thets.first <= thets.second);

        // for each range, duplicate the flipped range as well.
        res.push_back(thets);
        res.push_back(std::pair<float,float>
                      (thets.first + M_PI, thets.second + M_PI));
    }

    return res;
}

// make sure no values are more than 2pi while maintaining the range meaning.
std::vector< std::pair<float, float> > processRanges(std::vector< std::pair<float, float> > ranges){
    std::vector< std::pair<float, float> > res; res.clear();
    for (int i = 0; i < ranges.size(); i++) {
        std::pair<float,float> thets = ranges.at(i);

        // Range must be ordered, and free ranges can't be larger than 2pi
        assert(thets.first <= thets.second);
        assert(thets.second - thets.first <= 2 * M_PI + EPSILON);

        // Shift range so that the lower value is within 0..2pi
        while (thets.first > 2 * M_PI) {
            thets.first -= 2 * M_PI;
            thets.second -= 2 * M_PI;
        }

        // Check if the whole range lies within 0..2pi.
        if (thets.second <= 2 * M_PI) {
            res.push_back(thets);
        }
        // Otherwise, split the range into two pieces.
        else {
            assert(thets.second <= 4 * M_PI);
            // if the range is so large it encompasses all of s1, then just one range is enough.
            //            if (thets.second >= 4*M_PI){
            //                printf("*************RANGE_END > 4PI\n\n");
            //                res.push_back(std::pair<float,float>(0, 2*M_PI - .00001));
            //                return res;
            //            }
            // split range into two pieces.
            res.push_back(std::pair<float,float>
                          (thets.first, 2 * M_PI));
            res.push_back(std::pair<float,float>
                          (0, thets.second - 2 * M_PI));
        }
    }
    return res;
}

// given a set of angles, find if there are any angles in 0 to 2pi where the rod can move. Represented as
// the union of the ranges in the returned vector. ranges are represented as pairs where the latter is more than the former.
std::vector< std::pair<float, float> > getFreeAngles(std::vector<float> theta, int num_contacts, int p, bool front) // float theta[]
{
    //todo: implement
    std::vector< std::pair<float, float> > res;
    res.clear();

    // Set indices for range
    int startind = 0;
    int endind = num_contacts - 1;
    if (front) { endind = p; }
    else { startind = p + 1; }

    // If there are no contacts on this side of the pivot, all angles are free.
    if (startind > endind) {
        res.push_back(std::pair<float, float> (0, 2 * M_PI));
        return res;
    }

    // Create a vector of the angles on this side of the pivot.
    std::vector<float> fthetas;
    fthetas.clear();
    for (int i = startind; i <= endind; i++)
    {
        fthetas.push_back(theta[i]);
    }

    // this can't happen because of earlier index check
    assert(fthetas.size() != 0);

    // If there's a single contact on this side of the pivot, this end
    // can move out at any angle between theta + pi/2 and theta + 3pi/2.
    if (fthetas.size() == 1) {
        res.push_back(std::pair<float, float>
                      (fthetas[0] + M_PI / 2.0,
                       fthetas[0] + M_PI * 3.0 / 2.0));
        return res;
    }

    // 2+ contacts in the bunch. Sort the remaining angles.
    std::sort(fthetas.begin(), fthetas.end());

    assert(fthetas.size() >= 2);
    assert(fthetas[0] < fthetas[1]); // check the sorting order

    // look for gaps of more than pi in this bunch. implies free angle
    for (int i = 0; i < fthetas.size(); i++)
    {
        float tnext = fthetas[(i + 1) % fthetas.size()];
        float tcurr = fthetas[i];
        float dt = tnext - tcurr;

        // this should only happen once. at the last iteration of the loop.
        dt += (dt < 0 ? 2 * M_PI : 0);

        if (dt > M_PI)
        {
            float da = dt - M_PI; // size of free angle
            float la = tnext + M_PI/2.0;
            float ra = la + da;
            res.push_back(std::pair<float,float> (la,ra));
            return res; // dt can only be more than pi once.
        }
    }

    // res is empty. no free angles.
    return res;
}


//// select correct drawing functions

#ifdef dDOUBLE
#define dsDrawBox dsDrawBoxD
#define dsDrawSphere dsDrawSphereD
#define dsDrawCylinder dsDrawCylinderD
#define dsDrawCapsule dsDrawCapsuleD
#define dsDrawConvex dsDrawConvexD
#endif


// some constants

#define NUM 5000			// max number of objects
#define DENSITY (5.0)		// density of all objects
//#define GPB 3			// maximum number of geometries per body
#define MAX_CONTACTS 8          // maximum number of contact points per body
#define MAX_FEEDBACKNUM 20
#define GRAVITY         REAL(0.5)
#define USE_GEOM_OFFSET 1

// dynamics and collision objects
struct MyObject {
    dBodyID body;			// the body
    dGeomID geom;		// geometries representing this body
};

static int num=0;		// number of objects in simulation
static int num_stable = 0; // updated when 'v' is run
static int nextobj=0;		// next object to recycle if num==NUM
static dWorldID world;
static dSpaceID space;
//static MyObject obj[NUM];
// static
std::vector<MyObject> obj;
static dJointGroupID contactgroup;
static int selected = -1;	// selected object
static int show_aabb = 0;	// show geom AABBs?
static int show_contacts = 1;	// show contact points?
static int random_pos = 1;	// drop objects from random position?

static dReal MU = dInfinity;
static dReal MU2 = 0;

// global variables start now
int maxNumContactsSimulated = 0;
FILE * fp;
std::string savefilename = "savestate.txt";
std::string loadfilename = "savestate.txt";
float rad = .03; // rod radius

int simloopcount = 0;
int simloopmod = 70;
float AR = 50; // rod aspect ratio
int nwalls = 100; // number of walls. Should be 1 or 5 or 100
dGeomID *walls; // wall objects
float bound = rad*AR; //rad*AR/2.0+1.001; // shortest horz distance from walls to origin
float _MU;
float dropTheta;
bool hasboundary;

std::string contactdensityfilename;
std::string massdensityfilename;
std::string orientationdensityfilename;
std::string scalarsfilename;
std::string pythonplotter;
std::string parameters;
std::string filetypename;

struct MyFeedback {
    dJointFeedback fb;
    bool first;
};
static int doFeedback = 0;
static MyFeedback feedbacks[MAX_FEEDBACKNUM];
static int fbnum = 0;


// Return the spherical angles that correspond to the matrix R (of a rod that started aligned with the z axis)
std::vector<float> sphericalAnglesFromR(const dMatrix3 R, bool print) {
//    printf("SPHERICAL_ANGLES:\n");

    std::vector<float> sphericalAngles;
    dVector3 initial, final;
    initial[0] = 0;
    initial[1] = 0;
    initial[2] = 1;
    float theta, phi;

//    printf("\tR:\t%f, %f, %f\n\t\t%f,%f,%f\n\t\t%f,%f,%f\n",
//           R[0], R[1], R[2], R[4], R[5], R[6], R[8], R[9], R[10]);

    dMultiply0_331(final, R, initial);

//    printf("final vector: [%f, %f, %f]\n", final[0], final[1], final[2]);
//    printf("length: %f\n", sqrt(pow(final[0], 2) + pow(final[1], 2) + pow(final[2], 2)));

//    printf("\tfinal: [%f, %f, %f]\n", final[0], final[1], final[2]);

    // Enforce theta <= pi/2; we can do this because the rods are flippable
    theta = acos(std::abs(final[2]));
    theta = (theta <= M_PI / 2) ? theta : theta - M_PI / 2;
//    printf("\ttheta: %f\n\n", theta);
    sphericalAngles.push_back(theta);

    phi = atan(final[1] / final[0]);
    sphericalAngles.push_back(phi);

    if (print) {
        fprintf(fp, "theta: %f\n", theta);
        fprintf(fp, "phi:   %f\n", phi);
    }

    return sphericalAngles;
}

// Delete an object.
void removeObject(int index) {
    dBodyDestroy (obj[index].body);
    if (obj[index].geom) {
        dGeomDestroy(obj[index].geom);
        obj.erase(obj.begin() + index);
    }
    else {
        // This shouldn't happen
        throw std::runtime_error("No geom attached!");
    }
    return;
//    for (int k = 0; k < GPB; k++) {
//        if (obj[index].geom[k]) dGeomDestroy(obj[index].geom[k]);
//    }
    // memset(&obj[index], 0, sizeof(obj[index]));
}

std::vector<float> position(int index) {
    std::vector<float> pos;
    const dReal* pos_ = dGeomGetPosition(obj[index].geom);
    pos.push_back(pos_[0]);
    pos.push_back(pos_[1]);
    pos.push_back(pos_[2]);
    return pos;
}

std::vector<float> orientationAngles(int index) {
    std::vector<float> angles = sphericalAnglesFromR(dGeomGetRotation(obj[index].geom), false);

    return angles;
}

// this is called by dSpaceCollide when two objects in space are
// potentially colliding.
static void nearCallback (void *data, dGeomID o1, dGeomID o2)
{
    int i;
    // if (o1->body && o2->body) return;

    // this part confirms that position is in the center of the rod

    /*
     int interestedInd = selected > -1 ? selected : 0;
     dVector3 pos;
     dMatrix3 R;
     float halflen = AR*rad/2;
     dBodyCopyPosition(obj[interestedInd].body, pos);
     dBodyCopyRotation(obj[interestedInd].body, R);
     dMatrix3 RI2; dRSetIdentity (RI2); const dReal ss2[3] = {0.1,0.1,0.1};
     dsDrawBox (pos,RI2,ss2);

     dVector3 point; point[0]=0; point[1]=0; point[2] = halflen; //halflen ;
     dVector3 out;
     dMultiply0_331 (out,R,point);
     out[0]+=pos[0];
     out[1]+=pos[1];
     out[2]+=pos[2];
     dsDrawBox (out,RI2,ss2);
     */

    // exit without doing anything if the two bodies are connected by a joint
    dBodyID b1 = dGeomGetBody(o1);
    dBodyID b2 = dGeomGetBody(o2);
    if (b1 && b2 && dAreConnectedExcluding (b1,b2,dJointTypeContact)) return;

    dContact contact[MAX_CONTACTS];   // up to MAX_CONTACTS contacts per box-box
    for (i=0; i<MAX_CONTACTS; i++) {
        contact[i].surface.mode = dContactSoftCFM;
        contact[i].surface.mu = _MU;//dInfinity;
        contact[i].surface.rho = _MU;
        contact[i].surface.mu2 = 0;
        contact[i].surface.bounce = 0.0;
        contact[i].surface.bounce_vel = 10000000;
        contact[i].surface.soft_cfm = 0.01;
    }
    if (int numc = dCollide(o1, o2, MAX_CONTACTS, &contact[0].geom,
                             sizeof(dContact))) {
        dMatrix3 RI;
        dRSetIdentity (RI);
        const dReal ss[3] = {0.02,0.02,0.02};
        for (i=0; i<numc; i++) {
            dJointID c = dJointCreateContact (world,contactgroup,contact+i);
            dJointAttach (c,b1,b2);
            if (show_contacts) dsDrawBox (contact[i].geom.pos,RI,ss);

            if (doFeedback && (b1 == obj[selected].body || b2 == obj[selected].body)) {
                if (fbnum<MAX_FEEDBACKNUM) {
                    feedbacks[fbnum].first = (b1 == obj[selected].body);
                    dJointSetFeedback (c, &feedbacks[fbnum++].fb);
                }
                else fbnum++;
            }
        }
    }
}

void printInstructions(void) {
    printf("Objects:\n");
    printf ("  To drop another object, press:\n");
    printf ("     c for 1 capsule.\n");
    printf ("     x for 5 capsules.\n");
    printf ("     z for 50 capsules.\n");
    printf ("  To select an object, press space.\n");
    printf ("  To unselect objects, press u.\n");
    printf ("  To get the position/rotation of the selected object, press p.\n");
    printf ("  To toggle immobility of current objects, press k.\n\n");

    printf("Misc:\n");
    printf ("  To toggle showing the geom AABBs, press a.\n");
    printf ("  To toggle showing the contact points, press t.\n");
    printf ("  To show joint feedbacks of selected object, press f.\n\n");

    printf("Stability:\n");
    printf ("  To see which objects are stable, press v.\n");
    printf ("  To see if the selected object is stable, press b.\n\n");

    printf("Metrics:\n");
    printf("  To calculate the forces on all objects, press g.\n");
    printf("  To show the volumetric packing fraction, press n.\n");
    printf("  To show the highest midpoint value, press h.\n");
    printf("  To calculate mass density as a function of height, press m.\n");
    printf("  To calculate orientation density as a function of height, press o.\n\n");

    printf("Save/Load:\n");
    printf ("  To save the current state, press y.\n");
    printf ("  To load a state, press e.\n\n");

    printf ("To repeat the instructions, press i.\n\n");

    return;
}


// start simulation - set viewpoint

static void start()
{
  dAllocateODEDataForThread(dAllocateMaskAll);

  static float xyz[3] = {2.1640f,-1.3079f,1.7600f};
  static float hpr[3] = {125.5000f,-17.0000f,0.0000f};
  dsSetViewpoint (xyz,hpr);

    printInstructions();
}


char locase (char c)
{
    if (c >= 'A' && c <= 'Z') return c - ('a'-'A');
    else return c;
}

bool collidesWithWall(dGeomID g0){

    // get contacts with wall
    dContact contact_array[nwalls];
    for (int j = 0; j < nwalls; j++) {
      dGeomID g1 = walls[j];

      int numc = dCollide(g0, g1, MAX_CONTACTS, &contact_array[0].geom, sizeof(dContact));
      if (numc > 0) {
          return true;
      }
    }
    return false;
}

// copy-pasted code
#include "pile.hpp"


void drop(float hm = highest90percentileMidpoint(), float xyang = dRandReal()*2*M_PI) {
    size_t i;
    int j,k;
    dMass m;
    int setBody = 0;

    // Create and drop a rod.
    if (num < NUM) {
        i = num;
        num++;
    }
    else {
      i = nextobj;
      nextobj++;
      if (nextobj >= num) nextobj = 0;

      // destroy the body and geoms for slot i // this is useful for later. removing objects.
      dBodyDestroy (obj[i].body);
//      for (k=0; k < GPB; k++) {
//        if (obj[i].geom[k]) dGeomDestroy (obj[i].geom[k]);
//      }
      if (obj[i].geom) dGeomDestroy (obj[i].geom);
      memset (&obj[i], 0, sizeof(obj[i])); /****/
    }

    assert(i == obj.size());

    MyObject new_obj;

    new_obj.body = dBodyCreate (world);

    dBodySetLinearVel(new_obj.body, 0,0,-1);



    dBodySetData (new_obj.body,(void*) i);

    dMassSetCapsule (&m,DENSITY,3, rad, rad * AR);
    new_obj.geom = dCreateCapsule (space, rad, rad * AR);

    if (!setBody) {
        if (new_obj.geom) {
            dGeomSetBody(new_obj.geom, new_obj.body);
        }
    }
//      for (k=0; k < GPB; k++) {
//        if (new_obj.geom[k]) dGeomSetBody (new_obj.geom[k],new_obj.body);
//      }

    dBodySetMass (new_obj.body, &m);

    obj.push_back(new_obj);

    // keep changing pos/orientation until rod doesn't collide with wall.
    do{
        dMatrix3 R;

        // drops capsules from positions within [-1..1, -1..1, 2..3]
        dBodySetPosition (new_obj.body,
                          dRandReal()*2-1,dRandReal()*2-1, hm + rad*AR);

        // all orientations uniformly
        float xy_angle = xyang + dRandReal() * M_PI / 6 - M_PI / 12;
        float h_angle = M_PI/2;// dRandReal() * 2 * M_PI;
        float r_angle = dRandReal() * 2 * M_PI;

        r_angle = dropTheta + dRandReal() * M_PI/6-M_PI/12;

        //dRFromAxisAndAngle(R, 1, 0, 0, 0.01);
        dRFromAxisAndAngle (R, sin(xy_angle)*sin(h_angle), cos(xy_angle)*
                            sin(h_angle), cos(h_angle),r_angle);

        dBodySetRotation (new_obj.body, R);
    } while(collidesWithWall(new_obj.geom));

    return;
}

void dropBatch()
{
    float hm = highestMidpoint();
    float xyang = dRandReal()*2*M_PI;
    for (int j = 0; j < 5; j++) { drop(hm, xyang); }
}

bool CheckStable(int rodInd) {
    //fprintf(fp, "Checking stability of rod %d\n", rodInd);

    // get contacts
    dContact contact_array[num+nwalls];
    float contactPos[num+nwalls][3]; // add 5 for colliding with walls and floor
    int num_contacts = 0;
    dGeomID g0 = obj[rodInd].geom;
    for (int j = 0; j < num + nwalls; j++) {
      dGeomID g1;
      if(j < num){
        g1 = obj[j].geom;
      }
      else {
        g1 = walls[j-num];
      }

      int numc = dCollide(g0, g1, MAX_CONTACTS, &contact_array[0].geom, sizeof(dContact));
      //assert(numc < 2); // technically 2 cylinders can only collide in one location. but numerical error sometimes reports two. in which case we'll only take the first collision.
      if (numc > 0) {
          contactPos[num_contacts][0] = contact_array[0].geom.pos[0];
          contactPos[num_contacts][1] = contact_array[0].geom.pos[1];
          contactPos[num_contacts][2] = contact_array[0].geom.pos[2];

          //if(j >= num){
          //    fprintf(fp,"found contact with wall: %d\n",j-num);
          //}

          num_contacts++;
      }
    }
    // <= 3 contacts is always unstable.
    if (num_contacts <= 3) { return false; }

    dVector3 pos;
    dMatrix3 R, Rt;
    dBodyCopyPosition(obj[rodInd].body, pos);
    dBodyCopyRotation(obj[rodInd].body, R);

    // Rt is R transpose which is also the inverse of R
    // ODE matrices are represented 3x4 though, where the 4th col is useless which is why the transpose indices look wierd
    Rt[0] = R[0]; Rt[5] = R[5]; Rt[10] = R[10];
    Rt[1] = R[4]; Rt[2] = R[8]; Rt[6] = R[9];
    Rt[4] = R[1]; Rt[8] = R[2]; Rt[9] = R[6];

    // Organize all of the contacts
    std::vector<float> theta, z_order;
    for (int i = 0; i < num_contacts; i++) {
        // Recenter the contact
        for (int j = 0; j < 3; j++) {
            contactPos[i][j] -= pos[j];
        }

        // Unrotate
        dMultiply0_331(contactPos[i], Rt, contactPos[i]);

        float angle = fmodulo(atan2(contactPos[i][1],
                                    contactPos[i][0]), 2 * M_PI);
        theta.push_back(angle);
        z_order.push_back(contactPos[i][2]);
    }

    // sort the angles by their 'z_order' value
    theta = sortVectorBy(theta, z_order);

    //fprintf(fp, "\tOrdered angles: ");
    for (int i = 0; i < num_contacts; i++) {
        // re-center all angles so that the first angle is 0
        theta[i] -= theta[0];
        fprintf(fp, " %f", theta[i]);
    }
    fprintf(fp, "\n");

    // parametrize contact locations by l and theta, dist from start, and angle relative to y axis.
    /*   std::vector<float*> ltheta;
     ltheta.clear();
     float theta[num_contacts];
     for(int i = 0; i < num_contacts; i++) {
     // recenter
     contactPos[i][0] -= start[0];
     contactPos[i][1] -= start[1];
     contactPos[i][2] -= start[2];
     // reorient
     dVector3 rotCont;
     dMultiply0_331(rotCont, Rt, contactPos[i]);

     // calculate l and theta relative to y axis
     float * lthet = (float *)malloc(sizeof(float)*2);
     lthet[0] = rotCont[2];
     // float angle =
     // while(angle < 0) {angle += 2*M_PI;} angle = fmod(angle, 2*M_PI);
     lthet[1] = atan2(rotCont[1], rotCont[0]); //angle;
     ltheta.push_back(lthet);
     }
     std::sort(ltheta.begin(), ltheta.end(), orderbyfloat2);
     // fprintf(fp,"--- thetabegin\n");
     fprintf(fp, "\tthetas: ");
     for(int i = 0; i < num_contacts; i++){
     // recenter angles and assert necessary conditions
     theta[i] = (ltheta.at(i))[1] - (ltheta.at(0))[1];
     while(theta[i] < 0) {theta[i] += 8*M_PI;} theta[i] = fmod(theta[i], 2*M_PI);
     assert(theta[i] == theta[i]); // catch any nans. they pollute the batch.
     assert(theta[i] >= 0);
     fprintf(fp, "%f\t", theta[i]);
     }
     fprintf(fp, "\n");
     // fprintf(fp,"\n--- thetaend\n");
     */
    // calculate stability
    /*    for (int i = 0; i < num_contacts; i++)
     {
     free(ltheta.at(i)); // free unused memory.
     float dt = theta[(i+1)%num_contacts]-theta[i];
     if(std::abs(dt) > M_PI){ return false; } // a angle difference anywhere of more than pi means the rod can move in that dir.
     } */

    // Make a (deep) copy of the angles and sort the copy.
    std::vector<float> sorted_theta = theta;
    std::sort(sorted_theta.begin(), sorted_theta.end());

    // Append 2pi to the sorted angle vector
    sorted_theta.insert(sorted_theta.end(), 2 * M_PI);

    // Check adjacent angles for any pi-sized gaps.
    for (int i = 0; i < num_contacts; i++) {
        float dt = sorted_theta[(i + 1)] - sorted_theta[i];
        if (std::abs(dt) > M_PI) {
            fprintf(fp, "\tOpen hemispheres found.\n");
            return false;
        }
    }

    // Check all pivots.
    for (int p = 0; p < num_contacts; p++) {
        // p is a pivot. splits the contacts to [0:p] and [p+1:end]
        std::vector<std::pair<float, float> > front = processRanges(Suspend(getFreeAngles(theta, num_contacts, p, true)));
        std::vector<std::pair<float, float> > back = processRanges(Suspend(getFreeAngles(theta, num_contacts, p, false)));

        // get intersecting free angles
        if (hasIntersectingFreeAngles(front, back)) {
            return false;
        }
    }

    // no pivot revealed intersecting free angles. therefore stable
    return true;
}


// Save CSV file (delimited by spaces!) from a 2d vector
//      dataInCols == true  means that the first field values live in data[0][i]
//      dataInCols == false means that the first field values live in data[i][0]
template <class T>
void saveCSV(std::string filename,
             std::vector< std::string > fieldnames,
             std::vector< std::vector<T> > data,
             bool dataInCols=true) {

    if (dataInCols) { assert(fieldnames.size() == data.size()); }
    else            { assert(fieldnames.size() == data[0].size()); }

    std::ofstream outfile(filename.c_str());

    // Print field names
    for (int i = 0; i < fieldnames.size(); i++) {
        outfile << fieldnames[i] << " ";
    }
    outfile << std::endl;

    if (dataInCols) {
        // Print data row by row, its dumb but works
        for (int col = 0; col < data[0].size(); col++) {
            for (int row = 0; row < data.size(); row++) {
                outfile << data[row][col] << " ";
            }
            outfile << std::endl;
        }
    }
    else {
        for (int row = 0; row < data.size(); row++) {
            for (int col = 0; col < data[0].size(); col++) {
                outfile << data[row][col] << " ";
            }
            outfile << std::endl;
        }
    }

    outfile.close();
    return;
}


// Saves with the assumption of constant hardcoded mass, and walls
void SaveState(std::string filename){
    FILE * savefile = fopen(filename.c_str(),"w");
    if(savefile == NULL){
        fprintf(fp,"Save failed. Unable to open file. Maybe it doesn't exist? %s\n", filename.c_str());
        return;
    }

    fprintf(savefile, "numobj %d\n",num);
    // save radius length position orientation velocity rotationalV
    for(int i = 0; i < num; i ++)
    {
        dVector3 pos;
        dMatrix3 R, Rt;
        dBodyCopyPosition(obj[i].body, pos);
        dBodyCopyRotation(obj[i].body, R);
        // position rotation radius length
        fprintf(savefile, "ind %d pos %f %f %f ",i,pos[0],pos[1],pos[2]);
        fprintf(savefile, "R %f %f %f %f %f %f %f %f %f %f %f %f ",
                R[0],R[1],R[2],R[3],R[4],R[5],R[6],R[7],R[8],R[9],R[10],R[11]);
        fprintf(savefile, "Rad %f L %f ",rad,AR*rad);
        //vel rvel
        const dReal * lvel = dBodyGetLinearVel(obj[i].body);
        const dReal * rvel = dBodyGetAngularVel(obj[i].body);
        fprintf(savefile, "lvel %f %f %f rvel %f %f %f ",lvel[0],lvel[1],lvel[2],rvel[0],rvel[1],rvel[2]);
        fprintf(savefile, "\n");
    }
    fclose(savefile);
    fprintf(fp,"Save complete: %s\n", filename.c_str());
}

// Saves with the assumption of constant hardcoded mass, and walls. radius and length in loadifle are actually also ignored.
void LoadState(std::string filename){

    FILE * loadfile = fopen(filename.c_str(),"r");
    if(loadfile == NULL){
        fprintf(fp,"Load failed. Unable to open file. Maybe it doesn't exist? %s\n", filename.c_str());
        return;
    }

    fprintf(fp,"Load success. \n");
    while(num > 0)
    {
    //    fprintf(fp,"num %d \n",num);
        removeObject(0);
        num --;
    }
    fprintf(fp,"Cleared Current capsules. \n");

    dVector3 pos;
    dMatrix3 R, Rt;
    float lvl[3];
    float rvl[3];

    int numobjs, ind;
    float rad,l;
    fscanf(loadfile, "numobj %d\n", &numobjs);

    for(int i =0; i<numobjs;i ++)
    {
        fscanf(loadfile, "ind %d pos %f %f %f R %f %f %f %f %f %f %f %f %f %f %f %f Rad %f L %f lvel %f %f %f rvel %f %f %f \n",
               &ind,
               &(pos[0]),&(pos[1]),&(pos[2]),
               &(R[0]),&(R[1]),&(R[2]),&(R[3]),&(R[4]),&(R[5]),&(R[6]),&(R[7]),&(R[8]),&(R[9]),&(R[10]),&(R[11]),
               &rad,&l,
               &(lvl[0]),&(lvl[1]),&(lvl[2]),
               &(rvl[0]),&(rvl[1]),&(rvl[2])
               );

        //fprintf(fp,"dropping loaded rod. %d\n",ind);
        drop();
        dBodySetPosition (obj[num-1].body, pos[0],pos[1],pos[2]);
        dBodySetRotation (obj[num-1].body, R);
        dBodySetLinearVel(obj[num-1].body, lvl[0],lvl[1],lvl[2]);
        dBodySetAngularVel(obj[num-1].body, rvl[0],rvl[1],rvl[2]);
    }
    fclose(loadfile);
}


void save_mass_density(){
    if(num==0) {fprintf(fp,"nothing to save. no rods exist.\n"); return;}
    std::vector<float> mass_density = massDensityByHeight(rad / 2);
    printf("Mass density as a function of height: ");
    printFloatVector(mass_density); // print to terminal

    // ugh
    std::vector< std::vector<float> > data;
    data.push_back(mass_density);

    // only one field
    std::vector< std::string > fields;
    fields.push_back("mass_density");

    printf("Saving mass density plot to file...\n");
    saveCSV((massdensityfilename+parameters+filetypename).c_str(), fields, data, true);

    //printf("Plotting mass density. Simulation paused.\n");
    //system((pythonplotter+massdensityfilename+parameters+filetypename).c_str());
}

void save_contact_density(){
    if(num==0) {fprintf(fp,"nothing to save. no rods exist.\n"); return;}

    std::vector<float> cbh = ContactsByHeight(rad/2.0);

    printf("Number of contacts as a function of height: ");
    printFloatVector(cbh); // print to terminal

    // ugh
    std::vector< std::vector<float> > data;
    data.push_back(cbh);

    // only one field
    std::vector< std::string > fields;
    fields.push_back("contacts_by_height");

    printf("Saving contact density plot to file...\n");
    saveCSV((contactdensityfilename+parameters+filetypename).c_str(), fields, data, true);
    printf("Finished saving contact density plot to file...\n");

    //printf("Plotting contact density. Simulation paused.\n");
    //system((pythonplotter+contactdensityfilename+parameters+filetypename).c_str());
}
void save_orientation_density(){
    if(num==0) {fprintf(fp,"nothing to save. no rods exist.\n"); return;}
    std::vector<float> orientation_density = orientationDensityByHeight(rad / 2);
    printf("Orientation (polar angle) density as a function of height: ");
    printFloatVector(orientation_density); // print to terminal

    // ugh
    std::vector< std::vector<float> > data;
    data.push_back(orientation_density);

    // only one field
    std::vector< std::string > fields;
    fields.push_back("orientation_density");

    printf("Saving orientation density plot to file...\n");
    saveCSV((orientationdensityfilename+parameters+filetypename).c_str(), fields, data, true);

    //printf("Plotting orientation density. Simulation paused.\n");
    //system((pythonplotter+orientationdensityfilename+parameters+filetypename).c_str());

}

void save_scalars(){
    // save packing fraction and stable percentage.
    FILE * scalarfp = fopen((scalarsfilename+filetypename).c_str(),"w");
    num_stable = 0;
    for (int i = 0; i < num; i++) {
        bool isStable = CheckStable(i);
        if (isStable) {
            num_stable++;
        }
    }

    float pf = packingFraction(0);

    float stable_perc = ((float)num_stable)/((float)num);
    fprintf(scalarfp, "%f, %f\n", pf, stable_perc);
    fclose(scalarfp);
}

void save_all(){
    if(num==0) {fprintf(fp,"nothing to save. no rods exist.\n"); return;}
    save_mass_density();
    save_contact_density();
    save_orientation_density();
    save_scalars();
}


// called when a key pressed
static void command (int cmd)
{
    int i, j, k;
    cmd = locase(cmd);

    // utility function
    if (cmd == 'z') {
        for (i = 0; i < 50; i++) { drop(); }
        return;
    }

    else if (cmd == 'm') {
        save_mass_density();

        return;
    }
    else if (cmd == 'o') {

        save_contact_density();
        return;
    }

    else if (cmd == 'r') {
        save_orientation_density();
        return;
    }
      
    else if (cmd == 'i') {
        printInstructions();
        return;
    }

    else if (cmd == 'h') {
        printf("Highest midpoint value: %f\n", highestMidpoint());
        return;
    }

    else if (cmd == 'n') {
        printf("Volumetric packing fraction: %f\n", packingFraction(0));
        return;
    }

    // save current state
    else if (cmd == 'y') {
        SaveState(savefilename);
        return;
    }

    // load previous save state
    else if (cmd == 'e') {
        LoadState(loadfilename);
        return;
    }

    // Find if any rod is stable.
    else if(cmd == 'v')
    {
        num_stable = 0;
        for (int i = 0; i < num; i++) {
            bool isStable = CheckStable(i);
            if (isStable) {
                //fprintf(fp, "\t**************** Stable rod found! %d ****************\n", i);
                num_stable++;
            }
        }

        fprintf(fp, "Number of stable/total objects: %d/%d\n\n", num_stable, num);
    }

    // Check if the selected rod is stable.
    else if (cmd == 'b') {
        if (selected >= 0) { CheckStable(selected); }
    }

    // Drop a capsule (spherocylinder)
    else if (cmd == 'c') {
        drop();
    }

    // Drop lots of capsules
    else if (cmd == 'x') {
        dropBatch();
    }

    // toggle immobility
    else if (cmd == 'k') {
        if (dBodyIsKinematic(obj[0].body) == 0) {
            for (j = 0; j < num; j++) {
                dBodySetKinematic(obj[j].body);
                dBodySetLinearVel(obj[j].body, 0, 0, 0);
                dBodySetAngularVel(obj[j].body, 0, 0, 0);
            }
        }
        else {
            for (j = 0; j < num; j++) {
                dBodySetDynamic(obj[j].body);
            }
        }
    }

    else if (cmd == ' ') {
        selected++;
        if (selected >= num) selected = 0;
        if (selected < 0) selected = 0;
    }
    else if (cmd == 'u') {
        selected = -1;
    }
    else if (cmd == 'd') {
        if (selected >= 0) {
            removeObject(selected);
            printf("Object %d removed.\n", selected);
            selected--;
            num--;
        }
    }
    else if (cmd == 'a') {
        show_aabb ^= 1;
    }
    else if (cmd == 't') {
        show_contacts ^= 1;
    }
    else if (cmd == 'p' && selected >= 0)
    {
        std::vector<float> pos = position(selected);
        std::vector<float> angles = orientationAngles(selected);
        printf("Object %d:\n", selected);
        printf("\tposition:\t[%f, %f, %f]\n", pos[0], pos[1], pos[2]);
        printf("\trotation:\t[%f, %f]\t\t(theta, phi)\n\n",
               angles[0], angles[1]);
    }
    else if (cmd == 'f' && selected >= 0 && selected < num) {
        if (dBodyIsEnabled(obj[selected].body)) {
            doFeedback = 1;
        }
    }
}


// draw a geom
void drawGeom (dGeomID g, const dReal *pos, const dReal *R, int show_aabb)
{
    int i;

    if (!g) return;
    if (!pos) pos = dGeomGetPosition (g);
    if (!R) R = dGeomGetRotation (g);

    int type = dGeomGetClass (g);

    if (type == dCapsuleClass) {
        dReal radius, length;
        dGeomCapsuleGetParams (g,&radius,&length);
        dsDrawCapsule (pos,R,length,radius);
    }

    if (show_aabb) {
        // draw the bounding box for this geom
        dReal aabb[6];
        dGeomGetAABB (g,aabb);
        dVector3 bbpos;
        for (i=0; i<3; i++) bbpos[i] = 0.5*(aabb[i*2] + aabb[i*2+1]);
        dVector3 bbsides;
        for (i=0; i<3; i++) bbsides[i] = aabb[i*2+1] - aabb[i*2];
        dMatrix3 RI;
        dRSetIdentity (RI);
        dsSetColorAlpha (1,0,0,0.5);
        dsDrawBox (bbpos,RI,bbsides);
    }
}


// simulation loop

static void simLoop (int pause)
{
    simloopcount++;
    if(simloopcount % simloopmod == 0 && num < 300){
        dropBatch();
    }
    if(simloopcount > 4500){
        save_all();
        assert(0); // not the most elegant way to leave a program i know.
    }

    dsSetColor (0,0,2);
    dSpaceCollide (space,0,&nearCallback);
    if (!pause) { dWorldQuickStep (world, 0.02); }

    if (doFeedback)
    {
        if (fbnum>MAX_FEEDBACKNUM) {
            printf("joint feedback buffer overflow!\n");
        }
        else
        {
            // Sum the forces on the object
            dVector3 sum = {0, 0, 0};

            for (int i = 0; i < fbnum; i++) {
                dReal* f = feedbacks[i].first ? feedbacks[i].fb.f1 : feedbacks[i].fb.f2;
                printf("%f %f %f\n", f[0], f[1], f[2]);
                sum[0] += f[0];
                sum[1] += f[1];
                sum[2] += f[2];
            }
            printf("Sum: %f %f %f\n", sum[0], sum[1], sum[2]);

            // Gravitational force on the object
            dMass m;
            dBodyGetMass(obj[selected].body, &m);
            printf("Gravitational force: %f\n", GRAVITY*m.mass);
        }

        // Reset
        doFeedback = 0;
        fbnum = 0;
    }

  // remove all contact joints
  dJointGroupEmpty (contactgroup);

  dsSetColor (1,1,0);
  dsSetTexture (DS_WOOD);

  for (int i=0; i<num; i++) {

      if (i==selected) {
          dsSetColor (0,0.7,1);
      }
      else {
          dsSetColor (1,1,0);
      }
      drawGeom (obj[i].geom,0,0,show_aabb);
  }

  return;
}


int main(int argc, char **argv)
{

    fp = stdout;
    if(argc == 1)
    {
        _MU = 1000000;
        dropTheta = M_PI/2; // 1.57
        AR = 50;
        hasboundary = true;
    }
    else if(argc == 5){
        _MU = atof(argv[1]);
        dropTheta = atof(argv[2]);
        AR = atof(argv[3]);
        hasboundary = atoi(argv[4]) != 0;
    }
    else
    {
        fprintf(fp, "need 4 arguments! friction, drop angle, AR, boundary\n");
        return 0;
    }
    fprintf(fp, "nest sim started with parameters %f %f %f %d\n",_MU, dropTheta, AR, hasboundary);

    nwalls = 100;
    walls = (dGeomID *) malloc(sizeof(dGeomID)*nwalls);//[nwalls]
    bound = 1;//AR * rad / 2 + 1;
    if(!hasboundary){bound = 10000;} // there's some bug with not initializing the walls. so we can just secretly make them really far.
    fprintf(fp, "walls initiated\n");

    contactdensityfilename = "outputdata/contact_density";
    massdensityfilename = "outputdata/mass_density";
    orientationdensityfilename = "outputdata/orientation_density";
    scalarsfilename = "outputdata/scalar_data";
    pythonplotter = "python3 plotter/plotter.py ";
    filetypename = ".csv";
    parameters = "_mu_" + std::to_string(_MU) + "_dtheta_" + std::to_string(dropTheta) + "_AR_" + std::to_string(AR) + "_bounded_" + std::to_string(hasboundary);

    std::cout << parameters;

    // To check that saveCSV(...) works
//    std::vector<int> xvals;
//    xvals.push_back(0);
//    xvals.push_back(1);
//    xvals.push_back(2);
//    xvals.push_back(3);
//    std::vector<int> yvals;
//    yvals.push_back(2);
//    yvals.push_back(3);
//    yvals.push_back(1);
//    yvals.push_back(2);
//
//    std::vector< std::vector<int> > data;
//    data.push_back(xvals);
//    data.push_back(yvals);
//
//    std::vector< std::string > fields;
//    fields.push_back("x");
//    fields.push_back("y");
//
//    std::string filename = "plotter/example.csv";
//
//    saveCSV(filename, fields, data, true);


    // init global vars /***/
    maxNumContactsSimulated = 0;


    // setup pointers to drawstuff callback functions
    dsFunctions fn;
    fn.version = DS_VERSION;
    fn.start = &start;
    fn.step = &simLoop;
    fn.command = &command;
    fn.stop = 0;
    fn.path_to_textures = DRAWSTUFF_TEXTURE_PATH;
    
    // create world
    dInitODE2(0);
    world = dWorldCreate();
    space = dHashSpaceCreate (0);
    contactgroup = dJointGroupCreate (0);
    dWorldSetGravity (world,0,0,-GRAVITY);
    dWorldSetCFM (world,1e-5);
    dWorldSetAutoDisableFlag (world,1);
    
#if 1
    
    dWorldSetAutoDisableAverageSamplesCount( world, 10 );
    
#endif
  dWorldSetLinearDamping(world, 0.00001);
  dWorldSetAngularDamping(world, 0.005);
  dWorldSetMaxAngularSpeed(world, 200);

  dWorldSetContactMaxCorrectingVel (world,0.1);
  dWorldSetContactSurfaceLayer (world,0.001);

//  assert(nwalls == 1 || nwalls == 5); // need floor or floor+sidewalls
  dGeomID wall_D = dCreatePlane (space,0,0,1,EPSILON);
  walls[0] = wall_D;
    obj.clear();

  dThreadingImplementationID threading = dThreadingAllocateMultiThreadedImplementation();
  dThreadingThreadPoolID pool = dThreadingAllocateThreadPool(4, 0, dAllocateFlagBasicData, NULL);
  dThreadingThreadPoolServeMultiThreadedImplementation(pool, threading);
  // dWorldSetStepIslandsProcessingMaxThreadCount(world, 1);
  dWorldSetStepThreadingImplementation(world, dThreadingImplementationGetFunctions(threading), threading);

  // draw bounding box
//  float bound = rad*AR/2.0+1.001; // moved to globals
    if (nwalls == 5){
        dGeomID wall_N = dCreatePlane(space, 0, -1, 0, -bound);
        dGeomID wall_E = dCreatePlane(space, -1, 0, 0, -bound);
        dGeomID wall_S = dCreatePlane(space, 0, 1, 0, -bound);
        dGeomID wall_W = dCreatePlane(space, 1, 0, 0, -bound);
        walls[1]=wall_N;
        walls[2]=wall_E;
        walls[3]=wall_S;
        walls[4]=wall_W;
    }
    else if (nwalls > 5) {
        fprintf(fp, "making many walls\n");
        // put in like a ton of walls
        float angle_step = 2 * M_PI / (nwalls - 1.);
        float angle = 0;

        for (int i = 0; i < nwalls; i++) {
            walls[i] = dCreatePlane(space, cos(angle), sin(angle), 0, -bound);
            angle += angle_step;
        }
    }

  // run simulation
  dsSimulationLoop (argc,argv,800,600,&fn);

  dThreadingImplementationShutdownProcessing(threading);
  dThreadingFreeThreadPool(pool);
  dWorldSetStepThreadingImplementation(world, NULL, NULL);
  dThreadingFreeImplementation(threading);

  dJointGroupDestroy (contactgroup);
  dSpaceDestroy (space);
  dWorldDestroy (world);
  dCloseODE();
  return 0;
}
