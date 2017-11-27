/*************************************************************************
 *                                                                       *
 * Open Dynamics Engine, Copyright (C) 2001,2002 Russell L. Smith.       *
 * All rights reserved.  Email: russ@q12.org   Web: www.q12.org          *
 *                                                                       *
 * This library is free software; you can redistribute it and/or         *
 * modify it under the terms of EITHER:                                  *
 *   (1) The GNU Lesser General Public License as published by the Free  *
 *       Software Foundation; either version 2.1 of the License, or (at  *
 *       your option) any later version. The text of the GNU Lesser      *
 *       General Public License is included with this library in the     *
 *       file LICENSE.TXT.                                               *
 *   (2) The BSD-style license that is included with this library in     *
 *       the file LICENSE-BSD.TXT.                                       *
 *                                                                       *
 * This library is distributed in the hope that it will be useful,       *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the files    *
 * LICENSE.TXT and LICENSE-BSD.TXT for more details.                     *
 *                                                                       *
 *************************************************************************/

// ./nest -notex -noshadows

#include <ode/ode.h>
#include <drawstuff/drawstuff.h>
#include "texturepath.h"
#include <assert.h>
#include <vector>
#include <cmath>
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

//
////<---- Convex Object
//dReal planes[]= // planes for a cube, these should coincide with the face array
//{
//  1.0f ,0.0f ,0.0f ,0.25f,
//  0.0f ,1.0f ,0.0f ,0.25f,
//  0.0f ,0.0f ,1.0f ,0.25f,
//  -1.0f,0.0f ,0.0f ,0.25f,
//  0.0f ,-1.0f,0.0f ,0.25f,
//  0.0f ,0.0f ,-1.0f,0.25f
//  /*
//   1.0f ,0.0f ,0.0f ,2.0f,
//   0.0f ,1.0f ,0.0f ,1.0f,
//   0.0f ,0.0f ,1.0f ,1.0f,
//   0.0f ,0.0f ,-1.0f,1.0f,
//   0.0f ,-1.0f,0.0f ,1.0f,
//   -1.0f,0.0f ,0.0f ,0.0f
//   */
//};
//const unsigned int planecount=6;
//
//dReal points[]= // points for a cube
//{
//  0.25f,0.25f,0.25f,  //  point 0
//  -0.25f,0.25f,0.25f, //  point 1
//
//  0.25f,-0.25f,0.25f, //  point 2
//  -0.25f,-0.25f,0.25f,//  point 3
//
//  0.25f,0.25f,-0.25f, //  point 4
//  -0.25f,0.25f,-0.25f,//  point 5
//
//  0.25f,-0.25f,-0.25f,//  point 6
//  -0.25f,-0.25f,-0.25f,// point 7
//};
//const unsigned int pointcount=8;
//unsigned int polygons[] = //Polygons for a cube (6 squares)
//{
//  4,0,2,6,4, // positive X
//  4,1,0,4,5, // positive Y
//  4,0,1,3,2, // positive Z
//  4,3,1,5,7, // negative X
//  4,2,3,7,6, // negative Y
//  4,5,4,6,7, // negative Z
//};
////----> Convex Object
//
//// select correct drawing functions

#ifdef dDOUBLE
#define dsDrawBox dsDrawBoxD
#define dsDrawSphere dsDrawSphereD
#define dsDrawCylinder dsDrawCylinderD
#define dsDrawCapsule dsDrawCapsuleD
#define dsDrawConvex dsDrawConvexD
#endif


// some constants

#define NUM 1000			// max number of objects
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
static int write_world = 0;
static int show_body = 0;

// global variables start now
int maxNumContactsSimulated = 0;
FILE * fp;
float rad = .03; // rod radius
int AR = 50; // rod aspect ratio

struct MyFeedback {
    dJointFeedback fb;
    bool first;
};
static int doFeedback=0;
static MyFeedback feedbacks[MAX_FEEDBACKNUM];
static int fbnum=0;



// Return the spherical angles that correspond to the matrix R (of a rod that started aligned with the z axis)
std::vector<float> sphericalAnglesFromR(const dMatrix3 R, bool print) {
    std::vector<float> sphericalAngles;
    dVector3 initial, final;
    initial[2] = 1;
    float theta, phi;

    dMultiply0_331(final, R, initial);

    theta = acos(final[2]);
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
        throw std::runtime_error("No geom attached!");
    }
    return;
//    for (int k = 0; k < GPB; k++) {
//        if (obj[index].geom[k]) dGeomDestroy(obj[index].geom[k]);
//    }
    // memset(&obj[index], 0, sizeof(obj[index]));
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
        contact[i].surface.mode = dContactBounce | dContactSoftCFM;
        contact[i].surface.mu = dInfinity;
        contact[i].surface.mu2 = 0;
        contact[i].surface.bounce = 0.1;
        contact[i].surface.bounce_vel = 0.1;
        contact[i].surface.soft_cfm = 0.01;
    }
    if (int numc = dCollide (o1,o2,MAX_CONTACTS,&contact[0].geom,
                             sizeof(dContact))) {
        dMatrix3 RI;
        dRSetIdentity (RI);
        const dReal ss[3] = {0.02,0.02,0.02};
        for (i=0; i<numc; i++) {
            dJointID c = dJointCreateContact (world,contactgroup,contact+i);
            dJointAttach (c,b1,b2);
            if (show_contacts) dsDrawBox (contact[i].geom.pos,RI,ss);

            if (doFeedback && (b1==obj[selected].body || b2==obj[selected].body))
            {
                if (fbnum<MAX_FEEDBACKNUM)
                {
                    feedbacks[fbnum].first = b1==obj[selected].body;
                    dJointSetFeedback (c,&feedbacks[fbnum++].fb);
                }
                else fbnum++;
            }
        }
    }
}


// start simulation - set viewpoint

static void start()
{
  dAllocateODEDataForThread(dAllocateMaskAll);

  static float xyz[3] = {2.1640f,-1.3079f,1.7600f};
  static float hpr[3] = {125.5000f,-17.0000f,0.0000f};
  dsSetViewpoint (xyz,hpr);
  printf ("To drop another object, press:\n");
  printf ("   c for capsule.\n");
  printf ("   x for 50 capsules.\n");
  printf ("To select an object, press space.\n");
  printf ("To unselect objects, press u.\n");
  printf ("To toggle showing the geom AABBs, press a.\n");
  printf ("To toggle showing the contact points, press t.\n");
  printf ("To save the current state to 'state.dif', press 1.\n");
  printf ("To show joint feedbacks of selected object, press f.\n");
  printf ("To get the position/rotation of the selected object, press p.\n");

  printf ("To show/save total number of contacts, press w.\n");
  printf ("To show/save number of contacts for the selected object, press q.\n");
  printf ("To toggle immobility of current objects, press k.\n");
  printf ("To see if any objects are stable, press v.\n");
  printf ("To see if the selected object is stable, press b.\n\n");
}


char locase (char c)
{
    if (c >= 'A' && c <= 'Z') return c - ('a'-'A');
    else return c;
}

void drop(void) {
    size_t i;
    int j,k;
    dReal sides[3];
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
    for (k=0; k<3; k++) sides[k] = dRandReal()*0.5+0.1;

    dMatrix3 R;

    // drops capsules from positions within [-1..1, -1..1, 2..3]
    dBodySetPosition (new_obj.body,
                      dRandReal()*2-1,dRandReal()*2-1,dRandReal()+2);

    // all orientations uniformly
    float xy_angle = dRandReal() * 2 * M_PI;
    float h_angle = dRandReal() * 2 * M_PI;
    float r_angle = dRandReal() * 2 * M_PI;

    //dRFromAxisAndAngle(R, 1, 0, 0, 0.01);
    dRFromAxisAndAngle (R, sin(xy_angle)*sin(h_angle), cos(xy_angle)*
                        sin(h_angle), cos(h_angle),r_angle);

    dBodySetRotation (new_obj.body, R);

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

    return;
}

bool CheckStable(int rodInd) {
    fprintf(fp, "Checking stability of rod %d\n", rodInd);

    // get contacts
    dContact contact_array[num];
    float contactPos[num][3];
    int num_contacts = 0;
    dGeomID g0 = obj[rodInd].geom;
    for (int j = 0; j < num; j++) {
      dGeomID g1 = obj[j].geom;
      int numc = dCollide(g0, g1, MAX_CONTACTS, &contact_array[0].geom, sizeof(dContact));
      //assert(numc < 2); // technically 2 cylinders can only collide in one location. but numerical error sometimes reports two. in which case we'll only take the first collision.
      if (numc > 0) {
          contactPos[num_contacts][0] = contact_array[0].geom.pos[0];
          contactPos[num_contacts][1] = contact_array[0].geom.pos[1];
          contactPos[num_contacts][2] = contact_array[0].geom.pos[2];

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

    fprintf(fp, "\tOrdered angles: ");
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

// called when a key pressed
static void command (int cmd)
{
    int i, j, k;
    cmd = locase(cmd);

    // utility function
    if (cmd == 'z') {
        fprintf(fp, "num: %d\n", num);
        fprintf(fp, "Number of objects: %lu\n", obj.size());
        return;
    }

    // check all collisions
    else if (cmd == 'w')
    {
        dContact contact[MAX_CONTACTS];   // up to MAX_CONTACTS contacts per box-box
        for (i=0; i<MAX_CONTACTS; i++) {
            contact[i].surface.mode = dContactBounce | dContactSoftCFM;
            contact[i].surface.mu = dInfinity;
            contact[i].surface.mu2 = 0;
            contact[i].surface.bounce = 0.1;
            contact[i].surface.bounce_vel = 0.1;
            contact[i].surface.soft_cfm = 0.01;
        }

        int totalc = 0; // total number of contacts
        fprintf(fp, "num objs %d\n", num);

        // check all collisions between num x num objects
        for (i = 0; i < num; i ++){
            for (j = 0; j < num; j ++){
                dGeomID g1 = obj[i].geom;
                dGeomID g2 = obj[j].geom;

                totalc += dCollide (g1,g2,MAX_CONTACTS,&contact[0].geom,
                                    sizeof(dContact));
            }
        }
        // doublecounting collisions should always be even number
        assert(totalc%2==0);
        totalc/=2;

        if (totalc >= maxNumContactsSimulated )
        {
            maxNumContactsSimulated = totalc;
            fprintf(fp, "numContacts %d   maxcontacts %d \n\n", totalc, maxNumContactsSimulated);
        }
    }

    // Check collisions for selected object
    else if (cmd == 'q') {
        if (selected >= 0) {
            fprintf(fp, "\tselected object: %d\n", selected);

            // print object info (position, orientation, etc.)
            dVector3 pos;
            dMatrix3 R;
            dBodyCopyPosition(obj[selected].body, pos);
            dBodyCopyRotation(obj[selected].body, R);

            // get contacts
            dContact contact_array[num];
            int num_contacts = 0;
            dGeomID g0 = obj[selected].geom;
            for (j = 0; j < num; j++) {
                dGeomID g1 = obj[j].geom;
                if (dCollide(g0, g1, MAX_CONTACTS, &contact_array[0].geom, sizeof(dContact)) > 0) {
                    num_contacts++;
                }
            }

            fprintf(fp, "num contacts: %d\n\n", num_contacts);
        }
    }

    // Find if any rod is stable.
    else if(cmd == 'v')
    {
        num_stable = 0;
        for (int i = 0; i < num; i++) {
            bool isStable = CheckStable(i);
            if (isStable) {
                fprintf(fp, "\t**************** Stable rod found! %d ****************\n", i);
                num_stable++;
            }
        }

        fprintf(fp, "Number of stable/total objects: %d/%d\n", num_stable, num);
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
        for (j = 0; j < 50; j++) { drop(); }
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
            selected = -1;
            num--;
        }
    }
    else if (cmd == 'a') {
        show_aabb ^= 1;
    }
    else if (cmd == 't') {
        show_contacts ^= 1;
    }
    else if (cmd == '1') {
        write_world = 1;
    }
    else if (cmd == 'p' && selected >= 0)
    {
        const dReal* pos = dGeomGetPosition(obj[selected].geom);
        std::vector<float> angles = sphericalAnglesFromR(dGeomGetRotation(obj[selected].geom), false);
        printf("Object %d:\n", selected);
        printf("\tposition:\t[%f, %f, %f]\n", pos[0], pos[1], pos[2]);
        printf("\trotation:\t[%f, %f]\t\t(theta, phi)\n\n",
               angles[0], angles[1]);
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
        dReal radius,length;
        dGeomCapsuleGetParams (g,&radius,&length);
        dsDrawCapsule (pos,R,length,radius);
    }

    if (show_body) {
        dBodyID body = dGeomGetBody(g);
        if (body) {
            const dReal *bodypos = dBodyGetPosition (body);
            const dReal *bodyr = dBodyGetRotation (body);
            dReal bodySides[3] = { 0.1, 0.1, 0.1 };
            dsSetColorAlpha(0,1,0,1);
            dsDrawBox(bodypos,bodyr,bodySides);
        }
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
    dsSetColor (0,0,2);
    dSpaceCollide (space,0,&nearCallback);
    if (!pause) dWorldQuickStep (world, 0.02);

    if (write_world) {
        FILE *f = fopen ("state.dif","wt");
        if (f) {
            dWorldExportDIF (world,f,"X");
            fclose (f);
        }
        write_world = 0;
    }


  if (doFeedback)
  {
    if (fbnum>MAX_FEEDBACKNUM)
      printf("joint feedback buffer overflow!\n");
    else
    {
      dVector3 sum = {0, 0, 0};
      printf("\n");
      for (int i=0; i<fbnum; i++) {
        dReal* f = feedbacks[i].first?feedbacks[i].fb.f1:feedbacks[i].fb.f2;
        printf("%f %f %f\n", f[0], f[1], f[2]);
        sum[0] += f[0];
        sum[1] += f[1];
        sum[2] += f[2];
      }
      printf("Sum: %f %f %f\n", sum[0], sum[1], sum[2]);
      dMass m;
      dBodyGetMass(obj[selected].body, &m);
      printf("Object G=%f\n", GRAVITY*m.mass);
    }
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

//    for (int j=0; j < GPB; j++) {
//      if (i==selected) {
//        dsSetColor (0,0.7,1);
//      }
//      else {
//        dsSetColor (1,1,0);
//      }
//      drawGeom (obj[i].geom[j],0,0,show_aabb);
//    }
  }
}


int main (int argc, char **argv)
{
    // init global vars /***/
    maxNumContactsSimulated = 0;
    fp = stdout;//fopen("output.txt","w");
    AR = 50;

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
  dCreatePlane (space,0,0,1,0);
  // memset (obj,0,sizeof(obj));
    obj.clear();

  dThreadingImplementationID threading = dThreadingAllocateMultiThreadedImplementation();
  dThreadingThreadPoolID pool = dThreadingAllocateThreadPool(4, 0, dAllocateFlagBasicData, NULL);
  dThreadingThreadPoolServeMultiThreadedImplementation(pool, threading);
  // dWorldSetStepIslandsProcessingMaxThreadCount(world, 1);
  dWorldSetStepThreadingImplementation(world, dThreadingImplementationGetFunctions(threading), threading);

  // draw bounding box
//  dGeomID wall_N = dCreatePlane(space, 0, -1, 0, -2);
//  dGeomID wall_E = dCreatePlane(space, -1, 0, 0, -2);
//  dGeomID wall_S = dCreatePlane(space, 0, 1, 0, -2);
//  dGeomID wall_W = dCreatePlane(space, 1, 0, 0, -2);

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
