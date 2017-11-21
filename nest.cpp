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

#ifdef _MSC_VER
#pragma warning(disable:4244 4305)  // for VC++, no precision loss complaints
#endif

#include "icosahedron_geom.h"


//<---- Convex Object
dReal planes[]= // planes for a cube, these should coincide with the face array
{
  1.0f ,0.0f ,0.0f ,0.25f,
  0.0f ,1.0f ,0.0f ,0.25f,
  0.0f ,0.0f ,1.0f ,0.25f,
  -1.0f,0.0f ,0.0f ,0.25f,
  0.0f ,-1.0f,0.0f ,0.25f,
  0.0f ,0.0f ,-1.0f,0.25f
  /*
   1.0f ,0.0f ,0.0f ,2.0f,
   0.0f ,1.0f ,0.0f ,1.0f,
   0.0f ,0.0f ,1.0f ,1.0f,
   0.0f ,0.0f ,-1.0f,1.0f,
   0.0f ,-1.0f,0.0f ,1.0f,
   -1.0f,0.0f ,0.0f ,0.0f
   */
};
const unsigned int planecount=6;

dReal points[]= // points for a cube
{
  0.25f,0.25f,0.25f,  //  point 0
  -0.25f,0.25f,0.25f, //  point 1

  0.25f,-0.25f,0.25f, //  point 2
  -0.25f,-0.25f,0.25f,//  point 3

  0.25f,0.25f,-0.25f, //  point 4
  -0.25f,0.25f,-0.25f,//  point 5

  0.25f,-0.25f,-0.25f,//  point 6
  -0.25f,-0.25f,-0.25f,// point 7
};
const unsigned int pointcount=8;
unsigned int polygons[] = //Polygons for a cube (6 squares)
{
  4,0,2,6,4, // positive X
  4,1,0,4,5, // positive Y
  4,0,1,3,2, // positive Z
  4,3,1,5,7, // negative X
  4,2,3,7,6, // negative Y
  4,5,4,6,7, // negative Z
};
//----> Convex Object

// select correct drawing functions

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
#define GPB 3			// maximum number of geometries per body
#define MAX_CONTACTS 20          // maximum number of contact points per body
#define MAX_FEEDBACKNUM 20
#define GRAVITY         REAL(0.5)
#define USE_GEOM_OFFSET 1

// dynamics and collision objects

struct MyObject {
  dBodyID body;			// the body
  dGeomID geom[GPB];		// geometries representing this body
};

static int num=0;		// number of objects in simulation
static int nextobj=0;		// next object to recycle if num==NUM
static dWorldID world;
static dSpaceID space;
static MyObject obj[NUM];
static dJointGroupID contactgroup;
static int selected = -1;	// selected object
static int show_aabb = 0;	// show geom AABBs?
static int show_contacts = 1;	// show contact points?
static int random_pos = 1;	// drop objects from random position?
static int write_world = 0;
static int show_body = 0;

// global variables start now /***/
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

// this is called by dSpaceCollide when two objects in space are
// potentially colliding.

static void nearCallback (void *data, dGeomID o1, dGeomID o2)
{
  int i;
  // if (o1->body && o2->body) return;

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
  printf ("To select an object, press space.\n");
  printf ("To unselect objects, press u.\n");
  printf ("To dump transformation data for the selected object, press p.\n");
  printf ("To toggle showing the geom AABBs, press a.\n");
  printf ("To toggle showing the contact points, press t.\n");
  printf ("To toggle dropping from random position/orientation, press r.\n");
  printf ("To save the current state to 'state.dif', press 1.\n");
  printf ("To show joint feedbacks of selected object, press f.\n\n");

  printf ("To show/save total number of contacts, press w.\n");
  printf ("To show/save number of contacts for the selected object, press q.\n");
  printf ("To toggle immobility of current objects, press k.\n\n");
}


char locase (char c)
{
  if (c >= 'A' && c <= 'Z') return c - ('a'-'A');
  else return c;
}


// called when a key pressed

static void command (int cmd)
{
  size_t i;
  int j,k;
  dReal sides[3];
  dMass m;
  int setBody;

  cmd = locase (cmd);

  // /***/ check all collisions
  if (cmd == 'w')
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
        dGeomID g1 = obj[i].geom[0];
        dGeomID g2 = obj[j].geom[0]; // i->j

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
      fprintf(fp, "%d\n\n", maxNumContactsSimulated);
    }
  }

  // Check collisions for selected object
  if (cmd == 'q') {
    if (selected >=0) {
      fprintf(fp, "selected object: %d\n", selected);

      // print object info (position, orientation, etc.)
      // ...

      // get contacts
      dContact contact_array[num];
      int num_contacts = 0;
      dGeomID g0 = obj[selected].geom[0];
      for (j = 0; j < num; j++) {
        dGeomID g1 = obj[j].geom[0];
        if (dCollide(g0, g1, MAX_CONTACTS, &contact_array[0].geom, sizeof(dContact)) > 0) {
          num_contacts++;
        }
      }

      fprintf(fp, "num contacts: %d\n\n", num_contacts);
    }
  }

  // Drop a capsule (spherocylinder)
  if (cmd == 'c')
  {
    setBody = 0;
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
      for (k=0; k < GPB; k++) {
        if (obj[i].geom[k]) dGeomDestroy (obj[i].geom[k]);
      }
      memset (&obj[i],0,sizeof(obj[i]));
    }

    obj[i].body = dBodyCreate (world);
    for (k=0; k<3; k++) sides[k] = dRandReal()*0.5+0.1;

    dMatrix3 R;
    if (0==0) // if (random_pos) /**/
    {

      // [-1,1 -1,1 2]
      dBodySetPosition (obj[i].body,
                        dRandReal()*2-1,dRandReal()*2-1,dRandReal()+2);

      // all orientations uniformly /***/
      float xy_angle = dRandReal() * 2 * M_PI;
      float h_angle = dRandReal() * 2 * M_PI;
      float r_angle = dRandReal() * 2 * M_PI;

      dRFromAxisAndAngle (R, sin(xy_angle)*sin(h_angle), cos(xy_angle)*sin(h_angle),
                          cos(h_angle),r_angle);
    }

    dBodySetRotation (obj[i].body,R);
    dBodySetData (obj[i].body,(void*) i);

    if (cmd == 'c') { /***/
      dMassSetCapsule (&m,DENSITY,3, rad, rad * AR);
      obj[i].geom[0] = dCreateCapsule (space, rad, rad * AR);
    }

    if (!setBody)
      for (k=0; k < GPB; k++) {
        if (obj[i].geom[k]) dGeomSetBody (obj[i].geom[k],obj[i].body);
      }

    dBodySetMass (obj[i].body,&m);

    // Make immovable (infinite mass)
    // dBodySetKinematic(obj[i].body);
  }

  if (cmd == 'k') {
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

  if (cmd == ' ') {
    selected++;
    if (selected >= num) selected = 0;
    if (selected < 0) selected = 0;
  }
  else if (cmd == 'u') { /***/ // unselect body
    selected = -1;
  }
  else if (cmd == 'a') {
    show_aabb ^= 1;
  }
  else if (cmd == 't') {
    show_contacts ^= 1;
  }
  else if (cmd == 'r') {
    random_pos ^= 1;
  }
  else if (cmd == '1') {
    write_world = 1;
  }
  else if (cmd == 'p'&& selected >= 0)
  {
    const dReal* pos = dGeomGetPosition(obj[selected].geom[0]);
    const dReal* rot = dGeomGetRotation(obj[selected].geom[0]);
    printf("POSITION:\n\t[%f,%f,%f]\n\n",pos[0],pos[1],pos[2]);
    printf("ROTATION:\n\t[%f,%f,%f,%f]\n\t[%f,%f,%f,%f]\n\t[%f,%f,%f,%f]\n\n",
           rot[0],rot[1],rot[2],rot[3],
           rot[4],rot[5],rot[6],rot[7],
           rot[8],rot[9],rot[10],rot[11]);
  }
  else if (cmd == 'f' && selected >= 0 && selected < num) {
    if (dBodyIsEnabled(obj[selected].body))
      doFeedback = 1;
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
  if (type == dBoxClass) {
    dVector3 sides;
    dGeomBoxGetLengths (g,sides);
    dsDrawBox (pos,R,sides);
  }
  else if (type == dSphereClass) {
    dsDrawSphere (pos,R,dGeomSphereGetRadius (g));
  }
  else if (type == dCapsuleClass) {
    dReal radius,length;
    dGeomCapsuleGetParams (g,&radius,&length);
    dsDrawCapsule (pos,R,length,radius);
  }
  //<---- Convex Object
  else if (type == dConvexClass)
  {
#if 0
    dsDrawConvex(pos,R,planes,
                 planecount,
                 points,
                 pointcount,
                 polygons);
#else
    dsDrawConvex(pos,R,
                 Sphere_planes,
                 Sphere_planecount,
                 Sphere_points,
                 Sphere_pointcount,
                 Sphere_polygons);
#endif
  }
  //----> Convex Object
  else if (type == dCylinderClass) {
    dReal radius,length;
    dGeomCylinderGetParams (g,&radius,&length);
    dsDrawCylinder (pos,R,length,radius);
  }
  else if (type == dGeomTransformClass) {
    dGeomID g2 = dGeomTransformGetGeom (g);
    const dReal *pos2 = dGeomGetPosition (g2);
    const dReal *R2 = dGeomGetRotation (g2);
    dVector3 actual_pos;
    dMatrix3 actual_R;
    dMultiply0_331 (actual_pos,R,pos2);
    actual_pos[0] += pos[0];
    actual_pos[1] += pos[1];
    actual_pos[2] += pos[2];
    dMultiply0_333 (actual_R,R,R2);
    drawGeom (g2,actual_pos,actual_R,0);
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
  if (!pause) dWorldQuickStep (world,0.02);

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
    for (int j=0; j < GPB; j++) {
      if (i==selected) {
        dsSetColor (0,0.7,1);
      }
      else {
        dsSetColor (1,1,0);
      }
      drawGeom (obj[i].geom[j],0,0,show_aabb);
    }
  }
}


int main (int argc, char **argv)
{
  // init global vars /***/
  maxNumContactsSimulated = 0;
  fp = fopen("output.txt","w");
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
  memset (obj,0,sizeof(obj));

  dThreadingImplementationID threading = dThreadingAllocateMultiThreadedImplementation();
  dThreadingThreadPoolID pool = dThreadingAllocateThreadPool(4, 0, dAllocateFlagBasicData, NULL);
  dThreadingThreadPoolServeMultiThreadedImplementation(pool, threading);
  // dWorldSetStepIslandsProcessingMaxThreadCount(world, 1);
  dWorldSetStepThreadingImplementation(world, dThreadingImplementationGetFunctions(threading), threading);

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
