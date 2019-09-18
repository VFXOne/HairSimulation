#pragma once
// Functions that enable linking a nanogui viewer with projective dynamics

#include "projdyn_types.h"

constexpr int PROJDYN_FPS = 60;
constexpr int PROJDYN_NUM_ITS_INITIAL = 10;

class Viewer;

bool projdyn_setmesh(Viewer* viewer);
bool projdyn_setmesh(Viewer* viewer, bool add_tets);
bool projdyn_update(Viewer* viewer);
bool projdyn_start(Viewer* viewer);
void projdyn_stop();
void projdyn_add_constraints(Viewer* viewer, const std::vector<ProjDyn::ConstraintPtr>& constraints);
void init_projdyn_gui(Viewer* viewer);
void projdyn_grab(const std::vector<size_t>& grabbedVerts, const std::vector<Eigen::Vector3f>& grabPos);
void projdyn_release_grab();
bool projdyn_upload_positions(Viewer*);
