#include "projdyn_api.h"
#include <thread>
#include "projdyn.h"
#include <nanogui/slider.h>
#include "viewer.h"

// Global variables used by the callback functions
int projdyn_num_iterations = PROJDYN_NUM_ITS_INITIAL;
ProjDyn::Simulator sim;
Eigen::MatrixXf upload_pos;
std::thread projdyn_thread;
bool projdyn_active = false;

// Pass updated vertex positions to the viewer
bool projdyn_setmesh(Viewer* viewer) {
	return projdyn_setmesh(viewer, false);
}

// When a new mesh is loaded in the viewer, update the mesh in the simulator as well
bool projdyn_setmesh(Viewer* viewer, bool add_tets) {
	std::cout << "New mesh was loaded, re-initializing simulation..." << std::endl;

	// Stop the running simulation
	projdyn_stop();

    // Convert surface mesh to Eigen matrices
    surface_mesh::Surface_mesh* mesh = viewer->getMesh();
    ProjDyn::Positions vertices(mesh->n_vertices(), 3);
    ProjDyn::Triangles faces(mesh->n_faces(), 3);
    int j = 0;
    for (auto f : mesh->faces()) {
        int k = 0;
        for (auto v : mesh->vertices(f)) {
            faces(j, k) = (ProjDyn::Index)v.idx();
            ++k;
        }
        ++j;
    }
    j = 0;
    for (auto v : mesh->vertices()) {
        vertices.row(j) << (ProjDyn::Scalar)mesh->position(v).x,
                (ProjDyn::Scalar)mesh->position(v).y,
                (ProjDyn::Scalar)mesh->position(v).z;
        ++j;
    }

    sim.setMesh(vertices, faces);

	return true;
}

void init_projdyn_gui(Viewer* viewer) {
	Window* pd_win = new Window(viewer, "Simulation Controls");
	pd_win->setPosition(Vector2i(15, 230));
	pd_win->setLayout(new GroupLayout());

	Button* runsim_b = new Button(pd_win, "Run Simulation");
	runsim_b->setCallback([viewer]() {
		projdyn_start(viewer);
	});

	Button* stopsim_b = new Button(pd_win, "Stop Simulation");
	stopsim_b->setCallback([]() {
		projdyn_stop();
	});

	Button* reset_b = new Button(pd_win, "Reset Positions");
	reset_b->setCallback([viewer]() {
		bool was_active = projdyn_active;
		projdyn_stop();
		sim.resetPositions();
		if (was_active) {
			projdyn_start(viewer);
		} else {
			projdyn_upload_positions(viewer);
		}
	});

	Label* iterations_label = new Label(pd_win, "Num of Iterations: ");
	IntBox<int>* iterations_box = new IntBox<int>(pd_win, projdyn_num_iterations);
	iterations_box->setEditable(true);
	iterations_box->setCallback([viewer](int num_its) {
		projdyn_num_iterations = num_its;
	});

	viewer->performLayout();
}

bool projdyn_start(Viewer* viewer) {
	projdyn_stop();

	// Make sure the simulator is properly initialized
	if (!sim.isInitialized()) {
		if (!sim.initializeSystem())
			return false;
	}

	// Create a thread that runs the simulation
	// It calls a function that triggers a time-step every 1000/PROJDYN_FPS milliseconds
	projdyn_active = true;
	projdyn_thread = std::thread(
		[viewer]() {
		std::chrono::milliseconds time(1000 / PROJDYN_FPS);
		while (projdyn_active) {
			std::this_thread::sleep_for(time);
			projdyn_update(viewer);
			glfwPostEmptyEvent();
		}
	}
	);

	return true;
}

void projdyn_stop() {
	if (projdyn_active) {
		projdyn_active = false;
		projdyn_thread.join();
	}
}

// Will be called every frame:
// Performs a time step and updates the positions that are drawn in the shader window
bool projdyn_update(Viewer* viewer) {
	if (!sim.isInitialized()) return false;

	// Simulate one time step
	sim.step(projdyn_num_iterations);

	return projdyn_upload_positions(viewer);
}

// Extract positions, convert them to column-wise tripples of floats and
// upload them to the OpenGL buffer
bool projdyn_upload_positions(Viewer* viewer) {
	// In this function you need to extract the vertex positions of the simulation
	// and send them to the viewer.

	// This is done using the 3 x #Verts matrix upload_pos
	size_t num_verts = 1;
	upload_pos.resize(3, num_verts);
	upload_pos(0, 0) = 1.;
	upload_pos(1, 0) = 1.;
	upload_pos(2, 0) = 1.;

	// The matrix is sent to the viewer with this function
	viewer->updateShaderVertices(upload_pos);

	return true;
}

void projdyn_grab(const std::vector<size_t>& grabbedVerts, const std::vector<Vector3f>& grabPos) {
	if (!sim.isInitialized() || !projdyn_active) return;

	sim.setGrab(grabbedVerts, grabPos);
}

void projdyn_release_grab() {
	if (!sim.isInitialized()) return;

	sim.releaseGrab();
}
