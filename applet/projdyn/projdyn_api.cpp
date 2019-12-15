#include "projdyn_api.h"
#include <thread>
#include "projdyn.h"
#include <nanogui/slider.h>
#include "viewer.h"

// Global variables used by the callback functions
int projdyn_num_iterations = PROJDYN_NUM_ITS_INITIAL;
ProjDyn::Simulator sim;
Eigen::MatrixXf upload_pos, upload_rods, upload_rods_tan, upload_rods_norm;
std::thread projdyn_thread;
bool projdyn_active = false;
bool default_constraints = true;

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
    int j = 0;
    for (auto v : mesh->vertices()) {
        vertices.row(j) << (ProjDyn::Scalar)mesh->position(v).x,
                (ProjDyn::Scalar)mesh->position(v).y,
                (ProjDyn::Scalar)mesh->position(v).z;
        ++j;
    }

    if (viewer->is_using_rods()) {
        MatrixXf* pos = viewer->getRodsPos();
        std::vector<Positions> rods;
        std::vector<size_t> rod_indices = viewer->getRodIndices();

        for (size_t ind = 0; ind < rod_indices.size(); ind++) {
            size_t rod_index = rod_indices.at(ind);
            size_t next_index = ind == rod_indices.size()-1 ? pos->cols() : rod_indices.at(ind+1);

            Positions rodPos;
            rodPos.resize(next_index - rod_index, 3);

            for (size_t i = rod_index; i < next_index; i++) {
                Vector3f p = pos->col(i);
                rodPos.row(i-rod_index) << p.x(), p.y(), p.z();
            }
            rods.push_back(rodPos);
        }

        sim.setRods(rods);
        sim.setMesh(vertices);
    } else {
        sim.setMesh(vertices);
    }

    return true;
}

void setup_demo_scene(Viewer* viewer) {
    const float radius = 0.8;
    viewer->addRodsOnBall(radius, 5, 2, 0.5);
    default_constraints = false;
}

void init_projdyn_gui(Viewer* viewer) {
    //Simulation panel
	Window* pd_win = new Window(viewer, "Simulation Controls");
	pd_win->setPosition(Vector2i(15, 230));
	pd_win->setLayout(new GroupLayout());

	Button* sim_setup = new Button(pd_win, "Set up demo scene");
	sim_setup->setCallback([viewer]() {
		setup_demo_scene(viewer);
	});

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
	if (!sim.isInitialized())
		if (!sim.initializeSystem(default_constraints))
            return false;

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

	Positions* pos = sim.getPositions();

    // This is done using the 3 x #Verts matrix upload_pos
    size_t num_verts = pos->rows();
    upload_pos.resize(3, num_verts);

    //Copy the transpose of the matrix because different storage order
    for (size_t i = 0; i < num_verts; i++) {
        upload_pos(0, i) = pos->coeff(i, 0);
        upload_pos(1, i) = pos->coeff(i, 1);
        upload_pos(2, i) = pos->coeff(i,2);
    }

    std::vector<ProjDyn::Index> rod_indices;
    Positions* rods_pos = sim.getRodsPositions();
    size_t num_rods = rods_pos->rows();
    upload_rods.resize(3, num_rods);

    for (size_t i = 0; i < num_rods; i++) {
        upload_rods(0, i) = rods_pos->coeff(i, 0);
        upload_rods(1, i) = rods_pos->coeff(i, 1);
        upload_rods(2, i) = rods_pos->coeff(i,2);
    }

    upload_rods_tan.resize(3, num_rods);
    Positions* rods_tan = sim.getRodsTangents();

    for (size_t i = 0; i < num_rods; i++) {
        upload_rods_tan(0, i) = rods_tan->coeff(i, 0);
        upload_rods_tan(1, i) = rods_tan->coeff(i, 1);
        upload_rods_tan(2, i) = rods_tan->coeff(i,2);
    }

    upload_rods_norm.resize(3, num_rods);
    Positions* rods_normals = sim.getRodsNormals();

    for (size_t i = 0; i < num_rods; i++) {
        upload_rods_norm(0, i) = rods_normals->coeff(i, 0);
        upload_rods_norm(1, i) = rods_normals->coeff(i, 1);
        upload_rods_norm(2, i) = rods_normals->coeff(i,2);
    }

    rod_indices = sim.getRodIndices();

	// The matrix is sent to the viewer with this function
    viewer->updateShaderRods(upload_rods, upload_rods_tan, upload_rods_norm, rod_indices);
    viewer->updateShaderVertices(upload_pos);

	return true;
}

void projdyn_grab(const std::vector<ProjDyn::Index>& grabbedVerts, const std::vector<Vector3f>& grabPos) {
	if (!sim.isInitialized() || !projdyn_active) return;

	sim.setGrab(grabbedVerts, grabPos);
}

void projdyn_release_grab() {
	if (!sim.isInitialized()) return;

	sim.releaseGrab();
}
