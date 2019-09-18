//=============================================================================
//
// Skeleton framework for a projective dynamics simulation
//
//-----------------------------------------------------------------------------

#include "viewer.h"
#include "projdyn_api.h"
#include <thread>
#include <chrono>

// Entry point for the program
int main(int argc, char** argv) {
	try {
		// Initialize the GUI engine
		nanogui::init();

		// Create a new mesh viewer app, which will add a screen to nanogui
		// The callback function is triggered when loading a new mesh and (re-)initializes the
		// projective dynamics simulator with the new vertices, faces and tetrahedrons
		nanogui::ref<Viewer> app = new Viewer("Projective Dynamics", nullptr, projdyn_setmesh);
		app->setGrabCallbacks(projdyn_grab, projdyn_release_grab);

		init_projdyn_gui(app);
		app->drawAll();
		app->setVisible(true);

		// Start the main loop: keeps calling drawContents() of the viewer, until the window is closed
		nanogui::mainloop(-1);

		// Clean-up
		projdyn_stop();
		nanogui::shutdown();

		// Error handling:
	}
	catch (const std::runtime_error& e) {
		std::string error_msg = std::string("Caught a fatal error: ") + std::string(e.what());
#if defined(_WIN32)
		MessageBoxA(nullptr, error_msg.c_str(), NULL, MB_ICONERROR | MB_OK);
#else
		std::cerr << error_msg << std::endl;
#endif
		return -1;
	}

	return 0;
}
