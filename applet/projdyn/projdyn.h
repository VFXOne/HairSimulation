#pragma once

#include "projdyn_types.h"

#include <memory>
#include <iostream>

namespace ProjDyn {
	typedef Eigen::SimplicialLDLT<SparseMatrix> SparseSolver;

	// The simulator class, which handles objects represented as Eigen matrices
	// and the constraints above to run a Projective Dynamics simulation.
	class Simulator {
	public:

		Simulator()
		{
			// Construct the simulator object
			m_hasGrab = false;
			isReady = false;
		}

		void setMesh(Positions& pos, Triangles& tris) {
			// Pass the positions of the geometry
			// Here you will have to introduce edges instead of triangles
		}

		void resetPositions() {
			// Reset vertex positions to initial state and set velocities to 0
		}

		bool isInitialized() {
			return isReady;
		}

		bool initializeSystem() {
			// Setup simulation (collect constraints, set up and factorize global step system)

			//isReady = true
			return isReady;
		}

		bool step(int num_iterations) {
			// Perform a simulation step -> PD Algo
			/*
			 *s_n = q_n + h * v_n + h*h*M_inv*f_ext;
			 * q_n_1 = s_n
			 *
			 * int step = 0;
			 * while (step < num_iterations) {
			 *  for(int i = 0; i < C.length; i++) {//C -> Constraints
			 *   p[i] = ProjectOnConstraintSet(C[i], q_n_1) //local step
			 *  }
			 *  step++:
			 * }
			 * q_n_1 = SolveLinearSystem(s_n, p) //global step
			 * v_n_1 = (q_n_1 - q_n) / h;
			 */
			return false;
		}

		void setGrab(const std::vector<unsigned long>& grabVerts, const std::vector<Eigen::Vector3f> grabPos) {
			// This method will be called from the viewer if an alt-click has been performed and
			// vertices were grabbed.
			// You can use the list of grabbed vertices and the new positions to force vertices to move
			// in your step() method.

			m_grabVerts = grabVerts;
			m_grabPos = grabPos;
			m_hasGrab = true;
		}

		void releaseGrab() {
			m_hasGrab = false;
		}

	protected:
		// If m_hasGrab is true, vertices with the indices in m_grabVerts will be forced
		// to the positions set in m_grabPos
		bool m_hasGrab = false;
		std::vector<unsigned long> m_grabVerts;
		std::vector<Eigen::Vector3f> m_grabPos;

		const double h = 1.0; //Simulation step size

		//State variables
		bool isReady;

	};

}
