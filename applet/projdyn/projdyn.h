#pragma once

#include "projdyn_types.h"

#include <memory>
#include <iostream>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>

typedef OpenMesh::TriMesh_ArrayKernelT<> OMesh;

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
			m = pos.rows();
			q = pos;
			std::cout << "m = " << m << std::endl;
		}

		void resetPositions() {
			// Reset vertex positions to initial state and set velocities to 0
		}

		bool isInitialized() {
			return isReady;
		}

		bool initializeSystem() {
			// Setup simulation (collect constraints, set up and factorize global step system)

			masses_flat = Vector::Ones(m);
			masses = masses_flat.asDiagonal();
			masses_inv = masses.inverse();
			f_ext << 1, 1, 1;
			f_ext.replicate(m, 1);

			//isReady = true
			return isReady;
		}

		bool step(int num_iterations) {
			// Perform a simulation step -> PD Algo

			Positions s_n = q + h * velocities + h * h * masses_inv * f_ext;
			Positions q_n_1 = s_n;

			size_t step = 0;

			while (step < num_iterations) {
			    step++;
			    //What tf are constraints !?

			}

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
		Positions f_ext;
		size_t m;
		Positions q;
		Positions velocities;
		Vector masses_flat;
		Matrix masses;
		Matrix masses_inv;

		//State variables
		bool isReady;

		//----------------------------------------------

		float ProjectOnConstraintSet(float c_i, Positions q) {
		    return 0;
		}

		float SolveLinearSystem(Positions s_n, Positions p) {
		    return 0;
		}

	};

}
