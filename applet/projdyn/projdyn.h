
#include "projdyn_types.h"
#include "Edge.h"
#include "PDConstraint.h"

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
			m = pos.rows();
			q = pos;
			q_init = pos;
			createEdges(tris);
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
			f_ext << 0, 0, -1;
			f_ext.replicate(m, 1);
			A_i = masses.array().sqrt().matrix(); //take the square root element-wise
			B_i = A_i;

			addEdgeSpringConstraints();
			addGroundConstraints();

			factorizeLHS();

			//TODO: Factorize global step matrices

			//isReady = true
			return isReady;
		}

		bool step(int num_iterations) {
			// Perform a simulation step -> PD Algo

			size_t numConstraints = constraints->size();

			Positions s_n = q + h * velocities + h * h * masses_inv * f_ext;
			Positions q_n_1 = s_n;
			std::vector<Positions> p_i(numConstraints);

			size_t step = 0;

			while (step < num_iterations) {
			    step++;
			    //Local step
			    for (size_t c_ind = 0; c_ind < numConstraints; c_ind++) {
			        PDConstraint* c = constraints->at(c_ind);
                    p_i[c_ind] = (c->projectOnConstraintSet(q_n_1));
			    }
			    //Global step
			    /* sparseSolver.compute(A).solve(B); to solve Ax=b
			     * check state with:
			     * if (solver.info() != Success) //On both states (compute and solve)
			    */
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
        std::vector<Edge> edges;
		const double h = 1.0; //Simulation step size

        Positions f_ext;
        size_t m;
        Positions q_init;
        Positions q;
        Positions velocities;
        Vector masses_flat;
        Matrix masses;
        Matrix masses_inv;
        Matrix lhs;
        Matrix A_i;
        Matrix B_i;
		std::vector<PDConstraint*>* constraints;

		//State variables
		bool isReady;

		void createEdges(Triangles& triangles) {
		    for (size_t i = 0; i < triangles.size(); i++) { //iterate over all triangles
                size_t length = triangles.row(i).size(); //Find size of row

                for (size_t t = 0; t < length; t++) { //iterate over all coordinates
                    //Find the two indexes of an edge
		            Index t0 = triangles(i, t);
		            Index t1 = triangles(i, (t + 1) % length);
		            Edge edge = Edge(t0, t1);

		            if (std::find(edges.begin(), edges.end(), edge) == edges.end()) { //Check if the edge already exists
		                edges.push_back(edge);
		            }
		        }
		    }
		}

		void addGroundConstraints() {
            //Iterate over all vertices to add the ground constraint
            for (size_t i = 0; i < q.size(); i++) {
                GroundConstraint new_ground_c = GroundConstraint(m, i, 1);
                constraints->push_back(&new_ground_c);
            }
		}

		void addEdgeSpringConstraints() {
		    //Iterate over all edges to add the spring constraint
            for (auto e : edges) {
                EdgeSpringConstraint new_edge_c = EdgeSpringConstraint(m, &e, 1.0, 0.001, 1.5);
                constraints->push_back(&new_edge_c);
            }
		}

        void factorizeLHS() {
		    //Need to make sure that constraints exists
		    lhs = masses / (h*h);
		    for (size_t i = 0; i < constraints->size(); i++) {
                lhs += constraints->at(i)->getSelectionMatrix() * A_i.transpose() * A_i * constraints->at(i)->getSelectionMatrix();
		    }
        }

	};

}