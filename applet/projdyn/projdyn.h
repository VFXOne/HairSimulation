
#include "projdyn_types.h"
#include "Edge.h"
#include "PDConstraint.h"

#include<Eigen/Geometry>

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
            is_ready = false;
            is_simulating = false;
            m = 0;
		}

		void setMesh(Positions& pos, Triangles& tris) {
			// Pass the positions of the geometry
			// Here you will have to introduce edges instead of triangles
			m = pos.rows();
            positions = Positions(pos);
            positions_init = Positions(pos);
			velocities.setZero(pos.rows(), pos.cols());
			createEdges(tris);
		}

		void resetPositions() {
			// Reset vertex positions to initial state and set velocities to 0
			positions = Positions(positions_init);
			velocities.setZero(positions.rows(), positions.cols());
		}

		bool isInitialized() {
			return is_ready;
		}

        bool initializeSystem() {
			// Setup simulation (collect constraints, set up and factorize global step system)
			if (is_simulating) {
			    return false;
			}

			masses_flat = Vector::Ones(m) * 0.01;
			masses = masses_flat.asDiagonal();
			masses_inv = masses.cwiseInverse();

			f_ext.setZero(m, 3);
			//Set gravity on z axis m times
			for (size_t i = 0; i < m; i++) {
			    f_ext.row(i) << 0, 0, gravity;
			}

			//addEdgeSpringConstraints();
			//addGroundConstraints();

            lhs = masses / (h*h);
            for (size_t i = 0; i < constraints.size(); i++) {
                lhs += constraints.at(i)->getSelectionMatrixWeighted().transpose() * constraints.at(i)->getSelectionMatrix();
            }

			rhs.resize(m, 3);
			rhs.setZero();

			//TODO: Factorize global step matrices

			is_ready = true;
			return is_ready;
		}

		bool step(int num_iterations) {
			// Perform a simulation step -> PD Algo
            is_simulating = true;

            size_t numConstraints = constraints.size();

            Positions s_n = positions + h * velocities + masses_inv * f_ext * h * h;
            Positions positions_updated = Positions(s_n);
            std::vector<Positions> projections(numConstraints);

            size_t step = 0;

            solver.compute(lhs);

            while (step < num_iterations) {
                if (true) {
                    return true;
                }
                step++;
                /****************
                 ***Local step***
                 ****************/
                for (size_t c_ind = 0; c_ind < numConstraints; c_ind++) {
                    PDConstraint* c = constraints.at(c_ind);
                    projections[c_ind] = c->projectOnConstraintSet(&positions_updated);
                }

                /*****************
                 ***Global step***
                 *****************/
			    if (solver.info() != Eigen::Success) {
			        std::cout << "Unable to factorize left hand side" << std::endl;
			        return false;
			    }

                //rhs = masses_flat * s_n / (h*h);
                rhs.setZero();
                for (size_t v = 0; v < m; v++) {
                    for (int d = 0; d < 3; d++) {
                        rhs.coeffRef(v, d) += masses_flat(v) * s_n(v, d) / (h*h);
                    }
                }

			    //Factorize right hand side of the system
			    //TODO: Parallelize this
			    for (size_t i = 0; i < numConstraints; i++) {
                        SparseMatrix S_i_T = constraints.at(i)->getSelectionMatrixWeighted().transpose();
                        //rhs.col(d) += S_i_T * projections[i].col(d);
                    for (size_t d = 0; d < 3; d++) {
                        for (int k = 0; k < S_i_T.cols(); ++k) {
                            //for (SparseMatrix::InnerIterator it(S_i_T, k); it; ++it) {
                            //    rhs.coeffRef(it.row(), d) += it.value() * projections[i](it.row(), d);
                            //}
                            for (size_t j = 0; j < S_i_T.rows(); j++) {
                                rhs.coeffRef(j, d) += S_i_T.coeff(j, k) * projections[i](j, k);
                            }
                        }
                    }
                }

                for (size_t d = 0; d < 3; d++) {
                    positions_updated.col(d) = solver.solve(rhs.col(d));
                }

                is_simulating = false;

                if (solver.info() == Eigen::Success) {
                    //Update positions and velocities
                    velocities = (positions_updated - positions) / h;
                    positions = Positions(positions_updated);
                    return true;
                } else {
                    std::cout << "Unable to solve constraints" << std::endl;
                    return false;
                }
            }


			/*
			 *s_n = q_n + h * v_n + h*h*M_inv*f_ext;
			 * positions_updated = s_n
			 *
			 * int step = 0;
			 * while (step < num_iterations) {
			 *  for(int i = 0; i < C.length; i++) {//C -> Constraints
			 *   p[i] = ProjectOnConstraintSet(C[i], positions_updated) //local step
			 *  }
			 *  step++:
			 * }
			 * positions_updated = SolveLinearSystem(s_n, p) //global step
			 * v_n_1 = (positions_updated - q_n) / h;
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

		Positions* getPositions() {
		    return &positions;
		}

	protected:

        // If m_hasGrab is true, vertices with the indices in m_grabVerts will be forced
		// to the positions set in m_grabPos
		bool m_hasGrab = false;
        std::vector<unsigned long> m_grabVerts;
        std::vector<Eigen::Vector3f> m_grabPos;
        std::vector<Edge> edges;
		const double h = 1.0; //Simulation step size
		const float gravity = -1.0;

        Positions f_ext;
        size_t m;
        Positions positions_init;
        Positions positions;
        Positions velocities;
        Vector masses_flat;
        SparseMatrix masses;
        SparseMatrix masses_inv;
        SparseMatrix lhs;
        SparseMatrix rhs;
		std::vector<PDConstraint*> constraints;
		SparseSolver solver;

		//State variables
		bool is_ready;
		bool is_simulating;

		void createEdges(Triangles& triangles) {
		    for (size_t i = 0; i < triangles.rows(); i++) { //iterate over all triangles
                size_t length = triangles.rows();
                for (size_t t = 0; t < length; t++) { //iterate over all coordinates
                    //Find the two indexes of an edge
		            Index t0 = triangles(i, t);
		            Index t1 = triangles(i, (t + 1) % length);
		            //TODO: Use pointers
		            Edge edge = Edge(t0, t1);

                    edges.push_back(edge);
                    //if (std::find(edges.begin(), edges.end(), edge) == edges.end()) { //Check if the edge already exists
		            //}
		        }
		    }

		}

		void addGroundConstraints() {
            //Iterate over all vertices to add the ground constraint
            for (size_t i = 0; i < positions.rows(); i++) {
                GroundConstraint* new_ground_c = new GroundConstraint(m, i, 1.0);
                constraints.push_back(new_ground_c);
            }
		}

		void addEdgeSpringConstraints() {
		    //Iterate over all edges to add the spring constraint
            for (auto e : edges) {
                EdgeSpringConstraint* new_edge_c = new EdgeSpringConstraint(m, e, 1.0, 0.001, 1.5);
                constraints.push_back(new_edge_c);
            }
		}

	};

}