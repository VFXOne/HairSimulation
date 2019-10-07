
#include "projdyn_types.h"
#include "Edge.h"
#include "PDConstraint.h"
#include "Quaternion.h"

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
			m = pos.rows();
            m_positions = Positions(pos);
            m_positions_init = Positions(pos);
			m_velocities.setZero(pos.rows(), pos.cols());
			createEdges(tris);
			//TODO: Reset the correct variables (dimensions) to fit new mesh
		}

		void setRods(size_t length, Positions& pos) {
		    assert(pos.rows() % length == 0);
		    cr_m = pos.rows();
		    cr_positions_init = Positions(pos);
		    cr_positions = Positions(pos);

		    cr_velocities.setZero(pos.rows(), pos.cols());
		}

		void resetPositions() {
			// Reset vertex positions to initial state and set m_velocities to 0
			m_positions = Positions(m_positions_init);
			m_velocities.setZero(m_positions.rows(), m_positions.cols());
		}

		bool isInitialized() {
			return is_ready;
		}

        bool initializeSystem() {
			// Setup simulation (collect constraints, set up and factorize global step system)
			if (is_simulating) {
			    return false;
			}

            m_masses_flat = Vector::Ones(m) * 0.01;
            m_masses = m_masses_flat.asDiagonal();
            m_masses_inv = m_masses.cwiseInverse();

			f_ext.setZero(m, 3);
			//Set gravity on z axis m times
			for (size_t i = 0; i < m; i++) {
			    f_ext.row(i) << 0, 0, gravity;
			}

			//addEdgeSpringConstraints();
			addGroundConstraints();

            lhs = m_masses / (h * h);
            for (size_t i = 0; i < m_constraints.size(); i++) {
                lhs += m_constraints.at(i)->getSelectionMatrixWeighted().transpose() * m_constraints.at(i)->getSelectionMatrix();
            }

			rhs.resize(m, 3);
			rhs.setZero();

			//TODO: Factorize global step matrices

			is_ready = true;
			return is_ready;
		}

		bool step(int num_iterations) {
			// Perform a simulation step
            is_simulating = true;

            size_t numConstraints = m_constraints.size();

            Positions s_n = m_positions + h * m_velocities + m_masses_inv * f_ext * h * h;
            Positions positions_updated = Positions(s_n);
            std::vector<Positions> projections(numConstraints);

            size_t step = 0;

            m_solver.compute(lhs);

            while (step < num_iterations) {
                step++;
                /****************
                 ** Local step **
                 ****************/
                for (size_t c_ind = 0; c_ind < numConstraints; c_ind++) {
                    PDConstraint* c = m_constraints.at(c_ind);
                    projections[c_ind] = c->projectOnConstraintSet(positions_updated);
                    Positions p = c->projectOnConstraintSet(positions_updated);
                }

                /*****************
                 ** Global step **
                 *****************/
			    if (m_solver.info() != Eigen::Success) {
			        std::cout << "Unable to factorize left hand side" << std::endl;
			        return false;
			    }

                //rhs = m_masses_flat * s_n / (h*h);
                rhs.setZero();
                for (size_t v = 0; v < m; v++) {
                    for (int d = 0; d < 3; d++) {
                        rhs.coeffRef(v, d) += m_masses_flat(v) * s_n(v, d) / (h * h);
                    }
                }

			    //Factorize right hand side of the system
			    //TODO: Parallelize this
			    for (size_t i = 0; i < numConstraints; i++) {
                    SparseMatrix S_i_T = m_constraints.at(i)->getSelectionMatrixWeighted().transpose();
                    for (size_t d = 0; d < 3; d++) {
                        for (int k = 0; k < S_i_T.cols(); ++k) {
                            for (SparseMatrix::InnerIterator it(S_i_T, k); it; ++it) {
                                rhs.coeffRef(it.row(), d) += it.value() * projections.at(i)(0, d);
                            }
                        }
                    }
                }

                for (size_t d = 0; d < 3; d++) {
                    positions_updated.col(d) = m_solver.solve(rhs.col(d));
                }

                is_simulating = false;

                if (m_solver.info() == Eigen::Success) {
                    //Update positions and m_velocities
                    m_velocities = (positions_updated - m_positions) / h;
                    m_positions = Positions(positions_updated);
                    return true;
                } else {
                    std::cout << "Unable to solve constraints" << std::endl;
                    return false;
                }
            }

			return false;
		}

        bool CR_step(int num_iterations) {
		    /*
		     * s_x_t = x_t + h*v_t + h*h*M_inv*f_ext
		     * s_w_t = w_t + h*J_inv*(tau - w_t.cross(J*w_t))
		     * s_u_t = u_t + 1/2*h*(u_t {quaternion mult} s_w_t)
		     * s_t = [s_x_t, s_u_t].tranpose() //Optional because we might distinct the two (positions and quaternions)
		     * q_t_1 = s_t // original algo
		     * q_x_t_1 = s_x_t //Our implementation
		     * q_u_t_1 = s_u_t //Our implementation
		     *
		     * for num_iterations times:
		     *     for all constraints i:
		     *         p_i = ProjectOnConstraintSet(C_i, q_t_1) //Original algo
		     *         p_i = ProjectOnConstraintSet(C_i, q_x_t_1, q_u_t_1) //Our method
		     *     q_t_1 = SolveLinearSystem(s_t, p1, p2, ...)
		     * v_t_1 = 1/h*(x_t_1 - x_t) //Original algo
		     * v_t_1 = 1/h*(q_x_t_1 - x_t) //Our method
		     * w_t_1 = 2/h*(u_t {quaternion mult} u_t_1) //Original algo
		     * w_t_1 = 2/h*(u_t {quaternion mult} q_u_t_1) //Our method
		     */
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
		    return &m_positions;
		}

	protected:

        // If m_hasGrab is true, vertices with the indices in m_grabVerts will be forced
		// to the positions set in m_grabPos
		bool m_hasGrab = false;
        std::vector<unsigned long> m_grabVerts;
        std::vector<Eigen::Vector3f> m_grabPos;
        std::vector<Edge> edges;
		const double h = 0.005; //Simulation step size
		const float gravity = -1.0;

        Positions f_ext;

        /******************
         * Mesh variables *
         ******************/
        size_t m;
        Positions m_positions_init;
        Positions m_positions;
        Positions m_velocities;
        Vector m_masses_flat;
        SparseMatrix m_masses;
        SparseMatrix m_masses_inv;
        SparseMatrix lhs;
        SparseMatrix rhs;
		std::vector<PDConstraint*> m_constraints;
		SparseSolver m_solver;

        /***************************
         * Cosserat Rods variables *
         ***************************/
         size_t cr_m;
         Positions cr_positions_init;
         Positions cr_positions;
         Positions cr_velocities;

         std::vector<Quaternion<double>> cr_quaternions;

         std::vector<PDConstraint*> cr_constraints;
         SparseSolver cr_solver;

		//State variables
		bool is_ready;
		bool is_simulating;

		void createEdges(Triangles& triangles) {
		    for (size_t i = 0; i < triangles.rows(); i++) { //iterate over all triangles
                size_t length = triangles.row(i).size();
                for (size_t t = 0; t < length; t++) { //iterate over all coordinates
                    //Find the two indexes of an edge
		            Index t0 = triangles(i, t);
		            Index t1 = triangles(i, (t + 1) % length);
		            //TODO: Use pointers
		            Edge edge = Edge(t0, t1);

                    if (std::find(edges.begin(), edges.end(), edge) == edges.end()) { //Check if the edge already exists
                        edges.push_back(edge);
                    }
		        }
		    }

		}

		void addGroundConstraints() {
            //Iterate over all vertices to add the ground constraint
            for (size_t i = 0; i < m_positions.rows(); i++) {
                GroundConstraint* new_ground_c = new GroundConstraint(m, i, 1.0);
                m_constraints.push_back(new_ground_c);
            }
		}

		void addEdgeSpringConstraints() {
		    //Iterate over all edges to add the spring constraint
            for (auto e : edges) {
                EdgeSpringConstraint* new_edge_c = new EdgeSpringConstraint(m, e, 1.0, 0.001, 1.5);
                m_constraints.push_back(new_edge_c);
            }
		}

	};

}