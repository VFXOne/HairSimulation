
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
            use_cosserat_rods = false;
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
		    use_cosserat_rods = true;
		    assert(pos.rows() % length == 0);
            cr_num_positions = pos.rows();
		    cr_num_quaternions = cr_num_positions - 1;
		    cr_size = 3 * cr_num_positions + 4 * cr_num_quaternions;
		    cr_positions_init = Positions(pos);
		    cr_positions = Positions(pos);
		    //TODO: Compute initial quaternions

		    cr_velocities.setZero(cr_num_positions, 3);
		    cr_angular_velocities.setZero(cr_num_positions);
		}

		void resetPositions() {
			// Reset vertex positions to initial state and set m_velocities to 0
			m_positions = Positions(m_positions_init);
			m_velocities.setZero();

			if (use_cosserat_rods) {
                cr_positions = Positions(cr_positions_init);
                cr_velocities.setZero(cr_num_positions, 3);
                cr_angular_velocities.setZero(cr_num_positions, 3);
                //TODO: Recompute correct orientations
            }

        }

		bool isInitialized() {
			return is_ready;
		}

        bool initializeSystem() {
			// Setup simulation (collect constraints, set up and factorize global step system)
			if (is_simulating) {
			    return false;
			}


            /****************
             * PD Variables *
             ****************/
            const float unit_weight = 0.01;
            m_masses_flat = Vector::Ones(m) * unit_weight;
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
                //lhs += m_constraints.at(i)->computeLHS();
                auto c = m_constraints.at(i);
                lhs += c->getSelectionMatrixWeighted().transpose() * c->getSelectionMatrix();
            }
            m_solver.compute(lhs);

            rhs.resize(m, 3);
            rhs.setZero();

            //TODO: Initialize cosserat rods variables: M_star, torques, angular velocities,
            if (use_cosserat_rods) {
                /****************
                 * CR Variables *
                 ****************/
                Vector masses_flat = Vector::Ones(m) * cr_unit_weight;
                cr_masses = masses_flat.asDiagonal();
                cr_masses_inv = cr_masses.cwiseInverse();

                const float j_val = M_PI*cr_radius/4;
                Vector J_flat_vec;
                J_flat_vec << j_val, j_val, j_val*2;
                J_flat_vec *= cr_segment_length * cr_density;
                cr_J_vec = stackDiagonal(J_flat_vec, cr_num_positions);
                cr_J_vec_inv = cr_J_vec.cwiseInverse();

                Vector J_flat_quat;
                J_flat_quat << 0, j_val, j_val, j_val*2;
                J_flat_quat *= cr_segment_length * cr_density;
                cr_J_quat = stackDiagonal(J_flat_quat, cr_num_quaternions);
                cr_J_quat_inv = cr_J_quat.cwiseInverse();
            }

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

            while (step < num_iterations) {
                step++;
                /****************
                 ** Local step **
                 ****************/
                for (size_t c_ind = 0; c_ind < numConstraints; c_ind++) {
                    PDConstraint* c = m_constraints.at(c_ind);
                    Positions p_i = c->projectOnConstraintSet(positions_updated);
                    projections[c_ind] = p_i;
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
			    //TODO: Parallelize this badboi
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

			    //TODO: parallelize
                for (size_t d = 0; d < 3; d++) {
                    positions_updated.col(d) = m_solver.solve(rhs.col(d));
                }

                is_simulating = false;

                if (m_solver.info() == Eigen::Success) {
                    //Update positions and m_velocities
                    m_velocities = (positions_updated - m_positions) / h;
                    m_positions = Positions(positions_updated);
                } else {
                    std::cout << "Unable to solve constraints" << std::endl;
                    return false;
                }
            }

			return true;
		}

        bool CR_step(int num_iterations) {

		    is_simulating = true;

		    const size_t num_constraints = cr_constraints.size(); //TODO: global variable
		    std::vector<Vector> projections(num_constraints);

		    Positions s_x = cr_positions + h*cr_velocities + h*h * cr_masses_inv * f_ext;
		    Positions cross = Positions(cr_angular_velocities);
            vectorCross(cr_angular_velocities, cr_J_vec * cr_angular_velocities);

		    Positions s_w = cr_angular_velocities + h*cr_J_vec_inv * (cr_torques - cross);
            Orientations s_w_u = pos2quat(s_w);
		    Orientations s_u = pos2quat(quat2pos(cr_orientations) + h/2 * quat2pos(cr_orientations * s_w_u));

		    Vector s_t = posQuatConcat(s_x, s_u);
            cr_q_t = Vector(s_t);

		    size_t step = 0;
		    while (step < num_iterations) {
                step++;
                /****************
                 ** Local step **
                 ****************/
		        for (size_t i = 0; i < num_constraints; i++) {
		            projections[i] = cr_constraints.at(i)->projectOnConstraintSet(cr_q_t);
		        }

                /*****************
                 ** Global step **
                 *****************/
		        //Compute left-hand side of system (can be done earlier)
		        cr_lhs.resize(cr_size, cr_size);
		        cr_lhs.setZero();
                cr_lhs = cr_M_star / (h*h);

		        for (size_t i = 0; i < num_constraints; i++) {
		            auto c = cr_constraints.at(i);
		            auto Ai = c->getAiMatrix();
                    cr_lhs += c->getSelectionMatrixWeighted().transpose()
                            * Ai.transpose() * Ai * c->getSelectionMatrix();
		        }

		        //Compute right-hand side
		        cr_rhs = cr_M_star * s_t / (h*h);
		        for (size_t i = 0; i < num_constraints; i++) {
		            auto c = cr_constraints.at(i);
		            cr_rhs += c->getSelectionMatrixWeighted().transpose()
		                    * c->getAiMatrix().transpose() * c->getAiMatrix()
		                    * projections.at(i);
		        }

		    }

            cr_solver.compute(cr_lhs);
            cr_q_t = cr_solver.solve(cr_rhs);

            is_simulating = false;

            if (cr_solver.info() == Eigen::Success) {
                //Update velocities and angular velocities
                Orientations new_quat;
                Positions new_pos;
                separatePosQuat(&cr_q_t, new_pos, new_quat, cr_num_positions);
                cr_velocities = (new_pos - cr_positions) / h;
                cr_positions = Positions(new_pos);
                cr_angular_velocities = 2/h * quat2pos(conjugateQuat(cr_orientations) * new_quat);
            } else {
                std::cout << "Unable to solve constraints" << std::endl;
                return false;
            }

		    /*
		     * s_x_t = x_t + h*v_t + h*h*M_inv*f_ext
		     * s_w_t = w_t + h*J_inv*(tau - w_t.cross(J*w_t)) // tau = torques
		     * s_u_t = u_t + 1/2*h*(u_t {quaternion mult} s_w_t)
		     * s_t = [s_x_t, s_u_t].tranpose() //Optional because we might distinct the two (positions and quaternions)
		     * q_t_1 = s_t // original algo
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
		const double h = 0.05; //Simulation step size
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
        const float cr_unit_weight = 0.01;
        const float cr_radius = 0.1;
        const float cr_density = 0.01;
        float cr_segment_length;
        size_t cr_num_positions; //Number of positions
        size_t cr_num_quaternions;
        size_t cr_size; //Full size of the flat variables
        Positions cr_positions_init;
        Positions cr_positions;
        Vector cr_q_t; //The flat vector that holds every variables (positions and quaternions)
        Positions cr_velocities;
        Positions cr_angular_velocities;
        Positions cr_torques;

        SparseMatrix cr_masses;
        SparseMatrix cr_masses_inv;
        SparseMatrix cr_J_vec;
        SparseMatrix cr_J_vec_inv;
        SparseMatrix cr_J_quat;
        SparseMatrix cr_J_quat_inv;
        SparseMatrixRM cr_M_star; //Là on utilise J_quat

        Orientations cr_orientations;

        std::vector<CRConstraint*> cr_constraints;
        SparseSolver cr_solver;
        Vector cr_rhs;
        SparseMatrix cr_lhs;

		//State variables
		bool is_ready;
		bool is_simulating;
		bool use_cosserat_rods;

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

		//Create imaginary parts of quaternions from vectors. Real part is 0
		static Orientations pos2quat(Positions p) {
		    Orientations u(p.rows());
		    for (size_t i = 0; i < p.rows(); i++) {
		        Vector3 pos = p.row(i);
                Quaternion quat(0, pos[0], pos[1], pos[2]);
                u[i] = quat.normalized();
		    }
		    return u;
		}

		//Creates vectors from the imaginary part of the quaternions
		static Positions quat2pos(Orientations quat) {
		    Positions pos;
		    pos.resize(quat.size(), 3);
		    for (size_t i = 0; i < quat.size(); i++) {
		        auto q = quat[i].normalized();
		        //We asssume the real part is null
		        assert(q.w() == 0.0);
		        pos.row(i) << q.x(), q.y(), q.z();
		    }
		}

		//Can we normalize in place ?
		static Orientations conjugateQuat(Orientations quat) {
		    Orientations new_quat;
		    new_quat.resize(quat.size());
		    for (size_t i = 0; i < quat.size(); i++) {
		        new_quat[i] = quat[i].normalized();
		    }
		    return new_quat;
		}

		static Positions vectorCross(Positions a, Positions b) {
		    Positions cross;
		    assert(a.cols() == b.cols() && a.rows() == b.rows());
		    cross.setZero(a.rows(), a.cols());

		    for (size_t i = 0; i < cross.rows(); i++) {
		        cross.row(i) = a.row(i).cross(b.row(i));
		    }
		    return cross;
		}

		static Vector posQuatConcat(Positions pos, Orientations u) {
		    size_t m = pos.rows();
		    size_t dim = 3*pos.rows() + 4*u.rows();
		    Vector fp(dim);

		    for (size_t i = 0; i < m; i++) {
		        auto vec = pos.row(i);
		        fp(i*3) = vec.x();
		        fp(i*3+1) = vec.y();
		        fp(i*3+2) = vec.z();
		    }

		    for (size_t i = 0; i < u.rows(); i++) {
		        auto quat = u(i);
		        fp(m + i*4) = quat.w();
		        fp(m + i*4+1) = quat.x();
		        fp(m + i*4+2) = quat.y();
		        fp(m + i*4+3) = quat.z();
		    }

		    return fp;
		}

		static void separatePosQuat(Vector* concat, Positions& pos, Orientations& quat, size_t separation) {
		    size_t size = concat->size();
		    size_t pos_height = (size_t)(separation/3);
		    size_t num_quat = (size_t)((size-separation)/4);

		    pos.resize(pos_height, 3);
		    pos.setZero();

		    quat.resize(num_quat);

		    for (size_t i = 0; i < pos_height; i++) {
		        pos.row(i) << concat[i*3], concat[i*3+1], concat[i*3+2];
		    }
		    std::cout << "Concat size:" << size << " pos after: " << pos_height << "quat after: " << num_quat << std::endl;

		    for (size_t i = 0; i < num_quat; i++) {
		        size_t index = separation + i*4;
		        quat.coeffRef(i) = Eigen::Quaterniond(concat->coeff(index),
                                                    concat->coeff(index+1),
                                                    concat->coeff(index+2),
                                                    concat->coeff(index+3));
		    }
		}

		static SparseMatrix stackDiagonal(Vector x, size_t times) {
		    size_t height = x.size();
		    SparseMatrix concat;
		    concat.resize(times*height, times*height);
		    concat.setZero();
		    for (size_t i = 0; i < times; i++) {
		        for (size_t j = 0; j < height; j++) {
		            concat.coeffRef(i*height + j, i*height + j) = x.coeff(j);
		        }
		    }
		    return concat;
		}

	};

}