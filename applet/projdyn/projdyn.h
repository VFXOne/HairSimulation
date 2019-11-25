#define _USE_MATH_DEFINES
#include <math.h>

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

		void setRods(const std::vector<Positions> rods) {
		    assert(rods.size() >= 1); //At least one rod
            use_cosserat_rods = true;

            size_t rod_index = 0;
            for (size_t i = 0; i < rods.size(); i++) {
                rod_indices.push_back(rod_index);
                Positions pos = rods.at(i);
                assert(pos.rows() >= 2); //At least one segment in a rod

                cr_segments_length.push_back((pos.row(1) - pos.row(0)).norm()); //Assuming all segments have same length
                rod_index += pos.rows();
                cr_num_positions += pos.rows();

                Vector new_pos = pos2flatVec(pos);
                Vector temp_pos(cr_positions_init.size() + new_pos.size());
                temp_pos << cr_positions_init, new_pos;
                cr_positions_init = temp_pos;
            }

            cr_num_flat_pos = 3 * (cr_num_positions - rod_indices.size());
            cr_num_coord = 3 * cr_num_positions;
            cr_num_quaternions = cr_num_positions - rod_indices.size();
            cr_size = 3 * cr_num_positions + 4 * cr_num_quaternions;
            cr_positions = Vector(cr_positions_init);

		    cr_velocities.setZero(cr_num_coord);
		    cr_angular_velocities.setZero(cr_num_flat_pos);
		}

		void resetPositions() {
			// Reset vertex positions to initial state and set m_velocities to 0
			m_positions = Positions(m_positions_init);
			m_velocities.setZero();

			if (use_cosserat_rods) {
                cr_positions = cr_positions_init;
                cr_velocities.setZero(3 * cr_num_positions);
                cr_angular_velocities.setZero(cr_num_flat_pos);
                //TODO: Recompute correct orientations
            }

        }

		bool isInitialized() {
			return is_ready;
		}

		bool isUsingCR() {
		    return use_cosserat_rods;
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

            m_f_ext.setZero(m, 3);
            //Set gravity on z axis m times
            for (size_t i = 0; i < m; i++) {
                m_f_ext.row(i) << 0, gravity, 0;
            }

            //addEdgeSpringConstraints();
            //addGroundConstraints();

            lhs = m_masses / (h * h);
            for (auto c: m_constraints) {
                lhs += c->getSelectionMatrixWeighted().transpose() * c->getSelectionMatrix();
            }
            m_solver.compute(lhs);

            rhs.resize(m, 3);
            rhs.setZero();

            if (use_cosserat_rods) {
                /****************
                 * CR Variables *
                 ****************/

                //addFixedPosConstraint();
                addSSConstraints();
                addBTConstraints();

                cr_f_ext.setZero(cr_num_coord);
                //Set gravity on z axis m times
                for (size_t i = 0; i < cr_num_coord; i++) {
                    cr_f_ext.coeffRef(i) = 0;
                    cr_f_ext.coeffRef(++i) = 0; //gravity;
                    cr_f_ext.coeffRef(++i) = 10;
                }

                cr_orientations = quatFromPos(cr_positions);
                cr_torques.setZero(cr_num_flat_pos);

                Vector masses_flat = Vector::Ones(cr_num_coord) * cr_density;
                cr_masses = masses_flat.asDiagonal();
                cr_masses_inv = cr_masses.cwiseInverse();

                const float j_val = M_PI*cr_radius/4;
                Vector J_flat_vec;
                J_flat_vec.setZero(3);
                J_flat_vec << j_val, j_val, j_val*2;
                J_flat_vec *= cr_density;
                cr_J_vec = stackDiagonal(J_flat_vec, cr_num_positions - rod_indices.size());
                Vector J_flat_vec_inv;
                J_flat_vec_inv.setZero(3);
                J_flat_vec_inv << 1/j_val, 1/j_val, 1/(j_val*2);
                J_flat_vec_inv /= cr_density;
                cr_J_vec_inv = stackDiagonal(J_flat_vec_inv, cr_num_positions - rod_indices.size());

                Vector J_flat_quat;
                J_flat_quat.setZero(4);
                J_flat_quat << 0, j_val, j_val, j_val*2;
                J_flat_quat *= cr_density;
                cr_J_quat = stackDiagonal(J_flat_quat, cr_num_quaternions);
                Vector J_flat_quat_inv;
                J_flat_quat_inv.setZero(4);
                J_flat_quat_inv << 0, 1/j_val, 1/j_val, 1/(j_val*2);
                J_flat_quat_inv /= cr_density;
                cr_J_quat_inv = stackDiagonal(J_flat_quat_inv, cr_num_quaternions);

                //Scale by corresponding segment length
                for (size_t i = 0; i < rod_indices.size(); i++) {
                    Index rod_index = rod_indices.at(i);
                    Index next_index = i == rod_indices.size()-1 ? cr_num_positions : rod_indices.at(i+1);
                    for (size_t j = rod_index; j < next_index; j++) {
                        float seg_len = cr_segments_length.at(i);
                        cr_J_vec.row(j) *= seg_len;
                        cr_J_vec_inv.row(j) /= seg_len;
                        cr_J_quat.row(j) *= seg_len;
                        cr_J_quat_inv.row(j) /= seg_len;
                    }
                }

                cr_M_star.resize(cr_size, cr_size);
                cr_M_star.setZero();
                //Build M_star matrix
                for (size_t i = 0; i < cr_num_positions; i++) {
                    cr_M_star.row(i) = cr_masses.row(i);
                }
                for (size_t i = 0; i < cr_num_quaternions; i++) {
                    cr_M_star.row(i + cr_num_positions) = cr_J_quat.row(i);
                }
            }

			is_ready = true;
			return is_ready;
		}

		bool step(int num_iterations) {
			// Perform a simulation step
            is_simulating = true;

            size_t numConstraints = m_constraints.size();

            Positions s_n = m_positions + h * m_velocities + m_masses_inv * m_f_ext * h * h;
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
                        for (size_t k = 0; k < S_i_T.cols(); ++k) {
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
                    std::cout << "Unable to solve PD constraints" << std::endl;
                    return false;
                }
            }

            if (use_cosserat_rods) {
                CR_step(num_iterations);
            }

			return true;
		}

        bool CR_step(int num_iterations) {

		    is_simulating = true;

		    const size_t num_constraints = cr_constraints.size(); //TODO: global variable
		    std::vector<Vector> projections(num_constraints);

		    Vector s_x = cr_positions + h*cr_velocities + h*h * cr_masses_inv * cr_f_ext;
            Vector cross = vectorCross(cr_angular_velocities, cr_J_vec * cr_angular_velocities);

		    Vector s_w = cr_angular_velocities + h*cr_J_vec_inv * (cr_torques - cross);
            Orientations s_w_u = pos2quat(s_w);
		    Orientations s_u = pos2quat(quat2pos(cr_orientations) + h/2 * quat2pos(cr_orientations * s_w_u));

//		    std::cout << "cr_positions: " << cr_positions << std::endl;
//		    std::cout << "cr_orientations: " << quat2pos(cr_orientations) << std::endl;
//		    std::cout << "s_u: " << quat2pos(s_u) << std::endl;
//		    std::cout << "s_x: " << s_x << std::endl;

		    Vector s_t = posQuatConcat(s_x, s_u);
            cr_q_t = s_t;

            //std::cout << "cr_q_t before: " << cr_q_t << std::endl;

		    size_t step = 0;
		    while (step < num_iterations) {
                step++;
                /****************
                 ** Local step **
                 ****************/
		        for (size_t i = 0; i < num_constraints; i++) {
		            auto proj = cr_constraints.at(i)->projectOnConstraintSet(cr_q_t);
		            //std::cout << "Proj: " << proj << std::endl;
		            projections[i] = proj;
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

                    //Check:
                    //std::cout << "Projection check for constraint " << i << " (should be 0): " << Ai * c->getSelectionMatrix() * cr_q_t - c->getBiMatrix() * projections.at(i) << std::endl;
		        }

		        //Compute right-hand side
		        cr_rhs.setZero();
		        cr_rhs = cr_M_star * s_t / (h*h);
		        for (size_t i = 0; i < num_constraints; i++) {
		            auto c = cr_constraints.at(i);
		            cr_rhs += c->getSelectionMatrixWeighted().transpose()
		                    * c->getAiMatrix().transpose() * c->getBiMatrix()
		                    * projections.at(i);
		        }

		        //std::cout << "left hand side: " << cr_lhs << std::endl;
		        //std::cout << "right hand side: " << cr_rhs << std::endl;

                cr_solver.compute(cr_lhs);
                cr_q_t = cr_solver.solve(cr_rhs);

                //std::cout << "cr_qt_solved: " << cr_q_t << std::endl;

                if (cr_solver.info() == Eigen::Success) {
                    //Update velocities and angular velocities
                    Orientations new_quat;
                    Vector new_pos;
                    separatePosQuat(&cr_q_t, new_pos, new_quat, cr_num_positions * 3);
                    cr_velocities = (new_pos - cr_positions) / h;
                    cr_positions = new_pos;
                    cr_angular_velocities = 2/h * quat2pos(conjugateQuat(cr_orientations) * new_quat);
                    cr_orientations = new_quat;
                } else {
                    std::cout << "Unable to solve rods constraints" << std::endl;
                    return false;
                }
            }


            is_simulating = false;

		}

		void setGrab(const std::vector<Index>& grabVerts, const std::vector<Eigen::Vector3f> grabPos) {
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

		Positions* getRodsPositions() {
		    upload_pos = vec2pos(cr_positions);
		    return &upload_pos;
		}

		std::vector<Index> getRodIndices() {
		    return rod_indices;
		}

		Positions* getRodsTangents() {
		    upload_tan = vec2pos(cr_positions);
		    upload_tan.resize(cr_num_positions, 3);
		    Vector3 t(0,1,0);
		    size_t j = 0;
		    for (size_t i = 0; i < cr_num_quaternions; i++) {
		        Quaternion q = cr_orientations[i];
		        upload_tan.row(j++) = q.normalized().toRotationMatrix() * t;
                if (i != 0 and std::find(rod_indices.begin(), rod_indices.end(), i) != rod_indices.end()) {
                    upload_tan.row(j++) = q.normalized().toRotationMatrix() * t;
                }
		    }
		    upload_tan.row(upload_tan.rows()-1) = upload_tan.row(upload_tan.rows()-2);

		    return &upload_tan;
		}

		Positions* getRodsNormals() {
		    upload_normals.resize(cr_num_positions, 3);
		    Vector3 n(-1, 0, 0);
		    size_t j = 0;
		    for (size_t i = 0; i < cr_num_quaternions; i++) {
		        Quaternion q = cr_orientations[i];
                upload_normals.row(j++) = q.normalized().toRotationMatrix() * n;
                if (i != 0 and std::find(rod_indices.begin(), rod_indices.end(), i) != rod_indices.end()) {
                    upload_normals.row(j++) = q.normalized().toRotationMatrix() * n;
                }
            }
		    upload_normals.row(upload_normals.rows()-1) = upload_normals.row(upload_normals.rows()-2);

		    return &upload_normals;
		}

	protected:

        // If m_hasGrab is true, vertices with the indices in m_grabVerts will be forced
		// to the positions set in m_grabPos
		bool m_hasGrab = false;
        std::vector<Index> m_grabVerts;
        std::vector<Eigen::Vector3f> m_grabPos;
        std::vector<Edge> edges;
		const double h = 0.05; //Simulation step size
		const float gravity = -1;
        /******************
         * Mesh variables *
         ******************/
        size_t m;

        Positions m_f_ext;
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
        const float cr_unit_weight = 0.01f;
        const float cr_radius = CRConstraint::radius;
        const float cr_density = 1.3f;
        std::vector<float> cr_segments_length;
        std::vector<Index> rod_indices;

        //Position variables (dim: cr_num_coord)
        size_t cr_size; //Full size of the flat variables
        Vector cr_positions_init;
        Vector cr_positions;

        Vector cr_q_t; //The flat vector that holds every variables (positions and quaternions)
        size_t cr_num_coord;
        size_t cr_num_positions = 0; //Number of positions
        Vector cr_velocities;
        Vector cr_f_ext;
        SparseMatrix cr_masses;
        SparseMatrix cr_masses_inv;

        //Angular variables (dim: cr_num_flat_pos)
        size_t cr_num_flat_pos;
        Vector cr_angular_velocities;
        Vector cr_torques;
        SparseMatrix cr_J_vec;

        SparseMatrix cr_J_vec_inv; //TODO: Use the diagonal matrix type
        SparseMatrix cr_J_quat;
        SparseMatrix cr_J_quat_inv;
        SparseMatrixRM cr_M_star;

        size_t cr_num_quaternions;
        Orientations cr_orientations;

        std::vector<CRConstraint*> cr_constraints;
        SparseSolver cr_solver;
        Vector cr_rhs;
        SparseMatrix cr_lhs;

		//State variables
		bool is_ready;
		bool is_simulating;
		bool use_cosserat_rods;
		Positions upload_pos, upload_tan, upload_normals;

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
                auto new_ground_c = new GroundConstraint(m, i, 1000.0);
                m_constraints.push_back(new_ground_c);
            }
		}

		void addEdgeSpringConstraints() {
		    //Iterate over all edges to add the spring constraint
            for (auto e : edges) {
                auto new_edge_c = new EdgeSpringConstraint(m, e, 1.0, 0.001, 1.5);
                m_constraints.push_back(new_edge_c);
            }
		}

		void addSSConstraints() {
		    for (size_t ind = 0; ind < rod_indices.size(); ind++) {
		        Index rod_index = rod_indices.at(ind);
		        Index next_index = ind == rod_indices.size() - 1 ? cr_num_positions : rod_indices.at(ind + 1);

                for (size_t i = rod_index; i < next_index - 1; i++) {
                    auto new_c = new StretchShearConstraint(cr_size, 1.0, i*3,
                                                            cr_num_positions*3 + (i - ind) * 4, cr_segments_length.at(ind));
                    cr_constraints.push_back(new_c);
                }
            }
		}

		void addFixedPosConstraint() {
		    for (auto i: rod_indices) {
                Vector3 p;
                p << cr_positions.coeff(i, 0), cr_positions.coeff(i, 1), cr_positions.coeff(i, 2);
                auto fpc = new FixedPointConstraint(cr_size, 100.0,i, p);
                cr_constraints.push_back(fpc);
            }
		}

		void addBTConstraints() {
		    for (Index ind = 0; ind < rod_indices.size(); ind++) {
		        Index rod_index = rod_indices.at(ind);
		        Index next_index = ind == rod_indices.size() - 1 ? cr_num_positions : rod_indices.at(ind + 1);
                for (size_t i = rod_index; i < next_index - 2; i++) {
                    auto new_c = new BendTwistConstraint(cr_size, 1.0, cr_num_positions*3 + (i-ind)*4, cr_segments_length.at(ind));
                    cr_constraints.push_back(new_c);
                }
            }
		}

		//Create imaginary parts of quaternions from vectors. Real part is 0
		static Orientations pos2quat(const Vector p) {
		    //The size of the vector must a multiple of 3
		    assert(p.size() % 3 == 0);
		    size_t q_size = p.size() / 3;
		    Orientations u(q_size);
		    for (size_t i = 0; i < q_size; i++) {
                Quaternion quat(0, p[i*3], p[i*3 + 1], p[i*3 + 2]);
                u[i] = quat.normalized();
		    }
		    return u;
		}

		//Creates vectors from the imaginary part of the quaternions
		static Vector quat2pos(const Orientations quat) {
		    Vector pos;
		    pos.setZero(quat.size() * 3);
		    for (size_t i = 0; i < quat.size(); i++) {
		        auto q = quat[i].normalized();
		        //The real part should be null
		        //assert(q.w() == 0.0);
		        pos[i*3] = q.x();
		        pos[i*3 + 1] = q.y();
		        pos[i*3 + 2] = q.z();
		    }
		    return pos;
		}

		static Positions vec2pos(const Vector p) {
		    assert(p.size() % 3 == 0);
		    size_t q_size = p.size() / 3;
		    Positions q;
		    q.setZero(q_size, 3);

		    for (size_t i = 0; i < q_size; i++) {
		        Vector3 vec;
		        vec << p[i*3], p[i*3 + 1], p[i*3 + 2];
		        q.row(i) = vec;
		    }

		    return q;
		}

		static Vector pos2flatVec(const Positions pos) {
		    Vector v;
		    v.setZero(pos.rows() * 3);

		    for (size_t i = 0; i < pos.rows(); i++) {
		        v[i*3] = pos.coeff(i, 0);
		        v[i*3 + 1] = pos.coeff(i, 1);
		        v[i*3 + 2] = pos.coeff(i, 2);
		    }

		    return v;
		}

		Orientations quatFromPos(const Vector pos) {
		    assert(pos.size() % 3 == 0);
		    size_t num_pos = pos.size() / 3;
		    Orientations quat;
		    quat.resize(num_pos - rod_indices.size());

		    size_t j = 0;
		    for (size_t i = 0; i < num_pos-1; i++) {
		        size_t index = i*3;
		        size_t next_index = (i + 1) * 3;
                size_t previous_index = (i - 1) * 3;
                Vector3 v_i, v_next, v_previous;
                v_i << pos[index], pos[index + 1], pos[index + 2];
                v_next << pos[next_index], pos[next_index + 1], pos[next_index + 2];
                Quaternion q;
                if (std::find(rod_indices.begin(), rod_indices.end(), i) != rod_indices.end()) { //The first position always points forward
                    v_previous = Vector3(v_next);
                    q = Quaternion::FromTwoVectors(v_next - v_i, Vector3::UnitY());
                    if (i != 0) i++;
                } else {
                    v_previous << pos[previous_index], pos[previous_index + 1], pos[previous_index + 2];
                    //q = Quaternion::FromTwoVectors(v_i - v_previous, v_next - v_i);
                    q = Quaternion::FromTwoVectors(v_i - v_previous, Vector3::UnitY());
                }
                quat[j++] = q;
            }

		    return quat;
		}

		static Orientations conjugateQuat(const Orientations quat) {
		    Orientations new_quat;
		    new_quat.resize(quat.size());
		    for (size_t i = 0; i < quat.size(); i++) {
		        new_quat[i] = quat[i].conjugate();
		    }
		    return new_quat;
		}

		//Perform cross product column-wise of two vectors
		static Vector vectorCross(const Vector a, const Vector b) {
		    Vector cross;
		    assert(a.size() == b.size() && a.size() % 3 == 0);
		    cross.setZero(a.rows());
            Vector3 a_v, b_v, cross_v;

            for (size_t i = 0; i < cross.size(); i+=3) {
		        a_v << a[i], a[i+1], a[i+2];
		        b_v << b[i], b[i+1], b[i+2];
		        cross_v = a_v.cross(b_v);
		        cross[i] = cross_v[0];
		        cross[i+1] = cross_v[1];
		        cross[i+2] = cross_v[2];
		    }
		    return cross;
		}

		static Vector posQuatConcat(const Vector pos, const Orientations u) {
		    size_t m = pos.size();
		    size_t dim = m + 4*u.rows();
		    Vector fp(dim);

		    for (size_t i = 0; i < m; i++) {
		        fp(i) = pos[i];
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

		static void separatePosQuat(const Vector* concat, Vector& pos, Orientations& quat, size_t separation) {
		    size_t size = concat->size();
		    assert((size - separation) % 4 == 0);
		    size_t num_quat = (size-separation) / 4;

		    pos.setZero(separation);

		    quat.resize(num_quat);

		    for (size_t i = 0; i < separation; i++) {
		        pos.coeffRef(i) = concat->coeff(i);
		    }

		    for (size_t i = 0; i < num_quat; i++) {
		        size_t index = separation + i*4;
		        quat.coeffRef(i) = Eigen::Quaterniond(concat->coeff(index),
                                                    concat->coeff(index+1),
                                                    concat->coeff(index+2),
                                                    concat->coeff(index+3));
		    }
		}

		static SparseMatrix stackDiagonal(const Vector x, size_t times) {
		    size_t size = x.size();
		    SparseMatrix concat;
		    concat.resize(times * size, times * size);
		    concat.setZero();
		    for (size_t i = 0; i < times; i++) {
		        for (size_t j = 0; j < size; j++) {
		            concat.coeffRef(i*size + j, i*size + j) = x.coeff(j);
		        }
		    }
		    return concat;
		}

	};

}