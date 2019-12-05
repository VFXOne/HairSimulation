#pragma once

#define _USE_MATH_DEFINES

#include <math.h>
#include "projdyn_types.h"
#include "Edge.h"

#include "PDConstraint.h"

#include<Eigen/Geometry>
#include <memory>

#include <iostream>
#include <limits>

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

		void setMesh(Positions& pos) {
			// Pass the positions of the geometry
			m = pos.rows();
            m_positions = Positions(pos);
		}

		void setRods(const std::vector<Positions> rods) {
		    assert(rods.size() >= 1); //At least one rod
            use_cosserat_rods = true;
            rod_indices.clear();
            cr_num_positions = 0;
            cr_positions_init = Vector();

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
            cr_positions = cr_positions_init;

		    cr_velocities.setZero(cr_num_coord);
		    cr_angular_velocities.setZero(cr_num_flat_pos);
		}

		void addDefaultConstraints() {
		    cr_constraints.clear();
            addFixedPosConstraint();
            addSSConstraints();
            addBTConstraints();
		}

		void addMeshConstraints() {
		    cr_constraints.clear();
		    addSSConstraints();
		    addBTConstraints();
		    addMovingPosConstraints();
		}

		void resetPositions() {
			// Reset vertex positions to initial state and set m_velocities to 0
			if (use_cosserat_rods) {
                cr_positions = cr_positions_init;
                cr_velocities.setZero(3 * cr_num_positions);
                cr_angular_velocities.setZero(cr_num_flat_pos);
                cr_orientations = quatFromPos(cr_positions);
            }
        }

		bool isInitialized() {
			return is_ready;
		}

		bool isUsingCR() {
		    return use_cosserat_rods;
		}

        void initializeLHS() {
            //Compute left-hand side of system (can be done earlier)
            cr_lhs.resize(cr_size, cr_size);
            cr_lhs.setZero();
            cr_lhs = cr_M_star / (h*h);

            for (size_t i = 0; i < cr_constraints.size(); i++) {
                auto c = cr_constraints.at(i);
                auto Ai = c->getAiMatrix();
                cr_lhs += c->getSelectionMatrixWeighted().transpose()
                          * Ai.transpose() * Ai * c->getSelectionMatrix();
            }

            cr_solver.compute(cr_lhs);
		}

        bool initializeSystem(bool default_constraints) {
			// Setup simulation (collect constraints, set up and factorize global step system)
			if (is_simulating) {
			    return false;
			}

            if (use_cosserat_rods) {
                /****************
                 * CR Variables *
                 ****************/

                if (default_constraints) {
                    addDefaultConstraints();
                } else {
                    addMeshConstraints();
                }

                cr_f_ext.setZero(cr_num_coord);
                //Set gravity on z axis m times
                for (size_t i = 0; i < cr_num_coord; i++) {
                    cr_f_ext.coeffRef(i) = 0;
                    cr_f_ext.coeffRef(++i) = gravity;
                    cr_f_ext.coeffRef(++i) = 0;
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
                for (size_t i = 0; i < cr_masses.rows(); i++) {
                    cr_M_star.row(i) = cr_masses.row(i);
                }
                for (size_t i = 0; i < cr_J_quat.rows(); i++) {
                    for (size_t j = 0; j < cr_J_quat.cols(); j++) {
                        cr_M_star.coeffRef(i + cr_num_coord, j + cr_num_coord) = cr_J_quat.coeff(i, j);
                    }
                }

                //Finally, precompute the left-hand side matrix
                initializeLHS();
            }

			is_ready = true;
			return is_ready;
		}

        bool step(int num_iterations) {
		    if (!isInitialized()) return

		    is_simulating = true;

		    const size_t num_constraints = cr_constraints.size();
		    std::vector<Vector> projections(num_constraints);

		    Vector s_x = cr_positions + h*cr_velocities + h*h * cr_masses_inv * cr_f_ext;
            Vector cross = vectorCross(cr_angular_velocities, cr_J_vec * cr_angular_velocities);

		    Vector s_w = cr_angular_velocities + h*cr_J_vec_inv * (cr_torques - cross);
            Orientations s_w_u = pos2quat(s_w);
		    Orientations s_u = pos2quat(quat2pos(cr_orientations) + h/2 * quat2pos(cr_orientations * s_w_u));

		    Vector s_t = posQuatConcat(s_x, s_u);
            cr_q_t = s_t;

		    size_t step = 0;
		    while (step < num_iterations) {
                step++;
                /****************
                 ** Local step **
                 ****************/
#pragma omp parallel for
		        for (size_t i = 0; i < num_constraints; i++) {
		            auto proj = cr_constraints.at(i)->projectOnConstraintSet(cr_q_t);
		            projections[i] = proj;
		        }

                /*****************
                 ** Global step **
                 *****************/
		        //Compute right-hand side
		        cr_rhs.setZero();
		        cr_rhs = cr_M_star * s_t / (h*h);
		        for (size_t i = 0; i < num_constraints; i++) {
		            auto c = cr_constraints.at(i);
		            cr_rhs += c->getSelectionMatrixWeighted().transpose() //TODO: new function in constraints that computes this product in advance
		                    * c->getAiMatrix().transpose() * c->getBiMatrix()
		                    * projections.at(i);
		        }

                cr_q_t = cr_solver.solve(cr_rhs);

                if (cr_solver.info() == Eigen::Success) {
                    //Update velocities and angular velocities
                    Orientations new_quat;
                    Vector new_pos;
                    separatePosQuat(&cr_q_t, new_pos, new_quat, cr_num_coord);
					for (Index i = 0; i < new_quat.rows(); i++) {
						double phi = new_quat(i,0).angularDistance(s_u(i,0));
						if (std::abs(phi) > 1e-10) std::cout << phi << std::endl;
					}
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
		    return true;
		}

		void setGrab(const std::vector<Index>& grabVerts, const std::vector<Eigen::Vector3f> grabPos) {
			// This method will be called from the viewer if an alt-click has been performed and
			// vertices were grabbed.
			// You can use the list of grabbed vertices and the new positions to force vertices to move
			// in your step() method.

			m_grabVerts = grabVerts;
			m_grabPos = grabPos;
			m_hasGrab = true;

			assert(grabVerts.size() >= 1);
			Vector3 p = m_positions.row(grabVerts.at(0));
			Eigen::Vector3f d = grabPos.at(0);
			Vector3 f(d.x(), d.y(), d.z());
			Vector3 diff = f - p;
#pragma omp parallel for
			for (size_t i = 0; i < m_positions.rows(); i++) {
			    m_positions.row(i) += diff;
			}
		}

		void releaseGrab() {
			m_hasGrab = false;
		}

		void setTimestep(const float timestep) {
		    h = timestep;
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
		    Vector3 t = Vector3::UnitY();
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
		    Vector3 n = Vector3::UnitX();
		    size_t j = 0;
		    for (size_t i = 0; i < cr_num_quaternions; i++) {
		        Quaternion q = cr_orientations[i];
                upload_normals.row(j++) = q.normalized().toRotationMatrix() * n;
                std::cout << "uploaded normal: " << q.normalized().toRotationMatrix() * n << std::endl;
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
        double h = 0.05; //Simulation step size

        /******************
         * Mesh variables *
         ******************/
        size_t m;
        Positions m_positions;

        /***************************
         * Cosserat Rods variables *
         ***************************/
        std::vector<float> cr_segments_length;
        std::vector<Index> rod_indices;

        //Position variables (dim: cr_num_coord)
        size_t cr_size; //Full size of the flat variables
        Vector cr_positions_init;
        Vector cr_positions;

        Vector cr_q_t; //The flat vector that holds every variables (positions and quaternions)
        size_t cr_num_coord;
        size_t cr_num_positions; //Number of positions
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

		void addSSConstraints(double weight = 1.0) {
		    for (size_t ind = 0; ind < rod_indices.size(); ind++) {
		        Index rod_index = rod_indices.at(ind);
		        Index next_index = ind == rod_indices.size() - 1 ? cr_num_positions : rod_indices.at(ind + 1);

                for (size_t i = rod_index; i < next_index - 1; i++) {
                    auto new_c = new StretchShearConstraint(cr_size, weight, i*3,
                                                            cr_num_coord + (i - ind) * 4, cr_segments_length.at(ind));
                    cr_constraints.push_back(new_c);
                }
            }
		}

		void addFixedPosConstraint() {
		    for (auto i: rod_indices) {
                Vector3 p;
                p << cr_positions.coeff(i*3), cr_positions.coeff(i*3+1), cr_positions.coeff(i*3+2);
                auto fpc = new FixedPointConstraint(cr_size, 10000,i*3, p);
                cr_constraints.push_back(fpc);
            }
		}

		void addBTConstraints(double weight = 1.0) {
		    for (Index ind = 0; ind < rod_indices.size(); ind++) {
		        Index rod_index = rod_indices.at(ind);
		        Index next_index = ind == rod_indices.size() - 1 ? cr_num_positions : rod_indices.at(ind + 1);
                for (size_t i = rod_index; i < next_index - 2; i++) {
                    auto new_c = new BendTwistConstraint(cr_size, 1.0, cr_num_coord + (i-ind)*4, cr_segments_length.at(ind));
                    cr_constraints.push_back(new_c);
                }
            }
		}

		void addMovingPosConstraints() {
		    if (m_positions.size() == 0) return;
            for (auto i: rod_indices) {
                //Find closest point to the rod on the mesh
                size_t index  = 0;
                float min_distance = std::numeric_limits<float>::max();
                for (size_t j = 0; j < m_positions.rows(); j++) {
                    Vector3 p = m_positions.row(j);
                    Vector3 r(cr_positions.coeff(j), cr_positions.coeff(j+1), cr_positions.coeff(j+2));
                    float dist = (p-r).norm();
                    if (dist < min_distance) {
                        index = j;
                        min_distance = dist;
                    }
                }

                auto mpc = new MovingPointConstraint(cr_size, 1000, i, &m_positions, index);
                cr_constraints.push_back(mpc);
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

		    //for (size_t i = 0; i < separation; i++) {
		    //    pos.coeffRef(i) = concat->coeff(i);
		    //}

		    pos = concat->head(separation);

		    for (size_t i = 0; i < num_quat; i++) {
		        size_t index = separation + i*4;
		        quat.coeffRef(i) = Quaternion(concat->coeff(index),
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