//
// Created by lucas on 24.09.19.
//

#ifndef APPLET_PDCONSTRAINT_H
#define APPLET_PDCONSTRAINT_H

#include <iostream>
#include "Edge.h"
#include "projdyn_types.h"

using namespace ProjDyn;

class PDConstraint {
public:

    PDConstraint(Index numVertices, Scalar weight) {
        m_weight = weight;
        m = numVertices;
    }

    virtual Positions projectOnConstraintSet(Positions& q) = 0;

    virtual SparseMatrix getSelectionMatrix() = 0;

    virtual SparseMatrix getSelectionMatrixWeighted() = 0;

    SparseMatrix computeLHS() {
        return m_weight*m_selectionMatrix.transpose() * m_selectionMatrix;
    }

protected:

    SparseMatrix m_selectionMatrix;
    float m_weight;
    Index m;

    void initSM(Index rows, Index cols) {
        m_selectionMatrix = SparseMatrix(rows, cols);
        m_selectionMatrix.setZero();
    }

    double clamp(const double x, const double min, const double max) {
        if (x <= min) {
            return min;
        } else if (x >= max) {
            return max;
        } else {
            return x;
        }
    }
};

class EdgeSpringConstraint : public PDConstraint {
public:

    EdgeSpringConstraint(const Index m, Edge edge, const Scalar weight,
            const float rangeMin, const float rangeMax)
    : PDConstraint(m, weight) {
        m_rangeMax = rangeMax;
        m_rangeMin = rangeMin;
        m_edge = edge;

        initSM(1, m);
        Index i = m_edge.getFirstPos();
        Index j = m_edge.getSecondPos();
        m_selectionMatrix.coeffRef(0, m_edge.getFirstPos()) = 1;
        m_selectionMatrix.coeffRef(0, m_edge.getSecondPos()) = -1;
    }

    Positions projectOnConstraintSet(Positions& q_n) override {
        Index i1 = m_edge.getFirstPos();
        Index i2 = m_edge.getSecondPos();
        Positions edge_coord = q_n.row(m_edge.getFirstPos()) - q_n.row(m_edge.getSecondPos());
        Scalar edge_length = edge_coord.norm();
        edge_coord /= edge_length; //Normalize
        Scalar target_length = clamp(edge_length, m_rangeMin, m_rangeMax);
        return edge_coord *= target_length;
    }

    SparseMatrix getSelectionMatrix() override {
        return m_selectionMatrix;
    }

    SparseMatrix getSelectionMatrixWeighted() override {
        return m_weight * m_selectionMatrix;
    }

protected:
    Scalar m_rangeMin;
    Scalar m_rangeMax;
    Edge m_edge = Edge(0, 1);
};

class GroundConstraint : public PDConstraint {
public:

    GroundConstraint(const Index m, const Index vertexIndex, const Scalar weight,
            Scalar groundHeight = -1.0f, Index floorCoord = 1)
    : PDConstraint(m, weight){
        m_groundHeight = groundHeight;
        m_constrainedVertex = vertexIndex;
        m_groundCoord = floorCoord;

        initSM(1, m);
        m_selectionMatrix.coeffRef(0, m_constrainedVertex) = 1;
    }

    Positions projectOnConstraintSet(Positions& q) override {
        Scalar coord = q(m_constrainedVertex, m_groundCoord);
        Positions targetPos = q.row(m_constrainedVertex);

        if (coord < m_groundHeight) {
            targetPos(0, m_groundCoord) = m_groundHeight;
        }

        return targetPos;
    }

    SparseMatrix getSelectionMatrix() override {
        return m_selectionMatrix;
    }

    SparseMatrix getSelectionMatrixWeighted() override {
        return m_weight * m_selectionMatrix;
    }

protected:
    Scalar m_groundHeight;
    Index m_groundCoord;
    Index m_constrainedVertex;

};

class CRConstraint: public PDConstraint {
public:

    CRConstraint(Scalar numVertices, float weight)
    : PDConstraint(numVertices, weight) {
        m_weight = weight;
        m = numVertices;
    }

    Positions projectOnConstraintSet(Positions& q) override {
        throw std::logic_error("Function not implemented for Cosserat Rods. "
                               "Please use projectOnConstraintSet(FlatPos& fp)");
    }

    virtual Vector projectOnConstraintSet(Vector& fp) = 0;

    SparseMatrix getAiMatrix() {
        return A_i;
    }

    SparseMatrix getBiMatrix() {
        return B_i;
    }

    SparseMatrix getSelectionMatrix() override {
        return m_selectionMatrix;
    }

    SparseMatrix getSelectionMatrixWeighted() override {
        return m_weight * m_selectionMatrix;
    }


protected:
    SparseMatrix A_i;
    SparseMatrix B_i;
};

class StretchShearConstraint: public CRConstraint {
public:

    StretchShearConstraint(size_t num_coord, float weight_multiplier, size_t pos_index, size_t quat_index, float segment_length)
            : CRConstraint(num_coord, weight_multiplier) {

        seg_length = segment_length;
        p_index = pos_index;
        q_index = quat_index;

        m_weight = E * M_PI * cr_radius * cr_radius * seg_length * weight_multiplier;

        A_i.resize(7, 10);
        A_i.setZero();
        const float x = 1/seg_length;
        A_i.coeffRef(0,0) = x;
        A_i.coeffRef(1,1) = x;
        A_i.coeffRef(2,2) = x;
        A_i.coeffRef(0,3) = -x;
        A_i.coeffRef(1,4) = -x;
        A_i.coeffRef(2,5) = -x;
        A_i.coeffRef(3,6) = 1;
        A_i.coeffRef(4,7) = 1;
        A_i.coeffRef(5,8) = 1;
        A_i.coeffRef(6,9) = 1;

        B_i.resize(7, 7);
        B_i.setIdentity();

        initSM(10, num_coord);
        m_selectionMatrix.coeffRef(0,p_index + 3) = 1; //Select x_n_1
        m_selectionMatrix.coeffRef(1,p_index + 4) = 1;
        m_selectionMatrix.coeffRef(2,p_index + 5) = 1;
        m_selectionMatrix.coeffRef(3,p_index) = 1; //Select x_n
        m_selectionMatrix.coeffRef(4,p_index + 1) = 1;
        m_selectionMatrix.coeffRef(5,p_index + 2) = 1;
        m_selectionMatrix.coeffRef(6,q_index) = 1; //Select u_n
        m_selectionMatrix.coeffRef(7,q_index + 1) = 1;
        m_selectionMatrix.coeffRef(8,q_index + 2) = 1;
        m_selectionMatrix.coeffRef(9,q_index + 3) = 1;
    }

    Vector projectOnConstraintSet(Vector& q) override {
        Vector3 x_n, x_n_1, x_f, d_3;
        Quaternion u_n, diff_u_n, u_n_star;

        x_n << q.coeff(p_index), q.coeff(p_index + 1), q.coeff(p_index + 2);
        x_n_1 << q.coeff(p_index + 3), q.coeff(p_index + 4), q.coeff(p_index + 5);
        x_f = (x_n_1 - x_n) / seg_length;

        u_n = Quaternion(q.coeff(q_index), q.coeff(q_index+1), q.coeff(q_index+2), q.coeff(q_index+3));
        u_n.normalize();

        d_3 = u_n.toRotationMatrix() * Vector3::UnitY();
        //d_3 = u_n_star.toRotationMatrix() * x_f;
        d_3.normalize();

        diff_u_n = Quaternion::FromTwoVectors(d_3, x_f.normalized());
        u_n_star = u_n * diff_u_n;

        std::cout << "x_f: " << x_f << std::endl;
        std::cout << "d_3: " << d_3 << std::endl;
        std::cout << "diff_u_n: " << diff_u_n.toRotationMatrix() << std::endl;
        std::cout << "u_n_star: " << u_n_star.toRotationMatrix() << std::endl;

        Vector p_i;
        p_i.resize(7);

        p_i << d_3.coeff(0), d_3.coeff(1), d_3.coeff(2),
                u_n_star.w(), u_n_star.x(),u_n_star.y(), u_n_star.z();

        return p_i;
    }

protected:
    Scalar seg_length;
    Index q_index;
    Index p_index;
};

class BendTwistConstraint: public CRConstraint {
public:

    BendTwistConstraint(Index num_coord, Scalar weight_multiplier, Index quat_index, float segment_length)
    : CRConstraint(num_coord, weight_multiplier) {
        q_index = quat_index;

        m_weight = E * M_PI * cr_radius * cr_radius * cr_radius * cr_radius / ((1+poisson) * segment_length) * weight_multiplier;

        A_i.resize(8, 8);
        A_i.setIdentity();

        B_i.resize(8, 8);
        B_i.setIdentity();

        initSM(8, num_coord);
        m_selectionMatrix.coeffRef(0, q_index) = 1; // Select u_n
        m_selectionMatrix.coeffRef(1, q_index+1) = 1;
        m_selectionMatrix.coeffRef(2, q_index+2) = 1;
        m_selectionMatrix.coeffRef(3, q_index+3) = 1;

        m_selectionMatrix.coeffRef(4, q_index+4) = 1; //Select u_(n+1)
        m_selectionMatrix.coeffRef(5, q_index+5) = 1;
        m_selectionMatrix.coeffRef(6, q_index+6) = 1;
        m_selectionMatrix.coeffRef(7, q_index+7) = 1;

    }

    Vector projectOnConstraintSet(Vector& q) override {
        Quaternion u_n(q.coeff(q_index), q.coeff(q_index+1),q.coeff(q_index+2), q.coeff(q_index+3));
        Quaternion u_n_1(q.coeff(q_index+4), q.coeff(q_index+5),q.coeff(q_index+6), q.coeff(q_index+7));
        u_n.normalize();
        u_n_1.normalize();
        Quaternion r_curvature = u_n.conjugate() * u_n_1;
        r_curvature.normalize();
        //r_curvature.w() = 0;

        Quaternion r_curvature_c(r_curvature.conjugate());
        r_curvature_c.normalize();

        auto rotMatrix = r_curvature.toRotationMatrix();
        r_curvature = Quaternion(0.5 * rotMatrix);

        rotMatrix = r_curvature_c.toRotationMatrix();
        r_curvature_c = Quaternion(0.5 * rotMatrix);

        Quaternion u_n_star = u_n * r_curvature;
        Quaternion u_n_1_star = u_n_1 * r_curvature_c;

        Vector sol;
        sol.resize(8);
        sol << u_n_star.w(), u_n_star.x(), u_n_star.y(), u_n_star.z(),
                u_n_1_star.w(), u_n_1_star.x(), u_n_1_star.y(), u_n_1_star.z();

        return sol;
    }

private:
    Index q_index;
};

class FixedPointConstraint: public CRConstraint {
public:
    FixedPointConstraint(Index num_coord, Scalar weight, Index pos_index, const Vector3 fixed_pos)
    : CRConstraint(num_coord, weight) {
        p_index = pos_index;
        f_pos = fixed_pos;

        A_i.resize(3, 3);
        A_i.setIdentity();

        B_i.resize(3, 3);
        B_i.setIdentity();

        initSM(3, num_coord);
        m_selectionMatrix.coeffRef(0, p_index) = 1;
        m_selectionMatrix.coeffRef(1, p_index+1) = 1;
        m_selectionMatrix.coeffRef(2, p_index+2) = 1;
    }

    Vector projectOnConstraintSet(Vector& q) override {
        Vector p_i;
        p_i.resize(3);
        p_i << f_pos.x(), f_pos.y(), f_pos.z();
        return p_i;
    }

protected:
    Index p_index;
    Vector3 f_pos;
};

/*
 * Constraint describing a rod point sticking to another point (on a mesh for example)
 */
//TODO: Add normal parameter
class MovingPointConstraint: public CRConstraint {
public:
    MovingPointConstraint(Index num_coord, Scalar weight, Index pos_index, Positions* positions, Index moving_pos_index)
    : CRConstraint(num_coord, weight) {
        p_index = pos_index;
        m_pos = positions;
        m_pos_index = moving_pos_index;

        A_i.resize(3, 3);
        A_i.setIdentity();

        B_i.resize(3, 3);
        B_i.setIdentity();

        initSM(3, num_coord);
        m_selectionMatrix.coeffRef(0, p_index) = 1;
        m_selectionMatrix.coeffRef(1, p_index+1) = 1;
        m_selectionMatrix.coeffRef(2, p_index+2) = 1;
    }

    Vector projectOnConstraintSet(Vector& q) override {
        Vector p_i;
        p_i.resize(3);
        p_i << m_pos->coeff(m_pos_index, 0), m_pos->coeff(m_pos_index, 1), m_pos->coeff(m_pos_index, 2);
        return p_i;
    }

protected:
    Index p_index;
    Positions* m_pos;
    Index m_pos_index;
};

#endif //APPLET_PDCONSTRAINT_H
