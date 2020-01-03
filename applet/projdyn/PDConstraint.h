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

    virtual Vector projectOnConstraintSet(const Vector& fp) = 0;

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

    SparseMatrix getLHS() {
        return lhs;
    }


protected:

    void computeLHS()  {
        lhs = m_weight * m_selectionMatrix.transpose() * A_i.transpose() * B_i;
    }

    SparseMatrix A_i;
    SparseMatrix B_i;
    SparseMatrix lhs;
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

        computeLHS();
    }

    Vector projectOnConstraintSet(const Vector& q) override {
        assert(q.size() > q_index);
        Vector3 x_n, x_n_1, x_f, d_3;
        Quaternion u_n, diff_u_n, u_n_star;

        x_n << q.coeff(p_index), q.coeff(p_index + 1), q.coeff(p_index + 2);
        x_n_1 << q.coeff(p_index + 3), q.coeff(p_index + 4), q.coeff(p_index + 5);
        x_f = (x_n_1 - x_n) / seg_length;

        u_n = Quaternion(q.coeff(q_index), q.coeff(q_index+1), q.coeff(q_index+2), q.coeff(q_index+3));

        //TODO: Add constant default tangent variable
        d_3 = u_n.normalized().toRotationMatrix() * Vector3::UnitY();

        diff_u_n = Quaternion::FromTwoVectors(x_f, d_3);
        u_n_star = u_n * diff_u_n;

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

        computeLHS();
    }

    Vector projectOnConstraintSet(const Vector& q) override {
        assert(q.size() > q_index);
        Quaternion u_n(q.coeff(q_index), q.coeff(q_index+1),q.coeff(q_index+2), q.coeff(q_index+3));
        Quaternion u_n_1(q.coeff(q_index+4), q.coeff(q_index+5),q.coeff(q_index+6), q.coeff(q_index+7));
        Quaternion r_curvature = u_n.conjugate() * u_n_1;

        //Divide the quaternion by 2
        auto rotMatrix = r_curvature.normalized().toRotationMatrix();
        r_curvature = Quaternion(rotMatrix / 2);

        Quaternion u_n_star = u_n * r_curvature;
        Quaternion u_n_1_star = u_n_1 * r_curvature.conjugate();

        Vector sol;
        sol.resize(8);

        sol <<  u_n_star.w(), u_n_star.x(), u_n_star.y(), u_n_star.z(),
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

        computeLHS();
    }

    Vector projectOnConstraintSet(const Vector& q) override {
        assert(q.size() > p_index);
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

        computeLHS();
    }

    Vector projectOnConstraintSet(const Vector& q) override {
        assert(q.size() > p_index);
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

class NormalConstraint: public CRConstraint {
public:
    NormalConstraint(Index num_coord, Scalar weight, Index pos_index, Vector3 normal)
    : CRConstraint(num_coord, weight) {
        p_index = pos_index;
        m_normal = normal;

        A_i.resize(3, 3);
        A_i.setIdentity();

        B_i.resize(3, 3);
        B_i.setIdentity();

        initSM(3, num_coord);
//        m_selectionMatrix.coeffRef(0, p_index) = 1;
//        m_selectionMatrix.coeffRef(1, p_index+1) = 1;
//        m_selectionMatrix.coeffRef(2, p_index+2) = 1;
        m_selectionMatrix.coeffRef(0, p_index+3) = 1;
        m_selectionMatrix.coeffRef(1, p_index+4) = 1;
        m_selectionMatrix.coeffRef(2, p_index+5) = 1;

        computeLHS();
    }

    Vector projectOnConstraintSet(const Vector& q) override {
        assert(q.size() > p_index);
        Vector p_i;
        p_i.resize(3);
        Vector3 pos(q.coeff(p_index), q.coeff(p_index+1), q.coeff(p_index+2));
        Vector3 next(q.coeff(p_index+3), q.coeff(p_index+4), q.coeff(p_index+5));
        double dist = (pos - next).norm();
        Vector3 proj = pos + m_normal*dist;
        p_i << next.x(), next.y(), next.z();
        return p_i;
    }

protected:
    Index p_index;
    Vector3 m_normal;
};

class NormalConstraintQuat: public CRConstraint {
public:
    NormalConstraintQuat(Index num_coord, Scalar weight, Index quat_index, Quaternion normal)
    : CRConstraint(num_coord, weight) {
        q_index = quat_index;
        m_normal = normal;

        A_i.resize(4, 4);
        A_i.setIdentity();

        B_i.resize(4, 4);
        B_i.setIdentity();

        initSM(4, num_coord);
        m_selectionMatrix.coeffRef(0, q_index) = 1;
        m_selectionMatrix.coeffRef(1, q_index+1) = 1;
        m_selectionMatrix.coeffRef(2, q_index+2) = 1;
        m_selectionMatrix.coeffRef(3, q_index+3) = 1;

        computeLHS();
    }

    Vector projectOnConstraintSet(const Vector& q) override {
        assert(q.size() > q_index);
        Vector p_i;
        p_i.resize(4);
        p_i << m_normal.w(), m_normal.x(), m_normal.y(), m_normal.z();
        return p_i;
    }

protected:
    Index q_index;
    Quaternion m_normal;
};

class SphereConstraint : public CRConstraint {
public:
    SphereConstraint(Index num_coord, Scalar weight, Scalar radius, Index pos_index, Vector3 center_pos, Positions* positions, Index moving_pos_index, Scalar forceFactor = 1.)
            :
            CRConstraint(num_coord, weight)
    {
        m_center_pos = center_pos;
        m_radius = radius;
        m_pos_index = pos_index;
        m_force_factor = forceFactor;
        m_positions = positions;
        m_mov_pos_index = moving_pos_index;
        m_ref_point.x() = m_positions->coeff(moving_pos_index);
        m_ref_point.y() = m_positions->coeff(moving_pos_index+1);
        m_ref_point.z() = m_positions->coeff(moving_pos_index+2);

        A_i.resize(3, 3);
        A_i.setIdentity();

        B_i.resize(3, 3);
        B_i.setIdentity();

        initSM(3, num_coord);
        m_selectionMatrix.coeffRef(0, pos_index) = 1;
        m_selectionMatrix.coeffRef(1, pos_index+1) = 1;
        m_selectionMatrix.coeffRef(2, pos_index+2) = 1;

        computeLHS();
    }

    Vector projectOnConstraintSet(const Vector& q) override {
        // Check for correct size of the projection auxiliary variable;
        assert(q.size() > m_pos_index);
        //TODO: update sphere center
        Vector3 new_ref, diff, new_center;
        new_ref = m_positions->row(m_mov_pos_index);
        diff = new_ref - m_ref_point;
        new_center = m_center_pos + diff;

        std::cout << "new ref = " << new_ref << std::endl;
        std::cout << "new center = " << new_center << std::endl;
        std::cout << "diff between center and ref = " << (new_center - new_ref).norm() << std::endl;

        Vector3 pos;
        Vector p_i;
        p_i.resize(3);

        pos.x() = q(m_pos_index);
        pos.y() = q(m_pos_index + 1);
        pos.z() = q(m_pos_index + 2);
        Vector3 v = pos - new_center;
        float distance = v.norm();
        //TODO: compute distance with the rod radius too

        if (distance < m_radius) {
            v = v/v.norm();
            v = v*m_radius;

            //TODO: Add forceFactor!
            p_i.x() = (1 + m_force_factor) * new_center(0) - m_force_factor * v(0);
            p_i.y() = (1 + m_force_factor) * new_center(1) - m_force_factor * v(1);
            p_i.z() = (1 + m_force_factor) * new_center(2) - m_force_factor * v(2);
            std::cout << "INSIDE!" << std::endl;

        } else {
            p_i.x() = pos.x();
            p_i.y() = pos.y();
            p_i.z() = pos.z();
        }

        return p_i;
    }

private:
    Index m_pos_index;
    Scalar m_force_factor;
    Vector3 m_center_pos;
    Scalar m_radius;
    Positions* m_positions;
    Index m_mov_pos_index;
    Vector3 m_ref_point;
};

#endif //APPLET_PDCONSTRAINT_H
