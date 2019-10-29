//
// Created by lucas on 24.09.19.
//

#ifndef APPLET_PDCONSTRAINT_H
#define APPLET_PDCONSTRAINT_H

#include "projdyn_types.h"

#include "Edge.h"
typedef ProjDyn::Positions Positions;
typedef ProjDyn::Orientations Orientations;
typedef ProjDyn::Vector Vector;

typedef unsigned int Scalar;
class PDConstraint {
public:

    PDConstraint(Scalar numVertices, float weight) {
        m_weight = weight;
        m = numVertices;
    }

    virtual Positions projectOnConstraintSet(Positions& q) = 0;

    virtual ProjDyn::SparseMatrix getSelectionMatrix() = 0;

    virtual ProjDyn::SparseMatrix getSelectionMatrixWeighted() = 0;

    ProjDyn::SparseMatrix computeLHS() {
        return m_weight*m_selectionMatrix.transpose() * m_selectionMatrix;
    }

protected:

    ProjDyn::SparseMatrix m_selectionMatrix;
    float m_weight;
    Scalar m;

    void initSM(Scalar rows, Scalar cols) {
        m_selectionMatrix = ProjDyn::SparseMatrix(rows, cols);
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

    EdgeSpringConstraint(const Scalar m, Edge edge, const float weight,
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
        int i1 = m_edge.getFirstPos();
        int i2 = m_edge.getSecondPos();
        Positions edge_coord = q_n.row(m_edge.getFirstPos()) - q_n.row(m_edge.getSecondPos());
        double edge_length = edge_coord.norm();
        edge_coord /= edge_length; //Normalize
        double target_length = clamp(edge_length, m_rangeMin, m_rangeMax);
        return edge_coord *= target_length;
    }

    ProjDyn::SparseMatrix getSelectionMatrix() override {
        return m_selectionMatrix;
    }

    ProjDyn::SparseMatrix getSelectionMatrixWeighted() override {
        return m_weight * m_selectionMatrix;
    }

protected:
    float m_rangeMin;
    float m_rangeMax;
    Edge m_edge = Edge(0, 1);
};

class GroundConstraint : public PDConstraint {
public:

    GroundConstraint(const Scalar m, const Scalar vertexIndex, const float weight,
            float groundHeight = -1.0f, unsigned int floorCoord = 2)
    : PDConstraint(m, weight){
        m_groundHeight = groundHeight;
        m_constrainedVertex = vertexIndex;
        m_groundCoord = floorCoord;

        initSM(1, m);
        m_selectionMatrix.coeffRef(0, m_constrainedVertex) = 1;
    }

    Positions projectOnConstraintSet(Positions& q) override {
        float coord = q(m_constrainedVertex, m_groundCoord);
        Positions targetPos = q.row(m_constrainedVertex);

        if (coord < m_groundHeight) {
            targetPos(0, m_groundCoord) = m_groundHeight;
        }

        return targetPos;
    }

    ProjDyn::SparseMatrix getSelectionMatrix() override {
        return m_selectionMatrix;
    }

    ProjDyn::SparseMatrix getSelectionMatrixWeighted() override {
        return m_weight * m_selectionMatrix;
    }

protected:
    float m_groundHeight;
    unsigned int m_groundCoord;
    Scalar m_constrainedVertex;

};

class CRConstraint: public PDConstraint {
public:

    CRConstraint(Scalar m, Scalar weight)
    : PDConstraint(m, weight) {

    }

    Positions projectOnConstraintSet(Positions& q) override {
        throw std::logic_error("Function not implemented for Cosserat Rods. "
                               "Please use projectOnConstraintSet(ProjDyn::FlatPos& fp)");
    }

    virtual Vector projectOnConstraintSet(Vector& fp) = 0;

    virtual ProjDyn::SparseMatrix getAiMatrix() = 0;

    virtual ProjDyn::SparseMatrix getBiMatrix() = 0;

protected:
    ProjDyn::SparseMatrix A_i;
    ProjDyn::SparseMatrix B_i;
};

class StretchShearConstraint: public CRConstraint {
public:

    StretchShearConstraint(size_t num_coord, float weight, size_t pos_index, size_t quat_index, Scalar segment_length)
            : CRConstraint(num_coord,weight) {

        seg_length = segment_length;
        p_index = pos_index;
        q_index = quat_index;

        A_i.resize(7, 10);
        A_i.setZero();
        Scalar x = 1/seg_length;
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
        B_i.setZero();
        B_i.coeffRef(0,0) = 1;
        B_i.coeffRef(1,1) = 1;
        B_i.coeffRef(2,2) = 1;
        B_i.coeffRef(4, 4) = 1;
        B_i.coeffRef(5, 5) = 1;
        B_i.coeffRef(6, 6) = 1;

        initSM(10, num_coord);
        m_selectionMatrix.coeffRef(0,p_index + 3) = 1; //To get x_n_1
        m_selectionMatrix.coeffRef(1,p_index + 4) = 1;
        m_selectionMatrix.coeffRef(2,p_index + 5) = 1;
        m_selectionMatrix.coeffRef(3,p_index) = 1;
        m_selectionMatrix.coeffRef(4,p_index + 1) = 1;
        m_selectionMatrix.coeffRef(5,p_index + 2) = 1;
        m_selectionMatrix.coeffRef(6,q_index) = 1;
        m_selectionMatrix.coeffRef(7,q_index + 1) = 1;
        m_selectionMatrix.coeffRef(8,q_index + 2) = 1;
        m_selectionMatrix.coeffRef(9,q_index + 3) = 1;
    }

    Vector projectOnConstraintSet(Vector& q) override {
        ProjDyn::Vector3 x_n, x_n_1, x_f, d_3;
        ProjDyn::Quaternion u_n, diff_u_n, u_n_star;

        Vector p_i;
        p_i.resize(7);

        x_n << q.coeff(p_index), q.coeff(p_index + 1), q.coeff(p_index + 2);
        x_n_1 << q.coeff(p_index + 3), q.coeff(p_index + 4), q.coeff(p_index + 5);
        x_f = (x_n_1 - x_n) / seg_length;

        u_n = ProjDyn::Quaternion(q.coeff(q_index), q.coeff(q_index+1), q.coeff(q_index+2), q.coeff(q_index+3));
        d_3 = u_n.toRotationMatrix() * ProjDyn::Vector3(0,0,1);
        diff_u_n = ProjDyn::Quaternion::FromTwoVectors(d_3, x_f);

        u_n_star = u_n * diff_u_n;

        p_i << d_3.coeff(0), d_3.coeff(1), d_3.coeff(2),
                u_n_star.w(), u_n_star.x(),u_n_star.y(), u_n_star.z();

        return p_i;
    }

    ProjDyn::SparseMatrix getAiMatrix() override {
        return A_i;
    }

    ProjDyn::SparseMatrix getBiMatrix() override {
        return B_i;
    }

    ProjDyn::SparseMatrix getSelectionMatrix() override {
        return m_selectionMatrix;
    }

    ProjDyn::SparseMatrix getSelectionMatrixWeighted() override {
        return m_weight * m_selectionMatrix;
    }

protected:
    Scalar seg_length;
    size_t q_index;
    size_t p_index;
};

class BendingTwistingConstraint: public CRConstraint {
public:

    BendingTwistingConstraint(size_t num_coord, float weight)
    : CRConstraint(num_coord, weight) {

    }

    Vector projectOnConstraintSet(Vector& q) override {
        return q;
    }

};

class FixedPointConstraint: public CRConstraint {
public:
    FixedPointConstraint(size_t num_coord, float weight, size_t pos_index, ProjDyn::Vector3 fixed_pos)
    : CRConstraint(num_coord, weight) {
        p_index = pos_index;
        f_pos = fixed_pos;

        A_i.resize(3, 3);
        A_i.setIdentity();

        B_i.resize(3, 3);
        B_i.setIdentity();

        initSM(3, num_coord);
        m_selectionMatrix.coeffRef(0, p_index);
        m_selectionMatrix.coeffRef(1, p_index+1);
        m_selectionMatrix.coeffRef(2, p_index+2);
    }

    Vector projectOnConstraintSet(Vector& q) override {
        Vector p_i;
        p_i.resize(3);
        p_i << f_pos.x(), f_pos.y(), f_pos.z();
        return p_i;
    }

    ProjDyn::SparseMatrix getAiMatrix() override {
        return A_i;
    }

    ProjDyn::SparseMatrix getBiMatrix() override {
        return B_i;
    }

    ProjDyn::SparseMatrix getSelectionMatrix() override {
        return m_selectionMatrix;
    }

    ProjDyn::SparseMatrix getSelectionMatrixWeighted() override {
        return m_weight * m_selectionMatrix;
    }

protected:
    size_t p_index;
    ProjDyn::Vector3 f_pos;
};

#endif //APPLET_PDCONSTRAINT_H
