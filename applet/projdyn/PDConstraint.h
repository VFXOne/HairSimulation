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

    void initSM(Scalar rows, Scalar cols) {
        m_selectionMatrix = ProjDyn::SparseMatrix(rows, cols);
        m_selectionMatrix.setZero();
    }

    ProjDyn::SparseMatrix computeLHS() {
        return m_weight*m_selectionMatrix.transpose() * m_selectionMatrix;
    }

protected:
    ProjDyn::SparseMatrix m_selectionMatrix;
    float m_weight;
    Scalar m;

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
};

class BendingTwistingConstraint: public CRConstraint {
public:

    BendingTwistingConstraint(Scalar numVertices, float weight)
    : CRConstraint(numVertices, weight) {

    }

    Vector projectOnConstraintSet(Vector& q) override {
        return q;
    }

};

class StretchShearConstraint: public CRConstraint {
public:

    StretchShearConstraint(Scalar numVertices, float weight)
    : CRConstraint(numVertices,weight) {

    }

    Vector projectOnConstraintSet(Vector& q) override {
        return q;
    }
};

class FixedPointConstraint: public CRConstraint {
public:
    FixedPointConstraint(Scalar numVertices, float weight)
    : CRConstraint(numVertices, weight) {

    }

    Vector projectOnConstraintSet(Vector& q) override {
        return q;
    }
};

#endif //APPLET_PDCONSTRAINT_H
