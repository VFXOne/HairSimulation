//
// Created by lucas on 24.09.19.
//

#ifndef APPLET_CONSTRAINT_H
#define APPLET_CONSTRAINT_H

#include "projdyn_types.h"

#include "Edge.h"
typedef ProjDyn::Positions Positions;

typedef unsigned int Scalar;
class Constraint {
public:
    virtual Positions projectOnConstraintSet() = 0;

    void initSM(Scalar rows, Scalar cols) {
        m_selectionMatrix = ProjDyn::SparseMatrix(rows, cols);
    }

protected:
    ProjDyn::SparseMatrix m_selectionMatrix;

    double clamp(const double x, const double min, const double max) {
        if (x < min) {
            return min;
        } else if (x > max) {
            return max;
        } else {
            return x;
        }
    }
};

class EdgeSpringConstraint : public Constraint {

    EdgeSpringConstraint(const Scalar m, Edge* edge, const Scalar weight,
            const Scalar rangeMin, const Scalar rangeMax)
    : Constraint() {
        initSM(1, m);
        m_rangeMax = rangeMax;
        m_rangeMin = rangeMin;
        m_weight = weight;
        m_edge = edge;
    }

    Positions projectOnConstraintSet(Positions q_n) {
        Positions edge_coord = (q_n.row(m_edge->getFirstPos()) - q_n.row(m_edge->getSecondPos()));
        double edge_length = edge_coord.norm();
        edge_coord /= edge_length; //Normalize
        double target_length = clamp(edge_length, m_rangeMin, m_rangeMax);
        return edge_coord *= target_length;
    }

protected:
    Scalar m_rangeMin;
    Scalar m_rangeMax;
    Scalar m_weight;
    Edge* m_edge;
};

#endif //APPLET_CONSTRAINT_H
