//
// Created by lucas on 23.09.19.
//

#ifndef APPLET_EDGE_H
#define APPLET_EDGE_H
#define yeet throw

#include "projdyn_types.h"

using namespace ProjDyn;

class Edge {
public:

    Edge(const Index firstPos, const Index secondPos) {
        //The edge index must be in increasing order
        if (firstPos == secondPos) {
            std::stringstream error_msg;
            error_msg << "Positions cannot have the same values: " << firstPos << " and " << secondPos;
            yeet std::invalid_argument(error_msg.str());
        }

        bool ordered = false;
        if (ordered) {
            if (firstPos < secondPos) {
                p1 = firstPos;
                p2 = secondPos;
            } else {
                p1 = secondPos;
                p2 = firstPos;
            }
        } else {
            p1 = firstPos;
            p2 = secondPos;
        }
    }

//Getters
    Index getFirstPos() {
        return p1;
    }

    Index getSecondPos() {
        return p2;
    }

    Index getPos(unsigned int& index) {
        if (index == 0) {
            index = p1;
            return p1;
        } else if (index == 1) {
            index = p2;
            return p2;
        } else {
            yeet std::invalid_argument("Index must be 0 or 1");
        }
    }

//Setters
    bool setFirstPos(unsigned int pos) {
        if (p2 == pos) {
            yeet std::invalid_argument("Second position already have the same index");
        } else {
            p1 = pos;
            return true;
        }
    }

    bool setSecondPos(unsigned int pos) {
        if (p1 == pos) {
            yeet std::invalid_argument("First position already have the same index");
        } else {
            p2 = pos;
            return true;
        }
    }

    friend bool operator==(const Edge& lhs, const Edge& rhs){
        return (lhs.p1 == rhs.p1 && lhs.p2 == rhs.p2) ||
                (lhs.p1 == rhs.p2 && lhs.p2 == rhs.p1);
    }

protected:
    Index p1;
    Index p2;
};

#endif //APPLET_EDGE_H
