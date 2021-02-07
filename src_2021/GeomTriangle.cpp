/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "GeomTriangle.h"

GeomTriangle::GeomTriangle() {
}

GeomTriangle::~GeomTriangle() {
}

GeomTriangle::GeomTriangle(const GeomTriangle &copy) {
    fNodeIndices = copy.fNodeIndices;

}

GeomTriangle& GeomTriangle::operator=(const GeomTriangle& copy) {
    fNodeIndices = copy.fNodeIndices;

    return *this;
}

void GeomTriangle::Shape(const VecDouble& xi, VecDouble& phi, MatrixDouble& dphi) {
    double qsi = xi[0];
    double eta = xi[1];

    phi[0] = 1. - qsi - eta;
    phi[1] = qsi;
    phi[2] = eta;

    dphi(0, 0) = -1.;
    dphi(1, 0) = -1.;

    dphi(0, 1) = 1.;
    dphi(1, 1) = 0.;

    dphi(0, 2) = 0.;
    dphi(1, 2) = 1.;
}

void GeomTriangle::X(const VecDouble &xi, MatrixDouble &NodeCo, VecDouble &x) {
    VecDouble phi(3);
    MatrixDouble dphi(2, 3);

    Shape(xi, phi, dphi);
    int space = NodeCo.rows();

    for (int i = 0; i < space; i++) {
        x[i] = 0.0;
        for (int j = 0; j < 3; j++) {
            x[i] += phi[j] * NodeCo(i, j);
        }
    }
}

void GeomTriangle::GradX(const VecDouble &xi, MatrixDouble &NodeCo, VecDouble &x, MatrixDouble &gradx) {
    int space = NodeCo.rows();
    int ncol = NodeCo.cols();

    gradx.resize(space, 2);
    gradx.setZero();

    VecDouble phi(3);
    phi.setZero();
    MatrixDouble dphi(2, 3);
    Shape(xi, phi, dphi);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < space; j++) {
            gradx(j, 0) += NodeCo(j, i) * dphi(0, i);
            gradx(j, 1) += NodeCo(j, i) * dphi(1, i);
        }
    }
}

void GeomTriangle::SetNodes(const VecInt &nodes) {
    fNodeIndices = nodes;
}

void GeomTriangle::GetNodes(VecInt &nodes) {
    nodes = fNodeIndices;
}

int GeomTriangle::NodeIndex(int node) {
    return fNodeIndices[node];
}

int GeomTriangle::NumNodes() {
    return nCorners;
}

GeoElementSide GeomTriangle::Neighbour(int side) {
    return fNeighbours[side];
}

void GeomTriangle::SetNeighbour(int side, const GeoElementSide &neighbour) {
    fNeighbours[side] = neighbour;
}
