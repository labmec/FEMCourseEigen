/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include "GeomQuad.h"

GeomQuad::GeomQuad() {
}

GeomQuad::~GeomQuad() {
}

GeomQuad::GeomQuad(const GeomQuad &copy) {
    fNodeIndices = copy.fNodeIndices;
}

GeomQuad& GeomQuad::operator=(const GeomQuad& copy) {
    fNodeIndices = copy.fNodeIndices;
    return *this;
}

void GeomQuad::Shape(const VecDouble &xi, VecDouble &phi, MatrixDouble &dphi) {
    double qsi = xi[0];
    double eta = xi[1];

    phi[0] = 0.25 * (1. - qsi)*(1. - eta);
    phi[1] = 0.25 * (1. + qsi)*(1. - eta);
    phi[2] = 0.25 * (1. + qsi)*(1. + eta);
    phi[3] = 0.25 * (1. - qsi)*(1. + eta);

    dphi(0, 0) = 0.25 * (eta - 1.);
    dphi(1, 0) = 0.25 * (qsi - 1.);

    dphi(0, 1) = 0.25 * (1. - eta);
    dphi(1, 1) = -0.25 * (1. + qsi);

    dphi(0, 2) = 0.25 * (1. + eta);
    dphi(1, 2) = 0.25 * (1. + qsi);

    dphi(0, 3) = -0.25 * (1. + eta);
    dphi(1, 3) = 0.25 * (1. - qsi);

}

void GeomQuad::X(const VecDouble &xi, MatrixDouble &NodeCo, VecDouble &x) {
    VecDouble phi(4);
    MatrixDouble dphi(2, 4);

    Shape(xi, phi, dphi);
    int space = NodeCo.rows();

    for (int i = 0; i < space; i++) {
        x[i] = 0.0;
        for (int j = 0; j < 4; j++) {
            x[i] += phi[j] * NodeCo(i, j);
        }
    }
}

void GeomQuad::GradX(const VecDouble &xi, MatrixDouble &NodeCo, VecDouble &x, MatrixDouble &gradx) {

    int space = Dimension;
    int ncol = NodeCo.cols();

    gradx.resize(space, 2);
    gradx.setZero();

    VecDouble phi(4);
    MatrixDouble dphi(2, 4);
    Shape(xi, phi, dphi);
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < space; j++) {
            gradx(j, 0) += NodeCo(j, i) * dphi(0, i);
            gradx(j, 1) += NodeCo(j, i) * dphi(1, i);
        }
    }
}

void GeomQuad::SetNodes(const VecInt &nodes) {
    fNodeIndices = nodes;
}

void GeomQuad::GetNodes(VecInt &nodes) {
    nodes = fNodeIndices;
}

int GeomQuad::NodeIndex(int node) {
    return fNodeIndices[node];
}

int GeomQuad::NumNodes() {
    return nCorners;
}

GeoElementSide GeomQuad::Neighbour(int side) {
    return fNeighbours[side];
}

void GeomQuad::SetNeighbour(int side, const GeoElementSide &neighbour) {
    fNeighbours[side] = neighbour;
}
