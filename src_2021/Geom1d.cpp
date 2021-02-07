/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include "Geom1d.h"

Geom1d::Geom1d() {
}

Geom1d::~Geom1d() {
}

Geom1d::Geom1d(const Geom1d &copy) {
    fNodeIndices = copy.fNodeIndices;
}

Geom1d& Geom1d::operator=(const Geom1d& copy) {
    fNodeIndices = copy.fNodeIndices;
    return *this;
}

void Geom1d::Shape(const VecDouble &xi, VecDouble &phi, MatrixDouble &dphi) {
    double qsi = xi[0];
    
    phi[0] = (1.0 - qsi) / 2.;
    phi[1] = (1.0 + qsi) / 2.;
    
    dphi(0, 0) = -0.5;
    dphi(0, 1) = 0.5;
}

void Geom1d::X(const VecDouble &xi, MatrixDouble &NodeCo, VecDouble &x) {
    int nrow = NodeCo.rows();
    for (int i = 0; i < nrow; i++) {
        x[i] = NodeCo(i, 0)*(1. - xi[0])*0.5 + NodeCo(i, 1)*(1. + xi[0])*0.5;
    }
}

void Geom1d::GradX(const VecDouble &xi, MatrixDouble &NodeCo, VecDouble &x, MatrixDouble &gradx) {
    int nrow = NodeCo.rows();
    int ncol = NodeCo.cols();

    gradx.resize(nrow, 1);
    gradx.setZero();

    VecDouble phi(2);
    MatrixDouble dphi(2, 2);
    Shape(xi, phi, dphi);
    for (int i = 0; i < ncol; i++) {
        for (int j = 0; j < nrow; j++) {
            gradx(j, 0) += NodeCo(j, i) * dphi(0, i);
        }
    }
}

void Geom1d::SetNodes(const VecInt &nodes) {
    fNodeIndices = nodes;
}

void Geom1d::GetNodes(VecInt &nodes) {
    nodes = fNodeIndices;
}

int Geom1d::NodeIndex(int node) {
    return fNodeIndices[node];
}

int Geom1d::NumNodes() {
    return nCorners;    
}

GeoElementSide Geom1d::Neighbour(int side) {
    return fNeighbours[side];
}

void Geom1d::SetNeighbour(int side, const GeoElementSide &neighbour) {
    fNeighbours[side]=neighbour;
}
