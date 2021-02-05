/***************************************************************************
 *   Copyright (C) 2005 by Philippe R. B. Devloo                           *
 *   phil@corona                                                           *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include "telemento0d.h"
#include "tmaterial.h"
#include "tmalha.h"
#include "DataTypes.h"
#include "tpanic.h"

TElemento0d::TElemento0d()
{
}


TElemento0d::TElemento0d(int matid, int order, VectorXi & nodes): TElemento(matid, order, nodes)
{
}


TElemento0d::~TElemento0d()
{
}


std::string TElemento0d::TypeName(MElementType type)
{
    return TElemento::TypeName(type);
}

void TElemento0d::CalcStiff(TMalha &malha, MatrixXd &stiff, MatrixXd &rhs)
{
    //identificar o objeto material
    
    int matId =this->fMaterialId;
    TMaterial *mat=malha.getMaterial(matId);
    
    //inicializar as matrizes
    
    stiff.resize(1, 1);
    rhs.resize(1, 1);
    
    //Criar regra de integracao
    
    MatrixXd dphi(1,1);
    VectorXd phi(1), co(1);
    dphi.setZero();
    phi.setZero();
    co.setZero();
    double weight=1.0;
    Shape(co,phi,dphi);
    mat->Contribute(weight,phi,dphi,stiff,rhs);
        
    
    

}

void TElemento0d::Jacobian(VectorXd &point, MatrixXd &jacobian, MatrixXd &jacinv, double &detjac, TMalha &malha)
{
    jacobian.resize(0,0);
    jacinv.resize(0,0);
    detjac = 1;
}




//TElemento0d::getElType()

MElementType TElemento0d::getElType()
{
    return EPoint;
}

/*
 * Calcula os valores das funcoes de forma e suas derivadas
 * @point ponto onde calcular as funcoes de forma
 * @phi valores das funcoes de forma
 * @dphi valores das derivadas das funcoes de forma
 */

void TElemento0d::Shape(VectorXd &point, VectorXd &phi, MatrixXd &dphi)
{
    phi.resize(1);
    phi[0] = 1.;
    dphi.resize(1,1);
    dphi(0,0)=0;
}

/**
 * Calcula o a solucao e o gradiente do elemento
 * @param solution vetor de coeficientes multiplicadores (dof)
 * @param malha espaco de aproximacao
 * @param uhe combinacao linear de alpha_{i} phi_{i}
 * @param duhedx combinacao linear de alpha_{i} dphi_{i}
 */
void TElemento0d::uhe(MatrixXd &solution, TMalha &malha, MatrixXd &uhe, MatrixXd &duhedx){
    uhe.resize(0,0);
    duhedx.resize(0,0);

}