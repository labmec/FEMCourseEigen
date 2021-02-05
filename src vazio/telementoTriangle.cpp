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


#include "telementoTriangle.h"

#include "tno.h"
#include "tmalha.h"
#include "tmaterial.h"
#include "DataTypes.h"
#include "TIntRuleTriangle.h"
#include "TIntRule1d.h"
#include "tpanic.h"
#include <math.h>
//#include "pzquad.h"

TElementoTriangle::TElementoTriangle(int matid, int order, VectorXi &nodes) : TElemento(matid, order, nodes)
{
}


TElementoTriangle::TElementoTriangle() : TElemento()
{
}


TElementoTriangle::~TElementoTriangle()
{
}


MElementType TElementoTriangle::getElType()
{
    return ETriangle;
}

void TElementoTriangle::CalcStiff(TMalha &malha, MatrixXd &stiff, MatrixXd &rhs)
{
    
    
    //identificar o objeto material
    
    int matId =this->fMaterialId;
    TMaterial *mat=malha.getMaterial(matId);
    
    
    int nshape = fNodes.size();
    
    //inicializar as matrizes
    
    stiff.resize(nshape, nshape);
    rhs.resize(nshape, 1);
    
    //Criar regra de integracao
    
    TIntRuleTriangle iRuleTriangle(fPorder+fPorder);
    int np=iRuleTriangle.NPoints();
    
    
    MatrixXd dphi(1,nshape);
    dphi.setZero();
    VectorXd phi(nshape);
    
    for (int i=0; i<np; i++) {
        double weight=0.0;
        VectorXd co(1);
        iRuleTriangle.Point(i,co,weight);
        Shape(co,phi,dphi);
        
        MatrixXd jacobian(2,2);
        MatrixXd jacinv(2,2);
        double detjac=0.0;
        
        Jacobian(co,jacobian,jacinv,detjac,malha);
        
        weight = weight*detjac;
		dphi = jacinv*dphi;
        
        
        mat->Contribute(weight,phi,dphi,stiff,rhs);
       // stiff.Print("Pos Contribute");
    }
    
    
    
}


/**
 * Calculo do jacobiano
 * @point : ponto em qual calculamos o jacobiano do mapeamento
 * @jacobian : matriz jacobiana (a dimensao depende da dimensao do elemento
 * @jacinv : inverso da matriz jacobiana
 * @malha : objeto malha necessaria para relacionar os indices dos nos com os nos reais
 */

void TElementoTriangle::Jacobian(VectorXd &point, MatrixXd &jacobian, MatrixXd &jacinv, double &detjac, TMalha &malha)
{
    
    TVec<TNo> nos(fNodes.size());
    for (int ino=0; ino<nos.Size(); ino++) {
        nos[ino]=malha.getNode(fNodes[ino]);
    }
    
    
    
    MatrixXd Coord(2,fNodes.size());
    
    for (int i=0; i<fNodes.size(); i++) {
        Coord(0,i)=nos[i].Co(0);
        Coord(1,i)=nos[i].Co(1);
    }
    
    VectorXd phixi(fNodes.size());
    MatrixXd dphixi(2,fNodes.size());
    
    jacobian.resize(2, 2);
    jacobian.setZero();
    
    jacinv.resize(2, 2);
    jacinv.setZero();
    
    this->Shape(point, phixi, dphixi);
    
    for (int xi=0; xi<fNodes.size(); xi++) {
        jacobian(0,0)+=Coord(0,xi)*dphixi(0,xi);
        jacobian(0,1)+=Coord(0,xi)*dphixi(1,xi);
        jacobian(1,0)+=Coord(1,xi)*dphixi(0,xi);
        jacobian(1,1)+=Coord(1,xi)*dphixi(1,xi);
        
    }

	jacinv.resize(2, 2);
	jacinv.setZero();
	jacinv(0, 0) = 1.;
	jacinv(1, 1) = 1.;
    FullPivLU<MatrixXd> jacobianLU(jacobian);
    jacinv = jacobianLU.matrixLU();

	detjac = fabs(jacobian(0, 0)*jacobian(1, 1) - jacobian(1, 0)*jacobian(0, 1));


    
}

/**
 * Calcula os valores das funcoes de forma e suas derivadas
 * @point ponto onde calcular as funcoes de forma
 * @phi valores das funcoes de forma
 * @dphi valores das derivadas das funcoes de forma */


void TElementoTriangle::Shape(VectorXd &point, VectorXd &phi, MatrixXd &dphi)
{
    
	
    phi.resize(fNodes.size());
    
    dphi.resize(2, fNodes.size());
    
    
	if (fPorder==1) {
        
       
        
        phi[0]=1-point[0]-point[1];
        phi[1]=point[0];
        phi[2]=point[1];
        
    
        dphi(0,0)=-1;
        dphi(0,1)=1;
        dphi(0,2)=0;
        dphi(1,0)=-1;
        dphi(1,1)=0;
        dphi(1,2)=1;
        
        
    }
    
    if (fPorder==2) {
        
        VectorXd eps(fNodes.size());
        
        eps[0]=1-point[0]-point[1];
        eps[1]=point[0];
        eps[2]=point[1];
        
        
        
        phi[0]=2*eps[0]*(eps[0]-0.5);
        phi[1]=2*eps[1]*(eps[1]-0.5);
        phi[2]=2*eps[2]*(eps[2]-0.5);
        phi[3]=4*eps[0]*eps[1];
        phi[4]=4*eps[1]*eps[2];
        phi[5]=4*eps[2]*eps[0];
        
        dphi(0,0)=-2*(0.5 - point[1] - point[0]) - 2*(1 - point[1] - point[0]);
        dphi(0,1)=2*(-0.5+point[0])+2*point[0];
        dphi(0,2)=0;
        dphi(0,3)=4*(1 - point[1] - point[0]) - 4*point[0];
        dphi(0,4)=4*point[1];
        dphi(0,5)=-4*point[1];
        
        
        dphi(1,0)=-2*(0.5 - point[1] - point[0]) - 2*(1 - point[1] - point[0]);
        dphi(1,1)=0;
        dphi(1,2)=2*(-0.5 + point[1]) + 2*point[1];
        dphi(1,3)=-4*point[0];
        dphi(1,4)=4*point[0];
        dphi(1,5)=-4*point[1] + 4*(1 - point[1] - point[0]);
        
        
    }
    
    if (fPorder==3) {
        
        VectorXd eps(fNodes.size());
        
        eps[0]=1-point[0]-point[1];
        eps[1]=point[0];
        eps[2]=point[1];
        

        
        for (int i=0; i<=2; i++) {
            phi[i]=(9/2)*eps[i]*(eps[i]-2/3)*(eps[i]-1/3);
        }
        
        phi[3]=(27/2)*eps[0]*eps[1]*(eps[0]-1/3);
        phi[4]=(27/2)*eps[0]*eps[1]*(eps[1]-1/3);
        phi[5]=(27/2)*eps[1]*eps[2]*(eps[1]-1/3);
        
        phi[6]=(27/2)*eps[1]*eps[2]*(eps[1]-2/3);
        phi[7]=(27/2)*eps[0]*eps[2]*(eps[2]-1/3);
        phi[8]=(27/2)*eps[0]*eps[2]*(eps[0]-1/3);
        
        phi[9]=27*eps[0]*eps[1]*eps[2];
        
        
        dphi(0,0)=-(9/2)*(1/3 - point[0] - point[1])*(2/3 - point[0] - point[1]) -
        (9/2)*(1/3 - point[0] - point[1])*(1 - point[0] - point[1]) -
        (9/2)*(2/3 - point[0] - point[1])*(1 - point[0] - point[1]);
		dphi(0, 1) = (9 / 2)*(-(2 / 3) + point[0])*(-(1 / 3) + point[0]) +
			(9 / 2)*(-(2 / 3) + point[0])*point[0] + (9 / 2)*(-(1 / 3) + point[0])*point[0];
		dphi(0, 2) = 0;
		dphi(0, 3) = -(27 / 2)*point[0] * (2 / 3 - point[0] - point[1]) -
			(27 / 2)*point[0] * (1 - point[0] - point[1]) +
			(27 / 2)*(2 / 3 - point[0] - point[1])*(1 - point[0] - point[1]);
		dphi(0, 4) = -(27 / 2)*(-(1 / 3) + point[0])*point[0] +
			(27 / 2)*(-(1 / 3) + point[0])*(1 - point[0] - point[1]) +
			(27 / 2)*point[0] * (1 - point[0] - point[1]);
		dphi(0, 5) = (27 / 2)*(-(1 / 3) + point[0])*point[1] + (27 / 2)*point[0] * point[1];
		dphi(0, 6) = (27 / 2)*(-(2 / 3) + point[0])*point[1] + (27 / 2)*point[0] * point[1];
		dphi(0, 7) = -(27 / 2)*(-(1 / 3) + point[1])*point[1];
		dphi(0, 8) = -(27 / 2)*(2 / 3 - point[0] - point[1])*point[1] -
			(27 / 2)*(1 - point[0] - point[1])*point[1];
		dphi(0, 9) = -27 * point[0] * point[1] + 27 * (1 - point[0] - point[1])*point[1];
        
		dphi(1, 0) = -(9 / 2)*(1 / 3 - point[0] - point[1])*(2 / 3 - point[0] - point[1]) -
			(9 / 2)*(1 / 3 - point[0] - point[1])*(1 - point[0] - point[1]) -
			(9 / 2)*(2 / 3 - point[0] - point[1])*(1 - point[0] - point[1]);
		dphi(1, 1) = 0;
		dphi(1, 2) = (9 / 2)*(-(2 / 3) + point[1])*(-(1 / 3) + point[1]) +
			(9 / 2)*(-(2 / 3) + point[1])*point[1] + (9 / 2)*(-(1 / 3) + point[1])*point[1];
		dphi(1, 3) = -(27 / 2)*point[0] * (2 / 3 - point[0] - point[1]) -
			(27 / 2)*point[0] * (1 - point[0] - point[1]);
		dphi(1, 4) = -(27 / 2)*(-(1 / 3) + point[0])*point[0];
		dphi(1, 5) = (27 / 2)*(-(1 / 3) + point[0])*point[0];
		dphi(1, 6) = (27 / 2)*(-(2 / 3) + point[0])*point[0];
		dphi(1, 7) = (27 / 2)*(1 - point[0] - point[1])*(-(1 / 3) + point[1]) +
			(27 / 2)*(1 - point[0] - point[1])*point[1] -
			(27 / 2)*(-(1 / 3) + point[1])*point[1];
		dphi(1, 8) = (27 / 2)*(2 / 3 - point[0] - point[1])*(1 - point[0] - point[1]) -
			(27 / 2)*(2 / 3 - point[0] - point[1])*point[1] -
			(27 / 2)*(1 - point[0] - point[1])*point[1];
		dphi(1, 9) = 27 * point[0] * (1 - point[0] - point[1]) - 27 * point[0] * point[1];
        
        
    }

    
}

/**
 * Calcula o erro do elemento
 * @param exact funcao que calcula o valor exato
 * @param energy [out] erro na norma de energia
 * @param l2 [out] erro na norma l2
 */
void TElementoTriangle::Error(MatrixXd &solution, TMalha &malha, void (exact)(VectorXd &, double &, VectorXd &), double &energy, double &l2)
{

	int matId =	this->fMaterialId;
	TMaterial *mat = malha.getMaterial(matId);
    
	//calcluar as coordenadas , limites dos elementos
    
	TVec<TNo> nos(fNodes.size());
	for (int ino = 0; ino<nos.Size(); ino++) {
		nos[ino] = malha.getNode(fNodes[ino]);
	}
    
	int nnodes = fNodes.size();
    
	MatrixXd Coord(2, fNodes.size());
    
	for (int i = 0; i<fNodes.size(); i++) {
		Coord(0, i) = nos[i].Co(0);
		Coord(1, i) = nos[i].Co(1);
	}
    
    
	//calcular o n˙mero de pontos de integracao
	TIntRuleTriangle iRuleTriangle(19);
	int np = iRuleTriangle.NPoints();
	energy = 0.;
	l2 = 0.;
    
	for (int ip = 0; ip<np; ip++) {
		//calcualndo o psi-ponto de integracao
		double weight = 0.0;
		VectorXd co(2);
        co[0]=0.;
        co[1]=0.;
        
        
		iRuleTriangle.Point(ip, co, weight);
		//calculando o vetor point das cordenadas em funcao de psi (volta)
        
		VectorXd xPoint(2);
        xPoint[0]=0.;
        xPoint[1]=0.;
        
        VectorXd phixi(fNodes.size());
        MatrixXd dphixi(2,fNodes.size());
        dphixi.setZero();
        
        this->Shape(co, phixi, dphixi);
        
        for (int xi=0; xi<fNodes.size(); xi++) {
            xPoint[0]+=Coord(0,xi)*phixi[xi];
            xPoint[1]+=Coord(1,xi)*phixi[xi];
            
        }
        
        
		//calculo do jacobiano
        
		MatrixXd jacobian(2, 2);
        jacobian.setZero();
		MatrixXd jacinv(2, 2);
        jacinv.setZero();
		double detjac = 0.0;
        
		Jacobian(co, jacobian, jacinv, detjac, malha);
        
		//calculo das funcoes de forma
        
		MatrixXd dphi(2, fNodes.size());
        dphi.setZero();
        
		VectorXd phi(fNodes.size());
		Shape(co, phi, dphi);
        
		weight = weight*fabs(detjac);
        
		dphi = jacinv*dphi;
        
        
		//calculo dos vetor de solucao e alphas
        
		double uh = 0.;
		VectorXd duh(2);
		duh[0] = 0.;
        duh[1] = 0.;
        
		for (int i = 0; i<nnodes; i++) {
			uh += solution(fNodes[i], 0)*phi[i];
			duh[0] += solution(fNodes[i], 0)*dphi(0, i);
			duh[1] += solution(fNodes[i], 0)*dphi(1, i);
		}
        
		mat->ContributeErrorSquare(xPoint, weight, uh, duh,
                                   exact, energy, l2);
        
	}
    
	//    energy=sqrt(energy);
	//    l2=sqrt(l2);
    
}

void TElementoTriangle::uhe(MatrixXd &solution, TMalha &malha, MatrixXd &uhe, MatrixXd &duhedx){
	uhe.resize(0, 0);
	duhedx.resize(0, 0);
    
}