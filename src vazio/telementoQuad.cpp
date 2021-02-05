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

#include <math.h>
#include "telementoQuad.h"
#include "tmaterial.h"
#include "tmalha.h"
#include "DataTypes.h"
#include "tpanic.h"
#include "TIntRuleQuad.h"

TElementoQuad::TElementoQuad()
{
}


TElementoQuad::TElementoQuad(int matid, int order, VectorXi & nodes): TElemento(matid, order, nodes)
{
}


TElementoQuad::~TElementoQuad()
{
}

MElementType TElementoQuad::getElType()
{
    return EQuadrilateral;
}

void TElementoQuad::CalcStiff(TMalha &malha, MatrixXd &stiff, MatrixXd &rhs)
{

	//identificar o objeto material

	int matId = this->fMaterialId;
	TMaterial *mat = malha.getMaterial(matId);

	//inicializar as matrizes
	int nshape = fNodes.size();

	stiff.resize(nshape, nshape);
	rhs.resize(nshape, 2);

	//Criar regra de integracao

	TIntRuleQuad iRuleQuad(fPorder + fPorder);
	int np = iRuleQuad.NPoints();

	MatrixXd dphi(2, nshape);
	dphi.setZero();
	VectorXd phi(nshape);
    for (int i=0; i<nshape; i++) {
        phi[i]=0.0;
    }
    


	for (int i = 0; i<np; i++) {
		double weight = 0.0;
		VectorXd co(2);
        co[0]=0.0;
        co[1]=0.0;
        
		iRuleQuad.Point(i, co, weight);
		Shape(co, phi, dphi);

		MatrixXd jacobian(2, 2);
		jacobian.setZero();
		MatrixXd jacinv(2, 2);
		jacinv.setZero();
		double detjac = 0.0;

		Jacobian(co, jacobian, jacinv, detjac, malha);
		 
		weight = weight*detjac;
		dphi = jacinv*dphi;
		


		mat->Contribute(weight, phi, dphi, stiff, rhs);
		// stiff.Print("Pos Contribute");
	}
        
    
    

}

void TElementoQuad::Jacobian(VectorXd &point, MatrixXd &jacobian, MatrixXd &jacinv, double &detjac, TMalha &malha)
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


/*
 * Calcula os valores das funcoes de forma e suas derivadas
 * @point ponto onde calcular as funcoes de forma
 * @phi valores das funcoes de forma
 * @dphi valores das derivadas das funcoes de forma
 */

void TElementoQuad::Shape(VectorXd &point, VectorXd &phi, MatrixXd &dphi)
{
	
	if (fPorder==1) {
        
        int Indices[2][2]={{0,3},{1,2}};
        
        VectorXd coxi(1);
        coxi[0]=point[0];
        
        VectorXd coeta(1);
        coeta[0]=point[1];
        
        VectorXd phixi(fNodes.size()), phieta(fNodes.size());
        MatrixXd dphixi(1,fNodes.size()),dphieta(1,fNodes.size());
        
        phi.resize(4);
        dphi.resize(2, 4);
    
        
        for (int xi=0; xi<fPorder+1; xi++) {
            
            TElemento::Shape1d(fPorder, coxi, phixi, dphixi);
            
            for (int eta=0; eta<fPorder+1; eta++) {
                
                TElemento::Shape1d(fPorder, coeta, phieta, dphieta);
                
                phi[Indices[xi][eta]]=phixi[xi]*phieta[eta];
                
                dphi(0,Indices[xi][eta])=dphixi(0,xi)*phieta[eta];
                dphi(1,Indices[xi][eta])=dphieta(0,eta)*phixi[xi];
                
                }
            
        }

        
        
    }
    
    if (fPorder==2) {
        
        int Indices[3][3]={{0,7,3},{4,8,6},{1,5,2}};
        
        VectorXd coxi(1);
        coxi[0]=point[0];
        
        VectorXd coeta(1);
        coeta[0]=point[1];
        
        VectorXd phixi(fNodes.size()), phieta(fNodes.size());
        MatrixXd dphixi(1,fNodes.size()),dphieta(1,fNodes.size());
        
        phi.resize(9);
        dphi.resize(2, 9);
        
        
        for (int xi=0; xi<fPorder+1; xi++) {
            
            TElemento::Shape1d(fPorder, coxi, phixi, dphixi);
            
            for (int eta=0; eta<fPorder+1; eta++) {
                
                TElemento::Shape1d(fPorder, coeta, phieta, dphieta);
                
                phi[Indices[xi][eta]]=phixi[xi]*phieta[eta];
                
                dphi(0,Indices[xi][eta])=dphixi(0,xi)*phieta[eta];
                dphi(1,Indices[xi][eta])=dphieta(0,eta)*phixi[xi];
                
            }
            
        }
        
        
        
    }
    
    
    


}

void TElementoQuad::Error(MatrixXd &solution, TMalha &malha, void (exact)(VectorXd &, double &, VectorXd &), double &energy, double &l2)
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


	//calcular o n�mero de pontos de integracao
	TIntRuleQuad iRuleQuad(19);
	int np = iRuleQuad.NPoints();
	energy = 0.;
	l2 = 0.;

	for (int ip = 0; ip<np; ip++) {
		//calcualndo o psi-ponto de integracao
		double weight = 0.0;
		VectorXd co(2);

		iRuleQuad.Point(ip, co, weight);
		//calculando o vetor point das cordenadas em funcao de psi (volta)

		VectorXd xPoint(2);
        xPoint[0]=0.;
        xPoint[1]=0.;
        
        VectorXd phixi(fNodes.size());
        MatrixXd dphixi(2,fNodes.size());
        
        this->Shape(co, phixi, dphixi);
        
        for (int xi=0; xi<fNodes.size(); xi++) {
            xPoint[0]+=Coord(0,xi)*phixi[xi];
            xPoint[1]+=Coord(1,xi)*phixi[xi];
            
        }
        

		//calculo do jacobiano

		MatrixXd jacobian(2, 2);
		MatrixXd jacinv(2, 2);
		double detjac = 0.0;

		Jacobian(co, jacobian, jacinv, detjac, malha);

		//calculo das funcoes de forma

		MatrixXd dphi(2, fNodes.size());
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

	    //energy=sqrt(energy);
	    //l2=sqrt(l2);

}

/**
 * Calcula o a solucao e o gradiente do elemento
 * @param solution vetor de coeficientes multiplicadores (dof)
 * @param malha espaco de aproximacao
 * @param uhe combinacao linear de alpha_{i} phi_{i}
 * @param duhedx combinacao linear de alpha_{i} dphi_{i}
 */
void TElementoQuad::uhe(MatrixXd &solution, TMalha &malha, MatrixXd &uhe, MatrixXd &duhedx){
	uhe.resize(0, 0);
	duhedx.resize(0, 0);

}