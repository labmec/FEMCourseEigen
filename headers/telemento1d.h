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
#ifndef TELEMENTO1D_H
#define TELEMENTO1D_H

#include <telemento.h>
#include "TVec.h"

/**
@author Philippe R. B. Devloo
*/
class TElemento1d : public TElemento
{
public:
    TElemento1d(int matid, int order, VectorXi & nodes);
    
    TElemento1d();

    ~TElemento1d();

    virtual MElementType getElType();
    virtual void CalcStiff(TMalha &malha, MatrixXd& stiff, MatrixXd& rhs);
    virtual void Jacobian(VectorXd &point, MatrixXd& jacobian, MatrixXd& jacinv, double &detjac, TMalha& malha);
    
    
    /*
     * Calcula os valores das funcoes de forma e suas derivadas
     * @point ponto onde calcular as funcoes de forma
     * @phi valores das funcoes de forma
     * @dphi valores das derivadas das funcoes de forma
     */
    virtual void Shape(VectorXd &point, VectorXd &phi, MatrixXd &dphi);
    
    /**
     * Calcula o erro do elemento
     * @param exact funcao que calcula o valor exato
     * @param energy [out] erro na norma de energia
     * @param l2 [out] erro na norma l2
     */
virtual void Error(MatrixXd &solution, TMalha &malha, void (*f)(VectorXd &,double &, VectorXd &), double &energy, double &l2);
    
    /**
     * Calcula o a solucao e o gradiente do elemento
     * @param solution vetor de coeficientes multiplicadores (dof)
     * @param malha espaco de aproximacao
     * @param uhe combinacao linear de alpha_{i} phi_{i}
     * @param duhedx combinacao linear de alpha_{i} dphi_{i}
     */
virtual void uhe(MatrixXd &solution, TMalha &malha, MatrixXd &uhe, MatrixXd &duhedx);


};

#endif
