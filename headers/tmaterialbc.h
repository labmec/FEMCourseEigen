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
#ifndef TMATERIALBC_H
#define TMATERIALBC_H

#include <tmaterial.h>

#include "TVec.h"

/**
This class implements the contribution of a boundary condition

@author Philippe R. B. Devloo
*/
class TMaterialBC : public TMaterial
{
public:
    TMaterialBC(int id, int bctype, double contrstiff, double contrrhs);

    ~TMaterialBC();
  /**
   * Calcula o valor da contribui��o da equa��o variacional no ponto dado
   * na matriz de rigidez do elemento e no vetor de carga
   * @param pt [in]: ponto de integra��o de Gauss
   * @param weight [in]: peso de integra��o
   * @param phiVal [in] : valor da fun��o teste no ponto dado
   * @param dphi [in] : valor das derivadas da fun��o de forma no ponto de integra��o
   * @param elementK [inout]: matriz de rigidez do elemento
   * @param elementF [inout]: vetor de carga do elemento
   */
    virtual void Contribute (double  weight,
                             VectorXd & philVal,
                             MatrixXd & dphi,MatrixXd & elementK,
                             MatrixXd & elementF) const;
   /**
    * Calcula a contribuicao para o erro da solucao
    * @param x [in] localizacao do ponto
    * @param weight [in] peso do ponto de integracao
    * @param sol [in] valor da solucao
    * @param deriv [in] valor da derivada da solucao
    * @param function [in] ponteiro para funcao que calcula o valor exato
    * @param energy [in/out] contribuicao para norma da energia
    * @param l2 [in/out] contribuicao para norma em L2
    *
    */
   virtual void ContributeErrorSquare(VectorXd &x, double weight, double sol, VectorXd &deriv,
   	void (*function)(VectorXd& x, double &val, VectorXd&der), double &energy, double &l2)
	{
	}

    virtual void Print(std::ostream& out) const;

protected:
  /**
  * value to be contributed to the stiffness matrix
  */
    double fContrStiff;
    /**
    * value to be contributed to the rhs
    */
    double fContrRhs;
    /**
    * boundary condition type 0:Dirichlet 1:Neumann 2: Mixed
    */
    int fBCType;
};

#endif
