

#include <iostream>
#include <math.h>
#include "DataTypes.h"
#include "IntRule1d.h"
#include "IntRuleQuad.h"
#include "IntRuleTriangle.h"
#include "tpanic.h"
#include "DataTypes.h"
#include "GeoMesh.h"
#include "GeoElement.h"
#include "GeomQuad.h"
#include "GeoElementTemplate.h"
#include "MathStatement.h"
#include "Poisson.h"
#include "L2Projection.h"
#include "Analysis.h"
#include "IntRule1d.h"
#include "Shape1d.h"


#ifdef WIN32
#define __PRETTY_FUNCTION__ __FUNCTION__
#endif

using std::cout;
using std::endl;
using std::cin;

void test1D();

void testQuad();

void testTriangle();
//twodim Test
double Poli2D(VecDouble &coord, int orderx, int ordery);
double Poli2DQuad(VecDouble &coord, int orderx, int ordery);
void UXi(VecDouble &coord, VecDouble &uXi, VecDouble &gradu);
VecDouble X(VecDouble &coordXi);
void Jacobian(VecDouble &Coord, MatrixXd &jacobian, double &detjac);
double InnerVec(VecDouble &S , VecDouble &T);

int main ()
{
    //  test1D();
    //  testTriangle();
    //  testQuad();
    
    cout<<"---------------------------------------------------------------------------------"<< endl;
    cout<<"----------------- Regra de integração - Tabela de valores :  --------------------"<< endl;
    
    double weight;
    VecDouble CoordXi(2), uXi(1), gradu(2);
    
    //Teste 1 - Regra de integração - Tabela de valores;
    
    int order =1;
    double val1 = 0.;
    int NPoint;
    
    while(fabs(val1-12) >= 1.e-5) {
        IntRuleQuad TestInt1(order);
        NPoint = TestInt1.NPoints();
        for (int i=0; i<NPoint; i++) {
            TestInt1.Point(i, CoordXi, weight);
            UXi(CoordXi,uXi,gradu);
            val1 = val1 + weight*uXi[0];
        }
        order++;
    }
    
    cout<<"Resultado Integral u(xi,eta):  "<< val1  << "  -> Número de pontos de integração = "<< NPoint << endl;
    
    //Teste 2 - Regra de integração - Tabela de valores;
    
    order =17;
    double val2 = 0.;
    
    while(fabs(val2-5.03941698) >= 1.e-5) {
        IntRuleQuad TestInt2(order);
        NPoint = TestInt2.NPoints();
        for (int i=0; i<NPoint; i++) {
            TestInt2.Point(i, CoordXi, weight);
            UXi(CoordXi,uXi,gradu);
            val2 += weight*InnerVec(gradu, gradu);
        }
        val2=sqrt(val2);
    }
    cout<<"Resultado Integral ||Gradu(xi,eta)||:  "<< val2 << "  -> Número de pontos de integração = "<< NPoint << endl;
    
    
    GeoMesh mesh;
    CoordXi.setZero();
    MatrixXd jacobian, jacinv;
    double detjac;
    
    //Teste 3 - Regra de integração - Tabela de valores;
    
    order =1;
    weight=0.;
    double valA=0;
    CoordXi.setZero(),uXi.setZero(),gradu.setZero();
    
    while(fabs(valA-80) >= 1.e-5) {
        IntRuleQuad TestInt3(order);
        
        VecDouble coord;
        
        NPoint = TestInt3.NPoints();
        VectorXi nodes(NPoint);
        
        mesh.SetNumElements(1);
        
        for (int i=0; i<NPoint; i++) {
            nodes[i]=i;
        }
        mesh.SetElement(0, new GeoElementTemplate<GeomQuad>(nodes,1,&mesh,0));
        
        for (int i=0; i<NPoint; i++) {
            TestInt3.Point(i, CoordXi, weight);
            UXi(CoordXi,uXi,gradu);
            //coord = X(CoordXi);
            Jacobian(CoordXi, jacobian, detjac);
            valA = valA + detjac*weight;
        }
    }
    cout<<"Resultado Área quad. mapeado:  "<< valA << "  -> Número de pontos de integração = "<< NPoint << endl;
    
    
    //Teste 4 - Regra de integração - Tabela de valores;
    
    order =19;
    weight=0.;
    double val3 = 0.;
    CoordXi.setZero(),uXi.setZero(),gradu.setZero();
    
    //   while(fabs(val3-239.49661609) >= 1.e-5) {
    IntRuleQuad TestInt4(order);
    NPoint = TestInt4.NPoints();
    VectorXi nodes4(NPoint);
    mesh.SetNumElements(1);
    
    for (int i=0; i<NPoint; i++) {
        nodes4[i]=i;
    }
    mesh.SetElement(0, new GeoElementTemplate<GeomQuad>(nodes4,1,&mesh,0));
    
    for (int i=0; i<NPoint; i++) {
        TestInt4.Point(i, CoordXi, weight);
        UXi(CoordXi,uXi,gradu);
        //coord = X(CoordXi);
        Jacobian(CoordXi, jacobian, detjac);
        val3 = val3 + detjac*weight*uXi[0];
    }
    //   }
    cout<<"Resultado Integral u(x,y):  "<< val3 << "  -> Número de pontos de integração = "<< NPoint << " (Obs: é preciso mais pontos para a precisão 1.e-5)"<< endl;
    
    //Teste 5 - Regra de integração - Tabela de valores;
    
    order =18;
    weight=0.;
    CoordXi.setZero(),uXi.setZero(),gradu.setZero();
    
    double val4 = 0.;
    
    while(fabs(val4-22.53695786) >= 1.e-5) {
        IntRuleQuad TestInt5(order);
        NPoint = TestInt5.NPoints();
        VecInt nodes5(NPoint);
        
        for (int i=0; i<NPoint; i++) {
            nodes5[i]=i;
        }
        mesh.SetElement(0, new GeoElementTemplate<GeomQuad>(nodes5,1,&mesh,0));
        
        for (int i=0; i<NPoint; i++) {
            TestInt5.Point(i, CoordXi, weight);
            UXi(CoordXi,uXi,gradu);
            //coord = X(CoordXi);
            Jacobian(CoordXi, jacobian, detjac);
            val4 += detjac*weight*InnerVec(gradu, gradu);
            
        }
        
        val4=sqrt(val4);
    }
    cout<<"Resultado Integral ||Gradu(x,y)||:  "<< val4 << "  -> Número de pontos de integração = "<< NPoint << endl;
    
    
    
    
    cout<<"---------------------------------------------------------------------------------"<< endl;
    cout<<"---------------- Função Gauleg retirada do Numerical Recipes : ------------------"<< endl;
    
    
    //Teste 1 - Função Gauleg retirada do Numerical Recipes;
    
    int order1=1;
    val1=0.;
    int nPoints = 1;
    
    while(fabs(val1-12.) >= 1.e-5) {
        IntRuleQuad TestIntGauss(order1);
        VecDouble x(nPoints);
        VecDouble wvec(nPoints);
        
        TestIntGauss.gaulegQuad(-1, 1, x, wvec);
        for (int i=0; i<nPoints*nPoints; i++) {
            CoordXi[0]=x[i];
            CoordXi[1]=x[i+nPoints*nPoints];
            UXi(CoordXi,uXi,gradu);
            
            val1 += wvec[i]*uXi[0];
        }
        nPoints++;
    }
    
    cout<<"Resultado Integral u(xi,eta):  "<< val1 << "  -> Número de pontos de integração = "<< nPoints*nPoints << endl;
    
    //Teste 2 - Função Gauleg retirada do Numerical Recipes;
    val2 = 0.;
    nPoints = 9;
    
    while(fabs(val2-5.03941698) >= 1.e-5) {
        IntRuleQuad TestIntGauss2(order1);
        VecDouble x2(nPoints);
        VecDouble wvec2(nPoints);
        
        TestIntGauss2.gaulegQuad(-1, 1, x2, wvec2);
        for (int i=0; i<nPoints*nPoints; i++) {
            CoordXi[0]=x2[i];
            CoordXi[1]=x2[i+nPoints*nPoints];
            UXi(CoordXi,uXi,gradu);
            val2 += wvec2[i]*InnerVec(gradu, gradu);
        }
        
        val2=sqrt(val2);
        nPoints++;
    }
    cout<<"Resultado Integral ||Gradu(xi,eta)||:  "<< val2 << "  -> Número de pontos de integração = "<< nPoints*nPoints << endl;
    
    
    //Teste 3 - Função Gauleg retirada do Numerical Recipes;
    
    valA=0.;
    nPoints=1;
    while(fabs(valA-80.) >= 1.e-5) {
        IntRuleQuad TestIntGauss3(order1);
        
        VecDouble x3(nPoints);
        VecDouble wvec3(nPoints);
        
        TestIntGauss3.gaulegQuad(-1, 1, x3, wvec3);
        for (int i=0; i<nPoints*nPoints; i++) {
            CoordXi[0]=x3[i];
            CoordXi[1]=x3[i+nPoints*nPoints];
            UXi(CoordXi,uXi,gradu);
            //coord = X(CoordXi);
            Jacobian(CoordXi, jacobian, detjac);
            
            valA += detjac*wvec3[i];
        }
        nPoints++;
    }
    
    cout<<"Resultado Área quad. mapeado:  "<< valA << "  -> Número de pontos de integração = "<< nPoints*nPoints << endl;
    
    //Teste 4 - Função Gauleg retirada do Numerical Recipes;
    val3=0.;
    nPoints = 14;
    while(fabs(val3-239.49661609) >= 1.e-5) {
        IntRuleQuad TestIntGauss4(order1);
        VecDouble x4(nPoints);
        VecDouble wvec4(nPoints);
        
        TestIntGauss4.gaulegQuad(-1, 1, x4, wvec4);
        for (int i=0; i<nPoints*nPoints; i++) {
            CoordXi[0]=x4[i];
            CoordXi[1]=x4[i+nPoints*nPoints];
            UXi(CoordXi,uXi,gradu);
            //coord = X(CoordXi);
            Jacobian(CoordXi, jacobian, detjac);
            
            val3 += detjac*wvec4[i]*uXi[0];
        }
        nPoints++;
    }
    
    cout<<"Resultado Integral u(x,y):  "<< val3 << "  -> Número de pontos de integração = "<< nPoints*nPoints << endl;
    
    //Teste 5 - Função Gauleg retirada do Numerical Recipes;
    
    
    val4=0.;
    nPoints = 10;
    while(fabs(val4-22.53695786) >= 1.e-5) {
        IntRuleQuad TestIntGauss5(order1);
        VecDouble x5(nPoints);
        VecDouble wvec5(nPoints);
        
        TestIntGauss5.gaulegQuad(-1, 1, x5, wvec5);
        for (int i=0; i<nPoints*nPoints; i++) {
            CoordXi[0]=x5[i];
            CoordXi[1]=x5[i+nPoints*nPoints];
            UXi(CoordXi,uXi,gradu);
            //coord = X(CoordXi);
            Jacobian(CoordXi, jacobian, detjac);
            
            val4 += detjac*wvec5[i]*InnerVec(gradu, gradu);
        }
        val4=sqrt(val4);
        nPoints++;
    }
    
    cout<<"Resultado Integral ||Gradu(x,y)||:  "<< val4 << "  -> Número de pontos de integração = "<< nPoints*nPoints << endl;
    cout<<"-----------------------------------"<< endl;
    
    
    
    
    return 0;
    
}
//f[x_, y_] = x^2 + y^2 - 2 x y;
double Poli2D(VecDouble &coord, int orderx, int ordery)
{
    return pow(coord[0],orderx) * pow(1.-coord[1],ordery);
}

double Poli2DQuad(VecDouble &coord, int orderx, int ordery)
{
    return pow(coord[0],orderx) * pow(coord[1],ordery);
}

void UXi(VecDouble &coord, VecDouble &uxi, VecDouble &gradu)
{
    uxi[0] = 3.+cos(3.*coord[1])*sin(4.*coord[0]);
    gradu[0] = 4.*cos(3.*coord[1])* cos(4.*coord[0]);
    gradu[1] = -3.*sin(3.*coord[1])* sin(4.*coord[0]);
}


double InnerVec(VecDouble &S , VecDouble &T){
    
    double Val = 0;
    if(S.size()!=T.size()){
        DebugStop();
    }
    for(int i = 0; i < S.size(); i++){
        Val += S[i]*T[i];
    }
    return Val;
    
}





VecDouble X(VecDouble &CoordXi)
{
    
    VecDouble coord(2);
    
    coord[0] = 5.*CoordXi[0]+0.5*sin(3.*CoordXi[1]);
    coord[1] = 4.*CoordXi[1]+0.3*cos(10.*CoordXi[0]);
    
    return coord;
    
}

void Jacobian(VecDouble &Coord, MatrixXd &jacobian, double &detjac)
{
    
    jacobian.resize(2, 2);
    jacobian.setZero();
    
    jacobian(0,0)=5.;
    jacobian(0,1)=1.5*cos(3.*Coord[1]);
    jacobian(1,0)=-3.*sin(10.*Coord[0]);
    jacobian(1,1)=4.;
    
    detjac = fabs(jacobian(0, 0)*jacobian(1, 1) - jacobian(1, 0)*jacobian(0, 1));
    
}


/*
 * Calcula os valores das funcoes de forma e suas derivadas
 * @point ponto onde calcular as funcoes de forma
 * @phi valores das funcoes de forma
 * @dphi valores das derivadas das funcoes de forma
 */

void Shape(VecDouble &point, VecDouble &phi, MatrixDouble &dphi, int order)
{
    
    if (order==1) {
        
        int Indices[2][2]={{0,3},{1,2}};
        
        VecDouble coxi(1);
        coxi[0]=point[0];
        
        VecDouble coeta(1);
        coeta[0]=point[1];
        
        VecDouble phixi(2), phieta(2);
        MatrixXd dphixi(1,2),dphieta(1,2);
        VecInt ord(1);
        ord[0] = order;
        
        phi.resize(4);
        dphi.resize(2, 4);
        
        
        for (int xi=0; xi<order+1; xi++) {
            
            Shape1d::Shape(coxi,ord,phixi,dphixi);
//            TElemento::Shape1d(order, coxi, phixi, dphixi);
            
            for (int eta=0; eta<order+1; eta++) {
                
                Shape1d::Shape(coeta,ord,phieta,dphieta);
                
                phi[Indices[xi][eta]]=phixi[xi]*phieta[eta];
                
                dphi(0,Indices[xi][eta])=dphixi(0,xi)*phieta[eta];
                dphi(1,Indices[xi][eta])=dphieta(0,eta)*phixi[xi];
                
            }
            
        }
        
    }
    
    if (order==2) {
        
        int Indices[3][3]={{0,7,3},{4,8,6},{1,5,2}};
        
        VecDouble coxi(1);
        coxi[0]=point[0];
        
        VecDouble coeta(1);
        coeta[0]=point[1];
        
        VecInt ord(1);
        ord[0] = order;
        
        VecDouble phixi(2), phieta(1);
        MatrixXd dphixi(1,3),dphieta(1,3);
        
        phi.resize(9);
        dphi.resize(2, 9);
        
        
        for (int xi=0; xi<order+1; xi++) {
            
            Shape1d::Shape(coxi,ord,phixi,dphixi);
            
            for (int eta=0; eta<order+1; eta++) {
                
                Shape1d::Shape(coeta,ord,phieta,dphieta);
                
                phi[Indices[xi][eta]]=phixi[xi]*phieta[eta];
                
                dphi(0,Indices[xi][eta])=dphixi(0,xi)*phieta[eta];
                dphi(1,Indices[xi][eta])=dphieta(0,eta)*phixi[xi];
                
            }
            
        }
        
    }
    
}











double Poli1D(double coord, int polorder)
{
    double y;
    
    y = pow(coord, polorder);
    
    return y;
}

int TestTIntRule1d()
{
    cout <<__PRETTY_FUNCTION__<<endl;
    int result = 0;
    
    cout << "Teste1: Criacao de um TIntRule1d com ordem negativa" <<endl;
    try {
        IntRule1d Integral(-1);
    } catch (std::bad_exception) {
        result = 1;
    }
    if (result == 0) {
        cout << "Teste1: Errado." <<endl;
    }
    else
        cout << "Teste1: Ok." <<endl;
    
    result = 0;
    cout << "Teste2: Criacao de um TIntRulde1d com ordem > 19" <<endl;
    try {
        IntRule1d  Integral(20);
    } catch (std::bad_exception) {
        result = 1;
    }
    if (result == 0) {
        cout << "Teste2: Errado." <<endl;
    }
    else
        cout << "Teste2: Ok." <<endl;
    
    result = 0;
    cout << "Teste3: Verificacao se a ordem dada atribui corretamente o numero de pesos" <<endl;
    IntRule1d  Integral(7);
    IntRule1d  Integral2(18);
    
    if (Integral.NPoints() == 4) {
        result++;
    }
    else
        result = 0;
    
    if (Integral2.NPoints() == 10) {
        result++;
    }
    else
        result = 0;
    
    if (result == 2) {
        cout << "Teste3: Ok." <<endl;
    }
    else
        cout << "Teste3: Errado." <<endl;
    
    
    result = 0;
    cout << "Teste4: Verificacao dos pesos e coordenadas para uma ordem x." <<endl;
    VecDouble testcoord(1);
    double peso;
    int pontoteste = 3;
    
    Integral.Point(pontoteste,testcoord,peso);
    
    if (testcoord[0] == 0.33998104358485626 && peso == 0.65214515486254609) {
        result++;
    }
    else
        result = 0;
    
    pontoteste = 8;
    Integral2.Point(pontoteste,testcoord,peso);
    
    if (testcoord[0] == -0.14887433898163122 && peso == 0.29552422471475287) {
        result++;
    }
    else
        result = 0;
    
    if (result != 2) {
        cout << "Teste1: Errado." <<endl;
    }
    if (result == 2) {
        cout << "Teste4: Ok." <<endl;
        result = 1;
    }
    
    cout << "---------------//----------------" <<endl;
    
    return result;
}

int TestPointFunctionFull()
{
    cout <<__PRETTY_FUNCTION__<<endl;
    int result = 0;
    
    cout << "Teste1: Teste completo de pesos e coordenadas para uma ordem x." <<endl;
    IntRule1d  Teste(6);
    VecDouble weights(4);
    VecDouble valcoord(4);
    
    valcoord[0] = -0.86113631159405257;  weights[0] = 0.34785484513745385;
    valcoord[1] = 0.86113631159405257;   weights[1] = 0.34785484513745385;
    valcoord[2] = -0.33998104358485626;  weights[2] = 0.65214515486254609;
    valcoord[3] = 0.33998104358485626;   weights[3] = 0.65214515486254609;
    
    VecDouble coord(1);
    coord[0] = 0.;
    double peso;
    
    for(int i=0; i<Teste.NPoints() ;i++)
    {
        Teste.Point(i, coord, peso);
        
        if (coord[0] == valcoord[i] && peso == weights [i]) {
            result++;
        }
    }
    
    
    
    if (result != 4) {
        cout << "Teste1: Errado." <<endl;
    }
    if (result == 4) {
        cout << "Teste1: Ok." <<endl;
        result = 1;
    }
    
    cout << "---------------//----------------" <<endl;
    return result;
}

int TestPoliIntegration()
{
    cout <<__PRETTY_FUNCTION__<<endl;
    int result = 1;
    
    cout << "Teste1: Integral de um polinomio e verifica resultado." <<endl;
    
    for (int order = 0; order <20; order++)
    {
        IntRule1d TestInt(order);
        double val = 0;
        double weight;
        VecDouble Coord(1);
        for (int polorder = 0; polorder <= order; polorder++)
        {
            val = 0.;
            for (int i=0; i<TestInt.NPoints(); i++) {
                TestInt.Point(i, Coord, weight);
                val = val + weight*Poli1D(Coord[0], polorder);
            }
            double exact = 1./(polorder+1.)*(1-pow(-1, polorder+1));
            double error = val-exact;
            
            if (fabs(error) >= 1.e-15) {
                std::cout << "A regra de integracao de ordem " << order << " nao integrou um polinome de ordem " << polorder << " exato " << exact << " calculado " << val << " erro " << error << std::endl;
                result = 0;
            }
        }
        
    }
    if (result == 1) {
        cout << "Teste1: Ok." <<endl;
        result = 1;
    }
    else
        cout << "Teste1: Errado." <<endl;
    
    cout << "---------------//----------------" <<endl;
    return result;
}

int TestPoliIntegrationQuad()
{
    cout <<__PRETTY_FUNCTION__<<endl;
    int result = 1;
    
    cout << "Teste1: Integral de um polinomio e verifica resultado." <<endl;
    
    for (int order = 0; order <20; order++)
    {
        IntRuleQuad TestInt(order);
        
        //TestInt.Print(std::cout);
        
        double val = 0;
        double weight;
        VecDouble Coord(2);
        for (int polorderx = 0; polorderx <= order; polorderx++)
        {
            for (int polordery=0; polordery <= order; polordery++)
            {
                val = 0.;
                int NPoint = TestInt.NPoints();
                for (int i=0; i<NPoint; i++) {
                    TestInt.Point(i, Coord, weight);
                    val = val + weight*Poli2DQuad(Coord, polorderx,polordery);
                }
                double exact = 1./(polorderx+1.)*(1-pow(-1, polorderx+1))*1./(polordery+1.)*(1-pow(-1, polordery+1));
                double error = val-exact;
                
                if (fabs(error) >= 1.e-14) {
                    std::cout << "A regra de integracao quadrilatero de ordem " << order << " nao integrou um polinome de ordem " << polorderx << " , " << polordery << " exato " << exact << " calculado " << val << " erro " << error << std::endl;
                    result = 0;
                }
            }
        }
    }
    if (result == 1) {
        cout << "TesteQuad1: Ok." <<endl;
        result = 1;
    }
    else
        cout << "TesteQuad1: Errado." <<endl;
    
    cout << "---------------//----------------" <<endl;
    return result;
}

int TestPoliIntegrationTriangle()
{
    cout <<__PRETTY_FUNCTION__<<endl;
    int result = 1;
    
    cout << "Teste1: Integral de um polinomio e verifica resultado." <<endl;
    
    for (int order = 0; order <20; order++)
    {
        IntRuleTriangle TestInt(order);
        double val = 0;
        double weight;
        VecDouble Coord(2);
        for (int polorderx = 0; polorderx <= order; polorderx++)
        {
            for (int polordery=0; polordery <= order-polorderx; polordery++)
            {
                val = 0.;
                int npoints = TestInt.NPoints();
                for (int i=0; i<npoints; i++) {
                    TestInt.Point(i, Coord, weight);
                    val = val + weight*Poli2D(Coord, polorderx,polordery);
                }
                double exact = 1./(1.+polorderx)/(2.+polorderx+polordery);
                double error = val-exact;
                
                if (fabs(error) >= 1.e-15) {
                    std::cout << "A regra de integracao triangular de ordem " << order << " nao integrou um polinome de ordem " << polorderx << " , " << polordery << " exato " << exact << " calculado " << val << " erro " << error << std::endl;
                    result = 0;
                }
            }
        }
    }
    if (result == 1) {
        cout << "TesteTriangle1: Ok." <<endl;
        result = 1;
    }
    else
        cout << "TesteTriangle1: Errado." <<endl;
    
    cout << "---------------//----------------" <<endl;
    return result;
}


void test1D()
{
    int res = 0;
    int count = 0;
    
    res = TestTIntRule1d();
    if (res == 1) {
        count++;
    }
    
    res = TestPointFunctionFull();
    if (res == 1) {
        count++;
    }
    
    res = TestPoliIntegration();
    if (res == 1) {
        count++;
    }
    
    
    if (count == 1) {
        cout <<"\nSeu codigo passou no teste de Integracao para funcoes de uma dimensao. (:D)" << "\n---------------//----------------";
    }
    
}

void testQuad()
{
    int result = 0;
    result = TestPoliIntegrationQuad();
    if (result) {
        cout <<"\nSeu codigo passou no teste de Integracao sobre quadrilatero. (:D)" << "\n---------------//----------------";
        
    }
}

void testTriangle()
{
    int result = 0;
    result = TestPoliIntegrationTriangle();
    if (result) {
        cout <<"\nSeu codigo passou no teste de Integracao sobre triangulo. (:D)" << "\n---------------//----------------";
        
    }
}
