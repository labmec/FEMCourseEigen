

#include <iostream>
#include <math.h>
#include <Eigen/Core>
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

int main ()
{
    VectorXd vec1;

    // Testing Eigen matrix
    // Copied from: https://gitlab.com/libeigen/eigen/-/blob/master/doc/examples/QuickStart_example.cpp 
    MatrixXd m(2,2); // matrix composed by doubles.
    m(0,0) = 3;
    m(1,0) = 2.5;
    m(0,1) = -1;
    m(1,1) = m(1,0) + m(0,1);
    cout << m << endl;

    VectorXd v(2);
    v(0) = 4;
    v(1) = v(0) - 1;
    cout << "Here is the vector v:\n" << v << endl;
  
    return 0;
}