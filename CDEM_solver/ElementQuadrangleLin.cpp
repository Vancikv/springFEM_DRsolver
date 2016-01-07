#include "stdafx.h"
#include "ElementQuadrangleLin.h"


ElementQuadrangleLin::ElementQuadrangleLin()
{
}


ElementQuadrangleLin::~ElementQuadrangleLin()
{
}



// Calculate and store local matrices
void ElementQuadrangleLin::set_matrices()
{
}

// Calculate local stiffness matrix. Use reduced integration with hourglass stabilization.
void ElementQuadrangleLin::set_K_isoparametric()
{
}


// Calculate stress vector for element gauss points.
Eigen::MatrixXd ElementQuadrangleLin::get_stress()
{
	Eigen::MatrixXd sigma;
	return sigma;
}
