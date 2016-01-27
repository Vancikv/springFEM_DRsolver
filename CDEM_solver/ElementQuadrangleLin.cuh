#pragma once
#include "stdafx.h"
#include "Element.cuh"

class ElementQuadrangleLin :
	public Element
{
public:
	ElementQuadrangleLin(double, double, double, double, double, int *,int);
	ElementQuadrangleLin();
	~ElementQuadrangleLin();
	// Calculate local stiffness matrix. Use reduced integration with hourglass stabilization.
	void set_K_isoparametric();
	// Calculate stress vector for element gauss points.
	void set_matrices();
	Eigen::MatrixXd get_stress();
};

