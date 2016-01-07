#pragma once
#include "stdafx.h"

class Element
{
public:
	Element();
	~Element();
	double E;
	double nu;
	double density;
	double thickness;
	double alfaC;
	void * domain;
	int * nodes;
	Eigen::MatrixXd K_loc;
	Eigen::MatrixXd M_loc;
	// Calculate and store local matrices
	virtual void set_matrices();
	// Calculate the normal vector of each face and pass it to the adjacent nodes.
	void calc_normal_vector();
	// Perform one iteration of dynamic relaxation. Return the velocity norm.
	double iterate(double dt, double tau);
};
