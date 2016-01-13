#pragma once
#include "stdafx.h"
#include "Domain.h"

class Domain;
class Element
{
public:
	Element();
	virtual ~Element();
	double E;
	double nu;
	double density;
	double thickness;
	double alfaC;
	Domain * domain;
	int nnodes;
	int * nodes;
	int stiffness_dim;
	Eigen::MatrixXd K_loc;
	Eigen::MatrixXd * B_matrices;
	Eigen::MatrixXd M_loc;
	Eigen::MatrixXd M_loc_inv;
	Eigen::MatrixXd C_loc;
	double volume;
	// print data
	void print_self();
	// Calculate and store local matrices
	virtual void set_matrices();
	// Calculate the normal vector of each face and pass it to the adjacent nodes.
	void calc_normal_vectors();
	// Perform one iteration of dynamic relaxation. Return the velocity norm.
	double iterate(double dt, double tau);
};
