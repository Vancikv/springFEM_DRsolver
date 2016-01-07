#pragma once
#include "stdafx.h"

class Node
{
	int ndofs;
public:
	Node();
	~Node();
	double x;
	double y;
	Eigen::VectorXd v_disp;
	Eigen::VectorXd v_velo;
	Eigen::VectorXd v_acce;
	Eigen::VectorXd v_load;
	int * supports;
	int * neighbors;
	// Assign code numbers to node dofs. Increase the maxcode value accordingly.
	void set_codes(int maxcode);
	// Initiate nodal values prior to a dynamic relaxation calculation.
	void init_vals(double tau_0);
	// Set the displacement vector.
	void set_disp(double * val);
	// Set the velocity vector.
	void set_velo(double * val);
	// Set the acceleration vector.
	void set_acce(double * val);
};

