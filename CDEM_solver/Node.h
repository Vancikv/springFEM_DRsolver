#pragma once
#include "stdafx.h"

class Node
{
public:
	Node();
	Node(int, double, double, int *, int *, Eigen::VectorXd);
	~Node();
	int ndofs;
	double x;
	double y;
	Eigen::VectorXd v_disp;
	Eigen::VectorXd v_velo;
	Eigen::VectorXd v_acce;
	int * v_code;
	Eigen::Vector2d v_norm[2];
	Eigen::VectorXd v_load;
	int * supports;
	int * neighbors;
	// Assign code numbers to node dofs. Increase the maxcode value accordingly.
	void set_codes(int &maxcode);
	// Initiate nodal values prior to a dynamic relaxation calculation.
	void init_vals(double tau_0, double mass);
	// Set the displacement vector.
	void set_disp(Eigen::VectorXd);
	// Set the velocity vector.
	void set_velo(Eigen::VectorXd);
	// Set the acceleration vector.
	void set_acce(Eigen::VectorXd);
};

