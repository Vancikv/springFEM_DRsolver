#include "stdafx.h"
#include "Node.h"

Node::Node()
{}

Node::Node(int _ndofs, double _x, double _y, int * _supports, int * _neighbors, Eigen::VectorXd _v_load)
{
	ndofs = _ndofs;
	x = _x;
	y = _y;
	supports = _supports;
	neighbors = _neighbors;
	v_load = _v_load;
	v_disp = Eigen::VectorXd::Zero(ndofs);
	v_velo = Eigen::VectorXd::Zero(ndofs);
	v_acce = Eigen::VectorXd::Zero(ndofs);
	v_code = new double [ndofs];
	int i;
	for (i = 0; i < ndofs; i++)
	{
		v_code[i] = 0;
	}
}


Node::~Node()
{
	delete [] v_code;
}


// Assign code numbers to node dofs. Increase the maxcode value accordingly.
void Node::set_codes(int &maxcode)
{
	int i;
	for (i = 0; i < ndofs; i++)
	{
		if (supports[i] == 0)
		{
			maxcode++;
			v_code[i] = maxcode;
		}
	}
}


// Initiate nodal values prior to a dynamic relaxation calculation.
void Node::init_vals(double tau_0)
{
	int i;
	for (i = 0; i < ndofs; i++)
	{
		if (supports[i] == 0)
		{
			// v_acce(i) = load_function(tau_0) / mass * v_load(i)
			// define mass and load_function
		}
	}
}

// Set the displacement vector.
void Node::set_disp(Eigen::VectorXd val)
{
	v_disp = val;
}

// Set the velocity vector.
void Node::set_velo(Eigen::VectorXd val)
{
	v_velo = val;
}

// Set the acceleration vector.
void Node::set_acce(Eigen::VectorXd val)
{
	v_acce = val;
}
