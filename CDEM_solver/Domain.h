#pragma once
#include "stdafx.h"
#include "Element.h"
#include "ElementQuadrangleLin.h"
#include "Node.h"

class Domain
{
public:
	Domain(int, int, Node *, ElementQuadrangleLin *, double, double);
	Domain();
	~Domain();
	Eigen::Matrix2d m_contact_stiffness;
	int nelems, nnodes;
	ElementQuadrangleLin * elements;
	Node * nodes;
	// Create a domain structure defined by a text file
	void load_from_file(std::string filename);
	// Calculate the force acting on a node as a result of its relative displacement to the neighbor nodes.
	Eigen::Vector2d get_contact_force(int node_id);
	// Solve the system using the dynamic relaxation method.
	void solve(double t_load, double t_max, int maxiter);
};

