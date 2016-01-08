#pragma once
#include "stdafx.h"
#include "Element.h"
#include "Node.h"

class Domain
{
public:
	Domain(Node *, Element *, double, double);
	~Domain();
	Eigen::Matrix2d m_contact_stiffness;
	Element * elements;
	Node * nodes;
	// Calculate the force acting on a node as a result of its relative displacement to the neighbor nodes.
	Eigen::Vector2d get_contact_force(int node_id);
	// Solve the system using the dynamic relaxation method.
	void solve(double t_load, double t_max, int maxiter);
};

