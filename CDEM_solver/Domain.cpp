#include "stdafx.h"
#include "Domain.h"
#include "Element.h"
#include "Node.h"


Domain::Domain(Node * _nodes, Element * _elements, double c_n, double c_s)
{
	nodes = _nodes;
	elements = _elements;
	m_contact_stiffness << c_n, 0.0,
		0.0, c_s;
}


Domain::~Domain()
{
}


// Calculate the force acting on a node as a result of its relative displacement to the neighbor nodes.
Eigen::Vector2d Domain::get_contact_force(int node_id)
{
	Eigen::Vector2d F(0., 0.);
	return F;
}


// Solve the system using the dynamic relaxation method.
void Domain::solve(double t_load, double t_max, int maxiter)
{
}
