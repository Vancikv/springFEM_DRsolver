#include "stdafx.h"
#include "Domain.h"


Domain::Domain()
{
	//  c_n = 0.0;
	//  c_s = 0.0;
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
