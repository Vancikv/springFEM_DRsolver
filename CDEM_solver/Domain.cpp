#include "stdafx.h"
#include "Domain.h"
#include "Element.h"
#include "Node.h"


Domain::Domain(int _nelem, int _nnode, Node * _nodes, Element * _elements, double c_n, double c_s)
{
	nelem = _nelem;
	nnode = _nnode;
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
	Eigen::Vector2d du_g;
	Eigen::Matrix2d T(2, 2);
	Node n0 = nodes[node_id-1];
	int i;
	for (i = 0; i < 2; i++)
	{
		if (n0.neighbors[i] != 0)
		{
			T << n0.v_norm[i](0), n0.v_norm[i](1),
				-n0.v_norm[i](1), n0.v_norm[i](0);
			du_g = nodes[n0.neighbors[i] - 1].v_disp - n0.v_disp;
			F += T.transpose() * m_contact_stiffness * T * du_g;
		}
	}
	return F;
}


// Solve the system using the dynamic relaxation method.
void Domain::solve(double t_load, double t_max, int maxiter)
{
	double dt = t_max / maxiter;
	int i, j;
	for (i = 0; i < nelem; i++)
	{
		elements[i].set_matrices();
		elements[i].calc_normal_vectors();
	}
	for (i = 0; i < nnode; i++)
	{
		nodes[i].init_vals(dt);
	}
	int itercount = 0;
	for (i = 0; i < itercount; i++)
	{
		for (j = 0; j < nelem; j++)
		{
			elements[j].iterate(dt, itercount * dt / t_load);
		}
	}
}
