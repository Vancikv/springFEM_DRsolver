#include "stdafx.h"
#include "Domain.h"
#include "ElementQuadrangleLin.h"
#include "Node.h"
#include <fstream>
#include <istream>

Domain::Domain(int _nelems, int _nnodes, Node * _nodes, Element * _elements, double c_n, double c_s)
{
	nelems = _nelems;
	nnodes = _nnodes;
	nodes = _nodes;
	elements = _elements;
	m_contact_stiffness << c_n, 0.0,
		0.0, c_s;
}

Domain::Domain()
{
	m_contact_stiffness << 0.0, 0.0,
		0.0, 0.0;
}

Domain::~Domain()
{
	delete[] elements;
	delete[] nodes;
}

void Domain::write_state_to_file(std::string filename)
{

}

// Calculate the force acting on a node as a result of its relative displacement to the neighbor nodes.
Eigen::Vector2d Domain::get_contact_force(int node_id)
{
	Eigen::Vector2d F(0., 0.);
	Eigen::Vector2d du_g;
	Eigen::Matrix2d T;
	Node& n0 = nodes[node_id-1];
	int i;
	for (i = 0; i < 2; i++)
	{
		if (n0.neighbors[i] != 0)
		{
			T << n0.v_norm[i](0), n0.v_norm[i](1),
				-n0.v_norm[i](1), n0.v_norm[i](0);
			du_g = nodes[n0.neighbors[i] - 1].v_disp - n0.v_disp; // 2x1 - 2x1 = 2x1
			F += T.transpose() * m_contact_stiffness * T * du_g; // 2x2 * 2x2 * 2x2 * 2x1 = 2x1
		}
	}
	return F;
}


// Solve the system using the dynamic relaxation method.
void Domain::solve(double t_load, double t_max, int maxiter)
{	
	double dt = t_max / maxiter;
	int i, j;
	for (i = 0; i < nelems; i++)
	{
		elements[i].set_matrices();
		elements[i].calc_normal_vectors();
	}
	for (i = 0; i < nnodes; i++)
	{
		nodes[i].init_vals(dt);
	}

	for (i = 0; i < maxiter; i++)
	{
		for (j = 0; j < nelems; j++)
		{
			elements[j].iterate(dt, i * dt / t_load);
		}
	}
}

void Domain::load_from_file(std::string filename)
{
	std::ifstream input_file (filename);
	std::string line;
	int lncount=0;

	while (std::getline(input_file, line)) // Implement error catching!
	{
		lncount++;
		std::string entry;
		std::stringstream lss(line);
		if (lncount == 1) // First row - domain specs
		{
			while (std::getline(lss,entry,' '))
			{
				if (entry == "nnodes") lss >> nnodes;
				if (entry == "nelems") lss >> nelems;
				if (entry == "cn") lss >> m_contact_stiffness(0,0);
				if (entry == "cs") lss >> m_contact_stiffness(1,1);
			}
			elements = new ElementQuadrangleLin[nelems];
			nodes = new Node[nnodes];
		}
		else if (lncount <= 1 + nnodes) // Node records
		{
			Node& nd = nodes[lncount - 2];
			while (std::getline(lss, entry, ' '))
			{
				if (entry == "ndofs")
				{
					int ndf, i;
					lss >> ndf;
					nd.ndofs = ndf;
					nd.neighbors = new int[ndf];
					nd.supports = new int[ndf];
					nd.v_load = Eigen::VectorXd::Zero(ndf);
					nd.v_disp = Eigen::VectorXd::Zero(ndf);
					nd.v_velo = Eigen::VectorXd::Zero(ndf);
					nd.v_acce = Eigen::VectorXd::Zero(ndf);
					nd.v_code = new int[ndf];
					for (i = 0; i < ndf; i++)
					{
						nd.neighbors[i] = 0;
						nd.supports[i] = 0;
						nd.v_code[i] = 0;
					}
				}
				if (entry == "position")
				{
					lss >> nd.x;
					lss >> nd.y;
				}
				if (entry == "neighbors")
				{
					lss >> nd.neighbors[0];
					lss >> nd.neighbors[1];
				}
				if (entry == "supports")
				{
					int i;
					for (i = 0; i < nd.ndofs;i++) lss >> nd.supports[i];
				}
				if (entry == "load")
				{
					int i;
					for (i = 0; i < nd.ndofs; i++) lss >> nd.v_load(i);
				}
			}
		}
		else // Element records
		{
			Element& el = elements[lncount - 2 - nnodes];
			el.domain = this;
			while (std::getline(lss, entry, ' '))
			{
				if (entry == "nodes")
				{
					int nnds, i;
					lss >> nnds;
					el.nnodes = nnds;
					el.nodes = new int[nnds];
					for (i = 0; i < nnds; i++)
					{
						lss >> el.nodes[i];
					}
				}
				if (entry == "E") lss >> el.E;
				if (entry == "nu") lss >> el.nu;
				if (entry == "density") lss >> el.density;
				if (entry == "thickness") lss >> el.thickness;
				if (entry == "alfaC") lss >> el.alfaC;
			}
		}
	}

}