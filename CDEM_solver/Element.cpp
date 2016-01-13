#include "stdafx.h"
#include "Element.h"
#include "Node.h"
#include "aux_functions.h"

Element::Element()
{
	E = 0.0;
}

Element::~Element(){}

void Element::set_matrices(){}

void Element::print_self() 
{ 
	std::cout << "E = " << E << std::endl;
	std::cout << "nu = " << nu << std::endl;
	std::cout << "thickness = " << thickness << std::endl;
	std::cout << "density = " << density << std::endl;
	std::cout << "alfaC = " << alfaC << std::endl;
};

// Calculate the normal vector of each face and pass it to the adjacent nodes.
void Element::calc_normal_vectors()
{
	Eigen::VectorXd x(nnodes+1), y(nnodes+1);
	int i;
	for (i = 0; i < nnodes; i++)
	{
		x(i) = domain->nodes[nodes[i] - 1].x;
		y(i) = domain->nodes[nodes[i] - 1].y;
	}
	x(nnodes) = x(0);
	y(nnodes) = y(0);
	double dx, dy;
	Eigen::Vector2d vn;
	for (i = 0; i < nnodes; i++)
	{
		dx = x(i + 1) - x(i);
		dy = y(i + 1) - y(i);
		vn << dy, -dx;
		vn /= sqrt(dx*dx+dy*dy);
		domain->nodes[nodes[i % (nnodes - 1)] - 1].v_norm[1] = vn;
		domain->nodes[nodes[(i+1) % (nnodes - 1)] - 1].v_norm[0] = vn;
	}
}



// Perform one iteration of dynamic relaxation. Return the velocity norm.
double Element::iterate(double dt, double tau)
{
	Node ** nds;
	int i;
	int ndofs = stiffness_dim / nnodes;
	nds = new Node*[nnodes];
	for (i = 0; i < nnodes; i++)
	{
		nds[i] = &domain->nodes[nodes[i] - 1];
	}
	Eigen::VectorXd v_disp(stiffness_dim), v_velo(stiffness_dim), v_acce(stiffness_dim), v_supp(stiffness_dim);
	Eigen::VectorXd F_k_e(stiffness_dim), F_r(stiffness_dim), F_k_c(stiffness_dim), 
		F_c(stiffness_dim), F_ext(stiffness_dim), F_tot(stiffness_dim);
	// Get initial values.
	for (i = 0; i < stiffness_dim; i++)
	{
		v_disp(i) = nds[i / ndofs]->v_disp(i%ndofs);
		v_velo(i) = nds[i / ndofs]->v_velo(i%ndofs);
		v_acce(i) = nds[i / ndofs]->v_acce(i%ndofs);
		v_supp(i) = nds[i / ndofs]->supports[i%ndofs];
	}
	/*std::cout << "Displacement:" << std::endl << v_disp.transpose() << std::endl;*/

	// Get next values.
	v_disp = v_disp + dt * v_velo + 0.5 * dt * dt * v_acce;
	v_velo = v_velo + dt * v_acce;

	// Get forces.
	for (i = 0; i < stiffness_dim; i++)
	{
		F_k_c(i) = domain->get_contact_force(nodes[i / ndofs])(i%ndofs);
		F_ext(i) = load_function(tau) * nds[i / ndofs]->v_load[i%ndofs];
	}
	F_k_e = -1 * K_loc * v_disp;
	F_r = -1 * F_k_e.array() * v_supp.array();
	F_c = -1. * C_loc * v_velo;
	F_tot = F_k_e + F_k_c + F_c + F_ext + F_r;
	/*std::cout << "Stiffness force:" << std::endl << F_k_e.transpose() << std::endl;
	std::cout << "Reaction force:" << std::endl << F_r.transpose() << std::endl;
	std::cout << "Total force:" << std::endl << F_tot.transpose() << std::endl;*/

	v_acce = M_loc_inv * F_tot;
	for (i = 0; i < stiffness_dim; i++)
	{
		nds[i / ndofs]->v_disp(i%ndofs) = v_disp(i);
		nds[i / ndofs]->v_velo(i%ndofs) = v_velo(i);
		nds[i / ndofs]->v_acce(i%ndofs) = v_acce(i);
	}

	delete[] nds;
	return 0.0;
}
