#include "stdafx.h"
#include "Domain.h"
#include "ElementQuadrangleLin.h"
#include "Node.h"
#include <fstream>
#include <istream>
#include "aux_functions.h"
#include <cmath>
#include <ctime>

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

void Domain::write_state_to_file(std::string filename, double time)
{
	std::ofstream outfile(filename);
	int i;

	outfile << nnodes << " " << nelems << " " << time << std::endl;
	for (i = 0; i < nnodes; i++)
	{
		outfile << "node " << (i+1) << " " << nodes[i].x << " " << nodes[i].y << " " << nodes[i].v_disp(0) << " " << nodes[i].v_disp(1)
			<< " " << nodes[i].v_velo(0) << " " << nodes[i].v_velo(0) << " " << nodes[i].v_acce(0) << " " << nodes[i].v_acce(0) << std::endl;
	}
	for (i = 0; i < nelems; i++)
	{
		Element& e = elements[i];
		outfile << "element " << (i + 1) << " " << elements[i].nnodes;
		for (int j = 0; j < elements[i].nnodes; j++)
		{
			outfile << " " << elements[i].nodes[j];
		}
		outfile << std::endl;
	}
}

void Domain::write_state_to_file(std::string filename, double time, double loadfunc)
{
	std::ofstream outfile(filename);
	int i;

	outfile << nnodes << " " << nelems << " " << time << " " << loadfunc << std::endl;
	for (i = 0; i < nnodes; i++)
	{
		outfile << "node " << (i + 1) << " " << nodes[i].x << " " << nodes[i].y << " " << nodes[i].v_disp(0) << " " << nodes[i].v_disp(1)
			<< " " << nodes[i].v_velo(0) << " " << nodes[i].v_velo(0) << " " << nodes[i].v_acce(0) << " " << nodes[i].v_acce(0) << std::endl;
	}
	for (i = 0; i < nelems; i++)
	{
		Element& e = elements[i];
		outfile << "element " << (i + 1) << " " << elements[i].nnodes;
		for (int j = 0; j < elements[i].nnodes; j++)
		{
			outfile << " " << elements[i].nodes[j];
		}
		outfile << std::endl;
	}
}

void Domain::write_state_to_file(std::string filename, double time, double loadfunc, double *u_norm, double * v_norm, double *a_norm, int normsize)
{
	std::ofstream outfile(filename);
	int i;

	outfile << nnodes << " " << nelems << " " << time << " " << loadfunc << std::endl;
	for (i = 0; i < nnodes; i++)
	{
		outfile << "node " << (i + 1) << " " << nodes[i].x << " " << nodes[i].y << " " << nodes[i].v_disp(0) << " " << nodes[i].v_disp(1)
			<< " " << nodes[i].v_velo(0) << " " << nodes[i].v_velo(0) << " " << nodes[i].v_acce(0) << " " << nodes[i].v_acce(0) << std::endl;
	}
	for (i = 0; i < nelems; i++)
	{
		Element& e = elements[i];
		outfile << "element " << (i + 1) << " " << elements[i].nnodes;
		for (int j = 0; j < elements[i].nnodes; j++)
		{
			outfile << " " << elements[i].nodes[j];
		}
		outfile << std::endl;
	}
	outfile << "u_norm";
	for (i = 0; i < normsize; i++)
	{
		outfile << " " << u_norm[i];
	}
	outfile << std::endl;
	outfile << "v_norm";
	for (i = 0; i < normsize; i++)
	{
		outfile << " " << v_norm[i];
	}
	outfile << std::endl;
	outfile << "a_norm";
	for (i = 0; i < normsize; i++)
	{
		outfile << " " << a_norm[i];
	}
	outfile << std::endl;
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
void Domain::solve(double t_load, double t_max, int maxiter, std::string outfile, int output_frequency, int n_inner_steps)
{	
	int i, j, k;
	std::string timefile_name = outfile+ "_time.txt";
	clock_t begin_local_assembly = clock();
	for (i = 0; i < nelems; i++)
	{
		elements[i].set_matrices();
		elements[i].calc_normal_vectors();
	}
	clock_t end_local_assembly = clock();

	// Eigen matrices will be copied into arrays of doubles.
	// Using the Eigen::Map function defaults in a column by column layout.
	double * u, *v, *a,*z, *load, *supports, *K, *C, *Mi, *Kc, *n_vects;
	int * neighbors;
	int nnodedofs = 2, stiffdim = 8;
	int vdim = nnodes*nnodedofs;
	u = new double[vdim];
	v = new double[vdim];
	a = new double[vdim];
	z = new double[vdim];
	load = new double[vdim];
	supports = new double[vdim];
	K = new double[nelems*stiffdim*stiffdim];
	C = new double[nelems*stiffdim];
	Mi = new double[nelems*stiffdim];
	Kc = new double[4];
	n_vects = new double[4 * nnodes];
	neighbors = new int[vdim];

	for (i = 0; i < nnodes; i++)
	{
		int of = i*nnodedofs;
		Eigen::Map<Eigen::VectorXd>(u + of, nodes[i].v_disp.rows(), nodes[i].v_disp.cols()) = nodes[i].v_disp;
		Eigen::Map<Eigen::VectorXd>(v + of, nodes[i].v_velo.rows(), nodes[i].v_velo.cols()) = nodes[i].v_velo;
		Eigen::Map<Eigen::VectorXd>(a + of, nodes[i].v_acce.rows(), nodes[i].v_acce.cols()) = nodes[i].v_acce;
		Eigen::Map<Eigen::VectorXd>(load + of, nodes[i].v_load.rows(), nodes[i].v_load.cols()) = nodes[i].v_load;
		*(supports++) = nodes[i].supports[0];
		*(supports++) = nodes[i].supports[1];
		*(neighbors++) = nodes[i].neighbors[0];
		*(neighbors++) = nodes[i].neighbors[1];
		*(n_vects++) = nodes[i].v_norm[0](0);
		*(n_vects++) = nodes[i].v_norm[0](1);
		*(n_vects++) = nodes[i].v_norm[1](0);
		*(n_vects++) = nodes[i].v_norm[1](1);
	}

	supports -= vdim;
	neighbors -= vdim;
	n_vects -= 4 * nnodes;

	for (i = 0; i < nelems; i++)
	{
		int of1 = i*stiffdim*stiffdim, of2 = i*stiffdim;
		Eigen::Map<Eigen::MatrixXd>(K + of1, elements[i].K_loc.rows(), elements[i].K_loc.cols()) = elements[i].K_loc;
		Eigen::Map<Eigen::MatrixXd>(C + of2, elements[i].C_loc.diagonal().rows(), elements[i].C_loc.diagonal().cols()) = elements[i].C_loc.diagonal();
		Eigen::Map<Eigen::MatrixXd>(Mi + of2, elements[i].M_loc_inv.diagonal().rows(), elements[i].M_loc_inv.diagonal().cols()) = elements[i].M_loc_inv.diagonal();
	}

	Eigen::Map<Eigen::MatrixXd>(Kc, m_contact_stiffness.rows(), m_contact_stiffness.cols()) = m_contact_stiffness;

	clock_t end_eigen_to_double = clock();

	const double c = 0.01; // Load speed constant

	double F_k_e = 0;
	double F_k_c = 0;
	double F_c = 0;
	double F_r = 0;
	double * last_u, * last_v;
	last_u = new double[vdim];
	last_v = new double[vdim];
	double * u_norm, *v_norm, *a_norm;
	//int n_inner_steps=2;
	u_norm = new double[n_inner_steps+1];
	v_norm = new double[n_inner_steps+1];
	a_norm = new double[n_inner_steps+1];

	double loadfunc=0.0;

	// Estimate max time step length, take greatest stiffness and smallest mass:

	double dt = t_max / maxiter;
	double zT_Mi_p, pT_Mi_p;

	clock_t begin_dr = clock();
	for (k = 1; k <= maxiter; k++) // Loop of time steps
	{
		zT_Mi_p = 0.0, pT_Mi_p = 0.0;
		if (k*dt < t_load)
		{
			for (i = 0; i < nnodedofs*nnodes; i++) // Loop of dofs - load function
			{
				zT_Mi_p += z[i] * Mi[i] * load[i];
				pT_Mi_p += Mi[i] * load[i] * load[i];
			}
			loadfunc = zT_Mi_p / pT_Mi_p + c * (t_load - k*dt) / t_load;
		}else{
			loadfunc = 1.0;
		}


		for (i = 0; i < nnodedofs*nnodes; i++) // Loop of dofs - increment displacement
		{
			double buff = u[i];
			u[i] = 2*dt*v[i] + last_u[i];
			last_u[i] = buff;
			buff = v[i];
			v[i] = 2 * dt*a[i] + last_v[i];
			last_v[i] = buff;
		}

		for (i = 0; i < nnodedofs*nnodes; i++) // Loop of dofs - Calculate balance of forces and resulting acceleration
		{
			int eid = i / stiffdim; // global number of element
			int nid = (i / nnodedofs) * nnodedofs; // number of dof 1 of this node
			int ned = i % stiffdim; // number of dof within element
			int mdim = stiffdim*stiffdim; // number of elements of the stiffness matrix
			double kc11 = Kc[0];
			double kc21 = Kc[1];
			double kc12 = Kc[2];
			double kc22 = Kc[3];

			// Element stiffness force:
			F_k_e = 0;
			for (j = 0; j < stiffdim; j++)
			{
				F_k_e += -K[eid*mdim + j*stiffdim + ned] * u[eid*stiffdim + j];
			}
			// Contact stiffness force:
			F_k_c = 0;
			for (j = 0; j < 2; j++)
			{
				int nbr = neighbors[nid + j];
				if (nbr != 0)
				{
					double t11 = n_vects[4 * (i / nnodedofs) + 2 * j];
					double t12 = n_vects[4 * (i / nnodedofs) + 2 * j + 1];
					double t21 = -t12;
					double t22 = t11;
					double du_x = u[(nbr - 1)*nnodedofs] - u[nid];
					double du_y = u[(nbr - 1)*nnodedofs + 1] - u[nid + 1];
					if (i == nid) // X-component
					{
						F_k_c += du_x * (t11*(t11*kc11 + t21*kc21) + t21*(t11*kc12 + t21*kc22)) + du_y * (t12*(t11*kc11 + t21*kc21) + t22*(t11*kc12 + t21*kc22)); // T_T * Kc * T * du_g
					}
					else // Y-component
					{
						F_k_c += du_x * (t11*(t12*kc11 + t22*kc21) + t21*(t12*kc12 + t22*kc22)) + du_y * (t12*(t12*kc11 + t22*kc21) + t22*(t12*kc12 + t22*kc22)); // T_T * Kc * T * du_g
					}
				}
			}
			// Damping force:
			F_c = -C[i] * v[i];
			// Reaction force
			F_r = supports[i] * (-F_k_e - F_k_c - F_c - loadfunc*load[i]);
			z[i] = F_k_e + F_k_c + F_r + F_c;
			a[i] = Mi[i] * (z[i] + loadfunc*load[i]);
		}
		u_norm[0] = vector_norm(u, vdim);
		v_norm[0] = vector_norm(v, vdim);
		a_norm[0] = vector_norm(a, vdim);

		if ((k % output_frequency == 0) || (k == 1))
		{
			std::string num = std::to_string(k);
			std::cout << loadfunc << std::endl;
			num = num + ".txt";
			for (i = 0; i < nnodedofs*nnodes; i++) // Loop of dofs - update nodal values
			{
				int nid = i / nnodedofs, did = i%nnodedofs; // node id, nodal dof id
				nodes[nid].v_disp(did) = u[i];
				nodes[nid].v_velo(did) = v[i];
				nodes[nid].v_acce(did) = a[i];
			}
			write_state_to_file(outfile + num, k * dt, loadfunc, u_norm, v_norm, a_norm, n_inner_steps + 1);
		}
	}
	clock_t end_dr = clock();
	std::ofstream timefile(timefile_name);
	timefile << "Local assemblies: " << std::endl << double(end_local_assembly-begin_local_assembly)/CLOCKS_PER_SEC << " seconds" << std::endl << double(end_local_assembly - begin_local_assembly) << " clicks" << std::endl;
	timefile << "Eigen to double: " << std::endl << double(end_eigen_to_double - end_local_assembly) / CLOCKS_PER_SEC << " seconds" << std::endl << double(end_eigen_to_double - end_local_assembly) << " clicks" << std::endl;
	timefile << "Dynamic relaxation solution: " << std::endl << double(end_dr - begin_dr) / CLOCKS_PER_SEC << " seconds" << std::endl << double(end_dr - begin_dr) << " clicks" << std::endl;
}

void Domain::solve(double t_load, double t_max, int maxiter)
{
	double dt = t_max / maxiter;
	int i, j;
	for (i = 0; i < nelems; i++)
	{
		elements[i].set_matrices();
		elements[i].calc_normal_vectors();
		for (j = 0; j < elements[i].nnodes; j++)
		{
			nodes[elements[i].nodes[j] - 1].init_vals(dt / t_max, elements[i].M_loc(2 * j, 2 * j));
		}
	}

	// Eigen matrices will be copied into arrays of doubles.
	// Using the Eigen::Map function defaults in a column by column layout.
	double * u, *v, *a, *load, *supports, *K, *C, *Mi, *Kc, *n_vects;
	int * neighbors;
	int nnodedofs = 2, stiffdim = 8;
	int vdim = nnodes*nnodedofs;
	u = new double[vdim];
	v = new double[vdim];
	a = new double[vdim];
	load = new double[vdim];
	supports = new double[vdim];
	K = new double[nelems*stiffdim*stiffdim];
	C = new double[nelems*stiffdim];
	Mi = new double[nelems*stiffdim];
	Kc = new double[4];
	n_vects = new double[4 * nnodes];
	neighbors = new int[vdim];

	for (i = 0; i < nnodes; i++)
	{
		int of = i*nnodedofs;
		Eigen::Map<Eigen::VectorXd>(u + of, nodes[i].v_disp.rows(), nodes[i].v_disp.cols()) = nodes[i].v_disp;
		Eigen::Map<Eigen::VectorXd>(v + of, nodes[i].v_velo.rows(), nodes[i].v_velo.cols()) = nodes[i].v_velo;
		Eigen::Map<Eigen::VectorXd>(a + of, nodes[i].v_acce.rows(), nodes[i].v_acce.cols()) = nodes[i].v_acce;
		Eigen::Map<Eigen::VectorXd>(load + of, nodes[i].v_load.rows(), nodes[i].v_load.cols()) = nodes[i].v_load;
		*(supports++) = nodes[i].supports[0];
		*(supports++) = nodes[i].supports[1];
		*(neighbors++) = nodes[i].neighbors[0];
		*(neighbors++) = nodes[i].neighbors[1];
		*(n_vects++) = nodes[i].v_norm[0](0);
		*(n_vects++) = nodes[i].v_norm[0](1);
		*(n_vects++) = nodes[i].v_norm[1](0);
		*(n_vects++) = nodes[i].v_norm[1](1);

	}

	supports -= vdim;
	neighbors -= vdim;
	n_vects -= 4 * nnodes;

	for (i = 0; i < nelems; i++)
	{
		int of1 = i*stiffdim*stiffdim, of2 = i*stiffdim;
		Eigen::Map<Eigen::MatrixXd>(K + of1, elements[i].K_loc.rows(), elements[i].K_loc.cols()) = elements[i].K_loc;
		Eigen::Map<Eigen::MatrixXd>(C + of2, elements[i].C_loc.diagonal().rows(), elements[i].C_loc.diagonal().cols()) = elements[i].C_loc.diagonal();
		Eigen::Map<Eigen::MatrixXd>(Mi + of2, elements[i].M_loc_inv.diagonal().rows(), elements[i].M_loc_inv.diagonal().cols()) = elements[i].M_loc_inv.diagonal();
	}

	Eigen::Map<Eigen::MatrixXd>(Kc, m_contact_stiffness.rows(), m_contact_stiffness.cols()) = m_contact_stiffness;

	for (int k = 0; k < maxiter; k++) // Loop of time steps
	{
		double loadfunc = load_function(k*dt / t_load);
		for (i = 0; i < nnodedofs*nnodes; i++) // Loop of dofs - increment displacement
		{
			u[i] += dt*v[i] + 0.5*dt*dt*a[i];
			v[i] += dt*a[i];
		}
		for (i = 0; i < nnodedofs*nnodes; i++) // Loop of dofs - Calculate balance of forces and resulting acceleration
		{
			int eid = i / stiffdim; // global number of element
			int nid = (i / nnodedofs) * nnodedofs; // number of dof 1 of this node
			int ned = i % stiffdim; // number of dof within element
			int mdim = stiffdim*stiffdim; // number of elements of the stiffness matrix
			double kc11 = Kc[0];
			double kc21 = Kc[1];
			double kc12 = Kc[2];
			double kc22 = Kc[3];

			// Element stiffness force:
			double F_k_e = 0;
			for (j = 0; j < stiffdim; j++)
			{
				F_k_e += -K[eid*mdim + j*stiffdim + ned] * u[eid*stiffdim + j];
			}
			// Contact stiffness force:
			double F_k_c = 0;
			for (j = 0; j < 2; j++)
			{
				int nbr = neighbors[nid + j];
				if (nbr != 0)
				{
					double t11 = n_vects[4 * (i / nnodedofs) + 2 * j];
					double t12 = n_vects[4 * (i / nnodedofs) + 2 * j + 1];
					double t21 = -t12;
					double t22 = t11;
					double du_x = u[(nbr - 1)*nnodedofs] - u[nid];
					double du_y = u[(nbr - 1)*nnodedofs + 1] - u[nid + 1];
					if (i == nid) // X-component
					{
						F_k_c += du_x * (t11*(t11*kc11 + t21*kc21) + t21*(t11*kc12 + t21*kc22)) + du_y * (t12*(t11*kc11 + t21*kc21) + t22*(t11*kc12 + t21*kc22)); // T_T * Kc * T * du_g
					}
					else // Y-component
					{
						F_k_c += du_x * (t11*(t12*kc11 + t22*kc21) + t21*(t12*kc12 + t22*kc22)) + du_y * (t12*(t12*kc11 + t22*kc21) + t22*(t12*kc12 + t22*kc22)); // T_T * Kc * T * du_g
					}
				}
			}
			// Damping force:
			double F_c = -C[i] * v[i];
			// Reaction force
			double F_r = supports[i] * (-F_k_e - F_k_c - F_c - loadfunc*load[i]);
			a[i] = Mi[i] * (F_k_e + F_k_c + F_r + F_c + loadfunc*load[i]);
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