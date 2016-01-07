#include "stdafx.h"
#include "Node.h"


Node::Node()
{
	x = 0.0;
	y = 0.0;
	ndofs = 0;
}


Node::~Node()
{
}


// Assign code numbers to node dofs. Increase the maxcode value accordingly.
void Node::set_codes(int maxcode)
{
}


// Initiate nodal values prior to a dynamic relaxation calculation.
void Node::init_vals(double tau_0)
{
}


// Set the displacement vector.
void Node::set_disp(double * val)
{
}

// Set the velocity vector.
void Node::set_velo(double * val)
{
}

// Set the acceleration vector.
void Node::set_acce(double * val)
{
}
