#include "stdafx.h"
#include "Element.h"

Element::Element()
{
	E = 0.0;
}


Element::~Element()
{
}

// Calculate the normal vector of each face and pass it to the adjacent nodes.
void Element::calc_normal_vector()
{
}

// Calculate and store local matrices
void Element::set_matrices()
{
}

// Perform one iteration of dynamic relaxation. Return the velocity norm.
double Element::iterate(double dt, double tau)
{
	return 0.0;
}
