#include "stdafx.h"
#include "aux_functions.h"
#include <cmath>

double load_function(double tau)
{
	if (tau < 1.0)
	{
		return coef1*(tau*tau*tau) + coef2 * (tau*tau);
	}
	else
	{
		return 1.0;
	}
}

double vector_norm(double * vec, int size)
{
	double norm = 0.0;
	for (int i = 0; i < size; i++)
		norm += vec[i] * vec[i];
	norm = sqrt(norm);
	return norm;
}