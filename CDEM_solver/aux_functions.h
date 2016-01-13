#pragma once
const double parX = 0.8;
const double parY = 0.5;
const double coef1 = (parY-parX*parX) / (parX*parX*parX - parX*parX);
const double coef2 = 1.0 - coef1;

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

