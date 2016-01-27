#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdio.h>
//#include "stdafx.h"
#include "Eigen/Dense"

cudaError_t element_step_with_CUDA(Eigen::VectorXd * u, Eigen::VectorXd * v, Eigen::VectorXd * a,
	Eigen::VectorXd * load, Eigen::VectorXd * supports, int * neighbors, Eigen::MatrixXd * K, 
	Eigen::MatrixXd * C, Eigen::MatrixXd * M, Eigen::Matrix2d Kc, int n_elems, int n_nodes);

__global__ void element_step_kernel(double * u, double * v, double * a,
	double * load, double * supports, int * neighbors, double * K, double * C, double * M,
	double * Kc,int n_elems, int n_nodes);