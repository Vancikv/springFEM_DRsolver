#include "cuda_functions.cuh"

cudaError_t element_step_with_CUDA(Eigen::VectorXd * u, Eigen::VectorXd * v, Eigen::VectorXd * a,
	Eigen::VectorXd * load, Eigen::VectorXd * supports, int * neighbors, Eigen::MatrixXd * K,
	Eigen::MatrixXd * C, Eigen::MatrixXd * M, Eigen::Matrix2d Kc, int n_elems, int n_nodes)
{
	// Eigen structures will be copied into arrays of doubles
	double * u_;
	Eigen::Map<Eigen::VectorXd>(u_, u->rows(), u->cols()); // Layout is column by column
}
