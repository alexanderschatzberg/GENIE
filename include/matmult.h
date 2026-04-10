#ifndef PROPCA_MATMULT_H_
#define PROPCA_MATMULT_H_


#include "genotype.h"
#include "threadpool.h"

#include <Eigen/Dense>
#include <Eigen/Core>
using namespace Eigen;
//typedef Matrix<float, Dynamic, Dynamic, RowMajor> MatrixXdr;
//typedef Matrix<double, Dynamic, Dynamic, RowMajor> MatrixXdr;
#ifdef USE_DOUBLE
	typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixXdr;
#else
	typedef Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixXdr;
#endif

class MatMult {
 public:
 	genotype g;
	MatrixXdr geno_matrix;  // (p,n)

	bool debug; // = false;
	bool var_normalize; // = false;
	bool memory_efficient; // = false;
	bool missing; // = false;
	bool fast_mode; // = true;
	int nthreads = 1; // = 1;
	int Nindv_alloc = 0;  // allocated Nindv for y_e (to detect reuse)

	// How to batch columns:
	int blocksize = 0;  // k
	int hsegsize = 0;  // = log_3(n)
	int hsize = 0;
	int vsegsize = 0;  // = log_3(p)
	int vsize = 0;

	double **partialsums = nullptr;
	double *sum_op = nullptr;

	// Thread pool for reuse across multiply calls
	ForkJoinPool *pool = nullptr;

	// Intermediate computations in E-step.
	double **yint_e = nullptr;  // Size = nthreads X 3^(log_3(n)) * k
	double ***y_e = nullptr;    // nthreads X n X k (pointer arrays)
	double **y_e_data = nullptr; // nthreads contiguous blocks backing y_e

	// Intermediate computations in M-step.
	double **yint_m = nullptr;  // Size = nthreads X 3^(log_3(n)) * k
	double ***y_m = nullptr;    // nthreads X log_3(n) X k (pointer arrays)
	double **y_m_data = nullptr; // nthreads contiguous blocks backing y_m

	MatMult() {}
  	~MatMult();

	MatMult(const MatMult&) = delete;
	MatMult& operator=(const MatMult&) = delete;

	MatMult(MatMult&& other) noexcept;
	MatMult& operator=(MatMult&& other) noexcept;


	MatMult(genotype &xg,
			MatrixXdr &xgeno_matrix,
			bool xdebug,
			bool xvar_normalize,
			bool xmemory_efficient,
			bool xmissing,
			bool xfast_mode,
			int xnthreads,
			int xk);

	void multiply_y_pre_fast_thread(int begin, int end, MatrixXdr &op, int Ncol_op, double *yint_m, double **y_m, double *partialsums, MatrixXdr &res, double *y_m_data_flat);
	void multiply_y_post_fast_thread(int begin, int end, MatrixXdr &op, int Ncol_op, double *yint_e, double **y_e, double *partialsums, double *y_e_data_flat);

	/*
	 * Compute C = Y E 
	 * Y : p X n genotype matrix
	 * E : n K k matrix: X^{T} (XX^{T})^{-1}
	 * C = p X k matrix
	 *
	 * op : E
	 * Ncol_op : k
	 * res : C
	 * subtract_means :
	 */
	void multiply_y_pre_fast(MatrixXdr &op, int Ncol_op, MatrixXdr &res, bool subtract_means);

	/*
	 * Compute X = D Y 
	 * Y : p X n genotype matrix
	 * D : k X p matrix: (C^T C)^{-1} C^{T}
	 * X : k X n matrix
	 *
	 * op_orig : D
	 * Nrows_op : k
	 * res : X
	 * subtract_means :
	 */
	void multiply_y_post_fast(MatrixXdr &op_orig, int Nrows_op, MatrixXdr &res, bool subtract_means);
	void multiply_y_pre_naive_mem(MatrixXdr &op, int Ncol_op, MatrixXdr &res);
	void multiply_y_post_naive_mem(MatrixXdr &op, int Nrows_op, MatrixXdr &res);

	void multiply_y_post(MatrixXdr &op, int Nrows_op, MatrixXdr &res, bool subtract_means);

	void multiply_y_pre(MatrixXdr &op, int Ncol_op, MatrixXdr &res, bool subtract_means);

	void clean_up();

	// Reset for reuse: reuses existing buffers if dimensions match
	void reset(genotype &xg,
			MatrixXdr &xgeno_matrix,
			bool xdebug,
			bool xvar_normalize,
			bool xmemory_efficient,
			bool xmissing,
			bool xfast_mode,
			int xnthreads,
			int xk);

};


#endif  // PROPCA_MATMULT_H_
