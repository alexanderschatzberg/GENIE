#include "matmult.h"
#include "mailman.h"
#include "helper.h"
#include "profiler.h"

#include <iostream>
#include <thread>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Core>

#if defined(SSE_SUPPORT) && SSE_SUPPORT == 1
	#define fastmultiply fastmultiply_sse
	#define fastmultiply_pre fastmultiply_pre_sse
#else
	#define fastmultiply fastmultiply_normal
	#define fastmultiply_pre fastmultiply_pre_normal
#endif

MatMult::MatMult(genotype &xg,
			MatrixXdr &xgeno_matrix,
			bool xdebug,
			bool xvar_normalize,
			bool xmemory_efficient,
			bool xmissing,
			bool xfast_mode,
			int xnthreads,
			int xk) {
	ScopedTimer timer("matmult_setup");
	g = xg;
	geno_matrix = xgeno_matrix;

	debug = xdebug;
	var_normalize = xvar_normalize;
	memory_efficient = xmemory_efficient;
	missing = xmissing;
	fast_mode = xfast_mode;
	nthreads = xnthreads;

	blocksize = xk;

	hsegsize = g.segment_size_hori;  // = log_3(n)
	hsize = pow(3, hsegsize);
	vsegsize = g.segment_size_ver;  // = log_3(p)
	vsize = pow(3, vsegsize);
	Nindv_alloc = g.Nindv;

	if (nthreads > 1) {
		pool = new ForkJoinPool(nthreads);
	}

	partialsums = new double*[nthreads];
	yint_m = new double*[nthreads];
	for (int t = 0; t < nthreads; t++) {
		partialsums[t] = new double[blocksize];
		memset(partialsums[t], 0, blocksize * sizeof(double));
		yint_m[t] = new double[hsize*blocksize];
		memset (yint_m[t], 0, hsize*blocksize * sizeof(double));
	}

	sum_op = new double[blocksize];
	memset(sum_op, 0, blocksize * sizeof(double));

	yint_e = new double* [nthreads];
	for (int t = 0; t < nthreads; t++) {
		yint_e[t] = new double[hsize*blocksize];
		memset (yint_e[t], 0, hsize*blocksize * sizeof(double));
	}

	// Flat allocation for y_e: one contiguous block per thread
	y_e  = new double**[nthreads];
	y_e_data = new double*[nthreads];
	for (int t = 0 ; t < nthreads ; t++) {
		y_e_data[t] = new double[static_cast<size_t>(g.Nindv) * blocksize]();
		y_e[t] = new double*[g.Nindv];
		for (int i = 0 ; i < g.Nindv ; i++) {
			y_e[t][i] = y_e_data[t] + static_cast<size_t>(i) * blocksize;
		}
	}

	// Flat allocation for y_m: one contiguous block per thread
	y_m = new double**[nthreads];
	y_m_data = new double*[nthreads];
	for (int t = 0; t < nthreads; t++) {
		y_m_data[t] = new double[static_cast<size_t>(hsegsize) * blocksize]();
		y_m[t] = new double*[hsegsize];
		for (int i = 0; i < hsegsize; i++) {
			y_m[t][i] = y_m_data[t] + static_cast<size_t>(i) * blocksize;
		}
	}
}

void MatMult::multiply_y_pre_fast_thread(int begin, int end, MatrixXdr &op, int Ncol_op, double *yint_m, double **y_m, double *partialsums, MatrixXdr &res, double *y_m_data_flat) {
	if (debug) {
		cout << "begin = " << begin << endl;
		cout << "end = " << end << endl;
	}
	for (int seg_iter = begin; seg_iter < end; seg_iter++) {
		memset(partialsums, 0, blocksize * sizeof(double));
		memset(yint_m, 0, hsize * blocksize * sizeof(double));
		memset(y_m_data_flat, 0, static_cast<size_t>(g.segment_size_hori) * blocksize * sizeof(double));
		mailman::fastmultiply(g.segment_size_hori, g.Nindv, Ncol_op, g.p[seg_iter], op, yint_m, partialsums, y_m);
		int p_base = seg_iter * g.segment_size_hori;
		for (int p_iter = p_base; (p_iter < p_base + g.segment_size_hori) && (p_iter < g.Nsnp); p_iter++) {
			for (int k_iter = 0; k_iter < Ncol_op; k_iter++) {
				res(p_iter, k_iter) = y_m[p_iter - p_base][k_iter];
			}
		}
		if (debug)
			cout << "here = " << seg_iter << "\t" << res.sum () << endl;
	}
}

void MatMult::multiply_y_post_fast_thread(int begin, int end, MatrixXdr &op, int Ncol_op, double *yint_e, double **y_e, double *partialsums, double *y_e_data_flat) {
	// Zero the entire contiguous y_e block in one call instead of Nindv separate memsets
	std::memset(y_e_data_flat, 0, static_cast<size_t>(g.Nindv) * blocksize * sizeof(double));
	std::memset(yint_e,     0, static_cast<size_t>(hsize) * blocksize * sizeof(double));
  	std::memset(partialsums, 0, static_cast<size_t>(blocksize) * sizeof(double));

	for (int seg_iter = begin; seg_iter < end; seg_iter++) {
		mailman::fastmultiply_pre(g.segment_size_hori, g.Nindv, Ncol_op, seg_iter * g.segment_size_hori, g.p[seg_iter], op, yint_e, partialsums, y_e);
	}
}


void MatMult::multiply_y_pre_fast(MatrixXdr &op, int Ncol_op, MatrixXdr &res, bool subtract_means) {
	for (int k_iter = 0; k_iter < Ncol_op; k_iter++) {
		sum_op[k_iter] = op.col(k_iter).sum();
	}

	if (debug) {
		// print_time();
		std::cout << "Starting mailman on premultiply" << std::endl;
		std::cout << "Nops = " << Ncol_op << "\t" << g.Nsegments_hori << std::endl;
		std::cout << "Segment size = " << g.segment_size_hori << std::endl;
		std::cout << "Matrix size = " << g.segment_size_hori << "\t" << g.Nindv << std::endl;
		std::cout << "op = " <<  op.rows() << "\t" << op.cols() << std::endl;
	}

	// TODO: Memory Effecient SSE FastMultipy

	int eff_threads = (nthreads > g.Nsegments_hori) ? g.Nsegments_hori : nthreads;
	eff_threads = (eff_threads < 1) ? 1 : eff_threads;

	int perthread = g.Nsegments_hori / eff_threads;

	if (pool && eff_threads > 1) {
		pool->parallel_for(eff_threads, [&](int t) {
			int begin = t * perthread;
			int end = (t < eff_threads - 1) ? (t + 1) * perthread : g.Nsegments_hori - 1;
			multiply_y_pre_fast_thread(begin, end, op, Ncol_op, yint_m[t], y_m[t], partialsums[t], res, y_m_data[t]);
		});
	} else {
		multiply_y_pre_fast_thread(0, g.Nsegments_hori - 1, op, Ncol_op, yint_m[0], y_m[0], partialsums[0], res, y_m_data[0]);
	}

	int last_seg_size = (g.Nsnp % g.segment_size_hori != 0) ? g.Nsnp % g.segment_size_hori : g.segment_size_hori;
	memset(partialsums[0], 0, blocksize * sizeof(double));
	memset(yint_m[0], 0, hsize * blocksize * sizeof(double));
	memset(y_m_data[0], 0, static_cast<size_t>(g.segment_size_hori) * blocksize * sizeof(double));
	mailman::fastmultiply(last_seg_size, g.Nindv, Ncol_op, g.p[g.Nsegments_hori-1], op, yint_m[0], partialsums[0], y_m[0]);
	int p_base = (g.Nsegments_hori - 1) * g.segment_size_hori;
	for (int p_iter = p_base; (p_iter < p_base + g.segment_size_hori) && (p_iter < g.Nsnp); p_iter++) {
		for (int k_iter = 0; k_iter < Ncol_op; k_iter++) {
			res(p_iter, k_iter) = y_m[0][p_iter - p_base][k_iter];
		}
	}

	#if DEBUG == 1
		if (debug) {
			// print_time();
			std::cout << "Ending mailman on premultiply" << std::endl;
		}
	#endif

	if (!subtract_means) {
		return;
	}

	for (int p_iter = 0; p_iter < g.Nsnp; p_iter++) {
		for (int k_iter = 0; k_iter < Ncol_op; k_iter++) {
			res(p_iter, k_iter) = res(p_iter, k_iter) - (g.get_col_mean(p_iter) * sum_op[k_iter]);
			if (var_normalize) {
				double s = g.get_col_std(p_iter);
				if (s > 0.0)
					res(p_iter, k_iter) /= s;
				else
					res(p_iter, k_iter) = 0.0; // degenerate column
			}
		}
	}
}


void MatMult::multiply_y_post_fast(MatrixXdr &op_orig, int Nrows_op, MatrixXdr &res, bool subtract_means) {
	MatrixXdr op;
	op = op_orig.transpose();

	// safer standardization
	if (var_normalize && subtract_means) {
		for (int p_iter = 0; p_iter < g.Nsnp; p_iter++) {
			double s = g.get_col_std(p_iter);
			if (s > 0.0) {
				for (int k_iter = 0; k_iter < Nrows_op; k_iter++) {
					op(p_iter, k_iter) /= s;
				}
			}
		}
	}

	#if DEBUG == 1
		if (debug) {
			// print_time();
			std::cout << "Starting mailman on postmultiply" << std::endl;
		}
	#endif

	int Ncol_op = Nrows_op;

	int eff_threads = (nthreads > g.Nsegments_hori) ? g.Nsegments_hori : nthreads;
	eff_threads = (eff_threads < 1) ? 1 : eff_threads;

	int perthread = g.Nsegments_hori / eff_threads;

	if (pool && eff_threads > 1) {
		pool->parallel_for(eff_threads, [&](int t) {
			int begin = t * perthread;
			int end = (t < eff_threads - 1) ? (t + 1) * perthread : g.Nsegments_hori - 1;
			multiply_y_post_fast_thread(begin, end, op, Ncol_op, yint_e[t], y_e[t], partialsums[t], y_e_data[t]);
		});
	} else {
		multiply_y_post_fast_thread(0, g.Nsegments_hori - 1, op, Ncol_op, yint_e[0], y_e[0], partialsums[0], y_e_data[0]);
	}

	// int seg_iter;
	// for(seg_iter = 0; seg_iter < g.Nsegments_hori-1; seg_iter++){
	// 	mailman::fastmultiply_pre (g.segment_size_hori, g.Nindv, Ncol_op, seg_iter * g.segment_size_hori, g.p[seg_iter], op, yint_e, partialsums[0], y_e);
	// }

	// Reduce across threads using flat contiguous data
	{
		size_t total = static_cast<size_t>(g.Nindv) * Ncol_op;
		for (int t = 1; t < eff_threads; t++) {
			double *dst = y_e_data[0];
			double *src = y_e_data[t];
			for (size_t i = 0; i < total; i++) {
				dst[i] += src[i];
			}
		}
	}

	int last_seg_size = (g.Nsnp % g.segment_size_hori != 0) ? g.Nsnp % g.segment_size_hori : g.segment_size_hori;
	mailman::fastmultiply_pre(last_seg_size, g.Nindv, Ncol_op, (g.Nsegments_hori-1) * g.segment_size_hori,
							  g.p[g.Nsegments_hori-1], op, yint_e[0], partialsums[0], y_e[0]);

	for (int n_iter = 0; n_iter < g.Nindv; n_iter++)  {
		for (int k_iter = 0; k_iter < Ncol_op; k_iter++) {
			res(k_iter, n_iter) = y_e[0][n_iter][k_iter];
			y_e[0][n_iter][k_iter] = 0;
		}
	}

	#if DEBUG == 1
		if (debug) {
			// print_time();
			std::cout << "Ending mailman on postmultiply" << std::endl;
		}
	#endif


	if (!subtract_means) {
		return;
	}

	std::vector<double> sums_elements(Ncol_op, 0.0);

	for (int k_iter = 0; k_iter < Ncol_op; k_iter++) {
		double sum_to_calc = 0.0;
		for (int p_iter = 0; p_iter < g.Nsnp; p_iter++) {
			sum_to_calc += g.get_col_mean(p_iter) * op(p_iter, k_iter);
		}
		sums_elements[k_iter] = sum_to_calc;
	}
	for (int k_iter = 0; k_iter < Ncol_op; k_iter++) {
		for (int n_iter = 0; n_iter < g.Nindv; n_iter++) {
			res(k_iter, n_iter) = res(k_iter, n_iter) - sums_elements[k_iter];
		}
	}
}

void MatMult::multiply_y_pre_naive_mem(MatrixXdr &op, int Ncol_op, MatrixXdr &res) {
	for (int p_iter = 0; p_iter < g.Nsnp; p_iter++) {
		for (int k_iter = 0; k_iter < Ncol_op; k_iter++) {
			double temp = 0;
			for (int n_iter = 0; n_iter < g.Nindv; n_iter++) {
				temp+= g.get_geno(p_iter, n_iter, var_normalize) * op(n_iter, k_iter);
			}
			res(p_iter, k_iter) = temp;
		}
	}
}

void MatMult::multiply_y_post_naive_mem(MatrixXdr &op, int Nrows_op, MatrixXdr &res) {
	for (int n_iter = 0; n_iter < g.Nindv; n_iter++) {
		for (int k_iter = 0; k_iter < Nrows_op; k_iter++) {
			double temp = 0;
			for(int p_iter = 0; p_iter < g.Nsnp; p_iter++)
				temp += op(k_iter, p_iter) * (g.get_geno(p_iter, n_iter, var_normalize));
			res(k_iter, n_iter) = temp;
		}
	}
}

void MatMult::multiply_y_post(MatrixXdr &op, int Nrows_op, MatrixXdr &res, bool subtract_means) {
	ScopedTimer timer("Xt_v");
	if (fast_mode) {
		multiply_y_post_fast(op, Nrows_op, res, subtract_means);
	} else {
		if (memory_efficient)
			multiply_y_post_naive_mem(op, Nrows_op, res);
		else
			res = op * geno_matrix;
	}
}

void MatMult::multiply_y_pre(MatrixXdr &op, int Ncol_op, MatrixXdr &res, bool subtract_means) {
	ScopedTimer timer("Xv");
	if (fast_mode) {
		multiply_y_pre_fast(op, Ncol_op, res, subtract_means);
	} else {
		if (memory_efficient) {
			multiply_y_pre_naive_mem(op, Ncol_op, res);
		} else {
			res = geno_matrix * op;
		}
	}
}

void MatMult::clean_up() {
  if (y_m) {
    for (int t = 0; t < nthreads; ++t) {
      if (y_m[t]) delete[] y_m[t]; // pointer array only
    }
    delete[] y_m; y_m = nullptr;
  }
  if (y_m_data) {
    for (int t = 0; t < nthreads; ++t) delete[] y_m_data[t]; // contiguous data
    delete[] y_m_data; y_m_data = nullptr;
  }

  if (y_e) {
    for (int t = 0; t < nthreads; ++t) {
      if (y_e[t]) delete[] y_e[t]; // pointer array only
    }
    delete[] y_e; y_e = nullptr;
  }
  if (y_e_data) {
    for (int t = 0; t < nthreads; ++t) delete[] y_e_data[t]; // contiguous data
    delete[] y_e_data; y_e_data = nullptr;
  }

  if (yint_m) {
    for (int t = 0; t < nthreads; ++t) delete[] yint_m[t];
    delete[] yint_m; yint_m = nullptr;
  }
  if (partialsums) {
    for (int t = 0; t < nthreads; ++t) delete[] partialsums[t];
    delete[] partialsums; partialsums = nullptr;
  }

  if (yint_e) { for (int t = 0; t < nthreads; ++t) delete[] yint_e[t];
                delete[] yint_e; yint_e = nullptr; }

  if (sum_op) { delete[] sum_op; sum_op = nullptr; }
  if (pool) { delete pool; pool = nullptr; }

  // reset sizes to safe state
  nthreads = 0; blocksize = 0; hsegsize = 0; hsize = 0; vsegsize = 0; vsize = 0;
}


MatMult::~MatMult() { clean_up(); }

static void move_from(MatMult& dst, MatMult& src) noexcept {
  // trivial fields
  dst.g = src.g; // or std::move(src.g) if genotype is movable
  dst.geno_matrix = std::move(src.geno_matrix);
  dst.debug = src.debug;
  dst.var_normalize = src.var_normalize;
  dst.memory_efficient = src.memory_efficient;
  dst.missing = src.missing;
  dst.fast_mode = src.fast_mode;
  dst.nthreads = src.nthreads;
  dst.blocksize = src.blocksize;
  dst.hsegsize = src.hsegsize;
  dst.hsize = src.hsize;
  dst.vsegsize = src.vsegsize;
  dst.vsize = src.vsize;

  // steal pointers
  dst.sum_op      = src.sum_op;      src.sum_op = nullptr;
  dst.partialsums = src.partialsums; src.partialsums = nullptr;
  dst.yint_m      = src.yint_m;      src.yint_m = nullptr;
  dst.y_m         = src.y_m;         src.y_m = nullptr;
  dst.y_m_data    = src.y_m_data;    src.y_m_data = nullptr;
  dst.yint_e      = src.yint_e;      src.yint_e = nullptr;
  dst.y_e         = src.y_e;         src.y_e = nullptr;
  dst.y_e_data    = src.y_e_data;    src.y_e_data = nullptr;
  dst.pool        = src.pool;        src.pool = nullptr;
}

MatMult::MatMult(MatMult&& other) noexcept { move_from(*this, other); }

MatMult& MatMult::operator=(MatMult&& other) noexcept {
  if (this != &other) {
    clean_up();
    move_from(*this, other);
  }
  return *this;
}

void MatMult::reset(genotype &xg,
		MatrixXdr &xgeno_matrix,
		bool xdebug,
		bool xvar_normalize,
		bool xmemory_efficient,
		bool xmissing,
		bool xfast_mode,
		int xnthreads,
		int xk) {
	ScopedTimer timer("matmult_setup");

	int new_hsegsize = xg.segment_size_hori;
	int new_hsize = pow(3, new_hsegsize);
	int new_Nindv = xg.Nindv;

	// Check if we can reuse existing buffers
	bool can_reuse = (y_e != nullptr) &&
	                 (xnthreads == nthreads) &&
	                 (xk == blocksize) &&
	                 (new_Nindv == Nindv_alloc) &&
	                 (new_hsegsize == hsegsize);

	g = xg;
	geno_matrix = xgeno_matrix;
	debug = xdebug;
	var_normalize = xvar_normalize;
	memory_efficient = xmemory_efficient;
	missing = xmissing;
	fast_mode = xfast_mode;

	if (can_reuse) {
		// Just zero the existing buffers
		for (int t = 0; t < nthreads; t++) {
			memset(partialsums[t], 0, blocksize * sizeof(double));
			memset(yint_m[t], 0, static_cast<size_t>(hsize) * blocksize * sizeof(double));
			memset(yint_e[t], 0, static_cast<size_t>(hsize) * blocksize * sizeof(double));
			memset(y_e_data[t], 0, static_cast<size_t>(g.Nindv) * blocksize * sizeof(double));
			memset(y_m_data[t], 0, static_cast<size_t>(hsegsize) * blocksize * sizeof(double));
		}
		memset(sum_op, 0, blocksize * sizeof(double));
		return;
	}

	// Dimensions changed — fall back to full reallocation
	clean_up();

	nthreads = xnthreads;
	blocksize = xk;
	hsegsize = new_hsegsize;
	hsize = new_hsize;
	vsegsize = xg.segment_size_ver;
	vsize = pow(3, vsegsize);
	Nindv_alloc = new_Nindv;

	if (nthreads > 1) {
		pool = new ForkJoinPool(nthreads);
	}

	partialsums = new double*[nthreads];
	yint_m = new double*[nthreads];
	for (int t = 0; t < nthreads; t++) {
		partialsums[t] = new double[blocksize]();
		yint_m[t] = new double[static_cast<size_t>(hsize)*blocksize]();
	}

	sum_op = new double[blocksize]();

	yint_e = new double*[nthreads];
	for (int t = 0; t < nthreads; t++) {
		yint_e[t] = new double[static_cast<size_t>(hsize)*blocksize]();
	}

	y_e = new double**[nthreads];
	y_e_data = new double*[nthreads];
	for (int t = 0; t < nthreads; t++) {
		y_e_data[t] = new double[static_cast<size_t>(g.Nindv) * blocksize]();
		y_e[t] = new double*[g.Nindv];
		for (int i = 0; i < g.Nindv; i++) {
			y_e[t][i] = y_e_data[t] + static_cast<size_t>(i) * blocksize;
		}
	}

	y_m = new double**[nthreads];
	y_m_data = new double*[nthreads];
	for (int t = 0; t < nthreads; t++) {
		y_m_data[t] = new double[static_cast<size_t>(hsegsize) * blocksize]();
		y_m[t] = new double*[hsegsize];
		for (int i = 0; i < hsegsize; i++) {
			y_m[t][i] = y_m_data[t] + static_cast<size_t>(i) * blocksize;
		}
	}
}
