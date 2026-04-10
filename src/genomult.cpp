#include "genomult.h"
#include "profiler.h"

#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

// Note: multiply_y_pre/post take non-const MatrixXdr& but don't actually mutate op.
// We use const_cast to avoid unnecessary copies when passing const refs through.

//
// Compute y^T X X^T y : output is a scalar
// X = X0[1:Nindv,index:(index+num_snp-1)]
// X0 : genotype matrix of Nindv X Nsnp
// vec (matrix of dimension Nindv X 1) (usually phenotype)
double compute_yXXy (int num_snp, const MatrixXdr &vec){
	ScopedTimer timer("yXXy");
	MatrixXdr res = MatrixXdr::Zero (num_snp, 1);

	if (verbose >= 3){
		cout << "***In compute_yXXy***" << endl;
		cout << "res = " << res.rows() << "," << res.cols () << "\t" << res.sum()<<endl;
		cout << "means = " << means.rows() << "," << means.cols () << "\t" << means.sum()<<endl;
		cout << "stds = " << stds.rows() << "," << stds.cols () << "\t" << stds.sum()<<endl;
	}

	mm.multiply_y_pre (const_cast<MatrixXdr&>(vec), 1, res, false);

	if (verbose >= 4)
		cout << "res = " << res.rows() << "," << res.cols () << "\t" << res.sum()<<endl;

	res = res.cwiseProduct(stds);
	MatrixXdr resid(num_snp, 1);
	resid.noalias() = means.cwiseProduct(stds);
	resid *= vec.sum();

	if (verbose >= 4)
		cout << "resid = " << resid.rows() << "," << resid.cols () << "\t" << resid.sum()<<endl;

	res -= resid;  // Xy = res - resid, reuse res

	if (verbose >= 4)
		cout << "Xy = " << res.rows() << "," << res.cols () << "\t" << res.sum()<<endl;

	double yXXy = (res.array() * res.array()).sum();
	return yXXy;
}


double compute_yVXXVy(int num_snp){
	ScopedTimer timer("yXXy");
	MatrixXdr new_pheno_sum = new_pheno.colwise().sum();
	MatrixXdr res(num_snp, 1);

	mm.multiply_y_pre(new_pheno, 1, res, false);

	res = res.cwiseProduct(stds);
	MatrixXdr resid(num_snp, 1);
	resid.noalias() = means.cwiseProduct(stds);
	resid *= new_pheno_sum(0,0);
	res -= resid;
	double ytVXXVy = (res.array() * res.array()).sum();
	return ytVXXVy;
}

// Compute X X^T Z : Nindv X Nz matrix
// X = X0[1:Nindv,index:(index+num_snp-1)]
// X0 : genotype matrix of Nindv X Nsnp
// Z  : Zvec (matrix of dimension Nindv X Nz) (usually random vectors)
MatrixXdr  compute_XXz (int num_snp, const MatrixXdr &Zvec){
	ScopedTimer timer("XXz");
	MatrixXdr res = MatrixXdr::Zero (num_snp, Nz);

	if (verbose >= 3) {
		cout << "***In compute_XXz***" << endl;
		cout << "res [" << res.rows() << "," << res.cols () << "]\t" << res.sum()<<endl;
		cout << "Zvec [" << Zvec.rows() << "," << Zvec.cols () << "]\t" << Zvec.sum()<<endl;
	}

	mm.multiply_y_pre(const_cast<MatrixXdr&>(Zvec), Nz, res, false);

	if (verbose >= 4) {
		cout << "res [ " << res.rows() << ", " << res.cols () << "]\t" << res.sum()<<endl;
		cout << "means [ " << means.rows() << "," << means.cols () << "]\t" << means.sum()<<endl;
		cout << "stds [ " << stds.rows() << "," << stds.cols () << "]\t" << stds.sum()<<endl;
	}

	MatrixXdr zb_sum = Zvec.colwise().sum();

	for(int j = 0; j < num_snp; j++)
		for(int k = 0; k < Nz ; k++)
			res(j,k) = res(j,k) * stds(j,0);

	if (verbose >= 4)
		cout << "res [" << res.rows() << "," << res.cols () << "]\t" << res.sum()<<endl;

	// inter_zb = res - means.*stds * zb_sum, then scale by stds again
	// Fuse: compute inter = means.*stds once, subtract, scale
	MatrixXdr inter = means.cwiseProduct(stds);
	// res -= inter * zb_sum  (res becomes inter_zb)
	res.noalias() -= inter * zb_sum;

	if (verbose >= 4) {
		cout << "inter_zb [" << res.rows() << "," << res.cols () << "]\t" << res.sum()<<endl;
	}

	for(int k = 0; k < Nz; k++)
		for(int j = 0; j < num_snp ; j++)
			res(j,k) = res(j,k) * stds(j,0);

	// new_zb = inter_zb^T, reuse res as inter_zb
	MatrixXdr new_zb = res.transpose();
	MatrixXdr new_res(Nz, Nindv);

	mm.multiply_y_post(new_zb, Nz, new_res, false);
	if (verbose >= 4)
		cout << "new_res = " << new_res.rows() << "," << new_res.cols () << "\t" << new_res.sum()<<endl;

	// new_resid = (new_zb * means) * ones(1, Nindv)
	MatrixXdr zb_scale_sum;
	zb_scale_sum.noalias() = new_zb * means;
	// Subtract: new_res -= zb_scale_sum * ones
	for (int i = 0; i < Nz; i++) {
		double val = zb_scale_sum(i, 0);
		for (int j = 0; j < Nindv; j++) {
			new_res(i, j) = (new_res(i, j) - val) * mask(j, 0);
		}
	}

	if (verbose >= 3)
		cout << "temp = " << new_res.rows() << "," << new_res.cols () << "\t" << new_res.sum()<<endl;

	return new_res.transpose();
}

MatrixXdr  compute_XXUz (int num_snp){
	ScopedTimer timer("XXz");
	res.resize(num_snp, Nz);

	mm.multiply_y_pre(all_Uzb,Nz,res, false);

	MatrixXdr zb_sum = all_Uzb.colwise().sum();

	for(int j = 0; j < num_snp; j++)
		for(int k = 0; k < Nz ; k++)
			res(j,k) = res(j,k) * stds(j,0);

	MatrixXdr inter = means.cwiseProduct(stds);
	res.noalias() -= inter * zb_sum;

	for(int k = 0; k < Nz; k++)
		for(int j = 0; j < num_snp ; j++)
			res(j,k) = res(j,k) * stds(j,0);

	MatrixXdr new_zb = res.transpose();
	MatrixXdr new_res(Nz, Nindv);

	mm.multiply_y_post(new_zb, Nz, new_res, false);

	MatrixXdr zb_scale_sum;
	zb_scale_sum.noalias() = new_zb * means;
	for (int i = 0; i < Nz; i++) {
		double val = zb_scale_sum(i, 0);
		for (int j = 0; j < Nindv; j++) {
			new_res(i, j) = (new_res(i, j) - val) * mask(j, 0);
		}
	}

	return new_res.transpose();
}


MatrixXdr compute_yXXy_multi (int num_snp, const MatrixXdr &vec, int cur_pheno_count){
	ScopedTimer timer("yXXy");

	MatrixXdr pheno_sum=vec.colwise().sum();
	MatrixXdr res = MatrixXdr::Zero (num_snp, cur_pheno_count);

	if (verbose >= 3){
		cout << "***In compute_yXXy_multi***" << endl;
	}

	mm.multiply_y_pre (const_cast<MatrixXdr&>(vec), cur_pheno_count, res, false);

	for(int j=0; j<num_snp; j++)
		for(int k=0; k<cur_pheno_count;k++)
				res(j,k) = res(j,k)*stds(j,0);

	if (verbose >= 4)
		cout << "res = " << res.rows() << "," << res.cols () << "\t" << res.sum()<<endl;

	MatrixXdr resid(num_snp, cur_pheno_count);
	for (int j=0; j<num_snp; j++)
		for (int k=0; k<cur_pheno_count;k++)
			resid(j, k) = means(j, k)*stds(j, k);
	resid *= vec.sum();

	if (verbose >= 4)
		cout << "resid = " << resid.rows() << "," << resid.cols () << "\t" << resid.sum()<<endl;

	res -= resid;

	if (verbose >= 4)
		cout << "Xy = " << res.rows() << "," << res.cols () << "\t" << res.sum()<<endl;

	res = res.array() * res.array();
    MatrixXdr out_temp = res.colwise().sum();
	return out_temp;
}


MatrixXdr compute_yVXXVy_multi(int num_snp, const MatrixXdr &vec, int cur_pheno_count){
	ScopedTimer timer("yXXy");
	MatrixXdr new_pheno_sum = vec.colwise().sum();
	MatrixXdr res(num_snp, cur_pheno_count);

	mm.multiply_y_pre(new_pheno, cur_pheno_count, res, false);

	for(int j=0; j<num_snp; j++)
		for(int k=0; k<cur_pheno_count;k++)
				res(j,k) = res(j,k)*stds(j,0);

	MatrixXdr resid(num_snp, cur_pheno_count);

	for (int j=0; j<num_snp; j++)
		for (int k=0; k<cur_pheno_count;k++)
			resid(j, k) = means(j, k)*stds(j, k);

	resid *= new_pheno_sum(0,0);

	res -= resid;
	res = res.array() * res.array();
    MatrixXdr out_temp = res.colwise().sum();
    return out_temp;
}
