#ifndef __RECONSTRUC_H_
#define __RECONSTRUC_H_

#include "ele.h"
#include <vector>
#include "Eigen/Dense"

void reconstruct(std::vector<Element>& ele, 
    double a, 
    double b, 
    int n, 
    int n_p, 
    int m);

void find_quadrature_info(std::vector<Element>& ele, 
    double a, 
    double b, 
    int n, 
    int acc);

void lu_decompose(Eigen::MatrixXd &mtx,
    Eigen::MatrixXd& lmtx,
    Eigen::MatrixXd& umtx,
    Eigen::MatrixXd& mtx1);

void lu_inverse(Eigen::MatrixXd& lmtx,
    Eigen::MatrixXd& umtx,
    Eigen::MatrixXd& limtx,
    Eigen::MatrixXd& uimtx);

std::vector<std::vector<double>> get_basis(
    Eigen::MatrixXd& T,
    Eigen::MatrixXd& Tt);

double L2_erro(std::vector<Element>& ele);


#endif
