#define _USE_MATH_DEFINES
#include <cmath>
#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>
#include "config.h"

using namespace Eigen;

bool isZero(float value);

MatrixXf vecToso3(VectorXf w);

MatrixXf vecTose3(VectorXf V);

MatrixXf Tinv(MatrixXf T);

VectorXf so3Tovec(MatrixXf so3m);

VectorXf se3Tovec(MatrixXf se3m);

MatrixXf expso3(VectorXf w);

MatrixXf dexpso3(VectorXf w);

MatrixXf expse3(VectorXf V);

VectorXf logSO3(MatrixXf R);

VectorXf logSE3(MatrixXf T);

MatrixXf Adjoint(MatrixXf T);

MatrixXf adjop(VectorXf V);

MatrixXf E_list(4,4);

MatrixXf E_tilde_list(6,6);

MatrixXf E_tilde_dot_list(6,6);

MatrixXf lambda_list(6,6);

MatrixXf T0list(24,4);

MatrixXf Tlist(24,4);

MatrixXf TJlist(24,4);

MatrixXf Adlist(36,6);

MatrixXf AdJlist(36,6);

MatrixXf adlist(36,6);

MatrixXf Alist(36,6);

MatrixXf Astarlist(36,6);
