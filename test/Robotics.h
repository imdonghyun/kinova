#define _USE_MATH_DEFINES
#include <cmath>
#include <Eigen/Dense>


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

MatrixXf getTlist();

MatrixXf getAdlist();