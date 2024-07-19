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

MatrixXf E_tilde_list(6,6);

MatrixXf E_tilde_dot_list(6,6);

MatrixXf lambda_list(6,6);

MatrixXf A1(6,6);
MatrixXf A2(6,6);
MatrixXf A3(6,6);
MatrixXf A4(6,6);
MatrixXf A5(6,6);
MatrixXf A6(6,6);

Matrix4f T0_01;
Matrix4f T0_12;
Matrix4f T0_23;
Matrix4f T0_34;
Matrix4f T0_45;
Matrix4f T0_56;

Matrix4f T01;
Matrix4f T12;
Matrix4f T23;
Matrix4f T34;
Matrix4f T45;
Matrix4f T56;

Matrix4f T1J;
Matrix4f T2J;
Matrix4f T3J;
Matrix4f T4J;
Matrix4f T5J;
Matrix4f T6J;

MatrixXf Ad01(6,6);
MatrixXf Ad12(6,6);
MatrixXf Ad23(6,6);
MatrixXf Ad34(6,6);
MatrixXf Ad45(6,6);
MatrixXf Ad56(6,6);

MatrixXf Ad1J(6,6);
MatrixXf Ad2J(6,6);
MatrixXf Ad3J(6,6);
MatrixXf Ad4J(6,6);
MatrixXf Ad5J(6,6);
MatrixXf Ad6J(6,6);

MatrixXf ad01(6,6);
MatrixXf ad12(6,6);
MatrixXf ad23(6,6);
MatrixXf ad34(6,6);
MatrixXf ad45(6,6);
MatrixXf ad56(6,6);