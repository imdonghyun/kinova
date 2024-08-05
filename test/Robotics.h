#define _USE_MATH_DEFINES
#include <cmath>
#include <Eigen/Dense>
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

MatrixXf InvAd(MatrixXf Ad);

MatrixXf InvTransAd(MatrixXf Ad);

MatrixXf adjop(VectorXf V);

MatrixXf getInertia(int i);

void setConstant();

void updateFKList(VectorXf theta, VectorXf dtheta);

MatrixXf systemInertia();

MatrixXf systemBias(VectorXf dq);

VectorXf systemGravity();

void RNE(MatrixXf Fextlist, VectorXf V0, VectorXf dV0, VectorXf q, VectorXf dq, VectorXf ddq);

