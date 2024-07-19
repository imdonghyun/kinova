#include <iostream>
#include "Robotics.h"

bool isZero(float value)
{
    if (abs(value)<0.000001) return true;

    else return false;
}

MatrixXf Tinv(MatrixXf T)
{
    MatrixXf Tinv(4,4);
    MatrixXf Rinv(3,3);
    VectorXf p(3);

    Tinv.setIdentity();

    Rinv = T.block<3,3>(0,0).transpose();
    p << T(0,3), T(1,3), T(2,3);

    Tinv.block<3,3>(0,0) = Rinv;
    Tinv.block<3,1>(0,3) = -Rinv * p;
    return Tinv;
}

MatrixXf vecToso3(VectorXf w)
{
    float x = w(0);
    float y = w(1);
    float z = w(2);

    MatrixXf so3m(3,3);
    so3m << 0,-z, y,
            z, 0,-x,
            -y,x, 0;
    
    return so3m;
}

MatrixXf vecTose3(VectorXf V)
{
    MatrixXf se3m(4,4);
    
    se3m.block<3,3>(0,0) = vecToso3(V.head(3));
    se3m.block<3,1>(0,3) = V.tail(3);
    se3m.block<1,4>(3,0).setZero();

    return se3m;
}

VectorXf so3Tovec(MatrixXf so3m)
{
    VectorXf w(3);
    w << so3m(2,1), so3m(0,2), so3m(1,0);
    return w;
}

VectorXf se3Tovec(MatrixXf se3m)
{
    VectorXf V(6);
    V.head(3) = so3Tovec(se3m.block<3,3>(0,0));
    V.tail(3) = se3m.block<3,1>(0,3);
    return V;
}

MatrixXf expso3(VectorXf w)
{
    float theta = w.norm();
    MatrixXf R(3,3);
    R.setIdentity();

    if (!isZero(theta))
    {
        MatrixXf omg(3,3);
        omg = vecToso3(w/theta);
        R = R + sin(theta)*omg + (1-cos(theta))*omg*omg; 
    }
    return R;
}

MatrixXf dexpso3(VectorXf w)
{
    float theta = w.norm();
    MatrixXf dexp(3,3);
    dexp.setIdentity();

    if (!isZero(theta))
    {
        MatrixXf omg(3,3);
        omg = vecToso3(w/theta);
        dexp = dexp + (1-cos(theta))/theta*omg + (1-sin(theta)/theta)*omg*omg; 
    }
    return dexp;
}

MatrixXf dexpso3inv(VectorXf w)
{
    float theta = w.norm();
    MatrixXf dexpinv(3,3);
    dexpinv.setIdentity();

    if (!isZero(theta))
    {
        MatrixXf omg(3,3);
        float theta_2 = theta/2;
        omg = vecToso3(w/theta);
        dexpinv = dexpinv - theta_2*omg + (1-theta_2/tan(theta_2))*omg*omg; 
    }
    return dexpinv;
}

MatrixXf expse3(VectorXf V)
{
    MatrixXf T(4,4);
    VectorXf w = V.head(3);
    VectorXf v = V.tail(3);
    T.setIdentity();

    T.block<3,3>(0,0) = expso3(w);
    T.block<3,1>(0,3) = dexpso3(w)*v;
    return T;
}

VectorXf logSO3(MatrixXf R)
{
    float cos = (R.trace()-1)/2;
    VectorXf w(3);

    if (isZero(cos-1))
    {
        w << 0,0,0;
    }
    else if (isZero(cos+1))
    {
        if (!isZero(1+R(2,2)))
        {
            w << R(0,2), R(1,2), 1+R(2,2);
            w = w*M_PI/sqrt(2*(1+R(2,2)));
        }
        else if (!isZero(1+R(1,1)))
        {
            w << R(0,1), 1+R(1,1), R(2,1);
            w = w*M_PI/sqrt(2*(1+R(1,1)));
        }
        else
        {
            w << 1+R(0,0), R(1,0), R(2,0);
            w = w*M_PI/sqrt(2*(1+R(0,0)));
        }
    }
    else
    {
        float theta = acos(cos);
        w = so3Tovec(R-R.transpose())*theta/(2*sin(theta));
    }
    return w;
}

VectorXf logSE3(MatrixXf T)
{
    VectorXf V(6);
    VectorXf w(3);
    VectorXf v(3);
    VectorXf p(3);

    p << T(0,3), T(1,3), T(2,3);
    w = logSO3(T.block<3,3>(0,0));
    v = dexpso3inv(w) * p;
    V.head(3) = w;
    V.tail(3) = v;

    return V;
}

MatrixXf Adjoint(MatrixXf T)
{
    MatrixXf AdT(6,6);
    MatrixXf R(3,3);
    VectorXf p(3);

    
    R = T.block<3,3>(0,0);
    p << T(0,3), T(1,3), T(2,3);

    AdT.block<3,3>(0,0) = R;
    AdT.block<3,3>(0,3).setZero();
    AdT.block<3,3>(3,0) = vecToso3(p)*R;
    AdT.block<3,3>(3,3) = R;

    return AdT;
}

MatrixXf adjop(VectorXf V)
{
    MatrixXf adV(6,6);
    VectorXf w(3);
    VectorXf v(3);

    w = V.head(3);
    v = V.tail(3);

    adV.block<3,3>(0,0) = vecToso3(w);
    adV.block<3,3>(0,3).setZero();
    adV.block<3,3>(3,0) = vecToso3(v);
    adV.block<3,3>(3,3) = vecToso3(w);

    return adV;
}

void setConstant()
{
    VectorXf v1(6);
    VectorXf v2(6);
    VectorXf v3(6);
    VectorXf v4(6);
    VectorXf v5(6);
    VectorXf v6(6);

    v1 << 0,0,1,0,0,0;
    v2 << 0,0,1,-0.03,0,0;
    v3 << 0,0,1,0.28,0,0;
    v4 << 0,0,1,-0.14,0,0;
    v5 << 0,0,1,0,-0.285,0;
    v6 << 0,0,1,0,0.105,0;

    lambda_list << v1, v2, v3, v4, v5, v6;

    T0_01 << 1, 0, 0, CM[0][0]+link[0][0],
             0, 1, 0, CM[0][1]+link[0][1],
             0, 0, 1, CM[0][2]+link[0][2],
             0, 0, 0, 1;

    T0_12 << 1, 0, 0, CM[1][0]+link[1][0],
             0, 0,-1, CM[1][1]+link[1][1],
             0, 1, 0, CM[1][2]+link[1][2],
             0, 0, 0, 1;

    T0_23 << 1, 0, 0, CM[2][0]+link[2][0],
             0,-1, 0, CM[2][1]+link[2][1],
             0, 0,-1, CM[2][2]+link[2][2],
             0, 0, 0, 1;

    T0_34 << 1, 0, 0, CM[3][0]+link[3][0],
             0, 0,-1, CM[3][1]+link[3][1],
             0, 1, 0, CM[3][2]+link[3][2],
             0, 0, 0, 1;

    T0_45 << 0, 0, 1, CM[4][0]+link[4][0],
             0, 1, 0, CM[4][1]+link[4][1],
             -1,0, 0, CM[4][2]+link[4][2],
             0, 0, 0, 1;

    T0_56 << 0, 0,-1, CM[5][0]+link[5][0],
             0, 1, 0, CM[5][1]+link[5][1],
             1, 0, 0, CM[5][2]+link[5][2],
             0, 0, 0, 1;

    T1J << 1, 0, 0, -CM[0][0],
           0, 1, 0, -CM[0][1],
           0, 0, 1, -CM[0][2],
           0, 0, 0, 1;

    T2J << 1, 0, 0, -CM[1][0],
           0, 1, 0, -CM[1][1],
           0, 0, 1, -CM[1][2],
           0, 0, 0, 1;

    T3J << 1, 0, 0, -CM[2][0],
           0, 1, 0, -CM[2][1],
           0, 0, 1, -CM[2][2],
           0, 0, 0, 1;

    T4J << 1, 0, 0, -CM[3][0],
           0, 1, 0, -CM[3][1],
           0, 0, 1, -CM[3][2],
           0, 0, 0, 1;

    T5J << 1, 0, 0, -CM[4][0],
           0, 1, 0, -CM[4][1],
           0, 0, 1, -CM[4][2],
           0, 0, 0, 1;

    T6J << 1, 0, 0, -CM[5][0],
           0, 1, 0, -CM[5][1],
           0, 0, 1, -CM[5][2],
           0, 0, 0, 1;

    Ad1J = Adjoint(T1J);
    Ad2J = Adjoint(T2J);
    Ad3J = Adjoint(T3J);
    Ad4J = Adjoint(T4J);
    Ad5J = Adjoint(T5J);
    Ad6J = Adjoint(T6J);
    
    VectorXf Ei(6);
    Ei << 0,0,1,0,0,0;
    
    v1 = Ad1J*Ei;
    v2 = Ad2J*Ei;
    v3 = Ad3J*Ei;
    v4 = Ad4J*Ei;
    v5 = Ad5J*Ei;
    v6 = Ad6J*Ei;
    E_tilde_list << v1, v2, v3, v4, v5, v6;
    E_tilde_dot_list << 0,0,0,0,0,0,
                        0,0,0,0,0,0,
                        0,0,0,0,0,0,
                        0,0,0,0,0,0,
                        0,0,0,0,0,0,
                        0,0,0,0,0,0;

    A1 << Inertia[0][0], Inertia[0][1], Inertia[0][2], 0, 0, 0,
          Inertia[0][1], Inertia[0][3], Inertia[0][4], 0, 0, 0,
          Inertia[0][2], Inertia[0][4], Inertia[0][5], 0, 0, 0,
          0, 0, 0, Mass[0], 0, 0,
          0, 0, 0, 0, Mass[0], 0,
          0, 0, 0, 0, 0, Mass[0];

    A2 << Inertia[1][0], Inertia[1][1], Inertia[1][2], 0, 0, 0,
          Inertia[1][1], Inertia[1][3], Inertia[1][4], 0, 0, 0,
          Inertia[1][2], Inertia[1][4], Inertia[1][5], 0, 0, 0,
          0, 0, 0, Mass[1], 0, 0,
          0, 0, 0, 0, Mass[1], 0,
          0, 0, 0, 0, 0, Mass[1];
    
    A3 << Inertia[2][0], Inertia[2][1], Inertia[2][2], 0, 0, 0,
          Inertia[2][1], Inertia[2][3], Inertia[2][4], 0, 0, 0,
          Inertia[2][2], Inertia[2][4], Inertia[2][5], 0, 0, 0,
          0, 0, 0, Mass[2], 0, 0,
          0, 0, 0, 0, Mass[2], 0,
          0, 0, 0, 0, 0, Mass[2];

    A4 << Inertia[3][0], Inertia[3][1], Inertia[3][2], 0, 0, 0,
          Inertia[3][1], Inertia[3][3], Inertia[3][4], 0, 0, 0,
          Inertia[3][2], Inertia[3][4], Inertia[3][5], 0, 0, 0,
          0, 0, 0, Mass[3], 0, 0,
          0, 0, 0, 0, Mass[3], 0,
          0, 0, 0, 0, 0, Mass[3];

    A5 << Inertia[4][0], Inertia[4][1], Inertia[4][2], 0, 0, 0,
          Inertia[4][1], Inertia[4][3], Inertia[4][4], 0, 0, 0,
          Inertia[4][2], Inertia[4][4], Inertia[4][5], 0, 0, 0,
          0, 0, 0, Mass[4], 0, 0,
          0, 0, 0, 0, Mass[4], 0,
          0, 0, 0, 0, 0, Mass[4];

    A6 << Inertia[5][0], Inertia[5][1], Inertia[5][2], 0, 0, 0,
          Inertia[5][1], Inertia[5][3], Inertia[5][4], 0, 0, 0,
          Inertia[5][2], Inertia[5][4], Inertia[5][5], 0, 0, 0,
          0, 0, 0, Mass[5], 0, 0,
          0, 0, 0, 0, Mass[5], 0,
          0, 0, 0, 0, 0, Mass[5];
}

void updateTlist(VectorXf theta)
{
    T01 = expse3(lambda_list.col(0)*theta(0))*T0_01;
    T12 = expse3(lambda_list.col(1)*theta(1))*T0_12;
    T23 = expse3(lambda_list.col(2)*theta(2))*T0_23;
    T34 = expse3(lambda_list.col(3)*theta(3))*T0_34;
    T45 = expse3(lambda_list.col(4)*theta(4))*T0_45;
    T56 = expse3(lambda_list.col(5)*theta(5))*T0_56;
}

void updateAdlist()
{
    Ad01 = Adjoint(T01);
    Ad12 = Adjoint(T12);
    Ad23 = Adjoint(T23);
    Ad34 = Adjoint(T34);
    Ad45 = Adjoint(T45);
    Ad56 = Adjoint(T56);
}

void updateadlist(VectorXf dtheta)
{
    ad01 = adjop(E_tilde_list.col(0)*dtheta(0));
    ad12 = adjop(E_tilde_list.col(1)*dtheta(1));
    ad23 = adjop(E_tilde_list.col(2)*dtheta(2));
    ad34 = adjop(E_tilde_list.col(3)*dtheta(3));
    ad45 = adjop(E_tilde_list.col(4)*dtheta(4));
    ad56 = adjop(E_tilde_list.col(5)*dtheta(5));
}



// VectorXf systemGravity()
// {
//     VectorXf G(6);
//     VectorXf g_ref(6);
//     MatrixXf g(6,6);
//     g_ref << 0,0,0,0,0,-9.81;

//     for (int i=0; i<6; i++) 
//     {

//     }
// }
