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

MatrixXf getTlist()
{
    MatrixXf T(4,4);

}

VectorXf systemGravity()
{
    VectorXf G(6);
    VectorXf g_ref(6);
    MatrixXf g(6,6);
    g_ref << 0,0,0,0,0,-9.81;

    for (int i=0; i<6; i++) 
    {

    }
}