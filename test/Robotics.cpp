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

MatrixXf InvAd(MatrixXf Ad)
{
    MatrixXf tmp(6,6);

    tmp.setZero();

    tmp.block<3,3>(0,0) = Ad.block<3,3>(0,0).transpose();
    tmp.block<3,3>(3,0) = Ad.block<3,3>(3,0).transpose();
    tmp.block<3,3>(3,3) = Ad.block<3,3>(3,3).transpose();

    return tmp;
}

MatrixXf InvTransAd(MatrixXf Ad)
{
    MatrixXf tmp(6,6);

    tmp = Ad;

    tmp.block<3,3>(3,0) = Ad.block<3,3>(0,3);
    tmp.block<3,3>(0,3) = Ad.block<3,3>(3,0);

    return tmp;
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

    T0list << 1, 0, 0, 0.000025,
              0, 1, 0, 0.022135,
              0, 0, 1, 0.227677,
              0, 0, 0, 1, //{0}->{1}

              1, 0, 0, 0.029983,
              0, 0,-1,-0.075303,
              0, 1, 0, 0.484850,
              0, 0, 0, 1, //{0}->{2}

              1, 0, 0, 0.030156,
              0, 0, 1,-0.022644,
              0,-1, 0, 0.618322,
              0, 0, 0, 1, //{0}->{3}

              1, 0, 0, 0.005752,
              0, 1, 0, 0.000004,
              0, 0, 1, 0.750492,
              0, 0, 0, 1, //{0}->{4}

              0, 0, 1, 0.047228,
              0, 1, 0,-0.000196,
              -1,0, 0, 0.687735,
              0, 0, 0, 1, //{0}->{5}

              1, 0, 0, 0.066930,
              0, 1, 0,-0.000050,
              0, 0, 1, 0.934660,
              0, 0, 0, 1; //{0}->{6}

    TJlist << 1, 0, 0, -CM[0][0],
              0, 1, 0, -CM[0][1],
              0, 0, 1, -CM[0][2],
              0, 0, 0, 1, //{1}->{1J0}

              1, 0, 0, -CM[1][0],
              0, 1, 0, -CM[1][1],
              0, 0, 1, -CM[1][2],
              0, 0, 0, 1, //{2}->{2J1}

              1, 0, 0, -CM[2][0],
              0, 1, 0, -CM[2][1],
              0, 0, 1, -CM[2][2],
              0, 0, 0, 1, //{3}->{3J2}

              1, 0, 0, -CM[3][0],
              0, 1, 0, -CM[3][1],
              0, 0, 1, -CM[3][2],
              0, 0, 0, 1, //{4}->{4J3}

              1, 0, 0, -CM[4][0],
              0, 1, 0, -CM[4][1],
              0, 0, 1, -CM[4][2],
              0, 0, 0, 1, //{5}->{5J4}

              1, 0, 0, -CM[5][0],
              0, 1, 0, -CM[5][1],
              0, 0, 1, -CM[5][2],
              0, 0, 0, 1; //{6}->{6J5}

    MatrixXf Adtmp(6,6);
    VectorXf Ei(6);
    Ei << 0,0,1,0,0,0;

    for (int i=0; i<6; i++)
    {
        Adtmp = Adjoint(TJlist.block<4,4>(i*4, 0));

        AdJlist.block<6,6>(i*6 ,0) = Adtmp;

        E_tilde_list.col(i) = Adtmp*Ei;
    }
    
    E_tilde_dot_list.setZero();

    Alist << Inertia[0][0], Inertia[0][1], Inertia[0][2], 0, 0, 0,
             Inertia[0][1], Inertia[0][3], Inertia[0][4], 0, 0, 0,
             Inertia[0][2], Inertia[0][4], Inertia[0][5], 0, 0, 0,
             0, 0, 0, Mass[0], 0, 0,
             0, 0, 0, 0, Mass[0], 0,
             0, 0, 0, 0, 0, Mass[0],

             Inertia[1][0], Inertia[1][1], Inertia[1][2], 0, 0, 0,
             Inertia[1][1], Inertia[1][3], Inertia[1][4], 0, 0, 0,
             Inertia[1][2], Inertia[1][4], Inertia[1][5], 0, 0, 0,
             0, 0, 0, Mass[1], 0, 0,
             0, 0, 0, 0, Mass[1], 0,
             0, 0, 0, 0, 0, Mass[1];

             Inertia[2][0], Inertia[2][1], Inertia[2][2], 0, 0, 0,
             Inertia[2][1], Inertia[2][3], Inertia[2][4], 0, 0, 0,
             Inertia[2][2], Inertia[2][4], Inertia[2][5], 0, 0, 0,
             0, 0, 0, Mass[2], 0, 0,
             0, 0, 0, 0, Mass[2], 0,
             0, 0, 0, 0, 0, Mass[2];
             
             Inertia[3][0], Inertia[3][1], Inertia[3][2], 0, 0, 0,
             Inertia[3][1], Inertia[3][3], Inertia[3][4], 0, 0, 0,
             Inertia[3][2], Inertia[3][4], Inertia[3][5], 0, 0, 0,
             0, 0, 0, Mass[3], 0, 0,
             0, 0, 0, 0, Mass[3], 0,
             0, 0, 0, 0, 0, Mass[3];
             
             Inertia[4][0], Inertia[4][1], Inertia[4][2], 0, 0, 0,
             Inertia[4][1], Inertia[4][3], Inertia[4][4], 0, 0, 0,
             Inertia[4][2], Inertia[4][4], Inertia[4][5], 0, 0, 0,
             0, 0, 0, Mass[4], 0, 0,
             0, 0, 0, 0, Mass[4], 0,
             0, 0, 0, 0, 0, Mass[4];

             Inertia[5][0], Inertia[5][1], Inertia[5][2], 0, 0, 0,
             Inertia[5][1], Inertia[5][3], Inertia[5][4], 0, 0, 0,
             Inertia[5][2], Inertia[5][4], Inertia[5][5], 0, 0, 0,
             0, 0, 0, Mass[5], 0, 0,
             0, 0, 0, 0, Mass[5], 0,
             0, 0, 0, 0, 0, Mass[5];

}

void updateFKList(VectorXf theta, VectorXf dtheta)
{
    Matrix4f Ttmp;

    Ttmp.setIdentity();

    for (int i=0; i<6; i++)
    {
        if (i==0)
        {
            Ttmp = expse3(lambda_list.col(i)*theta(i))*T0list.block<4,4>(i*4, 0);
        }
        else
        {
            Ttmp = Tinv(T0list.block<4,4>((i-1)*4, 0))*expse3(lambda_list.col(i)*theta(i))*T0list.block<4,4>(i*4, 0);
        }

        Tlist.block<4,4>(i*4, 0) = Ttmp;

        Adlist.block<6,6>(i*6, 0) = Adjoint(Ttmp);

        adlist.block<6,6>(i*6, 0) = adjop(E_tilde_list.col(i)*dtheta(i));
    }
}

MatrixXf systemInertia()
{
    MatrixXf M(6,6);
    MatrixXf Ad(6,6);

    Astarlist.block<6,6>(30,0) = Alist.block<6,6>(30,0);

    for (int i=5; i>0; i--)
    {
        Ad = Adlist.block<6,6>(i*6, 0);
        Astarlist.block<6,6>((i-1)*6,0) = Alist.block<6,6>((i-1)*6,0) + InvTransAd(Ad)*Astarlist.block<6,6>(i*6,0)*InvAd(Ad);
    }

    for (int i=0; i<6; i++)
    {
        Ad.setIdentity();
        for (int j=i; j<6; j++)
        {
            if (j>=i+1) Ad = Ad * Adlist.block<6,6>(j*6, 0);
            M(i,j) = E_tilde_list.col(i).transpose() * InvTransAd(Ad) * Astarlist.block<6,6>(j*6,0) * E_tilde_list.col(j);

            if (i != j)
            {
                M(j,i) = M(i,j);
            }
        }
    }
    return M;
}

VectorXf systemGravity()
{
    VectorXf G(6);
    VectorXf gref(6);
    VectorXf gtmp(6);
    MatrixXf Ad(6,6);
    MatrixXf g(6,6);

    float tmp;

    gref << 0,0,0,0,0,-9.81;

    for (int i=0; i<6; i++)
    {
        Ad = Adlist.block<6,6>(i*6, 0);

        if (i==0) gtmp = InvAd(Ad) * gref;
        else gtmp = InvAd(Ad) * gtmp;
        g.col(i) = gtmp;
    }

    for (int i=0; i<6; i++) 
    {
        tmp = 0;
        Ad.setIdentity();
        for (int j=i; j<6; j++)
        {
            if (j>=i+1) Ad = Ad * Adlist.block<6,6>(j*6, 0);
            
            tmp = tmp + E_tilde_list.col(i).transpose() * InvTransAd(Ad) * Alist.block<6,6>(j*6, 0) * g.col(j);
        }
        G(i) = -tmp;
    }
    return G;
}
