#define _USE_MATH_DEFINES
#include <iostream>
#include "Robotics.h"

/*
R: rotation matrix
T: homogeneous transform matrix
w: angular velocity
v: linear velocity
V: body twist

{0}: base frame
idx: i=0~5
{n}: nth link CoM frame
{nJn-1}: nth joint frame
T0list: transform matrix from {0} to {i+1} at initial state
TJlist: transform matrix from {i+1} to {(i+1)Ji}
Tlist: transform matrix from {i} to {i+1} after rotation
M0: initial endeffector frame from {0}

E_list: joint matrix about joint
E_tilde_list: joint matrix about CoM
lambda_list: screw from {0} to {(i+1)Ji}
Adlist: Adjoint matrix from {i} to {i+1} 
AdJlist: Adjoint matrix from {i+1} to {(i+1)Ji}
adlist: adjoint operator from {i} to {i+1} 

idx: i:0~6
Alist: ith link CoM inertia matrix
Astarlist: accumulated inertia matrix
*/

int n=6; //dof

MatrixXf E_tilde_list(6,n);

MatrixXf E_tilde_dot_list(6,n);

MatrixXf lambda_list(6,n);

MatrixXf Teef0(4,4);

MatrixXf T0list(4*n,4);

MatrixXf Tlist(4*n,4);

MatrixXf TJlist(4*n,4);

MatrixXf Adlist(6*n,6);

MatrixXf AdJlist(6*n,6);

MatrixXf adlist(6*n,6);

MatrixXf Alist(6*(n+1),6);

VectorXf tau(n);

VectorXf Fc(n);

MatrixXf M(n,n);
MatrixXf C(n,n);
VectorXf c(n);
VectorXf G(n);


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
    p = T.block<3,1>(0,3);

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
    MatrixXf so3m(3,3);

    if (isZero(cos-1))
    {
        so3m << 0,0,0,
                0,0,0,
                0,0,0;
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
        so3m = vecToso3(w);
    }
    else
    {
        float theta = acos(cos);
        so3m = R-R.transpose()*theta/(2*sin(theta));
    }
    return so3m;
}

VectorXf logSE3(MatrixXf T)
{
    VectorXf w(3);
    VectorXf p(3);
    MatrixXf so3m(3,3);
    MatrixXf se3m(4,4);

    p = T.block<3,1>(0,3);
    so3m = logSO3(T.block<3,3>(0,0));
    w = so3Tovec(so3m);
    se3m.block<3,3>(0,0) = so3m;
    se3m.block<3,1>(0,3) = dexpso3inv(w) * p;

    return se3m;
}

MatrixXf Adjoint(MatrixXf T)
{
    MatrixXf AdT(6,6);
    MatrixXf R(3,3);
    VectorXf p(3);

    
    R = T.block<3,3>(0,0);
    p = T.block<3,1>(0,3);

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

MatrixXf getInertia(int i)
{
    Matrix3f Ib;

    Ib << Inertia[i][0], Inertia[i][1], Inertia[i][2],
          Inertia[i][1], Inertia[i][3], Inertia[i][4],
          Inertia[i][2], Inertia[i][4], Inertia[i][5];

    return Ib;
}

void setConstant()
{
    VectorXf v1(6);
    VectorXf v2(6);
    VectorXf v3(6);
    VectorXf v4(6);
    VectorXf v5(6);
    VectorXf v6(6);

    for (int i=0; i<n; i++)
    {
        for (int j=0; j<6; j++)
        {
            lambda_list(j,i) = lambdalist[i][j];
        }
    }

    for (int i=0; i<4; i++)
    {
        for (int j=0; j<4; j++)
        {
            Teef0(i,j) = M0[i][j];
        }
    }

    Matrix4f T;
    Matrix4f Ttmp;
    Vector4f CMtmp;
    T.setIdentity();

    for (int i=0; i<n; i++)
    {
        for (int j=0; j<4; j++)
        {
            for (int k=0; k<4; k++)
            {
                Ttmp(j,k) = tf_list[i][j][k];
            }
        }
        T = T*Ttmp;
        CMtmp << CM[i+1][0], CM[i+1][1], CM[i+1][2], 1;

        T0list.block<3,3>(i*4, 0) = T.block<3,3>(0,0);
        T0list.block<1,3>(i*4+3, 0).setZero();
        T0list.block<4,1>(i*4, 3) = T*CMtmp;

        TJlist.block<4,4>(i*4, 0).setIdentity();
        TJlist.block<3,1>(i*4, 3) << -CM[i+1][0], -CM[i+1][1], -CM[i+1][2];
    }

    MatrixXf tmp(6,6);
    VectorXf Ei(6);
    Ei << 0,0,1,0,0,0;

    for (int i=0; i<n; i++)
    {
        tmp = Adjoint(TJlist.block<4,4>(i*4, 0));

        AdJlist.block<6,6>(i*6 ,0) = tmp;

        E_tilde_list.col(i) = tmp*Ei;
    }

    E_tilde_dot_list.setZero();

    for (int i=0; i<=n; i++)
    {
        tmp.setZero();
        tmp.topLeftCorner(3,3) = getInertia(i);
        tmp.bottomRightCorner(3,3) = MatrixXf::Identity(3,3) * Mass[i];
        Alist.block<6,6>(i*6, 0) = tmp;
    }
}

void updateFKList(VectorXf q, VectorXf dq)
{
    Matrix4f Ttmp;

    for (int i=0; i<n; i++)
    {
        if (i==0)
        {
            Ttmp = expse3(lambda_list.col(i)*q(i))*T0list.block<4,4>(i*4, 0);
        }
        else
        {
            Ttmp = Tinv(T0list.block<4,4>((i-1)*4, 0))*expse3(lambda_list.col(i)*q(i))*T0list.block<4,4>(i*4, 0);
        }

        Tlist.block<4,4>(i*4, 0) = Ttmp;

        Adlist.block<6,6>(i*6, 0) = Adjoint(Ttmp);

        adlist.block<6,6>(i*6, 0) = adjop(E_tilde_list.col(i)*dq(i));
    }
}

MatrixXf systemInertia()
{
    MatrixXf M(n,n);
    MatrixXf Ad(6,6);
    MatrixXf Astarlist(6*(n+1),6);

    Astarlist.block<6,6>(6*n,0) = Alist.block<6,6>(6*n,0);

    for (int i=n; i>0; i--)
    {
        Ad = Adlist.block<6,6>((i-1)*6, 0);
        Astarlist.block<6,6>((i-1)*6,0) = Alist.block<6,6>((i-1)*6,0) + InvTransAd(Ad)*Astarlist.block<6,6>(i*6,0)*InvAd(Ad);
    }

    for (int i=0; i<n; i++)
    {
        Ad.setIdentity();
        for (int j=i; j<n; j++)
        {
            if (j>=i+1) Ad = Ad * Adlist.block<6,6>(j*6, 0);
            M(i,j) = E_tilde_list.col(i).transpose() * InvTransAd(Ad) * Astarlist.block<6,6>((j+1)*6,0) * E_tilde_list.col(j);

            if (i != j)
            {
                M(j,i) = M(i,j);
            }
        }
    }
    return M;
}

MatrixXf systemBias(VectorXf dq)
{
    MatrixXf C(n,n);
    MatrixXf ad(6,6);
    MatrixXf Blist(6*(n+1),6);
    MatrixXf Bstarlist(6*(n+1),6);
    MatrixXf Astarlist(6*(n+1),6);
    MatrixXf Ad(6,6);
    MatrixXf Astar(6,6);
    MatrixXf Bstar(6,6);
    MatrixXf tmp(6,6);
    VectorXf V(6);
    VectorXf Ei(6);
    VectorXf Ej(6);
    V.setZero();
    Blist.block<6,6>(0,0).setZero();

    for (int i=0; i<n; i++)
    {
        V = InvAd(Adlist.block<6,6>(i*6,0))*V + E_tilde_list.col(i)*dq(i);
        ad = adjop(V);
        Blist.block<6,6>((i+1)*6,0) = Alist.block<6,6>((i+1)*6,0)*ad - ad.transpose()*Alist.block<6,6>((i+1)*6,0);
    }

    Astarlist.block<6,6>(6*n,0) = Alist.block<6,6>(6*n,0);
    Bstarlist.block<6,6>(6*n,0) = Blist.block<6,6>(6*n,0);

    for (int i=n; i>0; i--)
    {
        Astar = Astarlist.block<6,6>(i*6,0);
        Ad = Adlist.block<6,6>((i-1)*6, 0);
        Astarlist.block<6,6>((i-1)*6,0) = Alist.block<6,6>((i-1)*6,0) + InvTransAd(Ad) * Astar * InvAd(Ad);
        Bstarlist.block<6,6>((i-1)*6,0) = Blist.block<6,6>((i-1)*6,0) + InvTransAd(Ad)*(Bstarlist.block<6,6>(i*6,0) - Astar*adlist.block<6,6>((i-1)*6,0))*InvAd(Ad);
    }

    for (int i=0; i<n; i++)
    {
        Ad.setIdentity();
        if (i != n-1)
        {
            ad = adlist.block<6,6>((i+1)*6,0);
        }
        Ei = E_tilde_list.col(i);

        for (int j=i; j<n; j++)
        {
            Astar = Astarlist.block<6,6>((j+1)*6,0);
            Bstar = Bstarlist.block<6,6>((j+1)*6,0);
            Ej = E_tilde_list.col(j);
            if (j>=i+1) 
            {
                tmp = Adlist.block<6,6>(j*6, 0);
                Ad = Ad * tmp;
                if (j != i+1)
                {
                    ad = adlist.block<6,6>(j*6,0) + InvAd(tmp) * ad * tmp;
                }
            }
            
            
            C(i,j) = Ei.transpose() * InvTransAd(Ad) * (Astar*E_tilde_dot_list.col(j) + Bstar*Ej);

            if (i != j)
            {
                C(j,i) = Ej.transpose() * ((Bstar - Astar*ad)*InvAd(Ad)*Ei + Astar*InvAd(Ad)*E_tilde_dot_list.col(i));
            }
        }
    }
    return C;
}

VectorXf systemGravity()
{
    VectorXf G(n);
    VectorXf gref(6);
    VectorXf gtmp(6);
    MatrixXf Ad(6,6);
    MatrixXf g(6,n);

    float tmp;

    gref << 0,0,0,0,0,-9.81;

    for (int i=0; i<n; i++)
    {
        Ad = Adlist.block<6,6>(i*6, 0);

        if (i==0) gtmp = InvAd(Ad) * gref;
        else gtmp = InvAd(Ad) * gtmp;
        g.col(i) = gtmp;
    }

    for (int i=0; i<n; i++) 
    {
        tmp = 0;
        Ad.setIdentity();
        for (int j=i; j<n; j++)
        {
            if (j>=i+1) Ad = Ad * Adlist.block<6,6>(j*6, 0);
            
            tmp = tmp + E_tilde_list.col(i).transpose() * InvTransAd(Ad) * Alist.block<6,6>((j+1)*6, 0) * g.col(j);
        }
        G(i) = -tmp;
    }
    return G;
}

void RNE(VectorXf V0, VectorXf dV0, VectorXf q, VectorXf dq, VectorXf ddq, MatrixXf Fextlist)
{
    //V0 = (0,0,0,0,0,0), dV0 = (0,0,0,0,0,-9.81) for fixed-grounded base
    Matrix4f Ttmp;
    MatrixXf Vlist(6,n);
    MatrixXf Vdotlist(6,n);
    MatrixXf Adtmp(6,6);
    MatrixXf adtmp(6,6);
    MatrixXf Teef(4,4);
    VectorXf F(6);
    VectorXf Vtmp(6);

    //FR
    for (int i=0; i<n; i++)
    {
        Vtmp = E_tilde_list.col(i)*dq(i);
        adlist.block<6,6>(i*6,0) = adjop(Vtmp);
        if (i==0)   
        {
            Ttmp = expse3(lambda_list.col(i)*q(i))*T0list.block<4,4>(i*4,0);
            Adtmp = Adjoint(Ttmp);
            Vlist.col(i) = InvAd(Adtmp)*V0 + Vtmp;
            Vdotlist.col(i) = InvAd(Adtmp)*dV0 - adlist.block<6,6>(i*6,0)*InvAd(Adtmp)*V0 + E_tilde_list.col(i)*ddq(i) + E_tilde_dot_list.col(i)*dq(i);
        }
        else
        {
            Ttmp = Tinv(T0list.block<4,4>((i-1)*4,0))*expse3(lambda_list.col(i)*q(i))*T0list.block<4,4>(i*4,0);
            Adtmp = Adjoint(Ttmp);
            Vlist.col(i) = InvAd(Adtmp)*Vlist.col(i-1) + Vtmp;
            Vdotlist.col(i) = InvAd(Adtmp)*Vdotlist.col(i-1) - adlist.block<6,6>(i*6,0)*InvAd(Adtmp)*Vlist.col(i-1) + E_tilde_list.col(i)*ddq(i) + E_tilde_dot_list.col(i)*dq(i);
        }        
        Tlist.block<4,4>(i*4, 0) = Ttmp;
        Adlist.block<6,6>(i*6, 0) = Adtmp;
    }
    Teef = Tinv(Tlist.block<4,4>((n-1)*4,0))*Teef0;
    Adtmp = Adjoint(Teef);

    //BR
    F = Fextlist.col(n+1);
    for (int i=n-1; i>=0; i--)
    {
        adtmp = adjop(Vlist.col(i)).transpose();
        F = Fextlist.col(i+1) - Alist.block<6,6>((i+1)*6,0)*Vdotlist.col(i) + adtmp*Alist.block<6,6>((i+1)*6,0)*Vlist.col(i) + InvTransAd(Adtmp)*F;
        tau(i) = -E_tilde_list.col(i).transpose() * F;
        Adtmp = Adlist.block<6,6>(i*6, 0);
    }
    adtmp = adjop(V0).transpose();
    Fc = Fextlist.col(0) - Alist.block<6,6>(0,0)*dV0 + adtmp*Alist.block<6,6>(0,0)*V0 + InvTransAd(Adtmp)*F;
}

void pRNE(VectorXf V0, VectorXf q, VectorXf dq, MatrixXf Fextlist, VectorXf V0ref, VectorXf dV0ref, VectorXf dqref, VectorXf ddqref)
{
    //V0 = (0,0,0,0,0,0), dV0 = (0,0,0,0,0,-9.81) for fixed-grounded base
    Matrix4f Ttmp;
    MatrixXf Vlist(6,n);
    MatrixXf Vreflist(6,n);
    MatrixXf dVreflist(6,n);
    MatrixXf Adtmp(6,6);
    MatrixXf adtmp(6,6);
    MatrixXf A(6,6);
    MatrixXf Teef(4,4);
    VectorXf F(6);
    VectorXf Ei(6);
    VectorXf Vtmp(6);

    //FR
    Ttmp.setIdentity();
    auto V = V0;
    auto Vref = V0ref;
    auto dVref = dV0ref;

    for (int i=0; i<n; i++)
    {
        Ei = E_tilde_list.col(i);
        Ttmp = Ttmp*expse3(lambda_list.col(i)*q(i))*T0list.block<4,4>(i*4,0);
        Adtmp = Adjoint(Ttmp);

        Vtmp = Ei*dq(i);
        adtmp = adjop(Vtmp);

        dVref = InvAd(Adtmp)*dVref - adtmp*InvAd(Adtmp)*Vref + Ei*ddqref(i) + E_tilde_dot_list.col(i)*dqref(i);
        Vref = InvAd(Adtmp)*Vref + Ei*dqref(i);
        V = InvAd(Adtmp)*V + Vtmp;
            
        Tlist.block<4,4>(i*4, 0) = Ttmp;
        Adlist.block<6,6>(i*6, 0) = Adtmp;
        Vlist.col(i) = V;
        Vreflist.col(i) = Vref;
        dVreflist.col(i) = dVref;

        Ttmp = Tinv(T0list.block<4,4>(i*4,0));
    }
    Teef = Tinv(T0list.block<4,4>((n-1)*4,0))*Teef0;
    Adtmp = Adjoint(Teef);

    //BR
    F = Fextlist.col(n+1);
    for (int i=n-1; i>=0; i--)
    {
        adtmp = adjop(Vlist.col(i));
        A = Alist.block<6,6>((i+1)*6,0);
        F = Fextlist.col(i+1) - A*dVreflist.col(i) + (adtmp.transpose()*A - A*adtmp)*Vreflist.col(i) + InvTransAd(Adtmp)*F;
        // std::cout << Adtmp << std::endl;
        // std::cout << adtmp << std::endl;

        tau(i) = -E_tilde_list.col(i).transpose() * F;
        Adtmp = Adlist.block<6,6>(i*6, 0);
    }
    
    adtmp = adjop(V0);
    A = Alist.block<6,6>(0,0);
    Fc = Fextlist.col(0) - A*dV0ref + (adtmp*A - A*adtmp)*V0ref + InvTransAd(Adtmp)*F;
}

void RNEdynamics(VectorXf V0, VectorXf q, VectorXf dq)
{
    VectorXf zeros(n);
    VectorXf e(n);
    VectorXf g(6);
    MatrixXf Fext(6,n+2); //F0, F1, ... , Fn, Feef

    zeros.setZero();
    Fext.setZero();
    g << 0,0,0,0,0,-9.81;

    for (int i=0; i<n; i++)
    {
        e.setZero();
        e(i) = 1;
        RNE(zeros, zeros, q, zeros, e, Fext);
        M.col(i) = tau;
    }
    RNE(V0, zeros, q, dq, zeros, Fext);
    c = tau;
    for (int i=0; i<n; i++)
    {
        e.setZero();
        e(i) = 1;
        pRNE(V0, q, dq, Fext, zeros, zeros, e, zeros);
        C.col(i) = tau;
    }
    RNE(zeros, g, q, zeros, zeros, Fext);
    G = tau;
}

int main()
{
    VectorXf q(6);
    VectorXf dq(6);
    VectorXf g(6);
    VectorXf zeros(6);
    Matrix4f T;
    Vector4f t;


    zeros.setZero();

    q << 0,0,0,0,0,0;
    dq << 1,1,1,1,1,1;

    

    setConstant();
    updateFKList(q, dq);

    // std::cout << adlist << std::endl;

    M = systemInertia();
    C = systemBias(dq);
    G = systemGravity();
    std::cout << M << std::endl << C << std::endl << G << std::endl;
    
    RNEdynamics(zeros, q, dq);
    std::cout << M << std::endl << C << std::endl << G << std::endl;
    // std::cout << c << std::endl;


    
}