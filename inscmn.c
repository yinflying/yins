/**
 * @file inscmn.c
 * @brief INS common functions, include basic vector,matrix,quaternion
 *      operations.
 * @author yinflying(yinfying@foxmail.com)
 * @version 0.0.1
 * @note
 *  2019-06-11  Created Files
 *  2019-06-14  Add earth_RN and earth_RE function
 *  2019-06-21  Add v3_normalize
 *  2019-06-27  Add m3_inv
 */
/*
 * Copyright (c) 2019 yinflying <yinflying@foxmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or any
 * later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see accompanying file LICENSE.txt
 * or <http://www.gnu.org/licenses/>.
 */

#include "ins.h"
#include <math.h>
#include <stdlib.h>

#define PI  3.14159265358979
#define EPS 1E-50
#define SQR(x)  ((x) * (x))
static const m3_t I3 = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
static const m3_t O3 = {0.0};

/**
 * @brief vector outer product/cross product
 * @param[in] v1    First vector
 * @param[in] v2    Second vector
 * @return cross product result of v1 and v2
 */
extern v3_t v3_cross(v3_t v1, v3_t v2)
{
    v3_t v;
    v.i = v1.j * v2.k - v1.k * v2.j;
    v.j = v1.k * v2.i - v1.i * v2.k;
    v.k = v1.i * v2.j - v1.j * v2.i;
    return v;
}
/**
 * @brief Sum of the correspoding elements of two vector
 * @param[in] v1    Fisrt vector
 * @param[in] v2    Second vector
 * @return sum of v1 and v2 (v1 + v2)
 */
extern v3_t v3_add(v3_t v1, v3_t v2)
{
    return (v3_t) { v1.i + v2.i, v1.j + v2.j, v1.k + v2.k };
}
/**
 * @brief Substraction of corresponding elements of two vector
 * @param[in] v1    Fisrt vector
 * @param[in] v2    Second vector
 * @return Substraction of v1 and v2 ( v1 - v2 )
 */
extern v3_t v3_del(v3_t v1, v3_t v2)
{
    return (v3_t) { v1.i - v2.i, v1.j - v2.j, v1.k - v2.k };
}
/**
 * @brief Scalar multiplication between number and vector
 * @param[in] s     Input number
 * @param[in] v     Input vector
 * @return Scalar muliplication result of s and v ( s x v )
 */
extern v3_t v3_scalar(double s, v3_t v)
{
    return (v3_t) { s * v.i, s * v.j, s * v.k };
}

/**
 * @brief 3D vector normalize(to an unit vector)
 * @param[in,out] v input & output vector
 * @return  0: OK, 1: Input vector are zero vector
 * @see v3_norm()
 * @warning  if input vector length less than EPS(Zero determination threshold),
 *      v3_normalize() will do nothing.
 */
extern int v3_normalize(v3_t *v){
    double norm = v3_norm(*v);
    if(fabs(norm) < EPS) return 1;
    double fac = 1.0 / norm;
    v->i *= fac; v->j *= fac; v->k *= fac;
    return 0;
}

/**
 * @brief L2-norm/Euclidean Metric/Euclidean Distance of vector
 * @param[in] v Input vector
 * @return L2-norm of v
 * @note Ref: http://mathworld.wolfram.com/VectorNorm.html
 */
extern double v3_norm(v3_t v)
{
    return sqrt(v.i * v.i + v.j * v.j + v.k * v.k);
}

/**
 * @brief row vector multiply by column vector, also called vector dot product
 * @param[in] v1    row vector
 * @param[in] v2    column vector
 * @see v3_mul_cxr()
 * @return product result number of v1 and v2
 */
extern double v3_mul_rxc(v3_t v1, v3_t v2)
{
    return v1.i * v2.i + v1.j * v2.j + v1.k * v2.k;
}

/**
 * @brief column vector multiply by row vector
 * @param[in] v1    column vector
 * @param[in] v2    row vector
 * @return product result matrix of v1 and v2
 */
extern m3_t v3_mul_cxr(v3_t v1, v3_t v2)
{
    return (m3_t) { v1.i * v2.i, v1.i * v2.j, v1.i * v2.k, v1.j * v2.i,
        v1.j * v2.j, v1.j * v2.k, v1.k * v2.i, v1.k * v2.j, v1.k * v2.k };
}
/**
 * @brief Form 3D diagonal matrix by using 3D vector
 * @param[in] v Input vector
 * @return  3D diagonal matrix
 */
extern m3_t v3_diag(v3_t v)
{
    return (m3_t) {v.i, 0.0, 0.0,   0.0, v.j, 0.0,   0.0, 0.0, v.k};
}
/**
 * @brief Power exponent of every elements in vector
 * @param[in] v     Input vector
 * @param[in] order Power order, e.g. 2.0 means square, 0.5 means root square
 * @return Power exponent result vector
 */
extern v3_t v3_pow(v3_t v, double order)
{
    return (v3_t){pow(v.i,order), pow(v.j, order), pow(v.k, order)};
}

/**
 * @brief Check two 3D vector are equal or not
 * @param[in] v1    Fisrt 3D vector
 * @param[in] v2    Second 3D vector
 * @param[in] eps   Zero threshold, e.g. 1E-10
 * @return true: v1 eqaul to v2, false: v1 not equal to v2
 * @note    two vector's corresponding elements difference should be less than
 *  zero threshold, then function return true
 */
extern bool v3_equal(const v3_t *v1, const v3_t *v2, double eps)
{
    double diff;
    const double *pv1 = (double *)v1;
    const double *pv2 = (double *)v2;
    for(int i = 0; i < 3; ++i) {
        if((diff = pv1[i] - pv2[i]) > eps || diff < -eps)
            return false;
    }
    return true;
}
/**
 * @brief 3D matrix transposition
 * @param[in] A Input 3D matrix
 * @return transposition matrix of A
 */
extern m3_t m3_transpose(m3_t A)
{
    m3_t dcm;
    dcm.m11 = A.m11, dcm.m12 = A.m21, dcm.m13 = A.m31;
    dcm.m21 = A.m12, dcm.m22 = A.m22, dcm.m23 = A.m32;
    dcm.m31 = A.m13, dcm.m32 = A.m23, dcm.m33 = A.m33;
    return dcm;
}

/**
 * @brief 3D matrix inverse
 * @param[in,out] A Input 3D matrix, output inverse matrix of A
 * @return 0: OK
 */
extern int m3_inv(m3_t* A)
{
    /* PLU decomposition */
    m3_t P, L, U;
    if(!m3_LU(A,&L,&U,&P)){
        /* inverse of L */
        L.m21 = L.m21 / ( L.m22 * L.m11 );
        L.m31 = L.m31 / ( L.m33 * L.m11 ) - L.m32 * L.m21 / L.m33;
        L.m22 = - 1.0/L.m22;
        L.m32 = - L.m32 * L.m22 / L.m33;
        L.m11 = - 1.0/L.m11;
        L.m33 = - 1.0/L.m33;

        /* inverse of U */
        U.m12 = U.m12 / ( U.m22 * U.m11 );
        U.m13 = U.m13 / ( U.m33 * U.m11 ) - U.m23 * U.m12 / U.m33;
        U.m22 = - 1.0/U.m22;
        U.m23 = - U.m23 * U.m22 / U.m33;
        U.m11 = - 1.0/U.m11;
        U.m33 = - 1.0/U.m33;
    }
    *A = m3_mul(U,m3_mul(L,P));
    return 0;
}

/**
 * @brief Sum of corresponding elements of two 3D matrix
 * @param[in] A Fisrt matrix
 * @param[in] B Second matrix
 * @return Sum of A and B (A + B)
 */
extern m3_t m3_add(m3_t A, m3_t B)
{
    return (m3_t){A.m11 + B.m11, A.m12 + B.m12, A.m13 + B.m13,
        A.m21 + B.m21, A.m22 + B.m22, A.m23 + B.m23,
        A.m31 + B.m31, A.m32 + B.m32, A.m33 + B.m33};
}

/**
 * @brief Substraction of corresponding elements of two matrix
 * @param[in] A Frist 3D matrix
 * @param[in] B Second 3D matrix
 * @return A minus B result (A - B)
 */
extern m3_t m3_del(m3_t A, m3_t B)
{
    return (m3_t){A.m11 - B.m11, A.m12 - B.m12, A.m13 - B.m13,
        A.m21 - B.m21, A.m22 - B.m22, A.m23 - B.m23,
        A.m31 - B.m31, A.m32 - B.m32, A.m33 - B.m33};
}

/**
 * @brief Scalar multiplication between number and 3D matrix
 * @param[in] s     Input number
 * @param[in] A     Input 3D matrix
 * @return Scalar muliplication result of s and A ( s x v )
 */
extern m3_t m3_scalar(double alpha, m3_t A)
{
    m3_t mat;
    mat.m11 = alpha * A.m11, mat.m12 = alpha * A.m12, mat.m13 = alpha * A.m13;
    mat.m21 = alpha * A.m21, mat.m22 = alpha * A.m22, mat.m23 = alpha * A.m23;
    mat.m31 = alpha * A.m31, mat.m32 = alpha * A.m32, mat.m33 = alpha * A.m33;
    return mat;
}

/**
 * @brief matrix multiplication
 * @param[in] A Fisrt 3D matrix
 * @param[in] B Second 3D matrix
 * @return matrix multiplication result of A and B (A * B)
 */
extern m3_t m3_mul(m3_t A, m3_t B)
{
    m3_t C;
    C.m11 = A.m11 * B.m11 + A.m12 * B.m21 + A.m13 * B.m31;
    C.m12 = A.m11 * B.m12 + A.m12 * B.m22 + A.m13 * B.m32;
    C.m13 = A.m11 * B.m13 + A.m12 * B.m23 + A.m13 * B.m33;
    C.m21 = A.m21 * B.m11 + A.m22 * B.m21 + A.m23 * B.m31;
    C.m22 = A.m21 * B.m12 + A.m22 * B.m22 + A.m23 * B.m32;
    C.m23 = A.m21 * B.m13 + A.m22 * B.m23 + A.m23 * B.m33;
    C.m31 = A.m31 * B.m11 + A.m32 * B.m21 + A.m33 * B.m31;
    C.m32 = A.m31 * B.m12 + A.m32 * B.m22 + A.m33 * B.m32;
    C.m33 = A.m31 * B.m13 + A.m32 * B.m23 + A.m33 * B.m33;
    return C;
}

/**
 * @brief 3D matrix multiply 3D column vector
 * @param[in] A     Input 3D matrix
 * @param[in] B     Input 3D column vector
 * @return result 3D vector (A * B)
 */
extern v3_t m3_mul_v3(m3_t A, v3_t B)
{
    v3_t C;
    C.i = A.m11 * B.i + A.m12 * B.j + A.m13 * B.k;
    C.j = A.m21 * B.i + A.m22 * B.j + A.m23 * B.k;
    C.k = A.m31 * B.i + A.m32 * B.j + A.m33 * B.k;
    return C;
}

/**
 * @biref main diagonal elements of a 3D matrix
 * @param[in] A     Input 3D matrix
 * @return 3D vector of main diagonal elements
 */
extern v3_t m3_diag(m3_t A)
{
    return (v3_t){A.m11, A.m22, A.m33};
}

/**
 * @brief Power exponent of every elements in vector
 * @param[in] A     Input vector
 * @param[in] order Power order, e.g. 2.0 means square, 0.5 means root square
 * @return Power exponent result vector
 */
extern m3_t m3_pow(m3_t A, double order)
{
    return (m3_t){pow(A.m11, order), pow(A.m12, order), pow(A.m13, order),
    pow(A.m21, order), pow(A.m22, order), pow(A.m23, order),
    pow(A.m31, order), pow(A.m32, order), pow(A.m33, order)};
}

/**
 * @brief judge the two matrix if equal or not
 * @param[in] A     First matrix
 * @param[in] B     Second matrix
 * @param[in] eps   the zero threshold e.g. 1e-20
 * @return true: A == B, false: A != B
 */
extern bool m3_equal(const m3_t *A, const m3_t *B, double eps)
{
    double diff;
    const double *pA = (double *)A;
    const double *pB = (double *)B;
    for(int i = 0; i < 9; ++i) {
        if((diff = pA[i] - pB[i]) > eps || diff < -eps)
            return false;
    }
    return true;
}

/**
 * @brief Singular Value Decomposition(SVD) of 3D Matrix
 *          A = U * diag(D) * V'
 * @param[in] A   3D matrix
 * @param[out] U  3D unit othogonal matrix
 * @param[out] D  3D vector
 * @param[out] V  3D unit othogonal matrix
 * @return 0: OK
 */
extern int m3_SVD(const m3_t *A, m3_t *U, v3_t *D, m3_t *V)
{
    m3_t B = m3_mul(m3_transpose(*A),*A);
    /* Jacobi iteration method to solve eigenvalue and eigenvector of B */
    const double thres = 1E-15;
    m3_t eVal = B, eVec = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
    double afa, sa, ca;
    m3_t Upq;
    for(int i = 0; i < 100; i++){
        double d12 = fabs(eVal.m12), d13 = fabs(eVal.m13), d23 = fabs(eVal.m23);
        if(d12 > thres && d12 >= d13 && d12 >= d23) {
            afa = atan2(2*eVal.m12, eVal.m11 - eVal.m22)/2.0;
            sa = sin(afa), ca = cos(afa);
            Upq = (m3_t){ca,-sa,0.0, sa,ca,0.0, 0.0, 0.0, 1.0};
        }
        else if(d13 > thres && d13 >= d23 && d13 >= d12){
            afa = atan2(2*eVal.m13, eVal.m33 - eVal.m11)/2.0;
            sa = sin(afa), ca = cos(afa);
            Upq = (m3_t){ca, 0.0, sa, 0.0, 1.0, 0.0, -sa, 0.0, ca};
        }
        else if(d23 > thres && d23 >= d12 && d23 >= d13){
            afa = atan2(2*eVal.m23, eVal.m22 - eVal.m33)/2.0;
            sa = sin(afa), ca = cos(afa);
            Upq = (m3_t){1.0, 0.0, 0.0, 0.0, ca, -sa, 0.0, sa, ca};
        }else{
            break;
        }
        eVec = m3_mul(eVec,Upq);
        eVal = m3_mul(m3_transpose(eVec),m3_mul(B,eVec));
    }
    /* Using eigenvalue&eigenvecotor to calculate U D V */
    *D = v3_pow(m3_diag(eVal),0.5);
    *V = eVec;
    *U = m3_mul(*A,m3_mul(eVec,v3_diag(v3_pow(*D,-1.0))));
    return 0;
}

/**
 * @brief swap 3D matrix row
 * @param[in,out]   A   Input and swap row matrix(swap row1 and row2)
 * @param[in]       r1  row1, range: 1~3
 * @param[in]       r2  row2, range: 1~3
 */
extern void m3_swap_row(m3_t *A, int r1, int r2)
{
    if(r1 == r2) return;
    double *pA = (double *)A;
    double tmp;
    int rr1 = r1 - 1, rr2 = r2 - 1;
    for(int i = 0; i < 3; ++i){
        tmp = pA[rr1*3 + i];
        pA[rr1*3 + i] = pA[rr2*3 +i];
        pA[rr2*3 + i] = tmp;
    }
}

/**
 * @brief swap 3D matrix column
 * @param[in,out]   A   Input and swap column matrix(swap column1 and column2)
 * @param[in]       c1  column1, range: 1~3
 * @param[in]       c2 column2, range: 1~3
 */
extern void m3_swap_clm(m3_t *A, int c1, int c2)
{
    if(c1 == c2) return;
    double *pA = (double *)A;
    double tmp;
    int cc1 = c1 - 1, cc2 = c2 - 1;
    for(int i = 0; i < 3; ++i){
        tmp = pA[i*3 + cc1];
        pA[i*3 + cc1] = pA[i*3 + cc2];
        pA[i*3 + cc2] = tmp;
    }
}

/**
 * @brief determinant value of a 3D matrix
 * @param[in] A Input matrix
 * @return determinant value
 */
extern double m3_det(const m3_t *A)
{
    return A->m11 * A->m22 * A->m33 + A->m12 * A->m23 * A->m31
        + A->m13 * A->m21 * A->m32 - A->m31 * A->m22 * A->m13
        - A->m21 * A->m12 * A->m33 - A->m11 * A->m32 * A->m23;
}

/**
 * @brief 3D matrix LU decomposition, express matrix A as the product of two
 *      essentially triangular matrices, and satisfying the equation:
 *      PA = LU, where L is lower triangular matrix, and U is Upper trianglur
 *      matrix, P is the row permutation matrix.
 * @param[in]   A   Input matrix
 * @param[out]  L   Lower triangular matrix
 * @param[out]  U   Upper triangular matrix
 * @param[out]  P   Row permutation matrix
 * @return 0: OK
 * @note  LU decompostion by Elimination with Maximal Column Pivoting
 */
extern int m3_LU(const m3_t *A, m3_t *L, m3_t *U, m3_t *P)
{
    *L = I3; *U = *A; *P = I3;
    double *pU = (double *)U, *pL = (double *)L;
    for(int i = 0; i < 2; ++i){
        /* choose the abs max column element as pivot */
        if(i == 0){
            if(fabs(U->m11) >= fabs(U->m21)){
                if(fabs(U->m11) < fabs(U->m31)){
                    m3_swap_row(U,1,3);
                    m3_swap_row(P,1,3);
                }
            }else{
                if(fabs(U->m21) >= fabs(U->m31)) {
                    m3_swap_row(U,1,2);
                    m3_swap_row(P,1,2);
                }else{
                    m3_swap_row(U,1,3);
                    m3_swap_row(P,1,3);
                }
            }
        }else if(i == 1){
            if(fabs(U->m22) < fabs(U->m32)){
                m3_swap_row(U,2,3);
                m3_swap_row(P,2,3);
                /* Should exchange the m31 and m21 of L */
                double tmp = L->m31; L->m31 = L->m21; L->m21 = tmp;
            }
        }
        /* solve pivot element equal to zero */
        if(fabs(pU[i*3 +i]) < 1E-64)
            continue;
        /* Gaussian elimination */
        for(int j = i+1; j < 3; ++j){
            pL[j*3 + i] = pU[j*3 + i] / pU[i*3 + i];
            for(int k = 0; k < 3; ++k){
                pU[j*3 + k] -= pL[j*3 + i] * pU[i*3 + k];
            }
        }
    }
    return 0;
}

earth_t wgs84 = { .wie = 7.292115E-5,
    .R0 = 6378137,
    .RP= 6356752.31425,
    .mu = 3.986004418E14,
    .J2 = 1.082627E-3,
    .e = 0.0818191908425,
    .f = 1.0 / 298.257223563
};

/**
 * @brief Caculate meridian radius of curvature(from the north to south)
 * @param[in] eth   earth parameters struct
 * @param[in] lat   latitude [rad]
 * @return meridian radius of curvature [m]
 * @see earth_RE()
 * @note
 *      Ref: Paul D. Groves, Principles of GNSS, Inertial, and Multisensor
 *          Integrated Navigation Systems(2nd Edition), 2013, P59
 */
extern double earth_RN(const earth_t *eth, double lat)
{
    double e2 = SQR(eth->e);
    return eth->R0 * (1-e2) * pow(1 - e2*SQR(sin(lat)), -1.5);
}

/**
 * @brief Caculate transverse radius of curvature(from the east to west)
 * @param[in] eth   earth parameters struct
 * @param[in] lat   latitude [rad]
 * @return transverse radius of curvature [m]
 * @see earth_RN()
 * @note Also called nomal radius of curvature, or prime vertical radius of
 *  curvature.
 *
 *      Ref: Paul D. Groves, Principles of GNSS, Inertial, and Multisensor
 *          Integrated Navigation Systems(2nd Edition), 2013, P59
 */
extern double earth_RE(const earth_t *eth, double lat)
{
    return eth->R0 / sqrt(1 - SQR(eth->e) * SQR(sin(lat)));
}

/**
 * @brief  get anti-symmetric matrix of a 3D vector
 * @param[in] v3    Input 3D vector
 * @param[in] mat   result anti-symmetric matrix
 * @return 0: OK
 */
int asymmetric_mat(const v3_t* v3, m3_t* mat)
{
    mat->m11 = 0, mat->m12 = -v3->k, mat->m13 = v3->j;
    mat->m21 = v3->k, mat->m22 = 0, mat->m23 = -v3->i;
    mat->m31 = -v3->j, mat->m32 = v3->i, mat->m33 = 0;
    return 0;
}

/**
 * @brief Convert calendar day/time to gtime_t struct
 * @param[in] ep    day/time {year,month,day,hour,min,sec}
 * @return gtime_t struct
 * @note proper in 1970-2037 or 1970-2099 (64bit time_t)
 */
extern gtime_t yins_epoch2time(const double* ep)
{
    const int doy[] = { 1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335 };
    gtime_t time = { 0 };
    int days, sec, year = (int)ep[0], mon = (int)ep[1], day = (int)ep[2];

    if (year < 1970 || 2099 < year || mon < 1 || 12 < mon)
        return time;

    /* leap year if year%4==0 in 1901-2099 */
    days = (year - 1970) * 365 + (year - 1969) / 4 + doy[mon - 1] + day - 2
        + (year % 4 == 0 && mon >= 3 ? 1 : 0);
    sec = (int)floor(ep[5]);
    time.time
        = (time_t)days * 86400 + (int)ep[3] * 3600 + (int)ep[4] * 60 + sec;
    time.sec = ep[5] - sec;
    return time;
}

/**
 * @brief Convert gps time struct to calendar day/time
 * @param[in] week  gps week
 * @param[in] sec   gps second of week[s]
 * @return gtime_t struct
 * @note proper in 1970-2037 or 1970-2099 (64bit time_t)
 */
extern gtime_t yins_gpst2time(int week, double sec)
{
    const double gpst0[] = { 1980, 1, 6, 0, 0, 0 }; /* gps time reference */
    gtime_t t = yins_epoch2time(gpst0);

    if (sec < -1E9 || 1E9 < sec)
        sec = 0.0;
    t.time += (time_t)86400 * 7 * week + (int)sec;
    t.sec = sec - (int)sec;
    return t;
}

/**
 * @brief convert gtime_t struct to calendar day/time
 * @param[in] t     gtime_t struct
 * @param[out] ep   day/time {year,month,day,hour,min,sec}
 * @note Proper in 1970-2037 or 1970-2099 (64bit time_t)
 */
extern void yins_time2epoch(gtime_t t, double* ep)
{
    const int mday[] = { /* # of days in a month */
        31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31, 31, 28, 31, 30, 31, 30,
        31, 31, 30, 31, 30, 31, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31,
        31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31
    };
    int days, sec, mon, day;

    /* leap year if year%4==0 in 1901-2099 */
    days = (int)(t.time / 86400);
    sec = (int)(t.time - (time_t)days * 86400);
    for (day = days % 1461, mon = 0; mon < 48; mon++) {
        if (day >= mday[mon])
            day -= mday[mon];
        else
            break;
    }
    ep[0] = 1970 + days / 1461 * 4 + mon / 12;
    ep[1] = mon % 12 + 1;
    ep[2] = day + 1;
    ep[3] = sec / 3600;
    ep[4] = sec % 3600 / 60;
    ep[5] = sec % 60 + t.sec;
}

/**
 * @brief Convert Euler attitude to Quaternion attitude(Eab => Qab)
 * @param[in] euler Input Euler attitude
 * @param[out] quat Ouput Quaternion attitude
 * @return 0: OK
 */
extern int euler2quat(const v3_t* euler, quat_t* quat)
{
    double si = sin(euler->i / 2), ci = cos(euler->i / 2);
    double sj = sin(euler->j / 2), cj = cos(euler->j / 2);
    double sk = sin(euler->k / 2), ck = cos(euler->k / 2);
    quat->q0 = ci * cj * ck + si * sj * sk;
    quat->q1 = si * cj * ck - ci * sj * sk;
    quat->q2 = ci * sj * ck + si * cj * sk;
    quat->q3 = ci * cj * sk - si * sj * ck;
    quat_conj(quat);
    return 0;
}

/**
 * @brief Convert Quaternion attitude to Euler attitude(Qab => Eab)
 * @param[in]  quat     Input Quaternion attitude
 * @param[out] euler    Output Euler attitue
 * @return  0: OK
 */
extern int quat2euler(const quat_t* quat, v3_t* euler)
{
    euler->i = atan2(2 * (-quat->q0 * quat->q1 + quat->q2 * quat->q3),
        1 - 2 * quat->q1 * quat->q1 - 2 * quat->q2 * quat->q2);
    euler->j = asin(2 * (-quat->q0 * quat->q2 - quat->q1 * quat->q3));
    euler->k = atan2(2 * (-quat->q0 * quat->q3 + quat->q1 * quat->q2),
        1 - 2 * quat->q2 * quat->q2 - 2 * quat->q3 * quat->q3);

    if (euler->i <= -PI) /* Limit Roll Angle to (-pi,pi] */
        euler->i += 2 * PI;
    else if (euler->i > PI)
        euler->i -= 2 * PI;
    if (euler->k < 0) /* Limit Heading Angle to [0,2pi) */
        euler->k += 2 * PI;
    /* Pitch Angle limit to [-pi/2,pi/2], asin return value range */
    return 0;
}

/**
 * @brief Convert DCM attitude to Euler attitude
 * @param[in]   dcm     Input DCM attitude
 * @param[out]  euler   Ouput Euler attitude
 * @return  0: OK
 */
int dcm2euler(const m3_t* dcm, v3_t* euler)
{
    euler->i = atan2(dcm->m23, dcm->m33);
    euler->j = -asin(dcm->m13);
    euler->k = atan2(dcm->m12, dcm->m11);

    if (euler->i <= -PI) /* Limit Roll Angle to (-pi,pi] */
        euler->i += 2 * PI;
    else if (euler->i > PI)
        euler->i -= 2 * PI;
    if (euler->k < 0) /* Limit Heading Angle to [0,2pi) */
        euler->k += 2 * PI;
    /* Pitch Angle limit to [-pi/2,pi/2], asin return value range */
    return 0;
}

/**
 * @brief Convert Euler to DCM attitude
 * @param[in] euler     Input Euler attitude
 * @param[out] dcm      Ouput DCM attitude
 * @return  0: OK
 */
int euler2dcm(const v3_t* euler, m3_t* dcm)
{
    double sin_phi = sin(euler->i);
    double cos_phi = cos(euler->i);
    double sin_theta = sin(euler->j);
    double cos_theta = cos(euler->j);
    double sin_psi = sin(euler->k);
    double cos_psi = cos(euler->k);

    /* Calculate coordinate transformation matrix using (2.22) */
    dcm->m11 = cos_theta * cos_psi;
    dcm->m12 = cos_theta * sin_psi;
    dcm->m13 = -sin_theta;
    dcm->m21 = -cos_phi * sin_psi + sin_phi * sin_theta * cos_psi;
    dcm->m22 = cos_phi * cos_psi + sin_phi * sin_theta * sin_psi;
    dcm->m23 = sin_phi * cos_theta;
    dcm->m31 = sin_phi * sin_psi + cos_phi * sin_theta * cos_psi;
    dcm->m32 = -sin_phi * cos_psi + cos_phi * sin_theta * sin_psi;
    dcm->m33 = cos_phi * cos_theta;
    return 0;
}

/**
 * @brief Convert DCM attitude to Quaternion attitude
 * @param[in]   dcm     DCM attitude
 * @param[out]  quat    Quaternion attitude
 * @return  O: OK
 */
extern int dcm2quat(const m3_t* dcm, quat_t* quat)
{
    double qq4;
    if (dcm->m11 >= dcm->m22 + dcm->m33) {
        quat->q1 = 0.5 * sqrt(1 + dcm->m11 - dcm->m22 - dcm->m33);
        qq4 = 4 * quat->q1;
        quat->q0 = (dcm->m32 - dcm->m23) / qq4;
        quat->q2 = (dcm->m12 + dcm->m21) / qq4;
        quat->q3 = (dcm->m13 + dcm->m31) / qq4;
    } else if (dcm->m22 >= dcm->m11 + dcm->m33) {
        quat->q2 = 0.5 * sqrt(1 - dcm->m11 + dcm->m22 - dcm->m33);
        qq4 = 4 * quat->q2;
        quat->q0 = (dcm->m13 - dcm->m31) / qq4;
        quat->q1 = (dcm->m12 + dcm->m21) / qq4;
        quat->q3 = (dcm->m23 + dcm->m32) / qq4;
    } else if (dcm->m33 >= dcm->m11 + dcm->m22) {
        quat->q3 = 0.5 * sqrt(1 - dcm->m11 - dcm->m22 + dcm->m33);
        qq4 = 4 * quat->q3;
        quat->q0 = (dcm->m21 - dcm->m12) / qq4;
        quat->q1 = (dcm->m13 + dcm->m31) / qq4;
        quat->q2 = (dcm->m23 + dcm->m32) / qq4;
    } else {
        quat->q0 = 0.5 * sqrt(1 + dcm->m11 + dcm->m22 + dcm->m33);
        qq4 = 4 * quat->q0;
        quat->q1 = (dcm->m32 - dcm->m23) / qq4;
        quat->q2 = (dcm->m13 - dcm->m31) / qq4;
        quat->q3 = (dcm->m21 - dcm->m12) / qq4;
    }
    quat_normalize(quat);
    return 0;
}

/**
 * @brief Convert Quaternion to DCM attitude
 * @param[in]   quat    Input Quaternion attitude
 * @param[out]  dcm     Oput DCM attitude
 * @return 0: OK
 */
extern int quat2dcm(const quat_t* quat, m3_t* dcm)
{
    double q11 = quat->q0 * quat->q0, q12 = quat->q0 * quat->q1,
           q13 = quat->q0 * quat->q2, q14 = quat->q0 * quat->q3,
           q22 = quat->q1 * quat->q1, q23 = quat->q1 * quat->q2,
           q24 = quat->q1 * quat->q3, q33 = quat->q2 * quat->q2,
           q34 = quat->q2 * quat->q3, q44 = quat->q3 * quat->q3;
    dcm->m11 = q11 + q22 - q33 - q44, dcm->m12 = 2 * (q23 - q14),
    dcm->m13 = 2 * (q24 + q13), dcm->m21 = 2 * (q23 + q14),
    dcm->m22 = q11 - q22 + q33 - q44, dcm->m23 = 2 * (q34 - q12),
    dcm->m31 = 2 * (q24 - q13), dcm->m32 = 2 * (q34 + q12),
    dcm->m33 = q11 - q22 - q33 + q44;
    return 0;
}

/**
 * @brief Quaternion normalization and make first element large than zero
 * @param[in,out] quat  Quaternion/unit Quaternion
 * @return  O: OK
 */
extern int quat_normalize(quat_t* quat)
{
    double nq = 1.0 / quat_norm(*quat);
    quat->q0 *= nq;
    quat->q1 *= nq;
    quat->q2 *= nq;
    quat->q3 *= nq;
    if (quat->q0 < 0) {
        quat->q0 = -quat->q0;
        quat->q1 = -quat->q1;
        quat->q2 = -quat->q2;
        quat->q3 = -quat->q3;
    }
    return 0;
}

/**
 * @brief Quaternion conjugation
 * @param[in,out] quat  Quaternion/Quaternion conjugation
 * @return 0: OK
 * @see quat_inv()
 * @note Unit Quaternion Conjugation is the same as inversion, but conjugation
 *  is much faster
 */
extern int quat_conj(quat_t *quat)
{
    quat->q1 = -quat->q1;
    quat->q2 = -quat->q2;
    quat->q3 = -quat->q3;
    return 0;
}

/**
 * @brief Quaternion inversion
 * @param[in,out] quat  Input Quaternion/Ouput Quaternion invertion
 * @see quat_dot()
 * @see quat_conj()
 * @return 0: OK, 1: input Quaternion norm too small, less than sqrt(EPS)
 * @note Unit Quaternion Conjugation is the same as inversion, but conjugation
 *  is much faster
 */
extern int quat_inv(quat_t *quat)
{
    double inv = 1.0 / quat_dot(*quat, *quat);
    if(fabs(inv) < EPS) return 1;

    quat_conj(quat);
    quat->q0 *= inv; quat->q1 *= inv; quat->q2 *= inv; quat->q3 *= inv;
    return 0;
}

/**
 * @brief Quaterion dot production, sum of corresponding elements product
 * @param[in] P     Fisrt Quaternion
 * @param[in] Q     Second Quaternion
 * @return dot product result number
 * @note dot production also called inner porduction
 */
extern double quat_dot(quat_t P, quat_t Q)
{
    return P.q0*Q.q0 + P.q1*Q.q1 + P.q2*Q.q2 + P.q3*Q.q3;
}

/**
 * @brief Quaternion cross product, result of the vector cross procdut of two
 *      quaternions' imaginary part(q1,q2,q3)
 * @param[in] P     Fisrt quaternion
 * @param[in] Q     Second quatenion
 * @return Quaternion cross product result vector
 */
extern v3_t quat_cross(quat_t P, quat_t Q)
{
    return v3_cross((v3_t){P.q1,P.q2,P.q3},(v3_t){Q.q1,Q.q2,Q.q3});
}

/**
 * @brief L2-norm of quaternion(also call modulus of quaternion)
 * @param[in] P Input quaternion
 * @return L2-norm of quaternion
 */
extern double quat_norm(quat_t P)
{
    return sqrt(P.q0*P.q0 + P.q1*P.q1 + P.q2*P.q2 + P.q3*P.q3);
}

/**
 * @brief Check two quaternions are equal or not
 * @param[in] P     Fisrt quaternion
 * @param[in] Q     Second quaternion
 * @param[in] eps   Zero threshold, e.g. 1E-10
 * @return true: P eqaul to Q, false: P not equal to Q
 * @warning minus quaternion do NOT equal to origin quaternion
 * @note    two quaternions' corresponding elements difference should be less
 *  than zero threshold, then function return true, however, minus quaterions
 *  are the same as origin if quaternion represent attitude, but this function
 *  can not identify.
 */
extern bool quat_equal(const quat_t *P, const quat_t *Q, double eps)
{
    double diff;
    const double *pP = (double *)P;
    const double *pQ = (double *)Q;
    for(int i = 0; i < 3; ++i) {
        if((diff = pP[i] - pQ[i]) > eps || diff < -eps)
            return false;
    }
    return true;
}

/**
 * @brief Rotation vector(Angular increment) to quaternion attitude
 *      transformation(New to Old)
 * @param[in]   dtheta  Rotation vector(Angular increment) [rad]
 * @param[out]  quat    Quaternion attitude transformation
 * @return O: OK
 * @see rv2dcm()
 * @note Quaternion attitude transforman can be expressed as Qbi+_bi-, And
 *      fellow the compluting: (Qbi+  = Qbi-  * Qbi+_bi- )
 *      Ref: Yan Gongming, 捷联惯导算法与组合导航原理讲义, 2016, P29
 */
extern int rv2quat(const v3_t* dtheta, quat_t* quat)
{
    const double F1 = 2 * 1; // define Fk = 2^k * k!
    const double F2 = F1 * 2 * 2;
    const double F3 = F2 * 2 * 3;
    const double F4 = F3 * 2 * 4;
    const double F5 = F4 * 2 * 5;
    double n2
        = dtheta->i * dtheta->i + dtheta->j * dtheta->j + dtheta->k * dtheta->k;
    double f;
    if (n2 < (PI / 180.0 * PI / 180.0)) {
        double n4 = n2 * n2;
        quat->q0 = 1.0 - n2 * (1.0 / F2) + n4 * (1.0 / F4);
        f = 0.5 - n2 * (1.0 / F3) + n4 * (1.0 / F5);
    } else {
        double n_2 = sqrt(n2) / 2.0;
        quat->q0 = cos(n_2);
        f = sin(n_2) / n_2 * 0.5;
    }
    quat->q1 = f * dtheta->i;
    quat->q2 = f * dtheta->j;
    quat->q3 = f * dtheta->k;
    return 0;
}

/**
 * @brief Rotation vector(Angular increment) to DCM attitude transformation
 *      (New to Old)
 * @param[in]   dtheta  Rotation vector(Angular increment) [rad]
 * @param[out]  dcm     DCM attitude transformation
 * @return  0: OK
 * @see rv2quat()
 * @note DCM attitude transforman can be expressed as Cbi+_bi-, And  fellow the
 *      compluting: (Cbi+  = Cbi-  * Cbi+_bi- )
 *      The DCM can also be expressed in "exp(skew-symmetric(dtheta))".
 *
 *      Ref: Paul D. Groves, Principles of GNSS, Inertial, and Multisensor
 *          Integrated Navigation Systems(2nd Edition), 2013, P184
 */
extern int rv2dcm(const v3_t* dtheta, m3_t* dcm)
{
    double norm_dtheta = v3_norm(*dtheta);
    m3_t Ax; asymmetric_mat(dtheta,&Ax);
    if(norm_dtheta > 1e-8){
        m3_t m1 = m3_scalar(sin(norm_dtheta)/norm_dtheta,Ax);
        m3_t m2 = m3_scalar((1-cos(norm_dtheta))/SQR(norm_dtheta),m3_mul(Ax,Ax));
        *dcm = m3_add(m3_add(I3, m1),m2);
    }else{
        *dcm = m3_add(I3, Ax);
    }
    return 0;
}

/**
 * @brief Hamiton product, is determined by the products of the basis elements
 *  and the ditributive law. Hamition product could represent the ratation
 *  quaternions
 * @param[in] P Fisrt quaternion
 * @param[in] Q Second quaternion
 * @return Hamiton product result quaternion (P * Q)
 */
extern quat_t quat_mul(quat_t P, quat_t Q)
{
    quat_t qtmp;
    qtmp.q0 = P.q0 * Q.q0 - P.q1 * Q.q1 - P.q2 * Q.q2 - P.q3 * Q.q3;
    qtmp.q1 = P.q0 * Q.q1 + P.q1 * Q.q0 + P.q2 * Q.q3 - P.q3 * Q.q2;
    qtmp.q2 = P.q0 * Q.q2 + P.q2 * Q.q0 + P.q3 * Q.q1 - P.q1 * Q.q3;
    qtmp.q3 = P.q0 * Q.q3 + P.q3 * Q.q0 + P.q1 * Q.q2 - P.q2 * Q.q1;
    return qtmp;
}

/**
 * @brief Rotate a 3D vector, to make (vn = Qbn * vb) work
 * @param[in] quat  Input Quaternion
 * @param[in] vec   Input 3D vector
 * @return result of rotation 3D vector ( quat * vec )
 */
extern v3_t quat_mul_v3(quat_t quat, v3_t vec)
{
    quat_t qtmp;
    v3_t vtmp;
    qtmp.q0 = -quat.q1 * vec.i - quat.q2 * vec.j - quat.q3 * vec.k;
    qtmp.q1 = quat.q0 * vec.i + quat.q2 * vec.k - quat.q3 * vec.j;
    qtmp.q2 = quat.q0 * vec.j + quat.q3 * vec.i - quat.q1 * vec.k;
    qtmp.q3 = quat.q0 * vec.k + quat.q1 * vec.j - quat.q2 * vec.i;
    vtmp.i = -qtmp.q0 * quat.q1 + qtmp.q1 * quat.q0 - qtmp.q2 * quat.q3
        + qtmp.q3 * quat.q2;
    vtmp.j = -qtmp.q0 * quat.q2 + qtmp.q2 * quat.q0 - qtmp.q3 * quat.q1
        + qtmp.q1 * quat.q3;
    vtmp.k = -qtmp.q0 * quat.q3 + qtmp.q3 * quat.q0 - qtmp.q1 * quat.q2
        + qtmp.q2 * quat.q1;
    return vtmp;
}

extern m3_t formCen_ned(double lat, double lon)
{
    double coslat = cos(lat), sinlat = sin(lat);
    double coslon = cos(lon), sinlon = sin(lon);
    m3_t Cen;
    Cen.m11 = -sinlat * coslon;
    Cen.m12 = -sinlat * sinlon;
    Cen.m13 = coslat;
    Cen.m21 = -sinlon;
    Cen.m22 = coslon;
    Cen.m23 = 0;
    Cen.m31 = -coslat * coslon;
    Cen.m32 = -coslat * sinlon;
    Cen.m33 = -sinlat;
    return Cen;
}

extern int ned2ecef(v3_t* pos, v3_t* vel, m3_t* dcm)
{
    double lat = pos->i, lon = pos->j, hgt = pos->k;
    double coslat = cos(lat), sinlat = sin(lat);
    double coslon = cos(lon), sinlon = sin(lon);

    double tmp = wgs84.e * sinlat;
    double Re = wgs84.R0 / sqrt(1 - tmp * tmp);

    pos->i = (Re + hgt) * coslat * coslon;
    pos->j = (Re + hgt) * coslat * sinlon;
    pos->k = ((1 - wgs84.e * wgs84.e) * Re + hgt) * sinlat;

    if (vel != NULL || dcm != NULL) {
        m3_t Cne;
        Cne.m11 = -sinlat * coslon;
        Cne.m21 = -sinlat * sinlon;
        Cne.m31 = coslat;
        Cne.m12 = -sinlon;
        Cne.m22 = coslon;
        Cne.m32 = 0;
        Cne.m13 = -coslat * coslon;
        Cne.m23 = -coslat * sinlon;
        Cne.m33 = -sinlat;
        if (vel != NULL)
            *vel = m3_mul_v3(Cne, *vel); /* Veb_n => Veb_e */
        if (dcm != NULL)
            *dcm = m3_mul(Cne, *dcm); /* Cb_n => Cb_e */
    }
    return 0;
}

extern int ecef2ned(v3_t* pos, v3_t* vel, m3_t* dcm)
{
    /* ref Pual 2012, C.29 - C.38 */
    double lon = atan2(pos->j, pos->i);

    double e2 = wgs84.e * wgs84.e;
    double k1 = sqrt(1.0 - e2) * fabs(pos->k);
    double k2 = e2 * wgs84.R0;
    double beta = sqrt(pos->i * pos->i + pos->j * pos->j);
    double E = (k1 - k2) / beta, F = (k1 + k2) / beta;
    double P = 4.0 / 3.0 * (E * F + 1.0);
    double Q = 2.0 * (E * E - F * F);
    double D = P * P * P + Q * Q;
    double V = pow(sqrt(D) - Q, 1.0 / 3.0) - pow(sqrt(D) + Q, 1.0 / 3.0);
    double G = 0.5 * (sqrt(E * E + V) + E);
    double T = sqrt(G * G + (F - V * G) / (2.0 * G - E)) - G;
    double signz = pos->k > 0 ? 1.0 : -1.0;
    double lat = signz * atan((1 - T * T) / (2 * T * sqrt(1 - e2)));

    double coslat = cos(lat), sinlat = sin(lat);
    double hgt = (beta - wgs84.R0 * T) * coslat
        + (pos->k - signz * wgs84.R0 * sqrt(1.0 - e2)) * sinlat;

    pos->i = lat, pos->j = lon, pos->k = hgt;

    if (vel != NULL || dcm != NULL) {
        double coslon = cos(lon), sinlon = sin(lon);
        m3_t Cen;
        Cen.m11 = -sinlat * coslon;
        Cen.m12 = -sinlat * sinlon;
        Cen.m13 = coslat;
        Cen.m21 = -sinlon;
        Cen.m22 = coslon;
        Cen.m23 = 0;
        Cen.m31 = -coslat * coslon;
        Cen.m32 = -coslat * sinlon;
        Cen.m33 = -sinlat;
        if (vel != NULL)
            *vel = m3_mul_v3(Cen, *vel); /* Veb_e => Veb_n */
        if (dcm != NULL)
            *dcm = m3_mul(Cen, *dcm); /* Cbe => Cbn */
    }
    return 0;
}

extern double yins_timediff(gtime_t t1, gtime_t t2)
{
    return difftime(t1.time, t2.time) + t1.sec - t2.sec;
}
