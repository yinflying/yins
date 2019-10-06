/**
 * @file inscmn.c
 * @brief INS common functions, include basic vector,matrix,quaternion
 *      operations.
 * @author yinflying(yinfying@foxmail.com)
 * @version 0.0.1
 * @note
 *  2019-06-11  Created Files \n
 *  2019-06-14  Add earth_RN and earth_RE function \n
 *  2019-06-21  Add v3_normalize \n
 *  2019-06-27  Add m3_inv \n
 *  2019-07-01  Full all functions' comments of the file \n
 *  2019-08-03  Add multiple attitude transform functions \n
 *  2019-08-10  Add cfg relatived functions \n
 *  2019-08-28  Add randn function
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
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <ctype.h>

#define PI  3.14159265358979
#define EPS 1E-50
#define SQR(x)  ((x) * (x))

const m3_t I3 = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0}; /**< unit 3D matrix */
const m3_t O3 = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; /**< zero 3D matrix */
const v3_t V0 = {0.0, 0.0, 0.0};            /**< zero 3D vector */
const v3_t V1 = {1.0, 1.0, 1.0};            /**< unit 3D vector */
v3_t WIE_E = {0.0, 0.0, 7.292115E-5};
cfg_t cfg = {
    .IPOS = 6,
    .IVEL = 3,
    .IATT = 0,
    .IBA = 9,
    .IBG = 12,
    .nx = 15,
    .isx_ba = true,
    .isx_bg = true,
    .isx_kod = false,
    .isx_eroll = false,
    .isx_epitch = false,
    .isx_eyaw = false,
    .isx_armgps = false,
    .isx_armcar = false,
    .max_ny = 3,
    .issol_header = true,
    .sol_refpos = 0,
    .feedratio = 1.0,
    .is_imu_samenoise = false,
    .log_path = "yins.log",
};

/**
 * @brief vector outer product/cross product
 * @param[in] v1    First vector
 * @param[in] v2    Second vector
 * @return cross product result of v1 and v2
 */
inline extern v3_t v3_cross(v3_t v1, v3_t v2)
{
    v3_t v;
    v.x = v1.y * v2.z - v1.z * v2.y;
    v.y = v1.z * v2.x - v1.x * v2.z;
    v.z = v1.x * v2.y - v1.y * v2.x;
    return v;
}
/**
 * @brief Sum of the correspoding elements of two vector
 * @param[in] v1    Fisrt vector
 * @param[in] v2    Second vector
 * @return sum of v1 and v2 (v1 + v2)
 */
inline extern v3_t v3_add(v3_t v1, v3_t v2)
{
    return (v3_t) { v1.x + v2.x, v1.y + v2.y, v1.z + v2.z };
}
/**
 * @brief Substraction of corresponding elements of two vector
 * @param[in] v1    Fisrt vector
 * @param[in] v2    Second vector
 * @return Substraction of v1 and v2 ( v1 - v2 )
 */
inline extern v3_t v3_del(v3_t v1, v3_t v2)
{
    return (v3_t) { v1.x - v2.x, v1.y - v2.y, v1.z - v2.z };
}
/**
 * @brief Scalar multiplication between number and vector
 * @param[in] s     Input number
 * @param[in] v     Input vector
 * @return Scalar muliplication result of s and v ( s x v )
 */
inline extern v3_t v3_scalar(double s, v3_t v)
{
    return (v3_t) { s * v.x, s * v.y, s * v.z };
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
    v->x *= fac; v->y *= fac; v->z *= fac;
    return 0;
}

/**
 * @brief L2-norm/Euclidean Metric/Euclidean Distance of vector
 * @param[in] v Input vector
 * @return L2-norm of v
 * @note Ref: http://mathworld.wolfram.com/VectorNorm.html
 */
inline extern double v3_norm(v3_t v)
{
    return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

/**
 * @brief row vector multiply by column vector, also called vector dot product
 * @param[in] v1    row vector
 * @param[in] v2    column vector
 * @see v3_mul_cxr() v3_mul_cxc() v3_cross()
 * @return product result number of v1 and v2
 */
inline extern double v3_mul_rxc(v3_t v1, v3_t v2)
{
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

/**
 * @brief column vector multiply by row vector
 * @param[in] v1    column vector
 * @param[in] v2    row vector
 * @return product result matrix of v1 and v2
 * @see v3_mul_rxc() v3_mul_cxc() v3_cross()
 */
inline extern m3_t v3_mul_cxr(v3_t v1, v3_t v2)
{
    return (m3_t) { v1.x * v2.x, v1.x * v2.y, v1.x * v2.z, v1.y * v2.x,
        v1.y * v2.y, v1.y * v2.z, v1.z * v2.x, v1.z * v2.y, v1.z * v2.z };
}

/**
 * @brief multiply coresponding elements of two vector
 * @param[in] v1    vector 1
 * @param[in] v2    vector 2
 * @return result vector
 * @see v3_mul_cxr() v3_mul_rxc() v3_cross()
 */
inline extern v3_t v3_mul_cxc(v3_t v1, v3_t v2)
{
    return (v3_t){v1.x*v2.x, v1.y*v2.y, v1.z*v2.z};
}

/**
 * @brief Form 3D diagonal matrix by using 3D vector
 * @param[in] v Input vector
 * @return  3D diagonal matrix
 */
inline extern m3_t v3_diag(v3_t v)
{
    return (m3_t) {v.x, 0.0, 0.0,   0.0, v.y, 0.0,   0.0, 0.0, v.z};
}
/**
 * @brief Power exponent of every elements in vector
 * @param[in] v     Input vector
 * @param[in] order Power order, e.g. 2.0 means square, 0.5 means root square
 * @return Power exponent result vector
 */
inline extern v3_t v3_pow(v3_t v, double order)
{
    return (v3_t){pow(v.x,order), pow(v.y, order), pow(v.z, order)};
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
    const double *pv1 = (const double *)v1;
    const double *pv2 = (const double *)v2;
    for(int i = 0; i < 3; ++i) {
        if((diff = pv1[i] - pv2[i]) > eps || diff < -eps)
            return false;
    }
    return true;
}
/**
 * @brief Square root of every elements in vector
 * @param[in] v     Input vector
 * @return Square root result vector
 */
inline extern v3_t v3_sqrt(v3_t v)
{
    return (v3_t){sqrt(v.x), sqrt(v.y), sqrt(v.z)};
}

/**
 * @brief sum of 3D vector list's corresponding elements
 * @param[in] v3_list   Input vector list
 * @param[in] n         number of vector
 * @return  sum vector
 */
extern v3_t v3_sum(const v3_t *v3_list, int n)
{
    v3_t vsum = V0;
    for(int i = 0; i < n; ++i){
        vsum.x += v3_list[i].x;
        vsum.y += v3_list[i].y;
        vsum.z += v3_list[i].z;
    }
    return vsum;
}

/**
 * @brief mean value of 3D vector list's corresponding elements
 * @param[in] v3_list   Input vector list
 * @param[in] n         number of vector
 * @return   mean value vector
 */
extern v3_t v3_mean(const v3_t *v3_list, int n)
{
    return v3_scalar(1.0/n, v3_sum(v3_list,n));
}

/**
 * @brief root-mean-square value of 3D vector list's corresponding elements
 * @param[in] v3_list   Input vector list
 * @param[in] n         number of vector
 * @return  root-mean-square vector
 */
v3_t v3_rms(const v3_t *v3_list, int n){
    v3_t vsum = V0;
    for(int i = 0; i < n; ++i){
        vsum.x += SQR(v3_list[i].x);
        vsum.y += SQR(v3_list[i].y);
        vsum.z += SQR(v3_list[i].z);
    }
    return v3_sqrt(v3_scalar(1.0/n, vsum));
}

/**
 * @brief standard deviation of 3D vector list's corresponding elements
 * @param[in] v3_list   Input vector list
 * @param[in] n         number of vector
 * @return  standard deviation vector
 */
v3_t v3_std(const v3_t *v3_list, int n)
{
    v3_t vmean = v3_mean(v3_list, n);
    v3_t vsum = V0;
    for(int i = 0; i < n; ++i){
        vsum.x += SQR(v3_list[i].x - vmean.x);
        vsum.y += SQR(v3_list[i].y - vmean.y);
        vsum.z += SQR(v3_list[i].z - vmean.z);
    }
    return v3_sqrt(v3_scalar(1.0/(n-1), vsum));
}

/**
 * @brief Convert 3D vector to 3D skew symmetric matrix.
 * @param[in] v     3D vector
 * @return skew symmetric matrix
 * @see m3_iaskew()
 */
inline extern m3_t v3_askew(v3_t v)
{
    m3_t Vx = O3;
    Vx.m12 = -v.z;  Vx.m13 = v.y;
    Vx.m21 = v.z;   Vx.m23 = -v.x;
    Vx.m31 = -v.y;  Vx.m32 = v.x;
    return Vx;
}

/**
 * @brief Convert 3D skew symmetric matrix to 3D vector, the inverse precedure
    of function v3_askew().
 * @param[in] m     3D skew symmetric matrix
 * @return correspoding 3D vector
 * @see v3_askew()
 */
inline extern v3_t m3_iaskew(m3_t m)
{
    v3_t v;
    v.x = (m.m32 - m.m23)/2.0;
    v.y = (m.m13 - m.m31)/2.0;
    v.z = (m.m21 - m.m12)/2.0;
    return v;
}
/**
 * @brief 3D matrix transposition
 * @param[in] A Input 3D matrix
 * @return transposition matrix of A
 */
inline extern m3_t m3_T(m3_t A)
{
    m3_t dcm;
    dcm.m11 = A.m11; dcm.m12 = A.m21; dcm.m13 = A.m31;
    dcm.m21 = A.m12; dcm.m22 = A.m22; dcm.m23 = A.m32;
    dcm.m31 = A.m13; dcm.m32 = A.m23; dcm.m33 = A.m33;
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
 * @param[in] alpha Input number
 * @param[in] A     Input 3D matrix
 * @return Scalar muliplication result of s and A ( s x v )
 */
extern m3_t m3_scalar(double alpha, m3_t A)
{
    m3_t mat;
    mat.m11 = alpha * A.m11; mat.m12 = alpha * A.m12; mat.m13 = alpha * A.m13;
    mat.m21 = alpha * A.m21; mat.m22 = alpha * A.m22; mat.m23 = alpha * A.m23;
    mat.m31 = alpha * A.m31; mat.m32 = alpha * A.m32; mat.m33 = alpha * A.m33;
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
    C.x = A.m11 * B.x + A.m12 * B.y + A.m13 * B.z;
    C.y = A.m21 * B.x + A.m22 * B.y + A.m23 * B.z;
    C.z = A.m31 * B.x + A.m32 * B.y + A.m33 * B.z;
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
    const double *pA = (const double *)A;
    const double *pB = (const double *)B;
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
    m3_t B = m3_mul(m3_T(*A),*A);
    /* Jacobi iteration method to solve eigenvalue and eigenvector of B */
    const double thres = 1E-15;
    m3_t eVal = B, eVec = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
    double afa, sa, ca;
    m3_t Upq;
    for(int i = 0; i < 100; i++){
        double d12 = fabs(eVal.m12), d13 = fabs(eVal.m13), d23 = fabs(eVal.m23);
        if(d12 > thres && d12 >= d13 && d12 >= d23) {
            afa = atan2(2*eVal.m12, eVal.m11 - eVal.m22)/2.0;
            sa = sin(afa); ca = cos(afa);
            Upq = (m3_t){ca,-sa,0.0, sa,ca,0.0, 0.0, 0.0, 1.0};
        }
        else if(d13 > thres && d13 >= d23 && d13 >= d12){
            afa = atan2(2*eVal.m13, eVal.m33 - eVal.m11)/2.0;
            sa = sin(afa); ca = cos(afa);
            Upq = (m3_t){ca, 0.0, sa, 0.0, 1.0, 0.0, -sa, 0.0, ca};
        }
        else if(d23 > thres && d23 >= d12 && d23 >= d13){
            afa = atan2(2*eVal.m23, eVal.m22 - eVal.m33)/2.0;
            sa = sin(afa); ca = cos(afa);
            Upq = (m3_t){1.0, 0.0, 0.0, 0.0, ca, -sa, 0.0, sa, ca};
        }else{
            break;
        }
        eVec = m3_mul(eVec,Upq);
        eVal = m3_mul(m3_T(eVec),m3_mul(B,eVec));
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

/**
 * @brief copy double array elements to 3D vector
 * @param[out] v3_dest  3D vector
 * @param[in]  src      first pointer of source double array
 */
extern void v3_copy(v3_t *v3_dest, const double *src)
{
    v3_dest->x = src[0]; v3_dest->y = src[1]; v3_dest->z = src[2];
}

/**
 * @brief copy double matrix(column-major) elements to 3D matrix
 * @param[out]  m3_dest 3D matrix
 * @param[in]   src     first pointer of source double matrix
 * @param[in]   n       row number of source double matrix
 */
extern void m3_copy(m3_t *m3_dest, const double *src, int n)
{
    double *imat = (double *)m3_dest;
    for(int j = 0; j < 3; ++j)
        for(int i = 0; i < 3; ++i)
            imat[i*3+j] =  src[j*n+i];
}

/**
 * @brief paste matrix to double array
 * @param[out]  dest    first pointer of double array
 * @param[in]   v3      3D vector
 * @see m3_copy() v3_copy() m3_paste()
 */
extern void v3_paste(double *dest, const v3_t *v3)
{
    dest[0] = v3->x; dest[1] = v3->y; dest[2] = v3->z;
}

/**
 * @brief paste matrix to double array matrix
 * @param[out]  dest    double array first address
 * @param[in]   n       matrix row number(n >= 3)
 * @param[in]   m3      3D matrix
 * @see m3_copy() v3_copy() v3_paste()
 */
extern void m3_paste(double *dest, int n, const m3_t *m3)
{
    const double *imat = (const double *)m3;
    for(int j = 0; j < 3; ++j)
        for(int i = 0; i < 3; ++i)
            dest[j*n+i] = imat[i*3+j];
}

earth_t wgs84 = {
    .wie = 7.292115E-5, /**< WGS84 earth parameters */
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
    mat->m11 = 0.0; mat->m12 = -v3->z; mat->m13 = v3->y;
    mat->m21 = v3->z; mat->m22 = 0.0; mat->m23 = -v3->x;
    mat->m31 = -v3->y; mat->m32 = v3->x; mat->m33 = 0.0;
    return 0;
}

#ifndef RTKLIB
/**
 * @brief Convert calendar day/time to gtime_t struct
 * @param[in] ep    day/time {year,month,day,hour,min,sec}
 * @return gtime_t struct
 * @note proper in 1970-2037 or 1970-2099 (64bit time_t)
 */
extern gtime_t epoch2time(const double* ep)
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
extern gtime_t gpst2time(int week, double sec)
{
    const double gpst0[] = { 1980, 1, 6, 0, 0, 0 }; /* gps time reference */
    gtime_t t = epoch2time(gpst0);

    if (sec < -1E9 || 1E9 < sec)
        sec = 0.0;
    t.time += (time_t)86400 * 7 * week + (int)sec;
    t.sec = sec - (int)sec;
    return t;
}

/**
 * @brief convert gtime_t struct to calendar day/time
 * @param[in]  t    gtime_t struct
 * @param[out] ep   day/time {year,month,day,hour,min,sec}
 * @note Proper in 1970-2037 or 1970-2099 (64bit time_t)
 */
extern void time2epoch(gtime_t t, double* ep)
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
 * @brief convert gtime_t struct to week and tow in gps time
 * @param[in]  t        gtime_t struct
 * @param[out] week     week number in gps time (NULL: no output)
 * @return time of week in gps time (s)
 */
extern double time2gpst(gtime_t t, int *week)
{
    const double gpst0[] = { 1980, 1, 6, 0, 0, 0 }; /* gps time reference */
    gtime_t t0 = epoch2time(gpst0);
    time_t sec=t.time-t0.time;
    int w=(int)(sec/(86400*7));

    if (week) *week=w;
    return (double)(sec-(double)w*86400*7)+t.sec;
}

/**
 * @brief time differenece between two gtime_t struct
 * @param[in] t1    First gtime_t struct
 * @param[in] t2    Second gtime_t struct
 * @return seconds of (t1 - t2)
 */
extern double timediff(gtime_t t1, gtime_t t2)
{
    return difftime(t1.time, t2.time) + t1.sec - t2.sec;
}

/**
 * @brief add time to gtime_t struct
 * @param[in] t     gtime_t struct
 * @param[in] sec   time to add[s]
 * @return gtime_t struct(t+sec)
 */
extern gtime_t timeadd(gtime_t t, double sec)
{
    double tt;

    t.sec+=sec; tt=floor(t.sec); t.time+=(int)tt; t.sec-=tt;
    return t;
}
#endif /* RTKLIB */

/**
 * @brief Convert Euler attitude to Quaternion attitude(Eab => Qab)
 * @param[in] euler Input Euler attitude
 * @param[out] quat Ouput Quaternion attitude
 * @return 0: OK
 */
extern int euler2quat(const v3_t* euler, quat_t* quat)
{
    double si = sin(euler->x / 2), ci = cos(euler->x / 2);
    double sj = sin(euler->y / 2), cj = cos(euler->y / 2);
    double sk = sin(euler->z / 2), ck = cos(euler->z / 2);
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
    euler->x = atan2(2 * (-quat->q0 * quat->q1 + quat->q2 * quat->q3),
        1 - 2 * quat->q1 * quat->q1 - 2 * quat->q2 * quat->q2);
    euler->y = asin(2 * (-quat->q0 * quat->q2 - quat->q1 * quat->q3));
    euler->z = atan2(2 * (-quat->q0 * quat->q3 + quat->q1 * quat->q2),
        1 - 2 * quat->q2 * quat->q2 - 2 * quat->q3 * quat->q3);

    if (euler->x <= -PI) /* Limit Roll Angle to (-pi,pi] */
        euler->x += 2 * PI;
    else if (euler->x > PI)
        euler->x -= 2 * PI;
    if (euler->z < 0) /* Limit Heading Angle to [0,2pi) */
        euler->z += 2 * PI;
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
    euler->x = atan2(dcm->m23, dcm->m33);
    euler->y = -asin(dcm->m13);
    euler->z = atan2(dcm->m12, dcm->m11);

    if (euler->x <= -PI) /* Limit Roll Angle to (-pi,pi] */
        euler->x += 2 * PI;
    else if (euler->x > PI)
        euler->x -= 2 * PI;
    if (euler->z < 0) /* Limit Heading Angle to [0,2pi) */
        euler->z += 2 * PI;
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
    double sin_phi = sin(euler->x);
    double cos_phi = cos(euler->x);
    double sin_theta = sin(euler->y);
    double cos_theta = cos(euler->y);
    double sin_psi = sin(euler->z);
    double cos_psi = cos(euler->z);

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
    dcm->m11 = q11 + q22 - q33 - q44; dcm->m12 = 2 * (q23 - q14);
    dcm->m13 = 2 * (q24 + q13); dcm->m21 = 2 * (q23 + q14);
    dcm->m22 = q11 - q22 + q33 - q44; dcm->m23 = 2 * (q34 - q12);
    dcm->m31 = 2 * (q24 - q13); dcm->m32 = 2 * (q34 + q12);
    dcm->m33 = q11 - q22 - q33 + q44;
    return 0;
}

/**
 * @brief Convert Euler attitude(roll, pitch, yaw) to DCM
 * @param[in]   att     Euler attitude(roll, pitch, yaw), Enb[rad]
 * @param[out]  dcm     DCM, Cbn
 * @see euler2dcm() dcm2att()
 * @return 0: OK
 */
extern int att2dcm(const v3_t *att, m3_t *dcm)
{
    euler2dcm(att, dcm);
    *dcm = m3_T(*dcm);
    return 0;
}

/**
 * @brief Convert DCM to Euler attitude(roll, pitch, yaw)
 * @param[in]   dcm     DCM, Cbn
 * @param[out]  att     Euler attitude[roll, pitch, yaw], Enb[rad]
 * @see dcm2euler() att2dcm()
 * @return 0: OK
 */
extern int dcm2att(const m3_t *dcm, v3_t *att)
{
    m3_t M = m3_T(*dcm);
    dcm2euler(&M, att);
    return 0;
}

/**
 * @brief Convert Euler attitude to Quaternion
 * @param[in]   att     Euler attitude(roll, pitch, yaw), Enb[rad]
 * @param[out]  quat    Quaternion, Qbn
 * @see euler2quat() quat2att()
 * @return O: OK
 */
extern int att2quat(const v3_t *att, quat_t *quat)
{
    euler2quat(att, quat);
    quat_conj(quat);
    return 0;
}

/**
 * @brief Convert Quaternion to Euler attitude
 * @param[in]   quat    Quaternion, Qbn
 * @param[out]  att     Euler attitude(roll, pitch, yaw), Enb[rad]
 * @see quat2euler() att2quat()
 * @return 0: OK
 */
extern int quat2att(const quat_t *quat, v3_t *att)
{
    quat_t Q = *quat; quat_conj(&Q);
    quat2euler(&Q, att);
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
    quat->q0 *= nq; quat->q1 *= nq; quat->q2 *= nq; quat->q3 *= nq;
    if (quat->q0 < 0) {
        quat->q0 = -quat->q0; quat->q1 = -quat->q1;
        quat->q2 = -quat->q2; quat->q3 = -quat->q3;
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
    quat->q1 = -quat->q1; quat->q2 = -quat->q2; quat->q3 = -quat->q3;
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
    const double *pP = (const double *)P;
    const double *pQ = (const double *)Q;
    for(int i = 0; i < 4; ++i) {
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
        = dtheta->x * dtheta->x + dtheta->y * dtheta->y + dtheta->z * dtheta->z;
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
    quat->q1 = f * dtheta->x;
    quat->q2 = f * dtheta->y;
    quat->q3 = f * dtheta->z;
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
 * @brief Convert Quaternion attitude to Rotation vector(Angular increment)
 * @param[in]  quat     Quaternion attitude trasnformation
 * @param[out] dtheta   Rotation vector(Angular increment) [rad]
 * @return 0: OK
 * @see dcm2rv()
 * @note  q = [ cos(|rv|/2); sin(|rv|/2)/|rv|*rv ]
 */
extern int quat2rv(const quat_t *quat, v3_t *dtheta)
{
    double n2 = acos(fabs(quat->q0));
    double k;
    if(n2 > 1e-40) k = 2 * n2 / sin(n2); else k = 2.0;
    if(quat->q0 < 0.0) k = -k;
    dtheta->x = k * quat->q1;
    dtheta->y = k * quat->q2;
    dtheta->z = k * quat->q3;
    return 0;
}

/**
 * @brief Convert DCM attitude to Rotation vector(Angular increment)
 * @param[in]   dcm     DCM attitude transformation
 * @param[out]  dtheta  Rotation vector(Angular increment) [rad]
 * @see quat2rv()
 * @return 0: OK
 * @note
 *      DCM = I + sin(|rv|)/|rv|*(rvx) + [1-cos(|rv|)]/|rv|^2*(rvx)^2
 *      rvx: askew matrix of rv
 */
extern int dcm2rv(const m3_t *dcm, v3_t *dtheta)
{

    quat_t quat; dcm2quat(dcm, &quat);
    quat2rv(&quat, dtheta);
    /* TODO: maybe need to iterate to get more accurate solution */
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
    qtmp.q0 = -quat.q1 * vec.x - quat.q2 * vec.y - quat.q3 * vec.z;
    qtmp.q1 = quat.q0 * vec.x + quat.q2 * vec.z - quat.q3 * vec.y;
    qtmp.q2 = quat.q0 * vec.y + quat.q3 * vec.x - quat.q1 * vec.z;
    qtmp.q3 = quat.q0 * vec.z + quat.q1 * vec.y - quat.q2 * vec.x;
    vtmp.x = -qtmp.q0 * quat.q1 + qtmp.q1 * quat.q0 - qtmp.q2 * quat.q3
        + qtmp.q3 * quat.q2;
    vtmp.y = -qtmp.q0 * quat.q2 + qtmp.q2 * quat.q0 - qtmp.q3 * quat.q1
        + qtmp.q1 * quat.q3;
    vtmp.z = -qtmp.q0 * quat.q3 + qtmp.q3 * quat.q0 - qtmp.q1 * quat.q2
        + qtmp.q2 * quat.q1;
    return vtmp;
}

/**
 * @brief form the DCM transform from e-frame to n-frame(Cbe)
 * @param[in] lat   Current postion's latitude[rad]
 * @param[in] lon   Current postion's longitude[rad]
 * @return Cen, DCM attitude transform from e-frame to n-frame
 */
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

/**
 * @brief Convert n-frame(NED) position/velocity/attitude to e-frame(ECEF)
 * @param[in,out] pos   Input (lat,lon,hgt)[rad,m] / Output ECEF position xyz[m]
 * @param[in,out] vel   Input velocity(vN, vE, vD)[m/s] / Output velocity ECEF[m/s]
 * @param[in,out] dcm   Input Cbn attitude / Ouput Cbe attitude
 * @return 0: OK
 * @see ecef2ned()
 * @note  vel and dcm could set to NULL pointer when do NOT interested.
 */
extern int ned2ecef(v3_t* pos, v3_t* vel, m3_t* dcm)
{
    double lat = pos->x, lon = pos->y, hgt = pos->z;
    double coslat = cos(lat), sinlat = sin(lat);
    double coslon = cos(lon), sinlon = sin(lon);

    double tmp = wgs84.e * sinlat;
    double Re = wgs84.R0 / sqrt(1 - SQR(tmp));

    pos->x = (Re + hgt) * coslat * coslon;
    pos->y = (Re + hgt) * coslat * sinlon;
    pos->z = ((1 - wgs84.e * wgs84.e) * Re + hgt) * sinlat;

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
            *vel = m3_mul_v3(Cne, *vel);    /* Veb_n => Veb_e */
        if (dcm != NULL)
            *dcm = m3_mul(Cne, *dcm);       /* Cb_n => Cb_e */
    }
    return 0;
}

/**
 * @brief Convert e-frame(ECEF) position/velocity/attitude to n-frame(NED)
 * @param[in,out] pos   Output ECEF position xyz[m] / Input (lat,lon,hgt)[rad,m]
 * @param[in,out] vel   Output velocity ECEF[m/s] / Input velocity(vN, vE, vD)[m/s]
 * @param[in,out] dcm   Ouput Cbe attitude / Input Cbn attitude
 * @see ned2ecef()
 * @return 0: OK
 * @note vel and dcm could set to NULL pointer when do NOT interested.
 *      Ref: Paul D. Groves, Principles of GNSS, Inertial, and Multisensor
 *          Integrated Navigation Systems(2nd Edition), 2013, C.29-C.38
 */
extern int ecef2ned(v3_t* pos, v3_t* vel, m3_t* dcm)
{
    double lon = atan2(pos->y, pos->x);

    double e2 = wgs84.e * wgs84.e;
    double k1 = sqrt(1.0 - e2) * fabs(pos->z);
    double k2 = e2 * wgs84.R0;
    double beta = sqrt(pos->x * pos->x + pos->y * pos->y);
    double E = (k1 - k2) / beta, F = (k1 + k2) / beta;
    double P = 4.0 / 3.0 * (E * F + 1.0);
    double Q = 2.0 * (E * E - F * F);
    double D = P * P * P + Q * Q;
    double V = pow(sqrt(D) - Q, 1.0 / 3.0) - pow(sqrt(D) + Q, 1.0 / 3.0);
    double G = 0.5 * (sqrt(E * E + V) + E);
    double T = sqrt(G * G + (F - V * G) / (2.0 * G - E)) - G;
    double signz = pos->z > 0 ? 1.0 : -1.0;
    double lat = signz * atan((1 - T * T) / (2 * T * sqrt(1 - e2)));

    double coslat = cos(lat), sinlat = sin(lat);
    double hgt = (beta - wgs84.R0 * T) * coslat
        + (pos->z - signz * wgs84.R0 * sqrt(1.0 - e2)) * sinlat;

    pos->x = lat; pos->y = lon; pos->z = hgt;

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

/**
 * @brief Convert Var-covariance matrix from ned frame to ecef frame
 * @param[in]      pos  Current geodetic position, BLH[rad,m]
 * @param[in,out]  Qpos pos covariance matrix, Q_ned to Q_xyz[m^2]
 * @param[in,out]  Qvel vel covariance matrix, Q_ned to Q_xyz[m^2/s^2]
 * @param[in,out]  Qatt att covariance matrix, Q_ned to Q_xyz[rad^2]
 * @return 0: OK
 * @note Qpos, Qvel and Qatt could be NULL pointer when NOT interested
 */
int ned2ecefQ(const v3_t *pos, m3_t *Qpos, m3_t *Qvel, m3_t *Qatt)
{
    m3_t Cen = formCen_ned(pos->x, pos->y);
    m3_t Cne = m3_T(Cen);
    if(Qpos != NULL) *Qpos = m3_mul(Cne, m3_mul(*Qpos, Cen));
    if(Qvel != NULL) *Qvel = m3_mul(Cne, m3_mul(*Qvel, Cen));
    if(Qatt != NULL) *Qatt = m3_mul(Cne, m3_mul(*Qatt, Cen));
    return 0;
}

/**
 * @brief Convert var-covariance matrix from ecef frame to ned frame
 * @param[in]       xyz  Current ecef position, xyz[m]
 * @param[in,out]   Qxyz pos covaraince matrix, Q_xyz to Q_ned[m^2]
 * @param[in,out]   Qvel vel covariance matrix, Q_xyz to Q_ned[m^2/s^2]
 * @param[in,out]   Qatt att covariance matrix, Q_xyz to Q_ned[rad^2]
 * @return 0: OK
 * @note Qpos, Qvel and Qatt could be NULL pointer when NOT interested
 */
int ecef2nedQ(const v3_t *xyz, m3_t *Qxyz, m3_t *Qvel, m3_t *Qatt)
{
    v3_t pos = *xyz;
    ecef2ned(&pos, NULL, NULL);
    m3_t Cen = formCen_ned(pos.x, pos.y);
    m3_t Cne = m3_T(Cen);
    if(Qxyz != NULL) *Qxyz = m3_mul(Cen, m3_mul(*Qxyz, Cne));
    if(Qvel != NULL) *Qvel = m3_mul(Cen, m3_mul(*Qvel, Cne));
    if(Qatt != NULL) *Qatt = m3_mul(Cen, m3_mul(*Qatt, Cne));
    return 0;
}

/**
 * @brief Convert attitude angle(Enb) to DCM transform from b-frame to
 *      e-frame(Cbe)
 * @param[in] pos   Current geodetic position, BLH[rad,m]
 * @param[in] Enb   Current attitue angle(roll, pitch, yaw)[rad]
 * @return DCM transform from b-frame to e-frame, Cbe
 */
m3_t att2Cbe(const v3_t *pos, const v3_t *Enb)
{
    m3_t Cbn; att2dcm(Enb, &Cbn);
    m3_t Cen = formCen_ned(pos->x, pos->y);
    return m3_mul(m3_T(Cen), Cbn);
}

/**
 * @brief Convert attitude angle(Enb) to quaternion transform from b-frame to
 *      e-frame(Qbe)
 * @param[in] pos   Current geodetic position, BLH[rad,m]
 * @param[in] Enb   Current attitue angle(roll, pitch, yaw)[rad]
 * @return quaternion transform from b-frame to e-frame, Qbe
 */
quat_t att2Qbe(const v3_t *pos, const v3_t *Enb)
{
    m3_t Cbe = att2Cbe(pos, Enb);
    quat_t Qbe; dcm2quat(&Cbe, &Qbe);
    return Qbe;
}

/**
 * @brief Convert attitude angle(Enb) to Euler transform from b-frame to
 *      e-frame(Ebe)
 * @param[in] pos   Current geodetic position, BLH[rad,m]
 * @param[in] Enb   Current attitue angle(roll, pitch, yaw)[rad]
 * @return Ebe, Euler transform from b-frame to e-frame
 */
v3_t att2Ebe(const v3_t *pos, const v3_t *Enb)
{
    m3_t Cbe = att2Cbe(pos, Enb);
    v3_t Ebe; dcm2euler(&Cbe, &Ebe);
    return Ebe;
}

/**
 * @brief Convert quaternion transform from b-frame to e-frame(Qbe) to attitude
 *      angle(Enb)
 * @param[in] xyz   Current geodetic position, BLH[rad,m]
 * @param[in] Cbe   Current DCM transform from b-frame to e-frame, Cbe
 * @return attitude angle(roll, pitch, yaw)[Enb]
 */
v3_t Cbe2att(const v3_t *xyz, const m3_t *Cbe)
{
    v3_t pos = *xyz; m3_t Cbn = *Cbe;
    ecef2ned(&pos, NULL, &Cbn);
    v3_t att; dcm2att(&Cbn, &att);
    return att;
}

/**
 * @brief Convert quaternion transform from b-frame to e-frame(Qbe) to attitude
 *      angle(Enb)
 * @param[in] xyz   Current geodetic position, BLH[rad,m]
 * @param[in] Qbe   Current quaternion transform from b-frame to e-frame, Qbe
 * @return Enb, attitude angle(roll, pitch, yaw) [rad]
 */
v3_t Qbe2att(const v3_t *xyz, const quat_t *Qbe)
{
    m3_t Cbe; quat2dcm(Qbe, &Cbe);
    return Cbe2att(xyz, &Cbe);
}

/**
 * @brief Convert Euler transform from b-frame to e-frame(Ebe) to attitude
 *      angle(Enb)
 * @param[in] xyz   Current geodetic position, BLH[rad,m]
 * @param[in] Ebe   Current Euler transform from b-frame to e-frame, Ebe[rad]
 * @return Enb, attitude angle(roll, pitch, yaw)[rad]
 */
v3_t Ebe2att(const v3_t *xyz, const v3_t *Ebe)
{
    m3_t Cbe; euler2dcm(Ebe, &Cbe);
    return Cbe2att(xyz, &Cbe);
}

/**
 * @brief Euler attitude add operation, the same as left-multiply DCM, and
 *      similar to vector add operation by superscript&subscript
 * @param[in] Eab   Euler attitude transform from a-frame to b-frame[rad]
 * @param[in] Ebc   Euler attitude transform from b-frame to c-frame[rad]
 * @return Eac, Euler attitude transofrom from a-frame to c-frame[rad]
 */
v3_t euler_add(v3_t Eab, v3_t Ebc)
{
    m3_t Cab; euler2dcm(&Eab, &Cab);
    m3_t Cbc; euler2dcm(&Ebc, &Cbc);
    m3_t Cac = m3_mul(Cbc, Cab);
    v3_t Eac; dcm2euler(&Cac, &Eac);
    return Eac;
}

/**
 * @brief Eulear attitue delete operation, get two Euler angel difference, and
 *      similar to vector delete operation by superscript&subscript
 * @param[in] Eac   Euler attitude transform from a-frame to c-frame
 * @param[in] Eab   Euler attitude transform from a-frame to b-frame
 * @return Eac, Euler attitude transform from a-frame to c-frame
 */
v3_t euler_del(v3_t Eac, v3_t Eab)
{
    m3_t Cab; euler2dcm(&Eab, &Cab);
    m3_t Cac; euler2dcm(&Eac, &Cac);
    m3_t Cbc = m3_mul(Cac, m3_T(Cab));
    v3_t Ebc; dcm2euler(&Cbc, &Ebc);
    return Ebc;
}

/**
 * @brief Euler angle attitude add a samll error angle
 * @param[in,out]   E       Euler angle attitude[rad]
 * @param[in]       phi     Small error angle at three axis[rad]
 * @return (I + phi x)E
 * @note phi must keep small, otherwise this function would lead to large
 *  calcuate error
 */
extern int euler_addphi(v3_t *E, const v3_t *phi)
{
    if(v3_norm(*phi) > 1.0*DEG2RAD){
        LOG_WARN("euler_addphi: error angel may be too large, %f %f %f",
                 phi->x, phi->y, phi->z);
    }
    m3_t dcm, phix;
    euler2dcm(E, &dcm);
    asymmetric_mat(phi, &phix);
    phix.m11 += 1.0; phix.m22 += 1.0; phix.m33 += 1.0;
    dcm = m3_mul(phix, dcm);
    dcm2euler(&dcm, E);
    return 0;
}

/**
 * @brief Euler angle attitude add a samll error angle(minus)
 * @param[in,out]   E       Euler angle attitude[rad]
 * @param[in]       phi     Small error angle at three axis[rad]
 * @return (I - phi x)E
 * @note phi must keep small, otherwise this function would lead to large
 *  calcuate error
 */
extern int euler_delphi(v3_t *E, const v3_t *phi)
{
    if(v3_norm(*phi) > 1.0*DEG2RAD){
        LOG_WARN("eluer_delphi: error angel may be too large: %f %f %f",
                 phi->x, phi->y, phi->z);
    }
    m3_t dcm, phix;
    euler2dcm(E, &dcm);
    asymmetric_mat(phi, &phix);
    phix = m3_scalar(-1.0, phix);
    phix.m11 += 1.0; phix.m22 += 1.0; phix.m33 += 1.0;
    dcm = m3_mul(phix, dcm);
    dcm2euler(&dcm, E);
    return 0;
}

/**
 * @brief yaw angel delete, get two yaw difference
 * @param[in] yaw1  Fisrt yaw angel[rad]
 * @param[in] yaw2  Second yaw angel[rad]
 * @return yaw1 - yaw2, angle differnece, range from -PI to PI[rad]
 */
extern double yaw_del(double yaw1, double yaw2)
{
    double dyaw = yaw1 - yaw2;
    while(dyaw > PI)    dyaw -= 2*PI;
    while(dyaw < -PI)   dyaw += 2*PI;
    return dyaw;
}

#ifndef RTKLIB
/**
 * @brief new matrix, allocate memeory of matrix
 * @param[in] n number of rows of matrix
 * @param[in] m number of columns of matrix
 * @return  matrix pointer(if n<=0 or m>=0, return NULL)
 */
extern double *mat(int n, int m)
{
    double *p;

    if (n<=0||m<=0) return NULL;
    if (!(p=(double *)malloc(sizeof(double)*n*m))) {
        LOG_FATAL("matrix memory allocation error: n=%d,m=%d",n,m);
    }
    return p;
}
/**
 * @brief new interger matrix, allocate memory of integer matrix
 * @param[in] n number of rows of matrix
 * @param[in] m number of columns of matrix
 * @return matrix pointer(if n<=0 or m<=0, return NULL)
 */
extern int *imat(int n, int m)
{
    int *p;

    if (n<=0||m<=0) return NULL;
    if (!(p=(int *)malloc(sizeof(int)*n*m))) {
        LOG_FATAL("integer matrix memory allocation error: n=%d,m=%d",n,m);
    }
    return p;
}
/**
 * @brief generate new zero matrix
 * @param[in] n nunmber of rows of matrix
 * @param[in] m number of colums of matrix
 * @return matrix pointer(if n<=0 or m<=0, return NULL)
 */
extern double *zeros(int n, int m)
{
    double *p;
#if NOCALLOC
    if ((p=mat(n,m))) for (n=n*m-1;n>=0;n--) p[n]=0.0;
#else
    if (n<=0||m<=0) return NULL;
    if (!(p=(double *)calloc(sizeof(double),n*m))) {
        LOG_FATAL("matrix memory allocation error: n=%d,m=%d",n,m);
    }
#endif
    return p;
}
/**
 * @brief generate new identify matrix
 * @param n number of rows and columns of matrix
 * @return  matrix pointer(if n<=0, return NULL)
 */
extern double *eye(int n)
{
    double *p;
    int i;

    if ((p=zeros(n,n))) for (i=0;i<n;i++) p[i+i*n]=1.0;
    return p;
}
/**
 * @brief inner product of vector
 * @param[in] a vector a(nx1)
 * @param[in] b vector b(nx1)
 * @param[in] n size of vector a,b
 * @return a'*b
 */
extern double dot(const double *a, const double *b, int n)
{
    double c=0.0;

    while (--n>=0) c+=a[n]*b[n];
    return c;
}
/**
 * @brief euclid norm of vector
 * @param[in] a vector a(nx1)
 * @param[in] n size of vector a
 * @return || a ||
 */
extern double norm(const double *a, int n)
{
    return sqrt(dot(a,a,n));
}

/**
 * @brief outer product of 3d vectors
 * @param[in] a     vector a(3x1)
 * @param[in] b     vector b(3x1)
 * @param[out] c    outer product result of a,b (3x1)
 */
extern void cross3(const double *a, const double *b, double *c)
{
    c[0]=a[1]*b[2]-a[2]*b[1];
    c[1]=a[2]*b[0]-a[0]*b[2];
    c[2]=a[0]*b[1]-a[1]*b[0];
}

/**
 * @brief normailize 3D vector
 * @param[in] 	a 	vector a (3x1)
 * @param[out] 	b	normlized vector(3x1) || b || = 1
 * @return status (1:ok, 0:error)
 */
extern int normv3(const double *a, double *b)
{
    double r;
    if ((r=norm(a,3))<=0.0) return 0;
    b[0]=a[0]/r;
    b[1]=a[1]/r;
    b[2]=a[2]/r;
    return 1;
}

/**
 * @brief copy matrix
 * @param[out] 	A	destination matrix A (n x m)
 * @param[in] 	B	source matrix B (n x m)
 * @param[in] 	n	number of rows of matrix
 * @param[in] 	m	number of columns of matrix
 */
extern void matcpy(double *A, const double *B, int n, int m)
{
    memcpy(A,B,sizeof(double)*n*m);
}
/* matrix routines -----------------------------------------------------------*/

#ifdef LAPACK /* with LAPACK/BLAS or MKL */
/**
 * @brief multiply matrix (wrapper of blas degmm)
 * 		( C = alpha * A * B + beta * C )
 * @param[in] 	tr		transpose flas ("N": normal, "T": transpose)
 * @param[in]  	n		size of (transposed) matrix A, B
 * @param[in] 	k		size of (transposed) matrix A, B
 * @param[in] 	m		size of (transposed) matrix A, B
 * @param[in] 	alpha	alpha
 * @param[in] 	A		(transposd) matrix A (n x m)
 * @param[in] 	B		(transposd) matrix B (m x k)
 * @param[in] 	beta	beta
 * @param[in,out] 	C	matrix C (n x k)
 */
extern void matmul(const char *tr, int n, int k, int m, double alpha,
                   const double *A, const double *B, double beta, double *C)
{
    int lda=tr[0]=='T'?m:n,ldb=tr[1]=='T'?k:m;

    dgemm_((char *)tr,(char *)tr+1,&n,&k,&m,&alpha,(double *)A,&lda,(double *)B,
           &ldb,&beta,C,&n);
}

/**
 * @brief inverse of matrix, ( A = A^-1 )
 * @param[in,out] 	A	matrix (n x n)
 * @param[in] 		n	size of matrix A
 * @return status (0: ok, >0: error)
 */
extern int matinv(double *A, int n)
{
    double *work;
    int info,lwork=n*16,*ipiv=imat(n,1);

    work=mat(lwork,1);
    dgetrf_(&n,&n,A,&n,ipiv,&info);
    if (!info) dgetri_(&n,A,&n,ipiv,work,&lwork,&info);
    free(ipiv); free(work);
    return info;
}

/**
 * @brief solve linear equation ( X = A\Y or X = A'\Y )
 * @param[in] 	tr	transpose flag( "N": normal, "T": transpose)
 * @param[in] 	A	input matrix A (n x n)
 * @param[in] 	Y	input matrix Y (n x n)
 * @param[in] 	n	rows of matrix A, Y
 * @param[in] 	m	columns of matrix A, Y
 * @param[out] 	X	X = A\Y or X = A'\Y (n x m)
 * @return status (0:ok, >0: error)
 * @note matrix stored by column-major order(fortran convertion), X can be same
 * 		as Y
 */
extern int solve(const char *tr, const double *A, const double *Y, int n,
                 int m, double *X)
{
    double *B=mat(n,n);
    int info,*ipiv=imat(n,1);

    matcpy(B,A,n,n);
    matcpy(X,Y,n,m);
    dgetrf_(&n,&n,B,&n,ipiv,&info);
    if (!info) dgetrs_((char *)tr,&n,&m,B,&n,ipiv,X,&n,&info);
    free(ipiv); free(B);
    return info;
}

#else /* without LAPACK/BLAS or MKL */

/* multiply matrix -----------------------------------------------------------*/
extern void matmul(const char *tr, int n, int k, int m, double alpha,
                   const double *A, const double *B, double beta, double *C)
{
    double d;
    int i,j,x,f=tr[0]=='N'?(tr[1]=='N'?1:2):(tr[1]=='N'?3:4);

    for (i=0;i<n;i++) for (j=0;j<k;j++) {
        d=0.0;
        switch (f) {
            case 1: for (x=0;x<m;x++) d+=A[i+x*n]*B[x+j*m]; break;
            case 2: for (x=0;x<m;x++) d+=A[i+x*n]*B[j+x*k]; break;
            case 3: for (x=0;x<m;x++) d+=A[x+i*m]*B[x+j*m]; break;
            case 4: for (x=0;x<m;x++) d+=A[x+i*m]*B[j+x*k]; break;
        }
        if (beta==0.0) C[i+j*n]=alpha*d; else C[i+j*n]=alpha*d+beta*C[i+j*n];
    }
}
/* LU decomposition ----------------------------------------------------------*/
static int ludcmp(double *A, int n, int *indx, double *d)
{
    double big,s,tmp,*vv=mat(n,1);
    int i,imax=0,j,k;

    *d=1.0;
    for (i=0;i<n;i++) {
        big=0.0; for (j=0;j<n;j++) if ((tmp=fabs(A[i+j*n]))>big) big=tmp;
        if (big>0.0) vv[i]=1.0/big; else {free(vv); return -1;}
    }
    for (j=0;j<n;j++) {
        for (i=0;i<j;i++) {
            s=A[i+j*n]; for (k=0;k<i;k++) s-=A[i+k*n]*A[k+j*n]; A[i+j*n]=s;
        }
        big=0.0;
        for (i=j;i<n;i++) {
            s=A[i+j*n]; for (k=0;k<j;k++) s-=A[i+k*n]*A[k+j*n]; A[i+j*n]=s;
            if ((tmp=vv[i]*fabs(s))>=big) {big=tmp; imax=i;}
        }
        if (j!=imax) {
            for (k=0;k<n;k++) {
                tmp=A[imax+k*n]; A[imax+k*n]=A[j+k*n]; A[j+k*n]=tmp;
            }
            *d=-(*d); vv[imax]=vv[j];
        }
        indx[j]=imax;
        if (A[j+j*n]==0.0) {free(vv); return -1;}
        if (j!=n-1) {
            tmp=1.0/A[j+j*n]; for (i=j+1;i<n;i++) A[i+j*n]*=tmp;
        }
    }
    free(vv);
    return 0;
}
/* LU back-substitution ------------------------------------------------------*/
static void lubksb(const double *A, int n, const int *indx, double *b)
{
    double s;
    int i,ii=-1,ip,j;

    for (i=0;i<n;i++) {
        ip=indx[i]; s=b[ip]; b[ip]=b[i];
        if (ii>=0) for (j=ii;j<i;j++) s-=A[i+j*n]*b[j]; else if (s) ii=i;
        b[i]=s;
    }
    for (i=n-1;i>=0;i--) {
        s=b[i]; for (j=i+1;j<n;j++) s-=A[i+j*n]*b[j]; b[i]=s/A[i+i*n];
    }
}
/* inverse of matrix ---------------------------------------------------------*/
extern int matinv(double *A, int n)
{
    double d,*B;
    int i,j,*indx;

    indx=imat(n,1); B=mat(n,n); matcpy(B,A,n,n);
    if (ludcmp(B,n,indx,&d)) {free(indx); free(B); return -1;}
    for (j=0;j<n;j++) {
        for (i=0;i<n;i++) A[i+j*n]=0.0;
        A[j+j*n]=1.0;
        lubksb(B,n,indx,A+j*n);
    }
    free(indx); free(B);
    return 0;
}
/* solve linear equation -----------------------------------------------------*/
extern int solve(const char *tr, const double *A, const double *Y, int n,
                 int m, double *X)
{
    double *B=mat(n,n);
    int info;

    matcpy(B,A,n,n);
    if (!(info=matinv(B,n))) matmul(tr[0]=='N'?"NN":"TN",n,m,n,1.0,B,Y,0.0,X);
    free(B);
    return info;
}
#endif
/* end of matrix routines ----------------------------------------------------*/

/**
 * @brief least square estimation by solving normal equation (x=(A*A')^-1*A*y)
 * @param[in] 	A	transpose of (weighted) design matrix (n x m)
 * @param[in] 	y 	(weighted) measurements (m x 1)
 * @param[in]	n 	number of parameters(n<=m)
 * @param[in] 	m	number of measurements(n<=m)
 * @param[out] 	x 	estmated parameters (n x 1)
 * @param[out] 	Q 	esimated parameters covariance matrix (n x n)
 * @return status (0: ok, >0: error)
 * @note for weighted least square, replace A and y by A*w and w*y (w=W^(1/2))
 *          matirix stored by column-major order (fortran convention)
 */
extern int lsq(const double *A, const double *y, int n, int m, double *x,
               double *Q)
{
    double *Ay;
    int info;

    if (m<n) return -1;
    Ay=mat(n,1);
    matmul("NN",n,1,m,1.0,A,y,0.0,Ay); /* Ay=A*y */
    matmul("NT",n,n,m,1.0,A,A,0.0,Q);  /* Q=A*A' */
    if (!(info=matinv(Q,n))) matmul("NN",n,1,n,1.0,Q,Ay,0.0,x); /* x=Q^-1*Ay */
    free(Ay);
    return info;
}

/**
 * @brief kalman filter
 * @details kalman filter state update as follows:
 *   		K=P*H*(H'*P*H+R)^-1, x=x+K*v, P=(I-K*H')*P
 * @param[in] 		x 	states vector (n x 1)
 * @param[in] 		P 	covariance matrix of states (n x n)
 * @param[in] 		H 	transpose of design matrix (n x m)
 * @param[in] 		v 	innovation (measurement - model) (m x 1)
 * @param[in] 		R 	covariance matrix of measurement error (m x m)
 * @param[in] 		n 	number of states
 * @param[in] 		m	number of measurements
 * @param[out] 		xp 	states vector after update (n x 1)
 * @param[out] 		Pp 	covariance matrix of states after update (n x n)
 * @return status (0:ok,<0:error)
 * @note  : matirix stored by column-major order (fortran convention)
 *          if state x[i]==0.0, not updates state x[i]/P[i+i*n]
 */
static int filter_(const double *x, const double *P, const double *H,
                   const double *v, const double *R, int n, int m,
                   double *xp, double *Pp)
{
    double *F=mat(n,m),*Q=mat(m,m),*K=mat(n,m),*I=eye(n);
    int info;

    matcpy(Q,R,m,m);
    matcpy(xp,x,n,1);
    matmul("NN",n,m,n,1.0,P,H,0.0,F);       /* Q=H'*P*H+R */
    matmul("TN",m,m,n,1.0,H,F,1.0,Q);
    if (!(info=matinv(Q,m))) {
        matmul("NN",n,m,m,1.0,F,Q,0.0,K);   /* K=P*H*Q^-1 */
        matmul("NN",n,1,m,1.0,K,v,1.0,xp);  /* xp=x+K*v */
        matmul("NT",n,n,m,-1.0,K,H,1.0,I);  /* Pp=(I-K*H')*P */
        matmul("NN",n,n,n,1.0,I,P,0.0,Pp);
    }
    free(F); free(Q); free(K); free(I);
    return info;
}

/**
 * @brief kalman filter
 * @details kalman filter state update as follows:
 *   		K=P*H*(H'*P*H+R)^-1, x=x+K*v, P=(I-K*H')*P
 * @param[in,out] 	x 	states vector (n x 1)
 * @param[in,out] 	P 	covariance matrix of states (n x n)
 * @param[in] 		H 	transpose of design matrix (n x m)
 * @param[in] 		v 	innovation (measurement - model) (m x 1)
 * @param[in] 		R 	covariance matrix of measurement error (m x m)
 * @param[in] 		n 	number of states
 * @param[in] 		m	number of measurements
 * @return status (0:ok,<0:error)
 * @note matirix stored by column-major order (fortran convention)
 *       if state x[i]==0.0, not updates state x[i]/P[i+i*n]
 */
extern int filter(double *x, double *P, const double *H, const double *v,
                  const double *R, int n, int m)
{
    double *x_,*xp_,*P_,*Pp_,*H_;
    int i,j,k,info,*ix;

    ix=imat(n,1);
    for(i=k=0; i<n; i++)
        if(x[i]!=0.0&&P[i+i*n]>0.0)
        ix[k++]=i;
    x_=mat(k,1); xp_=mat(k,1); P_=mat(k,k); Pp_=mat(k,k); H_=mat(k,m);
    for (i=0;i<k;i++) {
        x_[i]=x[ix[i]];
        for (j=0;j<k;j++) P_[i+j*k]=P[ix[i]+ix[j]*n];
        for (j=0;j<m;j++) H_[i+j*k]=H[ix[i]+j*n];
    }
    info=filter_(x_,P_,H_,v,R,k,m,xp_,Pp_);
    for (i=0;i<k;i++) {
        x[ix[i]]=xp_[i];
        for (j=0;j<k;j++) P[ix[i]+ix[j]*n]=Pp_[i+j*k];
    }
    free(ix); free(x_); free(xp_); free(P_); free(Pp_); free(H_);
    return info;
}

/**
 * @brief combine forward and backward filters by fixed-interval smoother
 * 	as follows:
*   		xs=Qs*(Qf^-1*xf+Qb^-1*xb), Qs=(Qf^-1+Qb^-1)^-1)
 * @param[in] 	xf 		forward solutions (n x 1)
 * @param[in] 	Qf 		forward solutions covariance matrix (n x n)
 * @param[in] 	xb 		backward solutions (n x 1)
 * @param[in] 	Qb		backward solutions covariance matrix (n x n)
 * @param[in] 	n 		number of solutions
 * @param[out] 	xs 		smoothed solutions (n x 1)
 * @param[out] 	Qs 		smoothed solutions covariance matrix (n x n)
 * @return status (0:ok,0>:error)
 * @note matirix stored by column-major order (fortran convention)
 */
extern int smoother(const double *xf, const double *Qf, const double *xb,
                    const double *Qb, int n, double *xs, double *Qs)
{
    double *invQf=mat(n,n),*invQb=mat(n,n),*xx=mat(n,1);
    int i,info=-1;

    matcpy(invQf,Qf,n,n);
    matcpy(invQb,Qb,n,n);
    if (!matinv(invQf,n)&&!matinv(invQb,n)) {
        for (i=0;i<n*n;i++) Qs[i]=invQf[i]+invQb[i];
        if (!(info=matinv(Qs,n))) {
            matmul("NN",n,1,n,1.0,invQf,xf,0.0,xx);
            matmul("NN",n,1,n,1.0,invQb,xb,1.0,xx);
            matmul("NN",n,1,n,1.0,Qs,xx,0.0,xs);
        }
    }
    free(invQf); free(invQb); free(xx);
    return info;
}

/**
 * @brief print matrix to file
 * @param[in] 	A		matrix A (n x m)
 * @param[in] 	n 		number of rows of A
 * @param[in] 	m 		number of columns of A
 * @param[in] 	p 		total columns, columns under decimal point
 * @param[in] 	q 		total columns, columns under decimal point
 * @param[in] 	fp 		output file pointer
 * @note  matirix stored by column-major order (fortran convention)
 */
static void matfprint(const double A[], int n, int m, int p, int q, FILE *fp)
{
    int i,j;

    for (i=0;i<n;i++) {
        for (j=0;j<m;j++) fprintf(fp," %*.*G",p,q,A[i+j*n]);
        fprintf(fp,"\n");
    }
}

/**
 * @brief print matrix to file
 * @param[in] 	A		matrix A (n x m)
 * @param[in] 	n 		number of rows of A
 * @param[in] 	m 		number of columns of A
 * @param[in] 	p 		total columns, columns under decimal point
 * @param[in] 	q 		total columns, columns under decimal point
 * @param[in] 	fp 		output file pointer
 * @note  matirix stored by column-major order (fortran convention)
 */
extern void matprint(const double A[], int n, int m, int p, int q)
{
    matfprint(A,n,m,p,q,stdout);
}

/**
 * @brief symmetrize matrix
 * @param[in,out] 	A	matrix A (n x )
 * @param[in] 		n 	size of matrix A
 * @return 0: OK
 */
extern int mat_symmetry(double *A, int n)
{
    for(int i =1; i < n; ++i ){
        for(int j = 0; j < i; ++j){
            A[i+j*n] = 0.5 * (A[i+j*n]  + A[j+i*n]);
            A[j+i*n] = A[i+j*n];
        }
    }
    return 0;
}
#endif  /* RTKLIB */

enum LOG_LEVEL FILE_LOG_LEVEL = LEVEL_TRACE;    /**< default file log level */
enum LOG_LEVEL STDOUT_LOG_LEVEL = LEVEL_DEBUG;  /**< default stdout log level */
FILE *LOG_FP = NULL;                /**< default log file pointer(stdout) */
gtime_t LOG_CURTIME = {0, 0.0};     /**< default log time */

/**
 * @brief Open log file, so that all logs could write into.
 * @param[in] file  log file name
 * @return 0: OK, 1: open file error
 */
extern int LOG_OPEN(const char *file)
{
    if((LOG_FP = fopen(file, "w")))
        return 0;
    else{
        printf("ERROR: LOG FILE (%s) failed to open!\n", file);
        LOG_FP = NULL;
        FILE_LOG_LEVEL = LEVEL_OFF;
        return 1;
    }
}

/**
 * @brief Close log file, and then all logs will write into stdout.
 * @return 0: OK, 1: filed to close log file
 */
extern int LOG_CLOSE(void)
{
    if(LOG_FP == NULL) return 1;
    if(!(fclose(LOG_FP))){ LOG_FP = NULL; return 1; }
    else{ printf("ERROR: LOG FILE failed to close!\n"); return 0; }
}

#define LOG_LEVEL_BODY(level) {                                             \
    va_list ap;                                                             \
    int week; double sec;                                                   \
    if(LOG_CURTIME.time == 0.0){                                            \
        week = 0; sec = LOG_CURTIME.sec;                                    \
    }                                                                       \
    else                                                                    \
        sec = time2gpst(LOG_CURTIME, &week);                                \
    if(FILE_LOG_LEVEL >= LEVEL_##level && LOG_FP != NULL){                  \
        if(week != 0 || sec != 0.0)                                         \
            fprintf(LOG_FP, #level"[%i_%.3f]: ", week, sec);                \
        else                                                                \
            fprintf(LOG_FP, #level": ");                                    \
        va_start(ap, format); vfprintf(LOG_FP, format, ap); va_end(ap);     \
        fprintf(LOG_FP, "\n");                                              \
        fflush(LOG_FP);                                                     \
    }                                                                       \
    if(STDOUT_LOG_LEVEL >= LEVEL_##level){                                  \
        if(week != 0 || sec != 0.0)                                         \
            fprintf(stdout, #level"[%i_%.3f]: ", week, sec);                \
        else                                                                \
            fprintf(stdout, #level": ");                                    \
        va_start(ap, format); vfprintf(stdout, format, ap); va_end(ap);     \
        fprintf(stdout, "\n");                                              \
        fflush(stdout);                                                     \
    }                                                                       \
}

/**
 * @brief LOG FATAL message, means that program would NOT work correctly at all,
 *      and if the program do not be terminated, may cause a disaster.
 * @param format    the same as "printf" function pararmeters
 * @note  This LOG will terminate current program with return value 1, and it
 *      will append a line break automatically
 */
extern void LOG_FATAL(const char *format, ...) {
    LOG_LEVEL_BODY(FATAL); LOG_CLOSE(); exit(1);
}
/**
 * @brief LOG ERROR message, meas that program may not work correctly, but in
 *      most situations, program will give a error result.
 * @param format    the same as "printf" function pararmeters
 * @note  the function will will append a line break automatically
 */
extern void LOG_ERROR(const char *format, ...) { LOG_LEVEL_BODY(ERROR); }
/**
 * @brief LOG WARNNING message, meas that program may not work as expected, but
 *      in most situations, program will give a not too bad result.
 * @param format    the same as "printf" function pararmeters
 * @note  the function will will append a line break automatically
 */
extern void LOG_WARN(const char *format, ...) { LOG_LEVEL_BODY(WARN); }
/**
 * @brief LOG INFORMATION message, means that current program situation,
 *      somthing important need to know.
 * @param format    the same as "printf" function pararmeters
 * @note  the function will will append a line break automatically
 */
extern void LOG_INFO(const char *format, ...) { LOG_LEVEL_BODY(INFO); }
/**
 * @brief LOG DEBUG message, more clear message of current program situation,
 *      this message will help you to locate the unexpected behaviors position.
 * @param format    the same as "printf" function pararmeters
 * @note  the function will will append a line break automatically
 */
extern void LOG_DEBUG(const char *format, ...) { LOG_LEVEL_BODY(DEBUG); }
/**
 * @brief LOG TRACE message, the most detailed data trace, in order to solve
 *      data error related problem.
 * @param format    the same as "printf" function pararmeters
 * @note  the function will will append a line break automatically
 */
extern void LOG_TRACE(const char *format, ...) { LOG_LEVEL_BODY(TRACE); }
#undef LOG_LEVEL_BODY

extern void SHOW_PROGRESS(const char *msg, int current, int total)
{
    printf("\r");
    char buff[1024];
    double percent = ((double)current)/((double)total)*100.0;
    sprintf(buff, "%s %i/%i[%3.2f%%]\n", msg, current, total, percent);
    if(strlen(buff) >= 1023){
        LOG_ERROR("SHOW_PROGRESS: string too long, out of range: %i", 1023);
    }
    printf("%s", buff);
    fflush(stdout);
}

/**
 * @brief determinate yaw angel by velocity(velocity yaw angel under n-frame)
 * @param[in]   veb_n   b-frame velocity under ECEF, project to NED frame[m/s]
 * @param[out]  yaw     velocity yaw angel[rad]
 * @param[in]   Qveb_n  velocity covarinace matirx under NED frame[m^2/s^2]
 * @param[out]  Qyaw    velocity yaw angle variance[rad^2]
 * @return 0: OK, 1: yaw angel can not be determiated
 * @note if velocity too small, function will stay remain and return 1
 *      Qveb_n and Qyaw could set to NULL pointer when not interested
 *
 *      ref: CHAI Yanju, OU Jikun, YUAN yunbin et al. The Adaptive Kalman
 *          Filtering for Single Antenna GPS/INS intergrated System with Heading
 *          Angle Constranint by selecting the parameter weights. 2011
 */
extern int vel2yaw(const v3_t *veb_n, double *yaw,
                    const m3_t *Qveb_n, double *Qyaw)
{
    double SQR_vel = SQR(veb_n->x) + SQR(veb_n->y);
    if(sqrt(SQR_vel) < 1.0){
        LOG_DEBUG("vel2yaw:  velocity too small to determinate attitude:(%f %f)",
                 veb_n->x, veb_n->y);
        return 1;
    }
    *yaw = atan2(veb_n->y, veb_n->x);
    if(Qveb_n != NULL && Qyaw != NULL){
        v3_t A = {- veb_n->x / SQR_vel, veb_n->y / SQR_vel, 0.0};
        *Qyaw = v3_mul_rxc(A, m3_mul_v3(*Qveb_n, A));
    }
    return 0;
}

/**
 * @brief calculate yaw difference between imu body velocity direction and
 *      carrier forward direction
 * @param[in] veb_b_x   imu body forward velocity[m/s]
 * @param[in] web_n_z   yaw angular rate[rad/s]
 * @param[in] dL        lever_arm_car.i[m]
 * @return swerve yaw error compensention[rad]
 * @see dyaw_swerve_Q()
 * @note
 *  1. dyaw_swerve will not be zero when carrier direction changing
 *  2. lever_arm_car.i is the distance between imu body reference point and
 *      carrier body refecence lateral axis at carrier forward direction.
 *  ref: Yan Gongming, Wen Jun. teaching materials. P160, 2016.9
 */
extern double dyaw_swerve(double veb_b_x, double web_n_z, double dL)
{
    if(fabs(veb_b_x) < 0.05){
        LOG_DEBUG("dyaw_swerve: velocity too small to calcuate swerve yaw"
                 "(return 0.0)");
        return 0.0;
    }
    double sinphi = dL * web_n_z / veb_b_x;
    if(fabs(sinphi) > 1.0){
        LOG_DEBUG("dyaw_swerve: dyaw_swerve can NOT be calculated(return 0.0)");
        return 0.0;
    }
    return asin(sinphi);
}

/**
 * @brief calcuate the dyaw_swerve variance
 * @param[in] dyaw      dyaw_sserve[rad]
 * @param[in] veb_b_x   imu body forward velocity[m/s]
 * @param[in] web_n_z   yaw angular rate[rad/s]
 * @param[in] dL        lever_arm_car.i[m]
 * @param[in] Qv        veb_b_x variance[m^2/s^2]
 * @param[in] Qw        web_n_z variance[rad^2]
 * @param[in] QL        dL variance[m^2]
 * @return variance of dyaw_swerve[rad^2]
 * @see dyaw_swerve()
 */
extern double dyaw_swerve_Q(double dyaw, double veb_b_x, double web_n_z,
                            double dL, double Qv, double Qw, double QL)
{
    if(fabs(veb_b_x) < 0.05){
        LOG_DEBUG("dyaw_swerve_Q: velocity too small to calculate converance"
                 "(return 999.9)");
        return 999.9;
    }
    double tmp = 1.0 / (veb_b_x  * cos(dyaw));
    v3_t H = { -web_n_z * tmp, -dL * tmp, dL * web_n_z * tmp / veb_b_x};
    m3_t Q = I3; Q.m11 = QL; Q.m22 = Qw; Q.m33 = Qv;
    return v3_mul_rxc(H, m3_mul_v3(Q, H));
}

/**
 * @brief update global variables cfg(cfg_t)
 * @return 0:OK
 * @note this function is use to control the index in cfg variable, when change
 *  cfg setting, this function should be called to fit the problem
 */
extern int updatecfg()
{
    /* configure for kf state vector */
    cfg.nx = 9;
    if(cfg.isx_ba){ cfg.IBA = cfg.nx; cfg.nx += 3; }
    if(cfg.isx_bg){ cfg.IBG = cfg.nx; cfg.nx += 3; }
    if(cfg.isx_kax){ cfg.IKAx = cfg.nx; cfg.nx += 1; }
    if(cfg.isx_kay){ cfg.IKAy = cfg.nx; cfg.nx += 1; }
    if(cfg.isx_kaz){ cfg.IKAz = cfg.nx; cfg.nx += 1; }
    if(cfg.isx_kgx){ cfg.IKGx = cfg.nx; cfg.nx += 1; }
    if(cfg.isx_kgy){ cfg.IKGy = cfg.nx; cfg.nx += 1; }
    if(cfg.isx_kgz){ cfg.IKGz = cfg.nx; cfg.nx += 1; }
    if(cfg.isx_kod){ cfg.IKOD = cfg.nx; cfg.nx += 1; }
    if(cfg.isx_eroll){ cfg.IEROLL = cfg.nx; cfg.nx += 1; }
    if(cfg.isx_epitch){ cfg.IEPITCH = cfg.nx; cfg.nx += 1; }
    if(cfg.isx_eyaw){ cfg.IEYAW = cfg.nx; cfg.nx += 1; }
    if(cfg.isx_armgps){ cfg.IARMGPS = cfg.nx; cfg.nx += 1; }
    if(cfg.isx_armcar){ cfg.IARMCAR = cfg.nx; cfg.nx += 1; }

    /* configure for intergral vector */
    cfg.nitg = 0;
    if(cfg.is_odincre){ cfg.ITGVEL_OD = cfg.nitg; cfg.nitg += 3; }
    if(cfg.is_odincre && cfg.isx_kod){ cfg.ITGC_KOD = cfg.nitg; cfg.nitg += 3; }
    if(cfg.is_odincre){
        cfg.ITGVEL_INS = cfg.nitg; cfg.nitg += 3;
        cfg.ITGF_INS = cfg.nitg; cfg.nitg += 3;
    }
    if(cfg.is_odincre && cfg.isx_eyaw){
        cfg.ITGC_EYAW = cfg.nitg; cfg.nitg += 3;
    }
    if(cfg.is_odincre && cfg.isx_eroll){
        cfg.ITGC_EROLL = cfg.nitg; cfg.nitg += 3;
    }
    if(cfg.is_odincre && cfg.isx_epitch){
        cfg.ITGC_EPITCH = cfg.nitg; cfg.nitg += 3;
    }

    /*  configure for sagehusa auto adative algorithm */
    cfg.nR = 0;
    if(cfg.iskf_itgdS_sagehusa){
        cfg.IR_itgdS = cfg.nR;  cfg.nR += 3;
    }
    return 0;
}

/**
 * @brief generate stanard normal distribution random number
 * @return radnom number
 */
extern double yins_randn(void)
{
    double N = 25.0, sum_x = 0.0;
    for(int i = 0; i < N; ++i){
        sum_x += ((double)(rand())) / RAND_MAX;
    }
    /* standard normal distrubution */
    return  (sum_x - 0.5*N) / sqrt(N/12.0);
}

/**
 * @brief generate an 3D vector filled with normal distribution random
 *      number
 * @param[in] mean      mean of normal distribution
 * @param[in] sigma     standard error of normal distribution
 * @return 3D random vector number
 * @see yins_randn()
 */
extern v3_t v3_randn(double mean, double sigma)
{
    return (v3_t){mean + sigma*yins_randn(), mean + sigma*yins_randn(),
                mean + sigma*yins_randn()};
}

/** FB UD LR */
void imu_orientation_adjust(imud_t *imud, const char *or1, const char *or2)
{
    imud_t old_imud = *imud;
    const char OR[3][2] = {"FB", "UD", "LR"};
    double *a = (double *)&imud->accel;
    double *old_a = (double *)&old_imud.accel;
    double *g = (double *)&imud->gyro;
    double *old_g = (double *)&old_imud.gyro;

    for(int i = 0; i < 3; i++){
        if(or1[i] == OR[i][0] || or1[i] == OR[i][1]){
            for(int j = 0; j < 3; j++){
                if(or2[j] ==  OR[i][0] || or2[j] == OR[i][1]){
                    if(or1[i] == or2[j] && i != j){
                        a[i] = old_a[j];
                        g[i] = old_g[j];
                    }
                    if(or1[i] != or2[j]){
                        a[i] = - old_a[j];
                        g[i] = - old_g[j];
                    }
                    break;
                }
            }
        }
    }
}

/**
 * @brief Add solution type to status
 * @param[in] status        old status
 * @param[in] SOL_TYPE      solution type, see enum SOL
 * @return  added SOL_TYPE new status
 */
inline extern unsigned int soltype_add(unsigned int status,
                                       unsigned int SOL_TYPE)
{
    return status|SOL_TYPE;
}

/**
 * @brief remove solution type to status
 * @param[in] status        old status
 * @param[in] SOL_TYPE      solution type, see enum SOL
 * @return  removed SOL_TYPE new status
 */
inline extern unsigned int soltype_remove(unsigned int status,
                                          unsigned int SOL_TYPE)
{
    return status&(~(SOL_TYPE));
}
/**
 * @brief judge if status contains SOL_TYPE or not
 * @param[in] status        current solution status
 * @param[in] SOL_TYPE      solution type, see enum SOL
 * @return  removed SOL_TYPE new status
 */
inline extern bool is_soltype(unsigned int status, unsigned int SOL_TYPE)
{
    return (SOL_TYPE == (status & SOL_TYPE));
}

inline extern bool is_blh(const v3_t *pos){
    if(fabs(pos->x) < PI+1e-6 && fabs(pos->y) < 2*PI+1e-6){
        if(fabs(pos->z) < 1e6) return true;
        else if(fabs(v3_norm(*pos)-wgs84.R0) < 1e6) return false;
        else{
            LOG_WARN("is_blh: pos type may be incorrect(%f, %f, %f)",
                     pos->x, pos->y, pos->z);
            return false;
        }
    }
    return false;
}

/**
 * @brief Convert Markov random process standard error(often named unstability)
 *      to it's driven white noise(random walk coefficient)
 * @param[in] std   standard error      [SI]
 * @param[in] T     correlative time    [SI]
 * @return random walk coefficient [SI]
 * @see markov_rw2std()
 * @note  ref Yan, IMU test & data analysis, 2012, P152
 */
inline extern double
markov_std2rw(double std, double T) { return sqrt(2.0*SQR(std))/T; }

/**
 * @brief Convert Markov random process driven white noise(random walk
 *  coefficient to standard error
 * @param[in] rw    random warlk coeffcient [SI]
 * @param[in] T     correlative time   [SI]
 * @return standard error
 * @see markov_std2rw()
 * @note  ref Yan, IMU test & data analysis, 2012, P152
 */
inline extern double
markov_rw2std(double rw, double T) { return sqrt(SQR(rw*T)/2.0); }
