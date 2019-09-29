/**
 * @file unittest.c
 * @brief ins.h functions unit tests
 * @author yinflying(yinflying@foxmail.com)
 * @version 0.0.1
 * @note
 *  2019-07-01  Add header comments
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
#include <criterion/criterion.h>
#include <math.h>
#include <stdio.h>

#define EPS 1E-6
#define PI 3.1415926
#define DEG2RAD 0.0174532925199433
#define RAD2DEG 57.2957795130823

static void m3_print(const m3_t* A)
{
    printf("-| %10.6f %10.6f %10.6f |-\n", A->m11, A->m12, A->m13);
    printf(" | %10.6f %10.6f %10.6f |\n", A->m21, A->m22, A->m23);
    printf("-| %10.6f %10.6f %10.6f |-\n", A->m31, A->m32, A->m33);
}
static void v3_print(const v3_t* V)
{
    printf("[ %10.6f %10.6f %10.6f ]\n", V->x, V->y, V->z);
}

Test(ned2ecef, real)
{
    v3_t pos = { 1.0, 1.0, 5.0 };
    v3_t vel = { 1.0, 2.0, 3.0 };
    m3_t att = { 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 };
    ned2ecef(&pos, &vel, &att);
    cr_expect_float_eq(pos.x, 1866377.86324143, EPS);
    cr_expect_float_eq(pos.y, 2906711.30133712, EPS);
    cr_expect_float_eq(pos.z, 5343772.65324308, EPS);
    cr_expect_float_eq(vel.x, -3.01337042820792, EPS);
    cr_expect_float_eq(vel.y, -0.991414946775814, EPS);
    cr_expect_float_eq(vel.z, -1.98411064855555, EPS);

    cr_expect_float_eq(att.m11, -0.454648713412841, EPS);
    cr_expect_float_eq(att.m12, -0.841470984807897, EPS);
    cr_expect_float_eq(att.m13, -0.291926581726429, EPS);
    cr_expect_float_eq(att.m21, -0.708073418273571, EPS);
    cr_expect_float_eq(att.m22, 0.540302305868140, EPS);
    cr_expect_float_eq(att.m23, -0.454648713412841, EPS);
    cr_expect_float_eq(att.m31, 0.540302305868140, EPS);
    cr_expect_float_eq(att.m32, 0.0, EPS);
    cr_expect_float_eq(att.m33, -0.841470984807897, EPS);
}

Test(v3,simple)
{
    v3_t v1 = {1.0, 2.0, 3.0};
    v3_t v2 = {-1.0, -2.0, -3.0};
    v3_t vadd = v3_add(v1,v2);
    v3_t vadd_result = V0;
    cr_expect(v3_equal(&vadd,&vadd_result,1e-14) ==  true);
    v3_t vdel = v3_del(v1,v2);
    v3_t vdel_result = {2.0, 4.0, 6.0};
    cr_expect(v3_equal(&vdel,&vdel_result,1e-14) == true);
    m3_t v1x = v3_askew(v1);
    v3_t v1_result = m3_iaskew(v1x);
    cr_expect(v3_equal(&v1,&v1_result,1e-14) == true);
}

Test(ecef2ned, real)
{
    v3_t pos = { -1890789.0, 5194902.0, 3170398.0 };
    v3_t vel = { 10.0, 15.0, 20.0 };
    m3_t att = { 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 };
    ecef2ned(&pos, &vel, &att);
    cr_expect_float_eq(pos.x, 0.523598749380288, EPS);
    cr_expect_float_eq(pos.y, 1.91986205856059, EPS);
    cr_expect_float_eq(pos.z, 48.8176567496266, EPS);
    cr_expect_float_eq(vel.x, 11.9829137792825, EPS);
    cr_expect_float_eq(vel.y, -14.5272270913126, EPS);
    cr_expect_float_eq(vel.z, -19.2449850713300, EPS);

    cr_expect_float_eq(att.m11, 0.171010008157753, EPS);
    cr_expect_float_eq(att.m12, -0.46984630934426, EPS);
    cr_expect_float_eq(att.m13, 0.866025416893444, EPS);
    cr_expect_float_eq(att.m21, -0.93969266136083, EPS);
    cr_expect_float_eq(att.m22, -0.34202003184695, EPS);
    cr_expect_float_eq(att.m23, 0.0, EPS);
    cr_expect_float_eq(att.m31, 0.296198040666165, EPS);
    cr_expect_float_eq(att.m32, -0.81379772880672, EPS);
    cr_expect_float_eq(att.m33, -0.49999997729453, EPS);
}

Test(ecef2ned_ned2ecef, real)
{
    v3_t pos, old_pos;
    v3_t vel = { 1.0, 2.0, 3.0 }, old_vel = vel;
    v3_t eular;
    m3_t att, old_att;
    double N = 5.0;
    for (int i1 = 0; i1 < N; ++i1) {
        for (int i2 = 0; i2 < N; ++i2) {
            for (int i3 = 0; i3 < N; ++i3) {
                eular.x = (double)i1 / N * 2 * PI - PI;
                eular.y = (double)i2 / N * PI - PI / 2;
                eular.z = (double)i3 / N * 2 * PI;
                euler2dcm(&eular, &att);
                old_att = att;
                for (int j1 = 0; j1 < N; ++j1) {
                    for (int j2 = 0; j2 < N; ++j2) {
                        for (int j3 = 0; j3 < N; ++j3) {
                            pos.x = (double)j1 / N * PI - PI / 2;
                            pos.y = (double)j2 / N * 2 * PI - PI;
                            pos.z = (double)j3 / N * 10000 - 5000;
                            ned2ecef(&pos, &vel, &att);
                            ecef2ned(&pos, &vel, &att);

                            cr_expect_float_eq(pos.x, pos.x, EPS);
                            cr_expect_float_eq(pos.y, pos.y, EPS);
                            cr_expect_float_eq(pos.z, pos.z, EPS);
                            cr_expect_float_eq(vel.x, vel.x, EPS);
                            cr_expect_float_eq(vel.y, vel.y, EPS);
                            cr_expect_float_eq(vel.z, vel.z, EPS);
                            cr_expect_float_eq(att.m11, att.m11, EPS);
                            cr_expect_float_eq(att.m12, att.m12, EPS);
                            cr_expect_float_eq(att.m13, att.m13, EPS);
                            cr_expect_float_eq(att.m21, att.m21, EPS);
                            cr_expect_float_eq(att.m22, att.m22, EPS);
                            cr_expect_float_eq(att.m23, att.m23, EPS);
                            cr_expect_float_eq(att.m31, att.m31, EPS);
                            cr_expect_float_eq(att.m32, att.m32, EPS);
                            cr_expect_float_eq(att.m33, att.m33, EPS);
                        }
                    }
                }
            }
        }
    }
}

Test(att_trans, real)
{
    v3_t euler = { 3.1415, 0.68, 3.1415 };
    /* forward */
    m3_t m1;
    euler2dcm(&euler, &m1);
    quat_t q1;
    dcm2quat(&m1, &q1);
    v3_t e1;
    quat2euler(&q1, &e1);

    cr_expect_float_eq(e1.x, euler.x, EPS);
    cr_expect_float_eq(e1.y, euler.y, EPS);
    cr_expect_float_eq(e1.z, euler.z, EPS);

    /* backward */
    quat_t q2;
    euler2quat(&euler, &q2);
    m3_t m2;
    quat2dcm(&q2, &m2);
    v3_t e2;
    dcm2euler(&m2, &e2);

    cr_expect_float_eq(e2.x, euler.x, EPS);
    cr_expect_float_eq(e2.y, euler.y, EPS);
    cr_expect_float_eq(e2.z, euler.z, EPS);

    v3_t rv1_q, rv1_d;
    quat2rv(&q1, &rv1_q);
    dcm2rv(&m1, &rv1_d);
    cr_expect(v3_equal(&rv1_q, &rv1_d, 1E-14) == true);
    quat_t q1_result;
    m3_t m1_result;
    rv2quat(&rv1_q, &q1_result);
    rv2dcm(&rv1_q, &m1_result);
    cr_expect(quat_equal(&q1, &q1_result, 1E-14) == true);
    cr_expect(m3_equal(&m1, &m1_result, 1E-14) == true);
}

Test(dblvec2att, simple)
{
    v3_t vn1 = { 0, 0, -9.79205619566115 };
    v3_t vn2 = { 0, 6.31515696436349e-05, 3.64605757335000e-05 };
    v3_t vb1 = { -0.0004890350025122, -0.0004890431881453, -9.79205617123735 };
    v3_t vb2
        = { 4.90421825303274e-08, 6.31426704787843e-05, 3.64759522280797e-05 };
    m3_t Cnb;
    dblvec2att(&vn1, &vn2, &vb1, &vb2, &Cnb);
    cr_expect_float_eq(Cnb.m11, 0.999999719107770, EPS);
    cr_expect_float_eq(Cnb.m12, 0.000747857056062642, EPS);
    cr_expect_float_eq(Cnb.m13, 4.99420134791420e-05, EPS);
    cr_expect_float_eq(Cnb.m21, -0.000747859550308376, EPS);
    cr_expect_float_eq(Cnb.m22, 0.999999719105863, EPS);
    cr_expect_float_eq(Cnb.m23, 4.99428494254436e-05, EPS);
    cr_expect_float_eq(Cnb.m31, -4.99046493383806e-05, EPS);
    cr_expect_float_eq(Cnb.m32, -4.99801850086273e-05, EPS);
    cr_expect_float_eq(Cnb.m33, 0.999999997505754, EPS);

    v3_t vb1t = m3_mul_v3(Cnb, vn1);
    v3_t vb2t = m3_mul_v3(Cnb, vn2);
    cr_expect_float_eq(vb1t.x, vb1.x, EPS);
    cr_expect_float_eq(vb1t.y, vb1.y, EPS);
    cr_expect_float_eq(vb1t.z, vb1.z, EPS);
    cr_expect_float_eq(vb2t.x, vb2.x, EPS);
    cr_expect_float_eq(vb2t.y, vb2.y, EPS);
    cr_expect_float_eq(vb2t.z, vb2.z, EPS);
}

Test(align_coarse_static_base, simple)
{
    imu_t imu;
    imu.property = NULL;
    imud_t imud;
    imu.n = 0;
    imu.nmax = 0;
    double ep[6] = { 2019, 1, 1, 0, 0, 0.0 };
    imud.time = epoch2time(ep);
    imud.accel = (v3_t) { 4.8898352e-04, 4.89016111e-04, -9.792545188 };
    imud.gyro = (v3_t) { 6.31998667e-05, 4.863401409e-08, -3.650871198e-05 };
    for (int i = 0; i < 10; i++) {
        imu_add(&imu, &imud);
        ep[5] += 1.0;
        imud.time = epoch2time(ep);
    }
    double lat = 0.523598775598299;
    m3_t Cnb;
    v3_t Enb;

    align_coarse_static_base(&imu, lat, &Cnb);
    dcm2euler(&Cnb, &Enb);
    printf("Align_coarse_static_base:simple %f %f %f\n", Enb.x, Enb.y, Enb.z);

    imu_free(&imu);
    cr_expect_float_eq(Enb.x, -4.994210169e-05, EPS);
    cr_expect_float_eq(Enb.y, 4.9942061337e-05, EPS);
    cr_expect_float_eq(Enb.z, 6.2824519208, 1E-3);
}

Test(align_coarse_inertial, simple)
{
    imu_t imus; imu_init(&imus);
    yins_readf("data/align_test_data.ycsv", FT_IMU_YCSV, &imus, NULL, NULL);
    double lat = 0.523598775598299;
    m3_t Cnb;
    v3_t Enb;
    int maxn = imus.n;
    align_coarse_inertial(&imus, lat, &Cnb);

    dcm2euler(&Cnb, &Enb);
    printf("Align_coarse_inertial:simple %f %f %f\n", Enb.x, Enb.y, Enb.z);

    imu_free(&imus);
    cr_expect_float_eq(Enb.x, -1.04880393840960e-05, 1E-5);
    cr_expect_float_eq(Enb.y, 5.04490437057051e-05, 1E-5);
    cr_expect_float_eq(Enb.z, 6.28242649768121, 1E-3);
}

Test(align_coarse_wuhba, simple)
{
    imu_t imus; imu_init(&imus);
    yins_readf("data/align_test_data.ycsv", FT_IMU_YCSV, &imus, NULL, NULL);
    double lat = 0.523598775598299;
    m3_t Cnb;
    v3_t Enb;
    int maxn = imus.n;

    const int N = 5;
    v3_t eb_n[N];
    for(int i = 0; i < N; ++i){
        eb_n[i] = (v3_t){0.0, 0.0, 0.0};
    }

    align_coarse_wuhba(&imus, lat, eb_n, N, &Cnb);

    dcm2euler(&Cnb, &Enb);
    printf("Align_coarse_wuhba:simple %f %f %f\n", Enb.x, Enb.y, Enb.z);

    imu_free(&imus);

    cr_expect_float_eq(Enb.x, 0.0, 1E-2);
    cr_expect_float_eq(Enb.y, 0.0, 1E-2);
    cr_expect_float_eq(Enb.z, 6.28242649768121, 1E-2);
}

Test(m3_SVD, simple)
{
    m3_t A = { 1.0, 2.0, 3.0, 2.0, 3.0, 6.0, 2.0, 1.0, 3.0 };
    m3_t U, V;
    v3_t D;
    m3_SVD(&A, &U, &D, &V);
    m3_t B = m3_mul(m3_mul(U, v3_diag(D)), m3_T(V));
    cr_expect_eq(m3_equal(&A, &B, 1E-14), true);
}

Test(m3_swap_row, simple)
{
    m3_t A = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0 };
    m3_t orignA = A;
    m3_swap_row(&A, 1, 2);
    m3_swap_row(&A, 2, 1);
    cr_expect(m3_equal(&A, &orignA, 1e-20) == true);

    m3_swap_row(&A, 1, 3);
    m3_swap_row(&A, 3, 1);
    cr_expect(m3_equal(&A, &orignA, 1e-20) == true);

    m3_swap_row(&A, 2, 3);
    m3_swap_row(&A, 3, 2);
    cr_expect(m3_equal(&A, &orignA, 1e-20) == true);

    m3_swap_row(&A, 1, 2);
    m3_swap_row(&A, 2, 3);
    m3_swap_row(&A, 1, 2);
    m3_t B = { 7.0, 8.0, 9.0, 4.0, 5.0, 6.0, 1.0, 2.0, 3.0 };
    cr_expect(m3_equal(&A, &B, 1e-20) == true);
}

Test(m3_swap_clm, simple)
{
    m3_t A = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0 };
    m3_t orignA = A;
    m3_swap_clm(&A, 1, 2);
    m3_swap_clm(&A, 2, 1);
    cr_expect(m3_equal(&A, &orignA, 1e-20) == true);

    m3_swap_clm(&A, 1, 3);
    m3_swap_clm(&A, 3, 1);
    cr_expect(m3_equal(&A, &orignA, 1e-20) == true);

    m3_swap_clm(&A, 2, 3);
    m3_swap_clm(&A, 3, 2);
    cr_expect(m3_equal(&A, &orignA, 1e-20) == true);

    m3_swap_clm(&A, 1, 2);
    m3_swap_clm(&A, 2, 3);
    m3_swap_clm(&A, 1, 2);
    m3_t B = { 3.0, 2.0, 1.0, 6.0, 5.0, 4.0, 9.0, 8.0, 7.0 };
    cr_expect(m3_equal(&A, &B, 1e-20) == true);
}

Test(m3_LU, simple)
{
    m3_t A, L, U, P, LU, PA;
    /* common matrix test */
    A = (m3_t) { 1.0, 2.0, 3.0, 22, 333.0, 1.0, 4.0, 5.0, 6.0 };
    m3_LU(&A, &L, &U, &P);
    cr_expect_float_eq(L.m12, 0.0, EPS);
    cr_expect_float_eq(L.m13, 0.0, EPS);
    cr_expect_float_eq(L.m23, 0.0, EPS);
    cr_expect_float_eq(U.m21, 0.0, EPS);
    cr_expect_float_eq(U.m31, 0.0, EPS);
    cr_expect_float_eq(U.m32, 0.0, EPS);
    LU = m3_mul(L, U);
    PA = m3_mul(P, A);
    cr_expect(m3_equal(&PA, &LU, EPS) == true);

    /* U matrix test */
    A = (m3_t) { 1.0, 2.0, 5.0, 0.0, 3.0, 1.0, 0.0, 0.0, 1.0 };
    m3_LU(&A, &L, &U, &P);
    cr_expect_float_eq(L.m12, 0.0, EPS);
    cr_expect_float_eq(L.m13, 0.0, EPS);
    cr_expect_float_eq(L.m23, 0.0, EPS);
    cr_expect_float_eq(U.m21, 0.0, EPS);
    cr_expect_float_eq(U.m31, 0.0, EPS);
    cr_expect_float_eq(U.m32, 0.0, EPS);
    LU = m3_mul(L, U);
    PA = m3_mul(P, A);
    cr_expect(m3_equal(&PA, &LU, EPS) == true);

    /* L matrix test */
    A = (m3_t) { 0.0, 0.0, 5.0, 0.0, 2.0, 7.0, 1.0, 3.0, 1.0 };
    m3_LU(&A, &L, &U, &P);
    cr_expect_float_eq(L.m12, 0.0, EPS);
    cr_expect_float_eq(L.m13, 0.0, EPS);
    cr_expect_float_eq(L.m23, 0.0, EPS);
    cr_expect_float_eq(U.m21, 0.0, EPS);
    cr_expect_float_eq(U.m31, 0.0, EPS);
    cr_expect_float_eq(U.m32, 0.0, EPS);
    LU = m3_mul(L, U);
    PA = m3_mul(P, A);
    cr_expect(m3_equal(&PA, &LU, EPS) == true);

    A = (m3_t){ 11.0, 2323.00, -3.343, 1e-3, 333.0, 1.0, 1.0, 3.0, 5.0 };
    m3_LU(&A, &L, &U, &P);
    cr_expect_float_eq(L.m12, 0.0, EPS);
    cr_expect_float_eq(L.m13, 0.0, EPS);
    cr_expect_float_eq(L.m23, 0.0, EPS);
    cr_expect_float_eq(U.m21, 0.0, EPS);
    cr_expect_float_eq(U.m31, 0.0, EPS);
    cr_expect_float_eq(U.m32, 0.0, EPS);
    LU = m3_mul(L, U);
    PA = m3_mul(P, A);
    cr_expect(m3_equal(&PA, &LU, EPS) == true);

    A = (m3_t) { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 };
    m3_LU(&A, &L, &U, &P);
    cr_expect_float_eq(L.m12, 0.0, EPS);
    cr_expect_float_eq(L.m13, 0.0, EPS);
    cr_expect_float_eq(L.m23, 0.0, EPS);
    cr_expect_float_eq(U.m21, 0.0, EPS);
    cr_expect_float_eq(U.m31, 0.0, EPS);
    cr_expect_float_eq(U.m32, 0.0, EPS);
    LU = m3_mul(L, U);
    PA = m3_mul(P, A);
    cr_expect(m3_equal(&PA, &LU, EPS) == true);

}

Test(m3_det, simple)
{
    m3_t A = { 1.0, 2.0, 3.0, 1.0, -3.0, 4.0, 1.0, 2.2, 77.0 };
    cr_expect_float_eq(m3_det(&A), -370.2, 1e-14);
}

Test(v3_normalize,simple)
{
    v3_t v = {3.0, 1e-10, -3};
    int ret = v3_normalize(&v);
    v3_t v_check = {0.707106781186548,2.35702260395516e-11,-0.707106781186548};
    cr_expect(v3_equal(&v,&v_check,1e-14) == true);
    cr_expect_eq(ret,0);

    v = (v3_t){1e-60, 1e-80, -1e-90};
    ret = v3_normalize(&v);
    v_check = (v3_t){1e-60, 1e-80, -1e-90};
    cr_expect(v3_equal(&v,&v_check,1e-14) == true);
    cr_expect_eq(ret,1);
}

Test(m3_inv,simple)
{
    m3_t A, A_check;
    A = (m3_t){1.2, 3, 2.5, -234, 1e-2, 0.0, 1.1, 2.2, 3.3};
    m3_inv(&A);
    A_check = (m3_t){3.20509053846590e-05, -0.00427345405128786, -2.4280988927771e-05,
        0.749991186001019, 0.00117519986410416, -0.568175140909863,
        -0.50000480763580, 0.00064101810769317, 0.681821823966521};
    cr_expect(m3_equal(&A, &A_check, EPS) == true);

    A = (m3_t){1.0, 2.0, 5.0, 0.0, 3.0, 1.0, 0.0, 0.0, 1.0};
    A_check = A;
    m3_inv(&A); m3_inv(&A);
    cr_expect(m3_equal(&A, &A_check, EPS) == true);

    A = (m3_t){1.0, 2.0, 5.0, 0.0, 3.0, 1.0, 5.0, 2.0, 1.0};
    A_check = A;
    m3_inv(&A); m3_inv(&A);
    cr_expect(m3_equal(&A, &A_check, EPS) == true);
}

Test(quat, simple)
{
    quat_t P,Q;

    P = (quat_t){1.0, 1.0, 4.0, 8.0};
    quat_normalize(&P);
    quat_t P_result = {0.1104, 0.1104, 0.4417, 0.8835};
    cr_expect(quat_equal(&P,&P_result,1e-3) == true);

    P = (quat_t){1.0, 1.0, 4.0, 8.0};
    quat_conj(&P);
    P_result = (quat_t){1.0, -1.0, -4.0,-8.0};
    cr_expect(quat_equal(&P,&P_result,1e-3) == true);

    P = (quat_t){1.0, 1.0, 4.0, 8.0};
    quat_t P_inv = P; quat_inv(&P_inv);
    P_result = quat_mul(P,P_inv);
    quat_t P_expect = {1.0, 0.0, 0.0, 0.0};
    cr_expect(quat_equal(&P_result,&P_expect,EPS));

    P = (quat_t){1.0, 1.0, 4.0, 8.0};
    Q = (quat_t){1.0, 3.0, 2.0, 9.0};
    cr_expect_float_eq(quat_dot(P,Q),84,EPS);
    cr_expect_float_eq(quat_norm(P),9.0554,1e-4);

    P = (quat_t){1.0, 1.0, 4.0, 8.0};
    Q = (quat_t){1.0, 3.0, 2.0, 9.0};
    v3_t v_result = quat_cross(P, Q);
    v3_t v_check = {20.0, 15.0, -10.0};
    cr_expect(v3_equal(&v_result,&v_check,EPS) == true);
}

Test(mat,simple)
{
    double *A = mat(3,3);
    A[0*3 + 0] = 1.2;
    A[1*3 + 0] = 3;
    A[2*3 + 0] = 2.5;

    A[0*3 + 1] = -234;
    A[1*3 + 1] = 1e-2;
    A[2*3 + 1] = 0.0;

    A[0*3 + 2] = 1.1;
    A[1*3 + 2] = 2.2;
    A[2*3 + 2] = 3.3;

    matinv(A, 3);

    m3_t A_check = (m3_t){3.20509053846590e-05, -0.00427345405128786, -2.4280988927771e-05,
        0.749991186001019, 0.00117519986410416, -0.568175140909863,
        -0.50000480763580, 0.00064101810769317, 0.681821823966521};
    m3_t *m3A = (m3_t *)A;
    m3_t B = m3_T(*m3A);
    cr_expect(m3_equal(&B, &A_check, EPS) == true);

    free(A);
}

Test(mat_symmetry, simple)
{
    double A[3*3] = {1.0, 3.0, 2.0, 3.0, 4.0, 4.0, 4.0, 5.0, 7.0};
    mat_symmetry(A, 3);
    cr_expect_float_eq(A[0],1.0,1e-32);
    cr_expect_float_eq(A[1],3.0,1e-32);
    cr_expect_float_eq(A[2],3.0,1e-32);
    cr_expect_float_eq(A[3],3.0,1e-32);
    cr_expect_float_eq(A[4],4.0,1e-32);
    cr_expect_float_eq(A[5],4.5,1e-32);
    cr_expect_float_eq(A[6],3.0,1e-32);
    cr_expect_float_eq(A[7],4.5,1e-32);
    cr_expect_float_eq(A[8],7.0,1e-32);
}

Test(v3_list_operation, simple)
{
    const int n = 3;
    v3_t v3_list[n];
    v3_list[0] = (v3_t){1.0, 2.0, 3.0};
    v3_list[1] = (v3_t){2.0, 3.0, 1.0};
    v3_list[2] = (v3_t){-1.0,1.0, 0.0};

    v3_t vsum = v3_sum(v3_list, n);
    cr_expect_float_eq(vsum.x,2.0,1e-32);
    cr_expect_float_eq(vsum.y,6.0,1e-32);
    cr_expect_float_eq(vsum.z,4.0,1e-32);

    v3_t vmean = v3_mean(v3_list,n);
    cr_expect_float_eq(vmean.x,2.0/3.0,1e-32);
    cr_expect_float_eq(vmean.y,2.0,1e-32);
    cr_expect_float_eq(vmean.z,1.0+1.0/3.0,1e-32);

    v3_t vrms = v3_rms(v3_list,n);
    cr_expect_float_eq(vrms.x,1.41421356237310,1e-14);
    cr_expect_float_eq(vrms.y,2.16024689946929,1e-14);
    cr_expect_float_eq(vrms.z,1.82574185835055,1e-14);

    v3_t vstd = v3_std(v3_list,n);
    printf("%f\n", vstd.x);
    cr_expect_float_eq(vstd.x,1.52752523165195,1e-14);
    cr_expect_float_eq(vstd.y,1.0,1e-14);
    cr_expect_float_eq(vstd.z,1.52752523165195,1e-14);
}
