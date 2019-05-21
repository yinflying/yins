#include "ins.h"
#include <criterion/criterion.h>
#include <math.h>

#define EPS 1E-6
#define PI 3.1415926

Test(ned2ecef, real)
{
    v3_t pos = { 1.0, 1.0, 5.0 };
    v3_t vel = { 1.0, 2.0, 3.0 };
    m3_t att = { 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 };
    ned2ecef(&pos, &vel, &att);
    cr_expect_float_eq(pos.i, 1866377.86324143, EPS);
    cr_expect_float_eq(pos.j, 2906711.30133712, EPS);
    cr_expect_float_eq(pos.k, 5343772.65324308, EPS);
    cr_expect_float_eq(vel.i, -3.01337042820792, EPS);
    cr_expect_float_eq(vel.j, -0.991414946775814, EPS);
    cr_expect_float_eq(vel.k, -1.98411064855555, EPS);

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

Test(ecef2ned, real)
{
    v3_t pos = { -1890789.0, 5194902.0, 3170398.0 };
    v3_t vel = { 10.0, 15.0, 20.0 };
    m3_t att = { 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 };
    ecef2ned(&pos, &vel, &att);
    cr_expect_float_eq(pos.i, 0.523598749380288, EPS);
    cr_expect_float_eq(pos.j, 1.91986205856059, EPS);
    cr_expect_float_eq(pos.k, 48.8176567496266, EPS);
    cr_expect_float_eq(vel.i, 11.9829137792825, EPS);
    cr_expect_float_eq(vel.j, -14.5272270913126, EPS);
    cr_expect_float_eq(vel.k, -19.2449850713300, EPS);

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
                eular.i = (double)i1 / N * 2 * PI - PI;
                eular.j = (double)i2 / N * PI - PI / 2;
                eular.k = (double)i3 / N * 2 * PI;
                euler2dcm(&eular, &att);
                old_att = att;
                for (int j1 = 0; j1 < N; ++j1) {
                    for (int j2 = 0; j2 < N; ++j2) {
                        for (int j3 = 0; j3 < N; ++j3) {
                            pos.i = (double)j1 / N * PI - PI / 2;
                            pos.j = (double)j2 / N * 2 * PI - PI;
                            pos.k = (double)j3 / N * 10000 - 5000;
                            ned2ecef(&pos, &vel, &att);
                            ecef2ned(&pos, &vel, &att);

                            cr_expect_float_eq(pos.i, pos.i, EPS);
                            cr_expect_float_eq(pos.j, pos.j, EPS);
                            cr_expect_float_eq(pos.k, pos.k, EPS);
                            cr_expect_float_eq(vel.i, vel.i, EPS);
                            cr_expect_float_eq(vel.j, vel.j, EPS);
                            cr_expect_float_eq(vel.k, vel.k, EPS);
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
    v3_t euler = {3.1415, 0.68, 3.1415};
    /* forward */
    m3_t m1;
    euler2dcm(&euler, &m1);
    quat_t q1;
    dcm2quat(&m1, &q1);
    v3_t e1;
    quat2euler(&q1, &e1);

    cr_expect_float_eq(e1.i,euler.i,EPS);
    cr_expect_float_eq(e1.j,euler.j,EPS);
    cr_expect_float_eq(e1.k,euler.k,EPS);

    /* backward */
    quat_t q2;
    euler2quat(&euler, &q2);
    m3_t m2;
    quat2dcm(&q2, &m2);
    v3_t e2;
    dcm2euler(&m2, &e2);

    cr_expect_float_eq(e2.i,euler.i,EPS);
    cr_expect_float_eq(e2.j,euler.j,EPS);
    cr_expect_float_eq(e2.k,euler.k,EPS);
}
