/**
 * @file ins.c
 * @brief ins core functions implementation
 * @author yinflying(yinflying@foxmail.com)
 * @version 0.0.1
 * @note
 *  2019-05-21  Created
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
#include <string.h>

#define PI 3.14159265358979
#define SQR(x) ((x) * (x))

/**
 * @brief Gravitational acceleration of Earth project to e-frame
 * @param[in]   r   position under e-axis
 * @param[out]  ge  Gravitational acceleration under ECEF(m s^-2)
 * @return 0: OK
 * @see gravity_ned()
 * @note ge do NOT contain centrifugal force
 */
extern int gravity_ecef(const v3_t* r, v3_t* ge)
{
    double mag_r = sqrt(r->x * r->x + r->y * r->y + r->z * r->z);
    if (fabs(mag_r) < 1e-32) {
        ge->x = 0.0; ge->y = 0.0; ge->z = 0.0;
    } else {
        /* Calculate gravitational acceleration using (2.142) */
        double f1 = -wgs84.mu / pow(mag_r, 3);
        double f2 = 1.5 * wgs84.J2 * pow(wgs84.R0 / mag_r, 2.0);
        double z_scale = 5.0 * pow((r->z / mag_r), 2.0);
        double g1 = f1 * (r->x + f2 * (1.0 - z_scale) * r->x);
        double g2 = f1 * (r->y + f2 * (1.0 - z_scale) * r->y);
        double g3 = f1 * (r->z + f2 * (3.0 - z_scale) * r->z);
        /* Add centripetal acceleration using (2.133) */
        ge->x = g1 + wgs84.wie * wgs84.wie * r->x;
        ge->y = g2 + wgs84.wie * wgs84.wie * r->y;
        ge->z = g3;
    }
    return 0;
}

/**
 * @brief Acceleration of gravity under n-axis(NED)
 * @param[in]   lat     latitude [rad]
 * @param[in]   hgt     ellipsoidal height [m]
 * @param[out]  gn      acceleration of gravity [m s^-2]
 * @return 0: OK
 * @see gravity_ecef()
 * @note gravity contains two part: gravitational and centrifugal accelration
 */
extern int gravity_ned(double lat, double hgt, v3_t* gn)
{
    double sinlat2 = sin(lat) * sin(lat);
    double e2 = SQR(wgs84.e);

    /* Calculate surface gravity using Somigliana model */
    double g0 = 9.7803253359 * (1.0 + 0.001931853 * sinlat2)
        / sqrt(1.0 - e2 * sinlat2);

    gn->x = -8.08E-9 * hgt * sin(2.0 * lat); /* North */
    gn->y = 0.0;                             /* East */
    /* Down */
    double tmp = 1.0 + wgs84.f * (1.0 - 2.0 * sinlat2)
        + (SQR(wgs84.wie) * SQR(wgs84.R0) * wgs84.RP / wgs84.mu);
    gn->z = g0 *
        (1.0 - (2.0 / wgs84.R0) * tmp * hgt + (3.0 * SQR(hgt) / SQR(wgs84.R0)));

    return 0;
}

/**
 * @brief Strapdown-INS equations under ECEF frame
 * @param[in]   dt      Time interval [s]
 * @param[in]   dtheta  Angular increment [rad]
 * @param[in]   dv      Velocity increment [m/s]
 * @param[in,out]   r   Start/End postion in ECEF(reb_e) [m]
 * @param[in,out]   v   Start/End velocity in ECEF(veb_e) [m]
 * @param[in,out]   q   Start/End attitude trans express by quaternion(qbe)
 * @return 0: OK
 * @see multisample()
 * @note This function do not contain conning&sculling error compensation,
 *      so it should work with multisample() function.
 *
 *      Ref: Paul. D. Groves. Principles of GNSS, Inertial, and Multisensor
 *      Integrated Navigation Systems(2nd Edition), P163
 */
extern int ins_nav_ecef(
    double dt, const v3_t* dtheta, const v3_t* dv, v3_t* r, v3_t* v, quat_t* q)
{
    /* Attitude update */
    v3_t dtheta_ie = {0.0, 0.0, -wgs84.wie * dt};
    quat_t q_earth, q_body, old_q = *q;
    rv2quat(&dtheta_ie, &q_earth);
    rv2quat(dtheta, &q_body);
    *q = quat_mul(quat_mul(q_earth, old_q), q_body);
    quat_normalize(q);
    /* Specific force transform(velocity form) */
    v3_t dtheta_ie_half = { 0.0, 0.0, -wgs84.wie * dt / 2.0 };
    rv2quat(&dtheta_ie_half, &q_earth);
    v3_t dv_rot = v3_scalar(0.5, v3_cross(*dtheta, *dv));
    v3_t dv_e = quat_mul_v3(quat_mul(q_earth, old_q), v3_add(*dv, dv_rot));
    /* Velocity update */
    v3_t ge; gravity_ecef(r, &ge);
    v3_t old_v = *v;
    v->x = old_v.x + dv_e.x + dt * (ge.x + 2.0 * wgs84.wie * old_v.y);
    v->y = old_v.y + dv_e.y + dt * (ge.y - 2.0 * wgs84.wie * old_v.x);
    v->z = old_v.z + dv_e.z + dt * ge.z;
    /* Position update */
    v3_t old_r = *r;
    *r = v3_add(old_r, v3_scalar(0.5 * dt, v3_add(old_v, *v)));
    return 0;
}

int ins_nav_ecef_back(double dt,const v3_t *dtheta, const v3_t *dv,
    v3_t *r, v3_t *v, quat_t *q)
{
    /* TODO: test this function */
    /* Attitude update */
    v3_t dtheta_ie = { 0.0, 0.0, wgs84.wie * dt };
    quat_t q_earth, q_body, old_q = *q;
    rv2quat(&dtheta_ie, &q_earth);
    rv2quat(dtheta,&q_body);
    quat_conj(&q_body);
    *q = quat_mul(quat_mul(q_earth, old_q), q_body);
    /* Specific force transform(velocity form) */
    v3_t dtheta_ie_half = { 0.0, 0.0, wgs84.wie * dt / 2.0 };
    rv2quat(&dtheta_ie_half, &q_earth);
    v3_t dv_rot = v3_scalar(-0.5, v3_cross(*dtheta, *dv));
    v3_t dv_e = quat_mul_v3(quat_mul(q_earth, old_q), v3_add(*dv, dv_rot));
    /* Velocity update */
    v3_t ge; gravity_ecef(r, &ge);
    v3_t old_v = *v;
    v->x = old_v.x + dv_e.x + dt * (ge.x + 2 * wgs84.wie * old_v.y);
    v->y = old_v.y + dv_e.y + dt * (ge.y - 2 * wgs84.wie * old_v.x);
    v->z = old_v.z + dv_e.z + dt * ge.z;

    return 0;
}

/**
 * @brief Use multi-subsample to compensate the conning&scull error
 * @param[in]   dtheta_list     Angular increment list,order by time,length:abs(N)
 * @param[in]   dv_list         Velocity increment list,order by time,length:abs(N)
 * @param[in]   N               Subsample number(1=<N<=5, -2: one-plus-previous)
 * @param[out]  dtheta          Sum of angular increment with conning error compensation
 * @param[out]  dv              Sum of velocity incrment with scull error compensation
 * @return 0: OK    1: Error
 * @note Ref:
 *  1. Paul. D. Groves. Principles of GNSS, Inertial, and Multisensor
 *      Integrated Navigation Systems(2nd Edition)
 *  2. Yan Gongming, 捷联惯导算法与组合导航讲义, 2016
 */
extern int multisample( const v3_t* dtheta_list, const v3_t* dv_list, int N,
        v3_t* dtheta, v3_t* dv)
{
    if (abs(N) == 1) {
        *dtheta = dtheta_list[0];
        *dv = dv_list[0];
        return 0;
    }
    /* Multi-subsample coning compensation */
    if (N >= 2 && N <= 5) {
        // coning error confficients, ref Yan2016(P38:table 2.5-2)
        static double conefactors[5][4] = {
            { 2. / 3 },                                         // 2
            { 9. / 20, 27. / 20 },                              // 3
            { 54. / 105, 92. / 105, 214. / 105 },               // 4
            { 250. / 504, 525. / 504, 650. / 504, 1375. / 504 } // 5
        };
        double* pcf = conefactors[N - 2];
        v3_t sum_c = {0.0, 0.0, 0.0}, sum_s = {0.0, 0.0, 0.0};
        v3_t sum_w = {0.0, 0.0, 0.0}, sum_v = {0.0, 0.0, 0.0};
        int i = 0;
        for (; i < N - 1; i++) {
            /* sum of cross_product factor */
            sum_c = v3_add(sum_c, v3_scalar(pcf[i], dtheta_list[i]));
            sum_s = v3_add(sum_s, v3_scalar(pcf[i], dv_list[i]));
            /* sum of angular increment and velocity increment */
            sum_w = v3_add(sum_w, dtheta_list[i]);
            sum_v = v3_add(sum_v, dv_list[i]);
        }
        sum_w = v3_add(sum_w, dtheta_list[i]);
        sum_v = v3_add(sum_v, dv_list[i]);
        /* coning error compensation for angular increment, ref Paul2013(5.97) */
        *dtheta = v3_add(sum_w, v3_cross(sum_c, dtheta_list[i]));
        /* sculling error compensation for velocity increment, ref Paul2013(5.98)*/
        v3_t scul = v3_add(
            v3_cross(sum_c, dv_list[i]), v3_cross(sum_s, dtheta_list[i]));
        *dv = v3_add(sum_v, scul);
        return 0;
    }
    if (N == -2) {
        v3_t sum_c = v3_scalar(1.0 / 12, dtheta_list[0]);
        v3_t sum_s = v3_scalar(1.0 / 12, dv_list[0]);
        /* ref Yan2016(P31:2.5-37) */
        *dtheta = v3_add(dtheta_list[1], v3_cross(sum_c, dtheta_list[1]));
        /* ref Yan2016(P73:4.1-36,P76:4.1-55) */
        v3_t scul = v3_add(v3_cross(sum_c, dv_list[1]),
                v3_cross(sum_s, dtheta_list[1]));
        *dv = v3_add(dv_list[1], scul);
        return 0;
    }
    return 1;
}

int dr_nav_ned()
{
    /* TODO: DR under ned */
    return 0;
}

/**
 * @brief Determinate attitude by using double vector
 * @param[in] vn1   First vector under n-frame
 * @param[in] vn2   Second vector under n-frame
 * @param[in] vb1   First vector under b-frame
 * @param[in] vb2   Second vector under b-frame
 * @param[out] Cnb  Attitude transform from n-frame to b-frame(or Attitude
 *      b-frame to n-frame)
 * @return 0: OK
 * @note Ref: Yan Gongming, 捷联惯导算法与组合导航讲义, 2016, P143:6.1-7
 */
extern int dblvec2att(const v3_t* vn1, const v3_t* vn2, const v3_t* vb1,
    const v3_t* vb2, m3_t* Cnb)
{
    v3_t vx = v3_scalar(1.0 / v3_norm(*vb1), *vb1);
    v3_t vtmp = v3_cross(*vb1, *vb2);
    v3_t vy = v3_scalar(1.0 / v3_norm(vtmp), vtmp);
    vtmp = v3_cross(vtmp, *vb1);
    v3_t vz = v3_scalar(1.0 / v3_norm(vtmp), vtmp);
    m3_t m3_b = { vx.x, vy.x, vz.x, vx.y, vy.y, vz.y, vx.z, vy.z, vz.z };

    vx = v3_scalar(1.0 / v3_norm(*vn1), *vn1);
    vtmp = v3_cross(*vn1, *vn2);
    vy = v3_scalar(1.0 / v3_norm(vtmp), vtmp);
    vtmp = v3_cross(vtmp, *vn1);
    vz = v3_scalar(1.0 / v3_norm(vtmp), vtmp);
    m3_t m3_n = { vx.x, vx.y, vx.z, vy.x, vy.y, vy.z, vz.x, vz.y, vz.z };

    *Cnb = m3_mul(m3_b, m3_n);
    return 0;
}

/**
 * @brief Coarse align on the static base
 * @param[in]   imu     static imu data
 * @param[in]   lat     imu latitude [rad]
 * @param[out]  Cnb     Output DCM attitude(a-axis refer to n-axis) on average
 * @return 0: Ok
 * @see align_coarse_inertial()
 * @see dblvec2att()
 * @note Yan Gongming, 捷联惯导算法与组合导航讲义, 2016, P144
 */
extern int align_coarse_static_base(const imu_t* imu, double lat, m3_t *Cnb)
{
    /* Fetch mean mean angular rate and specific force */
    v3_t mean_wib_b = { 0 }, mean_fib_b = { 0 };
    /* Start from the second(first epoch should not add in)*/
    for (unsigned int i = 1; i < imu->n; ++i) {
        mean_wib_b = v3_add(mean_wib_b, imu->data[i].gyro);
        mean_fib_b = v3_add(mean_fib_b, imu->data[i].accel);
    }
    double T = timediff(imu->data[imu->n - 1].time, imu->data[0].time);
    mean_wib_b = v3_scalar(1.0 / T, mean_wib_b);
    mean_fib_b = v3_scalar(1.0 / T, mean_fib_b);

    /* Double Vector to Attitude */
    v3_t gn; gravity_ned(lat,0.0,&gn);
    v3_t wie_n = {wgs84.wie * cos(lat),0.0,- wgs84.wie * sin(lat)};
    v3_t gib_b = v3_scalar(-1.0,mean_fib_b);

    dblvec2att(&gn,&wie_n,&gib_b,&mean_wib_b,Cnb);
    return 0;
}

/**
 * @brief Coarse align under inertial fame(anti-vibration method)
 * @param[in]   imu     (quasi-)static imu data(recommend imu->n/8 = 0)
 * @param[in]   lat     imu latitude [rad]
 * @param[out]  Cnb     Output DCM attitude(b-axis refer to n-axis) at last moment
 * @return 0: OK
 * @see align_coarse_static_base()
 * @see dblvec2att()
 * @note ref Qin Yunyuan, 惯性导航(2nd Edition), P319
 */
extern int align_coarse_inertial(const imu_t *imu, double lat, m3_t *Cnb)
{
    /* I-frame: inertial frame of n-frame at start moment */
    /* B-frame: inertial frame of b-frame at start moment */

    unsigned int sample_N = 4;
    double ts = timediff(imu->data[1].time,imu->data[0].time);
    /* double nts = sample_N * ts; */

    double sin_lat = sin(lat), cos_lat = cos(lat);
    v3_t gn; gravity_ned(lat,0.0,&gn);

    /* Calculate vib_B1,vib_B2*/
    v3_t dtheta[4],dv[4],sum_dtheta,sum_dv;
    v3_t vib_B = {0.0, 0.0, 0.0}, vib_B1 = {0.0, 0.0, 0.0};
    quat_t qb_B = {1.0, 0.0, 0.0, 0.0};   /* initial attitde*/
    quat_t qk_k1; /* trans from k to k+1 */
    unsigned int ind_mid = (imu->n/sample_N) / 2 * sample_N;
    for (unsigned int i = 0; i <= imu->n - sample_N ; i += sample_N) {
        for(unsigned int j = 0; j < sample_N; ++j){
            dtheta[j] = imu->data[i+j].gyro;
            dv[j] = imu->data[i+j].accel;
        }
        /* Calculate current fib_B */
        multisample(dtheta,dv,(int)sample_N,&sum_dtheta,&sum_dv);
        vib_B = v3_add(vib_B,quat_mul_v3(qb_B, sum_dv));
        /* qb_B attitude update uner inertial frame */
        rv2quat(&sum_dtheta,&qk_k1);
        qb_B = quat_mul(qb_B,qk_k1);

        /* record middle vib_B */
        if(i == ind_mid - sample_N) vib_B1 = vib_B;
    }
    /* Calculate vib_I1, vib_I2 */
    double total_t = ts * ind_mid * 2;
    double wie_dtheta = wgs84.wie * total_t;
    double gcl_wie = gn.z * cos_lat / wgs84.wie;
    v3_t vib_I1 = { gcl_wie * sin(wie_dtheta/2.0),
        gcl_wie * (1-cos(wie_dtheta/2.0)), total_t/2.0 * gn.z * sin_lat};
    v3_t vib_I2 = { gcl_wie * sin(wie_dtheta),
        gcl_wie * (1-cos(wie_dtheta)), total_t * gn.z * sin_lat};

    /* double vector to attitude */
    m3_t CB_I; dblvec2att(&vib_B1,&vib_B,&vib_I1,&vib_I2,&CB_I);

    /* Calculate Cnb */
    double cos_wie = cos(wie_dtheta), sin_wie = sin(wie_dtheta);
    m3_t CI_n = { -sin_lat * cos_wie, -sin_lat * sin_wie, cos_lat,
        -sin_wie, cos_wie, 0.0,
        -cos_lat * cos_wie, -cos_lat * sin_wie, -sin_lat };
    m3_t Cb_B; quat2dcm(&qb_B,&Cb_B);
    *Cnb = m3_T(m3_mul(m3_mul(CI_n,CB_I),Cb_B));

    return 0;
}

/**
 * @brief Coarse Alignment by solving Wuhba problem under inertial frame
 * @param[in]   imu     IMU data
 * @param[in]   lat     Imu latitude(Average Latitude) [rad]
 * @param[in]   veb_n   Imu velocity uner n-frame, length: Nveb_n [m/s]
 * @param[in]   Nveb_n  Imu velocity number
 * @param[out]  Cnb     Output DCM attitude at last moment
 * @return 0: Ok
 * @warning     1. (imu->n - 1)/(Nveb_n - 1) should be an interger
 *              2. Nveb_n >= 4(shouldn't be too small)
 *              3. imu data should be uniform sampling
 * @see align_coarse_static_base()
 * @see align_coarse_inertial()
 * @note Ref:
 *  Peter M.G. Silson, Coarse Align of Ship's Strapdown Inertial Attitude
 *  Reference System Using Velocity Loci, 2011
 *  F. Landis Markley, Attitude Determination using Vector Observations and the
 *  Singular Value Decompostion, 1988
 */
extern int align_coarse_wuhba(const imu_t *imu, double lat, const v3_t *veb_n,
        unsigned int Nveb_n, m3_t *Cnb)
{
    /* N-frame: inertial frame of n-frame at start moment */
    /* B-frame: inertial frame of b-frame at start moment */
    double ts = timediff(imu->data[1].time,imu->data[0].time);
    unsigned int len_dv = (imu->n - 1) / (unsigned int)(Nveb_n - 1);
    double nts = len_dv * ts;

    double sin_lat = sin(lat), cos_lat = cos(lat);
    v3_t gn; gravity_ned(lat,0.0,&gn);
    v3_t wie_n = {wgs84.wie*cos_lat, 0.0, -wgs84.wie*sin_lat};
    v3_t dtheta_ie_n = v3_scalar(nts,wie_n);

    /* Calculate vib_B1,vib_B2 */
    v3_t vib_B = {0.0, 0.0, 0.0};
    m3_t CbB = I3, CnN = I3;   /* initial attitde at start moment*/
    m3_t Ck_k1; /* trans from k to k+1 */
    /* int ind_mid = (imu->n) / 2; */

    v3_t *dv_N = (v3_t *)malloc(sizeof(v3_t) * (Nveb_n - 1));
    v3_t *dv_B = (v3_t *)malloc(sizeof(v3_t) * (Nveb_n - 1));

    /* TN, TN_last: wie_N x veb_N - gN */
    v3_t veb_N, veb_N_last = V0, TN, TN_last = V0;
    v3_t mean_v, wen_n, dtheta_Nn_n;

    for (unsigned int i = 1, n = 0; i < imu->n; ++i) {
        /* i start with 1: becase imu->data[0](first obs) should NOT be counted */
        /* fib_B integration, ref: Peter, 2011, eq.7 */
        vib_B = v3_add(vib_B,m3_mul_v3(CbB, imu->data[i].accel));
        /* CbB attitude update uner inertial frame */
        rv2dcm(&imu->data[i].gyro,&Ck_k1);
        CbB = m3_mul(CbB,Ck_k1);

        if(i%len_dv == 0){  /* Save vib_B and vib_N, then reset */
            dv_B[n] = vib_B; vib_B = (v3_t){0.0,0.0,0.0};

            if(n == 0){
                veb_N_last = veb_n[0];
                TN_last = v3_del(v3_cross(wie_n,veb_n[0]),gn);
            }
            mean_v = v3_scalar(0.5,v3_add(veb_n[n],veb_n[n+1]));
            wen_n = (v3_t){ mean_v.y / wgs84.R0,  - mean_v.x / wgs84.R0,
                - mean_v.y * tan(lat) / wgs84.R0};
            /* update CnN */
            dtheta_Nn_n = v3_add(dtheta_ie_n,v3_scalar(nts, wen_n));
            rv2dcm(&dtheta_Nn_n,&Ck_k1);
            CnN = m3_mul(CnN,Ck_k1);
            /* Calculate  dv_N  and save */
            veb_N = m3_mul_v3(CnN,veb_n[n+1]);
            TN = v3_del(v3_cross(m3_mul_v3(CnN,wie_n), m3_mul_v3(CnN,veb_n[n+1])),
                    m3_mul_v3(CnN,gn));
            dv_N[n] = v3_add(v3_del(veb_N, veb_N_last),
                    v3_scalar(0.5*nts, v3_add(TN, TN_last)));

            TN_last = TN; veb_N_last = veb_N;
            n++;
        }
    }

    /* interleave accumulate */
    len_dv = Nveb_n / 2;
    unsigned int N_dvsum = Nveb_n - len_dv;  /* Number of vectors to determiate attitude */
    v3_t *dv_N_sum = (v3_t *)malloc(sizeof(v3_t) * N_dvsum);
    v3_t *dv_B_sum = (v3_t *)malloc(sizeof(v3_t) * N_dvsum);
    for(unsigned int i = 0; i < N_dvsum; ++i){
        dv_N_sum[i] = (v3_t){0.0, 0.0, 0.0};
        dv_B_sum[i] = (v3_t){0.0, 0.0, 0.0};
        for(unsigned int j = 0; j < len_dv; ++j){
            dv_N_sum[i] = v3_add(dv_N_sum[i], dv_N[i+j]);
            dv_B_sum[i] = v3_add(dv_B_sum[i], dv_B[i+j]);
        }
        /* Normalize to unit vector */
        v3_normalize(&dv_N_sum[i]);
        v3_normalize(&dv_B_sum[i]);
    }

    /* Solve wuhba problem by using SVD solution */
    m3_t B = O3;
    for(unsigned int i = 0; i < N_dvsum; ++i){
        B = m3_add(B,v3_mul_cxr(dv_N_sum[i],dv_B_sum[i]));
    }
    m3_t U,V; v3_t D;
    m3_SVD(&B,&U,&D,&V);
    double d = m3_det(&U) * m3_det(&V); /* d = +-1 */
    m3_t CBN_opt = m3_mul(U,m3_mul(v3_diag((v3_t){1.0, 1.0, d}),m3_T(V)));

    /* Calculate variance-covariance */
    /* Beasue of sensors' bias, difference between CBN_opt and CBN_true is a
     * systemic bias, could not evaluating by RMSE */
    /* if(Q_Enb != NULL){
        v3_t diff_dv = {0.0}, SS_diff_dv = {0.0};
        m3_t AN = {0.0};
        for(int i = 0; i < N_dvsum; ++i){
            diff_dv = v3_del(dv_N_sum[i],m3_mul_v3(CBN_opt,dv_B_sum[i]));
            SS_diff_dv = v3_add(SS_diff_dv, v3_pow(diff_dv,2.0));
            AN = m3_add(AN,v3_mul_cxr(dv_N_sum[i],dv_N_sum[i]));
        }
        v3_t var_dv_v3 = v3_scalar(1.0/(Nveb_n - len_dv),SS_diff_dv);
        double var_dv = SQR(sqrt(var_dv_v3.i) + sqrt(var_dv_v3.j) + sqrt(var_dv_v3.k)) / 9.0;
        AN = m3_scalar(1.0/N_dvsum, AN);
        AN.m11 = 1 - AN.m11; AN.m22 = 1 - AN.m22; AN.m33 = 1 - AN.m33;
        m3_inv(&AN);
        *Q_Enb = m3_scalar(var_dv/N_dvsum, AN);
        *Q_Enb = m3_transpose(m3_mul(m3_mul(m3_transpose(CnN),*Q_Enb),CnN));
    }*/

    /* Caludate Cnb at last moment: Cbn = CNn * CBN * CbB */
    *Cnb = m3_T(m3_mul(m3_mul(m3_T(CnN),CBN_opt),CbB));

    free(dv_N); free(dv_B); free(dv_N_sum); free(dv_B_sum);

    return 0;
}

/**
 * @brief Form state transformation jacobi matrix from Ebe to Ebe
 * @param[out]  F   jacobi matrix
 * @param[in]   dt  time interval[s]
 * @return O: OK
 * @note This function relate to the earth rotation rate, use wgs84 parameter.
 *
 *      Ref: Paul. D. Groves. Principles of GNSS, Inertial, and Multisensor
 *      Integrated Navigation Systems(2nd Edition), P583, 14.48
 */
extern inline int jacobi_trans_Ebe2Ebe(m3_t *F, double dt)
{
    /* I3 - OMGie_e * dt */
    v3_t omg_ie = {0.0, 0.0, wgs84.wie};
    omg_ie = v3_scalar(-dt, omg_ie);
    asymmetric_mat(&omg_ie, F);
    F->m11 += 1; F->m22 += 1; F->m33 += 1;
    return 0;
}

/**
 * @brief Form state transformation jacobi matrix from Ebe to gryo bias
 * @param[out]  F   Jacobi matrix
 * @param[in]   dt  Time interval
 * @param[in]   Cbe Current average attitude
 * @return 0: OK
 * @note
 *      Ref: Paul. D. Groves. Principles of GNSS, Inertial, and Multisensor
 *      Integrated Navigation Systems(2nd Edition), P583, 14.48
 */
extern inline int jacobi_trans_Ebe2bg(m3_t *F, double dt, const m3_t *Cbe)
{
    /* - Cbe * dt */
    *F = m3_scalar(-dt, *Cbe);
    return 0;
}

/**
 * @brief Form state transformation jacobi matrix from veb_e to Ebe
 * @param[out]  F   Jacobi matrix
 * @param[in]   Cbe Current average attitude.
 * @param[in]   dv  Velocity increment[m/s]
 * @return 0: OK
 * @note
 *      Ref: Paul. D. Groves. Principles of GNSS, Inertial, and Multisensor
 *      Integrated Navigation Systems(2nd Edition), P583, 14.48
 */
extern inline int
jacobi_trans_veb_e2Ebe(m3_t *F, const m3_t *Cbe, const v3_t *dv)
{
    /* [- Cbe * fib_b ]x * dt */
    v3_t Tib_e = v3_scalar(-1.0, m3_mul_v3(*Cbe, *dv));
    asymmetric_mat(&Tib_e, F);
    return 0;
}

/**
 * @brief Form state transformation jacobi matrix from veb_e to veb_e
 * @param[out]  F   Jacobi matrix
 * @param[in]   dt  Time interval[s]
 * @return 0: OK
 * @note
 *      Ref: Paul. D. Groves. Principles of GNSS, Inertial, and Multisensor
 *      Integrated Navigation Systems(2nd Edition), P583, 14.48
 */
extern inline int jacobi_trans_veb_e2veb_e(m3_t *F, double dt)
{
    /* I3 - 2 * OMGie_e * dt */
    v3_t omg_ie = {0.0, 0.0, wgs84.wie};
    omg_ie = v3_scalar(-2 * dt, omg_ie);
    asymmetric_mat(&omg_ie, F);
    F->m11 += 1; F->m22 += 1; F->m33 += 1;
    return 0;
}

/**
 * @brief Form state transforamtion jacobi matrix from veb_e to reb_e
 * @param[out]  F       Jacobi matrix
 * @param[in]   dt      Time interval[s]
 * @param[in]   reb_e   ecef position[m]
 * @return 0: OK
 * @note
 *      Ref: Paul. D. Groves. Principles of GNSS, Inertial, and Multisensor
 *      Integrated Navigation Systems(2nd Edition), P583, 14.48
 */
extern inline int
jacobi_trans_veb_e2reb_e(m3_t *F, double dt, const v3_t *reb_e)
{
    /* F23 * dt */
    /* Ref Paul P71:2.137 */
    v3_t pos = *reb_e; ecef2ned(&pos, NULL, NULL);
    double RE = earth_RE(&wgs84,pos.x);
    double grtmp = SQR(cos(pos.x)) + SQR(1-SQR(wgs84.e)) * SQR(sin(pos.x));
    double geocentric_radius = RE * sqrt(grtmp);
    v3_t ge; gravity_ecef(reb_e, &ge);
    double fac = - dt * 2 / geocentric_radius / v3_norm(*reb_e);
    ge = v3_scalar(fac, ge);
    *F = v3_mul_cxr(ge, *reb_e);
    return 0;
}

/**
 * @brief Form state transformation jacobi matrix from veb_e to ba
 * @param[out]  F       Jacobi matrix
 * @param[in]   dt      Time interval[s]
 * @param[in]   Cbe     Current average attitude over time[m]
 * @return 0: OK
 * @note
 *      Ref: Paul. D. Groves. Principles of GNSS, Inertial, and Multisensor
 *      Integrated Navigation Systems(2nd Edition), P583, 14.48
 */
extern inline int jacobi_trans_veb_e2ba(m3_t *F, double dt, const m3_t *Cbe)
{
    *F = m3_scalar(-dt, *Cbe);
    return 0;
}

/**
 * @brief Form state transformation jacobi matrix from reb_e to veb_e
 * @param[out]  F   Jacobi matrix
 * @param[in]   dt  Time interval[s]
 * @return 0: OK
 * @note
 *      Ref: Paul. D. Groves. Principles of GNSS, Inertial, and Multisensor
 *      Integrated Navigation Systems(2nd Edition), P583, 14.48
 */
extern inline int jacobi_trans_reb_e2veb_e(m3_t *F, double dt) {
    F->m11 = dt; F->m12 = 0.0; F->m13 = 0.0;
    F->m21 = 0.0; F->m22 = dt; F->m23 = 0.0;
    F->m31 = 0.0; F->m32 = 0.0; F->m33 = dt;
    return 0.0;
}

/**
 * @brief Form state transformation jacobi matrix of
 *      1st order markov/random walk process
 * @param[out]   F      Jacobi matrix
 * @param[in]   dt      Time interval[s]
 * @param[in]   T       correlation time[s]
 * @warning dt should far less than T(dt << T)
 * @note if T == V0, this function will return random walk process factor(I3)
 *  markov process should not used if you can not comprehend what really it is.
 */
extern inline int jacobi_trans_markov(m3_t *F, double dt, const v3_t *T)
{
    F->m11 = T->x == 0.0 ? 1.0 : 1.0 - dt / T->x;
    F->m22 = T->y == 0.0 ? 1.0 : 1.0 - dt / T->y;
    F->m33 = T->z == 0.0 ? 1.0 : 1.0 - dt / T->z;
    F->m12 = 0.0; F->m13 = 0.0;
    F->m21 = 0.0; F->m23 = 0.0;
    F->m31 = 0.0; F->m32 = 0.0;
    return 0;
}

/**
 * @brief Form measurement jacobi matrix from reb_e to Ebe
 * @param[out]  H           Jacobi matrix
 * @param[in]   Cbe         Current average attitude
 * @param[in]   lever_arm   gnss to imu lever arm[m]
 * @return status(0: OK)
 * @note
 *      Ref: Paul. D. Groves. Principles of GNSS, Inertial, and Multisensor
 *      Integrated Navigation Systems(2nd Edition), P600, 14.111
 */
extern int jacobi_meas_reb_e2Ebe(m3_t *H, const m3_t *Cbe, const v3_t *lever_arm)
{
    v3_t v = m3_mul_v3(*Cbe, *lever_arm);
    asymmetric_mat(&v, H);
    return 0;
}

/**
 * @brief Form measurement jacobi matrix from veb_e to Ebe
 * @param[out]  H			Jacobi matrix
 * @param[in]   Cbe			Current average attitude
 * @param[in]   wib_b		imu output angluar rate[rad/s]
 * @param[in]   lever_arm	gnss to imu lever arm[m]
 * @return status(0: OK)
 * @note
 *      Ref: Paul. D. Groves. Principles of GNSS, Inertial, and Multisensor
 *      Integrated Navigation Systems(2nd Edition), P600, 14.111
 */
extern int jacobi_meas_veb_e2Ebe(m3_t *H, const m3_t *Cbe, const v3_t *wib_b,
    const v3_t *lever_arm)
{
    v3_t v1 = m3_mul_v3(*Cbe, v3_cross(*wib_b, *lever_arm));
    v3_t wie_e = {0.0, 0.0, wgs84.wie};
    v3_t v2 = v3_cross(wie_e, m3_mul_v3(*Cbe, *lever_arm));
    v3_t v = v3_del(v1, v2);
    asymmetric_mat(&v, H);
    return 0;
}

/**
 * @brief Form measurement jacobi matrix from veb_e to bg
 * @param[out] 		H			Jacobi matrix
 * @param[in] 		Cbe			Current average attitude
 * @param[in] 		lever_arm	gnss lever arm[m]
 * @return status(0: 0K)
 * @note
 *      Ref: Paul. D. Groves. Principles of GNSS, Inertial, and Multisensor
 *      Integrated Navigation Systems(2nd Edition), P600, 14.111
 */
extern int jacobi_meas_veb_e2bg(m3_t *H, const m3_t *Cbe, const v3_t *lever_arm)
{
    m3_t lx; asymmetric_mat(lever_arm, &lx);
    *H = m3_mul(*Cbe, lx);
    return 0;
}

/**
 * @brief Form state transformation matrix from veb_e to ka.x
 * @param[out] 	F		Jacobi matrix
 * @param[in] 	Cbe		current average attitude
 * @param[in] 	dv		imu output velocity increment
 * @return status(0: OK)
 */
extern inline int
jacobi_trans_veb_e2kax(double *F, const m3_t *Cbe, const v3_t *dv)
{
    F[cfg.IVEL  +cfg.IKAx*cfg.nx] = Cbe->m11 * dv->x;
    F[cfg.IVEL+1+cfg.IKAx*cfg.nx] = Cbe->m21 * dv->x;
    F[cfg.IVEL+2+cfg.IKAx*cfg.nx] = Cbe->m31 * dv->x;
    return 0;
}

/**
 * @brief Form state transformation matrix from veb_e to ka.y
 * @param[out] 	F		Jacobi matrix
 * @param[in] 	Cbe		current avearge attitude
 * @param[in] 	dv		imu output velocity increment
 * @return status(0: OK)
 */
extern inline int
jacobi_trans_veb_e2kay(double *F, const m3_t *Cbe, const v3_t *dv)
{
    F[cfg.IVEL  +cfg.IKAy*cfg.nx] = Cbe->m12 * dv->y;
    F[cfg.IVEL+1+cfg.IKAy*cfg.nx] = Cbe->m22 * dv->y;
    F[cfg.IVEL+2+cfg.IKAy*cfg.nx] = Cbe->m32 * dv->y;
    return 0;
}

/**
 * @brief Form state transformation matrix from veb_e to ka.z
 * @param[out] 	F		Jacobi matrix
 * @param[in] 	Cbe		current average attitude
 * @param[in] 	dv		imu output velocity increment
 * @return status(0: OK)
 */
extern inline int
jacobi_trans_veb_e2kaz(double *F, const m3_t *Cbe, const v3_t *dv)
{
    F[cfg.IVEL  +cfg.IKAz*cfg.nx] = Cbe->m13 * dv->z;
    F[cfg.IVEL+1+cfg.IKAz*cfg.nx] = Cbe->m23 * dv->z;
    F[cfg.IVEL+2+cfg.IKAz*cfg.nx] = Cbe->m33 * dv->z;
    return 0;
}

/**
 * @brief Form state transformation matrix from veb_e to kg.x
 * @param[out] 	F		Jacobi matrix
 * @param[in] 	Cbe		current average attitude
 * @param[in] 	dtheta	imu output angular increment
 * @return status(0: OK)
 */
extern inline int
jacobi_trans_Ebe2kgx(double *F, const m3_t *Cbe, const v3_t *dtheta)
{
    F[cfg.IATT  +cfg.IKGx*cfg.nx] = Cbe->m11 * dtheta->x;
    F[cfg.IATT+1+cfg.IKGx*cfg.nx] = Cbe->m21 * dtheta->x;
    F[cfg.IATT+2+cfg.IKGx*cfg.nx] = Cbe->m31 * dtheta->x;
    return 0;
}

/**
 * @brief Form state transformation matrix from veb_e to kg.y
 * @param[out] 	F		Jacobi matrix
 * @param[in] 	Cbe		current average attitude
 * @param[in] 	dtheta	imu output angular increment
 * @return status(0: OK)
 */
extern inline int
jacobi_trans_Ebe2kgy(double *F, const m3_t *Cbe, const v3_t *dtheta)
{
    F[cfg.IATT  +cfg.IKGy*cfg.nx] = Cbe->m12 * dtheta->y;
    F[cfg.IATT+1+cfg.IKGy*cfg.nx] = Cbe->m22 * dtheta->y;
    F[cfg.IATT+2+cfg.IKGy*cfg.nx] = Cbe->m32 * dtheta->y;
    return 0;
}

/**
 * @brief Form state transformation matrix from veb_e to kg.z
 * @param[out] 	F		Jacobi matrix
 * @param[in] 	Cbe		current average attitude
 * @param[in] 	dtheta	imu output angular increment
 * @return status(0: OK)
 */
extern inline int
jacobi_trans_Ebe2kgz(double *F, const m3_t *Cbe, const v3_t *dtheta)
{
    F[cfg.IATT  +cfg.IKGz*cfg.nx] =  Cbe->m13 * dtheta->z;
    F[cfg.IATT+1+cfg.IKGz*cfg.nx] =  Cbe->m23 * dtheta->z;
    F[cfg.IATT+2+cfg.IKGz*cfg.nx] =  Cbe->m33 * dtheta->z;
    return 0;
}

/**
 * @brief dead reckoning on navigation
 * @param[in] dt        time interval[s]
 * @param[in] dtheta    angular increment in the time interval dt [rad]
 * @param[in] dS        distance increment in the time interval dt [m]
 * @param[in,out] r     last/current position under e-frame[m]
 * @param[in,out] q     last/current attitude exepress by quaternion(Qbe)
 * @return status(0: OK)
 */
extern int dr_nav_ecef(double dt, const v3_t* dtheta, double dS, v3_t *r,
                       quat_t *q)
{
    /* position update */
    *r = v3_add(*r, quat_mul_v3(*q, (v3_t){dS, 0.0, 0.0}));
    /* Attitude update */
    v3_t dtheta_ie = {0.0, 0.0, -wgs84.wie * dt};
    quat_t q_earth, q_body, old_q = *q;
    rv2quat(&dtheta_ie, &q_earth);
    rv2quat(dtheta, &q_body);
    *q = quat_mul(quat_mul(q_earth, old_q), q_body);
    quat_normalize(q);
    return 0;
}
