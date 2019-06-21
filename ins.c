/**
 * @file ins.c
 * @brief ins core functions implementation
 * @author yinflying(yinflying@foxmail.com)
 */
#include "ins.h"

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define PI 3.14159265358979
#define SQR(x) (x) * (x)

static const m3_t I3 = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
static const m3_t O3 = {0.0};

/**
 * @brief Gravitational acceleration of Earth project to e-axis
 * @param[in]   r   postion under e-axis
 * @param[out]  ge  Gravitational acceleration under ECEF(m s^-2)
 * @return 0: OK
 * @see gravity_ned()
 * @note Do not contain centrifugal force
 */
int gravity_ecef(const v3_t* r, v3_t* ge)
{
    double mag_r = sqrt(r->i * r->i + r->j * r->j + r->k * r->k);
    if (mag_r == 0) {
        ge->i = 0, ge->j = 0, ge->k = 0;
    } else {
        /* Calulate gravitational acceleration using (2.142) */
        double f1 = -wgs84.mu / pow(mag_r, 3);
        double f2 = 1.5 * wgs84.J2 * pow(wgs84.R0 / mag_r, 2);
        double z_scale = 5 * pow((r->k / mag_r), 2);
        double g1 = f1 * (r->i + f2 * (1 - z_scale) * r->i);
        double g2 = f1 * (r->j + f2 * (1 - z_scale) * r->j);
        double g3 = f1 * (r->k + f2 * (3 - z_scale) * r->k);
        /* Add centripetal acceleration using (2.133) */
        ge->i = g1 + wgs84.wie * wgs84.wie * r->i;
        ge->j = g2 + wgs84.wie * wgs84.wie * r->j;
        ge->k = g3;
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
int gravity_ned(double lat, double hgt, v3_t* gn)
{
    double sinlat2 = sin(lat) * sin(lat);
    double e2 = SQR(wgs84.e);

    /* Calculate surface gravity using Somigliana model */
    double g0 = 9.7803253359 * (1.0 + 0.001931853 * sinlat2)
        / sqrt(1.0 - e2 * sinlat2);

    gn->i = -8.08E-9 * hgt * sin(2.0 * lat); /* North */
    gn->j = 0.0;                           /* East */
    /* Down */
    double tmp = 1.0 + wgs84.f * (1.0 - 2.0 * sinlat2)
        + (SQR(wgs84.wie) * SQR(wgs84.R0) * wgs84.RP / wgs84.mu);
    gn->k = g0 *
        (1.0 - (2.0 / wgs84.R0) * tmp * hgt + (3.0 * SQR(hgt) / SQR(wgs84.R0)));

    return 0;
}

/**
 * @brief Strapdown-INS equations under ECEF frame
 * @param[in]   dt      Time interval [s]
 * @param[in]   dtheta  Angular increment [rad]
 * @param[in]   dv      Velocity increment [m/s]
 * @param[in,out]   r   Start/End postion in ECEF [m]
 * @param[in,out]   v   Start/End velocity in ECEF [m]
 * @param[in,out]   q   Start/End attitude trans express by quaternion(qbe)
 * @return 0: OK
 * @see multisample()
 * @note This function do not contain conning&sculling error compensation, so
 *      so it should work with multisample() function.
 *
 *      Ref: Paul. D. Groves. Principles of GNSS, Inertial, and Multisensor
 *      Integrated Navigation Systems(2nd Edition), P163
 */
extern int nav_equations_ecef(
    double dt, const v3_t* dtheta, const v3_t* dv, v3_t* r, v3_t* v, quat_t* q)
{
    /* Attitude update */
    v3_t dtheta_ie = { 0, 0, -wgs84.wie * dt };
    quat_t q_earth;
    rv2quat(&dtheta_ie, &q_earth);
    quat_t q_body;
    rv2quat(dtheta, &q_body);
    quat_t old_q = *q;
    *q = quat_mul(quat_mul(q_earth, old_q), q_body);
    quat_normalize(q);
    /* Specific force transform(velocity form) */
    v3_t dtheta_ie_half = { 0, 0, -wgs84.wie * dt / 2 };
    rv2quat(&dtheta_ie_half, &q_earth);
    v3_t dv_rot = v3_dot(0.5, v3_cross(*dtheta, *dv));
    v3_t dv_e = quat_mul_v3(quat_mul(q_earth, old_q), v3_add(*dv, dv_rot));
    /* Velocity update */
    v3_t ge;
    gravity_ecef(r, &ge);
    v3_t old_v = *v;
    v->i = old_v.i + dv_e.i + dt * (ge.i + 2 * wgs84.wie * old_v.j);
    v->j = old_v.j + dv_e.j + dt * (ge.j - 2 * wgs84.wie * old_v.i);
    v->k = old_v.k + dv_e.k + dt * ge.k;
    /* Position update */
    v3_t old_r = *r;
    *r = v3_add(old_r, v3_dot(0.5 * dt, v3_add(old_v, *v)));
    return 0;
}

/**
 * @brief Use multi-subsample to compensate the conning&scull error
 * @param[in]   dtheta_list     Angular increment list,order by time,length:abs(N)
 * @param[in]   dv_list         Velocity increment list,order by time,length:abs(N)
 * @param[in]   N               Subsample number(1=<N<=5, -2: one-plus-previous)
 * @param[out]  dtheta          Sum of angular increment with conning error compensation
 * @param[out]  dv              Sum of velocity incrment with scull error compensation
 * @return 0: No error    1: Error
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
        v3_t sum_c = {}, sum_s = {}, sum_w = {}, sum_v = {};
        int i = 0;
        for (; i < N - 1; i++) {
            /* sum of cross_product factor */
            sum_c = v3_add(sum_c, v3_dot(pcf[i], dtheta_list[i]));
            sum_s = v3_add(sum_s, v3_dot(pcf[i], dv_list[i]));
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
        v3_t sum_c = v3_dot(1.0 / 12, dtheta_list[0]);
        v3_t sum_s = v3_dot(1.0 / 12, dv_list[0]);
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

/* Double vector to Atttitude
 * Cnb => the matix tranform n-axis to b-axis(or attitude b-axis to n-axis)
 * */
int dblvec2att(const v3_t* vn1, const v3_t* vn2, const v3_t* vb1,
    const v3_t* vb2, m3_t* Cnb)
{
    /* ref Yan2016(P143:6.1-7) */
    v3_t vx = v3_dot(1.0 / v3_norm(*vb1), *vb1);
    v3_t vtmp = v3_cross(*vb1, *vb2);
    v3_t vy = v3_dot(1.0 / v3_norm(vtmp), vtmp);
    vtmp = v3_cross(vtmp, *vb1);
    v3_t vz = v3_dot(1.0 / v3_norm(vtmp), vtmp);
    m3_t m3_b = { vx.i, vy.i, vz.i, vx.j, vy.j, vz.j, vx.k, vy.k, vz.k };

    vx = v3_dot(1.0 / v3_norm(*vn1), *vn1);
    vtmp = v3_cross(*vn1, *vn2);
    vy = v3_dot(1.0 / v3_norm(vtmp), vtmp);
    vtmp = v3_cross(vtmp, *vn1);
    vz = v3_dot(1.0 / v3_norm(vtmp), vtmp);
    m3_t m3_n = { vx.i, vx.j, vx.k, vy.i, vy.j, vy.k, vz.i, vz.j, vz.k };

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
    for (int i = 1; i < imu->n; ++i) {
        mean_wib_b = v3_add(mean_wib_b, imu->data[i].gryo);
        mean_fib_b = v3_add(mean_fib_b, imu->data[i].accel);
    }
    double T = yins_timediff(imu->data[imu->n - 1].time, imu->data[0].time);
    mean_wib_b = v3_dot(1.0 / T, mean_wib_b);
    mean_fib_b = v3_dot(1.0 / T, mean_fib_b);

    /* Double Vector to Attitude */
    v3_t gn; gravity_ned(lat,0.0,&gn);
    v3_t wie_n = {wgs84.wie * cos(lat),0.0,- wgs84.wie * sin(lat)};
    v3_t gib_b = v3_dot(-1.0,mean_fib_b);

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
    /* N-frame: inertial frame of n-frame at start moment */
    /* B-frame: inertial frame of b-frame at start moment */

    int sample_N = 4;
    double ts = yins_timediff(imu->data[1].time,imu->data[0].time);
    double nts = sample_N * ts;

    double sin_lat = sin(lat), cos_lat = cos(lat);
    v3_t gn; gravity_ned(lat,0.0,&gn);

    /* Calculate vib_B1,vib_B2*/
    v3_t dtheta[4],dv[4],sum_dtheta,sum_dv;
    v3_t vib_B = {0.0}, vib_B1 = {0.0};
    quat_t qb_B = {1.0, 0.0, 0.0, 0.0};   /* initial attitde*/
    quat_t qk_k1; /* trans from k to k+1 */
    int ind_mid = (imu->n/sample_N) / 2 * sample_N;
    for (int i = 0; i <= imu->n - sample_N ; i += sample_N) {
        for(int j = 0; j < sample_N; ++j){
            dtheta[j] = imu->data[i+j].gryo;
            dv[j] = imu->data[i+j].accel;
        }
        /* Calculate current fib_B */
        multisample(dtheta,dv,sample_N,&sum_dtheta,&sum_dv);
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
    double gcl_wie = gn.k * cos_lat / wgs84.wie;
    v3_t vib_I1 = { gcl_wie * sin(wie_dtheta/2.0),
        gcl_wie * (1-cos(wie_dtheta/2.0)), total_t/2.0 * gn.k * sin_lat};
    v3_t vib_I2 = { gcl_wie * sin(wie_dtheta),
        gcl_wie * (1-cos(wie_dtheta)), total_t * gn.k * sin_lat};

    /* double vector to attitude */
    m3_t CB_I; dblvec2att(&vib_B1,&vib_B,&vib_I1,&vib_I2,&CB_I);

    /* Calculate Cnb */
    double cos_wie = cos(wie_dtheta), sin_wie = sin(wie_dtheta);
    m3_t CI_n = { -sin_lat * cos_wie, -sin_lat * sin_wie, cos_lat,
        -sin_wie, cos_wie, 0.0,
        -cos_lat * cos_wie, -cos_lat * sin_wie, -sin_lat };
    m3_t Cb_B; quat2dcm(&qb_B,&Cb_B);
    *Cnb = m3_transpose(m3_mul(m3_mul(CI_n,CB_I),Cb_B));

    return 0;
}

/**
 * @brief Coarse Alignment by solving Wuhba problem under inertial frame
 * @param[in]   imu     Inertial IMU
 * @param[in]   lat     Imu latitude(Average Latitude) [rad]
 * @param[in]   veb_n   Imu velocity uner n-frame, length: Nveb_n [m/s]
 * @param[in]   Nveb_n  Imu velocity number
 * @param[out]  Cnb     Output DCM attitude at last moment
 * @return 0: Ok
 * @warning     1. (imu->n - 1)/(Nveb_n - 1) should be an interger
 *              2. Nveb_n >= 3(shouldn't be too small)
 * @note Ref:
 *  Peter M.G. Silson, Coarse Align of Ship's Strapdown Inertial Attitude
 *  Reference System Using Velocity Loci, 2011
 *  F. Landis Markley, Attitude Determination using Vector Observations and the
 *  Singular Value Decompostion, 1988
 */
extern int align_coarse_wuhba(
        const imu_t *imu, double lat, const v3_t *veb_n, int Nveb_n, m3_t *Cnb)
{
    /* N-frame: inertial frame of n-frame at start moment */
    /* B-frame: inertial frame of b-frame at start moment */
    double ts = yins_timediff(imu->data[1].time,imu->data[0].time);
    int len_dv = (imu->n - 1) / (Nveb_n - 1);
    double nts = len_dv * ts;

    double sin_lat = sin(lat), cos_lat = cos(lat);
    v3_t gn; gravity_ned(lat,0.0,&gn);
    v3_t wie_n = {wgs84.wie*cos_lat, 0.0, -wgs84.wie*sin_lat};
    v3_t dtheta_ie_n = v3_dot(nts,wie_n);

    /* Calculate vib_B1,vib_B2 */
    v3_t vib_B = {0.0}, vib_N = {0.0};
    m3_t CbB = I3, CnN = I3;   /* initial attitde at start moment*/
    m3_t Ck_k1; /* trans from k to k+1 */
    int ind_mid = (imu->n) / 2;

    v3_t *dv_N = (v3_t *)malloc(sizeof(v3_t) * (Nveb_n - 1));
    v3_t *dv_B = (v3_t *)malloc(sizeof(v3_t) * (Nveb_n - 1));

    /* TN, TN_last: wie_N x veb_N - gN */
    v3_t veb_N, veb_N_last, TN, TN_last;
    v3_t mean_v, wen_n, dtheta_Nn_n;

    for (int i = 1, n = 0; i < imu->n; ++i) {
        /* i start with 1: becase imu->data[0](first obs) should NOT be counted */
        /* Calculate current fib_B */
        vib_B = v3_add(vib_B,m3_mul_v3(CbB, imu->data[i].accel));
        /* qb_B attitude update uner inertial frame */
        rv2dcm(&imu->data[i].gryo,&Ck_k1);
        CbB = m3_mul(CbB,Ck_k1);

        if(i%len_dv == 0){  /* Save vib_B and vib_N reset */
            dv_B[n] = vib_B; vib_B = (v3_t){0.0,0.0,0.0};

            if(n == 0){
                veb_N_last = veb_n[0];
                TN_last = v3_del(v3_cross(wie_n,veb_n[0]),gn);
            }
            /* update CnN */
            mean_v = v3_dot(0.5,v3_add(veb_n[n],veb_n[n+1]));
            wen_n = (v3_t){ mean_v.j / wgs84.R0,  - mean_v.i / wgs84.R0,
                - mean_v.j * tan(lat) / wgs84.R0};
            dtheta_Nn_n = v3_add(dtheta_ie_n,v3_dot(nts, wen_n));
            rv2dcm(&dtheta_Nn_n,&Ck_k1);
            CnN = m3_mul(CnN,Ck_k1);
            /* Calculate  dv_N  and save */
            veb_N = m3_mul_v3(CnN,veb_n[n+1]);
            TN = v3_del(v3_cross(m3_mul_v3(CnN,wie_n), m3_mul_v3(CnN,veb_n[n+1])),
                    m3_mul_v3(CnN,gn));
            dv_N[n] = v3_add(v3_del(veb_N, veb_N_last),
                    v3_dot(0.5*nts, v3_add(TN, TN_last)));

            TN_last = TN; veb_N_last = veb_N;
            n++;
        }
    }

    /* interleave accumulate */
    len_dv = Nveb_n / 2;
    v3_t *dv_N_sum = (v3_t *)malloc(sizeof(v3_t) * (Nveb_n - len_dv));
    v3_t *dv_B_sum = (v3_t *)malloc(sizeof(v3_t) * (Nveb_n - len_dv));
    for(int i = 0; i < Nveb_n - len_dv; ++i){
        dv_N_sum[i] = (v3_t){0.0};
        dv_B_sum[i] = (v3_t){0.0};
        for(int j = 0; j < len_dv; ++j){
            dv_N_sum[i] = v3_add(dv_N_sum[i], dv_N[i+j]);
            dv_B_sum[i] = v3_add(dv_B_sum[i], dv_B[i+j]);
        }
        /* Nomalize to unit vector*/
        v3_normalize(&dv_N_sum[i]);
        v3_normalize(&dv_B_sum[i]);
    }

    /* Solve wuhba problem using SVD solution*/
    m3_t B = {0.0};
    for(int i = 0; i < Nveb_n - len_dv; ++i){
        B = m3_add(B,v3_mul_cxr(dv_N_sum[i],dv_B_sum[i]));
    }
    m3_t U,V; v3_t D;
    m3_SVD(&B,&U,&D,&V);
    double d = m3_det(&U) * m3_det(&V); /* d = +-1 */
    m3_t CBN_opt = m3_mul(U,m3_mul(v3_diag((v3_t){1.0, 1.0, d}),m3_transpose(V)));

    /* Caludate Cnb at last moment: Cbn = CNn * CBN * CbB */
    *Cnb = m3_transpose(m3_mul(m3_mul(m3_transpose(CnN),CBN_opt),CbB));

    return 0;
}
