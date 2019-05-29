#include "ins.h"

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define PI 3.14159265358979
#define SQR(x) (x) * (x)

extern v3_t v3_cross(v3_t v1, v3_t v2)
{
    v3_t v;
    v.i = v1.j * v2.k - v1.k * v2.j;
    v.j = v1.k * v2.i - v1.i * v2.k;
    v.k = v1.i * v2.j - v1.j * v2.i;
    return v;
}
extern v3_t v3_add(v3_t v1, v3_t v2)
{
    return (v3_t) { v1.i + v2.i, v1.j + v2.j, v1.k + v2.k };
}
extern v3_t v3_del(v3_t v1, v3_t v2)
{
    return (v3_t) { v1.i - v2.i, v1.j - v2.j, v1.k - v2.k };
}
extern v3_t v3_dot(double s, v3_t v)
{
    return (v3_t) { s * v.i, s * v.j, s * v.k };
}
extern double v3_norm(v3_t v)
{
    return sqrt(v.i * v.i + v.j * v.j + v.k * v.k);
}
extern double v3_mul_rxc(v3_t v1, v3_t v2)
{
    return v1.i * v2.i + v1.j * v2.j + v1.k * v2.k;
}
extern m3_t v3_mul_cxr(v3_t v1, v3_t v2)
{
    return (m3_t) { v1.i * v2.i, v1.i * v2.j, v1.i * v2.k, v1.j * v2.i,
        v1.j * v2.j, v1.j * v2.k, v1.k * v2.i, v1.k * v2.j, v1.k * v2.k };
}
extern m3_t v3_diag(v3_t v)
{
    return (m3_t) {v.i, 0.0, 0.0,   0.0, v.j, 0.0,   0.0, 0.0, v.k};
}

extern m3_t m3_transpose(m3_t A)
{
    m3_t dcm;
    dcm.m11 = A.m11, dcm.m12 = A.m21, dcm.m13 = A.m31;
    dcm.m21 = A.m12, dcm.m22 = A.m22, dcm.m23 = A.m32;
    dcm.m31 = A.m13, dcm.m32 = A.m23, dcm.m33 = A.m33;
    return dcm;
}
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
v3_t m3_mul_v3(m3_t A, v3_t B)
{
    v3_t C;
    C.i = A.m11 * B.i + A.m12 * B.j + A.m13 * B.k;
    C.j = A.m21 * B.i + A.m22 * B.j + A.m23 * B.k;
    C.k = A.m31 * B.i + A.m32 * B.j + A.m33 * B.k;
    return C;
}

v3_t m3_diag(m3_t A)
{
    return (v3_t){A.m11, A.m22, A.m33};
}

earth_t wgs84 = { .wie = 7.292115E-5,
    .R0 = 6378137,
    .RP= 6356752.31425,
    .mu = 3.986004418E14,
    .J2 = 1.082627E-3,
    .e = 0.0818191908425,
    .f = 1.0 / 298.257223563
};

int asymmetric_mat(const v3_t* v3, m3_t* mat)
{
    mat->m11 = 0, mat->m12 = -v3->k, mat->m13 = v3->j;
    mat->m21 = v3->k, mat->m22 = 0, mat->m23 = -v3->i;
    mat->m31 = -v3->j, mat->m32 = v3->i, mat->m33 = 0;
    return 0;
}

/* convert calendar day/time to time -------------------------------------------
 * convert calendar day/time to gtime_t struct
 * args   : double *ep       I   day/time {year,month,day,hour,min,sec}
 * return : gtime_t struct
 * notes  : proper in 1970-2037 or 1970-2099 (64bit time_t)
 *-----------------------------------------------------------------------------*/
extern gtime_t ins_epoch2time(const double* ep)
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

extern gtime_t ins_gpst2time(int week, double sec)
{
    const double gpst0[] = { 1980, 1, 6, 0, 0, 0 }; /* gps time reference */
    gtime_t t = ins_epoch2time(gpst0);

    if (sec < -1E9 || 1E9 < sec)
        sec = 0.0;
    t.time += (time_t)86400 * 7 * week + (int)sec;
    t.sec = sec - (int)sec;
    return t;
}

/* time to calendar day/time ---------------------------------------------------
 * convert gtime_t struct to calendar day/time
 * args   : gtime_t t        I   gtime_t struct
 *          double *ep       O   day/time {year,month,day,hour,min,sec}
 * return : none
 * notes  : proper in 1970-2037 or 1970-2099 (64bit time_t)
 *-----------------------------------------------------------------------------*/
extern void ins_time2epoch(gtime_t t, double* ep)
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

extern int euler2quat(const v3_t* euler, quat_t* quat)
{
    double si = sin(euler->i / 2), ci = cos(euler->i / 2);
    double sj = sin(euler->j / 2), cj = cos(euler->j / 2);
    double sk = sin(euler->k / 2), ck = cos(euler->k / 2);
    quat->q0 = ci * cj * ck + si * sj * sk;
    quat->q1 = si * cj * ck - ci * sj * sk;
    quat->q2 = ci * sj * ck + si * cj * sk;
    quat->q3 = ci * cj * sk - si * sj * ck;
    quat_inv(quat);
    return 0;
}

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

extern int quat_normalize(quat_t* quat)
{
    double nq = sqrt(quat->q0 * quat->q0 + quat->q1 * quat->q1
        + quat->q2 * quat->q2 + quat->q3 * quat->q3);
    quat->q0 /= nq;
    quat->q1 /= nq;
    quat->q2 /= nq;
    quat->q3 /= nq;
    if (quat->q0 < 0) {
        quat->q0 = -quat->q0;
        quat->q1 = -quat->q1;
        quat->q2 = -quat->q2;
        quat->q3 = -quat->q3;
    }
    return 0;
}

extern int quat_inv(quat_t* quat)
{
    quat->q1 = -quat->q1;
    quat->q2 = -quat->q2;
    quat->q3 = -quat->q3;
    return 0;
}

extern int dtheta2quat(const v3_t* dtheta, quat_t* quat)
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

extern quat_t quat_mul(quat_t P, quat_t Q)
{
    quat_t qtmp;
    qtmp.q0 = P.q0 * Q.q0 - P.q1 * Q.q1 - P.q2 * Q.q2 - P.q3 * Q.q3;
    qtmp.q1 = P.q0 * Q.q1 + P.q1 * Q.q0 + P.q2 * Q.q3 - P.q3 * Q.q2;
    qtmp.q2 = P.q0 * Q.q2 + P.q2 * Q.q0 + P.q3 * Q.q1 - P.q1 * Q.q3;
    qtmp.q3 = P.q0 * Q.q3 + P.q3 * Q.q0 + P.q1 * Q.q2 - P.q2 * Q.q1;
    return qtmp;
}

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

extern double ins_timediff(gtime_t t1, gtime_t t2)
{
    return difftime(t1.time, t2.time) + t1.sec - t2.sec;
}

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

extern int nav_equations_ecef(
    double dt, const v3_t* dtheta, const v3_t* dv, v3_t* r, v3_t* v, quat_t* q)
{
    /* Attitude update */
    v3_t dtheta_ie = { 0, 0, -wgs84.wie * dt };
    quat_t q_earth;
    dtheta2quat(&dtheta_ie, &q_earth);
    quat_t q_body;
    dtheta2quat(dtheta, &q_body);
    quat_t old_q = *q;
    *q = quat_mul(quat_mul(q_earth, old_q), q_body);
    quat_normalize(q);
    /* Specific force transform(velocity form) */
    v3_t dtheta_ie_half = { 0, 0, -wgs84.wie * dt / 2 };
    dtheta2quat(&dtheta_ie_half, &q_earth);
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

extern int multisample(
    const v3_t* dtheta_list, const v3_t* dv_list, int N, v3_t* dtheta, v3_t* dv)
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
        /* coning error compensation for angular increment, ref Paul2013(5.97)
         */
        *dtheta = v3_add(sum_w, v3_cross(sum_c, dtheta_list[i]));
        /* sculling error compensation for velocity increment, ref
         * Paul2013(5.98)*/
        /* v3_t rot = v3_dot(0.5,v3_cross(sum_w,sum_v)); */
        v3_t scul = v3_add(
            v3_cross(sum_c, dv_list[i]), v3_cross(sum_s, dtheta_list[i]));
        /* *dv = v3_add(sum_v, v3_add(rot,scul)); */
        *dv = v3_add(sum_v, scul);
    }
    if (N == -2) {
        v3_t sum_c = v3_dot(1.0 / 12, dtheta_list[0]);
        v3_t sum_s = v3_dot(1.0 / 12, dv_list[0]);
        /* ref Yan2016(P31:2.5-37) */
        *dtheta = v3_add(dtheta_list[1], v3_cross(sum_c, dtheta_list[1]));
        /* ref Yan2016(P73:4.1-36,P76:4.1-55) */
        v3_t scul = v3_add(
            v3_cross(sum_c, dv_list[1]), v3_cross(sum_s, dtheta_list[1]));
        /* v3_t rot = vec_dot(0.5,vec_cross(dtheta_list[1],dv_list[1])); */
        /* *dv = vec_add(dv_list[1], vec_add(rot,scul)); */
        *dv = v3_add(dv_list[1], scul);
    }
    return 0;
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

/* Coarse align on the static base */
extern int align_static_base(const imu_t* imu, double lat, m3_t *Cnb)
{
    /* Fetch mean mean angular rate and specific force */
    v3_t mean_wib_b = { 0 }, mean_fib_b = { 0 };
    /* Start from the second(first epoch should not add in)*/
    for (int i = 1; i < imu->n; ++i) {
        mean_wib_b = v3_add(mean_wib_b, imu->data[i].gryo);
        mean_fib_b = v3_add(mean_fib_b, imu->data[i].accel);
    }
    double T = ins_timediff(imu->data[imu->n - 1].time, imu->data[0].time);
    mean_wib_b = v3_dot(1.0 / T, mean_wib_b);
    mean_fib_b = v3_dot(1.0 / T, mean_fib_b);

    /* Double Vector to Attitude */
    v3_t gn; gravity_ned(lat,0.0,&gn);
    v3_t wie_n = {wgs84.wie * cos(lat),0.0,- wgs84.wie * sin(lat)};
    v3_t gib_b = v3_dot(-1.0,mean_fib_b);

    dblvec2att(&gn,&wie_n,&gib_b,&mean_wib_b,Cnb);
    return 0;
}
