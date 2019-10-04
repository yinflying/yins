/**
 * @file ins.h
 * @brief ins header file
 * @author yinflying(yinflying@foxmail.com)
 * @version 0.0.1
 * @note  All abbreviation in code
 * Cnb          C_n^b, trans matrix n-axis to b-axis(the same as Qnb)
 * Cnb,Qnb,Enb  trans from n-frame to b-axis(or attitude b-axis to n-axis)
 * n,b,e,i,c,d  reference frame, navitaion/imu body/ecef/eci/carrier/odometer
 * euler        Euler attitude, [roll, pitch, yaw]
 * g            gps(or gnss phase center),gravity
 * d            data, or DELTA, mean difference, e.g. dv, dtheta
 * theta        angle increment
 * v            velocity
 * SHF          Spherical Harmonics Function
 * dbl          double
 * enu,ned      East, North, Up, Down
 * v3           3D vector
 * m3           3D matrix
 * _t           data type
 * att          attitude(could express by DCM,Quat,Euler), most express as "Enb"
 * mul          multiply
 * dcm,ctm      Direct Cosine Matrix, Coordinate Transform Matrix
 * w            Omega, Rotation rate(wie_e mean w_ie^e, project to e-axis)
 * rv           rotation vector
 * lat          latitude
 * lon          longitude
 * hgt          height
 * pos          position, most represent under geodetic coordinate, (lat,lon,hgt)
 * xyz          position under ECEF, Cartesian coordinate
 * rad          radius
 * deg          degree, unit
 * nav          navigtaion
 * cfg          config/configuration
 * kf           kalman filter
 * od           odometer
 * sol          solution
 * ft           filetype
 * kod          odometer scalar factor, true/output
 * itg          intergral
 * ycsv         yins csv file, a self-explain extended csv file format for yins
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

#ifndef INS_H
#define INS_H

/* Solution status */
#ifndef YINS_CONVERT_PARAMETER
#define DEG2RAD     0.0174532925199433      /**< deg to rad */
#define RAD2DEG     57.2957795130823        /**< rad to deg */
#define DPH2RPS     4.84813681109536e-06    /**< deg/h to rad/s */
#define RPS2DPH     206264.806247096        /**< rad/s to deg/h */
#define DPSH2RPSS   2.90888208665722e-4     /**< deg/sqrt(h) to rad/sqrt(s) */
#define RPSS2DPSH   3437.74677078493        /**< rad/sqrt(s) to deg/sqrt(h) */
#define DPHPSHZ2RPSS 4.84813681109536e-06   /**< deg/h/sqrt(Hz) to rad/sqrt(s) */
#define RPSS2DPHPSHZ 206264.806247096       /**< rad/sqrt(s) to deg/h/sqrt(Hz) */

#define G2MPS2      9.80665                 /**< g0 to m/s^2 */
#define MPS22G      0.101971621297793       /**< m/s^2 to g0 */
#define MG2MPS2     9.80665e-3              /**< millo-g0(mg) to m/s^2 */
#define MPS22MG     101.971621297793        /**< m/s^2 to milli-g0(mg) */
#define UG2MPS2     9.80665e-6              /**< micro-g0(ug) to m/s^2 */
#define MPS22UG     101971.621297793        /**< m/s^2 to micro-g0(ug) */
#define GAL2MPS2    0.01                    /**< gal to m/s^2 */
#define MPS22GAL    100.0                   /**< m/s^2 to gal */
#endif  /* YINS_CONVERT_PARAMETER */

#include <time.h>
#include <stdbool.h>
#include <stdio.h>

#ifdef RTKLIB
#include "rtklib.h"
#endif /* RTKLIB */

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief file type definition
 */
enum FT{
    FT_IMU_CSV,   /**< user defined IMU CSV file */
    FT_PVA_CSV,   /**< user defined Pos,vel,att CSV file */
    FT_OD_CSV,    /**< user defined odometer output */
    FT_IMU_YCSV,  /**< ycsv for IMU data */
    FT_PVA_YCSV,  /**< ycsv for PVA data */
    FT_OD_YCSV,   /**< ycsv for OD data  */
    FT_CFG_YCSV,  /**< ycsv for yins configuration */
    FT_NVT_ASC,   /**< Novatel ascii log file */
    FT_YGM_IMU,   /**< Yan Gongming matlab simulation package's export imu data */
    FT_YGM_AVP,   /**< Yan Gongming matlab simulation package's export trajactory data */
    FT_YGM_OD,    /**< Yan Gongming matlab simulation package's export odometer data */
    FT_NMEA       /**< nmea 1803 */
};

/**
 * @brief output reference point type
 */
enum REFPOS{
    REFPOS_IMU,     /**< imu reference center */
    REFPOS_MANUAL,  /**< manual specified pointer uner b-frame, see imup.ref_point */
    REFPOS_GPS,     /**< gps/gnss phase center, see imup.lever_arm_gps */
    REFPOS_CAR      /**< car/carrier reference center, imup.lever_arm_car */
};

/**
 * @brief Solution type(status)
 */
enum SOL{
    SOL_NONE   =     0,         /**< unvalid solution */
    SOL_MANUAL =     1,         /**< manually input solution */
    SOL_INS    =     2,         /**< INS output solution */
    SOL_DR     =     4,         /**< Dead-recoken output */
    SOL_FIXED  =     8,         /**< GNSS fixed solution */
    SOL_FLOAT  =     16,        /**< GNSS float solution */
    SOL_DGNSS  =     32,        /**< GNSS pseudo-range difference solution */
    SOL_PPP    =     64,        /**< GNSS precise point positioning solution */
    SOL_SINGLE =     128,       /**< GNSS standard point positioning solution */
    SOL_CYAW   =     256,       /**< solution with yaw constraint */
    SOL_ZUPT   =     512,       /**< solution with zero velocity update */
    SOL_ZRAU   =     1024,      /**< solution with zero rangular rate update */
    SOL_MC     =     2048,      /**< solution with motion constraint */
    SOL_EYAW   =     4096       /**< solution with extra yaw input or output yaw valid */
};

/**
 * @brief Earth parameters struct
 */
typedef struct {
    double wie; /**< rotation rate(rad s^-1) */
    double R0;  /**< Equatorial radius(m) */
    double RP;  /**< Polar radius(m) */
    double mu;  /**< gravitational constant, GM(m^3 s^-2) */
    double J2;  /**< 2nd-order gravitational Spherical Harmonics Function coefficient */
    double e;   /**< Eccentricity */
    double f;   /**< Flattening */
} earth_t;

extern earth_t wgs84;
/* Earth functions */
double earth_RN(const earth_t *eth, double lat);
double earth_RE(const earth_t *eth, double lat);

#ifndef RTKLIB
/* matrix and vector functions */
double *mat  (int n, int m);
int    *imat (int n, int m);
double *zeros(int n, int m);
double *eye  (int n);
double dot (const double *a, const double *b, int n);
double norm(const double *a, int n);
void cross3(const double *a, const double *b, double *c);
int  normv3(const double *a, double *b);
void matcpy(double *A, const double *B, int n, int m);
void matmul(const char *tr, int n, int k, int m, double alpha, const double *A, const double *B, double beta, double *C);
int  matinv(double *A, int n);
int  solve (const char *tr, const double *A, const double *Y, int n, int m, double *X);
int  lsq   (const double *A, const double *y, int n, int m, double *x, double *Q);
int  filter(double *x, double *P, const double *H, const double *v, const double *R, int n, int m);
int  smoother(const double *xf, const double *Qf, const double *xb, const double *Qb, int n, double *xs, double *Qs);
void matprint(const double *A, int n, int m, int p, int q);
int mat_symmetry(double *A, int n);

typedef struct {
    time_t time;    /**< time (s) expressed by standard time_t */
    double sec;     /**< fraction of second under 1 s */
} gtime_t;          /**< time struct */

double  timediff (gtime_t t1, gtime_t t2);
gtime_t epoch2time(const double* ep);
void    time2epoch(gtime_t t, double* ep);
gtime_t gpst2time(int week, double sec);
double  time2gpst(gtime_t t, int *week);
gtime_t timeadd(gtime_t t, double sec);
#endif /* RTKLIB */

typedef struct {
    double x, y, z;
} v3_t;     /**< 3D vector sturct */

typedef struct {
    double q0, q1, q2, q3;
} quat_t;   /**< quaternion struct */

typedef struct {
    double m11, m12, m13, m21, m22, m23, m31, m32, m33;
} m3_t;     /**< 3D matrix struct */

typedef struct {
    gtime_t time;   /**< current time */
    v3_t gyro;      /**< gyro output, angular increment */
    v3_t accel;     /**< accelermeter output, velocity increment */
} imud_t;   /**< IMU epoch data struct */

typedef struct {
    unsigned int n, nmax;
    gtime_t *time;      /**< increment */
    double  *dS;        /**< increment */
} od_t;     /**< odometer data */

typedef struct {
    unsigned int freq_imu;  /**< IMU sample rate[Hz] */
    unsigned int freq_od;   /**< Odometer sample rate[Hz] */
    gtime_t tstart;     /**< first epoch */
    v3_t gyro_noise;    /**< Gyro output noise [rad/s] */
    v3_t accel_noise;   /**< Accelermeter output notput noise [m/s^2] */
    v3_t arw;           /**< Angular random walk [rad/sqrt(s)] */
    v3_t arrw;          /**< Angular rate(gryo bias) random walk [rad/s/sqrt(s)] */
    v3_t vrw;           /**< Velocity random walk [m/s/sqrt(s)] */
    v3_t vrrw;          /**< Acceleration(accel bias) random walk[m/s^2/sqrt(s)] */
    v3_t ka;            /**< Accelermeter scalar factor(or initial value) */
    v3_t ka_std;        /**< Standard error of accelermeter scalar factor(initial value) */
    v3_t kg;            /**< Gyro scalar factor(or initial value) */
    v3_t kg_std;        /**< Standard error of Gyro scalar factor(initial value) */
    v3_t Ta;            /**< Accel bias correlation time(1st order Markov) [s] */
    v3_t Tg;            /**< Gryo bias correlation time(1st order Markov) [s] */
    v3_t inita;         /**< initial attitude [rad] */
    m3_t initQa;        /**< initial attitude uncertainty [rad^2] */
    v3_t initv;         /**< initial velocity [m/s] */
    m3_t initQv;        /**< initial velocity uncertainty [m^2/s^2] */
    v3_t initr;         /**< initial position [m] */
    m3_t initQr;        /**< initial position uncertainty */
    v3_t ba;            /**< (initial) accel bias [m/s^2] */
    v3_t ba_std;        /**< (initial) accel bias stanadard error[m/s^2] */
    v3_t bg;            /**< (initial) gryo bias [rad/s] */
    v3_t bg_std;        /**< (initial) gryo bias stanadard error[m/s^2] */
    double kod;         /**< odometer scalar factor, true/output */
    double kod_std;         /**< initial odometer scalar factor uncertainty [^2] */
    v3_t lever_arm_gps;     /**< gnss phase center position under imu frame[m] */
    v3_t lever_arm_gps_std; /**< gnss lever arm uncertainty [m] */
    v3_t lever_arm_od;      /**< odometer reference center under imu frame[m]*/
    v3_t lever_arm_od_std;  /**< odometer lever arm uncertainty [m] */
    v3_t lever_arm_car;     /**< car tailing wheel center position under imu frame[m] */
    v3_t lever_arm_car_std; /**< car tailing whell center uncertainty [m] */
    v3_t err_angle_imu;     /**< IMU install error angle(car-imu, roll, pitch, yaw)[rad] */
    v3_t err_angle_imu_std; /**< IMU install error angle uncertainty[rad] */
    v3_t err_angle_imu_rw;  /**< IMU install error angle randon walk [rad/sqrt(s)] */
    v3_t Terr_angle_imu;    /**< IMU install error angle  correlation time(1st order Markov)[s] */
    v3_t err_angle_gps;     /**< GPS install error angle(gps-imu, roll, pitch, yaw)[rad] */
    v3_t err_angle_gps_std; /**< GPS install error angle uncertainty[rad] */
    v3_t ref_point;         /**< reference point under b-frame, use for solution output [m]*/
} imup_t;   /**< IMU property struct */

typedef struct {
    unsigned int n, nmax;   /**< number of data/allocated */
    imud_t* data;           /**< IMU observation data record */
    imup_t *property;       /**< IMU property */
} imu_t;    /**< IMU full data strct */

typedef struct{
    unsigned int n;         /**< current pva number */
    unsigned int nmax;      /**< pva mamory allocated max number */
    gtime_t *time;          /**< time */
    v3_t *pos;              /**< position(ecef, xyz) */
    v3_t *vel;              /**< velocity(ecef, xyz) */
    v3_t *att;              /**< attitude(Enb, n-frame), be careful! */
    unsigned int *status;   /**< pva status(see macro SOL_*) */

    bool is_cov;            /**< if contain variance-covrainace matrix or not*/
    m3_t *Qpos;             /**< variance matrix of position, control by pva_t.is_cov */
    m3_t *Qvel;             /**< variance matirx of velocity, control by pva_t.is_cov */
    m3_t *Qatt;             /**< variance matrix of attitude, control by pva_t.is_cov */

    bool is_ext;            /**< if contain extension variables or not */
    double *yaw2;           /**< extra yaw(such as double gnss antenna), control by pva_t.is_ext */
    double *std_yaw2;       /**< extra yaw standard error, control by pva_t.is_ext */
    double *pitch2;         /**< extra pitch(such as double gnss antenna), control by pva_t.is_ext */
    double *std_pitch2;     /**< extra pitch stanadard error, control by pva_t.is_ext */
    unsigned int *ext_status;   /**< extra variebles status, control by pva_t.is_ext */
} pva_t;

typedef struct{
    bool isx_ba;        /**< if acceleremter bias in state vector or not */
    bool isx_bg;        /**< if gryo bias in state vector or not */
    bool isx_kax;       /**< if accel scalar factor(x-axis) in state vector or not */
    bool isx_kay;       /**< if accel scalar factor(y-axis) in state vector or not */
    bool isx_kaz;       /**< if accel scalar factor(z-axis) in state vector or not */
    bool isx_kgx;       /**< if gyro scalar factor(x-axis) in state vector or not */
    bool isx_kgy;       /**< if gyro scalar factor(y-axis) in state vector or not */
    bool isx_kgz;       /**< if gyro scalar factor(z-axis) in state vector or not */
    bool isx_kod;       /**< if scalar factor of odometer in state vector or not */
    bool isx_eroll;     /**< if install roll error angel in state vector or not */
    bool isx_epitch;    /**< if install pitch error angle in state vector or not */
    bool isx_eyaw;      /**< if install yaw error angel in state vector or not */
    bool isx_armgps;    /**< if lever arm of imu to GPS in state vector or not */
    bool isx_armcar;    /**< if lever arm of imu to car in state vector or not */
    bool is_odincre;    /**< odometer/ins distance increment combine mode */
    bool iskf_itgdS_sagehusa;   /**< if using sage-husa self-adaption filter for itgdS */
    unsigned char max_ny;   /**< max observations when kalman filter measurement update */
    unsigned short nZST;    /**< zero speed test window size */

    bool issol_header;      /**< if output solution header or not(ycsv format) */
    enum REFPOS sol_refpos; /**< output solution reference point type */
    float feedratio;        /**< kalman filter feedback fatio */

    unsigned char IPOS;     /**< index position in state vector */
    unsigned char IVEL;     /**< index velocity in state vector */
    unsigned char IATT;     /**< index attitude in state vector */
    unsigned char IBA;      /**< index acceleremeter bias in state vector */
    unsigned char IBG;      /**< index gyro bias in state vector */
    unsigned char IKAx;     /**< index accel scalar factor(x-axis) in state vector */
    unsigned char IKAy;     /**< index accel scalar factor(y-axis) in state vector */
    unsigned char IKAz;     /**< index accel scalar factor(z-axis) in state vector */
    unsigned char IKGx;     /**< index gyro scalar factor(x-axis) in state vector */
    unsigned char IKGy;     /**< index gyro scalar factor(y-axis) in state vector */
    unsigned char IKGz;     /**< index gyro scalar factor(z-axis) in state vector */
    unsigned char IKOD;     /**< index scalar factor of odometer in state vector */
    unsigned char IEROLL;   /**< index imu to car install roll error in state vector */
    unsigned char IEPITCH;  /**< index imu to car install pitch error in state vector */
    unsigned char IEYAW;    /**< index imu to car install yaw error in state vector */
    unsigned char IARMGPS;  /**< index lever arm of imu to GPS in state vector */
    unsigned char IARMCAR;  /**< index lever arm of imu to car in state vector */
    unsigned char nx;       /**< numbner of states */

    unsigned char ITGVEL_OD;    /**< index odometer velocity integral in integral state vector */
    unsigned char ITGVEL_GPS;   /**< index gps velocity integral in integral state vector */
    unsigned char ITGVEL_INS;   /**< index ins velocity integral in integral state vector */
    unsigned char ITGF_INS;     /**< index ins specific force integral in integral state vector */
    unsigned char ITGC_EYAW;    /**< index yaw error angle coefficient in integral state vector */
    unsigned char ITGC_EPITCH;  /**< index pitch error angle coefficient in integral state vector */
    unsigned char ITGC_EROLL;   /**< index roll error angle coefficient in integral state vector */
    unsigned char ITGC_KOD;     /**< index odomemter scalar factor coefficient in intergral state vector */
    unsigned char nitg;         /**< numbner of integral states */

    unsigned char IR_itgdS;     /**< index intergral odometer distance increment measurement noise in R vector */
    unsigned char nR;           /**< necessary measurent noise number */

    bool is_imu_samenoise;      /**< if imu three axis have same noise(use for speed up kalman filter) */
} cfg_t;            /**< configure struct */

extern cfg_t cfg;   /**< global configuration */
int updatecfg(void);

typedef struct{
    gtime_t time;       /**< current solution time */
    unsigned int status;/**< solution status, see macro SOL_* */
    m3_t dcm;           /**< attitude in DCM */
    quat_t quat;        /**< attitude in quaternion */
    m3_t Qatt;          /**< var-covariance matrix of attitude */
    v3_t vel;           /**< velocity */
    m3_t Qvel;          /**< var-covariance matrix of velocity */
    v3_t pos;           /**< position */
    m3_t Qpos;          /**< var-covariance matrix of postion  */
    v3_t ba;            /**< accelermeter bias */
    v3_t ba_std;        /**< standard error of accelermeter bais */
    v3_t bg;            /**< gryo bias */
    v3_t bg_std;        /**< standard error of gyro bias */
    v3_t ka;            /**< accelermeter scalar factor */
    v3_t ka_std;        /**< standard error of accelermeter scalar factor */
    v3_t kg;            /**< gyro scalar scalar factor */
    v3_t kg_std;        /**< standard error of gyro scalar factor */
    double kod;         /**< scalar factor of odometer, true/output */
    double std_kod;     /**< standard error of odometer scalar factor */
    m3_t Cbc;           /**< install error angle */
    v3_t std_Cbc;       /**< standard error of install error angle */
} solins_t;             /**< ins solution struct */

typedef struct{
    gtime_t time;       /**< current time */
    double idt;         /**< time interval of imu */
    double odt;         /**< time interval of od */
    unsigned char nx;   /**< length of full state of x */
    unsigned char ny;   /**< length of full state of y */
    double *x;          /**< state vector */
    double *P;          /**< var-covariance matrix */
    double *F;          /**< transition matrix */
    double *Q;          /**< System noise covariance matrix */
    double *H;          /**< measurement matrix(transpose) */
    solins_t  *sol;     /**< solution of kalman fileter */
    gtime_t itg_start;  /**< intergral start time */
    double *itg;        /**< integral variables */
    double *R;          /**< necessary measurement noise(not for normal) */
    imud_t *imud;               /**< imu data list */
    unsigned short nimud;       /**< number of imud_t struct in kf_t.imud */
    unsigned short imudend;     /**< last imu data  */
    unsigned int ZST_count;     /**< zero speed test count */
} kf_t;     /**< kalman filter status struct */

/* Attitude transformation */
int rv2quat(const v3_t* dtheta, quat_t* quat);
int rv2dcm(const v3_t *dtheta, m3_t * dcm);
int dcm2rv(const m3_t *dcm, v3_t *dtheta);
int quat2rv(const quat_t *dcm, v3_t *dtheta);
int euler2quat(const v3_t* euler, quat_t* quat);
int quat2euler(const quat_t* quat, v3_t* euler);
int dcm2quat(const m3_t* dcm, quat_t* quat);
int quat2dcm(const quat_t* quat, m3_t* dcm);
int dcm2euler(const m3_t* dcm, v3_t* euler);
int euler2dcm(const v3_t* euler, m3_t* dcm);
int att2dcm(const v3_t *att, m3_t *dcm);
int dcm2att(const m3_t *dcm, v3_t *att);
int att2quat(const v3_t *att, quat_t *quat);
int quat2att(const quat_t *quat, v3_t *att);

/* 3D vector operation */
int asymmetric_mat(const v3_t* v3, m3_t* mat);
m3_t v3_askew(v3_t v3);
v3_t v3_cross(v3_t v1, v3_t v2);
v3_t v3_add(v3_t v1, v3_t v2);
v3_t v3_del(v3_t v1, v3_t v2);
v3_t v3_scalar(double s, v3_t v);
double v3_norm(v3_t v3);
int v3_normalize(v3_t *v3);
double v3_mul_rxc(v3_t v1, v3_t v2); /* row vector x column vector */
m3_t v3_mul_cxr(v3_t v1, v3_t v2);   /* column vector x row vector */
v3_t v3_mul_cxc(v3_t v1, v3_t v2);   /* multiply corresponding elements */
m3_t v3_diag(v3_t diag);
v3_t v3_pow(v3_t v, double order);
v3_t v3_sqrt(v3_t v);
bool v3_equal(const v3_t *v1, const v3_t *v2, double eps);
v3_t v3_sum(const v3_t *v3_list, int n);
v3_t v3_mean(const v3_t *v3_list, int n);
v3_t v3_rms(const v3_t *v3_list, int n);
v3_t v3_std(const v3_t *v3_list, int n);

/* 3D matrix operator */
v3_t m3_iaskew(m3_t m3);
m3_t m3_T(m3_t A);
int m3_inv(m3_t *A);
double m3_det(const m3_t *A);
m3_t m3_add(m3_t A, m3_t B);
m3_t m3_del(m3_t A, m3_t B);
m3_t m3_scalar(double alpha, m3_t A);
m3_t m3_mul(m3_t A, m3_t B);
v3_t m3_mul_v3(m3_t A, v3_t B);
v3_t m3_diag(m3_t diag);
m3_t m3_pow(m3_t A, double order);
bool m3_equal(const m3_t *A, const m3_t *B, double eps);
void m3_swap_row(m3_t *A, int r1, int r2);
void m3_swap_clm(m3_t *A, int c1, int c2);
int m3_SVD(const m3_t *A, m3_t *U, v3_t *D, m3_t *V);
int m3_LU(const m3_t *A, m3_t *L, m3_t *U, m3_t *P);

/* quaternion operation */
int quat_normalize(quat_t *quat);
int quat_conj(quat_t *quat);
int quat_inv(quat_t *quat);
double quat_dot(quat_t P, quat_t Q);
v3_t quat_cross(quat_t P, quat_t Q);
double quat_norm(quat_t P);
quat_t quat_mul(quat_t P, quat_t Q);
v3_t quat_mul_v3(quat_t quat, v3_t vec);
bool quat_equal(const quat_t *P, const quat_t *Q, double eps);

/* Attitude add and del operation */
v3_t euler_add(v3_t Eab, v3_t Ebc);
v3_t euler_del(v3_t Eab, v3_t Eac);
double yaw_del(double yaw1, double yaw2);
int euler_addphi(v3_t *Eab, const v3_t *phi_bc);
int euler_delphi(v3_t *Eab, const v3_t *phi_bc);

/* coordinate transformation */
m3_t formCen_ned(double lat, double lon);
int ned2ecef(v3_t* pos, v3_t* vel, m3_t* Cbn);
int ecef2ned(v3_t* xyz, v3_t* vel, m3_t* Cbe);
int ned2ecefQ(const v3_t *pos, m3_t *Qpos, m3_t *Qvel, m3_t *Qatt);
int ecef2nedQ(const v3_t *xyz, m3_t *Qxyz, m3_t *Qvel, m3_t *Qatt);
m3_t att2Cbe(const v3_t *pos, const v3_t *Enb);
quat_t att2Qbe(const v3_t *pos, const v3_t *Enb);
v3_t att2Ebe(const v3_t *pos, const v3_t *Enb);
v3_t Cbe2att(const v3_t *xyz, const m3_t *Cbe);
v3_t Qbe2att(const v3_t *xyz, const quat_t *Qbe);
v3_t Ebe2att(const v3_t *xyz, const v3_t *Ebe);
bool is_blh(const v3_t *pos);

/* Gravity model*/
int gravity_ecef(const v3_t *r, v3_t* ge);
int gravity_ned(double lat, double hgt, v3_t* gn);

/* INS Align */
int align_coarse_static_base(const imu_t *imu, double lat, m3_t *Cnb);
int dblvec2att(const v3_t *vn1, const v3_t *vn2, const v3_t *vb1, const v3_t*vb2, m3_t *Cnb);
int align_coarse_inertial(const imu_t *imu, double lat, m3_t *Cnb);
int align_coarse_wuhba(const imu_t *imu, double lat, const v3_t *veb_n,
                       unsigned int Nveb_n, m3_t *Cnb);

/* TODO: Static and dynamic judgment */

/* INS navigation */
int ins_nav_ecef(double dt, const v3_t* dtheta, const v3_t* dv, v3_t *r, v3_t *v, quat_t *q);
/* TODO NAV EQUATION BACKWARD */
int ins_nav_ecef_back(double dt,const v3_t *dtheta, const v3_t *dv, v3_t *r, v3_t *v, quat_t *q);
/* TODO NAV EQUATION NED */
int ins_nav_ned(double dt, const v3_t *dtheta, const v3_t *dv, v3_t *pos, v3_t *v, quat_t *q);
int multisample(const v3_t *dtheta_list, const v3_t *dv_list, int N, v3_t *dtheta, v3_t *dv);
int dr_nav_ecef(double dt, const v3_t* dtheta, double dS, v3_t *r, quat_t *q);

/* INS IO operation */
int yins_readf(const char * fname, enum FT ft, imu_t *imu, pva_t *pva, od_t *od);
int yins_writef(const char *fname, enum FT ft, const imu_t *imu, const pva_t *pva, const od_t *od);
int imu_init(imu_t *imu);
void imu_free(imu_t *imu);
int imu_add(imu_t *imu, const imud_t *data);
void od_init(od_t *od);
void od_add(od_t *od, const gtime_t *time, const double *dS);
void od_free(od_t *od);
void pva_init(pva_t *pva);
void pva_resize(pva_t *pva, unsigned int n);
void pva_free(pva_t *pva);
int outsolins(FILE *fp, const solins_t *sol, const imup_t *imup);

/* jacobi matrix */
int jacobi_trans_Ebe2Ebe(m3_t *F, double dt);
int jacobi_trans_Ebe2bg(m3_t *F, double dt, const m3_t *Cbe);
int jacobi_trans_veb_e2Ebe(m3_t *F, const m3_t *Cbe, const v3_t *dtheta);
int jacobi_trans_veb_e2veb_e(m3_t *F, double dt);
int jacobi_trans_veb_e2reb_e(m3_t *F, double dt, const v3_t *reb_e);
int jacobi_trans_veb_e2ba(m3_t *F, double dt, const m3_t *Cbe);
int jacobi_trans_reb_e2veb_e(m3_t *F, double dt);
int jacobi_trans_markov(m3_t *F, double dt, const v3_t *T);
int jacobi_meas_reb_e2Ebe(m3_t *H, const m3_t *Cbe, const v3_t *lever_arm);
int jacobi_meas_veb_e2Ebe(m3_t *H, const m3_t *Cbe, const v3_t *wib_b, const v3_t *lever_arm);
int jacobi_meas_veb_e2bg(m3_t *H, const m3_t *Cbe, const v3_t *lever_arm);
int jacobi_meas_Ebn2Ebe(m3_t *H, const m3_t *Cen);

int jacobi_trans_Ebe2kgx(double *F, const m3_t *Cbe, const v3_t *dtheta);
int jacobi_trans_Ebe2kgy(double *F, const m3_t *Cbe, const v3_t *dtheta);
int jacobi_trans_Ebe2kgz(double *F, const m3_t *Cbe, const v3_t *dtheta);

int jacobi_trans_veb_e2kax(double *F, const m3_t *Cbe, const v3_t *dv);
int jacobi_trans_veb_e2kay(double *F, const m3_t *Cbe, const v3_t *dv);
int jacobi_trans_veb_e2kaz(double *F, const m3_t *Cbe, const v3_t *dv);

/* determinate attitude by velocity */
int vel2yaw(const v3_t *veb_n, double *yaw, const m3_t *Qveb_n, double *Qyaw);
double dyaw_swerve(double veb_b_x, double web_n_z, double dL);
double dyaw_swerve_Q(double dyaw, double veb_b_x, double web_n_z, double dL, double Qv, double Qw, double QL);
int vel2pitch(const v3_t *veb_n, double *pitch, const m3_t *Qveb_n, double *Qpitch); /* TODO: vel2pithc() */

/* kalman filter operation for ins */
int inskf_init(kf_t *inskf, const imup_t *imup);
int inskf_udstate(kf_t *inskf, const imud_t *imud, const imup_t *imup);
int inskf_udmeasr(kf_t *inskf, const v3_t *reg_e, const m3_t *Qreg_e, const v3_t *lever_arm_gps);
int inskf_udmeasv(kf_t *inskf, const v3_t *veg_e, const m3_t *Qveg_e, const v3_t *lever_arm_gps, const v3_t *wib_b);
int inskf_udmeas_vod(kf_t *inskf, double vod, const v3_t *vod_std);
int inskf_constraint_yaw(kf_t *inskf, const v3_t *wib_b, const imup_t *imup);
int inskf_MC_car(kf_t *inskf, const v3_t * wib_b, const v3_t *lever_arm_car, double std_vy, double std_vz);
int inskf_ZST_pos(kf_t *inskf, const solins_t *last_sol);
bool inskf_ZST_imu(kf_t *inskf, double mean_dv, double std_dv);
bool inskf_ZST_auto(kf_t *inskf, const solins_t *last_sol, double *meandv, double *std_dv);
int inskf_ZRAU(kf_t *inskf, const quat_t *last_qbe, double dt, const m3_t *Qweb_b);
int inskf_ZUPT(kf_t *inskf, const m3_t *Qveb_e);
int inskf_ZAU(kf_t *inskf, double yaw, double Qyaw);
int inskf_udmeas_Cbe(kf_t *inskf, const m3_t *Cbe, const m3_t *QEbe);
int inskf_udmeas_yaw(kf_t *inskf, double yaw, double Qyaw);
int inskf_uditg_dS(kf_t *inskf, double dS, const v3_t *lever_arm_car, const v3_t *wib_b);
int inskf_udmeas_itgdS(kf_t *inskf, const v3_t *QitgdS);
int inskf_feedback(kf_t *inskf, unsigned int SOL_TYPE);
int inskf_norm_innov(kf_t *inskf, const double *R, const double *dz, double *ndz);
int inskf_close(kf_t *inskf);

/* logging */
enum LOG_LEVEL{
    LEVEL_OFF, LEVEL_FATAL, LEVEL_ERROR, LEVEL_WARN, LEVEL_INFO, LEVEL_DEBUG, LEVEL_TRACE
};

extern enum LOG_LEVEL FILE_LOG_LEVEL;       /**< log to file minimal level*/
extern enum LOG_LEVEL STDOUT_LOG_LEVEL;     /**< log to stdout minimal level */
extern FILE *LOG_FP;                        /**< log file pointer */
extern gtime_t LOG_CURTIME;                 /**< cureent log time */

int LOG_OPEN(const char *file);
int LOG_CLOSE(void);
void LOG_FATAL(const char *format, ...);
void LOG_ERROR(const char *format, ...);
void LOG_WARN(const char *format, ...);
void LOG_INFO(const char *format, ...);
void LOG_DEBUG(const char *format, ...);
void LOG_TRACE(const char *format, ...);
void SHOW_PROGRESS(const char *msg, int current, int total);

/* Constant defination */
extern const m3_t I3;       /**< Unit 3D matrix */
extern const m3_t O3;       /**< Zero 3D matirx */
extern const v3_t V0;       /**< Zero 3D vector */
extern const v3_t V1;       /**< Unit 3D vector */
extern v3_t WIE_E;          /**< Earth rotation rate under ECEF */

/* utilities */
double yins_randn(void);
v3_t v3_randn(double mean, double sigma);
void v3_paste(double *dest, const v3_t *v3);
void m3_paste(double *dest, int n, const m3_t *m3);
void v3_copy(v3_t *v3_dest, const double *src);
void m3_copy(m3_t *m3_dest, const double *src, int n);
void imu_orientation_adjust(imud_t *imud, const char *or1, const char *or2);
unsigned int soltype_add(unsigned int status, unsigned int SOL_TYPE);
unsigned int soltype_remove(unsigned int status, unsigned int SOL_TYPE);
bool is_soltype(unsigned int status, unsigned int SOL_TYPE);

double markov_std2rw(double std, double T);
double markov_rw2std(double rw, double T);

#ifdef __cplusplus
}
#endif

#endif /* ifndef INS_H */
