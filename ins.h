/* NOTE:
 * Cnb => C_n^b => trans matrix n-axis to b-axis(the same as qnb)
 * g=>gps,gravity
 * d=>data
 * SHF=>Spherical Harmonics Function
 * dbl=>double
 * enu,ned=>East, North, Up, Down
 * v3=> 3D vector
 * m3=> 3D matrix
 * _t => data type
 * att=> attitude(could express by DCM,Quat,Euler)
 * mul => multiply
 * */
#ifndef INS_H
#define INS_H

#define FT_CSV 0
#define FT_NVT 1
#define FT_RNX 2

#include <time.h>

typedef struct {
    double wie; /* rotation rate(rad s^-1) */
    double R0;  /* Equatorial radius(m) */
    double RP;  /* Polar radius(m) */
    double mu;  /* gravitational constant, GM(m^3 s^-2) */
    double J2;  /* 2nd-order gravitational Spherical Harmonics Function coefficient */
    double e;   /* Eccentricity */
    double f;   /* Flattening */
} earth_t;
earth_t wgs84;

#ifndef GTIME_T
#define GTIME_T
typedef struct { /* time struct */
    time_t time; /* time (s) expressed by standard time_t */
    double sec;  /* fraction of second under 1 s */
} gtime_t;
#endif /* ifndef GTIME_T */
double  ins_timediff (gtime_t t1, gtime_t t2);
gtime_t ins_epoch2time(const double* ep);
void ins_time2epoch(gtime_t t, double* ep);
gtime_t ins_gpst2time(int week, double sec);

typedef struct {
    double i, j, k;
} v3_t;

typedef struct {
    gtime_t time;
    v3_t gryo;
    v3_t accel;
} imud_t;

typedef struct {
    v3_t arw;    /**< Angular random walk */
    v3_t arrw;   /**< Augluar rate(gryo bias) random walk  */
    v3_t vrw;    /**< Velocity radom walk */
    v3_t vrrw;   /**< Acceleration(accel bias) radom walk */
    v3_t Ta;     /**< Accel bias correlation time */
    v3_t Tg;     /**< Gryo bias correlation time */
    v3_t initr;  /**< initial position */
    v3_t initQr; /**< initial position uncertinty */
    v3_t initv;  /**< initial velocity */
    v3_t initQv; /**< initial velocity uncertainty */
    v3_t inita;  /**< initial attitude */
    v3_t initQa; /**< initial attitude uncertainty */
    v3_t initQgb; /**< initial gryo bias uncertainty */
    v3_t initQab; /**< initial accel bias uncertainty*/
    v3_t lever_arm; /**< Reference position under imu frame */

    int n, nmax;   /**< number of data/allocated */
    imud_t* data;  /**< IMU observation data record */
} imu_t;

typedef struct {
    double q0, q1, q2, q3;
} quat_t;

typedef struct {
    double m11, m12, m13, m21, m22, m23, m31, m32, m33;
} m3_t;

/* Attitude transformation */
int dtheta2quat(const v3_t* dtheta, quat_t* quat);
int euler2quat(const v3_t* euler, quat_t* quat);
int quat2euler(const quat_t* quat, v3_t* euler);
int dcm2quat(const m3_t* dcm, quat_t* quat);
int quat2dcm(const quat_t* quat, m3_t* dcm);
int dcm2euler(const m3_t* dcm, v3_t* euler);
int euler2dcm(const v3_t* euler, m3_t* dcm);

/* 3D vector operation */
int asymmetric_mat(const v3_t* v3, m3_t* mat);
v3_t v3_cross(v3_t v1, v3_t v2);
v3_t v3_add(v3_t v1, v3_t v2);
v3_t v3_del(v3_t v1, v3_t v2);
v3_t v3_dot(double s, v3_t v);
double v3_norm(v3_t v3);
double v3_mul_rxc(v3_t v1, v3_t v2); /* row vector x column vector */
m3_t v3_mul_cxr(v3_t v1, v3_t v2);  /* column vector x row vector */
m3_t v3_diag(v3_t diag);

/* 3D matrix operator */
m3_t m3_transpose(m3_t A);
m3_t m3_mul(m3_t A, m3_t B);
v3_t m3_mul_v3(m3_t A, v3_t B);
v3_t m3_diag(m3_t diag);

/* quaternion operation */
int quat_normalize(quat_t* quat);
int quat_inv(quat_t* quat);
quat_t quat_mul(quat_t P, quat_t Q);
v3_t quat_mul_v3(quat_t quat, v3_t vec);

/* coordinate transformantion */
m3_t formCen_ned(double lat, double lon);
int ned2ecef(v3_t* pos, v3_t* vel, m3_t* att);
int ecef2ned(v3_t* pos, v3_t* vel, m3_t* att);

/* Gravity model*/
int gravity_ecef(const v3_t *r, v3_t* ge);
int gravity_ned(double lat, double hgt, v3_t* gn);

/* INS Align */
int align_static_base(const imu_t *imu, double lat, m3_t *Cnb);

int dblvec2att(const v3_t *vn1, const v3_t *vn2, const v3_t *vb1,
        const v3_t*vb2, m3_t *Cnb);

/* INS navgataion */
int nav_equations_ecef(double dt, const v3_t* dtheta, const v3_t* dv,
    v3_t* r, v3_t* v, quat_t* q);
int multisample(const v3_t* dtheta_list, const v3_t* dv_list, int N,
    v3_t* dtheta, v3_t* dv);

/* INS IO operation */
int readimu_file(const char* infile, imu_t* imu, int FILETYPE);
int imu_trans_rnx(const imu_t* imu, const char* outfile);
void freeimu(imu_t* imu);
int addimudata(imu_t* imu, const imud_t* data);

/* TODO INS Align */
/* TODO NAV EQUATION BACKWARD */
/* Static and dynamic judgement */

#endif /* ifndef INS_H */
