#ifndef INS_H
#define INS_H

#define FT_CSV 0
#define FT_NVT 1
#define FT_RNX 2

#include <time.h>

typedef struct {
    double wie;
    double R0;
    double mu;
    double J2;
    double e;
} earth_t;
earth_t wgs84;

#ifndef GTIME_T
#define GTIME_T
typedef struct { /* time struct */
    time_t time; /* time (s) expressed by standard time_t */
    double sec;  /* fraction of second under 1 s */
} gtime_t;
#endif /* ifndef GTIME_T */

typedef struct {
    double i, j, k;
} vec3_t;

typedef struct {
    gtime_t time;
    vec3_t gryo;
    vec3_t accel;
} imud_t;

typedef struct {
    vec3_t arw;    /**< Angular random walk */
    vec3_t arrw;   /**< Augluar rate(gryo bias) random walk  */
    vec3_t vrw;    /**< Velocity radom walk */
    vec3_t vrrw;   /**< Acceleration(accel bias) radom walk */
    vec3_t Ta;     /**< Accel bias correlation time */
    vec3_t Tg;     /**< Gryo bias correlation time */
    vec3_t initr;  /**< initial position */
    vec3_t initQr; /**< initial position uncertinty */
    vec3_t initv;  /**< initial velocity */
    vec3_t initQv; /**< initial velocity uncertainty */
    vec3_t inita;  /**< initial attitude */
    vec3_t initQa; /**< initial attitude uncertainty */
    vec3_t initQgb; /**< initial gryo bias uncertainty */
    vec3_t initQab; /**< initial accel bias uncertainty*/

    int n, nmax;   /**< number of data/allocated */
    imud_t* data;  /**< IMU observation data record */
} imu_t;

typedef struct {
    double q0, q1, q2, q3;
} quat_t;

typedef struct {
    double m11, m12, m13, m21, m22, m23, m31, m32, m33;
} dcm_t;

/* Attitude transformation */
int dtheta2quat(const vec3_t* dtheta, quat_t* quat);
int euler2quat(const vec3_t* euler, quat_t* quat);
int quat2euler(const quat_t* quat, vec3_t* euler);
int dcm2quat(const dcm_t* dcm, quat_t* quat);
int quat2dcm(const quat_t* quat, dcm_t* dcm);
int dcm2euler(const dcm_t* dcm, vec3_t* euler);
int euler2dcm(const vec3_t* euler, dcm_t* dcm);

/* 3D vector operation */
int asymmetric_mat(const vec3_t* v3, dcm_t* mat);
vec3_t vec_cross(vec3_t v1, vec3_t v2);
vec3_t vec_add(vec3_t v1, vec3_t v2);
vec3_t vec_del(vec3_t v1, vec3_t v2);
vec3_t vec_dot(double s, vec3_t v);

/* 3D matrix operator */
dcm_t m3_transpose(dcm_t A);
dcm_t m3_mul(dcm_t A, dcm_t B);
vec3_t m3_mul_vec(dcm_t A, vec3_t B);

/* quaternion operation */
int quat_normalize(quat_t* quat);
int quat_inv(quat_t* quat);
quat_t quat_mul(quat_t P, quat_t Q);
vec3_t quat_mul_vec(quat_t quat, vec3_t vec);

/* coordinate transformantion */
dcm_t formCen_ned(double lat, double lon);
int ned2ecef(vec3_t* pos, vec3_t* vel, dcm_t* att);
int ecef2ned(vec3_t* pos, vec3_t* vel, dcm_t* att);
int gravity_ecef(const vec3_t* r, vec3_t* ge);

/* INS navgataion */
int nav_equations_ecef(double dt, const vec3_t* dtheta, const vec3_t* dv,
    vec3_t* r, vec3_t* v, quat_t* q);
int multisample(const vec3_t* dtheta_list, const vec3_t* dv_list, int N,
    vec3_t* dtheta, vec3_t* dv);

/* INS IO operation */
int readimu_file(const char* infile, imu_t* imu, int FILETYPE);
int imu_trans_rnx(const imu_t* imu, const char* outfile);
void freeimu(imu_t* imu);
int addimudata(imu_t* imu, const imud_t* data);

#endif /* ifndef INS_H */
