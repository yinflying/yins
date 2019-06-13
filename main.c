#include "ins.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define DEG2RAD 0.0174532925199433
#define RAD2DEG 57.2957795130823
#define DPH2RPS 4.84813681109536e-06  /* deg/h to rad/s */
#define MG2MPS2 9.78046e-3             /* mg to m/s^2 */
#define DPSH2RPSS 2.90888208665722e-4 /* deg/sqrt(h) to rad/sqrt(s) */

void print_imu(const imu_t* imu)
{
    for (int i = 0; i < imu->n; ++i) {
        fprintf(stdout, "%16.10f,", imu->data[i].time.sec);
        fprintf(stdout, "%16.10f,%16.10f,%16.10f,", imu->data[i].accel.i,
            imu->data[i].accel.j, imu->data[i].accel.k);
        fprintf(stdout, "%16.10f,%16.10f,%16.10f\n", imu->data[i].gryo.i,
            imu->data[i].gryo.j, imu->data[i].gryo.k);
    }
}

static double str2num(const char* s, int i, int n)
{
    double value;
    char str[256], *p = str;

    if (i < 0 || (int)strlen(s) < i || (int)sizeof(str) - 1 < n)
        return 0.0;
    for (s += i; *s && --n >= 0; s++)
        *p++ = *s == 'd' || *s == 'D' ? 'E' : *s;
    *p = '\0';
    return sscanf(str, "%lf", &value) == 1 ? value : 0.0;
}

void nav_ins_ecef(const imu_t* imu, int N, v3_t r, v3_t v, quat_t q)
{
    v3_t euler;
    v3_t dtheta_list[5], dv_list[5], dtheta, dv;
    double dt = imu->data[1].time.sec - imu->data[0].time.sec;

    m3_t dcm; v3_t pos,veb_n,Enb;

    if (N < 0) /* update first abs(N)-1 when N<0*/
        for (int i = 0; i < abs(N) - 1; ++i) {
            nav_equations_ecef(
                dt, &imu->data[i].gryo, &imu->data[i].accel, &r, &v, &q);
            fprintf(stdout, "%6.3f ", imu->data[i].time.sec);
            fprintf(stdout, "%16.10f %16.10f %16.10f ", r.i, r.j, r.k);
            fprintf(stdout, "%16.10f %16.10f %16.10f ", v.i, v.j, v.k);
            quat2euler(&q, &euler);
            fprintf(
                stdout, "%16.10f %16.10f %16.10f\n", euler.i, euler.j, euler.k);
        }
    int nSample = N < 0 ? 1 : N;
    for (int i = abs(N) - 1; i < imu->n; i += nSample) {
        int start_ind = i - abs(N) + 1;
        for (int j = 0; j < abs(N); j++) {
            dtheta_list[j] = imu->data[start_ind + j].gryo;
            dv_list[j] = imu->data[start_ind + j].accel;
        }
        multisample(dtheta_list, dv_list, N, &dtheta, &dv);
        nav_equations_ecef(dt * nSample, &dtheta, &dv, &r, &v, &q);
        quat2dcm(&q,&dcm); pos=r; veb_n = v;
        ecef2ned(&pos,&veb_n,&dcm); dcm = m3_transpose(dcm); dcm2euler(&dcm,&Enb);
        fprintf(stdout, "%6.3f ", imu->data[i].time.sec);
        fprintf(stdout, "%16.10f %16.10f %16.10f ", pos.i*RAD2DEG, pos.j*RAD2DEG, pos.k);
        fprintf(stdout, "%16.10f %16.10f %16.10f ", veb_n.i, veb_n.j, veb_n.k);
        quat2euler(&q, &euler);
        fprintf(stdout, "%16.10f %16.10f %16.10f\n", Enb.i*RAD2DEG, Enb.j*RAD2DEG, Enb.k*RAD2DEG);
    }
}

void get_init_para(char* infile, v3_t* r, v3_t* v, quat_t* q)
{
    FILE* fp;
    if (!(fp = fopen(infile, "r"))) {
        fprintf(stderr, "trajectory file open error\n");
        return;
    }
    char buff[240];
    fgets(buff, 240, fp);
    v3_t e;
    r->i = str2num(buff, 17, 16);
    r->j = str2num(buff, 34, 16);
    r->k = str2num(buff, 51, 16);
    v->i = str2num(buff, 68, 16);
    v->j = str2num(buff, 85, 16);
    v->k = str2num(buff, 102, 16);
    e.i = str2num(buff, 119, 16);
    e.j = str2num(buff, 136, 16);
    e.k = str2num(buff, 153, 16);
    euler2quat(&e, q);
    fclose(fp);
    return;
}

void set_imu(imu_t* imus)
{
    imus->initQr = (v3_t) { 10.0, 10.0, 10.0 };
    imus->initQv = (v3_t) { 1.0, 1.0, 1.0 };
    imus->initQa = (v3_t) { 2.0 * DEG2RAD, 2.0 * DEG2RAD, 2.0 * DEG2RAD };
    imus->initQab = (v3_t) { 10.0 * MG2MPS2, 10.0 * MG2MPS2, 10.0 * MG2MPS2 };
    imus->initQgb = (v3_t) { 10.0 * DPH2RPS, 10.0 * DPH2RPS, 10.0 * DPH2RPS };
    imus->arw = (v3_t) { 0.0667 * DPSH2RPSS, 0.0667 * DPSH2RPSS,
        0.0667 * DPSH2RPSS };
    imus->arrw = (v3_t) { 1e-6, 1e-6, 1e-6 };
    imus->vrw = (v3_t) { 200 * 9.8e-6, 200 * 9.8e-6, 200 * 9.8e-6 };
    imus->vrrw = (v3_t) { 1e-4, 1e-4, 1e-4 };
    imus->Ta = (v3_t) { 1000.0, 1000.0, 1000.0 };
    imus->Tg = (v3_t) { 1000.0, 1000.0, 1000.0 };
    imus->lever_arm = (v3_t) {0.0, -0.18, -0.05};
}

void test1()
{
    imu_t imu;
    yins_readimu("./data/testIMUInter.csv",&imu,FT_CSV);
    v3_t r,v;
    v3_t att = {-0.5257575320*DEG2RAD, -2.1960440963*DEG2RAD, -123.503797619*DEG2RAD};
    v = (v3_t){ 0.0160980, -0.0429563, -0.0189353 };
    r = (v3_t){ (31+34/60.0+47.2719024/3600.0)*DEG2RAD,
        (104 +27/60.0+ 27.1144676/3600)*DEG2RAD, 569.916321};
    m3_t dcm; euler2dcm(&att,&dcm); dcm = m3_transpose(dcm); // Cnb => Cbn
    ned2ecef(&r,&v,&dcm); //blh=>reb_e, veb_n=>veb_e, Cbn => Cbe
    quat_t q; dcm2quat(&dcm,&q);
    v = (v3_t){ 0.0160980, -0.0429563, -0.0189353 };
    nav_ins_ecef(&imu,-2,r,v,q);

    freeimu(&imu);
}

int main()
{
    //imu_t imu;
    // yins_readimu("../testdata/20180416_rover_345_multi.ASC", &imu, FT_NVT);
    // set_imu(&imu);
    // yins_imu2rnx(&imu, "../testdata/1_level_rnx3/IMU.rnx");

    //v3_t r;
    //v3_t v;
    //quat_t q;
    // yins_readimu("./data/ECEF_IMU_meas_1.csv",&imu,FT_CSV);
    // print_imu(&imu);
    // get_init_para("./data/ECEF_trajectory_1.csv",&r,&v,&q);
    // nav_ins_ecef(&imu,-2, r, v, q);
    //

    //freeimu(&imu);

    test1();
    return 0;
}
