/**
 * @file inskf.c
 * @brief ins kalman filter related function
 * @author yinflying(yinflying@foxmail.com)
 * @note
 *  2019-05-21 Created
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

/* Simplfiy the filter usage */
#define KF_FILTER(kf, dz, Qdz)                                               \
    filter(kf->x, kf->P, kf->H, (dz), (Qdz), kf->nx, kf->ny)

/* set kf->ny and set kf->H to zero */
#define KF_HINIT(kf, _ny)                                                    \
{                                                                            \
    kf->ny = _ny;                                                            \
    if(cfg.max_ny < kf->ny)                                                  \
        LOG_FATAL("%f: cfg.max_ny should set large than %i",                 \
                __FUNCTION__, kf->ny);                                       \
    for(unsigned int i = 0; i < kf->ny*kf->nx; ++i)                          \
        kf->H[i] = 0.0;                                                      \
}

/**
 * @brief copy varainace(stanadard error) from inskf->P to inskf->sol.
 * @param   inskf   kalman filter struct
 */
static void kf_copyQ2sol(kf_t *inskf)
{
    #define GETSTD(x)   sqrt(inskf->P[(x)+(x)*inskf->nx])
    m3_copy(&inskf->sol->Qpos, &inskf->P[cfg.IPOS+cfg.IPOS*inskf->nx],
            inskf->nx);
    m3_copy(&inskf->sol->Qvel, &inskf->P[cfg.IVEL+cfg.IVEL*inskf->nx],
            inskf->nx);
    m3_copy(&inskf->sol->Qatt, &inskf->P[cfg.IATT+cfg.IATT*inskf->nx],
            inskf->nx);
    if(cfg.isx_ba){
        inskf->sol->ba_std.x = GETSTD(cfg.IBA);
        inskf->sol->ba_std.y = GETSTD(cfg.IBA+1);
        inskf->sol->ba_std.z = GETSTD(cfg.IBA+2);
    }
    if(cfg.isx_bg){
        inskf->sol->bg_std.x = GETSTD(cfg.IBG);
        inskf->sol->bg_std.y = GETSTD(cfg.IBG+1);
        inskf->sol->bg_std.z = GETSTD(cfg.IBG+2);
    }
    if(cfg.isx_kax)     inskf->sol->ka_std.x = GETSTD(cfg.IKAx);
    if(cfg.isx_kay)     inskf->sol->ka_std.y = GETSTD(cfg.IKAy);
    if(cfg.isx_kaz)     inskf->sol->ka_std.z = GETSTD(cfg.IKAz);
    if(cfg.isx_kgx)     inskf->sol->kg_std.x = GETSTD(cfg.IKGx);
    if(cfg.isx_kgy)     inskf->sol->kg_std.y = GETSTD(cfg.IKGy);
    if(cfg.isx_kgz)     inskf->sol->kg_std.z = GETSTD(cfg.IKGz);
    if(cfg.isx_eroll)   inskf->sol->std_Cbc.x = GETSTD(cfg.IEROLL);
    if(cfg.isx_epitch)  inskf->sol->std_Cbc.y = GETSTD(cfg.IEPITCH);
    if(cfg.isx_eyaw)    inskf->sol->std_Cbc.z = GETSTD(cfg.IEYAW);
    if(cfg.isx_kod)     inskf->sol->std_kod = GETSTD(cfg.IKOD);
    #undef GETSTD
}

/**
 * @brief form transfer noise matrix of kalman matrix
 *          P = F * P * F' + Q
 * @param[out]  Q           transer noise matrix
 * @param[in]   Cbe         Attitude transformation from b to e
 * @param[in]   dt_zoom     time interval[s] x zoom
 * @param[in]   imup        imu property struct
 * @return status(0: OK)
 */
static int
kf_formQ(double *Q, const m3_t *Cbe, double dt_zoom, const imup_t *imup)
{
    for(int i = 0; i < cfg.nx * cfg.nx; ++i)
        Q[i] = 0.0;
    /* Angular randon walk */
    m3_t Qphi = v3_diag(v3_scalar(dt_zoom, v3_pow(imup->arw, 2.0)));
    if(!cfg.is_imu_samenoise && Cbe != NULL)
        Qphi = m3_mul(m3_mul(*Cbe, Qphi), m3_T(*Cbe));
    m3_paste((Q+cfg.IATT+cfg.IATT*cfg.nx), cfg.nx, &Qphi);

    /* Velocity randon walk */
    m3_t Qv = v3_diag(v3_scalar(dt_zoom, v3_pow(imup->vrw, 2.0)));
    if(!cfg.is_imu_samenoise && Cbe != NULL)
        Qv = m3_mul(m3_mul(*Cbe, Qv), m3_T(*Cbe));
    m3_paste((Q+cfg.IVEL+cfg.IVEL*cfg.nx), cfg.nx, &Qv);

    if(cfg.isx_ba){
        /* Accel bias random walk */
        m3_t Qba = v3_diag(v3_scalar(dt_zoom, v3_pow(imup->vrrw, 2.0)));
        m3_paste((Q+cfg.IBA+cfg.IBA*cfg.nx), cfg.nx, &Qba);
    }
    if(cfg.isx_bg){
        /* Gyro bias random walk */
        m3_t Qbg = v3_diag(v3_scalar(dt_zoom, v3_pow(imup->arrw, 2.0)));
        m3_paste((Q+cfg.IBG+cfg.IBG*cfg.nx), cfg.nx, &Qbg);
    }

    /* gyro install error angle random walk */
    if(cfg.isx_eroll)
        Q[cfg.IEROLL+cfg.IEROLL*cfg.nx] = SQR(imup->err_angle_imu_rw.x)*dt_zoom;
    if(cfg.isx_epitch)
        Q[cfg.IEPITCH+cfg.IEPITCH*cfg.nx]=SQR(imup->err_angle_imu_rw.y)*dt_zoom;
    if(cfg.isx_eyaw)
        Q[cfg.IEYAW+cfg.IEYAW*cfg.nx] =   SQR(imup->err_angle_imu_rw.z)*dt_zoom;

    return 0;
}

/**
 * @brief initial ins kalman filter
 * @param[in,out]   inskf   kalman filter of ins
 * @param[in]       imup    imu property sturct
 * @return status(0: OK)
 */
extern int inskf_init(kf_t *inskf, const imup_t *imup)
{
    LOG_CURTIME = imup->tstart;
    inskf->time = imup->tstart;
    inskf->nx = cfg.nx;
    if(imup->freq_imu > 0)
        inskf->idt = 1.0 / (double)imup->freq_imu;
    else
        LOG_FATAL("%s: it is impossible that imup->freq_imu = 0", __FUNCTION__);
    LOG_INFO("kf_init: inskf->idt = %.3f", inskf->idt);
    if(imup->freq_od > 0)
        inskf->odt = 1.0 / (double)imup->freq_od;
    else
        inskf->odt = 0.0;
    LOG_INFO("kf_init: inskf->odt = %.3f", inskf->odt);

    if(!(inskf->x = (double *)malloc(inskf->nx*sizeof(double))))
        LOG_FATAL("%s: inskf->x memory allocation error", __FUNCTION__);
    if(!(inskf->sol = (solins_t *)malloc(sizeof(solins_t))))
        LOG_FATAL("%s: inskf->sol memory allocation error", __FUNCTION__);
    if(!(inskf->P = (double *)calloc(SQR(inskf->nx),sizeof(double))))
        LOG_FATAL("%s: inskf->P memory allocation error", __FUNCTION__);
    if(!(inskf->F = (double *)calloc(SQR(inskf->nx),sizeof(double))))
        LOG_FATAL("%s: inskf->F memory allocation error", __FUNCTION__);
    if(!(inskf->Q = (double *)calloc(SQR(inskf->nx), sizeof(double))))
        LOG_FATAL("%s: inskf->Q memory allocation error", __FUNCTION__);

    /* init x */
    for(int i = 0; i < inskf->nx; ++i)
        inskf->x[i] = 1e-32;    /* x[i] should NOT be 0.0, see filter() */
    /* init sol_t */
    inskf->sol->time = imup->tstart;
    inskf->sol->status = SOL_MANUAL;
    euler2dcm(&imup->inita, &inskf->sol->dcm);
    euler2quat(&imup->inita, &inskf->sol->quat);
    inskf->sol->vel = imup->initv;
    inskf->sol->pos = imup->initr;
    /* init P */
    m3_paste((inskf->P+cfg.IATT+cfg.IATT*cfg.nx), inskf->nx, &imup->initQa);
    m3_paste((inskf->P+cfg.IVEL+cfg.IVEL*cfg.nx), inskf->nx, &imup->initQv);
    m3_paste((inskf->P+cfg.IPOS+cfg.IPOS*cfg.nx), inskf->nx, &imup->initQr);
    /* init F (Some elements keep constant at equal time inteval */
    m3_t m3F;
    jacobi_trans_Ebe2Ebe(&m3F, inskf->idt);
    m3_paste(inskf->F, inskf->nx, &m3F);
    jacobi_trans_veb_e2veb_e(&m3F, inskf->idt);
    m3_paste((inskf->F+cfg.IVEL+cfg.IVEL*inskf->nx), inskf->nx, &m3F);
    jacobi_trans_reb_e2veb_e(&m3F, inskf->idt);
    m3_paste((inskf->F+cfg.IPOS+cfg.IVEL*inskf->nx), inskf->nx, &m3F);
    for(int i = cfg.IPOS; i < cfg.IPOS + 3; ++i)
        inskf->F[i + i*inskf->nx] = 1.0;

    if(cfg.isx_ba){
        inskf->sol->ba = imup->ba;
        m3_t initQba = v3_diag(v3_pow(imup->ba_std, 2.0));
        m3_paste((inskf->P+cfg.IBA+cfg.IBA*cfg.nx), inskf->nx, &initQba);
        jacobi_trans_markov(&m3F, inskf->idt,&imup->Ta);
        m3_paste((inskf->F+cfg.IBA+cfg.IBA*inskf->nx), inskf->nx, &m3F);
    }
    if(cfg.isx_bg){
        inskf->sol->bg = imup->bg;
        m3_t initQbg = v3_diag(v3_pow(imup->bg_std, 2.0));
        m3_paste((inskf->P+cfg.IBG+cfg.IBG*cfg.nx), inskf->nx, &initQbg);
        jacobi_trans_markov(&m3F, inskf->idt, &imup->Tg);
        m3_paste((inskf->F+cfg.IBG+cfg.IBG*inskf->nx), inskf->nx, &m3F);
    }

    if(v3_equal(&imup->ka, &V0, 1e-2)){
        LOG_WARN("%s: imup->ka shouldn't set to zero, use default: 1.0",
                  __FUNCTION__);
        inskf->sol->ka = V1;
    }else{
        inskf->sol->ka = imup->ka;
    }
    if(v3_equal(&imup->kg, &V0, 1e-2)){
        LOG_WARN("%s: imup->kg shouldn't set to zero, use default: 1.0",
                 __FUNCTION__);
        inskf->sol->kg = V1;
    }else{
        inskf->sol->kg = imup->kg;
    }
    if(cfg.isx_kax){
        inskf->F[cfg.IKAx + cfg.IKAx*cfg.nx] = 1.0;
        inskf->P[cfg.IKAx + cfg.IKAx*cfg.nx] = SQR(imup->ka_std.x);
    }
    if(cfg.isx_kay){
        inskf->F[cfg.IKAy + cfg.IKAy*cfg.nx] = 1.0;
        inskf->P[cfg.IKAy + cfg.IKAy*cfg.nx] = SQR(imup->ka_std.y);
    }
    if(cfg.isx_kaz){
        inskf->F[cfg.IKAz + cfg.IKAz*cfg.nx] = 1.0;
        inskf->P[cfg.IKAz + cfg.IKAz*cfg.nx] = SQR(imup->ka_std.z);
    }
    if(cfg.isx_kgx){
        inskf->F[cfg.IKGx + cfg.IKGx*cfg.nx] = 1.0;
        inskf->P[cfg.IKGx + cfg.IKGx*cfg.nx] = SQR(imup->kg_std.x);
    }
    if(cfg.isx_kgy){
        inskf->F[cfg.IKGy + cfg.IKGy*cfg.nx] = 1.0;
        inskf->P[cfg.IKGy + cfg.IKGy*cfg.nx] = SQR(imup->kg_std.y);
    }
    if(cfg.isx_kgz){
        inskf->F[cfg.IKGz + cfg.IKGz*cfg.nx] = 1.0;
        inskf->P[cfg.IKGz + cfg.IKGz*cfg.nx] = SQR(imup->kg_std.z);
    }
    if(imup->kod != 0.0){
        inskf->sol->kod = imup->kod;
    }else {
        LOG_WARN("%s: imup->kod shouldn't set to zero, use default: 1.0",
                 __FUNCTION__);
        inskf->sol->kod = 1.0;
    }
    if(cfg.isx_kod){
        inskf->P[cfg.IKOD+cfg.IKOD*cfg.nx] = SQR(imup->kod_std);
        inskf->F[cfg.IKOD+cfg.IKOD*cfg.nx] = 1.0;
    }
    LOG_INFO("kf_init: set install error angle: roll(%.2f deg), "
             "pitch(%.2f deg), yaw(%.2f deg)", imup->err_angle_imu.x*RAD2DEG,
             imup->err_angle_imu.y*RAD2DEG, imup->err_angle_imu.z*RAD2DEG);
    euler2dcm(&imup->err_angle_imu, &inskf->sol->Cbc);
    if(cfg.isx_eroll){
        LOG_INFO("kf_init: inskf->P of install roll error angle: %.2f deg",
                 imup->err_angle_imu_std.x*RAD2DEG);
        inskf->P[cfg.IEROLL+cfg.IEROLL*cfg.nx] = SQR(imup->err_angle_imu_std.x);
        inskf->F[cfg.IEROLL+cfg.IEROLL*cfg.nx] = 1.0;
//        if(imup->Terr_angle_imu.x == 0.0)
//            inskf->F[cfg.IEROLL+cfg.IEROLL*cfg.nx] = 1.0;
//        else
//            inskf->F[cfg.IEROLL+cfg.IEROLL*cfg.nx] =
//                    1.0 - inskf->idt / imup->Terr_angle_imu.x;
    }
    if(cfg.isx_epitch){
        LOG_INFO("kf_init: inskf->P of install pitch error angle: %.2f deg",
                 imup->err_angle_imu_std.y*RAD2DEG);
        inskf->P[cfg.IEPITCH+cfg.IEPITCH*cfg.nx] = SQR(imup->err_angle_imu_std.y);
        inskf->F[cfg.IEPITCH+cfg.IEPITCH*cfg.nx] = 1.0;
//        if(imup->Terr_angle_imu.y == 0.0)
//            inskf->F[cfg.IEPITCH+cfg.IEPITCH*cfg.nx] = 1.0;
//        else
//            inskf->F[cfg.IEPITCH+cfg.IEPITCH*cfg.nx] =
//                    1.0 - inskf->idt / imup->Terr_angle_imu.y;
    }
    if(cfg.isx_eyaw){
        LOG_INFO("kf_init: inskf->P of install yaw error angle: %.2f deg",
                 imup->err_angle_imu_std.z*RAD2DEG);
        inskf->P[cfg.IEYAW+cfg.IEYAW*cfg.nx] = SQR(imup->err_angle_imu_std.z);
        inskf->F[cfg.IEYAW+cfg.IEYAW*cfg.nx] = 1.0;
//        if(imup->Terr_angle_imu.z == 0.0)
//            inskf->F[cfg.IEYAW+cfg.IEYAW*cfg.nx] = 1.0;
//        else
//            inskf->F[cfg.IEYAW+cfg.IEYAW*cfg.nx] =
//                    1.0 - inskf->idt / imup->Terr_angle_imu.z;
    }
    if(cfg.is_odincre){
        if(!(inskf->itg = (double *)calloc(cfg.nitg, sizeof(double))))
            LOG_FATAL("%s: inskf->itg memory allocation error",__FUNCTION__);
        LOG_INFO("kf_init: start intergral variables ...");
        inskf->itg_start = inskf->time;
    }
    if(cfg.max_ny > 0){
        LOG_INFO("kf_init: inskf->H allocate memory to %i x %i",
                 cfg.max_ny, cfg.nx);
        if(!(inskf->H = (double *)malloc(cfg.max_ny*cfg.nx*sizeof(double))))
            LOG_FATAL("%s: inskf->H memory allocation error",__FUNCTION__);
    }
    if(cfg.nR > 0){
        LOG_INFO("kf_init: inskf->R allocate memory to %i", cfg.nR);
        if(!(inskf->R = (double *)malloc(cfg.nR*sizeof(double))))
            LOG_FATAL("%s: inskf->R memory allocation error",__FUNCTION__);
        if(cfg.iskf_itgdS_sagehusa){
            for(int i = 0; i < 3; ++i)
                inskf->R[cfg.IR_itgdS+i] = 0.1;
        }
    }
    if(cfg.nZST > 0){
        LOG_INFO("kf_init: inskf->imud allocate memory to %i", cfg.nZST);
        if(!(inskf->imud = (imud_t *)malloc(cfg.nZST*sizeof(imud_t))))
            LOG_FATAL("%s: inskf->imud memory allocation error",__FUNCTION__);
        inskf->nimud = 0;
        inskf->imudend = cfg.nZST - 1;
    }
    kf_copyQ2sol(inskf);

    if(!cfg.is_imu_samenoise){
        if( fabs(imup->arw.x - imup->arw.y) < 1e-32 &&
            fabs(imup->arw.x - imup->arw.z) < 1e-32 &&
            fabs(imup->vrw.x - imup->vrw.y) < 1e-32 &&
            fabs(imup->vrw.x - imup->vrw.z) < 1e-32)
            cfg.is_imu_samenoise = true;
    }
    if(cfg.is_imu_samenoise){
        LOG_INFO("kf_init: cfg.is_imu_samenoise = true, speed up for Q");
        kf_formQ(inskf->Q, NULL, inskf->idt*1.0, imup);
    }

    LOG_INFO("kf_init: kalman filter initialization end");
    return 0;
}

/**
 * @brief kalman filter update state
 *      x = Fx,   P = F * P * F' +  Q, and full state transfer
 * @param[in] inskf     kalman filter
 * @param[in] imud      imu data struct
 * @param[in] imup      imu property struct
 * @return status(0:OK)
 */
extern int inskf_udstate(kf_t *inskf, const imud_t *imud, const imup_t *imup)
{
    gtime_t cur_time = imud->time;
    double dt = timediff(cur_time, inskf->time);
    if(fabs(dt - inskf->idt) > 1E-8){
        LOG_WARN("%s: different time interval: %.3f(config: %.3f)",__FUNCTION__,
                 dt, inskf->idt);
        /* autofix for P2 time problem */
        if(fabs(dt - inskf->idt*2) < 1e-4){
            cur_time = timeadd(inskf->time, inskf->idt);
            dt = inskf->idt;
            LOG_WARN("%s: Autofix the time error", __FUNCTION__);
        }
    }
    /* corrent imu output */
    v3_t dv, dtheta;
    if(cfg.isx_ba) dv = v3_del(imud->accel, v3_scalar(dt, inskf->sol->ba));
    else dv = imud->accel;
    if(cfg.isx_bg) dtheta = v3_del(imud->gyro, v3_scalar(dt, inskf->sol->bg));
    else dtheta = imud->gyro;

    dv =        v3_mul_cxc(inskf->sol->ka, dv);
    dtheta =    v3_mul_cxc(inskf->sol->kg, dtheta);

    LOG_CURTIME = cur_time;
    inskf->time = cur_time;

    /* Update F */
    m3_t m3F;
    m3_t Cbe = inskf->sol->dcm;
    jacobi_trans_veb_e2Ebe(&m3F, &Cbe, &dv);
    m3_paste(inskf->F + cfg.IVEL + cfg.IATT*inskf->nx, inskf->nx, &m3F);
    jacobi_trans_veb_e2reb_e(&m3F, dt, &inskf->sol->pos);
    m3_paste(inskf->F + cfg.IVEL + cfg.IPOS*inskf->nx, inskf->nx, &m3F);
    if(cfg.isx_ba){
        jacobi_trans_veb_e2ba(&m3F, dt, &Cbe);
        m3_paste(inskf->F + cfg.IVEL + cfg.IBA*inskf->nx, inskf->nx, &m3F);
    }
    if(cfg.isx_bg){
        jacobi_trans_Ebe2bg(&m3F, dt, &Cbe);
        m3_paste(inskf->F + cfg.IATT + cfg.IBG*inskf->nx, inskf->nx, &m3F);
    }
    if(cfg.isx_kax) jacobi_trans_veb_e2kax(inskf->F, &Cbe, &dv);
    if(cfg.isx_kay) jacobi_trans_veb_e2kay(inskf->F, &Cbe, &dv);
    if(cfg.isx_kaz) jacobi_trans_veb_e2kaz(inskf->F, &Cbe, &dv);
    if(cfg.isx_kgx) jacobi_trans_Ebe2kgx(inskf->F, &Cbe, &dtheta);
    if(cfg.isx_kgy) jacobi_trans_Ebe2kgy(inskf->F, &Cbe, &dtheta);
    if(cfg.isx_kgz) jacobi_trans_Ebe2kgz(inskf->F, &Cbe, &dtheta);

    /* Change to Markov model after P is small enough */
    if(cfg.isx_eroll && imup->Terr_angle_imu.x != 0.0){
        double std=markov_rw2std(imup->err_angle_imu_rw.x,imup->Terr_angle_imu.x);
        if(sqrt(inskf->P[cfg.IEROLL + cfg.IEROLL*cfg.nx]) < std){
            inskf->F[cfg.IEROLL+cfg.IEROLL*cfg.nx] =
                    1.0 - inskf->idt / imup->Terr_angle_imu.x;
        }
    }
    if(cfg.isx_epitch && imup->Terr_angle_imu.y != 0.0){
        double std=markov_rw2std(imup->err_angle_imu_rw.y,imup->Terr_angle_imu.y);
        if(sqrt(inskf->P[cfg.IEPITCH + cfg.IEPITCH*cfg.nx]) < std){
            inskf->F[cfg.IEPITCH+cfg.IEPITCH*cfg.nx] =
                    1.0 - inskf->idt / imup->Terr_angle_imu.y;
        }
    }
    if(cfg.isx_eyaw && imup->Terr_angle_imu.z != 0.0){
        double std=markov_rw2std(imup->err_angle_imu_rw.z,imup->Terr_angle_imu.z);
        if(sqrt(inskf->P[cfg.IEYAW + cfg.IEYAW*cfg.nx]) < std){
            inskf->F[cfg.IEYAW+cfg.IEYAW*cfg.nx] =
                    1.0 - inskf->idt / imup->Terr_angle_imu.z;
        }
    }

//    matprint(inskf->F, inskf->nx, inskf->nx, 5, 2);
//    fflush(stdout);

    /* x = F*x */
    if((double)cfg.feedratio < 1.0){
        double *_x = mat(inskf->nx, 1);
        matmul("NN", inskf->nx, 1, inskf->nx, 1.0, inskf->F, inskf->x, 0.0, _x);
        matcpy(inskf->x, _x, inskf->nx, 1);
    }

    /* update var-covarinace: P = F*P*F' + Q */
    double *FP = mat(inskf->nx, inskf->nx);
    matmul("NN", inskf->nx, inskf->nx, inskf->nx, 1.0, inskf->F, inskf->P, 0.0, FP);
    matmul("NT", inskf->nx, inskf->nx, inskf->nx, 1.0, FP, inskf->F, 0.0, inskf->P);
    free(FP);

    /* process noise of INS */
    if(!cfg.is_imu_samenoise){
        double noise_zoom = 1.0;
        kf_formQ(inskf->Q, &Cbe, dt * noise_zoom, imup);
    }

    for(int i = 0; i < SQR(inskf->nx); ++i)
        inskf->P[i] += inskf->Q[i];

    v3_t last_pos = inskf->sol->pos;
    if(ins_nav_ecef(dt, &dtheta, &dv, &inskf->sol->pos, &inskf->sol->vel,
                          &inskf->sol->quat)){
        LOG_FATAL("%s: failed to update ins_nav_ecef", __FUNCTION__);
    }else{
        inskf->sol->time = inskf->time;
        inskf->sol->status = SOL_INS;
        quat2dcm(&inskf->sol->quat, &inskf->sol->dcm);
        kf_copyQ2sol(inskf);
    }

    if(cfg.is_odincre){
        v3_t dreb_e = v3_del(inskf->sol->pos, last_pos);
        inskf->itg[cfg.ITGVEL_INS  ] += dreb_e.x;
        inskf->itg[cfg.ITGVEL_INS+1] += dreb_e.y;
        inskf->itg[cfg.ITGVEL_INS+2] += dreb_e.z;
        v3_t vib_e = quat_mul_v3(inskf->sol->quat, imud->accel);
        inskf->itg[cfg.ITGF_INS] += vib_e.x;
        inskf->itg[cfg.ITGF_INS] += vib_e.y;
        inskf->itg[cfg.ITGF_INS] += vib_e.z;
    }

    if(cfg.nZST > 0){
        /* add imud data to inskf->imud queue */
        inskf->imudend ++;
        if(inskf->imudend >= cfg.nZST) inskf->imudend -= cfg.nZST;
        if(inskf->nimud < cfg.nZST) inskf->nimud ++;
        inskf->imud[inskf->imudend] = *imud;
    }

    return 0;
}

/**
 * @brief kalman filter position measment update
 * @param[in,out]   inskf       ins kalman filter struct
 * @param[in]       reg_e       (gnss) position under e-frame[m]
 * @param[in]       Qreg_e      (gnss) position variance matrix
 * @param[in]       lever_arm   lever arm of (gnss) postion(under b-frame)
 * @return  status(0: OK)
 */
extern int inskf_udmeasr(kf_t *inskf, const v3_t *reg_e, const m3_t *Qreg_e,
                      const v3_t *lever_arm)
{
    KF_HINIT(inskf, 3);
    m3_t Cbe = inskf->sol->dcm;

    /* reb_e measments update */
    v3_t reb_e_gnss = v3_del(*reg_e, m3_mul_v3(Cbe, *lever_arm));
    v3_t dz_r = v3_del(reb_e_gnss, inskf->sol->pos);

    LOG_TRACE("dz_rx %f, dz_ry %f, dz_rz %f", dz_r.x, dz_r.y, dz_r.z);

    m3_t m3H;
    jacobi_meas_reb_e2Ebe(&m3H, &Cbe, lever_arm);
    m3H = m3_T(m3H);
    m3_paste((inskf->H+cfg.IATT), inskf->nx, &m3H);
    inskf->H[cfg.IPOS              ] = -1.0;
    inskf->H[cfg.IPOS+1+1*inskf->nx] = -1.0;
    inskf->H[cfg.IPOS+2+2*inskf->nx] = -1.0;

    /* innovtaion check */
    double ndz_r[3];
    inskf_norm_innov(inskf, (const double *)&dz_r,(const double *)Qreg_e,
                     ndz_r);
    LOG_TRACE("ndz_rx %f, ndz_ry %f, ndz_rz %f", ndz_r[0], ndz_r[1], ndz_r[2]);
    for(unsigned int i = 0; i < inskf->ny; ++i){
        if(fabs(ndz_r[i]) > 5.0){
            const double *pdz_r = (const double *)&dz_r;
            LOG_WARN("%s: too small noise(%f m) or too large innovations(%f m)",
                     __FUNCTION__, pdz_r[i]/ndz_r[i], pdz_r[i]);
        }
    }

    /* note: Qreg_e is sysmmetry matrix */
    KF_FILTER(inskf, (const double *)&dz_r, (const double *)Qreg_e);

    return 0;
}

/**
 * @brief kalman filter velocity update
 * @param[in,out]   inskf       ins kalman filter struct
 * @param[in]       veg_e       (gnss) velocity under e-frame[m/s]
 * @param[in]       Qveg_e      (gnss) velocity variance matrix
 * @param[in]       lever_arm   lever arm of (gnss) position(unber b-frame)
 * @param[in]       wib_b       angular rate of imu
 * @return status(0: OK)
 */
extern int inskf_udmeasv(kf_t *inskf, const v3_t *veg_e, const m3_t *Qveg_e,
                      const v3_t *lever_arm, const v3_t *wib_b)
{
    KF_HINIT(inskf, 3);
    m3_t Cbe = inskf->sol->dcm;

    v3_t v1 = m3_mul_v3(Cbe, v3_cross(*wib_b, *lever_arm));
    v3_t wie_e = {0.0, 0.0, wgs84.wie};
    v3_t v2 = v3_cross(wie_e, m3_mul_v3(Cbe, *lever_arm));
    v3_t veb_e_gnss = v3_add(v3_del(*veg_e, v1), v2);
    v3_t dz_v = v3_del(veb_e_gnss, inskf->sol->vel);

    LOG_TRACE("dz_vx %f, dz_vy %f, dz_vz %f", dz_v.x, dz_v.y, dz_v.z);

    m3_t m3H;
    jacobi_meas_veb_e2Ebe(&m3H, &Cbe, wib_b, lever_arm);
    m3H = m3_T(m3H);
    m3_paste((inskf->H+cfg.IATT), inskf->nx, &m3H);
    if(cfg.isx_bg){
        jacobi_meas_veb_e2bg(&m3H, &Cbe, lever_arm);
        m3H = m3_T(m3H);
        m3_paste((inskf->H+cfg.IBG), inskf->nx, &m3H);
    }
    inskf->H[cfg.IVEL              ] = -1.0;
    inskf->H[cfg.IVEL+1+1*inskf->nx] = -1.0;
    inskf->H[cfg.IVEL+2+2*inskf->nx] = -1.0;

    double ndz_v[3];
    inskf_norm_innov(inskf, (const double *)&dz_v, (const double *)Qveg_e,
                     ndz_v);
    LOG_TRACE("ndz_rx %f, ndz_ry %f, ndz_rz %f", ndz_v[0], ndz_v[1], ndz_v[2]);
    for(unsigned int i = 0; i < inskf->ny; ++i){
        if(fabs(ndz_v[i]) > 5.0){
            const double *pdz_v = (const double *)&dz_v;
            LOG_WARN("%s: too small noise(%f m/s) or too large "
                     "innovations(%f m/s)",
                     __FUNCTION__, pdz_v[i]/ndz_v[i], pdz_v[i]);
        }
    }

    KF_FILTER(inskf, (const double *)&dz_v, (const double *)Qveg_e);
//    matprint(inskf->x, 1, inskf->nx, 8, 5);
//    fflush(stdout);

    return 0;
}

/**
 * @brief kalman filter yaw angle measument update
 * @param[in,out]   inskf   ins kalman fitler struct
 * @param[in]       yaw     yaw angle measuremnt(Enb.k) [rad]
 * @param[in]       Qyaw    uncertainty of yaw [rad^2]
 * @return status(0: OK, 1: failed)
 */
extern int inskf_udmeas_yaw(kf_t *inskf, double yaw, double Qyaw)
{
    KF_HINIT(inskf, 1);
    v3_t att = Cbe2att(&inskf->sol->pos, &inskf->sol->dcm);
    double dz_yaw = yaw_del(yaw, att.z);

    v3_t pos = inskf->sol->pos; ecef2ned(&pos, NULL, NULL);
    m3_t Cen = formCen_ned(pos.x, pos.y);

    inskf->H[cfg.IATT  ] = -Cen.m31;
    inskf->H[cfg.IATT+1] = -Cen.m32;
    inskf->H[cfg.IATT+2] = -Cen.m33;

    double ndz_yaw;
    inskf_norm_innov(inskf, &dz_yaw, &Qyaw, &ndz_yaw);

    LOG_TRACE("dz_yaw %f, ndz_yaw %f", dz_yaw*RAD2DEG, ndz_yaw);
    if(fabs(ndz_yaw) > 5.0){
        LOG_WARN("%s: disabled, too small yaw noise(%f deg) or large "
                 "innovations(%f deg)",
                 __FUNCTION__, (dz_yaw/ndz_yaw)*RAD2DEG, dz_yaw*RAD2DEG);
        return 1;
    }

    KF_FILTER(inskf, &dz_yaw, &Qyaw);
    return 0;
}

/**
 * @brief ins kalman filter Zero Angular Update
 * @param[in,out]   inskf   ins kalman filter struct
 * @param[in]       yaw     yaw measment at start of static moment[rad]
 * @param[in]       Qyaw    uncertainty of yaw[rad^2]
 * @return status(0: OK, 1: failed, current momoent may be not static)
 * @note this function nearly the same as inskf_udmeas_yaw(), but do not have
 *      same function.
 */
extern int inskf_ZAU(kf_t *inskf, double yaw, double Qyaw)
{
    KF_HINIT(inskf, 1);
    v3_t att = Cbe2att(&inskf->sol->pos, &inskf->sol->dcm);
    double dz_yaw = yaw_del(yaw, att.z);

    v3_t pos = inskf->sol->pos; ecef2ned(&pos, NULL, NULL);
    m3_t Cen = formCen_ned(pos.x, pos.y);

    inskf->H[cfg.IATT  ] = -Cen.m31;
    inskf->H[cfg.IATT+1] = -Cen.m32;
    inskf->H[cfg.IATT+2] = -Cen.m33;

    double ndz_yaw;
    inskf_norm_innov(inskf, &dz_yaw, &Qyaw, &ndz_yaw);

    LOG_TRACE("dz_yaw %f, ndz_yaw %f", dz_yaw*RAD2DEG, ndz_yaw);
    if(fabs(ndz_yaw) > 5.0){
        LOG_WARN("%s: disabled, too small yaw noise(%f deg) or large innovations(%f deg)",
                 __FUNCTION__, (dz_yaw/ndz_yaw)*RAD2DEG, dz_yaw*RAD2DEG);
        return 1;
    }

    KF_FILTER(inskf, &dz_yaw, &Qyaw);
    return 0;
}

/**
 * @brief ins kalman filter update attitude measment.
 * @param[in,out]   inskf   ins kalman filter struct
 * @param[in]       Cbe     DCM from b-frame to e-frame
 * @param[in]       QEbe    uncertainty of Ebe. [rad^2]
 * @return status(0: OK)
 */
extern int inskf_udmeas_Cbe(kf_t *inskf, const m3_t *Cbe, const m3_t *QEbe)
{

    KF_HINIT(inskf, 3);
//    m3_t dz_Cbe = m3_mul(*Cbe, m3_T(inskf->sol->dcm));
    m3_t dz_Cbe = m3_mul(inskf->sol->dcm, m3_T(*Cbe));
    v3_t dz_Ebe;  dcm2euler(&dz_Cbe, &dz_Ebe);
    if(dz_Ebe.z > PI) dz_Ebe.z -= 2*PI;

    LOG_TRACE("dz_Ebex %f, dz_Ebey %f, dz_Ebez %f",
              dz_Ebe.x*RAD2DEG, dz_Ebe.y*RAD2DEG, dz_Ebe.z*RAD2DEG);

    inskf->H[cfg.IATT           ] = - 1.0;
    inskf->H[cfg.IATT+1+1*cfg.nx] = - 1.0;
    inskf->H[cfg.IATT+2+2*cfg.nx] = - 1.0;

    KF_FILTER(inskf, (const double *)&dz_Ebe, (const double *)QEbe);

//    matprint((const double *)&dz_Ebe, 1,3,12,6);
//    matprint(inskf->x, 1, inskf->nx, 12, 6);
    return 0;
}

/**
 * @brief ins kalman filter feedback (move inskf->x to inskf.sol)
 * @param[in,out]   inskf       ins kalman filter
 * @param[in]       SOL_TYPE    solution type, see enum SOL
 * @return status(0: OK)
 */
extern int inskf_feedback(kf_t *inskf, unsigned int SOL_TYPE)
{
    if(cfg.feedratio < 1)
        LOG_FATAL("%s: Do not support cfg.feedratio < 1", __FUNCTION__);

    /* attitude feedback */
    inskf->sol->time = inskf->time;
    inskf->sol->status = SOL_TYPE;

    if(norm(inskf->x+cfg.IATT, 3) > 10.0*DEG2RAD){
        LOG_WARN("%s: Too large residuals for attitude update", __FUNCTION__);
    }

    m3_t phix = v3_askew(*(v3_t *)(inskf->x+cfg.IATT));
    inskf->sol->dcm = m3_mul(m3_del(I3, phix), inskf->sol->dcm);
    dcm2quat(&inskf->sol->dcm, &inskf->sol->quat);
    for(int i = cfg.IATT; i < 3+cfg.IATT; ++i)
        inskf->x[i] = 1e-32;
    inskf->sol->vel = v3_del(inskf->sol->vel, *(v3_t*)(inskf->x+cfg.IVEL));
    for(int i = cfg.IVEL; i < 3+cfg.IVEL; ++i)
        inskf->x[i] = 1e-32;
    inskf->sol->pos = v3_del(inskf->sol->pos, *(v3_t*)(inskf->x+cfg.IPOS));
    for(int i = cfg.IPOS; i < 3+cfg.IPOS; ++i)
        inskf->x[i] = 1e-32;

    if(cfg.isx_kod){
        inskf->sol->kod -= inskf->x[cfg.IKOD];
        inskf->x[cfg.IKOD] = 1e-32;
    }

    if(cfg.isx_eyaw || cfg.isx_eroll || cfg.isx_epitch){
        v3_t Ebc = V0;
        if(cfg.isx_eroll){
            Ebc.x = inskf->x[cfg.IEROLL];
            inskf->x[cfg.IEROLL] = 1e-32;
        }
        if(cfg.isx_epitch){
            Ebc.y = inskf->x[cfg.IEPITCH];
            inskf->x[cfg.IEPITCH] = 1e-32;
        }
        if(cfg.isx_eyaw){
            Ebc.z = inskf->x[cfg.IEYAW];
            inskf->x[cfg.IEYAW] = 1e-32;
        }
        inskf->sol->Cbc = m3_mul(inskf->sol->Cbc,m3_del(I3,v3_askew(Ebc)));
    }

    if(cfg.isx_kax){
        inskf->sol->ka.x -= inskf->x[cfg.IKAx];
        inskf->x[cfg.IKAx] = 1e-32;
    }
    if(cfg.isx_kay){
        inskf->sol->ka.y -= inskf->x[cfg.IKAy];
        inskf->x[cfg.IKAy] = 1e-32;
    }
    if(cfg.isx_kaz){
        inskf->sol->ka.z -= inskf->x[cfg.IKAz];
        inskf->x[cfg.IKAz] = 1e-32;
    }
    if(cfg.isx_kgy){
        inskf->sol->kg.x -= inskf->x[cfg.IKGx];
        inskf->x[cfg.IKGx] = 1e-32;
    }
    if(cfg.isx_kgy){
        inskf->sol->kg.y -= inskf->x[cfg.IKGy];
        inskf->x[cfg.IKGy] = 1e-32;
    }
    if(cfg.isx_kgz){
        inskf->sol->kg.z -= inskf->x[cfg.IKGz];
        inskf->x[cfg.IKGz] = 1e-32;
    }

    double feedratio = (double)cfg.feedratio;
    if(feedratio > 1.0 || feedratio < 0.0){
        LOG_WARN("set kalman filter feedback ratio to 1.0(origin=%f", feedratio);
        feedratio = 1.0;
    }

    if(fabs(feedratio - 1.0) < 1e-16){
        if(cfg.isx_ba){
            inskf->sol->ba = v3_del(inskf->sol->ba, *(v3_t*)(inskf->x+cfg.IBA));
            for(int i = cfg.IBA; i < 3+cfg.IBA; ++i)
                inskf->x[i] = 1e-32;
        }
        if(cfg.isx_bg){
            inskf->sol->bg = v3_del(inskf->sol->bg, *(v3_t*)(inskf->x+cfg.IBG));
            for(int i = cfg.IBG; i < 3+cfg.IBG; ++i)
                inskf->x[i] = 1e-32;
        }
    }else{
        v3_t x3;
        if(cfg.isx_ba){
            x3 = v3_scalar(feedratio, *(v3_t *)(inskf->x+cfg.IBA));
            inskf->sol->ba = v3_del(inskf->sol->ba, x3);
            for(int i = cfg.IBA; i < 3 + cfg.IBA; ++i)
                inskf->x[i] = (1 - feedratio) * inskf->x[i];
        }
        if(cfg.isx_bg){
            x3 = v3_scalar(feedratio, *(v3_t *)(inskf->x+cfg.IBG));
            inskf->sol->bg = v3_del(inskf->sol->bg, x3);
            for(int i = cfg.IBG; i < 3 + cfg.IBG; ++i)
                inskf->x[i] = (1 - feedratio) * inskf->x[i];
        }
    }
    return 0;
}

/**
 * @brief close ins kalman filter,  free all variable memory
 * @param[in]   inskf   ins kalman filter struct
 * @return status(0: OK)
 */
extern int inskf_close(kf_t *inskf)
{
    LOG_INFO("kalman filter closed");
    LOG_CURTIME = (gtime_t){0, 0.0};
    free(inskf->x);     inskf->x = NULL;
    free(inskf->F);     inskf->F = NULL;
    free(inskf->P);     inskf->P = NULL;
    free(inskf->Q);     inskf->Q = NULL;
    free(inskf->sol);   inskf->sol = NULL;
    if(cfg.nitg > 0)    { free(inskf->itg);   inskf->itg = NULL; }
    if(cfg.max_ny > 0)  { free(inskf->H);     inskf->H = NULL; }
    if(cfg.nR > 0)      { free(inskf->R);     inskf->R = NULL; }
    return 0;
}

/**
 * @brief constraint yaw by velocity(Use velocity to caculate yaw, and
 *      then update)
 * @param[in,out]   inskf   ins kalman filter struct
 * @param[in]       wib_b   currnet imu output angular rate
 * @param[in]       imup    imu property struct
 * @return status(0: OK)
 * @warning make sure all necessary imu proerties are well setted.
 * @note related imu property:
 *  imup->phi_install.z         imu install yaw angel error[rad]
 *  imup->phi_install_std.z     imu install yaw angel error variance[rad^2]
 *  imup->lever_arm_car.x       car reference point(forward direction)[m]
 *  imup->lever_arm_car_std.x   car reference point variance[m^2]
 *  imup->gyro_noise.z          gryo downward axis output error[rad/s]
 */
int inskf_constraint_yaw(kf_t *inskf, const v3_t *wib_b, const imup_t *imup)
{
    v3_t pos = inskf->sol->pos, vel = inskf->sol->vel;
    m3_t dcm = inskf->sol->dcm;

    m3_t Qveb_n;
    m3_copy(&Qveb_n, &inskf->P[cfg.IVEL+cfg.IVEL*inskf->nx], inskf->nx);
    ecef2nedQ(&pos, NULL, &Qveb_n, NULL);
    ecef2ned(&pos, &vel, &dcm); /* Cbe => Cbn */

    LOG_DEBUG("Qvebn(%f %f %f)",
              sqrt(Qveb_n.m11), sqrt(Qveb_n.m22),sqrt(Qveb_n.m33));

    double yaw = 0.0, Qyaw = 0.0;
    if(!vel2yaw(&vel, &yaw, &Qveb_n, &Qyaw)){
        m3_t Cen = formCen_ned(pos.x, pos.y);
        v3_t web_n = v3_del(m3_mul_v3(dcm, *wib_b),m3_mul_v3(Cen, WIE_E));
        v3_t veb_b = m3_mul_v3(m3_T(dcm), vel);

        if(veb_b.x < - 0.1){
            LOG_WARN("%s: disabled, car are backing", __FUNCTION__);
            return 1;
        }

        /* adust vecocity */
        m3_t Qveb_b = m3_mul(m3_T(dcm), m3_mul(Qveb_n, dcm));
        double Qv = Qveb_b.m11 < SQR(0.01) ? SQR(0.01) : Qveb_b.m11*2.0;

        double dL = imup->lever_arm_car.x;
        double dyaw_sw = dyaw_swerve(veb_b.x, web_n.z, dL);
        double Qyaw_sw = dyaw_swerve_Q(dyaw_sw, veb_b.x, web_n.z*10.0, dL,
            Qv,SQR(imup->gyro_noise.z), SQR(imup->lever_arm_car_std.x));

        LOG_TRACE("dyaw_swerve %f, dyaw_swerve_std %f, vQyaw %f"
                  ", web_nz %f, veb_bx %f",
                  dyaw_sw*RAD2DEG, sqrt(Qyaw_sw)*RAD2DEG, sqrt(Qyaw)*RAD2DEG,
                  web_n.z, veb_b.x);

        yaw = yaw - imup->err_angle_imu.z - dyaw_sw;
        Qyaw = Qyaw + SQR(imup->err_angle_imu_std.z) + Qyaw_sw;

        inskf_udmeas_yaw(inskf, yaw, Qyaw);
    }
    return 0;
}

/**
 * @brief  Zero angular rate update for ins kalman filter
 * @param[in]   inskf       ins kalman filter struct
 * @param[in]   last_qbe    last quaternion
 * @param[in]   dt          time interval[s]
 * @param[in]   Qweb_b      uncertainty of three axis angle
 * @return status(0: OK)
 */
int inskf_ZRAU(kf_t *inskf, const quat_t *last_qbe, double dt,
                const m3_t *Qweb_b)
{
    KF_HINIT(inskf, 3);
    quat_t last_qeb = *last_qbe; quat_conj(&last_qeb);
    quat_t qbb = quat_mul(last_qeb, inskf->sol->quat);
    v3_t rveb_b; quat2rv(&qbb, &rveb_b);
    v3_t web_b = v3_scalar(1.0/dt, rveb_b);

    inskf->H[cfg.IBG   + 0*inskf->nx] = -1.0;
    inskf->H[cfg.IBG+1 + 1*inskf->nx] = -1.0;
    inskf->H[cfg.IBG+2 + 2*inskf->nx] = -1.0;

    LOG_TRACE("ZRAU_dzx %f, ZRAU_dzy %f, ZRAU_dzz %f",
              web_b.x*RPS2DPH, web_b.y*RPS2DPH, web_b.z*RPS2DPH);

//    LOG_TRACE("ZRAU_ndzx %f, ZRAU_ndzy %f, ZRAU_ndzz %f",
//              ZRAU_ndz[0], ZRAU_ndz[1], ZRAU_ndz[2]);

    double ZRAU_ndz[3];
    inskf_norm_innov(inskf, (const double *)&Qweb_b, (const double *)&web_b,
                     ZRAU_ndz);
    for (int i = 0; i < 3; ++i) {
        if(fabs(ZRAU_ndz[i]) > 3.0){
            LOG_WARN("%s: disabled", __FUNCTION__);
            return 1;
        }
    }
    LOG_INFO("%s: enabled", __FUNCTION__);

    KF_FILTER(inskf, (const double *)&web_b, (const double *)Qweb_b);
    return 0;
}

/**
 * @brief Zero velocity update of ins kalman filter
 * @param[in,out]   inskf   ins kalman filter struct
 * @param[in]       Qveb_e  uncerainty of zero speed
 * @return Status(0: OK)
 */
extern int inskf_ZUPT(kf_t *inskf, const m3_t *Qveb_e)
{
    KF_HINIT(inskf,3);
    v3_t dz_v = inskf->sol->vel;

    inskf->H[cfg.IVEL              ] = 1.0;
    inskf->H[cfg.IVEL+1+1*inskf->nx] = 1.0;
    inskf->H[cfg.IVEL+2+2*inskf->nx] = 1.0;

    KF_FILTER(inskf, (const double*)&dz_v, (const double *)Qveb_e);
    return 0;
}

/**
 * @brief ins kalman filter of update integral varables of odomemter
 * @param[in,out]   inskf   ins kalman filter struct
 * @param[in]   dS      distance increment[s]
 * @param[in]   lever_arm_car   ref point uner b-frame
 * @param[in]   wib_b   imu output angular rate
 * @return status(0: OK)
 */
extern int inskf_uditg_dS(kf_t *inskf, double dS, const v3_t *lever_arm_car,
                             const v3_t *wib_b)
{
    if(!cfg.is_odincre)
        LOG_FATAL("%s: cfg.is_odincre = false, conflicted", __FUNCTION__);

    if(cfg.isx_kod)
        dS *= inskf->sol->kod;  /* dS scalar factor correct */

    v3_t dSec_e;
    m3_t Cce = m3_mul(inskf->sol->dcm, m3_T(inskf->sol->Cbc));
    dSec_e  = m3_mul_v3(Cce, (v3_t){dS, 0.0, 0.0});
    if(lever_arm_car != NULL && !v3_equal(lever_arm_car, &V0, 1e-4)){
        if(inskf->odt == 0.0)
            LOG_FATAL("%s: inskf->odt = 0.0, unexpected value", __FUNCTION__);
        /* NOTE: lever arm car need to do more test */
        v3_t v = m3_mul_v3(inskf->sol->dcm, v3_cross(*wib_b, *lever_arm_car));
        dSec_e = v3_del(dSec_e, v3_scalar(inskf->odt, v));
    }

    inskf->itg[cfg.ITGVEL_OD  ] += dSec_e.x;
    inskf->itg[cfg.ITGVEL_OD+1] += dSec_e.y;
    inskf->itg[cfg.ITGVEL_OD+2] += dSec_e.z;

    if(cfg.isx_kod){
        inskf->itg[cfg.ITGC_KOD  ] += inskf->sol->dcm.m11 * dS;
        inskf->itg[cfg.ITGC_KOD+1] += inskf->sol->dcm.m21 * dS;
        inskf->itg[cfg.ITGC_KOD+2] += inskf->sol->dcm.m31 * dS;
    }
    if(cfg.isx_epitch){
        inskf->itg[cfg.ITGC_EPITCH  ] += inskf->sol->dcm.m13 * dS;
        inskf->itg[cfg.ITGC_EPITCH+1] += inskf->sol->dcm.m23 * dS;
        inskf->itg[cfg.ITGC_EPITCH+2] += inskf->sol->dcm.m33 * dS;
    }
    if(cfg.isx_eyaw){
        inskf->itg[cfg.ITGC_EYAW  ] -= inskf->sol->dcm.m12 * dS;
        inskf->itg[cfg.ITGC_EYAW+1] -= inskf->sol->dcm.m22 * dS;
        inskf->itg[cfg.ITGC_EYAW+2] -= inskf->sol->dcm.m32 * dS;
    }
    return 0;
}

/**
 * @brief ins kalman filter, use odometer's integral variables as measument
 * @param[in,out]   inskf   ins kalman filter struct
 * @param[in]       itgdS_std   odometer's velocity intetral's stanadard error.
 * @return status(0: OK)
 */
extern int inskf_udmeas_itgdS(kf_t *inskf, const v3_t *itgdS_std)
{
    if(!cfg.is_odincre){
        LOG_FATAL("%s: cfg.is_odincre = false, confilcted", __FUNCTION__);
    }
    bool is_udmeas = true;

    KF_HINIT(inskf, 3);

    double dt = timediff(inskf->time, inskf->itg_start);
    v3_t tmp;
    tmp.x = - inskf->itg[cfg.ITGF_INS  ]*dt/2.0 + inskf->itg[cfg.ITGVEL_OD  ];
    tmp.y = - inskf->itg[cfg.ITGF_INS+1]*dt/2.0 + inskf->itg[cfg.ITGVEL_OD+1];
    tmp.z = - inskf->itg[cfg.ITGF_INS+2]*dt/2.0 + inskf->itg[cfg.ITGVEL_OD+2];
    m3_t mat =  m3_T(v3_askew(tmp));
    m3_paste((inskf->H+cfg.IATT), inskf->nx, &mat);

    inskf->H[cfg.IVEL] = dt;
    inskf->H[cfg.IVEL+1 + 1*inskf->nx] = dt;
    inskf->H[cfg.IVEL+2 + 2*inskf->nx] = dt;

    if(cfg.isx_kod){
        inskf->H[cfg.IKOD] = - inskf->itg[cfg.ITGC_KOD];
        inskf->H[cfg.IKOD + 1*inskf->nx] = - inskf->itg[cfg.ITGC_KOD+1];
        inskf->H[cfg.IKOD + 2*inskf->nx] = - inskf->itg[cfg.ITGC_KOD+2];
    }
    if(cfg.isx_epitch){
        inskf->H[cfg.IEPITCH] = - inskf->itg[cfg.ITGC_EPITCH];
        inskf->H[cfg.IEPITCH + 1*inskf->nx] = - inskf->itg[cfg.ITGC_EPITCH+1];
        inskf->H[cfg.IEPITCH + 2*inskf->nx] = - inskf->itg[cfg.ITGC_EPITCH+2];
    }
    if(cfg.isx_eyaw){
        inskf->H[cfg.IEYAW] = - inskf->itg[cfg.ITGC_EYAW];
        inskf->H[cfg.IEYAW + 1*inskf->nx] = - inskf->itg[cfg.ITGC_EYAW+1];
        inskf->H[cfg.IEYAW + 2*inskf->nx] = - inskf->itg[cfg.ITGC_EYAW+2];
    }

    double dz[3];
    for(int i = 0; i < 3; ++i)
        dz[i] = inskf->itg[cfg.ITGVEL_INS+i] - inskf->itg[cfg.ITGVEL_OD+i];
    LOG_TRACE("dz_dSx %f, dz_dSy %f, dz_dSz %f", dz[0], dz[1], dz[2]);

    LOG_DEBUG("dz %f %f %f (%f, %f)", dz[0], dz[1], dz[2],
           sqrt(SQR(dz[0])+SQR(dz[1])+SQR(dz[2])), inskf->sol->kod);
    m3_t Qdz = O3;
    if(cfg.iskf_itgdS_sagehusa && norm(inskf->itg+cfg.ITGVEL_OD, 3) > 1e-4){
        /* sage-husa auto filter */
        double sagehusa_itgdS_factor = 0.1;
        double sagehusa_itgdS_maxR = SQR(10.0);
        double sagehusa_itgdS_minR = SQR(5e-3);

        m3_t rho = v3_mul_cxr(*(v3_t *)dz,*(v3_t *) dz);
        double *tmpm = zeros(3, inskf->nx);
        matmul("TN", 3, inskf->nx, inskf->nx, 1.0, inskf->H, inskf->P, 0.0, tmpm);
        matmul("NN", 3, 3, inskf->nx, -1.0, tmpm, inskf->H, 1.0, (double *)&rho);
        free(tmpm);

        double *prho = (double *)&rho; double *pQdz = (double *)&Qdz;
        for(int i = 0; i < 3; ++i){
            /* limit the R range */
            if(prho[i+i*3] < sagehusa_itgdS_minR){
                pQdz[i+i*3] =
                    (1.0 - sagehusa_itgdS_factor) * inskf->R[cfg.IR_itgdS+i]
                    + sagehusa_itgdS_factor * sagehusa_itgdS_minR;
            }else if(prho[i+i*3] > sagehusa_itgdS_maxR){
                LOG_WARN("%s: sagehusa reach the max R", __FUNCTION__);
                pQdz[i+i*3] = sagehusa_itgdS_maxR;
                is_udmeas = false;
            }else{
                pQdz[i+i*3] =
                    (1.0 - sagehusa_itgdS_factor) * inskf->R[cfg.IR_itgdS+i]
                    + sagehusa_itgdS_factor * prho[i+i*3];
            }
        }
        if(is_udmeas){
            for(int i = 0; i < 3; ++i) inskf->R[cfg.IR_itgdS+i] = pQdz[i+i*3];
        }
    }else if(norm(inskf->itg+cfg.ITGVEL_OD, 3) < 1e-4){  /* Zero speed */
        Qdz = m3_scalar(SQR(1.58e-3), I3);
    }else{
        Qdz = v3_diag(v3_pow(*itgdS_std,2.0));
        Qdz = m3_mul(inskf->sol->dcm, m3_mul(Qdz, m3_T(inskf->sol->dcm)));
    }
    LOG_TRACE("Rx1 %f, Rx2 %f, Rx3 %f", sqrt(inskf->R[cfg.IR_itgdS]),
            sqrt(inskf->R[cfg.IR_itgdS+1]),
            sqrt(inskf->R[cfg.IR_itgdS+2]));

    m3_t QQ = Qdz;
    double *tmpm = zeros(3, inskf->nx);
    /** QQ =  H'*P*H + R*/
    matmul("TN", 3, inskf->nx, inskf->nx, 1.0, inskf->H, inskf->P, 0.0, tmpm);
    matmul("NN", 3, 3, inskf->nx, 1.0, tmpm, inskf->H, 1.0, (double *)&QQ);
    LOG_TRACE("std_odx %f, std_ody %f, std_odz %f",
              sqrt(QQ.m11), sqrt(QQ.m22), sqrt(QQ.m33));
    free(tmpm);

    if(is_udmeas) KF_FILTER(inskf, dz, (const double *)&Qdz);

    /* reset intergral variables */
    for(int i = 0; i < 3; ++i){
        inskf->itg[cfg.ITGVEL_INS+i]  = 0.0;
        inskf->itg[cfg.ITGF_INS+i]    = 0.0;
        inskf->itg[cfg.ITGVEL_OD+i]   = 0.0;
        inskf->itg[cfg.ITGC_KOD+i]    = 0.0;
        inskf->itg[cfg.ITGC_EYAW+i]   = 0.0;
        inskf->itg[cfg.ITGC_EPITCH+i] = 0.0;
    }
    inskf->itg_start = inskf->time;
    return 0;
}

/**
 * @brief ins kalman filter, use odometer output velocity as measurement
 * @param[in,out]   inskf       ins kalman filter struct
 * @param[in]       vod         odometer output velocity[m/s]
 * @param[in]       vod_std     odometer output velocity uncertainty[m/s]
 * @return status(0: OK)
 * @note vod_std.x is vod uncertainty, vod_std.y and vod_std.z is the Zero
 *      velocity noise
 */
int inskf_udmeas_vod(kf_t *inskf, double vod, const v3_t *vod_std)
{
    KF_HINIT(inskf, 3);
    double vod_correct = vod * inskf->sol->kod;
    m3_t Cce = m3_mul(inskf->sol->dcm, m3_T(inskf->sol->Cbc));
    v3_t vec_e  = m3_mul_v3(Cce, (v3_t){vod_correct, 0.0, 0.0});

    m3_t mat = m3_T(v3_askew(v3_scalar(-1.0, vec_e)));
    m3_paste((inskf->H+cfg.IATT), inskf->nx, &mat);

    inskf->H[cfg.IVEL] = 1.0;
    inskf->H[cfg.IVEL+ 1 + 1*inskf->nx] = 1.0;
    inskf->H[cfg.IVEL+ 2 + 2*inskf->nx] = 1.0;

    if(cfg.isx_kod){
        inskf->H[cfg.IKOD] =   -inskf->sol->dcm.m11 * vod_correct;
        inskf->H[cfg.IKOD + 1*inskf->nx] = - inskf->sol->dcm.m21 * vod_correct;
        inskf->H[cfg.IKOD + 2*inskf->nx] = - inskf->sol->dcm.m31 * vod_correct;
    }
    if(cfg.isx_eyaw){
        inskf->H[cfg.IEYAW] = inskf->sol->dcm.m12 * vod_correct;
        inskf->H[cfg.IEYAW + 1*inskf->nx] = inskf->sol->dcm.m22 * vod_correct;
        inskf->H[cfg.IEYAW + 2*inskf->nx] = inskf->sol->dcm.m32 * vod_correct;
    }
    if(cfg.isx_epitch){
        inskf->H[cfg.IEPITCH] = -inskf->sol->dcm.m13 * vod_correct;
        inskf->H[cfg.IEPITCH + 1*inskf->nx] =  -inskf->sol->dcm.m23*vod_correct;
        inskf->H[cfg.IEPITCH + 2*inskf->nx] =  -inskf->sol->dcm.m33*vod_correct;
    }

    v3_t dz = v3_del(inskf->sol->vel, vec_e);
    LOG_TRACE("dz_vodx %f, dz_vody %f, dz_vodz %f", dz.x, dz.y, dz.z);

    m3_t Qdz = v3_diag(v3_pow(*vod_std, 2.0));
    Qdz = m3_mul(m3_mul(inskf->sol->dcm, Qdz), m3_T(inskf->sol->dcm));

    m3_t QQ = Qdz;
    double *tmpm = zeros(3, inskf->nx);
    matmul("TN", 3, inskf->nx, inskf->nx, 1.0, inskf->H, inskf->P, 0.0, tmpm);
    matmul("NN", 3, 3, inskf->nx, 1.0, tmpm, inskf->H, 1.0, (double *)&QQ);
    LOG_TRACE("std_vodx %f, std_vody %f, std_vodz %f",
              sqrt(QQ.m11), sqrt(QQ.m22), sqrt(QQ.m33));
    free(tmpm);

    KF_FILTER(inskf, (const double *)&dz, (const double *)&Qdz);
    return 0;
}

/**
 * @brief ins kalman filter of car's motion constraint
 * @param[in,out]   inskf   ins kalman filter struct
 * @param[in]   wib_b       imu output angular rate[rad/s]
 * @param[in]   lever_arm_car car ref point under b-frame[m]
 * @param[in]   std_vy      zero velocity around y-axis under c-frame[m/s]
 * @param[in]   std_vz      zero veloctiy around z-axis under c-frame[m/s]
 * @return status(0: OK)
 */
extern int inskf_MC_car(kf_t *inskf, const v3_t * wib_b,
    const v3_t *lever_arm_car, double std_vy, double std_vz)
{
    KF_HINIT(inskf, 2);
    v3_t veb_b = m3_mul_v3(m3_T(inskf->sol->dcm), inskf->sol->vel);
    v3_t vec_b = v3_del(veb_b, v3_cross(*wib_b, *lever_arm_car));
    v3_t vec_c = m3_mul_v3(inskf->sol->Cbc, vec_b);

    double dz[2]; dz[0] = vec_c.y; dz[1] = vec_c.z;

    v3_t arm = m3_mul_v3(inskf->sol->Cbc, v3_cross(*wib_b, *lever_arm_car));
    LOG_TRACE("dz_MCy %f, dz_MCz %f", dz[0], dz[1]);
    LOG_TRACE("dz_MCarmy %f, dz_MCarmz %f", arm.y, arm.z);

    m3_t Cec    = m3_mul(inskf->sol->Cbc, m3_T(inskf->sol->dcm));
    m3_t Mphi   = m3_mul(Cec, v3_askew(inskf->sol->vel));

    inskf->H[cfg.IATT  ] = Mphi.m21;
    inskf->H[cfg.IATT+1] = Mphi.m22;
    inskf->H[cfg.IATT+2] = Mphi.m23;
    inskf->H[cfg.IATT+  cfg.nx] = Mphi.m31;
    inskf->H[cfg.IATT+1+cfg.nx] = Mphi.m32;
    inskf->H[cfg.IATT+2+cfg.nx] = Mphi.m33;

    inskf->H[cfg.IVEL  ] = Cec.m21;
    inskf->H[cfg.IVEL+1] = Cec.m22;
    inskf->H[cfg.IVEL+2] = Cec.m23;
    inskf->H[cfg.IVEL+  cfg.nx] = Cec.m31;
    inskf->H[cfg.IVEL+1+cfg.nx] = Cec.m32;
    inskf->H[cfg.IVEL+2+cfg.nx] = Cec.m33;

    if(cfg.isx_eyaw || cfg.isx_epitch){
        double veb_cx = Cec.m11*inskf->sol->vel.x + Cec.m12*inskf->sol->vel.y +
                Cec.m13*inskf->sol->vel.z;
        if(cfg.isx_epitch)  inskf->H[cfg.IEPITCH+cfg.nx] = -veb_cx;
        if(cfg.isx_eyaw)    inskf->H[cfg.IEYAW] =   veb_cx;
    }

//    matprint(inskf->H, inskf->nx, ny, 12, 6);
//    fflush(stdout);

//    std_vy = std_vy * (fabs(wib_b->z) > 0.01 ? fabs(wib_b->z)/0.01 : 1.0 );
//    LOG_INFO("std_vy: %f, ratio %f", std_vy, fabs(wib_b->z)/0.01);

    double Qdz[4] = {SQR(std_vy), 0.0, 0.0, SQR(std_vz)};
    double ndz[2]; inskf_norm_innov(inskf, dz, Qdz, ndz);
    if(fabs(ndz[0]) > 10.0 || fabs(ndz[1]) > 10.0){
        LOG_INFO("%s: disabled, large innovation", __FUNCTION__);
        return 1;
    }

    KF_FILTER(inskf, dz, Qdz);

    return 0;
}

/**
 * @brief Zero speed test by postion difference
 * @param[in]   inskf       ins kalman filter struct
 * @param[in]   last_sol    last solution
 * @return 0: kinematic  1: static  -1: can not be determniated
 */
extern int inskf_ZST_pos(kf_t *inskf, const solins_t *last_sol)
{
    double dt = timediff(inskf->sol->time, last_sol->time);
    double meanv = v3_norm(v3_del(inskf->sol->pos, last_sol->pos))/dt;

    v3_t fac1 = v3_scalar(1.0/v3_norm(inskf->sol->pos), inskf->sol->pos);
    v3_t fac2 = v3_scalar(1.0/v3_norm(last_sol->pos), last_sol->pos);
    double Qpos1 = v3_mul_rxc(fac1, m3_mul_v3(inskf->sol->Qpos, fac1));
    double Qpos2 = v3_mul_rxc(fac2, m3_mul_v3(last_sol->Qpos, fac2));
    double meanv_std = sqrt(Qpos1 + Qpos2)/dt;

    LOG_TRACE("meanv %f, meanv_std %f", meanv, meanv_std);
    if(meanv < 0.04 && meanv_std < 0.5)
        return 1;
    else if(meanv > 0.04 && meanv > 3*meanv_std)
        return 0;
    else
        return -1;
}

static int imu_mean_var_dv(const kf_t *inskf, double *mean_dv, double *var_dv)
{
    /* calculate varianace of accelration */
    double *acc = malloc(inskf->nimud*sizeof(double));
    double sum_acc = 0.0;
    for(unsigned short i = 0; i < inskf->nimud; ++i){
        acc[i] = v3_norm(inskf->imud[i].accel);
        sum_acc += acc[i];
    }
    *mean_dv = sum_acc / inskf->nimud;

    sum_acc = 0.0;
    for(unsigned short i = 0; i < inskf->nimud; ++i){
        sum_acc += SQR(acc[i] - *mean_dv);
    }
    *var_dv = sum_acc / (inskf->nimud -1);
    free(acc);
    return 0;
}

/**
 * @brief ins kalman filter to Zero Speed Test by imu data
 * @param[in,out]   inskf   ins kalman filter struct
 * @param[in]   mean_dv     mean of a list imu dv output
 * @param[in]   std_dv      stanadard error of a list of imu dv output
 * @return true: static, false: kinematic
 */
extern bool inskf_ZST_imu(kf_t *inskf, double mean_dv, double std_dv)
{
    if(cfg.nZST == 0){
        LOG_WARN("%s: cfg.nZST = 0, Zero Speed Test disabled", __FUNCTION__);
        inskf->ZST_count = 0;
        return false;
    }
    if(inskf->nimud < cfg.nZST){
        LOG_WARN("%s: Do not have enough imu data", __FUNCTION__);
        inskf->ZST_count = 0;
        return false;
    }
    if(fabs(mean_dv) < 1e-32 || std_dv < 1e-32){
        inskf->ZST_count = 0;
        return false;
    }

    double mean_dv_cur, var_dv_cur;
    imu_mean_var_dv(inskf, &mean_dv_cur, &var_dv_cur);

    /* F-test */
    double Ftest =  var_dv_cur / SQR(std_dv);
    LOG_TRACE("Ftest %f", Ftest);
    /* t-test */
    if(Ftest < 20.0){
        double t = (mean_dv_cur - mean_dv) / sqrt(var_dv_cur/(cfg.nZST-1));
        LOG_TRACE("Ttest %f", t);
        if(fabs(t) < 5.0){
            inskf->ZST_count++;
//            LOG_INFO("inskf->ZST_count: %i", inskf->ZST_count);
            if(inskf->ZST_count > 5) return true; else return false;
        }

    }

    inskf->ZST_count = 0;
    return false;
}

/**
 * @brief  Zero Speed Test automatically
 * @param[in,out]   inskf   ins kalman filter struct
 * @param[in]   last_sol    last solution(to check by position difference)
 * @param[in,out]   mean_dv     mean of a list imu dv output
 * @param[in,out]   std_dv      stanadard error of a list of imu dv output
 * @return true: static, false: kinematic
 */
extern bool inskf_ZST_auto(kf_t *inskf, const solins_t *last_sol,
                              double *mean_dv, double *std_dv)
{
    double mean_dv_cur, var_dv_cur;
    double fac_mean = 1.0/100;
    double fac_var = 1.0/1000;
    LOG_TRACE("mean_dv_thres %f, std_dv_thres %f", *mean_dv, *std_dv);
    switch (inskf_ZST_pos(inskf, last_sol)) {
    case 0:     /* determinate dynamic by pos */
        // inskf_ZST_imu(inskf, *mean_dv, *std_dv);
        inskf->ZST_count = 0;
        return false;
    case 1:     /* */
        /* estimate the mean_dv and std_dv by "static" imu data */
        imu_mean_var_dv(inskf, &mean_dv_cur, &var_dv_cur);
        if(fabs(*mean_dv) < 1e-32 || fabs(*std_dv) < 1e-32){
            *mean_dv = mean_dv_cur;
            *std_dv = sqrt(var_dv_cur);
        }else{
            if(sqrt(var_dv_cur) > 10.0 * *std_dv){
                LOG_WARN("%s: var_dv_cur to large: %f", __FUNCTION__,
                         sqrt(var_dv_cur));
                inskf->ZST_count = 0;
                return false;
            }
            *mean_dv = (1.0 - fac_mean) * (*mean_dv) + fac_mean * mean_dv_cur;
            *std_dv = (1.0 - fac_var) * (*std_dv) + fac_var * sqrt(var_dv_cur);
        }
        return inskf_ZST_imu(inskf, *mean_dv, *std_dv);
    case -1:    /* Can't determinate static by position */
        return inskf_ZST_imu(inskf, *mean_dv, *std_dv);
    default:
        break;
    }
    return false;
}

/**
 * @brief inskf_norm_innov
 * @param[in]   inskf ins kalman filter struct
 * @param[in]   dz    innovation(measment - model, inskf->ny x 1)
 * @param[in]   R     measment noise matrix(inskf->ny x inskf->ny)
 * @param[out]  ndz   ouput normialized innovations(inskf->ny x 1)
 * @return  0: OK
 */
extern int
inskf_norm_innov(kf_t *inskf, const double *dz, const double *R, double *ndz)
{
    /* H'*PH */
    double *tmpm = zeros(inskf->ny, inskf->nx);
    double *Rm = zeros(inskf->ny, inskf->ny);
    matmul("TN", inskf->ny, inskf->nx, inskf->nx, 1.0, inskf->H, inskf->P, 0.0,
           tmpm);
    matmul("NN", inskf->ny, inskf->ny, inskf->nx, 1.0, tmpm, inskf->H, 0.0, Rm);

    for(int i = 0; i < inskf->ny; ++i)
        ndz[i] = dz[i] / sqrt(Rm[i+i*inskf->ny] + R[i+i*inskf->ny]);

    free(tmpm); free(Rm);
    return 0;
}
