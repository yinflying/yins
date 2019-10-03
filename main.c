/**
 * @file main.c
 * @brief main function, provide CUI program
 * @author yinflying(yinflying@foxmail.com)
 * @version 0.0.1
 * @note
 *  2019-07-01  Add header comments
 */

/**
 ** This file is part of the qt_yins project.
 ** Copyright 2019 yinflying <yinflying@foxmail.com>.
 **
 ** This program is free software: you can redistribute it and/or modify
 ** it under the terms of the GNU General Public License as published by
 ** the Free Software Foundation, either version 3 of the License, or
 ** (at your option) any later version.
 **
 ** This program is distributed in the hope that it will be useful,
 ** but WITHOUT ANY WARRANTY; without even the implied warranty of
 ** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 ** GNU General Public License for more details.
 **
 ** You should have received a copy of the GNU General Public License
 ** along with this program.  If not, see <http://www.gnu.org/licenses/>.
 **/

#include "yins_core/ins.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define PI 3.14159265358979
#define SQR(x) ((x)*(x))

void test_ygm_insod(){
#define WPATH "../yins/yins_data/ygm_insod/"
#define NAME  "ygm_circle_"
    /* open a file to record solution */
    FILE *fp_sol = fopen(WPATH NAME "sol.ycsv","w");
    /* open a file to recoard log */
    LOG_OPEN(WPATH NAME "yins.log");

    /* read imu file  */
    imu_t imu;  imu_init(&imu);
    yins_readf(WPATH NAME "imu.txt",FT_YGM_IMU, &imu, NULL, NULL);
    /* read pva(position,velocity,attitude) file */
    pva_t pva;  pva_init(&pva);
    yins_readf(WPATH NAME "avp.txt",FT_YGM_AVP, NULL, &pva, NULL);
    /* read od(odometer) file */
    od_t  od;   od_init(&od);
    yins_readf(WPATH NAME "od.txt", FT_YGM_OD, NULL, NULL, &od);
    /* read configure file */
    yins_readf(WPATH NAME "cfg.ycsv", FT_CFG_YCSV, &imu, NULL, NULL);

    /* set kalman filter initial time and postion message */
    imu.property->tstart = pva.time[0];
    imu.property->initr = pva.pos[0];
    imu.property->initv = pva.vel[0];
    imu.property->inita = pva.att[0];
    v3_t initpos = imu.property->initr; ecef2ned(&initpos, NULL, NULL);
    imu.property->inita = att2Ebe(&initpos, &imu.property->inita);

    // fix var-covarinace setting
    m3_t Qr = m3_scalar(SQR(0.10), I3);
    m3_t Qv = m3_scalar(SQR(0.15), I3);
    v3_t lever_arm = V0;

    kf_t inskf;
    v3_t pos, vel;
    inskf_init(&inskf, imu.property);
    outsolins(fp_sol, inskf.sol, imu.property);
    for(unsigned int i = 1; i < imu.n; ++i){
        /* update odometer increment */
        inskf_uditg_dS(&inskf, od.dS[i], NULL, NULL);
        /* update kalman statue */
        inskf_udstate(&inskf, (imu.data+i), imu.property);
        /* update GNSS position and velocity every 0.2s */
        if(i % 20 == 0 && i != 0){
            if(fabs(timediff(pva.time[i], inskf.sol->time)) > 0.1){
                LOG_FATAL("pva.time differ from  inskf.sol->time");
            }
            /* measurement update of odometer increment */
            v3_t QitgdS = (v3_t){1e-2, 1e-2, 1e-2};
            inskf_udmeas_itgdS(&inskf, &QitgdS);
            inskf_feedback(&inskf, soltype_add(inskf.sol->status, SOL_DR));

            /* add random error to postion and velocity */
            pos = pva.pos[i]; vel = pva.vel[i];
            pos = v3_add(pos, v3_randn(0.0, 0.05));
            vel = v3_add(vel, v3_randn(0.0, 0.10));

            /* simulate GNSS lost after 300s */
            if(i < 100*300){
                /* measure update of position */
                inskf_udmeasr(&inskf, &pos, &Qr, &lever_arm);
                v3_t wib_b = v3_scalar(imu.property->freq_imu, imu.data[i].gyro);
                /* measure update of velocity */
                inskf_udmeasv(&inskf, &vel, &Qv, &lever_arm, &wib_b);
                inskf_feedback(&inskf, soltype_add(inskf.sol->status, SOL_SINGLE));
            }
        }
        /* output solution to file */
        outsolins(fp_sol, inskf.sol, imu.property);
    }
    fclose(fp_sol);
    imu_free(&imu); od_free(&od); pva_free(&pva); LOG_CLOSE();
#undef WPATH
#undef NAME
}


int main()
{
    FILE_LOG_LEVEL = LEVEL_TRACE;
    STDOUT_LOG_LEVEL = LEVEL_DEBUG;

    test_ygm_insod();

    LOG_CLOSE();
    return 0;
}
