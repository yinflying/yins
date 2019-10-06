#include "yinsapp.h"
#include <string.h>

extern int yinsapp_pureins(
        const char *fin, enum FT ft, const imup_t *imup, const char *fsol)
{
    cfg.isx_ba = 0; cfg.isx_bg = 0;
    solins_t sol; memset(&sol, 0, sizeof(solins_t));

    FILE *fp_sol;
    if(!(fp_sol = fopen(fsol, "w")))
        LOG_FATAL("Can't open file: %s", fsol);

    if(FILE_LOG_LEVEL != LEVEL_OFF) LOG_OPEN(cfg.log_path);

    imu_t imu; imu_init(&imu);
    yins_readf(fin, ft, &imu, NULL, NULL);

    v3_t reb_e = imup->initr;
    v3_t veb_e = imup->initv;
    ned2ecef(&reb_e, &veb_e, NULL);
    quat_t qbe =  att2Qbe(&imup->initr, &imup->inita);
    double dt = 1.0 / imup->freq_imu;

    /* search start time of IMU data */
    unsigned int istart = 0;
    if(imup->tstart.time == 0 && imup->tstart.sec == 0.0){
        LOG_WARN("imup->tstart may not set correctly");
    }else{
        double t = timediff(imup->tstart, imu.data[0].time);
        if(t < - 1.5*dt)
            LOG_FATAL("imup->tstart is ahead of imu.data[0].time %.2f sec", -t);
        else if(t < 0.5*dt)
            istart = 0;
        else {
            for(istart = 1 ;istart < imu.n ; ++istart)
                if(timediff(imup->tstart, imu.data[0].time) < 0.5*dt) break;
            if(istart == imu.n)
                LOG_FATAL("imup->tstart setting error");
        }
    }

    /* run ins_nav_ecef */
    sol.pos = reb_e; sol.vel = veb_e; sol.quat = qbe;
    quat2dcm(&sol.quat, &sol.dcm);
    sol.time = imup->tstart;
    outsolins(fp_sol, &sol, imup);
    for(unsigned int i = istart; i < imu.n; i++){
        ins_nav_ecef(dt, &imu.data[i].gyro, &imu.data[i].accel,
                     &reb_e, &veb_e, &qbe);

        sol.time = imu.data[i].time;
        sol.pos = reb_e; sol.vel = veb_e; sol.quat = qbe;
        quat2dcm(&sol.quat, &sol.dcm);
        outsolins(fp_sol, &sol, imup);
    }

    fclose(fp_sol);
    if(FILE_LOG_LEVEL != LEVEL_OFF) LOG_CLOSE();
    return 0;
}
