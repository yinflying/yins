#include "yinsapp.h"
#include "yins_core/insmacro.h"

typedef struct{
    bool isa_imu;
    bool isa_pps;
    imud_t imud;

    bool isa_gnss;
    gtime_t time;
    v3_t reg_e;
    m3_t Qreg_e;
    v3_t veg_e;
    m3_t Qveg_e;
    enum SOL stat;

    bool isa_gnss_yaw2;
    double yaw2;
    double yaw2_std;

    bool isa_od;
    double dS;
} prod_t;

void yinsapp_data_reset(yinsapp_data_t * appdata)
{
    appdata->week = 0;
    appdata->sow = 0.0;

    appdata->isa_imu = false;
    appdata->isa_pps = false;
    appdata->isa_gnss = false;
    appdata->isa_od = false;

    for(int i = 0; i < 3; ++i){
        appdata->gyro[i] = 0.0;
        appdata->accel[i] = 0.0;
        appdata->pos[i] = 0.0;
        appdata->pos_std[i] = 0.0;
        appdata->veg_n[i] = 0.0;
        appdata->veg_n_std[i] = 0.0;
    }
    appdata->gnss_stat = YINS_GNSS_STAT_NONE;
    appdata->yaw2ant = 0.0;
    appdata->yaw2ant_std = 0.0;
    appdata->yaw2ant_stat = YINS_YAW2ANT_STAT_NONE;
    appdata->dS = 0.0;
}


static void prod_reset_stat(prod_t *prod)
{
    prod->isa_od = false;
    prod->isa_imu = false;
    prod->isa_pps = false;
    prod->isa_gnss = false;
    prod->isa_gnss_yaw2 = false;
}

void yinsapp_data2prod(const yinsapp_data_t *appdata, prod_t *prod)
{
    prod_reset_stat(prod);
    prod->time = gpst2time(appdata->week, appdata->sow);
    if(appdata->isa_imu){
        prod->isa_imu = true;
        prod->imud.time = prod->time;
        v3_copy(&prod->imud.gyro, appdata->gyro);
        v3_copy(&prod->imud.accel, appdata->accel);
        return;
    }
    if(appdata->isa_gnss){
        v3_copy(&prod->reg_e, appdata->pos);
        v3_copy(&prod->veg_e, appdata->veg_n);

        prod->Qreg_e = O3;
        prod->Qreg_e.m11 = appdata->pos_std[0];
        prod->Qreg_e.m22 = appdata->pos_std[1];
        prod->Qreg_e.m33 = appdata->pos_std[2];
        prod->Qveg_e = O3;
        prod->Qveg_e.m11 = appdata->veg_n_std[0];
        prod->Qveg_e.m22 = appdata->veg_n_std[1];
        prod->Qveg_e.m33 = appdata->veg_n_std[2];

        ned2ecefQ(&prod->reg_e, &prod->Qreg_e, &prod->Qveg_e, NULL);
        ned2ecef(&prod->reg_e, &prod->veg_e, NULL);

        switch (appdata->gnss_stat) {
        case YINS_GNSS_STAT_NONE: 	prod->stat = SOL_NONE;
        case YINS_GNSS_STAT_SINGLE: prod->stat = SOL_SINGLE;
        case YINS_GNSS_STAT_DGNSS: 	prod->stat = SOL_DGNSS;
        case YINS_GNSS_STAT_PPP: 	prod->stat = SOL_PPP;
        case YINS_GNSS_STAT_FLOAT: 	prod->stat = SOL_FLOAT;
        case YINS_GNSS_STAT_FIX: 	prod->stat = SOL_FIXED;
        }

        if(appdata->yaw2ant_stat != YINS_YAW2ANT_STAT_NONE){
            prod->isa_gnss_yaw2 = true;
            prod->yaw2 = appdata->yaw2ant*DEG2RAD;
            prod->yaw2_std = appdata->yaw2ant_std*DEG2RAD;
        }
    }
    if(appdata->isa_od){
        prod->isa_od = true;
        prod->dS = appdata->dS;
    }
}

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

extern int yinsapp_process(const yinsapp_data_t *appdata,
                           yinsapp_result_t *result, const char *cfg_file)
{
    /* read cfg */
    bool is_cfg = false;
    static imup_t imup = {0};
    if(!is_cfg){
        FILE *fp = fopen(cfg_file, "r");
        readf_ycsv_header(fp, &imup);
        is_cfg = true;
    }

    /* align and initial kalman filter */
    static kf_t inskf;
    static bool is_align = false;
    static int align_count = 0;
    static v3_t align_att = {0.0, 0.0, 0.0};
    if(!is_align && !appdata->isa_gnss){
        return 1;
    }
    if(!is_align && appdata->gnss_stat != YINS_GNSS_STAT_NONE){
        double cur_yaw;
        double vel = sqrt(SQR(appdata->veg_n[0])  + SQR(appdata->veg_n[1]));

        /* use gnss 2-antenna yaw or velocity */
        if(appdata->yaw2ant_stat != YINS_YAW2ANT_STAT_NONE)
            cur_yaw = appdata->yaw2ant*DEG2RAD;
        else if(vel > 2.0){
            vel2yaw((const v3_t *)appdata->veg_n, &cur_yaw, NULL, NULL);
            align_count ++;
        } else
            align_count = 0;

        /* Check align yaw result consistent */
        if(align_count > 0){
            double dyaw = angle_to180(appdata->yaw2ant*DEG2RAD - align_att.z);
            if(dyaw < 5*DEG2RAD) align_count++; else align_count = 0;
        }
        /* Accept the aligning result */
        if(align_count > 3){
            LOG_INFO("Finish align");
            is_align = true;

            imup.tstart = gpst2time(appdata->week, appdata->sow);
            v3_copy(&imup.initv, appdata->veg_n);
            v3_copy(&imup.initr, appdata->pos);
            ned2ecef(&imup.initr, &imup.initv, NULL);
            inskf_init(&inskf, &imup);
        }
    }
    if(!is_align) return 1;

    static prod_t *prod; prod = malloc(sizeof(prod_t));
    yinsapp_data2prod(appdata, prod);

    /* run kalman filter */
    if(prod->isa_imu){
        inskf_udstate(&inskf, &prod->imud, &imup);
        return  2;
    }
    if(prod->isa_gnss){
        v3_t wib_b = v3_scalar(imup.freq_imu, prod->imud.gyro);
        inskf_udmeasv(&inskf, &prod->veg_e, &prod->Qveg_e, &imup.lever_arm_gps, &wib_b);
        inskf_udmeasr(&inskf, &prod->reg_e, &prod->Qreg_e, &imup.lever_arm_gps);
        inskf_feedback(&inskf, soltype_add(inskf.sol->status, prod->stat));
        if(prod->isa_gnss_yaw2){
            inskf_udmeas_yaw(&inskf, prod->yaw2, SQR(prod->yaw2_std));
            inskf_feedback(&inskf, soltype_add(inskf.sol->status, prod->stat));
        }
        return 3;
    }
    return 0;
}
