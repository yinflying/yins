/**
 * @file insio.c
 * @brief ins read and write function
 * @author yinflying(yinflying@foxmail.com)
 * @note
 *  2019-05-21 Created \n
 *  2019-09-03 Add ycsv read and output support
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
#include "insmacro.h"

#define MAXLINELEN  512     /**< char numbner limit of line when reading file */
#define MAXIMUOBS   100000  /**< max memory allocate for imu_t struct */
#define MAXODOBS    100000  /**< max memory allocate for od_t sturct */
#define MAXPVAOBS   100000  /**< max memory allocate for pva_t struct */

/** safe to call realloc */
#define REALLOC(pointer, type, sz)                                              \
    if(!((pointer) =                                                            \
        (type *)realloc((pointer), ((unsigned long)sz)*sizeof(type)))){         \
        LOG_FATAL(#pointer " memeory realloc failed");                          \
    }

/** safe to call malloc */
#define MALLOC(pointer, type, sz)                                               \
    if(!((pointer) = (type *)malloc(((unsigned long)sz)*sizeof(type)))){        \
        LOG_FATAL(#pointer " memory allocation failed");                        \
    }

/** safe to free pointer */
#define FREE(pointer) {free(pointer);(pointer) = NULL;}

/* string to number ------------------------------------------------------------
 * convert substring in string to number
 * args   : char   *s        I   string ("... nnn.nnn ...")
 *          int    i,n       I   substring position and width
 * return : converted number (0.0:error)
 *-----------------------------------------------------------------------------*/
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

static void fprintf_fixwidth(FILE* fp, double num, int fixwidth)
{
    char buf[32];
    sprintf(buf, "%.*f", fixwidth, num);
    buf[fixwidth] = '\0';
    fprintf(fp, "%s", buf);
}

static int writef_imup(FILE *fp, const imup_t *imup)
{
    fprintf(fp, "imup.freq_imu:   %12i  #%s\n", imup->freq_imu,
            "[Hz]");
    fprintf(fp, "imup.freq_od:    %12i  #%s\n", imup->freq_od,
            "[Hz]");
    fprintf(fp, "imup.accel_noise:%12.6G, %12.6G, %12.6G #%s\n",
        imup->accel_noise.x, imup->accel_noise.y, imup->accel_noise.z,
            "[m/s^2]");
    fprintf(fp, "imup.gyro_noise: %12.6G, %12.6G, %12.6G #%s\n",
        imup->gyro_noise.x*RAD2DEG, imup->gyro_noise.y*RAD2DEG,
        imup->gyro_noise.z*RAD2DEG,
            "[deg/s]");
    fprintf(fp, "imup.ba:         %12.6G, %12.6G, %12.6G #%s\n",
        imup->ba.x*MPS22MG, imup->ba.y*MPS22MG, imup->ba.z*MPS22MG,
            "[mg]");
    fprintf(fp, "imup.ba_std:     %12.6G, %12.6G, %12.6G #%s\n",
        imup->ba_std.x*MPS22MG, imup->ba_std.y*MPS22MG, imup->ba_std.z*MPS22MG,
            "[mg]");
    fprintf(fp, "imup.bg:         %12.6G, %12.6G, %12.6G #%s\n",
        imup->bg.x*RAD2DEG, imup->bg.y*RAD2DEG, imup->bg.z*RAD2DEG,
            "[deg/s]");
    fprintf(fp, "imup.bg_std:     %12.6G, %12.6G, %12.6G #%s\n",
        imup->bg_std.x*RAD2DEG, imup->bg_std.y*RAD2DEG, imup->bg_std.z*RAD2DEG,
            "[deg/s]");
    fprintf(fp, "imup.ka:         %12.6f, %12.6f, %12.6f #%s\n",
        imup->ka.x, imup->ka.y, imup->ka.z,
            "[]");
    fprintf(fp, "imup.ka_std:     %12.4G, %12.4G, %12.4G #%s\n",
        imup->ka_std.x, imup->ka_std.y, imup->ka_std.z,
            "[]");
    fprintf(fp, "imup.kg:         %12.6f, %12.6f, %12.6f #%s\n",
        imup->kg.x, imup->kg.y, imup->kg.z,
            "[]");
    fprintf(fp, "imup.kg_std:     %12.4G, %12.4G, %12.4G #%s\n",
        imup->kg_std.x, imup->kg_std.y, imup->kg_std.z,
            "[]");
    fprintf(fp, "imup.arw:        %12.4G, %12.4G, %12.4G #%s\n",
        imup->arw.x*RPSS2DPSH, imup->arw.y*RPSS2DPSH, imup->arw.z*RPSS2DPSH,
            "[deg/sqrt(h)]");
    fprintf(fp, "imup.arrw:       %12.4G, %12.4G, %12.4G #%s\n",
        imup->arrw.x, imup->arrw.y, imup->arrw.z,
            "[rad/s/sqrt(s)]");
    fprintf(fp, "imup.vrw:        %12.4G, %12.4G, %12.4G #%s\n",
        imup->vrw.x*MPS22MG, imup->vrw.y*MPS22MG, imup->vrw.z*MPS22MG,
            "[mg/sqrt(Hz)]");
    fprintf(fp, "imup.vrrw:       %12.4G, %12.4G, %12.4G #%s\n",
        imup->vrrw.x, imup->vrrw.y, imup->vrrw.z,
            "[m/s^2/sqrt(s)]");
    fprintf(fp, "imup.Ta:         %12.2f, %12.2f, %12.2f #%s\n",
            imup->Ta.x,  imup->Ta.y, imup->Ta.z,
            "[s]");
    fprintf(fp, "imup.Tg:         %12.2f, %12.2f, %12.2f #%s\n",
            imup->Tg.x,  imup->Tg.y, imup->Tg.z,
            "[s]");
    fprintf(fp, "imup.kod:        %12.4f #%s\n", imup->kod,
            "[]");
    fprintf(fp, "imup.kod_std:    %12.4f #%s\n", imup->kod_std,
            "[]");
    fprintf(fp, "imup.lever_arm_gps:     %12.3f, %12.3f, %12.3f #%s\n",
        imup->lever_arm_gps.x, imup->lever_arm_gps.y, imup->lever_arm_gps.z,
            "[m]");
    fprintf(fp, "imup.lever_arm_gps_std: %12.3f, %12.3f, %12.3f #%s\n",
            imup->lever_arm_gps_std.x,
            imup->lever_arm_gps_std.y,
            imup->lever_arm_gps_std.z,
            "[m]");
    fprintf(fp, "imup.lever_arm_od:      %12.3f, %12.3f, %12.3f #%s\n",
        imup->lever_arm_od.x, imup->lever_arm_od.y, imup->lever_arm_od.z,
            "[m]");
    fprintf(fp, "imup.lever_arm_od_std:  %12.3f, %12.3f, %12.3f #%s\n",
            imup->lever_arm_od_std.x,
            imup->lever_arm_od_std.y,
            imup->lever_arm_od_std.z,
            "[m]");
    fprintf(fp, "imup.lever_arm_car:     %12.3f, %12.3f, %12.3f #%s\n",
        imup->lever_arm_car.x, imup->lever_arm_car.y, imup->lever_arm_car.z,
            "[m]");
    fprintf(fp, "imup.lever_arm_car_std: %12.3f, %12.3f, %12.3f #%s\n",
            imup->lever_arm_car_std.x,
            imup->lever_arm_car_std.y,
            imup->lever_arm_car_std.z,
            "[m]");
    fprintf(fp, "imup.err_angle_imu:     %12.3f, %12.3f, %12.3f #%s\n",
            imup->err_angle_imu.x*RAD2DEG,
            imup->err_angle_imu.y*RAD2DEG,
            imup->err_angle_imu.z*RAD2DEG,
            "[deg]");
    fprintf(fp, "imup.err_angle_imu_std: %12.3f, %12.3f, %12.3f #%s\n",
            imup->err_angle_imu_std.x*RAD2DEG,
            imup->err_angle_imu_std.y*RAD2DEG,
            imup->err_angle_imu_std.z*RAD2DEG,
            "[deg]");
    fprintf(fp, "imup.err_angle_imu_rw:  %12.3G, %12.3G, %12.3G #%s\n",
            imup->err_angle_imu_rw.x,
            imup->err_angle_imu_rw.y,
            imup->err_angle_imu_rw.z,
            "[rad/sqrt(s)] IMU install error angle random walk");
    fprintf(fp, "imup.Terr_angle_imu:    %12.3f, %12.3f, %12.3f #%s\n",
            imup->Terr_angle_imu.x,
            imup->Terr_angle_imu.y,
            imup->Terr_angle_imu.z,
            "[s] IMU install error angle correlcation time(1st order Markov)");
    fprintf(fp, "imup.err_angle_gps:     %12.3f, %12.3f, %12.3f #%s\n",
            imup->err_angle_gps.x*RAD2DEG,
            imup->err_angle_gps.y*RAD2DEG,
            imup->err_angle_gps.z*RAD2DEG,
            "[deg]");
    fprintf(fp, "imup.err_angle_gps_std: %12.3f, %12.3f, %12.3f #%s\n",
            imup->err_angle_gps_std.x*RAD2DEG,
            imup->err_angle_gps_std.y*RAD2DEG,
            imup->err_angle_gps_std.z*RAD2DEG,
            "[deg]");
    fprintf(fp, "imup.ref_point:         %12.3f, %12.3f, %12.3f #%s\n",
        imup->ref_point.x, imup->ref_point.y, imup->ref_point.z,
            "[m]");
    return 0;
}

static int writef_cfg(FILE *fp){
    fprintf(fp, "cfg.isx_ba:    %i %16s#%s\n", cfg.isx_ba, "",
            "[bool] if acceleremter bias in kalman filter state vector or not");
    fprintf(fp, "cfg.isx_bg:    %i %16s#%s\n", cfg.isx_bg, "",
            "[bool] if gyro bias in kalman filter state vector or not");
    fprintf(fp, "cfg.isx_kax:   %i %16s#%s\n", cfg.isx_kax, "",
            "[bool]");
    fprintf(fp, "cfg.isx_kay:   %i %16s#%s\n", cfg.isx_kay, "",
            "[bool]");
    fprintf(fp, "cfg.isx_kaz:   %i %16s#%s\n", cfg.isx_kaz, "",
            "[bool]");
    fprintf(fp, "cfg.isx_kgx:   %i %16s#%s\n", cfg.isx_kgx, "",
            "[bool]");
    fprintf(fp, "cfg.isx_kgy:   %i %16s#%s\n", cfg.isx_kgy, "",
            "[bool]");
    fprintf(fp, "cfg.isx_kgz:   %i %16s#%s\n", cfg.isx_kgz, "",
            "[bool]");
    fprintf(fp, "cfg.isx_kod:   %i %16s#%s\n", cfg.isx_kod, "",
            "[bool]");
    fprintf(fp, "cfg.isx_eroll: %i %16s#%s\n", cfg.isx_eroll, "",
            "[bool]");
    fprintf(fp, "cfg.isx_epitch:    %i %12s#%s\n", cfg.isx_epitch, "",
            "[bool]");
    fprintf(fp, "cfg.isx_eyaw:      %i %12s#%s\n", cfg.isx_eyaw, "",
            "[bool]");
    fprintf(fp, "cfg.isx_armgps:    %i %12s#%s\n", cfg.isx_armgps,"",
            "[bool]");
    fprintf(fp, "cfg.isx_armcar:    %i %12s#%s\n", cfg.isx_armcar, "",
            "[bool]");
    fprintf(fp, "cfg.is_odincre:    %i %12s#%s\n", cfg.is_odincre, "",
            "[bool]");
    fprintf(fp, "cfg.iskf_itgdS_sagehusa:   %i %4s#%s\n",
            cfg.iskf_itgdS_sagehusa ,"",
            "[bool]");
    fprintf(fp, "cfg.max_ny:        %3i %10s#%s\n", cfg.max_ny, "",
            "[int]");
    fprintf(fp, "cfg.nZST:          %3i %10s#%s\n", cfg.nZST, "",
            "[int]");
    fprintf(fp, "cfg.issol_header:  %i %12s#%s\n", cfg.issol_header, "",
            "[bool]");
    fprintf(fp, "cfg.sol_refpos:    %i %12s#%s\n", cfg.sol_refpos, "",
            "[bool]");
    fprintf(fp, "cfg.feedratio:     %4.2f %9s#%s\n", (double)cfg.feedratio,"",
            "[float]");
    return 0;
}

/* read Novatel ascii file */
static int readf_imu_nvt(FILE* fp, imu_t* imu)
{
    if (!(imu->data = (imud_t*)malloc(sizeof(imud_t) * MAXIMUOBS))){
        LOG_FATAL("imu->data memory alloction error");
    }
    imu->n = 0;
    imu->nmax = MAXIMUOBS;

    char buff[MAXLINELEN];
    imud_t data;
    char seps[] = ",;*";
    char* token;
    int gpsw;

    double a_scale = 1.52587890625E-06;     /* 0.05*2^-15 */
    double g_scale = 1.085069444444444E-07; /* 0.1/(3600*256) */
    while (fgets(buff, MAXLINELEN, fp)) {
        token = strtok(buff, seps);
        if (!strcmp(token, "%RAWIMUSA")) {
            for (int i = 0; i < 3; ++i)
                token = strtok(NULL, seps); /* skip header */
            gpsw = atoi(token);
            token = strtok(NULL, seps);
            data.time = gpst2time(gpsw, atof(token));
            for (int i = 0; i < 2; ++i)
                token = strtok(NULL, seps); /* skip flag */
            /* velocity increment, m/s */
            data.accel.z = -atof(token) * a_scale;
            token = strtok(NULL, seps);
            data.accel.x = -atof(token) * a_scale;
            token = strtok(NULL, seps);
            data.accel.y = atof(token) * a_scale;
            token = strtok(NULL, seps);
            /* angle increment, rad */
            data.gyro.z = -atof(token) * g_scale;
            token = strtok(NULL, seps);
            data.gyro.x = -atof(token) * g_scale;
            token = strtok(NULL, seps);
            data.gyro.y = atof(token) * g_scale;
            imu_add(imu, &data);
        }
    }
    return 0;
}

static int readf_imu_csv(FILE *fp, imu_t *imu)
{
    char buff[MAXLINELEN];
    imud_t data;
    double week, sow;
    while (fgets(buff, MAXLINELEN, fp)) {
        sscanf(buff, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf",
            &week, &sow, &data.accel.x, &data.accel.y, &data.accel.z,
            &data.gyro.x, &data.gyro.y, &data.gyro.z);
        data.time = gpst2time((int)week,sow);
        imu_add(imu, &data);
    }
    return 0;
}

static int readf_ygm_imu(FILE *fp, imu_t *imu)
{
    char buff[MAXLINELEN];
    imud_t imud;
    double gx, gy, gz, ax, ay, az, t;
    while (fgets(buff, MAXLINELEN, fp)) {
        sscanf(buff, "%lf %lf %lf %lf %lf %lf %lf",
               &gx, &gy, &gz, &ax, &ay, &az, &t);
        imud.time.time = 0;
        imud.time.sec = t;
        imud.gyro = (v3_t){gy, gx, -gz};
        imud.accel = (v3_t){ay, ax, -az};
        imu_add(imu, &imud);
    }
    return 0;
}

static int readf_ygm_avp(FILE *fp, pva_t *pva)
{
    char buff[MAXLINELEN];
    double lat, lon, hgt;
    double vE, vN, vU;
    double pitch, roll, yaw;
    double t;
    while (fgets(buff, MAXLINELEN, fp)) {
        sscanf(buff, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
               &pitch, &roll, &yaw, &vE, &vN, &vU, &lat, &lon, &hgt, &t);

        if(pva->n > pva->nmax){
            pva->nmax *= 2;
            pva_resize(pva, pva->nmax);
        }

        pva->time[pva->n].time = 0;
        pva->time[pva->n].sec = t;
        pva->pos[pva->n] = (v3_t){lat, lon, hgt};
        pva->vel[pva->n] = (v3_t){vN, vE, -vU};
        pva->att[pva->n] = (v3_t){roll, pitch, -yaw};
        ned2ecef(&pva->pos[pva->n], &pva->vel[pva->n], NULL);
        pva->status[pva->n] = SOL_FIXED;
        pva->n ++;
    }
    return 0;
}

static int readf_ygm_od(FILE *fp, od_t *od)
{
    char buff[MAXLINELEN];
    double dS, t;
    gtime_t time; time.time = 0.0;
    while(fgets(buff, MAXLINELEN, fp)){
        sscanf(buff, "%lf %lf", &dS, &t);
        time.sec = t;
        od_add(od, &time, &dS);
    }
    return 0;
}

static int reads_nmea_GPGGA(char *buff, double *sec, v3_t *pos,
                      enum SOL *status)
{
    char seps[] = ",*";
    char *token = strtok(buff, seps);
    int deg; double arc;
    if(!strcmp(token, "$GPGGA")){
        if((token = strtok(NULL, seps))){
            int hh,mm; sscanf(token, "%2i%2i%lf", &hh, &mm, sec);
            *sec = hh*3600.0 + mm * 60.0 + *sec;
        }else{ return -1; }
        if((token = strtok(NULL, seps))){
            sscanf(token, "%2i%lf", &deg, &arc);
            if(deg < 0.0 || deg > 90.0 || arc < 0.0 || arc > 60.0)
                LOG_ERROR("%s: parse error", __FUNCTION__);
            pos->x = (deg + arc/60.0)*DEG2RAD;
        }else{ return -1; }
        if((token = strtok(NULL, seps))){
            if(token[0] == 'S') pos->x = - pos->x;
        }else{ return -1; }
        if((token = strtok(NULL, seps))){
            sscanf(token, "%3i%lf", &deg, &arc);
            if(deg < 0.0 || deg > 180.0 || arc < 0.0 || arc > 60.0){
                LOG_ERROR("%s: parse error", __FUNCTION__);
            }
            pos->y = (deg + arc/60.0)*DEG2RAD;
        }else{ return -1; }
        if((token = strtok(NULL, seps))){
            if(token[0] == 'W') pos->y = - pos->y;
        }
        if((token = strtok(NULL, seps)) && status != NULL){
            switch(*token) {
            case '0': *status = SOL_NONE;   break;      /* invalid */
            case '1': *status = SOL_SINGLE; break;      /* SPS fix */
            case '2': *status = SOL_DGNSS;  break;      /* DGPS fix */
            case '3': *status = SOL_SINGLE; break;      /* PPS fix */
            case '4': *status = SOL_FIXED;  break;
            case '5': *status = SOL_FLOAT;  break;
            case '6': *status = SOL_DR;     break;
            case '7': *status = SOL_MANUAL; break;
            default:
                LOG_ERROR("%s: Can not recognize GPS staus: %s", __FUNCTION__,
                          token);
                return -1;
            }
        }
        if((token = strtok(NULL, seps))){}      /* number of satellites */
        if((token = strtok(NULL, seps))){}      /* HDOP */
        if((token = strtok(NULL, seps))){
            pos->z = atof(token);               /* altitude(above sea level) */
            if((token = strtok(NULL, seps)) && strcmp(token, "M"))
                LOG_ERROR("%s: unit of altitude can not be identified",
                         __FUNCTION__);
        }
        if((token = strtok(NULL, seps))){  /* sea level above WGS84 ellipsoid */
            pos->z += atof(token);
            if((token = strtok(NULL, seps)) && strcmp(token, "M"))
                LOG_ERROR("%s: unit of sea level height can not be identified",
                          __FUNCTION__);
        }
        return  0;
    }else{ return -1; }
}

static int readf_nmea(FILE *fp, pva_t *pva)
{
    char buff[MAXLINELEN];
    while(fgets(buff, MAXLINELEN, fp)){
        if(pva->n >= pva->nmax){
            pva->nmax *= 2;
            pva_resize(pva, pva->nmax);
        }
        if(!reads_nmea_GPGGA(buff, &pva->time[pva->n].sec, &pva->pos[pva->n],
                             &pva->status[pva->n]))
            ned2ecef(&pva->pos[pva->n], &pva->vel[pva->n], NULL);
            pva->n++;
    }
    return 0;
}

inline static void reads_ycsv_header_convert(char *token, v3_t *imup){
    if((token = strtok(NULL, ","))) imup->x = atof(token);
    if((token = strtok(NULL, ","))) imup->y = atof(token);
    if((token = strtok(NULL, ","))) imup->z = atof(token);
}

static int reads_ycsv_header(char *buff, imup_t *imup)
{
    /* remove comment part from buff(after #) */
    for(unsigned int i = 0; i < strlen(buff); ++i){
        if(buff[i] == '#') buff[i] = '\0';
    }
    char *token = strtok(buff, ":");
    if(imup != NULL){
        if(!strcmp(token, "imup.freq_imu")){
            if((token = strtok(NULL, ",")))
                imup->freq_imu = (unsigned int)atoi(token);
        }else if(!strcmp(token, "imup.freq_od")){
            if((token = strtok(NULL, ",")))
                imup->freq_od = (unsigned int)atoi(token);
        }else if(!strcmp(token, "imup.accel_noise")){
            reads_ycsv_header_convert(token, &imup->accel_noise);
        }else if(!strcmp(token, "imup.gyro_noise")){    /* deg/s */
            reads_ycsv_header_convert(token, &imup->gyro_noise);
            imup->gyro_noise = v3_scalar(DEG2RAD, imup->gyro_noise);
        }else if(!strcmp(token, "imup.ba")){            /* mg */
            reads_ycsv_header_convert(token, &imup->ba);
            imup->ba = v3_scalar(MG2MPS2, imup->ba);
        }else if(!strcmp(token, "imup.ba_std")){        /* mg */
            reads_ycsv_header_convert(token, &imup->ba_std);
            imup->ba_std = v3_scalar(MG2MPS2, imup->ba_std);
        }else if(!strcmp(token, "imup.bg")){            /* deg/s */
            reads_ycsv_header_convert(token, &imup->bg);
            imup->bg = v3_scalar(DEG2RAD, imup->bg);
        }else if(!strcmp(token, "imup.bg_std")){        /* deg/s */
            reads_ycsv_header_convert(token, &imup->bg_std);
            imup->bg_std = v3_scalar(DEG2RAD, imup->bg_std);
        }else if(!strcmp(token, "imup.ka")){
            reads_ycsv_header_convert(token, &imup->ka);
        }else if(!strcmp(token, "imup.ka_std")){
            reads_ycsv_header_convert(token, &imup->ka_std);
        }else if(!strcmp(token, "imup.kg")){
            reads_ycsv_header_convert(token, &imup->kg);
        }else if(!strcmp(token, "imup.kg_std")){
            reads_ycsv_header_convert(token, &imup->kg_std);
        }else if(!strcmp(token, "imup.arw")){        /* deg/sqrt(h) */
            reads_ycsv_header_convert(token, &imup->arw);
            imup->arw = v3_scalar(DPSH2RPSS, imup->arw);
        }else if(!strcmp(token, "imup.arrw")){
            reads_ycsv_header_convert(token, &imup->arrw);
        }else if(!strcmp(token, "imup.vrw")){       /* mg/sqrt(Hz) */
            reads_ycsv_header_convert(token, &imup->vrw);
            imup->vrw = v3_scalar(MG2MPS2, imup->vrw);
        }else if(!strcmp(token, "imup.vrrw")){
            reads_ycsv_header_convert(token, &imup->vrrw);
        }else if(!strcmp(token, "imup.Ta")){
            reads_ycsv_header_convert(token, &imup->Ta);
        }else if(!strcmp(token, "imup.Tg")){
            reads_ycsv_header_convert(token, &imup->Tg);
        }else if(!strcmp(token, "imup.kod")){
            if((token = strtok(NULL, ","))) imup->kod = atof(token);
        }else if(!strcmp(token, "imup.kod_std")){
            if((token = strtok(NULL, ","))) imup->kod_std = atof(token);
        }else if(!strcmp(token, "imup.lever_arm_gps")){
            reads_ycsv_header_convert(token, &imup->lever_arm_gps);
        }else if(!strcmp(token, "imup.lever_arm_gps_std")){
            reads_ycsv_header_convert(token, &imup->lever_arm_gps_std);
        }else if(!strcmp(token, "imup.lever_arm_od")){
            reads_ycsv_header_convert(token, &imup->lever_arm_od);
        }else if(!strcmp(token, "imup.lever_arm_od_std")){
            reads_ycsv_header_convert(token, &imup->lever_arm_od_std);
        }else if(!strcmp(token, "imup.lever_arm_car")){
            reads_ycsv_header_convert(token, &imup->lever_arm_car);
        }else if(!strcmp(token, "imup.lever_arm_car_std")){
            reads_ycsv_header_convert(token, &imup->lever_arm_car_std);
        }else if(!strcmp(token, "imup.err_angle_imu")){     /* deg */
            reads_ycsv_header_convert(token, &imup->err_angle_imu);
            imup->err_angle_imu = v3_scalar(DEG2RAD, imup->err_angle_imu);
        }else if(!strcmp(token, "imup.err_angle_imu_std")){   /* deg */
            reads_ycsv_header_convert(token, &imup->err_angle_imu_std);
            imup->err_angle_imu_std =
                    v3_scalar(DEG2RAD, imup->err_angle_imu_std);
        }else if(!strcmp(token, "imup.err_angle_imu_rw")){
            reads_ycsv_header_convert(token, &imup->err_angle_imu_rw);
        }else if(!strcmp(token, "imup.Terr_angle_imu")){
            reads_ycsv_header_convert(token, &imup->Terr_angle_imu);
        }else if(!strcmp(token, "imup.err_angle_gps")){
            reads_ycsv_header_convert(token, &imup->err_angle_gps);
            imup->err_angle_gps = v3_scalar(DEG2RAD, imup->err_angle_gps);
        }else if(!strcmp(token, "imup.err_angle_gps_std")){
            reads_ycsv_header_convert(token, &imup->err_angle_gps_std);
            imup->err_angle_gps_std =
                    v3_scalar(DEG2RAD, imup->err_angle_gps_std);
        }else if(!strcmp(token, "imup.ref_point")){
            reads_ycsv_header_convert(token, &imup->ref_point);
        }
    }
    if(!strcmp(token, "cfg.isx_ba")){
        if((token = strtok(NULL, ","))) cfg.isx_ba = atoi(token);
    }else if(!strcmp(token, "cfg.isx_bg")){
        if((token = strtok(NULL, ","))) cfg.isx_bg = atoi(token);
    }else if(!strcmp(token, "cfg.isx_kax")){
        if((token = strtok(NULL, ","))) cfg.isx_kax = atoi(token);
    }else if(!strcmp(token, "cfg.isx_kay")){
        if((token = strtok(NULL, ","))) cfg.isx_kay = atoi(token);
    }else if(!strcmp(token, "cfg.isx_kaz")){
        if((token = strtok(NULL, ","))) cfg.isx_kaz = atoi(token);
    }else if(!strcmp(token, "cfg.isx_kgx")){
        if((token = strtok(NULL, ","))) cfg.isx_kgx = atoi(token);
    }else if(!strcmp(token, "cfg.isx_kgy")){
        if((token = strtok(NULL, ","))) cfg.isx_kgy = atoi(token);
    }else if(!strcmp(token, "cfg.isx_kgz")){
        if((token = strtok(NULL, ","))) cfg.isx_kgz = atoi(token);
    }else if(!strcmp(token, "cfg.isx_kod")){
        if((token = strtok(NULL, ","))) cfg.isx_kod = atoi(token);
    }else if(!strcmp(token, "cfg.isx_eroll")){
        if((token = strtok(NULL, ","))) cfg.isx_eroll = atoi(token);
    }else if(!strcmp(token, "cfg.isx_epitch")){
        if((token = strtok(NULL, ","))) cfg.isx_epitch = atoi(token);
    }else if(!strcmp(token, "cfg.isx_eyaw")){
        if((token = strtok(NULL, ","))) cfg.isx_eyaw = atoi(token);
    }else if(!strcmp(token, "cfg.isx_armgps")){
        if((token = strtok(NULL, ","))) cfg.isx_armgps = atoi(token);
    }else if(!strcmp(token, "cfg.isx_armcar")){
        if((token = strtok(NULL, ","))) cfg.isx_armcar = atoi(token);
    }else if(!strcmp(token, "cfg.is_odincre")){
        if((token = strtok(NULL, ","))) cfg.is_odincre = atoi(token);
    }else if(!strcmp(token, "cfg.iskf_itgdS_sagehusa")){
        if((token = strtok(NULL, ","))) cfg.iskf_itgdS_sagehusa = atoi(token);
    }else if(!strcmp(token, "cfg.max_ny")){
        if((token = strtok(NULL, ","))) cfg.max_ny = (unsigned char)atoi(token);
    }else if(!strcmp(token, "cfg.nZST")){
        if((token = strtok(NULL, ","))) cfg.nZST = (unsigned short)atoi(token);
    }else if(!strcmp(token, "cfg.issol_header")){
        if((token = strtok(NULL, ","))) cfg.issol_header = atoi(token);
    }else if(!strcmp(token, "cfg.sol_refpos")){
        if((token = strtok(NULL, ",")))
            cfg.sol_refpos = (unsigned int)atoi(token);
    }else if(!strcmp(token, "cfg.feedratio")){
        if((token = strtok(NULL, ","))) cfg.feedratio = (float)atof(token);
    }
    return 0;
}

extern int readf_ycsv_header(FILE *fp, imup_t *imup){
    char buff[MAXLINELEN];
    while (fgets(buff, MAXLINELEN, fp)) {
        if(buff[0] == '>') break;       /* skip comment */
        if(buff[0] == '#') continue;    /* skip comment */
        reads_ycsv_header(buff, imup);  /* */
    }
    return 0;
}

static int readf_ycsv(FILE *fp, imu_t *imu, pva_t *pva, od_t *od)
{
    char buff[MAXLINELEN];
    char seps[] = ",\n";
    char *token;
    unsigned char clm[128] = {0};
    unsigned char count_clm = 0;
    /* week sec */
    while (fgets(buff, MAXLINELEN, fp)) {
        if(buff[0] != '>'){
            if(imu != NULL)
                reads_ycsv_header(buff, imu->property);
            else
                reads_ycsv_header(buff, NULL);
        }else{
            token = strtok((buff+1), seps);
            do{
                while(true){
                    if(token[0] == ' ') token++; else break;
                }
                /* time msg, 4*/
                if(!strcmp(token, "week"))      clm[count_clm] = 1;
                else if(!strcmp(token, "sec"))  clm[count_clm] = 2;
                /* imu data, 8 */
                else if(imu != NULL){
                    if(!strcmp(token, "accel_x"))  clm[count_clm] = 5;
                    else if(!strcmp(token, "accel_y"))  clm[count_clm] = 6;
                    else if(!strcmp(token, "accel_z"))  clm[count_clm] = 7;
                    else if(!strcmp(token, "gyro_x"))   clm[count_clm] = 8;
                    else if(!strcmp(token, "gyro_y"))   clm[count_clm] = 9;
                    else if(!strcmp(token, "gyro_z"))   clm[count_clm] = 10;
                }
                /* od data, 8 */
                else if(od != NULL){
                    if(!strcmp(token, "dS"))       clm[count_clm] = 13;
                }
                /* pva data, 32*/
                else if(pva != NULL){
                    if(!strcmp(token, "lat"))      clm[count_clm] = 19;
                    else if(!strcmp(token, "lon"))      clm[count_clm] = 20;
                    else if(!strcmp(token, "hgt"))      clm[count_clm] = 21;
                    else if(!strcmp(token, "vN"))       clm[count_clm] = 22;
                    else if(!strcmp(token, "vE"))       clm[count_clm] = 23;
                    else if(!strcmp(token, "vD"))       clm[count_clm] = 24;
                    else if(!strcmp(token, "roll"))     clm[count_clm] = 25;
                    else if(!strcmp(token, "pitch"))    clm[count_clm] = 26;
                    else if(!strcmp(token, "yaw"))      clm[count_clm] = 27;
                    else if(!strcmp(token, "status"))   clm[count_clm] = 28;
                    if(pva->is_cov){
                        if(!strcmp(token, "std_lat"))       clm[count_clm] = 33;
                        else if(!strcmp(token, "std_lon"))  clm[count_clm] = 34;
                        else if(!strcmp(token, "std_hgt"))  clm[count_clm] = 35;
                        else if(!strcmp(token, "std_vN"))   clm[count_clm] = 36;
                        else if(!strcmp(token, "std_vE"))   clm[count_clm] = 37;
                        else if(!strcmp(token, "std_vD"))   clm[count_clm] = 38;
                        else if(!strcmp(token, "std_roll")) clm[count_clm] = 39;
                        else if(!strcmp(token, "std_pitch"))clm[count_clm] = 40;
                        else if(!strcmp(token, "std_yaw"))  clm[count_clm] = 41;
                    }
                    if(pva->is_ext){
                        if(!strcmp(token, "yaw2"))          clm[count_clm] = 49;
                        else if(!strcmp(token, "std_yaw2")) clm[count_clm] = 50;
                        else if(!strcmp(token, "pitch2"))   clm[count_clm] = 51;
                        else if(!strcmp(token, "std_pitch2"))clm[count_clm] = 52;
                        else if(!strcmp(token, "ext_status"))clm[count_clm] = 53;
                    }
                }
                count_clm ++;
            }while((token = strtok(NULL, seps)));
            break;
        }
    }
    count_clm = 0;
    imud_t imud; int week = 0; double sec = 0.0;
    double dS;
    bool is_imu = false, is_pva = false, is_od = false;
    while (fgets(buff, MAXLINELEN, fp)) {
        if(pva != NULL && pva->n == pva->nmax){
            pva->nmax *= 2;
            pva_resize(pva, pva->nmax);
        }
        /* init all */
        if(imu != NULL){
            imud.time.time = 0; imud.time.sec = 0.0;
            imud.accel = V0; imud.gyro = V0;
        }
        if(od != NULL){
            week = 0; sec = 0.0;
        }
        if(pva != NULL){
            pva->time[pva->n].time = 0; pva->time[pva->n].sec = 0;
            pva->pos[pva->n] = V0;
            pva->vel[pva->n] = V0;
            pva->att[pva->n] = V0;
            pva->status[pva->n] = SOL_NONE;
            if(pva->is_cov){
                pva->Qpos[pva->n] = O3;
                pva->Qvel[pva->n] = O3;
                pva->Qatt[pva->n] = O3;
            }
            if(pva->is_ext){
                pva->yaw2[pva->n] = 0.0;
                pva->std_yaw2[pva->n] = 0.0;
                pva->pitch2[pva->n] = 0.0;
                pva->std_pitch2[pva->n] = 0.0;
                pva->ext_status[pva->n] = SOL_NONE;
            }
        }
        /* reading */
        token = strtok(buff, seps);
        count_clm = 0;
        do{
            switch (clm[count_clm]) {
            case 1: week = atoi(token); break;
            case 2: sec = atof(token);  break;
            /* imu data */
            case 5: if(imu != NULL) imud.accel.x = atof(token);
                is_imu = true; break;
            case 6: if(imu != NULL) imud.accel.y = atof(token);
                is_imu = true; break;
            case 7: if(imu != NULL) imud.accel.z = atof(token);
                is_imu = true; break;
            case 8: if(imu != NULL) imud.gyro.x = atof(token);
                is_imu = true; break;
            case 9: if(imu != NULL) imud.gyro.y = atof(token);
                is_imu = true; break;
            case 10: if(imu != NULL) imud.gyro.z = atof(token);
                is_imu = true; break;
            /* od data */
            case 13: if(od != NULL) dS = atof(token);
                is_od = true; break;
            /* pva data */
            case 19: if(pva != NULL) pva->pos[pva->n].x = atof(token)*DEG2RAD;
                is_pva = true; break;
            case 20: if(pva != NULL) pva->pos[pva->n].y = atof(token)*DEG2RAD;
                is_pva = true; break;
            case 21: if(pva != NULL) pva->pos[pva->n].z = atof(token);
                is_pva = true; break;
            case 22: if(pva != NULL) pva->vel[pva->n].x = atof(token);
                is_pva = true; break;
            case 23: if(pva != NULL) pva->vel[pva->n].y = atof(token);
                is_pva = true; break;
            case 24: if(pva != NULL) pva->vel[pva->n].z = atof(token);
                is_pva = true; break;
            case 25: if(pva != NULL) pva->att[pva->n].x = atof(token)*DEG2RAD;
                is_pva = true; break;
            case 26: if(pva != NULL) pva->att[pva->n].y = atof(token)*DEG2RAD;
                is_pva = true; break;
            case 27: if(pva != NULL) pva->att[pva->n].z = atof(token)*DEG2RAD;
                is_pva = true; break;
            case 28: if(pva != NULL)
                pva->status[pva->n] = (unsigned int)atoi(token);
                is_pva = true; break;
            case 33:
                if(pva != NULL && pva->is_cov)
                    pva->Qpos[pva->n].m11 = SQR(atof(token));
                is_pva = true; break;
            case 34:
                if(pva != NULL && pva->is_cov)
                    pva->Qpos[pva->n].m22 = SQR(atof(token));
                is_pva = true; break;
            case 35:
                if(pva != NULL && pva->is_cov)
                    pva->Qpos[pva->n].m33 = SQR(atof(token));
                is_pva = true; break;
            case 36:
                if(pva != NULL && pva->is_cov)
                    pva->Qvel[pva->n].m11 = SQR(atof(token));
                is_pva = true; break;
            case 37:
                if(pva != NULL && pva->is_cov)
                    pva->Qvel[pva->n].m22 = SQR(atof(token));
                is_pva = true; break;
            case 38:
                if(pva != NULL && pva->is_cov)
                    pva->Qvel[pva->n].m33 = SQR(atof(token));
                is_pva = true; break;
            case 39:
                if(pva != NULL && pva->is_cov)
                    pva->Qatt[pva->n].m11 = SQR(atof(token)*DEG2RAD);
                is_pva = true; break;
            case 40:
                if(pva != NULL && pva->is_cov)
                    pva->Qatt[pva->n].m22 = SQR(atof(token)*DEG2RAD);
                is_pva = true; break;
            case 41:
                if(pva != NULL && pva->is_cov)
                    pva->Qatt[pva->n].m33 = SQR(atof(token)*DEG2RAD);
                is_pva = true; break;
            case 49:
                if(pva != NULL && pva->is_ext)
                    pva->yaw2[pva->n] = atof(token)*DEG2RAD;
                is_pva = true; break;
            case 50:
                if(pva != NULL && pva->is_ext)
                    pva->std_yaw2[pva->n] = atof(token)*DEG2RAD;
                is_pva = true; break;
            case 51:
                if(pva != NULL && pva->is_ext)
                    pva->pitch2[pva->n] = atof(token)*DEG2RAD;
                is_pva = true; break;
            case 52:
                if(pva != NULL && pva->is_ext)
                    pva->std_pitch2[pva->n] = atof(token)*DEG2RAD;
                is_pva = true; break;
            case 53:
                if(pva != NULL && pva->is_ext)
                    pva->ext_status[pva->n] = (unsigned int)atoi(token);
                is_pva = true; break;
            default: break;
            }
            count_clm ++;
        }while((token = strtok(NULL, seps)));
        gtime_t time = gpst2time(week, sec);
        if(is_imu){ is_imu = false;imud.time = time; imu_add(imu, &imud); }
        if(is_od) { is_od = false; od_add(od, &time, &dS);}
        if(is_pva){
            is_pva = false;
            /* convert NED to ECEF */
            if(pva != NULL && pva->is_cov)
                ned2ecefQ(&pva->pos[pva->n], &pva->Qpos[pva->n],
                        &pva->Qvel[pva->n], &pva->Qatt[pva->n]);
            if(pva != NULL)
                ned2ecef(&pva->pos[pva->n], &pva->vel[pva->n], NULL);
            pva->time[pva->n++] = time;
        }
    }
    return 0;
}

extern int imu_add(imu_t* imu, const imud_t* data)
{
    if(imu->n == 0 && imu->nmax == 0){
        imu->nmax = MAXIMUOBS;
        MALLOC(imu->data, imud_t, imu->nmax);
    }
    else if (imu->nmax <= imu->n) {
        imud_t* imu_data;
        imu->nmax *= 2;
        LOG_INFO("try to reallocate imu->data memory to %i", imu->nmax);
        if (!(imu_data
                = (imud_t*)realloc(imu->data, sizeof(imud_t) * imu->nmax))) {
            FREE(imu->data);
            imu->n = imu->nmax = 0;
            return -1;
        }
        imu->data = imu_data;
    }
    imu->data[imu->n++] = *data;
    return 0;
}

extern int imu_init(imu_t *imu)
{
    LOG_INFO("imu_init: try to allocate memory for imu: %i", MAXIMUOBS);
    MALLOC(imu->property, imup_t, 1);

    imu->n = 0;
    imu->nmax = MAXIMUOBS;
    MALLOC(imu->data, imud_t, imu->nmax);
    return 0;
}

extern void imu_free(imu_t *imu)
{
    if(imu->data != NULL){
        LOG_INFO("imu_free: try to free imu->data ...");
        FREE(imu->data);
    }
    if(imu->property != NULL){
        LOG_INFO("imu_free: try to free imu->property ...");
        FREE(imu->property);
    }
}

extern void od_init(od_t *od)
{
    od->n = 0;
    od->nmax = MAXODOBS;
    LOG_INFO("try to init od, allocate memory: %i", MAXODOBS);

    MALLOC(od->time, gtime_t, od->nmax);
    MALLOC(od->dS, double, od->nmax);
}

extern void od_init1(od_t *od)
{
    od->n = 1; od->nmax = 1;
    MALLOC(od->time, gtime_t, 1);
    MALLOC(od->dS, double, 1);
}

extern void od_add(od_t *od, const gtime_t *time, const double *dS)
{
    if(od->n == 0 && od->nmax == 0)
        od_init(od);
    else if(od->n + 1 > od->nmax){
        od->nmax *= 2;
        LOG_INFO("od_add: try to rellocate od memory to %i", od->nmax);
        REALLOC(od->time, gtime_t, od->nmax);
        REALLOC(od->dS, double, od->nmax);
    }
    od->time[od->n] = *time;
    od->dS[od->n] = *dS;
    od->n++;
}

extern void od_free(od_t *od)
{
    LOG_INFO("od_free: try to free od...");
    FREE(od->time);     FREE(od->dS);
    od->n = 0;          od->nmax = 0;
}

extern void pva_init(pva_t *pva)
{
    pva->n = 0;
    pva->nmax = MAXPVAOBS;
    LOG_INFO("pva_init: try to allocate pva memory: %i", MAXPVAOBS);
    MALLOC(pva->time, gtime_t, pva->nmax);
    MALLOC(pva->pos, v3_t, pva->nmax);
    MALLOC(pva->vel, v3_t, pva->nmax);
    MALLOC(pva->att, v3_t, pva->nmax);
    MALLOC(pva->status,    unsigned int, pva->nmax);

    LOG_INFO("pva_init: pva->is_cov = %s", pva->is_cov ? "true":"false");
    if(pva->is_cov){
        MALLOC(pva->Qpos, m3_t, pva->nmax);
        MALLOC(pva->Qvel, m3_t, pva->nmax);
        MALLOC(pva->Qatt, m3_t, pva->nmax);
    }else{
        LOG_WARN("pva->Qatt, pva->Qvel, pva->Qpos do NOT allocate memory");
    }

    LOG_INFO("pva_init: pva->is_ext = %s", pva->is_ext ? "true":"false");
    if(pva->is_ext){
        MALLOC(pva->yaw2,       double, pva->nmax);
        MALLOC(pva->std_yaw2,   double, pva->nmax);
        MALLOC(pva->pitch2,     double, pva->nmax);
        MALLOC(pva->std_pitch2, double, pva->nmax);
        MALLOC(pva->ext_status,unsigned int, pva->nmax);
    }else{
        LOG_WARN("pva_init: some pva properties do NOT allocate memory");
    }
}

extern void pva_resize(pva_t *pva, unsigned int nmax)
{
    pva->nmax = nmax;
    LOG_INFO("%s: try to reallocate pva memory to %i",__FUNCTION__, nmax);
    REALLOC(pva->time, gtime_t, nmax);
    REALLOC(pva->pos, v3_t, nmax);
    REALLOC(pva->vel, v3_t, nmax);
    REALLOC(pva->att, v3_t, nmax);
    REALLOC(pva->status, unsigned int, nmax);
    if(pva->is_cov){
        REALLOC(pva->Qpos, m3_t, nmax);
        REALLOC(pva->Qvel, m3_t, nmax);
        REALLOC(pva->Qatt, m3_t, nmax);
    }
    if(pva->is_ext){
        REALLOC(pva->yaw2,          double, nmax);
        REALLOC(pva->std_yaw2,      double, nmax);
        REALLOC(pva->pitch2,        double, nmax);
        REALLOC(pva->std_pitch2,    double, nmax);
        REALLOC(pva->ext_status,    unsigned int, nmax);
    }
}

extern void pva_free(pva_t *pva)
{
    LOG_INFO("%s: try to free pva...", __FUNCTION__);
    pva->n = 0;         pva->nmax = 0;
    FREE(pva->att); FREE(pva->att); FREE(pva->vel); FREE(pva->pos);
    FREE(pva->status);
    if(pva->is_cov){
        FREE(pva->Qatt); FREE(pva->Qvel); FREE(pva->Qpos);
    }
    if(pva->is_ext){
        FREE(pva->yaw2);    FREE(pva->std_yaw2);
        FREE(pva->pitch2);  FREE(pva->std_pitch2);
        FREE(pva->ext_status);
    }
}

extern int yins_readf(const char * fname, enum FT ft, imu_t *imu, pva_t *pva,
                      od_t *od)
{
    FILE* fp;
    if (!(fp = fopen(fname, "r"))) {
        LOG_FATAL("%s: Try to file open error: %s",__FUNCTION__, fname);
    }else{
        LOG_INFO("%s: Reading file: %s",__FUNCTION__, fname);
    }
    switch(ft){
    case FT_IMU_CSV:
        readf_imu_csv(fp, imu); break;
    case FT_NVT_ASC:
        readf_imu_nvt(fp, imu); break;
    case FT_YGM_IMU:
        readf_ygm_imu(fp, imu); break;
    case FT_YGM_AVP:
        readf_ygm_avp(fp, pva); break;
    case FT_YGM_OD:
        readf_ygm_od(fp, od); break;
    case FT_IMU_YCSV:
        readf_ycsv(fp, imu, NULL, NULL); break;
    case FT_PVA_YCSV:
        readf_ycsv(fp, imu, pva, NULL); break;
    case FT_OD_YCSV:
        readf_ycsv(fp, imu, NULL, od); break;
    case FT_NMEA:
        readf_nmea(fp,  pva); break;
    case FT_CFG_YCSV:
        if(imu == NULL) readf_ycsv_header(fp, NULL);
        else readf_ycsv_header(fp, imu->property);
        updatecfg();
        break;
    default:
        LOG_ERROR("%s: Can't identify FILETYPE %s",__FUNCTION__, fname);
        fclose(fp);
        return 1;
    }
    fclose(fp);
    return 0;
}

static int writef_imu_csv(FILE *fp, const imu_t *imu)
{
    for (unsigned int i = 0; i < imu->n; ++i) {
        if(imu->data[i].time.time == 0){
            fprintf(fp, "%5i,%10.3f",0, imu->data[i].time.sec);
        }else{
            int week;
            double sec = time2gpst(imu->data[i].time, &week);
            fprintf(fp, "%5i,%10.3f",week, sec);
        }
        fprintf(fp, ",%16.10f,%16.10f,%16.10f", imu->data[i].accel.x,
                imu->data[i].accel.y, imu->data[i].accel.z);
        fprintf(fp, ",%16.10f,%16.10f,%16.10f\n", imu->data[i].gyro.x,
                imu->data[i].gyro.y, imu->data[i].gyro.z);
    }
    fflush(fp);
    return 0;
}

static int writef_imu_ycsv(FILE *fp, const imu_t *imu)
{
    writef_imup(fp, imu->property);
    if(imu->n <= 0) return 0;
    fprintf(fp, ">");
    fprintf(fp, "week,%10s", "sec");
    fprintf(fp, ",%16s,%16s,%16s", "accel_x", "accel_y", "accel_z");
    fprintf(fp, ",%16s,%16s,%16s", "gyro_x", "gyro_y", "gyro_z");
    fprintf(fp, "\n");
    return writef_imu_csv(fp, imu);
}

static int writef_od_csv(FILE *fp, const od_t *od)
{
    for(unsigned int i = 0; i < od->n; ++i){
        if(od->time[i].time == 0){
            fprintf(fp, "%5i,%10.3f",0, od->time[i].sec);
        }else{
            int week;
            double sec = time2gpst(od->time[i], &week);
            fprintf(fp, "%5i,%10.3f",week, sec);
        }
        fprintf(fp, ",%10.4f", od->dS[i]);
        fprintf(fp, "\n");
    }
    fflush(fp);
    return  0;
}

static int writef_od_ycsv(FILE *fp, const od_t *od)
{
    fprintf(fp, ">");
    fprintf(fp, "week,%10s", "sec");
    fprintf(fp, ",%10s", "dS");
    fprintf(fp, "\n");
    return writef_od_csv(fp,od);
}

static int writef_pva_csv(FILE *fp, const pva_t *pva)
{
    v3_t pos, vel, att; m3_t Qpos, Qvel, Qatt;
    for(unsigned int i = 0; i < pva->n; ++i){
        if(pva->time[i].time == 0){
            fprintf(fp, "%5i,%10.3f",0, pva->time[i].sec);
        }else{
            int week;
            double sec = time2gpst(pva->time[i], &week);
            fprintf(fp, "%5i,%10.3f",week, sec);
        }
        pos = pva->pos[i]; vel = pva->vel[i];
        ecef2ned(&pos, &vel, NULL);
        fprintf(fp, ",%12.7f,%12.7f,%12.4f",
                pos.x*RAD2DEG, pos.y*RAD2DEG, pos.z);
        fprintf(fp, ",%12.4f,%12.4f,%12.4f",
                vel.x, vel.y, vel.z);
        fprintf(fp, ",%8.4f,%8.4f,%8.4f", pva->att[i].x*RAD2DEG,
                pva->att[i].y*RAD2DEG, pva->att[i].z*RAD2DEG);
        fprintf(fp, ",%6d", pva->status[i]);
        if(pva->is_cov){
            Qpos = pva->Qpos[i]; Qvel = pva->Qvel[i];
            ecef2nedQ(&pva->pos[i], &Qpos, &Qvel, NULL);
            fprintf(fp, ",%12.5G,%12.5G,%12.5G", sqrt(Qpos.m11), sqrt(Qpos.m22),
                    sqrt(Qpos.m33));
            fprintf(fp, ",%12.5G,%12.5G,%12.5G", sqrt(Qvel.m11), sqrt(Qvel.m22),
                    sqrt(Qvel.m33));
            fprintf(fp, ",%12.5G,%12.5G,%12.5G", sqrt(pva->Qatt[i].m11)*RAD2DEG,
                    sqrt(pva->Qatt[i].m22)*RAD2DEG,
                    sqrt(pva->Qatt[i].m33)*RAD2DEG);
        }
        if(pva->is_ext){
            fprintf(fp, ",%10.4f,%10.4f",
                    pva->yaw2[i]*RAD2DEG, pva->std_yaw2[i]*RAD2DEG);
            fprintf(fp, ",%10.4f,%10.4f",
                    pva->pitch2[i]*RAD2DEG, pva->std_pitch2[i]*RAD2DEG);
            fprintf(fp, ",%10d", pva->ext_status[i]);
        }
        fprintf(fp, "\n");
    }
    fflush(fp);
    return  0;
}

static int writef_pva_ycsv(FILE *fp, const pva_t *pva)
{
    writef_cfg(fp);
    fprintf(fp, ">");
    fprintf(fp, "week,%10s", "sec");
    fprintf(fp, ",%12s,%12s,%12s", "lat","lon","hgt");
    fprintf(fp, ",%12s,%12s,%12s", "vN","vE","vD");
    fprintf(fp, ",%8s,%8s,%8s", "roll","pitch","yaw");
    fprintf(fp,",%6s", "status");
    if(pva->is_cov){
        fprintf(fp, ",%12s,%12s,%12s", "std_lat","std_lon","std_hgt");
        fprintf(fp, ",%12s,%12s,%12s", "std_vN","std_vE","std_vD");
        fprintf(fp, ",%12s,%12s,%12s", "std_roll","std_pitch","std_yaw");
    }
    if(pva->is_ext){
        fprintf(fp, ",%10s,%10s", "yaw2", "std_yaw2");
        fprintf(fp, ",%10s,%10s", "pitch2", "std_pitch2");
        fprintf(fp,",%10s", "ext_status");
    }
    fprintf(fp, "\n");
    return writef_pva_csv(fp, pva);
}

extern int yins_writef(const char *fname, enum FT FILETYPE, const imu_t *imu,
                       const pva_t *pva, const od_t *od)
{
    FILE *fp;
    if (!(fp = fopen(fname, "w"))){
        LOG_FATAL("Try to open file error: %s", fname);
    }else{
        LOG_INFO("%s: Write to file: %s", __FUNCTION__ ,fname);
    }
    switch (FILETYPE) {
    case FT_IMU_CSV:
        writef_imu_csv(fp, imu); break;
    case FT_OD_CSV:
        writef_od_csv(fp, od); break;
    case FT_IMU_YCSV:
        writef_imu_ycsv(fp, imu); break;
    case FT_OD_YCSV:
        writef_od_ycsv(fp, od); break;
    case FT_PVA_CSV:
        writef_pva_csv(fp, pva); break;
    case FT_PVA_YCSV:
        writef_pva_ycsv(fp, pva); break;
    case FT_CFG_YCSV:
        if(imu != NULL && imu->property != NULL)
            writef_imup(fp, imu->property);
        writef_cfg(fp);
        break;
    default:
        LOG_ERROR("Can NOT support current filetype: %i", FILETYPE);
    }
    fclose(fp);
    return 0;
}

/**
 * @brief output ins solution struct to file
 * @param[in] 	fp	 	file pointer
 * @param[in] 	sol		solution struct
 * @param[in] 	imup 	imu propertty struct
 * @return status(0: OK)
 * @note this function control by global cfg varibale, related properties: \n
 * 		cfg.sol_refpos		solution ouput reference point \n
 * 		cfg.issol_header	output solution header or not(stanadard csv file) \n
 * 		cfg.isx_*			output correspoding property or not
 */
extern int outsolins(FILE *fp, const solins_t *sol, const imup_t *imup)
{
    solins_t outsol = *sol;
    switch (cfg.sol_refpos) {
    case REFPOS_IMU:    break;
    case REFPOS_MANUAL: outsol.pos =
            v3_add(outsol.pos, quat_mul_v3(outsol.quat, imup->ref_point));
        break;
    case REFPOS_GPS: outsol.pos =
            v3_add(outsol.pos, quat_mul_v3(outsol.quat, imup->lever_arm_gps));
        break;
    case REFPOS_CAR: outsol.pos =
            v3_add(outsol.pos, quat_mul_v3(outsol.quat, imup->lever_arm_car));
        break;
    }

    ecef2nedQ(&outsol.pos, &outsol.Qpos, &outsol.Qvel, &outsol.Qatt);
    ecef2ned(&outsol.pos, &outsol.vel, &outsol.dcm);
    v3_t Enb;  dcm2att(&outsol.dcm, &Enb);

    /* output header (ycsv format) */
    if(cfg.issol_header){
        writef_imup(fp, imup);
        writef_cfg(fp);
        fprintf(fp, ">");
        fprintf(fp, "week,%10s", "sec");
        fprintf(fp, ",%12s,%12s,%12s", "lat", "lon", "hgt");
        fprintf(fp, ",%10s,%10s,%10s", "vN", "vE", "vD");
        fprintf(fp, ",%8s,%8s,%8s", "roll", "pitch", "yaw");
        fprintf(fp, ",%10s,%10s,%10s", "std_lat", "std_lon", "std_hgt");
        fprintf(fp, ",%10s,%10s,%10s", "std_vN", "std_vE", "std_vU");
        fprintf(fp, ",%10s,%10s,%10s", "std_roll", "std_pitch", "std_yaw");
        if(cfg.isx_ba){
            fprintf(fp, ",%12s,%12s,%12s", "ba_x", "ba_y", "ba_z");
            fprintf(fp, ",%10s,%10s,%10s", "std_bax", "std_bay", "std_baz");
        }
        if(cfg.isx_bg){
            fprintf(fp, ",%12s,%12s,%12s", "bg_x", "bg_y", "bg_z");
            fprintf(fp, ",%10s,%10s,%10s", "std_bgx", "std_bgy", "std_bgz");
        }
        if(cfg.isx_kod){
            fprintf(fp, ",%8s", "kod");
            fprintf(fp, ",%10s", "std_kod");
        }
        if(cfg.isx_eroll){
            fprintf(fp, ",%8s", "eroll");
            fprintf(fp, ",%10s", "std_eroll");
        }
        if(cfg.isx_epitch){
            fprintf(fp, ",%8s", "epitch");
            fprintf(fp, ",%10s", "std_epitch");
        }
        if(cfg.isx_eyaw){
            fprintf(fp, ",%8s", "eyaw");
            fprintf(fp, ",%10s", "std_eyaw");
        }
        if(cfg.isx_kax){
            fprintf(fp, ",%8s", "kax");
            fprintf(fp, ",%10s", "std_kax");
        }
        if(cfg.isx_kay){
            fprintf(fp, ",%8s", "kay");
            fprintf(fp, ",%10s", "std_kay");
        }
        if(cfg.isx_kaz){
            fprintf(fp, ",%8s", "kaz");
            fprintf(fp,",%10s", "std_kaz");
        }
        if(cfg.isx_kgx){
            fprintf(fp, ",%8s", "kgx");
            fprintf(fp, ",%10s", "std_kgx");
        }
        if(cfg.isx_kgy){
            fprintf(fp, ",%8s", "kgy");
            fprintf(fp, ",%10s", "std_kgy");
        }
        if(cfg.isx_kgz){
            fprintf(fp, ",%8s", "kgz");
            fprintf(fp, ",%10s", "std_kgz");
        }
        fprintf(fp, ",%6s","status");
        fprintf(fp, "\n");
        cfg.issol_header = false;
    }

    if(outsol.time.time == 0.0){
        fprintf(fp, "%5i,%10.3f",0, outsol.time.sec);
    }else{
        int week;
        double sec = time2gpst(outsol.time, &week);
        fprintf(fp, "%5i,%10.3f",week, sec);
    }
    fprintf(fp, ",%12.8f,%12.8f,%12.8f",
            outsol.pos.x*RAD2DEG, outsol.pos.y*RAD2DEG, outsol.pos.z);
    fprintf(fp, ",%10.6f,%10.6f,%10.6f",
            outsol.vel.x, outsol.vel.y, outsol.vel.z);
    fprintf(fp, ",%8.4f,%8.4f,%8.4f",
            Enb.x*RAD2DEG, Enb.y*RAD2DEG, Enb.z*RAD2DEG);
    fprintf(fp, ",%10.5G,%10.5G,%10.5G", sqrt(outsol.Qpos.m11),
            sqrt(outsol.Qpos.m22), sqrt(outsol.Qpos.m33));
    fprintf(fp, ",%10.5G,%10.5G,%10.5G", sqrt(outsol.Qvel.m11),
            sqrt(outsol.Qvel.m22), sqrt(outsol.Qvel.m33));
    fprintf(fp, ",%10.5G,%10.5G,%10.5G", sqrt(outsol.Qatt.m11)*RAD2DEG,
            sqrt(outsol.Qatt.m22)*RAD2DEG, sqrt(outsol.Qatt.m33)*RAD2DEG);
    if(cfg.isx_ba){
        fprintf(fp, ",%12.6G,%12.6G,%12.6G",
            outsol.ba.x*MPS22MG, outsol.ba.y*MPS22MG, outsol.ba.z*MPS22MG);
        fprintf(fp, ",%10.5G,%10.5G,%10.5G", outsol.ba_std.x*MPS22MG,
            outsol.ba_std.y*MPS22MG, outsol.ba_std.z*MPS22MG);
    }
    if(cfg.isx_bg){
        fprintf(fp, ",%12.6G,%12.6G,%12.6G",
            outsol.bg.x*RPS2DPH, outsol.bg.y*RPS2DPH, outsol.bg.z*RPS2DPH);
        fprintf(fp,",%10.5G,%10.5G,%10.5G", outsol.bg_std.x*RPS2DPH,
            outsol.bg_std.y*RPS2DPH, outsol.bg_std.z*RPS2DPH);
    }
    if(cfg.isx_kod){
        fprintf(fp, ",%8.4f", 1.0/outsol.kod);
        fprintf(fp,",%10.5G", outsol.std_kod);
    }
    if(cfg.isx_eroll || cfg.isx_epitch  || cfg.isx_eyaw ){
        v3_t Ebc;  dcm2euler(&outsol.Cbc, &Ebc);
        if(cfg.isx_eroll){
            fprintf(fp, ",%8.4f", angle_to180(Ebc.x)*RAD2DEG);
            fprintf(fp, ",%10.5G", outsol.std_Cbc.x*RAD2DEG);
        }
        if(cfg.isx_epitch){
            fprintf(fp, ",%8.4f", angle_to180(Ebc.y)*RAD2DEG);
            fprintf(fp, ",%10.5G", outsol.std_Cbc.y*RAD2DEG);
        }
        if(cfg.isx_eyaw){
            fprintf(fp, ",%8.4f", angle_to180(Ebc.z)*RAD2DEG);
            fprintf(fp, ",%10.5G", outsol.std_Cbc.z*RAD2DEG);
        }
    }
    if(cfg.isx_kax){
        fprintf(fp, ",%8.6f", outsol.ka.x);
        fprintf(fp, ",%10.5G", outsol.ka_std.x);
    }
    if(cfg.isx_kay){
        fprintf(fp, ",%8.6f", outsol.ka.y);
        fprintf(fp, ",%10.5G", outsol.ka_std.y);
    }
    if(cfg.isx_kaz){
        fprintf(fp, ",%8.6f", outsol.ka.z);
        fprintf(fp, ",%10.5G", outsol.ka_std.z);
    }
    if(cfg.isx_kgx){
        fprintf(fp, ",%8.6f", outsol.kg.x);
        fprintf(fp, ",%10.5G", outsol.kg_std.x);
    }
    if(cfg.isx_kgy){
        fprintf(fp, ",%8.6f", outsol.kg.y);
        fprintf(fp, ",%10.5G", outsol.kg_std.y);
    }
    if(cfg.isx_kgz){
        fprintf(fp, ",%8.6f", outsol.kg.z);
        fprintf(fp, ",%10.5G", outsol.kg_std.z);
    }
    fprintf(fp, ",%6i", outsol.status);
    fprintf(fp, "\n");
    fflush(fp);

    return 0;
}
