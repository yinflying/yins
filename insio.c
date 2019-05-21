#include "ins.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXIMULEN 256
#define NINCOBS 10000
#define MAXIMU 10000
#define DEG2RAD 0.0174532925199433
#define RAD2DEG 57.2957795130823

static const double gpst0[] = { 1980, 1, 6, 0, 0, 0 }; /* gps time reference */
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
/* convert calendar day/time to time -------------------------------------------
 * convert calendar day/time to gtime_t struct
 * args   : double *ep       I   day/time {year,month,day,hour,min,sec}
 * return : gtime_t struct
 * notes  : proper in 1970-2037 or 1970-2099 (64bit time_t)
 *-----------------------------------------------------------------------------*/
static gtime_t epoch2time(const double* ep)
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
static gtime_t gpst2time(int week, double sec)
{
    gtime_t t = epoch2time(gpst0);

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
static void time2epoch(gtime_t t, double* ep)
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

static int readimu_csv(FILE* fp, imu_t* imu)
{
    if (!(imu->data = (imud_t*)malloc(sizeof(imud_t) * MAXIMU)))
        return 0;
    imu->n = 0;
    imu->nmax = MAXIMU;

    char buff[MAXIMULEN];
    imud_t data;
    while (fgets(buff, MAXIMULEN, fp)) {
        data.time.sec = str2num(buff, 0, 16);
        data.accel.i = str2num(buff, 17, 16);
        data.accel.j = str2num(buff, 34, 16);
        data.accel.k = str2num(buff, 51, 16);
        data.gryo.i = str2num(buff, 68, 16);
        data.gryo.j = str2num(buff, 85, 16);
        data.gryo.k = str2num(buff, 102, 16);
        addimudata(imu, &data);
    }
    return 0;
}

/* read Novatel ascii file*/
static int readimu_nvt(FILE* fp, imu_t* imu)
{
    if (!(imu->data = (imud_t*)malloc(sizeof(imud_t) * MAXIMU)))
        return 0;
    imu->n = 0;
    imu->nmax = MAXIMU;

    char buff[MAXIMULEN];
    imud_t data;
    char seps[] = ",;*";
    char* token;
    int gpsw;

    double a_scale = 1.52587890625E-06;     /* 0.05*2^-15 */
    double g_scale = 1.085069444444444E-07; /* 0.1/(3600*256) */

    /* Search the INSPVASA with INS_SOLUTION_GOOD status */
    vec3_t r, v, a;
    while (fgets(buff, MAXIMULEN, fp)) {
        token = strtok(buff, seps);
        if (!strcmp(token, "%INSPVASA")) {
            for (int i = 0; i < 5; ++i)
                token = strtok(NULL, seps); /* skip header */
            r.i = atof(token) * DEG2RAD;
            token = strtok(NULL, seps);
            r.j = atof(token) * DEG2RAD;
            token = strtok(NULL, seps);
            r.k = atof(token);
            token = strtok(NULL, seps);
            v.i = atof(token);
            token = strtok(NULL, seps);
            v.j = atof(token);
            token = strtok(NULL, seps);
            v.k = atof(token);
            token = strtok(NULL, seps);
            a.i = atof(token) * DEG2RAD;
            token = strtok(NULL, seps);
            a.j = atof(token) * DEG2RAD;
            token = strtok(NULL, seps);
            a.k = atof(token) * DEG2RAD;
            token = strtok(NULL, seps);
            if (!strcmp(token, "INS_SOLUTION_GOOD"))
                break;
        }
    }

    dcm_t dcm;
    euler2dcm(&a, &dcm);
    ned2ecef(&r, &v, &dcm);
    dcm2euler(&dcm, &a);
    imu->initr = r;
    imu->initv = v;
    imu->inita = a;

    while (fgets(buff, MAXIMULEN, fp)) {
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
            data.accel.k = -atof(token) * a_scale;
            token = strtok(NULL, seps);
            data.accel.i = -atof(token) * a_scale;
            token = strtok(NULL, seps);
            data.accel.j = atof(token) * a_scale;
            token = strtok(NULL, seps);
            /* angle increment, rad */
            data.gryo.k = -atof(token) * g_scale;
            token = strtok(NULL, seps);
            data.gryo.i = -atof(token) * g_scale;
            token = strtok(NULL, seps);
            data.gryo.j = atof(token) * g_scale;
            addimudata(imu, &data);
        }
    }
    return 0;
}

extern int addimudata(imu_t* imu, const imud_t* data)
{
    imud_t* imu_data;
    if (imu->nmax <= imu->n) {
        if (imu->nmax <= 0)
            imu->nmax = NINCOBS;
        else
            imu->nmax *= 2;
        if (!(imu_data
                = (imud_t*)realloc(imu->data, sizeof(imud_t) * imu->nmax))) {
            free(imu->data);
            imu->data = NULL;
            imu->n = imu->nmax = 0;
            return -1;
        }
        imu->data = imu_data;
    }
    imu->data[imu->n++] = *data;
    return 1;
}

extern void freeimu(imu_t* imu) { free(imu->data); }

extern int readimu_file(const char* infile, imu_t* imu, int FILETYPE)
{
    FILE* fp;
    if (!(fp = fopen(infile, "r"))) {
        fprintf(stderr, "imu file open error\n");
        return 1;
    }
    if (FILETYPE == FT_CSV)
        readimu_csv(fp, imu);
    else if (FILETYPE == FT_NVT)
        readimu_nvt(fp, imu);
    else
        fprintf(stderr, "Can't identify FILETYPE\n");
    fclose(fp);
    return 0;
}

static void fprintf_fixwidth(FILE* fp, double num, int fixwidth)
{
    char buf[30];
    sprintf(buf, "%.*f", fixwidth, num);
    buf[fixwidth] = '\0';
    fprintf(fp, "%s", buf);
}

extern int imu_trans_rnx(const imu_t* imu, const char* outfile)
{
    FILE* fp;
    double ep[6];
    if (!(fp = fopen(outfile, "w"))) {
        fprintf(stderr, "OUTPUT IMU RINEX open error\n");
        return 1;
    }
    fprintf(fp, "%9.2f%-11s%-20s%-20s%-20s\n", 3.02, "", "OBSERVATION DATA",
        "U", "RINEX VERSION / TYPE");
    fprintf(
        fp, "%-20.20s%-40.40s%-20s\n", "yinflying", "CUG", "OBSERVER / AGENCY");

    /* Position */
    fprintf(fp, "%14.4f%14.4f%14.4f%-18s%-20s\n", imu->initr.i, imu->initr.j,
        imu->initr.k, "", "INITIAL POSITION XYZ");
    fprintf(fp, "%14.6f%14.4f%14.6f%-18s%-20s\n", imu->initQr.i, imu->initQr.j,
        imu->initQr.k, "", "INITIAL POSITION STD");
    /*Velocity */
    fprintf(fp, "%14.6f%14.6f%14.6f%-18s%-20s\n", imu->initv.i, imu->initv.j,
        imu->initv.k, "", "INITIAL VELOCITY XYZ");
    fprintf(fp, "%14.6f%14.6f%14.6f%-18s%-20s\n", imu->initQv.i, imu->initQv.j,
        imu->initQv.k, "", "INITIAL VELOCITY STD");
    /* Attitude */
    fprintf(fp, "%14.6f%14.8f%14.6f%-18s%-20s\n", imu->inita.i, imu->inita.j,
        imu->inita.k, "", "INITIAL ATTITUDE CBE");
    fprintf(fp, "%14.6f%14.8f%14.6f%-18s%-20s\n", imu->initQa.i, imu->initQa.j,
        imu->initQa.k, "", "INITIAL ATTITUDE STD");

    /* Bias */
    fprintf(fp, "%14.7G%14.7G%14.7G%-18s%-20s\n", imu->initQab.i, imu->initQab.j,
        imu->initQab.k, "", "INITIAL ACCEBIAS STD");
    fprintf(fp, "%14.7G%14.7G%14.7G%-18s%-20s\n", imu->initQgb.i, imu->initQgb.j,
        imu->initQgb.k, "", "INITIAL GRROBIAS STD");

    /* In run uncertainty */
    fprintf(fp, "%14.7G%14.7G%14.7G%-18s%-20s\n", imu->arw.i, imu->arw.j,
        imu->arw.k, "", "GYRO ANGULAR RW");
    fprintf(fp, "%14.7G%14.7G%14.7G%-18s%-20s\n", imu->arrw.i, imu->arrw.j,
        imu->arrw.k, "", "GYRO ANGULAR RATE RW");
    fprintf(fp, "%14.7G%14.7G%14.7G%-18s%-20s\n", imu->vrw.i, imu->vrw.j,
        imu->vrw.k, "", "ACCE VELOCITY RW");
    fprintf(fp, "%14.7G%14.7G%14.7G%-18s%-20s\n", imu->vrrw.i, imu->vrrw.j,
        imu->vrrw.k, "", "ACCE ACCELERATE RW");

    /* Correlation Time */
    fprintf(fp, "%14.7G%14.7G%14.7G%-18s%-20s\n", imu->Ta.i, imu->Ta.j,
        imu->Ta.k, "", "ACCE CORR TIME");
    fprintf(fp, "%14.7G%14.7G%14.7G%-18s%-20s\n", imu->Tg.i, imu->Tg.j,
        imu->Tg.k, "", "GRYO CORR TIME");

    fprintf(fp, "%1s  %3i%4s%4s%4s%4s%4s%4s%30s%-20s\n", "U", 6, "A1X", "A1Y",
        "A1Z", "G1X", "G1Y", "G1Z", "", "SYS / # / OBS TYPES");
    fprintf(fp, "%-60s%-20s\n", "UNIT: GRYO(RAD) ACCEL(M/S)", "COMMENT");

    /* First&End epoch */
    time2epoch(imu->data[0].time, ep);
    fprintf(fp, "  %04.0f%6.0f%6.0f%6.0f%6.0f%13.7f     %-12s%-20s\n", ep[0],
        ep[1], ep[2], ep[3], ep[4], ep[5], "GPST", "TIME OF FIRST OBS");
    time2epoch(imu->data[imu->n - 1].time, ep);
    fprintf(fp, "  %04.0f%6.0f%6.0f%6.0f%6.0f%13.7f     %-12s%-20s\n", ep[0],
        ep[1], ep[2], ep[3], ep[4], ep[5], "GPST", "TIME OF LAST OBS");
    fprintf(fp, "%-60.60s%-20s\n", "", "END OF HEADER");

    for (int i = 0; i < imu->n; ++i) {
        time2epoch(imu->data[i].time, ep);
        fprintf(fp, "> %04.0f %2.0f %2.0f %2.0f %2.0f%11.7f  %d%3d%21s\n",
            ep[0], ep[1], ep[2], ep[3], ep[4], ep[5], 0, 1, "");
        fprintf(fp, "%1s%02i", "U", 1);
        fprintf_fixwidth(fp, imu->data[i].accel.i, 14);
        fprintf(fp, "  ");
        fprintf_fixwidth(fp, imu->data[i].accel.j, 14);
        fprintf(fp, "  ");
        fprintf_fixwidth(fp, imu->data[i].accel.k, 14);
        fprintf(fp, "  ");
        fprintf_fixwidth(fp, imu->data[i].gryo.i, 14);
        fprintf(fp, "  ");
        fprintf_fixwidth(fp, imu->data[i].gryo.j, 14);
        fprintf(fp, "  ");
        fprintf_fixwidth(fp, imu->data[i].gryo.k, 14);
        fprintf(fp, "  \n");
    }
    fclose(fp);
    return 0;
}
