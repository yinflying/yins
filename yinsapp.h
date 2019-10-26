/**
 * @file yinsapp.h
 * @brief ins application functions
 * @details
 * 	In order to hide nearly all details about yins core function,
 * 	yinsapp.h/yinapp.c provide a new protocol struct and and function
 * 	interface for users(It is easy to call under microcomputer,
 *  e.g. Raspberry Pi)
 *
 *  This file is an example how to write an example to develop by using
 * 	ins funtion.
 *
 * @author yinflying(yinflying@foxmail.com)
 * @version 0.0.1
 * @note
 *  2019-10-06  Created, Add yinsapp_pureins()
 * 	2019-10-16 	Add yinsapp_process() and yinsapp_data_reset()
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

#ifndef YINSAPP_H
#define YINSAPP_H
#include "yins_core/ins.h"

typedef enum{
    YINS_GNSS_STAT_NONE = 0,
    YINS_GNSS_STAT_SINGLE = 5,
    YINS_GNSS_STAT_DGNSS = 4,
    YINS_GNSS_STAT_PPP = 3,
    YINS_GNSS_STAT_FLOAT = 2,
    YINS_GNSS_STAT_FIX = 1
}YINS_GNSS_STAT;

typedef enum{
    YINS_YAW2ANT_STAT_NONE = 0,
    YINS_YAW2ANT_STAT_FIX = 1
}YINS_YAW2ANT_STAT;

typedef struct{
    int week;
    double sow;

    bool isa_imu;
    bool isa_pps;
    double gyro[3];
    double accel[3];

    bool isa_gnss;
    double pos[3];
    double pos_std[3];
    double veg_n[3];
    double veg_n_std[3];
    YINS_GNSS_STAT gnss_stat;
    double yaw2ant;
    double yaw2ant_std;
    YINS_YAW2ANT_STAT yaw2ant_stat;

    bool isa_od;
    double dS;
}yinsapp_data_t;

typedef struct{
    int week;
    double sow;
    double pos[3];
    double pos_std[3];
    double veb_n[3];
    double veb_n_std[3];
    double att[3];
    double att_std[3];
    unsigned int stat;
}yinsapp_result_t;

int yinsapp_process(const yinsapp_data_t *appdata,
                    yinsapp_result_t *result, const char *cfg_file);

/**
 * @brief reset yinsapp_data_t sturct to zero status
 * @param[in,out] appdata  reseted yinsapp_data_t struct
 */
void yinsapp_data_reset(yinsapp_data_t *appdata);

/**
 * @brief pure ins navigation
 * @param[in] 	fin	 	input imu data file
 * @param[in] 	ft		input imu data file type
 * @param[in] 	imup	imup property struct
 * @param[in] 	fsol 	output solution file
 * @return status(0: OK)
 * @note related imu property:
 * 		imup.freq_imu 	imu data sample frequanency[Hz]
 * 		imup.tstart		pure ins start time
 * 		imup.initr		pure ins start position, BLH[rad, m]
 * 		imup.initv		pure ins start velocity, vN, vE, vD[m/s]
 * 		imup.inita		pure ins start attitude, Ebn[rad]
 */
int yinsapp_pureins(const char *fin, enum FT ft, const imup_t *imup,
                    const char *fsol);

#endif // YINSAPP_H
