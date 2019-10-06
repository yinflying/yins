/**
 * @file yinsapp.h
 * @brief ins application functions
 * @author yinflying(yinflying@foxmail.com)
 * @version 0.0.1
 * @note
 *  2019-10-06  Created, Add yinsapp_pureins()
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
#include <yins_core/ins.h>

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
int yinsapp_pureins(const char *fin, enum FT ft, const imup_t *imup, const char *fsol);

#endif // YINSAPP_H
