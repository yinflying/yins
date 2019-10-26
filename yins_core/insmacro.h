/**
 * @file insmacro.h
 * @brief ins common include header and macro definition
 * @author yinflying(yinflying@foxmail.com)
 * @version 0.0.1
 * @date 2019-10-16
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
#ifndef INSMACRO_H
#define INSMACRO_H

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <ctype.h>

#define DEG2RAD     0.0174532925199433      /**< deg to rad */
#define RAD2DEG     57.2957795130823        /**< rad to deg */
#define DPH2RPS     4.84813681109536e-06    /**< deg/h to rad/s */
#define RPS2DPH     206264.806247096        /**< rad/s to deg/h */
#define DPSH2RPSS   2.90888208665722e-4     /**< deg/sqrt(h) to rad/sqrt(s) */
#define RPSS2DPSH   3437.74677078493        /**< rad/sqrt(s) to deg/sqrt(h) */
#define DPHPSHZ2RPSS 4.84813681109536e-06   /**< deg/h/sqrt(Hz) to rad/sqrt(s) */
#define RPSS2DPHPSHZ 206264.806247096       /**< rad/sqrt(s) to deg/h/sqrt(Hz) */

#define G2MPS2      9.80665                 /**< g0 to m/s^2 */
#define MPS22G      0.101971621297793       /**< m/s^2 to g0 */
#define MG2MPS2     9.80665e-3              /**< millo-g0(mg) to m/s^2 */
#define MPS22MG     101.971621297793        /**< m/s^2 to milli-g0(mg) */
#define UG2MPS2     9.80665e-6              /**< micro-g0(ug) to m/s^2 */
#define MPS22UG     101971.621297793        /**< m/s^2 to micro-g0(ug) */
#define GAL2MPS2    0.01                    /**< gal to m/s^2 */
#define MPS22GAL    100.0                   /**< m/s^2 to gal */

#define PI 3.14159265358979
#define SQR(x) ((x) * (x))
#define EPS 1E-50							/**< Zero determination threshold */

#endif // INSMACRO_H
