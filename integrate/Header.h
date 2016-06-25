//
//  Header.h
//  integrate
//
//  Created by Дмитрий on 15.06.16.
//  Copyright © 2016 Дмитрий. All rights reserved.
//

#ifndef Header_h
#define Header_h
/*!
 \file
 \brief Заголовочный файл с описанием функций
 
 Данный файл содержит в себе определения основных
 функций, используемых в демонстрационной программе
 */
typedef double(*function)(double);
typedef double(*dfunction)(double, double);
double integral_pram(function f, double a, double b, unsigned step_count);
double integral_pram_inf(function f, double a, double b, double h, double eps);
double integral_trap(function f, double a, double b, unsigned step_count);
double integral_trap_inf(function f, double a, double b, double h, double eps);
double integral_simp(function f, double a, double b, unsigned step_count);
double integral_simp_inf(function f, double a, double b, double h, double eps);
double integral_monte(function f, double a, double b, unsigned step_count);
void integral_runge4(dfunction f, double x0, double x1, double y0, double* x, double* y, double h);
void integral_runge5(dfunction f, double x0, double x1, double y0, double* x, double* y, double h);
void integral_eiler(dfunction f, double x0, double x1, double y0, double* x, double* y, double h);
void mnk(int n, double *x, double *y, double *a_res, double *b_res);
#endif /* Header_h */