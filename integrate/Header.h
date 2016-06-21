//
//  Header.h
//  integrate
//
//  Created by Дмитрий on 15.06.16.
//  Copyright © 2016 Дмитрий. All rights reserved.
//

#ifndef Header_h
#define Header_h
typedef double(*function)(double);
typedef double(*dfunction)(double, double);
double integral_pram(function f, double a, double b, unsigned step_count);
double integral_trap(function f, double a, double b, unsigned step_count);
double integral_simp(function f, double a, double b, unsigned step_count);
double integral_monte(function f, double a, double b, unsigned step_count);
void integral_runge4(dfunction f, double* x0, double* y0, double h);
void integral_runge5(dfunction f, double* x0, double* y0, double h);
#endif /* Header_h */