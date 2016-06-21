//
//  main.c
//  integrate
//
//  Created by Дмитрий on 15.06.16.
//  Copyright © 2016 Дмитрий. All rights reserved.
//

#include "Header.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
double integral_pram(function f, double a, double b, unsigned count_step)
{
    double h,S,x;
    int i;
    h=(b-a)/(1.0*count_step);
    S=0;
    for(i=0;i<count_step-1;i++)
    {
        x=a+i*h;
        S=S+f(x);
    }
    S=h*S;
    return S;
}

double integral_trap(function f, double a, double b, unsigned step_count)
{
    double sum = .0, step;
    int i;
    if (0 == step_count) return sum;
    
    step = (b - a)/(1.0*step_count);
    for (i = 1;i<step_count;++i) {
        sum += f (a + i * step);
    }
    sum += (f(a) + f(b)) / 2;
    sum *= step;
    return sum;
}

double integral_simp(function f, double a, double b, unsigned step_count)
{
    double h,sum,x0,x1;
    int i;
    h = (b - a)/(1.0*step_count);
    sum = 0;
    x0 = a;
    x1 = a + h;
    
    for (i=0; i<=step_count-1; i++) {
        sum += f(x0) + 4*f(x0 + h/2) + f(x1);
        
        x0 += h;
        x1 += h;
    }
    return (h/6)*sum;
}

double integral_monte(function f, double a, double b, unsigned step_count)
{
    double k,g,x,s = 0.0,aa;
    int i;
    k=b-a;
    for (i=0; i<=1.0*step_count; i++){
        srand(time(NULL));
        g=rand()%1000;
        x=a+g*(b-a)/1000;
        s=s+f(x);
    }
    aa=(1/(1.0*step_count))*k*s;
    return (aa);
}

void integral_runge4(dfunction f, double* x0, double* y0, double h)
{
    double k1,k2,k3,k4,d,i;
    for(i=*x0;i<=*y0;i+=h)
    {
        k1=h*f(*x0, *y0);
        k2=h*f(*x0+h/2, *y0+k1/2);
        k3=h*f(*x0+h/2, *y0+k2/2);
        k4=h*f(*x0+h, *y0+k3);
        d=(k1+2*k2+2*k3+k4)/6;
        *y0=*y0+d;
    }
    *x0=i;
}
