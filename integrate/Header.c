//
//  main.c
//  integrate
//
//  Created by Дмитрий on 15.06.16.
//  Copyright © 2016 Дмитрий. All rights reserved.
//

/*!
 \file
 \brief Файл с реализацией функций
 \author Dmitry
 \version 1.0
 Данный файл содержит в себе реализацию основных
 функций, используемых в демонстрационной программе
 */

#include "Header.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

double integral_pram(function f, double a, double b, unsigned step_count)
{
    double S=.0;
    if (a>b) {double c=b;b=a;a=c;}
    if (0 == step_count) return S;
    const double h=(b-a)/(1.0*step_count);
    for(unsigned int i=0;i<step_count-1;i++)
    {
        double x=a+i*h;
        S=S+f(x);
    }
    S=h*S;
    return S;
}
double integral_trap(function f, double a, double b, unsigned step_count)
{
    if (a>b) {double c=b;b=a;a=c;}
    double sum = .0;
    if (0 == step_count) return sum;
    
    const double step = (b - a)/(1.0*step_count);
    for (unsigned int i = 1;i<step_count;++i) {
        sum += f (a + i * step);
    }
    sum += (f(a) + f(b)) / 2;
    sum *= step;
    return sum;
}
double integral_simp(function f, double a, double b, unsigned step_count)
{
    if (a>b) {double c=b;b=a;a=c;}
    double sum,x0,x1;
    if (0 == step_count) return sum;
    const double h = (b - a)/(1.0*step_count);
    sum = 0;
    x0 = a;
    x1 = a + h;
    
    for (unsigned int i=0; i<=step_count-1; i++) {
        sum += f(x0) + 4*f(x0 + h/2) + f(x1);
        
        x0 += h;
        x1 += h;
    }
    return (h/6)*sum;
}
double integral_monte(function f, double a, double b, unsigned step_count)
{
    if (a>b) {double c=b;b=a;a=c;}
    double g,x,s = .0;
    const double k=b-a;
    srand(time(NULL));
    for (int i=0; i<=step_count; i++){
        g=rand()%1000;
        x=a+g*(b-a)/1000;
        s=s+f(x);
    }
    return ((1.0/step_count)*k*s);
}
void integral_runge4(dfunction f, double x0, double x1, double y0, double* x, double* y, double h)
{
    if (x0>x1) {double c=x1;x1=x0;x0=c;}
    const int n=(x1-x0)/h;
    int i=0;
    for(i=x0;i<n;i++)
    {
        const double xi=i*h;
        const double k1=h*f(xi, y0);
        const double k2=h*f(xi+h/2, y0+k1/2);
        const double k3=h*f(xi+h/2, y0+k2/2);
        const double k4=h*f(xi+h, y0+k3);
        const double d=(k1+2*k2+2*k3+k4)/6;
        y0=y0+d;
    }
    *x=i*h;
    *y=y0;
}
void integral_runge5(dfunction f, double x0, double x1, double y0, double* x, double* y, double h)
{
    if (x0>x1) {double c=x1;x1=x0;x0=c;}
    const int n=(x1-x0)/h;
    int i=0;
    for(i=x0;i<n;i++)
    {
        const double xi=i*h;
        const double k1=h*f(xi, y0)/3;
        const double k2=h*f(xi+h/3, y0+k1)/3;
        const double k3=h*f(xi+h/3, y0+k1/2+k2/2)/3;
        const double k4=h*f(xi+h/2, y0+3*k1/8+9*k3/8)/3;
        const double k5=h*f(xi+h, y0+3*k1/2-9*k3/2+6*k4)/3;
        const double d=(k1+4*k4+k5)/2;
        y0=y0+d;
    }
    *x=i*h;
    *y=y0;
}
void integral_eiler(dfunction f, double x0, double x1, double y0, double* x, double* y, double h)
{
    if (x0>x1) {double c=x1;x1=x0;x0=c;}
    const int n=(x1-x0)/h;
    for (int i=x0;i<n;i++)
    {
        y0 += h * f(x0, y0);
        x0 += h;
    }
    *x=x0;
    *y=y0;
}
double integral_pram_inf(function f, double a, double b, double h, double eps)
{
    double dS=eps,S=.0,s1=.0,a1=a;
    if (eps == 0) return S;
    int i=1;
    while ((dS>=eps)&&(a1>=a))
    {
        double x=a+i*h;
        s1=S;
        S=S+f(x);
        dS=fabs(S-s1);
        a1=a1+h;
        i++;
    }
    S=h*S;
    return S;
}
double integral_trap_inf(function f, double a, double b, double h, double eps)
{
    double sum = .0,a1=a,dS=eps, s1=.0;
    if (eps == 0) return sum;
    int i=1;
    double x = .0;
    while ((dS>=eps)&&(a1>=a))
    {
        x=a+i*h;
        s1=sum;
        sum += f (x);
        dS=fabs(sum-s1);
        a1=a1+h;
        i++;
    }
    sum += (f(a) + f(x)) / 2;
    sum *= h;
    return sum;
}
double integral_simp_inf(function f, double a, double b, double h, double eps)
{
    double sum = 0,x0 = a,x1,dS=eps, s1=.0;
    if (eps == 0) return sum;
    x1 = a + h;
    while (dS>=eps)
    {
        s1=sum;
        sum += f(x0) + 4*f(x0 + h/2) + f(x1);
        dS=fabs(sum-s1);
        x0 += h;
        x1 += h;
    }
    return (h/6)*sum;
}
void mnk(int s, int n, double *x, double *y, double *a_res[])
{
    double sums[n][n];
    double b[n],a[n];
    //упорядочиваем точки по возрастанию абсцисс
    for(int i=0; i<n; i++)
   	{
        for(int j=i;j>=1;j--)
            if(x[j]<x[j-1])
            {
                const double t=x[j-1]; x[j-1]=x[j]; x[j]=t;
                const double t1=y[j-1]; y[j-1]=y[j]; y[j]=t1;
            }
   	}
    //Заполняем коэффициенты системы уравнений
   	for(int i=0; i<s+1; i++)
    {
        for(int j=0; j<s+1; j++)
        {
            sums[i][j] = 0;
            for(int k=0; k<n; k++)
                sums[i][j]+=power1(x[k],i+j);
        }
   	}
    //Заполняем столбец свободных членов
    for(int i=0; i<s+1; i++)
    {
        b[i]=0;
        for(int k=0; k<n; k++)
            b[i] +=power1(x[k],i)*y[k];
    }
    //Применяем метод Гаусса для приведения матрицы системы к треугольному виду
    for(int k=0; k<s+1; k++)
    {
        for(int i=k+1; i<s+1; i++)
        {
            const double M=sums[i][k]/sums[k][k];
            for(int j=k; j<s+1; j++) sums[i][j] -= M * sums[k][j];
            b[i] -= M*b[k];
        }
    }
    //Вычисляем коэффициенты аппроксимирующего полинома
    for(int i=s;i>=0;i--)
    {
        double s1=0;
        for(int j=i; j<s+1; j++) s1+=sums[i][j]*a[j];
            a[i] = (b[i]-s1)/sums[i][i];
            printf("%f\n", a[i]);
    }
    for(int i=s;i>=0;i--) *a_res=&a[i];
}
// Вспомогательная функция возведения в степень
int power1(int t, int k)
{
    int res = 1;
    while (k)
    {
        if (k & 1)
            res *= t;
        t *= t;
        k >>= 1;
    }
    return res;
}