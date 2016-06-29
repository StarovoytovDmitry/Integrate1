//
//  test.c
//  integrate
//
//  Created by Дмитрий on 15.06.16.
//  Copyright © 2016 Дмитрий. All rights reserved.
//

/*!
 \file
 \brief Файл с тестами
 \author Dmitry
 \version 1.0
 Данный файл содержит в себе тесты основных
 функций
 */

#include <stdio.h>
#include "Integrate.h"
double f (double x) {
    return 2 * x;
}
double f1 (double x) {
    return x/(x*x*x*x+1);//сходится к pi/4=0.7853981 при a=0
}
double f2 (double x) {
    return 1/(x*x);//сходится к 1 при a=1
}
double df (double x, double y) {
    return (2 * (x*x + y));
}
int main(int argc, const char * argv[]) {
    double a,b,c,d;
    double a1=0;
    double x=0,y=0;
    // тест для мнк
    double x1[5]={-3,-1,0,1,3};
    double y1[5]={-4,-0.8,1.6,2.3,1.5};
    // тест для мнк 2
    //double x2[5]={-2,0,1,2,4};
    //double y2[5]={0.5,1,1.5,2,3};
    double aa[3]={0,0,0};
    a=integral_pram(f, 1, 5, 100);
    b=integral_monte(f, 1, 5, 100);
    c=integral_trap(f, 1, 5, 100);
    d=integral_simp(f, 1, 5, 100);
    printf("%f %f %f %f \n",a,b,c,d);
    x=0,y=0;
    integral_runge4(df, 0, 1, 1, &x, &y, 0.1);
    printf("runge4 x=%f y=%f \n",x,y);
    x=0,y=0;
    integral_runge5(df, 0, 1, 1, &x, &y, 0.1);
    printf("runge5 x=%f y=%f \n",x,y);
    x=0,y=0;
    integral_eiler(df, 0, 1, 1, &x, &y, 0.1);
    printf("eiler x=%f y=%f \n",x,y);
    a1=integral_pram_inf(f2, 1, 0.1, 0.001);
    printf("%f\n",a1);
    a1=integral_trap_inf(f2, 1, 0.1, 0.001);
    printf("%f\n",a1);
    a1=integral_simp_inf(f2, 1, 0.1, 0.001);
    printf("%f\n",a1);
    printf("mnk\n");
    mnk(2, 5, x1, y1, aa);
    for(int i=2;i>=0;i--) printf("%f \n", aa[i]);
}