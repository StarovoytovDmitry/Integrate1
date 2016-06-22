//
//  test.c
//  integrate
//
//  Created by Дмитрий on 15.06.16.
//  Copyright © 2016 Дмитрий. All rights reserved.
//

#include <stdio.h>
#include "Header.h"
double f (double x) {
    return 2 * x;
}
double df (double x, double y) {
    return (2 * (x*x + y));
}
int main(int argc, const char * argv[]) {
    double a,b,c;
    double x0=0,y0=1;
    double x=0,y=0;
    a=integral_pram(f, 1, 5, 100);
    b=integral_monte(f, 1, 5, 100);
    c=integral_trap(f, 1, 5, 100);
    printf("%f %f %f \n",a,b,c);
    x0=0,y0=1;
    x=0,y=0;
    integral_runge4(df, x0, 1, y0, &x, &y, 0.1);
    printf("x=%f y=%f \n",x,y);
    x0=0,y0=1;
    x=0,y=0;
    integral_runge5(df, x0, 1, y0, &x, &y, 0.1);
    printf("x=%f y=%f \n",x,y);
    x0=0,y0=1;
    x=0,y=0;
    integral_eiler(df, x0, 1, y0, &x, &y, 0.1);
    printf("x=%f y=%f \n",x,y);
}