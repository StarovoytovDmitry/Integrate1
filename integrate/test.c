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
    return (2 * x * y);
}
int main(int argc, const char * argv[]) {
    double a,b,c;
    double x0=0,y0=2;
    a=integral_pram(f, 1, 5, 100);
    b=integral_monte(f, 1, 5, 100);
    c=integral_trap(f, 1, 5, 100);
    printf("%f %f %f \n",a,b,c);
    integral_runge4(df, &x0, &y0, 0.1);
    printf("x=%f y=%f \n",x0,y0);
    integral_runge5(df, &x0, &y0, 0.1);
    printf("x=%f y=%f \n",x0,y0);
}