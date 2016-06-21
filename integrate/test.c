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
    double a,b;
    double x0=0,y0=1;
    a=integral_pram(f, 1, 5, 100);
    b=integral_monte(f, 1, 5, 100);
    printf("%f %f \n",a,b);
    integral_runge4(df, &x0, &y0, 0.1);
    //printf("x=%f y=%f \n",x0,y0);//izm
}