#ifndef __SSM0_INTERNAL_H__
#define __SSM0_INTERNAL_H__
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double rat(long int, long int);

double li22re(double, double);

double li2re(double);
double li3re(double);
double li4re(double);

double li2re2(double, double);
double li2im2(double, double);
double li3re2(double, double);

double SSm0_li2logA0(double);

void check_range(int, const char*);
void region_unimplemented(const char*);

struct PolylogConstants
{
  double epsdif;

  double li22x1c1;
  double li22x1c2;
  double li22x1c3;
  double li22x1c4;
  double li22x1c5;
  double li22x1c6;
  double li22x1c7;
  double li22x1c8;

  double xi0;
  double ccxzero;

  double Pi2;
  double Pi4;

  int LINLOGA0MAX;
  int LI2LOGA0MAX;
  int LINLOGA1MAX;
};

extern const struct PolylogConstants li_constants; 

#endif  /* __SSM0_INTERNAL_H__ */
