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

/* extern const struct PolylogConstants li_constants;  */

static const struct PolylogConstants li_constants =
  {
    .epsdif = 5e-14,            /* WTFepsdif */
    /* Li22 x->1 */
    .li22x1c1 = -2.,
    .li22x1c2 = 2.,
    .li22x1c3 = 0.5,
    .li22x1c4 = -0.166666666666666666667,
    .li22x1c5 = 0.333333333333333333333,
    .li22x1c6 = 1.64493406684822643647,  // Pi^2/6
    .li22x1c7 = -2.40411380631918857080, // -2*Zeta(3)
    .li22x1c8 = 2.16464646742227638303,  // Pi^4/45

    /* f1/f2 splitting for Li22 */
    .xi0 = 1,
    .ccxzero = 0.36787944117144232159552377016146, /* exp(-xi0); */

    .Pi2 = 9.8696044010893586188, /* Pi^2 */
    .Pi4 = 97.409091034002437236, /* Pi^4 */

    /* Recursion limits */
    .LI2LOGA0MAX = 8,
    .LINLOGA0MAX = 18,
    .LINLOGA1MAX = 8
  };

#endif  /* __SSM0_INTERNAL_H__ */
