#include "Point.h"
#include <cmath>

void Point::normalize()
{
    double d = sqrt(x * x + y * y);
    if (d == 0)
        return;
    x /= d;
    y /= d;
}


