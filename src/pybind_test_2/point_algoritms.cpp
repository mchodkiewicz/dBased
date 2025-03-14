#include "point_algoritms.h"

namespace point_algorithms {

    Point sum(
        const Point& p1,
        const Point& p2)
    {
        Point p;
        p.x = p1.x + p2.x;
        p.y = p1.y + p2.y;
        return p;
    }

    double product(const Point& p1, const Point& p2)
    {
        return p1.x * p2.x + p1.y * p2.y;
    }

}
