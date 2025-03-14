#include "line_algoritms.h"
#include <cmath>

namespace line_algorithms {

    Point find_cross_section(const Line& line1, const Line& line2)
    {
        if (!lines_cross(line1, line2))
        {
            Point p{ 0.0,0.0 };
            return p;
        }
        Point p;
        p.x = line1.contained_point.x + line2.contained_point.x;
        p.y = line1.contained_point.y + line2.contained_point.y;
    }

    bool lines_cross(
        const Line& line1,
        const Line& line2)
    {
        Point dir1_norm = line1.direction;
        dir1_norm.normalize();

        Point dir2_norm = line2.direction;
        dir2_norm.normalize();

        if (dir1_norm.x == dir2_norm.x && dir1_norm.y == dir2_norm.y)
            return false;

        if (dir1_norm.x == -dir2_norm.x && dir1_norm.y == -dir2_norm.y)
            return false;

        return true;
    }

}

