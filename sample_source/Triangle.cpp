#include <iomanip>
#include "Triangle.h"
#include "Helpers.h"

Triangle::Triangle(Vec3WithColor vid1, Vec3WithColor vid2, Vec3WithColor vid3)
{
    this->v1 = vid1;
    this->v2 = vid2;
    this->v3 = vid3;
}

Triangle::Triangle(const Triangle &other)
{
    this->v1 = other.v1;
    this->v2 = other.v2;
    this->v3 = other.v3;
}

std::ostream &operator<<(std::ostream &os, const Triangle &t)
{
    os << std::fixed << std::setprecision(0) << "Triangle with vertices (" << t.v1.vertexId << ", " << t.v2.vertexId << ", " << t.v3.vertexId << ")";
    return os;
}

Vec3 Triangle::triangleNormal(const Vec3 &a,
                                    const Vec3 &b,
                                    const Vec3 &c)
{
    Vec3 ab = subtractVec3(b, a);
    Vec3 ac = subtractVec3(c, a);
    Vec3 n  = crossProductVec3(ab, ac);
    return normalizeVec3(n);
}