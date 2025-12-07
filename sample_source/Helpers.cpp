#include <iostream>
#include <cmath>
#include "Helpers.h"
#include "Translation.h"
#include "Matrix4.h"
#include "Scaling.h"
#include "Rotation.h"
#include "Camera.h"

/*
 * Calculate cross product of vec3 a, vec3 b and return resulting vec3.
 */
Vec3 crossProductVec3(Vec3 a, Vec3 b)
{
    return Vec3(a.y * b.z - b.y * a.z, b.x * a.z - a.x * b.z, a.x * b.y - b.x * a.y);
}

/*
 * Calculate dot product of vec3 a, vec3 b and return resulting value.
 */
double dotProductVec3(Vec3 a, Vec3 b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

/*
 * Find length (|v|) of vec3 v.
 */
double magnitudeOfVec3(Vec3 v)
{
    return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

/*
 * Normalize the vec3 to make it unit vec3.
 */
Vec3 normalizeVec3(Vec3 v)
{
    double d = magnitudeOfVec3(v);
    return Vec3(v.x / d, v.y / d, v.z / d);
}

/*
 * Return -v (inverse of vec3 v)
 */
Vec3 inverseVec3(Vec3 v)
{
    return Vec3(-v.x, -v.y, -v.z);
}

/*
 * Add vec3 a to vec3 b and return resulting vec3 (a+b).
 */
Vec3 addVec3(Vec3 a, Vec3 b)
{
    return Vec3(a.x + b.x, a.y + b.y, a.z + b.z);
}

/*
 * Subtract vec3 b from vec3 a and return resulting vec3 (a-b).
 */
Vec3 subtractVec3(Vec3 a, Vec3 b)
{
    return Vec3(a.x - b.x, a.y - b.y, a.z - b.z);
}

/*
 * Multiply each element of vec3 with scalar.
 */
Vec3 multiplyVec3WithScalar(Vec3 v, double c)
{
    return Vec3(v.x * c, v.y * c, v.z * c);
}

/*
 * Prints elements in a vec3. Can be used for debugging purposes.
 */
void printVec3(Vec3 v)
{
    std::cout << "(" << v.x << "," << v.y << "," << v.z << ")" << std::endl;
}

/*
 * Check whether vec3 a and vec3 b are equal.
 * In case of equality, returns 1.
 * Otherwise, returns 0.
 */
int areEqualVec3(Vec3 a, Vec3 b)
{

    /* if x difference, y difference and z difference is smaller than threshold, then they are equal */
    if ((ABS((a.x - b.x)) < EPSILON) && (ABS((a.y - b.y)) < EPSILON) && (ABS((a.z - b.z)) < EPSILON))
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

/*
 * Returns an identity matrix (values on the diagonal are 1, others are 0).
 */
Matrix4 getIdentityMatrix()
{
    Matrix4 result;

    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            if (i == j)
            {
                result.values[i][j] = 1.0;
            }
            else
            {
                result.values[i][j] = 0.0;
            }
        }
    }

    return result;
}

/*
 * Multiply matrices m1 (Matrix4) and m2 (Matrix4) and return the result matrix r (Matrix4).
 */
Matrix4 multiplyMatrixWithMatrix(Matrix4 m1, Matrix4 m2)
{
    Matrix4 result;
    double total;

    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            total = 0;
            for (int k = 0; k < 4; k++)
            {
                total += m1.values[i][k] * m2.values[k][j];
            }

            result.values[i][j] = total;
        }
    }

    return result;
}

/*
 * Multiply matrix m (Matrix4) with vector v (Vec4WithColor) and store the result in vector r (Vec4WithColor).
 */
Vec4WithColor multiplyMatrixWithVec4WithColor(Matrix4 m, Vec4WithColor v)
{
    double values[4];
    double total;

    for (int i = 0; i < 4; i++)
    {
        total = 0;
        for (int j = 0; j < 4; j++)
        {
            total += m.values[i][j] * v.getNthComponent(j);
        }
        values[i] = total;
    }

    return Vec4WithColor(values[0], values[1], values[2], values[3], v.color);
}

// rotation matrix
// scaling matrix
// translation matrix

Matrix4 translationMatrix(Translation *t)
{
    Matrix4 result;
    result = getIdentityMatrix();
    result.values[0][3] = t->tx;
    result.values[1][3] = t->ty;
    result.values[2][3] = t->tz;
    return result;
}

//alternative method fro m slides
Matrix4 rotationMatrix(Rotation *r)
{
    // Axis from spec: rotation around line through (0,0,0) and (ux,uy,uz)
    Vec3 axis(r->ux, r->uy, r->uz);

    // 1) u: normalize axis (dönme yönü (ux,uy,uz) – spec ile uyumlu)
    Vec3 u = axis.unit();
    double ax = std::abs(u.x);
    double ay = std::abs(u.y);
    double az = std::abs(u.z);

    // 2) u ile paralel olmayan v seç (slaytlardaki “Alternative method”)
    Vec3 v;
    if (ax <= ay && ax <= az) {
        // x en küçük → x=0
        v = Vec3(0.0, -u.z, u.y);
    } else if (ay <= ax && ay <= az) {
        // y en küçük → y=0
        v = Vec3(-u.z, 0.0, u.x);
    } else {
        // z en küçük → z=0
        v = Vec3(-u.y, u.x, 0.0);
    }
    v = v.unit();

    // 3) w = u × v   (burada * cross ise böyle, değilse cross(u,v) kullan)
    Vec3 w = crossProductVec3(u, v).unit();

    // 4) M matrisi (xyz → uvw), tam 4×4
    Matrix4 M;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            M.values[i][j] = 0.0;

    // satırlar: u, v, w, (0 0 0 1)
    M.values[0][0] = u.x; M.values[0][1] = u.y; M.values[0][2] = u.z; M.values[0][3] = 0.0;
    M.values[1][0] = v.x; M.values[1][1] = v.y; M.values[1][2] = v.z; M.values[1][3] = 0.0;
    M.values[2][0] = w.x; M.values[2][1] = w.y; M.values[2][2] = w.z; M.values[2][3] = 0.0;
    M.values[3][0] = 0.0; M.values[3][1] = 0.0; M.values[3][2] = 0.0; M.values[3][3] = 1.0;

    // 5) M^{-1} = M^T (u,v,w ortonormal → ONB)
    Matrix4 M_inv = M.transpose();

    // 6) Rx(α): x-ekseni etrafında α derece CCW dönme
    double rad = r->angle * M_PI / 180.0;
    double c = std::cos(rad);
    double s = std::sin(rad);

    Matrix4 Rx;
    Rx = getIdentityMatrix();

    Rx.values[1][1] =  c;  
    Rx.values[1][2] = -s;
    Rx.values[2][1] =  s;
    Rx.values[2][2] =  c;

    // 7) Slayt formülü: R = M^{-1} * Rx * M
    // Bu, ekseni (ux,uy,uz) olan ve (0,0,0)'dan geçen doğru etrafında
    // (ux,uy,uz) yönü boyunca CCW dönmeyi verir (spec ile tam uyumlu).
    Matrix4 R = multiplyMatrixWithMatrix(multiplyMatrixWithMatrix(M_inv, Rx), M);

    return R;
}

Matrix4 scalingMatrix(Scaling *s)
{
    Matrix4 result;
    result = getIdentityMatrix();
    result.values[0][0] = s->sx;
    result.values[1][1] = s->sy;
    result.values[2][2] = s->sz;
    return result;
}


//kamera projection matrixleri
Matrix4 cameraTransformationMatrix(Camera *cam)
{
    Matrix4 result, T, R;
    T = getIdentityMatrix();
    R = getIdentityMatrix();

    Vec3 eye = inverseVec3(cam->position); 

    T.values[0][3] = eye.x;
    T.values[1][3] = eye.y;
    T.values[2][3] = eye.z;

    R.values[0][0] = cam->u.x;
    R.values[0][1] = cam->u.y;  
    R.values[0][2] = cam->u.z;
    R.values[1][0] = cam->v.x;
    R.values[1][1] = cam->v.y;
    R.values[1][2] = cam->v.z;
    R.values[2][0] = cam->w.x;
    R.values[2][1] = cam->w.y;
    R.values[2][2] = cam->w.z;        
    
    result = multiplyMatrixWithMatrix(R, T);

    return result;
}


Matrix4 perspectiveMatrix(Camera *c) {
    Matrix4 result = getIdentityMatrix();
    double r = c->right;
    double l = c->left;
    double t = c->top;
    double b = c->bottom;
    double n = c->near;
    double f = c->far;

    result.values[0][0] = 2 * n / (r - l);
    result.values[0][2] = (r + l) / (r - l);

    result.values[1][1] = 2 * n / (t - b);
    result.values[1][2] = (t + b) / (t - b);

    result.values[2][2] = -(f + n) / (f - n);
    result.values[2][3] = -2 * f * n / (f - n);

    result.values[3][2] = -1.0;
    result.values[3][3] = 0.0;

    return result;
}


Matrix4 orthographicMatrix(Camera *cam)
{
    Matrix4 result;

    result = getIdentityMatrix();

    result.values[0][0] = 2.0 / (cam->right - cam->left);
    result.values[1][1] = 2.0 / (cam->top - cam->bottom);
    result.values[2][2] = -2.0 / (cam->far - cam->near);
    result.values[0][3] = -(cam->right + cam->left) / (cam->right - cam->left);
    result.values[1][3] = -(cam->top + cam->bottom) / (cam->top - cam->bottom);
    result.values[2][3] = -(cam->far + cam->near) / (cam->far - cam->near);

    return result;
}

Matrix4 viewportMatrix(Camera *c) {
    Matrix4 result = getIdentityMatrix();
    double nx = c->horRes;
    double ny = c->verRes;

    result.values[0][0] = nx / 2.0;
    result.values[0][3] = (nx - 1) / 2.0;

    result.values[1][1] = ny / 2.0;
    result.values[1][3] = (ny - 1) / 2.0;

    result.values[2][2] = 0.5;
    result.values[2][3] = 0.5;

    return result;
}