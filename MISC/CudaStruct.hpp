#ifndef CUDA_DEFINICIONES_H
#define CUDA_DEFINICIONES_H

using uint = unsigned int;
using ushort = unsigned short;
using uchar = unsigned char;

typedef struct double4_t
{
    double x;
    double y;
    double z;
    double w;
}double4;

typedef struct double3_t
{
    double x;
    double y;
    double z;
}double3;

typedef struct double2_t
{
    double x;
    double y;
}double2;

typedef struct uint4_t
{
    uint x;
    uint y;
    uint z;
    uint w;
}uint4;

typedef struct uint3
{
    uint x;
    uint y;
    uint z;
}uint3;

template <typename T,typename K>
T InitV3(K vx,K vy,K vz)
{
    T var;
    var.x=vx;
    var.y=vy;
    var.z=vz;
    return var;
}

template <typename T,typename K>
T InitV4(K vx,K vy,K vz,K vw)
{
    T var;
    var.x=vx;
    var.y=vy;
    var.z=vz;
    var.w=vw;
    return var;
}

template <typename T,typename K>
T InitV2(K vx,K vy)
{
    T var;
    var.x=vx;
    var.y=vy;
    return var;
}

#endif