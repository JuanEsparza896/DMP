#ifndef DATA_TYPES_H
#define DATA_TYPES_H

enum termostatos_t{rescvel,berendsen,andersen,bdp,nosehoover};
using termostatos = enum termostatos_t;
enum algoconstr_t{fuerzasres,rattle};
using algoconstr = enum algoconstr_t;
enum ensamble_t{nve,nvt};
using ensamble = enum ensamble_t;
enum optimizaciones_t{no_opt,vec,cel,ambas};
using optimizaciones = enum optimizaciones_t;
enum coordenadas_t{cartesianas,esfericos};
using coordenadas = enum coordenadas_t;
enum potencial_t{LJ,Yuk};
using potencial = enum potencial_t;

ushort DefNparam(potencial pot)
{
    ushort nparam;
    switch (pot)
    {
    case LJ:
        nparam=2;
        break;
    
    case Yuk:
        nparam=6;
        break;
    }
    return nparam;

}

#endif