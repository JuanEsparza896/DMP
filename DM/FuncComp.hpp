#ifndef FUNCCOMP_H
#define FUNCCOMP_H

#include "../MISC/Definiciones.hpp"
#include "../MISC/CudaStruct.hpp"
#include <math.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <iomanip>

double3 Condper(double3 p1,double3 p2,double3 caja, double3 cajai)
{
    double3 dx;
    dx.x=p1.x-p2.x;
    dx.y=p1.y-p2.y;
    dx.z=p1.z-p2.z;
    dx.x-=rint(dx.x*cajai.x)*caja.x;
    dx.y-=rint(dx.y*cajai.y)*caja.y;
    dx.z-=rint(dx.z*cajai.z)*caja.z;

    return dx;
}

double3 CondperS(double3 dx,double3 caja,double3 cajai)
{
    dx.x-=rint(dx.x*cajai.x)*caja.x;
    dx.y-=rint(dx.y*cajai.y)*caja.y;
    dx.z-=rint(dx.z*cajai.z)*caja.z;

    return dx;
}

void CondperPos(uint ip,double *pos,double3 caja)
{
    if(pos[ip*nd] > caja.x){pos[ip*nd] -= caja.x;}
    if(pos[ip*nd] < 0){pos[ip*nd] += caja.x;}
    if(pos[1+ip*nd] > caja.y){pos[1+ip*nd] -= caja.y;}
    if(pos[1+ip*nd] < 0){pos[1+ip*nd] += caja.y;}
    if(pos[2+ip*nd] > caja.z){pos[2+ip*nd] -= caja.z;}
    if(pos[2+ip*nd] < 0){pos[2+ip*nd] += caja.z;}
}
double Dot(double3 dx,double3 dy)
{  
    return dx.x*dy.x+dx.y*dy.y+dx.z*dy.z;

}

double Distancia2(double3 var)
{
    return (var.x*var.x+var.y*var.y+var.z*var.z);
}
double Distancia2CP(double3 p1,double3 p2,double3 caja, double3 cajai)
{
    double3 var=Condper(p1,p2,caja,cajai);
    return (var.x*var.x+var.y*var.y+var.z*var.z);
}

void VV_Vel(uint np,uint *esp_p,double dt,double *arr,double *vel,double *acel,double *masa)
{
    double im;
    uint esp;
    for(int ip=0;ip<np;ip++){
        esp=esp_p[ip];
        im=1./masa[esp];
        for(int id=0;id<nd;id++)
            arr[id+ip*nd] =(vel[id+ip*nd] + im*acel[id+ip*nd]*dt*0.5);
    }
        
}

void CopiarArrMov(uint np, double *arr1,double *arr2)
{
    AllPart
        arr1[ip*nd+id]=arr2[ip*nd+id];
}

double EnergiaCinetica(uint np,double *vel,uint *esp_p,double *masa)
{
    double ec=0.0;
    uint esp;
    double m;
   for(int ip=0;ip<np;ip++){
        esp=esp_p[ip];
        m=masa[esp];
        for(int id=0;id<nd;id++)
        ec += m*vel[ip*nd+id] * vel[ip*nd+id];
    }
        
    return 0.5*ec; 
}

double Temperatura(double E_Kin,uint np,uint GdL)
{
    return (2.0*E_Kin)/((double)GdL*k_Boltzmann);
}

void ParticulasADisco(uint np,double *pos,double *vel,double *acel,std::ofstream &ofapos)
{
    for(uint i=0;i<np;i++){
        ofapos << std::fixed << std::setprecision(7) << i << "\t"
        << std::setw(10) << pos[i*nd] << " " << std::setw(10) << pos[i*nd+1] << " " << std::setw(10) << pos[i*nd+2] << "\t"
        << std::setw(10) << vel[i*nd] << " " << std::setw(10) << vel[i*nd+1] << " " << std::setw(10) << vel[i*nd+2] << "\t"
        << std::setw(10) << acel[i*nd] << " " << std::setw(10) << acel[i*nd+1] << " " << std::setw(10) << acel[i*nd+2] << std::endl;
    }
}

void AcumularProp(time_t ti,int ic,uint np,uint GdL,uint *esp_p,double dens,double epi,double &ecp,double &epp,double &etp,double &tempp,double *vel,double *masa,std::ofstream &ofasres)
{            
    time_t tf=clock();
    double dtt =((double)(tf - ti))/CLOCKS_PER_SEC;
    double eci=EnergiaCinetica(np,vel,esp_p,masa);
    epi/=static_cast<double>(np);
    double tempi=Temperatura(eci,np,GdL);
    eci/=np;
    double eti=eci+epi;
    
            
    ecp+=eci;
    epp+=epi;
    etp+=eti;
    tempp+=tempi;
           
    std::cout << std::fixed << std::setprecision(7) << ic << "\t"<< tempi << "\t" << dens << "\t"<< eci<<
    "\t"<< epi << "\t"<< eti << "\t" << dtt <<std::endl;

    ofasres << std::fixed << std::setprecision(7) << ic << "\t"<< tempi << "\t" << dens << "\t"<< eci<<
    "\t"<< epi << "\t"<< eti << "\t" << dtt <<std::endl;
}
#endif