#ifndef TERMOSTATOS_H
#define TERMOSTATOS_H

#include "../MISC/CudaStruct.hpp"
#include "../MISC/Definiciones.hpp"
#include "../MISC/DataTypes.hpp"
#include "../MISC/NumRand.h"
#include "FuncComp.hpp"
#include <math.h>
#include <stdlib.h>

void TermoReescVel(uint np,double temp_d,double temp_inst,double *vel)
{
    if(temp_inst<1e-6)printf("error en la temperatura\n");
    double lambda = sqrt(temp_d / temp_inst);
    
    AllPart
        vel[ip*nd+id]*=lambda;
}

void TermoBerendsen(uint np,double temp_d,double temp_inst,double p_termo,double dt,double *vel)
{
    if (temp_inst <= 1e-10) return;
    double lambda_2 = 1.0 + (dt / p_termo) * (temp_d / temp_inst - 1.0);
    if (lambda_2 < 0.0) lambda_2 = 0.0;
    double lambda = sqrt(lambda_2);
    //printf("lambda %lf\n",lambda);
    AllPart
        vel[ip*nd+id]*=lambda;
}

void TermoAndersen(uint np,uint *esp_p,double temp_d, double *vel,double *masas, double p_termo, double dt,uint64_t s[4]) 
{
    double prob_col = p_termo * dt; 
    
    double sigma_v;
    uint esp;
    double masa;
    
    for(int i=0;i<np;i++){
        esp=esp_p[i];
        masa=masas[esp];
        sigma_v = sqrt(k_Boltzmann * temp_d/masa);
        if(rand_d(s)<prob_col){
            vel[i*nd]=Num_Gaussiano(s,sigma_v);
            vel[i*nd+1]=Num_Gaussiano(s,sigma_v);
            vel[i*nd+2]=Num_Gaussiano(s,sigma_v);
        }
        
    }
}

void TermoBDP(uint np,uint GdL,double temp_d,double tempi,double p_termo,double dt,double *vel,uint64_t s[4])
{
    if(tempi<=0.)return;
    double fac1=exp(-dt/p_termo);
    double fac2= (1.-fac1)*(temp_d/tempi);
    if(fac1*fac2<0.)printf("Error\n");

    double n_g=Num_Gaussiano(s,1.);

    //Suma de cuadrados

    double sn_g2=S2_Gaussiano2(s,GdL);

    double alpha2=fac1+(fac2/GdL)*(sn_g2+n_g*n_g)+2.*n_g*sqrt(fac1*(fac2/GdL));

    if(alpha2<0.) alpha2=0.;

    double alpha=sqrt(alpha2);

    AllPart vel[ip*nd+id]*=alpha;
}

void ActualizarXiNoseHoover(uint np,uint GdL,uint *esp_p, double *vel,double *masa, double p_termo, double temp_d,double &xi,double f_dt)
{
    double Q=GdL*k_Boltzmann*temp_d*p_termo*p_termo;
    double eci=EnergiaCinetica(np,vel,esp_p,masa);
    double tempi=Temperatura(eci,np,GdL);

    xi+=f_dt*(GdL*k_Boltzmann/Q)*(tempi-temp_d);
}

void ReescalarNoseHoover(uint np,double xi,double f_dt,double *vel)
{
    double fa=exp(-xi*f_dt);

    AllPart 
        vel[id+ip*nd]*=fa;
}

void ApplyTermos(uint np,uint GdL,uint *esp_p, double *vel,double *masa, termostatos termo,double &xi,double dt,double temp_d,double p_termo,uint64_t s[4]){
    double eci=EnergiaCinetica(np,vel,esp_p,masa);
    double tempi=Temperatura(eci,np,GdL);
    //printf("tempi %lf\n",tempi);
    double d_dzeta;
    switch (termo)
    {
    case rescvel:
        TermoReescVel(np,temp_d,tempi,vel);
        break;
    case berendsen:
        TermoBerendsen(np,temp_d,tempi,p_termo,dt,vel);
        break;
    case andersen:
        TermoAndersen(np,esp_p,temp_d,vel,masa,p_termo,dt,s);
        break;
    case bdp:
        TermoBDP(np,GdL,temp_d,tempi,p_termo,dt,vel,s);
        break;
    case nosehoover:
        ReescalarNoseHoover(np,xi,0.5*dt,vel);
        ActualizarXiNoseHoover(np,GdL,esp_p,vel,masa,p_termo,temp_d,xi,0.5*dt);
        break;
    default:
        break;
    }
}
#endif