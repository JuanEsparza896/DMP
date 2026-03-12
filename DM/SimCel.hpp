#ifndef SIMCEL_H
#define SIMCEL_H

#include <stdio.h>
#include <time.h>
#include <fstream>
#include <iostream>

#include "../MISC/CudaStruct.hpp"
#include "../MISC/Definiciones.hpp"
#include "../MISC/DataTypes.hpp"
#include "FuncComp.hpp"
#include "PotencialesInt.hpp"
#include "Constricciones.hpp"
#include "Termostatos.hpp"
#include "Restricciones.hpp"
#include "Opt.hpp"

void AceleracionesC(uint np,ushort n_esp_p,uint ncv,uint mpc,uint ncel,uint *M_int,uint *esp_p,uint *vcl,uint *np_c,uint *p_c,uint3 *mad_de_p,double &pot,double *M_param,double *pos,double *acel,double3 itcel,double3 caja,bool nconf,potencial pote)
{
    pot=0.;
    double2 fuerza;
    double3 di,dj,dij,lcaja,lcajai;
    uint3 cel;
    uint mol1,mol2;
    double dis=0.;
    fuerza=InitV2<double2,double>(0.,0.);
    lcaja=InitV3<double3,double>(caja.x,caja.y,caja.z);
    lcajai=InitV3<double3,double>(1.0/caja.x,1.0/caja.y,1.0/caja.z);
    int fac1=caja.x;
    int fac2=caja.y; fac2*=fac1;
    //Se hace solo la diagonal superior de las fuerzas
    for(int i=0;i<np;i++){
        acel[i*nd]=0.;
        acel[i*nd+1]=0.;
        acel[i*nd+2]=0.;
    }

    for(int i=0;i<np;i++){
        mol1=mad_de_p[i].x;
        di=InitV3<double3,double>(pos[i*nd],pos[1+i*nd],pos[2+i*nd]);
        cel.x=di.x*itcel.x;
        cel.y=di.y*itcel.y;
        cel.z=di.z*itcel.z;
        if(di.x==caja.x&&di.y==caja.y&&di.z==caja.z){cel.x=cel.y=cel.z=0.;}
        int ce=cel.x+fac1*cel.y+fac2*cel.z;
        if(ce<0||ce>=ncel)printf("ncel %d",ce);

        for(int cv=0;cv<ncv;cv++){
            int clv=vcl[ce*ncv+cv];
            for(int qp=0;qp<np_c[clv];qp++){
                int j=p_c[mpc*clv+qp];
                mol2=mad_de_p[j].x;
                if(i==j)continue;
                if(mol1==mol2)continue;
                dj=InitV3<double3,double>(pos[j*nd],pos[1+j*nd],pos[2+j*nd]);
                dij=Condper(di,dj,lcaja,lcajai);
                dis=Distancia2(dij);
                fuerza=Interaccion(i,j,n_esp_p,dis,0.,M_param,M_int,esp_p,pote);
                acel[i*nd]+=fuerza.x*dij.x;
                acel[i*nd+1]+=fuerza.x*dij.y;
                acel[i*nd+2]+=fuerza.x*dij.z; 
                if(nconf){
                    pot+=fuerza.y;
                }
            }
        }                
    }
    pot*=0.5;
}


void SimulacionCel(float ncp,ushort n_esp_p,ushort n_esp_m,uint nc,uint np,uint max_p_en_esp_mr,uint GdL,uint *M_int, uint *esp_p,uint *n_p_esp_m,uint *n_m_esp_mr,uint *p_en_m,uint3 *mad_de_p,double rc,double dens,double dt,double temp_d,double p_termo,double *masa,
                    double *q_rat,double *pos,double *vel,double *acel,double *constr,double *param,double *dis_p_esp_mr_rep,double3 caja,ensamble ens,algoconstr vibrante,termostatos termo,potencial pot,std::ofstream &ofasres,std::ofstream &ofasat,std::ofstream &ofaeva)
{
    bool nconf=true;
    uint64_t semillaRand[4];
    ushort ncc=(float)nc*(ncp/100.),cont=0;
    double epi=0.,eci=0.,eti=0.,tempi=0.;
    double epp=0.,ecp=0.,etp=0.,tempp=0.;
    double *M_param = new double[n_esp_p*n_esp_p*DefNparam(pot)];
    //termo nosehoover
    double xi=0.;

    InicializarSemillas(semillaRand);

    time_t ti,tf;

    short dc=rc+0.5;
    double3 tamcel,itamcel;
    uint nveccel=pow(1+2*dc,nd);

    uint ncel=CrearCeldas(caja,tamcel,itamcel);
    uint mpc=NPCel(n_esp_p,DefNparam(pot),param,tamcel);
    uint *v_cel=new uint[ncel*nveccel];
    uint *np_c=new uint[ncel];
    uint *p_c=new uint[ncel*mpc];
    
    CalculoDeCeldasVecinas(caja,v_cel,dc,nveccel);
    AsignarCeldas(np,mpc,ncel,np_c,p_c,caja,itamcel,pos);

    printf("ic\tTemperatura\tDensidad\tE_k\t\tE_p\t\tE_Total\t\tTiempo\n");
    ofasres << "ic\tTemperatura\tDensidad\tE_k\t\tE_p\t\tE_Total\t\tTiempo\n";
    M_Parametros(n_esp_p,param,M_param,pot);
    AceleracionesC(np,n_esp_p,nveccel,mpc,ncel,M_int,esp_p,v_cel,np_c,p_c,mad_de_p,epi,M_param,pos,acel,itamcel,caja,nconf,pot);
    
    for(int i=0;i<np;i++)
        ofaeva << std::fixed << std::setprecision(7) << i << "\t"<< vel[i*nd] << "\t" << vel[i*nd+1] << "\t" << vel[i*nd+2] 
        << "\t" << acel[i*nd] << "\t" << acel[i*nd+1] << "\t" << acel[i*nd+2] << std::endl;
    ofaeva.close();

    ti=clock();
    for(int ic=0;ic<=nc;ic++)
    {
        nconf=false;
        if(ic%ncc==0)nconf=true;
        //printf("ic: %d",ic);
        
        if(ens==nvt&&termo==nosehoover){ActualizarXiNoseHoover(np,GdL,esp_p,vel,masa,p_termo,temp_d,xi,0.25*dt);ReescalarNoseHoover(np,xi,0.5*dt,vel);}

        //Esta es la actualizacion de posiciones y Velocidades
        VV_Vel(np,esp_p,dt,q_rat,vel,acel,masa);

        AllPart pos[id+ip*nd] += dt*q_rat[id+ip*nd];
            
        for(int i=0;i<np;i++)CondperPos(i,pos,caja);

        ParaMoleculas
            if(vibrante==rattle)
                RattlePos(n_esp_m,max_p_en_esp_mr,n_p_esp_m,n_m_esp_mr,p_en_m,esp_p,mad_de_p,dt,constr[1],constr[0],pos,q_rat,dis_p_esp_mr_rep,masa,caja);
        

        AsignarCeldas(np,mpc,ncel,np_c,p_c,caja,itamcel,pos);

        AceleracionesC(np,n_esp_p,nveccel,mpc,ncel,M_int,esp_p,v_cel,np_c,p_c,mad_de_p,epi,M_param,pos,acel,itamcel,caja,nconf,pot);
        
        ParaMoleculas
            if(vibrante==fuerzasres)epi+=PotencialesDeRestriccion(n_esp_m,max_p_en_esp_mr,n_m_esp_mr,n_p_esp_m,p_en_m,mad_de_p,constr[0],pos,acel,dis_p_esp_mr_rep,caja);    
        
        VV_Vel(np,esp_p,dt,vel,q_rat,acel,masa);

        ParaMoleculas
            if(vibrante==rattle)RattleVel(constr[1],n_esp_m,max_p_en_esp_mr,n_p_esp_m,n_m_esp_mr,p_en_m,esp_p,mad_de_p,constr[0],pos,vel,masa,caja);

        if(ens==nvt)ApplyTermos(np,GdL,esp_p,vel,masa,termo,xi,dt,temp_d,p_termo,semillaRand);

        if(!(ic%ncc)){
            AcumularProp(ti,ic,np,GdL,esp_p,dens,epi,ecp,epp,etp,tempp,vel,masa,ofasres);
            cont++;
            ParticulasADisco(np,pos,vel,acel,ofasat);
        }

    }
    ecp/=cont;
    epp/=cont;
    etp/=cont;
    tempp/=cont;

    printf("Resultados finales (Cantidades promedio)\n");
    printf("Temperatura\tE_Potencial\tE_Cinetica\tE_Total\n");
    printf("%.6f\t%.6f\t%.6f\t%.6f\n",tempp,epp,ecp,etp);
    
}

#endif