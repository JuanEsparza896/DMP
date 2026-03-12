#ifndef SIMVC_H
#define SIMVC_H

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

void AceleracionesCV(uint np,ushort n_esp_p,uint nmaxvec,uint *M_int,uint *nvec,uint *vec,uint *esp_p,double &pot,double rc,double *M_param,double *pos,double *acel,double3 caja,bool nconf,potencial pote)
{
    //cambiar para que en el calculo de vecinos este 
    //incluida la matriz de interaccion
    pot=0.;
    double2 fuerza;
    double3 di,dj,dij,lcaja,lcajai;
    double dis=0.;
    fuerza=InitV2<double2,double>(0.,0.);
    lcaja=InitV3<double3,double>(caja.x,caja.y,caja.z);
    lcajai=InitV3<double3,double>(1.0/caja.x,1.0/caja.y,1.0/caja.z);
    //Se hace solo la diagonal superior de las fuerzas
    for(int i=0;i<np;i++){
        acel[i*nd]=0.;
        acel[i*nd+1]=0.;
        acel[i*nd+2]=0.;
    }

    for(int i=0;i<np;i++){
        di=InitV3<double3,double>(pos[i*nd],pos[1+i*nd],pos[2+i*nd]);
        for(int k=0;k<nvec[i];k++){
            uint j=vec[i*nmaxvec+k];
            if(i>=j)continue;
            dj=InitV3<double3,double>(pos[j*nd],pos[1+j*nd],pos[2+j*nd]);
            dij=Condper(di,dj,lcaja,lcajai);
            dis=Distancia2(dij);
            fuerza=Interaccion(i,j,n_esp_p,dis,rc,M_param,M_int,esp_p,pote);
            if(fuerza.x>1000){printf("i%d,j%d,fuerza%lf\n",i,j,fuerza.x);}
            acel[i*nd]+=fuerza.x*dij.x;
            acel[i*nd+1]+=fuerza.x*dij.y;
            acel[i*nd+2]+=fuerza.x*dij.z;
            acel[j*nd]-=fuerza.x*dij.x;
            acel[j*nd+1]-=fuerza.x*dij.y;
            acel[j*nd+2]-=fuerza.x*dij.z;
            if(nconf){
                pot+=fuerza.y;
            }
        }           
    }
}


void SimulacionCV(float ncp,ushort n_esp_p,ushort n_esp_m,uint nc,uint np,uint GdL,uint max_p_en_esp_mr,uint *M_int, uint *esp_p,uint *n_p_esp_m,uint *n_m_esp_mr,uint *p_en_m,uint3 *mad_de_p,double dens,double dt,double temp_d,double p_termo,double rc,double *masa,
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

    time_t ti,tf;
    InicializarSemillas(semillaRand);
    
    uint nmaxvec=MaxNVec(n_esp_p,DefNparam(pot),np,param,rc,dens);
    uint *nvec=new uint[np];
    uint *vec=new uint[np*nmaxvec];
    short dc=rc+0.5;
    double3 tamcel,itamcel;
    uint nveccel=pow(1+2*dc,nd);

    uint ncel=CrearCeldas(caja,tamcel,itamcel);
    uint mpc=NPCel(n_esp_p,DefNparam(pot),param,tamcel);
    uint *v_cel=new uint[ncel*nveccel];
    uint *np_c=new uint[ncel];
    uint *p_c=new uint[ncel*mpc];

    //asignacion de radio de corte
    double sig=0;
    for(int i=0;i<n_esp_p;i++)if(param[i*n_esp_p]>=sig)sig=param[i*n_esp_p];
    rc*=sig;

    CalculoDeCeldasVecinas(caja,v_cel,dc,nveccel);
    AsignarCeldas(np,mpc,ncel,np_c,p_c,caja,itamcel,pos);
    VecinosConCel(np,nveccel,mpc,nmaxvec,nvec,vec,v_cel,np_c,p_c,mad_de_p,rc,pos,itamcel,caja);
    printf("nvec %d: %d\n",0,nvec[0]);

    printf("ic\tTemperatura\tDensidad\tE_k\t\tE_p\t\tE_Total\t\tTiempo\n");
    ofasres << "ic\tTemperatura\tDensidad\tE_k\t\tE_p\t\tE_Total\t\tTiempo\n";
    M_Parametros(n_esp_p,param,M_param,pot);
    AceleracionesCV(np,n_esp_p,nmaxvec,M_int,nvec,vec,esp_p,epi,rc,M_param,pos,acel,caja,nconf,pot);
    
    for(int i=0;i<np;i++)
        ofaeva << std::fixed << std::setprecision(7) << i << "\t"<< vel[i*nd] << "\t" << vel[i*nd+1] << "\t" << vel[i*nd+2] 
        << "\t" << acel[i*nd] << "\t" << acel[i*nd+1] << "\t" << acel[i*nd+2] << std::endl;
    ofaeva.close();

    
    ti=clock();
    for(int ic=0;ic<=nc;ic++)
    {
        nconf=false;
        if(ic%ncc==0)nconf=true;

        if(ens==nvt&&termo==nosehoover){ActualizarXiNoseHoover(np,GdL,esp_p,vel,masa,p_termo,temp_d,xi,0.5*dt);ReescalarNoseHoover(np,xi,0.5*dt,vel);}
        
        //Esta es la actualizacion de posiciones y Velocidades
        VV_Vel(np,esp_p,dt,q_rat,vel,acel,masa);

        ParaMoleculas
            if(vibrante==rattle)
                RattlePos(n_esp_m,max_p_en_esp_mr,n_p_esp_m,n_m_esp_mr,p_en_m,esp_p,mad_de_p,dt,constr[1],constr[0],pos,q_rat,dis_p_esp_mr_rep,masa,caja);
        
        AllPart pos[id+ip*nd] += dt*q_rat[id+ip*nd];
    
        for(int i=0;i<np;i++)CondperPos(i,pos,caja);


        if(ic%20==0){
            AsignarCeldas(np,mpc,ncel,np_c,p_c,caja,itamcel,pos);
            VecinosConCel(np,nveccel,mpc,nmaxvec,nvec,vec,v_cel,np_c,p_c,mad_de_p,rc,pos,itamcel,caja);
        }
        AceleracionesCV(np,n_esp_p,nmaxvec,M_int,nvec,vec,esp_p,epi,rc,M_param,pos,acel,caja,nconf,pot);
    
        ParaMoleculas
            if(vibrante==fuerzasres)epi+=PotencialesDeRestriccion(n_esp_m,max_p_en_esp_mr,n_m_esp_mr,n_p_esp_m,p_en_m,mad_de_p,constr[0],pos,acel,dis_p_esp_mr_rep,caja);    
        
        VV_Vel(np,esp_p,dt,vel,q_rat,acel,masa);

        //AllPart if(fabs(vel[ip*nd+id])>1000.)printf("ic %d,i %d,%lf\n",ic,ip,vel[ip*nd+id]);

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