#ifndef CONSTRICCIONES_H
#define CONSTRICCIONES_H

#include "../MISC/CudaStruct.hpp"
#include "FuncComp.hpp"

//Este algoritmo de rattle funciona cuando las moleculas solo tienen restriccion de longitud de enlace
void RattlePos(ushort n_esp_m,uint max_p_en_esp_mr,uint *n_p_esp_m,uint *n_m_esp_mr,uint *p_en_m,uint *esp_p,uint3 *mad_de_p,double dt,double m_it,double tol,double *pos,double *q_rat,double *dis_p_esp_mr_rep,double *masa,double3 caja)
{
    double dis,dif=0.,d_max=0.,grr=0.,srij=0.;
    uint esp1,esp2;
    double im1,im2;
    uint mol=0,part=0;
    double3 dx,s,rij,dxt,cajai;
    uint max_it = m_it;
    cajai.x=1./caja.x;
    cajai.y=1./caja.y;
    cajai.z=1./caja.z;
    for(uint it=0;it<max_it;it++){
        //Este algoritmo depende en que las partículas de una misma especie esten seguidas en la lista de moleculas, 
        //si no tendríamos que hacer un loop de moleculas y consultamos su especie. 
        mol=0;
        //Desde aquí
        d_max=0.;
        for(uint e_mol=0;e_mol<n_esp_m;e_mol++){
            for(uint n_emol=0;n_emol<n_m_esp_mr[e_mol];n_emol++){
                //hasta aquí,
                // lo único que simboliza es que se recorre cada molécula del sistema
                for(uint i=0;i<n_p_esp_m[e_mol];i++){
                    part=p_en_m[mol]+i;
                    esp1=esp_p[part];
                    im1=1./masa[esp1];
                    for(uint j=part+1;j<part+mad_de_p[part].z;j++){
                        if(part==j)continue;
                        esp2=esp_p[j];
                        im2=1./masa[esp2];                    
                        dx.x=pos[part*nd]+dt*q_rat[part*nd];
                        dxt.x=pos[j*nd]+dt*q_rat[j*nd];
                        dx.y=pos[part*nd+1]+dt*q_rat[part*nd+1];
                        dxt.y=pos[j*nd+1]+dt*q_rat[j*nd+1];
                        dx.z=pos[part*nd+2]+dt*q_rat[part*nd+2];
                        dxt.z=pos[j*nd+2]+dt*q_rat[j*nd+2];
                        s=Condper(dx,dxt,caja,cajai);
                        dis=Dot(s,s);
                        dif=dis - dis_p_esp_mr_rep[e_mol*max_p_en_esp_mr*max_p_en_esp_mr+i*max_p_en_esp_mr+j-p_en_m[mol]]*dis_p_esp_mr_rep[e_mol*max_p_en_esp_mr*max_p_en_esp_mr+i*max_p_en_esp_mr+j-p_en_m[mol]];  
                        d_max=fmax(d_max,fabs(dif));
                        if(fabs(dif)>tol){
                            /*
                              Si no se cumple la constriccion se busca una corrección g para ri y rj que haga que si se cumpla
                              la forma de ri_t y rj_t son:
                              ri_t=ri+dt*(q_rati-g*rij)
                              rj_t=rj+dt*(q_ratj+g*rij)
                              para la constricción de distancia se tiene la siguiente g
                            */
                            dx.x=pos[part*nd];
                            dxt.x=pos[j*nd];
                            dx.y=pos[part*nd+1];
                            dxt.y=pos[j*nd+1];
                            dx.z=pos[part*nd+2];
                            dxt.z=pos[j*nd+2];
                            rij=Condper(dx,dxt,caja,cajai);
                            srij = Dot(s,rij);
                            //si fueran distintos tipos de constricción aquí el valor de grr cambia, esto está en progreso
                            grr=dif/(2.*dt*srij*(im1+im2));    
                            q_rat[part*nd]-=grr*im1*rij.x;
                            q_rat[j*nd]+=grr*im2*rij.x;                         
                            q_rat[part*nd+1]-=grr*im1*rij.y;
                            q_rat[j*nd+1]+=grr*im2*rij.y;                         
                            q_rat[part*nd+2]-=grr*im1*rij.z;
                            q_rat[j*nd+2]+=grr*im2*rij.z;                         
                        }             
                    }
                }
                mol++;
            }
        }
        if(d_max<tol)break;
    }
}

void RattleVel(double m_it,uint n_esp_m,uint max_p_en_esp_mr,uint *n_p_esp_m,uint *n_m_esp_mr,uint *p_en_m,uint *esp_p,uint3 *mad_de_p,double tol,double *pos,double *vel,double *masa,double3 caja)
{
    uint esp1,esp2;
    double im1,im2;
    uint mol=0,part;
    double d_max=0.,rijv,ka=0.;
    double3 di,dj,dv,dij,cajai;
    uint max_it=m_it;
    cajai.x=1./caja.x;
    cajai.y=1./caja.y;
    cajai.z=1./caja.z;
    for(uint it=0;it<max_it;it++){
        mol=0;
        d_max=0.;
        for(uint emol=0;emol<n_esp_m;emol++){
            if(n_p_esp_m[emol]>1)
            for(uint n_emol=0;n_emol<n_m_esp_mr[emol];n_emol++){
                for(uint i=0;i<n_p_esp_m[emol];i++){
                    part=p_en_m[mol]+i;
                    esp1=esp_p[part];
                    im1=1./masa[esp1];
                    for(uint j=part+1;j<part+mad_de_p[part].z;j++){
                        if(part==j)continue;
                        esp2=esp_p[j];
                        im2=1./masa[esp2];
                        di.x=pos[part*nd];
                        dj.x=pos[j*nd];
                        di.y=pos[part*nd+1];
                        dj.y=pos[j*nd+1];
                        di.z=pos[part*nd+2];
                        dj.z=pos[j*nd+2];
                        dij=Condper(di,dj,caja,cajai);
                        dv.x=(vel[part*nd]-vel[j*nd]);
                        dv.y=(vel[part*nd+1]-vel[j*nd+1]);
                        dv.z=(vel[part*nd+2]-vel[j*nd+2]);
                        rijv=Dot(dij,dv);
                        d_max=fmax(d_max,fabs(rijv));
                        if(fabs(rijv)>tol){
                            ka=rijv/(Dot(dij,dij)*(im1+im2));        
                            vel[part*nd]-=ka*im1*dij.x;
                            vel[j*nd]+=ka*im2*dij.x;
                            vel[part*nd+1]-=ka*im1*dij.y;
                            vel[j*nd+1]+=ka*im2*dij.y;
                            vel[part*nd+2]-=ka*im1*dij.z;
                            vel[j*nd+2]+=ka*im2*dij.z;
                            
                        }    
                    }
                }
                mol++;
            }
        }
        if(d_max<tol)break;
    }
}


#endif