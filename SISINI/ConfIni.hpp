#ifndef CONF_INI_H
#define CONF_INI_H

#include "../MISC/CudaStruct.hpp"
#include "../MISC/DataTypes.hpp"
#include "../MISC/Definiciones.hpp"
#include <math.h>
#include <iostream>
#include <fstream>

double3 CreandoCeldaMinima(uint n_esp_m,double3 *pos_respecto_p_central,uint max_p_en_esp_mr,
                           double *param, ushort nparam,uint *esp_p_en_esp_mr,uint *n_p_esp_m)
{
    double3 lon_p,lon_n;
    double3 comparador_pos,comparador_neg;
    int k=0;
    
    lon_p=InitV3<double3,double>(0.0,0.0,0.0);
    lon_n=InitV3<double3,double>(0.0,0.0,0.0);
    
    
    for(int i=0;i<n_esp_m;i++){
        for(int j=0;j<n_p_esp_m[i];j++){
            k=esp_p_en_esp_mr[i*max_p_en_esp_mr+j];
            //estas expresiones asumen que el primer parametro para cada especie atomica siempre es su diametro
            comparador_pos=InitV3<double3,double>(pos_respecto_p_central[max_p_en_esp_mr*i+j].x+0.5*param[k*nparam],
                                                pos_respecto_p_central[max_p_en_esp_mr*i+j].y+0.5*param[k*nparam],
                                                pos_respecto_p_central[max_p_en_esp_mr*i+j].z+0.5*param[k*nparam]);

            comparador_neg=InitV3<double3,double>(pos_respecto_p_central[max_p_en_esp_mr*i+j].x-0.5*param[k*nparam],
                                                pos_respecto_p_central[max_p_en_esp_mr*i+j].y-0.5*param[k*nparam],
                                                pos_respecto_p_central[max_p_en_esp_mr*i+j].z-0.5*param[k*nparam]);
            if(comparador_pos.x>=lon_p.x)lon_p.x=comparador_pos.x;
            if(comparador_pos.y>=lon_p.y)lon_p.y=comparador_pos.y;
            if(comparador_pos.z>=lon_p.z)lon_p.z=comparador_pos.z;
            if(comparador_neg.x<=lon_n.x)lon_n.x=comparador_neg.x;
            if(comparador_neg.y<=lon_n.y)lon_n.y=comparador_neg.y;
            if(comparador_neg.z<=lon_n.z)lon_n.z=comparador_neg.z;
        }
    }
    return InitV3<double3,double>(lon_p.x-lon_n.x,lon_p.y-lon_n.y,lon_p.z-lon_n.z);
}

void CentrarMoleculas(double *centrar_m,ushort n_esp_m,uint *n_p_esp_m,
                      uint *esp_p_en_esp_mr,uint max_p_en_esp_mr,ushort nparam,
                      double *param,double3 *pos_respecto_p_central)
{

    double3 lon_n;
    double3 comparador_neg;
    int k=0;    

    for(int i=0;i<n_esp_m;i++){
        lon_n=InitV3<double3,double>(0.0,0.0,0.0);
        for(int j=0;j<n_p_esp_m[i];j++){
            k=esp_p_en_esp_mr[i*max_p_en_esp_mr+j];
            //estas expresiones asumen que el primer parametro para cada especie atomica siempre es su diametro
            comparador_neg=InitV3<double3,double>(pos_respecto_p_central[max_p_en_esp_mr*i+j].x-0.5*param[k*nparam],
                                                              pos_respecto_p_central[max_p_en_esp_mr*i+j].y-0.5*param[k*nparam],
                                                              pos_respecto_p_central[max_p_en_esp_mr*i+j].z-0.5*param[k*nparam]);
            if(comparador_neg.x<=lon_n.x)lon_n.x=comparador_neg.x;
            if(comparador_neg.y<=lon_n.y)lon_n.y=comparador_neg.y;
            if(comparador_neg.z<=lon_n.z)lon_n.z=comparador_neg.z;
        }
        centrar_m[nd*i]=-lon_n.x;
        centrar_m[nd*i+1]=-lon_n.y;
        centrar_m[nd*i+2]=-lon_n.z;
    }
}


void ConfiguracionCubica(ushort n_esp_m,uint *n_m_esp_mr,uint *n_p_esp_m,uint *p_en_m,double *pos,
                         double3 *pos_respecto_p_central,uint max_p_en_esp_mr,double3 &caja,double *centrar_m,
                         double3 celda_minima,double densidad,std::ofstream &ofapin,uint *esp_p_en_esp_mr,
                         uint nm,uint *esp_p,uint3 *m_de_p)
{
    /********************************************/
    uint3 particulas_por_lado;
    double3 cel;
    uint *moleculas_de_especie_acumuladas,contp=0;
    int k=0,part=0;
    /********************************************/
    particulas_por_lado.x=pow(nm,1.0/nd)+0.5;
    particulas_por_lado.y=nm/particulas_por_lado.x;
    particulas_por_lado.y=pow(particulas_por_lado.y,1.0/(nd-1))+0.5;
    particulas_por_lado.z=nm/(particulas_por_lado.x*particulas_por_lado.y);
    if(particulas_por_lado.z*particulas_por_lado.x*particulas_por_lado.y<nm)particulas_por_lado.z++;
    printf("\nMoleculas: %d\nMoleculas en cada direccion(inicialmente): (%d,%d,%d)\nMoleculas que caben en la caja de sumulacion %d\n",nm,particulas_por_lado.x,particulas_por_lado.y,particulas_por_lado.z,particulas_por_lado.x*particulas_por_lado.y*particulas_por_lado.z);
    
    cel=InitV3<double3,double>(celda_minima.x,celda_minima.y,celda_minima.z);
    
    caja.x = cel.x*particulas_por_lado.x;
    caja.y = cel.y*particulas_por_lado.y;
    caja.z = cel.z*particulas_por_lado.z;
    double di=0.0;
    di=pow(densidad,1.0/nd);
    caja.x/=di;
    caja.y/=di;
    caja.z/=di;
    cel.x=caja.x/particulas_por_lado.x;
    cel.y=caja.y/particulas_por_lado.y;
    cel.z=caja.z/particulas_por_lado.z;
    printf("\nTamano de la caja de simulacion: (%.3lf,%.3lf,%.3lf)\n",caja.x,caja.y,caja.z);

    int x=0,y=0,z=0;

    moleculas_de_especie_acumuladas=new uint[n_esp_m];
    for(int i=0;i<n_esp_m;i++)moleculas_de_especie_acumuladas[i]=n_m_esp_mr[i];
    
    ofapin << "particula molecula especie_molecular especie_atomica"<< std::endl;
    for(int i=0;i<nm;i++){
        p_en_m[i]=part;
        if(moleculas_de_especie_acumuladas[k]<=0)k++;
        if(k>=n_esp_m)return;
        if(x>=particulas_por_lado.x){
            x=0;
            y++;
        }
        if(y>=particulas_por_lado.y){
            y=0;
            z++;
        }
        contp=0;
        for(int j=0;j<n_p_esp_m[k];j++){
            pos[part*nd]=x*cel.x+pos_respecto_p_central[max_p_en_esp_mr*k+j].x+centrar_m[nd*k];
            pos[part*nd+1]=y*cel.y+pos_respecto_p_central[max_p_en_esp_mr*k+j].y+centrar_m[nd*k+1];
            pos[part*nd+2]=z*cel.z+pos_respecto_p_central[max_p_en_esp_mr*k+j].z+centrar_m[nd*k+2];
            esp_p[part] = esp_p_en_esp_mr[k*max_p_en_esp_mr+j];
            m_de_p[part].x=i;
            m_de_p[part].y=contp;
            m_de_p[part].z=n_p_esp_m[k]-contp;
            ofapin << part << "\t" << i << "\t" << k << "\t" << esp_p_en_esp_mr[k*max_p_en_esp_mr+j]<<
            "\t" << pos[part*nd] << "\t" <<  pos[part*nd+1] << "\t" <<  pos[part*nd+2] << std::endl;
            part++;
            contp++;
        }
        x++;
        moleculas_de_especie_acumuladas[k]--;
    }
}

void InicializarVelocidades(double temp,double *v,int np)
{
    double sumx=0.,sumy=0.,sumz=0.;
    double ek=0.;
    AllPart
        v[ip*nd+id]=((double)rand()/RAND_MAX)-0.5;
    for(int i=0;i<np;i++)
    {
        sumx+=v[i*nd];
        sumy+=v[i*nd+1];
        sumz+=v[i*nd+2];
    }
    for(int i=0;i<np;i++){
        v[i*nd]-=sumx/np;
        v[i*nd+1]-=sumy/np;
        v[i*nd+2]-=sumz/np;
    }
    AllPart
        ek += v[ip*nd+id] * v[ip*nd+id];
    double tin=ek/(nd*(double)np*k_Boltzmann);
    double fac= sqrt(temp/tin);
    AllPart
        v[ip*nd+id]*=fac;
    


}

void DistanciasEntreParticulasEnMoleculaIniciales(uint np,ushort n_esp_m,uint max_p_en_esp_mr,uint *n_p_esp_mr,double *dis_p_esp_mr_rep,double3 *pos_respecto_p_central)
{
    /*
    Sirve para RATTLE y en caso de usar potencial de Hooke como restricción.
    */
    for(int i=0;i<n_esp_m;i++){
        for(int j=0;j<n_p_esp_mr[i];j++){
            for(int k=0;k<n_p_esp_mr[i];k++){
                dis_p_esp_mr_rep[i*max_p_en_esp_mr*max_p_en_esp_mr+j*max_p_en_esp_mr+k] = 0.;
                dis_p_esp_mr_rep[i*max_p_en_esp_mr*max_p_en_esp_mr+j*max_p_en_esp_mr+k] += pow(( pos_respecto_p_central[max_p_en_esp_mr*i+j].x - pos_respecto_p_central[max_p_en_esp_mr*i+k].x ),2);
                dis_p_esp_mr_rep[i*max_p_en_esp_mr*max_p_en_esp_mr+j*max_p_en_esp_mr+k] += pow(( pos_respecto_p_central[max_p_en_esp_mr*i+j].y - pos_respecto_p_central[max_p_en_esp_mr*i+k].y ),2);
                dis_p_esp_mr_rep[i*max_p_en_esp_mr*max_p_en_esp_mr+j*max_p_en_esp_mr+k] += pow(( pos_respecto_p_central[max_p_en_esp_mr*i+j].z - pos_respecto_p_central[max_p_en_esp_mr*i+k].z ),2);
                dis_p_esp_mr_rep[i*max_p_en_esp_mr*max_p_en_esp_mr+j*max_p_en_esp_mr+k] =sqrt(dis_p_esp_mr_rep[i*max_p_en_esp_mr*max_p_en_esp_mr+j*max_p_en_esp_mr+k]);
            }    
        }
    }
}


#endif