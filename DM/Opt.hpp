#ifndef OPTIMIZACIONES_H
#define OPTIMIZACIONES_H

#include "../MISC/CudaStruct.hpp"
#include "../MISC/Definiciones.hpp"
#include "FuncComp.hpp"
/*
*   @brief
*   La densidad maxima de esferas es:
*   https://en.wikipedia.org/wiki/Close-packing_of_equal_spheres
*   la expresion para calcular el numero de esferas completas de 
*   radio sigma que caben en una esfera de radio rcc es:
*    
*   0.75*((rcc*sig)/(0.5*sig))³=6*rcc³
    
*   Para calcular vecinos se calcula la distancia entre centros de las particulas, entonces,
*   se pueden considerar los centros de particulas que no caben completamente en la esfera,
*   por eso consideraremos un radio de contenedor de rcc+(radio de la particula).
*   Ademas se debe multiplicar por la densidad para considerar la separacion inicial entre
*   particulas.
*   @param 
*   n_esp_p     numero de especies de particulas
*   nparam      numero de parametros
*   np          numero de particulas
*   param       lista de parametros por especie de particula (en particular el radio)
*   rc          radio de corte
*   dens        densidad del sistema
*/
uint MaxNVec(ushort n_esp_p,uint nparam,uint np,double *param,double rc,double dens)
{
    double sig=1000.;
    for(int i=0;i<n_esp_p;i++)if(param[i*nparam]<=sig)sig=param[i*nparam];
    double rcc=rc+0.5*sig;
    double calcvec = dens*6.*rcc*rcc*rcc;
    uint nmaxvec = calcvec;
    if(nmaxvec>=np)nmaxvec=np;
    printf("Numero de particulas %d\nproposicion del numero maximo de vecinos: %d\n",np,nmaxvec);ILinea
    return nmaxvec;
}

/*
*   @brief 
*   Se considera el lado de longitud mayor de las celdas (lcm) para asumirlas cúbicas.
*   https://doi.org/10.37236/1786 contiene un algoritmo para empacar esferas de un 
*   mismo tamaño dentro de una caja unitaria. En la tabla 2 d_n se indica la distancia
*   máxima entre 2 esferas, es decir su diametro (sigma).
*   El calculo del diametro de partículas que caben en celdas de lado (lcm) es:
*   new_sig = o_sig*lcm, los valores del articulo estan en el arreglo a.
*   @param 
*   n_esp_p         numero de especies de particulas
*   nparam          numero de parametros por especie en param
*   param           lista de parametros por especie de particula
*   tamcel          dimensiones de las celdas
*/

uint NPCel(uint n_esp_p,uint nparam,double *param, double3 tamcel)
{
    double sig=1000.;
    for(int i=0;i<n_esp_p;i++)if(param[i*nparam]<=sig)sig=param[i*nparam];
    uint num_part=1;
    
    //obteniendo las dimensiones para escalar las celdas
    double sigmax=0.0;
    sigmax=(tamcel.x>=sigmax)?tamcel.x:sigmax;
    sigmax=(tamcel.y>=sigmax)?tamcel.y:sigmax;
    sigmax=(tamcel.z>=sigmax)?tamcel.z:sigmax;
    
    double a[]=
    {
        2.0,sqrt(3.0),sqrt(2.0),sqrt(2.0),sqrt(5.0)/2.0,3*sqrt(2.0)/4.0,(-4.0+sqrt(10.0+4.0*sqrt(3.0)))/sqrt(3.0),1.0,
        sqrt(3.0)/2.0,3.0/4.0,0.710116382462,0.707106806467,sqrt(2.0)/2.0,sqrt(2.0)/2.0,5.0/8.0,0.606667120726,
        3.0*sqrt(2.0)/7.0,sqrt(13.0)/6.0,0.578209612716,0.554761174904,3.0/(2.0+2.0*sqrt(3.0)),3*sqrt(2.0)/8.0,0.523539214257,
        0.517638090205,0.505135865094,0.501074021252,1.0/2.0,0.471410634842,sqrt(2.0)/3.0,sqrt(2.0)/3.0,sqrt(2.0)/3.0,
    };
    //Reescalando los valores del articulo
    int tam_a=sizeof(a)/sizeof(a[0]);
    for(int i=0;i<tam_a;i++)a[i]*=sigmax;

    //Casos limite
    if(sig>a[0]){
        printf("En las celdas cabe 1 particula\n");ILinea
        return 1;
    }
    if(sig<a[tam_a-1]){
        printf("Las particulas son muy pequeñas\n");ILinea
        exit(EXIT_FAILURE);
    }

    for(int i=0;i<tam_a;i++)if(sig<=a[i])num_part++;
    
    printf("En la celda caben %d particulas\n",num_part);ILinea
    return num_part;
}

void VecinosCercanos(uint np,uint nmaxvec,uint *nvec,uint *vec,uint3 *m_d_p,double3 caja,double *pos,double rc)
{
    uint mol1,mol2;
    double3 dij,di,dj;
    double dis;
    double3 cajai=InitV3<double3,double>(1.0/caja.x,1.0/caja.y,1.0/caja.z);
    rc+=0.5;

    for(uint i=0;i<np;i++)
    {
        nvec[i]=0;
        di=InitV3<double3,double>(pos[i*nd],pos[1+i*nd],pos[2+i*nd]);
        mol1=m_d_p[i].x;
        for(int j=i+1;j<np;j++){
            mol2=m_d_p[j].x;
            if(mol1!=mol2){
                dj=InitV3<double3,double>(pos[j*nd],pos[1+j*nd],pos[2+j*nd]);
                dij=Condper(di,dj,caja,cajai);
                dis=Distancia2(dij);
                if(dis<=rc*rc){
                    vec[i*nmaxvec+nvec[i]]=j;
                    vec[j*nmaxvec+nvec[j]]=i;
                    nvec[i]++;
                    nvec[j]++;
                }
            }
        }
    }
}


uint CrearCeldas(double3 caja,double3 &tamcel,double3 &invtamcel)
{
    uint nceldas=0;
    uint3 celdas;
    celdas.x=caja.x; printf("Celdas en x: %d\n",celdas.x);
    celdas.y=caja.y; printf("Celdas en y: %d\n",celdas.y);
    celdas.z=caja.z; printf("Celdas en z: %d\n",celdas.z);
    tamcel.x=caja.x/celdas.x;
    tamcel.y=caja.y/celdas.y;
    tamcel.z=caja.z/celdas.z;
    printf("Tamaño de celda minima: (%.4lf,%.4lf,%.4lf)\n",tamcel.x,tamcel.y,tamcel.z);
    printf("-----------------------------------------------\n");
    invtamcel=InitV3<double3,double>(1./tamcel.x,1./tamcel.y,1./tamcel.z);
    nceldas=celdas.x*celdas.y*celdas.z;
    return nceldas;
}

//TODO mejorar este algoritmo para la lista vaya desde la celda mas interna a la externa (tratando de redicir error computacional).
void CalculoDeCeldasVecinas(double3 caja,uint *veccel,short dc,uint nveccel)
{
    int ncx=0,ncy=0,ncz=0;
    int celv;
    int cont=0;
    uint eLx=caja.x,eLy=caja.y,eLz=caja.z;
    uint fa=eLx*eLy;
    int cel;

    for(int celz=0;celz<eLz;celz++){
        for(int cely=0;cely<eLy;cely++){
            for(int celx=0;celx<eLx;celx++){
                /********************************************/
                cont=0;
                cel=celx + eLx*cely + fa*celz;
                for(int cz=-dc;cz<=dc;cz++){
                    ncz=celz+cz;
                    if(ncz<0) ncz+=eLz;
                    if(ncz>=eLz) ncz-=eLz;
                    for(int cy=-dc;cy<=dc;cy++){
                        ncy=cely+cy;
                        if(ncy<0) ncy+=eLy;
                        if(ncy>=eLy) ncy-=eLy;
                        for(int cx=-dc;cx<=dc;cx++){
                            ncx=celx+cx;
                            if(ncx<0) ncx+=eLx;
                            if(ncx>=eLx) ncx-=eLx;
                            if(ncx>=0&&ncy>=0&&ncz>=0){
                                celv=ncx + eLx*ncy + fa*ncz;
                                veccel[nveccel*cel+cont]=celv;
                                cont++;
                            }
                        }
                    }
                }
                /********************************************/
            }
        }
    }
}

void AsignarCeldas(int np,uint max_p_cel,uint nceldas,uint *np_c,uint *p_c,double3 caja,double3 itcel,double *pos)
{
    uint3 cel;
    for(int i=0;i<nceldas;i++)np_c[i]=0;

    for(int i=0;i<np;i++){
        int fac1=caja.x;
        int fac2=caja.y; fac2*=fac1;
        double px=pos[nd*i];
        double py=pos[nd*i+1];
        double pz=pos[nd*i+2];
        if(px==caja.x&&py==caja.y&&pz==caja.z){px=py=pz=0.;}
        cel.x=px*itcel.x;
        cel.y=py*itcel.y;
        cel.z=pz*itcel.z;
        int ncel=cel.x+fac1*cel.y+fac2*cel.z;
        //if(ncel<0||ncel>nceldas)printf("error celda: %d,part: %d posiciones(%lf,%lf,%lf)\n",ncel,i,px,py,pz);
        p_c[max_p_cel*ncel+np_c[ncel]]=i;
        np_c[ncel]++;
        //if(np_c[ncel]>max_p_cel)printf("exceso de particulas %d en celda %d, particulas %d\n",i,ncel,np_c[ncel]);    
    }
}

void VecinosConCel(uint np,uint ncv,uint max_p_cel,uint nmaxvec,uint *nvec,uint *vec,uint *vcl,uint *np_c,uint *p_c,uint3 *m_d_p,double rc,double *pos,double3 itcel,double3 caja)
{
    uint3 cel;
    uint mol1,mol2;
    double3 dij,di,dj;
    double dis;
    double3 cajai=InitV3<double3,double>(1.0/caja.x,1.0/caja.y,1.0/caja.z);
    rc+=0.5;
    int fac1=caja.x;
    int fac2=caja.y; fac2*=fac1;

    for(int i=0;i<np;i++)
    {
        nvec[i]=0;
        di=InitV3<double3,double>(pos[i*nd],pos[1+i*nd],pos[2+i*nd]);
        mol1=m_d_p[i].x;
        cel.x=di.x*itcel.x;
        cel.y=di.y*itcel.y;
        cel.z=di.z*itcel.z;
        if(di.x==caja.x&&di.y==caja.y&&di.z==caja.z){cel.x=cel.y=cel.z=0.;}
        int ce=cel.x+fac1*cel.y+fac2*cel.z;
        for(int cv=0;cv<ncv;cv++){
            int clv=vcl[ce*ncv+cv];
            for(int qp=0;qp<np_c[clv];qp++){
                int j=p_c[max_p_cel*clv+qp];
                mol2=m_d_p[j].x;
                if((i<j)&&(i!=j)&&(mol1!=mol2)){
                    dj=InitV3<double3,double>(pos[j*nd],pos[1+j*nd],pos[2+j*nd]);
                    dij=Condper(di,dj,caja,cajai);
                    dis=Distancia2(dij);
                    if(dis<=rc*rc){
                        vec[i*nmaxvec+nvec[i]]=j;
                        vec[j*nmaxvec+nvec[j]]=i;
                        nvec[i]++;
                        nvec[j]++;
                    }
                }
            }
        }
        
    }
}

#endif