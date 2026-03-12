#ifndef POTENCIALES_INT_H
#define POTENCIALES_INT_H

#include "../MISC/CudaStruct.hpp"
#include "../MISC/Math.hpp"
#include "../MISC/Definiciones.hpp"

void ParamLJ(ushort n_esp_p,double *param,double *M_param)
{
    //como se pueden tener multiples especies de particulas
    //se calculan los promedios para cuando se calcula la interaccion

    //El potencial de Lennard-Jones tiene 2 parametros
    int M_6,M_12;
    double sig,d2,d6,eps;
    for(int i=0;i<n_esp_p;i++)
        for(int j=0;j<n_esp_p;j++){
            M_6=i*n_esp_p+j;
            M_12=M_6+n_esp_p*n_esp_p;
            sig=PromedioSuma2(param[j*n_esp_p],param[i*n_esp_p]);
            d2=sig*sig;
            d6=d2*d2*d2;
            eps=PromedioRaiz2(param[j*n_esp_p+1],param[i*n_esp_p+1]);
            M_param[M_6]=4.0*eps*d6;
            M_param[M_12]=M_param[M_6]*d6;
        }
}

double2 InteraccionLJ(int i,int j,ushort n_esp_p,double dis,double rc,double *M_param,uint *M_int,uint *esp_p)
{
    uint esp1=esp_p[i];
    uint esp2=esp_p[j];
    uint M_6=esp1*n_esp_p+esp2;
    uint M_12=M_6+n_esp_p*n_esp_p;
    
    double pot=0.,fue=0.;
    double d2=1.0/dis; //se modifica esto cuando sea en paralelo para evitar i=j
    double dc=rc?1./rc:0.;
    double d2c=dc*dc;
    double d6=d2*d2*d2;
    double d6c=d2c*d2c*d2c;
    double d12=d6*d6;
    double d12c=d6c*d6c;
    if((dis<=rc*rc)||(rc<=0.)){
        fue=6.0*(M_param[M_12]*2.0*d12-M_param[M_6]*d6)*d2;
        pot=M_param[M_12]*(d12-d12c)-M_param[M_6]*(d6-d6c);
    }else{
        fue=0.;
        pot=0.;
    }
    
    pot*=M_int[esp1*n_esp_p+esp2];
    fue*=M_int[esp1*n_esp_p+esp2];
    double2 fuepot;
    fuepot.x=fue;
    fuepot.y=pot;
    return fuepot;
}

double2 InteraccionYukawa(int i, int j,uint n_esp_p,uint elem_M,double dis, double *M_param,bool nconf,bool eshift,double r2c=0.0)
{
    //[1]
    double apot=0.;
    double2 val;
    bool gidj=i-j;

    double d2=gidj?(1./dis):0.;
    double ee = exp(-M_param[4*n_esp_p*n_esp_p+elem_M] * M_param[n_esp_p*n_esp_p+elem_M] * (dis / M_param[n_esp_p*n_esp_p+elem_M] - 1));

    double  fue= M_param[5*n_esp_p*n_esp_p+elem_M] * dis * pow((M_param[n_esp_p*n_esp_p+elem_M] / dis),M_param[5*n_esp_p*n_esp_p+elem_M]);
            fue+= M_param[2*n_esp_p*n_esp_p+elem_M] * M_param[3*n_esp_p*n_esp_p+elem_M] * (1 + M_param[4*n_esp_p*n_esp_p+elem_M] * dis) * (M_param[elem_M] * M_param[n_esp_p*n_esp_p+elem_M]) * ee;
            fue /= dis*dis;
    
    if(nconf){
        apot=pow((M_param[n_esp_p*n_esp_p+elem_M] *d2),M_param[5*n_esp_p*n_esp_p+elem_M]) + M_param[2*n_esp_p*n_esp_p+elem_M] * M_param[3*n_esp_p*n_esp_p+elem_M] * (M_param[elem_M] * M_param[n_esp_p*n_esp_p+elem_M] *d2) * ee;
    }
    
    val.x=fue;
    val.y=apot;
    return val;
}

void M_Parametros(ushort n_esp_p,double *param,double *M_param,potencial pot)
{
    switch (pot)
    {
    case LJ:
        ParamLJ(n_esp_p,param,M_param);
        break;
    case Yuk:
        break;
    }
}

double2 Interaccion(int i,int j,ushort n_esp_p,double dis,double rc, double *M_param,uint *M_int,uint *esp_p,potencial pot)
{
    double2 en;
    en.x=en.y=0.;
    switch (pot)
    {
    case LJ:
        return InteraccionLJ(i,j,n_esp_p,dis,rc,M_param,M_int,esp_p);
    case Yuk:
        return en;
    default:
        printf("error en la ejecucion\n");
        exit(EXIT_FAILURE);
        return en;    
    }
}

#endif