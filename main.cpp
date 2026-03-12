#include <stdio.h>
#include <iostream>
#include <fstream>

#include "MISC/Definiciones.hpp"
#include "MISC/CudaStruct.hpp"
#include "MISC/DataTypes.hpp"
#include "SISINI/DatosIniciales.hpp"
#include "SISINI/ConfIni.hpp"
#include "DM/SimNoOpt.hpp"
#include "DM/SimVec.hpp"
#include "DM/SimCel.hpp"
#include "DM/SimCV.hpp"

int main()
{
    //************Datos Programa*************
    str dir = "/home/ornit/DMM/";       //directorio en el que se encuentra main.cpp.
    str dpsco;                          //directorio donde esta la carpeta de corridas.
    std::ofstream ofaedi;               //archivo de datos iniciales.
    std::ofstream ofapin;               //archivo de posiciones iniciales.
    std::ofstream ofasres;              //archivo de propiedades termodinámicas.
    std::ofstream ofasat;               //archivo de posiciones, velocidades y aceleraciones.
    std::ofstream ofaeva;               //archivo de velocidades y aceleraciones iniciales.
    //************Datos Corrida**************
    uint nc;                            //Número de configuraciones que se exploran en el sistema.
    float ncp;                          //Cada que porcentaje se imprimen propiedades termodinámicas.
    double rc;                          //Radio de corte para optimizaciones.
    double dt;                          //Tamano del paso de integración.
    optimizaciones opt;                 //Optimización que se utiliza para la simulación.
    coordenadas coord;                  //Tipo de coordenadas en el archivo DatosMolecula.txt
    termostatos termo;                  //Tipo de termostato.
    algoconstr vibrante;
    ensamble ens;
    //***********Datos Sistema***************
    double tempe;
    double dens;
    double p_termo;
    bool p_o_m;
    ushort n_esp_p;
    ushort n_esp_m;
    uint GdL;
    uint *M_int;
    double3 caja;
    double3 celda_min;
    //***********Datos Moleculas*************
    uint nm;
    uint n_constr;
    uint *n_m_esp_mr;
    uint *n_p_esp_m;
    uint max_p_en_esp_mr;
    uint *esp_p_en_esp_mr;
    uint *p_en_m;
    double *constr;
    //***********Datos Particulas************
    uint np;
    ushort nparam;
    uint *esp_p;
    uint3 *mad_de_p;
    double *pos,*vel,*acel,*q_rat;
    double *param;
    double3 *pos_respecto_p_central;
    potencial pot;
    double *masa;
    //***************************************
    uchar n_par_const=0,arr_temp=0;
    srand(time(NULL));
    
    std::cout << "dir: " << dir << std::endl;
    
    DatosIniciales1(dir,n_esp_m,n_esp_p);

    n_m_esp_mr = new uint[n_esp_m];
    n_p_esp_m = new uint[n_esp_m];
    
    DatosIniciales2(dir,n_esp_m,n_m_esp_mr,n_p_esp_m,np,nm);
    
    pos   = new double[np*nd];
    q_rat = new double[np*nd];
    vel   = new double[np*nd];
    acel  = new double[np*nd];
    p_en_m= new uint[nm];

    DatosCorrida(dir,nc,ncp,coord,ens,termo,pot,opt,dt,tempe,rc,dens,p_termo,vibrante);
    
    nparam=DefNparam(pot);

    switch (vibrante)
    {
    case fuerzasres:
        n_par_const=1;
        break;
    case rattle:
        n_par_const=2;
        break;
    }

    constr = new double[static_cast<int>(n_par_const)];
    param = new double[nparam*n_esp_p];
    M_int = new uint[n_esp_p*n_esp_p];
    masa = new double[n_esp_p];
    for(int i=0;i<n_esp_m;i++)
    if(n_p_esp_m[i]>=arr_temp)arr_temp=n_p_esp_m[i];
    max_p_en_esp_mr = arr_temp;arr_temp=0;
    esp_p_en_esp_mr = new uint[max_p_en_esp_mr*n_esp_m];
    pos_respecto_p_central = new double3[max_p_en_esp_mr*n_esp_m];
    
    DatosParticulas(dir,param,masa,pot,n_esp_p);
    DatosMoleculas(dir,n_esp_m,coord,n_p_esp_m,esp_p_en_esp_mr,max_p_en_esp_mr,pos_respecto_p_central);
    DatosInteraccion(dir,n_esp_p,M_int);
    DatosConstriccion(dir,vibrante,constr);
    AbrirArchivos(dir,dens,n_esp_m,n_esp_p,pot,n_m_esp_mr,ens,termo,nc,p_termo,param,ofaedi,ofapin,ofaeva,dpsco,vibrante);
    ImpresionDeDatos(nc,ncp,dt,tempe,rc,opt,pot,dens,n_esp_m,n_esp_p,ens,termo,n_m_esp_mr,n_p_esp_m,esp_p_en_esp_mr,max_p_en_esp_mr,pos_respecto_p_central);
    ImpresionDeDatosADisco(nc,ncp,dt,tempe,rc,opt,pot,dens,n_esp_m,n_esp_p,ens,termo,n_m_esp_mr,n_p_esp_m,esp_p_en_esp_mr,max_p_en_esp_mr,pos_respecto_p_central,ofaedi);
    
    celda_min=CreandoCeldaMinima(n_esp_m,pos_respecto_p_central,max_p_en_esp_mr,param,nparam,esp_p_en_esp_mr,n_p_esp_m);
    double *centrar_m=new double[nd*n_esp_m];

    CentrarMoleculas(centrar_m,n_esp_m,n_p_esp_m,esp_p_en_esp_mr,max_p_en_esp_mr,nparam,param,pos_respecto_p_central);

    esp_p = new uint[np];
    mad_de_p = new uint3[np];
    ConfiguracionCubica(n_esp_m,n_m_esp_mr,n_p_esp_m,p_en_m,pos,pos_respecto_p_central,max_p_en_esp_mr,caja,centrar_m,celda_min,dens,ofapin,esp_p_en_esp_mr,nm,esp_p,mad_de_p);
    
    double3 caja_i=InitV3<double3,double>(1.0/caja.x,1.0/caja.y,1.0/caja.z);
    double *dis_p_esp_mr_rep = new double[n_esp_m*max_p_en_esp_mr*max_p_en_esp_mr];

    InicializarVelocidades(tempe,vel,np);

    if(max_p_en_esp_mr-1)DistanciasEntreParticulasEnMoleculaIniciales(np,n_esp_m,max_p_en_esp_mr,n_p_esp_m,dis_p_esp_mr_rep,pos_respecto_p_central);
    
    double rbuf=0.5;
    ArchivosDeResultados(dpsco,ofasres,ofasat,opt);

    n_constr=CalculoConstricciones(n_esp_m,n_p_esp_m,n_m_esp_mr);
    printf("numero de constricciones: %d\n",n_constr);

    GdL=nd*np-n_constr;
    printf("Grados de libertad: %d\n",GdL);

    switch (opt)
    {
    case no_opt:
        SimulacionNOpt(ncp,n_esp_p,n_esp_m,nc,np,max_p_en_esp_mr,GdL,M_int,esp_p,n_p_esp_m,n_m_esp_mr,
                       p_en_m,mad_de_p,dens,dt,tempe,p_termo,masa,q_rat,pos,vel,acel,constr,param,dis_p_esp_mr_rep,
                       caja,ens,vibrante,termo,pot,ofasres,ofasat,ofaeva);
        break;
    case vec:
        SimulacionVec(ncp,n_esp_p,n_esp_m,nc,np,max_p_en_esp_mr,GdL,M_int,esp_p,n_p_esp_m,n_m_esp_mr,
                      p_en_m,mad_de_p,dens,dt,tempe,p_termo,rc,masa,q_rat,pos,vel,acel,constr,param,
                      dis_p_esp_mr_rep,caja,ens,vibrante,termo,pot,ofasres,ofasat,ofaeva);
        break;
    case cel:
        SimulacionCel(ncp,n_esp_p,n_esp_m,nc,np,max_p_en_esp_mr,GdL,M_int,esp_p,n_p_esp_m,n_m_esp_mr,
                      p_en_m,mad_de_p,rc,dens,dt,tempe,p_termo,masa,q_rat,pos,vel,acel,constr,param,
                      dis_p_esp_mr_rep,caja,ens,vibrante,termo,pot,ofasres,ofasat,ofaeva);
        break;
    case ambas:
        SimulacionCV(ncp,n_esp_p,n_esp_m,nc,np,GdL,max_p_en_esp_mr,M_int,esp_p,n_p_esp_m,n_m_esp_mr,
                      p_en_m,mad_de_p,dens,dt,tempe,p_termo,rc,masa,q_rat,pos,vel,acel,constr,param,
                      dis_p_esp_mr_rep,caja,ens,vibrante,termo,pot,ofasres,ofasat,ofaeva);
        break;
    }
    delete[] pos;
    delete[] vel;
    delete[] acel;
    
    ofaedi.close();
    ofapin.close();
    ofasres.close();
    ofasat.close();

    return 0;
}