#ifndef DATOS_INI_H
#define DATOS_INI_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <math.h>

#include "../MISC/Definiciones.hpp"
#include "../MISC/CudaStruct.hpp"
#include "../MISC/DataTypes.hpp"

void DatosIniciales1(str dir, ushort &n_esp_m, ushort &n_esp_p)
{
    str f = dir + "DATOS/DatosSistema.txt";
    str key;
    double value;
    std::ifstream iff(f);
    PAbrioArchivo(f,iff);

    while (iff >> key >> value) {
        if (key == "n_esp_m") n_esp_m = (ushort)value;
        else if (key == "n_esp_p") n_esp_p = (ushort)value;
    }

    iff.close();
}

void DatosIniciales2(str dir,ushort n_esp_m,uint *n_m_esp_mr,uint *n_p_esp_m,uint &np,uint &nm)
{
    str f=dir+"DATOS/DatosSistema.txt";
    str temp;
    np=nm=0;
    std::ifstream iff(f);
    PAbrioArchivo(f,iff);
    str line;

    while (std::getline(iff,line)){
        std::istringstream iss(line);
        
        iss >> temp;

        if(temp == "n_m_esp_mr"){
            for(uint i=0;i<n_esp_m;i++){
                if(!(iss >>n_m_esp_mr[i])){
                    std::cerr << "Error al leer n_m_esp_mr" << std::endl;
                    break;
                    exit(EXIT_FAILURE);
                }
            }
        }else if(temp == "n_p_esp_m"){
            for(uint i=0;i<n_esp_m;i++){
                if(!(iss >>n_p_esp_m[i])){
                    std::cerr << "Error al leer n_p_esp_m" << std::endl;
                    break;
                    exit(EXIT_FAILURE);
                }
            }
        }
    }

    for(int i=0;i< n_esp_m;i++)
        nm+=n_m_esp_mr[i];
    for(int i=0;i< n_esp_m;i++)
        np+=n_m_esp_mr[i]*n_p_esp_m[i];

    iff.close();
}

void DatosCorrida(str dir,uint &nc,float &ncp,coordenadas &coord,ensamble &ens,termostatos &termo,potencial &pot,optimizaciones &opt,double &dt,double &tempe,double &rc,double &dens,double &p_termo,algoconstr &vibrante)
{
    str f=dir+"DATOS/DatosCorrida.txt";
    str key;
    double value;
    uint tempo;
    std::ifstream iff(f);
    PAbrioArchivo(f,iff);
    uint cvec=0,ccel=0;

    while (iff >> key >> value) {
        if (key == "nc") nc = static_cast<uint>(value);
        else if (key == "ncp") ncp = static_cast<float>(value);
        else if (key == "dt") dt = value;
        else if (key == "temp") tempe = value;
        else if (key == "dens") dens = value;
        else if (key == "potencial") pot = static_cast<potencial>(static_cast<uint>(value));
        else if (key == "ensamble") ens = static_cast<ensamble>(static_cast<uint>(value));
        else if (key == "termo") termo = static_cast<termostatos>(static_cast<uint>(value));
        else if (key == "p_termo") p_termo = static_cast<double>(value);
        else if (key == "vibrante")vibrante = static_cast<algoconstr>(static_cast<uint>(value)); 
        else if (key == "rc") rc = value;
        else if (key == "ovec") cvec = static_cast<int>(value);
        else if (key == "ocel") ccel = static_cast<int>(value);
        else if (key == "coord") coord = static_cast<coordenadas>(static_cast<uint>(value));
    }
    iff.close();
    tempo=cvec+2*ccel;
    opt = static_cast<optimizaciones>(tempo);
}

void DatosParticulas(str dir,double *param,double *masa,potencial pot,ushort n_esp_p)
{
    /*esta rutina asume que siempre el dato que esta en la primera columna de este archivo es el diametro de la especie atomica */
    str f=dir+"DATOS/DatosParticulas.txt";
    str tp;
    ushort nparam;

    nparam = DefNparam(pot);
    
    double temporal;
    std::ifstream iff(f);
    PAbrioArchivo(f,iff);
    iff >> tp;iff >> tp;iff >> tp;
    for(int j=0;j<n_esp_p;j++){
        iff >> temporal;
        masa[j]=temporal;
        for(int i=0;i<nparam;i++){
            iff >> temporal;
            param[nparam*j+i]=temporal;
        }
    }
    iff.close();
}

void DatosMoleculas(str dir,ushort n_esp_m,coordenadas coord,uint *n_p_esp_m,uint *esp_p_en_esp_mr,uint max_p_en_esp_mr,double3 *pos_respecto_p_central)
{
    const double pi=2*acos(0.);
    uint tui;
    str f=dir+"DATOS/DatosMoleculas.txt";
    str tp;
    double temporal,t1,t2,t3;
    std::ifstream iff(f);
    PAbrioArchivo(f,iff);
    iff >> tp;iff >> tp;iff >> tp;
    for(int i=0;i<n_esp_m;i++){
        iff >> tui;
        if(tui != i){ std::cerr << "Los Datos no fueron bien introducidos, todavia faltan los datos de las particulas de la especie" << i << std::endl;}
        for(int j=0;j<n_p_esp_m[i];j++){
            iff >> tui;
            esp_p_en_esp_mr[max_p_en_esp_mr*i+j]=tui;
            switch (coord)
            {
            case cartesianas:
                iff >> temporal;pos_respecto_p_central[max_p_en_esp_mr*i+j].x=temporal;
                iff >> temporal;pos_respecto_p_central[max_p_en_esp_mr*i+j].y=temporal;
                iff >> temporal;pos_respecto_p_central[max_p_en_esp_mr*i+j].z=temporal;
                break;
            
            case esfericos:
                iff >>t1;iff >>t2;iff >>t3;
                t2*=pi/180.;
                t3*=pi/180.;
                pos_respecto_p_central[max_p_en_esp_mr*i+j].x=t1*sin(t2)*cos(t3);
                pos_respecto_p_central[max_p_en_esp_mr*i+j].y=t1*sin(t2)*sin(t3);
                pos_respecto_p_central[max_p_en_esp_mr*i+j].z=t1*cos(t2);
                break;
            }
        }
    }
    iff.close();
}

void DatosInteraccion(str dir,ushort n_esp_p,uint *M_int)
{
    str f=dir+"DATOS/DatosInteraccion.txt";
    str tp;
    std::ifstream iff(f);
    uint tmp;
    PAbrioArchivo(f,iff);
    for(int i=0;i<n_esp_p;i++)
    {
        for(int j=0;j<n_esp_p;j++)
        {
            iff >> tmp;
            M_int[n_esp_p*i+j]=tmp;
            printf("Interaccion de especies %d-%d: ",i,j);
            if(tmp)printf("Activa\n"); else printf("Inactiva\n"); 
        }
    }
    iff.close();
}

void DatosConstriccion(str dir,algoconstr vibrante,double *constr)
{
    /*esta rutina asume que siempre el dato que esta en la primera columna de este archivo es el diametro de la especie atomica */
    str f=dir+"DATOS/DatosConstricciones.txt";
    str line;
    std::ifstream iff(f);
    PAbrioArchivo(f,iff);

    switch (vibrante)
    {
    case fuerzasres:
        while (std::getline(iff,line)){
            std::istringstream iss(line);
            str temp;
            if(line.find("kres") != str::npos){
                double val;
                iss >> temp >> val;
                constr[0] = val;
            }
        }
        break;
    
    case rattle:
        while (std::getline(iff,line)){
            std::istringstream iss(line);
            str temp;
            if(line.find("tol") != str::npos){
                double val;
                iss >> temp >> val;
                constr[0] = val;
            }else if(line.find("maxit") != str::npos){
                double val;
                iss >> temp >> val;
                constr[1] = val;
            }
        }
        break;
    }

    iff.close();
}

void AbrirArchivos(str directorio,double dens, ushort nem, ushort nea,potencial pot, uint *nme,ensamble ens,termostatos termo,uint nc,double param_termo, double *param,std::ofstream &ofaedi,std::ofstream &ofapin,std::ofstream &ofaeva,str &dpsco,algoconstr vibrante)
{
    /*****************************************/
    ushort nparam;
    std::stringstream stream;
    std::stringstream strea;
    std::stringstream *stream1;
    str dpsc1,dpsc,ss,s,aedi,apin,aeva;
    /*****************************************/
    nparam=DefNparam(pot);
    printf("nparam: %d\n",nparam);
    stream1=new std::stringstream[nparam*nea];
    dpsc1 = directorio + "CORRIDAS";
    
    create_directory(dpsc1.c_str());
    
    dpsc = dpsc1;

    switch (ens)
    {
    case nve:
        dpsc+="/NVE_";
        break;
    
    case nvt:
        strea << std::fixed << std::setprecision(2) << param_termo;
        ss = strea.str();
        dpsc+="/NVT_";
        switch (termo)
        {
        case rescvel:
            dpsc+="VRes_"+ss+"_";
            break;
        case andersen:
            dpsc+="And_"+ss+"_";
            break;
        case berendsen:
            dpsc+="Ber_"+ss+"_";
            break;
        case bdp:
            dpsc+="BDP_"+ss+"_";
            break;
        case nosehoover:
            dpsc+="NH_"+ss+"_";
            break;
        
        }break;
    }

    switch (vibrante)
    {
    case rattle:
        dpsc+="rattle";
        break;
    
    case fuerzasres:
        dpsc+="fuerres";
        break;
    }
    
    stream << std::fixed << std::setprecision(2) << dens;
    s = stream.str();
    dpsc += "_dens_"+ s;
    dpsc += "_nd_" + std::to_string(nd) + "_nem_" + std::to_string(nem) + "_nea_" + std::to_string(nea);

    for(int i=0;i<nem;i++){
        dpsc+="_nme"+std::to_string(i)+"_"+std::to_string(nme[i]);
    }
    for(int j=0;j<nea;j++){
        dpsc += "_diametro"+std::to_string(j)+"_";
        stream1[j*nparam] << std::fixed << std::setprecision(2) << param[j*nparam];
        s = stream1[j*nparam].str();
        dpsc += s;
    }
    for(int j=0;j<nea;j++)
        for(int i=1;i<nparam;i++){
            dpsc += "_param_"+std::to_string(i)+"_"+std::to_string(j)+"_";
            stream1[nparam*j+i] << std::fixed << std::setprecision(2) << param[nparam*j+i];
            s = stream1[nparam*j+i].str();
            dpsc += s;    
        }
    dpsc += "_nc_"+std::to_string(nc);

    std::cout << "LOS ARCHIVOS DE ESTA CORRIDA SE GUARDAN EN: " << dpsc << std::endl;
    
    create_directory(dpsc.c_str());
    dpsco = dpsc+"/Resultados";
    create_directory(dpsco.c_str());

    aedi = dpsc + "/DatosIniciales.txt";
    aeva = dpsc + "/VelAcel.txt";
    apin = dpsc + "/Posiciones_Iniciales.txt";

    ofaedi.open(aedi.c_str());
    ofaeva.open(aeva.c_str());
    ofapin.open(apin.c_str());
}

void ImpresionDeDatos(uint nc,float ncp,double dt,double temp,double rc,optimizaciones opt,
                      potencial pot,double dens,ushort n_esp_m,ushort n_esp_p,
                      ensamble ens,termostatos termos,uint *n_m_esp_mr,uint *n_p_esp_m, 
                      uint *esp_de_p_en_m,uint max_p_en_esp_mr,double3 *pos_respecto_p_central)
{
    int ncc=(float)nc*(ncp/100.);
    printf("\n---------------------------------------------------------------------------\n\n");
    printf("Datos de Corrida:\n\n");
    printf("Configuraciones: %d\nCada cuantas configuraciones imprimimos propiedades:%d\n",nc,ncc);
    printf("Tamaño del paso de integracion: %.3lf\n",dt);
    printf("Ensamble: ");
    switch (ens)
    {
    case nve:
        printf("NVE\n");
        break;
    
    case nvt:
        printf("NVT\n");
        printf("Temperatura a la que se debe mantener el sistema: %.1lf\n",temp);
        printf("Termostato: ");
        switch (termos)
        {
        case rescvel:
            printf("Reescalamiento de velocidades\n");
            break;
        case andersen:
            printf("Andersen\n");
            break;
        case berendsen:
            printf("Berendsen\n");
            break;
        case bdp:
            printf("Bussi-Donaldio-Parinello\n");
            break;
        case nosehoover:
            printf("Nose-Hoover\n");
            break;
        }
        break;
    }
    printf("Radio de corte (en caso de optimizaciones): %.1lf\n",rc);
    switch (opt)
    {
    case no_opt:
        printf("Sin Optimizaciones\n");
        break;
    case vec:
        printf("Optimizacion de vecinos\n");
        break;
    case cel:
        printf("Optimizacion de celdas\n");
        break;
    case ambas:
        printf("Optimizacion de vecinos y celdas\n");
        break;
    }
    printf("Potencial: ");
    switch(pot)
    {
        case LJ:printf("Lennard-Jones\n");break;
        case Yuk:printf("Yukawa\n");break;
    }
    printf("\nDatos del Sistema:\n\n");
    printf("Densidad: %.2lf\n",dens);
    printf("Dimensiones: %d\n",nd);
    printf("\nDatos de Atomos y Moleculas\n\n");
    printf("Especies moleculares: %d\n",n_esp_m);
    printf("Especies atomicas: %d\n\n",n_esp_p);
    for(int i=0;i<n_esp_m;i++)printf("Moleculas de la especie %d: %d\n",i,n_m_esp_mr[i]);
    for(int i=0;i<n_esp_m;i++)printf("Atomos en la especie molecular %d: %d\n",i,n_p_esp_m[i]);
    printf("\n");
    for(int i=0;i<n_esp_m;i++)for(int j=0;j<n_p_esp_m[i];j++)printf("Especie de los atomos en la especie molecular %d: %d\n",i,esp_de_p_en_m[i*max_p_en_esp_mr+j]);
    printf("Respecto al atomo central de la especie molecular correspondiente: \n\n");
    for(int i=0;i<n_esp_m;i++)
        for(int j=0;j<n_p_esp_m[i];j++)
            printf("Posiciones de los atomos en la especie molecular %d: (%.4lf,%.4lf,%.4lf)\n",i,pos_respecto_p_central[max_p_en_esp_mr*i+j].x,pos_respecto_p_central[max_p_en_esp_mr*i+j].y,pos_respecto_p_central[max_p_en_esp_mr*i+j].z);
    printf("\n---------------------------------------------------------------------------\n\n");
}

void ImpresionDeDatosADisco(uint nc,float ncp,double dt,double temp,double rc,optimizaciones opt,
                            potencial pot,double dens,ushort n_esp_m,ushort n_esp_p,ensamble ens,
                            termostatos termos,uint *n_m_esp_mr,uint *n_p_esp_m,uint *esp_de_p_en_m,
                            uint max_p_en_esp_mr,double3 *pos_respecto_p_central,std::ofstream &ofaedi)
{
    if(!ofaedi){std::cerr << "error al abrir ofaedi" << std::endl; exit(EXIT_FAILURE);}
    /**************/
    int ncc=(float)nc*(ncp/100.);
    /**************/
    ofaedi << "Datos de Corrida:\n"<< std::endl;
    ofaedi << "Configuraciones: " << nc << "\nCada cuantas configuraciones imprimimos propiedades:" << ncc << std::endl;
    ofaedi << "Tamaño del paso de integracion:" << std::setprecision(3) << dt << std::endl;
    ofaedi << "Ensamble: ";
    switch (ens)
    {
    case nve:
        ofaedi<< "NVE"<<std::endl;
        break;
    
    case nvt:
        ofaedi << "NVT"<<std::endl;
        ofaedi << "Temperatura a la que se debe mantener el sistema (CASO NVT):" << std::setprecision(2) << temp << std::endl;
        ofaedi << "Termostato: ";
        switch (termos)
        {
        case rescvel:
            ofaedi << "Reescalamiento de velocidades"<<std::endl;
            break;
        case andersen:
            ofaedi << "Andersen"<<std::endl;
            break;
        case berendsen:
            ofaedi << "Berendsen"<<std::endl;
            break;
        case bdp:
            ofaedi << "Bussi-Donaldio-Parinello"<<std::endl;
            break;
        case nosehoover:
            ofaedi << "Nose-Hoover"<<std::endl;
            break;
        }
        break;
    }

    ofaedi << "Radio de corte (en caso de optimizaciones):" << std::setprecision(2) << rc <<std::endl;
    switch (opt)
    {
    case no_opt:
        ofaedi << "Sin Optimizaciones"<<std::endl;
        break;
    case vec:
        ofaedi << "Optimizacion de vecinos"<<std::endl;
        break;
    case cel:
        ofaedi << "Optimizacion de celdas"<<std::endl;
        break;
    case ambas:
        ofaedi << "Optimizacion de vecinos y celdas"<<std::endl;
        break;
    }

    switch(pot)
    {
        case LJ:ofaedi << "Lennard-Jones"<<std::endl;break;
        case Yuk:ofaedi << "Yukawa"<<std::endl;break;
    }
    ofaedi << "\nDatos del Sistema:\n"<<std::endl;
    ofaedi << "Densidad:" << std::setprecision(3) << dens <<std::endl;
    ofaedi << "Dimensiones:" << nd << std::endl;
    ofaedi << "\nDatos de Atomos y Moleculas\n"<<std::endl;
    ofaedi << "Especies moleculares:" << n_esp_m << std::endl;
    ofaedi << "Especies atomicas: " << n_esp_p << "\n" <<std::endl;
    for(int i=0;i<n_esp_m;i++)ofaedi << "Moleculas de la especie " << i << ": " << n_m_esp_mr[i] << std::endl;
    for(int i=0;i<n_esp_m;i++)ofaedi << "Atomos en la especie molecular " << i << ": " << n_p_esp_m[i] << std::endl;
    ofaedi << std::endl;
    for(int i=0;i<n_esp_m;i++)for(int j=0;j<n_p_esp_m[i];j++)ofaedi << "Especie de los atomos en la especie molecular " << i << ": " << esp_de_p_en_m[i*max_p_en_esp_mr+j] << std::endl;
    ofaedi << "Respecto al atomo central de la especie molecular correspondiente: \n"<<std::endl;
    for(int i=0;i<n_esp_m;i++)for(int j=0;j<n_p_esp_m[i];j++)ofaedi << "Posiciones de los atomos en la especie molecular " << i << ": (" << std::setprecision(4) << pos_respecto_p_central[max_p_en_esp_mr*i+j].x << "," << std::setprecision(4) << pos_respecto_p_central[max_p_en_esp_mr*i+j].y << "," << std::setprecision(4) << pos_respecto_p_central[max_p_en_esp_mr*i+j].z << ")"<<std::endl;
    
}

void ArchivosDeResultados(str dpsco,std::ofstream &ofasres,std::ofstream &ofasat,optimizaciones op)
{
    str k;
    switch(op)
    {
        case no_opt: k="SinOptimizaciones";break;
        case vec: k="Vecinos";break;
        case cel: k="Celdas";break;
        case ambas: k="AmbasOptimizaciones";break;
    }
    std::cout << "Optimizacion: " << k << std::endl;
    str asat = dpsco + "/Posiciones_"+k+".txt";
    str asres = dpsco + "/Resultados_"+k+".txt";

    ofasat.open(asat.c_str());
    ofasres.open(asres.c_str());
}

uint CalculoConstricciones(ushort n_esp_m,uint *n_p_esp_m,uint *n_m_esp_mr)
{
    //como ahorita tenemos moleculas completamente rigidas, para cada especie molecular se tendran
    //N(N-1)/2 consrticciones por molecula, donde N es el numero de particulas en esa especie molecular
    uint constr=0;
    for(int i=0;i<n_esp_m;i++)
        constr+=(n_p_esp_m[i]*(n_p_esp_m[i]-1)*0.5)*n_m_esp_mr[i];
    return constr;
}
#endif