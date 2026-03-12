#ifndef DEFINICIONES_H
#define DEFINICIONES_H

#include <string>

using str = std::string;

#define MAX_INT 0xFFFFFFFF 
#define nd 3
#define k_Boltzmann 1

#define AllPart for(int ip=0; ip<np; ip++)   \
                    for(int id=0; id<nd; id++)

#define ParaMoleculas if(max_p_en_esp_mr-1)
        
#define PAbrioArchivo(arch,farch)\
        if(!farch)\
        {std::cout<<"error al abrir archivo "<<arch<<std::endl;exit(EXIT_FAILURE);}

#define ILinea printf("-----------------------------------------------\n");

#ifdef _WIN32
    #include <direct.h>
    #define MKDIR(path) _mkdir(path)
#elif __linux__
    #include <sys/stat.h>
    #include <sys/types.h>
    #define MKDIR(path) mkdir(path,0777)
#else
    std::cout << "Sistema operativo no detectado\n";
    exit(0);
#endif

void create_directory(const str &path){
    MKDIR(path.c_str());
}

#endif