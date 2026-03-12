#ifndef RAND_HEADER
#define RAND_HEADER

//generador de numeros no lineales para evitar usar rand()
#include <stdint.h>
#include <time.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

//Inicializando los estados, este es para cuando se hace debug (semilla fija)
static uint64_t semilladebug[4] = {0x123,0x456,0x789,0xABC};
//para poder expandir una semilla en 4 distintas se utiliza el algoritmo de splitmix64
uint64_t splitmix64(uint64_t *state)
{
    //para dividir la semilla se utilizan numeros irracionales
    //Este primer numero hexadecimal corresponde al cociente de 2^64 y la proporcion aurea
    uint64_t z = (*state += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z>>30))*0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z>>27))*0x94D049BB133111EBULL;
    return z ^ (z>>31);
}

void InicializarSemillas(uint64_t s[4])
{
    uint64_t estado = (uint64_t)time(NULL);
    s[0] =splitmix64(&estado);
    s[1] =splitmix64(&estado);
    s[2] =splitmix64(&estado);
    s[3] =splitmix64(&estado);
}

static inline uint64_t rotl(const uint64_t x, int k)
{
    return (x << k) | (x >>(64-k));
}

//generacion de numero aleatorio utilizando el algoritmo de Xoshiro
//da un numero aleatorio y modifica la semilla para obtener el siguiente
uint64_t Rand_Xoshiro(uint64_t s[4])
{
    const uint64_t result = rotl(s[1]*5,7)*9;
    const uint64_t t = s[1] << 17;
    
    s[2] ^= s[0];
    s[3] ^= s[1];
    s[1] ^= s[2];
    s[2] ^= s[3];

    s[2] ^= t;
    s[3] = rotl(s[3],45);

    return result;
}

/**
*   @brief Genera numeros aleatorios tipo double.
*
*   @details Como Rand_Xoshiro da un numero de 64 bits, se debe transformar en
*   uno de 52 para ser double, por eso se recorre 11 bits, 9007199254740992.0
*   representa 2^53 por lo que el numero final es uno en el intervalo [0.0,1.0)
*   @param uint64_t s[4] semillas originales.
*   @return numero aleatorio double
*/
double rand_d(uint64_t s[4])
{
    return (Rand_Xoshiro(s) >> 11)*(1.0/9007199254740992.0);
}

/**
 * @brief Genera un número aleatorio Gaussiano (Normal) con media 0 y desviación estándar sigma.
 * @param s[4] Semillas para generar numeros aleatorios.
 * @param sigma Desviación estándar deseada.
 * @return double Número aleatorio con distribución Gaussiana.
 */
double Num_Gaussiano(uint64_t  s[4],double sigma)
{
    double u1; 
    do{
        u1= rand_d(s);
    }while(u1<=0.);
    double u2 = rand_d(s);
    if(u1<=0. )u1=1e-15;

    return sigma*sqrt(-2.0 * log(u1))*cos(2.0*u2*M_PI);
}

//suma de numeros aleatorios de una distribucion gaussiana
double S2_Gaussiano(u_int GdL,uint64_t s[4],double sigma)
{
    double sn_g2=0.;
    for(int i=0;i<GdL-1;i++){
        double num=Num_Gaussiano(s,sigma);
        sn_g2+=num*num;
    }
    return sn_g2;
}

//genera distribucion gamma(a,1)
//usa metodo de marsaglia y tsang
double distribucion_gamma(uint64_t s[4],double a)
{
    double d,c,x,v,u;
    d=a-1./3.;
    c=1./sqrt(9.*d);
    for(;;){
        do{
            x=Num_Gaussiano(s,1.);
            v=1.+c*x;
        }while (v<=0.);

        v=v*v*v;
        u=rand_d(s);

        if(u<1.-0.0331*(x*x)*(x*x))
            return d*v;
        if(log(u)<.5*x*x+d*(1.-v+log(v)))
            return d*v;        
    }
}

double S2_Gaussiano2(uint64_t s[4],uint GdL)
{
    return 2.*distribucion_gamma(s,(GdL-1)/2.);
}

#endif