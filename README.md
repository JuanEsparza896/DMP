# DMP
Programa de dinámica molecular con algoritmo de RATTLE y termostatos para simulación en NVT

El programa está planeado para ser usado en distribuciones de Linux.

Dentro del archivo main.cpp el string dir tiene que contener la ruta la carpeta donde se localiza el archivo main.cpp para que el programa realizar la simulación.

Para simular distintos sistemas se deben modificar los archivos dentro de la carpeta de DATOS, a continuación se muestra como se modifica cada archivo:

**DatosConstricciones.txt**

Las variables tol y maxit corresponden a la tolerancia y al máximo número de iteraciones para resolver las constricciones del algoritmo de RATTLE, el valor de kres corresponde al valor de la constante de restitución del potencial en caso de que se usa un resorte como fuerza de restricción para mantener partículas en una molécula unidas.

**DatosCorrida.txt**

Se muestra una tabla con la variable y lo que significa

| Variable | Significado|
| -------- | --- |
| nc |Número de configuraciones que se exploran en la simulación |
|ncp|Cada que porcentaje de la corrida se imprimen propiedades termodinámicas|
|dt|Tamaño del paso de integración|
|temp|Temperatura a la que se inicializa el sistema, o deseada en caso de NVT|
|dens|Densidad de partículas en el sistema|
|potencial|Potencial de interacción entre partículas|
|ensamble|ensamble que se simula|
|termo|termostato que se utiliza para simulación NVT|
|p_termo|parámetro asociado al termostato correspondiente|
|vibrante|si se utilizan potenciales de restricción o RATTLE para moléculas|
|rc|radio de corte para optimizaciones de celdas y vecinos cercanos|
|ovec|apagar o encender la optimización de vecinos cercanos|
|ocel|apagar o encender la optimización de celdas|
|coord|si las coordenadas en el archivo de moleculas son cartesianas o polares|

Para poder ver a que número le corresponde cada termostato u opciones de compilación se pueden ver como están organizados los enum en la carpeta MISC en el archivo DataTypes.hpp. los siguientes enum corresponden a los siguientes parámetros del archivo de datos

|Nombre del enum|Nombre en el archivo de datos|
|-----|---------|
|termostatos|termo|
|algoconstr|vibrante|
|ensamble|ensamble|
|optimizaciones|ocel y ovec|
|coordenadas|coord|
|potencial|potencial|

**DatosInteraccion.txt**

El archivo contiene una matriz que indica si distintas especies de partículas interactuan entre si

**DatosMoleculas.txt**

En este archivo se colocan las especies y las posiciones de las partículas en una molécula, se requiere de un orden específico para ingresar la información.

Primero se colocan todas las partículas que corresponden a la primera especie de molécula, las coordenadas son relativas a su centro de masa, de ahí se colocan la de las demás partículas correspondientes a la siguiente especie y de ahí en adelante.

Por ejemplo para un sistema de 2 especies de moleculas, la especie 0 con 3 moleculas y la segunda con 2 el archivo se ve de la siguiente manera

````
  especie_molecular   especies_particulas    posiciones_respecto_a_particula_central
  0  0  0.0 0.0 0.0
  0  1  0.0 0.5 0.0
  0  1  0.0 -0.5 0.0
  1  0  0.0 0.0 0.0
  1  0  0.0 0.0 0.4
````
Las especies de las partículas y las posiciones en este ejemplo no son relevantes, el objetivo es ilustrar como se vería este archivo.

**DatosParticulas.txt**

Cada fila corresponde a una esoecie de partícula, se coloca la masa, el diámetro y los parámetros que correspondan al potencial de interacción asociado.

**Datos_sistema.txt**

En este archivo se indican el número de especies de partículas, especies de moléculas, numero de moléculas de cada especie, numero de partículas en cada especie.

Los 2 ultimos parámetros son listas, un archivo con 2 especies de moléculas se vería de la siguiente manera

````
n_esp_m     1
n_esp_p     1
n_m_esp_mr  500  100
n_p_esp_m   2  1
````

Para compilar el archivo se ejecuta:

````
g++ main.cpp
````

Para ejecutarlo

````
./a.out
````
