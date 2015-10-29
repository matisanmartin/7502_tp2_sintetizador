#ifndef VALIDAR_H
#define VALIDAR_H

#define CANT_MIN_ARG 7
#define CANT_MAX_ARG 9
#define PARTITURA "-p"
#define SINTETIZADOR "-s"
#define AUDIO "-o"
#define FRECUENCIA "-f"
#define FRECUENCIA_DEFECTO 44000

#define MSJ_ERROR_PARTITURA_INVALIDO "No se encontro el archivo de partitura\n"
#define MSJ_ERROR_SINTETIZADOR_INVALIDO "No se encontro el archivo sintetizador\n"
#define MSJ_ERROR_ARGUMENTO_INVALIDO "Error en el argumento\n"
#define MSJ_ERROR_FRECUENCIA_INCORRECTA "Frecuencia incorrecta\n"

status_t validar(int, char**, unsigned int *frecuencia, FILE ** partitura, FILE ** sintetizador, FILE **audio); 

#endif
