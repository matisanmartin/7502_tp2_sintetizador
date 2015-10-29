#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "common.h"
#include "validar.h"

status_t validar (int argc, char * argv[], unsigned int *frecuencia, FILE ** partitura, FILE ** sintetizador, FILE **audio)
{
	char *perr;
	size_t i;
	
	*frecuencia = FRECUENCIA_DEFECTO;
	
	for(i = 1; i < argc; i += 2)
	{
		if(!strcmp(argv[i], FRECUENCIA))
		{
			*frecuencia = strtol(argv[i + 1], &perr, 10);
			if(*perr != '\0')
				return ERROR_FRECUENCIA_INVALIDA;
			else break;
		}
	}
	
	if(i < argc)/*hizo el break asi que encontro -f*/
	{
		if (argc != CANT_MAX_ARG)
			return ERROR_CANT_ARG;
	}
	else 
	{
		if (argc != CANT_MIN_ARG)
			return ERROR_CANT_ARG;
	}
	
	for(i = 1; i < argc; i += 2)
	{
		if (!strcmp(argv[i], PARTITURA))
		{ 
			if ((*partitura = fopen(argv[i + 1], "rt")) == NULL)
				return ERROR_PARTITURA_INVALIDA;
			break;
		}	
	}
	if (i >= argc)
		return ERROR_ARGUMENTO_INVALIDO;
		
	for(i = 1; i < argc; i += 2)
	{
		if (!strcmp(argv[i], SINTETIZADOR))
		{
			if ((*sintetizador = fopen(argv[i + 1], "rt")) == NULL)
				return ERROR_SINTETIZADOR_INVALIDO;
			break;
		}
	}
	if (i >= argc)
		return ERROR_ARGUMENTO_INVALIDO;
		
	for(i = 1; i < argc; i += 2)
	{
		if (!strcmp(argv[i], AUDIO)) 
		{
			if ((*audio=fopen(argv[i + 1], "wb")) == NULL)
				return ERROR_ARGUMENTO_INVALIDO;
			break;
		}
	}
	if (i >= argc)
		return ERROR_ARGUMENTO_INVALIDO;
	
	return OK;
}
