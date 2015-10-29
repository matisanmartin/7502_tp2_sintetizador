#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include "common.h"
#include "sintetizador.h"																

static const char *arr_funciones[CANT_FUNCIONES]={	"CONSTANT",
													"LINEAR",
													"INVLINEAR",
													"SIN",
													"EXP",
													"INVEXP",
													"QUARTCOS",
													"QUARTSIN",
													"HALFCOS",
													"HALFSIN",
													"LOG",
													"INVLOG",
													"TRI",
													"PULSES"};
													
static const char *funciones_ataque[CANT_FUNCIONES_ATAQUE]={"LINEAR",
															"EXP",
															"QUARTSIN",
															"HALFSIN",
															"LOG",
															"TRI"};
															
static const char *funciones_sostenido_decaimiento[CANT_FUNCIONES_SOST_DEC]={ "CONSTANT",
																		"INVLINEAR",
																		"SIN",
																		"INVEXP",
																		"QUARTCOS",
																		"HALFCOS",
																		"INVLOG",
																		"PULSES"};

double f_constant(double t, double *t0, double *param2, double *param3)
{
	return 1;
}

double f_linear(double t, double *t0, double *param2, double *param3)
{
	return t/(double)(*t0);
}

double f_invlinear(double t, double *t0, double *param2, double *param3)
{
	if((1-(t/(*t0)))>0)
		return (1-(t/(*t0)));
	else
		return 0;
}

double f_sin(double t, double *a, double *f, double *param1)
{
	return (1+(*a)*sin((*f)*t));
}

double f_exp(double t, double *t0, double *param2, double *param3)
{
	return exp((5-(t-(*t0)))/(*t0));
}

double f_invexp(double t, double *t0, double *param2, double *param3)
{
	return exp((-5)*t/(*t0));
}

double f_quartcos(double t, double *t0, double *param2, double *param3)
{
	return cos(PI*t/2*(*t0));
}

double f_quartsin(double t, double *t0, double *param2, double *param3)
{
	return sin(PI*t/2*(*t0));
}

double f_halfcos(double t, double *t0, double *param2, double *param3)
{
	return ((1+cos(PI*t/(*t0)))/2);
}

double f_halfsin(double t, double *t0, double *param2, double *param3)
{
	return ((1+sin(PI*((t/(*t0))-((double)1/2))))/2);
}

double f_log(double t, double *t0, double *param2, double *param3)
{
	return log10((9*t/(*t0))+1);
}

double f_invlog(double t, double *t0, double *param2, double *param3)
{
	if(t<(*t0))
		return log10(((-9)*t/(*t0))+10);
	else
		return 0;
}

double f_tri(double t, double *t0, double *param2, double *param3)
{
	if(t < (*param2))
		return ((t*(*param3))/(*param2));

	return (((((t-(*param2))/((*param2)-(*t0)))*((*param3)-1)))+((*param3)));
}

double f_pulses(double t, double *t0, double *param2, double *param3)
{

	double tprima;
	double aux;
	
	tprima=(((t/(*t0))-floor(t/(*t0)))*(*t0));
	
	aux=(((1-(*param3))/(*param2))*(tprima-(*t0)+(*param2)));
	
	if(aux<0)
		aux*=-1;
		
	if((aux+(*param3))>1)
		return 1;
	else
		return (aux+(*param3));
}

double (*pfunciones[CANT_FUNCIONES])(double , double * , double *, double *)=
{	f_constant,
	f_linear,
	f_invlinear,
	f_sin,
	f_exp,
	f_invexp,
	f_quartcos,
	f_quartsin,
	f_halfcos,
	f_halfsin,
	f_log,
	f_invlog,
	f_tri,
	f_pulses};

status_t TDA_sintetizador_creadora(TDA_sintetizador_t **sint)
{
	if(sint == NULL)
		return ERROR_PUNTERO_NULO;

	
	*sint = (TDA_sintetizador_t *) malloc(sizeof(TDA_sintetizador_t));
	if(*sint == NULL)
		return ERROR_MEMORIA;
		
	return OK;
}

status_t TDA_sintetizador_inicializar(TDA_sintetizador_t *sint)
{
	size_t i;
	
	if(sint == NULL)
		return ERROR_PUNTERO_NULO;
		
	sint -> contador_lineas = 0;
	sint -> lineas = NULL;
	sint -> num_armonico = NULL;
	sint -> intensidad = NULL;
	for(i = 0; i < CANT_MODULADORAS; i++)
	{
		((sint -> moduladoras)[i]).nombre = NULL;
		((sint -> moduladoras)[i]).param1 = NULL;
		((sint -> moduladoras)[i]).param2 = NULL;
		((sint -> moduladoras)[i]).param3 = NULL;
	}
	return OK;
}
	
status_t validar_nombre_y_uso(char *nombre_en_archivo,size_t i)
{
	/*se valida que las funciones sean correctas y que esten usadas
	  correctamente, es decir, que sea efectivamente de decaimiento, ataque
	  o sostenido */
	bool_t encontrado = FALSE;
	size_t j;
									
	if(i == INDICE_ATAQUE)
	{
		for(j=0;(j<CANT_FUNCIONES_ATAQUE) && encontrado == FALSE;j++)
		{
			if(!strcmp(nombre_en_archivo,funciones_ataque[j]))
			{
				encontrado = TRUE;
				return OK;
			}	
		}
		return ERROR_NOMBRE_FUNCION;
	}
	else
	{
		for(j=0;(j<CANT_FUNCIONES_SOST_DEC) &&  encontrado == FALSE;j++)
		{
			if(!strcmp(nombre_en_archivo,funciones_sostenido_decaimiento[j]))
			{
				encontrado = TRUE;
				return OK;
			}	
		}
		return ERROR_NOMBRE_FUNCION;
	}
		
	return OK;
}

status_t TDA_sintetizador_cargar(TDA_sintetizador_t * sint, FILE * archivo_sint)
{
	size_t used, i, j;
	char **paux;
	
	if(sint == NULL)
		return ERROR_PUNTERO_NULO;
	
	/*crea arreglo de punteros a cadenas*/
	(sint-> lineas) = (char **) malloc (sizeof(char*) * INIT_CHOP_SINTETIZADOR);
	if ((sint -> lineas) == NULL)
		return ERROR_MEMORIA;
		
	used = INIT_CHOP_SINTETIZADOR;
	
	/*crea strings*/
	for(j = 0; j < used; j++)
	{
		(sint -> lineas) [j] = (char*) malloc (sizeof(char) * MAX_LINEA);
		if ((sint -> lineas)[j] == NULL)
			return ERROR_MEMORIA;
	}
	
	/*carga los strings*/
	i = 0;
	while(fgets((sint -> lineas)[i], MAX_LINEA, archivo_sint) != NULL)
	{
		i++;
		if (i == used)
		{
			paux = (char**) realloc((sint -> lineas), sizeof(char*)*(used + CHOP_SIZE_SINTETIZADOR));
			if (paux == NULL)
				return ERROR_MEMORIA;
				
			used += CHOP_SIZE_SINTETIZADOR;
			(sint -> lineas) = paux;
			for (j = i; j < used; j++)
			{
				(sint -> lineas)[j] = (char*) malloc (sizeof(char) * MAX_LINEA);
				if ((sint -> lineas)[j] == NULL)
					return ERROR_MEMORIA;
			}	
		}
	}
	sint -> contador_lineas = i;
	return OK;	
}

status_t TDA_sintetizador_cargar_cant_armonicos (TDA_sintetizador_t * sint)
{
	char *perr;
	
	sint -> cant_armonicos = strtoul((sint -> lineas)[0], &perr, 10);
	if (((*perr) != '\0') && ((*perr) != '\n'))
		return ERROR_FORMATO_SINTETIZADOR;
	
	return OK;
}

status_t TDA_sintetizador_liberar_memoria(TDA_sintetizador_t **sint)
{
	size_t i;
	
	if(sint == NULL)
		return ERROR_PUNTERO_NULO;
	
	free((*sint) -> num_armonico);
	free((*sint) -> intensidad);
	
	for(i = 0; i < CANT_MODULADORAS; i++)
	{
		free((*sint) -> moduladoras[i].nombre);
		
		if((*sint) -> moduladoras[i].param1 != NULL)
		{
			free((*sint) -> moduladoras[i].param1);
			(*sint) -> moduladoras[i].param1 = NULL;
		}
		
		if((*sint) -> moduladoras[i].param2 != NULL)
		{
			free((*sint) -> moduladoras[i].param2);
			(*sint) -> moduladoras[i].param2 = NULL;
		}
		
		if((*sint) -> moduladoras[i].param3 != NULL)
		{
			free((*sint) -> moduladoras[i].param3);
			(*sint) -> moduladoras[i].param3 = NULL;
		}
	}
	if ((*sint) -> lineas != NULL)
		for(i = 0; i < (*sint) -> contador_lineas; i++)
			free(((*sint) -> lineas)[i]);
	free((*sint) -> lineas);
	free(*sint);
	*sint = NULL;
	return OK;
}
	
status_t TDA_sintetizador_cargar_amplitudes (TDA_sintetizador_t * sint)
{
	char *perr, *aux_intensidad;
	size_t i;

	/*Leer cada armonico y su amplitud hasta sint->cant_armonicos*/
	
	sint -> num_armonico = (float*) malloc (sizeof(float) * (sint -> cant_armonicos));
	if(sint -> num_armonico == NULL)
		return ERROR_MEMORIA;
		
	sint->intensidad = (float*) malloc (sizeof(float) * (sint->cant_armonicos));
	if(sint -> intensidad == NULL)
		return ERROR_MEMORIA;
		
	for(i = 0; i < sint -> cant_armonicos; i++)
	{
		(sint -> num_armonico)[i] = strtod((sint -> lineas)[i+1], &perr);
		
		if (!isspace(*perr))
			return ERROR_FORMATO_SINTETIZADOR;
	
		aux_intensidad = perr+1; /* perr queda apuntando al espacio, hay que incrementar aux_intensidad para que
		 						 apunte al char siguiente */
		(sint->intensidad[i]) = strtod(aux_intensidad, &perr);
	}
	return OK;
}

status_t TDA_sintetizador_cargar_moduladoras(TDA_sintetizador_t * sint)
{
	size_t i, j, k;
	char *perr;

	for(i = (sint -> cant_armonicos) + 1, k = 0; i < CANT_MODULADORAS + (sint -> cant_armonicos) + 1; i++, k++)
	{
		j=0;
	
		/*Obtener el nombre*/
		while(!isspace((sint -> lineas)[i][j]) &&(sint ->lineas)[i][j] != '\0')
			j++;
			
		((sint->moduladoras)[k]).nombre = (char*) malloc((sizeof(char) * j) + 1);
		if(((sint -> moduladoras)[k]).nombre == NULL)
			return ERROR_MEMORIA;

		memcpy(	(((sint->moduladoras)[k]).nombre), (sint -> lineas)[i], j);
		(((sint -> moduladoras)[k]).nombre)[j] = '\0';

		if(validar_nombre_y_uso(sint -> moduladoras[k].nombre, k) != OK)
			return ERROR_NOMBRE_FUNCION;

		/*Obtener parametros */
		
		sint->moduladoras[k].param1=NULL;
		sint->moduladoras[k].param2=NULL;
		sint->moduladoras[k].param3=NULL;
		
		((sint->moduladoras)[k]).param1=(double *) malloc(sizeof(double));
		if(((sint -> moduladoras)[k]).param1 == NULL)
			return ERROR_MEMORIA;

		*(((sint->moduladoras)[k]).param1)=strtod(((sint -> lineas)[i])+j, &perr);
			
		if (!isspace(*perr) && (*perr) != '\0')
			return ERROR_FORMATO_SINTETIZADOR;

		if(isspace(*perr))
		{
			((sint->moduladoras)[k]).param2=(double *) malloc(sizeof(double));
			if(((sint -> moduladoras)[k]).param2 == NULL)
				return ERROR_MEMORIA;
				
			*(sint->moduladoras[k].param2)=strtod(perr+1, &perr);	
		
			if (!isspace(*perr) && (*perr) != '\0')
				return ERROR_FORMATO_SINTETIZADOR;
		
			if(isspace(*perr))
			{
				((sint->moduladoras)[k]).param3=(double *) malloc(sizeof(double));
				
				if(((sint -> moduladoras)[k]).param1 == NULL)
					return ERROR_MEMORIA;
	
				*(sint->moduladoras[k].param3)=strtod(perr+1,&perr);
				
				if (!isspace(*perr) && (*perr) != '\0')
					return ERROR_FORMATO_SINTETIZADOR;
				
			}
		}
	}
	return OK;
}

funciones_t obtener_posicion_de_funcion(const TDA_sintetizador_t sint, size_t pos)
{
	/*La funcion se usa post-validacion de los nombres de las cadenas, por lo tanto
	  el strcmp deberia ser 0 en algun momento y se asigna el valor correspondiente.
	  Devuelve un valor que esta asociado al arreglo de punteros a funcion*/
	  
	size_t j=0;
	bool_t encontrado = FALSE;
	funciones_t num_funcion;
	
	for(j=0;(j<CANT_FUNCIONES)&&!encontrado;j++)
	{
		if(!strcmp(sint.moduladoras[pos].nombre,arr_funciones[j]))
		{
			encontrado=TRUE;
			num_funcion=j;
		}
	}
	return num_funcion;
}

status_t asignar_funciones_moduladoras(TDA_sintetizador_t * sint)
{
	funciones_t posicion[CANT_MODULADORAS];
	size_t i=0;
	
	while(i<CANT_MODULADORAS)
	{
		posicion[i]=obtener_posicion_de_funcion(*sint, i);
		i++;
	}
	
	sint->pataque=(*pfunciones[posicion[INDICE_ATAQUE]]);
	sint->psostenido=(*pfunciones[posicion[INDICE_SOSTENIDO]]);
	sint->pdecaimiento=(*pfunciones[posicion[INDICE_DECAIMIENTO]]);
	
	return OK;
}

