#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "common.h"
#include "validar.h"
#include "partitura.h" 
#include "sintetizador.h"
#include "wav.h"
#include "muestreo.h"

int main (int argc, char *argv[])
{	
	TDA_partitura_t *partitura;
	TDA_sintetizador_t *sintetizador;
	TDA_wav_t *wav;
	FILE *f_partitura, *f_sintetizador, *f_audio;
	uint frecuencia, cantidad_muestras;
	float *muestras;
	
	if(validar(argc, argv, &frecuencia, &f_partitura, &f_sintetizador, &f_audio) != OK) 
	{
		fprintf(stderr,"%s\n",MSJ_ERROR_ARGUMENTOS);
		return EXIT_FAILURE;
	}
	
	if(TDA_partitura_creadora(&partitura) != OK)
	{
		fprintf(stderr,"%s\n",MSJ_ERROR_CREAR_PARTITURA);
		return EXIT_FAILURE;
	}
	

	if(TDA_partitura_inicializar(partitura) != OK)
	{
		fprintf(stderr,"%s\n",MSJ_ERROR_INICIALIZAR_PARTITURA);
		TDA_partitura_liberar_memoria(&partitura);
		return EXIT_FAILURE;
	}
	

	if(TDA_partitura_cargar(partitura, f_partitura) != OK) 
	{
		fprintf(stderr,"%s\n",MSJ_ERROR_CARGAR_PARTITURA);
		TDA_partitura_liberar_memoria(&partitura);
		return EXIT_FAILURE;
	}
	
	if(TDA_partitura_cargar_inicio (partitura) != OK) 
	{
		fprintf(stderr,"%s\n",MSJ_ERROR_CARGAR_INICIO_PARTITURA);
		TDA_partitura_liberar_memoria(&partitura);
		return EXIT_FAILURE;
	}
	
	if(TDA_partitura_cargar_duracion (partitura) != OK) 
	{
		fprintf(stderr,"%s\n",MSJ_ERROR_CARGAR_DURACION);
		TDA_partitura_liberar_memoria(&partitura);
		return EXIT_FAILURE;
	}
	
	if(TDA_partitura_calcular_duracion (partitura) != OK) 
	{
		fprintf(stderr,"%s\n",MSJ_ERROR_CALCULAR_DURACION);
		TDA_partitura_liberar_memoria(&partitura);
		return EXIT_FAILURE;
	}
	
	if(TDA_partitura_calcular_frecuencias (partitura) != OK) 
	{
		fprintf(stderr,"%s\n",MSJ_ERROR_CALCULAR_FRECUENCIAS);
		TDA_partitura_liberar_memoria(&partitura);
		return EXIT_FAILURE;
	}
	
	if(TDA_sintetizador_creadora(&sintetizador) != OK)	
	{
		fprintf(stderr,"%s\n",MSJ_ERROR_CREAR_SINTETIZADOR);
		TDA_partitura_liberar_memoria(&partitura);
		return EXIT_FAILURE;
	}
	
	if((TDA_sintetizador_inicializar(sintetizador) != OK))	
	{
		fprintf(stderr,"%s\n",MSJ_ERROR_INICIALIZAR_SINTETIZADOR);
		TDA_sintetizador_liberar_memoria(&sintetizador);
		TDA_partitura_liberar_memoria(&partitura);
		return EXIT_FAILURE;
	}
	
	if(TDA_sintetizador_cargar(sintetizador,f_sintetizador) != OK)	
	{
		fprintf(stderr,"%s\n",MSJ_ERROR_CARGAR_SINTETIZADOR);
		TDA_sintetizador_liberar_memoria(&sintetizador);
		TDA_partitura_liberar_memoria(&partitura);
		return EXIT_FAILURE;
	}
	
	if(TDA_sintetizador_cargar_cant_armonicos (sintetizador) != OK) 
	{
		fprintf(stderr,"%s\n",MSJ_ERROR_CARGAR_CANT_ARMONICOS);
		TDA_sintetizador_liberar_memoria(&sintetizador);
		TDA_partitura_liberar_memoria(&partitura);
		return EXIT_FAILURE;
	}
	
	if(TDA_sintetizador_cargar_amplitudes (sintetizador) != OK)	
	{
		fprintf(stderr,"%s\n",MSJ_ERROR_CARGAR_AMPLITUDES);
		TDA_sintetizador_liberar_memoria(&sintetizador);
		TDA_partitura_liberar_memoria(&partitura);
		return EXIT_FAILURE;
	}
	
	if(TDA_sintetizador_cargar_moduladoras(sintetizador) != OK)	
	{
		fprintf(stderr,"%s\n",MSJ_ERROR_CARGAR_MODULADORAS);
		TDA_sintetizador_liberar_memoria(&sintetizador);
		TDA_partitura_liberar_memoria(&partitura);
		return EXIT_FAILURE;
	}
	
	if(asignar_funciones_moduladoras(sintetizador) != OK)	
	{
		fprintf(stderr,"%s\n",MSJ_ERROR_ASIGNAR_MODULADORAS);
		TDA_sintetizador_liberar_memoria(&sintetizador);
		TDA_partitura_liberar_memoria(&partitura);
		return EXIT_FAILURE;
	}
	
	cantidad_muestras=(unsigned int)((partitura->duracion_partitura + *((sintetizador->moduladoras[INDICE_DECAIMIENTO]).param1)) * frecuencia) +1;
	
	if(crear_vector_muestras(&muestras, cantidad_muestras) != OK)
	{
		fprintf(stderr,"%s\n",MSJ_ERROR_CARGAR_MUESTRAS);
		TDA_partitura_liberar_memoria(&partitura);
		TDA_sintetizador_liberar_memoria(&sintetizador);
		return EXIT_FAILURE;
	}
	
	if(muestrear(partitura -> contador_notas, partitura -> aux_frecuencia_nota, partitura -> aux_inicio_nota,
				partitura -> aux_duracion_nota, sintetizador -> cant_armonicos, sintetizador -> intensidad,
				sintetizador -> num_armonico, 1.0/frecuencia, frecuencia, muestras,
				sintetizador) != OK)
	{
		fprintf(stderr,"%s\n",MSJ_ERROR_MUESTREO);
		TDA_partitura_liberar_memoria(&partitura);
		TDA_sintetizador_liberar_memoria(&sintetizador);
		liberar_memoria_muestras(&muestras);
		return EXIT_FAILURE;
	}
	if(multiplicar_factor_amplitud(muestras, cantidad_muestras) != OK)
	{
		fprintf(stderr,"%s\n",MSJ_ERROR_FACTOR_AMPLITUD);
		TDA_partitura_liberar_memoria(&partitura);
		TDA_sintetizador_liberar_memoria(&sintetizador);
		liberar_memoria_muestras(&muestras);
		return EXIT_FAILURE;
	}
	
	if(TDA_wav_creadora(&wav) != OK)
	{
		fprintf(stderr,"%s\n",MSJ_ERROR_CREAR_WAV);
		TDA_partitura_liberar_memoria(&partitura);
		TDA_sintetizador_liberar_memoria(&sintetizador);
		liberar_memoria_muestras(&muestras);
		return EXIT_FAILURE;
	}
	
	if(TDA_wav_inicializar(wav, cantidad_muestras, frecuencia) != OK)
	{
		fprintf(stderr,"%s\n",MSJ_ERROR_INICIALIZAR_WAV);
		TDA_wav_liberar_memoria(&wav);
		TDA_partitura_liberar_memoria(&partitura);
		TDA_sintetizador_liberar_memoria(&sintetizador);
		liberar_memoria_muestras(&muestras);
		return EXIT_FAILURE;
	}
	
	wav -> data = muestras;

	if(crear_archivo_wav(wav, cantidad_muestras, f_audio) != OK)
	{
		fprintf(stderr,"%s\n",MSJ_ERROR_CREAR_ARCHIVO_WAV);
		TDA_wav_liberar_memoria(&wav);
		TDA_partitura_liberar_memoria(&partitura);
		TDA_sintetizador_liberar_memoria(&sintetizador);
		return EXIT_FAILURE;
	}
	
	TDA_partitura_liberar_memoria(&partitura);
	TDA_sintetizador_liberar_memoria(&sintetizador);
	TDA_wav_liberar_memoria(&wav);
	
	fclose(f_partitura); 
	fclose(f_sintetizador);
	/*no hace falta validarlos porque se abrieron en modo lectura*/
	
	if(fclose(f_audio) == EOF) 
		return EXIT_FAILURE;
	
	return EXIT_SUCCESS;
}
