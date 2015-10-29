#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include "common.h"
#include "partitura.h"
#include "wav.h"

status_t TDA_wav_creadora(TDA_wav_t **wav)
{
	if(wav == NULL)
		return ERROR_PUNTERO_NULO;
	
	*wav = (TDA_wav_t *) malloc(sizeof(TDA_wav_t));
	if(*wav == NULL)
		return ERROR_MEMORIA;
		
	return OK;
}

status_t TDA_wav_inicializar(TDA_wav_t * wav, const uint cantidad_muestras,
							const uint frecuencia_muestreo)
{
	
	if(wav==NULL)
		return ERROR_PUNTERO_NULO;
	
	if(strncpy(wav->chunk_id,CHUNK_ID,SIZE)==NULL)
		return ERROR_CARGAR_WAV;
		
	if(strncpy(wav->format,FORMAT,SIZE)==NULL)
		return ERROR_CARGAR_WAV;
		
	if(strncpy(wav->subchunk_id,SUBCHUNK_ID,SIZE)==NULL)
		return ERROR_CARGAR_WAV;
		
	if(strncpy(wav->subchunk2_id,SUBCHUNK2_ID,SIZE)==NULL)
		return ERROR_CARGAR_WAV;
		
	wav->chunk_size=CHUNK_SIZE_36 + 2*cantidad_muestras;
	wav->subchunk1_size=SUBCHUNK1_SIZE;
	wav->audio_format=AUDIO_FORMAT;
	wav->num_channels=NUM_CHANNELS;
	wav->sample_rate=frecuencia_muestreo;
	wav->byte_rate=(2*frecuencia_muestreo);
	wav->block_align=BLOCK_ALIGN;
	wav->bits_per_sample=BITS_PER_SAMPLE;
	wav->subchunk2_size=2*cantidad_muestras;
	
	wav -> data = NULL;
	
	return OK;
}

status_t TDA_wav_liberar_memoria(TDA_wav_t **wav)
{
	if(wav == NULL)
		return ERROR_PUNTERO_NULO;
	
	free((*wav) -> data);
	(*wav) -> data = NULL;
	free(*wav);
	*wav = NULL;
	return OK;
}
status_t convertir_a_uint32(uint valor, char aux[])
{
	aux[0]=valor&MASK;
	aux[1]=(valor>>SHIFT)&MASK;
	aux[2]=(valor>>(2*SHIFT))&MASK;
	aux[3]=(valor>>(3*SHIFT))&MASK;
	return OK;
}

status_t convertir_a_uint16(uint valor, char aux[])
{
	aux[0]=valor&MASK;
	aux[1]=(valor>>SHIFT)&MASK;
	return OK;
}

status_t crear_archivo_wav(TDA_wav_t *wav, const uint cantidad_muestras, FILE *archivo)
{
	char aux_32[SIZE];
	char aux_16[SIZE/2];
	size_t i;
	short int aux;
	
	fwrite(wav -> chunk_id,sizeof(char),SIZE,archivo);

	if(convertir_a_uint32((uint)(wav -> chunk_size), aux_32) != OK)
		return ERROR_CONVERTIR;
	
	fwrite(aux_32, sizeof(char), SIZE, archivo);
	
	fwrite(wav -> format, sizeof(char), SIZE, archivo);
	
	fwrite(wav -> subchunk_id, sizeof(char), SIZE, archivo);
	
	if(convertir_a_uint32((uint) (wav -> subchunk1_size), aux_32) != OK)
		return ERROR_CONVERTIR;
	
	fwrite(aux_32, sizeof(char), SIZE, archivo);
	
	if(convertir_a_uint16((uint) (wav -> audio_format), aux_16) != OK)
		return ERROR_CONVERTIR;
	
	fwrite(aux_16, sizeof(char), SIZE/2, archivo);
	
	if(convertir_a_uint16((uint) (wav -> num_channels), aux_16) != OK)
		return ERROR_CONVERTIR;
	
	fwrite(aux_16,sizeof(char),SIZE/2,archivo);
	
	if(convertir_a_uint32((uint) (wav -> sample_rate),aux_32)!=OK)
		return ERROR_CONVERTIR;
	
	fwrite(aux_32,sizeof(char),SIZE,archivo);
	
	if(convertir_a_uint32((uint)(wav -> byte_rate),aux_32)!=OK)
		return ERROR_CONVERTIR;
		
	fwrite(aux_32,sizeof(char),SIZE,archivo);
	
	if(convertir_a_uint16((uint)(wav -> block_align),aux_16)!=OK)
		return ERROR_CONVERTIR;
	
	fwrite(aux_16,sizeof(char),SIZE/2,archivo);
	
	if(convertir_a_uint16((uint)(wav -> bits_per_sample),aux_16)!=OK)
		return ERROR_CONVERTIR;
	
	fwrite(aux_16,sizeof(char),SIZE/2,archivo);
	
	fwrite(wav -> subchunk2_id, sizeof(char), 4, archivo);
	
	if(convertir_a_uint32((uint)(wav -> subchunk2_size),aux_32)!=OK)
		return ERROR_CONVERTIR;
	
	fwrite(aux_32,sizeof(char),SIZE,archivo);
	
	for(i = 0; i < cantidad_muestras; i++)
	{
		aux = (short int) (wav -> data)[i];
		fwrite(&aux,SIZE/2,1,archivo);
	}
	return OK;
}
