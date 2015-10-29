#include "sintetizador.h"

#ifndef MUESTREO_H
#define MUESTRO_H

#define MSJ_ERROR_CARGAR_MUESTRAS "No se pudieron cargar las muestras"
#define MSJ_ERROR_MUESTREO "No se pudo muestrear"
#define MSJ_ERROR_FACTOR_AMPLITUD "No se pudo multiplicar por el factor de amplitud"

#define INT_LIMIT 16384
/*pow(2, BITS_PER_SAMPLE) * FACTOR   */

status_t crear_vector_muestras (float **muestras,const uint cantidad_muestras);
status_t muestrear (const uint cant_notas, const float notas[], const float inicios[], const float duraciones[],
					const uint cant_armonicos, const float intensidades[], const float armonicos[],
					const float intervalo, const uint frecuencia, float muestras[],
					TDA_sintetizador_t *sint);					
status_t obtener_amplitud(const float delta_tiempo,const float duracion_nota,
						TDA_sintetizador_t *sint, float *muestra);
status_t multiplicar_factor_amplitud(float *muestras, const uint cant_muestras);
status_t obtener_factor_amplitud(float *factor_amplitud, const float *muestras, const uint cant_muestras);
status_t obtener_maximo(float *max, const float *muestras, const uint cant_muestras);
status_t obtener_minimo(float *min, const float *muestras, const uint cant_muestras);
status_t liberar_memoria_muestras(float **muestras);

#endif 
