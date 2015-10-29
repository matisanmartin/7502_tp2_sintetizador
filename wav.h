#ifndef WAV_H
#define WAV_H 

#define SIZE 4
#define SIZE_32_BITS 4
#define SIZE_16_BITS 2

#define MASK 0xFF
#define SHIFT 8

#define CHUNK_SIZE_36 36
#define CHUNK_ID "RIFF"
#define FORMAT "WAVE"
#define SUBCHUNK_ID "fmt "
#define SUBCHUNK1_SIZE 16
#define AUDIO_FORMAT 1
#define NUM_CHANNELS 1
#define BLOCK_ALIGN 2
#define BITS_PER_SAMPLE 16
#define SUBCHUNK2_ID "data"

#define MSJ_ERROR_CARGAR_WAV "No se pudo cargar el WAV\n"
#define MSJ_ERROR_CREAR_ARCHIVO_WAV "No se pudo crear el archivo WAV\n"
#define MSJ_ERROR_CONVERTIR "No se pudo convertir al formato de escritura\n"
#define MSJ_ERROR_INICIALIZAR_WAV "No se pudo inicializar el WAV"
#define MSJ_ERROR_CREAR_WAV "No se pudo crear el wav"

struct TDA_wav{
			char chunk_id[SIZE];
			unsigned int chunk_size; /*32 bytes*/
			char format[SIZE];
			char subchunk_id[SIZE];
			unsigned int subchunk1_size; /*32 bytes*/
			unsigned int audio_format; /*bytes 16*/
			unsigned int num_channels; /*16 bytes*/
			unsigned int sample_rate; /*32 bytes*/
			unsigned int byte_rate; /*32 bytes*/
			unsigned int block_align; /*16 bytes*/
			unsigned int bits_per_sample; /*16 bytes*/
			char subchunk2_id [SIZE];
			unsigned int subchunk2_size; /*32 bytes*/
			float *data; /*16 bytes*/
			};
			
typedef struct TDA_wav TDA_wav_t;

status_t TDA_wav_creadora(TDA_wav_t **wav);
status_t TDA_wav_inicializar(TDA_wav_t * wav, const uint cantidad_muestras,
							const uint frecuencia_muestreo);
status_t convertir_a_uint32(uint valor, char aux[]);
status_t convertir_a_uint16(uint valor, char aux[]);
status_t crear_archivo_wav(TDA_wav_t *wav, const uint cantidad_muestras, FILE *archivo);
status_t TDA_wav_liberar_memoria(TDA_wav_t **wav);

#endif
