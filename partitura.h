#ifndef PARTITURA_H
#define PARTITURA_H 

#define MSJ_ERROR_CALCULAR_DURACION "No se pudo calcular la duracion de la nota"
#define MSJ_ERROR_CALCULAR_FRECUENCIAS "No se pudo calcular la frecuencia de la nota"
#define MSJ_ERROR_CREAR_PARTITURA "No se pudo crear la partitura"
#define MSJ_ERROR_INICIALIZAR_PARTITURA "No se pudo inicializar partitura"
#define MSJ_ERROR_CARGAR_PARTITURA "No se pudo cargar partitura"
#define MSJ_ERROR_CARGAR_INICIO_PARTITURA "No se pudo cargar el inicio de la partitura"
#define MSJ_ERROR_CARGAR_DURACION "No se pudo cargar la duracion"

#define FREQ_C 65.4064
#define FREQ_CS 69.2957
#define FREQ_D 73.4162
#define FREQ_DS 77.7917
#define FREQ_E 82.4069
#define FREQ_F 87.3071
#define FREQ_FS 92.4986
#define FREQ_G 97.9989
#define FREQ_GS 103.8262
#define FREQ_A 110.0000
#define FREQ_AS 123.4708
#define FREQ_B 130.8128

#define DO 'C'
#define RE 'D'
#define MI 'E'
#define FA 'F'
#define SOL 'G'
#define LA 'A'
#define SI 'B'
#define SOSTENIDO 's'
#define BEMOL 'b'
#define PRIMER_ELEMENTO_ARREGLO 0
#define SEGUNDO_ELEMENTO_ARREGLO 1

#define INIT_SIZE_PARTITURA 100
#define CHOP_SIZE_PARTITURA 2	

#define NUMERO_O_ESPACIO " .1234567890"

struct TDA_partitura{
	size_t contador_notas;
	float duracion_partitura;
	float *aux_duracion_nota;
	float *aux_inicio_nota;
	float *aux_frecuencia_nota;
	char **lineas;
};

typedef struct TDA_partitura TDA_partitura_t; 

typedef enum {NOTA_DO = 0, NOTA_DO_S = 1, NOTA_RE = 2, NOTA_RE_S = 3, NOTA_MI = 4, NOTA_FA = 5, NOTA_FA_S = 6,
			NOTA_SOL = 7, NOTA_SOL_S = 8, NOTA_LA = 9, NOTA_LA_S = 10, NOTA_SI = 11} notas_t;


float obtener_frecuencia_nota (const char nota_americana[]);
status_t TDA_partitura_creadora(TDA_partitura_t **partitura);
status_t TDA_partitura_inicializar(TDA_partitura_t *partitura);
status_t TDA_partitura_cargar (TDA_partitura_t *partitura, FILE *archivo_partitura);
status_t TDA_partitura_cargar_inicio (TDA_partitura_t *partitura);
status_t TDA_partitura_cargar_duracion (TDA_partitura_t *partitura);
status_t TDA_partitura_calcular_duracion (TDA_partitura_t *partitura);
status_t TDA_partitura_calcular_frecuencias (TDA_partitura_t *partitura);
status_t TDA_partitura_liberar_memoria(TDA_partitura_t **partitura);

#endif
