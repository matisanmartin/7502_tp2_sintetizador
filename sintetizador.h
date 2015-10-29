#ifndef SINTETIZADOR_H
#define SINTETIZADOR_H 

#define INIT_CHOP_SINTETIZADOR 30
#define CHOP_SIZE_SINTETIZADOR 5

#define CANT_MODULADORAS 3
#define CANT_FUNCIONES 14
#define CANT_FUNCIONES_ATAQUE 6
#define CANT_FUNCIONES_SOST_DEC 8
#define CANT_PARAMETROS 2

#define MSJ_ERROR_CREAR_SINTETIZADOR "No se pudo crear el sintetizador"
#define MSJ_ERROR_INICIALIZAR_SINTETIZADOR "No se pudo inicializar el sintetizador"
#define MSJ_ERROR_CARGAR_SINTETIZADOR "No se pudo cargar el sintetizador"
#define MSJ_ERROR_CARGAR_CANT_ARMONICOS "No se pudo cargar la cantidad de armonicos"
#define MSJ_ERROR_CARGAR_AMPLITUDES "No se pudieron cargar las amplitudes"
#define MSJ_ERROR_CARGAR_MODULADORAS "No se pudieron cargar las moduladoras"
#define MSJ_ERROR_ASIGNAR_MODULADORAS "No se pudieron asignar las moduladoras"

typedef enum{INDICE_ATAQUE,INDICE_SOSTENIDO,INDICE_DECAIMIENTO}indice_t;

typedef enum{CONSTANT,LINEAR,INVLINEAR,SIN,EXP,INVEXP,QUARTCOS,QUARTSIN,
			 HALFCOS,HALFSIN,LOG,INVLOG,TRI,PULSES}funciones_t;
			 
typedef enum{FALSE,TRUE}bool_t;

struct moduladora{
	char *nombre;
	double *param1;
	double *param2;
	double *param3;
};

typedef struct TDA_sintetizador TDA_sintetizador_t;
typedef struct moduladora moduladora_t;


struct TDA_sintetizador{
	char **lineas;
	unsigned int cant_armonicos;
	float *num_armonico;
	float *intensidad;
	moduladora_t moduladoras[CANT_MODULADORAS];
	double (*pataque)(double , double *, double *, double *);
	double (*psostenido)(double , double *, double *, double *);
	double (*pdecaimiento)(double , double *, double *, double *);
	size_t contador_lineas;
};

double f_constant(double t, double *t0, double *param2, double *param3);
double f_linear(double t, double *t0, double *param2, double *param3);
double f_invlinear(double t, double *t0, double *param2, double *param3);
double f_sin(double t, double *a, double *f, double *param1);
double f_exp(double t, double *t0, double *param2, double *param3);
double f_invexp(double t, double *t0, double *param2, double *param3);
double f_quartcos(double t, double *t0, double *param2, double *param3);
double f_quartsin(double t, double *t0, double *param2, double *param3);
double f_halfcos(double t, double *t0, double *param2, double *param3);
double f_halfsin(double t, double *t0, double *param2, double *param3);
double f_log(double t, double *t0, double *param2, double *param3);
double f_invlog(double t, double *t0, double *param2, double *param3);
double f_tri(double t, double *t0, double *a1, double *t1);
double f_pulses(double t, double *t0, double *param2, double *param3);

status_t TDA_sintetizador_creadora(TDA_sintetizador_t **sint);
status_t TDA_sintetizador_inicializar(TDA_sintetizador_t *sint);
status_t TDA_sintetizador_cargar(TDA_sintetizador_t * sint, FILE * archivo_sint);
status_t validar_nombre_y_uso(char *nombre_en_archivo,size_t i);
status_t TDA_sintetizador_cargar_cant_armonicos (TDA_sintetizador_t * sint);
status_t TDA_sintetizador_liberar_memoria(TDA_sintetizador_t **sint);
status_t TDA_sintetizador_cargar_amplitudes (TDA_sintetizador_t * sint);
status_t TDA_sintetizador_cargar_moduladoras(TDA_sintetizador_t * sint);
funciones_t obtener_posicion_de_funcion(const TDA_sintetizador_t sint, size_t i);
status_t asignar_funciones_moduladoras(TDA_sintetizador_t * sint);

#endif
