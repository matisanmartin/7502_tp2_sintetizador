sintetizador: main.o validar.o partitura.o sintetizador.o muestreo.o wav.o
	gcc -o sintetizador main.o validar.o partitura.o sintetizador.o muestreo.o wav.o -lm

main.o: main.c common.h partitura.h sintetizador.h validar.h wav.h
	gcc -c main.c -Wall -ansi -pedantic

validar.o: validar.c validar.h common.h
	gcc -c validar.c -Wall -ansi -pedantic

partitura.o: partitura.c partitura.h common.h
	gcc -c -g partitura.c -Wall -ansi -pedantic

sintetizador.o: sintetizador.c sintetizador.h common.h
	gcc -c sintetizador.c -Wall -ansi -pedantic

muestreo.o: muestreo.c common.h muestreo.h sintetizador.h
	gcc -c muestreo.c -Wall -ansi -pedantic

wav.o: wav.c wav.h common.h
	gcc -c wav.c -Wall -ansi -pedantic
