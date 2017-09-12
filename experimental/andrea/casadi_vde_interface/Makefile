CC = gcc
CFLAGS = -O3
main: 	main.o model.o vde_forw_model.o
	$(CC) -o main main.o model.o vde_forw_model.o -lm -I.

clean:
	rm main.o model.o vde_forw_model.o
