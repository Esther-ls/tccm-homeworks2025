CC=gcc
SRC=functions.c main.c
OBJ=$(patsubst %.c, %.o, $(SRC))
TREXIO=-ltrexio -lm


main: $(OBJ)
	$(CC) -o main $(OBJ) $(TREXIO)

%.o: %.c
	$(CC) -c $<

clean:
	rm -f $(OBJ) main
