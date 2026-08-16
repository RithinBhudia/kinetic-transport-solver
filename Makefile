CC = gcc
CFLAGS = -Wall -Wextra -O2
LDLIBS = -llapack -lblas -lm

TARGET = solver
SRC = src/main.c src/config.c src/solver.c

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) $(SRC) -o $(TARGET) $(LDLIBS)

clean:
	rm -f $(TARGET) src/*.o