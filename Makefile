# Compiler and compilation flags
CC = gcc
CFLAGS = -ansi -Wall -Wextra -Werror -pedantic-errors
LINKER_FLAGS = -lm

# Executable name
EXEC = symnmf

# Source files
SRC = symnmf.c
HEADERS = symnmf.h

# Object files (compiled .o files)
OBJ = $(SRC:.c=.o)

# Default target: Compile the program
all: $(EXEC)

$(EXEC): $(OBJ)
	$(CC) $(CFLAGS) -o $(EXEC) $(OBJ) $(LINKER_FLAGS)

# Compile .c files into .o files
%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

# Clean target: Remove compiled files
clean:
	rm -f $(EXEC) $(OBJ)
