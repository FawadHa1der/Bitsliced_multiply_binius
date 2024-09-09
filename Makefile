# Makefile

# Compiler
CC = clang

# Compiler Flags
CFLAGS = -Wall -O3 -ftree-vectorize
# CFLAGS = -Wall 

# Output executable name
OUTPUT = bitsliced

# Source files
SOURCES = test.c bs_multiply_64.c bs.c 

# Object files
OBJECTS = $(SOURCES:.c=.o)

# Header files
HEADERS = bs_multiply_64.h bs.h

# Build all
all: $(OUTPUT)

# Link
$(OUTPUT): $(OBJECTS)
	$(CC) $(CFLAGS) -o $(OUTPUT) $(OBJECTS)

# Compile
%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -c  $<

# Clean
clean:
	rm -f $(OBJECTS) $(OUTPUT)
