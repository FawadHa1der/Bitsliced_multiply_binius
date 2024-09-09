# Makefile

# Compiler
CC = clang

# Compiler Flags
CFLAGS = -Wall -O3 -ftree-vectorize
# CFLAGS = -Wall 

# Output executable name
OUTPUT = unbitsliced

# Source files
SOURCES = binmul.c

# Object files
OBJECTS = $(SOURCES:.c=.o)

# # Header files
HEADERS = ./*.h

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
