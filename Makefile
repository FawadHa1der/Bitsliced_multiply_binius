# Makefile

# Compiler
CC = clang

# Compiler Flags
# CFLAGS = -Wall -O3 -ftree-vectorize
CFLAGS = -Wall -O2 -ftree-vectorize
# CFLAGS = -Wall

# Output executable name
OUTPUT = bitsliced

# Output static library name
LIBRARY = libbitsliced.a

# Source files (excluding test.c for static lib)
SOURCES = bs_multiply_64.c bs.c bs_multiply_128.c polynomial_byte_sliced_mul_16.c polynomial_byte_slice_mul_2.c

# Add test.c to SOURCES if not building static library
ifeq ($(STATIC),)
SOURCES += test.c
endif

# Object files
OBJECTS = $(SOURCES:.c=.o)

# Header files
HEADERS = bs_multiply_64.h bs.h bs_multiply_128.h polynomial_byte_sliced_mul_16.h polynomial_byte_slice_mul_2.h

# Build all
all: $(OUTPUT)

# Conditional static library build
ifeq ($(STATIC),1)
all: $(LIBRARY)
endif

# Link to create executable
$(OUTPUT): $(OBJECTS)
	$(CC) $(CFLAGS) -o $(OUTPUT) $(OBJECTS)

# Create the static library
$(LIBRARY): $(OBJECTS)
	# Use ar to create a static library
	ar rcs $(LIBRARY) $(OBJECTS)

# Compile object files
%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -c $<

# Clean
clean:
	rm -f $(OBJECTS) $(OUTPUT) $(LIBRARY)