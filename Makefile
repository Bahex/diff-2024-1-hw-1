.PHONY: clean

STD := -std=c99

# Specifies to GCC the required warnings
WARNS := -Wall -Wextra -Wpedantic -Wno-unused-function

# Debug options
DEBUG := -g3 -DDEBUG=1

OPTS := -O3

# Flags for compiling
CFLAGS := $(OPTS) $(STD) $(WARNS)

# Dependency libraries
CFLAGS += $(STD) $(DEBUG)
LDLIBS += -lm

main: main.c $(wildcard sensible/**/*.{h,c})

clean:
	rm -f main
