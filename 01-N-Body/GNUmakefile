OUT  = nbody
OBJS = vector.o io.o nbody.o integrators.o constants.o
HDRS = nbody.h integrators.h

CC 	 = gcc
ODIR = obj
SDIR = src
LIBS = -lmath

_OBJS = $(OBJS:%=$(ODIR)/%)
_HDRS = $(HDRS:%=$(SDIR)/%)

$(ODIR)/%.o: $(SDIR)/%.c $(_HDRS)
	$(CC) -c -o $@ $< $(CFLAGS)

$(OUT): $(_OBJS)
	$(CC) -o $(OUT) $(_OBJS)

.PHONY: all clean
all: $(OUT)
clean:
	rm -f $(ODIR)/*.o $(OUT)
