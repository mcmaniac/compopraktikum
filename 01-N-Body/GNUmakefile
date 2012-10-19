OUT  = nbody
OBJS = vector.o io.o nbody.o
HDRS = nbody.h

CC 	 = gcc
ODIR = obj
SDIR = src

_OBJS = $(OBJS:%=$(ODIR)/%)
_HDRS = $(HDRS:%=$(SDIR)/%)

$(ODIR)/%.o: $(SDIR)/%.c $(_HDRS)
	$(CC) -c $(CFLAGS) -o $@ $<

$(OUT): $(_OBJS)
	$(CC) -o $(OUT) $(_OBJS)

.PHONY: all clean
all: $(OUT)
clean:
	rm -f $(ODIR)/*.o $(OUT)
