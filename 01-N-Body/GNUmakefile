OUT  = nbody

CC 	 = gcc
ODIR = obj
SDIR = src
LIBS = -lm
CFLAGS = -Wall -std=c99

SRCS   = $(wildcard $(SDIR)/*.c)
HDRS   = $(wildcard $(SDIR)/*.h)
OBJS   = $(SRCS:$(SDIR)/%.c=$(ODIR)/%.o)

.PHONY: all clean
all: $(OUT)

clean:
	rm -f $(ODIR)/*.o $(OUT)

$(OUT): $(OBJS)
	$(CC) $(LIBS) -o $(OUT) $(OBJS)

$(ODIR)/%.o: $(SDIR)/%.c $(HDRS)
	$(CC) -c -o $@ $< $(CFLAGS)
