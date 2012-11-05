# name of the executable
EXEC = ssg
LIBS = -lm

#
# do not modify after this
#
CC     = gcc
CFLAGS = -Wall -std=c99 -pedantic

ODIR   = obj
SDIR   = src
SRCS   = $(wildcard $(SDIR)/*.c)
HDRS   = $(wildcard $(SDIR)/*.h)
OBJS   = $(SRCS:$(SDIR)/%.c=$(ODIR)/%.o)
DEPEND = $(ODIR)/GNUmakefile.dep

.PHONY: all clean

all: $(DEPEND) $(EXEC)

clean:
	rm -rf $(ODIR) $(EXEC) $(DEPEND)

# rule to build/link executable
$(EXEC): $(OBJS)
	$(CC) $(LIBS) -o $(EXEC) $(OBJS)

# rule to build .o files
$(OBJS):
	$(CC) -c -o $@ $< $(CFLAGS)

# make sure $(ODIR) exists
$(ODIR):
	mkdir $(ODIR)

# generate dependency makefile for all source files
$(DEPEND): $(ODIR)
	$(CC) $(CFLAGS) -MM $(SRCS) $(HDRS) | sed -E 's,^[[:alpha:]],$(ODIR)/&,' > $(DEPEND)

-include $(DEPEND)
