# name of the executable
LIB = libcompo.a

LIBF = -lm

#
# do not modify after this
#
CC     = gcc
CFLAGS = -Wall -std=c99 -pedantic
AR     = ar

ODIR   = obj
SDIR   = src
LDIR   = lib
SRCS   = $(wildcard $(SDIR)/*.c)
HDRS   = $(wildcard $(SDIR)/*.h)
OBJS   = $(SRCS:$(SDIR)/%.c=$(ODIR)/%.o)
LIB_   = $(LIB:%=$(LDIR)/%)
DEPEND = $(ODIR)/GNUmakefile.dep

.PHONY: all clean clean-all

all: $(DEPEND) $(LIB_)

clean:
	rm -rf $(ODIR) $(DEPEND)

clean-all: clean
	rm -rf $(LDIR)

# rule to build/link executable
#$(EXEC): $(OBJS)
#	$(CC) $(LIBF) -o $(EXEC) $(OBJS)

# rule to build .o files
$(OBJS):
	$(CC) -c -o $@ $< $(CFLAGS)

#rule to build .a files
$(LDIR)/$(LIB): $(OBJS)
	$(AR) rcs $@ $(OBJS)

# make sure $(ODIR) exists
$(ODIR):
	mkdir $(ODIR)
$(LDIR):
	mkdir $(LDIR)

# generate dependency makefile for all source files
$(DEPEND): $(ODIR) $(LDIR)
	$(CC) $(CFLAGS) -MM $(SRCS) $(HDRS) | sed -E 's,^[[:alpha:]],$(ODIR)/&,' > $(DEPEND)

-include $(DEPEND)
