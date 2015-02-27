srcdir     :=$(shell pwd)

CFLAGS= -O
fixed-flags = -I$(srcdir) -I$(srcdir)/box
all-C-flags:= -ansi $(fixed-flags) $(CFLAGS)
non-ansi-flags :=  $(fixed-flags) $(CFLAGS)

CC = gcc # some compiler for ANSI/ISO C

# These settings should also be used in subdirectories:
export CC all-C-flags fixed-flags CFLAGS

.SUFFIXES:

%.o: %.c
	$(CC) -c $(CPPFLAGS) $(all-C-flags) $<

common_objects=lexer.o parser.o\
 non-ANSI.o bigint.o binmat.o creatop.o gettype.o getvalue.o\
 init.o learn.o main.o mem.o node.o onoff.o output.o poly.o sym.o

objects=$(common_objects) print.o getl.o
GAP_objects=$(common_objects) gapprint.o gapgetl.o

# Global organisation (phony targets)

.PHONY: install all script finish no_readline
.PHONY:	math_functions binding_functions

# The first target makes everything to get an operational LiE program
install: all script INFO.a

# To 'make all', we first descend into the subdirectories, and afterwards
# return to finish here.

all:
	$(MAKE) math_functions binding_functions
	$(MAKE) finish

finish: Lie.exe LEARN.ind INFO.ind # do not call 'make finish' directly

math_functions:
	$(MAKE) -C box all

binding_functions:
	$(MAKE) -C static


# Real dependencies 

# The file non-ANSI.c should be compiled without -ansi flag

non-ANSI.o: non-ANSI.c
	$(CC) -c $(CPPFLAGS) $(non-ansi-flags) $<


# These object files depend on the data types declared in those .h files:

gettype.o getvalue.o init.o main.o mem.o node.o sym.o: memtype.h nodetype.h


# The parser is generated by bison (BSD-yacc will NOT do) and in compilation
# requires inclusion of parseaux.c (derived from parseaux.w).

parser.c parser.h: parser.y
	bison -d --output-file=parser.c parser.y
parser.o: parser.c parser.h parseaux.c

lexer.o: parser.h

# Binding to the GNU readline library is achieved by -Dpreprocessor below

getl.o:	getl.c
	$(CC) -c $(CPPFLAGS) -Dpreprocessor $(all-C-flags) $<
gapgetl.o: getl.c
	$(CC) -c $(CPPFLAGS) $(all-C-flags) -o gapgetl.o $<

# Though date.c never changes, it should be recompiled for each modifiaction.
# Since Liegap is a separate executable, it gets its own version of the date.
# Since global recompilation should be issued by 'make all' rather than by
# 'make Lie.exe', we may assume here that the '.last_compiled' dates have just
# been set to the most recent one of object files in the respective
# subdirectories, whence taking that dummy file as dependency suffices.

date.o: date.c $(objects) box/.last_compiled static/.last_compiled
	$(CC) -ansi -c date.c
gapdate.o: date.c $(GAP_objects) box/.last_compiled static/.last_compiled
	$(CC) -ansi -c -o gapdate.o date.c

Lie.exe: date.o
	$(CC) -o Lie.exe $(objects) date.o static/*.o box/*.o -lreadline
	chmod g+w Lie.exe
Liegap.exe: gapdate.o
	$(CC) -o Liegap.exe $(GAP_objects) gapdate.o static/*.o box/*.o
	chmod g+w Liegap.exe

noreadline: math_functions binding_functions $(common_objects) print.o
	$(CC) -c $(CPPFLAGS) $(all-C-flags) getl.c
	$(MAKE) date.o
	$(CC) -o Lie.exe $(objects) date.o static/*.o box/*.o
	chmod g+w Lie.exe
	$(MAKE) LEARN.ind INFO.ind script INFO.a

script:
	./make_lie

INFO.ind:	INFO.0 INFO.1 INFO.2 INFO.3 INFO.4 infoind
	./infoind

LEARN.ind:	LEARN learnind
	./learnind

infoind: util/infoind.c
	$(MAKE) -C util ../infoind
learnind: util/learnind.c
	$(MAKE) -C util ../learnind

INFO.a: progs/maxsub progs/maxsub0 progs/eqrank
	rm -f INFO.a
	./Lie.exe < progs/maxsub

.PHONY: clean
clean:
	$(MAKE) -C box clean
	$(MAKE) -C static clean
	rm -f *~ *.o parser.[ch] INFO.a LEARN.ind
	rm -f Lie.exe Liegap.exe infoind learnind util/*.o