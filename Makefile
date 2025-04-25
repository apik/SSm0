CFLAGS= -O2 -Wall -pedantic -I.
LDFLAGS= -L. -lSSm0 -lm

all:libSSm0.a ex_SSm0

.INTERMEDIATE:\
	SSm0_SijLi22.o \
	SSm0_SijLiR.o \
	SSm0_SijLiC.o \
	SSm0_IijLiR.o \
	SSm0_IijLiC.o \
	SSm0_tld.o \
	SSm0_li2.o \
	SSm0_li3.o \
	SSm0_li4.o \
	SSm0_li22.o \
	SSm0_internal.o

libSSm0.a: \
	SSm0_SijLi22.o \
	SSm0_SijLiR.o \
	SSm0_SijLiC.o \
	SSm0_IijLiR.o \
	SSm0_IijLiC.o \
	SSm0_tld.o \
	SSm0_li2.o \
	SSm0_li3.o \
	SSm0_li4.o \
	SSm0_li22.o \
	SSm0_internal.o
	ar rcs $@ $^

ex_SSm0:ex_SSm0.c
	gcc -I. ex_SSm0.c $(LDFLAGS) -o ex_SSm0


.PHONY:check

check:tests/test_SSm0.c
	gcc -I. tests/test_SSm0.c $(LDFLAGS) -o tests/test_SSm0


%.o:%.c
	gcc -c $(CFLAGS) $<
