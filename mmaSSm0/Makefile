
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

libmmaSSm0.so: \
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
	gcc -fPIC -shared -I../ mmaSSm0.c  $^ -o $@

%.o:../%.c
	gcc -c -fPIC $<
