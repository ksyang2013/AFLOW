CC = g++
LD = g++

.SUFFIXES: .o

CFLAGS = -c -O2 -s -Wall 
LFLAGS =

INCLUDE = -I ./

LIBS    =

SRC_XMLIB = apl2agr.cpp  xmaxis.cpp  xmgracefile.cpp  xmgraph.cpp  xmline.cpp  xmset.cpp
SRC_XMLIB_h = xmaxis.h  xmexcept.h  xmgracefile.h  xmgraph.h  xmlib.h  xmline.h  xmset.h

SRC = apl2agr.cpp 

OBJS_XMLIB = ${SRC_XMLIB:.cpp=.o}
OBJS = ${SRC:.cpp=.o}

all: xmlib apl2agr

.cpp.o:
	${CC} ${INCLUDE} ${CFLAGS} ${LIBS} -o $@ $<

xmlib:${OBJS_XMLIB}
	rm -f xmlib.a
	ar vr xmlib.a ${OBJS_XMLIB}
	#	touch apl2agr.cpp

apl2agr:${OBJS} $(SRC_XMLIB_h) $(SRC_XMLIB) Makefile
	${LD} ${LFLAGS} -o apl2agr ${OBJS} xmlib.a

clean:
	rm -f *.o xmlib.a apl2agr *~
realclean: clean
clear:
	rm *~
install:
	cp apl2agr /usr/local/bin
