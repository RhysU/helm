# GNU-like toolchain assumed
ifeq (icc,${CC})
    HOWSTRICT ?= -std=c99 -Wall
else
    HOWSTRICT ?= -std=c99 -pedantic -Wall -Wextra
endif
HOWFAST ?= -g -O2 -DNDEBUG
CFLAGS  ?= $(HOWSTRICT) $(HOWFAST)

all:       helm.o
helm.o:    helm.c helm.h

clean:
	rm -f *.o

###################################################################
# Build Graphviz-based block diagram
# Requires versions after 14 October 2011
###################################################################
DOT ?= $(shell which dot)
ifneq "$(DOT)" ""

all:      helm.png helm.svg helm.eps
helm.png: helm.dot ; dot $M -Tpng      -o $@
helm.svg: helm.dot ; dot $M -Tsvg      -o $@
helm.eps: helm.dot ; dot $M -Tps:cairo -o $@

endif
