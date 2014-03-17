# GNU-like toolchain assumed
ifeq (icc,${CC})
    HOWSTRICT ?= -std=c99 -Wall
else
    HOWSTRICT ?= -std=c99 -pedantic -Wall -Wextra
endif
HOWFAST ?= -g -O2 -DNDEBUG
CFLAGS  ?= $(HOWSTRICT) $(HOWFAST)

all:       helm.o step3
helm.o:    helm.c helm.h
step3.o:   step3.c helm.h

clean:
	rm -f *.o step3

###################################################################
# Build Graphviz-based block diagram
# FIXME Require version after 14 October 2011
###################################################################
DOT ?= $(shell which dot)
ifneq "$(DOT)" ""

all:      helm.png helm.svg helm.eps
helm.png: helm.dot ; dot $M -Tpng      -o $@
helm.svg: helm.dot ; dot $M -Tsvg      -o $@
helm.eps: helm.dot ; dot $M -Tps:cairo -o $@

endif

###################################################################
# Build Doxygen
# FIXME Require version 1.8.3+ for full Markdown support
###################################################################
DOXYGEN ?= $(shell which doxygen)
ifneq "$(DOXYGEN)" ""

all:  docs
docs: Doxyfile README.md
docs: $(wildcard *.c)   $(wildcard *.h)
docs: $(wildcard *.png) $(wildcard *.eps)
	doxygen Doxyfile
clean: clean-docs
clean-docs:
	rm -rf html latex

endif

