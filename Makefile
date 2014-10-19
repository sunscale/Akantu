# Configuration options.
cc         = not-set
cxx        = not-set
prefix     = not-set
shared     = not-set
assert     = not-set
verbose    = not-set
doc        = not-set
latex      = not-set
build      = not-set

# Basically proxies everything to the builddir cmake.

cputype = $(shell uname -m | sed "s/\\ /_/g")
systype = $(shell uname -s)

BUILDDIR = build/$(systype)-$(cputype)

# Process configuration options.
CONFIG_FLAGS = -DCMAKE_VERBOSE_MAKEFILE=0
ifneq ($(assert), not-set)
    CONFIG_FLAGS += -DNDEGUB=$(assert)
endif
ifneq ($(prefix), not-set)
    CONFIG_FLAGS += -DCMAKE_INSTALL_PREFIX=$(prefix)
endif
ifneq ($(shared), not-set)
    CONFIG_FLAGS += -DAKANTU_SHARED=$(shared)
endif
ifneq ($(verbose), not-set)
    CONFIG_FLAGS += -DAKANTU_VERBOSE=$(verbose)
endif
ifneq ($(cxx), not-set)
    CONFIG_FLAGS += -DCMAKE_CXX_COMPILER=$(cxx)
endif
ifneq ($(cc), not-set)
    CONFIG_FLAGS += -DCMAKE_C_COMPILER=$(cc)
endif
ifneq ($(doc), not-set)
    CONFIG_FLAGS += -DDOCUMENTATION=$(doc)
endif
ifneq ($(latex), not-set)
    CONFIG_FLAGS += -DLaTeX=$(latex)
endif
ifneq ($(build), not-set)
    CONFIG_FLAGS += -DCMAKE_BUILD_TYPE=$(build)
endif


define run-config
mkdir -p $(BUILDDIR)
cd $(BUILDDIR) && cmake $(CURDIR) $(CONFIG_FLAGS)
endef

all clean install doc check examples package package_source depend edit_cache install/local install/strip list_install_components rebuild_cache a.out:
	@if [ ! -f $(BUILDDIR)/Makefile ]; then \
		more README; \
	else \
	  	make -C $(BUILDDIR) $@ $(MAKEFLAGS); \
	fi

test:
	@if [ ! -f $(BUILDDIR)/Makefile ]; then \
		more README; \
	else \
	  	make -C $(BUILDDIR) test $(MAKEFLAGS); \
	fi

uninstall:
	xargs rm < $(BUILDDIR)/install_manifest.txt

config: distclean
	$(run-config)

distclean:
	rm -rf $(BUILDDIR)

remake:
	find . -name CMakeLists.txt -exec touch {} ';'


.PHONY: config distclean all clean install uninstall remake dist doc check examples package package_source depend edit_cache install/local install/strip list_install_components rebuild_cache a.out test