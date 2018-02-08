GNUCAP_CONF = $(shell which gnucap-conf$(SUFFIX))
PACKAGE_NAME = gnucap-sim

include Make.override

CXX = $(shell $(GNUCAP_CONF) --cxx)
GNUCAP_CPPFLAGS = $(shell $(GNUCAP_CONF) --cppflags) -DADD_VERSION -DPIC
GNUCAP_CXXFLAGS = $(shell $(GNUCAP_CONF) --cxxflags)
GNUCAP_LDFLAGS = $(shell $(GNUCAP_CONF) --ldflags)
GNUCAP_LIBDIR   = $(shell $(GNUCAP_CONF) --libdir)
GNUCAP_PKGLIBDIR = $(shell $(GNUCAP_CONF) --pkglibdir)
GNUCAP_EXEC_PREFIX = $(shell $(GNUCAP_CONF) --exec-prefix)
#	 GNUCAP_PREFIX = $(shell $(GNUCAP_CONF) --prefix)
# TODO complete gnucap-conf
GNUCAP_PREFIX   = $(shell $(GNUCAP_CONF) --exec-prefix)# BUG, should be prefix!
GNUCAP_PKGLIBDIR = $(GNUCAP_LIBDIR)/gnucap
GNUCAP_DOCDIR = $(GNUCAP_PREFIX)/share/doc

GNUCAP_CXXFLAGS+= -fPIC -shared
PLUGIN_FILES = s_sparam.so

s_sens_SOURCES = s_sens.cc
s_sparam_SOURCES = s_sparam.cc

CLEANFILES = ${EXEC_FILES} *.so *.o

all: $(EXEC_FILES) ${PLUGIN_FILES} ${CHECK_PLUGINS}

SPARAM_LIBS=-lgsl -lblas

%.so: %.cc
	$(CXX) $(CXXFLAGS) $(GNUCAP_CXXFLAGS) $(CPPFLAGS) $(GNUCAP_CPPFLAGS) -o $@ $< ${SPARAM_LIBS}

%.o: %.cc
	$(CXX) $(CXXFLAGS) $(GNUCAP_CXXFLAGS) $(CPPFLAGS) $(GNUCAP_CPPFLAGS) -o $@ -c $<

check: all
	${MAKE} $@ -C tests

d_sc2.so: $(s_sens_OBJS)
	$(CXX) -shared $+ -o $@ ${LIBS}

.PHONY: check

install: $(EXEC_FILES) ${PLUGIN_FILES}
	-install -d $(DESTDIR)/$(GNUCAP_PKGLIBDIR)
	install $(PLUGIN_FILES) $(DESTDIR)/$(GNUCAP_PKGLIBDIR)

	install -d $(DESTDIR)$(GNUCAP_DOCDIR)/$(PACKAGE_NAME)
	install README $(DESTDIR)$(GNUCAP_DOCDIR)/$(PACKAGE_NAME)

	install -d $(DESTDIR)$(GNUCAP_DOCDIR)/$(PACKAGE_NAME)/examples
	install $(EXAMPLES) $(DESTDIR)$(GNUCAP_DOCDIR)/$(PACKAGE_NAME)/examples

clean:
	rm -f $(CLEANFILES)

distclean: clean
	rm Make.override

Make.override:
	[ -e $@ ] || echo "# here you may override settings" > $@
