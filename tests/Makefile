GNUCAP=gnucap
GCENV=LD_LIBRARY_PATH=${LD_LIBRARY_PATH:%=%:}${CURDIR}/..

TESTS = \
	s_sparam.2.gc \
	s_sparam.3.gc

LOGS = ${TESTS:%=%.log}
OUTS = ${TESTS:%=%.out}

check: ${LOGS}

${OUTS}: %.out: %
	${GCENV} ${GNUCAP} < $< | sed -e '1,/^default\ plugins/d' > $@

${LOGS}: %.log: %.out %.ref
	diff $*.out $*.ref > $@

.PHONY: ${OUTS}
