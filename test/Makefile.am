AM_CPPFLAGS = @CPPFLAGS@ "-I@abs_top_srcdir@/src"
AM_LDFLAGS  = @LDFLAGS@

bin_PROGRAMS = test-traversal

test_traversal_CXXFLAGS = -W -Wall -Wextra -pedantic -std=c++11

test_traversal_SOURCES = test_traversal.cpp \
	../src/tree_template.cpp

check:
	./test-traversal "@abs_top_srcdir@/test/small_complex.tre"
	if test -z "$(OPEN_TREE_TAXONOMY_DIR)" ; \
		then echo "OPEN_TREE_TAXONOMY_DIR must defined for tests" ; \
		elif test -f "$(OPEN_TREE_TAXONOMY_DIR)/taxa_allCoL.tre" ; \
		then \
			time echo | "$(top_builddir)/src/big-tree-sim" -x -n 1300000 "$(OPEN_TREE_TAXONOMY_DIR)/taxa_allCoL.tre" ; \
		else echo "$(OPEN_TREE_TAXONOMY_DIR)/taxa_allCoL.tre not found. Test skipped" ; \
	fi

