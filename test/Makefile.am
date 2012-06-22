check:
	if test -z "$(OPEN_TREE_TAXONOMY_DIR)" ; \
		then echo "OPEN_TREE_TAXONOMY_DIR must defined for tests" ; \
		elif test -f "$(OPEN_TREE_TAXONOMY_DIR)/taxa_allCoL.tre" ; \
		then \
			time "$(top_builddir)/src/big-tree-sim" -n 1300000 "$(OPEN_TREE_TAXONOMY_DIR)/taxa_allCoL.tre" ; \
		else echo "$(OPEN_TREE_TAXONOMY_DIR)/taxa_allCoL.tre not found. Test skipped" ; \
	fi
