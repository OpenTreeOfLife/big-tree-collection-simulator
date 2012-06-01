big-tree-sim
============
This tool will be capable of simulating a large of a collection of big trees
based on a classification. The motivation is to provide a set of trees that
could be useful for testing the efficiency of the data store for the open tree
of life project. 

Building
--------
If you are using the bootstrap-open-tree-software system from 
    https://github.com/OpenTreeOfLife/bootstrap-open-tree-software
then you'll have environmental variables in 

    cd $OPEN_TREE_DEPENDENCY_DIR 
    python download-dev-resource.py install ncl
    cd $OPEN_TREE_SOURCE_DIR/big-tree-collection-simulator

    sh bootstrap.sh 
    mkdir "build$OPEN_TREE_BUILD_TAG"
    cd "build$OPEN_TREE_BUILD_TAG"
    ../configure --prefix=$OPEN_TREE_INSTALL_DIR --with-ncl=$OPEN_TREE_INSTALL_DIR


