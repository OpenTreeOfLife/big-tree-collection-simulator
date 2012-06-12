big-tree-collection-simulator
=============================
This tool will be capable of simulating a large of a collection of big trees
based on a classification. The motivation is to provide a set of trees that
could be useful for testing the efficiency of the data store for the open tree
of life project. A classification will be read in, substrees will be extracted
and topological changes will be introduced into them before they are written
out.

Building
--------
If you are using the bootstrap-open-tree-software system from 
    https://github.com/OpenTreeOfLife/bootstrap-open-tree-software
then you'll have environmental variables define such that you can install
the dependencies with:

    cd $OPEN_TREE_DEPENDENCY_DIR 
    python download-dev-resource.py install ncl
    cd $OPEN_TREE_SOURCE_DIR/big-tree-collection-simulator

Then from the top of this repository you can do something like this:

    sh bootstrap.sh 
    mkdir "build$OPEN_TREE_BUILD_TAG"
    cd "build$OPEN_TREE_BUILD_TAG"
    ../configure --prefix=$OPEN_TREE_INSTALL_DIR --with-ncl=$OPEN_TREE_INSTALL_DIR
    make

Note that "$OPEN_TREE_BUILD_TAG" is just a tag for the build (like "release" or
"debug" The <tt>--prefix</tt> argument should point to where you want to 
install the tool and <tt>--with-ncl</tt> should point to the prefix directory
specified when installing NCL.



