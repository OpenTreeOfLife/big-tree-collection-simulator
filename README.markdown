# big-tree-collection-simulator

This tool will be capable of simulating a large of a collection of big trees
based on a classification. The motivation is to provide a set of trees that
could be useful for testing the efficiency of the data store for the open tree
of life project. A classification will be read in, substrees will be extracted
and topological changes will be introduced into them before they are written
out.

## Building

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

## Usage

To use, specify a tree to serve as the base and then provde commands. For 
example:

    big-tree-sim -c test/basic-commands.txt test/poly.tre
    
or, equivalently:

    echo 'resolve ; sample ; spr ; print ; quit ;' | big-tree-sim test/poly.tre

The structure of the simulation algorithm is:
    1. read in a newick tree with all of the names (from unnamed command line arg).
    2. resolve the tree (recommended, but not required).
    3. sample a set of leaves from the tree.
    4. apply SPR operations to the subsampled tree if you want a set of trees
        that may conflict with the input taxonomy.
    5. print the subsampled tree.
    
The program works with two "slots" for a tree:
    1. the full tree, and
    2. the current, focal tree.
Initally the focal tree points to the full tree. The SAMPLE command creates a
subsampled tree as the focal tree. The PRINT, RESOLVE, and SPR commands act on
the  current focal tree (which will be the full tree if no SAMPLE command has been 
used).

### SAMPLE command

When using SAMPLE, the leaves are chosen in proportion to their leaf weights. 
This allows you to specify a set of leaves that will be preferentially sampled
to mimic sampling biases (e.g. frequent occurrence of model organisms in
phylogenetic studies, for example). By default all weights are 1.0, but with 
the WEIGHT command you can supply a weight and a filename all leaves with labels
found in the file will be assigned the specified weight. Note that using weight
differences so large that the result in rounding error when added can result 
in the SAMPLE command failing!
 
Note that the leaf weights apply to every part of the SAMPLE procedure, but not
other commands (e.g. highly weighted leaves are not more likely to be selected
in the SPR command).

With the SAMPLE command, you can specify a root depth, ingroup and ougroup size.
The outgroup is always a set of leaves selected from the sister of your ingroup.
The basic sample procedure is to repeat the following until a tree satisfies
the constraint:
    1. randomly choose a leaf (according to weights) from the full tree,
    2. randomly choose an ingroup root depth within the constraints,
    3. See if the ingroup and outgroup sizes can be met by the chosen root, if
        "yes" go to step 4. If "no" go to step 2 (or back to step 1 if there 
        are no other rooting points to try for this leaf choice).
    4. Randomly choose an ingroup size according within the constraints.
    5. Randomly select other ingroup leaves (the leaf chosen in step 1 is one 
        member of the ingroup) or go back to 4 if there are not enough leaves
        in the clade.
    6. Repeat steps 4 and 5 for the outgroup.
    7. Make a copy of the tree induced by the select leaves and store it as the
        focal tree.

### SPR command

The SPR command randomly chooses a (non-root) node from the sampled tree and 
moves it a number of edges that is that is affected by the user constraints. 
Specifically a desired ReconLimit is selected as a uniform from the 
[ReconMin, ReconMax] range that the user controls. The node selected is move 
the ReconLimit number of edges if the tree is big enough to accommodate that 
many movements. Note that if the tree is unbalanced and the clade of all of the 
leaves but one is selected, then the SPR will have no affect (there is no where
to move that clade).  So the operation can be a no-op.

If the node selected by an SPR is attached to a polytomy, it will be reattached
at a node (creating or adding to a polytomy). Thus the SPR operation does not
change the number of nodes in the tree.

### PRINT command

The prints a newick (or NEXUS) tree to the output stream specified by the 
last invocation of the OUT command (standard output is the default).

## Caveats

See the Issues.txt file for a gotcha with respect to parsing of newick trees and
disambiguation of repeated labels.

