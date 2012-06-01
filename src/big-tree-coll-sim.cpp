#include "ncl/ncl.h"
#include <iostream>
#include <climits>
#include <cassert>
#include "ncl/nxsmultiformat.h"

static long gStrictLevel = 1;
static long gTreeNum = -1;

void filepathToPatristic(const char * filename, std::ostream * os);
bool processClassification(const NxsFullTreeDescription &treeDesc, std::ostream * os, std::ostream * errStr, NxsTaxaBlockAPI &taxaB, const unsigned treeN);

bool processClassification(const NxsFullTreeDescription &treeDesc, std::ostream * os, std::ostream * errStr, NxsTaxaBlockAPI &taxaB, const unsigned treeN) {
    if (os) {
        *os << "Hi.  A tree was read\n";
    }
}


void filepathToPatristic(const char * filename, ostream *os)
{
	assert(filename);
	BlockReaderList blocks;
	try{
		MultiFormatReader nexusReader(-1, NxsReader::WARNINGS_TO_STDERR);
		if (gStrictLevel != 2)
			nexusReader.SetWarningToErrorThreshold((int)NxsReader::FATAL_WARNING + 1 - (int) gStrictLevel);
		NxsTreesBlock * treesB = nexusReader.GetTreesBlockTemplate();
		assert(treesB);
		treesB->SetAllowImplicitNames(true);
		if (gStrictLevel < 2)
			{
			NxsStoreTokensBlockReader *storerB =  nexusReader.GetUnknownBlockTemplate();
			assert(storerB);
			storerB->SetTolerateEOFInBlock(true);
			}
		cerr << "Executing" <<endl;
		nexusReader.ReadFilepath(filename, MultiFormatReader::RELAXED_PHYLIP_TREE_FORMAT);
        blocks = nexusReader.GetUsedBlocksInOrder();
        int prevTrees = 0;
        for (BlockReaderList::const_iterator bIt = blocks.begin(); bIt != blocks.end(); ++bIt) {
            NxsBlock * b = *bIt;
            if (b && b->GetID() == "TREES") {
                //*os << "TREES block found" << std::endl;

                NxsTreesBlock * treesBPtr = (NxsTreesBlock *) b;
                NxsTaxaBlockAPI * taxaBPtr =  treesBPtr->GetTaxaBlockPtr(NULL);
                if (!taxaBPtr)
                    throw NxsException("Trees block is not connected to a taxa block -- I don\'t know how that happened");
                const int nTreesThisBlock = (int) treesBPtr->GetNumTrees();
                if (gTreeNum < 0) {
                    for (unsigned i = 0; i < (unsigned) nTreesThisBlock; ++i) {
                        const NxsFullTreeDescription & treeDesc = treesBPtr->GetFullTreeDescription(i);
                        processClassification(treeDesc, os, &(std::cerr), *taxaBPtr, i + prevTrees);
                        *os << std::endl;
                    }
                }
                else if (prevTrees + nTreesThisBlock > gTreeNum) {
                    const NxsFullTreeDescription & treeDesc = treesBPtr->GetFullTreeDescription(gTreeNum - prevTrees);
                    processClassification(treeDesc, os, &(std::cerr), *taxaBPtr, gTreeNum);
                    break;
                }
            }
        }
		for (BlockReaderList::const_iterator bIt = blocks.begin(); bIt != blocks.end(); ++bIt) {
			NxsBlock * b = *bIt;
			if (b)
				delete b;
		}
	}
	catch (const NxsException &x) {
		cerr << "Error:\n " << x.msg << endl;
		if (x.line >=0)
			cerr << "at line " << x.line << ", column (approximately) " << x.col << " (and file position "<< x.pos << ")" << endl;
		exit(2);
	}
}

void readFilepathAsNEXUS(const char *filename) {
	cerr << "[Reading " << filename << "	 ]" << endl;
	try {
		ostream * outStream = &std::cout;
		filepathToPatristic(filename, outStream);
	}
	catch (...) {
		cerr << "Normalizing of " << filename << " failed (with an exception)" << endl;
		exit(1);
	}
}

void readFilesListedIsFile(const char *masterFilepath)
{
	ifstream masterStream(masterFilepath);
	if (masterStream.bad()) {
		cerr << "Could not open " << masterFilepath << "." << endl;
		exit(3);
	}
	char filename[1024];
	while ((!masterStream.eof())  && masterStream.good()) {
		masterStream.getline(filename, 1024);
		if (strlen(filename) > 0 && filename[0] != '#')
			readFilepathAsNEXUS(filename);
	}
}

void printHelp(ostream & out)
{
	out << "big-tree-collection-simulator reads a NEXUS file with a classification\n";
	out << "\nThe most common usage is simply:\n    big-tree-collection-simulator <path to NEXUS file>\n";
	out << "\nCommand-line flags:\n\n";
	out << "    -h on the command line shows this help message\n\n";
}

int main(int argc, char *argv[])
{
	//sleep(10);
	if (argc < 2) {
		cerr << "Expecting the path to NEXUS file as the only command line argument!\n" << endl;
		printHelp(cerr);
		return 1;
	}
	for (int i = 1; i < argc; ++i) {
		const char * filepath = argv[i];

		if (strlen(filepath) > 1 && filepath[0] == '-' && filepath[1] == 'h')
			printHelp(cout);
		else if (strlen(filepath) > 2 && filepath[0] == '-' && filepath[1] == 't') {
			if (!NxsString::to_long(filepath + 2, &gTreeNum) || gTreeNum < 1) {
				cerr << "Expecting a positive integer after -t\n" << endl;
				printHelp(cerr);
				return 2;
			}
			--gTreeNum;
		}
			readFilepathAsNEXUS(filepath);
	}
	return 0;
}

