#include <iostream>
#include <ios>
#include <fstream>
#include <cstring>
#include <string>
#include <cassert>
#include <list>
#include <stack>
#include <vector>
#include <stdexcept>
#include <unordered_map>
#include <sstream>
#include <cstdlib>

typedef unsigned long nnodes_t; // big enough to index the # of nodes




class Node;
class Tree;
class TaxonLabel {
    public:
        TaxonLabel()
            :_label(),
            _n_added(0) {
        }
        TaxonLabel(const std::string & label, unsigned char len_disambig_suffix=0)
           :_label(label),
           _n_added(len_disambig_suffix) {
        }
        const char * c_str() const {
            return this->_label.c_str();
        }
        bool empty() const {
            return this->_label.empty();
        }
        const std::string & get_label() const {
            return this->_label;
        }
        unsigned char get_num_char_added() const {
            return this->_n_added;
        }
    private:
        std::string _label;
        unsigned char _n_added;
            
};
class TaxonNameUniverse {
    public:
        TaxonNameUniverse();
        TaxonLabel add_label(const std::string & label, Node * nd, Tree * tree);
        nnodes_t get_num_labels() {
            return this->_str2TaxonLabel.size();
        }
    private:
        std::unordered_map<std::string, TaxonLabel> _str2TaxonLabel;

        TaxonLabel _register_new_name(const std::string & label);
};
class Node {
    public:
        Node()
            :_label(std::string(), 0),
            _left_child(nullptr),
            _right_sib(nullptr),
            _parent(nullptr) {
        }
        Node * get_left_child() const {
            return this->_left_child;
        }
        Node * get_new_sib(Tree &);
        Node * get_new_left_child(Tree &);
        Node * get_rightmost_sib();
        Node * get_parent() {
            return this->_parent;
        }
        void set_label(const TaxonLabel & tl) {
            this->_label = tl;
        }
        TaxonLabel get_label() const {
            return this->_label;
        }
    private:
            
        TaxonLabel _label;
        Node * _left_child;
        Node * _right_sib;
        Node * _parent;
        
        friend class Tree;
};
class Tree {
    public:
        static nnodes_t get_initial_node_store_size();
        static void set_initial_node_store_size(nnodes_t);

        Tree();
        Node * get_root() {
            return this->_root;
        }
        Node * get_new_node();
        void register_label2node(const TaxonLabel & tl, Node *nd) {
            this->_label2node[tl.get_label()] = nd;
        }
    private:
        void _replinish_node_store();

        static unsigned long initial_node_store_size;

        Node * _root;
        // node storage
        std::stack<Node *> _node_cache;
        std::list<std::vector<Node> > _node_alloc_storage;
        std::vector<Node> * _curr_node_alloc_ptr;
        nnodes_t _curr_place_in_nap;
        nnodes_t _last_ind_in_nap;
        std::unordered_map<std::string, Node *> _label2node;
};

class ParseExcept {
    public:
        ParseExcept(const char * msg, long line, long col, long pos):
            message(msg),
            fileline(line),
            filecol(col),
            filepos(pos) {
            }
        std::string message;
        long fileline;
        long filecol;
        long filepos;
};

////////////////////////////////////////////////////////////////////////////////
// IMPL
unsigned long gNumBuckets = 100;
TaxonNameUniverse::TaxonNameUniverse()
    :_str2TaxonLabel(gNumBuckets){
}
inline TaxonLabel TaxonNameUniverse::add_label(const std::string & label, Node * nd, Tree * tree) {
    TaxonLabel tl = this->_register_new_name(label);
    //std::cerr << "Adding label \"" << tl.get_label() << "\"\n";
    nd->set_label(tl);
    if (tree) {
        tree->register_label2node(tl, nd);
    }
    return tl;
}


inline TaxonLabel TaxonNameUniverse::_register_new_name(const std::string & label) {
    std::unordered_map<std::string, TaxonLabel>::const_iterator labIt = this->_str2TaxonLabel.find(label); 
    if (labIt == this->_str2TaxonLabel.end()) {
        TaxonLabel tl(label, 0);
        this->_str2TaxonLabel[label] = tl;
        return tl;
    }
    std::string prefix = label;
    unsigned dup_num = 2;
    for (;;) {
        std::stringstream disambig_label;
        disambig_label << prefix << "_DUPLICATE" << dup_num++;
        std::string new_label;
        new_label.assign(disambig_label.str());
        labIt = this->_str2TaxonLabel.find(new_label); 
        if (labIt == this->_str2TaxonLabel.end()) {
            TaxonLabel tl(new_label, 0);
            this->_str2TaxonLabel[new_label] = tl;
            return tl;
        }
    }
}

inline Node * Node::get_rightmost_sib() {
    Node * sib = this->_right_sib;
    if (sib) {
        while (sib->_right_sib) {
            sib = sib->_right_sib;
        }
        return sib;
    }
    return nullptr;
}
Node * Node::get_new_sib(Tree & tree) {
    if (this->_parent == nullptr) {
        throw std::range_error("root sib requested");
    }
    Node * nd = tree.get_new_node();
    nd->_parent = this->_parent;
    if (this->_right_sib) {
        Node * sib = this->get_rightmost_sib();
        sib->_right_sib = nd;
    }
    else {
        this->_right_sib = nd;
    }
    return nd;
}
Node * Node::get_new_left_child(Tree & tree) {
    if (this->_left_child != nullptr) {
        throw std::range_error("new left child requested on internal");
    }
    Node * nd = tree.get_new_node();
    nd->_parent = this;
    this->_left_child = nd;
    return nd;
}


nnodes_t Tree::initial_node_store_size = 100;
nnodes_t Tree::get_initial_node_store_size() {
    return Tree::initial_node_store_size;
}
void Tree::set_initial_node_store_size(unsigned long x) {
    if (x < 1) {
        throw std::range_error("Tree::initial_node_store_size");
    }
    Tree::initial_node_store_size = x;
}
Tree::Tree()
    :_label2node(gNumBuckets) {
    this->_node_alloc_storage.push_back(std::vector<Node>(Tree::get_initial_node_store_size()));
    this->_curr_node_alloc_ptr = &(*this->_node_alloc_storage.rbegin());
    this->_curr_place_in_nap = 0L;
    this->_last_ind_in_nap = this->_curr_node_alloc_ptr->size() - 1;
    this->_root = this->get_new_node();
}

void Tree::_replinish_node_store() {
    assert(this->_curr_node_alloc_ptr);
    nnodes_t last_alloc_size = this->_curr_node_alloc_ptr->size();
    try {
        this->_node_alloc_storage.push_back(std::vector<Node>(2*last_alloc_size));
    }
    catch (...) {
        this->_node_alloc_storage.push_back(std::vector<Node>(last_alloc_size));
    }
    this->_curr_node_alloc_ptr = &(*this->_node_alloc_storage.rbegin());
    this->_curr_place_in_nap = 0L;
    this->_last_ind_in_nap = this->_curr_node_alloc_ptr->size() - 1;
}

Node * Tree::get_new_node() {
    if (this->_curr_place_in_nap == this->_last_ind_in_nap) {
        if (! this->_node_cache.empty()) {
            Node * nd = this->_node_cache.top();
            this->_node_cache.pop();
            return nd;
        }
        this->_replinish_node_store();
    }
    assert(this->_curr_place_in_nap < this->_last_ind_in_nap);
    Node & nd = (*this->_curr_node_alloc_ptr)[this->_curr_place_in_nap];
    this->_curr_place_in_nap++;
    return &nd;

}

Tree * parse_from_newick_stream(std::istream & input, TaxonNameUniverse & taxa) {
    enum READING_MODE {IN_LABEL, OUT_OF_LABEL, IN_QUOTE};
    READING_MODE curr_mode = OUT_OF_LABEL;
    std::string label;
    std::string cache;
    long filepos = 0;
    long filecol = 0;
    long fileline = 1;
    int comment_level = 0;
    const bool preserve_comments = false;
    char prev_control_char = '\0';
    Node * curr_node = nullptr;
    Tree * tree = new Tree();
    curr_node = tree->get_root();
    std::streambuf * rdbuf = input.rdbuf();
    for(;;) {
        signed char c = rdbuf->sbumpc(); filepos++; filecol++;
        if (c == ',') {
            if (curr_mode == IN_LABEL and prev_control_char != ':') {
                taxa.add_label(label, curr_node, tree);
            }
            else if (curr_node->get_label().empty() and curr_node->get_left_child() != nullptr) {
                throw ParseExcept("Expecting every leaf to have a label. Found a comma.", fileline, filecol, filepos);
            }
            curr_node = curr_node->get_new_sib(*tree);
            prev_control_char = c;
            curr_mode = OUT_OF_LABEL;
        }
        else if (c == ')') {
            if (curr_mode == IN_LABEL and prev_control_char != ':') {
                taxa.add_label(label, curr_node, tree);
            }
            else if (curr_node->get_label().empty() and curr_node->get_left_child() != nullptr) {
                throw ParseExcept("Expecting every leaf to have a label. Found a closed parenthesis.", fileline, filecol, filepos);
            }
            curr_node = curr_node->get_parent();
            prev_control_char = c;
            curr_mode = OUT_OF_LABEL;
        }
        else if (c == '(') {
            if (curr_mode == IN_LABEL) {
                throw ParseExcept("Unexpected ( after taxon label. Expecting labels after clade", fileline, filecol, filepos);
            }
            curr_node = curr_node->get_new_left_child(*tree);
            prev_control_char = c;
            curr_mode = OUT_OF_LABEL;
        }
        else if (c == ':') {
            if (curr_mode == IN_LABEL and prev_control_char != ':') {
                taxa.add_label(label, curr_node, tree);
            }
            prev_control_char = c;
            curr_mode = OUT_OF_LABEL;
        }
        else if (c == ';') {
            assert(curr_node == tree->get_root());
            return tree;
        }
        else if(isgraph(c)) {
            if (curr_mode == OUT_OF_LABEL) {
                label.clear();
                if (c == '\'') {
                    if (!input.good()) {
                        throw ParseExcept("Quote started then EOF", fileline, filecol, filepos);
                    }
                    for (;;) {
                        if (!input.good()) {
                            throw ParseExcept("File reading before termination of quote", fileline, filecol, filepos);
                        }
                        char q = rdbuf->sbumpc(); filepos++; filecol++;
                        if (q == '\'') {
                            const signed char d = rdbuf->sgetc(); // peek
                            if (d == '\'') {
                                label.append(1, q);
                                q = rdbuf->sbumpc(); filepos++; filecol++;
                            }
                            else {
                                break;
                            }
                        }
                        else {
                            label.append(1, q);
                        }
                    }
                }
                else {
                    curr_mode = IN_LABEL;
                    cache.clear();
                }
            }
            else if (!cache.empty()) {
                label += cache;
                cache.clear();
            }
            if (c == '_') {
                label.append(1, ' ');
            }
            else if (c == '[') {
                if (preserve_comments) {
                    label.append(1, '[');
                }
                comment_level = 1;
                if (!input.good()) {
                    throw ParseExcept("Comment started then EOF", fileline, filecol, filepos);
                }
                while (comment_level > 0) {
                    if (!input.good()) {
                        throw ParseExcept("File reading before termination of comment", fileline, filecol, filepos);
                    }
                    char q = rdbuf->sbumpc(); filepos++; filecol++;
                    if (preserve_comments) {
                        label.append(1, q);
                    }
                    if (q == '[')
                        comment_level++;
                    else if (q == ']')
                        comment_level--;
                }
            }
            else {
                label.append(1, c);
            }
        }
        else {
            if (c == EOF) {
                assert(curr_node == tree->get_root());
                return tree;
            }
            if (curr_mode == IN_LABEL) {
                if (c == '\n') {
                    filecol = 0;
                    fileline = 1;
                }
                cache.append(1, ' '); // converts all non-graphical characters to spaces!
            }
        }
    }
    return tree;
}
int main(int argc, char *argv[]) {
    std::ios_base::sync_with_stdio(false);
    if (argc < 2) {
        std::cerr << "Expecting a file path to a newick formatted tree\n";
        return 1;
    }
    std::string tree_filename;
    char prev_flag = '\0';
    char * endptr;
    for (int i = 1; i < argc; ++i) {
        if (prev_flag == 'n') {
            long int n = std::strtol(argv[i], &endptr, 10);
            if (endptr == argv[i]) {
                std::cerr << "Expecting number after -n\n";
                return 1;
            }
            if (n < 4) {
                std::cerr << "Expecting number of leaves > 4\n";
                return 1;
            }
            gNumBuckets = (unsigned) n;
            prev_flag = '\0';
        }
        else if (prev_flag == '\0') {
            std::string arg(argv[i]);
            if (arg.length() > 1 and arg[0] == '-') {
                if (arg.length() > 2) {
                    std::cerr << "Expecting each flag to be specified as a separate argument\n";
                    return 1;
                }
                prev_flag = arg[1];
            }
            else {
                if (!tree_filename.empty()) {
                    std::cerr << "Expecting only one tree filename argument\n";
                    return 1;
                }
                tree_filename.assign(argv[i]);
            }
        }
    }

    if (tree_filename.empty()) {
        std::cerr << "Expecting only a tree filename argument\n";
        return 1;
    }
    std::ifstream inp(tree_filename.c_str());
    if (!inp.good()) {
        std::cerr << "Could not open " << tree_filename << "\n";
        return 2;
    }
    std::ostream & out(std::cout);
    TaxonNameUniverse taxa;
    const Tree * tree = nullptr;
    try {
        tree = parse_from_newick_stream(inp, taxa);
        out << '\n';
        std::cerr << taxa.get_num_labels() << " labels read.\n";
    }
    catch (ParseExcept & x) {
        std::cerr << "\nError:  " << x.message << "\nAt line = " << x.fileline << " at column = " << x.filecol << " at pos = " << x.filepos << "\n";
        return 3;
    }
    return 0;
}
