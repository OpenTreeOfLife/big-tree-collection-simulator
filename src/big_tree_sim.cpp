#include <iostream>
#include <fstream>
#include <cstring>
#include <string>
#include <cassert>
#include <list>
#include <stack>
#include <vector>
#include <stdexcept>


class Node;
class Tree;
class TaxonNameUniverse {
    public:
        unsigned long add_label(const char * label, Node * nd);
};
class Node {
    public:
        const char * get_label() const {
            return (this->_label.empty() ? nullptr : this->_label.c_str());
        }
        Node * get_left_child() const {
            return this->_left_child;
        }
        Node * get_new_sib(Tree *);
        Node * get_new_left_child(Tree *);
        Node * get_parent();
        Node()
            :_left_child(nullptr),
            _right_sib(nullptr),
            _parent(nullptr) {
        }
    private:
            
        std::string _label;
        Node * _left_child;
        Node * _right_sib;
        Node * _parent;
        
        friend class Tree;
};
class Tree {
    public:
        typedef unsigned long nnodes_t;
        static nnodes_t get_initial_node_store_size();
        static void set_initial_node_store_size(nnodes_t);

        Tree();
        Node * get_root() {
            return this->_root;
        }
        Node * get_new_node();
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
Tree::nnodes_t Tree::initial_node_store_size = 100;
Tree::nnodes_t Tree::get_initial_node_store_size() {
    return Tree::initial_node_store_size;
}
void Tree::set_initial_node_store_size(unsigned long x) {
    if (x < 1) {
        throw std::range_error("Tree::initial_node_store_size");
    }
    Tree::initial_node_store_size = x;
}
Tree::Tree() {
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

Tree * parse_from_newick_stream(std::istream & inp, TaxonNameUniverse & taxa) {
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
    while (inp.good()) {
        char c = inp.get(); filepos++; filecol++;
        if (std::isgraph(c)) {
            if (strchr("(),:;", c) != NULL) {
                if (c == ',') {
                    if (curr_mode == IN_LABEL) {
                        taxa.add_label(label.c_str(), curr_node);
                    }
                    else if (curr_node->get_label() == '\0' and curr_node->get_left_child() != nullptr) {
                        throw ParseExcept("Expecting every leaf to have a label. Found a comma.", fileline, filecol, filepos);
                    }
                    curr_node = curr_node->get_new_sib(tree);
                }
                else if (c == ')') {
                    if (curr_mode == IN_LABEL) {
                        taxa.add_label(label.c_str(), curr_node);
                    }
                    else if (curr_node->get_label() == '\0' and curr_node->get_left_child() != nullptr) {
                        throw ParseExcept("Expecting every leaf to have a label. Found a closed parenthesis.", fileline, filecol, filepos);
                    }
                    curr_node = curr_node->get_parent();
                }
                else if (c == '(') {
                    if (curr_mode == IN_LABEL) {
                        throw ParseExcept("Unexpected ( after taxon label. Expecting labels after clade", fileline, filecol, filepos);
                    }
                    curr_node = curr_node->get_new_left_child(tree);
                }
                else if (c == ':') {
                    if (curr_mode == IN_LABEL and prev_control_char != ':') {
                        taxa.add_label(label.c_str(), curr_node);
                    }
                }
                else if (c == ';') {
                    assert(curr_node == tree->get_root());
                    return tree;
                }
                prev_control_char = c;
                curr_mode = OUT_OF_LABEL;
            }
            else {
                if (curr_mode == OUT_OF_LABEL) {
                    label.clear();
                    if (c == '\'') {
                        if (!inp.good()) {
                            throw ParseExcept("Quote started then EOF", fileline, filecol, filepos);
                        }
                        for (;;) {
                            if (!inp.good()) {
                                throw ParseExcept("File reading before termination of quote", fileline, filecol, filepos);
                            }
                            char q = inp.get(); filepos++; filecol++;
                            if (q == '\'') {
                                const char d = inp.peek();
                                if (d == '\'') {
                                    label.append(1, q);
                                    q = inp.get(); filepos++; filecol++;
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
                if (strchr("\'(){}\"-]/\\,;_:=*`+<>", c) != NULL) {
                    if (c == '_') {
                        label.append(1, ' ');
                    }
                    else {
                        label.append(1, c);
                    }
                }
                else if (c == '[') {
                    if (preserve_comments) {
                        label.append(1, '[');
                    }
                    comment_level = 1;
                    if (!inp.good()) {
                        throw ParseExcept("Comment started then EOF", fileline, filecol, filepos);
                    }
                    while (comment_level > 0) {
                        if (!inp.good()) {
                            throw ParseExcept("File reading before termination of comment", fileline, filecol, filepos);
                        }
                        char q = inp.get(); filepos++; filecol++;
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
            //std::cerr << c << '\n';
        }
        else if (curr_mode == IN_LABEL) {
            if (c == '\n') {
                filecol = 0;
                fileline = 1;
            }
            cache.append(1, ' '); // converts all non-graphical characters to spaces!
        }
    }
    return tree;
}
int main(int argc, char *argv[]) {
    if (argc < 2) {
        std::cerr << "Expecting a file path to a newick formatted tree\n";
        return 1;
    }
    std::string tree_filename(argv[1]);
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
    }
    catch (ParseExcept & x) {
        std::cerr << "\nError:  " << x.message << "\nAt line = " << x.fileline << " at column = " << x.filecol << " at pos = " << x.filepos << "\n";
        return 3;
    }
    return 0;
}
