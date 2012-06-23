#ifndef TREE_TEMPLATE_HPP
#define TREE_TEMPLATE_HPP

#include <iostream>
#include <cassert>
#include <cstring>
#include <list>
#include <sstream>
#include <stack>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

typedef unsigned long nnodes_t; // big enough to index the # of nodes


/// Setting these appropriately can improve efficiency.
extern nnodes_t g_num_taxa_buckets; // 3* numbere of taxa works well
extern nnodes_t g_initial_node_store_size;

template<typename T> class Node;
template<typename T> class Tree;
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
        std::string newick() const {
            std::string _newick;
            if (true) {
                _newick = this->_create_newick_escaped();
            }
            else {
                _newick = this->_label;
            }
            return _newick;
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
        std::string _create_newick_escaped() const {
            bool need_quotes = false;
            for (auto c : this->_label) {
                if (std::strchr("\'[(){}\"-]/\\,;:=*`+<> \t\n\r", c) != nullptr) {
                    need_quotes = true;
                    break;
                }
            }
            if (need_quotes) {
                std::string _newick;
                _newick.reserve(this->_label.length() + 2);
                _newick = "\'";
                for (auto c : this->_label) {
                    if (c == '\'') {
                        _newick.append(1, '\'');
                    }
                    _newick.append(1, c);
                }
                _newick.append(1, '\'');
                return _newick;
            }
            else {
                return this->_label;
            }
        }

        std::string _label;
        unsigned char _n_added;
            
};

template<typename T>
class TaxonNameUniverse {
    public:
        TaxonNameUniverse();
        TaxonLabel add_label(const std::string & label, Node<T> * nd, Tree<T> * tree);
        nnodes_t get_num_labels() {
            return this->_str2TaxonLabel.size();
        }
    private:
        std::unordered_map<std::string, TaxonLabel> _str2TaxonLabel;

        TaxonLabel _register_new_name(const std::string & label);
};

template<typename T>
class Node {
    public:
        Node<T>()
            :_label(std::string(), 0),
            _left_child(nullptr),
            _right_sib(nullptr),
            _parent(nullptr) {
        }
        bool is_leaf() const {
            return (this->_left_child == nullptr);
        }
        bool is_internal() const {
            return (this->_left_child != nullptr);
        }
        const Node<T> * get_right_sib() const {
            return this->_right_sib;
        }
        const Node<T> * get_left_child() const {
            return this->_left_child;
        }
        const Node<T> * get_parent() const {
            return this->_parent;
        }
        Node<T> * get_right_sib() {
            return this->_right_sib;
        }
        Node<T> * get_left_child() {
            return this->_left_child;
        }
        Node<T> * get_parent() {
            return this->_parent;
        }
        Node<T> * get_new_sib(Tree<T> &);
        Node<T> * get_new_left_child(Tree<T> &);
        Node<T> * get_rightmost_sib();
        void set_label(const TaxonLabel & tl) {
            this->_label = tl;
        }
        TaxonLabel get_label() const {
            return this->_label;
        }
        const T & get_const_blob() const {
            return this->blob;
        }

        class const_preorder_iterator {
            public:
                const_preorder_iterator(const Node<T> * nd) {
                    this->_curr_nd = nd;
                }
                const_preorder_iterator(const const_preorder_iterator & other) {
                    this->_curr_nd = other._curr_nd;
                    this->_sib_stack = other._sib_stack;
                }
                const Node<T> & operator*() {
                    return *(this->_curr_nd);
                }
                const Node<T> * operator->() {
                    return &(operator*());
                }
                const_preorder_iterator & operator++() {
                    assert(this->_curr_nd);
                    const Node<T> * lc = this->_curr_nd->get_left_child();
                    if (lc) {
                        const Node<T> * lcrs = lc->get_right_sib();
                        if (lcrs) {
                            this->_sib_stack.push(lcrs);
                        }
                        this->_curr_nd = lc;
                    }
                    else if ( this->_sib_stack.empty()) {
                        this->_curr_nd = nullptr;
                    }
                    else {
                        this->_curr_nd = this->_sib_stack.top();
                        this->_sib_stack.pop();
                        const Node<T> * ps = this->_curr_nd->get_right_sib();
                        if (ps) {
                            this->_sib_stack.push(ps);                            
                        }
                    }
                    return *this;
                }
                const_preorder_iterator & operator++(int) {
                    const_preorder_iterator tmp(*this);
                    ++(*this);
                    return tmp;
                }
                bool operator==(const const_preorder_iterator & other) const {
                    return (this->_curr_nd == other._curr_nd);
                }
                bool operator!=(const const_preorder_iterator & other) const {
                    return (this->_curr_nd != other._curr_nd);
                }
            private:
                const Node<T> * _curr_nd;
                std::stack<const Node<T> *> _sib_stack;
        };

        class const_leaf_iterator {
            public:
                const_leaf_iterator(const Node<T> * nd) {
                    this->_curr_nd = nd;
                    while (this->_curr_nd && this->_curr_nd->is_internal()) {
                        this->_advance();
                    }
                }
                const_leaf_iterator(const const_leaf_iterator & other) {
                    this->_curr_nd = other._curr_nd;
                    this->_sib_stack = other._sib_stack;
                    while (this->_curr_nd && this->_curr_nd->is_internal()) {
                        this->_advance();
                    }
                }
                const Node<T> & operator*() {
                    return *(this->_curr_nd);
                }
                const Node<T> * operator->() {
                    return &(operator*());
                }
                const_leaf_iterator & operator++() {
                    this->_advance();
                    while (this->_curr_nd && this->_curr_nd->is_internal()) {
                        this->_advance();
                    }
                    return *this;
                }
                const_leaf_iterator & operator++(int) {
                    const_leaf_iterator tmp(*this);
                    ++(*this);
                    return tmp;
                }
                bool operator==(const const_leaf_iterator & other) const {
                    return (this->_curr_nd == other._curr_nd);
                }
                bool operator!=(const const_leaf_iterator & other) const {
                    return (this->_curr_nd != other._curr_nd);
                }
            private:
                const_leaf_iterator & _advance() {
                    assert(this->_curr_nd);
                    const Node<T> * lc = this->_curr_nd->get_left_child();
                    if (lc) {
                        const Node<T> * lcrs = lc->get_right_sib();
                        if (lcrs) {
                            this->_sib_stack.push(lcrs);
                        }
                        this->_curr_nd = lc;
                    }
                    else if ( this->_sib_stack.empty()) {
                        this->_curr_nd = nullptr;
                    }
                    else {
                        this->_curr_nd = this->_sib_stack.top();
                        this->_sib_stack.pop();
                        const Node<T> * ps = this->_curr_nd->get_right_sib();
                        if (ps) {
                            this->_sib_stack.push(ps);                            
                        }
                    }
                    return *this;
                }

                const Node<T> * _curr_nd;
                std::stack<const Node<T> *> _sib_stack;
        };


        const_preorder_iterator begin_preorder() const {
            return const_preorder_iterator(this);
        }        
        const_preorder_iterator end_preorder() const {
            return const_preorder_iterator(nullptr);
        }        
        const_leaf_iterator begin_leaf() const {
            return const_leaf_iterator(this);
        }        
        const_leaf_iterator end_leaf() const {
            return const_leaf_iterator(nullptr);
        }        

        void write_newick(std::ostream &o, bool edge_lengths) const {
            const_preorder_iterator it = this->begin_preorder();
            const_preorder_iterator end_it = this->end_preorder();
            for (; it != end_it; ++it) {
                //std::cerr << "nd\n";
                const Node<T> & nd = *it;
                if (nd.is_internal()) {
                    o << '(';
                }
                else {
                    o << nd.get_label().newick();
                    if (edge_lengths) {
                        o << ':' << nd.get_const_blob().edge_length();
                    }
                    if (nd.get_right_sib() != nullptr) {
                        o << ',';
                    }
                    else {
                        const Node<T> * c = &nd;
                        const Node<T> * a = nd.get_parent();
                        while (a and c != this and c->get_right_sib() == nullptr) {
                            o << ')' << a->get_label().newick();
                            if (edge_lengths) {
                                o << ':' << a->get_const_blob().edge_length();
                            }
                            c = a;
                            a = a->get_parent();
                        }
                        if (a and a != this) {
                            o << ',';
                        }
                    }
                }
            }
        }
    private:
        T blob;
        TaxonLabel _label;
        Node<T> * _left_child;
        Node<T> * _right_sib;
        Node<T> * _parent;
        
        friend class Tree<T>;
};

template<typename T>
class Tree {
    public:
        static nnodes_t get_initial_node_store_size();
        static void set_initial_node_store_size(nnodes_t);

        typedef Node<T> Node_T;

        Tree();
        Node_T * get_root() {
            return this->_root;
        }
        Node_T * get_new_node();
        void register_label2node(const TaxonLabel & tl, Node_T *nd) {
            this->_label2node[tl.get_label()] = nd;
        }

        typename Node_T::const_preorder_iterator begin_preorder() const {
            return Node_T::const_preorder_iterator(this->_root);
        }        
        typename Node_T::const_preorder_iterator end_preorder() const {
            return Node_T::const_preorder_iterator(nullptr);
        }        

        typename Node_T::const_leaf_iterator begin_leaf() const {
            return typename Node_T::const_leaf_iterator(this->_root);
        }        
        typename Node_T::const_leaf_iterator end_leaf() const {
            return typename Node_T::const_leaf_iterator(nullptr);
        }        
        nnodes_t get_num_leaves() const;
        void write_newick(std::ostream &o, bool edge_lengths) const {
            if (this->_root == nullptr) {
                return;
            }
            this->_root->write_newick(o, edge_lengths);
            o << ";";
        }
    private:
        void _replinish_node_store();

        Node_T * _root;
        // node storage
        std::stack<Node_T *> _node_cache;
        std::list<std::vector<Node_T > > _node_alloc_storage;
        std::vector<Node_T > * _curr_node_alloc_ptr;
        nnodes_t _curr_place_in_nap;
        nnodes_t _last_ind_in_nap;
        std::unordered_map<std::string, Node_T *> _label2node;
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

template<typename T> 
TaxonNameUniverse<T>::TaxonNameUniverse()
    :_str2TaxonLabel(g_num_taxa_buckets){
}

template<typename T> 
inline TaxonLabel TaxonNameUniverse<T>::add_label(const std::string & label, Node<T> * nd, Tree<T> * tree) {
    TaxonLabel tl = this->_register_new_name(label);
    //std::cerr << "Adding label \"" << tl.get_label() << "\"\n";
    nd->set_label(tl);
    if (tree) {
        tree->register_label2node(tl, nd);
    }
    return tl;
}


template<typename T> 
inline TaxonLabel TaxonNameUniverse<T>::_register_new_name(const std::string & label) {
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

template<typename T> 
inline Node<T> * Node<T>::get_rightmost_sib() {
    Node<T> * sib = this->_right_sib;
    if (sib) {
        while (sib->_right_sib) {
            sib = sib->_right_sib;
        }
        return sib;
    }
    return nullptr;
}

template<typename T>
Node<T> * Node<T>::get_new_sib(Tree<T> & tree) {
    if (this->_parent == nullptr) {
        throw std::range_error("root sib requested");
    }
    Node<T> * nd = tree.get_new_node();
    nd->_parent = this->_parent;
    if (this->_right_sib) {
        Node<T> * sib = this->get_rightmost_sib();
        sib->_right_sib = nd;
    }
    else {
        this->_right_sib = nd;
    }
    return nd;
}

template<typename T>
Node<T> * Node<T>::get_new_left_child(Tree<T> & tree) {
    if (this->_left_child != nullptr) {
        throw std::range_error("new left child requested on internal");
    }
    Node<T> * nd = tree.get_new_node();
    nd->_parent = this;
    this->_left_child = nd;
    return nd;
}



template<typename T>
nnodes_t Tree<T>::get_initial_node_store_size() {
    return g_initial_node_store_size;
}

template<typename T>
void Tree<T>::set_initial_node_store_size(unsigned long x) {
    if (x < 1) {
        throw std::range_error("Tree::initial_node_store_size");
    }
   g_initial_node_store_size = x;
}

template<typename T>
Tree<T>::Tree()
    :_label2node(g_num_taxa_buckets) {
    this->_node_alloc_storage.push_back(std::vector<Node_T>(Tree::get_initial_node_store_size()));
    this->_curr_node_alloc_ptr = &(*this->_node_alloc_storage.rbegin());
    this->_curr_place_in_nap = 0L;
    this->_last_ind_in_nap = this->_curr_node_alloc_ptr->size() - 1;
    this->_root = this->get_new_node();
}


template<typename T>
nnodes_t Tree<T>::get_num_leaves() const  {
    typename Node<T>::const_leaf_iterator s_it = this->begin_leaf();
    typename Node<T>::const_leaf_iterator e_it = this->end_leaf();
    nnodes_t num_nodes = 0;
    for (; s_it != e_it; ++s_it) {
        ++num_nodes;
    }
    return num_nodes;
}
template<typename T>
void Tree<T>::_replinish_node_store() {
    assert(this->_curr_node_alloc_ptr);
    nnodes_t last_alloc_size = this->_curr_node_alloc_ptr->size();
    try {
        this->_node_alloc_storage.push_back(std::vector<Node_T>(2*last_alloc_size));
    }
    catch (...) {
        this->_node_alloc_storage.push_back(std::vector<Node_T>(last_alloc_size));
    }
    this->_curr_node_alloc_ptr = &(*this->_node_alloc_storage.rbegin());
    this->_curr_place_in_nap = 0L;
    this->_last_ind_in_nap = this->_curr_node_alloc_ptr->size() - 1;
}

template<typename T>
Node<T> * Tree<T>::get_new_node() {
    if (this->_curr_place_in_nap == this->_last_ind_in_nap) {
        if (! this->_node_cache.empty()) {
            Node_T * nd = this->_node_cache.top();
            this->_node_cache.pop();
            return nd;
        }
        this->_replinish_node_store();
    }
    assert(this->_curr_place_in_nap < this->_last_ind_in_nap);
    Node_T & nd = (*this->_curr_node_alloc_ptr)[this->_curr_place_in_nap];
    this->_curr_place_in_nap++;
    return &nd;
}

template<typename T>
Tree<T> * parse_from_newick_stream(std::istream & input, TaxonNameUniverse<T> & taxa) {
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
    Node<T> * curr_node = nullptr;
    Tree<T> * tree = new Tree<T>();
    curr_node = tree->get_root();
    std::streambuf * rdbuf = input.rdbuf();
    for(;;) {
        signed char c = rdbuf->sbumpc(); filepos++; filecol++;
        if (c == ',') {
            if (curr_mode == IN_LABEL and prev_control_char != ':') {
                taxa.add_label(label, curr_node, tree);
            }
            else if (curr_node->get_label().empty() and curr_node->get_left_child() == nullptr) {
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
            else if (curr_node->get_label().empty() and curr_node->get_left_child() == nullptr) {
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

#endif // TREE_TEMPLATE_HPP
