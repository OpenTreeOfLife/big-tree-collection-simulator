#ifndef TREE_TEMPLATE_HPP
#define TREE_TEMPLATE_HPP

#include <iostream>
#include <cassert>
#include <cstring>
#include <iterator>
#include <list>
#include <set>
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
template<typename T, typename U> class Tree;
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
		void clear() {
			this->_label.clear();
			this->_n_added = 0;
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

class TaxonNameUniverse {
	public:
		TaxonNameUniverse();
		nnodes_t get_num_labels() {
			return this->_str2TaxonLabel.size();
		}
		TaxonLabel _register_new_name(const std::string & label);
		const TaxonLabel * find_name(const std::string & label) const {
			std::unordered_map<std::string, TaxonLabel>::const_iterator it;
			it = this->_str2TaxonLabel.find(label);
			return (it == this->_str2TaxonLabel.end() ? nullptr : &(it->second));
		}
	private:
		std::unordered_map<std::string, TaxonLabel> _str2TaxonLabel;

};

template<typename T, typename U>
TaxonLabel add_taxon_label(TaxonNameUniverse & taxa, const std::string & label, Node<T> * nd, Tree<T, U> * tree);

template<typename T>
class Node {
	public:
		Node<T>()
			:_label(std::string(), 0),
			_left_child(nullptr),
			_right_sib(nullptr),
			_parent(nullptr) {
		}

		bool is_polytomy() const {
			const Node<T> * c = this->_left_child;
			return (c and c->_right_sib and c->_right_sib->_right_sib);
		}
		std::vector<Node *> get_children() {
			std::vector<Node *> v;
			Node<T> * c = this->_left_child;
			while (c) {
				v.push_back(c);
				c = c->_right_sib;
			}
			return v;
		}
		std::vector<const Node *> get_children() const {
			std::vector<const  Node *> v;
			const Node<T> * c = this->_left_child;
			while (c) {
				v.push_back(c);
				c = c->_right_sib;
			}
			return v;
		}
		nnodes_t get_num_children() const {
			nnodes_t nc = 0;
			const Node<T> * c = this->_left_child;
			while (c) {
				++nc;
				c = c->_right_sib;
			}
			return nc;
		}
		/// introduces `nd` between this and its parent
		void bisect_edge_with_node(Node<T> * nd) {
			assert(nd);
			assert(this);
			assert(this->_parent);
			Node<T> * left_sib = this->find_left_sib();
			if (left_sib == nullptr) {
				this->_parent->_left_child = nd;
			}
			else {
				left_sib->_right_sib = nd;
			}
			nd->_parent = this->_parent;
			nd->_right_sib = this->_right_sib;
			this->_right_sib = nullptr;
			nd->_left_child = this;
			this->_parent = nd;
		}
		/// sets _left_child and _parent fields for a new left child. Does NOT alter previous _left_child Node
		void set_left_child(Node<T> * c) {
			this->_left_child = c;
			if (c) {
				c->_parent = this;
			}
		}
		bool is_leaf() const {
			return (this->_left_child == nullptr);
		}
		bool is_internal() const {
			return (this->_left_child != nullptr);
		}
		/// makes `r` the rsib of `this` and sets the parent of `r`
		void set_right_sib(Node<T> *r) {
			this->_right_sib = r;
			if (r != nullptr) {
				r->_parent = this->_parent;
			}
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
		typedef typename std::back_insert_iterator<std::list<const Node<T> *> > back_const_list_inserter;
		bool put_path_to_anc(const Node<T> * anc, back_const_list_inserter & b_it) const {
			Node<T> * curr_nd = this;
			for (;;) {
				b_it = curr_nd;
				if (curr_nd == anc) {
					return true;
				}
				curr_nd = curr_nd->get_parent();
				if (curr_nd == nullptr) {
					return false;
				}
			}
		}
		// if nullptr is returned then, *max_rank will be the rank of the root
		//	from `this` node.
		Node<T> * get_ancestor_by_edge_dist(unsigned anc_edge_dist, unsigned *max_rank) {
			unsigned curr_depth = 0;
			Node<T> * curr_nd = this;
			while (curr_depth < anc_edge_dist) {
				curr_nd = curr_nd->_parent;
				if (curr_nd == nullptr) {
					if (max_rank) {
						*max_rank = curr_depth;
					}
					return nullptr;
				}
				++curr_depth;
			}
			return curr_nd;
		}
		const Node<T> * find_rightmost_child() const {
			if (this->_left_child == nullptr) {
				return nullptr;
			}
			if (this->_left_child->_right_sib == nullptr) {
				return this->_left_child;
			}
			return this->_left_child->get_rightmost_sib();
		}
		Node<T> * find_left_sib(){
			if (this->_parent == nullptr) {
				return nullptr;
			}
			Node<T> * ls = this->_parent->_left_child;
			assert(ls);
			if (ls == this) {
				return nullptr;
			}
			assert(ls->_right_sib);
			while (ls->_right_sib != this) {
				assert(ls->_right_sib);
				ls = ls->_right_sib;
			}
			return ls;
		}
		const Node<T> * find_left_sib() const {
			if (this->_parent == nullptr) {
				return nullptr;
			}
			const Node<T> * ls = this->_parent->_left_child;
			assert(ls);
			if (ls == this) {
				return nullptr;
			}
			while (ls->_right_sib != this) {
				assert(ls->_right_sib);
				ls = ls->_right_sib;
			}
			return ls;
		}

		// Find the returns the node that is the leftmost child of the leftmost child of the leftmost child...
		Node<T> * find_furthest_left_des() {
			Node<T> * p = this->get_left_child();
			if (p == nullptr) {
				return nullptr;
			}
			Node<T> * n = p->get_left_child();
			while (n != nullptr) {
				 p = n;
				 n = p->get_left_child();
			}
			return p;
		}
		const Node<T> * find_furthest_left_des() const {
			const Node<T> * p = this->get_left_child();
			if (p == nullptr) {
				return nullptr;
			}
			const Node<T> * n = p->get_left_child();
			while (n != nullptr) {
				 p = n;
				 n = p->get_left_child();
			}
			return p;
		}

		// Find the returns the node that is the rightmost child of the rightmost child of the rightmost child...
		const Node<T> * find_furthest_right_des() const {
			const Node<T> * p = this->find_rightmost_child();
			if (p == nullptr) {
				return nullptr;
			}
			const Node<T> * n = p->find_rightmost_child();
			while (n != nullptr) {
				 p = n;
				 n = p->find_rightmost_child();
			}
			return p;
		}
		Node<T> * get_rightmost_sib();
		void set_label(const TaxonLabel & tl) {
			this->_label = tl;
		}
		TaxonLabel get_label() const {
			return this->_label;
		}

		// children before parents. Reverse of the preorder
		class const_child_iterator {
			public:
				const_child_iterator(const Node<T> * nd) {
					this->_curr_nd = (nd ? nd->_left_child : nullptr);
				}
				const_child_iterator(const const_child_iterator & other) {
					this->_curr_nd = other._curr_nd;
				}
				const Node<T> & operator*() {
					return *(this->_curr_nd);
				}
				const Node<T> * operator->() {
					return &(operator*());
				}
				const_child_iterator & operator++() {
					assert(this->_curr_nd);
					this->_curr_nd = this->_curr_nd->_right_sib;
					return *this;
				}
				const_child_iterator & operator++(int) {
					const_child_iterator tmp(*this);
					++(*this);
					return tmp;
				}
				bool operator==(const const_child_iterator & other) const {
					return (this->_curr_nd == other._curr_nd);
				}
				bool operator!=(const const_child_iterator & other) const {
					return (this->_curr_nd != other._curr_nd);
				}
			private:
				const Node<T> * _curr_nd;
		};
		// children before parents. Reverse of the preorder
		class child_iterator {
			public:
				child_iterator(Node<T> * nd) {
					this->_curr_nd = (nd ? nd->_left_child : nullptr);
				}
				child_iterator(const child_iterator & other) {
					this->_curr_nd = other._curr_nd;
				}
				Node<T> & operator*() {
					return (*(this->_curr_nd));
				}
				Node<T> * operator->() {
					return &(operator*());
				}
				child_iterator & operator++() {
					assert(this->_curr_nd);
					this->_curr_nd = this->_curr_nd->_right_sib;
					return *this;
				}
				child_iterator & operator++(int) {
					child_iterator tmp(*this);
					++(*this);
					return tmp;
				}
				bool operator==(const child_iterator & other) const {
					return (this->_curr_nd == other._curr_nd);
				}
				bool operator!=(const child_iterator & other) const {
					return (this->_curr_nd != other._curr_nd);
				}
			private:
				Node<T> * _curr_nd;
		};


		// children before parents. Reverse of the preorder
		class const_postorder_iterator {
			public:
				const_postorder_iterator(const Node<T> * nd) {
					this->_last_nd = nd;
					if (nd) {
						this->_curr_nd = nd->find_furthest_right_des();
						if (this->_curr_nd == nullptr) {
							this->_curr_nd = nd;
						}
					}
					else {
						this->_curr_nd = nullptr;
					}
				}
				const_postorder_iterator(const const_postorder_iterator & other) {
					this->_curr_nd = other._curr_nd;
					this->_last_nd = other._last_nd;
				}
				const Node<T> & operator*() {
					return *(this->_curr_nd);
				}
				const Node<T> * operator->() {
					return &(operator*());
				}
				const_postorder_iterator & operator++() {
					assert(this->_curr_nd);
					if (this->_curr_nd == this->_last_nd) {
						this->_curr_nd = nullptr;
						return *this;
					}
					const Node<T> * n = this->_curr_nd->find_left_sib();
					if (n == nullptr) {
						n = this->_curr_nd->get_parent();
						assert(n);
					}
					else if (n->is_internal()) {
						n = n->find_furthest_right_des();
					}
					this->_curr_nd = n;
					return *this;
				}
				const_postorder_iterator & operator++(int) {
					const_postorder_iterator tmp(*this);
					++(*this);
					return tmp;
				}
				bool operator==(const const_postorder_iterator & other) const {
					return (this->_curr_nd == other._curr_nd);
				}
				bool operator!=(const const_postorder_iterator & other) const {
					return (this->_curr_nd != other._curr_nd);
				}
			private:
				const Node<T> * _curr_nd;
				const Node<T> * _last_nd;
		};
		// children before parents. *Not* the reverse of the preorder (left first
		class const_fast_postorder_iterator {
			public:
				const_fast_postorder_iterator(const Node<T> * nd) {
					this->_last_nd = nd;
					if (nd) {
						this->_curr_nd = nd->find_furthest_left_des();
						if (this->_curr_nd == nullptr) {
							this->_curr_nd = nd;
						}
					}
					else {
						this->_curr_nd = nullptr;
					}
				}
				const_fast_postorder_iterator(const const_fast_postorder_iterator & other) {
					this->_curr_nd = other._curr_nd;
					this->_last_nd = other._last_nd;
				}
				const Node<T> & operator*() {
					return *(this->_curr_nd);
				}
				const Node<T> * operator->() {
					return &(operator*());
				}
				const_fast_postorder_iterator & operator++() {
					assert(this->_curr_nd);
					if (this->_curr_nd == this->_last_nd) {
						this->_curr_nd = nullptr;
						return *this;
					}
					const Node<T> * n = this->_curr_nd->get_right_sib();
					if (n == nullptr) {
						n = this->_curr_nd->get_parent();
						assert(n);
					}
					else if (n->is_internal()) {
						n = n->find_furthest_left_des();
					}
					this->_curr_nd = n;
					return *this;
				}
				const_fast_postorder_iterator & operator++(int) {
					const_fast_postorder_iterator tmp(*this);
					++(*this);
					return tmp;
				}
				bool operator==(const const_fast_postorder_iterator & other) const {
					return (this->_curr_nd == other._curr_nd);
				}
				bool operator!=(const const_fast_postorder_iterator & other) const {
					return (this->_curr_nd != other._curr_nd);
				}
			private:
				const Node<T> * _curr_nd;
				const Node<T> * _last_nd;
		};
		// children before parents. *Not* the reverse of the preorder (left first
		class fast_postorder_iterator {
			public:
				fast_postorder_iterator(Node<T> * nd) {
					this->_last_nd = nd;
					if (nd) {
						this->_curr_nd = nd->find_furthest_left_des();
						if (this->_curr_nd == nullptr) {
							this->_curr_nd = nd;
						}
					}
					else {
						this->_curr_nd = nullptr;
					}
				}
				fast_postorder_iterator(const fast_postorder_iterator & other) {
					this->_curr_nd = other._curr_nd;
					this->_last_nd = other._last_nd;
				}
				Node<T> & operator*() {
					return *(this->_curr_nd);
				}
				Node<T> * operator->() {
					return &(operator*());
				}
				fast_postorder_iterator & operator++() {
					assert(this->_curr_nd);
					if (this->_curr_nd == this->_last_nd) {
						this->_curr_nd = nullptr;
						return *this;
					}
					Node<T> * n = this->_curr_nd->get_right_sib();
					if (n == nullptr) {
						n = this->_curr_nd->get_parent();
						assert(n);
					}
					else if (n->is_internal()) {
						n = n->find_furthest_left_des();
					}
					this->_curr_nd = n;
					return *this;
				}
				fast_postorder_iterator & operator++(int) {
					fast_postorder_iterator tmp(*this);
					++(*this);
					return tmp;
				}
				bool operator==(const fast_postorder_iterator & other) const {
					return (this->_curr_nd == other._curr_nd);
				}
				bool operator!=(const fast_postorder_iterator & other) const {
					return (this->_curr_nd != other._curr_nd);
				}
			private:
				Node<T> * _curr_nd;
				Node<T> * _last_nd;
		};


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

		const_child_iterator begin_child() const {
			return const_child_iterator(this);
		}
		const_child_iterator end_child() const {
			return const_child_iterator(nullptr);
		}
		child_iterator begin_child() {
			return child_iterator(this);
		}
		child_iterator end_child() {
			return child_iterator(nullptr);
		}
		const_preorder_iterator begin_preorder() const {
			return const_preorder_iterator(this);
		}
		const_preorder_iterator end_preorder() const {
			return const_preorder_iterator(nullptr);
		}
		fast_postorder_iterator begin_fast_postorder() {
			return fast_postorder_iterator(this);
		}
		fast_postorder_iterator end_fast_postorder() {
			return fast_postorder_iterator(nullptr);
		}
		const_fast_postorder_iterator begin_fast_postorder() const {
			return const_fast_postorder_iterator(this);
		}
		const_fast_postorder_iterator end_fast_postorder() const {
			return const_fast_postorder_iterator(nullptr);
		}
		const_postorder_iterator begin_postorder() const {
			return const_postorder_iterator(this);
		}
		const_postorder_iterator end_postorder() const {
			return const_postorder_iterator(nullptr);
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
						o << ':' << nd.blob.get_edge_length();
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
								o << ':' << a->blob.get_edge_length();
							}
							c = a;
							a = a->get_parent();
						}
						if (a and c->get_right_sib() != nullptr) {
							o << ',';
						}
					}
				}
			}
		}

		void clear() {
			this->_label.clear();
			this->_left_child = nullptr;
			this->_right_sib = nullptr;
			this->_parent = nullptr;
		}
	void debug_check_subtree_nav_pointers() const {
#		if ! defined(NDEBUG)
			std::stack<const Node<T> *> c_stack;
			const Node * focal_nd = this;
			std::set<const Node *> checked;
			for (;;) {
				checked.insert(focal_nd);
				const Node * ls = focal_nd->find_left_sib();
				if (ls == nullptr) {
					assert(focal_nd->_parent == nullptr or focal_nd->_parent->_left_child == focal_nd);
				}
				else {
					assert(focal_nd == ls->_right_sib);
				}
				if (focal_nd->is_internal()) {
					typename std::vector<const Node<T> *> v = focal_nd->get_children();
					for (typename std::vector<const Node<T> *>::const_iterator i = v.begin(); i != v.end(); ++i) {
						if (checked.find(*i) != checked.end()) {
							c_stack.push(*i);
						}
						const Node<T> * c_nd = *i;
						assert(c_nd->_parent == focal_nd);
					}
					focal_nd = v[0];
				}
				else if (!c_stack.empty()) {
					focal_nd = c_stack.top();
					c_stack.pop();
				}
				else {
					break;
				}
			}
#		endif
	}

		void set_left_child_raw(Node<T> * n) {
			this->_left_child = n;
		}
		void set_parent_raw(Node<T> * n) {
			this->_parent = n;
		}
		void set_right_sib_raw(Node<T> * n) {
			this->_right_sib = n;
		}

		mutable T blob;

	private:
		TaxonLabel _label;
		Node<T> * _left_child;
		Node<T> * _right_sib;
		Node<T> * _parent;

		Node<T>(const Node<T> & other); // not defined, not copyable
		Node<T> & operator=(const Node<T> & other);// not defined, not copyable
};

template<typename T, typename U>
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
		void register_label2node(const TaxonLabel & tl, Node_T * nd) {
			this->_label2node[tl.get_label()] = nd;
		}
		Node_T * find_node_by_label(const TaxonLabel & tl) const {
			auto nd_it = this->_label2node.find(tl.get_label());
			if (nd_it == this->_label2node.end()) {
				return nullptr;
			}
			return nd_it->second;
		}

		typename Node_T::const_preorder_iterator begin_preorder() const {
			return typename Node_T::const_preorder_iterator(this->_root);
		}
		typename Node_T::const_preorder_iterator end_preorder() const {
			return typename Node_T::const_preorder_iterator(nullptr);
		}
		typename Node_T::fast_postorder_iterator begin_fast_postorder(){
			return typename Node_T::fast_postorder_iterator(this->_root);
		}
		typename Node_T::fast_postorder_iterator end_fast_postorder() {
			return typename Node_T::fast_postorder_iterator(nullptr);
		}
		typename Node_T::const_fast_postorder_iterator begin_fast_postorder() const {
			return typename Node_T::const_fast_postorder_iterator(this->_root);
		}
		typename Node_T::const_fast_postorder_iterator end_fast_postorder() const {
			return typename Node_T::const_fast_postorder_iterator(nullptr);
		}
		typename Node_T::const_postorder_iterator begin_postorder() const {
			return typename Node_T::const_postorder_iterator(this->_root);
		}
		typename Node_T::const_postorder_iterator end_postorder() const {
			return typename Node_T::const_postorder_iterator(nullptr);
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

		Node<T> * get_new_sib(Node<T> & old_nd);
		Node<T> * get_new_left_child(Node<T> & old_nd);

		void debug_check() const {
#			if ! defined(NDEBUG)
				this->_root->debug_check_subtree_nav_pointers();
#			endif
		}


		mutable U blob; //

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

		Tree<T, U>(const Tree<T, U> & other); // not defined, not copyable
		Tree<T, U> & operator=(const Tree<T, U> & other);// not defined, not copyable

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


template<typename T, typename U>
inline TaxonLabel add_taxon_label(TaxonNameUniverse & taxa, const std::string & label, Node<T> * nd, Tree<T, U> * tree) {
	TaxonLabel tl = taxa._register_new_name(label);
	//std::cerr << "Adding label \"" << tl.get_label() << "\"\n";
	nd->set_label(tl);
	if (tree) {
		tree->register_label2node(tl, nd);
	}
	return tl;
}

inline TaxonNameUniverse::TaxonNameUniverse()
	:_str2TaxonLabel(g_num_taxa_buckets){
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

template<typename T, typename U>
Node<T> * Tree<T, U>::get_new_sib(Node<T> & old_nd) {
	if (old_nd.get_parent() == nullptr) {
		throw std::range_error("root sib requested");
	}
	Node<T> * nd = this->get_new_node();
	nd->set_parent_raw(old_nd.get_parent());
	if (old_nd.get_right_sib()) {
		Node<T> * sib = old_nd.get_rightmost_sib();
		sib->set_right_sib_raw(nd);
	}
	else {
		old_nd.set_right_sib_raw(nd);
	}
	return nd;
}

template<typename T, typename U>
Node<T> * Tree<T, U>::get_new_left_child(Node<T> & old_nd) {
	if (old_nd.get_left_child() != nullptr) {
		throw std::range_error("new left child requested on internal");
	}
	Node<T> * nd = this->get_new_node();
	nd->set_parent_raw(&old_nd);
	old_nd.set_left_child_raw(nd);
	return nd;
}



template<typename T, typename U>
nnodes_t Tree<T, U>::get_initial_node_store_size() {
	return g_initial_node_store_size;
}

template<typename T, typename U>
void Tree<T, U>::set_initial_node_store_size(unsigned long x) {
	if (x < 1) {
		throw std::range_error("Tree::initial_node_store_size");
	}
   g_initial_node_store_size = x;
}

template<typename T, typename U>
Tree<T, U>::Tree()
	:_label2node(g_num_taxa_buckets) {
	this->_node_alloc_storage.push_back(std::vector<Node_T>(Tree::get_initial_node_store_size()));
	this->_curr_node_alloc_ptr = &(*this->_node_alloc_storage.rbegin());
	this->_curr_place_in_nap = 0L;
	this->_last_ind_in_nap = this->_curr_node_alloc_ptr->size() - 1;
	this->_root = this->get_new_node();
}


template<typename T, typename U>
nnodes_t Tree<T, U>::get_num_leaves() const {
	typename Node<T>::const_leaf_iterator s_it = this->begin_leaf();
	typename Node<T>::const_leaf_iterator e_it = this->end_leaf();
	nnodes_t num_nodes = 0;
	for (; s_it != e_it; ++s_it) {
		++num_nodes;
	}
	return num_nodes;
}
template<typename T, typename U>
void Tree<T, U>::_replinish_node_store() {
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

template<typename T, typename U>
Node<T> * Tree<T, U>::get_new_node() {
	if (this->_curr_place_in_nap == this->_last_ind_in_nap) {
		if (! this->_node_cache.empty()) {
			Node_T * nd = this->_node_cache.top();
			this->_node_cache.pop();
			nd->clear();
			return nd;
		}
		this->_replinish_node_store();
	}
	assert(this->_curr_place_in_nap < this->_last_ind_in_nap);
	Node_T & nd = (*this->_curr_node_alloc_ptr)[this->_curr_place_in_nap];
	this->_curr_place_in_nap++;
	return &nd;
}

std::vector<TaxonLabel> parse_labels_from_stream(std::istream & input,
												 TaxonNameUniverse & taxa);

template<typename T, typename U>
Tree<T, U> * parse_from_newick_stream(std::istream & input,
									  TaxonNameUniverse & taxa,
									  const T & nd_blob,
									  const U & tree_blob) {
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
	Tree<T, U> * tree = new Tree<T, U>();
	tree->blob = tree_blob;
	curr_node = tree->get_root();
	std::streambuf * rdbuf = input.rdbuf();
	for(;;) {
		signed char c = rdbuf->sbumpc(); filepos++; filecol++;
		if (c == ',') {
			if (curr_mode == IN_LABEL and prev_control_char != ':') {
				add_taxon_label(taxa, label, curr_node, tree);
			}
			else if (curr_node->get_label().empty() and curr_node->get_left_child() == nullptr) {
				throw ParseExcept("Expecting every leaf to have a label. Found a comma.", fileline, filecol, filepos);
			}
			curr_node = tree->get_new_sib(*curr_node);
			curr_node->blob = nd_blob;
			prev_control_char = c;
			curr_mode = OUT_OF_LABEL;
		}
		else if (c == ')') {
			if (curr_mode == IN_LABEL and prev_control_char != ':') {
				add_taxon_label(taxa, label, curr_node, tree);
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
			curr_node = tree->get_new_left_child(*curr_node);
			curr_node->blob = nd_blob;
			prev_control_char = c;
			curr_mode = OUT_OF_LABEL;
		}
		else if (c == ':') {
			if (curr_mode == IN_LABEL and prev_control_char != ':') {
				add_taxon_label(taxa, label, curr_node, tree);
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
					add_taxon_label(taxa, label, curr_node, tree); // handling this here works for typical trees, but it will also not reject internal node names before subtrees.
					cache.clear();
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

struct empty_node_blob_t {
	public:
		double get_edge_length() const {
			return 0.0;
		}
};
struct empty_tree_blob_t {
};

typedef Node<empty_node_blob_t> SlimNode;
typedef Tree<empty_node_blob_t, empty_tree_blob_t> SlimTree;

#endif // TREE_TEMPLATE_HPP
