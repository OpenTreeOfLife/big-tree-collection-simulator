

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <functional>
#include <fstream>
#include <ios>
#include <random>

#include "config.h"
#include "tree_template.hpp"
////////////////////////////////////////////////////////////////////////////////
// RNGImpl
////////////////////

class RandGen {
	public:
		typedef unsigned long uint_seed_t;
		typedef double rng_float_t;
		RandGen(uint_seed_t s)
			:_u_dist(0,1.0) {
			if (s == 0 || s == std::numeric_limits<uint_seed_t>::max()) {
				s = (uint_seed_t) time(nullptr);
			}
			this->seed(s);
		}
		void seed(uint_seed_t x) {
			this->_engine.seed(x);
			//this->_generator = std::bind(this->_u_dist, this->_engine);
		}
		rng_float_t uniform01(void) {
			const rng_float_t x = this->_u_dist(this->_engine); //this->_generator();
#			if ! defined(NDEBUG)
#				if defined(DEBUGGING_RNG)
					std::cerr << "UNIFORMO1 " << x << "\n";
#				endif
#			endif
			return x;
		}
		rng_float_t operator()(void) {
			return this->uniform01();
		}

		//// Returns a random integer from the range [beg_ind, end_ind)
		// @ not the greatest algorithm...
		uint_seed_t rand_range(uint_seed_t beg_ind, uint_seed_t end_ind) {
			assert(end_ind > beg_ind);
			const long diff = end_ind - beg_ind;
			const double real_diff = float(diff);
			const double rand_offset = real_diff*this->uniform01();
			const long int_rand_offset(rand_offset);
			uint_seed_t ret = (beg_ind + int_rand_offset);
			if (ret >= end_ind) {
				assert(ret == end_ind);
				return ret - 1;
			}
			return ret;
		}
	private:
		std::uniform_real_distribution<rng_float_t> _u_dist;
		std::minstd_rand _engine;
		//std::function<double ()> _generator;
};

template<typename T>
T & choose_element(T * arr, RandGen::uint_seed_t sz, RandGen & rng) {
	assert(sz > 0);
	RandGen::uint_seed_t ind = rng.rand_range(0, sz);
	assert(ind < sz);
	return arr[ind];
}

template<typename T>
T & choose_element(std::vector<T> & vec, RandGen & rng) {
	assert(!vec.empty());
	RandGen::uint_seed_t ind = rng.rand_range(0, vec.size());
	assert(ind < vec.size());
#	if defined(DEBUGGING_RNG)
		std::cerr << "Element " << ind << " chosen.\n";
#	endif
	return vec.at(ind);
}


////////////////////////////////////////////////////////////////////////////////
// Tree manipulation impl
////////////////////


/// Resolves each polytomy in the tree with a stepwise addition type procedure.
// For every polytomy, the specific Node instance of the polytomy will correspond
//	to the ancestor of the resolved portion.
// Returns then number nodes added.
template<typename T, typename U>
nnodes_t resolve_polytomies(Tree<T, U> & tree, RandGen & rng) {
	typename Node<T>::fast_postorder_iterator it = tree.begin_fast_postorder();
	nnodes_t num_nodes_added = 0;
	std::vector<Node<T> *> attachment_points;

#	if ! defined(NDEBUG)
		tree.debug_check();
#	endif
	for (; it != tree.end_fast_postorder(); ++it) {
		Node<T> & nd = *it;
		if (nd.is_polytomy()) {
#			if ! defined(NDEBUG)
				tree.debug_check();
#			endif
			std::vector<Node<T> *> children = nd.get_children();
			assert(children.size() > 2);
			//std::cerr << " resolving polytomy with " <<	 children.size() << " children\n";
			children[1]->set_right_sib(nullptr);
			attachment_points.clear();
			attachment_points.reserve(children.size());
			attachment_points.push_back(&nd);
			attachment_points.push_back(children[0]);
			attachment_points.push_back(children[1]);

			for (nnodes_t ch_ind = 2; ch_ind < children.size(); ++ch_ind) {
				Node<T> * child = children[ch_ind];
				child->set_right_sib(nullptr);
				Node<T> * new_sib = choose_element(attachment_points, rng);
				Node<T> * new_connector = tree.get_new_node();
				++num_nodes_added;
				child->set_right_sib(nullptr);
				if (new_sib == &nd) {
					// special handling so that nd remains the deepest node in the resolved polytomy
					Node<T> * left_of_base = nd.get_left_child();
					assert(left_of_base);
					Node<T> * right_of_base = left_of_base->get_right_sib();
					assert(right_of_base);
					assert(right_of_base->get_right_sib() == nullptr);
					left_of_base->set_right_sib(nullptr);
					left_of_base->bisect_edge_with_node(new_connector);
					left_of_base->set_right_sib(right_of_base);
					new_connector->set_right_sib(child);
				}
				else {
					new_sib->bisect_edge_with_node(new_connector);
					new_sib->set_right_sib(child);
				}
				attachment_points.push_back(child);
				attachment_points.push_back(new_connector);
#				if ! defined(NDEBUG)
					tree.debug_check();
#				endif
			}
		}
	}
	return num_nodes_added;
}


////////////////////////////////////////////////////////////////////////////////
// SimTree manipulations
////////////////////
class SimNdBlob {
	public:
		SimNdBlob()
			:_edge_length(0.0),
			_sum_leaf_weights(1.0),
			_num_leaves_below(1),
			_num_nodes_at_or_below(1) {
		}

		double get_sum_leaf_weights() const {
			return this->_sum_leaf_weights;
		}
		void set_sum_leaf_weights(double x) {
			this->_sum_leaf_weights = x;
		}
		double get_edge_length() const {
			return _edge_length;
		}
		nnodes_t get_num_leaves_below() const {
			return this->_num_leaves_below;
		}
		void set_num_leaves_below(nnodes_t x) {
			this->_num_leaves_below = x;
		}
		nnodes_t get_num_nodes_at_or_below() const {
			return this->_num_nodes_at_or_below;
		}
		void set_num_nodes_at_or_below(nnodes_t x) {
			this->_num_nodes_at_or_below = x;
		}
		void debug_dump(std::ostream &o) const {
			o << " sum_leaf_weights=" << this->_sum_leaf_weights << " _num_leaves_below=" << this->_num_leaves_below << " _num_nodes_at_or_below=" << this->_num_nodes_at_or_below;
		}
	private:
		double _edge_length;
		double _sum_leaf_weights;
		nnodes_t _num_leaves_below;
		nnodes_t _num_nodes_at_or_below;
};
class SimTreeBlob {
	public:
		SimTreeBlob()
			:_sum_leaf_weights(-1.0) {
		}
		void set_is_dirty(bool) {
			this->_sum_leaf_weights = -1; //@TEMP currently using the sum of leaf weights
		}
		bool is_dirty() const {
			return (this->_sum_leaf_weights < 0); //@TEMP
		}


		double get_sum_leaf_weights() {
			return this->_sum_leaf_weights;
		}
		void set_sum_leaf_weights(double x) {
			this->_sum_leaf_weights = x;
		}
	private:
		double _sum_leaf_weights; // <0 is a flag for "is_dirty".  @TEMP
};
typedef Node<SimNdBlob> SimNode;
typedef Tree<SimNdBlob, SimTreeBlob> SimTree;

////////////////////////////////////////////////////////////////////////////////
// Main impl
////////////////////

////////////////////////////////////////////////////////////////////////////////
// Sums the sum_leaf_weights and num_leaves_below attributues of the node blob
//	Sets the sum_leaf_weights for the tree.
////////
template <typename T, typename U>
void sum_leaf_weights_over_tree(Tree<T,U> & tree) {
	typedef typename Tree<T,U>::Node_T::fast_postorder_iterator nd_it_t;
	typedef typename Node<T>::child_iterator child_it;
	for (nd_it_t it = tree.begin_fast_postorder(); it != tree.end_fast_postorder(); ++it) {
		Node<T> & nd = *it;
		if (nd.is_internal()) {
			double x = 0.0;
			nnodes_t num_leaves_below = 0;
			nnodes_t num_nodes_at_or_below = 1;
			for (child_it c_it = nd.begin_child(); c_it != nd.end_child(); ++c_it) {
				const T & blob = c_it->blob;
				x += blob.get_sum_leaf_weights();
				num_leaves_below += blob.get_num_leaves_below();
				num_nodes_at_or_below += blob.get_num_nodes_at_or_below();
			}
			nd.blob.set_sum_leaf_weights(x);
			nd.blob.set_num_leaves_below(num_leaves_below);
			nd.blob.set_num_nodes_at_or_below(num_nodes_at_or_below);
		}
	}
	assert(tree.get_root());
	tree.blob.set_sum_leaf_weights(tree.get_root()->blob.get_sum_leaf_weights());
}

template <typename T, typename U>
void refresh_blob_data(Tree<T,U> & tree) {
	//std::cerr << "In refresh_blob_data. Before...\n"; tree.get_root()->debug_dump(std::cerr, true);
	sum_leaf_weights_over_tree<T, U>(tree);
	//std::cerr << "In refresh_blob_data. After...\n"; tree.get_root()->debug_dump(std::cerr, true);

}


template <typename T, typename U>
Node<T> * find_leaf_by_weight_range(Tree<T,U> & tree, double x) {
	if (tree.blob.is_dirty()) {
		sum_leaf_weights_over_tree(tree);
	}
	if (x > tree.blob.get_sum_leaf_weights()) {
		return nullptr;
	}
	Node<T> * curr_nd = tree.get_root();
	for (;;) {
		typedef typename Node<T>::child_iterator ch_it;
		assert(curr_nd->is_internal());
		for (ch_it c_it = curr_nd->begin_child(); c_it != curr_nd->end_child(); ++c_it) {
			const double subtree_wt = c_it->blob.get_sum_leaf_weights();
			if (subtree_wt < x) {
				x -= subtree_wt;
			}
			else {
				curr_nd = &(*c_it);
				if (c_it->is_leaf()) {
					return curr_nd;
				}
				break;
			}
		}
	}
}

template <typename T>
const Node<T> * find_leaf_by_weight_range_offsets(const Node <T> & root,
									      double x,
										  const std::unordered_map<const Node <T> *, double> & wt_offset) {
	double curr_nd_offset = get_or_def_const(wt_offset, &root, 0.0);
	if (x > root.blob.get_sum_leaf_weights() - curr_nd_offset) {
		return nullptr;
	}
	std::string indent = " ";
	//std::cerr << indent <<  "find_leaf_by_weight_range_offsets x = " << x << " root weight = " << root.blob.get_sum_leaf_weights() << " - " << curr_nd_offset << "\n";
	const Node<T> * curr_nd = &root;
	for (;;) {
		typedef typename Node<T>::const_child_iterator ch_it;
		assert(curr_nd->is_internal());
		const Node<T> * prev = curr_nd;
		indent.append(2, ' ');
		for (ch_it c_it = curr_nd->begin_child(); c_it != curr_nd->end_child(); ++c_it) {
			const double subtree_wt = c_it->blob.get_sum_leaf_weights();
			const double subtree_offset = get_or_def_const(wt_offset, &(*c_it), 0.0);
			//std::cerr << indent <<  "before x = " << x << " subtree_wt = " << subtree_wt << " - " << subtree_wt << "\n";
			if (subtree_offset < subtree_wt) {
				const double corrected_wt = subtree_wt - subtree_offset;
				if (corrected_wt < x) {
					x -= corrected_wt;
					//std::cerr << indent <<  "decremented x => " << x << "\n";
				}
				else {
					curr_nd = &(*c_it);
					//std::cerr << indent <<  "chose nd curr_nd =" << curr_nd << "\n";
					if (c_it->is_leaf()) {
						//std::cerr << indent <<  "It is a leaf!\n";
						return curr_nd;
					}
					break;
				}
			}
		}
		if (curr_nd == nullptr or prev == curr_nd) {
			return nullptr;
		}
	}
}




typedef std::vector<std::string> ProgCommand;
class ProgState;
bool rec_and_process_command(const ProgCommand & command_vec,
					 		 ProgState & prog_state,
					 		 bool record);

////////////////////////////////////////////////////////////////////////////////
// Wraps up the changeable state of the program (as it runs through the command
//	line execution.
///////
class ProgState {
	public:
		ProgState(TaxonNameUniverse & taxa_universe,
				  SimTree & tree, RandGen::uint_seed_t seed)
			:err_stream(std::cerr),
			strict_mode(false),
			taxa(taxa_universe),
			curr_repeat_list(0L),
			outp(&std::cout),
			full_tree(tree),
			current_tree(nullptr),
			last_seed(seed),
			rng(seed),
			verbose(false) {
			this->current_tree = &(this->full_tree);
		}
		void sample_new_focal(const SimNode * src_root,
		                      const std::set<const SimNode *> & leaves,
		                      const std::set<const SimNode *> & traversed) {
		    if (this->current_tree != &(this->full_tree)) {
		    	SimTree * t = this->current_tree;
		    	this->current_tree = 0L;
		    	delete t;

		    }
		    this->current_tree = create_new_subsampled_tree(src_root,
		    												leaves,
		    												traversed,
		    												this->nd_blob,
		    												this->tree_blob,
		    												2*src_root->blob.get_num_leaves_below());
		}

		void print_tree(bool edge_lengths, bool nexus) const {
			if (this->outp != nullptr) {
				assert(this->current_tree);
				if (nexus) {
					*this->outp << "#NEXUS\nBEGIN TAXA;\n  Dimensions ntax = ";
					*this->outp << this->current_tree->get_num_leaves() << ";\n	 Taxlabels";
					for (SimNode::const_leaf_iterator sn_it = this->current_tree->begin_leaf();
						sn_it != this->current_tree->end_leaf();
						++sn_it) {
						*this->outp << ' ' << sn_it->get_label().newick();
					}
				*this->outp << ";\nEND;\n\nBEGIN TREES;\n  Tree one = [&R] ";
				}
				//this->err_stream << "Calling write_newick...\n";
				this->current_tree->write_newick(*this->outp, edge_lengths);
				if (nexus) {
					*this->outp << "\nEND;";
				}
				*this->outp << std::endl;
			}
		}

		std::ostream * get_output_stream() {
			return this->outp;
		}

		void set_output_stream( std::ostream * new_out) {
			if (&(this->outp_obj) == this->outp) {
				this->outp_obj.close();
			}
			this->outp = new_out;
		}

		void set_output_file(const char * new_out) {
			if (&(this->outp_obj) == this->outp) {
				this->outp_obj.close();
			}
			if (new_out) {
				this->outp_obj.open(new_out, std::ios::app);
				this->outp = &(this->outp_obj);
			}
			else {
				this->outp = 0L;
			}
		}

		bool out_good() const {
			return this->outp and this->outp->good();
		}

		std::ostream & err_stream;
		bool strict_mode;
		SimTree & get_full_tree() const {
			return this->full_tree;
		}
		SimTree * get_focal_tree() const {
			return this->current_tree;
		}
		TaxonNameUniverse & taxa;
		std::vector<unsigned> repeat_count_vec;
		typedef std::pair<bool, std::list<ProgCommand> > cmd_recorder_t;
		cmd_recorder_t scratch_repeat_list;
		cmd_recorder_t *curr_repeat_list;
		std::list<cmd_recorder_t> command_list_list;
		unsigned max_tries = 100;

	private:
		std::ostream * outp;
		std::ofstream outp_obj;
		SimTree & full_tree;
		SimTree * current_tree;
		RandGen::uint_seed_t last_seed;
	public:
		RandGen rng;
		SimNdBlob nd_blob;
		SimTreeBlob tree_blob;
		bool verbose;

		ProgState(const ProgState & other); // not defined, not copyable
		ProgState & operator=(const ProgState & other);// not defined, not copyable
};





std::string capitalize(const std::string & x) {
	std::string cmd = x;
	std::transform(cmd.begin(), cmd.end(), cmd.begin(), ::toupper);
	return cmd;
}

bool unrecognize_arg(const char * cmd, const char * arg, ProgState & prog_state) {
	prog_state.err_stream << "Unrecognized argument \"" << arg << "\" to the \"" << cmd << "\" command.\n";
	return !prog_state.strict_mode;
}



std::pair<bool, long> parse_pos_int(ProgState & prog_state,
									unsigned arg_ind,
									const ProgCommand & command_vec,
									const char * arg_name) {
	std::pair<bool, long> r(false, 0L);
	if (arg_ind + 2 >= command_vec.size() or command_vec[1 + arg_ind] != "=") {
		prog_state.err_stream << "Expecting  = # after " << arg_name << "\n";
		return r;
	}
	char * end_ptr;
	std::string s = command_vec[2 + arg_ind];
	long count = std::strtol(s.c_str(), &end_ptr, 10);
	if (end_ptr != s.c_str() + s.length() or count < 1) {
		prog_state.err_stream << "Expected a positive number (the weight) after the ROOTMIN command. found " << s << ".\n";
		return r;
	}
	r.first = true;
	r.second = count;
	return r;
}

bool sample_leaves_and_traversed(const SimNode & root,
								 const nnodes_t num_leaves,
								 std::set<const SimNode *> * sampled,
								 std::set<const SimNode *> * traversed,
								 const SimNode * to_include,
								 const SimNode * avoid_node, // if we want to avoid a clade (e.g. avoid the ingroup) send in its pointer here...
								 RandGen & rng,
								 ProgState & prog_state) {
	assert(num_leaves <= root.blob.get_num_leaves_below());
	if (num_leaves == root.blob.get_num_leaves_below()) {
		// grab each leaf
		for (SimNode::const_fast_postorder_iterator nd_it = root.begin_fast_postorder(); nd_it != root.end_fast_postorder(); ++nd_it) {
			const SimNode * nd_ptr = &(*nd_it);
			assert(nd_ptr);
			if (nd_ptr->is_leaf()) {
				sampled->insert(nd_ptr);
			}
			else {
				traversed->insert(nd_ptr);
			}
		}
		return true;
	}
	if (num_leaves == 0) {
		return true;
	}
	//std::cerr << " Sampling ="<< num_leaves << " leaves from a clade of " << root.blob.get_num_leaves_below() << " leaves\n";
	std::unordered_map<const SimNode *, double> weight_offset;
	nnodes_t num_added = 0;
	std::list<const SimNode *> path_to_anc;
	std::list<const SimNode *>::const_iterator p2a_it;
	typedef std::back_insert_iterator<std::list<const SimNode *> > back_it_t;

	if (avoid_node) {
		//std::cerr << " Avoiding a clade of "<< avoid_node->blob.get_num_leaves_below() << " leaves.\n";
		back_it_t back_i = back_inserter(path_to_anc);
		if (!avoid_node->put_path_to_anc(&root, back_i)) {
			assert(false);
			//std::cerr << "Ancestor of avoid_node node not found!\n";
			return false;
		}
		const double wt = avoid_node->blob.get_sum_leaf_weights();
		for (p2a_it = path_to_anc.begin(); p2a_it != path_to_anc.end(); ++p2a_it) {
			const SimNode * nd_ptr2 = *p2a_it;
			weight_offset[nd_ptr2] += wt;
			//std::cerr << "   weight="<< nd_ptr2->blob.get_sum_leaf_weights() << " offset[" << nd_ptr2 << "] = " << weight_offset[nd_ptr2] << "\n";
		}
	}

	if (to_include) {
		assert(to_include->is_leaf());
		back_it_t back_i = back_inserter(path_to_anc);
		if (!to_include->put_path_to_anc(&root, back_i)) {
			assert(false);
			//std::cerr << "Ancestor of leaf node not found!\n";
			return false;
		}
		++num_added;

		const double wt = to_include->blob.get_sum_leaf_weights();
		for (p2a_it = path_to_anc.begin(); p2a_it != path_to_anc.end(); ++p2a_it) {
			const SimNode * nd_ptr2 = *p2a_it;
			if (nd_ptr2->is_leaf()) {
				sampled->insert(nd_ptr2);
			}
			else {
				traversed->insert(nd_ptr2);
			}
			weight_offset[nd_ptr2] += wt;
			//std::cerr << "   weight="<< nd_ptr2->blob.get_sum_leaf_weights() << " offset[" << nd_ptr2 << "] = " << weight_offset[nd_ptr2] << "\n";
		}
		//std::cerr << "  Orig leaf added\n";
	}

	while (num_added < num_leaves) {
		//std::cerr << "  Adding leaf " << 1 + num_added << "\n";
		double u = rng.uniform01();
		double x = u * (root.blob.get_sum_leaf_weights() - weight_offset[&root]);
		//std::cerr << "  u="<< u << " root.weight=" << root.blob.get_sum_leaf_weights() << " root offset=" << weight_offset[&root] << " x=" << x <<  "\n";
		const SimNode * next_leaf = nullptr;
		for (unsigned trial = 0; next_leaf == nullptr and trial < prog_state.max_tries; ++trial) {
			next_leaf = find_leaf_by_weight_range_offsets(root,
														  x,
														  weight_offset);
		}
		if (next_leaf == nullptr) {
			prog_state.err_stream << "Could not find a random leaf - this could be caused by rounding error if the leaf weights vary dramatically!\n";
			return false;
		}
		path_to_anc.clear();
		back_it_t back_i = back_inserter(path_to_anc);
		if (!next_leaf->put_path_to_anc(&root, back_i)) {
			assert(false);
			//std::cerr << "Ancestor of rand leaf node not found!\n";
			return false;
		}
		++num_added;
		const double wt = next_leaf->blob.get_sum_leaf_weights();
		for (p2a_it = path_to_anc.begin(); p2a_it != path_to_anc.end(); ++p2a_it) {
			const SimNode * nd_ptr3 = *p2a_it;
			if (nd_ptr3->is_leaf()) {
				sampled->insert(nd_ptr3);
			}
			else {
				traversed->insert(nd_ptr3);
			}
			weight_offset[nd_ptr3] += wt;
			//std::cerr << "   weight="<< nd_ptr2->blob.get_sum_leaf_weights() << " offset[" << nd_ptr2 << "] = " << weight_offset[nd_ptr2] << "\n";
		}
	}
	return true;
}

////////////////////////////////////////////////////////////////////////////////
// all nodes have an equal probability of being chosen. The root will not be
// returned
SimNode * choose_rand_node_below_exclusive(SimNode & root, ProgState & prog_state) {
	//std::cerr << "choose_rand_node_below_exclusive\n"; root.debug_dump(prog_state.err_stream, true);
	nnodes_t num_nodes_at_or_below_root = root.blob.get_num_nodes_at_or_below();
	if (num_nodes_at_or_below_root <= 1) {
		return nullptr;
	}
	// start at 1 to avoid root
	nnodes_t node_ind = prog_state.rng.rand_range(1, num_nodes_at_or_below_root);
	SimNode * curr_node = &root;
	std::stack<SimNode *> nd_stack;
	for (;;) {
		//std::cerr<< "node_ind = " << node_ind << "\ncurr_node info: "; curr_node->debug_dump(prog_state.err_stream, false);
		if (node_ind == 0) {
			return curr_node;
		}
		node_ind--; // subtract 1 for the current Node
		curr_node = curr_node->get_left_child();
		for (;;) {
			assert(curr_node);
			if (node_ind < curr_node->blob.get_num_nodes_at_or_below()) {
				break;
			}
			node_ind -= curr_node->blob.get_num_nodes_at_or_below();
			curr_node = curr_node->get_right_sib();
			//std::cerr<< "inner for loop: node_ind = " << node_ind << "\ncurr_node info: "; curr_node->debug_dump(prog_state.err_stream, false);
		}
	}
}


std::vector<SimNode *>	get_internal_children(SimNode * par, SimNode *avoid_node) {
	std::vector<SimNode *> child_vec = par->get_children();
	typedef std::vector<SimNode *>::const_iterator vec_node_ptr_it;
	std::vector<SimNode *> poss_up;
	poss_up.reserve(child_vec.size() - 1);
	for (vec_node_ptr_it c_it = child_vec.begin(); c_it != child_vec.end(); ++c_it) {
		SimNode * nd;
		nd = (*c_it);
		if (nd != avoid_node and nd->is_internal()) {
			poss_up.push_back(nd);
		}
	}
	return poss_up;
}
nnodes_t count_internal_children(SimNode * par, SimNode *avoid_node) {
	std::vector<SimNode *> child_vec = get_internal_children(par, avoid_node);
	return child_vec.size();
}

bool slide_node_lteq_num_edges(SimNode * moving_nd, nnodes_t max_recon_dist, SimNode * root, ProgState & prog_state) {
#	if ! defined(NDEBUG)
		root->debug_check_subtree_nav_pointers();
#	endif
	if (max_recon_dist == 0) {
		return true;
	}
	assert(moving_nd);
	SimNode * par =  moving_nd->get_parent();
	assert(par);
	SimNode * curr_attach = nullptr;
	SimNode * avoid_node = nullptr;
	bool moving_up;
	std::vector<SimNode *> child_vec;
	std::vector<SimNode *> poss_up;
	const bool attach_at_polytomy = par->is_polytomy();
	if (attach_at_polytomy) {
		// treat par as attachment, only move to internal nodes.
		curr_attach = par;
	}
	else {
		// treat sister as attachment, move to any node, attach below on exit, .
		curr_attach = moving_nd->get_right_sib();
		if (curr_attach == nullptr) {
			curr_attach = par->get_left_child();
			assert(curr_attach and curr_attach != moving_nd);
		}
		if (par == root) { // what a pain, we have to keep the root stable
			if (curr_attach->is_leaf()) {
				return true; // no move is possible, this will be a no-op swap
			}
		}
	}
	// first decide whether moving up or down...
	nnodes_t num_poss_up = (attach_at_polytomy ? count_internal_children(par, moving_nd) : curr_attach->get_num_children());
	const nnodes_t num_poss_down = 1;
	double prob_up = ((double) num_poss_up) / ((double)(num_poss_down + num_poss_up));
	moving_up = (prog_state.rng.uniform01() < prob_up ? true : false);

	// detach moving_nd from its polytomy
	SimNode * connector = nullptr;
	if (attach_at_polytomy) {
		if (par->get_left_child() == moving_nd) {
			par->set_left_child(moving_nd->get_right_sib());
		}
		else {
			SimNode * l_sib = moving_nd->find_left_sib();
			assert(l_sib);
			l_sib->set_right_sib(moving_nd->get_right_sib());
		}
	}
	else {
		// also detach the connector node to use when reattaching...
		if (par == root) { // what a pain, we have to keep the root stable
			assert(curr_attach->is_internal()); // we bale earlier if this is false
			connector = curr_attach;
			curr_attach = root;
			SimNode * n = connector->get_left_child();
			root->set_left_child_raw(n);
			while (n) {
				n->set_parent_raw(root);
				n = n->get_right_sib();
			}
			moving_up = true;
		}
		else {
			connector = par;
			SimNode * grand_parent = par->get_parent();
			curr_attach->set_parent_raw(grand_parent);
			assert(curr_attach->get_right_sib() == nullptr or curr_attach->get_right_sib() == moving_nd);
			curr_attach->set_right_sib_raw(par->get_right_sib());
			if (grand_parent->get_left_child() == connector) {
				grand_parent->set_left_child_raw(curr_attach);
				assert(curr_attach->get_right_sib() != nullptr and curr_attach->get_right_sib() != moving_nd);
			}
			else {
				SimNode * pls = par->find_left_sib();
				pls->set_right_sib_raw(curr_attach);
			}
		}
		connector->set_parent_raw(nullptr);
		moving_nd->set_parent_raw(connector);
		connector->set_left_child_raw(moving_nd);
		moving_nd->set_right_sib_raw(nullptr);
	}
#	if ! defined(NDEBUG)
		root->debug_check_subtree_nav_pointers();
#	endif

	// navigate to find its new attachment point

	while (max_recon_dist > 0) {
		if (moving_up) {
			poss_up = (attach_at_polytomy ? get_internal_children(curr_attach, avoid_node) : curr_attach->get_children());
			if (poss_up.empty()) {
				break;
			}
			curr_attach = avoid_node;
			if (poss_up.size() == 1) {
				assert(poss_up[0] and poss_up[0] != avoid_node);
			}
			else {
				while (curr_attach == avoid_node) {
					curr_attach = choose_element(poss_up, prog_state.rng);
				}
			}
			avoid_node = nullptr;
		}
		else {
			if (curr_attach == root) {
				moving_up = true;
				++max_recon_dist; // increment to cancel the decrement and make this pivot a no-op
			}
			else {
				avoid_node = curr_attach;
				par = curr_attach->get_parent();
				if (par == root and !attach_at_polytomy) {
					curr_attach = root;
					moving_up = true;
					++max_recon_dist; // increment to cancel the decrement and make this pivot a no-op
				}
				else {
					num_poss_up = (attach_at_polytomy ? count_internal_children(par, avoid_node) : (par->get_num_children() - 1));
					double prob_up = ((double) num_poss_up) / ((double)(num_poss_down + num_poss_up));
					moving_up = (prog_state.rng.uniform01() < prob_up ? true : false);
					curr_attach = par;
				}
			}
		}
		max_recon_dist--;
	}
	assert(curr_attach);
	assert(moving_nd);
	if (attach_at_polytomy) {
		curr_attach->add_new_child(moving_nd);
	}
	else {
		assert(connector);
		assert(curr_attach != root);
		assert(curr_attach->get_parent());
		curr_attach->bisect_edge_with_node(connector);
		connector->add_new_child(moving_nd);
	}
#	if ! defined(NDEBUG)
		root->debug_check_subtree_nav_pointers();
#	endif
	return true;
}

bool do_rand_spr(SimTree & tree, const nnodes_t recon_min, const nnodes_t recon_max, ProgState & prog_state) {
	if (tree.blob_is_dirty()) {
		refresh_blob_data(tree);
	}
	SimNode * moving_nd = choose_rand_node_below_exclusive(*(tree.get_root()), prog_state);
	nnodes_t recon_dist = (nnodes_t) (recon_min == recon_max ? recon_min : prog_state.rng.rand_range(recon_min, recon_max));
	bool success = slide_node_lteq_num_edges(moving_nd, recon_dist, tree.get_root(), prog_state);
	tree.flag_blob_as_dirty(); // temp should update blob, rather just flagging it as dirty...
	return success;
}

bool do_sample(ProgState & prog_state,
			  const unsigned arg_root_min,
			  const unsigned arg_root_max,
			  const nnodes_t in_min,
			  const nnodes_t arg_in_max,
			  const nnodes_t out_min,
			  const nnodes_t arg_out_max) {
	SimTree & tree = prog_state.get_full_tree();
	if (tree.blob_is_dirty()) {
		refresh_blob_data(tree);
	}
	const double w = tree.blob.get_sum_leaf_weights();
	// prog_state.err_stream << "tree.blob.get_sum_leaf_weights() = ";
// 	prog_state.err_stream.setf(std::ios::fixed);
// 	prog_state.err_stream.precision(5);
// 	prog_state.err_stream << w << '\n';
// 	prog_state.err_stream << "tree.blob.get_num_leaves_below() = ";
// 	prog_state.err_stream << tree.get_root()->blob.get_num_leaves_below() << '\n';
// 	prog_state.err_stream << "Will try to sample: root depth ["<< arg_root_min << ", " << arg_root_max << "]\n";
// 	prog_state.err_stream << "                    ingroup    ["<< in_min << ", " << arg_in_max << "]\n";
// 	prog_state.err_stream << "                    outgroup   ["<< out_min << ", " << arg_out_max << "]\n";

	for (unsigned trial = 0; trial < prog_state.max_tries; ++trial) {
		if (trial > 0) {
			prog_state.err_stream << " sample trial failed. Trying again...\n";
		}
		// Step 1: Choose a leaf (using the leaf weighting)
		double rand_x = prog_state.rng.uniform01() * w;
		const SimNode * leaf_nd = find_leaf_by_weight_range(tree, rand_x);
		//prog_state.err_stream << "Chose \"" << leaf_nd->get_label().c_str() << "\"\n";
		unsigned root_depth;
		unsigned root_min = arg_root_min;
		unsigned root_max = arg_root_max;

		// Step 2: choose the depth of the root
		const SimNode * ingroup_root = nullptr;
		while (ingroup_root == nullptr) {
			if (root_min < root_max) {
				root_depth = prog_state.rng.rand_range(root_min, root_max);
			}
			else {
				root_depth = root_min;
			}
			unsigned x = 0;
			ingroup_root = leaf_nd->get_ancestor_by_edge_dist(root_depth, &x);
			if (ingroup_root == nullptr or ingroup_root->get_parent() == nullptr) { // reached the root of the tree in the ingroup
				ingroup_root = nullptr;
				root_max = x - 1;
				if (prog_state.verbose) {
					prog_state.err_stream << "  Cropping root depth ["<< root_min << ", " << root_max << "]\n";
				}
				if (root_max < root_min) {
					break;
				}
			}
			if (ingroup_root) {
				nnodes_t potential_ingroup = ingroup_root->blob.get_num_leaves_below();
				if (potential_ingroup < in_min) { // not enough leaves at this depth, increase the min depth
					ingroup_root = nullptr;
					root_min = root_depth + 1;
					if (prog_state.verbose) {
						prog_state.err_stream << "  Bumping up root_min ["<< root_min << ", " << root_max << "]\n";
					}
					if (root_max < root_min) {
						break;
					}
				}
				if (ingroup_root) {
					nnodes_t potential_outgroup = ingroup_root->get_parent()->blob.get_num_leaves_below() - potential_ingroup;
					if (potential_outgroup < out_min) { // not enough in the outgroup. Reject this draw, but don't change params (this result depends on the size of the sister group, which can fluctuate freely with the depth of the tree.
						ingroup_root = nullptr;
					}
				}
			}
		}
		if (ingroup_root == nullptr) {
			continue;
		}
		const SimNode * sample_root = ingroup_root->get_parent();
		assert(sample_root);

		//Step 3: choose the ingroup
		const nnodes_t num_below_sr = sample_root->blob.get_num_leaves_below();
		const nnodes_t num_possible_ingroup = ingroup_root->blob.get_num_leaves_below();
		const nnodes_t num_possible_outgroup = num_below_sr - num_possible_ingroup;
		//prog_state.err_stream << " sample_root found num_possible_ingroup="<< num_possible_ingroup << ", num_below_sr=" << num_below_sr << "\n";
		std::set<const SimNode *> leaf_nodes;
		std::set<const SimNode *> traversed_internals_nodes;
		assert(num_possible_ingroup >= in_min);
		nnodes_t in_max = arg_in_max;
		if (num_possible_ingroup < in_max) {
			in_max = num_possible_ingroup;
		}
		nnodes_t in_size = (in_max == in_min ? in_min : prog_state.rng.rand_range(in_min, in_max));
		if (!sample_leaves_and_traversed(*ingroup_root,
										 in_size,
										 &leaf_nodes,
										 &traversed_internals_nodes,
										 leaf_nd,
										 nullptr,
										 prog_state.rng,
										 prog_state)) {
			continue;
		}

		traversed_internals_nodes.insert(ingroup_root);
		traversed_internals_nodes.insert(sample_root);

		//Step 4: choose the outgroup
		assert(num_possible_outgroup >= out_min);
		nnodes_t out_max = arg_out_max;
		if (num_possible_outgroup < out_max) {
			out_max = num_possible_outgroup;
		}
		nnodes_t out_size = (out_max == out_min ? out_min : prog_state.rng.rand_range(out_min, out_max));
		if (!sample_leaves_and_traversed(*sample_root,
										 out_size,
										 &leaf_nodes,
										 &traversed_internals_nodes,
										 nullptr,
										 ingroup_root,
										 prog_state.rng,
										 prog_state)) {
			continue;
		}

		//std::cerr << "Full Tree\n";
		//sample_root->debug_dump(std::cerr, true);
		//std::cerr << "Leaves\n";
		//for (auto nd_it : leaf_nodes) {
		//	nd_it->debug_dump(std::cerr, true);
		//}
		//std::cerr << "traversed_internals_nodes:\n";
		//for (auto nd_it : traversed_internals_nodes) {
		//	nd_it->debug_dump(std::cerr, false);
		//}
		prog_state.sample_new_focal(sample_root, leaf_nodes, traversed_internals_nodes);
		if (prog_state.get_focal_tree()) {
			prog_state.get_focal_tree()->flag_blob_as_dirty();
			return true;
		}
	}

	prog_state.err_stream << prog_state.max_tries << " sample trials failed. bailing out...\n";
	return !prog_state.strict_mode;
}


bool process_resolve_command(const ProgCommand & command_vec,
						   ProgState & prog_state) {
	if (command_vec.size() > 1) {
		if (!unrecognize_arg("RESOLVE", command_vec.at(1).c_str(), prog_state)) {
			return false;
		}
	}

	SimTree * focal_tree = prog_state.get_focal_tree();
	assert(focal_tree);
	resolve_polytomies(*focal_tree, prog_state.rng);
	if (prog_state.verbose) {
		prog_state.err_stream << "RESOLVE successful.\n";
	}
	focal_tree->flag_blob_as_dirty();
	return true;
}
bool process_status_command(const ProgCommand & command_vec,
						   ProgState & prog_state) {
	if (command_vec.size() > 1) {
		if (!unrecognize_arg("STATUS", command_vec.at(1).c_str(), prog_state)) {
			return false;
		}
	}
	SimTree & full_tree = prog_state.get_full_tree();
	SimTree * focal_tree = prog_state.get_focal_tree();

	if (full_tree.blob_is_dirty()) {
		refresh_blob_data(full_tree);
	}
	if (focal_tree and focal_tree->blob_is_dirty()) {
		refresh_blob_data(*focal_tree);
	}
	prog_state.err_stream << "Full tree has " << full_tree.get_root()->blob.get_num_leaves_below() << " leaves and " << full_tree.get_root()->blob.get_num_nodes_at_or_below() << " nodes.\n";
	prog_state.err_stream << "The sum of its leaf weights is " << full_tree.blob.get_sum_leaf_weights() << " .\n";
	if (focal_tree == nullptr or focal_tree == &full_tree) {
		prog_state.err_stream << "No subsampled tree\n";
	}
	else {
		prog_state.err_stream << "The subsampled tree has " << focal_tree->get_root()->blob.get_num_leaves_below() << " leaves and " << focal_tree->get_root()->blob.get_num_nodes_at_or_below() << " nodes.\n";
		prog_state.err_stream << "The sum of its leaf weights is " << focal_tree->blob.get_sum_leaf_weights() << " .\n";
	}

	return true;
}

bool process_repeat_command(const ProgCommand & command_vec,
						   ProgState & prog_state) {
	if (command_vec.size() != 2) {
		std::cerr << "Expected \"REPEAT #### ;\" where ### is an integer\n";
		return !prog_state.strict_mode;
	}
	if (prog_state.curr_repeat_list) {
		std::cerr << "Expected \"REPEAT\" commands cannot be nested! (in the current implementation)\n";
		return false;
	}
	std::string count_str = command_vec[1];
	char * e;
	long count = std::strtol(count_str.c_str(), &e, 10);
	if (e != count_str.c_str() + count_str.length()) {
		prog_state.err_stream << "Expected a number (the weight) after the REPEAT command. found " << count_str << ".\n";
		return !prog_state.strict_mode;
	}
	if (count < 1) {
		prog_state.err_stream << "Expected a positive number (the weight) after the REPEAT command. found " << count_str << ".\n";
		return !prog_state.strict_mode;
	}
	prog_state.repeat_count_vec.push_back((unsigned) count);
	if (prog_state.curr_repeat_list) {
		prog_state.command_list_list.push_back(*prog_state.curr_repeat_list);
		prog_state.curr_repeat_list->second.clear();
		prog_state.curr_repeat_list->first = true;
	}
	else {
		prog_state.curr_repeat_list = &(prog_state.scratch_repeat_list);
		prog_state.curr_repeat_list->second.clear();
		prog_state.curr_repeat_list->first = true;
	}
	return true;
}


bool process_end_repeat_command(const ProgCommand & command_vec,
						   ProgState & prog_state) {
	if (command_vec.size() != 1) {
		std::cerr << "Expected no arguments after \"ENDREPEAT ;\" found \"" <<command_vec.at(1) << "\"\n";
		return !prog_state.strict_mode;
	}
	if (!prog_state.curr_repeat_list) {
		std::cerr << "Expecting \"REPEAT\" command before \"ENDREPEAT\" \n";
		return !prog_state.strict_mode;
	}
	assert(!prog_state.repeat_count_vec.empty());
	assert(*prog_state.repeat_count_vec.rbegin() > 0);
	unsigned end_ind = 2;

	while (end_ind > 1) {
		end_ind = (*prog_state.repeat_count_vec.rbegin()) - 1;
		//std::cerr << "process_end_repeat_command loop end_ind = " << end_ind <<'\n';
		if (end_ind >= 1) {
			*prog_state.repeat_count_vec.rbegin() = end_ind;
			prog_state.curr_repeat_list->first = false;
			std::list<ProgCommand> to_repeat = prog_state.curr_repeat_list->second;
			// we repeat the loop that we are currently in. This could be recursive if loops are nested

			//std::cerr << "About to repeat: {\n";
			//for (std::list<ProgCommand>::const_iterator c_it = to_repeat.begin(); c_it != to_repeat.end(); ++c_it) {
			//	const ProgCommand & cmd = *c_it;
			//std::cerr << "  " << cmd[0] << '\n';
			//}
			//std::cerr << " }\n";

			for (std::list<ProgCommand>::const_iterator c_it = to_repeat.begin(); c_it != to_repeat.end(); ++c_it) {
				const ProgCommand & cmd = *c_it;
				rec_and_process_command(cmd, prog_state, false);
			}
		}
		else {
			prog_state.repeat_count_vec.pop_back();
			assert(prog_state.curr_repeat_list);
			if (prog_state.command_list_list.empty()) {
				prog_state.curr_repeat_list = nullptr;
			}
			else {
				*prog_state.curr_repeat_list = *prog_state.command_list_list.rbegin();
				prog_state.command_list_list.pop_back();
			}
		}
	}
	return true;
}

bool process_spr_command(const ProgCommand & command_vec,
						   ProgState & prog_state) {
	std::pair<bool, long> r;
	SimTree * tree = prog_state.get_focal_tree();
	assert(tree and tree->get_root());
	// assure that the tip counts are right... (side effect of sum_leaf_weights_over_tree)
	if (tree->blob_is_dirty()) {
		refresh_blob_data(*tree);
	}
	nnodes_t recon_min = 1;
	nnodes_t recon_max = 2*tree->get_root()->blob.get_num_leaves_below();
	unsigned num_moves = 1;
	for (unsigned arg_ind = 1; arg_ind < command_vec.size(); ++arg_ind) {
		const std::string cap = capitalize(command_vec[arg_ind]);
		if (cap == "REP") {
			r = parse_pos_int(prog_state, arg_ind, command_vec, "REP");
			if (!r.first) {
				return !prog_state.strict_mode;
			}
			num_moves = (unsigned)r.second;
			arg_ind += 2;
		}
		else if (cap == "RECONMIN") {
			r = parse_pos_int(prog_state, arg_ind, command_vec, "RECONMAX");
			if (!r.first) {
				return !prog_state.strict_mode;
			}
			recon_min = (nnodes_t)r.second;
			arg_ind += 2;
		}
		else if (cap == "RECONMAX") {
			r = parse_pos_int(prog_state, arg_ind, command_vec, "RECONMAX");
			if (!r.first) {
				return !prog_state.strict_mode;
			}
			recon_max = (nnodes_t)r.second;
			arg_ind += 2;
		}
		else {
			if (!unrecognize_arg("SPR", command_vec.at(arg_ind).c_str(), prog_state)) {
				return false;
			}
		}
	}
	if (recon_min > recon_max) {
		prog_state.err_stream << "ReconMax cannot be below ReconMin!\n";
		return !prog_state.strict_mode;
	}
	for (unsigned i = 0; i < num_moves; ++i) {
		if (!do_rand_spr(*tree, recon_min, recon_max, prog_state)) {
			return false; // failure in tree manip is fatal.
		}
	}
	if (prog_state.verbose) {
		prog_state.err_stream << "RESOLVE successful.\n";
	}
	return true;
}

bool process_set_command(const ProgCommand & command_vec,
						   ProgState & prog_state) {
	std::pair<bool, long> r;
	for (unsigned arg_ind = 1; arg_ind < command_vec.size(); ++arg_ind) {
		const std::string cap = capitalize(command_vec[arg_ind]);
		if (cap == "SEED") {
			r = parse_pos_int(prog_state, arg_ind, command_vec, "SEED");
			if (!r.first) {
				return !prog_state.strict_mode;
			}
			if (prog_state.verbose) {
				prog_state.err_stream << "Setting seed to " << r.second << "\n";

			}
			prog_state.rng.seed(r.second);
			arg_ind += 2;
		}
		else if (cap == "MAXTRIES") {
			r = parse_pos_int(prog_state, arg_ind, command_vec, "MAXTRIES");
			if (!r.first) {
				return !prog_state.strict_mode;
			}
			if (prog_state.verbose) {
				prog_state.err_stream << "Setting MaxTries to " << r.second << "\n";

			}
			prog_state.max_tries = r.second;
			arg_ind += 2;
		}
		else if (cap == "STRICT") {
			if (prog_state.verbose) {
				prog_state.err_stream << "Setting execution to strict mode.\n";
			}
			prog_state.strict_mode = true;
		}
		else if (cap == "VERBOSE") {
			prog_state.err_stream << "Verbose mode enabled.\n";
			prog_state.verbose = true;
		}
		else {
			if (!unrecognize_arg("SET", command_vec.at(arg_ind).c_str(), prog_state)) {
				return false;
			}
		}
	}
	return true;
}
bool process_weight_command(const ProgCommand & command_vec,
						   ProgState & prog_state) {
	if (command_vec.size() <= 2) {
		std::cerr << "Expected a weight then a filepath of names after the WEIGHT command.\n";
		return !prog_state.strict_mode;
	}
	std::string weight_string = command_vec[1];
	char * e;
	double wt = std::strtod(weight_string.c_str(), &e);
	if (e != weight_string.c_str() + weight_string.length()) {
		prog_state.err_stream << "Expected a number (the weight) after the WEIGHT command. found " << weight_string << ".\n";
		return !prog_state.strict_mode;
	}
	if (wt < 0.0) {
		prog_state.err_stream << "Expected a non-negative number (the weight) after the WEIGHT command. found " << weight_string << ".\n";
		return !prog_state.strict_mode;
	}
	if (command_vec.size() > 3) {
		if (!unrecognize_arg("WEIGHT", command_vec.at(3).c_str(), prog_state)) {
			return false;
		}
	}


	std::string fn = command_vec[2];
	std::ifstream inp(fn);
	if (!inp.good()) {
		prog_state.err_stream << "Could not open the file " << fn << ".\n";
		return !prog_state.strict_mode;
	}
	const std::vector<TaxonLabel> labels = parse_labels_from_stream(inp, prog_state.taxa);
	SimTree & tree = prog_state.get_full_tree();
	std::vector<TaxonLabel>::const_iterator lab_it = labels.begin();
	tree.blob.set_sum_leaf_weights(-1);
	for (; lab_it != labels.end(); ++lab_it) {
		SimNode * nd = tree.find_node_by_label(*lab_it);
		if (nd == nullptr) {
			prog_state.err_stream << "Label \"" << lab_it->get_label() << "\" not found in tree. Ignored...\n";
		}
		else if (nd->is_internal()) {
			for (auto l_it = nd->begin_leaf(); l_it != nd->end_leaf(); ++l_it) {
				l_it->blob.set_sum_leaf_weights(wt);
			}
		}
		else {
			nd->blob.set_sum_leaf_weights(wt);
		}
	}
	tree.flag_blob_as_dirty();
	return true;
}


bool process_sample_command(const ProgCommand & command_vec,
						   ProgState & prog_state) {
	int root_min = 1;
	int root_max = std::numeric_limits<int>::max();
	int in_min = 2;
	int in_max = 2;
	int out_min = 1;
	int out_max = 1;
	std::pair<bool, long> r;
	for (unsigned arg_ind = 1; arg_ind < command_vec.size(); ++arg_ind) {
		const std::string cap = capitalize(command_vec[arg_ind]);
		if (cap == "ROOTMIN") {
			r = parse_pos_int(prog_state, arg_ind, command_vec, "ROOTMIN");
			if (!r.first) {
				return !prog_state.strict_mode;
			}
			root_min = r.second;
			arg_ind += 2;
		}
		else if (cap == "ROOTMAX") {
			r = parse_pos_int(prog_state, arg_ind, command_vec, "ROOTMAX");
			if (!r.first) {
				return !prog_state.strict_mode;
			}
			root_max = r.second;
			arg_ind += 2;
		}
		else if (cap == "INMIN") {
			r = parse_pos_int(prog_state, arg_ind, command_vec, "INMIN");
			if (!r.first) {
				return !prog_state.strict_mode;
			}
			in_min = r.second;
			arg_ind += 2;
		}
		else if (cap == "INMAX") {
			r = parse_pos_int(prog_state, arg_ind, command_vec, "INMAX");
			if (!r.first) {
				return !prog_state.strict_mode;
			}
			in_max = r.second;
			arg_ind += 2;
		}
		else if (cap == "OUTMIN") {
			r = parse_pos_int(prog_state, arg_ind, command_vec, "OUTMIN");
			if (!r.first) {
				return !prog_state.strict_mode;
			}
			out_min = r.second;
			arg_ind += 2;
		}
		else if (cap == "OUTMAX") {
			r = parse_pos_int(prog_state, arg_ind, command_vec, "OUTMAX");
			if (!r.first) {
				return !prog_state.strict_mode;
			}
			out_max = r.second;
			arg_ind += 2;
		}
		else {
			if (!unrecognize_arg("SAMPLE", command_vec.at(arg_ind).c_str(), prog_state)) {
				return false;
			}
		}
	}
	if (root_min > root_max) {
		std::cerr << "ROOTMAX must be at least as large as ROOTMIN\n";
		return !prog_state.strict_mode;
	}
	if (in_min > in_max) {
		std::cerr << "INMAX must be at least as large as INMIN\n";
		return !prog_state.strict_mode;
	}
	if (out_min > out_max) {
		std::cerr << "OUTMAX must be at least as large as OUTMIN\n";
		return !prog_state.strict_mode;
	}
	if (do_sample(prog_state, root_min, root_max, in_min, in_max, out_min, out_max)) {
		if (prog_state.get_focal_tree()) {
			prog_state.get_focal_tree()->flag_blob_as_dirty(); // temp should update blob, rather just flagging it as dirty...
		}
		if (prog_state.verbose) {
			prog_state.err_stream << "SAMPLE successful.\n";
		}

		return true;
	}
	if (prog_state.verbose) {
		prog_state.err_stream << "SAMPLE unsuccessful.\n";
	}
	return !prog_state.strict_mode;
}



bool process_print_command(const ProgCommand & command_vec,
						   ProgState & prog_state) {
	bool nexus = false;
	if (command_vec.size() > 1) {
		if (capitalize(command_vec[1]) == "NEXUS") {
			nexus = true;
		}
		else {
			if (!unrecognize_arg("PRINT", command_vec.at(1).c_str(), prog_state)) {
				return false;
			}
		}
	}
	prog_state.print_tree(false, nexus);
	if (prog_state.verbose) {
		prog_state.err_stream << "PRINT successful.\n";
	}
	return true;
}
bool process_out_command(const ProgCommand &  command_vec,
						 ProgState & prog_state) {
	if (command_vec.size() == 1) {
		prog_state.err_stream << "Expecting an argument (FILE, STD, or STOP) after OUT command.";
		return !prog_state.strict_mode;
	}
	std::string arg = capitalize(command_vec[1]);
	if (arg == "FILE") {
		if (command_vec.size() < 4 or command_vec.at(2) != "=") {
			prog_state.err_stream << "Expecting filepaths to be specified as OUT FILE = filepath\n";
			return !prog_state.strict_mode;
		}
		const char * filepath = command_vec.at(3).c_str();
		prog_state.set_output_file(filepath);
		if (! prog_state.out_good()) {
			prog_state.err_stream << "Could not open the file \"" << filepath << "\". Output set to standard out!\n";
			prog_state.set_output_stream(&std::cout);
			return !prog_state.strict_mode;
		}
		if (prog_state.verbose) {
			prog_state.err_stream << "OUT set to \"" << filepath <<"\"\n";
		}
		return true;
	}
	if (arg == "STD") {
		prog_state.set_output_stream(&std::cout);
		if (prog_state.verbose) {
			prog_state.err_stream << "OUT stopped.\n";
		}
		return true;
	}
	if (arg == "STOP") {
		prog_state.set_output_file(nullptr);
		if (prog_state.verbose) {
			prog_state.err_stream << "OUT stopped.\n";
		}
		return true;
	}
	return unrecognize_arg("OUT", command_vec.at(1).c_str(), prog_state);
}


void print_help(std::ostream & out) {
	out << PACKAGE_STRING << '\n';
	out << "\nUsage:\n	"<< PACKAGE << " [-chns] <filename-with-taxonomy-in-newick>\n";
	out << "Command line options:\n";
	out << "   -c <filename>   Specifies the path to a command file \n";
	out << "   -h              Prints this help message\n";
	out << "   -n ####         Gives a hint for the number of leaves\n";
	out << "   -s ####         Sets the initial random number seed\n";
	out << "   -x              Do not warn about duplicate labels\n";
	out << "\nAfter reading the taxonomy file commands will be read.\n";
	out << "  if a command file is not specified with the -c flag, then commands\n";
	out << "  will be read from standard input.\nThe commands are:\n";
	out << "\n   HELP     shows this message ;\n";
	out << "\n   OUTPUT [STD|STOP|FILE filename] ;\n";
	out << "             Specifies an output stream\n";
	out << "\n   PRINT ;    Writes the current tree to the output stream (in newick).\n";
	out << "\n   QUIT ;     Quits the program (surprise!)\n";
	out << "\n   REPEAT # ; CMD ; ENDREPEAT ;\n";
	out << "             Repeats the command(s) # times.\n";
	out << "\n   RESOLVE ; Randomly resolves all polytomies in the current tree.\n";
	out << "\n   SAMPLE RootMin = # RootMax = # InMin = # InMax = # OutMin = # OutMax = # ;\n";
	out << "             Causes the focal tree to be created by subsampling the full tree with\n";
	out << "               the specified root depth, ingroup and outgroup size. May fail if\n";
	out << "               relatively few leaves in the tree satisfy the constraints.\n";
	out << "\n   SET Seed = # Strict Verbose MaxTries = #;\n";
	out << "             Used to change program behaviors. \"Seed\" control the random number\n";
	out << "               generator. \"Strict\" causes the program to exit on any failed command.\n";
	out << "               \"Verbose\" produces more status updates in the standard error stream.\n";
	out << "               \"MaxTries\" controls the maximum number of times the SAMPLE command will\n";
	out << "               repeat its loop of trying to find a subsample that agrees with our constraints.\n";
	out << "\n   SPR Rep = # ReconMin = # ReconMax = # ;\n";
	out << "             Applies the specified number of random SPR edits to the focal tree\n";
	out << "               ReconMin and ReconMax limit the distance the moving subtree can travel.\n";
	out << "               The default ReconMin is 1. The movement of a subtree is similar to the.\n";
	out << "               extending TBR of MrBayes - a decision of which way to move is made at.\n";
	out << "               every fork. Note that the SPR will stop if the moving node \"runs out\"\n";
	out << "               of tree to traverse. Thus the movement can be smaller than reconmin; in\n";
	out << "               fact it is possible that a move will not change the topology.\n";
	out << "\n   STATUS ; Reports summary statistics about the current tree(s).\n";
	out << "\n   WEIGHT #.#  filename ;\n";
	out << "             Assigns tip weigth of #.# to every taxon listed in filename (a newline-\n";
	out << "               separated file). These weights are used to select the tips in the\n";
	out << "               SAMPLE command. For each draw, a leaf's probability of being chosen\n";
	out << "               is equal to its weight divided by the sum of all \"eligible\" leaves.\n";
	out << "               weights must be positive. The default weight is 1.0.\n";
	out << "\nCommands must be separated by semicolons !\n";
}

bool process_command(const ProgCommand & command_vec,
					 ProgState & prog_state) {
	if (command_vec.empty()) {
		return true;
	}

	//prog_state.err_stream << "Command:	";
	std::string cmd = capitalize(command_vec[0]);
	//prog_state.err_stream << cmd << "  \n";
	//for (const std::string & word : command_vec) {
	//		prog_state.err_stream << '"' << word << "\" ";
	//}
	//prog_state.err_stream << "\n";
	if (cmd == "QUIT" or cmd == "EXIT" or cmd == "Q") {
		return false;
	}
	if (cmd == "HELP") {
		print_help(std::cerr);
		return true;
	}
	else if (cmd == "ENDREPEAT") {
		return process_end_repeat_command(command_vec, prog_state);
	}
	else if (cmd == "OUT") {
		return process_out_command(command_vec, prog_state);
	}
	else if (cmd == "PRINT") {
		return process_print_command(command_vec, prog_state);
	}
	else if (cmd == "REPEAT") {
		return process_repeat_command(command_vec, prog_state);
	}
	else if (cmd == "RESOLVE") {
		return process_resolve_command(command_vec, prog_state);
	}
	else if (cmd == "SAMPLE") {
		return process_sample_command(command_vec, prog_state);
	}
	else if (cmd == "SET") {
		return process_set_command(command_vec, prog_state);
	}
	else if (cmd == "SPR") {
		return process_spr_command(command_vec, prog_state);
	}
	else if (cmd == "STATUS") {
		return process_status_command(command_vec, prog_state);
	}
	else if (cmd == "WEIGHT") {
		return process_weight_command(command_vec, prog_state);
	}
	else {
		prog_state.err_stream << PACKAGE << ": Command \"" << command_vec[0] << "\" is not recognized (use \"QUIT ;\" to exit)\n";
		return !prog_state.strict_mode;
	}
	return true;
}

bool rec_and_process_command(const ProgCommand & command_vec,
					 		 ProgState & prog_state,
					 		 bool record) {
	if (command_vec.empty()) {
		return true;
	}
	if (prog_state.curr_repeat_list) {
		//std::cerr << "prog_state.curr_repeat_list has " <<prog_state.curr_repeat_list->second.size() << " elements.\n";
	}
	else {
		//std::cerr << "prog_state.curr_repeat_list is NULL\n";
	}
	if (prog_state.command_list_list.empty()) {
		//std::cerr << "prog_state.command_list_list is empty\n";
	}
	else {
		//std::cerr << "prog_state.command_list_list has " <<prog_state.command_list_list.size() << " elements (first with " << prog_state.command_list_list.begin()->second.size()<< " commands)\n";
	}
	std::string cmd = capitalize(command_vec[0]);
	typedef std::list<ProgState::cmd_recorder_t>::iterator cmd_stack_it;
	if (record and cmd == "ENDREPEAT") {
		for (cmd_stack_it cll_it = prog_state.command_list_list.begin(); cll_it != prog_state.command_list_list.end(); ++cll_it) {
			ProgState::cmd_recorder_t & com_list = *cll_it;
			if (com_list.first) {
				com_list.second.push_back(command_vec);
			}
		}
	}
	const bool ret = process_command(command_vec, prog_state);
	if (cmd != "REPEAT" and prog_state.curr_repeat_list) {
		if (prog_state.curr_repeat_list->first) {
			prog_state.curr_repeat_list->second.push_back(command_vec);
		}
	}
	if (ret and record and cmd != "ENDREPEAT") {
		for (cmd_stack_it cll_it = prog_state.command_list_list.begin(); cll_it != prog_state.command_list_list.end(); ++cll_it) {
			ProgState::cmd_recorder_t & com_list = *cll_it;
			if (com_list.first) {
				com_list.second.push_back(command_vec);
			}
		}
	}
	return ret;
}

ProgCommand parse_command(std::istream & inp) {
	std::stringbuf str_buf;
	inp.get(str_buf, ';');
	std::string x = str_buf.str();
	if (inp.good()) {
		inp.get(); // skip the delimiting character
	}
	ProgCommand command_vec;
	if (x.empty()) {
		return command_vec;
	}
	bool in_word = false;
	unsigned word_start = 0;
	unsigned i = 0;
	for (; i < x.length(); ++i) {
		if (isgraph(x[i])) {
			if (std::strchr("(),=[]{}<>", x[i]) != nullptr) {
				if (in_word) {
					command_vec.push_back(x.substr(word_start, i - word_start));
					in_word = false;
				}
				command_vec.push_back(std::string(1, x[i]));
			}
			else {
				if (!in_word) {
					in_word = true;
					word_start = i;
				}
			}
		}
		else {
			if (in_word) {
				command_vec.push_back(x.substr(word_start, i - word_start));
			}
			in_word = false;
		}
	}
	if (in_word) {
		command_vec.push_back(x.substr(word_start, i - word_start));
	}
	return command_vec;
}
void run_command_interpreter(std::istream & command_stream,
							 ProgState & prog_state) {
	std::string command;
	bool keep_executing = true;
	while (keep_executing and command_stream.good()) {
		ProgCommand command_vec = parse_command(command_stream);
		keep_executing = rec_and_process_command(command_vec,
										 		 prog_state,
										 		 true);
	}
}

void warn_dup_label_callback(const std::string & original_label,
					const std::string & disambiguated_label,
					unsigned dup_num,
					void * /*nd_ptr*/) {
	std::cerr << "Duplicate label \"" << original_label << "\" instance #" << dup_num << " found. Label changed to \"" << disambiguated_label << "\"\n";
}
int main(int argc, char *argv[]) {
	std::ios_base::sync_with_stdio(false);
	if (argc < 2) {
		std::cerr << PACKAGE << ": Expecting a file path to a newick formatted tree\n";
		return 1;
	}
	std::string tree_filename;
	std::string cmd_filename;
	char prev_flag = '\0';
	char * endptr;
	RandGen::uint_seed_t seed = 0;
	bool warn_for_dup_label = true;
	for (int i = 1; i < argc; ++i) {
		if (prev_flag == 'n') {
			long int n = std::strtol(argv[i], &endptr, 10);
			if (endptr == argv[i]) {
				std::cerr << PACKAGE << ": Expecting number after -n\n";
				return 1;
			}
			if (n < 4) {
				std::cerr << PACKAGE << ": Expecting number of leaves > 4\n";
				return 1;
			}
			g_num_taxa_buckets = (unsigned long) 3*n;
			SimTree::set_initial_node_store_size(1.2*n);
			prev_flag = '\0';
		}
		else if (prev_flag == 's') {
			long int n = std::strtol(argv[i], &endptr, 10);
			if (endptr == argv[i]) {
				std::cerr << PACKAGE << ": Expecting number after -s\n";
				return 1;
			}
			if (n < 1) {
				std::cerr << PACKAGE << ": Expecting the random number seed to be > 1\n";
				return 1;
			}
			seed = (RandGen::uint_seed_t) n; //safe, positivity checked
			prev_flag = '\0';
		}
		else if (prev_flag == 'c') {
			cmd_filename = argv[i];
			if (cmd_filename.empty()) {
				std::cerr << PACKAGE << ": Expecting a filepath after the -c flag.\n";
				return 1;
			}
			prev_flag = '\0';
		}
		else if (prev_flag == '\0') {
			std::string arg(argv[i]);
			if (arg.length() > 1 and arg[0] == '-') {
				if (arg.length() > 2) {
					std::cerr << PACKAGE << ": Expecting each flag to be specified as a separate argument\n";
					return 1;
				}
				if (arg[1] == 'h') {
					print_help(std::cerr);
					return 0;
				}
				else if (arg[1] == 'x') {
					warn_for_dup_label = false;
				}
				else {
					prev_flag = arg[1];
					if (std::strchr("csn", prev_flag) == nullptr) {
						std::cerr << PACKAGE << ": Unknown flag -" << prev_flag << "\n";
						return 1;
					}
				}

			}
			else {
				if (!tree_filename.empty()) {
					std::cerr << PACKAGE << ": Expecting only one tree filename argument\n";
					return 1;
				}
				tree_filename.assign(argv[i]);
			}
		}
	}

	if (tree_filename.empty()) {
		std::cerr << PACKAGE << ": Expecting only a tree filename argument\n";
		return 1;
	}
	std::ifstream inp(tree_filename.c_str());
	if (!inp.good()) {
		std::cerr << PACKAGE << ": Could not open " << tree_filename << "\n";
		return 2;
	}

	std::ifstream cmd_inp;
	if (!cmd_filename.empty()) {
		cmd_inp.open(cmd_filename);
		if (!cmd_inp.good()) {
			std::cerr << PACKAGE << ": Could not open " << cmd_filename << "\n";
			return 2;
		}
	}

	TaxonNameUniverse taxa;
	SimTree * tree = nullptr;
	SimNdBlob nd_blob;
	SimTreeBlob tree_blob;
	try {
		dup_label_callback_ptr_t w_ptr = (warn_for_dup_label ? warn_dup_label_callback : nullptr);
		tree = parse_from_newick_stream<SimNdBlob, SimTreeBlob>(inp,
																taxa,
																nd_blob,
																tree_blob,
																w_ptr);
		if (tree == nullptr) {
			std::cerr << PACKAGE << ": No tree found!\n";
			return 1;
		}
		std::cerr << taxa.get_num_labels() << " labels read.\n";
	}
	catch (ParseExcept & x) {
		std::cerr << PACKAGE << " Error:	" << x.message << "\nAt line = " << x.fileline << " at column = " << x.filecol << " at pos = " << x.filepos << "\n";
		return 3;
	}

	ProgState prog_state(taxa, *tree, seed);
	prog_state.nd_blob = nd_blob;
	prog_state.tree_blob = tree_blob;
	prog_state.set_output_stream(&std::cout);
	std::istream & command_stream = std::cin;
	if (cmd_filename.empty()) {
		run_command_interpreter(command_stream, prog_state);
	}
	else {
		run_command_interpreter(cmd_inp, prog_state);
	}
	prog_state.set_output_file(nullptr);

	return 0;
}

////////////////////////////////////////////////////////////////////////////////
// Written by Mark T. Holder. 2012. Available for use under the terms of either
//  the BSD or GPL (v3) license. (see BSDLicense.txt GPL.txt). Take your pick.
// No WARRANTY!
////////////////////////////////////////////////////////////////////////////////
