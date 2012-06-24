

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <functional>
#include <fstream>
#include <ios>
#include <random>

#include <config.h>
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

	tree.debug_check();
	for (; it != tree.end_fast_postorder(); ++it) {
		Node<T> & nd = *it;
		if (nd.is_polytomy()) {
			tree.debug_check();
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
				tree.debug_check();
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
			_sum_leaf_weights(1.0) {
		}

		double get_sum_leaf_weights() {
			return this->_sum_leaf_weights;
		}
		void set_sum_leaf_weights(double x) {
			this->_sum_leaf_weights = x;
		}
		double get_edge_length() const {
			return _edge_length;
		}

	private:
		double _edge_length;
		double _sum_leaf_weights;
};
class SimTreeBlob {
	public:
		SimTreeBlob()
			:sum_leaf_weights(-1.0) {
		}

		double get_sum_leaf_weights() {
			return this->sum_leaf_weights;
		}
		void set_sum_leaf_weights(double x) {
			this->sum_leaf_weights = x;
		}
	private:
		double sum_leaf_weights;
};
typedef Node<SimNdBlob> SimNode;
typedef Tree<SimNdBlob, SimTreeBlob> SimTree;

////////////////////////////////////////////////////////////////////////////////
// Main impl
////////////////////
template <typename T, typename U>
void sum_leaf_weights_over_tree(Tree<T,U> & tree) {
	typedef typename Tree<T,U>::Node_T::fast_postorder_iterator nd_it_t;
	typedef typename Node<T>::const_child_iterator const_child_it;
	for (nd_it_t it = tree.begin_fast_postorder(); it != tree.end_fast_postorder(); ++it) {
		Node<T> & nd = *it;
		if (nd.is_internal()) {
			double x = 0.0;
			for (const_child_it c_it = nd.begin_child(); c_it != nd.end_child(); ++c_it) {
				x += c_it->blob.get_sum_leaf_weights();
			}
			nd.blob.set_sum_leaf_weights(x);
		}
	}
	assert(tree.get_root());
	tree.blob.set_sum_leaf_weights(tree.get_root()->blob.get_sum_leaf_weights());
}


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
			outp(&std::cout),
			full_tree(tree),
			current_tree(nullptr),
			last_seed(seed),
			rng(seed) {
			this->current_tree = &(this->full_tree);
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
				this->err_stream << "Calling write_newick...\n";
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
		SimTree * get_focal_tree() const {
			return this->current_tree;
		}
		TaxonNameUniverse & taxa;
	private:
		std::ostream * outp;
		std::ofstream outp_obj;
		SimTree & full_tree;
		SimTree * current_tree;
		RandGen::uint_seed_t last_seed;
	public:
		RandGen rng;

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

bool process_resolve_command(const std::vector<std::string> ,
						   ProgState & prog_state) {
	SimTree * focal_tree = prog_state.get_focal_tree();
	assert(focal_tree);
	resolve_polytomies(*focal_tree, prog_state.rng);
	return true;
}

bool process_weight_command(const std::vector<std::string> command_vec,
						   ProgState & prog_state) {
	if (command_vec.size() <= 2) {
		std::cerr << "Expected a weight then a filepath of names after the WEIGHT command.\n";
		return prog_state.strict_mode;
	}
	std::string weight_string = command_vec[1];
	char * e;
	double wt = std::strtod(weight_string.c_str(), &e);
	if (e != weight_string.c_str() + weight_string.length()) {
		prog_state.err_stream << "Expected a number (the weight) after the WEIGHT command. found " << weight_string << ".\n";
		return prog_state.strict_mode;
	}

	std::string fn = command_vec[2];
	std::ifstream inp(fn);
	if (!inp.good()) {
		prog_state.err_stream << "Could not open the file " << fn << ".\n";
		return prog_state.strict_mode;
	}
	const std::vector<TaxonLabel> labels = parse_labels_from_stream(inp, prog_state.taxa);
	assert(prog_state.get_focal_tree());
	SimTree & tree = *(prog_state.get_focal_tree());
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
	return true;
}

bool process_sample_command(const std::vector<std::string> command_vec,
						   ProgState & prog_state) {
	SimTree * focal_tree = prog_state.get_focal_tree();
	assert(focal_tree);
	if (focal_tree->blob.get_sum_leaf_weights() < 0.0) {
		sum_leaf_weights_over_tree(*focal_tree);
	}
	prog_state.err_stream << "focal_tree->blob.get_sum_leaf_weights() = " << focal_tree->blob.get_sum_leaf_weights() << '\n';
	return true;
}

bool process_print_command(const std::vector<std::string> command_vec,
						   ProgState & prog_state) {
	bool nexus = false;
	if (command_vec.size() > 1 && capitalize(command_vec[1]) == "NEXUS") {
		nexus = true;
	}
	prog_state.print_tree(false, nexus);
	return true;
}
bool process_out_command(const std::vector<std::string> command_vec,
						 ProgState & prog_state) {
	if (command_vec.size() == 1) {
		prog_state.err_stream << "Expecting an argument (FILE, STD, or STOP) after OUT command.";
		return !prog_state.strict_mode;
	}
	std::string arg = capitalize(command_vec[1]);
	if (arg == "FILE") {
		if (command_vec.size() == 2) {
			prog_state.err_stream << "Expecting a file path after OUT FILE";
			return !prog_state.strict_mode;
		}
		const char * filepath = command_vec.at(2).c_str();
		prog_state.set_output_file(filepath);
		if (! prog_state.out_good()) {
			prog_state.err_stream << "Could not open the file \"" << filepath << "\". Output set to standard out!\n";
			prog_state.set_output_stream(&std::cout);
			return !prog_state.strict_mode;
		}
		return true;
	}
	if (arg == "STD") {
		prog_state.set_output_stream(&std::cout);
		return true;
	}
	if (arg == "STOP") {
		prog_state.set_output_file(nullptr);
		return true;
	}
	return unrecognize_arg("OUT", command_vec.at(1).c_str(), prog_state);
}

bool process_command(const std::vector<std::string> & command_vec,
					 const TaxonNameUniverse & ,
					 const SimTree & ,
					 ProgState & prog_state) {
	if (command_vec.empty()) {
		return true;
	}

	prog_state.err_stream << "Command:	";
	std::string cmd = capitalize(command_vec[0]);
	prog_state.err_stream << cmd << "  \n";
	for (const std::string & word : command_vec) {
		prog_state.err_stream << '"' << word << "\" ";
	}
	prog_state.err_stream << "\n";
	if (cmd == "QUIT") {
		return false;
	}
	else if (cmd == "OUT") {
		return process_out_command(command_vec, prog_state);
	}
	else if (cmd == "PRINT") {
		return process_print_command(command_vec, prog_state);
	}
	else if (cmd == "RESOLVE") {
		return process_resolve_command(command_vec, prog_state);
	}
	else if (cmd == "SAMPLE") {
		return process_sample_command(command_vec, prog_state);
	}
	else if (cmd == "WEIGHT") {
		return process_weight_command(command_vec, prog_state);
	}
	else {
		prog_state.err_stream << "Command \"" << command_vec[0] << "\" is not recognized (use \"QUIT ;\" to exit)\n";
		return !prog_state.strict_mode;
	}
	return true;

}

std::vector<std::string> parse_command(std::istream & inp) {
	std::stringbuf str_buf;
	inp.get(str_buf, ';');
	std::string x = str_buf.str();
	if (inp.good()) {
		inp.get(); // skip the delimiting character
	}
	std::vector<std::string> command_vec;
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
							 const TaxonNameUniverse & taxa,
							 const SimTree & tree,
							 ProgState & prog_state) {
	std::string command;
	bool keep_executing = true;
	while (keep_executing && command_stream.good()) {
		std::vector<std::string> command_vec = parse_command(command_stream);
		keep_executing = process_command(command_vec,
										 taxa,
										 tree,
										 prog_state);
	}
}
void print_help(std::ostream & out) {
	out << PACKAGE_STRING << '\n';
	out << "\nUsage:\n	"<< PACKAGE << " [-ns] <taxonomy-file>\n";
	out << "Command line options:\n";
	out << "   -h		  this help message\n";
	out << "   -i		  interactive mode\n";
	out << "   -n ####	  number of leaves\n";
	out << "   -s ####	  random number seed\n";
	out << "\n\nIf you enter interactive mode, then commands will be read from standard input.\nThe commands are:\n";
	out << "   OUTPUT [STD|STOP|FILE fn]   Specifies an output stream\n";
	out << "   PRINT	 Writes the current tree to the output stream (in newick).\n";
	out << "   RESOLVE	 Randomly resolves all polytomies in the current tree.\n";
	out << "   QUIT		 Quits the program (surprise!)\n";
	out << "\nCommands must be separated by semicolons !\n";
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
	bool interactive = false;
	RandGen::uint_seed_t seed = 0;
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
			g_num_taxa_buckets = (unsigned long) 3*n;
			SimTree::set_initial_node_store_size(1.2*n);
			prev_flag = '\0';
		}
		else if (prev_flag == 's') {
			long int n = std::strtol(argv[i], &endptr, 10);
			if (endptr == argv[i]) {
				std::cerr << "Expecting number after -s\n";
				return 1;
			}
			if (n < 1) {
				std::cerr << "Expecting the random number seed to be > 1\n";
				return 1;
			}
			seed = (RandGen::uint_seed_t) n; //safe, positivity checked
			prev_flag = '\0';
		}
		else if (prev_flag == '\0') {
			std::string arg(argv[i]);
			if (arg.length() > 1 and arg[0] == '-') {
				if (arg.length() > 2) {
					std::cerr << "Expecting each flag to be specified as a separate argument\n";
					return 1;
				}
				if (arg[1] == 'h') {
					print_help(std::cerr);
					return 0;
				}
				else if (arg[1] == 'i') {
					interactive = true;
				}
				else {
					prev_flag = arg[1];
				}
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
	TaxonNameUniverse taxa;
	SimTree * tree = nullptr;
	try {
		SimNdBlob nd_blob;
		SimTreeBlob tree_blob;
		tree = parse_from_newick_stream<SimNdBlob, SimTreeBlob>(inp, taxa, nd_blob, tree_blob);
		if (tree == nullptr) {
			std::cerr << "No tree found!\n";
			return 1;
		}
		std::cerr << taxa.get_num_labels() << " labels read.\n";
	}
	catch (ParseExcept & x) {
		std::cerr << "\nError:	" << x.message << "\nAt line = " << x.fileline << " at column = " << x.filecol << " at pos = " << x.filepos << "\n";
		return 3;
	}

	ProgState prog_state(taxa, *tree, seed);
	prog_state.set_output_stream(&std::cout);
	std::istream & command_stream = std::cin;
	if (interactive) {
		run_command_interpreter(command_stream, taxa, *tree, prog_state);
	}
	prog_state.set_output_file(nullptr);

	return 0;
}

