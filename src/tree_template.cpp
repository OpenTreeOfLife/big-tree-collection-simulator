#include <unordered_map>
#include <vector>
#include <fstream>
#include <iostream>
#include "tree_template.hpp"

nnodes_t g_num_taxa_buckets = 100;
nnodes_t g_initial_node_store_size = 100;

std::vector<TaxonLabel> parse_labels_from_stream(std::istream & input,
												 TaxonNameUniverse & taxa) {
	enum READING_MODE {IN_LABEL, OUT_OF_LABEL, IN_QUOTE};
	READING_MODE curr_mode = OUT_OF_LABEL;
	std::string label;
	std::string cache;
	long filepos = 0;
	long filecol = 0;
	long fileline = 1;
	int comment_level = 0;
	const bool preserve_comments = false;
	std::streambuf * rdbuf = input.rdbuf();
	std::vector<TaxonLabel> vec;
	for(;;) {
		signed char c = rdbuf->sbumpc(); filepos++; filecol++;
		if (c == ',' or c == '\n') {
			if (curr_mode == IN_LABEL) {
				const TaxonLabel * tl = taxa.find_name(label);
				if (tl == nullptr) {
					std::cerr << "label \"" << label << "\" not found. Ignored...\n";
				}
				vec.push_back(*tl);
			}
			if (c == '\n') {
				filecol = 0;
				fileline = 1;
			}
			curr_mode = OUT_OF_LABEL;
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
					const TaxonLabel * tl = taxa.find_name(label);
					if (tl == nullptr) {
						std::cerr << "label \"" << label << "\" not found. Ignored...\n";
					}
					vec.push_back(*tl);
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
				return vec;
			}
			if (curr_mode == IN_LABEL) {
				cache.append(1, ' '); // converts all non-graphical characters to spaces!
			}
		}
	}
	return vec;
}
