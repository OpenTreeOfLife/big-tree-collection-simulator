#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "tree_template.hpp"

int main(int argc, char * argv[]) {
    if (argc != 2) {
        std::cerr << "Expecting exactly one argument - a filepath to a newick file\n";
        return 1;
    }
    std::string fn = argv[1];
    std::ifstream inp(fn);
    SlimTaxonNameUniverse taxa;
    const SlimTree * tree = parse_from_newick_stream(inp, taxa);
    if (tree == nullptr) {
        std::cerr << "No tree found in " << fn << "\n";
        return 1;
    }
    std::vector<const SlimNode *> pre_vec, post_vec, leaf_vec;
    SlimNode::const_preorder_iterator it = tree->begin_preorder();
    for (; it != tree->end_preorder(); ++it) {
        pre_vec.push_back(&(*it));
    }
    SlimNode::const_postorder_iterator po_it = tree->begin_postorder();
    for (; po_it != tree->end_postorder(); ++po_it) {
        post_vec.push_back(&(*po_it));
    }
    SlimNode::const_leaf_iterator l_it = tree->begin_leaf();
    for (; l_it != tree->end_leaf(); ++l_it) {
        leaf_vec.push_back(&(*l_it));
    }
    if (post_vec.size() != pre_vec.size()) {
        std::cerr << "post_vec.size() == " << post_vec.size() << "\npre_vec.size() == " << pre_vec.size() << '\n';
        return 1;
    }
    std::vector<const SlimNode *>::const_iterator pre_it = pre_vec.begin();
    std::vector<const SlimNode *>::const_iterator leaf_it = leaf_vec.begin();
    std::vector<const SlimNode *>::const_reverse_iterator pre_rev_it = post_vec.rbegin();
    for (; pre_it != pre_vec.end(); ++pre_rev_it, ++pre_it) {
        const SlimNode * p = *pre_it;
        if (p->is_leaf()) {
            const SlimNode * pl = *leaf_it++;
            if (pl != p) {
                std::cerr << "leaf != preorder leaf\n";
                return 1;
            }
        }
        const SlimNode * q = *pre_rev_it;
        if (p != q) {
            std::cerr << "p != q\n";
            return 1;
        }
    }
    return 0;
}
