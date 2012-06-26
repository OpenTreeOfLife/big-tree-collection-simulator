#include <iostream>
#include <fstream>
#include <string>
#include <set>
#include <vector>
#include "tree_template.hpp"

int main(int argc, char * argv[]) {
    if (argc != 2) {
        std::cerr << "Expecting exactly one argument - a filepath to a newick file\n";
        return 1;
    }
    std::string fn = argv[1];
    std::ifstream inp(fn);
    TaxonNameUniverse taxa;
    empty_node_blob_t nd_blob;
    empty_tree_blob_t tree_blob;
    const SlimTree * tree = parse_from_newick_stream(inp, taxa, nd_blob, tree_blob, nullptr);
    if (tree == nullptr) {
        std::cerr << "No tree found in " << fn << "\n";
        return 1;
    }
    std::vector<const SlimNode *> pre_vec, post_vec, leaf_vec, fast_post_vec;
    SlimNode::const_preorder_iterator it = tree->begin_preorder();
    for (; it != tree->end_preorder(); ++it) {
        pre_vec.push_back(&(*it));
    }
    SlimNode::const_postorder_iterator po_it = tree->begin_postorder();
    for (; po_it != tree->end_postorder(); ++po_it) {
        post_vec.push_back(&(*po_it));
    }
    SlimNode::const_fast_postorder_iterator fast_post_it = tree->begin_fast_postorder();
    for (; fast_post_it != tree->end_fast_postorder(); ++fast_post_it) {
        fast_post_vec.push_back(&(*fast_post_it));
    }
    SlimNode::const_leaf_iterator l_it = tree->begin_leaf();
    for (; l_it != tree->end_leaf(); ++l_it) {
        leaf_vec.push_back(&(*l_it));
    }
    if (post_vec.size() != pre_vec.size()) {
        std::cerr << "post_vec.size() == " << post_vec.size() << "\npre_vec.size() == " << pre_vec.size() << '\n';
        return 1;
    }
    if (post_vec.size() != fast_post_vec.size()) {
        std::cerr << "post_vec.size() == " << post_vec.size() << "\nfast_post_vec.size() == " << pre_vec.size() << '\n';
        return 1;
    }
    std::set<const SlimNode *> left_seen;
    std::vector<const SlimNode *>::const_iterator pre_it = pre_vec.begin();
    std::vector<const SlimNode *>::const_iterator leaf_it = leaf_vec.begin();
    std::vector<const SlimNode *>::const_iterator fast_it = fast_post_vec.begin();
    std::vector<const SlimNode *>::const_reverse_iterator pre_rev_it = post_vec.rbegin();
    for (; pre_it != pre_vec.end(); ++pre_rev_it, ++pre_it, ++fast_it) {
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
        const SlimNode * n = *fast_it;
        if (left_seen.find(n) != left_seen.end()) {
            std::cerr << "repeated node in fast_postorder\n";
            return 1;
        }
        left_seen.insert(n);
        if (left_seen.find(n->get_parent()) != left_seen.end()) {
            std::cerr << "parent before child in fast_postorder\n";
            return 1;
        }
    }
    return 0;
}

////////////////////////////////////////////////////////////////////////////////
// Written by Mark T. Holder. 2012. Available for use under the terms of either
//  the BSD or GPL (v3) license. (see BSDLicense.txt GPL.txt). Take your pick.
// No WARRANTY!
////////////////////////////////////////////////////////////////////////////////
