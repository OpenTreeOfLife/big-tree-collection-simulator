#include <iostream>
#include <fstream>
#include <cstring>
#include <string>

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

void cleantree(std::istream & inp, std::ostream & out) {
    enum READING_MODE {IN_LABEL, OUT_OF_LABEL, IN_QUOTE};
    READING_MODE curr_mode = OUT_OF_LABEL;
    std::string label;
    std::string cache;
    bool needsQuoting = false;
    long filepos = 0;
    long filecol = 0;
    long fileline = 1;
    while (inp.good()) {
        char c = inp.get(); filepos++; filecol++;
        if (std::isgraph(c)) {
            if (strchr("(),:;", c) != NULL) {
                if (curr_mode == IN_LABEL) {
                    curr_mode = OUT_OF_LABEL;
                    if (needsQuoting) {
                        out << '\'' << label << '\'';
                    }
                    else {
                        out << label;
                    }
                    needsQuoting = false;
                }
                out << c;
                if (c == ';')
                    out << "\n";
            }
            else {
                if (curr_mode == OUT_OF_LABEL) {
                    label.clear();
                    needsQuoting = false;
                    curr_mode = IN_LABEL;
                    cache.clear();
                    if (c == '\'') {
                        if (!inp.good()) {
                            throw ParseExcept("Quote started then EOF", fileline, filecol, filepos);
                        }
                        char q = inp.get(); filepos++; filecol++;
                        for (;;) {
                            if (!inp.good()) {
                                throw ParseExcept("File reading before termination of quote", fileline, filecol, filepos);
                            }
                            label.append(1, q);
                            if (q == '\'') {
                                const char d = inp.peek();
                                if (d == '\'') {
                                    label.append(1, q);
                                    q = inp.get(); filepos++; filecol++;
                                }
                            }
                            else {
                                break;
                            }
                        }
                    }
                }
                else if (!cache.empty()) {
                    label += cache;
                    cache.clear();
                    needsQuoting = true;
                }
                label.append(1, c);
                if (strchr("\'[(){}\"-]/\\,;:=*`+<>", c) != NULL) {
                    needsQuoting = true;
                    if (c == '\'') {
                        label.append(1, '\'');
                    }
                }
            }
            //std::cerr << c << '\n';
        }
        else if (curr_mode == IN_LABEL) {
            if (c == '\n') {
                filecol = 0;
                fileline = 1;
            }
            cache.append(1, ' ');
        }
    }
}
int main(int argc, char *argv[]) {
    if (argc < 2) {
        std::cerr << "Expecting a file path to a newick formatted tree\n";
        return 1;
    }
    std::string filename(argv[1]);
    std::ifstream inp(filename.c_str());
    if (!inp.good()) {
        std::cerr << "Could not open " << filename << "\n";
        return 2;
    }
    std::ostream & out(std::cout);
    try {
        cleantree(inp, out);
        out << '\n';
    }
    catch (ParseExcept & x) {
        std::cerr << "Error:  " << x.message << "\nAt line = " << x.fileline << " at column = " << x.filecol << " at pos = " << x.filepos << "\n";
        return 3;
    }
    return 0;
}
