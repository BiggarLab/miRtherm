// Command line front end.
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "../core/mirtherm.h"

using namespace mirtherm;

static void usage() {
  std::cerr <<
      "Usage: \n\n"
      "mirtherm mirnafile utrfile [-o outfile] [-tab|-text]\n"
      "     [-at at_weight] [-gc gc_weight] [-gt gt_weight]\n"
      "     [-mis mis_weight] [-go gap_open] [-ge gap_elong]\n"
      "     [-energy energy] [-score score] [-loop loopscore]\n"
      "     [-T temperature] [-noLP] [-noGU] [-noCloseGU]\n\n"
      "   mirnafile and utrfile are miRNA and UTR sequences in FASTA format.\n"
      "   -energy  keep sites with duplex energy <= energy (kcal/mol)\n"
      "   -score   keep sites with alignment score >= score\n"
      "   -loop    keep sites with loop score >= loopscore\n"
      "   -T       folding temperature in degrees C (default 37)\n";
  std::exit(1);
}

static void bad(const char *a) {
  std::cerr << "Unknown or incomplete option: " << a << "\n";
  usage();
}

int main(int argc, char **argv) {
  Params p;
  std::vector<std::string> files;
  std::string outfile;
  for (int i = 1; i < argc; ++i) {
    std::string a = argv[i];
    auto need = [&]() -> const char * {
      if (i + 1 >= argc) bad(argv[i]);
      return argv[++i];
    };
    if (a.empty() || a[0] != '-') files.push_back(a);
    else if (a == "-h" || a == "--help") usage();
    else if (a == "--version") { std::cerr << "mirtherm :: " << kVersion << std::endl; return 1; }
    else if (a == "-o") outfile = need();
    else if (a == "-tab") p.format = 'g';
    else if (a == "-text") p.format = 't';
    else if (a == "-at") p.at_weight = (int)std::strtol(need(), nullptr, 10);
    else if (a == "-gc") p.gc_weight = (int)std::strtol(need(), nullptr, 10);
    else if (a == "-gt") p.gt_weight = (int)std::strtol(need(), nullptr, 10);
    else if (a == "-mis") p.mis_weight = (int)std::strtol(need(), nullptr, 10);
    else if (a == "-go") p.gap_open = (int)std::strtol(need(), nullptr, 10);
    else if (a == "-ge") p.gap_ext = (int)std::strtol(need(), nullptr, 10);
    else if (a == "-se") p.seed_end = (int)std::strtol(need(), nullptr, 10);
    else if (a == "-energy") p.filter_energy = (float)std::strtod(need(), nullptr);
    else if (a == "-score") p.filter_score = (int)std::strtol(need(), nullptr, 10);
    else if (a == "-loop") p.filter_loop = (float)std::strtod(need(), nullptr);
    else if (a == "-T") p.temperature = std::strtod(need(), nullptr);
    else if (a == "-noLP") p.no_lonely_pairs = true;
    else if (a == "-noGU") p.no_gu = true;
    else if (a == "-noCloseGU") p.no_closing_gu = true;
    else bad(argv[i]);
  }
  if (files.size() != 2) usage();

  std::vector<Sequence> mirs, utrs;
  {
    std::ifstream f(files[0].c_str());
    if (!f) { std::cerr << "Cannot open " << files[0] << "\n"; return 1; }
    read_fasta(f, mirs);
  }
  {
    std::ifstream f(files[1].c_str());
    if (!f) { std::cerr << "Cannot open " << files[1] << "\n"; return 1; }
    read_fasta(f, utrs);
  }
  if (mirs.empty() || utrs.empty()) { std::cerr << "No input sequences!\n"; return 1; }

  std::ofstream of;
  std::ostream *os = &std::cout;
  if (!outfile.empty()) {
    std::cerr << "Output results to " << outfile << std::endl;
    of.open(outfile.c_str());
    if (!of) { std::cerr << "Cannot create output file.\n"; return 1; }
    os = &of;
  }
  if (!run(mirs, utrs, p, *os)) {
    os->flush();
    std::cerr << "Not a nucleotide!" << std::endl;
    return 1;
  }
  if (!outfile.empty()) std::cerr << "Store successfully." << std::endl;
  return 0;
}
