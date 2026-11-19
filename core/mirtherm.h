// miRtherm: microRNA target prediction with temperature-adjusted duplex energy.
//
// Pipeline, per miRNA / UTR pair:
//   1. seed scan      6-mer seed matches for miRNA nt 1-6, 2-7 and 3-8, allowing
//                     one G:U wobble
//   2. gapped align   affine-gap alignment of the miRNA 3' part against the
//                     upstream target window (max 40 nt)
//   3. duplex energy  ViennaRNA 1.5 cofold (Turner 1999 parameters) at the
//                     requested temperature, extending the target if the miRNA
//                     3' end overhangs
//   4. loop score     position/size score of the central loop
//   5. filters        energy, alignment score, loop score
#ifndef MIRTHERM_H
#define MIRTHERM_H

#include <ostream>
#include <string>
#include <vector>

namespace mirtherm {

extern const char *const kVersion;  // version string printed in the output header

struct Params {
  int at_weight = 5;
  int gc_weight = 5;
  int gt_weight = 2;
  int mis_weight = -3;
  int gap_open = -8;
  int gap_ext = -2;
  int seed_end = 8;
  float filter_energy = 1000;  // keep sites with energy <= filter_energy
  int filter_score = -1000;    // keep sites with score  >= filter_score
  float filter_loop = -1000;   // keep sites with loop   >= filter_loop
  double temperature = 37;     // degrees C
  bool no_lonely_pairs = false;
  bool no_gu = false;
  bool no_closing_gu = false;
  char format = 'g';           // 'g' tab-delimited (default), 't' text
};

struct Sequence {
  std::string id;
  std::string seq;
};

struct Align {
  std::string ds;   // miRNA row (3'->5')
  std::string aln;  // match row
  std::string cs;   // target row (5'->3')
};

struct Result {
  std::string utrid, mirid;
  Align dynaln;  // dynamic-programming alignment
  Align eryaln;  // alignment derived from the cofold structure
  std::string structure;
  int pos[2] = {1, 1};
  int score = 0;
  float seedery = 0;
  int seedscore = 0;
  float energy = 0;
  float loopscore = 0;
  int loop[2] = {0, 0};
  float tr_score = 0;
  float ol_score = 0;
};

class Scanner {
 public:
  explicit Scanner(const Params &p);

  // Scan one UTR with one miRNA. Every site passing the filters is passed to
  // `emit`. Returns the number of candidate sites examined, or -1 on
  // a non-nucleotide character or a miRNA that is too short or too long.
  template <class F>
  int scan(const Sequence &utr, const Sequence &mir, F emit) {
    return scan_impl(utr, mir, &thunk<F>, &emit);
  }

  // Write one result (tab or text format).
  void write(std::ostream &os, const Result &r) const;
  // Header block ("# version ..." etc.).
  void write_header(std::ostream &os) const;

  const Params &params() const { return p_; }

 private:
  template <class F>
  static void thunk(void *ctx, const Result &r) { (*static_cast<F *>(ctx))(r); }
  int scan_impl(const Sequence &utr, const Sequence &mir,
                void (*emit)(void *, const Result &), void *ctx);
  bool get_result(const std::string &utrsub, const std::string &mir, Result &r);
  int align(const char *a, const char *b, Align &out) const;
  float get_energy(Result &r, const std::string &utrsub, const std::string &mir);
  void score_loop(Result &r) const;

  Params p_;
  int matrix_[5][5];
};

// FASTA reader (ids keep the full header line,
// sequences are upper-cased, text after a '\r' is dropped).
void read_fasta(std::istream &in, std::vector<Sequence> &out);

// Whole run, written as the command line tool prints it. Returns
// false if it stopped on a non-nucleotide character.
bool run(const std::vector<Sequence> &mirs, const std::vector<Sequence> &utrs,
         const Params &p, std::ostream &os);

}  // namespace mirtherm

#endif
