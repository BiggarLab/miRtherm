// WebAssembly entry points used by the browser dashboard.
#include <sstream>
#include <string>

#include "../core/mirtherm.h"

using namespace mirtherm;

static Scanner *g_scanner = nullptr;
static std::string g_out;

extern "C" {

void mt_configure(int at, int gc, int gt, int mis, int go, int ge, double temp, double f_energy,
                  int f_score, double f_loop, int no_lp, int no_gu, int no_close_gu) {
  Params p;
  p.at_weight = at;
  p.gc_weight = gc;
  p.gt_weight = gt;
  p.mis_weight = mis;
  p.gap_open = go;
  p.gap_ext = ge;
  p.temperature = temp;
  p.filter_energy = (float)f_energy;
  p.filter_score = f_score;
  p.filter_loop = (float)f_loop;
  p.no_lonely_pairs = no_lp;
  p.no_gu = no_gu;
  p.no_closing_gu = no_close_gu;
  delete g_scanner;
  g_scanner = new Scanner(p);
}

// Header block of the output file for the current configuration.
const char *mt_header() {
  std::ostringstream os;
  if (g_scanner) g_scanner->write_header(os);
  g_out = os.str();
  return g_out.c_str();
}

// Scan one UTR with one miRNA. Returns the result lines (tab format,
// one per site) or NULL if a sequence holds a non-nucleotide where
// the command line program stops.
const char *mt_scan(const char *mir_id, const char *mir_seq, const char *utr_seq) {
  if (!g_scanner) return nullptr;
  Sequence mir{mir_id, mir_seq}, utr{"", utr_seq};
  std::ostringstream os;
  int n = g_scanner->scan(utr, mir, [&](const Result &r) { g_scanner->write(os, r); });
  if (n < 0) return nullptr;
  g_out = os.str();
  return g_out.c_str();
}

}  // extern "C"
