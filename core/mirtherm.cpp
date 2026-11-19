// miRtherm core. See mirtherm.h for an overview.
//
// A few steps look odd on purpose (noted where they occur). They are kept so
// that output stays identical to the reference outputs used for testing.
#include "mirtherm.h"

#include <algorithm>
#include <cstring>
#include <istream>

extern "C" {
// Bundled ViennaRNA 1.5 (see core/vienna).
float cofold(const char *sequence, char *structure);
void update_cofold_params(void);
extern int cut_point;
extern double temperature;
extern int noLonelyPairs;
extern int noGU;
extern int no_closingGU;
}

namespace mirtherm {

const char *const kVersion = "3.11.12";

namespace {

const float kTypeScore[9] = {0.5f, 0.75f, 1, 1, 1, 1, 0.75f, 0.75f, 0.5f};
const int kLenScore[5] = {10, 20, 25, 25, 20};
const int kWindow = 40;   // target nucleotides upstream of the miRNA 5' end
const int kMinWindow = 10;
const int kMaxAln = 50;   // DP matrix dimension

inline void reverse(std::string &s) { std::reverse(s.begin(), s.end()); }

// Complement; returns false on a non-nucleotide character.
bool comp_seq(const std::string &in, std::string &out) {
  out.clear();
  for (char ch : in) {
    switch (ch) {
      case 'A': out += 'T'; break;
      case 'C': out += 'G'; break;
      case 'G': out += 'C'; break;
      case 'T':
      case 'U': out += 'A'; break;
      default: return false;
    }
  }
  return true;
}

// Target-side patterns for a 6-nt miRNA seed: the reverse complement, plus
// every single G:U wobble variant (target C->T and target A->G).
bool get_pattern(const std::string &seed, std::string pat[7]) {
  std::string rc;
  if (!comp_seq(seed, rc)) return false;
  reverse(rc);
  pat[0] = rc;
  for (size_t i = 0; i < 6 && i < rc.size(); ++i) {
    if (rc[i] == 'C') {
      pat[i + 1] = rc;
      pat[i + 1][i] = 'T';
    } else if (rc[i] == 'A') {
      pat[i + 1] = rc;
      pat[i + 1][i] = 'G';
    } else {
      pat[i + 1].clear();
    }
  }
  return true;
}

// Match row and score for two aligned rows. Uses fixed 5 / 2 / -3 scores,
// independent of the alignment weights.
int get_aln(const std::string &a, const std::string &b, std::string &out) {
  out.clear();
  int score = 0;
  for (size_t i = 0; i < a.size(); ++i) {
    switch (a[i] + b[i]) {
      case 'C' + 'G':
      case 'A' + 'T':
      case 'A' + 'U':
        out += '|';
        score += 5;
        break;
      case 'G' + 'T':
      case 'G' + 'U':
        out += ':';
        score += 2;
        break;
      default:
        out += '*';
        score -= 3;
    }
  }
  return score;
}

// Left-pad the shorter string with '-'.
void costring(std::string &a, std::string &b) {
  if (a.size() > b.size())
    b.insert(0, a.size() - b.size(), '-');
  else if (b.size() > a.size())
    a.insert(0, b.size() - a.size(), '-');
}

// Index of the first character differing from s[0] if s starts with c;
// 0 if s does not start with c; npos if s is empty or uniform.
size_t leading_run(const std::string &s, char c) {
  if (s.empty()) return std::string::npos;
  if (s[0] != c) return 0;
  for (size_t i = 1; i < s.size(); ++i)
    if (s[i] != s[0]) return i;
  return std::string::npos;
}

// 1-based pair table from a dot-bracket string; 0 = unpaired.
void pair_table(const char *structure, std::vector<int> &pt) {
  size_t n = std::strlen(structure);
  pt.assign(n + 2, 0);
  std::vector<int> stack;
  for (size_t i = 1; i <= n; ++i) {
    char ch = structure[i - 1];
    if (ch == '(') {
      stack.push_back((int)i);
    } else if (ch == ')' && !stack.empty()) {
      int j = stack.back();
      stack.pop_back();
      pt[i] = j;
      pt[j] = (int)i;
    }
  }
}

// Turn "miRNA&target" plus its dot-bracket structure into three alignment
// rows: miRNA 3'->5', match symbols, target 5'->3'. Only inter-molecular pairs
// are drawn as pairs.
void format_aln(const std::string &seq, const char *structure, Align &out) {
  std::vector<int> pt;
  pair_table(structure, pt);
  const int n = (int)std::strlen(structure);
  out.ds.clear();
  out.cs.clear();
  out.aln.clear();

  int mir_hi = cut_point - 1;   // last unpaired miRNA position still to draw
  int tar_lo = cut_point + 1;   // first unpaired target position to draw
  int last_j = 0;               // miRNA partner of the last drawn pair
  bool prev_paired = false;
  std::string ms, ts;

  auto flush = [&]() {
    costring(ms, ts);
    out.ds += ms;
    out.cs += ts;
    out.aln.append(ms.size(), '*');
  };

  for (int i = cut_point + 1; i <= n; ++i) {
    int j = pt[i];
    if (j != 0 && j < cut_point) {
      if (prev_paired) {
        int d = pt[i - 1] - j;
        if (d != 1) {
          ms = seq.substr(j, d - 1);
          reverse(ms);
          ts.clear();
          flush();
        }
      } else {
        ms = seq.substr(j, mir_hi - j);
        reverse(ms);
        ts = seq.substr(tar_lo - 1, i - tar_lo);
        flush();
        prev_paired = true;
      }
      char a = seq[j - 1], b = seq[i - 1];
      out.ds += a;
      out.cs += b;
      out.aln += (a + b == 'G' + 'T' || a + b == 'G' + 'U') ? ':' : '|';
      last_j = j;
    } else if (prev_paired) {
      tar_lo = i;
      prev_paired = false;
      mir_hi = pt[i - 1] - 1;
    }
  }
  if (prev_paired) {
    if (last_j != 1) {
      ms = seq.substr(0, last_j - 1);
      reverse(ms);
      ts.clear();
      flush();
    }
  } else {
    ms = seq.substr(0, mir_hi);
    reverse(ms);
    ts = (size_t)(tar_lo - 1) <= seq.size() ? seq.substr(tar_lo - 1) : std::string();
    flush();
  }
}

float fold_pair(const std::string &a, const std::string &b, std::vector<char> &buf,
                std::string &seq) {
  seq = a;
  seq += '&';
  seq += b;
  // The '&' stays in the sequence handed to cofold (as an unpairable base at
  // the first position of the second strand). Reference outputs depend on it.
  cut_point = (int)a.size() + 1;
  buf.assign(seq.size() + 10, 0);
  float e = cofold(seq.c_str(), buf.data());
  buf[cut_point - 1] = '&';
  return e;
}

std::ostream &operator<<(std::ostream &os, const Align &a) {
  return os << "miRNA: " << a.ds << "\n       " << a.aln << "\n3'UTR: " << a.cs << "\n";
}

}  // namespace

Scanner::Scanner(const Params &p) : p_(p) {
  for (int i = 0; i < 5; ++i)
    for (int j = 0; j < 5; ++j) matrix_[i][j] = 0;
  // codes: A=1, T/U=2, C=3, G=4
  for (int i = 1; i <= 4; ++i)
    for (int j = 1; j <= 4; ++j) matrix_[i][j] = p_.mis_weight;
  matrix_[1][2] = matrix_[2][1] = p_.at_weight;
  matrix_[2][4] = matrix_[4][2] = p_.gt_weight;
  matrix_[3][4] = matrix_[4][3] = p_.gc_weight;

  temperature = p_.temperature;
  noLonelyPairs = p_.no_lonely_pairs;
  noGU = p_.no_gu;
  no_closingGU = p_.no_closing_gu;
  update_cofold_params();
}

static inline int base_code(char c) {
  switch (c) {
    case 'A': return 1;
    case 'T':
    case 'U': return 2;
    case 'C': return 3;
    case 'G': return 4;
  }
  return 0;
}

// Affine-gap alignment of a (miRNA beyond the seed, 5'->3') against b
// (target upstream of the seed, read 3'->5'). All of a is aligned; the far
// end of b is free. Returns the score, or INT_MIN on invalid input.
int Scanner::align(const char *a, const char *b, Align &out) const {
  const int n = (int)std::strlen(a), m = (int)std::strlen(b);
  if (n + 1 > kMaxAln - 1 || m + 1 > kMaxAln - 1) return -0x7fffffff - 1;
  static thread_local int M[kMaxAln][kMaxAln], E[kMaxAln][kMaxAln], F[kMaxAln][kMaxAln];
  static thread_local int Hl[kMaxAln][kMaxAln], Vl[kMaxAln][kMaxAln];
  static thread_local int tb[kMaxAln][kMaxAln][2];
  const int go = p_.gap_open, ge = p_.gap_ext;

  M[0][0] = 0;
  for (int i = 1; i <= n; ++i) {
    M[i][0] = go + (i - 1) * ge;
    E[i][0] = 2 * go - ge + (i - 1) * ge;
    Hl[i][0] = -1;
    tb[i][0][0] = -1;
    tb[i][0][1] = 0;
  }
  for (int j = 1; j <= m; ++j) {
    F[0][j] = 2 * go - ge + (j - 1) * ge;
    Vl[0][j] = -1;
    tb[0][j][0] = 0;
    M[0][j] = go + (j - 1) * ge;
    tb[0][j][1] = -1;
  }
  for (int i = 1; i <= n; ++i) {
    int ca = base_code(a[i - 1]);
    if (!ca) return -0x7fffffff - 1;
    for (int j = 1; j <= m; ++j) {
      int ext = ge + E[i][j - 1], open = go + M[i][j - 1];
      if (open > ext) {
        E[i][j] = open;
        Hl[i][j] = -1;
      } else {
        E[i][j] = ext;
        Hl[i][j] = Hl[i][j - 1] - 1;
      }
      open = go + M[i - 1][j];
      ext = ge + F[i - 1][j];
      if (open > ext) {
        F[i][j] = open;
        Vl[i][j] = -1;
      } else {
        F[i][j] = ext;
        Vl[i][j] = Vl[i - 1][j] - 1;
      }
      int cb = base_code(b[j - 1]);
      if (!cb) return -0x7fffffff - 1;
      int diag = M[i - 1][j - 1] + matrix_[ca][cb];
      int best = E[i][j];
      bool from_e = true;
      if (diag >= best) {
        best = diag;
        from_e = false;
      }
      if (F[i][j] <= best) {
        M[i][j] = best;
        if (from_e) {
          tb[i][j][0] = 0;
          tb[i][j][1] = Hl[i][j];
        } else {
          tb[i][j][0] = -1;
          tb[i][j][1] = -1;
        }
      } else {
        M[i][j] = F[i][j];
        tb[i][j][0] = Vl[i][j];
        tb[i][j][1] = 0;
      }
    }
  }
  int bj = 0, bs = M[n][0];
  for (int j = 0; j <= m; ++j)
    if (M[n][j] > bs) {
      bs = M[n][j];
      bj = j;
    }

  // traceback; rows are built from the far end, i.e. miRNA 3'->5'
  int i = n, j = bj;
  while (j > 0 && i > 0) {
    const int *t = tb[i][j];
    if (t[0] == 0) {
      int k = -t[1];
      if (k <= 0) break;
      for (int c = 0; c < k; ++c) {
        out.ds += '-';
        out.cs += b[j - 1 - c];
      }
      j -= k;
    } else if (t[1] == 0) {
      int k = -t[0];
      if (k <= 0) break;
      for (int c = 0; c < k; ++c) {
        out.ds += a[i - 1 - c];
        out.cs += '-';
      }
      i -= k;
    } else {
      --i;
      --j;
      out.ds += a[i];
      out.cs += b[j];
    }
  }
  for (; i > 0; --i) {
    out.ds += a[i - 1];
    out.cs += '-';
  }
  for (; j > 0; --j) {
    out.ds += '-';
    out.cs += b[j - 1];
  }
  get_aln(out.ds, out.cs, out.aln);
  return bs;
}

// Duplex energy of the miRNA with the aligned target stretch. If the cofold
// alignment leaves the miRNA 3' end hanging over the target start, the target
// is extended upstream by that many nucleotides and folded again.
float Scanner::get_energy(Result &r, const std::string &utrsub, const std::string &mir) {
  std::string tar;
  for (char ch : r.dynaln.cs)
    if (ch != '-') tar += ch;

  static thread_local std::vector<char> buf;
  std::string seq;
  Align out;
  out.ds = "-";
  out.cs = "-";
  float energy = 0;
  for (int guard = 0; guard < 1000; ++guard) {
    if (out.cs.empty() || out.cs[0] != '-') break;
    energy = fold_pair(mir, tar, buf, seq);
    format_aln(seq, buf.data(), out);
    if (out.ds.empty()) break;
    if (out.ds[0] != '-') {
      if (tar.size() == utrsub.size()) break;
      size_t p = leading_run(out.cs, '-');
      if (p == 0) break;
      size_t n = utrsub.size() - tar.size();
      size_t start, count;
      if ((long)(n - p) < 0 || p == std::string::npos) {
        // Not enough upstream sequence: take all of it. (p == npos cannot
        // happen for a non-empty target; treated the same way to stay safe.)
        start = 0;
        count = n;
      } else {
        start = n - p;
        count = p;
      }
      tar = utrsub.substr(start, count) + tar;
    } else {
      size_t p = leading_run(out.ds, '-');
      out.ds.erase(0, std::min(p, out.ds.size()));
      out.cs.erase(0, std::min(p, out.cs.size()));
      out.aln.erase(0, std::min(p, out.aln.size()));
      break;
    }
  }
  r.eryaln = out;
  r.structure = buf.data();
  return energy;
}

void Scanner::score_loop(Result &r) const {
  const Align &e = r.eryaln;
  std::string s;
  size_t p = e.ds.find_last_not_of('-');
  if (p != e.ds.size() - 1 || (p = e.cs.find_last_not_of('-')) != e.cs.size() - 1)
    s = e.aln.substr(0, p + 1);
  else
    s = e.aln;
  reverse(s);  // now index 0 is the miRNA 5' end

  // Unpaired runs ("loops") as [start (1-based), length], beginning with the
  // run that contains the first '*' at or after index 6.
  struct Loop {
    int start, len;
  };
  std::vector<Loop> L;
  size_t pos = s.find_last_not_of('*', s.find('*', 6)) + 1;
  for (;;) {
    size_t q = s.find_first_not_of('*', pos);
    if ((int)q <= 0) {
      L.push_back({(int)pos + 1, (int)s.size() - (int)pos});
      break;
    }
    L.push_back({(int)pos + 1, (int)(q - pos)});
    pos = s.find('*', q + 1);
    if (pos == std::string::npos) break;
  }
  const int nl = (int)L.size();
  // Sentinel entry so the merge below can look one past the last loop.
  L.push_back({0x40000000, 0});
  // Merge two loops separated by a single pair.
  for (int k = 0; k < nl; ++k) {
    if (L[k + 1].start - L[k].start == 2 && L[k + 1].len + L[k].len != 0) {
      L[k + 1].start = L[k].start;
      L[k + 1].len = L[k + 1].len + L[k].len + 1;
      L[k].start = -1;
    }
  }
  int sel = 0;
  for (int k = 0; k < nl; ++k)
    if (L[k].start > 0 && L[k].len > 1) {
      sel = k;
      break;
    }
  const int start = L[sel].start, len = L[sel].len;
  float score;
  if (start <= 6)
    score = 0;
  else if (len > 0)
    score = (float)((double)kLenScore[len < 6 ? len - 1 : 4] * (double)kTypeScore[start < 16 ? start - 7 : 8]);
  else
    score = -1.0f;
  r.loopscore = score;
  r.loop[0] = start;
  r.loop[1] = start + len - 1;
  float extra = 0;
  if ((unsigned)r.loop[1] < (unsigned)(s.size() - 1)) {
    // perfectly paired 3' part after a small loop
    if (s.find('*', (size_t)(unsigned)r.loop[1]) == std::string::npos && len <= 3) extra = -7.5f;
  }
  r.loopscore = (float)((double)r.loopscore + (double)extra);
  r.tr_score = extra;
  r.ol_score = 0;
}

bool Scanner::get_result(const std::string &utrsub, const std::string &mir, Result &r) {
  const int se = p_.seed_end;
  std::string t = utrsub;
  reverse(t);

  r.dynaln = Align();
  r.score = align(mir.c_str() + se, t.c_str() + se, r.dynaln);
  if (r.score == -0x7fffffff - 1) return false;

  std::string s6 = mir.substr(0, se);
  std::string s5 = t.substr(0, se);
  reverse(s5);
  {
    static thread_local std::vector<char> buf;
    std::string seq;
    r.seedery = fold_pair(s6, s5, buf, seq);
  }
  reverse(s6);
  std::string m;
  r.seedscore = get_aln(s6, s5, m);
  r.dynaln.ds += s6;
  r.dynaln.aln += m;
  r.dynaln.cs += s5;

  int used = 0;
  for (char ch : r.dynaln.cs) used += ch != '-';
  r.pos[0] = r.pos[1] + 1 - used;

  r.energy = get_energy(r, utrsub, mir);
  score_loop(r);
  return true;
}

int Scanner::scan_impl(const Sequence &utr, const Sequence &mir,
                       void (*emit)(void *, const Result &), void *ctx) {
  const std::string &u = utr.seq, &m = mir.seq;
  if ((int)m.size() < p_.seed_end || p_.seed_end < 8) return -1;

  // Seed scan: hit[i] > 0 when a seed match puts miRNA nt 1 opposite target
  // position i (1-based), i.e. the candidate site ends at i.
  std::vector<int> hit(u.size() + 10, 0);
  std::string pat[7];
  for (int k = 0; k < 3; ++k) {
    if (!get_pattern(m.substr(k, 6), pat)) return -1;
    for (int v = 0; v < 7; ++v) {
      if (pat[v].empty()) continue;
      for (size_t at = u.find(pat[v]); at != std::string::npos; at = u.find(pat[v], at + 1))
        ++hit[at + 6 + k];
    }
  }

  int count = 0;
  Result r;
  r.utrid = utr.id;
  r.mirid = mir.id;
  for (size_t i = 0; i < u.size(); ++i) {
    if (hit[i] <= 0) continue;
    std::string sub;
    if ((int)i >= kWindow)
      sub = u.substr(i - kWindow, kWindow);
    else if ((int)i >= kMinWindow)
      sub = u.substr(0, i);
    else
      continue;
    r.pos[1] = (int)i;
    if (!get_result(sub, m, r)) return -1;
    ++count;
    if (r.energy <= p_.filter_energy && r.score >= p_.filter_score && p_.filter_loop <= r.loopscore)
      emit(ctx, r);
  }
  return count;
}

void Scanner::write(std::ostream &os, const Result &r) const {
  if (p_.format != 't') {
    os << r.mirid << "\t" << r.score << "\t" << (double)r.energy << "\t" << r.pos[0] << "\t"
       << r.pos[1] << "\t" << r.dynaln.ds << "," << r.dynaln.aln << "," << r.dynaln.cs << "\t"
       << r.eryaln.ds << "," << r.eryaln.aln << "," << r.eryaln.cs << "\t" << r.seedscore << "\t"
       << (double)r.seedery << "\t" << (double)r.loopscore << "[" << r.loop[0] << "," << r.loop[1]
       << "|" << (double)r.tr_score << "," << (double)r.ol_score << "]"
       << "\t" << r.structure << "\n";
  } else {
    os << "Position: " << r.pos[0] << ", " << r.pos[1] << "   Score: " << r.score
       << "    Energy: " << (double)r.energy << "   seed: " << r.seedscore << ", "
       << (double)r.seedery << "\n\n"
       << "dynamic programming structure:\n"
       << r.dynaln << "\nduplex structure:\n"
       << r.eryaln << "\nSecondary structure:" << r.structure << "\n\n";
  }
}

void Scanner::write_header(std::ostream &os) const {
  os << "# version     " << kVersion << "\n"
     << "# ATWeight    " << p_.at_weight << "\n"
     << "# GCWeight    " << p_.gc_weight << "\n"
     << "# GTWeight    " << p_.gt_weight << "\n"
     << "# MisWeight   " << p_.mis_weight << "\n"
     << "# GapOpen     " << p_.gap_open << "\n"
     << "# GapExt      " << p_.gap_ext << "\n"
     << "# SeedEnd     " << p_.seed_end << "\n";
  if (p_.filter_energy != 1000) os << "# Energy      " << (double)p_.filter_energy << "\n";
  if (p_.filter_score != -1000) os << "# Score       " << p_.filter_score << "\n";
  if (p_.filter_loop != -1000) os << "# Loop Score  " << (double)p_.filter_loop << "\n";
  os << "# Temperature " << p_.temperature << "\n";
  os << "# Fields: mirna score energy start end align energy_align seed_score seed_energy loop "
        "structure\n";
}

void read_fasta(std::istream &in, std::vector<Sequence> &out) {
  std::string line;
  auto chomp = [](std::string &s) {
    size_t p = s.find('\r');
    if (p != std::string::npos) s.erase(p);
  };
  bool have = false;
  while (std::getline(in, line)) {
    chomp(line);
    if (!line.empty() && line[0] == '>') {
      have = true;
      break;
    }
  }
  while (have) {
    Sequence s;
    s.id = line.substr(1);
    have = false;
    while (std::getline(in, line)) {
      chomp(line);
      if (!line.empty() && line[0] == '>') {
        have = true;
        break;
      }
      s.seq += line;
    }
    for (char &ch : s.seq) ch = (char)std::toupper((unsigned char)ch);
    out.push_back(s);
  }
}

bool run(const std::vector<Sequence> &mirs, const std::vector<Sequence> &utrs, const Params &p,
         std::ostream &os) {
  Scanner sc(p);
  sc.write_header(os);
  const char *bar = "=============================================================";
  for (size_t ui = 0; ui < utrs.size(); ++ui) {
    const Sequence &utr = utrs[ui];
    os << utr.id << "\n";
    if (!utr.seq.empty()) {
      for (const Sequence &mir : mirs) {
        if (mir.seq.empty()) continue;
        if (p.format == 't') os << bar << "\n" << mir.id << " " << utr.id << "\n";
        int n = sc.scan(utr, mir, [&](const Result &r) { sc.write(os, r); });
        if (n < 0) return false;  // non-nucleotide character
        if (p.format == 't')
          os << "<" << utr.id << ">" << mir.id << ": " << n << " targets sites\n" << bar << "\n\n";
      }
    }
    if (ui + 1 < utrs.size()) os << "\n";  // blank line between UTRs, none after the last
  }
  return true;
}

}  // namespace mirtherm
