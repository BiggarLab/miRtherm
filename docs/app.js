'use strict';
// miRtherm dashboard. All computation happens in Web Workers running the
// WebAssembly build of the core; nothing leaves the browser.

const $ = (id) => document.getElementById(id);
const PAGE = 40;
const SCALE_MIN = -5, SCALE_MAX = 50;
const COLD = [74, 79, 196], MID = [184, 65, 138], WARM = [224, 112, 31];

const state = {
  utrs: null,
  utrLabel: '',
  run: null,     // {temps, mirs, text, rows, cutoff}
  view: { rows: [], page: 0, sort: 'e0', dir: 1, sel: null, bin: null },
  workers: [],
  cancelled: false,
};

// ---------- input parsing ----------

// The id is the whole header line, text after a carriage return is dropped,
// sequences are upper-cased. The command line reader follows the same rules.
function parseFasta(text) {
  const out = [];
  let cur = null;
  for (let line of text.split('\n')) {
    const cr = line.indexOf('\r');
    if (cr >= 0) line = line.slice(0, cr);
    if (line[0] === '>') {
      cur = { id: line.slice(1), seq: '' };
      out.push(cur);
    } else if (cur) {
      cur.seq += line;
    }
  }
  for (const s of out) s.seq = s.seq.toUpperCase();
  return out;
}

function parseMirs(text) {
  let mirs;
  if (text.includes('>')) mirs = parseFasta(text);
  else mirs = text.split(/\s+/).filter(Boolean).map((s, i) => ({ id: 'miRNA-' + (i + 1), seq: s.toUpperCase() }));
  mirs = mirs.filter((m) => m.seq.length);
  if (!mirs.length) throw new Error('Enter at least one miRNA sequence.');
  for (const m of mirs) {
    m.seq = m.seq.replace(/\s+/g, '');
    if (!/^[ACGTU]+$/.test(m.seq)) throw new Error(`“${m.id}” has characters other than A, C, G, T or U.`);
    if (m.seq.length < 9 || m.seq.length > 56) throw new Error(`“${m.id}” is ${m.seq.length} nt long. miRNAs must be 9 to 56 nt.`);
  }
  return mirs;
}

function setUtrs(text, label) {
  const utrs = parseFasta(text);
  if (!utrs.length) throw new Error('No FASTA records found in that file.');
  for (const u of utrs) {
    u.seq = u.seq.replace(/\s+/g, '');
    u.name = u.id.split(/\s+/)[0] || u.id;
    const g = /(?:parent|gene|gene_id|gene_symbol)[=:]([^;\s,]+)/i.exec(u.id);
    u.gene = g ? g[1] : '';
  }
  state.utrs = utrs;
  state.utrLabel = label;
  const nt = utrs.reduce((a, u) => a + u.seq.length, 0);
  $('utrinfo').textContent = `${label}: ${utrs.length.toLocaleString()} sequences, ${(nt / 1e6).toFixed(2)} Mb.`;
}

function readTemps() {
  const parts = $('temps').value.split(/[,;\s]+/).filter(Boolean).map(Number);
  if (!parts.length || parts.some((t) => !Number.isFinite(t) || t < -20 || t > 100)) return null;
  return [...new Set(parts)].sort((a, b) => a - b);
}

function readSettings() {
  const num = (id, fallback) => {
    const v = $(id).value.trim();
    if (v === '') return fallback;
    const n = Number(v);
    if (!Number.isFinite(n)) throw new Error(`“${v}” is not a number.`);
    return n;
  };
  const temps = readTemps();
  if (!temps) throw new Error('Temperatures must be numbers between -20 and 100 °C.');
  if (temps.length > 6) throw new Error('Compare at most six temperatures at a time.');
  const base = {
    at: Math.trunc(num('at', 5)), gc: Math.trunc(num('gc', 5)), gt: Math.trunc(num('gt', 2)),
    mis: Math.trunc(num('mis', -3)), go: Math.trunc(num('go', -8)), ge: Math.trunc(num('ge', -2)),
    energy: num('energy', 1000), score: Math.trunc(num('score', -1000)), loop: num('loop', -1000),
    noLP: $('noLP').checked, noGU: $('noGU').checked, noCloseGU: $('noCloseGU').checked,
  };
  return { temps, base, threads: Math.max(1, Math.min(64, Math.trunc(num('threads', 1)))) };
}

// ---------- temperature scale ----------

const xOf = (t) => (Math.min(SCALE_MAX, Math.max(SCALE_MIN, t)) - SCALE_MIN) / (SCALE_MAX - SCALE_MIN) * 100;

function drawScale() {
  if (!$('ticks').childElementCount) {
    let h = '';
    for (let t = SCALE_MIN; t <= SCALE_MAX; t += 5) h += `<span style="left:${xOf(t)}%">${t}</span>`;
    $('ticks').innerHTML = h;
  }
  const temps = readTemps() || [];
  $('pins').innerHTML = temps.map((t) =>
    `<button type="button" class="pin" style="left:${xOf(t)}%" data-t="${t}" aria-label="Remove ${t} °C">${t}°</button>`).join('');
}

function writeTemps(list) {
  $('temps').value = [...new Set(list)].sort((a, b) => a - b).join(', ');
  drawScale();
}

// ---------- worker pool ----------

function makePool(n) {
  const workers = [];
  for (let i = 0; i < n; i++) workers.push(new Worker('worker.js'));
  state.workers = workers;
  let seq = 0;
  const pending = new Map(), idle = workers.slice(), queue = [];
  const pump = () => {
    while (idle.length && queue.length) {
      const w = idle.pop(), job = queue.shift();
      pending.set(job.msg.id, job);
      w.postMessage(job.msg);
    }
  };
  for (const w of workers) {
    w.onmessage = (e) => {
      const job = pending.get(e.data.id);
      pending.delete(e.data.id);
      idle.push(w);
      job.resolve(e.data);
      pump();
    };
    w.onerror = (e) => {
      const err = new Error(e.message || 'A worker failed to start.');
      for (const j of pending.values()) j.reject(err);
      for (const j of queue) j.reject(err);
    };
  }
  return {
    post(msg) {
      return new Promise((resolve, reject) => {
        msg.id = ++seq;
        queue.push({ msg, resolve, reject });
        pump();
      });
    },
    close() {
      for (const w of workers) w.terminate();
      state.workers = [];
    },
  };
}

function chunkUtrs(utrs) {
  const chunks = [];
  let start = 0, len = 0;
  for (let i = 0; i < utrs.length; i++) {
    len += utrs[i].seq.length;
    if (len > 60000 || i - start + 1 >= 150) {
      chunks.push([start, i + 1]);
      start = i + 1;
      len = 0;
    }
  }
  if (start < utrs.length) chunks.push([start, utrs.length]);
  return chunks;
}

// ---------- run ----------

function progress(frac, msg) {
  $('meter').hidden = false;
  $('meterfill').style.width = (frac * 100).toFixed(1) + '%';
  $('status').textContent = msg;
}

async function runAll() {
  const mirs = parseMirs($('mirs').value);
  if (!state.utrs) throw new Error('Choose a FASTA file of target sequences first.');
  const { temps, base, threads } = readSettings();
  const utrs = state.utrs;
  const pool = makePool(threads);
  state.cancelled = false;
  const chunks = chunkUtrs(utrs);
  const total = chunks.length * temps.length;
  let done = 0;
  const t0 = performance.now();
  const text = [], blocksByTemp = [], skipped = new Set();
  progress(0, 'Starting');
  try {
    for (const T of temps) {
      const cfg = Object.assign({ temperature: T }, base);
      const header = (await pool.post({ type: 'header', cfg })).header;
      const blocks = new Array(utrs.length);
      await Promise.all(chunks.map(async ([a, b]) => {
        const res = await pool.post({ type: 'scan', cfg, mirs, utrs: utrs.slice(a, b).map((u) => u.seq) });
        for (let i = a; i < b; i++) {
          if (res.out[i - a] === null) skipped.add(i);
          blocks[i] = res.out[i - a] || '';
        }
        done++;
        const left = (performance.now() - t0) / 1000 / done * (total - done);
        progress(done / total, `Folding at ${T} °C, about ${Math.ceil(left)} s left`);
      }));
      blocksByTemp.push(blocks);
      // Same layout as the command line output file.
      text.push(header + utrs.map((u, i) => u.id + '\n' + blocks[i]).join('\n'));
    }
  } finally {
    pool.close();
  }
  if (state.cancelled) return;
  const secs = (performance.now() - t0) / 1000;
  state.run = { temps, mirs, text, rows: mergeRows(blocksByTemp), cutoff: base.energy };
  progress(1, `Done in ${secs.toFixed(1)} s on ${threads} thread${threads > 1 ? 's' : ''}.` +
    (skipped.size ? ` Skipped ${skipped.size} sequence${skipped.size > 1 ? 's' : ''} with characters other than A, C, G, T or U.` : ''));
  renderResults();
}

// One row per site (transcript, miRNA, end position) with a slot per temperature.
function mergeRows(blocksByTemp) {
  const map = new Map();
  blocksByTemp.forEach((blocks, ti) => {
    blocks.forEach((blk, ui) => {
      if (!blk) return;
      for (const line of blk.split('\n')) {
        if (!line) continue;
        const f = line.split('\t');
        const key = ui + '|' + f[0] + '|' + f[4];
        let row = map.get(key);
        if (!row) {
          row = { u: ui, mir: f[0], start: +f[3], end: +f[4], score: +f[1], dyn: f[5], seedScore: +f[7], at: [] };
          map.set(key, row);
        }
        row.at[ti] = { energy: +f[2], aln: f[6], seedEnergy: +f[8], loop: f[9], structure: f[10] };
      }
    });
  });
  return [...map.values()];
}

// ---------- rendering ----------

function tempColor(i, n, alpha) {
  let c;
  if (n < 2) c = COLD;
  else {
    const t = i / (n - 1);
    const [a, b, f] = t <= 0.5 ? [COLD, MID, t * 2] : [MID, WARM, (t - 0.5) * 2];
    c = n === 2 ? (i ? WARM : COLD) : a.map((v, k) => Math.round(v + (b[k] - v) * f));
  }
  return alpha === undefined ? `rgb(${c.join(',')})` : `rgba(${c.join(',')},${alpha})`;
}
const deg = (t) => `${t} °C`;
const esc = (s) => String(s).replace(/[&<>"]/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
const minus = (x) => String(x).replace('-', '−');

function renderResults() {
  const { temps, rows } = state.run;
  const utrs = state.utrs;
  $('results').hidden = false;
  const n = temps.length;

  // tally
  const tx = temps.map(() => new Set());
  const sites = temps.map(() => 0);
  for (const r of rows) for (let i = 0; i < n; i++) if (r.at[i]) { sites[i]++; tx[i].add(r.u); }
  let h = temps.map((T, i) =>
    `<div class="tally-row"><div class="tally-temp"><span class="swatch" style="background:${tempColor(i, n)}"></span>${deg(T)}</div>` +
    `<div class="tally-big">${tx[i].size.toLocaleString()}</div>` +
    `<div class="tally-sub">transcripts targeted, ${sites[i].toLocaleString()} sites</div>` +
    `<button type="button" data-dl="${i}">Download ${deg(T)} output</button></div>`).join('');
  if (n >= 2) {
    const a = 0, b = n - 1;
    let lost = 0, gained = 0, kept = 0;
    for (const u of tx[a]) tx[b].has(u) ? kept++ : lost++;
    for (const u of tx[b]) if (!tx[a].has(u)) gained++;
    h += `<div class="shift"><p>Warming from ${deg(temps[a])} to ${deg(temps[b])}:</p>` +
      `<p><b>${lost.toLocaleString()}</b> transcripts drop out</p><p><b>${kept.toLocaleString()}</b> stay</p><p><b>${gained.toLocaleString()}</b> appear</p></div>`;
  }
  $('tally').innerHTML = h;

  drawLanes();

  let opts = '<option value="all">All sites</option>';
  if (n > 1) {
    opts += '<option value="every">Sites at every temperature</option>';
    temps.forEach((T, i) => { opts += `<option value="only${i}">Sites only at ${deg(T)}</option>`; });
  }
  $('show').innerHTML = opts;
  $('show').hidden = n < 2;

  const hasGene = utrs.some((u) => u.gene);
  state.cols = [
    { k: 'name', label: 'Transcript', text: true, get: (r) => utrs[r.u].name },
    ...(hasGene ? [{ k: 'gene', label: 'Gene', text: true, get: (r) => utrs[r.u].gene }] : []),
    ...(state.run.mirs.length > 1 ? [{ k: 'mir', label: 'miRNA', text: true, get: (r) => r.mir }] : []),
    { k: 'start', label: 'Start', get: (r) => r.start },
    { k: 'end', label: 'End', get: (r) => r.end },
    ...temps.map((T, i) => ({ k: 'e' + i, label: `ΔG at ${deg(T)}`, energy: i, get: (r) => (r.at[i] ? r.at[i].energy : null) })),
    { k: 'score', label: 'Alignment score', get: (r) => r.score },
  ];
  $('sites').querySelector('thead').innerHTML = '<tr>' + state.cols.map((c) =>
    `<th scope="col" data-k="${c.k}"${c.text ? ' class="t"' : ''}>${esc(c.label)}</th>`).join('') + '</tr>';
  state.view = { rows: [], page: 0, sort: 'e0', dir: 1, sel: null, bin: null };
  $('search').value = '';
  $('detail').hidden = true;
  applyView();
  $('results').scrollIntoView({ block: 'start', behavior: 'smooth' });
}

// Small multiples: one horizontal histogram of site energies per temperature,
// sharing both axes so lanes can be compared directly.
function drawLanes() {
  const { temps, rows } = state.run;
  const n = temps.length;
  let lo = Infinity, hi = -Infinity;
  for (const r of rows) for (const s of r.at) if (s) { lo = Math.min(lo, s.energy); hi = Math.max(hi, s.energy); }
  if (!Number.isFinite(lo)) { $('lanes').innerHTML = '<p class="note">No sites passed the filters.</p>'; state.bins = null; return; }
  const top = Math.ceil(hi), bottom = Math.floor(lo);          // top = weakest
  const nb = Math.max(1, top - bottom);
  const counts = temps.map(() => new Array(nb).fill(0));
  for (const r of rows) r.at.forEach((s, i) => { if (s) counts[i][Math.min(nb - 1, Math.floor(top - s.energy))]++; });
  const max = Math.max(1, ...counts.flat());
  const rowH = Math.max(6, Math.min(16, Math.floor(340 / nb))), gap = rowH > 8 ? 2 : 1;
  const laneW = 190, axisW = 46, headH = 26, padB = 22;
  const W = axisW + n * laneW, H = headH + nb * rowH + padB;
  const step = nb > 30 ? 10 : 5;
  let s = `<svg width="${W}" height="${H}" viewBox="0 0 ${W} ${H}" role="img" aria-label="Histograms of duplex energy at each temperature">`;
  for (let e = top; e >= bottom; e--) {
    if (e % step) continue;
    const y = headH + (top - e) * rowH;
    s += `<line class="grid" x1="${axisW}" x2="${W}" y1="${y}" y2="${y}"/><text x="${axisW - 8}" y="${y + 4}" text-anchor="end">${minus(e)}</text>`;
  }
  s += `<text x="${axisW - 8}" y="14" text-anchor="end">ΔG</text>`;
  temps.forEach((T, i) => {
    const x0 = axisW + i * laneW;
    const col = tempColor(i, n);
    s += `<text class="lane-title" x="${x0 + 6}" y="14">${deg(T)}</text><line class="base" x1="${x0}" x2="${x0}" y1="${headH}" y2="${headH + nb * rowH}"/>`;
    for (let b = 0; b < nb; b++) {
      const c = counts[i][b];
      if (!c) continue;
      const w = Math.max(2, c / max * (laneW - 34)), y = headH + b * rowH;
      const r = Math.min(4, (rowH - gap) / 2);
      s += `<rect class="hit" x="${x0}" y="${y}" width="${laneW - 6}" height="${rowH}" data-t="${i}" data-b="${b}" data-c="${c}"/>` +
        `<path class="bar" data-key="${i}:${b}" fill="${col}" d="M${x0 + 1},${y + gap / 2}h${w - r}a${r},${r} 0 0 1 ${r},${r}v${rowH - gap - 2 * r}a${r},${r} 0 0 1 ${-r},${r}h${-(w - r)}z"/>`;
    }
    s += `<text x="${x0 + 6}" y="${H - 5}">most in one bin: ${Math.max(...counts[i]).toLocaleString()}</text>`;
  });
  $('lanes').innerHTML = s + '</svg>';
  state.bins = { top, nb };
}

function applyView() {
  const { rows, temps } = state.run;
  const utrs = state.utrs, v = state.view;
  const q = $('search').value.trim().toLowerCase();
  const show = $('show').value || 'all';
  const out = rows.filter((r) => {
    if (show === 'every') { for (let i = 0; i < temps.length; i++) if (!r.at[i]) return false; }
    else if (show.startsWith('only')) {
      const k = +show.slice(4);
      for (let i = 0; i < temps.length; i++) if (!!r.at[i] !== (i === k)) return false;
    }
    if (v.bin) {
      const s = r.at[v.bin.t];
      if (!s || s.energy > v.bin.hi || s.energy <= v.bin.lo) return false;
    }
    if (q && !(utrs[r.u].id.toLowerCase().includes(q) || r.mir.toLowerCase().includes(q))) return false;
    return true;
  });
  const col = state.cols.find((c) => c.k === v.sort);
  if (col) {
    out.sort((a, b) => {
      const x = col.get(a), y = col.get(b);
      if (x === null || y === null) return (x === null) - (y === null);
      return (x < y ? -1 : x > y ? 1 : 0) * v.dir || a.u - b.u || a.end - b.end;
    });
  }
  v.rows = out;
  v.page = Math.min(v.page, Math.max(0, Math.ceil(out.length / PAGE) - 1));
  for (const th of $('sites').querySelectorAll('th')) {
    if (th.dataset.k === v.sort) th.setAttribute('aria-sort', v.dir > 0 ? 'ascending' : 'descending');
    else th.removeAttribute('aria-sort');
  }
  const chip = $('binchip');
  chip.hidden = !v.bin;
  if (v.bin) chip.textContent = `${deg(temps[v.bin.t])}, ${minus(v.bin.lo)} to ${minus(v.bin.hi)} kcal/mol`;
  for (const p of $('lanes').querySelectorAll('.bar')) p.classList.toggle('on', !!v.bin && p.dataset.key === v.bin.t + ':' + v.bin.b);
  renderPage();
}

function renderPage() {
  const v = state.view, n = state.run.temps.length;
  const span = state.bins ? Math.max(1, state.bins.nb) : 1;
  const slice = v.rows.slice(v.page * PAGE, (v.page + 1) * PAGE);
  $('sites').querySelector('tbody').innerHTML = slice.map((r, i) =>
    `<tr data-i="${v.page * PAGE + i}"${r === v.sel ? ' class="sel"' : ''}>` + state.cols.map((c) => {
      const val = c.get(r);
      if (val === null) return '<td class="none">not found</td>';
      if (c.energy !== undefined) {
        // small bar under the number: how far below the weakest site
        const w = Math.max(4, Math.min(100, (state.bins.top - val) / span * 100));
        return `<td class="e">${minus(val)}<i style="width:calc(${w}% - 1.8rem);background:${tempColor(c.energy, n)}"></i></td>`;
      }
      return `<td${c.text ? ' class="t"' : ''}>${esc(c.text ? val : minus(val))}</td>`;
    }).join('') + '</tr>').join('');
  const pages = Math.max(1, Math.ceil(v.rows.length / PAGE));
  $('count').textContent = `${v.rows.length.toLocaleString()} site${v.rows.length === 1 ? '' : 's'}`;
  $('pageinfo').textContent = `Page ${v.page + 1} of ${pages.toLocaleString()}`;
  $('prev').disabled = v.page === 0;
  $('next').disabled = v.page >= pages - 1;
}

// A duplex drawn column by column: miRNA on top (3′ to 5′), target below.
function duplexHtml(aln, tint) {
  const [top, mid, bot] = aln.split(',');
  let cells = '';
  for (let i = 0; i < top.length; i++) {
    const m = mid[i], paired = m === '|' || m === ':';
    const nt = (ch) => ch === '-' ? '<span class="gap">·</span>' : `<span class="nt${paired ? ' p' : ''}">${ch}</span>`;
    cells += nt(top[i]) + (paired ? `<span class="bond${m === ':' ? ' w' : ''}"></span>` : '<span></span>') + nt(bot[i]);
  }
  return `<div class="strand" style="--tint:${tint}">` +
    `<span class="rowlabel">miRNA <span class="end">3′</span></span><div class="cells">${cells}</div><span class="end">5′</span>` +
    `<span></span><span></span>` +
    `<span class="rowlabel">target <span class="end">5′</span></span><span class="end">3′</span></div>`;
}

function renderDetail(r) {
  const { temps } = state.run;
  const n = temps.length, u = state.utrs[r.u];
  let h = `<h2>${esc(r.mir)} on ${esc(u.name)}, nucleotides ${r.start} to ${r.end}</h2><p class="where">${esc(u.id)}</p><div class="duplexes">`;
  temps.forEach((T, i) => {
    const s = r.at[i];
    h += `<div class="dx"><h3><span class="swatch" style="background:${tempColor(i, n)}"></span>${deg(T)}</h3>`;
    if (!s) { h += '<p class="absent">Not kept at this temperature. The duplex is weaker than the cutoff.</p></div>'; return; }
    h += `<p class="facts">ΔG ${minus(s.energy)} kcal/mol, seed ${minus(s.seedEnergy)}, loop ${esc(s.loop)}</p>` +
      duplexHtml(s.aln, tempColor(i, n, 0.22)) + `<p class="struct">${esc(s.structure)}</p></div>`;
  });
  h += `<div class="dx"><h3>Sequence alignment alone</h3><p class="facts">Score ${minus(r.score)}, seed score ${r.seedScore}</p>` +
    duplexHtml(r.dyn, 'rgba(36,26,53,0.12)') + '</div></div>';
  const d = $('detail');
  d.innerHTML = h;
  d.hidden = false;
}

// ---------- downloads ----------

function download(name, text, type) {
  const url = URL.createObjectURL(new Blob([text], { type: type || 'text/plain' }));
  const a = document.createElement('a');
  a.href = url;
  a.download = name;
  document.body.appendChild(a);
  a.click();
  a.remove();
  setTimeout(() => URL.revokeObjectURL(url), 5000);
}

function csvOf(rows) {
  const { temps } = state.run;
  const q = (s) => '"' + String(s).replace(/"/g, '""') + '"';
  const lines = [['transcript', 'gene', 'mirna', 'start', 'end', 'score', 'seed_score',
    ...temps.flatMap((T) => [`energy_${T}C`, `seed_energy_${T}C`, `loop_${T}C`, `alignment_${T}C`, `structure_${T}C`]),
    'sequence_alignment', 'target_id'].join(',')];
  for (const r of rows) {
    const u = state.utrs[r.u];
    const cells = [q(u.name), q(u.gene), q(r.mir), r.start, r.end, r.score, r.seedScore];
    temps.forEach((T, i) => {
      const s = r.at[i];
      if (s) cells.push(s.energy, s.seedEnergy, q(s.loop), q(s.aln), q(s.structure));
      else cells.push('', '', '', '', '');
    });
    cells.push(q(r.dyn), q(u.id));
    lines.push(cells.join(','));
  }
  return lines.join('\n') + '\n';
}

// ---------- events ----------

function fail(err) {
  $('error').textContent = err.message || String(err);
  $('error').hidden = false;
}

$('threads').value = Math.max(1, (navigator.hardwareConcurrency || 2) - 1);
drawScale();

$('temps').addEventListener('input', drawScale);
$('scale').addEventListener('click', (e) => {
  const pin = e.target.closest('.pin');
  const temps = readTemps() || [];
  if (pin) { writeTemps(temps.filter((t) => t !== +pin.dataset.t)); return; }
  const box = $('track').getBoundingClientRect();
  const t = Math.round(SCALE_MIN + (e.clientX - box.left) / box.width * (SCALE_MAX - SCALE_MIN));
  if (t >= SCALE_MIN && t <= SCALE_MAX) writeTemps(temps.concat(t));
});

$('more').addEventListener('click', () => {
  const adv = $('advanced');
  adv.hidden = !adv.hidden;
  $('more').setAttribute('aria-expanded', String(!adv.hidden));
  $('more').textContent = adv.hidden ? 'More settings' : 'Fewer settings';
});

$('utrfile').addEventListener('change', async (e) => {
  const f = e.target.files[0];
  if (!f) return;
  $('error').hidden = true;
  try { setUtrs(await f.text(), f.name); } catch (err) { fail(err); }
});

$('example').addEventListener('click', async () => {
  $('error').hidden = true;
  try {
    const res = await fetch('example_utrs.fasta');
    if (!res.ok) throw new Error('Could not load the example sequences.');
    setUtrs(await res.text(), 'Fly example');
  } catch (err) { fail(err); }
});

$('setup').addEventListener('submit', async (e) => {
  e.preventDefault();
  $('error').hidden = true;
  $('run').disabled = true;
  $('cancel').hidden = false;
  try {
    await runAll();
  } catch (err) {
    if (!state.cancelled) { fail(err); $('meter').hidden = true; $('status').textContent = ''; }
  } finally {
    $('run').disabled = false;
    $('cancel').hidden = true;
  }
});

$('cancel').addEventListener('click', () => {
  state.cancelled = true;
  for (const w of state.workers) w.terminate();
  state.workers = [];
  $('run').disabled = false;
  $('cancel').hidden = true;
  progress(0, 'Stopped.');
});

$('tally').addEventListener('click', (e) => {
  const ti = e.target.dataset.dl;
  if (ti === undefined) return;
  const T = state.run.temps[+ti];
  const mir = state.run.mirs.length === 1 ? state.run.mirs[0].id.split(/\s+/)[0].replace(/[^\w.-]+/g, '_') : 'targets';
  download(`${mir}-${T}.txt`, state.run.text[+ti]);
});

const tip = $('tip');
$('lanes').addEventListener('mousemove', (e) => {
  const hit = e.target.closest('.hit');
  if (!hit || !state.bins) { tip.hidden = true; return; }
  const hi = state.bins.top - +hit.dataset.b, T = state.run.temps[+hit.dataset.t];
  tip.innerHTML = `<b>${(+hit.dataset.c).toLocaleString()}</b> sites at ${deg(T)}<br>${minus(hi - 1)} to ${minus(hi)} kcal/mol`;
  tip.hidden = false;
  tip.style.left = Math.min(window.innerWidth - 200, e.clientX + 14) + 'px';
  tip.style.top = e.clientY + 14 + 'px';
});
$('lanes').addEventListener('mouseleave', () => { tip.hidden = true; });
$('lanes').addEventListener('click', (e) => {
  const hit = e.target.closest('.hit');
  if (!hit) return;
  const b = +hit.dataset.b, t = +hit.dataset.t, hi = state.bins.top - b;
  const same = state.view.bin && state.view.bin.t === t && state.view.bin.b === b;
  state.view.bin = same ? null : { t, b, hi, lo: hi - 1 };
  state.view.sort = 'e' + t;
  state.view.dir = 1;
  state.view.page = 0;
  applyView();
});
$('binchip').addEventListener('click', () => { state.view.bin = null; state.view.page = 0; applyView(); });

$('csv').addEventListener('click', () => download('mirtherm-sites.csv', csvOf(state.view.rows), 'text/csv'));
$('search').addEventListener('input', () => { state.view.page = 0; applyView(); });
$('show').addEventListener('change', () => { state.view.page = 0; applyView(); });
$('prev').addEventListener('click', () => { state.view.page--; renderPage(); });
$('next').addEventListener('click', () => { state.view.page++; renderPage(); });

$('sites').querySelector('thead').addEventListener('click', (e) => {
  const th = e.target.closest('th');
  if (!th) return;
  const v = state.view;
  if (v.sort === th.dataset.k) v.dir = -v.dir;
  else { v.sort = th.dataset.k; v.dir = 1; }
  v.page = 0;
  applyView();
});

$('sites').querySelector('tbody').addEventListener('click', (e) => {
  const tr = e.target.closest('tr');
  if (!tr) return;
  const r = state.view.rows[+tr.dataset.i];
  state.view.sel = r;
  renderPage();
  renderDetail(r);
  $('detail').scrollIntoView({ block: 'nearest', behavior: 'smooth' });
});
