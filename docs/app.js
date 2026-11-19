'use strict';
// miRtherm dashboard. All computation happens in Web Workers running the
// WebAssembly build of the core; nothing leaves the browser.

const $ = (id) => document.getElementById(id);
const PAGE = 40;
const SCALE_MIN = -5, SCALE_MAX = 50;
const COLD = [74, 79, 196], MID = [184, 65, 138], WARM = [224, 112, 31];

const state = {
  input: { mir: null, utr: null },   // parsed records when a file is loaded instead of pasted text
  run: null,     // {temps, cutoff, utrs, mirs: [{id, seq, name, status, stats, files}]}
  detail: null,  // {mi, rows} for the miRNA currently shown
  view: { rows: [], page: 0, sort: 'e0', dir: 1, sel: null, bin: null },
  workers: [],
  cancelled: false,
};
const PASTE_LIMIT = 200000;   // larger files are kept out of the text box

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

function mirRecords(text) {
  let mirs;
  if (text.includes('>')) mirs = parseFasta(text);
  else mirs = text.split(/\s+/).filter(Boolean).map((q, i) => ({ id: 'miRNA-' + (i + 1), seq: q.toUpperCase() }));
  for (const m of mirs) m.seq = m.seq.replace(/\s+/g, '');
  return mirs.filter((m) => m.seq.length);
}

function checkMirs(mirs) {
  if (!mirs.length) throw new Error('Enter at least one miRNA sequence.');
  for (const m of mirs) {
    if (!/^[ACGTU]+$/.test(m.seq)) throw new Error(`miRNA “${m.id}” has characters other than A, C, G, T or U.`);
    if (m.seq.length < 9 || m.seq.length > 56) throw new Error(`miRNA “${m.id}” is ${m.seq.length} nt long. miRNAs must be 9 to 56 nt.`);
  }
  return mirs;
}

function utrRecords(text) {
  const utrs = parseFasta(text);
  for (const u of utrs) {
    u.seq = u.seq.replace(/\s+/g, '');
    u.name = u.id.split(/\s+/)[0] || u.id;
    const g = /(?:parent|gene|gene_id|gene_symbol)[=:]([^;\s,]+)/i.exec(u.id);
    u.gene = g ? g[1] : '';
  }
  return utrs;
}

const getMirs = () => checkMirs(state.input.mir ? state.input.mir.records : mirRecords($('mirs').value));
function getUtrs() {
  const utrs = state.input.utr ? state.input.utr.records : utrRecords($('utrs').value);
  if (!utrs.length) throw new Error('Paste target sequences in FASTA format, or upload a file.');
  return utrs;
}

// Each sequence input is a paste box that turns into a file summary when a
// large file is loaded.
function showInput(kind) {
  const f = state.input[kind];
  const box = $(kind === 'mir' ? 'mirs' : 'utrs');
  $(kind + 'chip').hidden = !f;
  box.hidden = !!f;
  const info = $(kind + 'info');
  info.classList.remove('bad');
  let recs;
  if (f) {
    $(kind + 'chiptext').textContent = f.label;
    recs = f.records;
  } else {
    recs = kind === 'mir' ? mirRecords(box.value) : utrRecords(box.value);
  }
  if (!recs.length) { info.textContent = ''; return; }
  const nt = recs.reduce((a, r) => a + r.seq.length, 0);
  info.textContent = `${recs.length.toLocaleString()} sequence${recs.length === 1 ? '' : 's'}` +
    (kind === 'utr' ? `, ${nt >= 1e6 ? (nt / 1e6).toFixed(2) + ' Mb' : nt.toLocaleString() + ' nt'}` : '');
  if (kind === 'mir' && recs.length > 100) info.textContent += '. A batch this size can take a long time in the browser.';
}

function loadText(kind, text, label) {
  const recs = kind === 'mir' ? mirRecords(text) : utrRecords(text);
  if (!recs.length) throw new Error(`No FASTA records found in ${label}.`);
  if (text.length <= PASTE_LIMIT) {
    state.input[kind] = null;
    $(kind === 'mir' ? 'mirs' : 'utrs').value = text;
  } else {
    state.input[kind] = { label, records: recs };
  }
  showInput(kind);
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
    // Stop everything; calls still waiting for a worker fail with this error.
    abort() {
      this.close();
      const err = new Error('Stopped.');
      for (const j of pending.values()) j.reject(err);
      for (const j of queue) j.reject(err);
      pending.clear();
      queue.length = 0;
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
  const mirList = getMirs();
  const utrs = getUtrs();
  const { temps, base, threads } = readSettings();
  const names = new Set();
  const mirs = mirList.map((m) => {
    let name = (m.id.split(/\s+/)[0] || 'miRNA').replace(/[^\w.+-]+/g, '_'), k = 2;
    const stem = name;
    while (names.has(name)) name = stem + '_' + k++;
    names.add(name);
    return { id: m.id, seq: m.seq, name, status: 'queued', stats: null, files: [] };
  });
  state.run = { temps, cutoff: base.energy, utrs, mirs };
  state.detail = null;
  state.cancelled = false;
  $('results').hidden = false;
  $('one').hidden = true;
  renderBatch();

  const pool = makePool(threads);
  state.pool = pool;
  const chunks = chunkUtrs(utrs);
  const total = chunks.length * temps.length * mirs.length;
  let done = 0, skippedAll = 0;
  const t0 = performance.now();
  progress(0, 'Starting');
  try {
    for (let mi = 0; mi < mirs.length; mi++) {
      const M = mirs[mi];
      M.status = 'running';
      renderBatch();
      const blocksByTemp = [];
      const skipped = new Set();
      for (const T of temps) {
        const cfg = Object.assign({ temperature: T }, base);
        const header = (await pool.post({ type: 'header', cfg })).header;
        const blocks = new Array(utrs.length);
        await Promise.all(chunks.map(async ([a, b]) => {
          const res = await pool.post({ type: 'scan', cfg, mirs: [{ id: M.id, seq: M.seq }], utrs: utrs.slice(a, b).map((u) => u.seq) });
          for (let i = a; i < b; i++) {
            if (res.out[i - a] === null) skipped.add(i);
            blocks[i] = res.out[i - a] || '';
          }
          done++;
          const left = (performance.now() - t0) / 1000 / done * (total - done);
          const eta = left > 90 ? `about ${Math.ceil(left / 60)} min left` : `about ${Math.ceil(left)} s left`;
          progress(done / total, (mirs.length > 1 ? `miRNA ${mi + 1} of ${mirs.length}, ` : '') + `folding at ${T} °C, ${eta}`);
        }));
        blocksByTemp.push(blocks);
        // Same layout as the command line output file.
        M.files.push(await packFile(`${M.name}-${T}C.txt`, header + utrs.map((u, i) => u.id + '\n' + blocks[i]).join('\n')));
      }
      M.stats = statsOf(blocksByTemp);
      M.status = 'done';
      skippedAll += skipped.size;
      renderBatch();
      if (!state.detail) openMir(mi, false);
    }
  } finally {
    pool.close();
    for (const M of mirs) if (M.status !== 'done') M.status = state.cancelled ? 'stopped' : 'failed';
    renderBatch();
  }
  const secs = (performance.now() - t0) / 1000;
  progress(1, `Done in ${secs > 90 ? (secs / 60).toFixed(1) + ' min' : secs.toFixed(1) + ' s'} on ${threads} thread${threads > 1 ? 's' : ''}.` +
    (skippedAll ? ` Skipped ${skippedAll} target sequence${skippedAll > 1 ? 's' : ''} with characters other than A, C, G, T or U.` : ''));
}

// Counts for the summary table, taken straight from the per-UTR result blocks.
function statsOf(blocksByTemp) {
  const n = blocksByTemp.length;
  const sites = [], tx = [];
  for (const blocks of blocksByTemp) {
    let s = 0, t = 0;
    for (const blk of blocks) if (blk) { t++; for (let i = 0; i < blk.length; i++) if (blk.charCodeAt(i) === 10) s++; }
    sites.push(s);
    tx.push(t);
  }
  let lost = 0, kept = 0, gained = 0;
  if (n >= 2) {
    const A = blocksByTemp[0], B = blocksByTemp[n - 1];
    for (let i = 0; i < A.length; i++) {
      if (A[i] && B[i]) kept++;
      else if (A[i]) lost++;
      else if (B[i]) gained++;
    }
  }
  return { sites, tx, lost, kept, gained };
}

// ---------- result files (kept compressed) and zip ----------

const CRC_TABLE = (() => {
  const t = new Uint32Array(256);
  for (let n = 0; n < 256; n++) {
    let c = n;
    for (let k = 0; k < 8; k++) c = c & 1 ? 0xedb88320 ^ (c >>> 1) : c >>> 1;
    t[n] = c >>> 0;
  }
  return t;
})();
function crc32(bytes) {
  let c = 0xffffffff;
  for (let i = 0; i < bytes.length; i++) c = CRC_TABLE[(c ^ bytes[i]) & 0xff] ^ (c >>> 8);
  return (c ^ 0xffffffff) >>> 0;
}

async function packFile(name, text) {
  const bytes = new TextEncoder().encode(text);
  const file = { name, crc: crc32(bytes), size: bytes.length, method: 0, blob: null };
  if (typeof CompressionStream === 'function') {
    try {
      file.blob = await new Response(new Blob([bytes]).stream().pipeThrough(new CompressionStream('deflate-raw'))).blob();
      file.method = 8;
    } catch (e) { file.blob = null; }
  }
  if (!file.blob) file.blob = new Blob([bytes]);
  return file;
}

function unpackFile(file) {
  if (file.method === 0) return file.blob.text();
  return new Response(file.blob.stream().pipeThrough(new DecompressionStream('deflate-raw'))).text();
}

function zipBlob(files) {
  const enc = new TextEncoder();
  const d = new Date();
  const time = (d.getHours() << 11) | (d.getMinutes() << 5) | (d.getSeconds() >> 1);
  const date = ((d.getFullYear() - 1980) << 9) | ((d.getMonth() + 1) << 5) | d.getDate();
  const parts = [], central = [];
  let offset = 0;
  for (const f of files) {
    const name = enc.encode(f.name);
    const lh = new DataView(new ArrayBuffer(30));
    lh.setUint32(0, 0x04034b50, true); lh.setUint16(4, 20, true); lh.setUint16(6, 0x0800, true);
    lh.setUint16(8, f.method, true); lh.setUint16(10, time, true); lh.setUint16(12, date, true);
    lh.setUint32(14, f.crc, true); lh.setUint32(18, f.blob.size, true); lh.setUint32(22, f.size, true);
    lh.setUint16(26, name.length, true); lh.setUint16(28, 0, true);
    parts.push(lh.buffer, name, f.blob);
    const ch = new DataView(new ArrayBuffer(46));
    ch.setUint32(0, 0x02014b50, true); ch.setUint16(4, 20, true); ch.setUint16(6, 20, true); ch.setUint16(8, 0x0800, true);
    ch.setUint16(10, f.method, true); ch.setUint16(12, time, true); ch.setUint16(14, date, true);
    ch.setUint32(16, f.crc, true); ch.setUint32(20, f.blob.size, true); ch.setUint32(24, f.size, true);
    ch.setUint16(28, name.length, true); ch.setUint32(42, offset, true);
    central.push(ch.buffer, name);
    offset += 30 + name.length + f.blob.size;
  }
  const cdSize = central.reduce((a, p) => a + p.byteLength, 0);
  if (files.length > 65535 || offset + cdSize > 0xffffffff) throw new Error('The results are too large for one zip file. Download them per miRNA instead.');
  const end = new DataView(new ArrayBuffer(22));
  end.setUint32(0, 0x06054b50, true); end.setUint16(8, files.length, true); end.setUint16(10, files.length, true);
  end.setUint32(12, cdSize, true); end.setUint32(16, offset, true);
  return new Blob([...parts, ...central, end.buffer], { type: 'application/zip' });
}

// Read an output file back into site rows. `map` collects one row per site
// (transcript and end position) with a slot per temperature.
function parseOutput(text, ti, map) {
  const lines = text.split('\n');
  let i = 0;
  while (i < lines.length && lines[i][0] === '#') i++;
  let ui = -1, expectId = true;
  for (; i < lines.length; i++) {
    const line = lines[i];
    if (expectId) { ui++; expectId = false; continue; }
    if (line === '') { expectId = true; continue; }
    const f = line.split('\t');
    const key = ui + '|' + f[4];
    let row = map.get(key);
    if (!row) {
      row = { u: ui, mir: f[0], start: +f[3], end: +f[4], score: +f[1], dyn: f[5], seedScore: +f[7], at: [] };
      map.set(key, row);
    }
    row.at[ti] = { energy: +f[2], aln: f[6], seedEnergy: +f[8], loop: f[9], structure: f[10] };
  }
}

async function openMir(mi, scroll) {
  const M = state.run.mirs[mi];
  if (!M || M.status !== 'done') return;
  const map = new Map();
  for (let ti = 0; ti < M.files.length; ti++) parseOutput(await unpackFile(M.files[ti]), ti, map);
  state.detail = { mi, rows: [...map.values()] };
  renderResults(scroll);
  renderBatch();
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

function summaryCols() {
  const { temps } = state.run;
  const cols = [{ label: 'miRNA', text: true, get: (M) => M.name }];
  temps.forEach((T, i) => cols.push({ label: `Transcripts at ${deg(T)}`, get: (M) => M.stats.tx[i] }));
  temps.forEach((T, i) => cols.push({ label: `Sites at ${deg(T)}`, get: (M) => M.stats.sites[i] }));
  if (temps.length >= 2) {
    cols.push({ label: 'Drop out on warming', get: (M) => M.stats.lost });
    cols.push({ label: 'Stay', get: (M) => M.stats.kept });
    cols.push({ label: 'Appear', get: (M) => M.stats.gained });
  }
  return cols;
}

// Summary table, one row per miRNA. Only shown for batches.
function renderBatch() {
  const { mirs } = state.run;
  const sec = $('batch');
  sec.hidden = mirs.length < 2;
  if (sec.hidden) return;
  const cols = summaryCols();
  const tbl = $('batchtable');
  tbl.querySelector('thead').innerHTML = '<tr>' + cols.map((c) => `<th scope="col"${c.text ? ' class="t"' : ''}>${esc(c.label)}</th>`).join('') + '</tr>';
  const words = { queued: 'waiting', running: 'running', stopped: 'stopped', failed: 'failed' };
  tbl.querySelector('tbody').innerHTML = mirs.map((M, i) => {
    if (M.status !== 'done')
      return `<tr class="wait"><td class="t">${esc(M.name)}</td><td class="t" colspan="${cols.length - 1}">${words[M.status]}</td></tr>`;
    return `<tr data-mi="${i}"${state.detail && state.detail.mi === i ? ' class="sel"' : ''}>` +
      cols.map((c) => `<td${c.text ? ' class="t"' : ''}>${c.text ? esc(c.get(M)) : c.get(M).toLocaleString()}</td>`).join('') + '</tr>';
  }).join('');
  const nd = mirs.filter((M) => M.status === 'done').length;
  $('batchcount').textContent = `${nd} of ${mirs.length} finished`;
  $('zip').disabled = $('sumcsv').disabled = nd === 0;
}

function renderResults(scroll) {
  const { temps, utrs, mirs } = state.run;
  const rows = state.detail.rows;
  $('results').hidden = false;
  $('one').hidden = false;
  $('onetitle').textContent = mirs.length > 1 ? mirs[state.detail.mi].id : '';
  $('onetitle').hidden = mirs.length < 2;
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
  if (scroll !== false) $('one').scrollIntoView({ block: 'start', behavior: 'smooth' });
  else if (mirs.length < 2) $('results').scrollIntoView({ block: 'start', behavior: 'smooth' });
}

// Small multiples: one horizontal histogram of site energies per temperature,
// sharing both axes so lanes can be compared directly.
function drawLanes() {
  const { temps } = state.run;
  const rows = state.detail.rows;
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
  const { temps, utrs } = state.run;
  const rows = state.detail.rows, v = state.view;
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
    if (q && !utrs[r.u].id.toLowerCase().includes(q)) return false;
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
  const n = temps.length, u = state.run.utrs[r.u];
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
  downloadBlob(name, new Blob([text], { type: type || 'text/plain' }));
}

function downloadBlob(name, blob) {
  const url = URL.createObjectURL(blob);
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
    const u = state.run.utrs[r.u];
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

for (const kind of ['mir', 'utr']) {
  $(kind + 'file').addEventListener('change', async (e) => {
    const f = e.target.files[0];
    e.target.value = '';
    if (!f) return;
    $('error').hidden = true;
    try { loadText(kind, await f.text(), f.name); } catch (err) { fail(err); }
  });
  $(kind + 'clear').addEventListener('click', () => { state.input[kind] = null; showInput(kind); });
  $(kind === 'mir' ? 'mirs' : 'utrs').addEventListener('input', () => showInput(kind));
  showInput(kind);
}

$('example').addEventListener('click', async () => {
  $('error').hidden = true;
  try {
    const res = await fetch('example_utrs.fasta');
    if (!res.ok) throw new Error('Could not load the example sequences.');
    const text = await res.text();
    state.input.utr = { label: 'Fly example, 400 3′ UTRs', records: utrRecords(text) };
    showInput('utr');
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
  if (state.pool) state.pool.abort();
  progress(0, 'Stopped. Finished miRNAs are kept.');
});

$('tally').addEventListener('click', async (e) => {
  const ti = e.target.dataset.dl;
  if (ti === undefined) return;
  const f = state.run.mirs[state.detail.mi].files[+ti];
  download(f.name, await unpackFile(f));
});

$('batchtable').querySelector('tbody').addEventListener('click', (e) => {
  const tr = e.target.closest('tr[data-mi]');
  if (tr) openMir(+tr.dataset.mi, true);
});

function summaryCsv() {
  const cols = summaryCols();
  const q = (v) => '"' + String(v).replace(/"/g, '""') + '"';
  const lines = [['mirna_id', 'sequence', ...cols.slice(1).map((c) => c.label)].map(q).join(',')];
  for (const M of state.run.mirs) if (M.status === 'done')
    lines.push([q(M.id), M.seq, ...cols.slice(1).map((c) => c.get(M))].join(','));
  return lines.join('\n') + '\n';
}

$('sumcsv').addEventListener('click', () => download('mirtherm-summary.csv', summaryCsv(), 'text/csv'));
$('zip').addEventListener('click', async () => {
  try {
    const files = state.run.mirs.filter((M) => M.status === 'done').flatMap((M) => M.files);
    files.push(await packFile('summary.csv', summaryCsv()));
    downloadBlob('mirtherm-results.zip', zipBlob(files));
  } catch (err) { fail(err); }
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
