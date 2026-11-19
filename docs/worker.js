// Runs the miRtherm core (WebAssembly) on chunks of UTRs. One instance per CPU thread.
importScripts('mirtherm_wasm.js');

let api = null;
let currentKey = '';

const ready = MirthermModule().then((m) => {
  api = {
    configure: m.cwrap('mt_configure', null, Array(13).fill('number')),
    header: m.cwrap('mt_header', 'string', []),
    // returns a pointer; 0 means the sequence could not be scanned
    scanPtr: m.cwrap('mt_scan', 'number', ['string', 'string', 'string']),
    str: m.UTF8ToString,
  };
});

function configure(c) {
  const key = JSON.stringify(c);
  if (key === currentKey) return;
  api.configure(c.at, c.gc, c.gt, c.mis, c.go, c.ge, c.temperature, c.energy, c.score, c.loop,
                c.noLP ? 1 : 0, c.noGU ? 1 : 0, c.noCloseGU ? 1 : 0);
  currentKey = key;
}

onmessage = async (e) => {
  const msg = e.data;
  await ready;
  configure(msg.cfg);
  if (msg.type === 'header') {
    postMessage({ id: msg.id, header: api.header() });
    return;
  }
  // msg.type === 'scan': one output block per UTR, null if it could not be scanned
  const out = msg.utrs.map((seq) => {
    let block = '';
    for (const mir of msg.mirs) {
      const ptr = api.scanPtr(mir.id, mir.seq, seq);
      if (!ptr) return null;
      block += api.str(ptr);
    }
    return block;
  });
  postMessage({ id: msg.id, out });
};
