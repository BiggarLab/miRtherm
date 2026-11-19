# miRtherm

microRNA target prediction with duplex free energy calculated at a chosen
temperature. There are two optional user interfaces:

* a browser dashboard (`docs/`) that runs entirely on the user's computer
  through WebAssembly
* a command line program (`mirtherm`)

## Dashboard

Enter one or more mature miRNA sequences, choose a FASTA file of target
sequences (3' UTRs), set the temperatures on the scale and press "Find
targets". Sequences are processed locally. The page shows
how many transcripts are targeted at each temperature, the distribution of
duplex energies, a sortable site table and the duplex at each temperature.
Each temperature can be downloaded in the command line output format, and the
site table as CSV.

miRNA and target sequences can be pasted or uploaded as FASTA files. With
several miRNAs the page adds a summary table, one row per miRNA, that fills in
as each one finishes; click a row to see that miRNA's sites. All output files
can be downloaded together as a zip. Batches run one miRNA after another, and
a batch of several hundred is better suited to the command line program.

To try it locally: `cd docs && python3 -m http.server`, then open
http://localhost:8000. Opening `index.html` straight from disk does not work
because browsers block Web Workers on `file://` pages.

## Command line

    make
    ./mirtherm mirna.fa utrs.fasta -energy -20 -T 25 -o out.txt

Options: `-at -gc -gt -mis -go -ge` (alignment weights), `-energy -score
-loop` (filters), `-T` (temperature in °C, decimals allowed), `-noLP -noGU
-noCloseGU`, `-tab` or `-text`, `-o file`.

Output columns (tab format): miRNA, alignment score, duplex energy, start,
end, sequence alignment, energy alignment, seed score, seed energy, loop
score, dot-bracket structure.

## How a site is scored

1. Seed scan. 6-mer matches to miRNA nucleotides 1-6, 2-7 or 3-8, allowing one
   G:U wobble.
2. Alignment. The rest of the miRNA is aligned (affine gaps) to the 40 nt
   upstream of the seed match. This gives the score.
3. Duplex energy. The miRNA and the aligned target stretch are folded together
   with the ViennaRNA 1.5 `cofold` routine at the requested temperature. If the
   miRNA 3' end overhangs the target, the target is extended and refolded.
   Pairs within the target count toward the energy.
4. Loop score. Position and size of the central loop.
5. Filters on energy, score and loop score.

## Testing

`tests/compare.sh` runs the program and compares the result with a reference
output file, byte for byte:

    tests/compare.sh mirna.fa utrs.fasta reference.txt -energy -20 -T 25

The command line and WebAssembly builds give identical output.

## Building the WebAssembly files

    make wasm        # needs Emscripten; writes docs/mirtherm_wasm.js and .wasm

## Credits and licence

Duplex folding uses a subset of the ViennaRNA Package 1.5 (Ivo Hofacker and
colleagues, Institute for Theoretical Chemistry, University of Vienna) in
`core/vienna`, unchanged apart from build fixes. Its licence
(`core/vienna/COPYING.ViennaRNA`) allows research, educational and commercial
use and modification, provided the package and derived works are not
redistributed for a fee and credit is given. That condition applies to this
repository as a whole.

