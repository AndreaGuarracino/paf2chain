# paf2chain
convert PAF format to CHAIN format

## installation

`paf2chain` is built with rust, and so we install using `cargo`:

```
https://github.com/AndreaGuarracino/paf2chain
cd paf2chain
cargo install --force --path .
```

## usage

Generate alignments with the cigar string attached to the `cg:Z:` tag.
These can be made by several aligners, including `minimap2 -c`, `wfmash`, or `lastz --format=paf:wfmash`.
With alignments in `aln.paf`, we would convert it into a CHAIN format file using this call:

```
paf2chain -i aln.paf > aln.chain
```
