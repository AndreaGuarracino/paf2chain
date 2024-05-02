# paf2chain
convert PAF format to CHAIN format

## citation
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8108447.svg)](https://doi.org/10.5281/zenodo.8108447)

## installation

`paf2chain` is built with rust, and so we install using `cargo`:

```
git clone https://github.com/AndreaGuarracino/paf2chain
cd paf2chain
cargo install --force --path .
```

## usage

Generate alignments with the cigar string attached to the `cg:Z:` tag.
These can be made by several aligners, including `minimap2 -c`, `wfmash`, or `lastz --format=paf:wfmash`.
With alignments in `aln.paf` (or `aln.paf.gz` if gzip-compressed), we would convert it into a CHAIN format file using this call:

```
paf2chain -i aln.paf > aln.chain
```

## info

`paf2chain` performs the reverse operation of [chain2paf](https://github.com/AndreaGuarracino/chain2paf).
