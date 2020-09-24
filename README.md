# Pangenome Graph Variation Format (PGVF)

## Pangenomes

Pangenomes represent the total genomic sequence of a collection of genomes of the same species or clade.
Pangenome *models* allow us to understand the relationship between genomes and their component sequences in the context of the pangenome.
Pangenome *graphs* are representative graphical models of pangenomes.

While pangenome graphs let us represent differences between genomes, we have no mechanism to represent *differences between pangenome graphs*, or to combine multiple pangenome graphs into one structure without losing information.
This motivates the development of a new biological data format.

## Graph-to-graph mappings in PGVF

PGVF is a hard fork of the GFAv1 format that allows the description of graph-to-graph alignments.
It represents a collection of aligned graphs as a network of walks through an underlying merged sequence graph.
In principle, the structure of the underlying graph is unstable, while graphs mapped into it can always be extracted losslessly and are maintained across updates and modifications to the graph.
Subsets and projections of PGVF into GFAv1 are trivial, but precise compatibility is avoided to avoid confusion and eliminate features of GFA that cause problems for pangenome graphs (such as overlaps).

## Approach

Linear sequences are simple graphs, and their alignment to the graph is can be defined as a walk with edits.
When they are perfectly embedded in the underlying graph, a pure walk through the nodes of the graph is sufficient to describe them.
Their embedding makes their relationship to other embedded paths precise.
(We can already represent these in GFAv1.)

A sequence graph is a collection of sequences and links between their ends.
It an be represented within a new graph as a set of sequence walks (S-records) connected by a special type of graph edge (a L-record "link" in PGVF).
The addition of this additional link concept thus lets us express the alignment of graphs to each other by their mutual embedding in an underlying graph.
PGVF thus represents collections of graph-to-graph mappings and any simple structures, including sequence-to-graph and sequence-to-sequence relationships.

## Purpose

Graph-to-graph mapping can describe a number of concepts that are important in pangenomics.
A variable site is defined by set of walks representing its alleles, and genotypes over such alleles.
In pangenome graphs,  annotations on linear genomes (intervals of exons or protein binding motifs) are simply defined as walks, but genes and regions of interest must be defined as subgraphs.

Due to incomplete data (such as short reads), many genomes cannot be completely phased.
But, it may be possible to represent them as diploid assembly graphs.
We need graph-to-graph mapping to compare them to each other and to build pangenomes from collections of such assemblies.

## PGVF

We are aware of no format or existing schema (e.g. in RDF) that allows us to represent these relationships.
This specification defines such a format based on a collection of fixed record types with optional tags.

Our model describes variation between pangenome graphs, thus we call it the Pangenome Graph Variation Format, or PGVF.

## Grammar

```
<pgvf>     <- <line>+
<line>     <- <S-line> | <L-line> | <N-line> | <E-line> | <W-line> | <A-line> | <V-line> | <G-line> | <C-line>
<seq>      <- string
<segId>    <- string
<nodeId>   <- integer
<strand>   <- '>' | '<'
<walk>     <- (<strand>(<nodeId>|<segId>))+
<pathId>   <- string
<sampleId> <- string
<phaseId>  <- integer
<seqId>    <- string
<fragId>   <- integer
<varId>    <- string
<phaseId>  <- integer
<segPos>   <- integer
<walkPos>  <- integer
<cigar>    <- extended CS tag alignment encoding
<offset>   <- integer
<genotype> <- VCF-style genotype format
<comment>  <- string
<S-line>   <- 'S' <segId> ( <seq> | '*' ) ( <walk> | '*' ) ( <cigar> | '*' ) <tag>*
<L-line>   <- 'L' <strand><segId> <strand><segId> <tag>*
<N-line>   <- 'N' <nodeId> <seq> <tag>*
<E-line>   <- 'E' <strand><nodeId> <strand><nodeId> <tag>*
<W-line>   <- 'W' <sampleId> <phaseId> <seqId> <fragId> <walk> <tag>*
<A-line>   <- 'A' <segId> <walk> <cigar> <tag>*
<V-line>   <- 'V' <varId> <segId>+ <tag>*
<G-line>   <- 'G' <varId> <sampleId> <genotype> <tag>*
<C-line>   <- 'C' ( '*' | 'S' | 'L' | 'N' | 'E' | 'W' | 'A' | 'V' | 'G' ) <tagId> '"'<comment>'"'
```

## The base graph

The base graph is represented with `N` and `E` lines.
It's conceptually the same as GFA, just using different letters and the `>` and `<` characters to represent node strands.
The base graph should be stored in using a compact range of integer IDs.
These are not meant to be stable, and the sructure of nodes and edges can change with the addition or subtraction of sequences and graphs from the PGVF.

```
N 1 GAT
N 2 C
N 3 T
N 4 ACA
E >1 >2
E >2 >4
E >3 >4
```

## Overlay graphs

The overlay graphs are very similar to GFAv1, and use the same characters for the same concepts (S=sequences, L=links).
In PGVF, they represent stable coordinate and sequence systems projected into the base graph.

`S` records link sequences to walks through the graph.
They have optional cigars that allow the representation of variant alleles relative to the base graph
This facility allows us to establish overlay graphs containing rare or novel variation without disrupting the base graph. 
Either the sequence or the walk may be omitted depending on use, but the cigar requires the presence of the walk field.
If the cigar is omitted, the walk is assumed to be pure.

`L` records define links between the end of `S` records.
These let us embed the full graph topology of the overlay within our base graph.

### Embedding an overlay graph

This overlay is itself a valid graph, but it has not been mapped into a base graph.

```
S ctg1 GATTACA * *
S ctg2 GAT * *
S ctg3 T * *
S ctg4 C * *
S ctg5 ACA * *
S ctg6 ATGACA * *
L >ctg2 >ctg3
L >ctg2 >ctg4
L >ctg3 >ctg5
L >ctg4 >ctg5
```

We can represent the mapping of the overlay graph onto the base graph given above using this PGVF:

```
S ctg1 * >1>3>5 *
S ctg2 * >1 *
S ctg3 * >3 *
S ctg4 * >2 *
S ctg5 * >4 *
S ctg6 * >1>3>5 -1:3*tg:3
L >ctg2 >ctg3
L >ctg2 >ctg4
L >ctg3 >ctg5
L >ctg4 >ctg5
N 1 GAT
N 2 C
N 3 T
N 4 ACA
E >1 >2
E >2 >4
E >3 >4
```

Because the `S` records can be reconstructed by walking through the base graph, we can drop their sequences and represent them as walks.
But, we've left ctg6 out of the base graph, representing it with a walk and a cigar.
The particular format of the cigar allows us to indicate the starting position of the alignment on the target walk and the sequence of the `S`-line.

## CIGARs

Alignment descriptions are used in several parts of PGVF to represent differences from the base or overlay graphs.
They have a particular set of features that distinguish them from several other implementations.
In general, we follow the cs SAM/PAF tag format, beacuse it supports the lossless representation of the differences.
This lets us optionally drop the literal sequence of the query.
However, we extend the format slightly to encode the full alignment without additional offset information, and allow for deletions to be represented without their sequence, and SNPs without using the '*' character that indicates an empty field in PGVF.

The PGVF cigar matches the regular expression /(^.[0-9]+|,[0-9]+|:[0-9]+|\![a-z]|[-\+][a-z]+|=[A-Z]+)+/.
It consists of series of operations.
Each leading character specifies the operation, while the following integer or sequence defines the outcome of the operation.

This alignment:

```
ATGCGATCGATAAATAGAGTAG---GAATAGCA
   ||||||   ||||||||||   |||| |||
   CGATCG---AATAGAGTAGGTCGAATTGCAATT
```

... is represented as `_3:6_3:10+gtc:4!at:3:?3`

- `_3` indicates that we start at offset=3 (the fourth position) in the target sequence ("deletions" at the beginning of the target are treated as offsets)
- `:[0-9]+` represents an identical block
- `_3` represents a deltion
- `+gtc` an insertion
- `!at` indicates reference base A is substituted with a query base T
- `?3` indicates that we have a softclip of 3bp at the end of the mapping ("insertions" at the beginning and end of the query are treated as softclips)

Using `=`, `+`, and `-` rather than `:`, `?`, and `_`, the matching sequences as well as sequences for insertions and deletions may be represented.
If so, the previous example would become:

```
_3=CGATCG-ata=AATAGAGTAG+gtc=GAAT!at=GCA?3
```


## Describing genomes in the graph: walks, alignments, and variants

Pangenome graphs have unique applications that differ from assembly graphs.
If they serve as a reference data structure, we will use them to compare other samples to each other.
This encourages us to define variation across alleles, and phased haplotypes through those variable sites.
As with other features in PGVF, we represent these using walks through the graph.

### Walks (haplotypes)

A haplotype record is stored in W, which lets us collect together and order a set of haplotypes related to one individual.
Each `W`-line has a numeric index describing where it sits in the given phase set, allowing us to scaffold a collection of haplotypes together even if we cannot reconstruct them completely.
Here, phase 0 on chrz is unified, but the second phase breaks at the SNP and continues after.
We also represent a segment that differs from the base graph, which is used by the walk made by HG999.

```
S seg1 * >3 !tg
W HG002 0 chrz 0 >1>3>4
W HG002 1 chrz 0 >1
W HG002 1 chrz 1 >4
W HG999 1 chrz 0 >1>seg1>4
N 1 GAT
N 2 C
N 3 T
N 4 ACA
E >1 >2
E >2 >4
E >3 >4
```

### Alignments

Alignments (`A`-lines) provide a flexible way to describe the mapping between a sequence and the graph.
Unlike `S`-lines, they allow us to represent multiple alignments.
But, they are not meant to be referred to by other entities (like walks in `W`-lines or variable sites in `V`-lines).
The key objective of `A`-lines is to allow us to use a single input file in a pangenome construction pipeline.

`A`-lines are conceptually equivalent to the GAF format, but have fewer fields.
They must refer to sequences that are represented within the `S`-lines of the given PGVF.
In the following example, we show how to describe the mapping of `A`-lines.
Note that two mappings are shown for `f9b5e7d`.

```
S f9b5e7d TGACA * *
S 041eb7a TCACA * *
A f9b5e7d >1>3>4 _2:1!tg:3 as:Z:primary
A f9b5e7d >1 ?1:2?2 as:Z:secondary
A 041eb7a >1>2>4 _2:5 as:Z:primary
N 1 GAT
N 2 C
N 3 T
N 4 ACA
E >1 >2
E >2 >4
E >3 >4
```

Any PGVF file containing alignments can be converted into one with only `S`, `N`, and `E`-lines by embedding any novel sequences represented by the alignments in the graph so that they form pure walks through the `N` records.
For instance (using only the primary alignments shown above), the previous PGVF could be reduced to this one:

```
S f9b5e7d * >1>2>4>5
S 041eb7a * >2>4>5
N 1 GA
N 2 T
N 3 C
N 4 G
N 5 T
N 6 ACA
E >1 >2
E >2 >4
E >3 >4
```

### Variants and Genotypes

Variable sites collect mutually-exclusive alleles across a locus in the pangenome.
Optional genotype records provide a method to describe 
These genotypes are unphased.
The alleles are given as tuples of walks, an offset from the beginning of the first step, the offset to the end of the last step, and the alternative sequence of the allele.

```
S allele1 * >6 _1:1
S allele2 * >6 _1!cg
V var1 allele1,allele2 qual:f:89.3
G var1 sample1 0/1 gq:f:38.2 rd:i:23
G var1 sample2 1/1 gq:f:26.9 rd:i:18
G var1 sample3 0/0 gq:f:57.1 rd:i:37
```

In effect, the `V`-line can be used to describe bubble structures found in the graph.
However, there is no requirement that the `V`-line describe a formal bubble, ultrabbuble, or snarl.
It merely defines the allele scope for a set of genotype calls.

In the current form of this specification, phase blocks should be represented with `W`-lines, but future updates may allow `G`-lines to represent phase.

### Comments

Inline comments allow us to describe the meaning of tags for given line types.
For instance, this comment might explain the variant call set in the previous example:

```
C * * "a variant call set for samples 1, 2, and 3"
C V qual:f "PHRED-scaled variant call quality"
C G gq:f "PHRED-scaled genotype quality"
C G rd:i "median read depth over bases in the locus"
```

Because tags define their type, this is not necessary for parsing, and is optional.
Generic comments replace the line and field types with `*`.

## Notes / TODO

*This specification is a work in progress.*
Comments, suggestions, and improvements are welcome.

Provided consensus about its basic form, we will move to implement key required transformations via the `seqwish` graph induction kernel.
Basic projections and manipulations of the graph can be implemented in a particular "tools" library.

The representation of the PGVF depends heavily on walks, which we can represented efficiently in a dynamic form using `libbdsg` and other HandleGraph models, and in ultra-compact queryable form in the GBWT.

An initial application of this approach is to layer small variants on base "consensus" pangenome reference graphs built by `minigraph` and `smoothxg`.

The overlay graphs will allow us to use the hierarchical coordinate systems elaborated in rGFA while still representing dense, small variation.

Although graph-to-graph alignment has been used to generate multiple sequence alignments (in the PO-POA algorithm), no current method supports the heuristic alignment of large sequence graphs.
