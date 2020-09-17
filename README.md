# pangenome graph variation format (PGVF)

## pangenomes

Pangenomes represent the total genomic sequence of a collection of genomes of the same species or clade.
Pangenome *models* allow us to understand the relationship between genomes and their component sequences in the context of the pangenome.
Pangenome *graphs* are representative graphical models of pangenomes.

While pangenome graphs let us represent differences between genomes, we have no mechanism to represent *differences between pangenome graphs*, or to combine multiple pangenome graphs into one structure without losing information.
This motivates the development of a new biological data format.

## graph-to-graph mappings in PGVF

PGVF is a hard fork of the GFAv1 format that allows the description of graph-to-graph alignments.
It represents a collection of aligned graphs as a network of walks through an underlying merged sequence graph.
In principle, the structure of the underlying graph is unstable, while graphs mapped into it can always be extracted losslessly and are maintained across updates and modifications to the graph.
Subsets and projections of PGVF into GFAv1 are trivial, but precise compatibility is avoided to avoid confusion and eliminate features of GFA that cause problems for pangenome graphs (such as overlaps).

## approach

Linear sequences are simple graphs, and their alignment to the graph is can be defined as a walk with edits.
When they are perfectly embedded in the underlying graph, a pure walk through the nodes of the graph is sufficient to describe them.
Their embedding makes their relationship to other embedded paths precise.
(We can already represent these in GFAv1.)

A sequence graph is a collection of sequences and links between their ends.
It an be represented within a new graph as a set of sequence walks (S-records) connected by a special type of graph edge (a L-record "link" in PGVF).
The addition of this additional link concept thus lets us express the alignment of graphs to each other by their mutual embedding in an underlying graph.
PGVF thus represents collections of graph-to-graph mappings and any simple structures, including sequence-to-graph and sequence-to-sequence relationships.

## purpose

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

## grammar

```
<pgvf>     <- <line>+
<line>     <- <S-line> | <L-line> | <N-line> | <E-line> | <W-line> | <A-line> | <V-line> | <G-line> | <C-line>
<seq>      <- string
<segId>    <- string
<nodeId>   <- integer
<strand>   <- '>' | '<'
<walk>     <- (<strand><nodeId>)+
<pathId>   <- string
<sampleId> <- string
<hapId>    <- string
<varId>    <- string
<phaseId>  <- integer
<segPos>   <- integer
<walkPos>  <- integer
<cigar>    <- extended cigar
<offset>   <- integer
<allele>   <- <walk>':'<offset>':'<offset>':'<seq>
<genotype> <- VCF-style genotype format
<comment>  <- string
<S-line>   <- 'S' <segId> ( <seq> | '*' ) ( <walk> | '*' ) <tag>*
<L-line>   <- 'L' <strand><segId> <strand><segId> <tag>*
<N-line>   <- 'N' <nodeId> <seq> <tag>*
<E-line>   <- 'E' <strand><nodeId> <strand>nodeId> <tag>*
<W-line>   <- 'W' <sampleId> <hapId> <phaseId> <walk> <tag>*
<A-line>   <- 'A' <segId> <segPos> <strand> <walkpos> <walk> <cigar> <tag>*
<V-line>   <- 'V' <varId> <allele>+ (<sampleId>:<genotype>)* <tag>*
<G-line>   <- 'G' <varId> <sampleId> <genotype> <tag>*
<C-line>   <- 'C' ( '*' | 'S' | 'L' | 'N' | 'E' | 'W' | 'A' | 'V' | 'G' ) <tagId> '"'<comment>'"'
```

## example

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

The overlay graphs are nearly equivalent to GFA, using the same characters.
They represent stable coordinate and sequence systems in the PGVF.

`S` records link sequences to walks through the graph.
Either the sequence or the walk may be omitted depending on use.
This overlay has not been mapped into the graph.

```
S ctg1 GATTACA *
S ctg2 GAT *
S ctg3 T *
S ctg4 C *
S ctg5 ACA *
L >ctg2 >ctg3
L >ctg2 >ctg4
L >ctg3 >ctg5
L >ctg4 >ctg5
```

We can represent the mapping of the overlay graph onto the first base graph in PGVF:

```
S ctg1 * >1>3>5
S ctg2 * >1
S ctg3 * >3
S ctg4 * >2
S ctg5 * >4
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

## walks, alignments, and variants

Pangenome graphs have unique applications that differ from assembly graphs.
If they serve as a reference data structure, we will use them to compare other samples to each other.
This encourages us to define variation across alleles, and phased haplotypes through those variable sites.
As with other features in PGVF, we represent these using walks through the graph.

### walks (haplotype records)

A haplotype record is stored in W, which lets us collect together and order a set of haplotypes related to one individual.
Each W line has a numeric index describing where it sits in the given haplotype, allowing us to scaffold a collection of haplotypes together even if we cannot reconstruct them completely.
Here, hap1 is unified, but hap2 breaks at the SNP and continues after.

```
W HG002 hap1 1 >1>3>4
W HG002 hap2 1 >1
W HG002 hap2 2 >4
N 1 GAT
N 2 C
N 3 T
N 4 ACA
E >1 >2
E >2 >4
E >3 >4
```

### alignments

Alignments (`A` lines) allow us to describe the mapping between a sequence and part of the graph.
They are conceptually equivalent to the GAF format, but have fewer fields.
They refer to sequences that are represented within the `S` lines of the given PGVF.

```
S f9b5e7d TGACA *
S 041eb7a TCACA *
A f9b5e7d 0 + 2 >1>3>4 1=1X3=
A 041eb7a 2 + 1 >1>2>4 5=
N 1 GAT
N 2 C
N 3 T
N 4 ACA
E >1 >2
E >2 >4
E >3 >4
```

Any PGVF file containing alignments can be converted into one with only `S`, `N`, and `E` lines by embedding any novel sequences represented by the alignments in the graph so that they form pure walks through the `N` records.
For instance, the previous PGVF could be reduced to this one:

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

### variants and genotypes

Variable sites collect mutually-exclusive alleles across a locus in the pangenome.
Optional genotype records provide a method to describe 
These genotypes are unphased.
(Phase blocks should be represented with `W` lines.)
The alleles are given as tuples of walks, an offset from the beginning of the first step, the offset to the end of the last step, and the alternative sequence of the allele.

```
V var1 >6:1:1:C,>6:1:1:G qual:f:89.3
G var1 sample1 0/1 gq:f:38.2 rd:i:23
G var1 sample2 1/1 gq:f:26.9 rd:i:18
G var1 sample3 0/0 gq:f:57.1 rd:i:37
```

### comments

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

## notes / todo

*This specification is a work in progress.*
Comments, suggestions, and improvements are welcome.

Provided consensus about its basic form, we will move to implement key required transformations via the `seqwish` graph induction kernel.
Basic projections and manipulations of the graph can be implemented in a particular "tools" library.

The representation of the PGVF depends heavily on walks, which we can represented efficiently in a dynamic form using `libbdsg` and other HandleGraph models, and in ultra-compact queryable form in the GBWT.

An initial application of this approach is to layer small variants on base "consensus" pangenome reference graphs built by `minigraph` and `smoothxg`.

The overlay graphs will allow us to use the hierarchical coordinate systems elaborated in rGFA while still representing dense, small variation.

Although graph-to-graph alignment has been used to generate multiple sequence alignments (in the PO-POA algorithm), no current method supports the heuristic alignment of large sequence graphs.
