bigly: a pileup library that embraces the huge
==============================================

[![Build Status](https://travis-ci.org/brentp/bigly.svg?branch=master)](https://travis-ci.org/brentp/bigly)
[![GoDoc] (https://godoc.org/github.com/brentp/bigly?status.png)](https://godoc.org/github.com/brentp/bigly)

`bigly` is an API and a command-line app ([binaries here](https://github.com/brentp/bigly/releases)). It is similar to `samtools mpileup` but it reports a *huge* number of
additional variables that are useful for structural variant calling and visualization.

Each for each requested position, the struct below is filled with by the appropriate position in any overlapping
alignment that meets the requested filters:

```Go
// Pile holds the information about a single base.
type Pile struct {
    Chrom                  string
    Pos                    int
    Depth                  int    // count of reads passing filters.
    RefBase                byte   // Reference base a this position.
    MisMatches             uint32 // number of mismatches .
    ProperPairs            int    // count of reads with paired flag
    SoftStarts             uint32 // counts of base preceding an 'S' cigar op
    SoftEnds               uint32 // ...  following ...
    HardStarts             uint32 // counts of base preceding an 'H' cigar op
    HardEnds               uint32
    InsertionStarts        uint32 // counts of base preceding an 'I' cigar op
    InsertionEnds          uint32
    Deletions              uint32  // counts of deletions 'D' at this base
    Heads                  uint32  // counts of starts of reads at this base
    Tails                  uint32  // counts of ends of reads at this base
    Splitters              uint32  // count of non-secondary reads with SA tags.
    Splitters1             uint32  // count of non-secondary reads with exactly 1 SA tag.
    Bases                  []byte  // All bases from reads covering this position
    Quals                  []uint8 // All quals from reads covering this position
    MeanInsertSizeLP       uint32  // Calculated with left-most of pair
    MeanInsertSizeRM       uint32  // Calculated with right-most of pair
    WeirdCount             uint32  // Calculated with right-most of pair
    Discordant             uint32  // Number of reads with insert size > ConcordantCutoff
    DiscordantChrom        uint32  // Number of reads mapping on different chroms
    DiscordantChromEntropy float32 // high value means all discordants came from same chrom.
    GC65                   uint32  // GC content in a 65 base window centered on the current base.
    GC257                  uint32  
    Duplicity65            float32 // measure of lack of sequence entropy.
    Duplicity257           float32 // measure of lack of sequence entropy.
    SplitterPositions      []int   // Start, End Positions of splitters for reads overlapping this base.
    SplitterStrings        []string
}

```

The program in `cmd/bigly/main.go` is distributed as an example program of what one can do with this
library--namely make an enhanced pileup in a few lines of code.

Usage
-----

At this time, the usage of the example program is very simple.
Default exclude flags are `(sam.Unmapped | sam.QCFail | sam.Duplicate)`

```
bigly $bam $chrom:$start-$end > o
```

if a reference is specified with `-r` it will report statistics about GC content in windows
surrounding each base.

help:
```
bigly 0.2.0
usage: bigly [--minbasequality MINBASEQUALITY] [--minmappingquality MINMAPPINGQUALITY] [--excludeflag EXCLUDEFLAG] [--includeflag INCLUDEFLAG] [--mincliplength MINCLIPLENGTH] [--includebases] [--splitterverbosity SPLITTERVERBOSITY] [--reference REFERENCE] BAMPATH REGION

positional arguments:
  bampath
  region

options:
  --minbasequality MINBASEQUALITY, -q MINBASEQUALITY
                         base quality threshold [default: 10]
  --minmappingquality MINMAPPINGQUALITY, -Q MINMAPPINGQUALITY
                         mapping quality threshold [default: 5]
  --excludeflag EXCLUDEFLAG, -F EXCLUDEFLAG
  --includeflag INCLUDEFLAG, -f INCLUDEFLAG
  --mincliplength MINCLIPLENGTH, -c MINCLIPLENGTH
                         only count H/S clips of at least this length [default: 15]
  --includebases, -b     output each base and base quality score
  --splitterverbosity SPLITTERVERBOSITY, -s SPLITTERVERBOSITY
                         0-only count; 1:count and single most frequent; 2:all SAs; 3:dont shorten positions
  --reference REFERENCE, -r REFERENCE
                         optional path to reference fasta.
  --help, -h             display this help and exit
  --version              display version and exit

```

API
---

GoDoc Documentation is here: [![GoDoc] (https://godoc.org/github.com/brentp/bigly?status.png)](https://godoc.org/github.com/brentp/bigly)

With a supporting library easing regional bam queries here: [![GoDoc] (https://godoc.org/github.com/brentp/bigly/bamat?status.png)](https://godoc.org/github.com/brentp/bigly/bamat)

View the tests for more examples. `bigly` uses [biogo](https://github.com/biogo/hts) library for bam access.

Plotting
--------

```
python scripts/plotter.py bigly.output
```

will make a plot with python+matplotlib on output from bigly.

An example image looks like: 
![Example](https://cloud.githubusercontent.com/assets/1739/20151721/7a23b46a-a678-11e6-87be-0d4666faffdd.png "Bigly Deletion")

This is a homozygous deletion easily seen from the depth, but we can see that the deletion is nicely
delineated by soft-clips and we see aberrant insert sizes bounding the deletion.
In most libraries, we would also see splitters flanking the region.


TODO
----

+ Track strand of bases.
+ Report the 5th and 95th percentile of insert size.
+ More efficient insert-size calc

Credits
-------

Ryan Layer's [svv](https://github.com/ryanlayer/svv) was the inspiration for the plotter.

Obviously, [samtools](https://github.com/samtools/samtools) is the reference pileup implementation.

