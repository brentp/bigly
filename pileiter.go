package bigly

import (
	"io"
	"log"

	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"
	"github.com/brentp/bigly/bamat"
	"github.com/brentp/faidx"
)

// Iterator generates piles the order they were queried.
// Internally, it holds a cache of alignment objects that overlap
// the current position. It requests new ones and drops old ones
// as they start or end overlapping with the (updated) pile position.
type Iterator struct {
	// underlying iterator from which to pull *sam.Records.
	bit   *bam.Iterator
	bamat *bamat.BamAt
	chrom string
	err   error
	// pile gets created and set on each call to Next()
	pile *Pile
	opts Options
	// track current genomic position
	pos int
	// hold all the alignments we (might) need to look at given the current pos
	cache   []*Align
	hasMore bool
	end     int
	fai     *faidx.Faidx
	gcs     [2]*faidx.FaPos
}

// Position indicates where the pileup is performed
type Position struct {
	Chrom string
	Start int
	End   int
}

// Up is the user-facing function that performs the pileup.
func Up(bampath string, opts Options, pos Position, fai *faidx.Faidx) *Iterator {
	b, err := bamat.New(bampath)
	if err != nil {
		return &Iterator{err: err}
	}

	bit, err := b.Query(pos.Chrom, pos.Start, pos.End)
	if err != nil {
		b.Close()
		return &Iterator{err: err}
	}

	it := &Iterator{bit: bit, bamat: b, pos: pos.Start, chrom: pos.Chrom, opts: opts, hasMore: true, end: pos.End}
	it.cache = make([]*Align, 0, 32)

	// prime the cache and potentially advance it to the start of the first read.
	for bit.Next() {
		rec := it.bit.Record()
		if !passes(rec, opts) {
			continue
		}
		it.cache = append(it.cache, &Align{Record: rec})
		if it.cache[0].Start() > it.pos {
			it.pos = it.cache[0].Start()
		}
		break
	}
	if bit.Error() != nil {
		log.Fatal(err)
	}
	if fai != nil {
		it.fai = fai
		it.gcs = [2]*faidx.FaPos{
			&faidx.FaPos{Chrom: pos.Chrom}, // 65
			&faidx.FaPos{Chrom: pos.Chrom}, // 256
		}
	}

	return it
}

// Error returns any error encountered by the Iterator
func (it *Iterator) Error() error {
	if it.err == io.EOF {
		return nil
	}
	return it.err
}

func passes(r *sam.Record, o Options) bool {
	if uint16(r.Flags)&o.IncludeFlag != o.IncludeFlag {
		return false
	}
	if uint16(r.Flags)&o.ExcludeFlag != 0 {
		return false
	}
	return r.MapQ >= o.MinMappingQuality
}

// Next returns true as long as any remaning pileups are available.
func (it *Iterator) Next() bool {
	if it.err != nil || it.pos >= it.end {
		return false
	}
	// drop from the cache where the end is < the current position. this could leave
	// on some unneeded alignments when a longer alignment precedes a shorter one,
	// that will just cause pile.Update to do a bit more work but wont affect they
	// result.
	var i int
	// we have to copy the items in the cache that we need to the start to avoid
	// a memory leak.
	for i = 0; i < len(it.cache) && it.cache[i].End() < it.pos; i++ {
	}
	if i > 0 {
		copy(it.cache, it.cache[i:])
		for j := i; j > 0; j-- {
			// nil out the items so we don't leak.
			it.cache[len(it.cache)-j] = nil
		}
		// then truncate.
		it.cache = it.cache[:len(it.cache)-i]
	}

	// add to the cache as long until the start of the most recently added record
	// is greater than the current position.
	hasMore := true
	for len(it.cache) == 0 || it.cache[len(it.cache)-1].Start() <= it.pos {
		if it.bit.Next() {
			rec := it.bit.Record()
			if !passes(rec, it.opts) {
				continue
			}
			it.cache = append(it.cache, &Align{Record: rec})
		} else {
			hasMore = false
			it.err = it.bit.Error()
			break
		}
	}
	if len(it.cache) == 0 && !hasMore {
		return false
	}
	if it.err == nil {
		it.pile = &Pile{Chrom: it.chrom, Pos: it.pos, RefBase: 'N'}
		if it.fai != nil {
			it.faiUpdate()
		}
		it.pile.Update(it.opts, it.cache)
		it.pos++
		// skip missing regions.
		if it.pile.Depth == 0 && len(it.cache) > 0 && it.cache[0].Start() > it.pos {
			it.pos = it.cache[0].Start()
		}
		return true
	}
	it.pile = nil
	return false
}

// update the stuff that relies on a fasta.
func (it *Iterator) faiUpdate() {
	it.gcs[0].Start, it.gcs[0].End = it.pos-32, it.pos+32
	it.gcs[1].Start, it.gcs[1].End = it.pos-128, it.pos+128
	var err error
	it.pile.GC65, err = it.fai.Q(it.gcs[0])
	if err != nil {
		log.Fatal(err)
	}
	it.pile.GC257, err = it.fai.Q(it.gcs[1])
	if err != nil {
		log.Fatal(err)
	}
	it.pile.RefBase, err = it.fai.At(it.chrom, it.pos)
	it.pile.Duplicity65 = it.gcs[0].Duplicity()
	it.pile.Duplicity257 = it.gcs[1].Duplicity()
}

// Pile returns the next pile from the iterator.
func (it *Iterator) Pile() *Pile { return it.pile }

// Close the underlying bam iterator and bam file.
func (it *Iterator) Close() error {
	it.bamat.Close()
	it.cache = it.cache[:0]
	if it.bit != nil {
		return it.bit.Close()
	}
	return nil
}
