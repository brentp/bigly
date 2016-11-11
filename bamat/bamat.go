package bamat

import (
	"os"

	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"
)

// BamAt holds bam index and the Bam Reader.
// Because BamAt holds the underlying os.File open, it is not
// safe to query from multiple go routines.
type BamAt struct {
	*bam.Reader
	idx  *bam.Index
	fh   *os.File
	Refs map[string]*sam.Reference
}

// New returns a BamAt from the given path to the indexed bam
func New(path string) (*BamAt, error) {
	f, err := os.Open(path + ".bai")
	if err != nil && len(path) > 4 {
		f, err = os.Open(path[:len(path)-4] + ".bai")
	}
	if err != nil {
		return nil, err
	}
	defer f.Close()
	idx, err := bam.ReadIndex(f)
	if err != nil {
		return nil, err
	}
	bamat := &BamAt{idx: idx}

	bamat.fh, err = os.Open(path)
	if err != nil {
		return nil, err
	}

	br, err := bam.NewReader(bamat.fh, 2)
	if err != nil {
		bamat.fh.Close()
		return nil, err
	}
	bamat.Reader = br
	hdr := br.Header()
	bamat.Refs = make(map[string]*sam.Reference, 40)
	for _, r := range hdr.Refs() {
		bamat.Refs[r.Name()] = r
	}
	return bamat, nil
}

// Query the BamAt with 0-base half-open interval.
func (b *BamAt) Query(chrom string, start int, end int) (*bam.Iterator, error) {
	ref := b.Refs[chrom]
	chunks, err := b.idx.Chunks(ref, start, end)
	if err != nil {
		return nil, err
	}
	it, err := bam.NewIterator(b.Reader, chunks)
	if err != nil {
		return nil, err
	}
	return it, err
}

// Close closes the underlying file and the Bam.Reader
func (b *BamAt) Close() error {
	if b == nil {
		return nil
	}
	if b.Reader != nil {
		b.Reader.Close()
	}
	if b.fh != nil {
		return b.fh.Close()
	}
	return nil
}
