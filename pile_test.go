package bigly_test

import (
	"testing"

	"github.com/biogo/hts/sam"
	"github.com/brentp/bigly"

	. "gopkg.in/check.v1"
)

var records = []*sam.Record{
	{Name: "r001/1", Pos: 6, Cigar: []sam.CigarOp{
		sam.NewCigarOp(sam.CigarMatch, 8),
		sam.NewCigarOp(sam.CigarInsertion, 2),
		sam.NewCigarOp(sam.CigarMatch, 4),
		sam.NewCigarOp(sam.CigarDeletion, 1),
		sam.NewCigarOp(sam.CigarMatch, 3),
	},
		Flags: sam.Paired | sam.ProperPair | sam.MateReverse | sam.Read1,
		Seq:   sam.NewSeq([]byte("TTAGATAAAGGATACTG")),
		Qual:  []uint8{0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff},
	},
}

type PileTest struct{}

var _ = Suite(&PileTest{})

func Test(t *testing.T) { TestingT(t) }

func (t *PileTest) TestAt(c *C) {
	r := bigly.Align{Record: records[0], CursorPos: 0, CursorCigar: 0}
	cig := r.At(0)
	c.Assert(cig, IsNil)

	cig = r.At(8)
	c.Assert(cig.At.Type(), Equals, sam.CigarMatch)

	for off := 13; off < 18; off++ {
		cig = r.At(off)
		c.Assert(cig.At.Type(), Equals, sam.CigarMatch)
	}
	cig = r.At(18)
	c.Assert(cig.At.Type(), Equals, sam.CigarDeletion)
	c.Assert(cig.Left.Type(), Equals, sam.CigarMatch)
	cig = r.At(19)
	c.Assert(cig.At.Type(), Equals, sam.CigarMatch)

}
func (t *PileTest) TestLeft(c *C) {
	r := bigly.Align{Record: records[0], CursorPos: 0, CursorCigar: 0}
	cig := r.At(18)
	c.Assert(cig.Left.Type(), Equals, sam.CigarMatch)
	cig = r.At(19)
	c.Assert(cig.Left.Type(), Equals, sam.CigarDeletion)
	cig = r.At(20)
	c.Assert(cig.Left.Type(), Equals, sam.CigarMatch)
}

func (t *PileTest) TestRight(c *C) {
	r := bigly.Align{Record: records[0], CursorPos: 0, CursorCigar: 0}
	cig := r.At(17)
	c.Assert(cig.Right.Type(), Equals, sam.CigarDeletion)
	cig = r.At(18)
	c.Assert(cig.Right.Type(), Equals, sam.CigarMatch)
	c.Assert(cig.Tail, Equals, false)

	cig = r.At(21)
	c.Assert(cig.Tail, Equals, true)
}

func (t *PileTest) TestHead(c *C) {
	r := bigly.Align{Record: records[0], CursorPos: 0, CursorCigar: 0}
	cig := r.At(r.Pos)
	c.Assert(cig.Head, Equals, true)

	cig = r.At(r.Pos + 1)
	c.Assert(cig.Head, Equals, false)

}

// NOTE: this isn't really a test. Just reminds me how to get bases from an Aln.
func (t *PileTest) TestExBases(c *C) {
	r := bigly.Align{Record: records[0], CursorCigar: 0, CursorPos: 0}
	r.Sequence = r.Record.Seq.Expand()

	var seqs []byte
	for off := 0; off < len(r.Sequence); off++ {
		seqs = append(seqs, r.Sequence[off])
	}
	c.Assert(string(seqs), Equals, "TTAGATAAAGGATACTG")

}

func (t *PileTest) TestBases(c *C) {
	r := bigly.Align{Record: records[0]}
	r.Sequence = r.Record.Seq.Expand()
	var cig *bigly.CigarSummary

	cig = r.At(r.Pos)
	c.Assert(string(cig.Base), Equals, "T")

	cig = r.At(r.Pos + 1)
	c.Assert(string(cig.Base), Equals, "T")

	cig = r.At(r.Pos + 2)
	c.Assert(string(cig.Base), Equals, "A")

	cig = r.At(r.Pos + 7)
	c.Assert(string(cig.Base), Equals, "A")
	c.Assert(cig.Right.Type(), Equals, sam.CigarInsertion)
	c.Assert(cig.Insertion, DeepEquals, []byte("AG"))

	// -- insertion between 7 and 8

	cig = r.At(r.Pos + 8)
	c.Assert(string(cig.Base), Equals, "G")
	c.Assert(cig.Left.Type(), Equals, sam.CigarInsertion)

	cig = r.At(r.Pos + 9)
	c.Assert(string(cig.Base), Equals, "A")

	cig = r.At(r.Pos + 10)
	c.Assert(string(cig.Base), Equals, "T")

	cig = r.At(r.Pos + 11)
	c.Assert(string(cig.Base), Equals, "A")

	// TODO: fix this.handle deletion
	cig = r.At(r.Pos + 12)
	c.Assert(string(cig.Base), Equals, "*")

	cig = r.At(r.Pos + 15)
	c.Assert(string(cig.Base), Equals, "G")

}

func (t *PileTest) TestIndelBases(c *C) {
	r := bigly.Align{Record: records[0]}
	r.Sequence = r.Record.Seq.Expand()
	cig := r.At(18)
	c.Assert(cig.At.Type(), Equals, sam.CigarDeletion)
	c.Assert(cig.Base, Equals, byte('*'))
}
