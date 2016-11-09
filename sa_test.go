package bigly_test

import (
	"github.com/brentp/bigly"
	. "gopkg.in/check.v1"
)

type SATest struct{}

var _ = Suite(&SATest{})

func (t *SATest) TestSA(c *C) {
	sa := bigly.ParseSA([]byte("7,70999871,+,117S83M50S,42,8"))
	c.Assert(sa.Chrom, DeepEquals, []byte("7"))
	c.Assert(sa.Pos, Equals, 70999870)
	c.Assert(sa.Strand, Equals, int8(1))
	c.Assert(sa.Cigar, DeepEquals, []byte("117S83M50S"))
	c.Assert(sa.MapQ, Equals, uint8(42))
	c.Assert(sa.NM, Equals, uint16(8))
	c.Assert(sa.End(), Equals, sa.Pos+83)
}
