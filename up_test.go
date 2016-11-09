package bigly_test

import (
	"github.com/biogo/hts/sam"
	"github.com/brentp/bigly"
	. "gopkg.in/check.v1"
)

// Setup structs for tests copied from hts/sam/sam_tests.go
func mustAux(a sam.Aux, err error) sam.Aux {
	if err != nil {
		panic(err)
	}
	return a
}

var precords = []*sam.Record{
	{
		Name: "r001",
		Pos:  6,
		MapQ: 30,
		Cigar: sam.Cigar{
			sam.NewCigarOp(sam.CigarMatch, 8),
			sam.NewCigarOp(sam.CigarInsertion, 2),
			sam.NewCigarOp(sam.CigarMatch, 4),
			sam.NewCigarOp(sam.CigarDeletion, 1),
			sam.NewCigarOp(sam.CigarMatch, 3),
		},
		Flags:   sam.Paired | sam.ProperPair | sam.MateReverse | sam.Read1,
		MatePos: 36,
		TempLen: 39,
		Seq:     sam.NewSeq([]byte("TTAGATAAAGGATACTG")),
		Qual:    []uint8{0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff},
	},
	{
		Name: "r002",
		Pos:  8,
		MapQ: 30,
		Cigar: sam.Cigar{
			sam.NewCigarOp(sam.CigarSoftClipped, 3),
			sam.NewCigarOp(sam.CigarMatch, 6),
			sam.NewCigarOp(sam.CigarPadded, 1),
			sam.NewCigarOp(sam.CigarInsertion, 1),
			sam.NewCigarOp(sam.CigarMatch, 4),
		},
		MatePos: -1,
		TempLen: 0,
		Seq:     sam.NewSeq([]byte("AAAAGATAAGGATA")),
		Qual:    []uint8{0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff},
	},
	{
		Name: "r003",
		Pos:  8,
		MapQ: 30,
		Cigar: sam.Cigar{
			sam.NewCigarOp(sam.CigarSoftClipped, 5),
			sam.NewCigarOp(sam.CigarMatch, 6),
		},
		MatePos: -1,
		TempLen: 0,
		Seq:     sam.NewSeq([]byte("GCCTAAGCTAA")),
		Qual:    []uint8{0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff},
		AuxFields: []sam.Aux{
			mustAux(sam.NewAux(sam.NewTag("SA"), "ref,29,-,6H5M,17,0;")),
		},
	},
	{
		Name: "r004",
		Pos:  15,
		MapQ: 30,
		Cigar: sam.Cigar{
			sam.NewCigarOp(sam.CigarMatch, 6),
			sam.NewCigarOp(sam.CigarSkipped, 14),
			sam.NewCigarOp(sam.CigarMatch, 5),
		},
		MatePos: -1,
		TempLen: 0,
		Seq:     sam.NewSeq([]byte("ATAGCTTCAGC")),
		Qual:    []uint8{0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff},
	},
	{
		Name: "r003",
		Pos:  28,
		MapQ: 17,
		Cigar: sam.Cigar{
			sam.NewCigarOp(sam.CigarHardClipped, 6),
			sam.NewCigarOp(sam.CigarMatch, 5),
		},
		Flags:   sam.Reverse | sam.Supplementary,
		MatePos: -1,
		TempLen: 0,
		Seq:     sam.NewSeq([]byte("TAGGC")),
		Qual:    []uint8{0xff, 0xff, 0xff, 0xff, 0xff},
		AuxFields: []sam.Aux{
			mustAux(sam.NewAux(sam.NewTag("SA"), "ref,9,+,5S6M,30,1;")),
		},
	},
	{
		Name: "r001",
		Pos:  36,
		MapQ: 30,
		Cigar: sam.Cigar{
			sam.NewCigarOp(sam.CigarMatch, 9),
		},
		Flags:   sam.Paired | sam.ProperPair | sam.Reverse | sam.Read2,
		MatePos: 6,
		TempLen: -39,
		Seq:     sam.NewSeq([]byte("CAGCGGCAT")),
		Qual:    []uint8{0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff},
		AuxFields: []sam.Aux{
			mustAux(sam.NewAux(sam.NewTag("NM"), uint(1))),
		},
	},
}

type UpTest struct {
	alns []*bigly.Align
}

var _ = Suite(&UpTest{})

func (t *UpTest) SetUpTest(c *C) {
	t.alns = make([]*bigly.Align, 0)
	for _, r := range precords {
		t.alns = append(t.alns, &bigly.Align{Record: r})
	}
}

/*
@SQ SN:ref  LN:45
@CO --------------------------------------------------------
@CO Coor     12345678901234  5678901234567890123456789012345
@CO ref      AGCATGTTAGATAA**GATAGCTGTGCTAGTAGGCAGTCAGCGCCAT
@CO --------------------------------------------------------
@CO +r001/1        TTAGATAAAGGATA*CTG
@CO +r002         aaaAGATAA*GGATA
@CO +r003       gcctaAGCTAA
@CO +r004                     ATAGCT..............TCAGC
@CO -r003                            ttagctTAGGC
@CO -r001/2                                        CAGCGGCAT
@CO --------------------------------------------------------
r001    99  ref 7   30  8M2I4M1D3M  =   37  39  TTAGATAAAGGATACTG   *
r002    0   ref 9   30  3S6M1P1I4M  *   0   0   AAAAGATAAGGATA  *
r003    0   ref 9   30  5S6M    *   0   0   GCCTAAGCTAA *   SA:Z:ref,29,-,6H5M,17,0;
r004    0   ref 16  30  6M14N5M *   0   0   ATAGCTTCAGC *
r003    2064    ref 29  17  6H5M    *   0   0   TAGGC   *   SA:Z:ref,9,+,5S6M,30,1;
r001    147 ref 37  30  9M  =   7   -39 CAGCGGCAT   *   NM:i:1
*/

func (t *UpTest) TestCigarFirstMatch(c *C) {

	c1, _ := sam.ParseCigar([]byte("85S36M129S"))
	c.Assert(bigly.FirstMatch(c1), Equals, 85)

	c2, _ := sam.ParseCigar([]byte("128S36M86S"))
	c.Assert(bigly.FirstMatch(c2), Equals, 128)

	c3, _ := sam.ParseCigar([]byte("165S33M52S"))
	c.Assert(bigly.FirstMatch(c3), Equals, 165)

	c4, _ := sam.ParseCigar([]byte("183S67M"))
	c.Assert(bigly.FirstMatch(c4), Equals, 183)
}

func (t *UpTest) TestReadPieces(c *C) {
	r := bigly.Align{Record: precords[0]}
	c.Assert(bigly.ReadPieces(r.Cigar), DeepEquals, []int{0, 8, 10, 17})

	cig, _ := sam.ParseCigar([]byte("1S22M"))
	c.Assert(bigly.ReadPieces(cig), DeepEquals, []int{1, 23})

	cig, _ = sam.ParseCigar([]byte("1S22M1S"))
	c.Assert(bigly.ReadPieces(cig), DeepEquals, []int{1, 23})

}

func (t *UpTest) TestRefPieces(c *C) {
	r := bigly.Align{Record: precords[0]}
	c.Assert(bigly.RefPieces(r.Pos, r.Cigar), DeepEquals, []int{6, 18, 19, 22})

	r = bigly.Align{Record: precords[1]}
	c.Assert(bigly.RefPieces(r.Pos, r.Cigar), DeepEquals, []int{8, 18})

	r = bigly.Align{Record: precords[2]}
	c.Assert(bigly.RefPieces(r.Pos, r.Cigar), DeepEquals, []int{8, 14})

	r = bigly.Align{Record: precords[3]}
	c.Assert(bigly.RefPieces(r.Pos, r.Cigar), DeepEquals, []int{15, 21, 35, 40})

	r = bigly.Align{Record: precords[4]}
	c.Assert(bigly.RefPieces(r.Pos, r.Cigar), DeepEquals, []int{28, 33})

	r = bigly.Align{Record: precords[5]}
	c.Assert(bigly.RefPieces(r.Pos, r.Cigar), DeepEquals, []int{36, 45})

}

func (t *UpTest) TestSimple(c *C) {
	alns := t.alns

	opts := bigly.Options{IncludeBases: true}

	p := &bigly.Pile{Chrom: "ref", Pos: 8}

	p.Update(opts, []*bigly.Align{alns[0]})
	c.Assert(p.Depth, Equals, 1)
	p.Update(opts, []*bigly.Align{alns[1]})
	c.Assert(p.Depth, Equals, 2)

	p.Update(opts, []*bigly.Align{alns[3]})
	c.Assert(p.Depth, Equals, 2)

	p.Update(opts, []*bigly.Align{alns[2]})
	c.Assert(string(p.Bases), Equals, "AAA")

	c.Assert(p.SoftEnds, Equals, uint32(2))
}

func (t *UpTest) TestInsertionEnds(c *C) {

	opts := bigly.Options{IncludeBases: true}
	p := &bigly.Pile{Chrom: "ref", Pos: 14}
	p.Update(opts, t.alns)
	c.Assert(p.InsertionEnds, Equals, uint32(2))
	c.Assert(string(p.Bases), Equals, "GG")
}

func (t *UpTest) TestInsertionStarts(c *C) {
	opts := bigly.Options{IncludeBases: true}
	p := &bigly.Pile{Chrom: "ref", Pos: 13}
	p.Update(opts, t.alns)

	c.Assert(p.InsertionStarts, Equals, uint32(1))
	c.Assert(string(p.Bases), Equals, "AAA")
	c.Assert(p.Tails, Equals, uint32(1))
}

func (t *UpTest) TestMinMapQ(c *C) {
	opts := bigly.Options{MinMappingQuality: 100, IncludeBases: true}
	p := &bigly.Pile{Chrom: "ref", Pos: 13}
	p.Update(opts, t.alns)
	c.Assert(string(p.Bases), Equals, "")
	c.Assert(p.Depth, Equals, 0)
}

func (t *UpTest) TestDel(c *C) {
	opts := bigly.Options{IncludeBases: true}
	p := &bigly.Pile{Chrom: "ref", Pos: 18}
	p.Update(opts, t.alns)
	c.Assert(string(p.Bases), Equals, "*G")
	c.Assert(p.Deletions, Equals, uint32(1))
}

func (t *UpTest) TestHard(c *C) {
	opts := bigly.Options{IncludeBases: true}
	p := &bigly.Pile{Chrom: "ref", Pos: 28}
	p.Update(opts, t.alns)
	c.Assert(string(p.Bases), Equals, ".T")
	c.Assert(p.HardEnds, Equals, uint32(1))
}
