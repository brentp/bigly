package bigly

import (
	"bytes"
	"fmt"
	"log"
	"strconv"
	"strings"

	"github.com/biogo/hts/sam"
)

// SkipBase holds the byte used for a Skip Operation.
var SkipBase byte = '.'

func min(a, b int) int {
	if a > b {
		return b
	}
	return a
}

func max(a, b int) int {
	if a < b {
		return b
	}
	return a
}

// Options holds information about what is reported in the Pile
type Options struct {
	MinBaseQuality    uint8  `arg:"-q,help:base quality threshold"`
	MinMappingQuality uint8  `arg:"-Q,help:mapping quality threshold"`
	ExcludeFlag       uint16 `arg:"-F"`
	IncludeFlag       uint16 `arg:"-f"`
	MinClipLength     int    `arg:"-c,help:only count H/S clips of at least this length"`
	IncludeBases      bool   `arg:"-b,help:output each base and base quality score"`
	SplitterVerbosity int    `arg:"-s,help:0-only count; 1:count and single most frequent; 2:all SAs; 3:dont shorten positions"`
	ConcordantCutoff  int    `arg:"-o,help:distance beyond which mates are called discordant"`
}

// Pile holds the information about a single base.
type Pile struct {
	Chrom                 string
	Pos                   int
	Depth                 int    // count of reads passing filters.
	RefBase               byte   // Reference base a this position.
	MisMatches            uint32 // number of mismatches .
	ProperPairs           int    // count of reads with paired flag
	SoftStarts            uint32 // counts of base preceding an 'S' cigar op
	SoftEnds              uint32 // ...  following ...
	HardStarts            uint32 // counts of base preceding an 'H' cigar op
	HardEnds              uint32
	InsertionStarts       uint32 // counts of base preceding an 'I' cigar op
	InsertionEnds         uint32
	Deletions             uint32  // counts of deletions 'D' at this base
	Heads                 uint32  // counts of starts of reads at this base
	Tails                 uint32  // counts of ends of reads at this base
	Splitters             uint32  // count of non-secondary reads with SA tags.
	Splitters1            uint32  // count of non-secondary reads with exactly 1 SA tag.
	Bases                 []byte  // All bases from reads covering this position
	Quals                 []uint8 // All quals from reads covering this position
	InsertSizeLPs         []int32 // All insert size for left-most of pair
	InsertSizeRMs         []int32 // ...                 right-most of pair
	OrientationPlusPlus   uint32  // Paired reads mapped in +/+ orientation
	OrientationMinusMinus uint32  // Paired reads mapped in -/- orientation
	OrientationMinusPlus  uint32  // Paired reads mapped in -/+ orientation
	// count of reads where splitters were in different orientation than read.
	// only set if SplitterVerbosity > 1
	OrientationSplitter uint32
	/*
		TODO: figure out how to do this because these are excluded by ExcludeFlag.
			   could try to calculate and it will ony be non-zero if these are not excluded.C
		Duplicates             uint32  // reads counted as duplicates
		Supplementary          uint32  // reads counted as supplementary
	*/
	Discordant             uint32  // Number of reads with insert size > ConcordantCutoff
	DiscordantChrom        uint32  // Number of reads mapping on different chroms
	DiscordantChromEntropy float32 // high value means all discordants came from same chrom.
	GC65                   uint32
	GC257                  uint32
	Duplicity65            float32 // measure of lack of sequence entropy.
	Duplicity257           float32 // measure of lack of sequence entropy.
	SplitterPositions      []Position
}

// from biogo/hts
func formatQual(q []byte) []byte {
	for _, v := range q {
		if v != 0xff {
			a := make([]byte, len(q))
			for i, p := range q {
				a[i] = p + 33
			}
			return a
		}
	}
	return []byte{}
}

func formatInts(a []int) string {
	s := ""
	for _, v := range a {
		s += strconv.Itoa(v) + ","
	}
	if len(s) > 0 {
		return s[:len(s)-2]
	}
	return s
}

func mean(arr []int32) int {
	if len(arr) == 0 {
		return 0
	}
	var s float64
	for _, a := range arr {
		s += float64(a)
	}
	return int(0.5 + s/float64(len(arr)))
}

// TabString prints a tab-delimited version of the Pile
func (p Pile) TabString(o Options) string {
	spl := ""
	if o.SplitterVerbosity > 0 && len(p.SplitterPositions) > 0 {
		if o.SplitterVerbosity == 1 {
			m, c := pMode(p.SplitterPositions)
			spl = fmt.Sprintf("%d/%d/%d", m, c, len(p.SplitterPositions))
		} else {
			aspl := make([]string, len(p.SplitterPositions))
			for i, s := range p.SplitterPositions {
				aspl[i] = s.String()
			}
			spl = strings.Join(aspl, ",")
		}
	}
	return fmt.Sprintf("%s\t%d\t%d\t%c\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.2f\t%d\t%d\t%.2f\t%.2f\t%s",
		p.Chrom,
		p.Pos+1, p.Depth, p.RefBase, p.MisMatches,
		p.ProperPairs, p.SoftStarts, p.SoftEnds,
		p.HardStarts, p.HardEnds, p.InsertionStarts, p.InsertionEnds, p.Deletions,
		//p.Heads, p.Tails,
		p.Splitters,
		p.Splitters1,
		// string(p.RefBase), string(p.Bases), string(formatQual(p.Quals)),
		mean(p.InsertSizeLPs),
		mean(p.InsertSizeRMs),
		p.OrientationPlusPlus+p.OrientationMinusPlus+p.OrientationMinusMinus+p.OrientationSplitter,
		p.Discordant,
		p.DiscordantChrom,
		p.DiscordantChromEntropy,
		p.GC65,
		p.GC257,
		p.Duplicity65,
		p.Duplicity257,
		spl,
	)
}

func abs(a int) int {
	if a < 0 {
		return -a
	}
	return a
}

// Update the Pile with info from the Alignment if it meets the requirements in Options
func (p *Pile) Update(o Options, alns []*Align) {
	// the Total values hold the sum of 1 / val so we can calc harmonic Mean
	// with less susceptiblity to outliers.
	var discMates []int
	for _, a := range alns {
		if a.MapQ < o.MinMappingQuality {
			continue
		}
		s := a.At(p.Pos)
		if s == nil || s.Qual < o.MinBaseQuality {
			continue
		}

		if a.Flags&sam.Paired == sam.Paired {
			if a.Flags&sam.ProperPair == sam.ProperPair {
				p.ProperPairs++
			}
			if abs(a.Start()-a.MatePos) > o.ConcordantCutoff {
				p.Discordant++
			}
			if a.MateRef.ID() != a.Ref.ID() {
				p.DiscordantChrom++
				discMates = append(discMates, a.MateRef.ID())
			} else {
				// same chromosome.
				if a.Start() < a.MatePos && a.Flags&sam.Reverse != sam.Reverse {
					p.InsertSizeLPs = append(p.InsertSizeLPs, int32(a.MatePos-a.Start()))
				} else if a.Start() > a.MatePos && a.Flags&sam.Reverse == sam.Reverse {
					p.InsertSizeRMs = append(p.InsertSizeRMs, int32(float64(a.Start()-a.MatePos)))
				} else {
					if a.Flags&sam.Reverse == sam.Reverse && a.Flags&sam.MateReverse == sam.MateReverse {
						p.OrientationMinusMinus++
					} else if 0 == a.Flags&sam.Reverse && 0 == a.Flags&sam.MateReverse {
						p.OrientationPlusPlus++
					} else {
						p.OrientationMinusPlus++
					}
				}
			}

		}
		if a.Flags&sam.Secondary == 0 {
			if tags, ok := a.Record.Tag([]byte{'S', 'A'}); ok {
				p.Splitters++
				if bytes.Count([]byte(tags), []byte{';'}) <= 1 {
					p.Splitters1++
				}

				p.updateSplitters(o, tags, a.Strand() != -1)
			}
		}

		p.Depth++

		switch s.Right.Type() {
		case sam.CigarMatch:
			break
		case sam.CigarInsertion:
			p.InsertionStarts++
		case sam.CigarSoftClipped:
			if s.Right.Len() >= o.MinClipLength {
				p.SoftStarts++
			}
		case sam.CigarHardClipped:
			if s.Right.Len() >= o.MinClipLength {
				p.HardStarts++
			}
		}

		switch s.Left.Type() {
		case sam.CigarMatch:
			break
		case sam.CigarInsertion:
			p.InsertionEnds++
		case sam.CigarSoftClipped:
			if s.Left.Len() >= o.MinClipLength {
				p.SoftEnds++
			}
		case sam.CigarHardClipped:
			if s.Left.Len() >= o.MinClipLength {
				p.HardEnds++
			}
		}

		if s.Head {
			p.Heads++
		} else if s.Tail {
			p.Tails++
		}

		if s.At.Type() == sam.CigarDeletion {
			p.Deletions++
		}
		if o.IncludeBases {
			p.Bases = append(p.Bases, s.Base)
			p.Quals = append(p.Quals, s.Qual)
		}
		if s.Base != p.RefBase {
			p.MisMatches++
		}
	}

	if p.DiscordantChrom > 1 {
		p.DiscordantChromEntropy = float32(entropy(discMates))
	}
	// don't set this if we don't know the reference base.
	if p.RefBase == 'N' {
		p.MisMatches = 0
	}
}

// track the actual positions of the splitters and check the orientation.
func (p *Pile) updateSplitters(o Options, tags []byte, readStrand bool) {
	if o.SplitterVerbosity == 0 {
		return
	}
	sas := ParseSAs(tags)
	// if there is an orientation change, we only want to count it once.
	var orientationChange bool
	for _, sa := range sas {
		if sa.MapQ >= o.MinMappingQuality {
			p.SplitterPositions = append(p.SplitterPositions, Position{Chrom: string(sa.Chrom), Start: sa.Pos, End: sa.End(), Strand: sa.Strand})
			if readStrand != sa.Strand {
				orientationChange = true
			}
		}
	}
	if orientationChange {
		p.OrientationSplitter++
	}
}

// Align is a sam.Record with a cursor to track position in the read and reference.
type Align struct {
	*sam.Record
	// track the into the Cigar Ops.
	CursorCigar int
	// track basepairs into the Reference
	CursorPos int
	// track basepairs into the read
	CursorRead int

	// Sequence holds the expanded sequence.
	Sequence []byte
	// track that we're alwasy moving forward.
	lastPos int
}

// CigarSummary is a summary of the context of a position for 1 alignment
type CigarSummary struct {
	At        sam.CigarOp
	Left      sam.CigarOp
	Right     sam.CigarOp
	Qual      uint8
	Head      bool
	Tail      bool
	Base      byte
	Insertion []byte
}

// At returns the CigarOp for a particular genomic position of the given read.
func (a *Align) At(pos0 int) *CigarSummary {
	if pos0+1 <= a.lastPos {
		log.Fatal("can't use align on the same position or lower position's that previous calls")
	}
	a.lastPos = pos0 + 1
	pos := a.Pos + a.CursorPos
	if pos0 < pos || len(a.Cigar) == 0 {
		return nil
	}

	if a.Sequence == nil {
		a.Sequence = a.Seq.Expand()
	}

	res := &CigarSummary{Left: a.Cigar[max(a.CursorCigar-1, 0)]}
	res.Head = pos0 == pos && a.CursorPos == 0

	for _, co := range a.Cigar[a.CursorCigar:] {
		t := co.Type()
		con := t.Consumes()
		l := co.Len()
		lr := l * con.Reference
		lq := l * con.Query

		// we can update the Left if the current position is still less
		// then the requested position.
		if pos0 > pos {
			res.Left = co
		}
		if pos+lr > pos0 {
			res.At = co
			right := co
			res.Right = co
			// add distance into read. Subtract how far before la we are.
			readi := a.CursorRead + (pos0 - pos)
			if lq > 0 {
				// if we consume bases from the read, then grab them.
				res.Base = a.Sequence[readi]
				res.Qual = a.Qual[readi]
			} else {
				res.Base = '*'
				if co.Type() == sam.CigarSkipped {
					res.Base = SkipBase
				}
			}
			if pos+lr-1 == pos0 {
				idx := a.CursorCigar + 1
				if idx == len(a.Cigar) {
					res.Tail = true
				} else {
					res.Right = a.Cigar[idx]
					right = a.Cigar[idx]
				}
			}
			if res.Right.Type() == sam.CigarInsertion {
				res.Insertion = a.Sequence[readi+1 : readi+1+right.Len()]
			}

			return res
		}
		res.Left = co
		a.CursorCigar++
		a.CursorPos += lr
		a.CursorRead += lq
		pos += lr
	}
	return nil
}
