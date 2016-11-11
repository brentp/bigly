package bigly

import (
	"bytes"
	"fmt"
	"os"
	"strconv"

	"github.com/biogo/hts/sam"
)

// SA holds the attributes in a SAM SA:Z tag.
type SA struct {
	Chrom  []byte
	Pos    int
	Strand int8
	Cigar  []byte
	MapQ   uint8
	NM     uint16
	end    int
	Parsed sam.Cigar // parsed Cigar
}

func RecordToSA(r *sam.Record) *SA {
	return &SA{Chrom: []byte(r.Ref.Name()), Pos: r.Start(), Strand: r.Strand(), Parsed: r.Cigar}
}

func AsSAs(r *sam.Record, cigs []byte) []*SA {
	sas := ParseSAs(cigs)
	sas = append(sas, RecordToSA(r))
	return sas
}

// End returns the end of the Cigar string
func (s *SA) End() int {
	if s.end != 0 {
		return s.end
	}
	var err error
	if s.Parsed == nil {
		s.Parsed, err = sam.ParseCigar(s.Cigar)
	}
	if err != nil {
		fmt.Fprintf(os.Stderr, "warning bad cigar in SA tag: %s\n", s.Cigar)
		return s.Pos
	}
	s.end = s.Pos
	for _, co := range s.Parsed {
		if t := co.Type(); t != sam.CigarBack {
			s.end += co.Len() * t.Consumes().Reference
		}
	}
	return s.end
}

func ParseSAs(s []byte) []*SA {
	if s[len(s)-1] == ';' {
		s = s[:len(s)-1]
	}
	if len(s) > 3 && s[0] == 'S' && s[1] == 'A' && s[2] == 'Z' {
		s = s[3:]
	}
	ss := bytes.Split(s, []byte{';'})

	sas := make([]*SA, len(ss), len(ss)+1)
	for i, sa := range ss {
		tmp := ParseSA(sa)
		sas[i] = &tmp
	}
	return sas
}

// ParseSA returns an SA struct from the bytes
func ParseSA(sa []byte) SA {
	// "7,70999871,+,117S83M50S,42,8"
	//parts := bytes.SplitN(sa, []byte{','}, 6)
	var s SA
	off := 0
	for i := 0; i < 6; i++ {
		next := off + bytes.Index(sa[off:], []byte{','})
		if i == 5 && next < off {
			next = len(sa)
		}
		switch i {
		case 0:
			s.Chrom = sa[:next]
		case 1, 4, 5:
			p, err := strconv.Atoi(string(sa[off:next]))
			if err != nil {
				panic(err)
			}
			if i == 1 {
				s.Pos = p - 1
			} else if i == 4 {
				s.MapQ = uint8(p)
			} else {
				s.NM = uint16(p)
			}
		case 3:
			s.Cigar = sa[off:next]
		case 2:
			if sa[off] == '-' {
				s.Strand = -1
			} else {
				s.Strand = 1
			}

		}
		off = 1 + next
	}
	return s
}
