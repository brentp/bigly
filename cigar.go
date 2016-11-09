package bigly

import (
	"log"

	"github.com/biogo/hts/sam"
)

// RefPieces returns an array of start, end pairs of places where the given Cigar covers the reference.
// E.g. 6, 8M2I4M1D3M will return {6, 18, 19, 22} meaning 6-18, 19-22.
func RefPieces(pos int, c sam.Cigar) []int {
	if len(c) == 1 && c[0].Type() == sam.CigarMatch {
		return []int{pos, pos + c[0].Len()}
	}
	p := make([]int, 0, 4)

	for _, co := range c {
		con := co.Type().Consumes()
		if con.Reference > 0 {
			if con.Query != 0 {
				if len(p) == 0 || pos != p[len(p)-1] {
					p = append(p, pos)
					p = append(p, pos+co.Len())
				} else {
					p[len(p)-1] = pos + co.Len()
				}
			}
			pos += co.Len()
		}
	}
	return p
}

// ReadPieces gives the offsets into the cigar that consome reference bases.
func ReadPieces(c sam.Cigar) []int {
	if len(c) == 1 && c[0].Type() == sam.CigarMatch {
		return []int{0, c[0].Len()}
	}
	p := make([]int, 0, 4)
	off := 0
	for _, co := range c {
		con := co.Type().Consumes()
		if con.Query != 0 && con.Reference != 0 {
			if len(p) == 0 || off != p[len(p)-1] {
				p = append(p, off)
				p = append(p, off+co.Len())
			} else {
				p[len(p)-1] = off + co.Len()
			}
		}
		if con.Query != 0 {
			off += co.Len()
		}
	}
	return p
}

// FirstMatch reports the first base in the read that matches the reference.
func FirstMatch(c sam.Cigar) int {
	start := 0
	if c == nil {
		log.Fatal("no cigar to parse")
	}
	for _, co := range c {
		if co.Type() == sam.CigarMatch {
			return start
		}
		con := co.Type().Consumes()
		if con.Query > 0 {
			start += co.Len()
		}

	}
	return start
}
