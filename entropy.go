package bigly

import (
	"math"
	"sort"
)

// i should be the refids of the chroms and nstates should be the number of chroms in the genomes.
// the value returned will be 0 if all chroms are the same and 1 if all chroms are different.
func entropy(i []int) float64 {
	// since this is for chroms from mates, we know they should fall in fixed range so we can use
	// a slice instead of a map to track frequncies:
	if len(i) < 2 {
		return 0
	}
	if len(i) == 2 {
		if i[0] == i[1] {
			return 0.0
		}
		return 0.6309297535714575
	} else if len(i) == 3 && i[0] == i[1] && i[0] == i[2] {
		return 0.0
	}
	counts := make([]int, 92)
	for _, v := range i {
		counts[v%92]++
	}
	var s float64
	k, n := 0, float64(len(i))

	for _, c := range counts {
		if c != 0 {
			k++
			s += float64(c) / n * math.Log(float64(c)/n)
		}
	}
	return -s / math.Log(float64(k+1))
}

// mode returns the most frequent value and the count of times it was seen.
func mode(a []int) (int, int) {
	if len(a) < 1 {
		return 0, 0
	}
	sort.Ints(a)
	m := a[0]
	imode := m
	most := 1
	current := 0

	for _, v := range a {
		if v == m {
			current++
		} else {
			if current > most {
				most = current
				imode = m
			}
			current = 1
			m = v
		}
	}
	if current > most {
		imode = m
		most = current
	}
	return imode, most
}
