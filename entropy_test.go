package bigly

import . "gopkg.in/check.v1"

type EntropyTest struct{}

var _ = Suite(&EntropyTest{})

func (t *EntropyTest) TestEven(c *C) {
	e := entropy([]int{1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 12})
	c.Assert(e > 0.9, Equals, true)

	e = entropy([]int{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1})
	c.Assert(e < 0.01, Equals, true)

	/*
		e = entropy([]int{1, 1})

		e = entropy([]int{0, 1})
	*/

}
