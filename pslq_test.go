package pslq

import (
	"testing"

	fp "fixedpoint"
	"fmt"
)

func compareResult(t *testing.T, actual []int64, expected ...int64) {
	if len(actual) != len(expected) {
		t.Fatalf("lengths wrong of answers got %d expecting %d", len(actual), len(expected))
	}
	for i := range actual {
		if actual[i] != expected[i] {
			t.Errorf("actual[%d]=%d != expected[%d]=%d", i, actual[i], i, expected[i])
		}
	}
}

func TestPslq(t *testing.T) {
	// assert pslq([3*pi+4*e/7, pi, e, log(2)]) == [7, -21, -4, 0]
	// assert pslq([4.9999999999999991, 1]) == [1, -5]
	// assert pslq([2,1]) == [1, -2]

	env := fp.NewEnvironment(63)
	//one := float64(1<<60)

	in := make([]fp.FixedPoint, 2)
	in[0].SetInt64(env, 1)
	in[1].SetInt64(env, -2)

	out, err := Pslq(env, in, 0, 0, true)
	if err != nil {
		t.Error("Got error", err)
	}
	fmt.Printf("out = %v\n", out)
	compareResult(t, out, 2, 1)
}
