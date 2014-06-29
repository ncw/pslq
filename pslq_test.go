package pslq

import (
	"fmt"
	"fp"
	"math"
	"testing"
)

func compareResult(t *testing.T, actual []int64, expected ...int64) {
	if len(actual) < len(expected) {
		t.Fatalf("lengths wrong of answers got %d expecting %d", len(actual), len(expected))
	}
	for i := range actual {
		var e int64
		if i >= len(expected) {
			e = 0
		} else {
			e = expected[i]

		}
		if actual[i] != e {
			t.Errorf("actual[%d]=%d != expected[%d]=%d", i, &actual[i], i, e)
		}
	}
}

// assert pslq([3*pi+4*e/7, pi, e, log(2)]) == [7, -21, -4, 0]
// assert pslq([4.9999999999999991, 1]) == [1, -5]
// assert pslq([2,1]) == [1, -2]

func TestPslqSimple(t *testing.T) {
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

func TestPslq2(t *testing.T) {
	env := fp.NewEnvironment(63)

	inFloat := []float64{
		3*math.Pi + 4*math.E/7,
		math.Pi,
		math.E,
		math.Log(2),
	}
	in := make([]fp.FixedPoint, len(inFloat))
	for i := range inFloat {
		in[i].SetFloat64(env, inFloat[i])
	}
	out, err := Pslq(env, in, 0, 0, true)
	if err != nil {
		t.Error("Got error", err)
	}
	fmt.Printf("out = %v\n", out)
	compareResult(t, out, 7, -21, -4, 0)
}

func TestPslq3(t *testing.T) {
	env := fp.NewEnvironment(63)

	inFloat := []float64{
		3*math.Pi + 4*math.E/7,
		math.Pi,
		math.E,
		math.Log(2),
		0.28917320090206799499,
		0.57591529756863646394,
		0.55698607277729539344,
		0.54073048514703925260,
		0.99835889431176827458,
		0.11551877481656358526,
	}
	in := make([]fp.FixedPoint, len(inFloat))
	for i := range inFloat {
		in[i].SetFloat64(env, inFloat[i])
	}
	out, err := Pslq(env, in, 0, 1000, true)
	if err != nil {
		t.Error("Got error", err)
	}
	fmt.Printf("out = %v\n", out)
	compareResult(t, out, 7, -21, -4, 0, 0, 0, 0, 0, 0, 0)
}

// Evaluates a BBP term
//
// sum(k=0->inf)(1/base**k * (1/a*k + b))
func bbp(env *fp.Environment, base, a, b int64, result *fp.FixedPoint) {
	var term, power, aFp, bFp, _1 fp.FixedPoint
	power.SetInt64(env, 1)
	result.SetInt64(env, 0)
	aFp.SetInt64(env, a)
	bFp.SetInt64(env, b)
	_1.SetInt64(env, 1)
	for k := int64(0); ; k++ {
		term.MulInt64(&aFp, k)
		term.Add(&term, &bFp)
		term.Div(&_1, &term)
		term.Mul(&term, &power)
		if term.Sign() == 0 {
			break
		}
		result.Add(result, &term)
		power.DivInt64(&power, base)
	}
}

func TestPslq4(t *testing.T) {
	env := fp.NewEnvironment(64)

	in := make([]fp.FixedPoint, 8)
	for i := range in {
		if i == 0 {
			in[i].SetFloat64(env, math.Pi)
		} else {
			bbp(env, 16, 8, int64(i), &in[i])
		}
		fmt.Printf("in[%d] = %d\n", i, &in[i])
	}
	out, err := Pslq(env, in, 0, 1000, true)
	if err != nil {
		t.Error("Got error", err)
	}
	fmt.Printf("out = %v\n", out)
	compareResult(t, out, 1, -4, 0, 0, 2, 1, 1, 0)
}

// Returns acot(x) in result
func acot(env *fp.Environment, x int64, result *fp.FixedPoint) {
	var term, power fp.FixedPoint
	power.SetInt64(env, 1)
	power.DivInt64(&power, x) // 1/x
	x2 := x * x
	result.SetInt64(env, 0)
	positive := true
	for k := int64(1); ; k += 2 {
		kp := k
		if !positive {
			kp = -k
		}
		positive = !positive
		term.DivInt64(&power, kp)
		if term.Sign() == 0 {
			break
		}
		result.Add(result, &term)
		power.DivInt64(&power, x2)
	}
}

// Returns pi using Machin's formula
func pi(env *fp.Environment, result *fp.FixedPoint) {
	var tmp fp.FixedPoint
	acot(env, 5, &tmp)
	tmp.Lsh(&tmp, 2)
	acot(env, 239, result)
	result.Sub(&tmp, result)
	result.Lsh(result, 2)
}

func TestPslq5(t *testing.T) {
	env := fp.NewEnvironment(64)

	in := make([]fp.FixedPoint, 8)
	in[0].SetFloat64(env, math.Pi)
	acot(env, 2, &in[1])
	acot(env, 4, &in[2])
	acot(env, 6, &in[3])
	acot(env, 7, &in[4])
	acot(env, 8, &in[5])
	acot(env, 9, &in[6])
	acot(env, 10, &in[7])
	out, err := Pslq(env, in, 0, 1000, true)
	if err != nil {
		t.Error("Got error", err)
	}
	fmt.Printf("out = %v\n", out)
	compareResult(t, out, 1, -8, 0, 0, 4, 0, 0, 0)
}

func TestPslq6(t *testing.T) {
	env := fp.NewEnvironment(64)

	in := make([]fp.FixedPoint, 3)
	in[0].SetFloat64(env, math.Pi/4)
	acot(env, 5, &in[1])
	acot(env, 239, &in[2])
	out, err := Pslq(env, in, 0, 1000, true)
	if err != nil {
		t.Error("Got error", err)
	}
	fmt.Printf("out = %v\n", out)
	compareResult(t, out, 1, -4, 1)
}

func TestPslq7(t *testing.T) {
	env := fp.NewEnvironment(64)

	in := make([]fp.FixedPoint, 5)
	in[0].SetFloat64(env, math.Pi/4)
	acot(env, 49, &in[1])
	acot(env, 57, &in[2])
	acot(env, 239, &in[3])
	acot(env, 110443, &in[4])
	out, err := Pslq(env, in, 0, 1000, true)
	if err != nil {
		t.Error("Got error", err)
	}
	fmt.Printf("out = %v\n", out)
	compareResult(t, out, 1, -12, -32, 5, -12)
}

func TestPslq8(t *testing.T) {
	env := fp.NewEnvironment(256)

	in := make([]fp.FixedPoint, 16)
	for i := range in {
		if i == 0 {
			pi(env, &in[i])
		} else {
			bbp(env, 16, 8, int64(i), &in[i])
		}
		fmt.Printf("in[%d] = %d\n", i, &in[i])
	}
	out, err := Pslq(env, in, 0, 10000, true)
	if err != nil {
		t.Error("Got error", err)
	}
	fmt.Printf("out = %v\n", out)
	compareResult(t, out, 1, -4, 0, 0, 2, 1, 1)
}
