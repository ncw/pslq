package pslq

import (
	"fp"
	"math"
	"math/big"
	"math/rand"
	"testing"
)

const verbose = false

func compareResult(t *testing.T, actual []big.Int, expected ...int64) {
	if len(actual) < len(expected) {
		t.Fatalf("lengths wrong of answers got %d expecting %d", len(actual), len(expected))
	}
	for i := range actual {
		var e big.Int
		if i >= len(expected) {
			e.SetInt64(0)
		} else {
			e.SetInt64(expected[i])
		}
		if actual[i].Cmp(&e) != 0 {
			t.Errorf("actual[%d]=%d != expected[%d]=%d", i, &actual[i], i, &e)
		}
	}
}

func TestErrorBadArguments(t *testing.T) {
	env := fp.NewEnvironment(64)
	in := make([]fp.FixedPoint, 1)
	_, err := Pslq(env, in, nil, 0, verbose)
	if err != ErrorBadArguments {
		t.Error("Expecting", ErrorBadArguments, "but got", err)
	}
}

func TestErrorPrecisionTooLow(t *testing.T) {
	env := fp.NewEnvironment(63)
	in := make([]fp.FixedPoint, 2)
	in[0].SetInt64(env, 1)
	in[1].SetInt64(env, -2)
	_, err := Pslq(env, in, nil, 0, verbose)
	if err != ErrorPrecisionTooLow {
		t.Error("Expecting", ErrorPrecisionTooLow, "but got", err)
	}
}

func TestErrorToleranceRoundsToZero(t *testing.T) {
	// FIXME Can't test this until can pass tolerance in
	// if err != ErrorToleranceRoundsToZero {
	// 	t.Error("Expecting", ErrorToleranceRoundsToZero, "but got", err)
	// }
}

func TestErrorZeroArguments(t *testing.T) {
	env := fp.NewEnvironment(64)
	in := make([]fp.FixedPoint, 2)
	in[0].SetInt64(env, 0)
	in[1].SetInt64(env, -2)
	_, err := Pslq(env, in, nil, 0, verbose)
	if err != ErrorZeroArguments {
		t.Error("Expecting", ErrorZeroArguments, "but got", err)
	}
}

func TestErrorArgumentTooSmall(t *testing.T) {
	env := fp.NewEnvironment(64)
	in := make([]fp.FixedPoint, 2)
	in[0].Init(env)
	in[1].SetInt64(env, -2)
	const minOkValue = 1 << (64/4 - 7)
	in[0].Int.SetInt64(minOkValue - 1)
	_, err := Pslq(env, in, nil, 0, verbose)
	if err != ErrorArgumentTooSmall {
		t.Error("Expecting", ErrorArgumentTooSmall, "but got", err)
	}
	in[0].Int.SetInt64(minOkValue)
	_, err = Pslq(env, in, nil, 0, verbose)
	if err == ErrorArgumentTooSmall {
		t.Error("Not expecting", err)
	}
}

// assert pslq([3*pi+4*e/7, pi, e, log(2)]) == [7, -21, -4, 0]
// assert pslq([4.9999999999999991, 1]) == [1, -5]
// assert pslq([2,1]) == [1, -2]

func TestPslqSimple(t *testing.T) {
	env := fp.NewEnvironment(64)
	//one := float64(1<<60)

	in := make([]fp.FixedPoint, 2)
	in[0].SetInt64(env, 1)
	in[1].SetInt64(env, -2)

	out, err := Pslq(env, in, nil, 0, verbose)
	if err != nil {
		t.Error("Got error", err)
	}
	if verbose {
		printBigIntVector(t, "out", out)
	}
	compareResult(t, out, 2, 1)
}

func TestPslq2(t *testing.T) {
	env := fp.NewEnvironment(64)

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
	out, err := Pslq(env, in, nil, 0, verbose)
	if err != nil {
		t.Error("Got error", err)
	}
	if verbose {
		printBigIntVector(t, "out", out)
	}
	compareResult(t, out, 7, -21, -4, 0)
}

func TestPslq3(t *testing.T) {
	env := fp.NewEnvironment(64)

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
	out, err := Pslq(env, in, nil, 1000, verbose)
	if err != nil {
		t.Error("Got error", err)
	}
	if verbose {
		printBigIntVector(t, "out", out)
	}
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
		if verbose {
			t.Logf("in[%d] = %d\n", i, &in[i])
		}
	}
	out, err := Pslq(env, in, nil, 1000, verbose)
	if err != nil {
		t.Error("Got error", err)
	}
	if verbose {
		printBigIntVector(t, "out", out)
	}
	compareResult(t, out, 1, -4, 0, 0, 2, 1, 1, 0)
}

// Print a vector
func printBigIntVector(t *testing.T, name string, x []big.Int) {
	for i := range x {
		t.Logf("%s[%d] = %d\n", name, i, &x[i])
	}
}

func piFailTest(t *testing.T, prec uint, iterations int, logMaxCoeff int64, expectedError error) {
	if prec >= 2048 && testing.Short() {
		t.Skip("skipping test in short mode.")
	}
	env := fp.NewEnvironment(prec)
	_10 := big.NewInt(10)
	maxCoeff := big.NewInt(logMaxCoeff)
	maxCoeff.Exp(_10, maxCoeff, nil)

	in := make([]fp.FixedPoint, 5)
	for i := range in {
		if i == 0 {
			pi(env, &in[i])
		} else {
			bbp(env, 10, 5, int64(i), &in[i])
		}
		if verbose {
			t.Logf("in[%d] = %d\n", i, &in[i])
		}
	}
	out, err := Pslq(env, in, maxCoeff, iterations, verbose)
	if err == nil {
		t.Errorf("Expecting error but got none")
	} else {
		if expectedError != nil && err != expectedError {
			t.Errorf("Wrong error: %v != %v", err, expectedError)
		}
	}
	if verbose {
		printBigIntVector(t, "out", out)
	}
	if out != nil {
		t.Errorf("Expecting nil out, got")
		t.Errorf("piFailText(prec=%d, iterations=%d, logMaxCoeff=%d) = %v\n", prec, iterations, logMaxCoeff, err)
		printBigIntVector(t, "out", out)
		t.Fatalf("finishing")
	}
}

func TestPslqPiFail01(t *testing.T) { piFailTest(t, 64, 75, 0, ErrorNoRelationFound) }
func TestPslqPiFail02(t *testing.T) { piFailTest(t, 64, 75, 3, ErrorPrecisionExhausted) }
func TestPslqPiFail03(t *testing.T) { piFailTest(t, 128, 75, 3, ErrorIterationsExceeded) }
func TestPslqPiFail04(t *testing.T) { piFailTest(t, 128, 150, 3, ErrorNoRelationFound) }
func TestPslqPiFail05(t *testing.T) { piFailTest(t, 128, 200, 10, ErrorPrecisionExhausted) }
func TestPslqPiFail06(t *testing.T) { piFailTest(t, 256, 200, 10, ErrorIterationsExceeded) }
func TestPslqPiFail07(t *testing.T) { piFailTest(t, 256, 1E8, 10, ErrorNoRelationFound) }
func TestPslqPiFail08(t *testing.T) { piFailTest(t, 256, 1E8, 13, ErrorPrecisionExhausted) }
func TestPslqPiFail09(t *testing.T) { piFailTest(t, 512, 1E8, 13, ErrorNoRelationFound) }
func TestPslqPiFail10(t *testing.T) { piFailTest(t, 512, 1E8, 22, ErrorNoRelationFound) }
func TestPslqPiFail11(t *testing.T) { piFailTest(t, 512, 1E8, 25, ErrorPrecisionExhausted) }
func TestPslqPiFail12(t *testing.T) { piFailTest(t, 1024, 1E8, 25, ErrorNoRelationFound) }
func TestPslqPiFail13(t *testing.T) { piFailTest(t, 1024, 1E8, 46, ErrorNoRelationFound) }
func TestPslqPiFail14(t *testing.T) { piFailTest(t, 1024, 1E8, 49, ErrorPrecisionExhausted) }
func TestPslqPiFail15(t *testing.T) { piFailTest(t, 2048, 1E8, 49, ErrorNoRelationFound) }
func TestPslqPiFail16(t *testing.T) { piFailTest(t, 2048, 1E8, 98, ErrorNoRelationFound) }
func TestPslqPiFail17(t *testing.T) { piFailTest(t, 2048, 1E8, 99, ErrorPrecisionExhausted) }
func TestPslqPiFail18(t *testing.T) { piFailTest(t, 2048+64, 1E8, 99, ErrorNoRelationFound) }

func TestPslqPiFail20(t *testing.T) { piFailTest(t, 67, 721, 12, ErrorPrecisionExhausted) }

// Run lots of tests to check we always get sensible errors and not an invalid solution
func TestPslqPiFailRandom(t *testing.T) {
	iterations := 100
	if testing.Short() {
		iterations = 10
	}
	for i := 0; i < iterations; i++ {
		precision := rand.Intn(512) + 64
		iterations := rand.Intn(1000) + 64
		logMaxCoeff := rand.Intn(28) + 2
		piFailTest(t, uint(precision), iterations, int64(logMaxCoeff), nil)
	}
}

func piFailBench(b *testing.B, prec uint, iterations int, logMaxCoeff int64) {
	env := fp.NewEnvironment(prec)
	_10 := big.NewInt(10)
	maxCoeff := big.NewInt(logMaxCoeff)
	maxCoeff.Exp(_10, maxCoeff, nil)

	in := make([]fp.FixedPoint, 5)
	for i := range in {
		if i == 0 {
			pi(env, &in[i])
		} else {
			bbp(env, 10, 5, int64(i), &in[i])
		}
	}

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		Pslq(env, in, maxCoeff, iterations, verbose)
	}
}

func BenchmarkPslqPiFail128(b *testing.B)  { piFailBench(b, 128, 150, 3) }
func BenchmarkPslqPiFail256(b *testing.B)  { piFailBench(b, 256, 1E8, 10) }
func BenchmarkPslqPiFail512(b *testing.B)  { piFailBench(b, 512, 1E8, 22) }
func BenchmarkPslqPiFail1024(b *testing.B) { piFailBench(b, 1024, 1E8, 46) }

// FIXME try bpp base 100 as well? and 16?

// out = [44 120 -359 -665 431 -138 248 -166 146 -22 -5 20 339 -563 -606 -89 391 201 351 -31 -5 588 235 -663 183 646 -130 -73 11 167 -31 -788 666 -645 580 -15 -145 -523 -519 532 -169 686 43 80 -387 -234 560 486 285 -318]

func XXXXTestPslq4b(t *testing.T) {
	env := fp.NewEnvironment(1024 * 2)

	in := make([]fp.FixedPoint, 50)
	for i := range in {
		if i == 0 {
			pi(env, &in[i])
		} else {
			bbp(env, 100, 50, int64(i), &in[i])
		}
		if verbose {
			t.Logf("in[%d] = %d\n", i, &in[i])
		}
	}
	out, err := Pslq(env, in, big.NewInt(1E18), 1E6, verbose)
	if err == nil || err.Error() != "could not find an integer relation" {
		t.Errorf("Wrong error %v", err)
	}
	if verbose {
		printBigIntVector(t, "out", out)
	}
	if out != nil {
		t.Errorf("Expecting nil out, got %v", out)
	}
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
	out, err := Pslq(env, in, nil, 1000, verbose)
	if err != nil {
		t.Error("Got error", err)
	}
	if verbose {
		printBigIntVector(t, "out", out)
	}
	compareResult(t, out, 1, -8, 0, 0, 4, 0, 0, 0)
}

func TestPslq6(t *testing.T) {
	env := fp.NewEnvironment(64)

	in := make([]fp.FixedPoint, 3)
	in[0].SetFloat64(env, math.Pi/4)
	acot(env, 5, &in[1])
	acot(env, 239, &in[2])
	out, err := Pslq(env, in, nil, 1000, verbose)
	if err != nil {
		t.Error("Got error", err)
	}
	if verbose {
		printBigIntVector(t, "out", out)
	}
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
	out, err := Pslq(env, in, nil, 1000, verbose)
	if err != nil {
		t.Error("Got error", err)
	}
	if verbose {
		printBigIntVector(t, "out", out)
	}
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
		if verbose {
			t.Logf("in[%d] = %d\n", i, &in[i])
		}
	}
	out, err := Pslq(env, in, nil, 10000, verbose)
	if err != nil {
		t.Error("Got error", err)
	}
	if verbose {
		printBigIntVector(t, "out", out)
	}
	compareResult(t, out, 1, -4, 0, 0, 2, 1, 1)
}
