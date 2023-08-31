package pslq

import (
	"flag"
	"fmt"
	"math"
	"math/big"
	"math/rand"
	"testing"
)

var verbose = false

func init() {
	flag.BoolVar(&verbose, "verbose", false, "set for more output in the tests")
}

// Print a vector
func printBigIntVector(t *testing.T, name string, x []big.Int) {
	for i := range x {
		t.Logf("%s[%d] = %d\n", name, i, &x[i])
	}
}

// Returns acot(x) in result
func acot(prec uint, x int64, result *big.Float) {
	var term, power, _x, _kp, x2, oldresult big.Float
	_x.SetPrec(prec).SetInt64(x)
	power.SetPrec(prec).SetInt64(1)
	power.Quo(&power, &_x) // 1/x
	x2.Mul(&_x, &_x)
	result.SetPrec(prec).SetInt64(0)
	positive := true
	for k := int64(1); ; k += 2 {
		oldresult.Set(result)
		kp := k
		if !positive {
			kp = -k
		}
		positive = !positive
		_kp.SetPrec(prec).SetInt64(kp)
		term.Quo(&power, &_kp)
		result.Add(result, &term)
		if oldresult.Cmp(result) == 0 {
			break
		}
		power.Quo(&power, &x2)
	}
}

// Returns pi using Machin's formula
func pi(prec uint, result *big.Float) {
	var tmp, _4 big.Float
	_4.SetPrec(prec).SetInt64(4)
	acot(prec, 5, &tmp)
	tmp.SetPrec(prec).Mul(&tmp, &_4)
	acot(prec, 239, result)
	result.Sub(&tmp, result)
	result.SetPrec(prec).Mul(result, &_4)
}

// Evaluates a BBP term
//
// sum(k=0->inf)(1/base**k * (1/a*k + b))
func bbp(prec uint, base, a, b int64, result *big.Float) {
	var term, power, aFp, bFp, _1, k, _base, oldresult big.Float
	power.SetPrec(prec).SetInt64(1)
	result.SetPrec(prec).SetInt64(0)
	aFp.SetPrec(prec).SetInt64(a)
	bFp.SetPrec(prec).SetInt64(b)
	_1.SetPrec(prec).SetInt64(1)
	k.SetPrec(prec).SetInt64(0)
	_base.SetPrec(prec).SetInt64(base)
	for {
		oldresult.Set(result)
		term.Mul(&aFp, &k)
		term.Add(&term, &bFp)
		term.Quo(&_1, &term)
		term.Mul(&term, &power)
		result.Add(result, &term)
		if oldresult.Cmp(result) == 0 {
			break
		}
		power.Quo(&power, &_base)
		k.Add(&k, &_1)
	}
}

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

// check that the PSLQ actually found a relation
func checkResult(t *testing.T, in []big.Float, out []big.Int) {
	if len(in) != len(out) {
		t.Fatalf("lengths wrong of answers got %d expecting %d", len(in), len(out))
	}
	prec := in[0].Prec()
	var sum big.Float
	sum.SetPrec(prec)
	expectedPrecision := int(prec) - 1
	for i := range in {
		if out[i].Sign() == 0 {
			continue
		}
		var tmp big.Float
		tmp.SetPrec(prec)
		tmp.SetInt(&out[i])
		tmp.Mul(&tmp, &in[i])
		sum.Add(&sum, &tmp)
		expectedPrecision -= 1
	}
	if sum.Sign() != 0 {
		sumExponent := -sum.MantExp(nil)
		if sumExponent < expectedPrecision {
			t.Errorf("Expecting to sum to 0 got %g instead with precision %d bits but wanted %d bits", &sum, sumExponent, expectedPrecision)
		}
	}
}

func TestAlgorithm1(t *testing.T) {
	testPSLQ(t, 1, 64)
}

func TestAlgorithm2(t *testing.T) {
	// algorithm 2 needs more precision to come to an answer
	testPSLQ(t, 2, 256)
}

func testPSLQ(t *testing.T, algorithm int, minPrec uint) {
	new := func(prec uint) *Pslq {
		return New(prec).SetVerbose(verbose).SetAlgorithm(algorithm)
	}

	getPrec := func(prec uint) uint {
		if prec >= minPrec {
			return prec
		}
		return minPrec
	}

	t.Run("ErrorBadArguments", func(t *testing.T) {
		prec := minPrec
		in := make([]big.Float, 1)
		pslq := new(prec)
		_, err := pslq.Run(in)
		if err != ErrorBadArguments {
			t.Error("Expecting", ErrorBadArguments, "but got", err)
		}
	})

	t.Run("ErrorPrecisionTooLow", func(t *testing.T) {
		prec := uint(63)
		in := make([]big.Float, 2)
		in[0].SetPrec(prec).SetInt64(1)
		in[1].SetPrec(prec).SetInt64(-2)
		pslq := new(prec)
		_, err := pslq.Run(in)
		if err != ErrorPrecisionTooLow {
			t.Error("Expecting", ErrorPrecisionTooLow, "but got", err)
		}
	})

	t.Run("ErrorToleranceRoundsToZero", func(t *testing.T) {
		// FIXME Can't test this until can pass tolerance in
		// if err != ErrorToleranceRoundsToZero {
		// 	t.Error("Expecting", ErrorToleranceRoundsToZero, "but got", err)
		// }
	})

	t.Run("ErrorZeroArguments", func(t *testing.T) {
		prec := minPrec
		in := make([]big.Float, 2)
		in[0].SetPrec(prec).SetInt64(0)
		in[1].SetPrec(prec).SetInt64(-2)
		pslq := new(prec)
		_, err := pslq.Run(in)
		if err != ErrorZeroArguments {
			t.Error("Expecting", ErrorZeroArguments, "but got", err)
		}
	})

	t.Run("ErrorArgumentTooSmall", func(t *testing.T) {
		prec := minPrec
		in := make([]big.Float, 2)
		tol := math.Pow(2, -float64(3*prec/4))
		in[0].SetPrec(prec).SetFloat64(tol / 129)
		in[1].SetPrec(prec).SetInt64(-2)
		pslq := new(prec)
		_, err := pslq.Run(in)
		if err != ErrorArgumentTooSmall {
			t.Error("Expecting", ErrorArgumentTooSmall, "but got", err)
		}
		in[0].SetFloat64(tol / 127)
		_, err = pslq.Run(in)
		if err == ErrorArgumentTooSmall {
			t.Error("Not expecting", err)
		}
	})

	// assert pslq([3*pi+4*e/7, pi, e, log(2)]) == [7, -21, -4, 0]
	// assert pslq([4.9999999999999991, 1]) == [1, -5]
	// assert pslq([2,1]) == [1, -2]

	t.Run("PslqSimple", func(t *testing.T) {
		prec := minPrec
		//one := float64(1<<60)

		in := make([]big.Float, 2)
		in[0].SetPrec(prec).SetInt64(1)
		in[1].SetPrec(prec).SetInt64(-2)

		pslq := new(prec)
		out, err := pslq.Run(in)
		if err != nil {
			t.Error("Got error", err)
		}
		if verbose {
			printBigIntVector(t, "out", out)
		}
		compareResult(t, out, 2, 1)
		checkResult(t, in, out)
	})

	t.Run("Pslq2", func(t *testing.T) {
		if minPrec > 64 {
			t.Skip("Can't do float64 test when minPrec is above 64")
		}
		prec := minPrec

		inFloat := []float64{
			3*math.Pi + 4*math.E/7,
			math.Pi,
			math.E,
			math.Log(2),
		}
		in := make([]big.Float, len(inFloat))
		for i := range inFloat {
			in[i].SetPrec(prec).SetFloat64(inFloat[i])
		}
		pslq := new(prec)
		out, err := pslq.Run(in)
		if err != nil {
			t.Error("Got error", err)
		}
		if verbose {
			printBigIntVector(t, "out", out)
		}
		compareResult(t, out, 7, -21, -4, 0)
		// checkResult(t, in, out) - float64's above have incorrect precision
	})

	t.Run("Pslq3", func(t *testing.T) {
		if minPrec > 64 {
			t.Skip("Can't do float64 test when minPrec is above 64")
		}
		prec := minPrec

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
		in := make([]big.Float, len(inFloat))
		for i := range inFloat {
			in[i].SetPrec(prec).SetFloat64(inFloat[i])
		}
		pslq := new(prec).SetMaxSteps(1000)
		out, err := pslq.Run(in)
		if err != nil {
			t.Error("Got error", err)
		}
		if verbose {
			printBigIntVector(t, "out", out)
		}
		compareResult(t, out, 7, -21, -4, 0, 0, 0, 0, 0, 0, 0)
		// checkResult(t, in, out) - float64's above have incorrect precision
	})

	// Original BBP series
	t.Run("PslqOriginalBBP", func(t *testing.T) {
		prec := minPrec

		in := make([]big.Float, 8)
		for i := range in {
			if i == 0 {
				pi(prec, &in[i])
			} else {
				bbp(prec, 16, 8, int64(i), &in[i])
			}
			if verbose {
				t.Logf("in[%d] = %g\n", i, &in[i])
			}
		}
		pslq := new(prec).SetMaxSteps(1000)
		out, err := pslq.Run(in)
		if err != nil {
			t.Error("Got error", err)
		}
		if verbose {
			printBigIntVector(t, "out", out)
		}
		compareResult(t, out, 1, -4, 0, 0, 2, 1, 1, 0)
		checkResult(t, in, out)
	})

	// New BBP series
	t.Run("PslqBBP64", func(t *testing.T) {
		prec := getPrec(128)

		in := make([]big.Float, 12)
		for i := range in {
			if i == 0 {
				pi(prec, &in[i])
			} else {
				bbp(prec, -64, 12, int64(i), &in[i])
			}
			if verbose {
				t.Logf("in[%d] = %g\n", i, &in[i])
			}
		}
		pslq := new(prec).SetMaxSteps(1000)
		out, err := pslq.Run(in)
		if err != nil {
			t.Error("Got error", err)
		}
		if verbose {
			printBigIntVector(t, "out", out)
		}
		// This isn't a Pi relation but it is what we find!
		compareResult(t, out, 0, 0, 16, -16, -24, 0, 8, 4, 6, 6, 1, -1)
		checkResult(t, in, out)
	})

	// Attempt to find Bellard's BBP series
	t.Run("PslqBellard", func(t *testing.T) {
		prec := getPrec(1024)

		in := make([]big.Float, 20)
		for i := range in {
			if i == 0 {
				pi(prec, &in[i])
			} else {
				bbp(prec, -1024, 20, int64(i), &in[i])
			}
			if verbose {
				t.Logf("in[%d] = %g\n", i, &in[i])
			}
		}
		pslq := new(prec).SetMaxSteps(10000).SetMaxCoeff(big.NewInt(1000000))
		out, err := pslq.Run(in)
		if err != nil {
			t.Error("Got error", err)
		}
		if verbose {
			printBigIntVector(t, "out", out)
		}
		// This isn't Bellard's formula but it is a valid formula for pi
		compareResult(t, out, 192, -512, 0, -256, 0, -32, 0, 64, 0, -32, -40, -16, 0, 8, 0, -1, 0, -2, 0, -1)
		checkResult(t, in, out)
	})

	piFailTest := func(t *testing.T, prec uint, iterations int, logMaxCoeff int64, expectedError error) {
		if prec >= 2048 && testing.Short() {
			t.Skip("skipping test in short mode.")
		}
		_10 := big.NewInt(10)
		maxCoeff := big.NewInt(logMaxCoeff)
		maxCoeff.Exp(_10, maxCoeff, nil)

		in := make([]big.Float, 5)
		for i := range in {
			if i == 0 {
				pi(prec, &in[i])
			} else {
				bbp(prec, 10, 5, int64(i), &in[i])
			}
			if verbose {
				t.Logf("in[%d] = %g\n", i, &in[i])
			}
		}
		pslq := new(prec).SetMaxSteps(iterations).SetMaxCoeff(maxCoeff)
		out, err := pslq.Run(in)
		if err == nil {
			t.Errorf("Expecting error but got none")
		} else {
			if expectedError != nil && err != expectedError {
				t.Errorf("Wrong error: want: %v, got: %v", expectedError, err)
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

	t.Run("Fail", func(t *testing.T) {
		for _, test := range []struct {
			prec          uint
			iterations    int
			logMaxCoeff   int64
			expectedError error
		}{
			{64, 75, 0, ErrorNoRelationFound},
			{64, 75, 3, ErrorPrecisionExhausted},
			{128, 75, 3, ErrorIterationsExceeded},
			{128, 150, 3, ErrorNoRelationFound},
			{128, 200, 10, ErrorPrecisionExhausted},
			{256, 50, 10, ErrorIterationsExceeded},
			{256, 1e8, 3, ErrorNoRelationFound},
			{256, 1e8, 13, ErrorPrecisionExhausted},
			{512, 1e8, 13, ErrorNoRelationFound},
			{512, 1e8, 18, ErrorNoRelationFound},
			{512, 1e8, 25, ErrorPrecisionExhausted},
			{1024, 1e8, 25, ErrorNoRelationFound},
			{1024, 1e8, 46, ErrorNoRelationFound},
			{1024, 1e8, 51, ErrorPrecisionExhausted},
			{2048, 1e8, 49, ErrorNoRelationFound},
			{2048, 1e8, 98, ErrorNoRelationFound},
			{2048, 1e8, 115, ErrorPrecisionExhausted},
			{2048 + 64, 1e8, 99, ErrorNoRelationFound},
			{67, 721, 12, ErrorPrecisionExhausted},
		} {
			if test.prec < minPrec {
				continue
			}
			t.Run(fmt.Sprintf("prec=%d,iters=%d,logMaxCoeff=%d", test.prec, test.iterations, test.logMaxCoeff), func(t *testing.T) {
				piFailTest(t, test.prec, test.iterations, test.logMaxCoeff, test.expectedError)
			})
		}
	})

	// Run lots of tests to check we always get sensible errors and not an invalid solution
	t.Run("FailRandom", func(t *testing.T) {
		iterations := 100
		if testing.Short() {
			iterations = 10
		}
		for i := 0; i < iterations; i++ {
			precision := rand.Intn(512) + 64
			iterations := rand.Intn(1000) + 64
			logMaxCoeff := rand.Intn(28) + 2
			piFailTest(t, getPrec(uint(precision)), iterations, int64(logMaxCoeff), nil)
		}
	})

	t.Run("Pslq4b", func(t *testing.T) {
		t.Skip("Takes too long")
		prec := getPrec(1024 * 2)

		in := make([]big.Float, 50)
		for i := range in {
			if i == 0 {
				pi(prec, &in[i])
			} else {
				bbp(prec, 100, 50, int64(i), &in[i])
			}
			if verbose {
				t.Logf("in[%d] = %g\n", i, &in[i])
			}
		}
		pslq := new(prec).SetMaxCoeff(big.NewInt(1e18)).SetMaxSteps(1e6)
		out, err := pslq.Run(in)
		if err == nil || err.Error() != "could not find an integer relation" {
			t.Errorf("Wrong error %v", err)
		}
		if verbose {
			printBigIntVector(t, "out", out)
		}
		if out != nil {
			t.Errorf("Expecting nil out, got %v", out)
		}
	})

	t.Run("Pslq5", func(t *testing.T) {
		prec := minPrec

		in := make([]big.Float, 8)
		pi(prec, &in[0])
		acot(prec, 2, &in[1])
		acot(prec, 4, &in[2])
		acot(prec, 6, &in[3])
		acot(prec, 7, &in[4])
		acot(prec, 8, &in[5])
		acot(prec, 9, &in[6])
		acot(prec, 10, &in[7])
		pslq := new(prec).SetMaxSteps(1000)
		out, err := pslq.Run(in)
		if err != nil {
			t.Error("Got error", err)
		}
		if verbose {
			printBigIntVector(t, "out", out)
		}
		compareResult(t, out, 1, -8, 0, 0, 4, 0, 0, 0)
		checkResult(t, in, out)
	})

	t.Run("Pslq6", func(t *testing.T) {
		prec := minPrec

		in := make([]big.Float, 3)
		pi(prec, &in[0])
		acot(prec, 5, &in[1])
		acot(prec, 239, &in[2])
		pslq := new(prec).SetMaxSteps(1000)
		out, err := pslq.Run(in)
		if err != nil {
			t.Error("Got error", err)
		}
		if verbose {
			printBigIntVector(t, "out", out)
		}
		compareResult(t, out, 1, -16, 4)
		checkResult(t, in, out)
	})

	t.Run("Pslq7", func(t *testing.T) {
		prec := minPrec

		in := make([]big.Float, 5)
		pi(prec, &in[0])
		acot(prec, 49, &in[1])
		acot(prec, 57, &in[2])
		acot(prec, 239, &in[3])
		acot(prec, 110443, &in[4])
		pslq := new(prec).SetMaxSteps(1000)
		out, err := pslq.Run(in)
		if err != nil {
			t.Error("Got error", err)
		}
		if verbose {
			printBigIntVector(t, "out", out)
		}
		compareResult(t, out, 1, -48, -128, 20, -48)
		checkResult(t, in, out)
	})

	t.Run("Pslq8", func(t *testing.T) {
		prec := getPrec(256)

		in := make([]big.Float, 16)
		for i := range in {
			if i == 0 {
				pi(prec, &in[i])
			} else {
				bbp(prec, 16, 8, int64(i), &in[i])
			}
			if verbose {
				t.Logf("in[%d] = %g\n", i, &in[i])
			}
		}
		pslq := new(prec).SetMaxSteps(10000)
		out, err := pslq.Run(in)
		if err != nil {
			t.Error("Got error", err)
		}
		if verbose {
			printBigIntVector(t, "out", out)
		}
		compareResult(t, out, 1, -4, 0, 0, 2, 1, 1)
		checkResult(t, in, out)
	})
}

func piFailBench(b *testing.B, prec uint, iterations int, logMaxCoeff int64) {
	_10 := big.NewInt(10)
	maxCoeff := big.NewInt(logMaxCoeff)
	maxCoeff.Exp(_10, maxCoeff, nil)

	pslq := New(prec).SetVerbose(verbose).SetMaxCoeff(maxCoeff).SetMaxSteps(iterations)

	in := make([]big.Float, 5)
	for i := range in {
		if i == 0 {
			pi(prec, &in[i])
		} else {
			bbp(prec, 10, 5, int64(i), &in[i])
		}
	}

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		pslq.Run(in)
	}
}

func BenchmarkPslqPiFail128(b *testing.B)  { piFailBench(b, 128, 150, 3) }
func BenchmarkPslqPiFail256(b *testing.B)  { piFailBench(b, 256, 1e8, 10) }
func BenchmarkPslqPiFail512(b *testing.B)  { piFailBench(b, 512, 1e8, 22) }
func BenchmarkPslqPiFail1024(b *testing.B) { piFailBench(b, 1024, 1e8, 46) }
