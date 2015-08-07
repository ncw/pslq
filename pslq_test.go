package pslq

import (
	"fmt"
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

func TestSqrt(t *testing.T) {
	tests := []struct {
		prec uint
		in   float64
	}{
		{16, 0},
		{16, 1},
		{16, 4},
		{16, 10000},
		{16, 2},
		{64, 2},
		{256, 2},
		{1024, 1.5},
	}
	for _, test := range tests {
		x := new(big.Float).SetPrec(test.prec)
		x.SetFloat64(test.in)
		var got, got2, diff big.Float
		Sqrt(x, &got)
		got2.SetPrec(test.prec).Mul(&got, &got)
		diff.Sub(&got2, x)
		if diff.MinPrec() > 1 {
			t.Errorf("sqrt(%f) prec %d wrong got %.20f square %.20f expecting %f diff %g minprec %d", test.in, test.prec, &got, &got2, x, &diff, diff.MinPrec())
		}
	}
}

func TestSqrt2HighPrecision(t *testing.T) {
	prec := uint(100/math.Log10(2)) + 1 // 100 decimal digits
	x := new(big.Float).SetPrec(prec).SetInt64(2)
	out := new(big.Float).SetPrec(prec)
	Sqrt(x, out)
	result := fmt.Sprintf("%.100f", out)
	expected := "1.4142135623730950488016887242096980785696718753769480731766797379907324784621070388503875343276415727"
	if result != expected {
		t.Errorf("100 digits of Sqrt(2) wrong\nExpected %q\nActual   %q", expected, result)
	}
}

func TestSqrt3HighPrecision(t *testing.T) {
	prec := uint(1000/math.Log10(2)) + 1 // 1000 decimal digits
	x := new(big.Float).SetPrec(prec).SetInt64(3)
	out := new(big.Float).SetPrec(prec)
	Sqrt(x, out)
	result := fmt.Sprintf("%.1000f", out)
	expected := "1.7320508075688772935274463415058723669428052538103806280558069794519330169088000370811461867572485756756261414154067030299699450949989524788116555120943736485280932319023055820679748201010846749232650153123432669033228866506722546689218379712270471316603678615880190499865373798593894676503475065760507566183481296061009476021871903250831458295239598329977898245082887144638329173472241639845878553976679580638183536661108431737808943783161020883055249016700235207111442886959909563657970871684980728994932964842830207864086039887386975375823173178313959929830078387028770539133695633121037072640192491067682311992883756411414220167427521023729942708310598984594759876642888977961478379583902288548529035760338528080643819723446610596897228728652641538226646984200211954841552784411812865345070351916500166892944154808460712771439997629268346295774383618951101271486387469765459824517885509753790138806649619119622229571105552429237231921977382625616314688420328537166829386496119170497388363954959381"
	if result != expected {
		t.Errorf("1000 digits of Sqrt(3) wrong\nExpected %q\nActual   %q", expected, result)
	}
}

func checkErr(t *testing.T, err interface{}, expectedMsg string) {
	if err != nil {
		msg, ok := err.(string)
		if !ok {
			t.Error("expecting string error")
			return
		}
		if msg != expectedMsg {
			t.Errorf("wrong error returned %q expecting %q", msg, expectedMsg)
		}
	} else {
		t.Error("expecting panic but didn't get one")
	}
}

func TestSqrtNegative(t *testing.T) {
	defer func() { checkErr(t, recover(), "Sqrt of negative number") }()
	x := new(big.Float).SetInt64(-2)
	var z big.Float
	Sqrt(x, &z)
}

func TestErrorBadArguments(t *testing.T) {
	prec := uint(64)
	in := make([]big.Float, 1)
	_, err := Pslq(prec, in, nil, 0, verbose)
	if err != ErrorBadArguments {
		t.Error("Expecting", ErrorBadArguments, "but got", err)
	}
}

func TestErrorPrecisionTooLow(t *testing.T) {
	prec := uint(63)
	in := make([]big.Float, 2)
	in[0].SetPrec(prec).SetInt64(1)
	in[1].SetPrec(prec).SetInt64(-2)
	_, err := Pslq(prec, in, nil, 0, verbose)
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
	prec := uint(64)
	in := make([]big.Float, 2)
	in[0].SetPrec(prec).SetInt64(0)
	in[1].SetPrec(prec).SetInt64(-2)
	_, err := Pslq(prec, in, nil, 0, verbose)
	if err != ErrorZeroArguments {
		t.Error("Expecting", ErrorZeroArguments, "but got", err)
	}
}

func TestErrorArgumentTooSmall(t *testing.T) {
	prec := uint(64)
	in := make([]big.Float, 2)
	tol := math.Pow(2, -float64(3*prec/4))
	in[0].SetPrec(prec).SetFloat64(tol / 129)
	in[1].SetPrec(prec).SetInt64(-2)
	_, err := Pslq(prec, in, nil, 0, verbose)
	if err != ErrorArgumentTooSmall {
		t.Error("Expecting", ErrorArgumentTooSmall, "but got", err)
	}
	in[0].SetFloat64(tol / 127)
	_, err = Pslq(prec, in, nil, 0, verbose)
	if err == ErrorArgumentTooSmall {
		t.Error("Not expecting", err)
	}
}

// assert pslq([3*pi+4*e/7, pi, e, log(2)]) == [7, -21, -4, 0]
// assert pslq([4.9999999999999991, 1]) == [1, -5]
// assert pslq([2,1]) == [1, -2]

func TestPslqSimple(t *testing.T) {
	prec := uint(64)
	//one := float64(1<<60)

	in := make([]big.Float, 2)
	in[0].SetPrec(prec).SetInt64(1)
	in[1].SetPrec(prec).SetInt64(-2)

	out, err := Pslq(prec, in, nil, 0, verbose)
	if err != nil {
		t.Error("Got error", err)
	}
	if verbose {
		printBigIntVector(t, "out", out)
	}
	compareResult(t, out, 2, 1)
}

func TestPslq2(t *testing.T) {
	prec := uint(64)

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
	out, err := Pslq(prec, in, nil, 0, verbose)
	if err != nil {
		t.Error("Got error", err)
	}
	if verbose {
		printBigIntVector(t, "out", out)
	}
	compareResult(t, out, 7, -21, -4, 0)
}

func TestPslq3(t *testing.T) {
	prec := uint(64)

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
	out, err := Pslq(prec, in, nil, 1000, verbose)
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

func TestPslq4(t *testing.T) {
	prec := uint(64)

	in := make([]big.Float, 8)
	for i := range in {
		if i == 0 {
			in[i].SetPrec(prec).SetFloat64(math.Pi)
		} else {
			bbp(prec, 16, 8, int64(i), &in[i])
		}
		if verbose {
			t.Logf("in[%d] = %g\n", i, &in[i])
		}
	}
	out, err := Pslq(prec, in, nil, 1000, verbose)
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
	out, err := Pslq(prec, in, maxCoeff, iterations, verbose)
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
	}

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		Pslq(prec, in, maxCoeff, iterations, verbose)
	}
}

func BenchmarkPslqPiFail128(b *testing.B)  { piFailBench(b, 128, 150, 3) }
func BenchmarkPslqPiFail256(b *testing.B)  { piFailBench(b, 256, 1E8, 10) }
func BenchmarkPslqPiFail512(b *testing.B)  { piFailBench(b, 512, 1E8, 22) }
func BenchmarkPslqPiFail1024(b *testing.B) { piFailBench(b, 1024, 1E8, 46) }

// FIXME try bpp base 100 as well? and 16?

// out = [44 120 -359 -665 431 -138 248 -166 146 -22 -5 20 339 -563 -606 -89 391 201 351 -31 -5 588 235 -663 183 646 -130 -73 11 167 -31 -788 666 -645 580 -15 -145 -523 -519 532 -169 686 43 80 -387 -234 560 486 285 -318]

func XXXXTestPslq4b(t *testing.T) {
	prec := uint(1024 * 2)

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
	out, err := Pslq(prec, in, big.NewInt(1E18), 1E6, verbose)
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

func TestPslq5(t *testing.T) {
	prec := uint(64)

	in := make([]big.Float, 8)
	in[0].SetPrec(prec).SetFloat64(math.Pi)
	acot(prec, 2, &in[1])
	acot(prec, 4, &in[2])
	acot(prec, 6, &in[3])
	acot(prec, 7, &in[4])
	acot(prec, 8, &in[5])
	acot(prec, 9, &in[6])
	acot(prec, 10, &in[7])
	out, err := Pslq(prec, in, nil, 1000, verbose)
	if err != nil {
		t.Error("Got error", err)
	}
	if verbose {
		printBigIntVector(t, "out", out)
	}
	compareResult(t, out, 1, -8, 0, 0, 4, 0, 0, 0)
}

func TestPslq6(t *testing.T) {
	prec := uint(64)

	in := make([]big.Float, 3)
	in[0].SetPrec(prec).SetFloat64(math.Pi / 4)
	acot(prec, 5, &in[1])
	acot(prec, 239, &in[2])
	out, err := Pslq(prec, in, nil, 1000, verbose)
	if err != nil {
		t.Error("Got error", err)
	}
	if verbose {
		printBigIntVector(t, "out", out)
	}
	compareResult(t, out, 1, -4, 1)
}

func TestPslq7(t *testing.T) {
	prec := uint(64)

	in := make([]big.Float, 5)
	in[0].SetPrec(prec).SetFloat64(math.Pi / 4)
	acot(prec, 49, &in[1])
	acot(prec, 57, &in[2])
	acot(prec, 239, &in[3])
	acot(prec, 110443, &in[4])
	out, err := Pslq(prec, in, nil, 1000, verbose)
	if err != nil {
		t.Error("Got error", err)
	}
	if verbose {
		printBigIntVector(t, "out", out)
	}
	compareResult(t, out, 1, -12, -32, 5, -12)
}

func TestPslq8(t *testing.T) {
	prec := uint(256)

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
	out, err := Pslq(prec, in, nil, 10000, verbose)
	if err != nil {
		t.Error("Got error", err)
	}
	if verbose {
		printBigIntVector(t, "out", out)
	}
	compareResult(t, out, 1, -4, 0, 0, 2, 1, 1)
}
