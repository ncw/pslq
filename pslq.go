// Implements PSLQ algorithm for integer relation detection.
//
// This code was originally ported from the sympy identification.py module to Go
//
// It was subsequently modified to improve correctness (in particular
// implementing the termination condition using the A matrix) and
// performance
package pslq

// Original code: Copyright (c) 2006-2014 SymPy Development Team
// Modifications: Copyright (c) 2014-2015 Nick Craig-Wood

import (
	"errors"
	"fmt"
	"log"
	"math"
	"math/big"
)

// Set to true to print a lot of debugging info
const debug = false

// Errors
var (
	ErrorPrecisionExhausted    = errors.New("precision exhausted")
	ErrorBadArguments          = errors.New("bad arguments: need at least 2 items")
	ErrorPrecisionTooLow       = errors.New("precision of input is too low")
	ErrorToleranceRoundsToZero = errors.New("tolerance is zero")
	ErrorZeroArguments         = errors.New("all input numbers must be non zero")
	ErrorArgumentTooSmall      = errors.New("one or more arguments are too small")
	ErrorNoRelationFound       = errors.New("could not find an integer relation")
	ErrorIterationsExceeded    = errors.New("ran out of iterations looking for relation")
)

// Pslq environment - may be reused and used concurrently
type Pslq struct {
	prec        uint
	target      uint
	tol         big.Float
	maxcoeff    big.Int
	maxcoeff_fp big.Float
	maxsteps    int
	verbose     bool
	one         big.Float
	half        big.Float
}

func max(a, b int) int {
	if a >= b {
		return a
	}
	return b
}

func min(a, b int) int {
	if a <= b {
		return a
	}
	return b
}

// Make a new matrix with that many rows and that many cols
func newMatrix(rows, cols int) [][]big.Float {
	M := make([][]big.Float, rows)
	for i := 0; i < cols; i++ {
		M[i] = make([]big.Float, cols)
	}
	return M
}

// Make a new matrix with that many rows and that many cols
func newBigIntMatrix(rows, cols int) [][]big.Int {
	M := make([][]big.Int, rows)
	for i := 0; i < cols; i++ {
		M[i] = make([]big.Int, cols)
	}
	return M
}

// Print a matrix
func printMatrix(name string, X [][]big.Float) {
	n := len(X) - 1
	for i := 1; i <= n; i++ {
		for j := 1; j <= n; j++ {
			fmt.Printf("%s[%d,%d] = %f\n", name, i, j, &X[i][j])
		}
		fmt.Printf("\n")
	}
}

// Print a matrix
func printBigIntMatrix(name string, X [][]big.Int) {
	n := len(X) - 1
	for i := 1; i <= n; i++ {
		for j := 1; j <= n; j++ {
			fmt.Printf("%s[%d,%d] = %d\n", name, i, j, &X[i][j])
		}
		fmt.Printf("\n")
	}
}

// Print a vector
func printVector(name string, x []big.Float) {
	for i := range x {
		if i == 0 {
			continue
		}
		fmt.Printf("%s[%d] = %f\n", name, i, &x[i])
	}
}

// Create a new environment for evaluating Pslq at the given
// precision.  This can be re-used multiple times and used
// concurrently after it has been set up.
func New(Prec uint) *Pslq {
	e := &Pslq{
		prec:     Prec,
		maxsteps: 100,
	}
	e.one.SetPrec(Prec).SetFloat64(1)
	e.half.SetPrec(Prec).SetFloat64(0.5)
	e.SetMaxCoeff(big.NewInt(1000))
	e.SetTarget((Prec * 3) / 4)
	return e
}

// SetVerbose if passed a true parameter then Run will log its
// progress
func (e *Pslq) SetVerbose(verbose bool) *Pslq {
	e.verbose = verbose
	return e
}

// SetMaxCoeff sets the maximum size of the parameter to be searched
// for - default 1000
func (e *Pslq) SetMaxCoeff(maxcoeff *big.Int) *Pslq {
	e.maxcoeff.Set(maxcoeff)
	e.maxcoeff_fp.SetPrec(e.prec).SetInt(&e.maxcoeff)
	return e
}

// SetMaxSteps sets the maximum number of steps of the algorithm to be
// run default 100
func (e *Pslq) SetMaxSteps(maxsteps int) *Pslq {
	e.maxsteps = maxsteps
	return e
}

// SetTarget sets the target precision of the result.
//
// By default this is 3/4 of the precision
func (e *Pslq) SetTarget(target uint) *Pslq {
	e.target = target
	e.tol.SetPrec(e.prec).SetMantExp(&e.one, -int(e.target))
	return e
}

// Compute the square root of n using Newton's Method. We start with
// an initial estimate for sqrt(n), and then iterate
//     x_{i+1} = 1/2 * ( x_i + (n / x_i) )
// Result is returned in x
func (e *Pslq) Sqrt(n, x *big.Float) {
	if n == x {
		panic("need distinct input and output")
	}
	if n.Sign() == 0 {
		x.Set(n)
		return
	} else if n.Sign() < 0 {
		panic("Sqrt of negative number")
	}
	prec := n.Prec()

	// Use the floating point square root as initial estimate
	nFloat64, _ := n.Float64()
	x.SetPrec(prec).SetFloat64(math.Sqrt(nFloat64))

	// We use t as a temporary variable. There's no need to set its precision
	// since big.Float values with unset (== 0) precision automatically assume
	// the largest precision of the arguments when used as the result (receiver)
	// of a big.Float operation.
	var t big.Float

	// Iterate.
	for {
		t.Quo(n, x)        // t = n / x_i
		t.Add(x, &t)       // t = x_i + (n / x_i)
		t.Mul(&e.half, &t) // x_{i+1} = 0.5 * t
		if x.Cmp(&t) == 0 {
			// Exit loop if no change to result
			break
		}
		x.Set(&t)
	}
}

// Sets res to nearest_int(x)
func (e *Pslq) NearestInt(x *big.Float, res *big.Int) {
	prec := x.Prec()
	var tmp big.Float
	tmp.SetPrec(prec)
	if x.Sign() >= 0 {
		tmp.Add(x, &e.half)
	} else {
		tmp.Sub(x, &e.half)
	}
	tmp.Int(res)
}

// Given a vector of real numbers x = [x_0, x_1, ..., x_n], Pslq(x)
// uses the PSLQ algorithm to find a list of integers
// [c_0, c_1, ..., c_n] such that
//
//     |c_1 x_1 + c_2 x_2 + ... + c_n x_n| < tolerance
//
// and such that max |c_k| < maxcoeff. If no such vector exists, Pslq
// returns one of the errors in this package depending on whether it
// has run out of iterations, precision or explored up to the
// maxcoeff. The tolerance defaults to 3/4 of the precision of the
// Environment passed in.
//
// This is a fairly direct translation of the pseudocode given by
// David Bailey, "The PSLQ Integer Relation Algorithm":
// http://www.cecm.sfu.ca/organics/papers/bailey/paper/html/node3.html
//
// This implementation uses fixed-point instead of floating-point
// arithmetic, since this is significantly faster.  This does
// introduce some additional failure modes over the original paper
// which should hopefully be covered correctly.
//
// If a result is returned, the first non-zero element will be positive
func (e *Pslq) Run(x []big.Float) ([]big.Int, error) {
	n := len(x)
	if n <= 1 {
		return nil, ErrorBadArguments
	}

	// At too low precision, the algorithm becomes meaningless
	if e.prec < 64 {
		return nil, ErrorPrecisionTooLow
	}

	if e.verbose && int(e.prec)/max(2, int(n)) < 5 {
		log.Printf("Warning: precision for PSLQ may be too low")
	}

	if e.verbose {
		log.Printf("PSLQ using prec %d and tol %g", e.prec, e.tol)
	}

	if e.tol.Sign() == 0 {
		return nil, ErrorToleranceRoundsToZero
	}

	// Temporary variables
	tmp0 := new(big.Float).SetPrec(e.prec)
	tmp1 := new(big.Float).SetPrec(e.prec)
	bigTmp := new(big.Int)

	// Convert to use 1-based indexing to allow us to be
	// consistent with Bailey's indexing.
	xNew := make([]big.Float, len(x)+1)
	minx := new(big.Float).SetPrec(e.prec)
	minxFirst := true
	for i, xk := range x {
		p := &xNew[i+1]
		p.Set(&xk)
		tmp0.Abs(p)
		if minxFirst || tmp0.Cmp(minx) < 0 {
			minxFirst = false
			minx.Set(tmp0)
		}
	}
	x = xNew
	if debug {
		printVector("x", x)
	}

	// Sanity check on magnitudes
	if minx.Sign() == 0 {
		return nil, ErrorZeroArguments
	}
	tmp1.SetInt64(128)
	tmp0.Quo(&e.tol, tmp1)
	if minx.Cmp(tmp0) < 0 { //  minx < tol/128
		return nil, ErrorArgumentTooSmall
	}

	tmp0.SetInt64(4)
	tmp1.SetInt64(3)
	tmp0.Quo(tmp0, tmp1)
	var γ big.Float
	e.Sqrt(tmp0, &γ) // sqrt(4<<prec)/3)
	if debug {
		fmt.Printf("γ = %f\n", &γ)
	}
	A := newBigIntMatrix(n+1, n+1)
	B := newBigIntMatrix(n+1, n+1)
	H := newMatrix(n+1, n+1)
	// Initialization Step 1
	//
	// Set the n×n matrices A and B to the identity.
	for i := 1; i <= n; i++ {
		for j := 1; j <= n; j++ {
			if i == j {
				A[i][j].SetInt64(1)
				B[i][j].SetInt64(1)
			} else {
				A[i][j].SetInt64(0)
				B[i][j].SetInt64(0)
			}
			H[i][j].SetInt64(0)
		}
	}
	if debug {
		printBigIntMatrix("A", A)
		printBigIntMatrix("B", B)
		printMatrix("H", H)
	}
	// Initialization Step 2
	//
	// For k := 1 to n
	//     compute s_k := sqrt( sum_j=k^n x_j^2 )
	// endfor.
	// Set t = 1/s1.
	// For k := 1 to n:
	//     y_k := t * x_k
	//     s_k := t * s_k
	// endfor.
	s := make([]big.Float, n+1)
	for i := 1; i <= n; i++ {
		s[i].SetInt64(0)
	}
	for k := 1; k <= n; k++ {
		var t big.Float
		t.SetInt64(0)
		for j := k; j <= n; j++ {
			tmp0.Mul(&x[j], &x[j])
			t.Add(&t, tmp0)
		}
		e.Sqrt(&t, &s[k])
	}
	if debug {
		fmt.Println("Init Step 2")
		printVector("s", s)
	}
	var t big.Float
	t.Set(&s[1])
	y := make([]big.Float, len(x))
	copy(y, x)
	for k := 1; k <= n; k++ {
		// y[k] = (x[k] << prec) / t
		y[k].Quo(&x[k], &t)
		// s[k] = (s[k] << prec) / t
		s[k].Quo(&s[k], &t)
	}
	if debug {
		printVector("y", y)
		printVector("s", s)
	}
	// Init Step 3
	//
	// Compute the n×(n−1) matrix H as follows:
	// For i := 1 to n:
	//     for j := i + 1 to n − 1:
	//         set Hij := 0
	//     endfor
	//     if i ≤ n − 1 then set Hii := s_(i+1)/s_i
	//     for j := 1 to i−1:
	//         set Hij := −y_i * y_j / (s_j * s_(j+1))
	//     endfor
	// endfor
	for i := 1; i <= n; i++ {
		for j := i + 1; j < n; j++ {
			H[i][j].SetInt64(0)
		}
		if i <= n-1 {
			if s[i].Sign() == 0 {
				// Precision probably exhausted
				return nil, ErrorPrecisionExhausted
			}
			// H[i][i] = (s[i+1] << prec) / s[i]
			H[i][i].Quo(&s[i+1], &s[i])
		}
		for j := 1; j < i; j++ {
			var sjj1 big.Float
			sjj1.Mul(&s[j], &s[j+1])
			if debug {
				fmt.Printf("sjj1 = %f\n", &sjj1)
			}
			if sjj1.Sign() == 0 {
				// Precision probably exhausted
				return nil, ErrorPrecisionExhausted
			}
			// H[i][j] = ((-y[i] * y[j]) << prec) / sjj1
			tmp0.Mul(&y[i], &y[j])
			tmp0.Neg(tmp0)
			H[i][j].Quo(tmp0, &sjj1)
		}
	}
	if debug {
		fmt.Println("Init Step 3")
		printMatrix("H", H)
	}
	// Init Step 4
	//
	// Perform full reduction on H, simultaneously updating y, A and B:
	//
	// For i := 2 to n:
	//     for j := i−1 to 1 step−1:
	//         t := nint(Hij/Hjj)
	//         y_j := y_j + t * y_i
	//         for k := 1 to j:
	//             Hik := Hik − t * Hjk
	//         endfor
	//         for k := 1 to n:
	//             Aik := Aik − t * Ajk
	//             Bkj := Bkj + t * Bki
	//         endfor
	//     endfor
	// endfor
	for i := 2; i <= n; i++ {
		for j := i - 1; j > 0; j-- {
			//t = floor(H[i][j]/H[j,j] + 0.5)
			var t big.Int
			var tFloat big.Float
			if H[j][j].Sign() == 0 {
				// Precision probably exhausted
				return nil, ErrorPrecisionExhausted
			}
			tmp0.Quo(&H[i][j], &H[j][j])
			e.NearestInt(tmp0, &t)
			tFloat.SetInt(&t).SetPrec(e.prec)
			if debug {
				fmt.Printf("H[i][j]=%f\n", &H[i][j])
				fmt.Printf("H[j][j]=%f\n", &H[j][j])
				fmt.Printf("tmp=%f\n", tmp0)
				fmt.Printf("t=%d\n", &t)
			}
			// y[j] = y[j] + (t * y[i] >> prec)
			tmp0.Mul(&y[i], &tFloat)
			y[j].Add(&y[j], tmp0)
			for k := 1; k <= j; k++ {
				// H[i][k] = H[i][k] - (t * H[j][k] >> prec)
				tmp0.Mul(&H[j][k], &tFloat)
				H[i][k].Sub(&H[i][k], tmp0)
			}
			for k := 1; k <= n; k++ {
				bigTmp.Mul(&t, &A[j][k])
				A[i][k].Sub(&A[i][k], bigTmp)
				bigTmp.Mul(&t, &B[k][i])
				B[k][j].Add(&B[k][j], bigTmp)
			}
		}
	}
	if debug {
		fmt.Println("Init Step 4")
		printBigIntMatrix("A", A)
		printBigIntMatrix("B", B)
		printMatrix("H", H)
	}
	// Main algorithm
	var REP int
	var norm big.Int
	vec := make([]big.Int, n)
	for REP = 0; REP < e.maxsteps; REP++ {
		// Step 1
		//
		// Select m such that γ^i * |Hii| is maximal when i = m.
		m := -1
		var szmax big.Float
		szmax.SetInt64(-1)
		var γPower big.Float
		γPower.Set(&γ)
		for i := 1; i < n; i++ {
			var absH big.Float
			absH.Abs(&H[i][i])
			var sz big.Float
			sz.Mul(&γPower, &absH)
			// sz := (g**i * abs(h)) >> (prec * (i - 1))
			if sz.Cmp(&szmax) > 0 {
				m = i
				szmax.Set(&sz)
			}
			γPower.Mul(&γPower, &γ)
		}
		if debug {
			fmt.Println("Step 1")
			fmt.Printf("szmax=%f\n", &szmax)
			fmt.Printf("m=%d\n", m)
		}
		// Step 2
		//
		// Exchange entries m and m+1 of y, corresponding rows
		// of A and H, and corresponding columns of B.
		y[m], y[m+1] = y[m+1], y[m]
		for i := 1; i < n+1; i++ {
			H[m][i], H[m+1][i] = H[m+1][i], H[m][i]
		}
		for i := 1; i < n+1; i++ {
			A[m][i], A[m+1][i] = A[m+1][i], A[m][i]
		}
		for i := 1; i < n+1; i++ {
			B[i][m], B[i][m+1] = B[i][m+1], B[i][m]
		}
		if debug {
			fmt.Println("Step 2")
			printVector("y", y)
			printBigIntMatrix("A", A)
			printBigIntMatrix("B", B)
			printMatrix("H", H)
		}
		// Step 3
		//
		// If m ≤ n−2 then update H as follows:
		//
		// t0 := sqrt( Hmm^2 + H(m,m+1)^2 )
		// t1 := Hmm/t0
		// t2 := H(m,m+1)/t0.
		// for i := m to n:
		//     t3 := Him
		//     t4 := Hi,m+1
		//     Him := t1t3 +t2t4
		//     Hi,m+1 := −t2t3 +t1t4
		// endfor.
		if m <= n-2 {
			tmp0.Mul(&H[m][m], &H[m][m])
			tmp1.Mul(&H[m][m+1], &H[m][m+1])
			tmp0.Add(tmp0, tmp1)
			var t0 big.Float
			e.Sqrt(tmp0, &t0)
			// Precision probably exhausted
			if t0.Sign() == 0 {
				return nil, ErrorPrecisionExhausted
			}
			var t1, t2 big.Float
			t1.Quo(&H[m][m], &t0)
			t2.Quo(&H[m][m+1], &t0)
			for i := m; i <= n; i++ {
				var t3, t4 big.Float
				t3.Set(&H[i][m])
				t4.Set(&H[i][m+1])
				// H[i][m] = (t1*t3 + t2*t4) >> prec
				tmp0.Mul(&t1, &t3)
				tmp1.Mul(&t2, &t4)
				H[i][m].Add(tmp0, tmp1)
				// H[i][m+1] = (-t2*t3 + t1*t4) >> prec
				tmp0.Mul(&t2, &t3)
				tmp1.Mul(&t1, &t4)
				H[i][m+1].Sub(tmp1, tmp0)
			}
		}
		if debug {
			fmt.Println("Step 3")
			printMatrix("H", H)
		}
		// Step 4
		// Perform block reduction on H, simultaneously updating y, A and B:
		//
		// For i := m+1 to n:
		//     for j := min(i−1, m+1) to 1 step −1:
		//         t := nint(Hij/Hjj)
		//         yj := yj + t * yi
		//         for k := 1 to j:
		//             Hik := Hik − tHjk
		//         endfor
		//         for k := 1 to n:
		//             Aik := Aik −tAjk
		//             Bkj := Bkj +tBki
		//         endfor
		//     endfor
		// endfor.
		for i := m + 1; i <= n; i++ {
			var t big.Int
			var tFloat big.Float
			for j := min(i-1, m+1); j > 0; j-- {
				if H[j][j].Sign() == 0 {
					// Precision probably exhausted
					return nil, ErrorPrecisionExhausted
				}
				tmp0.Quo(&H[i][j], &H[j][j])
				e.NearestInt(tmp0, &t)
				tFloat.SetInt(&t).SetPrec(e.prec)
				// y[j] = y[j] + ((t * y[i]) >> prec)
				tmp0.Mul(&y[i], &tFloat)
				y[j].Add(&y[j], tmp0)
				for k := 1; k <= j; k++ {
					// H[i][k] = H[i][k] - (t * H[j][k] >> prec)
					tmp0.Mul(&H[j][k], &tFloat)
					H[i][k].Sub(&H[i][k], tmp0)
				}
				for k := 1; k <= n; k++ {
					bigTmp.Mul(&t, &A[j][k])
					A[i][k].Sub(&A[i][k], bigTmp)
					bigTmp.Mul(&t, &B[k][i])
					B[k][j].Add(&B[k][j], bigTmp)
				}
			}
		}
		if debug {
			fmt.Println("Step 4")
			printBigIntMatrix("A", A)
			printBigIntMatrix("B", B)
			printMatrix("H", H)
		}

		// Step 6
		//
		// Termination test: If the largest entry of A exceeds
		// the level of numeric precision used, then precision
		// is exhausted. If the smallest entry of the y vector
		// is less than the detection threshold, a relation
		// has been detected and is given in the corresponding
		// column of B.
		//
		// Until a relation is found, the error typically decreases
		// slowly (e.g. a factor 1-10) with each step TODO: we could
		// compare err from two successive iterations. If there is a
		// large drop (several orders of magnitude), that indicates a
		// "high quality" relation was detected. Reporting this to
		// the user somehow might be useful.

		maxAPrecision := 0
		for i := 1; i <= n; i++ {
			for j := 1; j <= n; j++ {
				precision := A[i][j].BitLen()
				if precision > maxAPrecision {
					maxAPrecision = precision
				}
			}
		}
		if debug {
			log.Printf("Max A precision = %d, precision = %d, tolerance %d, ratio = %.3f\n", maxAPrecision, e.prec, e.target, float64(maxAPrecision)/float64(e.target))
		}
		if float64(maxAPrecision)/float64(e.target) > 0.85 {
			if e.verbose {
				log.Printf("CANCELLING after step %d/%d.", REP, e.maxsteps)
			}
			return nil, ErrorPrecisionExhausted
		}

		var best_err big.Float
		best_err.Set(&e.maxcoeff_fp)
		for i := 1; i <= n; i++ {
			var err big.Float
			err.Abs(&y[i])
			// Maybe we are done?
			if err.Cmp(&e.tol) < 0 {
				// We are done if the coefficients are acceptable
				var maxc big.Int
				for j := 1; j <= n; j++ {
					if debug {
						fmt.Printf("vec[%d]=%d\n", j-1, &B[j][i])
					}
					t := B[j][i]
					if debug {
						fmt.Printf("vec[%d]=%d\n", j-1, t)
					}
					vec[j-1] = t
					if t.Sign() < 0 {
						t.Neg(&t)
					}
					if t.Cmp(&maxc) > 0 {
						maxc.Set(&t)
					}
				}
				if debug {
					fmt.Printf("maxc = %d, maxcoeff = %d\n", maxc, e.maxcoeff)
				}
				if maxc.Cmp(&e.maxcoeff) < 0 {
					if e.verbose {
						log.Printf("FOUND relation at iter %d/%d, error: %g", REP, e.maxsteps, &err)
					}
					// Find sign of first non zero item
					sign := 0
					for i := range vec {
						sign = vec[i].Sign()
						if sign != 0 {
							break
						}
					}
					// Normalise vec making first non-zero argument positive
					if sign < 0 {
						for i := range vec {
							vec[i].Neg(&vec[i])
						}
					}
					return vec, nil
				}
			}
			if err.Cmp(&best_err) < 0 {
				best_err = err
			}
		}
		// Step 5
		//
		// Norm bound: Compute M := 1/maxj |Hj|, where Hj
		// denotes the j-th row of H.
		//
		// Then there can exist no relation vector whose
		// Euclidean norm is less than M.
		//
		// Calculate a lower bound for the norm. We could do this
		// more exactly (using the Euclidean norm) but there is probably
		// no practical benefit.
		var recnorm big.Float
		recnorm.SetInt64(0)
		for i := 1; i <= n; i++ {
			for j := 1; j <= n; j++ {
				tmp0.Abs(&H[i][j])
				if tmp0.Cmp(&recnorm) > 0 {
					recnorm.Set(tmp0)
				}
			}
		}
		norm.Set(&e.maxcoeff)
		if recnorm.Sign() != 0 {
			// norm = ((1 << (2 * prec)) / recnorm) >> prec
			tmp0.Quo(&e.one, &recnorm)
			tmp0.Int(&norm)
		}
		if e.verbose {
			log.Printf("%2d/%2d:  Error: %g   Norm: %d", REP, e.maxsteps, &best_err, &norm)
		}
		if norm.Cmp(&e.maxcoeff) >= 0 {
			if e.verbose {
				log.Printf("CANCELLING after step %d/%d.", REP, e.maxsteps)
				log.Printf("Could not find an integer relation. Norm bound: %d", &norm)
			}
			return nil, ErrorNoRelationFound
		}
	}
	if e.verbose {
		log.Printf("CANCELLING after step %d/%d.", REP, e.maxsteps)
		log.Printf("Could not find an integer relation. Norm bound: %d", &norm)
	}
	return nil, ErrorIterationsExceeded
}
