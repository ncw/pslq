// Implens PSLQ algorithm for iteger releation detection.
package pslq

import (
	"errors"
	"log"

	fp "fixedpoint"
	"fmt"
)

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
func newMatrix(rows, cols int) [][]fp.FixedPoint {
	M := make([][]fp.FixedPoint, rows)
	for i := 0; i < cols; i++ {
		M[i] = make([]fp.FixedPoint, cols)
	}
	return M
}

// Print a matrix
func printMatrix(name string, X [][]fp.FixedPoint) {
	n := len(X) - 1
	for i := 1; i <= n; i++ {
		for j := 1; j <= n; j++ {
			fmt.Printf("%s[%d,%d] = %d\n", name, i, j, &X[i][j])
		}
		fmt.Printf("\n")
	}
}

// Print a vector
func printVector(name string, x []fp.FixedPoint) {
	for i := range x {
		if i == 0 {
			continue
		}
		fmt.Printf("%s[%d] = %d\n", name, i, &x[i])
	}
}

// Given a vector of real numbers `x = [x_0, x_1, ..., x_n]`, ``pslq(x)``
// uses the PSLQ algorithm to find a list of integers
// `[c_0, c_1, ..., c_n]` such that
//
// .. math ::
//
//     |c_1 x_1 + c_2 x_2 + ... + c_n x_n| < \mathrm{tol}
//
// and such that `\max |c_k| < \mathrm{maxcoeff}`. If no such vector
// exists, :func:`~mpmath.pslq` returns ``None``. The tolerance defaults to
// 3/4 of the working precision.
//
// **Algorithm**
//
// This is a fairly direct translation to Python of the pseudocode given by
// David Bailey, "The PSLQ Integer Relation Algorithm":
// http://www.cecm.sfu.ca/organics/papers/bailey/paper/html/node3.html
//
// The present implementation uses fixed-point instead of floating-point
// arithmetic, since this is significantly (about 7x) faster.
//
// prec is the number of bits of precision each fp.FixedPoint has
func Pslq(startEnv *fp.Environment, x []fp.FixedPoint, maxcoeff int64, maxsteps int, verbose bool) ([]int64, error) {
	if maxcoeff == 0 {
		maxcoeff = 1000
	}
	if maxsteps == 0 {
		maxsteps = 100
	}

	n := len(x)
	if n <= 1 {
		return nil, errors.New("need at least 2 items")
	}

	// At too low precision, the algorithm becomes meaningless
	if startEnv.Prec < 53 {
		return nil, errors.New("prec too low")
	}

	if verbose && int(startEnv.Prec)/max(2, int(n)) < 5 {
		log.Printf("Warning: precision for PSLQ may be too low")
	}

	target := (startEnv.Prec * 3) / 4

	extra := uint(60) // was 60 but made it 2**n FIXME make 64 again
	env := fp.NewEnvironment(startEnv.Prec + extra)

	// FIXME make it so you can pass tol in
	tol := env.NewInt(1)
	tol.Rsh(tol, target)

	if verbose {
		log.Printf("PSLQ using prec %d and tol %s", env.Prec, tol)
	}

	if tol.Sign() == 0 {
		return nil, errors.New("tol must be > 0 here")
	}

	// Useful constants
	_1 := env.NewInt(1)
	_100 := env.NewInt(100)

	// Temporary variables
	tmp0 := env.New()
	tmp1 := env.New()

	// Convert to fixed-point numbers. The dummy None is added so we can
	// use 1-based indexing. (This just allows us to be consistent with
	// Bailey's indexing. The algorithm is 100 lines long, so debugging
	// a single wrong index can be painful.)
	xNew := make([]fp.FixedPoint, len(x)+1)
	minx := env.New()
	minxFirst := true
	for i, xk := range x {
		p := &xNew[i+1]
		p.Convert(env, &xk)
		tmp0.Abs(p)
		if minxFirst || tmp0.Cmp(minx) < 0 {
			minxFirst = false
			minx.Set(tmp0)
		}
	}
	x = xNew
	printVector("x", x)

	// Sanity check on magnitudes
	if minx.Sign() == 0 {
		return nil, errors.New("PSLQ requires a vector of nonzero numbers")
	}
	tmp0.Div(tol, _100)
	if minx.Cmp(tmp0) < 0 { //  minx < tol/100
		return nil, errors.New("STOPPING: (one number is too small)")
	}

	tmp0.SetInt64(env, 4)
	tmp0.DivInt64(tmp0, 3)
	var g fp.FixedPoint
	g.Sqrt(tmp0) /// sqrt(4<<prec)/3)
	fmt.Printf("g = %d\n", &g)
	A := newMatrix(n+1, n+1)
	B := newMatrix(n+1, n+1)
	H := newMatrix(n+1, n+1)
	// Initialization
	// step 1
	for i := 1; i <= n; i++ {
		for j := 1; j <= n; j++ {
			if i == j {
				A[i][j].SetInt64(env, 1)
				B[i][j].SetInt64(env, 1)
			} else {
				A[i][j].SetInt64(env, 0)
				B[i][j].SetInt64(env, 0)
			}
			H[i][j].SetInt64(env, 0)
		}
	}
	printMatrix("A", A)
	printMatrix("B", B)
	printMatrix("H", H)
	// step 2
	s := make([]fp.FixedPoint, n+1)
	for i := 1; i <= n; i++ {
		s[i].SetInt64(env, 0)
	}
	for k := 1; k <= n; k++ {
		var t fp.FixedPoint
		t.SetInt64(env, 0)
		for j := k; j <= n; j++ {
			tmp0.Mul(&x[j], &x[j])
			t.Add(&t, tmp0)
		}
		s[k].Sqrt(&t)
	}
	printVector("s", s)
	var t fp.FixedPoint
	t.Set(&s[1])
	y := make([]fp.FixedPoint, len(x))
	copy(y, x)
	for k := 1; k <= n; k++ {
		// y[k] = (x[k] << prec) / t
		y[k].Div(&x[k], &t)
		// s[k] = (s[k] << prec) / t
		s[k].Div(&s[k], &t)
	}
	printVector("y", y)
	printVector("s", s)
	// step 3
	for i := 1; i <= n; i++ {
		for j := i + 1; j < n; j++ {
			H[i][j].SetInt64(env, 0)
		}
		if i <= n-1 {
			if s[i].Sign() != 0 {
				// H[i][i] = (s[i+1] << prec) / s[i]
				H[i][i].Div(&s[i+1], &s[i])
			} else {
				H[i][i].SetInt64(env, 0)
			}
		}
		for j := 1; j < i; j++ {
			var sjj1 fp.FixedPoint
			sjj1.Mul(&s[j], &s[j+1])
			fmt.Printf("sjj1 = %d\n", &sjj1)
			if sjj1.Sign() != 0 {
				// H[i][j] = ((-y[i] * y[j]) << prec) / sjj1
				tmp0.Mul(&y[i], &y[j])
				tmp0.Neg(tmp0)
				H[i][j].Div(tmp0, &sjj1)
			} else {
				H[i][j].SetInt64(env, 0)
			}
		}
	}
	printMatrix("H", H)
	// step 4
	for i := 2; i <= n; i++ {
		for j := i - 1; j > 0; j-- {
			//t = floor(H[i][j]/H[j,j] + 0.5)
			if H[j][j].Sign() != 0 {
				tmp0.Div(&H[i][j], &H[j][j])
				t.Round(tmp0)
			} else {
				//t = 0
				continue
			}
			// y[j] = y[j] + (t * y[i] >> prec)
			tmp0.Mul(&t, &y[i])
			y[j].Add(&y[j], tmp0)
			for k := 1; k <= j; k++ {
				// H[i][k] = H[i][k] - (t * H[j][k] >> prec)
				tmp0.Mul(&t, &H[j][k])
				H[i][k].Sub(&H[i][k], tmp0)
			}
			for k := 1; k <= n; k++ {
				// A[i][k] = A[i][k] - (t * A[j][k] >> prec)
				tmp0.Mul(&t, &A[j][k])
				A[i][k].Sub(&A[i][k], tmp0)
				// B[k][j] = B[k][j] + (t * B[k][i] >> prec)
				tmp0.Mul(&t, &B[k][i])
				B[k][j].Add(&B[k][j], tmp0)
			}
		}
	}
	// Main algorithm
	var REP int
	var norm int64
	vec := make([]int64, n) // FIXME big.Int?
	for REP = 0; REP < maxsteps; REP++ {
		// Step 1
		m := -1
		var szmax fp.FixedPoint
		szmax.SetInt64(env, -1)
		var gPower fp.FixedPoint
		gPower.Set(&g)
		for i := 1; i < n; i++ {
			h := H[i][i]
			var absH fp.FixedPoint
			absH.Abs(&h)
			var sz fp.FixedPoint
			sz.Mul(&gPower, &absH)
			// sz := (g**i * abs(h)) >> (prec * (i - 1))
			if sz.Cmp(&szmax) > 0 {
				m = i
				szmax.Set(&sz)
			}
			gPower.Mul(&gPower, &g)
		}
		// Step 2
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
		// Step 3
		if m <= n-2 {
			tmp0.Mul(&H[m][m], &H[m][m])
			tmp1.Mul(&H[m][m+1], &H[m][m+1])
			tmp0.Add(tmp0, tmp1)
			var t0 fp.FixedPoint
			t0.Sqrt(tmp0)
			// A zero element probably indicates that the precision has
			// been exhausted. XXX: this could be spurious, due to
			// using fixed-point arithmetic
			if t0.Sign() == 0 {
				break
			}
			var t1, t2 fp.FixedPoint
			t1.Div(&H[m][m], &t0)
			t2.Div(&H[m][m+1], &t0)
			for i := m; i <= n; i++ {
				var t3, t4 fp.FixedPoint
				t3 = H[i][m]
				t4 = H[i][m+1]
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
		// Step 4
		for i := m + 1; i <= n; i++ {
			for j := min(i-1, m+1); j > 0; j-- {
				if H[j][j].Sign() == 0 {
					// Precision probably exhausted
					break
				}
				tmp0.Div(&H[i][j], &H[j][j])
				t.Round(tmp0)
				// y[j] = y[j] + ((t * y[i]) >> prec)
				tmp0.Mul(&t, &y[i])
				y[j].Add(&y[j], tmp0)
				for k := 1; k <= j; k++ {
					// H[i][k] = H[i][k] - (t * H[j][k] >> prec)
					tmp0.Mul(&t, &H[j][k])
					H[i][k].Sub(&H[i][k], tmp0)
				}
				for k := 1; k <= n; k++ {
					// A[i][k] = A[i][k] - (t * A[j][k] >> prec)
					tmp0.Mul(&t, &A[j][k])
					A[i][k].Sub(&A[i][k], tmp0)
					// B[k][j] = B[k][j] + (t * B[k][i] >> prec)
					tmp0.Mul(&t, &B[k][i])
					B[k][j].Add(&B[k][j], tmp0)
				}
			}
		}
		// Until a relation is found, the error typically decreases
		// slowly (e.g. a factor 1-10) with each step TODO: we could
		// compare err from two successive iterations. If there is a
		// large drop (several orders of magnitude), that indicates a
		// "high quality" relation was detected. Reporting this to
		// the user somehow might be useful.
		var best_err fp.FixedPoint
		best_err.SetInt64(env, int64(maxcoeff))
		for i := 1; i <= n; i++ {
			var err fp.FixedPoint
			err.Abs(&y[i])
			// Maybe we are done?
			if err.Cmp(tol) < 0 {
				// We are done if the coefficients are acceptable
				// FIXME use big.Int here?
				maxc := int64(0)
				for j := 1; i <= n; i++ {
					vec[j-1] = B[j][i].RoundInt64()
					t := vec[j-1]
					if t < 0 {
						t = -t
					}
					if t > maxc {
						maxc = t
					}
				}
				if maxc < maxcoeff {
					if verbose {
						log.Printf("FOUND relation at iter %d/%d, error: %d", REP, maxsteps, &err)
					}
					return vec, nil
				}
			}
			if err.Cmp(&best_err) < 0 {
				best_err = err
			}
		}
		// Calculate a lower bound for the norm. We could do this
		// more exactly (using the Euclidean norm) but there is probably
		// no practical benefit.
		var recnorm fp.FixedPoint
		recnorm.SetInt64(env, 0)
		for i := 1; i <= n; i++ {
			for j := 1; j <= n; j++ {
				tmp0.Abs(&H[i][j])
				if tmp0.Cmp(&recnorm) > 0 {
					recnorm.Set(tmp0)
				}
			}
		}
		norm = maxcoeff
		if recnorm.Sign() != 0 {
			// norm = ((1 << (2 * prec)) / recnorm) >> prec
			tmp0.Div(_1, &recnorm)
			norm = tmp0.Int64()
			norm /= 100
		}
		if verbose {
			log.Printf("%2d/%2d:  Error: %d   Norm: %d", REP, maxsteps, &best_err, norm)
		}
		if norm >= maxcoeff {
			break
		}
	}
	if verbose {
		log.Printf("CANCELLING after step %d/%d.", REP, maxsteps)
		log.Printf("Could not find an integer relation. Norm bound: %d", norm)
	}
	return nil, nil
}
