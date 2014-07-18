// Implens PSLQ algorithm for iteger releation detection.
package pslq

// Try raising internal precision so that we can eliminate false
// positives.  Using twice the input precision should emulate floating
// point arithmetic pretty well.

// Make an adaptive Pslq which raises precision

// Can we detect false positives?

// FIXME A isn't used in the result?
//
// Mising one of the terminaton tests: If the largest entry of A
// exceeds the level of numeric precision used, then precision is
// exhausted.

// FIXME where did / 100s come from?

import (
	"errors"
	"fmt"
	"fp"
	"log"
	"math/big"
)

const debug = false

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

// Make a new matrix with that many rows and that many cols
func newInt64Matrix(rows, cols int) [][]int64 {
	M := make([][]int64, rows)
	for i := 0; i < cols; i++ {
		M[i] = make([]int64, cols)
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

// Print a matrix
func printInt64Matrix(name string, X [][]int64) {
	n := len(X) - 1
	for i := 1; i <= n; i++ {
		for j := 1; j <= n; j++ {
			fmt.Printf("%s[%d,%d] = %d\n", name, i, j, X[i][j])
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
	_100big := big.NewInt(100)
	var maxcoeff_fp fp.FixedPoint
	maxcoeff_fp.SetInt64(env, maxcoeff)
	var maxcoeff_big big.Int
	maxcoeff_big.SetInt64(maxcoeff)

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
	if debug {
		printVector("x", x)
	}

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
	var γ fp.FixedPoint
	γ.Sqrt(tmp0) /// sqrt(4<<prec)/3)
	if debug {
		fmt.Printf("γ = %d\n", &γ)
	}
	A := newInt64Matrix(n+1, n+1)
	B := newInt64Matrix(n+1, n+1)
	H := newMatrix(n+1, n+1)
	// Initialization Step 1
	//
	// Set the n×n matrices A and B to the identity.
	for i := 1; i <= n; i++ {
		for j := 1; j <= n; j++ {
			if i == j {
				A[i][j] = 1
				B[i][j] = 1
			} else {
				A[i][j] = 0
				B[i][j] = 0
			}
			H[i][j].SetInt64(env, 0)
		}
	}
	if debug {
		printInt64Matrix("A", A)
		printInt64Matrix("B", B)
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
	if debug {
		fmt.Println("Init Step 2")
		printVector("s", s)
	}
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
			if debug {
				fmt.Printf("sjj1 = %d\n", &sjj1)
			}
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
			var t int64
			if H[j][j].Sign() != 0 {
				tmp0.Div(&H[i][j], &H[j][j])
				// FIXME div is 1 different in py - temp for matching up
				one := big.NewInt(1)         // FIXME
				tmp0.Int.Sub(&tmp0.Int, one) // FIXME
				t = tmp0.RoundInt64()
			} else {
				//t = 0
				continue
			}
			if debug {
				fmt.Printf("H[i][j]=%d\n", &H[i][j])
				fmt.Printf("H[j][j]=%d\n", &H[j][j])
				fmt.Printf("tmp=%d\n", &tmp0.Int)
				fmt.Printf("tmp=%d\n", tmp0)
				fmt.Printf("t=%d\n", t)
			}
			// y[j] = y[j] + (t * y[i] >> prec)
			tmp0.MulInt64(&y[i], t)
			y[j].Add(&y[j], tmp0)
			for k := 1; k <= j; k++ {
				// H[i][k] = H[i][k] - (t * H[j][k] >> prec)
				tmp0.MulInt64(&H[j][k], t)
				H[i][k].Sub(&H[i][k], tmp0)
			}
			for k := 1; k <= n; k++ {
				A[i][k] -= t * A[j][k]
				B[k][j] += t * B[k][i]
			}
		}
	}
	if debug {
		fmt.Println("Init Step 4")
		printInt64Matrix("A", A)
		printInt64Matrix("B", B)
		printMatrix("H", H)
	}
	// Main algorithm
	var REP int
	var norm big.Int
	vec := make([]int64, n)
	for REP = 0; REP < maxsteps; REP++ {
		// Step 1
		//
		// Select m such that γ^i * |Hii| is maximal when i = m.
		m := -1
		var szmax fp.FixedPoint
		szmax.SetInt64(env, -1)
		var γPower fp.FixedPoint
		γPower.Set(&γ)
		for i := 1; i < n; i++ {
			var absH fp.FixedPoint
			absH.Abs(&H[i][i])
			var sz fp.FixedPoint
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
			fmt.Printf("szmax=%d\n", &szmax)
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
			printInt64Matrix("A", A)
			printInt64Matrix("B", B)
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
			var t0 fp.FixedPoint
			t0.Sqrt(tmp0)
			// A zero element probably indicates that the precision has
			// been exhausted. FIXME: this could be spurious, due to
			// using fixed-point arithmetic
			if t0.Sign() == 0 {
				break
			}
			var t1, t2 fp.FixedPoint
			t1.Div(&H[m][m], &t0)
			t2.Div(&H[m][m+1], &t0)
			for i := m; i <= n; i++ {
				var t3, t4 fp.FixedPoint
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
			var t int64
			for j := min(i-1, m+1); j > 0; j-- {
				if H[j][j].Sign() == 0 {
					// Precision probably exhausted
					break
				}
				tmp0.Div(&H[i][j], &H[j][j])
				t = tmp0.RoundInt64()
				// y[j] = y[j] + ((t * y[i]) >> prec)
				tmp0.MulInt64(&y[i], t)
				y[j].Add(&y[j], tmp0)
				for k := 1; k <= j; k++ {
					// H[i][k] = H[i][k] - (t * H[j][k] >> prec)
					tmp0.MulInt64(&H[j][k], t)
					H[i][k].Sub(&H[i][k], tmp0)
				}
				for k := 1; k <= n; k++ {
					A[i][k] -= t * A[j][k]
					B[k][j] += t * B[k][i]
				}
			}
		}
		if debug {
			fmt.Println("Step 4")
			printInt64Matrix("A", A)
			printInt64Matrix("B", B)
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
		var best_err fp.FixedPoint
		best_err.Set(&maxcoeff_fp)
		for i := 1; i <= n; i++ {
			// FIXME what if there is more than one value
			// of |y| < tol?
			var err fp.FixedPoint
			err.Abs(&y[i])
			// Maybe we are done?
			if err.Cmp(tol) < 0 {
				// We are done if the coefficients are acceptable
				var maxc int64
				for j := 1; j <= n; j++ {
					if debug {
						fmt.Printf("vec[%d]=%d\n", j-1, &B[j][i])
					}
					t := B[j][i]
					if debug {
						fmt.Printf("vec[%d]=%d\n", j-1, t)
					}
					vec[j-1] = t
					if t < 0 {
						t = -t
					}
					if t > maxc {
						maxc = t
					}
				}
				if debug {
					fmt.Printf("maxc = %d, maxcoeff = %d\n", maxc, maxcoeff)
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
		norm.SetInt64(maxcoeff)
		if recnorm.Sign() != 0 {
			// norm = ((1 << (2 * prec)) / recnorm) >> prec
			tmp0.Div(_1, &recnorm)
			tmp0.BigInt(&norm)
			norm.Div(&norm, _100big)
		}
		if verbose {
			log.Printf("%2d/%2d:  Error: %d   Norm: %d", REP, maxsteps, &best_err, &norm)
		}
		if norm.Cmp(&maxcoeff_big) >= 0 {
			break
		}
	}
	if verbose {
		log.Printf("CANCELLING after step %d/%d.", REP, maxsteps)
		log.Printf("Could not find an integer relation. Norm bound: %d", &norm)
	}
	return nil, errors.New("could not find an integer relation")
}
