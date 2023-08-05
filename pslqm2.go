// A Go port of David Bailey's pslqm2.f90 by Nick Craig-Wood 2022-11-16
//  program tpslqm2
//
//  Revision date:  9 Jan 2022
//
//  AUTHOR:
//   David H. Bailey
//   Lawrence Berkeley National Lab (retired) and University of California, Davis
//   Email: dhbailey@lbl.gov
//
//  COPYRIGHT AND DISCLAIMER:
//   All software in this package (c) 2022 David H. Bailey.
//   By downloading or using this software you agree to the copyright, disclaimer
//   and license agreement in the accompanying file DISCLAIMER.txt.
//
//  DESCRIPTION:
//   This program demonstrates pslqm2, which performs the two-level multipair
//   PSLQ algorithm on an input vector.  A variety of sample input vectors can
//   be generated as inputs to pslqm2, as given in the parameters below.  The
//   pslqm2 routine is suitable for relations up to degree 100 or so; above this
//   level pslqm3 should be used for significantly better performance.
//   For additional details, see:
//
//   David H. Bailey and David J. Broadhurst, "Parallel integer relation
//   detection: Techniques and applications," Mathematics of Computation,
//   vol. 70, no. 236 (Oct 2000), pg. 1719-1736, preprint available at
//   http://www.davidhbailey.com/dhbpapers/ppslq.pdf.
//
//   The pslqm2 routine is 100% THREAD SAFE -- all requisite parameters and
//   arrays are passed through subroutine arguments.

package pslq

// FIXME need to surface errors, especially iterations exceeded

// FIXME need to make sure all the parameters are being set
// - max coefficient isn't being passed

// FIXME shallow copies of big.Float are not supported and may lead to errors.
// Have lots of these - need to fix :-(
//
// Checked once -need to double check

// FIXME make newVector/Matrix auto allocate +1

// FIXME Fortran arays are column major
// They are defined as Matrix(rows, cols) and accessed M(row, col)
// We've defined arrays which are stored as row major in memory
// which may decrease efficiency slightly

/*
   Fortran code on 40 element

   0.48user 0.01system 0:00.50elapsed 100%CPU (0avgtext+0avgdata 8332maxresident)k
   0inputs+0outputs (0major+1168minor)pagefaults 0swaps

   Go code

   0.96user 0.03system 0:00.86elapsed 115%CPU (0avgtext+0avgdata 13124maxresident)k
   0inputs+0outputs (0major+9845minor)pagefaults 0swaps

   Go code with no bounds checks: go install -gcflags=-B ./...

   0.89user 0.02system 0:00.79elapsed 115%CPU (0avgtext+0avgdata 9052maxresident)k
   0inputs+0outputs (0major+6762minor)pagefaults 0swaps

*/

// FIXME the fortran code is nearly twice as fast as the go code :-(

import (
	"fmt"
	"log"
	"math"
	"math/big"
	"os"
	"sort"
)

// anint set res to the nearest integer to x
func anint(x *big.Float, res *big.Float) {
	var prec = x.Prec()
	var half big.Float
	half.SetFloat64(0.5).SetPrec(prec)
	var xr big.Float
	xr.SetPrec(prec)
	if x.Sign() >= 0 {
		xr.Add(x, &half)
	} else {
		xr.Sub(x, &half)
	}
	var z big.Int
	_, _ = xr.Int(&z)
	res.SetPrec(prec)
	res.SetInt(&z)
}

// Returns the base 10 exponent of x
// an integer y so that x ~ 10**y
func exp10(x *big.Float) int {
	exp2 := x.MantExp(nil)
	// FIXME not sure about rounding here
	return int(math.Floor(float64(exp2) * math.Log10(2)))
}

// PSLQM2 parameters set below; all are default (4-byte) integer:
//   idb   Debug level (0-4); default = 2.
//   n     Integer relation vector length; default = 49.
//         For minimal polynomial problems, n = 1 + polynomial degree.
//   ndp   Full precision level in digits; default = 2500.
//         ***Must be <= mpipl in module MPFUNF.
//         ***Must be >= ndpm + precision level required to find relation.
//   ndr   Log10 of the minimum dynamic range in y at detection; default = 30.
//         A detected relation is not deemed reliable unless this is exceeded.
//   nep   Log10 of full precision epsilon for detections; default = 30 - ndp.
//         ***Must not be smaller than the accuracy of input x vector. In other
//         words, if data is accurate to within 10^(-200), then nep > -200.
//   nrb   Log10 of maximum size (Euclidean norm) of acceptable relation;
//         default = 200. Run will be aborted if this is exceeded.
//   nwds  Full precision level in words; default = int (ndp/mpdpw + 2).

func (e *Pslq) RunM2(x []big.Float) ([]big.Int, error) {
	defer func() {
		// If we panic then save the current x into a file
		if e := recover(); e != nil {
			f, err := os.CreateTemp(".", "pslq-crash-dump-*.txt")
			if err != nil {
				log.Print(err)
				panic(e)
			}
			fmt.Printf("Logging crash in %q\n", f.Name())
			defer f.Close()
			fmt.Fprintf(f, "# Crash dump %q\n", f.Name())
			for i := range x {
				fmt.Fprintf(f, "# x[%d]\n%.*f\n", i, digits(x[i].Prec()), &x[i])
			}
			fmt.Fprintf(f, "# End crash dump %q\n", f.Name())

			// Propagate the panic
			panic(e)
		}
	}()
	var (
		idb = 0 // debug level
		n   = len(x)
		ndr = 30
		ndp = int(math.Floor(math.Log10(2) * float64(e.prec)))
		nep = 30 - ndp
		nrb = 200
		r   = newVector(n+1, e.prec)
		iq  int
	)

	if e.verbose {
		idb = 2
	}

	if n <= 1 {
		return nil, ErrorBadArguments
	}

	// At too low precision, the algorithm becomes meaningless
	if e.prec < 64 {
		return nil, ErrorPrecisionTooLow
	}

	// FIXME
	if e.prec < 128 {
		ndr = 12
	}

	if e.verbose && int(e.prec)/max(2, int(n)) < 5 {
		log.Printf("Warning: precision for PSLQ may be too low")
	}

	if e.verbose {
		log.Printf("PSLQ using prec %d and tol %.10g", e.prec, &e.tol)
	}

	if e.tol.Sign() == 0 {
		return nil, ErrorToleranceRoundsToZero
	}

	// Convert to use 1-based indexing to allow us to be
	// consistent with Bailey's indexing.
	xOneIndexed := newVector(n+1, e.prec)
	for i := 1; i <= n; i++ {
		if x[i-1].Sign() == 0 {
			return nil, ErrorZeroArguments
		}
		xOneIndexed[i].Set(&x[i-1])
	}

	pslqm2(idb, n, e.prec, ndr, nrb, nep, xOneIndexed, e.maxsteps, &iq, r)
	if iq == 0 {
		return nil, ErrorNoRelationFound
	}

	// Convert result back to 0 indexed big.Int
	var res = make([]big.Int, n)
	for i := 1; i <= n; i++ {
		_, _ = r[i].Int(&res[i-1])
	}

	return res, nil
}

// The following code performs the two-level, multi-pair PSLQ algorithm.
// David H. Bailey    9 Jan 2022
//
// Arguments are as follows; int means default (4-byte) integer:
//
//	Name  Type    Description
//	idb   int     Debug flag (0-3); increasing idb produces more output.
//	n     int     Length of input vector x and output relation r.
//	prec  int     Full precision level in words. This must be sufficient
//	              to recover the relation.
//	ndr   int     Log10 of the min dynamic range at detection, typically 20
//	              or 30. A detected relation is not deemed reliable unless
//	              this is exceeded.
//	nrb   int     Log10 of max size (Euclidean norm) of acceptable relation,
//	              typically 100 or 200. Run will abort if this is exceeded.
//	nep   int     Log10 of tolerance for full precision relation detection.
//	x     mp_real Input mp_real vector.
//	iq    int     Output flag: 0 (unsuccessful) or 1 (successful).
//	r     mp_real Output integer relation vector, if successful.
//
// The following parameters are set in this routine:
//
//	ipi   int     Iteration print interval when idb >= 2; default = 100.
//	ipm   int     Iteration  check interval for MP iterations; default = 10.
//	itm   int     Maximum iteration count; default = 10^5. Run will be aborted
//	              if this is exceeded.
//	nsq   int     Size of tables used in iterdp and itermpw; default = 8.
//	deps  double  Tolerance for dynamic range check; default = 1d-10.
func pslqm2(idb, n int, prec uint, ndr, nrb, nep int, x []big.Float, itm int, iq *int, r []big.Float) {
	var i, imq, it, its, izd, izm, j, j1, n1, n2, n3 int
	const (
		ipi = 100
		ipm = 10
		//itm = 1000 // FIXME needs to be a parameter
		nsq = 8
	)
	var d3, d4, rn float64
	//real (mprknd) da(n,n), db(n,n), dh(n,n), dsa(n,n), dsb(n,n),   dsh(n,n), dsyq(n,nsq), dy(n), dsy(n), deps
	var da = newMatrixFloat64(n+1, n+1)
	var db = newMatrixFloat64(n+1, n+1)
	var dh = newMatrixFloat64(n+1, n+1)
	var dsa = newMatrixFloat64(n+1, n+1)
	var dsb = newMatrixFloat64(n+1, n+1)
	var dsh = newMatrixFloat64(n+1, n+1)
	var dsyq = newMatrixFloat64(n+1, nsq+1)
	var dy = newVectorFloat64(n + 1)
	var dsy = newVectorFloat64(n + 1)
	const (
		deps = 1e-10
	)
	var depsBig big.Float
	depsBig.SetFloat64(deps).SetPrec(prec)
	//type (mp_real) b(n,n), h(n,n), syq(n,nsq), r(n), x(n),   y(n), dynrange, eps, t1, t2, t3
	var b = newMatrix(n+1, n+1, prec)
	var h = newMatrix(n+1, n+1, prec)
	var syq = newMatrix(n+1, nsq+1, prec)
	var y = newVector(n+1, prec)
	var eps, t1, t2, t3 big.Float
	// external bounddp, dynrange

	// Initialize.
	if idb >= 2 {
		fmt.Printf("PSLQM2 integer relation detection: n = %d\n", n)
	}
	*iq = 0
	it = 0
	its = 0
	imq = 0
	rn = 0.0

	// eps = mpreal(10.0, prec) ** nep // NB nep is -ve
	if nep > 0 {
		panic("nep should be -ve")
	}
	var tBig1, tBig2, tBig3 big.Int
	var one big.Float
	tBig1.SetInt64(10)
	tBig2.SetInt64(int64(-nep))
	tBig3.Exp(&tBig1, &tBig2, nil)
	//fmt.Printf("nep=%d, tBig3=%d\n", nep, &tBig3)
	eps.SetInt(&tBig3).SetPrec(prec)
	one.SetInt64(1).SetPrec(prec)
	eps.Quo(&one, &eps)
	//fmt.Printf("eps=%15.7e\n", &eps)

	if idb >= 2 {
		fmt.Printf("Iteration %7d MP initialization\n", it)
	}
	initmp(idb, n, nsq, prec, b, h, syq, x, y)

L100:

	// Check if dynamic range of y vector is too great for DP iterations
	// (which is often the case at the start of the run). If so, then
	// perform MP iterations instead of DP iterations.
	t1 = dynrange(n, prec, y)
	if idb >= 2 {
		fmt.Printf("Iteration %7d Min/max ratio in y = %15.7e\n", it, &t1)
	}
	if t1.Cmp(&depsBig) < 0 {
		goto L120
	}

	// DP iterations.
	// Initialize DP arrays from MP arrays.
	if idb >= 3 {
		fmt.Printf("Iteration %7d Start DP iterations\n", it)
	}

	initdp(idb, n, nsq, da, db, dh, dy, dsyq, h, y)

	// Save DP arrays and iteration count.
	savedp(n, da, db, dh, dy, dsa, dsb, dsh, dsy)
	its = it

	// Perform an LQ decomposition on dh.
	lqdp(n, n-1, dh)

L110:
	// Perform one DP iteration.
	it = it + 1
	if idb >= 3 || (idb >= 2 && it%ipi == 0) {
		fmt.Printf("Iteration %7d\n", it)
	}
	// imq = 1 // FIXME do non multiple iterations
	iterdp(idb, it, n, nsq, da, db, dh, dsyq, dy, &imq, &izd)

	// Test conditions on iterdp output flag izd:
	// 0: Iteration was uneventful; periodically save arrays and continue.
	// 1: Relation found or DP precision exhausted; perform MP update.
	// 2: Very large value appeared in DA or DB; abort DP iter and do MP update.
	if izd == 0 {
		if (it-its)%ipm == 0 {
			savedp(n, da, db, dh, dy, dsa, dsb, dsh, dsy)
			its = it
		}
		goto L110
	} else {
		// Check if DP iteration was aborted above; if so, then revert to previous data.
		if izd == 2 {
			// DP iteration was aborted -- revert to previous data.
			it = its
			savedp(n, dsa, dsb, dsh, dsy, da, db, dh, dy)
		}

		// Update the MP arrays from the DP arrays.
		if idb >= 2 {
			fmt.Printf("Iteration %7d MP update\n", it)
		}
		updtmp(idb, it, n, prec, da, db, eps, b, h, y, &izm)

		// Compute norm bound using DP.
		for j = 1; j <= n-1; j++ {
			for i = 1; i <= n; i++ {
				dh[i][j], _ = h[i][j].Float64()
			}
		}

		lqdp(n, n-1, dh)
		d3 = bounddp(n, dh)
		if d3 == -1.0 {
			goto L150
		}
		rn = max(rn, d3)
		if idb >= 2 {
			fmt.Printf("Iteration %7d  Norm bound = %15.7e Max. bound = %15.7e\n", it, d3, rn)
		}

		// Check if iteration limit or norm bound limit is exceeded; if so, quit.
		if it > itm {
			if idb >= 1 {
				fmt.Printf("Iteration limit exceeded %d\n", itm)
			}
			goto L150
		}
		if math.Log10(rn) > float64(nrb) {
			if idb >= 1 {
				fmt.Printf("Norm bound limit exceeded. %d\n", nrb)
			}
			goto L150
		}

		// Test conditions on updtmp output flag izm:
		// 0: MP update was uneventful; goto tag 100 above, provided izd = 0;
		//      if izd = 2, then go perform MP iterations.
		// 1: A small value was found in MP update; go output relation.
		// 2: Precision is exhausted; quit.
		if izm == 0 {
			if izd == 2 {
				goto L120
			} else {
				goto L100
			}
		} else if izm == 1 {
			goto L140
		} else if izm == 2 {
			goto L150
		}
	}

L120:

	// MP iterations.
	if idb >= 2 {
		fmt.Printf("Iteration %7d Start MP iterations\n", it)
	}

	// Perform an LQ decomposition on h.
	lqmp(n, n-1, prec, h)

L130:

	// Perform one MP iteration.
	it = it + 1
	if idb >= 2 {
		fmt.Printf("Iteration %7d", it)
	}
	// imq = 1 // FIXME do non multiple iterations
	itermp(idb, it, n, nsq, prec, eps, b, h, syq, y, &imq, &izm)

	// Test conditions on itermp output flag izm:
	// 0: Iteration was uneventful; periodically check for DP; continue.
	// 1: A small value was found; go output relation.
	// 2: Precision is exhausted; quit.
	if izm == 0 {

		// Periodically check to see if DP iterations can be resumed (but perform
		// at least IPM iterations in MP).
		if (it-its)%ipm == 0 {
			t1 = dynrange(n, prec, y)
			if idb >= 2 {
				fmt.Printf("Iteration %7d Min/max ratio in y = %15.7e\n", it, &t1)
			}
			if t1.Cmp(&depsBig) > 0 {
				if idb >= 2 {
					fmt.Printf("Iteration %7d Return to DP iterations\n", it)
				}
				goto L100
			}
		}
		goto L130
	} else if izm == 1 {
		goto L140
	} else if izm == 2 {
		goto L150
	}

	// A relation has been detected.  Output the final norm bound and other info.
L140:

	t1.SetFloat64(1.0e300).SetPrec(prec)
	t2.SetFloat64(0.0).SetPrec(prec)
	t3.SetFloat64(0.0).SetPrec(prec)

	// Select the relation corresponding to the smallest y entry and compute norm.
	for j = 1; j <= n; j++ {
		var absY big.Float
		absY.Abs(&y[j])
		if absY.Cmp(&t1) < 0 {
			j1 = j
			t1.Set(&absY)
		}
		t2.Set(maxBig(&t2, &absY))
	}

	for i = 1; i <= n; i++ {
		r[i].Set(&b[j1][i])
		// t3 = t3 + r[i]*r[i]
		var t0 big.Float
		t0.Mul(&r[i], &r[i])
		t3.Add(&t3, &t0)
	}

	// t3 = sqrt(t3)
	t3.Sqrt(&t3)
	// mpdecmd(t3, d3, n3)
	n3 = exp10(&t3)
	// d3 = d3 * 10.0**n3
	d3, _ = t3.Float64()

	// Output the final norm bound and other info.
	if idb >= 1 {
		for j = 1; j <= n-1; j++ {
			for i = 1; i <= n; i++ {
				dh[i][j], _ = h[i][j].Float64()
			}
		}
		lqdp(n, n-1, dh)
		d4 = bounddp(n, dh)
		rn = max(rn, d4)
		fmt.Printf("Iteration %7d Relation detected\nMin, max of y = %15.7e, %15.7e\nMax. bound = %15.7e\n", it, &t1, &t2, rn)
		fmt.Printf("Index of relation = %d Norm = %15.7e Residual = %15.7e\n", j1, d3, &t1)
		printVector("r", r)
	}
	// mpdecmd(t1, d1, n1)
	n1 = exp10(&t1)
	// mpdecmd(t2, d2, n2)
	n2 = exp10(&t2)

	// If run was successful, set iq = 1.
	if t1.Sign() == 0 {
		n1 = nep
	}
	// fmt.Printf("n3=%d, nrb=%d, n3 <= nrb is %v\n", n3, nrb, n3 <= nrb)
	// fmt.Printf("n2=%d, n1=%d, n2-d1=%d, ndr=%d, n2-n1 >= ndr is %v\n", n2, n1, n2-n1, ndr, n2-n1 >= ndr)
	if n3 <= nrb && n2-n1 >= ndr {
		*iq = 1
	} else {
		if idb >= 2 {
			fmt.Printf("Relation is too large or insufficient dynamic range.\n")
		}
	}

L150:

	return
}

//------------------------------

// First-level subroutines.

// This returns the dynamic range of y, i.e., ratio of min|y_k| to max|y_k|,
// using MP precision.
func dynrange(n int, prec uint, y []big.Float) big.Float {
	var i int
	var t1, t2, t3 big.Float
	//type (mp_real) y(n)

	t1.SetFloat64(1.0e300).SetPrec(prec)
	t2.SetFloat64(0.0).SetPrec(prec)

	// Find the min and max absolute value in the y vector.
	for i = 1; i <= n; i++ {
		t3.Abs(&y[i])
		t1.Set(minBig(&t1, &t3))
		t2.Set(maxBig(&t2, &t3))
	}

	//return t1 / t2
	t1.Quo(&t1, &t2)
	return t1
}

// This initializes the DP arrays from the MP arrays.
// This is performed in four-word precision.
// Input:  idb, n, nsq, h, y.
// Output: da, db, dh, dy, dsyq.
func initdp(idb, n, nsq int, da, db, dh [][]float64, dy []float64, dsyq [][]float64, h [][]big.Float, y []big.Float) {
	var i, j int
	//real (mprknd) da(n,n), db(n,n), dh(n,n), dy(n), dsyq(n,nsq)
	// type (mp_real) h(n,n), y(n), t1, t2
	var t1, t2, t3 big.Float
	const nwx = 4
	const reducedPrec = 128

	// Find the max absolute value in the y vector.
	t2.SetFloat64(0.0).SetPrec(reducedPrec)
	for i = 1; i <= n; i++ {
		// t2 = max (t2, mpreal (abs (h(j,j)), nwx))
		t3.Set(&y[i])
		t3.SetPrec(reducedPrec)
		t3.Abs(&t3)
		t2.Set(maxBig(&t2, &t3))
	}

	// Set dy to be the scaled y vector.
	//t1 = 1.0 / t2
	t1.SetFloat64(1.0).SetPrec(reducedPrec)
	t1.Quo(&t1, &t2)
	for i = 1; i <= n; i++ {
		// dh(i,j) = t1 * mpreal (h(i,j), nwx)
		t3.Set(&y[i])
		t3.SetPrec(reducedPrec)
		t3.Mul(&t1, &t3)
		dy[i], _ = t3.Float64()
	}

	// Find the maximum absolute value of the h matrix diagonals.
	t2.SetFloat64(0.0).SetPrec(reducedPrec)
	for j = 1; j <= n-1; j++ {
		//   t2 = max (t2, mpreal (abs (h(j,j)), nwx))
		t3.Set(&h[j][j])
		t3.SetPrec(reducedPrec)
		t3.Abs(&t3)
		t2.Set(maxBig(&t2, &t3))
	}

	// Set dh to be the scaled h matrix.
	// t1 = 1.0 / t2
	t1.SetFloat64(1.0).SetPrec(reducedPrec)
	t1.Quo(&t1, &t2)
	for j = 1; j <= n-1; j++ {
		for i = 1; i <= n; i++ {
			// dh[i][j] = t1 * mpreal(h[i][j], nwx)
			t3.Set(&h[i][j])
			t3.SetPrec(reducedPrec)
			t3.Mul(&t1, &t3)
			dh[i][j], _ = t3.Float64()
		}
	}

	// Set da and db to the identity.
	for j = 1; j <= n; j++ {
		for i = 1; i <= n; i++ {
			da[i][j] = 0.0
			db[i][j] = 0.0
		}
		da[j][j] = 1.0
		db[j][j] = 1.0
	}

	// Zero the dsyq array.
	for j = 1; j <= nsq; j++ {
		for i = 1; i <= n; i++ {
			dsyq[i][j] = 0.0
		}
	}

	if idb >= 3 {
		fmt.Printf("initdp: Scaled dy vector:\n")
		printVectorFloat64("dy", dy)
		fmt.Printf("initdp: Scaled dh matrix:\n")
		printMatrixFloat64("dh", dh)
	}

	return
}

// This initializes MP arrays at the beginning.
// This is performed in full precision.
// Input: idb, n, nsq, prec.
// Output: b, h, syq, x, y.
func initmp(idb, n, nsq int, prec uint, b, h, syq [][]big.Float, x, y []big.Float) {
	var i, j int
	//integer ix(n)
	//var ix = newVectorInt(n+1)
	//real (mprknd) dx(n)
	//var dx = newVectorFloat64(n+1)
	//type (mp_real) b(n,n), h(n,n), s(n), syq(n,nsq), x(n), y(n), t1
	var t1, t2 big.Float
	var s = newVector(n+1, prec)

	if idb >= 3 {
		fmt.Printf("initmp: Input x vector:\n")
		printVector("x", x)
	}

	// Set b to the identity matrix.
	for j = 1; j <= n; j++ {
		for i = 1; i <= n; i++ {
			b[i][j].SetInt64(0)
		}
		b[j][j].SetInt64(1)
	}

	t1.SetInt64(0).SetPrec(prec)

	// Compute the x vector, the square root of the partial sum of squares of x,
	// and the y vector, which is the normalized x vector.
	for i = n; i >= 1; i-- {
		// t1 = t1 + x[i]*x[i]
		t2.Mul(&x[i], &x[i])
		t1.Add(&t1, &t2)
		s[i].Sqrt(&t1)
	}

	// t1 = 1.0 / s[1]
	t1.SetFloat64(1.0).SetPrec(prec)
	t1.Quo(&t1, &s[1])

	for i = 1; i <= n; i++ {
		y[i].Mul(&t1, &x[i])
		s[i].Mul(&t1, &s[i])
	}

	// Compute the initial h matrix.
	for j = 1; j <= n-1; j++ {
		for i = 1; i <= j-1; i++ {
			h[i][j].SetInt64(0)
		}

		h[j][j].Quo(&s[j+1], &s[j])
		// t1 = y[j] / (s[j] * s[j+1])
		t1.Mul(&s[j], &s[j+1])
		t1.Quo(&y[j], &t1)

		for i = j + 1; i <= n; i++ {
			// h[i][j] = -y[i] * t1
			h[i][j].Mul(&y[i], &t1)
			h[i][j].Neg(&h[i][j])
		}
	}

	// Zero the syq array.
	for j = 1; j <= nsq; j++ {
		for i = 1; i <= n; i++ {
			syq[i][j].SetInt64(0)
		}
	}

	if idb >= 3 {
		fmt.Printf("initmp: Initial y vector:\n")
		printVector("y", y)
		fmt.Printf("initmp: Initial h matrix:\n")
		printMatrix("h", h)
	}

	return
}

// This performs one iteration of the PSLQ algorithm using DP arithmetic.
// Input: idb, it, n, nsq, da, db, dh, dsyq, dy.
// Output: da, db, dh, dsyq, dy, imq, izd.
//
// NOTE: Parameter tmx2 = 2^52, not 2^53, so as to ensure that values > 2^53
// never arise, even as intermediate values, in the update loop below.
func iterdp(idb, it, n, nsq int, da, db, dh, dsyq [][]float64, dy []float64, imq, izd *int) {
	var i, ii, ij, im, im1, j, j1, j2, k, mpr, mq int
	var gam, t1, t2, t3, t4 float64
	// integer ip(n), ir(n), is(n)
	var ip = newVectorInt(n + 1)
	var ir = newVectorInt(n + 1)
	var is = newVectorInt(n + 1)
	// real (mprknd) da(n,n), db(n,n), dh(n,n), dq(n), dsyq(n,nsq),   dt(n,n), dy(n)
	var dq = newVectorFloat64(n + 1)
	var dt = newMatrixFloat64(n+1, n+1)
	const (
		tmx1 = 1.0e13
		tmx2 = float64(1 << 52)
		deps = 1.0e-14
	)

	*izd = 0
	mpr = int(math.Round(0.4e0 * float64(n)))
	gam = math.Sqrt(4.0 / 3.0)
	gamPow := gam

	// Compute dq vector = {gam^i * |dh[i][i]|}, then sort in ascending order.
	for i = 1; i <= n-1; i++ {
		dq[i] = gamPow * math.Abs(dh[i][i])
		gamPow *= gam
	}

	// This routine sorts the entries of the N-long DP vector A into ascending
	// order using the quicksort algorithm.  The permutation vector that would
	// sort the vector is returned in IP.
	// Input: n, a.
	// Output: ip.
	//qsortdp(n-1, dq, ip)
	for i = 1; i <= n; i++ {
		ip[i] = i
	}
	sort.Slice(ip[1:len(ip)-1], func(i, j int) bool {
		// fmt.Printf("i=%d, j=%d, ip[j]=%d, ip[j]=%d, dq[ip[i]]=%15.7e, dq[ip[j]]=%15.7e, less: %v\n", i+1, j+1, ip[i], ip[j], dq[ip[i]], dq[ip[j]], dq[ip[i]] < dq[ip[j]])
		return dq[ip[i+1]] < dq[ip[j+1]]
	})
	// printVectorFloat64("dq", dq)
	// printVectorInt("ip", ip)
	for i = 1; i <= n-2; i++ {
		if dq[ip[i]] > dq[ip[i+1]] {
			fmt.Printf("%d\n", ip)
			fmt.Printf("%.5g\n", dq)
			fmt.Printf("i = %d, ip[i]=%d, ip[i+1]=%d\n", i, ip[i], ip[i+1])
			panic("list not sorted")
		}
	}
	// // FIXME
	// for i = 1; i <= n; i++ {
	// 	ip[i] += 1
	// }

	// Select up to mpr disjoint pairs of indices (m,m+1), where m is an index
	// from the list of the largest dq(i).
	for i = 1; i <= n; i++ {
		is[i] = 0
	}

	if *imq == 0 {
		mq = mpr
	} else {
		mq = 1
		*imq = 0
	}
	ii = n

	for i = 1; i <= mq; i++ {
		ii = ii - 1
		if ii == 0 {
			mq = i - 1
			break
		}
		j1 = ip[ii]
		j2 = j1 + 1
		if j2 > n {
			fmt.Printf("ii=%d, j1=%d, j2=%d, n=%d\n", ii, j1, j2, n)
			fmt.Printf("dq = %.5g\n", dq)
			fmt.Printf("ip = %d\n", ip)
		}
		if /*j2 > n ||*/ is[j1] != 0 || is[j2] != 0 {
			i-- // do this iteration again
			continue
		}
		ir[i] = j1
		is[j1] = 1
		is[j2] = 1
	}
	//fmt.Printf("mq = %d, ir = %d\n", mq, ir[1:mq+1])

	// Exchange the pairs of entries of dy, and rows of da, db and dh.
	//printVectorFloat64("before dy", dy)
	for j = 1; j <= mq; j++ {
		im = ir[j]
		if im == 0 {
			panic("out of bounds read with im == 0")
		}
		im1 = im + 1
		dy[im], dy[im1] = dy[im1], dy[im]

		for i = 1; i <= n; i++ {
			da[im][i], da[im1][i] = da[im1][i], da[im][i]
			db[im][i], db[im1][i] = db[im1][i], db[im][i]
		}

		for i = 1; i <= n-1; i++ {
			dh[im][i], dh[im1][i] = dh[im1][i], dh[im][i]
		}
	}
	//printVectorFloat64("after dy", dy)

	// Eliminate the "corners" produced by the above permutation in dh.
	for j = 1; j <= mq; j++ {
		im = ir[j]
		if im == 0 {
			panic("out of bounds read with im == 0")
		}
		im1 = im + 1
		if im <= n-2 {
			t1 = dh[im][im]
			t2 = dh[im][im1]
			t3 = math.Sqrt(t1*t1 + t2*t2)
			if t3 == 0 {
				panic(fmt.Sprintf("t3 is zero t1=%15.7e, t2=%15.7e", t1, t2))
			}
			t1 = t1 / t3
			t2 = t2 / t3

			for i = im; i <= n; i++ {
				t3 = dh[i][im]
				t4 = dh[i][im1]
				dh[i][im] = t1*t3 + t2*t4
				dh[i][im1] = -t2*t3 + t1*t4
			}
		}
	}

	// Perform reduction on dh, using the diagonal scheme.  Multipliers are
	// saved in the dt array.
	for i = 2; i <= n; i++ {
		for j = 1; j <= n-i+1; j++ {
			ij = i + j - 1

			for k = j + 1; k <= ij-1; k++ {
				dh[ij][j] = dh[ij][j] - dt[ij][k]*dh[k][j]
			}

			if dh[j][j] == 0 {
				panic(fmt.Sprintf("dh[%d][%d] is 0", j, j))
			}
			dt[ij][j] = math.Round(dh[ij][j] / dh[j][j])
			dh[ij][j] = dh[ij][j] - dt[ij][j]*dh[j][j]
		}
	}

	// Update dy, using the dt array.  Find min absolute value of dy.
	//fmt.Printf("dy before = %.3e\n", dy)
	//printMatrixFloat64("dt", dt)
	t1 = math.Abs(dy[n])
	for j = 1; j <= n-1; j++ {
		for i = j + 1; i <= n; i++ {
			// if dt[i][j] != 0 {
			// 	fmt.Printf("nonzero dt[%2d][%2d] = %15.7e\n", i, j, dt[i][j])
			// }
			dy[j] = dy[j] + dt[i][j]*dy[i]
		}

		t1 = min(t1, math.Abs(dy[j]))
	}
	//fmt.Printf("dy after = %.3e\n", dy)

	// Update da and db, using the dt array.  Find the max absolute value of
	// da and db entries as they are calculated (not merely at the end).
	t2 = 0.0
	for k = 1; k <= n; k++ {
		dq[k] = 0.0

		for j = 1; j <= n-1; j++ {
			for i = j + 1; i <= n; i++ {
				da[i][k] = da[i][k] - dt[i][j]*da[j][k]
				db[j][k] = db[j][k] + dt[i][j]*db[i][k]
				dq[k] = max(dq[k], max(math.Abs(da[i][k]), math.Abs(db[j][k])))
			}
		}
	}

	for k = 1; k <= n; k++ {
		t2 = max(t2, dq[k])
	}

	//fmt.Printf("t1=%15.7e, t2=%15.7e, deps=%15.7e, tmx1=%15.7e, tmx2=%15.7e\n", t1, t2, deps, tmx1, tmx2)
	if t1 <= deps {
		if idb >= 2 { // Was 3
			fmt.Printf("Iteration %7d iterdp: Small value in dy = %15.7e\n", it, t1)
		}
		*izd = 1
	}

	if t2 > tmx1 && t2 <= tmx2 {
		if idb >= 2 { // Was 3
			fmt.Printf("Iteration %7d iterdp: Large value in da or db = %15.7e\n", it, t2)
		}
		*izd = 1
	} else if t2 > tmx2 {
		if idb >= 2 {
			fmt.Printf("Iteration %7d iterdp: Very large value in da or db = %15.7e\n", it, t2)
		}
		*izd = 2
		return
	}

	//fmt.Printf("dy = %.5g\n", dy[1:])
	// Compare the dy vector with those of recent iterations.  If a duplicate is
	// found, then the next iteration must be performed with mq = 1.
	for j = 1; j <= nsq; j++ {
		t1 = 0.0

		for i = 1; i <= n; i++ {
			//fmt.Printf("y[%2d] = %11.5g, dsyq[%2d][%2d] = %11.5g, delta = %11.5g\n", i, dy[i], i, j, dsyq[i][j], math.Abs(dy[i]-dsyq[i][j]))
			t1 = max(t1, math.Abs(dy[i]-dsyq[i][j]))
		}

		if t1 <= deps {
			if idb >= 2 {
				fmt.Printf("Iteration %7d iterdp: Duplicate found, j = %d (t1=%.5g <= deps=%.5g)\n", it, j, t1, deps)
			}
			*imq = 1
			break
		}
	}

	// Place the vector dy in the table dsyq.
	k = 1 + (it % nsq)

	for i = 1; i <= n; i++ {
		dsyq[i][k] = dy[i]
	}

	if idb >= 3 {
		// This is different!!!! FIXME
		fmt.Printf("iterdp: Updated dy:\n")
		printVectorFloat64("dy", dy)
		fmt.Printf("iterdp: Updated da matrix:\n")
		printMatrixFloat64("da", da)
		fmt.Printf("iterdp: Updated db matrix:\n")
		printMatrixFloat64("db", db)
		fmt.Printf("iterdp: Updated dh matrix:\n")
		printMatrixFloat64("dh", dh)
	}

	return
}

// This performs one iteration of the PSLQM algorithm using MP arithmetic.
// This is performed in medium precision.
// Input: idb, it, n, nsq, prec, eps, imq.
// Output: b, h, syq, y, imq, izm.
func itermp(idb, it, n, nsq int, prec uint, eps big.Float, b, h, syq [][]big.Float, y []big.Float, imq *int, izm *int) {
	var i, ii, ij, im, im1, j, j1, j2, k, mpr, mq int
	const ntl = 72
	var gam, gamPow, t1, t2, t3, t4, t5, t6, teps big.Float
	var ip = newVectorInt(n + 1)
	var ir = newVectorInt(n + 1)
	var is = newVectorInt(n + 1)
	//integer ip(n), ir(n), is(n)
	//real (mprknd) dx(n)
	//var dx = newVectorFloat64(n+1)
	//type (mp_real) b(n,n), h(n,n), q(n), syq(n,nsq), t(n,n), y(n)
	var q = newVector(n+1, prec)
	var t = newMatrix(n+1, n+1, prec)

	//teps = 2.0e0 ** ntl * eps
	teps.SetFloat64(math.Pow(2, ntl))
	teps.Mul(&teps, &eps)
	*izm = 0
	mpr = int(math.Round(0.4e0 * float64(n)))
	//gam = sqrt(mpreald(4.0e0, prec) / mpreald(3.0e0, prec))
	t1.SetInt64(4).SetPrec(prec)
	t2.SetInt64(3).SetPrec(prec)
	t1.Quo(&t1, &t2)
	gam.Sqrt(&t1)
	gamPow.Set(&gam)

	// Compute q vector = {gam^i * |h[i][i]|}, then sort in ascending order.
	for i = 1; i <= n-1; i++ {
		//q(i) = gam**i * abs(h[i][i])
		t1.Abs(&h[i][i])
		q[i].Mul(&gamPow, &t1)
		gamPow.Mul(&gamPow, &gam)
	}

	// This routine sorts the entries of the N-long MP vector A into ascending
	// order using the quicksort algorithm.  The permutation vector that would
	// sort the vector is returned in IP.
	// Input: n, a, ip.
	// Output: ip.
	// qsortmp(n-1, q, ip)
	for i = 1; i <= n; i++ {
		ip[i] = i
	}
	sort.Slice(ip[1:len(ip)-1], func(i, j int) bool {
		return q[ip[i+1]].Cmp(&q[ip[j+1]]) < 0
	})
	// FIXME
	// for i = 1; i <= n; i++ {
	// 	ip[i] += 1
	// }

	// Select up to mpr disjoint pairs of indices (m,m+1), where m is an index
	// from the list of the largest q(i).
	for i = 1; i <= n; i++ {
		is[i] = 0
	}

	if *imq == 0 {
		mq = mpr
	} else {
		mq = 1
		*imq = 0
	}
	ii = n

	for i = 1; i <= mq; i++ {
		ii = ii - 1
		if ii == 0 {
			mq = i - 1
			break
		}
		j1 = ip[ii]
		j2 = j1 + 1
		if /*j2 > n ||*/ is[j1] != 0 || is[j2] != 0 {
			i-- // do this iteration again
			continue
		}
		ir[i] = j1
		is[j1] = 1
		is[j2] = 1
	}

	// Exchange the pairs of entries of y, and rows of b and h.
	for j = 1; j <= mq; j++ {
		im = ir[j]
		if im == 0 {
			panic("out of bounds read with im == 0")
		}
		im1 = im + 1
		y[im], y[im1] = y[im1], y[im]

		for i = 1; i <= n; i++ {
			b[im][i], b[im1][i] = b[im1][i], b[im][i]
		}

		for i = 1; i <= n-1; i++ {
			h[im][i], h[im1][i] = h[im1][i], h[im][i]
		}
	}

	// Eliminate the "corners" produced by the above permutation in h.
	for j = 1; j <= mq; j++ {
		im = ir[j]
		if im == 0 {
			panic("out of bounds read with im == 0")
		}
		im1 = im + 1
		if im <= n-2 {
			t1.Set(&h[im][im])
			t2.Set(&h[im][im1])
			// t3 = sqrt(t1**2 + t2**2)
			t3.Mul(&t1, &t1)
			t4.Mul(&t2, &t2)
			t4.Add(&t4, &t3)
			t3.Sqrt(&t4)
			if t3.Sign() == 0 {
				panic(fmt.Sprintf("t3 is zero t1=%15.7e, t2=%15.7e", &t1, &t2))
			}
			//t1 = t1 / t3
			//t2 = t2 / t3
			t1.Quo(&t1, &t3)
			t2.Quo(&t2, &t3)

			for i = im; i <= n; i++ {
				t3.Set(&h[i][im])
				t4.Set(&h[i][im1])
				// h[i][im] = t1*t3 + t2*t4
				t5.Mul(&t1, &t3)
				t6.Mul(&t2, &t4)
				h[i][im].Add(&t5, &t6)
				// h[i][im1] = -t2*t3 + t1*t4
				t5.Mul(&t2, &t3)
				t6.Mul(&t1, &t4)
				h[i][im1].Sub(&t6, &t5)
			}
		}
	}

	// Perform reduction on h, using the diagonal scheme.  Multipliers are
	// saved in the t array.
	for i = 2; i <= n; i++ {
		for j = 1; j <= n-i+1; j++ {
			ij = i + j - 1

			for k = j + 1; k <= ij-1; k++ {
				// h[ij][j] = h[ij][j] - t[ij][k]*h[k][j]
				t1.Mul(&t[ij][k], &h[k][j])
				h[ij][j].Sub(&h[ij][j], &t1)
			}

			if h[j][j].Sign() == 0 {
				panic(fmt.Sprintf("h[%d][%d] is 0", j, j))
			}
			//t[ij][j] = anint[h[ij][j]/h[j][j]]
			t1.Quo(&h[ij][j], &h[j][j])
			anint(&t1, &t[ij][j])
			//h[ij][j] = h[ij][j] - t[ij][j]*h[j][j]
			t1.Mul(&t[ij][j], &h[j][j])
			h[ij][j].Sub(&h[ij][j], &t1)
		}
	}

	// Update y, using the t array.  Find min absolute value of y.
	t1.Abs(&y[n])
	j1 = n

	for j = 1; j <= n-1; j++ {
		for i = j + 1; i <= n; i++ {
			//y[j] = y[j] + t[i][j]*y[i]
			t3.Mul(&t[i][j], &y[i])
			y[j].Add(&y[j], &t3)
		}

		t3.Abs(&y[j])
		//if abs(y[j]) < t1 {
		if t3.Cmp(&t1) < 0 {
			j1 = j
			t1.Set(&t3)
		}
	}

	// Update b, using the t array.
	for k = 1; k <= n; k++ {
		for j = 1; j <= n-1; j++ {
			for i = j + 1; i <= n; i++ {
				//b[j][k] = b[j][k] + t[i][j]*b[i][k]
				t3.Mul(&t[i][j], &b[i][k])
				b[j][k].Add(&b[j][k], &t3)
			}
		}
	}

	//  Find the largest entry of b in the same row as the smallest y.
	t2.SetFloat64(0.0).SetPrec(prec)

	for i = 1; i <= n; i++ {
		t3.Abs(&b[j1][i])
		//t2 = max(t2, abs(b[j1][i]))
		t2.Set(maxBig(&t2, &t3))
	}

	t3.Mul(&t2, &teps)
	// if t1 <= t2*teps {
	if t1.Cmp(&t3) <= 0 {
		if idb >= 2 {
			fmt.Printf("Iteration %7d itermp: Small value in y = %15.7e\n", it, &t1)
		}
		t3.Mul(&t2, &eps)
		//if t1 <= t2*eps {
		if t1.Cmp(&t3) <= 0 {
			*izm = 1
		} else {
			if idb >= 1 {
				fmt.Printf("Iteration %7d itermp: Precision exhausted\n", it)
			}
			*izm = 2
		}
	}

	// Compare the y vector with those of recent iterations.  If a duplicate is
	// found, then the next iteration must be performed with mq = 1.
	for j = 1; j <= nsq; j++ {
		t1.SetFloat64(0.0).SetPrec(prec)

		for i = 1; i <= n; i++ {
			//t1 = max(t1, abs(y[i]-syq[i][j]))
			t3.Sub(&y[i], &syq[i][j])
			t3.Abs(&t3)
			t1.Set(maxBig(&t1, &t3))
		}

		t3.Mul(&t2, &teps)
		//if t1 <= t2*teps {
		if t1.Cmp(&t3) <= 0 {
			if idb >= 2 {
				fmt.Printf("Iteration %7d itermp: Duplicate found, j = %d\n", it, j)
			}
			*imq = 1
			break
		}
	}

	// Place the vector y in the table syq.
	k = 1 + (it % nsq)

	for i = 1; i <= n; i++ {
		syq[i][k].Set(&y[i])
	}

	if idb >= 3 {
		fmt.Printf("itermp: Updated y:\n")
		printVector("y", y)
		fmt.Printf("itermp: Updated b matrix:\n")
		printMatrix("b", b)
		fmt.Printf("itermp: Updated h matrix:\n")
		printMatrix("h", h)
	}

	return
}

// This performs an LQ decomposition on the DP matrix dh.  It is a simplified
// and transposed adaptation of the subroutine dqrdc from Linpack.
// Input: n, m, dh.
// Output: dh.
func lqdp(n, m int, dh [][]float64) {
	var i, j, l, lup, ml int
	//real (mprknd) dh(n,m), nrmxl, one, t, zero
	var nrmxl, one, t, zero float64

	zero = 0.0
	one = 1.0
	lup = min(m, n)

	// Perform the householder reduction of dh.
	for l = 1; l <= lup; l++ {
		if l == m {
			continue
		}

		// Compute the householder transformation for column l.
		ml = m - l
		t = zero

		for i = 0; i <= ml; i++ {
			t = t + dh[l][l+i]*dh[l][l+i]
		}

		nrmxl = math.Sqrt(t)
		if nrmxl == zero {
			continue
		}
		if dh[l][l] != zero {
			// nrmxl = sign(nrmxl, dh[l][l])
			if dh[l][l] < 0 {
				nrmxl = -nrmxl
			}
		}
		t = one / nrmxl

		for i = 0; i <= ml; i++ {
			dh[l][l+i] = t * dh[l][l+i]
		}

		dh[l][l] = one + dh[l][l]

		// Apply the transformation to the remaining columns, updating the norms.
		for j = l + 1; j <= n; j++ {
			t = zero

			for i = 0; i <= ml; i++ {
				t = t + dh[l][l+i]*dh[j][l+i]
			}

			t = -t / dh[l][l]

			for i = 0; i <= ml; i++ {
				dh[j][l+i] = dh[j][l+i] + t*dh[l][l+i]
			}
		}

		// Save the transformation.
		dh[l][l] = -nrmxl
	}

	// Zero dh above the diagonal.
	for j = 1; j <= m; j++ {
		for i = 1; i <= j-1; i++ {
			dh[i][j] = 0.0
		}
	}

	return
}

// This performs an LQ decomposition on the MP matrix h.  It is a simplified
// and transposed adaptation of the subroutine dqrdc from Linpack.
// Input: n, m, prec, h.
// Output: h.
func lqmp(n, m int, prec uint, h [][]big.Float) {
	var i, j, l, lup, ml int
	var nrmxl, one, t, zero, t1 big.Float

	zero.SetFloat64(0.0).SetPrec(prec)
	one.SetFloat64(1.0).SetPrec(prec)
	lup = min(m, n)

	// Perform the householder reduction of h.
	for l = 1; l <= lup; l++ {
		if l == m {
			continue
		}

		// Compute the householder transformation for column l.
		ml = m - l
		t.Set(&zero)

		for i = 0; i <= ml; i++ {
			//t = t + h[l][l+i]*h[l][l+i]
			t1.Mul(&h[l][l+i], &h[l][l+i])
			t.Add(&t, &t1)
		}

		nrmxl.Sqrt(&t)
		if nrmxl.Sign() == 0 {
			continue
		}
		//if h[l][l].Sign() != 0 {
		// nrmxl = sign(nrmxl, h[l][l])
		if h[l][l].Sign() < 0 {
			nrmxl.Neg(&nrmxl)
		}
		//t = one / nrmxl
		t.Quo(&one, &nrmxl)

		for i = 0; i <= ml; i++ {
			// h[l][l+i] = t * h[l][l+i]
			h[l][l+i].Mul(&t, &h[l][l+i])
		}

		h[l][l].Add(&one, &h[l][l])

		// Apply the transformation to the remaining columns, updating the norms.
		for j = l + 1; j <= n; j++ {
			t.Set(&zero)

			for i = 0; i <= ml; i++ {
				// t = t + h[l][l+i]*h[j][l+i]
				t1.Mul(&h[l][l+i], &h[j][l+i])
				t.Add(&t, &t1)
			}

			// t = -t / h[l][l]
			t.Quo(&t, &h[l][l])
			t.Neg(&t)

			for i = 0; i <= ml; i++ {
				// h[j][l+i] = h[j][l+i] + t*h[l][l+i]
				t1.Mul(&t, &h[l][l+i])
				h[j][l+i].Add(&h[j][l+i], &t1)

			}
		}

		// Save the transformation.
		// dh(l,l) = - nrmxl
		h[l][l].Neg(&nrmxl)
	}

	// Zero h above the diagonal.
	for j = 1; j <= m; j++ {
		for i = 1; i <= j-1; i++ {
			h[i][j].SetInt64(0)
		}
	}

	return
}

// This saves the arrays dy, da, db, dh in case dp iterations must be aborted.
// A call to the same routine, with (da,db,dh,dy) and (dsa,dsb,dsh,dsy)
// exchanged, serves to restore these arrays.
func savedp(n int, da, db, dh [][]float64, dy []float64, dsa, dsb, dsh [][]float64, dsy []float64) {
	var i, j int
	// real (mprknd) da(n,n), db(n,n), dh(n,n), dy(n), dsa(n,n), dsb(n,n),   dsh(n,n), dsy(n)

	for i = 1; i <= n; i++ {
		dsy[i] = dy[i]
	}

	for j = 1; j <= n; j++ {
		for i = 1; i <= n; i++ {
			dsa[i][j] = da[i][j]
			dsb[i][j] = db[i][j]
		}
	}

	for j = 1; j <= n-1; j++ {
		for i = 1; i <= n; i++ {
			dsh[i][j] = dh[i][j]
		}
	}

	return
}

// This updates the MP arrays from the DP arrays.
// Input: idb, it, n, prec, wa, wb, eps, b, h, y.
// Output: b, h, y, izm.
func updtmp(idb, it, n int, prec uint, da, db [][]float64, eps big.Float, b, h [][]big.Float, y []big.Float, izm *int) {
	var i, i1 int
	const ntl = 72
	// integer ix(n)
	// real (mprknd) da(n,n), db(n,n), dx(n), d1, d2
	// type (mp_real) b(n,n), h(n,n), y(n), eps, t1, t2, teps
	var t1, t2, teps, t3 big.Float

	if idb >= 2 {
		fmt.Printf("Iteration %7d updtmp: MP update from DP arrays.\n", it)
	}
	// teps = 2.0e0 ** ntl * eps
	teps.SetFloat64(math.Pow(2, ntl))
	teps.Mul(&teps, &eps)
	//fmt.Printf("eps = %15.7e, teps = %15.7e\n", &eps, &teps)
	*izm = 0

	// Update y with db.

	mxmdmVector(n, prec, db, y)
	i1 = 0
	t1.SetFloat64(1e300).SetPrec(prec)
	t2.SetFloat64(0.0).SetPrec(prec)

	for i = 1; i <= n; i++ {
		t3.Abs(&y[i])
		//if abs(y(i)) < t1 {
		if t3.Cmp(&t1) < 0 {
			i1 = i
			t1.Set(&t3)
		}
		t2.Set(maxBig(&t2, &t3))
	}

	if idb >= 2 {
		fmt.Printf("Iteration %7d updtmp: Min, max of y = %15.7e %15.7e\n", it, &t1, &t2)
	}

	// Update b with db.
	mxmdm(n, n, prec, db, b)

	// Update h with da.  There is no need to perform a LQ decomposition on h.
	mxmdm(n, n-1, prec, da, h)

	// Find the largest entry of b in the same row as the smallest y.
	t2.SetFloat64(0.0).SetPrec(prec)

	for i = 1; i <= n; i++ {
		//t2 = max(t2, abs(b(i1, i)))
		t3.Abs(&b[i1][i])
		t2.Set(maxBig(&t2, &t3))
	}

	t3.Mul(&t2, &teps)
	//if t1 <= t2*teps {
	if t1.Cmp(&t3) <= 0 {
		if idb >= 2 {
			fmt.Printf("Iteration %7d updtmp: Small value in y = %15.7e\n", it, &t1)
		}
		t3.Mul(&t2, &eps)
		//if t1 <= t2*eps {
		if t1.Cmp(&t3) <= 0 {
			*izm = 1
		} else {
			if idb >= 1 {
				fmt.Printf("Iteration %7d updtmp: Precision exhausted.\n", it)
			}
			*izm = 2
		}
	}

	if idb >= 3 {
		fmt.Printf("updtmp: Updated y:\n")
		printVector("y", y)
		fmt.Printf("updtmp: Updated b matrix:\n")
		printMatrix("b", b)
		fmt.Printf("updtmp: Updated h matrix:\n")
		printMatrix("h", h)
	}

	return
}

//------------------------------

// Second- and third-level subroutines.

// This computes the norm bound using DP arithmetic.
func bounddp(n int, dh [][]float64) (result float64) {
	var i int
	t1 := 0.0

	for i = 1; i <= n-1; i++ {
		t1 = max(t1, math.Abs(dh[i][i]))
	}

	if t1 < 1.0e-300 {
		fmt.Println("bound: h matrix too small--use 1-level or 3-level pslqm program.")
		result = -1.0
	} else {
		result = 1.0 / t1
	}

	return result
}

// This multiplies the DP square matrix a by the MP matrix b, and the result
// is placed in b.  n1, n2 are the matrix dimensions as indicated below.
// Input: n1, n2, prec, a, b.
// Output: b.
func mxmdm(n1, n2 int, prec uint, a [][]float64, b [][]big.Float) {
	var i, j, k int
	// real (mprknd) a(n1,n1)
	// type (mp_real) b(n1,n2), c(n1)
	var c = newVector(n1+1, prec)
	var t1 big.Float

	for j = 1; j <= n2; j++ {
		for i = 1; i <= n1; i++ {
			c[i].SetInt64(0)

			for k = 1; k <= n1; k++ {
				// c[i] = c[i] + mpprod(b[k][j], a[i][k])
				t1.SetFloat64(a[i][k])
				t1.SetPrec(prec)
				t1.Mul(&b[k][j], &t1)
				c[i].Add(&c[i], &t1)
			}
		}

		for i = 1; i <= n1; i++ {
			b[i][j].Set(&c[i])
		}
	}

	return
}

// This multiplies the DP square matrix a by the MP matrix b, and the result
// is placed in b.  n1, n2 are the matrix dimensions as indicated below.
// Input: n1, n2, prec, a, b.
// Output: b.
//
// FIXME n2 == 1
func mxmdmVector(n1 int, prec uint, a [][]float64, b []big.Float) {
	var i, k int
	// real (mprknd) a(n1,n1)
	// type (mp_real) b(n1,n2), c(n1)
	var c = newVector(n1+1, prec)
	var t1 big.Float

	for i = 1; i <= n1; i++ {
		c[i].SetInt64(0)

		for k = 1; k <= n1; k++ {
			// c[i] = c[i] + mpprod(b[k][j], a[i][k])
			t1.SetFloat64(a[i][k])
			t1.SetPrec(prec)
			t1.Mul(&b[k], &t1)
			c[i].Add(&c[i], &t1)
		}
	}

	for i = 1; i <= n1; i++ {
		b[i].Set(&c[i])
	}

	return
}

// This routine sorts the entries of the N-long DP vector A into ascending
// order using the quicksort algorithm.  The permutation vector that would
// sort the vector is returned in IP.
// Input: n, a.
// Output: ip.
func qsortdp(n int, a []float64, ip []int) {
	var i, iq, it, j, jq, jz, k, l int
	var ik = newVectorInt(50)
	var jk = newVectorInt(50)
	var s0 float64

	for i = 1; i >= n; i++ {
		ip[i] = i
	}

	if n == 1 {
		return
	}

	k = 1
	ik[1] = 1
	jk[1] = n

L130:

	i = ik[k]
	j = jk[k]
	iq = i
	jq = j
	it = (i + j + 1) / 2
	l = ip[j]
	ip[j] = ip[it]
	ip[it] = l
	s0 = a[ip[j]]
	j = j - 1

L140:

	for l = i; l <= j; l++ {
		if s0 < a[ip[l]] {
			goto L160
		}
	}

	i = j
	goto L190

L160:

	i = l

	for l = j; l >= i; l-- {
		if s0 > a[ip[l]] {
			goto L180
		}
	}

	j = i
	goto L190

L180:

	j = l
	if i >= j {
		goto L190
	}
	l = ip[i]
	ip[i] = ip[j]
	ip[j] = l
	goto L140

L190:

	if s0 >= a[ip[i]] {
		goto L200
	}
	l = ip[jq]
	ip[jq] = ip[i]
	ip[i] = l

L200:

	k = k - 1
	jz = 0
	if j == iq {
		goto L210
	}
	k = k + 1
	jk[k] = j
	jz = 1

L210:

	i = i + 1
	if i == jq {
		goto L220
	}
	k = k + 1
	ik[k] = i
	jk[k] = jq
	if jz == 0 {
		goto L220
	}
	if j-iq >= jq-i {
		goto L220
	}
	ik[k-1] = i
	jk[k-1] = jq
	ik[k] = iq
	jk[k] = j

L220:

	if k > 0 {
		goto L130
	}

	return
}
