// Maths utilities

package pslq

import (
	"fmt"
	"math"
	"math/big"
	"strings"
)

type num interface {
	~int | ~int64 | ~float64
}

func max[T num](a, b T) T {
	if a >= b {
		return a
	}
	return b
}

func min[T num](a, b T) T {
	if a <= b {
		return a
	}
	return b
}

func maxBig(a, b *big.Float) *big.Float {
	if a.Cmp(b) >= 0 {
		return a
	}
	return b
}

func minBig(a, b *big.Float) *big.Float {
	if a.Cmp(b) <= 0 {
		return a
	}
	return b
}

// Make a new matrix with that many rows and that many cols with the
// provided precision
func newMatrix(rows, cols int, prec uint) [][]big.Float {
	// Make a continuously allocated row with all the items
	U := make([]big.Float, rows*cols)
	M := make([][]big.Float, rows)
	for i := 0; i < rows; i++ {
		M[i] = U[cols*i : cols*(i+1) : cols*(i+1)]
		// M[i] = make([]big.Float, cols)
		for j := range M[i] {
			M[i][j].SetPrec(prec)
		}
	}
	return M
}

// Make a new matrix with that many rows and that many cols
func newMatrixFloat64(rows, cols int) [][]float64 {
	// Make a continuously allocated row with all the items
	U := make([]float64, rows*cols)
	M := make([][]float64, rows)
	for i := 0; i < rows; i++ {
		M[i] = U[cols*i : cols*(i+1) : cols*(i+1)]
	}
	return M
}

// Make a new matrix with that many rows and that many cols
func newBigIntMatrix(rows, cols int) [][]big.Int {
	// Make a continuously allocated row with all the items
	U := make([]big.Int, rows*cols)
	M := make([][]big.Int, rows)
	for i := 0; i < rows; i++ {
		M[i] = U[cols*i : cols*(i+1) : cols*(i+1)]
		// M[i] = make([]big.Int, cols)
	}
	return M
}

// Make a new vector with n items and the provided precision
func newVector(n int, prec uint) []big.Float {
	V := make([]big.Float, n)
	for i := range V {
		V[i].SetPrec(prec)
	}
	return V
}

// Make a new vector with n items
func newVectorFloat64(n int) []float64 {
	return make([]float64, n)
}

// Make a new vector with n items
func newVectorInt(n int) []int {
	return make([]int, n)
}

// Return how many decimal digits we should print given a given binary precision
func digits(prec uint) uint {
	return uint(math.Ceil(math.Log10(2) * float64(prec)))
}

// Print a matrix
func printMatrix(name string, X [][]big.Float) {
	n := len(X) - 1
	for i := 1; i <= n; i++ {
		for j := 1; j <= n; j++ {
			fmt.Printf("%s[%d,%d] = %.*f (prec = %d)\n", name, i, j, digits(X[i][j].Prec()), &X[i][j], X[i][j].Prec())
		}
		fmt.Printf("\n")
	}
}

// Print a matrix
// Row  1
//
//		1.00000D+00    0.00000D+00    0.00000D+00    0.00000D+00    0.00000D+00
//	    -1.86301D-01   -2.00193D-01   -2.16280D-01    9.96727D-01    0.00000D+00
//		0.00000D+00    0.00000D+00    0.00000D+00    0.00000D+00    0.00000D+00
//		0.00000D+00    0.00000D+00    0.00000D+00    0.00000D+00    0.00000D+00
/*
   4.13127D-01    0.00000D+00    0.00000D+00    0.00000D+00    0.00000D+00
   0.00000D+00    0.00000D+00    0.00000D+00    0.00000D+00    0.00000D+00
   0.00000D+00    0.00000D+00    0.00000D+00    0.00000D+00    0.00000D+00
   0.00000D+00    0.00000D+00    0.00000D+00    0.00000D+00

   4.13127E-01    0.00000E+00    0.00000E+00    0.00000E+00    0.00000E+00
   0.00000E+00    0.00000E+00    0.00000E+00    0.00000E+00    0.00000E+00
   0.00000E+00    0.00000E+00    0.00000E+00    0.00000E+00    0.00000E+00
   0.00000E+00    0.00000E+00    0.00000E+00    0.00000E+00
*/
func printMatrixFloat64(name string, X [][]float64) {
	// n := len(X) - 1
	// for i := 1; i <= n; i++ {
	// 	for j := 1; j <= n; j++ {
	// 		fmt.Printf("%s[%d,%d] = %.*f\n", name, i, j, digits(53), X[i][j])
	// 	}
	// 	fmt.Printf("\n")
	// }

	//fmt.Println(name)
	n := len(X) - 1
	colMax := n
	if name == "dh" {
		colMax -= 1
	}
	var out strings.Builder
	for i := 1; i <= n; i++ {
		fmt.Fprintf(&out, "Row%3d\n", i)
		for j := 1; j <= colMax; j += 5 {
			jEnd := min(j+5, colMax+1)
			for _, v := range X[i][j:jEnd] {
				fmt.Fprintf(&out, "%15.5E", v)
			}
			//fmt.Printf("%s[%d,%d] = %.*f\n", name, i, j, digits(53), X[i][j])
			fmt.Fprintf(&out, "\n")
		}
	}
	fmt.Print(strings.ReplaceAll(out.String(), "E", "D"))
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
		fmt.Printf("%s[%d] = %.*f (prec = %d)\n", name, i, digits(x[i].Prec()), &x[i], x[i].Prec())
	}
}

// Print a vector
func printVectorFloat64(name string, x []float64) {
	// for i := range x {
	// 	if i == 0 {
	// 		continue
	// 	}
	// 	fmt.Printf("%s[%d] = %.*f\n", name, i, digits(53), x[i])
	// }

	// fmt.Println(name)
	n := len(x) - 1
	var out strings.Builder
	for i := 1; i <= 1; i++ {
		fmt.Fprintf(&out, "Row%3d\n", i)
		for j := 1; j <= n; j += 5 {
			jEnd := min(j+5, n+1)
			for _, v := range x[j:jEnd] {
				fmt.Fprintf(&out, "%15.5E", v)
			}
			//fmt.Printf("%s[%d,%d] = %.*f\n", name, i, j, digits(53), X[i][j])
			fmt.Fprintf(&out, "\n")
		}
	}
	fmt.Print(strings.ReplaceAll(out.String(), "E", "D"))
}

// Print a vector
func printVectorInt(name string, x []int) {
	for i := range x {
		if i == 0 {
			continue
		}
		fmt.Printf("%s[%d] = %d\n", name, i, x[i])
	}
}
