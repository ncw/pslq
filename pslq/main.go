// Read a file of real numbers (one per line) and run PSLQ on them
package main

import (
	"bufio"
	"errors"
	"flag"
	"fmt"
	"io"
	"log"
	"math"
	"math/big"
	"math/bits"
	"os"
	"runtime"
	"sort"
	"strings"
	"sync"
	"sync/atomic"
	"time"

	"github.com/ncw/pslq"
)

var (
	verbose                   = flag.Bool("verbose", false, "Print lots of stuff while running")
	iterations                = flag.Int("iterations", 1000, "Number of iterations to use max")
	prec                      = flag.Uint("prec", 64, "Precision to use (bits)")
	logMaxCoeff               = flag.Uint("log-max-coeff", 64, "log2(max coefficient size)")
	needFirst                 = flag.Bool("need-first", false, "Retry if first entry is not used")
	targetPrecision           = flag.Float64("target-precision", 0.75, "Target precision of the result as fraction of prec")
	tryAll                    = flag.Bool("try-all", false, "Try all combinations of input until solution found")
	workers                   = flag.Int("workers", runtime.NumCPU(), "Use this many threads in -try-all")
	stdin           io.Reader = os.Stdin
	stdout          io.Writer = os.Stdout
	digits          int
)

// syntaxError prints the syntax
func syntaxError() {
	fmt.Fprintf(os.Stderr, `pslq - find integer relations

Usage pslq [Options] <file>

Where file should contain decimal numbers, one per line.  White space
is ignored, as are comment lines starting with '#'.

If more than one file is passed in then they are concatenated

If file is '-' then stdin will be read

If -need-first is set then it will retry without items in the
list of numbers until it finds a match with the first item included.

If -try-all is set it will try all possible combinations of the input
items as the PSLQ algorithm seems sensitive to exactly what and how
many items are presented.

Options:
`)
	flag.PrintDefaults()
}

// Exit with the message
func fatalf(message string, args ...interface{}) {
	syntaxError()
	fmt.Fprintf(os.Stderr, message, args...)
	os.Exit(1)
}

// Read lines from in as big.Float
func read(in io.Reader, xs []big.Float) []big.Float {
	scanner := bufio.NewScanner(in)
	for scanner.Scan() {
		var x big.Float
		x.SetPrec(*prec)
		text := strings.TrimSpace(scanner.Text())
		if len(text) == 0 || text[0] == '#' {
			continue
		}
		_, ok := x.SetString(text)
		if !ok {
			log.Fatalf("Failed to parse line %q", text)
		}
		xs = append(xs, x)
	}
	if err := scanner.Err(); err != nil {
		log.Fatalf("Error reading input: %v", err)
	}
	return xs
}

// Read lines from the file name as big.Float
//
// if name is '-' then reads from stdin
func readFile(name string, xs []big.Float) []big.Float {
	if name == "-" {
		return read(stdin, xs)
	}
	in, err := os.Open(name)
	if err != nil {
		log.Fatalf("Error opening file %q: %v", name, err)
	}
	defer in.Close()
	return read(in, xs)
}

// Sum the relation PSLQ found
func sumResult(xs []big.Float, result []big.Int) *big.Float {
	prec := xs[0].Prec()
	var sum = new(big.Float)
	sum.SetPrec(prec)
	for i := range xs {
		if result[i].Sign() == 0 {
			continue
		}
		var tmp big.Float
		tmp.SetPrec(prec)
		tmp.SetInt(&result[i])
		tmp.Mul(&tmp, &xs[i])
		sum.Add(sum, &tmp)
	}
	return sum
}

var (
	printedResults = make(map[string]struct{})
	printResultsMu sync.Mutex
)

// Returns true if the result hasn't been seen before
func printResults(p *pslq.Pslq, xs []big.Float, result []big.Int) bool {
	printResultsMu.Lock()
	defer printResultsMu.Unlock()
	var out strings.Builder
	terms := 0
	for i := range result {
		d := &result[i]
		if d.Sign() != 0 {
			terms++
		}
	}
	fmt.Fprintf(&out, "Result with %d terms is:\n", terms)
	for i := range result {
		d := &result[i]
		if d.Sign() == 0 {
			continue
		}
		fmt.Fprintf(&out, "%d * %.*f\n", d, digits, &xs[i])
	}
	if _, found := printedResults[out.String()]; !found {
		fmt.Fprintf(stdout, "%s", out.String())
		fmt.Fprintf(stdout, "Result is accurate to %.5g\n", sumResult(xs, result))
		printedResults[out.String()] = struct{}{}
	} else {
		//fmt.Fprintf(stdout, "Duplicate found\n")
		return false
	}
	return true
}

// Do a single run of pslq with xs
func run(p *pslq.Pslq, xs []big.Float) error {
	result, err := p.Run(xs)
	if err != nil {
		return err
	}
	if *needFirst && result[0].Sign() == 0 {
		// Need the first item in the results, so retry
		// without each item in result in turn
		for i := range result {
			d := &result[i]
			if d.Sign() != 0 {
				// xs without i
				xsCopy := append([]big.Float(nil), xs[:i]...)
				xsCopy = append(xsCopy, xs[i+1:]...)
				err := run(p, xsCopy)
				fmt.Printf("xs[%d] %v\n", len(xsCopy), err)
				if err == nil {
					// Have printed a result already so return
					return nil
				}
				if err == pslq.ErrorPrecisionExhausted {
					return err
				}
			}
		}
		return errors.New("couldn't find solution with the first item")
	}
	printResults(p, xs, result)
	return nil
}

// Do an All run of pslq with xs
func runTryAll(p *pslq.Pslq, xs []big.Float) error {
	const statsPrintTime = 10 * time.Second
	start := time.Now()
	nextStat := time.Now().Add(statsPrintTime)
	if len(xs) >= 64 {
		return errors.New("Can't have 64 or more items with -try-all")
	}
	trials := uint64(1) << len(xs)
	if *needFirst {
		trials >>= 1
	}
	found := uint64(0)

	// Create worker routines
	in := make(chan uint64, *workers)
	var wg sync.WaitGroup
	for k := 0; k < *workers; k++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			xsCopy := make([]big.Float, len(xs))
			for i := range in {
				xsCopy = xsCopy[:0]
				for j := range xs {
					if ((1 << j) & i) == 0 {
						xsCopy = append(xsCopy, xs[j])
					}
				}
				result, err := p.Run(xsCopy)
				if err == pslq.ErrorPrecisionExhausted {
					log.Fatal(err)
				}
				if err != nil {
					continue
				}
				if *needFirst && result[0].Sign() == 0 {
					continue
				}
				if printResults(p, xsCopy, result) {
					atomic.AddUint64(&found, 1)
				}
			}
		}()
	}

	// Make a bitmask for each trial
	// If *needFirst then lowest bit is always 0 (included)
	masks := make([]uint64, trials)
	for i := uint64(0); i < trials; i++ {
		if *needFirst {
			masks[i] = i << 1
		} else {
			masks[i] = i
		}
	}

	// Sort these bitmasks by number of bits set
	//
	// There is probably a cool way to produce this sequence
	// directly without writing it out and sorting but I couldn't
	// think of it!
	sort.Slice(masks, func(i, j int) bool {
		mi, mj := masks[i], masks[j]
		diffOnes := bits.OnesCount64(mi) - bits.OnesCount64(mj)
		if diffOnes != 0 {
			return diffOnes < 0
		}
		return mi < mj
	})

	for i, mask := range masks {
		now := time.Now()
		if now.After(nextStat) {
			dt := time.Since(start)
			iterationsPerSecond := float64(i) / float64(dt) * float64(time.Second)
			eta := time.Duration(float64(trials-uint64(i))/iterationsPerSecond) * time.Second
			fmt.Fprintf(stdout, "Iteration %d/%d iterations/second %.2f eta %v\n", i, trials, iterationsPerSecond, eta)
			nextStat = nextStat.Add(statsPrintTime)
		}
		// Get the workers to calculate the mask
		in <- mask
	}
	// Signal to workers they are finished
	close(in)
	// Wait for workers to exit
	wg.Wait()
	fmt.Fprintf(stdout, "Found %d results\n", found)
	if found == 0 {
		return errors.New("No results found with -try-all")
	}
	return nil
}

func main() {
	flag.Usage = syntaxError
	flag.Parse()
	args := flag.Args()
	if len(args) < 1 {
		fatalf("No input supplied\n")
	}
	var xs []big.Float
	for _, arg := range args {
		xs = readFile(arg, xs)
	}
	// Set the max precision
	if *prec == 0 {
		for i := range xs {
			x := &xs[i]
			if x.Prec() > *prec {
				*prec = x.Prec()
			}
		}
	}
	digits = int(math.Log10(2)*float64(*prec) + 1)
	fmt.Fprintf(stdout, "Using precision %d\n", *prec)
	for i := range xs {
		x := &xs[i]
		fmt.Fprintf(stdout, "x[%d] = %*.*f\n", i, digits+6, digits, x)
	}

	// max coefficient is 2^logMaxCoeff
	maxCoeff := big.NewInt(1)
	maxCoeff.Lsh(maxCoeff, *logMaxCoeff)

	pslq := pslq.New(*prec).SetMaxSteps(*iterations).SetVerbose(*verbose).SetMaxCoeff(maxCoeff).SetTarget(uint(float64(*prec) * (*targetPrecision)))
	var err error
	if *tryAll {
		err = runTryAll(pslq, xs)
	} else {
		err = run(pslq, xs)
	}
	if err != nil {
		fmt.Fprintf(os.Stderr, "PSLQ failed: %v\n", err)
		os.Exit(1)
	}
}
