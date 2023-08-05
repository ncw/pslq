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
	"math/rand"
	"os"
	"runtime"
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
	tryAllRandom              = flag.Bool("try-all-random", false, "If set, uses random masks for -try-all")
	workers                   = flag.Int("workers", runtime.NumCPU(), "Use this many threads in -try-all")
	algorithm                 = flag.Int("algorithm", 1, "Which algorithm to use. 1: orig, 2: pslqm2")
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

The comment immediately before a the number will be taken as its name
if it is present and printed in the results.

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

// Read lines from in as big.Float with a name
func read(in io.Reader, xs []big.Float, names []string) ([]big.Float, []string) {
	scanner := bufio.NewScanner(in)
	name := ""
	for scanner.Scan() {
		var x big.Float
		x.SetPrec(*prec)
		text := strings.TrimSpace(scanner.Text())
		if len(text) == 0 {
			continue
		}
		if text[0] == '#' {
			name = strings.TrimSpace(text[1:])
			continue
		}
		_, ok := x.SetString(text)
		if !ok {
			log.Fatalf("Failed to parse line %q", text)
		}
		xs = append(xs, x)
		if name == "" {
			name = fmt.Sprintf("x[%d]", len(xs)-1)
		}
		names = append(names, name)
		name = ""
	}
	if err := scanner.Err(); err != nil {
		log.Fatalf("Error reading input: %v", err)
	}
	return xs, names
}

// Read lines from the file name as big.Float
//
// if name is '-' then reads from stdin
func readFile(name string, xs []big.Float, names []string) ([]big.Float, []string) {
	if name == "-" {
		return read(stdin, xs, names)
	}
	in, err := os.Open(name)
	if err != nil {
		log.Fatalf("Error opening file %q: %v", name, err)
	}
	defer in.Close()
	return read(in, xs, names)
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
func printResults(p *pslq.Pslq, xs []big.Float, names []string, result []big.Int) bool {
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
	fmt.Fprintf(&out, "\nResult with %d terms is:\n", terms)
	// Print numerically
	for i := range result {
		d := &result[i]
		if d.Sign() == 0 {
			continue
		}
		fmt.Fprintf(&out, "%+d * %.*f\n", d, digits, &xs[i])
	}
	// Print symbolically
	for i := range result {
		d := &result[i]
		if d.Sign() == 0 {
			continue
		}
		if i != 0 && *needFirst {
			tmp := new(big.Int)
			tmp.Neg(d)
			d = tmp
		}
		fmt.Fprintf(&out, "%+d * %s ", d, names[i])
		if i == 0 && *needFirst {
			fmt.Fprintf(&out, "= ")
		}
	}
	if !*needFirst {
		fmt.Fprintf(&out, "= 0")
	}
	fmt.Fprintln(&out)
	// Print vector
	fmt.Fprintf(&out, "Vector = [ ")
	for i := range result {
		d := &result[i]
		fmt.Fprintf(&out, "%d", d)
		if i != len(result)-1 {
			fmt.Fprintf(&out, ", ")
		}
	}
	fmt.Fprintf(&out, " ]\n")
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
func run(p *pslq.Pslq, xs []big.Float, names []string) error {
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
				namesCopy := append([]string(nil), names[:i]...)
				namesCopy = append(namesCopy, names[i+1:]...)
				err := run(p, xsCopy, namesCopy)
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
	printResults(p, xs, names, result)
	return nil
}

// Iterate through numbers < 1<<b in number of bits set order
//
// # When the iteration is finished it returns 0
//
// Given an n (current number) and b (highest bit to be set) this
// returns the next number. This will have the same number of bits set
// as n unless there are no more numbers with that many bits.
func next(n uint64, b int) uint64 {
	if n == 0 {
		return 1
	}
	lo := n & -n           // lowest one bit
	lz := (n + lo) & ^n    // lowest zero bit above lo
	n |= lz                // add lz to the set
	n &= ^(lz - 1)         // reset bits below lz
	n |= (lz / lo / 2) - 1 // put back right number of bits at end
	// If we exceed bits, start next count
	if n > ((1 << b) - 1) {
		newBits := bits.OnesCount64(n) + 1
		if newBits > b {
			// reset
			return 0
		}
		n = (1 << newBits) - 1
	}
	return n
}

// Do an All run of pslq with xs
func runTryAll(p *pslq.Pslq, xs []big.Float, names []string) error {
	const statsPrintTime = 10 * time.Second
	start := time.Now()
	nextStat := time.Now().Add(statsPrintTime)
	if len(xs) >= 64 {
		return errors.New("Can't have 64 or more items with -try-all")
	}
	trials := uint64(1) << len(xs)
	trialsBits := len(xs)
	if *needFirst {
		trials >>= 1
		trialsBits -= 1
	}
	found := uint64(0)
	total := uint64(0)
	badRelation := uint64(0)
	random := rand.New(rand.NewSource(time.Now().UnixNano()))

	// Create worker routines
	in := make(chan uint64, *workers)
	var wg sync.WaitGroup
	for k := 0; k < *workers; k++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			xsCopy := make([]big.Float, len(xs))
			namesCopy := make([]string, len(names))
			index := make([]int, len(xs))
			unscrambledResults := make([]big.Int, len(xs))
			for i := range in {
				// fmt.Printf("Starting run with mask %0*b\n", trialsBits, i)
				xsCopy = xsCopy[:0]
				namesCopy = namesCopy[:0]
				index = index[:0]
				for j := range xs {
					if ((1 << j) & i) == 0 {
						xsCopy = append(xsCopy, xs[j])
						namesCopy = append(namesCopy, names[j])
						index = append(index, j)
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
					atomic.AddUint64(&badRelation, 1)
					continue
				}
				for i := range unscrambledResults {
					unscrambledResults[i] = big.Int{}
				}
				for i, j := range index {
					unscrambledResults[j] = result[i]
				}
				atomic.AddUint64(&total, 1)
				//if printResults(p, xsCopy, namesCopy, result) {
				if printResults(p, xs, names, unscrambledResults) {
					atomic.AddUint64(&found, 1)
				}
			}
		}()
	}

	// Bitmask for each trial
	mask := uint64(0)
	for i := uint64(0); i < trials; i++ {
		now := time.Now()
		if now.After(nextStat) {
			dt := time.Since(start)
			iterationsPerSecond := float64(i) / float64(dt) * float64(time.Second)
			eta := time.Duration(float64(trials-uint64(i))/iterationsPerSecond) * time.Second
			fmt.Fprintf(stdout, "Iteration %d/%d iterations/second %.2f eta %v, Unique %d/%d, Bad %d\n", i, trials, iterationsPerSecond, eta, atomic.LoadUint64(&found), atomic.LoadUint64(&total), atomic.LoadUint64(&badRelation))
			nextStat = nextStat.Add(statsPrintTime)
		}
		// If needFirst is set we always want the first item in the mask
		workerMask := mask
		if *needFirst {
			workerMask = mask << 1
		}
		// Get the workers to calculate the mask
		in <- workerMask
		if *tryAllRandom {
			mask = random.Uint64() & (trials - 1)
		} else {
			mask = next(mask, trialsBits)
		}
	}
	// Signal to workers they are finished
	close(in)
	// Wait for workers to exit
	wg.Wait()
	fmt.Fprintf(stdout, "\nFound %d unique results out of %d total results\n", found, total)
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
	var names []string
	for _, arg := range args {
		xs, names = readFile(arg, xs, names)
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
		fmt.Fprintf(stdout, "%s = %.*f\n", names[i], digits, x)
	}

	// max coefficient is 2^logMaxCoeff
	maxCoeff := big.NewInt(1)
	maxCoeff.Lsh(maxCoeff, *logMaxCoeff)

	pslq := pslq.New(*prec).SetMaxSteps(*iterations).SetVerbose(*verbose).SetMaxCoeff(maxCoeff).SetTarget(uint(float64(*prec) * (*targetPrecision))).SetAlgorithm(*algorithm)
	var err error
	if *tryAll {
		err = runTryAll(pslq, xs, names)
	} else {
		err = run(pslq, xs, names)
	}
	if err != nil {
		fmt.Fprintf(os.Stderr, "PSLQ failed: %v\n", err)
		os.Exit(1)
	}
}
