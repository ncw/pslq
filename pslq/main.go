// Read a file of real numbers (one per line) and run PSLQ on them
package main

import (
	"bufio"
	"flag"
	"fmt"
	"io"
	"log"
	"math"
	"math/big"
	"os"
	"strings"

	"github.com/ncw/pslq"
)

var (
	verbose              = flag.Bool("verbose", false, "Print lots of stuff while running")
	iterations           = flag.Int("iterations", 1000, "Number of iterations to use max")
	prec                 = flag.Uint("prec", 64, "Precision to use")
	stdin      io.Reader = os.Stdin
	stdout     io.Writer = os.Stdout
)

// syntaxError prints the syntax
func syntaxError() {
	fmt.Fprintf(os.Stderr, `pslq - find integer relations

Usage pslq [Options] <file>

Where file should contain decimal numbers, one per line.  White space
is ignored, as are comment lines starting with '#'.

If more than one file is passed in then they are concatenated

If file is '-' then stdin will be read

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
	digits := int(math.Log10(2)*float64(*prec) + 1)
	fmt.Fprintf(stdout, "Using precision %d\n", *prec)
	for i := range xs {
		x := &xs[i]
		fmt.Fprintf(stdout, "x[%d] = %*.*f\n", i, digits+6, digits, x)
	}

	pslq := pslq.New(*prec).SetMaxSteps(*iterations).SetVerbose(*verbose)
	result, err := pslq.Run(xs)
	if err != nil {
		fmt.Fprintf(os.Stderr, "PSLQ failed: %v\n", err)
		os.Exit(1)
	}
	fmt.Fprintf(stdout, "Result is\n")
	for i := range result {
		d := &result[i]
		if d.Sign() == 0 {
			continue
		}
		fmt.Fprintf(stdout, "%d * %.*f\n", d, digits, &xs[i])
	}
}
