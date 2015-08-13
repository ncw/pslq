// Example for psql module
package pslq_test

import (
	"fmt"
	"log"
	"math"
	"math/big"

	"github.com/ncw/pslq"
)

func Example() {
	in := make([]big.Float, 4)
	in[0].SetFloat64(3*math.Pi + 4*math.E/7)
	in[1].SetFloat64(math.Pi)
	in[2].SetFloat64(math.E)
	in[3].SetFloat64(math.Log(2))
	pslq := pslq.New(64).SetMaxSteps(1000)
	out, err := pslq.Run(in)
	if err != nil {
		log.Fatal(err)
	}
	for i := range out {
		d := &out[i]
		if d.Sign() != 0 {
			fmt.Printf("%d * %.10f\n", d, &in[i])
		}
	}
	// Output: 7 * 10.9780818627
	// -21 * 3.1415926536
	// -4 * 2.7182818285
}
