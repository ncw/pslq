package main

import (
	"bytes"
	"math/bits"
	"os"
	"sort"
	"strings"
	"testing"
)

var (
	testMainIn = `# (3*math.Pi + 4*math.E/7)
10.978081862745977
# math.Pi
3.141592653589793
# math.E
2.718281828459045
# math.Log(2)
0.6931471805599453
0.28917320090206799499
0.57591529756863646394
0.55698607277729539344
0.54073048514703925260
0.99835889431176827458
0.11551877481656358526
`
	testMainWant = `Using precision 64
(3*math.Pi + 4*math.E/7) = 10.97808186274597699959
math.Pi = 3.14159265358979299999
math.E = 2.71828182845904499994
math.Log(2) = 0.69314718055994530000
x[4] = 0.28917320090206799499
x[5] = 0.57591529756863646394
x[6] = 0.55698607277729539344
x[7] = 0.54073048514703925260
x[8] = 0.99835889431176827458
x[9] = 0.11551877481656358526

Result with 3 terms is:
+7 * 10.97808186274597699959
-21 * 3.14159265358979299999
-4 * 2.71828182845904499994
+7 * (3*math.Pi + 4*math.E/7) -21 * math.Pi -4 * math.E = 0
Vector = [ 7, -21, -4, 0, 0, 0, 0, 0, 0, 0 ]
Result is accurate to 5.9978e-15
`
)

func TestMain(t *testing.T) {
	oldArgs := os.Args
	defer func() {
		os.Args = oldArgs
	}()
	stdin = bytes.NewBufferString(testMainIn)
	buf := new(bytes.Buffer)
	stdout = buf
	os.Args = []string{"pslq", "-prec", "64", "-"}
	main()
	got := strings.TrimSpace(buf.String())
	testMainWant = strings.TrimSpace(testMainWant)
	if got != testMainWant {
		t.Error("expecting", testMainWant)
		t.Error("got", got)
	}
}

func testNext(t *testing.T, b int, expected []uint64) {
	n := uint64(0)
	for i, want := range expected {
		if want != n {
			t.Errorf("Want[%d] %06b got %06b", i, want, n)
		}
		n = next(n, b)
	}
}

func TestNext1(t *testing.T) {
	testNext(t, 1, []uint64{
		0b0,
		0b1,
		0b0,
		0b1,
		0b0,
	})
}

func TestNext2(t *testing.T) {
	testNext(t, 2, []uint64{
		0b00,
		0b01,
		0b10,
		0b11,
		0b00,
		0b01,
		0b10,
		0b11,
	})
}

func TestNext3(t *testing.T) {
	testNext(t, 3, []uint64{
		0b000,
		0b001,
		0b010,
		0b100,
		0b011,
		0b101,
		0b110,
		0b111,
		0b000,
	})
}

func TestNext6(t *testing.T) {
	testNext(t, 6, []uint64{
		0b000000,
		0b000001,
		0b000010,
		0b000100,
		0b001000,
		0b010000,
		0b100000,
		0b000011,
		0b000101,
		0b000110,
		0b001001,
		0b001010,
		0b001100,
		0b010001,
		0b010010,
		0b010100,
		0b011000,
		0b100001,
		0b100010,
		0b100100,
		0b101000,
		0b110000,
		0b000111,
		0b001011,
		0b001101,
		0b001110,
		0b010011,
		0b010101,
		0b010110,
		0b011001,
		0b011010,
		0b011100,
		0b100011,
		0b100101,
		0b100110,
		0b101001,
		0b101010,
		0b101100,
		0b110001,
		0b110010,
		0b110100,
		0b111000,
		0b001111,
		0b010111,
		0b011011,
		0b011101,
		0b011110,
		0b100111,
		0b101011,
		0b101101,
		0b101110,
		0b110011,
		0b110101,
		0b110110,
		0b111001,
		0b111010,
		0b111100,
		0b011111,
		0b101111,
		0b110111,
		0b111011,
		0b111101,
		0b111110,
		0b111111,
		0b000000,
	})
}

func TestNext20(t *testing.T) {
	ns := make([]uint64, 0, 1<<20)
	for j := 1; j <= 20; j++ {
		n := uint64(0)
		ns := ns[:0]
		for {
			ns = append(ns, n)
			n = next(n, j)
			if n == 0 {
				break
			}
		}
		if len(ns) != 1<<j {
			t.Errorf("Expecting %d iterations but got %d", 1<<j, len(ns))
		}
		// Check the bit count keeps going up
		bitCount := 0
		for i, n := range ns {
			newBitCount := bits.OnesCount64(n)
			if newBitCount < bitCount {
				t.Errorf("Bitcount got smaller at %d from %d to %d", i, bitCount, newBitCount)
			}
			bitCount = newBitCount
		}
		sort.Slice(ns, func(i, j int) bool {
			return ns[i] < ns[j]
		})
		// Check we covered all the numbers
		for i, n := range ns {
			if uint64(i) != n {
				t.Errorf("Expecting %d in position %d for j=%d", n, i, j)
			}
		}
	}
}
