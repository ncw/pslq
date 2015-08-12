package main

import (
	"bytes"
	"os"
	"testing"
)

const (
	in = `# 3*math.Pi + 4*math.E/7
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
	out = `Using precision 64
x[0] =    10.97808186274597699959
x[1] =     3.14159265358979299999
x[2] =     2.71828182845904499994
x[3] =     0.69314718055994530000
x[4] =     0.28917320090206799499
x[5] =     0.57591529756863646394
x[6] =     0.55698607277729539344
x[7] =     0.54073048514703925260
x[8] =     0.99835889431176827458
x[9] =     0.11551877481656358526
Result is
7 * 10.97808186274597699959
-21 * 3.14159265358979299999
-4 * 2.71828182845904499994
`
)

func TestMain(t *testing.T) {
	oldArgs := os.Args
	defer func() {
		os.Args = oldArgs
	}()
	stdin = bytes.NewBufferString(in)
	buf := new(bytes.Buffer)
	stdout = buf
	os.Args = []string{"pslq", "-prec", "64", "-"}
	main()
	got := buf.String()
	if got != out {
		t.Error("expecting", out)
		t.Error("got", got)
	}
}
