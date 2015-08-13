[![GoDoc](https://godoc.org/github.com/ncw/rclone?status.svg)](https://godoc.org/github.com/ncw/pslq)
[![Build Status](https://travis-ci.org/ncw/pslq.png)](https://travis-ci.org/ncw/pslq)

This is a library and a command line tool for using the
[PSLQ](http://mathworld.wolfram.com/PSLQAlgorithm.html) [integer
relation algorithm](https://en.wikipedia.org/wiki/Integer_relation_algorithm).

It is used to find relations of the form

    a[0]*x[0] + a[1]*x[1] + a[2]*x[2] + ... a[n]*x[n] = 0

where `a[]` are floating point `big.Float` values to a given precision
and `x[]` are integer `big.Int` values.

The algorithm will tell you that a relation exists, or not within
certain bounds.

NB This requires Go >= 1.5 for `big.Float` support.

Install
-------

Install using go get

    go get github.com/ncw/pslq/...

and this will install the library and build the `pslq` binary in
`$GOPATH/bin`.

Using the library
-----------------

See the [go doc](https://godoc.org/github.com/ncw/pslq) for details,
in particular the example.

Using the binary
----------------

Usage pslq [Options] <file>

Where file should contain decimal numbers, one per line.  White space
is ignored, as are comment lines starting with '#'.

If more than one file is passed in then they are concatenated

If file is '-' then stdin will be read

Options:
```
  -iterations int
    	Number of iterations to use max (default 1000)
  -prec uint
    	Precision to use (default 64)
  -verbose
    	Print lots of stuff while running
```

License
-------

This is free software under the terms of the MIT license (check the
LICENSE file included in this package).

Portions of the code have been ported from the SymPy source code.
These are identified in the relevant files and these parts are
additionally licensed under the SymPy Licence (see the SymPy-LICENSE
file).

Contact and support
-------------------

The project website is at:

  * https://github.com/ncw/pslq

There you can file bug reports, ask for help or contribute patches.

Authors
-------

  * SymPy Development Team - original code that was ported to Go
  * Nick Craig-Wood <nick@craig-wood.com>
