PSLQ
====

This is a library and a command line tool for using the PSLQ
integer relation algorithm.

It is used to find relations of the form

    a[0]*x[0] + a[1]*x[1] + a[2]*x[2] + ... a[n]*x[n] = 0

where `a[]` are floating point `big.Float` values to a given precision
and `x[]` are integer `big.Int` values.

The algorithm will tell you that a relation exists, or not within
certain bounds.

[![GoDoc](https://godoc.org/github.com/ncw/rclone?status.svg)](https://godoc.org/github.com/ncw/pslq)
[![Build Status](https://travis-ci.org/ncw/pslq.png)](https://travis-ci.org/ncw/pslq)


Install
-------

Install using go get

    go get github.com/ncw/pslq/...

and this will build the `pslq` binary in `$GOPATH/bin`.

Using the library
-----------------

FIXME

Using the binary
----------------

FIXME

License
-------

This is free software under the terms of MIT the license (check the
COPYING file included in this package).

Portions of the code have been copied from the Go source.  These are
identified by comments at the head of each file and these are
Copyright (c) The Go Authors.  See the GO-LICENSE file for full details.

Contact and support
-------------------

The project website is at:

  * https://github.com/ncw/pslq

There you can file bug reports, ask for help or contribute patches.

Authors
-------

  * Nick Craig-Wood <nick@craig-wood.com>
