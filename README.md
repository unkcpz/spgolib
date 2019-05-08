# spgolib: spglib Go interface [![GoDoc](https://godoc.org/github.com/unkcpz/spgolib?status.svg)](https://godoc.org/github.com/unkcpz/spgolib)   [![Build Status](https://travis-ci.org/unkcpz/spgolib.svg)](https://travis-ci.org/unkcpz/spgolib)

Go wrapper for spglib C-API, dependent on "gonum.org/v1/gonum/mat". Spglib's C code
is embedded in this package, resulting binaries are statically linked.
You can just `go get -x github.com/unkcpz/spgolib` (-x to show what it's doing, compilation takes long due to C files).

## Copyright

This package using MIT, and C code from https://github.com/atztogo/spglib with BSD-3-Clause below.

Copyright (c) 2014, Atsushi Togo
All rights reserved.
