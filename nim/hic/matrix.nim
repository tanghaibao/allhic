#
#
#            Nim Simple Matrix Module
#        (c) Copyright 2015 Tom Krauss
#
# This module implements a linear algebra (matrix) class.
#
# The algorithms in this module are not necessarily the fastest
# possible nor the most memory efficient.  Some optimizations have been
# made, however.
#
# Limitations
# -----------
#   -  supports only dense matrices
#
#   -  not the most memory efficient algorithms (returns new
#      matrices rather than performing operations in place)
#
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
import
  math,
  complex


type
  Matrix*[T] = object
    numRows:   int
    numCols:   int
    data:       seq[T]

{.deprecated: [TMatrix: Matrix].}

const
  EPS = 5.0e-10


proc index[T](x: Matrix[T], r,c: int): int {.inline.} =
  ## Internal element access.  Returns the matrix element at
  ## row=r, column=c.  Element access is _unchecked_.
  return r*x.numCols+c


proc rows*[T](x: Matrix[T]): int {.inline.} =
  ## Returns the number of rows in the matrix `x`.
  result = x.numRows

proc cols*[T](x: Matrix[T]): int {.inline.}  =
  ## Returns the number of columns in the matrix `x`.
  result = x.numCols



proc newMatrix*[T](rows, cols: int, d: openarray[T]): Matrix[T] =
  ## Constructor.  Initializes the matrix by allocating memory
  ## for the data and setting the number of rows and columns
  ## Sets the data to the values specified in `d`.
  result.numRows = rows
  result.numCols = cols
  newSeq(result.data, rows*cols)
  if len(d)>0:
    if len(d)<(rows*cols):
      raise newException(IndexError, "insufficient data supplied in matrix constructor")

    for i in countup(0,rows*cols-1):
      result.data[i] = d[i]


proc newMatrix*[T](rows, cols: int): Matrix[T] =
  ## Constructor.  Initializes the matrix by allocating memory
  ## for the data and setting the number of rows and columns.
  ## Initially populated with 0.
  result.numRows = rows
  result.numCols = cols
  newSeq(result.data, rows*cols)
  when not (T is Complex):
    for i in countup(0,rows*cols-1):
      result.data[i] = T(0)
  when (T is Complex):
    for i in countup(0,rows*cols-1):
      result.data[i] = (0.0, 0.0)



proc eye*[T](N: int): Matrix[T] =
  ## Returns the NxN square identity matrix.
  result = newMatrix[T](N,N)
  when not (T is Complex):
    for i in countup(0,N-1):
      result[i,i] = T(1.0)
  when (T is Complex):
    for i in countup(0,N-1):
      result[i,i] = (1.0,0.0)


proc ones*[T](r,c: int): Matrix[T] =
  ## Returns a matrix with `r` rows of `c` columns which
  ## has all elements set to 1.
  result = newMatrix[T](r,c)
  when not (T is Complex):
    for i in countup(0,r*c-1):
      result.data[i] = T(1)
  when (T is Complex):
    for i in countup(0,N-1):
      result[i,i] = (1.0,0.0)


proc rand*[T](r,c: int, maxVal: float): Matrix[T] =
  ## Returns a matrix with `r` rows of `c` columns which
  ## has random elements between 0 and `maxVal`.
  result = newMatrix[T](r,c)
  when not (T is Complex):
    for i in countup(0,r*c-1):
      result.data[i] = T(random(maxVal))
  when (T is Complex):
    for i in countup(0,N-1):
      result[i,i] = (random(maxVal),0.0)



proc zeros*[T](r,c: int): Matrix[T] = newMatrix[T](r, c)
  ## Returns a matrix with `r` rows of `c` columns which
  ## has all elements set to 0.



proc setSize*[T](x: var Matrix[T], rows, cols: int) =
  ## Initializes the matrix by allocating memory
  ## for the data and setting the number of rows and columns.
  x.numRows = rows
  x.numCols = cols
  newSeq(x.data, rows*cols)




proc `[]`*[T](x: Matrix[T], r,c: int): T =
  ## Element access.  Returns the element at row `r` column `c`.
  if r<0  or  r>(x.rows()-1):
    raise newException(IndexError, "matrix index out of range")
  if c<0  or  c>(x.cols()-1):
    raise newException(IndexError, "matrix index out of range")

  result = x.data[x.index(r,c)]


proc `[]=`*[T](x: var Matrix[T], r,c: int, a: T) =
  ## Sets the value of the element at row `r` column `c` to
  ## the value supplied in `a`.
  if r<0  or  r>(x.rows()-1):
    raise newException(IndexError, "matrix index out of range")
  if c<0  or  c>(x.cols()-1):
    raise newException(IndexError, "matrix index out of range")

  x.data[x.index(r,c)] = a


proc `$`*[T](x: Matrix[T]): string =
  ## "Pretty prints" the matrix.  All elements of the matrix are
  ## included so large matrices will result in large strings.
  result = ""
  for r in countup(0,x.rows()-1):
    result = result & "|"
    for c in countup(0,x.cols()-1):
      result = result & $x[r,c]
      if c != (x.cols()-1):
        result = result & ", "
    result = result & "|\n"


proc transpose*[T](x: Matrix[T]): Matrix[T] =
  ## Transpose matrix
  result = newMatrix[T](x.cols, x.rows)
  for r in countup(0,x.rows()-1):
    for c in countup(0,x.cols()-1):
      result.data[ result.index(c,r) ] = x.data[ x.index(r,c) ]


proc `==`*[T](a: Matrix[T], b: Matrix[T]): bool =
  ## Compare two matrices `x` and `y` for equality.  The
  ## matrices must be of the same size and all elements
  ## must have the same value or the return result will
  ## be `false`.
  if (a.rows()==b.rows()) and (a.cols()==b.cols()):
    result = true
    for i in low(a.data)..high(a.data):
      if a.data[i] != b.data[i]:
        result = false
        break
  else:
    result=false



proc `=~`*[T](a, b: Matrix[T]): bool =
  ## Compare two matrices `a` and `b` approximately.
  if (a.rows()==b.rows()) and (a.cols()==b.cols()):
    result = true
    for i in low(a.data)..high(a.data):
      if abs(a.data[i] - b.data[i]) > EPS:
        result = false
        break
  else:
    result=false



proc `*`*[T](a: Matrix[T], b: Matrix[T]): Matrix[T] =
  ## Matrix multiply
  assert( a.cols()==b.rows() )
  result.setSize(a.rows(), b.cols())
  for i in countup(0,a.rows()-1):
    for j in countup(0,b.cols()-1):
      result[i,j] = 0.0
      for k in countup(0,a.cols()-1):
        result.data[ result.index(i,j) ] =
              result.data[ result.index(i,j) ] +
              a.data[ a.index(i,k) ]*b.data[ b.index(k,j) ]


proc `*.`*[T](a: Matrix[T], b: Matrix[T]): Matrix[T] =
  ## Element-by-element multiply
  assert( a.rows()==b.rows() )
  assert( a.cols()==b.cols() )
  result.setSize(a.rows(), a.cols())
  for r in countup(0,a.rows()-1):
    for c in countup(0,a.cols()-1):
      result.data[ result.index(r,c) ] =
         a.data[ a.index(r,c) ]*b.data[ b.index(r,c) ]

proc `*`*[T](a: Matrix[T], b: T): Matrix[T] =
  ## Multiplies a matrix `a` by type T `b` (a*b)
  result.setSize(a.rows(),a.cols())
  for i in low(result.data)..high(result.data):
    result.data[i] = a.data[i] * b

proc `*`*[T](a: T, b: Matrix[T]): Matrix[T] =
  ## Multiplies a matrix `b` by type T `a` (a*b)
  result.setSize(b.rows(),b.cols())
  for i in low(result.data)..high(result.data):
    result.data[i] = a * b.data[i]


proc `-`*[T](a: Matrix[T], b: Matrix[T]): Matrix[T] =
  ## Element-by-element subtraction (a-b)
  assert( a.rows()==b.rows() )
  assert( a.cols()==b.cols() )
  result.setSize(a.rows(),a.cols())
  for i in low(result.data)..high(result.data):
    result.data[i] = a.data[i] - b.data[i]

proc `-`*[T](a: Matrix[T], b: T): Matrix[T] =
  ## Subtraction of type T `b` from matrix `a` (a-b)
  result.setSize(a.rows(),a.cols())
  for i in low(result.data)..high(result.data):
    result.data[i] = a.data[i] - b

proc `-`*[T](a: T, b: Matrix[T]): Matrix[T] =
  ## Subtraction of matrix `a` from type T `b` (a-b)
  result.setSize(b.rows(), b.cols())
  for i in low(result.data)..high(result.data):
    result.data[i] = a-b.data[i]


proc `+`*[T](a: Matrix[T], b: Matrix[T]): Matrix[T] =
  ## Element-by-element addition (a+b)
  assert( a.rows()==b.rows() )
  assert( a.cols()==b.cols() )
  result.setSize(a.rows(),a.cols())
  for i in low(result.data)..high(result.data):
    result.data[i] = a.data[i] + b.data[i]

proc `+`*[T](a: Matrix[T], b: T): Matrix[T] =
  ## Addition of type T to matrix (a+b)
  result.setSize(a.rows(), a.cols())
  for i in low(result.data)..high(result.data):
    result.data[i] = a.data[i] + b

proc `+`*[T](a: T, b: Matrix[T]): Matrix[T] =
  ## Addition of type T to matrix (a+b)
  result.setSize(b.rows(), b.cols())
  for i in low(result.data)..high(result.data):
    result.data[i] = a + b.data[i]


proc `/`*[T](a: Matrix[T], b: T): Matrix[T] =
  ## Division of matrix `a` by type T `b` (a/b)
  result.setSize(a.rows(), a.cols())
  for i in low(result.data)..high(result.data):
    result.data[i] = a.data[i] / b


# Apply a math function to each element of the matrix.
template apply (a, b: typed): untyped =
  result.setSize(a.rows(), a.cols())
  for i in countup(0,len(a.data)-1):
    result.data[i] = b(a.data[i])

proc sin*[T](a: Matrix[T]):    Matrix[T] = apply(a, sin)
proc cos*[T](a: Matrix[T]):    Matrix[T] = apply(a, cos)
proc tan*[T](a: Matrix[T]):    Matrix[T] = apply(a, tan)
proc arcsin*[T](a: Matrix[T]): Matrix[T] = apply(a, arcsin)
proc arccos*[T](a: Matrix[T]): Matrix[T] = apply(a, arccos)
proc arctan*[T](a: Matrix[T]): Matrix[T] = apply(a, arctan)
proc sinh*[T](a: Matrix[T]):   Matrix[T] = apply(a, sinh)
proc cosh*[T](a: Matrix[T]):   Matrix[T] = apply(a, cosh)
proc tanh*[T](a: Matrix[T]):   Matrix[T] = apply(a, tanh)
proc sqrt*[T](a: Matrix[T]):   Matrix[T] = apply(a, sqrt)
proc ln*[T](a: Matrix[T]):     Matrix[T] = apply(a, ln)
proc log10*[T](a: Matrix[T]):  Matrix[T] = apply(a, log10)
proc log2*[T](a: Matrix[T]):   Matrix[T] = apply(a, log2)
proc exp*[T](a: Matrix[T]):    Matrix[T] = apply(a, exp)
proc abs*[T](a: Matrix[T]):    Matrix[T] = apply(a, abs)




proc qr*[T](AM: Matrix[T], Q: var Matrix[T], R: var Matrix[T]) =
  var A = AM
  var n = A.rows()
  var m = A.cols()

  var d: Matrix[T] = newMatrix[T](n,m)
  for j in countup(0,n-1):
    var s = T(0)
    for i in countup(j,m-1):
      s = s + A[i,j]*A[i,j]

    s = sqrt(s)
    if A[j,j]>T(0):
      d[j,0] = -s
    else:
      d[j,0] = s

    var fak = sqrt(s * (s+abs(A[j,j])))
    A[j,j] = A[j,j] - d[j,0]
    for k in countup(j,m-1):
      A[k,j] = A[k,j]/fak

    for i in countup((j+1),n-1):
      s = T(0)
      for k in countup(j,m-1):
        s = s + A[k,j]*A[k,i]

      for k in countup(j,m-1):
        A[k,i] = A[k,i] - A[k,j]*s

  # Reconstruct Q and R.  The diagonals of R are in d, the upper
  # triangle of R is in the upper triangle of A.  The lower
  # triangle of A holds the Householder vectors which we use to
  # rebuild Q (via Q transpose).
  var Qt: Matrix[T] = eye[T](n)
  for i in countup(0,min(n,m)-1):
    R[i,i] = d[i,0]

    var w: Matrix[T] = newMatrix[T](n,1)

    for k in countup(i,n-1):
      w[k,0] = A[k,i]

    Qt = (eye[T](n) - w*transpose(w))*Qt
    for j in countup((i+1),n-1):
      R[i,j] = A[i,j]

  Q = transpose(Qt)
