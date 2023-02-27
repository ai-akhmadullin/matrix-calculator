from __future__ import annotations
from typing import List, Tuple

operations = ["Addition",
              "Subtraction",
              "Multiplication",
              "Multiplication by a scalar",
              "Transposition",
              "Matrix inverse",
              "Determinant",
              "LU decomposition",
              "PLU decomposition",
              "RREF"]


class Matrix:
    """A class used to represent matrices, 'r' for the number of rows, 'c' for the number of columns,
    'val' will contain actual data of the matrix."""

    def __init__(self, val: List[List[float]]) -> None:
        self.val = val
        self.r = len(self.val)
        self.c = len(self.val[0])

    def print(self) -> None:
        """Print the matrix in a square form."""
        print("[", end="")
        for i in range(self.r):
            if i == self.r - 1:
                print(self.val[i], end=']')
            else:
                print(self.val[i])
        print()

    def add(self, other: Matrix) -> Matrix:
        """Add two matrices."""
        if self.r != other.r or self.c != other.c:
            raise ValueError("These matrices cannot be added. You can try again or change the operation.")

        result = Matrix([[0] * self.c for _ in range(self.r)])
        for i in range(self.r):
            for j in range(self.c):
                result.val[i][j] = self.val[i][j] + other.val[i][j]
        return result

    def subtract(self, other: Matrix) -> Matrix:
        """Subtract two matrices."""
        if self.r != other.r or self.c != other.c:
            raise ValueError("These matrices cannot be subtracted. You can try again or change the operation.")

        result = Matrix([[0] * self.c for _ in range(self.r)])
        for i in range(self.r):
            for j in range(self.c):
                result.val[i][j] = self.val[i][j] - other.val[i][j]
        return result

    def scalar(self, k: float) -> Matrix:
        """Multiply a matrix by a scalar."""
        for i in range(self.r):
            for j in range(self.c):
                self.val[i][j] *= k
        return self

    def multiply(self, other: Matrix) -> Matrix:
        """Multiply two matrices."""
        if self.c != other.r:
            raise ValueError("These matrices cannot be multiplied. You can try again or change the operation.")

        result = Matrix([[0] * other.c for _ in range(self.r)])
        for i in range(result.r):
            for j in range(result.c):
                for k in range(self.c):
                    result.val[i][j] += (self.val[i][k] * other.val[k][j])
        return result

    def transpose(self) -> Matrix:
        """Transpose a matrix."""
        temp = Matrix([[0] * self.r for _ in range(self.c)])
        for i in range(temp.r):
            for j in range(temp.c):
                temp.val[i][j] = self.val[j][i]
        self.val, self.r, self.c = temp.val, temp.r, temp.c
        return self

    def rref(self) -> Matrix:
        """Find a reduced row echelon form of a matrix."""
        n, m = 0, 0  # the indexes of the pivot we will be working with, 'n' for a row and 'm' for a column.

        # Stop if all the elements to the right and below equal 0, otherwise find the pivot,
        # leaving behind zero sub-columns.
        resume = True
        while resume:
            resume = False
            k = n
            for j in range(m, self.c):
                for i in range(n, self.r):
                    if self.val[i][j] != 0:
                        k = i
                        m = j
                        resume = True
                        break
                if resume:
                    break
            if not resume:
                return self
            self.val[k], self.val[n] = self.val[n], self.val[k]

            # Divide the row by the value of the pivot, making the pivot equal 1.
            a = self.val[n][m]
            for j in range(self.c):
                self.val[n][j] /= a

            # Subtract the row with the current pivot from other rows above and below it,
            # multiplying it by the value of the element above/below the pivot, thus making
            # all these elements above/below the pivot equal 0.
            for i in range(self.r):
                if i != n:
                    b = self.val[i][m]
                    for j in range(self.c):
                        self.val[i][j] -= (b * self.val[n][j])

            # Check whether we reached the right or bottom edge and update the indexes of
            # the next potential pivot if not.
            if n == (self.r - 1) or m == (self.c - 1):
                return self
            n += 1
            m += 1

    def invert(self) -> Matrix:
        """Find an inverse of a matrix."""
        if self.r != self.c:
            raise ValueError("An inverse matrix exists only for a square matrix. You can try again or change "
                             "the operation.")

        # If we add an identity matrix to the right side of our matrix and apply RREF, we will get an
        # identity matrix on the left side and the inverse of our matrix on the right side.
        for i in range(self.r):
            a = [0] * (self.c - 1)
            a.insert(i, 1)
            self.val[i].extend(a)
        self.c *= 2
        self.rref()
        self.c //= 2

        # If we didn't get an identity matrix on the left side, it means that our matrix is not regular and
        # does not have an inverse. If there is an identity matrix on the left, we want to delete it and return
        # only our inverse.
        for i in range(self.r):
            for j in range(self.c):
                if (i == j and self.val[i][j] != 1) or (i != j and self.val[i][j] != 0):
                    raise ValueError("The matrix does not have an inverse matrix. You can try again or "
                                     "change the operation.")
            del self.val[i][0:self.c]
        return self

    def lu(self) -> Tuple[Matrix, Matrix]:
        """LU decomposition of a matrix."""
        if self.r != self.c:
            raise ValueError("LU decomposition exists only for a square matrix. You can try again or "
                             "change the operation.")

        lower = Matrix([[0] * self.c for _ in range(self.r)])
        upper = Matrix([[0] * self.c for _ in range(self.r)])
        for i in range(self.r):

            # We fill the upper matrix from left to right. 's' will represent the value by which our element has
            # changed during previous REF transformations.
            for j in range(i, self.c):
                s = 0
                for k in range(i):
                    s += (lower.val[i][k] * upper.val[k][j])
                upper.val[i][j] = self.val[i][j] - s

            # We fill the lower matrix from top to bottom, it is filled with -factors by which the pivoting rows were
            # multiplied before subtracting them from the processed row.
            for j in range(i, self.r):
                if i == j:
                    lower.val[i][i] = 1
                else:
                    s = 0
                    for k in range(i):
                        s += (lower.val[j][k] * upper.val[k][i])
                    lower.val[j][i] = self.val[j][i] - s / upper.val[i][i]

        return lower, upper

    def plu(self) -> Tuple[Matrix, Matrix, Matrix]:
        """PLU decomposition of a matrix."""
        if self.r != self.c:
            raise ValueError("PLU decomposition exists only for a square matrix. You can try again or "
                             "change the operation.")

        # We will start with an identity matrix as our permutation matrix.
        perm = Matrix([[0] * self.c for _ in range(self.r)])
        for i in range(self.r):
            perm.val[i][i] = 1

        # We want to find a pivot which will not equal 0 (and find the maximum of such), and interchange rows in our
        # permutation matrix. However, if all the elements in a column equal 0, then our matrix is singular and PLU
        # decomposition does not exist.
        for k in range(self.r - 1):
            p = 0
            k_ = k  # k_ is the index of the row that will be interchanged.
            for i in range(k, self.r):
                if abs(self.val[i][k]) > p:
                    p = self.val[i][k]
                    k_ = i
            if p == 0:
                raise ValueError("PLU decomposition exists only for a regular matrix. You can try again or "
                                 "change the operation.")
            perm.val[k], perm.val[k_] = perm.val[k_], perm.val[k]
            self.val[k], self.val[k_] = self.val[k_], self.val[k]

            # We will prepare our matrix for the decomposition, dividing the rows by the pivot and
            # subtracting the pivoting row from the other rows.
            for i in range(k + 1, self.r):
                self.val[i][k] /= self.val[k][k]
                for j in range(k + 1, self.c):
                    self.val[i][j] -= (self.val[i][k] * self.val[k][j])

        perm.invert()
        lower = Matrix([[0] * self.c for _ in range(self.r)])
        upper = Matrix([[0] * self.c for _ in range(self.r)])

        # We are ready to form the lower and the upper matrices.
        for i in range(self.r):
            for j in range(self.c):
                if i == j:
                    lower.val[i][j] = 1
                    upper.val[i][j] = self.val[i][j]
                elif i > j:
                    lower.val[i][j] = self.val[i][j]
                elif i < j:
                    upper.val[i][j] = self.val[i][j]

        return perm, lower, upper

    def determinant(self) -> float:
        """Find a determinant of a matrix."""
        if self.r != self.c:
            raise ValueError("A determinant can be calculated only for a square matrix. You can try again or "
                             "change the operation.")

        perm, lower, upper = self.plu()
        det = 1

        # The determinant of a permutation matrix is either 1 or -1, depending on the number of interchanged rows.
        sign = 1
        for i in range(self.r):
            if perm.val[i][i] != 1:
                sign *= -1
                k = perm.val[i].index(1)
                perm.val[i], perm.val[k] = perm.val[k], perm.val[i]

        # det(A) = det(PLU) = det(P)det(L)det(U). The determinant of a triangular matrix is
        # the product of elements on the diagonal.
        for i in range(self.r):
            det *= (lower.val[i][i] * upper.val[i][i])
        return sign * det


def addition(tensor: List[List[List[float]]]) -> Matrix:
    """Take a list of matrices, add them, and return the resulting matrix."""
    n = len(tensor)
    _matrix = Matrix(tensor[0])
    for i in range(1, n):
        _matrix = _matrix.add(Matrix(tensor[i]))
    return _matrix


def subtraction(tensor: List[List[List[float]]]) -> Matrix:
    """Take a list of matrices, subtract them, and return the resulting matrix."""
    n = len(tensor)
    _matrix = Matrix(tensor[0])
    for i in range(1, n):
        _matrix = _matrix.subtract(Matrix(tensor[i]))
    return _matrix


def multiply(tensor: List[List[List[float]]]) -> Matrix:
    """Take a list of matrices, multiply them, and return the resulting matrix."""
    n = len(tensor)
    _matrix = Matrix(tensor[0])
    for i in range(1, n):
        _matrix = _matrix.multiply(Matrix(tensor[i]))
    return _matrix


def scalar(tensor: List[List[List[float]]]) -> Matrix:
    """Take a list of matrices with one matrix, multiply it by a scalar, and return the resulting matrix."""
    k = float(input("Enter the scalar -> "))
    print()
    _matrix = Matrix(tensor[0])
    return _matrix.scalar(k)


def transpose(tensor: List[List[List[float]]]) -> Matrix:
    """Take a list of matrices with one matrix, transpose it, and return the resulting matrix."""
    _matrix = Matrix(tensor[0])
    return _matrix.transpose()


def rref(tensor: List[List[List[float]]]) -> Matrix:
    """Take a list of matrices with one matrix, perform RREF, and return the resulting matrix."""
    _matrix = Matrix(tensor[0])
    return _matrix.rref()


def invert(tensor: List[List[List[float]]]) -> Matrix:
    """Take a list of matrices with one matrix and return the the inverse of the matrix."""
    _matrix = Matrix(tensor[0])
    return _matrix.invert()


def lu(tensor: List[List[List[float]]]) -> Tuple[Matrix, Matrix]:
    """Take a list of matrices with one matrix, apply LU decomposition, and return the resulting matrices."""
    _matrix = Matrix(tensor[0])
    return _matrix.lu()


def plu(tensor: List[List[List[float]]]) -> Tuple[Matrix, Matrix, Matrix]:
    """Take a list of matrices with one matrix, apply PLU decomposition, and return the resulting matrices."""
    _matrix = Matrix(tensor[0])
    return _matrix.plu()


def determinant(tensor: List[List[List[float]]]) -> int:
    """Take a list of matrices with one matrix and return its determinant."""
    _matrix = Matrix(tensor[0])
    return _matrix.determinant()
