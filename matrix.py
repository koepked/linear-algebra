# Copyright 2014 Dan Koepke koepked@gmail.com
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

class Matrix(object):
    """Perform the mathematical functions of a martix.
    """
    @classmethod
    def identity(cls, x):
        """Return an x-by-x identity matrix.

        Args:
            x (int): Square dimension of the returned identity matrix.
        """
        M = [[0 for y in range(x)] for z in range(x)]
        for i in range(x):
            M[i][i] = 1
        return Matrix(M)

    @classmethod
    def zero(cls, rows, cols):
        """Return a matrix of given dimension populated with zeros.

        Args:
            rows (int): Number of rows in returned matrix.
            cols (int): Number of columns in returned matrix.
        """
        M = [[0 for y in range(cols)] for z in range(rows)]
        return Matrix(M)

    def __add__(self, y):
        """x.__add__(y) <==> x + y

        Raises:
            ValueError if x + y is undefined.
        """
        if len(self._data) != len(y._data) or len(self._data[0]) != len(y._data[0]):
            raise ValueError("Improper dimensions for operation.")

        return Matrix( \
            [[w+x for w,x in zip(y,z)] for y,z in zip(self._data, y._data)])

    def __eq__(self, y):
        """x.__eq__(y) <==> x == y
        """
        return self._data == y._data

    def __init__(self, values=[]):
        """Populate a Matrix instance with the given values.

        Args:
            values (list of list of int): List of values to be entered into the
                matrix. List format is [[a00, a01, ...],[a10, a11, ...], ...]
        """
        for i in range(1, len(values)):
            if len(values[i]) != len(values[0]):
                raise ValueError("")

        self._data = tuple([tuple(x) for x in values])


    def __mul__(self, y):
        """x.__mul__(y) <==> x * y

        Raises:
            ValueError if x * y is undefined.
        """
        if len(self._data[0]) != len(y._data):
            raise ValueError("Improper dimensions for operation.")

        new = [[None for x in range(len(y._data[0]))] for z in range(len(self._data))]

        # Prolly can do the following as list comp. and remove line above this.
        for i in range(0, len(self._data)):
            for j in range(0, len(y._data[0])):
                new[i][j] = sum([ self._data[i][x] * y._data[x][j] for x in \
                    range(0, len(y._data))])

        return Matrix(new)

    def __pow__(self, x):
        """x.__pow__(y) <==> x ** y

        Raises:
            ValueError if x ** y is undefined.
        """
        if len(self._data) != len(self._data[0]):
            raise ValueError("Improper dimensions for operation.")

        if x == -1:
            M = self.augment(Matrix.identity(len(self._data))).rref()
            N = M.get_col(len(self._data))
            for i in range(len(self._data) + 1, 2 * len(self._data)):
                N = N.augment(M.get_col(i))
            return N

        if x < -1:
            M = self ** -1
            return M**(-x)

        M = self.copy()
        for i in range(x-1):
            M *= self

        return M

    def __repr__(self):
        """Return a string representation of the matrix.
        """
        return self._data.__repr__()

    def __sub__(self, x):
        """x.__sub__(y) <==> x - y

        Raises:
            ValueError if x - y is undefined.
        """
        return self + x.mult_by_scalar(-1)

    def _sub_matrix(self, nixedRow, nixedCol):
        # Return the matrix given by removing nixedRow and nixedCol from this
        # matrix.
        if nixedRow >= len(self._data) or nixedCol >= len(self._data[0]):
            raise ValueError("Improper dimensions for operation.")

        l= [x[:nixedCol] + x[nixedCol+1:] for x in self._data[:nixedRow]
                        + self._data[nixedRow+1:]]
        return Matrix(l)

    def augment(self, M):
        """Return and augmented matrix of the form [self|M].

        Args:
            M (Matrix): Matrix self will be augmented with.

        Raises:
            ValueError if self.num_rows() != M.num_rows()
        """
        if len(self._data) != len(M._data):
            raise ValueError("Improper dimensions for operation.")

        m1 = [[x for x in y] for y in self._data]
        m2 = [[x for x in y] for y in M._data]

        return Matrix([x+y for x,y in zip(m1, m2)])

    def copy(self):
        """Return a deep copy of the matrix.
        """
        return Matrix(list(self._data))

    def det(self):
        """Return the determinant of the matrix.

        Raises:
            ValueError if the determinant for the matrix is undefined.
        """
        if len(self._data) != len(self._data[0]):
            raise ValueError("Improper dimensions for operation.")

        if len(self._data) == 2 and len(self._data[0]) == 2:
            return self._data[0][0] * self._data[1][1] \
                    - self._data[0][1] * self._data[1][0]

        total = 0
        for i in range(len(self._data)):
            total += self._data[i][0] * self._sub_matrix(i, 0).det() * (-1)**i

        return total

    def get_col(self, col):
        """Return a matrix consisting of the given column.

        Args:
            col (int): The column number to return.
        """
        return Matrix([[y[col]] for y in self._data])

    def get_row(self, row):
        """Return a matrix consisting of the given row.

        Args:
            row (int): The row number to return.
        """
        return Matrix([self._data[row]])

    def mult_by_scalar(self, x):
        """Return the given scalar multiple of the matrix.

        Args:
            x (number): Multiply martix by x
        """
        return Matrix([[x * y for y in z] for z in self._data])

    def rref(self):
        """Return the reduced row echelon format of the matrix.
        """
        M = list(self._data)
        rowCount = len(M)
        colCount = len(M[0])
        lead = 0
    
        # For each row
        for r in range(rowCount):
            if lead >= colCount:
                return Matrix(M)

            i = r
    
            # Find the next row that has leading nonzero val (left-to-right)
            while M[i][lead] == 0:
                i += 1
                if i >= rowCount:
                    i = r
                    lead += 1
                    if lead >= colCount:
                        return Matrix(M)
    
            # Swap row found above with current row
            M[i], M[r] = M[r], M[i]
    
            # Nomralize row (leading 1)
            lv = M[r][lead]
            M[r] = [float(Mval) / lv for Mval in M[r]]
    
            # Subtract multiple of current row from all other rows so other rows
            # have zeroes in column of current row's leading 1
            for i in range(rowCount):
                if i != r:
                    lv = M[i][lead]
                    M[i] = [iv - lv*rv for iv,rv in zip(M[i], M[r])]
    
            lead += 1
        return Matrix(M)

    def transpose(self):
        """Return the transpose of this matrix.
        """
        l = [[row[x] for row in self._data] for x in range(len(self._data[0]))]
        return Matrix(l)


if __name__ == "__main__":
    """Perform unit tests.
    """
    # Initialization one row
    t = ((1,2,3),)
    l = [[x for x in y] for y in t]
    m = Matrix(l)
    assert m._data == t

    # Initialization one col
    t = ((1,),(2,),(3,))
    l = [[x for x in y] for y in t]
    n = Matrix(l)
    assert n._data == t

    # Initialization square
    t = ((1,2,3),(4,5,6),(7,8,9))
    l = [[x for x in y] for y in t]
    o = Matrix(l)
    assert o._data == t

    # Initialization numRows > numCols
    t = ((1,2),(4,5),(7,8))
    l = [[x for x in y] for y in t]
    p = Matrix(l)
    assert p._data == t

    # Initialization numCols > numRows
    t = ((1,2,3),(4,5,6))
    l = [[x for x in y] for y in t]
    q = Matrix(l)
    assert q._data == t

    # Initialization list arg has varying row lengths
    t = ((1,2,3),(4,5))
    l = [[x for x in y] for y in t]
    proper = False
    try:
        r = Matrix(l)
    except ValueError:
        proper = True
    assert proper

    # Zero
    assert Matrix.zero(1,1) == Matrix([[0]])
    assert Matrix.zero(3,3) == Matrix([[0,0,0],[0,0,0],[0,0,0]])
    assert Matrix.zero(2,5) == Matrix([[0,0,0,0,0],[0,0,0,0,0]])

    # Add/scalar multiply all of the above and compare
    assert m+m == m.mult_by_scalar(2)
    assert n+n == n.mult_by_scalar(2)
    assert o+o == o.mult_by_scalar(2)
    assert p+p == p.mult_by_scalar(2)
    assert q+q == q.mult_by_scalar(2)
    assert m+m+m == m.mult_by_scalar(3)
    assert n+n+n == n.mult_by_scalar(3)
    assert o+o+o == o.mult_by_scalar(3)
    assert p+p+p == p.mult_by_scalar(3)
    assert q+q+q == q.mult_by_scalar(3)
    
    # Add matrices w/ improper dimensions
    proper = False
    try:
        m+n
    except ValueError:
        proper = True
    assert proper

    # Subtract
    assert m-m == Matrix.zero(1,3)
    assert n-n == Matrix.zero(3,1)
    assert o-o == Matrix.zero(3,3)
    assert p-p == Matrix.zero(3,2)
    assert q-q == Matrix.zero(2,3)
    assert m-m-m == m.mult_by_scalar(-1)
    assert n-n-n == n.mult_by_scalar(-1)
    assert o-o-o == o.mult_by_scalar(-1)
    assert p-p-p == p.mult_by_scalar(-1)
    assert q-q-q == q.mult_by_scalar(-1)

    # Subtract matrices w/ improper dimensions
    proper = False
    try:
        m-n
    except ValueError:
        proper = True
    assert proper

    # Initalization empty list
    l = []
    m = Matrix(l)
    assert m._data == ()

    # Equality check
    t = ((1,2,3),(4,5,6),(7,8,9))
    l1 = [[x for x in y] for y in t]
    l2 = [[x for x in y] for y in t]
    m1 = Matrix(l1)
    m2 = Matrix(l2)
    t = ((1,2,3),(4,1,6),(7,8,9))
    l3 = [[x for x in y] for y in t]
    m3 = Matrix(l3)
    assert m1 == m2
    assert m1 != m3
    assert m2 != m3

    # Identity
    t = ((1,0,0),(0,1,0),(0,0,1))
    l = [[x for x in y] for y in t]
    m = Matrix(l)
    assert m == Matrix.identity(3)
    t = ((1,0),(0,1))
    l = [[x for x in y] for y in t]
    m = Matrix(l)
    assert m == Matrix.identity(2)

    # getCol
    m = Matrix([[1,2,3],[4,5,6],[7,8,9]])
    assert m.get_col(0) == Matrix([[1],[4],[7]])
    assert m.get_col(1) == Matrix([[2],[5],[8]])
    assert m.get_col(2) == Matrix([[3],[6],[9]])

    # getRow
    m = Matrix([[1,2,3],[4,5,6],[7,8,9]])
    assert m.get_row(0) == Matrix([[1,2,3]])
    assert m.get_row(1) == Matrix([[4,5,6]])
    assert m.get_row(2) == Matrix([[7,8,9]])

    # copy
    m = Matrix([[1,2,3],[4,5,6],[7,8,9]])
    m_copy = m.copy()
    assert m == m_copy
    l = [[x for x in y] for y in m._data]
    l[1][1] = 99
    m._data = (tuple(l[0]), tuple(l[1]), tuple(l[2]))
    assert m != m_copy

    # _sub_matrix square
    m = Matrix([[1,2,3],[4,5,6],[7,8,9]])
    assert m._sub_matrix(0,0) == Matrix([[5,6],[8,9]])
    assert m._sub_matrix(0,1) == Matrix([[4,6],[7,9]])
    assert m._sub_matrix(0,2) == Matrix([[4,5],[7,8]])
    assert m._sub_matrix(1,0) == Matrix([[2,3],[8,9]])
    assert m._sub_matrix(1,1) == Matrix([[1,3],[7,9]])
    assert m._sub_matrix(1,2) == Matrix([[1,2],[7,8]])
    assert m._sub_matrix(2,0) == Matrix([[2,3],[5,6]])
    assert m._sub_matrix(2,1) == Matrix([[1,3],[4,6]])
    assert m._sub_matrix(2,2) == Matrix([[1,2],[4,5]])

    # _sub_matrix rectangular
    m = Matrix([[1,2,3],[4,5,6]])
    assert m._sub_matrix(0,0) == Matrix([[5,6]])
    assert m._sub_matrix(0,1) == Matrix([[4,6]])
    assert m._sub_matrix(0,2) == Matrix([[4,5]])
    assert m._sub_matrix(1,0) == Matrix([[2,3]])
    assert m._sub_matrix(1,1) == Matrix([[1,3]])
    assert m._sub_matrix(1,2) == Matrix([[1,2]])

    # determinant square
    assert Matrix([[3,4,5],[1,1,4],[8,2,3]]).det() == 71
    assert Matrix([[-3,-4,-5],[1,1,4],[8,2,3]]).det() == -71
    assert Matrix([[1,2,3],[2,4,6],[1,1,1]]).det() == 0
    # determinant nonsquare
    proper = False
    try:
        Matrix([[1,2,3],[1,2,3]]).det()
    except ValueError:
        proper = True
    assert proper

    # rref square invertible
    t = ((1,2,3),(4,5,6),(7,8,9))
    l = [[x for x in y] for y in t]
    m = Matrix(l)
    # Have to check like this for now due to floating point implementation:
    assert abs((m.rref() - Matrix.identity(3)).det()) < .00000000001

    # rref square noninvertible
    m = Matrix([[1,2,3],[2,4,6],[1,1,1]])
    assert m.rref() == Matrix([[1,0,-1],[0,1,2],[0,0,0]])

    # rref rectangular
    m = Matrix([[1,2,3],[1,1,1]])
    assert m.rref() == Matrix([[1,0,-1],[0,1,2]])

    # Multipy square and rectangular matrices
    m1 = Matrix([[3,4,5],[1,1,4],[8,2,3]])
    m2 = Matrix([[1,2,3],[2,4,6],[3,1,1]])
    assert m1*m2 == Matrix([[26,27,38],[15,10,13],[21,27,39]])
    assert m2*m1 == Matrix([[29,12,22],[58,24,44],[18,15,22]])
    m1 = Matrix([[1,2,3]])
    m2 = Matrix([[1,2],[2,4],[3,1]])
    assert m1*m2 == Matrix([[14, 13]])

    # Multiply matricies w/ improper dimensions
    proper = False
    try:
        m2*m1
    except ValueError:
        proper = True
    assert proper

    # Inverse of square multiplied equals proper identity
    #
    # Won't work with floating point implementation
    #
    #assert m1 * (m1 ** -1) == Matrix.identity(3)
    #assert (m1 ** -1) * m1 == Matrix.identity(3)

    # Inverse of nonsquare
    proper = False
    try:
        Matrix([[1,2,3],[4,5,6]]) ** -1
    except ValueError:
        proper = True
    assert proper

    # Transpose square
    assert Matrix([[1,2,3],[4,5,6],[7,8,9]]).transpose() \
            == Matrix([[1,4,7],[2,5,8],[3,6,9]])

    # Transpose nonsquare
    assert Matrix([[1,2,3],[4,5,6]]).transpose() \
            == Matrix([[1,4],[2,5],[3,6]])

    # Augment
    assert Matrix([[1],[2],[3]]).augment(Matrix([[4],[5],[6]])) \
            == Matrix([[1,4],[2,5],[3,6]])

    # Augment improper dimensions
    proper = False
    try:
        Matrix([[1],[2]]).augment(Matrix([[1]]))
    except ValueError:
        proper = True
    assert proper
