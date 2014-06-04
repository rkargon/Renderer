import java.awt.Dimension;
import java.util.Arrays;

/**
 * Matrix class. Mainly used for multiplying matrices while rotating Cameras.
 * This was copied from a larger, more complete Matrix library, so there lots of
 * extra utiltity functions that may not be used at the moment, but will
 * probably come in handy in the future.
 * 
 * @author raphaelkargon
 * 
 */
public class Matrix {
	public static double DEFAULT_EPSILON = 1E-14;

	private int m, n; //rows, cols
	private double[] data; //stored as m blocks of n numbers, row 1 | row 2 | ... | row m

	public Matrix(int m, int n) {
		if (m < 0 || n < 0)
			throw new IllegalArgumentException("Negative Matrix Size!");
		this.m = m;
		this.n = n;
		this.data = new double[m * n];
	}

	public Matrix(double[] data, int m, int n) {
		this(m, n);
		this.data = Arrays.copyOf(data, data.length);
	}

	public Matrix(double[][] A) {
		this(A.length, A[0].length);

		for (int i = 0; i < A.length; i++) {
			System.arraycopy(A[i], 0, this.data, i * this.n, A[i].length);
		}
	}

	public Matrix(Matrix A) {
		this(A.data, A.m, A.n);
	}

	/**
	 * Copies a matrix to a new size, either adding 0's or discarding data that
	 * does not fit
	 * 
	 * @param A
	 * @param m
	 * @param n
	 */
	public Matrix(Matrix A, int m, int n) {
		this(m, n);

		for (int i = 0; i < m; i++) {
			if (i >= A.m) return;
			System.arraycopy(A.data, i * A.n, data, i * n, Math.min(A.n, n));
		}
	}

	/**
	 * 
	 * @return a COPY of the data of the matrix
	 */
	public double[] data() {
		return Arrays.copyOf(data, data.length);
	}

	public double[][] data2D() {
		double[][] d = new double[m][n];
		for (int i = 0; i < m; i++) {
			d[i] = Arrays.copyOfRange(data, i * n, (i + 1) * n);
		}

		return d;
	}

	public int rows() {
		return m;
	}

	public int cols() {
		return n;
	}

	public Dimension bounds() {
		return new Dimension(m, n);
	}

	/**
	 * Returns an item from the matrix
	 * Uses array indices, ie from [0...n-1/m-1];
	 * 
	 * @param i
	 * @param j
	 * @return
	 */
	public double get(int i, int j) {
		if (i < 0 || j < 0 || i >= m || j >= n)
			throw new IllegalArgumentException("Indices are out of bounds.");
		return this.data[i * n + j];
	}

	public void set(int i, int j, double x) {
		if (i < 0 || j < 0 || i >= m || j >= n)
			throw new IllegalArgumentException("Indices are out of bounds.");
		this.data[i * n + j] = x;
	}

	/**
	 * No boundary checks, directly accesses data array. get(i,j) ~
	 * fastget(i*n+j)
	 * 
	 * @param i
	 *            The index of the data array
	 * @return The value at that index
	 */
	public double fastget(int i) {
		return data[i];
	}

	/**
	 * No boundary checks, directly modifies data array. set(i,j, x) ~
	 * fastet(i*n+j, x)
	 * 
	 * @param i
	 *            The index of the data array
	 * @return The value at that index
	 */
	public void fastset(int i, double x) {
		data[i] = x;
	}

	public Matrix getRow(int i) {
		if (i < 0 || i >= n)
			throw new IllegalArgumentException("Row index " + i
					+ " is out of bounds.");
		double[] row = new double[n];
		for (int j = 0; j < m; j++) {
			row[j] = data[i * n + j];
		}
		return new Matrix(row, 1, n);
	}

	public Matrix getColumn(int j) {
		if (j < 0 || j >= n)
			throw new IllegalArgumentException("Column index " + j
					+ " is out of bounds.");
		double[] col = new double[m];
		for (int i = 0; i < m; i++) {
			col[i] = data[i * n + j];
		}
		return new Matrix(col, m, 1);
	}

	/**
	 * Returns a segment of the current matrix
	 * 
	 * Array indices (ie start at 0)
	 * Indices are inclusive, ie item (i2, j2) is returned in new array
	 * 
	 * @param i1
	 *            Top row
	 * @param j1
	 *            Left row
	 * @param i2
	 *            Bottom row (inclusive)
	 * @param j2
	 *            Right row (inclusive)
	 * @return A matrix that is a subsegment of this matrix
	 */
	public Matrix segment(int i1, int j1, int i2, int j2) {
		if (i1 > i2 || j1 > j2 || i1 < 0 || j1 < 0 || i2 >= m || j2 >= n)
			throw new IllegalArgumentException("Indices are out of bounds.");

		Matrix sub = new Matrix(1 + i2 - i1, 1 + j2 - j1);
		for (int i = i1; i <= i2; i++) {
			System.arraycopy(data, i * n + j1, sub.data, (i - i1) * sub.n, 1
					+ j2 - j1);
		}

		return sub;
	}

	public void swapRows(int i, int k) {
		if (i == k) return;
		if (i < 0 || k < 0 || i > m || k > m)
			throw new IllegalArgumentException("Indices are out of bounds.");
		double[] i_tmp = new double[n];
		System.arraycopy(data, i * n, i_tmp, 0, n); //row i to i_tmp
		System.arraycopy(data, k * n, data, i * n, n);//row k to i
		System.arraycopy(i_tmp, 0, data, k * n, n);//i_tmp to row k
	}

	public Matrix removeColumns(int[] cols) {
		if (cols.length == 0) return this.clone();

		Matrix A = new Matrix(m, n - cols.length);
		Arrays.sort(cols);
		int colsindex = 0;
		for (int i = 0; i < m; i++) {
			colsindex = 0;
			for (int j = 0; j < n; j++) {
				if (cols[colsindex] == j) {
					//skip column if it is to be omitted
					colsindex++;
					continue;
				}
				A.data[i * A.n + j - colsindex] = data[i * n + j]; //-colsindex shifts items down in new matrix when a column is skipped
			}
		}

		return A;
	}

	public Matrix removeRows(int[] rows) {
		if (rows.length == 0) return this.clone();

		Matrix A = new Matrix(m - rows.length, n);
		Arrays.sort(rows);
		int rowsindex = 0;

		for (int i = 0; i < m; i++) {
			if (rows[rowsindex] == i) {
				rowsindex++;
				continue;
			}
			System.arraycopy(data, i * n, A.data, (i - rowsindex) * n, n);
		}

		return A;

	}

	public void scale(double a) {
		for (int i = 0; i < data.length; i++) {
			data[i] *= a;
		}
	}

	public Matrix transpose() {
		Matrix tr = new Matrix(n, m);

		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				tr.data[j * m + i] = data[i * n + j];
			}
		}

		return tr;
	}

	/**
	 * Matrix multiplication.
	 * Does this naively, O(n^3) for nxn square matrices
	 * 
	 * @param B
	 *            Second matrix to multiply
	 * @return If this=M(mxn), B=M(nxp), then it returns this*B = M(mxp)
	 */
	public Matrix product(Matrix B) {
		if (n != B.m)
			throw new IllegalArgumentException("Matrix has noncompatible number of rows. This=("
					+ m + ", " + n + "), B=(" + B.m + ", " + B.n + ")");
		if (m == 1 && n == 1 && B.n == 1)
			return new Matrix(new double[] { data[0] * B.data[0] }, 1, 1);
		if (m == 2 && n == 2 && B.n == 2) return product_2x2(B);
		int p = B.n, n = B.m; //avoid confusion with m's and n's.
		Matrix prod = new Matrix(m, p);

		double sum = 0;
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < p; j++) {
				sum = 0;
				for (int k = 0; k < n; k++) {
					sum += data[i * n + k] * B.data[k * p + j];
				}
				prod.data[i * p + j] = sum;
			}
		}

		return prod;
	}

	public Matrix product_2x2(Matrix B) {
		if (B.m != 2 || B.n != 2 || m != 2 || n != 2)
			throw new IllegalArgumentException("Not 2x2 matrices.");
		double[] proddata = { data[0] * B.data[0] + data[1] * B.data[2],
				data[0] * B.data[1] + data[1] * B.data[3],
				data[2] * B.data[0] + data[3] * B.data[2],
				data[2] * B.data[1] + data[3] * B.data[3] };
		Matrix prod = new Matrix(proddata, 2, 2);
		return prod;
	}

	/**
	 * Fixes rounding errors when manipulating matrices of whole numbers
	 * represented as <code>doubles</code>.
	 * eg. 7.2875825470845E-15 is replaced by 0, 2.9999999999745 is replaced by
	 * 3, etc. <code>EPSILON</code> is used as an error threshold, differences
	 * above this value will be regarded as separate values and not rounding
	 * errors
	 */
	public void fixRoundingErrors(double EPSILON) {
		for (int i = 0; i < data.length; i++) {
			if (Math.abs(data[i] - Math.round(data[i])) < EPSILON)
				data[i] = Math.round(data[i]);
		}
	}

	public void fixRoundingErrors() {
		fixRoundingErrors(DEFAULT_EPSILON);
	}

	@Override
	public Matrix clone() {
		return new Matrix(data, m, n);
	}

	/**
	 * Comapres two matrices, using a given error range
	 */
	public boolean equals(Matrix M, double EPSILON) {
		if (this.m != M.m || this.n != M.n) return false;
		for (int i = 0; i < data.length; i++) {
			if (Math.abs(data[i] - M.data[i]) > EPSILON) return false;
		}
		return true;
	}

	@Override
	public boolean equals(Object o) {
		if (o instanceof Matrix) {
			Matrix M = (Matrix) o;
			if (this.m != M.m || this.n != M.n) return false;
			return Arrays.equals(data, M.data);
		}

		return false;
	}

	//Jenkin's hash function
	public int hashCode() {
		int result = 1;

		result = 31 * result + m;
		result = 31 * result + n;
		result = 31 * result + Arrays.hashCode(data);
		return result;
	}

	@Override
	public String toString() {
		return toString(true);
	}

	public String toString(boolean multiline) {
		if (n == 0 || m == 0) return "[ ]";
		String s = "";
		for (int i = 0; i < m; i++) {
			s += Arrays.toString(Arrays.copyOfRange(data, i * n, i * n + n));

			if (i < m - 1) s += ", ";
			if (multiline) s += "\n";
		}

		return s;
	}

	/**
	 * Returns an identity matrix of size n
	 * 
	 * @param n
	 *            The number of rows and columns in the matrix
	 */
	public static Matrix identityMatrix(int n) {
		Matrix I = new Matrix(n, n);
		for (int i = 0; i < n; i++) {
			I.data[i * n + i] = 1;
		}
		return I;
	}

	/**
	 * Creates a block matrix by combining A (mxn) and B (mxp) to form matrix
	 * |A B| (mx(n+p))
	 * 
	 * @param B
	 * @return
	 */
	public static Matrix blockMatrix(Matrix A, Matrix B) {
		if (A.m != B.m)
			throw new IllegalArgumentException("Matrices have unequal number of rows.");
		Matrix block = new Matrix(A.m, A.n + B.n);
		for (int i = 0; i < A.m; i++) {
			System.arraycopy(A.data, i * A.n, block.data, 2 * i * A.n, A.n);
			System.arraycopy(B.data, i * A.n, block.data, (2 * i + 1) * A.n, A.n);
		}
		return block;
	}
}
