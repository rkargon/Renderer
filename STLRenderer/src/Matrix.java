import java.awt.Dimension;
import java.util.Arrays;

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
	 * 
	 * @return a COPY of the data of the matrix
	 */
	public double[] data() {
		return Arrays.copyOf(data, data.length);
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
		if(cols.length==0) return this.clone();
		
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
		if(rows.length==0) return this.clone();
		
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

	/**
	 * Returns the adjoint (conjugate transpose) of a matrix
	 * Because complex numbers are not used, it is equivalent ot the transpose
	 * 
	 * @return The conjugate transpose of a matrix
	 */
	public Matrix adjoint() {
		//not dealing with complex numbers here
		return transpose();
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
	 * 
	 * @return The matrix of cofacors of this matrix
	 * @throws NonSquareMatrixException
	 */
	public Matrix cofactorMatrix() throws NonSquareMatrixException {
		Matrix cof = new Matrix(n, m); //matrix of cofactors
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				cof.data[i * n + j] = this.cofactor(i, j);
			}
		}

		return cof;
	}

	/**
	 * Returns the adjugate, or transpose of the matrix of cofactors, of the
	 * matrix
	 * 
	 * @return The transpose of the matrix of cofactors of the matrix
	 * @throws NonSquareMatrixException
	 */
	public Matrix adjugate() throws NonSquareMatrixException {
		Matrix adj = new Matrix(n, m); //matrix of cofactors
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				adj.data[j * m + i] = this.cofactor(i, j);
			}
		}

		return adj;
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

	public Matrix minorMatrix(int row, int col) {
		if (row < 0 || col < 0 || row >= m || col >= n)
			throw new IllegalArgumentException("Matrix indices are out of bounds.");
		Matrix minor = new Matrix(m - 1, n - 1);
		for (int i = 0; i < row; i++) {
			System.arraycopy(data, i * n, minor.data, i * minor.n, col);
			System.arraycopy(data, i * n + col + 1, minor.data, i * minor.n
					+ col, minor.n - col);
		}
		for (int i = row + 1; i < n; i++) {
			System.arraycopy(data, i * n, minor.data, (i - 1) * minor.n, col);
			System.arraycopy(data, i * n + col + 1, minor.data, (i - 1)
					* minor.n + col, minor.n - col);
		}
		return minor;
	}

	/**
	 * Returns the minor of the matrix
	 * 
	 * @param row
	 *            The row to exclude
	 * @param col
	 *            The column to exclude
	 * @return The minor of the matrix, with the given row and colum excluded
	 * @throws NonSquareMatrixException
	 */
	public double minor(int row, int col) throws NonSquareMatrixException {
		return minorMatrix(row, col).determinant();
	}

	/**
	 * Returns the cofactor, or signed minor, of the matrix
	 * 
	 * @param row
	 *            The row to exclude
	 * @param col
	 *            The column to exclude
	 * @return minor * (-1)^(row+col)
	 * @throws NonSquareMatrixException
	 */
	public double cofactor(int row, int col) throws NonSquareMatrixException {
		double det = minor(row, col);
		if ((row + col) % 2 == 1) det *= -1;
		return det;

	}

	/**
	 * Calculates trace, ie the sum of the diagonal values, of a matrix. If
	 * matrix is not
	 * square, adds up first min(m, n) diagonal items.
	 * 
	 * @return The trace of the matrix
	 */
	public double trace() {
		int range = Math.min(m, n);
		int sum = 0;
		for (int i = 0; i < range; i++) {
			sum += data[i * n + i];
		}
		return sum;
	}

	public double diagonalProduct() {
		int range = Math.min(m, n);
		int prod = 1;
		for (int i = 0; i < range; i++) {
			prod *= data[i * n + i];
		}

		return prod;
	}

	public double determinant() throws NonSquareMatrixException {
		if (n != m) throw new NonSquareMatrixException(m, n);
		Matrix tmp = this.clone();
		double det_coeff = tmp.GaussianElimination();
		return det_coeff * tmp.diagonalProduct();
	}

	public Matrix inverse() throws SingularMatrixException,
			NonSquareMatrixException {
		if (m != n) throw new NonSquareMatrixException(m, n);
		//create matrix |A I|
		Matrix block = Matrix.blockMatrix(this, Matrix.identityMatrix(m));
		double f = block.GaussianElimination();
		if (f * block.diagonalProduct() == 0)
			throw new SingularMatrixException(this);
		return block.segment(0, block.n / 2, block.m - 1, block.n - 1);
	}

	/**
	 * Calculates the Reduced Row Echelon Form (RREF) of a matrix through
	 * Gaussian elimination.
	 * NOTE: this algorithm is destructive, it modifies the replaces the current
	 * matrix with the RREF.
	 * 
	 * The determinant coefficient. This is multiplied by the diagonal product
	 * of
	 * the matrix in order to calculat the determinant
	 */
	public double GaussianElimination() {
		double det_coeff = 1; //coefficient for calculating determinant

		//k is stored as two variables because in some cases, pivot is not on diagonal
		for (int k_i = 0, k_j = 0; k_i < m && k_j < n; k_i++, k_j++) {
			int i_max = -1;
			double val_max = -1;
			while (val_max <= 0) {
				//if all remaining columns are 0s below row k_i
				if (k_j >= n) return det_coeff;

				for (int i = k_i; i < m; i++) {
					//find row with  largest (absolute) value in column
					if (val_max < Math.abs(data[i * n + k_j])) {
						val_max = Math.abs(data[i * n + k_j]);
						i_max = i;
					}
				}
				//skip this column if it has all 0s
				if (val_max == 0) {
					k_j++;
				}
			}

			//swap rows i_max and k
			swapRows(k_i, i_max);
			if (k_i != i_max) det_coeff *= -1;

			//Make this pivot 1
			double pv = data[k_i * n + k_j], pvInv = 1.0 / pv;
			data[k_i * n + k_j] = 1;
			for (int i = k_j + 1; i < n; i++) {
				data[k_i * n + i] *= pvInv;
			}
			det_coeff *= pv;

			//for all rows
			for (int i = 0; i < m; i++) {
				//if this row has a 0 at column k
				if (data[i * n + k_j] == 0) continue;

				if (i != k_i) {
					double factor = data[i * n + k_j];

					//for each column in this row
					for (int j = k_j; j < n; j++) {
						data[i * n + j] -= data[k_i * n + j] * factor;
					}
				}
			}
		}

		return det_coeff;
	}

	/**
	 * Rank-decomposes a matrix A (m x n) into matrices C (m x r) and F (r x n)
	 * such that A = C*F
	 * 
	 * @return An array containing C and F, in that order.
	 */
	public Matrix[] rankDecomposition() {
		Matrix B = this.clone();
		B.GaussianElimination();

		Matrix C = this.removeColumns(this.getNonPivotColumns());
		Matrix F = B.removeRows(B.getZeroRows());
		
		return new Matrix[] { C, F };
	}

	/**
	 * Returns the non-pivot columns of a matrix. Those are the columns that
	 * contain
	 * the leading nonzero numbers for each row. Usually this is done for
	 * row-reduced matrices, so all leading numbers are 1's. However, this
	 * function checks for any nonzero number as a leading entry.
	 * 
	 * NOTE: It is assumed that the leading coefficient of each row is to the
	 * right of the coefficient of the row above it. This condition is met for
	 * row echelon form matrices.
	 * 
	 * @return An array of the column indices for non-pivot columns
	 */
	public int[] getNonPivotColumns() {
		int row, col, ncols = 0;//ncols is the current index in the array of non=pivot columns found
		int[] nonpivots = new int[n];

		for (row = 0, col = 0; row < m && col < n; row++, col++) {
			if (data[row * n + col] == 0) {
				//add column to array
				nonpivots[ncols] = col;
				ncols++;

				//skip column for next iteration
				col++;
			}
		}

		return Arrays.copyOf(nonpivots, ncols);
	}

	public int[] getZeroRows() {
		int[] zerorows = new int[m];
		int nrows = 0; //index current of zerorows
		boolean isempty = true;

		for (int i = 0; i < m; i++) {
			isempty = true;
			for (int j = 0; j < n; j++) {
				if (data[i * n + j] != 0) {
					isempty = false;
					break;
				}
			}
			if (isempty) {
				zerorows[nrows] = i;
				nrows++;
			}
		}

		return Arrays.copyOf(zerorows, nrows);
	}

	public Matrix MoorePenrosePesudoinverse() {
		Matrix B, Btr, C, Ctr, RD[] = this.rankDecomposition();
		B = RD[0];
		C = RD[1];
		Btr = B.transpose();
		Ctr = C.transpose();
		
		try {
			Matrix B_plus = Btr.product(B).inverse().product(Btr);
			Matrix C_plus = Ctr.product(C.product(Ctr).inverse());
			return C_plus.product(B_plus);
		}
		catch (SingularMatrixException | NonSquareMatrixException e) {
			//Should not occur, every matrix has a pseudo-inverse
			e.printStackTrace();
			return null;
		}
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
		int hash = 0;

		//TODO use Arrays's hash function
		hash += m;
		hash += (hash << 10);
		hash ^= (hash >>> 6);

		hash += n;
		hash += (hash << 10);
		hash ^= (hash >>> 6);

		hash += Arrays.hashCode(data);
		hash += (hash << 10);
		hash ^= (hash >>> 6);

		hash += (hash << 3);
		hash += (hash >>> 11);
		hash += (hash << 15);
		return hash;
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

	/**
	 * Uses the Moore-Penrose pseudoinverse to calculate a linear fit for a set
	 * of data points.
	 * 
	 * @param X
	 *            Matrix of coefficients of linear equations, size m x n, m>n
	 * @param y
	 *            Matrix of constants, size m x 1
	 * @return
	 * @throws SingularMatrixException
	 * @throws NonSquareMatrixException
	 */
	public static Matrix linearLeastSquares(Matrix X, Matrix y)
			throws SingularMatrixException, NonSquareMatrixException {
		if (X.m != y.m || X.m <= X.n)
			throw new IllegalArgumentException("invalid matrix sizes!");

		return X.MoorePenrosePesudoinverse().product(y);
	}

	public static void main(String[] args) {

		try {
			Matrix X = new Matrix(new double[][] { { 1, 1 }, { 1, 2 },
					{ 1, 3 }, { 1, 4 } });
			Matrix Y = new Matrix(new double[][] { { 6 }, { 5 }, { 7 }, { 10 } });
			System.out.println(linearLeastSquares(X, Y));
		}
		catch (Exception e) {
			// Auto-generated catch block
			e.printStackTrace();
		}
	}
}
