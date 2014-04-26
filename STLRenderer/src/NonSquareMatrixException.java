
public class NonSquareMatrixException extends Exception {
	
	public NonSquareMatrixException(int m, int n){
		this(m, n, "");
	}
	
	public NonSquareMatrixException(int m, int n, String s){
		super("Matrix ("+m+" rows, "+n+" columns) is not square. "+s);
	}
}
