
public class SingularMatrixException extends Exception {
	
	public SingularMatrixException(Matrix M){
		this(M, "");
	}
	
	public SingularMatrixException(Matrix M, String s){
		super(M+" is a singular matrix. "+s);
	}
}
