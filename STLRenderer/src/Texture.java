import java.awt.Color;


public interface Texture {
	
	//get color from texture given a face barycentric coordinates of vertex
	public Color getCol(Face f, double w1, double w2);
	
	//get value from texture given a face barycentric coordinates of vertex
	public double getVal(Face f, double w1, double w2);
}
