import java.util.ArrayList;


public class MeshVertex extends Vertex {
	
	public ArrayList<Face> faces;
	
	public MeshVertex(double x, double y, double z) {
		super(x, y, z);
		faces = new ArrayList<Face>();
	}
	
	public Vertex vertexNormal(){
		Vertex n = new Vertex(0, 0, 0);
		
		for(Face f : faces){
			n=n.add(f.normal);
		}
		n.normalize();
		
		return n;
	}

}
