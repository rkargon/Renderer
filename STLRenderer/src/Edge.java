
public class Edge {
	MeshVertex v1, v2;

	public Edge(MeshVertex v1, MeshVertex v2) {
		this.v1 = v1;
		this.v2 = v2;
	}

	public Vertex getVector(){
		return v2.subtract(v1);
	}
	
	public double length(){
		return v2.subtract(v1).length();
	}
	
	public int hashCode() {
		int hash1 = v1.hashCode();
		int hash2 = v2.hashCode();

		int result;
		//hash of Edge(v1, v2) should equal hash of Edge(v2, v1)
		if (hash1 > hash2) result = 31 * (31 + hash1) + hash2;
		else result = 31 * (31 + hash2) + hash1;
		return result;
	}
	

	public boolean equals(Object o) {
		if (o instanceof Edge) {
			Edge e = (Edge) o;
			if(v1.equals(e.v1)){
				return v2.equals(e.v2);
			}
			else if (v1.equals(e.v2)){
				return v2.equals(e.v1);
			}
		}
		return false;
	}
}
