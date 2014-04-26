public class Face {
	public Vertex normal;
	public Vertex[] vertices;

	/**
	 * Creates a face based off of three vertices and a normal. If the normal
	 * does not correspond to the triangle, a new normal will be
	 * created, based on the order of the vertices.
	 * 
	 * The new normal will be based off of (v2-v1) X (v3-v2)
	 * 
	 * @param normal
	 * @param v1
	 * @param v2
	 * @param v3
	 */
	public Face(Vertex normal, Vertex v1, Vertex v2, Vertex v3) {
		vertices = new Vertex[] { v1, v2, v3 };
		if (isPerpendicular(normal)) {
			normal.normalize();//normalizin' the normal!
			this.normal = normal;
		}
		else this.normal = generateNormal();
	}

	/**
	 * Creates a face based on a normal and an array of 3 vertices.
	 * 
	 * @see {@link #Face(Vertex, Vertex, Vertex, Vertex)}
	 * 
	 * @param normal
	 *            the normal to the face
	 * @param v
	 *            the array of vertices
	 */
	public Face(Vertex normal, Vertex[] v) {
		this(normal, v[0], v[1], v[2]);
	}

	public boolean isPerpendicular(Vertex normal) {
		Vertex side1 = vertices[1].subtract(vertices[0]);
		Vertex side2 = vertices[2].subtract(vertices[1]);

		return (normal.dotproduct(side1) == 0 && normal.dotproduct(side2) == 0 && normal.lensquared()>0);
	}

	public Vertex generateNormal() {
		Vertex side1 = vertices[1].subtract(vertices[0]);
		Vertex side2 = vertices[2].subtract(vertices[1]);

		Vertex n = side1.crossproduct(side2);
		n.normalize();
		return n;
	}
	
	public Vertex center(){
		double x = vertices[0].x+vertices[1].x+vertices[2].x;
		double y = vertices[0].y+vertices[1].y+vertices[2].y;
		double z = vertices[0].z+vertices[1].z+vertices[2].z;
		
		return new Vertex(x/3, y/3, z/3);
	}
	
	public String toString(){
		return "Face: [Normal: "+normal+", V1: "+vertices[0]+", V2: "+vertices[1]+", V3: "+vertices[2]+"]";
	}
}
