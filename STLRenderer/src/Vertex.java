import java.util.Comparator;
import java.util.List;

public class Vertex {
	public double x, y, z;
	
	public Vertex(double x, double y, double z) {
		super();
		this.x = x;
		this.y = y;
		this.z = z;
	}

	public void setVertex(Vertex v) {
		this.x = v.x;
		this.y = v.y;
		this.z = v.z;
	}

	public double get(int axis) {
		if (axis == 0) return x;
		else if (axis == 1) return y;
		else if (axis == 2) return z;
		else return get(Math.abs(axis) % 3);
	}

	public void set(int axis, double val) {
		if (axis == 0) x = val;
		else if (axis == 1) y = val;
		else if (axis == 2) z = val;
		else set(Math.abs(axis) % 3, val);
	}

	public Vertex subtract(Vertex v) {
		return new Vertex(x - v.x, y - v.y, z - v.z);
	}

	public Vertex add(Vertex v) {
		return new Vertex(x + v.x, y + v.y, z + v.z);
	}

	public Vertex scalarproduct(double a) {
		return new Vertex(x * a, y * a, z * a);
	}

	public double dotproduct(Vertex v) {
		return x * v.x + y * v.y + z * v.z;
	}

	public Vertex crossproduct(Vertex v) {
		return new Vertex(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y
				* v.x);
	}

	public double lensquared() {
		return x * x + y * y + z * z;
	}

	public double length() {
		return Math.sqrt(x * x + y * y + z * z);
	}

	/**
	 * Returns a unit vector in the direction of this vector.
	 * If the vector is zero, returns a zero vector
	 * 
	 * @return a unit vector in the direction of this vector, or a zero vector
	 */
	public Vertex getUnitVector() {
		double len = x * x + y * y + z * z;
		if (len == 0) {
			return clone();
		}
		if (len != 1 && len > 0) {
			len = Math.sqrt(len);
			return new Vertex(x / len, y / len, z / len);
		}
		else return clone();
	}

	public void normalize() {
		double len = x * x + y * y + z * z;
		if (len != 1 && len > 0) {
			len = Math.sqrt(len);
			x /= len;
			y /= len;
			z /= len;
		}
	}

	public Vertex reflection(Vertex normal) {
		//V- N * (2 * (V.N))
		return this.subtract(normal.scalarproduct(2 * this.dotproduct(normal)));
	}

	/**
	 * Generates a hash code based on the values of the x,y,z coordinates
	 * Two Vertex objects have the same hash if they have the same coordinates
	 */
	public int hashCode() {
		int result = 1;

		long bits = Double.doubleToLongBits(x);
		result = 31 * result + (int) (bits ^ (bits >>> 32));
		bits = Double.doubleToLongBits(y);
		result = 31 * result + (int) (bits ^ (bits >>> 32));
		bits = Double.doubleToLongBits(z);
		result = 31 * result + (int) (bits ^ (bits >>> 32));

		return result;
	}

	public boolean equals(Object o) {
		if (o instanceof Vertex) {
			Vertex v = (Vertex) o;
			return x == v.x && y == v.y && z == v.z;
		}
		return false;
	}

	public Vertex clone() {
		return new Vertex(x, y, z);
	}

	public String toString() {
		return "(" + x + ", " + y + ", " + z + ")";
	}

	public static Comparator<Vertex> vertexsorter = new Comparator<Vertex>(){
		@Override
		public int compare(Vertex v1, Vertex v2) {
			int c = Double.compare(v1.x, v2.x);
			if(c==0){
				c = Double.compare(v1.y, v2.y);
				if(c == 0) c = Double.compare(v1.z, v2.z);
			}
			return c;
		}
	};
	
	// Test method for math functions
	public static void main(String[] args) {
		System.out.println(new Vertex(2, 1, -1)
				.crossproduct(new Vertex(-3, 4, 1)));
	}

	public static Vertex ORIGIN() {
		return new Vertex(0, 0, 0);
	}

	public static Matrix MatrixRows(Vertex... verts) {
		Matrix m = new Matrix(verts.length, 3);

		for (int i = 0; i < verts.length; i++) {
			m.set(i, 0, verts[i].x);
			m.set(i, 1, verts[i].y);
			m.set(i, 2, verts[i].z);
		}

		return m;
	}

	/**
	 * Interpolates three vector values using the given barycentric coordinates.
	 * This can be used to interpolate the coordinates of a point on a face, or
	 * to interpolate normals
	 * 
	 * @return
	 */
	public static Vertex lerp(Vertex v1, Vertex v2, Vertex v3, double w1,
			double w2, double w3) {
		double x = w1 * v1.x + w2 * v2.x + w3 * v3.x;
		double y = w1 * v1.y + w2 * v2.y + w3 * v3.y;
		double z = w1 * v1.z + w2 * v2.z + w3 * v3.z;
		return new Vertex(x, y, z);
	}

	public static Vertex lerp(Vertex v1, Vertex v2, double r) {
		//v1 + (v2-v1)*r
		return new Vertex(v1.x + (v2.x - v1.x) * r, v1.y + (v2.y - v1.y) * r, v1.z
				+ (v2.z - v1.z) * r);
	}

	public static Vertex[] listBounds(List<Vertex> vertices) {
		double minx, miny, minz, maxx, maxy, maxz;
		minx = miny = minz = Double.POSITIVE_INFINITY;
		maxx = maxy = maxz = Double.NEGATIVE_INFINITY;
		for (Vertex v : vertices) {
			if (v.x < minx) minx = v.x;
			if (v.x > maxx) maxx = v.x;
			if (v.y < miny) miny = v.y;
			if (v.y > maxy) maxy = v.y;
			if (v.z < minz) minz = v.z;
			if (v.z > maxz) maxz = v.z;
		}
		return new Vertex[] { new Vertex(minx, miny, minz),
				new Vertex(maxx, maxy, maxz) };
	}

	public static Vertex[] listBounds(Vertex... vertices) {
		double minx, miny, minz, maxx, maxy, maxz;
		minx = miny = minz = Double.POSITIVE_INFINITY;
		maxx = maxy = maxz = Double.NEGATIVE_INFINITY;
		for (Vertex v : vertices) {
			if (v.x < minx) minx = v.x;
			if (v.x > maxx) maxx = v.x;
			if (v.y < miny) miny = v.y;
			if (v.y > maxy) maxy = v.y;
			if (v.z < minz) minz = v.z;
			if (v.z > maxz) maxz = v.z;
		}
		return new Vertex[] { new Vertex(minx, miny, minz),
				new Vertex(maxx, maxy, maxz) };
	}

	public static Vertex[] intersectBoundingBoxes(Vertex[] b1, Vertex[] b2){
		Vertex[] newbounds = new Vertex[]{new Vertex(0,0,0), new Vertex(0,0,0)};
		newbounds[0] = max3(b1[0], b2[0]);
		newbounds[1] = min3(b1[1], b2[1]);
		return newbounds;
	}

	public static Vertex max3(Vertex v1, Vertex v2) {
		return new Vertex(Math.max(v1.x, v2.x), Math.max(v1.y, v2.y), Math.max(v1.z, v2.z));
	}

	public static Vertex min3(Vertex v1, Vertex v2) {
		return new Vertex(Math.min(v1.x, v2.x), Math.min(v1.y, v2.y), Math.min(v1.z, v2.z));
	}
}
