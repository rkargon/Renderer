public class Vertex {
	public double x, y, z;

	public Vertex(double x, double y, double z) {
		super();
		this.x = x;
		this.y = y;
		this.z = z;
	}

	public Vertex subtract(Vertex v) {
		return new Vertex(x - v.x, y - v.y, z - v.z);
	}

	public Vertex add(Vertex v) {
		return new Vertex(x + v.x, y + v.y, z + v.y);
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
		return Math.sqrt(lensquared());
	}

	/**
	 * Returns a unit vector in the direction of this vector. 
	 * If the vector is zero, returns a zero vector
	 * @return a unit vector in the direction of this vector, or a zero vector
	 */
	public Vertex getUnitVector(){
		double len = x * x + y * y + z * z;
		if(len==0){
			return this;
		}
		if (len != 1 && len > 0) {
			len = Math.sqrt(len);
			return new Vertex(x/len, y/len, z/len);
		}
		else return this;
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

	public Vertex clone(){
		return new Vertex(x, y, z);
	}
	
	public String toString() {
		return "<" + x + ", " + y + ", " + z + ">";
	}

	// Test method for math functions
	public static void main(String[] args) {
		System.out.println(new Vertex(2, 1, -1)
				.crossproduct(new Vertex(-3, 4, 1)));
	}
	
	public static Vertex ORIGIN(){
		return new Vertex(0, 0, 0);
	}

}
