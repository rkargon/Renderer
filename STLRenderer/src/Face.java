import java.util.Comparator;
import java.util.List;

public class Face {
	public Vertex normal;
	public MeshVertex[] vertices;
	public Object3D obj;
	
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
	public Face(Vertex normal, MeshVertex v1, MeshVertex v2, MeshVertex v3, Object3D obj) {
		vertices = new MeshVertex[] { v1, v2, v3 };
		v1.faces.add(this);
		v2.faces.add(this);
		v3.faces.add(this);

		if (normal != null && isPerpendicular(normal)) {
			normal.normalize();//normalizin' the normal!
			this.normal = normal;
		}
		else this.normal = generateNormal();

		this.obj = obj;
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
	public Face(Vertex normal, MeshVertex[] v, Object3D obj) {
		this(normal, v[0], v[1], v[2], obj);
	}

	public boolean isPerpendicular(Vertex normal) {
		Vertex side1 = vertices[1].subtract(vertices[0]);
		Vertex side2 = vertices[2].subtract(vertices[1]);

		return (normal.dotproduct(side1) == 0 && normal.dotproduct(side2) == 0 && normal
				.lensquared() > 0);
	}

	public Vertex generateNormal() {
		Vertex side1 = vertices[1].subtract(vertices[0]);
		Vertex side2 = vertices[2].subtract(vertices[1]);

		Vertex n = side1.crossproduct(side2);
		n.normalize();
		return n;
	}

	public Vertex center() {
		double x = vertices[0].x + vertices[1].x + vertices[2].x;
		double y = vertices[0].y + vertices[1].y + vertices[2].y;
		double z = vertices[0].z + vertices[1].z + vertices[2].z;

		return new Vertex(x / 3, y / 3, z / 3);
	}

	/**
	 * Determines if a face is in a certain range, along a given axis.
	 * 
	 * @param min
	 * @param max
	 * @param axis
	 *            0 = x axis, 1 = y axis, 2+ = z axis
	 * @return Whether or not the face overlaps the range [min, max] along the
	 *         given axis
	 */
	boolean inRange(double min, double max, int axis) {
		if (max < min) {
			double tmp = min;
			min = max;
			max = tmp;
		}

		double fmin = minCoord(axis);
		double fmax = maxCoord(axis);

		return (fmin <= max && fmax >= min);
	}

	boolean inBounds(Vertex[] b) {
		return inRange(b[0].x, b[1].x, 0) && inRange(b[0].y, b[1].y, 1)
				&& inRange(b[0].z, b[1].z, 2);
	}

	double minCoord(int axis) {
		return Math
				.min(Math.min(vertices[0].get(axis), vertices[1].get(axis)), vertices[2]
						.get(axis));

	}

	double maxCoord(int axis) {
		return Math
				.max(Math.max(vertices[0].get(axis), vertices[1].get(axis)), vertices[2]
						.get(axis));

	}

	// Fast, Minimum Storage Ray-Triangle Intersection
	//
	// Tomas Moller
	// Prosolvia Clarus AB
	// Sweden
	// tompa@clarus.se
	//
	// Ben Trumbore
	// Cornell University
	// Ithaca, New York
	// wbt@graphics.cornell.edu
	boolean intersectRayTriangle(Vertex orig, Vertex dir, Vertex tuv) {
		Vertex edge1, edge2, tvec, pvec, qvec;
		double det, inv_det;
		double epsilon = 0.00000001;
		double t, u, v;

		// find vectors for two edges sharing v0
		edge1 = vertices[1].subtract(vertices[0]);
		edge2 = vertices[2].subtract(vertices[0]);

		// begin calculating determinant - also used to calculate U parameter
		pvec = dir.crossproduct(edge2);

		// if determinant is near zero, ray lies in plane of triangle
		det = edge1.dotproduct(pvec);

		// calculate distance from v0 to ray origin
		tvec = orig.subtract(vertices[0]);
		inv_det = 1.0 / det;

		qvec = tvec.crossproduct(edge1);

		if (det > epsilon) {
			u = tvec.dotproduct(pvec);
			if (u < 0.0 || u > det) return false;

			// calculate V parameter and test bounds
			v = dir.dotproduct(qvec);
			if (v < 0.0 || u + v > det) return false;
		}
		else if (det < -epsilon) {
			// calculate U parameter and test bounds
			u = tvec.dotproduct(pvec);
			if (u > 0.0 || u < det) return false;

			// calculate V parameter and test bounds
			v = dir.dotproduct(qvec);
			if (v > 0.0 || u + v < det) return false;
		}
		else {
			return false;  // ray is parallel to the plane of the triangle
		}

		u *= inv_det;
		v *= inv_det;
		t = edge2.dotproduct(qvec) * inv_det;

		if (tuv != null) {
			tuv.x = t;
			tuv.y = u;
			tuv.z = v;
		}

		return (t > epsilon);
	}

	public Vertex[] bounds(){
		double xmin = minCoord(0);
		double xmax = maxCoord(0);
		double ymin = minCoord(1);
		double ymax = maxCoord(1);
		double zmin = minCoord(2);
		double zmax = maxCoord(2);
		return new Vertex[]{new Vertex(xmin, ymin, zmin), new Vertex(xmax, ymax, zmax)};
	}
	
	/**
	 * Checks for intersection between a ray and an axis-aligned bounding box
	 * (AABB)
	 * 
	 * @param bounds
	 *            The bounds of the bounding box
	 * @param origin
	 *            The origin of the ray
	 * @param dir
	 *            The direction of the ray
	 * @return Whether the ray intersects the bounding box
	 */
	public static boolean rayAABBIntersect(Vertex[] bounds, Vertex origin, Vertex dir) {
		//System.out.println("AABB intersect: "+bounds[0]+", "+bounds[1]);
		double tmp;
		double tmin = (bounds[0].x - origin.x) / dir.x;
		double tmax = (bounds[1].x - origin.x) / dir.x;
		if (tmin > tmax) {
			tmp = tmin;
			tmin = tmax;
			tmax = tmp;
		}

		double tymin = (bounds[0].y - origin.y) / dir.y;
		double tymax = (bounds[1].y - origin.y) / dir.y;
		if (tymin > tymax) {
			tmp = tymin;
			tymin = tymax;
			tymax = tmp;
		}

		if (tmin > tymax || tmax < tymin) return false;

		if (tymin > tmin) tmin = tymin;
		if (tymax < tmax) tmax = tymax;

		double tzmin = (bounds[0].z - origin.z) / dir.z;
		double tzmax = (bounds[1].z - origin.z) / dir.z;
		if (tzmin > tzmax) {
			tmp = tzmin;
			tzmin = tzmax;
			tzmax = tmp;
		}

		if (tmin > tzmax || tmax < tzmin) return false;

		if (tzmin > tmin) tmin = tzmin;
		if (tzmax < tmax) tmax = tzmax;
		
		if(tmin<=0 && tmax<=0) return false;//ray should not intersect bounding box behind it

		return true;
	}

	public String toString() {
		return "Face: [Normal: " + normal + ", V1: " + vertices[0] + ", V2: "
				+ vertices[1] + ", V3: " + vertices[2] + "]";
	}

	public static Vertex[] calcBoundingBox(List<Face> faces) {
		double minx, miny, minz, maxx, maxy, maxz, tmp;
		minx = miny = minz = maxx = maxy = maxz = Double.NaN;

		for (Face f : faces) {
			tmp = f.minCoord(0);
			if (minx != minx || minx > tmp) minx = tmp;
			tmp = f.minCoord(1);
			if (miny != miny || miny > tmp) miny = tmp;
			tmp = f.minCoord(2);
			if (minz != minz || minz > tmp) minz = tmp;

			tmp = f.maxCoord(0);
			if (maxx != maxx || maxx < tmp) maxx = tmp;
			tmp = f.maxCoord(1);
			if (maxy != maxy || maxy < tmp) maxy = tmp;
			tmp = f.maxCoord(2);
			if (maxz != maxz || maxz < tmp) maxz = tmp;
		}

		Vertex min = new Vertex(minx, miny, minz);
		Vertex max = new Vertex(maxx, maxy, maxz);

		return new Vertex[] { min, max };
	}

	public static double surfaceArea(Vertex[] b) {
		return surfaceArea(b[0], b[1]);
	}

	public static double surfaceArea(Vertex low, Vertex high) {
		double dx = Math.abs(high.x - low.x), dy = Math.abs(high.y - low.y), dz = Math
				.abs(high.z - low.z);
		double area = 2 * (dx * dy + dx * dz + dy * dz);
		return area;
	}

	public static void main(String[] args) {
		System.out.println(surfaceArea(new Vertex[] { new Vertex(1, 2, 3),
				new Vertex(11, 12, 3.1) }));
	}
}
