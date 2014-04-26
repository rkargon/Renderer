import java.util.Comparator;

public class Camera {
	//coordinates for the point the camera focuses on
	public Vertex center;

	//Normal vector, vector pointing 'up' relative to camera, and horizontal vector, pointing 'left' relative to camera
	private Matrix axes;
	private double theta, rho, psi;

	//horizontal field of view (radians)
	public double fov;

	public double mindist;//minimum distance from the center of the camera for something to be rendered

	/* COMPARATORS */
	//sort vertices and faces by their distance to the camera (in REVERSE order, so iterating goes from far to near

	public Comparator<Face> zsortFaces = new Comparator<Face>() {

		@Override
		public int compare(Face f1, Face f2) {
			Vertex f1v1 = f1.vertices[0].subtract(center), f1v2 = f1.vertices[1]
					.subtract(center), f1v3 = f1.vertices[2].subtract(center);
			Vertex f2v1 = f2.vertices[0].subtract(center), f2v2 = f2.vertices[1]
					.subtract(center), f2v3 = f2.vertices[2].subtract(center);

			double f1d = (f1v1.length() + f1v2.length() + f1v3.length()) / 3;
			double f2d = (f2v1.length() + f2v2.length() + f2v3.length()) / 3;
			return Double.compare(f2d, f1d);
		}

	};

	public Comparator<Vertex> zsortVerts = new Comparator<Vertex>() {

		@Override
		public int compare(Vertex v1, Vertex v2) {
			v1 = v1.subtract(center);
			v2 = v2.subtract(center);
			double d1 = v1.dotproduct(normal());
			double d2 = v2.dotproduct(normal());
			return Double.compare(d2, d1);
		}

	};

	public Camera() {
		this(0, 0, 0, 0, 0, 0, 1, 0.01);
	}

	/**
	 * Creates a camera object with the given center, rotation, field of view,
	 * and minimum clipping distance.
	 * 
	 * @param x
	 *            The x-coordinate of the center of the camera
	 * @param y
	 *            The y-coordinate of the center of the camera
	 * @param z
	 *            The z-coordinate of the center of the camera
	 * @param theta
	 *            The (global) rotation along the vertical axis of the camera
	 * @param rho
	 *            The (global) rotation along the horizontal axis of the camera
	 * @param psi
	 *            The (global) rotation along the normal axis of the camera
	 * @param fov
	 *            The horizontal field of view, in radians, of the camera. View
	 *            extends
	 *            from -fov/2 to fov/2 radians from the camera's normal along
	 *            the horizontal. (Field of view may be smaller along vertical,
	 *            depending on aspect ratio of image). If
	 *            fov>180, images will not repeat
	 * @param mindist
	 */
	public Camera(double x, double y, double z, double theta, double rho, double psi, double fov, double mindist) {
		this.center = new Vertex(x, y, z);
		this.theta = theta;
		this.rho = rho;
		this.psi = psi;

		calcNormals();

		if (fov <= 0)
			throw new IllegalArgumentException("Field of view (" + fov
					+ ") must be positive.");
		this.fov = fov;
		this.mindist = Math.max(mindist, 0);
	}

	/**
	 * Sets up axes matrix based on global rotations
	 * All I know is that theta, rho ,and psi apply to GLOBAL Z, Y, and X axes
	 */
	public void calcNormals() {
		double st = Math.sin(theta), ct = Math.cos(theta), sr = Math.sin(rho), cr = Math
				.cos(rho), sp = Math.sin(psi), cp = Math.cos(psi);

		axes = new Matrix(new double[][] {
				//horizontal, ie local X
				{ cr * st, ct * cp + st * sr * sp, cp * st * sr - ct * sp },
				//vertical, ie local Y
				{ -sr, cr * sp, cr * cp },
				//normal, ie local Z
				{ ct * cr, ct * sr * sp - cp * st, st * sp + ct * cp * sr } });
	}

	public Vertex horizontal() {
		return new Vertex(axes.get(0, 0), axes.get(0, 1), axes.get(0, 2));
	}

	public Vertex vertical() {
		return new Vertex(axes.get(1, 0), axes.get(1, 1), axes.get(1, 2));
	}

	public Vertex normal() {
		return new Vertex(axes.get(2, 0), axes.get(2, 1), axes.get(2, 2));
	}

	public double theta() {
		return theta;
	}

	public double rho() {
		return rho;
	}

	public double psi() {
		return psi;
	}

	public void setGlobalRotations(double theta, double rho, double psi) {
		this.theta = theta;
		this.rho = rho;
		this.psi = psi;

		calcNormals();

	}

	public void rotateGlobalZ(double dtheta) {
		theta += dtheta;
		double st = Math.sin(dtheta), ct = Math.cos(dtheta);
		Matrix rotZ = new Matrix(new double[] { ct, -st, 0, st, ct, 0, 0, 0, 1 }, 3, 3);
		axes = axes.product(rotZ);
	}

	public void rotateGlobalY(double drho) {
		rho += drho;
		double sr = Math.sin(drho), cr = Math.cos(drho);
		Matrix rotY = new Matrix(new double[] { cr, 0, sr, 0, 1, 0, -sr, 0, cr }, 3, 3);
		axes = axes.product(rotY);
	}

	public void rotateGlobalX(double dpsi) {
		psi += dpsi;
		double sp = Math.sin(dpsi), cp = Math.cos(dpsi);
		Matrix rotX = new Matrix(new double[] { 1, 0, 0, 0, cp, -sp, 0, sp, cp }, 3, 3);
		axes = axes.product(rotX);
	}

	public void rotateLocalZ(double dtheta){
		rotateAxis(normal(), dtheta);
	}
	
	public void rotateLocalY(double drho){
		rotateAxis(vertical(), drho);
	}
	
	public void rotateLocalX(double dpsi){
		rotateAxis(horizontal(), dpsi);
	}

	public void rotateAxis(Vertex a, double dtheta) {
		a.normalize();
		double l = a.x, m = a.y, n = a.z;
		double s = Math.sin(dtheta), c = Math.cos(dtheta);
		Matrix rot = new Matrix(new double[][] {
				{ l * l * (1 - c) + c, m * l * (1 - c) - n * s,
						n * l * (1 - c) + m * s },
				{ l * m * (1 - c) + n * s, m * m * (1 - c) + c,
						n * m * (1 - c) - l * s },
				{ l * n * (1 - c) - m * s, m * n * (1 - c) + l * s,
						n * n * (1 - c) + c } });
		axes=axes.product(rot);
	}

	public void centerOrigin() {
		center = normal().scalarproduct(-center.length());
	}

	//TODO centerVertex(Vertex v)

	public boolean verifyVectors() {
		Vertex normal = normal(), horizontal = horizontal(), vertical = vertical();

		System.out.println("normal: " + normal);
		System.out.println("horiz: " + horizontal);
		System.out.println("vert: " + vertical);
		System.out.println(normal().dotproduct(vertical));
		System.out.println(normal().dotproduct(horizontal));
		System.out.println(horizontal().dotproduct(vertical));
		return (normal.dotproduct(vertical) == 0
				&& normal.dotproduct(horizontal) == 0 && horizontal
					.dotproduct(vertical) == 0);
	}

	public double vertexDistance(Vertex v) {
		v = v.subtract(center);
		return v.dotproduct(normal());
	}

	//returns the distance of a face from the camera, based on the average distance of the vertices
	public double faceDistance(Face f) {

		Vertex v1 = f.vertices[0].subtract(center), v2 = f.vertices[1]
				.subtract(center), v3 = f.vertices[2].subtract(center);
		Vertex normal = normal();

		return (v1.dotproduct(normal) + v2.dotproduct(normal) + v3
				.dotproduct(normal)) / 3;
	}

}
