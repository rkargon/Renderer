import java.awt.Point;
import java.util.Comparator;

public class Camera {
	public Vertex center;
	public Vertex focus;
	public Vertex normal;
	public Vertex vert;
	public double fov;
	public double mindist, maxdist;
	public boolean ortho;

	//temporary vars used for ray casting.
	private Vertex cx, cy; //horizontal and vertical vectors of image plane

	public Camera() {
		this(new Vertex(-5, 0, 0), Vertex.ORIGIN(), new Vertex(1, 0, 0), new Vertex(0, 0, 1), 0.75, 0.01, 100, false);
	}

	/**
	 * 
	 * @param center
	 *            Location of camera
	 * @param normal
	 *            Direction of view
	 * @param vert
	 *            Orientation (up direction) of view
	 * @param fov
	 *            Field of view (radians)
	 * @param asp
	 *            the aspect ratio
	 * @param ortho
	 *            whether view is perspective or orthogonal
	 */
	public Camera(Vertex center, Vertex focus, Vertex normal, Vertex vert, double fov, double mindist, double maxdist, boolean ortho) {
		super();

		if (fov == 0)
			throw new IllegalArgumentException("Field of view cannot be 0.");
		if (mindist <= 0)
			throw new IllegalArgumentException("Minimum clipping distance must be greater than 0.");
		if (normal.dotproduct(vert) != 0)
			throw new IllegalArgumentException("Normal and vertical vectors of camera must be perpendicular");

		this.center = center;
		this.focus = focus;
		this.normal = normal;
		this.vert = vert;
		this.fov = fov;
		this.mindist = mindist;
		this.maxdist = maxdist > 0 ? maxdist : Double.POSITIVE_INFINITY;
		this.ortho = ortho;

		calcImageVectors();
	}

	private void calcImageVectors() {
		cx = normal.crossproduct(vert);
		cy = cx.crossproduct(normal);
		cx.normalize();
		cy.normalize();
	}

	public void centerFocus() {
		center = focus.add(normal.scalarproduct(-focus.subtract(center)
				.length()));
	}

	public void zoom(double zoomfactor) {
		center = focus.add(center.subtract(focus).scalarproduct(zoomfactor));
	}

	public void shiftFocus(double dx, double dy) {
		focus = focus.add(getImagePlaneVector(dx, dy));
	}

	public Vertex getImagePlaneVector(double dx, double dy) {
		return cx.scalarproduct(dx).add(cy.scalarproduct(dy));
	}

	public void rotateAxis(Vertex a, double dtheta) {
		a.normalize();
		double l = a.x, m = a.y, n = a.z;
		double s = Math.sin(dtheta), c = Math.cos(dtheta);
		Matrix axes = Vertex.MatrixRows(vert, normal);
		Matrix rot = new Matrix(new double[][] {
				{ l * l * (1 - c) + c, m * l * (1 - c) - n * s,
						n * l * (1 - c) + m * s },
				{ l * m * (1 - c) + n * s, m * m * (1 - c) + c,
						n * m * (1 - c) - l * s },
				{ l * n * (1 - c) - m * s, m * n * (1 - c) + l * s,
						n * n * (1 - c) + c } });
		axes = axes.product(rot);

		vert = new Vertex(axes.fastget(0), axes.fastget(1), axes.fastget(2));
		normal = new Vertex(axes.fastget(3), axes.fastget(4), axes.fastget(5));
		calcImageVectors();
	}

	public void rotateLocalZ(double dtheta) {
		rotateAxis(normal, dtheta);
	}

	public void rotateLocalY(double drho) {
		rotateAxis(vert, drho);
	}

	public void rotateLocalX(double dpsi) {
		rotateAxis(normal.crossproduct(vert), dpsi);
	}

	public void setGlobalRotation(double theta, double rho, double psi) {
		double st = Math.sin(theta), ct = Math.cos(theta), sr = Math.sin(rho), cr = Math
				.cos(rho), sp = Math.sin(psi), cp = Math.cos(psi);
		vert = new Vertex(-sr, cr * sp, cr * cp);
		normal = new Vertex(ct * cr, ct * sr * sp - cp * st, st * sp + ct * cp
				* sr);
		calcImageVectors();
	}

	public Point projectVertex(Vertex v, int w, int h) {
		Point p;
		Vertex dv = v.subtract(center);
		double x, y;

		if (ortho) {
			x = dv.dotproduct(cx);
			y = dv.dotproduct(cy);
		}
		else {

			double len = dv.length();
			if (dv.dotproduct(normal) <= 0 || len < mindist || len > maxdist) {
				return null;
			}
			Vertex proj_horiz = dv
					.subtract(cy.scalarproduct(dv.dotproduct(cy))); //projection onto horizontal plane
			Vertex proj_vert = dv.subtract(cx.scalarproduct(dv.dotproduct(cx))); //projection onto horizontal plane

			x = Math.asin(proj_horiz.dotproduct(cx) / proj_horiz.length());
			y = Math.asin(proj_vert.dotproduct(cy) / proj_vert.length());
		}

		double px = (0.5 + x/fov) * w;
		double py = (0.5 - y*w/(h*fov)) * h;
		p = new Point((int) px, (int) py);
		return p;
	}

	public Comparator<Face> zsortFaces = new Comparator<Face>() {
		public int compare(Face f1, Face f2) {
			return Double.compare(faceDepth(f1), faceDepth(f2));
		};
	};
	public Comparator<Vertex> zsortVertices = new Comparator<Vertex>() {
		public int compare(Vertex v1, Vertex v2) {
			return Double.compare(vertexDepth(v1), vertexDepth(v2));
		};
	};

	public double vertexDepth(Vertex v) {
		v = v.subtract(center);

		if (ortho) {
			return v.dotproduct(normal);
		}
		else {
			double d = v.length();
			if (v.dotproduct(normal) < 0) d = -d;
			return d;
		}
	}

	public double faceDepth(Face f) {
		return (vertexDepth(f.vertices[0]) + vertexDepth(f.vertices[1]) + vertexDepth(f.vertices[2])) / 3;
	}

	public Vertex[] castRay(double px, double py, int w, int h) {
		double img_w = 2 * Math.tan(fov / 2);
		double img_h = (img_w * h / w);

		double x = px / w - 0.5;
		double y = py / h - 0.5;

		Vertex[] ray = new Vertex[2];

		if (ortho) {
			ray[1] = normal;
			//center + cx * (x*img_w) + cy * (y*img_h)
			ray[0] = center.add(cx.scalarproduct(x * img_w))
					.add(cy.scalarproduct(-y * img_h)); //-y because of java panel coordinates
			return ray;
		}
		else {
			ray[0] = center;
			//normal + cx * (x*img_w) + cy * (y*img_h)
			ray[1] = normal.add(cx.scalarproduct(x * img_w))
					.add(cy.scalarproduct(-y * img_h)).getUnitVector();
			return ray;
		}
	}

	public Vertex viewVector(Vertex v) {
		if (ortho) {
			return normal.scalarproduct(v.dotproduct(normal));
		}
		else {
			return v.subtract(center);
		}
	}

	public static void main(String[] args) {
		Camera cam = new Camera();
		Vertex[] ray = cam.castRay(0, 0, 100, 100);
		Vertex v = ray[0].add(ray[1]);
		System.out.println(cam.projectVertex(v, 100, 100));
	}

}
