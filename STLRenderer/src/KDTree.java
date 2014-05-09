import java.awt.Graphics2D;
import java.awt.geom.Line2D;
import java.awt.geom.Line2D.Double;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

public class KDTree {

	//current axis for this tree
	public int dir;
	//position of division along axis
	public double pos;

	//subtrees
	public KDTree lower;
	public KDTree upper;

	//leaf nodes store faces
	public List<Face> faces;

	public static final int MAXDEPTH = 10;
	public static final int MINTRIS = 64;

	public KDTree(List<Face> faces) {
		this(faces, calcBoundingBox(faces), 0);
	}

	private KDTree(List<Face> faces, Vertex[] bounds, int d) {
		if (faces.size() <= MINTRIS || d >= MAXDEPTH) {
			this.faces = faces;
			return;
		}

		//sort faces by current axis, get min and max values
		this.dir = d % 3;
		CompareFacesByAxis fcomp = new CompareFacesByAxis(dir);
		Collections.sort(faces, fcomp);

		//use median of face centers by default. If this is not in the bounding box, just use midpoint of bounding box
		this.pos = faces.get(faces.size() / 2).center().get(dir);
		double min = bounds[0].get(dir), max = bounds[1].get(dir);
		if (pos <= min || pos >= max) {
			pos = (min + max) / 2;
		}

		List<Face> lowerfaces = new ArrayList<Face>();
		List<Face> upperfaces = new ArrayList<Face>();
		for (Face f : faces) {
			if (f.inRange(min, pos, dir)) {
				lowerfaces.add(f);
			}
			if (f.inRange(pos, max, dir)) {
				upperfaces.add(f);
			}

		}

		Vertex[] bounds_cpy = new Vertex[] { bounds[0].clone(),
				bounds[1].clone() };
		bounds[1].set(dir, pos);
		bounds_cpy[0].set(dir, pos);
		lower = new KDTree(lowerfaces, bounds, d + 1);
		upper = new KDTree(upperfaces, bounds_cpy, d + 1);
	}

	public String toString() {
		return toString("", true);
	}

	public String toString(String prefix, boolean isTail) {
		String s = "";

		s += prefix + (isTail ? "\\-- " : "|-- ")
				+ (faces != null ? "leaf " + faces.size() : "") + "\n";

		if (faces == null) {
			s += lower.toString(prefix + (isTail ? "    " : "|   "), false);
			s += upper.toString(prefix + (isTail ? "    " : "|   "), true);
		}

		return s;
	}

	public static Vertex[] calcBoundingBox(List<Face> faces) {
		double minx, miny, minz, maxx, maxy, maxz;
		minx = miny = minz = maxx = maxy = maxz = java.lang.Double.NaN;

		double min, max;
		for (Face f : faces) {
			min = f.minCoord(0);
			max = f.minCoord(0);
			if (minx != minx || minx > min) minx = min;
			if (maxx != maxx || maxx < max) maxx = max;

			min = f.minCoord(1);
			max = f.minCoord(1);
			if (miny != miny || miny > min) miny = min;
			if (maxy != maxy || maxy < max) maxy = max;

			min = f.minCoord(2);
			max = f.minCoord(2);
			if (minz != minz || minz > min) minz = min;
			if (maxz != maxz || maxz < max) maxz = max;
		}

		return new Vertex[] { new Vertex(minx, miny, minz),
				new Vertex(maxx, maxy, maxz) };
	}

	public class CompareFacesByAxis implements Comparator<Face> {
		public int axis = 0;

		public CompareFacesByAxis(int axis) {
			this.axis = axis % 3;//ensure axis is from 0 to 2
		}

		@Override
		public int compare(Face f1, Face f2) {
			return java.lang.Double.compare(f1.center().get(axis), f2.center()
					.get(axis));
		}
	}
}
