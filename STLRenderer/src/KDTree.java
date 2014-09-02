import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

/**
 * The terms 'upper' and 'lower' are used interchangeably with 'right'
 * and 'left', respectively to indicate different sides of the splitting planes.
 * Apologies to Hebrew- and Japanese-speakers, among
 * others.
 * 
 * @author raphaelkargon
 * 
 */
public class KDTree {
	private static final int MIN_FACES = 4;
	private static final int MAX_DEPTH = 20;
	//used for SAH calculations
	private static final double C_TRAVERSAL = 0; //0 prevents bias against flat, empty cells, which can be advantageous in architectural scenes, etc.
	private static final double C_INTERSECT = 1;

	public KDTree lower, upper;
	public List<Face> faces; //stores faces in leaf nodes
	public double pos;
	public int axis;
	public Vertex[] bounds;

	public static double completion = 0;

	public static KDTree buildTree(List<Face> faces) {
		Vertex[] bounds = Face.calcBoundingBox(faces);
		List<FaceWrapper> faceswrapper = FaceWrapper.wrapperList(faces); //wrapper object used to store "classification" values along with faces during tree building
		List<PlanarEvent> events = buildEventList(faceswrapper);
		Collections.sort(events);
		return new KDTree(events, bounds, 0);
	}

	private KDTree(List<PlanarEvent> events, Vertex[] bounds, int depth) {
		if (depth == 0) completion = 0;
		this.bounds = bounds;
		int nfaces = 0;
		//face classifications are initially not 0, and then are not 0 when passed from parent nodes.
		for (PlanarEvent e : events) {

			if (e.fw.classification != 0) {
				nfaces++;
				e.fw.classification = 0;
			}
		}
		//exit if max recursion depth is reached or too few faces are remaining
		if (depth > MAX_DEPTH || nfaces < MIN_FACES) {
			List<Face> faces = new ArrayList<Face>();
			//all face classifications are zero due to first pass
			for (PlanarEvent e : events) {
				if (e.fw.classification == 0) {
					faces.add(e.fw.f);
					e.fw.classification = 1;
				}
			}
			this.faces = faces;
			completion += 1.0 / (1 << depth);
			//System.out.println("depth: " + depth + " nfaces: " + this.faces.size());
			return;

		}

		double cost = C_TRAVERSAL + C_INTERSECT * nfaces;
		double[] planedata = findOptimalPlane(nfaces, events, bounds);
		double plane_pos = planedata[0];
		int plane_axis = (int) planedata[1];
		double splitcost = planedata[2];
		double planarside = planedata[3];

		//only split if SAH finds it would reduce cost
		if (cost <= splitcost) {
			List<Face> faces = new ArrayList<Face>();
			for (PlanarEvent e : events) {
				if (e.fw.classification == 0) {
					faces.add(e.fw.f);
					e.fw.classification = 1;
				}
			}
			this.faces = faces;
			completion += 1.0 / (1 << depth);
			return;
		}

		this.pos = plane_pos;
		this.axis = plane_axis;
		Vertex newmin = bounds[0].clone();
		newmin.set(axis, pos);
		Vertex newmax = bounds[1].clone();
		newmax.set(axis, pos);
		Vertex[] lowerbounds = new Vertex[] { bounds[0], newmax };
		Vertex[] upperbounds = new Vertex[] { newmin, bounds[1] };

		//current events grouped into new subvoxels
		List<PlanarEvent> events_left = new ArrayList<PlanarEvent>();
		List<PlanarEvent> events_right = new ArrayList<PlanarEvent>();
		//generate new events from splitting overlapping triangles
		List<FaceWrapper> faces_bothsides = new ArrayList<FaceWrapper>();
		List<PlanarEvent> newevents_left = new ArrayList<PlanarEvent>();
		List<PlanarEvent> newevents_right = new ArrayList<PlanarEvent>();

		//determine which side of split plane faces fall into 
		for (PlanarEvent e : events)
			e.fw.classification = 3;
		for (PlanarEvent e : events) {
			if (e.axis == axis) {
				if (e.type == 0 && e.pos <= pos) {
					e.fw.classification = 2; //left only
				}
				else if (e.type == 2 && e.pos >= pos) {
					e.fw.classification = 1; //right only
				}
				else if (e.type == 1) {
					if (e.pos == pos) e.fw.classification = (planarside == 0) ? 2
							: 1;
					else if (e.pos > pos) e.fw.classification = 1;
					else e.fw.classification = 2;
				}
			}
		}

		//second pass necessary because only events of one axis are looked at. 
		for (PlanarEvent e : events) {
			if (e.fw.classification == 1) events_right.add(e);
			else if (e.fw.classification == 2) events_left.add(e);
			else if (e.fw.classification == 3) {
				faces_bothsides.add(e.fw);
				e.fw.classification = 4;//prevents face from being added twice.
			}
		}
		events = null;

		for (FaceWrapper fw : faces_bothsides) {
			generateClippedEvents(fw, pos, axis, bounds, newevents_left, newevents_right);
		}

		faces_bothsides = null;

		Collections.sort(newevents_left);
		Collections.sort(newevents_right);
		events_left = PlanarEvent.mergeLists(events_left, newevents_left);
		events_right = PlanarEvent.mergeLists(events_right, newevents_right);
		newevents_left = newevents_right = null;
		this.lower = new KDTree(events_left, lowerbounds, depth + 1);
		events_left = null;
		this.upper = new KDTree(events_right, upperbounds, depth + 1);
	}

	private static List<PlanarEvent> buildEventList(List<FaceWrapper> faces) {
		List<PlanarEvent> events = new ArrayList<PlanarEvent>();

		double minx, maxx, miny, maxy, minz, maxz;
		for (FaceWrapper fw : faces) {
			minx = fw.f.minCoord(0);
			maxx = fw.f.maxCoord(0);
			miny = fw.f.minCoord(1);
			maxy = fw.f.maxCoord(1);
			minz = fw.f.minCoord(2);
			maxz = fw.f.maxCoord(2);

			if (minx == maxx) events.add(new PlanarEvent(fw, minx, 0, 1));
			else {
				events.add(new PlanarEvent(fw, minx, 0, 2));
				events.add(new PlanarEvent(fw, maxx, 0, 0));
			}
			if (miny == maxy) events.add(new PlanarEvent(fw, miny, 1, 1));
			else {
				events.add(new PlanarEvent(fw, miny, 1, 2));
				events.add(new PlanarEvent(fw, maxy, 1, 0));
			}
			if (minz == maxz) events.add(new PlanarEvent(fw, minz, 2, 1));
			else {
				events.add(new PlanarEvent(fw, minz, 2, 2));
				events.add(new PlanarEvent(fw, maxz, 2, 0));
			}
		}
		return events;
	}

	/**
	 * Finds optimal splitting plane in O(N) time.
	 * From: http://www.eng.utah.edu/~cs6965/papers/kdtree.pdf
	 * 
	 * @param nfaces
	 *            The number of faces in the current bounding box
	 * @param events
	 *            A list of splitting plane candidates
	 * @param bounds
	 *            THe current bounding box
	 * @return
	 */
	private static double[] findOptimalPlane(int nfaces,
			List<PlanarEvent> events, Vertex[] bounds) {
		double area = Face.surfaceArea(bounds);

		int[] Nleft = new int[3], Nplanar = new int[3], Nright = new int[] {
				nfaces, nfaces, nfaces };
		double pos_best = 0, axis_best = 0, cost_best = Double.POSITIVE_INFINITY, p_side = 0;
		//temporary variables
		int plane_axis, start_count, planar_count, end_count;
		double plane_pos, planedata[];//planedata = {cost, pside}
		for (int i = 0; i < events.size();) {
			PlanarEvent e = events.get(i);
			plane_pos = e.pos;
			plane_axis = e.axis;
			start_count = planar_count = end_count = 0;

			while (i < events.size() && e.axis == plane_axis
					&& e.pos == plane_pos) {
				e = events.get(i);
				if (e.type == 0) end_count++;
				else if (e.type == 1) planar_count++;
				else start_count++;
				i++;
			}

			Nplanar[plane_axis] = planar_count;
			Nright[plane_axis] -= planar_count;
			Nright[plane_axis] -= end_count;

			planedata = SAH(bounds, plane_pos, plane_axis, area, Nleft[plane_axis], Nplanar[plane_axis], Nright[plane_axis]);
			if (planedata[0] < cost_best) {
				cost_best = planedata[0];
				pos_best = plane_pos;
				axis_best = plane_axis;
				p_side = planedata[1];
			}
			Nleft[plane_axis] += start_count;
			Nleft[plane_axis] += planar_count;
			Nplanar[plane_axis] = 0;
		}
		return new double[] { pos_best, axis_best, cost_best, p_side };
	}

	/**
	 * Computes approximate cost of traversing a ray through a voxel with a
	 * given split plane, using the Surface Area Heuristic
	 * 
	 * @param bounds
	 *            The bounds of the whole voxel
	 * @param position
	 *            The position of the splitting plane
	 * @param axis
	 *            The axis of the splitting plane
	 * @param Nl
	 *            The number of faces to the left (below) the splitting plane
	 * @param Np
	 *            The number of faces exactly (flat) on the splitting plane
	 * @param Nr
	 *            The number of faces to the right (above) the splitting plane
	 * @return a double array {cost, pside}. pside = 0 or 1, depending if planar
	 *         faces should be grouped in the left or right voxel, respectively
	 */
	private static double[] SAH(Vertex[] bounds, double position, int axis,
			double area, int Nl, int Np, int Nr) {
		Vertex newmin = bounds[0].clone();
		newmin.set(axis, position);
		Vertex newmax = bounds[1].clone();
		newmax.set(axis, position);

		if (area == 0) return new double[] { C_TRAVERSAL, 0 };
		double lowerarea = Face.surfaceArea(bounds[0], newmax);
		double upperarea = Face.surfaceArea(newmin, bounds[1]);
		lowerarea /= area;
		upperarea /= area;

		double costleftbias = lowerarea * (Nl + Np) + upperarea * Nr;
		double costrightbias = lowerarea * Nl + upperarea * (Np + Nr);
		if (costleftbias <= costrightbias) return new double[] {
				costleftbias * C_INTERSECT + C_TRAVERSAL, 0 };

		else return new double[] { costrightbias * C_INTERSECT + C_TRAVERSAL, 1 };
	}

	private static void generateClippedEvents(FaceWrapper fw, double pos,
			int axis, Vertex[] bounds, List<PlanarEvent> lowerevents,
			List<PlanarEvent> upperevents) {

		//0 if vertex is below position
		int sideflags = (fw.f.vertices[0].get(axis) < pos ? 0 : 1)
				| (fw.f.vertices[1].get(axis) < pos ? 0 : 2)
				| (fw.f.vertices[2].get(axis) < pos ? 0 : 4);

		//this checks if sideflags is a power of two, ie there is only one vertex on the right.
		boolean oneUpperSide = (sideflags & (sideflags - 1)) == 0;

		//For whichever side has only one vertex, check with vertex it is
		//Assign the 'lonely' vertex to v1
		Vertex v1 = null, v2 = null, v3 = null;
		switch (oneUpperSide ? sideflags : (~sideflags) & 7) {
		case 1:
			v1 = fw.f.vertices[0];
			v2 = fw.f.vertices[1];
			v3 = fw.f.vertices[2];
			break;
		case 2:
			v1 = fw.f.vertices[1];
			v2 = fw.f.vertices[0];
			v3 = fw.f.vertices[2];
			break;
		case 4:
			v1 = fw.f.vertices[2];
			v2 = fw.f.vertices[0];
			v3 = fw.f.vertices[1];
			break;
		}

		Vertex v12cut, v13cut; //points on split plane along triangle edges
		double v1pos = v1.get(axis);
		double r12 = (pos - v1pos) / (v2.get(axis) - v1pos);
		double r13 = (pos - v1pos) / (v3.get(axis) - v1pos);
		v12cut = Vertex.lerp(v1, v2, r12);
		v13cut = Vertex.lerp(v1, v3, r13);

		int ax1 = (axis + 2) % 3, ax2 = (axis + 1) % 3;
		//clamp values of intersecting points to the bounding box
		v12cut.set(ax1, MathUtils.clamp(v12cut.get(ax1), bounds[0].get(ax1), bounds[1]
				.get(ax1)));
		v12cut.set(ax2, MathUtils.clamp(v12cut.get(ax2), bounds[0].get(ax2), bounds[1]
				.get(ax2)));
		v13cut.set(ax1, MathUtils.clamp(v13cut.get(ax1), bounds[0].get(ax1), bounds[1]
				.get(ax1)));
		v13cut.set(ax2, MathUtils.clamp(v13cut.get(ax2), bounds[0].get(ax2), bounds[1]
				.get(ax2)));

		Vertex[] lowerbounds, upperbounds;
		if (oneUpperSide) {
			upperbounds = Vertex.listBounds(v1, v12cut, v13cut);
			lowerbounds = Vertex.listBounds(v2, v3, v12cut, v13cut);
		}
		else {
			lowerbounds = Vertex.listBounds(v1, v12cut, v13cut);
			upperbounds = Vertex.listBounds(v2, v3, v12cut, v13cut);

		}

		//intersecting triangle bounding boxes with bounding box of volume
		lowerbounds = Vertex.intersectBoundingBoxes(lowerbounds, bounds);
		upperbounds = Vertex.intersectBoundingBoxes(upperbounds, bounds);

		//lower events
		if (lowerbounds[0].x == lowerbounds[1].x) {
			lowerevents.add(new PlanarEvent(fw, lowerbounds[0].x, 0, 1));
		}
		else {
			lowerevents.add(new PlanarEvent(fw, lowerbounds[0].x, 0, 2));
			lowerevents.add(new PlanarEvent(fw, lowerbounds[1].x, 0, 0));
		}
		if (lowerbounds[0].y == lowerbounds[1].y) {
			lowerevents.add(new PlanarEvent(fw, lowerbounds[0].y, 1, 1));
		}
		else {
			lowerevents.add(new PlanarEvent(fw, lowerbounds[0].y, 1, 2));
			lowerevents.add(new PlanarEvent(fw, lowerbounds[1].y, 1, 0));
		}
		if (lowerbounds[0].z == lowerbounds[1].z) {
			lowerevents.add(new PlanarEvent(fw, lowerbounds[0].z, 2, 1));
		}
		else {
			lowerevents.add(new PlanarEvent(fw, lowerbounds[0].z, 2, 2));
			lowerevents.add(new PlanarEvent(fw, lowerbounds[1].z, 2, 0));
		}

		//upper events
		if (upperbounds[0].x == upperbounds[1].x) {
			upperevents.add(new PlanarEvent(fw, upperbounds[0].x, 0, 1));
		}
		else {
			upperevents.add(new PlanarEvent(fw, upperbounds[0].x, 0, 2));
			upperevents.add(new PlanarEvent(fw, upperbounds[1].x, 0, 0));
		}
		if (upperbounds[0].y == upperbounds[1].y) {
			upperevents.add(new PlanarEvent(fw, upperbounds[0].y, 1, 1));
		}
		else {
			upperevents.add(new PlanarEvent(fw, upperbounds[0].y, 1, 2));
			upperevents.add(new PlanarEvent(fw, upperbounds[1].y, 1, 0));
		}
		if (upperbounds[0].z == upperbounds[1].z) {
			upperevents.add(new PlanarEvent(fw, upperbounds[0].z, 2, 1));
		}
		else {
			upperevents.add(new PlanarEvent(fw, upperbounds[0].z, 2, 2));
			upperevents.add(new PlanarEvent(fw, upperbounds[1].z, 2, 0));
		}
	}

	/**
	 * Represents a "planar event", used in sorting triangles for SAH. It stores
	 * a reference to a face, the position of the event, the dimension the event
	 * corresponds to, and whether the triangle starts, ends or is planar
	 */
	private static class PlanarEvent implements Comparable<PlanarEvent> {
		public FaceWrapper fw;
		public double pos;

		// 0 - x
		// 1 - y
		// 2 - z
		public int axis;

		// 0 - f ends at p
		// 1 - f is planar on p
		// 2 - f starts at p
		public int type;

		public PlanarEvent(FaceWrapper fw, double p, int k, int type) {
			this.fw = fw;
			this.pos = p;
			this.axis = k;
			this.type = type;
		}

		/**
		 * Combines two sorted lists
		 * 
		 * @param events_leftonly
		 * @param newevents_left
		 */
		public static List<PlanarEvent> mergeLists(List<PlanarEvent> l1,
				List<PlanarEvent> l2) {
			int l1size = l1.size(), l2size = l2.size();
			List<PlanarEvent> outlist = new ArrayList<PlanarEvent>(l1size
					+ l2size);
			int p1 = 0, p2 = 0;
			while (p1 < l1size && p2 < l2size) {
				if (l2.get(p2).compareTo(l1.get(p1)) < 0) {
					outlist.add(l2.get(p2));
					p2++;
				}
				else {
					outlist.add(l1.get(p1));
					p1++;
				}
			}
			if (p1 >= l1size) {
				for (; p2 < l2size; p2++)
					outlist.add(l2.get(p2));
			}
			else if (p2 >= l2size) {
				for (; p1 < l1size; p1++)
					outlist.add(l1.get(p1));
			}
			return outlist;
		}

		@Override
		public int compareTo(PlanarEvent e) {
			int c = Double.compare(pos, e.pos);
			if (c == 0) {
				c = Integer.compare(axis, e.axis);
				if (c == 0) {
					c = Integer.compare(type, e.type);
				}
			}
			return c;
		}
	}

	/**
	 * Wrapper class that pairs a face with a classification. This is used to
	 * group faces into the left/right sides of a splitting plane during KD Tree
	 * building.
	 */
	private static class FaceWrapper {
		public Face f;

		//last two bits of int correspond to left, right subset during classification.
		//00 = 0 = none
		//01 = 1 = right
		//10 = 2 = left
		//11 = 3 = both
		public int classification;

		public FaceWrapper(Face f, int classification) {
			this.f = f;
			this.classification = classification;
		}

		public FaceWrapper(Face f) {
			this.f = f;
		}

		public static List<Face> toFaceList(List<FaceWrapper> wrappedfaces) {
			List<Face> faces = new ArrayList<Face>();
			for (FaceWrapper w : wrappedfaces)
				faces.add(w.f);
			return faces;
		}

		public static List<FaceWrapper> wrapperList(List<Face> faces) {
			List<FaceWrapper> wf = new ArrayList<FaceWrapper>();
			for (Face f : faces) {
				wf.add(new FaceWrapper(f, 1));
			}
			return wf;
		}
	}

	public List<Edge> wireframe() {
		return wireframe(true);
	}

	private List<Edge> wireframe(boolean isRoot) {
		List<Edge> edges = new ArrayList<Edge>();

		if (isRoot) {
			//create corner vertices
			Vertex xyz = bounds[0];
			Vertex XYZ = bounds[1];
			Vertex xyZ = new Vertex(xyz.x, xyz.y, XYZ.z);
			Vertex xYz = new Vertex(xyz.x, XYZ.y, xyz.z);
			Vertex xYZ = new Vertex(xyz.x, XYZ.y, XYZ.z);
			Vertex Xyz = new Vertex(XYZ.x, xyz.y, xyz.z);
			Vertex XyZ = new Vertex(XYZ.x, xyz.y, XYZ.z);
			Vertex XYz = new Vertex(XYZ.x, XYZ.y, xyz.z);

			//from bottom corner
			edges.add(new Edge(xyz, xyZ));
			edges.add(new Edge(xyz, xYz));
			edges.add(new Edge(xyz, Xyz));

			//from top corner
			edges.add(new Edge(XYZ, XYz));
			edges.add(new Edge(XYZ, XyZ));
			edges.add(new Edge(XYZ, xYZ));

			//connect remaining edges
			edges.add(new Edge(xyZ, xYZ));
			edges.add(new Edge(xyZ, XyZ));
			edges.add(new Edge(xYz, XYz));
			edges.add(new Edge(xYz, xYZ));
			edges.add(new Edge(Xyz, XYz));
			edges.add(new Edge(Xyz, XyZ));
		}

		if (faces == null) {
			//add edges for split plane
			switch (axis) {
			case 0:
				edges.add(new Edge(new Vertex(pos, bounds[0].y, bounds[0].z), new Vertex(pos, bounds[0].y, bounds[1].z)));
				edges.add(new Edge(new Vertex(pos, bounds[0].y, bounds[0].z), new Vertex(pos, bounds[1].y, bounds[0].z)));
				edges.add(new Edge(new Vertex(pos, bounds[0].y, bounds[1].z), new Vertex(pos, bounds[1].y, bounds[1].z)));
				edges.add(new Edge(new Vertex(pos, bounds[1].y, bounds[0].z), new Vertex(pos, bounds[1].y, bounds[1].z)));
				break;
			case 1:
				edges.add(new Edge(new Vertex(bounds[0].x, pos, bounds[0].z), new Vertex(bounds[0].x, pos, bounds[1].z)));
				edges.add(new Edge(new Vertex(bounds[0].x, pos, bounds[0].z), new Vertex(bounds[1].x, pos, bounds[0].z)));
				edges.add(new Edge(new Vertex(bounds[0].x, pos, bounds[1].z), new Vertex(bounds[1].x, pos, bounds[1].z)));
				edges.add(new Edge(new Vertex(bounds[1].x, pos, bounds[0].z), new Vertex(bounds[1].x, pos, bounds[1].z)));
				break;
			case 2:
				edges.add(new Edge(new Vertex(bounds[0].x, bounds[0].y, pos), new Vertex(bounds[0].x, bounds[1].y, pos)));
				edges.add(new Edge(new Vertex(bounds[0].x, bounds[0].y, pos), new Vertex(bounds[1].x, bounds[0].y, pos)));
				edges.add(new Edge(new Vertex(bounds[1].x, bounds[0].y, pos), new Vertex(bounds[1].x, bounds[1].y, pos)));
				edges.add(new Edge(new Vertex(bounds[0].x, bounds[1].y, pos), new Vertex(bounds[1].x, bounds[1].y, pos)));
				break;
			}

			edges.addAll(lower.wireframe(false));
			edges.addAll(upper.wireframe(false));
		}
		return edges;
	}
}
