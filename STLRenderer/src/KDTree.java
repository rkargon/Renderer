import java.sql.PseudoColumnUsage;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

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

	//TODO save memory by settings objects to null right before recursion?
	public static KDTree buildTree(List<Face> faces) {
		Vertex[] bounds = Face.calcBoundingBox(faces);
		List<FaceWrapper> faceswrapper = FaceWrapper.wrapperList(faces); //wrapper object used to store "classification" values along with faces during tree building
		List<PlanarEvent> events = buildEventList(faceswrapper);
		Collections.sort(events);
		return new KDTree(faceswrapper, events, bounds, 0);
	}

	private KDTree(List<FaceWrapper> faces, List<PlanarEvent> events, Vertex[] bounds, int depth) {
		if (depth == 0) completion = 0;
		this.bounds = bounds;

		if (depth > MAX_DEPTH || faces.size() < MIN_FACES) {
			this.faces = FaceWrapper.toFaceList(faces);
			completion += 1.0 / (1 << depth);
			return;
		}

		double cost = C_TRAVERSAL + C_INTERSECT * faces.size();
		double[] planedata = findOptimalPlane(faces.size(), events, bounds);
		double plane_pos = planedata[0];
		int plane_axis = (int) planedata[1];
		double splitcost = planedata[2];
		double planarside = planedata[3];

		//only split if SAH finds it would reduce cost
		if (cost < splitcost) {
			this.faces = FaceWrapper.toFaceList(faces);
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

		List<FaceWrapper> faces_left = new ArrayList<FaceWrapper>();
		List<FaceWrapper> faces_right = new ArrayList<FaceWrapper>();
		List<FaceWrapper> faces_bothsides = new ArrayList<FaceWrapper>();
		//determine which side of split plane faces fall into 
		for (FaceWrapper fw : faces)
			fw.classification = 3;
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

		//TODO second pass unnecessary? (still classify 'bothsides' faces)
		for (FaceWrapper fw : faces) {
			if (fw.classification == 1) faces_right.add(fw);
			else if (fw.classification == 2) faces_right.add(fw);
			else if (fw.classification == 3) {
				faces_right.add(fw);
				faces_left.add(fw);
				faces_bothsides.add(fw); //used to generate new split events
			}
		}

		//current events grouped into new subvoxels
		List<PlanarEvent> events_leftonly = new ArrayList<PlanarEvent>();
		List<PlanarEvent> events_rightonly = new ArrayList<PlanarEvent>();
		//generate new events from splitting overlapping triangles
		List<PlanarEvent> newevents_left = new ArrayList<PlanarEvent>();
		List<PlanarEvent> newevents_right = new ArrayList<PlanarEvent>();

		for (PlanarEvent e : events) {
			if (e.fw.classification == 1) events_rightonly.add(e);
			else if (e.fw.classification == 2) events_leftonly.add(e);
		}

		for (FaceWrapper fw : faces_bothsides) {
			generateClippedEvents(fw, lowerbounds, newevents_left);
			generateClippedEvents(fw, upperbounds, newevents_right);
		}
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
		for (int i = 0; i < events.size(); i++) {
			PlanarEvent e = events.get(i);
			plane_pos = e.pos;
			plane_axis = e.axis;
			start_count = planar_count = end_count = 0;

			while (i < events.size() && e.axis == plane_axis
					&& e.pos == plane_pos) {
				if (e.type == 0) end_count++;
				else if (e.type == 1) planar_count++;
				else start_count++;
				i++;
				e = events.get(i);
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

	/**
	 * Generates split position events for a triangle clipped to a bounding box
	 * 
	 * @param fw
	 *            The face to be used
	 * @param bounds
	 *            The bounds the face is to be clipped to
	 * @param eventlist
	 *            The event list to which to add events
	 */
	private static void generateClippedEvents(FaceWrapper fw, Vertex[] bounds,
			List<PlanarEvent> events) {
		Face f = fw.f;
		List<PlanarEvent> newevents = new ArrayList<PlanarEvent>();
		List<Vertex> clippedPoints = new ArrayList<Vertex>();
		Vertex vcurr, vnext, v;
		clippedPoints.add(f.vertices[0]);
		clippedPoints.add(f.vertices[1]);
		clippedPoints.add(f.vertices[2]);
		double a, min, max, pcurr, pnext; //ratio used to interpolating vertices

		//for each edge
		for (int i = 0; i <= 2; i++) {
			vcurr = f.vertices[i];
			vnext = f.vertices[(i + 1) % 3];
			//for each axis
			for (int k = 0; k <= 2; k++) {
				pcurr = vcurr.get(k);
				pnext = vnext.get(k);
				min = bounds[0].get(k);
				max = bounds[1].get(k);

				//if two vertices are on opposite sides of bounding plane, clip them to the plane and add that new point to clippedPoints
				if ((pcurr < min && pnext >= min)
						|| (pcurr >= min && pnext < min)) {
					a = (min - pcurr) / (pnext - pcurr);
					clippedPoints.add(Vertex.lerp(vcurr, vnext, a));
				}
				if ((pcurr > max && pnext <= max)
						|| (pcurr <= max && pnext > max)) {
					a = (max - pcurr) / (pnext - pcurr);
					clippedPoints.add(Vertex.lerp(vcurr, vnext, a));
				}
			}
		}

		//TODO maybe add bounds check when each point is added to clippedPoints in the first place?
		for (int i = clippedPoints.size() - 1; i >= 0; i--) {
			v = clippedPoints.get(i);
			//if point is out of bounds, remove it
			if (v.x < bounds[0].x || v.x > bounds[1].x || v.y < bounds[0].y
					|| v.y > bounds[1].y || v.z < bounds[0].z
					|| v.z > bounds[1].z) {
				clippedPoints.remove(i);
			}
		}

		//get bounds of new, clipped points
		double minx, miny, minz, maxx, maxy, maxz;
		minx = miny = minz = Double.POSITIVE_INFINITY;
		maxx = maxy = maxz = Double.NEGATIVE_INFINITY;
		for (Vertex vtmp : clippedPoints) {
			if (vtmp.x < minx) minx = vtmp.x;
			if (vtmp.x > maxx) maxx = vtmp.x;
			if (vtmp.y < miny) miny = vtmp.y;
			if (vtmp.y > maxy) maxy = vtmp.y;
			if (vtmp.z < minz) minz = vtmp.z;
			if (vtmp.z > maxz) maxz = vtmp.z;
		}

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
	 * Wrapper class that pairs a face with a classifcation. This is used to
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
				wf.add(new FaceWrapper(f));
			}
			return wf;
		}
	}

	/** OLD KD TREE IMPLEMENTATION XXX **/

	public KDTree(List<Face> faces) {
		this(faces, Face.calcBoundingBox(faces), 0);
	}

	private KDTree(List<Face> facelist, Vertex[] bounds, int depth) {
		if (depth == 0) completion = 0;
		this.axis = depth % 3;
		this.bounds = bounds;

		if (depth > MAX_DEPTH || facelist.size() < MIN_FACES) {
			this.faces = facelist;
			completion += 1.0 / (1 << depth);
			return;
		}

		//cost of tracing a ray through this node
		double cost = C_TRAVERSAL + C_INTERSECT * facelist.size()
				* Face.surfaceArea(bounds);

		double[] poscost = splitPosition(facelist, axis);
		if (poscost[1] > cost) {
			this.faces = facelist;
			completion += 1.0 / (1 << depth);
			return;
		}
		this.pos = poscost[0];

		ArrayList<Face> lowerfaces = new ArrayList<Face>(), upperfaces = new ArrayList<Face>();
		Vertex[] lowerbounds = new Vertex[] { bounds[0].clone(),
				bounds[1].clone() };
		Vertex[] upperbounds = new Vertex[] { bounds[0].clone(),
				bounds[1].clone() };
		lowerbounds[1].set(axis, pos);
		upperbounds[0].set(axis, pos);
		for (Face f : facelist) {
			//inBounds used instead of inRange, (ie using 3 dimensions to calculate bounds) since faces can overlap divisions, 
			//and it's possible for a face to be in range in one axis but not be in bounds
			if (f.inBounds(lowerbounds)) lowerfaces.add(f);
			if (f.inBounds(upperbounds)) upperfaces.add(f);
		}
		this.lower = new KDTree(lowerfaces, lowerbounds, depth + 1);
		this.upper = new KDTree(upperfaces, upperbounds, depth + 1);
	}

	public double[] splitPosition(List<Face> facelist, int axis) {
		double bestpos = 0;
		double bestcost = Double.NaN;
		double min = bounds[0].get(axis), max = bounds[1].get(axis);

		double costtmp;

		//just check 8 possible locations, much faster than checking face boundaries
		for (double i = min; i < max; i += (max - min) / 8) {
			costtmp = splitCost(facelist, i, axis);
			if (costtmp < bestcost || bestcost != bestcost) {
				bestcost = costtmp;
				bestpos = i;
			}
		}
		return new double[] { bestpos, bestcost };
	}

	public double splitCost(List<Face> facelist, double pos, int axis) {
		double lowerfacecount = 0, upperfacecount = 0;

		Vertex[] lowerbounds = new Vertex[] { bounds[0].clone(),
				bounds[1].clone() };
		Vertex[] upperbounds = new Vertex[] { bounds[0].clone(),
				bounds[1].clone() };
		lowerbounds[1].set(axis, pos);
		upperbounds[0].set(axis, pos);
		for (Face f : facelist) {
			if (f.inBounds(lowerbounds)) lowerfacecount++;
			if (f.inBounds(upperbounds)) upperfacecount++;
		}
		double lowerarea = Face.surfaceArea(lowerbounds);
		double upperarea = Face.surfaceArea(upperbounds);

		double cost = C_TRAVERSAL + C_INTERSECT
				* (lowerarea * lowerfacecount + upperarea * upperfacecount);
		return cost;
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
