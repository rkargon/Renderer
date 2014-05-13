import java.util.ArrayList;
import java.util.List;

public class KDTree {
	//TODO implement SAH

	private static final int MIN_FACES = 4;
	private static final int MAX_DEPTH = 20;
	//used for SAH calculations
	//values derived from code samle from http://www.flipcode.com/archives/Raytracing_Topics_Techniques-Part_7_Kd-Trees_and_More_Speed.shtml
	private static final double C_TRAVERSAL = 0.3;
	private static final double C_INTERSECT = 1;

	public KDTree lower, upper;
	public List<Face> faces; //stores faces in leaf nodes
	public double loc;
	public int axis;
	public Vertex[] bounds;

	public KDTree(List<Face> faces) {
		this(faces, Face.calcBoundingBox(faces), 0);
	}

	private KDTree(List<Face> facelist, Vertex[] bounds, int depth) {
		this.axis = depth % 3;
		this.bounds = bounds;

		if (depth > MAX_DEPTH || facelist.size() < MIN_FACES) {
			this.faces = facelist;
			return;
		}

		//cost of tracing a ray through this node
		double cost = C_TRAVERSAL + C_INTERSECT * facelist.size()
				* Face.surfaceArea(bounds);

		//this.loc=optimalsplitposition();
		double[] poscost = splitPosition(facelist, axis);
		if (poscost[1] > cost) {
			this.faces = facelist;
			return;
		}
		this.loc = poscost[0];

		ArrayList<Face> lowerfaces = new ArrayList<Face>(), upperfaces = new ArrayList<Face>();
		Vertex[] lowerbounds = new Vertex[] { bounds[0].clone(),
				bounds[1].clone() };
		Vertex[] upperbounds = new Vertex[] { bounds[0].clone(),
				bounds[1].clone() };
		lowerbounds[1].set(axis, loc);
		upperbounds[0].set(axis, loc);
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
		//TODO optimize
		// eg. faces already arranged when checking cost. Use this in constructor to build children
		// Generate split positions and sort. If a face is already below the current position, it will remain below
		//TODO finish
		double bestpos = 0;
		double bestcost = Double.NaN;
		double min = bounds[0].get(axis), max = bounds[1].get(axis);

		double costtmp;
		//		for (Face f : facelist) {
		//			postmp = f.minCoord(axis);
		//			//check that position is in bounds
		//			if (postmp > bounds[0].get(axis) && postmp < bounds[1].get(axis)) {
		//				//get cost of this split
		//				costtmp = splitCost(facelist, postmp, axis);
		//				//check if cost is best so far
		//				if (costtmp < bestcost || bestcost != bestcost) {
		//					bestcost = costtmp;
		//					bestpos = postmp;
		//				}
		//			}
		//
		//			postmp = f.maxCoord(axis);
		//			if (postmp > bounds[0].get(axis) && postmp < bounds[1].get(axis)) {
		//				costtmp = splitCost(facelist, postmp, axis);
		//				if (costtmp < bestcost || bestcost != bestcost) {
		//					bestcost = costtmp;
		//					bestpos = postmp;
		//				}
		//			}
		//		}

		//just check 8 possible locations, much faster than checking face boundaries
		for (double i = min; i < max; i += (max - min) / 8) {
			costtmp = splitCost(facelist, i, axis);
			if (costtmp < bestcost || bestcost != bestcost) {
				bestcost = costtmp;
				bestpos = i;
			}
		}
		return new double[] { bestpos, bestcost};
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
			switch(axis){
			case 0:
				edges.add(new Edge(new Vertex(loc, bounds[0].y, bounds[0].z), new Vertex(loc, bounds[0].y, bounds[1].z)));
				edges.add(new Edge(new Vertex(loc, bounds[0].y, bounds[0].z), new Vertex(loc, bounds[1].y, bounds[0].z)));
				edges.add(new Edge(new Vertex(loc, bounds[0].y, bounds[1].z), new Vertex(loc, bounds[1].y, bounds[1].z)));
				edges.add(new Edge(new Vertex(loc, bounds[1].y, bounds[0].z), new Vertex(loc, bounds[1].y, bounds[1].z)));
				break;
			case 1:
				edges.add(new Edge(new Vertex(bounds[0].x, loc, bounds[0].z), new Vertex(bounds[0].x, loc, bounds[1].z)));
				edges.add(new Edge(new Vertex(bounds[0].x, loc, bounds[0].z), new Vertex(bounds[1].x, loc, bounds[0].z)));
				edges.add(new Edge(new Vertex(bounds[0].x, loc, bounds[1].z), new Vertex(bounds[1].x, loc, bounds[1].z)));
				edges.add(new Edge(new Vertex(bounds[1].x, loc, bounds[0].z), new Vertex(bounds[1].x, loc, bounds[1].z)));
				break;
			case 2:
				edges.add(new Edge(new Vertex(bounds[0].x, bounds[0].y, loc), new Vertex(bounds[0].x, bounds[1].y, loc)));
				edges.add(new Edge(new Vertex(bounds[0].x, bounds[0].y, loc), new Vertex(bounds[1].x, bounds[0].y, loc)));
				edges.add(new Edge(new Vertex(bounds[1].x, bounds[0].y, loc), new Vertex(bounds[1].x, bounds[1].y, loc)));
				edges.add(new Edge(new Vertex(bounds[0].x, bounds[1].y, loc), new Vertex(bounds[1].x, bounds[1].y, loc)));
				break;
			}
			
			edges.addAll(lower.wireframe(false));
			edges.addAll(upper.wireframe(false));
		}
		return edges;
	}
}
