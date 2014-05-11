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

	public KDTree(List<Face> faces) {
		this(faces, Face.calcBoundingBox(faces), 0);
	}

	public KDTree(List<Face> facelist, Vertex[] bounds, int depth) {
		if (depth > MAX_DEPTH || facelist.size() < MIN_FACES) {
			this.faces = facelist;
			return;
		}
		this.axis = depth % 3;

		//cost of tracing a ray through this node
		double cost = C_TRAVERSAL + C_INTERSECT * facelist.size()
				* Face.surfaceArea(bounds);

		//TODO this.loc=optimalsplitposition();
		this.loc = splitPosition(facelist, axis);
		//TODO for each face
		//TODO if face is in lower region, add to lower array
		//TODO if face is in upper region, add to upper array
		//this.lower = new KDTree(lowerfaces, depth+1);
		//this.upper = new KDTree(upperfaces, depth+1);
	}

	public double splitPosition(List<Face> facelist, int axis) {
		//optimize
		double bestpos = Double.NaN;
		double bestcost = Double.NaN;

		for (Face f : facelist){
			
		}
		
		return 0;
	}
}
