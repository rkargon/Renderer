import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.Collections;
import java.util.HashMap;

import javax.swing.Timer;

public class Object3D {
	public MeshVertex[] vertices;
	public Edge[] edges;
	public Face[] faces;
	public String name;

	public Material mat;

	public boolean smooth;


	public Object3D(File f) throws IOException {
		BufferedInputStream in = new BufferedInputStream(new FileInputStream(f));

		//header, 80 bytes
		byte[] headerBytes = new byte[80];
		in.read(headerBytes);

		byte[] face_tmp = new byte[50]; //50 bytes, stores one mesh face

		//read vertex count
		in.read(face_tmp, 0, 4);
		int nfaces = littleEndianInt(face_tmp, 0);
		//set up faces and vertices data 
		HashMap<Vertex, MeshVertex> vertices = new HashMap<Vertex, MeshVertex>();
		HashMap<Edge, Edge> edges = new HashMap<Edge, Edge>();
		faces = new Face[nfaces]; //set face array initial capacity to nfaces for better performance

		//read each triangle
		float x_tmp, y_tmp, z_tmp;
		MeshVertex norm, v1, v2, v3;
		Edge e12, e23, e31; //eij is edge between vertices i and j
		int i = 0; //face counter
		while (in.read(face_tmp) == face_tmp.length) {

			//read normal vector
			x_tmp = littleEndianFloat(face_tmp, 0);
			y_tmp = littleEndianFloat(face_tmp, 4);
			z_tmp = littleEndianFloat(face_tmp, 8);
			norm = new MeshVertex(x_tmp, y_tmp, z_tmp);

			//read vertex 1
			x_tmp = littleEndianFloat(face_tmp, 12);
			y_tmp = littleEndianFloat(face_tmp, 16);
			z_tmp = littleEndianFloat(face_tmp, 20);
			v1 = new MeshVertex(x_tmp, y_tmp, z_tmp);
			//if the vertex does not already exist, add it to the hash. Otherwise, set v1 to the existing vertex
			if (vertices.containsKey(v1)) v1 = vertices.get(v1);
			else vertices.put(v1, v1);

			//read vertex 2
			x_tmp = littleEndianFloat(face_tmp, 24);
			y_tmp = littleEndianFloat(face_tmp, 28);
			z_tmp = littleEndianFloat(face_tmp, 32);
			v2 = new MeshVertex(x_tmp, y_tmp, z_tmp);
			//if the vertex does not already exist, add it to the hash. Otherwise, set v1 to the existing vertex
			if (vertices.containsKey(v2)) v2 = vertices.get(v2);
			else vertices.put(v2, v2);

			//read vertex 3
			x_tmp = littleEndianFloat(face_tmp, 36);
			y_tmp = littleEndianFloat(face_tmp, 40);
			z_tmp = littleEndianFloat(face_tmp, 44);
			v3 = new MeshVertex(x_tmp, y_tmp, z_tmp);
			//if the vertex does not already exist, add it to the hash. Otherwise, set v1 to the existing vertex
			if (vertices.containsKey(v3)) v3 = vertices.get(v3);
			else vertices.put(v3, v3);

			e12 = new Edge(v1, v2);
			e23 = new Edge(v2, v3);
			e31 = new Edge(v3, v1);

			if (edges.containsKey(e12)) e12 = edges.get(e12);
			else edges.put(e12, e12);
			if (edges.containsKey(e23)) e23 = edges.get(e23);
			else edges.put(e23, e23);
			if (edges.containsKey(e31)) e31 = edges.get(e31);
			else edges.put(e31, e31);

			faces[i] = new Face(norm, v1, v2, v3, this);
			i++;

			//ignore attribute byte count, should be 0
		}
		System.out.println(f.getName() + " has " + faces.length + " faces, "
				+ edges.size() + " edges, and " + vertices.size()
				+ " vertices.");
		in.close();

		this.name=f.getName();
		this.vertices = vertices.values().toArray(new MeshVertex[0]);
		this.edges = edges.values().toArray(new Edge[0]);
		smooth = true;
	}

	public Object3D(String name, Face... facearr) {
		this.name=name;
		smooth = true;
		this.faces = facearr;

		HashMap<Vertex, MeshVertex> vertices = new HashMap<Vertex, MeshVertex>();
		HashMap<Edge, Edge> edges = new HashMap<Edge, Edge>();

		MeshVertex v1, v2, v3;
		Edge e12, e23, e31;
		for (Face f : this.faces) {
			f.obj = this;

			v1 = f.vertices[0];
			v2 = f.vertices[1];
			v3 = f.vertices[2];

			if (vertices.containsKey(v1)) v1 = vertices.get(v1);
			else vertices.put(v1, v1);
			if (vertices.containsKey(v2)) v2 = vertices.get(v2);
			else vertices.put(v2, v2);
			if (vertices.containsKey(v3)) v3 = vertices.get(v3);
			else vertices.put(v3, v3);
			
			e12 = new Edge(v1, v2);
			e23 = new Edge(v2, v3);
			e31 = new Edge(v3, v1);

			if (edges.containsKey(e12)) e12 = edges.get(e12);
			else edges.put(e12, e12);
			if (edges.containsKey(e23)) e23 = edges.get(e23);
			else edges.put(e23, e23);
			if (edges.containsKey(e31)) e31 = edges.get(e31);
			else edges.put(e31, e31);
		}

		this.vertices = vertices.values().toArray(new MeshVertex[0]);
		this.edges = edges.values().toArray(new Edge[0]);
	}

	public Vertex center(){
		if(faces.length==0){
			return Vertex.ORIGIN();
		}
		
		Vertex center=Vertex.ORIGIN();
		for(Face f : faces){
			center.add(f.center());
		}
		return center.scalarproduct(1.0/faces.length);
	}
	
	public void move(Vertex dv){
		for(Vertex v : vertices){
			v.setVertex(v.add(dv));
		}
	}
	
	public void moveCenterTo(Vertex v){
		Vertex dv = v.subtract(center());
		move(dv);
	}
	
	/* LITTLE ENDIAN IO */
	public int littleEndianInt(byte[] b, int start) {
		if (start < 0 || start > b.length - 4)
			throw new IllegalArgumentException("start index " + start
					+ " out of bounds");
		return (b[start] & 0xFF) | (b[start + 1] & 0xFF) << 8
				| (b[start + 2] & 0xFF) << 16 | (b[start + 3] & 0xFF) << 24;
	}

	public float littleEndianFloat(byte[] b, int start) {
		return Float.intBitsToFloat(littleEndianInt(b, start));
	}

	public byte[] intToLittleEndianBytes(int i) {

		byte[] b = new byte[4];
		b[0] = (byte) i;
		b[1] = (byte) ((i >>> 8) & 0xFF);
		b[2] = (byte) ((i >>> 16) & 0xFF);
		b[3] = (byte) ((i >>> 24) & 0xFF);

		return b;
	}

	public byte[] floatToLittleEndianBytes(float f) {
		return intToLittleEndianBytes(Float.floatToIntBits(f));
	}

	/* Output */

	public void writeSTL(File f) throws IOException {
		BufferedOutputStream out = new BufferedOutputStream(new FileOutputStream(f));
		byte[] emptyshort = new byte[2];

		out.write(new byte[80]); //header is empty
		out.write(intToLittleEndianBytes(faces.length)); //number of faces

		for (Face face : faces) {
			out.write(floatToLittleEndianBytes((float) face.normal.x));
			out.write(floatToLittleEndianBytes((float) face.normal.y));
			out.write(floatToLittleEndianBytes((float) face.normal.z));

			out.write(floatToLittleEndianBytes((float) face.vertices[0].x));
			out.write(floatToLittleEndianBytes((float) face.vertices[0].y));
			out.write(floatToLittleEndianBytes((float) face.vertices[0].z));

			out.write(floatToLittleEndianBytes((float) face.vertices[1].x));
			out.write(floatToLittleEndianBytes((float) face.vertices[1].y));
			out.write(floatToLittleEndianBytes((float) face.vertices[1].z));

			out.write(floatToLittleEndianBytes((float) face.vertices[2].x));
			out.write(floatToLittleEndianBytes((float) face.vertices[2].y));
			out.write(floatToLittleEndianBytes((float) face.vertices[2].z));

			out.write(emptyshort);
		}
		out.flush();
		out.close();
	}

	public String toString() {
		String s = "";
		s += "STL Mesh: { Name=\"" + name + "\", ";
		s += (faces.length + " faces, " + vertices.length + " vertices.\n");
		for (Face f : faces) {
			s += (f + "\n");
		}
		s += "}";

		return s;
	}

	public static void main(String[] args) {
		try {
			Object3D monkey = new Object3D(new File("/Users/raphaelkargon/Documents/Programming/STL Renderer/suzanne.stl"));
			monkey.writeSTL(new File("/Users/raphaelkargon/Documents/Programming/STL Renderer/suzanne_out.stl"));
		}
		catch (IOException e) {
			//Auto-generated catch block
			e.printStackTrace();
		}
	}
}
