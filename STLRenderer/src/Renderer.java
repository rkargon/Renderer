import java.awt.Color;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Point;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.MouseWheelEvent;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;

import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.Timer;

/**
 * TODO set up editing/manipulation view and render view
 * TODO object selection
 * TODO raytracing
 * 
 * @author raphaelkargon
 * 
 */
public class Renderer extends JPanel {
	ArrayList<STLObject> meshes;
	ArrayList<Lamp> lamps;
	Camera cam;

	boolean wireframe = false;

	//animation
	Timer t;

	public Renderer(File f) {
		STLObject mesh;

		try {
			mesh = new STLObject(f);
			mesh.mat = new Material(Color.white);
		}
		catch (IOException e) {
			mesh = null;
			e.printStackTrace();
		}

		meshes = new ArrayList<STLObject>();
		meshes.add(mesh);
		lamps = new ArrayList<Lamp>();
		lamps.add(new Lamp(5, new Vertex(0, -3, 0), new Color(255, 100, 100)));
		lamps.add(new Lamp(5, new Vertex(0, 0, 3), new Color(100, 255, 100)));
		lamps.add(new Lamp(5, new Vertex(3, 0, 0), new Color(100, 100, 255)));
		lamps.add(new Lamp(5, new Vertex(0, 3, 0), Color.WHITE));

		this.setBackground(Color.WHITE);
		cam = new Camera(-10, 0, 0, 0, 0, 0, 1, 0.1);

		this.addListeners();

		t = new Timer(20, new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				for (Vertex v : meshes.get(0).vertices.values()) {

					v.x += (Math.random() - 0.5) * 0.01;
					v.y += (Math.random() - 0.5) * 0.01;
					v.z += (Math.random() - 0.5) * 0.01;

				}

				repaint();
			}
		});

		//t.start();
	}

	public void addListeners() {
		MouseAdapter m = new MouseAdapter() {

			//keep's track of mouse's original position while dragging mouse
			public Point mousePt = null;

			@Override
			public void mousePressed(MouseEvent e) {
				if (e.getButton() == MouseEvent.BUTTON2) {
					cam.setGlobalRotations(0, 0, 0);
					cam.centerOrigin();
					repaint();
				}

				mousePt = e.getPoint();
			}

			@Override
			public void mouseDragged(MouseEvent e) {
				Point p = e.getPoint();
				cam.rotateLocalY((p.x - mousePt.x) / 100.0);
				cam.rotateLocalX(-(p.y - mousePt.y) / 100.0);
				cam.centerOrigin();

				mousePt = p;
				repaint();
			}

			@Override
			public void mouseWheelMoved(MouseWheelEvent e) {
				double zoomfactor = Math.pow(1.05, e.getWheelRotation());
				cam.center = cam.center.scalarproduct(zoomfactor);
				repaint();
			}
		};

		this.addMouseListener(m);
		this.addMouseMotionListener(m);
		this.addMouseWheelListener(m);

		KeyAdapter k = new KeyAdapter() {
			@Override
			public void keyPressed(KeyEvent e) {
				if (e.getKeyCode() == KeyEvent.VK_UP) {
					if (e.isShiftDown()) cam.rotateLocalZ(-0.1);
					else cam.rotateLocalX(0.1);
				}
				else if (e.getKeyCode() == KeyEvent.VK_DOWN) {
					if (e.isShiftDown()) cam.rotateLocalZ(0.1);
					else cam.rotateLocalX(-0.1);
				}
				else if (e.getKeyCode() == KeyEvent.VK_RIGHT) {
					cam.rotateLocalY(0.1);
				}
				else if (e.getKeyCode() == KeyEvent.VK_LEFT) {
					cam.rotateLocalY(-0.1);
				}
				else if (e.getKeyCode() == KeyEvent.VK_SPACE) {
					cam.setGlobalRotations(0, 0, 0);
				}
				else if (e.getKeyCode()==KeyEvent.VK_Z){
					wireframe=!wireframe;
				}

				cam.centerOrigin();
				repaint();
			}
		};

		this.addKeyListener(k);
	}

	/**
	 * Projects a 3d vertex onto the 2d panel's plane, returning the pixel
	 * coordinates on which to draw the vertex.
	 * 
	 * @param v
	 *            The vertex to be drawn
	 * @return The pixel coordinates of the vertex drawn on the panel. Returns
	 *         <code>(-1,-1)</code> if vertex is too close to camera (based on
	 *         <code>cam.mindist</code>, the camera's clipping distance)
	 */
	public Point projectVertex(Vertex v) {
		Vertex norm = cam.normal(), horiz = cam.horizontal(), vert = cam
				.vertical();
		v = v.subtract(cam.center); //set v to the vector from the camera to the point

		//is point too close to camera? (This is to avoid points on the camera's center, which don't have an angle of incidence)
		if (v.dotproduct(norm) <= cam.mindist) return null;

		Vertex v_nv = v.subtract(horiz.scalarproduct(v.dotproduct(horiz))); //projection onto plane of normal and vert ('vertical plane')
		Vertex v_nh = v.subtract(vert.scalarproduct(v.dotproduct(vert))); //projection onto plane of normal and horiz ('horizontal plane')

		//angles of incidence
		double theta_v = Math.PI / 2
				- Math.acos(v_nv.dotproduct(vert) / v_nv.length()); //vertical angle of incidence
		double theta_h = Math.PI / 2
				- Math.acos(v_nh.dotproduct(horiz) / v_nh.length()); //horizontal angle of incidence

		Point p = new Point((int) ((-theta_h / cam.fov) * getWidth() + getWidth() / 2), (int) ((-theta_v / cam.fov)
				* getWidth() + getHeight() / 2));

		return p;
	}

	public Point projectVertexOrthographic(Vertex v) {
		Vertex horiz = cam.horizontal(), vert = cam.vertical();
		v = v.subtract(cam.center); //set v to the vector from the camera to the point

		//project vertex onto horiz, vert vectors
		Vertex v_h = horiz.scalarproduct(v.dotproduct(horiz));
		Vertex v_v = vert.scalarproduct(v.dotproduct(vert));

		//get x,y coordinates of point on image plane (relative to horiz and vert vectors)
		double x = v_h.dotproduct(horiz);
		double y = v_v.dotproduct(vert);

		Point p = new Point((int) ((-x / cam.fov) * getWidth() + getWidth() / 2), (int) ((-y / cam.fov)
				* getWidth() + getHeight() / 2));

		return p;
	}

	public void paintComponent(Graphics gr) {
		super.paintComponent(gr);
		Graphics2D g2d = (Graphics2D) gr;
		//g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
		g2d.setColor(Color.BLACK);

		for (STLObject mesh : meshes) {
			if (mesh == null) continue;

			//generate set of projected points
			HashMap<Vertex, Point> projectedPoints = new HashMap<Vertex, Point>();
			for (Vertex v : mesh.vertices.values()) {
				projectedPoints.put(v, projectVertex(v));
			}
			List<Face> faces = Arrays.asList(mesh.faces);
			Collections.sort(faces, cam.zsortFaces);

			Point v1, v2, v3;

			if (wireframe) {
				for (Edge e : mesh.edges.values()) {
					v1 = projectedPoints.get(e.v1);
					v2 = projectedPoints.get(e.v2);
					if (v1 == null || v2 == null) continue;
					g2d.drawLine(v1.x, v1.y, v2.x, v2.y);
				}
			}
			else {
				for (Face f : faces) {

					//optional backface culling
					//if(f.normal.dotproduct(cam.normal()) > 0.01) continue;

					v1 = projectedPoints.get(f.vertices[0]);
					v2 = projectedPoints.get(f.vertices[1]);
					v3 = projectedPoints.get(f.vertices[2]);

					if (v1 == null || v2 == null || v3 == null) continue;

					//flat shading, no distance falloff
					/*
					 * For each lamp, increments r,g,b intensity values based on
					 * direction, distance, and intensity of light.
					 * These values are summed for all lamps, and then
					 * multiplied by the color of the material (and clamped to
					 * [0, 255])
					 * Color values are handles separately.
					 */
					double r = 0, g = 0, b = 0;
					for (Lamp l : lamps) {
						Vertex lampvect = l.loc.subtract(f.center());

						//determines angle between normal and lamp.
						double dotprod = lampvect.getUnitVector()
								.dotproduct(f.normal);
						double dstsqr = lampvect.lensquared();

						//assumes well formed solids, no weird normals
						//if solid is well formed, then normals facing away from lamp are in shadow. 
						if (dotprod > 0) {
							if (dstsqr == 0) {
								//if lamp is on the face's center, add full brightness to each color that the lamp emits
								if (l.col.getRed() != 0) r += 1;
								if (l.col.getGreen() != 0) g += 1;
								if (l.col.getBlue() != 0) b += 1;
								continue;
							}
							double intensity = l.intensity * dotprod / dstsqr; //square falloff
							r += (l.col.getRed() / 255.0) * intensity;
							g += (l.col.getGreen() / 255.0) * intensity;
							b += (l.col.getBlue() / 255.0) * intensity;
						}
					}
					Color col = colMultiply(mesh.mat.col, r, g, b);
					g2d.setColor(col);
					//g2d.setColor(new Color(f.hashCode()));

					g2d.fillPolygon(new int[] { v1.x, v2.x, v3.x }, new int[] {
							v1.y, v2.y, v3.y }, 3);
				}
			}
		}

		for (Lamp l : lamps) {
			Point lampcenter = projectVertex(l.loc);
			if (lampcenter == null) continue;

			g2d.setColor(l.col);
			g2d.fillOval(lampcenter.x - 5, lampcenter.y - 5, 10, 10);

			g2d.setColor(Color.BLACK);
			g2d.drawOval(lampcenter.x - 5, lampcenter.y - 5, 10, 10);
			g2d.drawLine(lampcenter.x - 14, lampcenter.y, lampcenter.x - 10, lampcenter.y);// - .
			g2d.drawLine(lampcenter.x, lampcenter.y - 14, lampcenter.x, lampcenter.y - 10);// |.
			g2d.drawLine(lampcenter.x + 14, lampcenter.y, lampcenter.x + 10, lampcenter.y);// . -
			g2d.drawLine(lampcenter.x, lampcenter.y + 14, lampcenter.x, lampcenter.y + 10);// '|
			g2d.drawLine(lampcenter.x - 11, lampcenter.y - 11, lampcenter.x - 8, lampcenter.y - 8);// \.
			g2d.drawLine(lampcenter.x + 11, lampcenter.y + 11, lampcenter.x + 8, lampcenter.y + 8);// '\
			g2d.drawLine(lampcenter.x + 11, lampcenter.y - 11, lampcenter.x + 8, lampcenter.y - 8);// ./
			g2d.drawLine(lampcenter.x - 11, lampcenter.y + 11, lampcenter.x - 8, lampcenter.y + 8);// /`
		}
		//stop=System.nanoTime();
		//System.out.println(stop-start);

	}

	public Color colMultiply(Color c, double r, double g, double b) {
		double col_r = c.getRed() / 255.0;
		double col_g = c.getGreen() / 255.0;
		double col_b = c.getBlue() / 255.0;

		double out_r = Math.min(1, r * col_r); //multiply, clamp
		double out_g = Math.min(1, g * col_g); //multiply, clamp
		double out_b = Math.min(1, b * col_b); //multiply, clamp

		return new Color((int) (out_r * 255), (int) (out_g * 255), (int) (out_b * 255));
	}

	public static void main(String[] args) {
		JFrame f = new JFrame("Raph's ghetto-ass renderer");
		f.setSize(800, 800);
		f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

		Renderer p = new Renderer(new File("/Users/raphaelkargon/Documents/Programming/STL Renderer/suzanneplus.stl"));
		f.getContentPane().add(p);
		p.setFocusable(true);
		p.requestFocusInWindow();

		f.setVisible(true);
	}
}
