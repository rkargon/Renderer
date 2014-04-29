import java.awt.Color;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Point;
import java.awt.Rectangle;
import java.awt.RenderingHints;
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
	World world;

	//0 - wireframe
	//1 - painter's algorithm
	//2 - z-buffering
	int rendermode = 0;

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
		lamps.add(new Lamp(5, new Vertex(0, -3, 0), new Color(255, 100, 100), 2));
		lamps.add(new Lamp(5, new Vertex(0, 0, 3), new Color(100, 255, 100), 2));
		lamps.add(new Lamp(5, new Vertex(3, 0, 0), new Color(100, 100, 255), 2));
		lamps.add(new Lamp(5, new Vertex(0, 3, 0), Color.WHITE, 2));

		this.setBackground(Color.WHITE);
		cam = new Camera(-10, 0, 0, 0, 0, 0, 1, 0.1, 1000.0, false);
		world = new World(Color.WHITE);

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
				else if (e.getKeyCode() == KeyEvent.VK_Z) {
					rendermode++;
					if (rendermode > 2) rendermode = 0;
				}
				else if (e.getKeyCode() == KeyEvent.VK_5) {
					cam.ortho = !cam.ortho;
				}

				cam.centerOrigin();
				repaint();
			}
		};

		this.addKeyListener(k);
	}

	/**
	 * For each lamp, increments r,g,b intensity values based on
	 * direction, distance, and intensity of light.
	 * These values are summed for all lamps, and then
	 * multiplied by the color of the material (and clamped to
	 * [0, 255])
	 * Color values are handles separately.
	 * 
	 * @param v
	 *            the vertex for which to calculate lighting
	 * @param n
	 *            the normal of the vertex, normalized so ||n||=1
	 * @param m
	 *            the material of the face
	 * 
	 * @return the resulting color of the face
	 */
	public Color calcVertexLighting(Vertex v, Vertex n, Material m) {
		//flat shading, no distance falloff

		double r = 0, g = 0, b = 0, spr=0, spg=0, spb=0;
		for (Lamp l : lamps) {
			Vertex lampvect = l.loc.subtract(v);
			Vertex lampvnorm = lampvect.getUnitVector(); //normalize lamp vector
			
			double dotprod = lampvnorm.dotproduct(n);
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

				//n is normalized, so
				//r = d- 2(d . n)n, where r = reflection, d = lampvector (normalized), n = normal
				//Phong relfection: spec intensity = (relf*view)^sphardness * intensity
				Vertex refl = lampvnorm.subtract(n.scalarproduct(2 * n
						.dotproduct(lampvnorm)));
				Vertex view = l.loc.subtract(cam.center).getUnitVector();
				double spec_intensity = l.intensity * m.spintensity * Math.abs(Math.pow(view.dotproduct(refl), m.sphardness));
				double diff_intensity = l.intensity*dotprod*m.diff_intensity;
				
				//calc falloff
				switch (l.falloff) {
				case 2:
					diff_intensity *= 1.0 / dstsqr;
					spec_intensity *= 1.0/dstsqr;
					break;
				case 1:
					diff_intensity *= 1.0 / Math.sqrt(dstsqr);
					spec_intensity *= 1.0/Math.sqrt(dstsqr);
					break;
				case 0:
					break;
				}

				r += (l.col.getRed() / 255.0) * diff_intensity;
				g += (l.col.getGreen() / 255.0) * diff_intensity;
				b += (l.col.getBlue() / 255.0) * diff_intensity;

				spr += (l.col.getRed() / 255.0) * spec_intensity;
				spg += (l.col.getGreen() / 255.0) * spec_intensity;
				spb += (l.col.getBlue() / 255.0) * spec_intensity;
				
			}
		}
		Color col = colMultiply(m.diff_col, r, g, b); //diffuse color
		Color spcol = colMultiply(m.spcol, spr, spg, spb);//specular color
		return colAdd(col, spcol);
	}

	public void paintComponent(Graphics gr) {
		super.paintComponent(gr);
		Graphics2D g2d = (Graphics2D) gr;
		g2d.setColor(Color.BLACK);

		long start = System.nanoTime();

		//create zbuffer for screen, fill with 1's 
		double[][] zbuffer = new double[getHeight()][getWidth()];
		for (double[] row : zbuffer)
			Arrays.fill(row, 1);

		for (STLObject mesh : meshes) {
			if (mesh == null) continue;

			//generate set of projected points
			HashMap<Vertex, Point> projectedPoints = new HashMap<Vertex, Point>();
			for (Vertex v : mesh.vertices.values()) {
				projectedPoints.put(v, cam
						.projectVertex(v, getWidth(), getHeight()));
			}
			//sort faces
			List<Face> faces = Arrays.asList(mesh.faces);
			Collections.sort(faces, cam.zsortFaces);

			Point v1, v2, v3;

			if (rendermode == 0) {
				g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
				for (Edge e : mesh.edges.values()) {
					v1 = projectedPoints.get(e.v1);
					v2 = projectedPoints.get(e.v2);
					if (v1 == null || v2 == null) continue;
					g2d.drawLine(v1.x, v1.y, v2.x, v2.y);
				}
			}
			else if (rendermode == 1) {
				Collections.reverse(faces);
				for (Face f : faces) {

					//optional backface culling
					//if(f.normal.dotproduct(cam.normal()) > 0.01) continue;

					v1 = projectedPoints.get(f.vertices[0]);
					v2 = projectedPoints.get(f.vertices[1]);
					v3 = projectedPoints.get(f.vertices[2]);

					if (v1 == null || v2 == null || v3 == null) continue;
					g2d.setColor(calcVertexLighting(f.center(), f.normal, mesh.mat));
					g2d.fillPolygon(new int[] { v1.x, v2.x, v3.x }, new int[] {
							v1.y, v2.y, v3.y }, 3);
				}
			}
			else if (rendermode == 2) {

				for (Face f : faces) {
					//get triangle corners
					v1 = projectedPoints.get(f.vertices[0]);
					v2 = projectedPoints.get(f.vertices[1]);
					v3 = projectedPoints.get(f.vertices[2]);

					if (v1 == null || v2 == null || v3 == null) continue;

					//zdepth values for vertices. For other pixels, lerp between these values 
					double z1 = (cam.vertexDistance(f.vertices[0]) - cam.mindist)
							/ (cam.maxdist - cam.mindist);
					double z2 = (cam.vertexDistance(f.vertices[1]) - cam.mindist)
							/ (cam.maxdist - cam.mindist);
					double z3 = (cam.vertexDistance(f.vertices[2]) - cam.mindist)
							/ (cam.maxdist - cam.mindist);
					//store difference values. Makes interpolation later on slightly faster 
					double dz21 = z2 - z1;
					double dz31 = z3 - z1;

					//triangle bounding box
					int minX = Math.min(Math.min(v1.x, v2.x), v3.x);
					int maxX = Math.max(Math.max(v1.x, v2.x), v3.x);
					int minY = Math.min(Math.min(v1.y, v2.y), v3.y);
					int maxY = Math.max(Math.max(v1.y, v2.y), v3.y);

					//clip against screen bounds
					minX = Math.max(minX, 0);
					maxX = Math.min(maxX, getWidth() - 1);
					minY = Math.max(minY, 0);
					maxY = Math.min(maxY, getHeight() - 1);

					//triangle edge setup
					int A12 = v1.y - v2.y, B12 = v2.x - v1.x;
					int A23 = v2.y - v3.y, B23 = v3.x - v2.x;
					int A31 = v3.y - v1.y, B31 = v1.x - v3.x;

					//initial barycentric coordinates at corner
					Point p = new Point(minX, minY);
					int w1_row = orient2D(v2, v3, p);
					int w2_row = orient2D(v3, v1, p);
					int w3_row = orient2D(v1, v2, p);

					int w = orient2D(v1, v2, v3);
					if (w == 0) continue;
					int wsgn = Integer.signum(w);

					double z;
					Color col = null;

					//rasterize
					for (p.y = minY; p.y <= maxY; p.y++) {
						//barycentric coordinates at the start of the row
						int w1 = w1_row;
						int w2 = w2_row;
						int w3 = w3_row;

						for (p.x = minX; p.x <= maxX; p.x++) {

							//if w123 have the same sign as w, or are 0, then point is inside triangle
							if ((Integer.signum(w1) == wsgn || w1 == 0)
									&& (Integer.signum(w2) == wsgn || w2 == 0)
									&& (Integer.signum(w3) == wsgn || w3 == 0)) {

								//interpolate z value
								z = z1 + w2 * dz21 / w + w3 * dz31 / w;

								if (z < zbuffer[p.y][p.x]) {
									zbuffer[p.y][p.x] = z;

									//only calculate face color once
									if (col == null)
										col = calcVertexLighting(f.center(), f.normal, mesh.mat);

									g2d.setColor(col);
									g2d.drawLine(p.x, p.y, p.x, p.y);
								}
							}

							//move one pixel to the right
							w1 += A23;
							w2 += A31;
							w3 += A12;
						}

						//move one row down
						w1_row += B23;
						w2_row += B31;
						w3_row += B12;
					}
				}
			}
		}
		//System.out.println(System.nanoTime() - start);

		for (Lamp l : lamps) {
			Point lampcenter = cam
					.projectVertex(l.loc, getWidth(), getHeight());
			if (lampcenter == null) continue;

			//obscure lamp if zbuffer is generated and lamp is occluded
			if (rendermode == 2 && lampcenter.x >= 0
					&& lampcenter.x < getWidth() && lampcenter.y > 0
					&& lampcenter.y < getHeight()) {
				double zlamp = (l.loc.subtract(cam.center).length() - cam.mindist)
						/ (cam.maxdist - cam.mindist);
				if (zbuffer[lampcenter.y][lampcenter.x] < zlamp) continue;
			}

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

	}

	/**
	 * Basically cross product, gives twice signed area of triangle abc
	 * Can also determine if a point p is inside a triangle abc
	 * 
	 * @param a
	 *            Point A
	 * @param b
	 *            Point B
	 * @param c
	 *            Point C
	 * @return (B-A) X (C-A)
	 */
	public int orient2D(Point a, Point b, Point c) {
		return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
	}

	/**
	 * Finds the location of a point <code>p</code> on a face by linearly
	 * interpolating a face's vertices using barycentric coordinates (normalized
	 * to [0,1]).
	 * 
	 * @param f
	 *            The face to be interpolated
	 * @param w1
	 *            the barycentric coordinate corresponding to the first vertex
	 * @param w2
	 *            the barycentric coordinate corresponding to the second vertex
	 * @param w3
	 *            the barycentric coordinate corresponding to the third vertex
	 * @return
	 */
	public Vertex lerpVertex(Face f, double w1, double w2, double w3) {
		double x = w1 * f.vertices[0].x + w2 * f.vertices[1].x + w3
				* f.vertices[2].x;
		double y = w1 * f.vertices[0].y + w2 * f.vertices[1].y + w3
				* f.vertices[2].y;
		double z = w1 * f.vertices[0].z + w2 * f.vertices[1].z + w3
				* f.vertices[2].z;
		return new Vertex(x, y, z);
	}

	public Color colAdd(Color a, Color b){
		return new Color(Math.min(255, a.getRed()+b.getRed()), Math.min(255, a.getGreen()+b.getGreen()), Math.min(255, a.getBlue()+b.getBlue()));
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

		Renderer p = new Renderer(new File("/Users/raphaelkargon/Documents/Programming/STL Renderer/suzanneplusplus.stl"));
		f.getContentPane().add(p);
		p.setFocusable(true);
		p.requestFocusInWindow();

		f.setVisible(true);
	}
}
