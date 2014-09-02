import java.awt.BasicStroke;
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
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;

import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.Timer;

/**
 * 
 * @author raphaelkargon
 * 
 */
public class Renderer extends JPanel {

	/* FIELDS */

	ArrayList<Object3D> meshes;
	ArrayList<Lamp> lamps;
	Camera cam;
	World world;

	Object3D selected = null;
	List<Face> faces;
	JLabel status;

	boolean isManipulating = false;
	int manipMode = 0; //0-grab, 1-scale, 2-rotate
	Point oldmouseloc;
	Vertex oldcenter;

	//0 - wireframe
	//1 - painter's algorithm
	//2 - z-buffering
	int rendermode = 0;
	boolean rendering;
	final int RAY_DEPTH = 10;
	final int TILESIZE = 32;
	BufferedImage render;
	KDTree renderkdt;
	List<Edge> kdedges;
	Thread renderthread = new Thread() {

		public void run() {
			raytraceRender();
			this.interrupt();
		}
	};
	//real time update of rendering process
	Timer t = new Timer(100, new ActionListener() {
		public void actionPerformed(ActionEvent e) {
			if (!rendering) t.stop();
			repaint();
		}
	});
	private long raysTraced = 0;//used for performance stats

	/* METHODS */

	public Renderer() {
		meshes = new ArrayList<Object3D>();

		Object3D mesh;

		File f = new File("/Users/raphaelkargon/Documents/Programming/STL Renderer/suzanneplusplus.stl");
		try {
			mesh = new Object3D(f);
			mesh.mat = new Material(Color.WHITE);
			meshes.add(mesh);
		}
		catch (IOException e) {
			e.printStackTrace();
		}

		/* Compile list of faces. Add meshes before this point */
		faces = new ArrayList<Face>();
		for (Object3D meshtmp : meshes) {
			Collections.addAll(faces, meshtmp.faces);
		}

		lamps = new ArrayList<Lamp>();
		lamps.add(new Lamp(12, new Vertex(-3, -3, 0), new Color((int) ((Math
				.random() - 0.5) * Integer.MAX_VALUE)), 2));
		lamps.add(new Lamp(12, new Vertex(3, -3, 0), new Color((int) ((Math
				.random() - 0.5) * Integer.MAX_VALUE)), 2));
		lamps.add(new Lamp(12, new Vertex(0, 4.2, 0), new Color((int) ((Math
				.random() - 0.5) * Integer.MAX_VALUE)), 2));

		cam = new Camera();
		world = new World(Color.WHITE, new Color(100, 100, 255));
		this.setBackground(Color.WHITE);
		this.setLayout(null);
		status = new JLabel();
		status.setLocation(10, 0);
		status.setSize(300, 30);
		this.add(status);
		status.setVisible(true);

		this.addListeners();
		
		renderkdt = KDTree.buildTree(faces);
		long begin = System.nanoTime();
		for(int i=1; i<=1_000_000; i++){
			rayIntersect(renderkdt, Vertex.ORIGIN(), new Vertex(Math.random(),  Math.random(), Math.random()), false, null);
			//rayIntersect(renderkdt, Vertex.ORIGIN(), new Vertex(1,1,1), false, null);
			//rayIntersect(faces, Vertex.ORIGIN(), new Vertex(1,1,1), false, null);
			//faces.get(0).intersectRayTriangle(Vertex.ORIGIN(), new Vertex(1,1,1), null);
		}
		long end = System.nanoTime();
		System.out.println("1,000,000 kdt intersects: " + (double)(end-begin)/1_000_000_000);
		kdedges = renderkdt.wireframe();
		System.out.println(kdedges.size() + " edges in kdt wireframe.");
		rayIntersect(renderkdt, Vertex.ORIGIN(), new Vertex(1,1,1), false, null);
	}

	public void addListeners() {
		MouseAdapter m = new MouseAdapter() {

			//keep's track of mouse's original position while dragging mouse
			public Point mousePt = null;

			@Override
			public void mousePressed(MouseEvent e) {
				if (isManipulating) {
					isManipulating = false;
				}
				if (e.getButton() == MouseEvent.BUTTON2) {
					//if camera is already centered on origin, reset rotation
					if (cam.focus.lensquared() == 0) {
						cam.setGlobalRotation(0, 0, 0);
					}
					else {
						cam.focus = Vertex.ORIGIN();
					}
					cam.centerFocus();
				}
				else if (e.getButton() == MouseEvent.BUTTON3) {
					Vertex[] ray = cam
							.castRay(e.getX(), e.getY(), getWidth(), getHeight());
					Face f = rayIntersect(faces, ray[0], ray[1], false, null);
					if (f == null || f.obj == selected) selected = null;
					else selected = f.obj;
				}

				mousePt = e.getPoint();
				repaint();
			}

			@Override
			public void mouseDragged(MouseEvent e) {
				Point p = e.getPoint();
				double dx = p.x - mousePt.x;
				double dy = p.y - mousePt.y;

				if (e.isShiftDown()) {
					cam.shiftFocus(-dx / 100, dy / 100);
				}
				else {
					cam.rotateLocalY(dx / 100.0);
					cam.rotateLocalX(dy / 100.0);
				}
				cam.centerFocus();
				mousePt = p;
				repaint();
			}

			@Override
			public void mouseMoved(MouseEvent e) {
				if (isManipulating) {
					if (oldmouseloc == null) {
						oldmouseloc = e.getPoint();
					}
					Point loc = e.getPoint();
					double dx = loc.x - oldmouseloc.x;
					double dy = loc.y - oldmouseloc.y;
					switch (manipMode) {
					case 0:
						selected.move(cam
								.getImagePlaneVector(dx / 100, -dy / 100));
						break;
					}

					oldmouseloc = e.getPoint();
					repaint();
				}
			}

			@Override
			public void mouseWheelMoved(MouseWheelEvent e) {
				double zoomfactor = Math.pow(1.05, e.getWheelRotation());
				cam.zoom(zoomfactor);
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
					cam.setGlobalRotation(0, 0, 0);
				}
				else if (e.getKeyCode() == KeyEvent.VK_G) {
					if (selected != null && !isManipulating) {
						isManipulating = true;
						manipMode = 0;
						oldcenter = selected.center();
						oldmouseloc = null;
					}
				}
				else if (e.getKeyCode() == KeyEvent.VK_R) {
					rendering = !rendering;
					if (rendering) {
						status.setVisible(false);
						t.start();
						if (!renderthread.isAlive()) renderthread.start();
					}
					else {
						status.setVisible(true);
					}
				}
				else if (e.getKeyCode() == KeyEvent.VK_S) {
					if (selected != null) selected.smooth = !selected.smooth;
					else {
						for (Object3D mesh : meshes)
							mesh.smooth = !mesh.smooth;
					}
				}
				else if (e.getKeyCode() == KeyEvent.VK_Z) {
					if (e.isShiftDown()) {
						rendermode--;
						if (rendermode < 0) rendermode = 2;
					}
					else {
						rendermode++;
						if (rendermode > 2) rendermode = 0;
					}
				}
				else if (e.getKeyCode() == KeyEvent.VK_5) {
					cam.ortho = !cam.ortho;
				}

				cam.centerFocus();
				repaint();
			}
		};

		this.addKeyListener(k);
	}

	/**
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
	public Color calcPointLighting(Vertex v, Vertex n, Material m) {
		//flat shading, no distance falloff

		double r = 0, g = 0, b = 0, spr = 0, spg = 0, spb = 0;
		for (Lamp l : lamps) {
			Vertex lampvect = l.loc.subtract(v);
			Vertex lampvnorm = lampvect.getUnitVector(); //normalize lamp vector

			Vertex view = cam.viewVector(v).getUnitVector();
			//if lamp and view are on different sides of face, skip this lamp.
			if (n.dotproduct(lampvect) * n.dotproduct(view) > 0) continue;

			double dotprod = lampvnorm.dotproduct(n);
			double dstsqr = lampvect.lensquared();
 
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
			Vertex refl = lampvnorm.reflection(n);
			double spec_intensity = l.intensity
					* m.spintensity
					* Math.max(0, Math.pow(view.dotproduct(refl), m.sphardness));
			double diff_intensity = l.intensity * Math.abs(dotprod)
					* m.diff_intensity;

			//calc falloff
			switch (l.falloff) {
			case 2:
				diff_intensity *= 1.0 / dstsqr;
				spec_intensity *= 1.0 / dstsqr;
				break;
			case 1:
				diff_intensity *= 1.0 / Math.sqrt(dstsqr);
				spec_intensity *= 1.0 / Math.sqrt(dstsqr);
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
		Color col = colMultiply(m.diff_col, r, g, b); //diffuse color
		Color spcol = colMultiply(m.spcol, spr, spg, spb);//specular color
		return colAdd(col, spcol);
	}

	/**
	 * Calculates the color of a point by casting shadow rays, reflection rays,
	 * and transmission rays.
	 * 
	 * @param v
	 *            The point to be lit
	 * @param n
	 *            The normal of the point
	 * @param m
	 *            The material of the point
	 * @param faces
	 *            A list of the scene's faces
	 * @return The color of the point
	 */
	public Color traceRay(Vertex origin, Vertex ray, KDTree kdt, int depth) {

		raysTraced++;//For performance stats

		//tuv is used for current vertex, while tuvtmp is used for lamp intersections and whatnot
		Vertex tuv = new Vertex(0, 0, 0);
		Face f = rayIntersect(kdt, origin, ray, false, tuv);
		if (f == null) return world.getColor(ray);
		else {
			Material m = f.obj.mat;

			//calculate vertex location based on tuv
			Vertex v = origin.add(ray.scalarproduct(tuv.x));
			//calculate normal
			Vertex n;
			if (f.obj.smooth) {
				Vertex n0 = f.vertices[0].vertexNormal();
				Vertex n1 = f.vertices[1].vertexNormal();
				Vertex n2 = f.vertices[2].vertexNormal();
				n = n0.scalarproduct(1 - tuv.y - tuv.z)
						.add(n1.scalarproduct(tuv.y))
						.add(n2.scalarproduct(tuv.z));
				n.normalize();
			}
			else {
				n = f.normal;
			}

			double ndotray = n.dotproduct(ray);

			double r = 0, g = 0, b = 0, spr = 0, spg = 0, spb = 0;
			for (Lamp l : lamps) {
				Vertex lampvect = l.loc.subtract(v);

				//if lamp and view are on different sides of face, skip this lamp.
				if (n.dotproduct(lampvect) * ndotray > 0) continue;

				if (rayIntersect(kdt, v, lampvect, true, null) != null)
					continue;

				Vertex lampvnorm = lampvect.getUnitVector();
				double dotprod = lampvnorm.dotproduct(n);
				double dstsqr = lampvect.lensquared();
				if (dstsqr == 0) {
					//if lamp is on the face's center, add full brightness to each color that the lamp emits
					if (l.col.getRed() != 0) r += 1;
					if (l.col.getGreen() != 0) g += 1;
					if (l.col.getBlue() != 0) b += 1;
					continue;
				}

				//n is normalized, so
				//r = d- 2(d . n)n, where r = reflection, d = lampvector (normalized), n = normal
				//Phong reflection: spec intensity = (relf*view)^sphardness * intensity
				Vertex lamprefl = lampvnorm.reflection(n);
				double spec_intensity = l.intensity
						* m.spintensity
						* Math.max(0, Math.pow(ray.dotproduct(lamprefl), m.sphardness));
				double diff_intensity = l.intensity * Math.abs(dotprod)
						* m.diff_intensity;

				//calc falloff
				switch (l.falloff) {
				case 2:
					diff_intensity *= 1.0 / dstsqr;
					spec_intensity *= 1.0 / dstsqr;
					break;
				case 1:
					diff_intensity *= 1.0 / Math.sqrt(dstsqr);
					spec_intensity *= 1.0 / Math.sqrt(dstsqr);
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

			Color diffcol = colMultiply(m.col(f, tuv.y, tuv.z), r, g, b); //diffuse color
			Color spcol = colMultiply(m.spcol, spr, spg, spb);//specular color	
			Color totcol = colAdd(diffcol, spcol);

			//an approximation. Assumes normals point outside and doesn't really deal with with concentric/intersecting objects.Currently assumes 'outside' of every object is air.
			//also doesn't do fresnel formula, reflection and refraction are handled separately, except for total internal reflection
			if (depth < RAY_DEPTH) {
				//used both for reflection and refraction (since refraction sometimes reflects)
				Vertex refl = ray.reflection(n);

				if (m.alpha < 1) {
					double n1, n2;
					Vertex transray = Vertex.ORIGIN();

					if (ndotray < 0) {
						n1 = 1;
						n2 = m.ior;
					}
					else {
						n1 = m.ior;
						n2 = 1;
					}

					Vertex raynorm = n.scalarproduct(ndotray);
					Vertex raytang = ray.subtract(raynorm);
					Vertex transtang = raytang.scalarproduct(n1 / n2);
					double transsinsquared = transtang.lensquared();
					if (transsinsquared > 1) transray = refl;
					else {
						Vertex transnorm = n.scalarproduct(Math.signum(ndotray)
								* Math.sqrt(1 - transsinsquared));
						transray = transnorm.add(transtang);
					}

					Color transcol = traceRay(v, transray, kdt, depth + 1);
					totcol = colMix(transcol, totcol, m.alpha);
				}

				if (m.refl_intensity > 0) {
					Color refcol = traceRay(v, refl, kdt, depth + 1);
					totcol = colMix(totcol, refcol, m.refl_intensity);
				}
			}
			return totcol;
		}
	}

	public void paintComponent(Graphics gr) {
		super.paintComponent(gr);
		Graphics2D g2d = (Graphics2D) gr;
		g2d.setColor(Color.BLACK);

		if (rendering && render != null) {
			g2d.drawImage(render, null, 0, 0);
			if (KDTree.completion < 1) {
				status.setText("Building KD Tree... " + KDTree.completion * 100
						+ "%");
				status.setVisible(true);
			}
			else status.setVisible(false);
			return;
		}

		//create zbuffer for screen
		double[][] zbuffer = new double[getHeight()][getWidth()];
		for (double[] row : zbuffer)
			Arrays.fill(row, 1);

		//combine data of all meshes
		List<Edge> edges = new ArrayList<Edge>();
		HashMap<MeshVertex, Point> vertexPixels = new HashMap<MeshVertex, Point>(); //vertex pixel locations
		HashMap<MeshVertex, Double> vertexZvalues = new HashMap<MeshVertex, Double>(); //vertex Z values
		HashMap<MeshVertex, Color> vertexColors = new HashMap<MeshVertex, Color>(); //vertex colors
		for (Object3D mesh : meshes) {
			Collections.addAll(edges, mesh.edges);

			for (MeshVertex v : mesh.vertices) {
				vertexPixels.put(v, cam
						.projectVertex(v, getWidth(), getHeight()));
				if (rendermode > 1) {
					vertexZvalues.put(v, (cam.vertexDepth(v) - cam.mindist)
							/ (cam.maxdist - cam.mindist));
					if (mesh.smooth) {
						vertexColors.put(v, calcPointLighting(v, v
								.vertexNormal(), mesh.mat));
					}
				}
			}
		}
		Collections.sort(faces, cam.zsortFaces);

		Point v1, v2, v3;
		double z1, z2, z3;
		Color v1col, v2col, v3col;

//		g2d.setColor(Color.RED);
//		for (Edge e : kdedges) {
//			Point p1 = cam.projectVertex(e.v1, getWidth(), getHeight());
//			Point p2 = cam.projectVertex(e.v2, getWidth(), getHeight());
//			if (p1 != null && p2 != null) g2d.drawLine(p1.x, p1.y, p2.x, p2.y);
//		}

		//draw wireframe of selected
		if (selected != null) {
			g2d.setStroke(new BasicStroke(4));
			g2d.setColor(new Color(255, 100, 0));
			for (Edge e : selected.edges) {
				v1 = vertexPixels.get(e.v1);
				v2 = vertexPixels.get(e.v2);
				if (v1 == null || v2 == null) continue;

				if (rendermode == 2) {
					if (this.getBounds().contains(v1)
							&& vertexZvalues.get(e.v1) > zbuffer[v1.y][v1.x])
						continue;
					if (this.getBounds().contains(v2)
							&& vertexZvalues.get(e.v2) > zbuffer[v2.y][v2.x])
						continue;
				}
				g2d.drawLine(v1.x, v1.y, v2.x, v2.y);
			}
		}

		g2d.setStroke(new BasicStroke(1));
		if (rendermode == 0) {
			status.setText("Wireframe");
			g2d.setColor(Color.BLACK);
			for (Edge e : edges) {
				v1 = vertexPixels.get(e.v1);
				v2 = vertexPixels.get(e.v2);
				if (v1 == null || v2 == null) continue;
				g2d.drawLine(v1.x, v1.y, v2.x, v2.y);
			}
		}
		else if (rendermode == 1) {
			status.setText("Painter's Algorithm");
			Collections.reverse(faces);
			for (Face f : faces) {

				//optional backface culling
				//if(f.normal.dotproduct(cam.normal()) > 0.01) continue;

				v1 = vertexPixels.get(f.vertices[0]);
				v2 = vertexPixels.get(f.vertices[1]);
				v3 = vertexPixels.get(f.vertices[2]);

				if (v1 == null || v2 == null || v3 == null) continue;
				g2d.setColor(calcPointLighting(f.center(), f.normal, f.obj.mat));
				g2d.fillPolygon(new int[] { v1.x, v2.x, v3.x }, new int[] {
						v1.y, v2.y, v3.y }, 3);
			}
		}
		else if (rendermode == 2) {
			status.setText("Z-buffer");
			for (Face f : faces) {
				//get triangle corners
				v1 = vertexPixels.get(f.vertices[0]);
				v2 = vertexPixels.get(f.vertices[1]);
				v3 = vertexPixels.get(f.vertices[2]);

				if (v1 == null || v2 == null || v3 == null) continue;

				//zdepth values for vertices. For other pixels, lerp between these values 
				z1 = vertexZvalues.get(f.vertices[0]);
				z2 = vertexZvalues.get(f.vertices[1]);
				z3 = vertexZvalues.get(f.vertices[2]);

				//store difference values. Makes interpolation later on slightly faster 
				double dz21 = z2 - z1;
				double dz31 = z3 - z1;

				Vertex fcenter = f.center();//face center

				//vertex colors, for smooth shading (Gouraud shading, faster than Phong)
				v1col = vertexColors.get(f.vertices[0]);
				v2col = vertexColors.get(f.vertices[1]);
				v3col = vertexColors.get(f.vertices[2]);

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
				int w1_row = MathUtils.orient2D(v2, v3, p);
				int w2_row = MathUtils.orient2D(v3, v1, p);
				int w3_row = MathUtils.orient2D(v1, v2, p);

				int w = MathUtils.orient2D(v1, v2, v3);
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

								//only calculate face color once, unless smooth shading is enabled
								if (f.obj.smooth) col = lerpColor(v1col, v2col, v3col, (double) w1
										/ w, (double) w2 / w, (double) w3 / w);
								else {
									if (col == null)
										col = calcPointLighting(fcenter, f.normal, f.obj.mat);
								}
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

		if (selected != null) {
			status.setText(selected.name);
		}
		if (isManipulating) {
			String s = "";
			switch (manipMode) {
			case 0:
				s = "Moving";
				break;
			case 1:
				s = "Scaling";
				break;
			case 2:
				s = "Rotating";
				break;
			}
			status.setText(s + " " + selected.name + "...");
		}
	}

	public void raytraceRender() {
		//set up graphics
		render = new BufferedImage(getWidth(), getHeight(), BufferedImage.TYPE_3BYTE_BGR);
		Graphics g = render.getGraphics();
		g.setColor(Color.WHITE);
		g.fillRect(0, 0, getWidth(), getHeight());

		renderkdt = KDTree.buildTree(faces);
		long start = System.nanoTime();
		for (int y = 0; y < getHeight(); y += TILESIZE) {
			for (int x = 0; x < getWidth(); x += TILESIZE) {
				raytraceTile(x, Math.min(getWidth() - 1, x + TILESIZE - 1), y, Math
						.min(getHeight() - 1, y + TILESIZE - 1), g);
			}
		}
		long stop = System.nanoTime();
		long time = stop - start;
		System.out.println(faces.size() + " faces. Traced " + raysTraced
				+ " rays in " + time + "ns, or " + raysTraced
				/ (time / 1_000_000_000.0) + " rays per second.");

	}

	public void raytraceTile(int xmin, int xmax, int ymin, int ymax, Graphics g) {
		int i, j;
		Vertex[] ray;
		Color col;
		for (i = ymin; i <= ymax; i++) {
			for (j = xmin; j <= xmax; j++) {
				g.setColor(Color.BLACK);
				g.drawLine(j, i, j, i); //to show status
				ray = cam.castRay(j, i, getWidth(), getHeight());
				col = traceRay(ray[0], ray[1], renderkdt, 1);
				//col = (rayIntersect(renderkdt, ray[0], ray[1], true, null)!=null) ? Color.BLACK : Color.WHITE;
				g.setColor(col);
				g.drawLine(j, i, j, i);
			}
		}
	}

	public Face rayIntersect(KDTree kdt, Vertex origin, Vertex ray,
			boolean lazy, Vertex tuv) {
		//System.out.println("pos: "+kdt.pos + ", axis: "+kdt.axis);
		Vertex tuvtmp1 = Vertex.ORIGIN(), tuvtmp2 = Vertex.ORIGIN();
		if (tuv == null) tuv = Vertex.ORIGIN();
		Face f1, f2;

		//does ray intersect this bounding box?
		if (Face.rayAABBIntersect(kdt.bounds, origin, ray)) {
			//is this a leaf node?
			if (kdt.faces != null) {
				//System.out.println("intersect "+kdt.faces.size()+" faces.");
				f1 = rayIntersect(kdt.faces, origin, ray, lazy, tuvtmp1);
				tuv.setVertex(tuvtmp1);
				return f1;
			}
			else {
				f1 = rayIntersect(kdt.lower, origin, ray, lazy, tuvtmp1);
				if (lazy && f1 != null) {
					tuv.setVertex(tuvtmp1);
					return f1;
				}
				f2 = rayIntersect(kdt.upper, origin, ray, lazy, tuvtmp2);
				if (lazy && f2 != null) {
					tuv.setVertex(tuvtmp2);
					return f2;
				}

				//non-lazy compare returned face from each child
				if (f1 == null) {
					tuv.setVertex(tuvtmp2);
					return f2;
				}
				else if (f2 == null) {
					tuv.setVertex(tuvtmp1);
					return f1;
				}
				else {
					if (tuvtmp1.x < tuvtmp2.x) {
						tuv.setVertex(tuvtmp1);
						return f1;
					}
					else {
						tuv.setVertex(tuvtmp2);
						return f2;
					}
				}

			}
		}
		return null;
	}

	public Face rayIntersect(List<Face> faces, Vertex origin, Vertex ray,
			boolean lazy, Vertex tuv) {
		Face f = null;
		if (tuv == null) tuv = new Vertex(0, 0, 0);
		Vertex tuvtmp = new Vertex(0, 0, 0);
		double zmin = Double.NaN;

		for (Face ftmp : faces) {
			if (ftmp.intersectRayTriangle(origin, ray, tuvtmp)
					&& (tuvtmp.x < zmin || zmin != zmin)) {
				zmin = tuvtmp.x;
				f = ftmp;
				tuv.setVertex(tuvtmp);
				if (lazy) return f;
			}
		}

		return f;
	}

	/**
	 * Linearly interpolates colors based on the rgb color space, using
	 * (normalized) barycentric coordinates
	 * 
	 * @param c1
	 *            color
	 * @param c2
	 *            color
	 * @param c3
	 *            color
	 * @param w1
	 *            Barycentric coordinate
	 * @param w2
	 *            Barycentric coordinate
	 * @param w3
	 *            Barycentric coordinate
	 * @return
	 */
	public static Color lerpColor(Color c1, Color c2, Color c3, double w1,
			double w2, double w3) {
		double r = c1.getRed() * w1 + c2.getRed() * w2 + c3.getRed() * w3;
		double g = c1.getGreen() * w1 + c2.getGreen() * w2 + c3.getGreen() * w3;
		double b = c1.getBlue() * w1 + c2.getBlue() * w2 + c3.getBlue() * w3;

		Color col = new Color((int) MathUtils.clamp(r, 0, 255), (int) MathUtils.clamp(g, 0, 255), (int) MathUtils
				.clamp(b, 0, 255));
		return col;
	}

	/**
	 * Linearly interpolates two colors
	 * 
	 * @param a
	 * @param b
	 * @param v
	 *            between 0 and 1. 0 produces a, 1 produces b. 0.5 produces an
	 *            even mix, etc.
	 * @return
	 */
	public static Color colMix(Color a, Color b, double v) {
		int red = (int) (a.getRed() + v * (b.getRed() - a.getRed()));
		int green = (int) (a.getGreen() + v * (b.getGreen() - a.getGreen()));
		int blue = (int) (a.getBlue() + v * (b.getBlue() - a.getBlue()));

		return new Color(red, green, blue);
	}

	public static Color colBlend(Color... cols) {
		double red = 0, green = 0, blue = 0;
		for (Color c : cols) {
			red += c.getRed();
			green += c.getGreen();
			blue += c.getBlue();
		}
		red /= cols.length;
		green /= cols.length;
		blue /= cols.length;

		return new Color((int) red, (int) green, (int) blue);
	}

	public static Color colAdd(Color a, Color b) {
		return new Color(Math.min(255, a.getRed() + b.getRed()), Math.min(255, a
				.getGreen() + b.getGreen()), Math.min(255, a.getBlue()
				+ b.getBlue()));
	}

	public static Color colMultiply(Color c, double r, double g, double b) {
		double col_r = c.getRed() / 255.0;
		double col_g = c.getGreen() / 255.0;
		double col_b = c.getBlue() / 255.0;

		double out_r = MathUtils.clamp(r * col_r, 0, 1); //multiply, clamp
		double out_g = MathUtils.clamp(g * col_g, 0, 1); //multiply, clamp
		double out_b = MathUtils.clamp(b * col_b, 0, 1); //multiply, clamp

		return new Color((float) out_r, (float) out_g, (float) out_b);
	}

	class RenderThread extends Thread {
		public void run() {
			raytraceRender();
			this.interrupt();
		}
	}

	public static void main(String[] args) {
		JFrame f = new JFrame("Raph's renderer");
		f.setSize(800, 800);
		f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

		Renderer p = new Renderer();
		f.getContentPane().add(p);
		p.setFocusable(true);
		p.requestFocusInWindow();

		f.setVisible(true);
	}
}
