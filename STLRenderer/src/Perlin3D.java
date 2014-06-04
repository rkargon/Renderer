import java.awt.Color;

public class Perlin3D implements Texture {

	@Override
	public Color getCol(Face f, double w1, double w2) {
		float col = (float) getVal(f, w1, w2);
		col = (float) (col / 1.8 + 0.5);
		return new Color(col, col, col);
	}

	@Override
	public double getVal(Face f, double w1, double w2) {
		Vertex v = f.vertices[0].scalarproduct(1 - w1 - w2)
				.add(f.vertices[1].scalarproduct(w1))
				.add(f.vertices[2].scalarproduct(w2));

		double val=0, valtmp, scale=1;
		for(int i=1; i<=depth; i++){
			valtmp = PerlinNoise3D(v.x, v.y, v.z, size/scale);
			
			val+=valtmp/scale;
			scale*=2;
		}
		
		//used to scale output to [0, 1]
		double rangefactor = 2*(1-1/scale)*1.8;
		val = val/rangefactor+0.5;
		return val;
	}

	private double size, depth;
	private Vertex[] grads3d;
	private int[] perms;

	/**
	 * Initializes a Perlin noise generator by setting up gradient and
	 * permutation arrays.
	 */
	public Perlin3D(double size, int depth) {
		this.size = size;
		this.depth = Math.max(1, depth);
		initializePerlin3D();
	}

	public int[] initializePerlin3D() {

		// spherical coordinates, theta is plane angle, phi is "lattitude"
		double theta = 0, phi = Math.PI / (8 + 1); // phi is set up to avoid
													// poles; poles added
													// manually

		int k, tmp, loc;
		grads3d = new Vertex[64 + 2];
		perms = new int[256];

		// manually assign poles, so they show up only once (when phi=0, all
		// theta map to same point)
		grads3d[0] = new Vertex(0, 0, 1);
		grads3d[1] = new Vertex(0, 0, -1);

		// set up Perlin gradients
		for (int i = 0; i < 8; i++) {
			for (int j = 0; j < 8; j++) {
				loc = i * 8 + j + 2;
				grads3d[loc] = new Vertex(Math.sin(phi) * Math.cos(theta), Math.sin(phi)
						* Math.sin(theta), Math.cos(phi));

				theta += 2 * Math.PI / 8;
			}
			phi += Math.PI / (8 + 1);
		}

		// fill permutations array
		for (int i = 0; i < 256; i++)
			perms[i] = i;
		// shuffle permutations
		for (int i = 0; i < 256; i++) {
			k = (int) (Math.random() * 256);
			tmp = perms[i];
			perms[i] = perms[k];
			perms[k] = tmp;
		}
		return perms;
	}

	public double PerlinNoise3D(double x, double y, double z, double del) {
		double sx = x / del - Math.floor(x / del);
		double sy = y / del - Math.floor(y / del);
		double sz = z / del - Math.floor(z / del);

		// indices for nearest grid points
		// w = west, e = east, n = north, s = south, b= bottom, t = top
		int w = (int) Math.floor(x / del) & 255;
		int e = (w + 1) & 255;
		int n = (int) Math.floor(y / del) & 255;
		int s = (n + 1) & 255;
		int b = (int) Math.floor(z / del) & 255;
		int t = (b + 1) & 255;

		// gradient values at corners of cube around point
		Vertex grad_nwb = grads3d[perms[(w + perms[(n + perms[b]) % 256]) % 256]
				% grads3d.length];
		Vertex grad_neb = grads3d[perms[(e + perms[(n + perms[b]) % 256]) % 256]
				% grads3d.length];
		Vertex grad_swb = grads3d[perms[(w + perms[(s + perms[b]) % 256]) % 256]
				% grads3d.length];
		Vertex grad_seb = grads3d[perms[(e + perms[(s + perms[b]) % 256]) % 256]
				% grads3d.length];
		Vertex grad_nwt = grads3d[perms[(w + perms[(n + perms[t]) % 256]) % 256]
				% grads3d.length];
		Vertex grad_net = grads3d[perms[(e + perms[(n + perms[t]) % 256]) % 256]
				% grads3d.length];
		Vertex grad_swt = grads3d[perms[(w + perms[(s + perms[t]) % 256]) % 256]
				% grads3d.length];
		Vertex grad_set = grads3d[perms[(e + perms[(s + perms[t]) % 256]) % 256]
				% grads3d.length];

		// difference vectors, between point x,y and grid points
		Vertex diff_nwb = new Vertex(sx, sy, sz);
		Vertex diff_neb = new Vertex(sx - 1, sy, sz);
		Vertex diff_swb = new Vertex(sx, sy - 1, sz);
		Vertex diff_seb = new Vertex(sx - 1, sy - 1, sz);
		Vertex diff_nwt = new Vertex(sx, sy, sz - 1);
		Vertex diff_net = new Vertex(sx - 1, sy, sz - 1);
		Vertex diff_swt = new Vertex(sx, sy - 1, sz - 1);
		Vertex diff_set = new Vertex(sx - 1, sy - 1, sz - 1);

		// dot products of gradient and different values
		double dotp_nwb = grad_nwb.dotproduct(diff_nwb);
		double dotp_neb = grad_neb.dotproduct(diff_neb);
		double dotp_swb = grad_swb.dotproduct(diff_swb);
		double dotp_seb = grad_seb.dotproduct(diff_seb);
		double dotp_nwt = grad_nwt.dotproduct(diff_nwt);
		double dotp_net = grad_net.dotproduct(diff_net);
		double dotp_swt = grad_swt.dotproduct(diff_swt);
		double dotp_set = grad_set.dotproduct(diff_set);

		// interpolate along y axis to get four values
		double inter_y_wb = MathUtils.smootherstep(dotp_nwb, dotp_swb, sy);
		double inter_y_eb = MathUtils.smootherstep(dotp_neb, dotp_seb, sy);
		double inter_y_wt = MathUtils.smootherstep(dotp_nwt, dotp_swt, sy);
		double inter_y_et = MathUtils.smootherstep(dotp_net, dotp_set, sy);

		// interpolate along x axis to get two values
		double inter_x_b = MathUtils.smootherstep(inter_y_wb, inter_y_eb, sx);
		double inter_x_t = MathUtils.smootherstep(inter_y_wt, inter_y_et, sx);

		// interpolate along z axis to get one, final value
		double inter_z = MathUtils.smootherstep(inter_x_b, inter_x_t, sz);

		return inter_z;
	}

}
