import java.awt.Color;

public class World {
	public Color horizcol, zenithcol;
	boolean isFlat;

	public World(Color bgcol) {
		this.horizcol = bgcol;
		isFlat = true;
	}

	public World(Color bgcol, Color zenithcol) {
		this.horizcol = bgcol;
		this.zenithcol = zenithcol;
		isFlat = false;
	}

	public Color getColor(Vertex dir) {
		if (isFlat) {
			return horizcol;
		}
		else{
			double uplen = dir.dotproduct(new Vertex(0, 0, 1));
			uplen*=uplen;
			double r = uplen/dir.lensquared();
			return Renderer.colMix(horizcol, zenithcol, r);
		}
	}
}
