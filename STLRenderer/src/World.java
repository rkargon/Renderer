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
			return Renderer.colMix(horizcol, zenithcol, Math.abs(dir.z));
		}
	}
}
