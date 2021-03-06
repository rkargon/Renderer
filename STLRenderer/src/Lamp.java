import java.awt.Color;

public class Lamp {
	public double intensity;
	public Vertex loc;
	public Color col;
	
	/**
	 * Falloff modes:
	 * 0 - constant (no falloff)
	 * 1 - inv. linear 
	 * 2 - inv. square (default)
	 */
	public int falloff;

	public Lamp(){
		this(1, Vertex.ORIGIN(), Color.WHITE, 2);
	}
	
	public Lamp(double i, Vertex p, Color c, int falloff) {
		this.intensity = i;
		this.loc = p;
		this.col = c;
		this.falloff=falloff;
	}

}
