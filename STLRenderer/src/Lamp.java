import java.awt.Color;

public class Lamp {
	public double intensity;
	public Vertex loc;
	public Color col;

	public Lamp(){
		this(1, Vertex.ORIGIN(), Color.WHITE);
	}
	
	public Lamp(double i, Vertex p, Color c) {
		this.intensity = i;
		this.loc = p;
		this.col = c;
	}
}
