import java.awt.Color;

public class World {
	public Color bgcol;
	
	public World(Color bgcol){
		this.bgcol = bgcol;
	}
	
	public Color getColor(Vertex dir){
		//TODO add horizons etc
		return bgcol;
	}
}
