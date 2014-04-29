import java.awt.Color;

public class Material {

	public Color diff_col = Color.WHITE;
	public double diff_intensity=1;
	
	public Color spcol = Color.WHITE;
	public double spintensity = 1;
	public double sphardness = 128;

	public Material(Color col) {
		this.diff_col = col;
	}


	public Material(Color col, double diff_intensity, Color spcol, double spintensity, double sphardness) {
		super();
		this.diff_col = col;
		this.spcol = spcol;
		this.spintensity = spintensity;
		this.sphardness = sphardness;
	}

}
