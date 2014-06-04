import java.awt.Color;

public class Material {

	public Color diff_col = Color.WHITE;
	public double diff_intensity = 1;
	Texture col_tex;
	Texture bump_tex;

	public Color spcol = Color.WHITE;
	public double spintensity = 1;
	public double sphardness = 128;

	public double relf_intensity = 0;
	public double alpha = 0.3;
	public double ior = 1.7;

	public Material(Color col) {
		this.diff_col = col;
	}

	public Material(Color col, double diff_intensity, Color spcol, double spintensity, double sphardness, double refl_intensity, double trans_alpha, double ior) {
		super();
		this.diff_col = col;
		this.spcol = spcol;
		this.spintensity = spintensity;
		this.sphardness = sphardness;
		this.relf_intensity = refl_intensity;
		this.alpha = trans_alpha;
		this.ior = ior;
	}
	
	public Color col(Face f, double w1, double w2){
		if(col_tex==null) return diff_col;
		else return col_tex.getCol(f, w1, w2);
	}

}
