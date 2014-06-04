import java.awt.Point;

public class MathUtils {
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
	public static int orient2D(Point a, Point b, Point c) {
		return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
	}

	public static double lerp(double a, double b, double t) {
		return a + (b - a) * t;
	}

	/**
	 * Smoothly interpolates between two points, so that 1st derivative of curve
	 * is 0 at endpoints
	 * 
	 * @param a
	 * @param b
	 * @param t
	 * @return <code>lerp(a, b, t * t * (3 - 2 * t))</code>
	 */
	public static double smoothstep(double a, double b, double t) {
		double r = t * t * (3 - 2 * t);
		return a + r * (b - a);

	}

	/**
	 * Smoothly interpolates between two points, so that 1st and 2nd derivatives
	 * are 0 at endpoints
	 * 
	 * @param a
	 * @param b
	 * @param t
	 * @return
	 */
	public static double smootherstep(double a, double b, double t) {
		double r = t * t * t * (t * (t * 6 - 15) + 10);
		return a + r * (b - a);
	}

	public static double clamp(double v, double min, double max) {
		if (v < min) return min;
		if (v > max) return max;
		return v;
	}
}
