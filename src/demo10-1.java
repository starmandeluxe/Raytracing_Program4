//                      File demo10.java

// JFrame with Canvas. Sphere, ellipsoid, box, cylinder, cone,
// torus. 
// Shadows, no reflection, highlights. Arbitrary position, rotation 
// and scale of solids. Background is an object.
// One global normal Vector - no calls to "new" in normal().
// illum() sets global variables in MainThread and returns void.
// Uses thread for the main loop. IntersectAll() and illum() are 
// in the main loop thread.

import javax.swing.*;
import java.awt.*;
import java.awt.image.*;
import java.awt.event.*;
import java.lang.reflect.Method;


class Main
{  public static void main(String[] arg)
   {  new RFrame(800, 800);
   }
}

class RFrame extends JFrame
{  public RFrame(int w, int h)
   {  RCanvas canvas = new RCanvas();
      getContentPane().add("Center", canvas);
      setSize(w, h);
      setVisible(true);
      setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
   }
}

class RCanvas extends Canvas
{   Image image;
    Light lights;
    Solid all;
    private RCanvas canvas = this;

   public RCanvas() 
   {  lights = 
       new Light(0f,0f,0f, 0.4f,0.4f,0.4f,            // ambient
       new Light(4f,  6f,  3f,     0.1f,0.1f,0.1f, 
       new Light(4f,6.5f,  3.4f,   0.1f,0.1f,0.1f, 
       new Light(4.3f,6.1f,3.8f,   0.1f,0.1f,0.1f, 
       new Light(4.8f,6.2f,4.2f,   0.1f,0.1f,0.1f, 
       new Light(5f,  6.3f,4.6f,   0.1f,0.1f,0.1f, 
       new Light(3.3f,6.8f,5f,     0.1f,0.1f,0.1f, null)))))));    // point

      all = 
       new Box(
        new Point(0.1, 0.3, 0.1), 
        new Point(0.2,0.2,-0.4),
        new Point(0.7, -0.4, 0.1), 
        Red.class.getMethods()[0],

       //new Torus(
       // 0.4f, 0.03f, 0.03f,    // big radius, hor radius, ver radius
       // new Point(1.3,0.0,0),                      // rotation
       // new Point(-0.2, -0.3, 0.3),                // center
       // Red.class.getMethods()[0],
       
       new Ellipsoid(
        new Point(0.5, 0.1, 0.7),                // axes
        new Point(.6,-.4,0),                  // rotation
        new Point(-0.4, 0.6, 0.0),                // center 
        Orange.class.getMethods()[0],         // texture

       new Sphere(
        0.4,                                        // radius
        new Point(0.7,-0.15,0),                     // rotation
        new Point(0.2, -0.2, -0.1),                 // center
        Blue.class.getMethods()[0],
        
       new Sphere(.25, 
       new Point(0.05, 0.5, 0.0), 
       new Point(0.5, 0.6, 0.0), 
       Green.class.getMethods()[0],
       
       new Sphere(.1, 
       new Point(0.05, 0.5, 0.0), 
       new Point(0.45, 0.0, 0.1), 
       Blue.class.getMethods()[0],
       
       new Sphere(.1, 
       new Point(0.05, 0.5, 0.0), 
       new Point(0.1, 0.15, 0.1), 
       Blue.class.getMethods()[0],
       
       new Sphere(.1, 
       new Point(0.05, 0.5, 0.0), 
       new Point(-.08, -0.35, 0.18), 
       Blue.class.getMethods()[0],

       new Ellipsoid(
        new Point(0.4, 0.2, 0.2),                // axes
        new Point(0,-0.3,-0.2),                  // rotation
        new Point(0.5, 0.6, 0.0),                // center 
        Green.class.getMethods()[0],         // texture
 
       new Cylinder(
        new Point(0.2, 0.4, 0.3),            // axes
        new Point(0,0,.3),                  // rotation
        new Point(-0.4, -.35, 0.4),            // center
        LightBlue.class.getMethods()[0],         // texture

       //new Cone(
        //new Point(0.15, 0.3, 0.15),            // axes
        //new Point( 0.2,0,0.3),                 // rotation
        //new Point(0.25, 0.01, 1.0),            // center
        //OffWhite.class.getMethods()[0],        // texture
        
        new Box(
        new Point(.8,.09,2.3), 
        new Point(0,0,0),
        new Point(0.0, -.6, -1), 
        Purple.class.getMethods()[0],
            
       new Background(
        BackgroundColor.class.getMethods()[0],
        null
       )))))))))));

      addComponentListener(new ComponentAdapter()
      {  MainLoopThread mainLoopThread;
         public void componentResized(ComponentEvent ce)
         {  mainLoopThread = new MainLoopThread(canvas);
            mainLoopThread.start();
         }
      });
   } 

   public void paint(Graphics g)
   {  g.drawImage(image, 0, 0, this); 
   }
}


class MainLoopThread extends Thread
{  private RCanvas canvas;
   private Light lights;
   private Solid all;
   private Graphics gr;
   private Vector normal = new Vector(0,0,0);
   private RGB Refrgb = new RGB();;
   private RGB rgb = new RGB();
   //private double illumr;
   //private double illumg;
   //private double illumb;
   private SecInf u = new SecInf(0, null, new SecInf(0,null,
      new SecInf(0, null, new SecInf(0, null,null))));

   public MainLoopThread(RCanvas canvas)
   {  this.canvas = canvas;
      lights = canvas.lights;
      all = canvas.all;
   }

   public void run()
   {  mainLoop(canvas.getWidth(), canvas.getHeight());
   }

   private void mainLoop(int wid, int hei)
   {  gr = canvas.getGraphics();
      int[] pixels = new int[wid*hei];

      int swath = 30;
      int index = 0;
      //double[] u = new double[4];
      short r, g, b;
      MemoryImageSource source = null;
      
      //Sec Inf object
      //SecInf u = new SecInf(0, null, new SecInf(0,null,
      //new SecInf(0, null, new SecInf(0, null,null))));

      for (int y = 0; y < hei; y++)                     // main loop
      {  for (int x = 0; x < wid; x++)
         {  double xw = (2*x - wid)/((double)wid);
            double yw = (hei - 2*y)/((double)wid);
            Ray ray = new Ray(0, 0, 4, xw, yw, -4);  // original ray
            ray.normalize();
            
            Solid solid = intersectAll(all, ray, u);
            u.solid = solid;
            u.next.solid = solid;
            u.next.next.solid = solid;
            u.next.next.next.solid = solid;
            Point hit = new Point(
             ray.sx + u.u*ray.dx,
             ray.sy + u.u*ray.dy,
             ray.sz + u.u*ray.dz);
            illum(solid, all, hit, ray.getDir(), lights);
            r = (short)(rgb.r*255);
            g = (short)(rgb.g*255);
            b = (short)(rgb.b*255);
            pixels[index++] = 255<<24 | r<<16 | g<<8 | b;
         }
         if ((y+1)%swath == 0)
         {  source = new MemoryImageSource(wid,swath,pixels,(y-swath+1)*wid,wid); 
            canvas.image = canvas.createImage(source);
            gr.drawImage(canvas.image,0,y-swath+1,wid,swath,null);
         }
      }
      int rest = hei%swath;
      source = new MemoryImageSource(wid,rest,pixels,(hei-rest)*wid,wid);
      canvas.image = canvas.createImage(source);
      gr.drawImage(canvas.image,0,hei-rest,wid,rest,null);
                                       // prepare full image for paint()
      source = new MemoryImageSource(wid,hei,pixels,0,wid);
      canvas.image = canvas.createImage(source);   
   }

   private Solid intersectAll(Solid all, Ray ray, SecInf u)
   {  double closest = 1000000;
      Solid h = null;
      
      for (Solid s = all; s != null; s = s.next)
      {  if (s.intersect(ray, u) == true)
            if (u.u < closest)
            {  closest = u.u;
               h = s;
            }
      }
      if (h != null) u.u = closest;
      return h;
   }
   
   public float VectorDot(Vector a, Vector b)
    {
        return (float)((a.x * b.x) + (a.y * b.y) + (a.z * b.z));
    }
   
   public Ray CalculateReflectedRay(Vector dir, Point hit, Vector normal)
    {
	    // R = I - 2(N.I)*N
	    
	    //r = 2 * (l dot n ) * n - l
	    Ray reflectedRay;
	    //ray.normalize();
	    float nDotI = 2 * VectorDot(normal, dir); 
	    double x = dir.x - nDotI * normal.x;					
	    double y = dir.y - nDotI * normal.y; 
	    double z = dir.z - nDotI * normal.z;
	    
	    reflectedRay = new Ray(hit.x, hit.y, hit.z, x,y,z);
	    
	    return reflectedRay;
    }

   private RGB illum(Solid solid, Solid all, Point hp, Vector hv, Light lights)
   {  
   	  Refrgb = new RGB(); 	
   	  rgb = new RGB();
   	  solid.normal(hp, normal);
      Light amb = lights;
      double totptlr = 0;
      double totptlg = 0;
      double totptlb = 0;
      double tothighr = 0;
      double tothighg = 0;
      double tothighb = 0;

      Texture col = null;
      try                                // get texture at hit point
      {  col = (Texture)solid.texture.invoke(null, solid, hp);
      }
      catch(Exception e) 
      {  System.out.println(e);
      }

      Vector flecVec;
      if (col.spec > 0)
         flecVec = normal.reflect(hv);
      else 
         flecVec = null;
         
      Vector normal = new Vector(0,0,0); 
      solid.normal(hp, normal);
         Ray refl = CalculateReflectedRay(hv,hp,normal);
         if (col.ref > 0)
         {
	         refl.normalize(); //need this?
             Ray ray2light = new Ray(hp.x+.00001*normal.x, hp.y+.00001*normal.y, hp.z+.00001*normal.z, 
    			refl.getDir().x, refl.getDir().y, refl.getDir().z);
    		 ray2light.normalize();
    		 Solid hits = intersectAll(all,ray2light,u);
    		 if (hits != null)
    		 {
    		 	Point hitp = new Point(
	             ray2light.sx + u.u*ray2light.dx,
	             ray2light.sy + u.u*ray2light.dy,
	             ray2light.sz + u.u*ray2light.dz);
	             Refrgb = illum(hits, all, hitp, ray2light.getDir(), lights);
    		 }
    		 
	     }
	     else
	     	Refrgb = new RGB(0,0,0);
	     	
                                                        // light sources loop
      for (Light ptl = lights.next; ptl != null; ptl = ptl.next)
      {  double vtlx = ptl.x-hp.x;                      // vector to light
         double vtly = ptl.y-hp.y;
         double vtlz = ptl.z-hp.z;
         float distToLight = (float)Math.sqrt(vtlx*vtlx+vtly*vtly+vtlz*vtlz);
         vtlx /= distToLight;                    // normalize vector to light
         vtly /= distToLight;
         vtlz /= distToLight;
         
         
         
      
         double dot = vtlx*normal.x + vtly*normal.y + vtlz*normal.z;
         if (dot < 0) dot = 0;

         if (dot > 0)                                 // test for shadows
         {  
         	Ray rayToLight = new Ray(
             hp.x+normal.x*0.00001,
             hp.y+normal.y*0.00001,
             hp.z+normal.z*0.00001,
             vtlx, vtly, vtlz); 
            Solid test = intersectAll(all, rayToLight, u); // shadow test
            if (test != null && u.u < distToLight)
               dot = 0;
         }
         totptlr += dot*ptl.r;
         totptlg += dot*ptl.g;
         totptlb += dot*ptl.b;

         double high = 0;
         if (col.spec > 0 && dot > 0)                 // soft highlights
         {  high = vtlx*flecVec.x
                 + vtly*flecVec.y
                 + vtlz*flecVec.z;
            high = Math.pow(high, (int)(col.spec*500));
         }
         tothighr += high*ptl.r;
         tothighg += high*ptl.g;
         tothighb += high*ptl.b;
      }
      
      

      //rgb.r += ((amb.r+totptlr)*col.r+tothighr);
      //rgb.g += ((amb.g+totptlg)*col.g+tothighg);
      //rgb.b += ((amb.b+totptlb)*col.b+tothighb);
      
      rgb.r = ((1.0 - col.ref) * ((amb.r+totptlr)*col.r+tothighr)) 
                       + (col.ref * Refrgb.r); 
      rgb.g = ((1.0 - col.ref) * ((amb.g+totptlg)*col.g+tothighg)) 
                         + (col.ref * Refrgb.g); 
      rgb.b = ((1.0 - col.ref) * ((amb.b+totptlb)*col.b+tothighb)) 
                       + (col.ref * Refrgb.b); 

      if (rgb.r > 1) rgb.r = 1;                 // trimming
      if (rgb.g > 1) rgb.g = 1;
      if (rgb.b > 1) rgb.b = 1;
      
      return rgb;
   }
}                                 // end of class MainLoopThread


class Ray
{  public double sx, sy, sz, dx, dy, dz;
   public Ray(double sx, double sy, double sz,
              double dx, double dy, double dz)
   { this.sx = sx;
     this.sy = sy;
     this.sz = sz;
     this.dx = dx;
     this.dy = dy;
     this.dz = dz;
   }
   public Ray(Ray r)                         // copy constructor
   {  sx = r.sx;
      sy = r.sy;
      sz = r.sz;
      dx = r.dx;
      dy = r.dy;
      dz = r.dz;
   }
   public Ray(double sx,double sy, double sz, Vector dir)
   {  this(sx, sy, sz, dir.x, dir.y, dir.z);
   }
   public void normalize()
   {  double l = Math.sqrt(dx*dx+dy*dy+dz*dz);
      if (l != 0)
      {  dx /= l;
         dy /= l;
         dz /= l;
      }
   }
   public Vector getDir()
   {  return new Vector(dx, dy, dz);
   }
}

class Point
{  public double x, y, z;
   public Point(double x, double y, double z)
   {  this.x = x;
      this.y = y;
      this.z = z;
   }
}

class Vector
{  public double x, y, z;
   public Vector(double x, double y, double z)
   {  this.x = x;
      this.y = y;
      this.z = z;
   }
   public float normalize()
   {  float l = (float)Math.sqrt(x*x+y*y+z*z);
      if (l != 0)
      {  x /= l;
         y /= l;
         z /= l;
      }
      return l;
   }
   public Vector reflect(Vector in)      // reflect "in" about "this"
   {  double dot = 2*(x*in.x + y*in.y + z*in.z);
      return new Vector(                          // is normalized
       in.x - x*dot,
       in.y - y*dot,
       in.z - z*dot);
   }
}

class Texture 
{  public float r, g, b, spec, ref; 
   public Texture() {}
   public Texture(float r, float g, float b, float spec, float ref)
   {  this.r = r;
      this.g = g;
      this.b = b;
      this.spec = spec;
      this.ref = ref;
   }
}

class Light
{  public float x, y, z;
   public float r, g, b;
   public Light next;

   public Light(float x, float y, float z,
                float r, float g, float b, Light next)
   {  this.x = x;
      this.y = y;
      this.z = z;
      this.r = r;
      this.g = g;
      this.b = b;
      this.next = next;
   }
}

class BackgroundColor
{  static Texture backBlue = new Texture(1.5f,1.5f,2.1f,0,0);
   static public Texture texture(Solid s, Point hit)
   {  return backBlue;
   }
}
class LatLong                          // latitude - logitude lattice
{  static Texture black = new Texture(0f, 0f, 0f, 0,0);
   static Texture pink = new Texture(1f, 0.85f, 0.85f, 0f,.5f);
   static public Texture texture(Solid solid, Point h)
   {  Sphere s = (Sphere)solid;
      double[] m = s.u2a;
      double x, y, z;      // (x,y,z) is set to hitpoint on unit sphere

      x = ((h.x-s.cx)*m[0]+(h.y-s.cy)*m[3]+(h.z-s.cz)*m[6])/s.radius;
      y = ((h.x-s.cx)*m[1]+(h.y-s.cy)*m[4]+(h.z-s.cz)*m[7])/s.radius;
      z = ((h.x-s.cx)*m[2]+(h.y-s.cy)*m[5]+(h.z-s.cz)*m[8])/s.radius;

      double w;                    // compute longitude - pole to pole
      if (Math.abs(x) < 0.00001) 
         w = (z > 0) ? -Math.PI/2 : Math.PI/2;
      else
         w = Math.atan(-z/x);  // right: front -PI/2, back +PI/2
      w += 3*Math.PI/2;           // right: front PI,  back 2*PI
      if (x < -0.00001)
         w -= Math.PI;            // left: back 0, front PI
      w /= 2*Math.PI;             // 0 < w < 1: back ccw to back
      w *= 15;                         // scaling, number of stripes
      double v;                        // compute latitude
      v = (Math.asin(y)+Math.PI/2)/Math.PI; // 0 < v < 1, top = 1
      v *= 15;                         // scaling, number of stripes

      if (Math.abs((int)(w+0.5)-w) < 0.05
          ||
          Math.abs((int)(v+0.5)-v) < 0.05)
         return black;
      else
         return pink;
   }
}
class Red
{  static Texture red = new Texture(1f, 0f, 0f, 0.1f,.3f);
   static public Texture tex(Solid s, Point hit)
   {  return red;
   }
}
class Yellow
{  static Texture yellow = new Texture(1f, 1f, 0.4f, 0.05f,.3f);
   static public Texture tex(Solid s, Point hit)
   {  return yellow;
   }
}
class Green
{  static Texture green = new Texture(0f, 1f, 0f, 0.05f,.3f);
   static public Texture tex(Solid s, Point hit)
   {  return green;
   }
}
class Blue
{  static Texture blue = new Texture(0f, 0f, 1f, 0.05f,.4f);
   static public Texture tex(Solid s, Point hit)
   {  return blue;
   }
}
class Purple
{  static Texture purple = new Texture(1f, 0f, 1f, 0.05f,.7f);
   static public Texture tex(Solid s, Point hit)
   {  return purple;
   }
}
class Orange
{  static Texture orange = new Texture(1f, .8f, .2f, 0.05f,.2f);
   static public Texture tex(Solid s, Point hit)
   {  return orange;
   }
}

class LightBlue
{  static Texture lightBlue = new Texture(0.7f, 0.9f, 1f, 0.2f,.3f);
   static public Texture tex(Solid s, Point hit)
   {  return lightBlue;
   }
}
class OffWhite
{  static Texture offWhite = new Texture(1f, 1f, 0.92f, 0,.3f);
   static public Texture tex(Solid s, Point hit)
   {  return offWhite;
   }
}

//************  Base class Solid  *************

abstract class Solid
{  public double cx, cy, cz;      // center
   public Method texture;         // texture
   public double[] u2a = new double[9]; // rotation matrix
   public Solid next;             // link to next solid
   public Ray r = new Ray(0,0,0,0,0,0);

   public void transRot(Ray ray) //untranslate and unrotate - no unscaling
   {  r.sx = (ray.sx-cx)*u2a[0]+(ray.sy-cy)*u2a[3]+(ray.sz-cz)*u2a[6];
      r.sy = (ray.sx-cx)*u2a[1]+(ray.sy-cy)*u2a[4]+(ray.sz-cz)*u2a[7];
      r.sz = (ray.sx-cx)*u2a[2]+(ray.sy-cy)*u2a[5]+(ray.sz-cz)*u2a[8];
      r.dx = ray.dx*u2a[0]+ray.dy*u2a[3]+ray.dz*u2a[6];
      r.dy = ray.dx*u2a[1]+ray.dy*u2a[4]+ray.dz*u2a[7];
      r.dz = ray.dx*u2a[2]+ray.dy*u2a[5]+ray.dz*u2a[8];
   }

   public Solid(Point rot, Point center, Method texture, Solid next)
   {  cx = center.x;
      cy = center.y;
      cz = center.z;
      double cx = Math.cos(rot.x);         // original rotations
      double sx = Math.sin(rot.x);         // of the arbitrarily  
      double cy = Math.cos(rot.y);         // positioned solid
      double sy = Math.sin(rot.y);         // in the order x-y-z
      double cz = Math.cos(rot.z);
      double sz = Math.sin(rot.z);

      u2a[0] =  cy*cz;          // unit-to-arbitray without scale
      u2a[1] = -cx*sz - sx*sy*cz;   // order x-y-z:   Rz*Ry*Rx
      u2a[2] =  sx*sz - cx*sy*cz;   // matrix left - vector right
      u2a[3] =  cy*sz; 
      u2a[4] =  cx*cz - sx*sy*sz;
      u2a[5] = -sx*cz - cx*sy*sz;
      u2a[6] =  sy; 
      u2a[7] =  sx*cy; 
      u2a[8] =  cx*cy;

      this.texture = texture;
      this.next = next;
   }
   abstract public boolean intersect(Ray ray, SecInf u);
   abstract public void normal(Point hit, Vector v);
}

class Sphere extends Solid
{  public double radius;

   public Sphere(double radius, Point rot, Point center, 
    Method texture, Solid next)
   {  super(rot, center, texture, next);
      this.radius = radius;
   }
   
   public boolean intersect(Ray ray, SecInf u)      
   {  r.sx=ray.sx-cx; r.sy=ray.sy-cy; r.sz=ray.sz-cz;    // no rotation
      r.dx=ray.dx;    r.dy=ray.dy;    r.dz=ray.dz;       // no scaling

      double a1 = r.sx*r.dx + r.sy*r.dy + r.sz*r.dz;
      double a0 = r.sx*r.sx + r.sy*r.sy + r.sz*r.sz - radius*radius;
      double disc = a1*a1 - a0;
      if (disc < 0.00001)          // no solution to quadratic equation 
         return false;

      disc = Math.sqrt(disc);
      u.u = -a1 - disc;
      u.next.u = -a1 + disc;
      if (u.next.u <= 0)                           // both u-values < 0
         return false;
      
      return true;
   }

   public void normal(Point hit, Vector v)
   {  v.x = (hit.x-cx)/radius; 
      v.y = (hit.y-cy)/radius; 
      v.z = (hit.z-cz)/radius;
   }
}


class Ellipsoid extends Solid
{  private double ax, ay, az;

   public Ellipsoid(Point axes, Point rot, Point center, 
    Method texture, Solid next)
   {  super(rot, center, texture, next);
      ax = axes.x;
      ay = axes.y;
      az = axes.z;
   }
   
   public boolean intersect(Ray ray, SecInf u)
   {  transRot(ray);                      // untranslate and unrotate
      r.sx /= ax; r.sy /= ay; r.sz /= az;   // r is ray member 
      r.dx /= ax; r.dy /= ay; r.dz /= az;   // in base class

      double a2 = r.dx*r.dx + r.dy*r.dy + r.dz*r.dz;      // need full 
      double a1 = r.sx*r.dx + r.sy*r.dy + r.sz*r.dz;      // quadratic
      double a0 = r.sx*r.sx + r.sy*r.sy + r.sz*r.sz - 1;  // equation
      double disc = a1*a1 - a2*a0;
      if (disc < 0.00001) return false; // no solution to quadratic equation 
         
      disc = Math.sqrt(disc);
      u.u = (-a1 - disc)/a2;        // u[0] is always smaller
      u.next.u = (-a1 + disc)/a2;        // u[1] is always bigger
      if (u.next.u <= 0) return false;                 // both u-values < 0
         
      return true;
   }

   public void normal(Point hit, Vector v)
   {  double 
      x = ((hit.x-cx)*u2a[0]+(hit.y-cy)*u2a[3]+(hit.z-cz)*u2a[6])/(ax*ax),
      y = ((hit.x-cx)*u2a[1]+(hit.y-cy)*u2a[4]+(hit.z-cz)*u2a[7])/(ay*ay),
      z = ((hit.x-cx)*u2a[2]+(hit.y-cy)*u2a[5]+(hit.z-cz)*u2a[8])/(az*az);

      v.x = (x*u2a[0]+y*u2a[1]+z*u2a[2]);              // rotate back out
      v.y = (x*u2a[3]+y*u2a[4]+z*u2a[5]);
      v.z = (x*u2a[6]+y*u2a[7]+z*u2a[8]);

      v.normalize();
   }
} 


class Box extends Solid
{  private double ax, ay, az;

   public Box(Point axes, Point rot, Point center, 
    Method texture, Solid next)
   {  super(rot, center, texture, next);
      ax = axes.x;
      ay = axes.y;
      az = axes.z;
   }

   public boolean intersect(Ray ray, SecInf u)
   {  double v; 
      double uen = -1000000.0;
      double uex =  1000000.0;

      transRot(ray);                 // this computes ray member r
      r.sx /= ax; r.sy /= ay; r.sz /= az;       // r is ray member 
      r.dx /= ax; r.dy /= ay; r.dz /= az;       // in base class

      if (Math.abs(r.dx) < 0.00001)                 // parallel
      {  if (0.99999 < -r.sx)
            return false; 
         if (0.99999 < r.sx)
            return false;      
      }
      else 
      {  v = (-r.sx - 1)/r.dx;
         if (0 < r.dx) 
         {  if (uen < v) uen = v;
         }
         else  
            if (v < uex) uex = v;
         v = (-r.sx + 1)/r.dx;
         if (0 < -r.dx)  
         {  if (uen < v) uen = v;
         }
         else 
            if (v < uex) uex = v; 
      }

      if (Math.abs(r.dy) < 0.00001)                  // parallel
      {  if (0.99999 < -r.sy)
            return false; 
         if (0.99999 <  r.sy)
            return false;
      }
      else 
      {  v = (-r.sy - 1)/r.dy;
         if (0 < r.dy) 
         {  if (uen < v) uen = v; 
         }
         else 
            if (v < uex) uex = v; 
         v = (-r.sy + 1)/r.dy;
         if (0 < -r.dy) 
         {  if (uen < v) uen = v; 
         }
         else 
            if (v < uex) uex = v; 
      }

      if (Math.abs(r.dz) < 0.00001)                   // parallel
      {  if (0.99999 < -r.sz)
            return false;
         if (0.99999 < r.sz)
            return false; 
      }
      else 
      {  v = (-r.sz - 1)/r.dz;
         if (0 < r.dz) 
         {  if (uen < v) uen = v; 
         }
         else 
            if (v < uex) uex = v;
         v = (-r.sz + 1)/r.dz;
         if (0 < -r.dz) 
         {  if (uen < v) uen = v; 
         }
         else 
            if (v < uex) uex = v;
      }

      if (uen > uex || uex < 0) return false;

      u.u = uen;
      u.next.u = uex;
      return true;
   }

   public void normal(Point hit, Vector v)
   {  double 
      x = ((hit.x-cx)*u2a[0]+(hit.y-cy)*u2a[3]+(hit.z-cz)*u2a[6])/ax,
      y = ((hit.x-cx)*u2a[1]+(hit.y-cy)*u2a[4]+(hit.z-cz)*u2a[7])/ay,
      z = ((hit.x-cx)*u2a[2]+(hit.y-cy)*u2a[5]+(hit.z-cz)*u2a[8])/az;

      if      (x >  0.9999)             // test which side hit is on
      {  v.x = u2a[0]; v.y = u2a[3]; v.z = u2a[6]; 
      }
      else if (x < -0.9999)
      {  v.x = -u2a[0]; v.y = -u2a[3]; v.z = -u2a[6];
      }
      else if (y >  0.9999)
      {  v.x = u2a[1]; v.y = u2a[4]; v.z = u2a[7]; 
      }
      else if (y < -0.9999)
      {  v.x = -u2a[1]; v.y = -u2a[4]; v.z = -u2a[7]; 
      }
      else if (z >  0.9999)
      {  v.x = u2a[2]; v.y = u2a[5]; v.z = u2a[8]; 
      }      
      else if (z < -0.99999)
      {  v.x = -u2a[2]; v.y = -u2a[5]; v.z = -u2a[8]; 
      }
   }
}


class Cylinder extends Solid
{  private double ax, ay, az;

   public Cylinder(Point axes, Point rot, Point center, 
    Method texture, Solid next)
   {  super(rot, center, texture, next);
      ax = axes.x;
      ay = axes.y;
      az = axes.z;
   }

   public boolean intersect(Ray ray, SecInf u)
   {  transRot(ray);                      // this computes ray member r
      r.sx /= ax; r.sy /= ay; r.sz /= az; // unscale
      r.dx /= ax; r.dy /= ay; r.dz /= az;

      if (  r.sy >  0.99999 && r.dy > -0.000001
         || r.sy < -0.99999 && r.dy <  0.000001)
         return false;

      double a0 = r.sx*r.sx + r.sz*r.sz - 1,
             a1 = r.sx*r.dx + r.sz*r.dz;
      if (a0 > 0 && a1 > 0) return false;  // start outside, point away 

      double a2 = r.dx*r.dx + r.dz*r.dz;
      double dis = a1*a1 - a2*a0;
      if (dis < 0) return false;     // no intersect with unit cylinder 

      dis = Math.sqrt(dis);
      double uex = (-a1 + dis)/a2;
      if (uex < 0) return false;           // cylinder behind ray start

      double uen = (-a1 - dis)/a2;
      double pen, pex;                   // plane entry and exit values
      if (Math.abs(r.dy) < 0.000001)          // ray parallel to planes
      {  pen = -100000;
         pex =  100000;
      }
      else
      {  pen = r.dy < 0 ? ( 1 - r.sy)/r.dy : (-1 - r.sy)/r.dy;
         pex = r.dy < 0 ? (-1 - r.sy)/r.dy : ( 1 - r.sy)/r.dy;
      }

      uen = uen > pen ? uen : pen;              // largest entry point 
      uex = uex < pex ? uex : pex;              // smallest exit point 
      if (uen > uex || uex < 0) return false;

      u.u = uen;
      u.next.u = uex;
      return true;
   }

   public void normal(Point h, Vector v)
   {  double                 // (x y z) is set to hit on unit cylinder
      x = ((h.x-cx)*u2a[0]+(h.y-cy)*u2a[3]+(h.z-cz)*u2a[6])/ax,
      y = ((h.x-cx)*u2a[1]+(h.y-cy)*u2a[4]+(h.z-cz)*u2a[7])/ay,
      z = ((h.x-cx)*u2a[2]+(h.y-cy)*u2a[5]+(h.z-cz)*u2a[8])/az;

      if (0.99999 < y)                               // on top plane
      {  v.x = u2a[1];                               // u2a*(0 1 0)
         v.y = u2a[4]; 
         v.z = u2a[7];                    // no normalization needed
      }
      else if (y < -0.99999)                      // on bottom plane
      {  v.x = -u2a[1];                              // u2a*(0 -1 0)
         v.y = -u2a[4]; 
         v.z = -u2a[7];                   // no normalization needed
      }
      else                                                // on body
      {  x /= ax; z /= az;             // for elliptic cross section   
         v.x = u2a[0]*x + u2a[2]*z;                   // u2a*(x 0 z)
         v.y = u2a[3]*x + u2a[5]*z;
         v.z = u2a[6]*x + u2a[8]*z;
         v.normalize();
      }
   }
}

class Background extends Solid
{  public Background(Method texture, Solid next)
   {  super(new Point(0,0,0), new Point(0,0,0), texture, next);
   }
   public boolean intersect(Ray ray, SecInf u)
   {  u.u = 100000;                     // both very far away
      u.next.u = 100001;
      return true;
   }
   public void normal(Point h, Vector v)
   {  v.x = 0; v.y = 0; v.z = 0;
   }
}

class Cone extends Solid
{  private double ax, ay, az;

   public Cone(Point axes, Point rot, Point center, 
    Method texture, Solid next)
   {  super(rot, center, texture, next);
      ax = axes.x;
      ay = axes.y;
      az = axes.z;
   }
   
   public boolean intersect(Ray ray, SecInf u)
   {  transRot(ray);
      r.sx /= ax; r.sy /= ay; r.sz /= az;
      r.dx /= ax; r.dy /= ay; r.dz /= az;

      if (   r.sy > -0.00001 && r.dy > -0.000001   // outside planes  
          || r.sy < -0.99999 && r.dy <  0.000001)  // pointing away  
         return false;

      double a0 = r.sx*r.sx - r.sy*r.sy + r.sz*r.sz,
             a1 = r.sx*r.dx - r.sy*r.dy + r.sz*r.dz,
             a2 = r.dx*r.dx - r.dy*r.dy + r.dz*r.dz;
      if (a0 > 0 && a1 > 0 && a2 > 0) return false;
      if (a0 < 0 && r.sy > 0 && a2 > 0) return false;

      double dis = a1*a1 - a2*a0;
      if (dis < 0.000001) return false; // no intersect with unit cone 

      dis = Math.sqrt(dis);
      double u0 = (-a1 - dis)/a2,
             u1 = (-a1 + dis)/a2;
      if (u0 < 0 && u1 < 0)           // both intersections negative 
         return false;                     // cone behind startpoint 

      double p0, p1;
      if (Math.abs(r.dy) < 0.000001)
      {  p0 = -100000;
         p1 =  100000;
      }
      else
      {  p0 = r.dy < 0 ?       -r.sy/r.dy : (-1 - r.sy)/r.dy;
         p1 = r.dy < 0 ? (-1 - r.sy)/r.dy :       -r.sy/r.dy;
      }

      if (a2 < 0)                             // concave case, u1 < u0 
      {  if (r.dy > 0) u0 = p0;                            // steep up 
         else  u1 = p1;                                  // steep down 
      }
      else                                     // convex case, u0 < u1 
      {  if (u0 < p0) u0 = p0;                  // largest entry point 
         if (u1 > p1) u1 = p1;                  // smallest exit point 
      }

      if (u0 > u1 || u1 < 0) return false;

      u.u = u0;
      u.next.u = u1;
      return true;
   }

   public void normal(Point h, Vector v)
   {  double                    // (x y z) is set to hit on unit cone
      x = ((h.x-cx)*u2a[0]+(h.y-cy)*u2a[3]+(h.z-cz)*u2a[6])/ax,
      y = ((h.x-cx)*u2a[1]+(h.y-cy)*u2a[4]+(h.z-cz)*u2a[7])/ay,
      z = ((h.x-cx)*u2a[2]+(h.y-cy)*u2a[5]+(h.z-cz)*u2a[8])/az;

      if (y < -0.99999)            // bottom plane like in cylinder
      {  v.x = -u2a[1];                             // u2a*(0 -1 0)
         v.y = -u2a[4]; 
         v.z = -u2a[7];                  // no normalization needed
      }
      else                                     // hit point on body
      {  x /= ax; y /= ay; z /= az;   // for elliptic cross section
         v.x = u2a[0]*x-u2a[1]*y+u2a[2]*z; 
         v.y = u2a[3]*x-u2a[4]*y+u2a[5]*z;
         v.z = u2a[6]*x-u2a[7]*y+u2a[8]*z;
         v.normalize();
      }
   }
}

class Torus extends Solid
{  private double R;
   private double hr;           // normalized small horizontal radius
   private double rsqu;                            // square of above
   private double vr;                              // vertical radius
   private double[] x = new double[4];
   private double vh;                               // ratio of vr/hr

   public Torus(
    double R,                              // large radius
    double r,                              // horizontal small radius
    double v,                              // vertical small radius
    Point rot, Point center, Method texture, Solid next)
   {  super(rot, center, texture, next);
      this.R = R;
      hr =  r/R;                 // normalize horizontal small radius
      rsqu = hr*hr;
      vr = v;
      vh = vr/hr;
   }
          // returns real part of cube root of a complex number (a b)
   static double realcbrt(double a, double b)                            
   {  double c = 0.866025403784438;
      double r = a*a + b*b;

      if (r > 0) r = Math.exp(Math.log(r)/6);           // sixth root
      if (b == 0)                                 // argument is real 
         return a < 0 ? -r : r;

      if (Math.abs(a) > 1.0e-24)               // argument is complex
         c = Math.cos(Math.atan(b/a)/3+(a<0 ? 1.047197551196598 : 0));
      return r*c;
   }
             // real solution of cubic equation x^3 + p*x^2 + q*x + r
   static double cubic(double p, double q, double r)
   {  double s = p/3,              
             t = s*s,
             a = q/3 - t,
             b = s*(q/2 - t) - r/2,
             y = b*b + a*a*a,
             rad = Math.sqrt(Math.abs(y));

      if (y < 0)                                      // complex case
         y = 2*realcbrt(b, rad) - s;
      else                                               // real case
         y = realcbrt(b+rad, 0.0) + realcbrt(b-rad, 0.0) - s;
      return y;
   }

   static void quartic(    // four real solutions of quartic equation
      double a, double b, double c, double d, double[] x)
   {
      double y = cubic(-b, a*c-4*d, -a*a*d+4*b*d-c*c);
      double res = 0.25*a*a - b + y;

      if (res < 0)
      {  x[0] = x[1] = x[2] = x[3] = -1;
         return;
      }

      double R  = Math.sqrt(res);
      double t1 = 0.75*a*a - R*R - 2*b;
      double t2 = (a*(b - 0.25*a*a) - 2*c)/R;

      if (t1 - t2 >= 0)
      {  double root = Math.sqrt(t1 - t2);
         x[0] = -a/2 - R - root;
         x[1] = -a/2 - R + root;
      }
      else
      {  x[0] = x[1] = -1;
      }

      if (t1 + t2 >= 0)
      {  double root = Math.sqrt(t1 + t2);
         x[2] = -a/2 + R - root;
         x[3] = -a/2 + R + root;
      }
      else
      {  x[2] = x[3] = -1;
      }
   }

   public boolean intersect(Ray ray, SecInf u)
   {  double t2, t1, t0;   // unit torus: R = 1, 0 < hr = vr < 1
      double a4, a3, a2, a1, a0;

      transRot(ray);
                                   //(vr/(hr*R))  is ratio of v to r
      r.sx /= R; r.sy /= vh; r.sz /= R;              //R*(vr/(hr*R))
      r.dx /= R; r.dy /= vh; r.dz /= R;             // R cancels out

      if (r.sy >  hr && r.dy > 0
      ||  r.sy < -hr && r.dy < 0)
         return false;

      t2 = r.dx*r.dx+r.dy*r.dy+r.dz*r.dz; 
      t1 = r.dx*r.sx+r.dy*r.sy+r.dz*r.sz; 
      t0 = r.sx*r.sx+r.sy*r.sy+r.sz*r.sz - 1 - rsqu; 
      a4 = 1/(t2*t2);
      a3 = 2*t1/t2;
      a2 = (t1*t1 + 0.5*t2*t0 + r.dy*r.dy)*a4;
      a1 = (0.5*t1*t0 + r.sy*r.dy)*a4;
      a0 = (0.0625*t0*t0 - 0.25*(rsqu - r.sy*r.sy))*a4;

      quartic(a3, a2, a1, a0, x);

      if (x[0] == Double.NEGATIVE_INFINITY
       || x[2] == Double.NEGATIVE_INFINITY)
         return false;

      boolean first = 0.0001 < Math.abs(x[2]-x[3]) && 0 < x[3];
      boolean second = 0.0001 < Math.abs(x[0]-x[1]) && 0 < x[1];

      if (first && second)
      {  u.u = x[0];
         u.next.u = x[1]; 
         u.next.next.u = x[2];
         u.next.next.next.u = x[3]; 
         return true;
      }
      if (first && !second)
      {  u.u = x[2];
         u.next.u = x[3]; 
         return true;
      }
      if (!first && second)
      {  u.u = x[0];
         u.next.u = x[1]; 
         return true;
      }

      return false;
   }


   public void normal(Point h, Vector n)
   {  double[] m = u2a;
      double x = (m[0]*(h.x-cx)+m[3]*(h.y-cy)+m[6]*(h.z-cz))/R;
      double y = (m[1]*(h.x-cx)+m[4]*(h.y-cy)+m[7]*(h.z-cz))/vh;
      double z = (m[2]*(h.x-cx)+m[5]*(h.y-cy)+m[8]*(h.z-cz))/R;

      double px = x, pz = z, len = Math.sqrt(x*x + z*z);
      px /= len; pz /= len;

      n.x = (m[0]*(x-px)+m[1]*y+m[2]*(z-pz))/R; 
      n.y = (m[3]*(x-px)+m[4]*y+m[5]*(z-pz))/vh;
      n.z = (m[6]*(x-px)+m[7]*y+m[8]*(z-pz))/R;
      n.normalize();
   }   
}
