#include "glm-0.9.7.1/glm/glm.hpp"
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <deque>
#include <stack>

#define M_RAD 0.017453292519943295769236907684
//#define M_PI 3.141592653589793238462643383279
using namespace glm;
using namespace std;

static const vec3 zero = vec3(0.0);

class Ray{
    public:
    vec3 origin;
    vec3 direc;
    float t_min, t_max;
    Ray(vec3 o, vec3 d){
        origin = o;
        direc = d;
        t_min = -1e34;
        t_max = 1e34;
    }
};

class LocalGeo{
    public:
    vec3 position;
    vec3 normal;
    LocalGeo(){}
    LocalGeo(vec3 p, vec3 n){
        position = p;
        normal = n;
    }
};

class BRDF{
    public:
        vec3 kd;
        vec3 ks;
        vec3 ka;
        vec3 ke;
        float kr;
        float ksh;
        float kt;
        BRDF(){
            kd = vec3(1.0, 0.0, 0.0);
            ks = vec3(1.0, 1.0, 1.0);
            ka = vec3(0.0, 0.0, 0.0);
            kr = 0.3;
            kt = 0.0;
            ke = ka;
            ksh = 50;
        }
        BRDF(vec3 d, vec3 s, vec3 a, vec3 e, float r, float t, float sh){
            kd = d;
            ks = s;
            ka = a;
            kr = r;
            kt = t;
            ke = e;
            ksh = sh;
        }
        BRDF(vec3 d ){
            kd = d;
            ks = vec3(1.0, 1.0, 1.0);
            ka = vec3(0.0, 0.0, 0.0);
            kr = 0.3;
            kt = 0.0;
            ke = ka;
            ksh = 50;
        }
};
class Primitive;

class Intersection{
    public:
    LocalGeo localgeo;
    Primitive * prim;
    Intersection(LocalGeo loc, Primitive* p){
        localgeo = loc;
        prim = p;
    }
};
class Material;
class Primitive{
    public:
    Material * mat;
    virtual bool intersect(Ray& ray, float& t,Intersection ** in) = 0;
    virtual bool intersectP(Ray& ray) = 0;
    virtual BRDF getBRDF() = 0;
};

class Camera{
    public:
    vec3 eyepos;
    vec3 center;
    vec3 up;
    vec3 lookat;
    vec3 u;
    vec3 v;
    vec3 w;
    float fovy;
    float width, height;
    float invwidth, invheight, aspectratio;
    Camera(vec3 eye, vec3 cen, vec3 upvec, float fov, unsigned x, unsigned y){
	    eyepos = eye;
        center = cen;
        up = upvec;
        fovy = (M_PI/180)*fov;
        width = x;
        height = y;
        invwidth = 1/(float)width;
        invheight = 1/(float)height;
        aspectratio = width*invheight;
        lookat = (eyepos - center);
        //if(lookat[0] != 0.0f && lookat[1] != 0.0f && lookat[2] != 0.0f){
            w = normalize(lookat);
        //}else{
        //    w = vec3(0.0);
        //}
        //vec3 norm = cross(up, w);

        //if( norm[0] != 0.0f && norm[1] != 0.0f && norm[2] != 0.0f){
            u = normalize(cross(up, w));
        //}else{
        //    u = vec3(0.0);
        //}
        v = cross(w, u);
    }
    Camera(){}
    Ray generateRay(int x, int y){
        vec3 origin = eyepos;
        float a = tan(fovy*0.5)*aspectratio * ((x - width*0.5)/(width*0.5));
        float b = tan(fovy*0.5)*( (height*0.5 - y)/(height*0.5));
        vec3 norm = a*u + b*v - w;
        vec3 direction = vec3(0.0);
        //if( norm[0] != 0.0f && norm[1] != 0.0f && norm[2] != 0.0f){
            direction  = normalize(a*u + b*v - w); 
        //} 
        return Ray(origin, direction);
    }
};

class Shape{
    public:
    virtual bool intersect(Ray& ray,float& t_hit, LocalGeo* local) = 0;
    virtual bool intersectP(Ray& ray) = 0;
};

class Sphere:public Shape{
    public:
    vec3 center;
    float radius;
    float radius2;
    Sphere();
    Sphere(vec3 c, float r){
        center = c;
        radius = r;
        radius2 = radius*radius;
    }
    
    bool intersect(Ray& ray, float& t_hit, LocalGeo* local){
        vec3 co = center - ray.origin;
        float tca = dot(co, ray.direc);
        if( tca < 0 ) return false;
        float dist = dot(co, co) - tca*tca;
        if( dist > radius2) return false;
        float thc = sqrt(radius2 - dist);
        float t0 = tca - thc;
        float t1 = tca + thc;

        if(t0 > t1)
            swap(t0, t1);

        if(t0 < 0){
            t0 = t1;
            if( t0 < 0)
                return false;
        }
        t_hit = t0;
        vec3 p = ray.origin + t0*ray.direc; 
        local->position = p;
        vec3 norm = p - center;
        //if( norm[0] != 0.0f && norm[1] != 0.0f && norm[2] != 0.0f){
            local->normal = normalize( p - center );
        //}else{
        //    local->normal = vec3(0.0);
        //    cout << "wtf" << endl;
        //}
        return true;
    }

    bool intersectP(Ray& ray){
        vec3 co = center - ray.origin;
        float tca = dot(co, ray.direc);
        if( tca < 0 ) return false;
        float dist = dot(co, co) - tca*tca;
        if( dist > radius2) return false;
        float thc = sqrt(radius2 - dist);
        float t0 = tca - thc;
        float t1 = tca + thc;
        if(t0 > t1)
            swap(t0, t1);

        if(t0 < 0){
            t0 = t1;
            if( t0 < 0)
                return false;
        }
        return true;
    }
};

class Triangle:public Shape{
    public:
    vec3 v1;
    vec3 v2;
    vec3 v3;
    Triangle( vec3 p1, vec3 p2, vec3 p3){
        v1 = p1;
        v2 = p2;
        v3 = p3;
    }

    bool intersect(Ray& ray, float& t_hit, LocalGeo* local){
        vec3 v2v1 = v2 - v1;
        vec3 v3v1 = v3 - v1;
        vec3 pvec = cross(ray.direc, v3v1);
        float det = dot(v2v1, pvec);
        
        if(det < 1e-8) return false;

        float invdet = 1/det;

        vec3 tvec = ray.origin - v1;
        float u = dot(tvec, pvec)*invdet;
        if(u < 0 || u > 1)
            return false;

        vec3 qvec = cross(tvec, v2v1);
        float v = dot(ray.direc, qvec)*invdet;
        if( v < 0 || u + v > 1) 
            return false;

        t_hit = dot(v3v1, qvec)*invdet;
        
        local->position = ray.origin + t_hit*ray.direc;
        vec3 norm  = cross(v2v1, v3v1); 
        //if( norm[0] != 0.0f && norm[1] != 0.0f && norm[2] != 0.0f){
            local->normal = normalize(cross(v2v1, v3v1));
        //}else{
        //    local->normal = vec3(0.0);
        //}
        return true;
    }
    bool intersectP(Ray& ray){
        vec3 v2v1 = v2 - v1;
        vec3 v3v1 = v3 - v1;
        vec3 pvec = cross(ray.direc, v3v1);
        float det = dot(v2v1, pvec);
        
        if(det < 1e-8) return false;

        float invdet = 1/det;

        vec3 tvec = ray.origin - v1;
        float u = dot(tvec, pvec)*invdet;
        if(u < 0 || u > 1)
            return false;

        vec3 qvec = cross(tvec, v2v1);
        float v = dot(ray.direc, qvec)*invdet;
        if( v < 0 || u + v > 1) 
            return false;

        return true;
    }
};

class ShaderProgram{
    public:
    virtual vec3 shader(LocalGeo &, Ray&, vec3 &, BRDF &, vec3 & ) = 0;
};

class BlinnPhong:public ShaderProgram{
    public:
    BlinnPhong(){}
    vec3 shader(LocalGeo & loc, Ray& lray, vec3 & lcolor, BRDF & brdf, vec3 & eye){
       vec3 diffuse = lcolor*brdf.kd*glm::max(float(0), dot(loc.normal, lray.direc));
       vec3 half = lray.direc + eye;
       //if( half[0] != 0.0f && half[1] != 0.0f && half[2] != 0.0f)
            half = normalize(lray.direc - eye);
       float nDotH = dot(loc.normal, half);
       vec3 specular = lcolor*brdf.ks*pow(glm::max(nDotH, 0.0f), brdf.ksh);
       //cout << specular[0] << " "  << specular[1] << " " << specular[2] << endl;
       return diffuse + specular;
    }
};

class CellShader:public ShaderProgram{
    public:
    CellShader(){}
    vec3 shader(LocalGeo & loc, Ray& lray, vec3 & lcolor, BRDF & brdf, vec3 & eye){
        float nDotL = dot(loc.normal, lray.direc);
        if( nDotL > 0.7 && nDotL < 1.0){
            nDotL = 0.8;
        }else if(nDotL > 0.3 && nDotL < 0.7){
            nDotL = 0.4;

        }else{
            nDotL = 0.0;
        }
       vec3 diffuse = lcolor*brdf.kd*glm::max(float(0), nDotL);
       vec3 half = lray.direc + eye;
       //if( half[0] != 0.0f && half[1] != 0.0f && half[2] != 0.0f)
            half = normalize(lray.direc - eye);
       float nDotH = dot(loc.normal, half);
       vec3 specular = lcolor*brdf.ks*pow(glm::max(nDotH, 0.0f), brdf.ksh);
       //cout << specular[0] << " "  << specular[1] << " " << specular[2] << endl;
       return diffuse + specular;
    }
};
class Material{
    public:
    ShaderProgram * shaderprog;
    BRDF constantBRDF;
    Material(){
        constantBRDF = BRDF();
        shaderprog = new BlinnPhong();
    }
    Material(BRDF brdf){
        constantBRDF = brdf;
        shaderprog = new BlinnPhong();
    }
    BRDF getBRDF(){
        return constantBRDF;
    }
};

class Transformation{
    public:
    mat4 m, minvt;
    Transformation(){}
    Transformation(mat4 mat){
        m = mat;
        minvt = transpose(inverse(mat));
    }
    LocalGeo operator*(const LocalGeo & local){
       vec4 npos =  m*vec4(local.position, 1.0);
       vec4 norm = (minvt*vec4(local.normal, 0.0));
       npos = npos/npos.w;
       //if( norm[0] != 0.0f && norm[1] != 0.0f && norm[2] != 0.0f){
           return LocalGeo(vec3(npos), normalize(vec3(norm)));
       //}else{
       //    return LocalGeo(vec3(npos), vec3(norm));
       //}
    }

    Ray operator*(const Ray & ray){ 
       vec4 norig =  m*vec4(ray.origin, 1.0);
       vec4 norm = m*vec4(ray.direc, 0.0);
       norig = norig/norig.w;
       //if( norm[0] != 0.0f && norm[1] != 0.0f && norm[2] != 0.0f){
           return Ray(vec3(norig), normalize(vec3(norm)));
        //}else{
        //   return Ray(vec3(norig), vec3(norm));
        //}
    }
};

class GeometricPrimitive:public Primitive{
    public:
    Transformation objToWorld, worldToObj;
    Shape* shape;
    GeometricPrimitive(Shape* shpe){
        shape = shpe;
        mat = new Material();
    }
    GeometricPrimitive(Shape* shpe, Material* m){
        shape = shpe;
        mat = m;
    }
    bool intersect(Ray& ray, float& thit, Intersection **  in){
        Ray oray = worldToObj*ray;
        LocalGeo olocal, wlocal;
        float thit2;
        if(!shape->intersect(oray, thit, &olocal)){
            return false;
        }
        wlocal = objToWorld*olocal;
        
        thit = length(ray.origin - wlocal.position);
        *in = new Intersection(wlocal, this);
        return true;
    }
    bool intersectP(Ray& ray){
        Ray oray = worldToObj*ray;
        return shape->intersectP(oray);
    }

    BRDF getBRDF(){
        return mat->getBRDF();
    }
};

class Light{
    public:
        vec3 color;
        virtual Ray generateLightRay(LocalGeo & local) = 0;
};

class DirectionalLight:public Light{
    public:
        vec3 position;
        DirectionalLight( vec3 d ){
            position = d;
            color = vec3(1.0);
        }
        DirectionalLight( vec3 d, vec3 c){
            position = d;
            color = c;
        }
        Ray generateLightRay(LocalGeo & local){
            vec3 norm = position;
            //if( norm[0] != 0.0f && norm[1] != 0.0f && norm[2] != 0.0f){
                return Ray(local.position + float(1e-2)*local.normal, normalize(position)); 
            //}else{
            //    return Ray(local.position + local.normal*vec3(1e-8), position); 
            //}
        }
};

class PositionalLight:public Light{
    public:
        vec3 position;
        PositionalLight(){}
        PositionalLight(vec3 p){
            position = p;
            color = vec3(1.0);
        }
        PositionalLight(vec3 p, vec3 c){
            position = p;
            color = c;
        }
        Ray generateLightRay(LocalGeo & local){
            //vec3 norm = position - local.position;
            //if( norm[0] != 0.0f && norm[1] != 0.0f && norm[2] != 0.0f){
                return Ray(local.position + float(1e-2)*local.normal, normalize(position - (local.position + local.normal)*float(1e-2)));
            //}else{
            //    return Ray(local.position, position - local.position);
            //}
        }
};

class Raytracer{
    public:
    vector<Primitive*> * prims;
    vector<Light*> * lights;
    int maxdepth;
    Raytracer(){}
    
    Raytracer(vector<Primitive*> *primitives, vector<Light*> * l, int maxd){
        prims = primitives;
        lights = l;
        maxdepth = maxd;
    }

    vec3 trace(Ray& ray, int depth){
       Primitive * primitive = NULL;

       if(depth > maxdepth )
           return vec3(0.0);
       float dist_hit = 1e34; 
       float dist_near = 1e34;
       Intersection * in_test = NULL;
       Intersection * in = NULL;
       for( int i = 0; i < prims->size(); ++i){
            if((*prims)[i]->intersect(ray, dist_hit, &in_test)){
                if( dist_hit < dist_near ){
                    dist_near = dist_hit;
                    in = in_test;
                }
            }
        }
       //cout << "dude" << endl;
       if(!in)
           return vec3(0.0);
       //cout << "deg" << endl; 
       Material * mat = in->prim->mat;
       //cout << " " << color.x << color.y << color.z << endl; 
       LocalGeo localgeo = in->localgeo;
       vec3 color = vec3(0.0);
       for(int i = 0; i < lights->size(); ++i){
            Ray lray = (*lights)[i]->generateLightRay(localgeo);
            bool intersect = false;
            for( int j = 0; j < prims->size(); ++j){
                if((*prims)[j]->intersectP(lray) && (*prims)[j] != in->prim){
                    intersect = true;
                }
            }
           if( !intersect){
                 color += mat->shaderprog->shader(localgeo, lray, (*lights)[i]->color, mat->constantBRDF, ray.direc);
            }
       }
       if(mat->constantBRDF.kr > 0.0f){
           Ray reflection = Ray(localgeo.position + float(1e-2)*localgeo.normal, reflect(ray.direc, localgeo.normal)); 
           color += trace(reflection, depth+1)*(mat->constantBRDF.kr);
       }
       return color;
    }
};

class Scene{
    public:
    Camera cam;
    Raytracer raytracer;
    vector<Primitive*> prim;
    vector<Light*> lights;
    Scene(Camera camera, vector<Primitive*> primitive, vector<Light*> l, int maxdepth){
        cam = camera;
        prim = primitive;
        lights = l;
        raytracer = Raytracer( &prim, &lights, maxdepth);
    }

    void render(){
        vec3 *image = new vec3[(int)cam.width *  (int)cam.height];
        vec3 *pixel = image;
        int depth = 5;
        for( unsigned x = 0; x < cam.width; ++x){
            for( unsigned y = 0; y < cam.height; ++y, ++pixel){
                Ray p = cam.generateRay(x, y);
                *pixel = raytracer.trace(p, 0);
            }
            cout << x << endl;
        }
        std::ofstream ofs("./untitled.ppm", std::ios::out | std::ios::binary); 
    	ofs << "P6\n" << cam.width << " " << cam.height << "\n255\n"; 
    	for (unsigned i = 0; i < cam.width * cam.height; ++i) { 
    	    ofs << (unsigned char)(std::min(float(1), image[i].x) * 255) << 
    	           (unsigned char)(std::min(float(1), image[i].y) * 255) << 
    	           (unsigned char)(std::min(float(1), image[i].z) * 255); 
    	} 
    	ofs.close();  

    }
};



void matransform(stack<mat4> &transfstack, float* values) 
{
  mat4 transform = transfstack.top(); 
  vec4 valvec = vec4(values[0],values[1],values[2],values[3]); 
  vec4 newval = transform * valvec; 
  for (int i = 0; i < 4; i++) values[i] = newval[i]; 
}

void rightmultiply(const mat4 & M, stack<mat4> &transfstack) 
{
  mat4 &T = transfstack.top(); 
  T = T * M; 
}

bool readvals(stringstream &s, const int numvals, float* values) 
{
  for (int i = 0; i < numvals; i++) {
    s >> values[i]; 
    if (s.fail()) {
      cout << "Failed reading value " << i << " will skip\n"; 
      return false;
    }
  }
  return true; 
}

mat3 rotate(float degrees, vec3 axist) 
{
  // YOUR CODE FOR HW2 HERE
  // Please implement this.  Likely the same as in HW 1.
	float rads = degrees*M_RAD;
    vec3 axis = vec3(0.0);
    vec3 norm = axis;
//    if( norm[0] != 0.0f && norm[1] != 0.0f && norm[2] != 0.0f){
	    axis = normalize(axist);
  //  }
	float cosine = cos(rads);
	float sine = sin(rads);
	mat3 unchcos = mat3(cosine, 0.0, 0.0,
			    0.0, cosine, 0.0, 
			    0.0, 0.0, cosine);
	mat3 unchpar = mat3(axis[0]*axis[0], axis[0]*axis[1], axis[0]*axis[2],
			axis[1]*axis[0], axis[1]*axis[1], axis[1]*axis[2],
			axis[2]*axis[0], axis[2]*axis[1], axis[2]*axis[2]);
	mat3 chperp = mat3(0.0, axis[2]*sine, -axis[1]*sine, 
			   -axis[2]*sine, 0.0, axis[0]*sine, 
			    axis[1]*sine, -axis[0]*sine, 0.0);

	// You will change this return call
	return unchcos + (1-cosine)*unchpar + chperp;  
}

void readfile(const char* filename, Camera & cam, vector<Primitive*> & prim, vector<Light*> & lights) 
{
  string str, cmd; 
  ifstream in;
  in.open(filename);
  vec3 *verts;
  int w = 500, h = 500, maxdepth = 5, index = 0, maxverts = 0;
  vec3 diffuse = vec3(0.5), ambient = vec3(0.2), specular = vec3(1.0), emission = vec3(0.0), eyeinit, center, upinit;
  
  float shininess, fovy = 90;
  if (in.is_open()) {

    // I need to implement a matrix stack to store transforms.  
    // This is done using standard STL Templates 
    stack <mat4> transfstack; 
    transfstack.push(mat4(1.0));  // identity

    getline (in, str); 
    while (in) {
      if ((str.find_first_not_of(" \t\r\n") != string::npos) && (str[0] != '#')) {
        // Ruled out comment and blank lines 

        stringstream s(str);
        s >> cmd; 
        int i; 
        float values[10]; // Position and color for light, colors for others
        // Up to 10 params for cameras.  
        bool validinput; // Validity of input 

        // Process the light, add it to database.
        // Lighting Command
        if (cmd == "directional") {
            validinput = readvals(s, 6, values); // Position/color for lts.
            if (validinput) {
                DirectionalLight * dl = new DirectionalLight( vec3(values[0], values[1], values[2]),
                                            vec3(values[3], values[4], values[5]));
                lights.push_back(dl);
	      	}
        }else if (cmd == "point"){
             validinput = readvals(s, 6, values); // Position/color for lts.
            if (validinput) {
                 PositionalLight * dl = new PositionalLight( vec3(values[0], values[1], values[2]),
                                            vec3(values[3], values[4], values[5]));
                lights.push_back(dl);
	      	}
        // Material Commands 
        // Ambient, diffuse, specular, shininess properties for each object.
        // Filling this in is pretty straightforward, so I've left it in 
        // the skeleton, also as a hint of how to do the more complex ones.
        // Note that no transforms/stacks are applied to the colors. 

        }else if (cmd == "ambient") {
          validinput = readvals(s, 3, values); // colors 
          if (validinput) {
            for (i = 0; i < 3; i++) {
              ambient[i] = values[i]; 
            }
          }
        } else if (cmd == "diffuse") {
          validinput = readvals(s, 3, values); 
          if (validinput) {
            for (i = 0; i < 3; i++) {
              diffuse[i] = values[i]; 
            }
          }
        } else if (cmd == "specular") {
          validinput = readvals(s, 3, values); 
          if (validinput) {
            for (i = 0; i < 3; i++) {
              specular[i] = values[i]; 
            }
          }
        } else if (cmd == "emission") {
          validinput = readvals(s, 3, values); 
          if (validinput) {
            for (i = 0; i < 3; i++) {
              emission[i] = values[i]; 
            }
          }
        } else if (cmd == "shininess") {
          validinput = readvals(s, 1, values); 
          if (validinput) {
            shininess = values[0]; 
          }
        } else if (cmd == "size") {
          validinput = readvals(s,2,values); 
          if (validinput) { 
            w = (int) values[0]; h = (int) values[1]; 
          } 
        }else if (cmd == "maxdepth") {
            validinput = readvals(s, 1, values);
            if(validinput) {
                maxdepth = (int) values[0];
            }

        }else if ( cmd == "maxverts" ){
            validinput = readvals(s, 1, values);
            if(validinput){
                verts = new vec3[(int)values[0]];
                maxverts = (int)values[0];
            }

        }else if ( cmd == "vertex"){
            validinput = readvals(s, 3, values);
            if(validinput){
                if(index < maxverts){ 
                    verts[index] = vec3(values[0], values[1], values[2]);    
                    index++;
                }
            }
        } else if (cmd == "camera") {
          validinput = readvals(s,10,values); // 10 values eye cen up fov
          if (validinput) {

            // YOUR CODE FOR HW 2 HERE
            // Use all of values[0...9]
            // You may need to use the upvector fn in Transform.cpp
            // to set up correctly. 
            // Set eyeinit upinit center fovy in variables.h
	    eyeinit = vec3(values[0], values[1], values[2]);
	    center = vec3(values[3], values[4], values[5]);
	    upinit = vec3(values[6], values[7], values[8]);
	    fovy = values[9];
	    cam = Camera(eyeinit, center, upinit, fovy, w, h);

          }
        }

        // I've left the code for loading objects in the skeleton, so 
        // you can get a sense of how this works.  
        // Also look at demo.txt to get a sense of why things are done this way.
        else if (cmd == "sphere") {
            validinput = readvals(s, 4, values); 
            if (validinput) {

             Shape * shape = new Sphere(vec3(values[0], values[1], values[2]), values[3]);
             BRDF brdf;
             brdf.kd = diffuse;
             brdf.ks = specular;
             brdf.ksh = shininess;
             brdf.ke = emission;
             brdf.ka = ambient;
             Material * mat = new Material( brdf);
	        GeometricPrimitive * geoprim = new GeometricPrimitive( shape, mat);
              // Set the object's transform
            geoprim->worldToObj = Transformation(inverse(transfstack.top()));
            geoprim->objToWorld = Transformation(transfstack.top());
            prim.push_back(geoprim);
          }
        }else if(cmd == "tri"){
            validinput = readvals(s, 3, values); 
            if (validinput) {

             Shape * shape = new Triangle(verts[(int)values[0]], verts[(int)values[1]], verts[(int)values[2]]);
             BRDF brdf;
             brdf.kd = diffuse;
             brdf.ks = specular;
             brdf.ksh = shininess;
             brdf.ke = emission;
             brdf.ka = ambient;
             Material * mat = new Material( brdf);
	        GeometricPrimitive * geoprim = new GeometricPrimitive( shape, mat);
              // Set the object's transform
            geoprim->worldToObj = Transformation(inverse(transfstack.top()));
            geoprim->objToWorld = Transformation(transfstack.top());
            prim.push_back(geoprim);
          }

        }else if (cmd == "translate") {
          validinput = readvals(s,3,values); 
          if (validinput) {

            // YOUR CODE FOR HW 2 HERE.  
            // Think about how the transformation stack is affected
            // You might want to use helper functions on top of file. 
            // Also keep in mind what order your matrix is!
	    mat4 translation = mat4(1.0, 0.0, 0.0, 0.0,
				    0.0, 1.0, 0.0, 0.0,
				    0.0, 0.0, 1.0, 0.0,
				    values[0], values[1], values[2], 1.0);	
	    rightmultiply(translation, transfstack);

          }
        }
        else if (cmd == "scale") {
          validinput = readvals(s,3,values); 
          if (validinput) {

            // YOUR CODE FOR HW 2 HERE.  
            // Think about how the transformation stack is affected
            // You might want to use helper functions on top of file.  
            // Also keep in mind what order your matrix is!
	     mat4 scale = mat4(values[0], 0.0, 0.0, 0.0,
				    0.0, values[1], 0.0, 0.0,
				    0.0, 0.0, values[2], 0.0,
				    0.0, 0.0, 0.0, 1.0);	
	    rightmultiply(scale, transfstack);

          }
        }
        else if (cmd == "rotate") {
          validinput = readvals(s,4,values); 
          if (validinput) {

            // YOUR CODE FOR HW 2 HERE. 
            // values[0..2] are the axis, values[3] is the angle.  
            // You may want to normalize the axis (or in Transform::rotate)
            // See how the stack is affected, as above.  
            // Note that rotate returns a mat3. 
            // Also keep in mind what order your matrix is!

	    mat4 rot = mat4(rotate(values[3], vec3(values[0], values[1], values[2])));
	    rightmultiply(rot, transfstack);	    

          }
        }

        // I include the basic push/pop code for matrix stacks
        else if (cmd == "pushTransform") {
          transfstack.push(transfstack.top()); 
        } else if (cmd == "popTransform") {
          if (transfstack.size() <= 1) {
            cerr << "Stack has no elements.  Cannot Pop\n"; 
          } else {
            transfstack.pop(); 
          }
        }

        else {
          cerr << "Unknown Command: " << cmd << " Skipping \n"; 
        }
      }
      getline (in, str); 
    }

  } else {
    cerr << "Unable to Open Input Data File " << filename << "\n"; 
    throw 2; 
  }
}






int main(int argc, char* argv[]){
    vector<Light*> lights;

    vector<Primitive*> prim;
    Camera cam;
    readfile( argv[1], cam, prim,lights);
    Scene scene(cam, prim, lights, 5);
    scene.render();
}
