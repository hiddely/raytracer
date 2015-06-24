#include <stdio.h>
#ifdef WIN32
#include <windows.h>
#include <GL/glut.h>
#endif
#ifdef __APPLE__
#include <GLUT/glut.h>
#endif
#include "raytracing.h"
#include "paths.h"


//temporary variables
//these are only used to illustrate 
//a simple debug drawing. A ray 
Vec3Df testRayOrigin;
Vec3Df testRayDestination;
Vec3Df moveLightRayOrigin; 
Vec3Df moveLightRayDestination;
Vec3Df moveLightDirection = Vec3Df(0.1, 0, 0);
double red;
double green;
double blue;
int selectedLight = 0;
Vec3Df cameraOrigin;
int maxLevel = 2;

//use this function for any preprocessing of the mesh.
void init()
{
	//load the mesh file
	//please realize that not all OBJ files will successfully load.
	//Nonetheless, if they come from Blender, they should, if they 
	//are exported as WavefrontOBJ.
	//PLEASE ADAPT THE LINE BELOW TO THE FULL PATH OF THE dodgeColorTest.obj
	//model, e.g., "C:/temp/myData/GraphicsIsFun/dodgeColorTest.obj", 
	//otherwise the application will not load properly
    MyMesh.loadMesh(MESH_PATH, true);
	MyMesh.computeVertexNormals();

	//one first move: initialize the first light source
	//at least ONE light source has to be in the scene!!!
	//here, we set it to the current location of the camera
	MyLightPositions.push_back(MyCameraPosition);
}

//return the color of your pixel.
Vec3Df performRayTracing(const Vec3Df & origin, const Vec3Df & dest)
{
    cameraOrigin = origin;
    int triangleIndex;
    Vec3Df hit;
    intersect(origin, dest, triangleIndex, hit);
    if (triangleIndex != -1) {
        // we have a hit
        Vec3Df ray = dest-origin;
        return shade(0, triangleIndex, hit, ray);
    }
    
	return Vec3Df(0, 0, 0);
}

/**
 Checks if there is an intersection between the ray and the triangles, returns closest triangleIndex
 **/
void intersect(const Vec3Df & origin, const Vec3Df & dest, int & triangleIndex, Vec3Df & hit) {
    std::vector<Vertex> vertices = MyMesh.vertices;
    float lastDistance = 100000000; // big number
    triangleIndex = -1; // -1 means no hit
    for(std::vector<int>::size_type i = 0; i != MyMesh.triangles.size(); i++) {
        /* std::cout << *it; ... */
        // single triangle
        Triangle triangle = MyMesh.triangles[i];
        Vertex v0 = vertices[triangle.v[0]];
        Vertex v1 = vertices[triangle.v[1]];
        Vertex v2 = vertices[triangle.v[2]];
        
        // d in n
        
        Vec3Df normal = surfaceNormalTriangle(v0, v1, v2);

        float ndotd = normal.dotProduct(normal, dest);
        
        // calculate if our ray has a non-zero dot product with the normal
        if (ndotd != 0) {
            Vec3Df D = normal.projectOntoVector(vertices[triangle.v[0]].p, normal);
            float odotn = normal.dotProduct(normal, origin);
            
            float t = (D.getLength()-odotn)/ndotd;
            
            if (t < 0.0005f) {
                continue;
            }
            
            // now we have t, check if we are inside the triangle
            Vec3Df p = origin + t * dest;
            // but p is also ... = a * v0 + b * v1 + (1-a-b) * v2
            
            float a, b, c;
            computeBarycentric(p, v0.p, v1.p, v2.p, a, b, c);
            
            if (!(a < 0 || a > 1 || b < 0 || a+b > 1)) {
                // we are inside triangle

                if (lastDistance > t) {
                    triangleIndex = i;
                    lastDistance = t;
                    hit = p;
                }
            }
        }
        
    }
}

/**
 Calculates shading color
 **/
Vec3Df shade(unsigned int level, const unsigned int triangleIndex, Vec3Df & hit, Vec3Df ray) {
    
    if (level > maxLevel) {
        return Vec3Df(0, 0, 0);
    }
    
    Vec3Df directLight = Vec3Df(0, 0, 0);
    
    Material m = getTriangleMaterial(triangleIndex);
    
    Triangle triangle = MyMesh.triangles[triangleIndex];
    
    float lightintensity_ambient = 1.1;
    float lightintensity_specular = 0.8;
    
    Vec3Df ambient = 0.5 * m.Kd();
    
    // for each light
    for(std::vector<int>::size_type i = 0; i != MyLightPositions.size(); i++) {
        Vec3Df lightsource = MyLightPositions[i];
        
        Vec3Df direction = hit - lightsource;
        
        int closestTriangleIndex;
        Vec3Df closestHit;
        intersect(lightsource, direction, closestTriangleIndex, closestHit);
        
        Vec3Df n = getNormal(triangle);
        
        // calculate ambient term
        Vec3Df diffuse = Vec3Df(0, 0, 0);
        Vec3Df specular = Vec3Df(0, 0, 0);
        Vec3Df reflectedColor = Vec3Df(0, 0, 0);
        Vec3Df refractedColor = Vec3Df(0, 0, 0);
        
        if (triangleIndex == closestTriangleIndex) {
            // let there be light
            
            // calculate diffuse term
            direction.normalize();
            float costheta = Vec3Df::dotProduct(n, direction);
            diffuse = lightintensity_ambient * fabs(powf(costheta, 1)) * m.Kd();
        
            // specular
            float n_inc = 8;
            Vec3Df link = ((cameraOrigin - hit) - direction);
            link.normalize();
            specular = lightintensity_specular * powf(fabsf(link.dotProduct(link, n)), n_inc) * m.Ks();
            
        }
        
        // compute reflected ray
        if (m.Ns() > 25) {
            // we shine
            reflectedColor = traceReflectedRay(level, n, hit, ray);
        }
        
        directLight += ambient + diffuse + specular + reflectedColor + refractedColor;

    }
    
    //return getTriangleColor(triangleIndex);

    return directLight / MyLightPositions.size() * powf(1.5, MyLightPositions.size()-1);
}

bool set = false;

//trace the reflection of the ray in p with normal n
Vec3Df traceReflectedRay(unsigned int level, const Vec3Df n, const Vec3Df p, const Vec3Df ray){
    Vec3Df v = ray;
    v.normalize();
    Vec3Df reflectedRay = v - 2 * (Vec3Df::dotProduct(n, v)) * n;
    Vec3Df dest = p + reflectedRay;
    if (!set) {
        testRayOrigin = p;
        testRayDestination = dest;
        set = true;
    }
    int reflectedTriangleIndex;
    Vec3Df reflectedHit;
    yourDebugDraw();
    
    //recursively trace the reflected ray
    intersect(p, dest, reflectedTriangleIndex, reflectedHit);

    if (reflectedTriangleIndex != -1) {
        return shade(level+1, reflectedTriangleIndex, reflectedHit, reflectedHit-p);
    }
    return Vec3Df(0, 0, 0);
}

Vec3Df getTriangleColor(const unsigned int triangleIndex) {
    Material m = MyMesh.materials[MyMesh.triangleMaterials[triangleIndex]];
    
    // for now return diffuse value
    return m.Kd();
}

Material getTriangleMaterial(const unsigned int triangleIndex) {
    return MyMesh.materials[MyMesh.triangleMaterials[triangleIndex]];
}

// calculates the surface normal vector n
Vec3Df surfaceNormalTriangle(const Vertex & v0, const Vertex & v1, const Vertex & v2) {
    
    Vec3Df product = v0.p.crossProduct((v0.p-v2.p), (v1.p-v2.p));
    product.normalize();
    
    if (product.dotProduct(product, v0.p) < 0) {
        product = -1 * product;
    }
    
    return product;
}

Vec3Df getNormal(const Triangle & triangle)
{
    Vec3Df edge01 = MyMesh.vertices[triangle.v[1]].p - MyMesh.vertices[triangle.v[0]].p;
    Vec3Df edge02 = MyMesh.vertices[triangle.v[2]].p - MyMesh.vertices[triangle.v[0]].p;
    Vec3Df n = Vec3Df::crossProduct(edge01, edge02);
    n.normalize();
    return n;
}

void computeBarycentric(Vec3Df p, Vec3Df a, Vec3Df b, Vec3Df c, float &u, float &v, float &w)
{
    Vec3Df v0 = b - a, v1 = c - a, v2 = p - a;
    float d00 = v0.dotProduct(v0, v0);
    float d01 = v0.dotProduct(v0, v1);
    float d11 = v0.dotProduct(v1, v1);
    float d20 = v0.dotProduct(v2, v0);
    float d21 = v0.dotProduct(v2, v1);
    float denom = d00 * d11 - d01 * d01;
    v = (d11 * d20 - d01 * d21) / denom;
    w = (d00 * d21 - d01 * d20) / denom;
    u = 1.0f - v - w;
}

bool equals(const Vec3Df & one, const Vec3Df & two) {
    return fabs(one[0] - two[0]) < 1 && fabs(one[1] - two[1]) < 1 && fabs(one[2] - two[2]) < 1;
}

/**
 Blends two colors with gradient a
 **/
Vec3Df blendColors(const Vec3Df & c1, const Vec3Df & c2, const float a)
{
    if (a >= 0 && a <= 1) {
        return (c1 * a) + (c2 * (1 - a));
    }
    else if (a > 0) {
        return c2;
    }
    return c1;
}

void yourDebugDraw()
{
	//draw open gl debug stuff
	//this function is called every frame

	//let's draw the mesh
	MyMesh.draw();
	
	//let's draw the lights in the scene as points
	glPushAttrib(GL_ALL_ATTRIB_BITS); //store all GL attributes
	glDisable(GL_LIGHTING);
	glColor3f(1,1,1);
	glPointSize(10);
	glBegin(GL_POINTS);
	for (int i=0;i<MyLightPositions.size();++i)
		glVertex3fv(MyLightPositions[i].pointer());
	glEnd();
	glPopAttrib();//restore all GL attributes
	//The Attrib commands maintain the state. 
	//e.g., even though inside the two calls, we set
	//the color to white, it will be reset to the previous 
	//state after the pop.


	//as an example: we draw the test ray, which is set by the keyboard function
	glPushAttrib(GL_ALL_ATTRIB_BITS);
	glDisable(GL_LIGHTING);
	glBegin(GL_LINES);
	glColor3f(0,1,1);
	glVertex3f(testRayOrigin[0], testRayOrigin[1], testRayOrigin[2]);
	glColor3f(0,0,1);
	glVertex3f(testRayDestination[0], testRayDestination[1], testRayDestination[2]);
	glEnd();
	glBegin(GL_LINES);
	glColor3f(red, green, blue);
	glVertex3f(moveLightRayOrigin[0], moveLightRayOrigin[1], moveLightRayOrigin[2]);
	glColor3f(red, green, blue);
	glVertex3f(moveLightRayDestination[0], moveLightRayDestination[1], moveLightRayDestination[2]);
	glEnd();
	glPointSize(10);
	glBegin(GL_POINTS);
	glVertex3fv(MyLightPositions[0].pointer());
	glEnd();
	glPopAttrib();
	
	//draw whatever else you want...
	////glutSolidSphere(1,10,10);
	////allows you to draw a sphere at the origin.
	////using a glTranslate, it can be shifted to whereever you want
	////if you produce a sphere renderer, this 
	////triangulated sphere is nice for the preview
}


//yourKeyboardFunc is used to deal with keyboard input.
//t is the character that was pressed
//x,y is the mouse position in pixels
//rayOrigin, rayDestination is the ray that is going in the view direction UNDERNEATH your mouse position.
//
//A few keys are already reserved: 
//'L' adds a light positioned at the camera location to the MyLightPositions vector
//'l' modifies the last added light to the current 
//    camera position (by default, there is only one light, so move it with l)
//    ATTENTION These lights do NOT affect the real-time rendering. 
//    You should use them for the raytracing.
//'r' calls the function performRaytracing on EVERY pixel, using the correct associated ray. 
//    It then stores the result in an image "result.ppm".
//    Initially, this function is fast (performRaytracing simply returns 
//    the target of the ray - see the code above), but once you replaced 
//    this function and raytracing is in place, it might take a 
//    while to complete...
void yourKeyboardFunc(char t, int x, int y, const Vec3Df & rayOrigin, const Vec3Df & rayDestination)
{

	//here, as an example, I use the ray to fill in the values for my upper global ray variable
	//I use these variables in the debugDraw function to draw the corresponding ray.
	//try it: Press a key, move the camera, see the ray that was launched as a line.
	//testRayOrigin=rayOrigin;
	//testRayDestination=rayDestination;
    set = false;
	Vec3Df offset;
	
	switch (t) {
	case 'L': selectedLight = MyLightPositions.size() - 1; break;
	case 'f': selectedLight = (selectedLight + 1) % MyLightPositions.size(); break;
	case 'x': {
		moveLightDirection[0] = 0.1;
		moveLightDirection[1] = 0;
		moveLightDirection[2] = 0;
		red = 0; green = 0; blue = 1;
		setMoveLightRay("x"); break;
	}
	case 'y': {
		moveLightDirection[0] = 0;
		moveLightDirection[1] = 0.1;
		moveLightDirection[2] = 0;
		red = 1; green = 0; blue = 0;
		setMoveLightRay("y"); break;
	}
	case 'z': { 
		moveLightDirection[0] = 0;
		moveLightDirection[1] = 0;
		moveLightDirection[2] = 0.1;
		red = 0; green = 1; blue = 0;
		setMoveLightRay("z"); break;
	}
	case 'w': MyLightPositions[selectedLight] = MyLightPositions[selectedLight] + moveLightDirection; break;
	case 's': MyLightPositions[selectedLight] = MyLightPositions[selectedLight] + moveLightDirection * -1; break;
	}	
	std::cout<<t<<" pressed! The mouse was in location "<<x<<","<<y<<"!"<<std::endl;	
}

void setMoveLightRay(std::string dir) {
	Vec3Df lastLight = MyLightPositions[selectedLight];
	moveLightRayOrigin = moveLightDirection * 100 + lastLight;
	moveLightRayDestination = moveLightDirection * -100 + lastLight;
	
	std::cout << "Move lightsource along "<< dir << " axis by pressing w and s." << std::endl;
}
