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
    std::vector<Vertex> vertices = MyMesh.vertices;
    float lastDistance = 100000000; // big number
    Vec3Df lastColor = Vec3Df(0, 0, 0);
    for(std::vector<int>::size_type i = 0; i != MyMesh.triangles.size(); i++) {
        /* std::cout << *it; ... */
        // single triangle
        Triangle triangle = MyMesh.triangles[i];
        Vertex v0 = vertices[triangle.v[0]];
        Vertex v1 = vertices[triangle.v[1]];
        Vertex v2 = vertices[triangle.v[2]];

        
        // d in n
        
        Vec3Df normal = surfaceNormalTriangle(v0, v1, v2);
        if (normal.dotProduct(normal, v0.p) < 0) {
            normal = -1 * normal;
        }
        float ndotd = normal.dotProduct(normal, dest);
        
        // calculate if our ray has a non-zero dot product with the normal
        if (ndotd != 0) {
            Vec3Df D = normal.projectOntoVector(vertices[triangle.v[0]].p, normal);
            float odotn = normal.dotProduct(normal, origin);
            
            float t = (D.getLength()-odotn)/ndotd;
            
            // now we have t, check if we are inside the triangle
            Vec3Df p = origin + t * dest;
            // but p is also ... = a * v0 + b * v1 + (1-a-b) * v2
            
            float a, b, c;
            computeBarycentric(p, v0.p, v1.p, v2.p, a, b, c);
            
            if (!(a < 0 || a > 1 || b < 0 || a+b > 1)) {
                // we are inside triangle
                
                if (lastDistance > t) {
                    lastColor = getTriangleColor(i);
                    lastDistance = t;
                }
            }
            
            //std::cout << "T:"<<t<<"\n";
            
            //return Vec3Df;
        }
        
    }
    
	return lastColor;
}


// calculates the surface normal vector n
Vec3Df surfaceNormalTriangle(const Vertex & v0, const Vertex & v1, const Vertex & v2) {
    
    Vec3Df product = v0.p.crossProduct((v0.p-v2.p), (v1.p-v2.p));
    product.normalize();
    
    return product;
}

Vec3Df getTriangleColor(const unsigned int triangleIndex) {
    Material m = MyMesh.materials[MyMesh.triangleMaterials[triangleIndex]];
    
    // for now return ambient value
    return m.Kd();
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
	testRayOrigin=rayOrigin;	
	testRayDestination=rayDestination;
	
	// do here, whatever you want with the keyboard input t.
	
	//...
	
	
	std::cout<<t<<" pressed! The mouse was in location "<<x<<","<<y<<"!"<<std::endl;	
}
