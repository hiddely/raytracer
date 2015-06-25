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

// Bounding box class, used for KD tree. Has a min vector and a max vector, containing
// the boundaries of the box.
class BBox {
public:
    Vec3Df min;
    Vec3Df max;
    
    void check(Triangle triangle);
    bool hit(Vec3Df rayOrigin, Vec3Df rayDestination);
    int longestAxis() const;
};

// Returns a number which represents the longest axis of a boundingbox:
// 0 for the x-axis, 1 for the y-axis and 2 for the z-axis
int BBox::longestAxis() const {
	int longest = 0;
	float longestLength = max[0] - min[0];
	float lengthSecond = max[1] - min[1];
	float lengthThird = max[2] - min[2];

	if (lengthSecond > longestLength) {
		longestLength = lengthSecond;
		longest = 1;
	}

	if (lengthThird > longestLength) {
		longestLength = lengthThird;
		longest = 2;
	}
	return longest;
}

// Check of the coordianates of a triangle are outside the boundingbox. If so, make the boundigbox bigger.
void BBox::check(Triangle triangle) {
    for (int i = 0; i < 3; i++) {
        Vec3Df v = MyMesh.vertices[triangle.v[i]].p;
		if (v[0] < min[0]) {
			min[0] = v[0];
		}
		else if (v[0] > max[0]) {
			max[0] = v[0];
		}

		if (v[1] < min[1]) {
			min[1] = v[1];
		}
		else if (v[1] > max[1]) {
			max[1] = v[1];
		}

		if (v[2] < min[2]) {
			min[2] = v[2];
		}            
		else if (v[2] > max[2]) {
			max[2] = v[2];
		}
    }
}	

// Method to check if a ray hit a bounding box (ray/box intersection)
bool BBox::hit(Vec3Df rayOrigin, Vec3Df rayDest) {
    
    Vec3Df rayDir = rayDest - rayOrigin;
	float tmin = 0;
	float tmax = 1000000;
	
	//Check if x direction of the ray goes between the minimum and maximum x-values of the boundingbox
    if (rayDir[0] != 0.f) {
        float tx1 = (min[0] - rayOrigin[0]) / rayDir[0];
        float tx2 = (max[0] - rayOrigin[0]) / rayDir[0];
        
        tmin = std::fmaxf(tmin, std::fminf(tx1, tx2));
        tmax = std::fminf(tmax, std::fmaxf(tx1, tx2));
    }
    
	//Check if y direction of the ray goes between the minimum and maximum x-values of the boundingbox
    if (rayDir[1] != 0.f) {
        float tx1 = (min[1] - rayOrigin[1]) / rayDir[1];
        float tx2 = (max[1] - rayOrigin[1]) / rayDir[1];
        
        tmin = std::fmaxf(tmin, std::fminf(tx1, tx2));
        tmax = std::fminf(tmax, std::fmaxf(tx1, tx2));
    }
    
	//Check if z direction of the ray goes between the minimum and maximum x-values of the boundingbox
    if (rayDir[2] != 0.f) {
        float tx1 = (min[2] - rayOrigin[2]) / rayDir[2];
        float tx2 = (max[2] - rayOrigin[2]) / rayDir[2];
        
        tmin = std::fmaxf(tmin, std::fminf(tx1, tx2));
        tmax = std::fminf(tmax, std::fmaxf(tx1, tx2));
    }
	else if (rayOrigin[0] < min[0] || rayOrigin[0] > max[0] || rayOrigin[1] < min[1] || 
		rayOrigin[1] > max[1] || rayOrigin[2] < min[2] || rayOrigin[2] > max[2]) {

        return false;
    }
    
	// return if the ray intersects the boundingbox
    return tmax >= tmin;
}

// Initializes a bounding box with one triangle
BBox initBBox(Triangle triangle) {
    BBox bbox;
    
    Vec3Df v = MyMesh.vertices[triangle.v[0]].p;
    
    Vec3Df min = Vec3Df(v[0], v[1], v[2]);
    Vec3Df max = Vec3Df(v[0], v[1], v[2]);
    
	// find the smallest x, y and z out of all three vertices of the triangle
    for (int i = 1; i < 3; i++) {
        Vec3Df vertex = MyMesh.vertices[triangle.v[i]].p;
        if (vertex[0] < min[0])
            min[0] = vertex[0];
        else if (vertex[0] > max[0])
            max[0] = vertex[0];
        if (vertex[1] < min[1])
            min[1] = vertex[1];
        else if (vertex[1] > max[1])
            max[1] = vertex[1];
        if (vertex[2] < min[2])
            min[2] = vertex[2];
        else if (vertex[2] > max[2])
            max[2] = vertex[2];
    }
    
    bbox.min = min;
    bbox.max = max;
    
    return bbox;
}

// Gets the middlepoint of a triangle
Vec3Df findMiddle(Triangle tr) {
	std::vector<Vertex> vertices = MyMesh.vertices;
	Vertex v0 = vertices[tr.v[0]];
	Vertex v1 = vertices[tr.v[1]];
	Vertex v2 = vertices[tr.v[2]];

	float xMid = (v0.p[0] + v1.p[0] + v2.p[0]) / 3;
	float yMid = (v0.p[1] + v1.p[1] + v2.p[1]) / 3;
	float zMid = (v0.p[2] + v1.p[2] + v2.p[2]) / 3;

	Vec3Df result = Vec3Df(xMid, yMid, zMid);
	return result;
}

// Class which represents a kdtree: each KDNode has a left and right child which are KDNodes, has a BoundingBox which contain the minima
// and maxima of all the triangles in the node and contains the indices of all the triangles in the node
class KDNode {
public:
    BBox bbox;
    KDNode* left;
    KDNode* right;
    std::vector<int> triangles;
    
    KDNode* build(std::vector<int>& triangles, int depth) const;
};

// Builds a kdtree
KDNode* KDNode::build(std::vector<int>& tris, int depth) const {
	// initialize nodes
    KDNode* node = new KDNode();
    node->triangles = tris;
    node->left = NULL;
    node->right = NULL;
    node->bbox = BBox();
    
	// If no triangles the tree is empty
    if (tris.size() == 0)
        return node;
	// Handle 1 triangle corner scenario
    if (tris.size() == 1) {
        node->bbox = initBBox(MyMesh.triangles[tris[0]]);
        node->left = new KDNode();
        node->right = new KDNode();
        node->left->triangles = std::vector<int>();
        node->right->triangles = std::vector<int>();
        return node;
    }
    
    // get a bounding box surrounding all the triangles
    node->bbox = initBBox(MyMesh.triangles[tris[0]]);
    
	// build the tree
    for (int i = 1; i < tris.size(); i++) {
        node->bbox.check(MyMesh.triangles[tris[i]]);
    }
    
    Vec3Df midpt(0, 0, 0);
    for (int i = 0; i < tris.size(); i++) {
        // find midpoint of all triangles
        midpt = midpt + findMiddle(MyMesh.triangles[tris[i]]) * (1.0 / tris.size());
    }
    
    std::vector<int> left_triangles;
    std::vector<int> right_triangles;
    int axis = node->bbox.longestAxis();
    for (int i = 0; i < tris.size(); i++) {
        // split triangles based on their midpoints side of avg in longest axis
        switch (axis) {
            case 0:
                midpt[0] >= findMiddle(MyMesh.triangles[tris[i]])[0] ? right_triangles.push_back(tris[i]) : left_triangles.push_back(tris[i]);
                break;
            case 1:
                midpt[1] >= findMiddle(MyMesh.triangles[tris[i]])[1] ? right_triangles.push_back(tris[i]) : left_triangles.push_back(tris[i]);
                break;
            case 2:
                midpt[2] >= findMiddle(MyMesh.triangles[tris[i]])[2] ? right_triangles.push_back(tris[i]) : left_triangles.push_back(tris[i]);
                break;
        }
    }
    
    bool done = false;
    // If either of the childs has no triangles, stop
    if (left_triangles.size() == 0 || right_triangles.size() == 0) 
        done = true;
    
    if (!done) {
        // recurse down left and right sides
        node->left = build(left_triangles, depth + 1);
        node->right = build(right_triangles, depth + 1);
    }
    
    return node;
}


KDNode* root;// = new KDNode();

bool hit(KDNode* node, const Vec3Df& rayOrigin, const Vec3Df& rayDestination, const Vec3Df& rayDirection, Vec3Df& _p, float & lastDistance, int & triangleIndex) {
    // check if ray intersects bounding box of given node
    if (node->bbox.hit(rayOrigin, rayDestination)) {
        bool hit_triangle = false;
        bool hasChildren = false;
        bool hitleft;
        bool hitright;
        
        // if either child still has triangles, recurse down both sides and check for intersections
        if (node->left != NULL) {
            if (node->left->triangles.size() > 0) {
                hasChildren = true;
                hitleft = hit(node->left, rayOrigin, rayDestination, rayDirection, _p, lastDistance, triangleIndex);
            }
        }
        if (node->right != NULL) {
            if (node->right->triangles.size() > 0) {
                hasChildren = true;
                hitright = hit(node->right, rayOrigin, rayDestination, rayDirection, _p, lastDistance, triangleIndex);
            }
        }
        if (hasChildren)
            return hitleft || hitright;
        else {
            // reached a leaf
            for (int i = 0; i < node->triangles.size(); i++) {
                int triangleIn = node->triangles[i];
                Triangle triangle = MyMesh.triangles[triangleIn];
                float t = intersectPlane(rayOrigin, rayDestination, triangle);
                
                //correct for small floating point errors in recursive steps
                if (t > 0.0005f)
                {
                    Vec3Df p = rayOrigin + (t * rayDirection);
                    // check if there is a hit with the triangle
                    if (pointInTriangle(p, triangle))
                    {
                        if (t < lastDistance)
                        {
                            hit_triangle = true;
                            lastDistance = t;
                            triangleIndex = triangleIn;
                            _p = p;
                        }
                    }
                }
            }
            if (hit_triangle)
                return true;
            return false;
        }
    }
    return false;
}

//return the t for which the ray intersects with the triangle.
float intersectPlane(const Vec3Df & rayOrigin, const Vec3Df & rayDest, const Triangle & tr)
{
    Vec3Df vec0 = MyMesh.vertices[tr.v[0]].p;
    Vec3Df vec1 = MyMesh.vertices[tr.v[1]].p;
    Vec3Df vec2 = MyMesh.vertices[tr.v[2]].p;
    
    Vec3Df normal = getNormal(tr);
    
	// If necessary turn around the normal
	if (Vec3Df::dotProduct(normal, vec0) < 0 || Vec3Df::dotProduct(normal, vec1) < 0 || Vec3Df::dotProduct(normal, vec2) < 0)
    {
		//get the right normal by making it positive.
        normal *= -1;
    }
    
    Vec3Df rayDir = rayDest - rayOrigin;
    
	float Dot = Vec3Df::dotProduct(vec2, normal);
    float distRayNorm = Vec3Df::dotProduct(rayDir, normal);
    
	// compute actual  ray triangle intersection point
	if (!(distRayNorm < 0.0001f && distRayNorm > -0.0001f)){
		float t = (Dot - (Vec3Df::dotProduct(rayOrigin, normal))) / distRayNorm;
        return t;
    }
    else
    {
        return 0.0f;
    }
}

//use this function for any preprocessing of the mesh.
void init()
{
	//load the mesh file, change path in paths.h
    MyMesh.loadMesh(MESH_PATH, true);
	MyMesh.computeVertexNormals();

	//one first move: initialize the first light source
	//at least ONE light source has to be in the scene!!!
	//here, we set it to the current location of the camera
	MyLightPositions.push_back(MyCameraPosition);
    
    //make an array of the indexes of the triangles
    std::vector<int> triangleIndexes;
    for (int i = 0; i < MyMesh.triangles.size(); i++) {
        triangleIndexes.push_back(i);
    }
    
    root = root->build(triangleIndexes, 0);
}

//return the color of your pixel.
Vec3Df performRayTracing(const Vec3Df & origin, const Vec3Df & dest)
{
    cameraOrigin = origin;
    int triangleIndex = -1;
    Vec3Df p;
    //intersect(origin, dest, triangleIndex, p);
    Vec3Df ray = dest-origin;
    float lastDistance = 10000;

    if (hit(root, origin, dest, ray, p, lastDistance, triangleIndex)) {
    //if (triangleIndex != -1) {
        // we have a hit
        //std::cout << "Ind: " << triangleIndex;
        return shade(0, triangleIndex, p, ray);
    }
    
	return Vec3Df(0.1, 0.1, 0.1);
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

//return if point p is within the triangle
bool pointInTriangle(const Vec3Df & p, const Triangle & triangle)
{
    Vec3Df v0 = MyMesh.vertices[triangle.v[0]].p;
    Vec3Df v1 = MyMesh.vertices[triangle.v[1]].p;
    Vec3Df v2 = MyMesh.vertices[triangle.v[2]].p;
    
    float a;
    float b;
	float c;
    computeBarycentric(p, v0, v1, v2, c, a, b);
    
	// compute if the point is within the triangle borders
    if (a >= 0 && a <= 1 && b >= 0 && c>= 0)
    {
        return true;
    }
    else
    {
        return false;
    }
}

/**
 Calculates shading color
 **/
Vec3Df shade(unsigned int level, const int triangleIndex, Vec3Df & p, Vec3Df ray) {
    
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
        
        Vec3Df direction = p - lightsource;
        
        int closestTriangleIndex = 1000;
        Vec3Df closestHit;
        //intersect(lightsource, direction, closestTriangleIndex, closestHit);
        float lastDistance = 10000;
        hit(root, lightsource, p, direction, closestHit, lastDistance, closestTriangleIndex);
        
        //std::cout << "Index: " << closestTriangleIndex << " and n " << triangleIndex << "Hit; " << h << std::endl;
        
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
            diffuse = lightintensity_ambient * fabs(powf(costheta, 5)) * m.Kd();
        
            // specular
            float n_inc = 8;
            Vec3Df link = ((cameraOrigin - p) - direction);
            link.normalize();
            specular = lightintensity_specular * powf(fabsf(link.dotProduct(link, n)), n_inc) * m.Ks();
            
        }
        
        // compute reflected ray
        if (m.Ns() > 25) {
            // we shine
            reflectedColor = (m.Ns() / 100.0) * traceReflectedRay(level+1, n, p, ray);
        }
        
        directLight += ambient + diffuse + specular + reflectedColor + refractedColor;

    }
    
    //return getTriangleColor(triangleIndex);

    return directLight;
    //return directLight / MyLightPositions.size() * powf(1.5, MyLightPositions.size()-1);
}

bool set = false;

//trace the reflection of the ray in p with normal n
Vec3Df traceReflectedRay(unsigned int level, const Vec3Df n, const Vec3Df p, const Vec3Df ray){
    if (level > maxLevel) {
        return Vec3Df(0, 0, 0);
    }
    
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
    float lastDistance = 10000;
    Vec3Df reflectedHit;
    yourDebugDraw();
    
    //recursively trace the reflected ray
    //intersect(p, dest, reflectedTriangleIndex, reflectedHit);
    
    if ( hit(root, p, dest, dest-p, reflectedHit, lastDistance, reflectedTriangleIndex) ) {
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

// Return the normal of a triangle
Vec3Df getNormal(const Triangle & triangle)
{
    Vec3Df edge01 = MyMesh.vertices[triangle.v[1]].p - MyMesh.vertices[triangle.v[0]].p;
    Vec3Df edge02 = MyMesh.vertices[triangle.v[2]].p - MyMesh.vertices[triangle.v[0]].p;
    Vec3Df n = Vec3Df::crossProduct(edge01, edge02);
    n.normalize();
    return n;
}

// Computes the hitpoint within vectors a, b and c
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

// Equals method for two vectors
bool equals(const Vec3Df & one, const Vec3Df & two) {
    return fabs(one[0] - two[0]) < 1 && fabs(one[1] - two[1]) < 1 && fabs(one[2] - two[2]) < 1;
}

// Method to draw the debugger screen every frame
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
	// Draw the ray for moving a light
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
	
	// t is the key character
	switch (t) {
	// Set new light
	case 'L': selectedLight = MyLightPositions.size() - 1; break;
	// switch between selected light to move 
	case 'f': selectedLight = (selectedLight + 1) % MyLightPositions.size(); break;
	// set x-axis to be the ray along the light will move
	case 'x': {
		moveLightDirection[0] = 0.1;
		moveLightDirection[1] = 0;
		moveLightDirection[2] = 0;
		red = 0; green = 0; blue = 1;
		setMoveLightRay("x"); break;
	}
	// set y-axis to be the ray along the light will move
	case 'y': {
		moveLightDirection[0] = 0;
		moveLightDirection[1] = 0.1;
		moveLightDirection[2] = 0;
		red = 1; green = 0; blue = 0;
		setMoveLightRay("y"); break;
	}
    // set z-axis to be the ray along the light will move
	case 'z': { 
		moveLightDirection[0] = 0;
		moveLightDirection[1] = 0;
		moveLightDirection[2] = 0.1;
		red = 0; green = 1; blue = 0;
		setMoveLightRay("z"); break;
	}
	// Move light position
	case 'w': MyLightPositions[selectedLight] = MyLightPositions[selectedLight] + moveLightDirection; break;
	case 's': MyLightPositions[selectedLight] = MyLightPositions[selectedLight] + moveLightDirection * -1; break;
	}	
	std::cout<<t<<" pressed! The mouse was in location "<<x<<","<<y<<"!"<<std::endl;	
}

// Computes the ray dependent move light direction  
void setMoveLightRay(std::string dir) {
	Vec3Df lastLight = MyLightPositions[selectedLight];
	moveLightRayOrigin = moveLightDirection * 100 + lastLight;
	moveLightRayDestination = moveLightDirection * -100 + lastLight;
	
	std::cout << "Move lightsource along "<< dir << " axis by pressing w and s." << std::endl;
}
