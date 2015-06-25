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
class BoundingBox {
public:
    Vec3Df min;
    Vec3Df max;
    
    void expand(Triangle triangle);
    bool hit(Vec3Df rayOrigin, Vec3Df rayDestination);
    int longest_axis() const;
};

// Returns a number which represents the longest axis of a boundingbox:
// 0 for the x-axis, 1 for the y-axis and 2 for the z-axis
int BoundingBox::longest_axis() const {
    int longest = 0;
    float lengthLongest = max[0] - min[0];
    
    if (max[1] - min[1] > lengthLongest) {
        lengthLongest = max[1] - min[1];
        longest = 1;
    }
    
    if (max[2] - min[2] > lengthLongest) {
        longest = 2;
    }
    
    return longest;
}

// Expand a boundingbax with a triangle (check for every vertex of the triangle if it extends the
// boundaries of the box)
void BoundingBox::expand(Triangle triangle) {
    for (int i = 0; i < 3; i++) {
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
}

// Method to check if a ray hit a bounding box
bool BoundingBox::hit(Vec3Df rayOrigin, Vec3Df rayDestination) {
    float tmin = std::numeric_limits<float>::min(),
    tmax = std::numeric_limits<float>::max();
    
    Vec3Df rayDir = rayDestination - rayOrigin;
    
    if (rayDir[0] != 0.0) {
        float tx1 = (min[0] - rayOrigin[0]) / rayDir[0];
        float tx2 = (max[0] - rayOrigin[0]) / rayDir[0];
        
        tmin = std::max(tmin, std::min(tx1, tx2));
        tmax = std::min(tmax, std::max(tx1, tx2));
    }
    else if (rayOrigin[0] < min[0] || rayOrigin[0] > max[0]) {
        return false;
    }
    
    if (rayDir[1] != 0.0) {
        float tx1 = (min[1] - rayOrigin[1]) / rayDir[1];
        float tx2 = (max[1] - rayOrigin[1]) / rayDir[1];
        
        tmin = std::max(tmin, std::min(tx1, tx2));
        tmax = std::min(tmax, std::max(tx1, tx2));
    }
    else if (rayOrigin[1] < min[1] || rayOrigin[1] > max[1]) {
        return false;
    }
    
    if (rayDir[2] != 0.0) {
        float tx1 = (min[2] - rayOrigin[2]) / rayDir[2];
        float tx2 = (max[2] - rayOrigin[2]) / rayDir[2];
        
        tmin = std::max(tmin, std::min(tx1, tx2));
        tmax = std::min(tmax, std::max(tx1, tx2));
    }
    else if (rayOrigin[2] < min[2] || rayOrigin[2] > max[2]) {
        return false;
    }
    
    return tmax >= tmin;
}

// Initializes a bounding box with one triangle
BoundingBox getBB(Triangle triangle) {
    BoundingBox bb;
    
    Vec3Df vertex = MyMesh.vertices[triangle.v[0]].p;
    
    Vec3Df min = Vec3Df(vertex[0], vertex[1], vertex[2]);
    Vec3Df max = Vec3Df(vertex[0], vertex[1], vertex[2]);
    
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
    
    bb.min = min;
    bb.max = max;
    
    return bb;
}

// Gets the middlepoint of a triangle
Vec3Df getMidPoint(Triangle triangle) {
    Vec3Df v0 = MyMesh.vertices[triangle.v[0]].p;
    Vec3Df v1 = MyMesh.vertices[triangle.v[1]].p;
    Vec3Df v2 = MyMesh.vertices[triangle.v[2]].p;
    
    float xMiddle = (v0[0] + v1[0] + v2[0]) / 3;
    float yMiddle = (v0[1] + v1[1] + v2[1]) / 3;
    float zMiddle = (v0[2] + v1[2] + v2[2]) / 3;
    
    return Vec3Df(xMiddle, yMiddle, zMiddle);
}

// Class which represents a kdtree: each KDNode has a left and right child which are KDNodes, has a BoundingBox which contain the minima
// and maxima of all the triangles in the node and contains the indices of all the triangles in the node
class KDNode {
public:
    BoundingBox bbox;
    KDNode* left;
    KDNode* right;
    std::vector<int> triangles;
    
    KDNode* build(std::vector<int>& triangles, int depth) const;
};

// Builds a kdtree
KDNode* KDNode::build(std::vector<int>& triangles, int depth) const {
    KDNode* node = new KDNode();
    node->triangles = triangles;
    node->left = NULL;
    node->right = NULL;
    node->bbox = BoundingBox();
    
    if (triangles.size() == 0)
        return node;
    
    if (triangles.size() == 1) {
        node->bbox = getBB(MyMesh.triangles[triangles[0]]);
        node->left = new KDNode();
        node->right = new KDNode();
        node->left->triangles = std::vector<int>();
        node->right->triangles = std::vector<int>();
        return node;
    }
    
    // get a bounding box surrounding all the triangles
    node->bbox = getBB(MyMesh.triangles[triangles[0]]);
    
    for (int i = 1; i < triangles.size(); i++) {
        node->bbox.expand(MyMesh.triangles[triangles[i]]);
    }
    
    Vec3Df midpt(0, 0, 0);
    for (int i = 0; i < triangles.size(); i++) {
        // find midpoint of all triangles
        midpt = midpt + (getMidPoint(MyMesh.triangles[triangles[i]]) * (1.0 / triangles.size()));
    }
    
    std::vector<int> left_triangles;
    std::vector<int> right_triangles;
    int axis = node->bbox.longest_axis();
    for (int i = 0; i < triangles.size(); i++) {
        // split triangles based on their midpoints side of avg in longest axis
        switch (axis) {
            case 0:
                midpt[0] >= getMidPoint(MyMesh.triangles[triangles[i]])[0] ? right_triangles.push_back(triangles[i]) : left_triangles.push_back(triangles[i]);
                break;
            case 1:
                midpt[1] >= getMidPoint(MyMesh.triangles[triangles[i]])[1] ? right_triangles.push_back(triangles[i]) : left_triangles.push_back(triangles[i]);
                break;
            case 2:
                midpt[2] >= getMidPoint(MyMesh.triangles[triangles[i]])[2] ? right_triangles.push_back(triangles[i]) : left_triangles.push_back(triangles[i]);
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

// Data for the recursive calls for hit
struct intersection_state {
    float intersectionDepth = 1000000;
    int triangleIndex = -1;
};
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

//return the position along the ray where the plane spanned by the triangle is intersected
float intersectPlane(const Vec3Df & rayOrigin, const Vec3Df & rayDestination, const Triangle & triangle)
{
    Vec3Df v0 = MyMesh.vertices[triangle.v[0]].p;
    Vec3Df v1 = MyMesh.vertices[triangle.v[1]].p;
    Vec3Df v2 = MyMesh.vertices[triangle.v[2]].p;
    
    Vec3Df n = getNormal(triangle);
    
    if (Vec3Df::dotProduct(n, v0) < 0 || Vec3Df::dotProduct(n, v1) < 0 || Vec3Df::dotProduct(n, v2) < 0)
    {
        n *= -1;
    }
    
    Vec3Df rayDirection = rayDestination - rayOrigin;
    
    float D = Vec3Df::dotProduct(v2, n);
    float distanceFromRayToNormal = Vec3Df::dotProduct(rayDirection, n);
    
    if (!(distanceFromRayToNormal < 0.0001f && distanceFromRayToNormal > -0.0001f)){
        float t = (D - (Vec3Df::dotProduct(rayOrigin, n))) / distanceFromRayToNormal;
        return t;
    }
    else
    {
        return 0.0f;
    }
}

// Compute barycentric coordinates (u, v, w) for
// point p with respect to triangle (a, b, c)
void Barycentric(const Vec3Df & p, const Vec3Df & a, const Vec3Df & b, const Vec3Df & c, float &v, float &w)
{
    Vec3Df v0 = b - a, v1 = c - a, v2 = p - a;
    float d00 = Vec3Df::dotProduct(v0, v0);
    float d01 = Vec3Df::dotProduct(v0, v1);
    float d11 = Vec3Df::dotProduct(v1, v1);
    float d20 = Vec3Df::dotProduct(v2, v0);
    float d21 = Vec3Df::dotProduct(v2, v1);
    float denom = d00 * d11 - d01 * d01;
    v = (d11 * d20 - d01 * d21) / denom;
    w = (d00 * d21 - d01 * d20) / denom;
}

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

//return if point p is within the triangle
bool pointInTriangle(const Vec3Df & p, const Triangle & triangle)
{
    Vec3Df v0 = MyMesh.vertices[triangle.v[0]].p;
    Vec3Df v1 = MyMesh.vertices[triangle.v[1]].p;
    Vec3Df v2 = MyMesh.vertices[triangle.v[2]].p;
    
    float a;
    float b;
    Barycentric(p, v0, v1, v2, a, b);
    
    if (a >= 0 && a <= 1 && b >= 0 && (a + b) <= 1)
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
