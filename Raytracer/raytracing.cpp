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


/*class BBox {
public:
	Vec3Df min;
	Vec3Df max;

	//BBox();
	void check(Triangle tr);
	bool hit(Vec3Df rayOrigin, Vec3Df rayDest);
	int longestAxis() const;
	void init(Triangle tr);
};*/

/*BBox::BBox(void) {
	std::cout << "New Bbox object created" << std::endl;
}*/

bool BBox::hit(Vec3Df rayOrigin, Vec3Df rayDest) {
	//X0 is origin van de ray en n is de direction
	Vec3Df rayDirection = rayDest - rayOrigin;
	float tmin, tmax;

	//divide by 0 is not possible.
	if (rayDirection[0] != 0.f) {
		float tx1 = (min[0] - rayOrigin[0]) / (rayDirection)[0];
		float tx2 = (max[0] - rayOrigin[0]) / (rayDirection)[0];

		tmin = min(tx1, tx2);
		tmax = max(tx1, tx2);
	}

	//if one of the coordinates is outside the coordinates of the box, there is no hit.
	//I do not understand why the coordinates of the origin are taken though and not those of e.g. the destination.
	else if (rayOrigin[0] < min[0] || rayOrigin[0] > max[0]) return false;

	//divide by 0 is not possible.
	if (rayDirection[1] != 0.f) {
		float ty1 = (min[1] - rayOrigin[1]) / rayDirection[1];
		float ty2 = (max[1] - rayOrigin[1]) / (rayDirection)[1];

		tmin = max(tmin, min(ty1, ty2));
		tmax = min(tmax, max(ty1, ty2));
	}

	//coordinates are outside the box.
	else if (rayOrigin[1] < min[1] || rayOrigin[1] > max[1]) return false;

	if (rayDirection[2] != 0.f) {

		float tz1 = (min[2] - rayOrigin[2]) / (rayDirection)[2];
		float tz2 = (max[2] - rayOrigin[2]) / (rayDirection)[2];

		tmin = max(tmin, min(tz1, tz2));
		tmax = min(tmax, max(tz1, tz2));
		std::cout << tmax << std::endl;
	}

	else if (rayOrigin[2] < min[2] || rayOrigin[2] > max[2]) return false;
	
	bool result = tmax >= tmin;
	//std::cout << "Boolean box is hit: "<< result << std::endl;
	return result;
}

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

BBox initBB(Triangle tr) {
	BBox bb;

	Vec3Df vertex = MyMesh.vertices[tr.v[0]].p;

	Vec3Df mini = Vec3Df(vertex[0], vertex[1], vertex[2]);
	Vec3Df maxi = Vec3Df(vertex[0], vertex[1], vertex[2]);

	for (int i = 1; i < 3; i++) {
		Vec3Df vertex = MyMesh.vertices[tr.v[i]].p;
		if (vertex[0] < mini[0])
			mini[0] = vertex[0];
		else if (vertex[0] > maxi[0])
			maxi[0] = vertex[0];
		if (vertex[1] < mini[1])
			mini[1] = vertex[1];
		else if (vertex[1] > maxi[1])
			maxi[1] = vertex[1];
		if (vertex[2] < mini[2])
			mini[2] = vertex[2];
		else if (vertex[2] > maxi[2])
			maxi[2] = vertex[2];
	}

	bb.min = mini;
	bb.max = maxi;

	return bb;
}

void BBox::check(Triangle tr) {

	std::vector<Vertex> vertices = MyMesh.vertices;
	Vertex v0 = vertices[tr.v[0]];
	Vertex v1 = vertices[tr.v[1]];
	Vertex v2 = vertices[tr.v[2]];

	float minX = v0.p[0];
	float maxX = v0.p[0];

	float minY = v0.p[1];
	float maxY = v0.p[1];

	float minZ = v0.p[2];
	float maxZ = v0.p[2];

	float v1X = v1.p[0];
	float v1Y = v1.p[1];
	float v1Z = v1.p[2];

	float v2X = v2.p[0];
	float v2Y = v2.p[1];
	float v2Z = v2.p[2];

	minX = min(minX, v1X);
	minX = min(minX, v2X);

	minY = min(minY, v1Y);
	minY = min(minY, v2Y);

	minZ = min(minZ, v1Z);
	minZ = min(minZ, v2Z);

	maxX = max(maxX, v1X);
	maxX = max(maxX, v2X);

	maxY = max(maxY, v1Y);
	maxY = max(maxY, v2Y);

	maxZ = max(maxZ, v1Z);
	maxZ = max(maxZ, v2Z);

	Vec3Df min = min;
	Vec3Df max = max;

	float oldMinX = min.p[0];
	float oldMinY = min.p[1];
	float oldMinZ = min.p[2];

	min[0] = min(oldMinX, minX);
	min[1] = min(oldMinY, minY);
	min[2] = min(oldMinZ, minZ);

	float oldMaxX = max.p[0];
	float oldMaxY = max.p[1];
	float oldMaxZ = max.p[2];

	max[0] = max(oldMaxX, maxX);
	max[1] = max(oldMaxY, maxY);
	max[2] = max(oldMaxZ, maxZ);
}

class KDNode {
public: 
	BBox bbox;
	KDNode* left;
	KDNode* right;
	std::vector<int> triangles;

	//KDNode();
	KDNode* build(std::vector<int> tris, int depth) const;
};

//KDNode::KDNode(void) {
	//std::cout << "KDNode created" << std::endl;
//}

Vec3Df findMidPoint(Triangle tr) {
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

KDNode* KDNode::build(std::vector<int> tris, int depth) const{
	KDNode* node = new KDNode();
	node->triangles = tris;
	node->left = NULL;
	node->right = NULL;
	node->bbox = BBox();

	if (tris.size() == 0) {
		return node;
	}

	if (tris.size() == 1) {
		node->bbox.check(MyMesh.triangles[tris[0]]);
		node->left = new KDNode();
		node->right = new KDNode();
		node->left->triangles = std::vector<int>();
		node->right->triangles = std::vector<int>();
		return node;
	}

	Triangle t1 = MyMesh.triangles[tris[0]];
	node->bbox.check(t1);
	std::cout << node->bbox.min[0] << std::endl;

	for (int i = 1; i < tris.size(); i++) {
		node->bbox.check(MyMesh.triangles[tris[i]]);
	}

	Vec3Df midPoint = Vec3Df(255, 0, 0);
	for (int i = 0; i < tris.size(); i++) {
		//find midpoint 
		midPoint = midPoint + (findMidPoint(MyMesh.triangles[tris[i]]) * (1.0 / tris.size()));
	}

	std::vector<int> left_tris;
	std::vector<int> right_tris;
	float midX, midY, midZ, midXTriangle, midYTriangle, midZTriangle;
	int longestAxis = node->bbox.longestAxis();
	for (int i = 0; i < tris.size(); i++) {
		switch (longestAxis) {
		case 0:
			midX = midPoint.p[0];
			midXTriangle = findMidPoint(MyMesh.triangles[tris[i]]).p[0];
			midX >= midXTriangle ? right_tris.push_back(tris[i]) : left_tris.push_back(tris[i]);
			break;
		case 1:
			midY = midPoint.p[1];
			midYTriangle = findMidPoint(MyMesh.triangles[tris[i]]).p[1];
			midY >= midYTriangle ? right_tris.push_back(tris[i]) : left_tris.push_back(tris[i]);
			break;
		case 2:
			midZ = midPoint.p[2];
			midZTriangle = findMidPoint(MyMesh.triangles[tris[i]]).p[2];
			midZ >= midZTriangle ? right_tris.push_back(tris[i]) : left_tris.push_back(tris[i]);
			break;
		}
	}

	if (left_tris.size() == 0 && right_tris.size() > 0) left_tris = right_tris;
	if (right_tris.size() == 0 && left_tris.size() > 0) right_tris = left_tris;

	int matches = 0;
	for (int i = 0; i < left_tris.size(); i++) {
		for (int k = 0; k < right_tris.size(); k++) {
			if (left_tris[i] == right_tris[k]) {
				matches++;
			}
		}
	}

	//If 50% of the triangles match, do not subdivide anymore.
	if ((float)matches / left_tris.size() < 0.5 && (float)matches / right_tris.size() < 0.5) {
		//recurse down left and right sides
		node->left = build(left_tris, depth + 1);
		node->right = build(right_tris, depth + 1);
	}

	//else might not be necessary.
	else {
		node->left = new KDNode();
		node->right = new KDNode();
		node->left->triangles = std::vector<int>();
		node->right->triangles = std::vector<int>();
	}

	return node;
}

bool hit(KDNode* node, const Vec3Df& rayOrigin, const Vec3Df& rayDest, Vec3Df hitPoint, int triangleIndex) {
	if (node->bbox.hit(rayOrigin, rayDest)) {
		bool hit_tri = false;
		bool hasChildren = false;
		bool hitleft;
		bool hitright;
		Vec3Df hit_pt;

		if (node->left != NULL) {
			if (node->left->triangles.size() > 0) {
				hasChildren = true;
				hitleft = hit(node->left, rayOrigin, rayDest, hitPoint, triangleIndex);
			}
		}

		if (node->right != NULL) {
			if (node->right->triangles.size() > 0) {
				hasChildren = true;
				hitright = hit(node->right, rayOrigin, rayDest, hitPoint, triangleIndex);
			}
		}

		if (hasChildren) { return hitleft || hitright; }
		/*//If either child still has triangles, recurse down both sides and check for intersections
		if (node->left->triangles.size() > 0 || node->right->triangles.size() > 0) {
			bool hitleft = hit(node->left, rayOrigin, rayDest, hitPoint, triangleIndex);
			bool hitright = hit(node->right, rayOrigin, rayDest, hitPoint, triangleIndex);
			return hitleft || hitright;
		} */
		else {
			//we reached a leaf
			for (int i = 0; i < node->triangles.size(); i++) {
				int triangleInputIndex = node->triangles[i];
				std::cout << triangleInputIndex << std::endl;
				Triangle tr = MyMesh.triangles[triangleInputIndex];
				float intersectionDepth = 0;
				float intersectionDepthTemp = 10000000;
				Vec3Df hitPointTemp = Vec3Df(0, 0, 0);
				if (intersectTriangle(tr, rayOrigin, rayDest, hitPointTemp, intersectionDepth)) {
					if (intersectionDepth < intersectionDepthTemp) {
						hit_tri = true;
						hitPoint = hitPointTemp;
						triangleIndex = triangleInputIndex;
					}

				}
			}
			if (hit_tri) {
				return true;
			}
			return false;
		}
	}
	return false;
}

bool intersectTriangle(Triangle tr, const Vec3Df & rayOrigin, const Vec3Df & rayDest, Vec3Df& hitPoint, float intersectionDepth) {
	std::vector<Vertex> vertices = MyMesh.vertices;
	bool hit = false;
	float lastDistance = 10000000; //big number
	Vertex v0 = vertices[tr.v[0]];
	Vertex v1 = vertices[tr.v[1]];
	Vertex v2 = vertices[tr.v[2]];

	// d in n

	Vec3Df normal = surfaceNormalTriangle(v0, v1, v2);
	
	float ndotd = normal.dotProduct(normal, rayDest);

	// calculate if our ray has a non-zero dot product with the normal
	if (ndotd != 0) {

		Vec3Df D = normal.projectOntoVector(vertices[tr.v[0]].p, normal);
		float odotn = normal.dotProduct(normal, rayOrigin);

		float t = (D.getLength() - odotn) / ndotd;

		// now we have t, check if we are inside the triangle
		Vec3Df p = rayOrigin + t * rayDest;
		// but p is also ... = a * v0 + b * v1 + (1-a-b) * v2

		float a, b, c;
		computeBarycentric(p, v0.p, v1.p, v2.p, a, b, c);
		
		if (!(a < 0 || a > 1 || b < 0 || a + b > 1)) {

			if (lastDistance > t) {
				hit = true;
				intersectionDepth = t;				
				hitPoint = p;
			}
		}
	}
	return hit;
}

KDNode* root;

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
int maxLevel = 1;

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

	std::vector<int> tris;
	for (int i = 0; i < MyMesh.triangles.size(); i++) {
		tris.push_back(i);
	}

	root = root->build(tris, 0);
	//one first move: initialize the first light source
	//at least ONE light source has to be in the scene!!!
	//here, we set it to the current location of the camera
	MyLightPositions.push_back(MyCameraPosition);
}

//return the color of your pixel.
Vec3Df performRayTracing(const Vec3Df & origin, const Vec3Df & dest)
{
	return trace(0, origin, dest);
}

Vec3Df trace(int level, const Vec3Df& origin, const Vec3Df& dest) {
	int triangleIndex = -1;
	Vec3Df hit = Vec3Df(0, 0, 0);
	intersect(origin, dest, triangleIndex, hit);
	if (triangleIndex != -1) {
		// we have a hit
		return shade(level, triangleIndex, hit);
	}
	return Vec3Df(0, 0, 0);
}

/**
 Checks if there is an intersection between the ray and the triangles, returns closest triangleIndex
 **/
void intersect(const Vec3Df & origin, const Vec3Df & dest, int triangleIndex, Vec3Df & hitPoint) {
	//Vec3Df rayDirection = dest - origin;
	//int triangleIn = -1;
	hit(root, origin, dest, hitPoint, triangleIndex);
	//triangleIndex = triangleIn;
}

/**
 Calculates shading color
 **/
Vec3Df shade(unsigned int level, const unsigned int triangleIndex, Vec3Df & hit) {
    
    if (level > maxLevel) {
        return Vec3Df(0, 0, 0);
    }
	Triangle tr = MyMesh.triangles[triangleIndex];
    
    Vec3Df directLight = Vec3Df(0, 0, 0);
    
    float lightintensity_ambient = 1.1;
    float lightintensity_specular = 0.8;
    
    // for each light
    for(std::vector<int>::size_type i = 0; i != MyLightPositions.size(); i++) {
		Vec3Df lightsource = MyLightPositions[i];

        Vec3Df direction = hit - lightsource;
        
        int closestTriangleIndex = -1;
        Vec3Df closestHit = Vec3Df(0, 0, 0);
        intersect(lightsource, direction, closestTriangleIndex, closestHit);
        
        Material m = getTriangleMaterial(triangleIndex);
        
        // calculate ambient term
        Vec3Df ambient = 0.5 * m.Kd();
        
        if (triangleIndex == closestTriangleIndex) {
            // let there be light
        
            
            // calculate diffuse term
            
            Triangle triangle = MyMesh.triangles[triangleIndex];
            Vertex v0 = MyMesh.vertices[triangle.v[0]];
            Vertex v1 = MyMesh.vertices[triangle.v[1]];
            Vertex v2 = MyMesh.vertices[triangle.v[2]];
            
            Vec3Df surfaceNormal = -1 * surfaceNormalTriangle(v0, v1, v2);
            
            direction.normalize();
            
            float costheta = surfaceNormal.dotProduct(surfaceNormal, direction);
            
            // diffuse
            Vec3Df diffuse = lightintensity_ambient * fabs(powf(costheta, 1)) * m.Kd();
        
            //std::cout << "Cos theta: " << surfaceNormal << " and " << diffuse << std::endl;

            // specular
            float n = 8;
            Vec3Df link = ((cameraOrigin - hit) - direction);
            link.normalize();
            Vec3Df specular = lightintensity_specular * powf(fabsf(link.dotProduct(link, surfaceNormal)), n) * m.Ks();
            
            //std::cout << "Pwo: " << powf(link.dotProduct(link, surfaceNormal), n) << std::endl;
            
            /*// compute reflected ray
            Vec3Df reflectedColor;
            Vec3Df reflectedRay = hit - ( (2 * hit.dotProduct(surfaceNormal, hit) ) * surfaceNormal);
            int reflectedTriangleIndex;
            Vec3Df reflectedHit;
            intersect(hit, reflectedRay, reflectedTriangleIndex, reflectedHit);
            if (reflectedTriangleIndex != -1) {
                // we have a hit
                reflectedColor = shade(level+1, reflectedTriangleIndex, reflectedHit);
            }*/
            
            directLight += ambient + diffuse + specular /*+ reflectedColor*/;
        } else {
            // shadow
            //directLight = Vec3Df(1, 1, 0);
            directLight += ambient; // just show ambient
        }

    }
    
    //return getTriangleColor(triangleIndex);

    return directLight;
}

Vec3Df getTriangleColor(const unsigned int triangleIndex) {
    Material m = MyMesh.materials[MyMesh.triangleMaterials[triangleIndex]];
    
    // for now return diffuse value
    return m.Kd();
}

Material getTriangleMaterial(const unsigned int triangleIndex) {
	if (triangleIndex >300) { std::cout << "index was: " << triangleIndex << std::endl; }
		
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
	testRayOrigin=rayOrigin;	
	testRayDestination=rayDestination;
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
