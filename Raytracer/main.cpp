#ifdef WIN32
#include <windows.h>
#include <GL/glut.h>
#endif
#ifdef __APPLE__
#include <GLUT/glut.h>
#include <pthread.h>
#endif
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "raytracing.h"
#include "mesh.h"
#include "traqueboule.h"
#include "imageWriter.h"
#include "paths.h"
#include <ctime>
#include <chrono>
#include <thread>
#include <atomic>
#include <future>

//This is the main application
//Most of the code in here, does not need to be modified.
//It is enough to take a look at the function "drawFrame",
//in case you want to provide your own different drawing functions



Vec3Df MyCameraPosition;

//MyLightPositions stores all the light positions to use
//for the ray tracing. Please notice, the light that is 
//used for the real-time rendering is NOT one of these, 
//but following the camera instead.
std::vector<Vec3Df> MyLightPositions;

//Main mesh 
Mesh MyMesh;

unsigned int WindowSize_X = OUTPUT_PIXELS;  // resolution X
unsigned int WindowSize_Y = OUTPUT_PIXELS;  // resolution Y




/**
 * Main function, which is drawing an image (frame) on the screen
*/
void drawFrame( )
{
	yourDebugDraw();
}

//animation is called for every image on the screen once
void animate()
{
	MyCameraPosition=getCameraPosition();
	glutPostRedisplay();
}



void display(void);
void reshape(int w, int h);
void keyboard(unsigned char key, int x, int y);

/**
 * Main Programme
 */
int main(int argc, char** argv)
{
    glutInit(&argc, argv);

    //framebuffer setup
    glutInitDisplayMode( GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH );

    // positioning and size of window
    glutInitWindowPosition(200, 100);
    glutInitWindowSize(WindowSize_X,WindowSize_Y);
    glutCreateWindow(argv[0]);	

    //initialize viewpoint
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glTranslatef(0,0,-4);
    tbInitTransform();     // This is for the trackball, please ignore
    tbHelp();             // idem
	MyCameraPosition=getCameraPosition();

	//activate the light following the camera
    glEnable( GL_LIGHTING );
    glEnable( GL_LIGHT0 );
    glEnable(GL_COLOR_MATERIAL);
    int LightPos[4] = {0,0,2,0};
    int MatSpec [4] = {1,1,1,1};
    glLightiv(GL_LIGHT0,GL_POSITION,LightPos);

	//normals will be normalized in the graphics pipeline
	glEnable(GL_NORMALIZE);
    //clear color of the background is black.
	glClearColor (0.0, 0.0, 0.0, 0.0);

	
	// Activate rendering modes
    //activate depth test
	glEnable( GL_DEPTH_TEST ); 
    //draw front-facing triangles filled
	//and back-facing triangles as wires
    glPolygonMode(GL_FRONT,GL_FILL);
    glPolygonMode(GL_BACK,GL_LINE);
    //interpolate vertex colors over the triangles
	glShadeModel(GL_SMOOTH);


	// glut setup... to ignore
    glutReshapeFunc(reshape);
    glutKeyboardFunc(keyboard);
    glutDisplayFunc(display);
    glutMouseFunc(tbMouseFunc);    // trackball
    glutMotionFunc(tbMotionFunc);  // uses mouse
    glutIdleFunc( animate);


	init();

    
	//main loop for glut... this just runs your application
    glutMainLoop();
        
    return 0;  // execution never reaches this point
}











/**
 * OpenGL setup - functions do not need to be changed! 
 * you can SKIP AHEAD TO THE KEYBOARD FUNCTION
 */
//what to do before drawing an image
 void display(void)
{
	glPushAttrib(GL_ALL_ATTRIB_BITS);//store GL state
    // Effacer tout
    glClear( GL_COLOR_BUFFER_BIT  | GL_DEPTH_BUFFER_BIT); // clear image
    
    glLoadIdentity();  

    tbVisuTransform(); // init trackball

    drawFrame( );    //actually draw

    glutSwapBuffers();//glut internal switch
	glPopAttrib();//return to old GL state
}
//Window changes size
void reshape(int w, int h)
{
    glViewport(0, 0, (GLsizei) w, (GLsizei) h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    //glOrtho (-1.1, 1.1, -1.1,1.1, -1000.0, 1000.0);
    gluPerspective (50, (float)w/h, 0.01, 10);
    glMatrixMode(GL_MODELVIEW);
}


//transform the x, y position on the screen into the corresponding 3D world position
void produceRay(int x_I, int y_I, Vec3Df * origin, Vec3Df * dest)
{
		int viewport[4];
		double modelview[16];
		double projection[16];
		glGetDoublev(GL_MODELVIEW_MATRIX, modelview); //recuperer matrices
		glGetDoublev(GL_PROJECTION_MATRIX, projection); //recuperer matrices
		glGetIntegerv(GL_VIEWPORT, viewport);//viewport
		int y_new = viewport[3] - y_I;

		double x, y, z;
		
		gluUnProject(x_I, y_new, 0, modelview, projection, viewport, &x, &y, &z);
		origin->p[0]=float(x);
		origin->p[1]=float(y);
		origin->p[2]=float(z);
		gluUnProject(x_I, y_new, 1, modelview, projection, viewport, &x, &y, &z);
		dest->p[0]=float(x);
		dest->p[1]=float(y);
		dest->p[2]=float(z);
}

// Raytraces a section from yStart to yEnd
void rayTraceSection(
	int yStart,
	int yEnd,
	Image* result, 
	Vec3Df& origin, Vec3Df& dest,
	Vec3Df& origin00, Vec3Df& dest00,
	Vec3Df& origin01, Vec3Df& dest01,
	Vec3Df& origin10, Vec3Df& dest10,
	Vec3Df& origin11, Vec3Df& dest11,
	int threadID
	)
{
	cout << "Thread " << threadID << " starting... \n" <<
		"from " << yStart << " to " << yEnd << endl;
	for (double y = yStart; y < yEnd; ++y) {
		//std::cout << WindowSize_Y << std::endl;
		for (unsigned int x = 0; x < WindowSize_X; ++x)
		{
			//produce the rays for each pixel, by interpolating 
			//the four rays of the frustum corners.
			float xscale = 1.0f - float(x) / (WindowSize_X - 1);
			float yscale = 1.0f - float(y) / (WindowSize_Y - 1);

			origin = yscale*(xscale*origin00 + (1 - xscale)*origin10) +
				(1 - yscale)*(xscale*origin01 + (1 - xscale)*origin11);
			dest = yscale*(xscale*dest00 + (1 - xscale)*dest10) +
				(1 - yscale)*(xscale*dest01 + (1 - xscale)*dest11);

			//launch raytracing for the given ray.
			Vec3Df rgb = performRayTracing(origin, dest);
			//store the result in an image 
			(*result).setPixel(x, y, RGBValue(rgb[0], rgb[1], rgb[2]));
		}
	}

	cout << "Thread " << threadID << " done" << endl;
}








// react to keyboard input
void keyboard(unsigned char key, int x, int y)
{
    printf("key %d pressed at %d,%d\n",key,x,y);
    fflush(stdout);
    switch (key)
    {
	//add/update a light based on the camera position.
	case 'L':
		MyLightPositions.push_back(getCameraPosition());
		break;
	case 'l':
		MyLightPositions[MyLightPositions.size()-1]=getCameraPosition();
		break;
	case 'r':
	{
		//Pressing r will launch the raytracing.
		std::cout<<"Raytracing"<<endl;
		auto begin = std::chrono::high_resolution_clock::now();
		//time_t startTime = time(0);

		//Setup an image with the size of the current image.
		Image result(WindowSize_X,WindowSize_Y);
		
		//produce the rays for each pixel, by first computing
		//the rays for the corners of the frustum.
		Vec3Df origin00, dest00;
		Vec3Df origin01, dest01;
		Vec3Df origin10, dest10;
		Vec3Df origin11, dest11;
		Vec3Df origin, dest;


		produceRay(0,0, &origin00, &dest00);
		produceRay(0,WindowSize_Y-1, &origin01, &dest01);
		produceRay(WindowSize_X-1,0, &origin10, &dest10);
		produceRay(WindowSize_X-1,WindowSize_Y-1, &origin11, &dest11);

		// True for multithreaded tracing
#ifdef __APPLE__
		bool MULTITHREAD = false;
#else
		boolean MULTITHREAD = false;
#endif

		if (MULTITHREAD){
			unsigned static const int nthreads = std::thread::hardware_concurrency();
			cout << "Amount of Cores: " << nthreads << endl;
			int rowsPerThread = WindowSize_Y / nthreads;
			cout << "Amount of rows per thread: " << rowsPerThread << endl;
			int leftOverRows = WindowSize_Y % (rowsPerThread * nthreads);
			cout << "Amount of leftover rows: " << leftOverRows << endl;

			std::thread *t = new std::thread[nthreads];

			// test stuff which is not dependant on the amount of cores the machine has, but is there to see if the multithreading works in the first place
			/*std::thread t1(rayTraceSection, 0, WindowSize_Y / 3, &result, origin, dest, origin00, dest00, origin01, dest01, origin10, dest10, origin11, dest11, 1);
			std::thread t2(rayTraceSection, (WindowSize_Y / 3), (WindowSize_Y / 3) * 2, &result, origin, dest, origin00, dest00, origin01, dest01, origin10, dest10, origin11, dest11, 2);
			std::thread t3(rayTraceSection, (WindowSize_Y / 3) * 2, WindowSize_Y, &result, origin, dest, origin00, dest00, origin01, dest01, origin10, dest10, origin11, dest11, 3);

			t1.join();
			t2.join();
			t3.join();*/

			// this is the row that the thread starts on
			int currStart = 0;
			int currEnd = rowsPerThread;

			//Launch a group of threads
			for (int i = 0; i < nthreads; ++i) {
				// start the thread
				//t[i] = std::thread(rayTraceSection, currStart, currEnd, &result, origin, dest, origin00, dest00, origin01, dest01, origin10, dest10, origin11, dest11, i);
				// set the start for the next loop
				currStart += rowsPerThread;

				// set the end for the next loop, if the next loop is the last loop tack on the remaining lines
				if (i == nthreads - 2)
				{
					currEnd = WindowSize_Y;
				}
				else
				{
					currEnd += rowsPerThread;
				}
			}

			for (int i = 0; i < nthreads; ++i)
			{
				t[i].join();
			}
		}
		else
		{
			cout << "Single thread Raytracing" << endl;
			for (double y = 0; y < WindowSize_Y; ++y) {
				double perc = round((y / WindowSize_Y) * 100);
				std::cout << "[Raytracing with " << WindowSize_Y << " pixels running]  [" << perc << "%]" << '\r' << std::flush;
				//std::cout << WindowSize_Y << std::endl;
				for (unsigned int x = 0; x < WindowSize_X; ++x)
				{
					//produce the rays for each pixel, by interpolating 
					//the four rays of the frustum corners.
					float xscale = 1.0f - float(x) / (WindowSize_X - 1);
					float yscale = 1.0f - float(y) / (WindowSize_Y - 1);

					origin = yscale*(xscale*origin00 + (1 - xscale)*origin10) +
						(1 - yscale)*(xscale*origin01 + (1 - xscale)*origin11);
					dest = yscale*(xscale*dest00 + (1 - xscale)*dest10) +
						(1 - yscale)*(xscale*dest01 + (1 - xscale)*dest11);

					//launch raytracing for the given ray.
					Vec3Df rgb = performRayTracing(origin, dest);
					//store the result in an image 
					result.setPixel(x, y, RGBValue(rgb[0], rgb[1], rgb[2]));
				}
			}
		}

		std::cout << "\n" << endl;

		result.writeImage(IMAGE_PATH);
		//time_t endTime = time(0);
		auto end = std::chrono::high_resolution_clock::now();
		auto dur = end - begin;
		auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(dur).count();
		std::cout << "The tracing took " << float(ms/1000.0) << " seconds." << std::endl;
		break;
	}
	case 27:     // touche ESC
        exit(0);
    }

	
	//produce the ray for the current mouse position
	Vec3Df testRayOrigin, testRayDestination;
	produceRay(x, y, &testRayOrigin, &testRayDestination);

	yourKeyboardFunc(key,x,y, testRayOrigin, testRayDestination);
}
