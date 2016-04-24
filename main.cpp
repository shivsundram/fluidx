#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <GLUT/glut.h>
#include <OpenGL/gl.h>
#include <math.h>       /* fmod */
#include <iostream>
#define DIM 100
#define DS (DIM*DIM)

using namespace std;

struct float2 
{
   float x;
   float y;
};


//UI and intialization is taken straight from stable fluids
//ego
// the opengl functionality, initiliaze methods
// and UI/click system have been taken from the nvidia stable fluids 
//However, the rest of the code, specifically the computaitonal and algorithmic 
//portions eg the fluid solvers parts are my own


// Particle data
static GLuint vbo = 0;                 // OpenGL vertex buffer object
static float2 *particles = NULL; // particle positions in host memory

//window data
static int wWidth  = DIM;
static int wHeight = DIM;

//click data
static int clicked  = 0;
static int lastx = 0, lasty = 0;


//velocity data
float2 *vfield = NULL;
float2 *vfield_temp = NULL;
float2 *rhs = NULL; 
float2 *divfield=NULL;
static float dt=.01f;

//random float between -.5f and .5f
float randhalf()
{
	return .5f-((float) rand())/RAND_MAX;
}

//initialize DIM*DIM particles, and small
//random offset as defined by randhalf
void initParticles(float2 *p, int dx, int dy)
{
    int i, j;

    for (i = 0; i < dy; i++)
    {
        for (j = 0; j < dx; j++)
        {
        	//printf("%f %f \n",0.0+j,0.0+i );
            p[j*dx+i].x = randhalf()/dx + (0.0f+i)/dx;
            p[j*dx+i].y = randhalf()/dy + (0.0f+j)/dy;
        }
    }
}

//init a constant velocity field
void initvfield(float2* field, float initial_x, float initial_y )
{
    for (int i = 0; i < DIM; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
        	//printf("%f %f \n",0.0+j,0.0+i );
            field[j*DIM+i].x = initial_x;
            field[j*DIM+i].y = initial_y;
        }
    }
}


//advect particles using velocity field and DT
//uses periodic boundary conditions
void advectParticles()
{
	for (int i = 0; i < DIM; i++)
	{
	    for (int j = 0; j < DIM; j++)
	    {
	    	//printf("%f %f \n",0.0+j,0.0+i );
	    	int x = particles[j*DIM+i].x*DIM;
	    	int y = particles[j*DIM+i].y*DIM;
	        particles[j*DIM+i].x =fmod(particles[j*DIM+i].x + dt*vfield[y*DIM+x].x,1.0f);
	        particles[j*DIM+i].y =fmod(particles[j*DIM+i].y + dt*vfield[y*DIM+x].y,1.0f);
	    }
	}

}


void printGridX(float2* u, int n )
{
    for (int y = 0 ; y<n; y++){
        for (int x =0 ; x<n; x++){
            cout<<u[y*n+x].x<<", ";
        }
        cout<<endl;
    }
}



float2 bilinterp(float2 * u, float x, float y, int dim)
{
    float2 result;
    //grab points: upper left, upper right, bottom left, bottom right

    //convert from gl space to 
    x=(int)(x*dim);
    y=(int)(y*dim);

    float ul = u[(int) ((floor(y)+1) *dim + floor(x))].x;
    float ur = u[(int) ((floor(y)+1) *dim + floor(x)+1)].x;
    float bl = u[(int) ((floor(y))   *dim + floor(x))].x;
    float br = u[(int) ((floor(y))   *dim + floor(x)+1)].x;

    //printf("%f %f\n",x,y );
    //printf("%f %f %f %f \n",ul, ur, bl, br);
    //interpolate in horizontal (x) direction 
    float fu = (x-floor(x))*ur + (floor(x)+1-x)*ul ; 
    float fb = (x-floor(x))*br + (floor(x)+1-x)*bl ; 

    //interpolate in vertical (j) direction 
    result.x = .99*(y-floor(y))*fu + .99*(floor(y)+1-y)*fb;   

    ul = u[(int) ((floor(y)+1) *dim + floor(x))].y;
    ur = u[(int) ((floor(y)+1) *dim + floor(x)+1)].y;
    bl = u[(int) ((floor(y))   *dim + floor(x))].y;
    br = u[(int) ((floor(y))   *dim + floor(x)+1)].y;

    //printf("%f %f %f %f \n",ul, ur, bl, br);
    //interpolate in horizontal (x) direction 
    fu = (x-floor(x))*ur + (floor(x)+1-x)*ul ; 
    fb = (x-floor(x))*br + (floor(x)+1-x)*bl ; 

    //interpolate in vertical (j) direction 
    result.y = .99*(y-floor(y))*fu + .99*(floor(y)+1-y)*fb;   
    //printf("%f \n",result.x );

    return result;
}

void vswap(){
    float2* temp = vfield;
    vfield=vfield_temp;
    vfield_temp =temp;
}

void advectVelocity(float time_step)
{
    float x_back = 0.0f; 
    float y_back = 0.0f; 
    for (int x = 1; x < DIM; x++){
        for (int y = 1 ;y < DIM; y++){
            x_back = ((1.0f*x)/DIM) - time_step*vfield[y*DIM+x].x;
            //printf("%f %f\n",x_back, fmod(x_back,1.0) );
            x_back = fmod(x_back, 1.0);
            y_back = ((1.0f*y)/DIM) - time_step*vfield[y*DIM+x].y;
            //x_back = fmod(y_back, 1.0);
            vfield_temp[y*DIM+x] =  bilinterp(vfield, x_back, y_back, DIM); 
            //printf("%i %i\n",x,y );
            //printf("%f\n",vfield[y*DIM+x].x  );
           //printf("%f %f %f \n",x_back, 1.0f*x/DIM, bilinterp(vfield, x_back, y_back, DIM).x);
        }
    }
    swap(vfield, vfield_temp); 
    return ;
}


void exchangeboundary(float2* u, int n){
    //exchange top and bottom
    int top = n-1;
    for (int i = 0; i<n; i++){
        u[i].x=u[(n)*(top-1)+i].x;
        u[i].y=u[(n)*(top-1)+i].y;

        u[n*(top)+i].x=u[n+i].x;
        u[n*(top)+i].y=u[n+i].y;

    }

    //exchange left and right
    for (int i = 1; i<n-1; i++){
        u[i*n].x=u[(n)*i+top-1].x;
        u[i*n].y=u[(n)*i+top-1].y;

        u[n*i+top].x=u[n*i+1].x;
        u[n*i+top].y=u[n*i+1].y;
    }
}

//gauss siedel preconditioned by the L1 inverse
//using L1 inverse requires only 1 array, eliminating
//the need for buffer swapping
void gsl1(float2* u, float2* b, float alpha, float rbeta,int n, int iters)
{
    for (int i = 0; i <iters; i ++)
    {
        for (int x = 1 ; x<n-1; x++)
        {
            for (int y =1; y<n-1; y++)
            {
                float up    = u[(y+1)*n+x].x;
                float down  = u[(y-1)*n+x].x;
                float left  = u[(y)*n+x-1].x;
                float right = u[(y)*n+x+1].x;
                u[y*n+x].x=rbeta*(up+down+left+right+alpha*b[y*n+x].x);
            
                up    = u[(y+1)*n+x].y;
                down  = u[(y-1)*n+x].y;
                left  = u[(y)*n+x-1].y;
                right = u[(y)*n+x+1].y;
                u[y*n+x].y=rbeta*(up+down+left+right+alpha*b[y*n+x].y);
            }
        }
        exchangeboundary(u, n);
    }
}


void diffuse()
{
    //todo, correctly derive jacobi weights
    int n = DIM;
    float dx2 = (1.0f/DIM)*(1.0f/DIM);
    float alpha= dx2/(dt*.0001f);
    float rbeta = 1.0f/(4.0f+alpha);
    gsl1(vfield, vfield, alpha, rbeta, DIM, 40);  
}

void divergence(float2 *u, float2 *d, int n){

    float scale = ((1.0)/n)*.5f;
    for (int x =1; x<n-1; x++){
        for (int y=1; y<n-1;y++){
            float2 up    = u[(y+1)*n+x];
            float2 down  = u[(y-1)*n+x];
            float2 left  = u[(y)*n+x-1];
            float2 right = u[(y)*n+x+1];

            d[y*n+x].x=((right.x-left.x)+(up.x-down.x))*scale;
            d[y*n+x].y=((right.y-left.y)+(up.y-down.y))*scale;


        }
    }

}



void subtractpressure(){


}

void project()
{


}


void simulate()
{
    advectVelocity(dt);
    advectParticles();
    diffuse();
}

void memcopy()
{
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float2) * DS,
                particles, GL_DYNAMIC_DRAW);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
}

void display(void)
{

	simulate();
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0, 1, 0, 1, 0, 1);

    // render points
    glClear(GL_COLOR_BUFFER_BIT);
    glClearColor(0, 0, 0, 1.0f);
    glColor4f(0,1,0,1.0f);
    glPointSize(1);
    glEnable(GL_POINT_SMOOTH);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnableClientState(GL_VERTEX_ARRAY);
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_CULL_FACE);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glVertexPointer(2, GL_FLOAT, 0, NULL);
    glDrawArrays(GL_POINTS, 0, DS);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glDisableClientState(GL_VERTEX_ARRAY);
    glDisableClientState(GL_TEXTURE_COORD_ARRAY);
    glDisable(GL_TEXTURE_2D);

    // Finish timing before swap buffers to avoid refresh sync
	memcopy();
	glutSwapBuffers();
	glutPostRedisplay();
}

void reshape(int x, int y)
{
    wWidth = x;
    wHeight = y;
    glViewport(0, 0, x, y);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0, 1, 0, 1, 0, 1);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glutPostRedisplay();
}


void click(int button, int updown, int x, int y)
{
    lastx = x;
    lasty = y;
    clicked = !clicked;
}

void motion(int x, int y)
{
    // Convert motion coordinates to domain
    float fx = (x / (float)wWidth);
    float fy = ((wHeight-y) / (float)wHeight);
    int nx = (int)(fx * DIM);
    int ny = (int)(fy * DIM);

    if (clicked)
    {
    		//printf("%f %f %i %i \n",fx,fy,nx, ny );
    		for (int i = -10 ; i<10 ;i ++){
    			for (int j = -10 ; j<10 ;j ++){
                    //TODO, use modulo function for indices of vfield
                    if(ny+i<DIM && ny+i>=0 && nx+j<DIM && nx+j>=0){
		              vfield[(ny+i)*DIM+nx+j].x +=.8f;
		              vfield[(ny+i)*DIM+nx+j].y +=.8f;
                    }
    			}
    		}

            //printf("%f\n",vfield[ny*DIM+nx].x );
 
    }

    glutPostRedisplay();
}



int initGL(int *argc, char **argv)
{
    glutInit(argc, argv);
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
    glutInitWindowSize(wWidth*2, wHeight*2);
    glutCreateWindow("Compute Stable Fluids");
    glutDisplayFunc(display);
    glutMouseFunc(click);
    glutMotionFunc(motion);
    glutReshapeFunc(reshape);
    return true;
}

void fluidx(int argc, char **argv){
    vfield = (float2 *)malloc(sizeof(float2) * DS);
    vfield_temp = (float2 *)malloc(sizeof(float2) * DS);

    initvfield(vfield, 0.01f, 0.01f);
    initvfield(vfield_temp , 0.01f, 0.01f);
   // initvfield(divfield, 0.0f, 0.0f);

    particles = (float2 *)malloc(sizeof(float2) * DS);
    initParticles(particles, DIM, DIM);
    initGL(&argc,  argv);


    glGenBuffers(1, &vbo);
    memcopy();

    glutMainLoop();
}


void testBoundary(){
    int n =4;
    float2* z = (float2*) malloc(n*n*sizeof(float2));
    for (int i =0 ; i<n*n; i++){
        z[i].x=0.0f;
        z[i].y=0.0f;
    }   
    printGridX(z, n);
    float c = 1.0f;
    cout<<endl;
    for (int i =1 ; i <n-1; i++){
        for (int j=1;j<n-1;j++){
            z[n*j+i].x=c;
            c=c+1.0f;
        }
    }
    printGridX(z, n);
    exchangeboundary(z,n);
    cout<<endl;
    printGridX(z, n);

}

void testswap(){
    vfield = (float2*) malloc(4*sizeof(float2));
    vfield_temp = (float2*) malloc(4*sizeof(float2));
    for (int i =0 ; i<4; i++){
        vfield[i].x=1.0f*i;
        vfield[i].y=1.0f*i;
    }   
    for (int i =0 ; i<4; i++){
        vfield_temp[i].x=0.0f*i;
        vfield_temp[i].y=0.0f*i;
    }
    printf("%f %f\n",vfield[1].x, vfield_temp[1].x );
    vswap();      
    printf("%f %f\n",vfield[1].x, vfield_temp[1].x );
    vswap();
    printf("%f %f\n",vfield[1].x, vfield_temp[1].x );
}

void testbilinterp(){
    float x = .25f; 
    float y = .5f;
    float2 * u = (float2*) malloc(4*sizeof(float2));
    for (int i =0 ; i<4; i++){
        u[i].x=1.0f*i;
        u[i].y=0.0f;
    }
    printf("u(%f %f).x = %f\n",x,y, bilinterp(u,x,y,2).x);
    printf("u(%f %f).y = %f\n",x,y, bilinterp(u,x,y,2).y);
    free(u);
}

void testgsl1(){
    //formula
    //u[y*n+x].y=rbeta*(up+down+left+right+alpha*b[y*n+x].y);
    int n =4;
    float2* z = (float2*) malloc(n*n*sizeof(float2));
    for (int i =0 ; i<n*n; i++){
        z[i].x=0.0f;
        z[i].y=0.0f;
    }   
    float c = 1.0f;
    cout<<endl;
    for (int i =1 ; i <n-1; i++){
        for (int j=1;j<n-1;j++){
            z[n*j+i].x=c;
            c=c+1.0f;
        }
    }
    printGridX(z, n);
    gsl1(z, z, 1.0f, 1.0f, n,1);
    cout<<endl;
    printGridX(z,n);
}

void testdivergence(){


}

int main(int argc, char **argv)
{
    fluidx(argc, argv);
    //testbilinterp();
    //testswap();
    //testjacobi();
    //testBoundary();
    //testgsl1();
	return 0;
}