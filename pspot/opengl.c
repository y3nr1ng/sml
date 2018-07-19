/*--------------------------------------------------------------------------
 *
 *	Use OpenGL to display the image on the screen
 *
 *-------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#ifdef __APPLE__
#  include <GLUT/glut.h>
#else
#  include <GL/glut.h>
#endif
#include "glfont/glfont.h"
#include "pspot.h"

static para_t *p;

typedef struct {
    int    press;
    double x1, y1, x2, y2;
} mouse_t;

static GLuint         texName;
static GLFONT         font;
static unsigned char *Timg;
static mouse_t        ms;
static rec_t         *markrec;

/*--------------------------------------------------------------------------
 *
 *	Initialization for different displaying modes.
 *
 *-------------------------------------------------------------------------*/

void texture_init(void)
{
    int  i, np;

    np = p->imgWidth * p->imgHeight;
    if ((Timg = malloc(np*3)) == NULL)
	pstop("!!! texture_init: not enough memory.\n");
    for (i=0; i < np; i++) {
	Timg[i*3  ] = p->gimg[i];
	Timg[i*3+1] = p->gimg[i];
	Timg[i*3+2] = p->gimg[i];
    }

    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glGenTextures(1, &texName);
    glBindTexture(GL_TEXTURE_2D, texName);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, p->imgWidth, p->imgHeight,
		 0, GL_RGB, GL_UNSIGNED_BYTE, Timg);
}

/*--------------------------------------------------------------------------
 *
 *	Output the framebuffer to JPEG file.
 *
 *-------------------------------------------------------------------------*/

static void output_framebuffer(para_t *p)
{
    int            i, j, x, y, iW, iH;
    unsigned char *img;

    iW = p->imgWidth;
    iH = p->imgHeight;
    if ((img = malloc(3*iW*iH)) == NULL)
	pstop("!!! output_framebuffer: not enough memory.\n");
    glReadPixels(0, 0, iW, iH, GL_RGB, GL_UNSIGNED_BYTE, img);
    for (j=0; j < iH; j++) {
    for (i=0; i < iW;  i++) {
	x = i + j*iW;
	y = i + (iH-j-1)*iW;
	p->gimg[x*3  ] = img[y*3  ];
	p->gimg[x*3+1] = img[y*3+1];
	p->gimg[x*3+2] = img[y*3+2];
    }}
    output_JPEG(p, 0, p->inpfn);
    free(img);
}

/*--------------------------------------------------------------------------
 *
 *	Plot the image spots and rectangles.
 *
 *-------------------------------------------------------------------------*/

static void glfontInit(void)
{
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    if (! glFontCreate(&font, p->fontfn, 0)) {
	printf("!!! glfontInit: cannot create font: %s\n", p->fontfn);
	return;
    }
}

static void DisplayDraw(void)
{
    int    i;
    char   rid[16];
    float  fontscale, fontshift;

    glClearColor(1.0, 1.0, 1.0, 1.0);   // white background
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glColor3f( 0, 0, 1);
    glBegin(GL_POINTS);
    for (i=0; i < p->np; i++)
        glVertex3d(p->x[i], p->y[i], 0.0);
    glEnd();
    if (markrec != NULL) {
	int  idx;
        glColor3f( 0, 0.4, 0);
        glBegin(GL_POINTS);
	for (i=0; i < markrec->n_sp; i++) {
	    idx = markrec->sp[i];
	    glVertex3d(p->x[idx], p->y[idx], 0.0);
	}
	glEnd();
    }
    if (ms.press == 1) {
	glColor3f( 1.0, 0.7, 0);
	glBegin(GL_LINE_STRIP);
	glVertex3d(ms.x1, ms.y1, 0.0);
	glVertex3d(ms.x1, ms.y2, 0.0);
	glVertex3d(ms.x2, ms.y2, 0.0);
	glVertex3d(ms.x2, ms.y1, 0.0);
	glVertex3d(ms.x1, ms.y1, 0.0);
	glEnd();
    }
    glColor3f( 1.0, 0, 0);
    for (i=0; i < p->n_rec; i++) {
        glBegin(GL_LINE_STRIP);
	glVertex3d(p->rec[i].x1, p->rec[i].y1, 0.0);
	glVertex3d(p->rec[i].x1, p->rec[i].y2, 0.0);
	glVertex3d(p->rec[i].x2, p->rec[i].y2, 0.0);
	glVertex3d(p->rec[i].x2, p->rec[i].y1, 0.0);
	glVertex3d(p->rec[i].x1, p->rec[i].y1, 0.0);
        glEnd();
    }
    if (font.TexWidth > 0 && font.TexHeight > 0) {
	fontscale = p->iW / ((double)(p->imgHeight)) * 9.0;
	fontshift = p->iW / ((double)(p->imgHeight)) * 2.0;
	glEnable(GL_TEXTURE_2D);
	glFontBegin(&font);
	for (i=0; i < p->n_rec; i++) {
	    sprintf(rid, "%d", i);
	    glLoadIdentity();
	    glTranslatef(p->rec[i].x1+fontshift, p->rec[i].y2-fontshift, 0);
	    glScalef(fontscale, fontscale, 1.0);
	    glFontTextOut(rid, 0.0, 0.0, 0.0);
	}
	glFontEnd();
	glDisable(GL_TEXTURE_2D);
    }
    glFlush();
}

void DisplayTexture(void)
{
   glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
   glEnable(GL_TEXTURE_2D);
   glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
   glBindTexture(GL_TEXTURE_2D, texName);
   glBegin(GL_QUADS);
   glTexCoord2d(0.0, 1.0); glVertex3d(0.0,   0.0,   0.0);
   glTexCoord2d(0.0, 0.0); glVertex3d(0.0,   p->imgHeight, 0.0);
   glTexCoord2d(1.0, 0.0); glVertex3d(p->imgWidth, p->imgHeight, 0.0);
   glTexCoord2d(1.0, 1.0); glVertex3d(p->imgWidth, 0.0,   0.0);
   glEnd();
   glFlush();
   glDisable(GL_TEXTURE_2D);
}

static void Keyboard(unsigned char key, int x, int y)
{
    double  xx, yy;

    if (key == 0x71) {				// 'q'
	if (p->datafmt == 0) {
	    recspot(p);
	    xcor_output(p);
	    output_framebuffer(p);
	}
	exit(0);
    }
    else if (key == 0x1b) {			// ESC
	if (p->datafmt == 0) {
	    xx = ((double)x) / ((double)p->imgWidth) * p->iW;
	    yy =-((double)y) / ((double)p->imgHeight)* p->iH + p->iH;
	    del_rec(p, xx, yy);
	    glutPostRedisplay();
	}
    }
    else if (key == 0x76) {			// 'v'
	if (p->datafmt == 0) {
	    xx = ((double)x) / ((double)p->imgWidth) * p->iW;
	    yy =-((double)y) / ((double)p->imgHeight)* p->iH + p->iH;
	    markrec = search_rec(p, xx, yy);
	    glutPostRedisplay();
	}
    }
}

static void MouseButton(int button, int state, int x, int y)
{
    if (button != GLUT_LEFT_BUTTON) return;

    if (state == GLUT_DOWN) {
	ms.press = 1;
	ms.x1    = ((double)x) / ((double)p->imgWidth) * p->iW;
	ms.y1    =-((double)y) / ((double)p->imgHeight)* p->iH + p->iH;
	ms.x2    = ((double)x) / ((double)p->imgWidth) * p->iW;
	ms.y2    =-((double)y) / ((double)p->imgHeight)* p->iH + p->iH;
    }
    else {
	ms.press = 0;
	ms.x2    = ((double)x) / ((double)p->imgWidth) * p->iW;
	ms.y2    =-((double)y) / ((double)p->imgHeight)* p->iH + p->iH;
	push_rec(p, ms.x1, ms.x2, ms.y1, ms.y2);
	glutPostRedisplay();
    }
}

static void MouseMotion(int x, int y)
{
    if (ms.press == 1) {
	ms.x2 = ((double)x) / ((double)p->imgWidth) * p->iW;
	ms.y2 =-((double)y) / ((double)p->imgHeight)* p->iH + p->iH;
	glutPostRedisplay();
    }
}

static void WindowSize(int w, int h)
{
    glViewport(0, 0, w, h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    if (p->datafmt == 0)
	glOrtho(0.0, p->iW, 0.0, p->iH, -1, 1);
    else
	glOrtho(0.0, p->imgWidth, 0.0, p->imgHeight, -1, 1);
}

void display_image(int argc, char **argv, para_t *pp)
{
    p = pp;

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
    glutInitWindowSize(p->imgWidth, p->imgHeight);
    glutInitWindowPosition(0, 0);
    glutCreateWindow(p->inpfn);
    glfontInit();
   
    if (p->datafmt == 0) {
	xcor_input(p);
	glutDisplayFunc(DisplayDraw);
	glutMouseFunc(MouseButton);
	glutMotionFunc(MouseMotion);
    }
    else {
	texture_init();
	glutDisplayFunc(DisplayTexture);
    }
    glutReshapeFunc(WindowSize);
    glutKeyboardFunc(Keyboard);
    glutMainLoop();
}

