#include <stdlib.h>
#include <memory.h>
#include <math.h>
#include <assert.h>
#include <stdbool.h>
#include <stdio.h>

#include <GL/gl.h>
#include <GL/glut.h>

#define STB_IMAGE_IMPLEMENTATION
#include <stb_image.h>

#define min(a, b) ((a) < (b) ? (a) : (b))
#define max(a, b) ((b) < (a) ? (a) : (b))

double *get_diff_polys(const double *x, const double *y, int n)
{
    double *y2 = malloc(sizeof(double) * n);
    for (int i = 0; i < n - 1; i++)
    {
        y2[i] = (y[i] - y[i + 1]) / (x[i + 1] - x[i]);
    }
    for (int i = 0; i < n - 2; i++)
    {
        y2[i] = (y2[i] - y2[i + 1]) / (x[i + 2] - x[i]);
    }
    return y2;
}

typedef struct
{
    double *x;
    double *y;
    double *M;
    int n;
} nifs3_t;

nifs3_t *nifs3_init(const double *x, const double *y, int n)
{
    nifs3_t *interp = malloc(sizeof(nifs3_t));

    double *M = malloc(sizeof(double) * n);

    double *q = malloc(sizeof(double) * (n - 1));
    double *u = malloc(sizeof(double) * (n - 1));

    double *d = get_diff_polys(x, y, n);

    q[0] = u[0] = 0;

    for (int i = 1; i < n - 1; i++)
    {
        double h_i = x[i] - x[i - 1];
        double h_i1 = x[i + 1] - x[i];
        double lam = h_i / (h_i + h_i1);

        double p = lam * q[i - 1] + 2;
        q[i] = (lam - 1) / p;
        u[i] = (d[i - 1] - lam * u[i - 1]) / p;
    }

    M[0] = M[n] = 0;
    M[n - 1] = u[n - 1];

    for (int i = n - 2; i > 0; i--)
        M[i] = u[i] - q[i] * M[i + 1];

    free(u);
    free(q);
    free(d);

    interp->x = malloc(sizeof(double) * n);
    interp->y = malloc(sizeof(double) * n);
    interp->M = M;
    interp->n = n;

    memcpy(interp->x, x, sizeof(double) * n);
    memcpy(interp->y, y, sizeof(double) * n);

    return interp;
}

void nifs3_free(nifs3_t *interp)
{
    if (interp == NULL)
        return;
    free(interp->x);
    free(interp->y);
    free(interp->M);
    free(interp);
}

double nifs3_get(nifs3_t *interp, double x)
{
    if (x < interp->x[0] || x > interp->x[interp->n - 1])
        return 0;

    int i = 0;
    while (x > interp->x[i])
        i++;

    double h = interp->x[i] - interp->x[i - 1];
    double t1 = interp->x[i] - x;
    double t2 = x - interp->x[i - 1];

    return (1 / h) * (interp->M[i - 1] / 6 * t1 * t1 * t1 +
                      interp->M[i] / 6 * t2 * t2 * t2 +
                      (interp->y[i - 1] - interp->M[i - 1] / 6 * h * h) * t1 +
                      (interp->y[i] - interp->M[i] / 6 * h * h) * t2);
}

typedef struct
{
    nifs3_t *iX, *iY;

    double xMax, xMin, yMax, yMin;

} nifs3_2d_t;

#define MAX_INTERPOLATORS 128
nifs3_2d_t interp[MAX_INTERPOLATORS];

void draw_interpolator(int inp, double *us, int n)
{
    glBegin(GL_LINE_STRIP);
    for (int i = 0; i < n; i++)
    {
        double x = nifs3_get(interp[inp].iX, us[i]);
        double y = nifs3_get(interp[inp].iY, us[i]);
        glVertex2d(x, y);
    }
    glEnd();
}

void free_interpolator2d(int i)
{
    nifs3_free(interp[i].iX);
    nifs3_free(interp[i].iY);
    memset(&interp[i], 0, sizeof(nifs3_2d_t));
}

void free_interpolator2ds()
{
    for (int i = 0; i < MAX_INTERPOLATORS; i++)
        free_interpolator2d(i);
}

int create_interpolator2d(const double *x, const double *y, const double *t, int n)
{
    for (int i = 0; i < MAX_INTERPOLATORS; i++)
    {
        if (interp[i].iX == NULL)
        {
            interp[i].iX = nifs3_init(t, x, n);
            interp[i].iY = nifs3_init(t, y, n);

            interp[i].xMin = INFINITY;
            interp[i].xMax = -INFINITY;
            interp[i].yMin = INFINITY;
            interp[i].yMax = -INFINITY;

            for (int j = 0; j < n; j++)
            {
                interp[i].xMin = min(interp[i].xMin, x[j]);
                interp[i].xMax = max(interp[i].xMax, x[j]);
                interp[i].yMin = min(interp[i].yMin, y[j]);
                interp[i].yMax = max(interp[i].yMax, y[j]);
            }

            return i;
        }
    }
    return -1;
}

const double Xs[] = {5.5, 8.5, 10.5, 13, 17, 20.5, 24.5, 28, 32.5, 37.5, 40.5, 42.5, 45, 47,
                     49.5, 50.5, 51, 51.5, 52.5, 53, 52.8, 52, 51.5, 53, 54, 55, 56, 55.5, 54.5, 54, 55, 57, 58.5,
                     59, 61.5, 62.5, 63.5, 63, 61.5, 59, 55, 53.5, 52.5, 50.5, 49.5, 50, 51, 50.5, 49, 47.5, 46,
                     45.5, 45.5, 45.5, 46, 47.5, 47.5, 46, 43, 41, 41.5, 41.5, 41, 39.5, 37.5, 34.5, 31.5, 28, 24,
                     21, 18.5, 17.5, 16.5, 15, 13, 10, 8, 6, 6, 6, 5.5, 3.5, 1, 0, 0, 0.5, 1.5, 3.5, 5, 5, 4.5, 4.5, 5.5,
                     6.5, 6.5, 5.5};
const double Ys[] = {41, 40.5, 40, 40.5, 41.5, 41.5, 42, 42.5, 43.5, 45, 47, 49.5, 53, 57, 59,
                     59.5, 61.5, 63, 64, 64.5, 63, 61.5, 60.5, 61, 62, 63, 62.5, 61.5, 60.5, 60, 59.5, 59, 58.5,
                     57.5, 55.5, 54, 53, 51.5, 50, 50, 50.5, 51, 50.5, 47.5, 44, 40.5, 36, 30.5, 28, 25.5, 21.5,
                     18, 14.5, 10.5, 7.50, 4, 2.50, 1.50, 2, 3.50, 7, 12.5, 17.5, 22.5, 25, 25, 25, 25.5, 26.5,
                     27.5, 27.5, 26.5, 23.5, 21, 19, 17, 14.5, 11.5, 8, 4, 1, 0, 0.5, 3, 6.50, 10, 13, 16.5, 20.5,
                     25.5, 29, 33, 35, 36.5, 39, 41};

#define countof(x) (sizeof(x) / sizeof(x[0]))

#define N 1000
double Us[N];

void linspace(double start, double end, int n, double *data)
{
    double step = (end - start) / (n - 1);
    for (int i = 0; i < n; i++)
        data[i] = start + i * step;
}

struct scene_data_t
{
    double xCenter, yCenter;

    // units per pixel
    double scale;

    int w, h;

    bool dragging;
    int lastX, lastY;

    double xMin, xMax, yMin, yMax;

    bool showImage;
    int imageW, imageH;
    int texture;

} scene_data;

void init()
{
    assert(countof(Xs) == countof(Ys));

    const int n = countof(Xs);
    double Ts[n];
    linspace(0, 1, n, Ts);
    linspace(0, 1, N, Us);

    create_interpolator2d(Xs, Ys, Ts, n);

    // create window
    glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
    glutInitWindowSize(500, 500);
    glutCreateWindow("Interpolation");

    // set up coordinate system
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    double xMin = INFINITY;
    double xMax = -INFINITY;
    double yMin = INFINITY;
    double yMax = -INFINITY;

    for (int i = 0; i < MAX_INTERPOLATORS; i++)
    {
        if (interp[i].iX == NULL)
            continue;

        xMin = min(xMin, interp[i].xMin);
        xMax = max(xMax, interp[i].xMax);
        yMin = min(yMin, interp[i].yMin);
        yMax = max(yMax, interp[i].yMax);
    }

    scene_data.w = glutGet(GLUT_WINDOW_WIDTH);
    scene_data.h = glutGet(GLUT_WINDOW_HEIGHT);

    scene_data.xCenter = (xMax + xMin) / 2;
    scene_data.yCenter = (yMax + yMin) / 2;

    scene_data.scale = min((xMax - xMin) / scene_data.w,
                           (yMax - yMin) / scene_data.h);
    scene_data.scale *= 1.1;

    scene_data.showImage = false;

    stbi_set_flip_vertically_on_load(true);

    int components = 4;
    stbi_uc *data = stbi_load("image.png", &scene_data.imageW, &scene_data.imageH, &components, 4);
    if (data == NULL)
    {
        printf("Failed to load image\n");
        exit(1);
    }

    glEnable(GL_TEXTURE_2D);
    glGenTextures(1, &scene_data.texture);
    glBindTexture(GL_TEXTURE_2D, scene_data.texture);
    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_GENERATE_MIPMAP, GL_TRUE);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, scene_data.imageW, scene_data.imageH, 0,
                 components == 4 ? GL_RGBA : GL_RGB, GL_UNSIGNED_BYTE, data);
    stbi_image_free(data);
    glDisable(GL_TEXTURE_2D);
}

void mouse(int button, int state, int x, int y)
{
    if (button == GLUT_LEFT_BUTTON)
    {
        scene_data.dragging = (state == GLUT_DOWN);

        scene_data.lastX = x;
        scene_data.lastY = y;
    }
    else if (button == 3)
    {
        scene_data.scale *= 1.1;
    }
    else if (button == 4)
    {
        scene_data.scale /= 1.1;
    }

    glutPostRedisplay();
}

void keyboard(unsigned char c, int x, int y)
{
    if (c == 'i')
        scene_data.showImage = !scene_data.showImage;

    glutPostRedisplay();
}

void motion(int x, int y)
{
    if (scene_data.dragging)
    {
        double dx = x - scene_data.lastX;
        double dy = y - scene_data.lastY;
        scene_data.xCenter -= dx * scene_data.scale;
        scene_data.yCenter += dy * scene_data.scale;
    }

    scene_data.lastX = x;
    scene_data.lastY = y;

    glutPostRedisplay();
}

void reshape(int w, int h)
{
    scene_data.scale *= (double)scene_data.w / w;

    scene_data.w = w;
    scene_data.h = h;

    glViewport(0, 0, w, h);
    glutPostRedisplay();
}

void drawText()
{
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    glOrtho(0, scene_data.w, 0, scene_data.h, -1, 1);

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    glRasterPos2i(10, scene_data.h - 40); // move in 10 pixels from the left and bottom edges

    char str[1024];
    sprintf(str, "X: [%f %f], Y: [%f %f]",
            scene_data.xMin, scene_data.xMax,
            scene_data.yMin, scene_data.yMax);

    for (char *s = str; *s; s++)
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, *s);

    glPopMatrix();

    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
}

void display()
{
    glClearColor(0, 0, 0, 0);
    glClear(GL_COLOR_BUFFER_BIT);

    double xCenter = scene_data.xCenter;
    double yCenter = scene_data.yCenter;

    double w = scene_data.w;
    double h = scene_data.h;

    double scale = scene_data.scale;

    double xOffset = scale * w / 2;
    double yOffset = scale * h / 2;

    scene_data.xMin = xCenter - xOffset;
    scene_data.xMax = xCenter + xOffset;
    scene_data.yMin = yCenter - yOffset;
    scene_data.yMax = yCenter + yOffset;

    if (scene_data.showImage)
    {
        glEnable(GL_TEXTURE_2D);
        glBindTexture(GL_TEXTURE_2D, scene_data.texture);
        glBegin(GL_QUADS);
        glTexCoord2f(0, 0);
        glVertex2f(scene_data.imageW / -2, scene_data.imageH / -2);
        glTexCoord2f(1, 0);
        glVertex2f(scene_data.imageW / 2, scene_data.imageH / -2);
        glTexCoord2f(1, 1);
        glVertex2f(scene_data.imageW / 2, scene_data.imageH / 2);
        glTexCoord2f(0, 1);
        glVertex2f(scene_data.imageW / -2, scene_data.imageH / 2);
        glEnd();
        glDisable(GL_TEXTURE_2D);
    }

    drawText();

    glLoadIdentity();
    glOrtho(scene_data.xMin, scene_data.xMax,
            scene_data.yMin, scene_data.yMax,
            -1, 1);

    for (int i = 0; i < MAX_INTERPOLATORS; i++)
    {
        if (interp[i].iX == NULL)
            continue;
        glColor3f(1, 1, 1);
        draw_interpolator(i, Us, N);
    }

    glutSwapBuffers();
}

int main(int argc, char **argv)
{
    glutInit(&argc, argv);

    init();

    // set up callback functions
    glutMouseFunc(mouse);
    glutKeyboardFunc(keyboard);
    glutMotionFunc(motion);
    glutReshapeFunc(reshape);
    glutDisplayFunc(display);
    glutMainLoop();

    free_interpolator2ds();

    return 0;
}