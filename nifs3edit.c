#include <stdlib.h>
#include <memory.h>
#include <math.h>
#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <ctype.h>

#include <GL/gl.h>
#include <GL/glut.h>

#define STB_IMAGE_IMPLEMENTATION
#include <stb_image.h>

#define min(a, b) ((a) < (b) ? (a) : (b))
#define max(a, b) ((b) < (a) ? (a) : (b))

//////////// Natural cubic spline interpolation //////////////
double *get_diff_polys(const double *x, const double *y, int n)
{
    double *y2 = malloc(sizeof(double) * (n - 1));
    for (int i = 0; i < n - 1; i++)
    {
        y2[i] = (y[i] - y[i + 1]) / (x[i + 1] - x[i]);
    }

    double *y3 = malloc(sizeof(double) * (n - 2));
    for (int i = 0; i < n - 2; i++)
    {
        y3[i] = (y2[i] - y2[i + 1]) / (x[i + 2] - x[i]);
    }

    free(y2);
    return y3;
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

    double *M = malloc(sizeof(double) * (n + 1));

    double *q = malloc(sizeof(double) * n);
    double *u = malloc(sizeof(double) * n);

    double *d = get_diff_polys(x, y, n);
    for (int i = 0; i < n - 2; i++)
        d[i] *= 6;

    q[0] = u[0] = 0;

    for (int i = 1; i <= n - 1; i++)
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
        M[i] = u[i] + q[i] * M[i + 1];

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
    {
        printf("x out of range\n");
        exit(-1);
    }

    int i = 0;
    while (i < interp->n && x >= interp->x[i])
        i++;

    i = min(max(1, i), interp->n);

    double h = interp->x[i] - interp->x[i - 1];
    double t1 = interp->x[i] - x;
    double t2 = x - interp->x[i - 1];

    return (1 / h) * (interp->M[i - 1] / 6 * t1 * t1 * t1 +
                      interp->M[i] / 6 * t2 * t2 * t2 +
                      (interp->y[i - 1] - interp->M[i - 1] / 6 * h * h) * t1 +
                      (interp->y[i] - interp->M[i] / 6 * h * h) * t2);
}

///////////// 2D Interpolation //////////////
typedef struct
{
    nifs3_t *iX, *iY;
    double xMax, xMin, yMax, yMin;

    int n;
    double *u; // interpolation points
} nifs3_2d_t;

#define MAX_INTERPOLATORS 128
nifs3_2d_t interp[MAX_INTERPOLATORS];

void draw_nifs3_2d(int inp)
{
    for (int i = 0; i < interp[inp].n; i++)
    {
        double x = nifs3_get(interp[inp].iX, interp[inp].u[i]);
        double y = nifs3_get(interp[inp].iY, interp[inp].u[i]);
        glVertex2d(x, y);
    }
}

void free_nifs3_2d(int i)
{
    nifs3_free(interp[i].iX);
    nifs3_free(interp[i].iY);
    free(interp[i].u);
    memset(&interp[i], 0, sizeof(nifs3_2d_t));
}

void cleanup_nifs3_2d()
{
    for (int i = 0; i < MAX_INTERPOLATORS; i++)
        free_nifs3_2d(i);
}

void set_nifs3_2d_interpolation_pts(int i, const double *u, int n)
{
    assert(interp[i].iX != NULL);
    interp[i].n = n;
    free(interp[i].u);
    interp[i].u = malloc(sizeof(double) * n);
    memcpy(interp[i].u, u, sizeof(double) * n);
}

int create_nifs3_2d(const double *x, const double *y, const double *t, int n)
{
    for (int i = 0; i < MAX_INTERPOLATORS; i++)
    {
        if (interp[i].iX != NULL)
            continue;

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

        set_nifs3_2d_interpolation_pts(i, t, n);

        return i;
    }

    assert(false && "No free 2d interpolators");
    return -1;
}

///////////// Loading 2d interpolators from file //////////////
double *get_line_array(FILE *fh, int *count)
{
    int real_count = 0;
    int capacity = 32;
    double *arr = malloc(sizeof(double) * capacity);

    while (fscanf(fh, "%lf", &arr[real_count]) > 0)
    {
        real_count++;

        int c;
        while (isspace(c = fgetc(fh)) && c != '\n' && c != EOF)
        {
        }

        if (c == '\n' || c == EOF)
            break;

        ungetc(c, fh);

        if (real_count < capacity)
            continue;

        capacity *= 2;
        arr = realloc(arr, sizeof(double) * capacity);
    }

    arr = realloc(arr, sizeof(double) * real_count);
    *count = real_count;
    return arr;
}

bool load_from_file(const char *path)
{
    cleanup_nifs3_2d();

    FILE *fh = fopen(path, "r");
    if (fh == NULL)
    {
        printf("Failed to open file %s\n", path);
        return false;
    }

    while (!feof(fh))
    {
        // get to next non-whitespace character
        int c;
        while (isspace(c = fgetc(fh)) && c != EOF)
        {
        }
        if (c == EOF)
            break;

        ungetc(c, fh);

        int nx, ny, nt, nu;
        double *x = get_line_array(fh, &nx);
        double *y = get_line_array(fh, &ny);
        double *t = get_line_array(fh, &nt);
        double *u = get_line_array(fh, &nu);

        if (nx == 0 && ny == 0 && nt == 0 && nu == 0)
            break;

        printf("nx: %d, ny: %d, nt: %d, nu: %d\n", nx, ny, nt, nu);

        if (nx != ny || nx != nt || nx == 0 || nu == 0)
        {
            printf("Invalid file format\n");
            cleanup_nifs3_2d();
            return false;
        }

        int i = create_nifs3_2d(x, y, t, nt);
        set_nifs3_2d_interpolation_pts(i, u, nu);
    }

    return true;
}

void save_to_file(const char *path)
{
    FILE *fh = fopen(path, "w");
    if (fh == NULL)
    {
        printf("Failed to open file %s\n", path);
        return;
    }

    for (int i = 0; i < MAX_INTERPOLATORS; i++)
    {
        if (interp[i].iX == NULL)
            continue;

        for (int j = 0; j < interp[i].iX->n; j++)
            fprintf(fh, "%lf ", interp[i].iX->y[j]); // Xs
        fprintf(fh, "\n");

        for (int j = 0; j < interp[i].iY->n; j++)
            fprintf(fh, "%lf ", interp[i].iY->y[j]); // Ys
        fprintf(fh, "\n");

        for (int j = 0; j < interp[i].iX->n; j++)
            fprintf(fh, "%lf ", interp[i].iX->x[j]); // Ts
        fprintf(fh, "\n");

        for (int j = 0; j < interp[i].n; j++)
            fprintf(fh, "%lf ", interp[i].u[j]); // Us
        fprintf(fh, "\n");

        fprintf(fh, "\n");
    }
}

///////////// APPLICATION //////////////
void linspace(double start, double end, int n, double *data)
{
    double step = (end - start) / (n - 1);
    for (int i = 0; i < n; i++)
        data[i] = start + i * step;
}

double *alloc_linspace(double start, double end, int n)
{
    double *data = malloc(sizeof(double) * n);
    linspace(start, end, n, data);
    return data;
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
    bool ok = load_from_file("zadanie7.data");
    if (!ok)
        exit(1);

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
    stbi_uc *data = stbi_load("image_transparent.png", &scene_data.imageW, &scene_data.imageH, &components, 4);
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

    glColor3f(0, 0, 1);
    for (char *s = str; *s; s++)
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, *s);
    glColor3f(1, 1, 1);

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

        glColor3f(1, 0, 0);
        glBegin(GL_LINE_STRIP);
        draw_nifs3_2d(i);
        glEnd();

        glColor3f(0, 1, 0);
        glPointSize(2);
        glBegin(GL_POINTS);
        draw_nifs3_2d(i);
        glEnd();

        glColor3f(0, 0, 1);
        glPointSize(4);
        glBegin(GL_POINTS);
        for (int j = 0; j < interp[i].iX->n; j++)
            glVertex2d(interp[i].iX->y[j], interp[i].iY->y[j]);
        glEnd();
        glColor3f(1, 1, 1);
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

    cleanup_nifs3_2d();

    return 0;
}