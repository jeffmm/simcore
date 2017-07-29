#ifndef _GRAPHICS_H
#define _GRAPHICS_H

#ifndef NOGRAPH

#include <vector>
#include <GL/glew.h>
#define GLEW_STATIC
#include <GLFW/glfw3.h>
#include <string>
#define SQR(x)             ((x) * (x))
#define ABS(x)             ((x) < 0 ? -(x) : (x))
#define MAX(x,y)           ((x) > (y) ? (x) : (y))
#include "definitions.h"

// Structure that holds all of the relevant drawing data
struct point {
	GLfloat x;
	GLfloat y;
	GLfloat s;
	GLfloat t;
};

class GraphicsPrimitive {
 public:
    GLuint vertex_buffer_, element_buffer_;
    GLuint vertex_shader_, fragment_shader_, program_;
    GLfloat *vertex_buffer_data_;
    GLushort *element_buffer_data_;

    struct {
        GLint half_l;
        GLint diameter;
    } uniforms_;

    struct {
        GLint position;
    } attributes_;

    int n_coords_;
    int n_triangles_;

    void MakeProgram();
    GLuint MakeBuffer(GLenum target, const void *buffer_data, GLsizei buffer_size);
    void MakeVertexBuffer();
    void MakeElementBuffer();
    void ShowInfoLog(GLuint object, PFNGLGETSHADERIVPROC glGet__iv,
                     PFNGLGETSHADERINFOLOGPROC glGet__InfoLog);
};

/* The <graphics_parameters> structure contains a variety of graphics parameters. */
class Graphics {
 public:
    GLFWwindow *window_; // Our actual window
    int windx_; // window x dimension in pixels
    int windy_; // window y dimension in pixels

    int n_dim_; // number of physical dimensions we're drawing
    double *unit_cell_;
    double z_correct_; // Used to recenter graphics
    space_struct * space_;
    std::vector<graph_struct*> * graph_array_;

    GraphicsPrimitive discorectangle_; // 2d spherocylinder
    GraphicsPrimitive spherocylinder_; // actual spherocylinder

    GLfloat xAngle_, yAngle_, zAngle_; // system rotation
    GLfloat xyzScale_; // zoom
    GLfloat xTrans_, yTrans_, zTrans_; // translation of system
    GLfloat zNear_, zFar_, cNear_, cFar_; // camera positions and clipping planes
    GLfloat cxAngle_, cyAngle_; // clipping plane angles
    GLfloat sphere_color_[3], box_color_[3], background_color_[3];
    GLfloat alpha_; // transparency (1 opaque, 0 invisible)

    static const int n_rgb_ = 384;
    static const int shift_ = 84;
    GLfloat colormap_[n_rgb_*4];

    int color_switch_; // which coloring algorithm to use
    int keep_going_; // Allow to draw same configuration in loop
    std::string boundary_type_; // boundary type to draw. currently supports
                                // "box" and "sphere"

 public:
    void Init(std::vector<graph_struct*> * const graph_array, space_struct * s_struct, double background); // Init. Must always be called.
    void Clear();

    void DrawLoop();
    void Draw();
    void SetBoundaryType(std::string boundary_type); //Set the boundary type to draw
    
 private:
    void InitColormap(); 
    void InitDiscoRectangle();
    void InitSpheroCylinder();
    void InitWindow();
    void Init2dWindow();
    void Init3dWindow();
    void InitText();
    void KeyInteraction();
    void DrawBox(); // Wrapper for rectangular boundaries
    void DrawSquare(); //boundary square
    void DrawCube(); // boundary cube
    void DrawWireSphere(double r, int lats, int longs); // Slow, not for mass drawing, just for boundary
    void DrawBoundary(); 
    void DrawBudding();
    void Draw2dBudding(); // Draws budding boundary
    void Draw3dBudding(); // Draws budding boundary (3d)
    void DrawSpheros(); // draw spherocylinders (3d)
    void DrawDiscorectangles(); // draw solid discorectangles (2d spherocylinders)
    void UpdateWindow(); // update window parameters in case of resize
    void Draw2d();
    void Draw3d();
    void DrawText();
};

#endif

#endif
