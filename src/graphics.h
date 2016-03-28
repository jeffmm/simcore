#ifndef _GRAPHICS_H
#define _GRAPHICS_H

#include <vector>
#include <GL/glew.h>
#define GLEW_STATIC
#include <GLFW/glfw3.h>
#include <string>
#include "auxiliary.h"

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
                                // "cube" and "sphere"

 public:
    void Init(int n_dim, double **h, double background); // Init. Must always be called.
    // n_dim: number of dimensions (2 or 3 currently supported)
    // h: n_dim x n_dim unit cell matrix 
    void Clear();

    void DrawLoop(graph_struct g_struct);
    void DrawLoop(int n_bonds, double **h,
                  double **r, double **u, double *length,
                  int n_sites, double **r_site,
                  double *sphero_diameter,
                  double sphere_diameter); // Same as draw, but loop indefinitely
    void DrawLoop(int n_bonds, double **h,
                  double **r, double **u, double *length, double *sphero_diameter);
    void DrawLoop(int n_bonds, double **h,
                        double **r, double **u, double *length, double *sphero_diameter,
                        double m_rad, double d_rad, double m_d_dist);
    void Draw(graph_struct g_struct);
    void Draw(int n_bonds, double **h,
              double **r, double **u, double *length, double *sphero_diameter); // Draw current simulation frame
    // n_bonds: number of spherocylinders
    // h: unit cell matrix
    // r: center of mass positions of spherocylinders
    // u: orientation vectors for spherocylinders
    // length: length of spherocylinders
    void Draw(int n_bonds, double **h,
              double **r, double **u, double *length,
              int n_sites, double **r_site,
              double *sphero_diameter,
              double sphere_diameter); // Draw current simulation frame with spheres too!
    // n_bonds: number of spherocylinders
    // h: unit cell matrix
    // r: center of mass positions of spherocylinders
    // u: orientation vectors for spherocylinders
    // length: length of spherocylinders
    // n_sites: actual number of sites (2 * n_bonds + n_spheres) FIXME
    // r_site: site positions matrix
    // sphere_diameter: diameter of spheres! //fixme, could be carried along 
    void Draw(int n_bonds, double **h,
                    double **r, double **u, double *length, double *sphero_diameter,
                    double m_rad, double d_rad, double m_d_dist); 
    void SetBoundaryType(std::string boundary_type); //Set the boundary type to draw
    
 private:
    void InitColormap(); 
    void InitDiscoRectangle();
    void InitSpheroCylinder();
    void InitWindow(int n_dim, double **h);
    void Init2dWindow(double **h);
    void Init3dWindow(double **h);
    void KeyInteraction();
    void DrawBox(double **h); // Wrapper for rectangular boundaries
    void DrawSquare(double **h); //boundary square
    void DrawCube(double **h); // boundary cube
    void DrawWireSphere(double r, int lats, int longs); // Slow, not for mass drawing, just for boundary
    void DrawBoundary(double **h); 
    void DrawBoundary(double **h, double m_rad, double d_rad, double m_d_dist); // used for snowman
    void DrawSnowman(double m_rad, double d_rad, double m_d_dist); // Draws snowman boundary
    void Draw2dSnowman(double m_rad, double d_rad, double m_d_dist); // Draws snowman boundary
    void Draw3dSnowman(double m_rad, double d_rad, double m_d_dist); // Draws snowman boundary
    void DrawSpheros(int n_spheros, double **r, double **u, double *length,
                     double *sphero_diameter); // draw spherocylinders (3d)
    void DrawDiscorectangles(int n_spheros, double **r, double **u, double *length, double *sphero_diameter); // draw solid discorectangles (2d spherocylinders)
    void DrawSpheres(int n_spheres, double **r, double sphere_diameter); // Draw solid spheres
    void UpdateWindow(); // update window parameters in case of resize
    void Draw2d(int n_bonds, double **h, double **r, double **u, double *length, double *sphero_diameter);
    void Draw2d(int n_bonds, double **h, double **r, double **u, double *length, 
                      double *sphero_diameter, double m_rad, double d_rad, double m_d_dist);
    void Draw3d(int n_bonds, double **h, double **r, double **u, double *length, double *sphero_diameter, double m_rad,
                double d_rad, double m_d_dist);
    void Draw3d(int n_bonds, double **h, double **r, double **u, double *length,
                int n_sites, double **r_site, double *sphero_diameter, double sphere_diameter);
    void Draw3d(int n_bonds, double **h, double **r, double **u, double *length, double *sphero_diameter);
};

#endif
