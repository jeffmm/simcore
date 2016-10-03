#ifndef _GRAPHICS_H
#define _GRAPHICS_H

#ifndef NOGRAPH

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
    double **unit_cell_;
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
                                // "cube" and "sphere"

 public:
    void Init(std::vector<graph_struct*> * const graph_array, space_struct * s_struct, double background); // Init. Must always be called.
    // n_dim: number of dimensions (2 or 3 currently supported)
    // h: n_dim x n_dim unit cell matrix 
    void Clear();

    void DrawLoop();
    //void DrawLoop(graph_struct g_struct);
    //void DrawLoop(int n_bonds,
                  //double **r, double **u, double *length,
                  //int n_sites, double **r_site,
                  //double *sphero_diameter,
                  //double sphere_diameter); // Same as draw, but loop indefinitely
    //void DrawLoop(int n_bonds,
                  //double **r, double **u, double *length, double *sphero_diameter);
    //void DrawLoop(int n_bonds,
                        //double **r, double **u, double *length, double *sphero_diameter,
                        //double m_rad, double d_rad, double m_d_dist);
    void Draw();
    //void Draw(graph_struct g_struct);
    //void Draw(int n_bonds,
              //double **r, double **u, double *length, double *sphero_diameter); // Draw current simulation frame
    // n_bonds: number of spherocylinders
    // h: unit cell matrix
    // r: center of mass positions of spherocylinders
    // u: orientation vectors for spherocylinders
    // length: length of spherocylinders
    //void Draw(int n_bonds,
              //double **r, double **u, double *length,
              //int n_sites, double **r_site,
              //double *sphero_diameter,
              //double sphere_diameter); // Draw current simulation frame with spheres too!
    // n_bonds: number of spherocylinders
    // h: unit cell matrix
    // r: center of mass positions of spherocylinders
    // u: orientation vectors for spherocylinders
    // length: length of spherocylinders
    // n_sites: actual number of sites (2 * n_bonds + n_spheres) FIXME
    // r_site: site positions matrix
    // sphere_diameter: diameter of spheres! //fixme, could be carried along 
    //void Draw(int n_bonds,
                    //double **r, double **u, double *length, double *sphero_diameter,
                    //double m_rad, double d_rad, double m_d_dist); 
    void SetBoundaryType(std::string boundary_type); //Set the boundary type to draw
    
 private:
    void InitColormap(); 
    void InitDiscoRectangle();
    void InitSpheroCylinder();
    void InitWindow();
    void Init2dWindow();
    void Init3dWindow();
    void KeyInteraction();
    void DrawBox(); // Wrapper for rectangular boundaries
    void DrawSquare(); //boundary square
    void DrawCube(); // boundary cube
    void DrawWireSphere(double r, int lats, int longs); // Slow, not for mass drawing, just for boundary
    void DrawBoundary(); 
    //void DrawBoundary(, double m_rad, double d_rad, double m_d_dist); // used for snowman
    //void DrawSnowman(double m_rad, double d_rad, double m_d_dist); // Draws snowman boundary
    void DrawSnowman();
    void Draw2dSnowman(); // Draws snowman boundary
    void Draw3dSnowman(); // Draws snowman boundary
    void DrawSpheros(); // draw spherocylinders (3d)
    void DrawDiscorectangles(); // draw solid discorectangles (2d spherocylinders)
    //void DrawSpheres(int n_spheres, double **r, double sphere_diameter); // Draw solid spheres
    void UpdateWindow(); // update window parameters in case of resize
    void Draw2d();
    void Draw3d();
    //void Draw2d(int n_bonds, double **r, double **u, double *length, double *sphero_diameter);
    //void Draw2d(int n_bonds, double **r, double **u, double *length, 
                      //double *sphero_diameter, double m_rad, double d_rad, double m_d_dist);
    //void Draw3d(int n_bonds, double **r, double **u, double *length, double *sphero_diameter, double m_rad,
                //double d_rad, double m_d_dist);
    //void Draw3d(int n_bonds, double **r, double **u, double *length,
                //int n_sites, double **r_site, double *sphero_diameter, double sphere_diameter);
    //void Draw3d(int n_bonds, double **r, double **u, double *length, double *sphero_diameter);
};

#endif

#endif
