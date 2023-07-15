#include "rasterizer.h"
#include <cmath>

using namespace std;

namespace CGL {

  RasterizerImp::RasterizerImp(PixelSampleMethod psm, LevelSampleMethod lsm,
    size_t width, size_t height,
    unsigned int sample_rate) {
    this->psm = psm;
    this->lsm = lsm;
    this->width = width;
    this->height = height;
    this->sample_rate = sample_rate;

    sample_buffer.resize(width * height * sample_rate, Color::White);
  }

  // Used by rasterize_point and rasterize_line
  void RasterizerImp::fill_pixel(size_t x, size_t y, Color c) {
    // TODO: Task 2: You might need to this function to fix points and lines (such as the black rectangle border in test4.svg)
    // NOTE: You are not required to implement proper supersampling for points and lines
    // It is sufficient to use the same color for all supersamples of a pixel for points and lines (not triangles)
    

    sample_buffer[y * width + x] = c;
  }

  // Rasterize a point: simple example to help you start familiarizing
  // yourself with the starter code.
  //
  void RasterizerImp::rasterize_point(float x, float y, Color color) {
    // fill in the nearest pixel
    int sx = (int)floor(x);
    int sy = (int)floor(y);

    // check bounds
    if (sx < 0 || sx >= width) return;
    if (sy < 0 || sy >= height) return;

    fill_pixel(sx, sy, color);
    return;
  }

  // Rasterize a line.
  void RasterizerImp::rasterize_line(float x0, float y0,
    float x1, float y1,
    Color color) {
    if (x0 > x1) {
      swap(x0, x1); swap(y0, y1);
    }

    float pt[] = { x0,y0 };
    float m = (y1 - y0) / (x1 - x0);
    float dpt[] = { 1,m };
    int steep = abs(m) > 1;
    if (steep) {
      dpt[0] = x1 == x0 ? 0 : 1 / abs(m);
      dpt[1] = x1 == x0 ? (y1 - y0) / abs(y1 - y0) : m / abs(m);
    }

    while (floor(pt[0]) <= floor(x1) && abs(pt[1] - y0) <= abs(y1 - y0)) {
      rasterize_point(pt[0], pt[1], color);
      pt[0] += dpt[0]; pt[1] += dpt[1];
    }
  }

  // Rasterize a triangle.
  void RasterizerImp::rasterize_triangle(float x0, float y0,
    float x1, float y1,
    float x2, float y2,
    Color color) {

    // boundary
    float minx = floor(min(min(x0, x1), x2));
    float miny = floor(min(min(y0, y1), y2));
    float maxx = ceil(max(max(x0, x1), x2));
    float maxy = ceil(max(max(y0, y1), y2));

    if(minx < 0 || maxx > width) return;
    if(miny < 0 || maxy > height) return;

    Vector2D v0(x1-x0,y1-y0);
    Vector2D v1(x2-x1,y2-y1);
    Vector2D v2(x0-x2,y0-y2);
    Vector2D p0, p1, p2;
    float x, y;
    float c0, c1, c2; // cross product
    float inside;

    // TODO: Task 1: Implement basic triangle rasterization here, no supersampling
    if(sample_rate == 1){

      for(int i = (int)minx; i < (int)maxx; ++i){
        x = (float)i + 0.5;
        for(int j = (int)miny; j < (int)maxy; ++j){
          y = (float)j + 0.5;
          inside = 0;
          p0.x = x - x0; 
          p0.y = y - y0;
          c0 = v0.x * p0.y - p0.x * v0.y;
          p1.x = x - x1;
          p1.y = y - y1;
          c1 = v1.x * p1.y - p1.x * v1.y;
          p2.x = x - x2;
          p2.y = y - y2;
          c2 = v2.x * p2.y - p2.x * v2.y;
          if((c0 >= 0 && c1 >= 0 && c2 >= 0) || (c0 <= 0 && c1 <= 0 && c2 <= 0)) inside = 1;
          if(inside != 0) rasterize_point(x, y, color);
        }
      }

    }else{
    // TODO: Task 2: Update to implement super-sampled rasterization

      float step = sqrt(sample_rate);
      for(int i = (int)minx; i < (int)maxx; ++i){
        for(int j = (int)miny; j < (int)maxy; ++j){
          inside = 0;
          // super-sampling within one pixel:
          for(int a = 0; a < step; ++a){
            for(int b = 0; b < step; ++b){
              
              x = (float)i + (a + 0.5)/step;
              y = (float)j + (b + 0.5)/step;

              p0.x = x - x0; 
              p0.y = y - y0;
              c0 = v0.x * p0.y - p0.x * v0.y;
              p1.x = x - x1;
              p1.y = y - y1;
              c1 = v1.x * p1.y - p1.x * v1.y;
              p2.x = x - x2;
              p2.y = y - y2;
              c2 = v2.x * p2.y - p2.x * v2.y;

              if((c0 >= 0 && c1 >= 0 && c2 >= 0) || (c0 <= 0 && c1 <= 0 && c2 <= 0)) inside++;
            }
          }

          if(inside != 0){
            inside /= sample_rate;
            Color fill_color = (1 - inside) * Color::White + inside * color;
            rasterize_point(x, y, fill_color);
          }

        }
      }

    }
  }


  void RasterizerImp::rasterize_interpolated_color_triangle(float x0, float y0, Color c0,
    float x1, float y1, Color c1,
    float x2, float y2, Color c2)
  {
    // TODO: Task 4: Rasterize the triangle, calculating barycentric coordinates and using them to interpolate vertex colors across the triangle
    // Hint: You can reuse code from rasterize_triangle

    // boundary
    float minx = floor(min(min(x0, x1), x2));
    float miny = floor(min(min(y0, y1), y2));
    float maxx = ceil(max(max(x0, x1), x2));
    float maxy = ceil(max(max(y0, y1), y2));

    if(minx < 0 || maxx > width) return;
    if(miny < 0 || maxy > height) return;

    float x, y;
    float l0, l1, l2;
    float inside;

    if(sample_rate == 1){

      for(int i = (int)minx; i < (int)maxx; ++i){
        x = (float)i + 0.5;
        for(int j = (int)miny; j < (int)maxy; ++j){
          y = (float)j + 0.5;

          l0 = ((y1-y2)*(x-x2)+(x2-x1)*(y-y2))/((y1-y2)*(x0-x2)+(x2-x1)*(y0-y2));
          l1 = ((y2-y0)*(x-x2)+(x0-x2)*(y-y2))/((y1-y2)*(x0-x2)+(x2-x1)*(y0-y2));
          l2 = 1 - l0 - l1;
          
          if(l0 >= 0 && l1 >= 0 && l2 >= 0){
            Color fill_color;
            fill_color = l0 * c0 + l1 * c1 + l2 * c2;
            rasterize_point(x, y, fill_color);
          }
        }
      }
    }else{

      float step = sqrt(sample_rate);
      for(int i = (int)minx; i < (int)maxx; ++i){
        for(int j = (int)miny; j < (int)maxy; ++j){
          // super-sampling within one pixel:
          Color fill_color(0, 0, 0);
          inside = 0;
          for(int a = 0; a < step; ++a){
            for(int b = 0; b < step; ++b){
              
              x = (float)i + (a + 0.5)/step;
              y = (float)j + (b + 0.5)/step;

          
              l0 = ((y1-y2)*(x-x2)+(x2-x1)*(y-y2))/((y1-y2)*(x0-x2)+(x2-x1)*(y0-y2));
              l1 = ((y2-y0)*(x-x2)+(x0-x2)*(y-y2))/((y1-y2)*(x0-x2)+(x2-x1)*(y0-y2));
              l2 = 1 - l0 - l1;
          
              if(l0 >= 0 && l1 >= 0 && l2 >= 0){
                inside++;
                fill_color += l0 * c0 + l1 * c1 + l2 * c2;
              }
            }
          }
          if(inside != 0){
            fill_color *= 1/inside;
            rasterize_point(x, y, fill_color);
          }
        }
      }
    }
  }


  void RasterizerImp::rasterize_textured_triangle(float x0, float y0, float u0, float v0,
    float x1, float y1, float u1, float v1,
    float x2, float y2, float u2, float v2,
    Texture& tex)
  {
    // TODO: Task 5: Fill in the SampleParams struct and pass it to the tex.sample function.
    // TODO: Task 6: Set the correct barycentric differentials in the SampleParams struct.
    // Hint: You can reuse code from rasterize_triangle/rasterize_interpolated_color_triangle

    // boundary
    float minx = floor(min(min(x0, x1), x2));
    float miny = floor(min(min(y0, y1), y2));
    float maxx = ceil(max(max(x0, x1), x2));
    float maxy = ceil(max(max(y0, y1), y2));

    if(minx < 0 || maxx > width) return;
    if(miny < 0 || maxy > height) return;

    float x, y;
    // barycentric coefficients
    float l0, l1, l2;
    float a, b, c;
    float inside;

    SampleParams sp;

    if(sample_rate == 1){

      for(int i = (int)minx; i < (int)maxx; ++i){
        x = (float)i + 0.5;
        for(int j = (int)miny; j < (int)maxy; ++j){
          y = (float)j + 0.5;
          
          // barycentric of (x, y)
          l0 = ((y1-y2)*(x-x2)+(x2-x1)*(y-y2))/((y1-y2)*(x0-x2)+(x2-x1)*(y0-y2));
          l1 = ((y2-y0)*(x-x2)+(x0-x2)*(y-y2))/((y1-y2)*(x0-x2)+(x2-x1)*(y0-y2));
          l2 = 1 - l0 - l1;
          
          if(l0 >= 0 && l1 >= 0 && l2 >= 0){

            sp.p_uv.x = l0 * u0 + l1 * u1 + l2 * u2;
            sp.p_uv.y = l0 * v0 + l1 * v1 + l2 * v2;

            // barycentric of (x+1, y)
            a = ((y1-y2)*(x+1-x2)+(x2-x1)*(y-y2))/((y1-y2)*(x0-x2)+(x2-x1)*(y0-y2));
            b = ((y2-y0)*(x+1-x2)+(x0-x2)*(y-y2))/((y1-y2)*(x0-x2)+(x2-x1)*(y0-y2));
            c = 1 - a - b;

            sp.p_dx_uv.x = a * u0 + b * u1 + c * u2;
            sp.p_dx_uv.y = a * v0 + b * v1 + c * v2;

            // barycentric of (x, y+1)
            a = ((y1-y2)*(x-x2)+(x2-x1)*(y+1-y2))/((y1-y2)*(x0-x2)+(x2-x1)*(y0-y2));
            b = ((y2-y0)*(x-x2)+(x0-x2)*(y+1-y2))/((y1-y2)*(x0-x2)+(x2-x1)*(y0-y2));
            c = 1 - a - b;
            
            sp.p_dy_uv.x = a * u0 + b * u1 + c * u2;
            sp.p_dy_uv.y = a * v0 + b * v1 + c * v2;

            rasterize_point(x, y, tex.sample(sp));
          }
        }
      }

    }else{

      float step = sqrt(sample_rate);
      for(int i = (int)minx; i < (int)maxx; ++i){
        for(int j = (int)miny; j < (int)maxy; ++j){
          // super-sampling within one pixel:
          Color fill_color(0, 0, 0);
          inside = 0;
          for(int m = 0; m < step; ++m){
            for(int n = 0; n < step; ++n){
              
              x = (float)i + (m + 0.5)/step;
              y = (float)j + (n + 0.5)/step;
          
              l0 = ((y1-y2)*(x-x2)+(x2-x1)*(y-y2))/((y1-y2)*(x0-x2)+(x2-x1)*(y0-y2));
              l1 = ((y2-y0)*(x-x2)+(x0-x2)*(y-y2))/((y1-y2)*(x0-x2)+(x2-x1)*(y0-y2));
              l2 = 1 - l0 - l1;

              if(l0 >= 0 && l1 >= 0 && l2 >= 0){

                sp.p_uv.x = l0 * u0 + l1 * u1 + l2 * u2;
                sp.p_uv.y = l0 * v0 + l1 * v1 + l2 * v2;

                // barycentric of (x+1, y)
                a = ((y1-y2)*(x+1-x2)+(x2-x1)*(y-y2))/((y1-y2)*(x0-x2)+(x2-x1)*(y0-y2));
                b = ((y2-y0)*(x+1-x2)+(x0-x2)*(y-y2))/((y1-y2)*(x0-x2)+(x2-x1)*(y0-y2));
                c = 1 - a - b;

                sp.p_dx_uv.x = a * u0 + b * u1 + c * u2;
                sp.p_dx_uv.y = a * v0 + b * v1 + c * v2;

                // barycentric of (x, y+1)
                a = ((y1-y2)*(x-x2)+(x2-x1)*(y+1-y2))/((y1-y2)*(x0-x2)+(x2-x1)*(y0-y2));
                b = ((y2-y0)*(x-x2)+(x0-x2)*(y+1-y2))/((y1-y2)*(x0-x2)+(x2-x1)*(y0-y2));
                c = 1 - a - b;
                
                sp.p_dy_uv.x = a * u0 + b * u1 + c * u2;
                sp.p_dy_uv.y = a * v0 + b * v1 + c * v2;

                fill_color += tex.sample(sp);
                inside++;
              }
            }
          }
          if(inside != 0){
            fill_color *= 1/inside;
            rasterize_point(x, y, fill_color);
          }
          
        }
      }

    }

  }

  void RasterizerImp::set_sample_rate(unsigned int rate) {
    // TODO: Task 2: You may want to update this function for supersampling support

    this->sample_rate = rate;


    this->sample_buffer.resize(width * height, Color::White);
  }


  void RasterizerImp::set_framebuffer_target(unsigned char* rgb_framebuffer,
    size_t width, size_t height)
  {
    // TODO: Task 2: You may want to update this function for supersampling support

    this->width = width;
    this->height = height;
    this->rgb_framebuffer_target = rgb_framebuffer;


    this->sample_buffer.resize(width * height, Color::White);
  }


  void RasterizerImp::clear_buffers() {
    std::fill(rgb_framebuffer_target, rgb_framebuffer_target + 3 * width * height, 255);
    std::fill(sample_buffer.begin(), sample_buffer.end(), Color::White);
  }


  // This function is called at the end of rasterizing all elements of the
  // SVG file.  If you use a supersample buffer to rasterize SVG elements
  // for antialising, you could use this call to fill the target framebuffer
  // pixels from the supersample buffer data.
  //
  void RasterizerImp::resolve_to_framebuffer() {
    // TODO: Task 2: You will likely want to update this function for supersampling support


    for (int x = 0; x < width; ++x) {
      for (int y = 0; y < height; ++y) {
        Color col = sample_buffer[y * width + x];
      
        for (int k = 0; k < 3; ++k) {
          this->rgb_framebuffer_target[3 * (y * width + x) + k] = (&col.r)[k] * 255;
        }
      }
    }

  }

  Rasterizer::~Rasterizer() { }


}// CGL
