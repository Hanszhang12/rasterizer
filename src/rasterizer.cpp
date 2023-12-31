#include "rasterizer.h"

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

    sample_buffer[y * width * this->sample_rate + (x * sqrt(this->sample_rate))] = c;

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
      // TODO: Task 1: Implement basic triangle rasterization here, no supersampling
      x0 = x0 * (float)sqrt(this->sample_rate);
      x1 = x1 * (float)sqrt(this->sample_rate);
      x2 = x2 * (float)sqrt(this->sample_rate);

      y0 = y0 * (float)sqrt(this->sample_rate);
      y1 = y1 * (float)sqrt(this->sample_rate);
      y2 = y2 * (float)sqrt(this->sample_rate);

      float dx_0 = x1 - x0;
      float dy_0 = y1 - y0;

      float dx_1 = x2 - x1;
      float dy_1 = y2 - y1;

      float dx_2 = x0 - x2;
      float dy_2 = y0 - y2;


      for (int x = 0; x < width; ++x) {
        for (int y = 0; y < height; ++y) {

          for (int mini_x = 0; mini_x < (int)sqrt(this->sample_rate); ++mini_x) {
            for (int mini_y = 0; mini_y < (int)sqrt(this->sample_rate); ++mini_y) {
              float point_x = x * sqrt(this->sample_rate) + mini_x + 0.5;
              float point_y = y * sqrt(this->sample_rate) + mini_y + 0.5;
              float l0 = -(point_x - x0) * dy_0 + (point_y - y0) * dx_0;
              float l1 = -(point_x - x1) * dy_1 + (point_y - y1) * dx_1;
              float l2 = -(point_x - x2) * dy_2 + (point_y - y2) * dx_2;
              if ((l0 >= 0 && l1 >= 0 && l2 >= 0) || (l0 <= 0 && l1 <= 0 && l2 <= 0)) {
                sample_buffer[((y * sqrt(this->sample_rate) + mini_y) * width * sqrt(this->sample_rate)) + (x * sqrt(this->sample_rate) + mini_x)] = color;
              }
            }
          }

        }
      }

      // TODO: Task 2: Update to implement super-sampled rasterization
    }


    void RasterizerImp::rasterize_interpolated_color_triangle(float x0, float y0, Color c0,
    float x1, float y1, Color c1,
    float x2, float y2, Color c2)
  {
    // TODO: Task 4: Rasterize the triangle, calculating barycentric coordinates and using them to interpolate vertex colors across the triangle
    // Hint: You can reuse code from rasterize_triangle

    x0 = x0 * (float)sqrt(this->sample_rate);
    x1 = x1 * (float)sqrt(this->sample_rate);
    x2 = x2 * (float)sqrt(this->sample_rate);

    y0 = y0 * (float)sqrt(this->sample_rate);
    y1 = y1 * (float)sqrt(this->sample_rate);
    y2 = y2 * (float)sqrt(this->sample_rate);

    float dx_0 = x1 - x0;
    float dy_0 = y1 - y0;

    float dx_1 = x2 - x1;
    float dy_1 = y2 - y1;

    float dx_2 = x0 - x2;
    float dy_2 = y0 - y2;


    for (int x = 0; x < width; ++x) {
      for (int y = 0; y < height; ++y) {

        for (int mini_x = 0; mini_x < (int)sqrt(this->sample_rate); ++mini_x) {
          for (int mini_y = 0; mini_y < (int)sqrt(this->sample_rate); ++mini_y) {
            float point_x = x * sqrt(this->sample_rate) + mini_x + 0.5;
            float point_y = y * sqrt(this->sample_rate) + mini_y + 0.5;
            float l0 = -(point_x - x0) * dy_0 + (point_y - y0) * dx_0;
            float l1 = -(point_x - x1) * dy_1 + (point_y - y1) * dx_1;
            float l2 = -(point_x - x2) * dy_2 + (point_y - y2) * dx_2;
            if ((l0 >= 0 && l1 >= 0 && l2 >= 0) || (l0 <= 0 && l1 <= 0 && l2 <= 0)) {
              float alpha =  (-(point_x - x1) * dy_1 + (point_y - y1) * dx_1) /(-(x0 - x1) * dy_1 + (y0 - y1) * dx_1);
              float beta =   (-(point_x - x2) * dy_2 + (point_y - y2) * dx_2)/ (-(x1 - x2) * dy_2 + (y1 - y2) * dx_2);
              float gamma = 1 - alpha - beta;
              Color col = c0 * alpha + c1 * beta + c2 * gamma;
              sample_buffer[((y * sqrt(this->sample_rate) + mini_y) * width * sqrt(this->sample_rate)) + (x * sqrt(this->sample_rate) + mini_x)] = col;
            }
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
    int sqrt_sampr = sqrt(get_sample_rate());

    x0 = x0 * sqrt_sampr;
    y0 = y0 * sqrt_sampr;
    x1 = x1 * sqrt_sampr;
    y1 = y1 * sqrt_sampr;
    x2 = x2 * sqrt_sampr;
    y2 = y2 * sqrt_sampr;

    int a_x = (int) (floor(x0));
    int a_y = (int) (floor(y0));
    int b_x = (int) (floor(x1));
    int b_y = (int) (floor(y1));
    int c_x = (int) (floor(x2));
    int c_y = (int) (floor(y2));

    //first, you have to ensure that the coordinates are within appropriate bounds
    if (a_x < 0 || a_x >= width * sqrt_sampr || a_y < 0 || a_y >= height * sqrt_sampr ||
    b_x < 0 || b_x >= width * sqrt_sampr || a_y < 0 || a_y >= height * sqrt_sampr ||
    c_x < 0 || c_x >= width * sqrt_sampr || c_y < 0 || c_y >= height * sqrt_sampr) {
      return;
    }
    //initialize sampling parameters
    int lower_x = std::min({a_x, b_x, c_x});
    int upper_x = std::max({a_x, b_x, c_x});
    int lower_y = std::min({a_y, b_y, c_y});
    int upper_y = std::max({a_y, b_y, c_y});

    for (int y = lower_y; y <= upper_y; y++) {
      for (int x = lower_x; x <= upper_x; x++) {
        //construct line equations
        float ac = -1*(x-a_x + 0.5)*(c_y-a_y) + (y-a_y + 0.5)*(c_x-a_x);
        float cb = -1*(x-c_x+ 0.5)*(b_y-c_y) + (y-c_y)*(b_x-c_x);
        float ba = -1*(x-b_x + 0.5)*(a_y-b_y) + (y-b_y + 0.5)*(a_x-b_x);

        float alpha = (-1*(x-x1 + 0.5)*(y2-y1)+(y-y1 + 0.5)*(x2-x1))/(-1*(x0-x1)*(y2-y1)+(y0-y1)*(x2-x1));
        float beta = (-1*(x-x2 + 0.5)*(y0-y2)+(y-y2 + 0.5)*(x0-x2))/(-1*(x1-x2)*(y0-y2)+(y1-y2)*(x0-x2));
        float gamma = 1-alpha-beta;


        float x_and1 = x == upper_x ? x : x + 1;
        float alpha_xand1 = (-1 * (x_and1 - x1 + 0.5) * (y2 - y1) + (y - y1 + 0.5) * (x2 - x1)) / (-1 * (x0 - x1) * (y2 - y1) + (y0 - y1) * (x2 - x1));
        float beta_xand1 = (-1 * (x_and1 - x2 + 0.5) * (y0 - y2) + (y - y2 + 0.5) * (x0 - x2)) / (-1 * (x1 - x2) * (y0 - y2) + (y1 - y2) * (x0 - x2));
        float gamma_xand1 = 1 - alpha_xand1 - beta_xand1;

        float y_and1 = y == upper_y ? y: y + 1;
        float alpha_yand1 = (-1 * (x - x1 + 0.5) * (y2 - y1) + (y_and1 - y1 + 0.5) * (x2 - x1)) / (-1 * (x0 - x1) * (y2 - y1) + (y0 - y1) * (x2 - x1));
        float beta_yand1 = (-1 * (x - x2 + 0.5) * (y0 - y2) + (y_and1 - y2 + 0.5) * (x0 - x2)) / (-1 * (x1 - x2) * (y0 - y2) + (y1 - y2) * (x0 - x2));
        float gamma_yand1 = 1 - alpha_yand1 - beta_yand1;

        Vector2D uv = alpha * Vector2D(u0, v0) + beta * Vector2D(u1, v1) + gamma * Vector2D(u2, v2);
        Vector2D uv_deltax = alpha_xand1 * Vector2D(u0, v0) + beta_xand1 * Vector2D(u1, v1) + gamma_xand1 * Vector2D(u2, v2);
        Vector2D uv_deltay = alpha_yand1 * Vector2D(u0, v0) + beta_yand1 * Vector2D(u1, v1) + gamma_yand1 * Vector2D(u2, v2);

        //form params struct and fill in approproate arguments
        SampleParams params;
        params.p_uv = uv;
        params.p_dx_uv = uv_deltax;
        params.p_dy_uv = uv_deltay;
        params.lsm = this->lsm;
        params.psm = this->psm;

        Color col = tex.sample(params);

        //check to see if in all three planes -> in triangle, or for some reason
        //in none of the three planes -> contradiction -> must be in triangle
        if ((ac >= 0 && ba >= 0 && cb >= 0) || (ac <= 0 && ba <= 0 && cb <= 0)){
          sample_buffer[y * width * sqrt(this->sample_rate) + x] = col;
        }
      }
    }
    return;
  }

  void RasterizerImp::set_sample_rate(unsigned int rate) {
    // TODO: Task 2: You may want to update this function for supersampling support

    this->sample_rate = rate;

    //added rate in resizing
    this->sample_buffer.resize(width * height * this->sample_rate, Color::White);
    clear_buffers();

  }


  void RasterizerImp::set_framebuffer_target(unsigned char* rgb_framebuffer,
    size_t width, size_t height)
  {
    // TODO: Task 2: You may want to update this function for supersampling support

    this->width = width;
    this->height = height;
    this->rgb_framebuffer_target = rgb_framebuffer;

    //added rate in resizing
    this->sample_buffer.resize(width * height * this->sample_rate, Color::White);
    clear_buffers();

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

    for (int x = 0; x < width; x++) {
      for (int y = 0; y < height; y++) {
        float temp_r = 0.0;
        float temp_g = 0.0;
        float temp_b = 0.0;
        for (int mini_x = 0; mini_x < sqrt(this->sample_rate); ++mini_x) {
          for (int mini_y = 0; mini_y < sqrt(this->sample_rate); ++mini_y) {
            Color col = sample_buffer[((y * sqrt(this->sample_rate) + mini_y) * width * sqrt(this->sample_rate)) + (x * sqrt(this->sample_rate) + mini_x)];
            temp_r += (float)(&col.r)[0];
            temp_g += (float)(&col.r)[1];
            temp_b += (float)(&col.r)[2];

          }
        }
        temp_r = temp_r / (float)this->sample_rate;
        temp_g = temp_g / (float)this->sample_rate;
        temp_b = temp_b / (float)this->sample_rate;

        Color col = Color(temp_r, temp_g, temp_b);

        for (int k = 0; k < 3; ++k) {
          this->rgb_framebuffer_target[3 * (y * width + x) + k] = (&col.r)[k] * 255;
        }
      }
    }
  }
  Rasterizer::~Rasterizer() { }


}// CGL
