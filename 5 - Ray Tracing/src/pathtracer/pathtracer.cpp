#include "pathtracer.h"

#include "scene/light.h"
#include "scene/sphere.h"
#include "scene/triangle.h"


using namespace CGL::SceneObjects;

namespace CGL {

PathTracer::PathTracer() {
  gridSampler = new UniformGridSampler2D();
  hemisphereSampler = new UniformHemisphereSampler3D();

  tm_gamma = 2.2f;
  tm_level = 1.0f;
  tm_key = 0.18;
  tm_wht = 5.0f;
}

PathTracer::~PathTracer() {
  delete gridSampler;
  delete hemisphereSampler;
}

void PathTracer::set_frame_size(size_t width, size_t height) {
  sampleBuffer.resize(width, height);
  sampleCountBuffer.resize(width * height);
}

void PathTracer::clear() {
  bvh = NULL;
  scene = NULL;
  camera = NULL;
  sampleBuffer.clear();
  sampleCountBuffer.clear();
  sampleBuffer.resize(0, 0);
  sampleCountBuffer.resize(0, 0);
}

void PathTracer::write_to_framebuffer(ImageBuffer &framebuffer, size_t x0,
                                      size_t y0, size_t x1, size_t y1) {
  sampleBuffer.toColor(framebuffer, x0, y0, x1, y1);
}

Vector3D
PathTracer::estimate_direct_lighting_hemisphere(const Ray &r,
                                                const Intersection &isect) {
  // Estimate the lighting from this intersection coming directly from a light.
  // For this function, sample uniformly in a hemisphere.

  // Note: When comparing Cornel Box (CBxxx.dae) results to importance sampling, you may find the "glow" around the light source is gone.
  // This is totally fine: the area lights in importance sampling has directionality, however in hemisphere sampling we don't model this behaviour.

  // make a coordinate system for a hit point
  // with N aligned with the Z direction.
  Matrix3x3 o2w;
  make_coord_space(o2w, isect.n);
  Matrix3x3 w2o = o2w.T();

  // w_out points towards the source of the ray (e.g.,
  // toward the camera if this is a primary ray)
  const Vector3D hit_p = r.o + r.d * isect.t;
  const Vector3D w_out = w2o * (-r.d);

  // This is the same number of total samples as
  // estimate_direct_lighting_importance (outside of delta lights). We keep the
  // same number of samples for clarity of comparison.
  int num_samples = scene->lights.size() * ns_area_light;
  Vector3D L_out;

  // TODO (Part 3): Write your sampling loop here
  // TODO BEFORE YOU BEGIN
  // UPDATE `est_radiance_global_illumination` to return direct lighting instead of normal shading 

  for (int i = 0; i < num_samples; ++i) {
		Vector3D wi = hemisphereSampler->get_sample(); // sample wi in local coordinate
		Vector3D wi_world = o2w * wi;
		double pdf = 1 / (2 * PI); // uniform pdf in hemisphere
		Ray r_(hit_p + EPS_D * wi_world, wi_world, r.max_t);
		Intersection is;
		if (bvh->intersect(r_, &is)){
		  L_out += is.bsdf->get_emission() * isect.bsdf->f(w_out, wi) * wi.z / pdf; // note: wi.z direction is aligned with the normal vector
    }
	}
	L_out /= num_samples;

  return L_out;

}

Vector3D
PathTracer::estimate_direct_lighting_importance(const Ray &r,
                                                const Intersection &isect) {
  // Estimate the lighting from this intersection coming directly from a light.
  // To implement importance sampling, sample only from lights, not uniformly in
  // a hemisphere.

  // make a coordinate system for a hit point
  // with N aligned with the Z direction.
  Matrix3x3 o2w;
  make_coord_space(o2w, isect.n);
  Matrix3x3 w2o = o2w.T();

  // w_out points towards the source of the ray (e.g.,
  // toward the camera if this is a primary ray)
  const Vector3D hit_p = r.o + r.d * isect.t;
  const Vector3D w_out = w2o * (-r.d);
  Vector3D L_out;

  for (auto light : scene->lights) {
		if (light->is_delta_light()) { // sample point light once
			Vector3D wi;
			double pdf, distToLight;
			Vector3D s = light->sample_L(hit_p, &wi, &distToLight, &pdf);
			Vector3D w_in = w2o * wi;
			if (w_in.z < 0) continue; // only hemisphere

			Ray r_(hit_p + EPS_D * wi, wi, distToLight);
			Intersection is;
			if (bvh->intersect(r_, &is)) continue; // blocked

			L_out += s * isect.bsdf->f(w_out, w_in) * w_in.z / pdf;
		}
		else {
			for (int i = 0; i < ns_area_light; ++i) {
				Vector3D wi;
				double pdf, distToLight;
				Vector3D s = light->sample_L(hit_p, &wi, &distToLight, &pdf);
				Vector3D w_in = w2o * wi;
				if (w_in.z < 0) continue;

				Ray r_(hit_p + EPS_D * wi, wi, distToLight);
				Intersection is;
				if (bvh->intersect(r_, &is)) continue;

				L_out += s * isect.bsdf->f(w_out, w_in) * w_in.z / pdf;
			}
      L_out /= ns_area_light;
		}
	}
	
	return L_out;

}

Vector3D PathTracer::zero_bounce_radiance(const Ray &r,
                                          const Intersection &isect) {
  // TODO: Part 3, Task 2
  // Returns the light that results from no bounces of light

  return isect.bsdf->get_emission();
}

Vector3D PathTracer::one_bounce_radiance(const Ray &r,
                                         const Intersection &isect) {
  // TODO: Part 3, Task 3
  // Returns either the direct illumination by hemisphere or importance sampling
  // depending on `direct_hemisphere_sample`

  return (direct_hemisphere_sample) ?
    estimate_direct_lighting_hemisphere(r, isect) : estimate_direct_lighting_importance(r, isect);
}

Vector3D PathTracer::at_least_one_bounce_radiance(const Ray &r,
                                                  const Intersection &isect) {
  Matrix3x3 o2w;
  make_coord_space(o2w, isect.n);
  Matrix3x3 w2o = o2w.T();

  Vector3D hit_p = r.o + r.d * isect.t;
  Vector3D w_out = w2o * (-r.d);

  Vector3D L_out;
  if (!isect.bsdf->is_delta()) {
    L_out += one_bounce_radiance(r, isect);
  }

  // TODO: Part 4, Task 2
  // Returns the one bounce radiance + radiance from extra bounces at this point.
  // Should be called recursively to simulate extra bounces.

  Vector3D w_in;
	double pdf;
	Vector3D bsdf = isect.bsdf->sample_f(w_out, &w_in, &pdf);
	Vector3D wi = o2w * w_in;
	if (pdf == 0) return L_out;

  double p_RR = 0.3; // Russian Roulette probability
  if (max_ray_depth > 1 && r.depth == 0) {
    p_RR = 0;
  }
  if (coin_flip(1.0 - p_RR) && r.depth < max_ray_depth - 1) {
    Ray r_(hit_p + EPS_D * wi, wi, r.max_t, r.depth + 1);
    Intersection is;
    if (bvh->intersect(r_, &is)) {
      Vector3D L = at_least_one_bounce_radiance(r_, is);
      if(isect.bsdf->is_delta()) 
        L += zero_bounce_radiance(r_, is);
      
      L_out += L * bsdf * w_in.z / pdf / (1.0 - p_RR);
    }
  }

  return L_out;
}

Vector3D PathTracer::est_radiance_global_illumination(const Ray &r) {
  Intersection isect;
  Vector3D L_out;

  // You will extend this in assignment 3-2.
  // If no intersection occurs, we simply return black.
  // This changes if you implement hemispherical lighting for extra credit.

  // The following line of code returns a debug color depending
  // on whether ray intersection with triangles or spheres has
  // been implemented.
  //
  // REMOVE THIS LINE when you are ready to begin Part 3.
  
  if (!bvh->intersect(r, &isect))
    return envLight ? envLight->sample_dir(r) : L_out;

  // NORMAL SHADING
  // return (isect.t == INF_D) ? debug_shading(r.d) : normal_shading(isect.n);

  // TODO (Part 3): Return the direct illumination.

  // TODO (Part 4): Accumulate the "direct" and "indirect"
  // parts of global illumination into L_out rather than just direct


  L_out += zero_bounce_radiance(r, isect);
  if (max_ray_depth > 0) {
    L_out += at_least_one_bounce_radiance(r, isect);
  }
  
  // Vector3D L_direct, L_indirect;
  // L_direct = one_bounce_radiance(r, isect);
  // L_indirect = L_out - L_direct;

  return L_out;
}

void PathTracer::raytrace_pixel(size_t x, size_t y) {
  // TODO (Part 1.2):
  // Make a loop that generates num_samples camera rays and traces them
  // through the scene. Return the average Vector3D.
  // You should call est_radiance_global_illumination in this function.
  
  // TODO (Part 5):
  // Modify your implementation to include adaptive sampling.
  // Use the command line parameters "samplesPerBatch" and "maxTolerance"

  size_t num_samples = ns_aa;          // total samples to evaluate
  Vector2D origin = Vector2D(x, y); // bottom left corner of the pixel

  double s1, s2;

  Vector3D sample_sum;
  size_t count = 0;

  while(count < num_samples){
    count++;
    Vector2D sample = origin + gridSampler->get_sample();
    if(num_samples == 1){
      sample = origin + Vector2D(0.5, 0.5);
    }
    Ray r = camera->generate_ray(sample.x / sampleBuffer.w, sample.y / sampleBuffer.h);
    Vector3D s = est_radiance_global_illumination(r);
    sample_sum += s;

    float illum = s.illum();
    s1 += illum;
    s2 += illum * illum;

    // Check for convergence
    if (count % samplesPerBatch == 0) {
      double mean = s1 / count;
      double var = (1.0 / (count - 1.0)) * (s2 - (s1 * s1 / count));
      if (1.96 * sqrt(var) / sqrt(count) <= maxTolerance * mean) break; // converged
    }
  }

  sampleBuffer.update_pixel(sample_sum/count, x, y);
  sampleCountBuffer[x + y * sampleBuffer.w] = count;

}

void PathTracer::autofocus(Vector2D loc) {
  Ray r = camera->generate_ray(loc.x / sampleBuffer.w, loc.y / sampleBuffer.h);
  Intersection isect;

  bvh->intersect(r, &isect);

  camera->focalDistance = isect.t;
}

} // namespace CGL
