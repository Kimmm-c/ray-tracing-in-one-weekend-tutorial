/*
 * This ray tracing project follows the tutorial https://raytracing.github.io/books/RayTracingInOneWeekend.html#overview
 * Refer to this source for math concepts such as linear interpolation, ray-sphere intersection calculation, etc..
 */

#include "src/color.h"
#include "src/vec3.h"
#include "src/Ray.h"

#include <iostream>

/**
 * @brief Determine if the ray hits the sphere
 *
 * Notes:
 * 1. The ray equation: P(t) = Q + t*d. Where Q is the origin, d is the direction, t is a variable.
 * 2. The sphere with center C(cx, cy, cz) equation: (x - cx)^2 + (y - cy)^2 + (z - cz)^2 = r^2. With (x, y, z) is a point
 * on the sphere. The distance from this point to the center of the sphere is equal to the radius.
 * 3. For a ray to hit the sphere, there has to exist at least one point (x, y ,z) that satisfies the 2 equations above.
 * As t is the only variable, we have to solve for t. If t exists => the ray hits the sphere. Otherwise, the ray doesn't
 * intersect the sphere.
 *
 * @param center The center of the sphere.
 * @param radius The radius of the sphere.
 * @param r The ray.
 * @return True if the ray intersects the sphere at some points. Otherwise, return False.
 */
bool hit_sphere(const point3 &center, double radius, const ray &r) {
    // Refer to the tutorial link for the math behind this.
    vec3 oc = center - r.origin();
    auto a = dot(r.direction(), r.direction());
    auto b = -2.0 * dot(r.direction(), oc);
    auto c = dot(oc, oc) - radius * radius;
    auto discriminant = b * b - 4 * a * c;
    return (discriminant >= 0);
}

/**
 * @brief Return the color of the given scene ray. We want to create an image of a gradient running from white (at the
 * bottom) to sky blue (at the top).
 *
 * Notes:
 * 1. Why normalizing the ray direction?
 * Normalizing makes ray direction a unit vector (vector with length = 1). This ensures the mapping between 'y' and 'a'
 * works properly, which is critical in achieving the desired color interpolation (blending). Additionally, in ray
 * tracing, we only care about the direction of the ray, not its magnitude (length).
 *
 * 2. Why converting y to a?
 * y represents the ray direction with respective to the viewport. Therefore, it may hold negative values, which we
 * can't use to compute color values, which are always positive.
 *
 * @param r A scene ray.
 * @return The color of the given scene ray.
 */
color ray_color(const ray &r) {
    color startColor(1.0, 1.0, 1.0);    // White at a=0
    color endColor(0.5, 0.7, 1.0);      // Sky blue at a=1

    // Return red for any ray that hits the sphere
    if (hit_sphere(point3(0, 0, -1), 0.5, r)) {
        return color(1, 0, 0); // Pure red
    }


    // Normalize the ray direction
    vec3 unit_direction = unit_vector(r.direction());

    // Convert y interval from [-1, 1] to [0, 1]
    auto a = 0.5 * (unit_direction.y() + 1.0);

    // Calculate the blend color
    color blendColor = (1.0 - a) * startColor + a * endColor;

    return blendColor;
}

int main() {

    // Image

//    auto aspect_ratio = 16.0 / 9.0;
//    int image_width = 400;
//
//    // Calculate the image height, and ensure that it's at least 1.
//    int image_height = int(image_width / aspect_ratio);
//    image_height = (image_height < 1) ? 1 : image_height;
//
//    // Camera
//    auto focal_length = 1.0; // Focal length: Length from the eye point/camera to the center of the viewport orthogonally.
//    auto viewport_height = 2.0;
//    auto viewport_width = viewport_height * (double(image_width) / image_height);
//    auto camera_center = point3(0, 0, 0); // Initialize the eye point at the origin.
//
//    // Calculate the vectors across the horizontal and down the vertical viewport edges.
//    auto viewport_u = vec3(viewport_width, 0, 0); // This vector represents the top edge of the viewport.
//    // Note: We negate the viewport height because the pixel grid's y coordinate is the inverse of the viewport's y coordinate.
//    // Pixel grid: y starts from 0 at the top left and increases towards the bottom.
//    // Viewport: y starts from 0 at the bottom left and increases towards the top.
//    auto viewport_v = vec3(0, -viewport_height, 0); // This vector represents the left edge of the viewport.
//
//    // Calculate the horizontal and vertical delta vectors from pixel to pixel
//    // delta vector: the pixel grid's width OR the distance between 2 pixels.
//    auto pixel_delta_u = viewport_u / image_width;
//    auto pixel_delta_v = viewport_v / image_height;
//
//    // Calculate the location of the upper left pixel
//    // Notes:
//    // 1. The viewport upper left corner's position is NOT the same as the upper left pixel's position.
//    // 2. In 3D graphics, we use right-handed coordinate. This means z gets more negative moving away from the camera.
//    // Therefore, to move towards the viewport, we have to negate the focal length.
//    auto viewport_upper_left = camera_center + vec3(0, 0, -focal_length) - viewport_u / 2 - viewport_v / 2;
//    auto pixel_upper_left = viewport_upper_left + 0.5 * (pixel_delta_u + pixel_delta_v);
//
//    // Render
//    std::cout << "P3\n" << image_width << ' ' << image_height << "\n255\n";
//
//    for (int j = 0; j < image_height; j++) {
//        std::clog << "\rScanlines remaining: " << (image_height - j) << ' ' << std::flush;
//        for (int i = 0; i < image_width; i++) {
//            // Get the location of the current pixel
//            auto pixel_center = pixel_upper_left + (i * pixel_delta_u) + (j * pixel_delta_v);
//            // Get the position of the current pixel relative to the camera
//            auto ray_direction = pixel_center - camera_center;
//            ray r(camera_center, ray_direction);
//
//            color pixel_color = ray_color(r);
//            write_color(std::cout, pixel_color);
//        }
//    }
//
//    std::clog << "\rDone.                 \n";

    point3 c(0, 0, 0);
    vec3 d(-0.25, -0.25, -1);
    ray r(c,d);
    hit_sphere(c, 0.5, r);
}