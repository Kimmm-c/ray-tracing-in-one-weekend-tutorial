#ifndef UTILS_H
#define UTILS_H

#include <cmath>
#include <iostream>
#include <limits>
#include <memory>

// C++ Std Usings

// make_shared: allocate a new instance of the specified type and return a shared_ptr.
// Eg: shared_ptr<vec3>   vec3_ptr   = make_shared<vec3>(1.414214, 2.718281, 1.618034);
using std::make_shared;
// shared_ptr: reference count increases by 1 everytime the value is referenced to by another shared_ptr.
// The object is safely deleted once the count reaches 0.
using std::shared_ptr;

// Constants

const double infinity = std::numeric_limits<double>::infinity();
const double pi = 3.1415926535897932385;

// Utility Functions

inline double degrees_to_radians(double degrees) {
    return degrees * pi / 180.0;
}

// Common Headers

#include "color.h"
#include "ray.h"
#include "vec3.h"

#endif