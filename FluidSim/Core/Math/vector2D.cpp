#include "vector2D.h"
#include <ostream>

std::ostream& operator<<(std::ostream& os, const Vector2D& v) {
	os << "(" << v.x << "," << v.y << ")";
	return os;
}