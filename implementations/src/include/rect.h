#pragma once

#include <type_traits>
#include <stdexcept>
#include <cmath>

// Put in a namespace because 'rect' is such an exceedingly common name
// that collisions may easily occur
namespace custom {

template <typename T>
class rect
{
public:
    static_assert(std::is_arithmetic<T>::value, "type T must be an arithmetic type");

    rect(T x, T y, T width, T height);

    T x() const { return _x; }
    T y() const { return _y; }
    T width() const { return _w; }
    T height() const { return _h; }

    double diagonal() const;

    static rect<T> from_corners(T x1, T x2, T y1, T y2);

private:
    T _x;
    T _y;
    T _w;
    T _h;
};

template <typename T>
inline rect<T>::rect(T x, T y, T width, T height)
    :   _x(x), _y(y), _w(width), _h(height)
{
    if (_w < 0 || _h < 0)
        throw std::logic_error("rect width and height must be non-negative.");
}

template <typename T>
inline double rect<T>::diagonal() const
{
    return sqrt(width() * width() + height() * height());
}

template <typename T>
inline rect<T> rect<T>::from_corners(T x1, T x2, T y1, T y2)
{
    return rect(std::min(x1, x2),
                std::min(y1, y2),
                std::abs(x2 - x1),
                std::abs(y2 - y1));
}

}
