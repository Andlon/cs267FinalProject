#pragma once

#include <cstdlib>
#include <stdexcept>

namespace pev {

template <typename T>
class interval
{
public:
    interval(T a, T b)
        : _a(a), _b(b)
    {
        if (b < a)
            throw std::logic_error("b must be greater than or equal to a");
    }

    T a() const { return _a; }
    T b() const { return _b; }

private:
    T _a;
    T _b;
};

/**
 * @brief The uniform_distribution class represents the semi-uniform distribution
 * of n objects into p buckets, with the last bucket potentially containing
 * fewer objects than the other buckets, which hold an equal number of objects.
 *
 * The objects are labelled { 0, 1, ..., n - 1 } and the buckets are numbered
 * { 0, 1, ..., p - 1 }, in line with usual C++ indexing.
 *
 * The class is a useful tool when dealing with resource allocation as it
 * abstracts away the details of dealing with error-prone integer divisions.
 */
class uniform_distribution
{
public:
    uniform_distribution(size_t object_count, size_t bucket_count);

    size_t object_count() const { return _object_count; }
    size_t bucket_count() const { return _bucket_count; }

    size_t objects_in_bucket(size_t bucket) const;
    interval<size_t> interval_for_bucket(size_t bucket) const;

private:
    size_t _object_count;
    size_t _bucket_count;
    size_t _divider;
    size_t _lower;
    size_t _upper;
};

inline uniform_distribution::uniform_distribution(size_t object_count, size_t bucket_count)
    :   _object_count(object_count), _bucket_count(bucket_count)
{
    if (bucket_count == 0u)
        throw std::logic_error("bucket_count must be non-zero");

    _lower = object_count / bucket_count;
    _upper = (object_count + bucket_count - 1) / bucket_count;

    // Buckets with index < _divider have _upper items in their buckets,
    // while buckets with index >= _divider have _lower items in their buckets
    if (_lower == _upper)
        _divider = object_count;
    else
        _divider = object_count - bucket_count * _lower;
}

inline size_t uniform_distribution::objects_in_bucket(size_t bucket) const
{
    if (bucket >= bucket_count())
        throw std::logic_error("bucket must not be larger than or equal to bucket_count");

    return bucket < _divider ? _upper : _lower;
}

inline interval<size_t> uniform_distribution::interval_for_bucket(size_t bucket) const
{
    if (bucket >= bucket_count())
        throw std::logic_error("bucket must not be larger than or equal to bucket_count");

    size_t a, b;
    if (bucket < _divider)
    {
        a = bucket * _upper;
        b = a + _upper;
    } else
    {
        a = _divider * _upper + (bucket - _divider) * _lower;
        b = a + _lower;
    }

    if (a >= b) {
        std::cout << "(a, b) = (" << a << ", " << b << ") for bucket " << bucket
                   << std::endl;
    }

    return interval<size_t>(a, b);
}

}
