#pragma once

#include <cstdlib>
#include <stdexcept>

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
    size_t _objects_per_bucket;
    size_t _objects_last_bucket;
};

inline uniform_distribution::uniform_distribution(size_t object_count, size_t bucket_count)
    :   _object_count(object_count), _bucket_count(bucket_count)
{
    if (bucket_count == 0u)
        throw std::logic_error("bucket_count must be non-zero");

    _objects_per_bucket = (object_count + bucket_count - 1) / (bucket_count);
    _objects_last_bucket = object_count - (bucket_count - 1) * _objects_per_bucket;
}

inline size_t uniform_distribution::objects_in_bucket(size_t bucket) const
{
    if (bucket >= bucket_count())
        throw std::logic_error("bucket must not be larger than or equal to bucket_count");

    return std::min(_objects_per_bucket, object_count() - bucket * _objects_per_bucket);
}

inline interval<size_t> uniform_distribution::interval_for_bucket(size_t bucket) const
{
    size_t a = bucket * _objects_per_bucket;
    size_t b = std::min(a + _objects_per_bucket, object_count());
    return interval<size_t>(a, b);
}
