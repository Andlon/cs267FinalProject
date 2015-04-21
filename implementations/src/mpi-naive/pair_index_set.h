#pragma once

#include <cstdint>
#include <cmath>
#include <utility>

typedef uint_least64_t index_type;
typedef std::pair<index_type, index_type> index_pair;

/**
 * @brief The pair_index_set class represents an index set for a set of
 * pairs (i, j) for i, j = 1, ..., n with j > i.
 *
 * In other words, it is an index set for the set of all combinations of integers
 * i, j with i != j and i, j = 1, ..., n
 *
 * The index set itself is given by { 1, ..., count }.
 *
 */
class pair_index_set {
public:
    explicit pair_index_set(index_type n)
        : _n(n),
          _n2(n * n),
          _count( (n - 1) * n / 2)
    { }

    /**
     * @brief map_to_pair Maps the given index into the set to the corresponding (i, j) pair.
     * If k is not in the range [1, count()], the result is undefined.
     *
     * @param k The index into the set.
     * @return The corresponding (i, j) pair.
     */
    index_pair map_to_pair(index_type k) const;

    /**
     * @brief count Returns the number of pairs in the index set.
     * @return The number of pairs in the index set.
     */
    index_type count() const { return _count; }

private:
    index_type _n;
    index_type _n2;
    index_type _count;
};

inline index_pair pair_index_set::map_to_pair(index_type k) const
{
    // Note: need to work more on proving/testing this formula. Since
    // we're dealing with enormous numbers the cast to double might lose
    // integer precision and we might end up with the wrong i.
    // Further testing required.

    double a = round(_n) + 0.5;
    double b = round(_n2 - _n) + (9.0 / 4.0) - round(2 * k);

    index_type i = floor(a - sqrt(b));
    index_type j = k + i + (i - 1) * i / 2 - _n * (i - 1);
    return std::make_pair(i, j);
}
