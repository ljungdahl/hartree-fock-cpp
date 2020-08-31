#pragma once
/* NOTE(anton):
 * At the moment this is just a thin template wrapper around the std::vector class.
 * Hopefully this can help in switching underlying vector type later.
 * TODO(anton): Make this a proper mathematical object as well.
 */
#include "logger.h"
#include "typedefs.h"

namespace LA {

template<class T>
class Vector {
public:
    Vector();
    Vector(u32 size);
    Vector(u32 size, T initial_value);

    u32 size() const;
    u32 last_index();
    void resize(u32 size);
    void resize(u32 size, T value);
    void fill(T value);
    void push_back(T value);
    void clear();


    void operator=(T value);
    T& operator[](u32 i);
    const T&operator[](u32 i) const;

    std::vector<T>& data();
    const std::vector<T>& data() const;
    T* dataPtr();
    const T* dataPtr() const;
public:

private:
    u32 m_size;
    std::vector<T> m_data;
private:

};

}
