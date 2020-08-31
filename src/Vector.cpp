#include "custom_asserts.h"
#include "Vector.h"

// We have to forward declare possible instantiations of the template class here so the compiler knows what's up.
// These two are by far the most likely for computational physics applications.
// If we don't do this then we'll have trouble linking while keeping the function definitions in the .cpp-file.
// An alternative is to keep the definitions in the header.
template
class LA::Vector<f64>;
template
class LA::Vector<Complex>;

// Definitions
template<class T>
LA::Vector<T>::Vector() : m_size(0) {
    m_data.clear();
}

template<class T>
LA::Vector<T>::Vector(u32 size, T initial_value) :
    m_size(size) {

    m_data.resize(size, initial_value);
}

template<class T>
u32 LA::Vector<T>::size() const {

    return m_data.size();
}

template<class T>
void LA::Vector<T>::fill(T value) {
    ASSERT(m_data.size() > 0);
    std::fill(m_data.begin(), m_data.end(), value);
}

template<class T>
void LA::Vector<T>::resize(u32 size) {
    if (m_data.size() > 0) {
        m_data.clear();
    }
    m_data.resize(size);
    m_size = m_data.size();
}

template<class T>
void LA::Vector<T>::resize(u32 size, T value) {
    m_data.clear();
    m_data.resize(size, value);
    m_size = m_data.size();
}
template<class T>
void LA::Vector<T>::operator=(T value) {
    fill(value);
}
template<class T>
T &LA::Vector<T>::operator[](u32 i) {
    ASSERT( i < m_data.size());
    return m_data[i];
}

template<class T>
const T &LA::Vector<T>::operator[](u32 i) const {
    ASSERT( i < m_data.size());
    return m_data[i];
}

template<class T>
LA::Vector<T>::Vector(u32 size) :
    m_size(size) {

    m_data.resize(size, T());
}
template<class T>
u32 LA::Vector<T>::last_index() {
    if (m_data.size() > 0) {
        return m_data.size() - 1;
    } else {
        return 0;
    }
}
template<class T>
void LA::Vector<T>::push_back(T value) {
    m_data.push_back(value);
}

template<class T>
std::vector<T> &LA::Vector<T>::data() {
    return m_data;
}

template<class T>
T *LA::Vector<T>::dataPtr() {
    return m_data.data();
}

template<class T>
const T *LA::Vector<T>::dataPtr() const {
    return m_data.data();
}


template<class T>
const std::vector<T> &LA::Vector<T>::data() const {
    return m_data;
}
template<class T>
void LA::Vector<T>::clear() {
m_data.clear();
}
