#include "Matrix.h"

// We have to forward declare possible instantiations of the template class here so the compiler knows what's up.
// These two are by far the most likely for computational physics applications.
// If we don't do this then we'll have trouble linking while keeping the function definitions in the .cpp-file.
// An alternative is to keep the definitions in the header.
template
class LA::Matrix<f64>;
template
class LA::Matrix<Complex>;

template<class T>
LA::Matrix<T>::Matrix(u32 rows, u32 cols)
    : m_rows(rows),
      m_cols(cols),
      m_data(rows * cols) {
    // The std::vector<T> class has a constructor explicit vector (size_t count)
    // which is as follows:ma
    // 4) Constructs the container with count default-inserted instances of T. No copies are made.
    // ref: https://en.cppreference.com/w/cpp/container/vector/vector
    std::fill(m_data.begin(), m_data.end(), T());
}

template<class T>
void LA::Matrix<T>::zero() {
    std::fill(m_data.begin(), m_data.end(), T());
}

template<class T>
T &LA::Matrix<T>::operator()(u32 i, u32 j) {
    u32 accessIndex = i * m_cols + j;
    ASSERT(accessIndex < m_data.size());
    return m_data[accessIndex];
}

template<class T>
T LA::Matrix<T>::operator()(u32 i, u32 j) const {
    u32 accessIndex = i * m_cols + j;
    ASSERT(accessIndex < m_data.size());
    return m_data[accessIndex];
}

template<class T>
T* LA::Matrix<T>::GetDataPtr() {
    return m_data.data();
}

template<class T>
void LA::Matrix<T>::copyToVector(std::vector<T>& result) const {
    ASSERT(result.size() == m_data.size());
    result = m_data;
}

template<class T>
u32 LA::Matrix<T>::totalSize() const {
    return m_rows*m_cols;
}

template<class T>
u32 LA::Matrix<T>::numRows() const {
    return m_rows;
}

template<class T>
u32 LA::Matrix<T>::numCols() const {
    return m_cols;
}