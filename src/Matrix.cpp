#include "Matrix.h"

#include "custom_asserts.h"

template<class T>
T &LA::Matrix<T>::operator()(u32 i, u32 j) {
    u32 accessIndex = i * m_cols + j;
    ASSERT(accessIndex < m_data.size());
    return m_data
}

template<class T>
T LA::Matrix<T>::operator()(u32 i, u32 j) const {
    u32 accessIndex = i * m_cols + j;
    ASSERT(accessIndex < m_data.size());
    return m_data
}