#pragma once

#include <vector>

#include "typedefs.h"
#include "logger.h"
#include "custom_asserts.h"

namespace LA {

    template<class T>
    class Matrix {
    public:

        Matrix(u32 rows, u32 cols)
                : m_rows(rows),
                  m_cols(cols),
                  m_data(rows * cols) {
            // The std::vector<T> class has a constructor explicit vector (size_t count)
            // which is as follows:ma
            // 4) Constructs the container with count default-inserted instances of T. No copies are made.
            // ref: https://en.cppreference.com/w/cpp/container/vector/vector
            std::fill(m_data.begin(), m_data.end(), T());
        }
//
//        T& operator()(u32 i, u32 j);
//        T operator()(u32 i, u32 j) const;

        void setToZero() {
            std::fill(m_data.begin(), m_data.end(), T());
        }


        T &operator()(u32 i, u32 j) {
            u32 accessIndex = i * m_cols + j;
            ASSERT(accessIndex < m_data.size());
            return m_data[accessIndex];
        }

        T operator()(u32 i, u32 j) const {
            u32 accessIndex = i * m_cols + j;
            ASSERT(accessIndex < m_data.size());
            return m_data[accessIndex];
        }

        void set(u32 i, u32 j, T value) {
            u32 accessIndex = i * m_cols + j;
            ASSERT(accessIndex < m_data.size());
            m_data[accessIndex] = value;
        }

        T *GetDataPtr() {
            return m_data.data();
        }

        void copyToVector(std::vector<T>& result) const {
            ASSERT(result.size() == m_data.size());
            result = m_data;
        }

        u32 TotalSize() const {
            return m_rows*m_cols;
        }

        u32 numRows() const {
            return m_rows;
        }

        u32 numCols() const {
            return m_cols;
        }

    private:
        u32 m_rows, m_cols;
        std::vector<T> m_data;

    private:
    };

}


