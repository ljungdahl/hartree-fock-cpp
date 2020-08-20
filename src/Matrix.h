#pragma once
#include <vector>

#include "typedefs.h"
#include "logger.h"

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

        T& operator()(u32 i, u32 j);
        T operator()(u32 i, u32 j) const;

        void setToZero() {
            std::fill(m_data.begin(), m_data.end(), T());
        }

    private:
    u32 m_rows, m_cols;
    std::vector<T> m_data;

    private:
    };

}


