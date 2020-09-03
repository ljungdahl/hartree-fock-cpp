#pragma once

#include <vector>

#include "typedefs.h"
#include "logger.h"
#include "custom_asserts.h"

namespace LA {

    template<class T>
    class Matrix {
    public:

        Matrix();
        Matrix(u32 rows, u32 cols);

        void zero();

        T* GetDataPtr();

        T &operator()(u32 i, u32 j);

        T operator()(u32 i, u32 j) const;

        void copyToVector(std::vector<T>& result) const;

        u32 totalSize() const;
        u32 numRows() const;
        u32 numCols() const;

        void resize(u32 rows, u32 cols);
    private:
        u32 m_rows, m_cols;
        std::vector<T> m_data;

    private:
    };

}


