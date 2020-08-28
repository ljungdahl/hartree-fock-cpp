#include <cmath>
#include "custom_asserts.h"
#include "Bsplines.h"


Atom::Bsplines::Bsplines(u32 numKnotPoints_, u32 bsplineOrder_)
        : m_numKnotPoints(numKnotPoints_),
          m_order(bsplineOrder_),
          m_numBsplines(numKnotPoints_ - bsplineOrder_) {

    m_knotPoints.resize(m_numKnotPoints, Complex(0.0));

    m_bsplvb_dL.resize(m_order, Complex(0.0));
    m_bsplvb_dR.resize(m_order, Complex(0.0));
    m_Sp.resize(m_order, Complex(0.0));
}

void Atom::Bsplines::setupKnotPoints(const std::vector<Complex> &gridPoints, Atom::knotSequenceType sequenceType) {
    u32 numberOfGhostPointsInEachEnd = (m_order - 1);
    u32 numberOfPhysicalKnotPoints = m_numKnotPoints - 2 * numberOfGhostPointsInEachEnd;

    // For proper indexing
    u32 offset = numberOfGhostPointsInEachEnd;


    if (sequenceType == Atom::knotSequenceType::Linear) {
        // Linear knot sequence
        {
            // Linear stride
            u32 stride = (u32) (gridPoints.size() / (numberOfPhysicalKnotPoints - 1));

            // set ghost knot points
            for (u32 i = 0; i < numberOfGhostPointsInEachEnd; ++i) {
                // Set first k points to first grid point value.
                m_knotPoints[i] = gridPoints[0];

                // Set last k points to last grid point value.
                u32 backwardsCountIndex = (m_knotPoints.size() - i) - 1;
                m_knotPoints[backwardsCountIndex] = gridPoints[gridPoints.size() - 1];
            }

            // Set starting physical knotpoint
            m_knotPoints[offset] = gridPoints[0];

            // Set physical knotpoints except start and end
            for (u32 i = 1; i < numberOfPhysicalKnotPoints - 1; ++i) {
                u32 knotPointIndex = i + offset;
                u32 gridIndex = i * stride;
                ASSERT(gridIndex < gridPoints.size() - 1);
                Complex knotPointValue = gridPoints[gridIndex];
                m_knotPoints[knotPointIndex] = knotPointValue;
            }

            // Set physical endpoint
            m_knotPoints[m_knotPoints.size() - offset - 1] = gridPoints[gridPoints.size() - 1];
        }
    }

    if (sequenceType == Atom::knotSequenceType::firstPointsCloserThenLinear) {
        {
            // set ghost knot points
            for (u32 i = 0; i < numberOfGhostPointsInEachEnd; ++i) {
                // Set first k points to first grid point value.
                m_knotPoints[i] = gridPoints[0];

                // Set last k points to last grid point value.
                u32 backwardsCountIndex = (m_knotPoints.size() - i) - 1;
                m_knotPoints[backwardsCountIndex] = gridPoints[gridPoints.size() - 1];
            }

            // Set starting physical knotpoint
            m_knotPoints[offset] = gridPoints[0];


            // Set tighter cluster close to zero
            u32 clusterEndIndex = 12;
            auto clusterPointDiff = gridPoints[5] - gridPoints[0] / 3.0;
            auto lastKnotPoint = Complex(0.0);
            for (u32 i = 1; i <= clusterEndIndex; i++) {
                u32 knotPointIndex = i + offset;
                Complex knotPointValue = gridPoints[0] + clusterPointDiff * ((f64) i);
                m_knotPoints[knotPointIndex] = knotPointValue;
                lastKnotPoint = knotPointValue;
            }

            // find index for lastKnotPoint on grid:
            u32 gridPointIndex = 0;
            for (auto point : gridPoints) {
                if (point.real() > lastKnotPoint.real()) {
                    break;
                }
                gridPointIndex += 1;
            }
//            Logger::Trace("gridPointIndex for lastKnotPoint: %i", gridPointIndex);

            u32 numberOfPhysicalPointsAfterCluster = numberOfPhysicalKnotPoints - clusterEndIndex;
            u32 stride = (u32) (gridPoints.size() - gridPointIndex + 1) / (numberOfPhysicalPointsAfterCluster - 1);
            // Set physical knotpoints between cluster and end point
            for (u32 i = 0; i < numberOfPhysicalPointsAfterCluster; ++i) {
                u32 knotPointIndex = i + offset + clusterEndIndex;
                u32 gridIndex = gridPointIndex + i * stride;
//                Logger::Trace("gridIndex %i", gridIndex);
                ASSERT(gridIndex < gridPoints.size() - 1);
                Complex knotPointValue = gridPoints[gridIndex];
                m_knotPoints[knotPointIndex] = knotPointValue;
            }

            // Set physical endpoint
            m_knotPoints[m_knotPoints.size() - offset - 1] = gridPoints[gridPoints.size() - 1];
        }
    }

}

void Atom::Bsplines::DebugLogKnotSequence() {
    u32 index = 0;
    for (auto point : m_knotPoints) {
        Logger::Log("Knotpoint %i: (%E, %E)", index, point.real(), point.imag());
        index += 1;
    }
}

u32 Atom::Bsplines::numberOfBsplines() {
    return m_numKnotPoints - m_order;
}

void Atom::Bsplines::bsplvb_Complex(Complex coordinate, u32 left, u32 k_order) {
/*
 * NOTE(anton):
 * This is a C++ version of the de Boor Fortran routine bsplvb, for Complex valued Bsplines.
 * If the bspline order is k, it calculates the k contributing (ie non-zero bsplines) at the
 * specified coordinate. For simplicity (and since we rarely deal with very large k), we always
 * calculate all necessary orders in the recursion (amounting to assuming index = 1 below).
 *
 * We also keep the necessary arrays as private member variables so we can allocate them once
 * at initialisation, since the bspline order is constant for each instance of the Bsplines class.
 * TODO(anton): Test if this is an impactful optimisation.
 *
 * TODO(anton): More detailed comments.
 * For now refer to the de Boor-routine bsplvb.f supplied commented at the end of this file.
 */
    Complex saved;
    Complex term;
    Complex x = coordinate;

    // Reset arrays
    std::fill(m_bsplvb_dR.begin(), m_bsplvb_dR.end(), Complex(0.0));
    std::fill(m_bsplvb_dL.begin(), m_bsplvb_dL.end(), Complex(0.0));
    std::fill(m_Sp.begin(), m_Sp.end(), Complex(0.0));

    // Assume index == 1, ie starting order is j = 0;
    m_Sp[0] = Complex(1.0, 0.0);

    u32 k = k_order;
    auto &t = m_knotPoints;
    auto &dR = m_bsplvb_dR;
    auto &dL = m_bsplvb_dL;

    for (u32 j = 0; j < k - 1; j += 1) {

        dR[j] = t[left + j + 1] - x;
        dL[j] = x - t[left - j];
        saved = Complex(0.0);

        for (u32 i = 0; i <= j; i++) {
            term = m_Sp[i] / (dR[i] + dL[j - i]);
            m_Sp[i] = saved + dR[i] * term;
            saved = dL[j - i] * term;
        }

        m_Sp[j + 1] = saved;

    }

}

Complex Atom::Bsplines::GetBsplineAtCoordinate(Complex coordinate, u32 bsplineIndex) {
    //TODO(anton): How to get proper first Bspline? How to think about it?
    /*
     * NOTE(anton):
     * Based on Fortran routine bget_cmplx. We get the bsplineIndex:th Bspline at x = coordinate.
     */

    auto x_real = coordinate.real();
    auto knotPointsStart_real = m_knotPoints[0].real();
    auto knotPointsEnd_real = m_knotPoints[m_knotPoints.size() - 1].real();

    bool isCoordinateOutsideGrid = (x_real > knotPointsEnd_real || x_real < knotPointsStart_real);

    ASSERT(!isCoordinateOutsideGrid);

    // NOTE(anton):
    // If we're close to the last grid point (or at the last grid point),
    // and we're looking for the value of the last Bspline at that point,
    // we know that the value should be 1. The indexing for this last point
    // gets messy in bsplvb_Complex, so we just force it to one here.
    // In the collocation application for our usual boundary conditions this Bspline has
    // a zero coefficient anyway.
    if (bsplineIndex == m_numBsplines - 1) {
        f64 tolerance = 1e-8;
        if (std::abs(coordinate - m_knotPoints[m_knotPoints.size() - 1]) < tolerance) {
            return Complex(1.0, 0.0);
        }
    }


    // Find the left knot point index left_knotPoint_index such that
    // m_knotPoints[left_knotPoint_index].real() < coordinate.real() m_knotPoints[left_knotPoint_index+1].real()
    u32 left_knotPoint_index = 0;
    for (int t = 0; t < m_numKnotPoints; t++) {
        if (x_real >= m_knotPoints[t].real()) {
            left_knotPoint_index = t;
        }
    }

    // The index used to get the bspline from bsplvb_Complex().
    u32 accessIndex = bsplineIndex - left_knotPoint_index + m_order - 1;
    if (accessIndex < 0 || accessIndex > m_Sp.size() - 1) {
        // We're not in a bspline that is nonzero on this coordinate. So return zero.
        return Complex(0.0);
    }

    // Calculate non-zero Bsplines on the coordinate.
    bsplvb_Complex(coordinate, left_knotPoint_index, m_order);

    // Get the correct Bspline from the k = bsplineOrder non-zero ones.
    Complex return_value = m_Sp[accessIndex];

    return return_value;
}

Complex Atom::Bsplines::GetBsplineAtCoordinate(Complex coordinate, u32 bsplineIndex, u32 order) {
    //TODO(anton): How to get proper first Bspline? How to think about it?
    /*
     * NOTE(anton):
     * Based on Fortran routine bget_cmplx. We get the bsplineIndex:th Bspline at x = coordinate.
     */

    auto x_real = coordinate.real();
    auto knotPointsStart_real = m_knotPoints[0].real();
    auto knotPointsEnd_real = m_knotPoints[m_knotPoints.size() - 1].real();

    bool isCoordinateOutsideGrid = (x_real > knotPointsEnd_real || x_real < knotPointsStart_real);

    ASSERT(!isCoordinateOutsideGrid);

    // NOTE(anton):
    // If we're close to the last grid point (or at the last grid point),
    // and we're looking for the value of the last Bspline at that point,
    // we know that the value should be 1. The indexing for this last point
    // gets messy in bsplvb_Complex, so we just force it to one here.
    // In the collocation application for our usual boundary conditions this Bspline has
    // a zero coefficient anyway.
    if (bsplineIndex == m_numBsplines - 1) {
        f64 tolerance = 1e-8;
        if (std::abs(coordinate - m_knotPoints[m_knotPoints.size() - 1]) < tolerance) {
            return Complex(1.0, 0.0);
        }
    }


    // Find the left knot point index left_knotPoint_index such that
    // m_knotPoints[left_knotPoint_index].real() < coordinate.real() m_knotPoints[left_knotPoint_index+1].real()
    u32 left_knotPoint_index = 0;
    for (int t = 0; t < m_numKnotPoints; t++) {
        if (x_real >= m_knotPoints[t].real()) {
            left_knotPoint_index = t;
        }
    }

    // The index used to get the bspline from bsplvb_Complex().
    u32 accessIndex = bsplineIndex - left_knotPoint_index + order - 1;
    if (accessIndex < 0 || accessIndex > m_Sp.size() - 1) {
        // We're not in a bspline that is nonzero on this coordinate. So return zero.
        return Complex(0.0);
    }

    // Calculate non-zero Bsplines on the coordinate.
    bsplvb_Complex(coordinate, left_knotPoint_index, order);

    // Get the correct Bspline from the k = bsplineOrder non-zero ones.
    Complex return_value = m_Sp[accessIndex];

    return return_value;
}


Complex Atom::Bsplines::GetDerivativeAtCoordinate(Complex coordinate, u32 bsplineIndex) {
    Complex x = coordinate;
    u32 i = bsplineIndex;
    u32 k = m_order;
    auto dB = GetDerivativeAtCoordinate(x, i, k);
//    auto B_i_k_min_1 = GetBsplineAtCoordinate(x, i, k - 1);
//    auto B_i_plus_1_k_min_1 = GetBsplineAtCoordinate(x, i + 1, k - 1);
//    Complex term1 = B_i_k_min_1 / (m_knotPoints[i + k - 1] - m_knotPoints[i]);
//    Complex term2 = B_i_plus_1_k_min_1 / (m_knotPoints[i + k] - m_knotPoints[i + 1]);
//    Complex dB = Complex((f64) (k - 1)) * (term1 - term2);

    return dB;
}

Complex Atom::Bsplines::GetDerivativeAtCoordinate(Complex coordinate, u32 bsplineIndex, u32 order) {
    Complex x = coordinate;
    u32 i = bsplineIndex;
    u32 k = order;
    Complex dB = Complex(0.0);

    // Check if we're on the grid
    if (std::abs(x) > std::abs(m_knotPoints[m_numKnotPoints - 1])) {
        return dB;
    }
    if (std::abs(x) < std::abs(m_knotPoints[0])) {
        return dB;
    }

    // Check if we're in the last knotpoint.
    if (std::abs(x - m_knotPoints[m_numKnotPoints - 1]) < 1e-8) {
        if (i < m_numBsplines - 2) { // Only the last and second to last Bspline is nonzero in the last knotpoint.
            return dB;
        }
        if (i == m_numBsplines - 1) { // If we are in the last Bspline in the last knotpoint the derivative will be
            dB = Complex((f64) (k - 1)) * Complex(1.0) /
                 (m_knotPoints[m_numKnotPoints - 1] - m_knotPoints[m_numKnotPoints - 1 - k]);
            return dB;
        }
        if (i == m_numBsplines - 2) { // If we are in the second to last Bspline we get this derivative
            dB = Complex((f64) (k - 1)) * Complex(-1.0) /
                 (m_knotPoints[m_numKnotPoints - 1] - m_knotPoints[m_numKnotPoints - 1 - k]);
            return dB;
        }
    }

    u32 left_knotPoint_index = 0;
    f64 x_real = x.real();
    for (int t = 0; t < m_numKnotPoints; t++) {
        if (x_real >= m_knotPoints[t].real()) {
            left_knotPoint_index = t;
        }
    }

    u32 accessIndex = i - left_knotPoint_index + k - 1;
    if (accessIndex < 0 || accessIndex > k - 1) {
        // We're not in a bspline that is nonzero on this coordinate. So return zero.
        return Complex(0.0);
    }

    bsplvb_Complex(x, left_knotPoint_index, k - 1);

    dB = Complex(0.0);
    if (accessIndex == 0) {
        Complex term1 = -m_Sp[accessIndex] / (m_knotPoints[i + k] - m_knotPoints[i + 1]);
        dB = Complex((f64) (k - 1)) * term1;
    } else if (accessIndex == k - 1) {
        Complex term1 = m_Sp[accessIndex - 1] / (m_knotPoints[i + k - 1] - m_knotPoints[i]);
        dB = Complex((f64) (k - 1)) * term1;
    } else {
        Complex term1 = m_Sp[accessIndex - 1] / (m_knotPoints[i + k - 1] - m_knotPoints[i]);
        Complex term2 = m_Sp[accessIndex] / (m_knotPoints[i + k] - m_knotPoints[i + 1]);
        dB = Complex((f64) (k - 1)) * (term1 - term2);
    }
//    auto B_i_k_min_1 = GetBsplineAtCoordinate(x, i, k - 1);
//    auto B_i_plus_1_k_min_1 = GetBsplineAtCoordinate(x, i + 1, k - 1);
//    Complex term1 = B_i_k_min_1 / (m_knotPoints[i + k - 1] - m_knotPoints[i]);
//    Complex term2 = B_i_plus_1_k_min_1 / (m_knotPoints[i + k] - m_knotPoints[i + 1]);
//    dB = Complex((f64) (k - 1)) * (term1 - term2);

    return dB;
}

Complex Atom::Bsplines::GetSecondDerivativeAtCoordinate(Complex coordinate, u32 bsplineIndex) {
    Complex x = coordinate;
    u32 i = bsplineIndex;
    u32 k = m_order;

    auto dB2 = GetSecondDerivativeAtCoordinate(x, i, k);

    return dB2;
}

Complex Atom::Bsplines::GetSecondDerivativeAtCoordinate(Complex coordinate, u32 bsplineIndex, u32 order) {
    Complex x = coordinate;
    u32 i = bsplineIndex;
    u32 k = order;
    ASSERT(k > 2); // Else we can't have second derivative.

    Complex dB2 = Complex(0.0);

    // Check if we're on the grid
    if (std::abs(x) > std::abs(m_knotPoints[m_numKnotPoints - 1])) {
        return dB2;
    }
    if (std::abs(x) < std::abs(m_knotPoints[0])) {
        return dB2;
    }

    auto prefactor = Complex((f64) (k - 1) * (k - 2));

    u32 np = m_numKnotPoints - 1; // Last knotpoint index, to match with fortran bder2 style.
    // Check if we're in the last knotpoint.
    if (std::abs(x - m_knotPoints[m_numKnotPoints - 1]) < 1e-8) {

        if (i < m_numBsplines - 3) { // Only the last,second and third to last Bspline is nonzero in the last knotpoint.
            return dB2;
        }

        if (i == m_numBsplines - 1) { // If we are in the last Bspline in the last knotpoint the derivative will be
            auto denominator =
                    (m_knotPoints[np - 1] - m_knotPoints[np - k]) * (m_knotPoints[np - 2] - m_knotPoints[np - k]);
            dB2 = prefactor * Complex(1.0) / denominator;
            return dB2;
        }

        if (i == m_numBsplines - 2) { // if we are in the second to last Bspline we get this
            auto denominator1 =
                    (m_knotPoints[np - 2] - m_knotPoints[np - k - 1]) * (m_knotPoints[np - 2] - m_knotPoints[np - k]);
            auto term1 = prefactor * Complex(1.0) / denominator1;
            auto denominator2 =
                    (m_knotPoints[np - 1] - m_knotPoints[np - k]) * (m_knotPoints[np - 2] - m_knotPoints[np - k]);
            auto term2 = prefactor * Complex(1.0) / denominator2;
            dB2 = -term1 - term2;
            return dB2;
        }

        if (i == m_numBsplines - 3) { // If we are in the third to last Bspline we get this derivative
            auto denominator =
                    (m_knotPoints[np - 1] - m_knotPoints[np - k - 1]) * (m_knotPoints[np - 2] - m_knotPoints[np - k]);
            dB2 = prefactor * Complex(1.0) / denominator;
            return dB2;
        }


    }

    // Find the left knot point index left_knotPoint_index such that
    // m_knotPoints[left_knotPoint_index].real() < coordinate.real() m_knotPoints[left_knotPoint_index+1].real()
    u32 left_knotPoint_index = 0;
    for (int t = 0; t < m_numKnotPoints; t++) {
        if (x.real() >= m_knotPoints[t].real()) {
            left_knotPoint_index = t;
        }
    }

    u32 accessIndex = i - left_knotPoint_index + k - 1;
    if (accessIndex < 0 || accessIndex > k - 1) {
        // We're not in a bspline that is nonzero on this coordinate. So return zero.
        return Complex(0.0);
    }

    bsplvb_Complex(x, left_knotPoint_index, k - 2);

    dB2 = Complex(0.0);
    if (accessIndex > 1) {
        auto denominator = (m_knotPoints[i + k - 2] - m_knotPoints[i]) * (m_knotPoints[i + k - 1] - m_knotPoints[i]);
        dB2 += prefactor * m_Sp[accessIndex - 2] / denominator;
    }
    if (accessIndex > 0 && accessIndex < k - 1) {
        auto denominator1 =
                (m_knotPoints[i + k - 1] - m_knotPoints[i + 1]) * (m_knotPoints[i + k - 1] - m_knotPoints[i]);
        auto term1 = prefactor * m_Sp[accessIndex - 1] / denominator1;
        auto denominator2 =
                (m_knotPoints[i + k - 1] - m_knotPoints[i + 1]) * (m_knotPoints[i + k] - m_knotPoints[i + 1]);
        auto term2 = prefactor * m_Sp[accessIndex - 1] / denominator2;
        dB2 += - term1 - term2;
    }
    if (accessIndex < k - 2) {
        auto denominator = (m_knotPoints[i + k] - m_knotPoints[i + 2]) * (m_knotPoints[i + k] - m_knotPoints[i + 1]);
        dB2 += prefactor * m_Sp[accessIndex] / denominator;
    }
//    auto B_i_k_min_1 = GetDerivativeAtCoordinate(x, i, k - 1);
//    auto B_i_plus_1_k_min_1 = GetDerivativeAtCoordinate(x, i + 1, k - 1);
//    Complex term1 = B_i_k_min_1 / (m_knotPoints[i + k - 1] - m_knotPoints[i]);
//    Complex term2 = B_i_plus_1_k_min_1 / (m_knotPoints[i + k] - m_knotPoints[i + 1]);
//    Complex dB = Complex((f64) (k - 1)) * (term1 - term2);

    return dB2;
//    Complex x = coordinate;
//    u32 i = bsplineIndex;
//    u32 k = m_order;
//    auto B_i_k_min_2 = GetBsplineAtCoordinate(x, i, k - 2);
//    auto B_i_plus_1_k_min_2 = GetBsplineAtCoordinate(x, i + 1, k - 2);
//    auto B_i_plus_2_k_min_2 = GetBsplineAtCoordinate(x, i + 2, k - 2);
//
//    f64 kk = (f64) k;
//    if(i > 2) {
//        Complex term1 =
//                B_i_k_min_2 /
//                ((m_knotPoints[i + k - 1] - m_knotPoints[i]) * (m_knotPoints[i + k - 2] - m_knotPoints[i]));
//    }
//    Complex term2 = B_i_plus_1_k_min_2 /
//                    ((m_knotPoints[i + k - 1] - m_knotPoints[i]) * (m_knotPoints[i + k - 1] - m_knotPoints[i + 1]));
//    Complex term3 = B_i_plus_1_k_min_2 /
//                    ((m_knotPoints[i + k] - m_knotPoints[i + 1]) * (m_knotPoints[i + k - 1] - m_knotPoints[i + 1]));
//    Complex term4 = B_i_plus_2_k_min_2 /
//                    ((m_knotPoints[i + k] - m_knotPoints[i + 1]) * (m_knotPoints[i + k] - m_knotPoints[i + 2]));
//    Complex dB2 = (kk - 1) * (kk - 2) * (term1 - term2 - term3 + term4);
//
//    return dB2;
}

void Atom::Bsplines::SetBoundaryConditionBsplineIndices(std::vector<u32> &indices) {
    for (auto index : indices) {
        m_usedBsplineIndices.push_back(index);
    }
}


// de Boor's bsplvb for reference.
/*

    subroutine bsplvb ( t, jhigh, index, x, left, biatx )
    !  from  * a practical guide to splines *  by c. de boor
    !calculates the value of all possibly nonzero b-splines at  x  of order
    !
    !               jout  =  max( jhigh , (j+1)*(index-1) )
    !c
    !  with knot sequence  t .
    !c
    !c******  i n p u t  ******
    !  t.....knot sequence, of length  left + jout  , assumed to be nonde-
    !        creasing.  a s s u m p t i o n . . . .
    !                       t(left)  .lt.  t(left + 1)   .
    !   d i v i s i o n  b y  z e r o  will result if  t(left) = t(left+1)
    !  jhigh,
    !  index.....integers which determine the order  jout = max(jhigh,
                                                                !        (j+1)*(index-1))  of the b-splines whose values at  x  are to
    !        be returned.  index  is used to avoid recalculations when seve-
    !        ral columns of the triangular array of b-spline values are nee-
    !        ded (e.g., in  bsplpp  or in  bsplvd ). precisely,
    !                     if  index = 1 ,
    !        the calculation starts from scratch and the entire triangular
    !        array of b-spline values of orders 1,2,...,jhigh  is generated
    !        order by order , i.e., column by column .
    !                     if  index = 2 ,
    !        only the b-spline values of order  j+1, j+2, ..., jout  are ge-
    !        nerated, the assumption being that  biatx , j , deltal , deltar
    !        are, on entry, as they were on exit at the previous call.
    !           in particular, if  jhigh = 0, then  jout = j+1, i.e., just
    !        the next column of b-spline values is generated.
    !c
    !  w a r n i n g . . .  the restriction   jout .le. jmax (= 20)  is im-
    !        posed arbitrarily by the dimension statement for  deltal  and
    !        deltar  below, but is  n o w h e r e  ! h e ! k e d  for .
    !c
    !  x.....the point at which the b-splines are to be evaluated.
    !  left.....an integer chosen (usually) so that
    !                  t(left) .le. x .le. t(left+1)  .
    !c
    !c******  o u t p u t  ******
    !  biatx.....array of length  jout , with  biatx(i)  containing the val-
    !        ue at  x  of the polynomial of order  jout  which agrees with
    !        the b-spline  b(left-jout+i,jout,t)  on the interval (t(left),
    !        t(left+1)) .
    !c!
    !c******  m e t h o d  ******
    !  the recurrence relation
    !c
    !                       x - t(i)              t(i+j+1) - x
    !     b(i,j+1)(x)  =  -----------b(i,j)(x) + ---------------b(i+1,j)(x)
    !                     t(i+j)-t(i)            t(i+j+1)-t(i+1)
    !c
    !  is used (repeatedly) to generate the (j+1)-vector  b(left-j,j+1)(x),
    !  ...,b(left,j+1)(x)  from the j-vector  b(left-j+1,j)(x),...,
    !  b(left,j)(x), storing the new values in  biatx  over the old. the
    !  facts that
    !            b(i,1) = 1  if  t(i) .le. x .lt. t(i+1)
    !  and that
    !            b(i,j)(x) = 0  unless  t(i) .le. x .lt. t(i+j)
    !  are used. the particular organization of the calculations follows al-
    !  gorithm  (8)  in chapter x of the text.
    !c
            integer index,jhigh,left,   i,j,jmax,jp1
    parameter (jmax = 20)
    real biatx(jhigh),t(1),x,   deltal(jmax),deltar(jmax),saved,term
    !     real biatx(jhigh),t(1),x,   deltal(20),deltar(20),saved,term
    !     dimension biatx(jout), t(left+jout)
    !current fortran standard makes it impossible to specify the length of
    !  t  and of  biatx  precisely without the introduction of otherwise
    !  superfluous additional arguments.
    data j/1/
    save j,deltal,deltar
    !c
            go to (10,20), index
    10 j = 1
    biatx(1) = 1.
    if (j .ge. jhigh)                 go to 99
    !c
    20    jp1 = j + 1
    deltar(j) = t(left+j) - x
    deltal(j) = x - t(left+1-j)
    saved = 0.
    do 26 i=1,j
    term = biatx(i)/(deltar(i) + deltal(jp1-i))
    biatx(i) = saved + deltar(i)*term
    26       saved = deltal(jp1-i)*term
    biatx(jp1) = saved
    j = jp1
    if (j .lt. jhigh)              go to 20
    !c
    99                                   return
    end

 */