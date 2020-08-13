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

void Atom::Bsplines::setupKnotPoints(const std::vector<Complex> &gridPoints) {
    u32 numberOfGhostPointsInEachEnd = (m_order - 1);
    u32 numberOfPhysicalKnotPoints = m_numKnotPoints - 2 * numberOfGhostPointsInEachEnd;

    // For proper indexing
    u32 offset = numberOfGhostPointsInEachEnd;
    u32 stride = (u32) (gridPoints.size() / (numberOfPhysicalKnotPoints - 1));

    // Linear knot sequence
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

void Atom::Bsplines::LogKnotSequence() {
    u32 index = 0;
    for (auto point : m_knotPoints) {
        Logger::Log("Knotpoint %i: (%E, %E)", index, point.real(), point.imag());
        index += 1;
    }
}

u32 Atom::Bsplines::numberOfBsplines() {
    return m_numKnotPoints - m_order;
}

void Atom::Bsplines::bsplvb_Complex(Complex coordinate, u32 left) {
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

    u32 k = m_order;
    auto &t = m_knotPoints;
    auto &dR = m_bsplvb_dR;
    auto &dL = m_bsplvb_dL;

    for (u32 j = 0; j < k - 1; ++j) {
        dR[j] = t[left + j + 1] - x;
        dL[j] = x - t[left - j];

//        Logger::Trace("dR[%i]: %f, dL[%i]: %f", j, dR[j].Re, j, dL[j].Re);
        saved = Complex(0.0);

        for (int i = 0; i <= j; i++) {
            term = m_Sp[i] / (dR[i] + dL[j - i]);
//            Logger::Trace("term: %f", term.Re);
            m_Sp[i] = saved + dR[i] * term;
            saved = dL[j - i] * term;
        }

//        Logger::Trace("saved: %f", saved.Re);
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
    if(bsplineIndex == m_numBsplines-1) {
        f64 tolerance = 1e-6;
        if (std::abs(coordinate-m_knotPoints[m_knotPoints.size()-1]) < tolerance ) {
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
    u32 accessIndex = bsplineIndex - left_knotPoint_index + m_order;
    if(accessIndex < 0 || accessIndex > m_Sp.size()-1) {
        // We're not in a bspline that is nonzero on this coordinate. So return zero.
        return Complex(0.0);
    }

    // Calculate non-zero Bsplines on the coordinate.
    bsplvb_Complex(coordinate, left_knotPoint_index);

    // Get the correct Bspline from the k = bsplineOrder non-zero ones.
    Complex return_value = m_Sp[accessIndex];
    return return_value;
}


//
//void Bsplines::bsplvb_Complex(const Complex *knotPoints, const u32 bsplineOrder,
//                              const Complex x, const u32 left,
//                              Complex *outputArray) {
//    Complex saved;
//    Complex term;
//
//    // Assume index == 1, ie starting order is j = 0;
//    outputArray[0] = Complex(1.0, 0.0);
//
//    u32 k = bsplineOrder;
//    auto *t = knotPoints;
//    auto *dR = new Complex[k];
//    auto *dL = new Complex[k];
//
//    for (u32 j = 0; j < k - 1; ++j) {
//        dR[j] = t[left + j + 1] - x;
//        dL[j] = x - t[left - j];
//
////        Logger::Trace("dR[%i]: %f, dL[%i]: %f", j, dR[j].Re, j, dL[j].Re);
//        saved = Complex(0.0);
//
//        for (int i = 0; i <= j; i++) {
//            term = outputArray[i] / (dR[i] + dL[j - i]);
////            Logger::Trace("term: %f", term.Re);
//            outputArray[i] = saved + dR[i] * term;
//            saved = dL[j - i] * term;
//        }
//
////        Logger::Trace("saved: %f", saved.Re);
//        outputArray[j + 1] = saved;
//    }
//
//    delete[] dR;
//    delete[] dL;
//}
//

//Complex Bsplines::GetValueAtCoordinate(const Complex x_value, const Complex *knotPoints,
//                                       const u32 numKnotPoints, const u32 bsplineOrder,
//                                       const u32 bsplineIndex, const u32 numberOfBsplines) {
//    f64 difference_tolerance = 1e-8;
//
//    if (x_value.Re < knotPoints[0].Re) {
////        Logger::Trace("x_val is less than the first grid point. Return Complex(0.0, 0.0)");
//        return Complex(0.0, 0.0);
//    }
//
//    if (x_value.Re > knotPoints[numKnotPoints - 1].Re) {
////        Logger::Trace("x_val is greater than the last grid point. Return Complex(0.0, 0.0)");
//        return Complex(0.0, 0.0);
//    }
//
////    if (bsplineIndex == 0) {
////        // We're looking at the first bspline
////        if (abs(x_value.Re - knotPoints[0].Re) < difference_tolerance) {
//////            Logger::Trace("Difference too small between x_val and last knot point. Return Complex(1.0, 0.0)");
////            return Complex(1.0, 0.0);
////        }
////    }
//
//    if (bsplineIndex == numberOfBsplines - 1) {
//        // We're looking at the last bspline
//        if (abs(x_value.Re - knotPoints[numKnotPoints - 1].Re) < difference_tolerance) {
////            Logger::Trace("Difference too small between x_val and last knot point. Return Complex(1.0, 0.0)");
//            return Complex(1.0, 0.0);
//        }
//    }
//
//
//
//    Logger::Trace("left_knotPoint_index %i", left_knotPoint_index);
//    u32 accessIndex = bsplineIndex - left_knotPoint_index + bsplineOrder;
//    if (
//            accessIndex < 0 || accessIndex > bsplineOrder-1) {
////        Logger::Trace("Index condition not fullfilled for:");
////        Logger::Trace("indexCondition: %i, bsplineIndex: %i, left_knotPoint_index: %i, bsplineOrder: %i",
////                      indexCondition, bsplineIndex, left_knotPoint_index, bsplineOrder);
//        return Complex(0.0, 0.0);
//    }
//
////    Logger::Trace("Calculating Bspline value for");
////    Logger::Trace("bsplineIndex: %i, left_knotPoint_index: %i, bsplineOrder: %i",
////                  bsplineIndex, left_knotPoint_index, bsplineOrder);
//
//    auto *Sp = new Complex[bsplineOrder]; // TODO(anton): Explain this.
//    Bsplines::bsplvb_Complex(knotPoints, bsplineOrder, x_value, left_knotPoint_index, Sp);
//
//    Complex ret = Sp[accessIndex];
//    delete[] Sp;
//    return ret;
//}

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