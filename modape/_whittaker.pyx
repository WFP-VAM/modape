#cython: boundscheck=True
#cython: wraparound=False
#cython: cdivision=True
"""
Core whittaker functions

Author: Valentin Pesendorfer, April 2019

updated: September 2019
"""

from cpython.array cimport array, clone
from libc.math cimport log, pow, sqrt
cimport numpy as np
import numpy as np

tFloat = np.double
ctypedef np.double_t dtype_t

__all__ = ["lag1corr", "ws2d", "ws2dp", "ws2doptv", "ws2doptvp"]

cpdef lag1corr(np.ndarray[dtype_t] data1, np.ndarray[dtype_t] data2, double nd):
    """Calculates Lag-1 autocorrelation.

    Adapted from https://stackoverflow.com/a/29194624/5997555

    Args:
        data1: fist data series
        data2: second data series
        nd: no-data value (will be exluded from calulation)

    Returns:
        Lag-1 autocorrelation value
    """

    cdef int M, sub
    cdef double sum1, sum2, var_sum1, var_sum2, cross_sum, std1, std2, cross_mean

    M = data1.size

    sum1 = 0.
    sum2 = 0.
    sub = 0
    for i in range(M):
        if data1[i] != nd and data2[i] != nd:
            sum1 += data1[i]
            sum2 += data2[i]
        else:
            sub += 1
    mean1 = sum1 / (M-sub)
    mean2 = sum2 / (M-sub)

    var_sum1 = 0.
    var_sum2 = 0.
    cross_sum = 0.
    for i in range(M):
        if data1[i] != nd and data2[i] != nd:
            var_sum1 += (data1[i] - mean1) ** 2
            var_sum2 += (data2[i] - mean2) ** 2
            cross_sum += (data1[i] * data2[i])

    std1 = (var_sum1 / (M-sub)) ** .5
    std2 = (var_sum2 / (M-sub)) ** .5
    cross_mean = cross_sum / (M-sub)
    return (cross_mean - mean1 * mean2) / (std1 * std2)

cpdef ws2d(np.ndarray[dtype_t] y, double lmda, np.ndarray[dtype_t] w):
    cdef array dbl_array_template = array("d", [])
    cdef int i, i1, i2, m, n
    cdef array z, d, c, e

    n = y.shape[0]
    m = n - 1

    z = clone(dbl_array_template, n, zero=False)
    d = clone(dbl_array_template, n, zero=False)
    c = clone(dbl_array_template, n, zero=False)
    e = clone(dbl_array_template, n, zero=False)

    d.data.as_doubles[0] = w[0] + lmda
    c.data.as_doubles[0] = (-2 * lmda) / d.data.as_doubles[0]
    e.data.as_doubles[0] = lmda /d.data.as_doubles[0]
    z.data.as_doubles[0] = w[0] * y[0]
    d.data.as_doubles[1] = w[1] + 5 * lmda - d.data.as_doubles[0] * (c.data.as_doubles[0] * c.data.as_doubles[0])
    c.data.as_doubles[1] = (-4 * lmda - d.data.as_doubles[0] * c.data.as_doubles[0] * e.data.as_doubles[0]) / d.data.as_doubles[1]
    e.data.as_doubles[1] =  lmda / d.data.as_doubles[1]
    z.data.as_doubles[1] = w[1] * y[1] - c.data.as_doubles[0] * z.data.as_doubles[0]
    for i in range(2, m-1):
        i1 = i - 1
        i2 = i - 2
        d.data.as_doubles[i]= w[i] + 6 *  lmda - (c.data.as_doubles[i1] * c.data.as_doubles[i1]) * d.data.as_doubles[i1] - (e.data.as_doubles[i2] * e.data.as_doubles[i2]) * d.data.as_doubles[i2]
        c.data.as_doubles[i] = (-4 *  lmda - d.data.as_doubles[i1] * c.data.as_doubles[i1] * e.data.as_doubles[i1])/ d.data.as_doubles[i]
        e.data.as_doubles[i] =  lmda / d.data.as_doubles[i]
        z.data.as_doubles[i] = w[i] * y[i] - c.data.as_doubles[i1] * z.data.as_doubles[i1] - e.data.as_doubles[i2] * z.data.as_doubles[i2]
    i1 = m - 2
    i2 = m - 3
    d.data.as_doubles[m - 1] = w[m - 1] + 5 *  lmda - (c.data.as_doubles[i1] * c.data.as_doubles[i1]) * d.data.as_doubles[i1] - (e.data.as_doubles[i2] * e.data.as_doubles[i2]) * d.data.as_doubles[i2]
    c.data.as_doubles[m - 1] = (-2 *  lmda - d.data.as_doubles[i1] * c.data.as_doubles[i1] * e.data.as_doubles[i1]) / d.data.as_doubles[m - 1]
    z.data.as_doubles[m - 1] = w[m - 1] * y[m - 1] - c.data.as_doubles[i1] * z.data.as_doubles[i1] - e.data.as_doubles[i2] * z.data.as_doubles[i2]
    i1 = m - 1
    i2 = m - 2
    d.data.as_doubles[m] = w[m] +  lmda - (c.data.as_doubles[i1] * c.data.as_doubles[i1]) * d.data.as_doubles[i1] - (e.data.as_doubles[i2] * e.data.as_doubles[i2]) * d.data.as_doubles[i2]
    z.data.as_doubles[m] = (w[m] * y[m] - c.data.as_doubles[i1] * z.data.as_doubles[i1] - e.data.as_doubles[i2] * z.data.as_doubles[i2]) / d.data.as_doubles[m]
    z.data.as_doubles[m - 1] = z.data.as_doubles[m - 1] / d.data.as_doubles[m - 1] - c.data.as_doubles[m - 1] * z.data.as_doubles[m]
    for i in range(m-2, -1, -1):
        z.data.as_doubles[i] = z.data.as_doubles[i] / d.data.as_doubles[i] - c.data.as_doubles[i] * z.data.as_doubles[i + 1] - e.data.as_doubles[i] * z.data.as_doubles[i + 2]
    return z

cdef _ws2d(np.ndarray[dtype_t] y, double lmda, array[double] w):
    """Internal whittaker function for use in asymmetric smoothing.

    Args:
      y: time-series numpy array
      lmbda: lambda (s) value
      w: weights numpy array

    Returns:
        smoothed time-series array z
    """

    cdef array dbl_array_template = array("d", [])
    cdef int i, i1, i2, m, n
    cdef array z, d, c, e

    n = y.shape[0]
    m = n - 1

    z = clone(dbl_array_template, n, zero=False)
    d = clone(dbl_array_template, n, zero=False)
    c = clone(dbl_array_template, n, zero=False)
    e = clone(dbl_array_template, n, zero=False)

    d.data.as_doubles[0] = w.data.as_doubles[0] + lmda
    c.data.as_doubles[0] = (-2 * lmda) / d.data.as_doubles[0]
    e.data.as_doubles[0] = lmda /d.data.as_doubles[0]
    z.data.as_doubles[0] = w.data.as_doubles[0] * y[0]
    d.data.as_doubles[1] = w.data.as_doubles[1] + 5 * lmda - d.data.as_doubles[0] * (c.data.as_doubles[0] * c.data.as_doubles[0])
    c.data.as_doubles[1] = (-4 * lmda - d.data.as_doubles[0] * c.data.as_doubles[0] * e.data.as_doubles[0]) / d.data.as_doubles[1]
    e.data.as_doubles[1] =  lmda / d.data.as_doubles[1]
    z.data.as_doubles[1] = w.data.as_doubles[1] * y[1] - c.data.as_doubles[0] * z.data.as_doubles[0]
    for i in range(2, m-1):
        i1 = i - 1
        i2 = i - 2
        d.data.as_doubles[i]= w.data.as_doubles[i] + 6 *  lmda - (c.data.as_doubles[i1] * c.data.as_doubles[i1]) * d.data.as_doubles[i1] - (e.data.as_doubles[i2] * e.data.as_doubles[i2]) * d.data.as_doubles[i2]
        c.data.as_doubles[i] = (-4 *  lmda - d.data.as_doubles[i1] * c.data.as_doubles[i1] * e.data.as_doubles[i1])/ d.data.as_doubles[i]
        e.data.as_doubles[i] =  lmda / d.data.as_doubles[i]
        z.data.as_doubles[i] = w.data.as_doubles[i] * y[i] - c.data.as_doubles[i1] * z.data.as_doubles[i1] - e.data.as_doubles[i2] * z.data.as_doubles[i2]
    i1 = m - 2
    i2 = m - 3
    d.data.as_doubles[m - 1] = w.data.as_doubles[m - 1] + 5 *  lmda - (c.data.as_doubles[i1] * c.data.as_doubles[i1]) * d.data.as_doubles[i1] - (e.data.as_doubles[i2] * e.data.as_doubles[i2]) * d.data.as_doubles[i2]
    c.data.as_doubles[m - 1] = (-2 *  lmda - d.data.as_doubles[i1] * c.data.as_doubles[i1] * e.data.as_doubles[i1]) / d.data.as_doubles[m - 1]
    z.data.as_doubles[m - 1] = w.data.as_doubles[m - 1] * y[m - 1] - c.data.as_doubles[i1] * z.data.as_doubles[i1] - e.data.as_doubles[i2] * z.data.as_doubles[i2]
    i1 = m - 1
    i2 = m - 2
    d.data.as_doubles[m] = w.data.as_doubles[m] +  lmda - (c.data.as_doubles[i1] * c.data.as_doubles[i1]) * d.data.as_doubles[i1] - (e.data.as_doubles[i2] * e.data.as_doubles[i2]) * d.data.as_doubles[i2]
    z.data.as_doubles[m] = (w.data.as_doubles[m] * y[m] - c.data.as_doubles[i1] * z.data.as_doubles[i1] - e.data.as_doubles[i2] * z.data.as_doubles[i2]) / d.data.as_doubles[m]
    z.data.as_doubles[m - 1] = z.data.as_doubles[m - 1] / d.data.as_doubles[m - 1] - c.data.as_doubles[m - 1] * z.data.as_doubles[m]
    for i in range(m-2, -1, -1):
        z.data.as_doubles[i] = z.data.as_doubles[i] / d.data.as_doubles[i] - c.data.as_doubles[i] * z.data.as_doubles[i + 1] - e.data.as_doubles[i] * z.data.as_doubles[i + 2]
    return z

cpdef ws2dp(np.ndarray[dtype_t] y, double lmda, np.ndarray[dtype_t] w, double p):
  """Whittaker smoother with asymmetric smoothing and fixed lambda (S).

  Args:
      y: time-series numpy array
      l: smoothing parameter lambda (S)
      w: weights numpy array
      p: "Envelope" value

  Returns:
      Smoothed time-series array z
  """
  cdef array template = array("d", [])
  cdef int m, i, j
  cdef double y_tmp, z_tmp, p1

  m = y.shape[0]
  i = 0
  j = 0
  p1 = 1-p

  template = array("d", [])
  z = clone(template, m, True)
  znew = clone(template, m, True)
  wa = clone(template, m, False)
  ww = clone(template, m, False)

  # Calculate weights

  for i in range(10):
    for j in range(m):
      y_tmp = y[j]
      z_tmp = z.data.as_doubles[j]

      if y_tmp > z_tmp:
        wa.data.as_doubles[j] = p
      else:
        wa.data.as_doubles[j] = p1
      ww.data.as_doubles[j] = w[j] * wa.data.as_doubles[j]

    znew[0:m] = _ws2d(y, lmda, ww)
    z_tmp = 0.0
    j = 0
    for j in range(m):
      z_tmp += abs(znew.data.as_doubles[j] - z.data.as_doubles[j])

    if z_tmp == 0.0:
      break

    z[0:m]= znew[0:m]

  z[0:m] = _ws2d(y, lmda, ww)
  return z

cpdef ws2doptv(np.ndarray[dtype_t] y, np.ndarray[dtype_t] w, array[double] llas):
    """Whittaker smoother with normal V-curve optimization of lambda (S).

    Args:
        y: time-series numpy array
        w: weights numpy array
        llas: array with lambda values to iterate (S-range)

    Returns:
        Smoothed time-series array z and optimized lambda (S) value lopt
    """
    cdef array template = array("d", [])
    cdef array fits, pens, diff1, lamids, v, z
    cdef int m, m1, m2, nl, nl1, lix, i, k
    cdef double w_tmp, y_tmp, z_tmp, z2, llastep, f1, f2, p1, p2, l, l1, l2, vmin, lopt

    m = y.shape[0]
    m1 = m - 1
    m2 = m - 2
    nl = len(llas)
    nl1 = nl - 1
    i = 0
    k = 0

    template = array("d", [])

    fits = clone(template, nl, True)
    pens = clone(template, nl, True)
    z = clone(template, m, False)
    diff1 = clone(template, m1, True)
    lamids = clone(template, nl1, False)
    v = clone(template, nl1, False)

    # Compute v-curve
    for lix in range(nl):
        l = pow(10,llas.data.as_doubles[lix])
        z[0:m] = ws2d(y, l, w)
        for i in range(m):
            w_tmp = w[i]
            y_tmp = y[i]
            z_tmp = z.data.as_doubles[i]
            fits.data.as_doubles[lix] += pow(w_tmp * (y_tmp - z_tmp),2)
        fits.data.as_doubles[lix] = log(fits.data.as_doubles[lix])

        for i in range(m1):
            z_tmp = z.data.as_doubles[i]
            z2 = z.data.as_doubles[i+1]
            diff1.data.as_doubles[i] = z2 - z_tmp
        for i in range(m2):
            z_tmp = diff1.data.as_doubles[i]
            z2 = diff1.data.as_doubles[i+1]
            pens.data.as_doubles[lix] += pow(z2 - z_tmp,2)
        pens.data.as_doubles[lix] = log(pens.data.as_doubles[lix])

    # Construct v-curve
    llastep = llas[1] - llas[0]

    for i in range(nl1):
        l1 = llas.data.as_doubles[i]
        l2 = llas.data.as_doubles[i+1]
        f1 = fits.data.as_doubles[i]
        f2 = fits.data.as_doubles[i+1]
        p1 = pens.data.as_doubles[i]
        p2 = pens.data.as_doubles[i+1]
        v.data.as_doubles[i] = sqrt(pow(f2 - f1,2) + pow(p2 - p1,2)) / (log(10) * llastep)
        lamids.data.as_doubles[i] = (l1+l2) / 2

    vmin = v.data.as_doubles[k]
    for i in range(1, nl1):
        if v.data.as_doubles[i] < vmin:
            vmin = v.data.as_doubles[i]
            k = i

    lopt = pow(10, lamids.data.as_doubles[k])

    z[0:m] = ws2d(y, lopt, w)

    return z, lopt


cpdef ws2doptvp(np.ndarray[dtype_t] y, np.ndarray[dtype_t] w, array[double] llas, double p):
    """Whittaker smoother with asymmetric V-curve optimization of lambda (S).

    Args:
        y: time-series numpy array
        w: weights numpy array
        llas: array with lambda values to iterate (S-range)
        p: "Envelope" value

    Returns:
        Smoothed time-series array z and optimized lambda (S) value lopt
    """
    cdef array template = array("d", [])
    cdef array fits, pens, diff1, lamids, v, z
    cdef int m, m1, m2, nl, nl1, lix, i, j, k
    cdef double w_tmp, y_tmp, z_tmp, z2, llastep, fit1, fit2, pen1, pen2, l, l1, l2, vmin, lopt, p1

    m = y.shape[0]
    m1 = m - 1
    m2 = m - 2
    nl = len(llas)
    nl1 = nl - 1
    i = 0
    k = 0
    j = 0
    p1 = 1-p

    template = array("d", [])
    fits = clone(template, nl, True)
    pens = clone(template, nl, True)
    z = clone(template, m, True)
    znew = clone(template, m, True)
    diff1 = clone(template, m1, True)
    lamids = clone(template, nl1, False)
    v = clone(template, nl1, False)
    wa = clone(template, m, False)
    ww = clone(template, m, False)

    # Compute v-curve
    for lix in range(nl):
        l = pow(10,llas.data.as_doubles[lix])

        for i in range(10):
          for j in range(m):
            y_tmp = y[j]
            z_tmp = z.data.as_doubles[j]
            if y_tmp > z_tmp:
              wa.data.as_doubles[j] = p
            else:
              wa.data.as_doubles[j] = p1
            ww.data.as_doubles[j] = w[j] * wa.data.as_doubles[j]

          znew[0:m] = _ws2d(y, l, ww)
          z_tmp = 0.0
          j = 0
          for j in range(m):
            z_tmp += abs(znew.data.as_doubles[j] - z.data.as_doubles[j])

          if z_tmp == 0.0:
            break

          z[0:m]= znew[0:m]

        for i in range(m):
            w_tmp = w[i]
            y_tmp = y[i]
            z_tmp = z.data.as_doubles[i]
            fits.data.as_doubles[lix] += pow(w_tmp * (y_tmp - z_tmp),2)
        fits.data.as_doubles[lix] = log(fits.data.as_doubles[lix])

        for i in range(m1):
            z_tmp = z.data.as_doubles[i]
            z2 = z.data.as_doubles[i+1]
            diff1.data.as_doubles[i] = z2 - z_tmp
        for i in range(m2):
            z_tmp = diff1.data.as_doubles[i]
            z2 = diff1.data.as_doubles[i+1]
            pens.data.as_doubles[lix] += pow(z2 - z_tmp,2)
        pens.data.as_doubles[lix] = log(pens.data.as_doubles[lix])

    # Construct v-curve
    llastep = llas[1] - llas[0]

    for i in range(nl1):
        l1 = llas.data.as_doubles[i]
        l2 = llas.data.as_doubles[i+1]
        fit1 = fits.data.as_doubles[i]
        fit2 = fits.data.as_doubles[i+1]
        pen1 = pens.data.as_doubles[i]
        pen2 = pens.data.as_doubles[i+1]
        v.data.as_doubles[i] = sqrt(pow(fit2 - fit1,2) + pow(pen2 - pen1,2)) / (log(10) * llastep)
        lamids.data.as_doubles[i] = (l1+l2) / 2

    vmin = v.data.as_doubles[k]
    for i in range(1, nl1):
        if v.data.as_doubles[i] < vmin:
            vmin = v.data.as_doubles[i]
            k = i

    lopt = pow(10, lamids.data.as_doubles[k])

    del z
    z = clone(template, m, True)

    for i in range(10):
      for j in range(m):
        y_tmp = y[j]
        z_tmp = z.data.as_doubles[j]

        if y_tmp > z_tmp:
          wa.data.as_doubles[j] = p
        else:
          wa.data.as_doubles[j] = p1
        ww.data.as_doubles[j] = w[j] * wa.data.as_doubles[j]

      znew[0:m] = _ws2d(y, lopt, ww)
      z_tmp = 0.0
      j = 0
      for j in range(m):
        z_tmp += abs(znew.data.as_doubles[j] - z.data.as_doubles[j])

      if z_tmp == 0.0:
        break

      z[0:m]= znew[0:m]

    z[0:m] = _ws2d(y, lopt, ww)
    return z, lopt
