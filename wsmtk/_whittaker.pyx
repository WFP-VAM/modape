#cython: boundscheck=False
#cython: wraparound=False
#cython: cdivision=True
from cpython.array cimport array, clone

cimport numpy as np
import numpy as np
tFloat = np.float32
ctypedef np.float32_t dtype_t
from libc.math cimport log, pow, sqrt

cpdef ws2d(np.ndarray[dtype_t] y, float lmda, np.ndarray[dtype_t] w):
    cdef array flt_array_template = array('f', [])
    cdef int i, i1, i2, m, n
    cdef array z, d, c, e

    n = y.shape[0]
    m = n - 1

    z = clone(flt_array_template,n,zero=False)
    d = clone(flt_array_template,n,zero=False)
    c = clone(flt_array_template,n,zero=False)
    e = clone(flt_array_template,n,zero=False)

    d.data.as_floats[0] = w[0] + lmda
    c.data.as_floats[0] = (-2 * lmda) / d.data.as_floats[0]
    e.data.as_floats[0] = lmda /d.data.as_floats[0]
    z.data.as_floats[0] = w[0] * y[0]
    d.data.as_floats[1] = w[1] + 5 * lmda - d.data.as_floats[0] * (c.data.as_floats[0] * c.data.as_floats[0])
    c.data.as_floats[1] = (-4 * lmda - d.data.as_floats[0] * c.data.as_floats[0] * e.data.as_floats[0]) / d.data.as_floats[1]
    e.data.as_floats[1] =  lmda / d.data.as_floats[1]
    z.data.as_floats[1] = w[1] * y[1] - c.data.as_floats[0] * z.data.as_floats[0]
    for i in range(2,m-1):
        i1 = i - 1
        i2 = i - 2
        d.data.as_floats[i]= w[i] + 6 *  lmda - (c.data.as_floats[i1] * c.data.as_floats[i1]) * d.data.as_floats[i1] - (e.data.as_floats[i2] * e.data.as_floats[i2]) * d.data.as_floats[i2]
        c.data.as_floats[i] = (-4 *  lmda - d.data.as_floats[i1] * c.data.as_floats[i1] * e.data.as_floats[i1])/ d.data.as_floats[i]
        e.data.as_floats[i] =  lmda / d.data.as_floats[i]
        z.data.as_floats[i] = w[i] * y[i] - c.data.as_floats[i1] * z.data.as_floats[i1] - e.data.as_floats[i2] * z.data.as_floats[i2]
    i1 = m - 2
    i2 = m - 3
    d.data.as_floats[m - 1] = w[m - 1] + 5 *  lmda - (c.data.as_floats[i1] * c.data.as_floats[i1]) * d.data.as_floats[i1] - (e.data.as_floats[i2] * e.data.as_floats[i2]) * d.data.as_floats[i2]
    c.data.as_floats[m - 1] = (-2 *  lmda - d.data.as_floats[i1] * c.data.as_floats[i1] * e.data.as_floats[i1]) / d.data.as_floats[m - 1]
    z.data.as_floats[m - 1] = w[m - 1] * y[m - 1] - c.data.as_floats[i1] * z.data.as_floats[i1] - e.data.as_floats[i2] * z.data.as_floats[i2]
    i1 = m - 1
    i2 = m - 2
    d.data.as_floats[m] = w[m] +  lmda - (c.data.as_floats[i1] * c.data.as_floats[i1]) * d.data.as_floats[i1] - (e.data.as_floats[i2] * e.data.as_floats[i2]) * d.data.as_floats[i2]
    z.data.as_floats[m] = (w[m] * y[m] - c.data.as_floats[i1] * z.data.as_floats[i1] - e.data.as_floats[i2] * z.data.as_floats[i2]) / d.data.as_floats[m]
    z.data.as_floats[m - 1] = z.data.as_floats[m - 1] / d.data.as_floats[m - 1] - c.data.as_floats[m - 1] * z.data.as_floats[m]
    for i in range(m-2,-1,-1):
        z.data.as_floats[i] = z.data.as_floats[i] / d.data.as_floats[i] - c.data.as_floats[i] * z.data.as_floats[i + 1] - e.data.as_floats[i] * z.data.as_floats[i + 2]

    return z


cdef ws2d_internal(np.ndarray[dtype_t] y, float lmda, array[float] w):
    ## internal whittaker function for asymmetric smoothing
    cdef array flt_array_template = array('f', [])
    cdef int i, i1, i2, m, n
    cdef array z, d, c, e

    n = y.shape[0]
    m = n - 1

    z = clone(flt_array_template,n,zero=False)
    d = clone(flt_array_template,n,zero=False)
    c = clone(flt_array_template,n,zero=False)
    e = clone(flt_array_template,n,zero=False)

    d.data.as_floats[0] = w.data.as_floats[0] + lmda
    c.data.as_floats[0] = (-2 * lmda) / d.data.as_floats[0]
    e.data.as_floats[0] = lmda /d.data.as_floats[0]
    z.data.as_floats[0] = w.data.as_floats[0] * y[0]
    d.data.as_floats[1] = w.data.as_floats[1] + 5 * lmda - d.data.as_floats[0] * (c.data.as_floats[0] * c.data.as_floats[0])
    c.data.as_floats[1] = (-4 * lmda - d.data.as_floats[0] * c.data.as_floats[0] * e.data.as_floats[0]) / d.data.as_floats[1]
    e.data.as_floats[1] =  lmda / d.data.as_floats[1]
    z.data.as_floats[1] = w.data.as_floats[1] * y[1] - c.data.as_floats[0] * z.data.as_floats[0]
    for i in range(2,m-1):
        i1 = i - 1
        i2 = i - 2
        d.data.as_floats[i]= w.data.as_floats[i] + 6 *  lmda - (c.data.as_floats[i1] * c.data.as_floats[i1]) * d.data.as_floats[i1] - (e.data.as_floats[i2] * e.data.as_floats[i2]) * d.data.as_floats[i2]
        c.data.as_floats[i] = (-4 *  lmda - d.data.as_floats[i1] * c.data.as_floats[i1] * e.data.as_floats[i1])/ d.data.as_floats[i]
        e.data.as_floats[i] =  lmda / d.data.as_floats[i]
        z.data.as_floats[i] = w.data.as_floats[i] * y[i] - c.data.as_floats[i1] * z.data.as_floats[i1] - e.data.as_floats[i2] * z.data.as_floats[i2]
    i1 = m - 2
    i2 = m - 3
    d.data.as_floats[m - 1] = w.data.as_floats[m - 1] + 5 *  lmda - (c.data.as_floats[i1] * c.data.as_floats[i1]) * d.data.as_floats[i1] - (e.data.as_floats[i2] * e.data.as_floats[i2]) * d.data.as_floats[i2]
    c.data.as_floats[m - 1] = (-2 *  lmda - d.data.as_floats[i1] * c.data.as_floats[i1] * e.data.as_floats[i1]) / d.data.as_floats[m - 1]
    z.data.as_floats[m - 1] = w.data.as_floats[m - 1] * y[m - 1] - c.data.as_floats[i1] * z.data.as_floats[i1] - e.data.as_floats[i2] * z.data.as_floats[i2]
    i1 = m - 1
    i2 = m - 2
    d.data.as_floats[m] = w.data.as_floats[m] +  lmda - (c.data.as_floats[i1] * c.data.as_floats[i1]) * d.data.as_floats[i1] - (e.data.as_floats[i2] * e.data.as_floats[i2]) * d.data.as_floats[i2]
    z.data.as_floats[m] = (w.data.as_floats[m] * y[m] - c.data.as_floats[i1] * z.data.as_floats[i1] - e.data.as_floats[i2] * z.data.as_floats[i2]) / d.data.as_floats[m]
    z.data.as_floats[m - 1] = z.data.as_floats[m - 1] / d.data.as_floats[m - 1] - c.data.as_floats[m - 1] * z.data.as_floats[m]
    for i in range(m-2,-1,-1):
        z.data.as_floats[i] = z.data.as_floats[i] / d.data.as_floats[i] - c.data.as_floats[i] * z.data.as_floats[i + 1] - e.data.as_floats[i] * z.data.as_floats[i + 2]

    return z


cpdef ws2d_vc(np.ndarray[dtype_t] y, np.ndarray[dtype_t] w, array[float] llas):
    ## vcurve
    cdef array template = array('f', [])

    cdef array fits, pens, diff1, lamids, v, z
    cdef int m, m1, m2, nl, nl1, lix, i, k
    cdef float w_tmp, y_tmp, z_tmp, z2, llastep, f1, f2, p1, p2, l, l1, l2, vmin, lopt

    m = y.shape[0]
    m1 = m - 1
    m2 = m - 2
    nl = len(llas)
    nl1 = nl - 1
    i = 0
    k = 0

    template = array('f',[])

    fits = clone(template, nl, True)
    pens = clone(template,nl,True)
    z = clone(template,m,False)
    diff1 = clone(template,m1,True)
    lamids = clone(template,nl1,False)
    v = clone(template,nl1,False)

    # Compute v-curve

    for lix in range(nl):
        l = pow(10,llas.data.as_floats[lix])
        z[0:m] = ws2d(y,l,w)
        for i in range(m):
            w_tmp = w[i]
            y_tmp = y[i]
            z_tmp = z.data.as_floats[i]
            fits.data.as_floats[lix] += pow(w_tmp * (y_tmp - z_tmp),2)
        fits.data.as_floats[lix] = log(fits.data.as_floats[lix])

        for i in range(m1):
            z_tmp = z.data.as_floats[i]
            z2 = z.data.as_floats[i+1]
            diff1.data.as_floats[i] = z2 - z_tmp
        for i in range(m2):
            z_tmp = diff1.data.as_floats[i]
            z2 = diff1.data.as_floats[i+1]
            pens.data.as_floats[lix] += pow(z2 - z_tmp,2)
        pens.data.as_floats[lix] = log(pens.data.as_floats[lix])


    # Construct v-curve

    llastep = llas[1] - llas[0]

    for i in range(nl1):
        l1 = llas.data.as_floats[i]
        l2 = llas.data.as_floats[i+1]
        f1 = fits.data.as_floats[i]
        f2 = fits.data.as_floats[i+1]
        p1 = pens.data.as_floats[i]
        p2 = pens.data.as_floats[i+1]
        v.data.as_floats[i] = sqrt(pow(f2 - f1,2) + pow(p2 - p1,2)) / (log(10) * llastep)
        lamids.data.as_floats[i] = (l1+l2) / 2


    vmin = v.data.as_floats[k]
    for i in range(1,nl1):
        if v.data.as_floats[i] < vmin:
            vmin = v.data.as_floats[i]
            k = i

    lopt = pow(10,lamids.data.as_floats[k])

    z[0:m] = ws2d(y,lopt,w)

    return z, lopt


cpdef ws2d_vc_asy(np.ndarray[dtype_t] y, np.ndarray[dtype_t] w, array[float] llas, float p):
    ## vcurve with asymmetric smoothing

    cdef array template = array('f', [])

    cdef array fits, pens, diff1, lamids, v, z
    cdef int m, m1, m2, nl, nl1, lix, i, j, k
    cdef float w_tmp, y_tmp, z_tmp, z2, llastep, fit1, fit2, pen1, pen2, l, l1, l2, vmin, lopt, p1

    m = y.shape[0]
    m1 = m - 1
    m2 = m - 2
    nl = len(llas)
    nl1 = nl - 1
    i = 0
    k = 0
    j = 0
    p1 = 1-p

    template = array('f',[])

    fits = clone(template, nl, True)
    pens = clone(template,nl,True)
    z = clone(template,m,True)
    znew = clone(template,m,True)
    diff1 = clone(template,m1,True)
    lamids = clone(template,nl1,False)
    v = clone(template,nl1,False)
    wa = clone(template,m,False)
    ww = clone(template,m,False)

    # Compute v-curve

    for lix in range(nl):
        l = pow(10,llas.data.as_floats[lix])

        for i in range(10):
          for j in range(m):
            y_tmp = y[j]
            z_tmp = z.data.as_floats[j]
            if y_tmp > z_tmp:
              wa.data.as_floats[j] = p
            else:
              wa.data.as_floats[j] = p1
            ww.data.as_floats[j] = w[j] * wa.data.as_floats[j]

          znew[0:m] = ws2d_internal(y,l,ww)

          z_tmp = 0.0
          j = 0
          for j in range(m):
            z_tmp += abs(znew.data.as_floats[j] - z.data.as_floats[j])

          if z_tmp == 0.0:
            break

          z [0:m]= znew[0:m]

        for i in range(m):
            w_tmp = w[i]
            y_tmp = y[i]
            z_tmp = z.data.as_floats[i]
            fits.data.as_floats[lix] += pow(w_tmp * (y_tmp - z_tmp),2)
        fits.data.as_floats[lix] = log(fits.data.as_floats[lix])

        for i in range(m1):
            z_tmp = z.data.as_floats[i]
            z2 = z.data.as_floats[i+1]
            diff1.data.as_floats[i] = z2 - z_tmp
        for i in range(m2):
            z_tmp = diff1.data.as_floats[i]
            z2 = diff1.data.as_floats[i+1]
            pens.data.as_floats[lix] += pow(z2 - z_tmp,2)
        pens.data.as_floats[lix] = log(pens.data.as_floats[lix])


    # Construct v-curve

    llastep = llas[1] - llas[0]

    for i in range(nl1):
        l1 = llas.data.as_floats[i]
        l2 = llas.data.as_floats[i+1]
        fit1 = fits.data.as_floats[i]
        fit2 = fits.data.as_floats[i+1]
        pen1 = pens.data.as_floats[i]
        pen2 = pens.data.as_floats[i+1]
        v.data.as_floats[i] = sqrt(pow(fit2 - fit1,2) + pow(pen2 - pen1,2)) / (log(10) * llastep)
        lamids.data.as_floats[i] = (l1+l2) / 2


    vmin = v.data.as_floats[k]
    for i in range(1,nl1):
        if v.data.as_floats[i] < vmin:
            vmin = v.data.as_floats[i]
            k = i

    lopt = pow(10,lamids.data.as_floats[k])

    z[0:m] = ws2d_internal(y,lopt,ww)

    return z, lopt
