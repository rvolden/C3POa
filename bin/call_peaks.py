#!/usr/bin/env python3
# Roger Volden

from savitzky_golay import savitzky_golay
from scipy.signal import find_peaks
import numpy as np

def call_peaks(scores, min_dist, iters, window, order):
    peaks = []
    for i in range(iters):
        scores = savitzky_golay(scores, window, order, deriv=0, rate=1)
    med_score = np.median(scores)
    if max(scores) < 6 * med_score:
        return peaks
    peaks, _ = find_peaks(scores, distance=min_dist, height=med_score * 3)
    return peaks
