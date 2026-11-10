import numpy as np
from scipy.fft import fft, fftfreq
from scipy.signal import find_peaks

class TAD_X:
    """
    TAD-X
    Data-Driven FFT Scale Discovery + DP Tree + Domain Inference + Overlap Resolution
    """
    def __init__(self, resolution, prominence_factor=0.25, _min_tad_bp=120000, _max_tad_bp=2500000):
        self.resolution = resolution
        self.min_bin = max(3, _min_tad_bp // resolution)
        self.max_bin = _max_tad_bp // resolution
        self.prominence_factor = prominence_factor

    def _get_1d_boundaries(self, matrix, min_s):
        n = matrix.shape[0]
        base_signal = np.zeros(n)
        for i in range(min_s, n - min_s):
            base_signal[i] = np.mean(matrix[i-min_s:i, i:i+min_s])

        base_signal[:min_s] = base_signal[min_s]
        base_signal[n-min_s:] = base_signal[n-min_s-1]

        sig_centered = base_signal - np.mean(base_signal)
        yf = fft(sig_centered)
        xf = fftfreq(n, d=1.0)

        power = np.abs(yf)**2
        positive_freqs = xf[xf > 0]
        positive_power = power[xf > 0]

        freq_peaks, _ = find_peaks(positive_power, prominence=np.std(positive_power)*0.5)
        dominant_freqs = positive_freqs[freq_peaks]

        discovered_scales = []
        for f in dominant_freqs:
            scale = int(round(1.0 / f))
            if self.min_bin <= scale <= self.max_bin:
                discovered_scales.append(scale)

        if not discovered_scales:
            discovered_scales = [min_s, min_s + 2, min_s + 4]
        else:
            discovered_scales = sorted(list(set(discovered_scales)))[:4]

        combined_score = np.zeros(n)
        for w in discovered_scales:
            score = np.zeros(n)
            for i in range(w, n - w):
                left_TAD = np.mean(matrix[i-w:i, i-w:i])
                right_TAD = np.mean(matrix[i:i+w, i:i+w])
                cross_bg = np.mean(matrix[i-w:i, i:i+w])

                if cross_bg > 1e-10:
                    score[i] = (left_TAD + right_TAD) / (2 * cross_bg)

            if np.std(score) > 0:
                combined_score += (score - np.mean(score)) / np.std(score)

        prom = np.std(combined_score) * self.prominence_factor
        peaks, _ = find_peaks(combined_score, distance=min_s, prominence=prom)

        boundaries = sorted(list(peaks))
        if 0 not in boundaries: boundaries.insert(0, 0)
        if n not in boundaries: boundaries.append(n)

        return boundaries

    def detect(self, matrix):
        n = matrix.shape[0]

        # Stage 0: O/E matrix computation
        OE = np.zeros_like(matrix, dtype=float)
        for d in range(n):
            diag_idx = (np.arange(n-d), np.arange(d, n))
            diag_mean = np.mean(matrix[diag_idx])
            if diag_mean > 1e-10:
                OE[diag_idx] = matrix[diag_idx] / diag_mean
                if d > 0:
                    lower_idx = (np.arange(d, n), np.arange(n-d))
                    OE[lower_idx] = matrix[lower_idx] / diag_mean

        B_base = self._get_1d_boundaries(matrix, self.min_bin)
        N_b = len(B_base)
        is_tad = np.zeros((N_b, N_b), dtype=bool)

        for i in range(N_b-1):
            length = B_base[i+1] - B_base[i]
            if self.min_bin <= length <= self.max_bin:
                is_tad[i, i+1] = True

        def eval_corner(s1, e1, s2, e2):
            corner = OE[s1:e1+1, s2:e2+1]
            c_mean = np.mean(corner) if corner.size else 0.0
            L_left = e1 - s1 + 1
            L_right = e2 - s2 + 1
            flank_L = OE[max(0, s1-L_left):s1, s2:e2+1]
            flank_R = OE[s1:e1+1, e2+1:min(n, e2+1+L_right)]

            fL = np.mean(flank_L) if flank_L.size else 1.0
            fR = np.mean(flank_R) if flank_R.size else 1.0

            ext_contrast = c_mean / (max(fL, fR) + 1e-10)
            b1 = OE[s1:e1+1, s1:e1+1]
            b2 = OE[s2:e2+1, s2:e2+1]
            b1m = np.mean(b1) if b1.size else 1.0
            b2m = np.mean(b2) if b2.size else 1.0
            retention = c_mean / ((b1m+b2m)/2.0 + 1e-10)

            return ext_contrast, retention, c_mean

        # Stage 1: DP Tree Construction
        for span in range(2, N_b):
            for i in range(N_b-span):
                j = i + span
                length = B_base[j] - B_base[i]
                if length > self.max_bin: continue
                for k in range(i+1, j):
                    if not (is_tad[i, k] and is_tad[k, j]): continue
                    s1, e1 = B_base[i], B_base[k]-1
                    s2, e2 = B_base[k], B_base[j]-1
                    _, retention, cmean = eval_corner(s1, e1, s2, e2)
                    if retention >= 0.40 and cmean >= 1.05:
                        is_tad[i, j] = True
                        break

        # Stage 2: Domain Inference
        cand_tads = []
        for i in range(N_b-1):
            if is_tad[i, i+1]:
                cand_tads.append([B_base[i], B_base[i+1]-1])

        for span in range(2, N_b):
            for i in range(N_b-span):
                j = i+span
                if not is_tad[i, j]: continue
                valid = False
                for k in range(i+1, j):
                    if not (is_tad[i, k] and is_tad[k, j]): continue
                    s1, e1 = B_base[i], B_base[k]-1
                    s2, e2 = B_base[k], B_base[j]-1
                    ext_c, _, cmean = eval_corner(s1, e1, s2, e2)
                    tip_s = min(3, e1-s1+1)
                    tip_e = min(3, e2-s2+1)
                    tip = OE[s1:s1+tip_s, e2-tip_e+1:e2+1]
                    tip_mean = np.mean(tip) if tip.size else 0

                    if (ext_c >= 1.15 and cmean >= 1.10 and tip_mean >= 1.10):
                        valid = True
                        break
                if valid:
                    cand_tads.append([B_base[i], B_base[j]-1])

        # Stage 3: Overlap Resolution
        final_tads = list(cand_tads)
        overlap_tads = []
        for i in range(len(cand_tads)):
            for j in range(i+1, len(cand_tads)):
                s1, e1 = cand_tads[i]
                s2, e2 = cand_tads[j]
                if s1 < s2 <= e1 < e2:
                    ext_c, _, ctop = eval_corner(s1, s2-1, e1+1, e2)
                    if (ext_c >= 1.20 and ctop >= 1.15):
                        new_tad = [s1, e2]
                        if (new_tad not in final_tads and (e2-s1+1) <= self.max_bin):
                            final_tads.append(new_tad)
                    else:
                        if [s1, e1] not in overlap_tads: overlap_tads.append([s1, e1])
                        if [s2, e2] not in overlap_tads: overlap_tads.append([s2, e2])

        unique = []
        for t in final_tads:
            if t not in unique: unique.append(t)
        unique = sorted(unique, key=lambda x: (x[0], x[1]-x[0]))
        overlap_tads = sorted(overlap_tads, key=lambda x: (x[0], x[1]-x[0]))

        return {"all_tads": unique, "overlap_tads": overlap_tads, "base_boundaries": B_base}
