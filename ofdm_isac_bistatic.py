#!/usr/bin/env python3
"""
ofdm_isac_bistatic.py — Phase 1 ISAC bistatic extension dari ofdm_robust_multi.py

Status:
  Phase 1A — IMPLEMENTED  : FS sweep mode, tentukan max FS hardware-feasible
  Phase 1B — IMPLEMENTED  : V2I "STEI" comm payload (QPSK center 8 SC)
  Phase 1C — IMPLEMENTED  : Range estimation (CIR + delta-range, in-band sync)
  Phase 1D — IMPLEMENTED  : Real-time matplotlib live plot (--plot flag)
  Phase 1E — PENDING      : Lab validation tests (moving target scenarios)

Frame structure (Nfft=64 fixed, FS-scalable):
  [STF | LTF1 | LTF2 | DATA × 26]
  Total samples = 240 + 26*80 = 2320

Phase 0/1A frame layout (default):
  All DATA_REL (46 SC) → BPSK known sequence (deterministic seed=42)
  Pilot SC (4)         → 1+0j (CPE reference)

Phase 1B frame layout (--phase1b):
  COMM_SC  (8 center)  → QPSK packet "STEI"+ctr+CRC16 + random pad (416 bit cap)
  SENSE_SC (38)        → BPSK known (sense reference, seed=123)
  Pilot SC (4)         → 1+0j
  TX cycles 256 frames (counter 0..255), each ~50 μs @ FS=40 MHz

V2I packet (56 bit):
  [ASCII text 32b | counter 8b | CRC16-CCITT 16b]

Phase 1C range estimation:
  CIR = IFFT(Hanning_window(H_est), 4× zero-pad) → 256-bin
  Range bin = c / (FS × oversample_factor)
  @ FS=40 MHz, 4× pad → 1.875 m/bin (interpolated; theoretical δR = 4.8 m)
  Detect: direct path (max) + echoes (10 dB above noise floor, sidelobe skip)
  Output: list of (delta_range_m, peak_db) per frame

CFO unambiguous range = fs/N = fs/64 (Schmidl-Cox dengan delay=32):
  @  2 MHz : ±31.25 kHz
  @ 40 MHz : ±625 kHz
B210 ±2 ppm @ 5.9 GHz → CFO actual ±11.8 kHz. Aman semua FS.

Usage:
  # Phase 1A — FS sweep (sudah validated, fs_winner = 40 MHz untuk hardware ini)
  python3 ofdm_isac_bistatic.py --fs-sweep \\
      --fs-candidates 5e6,10e6,20e6,30e6,40e6 \\
      --frames-per-fs 200 --tx-gain 80 --rx-gain 70

  # Phase 1B+1C — V2I "STEI" comm + range estimation @ FS_winner
  python3 ofdm_isac_bistatic.py --phase1b --fs 40e6 --frames 500 \\
      --tx-gain 80 --rx-gain 70 --text "STEI"

  # Phase 1B+1C+1D — Tambah live plot untuk demo BRIN
  python3 ofdm_isac_bistatic.py --phase1b --plot --fs 40e6 --frames 500 \\
      --tx-gain 80 --rx-gain 70

  # Plot dengan echo threshold lebih sensitif (default 6 dB)
  python3 ofdm_isac_bistatic.py --phase1b --plot --fs 40e6 --frames 500 \\
      --tx-gain 80 --rx-gain 70 --echo-threshold-db 4

  # Phase 1B dengan custom CSV log
  python3 ofdm_isac_bistatic.py --phase1b --fs 40e6 --frames 1000 \\
      --tx-gain 80 --rx-gain 70 --log-csv lab_test_1.csv

  # Single FS run Phase 0-style (no comm split)
  python3 ofdm_isac_bistatic.py --fs 40e6 --frames 100 \\
      --tx-gain 80 --rx-gain 70

  # AWGN self-test (no hardware)
  python3 ofdm_isac_bistatic.py --simulate --fs 20e6
  python3 ofdm_isac_bistatic.py --simulate --fs 40e6 --phase1b
"""
import multiprocessing as mp
import threading
import time
import os
import sys
import argparse
import csv
import json
import numpy as np
from collections import deque

# ═══════════════════════════════════════════════════════════════════
# HARDWARE CONFIG (sama dengan Phase 0)
# ═══════════════════════════════════════════════════════════════════
TX_SERIAL    = "000000037"
RX_SERIAL    = "HQHGTFH"
TX_IMAGE_DIR = "/home/telmat/uhd_images/asli"
RX_IMAGE_DIR = "/home/telmat/uhd_images/libre"
TX_ANT       = "TX/RX"
RX_ANT       = "TX/RX"
FC           = 5.9e9


def _find_fpga(directory):
    import glob
    for name in ("usrp_b210_fpga.bin", "usrp_b210_fpga.bit",
                 "usrp_b200_fpga.bin", "usrp_b200_fpga.bit"):
        full = os.path.join(directory, name)
        if os.path.isfile(full):
            return os.path.abspath(full)
    for pat in ("*b210*.bin", "*b210*.bit", "*b200*.bin", "*b200*.bit"):
        hits = sorted(glob.glob(os.path.join(directory, pat)))
        if hits:
            return os.path.abspath(hits[0])
    return None


TX_FPGA = _find_fpga(TX_IMAGE_DIR) if os.path.isdir(TX_IMAGE_DIR) else None
RX_FPGA = _find_fpga(RX_IMAGE_DIR) if os.path.isdir(RX_IMAGE_DIR) else None


# ═══════════════════════════════════════════════════════════════════
# OFDM PARAMETERS — Nfft fixed, FS variable (set via init_params)
# ═══════════════════════════════════════════════════════════════════
NSC     = 64
NCP     = 16
LSYM    = NSC + NCP            # 80
NSYM    = 26
N_HALF  = NSC // 2             # 32 = S&C delay

DATA_REL  = [i for i in range(-NSC // 2, NSC // 2)
             if 1 <= abs(i) <= 25 and i not in (-21, -7, 7, 21)]
PILOT_REL = [-21, -7, 7, 21]
N_DATA    = len(DATA_REL)
ACTIVE_REL = sorted(set(DATA_REL + PILOT_REL))
N_ACTIVE   = len(ACTIVE_REL)

STF_OFF   = 0
LTF1_OFF  = LSYM
LTF2_OFF  = 2 * LSYM
DATA_OFF  = 3 * LSYM
FRAME_LEN = DATA_OFF + NSYM * LSYM   # 2320

# ── Runtime-variable globals (set by init_params) ─────────────────
FS           = None
BW           = None
MAX_CFO      = None
KNOWN_BITS   = None
STF          = None
X_STF        = None
LTF_SYM      = None
X_LTF        = None
LTF_TD       = None
LTF_TD_NO_CP = None
TX_FRAME     = None

# ═══════════════════════════════════════════════════════════════════
# PHASE 1B/1C — V2I COMM + SENSING
# ═══════════════════════════════════════════════════════════════════
# Subcarrier split:
#   COMM_SC  = 8 center SC (avoid DC, |k|≤4)  → QPSK packet payload
#   SENSE_SC = sisanya dari DATA_REL          → BPSK known (sense reference)
#   PILOT    = [-21, -7, 7, 21]               → 1+0j
COMM_SC  = [-4, -3, -2, -1, 1, 2, 3, 4]
SENSE_SC = sorted([sc for sc in DATA_REL if sc not in COMM_SC])
N_COMM   = len(COMM_SC)        # 8
N_SENSE  = len(SENSE_SC)       # 38
COMM_BITS_PER_FRAME  = N_COMM * NSYM * 2   # 8 × 26 × 2 = 416 bit/frame (QPSK)
SENSE_BITS_PER_FRAME = N_SENSE * NSYM      # 38 × 26 = 988 bit/frame (BPSK)

# V2I packet format: [text 32b | counter 8b | CRC16 16b] = 56 bit
PKT_TEXT_BITS  = 32
PKT_CTR_BITS   = 8
PKT_CRC_BITS   = 16
PKT_TOTAL_BITS = PKT_TEXT_BITS + PKT_CTR_BITS + PKT_CRC_BITS  # 56

N_FRAMES_CYCLE = 256   # TX cycles through 256 pre-built frames (counter 0..255)
CIR_OVERSAMPLE = 4     # IFFT zero-pad factor untuk smooth CIR (Phase 1C)
CIR_NFFT       = NSC * CIR_OVERSAMPLE  # 256

# Phase 1B/1C runtime globals
PHASE1B_MODE       = False
TX_TEXT            = "STEI"
SENSE_KNOWN_BITS   = None
TX_FRAMES_P1B      = None         # list of N_FRAMES_CYCLE pre-built frames
LTF_FREQ_WINDOW    = None         # Hanning window for CIR sidelobe suppression


# ─────────────────────────────────────────────────────────────────
# CRC16-CCITT-FALSE (poly 0x1021, init 0xFFFF, no reflect, no xorout)
# ─────────────────────────────────────────────────────────────────
def crc16_ccitt(data_bytes):
    crc = 0xFFFF
    for b in data_bytes:
        crc ^= (b << 8)
        for _ in range(8):
            if crc & 0x8000:
                crc = ((crc << 1) ^ 0x1021) & 0xFFFF
            else:
                crc = (crc << 1) & 0xFFFF
    return crc


def encode_packet(text, counter):
    """Build 56-bit packet: ASCII(4) + counter(1) + CRC16(2)."""
    text_bytes = text.encode('ascii', errors='replace')[:4].ljust(4, b'\x00')
    ctr_byte = bytes([counter & 0xFF])
    payload = text_bytes + ctr_byte           # 5 bytes = 40 bit
    crc = crc16_ccitt(payload)
    crc_bytes = bytes([(crc >> 8) & 0xFF, crc & 0xFF])
    full_bytes = payload + crc_bytes          # 7 bytes = 56 bit
    bits = np.zeros(56, dtype=np.uint8)
    for i, byte in enumerate(full_bytes):
        for j in range(8):
            bits[i*8 + j] = (byte >> (7 - j)) & 1
    return bits


def decode_packet(bits56):
    """Decode 56-bit array → (text, counter, crc_ok)."""
    if len(bits56) < 56:
        return ("????", 0, False)
    bytes_arr = bytearray(7)
    for i in range(7):
        b = 0
        for j in range(8):
            b = (b << 1) | int(bits56[i*8 + j])
        bytes_arr[i] = b
    text = bytes_arr[:4].decode('ascii', errors='replace')
    counter = bytes_arr[4]
    rx_crc = (bytes_arr[5] << 8) | bytes_arr[6]
    expected = crc16_ccitt(bytes(bytes_arr[:5]))
    return (text, counter, rx_crc == expected)


# ─────────────────────────────────────────────────────────────────
# QPSK map/demap (Gray-coded, unit power)
# ─────────────────────────────────────────────────────────────────
# bit pair (b0,b1):  (0,0)->+1+1j, (0,1)->+1-1j, (1,0)->-1+1j, (1,1)->-1-1j  / sqrt(2)
QPSK_TABLE = np.array([1+1j, 1-1j, -1+1j, -1-1j], dtype=complex) / np.sqrt(2)


def qpsk_map(bits_pairs):
    """bits_pairs: shape (N, 2) → N complex symbols."""
    idx = bits_pairs[:, 0] * 2 + bits_pairs[:, 1]
    return QPSK_TABLE[idx]


def qpsk_demap(symbols):
    """N complex → shape (N, 2) bits via hard sign decision."""
    bits = np.zeros((len(symbols), 2), dtype=np.uint8)
    bits[:, 0] = (np.real(symbols) < 0).astype(np.uint8)
    bits[:, 1] = (np.imag(symbols) < 0).astype(np.uint8)
    return bits


def _build_phase1b_frame(comm_payload_bits_416, sense_known_bits_988):
    """Construct Phase 1B time-domain frame.

    comm_payload_bits_416: 416 bit (QPSK on 8 center SC × 26 sym × 2 bit)
    sense_known_bits_988:  988 bit (BPSK on 38 sense SC × 26 sym × 1 bit)
    """
    # QPSK comm symbols
    comm_pairs = comm_payload_bits_416.reshape(-1, 2)         # (208, 2)
    comm_syms = qpsk_map(comm_pairs).reshape(NSYM, N_COMM)    # (26, 8)

    # BPSK sense symbols
    sense_syms = (1.0 - 2.0 * sense_known_bits_988.astype(float)).astype(complex)
    sense_syms = sense_syms.reshape(NSYM, N_SENSE)            # (26, 38)

    parts = [STF, LTF_SYM, LTF_SYM]
    for m in range(NSYM):
        f = np.zeros(NSC, dtype=complex)
        for i, sc in enumerate(COMM_SC):
            f[sc % NSC] = comm_syms[m, i]
        for i, sc in enumerate(SENSE_SC):
            f[sc % NSC] = sense_syms[m, i]
        for p in PILOT_REL:
            f[p % NSC] = 1.0 + 0j
        td = np.fft.ifft(f) * np.sqrt(NSC)
        parts.append(np.concatenate([td[-NCP:], td]).astype(np.complex64))
    return np.concatenate(parts).astype(np.complex64)


def _build_all_phase1b_frames(text):
    """Pre-build N_FRAMES_CYCLE frames (counter 0..255). Returns list."""
    np.random.seed(123)  # different seed dari Phase 0 KNOWN_BITS
    sense_known = np.random.randint(0, 2, SENSE_BITS_PER_FRAME).astype(np.uint8)
    pad_bits = np.random.randint(0, 2, COMM_BITS_PER_FRAME - PKT_TOTAL_BITS).astype(np.uint8)

    frames = []
    for ctr in range(N_FRAMES_CYCLE):
        pkt_bits = encode_packet(text, ctr)
        # Layout: [packet 56b | random pad 360b] = 416 bit
        comm_bits = np.concatenate([pkt_bits, pad_bits])
        raw = _build_phase1b_frame(comm_bits, sense_known)
        max_amp = float(np.max(np.abs(raw)))
        frame = (raw / max_amp * 0.95).astype(np.complex64)
        frames.append(frame)
    return frames, sense_known


# ─────────────────────────────────────────────────────────────────
# Phase 1C: CIR + range estimation
# ─────────────────────────────────────────────────────────────────
def cir_from_h_est(H_est, fs, n_pad=CIR_NFFT, window=True):
    """Return oversampled CIR magnitude (linear) + bin-to-meter scale.

    Window in freq domain (Hanning over active SC) suppress sidelobe.
    Zero-pad to n_pad untuk smooth peaks (interpolation, NOT extra resolution).
    """
    H = H_est.copy()
    if window:
        # Apply Hanning over active SC only (preserve spectral mask)
        win = np.zeros(NSC, dtype=float)
        win_active = np.hanning(N_ACTIVE)
        for i, k in enumerate(ACTIVE_REL):
            win[k % NSC] = win_active[i]
        H = H * win

    # Zero-pad in freq domain (FFT-shift convention)
    H_shifted = np.fft.fftshift(H)            # DC at center NSC/2
    pad_lo = (n_pad - NSC) // 2
    pad_hi = n_pad - NSC - pad_lo
    H_padded = np.concatenate([
        np.zeros(pad_lo, dtype=complex),
        H_shifted,
        np.zeros(pad_hi, dtype=complex),
    ])
    H_back = np.fft.ifftshift(H_padded)
    cir = np.fft.ifft(H_back) * (n_pad / NSC)  # preserve amplitude

    # Each CIR bin = 1/(fs × oversample) sec → range bin = c/(fs × oversample)
    bin_to_meter = 3e8 / (fs * (n_pad / NSC))
    return np.abs(cir), bin_to_meter


def estimate_ranges(cir_mag, bin_to_meter, snr_threshold_db=6,
                    direct_skip_bins=2, max_echoes=8): # turun dari 4 ke 1, karena ruangan kecil (1)-> blind zone ~1.9 m
    """Detect direct path + echo peaks, return delta-range list.

    Default threshold 6 dB (lebih sensitif untuk indoor lab realistic).
    Noise floor pakai 25th percentile (robust thd clutter peaks).
    direct_skip_bins=4 @ 4× oversample = ~7.5 m blind zone (FS=40 MHz).

    Returns dict:
      direct_bin, direct_db, noise_floor_db, echoes=[(delta_range_m, peak_db), ...]
    """
    cir_db = 20 * np.log10(cir_mag + 1e-12)
    n = len(cir_mag)

    # Direct path = global max
    direct_bin = int(np.argmax(cir_mag))
    direct_db = float(cir_db[direct_bin])

    # Noise floor: 25th percentile of bins far from direct (more robust than median)
    far_lo = max(0, direct_bin - n // 3)
    far_hi = min(n, direct_bin + n // 3)
    far_region = np.concatenate([cir_db[:far_lo], cir_db[far_hi:]])
    if len(far_region) >= 8:
        noise_floor_db = float(np.percentile(far_region, 25))
    else:
        noise_floor_db = direct_db - 30
    threshold_db = noise_floor_db + snr_threshold_db

    # Echo search: forward bins after direct (skip sidelobe region)
    echoes = []
    search_start = direct_bin + direct_skip_bins
    search_end = min(n, direct_bin + n // 2)  # forward only, half range
    for b in range(search_start + 1, search_end - 1):
        if cir_db[b] > threshold_db:
            # Local maximum check
            if cir_mag[b] > cir_mag[b-1] and cir_mag[b] > cir_mag[b+1]:
                delta_bin = b - direct_bin
                delta_range = delta_bin * bin_to_meter
                echoes.append((float(delta_range), float(cir_db[b])))

    # Keep top-K strongest echoes
    echoes.sort(key=lambda x: x[1], reverse=True)
    echoes = echoes[:max_echoes]
    echoes.sort(key=lambda x: x[0])  # sort by range for display

    return {
        "direct_bin": direct_bin,
        "direct_db": direct_db,
        "noise_floor_db": noise_floor_db,
        "threshold_db": threshold_db,
        "echoes": echoes,
        "n_echoes": len(echoes),
    }



# ─────────────────────────────────────────────────────────────────
# Phase 1F: Forward-Scattering / LoS disruption detector
# ─────────────────────────────────────────────────────────────────
class ForwardScatterDetector:
    """One-look object detector for 2-USRP forward-scattering JCAS.

    Prinsip:
      1) Ambil fitur LoS dari amplitude dan phase channel estimate.
      2) Hilangkan komponen statis dengan moving median.
      3) Hitung score gangguan LoS.
      4) Pakai adaptive CFAR + debounce agar status tidak flicker.

    Output utama untuk demo:
      - KALIBRASI : baseline masih dikumpulkan, jangan lewatkan objek dulu.
      - CLEAR     : tidak ada benda di antara TX dan RX.
      - ADA_OBJEK : ada gangguan/movement pada jalur LoS TX-RX.
    """

    def __init__(self, fs_frame, ma_len=30, cfar_len=160, threshold_k=5.0,
                 min_score=2.2, amp_weight=1.0, doppler_weight=0.55,
                 baseline_frames=80, consec_on=3, consec_off=12,
                 hold_frames=18):
        self.fs_frame = float(fs_frame)
        self.ma_len = int(max(3, ma_len))
        self.cfar_len = int(max(self.ma_len + 10, cfar_len))
        self.threshold_k = float(threshold_k)
        self.min_score = float(min_score)
        self.amp_weight = float(amp_weight)
        self.doppler_weight = float(doppler_weight)
        self.baseline_frames = int(max(10, baseline_frames))
        self.consec_on = int(max(1, consec_on))
        self.consec_off = int(max(1, consec_off))
        self.hold_frames = int(max(0, hold_frames))

        self.amp_db_hist = deque(maxlen=self.cfar_len)
        self.phase_unwrapped_hist = deque(maxlen=self.cfar_len)
        self.doppler_hist = deque(maxlen=self.cfar_len)
        self.score_hist = deque(maxlen=self.cfar_len)
        self.prev_phase = None
        self.prev_unwrapped_phase = None
        self.frame_count = 0
        self.on_count = 0
        self.off_count = 0
        self.hold_count = 0
        self.object_state = False

    @staticmethod
    def _robust_median(x, default=0.0):
        arr = np.asarray(list(x), dtype=float)
        arr = arr[np.isfinite(arr)]
        if arr.size == 0:
            return float(default)
        return float(np.median(arr))

    @staticmethod
    def _mad_sigma(x):
        arr = np.asarray(list(x), dtype=float)
        arr = arr[np.isfinite(arr)]
        if arr.size < 5:
            return 0.0
        med = np.median(arr)
        return float(1.4826 * np.median(np.abs(arr - med)) + 1e-12)

    @staticmethod
    def _wrap_phase_delta(phi_now, phi_prev):
        return float(np.angle(np.exp(1j * (phi_now - phi_prev))))

    def update(self, avg_amp=None, H_est=None, cir_mag=None):
        self.frame_count += 1

        # 1) LoS amplitude metric. Phase-0 uses avg_amp; phase1b may use direct CIR peak.
        if cir_mag is not None and len(cir_mag) > 0:
            amp_lin = float(np.max(np.asarray(cir_mag)))
        elif avg_amp is not None:
            amp_lin = float(avg_amp)
        else:
            amp_lin = 0.0
        amp_db = float(20 * np.log10(max(amp_lin, 1e-12)))

        # 2) Phase metric from weighted mean channel estimate.
        if H_est is not None:
            active_h = np.asarray([H_est[k % NSC] for k in ACTIVE_REL], dtype=complex)
            mag = np.abs(active_h)
            good = np.isfinite(mag) & (mag > np.percentile(mag, 35) if mag.size else False)
            if np.any(good) and np.sum(mag[good]) > 1e-12:
                los_complex = np.sum(active_h[good] * mag[good]) / np.sum(mag[good])
            elif np.sum(mag) > 1e-12:
                los_complex = np.sum(active_h * mag) / np.sum(mag)
            else:
                los_complex = np.mean(active_h)
            phase_raw = float(np.angle(los_complex))
        else:
            phase_raw = 0.0

        if self.prev_phase is None:
            phase_unwrapped = phase_raw
            phase_delta = 0.0
        else:
            phase_delta = self._wrap_phase_delta(phase_raw, self.prev_phase)
            phase_unwrapped = self.prev_unwrapped_phase + phase_delta
        self.prev_phase = phase_raw
        self.prev_unwrapped_phase = phase_unwrapped

        # 3) Doppler proxy: f_D = Δphase/(2π) × frame_rate.
        doppler_hz = float((phase_delta / (2 * np.pi)) * self.fs_frame)

        # 4) High-pass by moving median background.
        amp_bg = self._robust_median(list(self.amp_db_hist)[-self.ma_len:], amp_db)
        phase_bg = self._robust_median(list(self.phase_unwrapped_hist)[-self.ma_len:], phase_unwrapped)
        dopp_bg = self._robust_median(list(self.doppler_hist)[-self.ma_len:], doppler_hz)

        amp_hp_db = float(amp_db - amp_bg)
        phase_hp_rad = float(phase_unwrapped - phase_bg)
        doppler_hp_hz = float(doppler_hz - dopp_bg)

        # 5) Robust score. Doppler is down-weighted because phase can be jumpy indoors.
        amp_sigma = max(self._mad_sigma(self.amp_db_hist), 0.08)       # dB
        dopp_sigma = max(self._mad_sigma(self.doppler_hist), 1.5)      # Hz
        score = float(np.sqrt(
            (self.amp_weight * abs(amp_hp_db) / amp_sigma) ** 2 +
            (self.doppler_weight * abs(doppler_hp_hz) / dopp_sigma) ** 2
        ))

        # 6) Adaptive CFAR threshold from stable baseline/history.
        baseline_ready = len(self.score_hist) >= self.baseline_frames
        if baseline_ready:
            score_floor = self._robust_median(self.score_hist, 0.0)
            score_sigma = self._mad_sigma(self.score_hist)
            threshold = max(self.min_score, float(score_floor + self.threshold_k * score_sigma))
        else:
            threshold = self.min_score

        raw_hit = bool(baseline_ready and score > threshold)

        # 7) Debounce + hold. This makes the display clear and prevents flicker.
        if raw_hit:
            self.on_count += 1
            self.off_count = 0
        else:
            self.off_count += 1
            self.on_count = 0

        if self.on_count >= self.consec_on:
            self.object_state = True
            self.hold_count = self.hold_frames
        elif self.object_state:
            if self.hold_count > 0:
                self.hold_count -= 1
            elif self.off_count >= self.consec_off:
                self.object_state = False

        if not baseline_ready:
            status = "KALIBRASI"
            object_detected = False
        else:
            object_detected = bool(self.object_state)
            status = "ADA_OBJEK" if object_detected else "CLEAR"

        # Update histories after decision so current event does not train threshold first.
        self.amp_db_hist.append(amp_db)
        self.phase_unwrapped_hist.append(phase_unwrapped)
        self.doppler_hist.append(doppler_hz)
        self.score_hist.append(score)

        ratio = float(score / max(threshold, 1e-9))
        progress = float(min(1.0, len(self.score_hist) / max(1, self.baseline_frames)))

        return {
            "amp_db": amp_db,
            "amp_hp_db": amp_hp_db,
            "phase_rad": phase_raw,
            "phase_unwrapped_rad": float(phase_unwrapped),
            "phase_hp_rad": phase_hp_rad,
            "doppler_hz": doppler_hz,
            "doppler_hp_hz": doppler_hp_hz,
            "score": score,
            "threshold": float(threshold),
            "score_ratio": ratio,
            "baseline_ready": baseline_ready,
            "baseline_progress": progress,
            "raw_hit": raw_hit,
            "object_detected": object_detected,
            "status": status,
        }

def _build_stf():
    """STF: only even subcarriers populated → time-domain x[n] = x[n+N/2]."""
    X = np.zeros(NSC, dtype=complex)
    even_active = [k for k in ACTIVE_REL if k % 2 == 0]
    L = len(even_active)
    n_zc = np.arange(L)
    zc = np.exp(-1j * np.pi * 5 * n_zc * (n_zc + 1) / L)
    for i, k in enumerate(even_active):
        X[k % NSC] = zc[i] * np.sqrt(2.0)
    td = np.fft.ifft(X) * np.sqrt(NSC)
    cp = td[-NCP:]
    return np.concatenate([cp, td]).astype(np.complex64), X


def _build_ltf():
    """LTF: full-band ZC, 2 identical symbols for averaging."""
    X = np.zeros(NSC, dtype=complex)
    L = N_ACTIVE
    n_zc = np.arange(L)
    zc = np.exp(-1j * np.pi * 25 * n_zc * (n_zc + 1) / L)
    for i, k in enumerate(ACTIVE_REL):
        X[k % NSC] = zc[i]
    td = np.fft.ifft(X) * np.sqrt(NSC)
    cp = td[-NCP:]
    return np.concatenate([cp, td]).astype(np.complex64), X, td


def _build_frame_signal(bits):
    """Construct time-domain frame from bits (BPSK on data SC, 1+0j on pilots)."""
    syms = (1.0 - 2.0 * bits.astype(float)).astype(complex).reshape(NSYM, N_DATA)
    parts = [STF, LTF_SYM, LTF_SYM]
    for m in range(NSYM):
        f = np.zeros(NSC, dtype=complex)
        for i, sc in enumerate(DATA_REL):
            f[sc % NSC] = syms[m, i]
        for p in PILOT_REL:
            f[p % NSC] = 1.0 + 0j
        td = np.fft.ifft(f) * np.sqrt(NSC)
        parts.append(np.concatenate([td[-NCP:], td]).astype(np.complex64))
    return np.concatenate(parts).astype(np.complex64)


def init_params(fs, phase1b=False, text="STEI"):
    """Initialize FS-dependent globals. Called once per FS (main + each TX worker).

    phase1b=True akan additionally pre-build N_FRAMES_CYCLE V2I packet frames.
    """
    global FS, BW, MAX_CFO
    global KNOWN_BITS, STF, X_STF, LTF_SYM, X_LTF, LTF_TD, LTF_TD_NO_CP, TX_FRAME
    global PHASE1B_MODE, TX_TEXT, SENSE_KNOWN_BITS, TX_FRAMES_P1B

    FS = float(fs)
    BW = FS  # Effective signal BW = sample rate (B210 sets analog BW = fs)

    # CFO threshold: scale dengan subcarrier spacing.
    # Schmidl-Cox unambiguous = fs/N. B210 ±2 ppm @ 5.9 GHz = ±11.8 kHz absolute.
    # Threshold = min(30 kHz, 30% subcarrier spacing) untuk avoid edge ambiguity & ICI.
    sc_spacing = FS / NSC
    MAX_CFO = min(30000.0, 0.3 * sc_spacing)

    np.random.seed(42)
    KNOWN_BITS = np.random.randint(0, 2, NSYM * N_DATA).astype(np.uint8)

    STF, X_STF = _build_stf()
    LTF_SYM, X_LTF, LTF_TD = _build_ltf()
    LTF_TD_NO_CP = LTF_TD.astype(np.complex64)

    raw = _build_frame_signal(KNOWN_BITS)
    max_amp = float(np.max(np.abs(raw)))
    TX_FRAME = (raw / max_amp * 0.95).astype(np.complex64)
    assert len(TX_FRAME) == FRAME_LEN

    PHASE1B_MODE = bool(phase1b)
    TX_TEXT = text
    if PHASE1B_MODE:
        TX_FRAMES_P1B, SENSE_KNOWN_BITS = _build_all_phase1b_frames(text)
    else:
        TX_FRAMES_P1B, SENSE_KNOWN_BITS = None, None


# ═══════════════════════════════════════════════════════════════════
# SYNCHRONIZATION + DEMOD (sama logic dengan Phase 0, tinggal pakai globals)
# ═══════════════════════════════════════════════════════════════════
def schmidl_cox_metric(buf):
    L = len(buf) - N_HALF
    if L < N_HALF:
        return None, None
    mult = np.conj(buf[:L]) * buf[N_HALF:N_HALF + L]
    pwr  = np.abs(buf[N_HALF:N_HALF + L]) ** 2
    P_cum = np.cumsum(mult)
    R_cum = np.cumsum(pwr)
    nw = L - N_HALF + 1
    P = P_cum[N_HALF - 1:N_HALF - 1 + nw].copy()
    R = R_cum[N_HALF - 1:N_HALF - 1 + nw].copy()
    if nw > 1:
        P[1:] -= P_cum[:nw - 1]
        R[1:] -= R_cum[:nw - 1]
    M = (np.abs(P) ** 2) / (R ** 2 + 1e-18)
    return M, P


def find_all_plateaus(M, threshold=0.7, min_width=4):
    above = M > threshold
    if not np.any(above):
        return []
    edges = np.diff(np.concatenate([[0], above.astype(int), [0]]))
    starts = np.where(edges == 1)[0]
    ends   = np.where(edges == -1)[0]
    plateaus = []
    for s, e in zip(starts, ends):
        w = int(e - s)
        if w < min_width:
            continue
        c = int((s + e - 1) // 2)
        plateaus.append((c, w, float(M[c])))
    plateaus.sort(key=lambda t: t[1] * t[2], reverse=True)
    return plateaus


def validate_sync(Y_eq_pilots, evm_db, cfo_hz, max_cfo=None):
    if max_cfo is None:
        max_cfo = MAX_CFO
    if abs(cfo_hz) > max_cfo:
        return False, f"CFO out of range ({cfo_hz:.0f} Hz)"
    if evm_db > 0:
        return False, f"EVM too high ({evm_db:.1f} dB)"
    p = np.asarray(Y_eq_pilots)
    pilot_mean = np.mean(p)
    pilot_std  = np.std(p)
    pilot_snr_db = 20 * np.log10(np.abs(pilot_mean) / (pilot_std + 1e-9))
    if pilot_snr_db < 3:
        return False, f"pilot scatter (pSNR={pilot_snr_db:.1f} dB)"
    return True, "ok"


def fine_timing_ltf(buf, coarse_start, search_radius=24):
    expected_ltf1 = coarse_start + LTF1_OFF + NCP
    lo = max(0, expected_ltf1 - search_radius)
    hi = expected_ltf1 + search_radius + len(LTF_TD_NO_CP)
    if hi > len(buf):
        return None
    seg = buf[lo:hi]
    mf = np.conj(LTF_TD_NO_CP[::-1])
    corr = np.convolve(seg, mf, mode='valid')
    if len(corr) == 0:
        return None
    peak_local = int(np.argmax(np.abs(corr)))
    refined_ltf1_start = lo + peak_local
    return refined_ltf1_start - (LTF1_OFF + NCP)


def _demod_at(buf_corr, frame_start):
    """Phase 0/1A demod: BPSK seluruh DATA_REL.
    Returns (bits, evm_db, sample_eq, pilot_eq, avg_amp, H_est)."""
    if frame_start < 0 or frame_start + FRAME_LEN > len(buf_corr):
        return None
    fb = buf_corr[frame_start:frame_start + FRAME_LEN]
    avg_amp = float(np.mean(np.abs(fb)))
    Y_ltf1 = np.fft.fft(fb[LTF1_OFF + NCP:LTF1_OFF + NCP + NSC]) / np.sqrt(NSC)
    Y_ltf2 = np.fft.fft(fb[LTF2_OFF + NCP:LTF2_OFF + NCP + NSC]) / np.sqrt(NSC)
    Y_ltf  = 0.5 * (Y_ltf1 + Y_ltf2)
    H_est = np.ones(NSC, dtype=complex)
    active_mask = np.abs(X_LTF) > 1e-9
    H_est[active_mask] = Y_ltf[active_mask] / X_LTF[active_mask]

    decoded = np.empty(NSYM * N_DATA, dtype=np.uint8)
    sample_eq_first = None
    evm_acc = 0.0
    evm_n = 0
    pilots_collected = []

    for m in range(NSYM):
        s0 = DATA_OFF + m * LSYM
        td = fb[s0 + NCP:s0 + LSYM]
        Y  = np.fft.fft(td) / np.sqrt(NSC)
        Y_eq = np.where(np.abs(H_est) > 1e-9, Y / H_est, Y)
        pilot_vals = np.array([Y_eq[p % NSC] for p in PILOT_REL])
        cpe = np.angle(np.mean(pilot_vals))
        Y_eq *= np.exp(-1j * cpe)
        pilots_collected.extend([Y_eq[p % NSC] for p in PILOT_REL])
        data_eq = np.array([Y_eq[sc % NSC] for sc in DATA_REL])
        bits = (np.real(data_eq) < 0).astype(np.uint8)
        decoded[m * N_DATA:(m + 1) * N_DATA] = bits
        ideal = np.where(bits == 0, 1.0 + 0j, -1.0 + 0j)
        evm_acc += float(np.sum(np.abs(data_eq - ideal) ** 2))
        evm_n   += N_DATA
        if m == 0:
            sample_eq_first = data_eq[:10].copy()

    evm_rms = np.sqrt(evm_acc / max(evm_n, 1))
    evm_db  = 20 * np.log10(evm_rms + 1e-12)
    return decoded, evm_db, sample_eq_first, pilots_collected, avg_amp, H_est


def _demod_at_phase1b(buf_corr, frame_start):
    """Phase 1B demod: split SC. QPSK comm + BPSK sense + H_est utk Phase 1C.

    Returns dict atau None.
    """
    if frame_start < 0 or frame_start + FRAME_LEN > len(buf_corr):
        return None
    fb = buf_corr[frame_start:frame_start + FRAME_LEN]
    avg_amp = float(np.mean(np.abs(fb)))

    # Channel estimation
    Y_ltf1 = np.fft.fft(fb[LTF1_OFF + NCP:LTF1_OFF + NCP + NSC]) / np.sqrt(NSC)
    Y_ltf2 = np.fft.fft(fb[LTF2_OFF + NCP:LTF2_OFF + NCP + NSC]) / np.sqrt(NSC)
    Y_ltf  = 0.5 * (Y_ltf1 + Y_ltf2)
    H_est = np.ones(NSC, dtype=complex)
    active_mask = np.abs(X_LTF) > 1e-9
    H_est[active_mask] = Y_ltf[active_mask] / X_LTF[active_mask]

    comm_bits  = np.empty(COMM_BITS_PER_FRAME, dtype=np.uint8)   # 416
    sense_bits = np.empty(SENSE_BITS_PER_FRAME, dtype=np.uint8)  # 988
    pilots_collected = []
    comm_evm_acc = 0.0
    comm_evm_n = 0
    sense_evm_acc = 0.0
    sense_evm_n = 0
    comm_sample_eq = None

    for m in range(NSYM):
        s0 = DATA_OFF + m * LSYM
        td = fb[s0 + NCP:s0 + LSYM]
        Y = np.fft.fft(td) / np.sqrt(NSC)
        Y_eq = np.where(np.abs(H_est) > 1e-9, Y / H_est, Y)
        # Pilot CPE
        pilot_vals = np.array([Y_eq[p % NSC] for p in PILOT_REL])
        cpe = np.angle(np.mean(pilot_vals))
        Y_eq *= np.exp(-1j * cpe)
        pilots_collected.extend([Y_eq[p % NSC] for p in PILOT_REL])

        # COMM: QPSK demap dari 8 center SC
        comm_eq = np.array([Y_eq[sc % NSC] for sc in COMM_SC])
        c_pairs = qpsk_demap(comm_eq)
        comm_bits[m * N_COMM * 2:(m + 1) * N_COMM * 2] = c_pairs.flatten()
        # Comm EVM (vs ideal QPSK)
        ideal_q = qpsk_map(c_pairs)
        comm_evm_acc += float(np.sum(np.abs(comm_eq - ideal_q) ** 2))
        comm_evm_n   += N_COMM

        # SENSE: BPSK hard slice
        sense_eq = np.array([Y_eq[sc % NSC] for sc in SENSE_SC])
        s_bits = (np.real(sense_eq) < 0).astype(np.uint8)
        sense_bits[m * N_SENSE:(m + 1) * N_SENSE] = s_bits
        ideal_s = np.where(s_bits == 0, 1.0 + 0j, -1.0 + 0j)
        sense_evm_acc += float(np.sum(np.abs(sense_eq - ideal_s) ** 2))
        sense_evm_n   += N_SENSE

        if m == 0:
            comm_sample_eq = comm_eq[:8].copy()

    comm_evm_rms = np.sqrt(comm_evm_acc / max(comm_evm_n, 1))
    comm_evm_db  = 20 * np.log10(comm_evm_rms + 1e-12)
    sense_evm_rms = np.sqrt(sense_evm_acc / max(sense_evm_n, 1))
    sense_evm_db  = 20 * np.log10(sense_evm_rms + 1e-12)

    return {
        "comm_bits": comm_bits, "sense_bits": sense_bits,
        "comm_evm_db": comm_evm_db, "sense_evm_db": sense_evm_db,
        "comm_sample_eq": comm_sample_eq,
        "pilots": pilots_collected, "avg_amp": avg_amp,
        "H_est": H_est,
    }


def sync_and_demod(buf, sc_threshold=0.7, max_candidates=3,
                   echo_threshold_db=6.0):
    """Full pipeline. Returns dict including H_est for sensing extension.

    echo_threshold_db: passed ke estimate_ranges (Phase 1C).
    """
    if len(buf) < FRAME_LEN + 64:
        return {"bits": None, "consume": 0}

    search_len = min(len(buf), FRAME_LEN + 256)
    M, P = schmidl_cox_metric(buf[:search_len])
    if M is None:
        return {"bits": None, "consume": 0}

    plateaus = find_all_plateaus(M, threshold=sc_threshold, min_width=3)
    if not plateaus:
        return {"bits": None, "consume": min(LSYM, len(buf) - FRAME_LEN)}

    candidates = plateaus[:max_candidates]
    best_invalid = None

    for center, width, m_peak in candidates:
        epsilon = np.angle(P[center]) / (2 * np.pi * N_HALF)
        frac_cfo_rad = 2 * np.pi * epsilon
        cfo_hz = epsilon * FS

        if abs(cfo_hz) > MAX_CFO:
            continue

        n_idx = np.arange(len(buf))
        buf_corr = buf * np.exp(-1j * frac_cfo_rad * n_idx)
        coarse_frame_start = center - NCP // 2

        frame_start = fine_timing_ltf(buf_corr, coarse_frame_start, search_radius=24)
        if frame_start is None or frame_start < 0:
            continue
        if frame_start + FRAME_LEN > len(buf):
            continue

        out = _demod_at(buf_corr, frame_start)
        if out is None:
            continue
        bits, evm_db, sample_eq, pilot_eq, avg_amp, H_est = out

        valid, reason = validate_sync(pilot_eq, evm_db, cfo_hz)
        result = {
            "bits": bits, "cfo_hz": float(cfo_hz),
            "frame_start": int(frame_start),
            "consume": int(frame_start + FRAME_LEN),
            "avg_amp": avg_amp, "evm_db": float(evm_db),
            "sample_eq": sample_eq, "m_peak": float(m_peak),
            "plateau_w": int(width), "valid": valid, "reason": reason,
            "H_est": H_est,  # untuk Phase 1C (CIR/range estimation)
        }
        # ── Phase 1B: extract comm + sense bila aktif ────────────
        if PHASE1B_MODE:
            p1b = _demod_at_phase1b(buf_corr, frame_start)
            if p1b is not None:
                # Decode V2I packet (first 56 bit)
                text, ctr, crc_ok = decode_packet(p1b["comm_bits"][:PKT_TOTAL_BITS])
                # Sense BER vs known
                if SENSE_KNOWN_BITS is not None:
                    n_sb = min(len(p1b["sense_bits"]), len(SENSE_KNOWN_BITS))
                    sense_ber = float(np.sum(
                        p1b["sense_bits"][:n_sb] != SENSE_KNOWN_BITS[:n_sb]
                    )) / max(n_sb, 1)
                else:
                    sense_ber = float('nan')
                # Range estimation (Phase 1C)
                cir_mag, bin2m = cir_from_h_est(H_est, FS)
                rng = estimate_ranges(cir_mag, bin2m,
                                      snr_threshold_db=echo_threshold_db)
                result.update({
                    "p1b_text": text, "p1b_counter": ctr, "p1b_crc_ok": crc_ok,
                    "p1b_comm_evm_db": float(p1b["comm_evm_db"]),
                    "p1b_sense_evm_db": float(p1b["sense_evm_db"]),
                    "p1b_sense_ber": sense_ber,
                    "p1b_comm_sample_eq": p1b["comm_sample_eq"],
                    "p1c_direct_db": rng["direct_db"],
                    "p1c_noise_floor_db": rng["noise_floor_db"],
                    "p1c_echoes": rng["echoes"],
                    "p1c_n_echoes": rng["n_echoes"],
                    "p1c_cir_mag": cir_mag,    # 256-bin CIR (untuk plot/log)
                    "p1c_bin_to_meter": bin2m,
                })
        if valid:
            return result
        if best_invalid is None:
            best_invalid = result

    if best_invalid is not None:
        best_invalid["consume"] = LSYM
        best_invalid["bits"] = None
        return best_invalid
    return {"bits": None, "consume": LSYM}


def calc_ber(rx_bits):
    n = min(len(rx_bits), len(KNOWN_BITS))
    return float(np.sum(rx_bits[:n] != KNOWN_BITS[:n])) / n if n else 0.5


# ═══════════════════════════════════════════════════════════════════
# AWGN SELF-TEST
# ═══════════════════════════════════════════════════════════════════
def run_simulation(snr_db_list=(0, 5, 10, 15, 20), n_trials=200, cfo_hz=8e3,
                   timing_offset=37):
    print(f"\n{'═' * 60}")
    print(f"  AWGN SELF-TEST   FS={FS/1e6:.1f} MHz  "
          f"(CFO={cfo_hz/1e3:.1f} kHz, timing_offset={timing_offset})")
    print(f"{'═' * 60}")
    print(f"  {'SNR(dB)':>8} | {'BER':>10} | {'EVM(dB)':>8} | {'CFO_err(Hz)':>11}")
    print(f"  {'-' * 50}")

    sig_pwr = float(np.mean(np.abs(TX_FRAME) ** 2))

    for snr_db in snr_db_list:
        ber_acc = 0; ber_n = 0; evm_acc = 0.0
        cfo_err_acc = 0.0; cfo_err_n = 0
        snr_lin = 10 ** (snr_db / 10)
        n0 = sig_pwr / snr_lin

        for _ in range(n_trials):
            pad_pre  = np.zeros(timing_offset, dtype=np.complex64)
            pad_post = np.zeros(256, dtype=np.complex64)
            tx = np.concatenate([pad_pre, TX_FRAME, pad_post])
            n_idx = np.arange(len(tx))
            tx = tx * np.exp(1j * 2 * np.pi * cfo_hz / FS * n_idx)
            noise = (np.random.randn(len(tx)) + 1j * np.random.randn(len(tx)))
            noise = noise.astype(np.complex64) * np.sqrt(n0 / 2)
            rx = (tx + noise).astype(np.complex64)
            res = sync_and_demod(rx)
            if res["bits"] is None:
                continue
            ber_acc += int(np.sum(res["bits"] != KNOWN_BITS))
            ber_n   += len(KNOWN_BITS)
            evm_acc += res["evm_db"]
            cfo_err_acc += abs(res["cfo_hz"] - cfo_hz)
            cfo_err_n   += 1

        if ber_n == 0:
            print(f"  {snr_db:>8.0f} | {'(no sync)':>10} | {'-':>8} | {'-':>11}")
            continue
        ber = ber_acc / ber_n
        evm = evm_acc / max(cfo_err_n, 1)
        cfe = cfo_err_acc / max(cfo_err_n, 1)
        print(f"  {snr_db:>8.0f} | {ber:>10.3e} | {evm:>8.2f} | {cfe:>11.0f}")
    print(f"{'═' * 60}\n")


# ═══════════════════════════════════════════════════════════════════
# USRP INIT & WORKERS
# ═══════════════════════════════════════════════════════════════════
def _init_usrp(serial, fpga, image_dir, is_tx, gain, ant, fs):
    sys.path.append("/usr/local/lib/python3.12/site-packages")
    import uhd
    fpga_suffix = f",fpga={fpga}" if fpga else ""
    strategies = [f"serial={serial}{fpga_suffix}"] if is_tx else \
                 [f"serial={serial}{fpga_suffix}",
                  f"name=LibreSDR_B210mini{fpga_suffix}"]
    old_env = os.environ.get("UHD_IMAGES_DIR")
    os.environ["UHD_IMAGES_DIR"] = image_dir
    usrp = None
    for args in strategies:
        try:
            usrp = uhd.usrp.MultiUSRP(args)
            break
        except RuntimeError:
            time.sleep(1)
    if old_env is None:
        os.environ.pop("UHD_IMAGES_DIR", None)
    else:
        os.environ["UHD_IMAGES_DIR"] = old_env
    if usrp is None:
        raise RuntimeError(f"Device {serial} tidak ditemukan")

    if is_tx:
        usrp.set_tx_rate(fs)
        usrp.set_tx_freq(uhd.libpyuhd.types.tune_request(FC), 0)
        usrp.set_tx_gain(gain, 0)
        usrp.set_tx_antenna(ant, 0)
        usrp.set_tx_bandwidth(fs, 0)
        actual_gain = float(usrp.get_tx_gain(0))
    else:
        usrp.set_rx_rate(fs)
        usrp.set_rx_freq(uhd.libpyuhd.types.tune_request(FC), 0)
        usrp.set_rx_gain(gain, 0)
        usrp.set_rx_antenna(ant, 0)
        usrp.set_rx_bandwidth(fs, 0)
        actual_gain = float(usrp.get_rx_gain(0))

    role = "TX (Lutetia)" if is_tx else "RX (LibreSDR)"
    if abs(actual_gain - gain) > 0.5:
        print(f"[{role}] serial={serial} | gain={gain} dB REQUESTED → "
              f"{actual_gain:.2f} dB CLAMPED (hardware max)")
        print(f"       Max valid: TX≈89.75 dB, RX≈76 dB untuk B210/AD9361")
    else:
        print(f"[{role}] serial={serial} | gain={actual_gain:.2f} dB | "
              f"fc={FC/1e9:.3f} GHz | fs={fs/1e6:.2f} MHz")
    return usrp


def init_tx(fs, gain):
    return _init_usrp(TX_SERIAL, TX_FPGA, TX_IMAGE_DIR, True, gain, TX_ANT, fs)


def init_rx(fs, gain):
    return _init_usrp(RX_SERIAL, RX_FPGA, RX_IMAGE_DIR, False, gain, RX_ANT, fs)


def tx_worker(stop_event, gain_val, fs, frame_delay_s=0.0,
              phase1b=False, text="STEI"):
    """TX worker process. Re-init params for fresh process (mp 'spawn').

    phase1b=True: cycle through TX_FRAMES_P1B (counter 0..255) per send.
    """
    init_params(fs, phase1b=phase1b, text=text)
    try:
        usrp = init_tx(fs, gain_val.value)
    except Exception as e:
        print(f"[TX] FATAL: {e}")
        return
    import uhd
    st_args = uhd.usrp.StreamArgs("fc32", "sc16")
    st_args.args = "num_send_frames=1000"
    st = usrp.get_tx_stream(st_args)
    md = uhd.types.TXMetadata()
    md.start_of_burst = True
    md.end_of_burst   = False

    pad = (-FRAME_LEN) % 1024
    silence_samples = int(frame_delay_s * fs) if frame_delay_s > 0 else 0
    silence_pad = np.zeros(pad + silence_samples, dtype=np.complex64)

    # Pre-build padded frames once (avoid per-loop concat)
    if phase1b:
        padded_frames = [
            np.concatenate([f, silence_pad]).astype(np.complex64)
            for f in TX_FRAMES_P1B
        ]
        print(f"[TX] PHASE 1B: {len(padded_frames)} frame cycle | text='{text}' | "
              f"frame={FRAME_LEN}+{pad}pad samples | dur={FRAME_LEN/fs*1e6:.1f} μs")
    else:
        padded_frames = [np.concatenate([TX_FRAME, silence_pad]).astype(np.complex64)]
        print(f"[TX] frame={FRAME_LEN}+{pad}pad+{silence_samples}silence samples | "
              f"frame_dur={FRAME_LEN/fs*1e6:.1f} μs")

    last_g, ctr = gain_val.value, 0
    n_frames_total = len(padded_frames)

    while not stop_event.is_set():
        if ctr % 50 == 0:
            g = gain_val.value
            if g != last_g:
                usrp.set_tx_gain(g, 0)
                last_g = g
        st.send(padded_frames[ctr % n_frames_total], md)
        md.start_of_burst = False
        ctr += 1
    md.end_of_burst = True
    st.send(np.zeros(256, dtype=np.complex64), md)
    print("[TX] Stop.")


def amplitude_probe(rx_usrp, fs, duration_s=2.0):
    """Capture raw IQ for `duration_s` and report signal stats."""
    import uhd
    st  = rx_usrp.get_rx_stream(uhd.usrp.StreamArgs("fc32", "sc16"))
    cmd = uhd.types.StreamCMD(uhd.types.StreamMode.start_cont)
    cmd.stream_now = True
    st.issue_stream_cmd(cmd)

    chunk = np.zeros(8192, dtype=np.complex64)
    md_rx = uhd.types.RXMetadata()
    samples = []
    target = int(duration_s * fs)
    got = 0
    while got < target:
        n = st.recv(chunk, md_rx)
        if md_rx.error_code != uhd.types.RXMetadataErrorCode.none:
            continue
        samples.append(chunk[:n].copy())
        got += n
    st.issue_stream_cmd(uhd.types.StreamCMD(uhd.types.StreamMode.stop_cont))

    s = np.concatenate(samples)[:target]
    mean_amp = float(np.mean(np.abs(s)))
    rms      = float(np.sqrt(np.mean(np.abs(s) ** 2)))
    peak     = float(np.max(np.abs(s)))
    dc_re    = float(np.mean(np.real(s)))
    dc_im    = float(np.mean(np.imag(s)))
    dc_mag   = np.hypot(dc_re, dc_im)
    crest_db = 20 * np.log10(peak / (rms + 1e-12))

    NF = 4096
    psd = np.zeros(NF)
    n_seg = 0
    for i in range(0, len(s) - NF, NF):
        seg = s[i:i + NF] * np.hanning(NF)
        psd += np.abs(np.fft.fftshift(np.fft.fft(seg))) ** 2
        n_seg += 1
    psd /= max(n_seg, 1)
    psd_db = 10 * np.log10(psd + 1e-18)
    psd_peak = float(np.max(psd_db))
    psd_med  = float(np.median(psd_db))
    spur_dr  = psd_peak - psd_med

    print(f"\n{'─' * 60}")
    print(f"  AMPLITUDE PROBE  (capture {duration_s}s @ {fs/1e6:.1f} MS/s)")
    print(f"{'─' * 60}")
    flag_amp = '⚠ TOO WEAK' if mean_amp < 0.005 else 'ok' if mean_amp < 0.5 else '⚠ NEAR SAT'
    flag_pk  = '⚠ SATURATED' if peak > 0.95 else 'ok'
    flag_dc  = '⚠ HIGH DC' if dc_mag > 0.01 else 'ok'
    flag_sd  = '⚠ NO BAND ENERGY' if spur_dr < 6 else 'signal present' if spur_dr > 15 else 'marginal'
    print(f"  mean|x|     : {mean_amp:.5f}      ({flag_amp})")
    print(f"  RMS         : {rms:.5f}")
    print(f"  peak        : {peak:.5f}      ({flag_pk})")
    print(f"  crest       : {crest_db:.1f} dB  (OFDM expected ≈ 8–12 dB)")
    print(f"  DC offset   : {dc_mag:.5f}      ({flag_dc})")
    print(f"  PSD peak    : {psd_peak:.1f} dB")
    print(f"  PSD median  : {psd_med:.1f} dB")
    print(f"  peak/median : {spur_dr:.1f} dB     ({flag_sd})")
    print(f"{'─' * 60}\n")
    return {"mean_amp": mean_amp, "rms": rms, "peak": peak,
            "dc_mag": dc_mag, "spur_dr": spur_dr}



# ═══════════════════════════════════════════════════════════════════
# PHASE 1D — Real-time visualization (matplotlib)
# ═══════════════════════════════════════════════════════════════════
def setup_live_plot(fs, max_history=200, max_range_m=8.0): #turun dari 80 ke 8, karena ruangan kecil
    """Setup 4-panel live plot. Return handles dict for update.

    Panels:
      (0,0) Range Profile   : CIR magnitude vs delta-range from direct
      (0,1) Constellation   : QPSK comm SC scatter
      (1,0) EVM trend       : Last `max_history` frames
      (1,1) PRR trend + log : Running PRR + last 6 packet text/ctr
    """
    # Let matplotlib auto-select backend (user's desktop env akan punya GUI backend)
    import matplotlib.pyplot as plt

    plt.ion()
    fig, ((ax_range, ax_const), (ax_evm, ax_log)) = plt.subplots(
        2, 2, figsize=(13, 8.5))
    fig.suptitle(f"ISAC Bistatic Live  |  FS={fs/1e6:.1f} MS/s  |  "
                 f"V2I 'STEI' + Range Estimation",
                 fontsize=12, fontweight='bold')

    # ── Range profile ────────────────────────────────────────────
    bin_max_m = 3e8 / (fs * (CIR_NFFT / NSC)) * (CIR_NFFT // 2)
    range_axis = np.arange(CIR_NFFT // 2) * 3e8 / (fs * CIR_NFFT / NSC)
    line_cir, = ax_range.plot(range_axis, np.zeros(CIR_NFFT // 2),
                               'b-', linewidth=1.0, label='CIR (current)')
    line_bg,  = ax_range.plot(range_axis, np.zeros(CIR_NFFT // 2),
                               'g--', linewidth=0.7, alpha=0.6, label='Background median')
    line_sub, = ax_range.plot(range_axis, np.zeros(CIR_NFFT // 2),
                               'r-', linewidth=1.2, alpha=0.8, label='CIR − background')
    direct_marker = ax_range.axvline(0, color='k', linestyle=':', linewidth=0.8, alpha=0.5)
    echo_markers = ax_range.scatter([], [], marker='v', color='red', s=60,
                                     zorder=5, label='detected echoes')
    ax_range.set_xlabel("Delta-range from direct path (m)")
    ax_range.set_ylabel("CIR magnitude (dB)")
    ax_range.set_title("Range Profile (window: Hanning, IFFT zero-pad 4×)")
    ax_range.set_xlim(0, min(max_range_m, bin_max_m))
    ax_range.set_ylim(-80, 5)
    ax_range.grid(True, alpha=0.3)
    ax_range.legend(loc='upper right', fontsize=8)

    # ── Comm constellation ───────────────────────────────────────
    scatter_const = ax_const.scatter([], [], s=20, c='cyan', alpha=0.6, edgecolors='none')
    # Ideal QPSK ref points
    for s in QPSK_TABLE:
        ax_const.plot(s.real, s.imag, '+', color='red', markersize=14,
                      markeredgewidth=2)
    ax_const.set_xlim(-1.6, 1.6)
    ax_const.set_ylim(-1.6, 1.6)
    ax_const.set_aspect('equal')
    ax_const.set_xlabel("In-phase")
    ax_const.set_ylabel("Quadrature")
    ax_const.set_title(f"Comm constellation (QPSK, {N_COMM} center SC)")
    ax_const.grid(True, alpha=0.3)
    ax_const.axhline(0, color='gray', linewidth=0.5)
    ax_const.axvline(0, color='gray', linewidth=0.5)

    # ── EVM trend ────────────────────────────────────────────────
    line_evm_comm,  = ax_evm.plot([], [], 'b-', label='Comm EVM (QPSK)', linewidth=1.0)
    line_evm_sense, = ax_evm.plot([], [], 'g-', label='Sense EVM (BPSK)', linewidth=1.0)
    ax_evm.axhline(-10, color='r', linestyle='--', linewidth=0.8, alpha=0.5,
                   label='QPSK threshold')
    ax_evm.set_xlabel("Frame index")
    ax_evm.set_ylabel("EVM (dB)")
    ax_evm.set_title("EVM trend")
    ax_evm.set_ylim(-25, 5)
    ax_evm.grid(True, alpha=0.3)
    ax_evm.legend(loc='lower right', fontsize=8)

    # ── PRR + Text log ───────────────────────────────────────────
    ax_log.axis('off')
    log_title = ax_log.text(0.02, 0.96, "Recent packets",
                            transform=ax_log.transAxes,
                            fontsize=11, fontweight='bold')
    log_text = ax_log.text(0.02, 0.85, "",
                           transform=ax_log.transAxes,
                           fontsize=9, family='monospace',
                           verticalalignment='top')
    prr_text = ax_log.text(0.02, 0.10, "",
                           transform=ax_log.transAxes,
                           fontsize=11, fontweight='bold',
                           color='green')

    fig.tight_layout(rect=[0, 0, 1, 0.96])
    plt.show(block=False)
    plt.pause(0.1)

    # ── JCAS forward-scattering radar scope (second live window) ─────────
    fig_jcas, (ax_amp, ax_dopp, ax_score) = plt.subplots(3, 1, figsize=(12, 7.5), sharex=True)
    fig_jcas.suptitle("JCAS Forward-Scattering Detector | LoS amplitude + phase/Doppler", fontsize=12, fontweight='bold')
    line_amp_hp, = ax_amp.plot([], [], linewidth=1.1, label='Amplitude HP (dB)')
    ax_amp.axhline(0, linestyle=':', linewidth=0.8, alpha=0.5)
    ax_amp.set_ylabel("Amp HP (dB)")
    ax_amp.set_title("LoS disruption / shadowing metric")
    ax_amp.grid(True, alpha=0.3)
    ax_amp.legend(loc='upper right', fontsize=8)

    line_dopp_hp, = ax_dopp.plot([], [], linewidth=1.1, label='Doppler HP (Hz)')
    ax_dopp.axhline(0, linestyle=':', linewidth=0.8, alpha=0.5)
    ax_dopp.set_ylabel("Doppler HP (Hz)")
    ax_dopp.set_title("Phase-derived Doppler proxy")
    ax_dopp.grid(True, alpha=0.3)
    ax_dopp.legend(loc='upper right', fontsize=8)

    line_score, = ax_score.plot([], [], linewidth=1.2, label='Detection score')
    line_threshold, = ax_score.plot([], [], linestyle='--', linewidth=1.0, label='Adaptive threshold')
    object_markers = ax_score.scatter([], [], s=45, marker='o', label='ADA OBJEK', zorder=5)
    status_text = ax_score.text(0.02, 0.86, "STATUS: CLEAR", transform=ax_score.transAxes,
                                fontsize=13, fontweight='bold')
    ax_score.set_xlabel("Frame index")
    ax_score.set_ylabel("Score")
    ax_score.set_title("Simple CFAR / adaptive threshold")
    ax_score.grid(True, alpha=0.3)
    ax_score.legend(loc='upper right', fontsize=8)
    fig_jcas.tight_layout(rect=[0, 0, 1, 0.95])
    plt.show(block=False)
    plt.pause(0.1)

    return {
        "fig": fig,
        "fig_jcas": fig_jcas,
        "ax_range": ax_range,
        "ax_const": ax_const,
        "ax_evm": ax_evm,
        "ax_log": ax_log,
        "ax_amp": ax_amp,
        "ax_dopp": ax_dopp,
        "ax_score": ax_score,
        "line_cir": line_cir,
        "line_bg": line_bg,
        "line_sub": line_sub,
        "direct_marker": direct_marker,
        "echo_markers": echo_markers,
        "scatter_const": scatter_const,
        "line_evm_comm": line_evm_comm,
        "line_evm_sense": line_evm_sense,
        "line_amp_hp": line_amp_hp,
        "line_dopp_hp": line_dopp_hp,
        "line_score": line_score,
        "line_threshold": line_threshold,
        "object_markers": object_markers,
        "status_text": status_text,
        "log_text": log_text,
        "prr_text": prr_text,
        "range_axis": range_axis,
        "max_history": max_history,
    }


def update_live_plot(handles, state, last_n_displayed):
    """Update plot from latest state. Returns # new frames shown.
    state["cir_history"] = deque of dicts: {cir_mag, bin_to_meter, frame_idx}.
    """
    import matplotlib.pyplot as plt

    if not state["p1b_records"]:
        return last_n_displayed
    cur_n = len(state["p1b_records"])
    if cur_n == last_n_displayed:
        return cur_n

    rec_latest = state["p1b_records"][-1]

    # ── Range profile (latest CIR + background subtraction) ─────
    if state["cir_history"]:
        latest = state["cir_history"][-1]
        cir = np.maximum(latest["cir_mag"], 0)  # defensive: ensure non-negative
        bin2m = latest["bin_to_meter"]

        # Take half (forward range only)
        n_half = len(cir) // 2
        cir_fwd = cir[:n_half]
        cir_db = 20 * np.log10(cir_fwd + 1e-12)
        # Find direct path bin in this CIR
        direct_bin = int(np.argmax(cir_fwd))
        direct_db = float(cir_db[direct_bin])
        # Shift: delta-range = 0 at direct path
        cir_shifted = cir_fwd[direct_bin:]
        n_show = len(cir_shifted)
        delta_range_x = np.arange(n_show) * bin2m
        cir_db_shifted = 20 * np.log10(cir_shifted + 1e-12) - direct_db

        handles["line_cir"].set_data(delta_range_x, cir_db_shifted)

        # Background: median of last K CIRs
        if len(state["cir_history"]) >= 10:
            cir_stack = np.stack([np.maximum(h["cir_mag"][:n_half], 0)
                                  for h in state["cir_history"]])
            bg = np.median(cir_stack, axis=0)
            bg_shifted = bg[direct_bin:]
            bg_db_shifted = 20 * np.log10(bg_shifted + 1e-12) - direct_db
            handles["line_bg"].set_data(delta_range_x, bg_db_shifted)
            # Subtracted (positive part only)
            sub = np.maximum(cir_shifted - bg_shifted, 1e-12)
            sub_db_shifted = 20 * np.log10(sub) - direct_db
            handles["line_sub"].set_data(delta_range_x, sub_db_shifted)
        else:
            handles["line_bg"].set_data([], [])
            handles["line_sub"].set_data([], [])

        # Echo markers from latest record (already delta-range)
        if rec_latest["echoes"]:
            ranges = [e[0] for e in rec_latest["echoes"]]
            db_vals = [e[1] - rec_latest["direct_db"] for e in rec_latest["echoes"]]
            handles["echo_markers"].set_offsets(np.column_stack([ranges, db_vals]))
        else:
            handles["echo_markers"].set_offsets(np.empty((0, 2)))

    # ── Constellation ────────────────────────────────────────────
    if "comm_sample_eq" in rec_latest and rec_latest.get("comm_sample_eq") is not None:
        pts = rec_latest["comm_sample_eq"]
        handles["scatter_const"].set_offsets(
            np.column_stack([np.real(pts), np.imag(pts)]))

    # ── EVM trend ────────────────────────────────────────────────
    n_show = min(handles["max_history"], cur_n)
    recs_show = state["p1b_records"][-n_show:]
    xs = [r["frame_idx"] for r in recs_show]
    comm_y = [r["comm_evm_db"] for r in recs_show]
    sense_y = [r["sense_evm_db"] for r in recs_show]
    handles["line_evm_comm"].set_data(xs, comm_y)
    handles["line_evm_sense"].set_data(xs, sense_y)
    if xs:
        handles["ax_evm"].set_xlim(min(xs), max(xs) + 1)

    # ── PRR + log ───────────────────────────────────────────────
    n_total = len(state["p1b_records"])
    n_ok = state["p1b_crc_ok_count"]
    prr = n_ok / max(n_total, 1) * 100
    handles["prr_text"].set_text(
        f"PRR: {prr:.2f}%  ({n_ok}/{n_total} packets)\n"
        f"OVF: {state['ovf']}  |  Recent comm EVM: {rec_latest['comm_evm_db']:.2f} dB"
    )

    # Last 8 packets
    recent = state["p1b_records"][-8:]
    log_lines = [f"{'frm':>4} {'text':>4} {'ctr':>4} {'CRC':>4} "
                 f"{'commEVM':>8} {'echoes':<14}"]
    log_lines.append("-" * 50)
    for r in recent:
        crc_str = "OK" if r["crc_ok"] else "FAIL"
        echo_str = ", ".join(f"{e[0]:.0f}m" for e in r["echoes"][:3]) or "-"
        log_lines.append(
            f"{r['frame_idx']:>4} {r['text']:>4} {r['counter']:>4} {crc_str:>4} "
            f"{r['comm_evm_db']:>8.2f} {echo_str:<14}"
        )
    handles["log_text"].set_text("\n".join(log_lines))

    handles["fig"].canvas.draw_idle()
    handles["fig"].canvas.flush_events()
    return cur_n


def setup_jcas_comm_plot(fs, max_history=180, sweep_frames=90):
    """Very simple demo monitor: tinggal lihat CLEAR / ADA OBJEK.

    Panel dibuat untuk orang awam:
      - kiri: status super besar,
      - kanan: radar sweep visual,
      - bawah: confidence meter + score/threshold sebagai bukti teknis.
    """
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    plt.ion()
    fig = plt.figure(figsize=(13.5, 8.2))
    fig.suptitle(
        f"JCAS Object Monitor — USRP TX ↔ RX Line-of-Sight | FS={fs/1e6:.1f} MS/s",
        fontsize=15, fontweight='bold', color='white'
    )
    fig.patch.set_facecolor('#0b0b0b')

    gs = fig.add_gridspec(3, 2, width_ratios=[1.25, 1.0], height_ratios=[2.25, 0.42, 1.0])
    ax_status = fig.add_subplot(gs[0:2, 0])
    ax_radar = fig.add_subplot(gs[0:2, 1], projection='polar')
    ax_meter = fig.add_subplot(gs[2, 0])
    ax_score = fig.add_subplot(gs[2, 1])

    # ── Big status panel ─────────────────────────────────────────
    ax_status.set_axis_off()
    ax_status.set_facecolor('#001b0a')
    status_box = Rectangle((0.02, 0.04), 0.96, 0.92, transform=ax_status.transAxes,
                           facecolor='#003d16', edgecolor='#00ff66', linewidth=5)
    ax_status.add_patch(status_box)
    status_text = ax_status.text(
        0.5, 0.68, 'KALIBRASI', transform=ax_status.transAxes,
        ha='center', va='center', fontsize=48, fontweight='bold', color='#eeeeee'
    )
    subtitle_text = ax_status.text(
        0.5, 0.49, 'Diamkan area TX–RX sebentar untuk ambil baseline',
        transform=ax_status.transAxes, ha='center', va='center',
        fontsize=16, fontweight='bold', color='#eeeeee'
    )
    object_text = ax_status.text(
        0.5, 0.33, 'TX  ━━━━━━━━━━━━━  LoS  ━━━━━━━━━━━━━  RX',
        transform=ax_status.transAxes, ha='center', va='center',
        fontsize=16, family='monospace', color='#eaffef'
    )
    metrics_text = ax_status.text(
        0.5, 0.15, 'Frame: -   Confidence: -   BER: -   EVM: -',
        transform=ax_status.transAxes, ha='center', va='center',
        fontsize=12, family='monospace', color='#d0d0d0'
    )

    # ── Realtime communication text from TX to RX ────────────────
    # Ini khusus untuk Phase 1B: payload teks dari TX akan muncul langsung
    # di panel radar, misalnya: TX → RX MESSAGE: "STEI" | CTR=104 | CRC OK
    comm_text = ax_status.text(
        0.5, 0.245, 'TX → RX MESSAGE: waiting...',
        transform=ax_status.transAxes, ha='center', va='center',
        fontsize=18, fontweight='bold', family='monospace', color='#00ffcc'
    )
    comm_history_text = ax_status.text(
        0.5, 0.075, 'Recent RX packets: -',
        transform=ax_status.transAxes, ha='center', va='center',
        fontsize=10, family='monospace', color='#d7ffe7'
    )

    # ── Radar scope visual ───────────────────────────────────────
    ax_radar.set_facecolor('#00140b')
    ax_radar.set_theta_zero_location('N')
    ax_radar.set_theta_direction(-1)
    ax_radar.set_rlim(0, 10)
    ax_radar.set_rticks([2, 4, 6, 8, 10])
    ax_radar.set_yticklabels(['', '', '', '', ''])
    ax_radar.set_thetagrids(range(0, 360, 45), labels=[''] * 8)
    ax_radar.grid(True, color='#00aa44', alpha=0.45, linewidth=0.9)
    ax_radar.spines['polar'].set_color('#00ff66')
    ax_radar.set_title('Radar visual', color='#8cffb2', pad=16, fontweight='bold')

    theta = np.linspace(0, 2*np.pi, 361)
    for rr in [2, 4, 6, 8, 10]:
        ax_radar.plot(theta, np.full_like(theta, rr), color='#00aa44', alpha=0.22, linewidth=0.8)
    for deg in range(0, 360, 45):
        t = np.deg2rad(deg)
        ax_radar.plot([t, t], [0, 10], color='#00aa44', alpha=0.18, linewidth=0.8)

    sweep_line, = ax_radar.plot([0, 0], [0, 10], color='#00ff66', linewidth=3.0, alpha=0.95)
    sweep_tail_lines = []
    for i in range(1, 7):
        ln, = ax_radar.plot([0, 0], [0, 10], color='#00ff66', linewidth=2.0, alpha=max(0.06, 0.45 - i*0.055))
        sweep_tail_lines.append(ln)
    blip_scatter = ax_radar.scatter([], [], s=[], c='#ff2222', edgecolors='#ffff66', linewidths=1.4, alpha=0.98, zorder=6)
    ax_radar.scatter([0], [0], s=95, c='#00ff66', edgecolors='#ccffdd', linewidths=1.0, zorder=7)
    radar_status_text = ax_radar.text(
        0.5, -0.08, 'Tidak ada blip = tidak ada objek', transform=ax_radar.transAxes,
        ha='center', va='top', color='#d8ffe4', fontsize=11, fontweight='bold'
    )

    # ── Confidence meter ─────────────────────────────────────────
    ax_meter.set_facecolor('#181818')
    ax_meter.set_xlim(0, 100)
    ax_meter.set_ylim(0, 1)
    ax_meter.set_yticks([])
    ax_meter.set_xticks([0, 25, 50, 75, 100])
    ax_meter.tick_params(colors='white')
    ax_meter.set_title('Object confidence meter', color='white', fontsize=11, fontweight='bold')
    ax_meter.grid(True, axis='x', alpha=0.18)
    meter_bg = Rectangle((0, 0.25), 100, 0.5, facecolor='#303030', edgecolor='#aaaaaa', linewidth=1.2)
    meter_fill = Rectangle((0, 0.25), 0, 0.5, facecolor='#00cc55', edgecolor='none')
    ax_meter.add_patch(meter_bg)
    ax_meter.add_patch(meter_fill)
    meter_line = ax_meter.axvline(100, color='#ffcc00', linestyle='--', linewidth=2)
    meter_text = ax_meter.text(50, 0.50, 'Baseline...', ha='center', va='center', color='white', fontsize=13, fontweight='bold')
    for spine in ax_meter.spines.values():
        spine.set_color('#aaaaaa')

    # ── Small technical score plot ───────────────────────────────
    ax_score.set_facecolor('#181818')
    line_score, = ax_score.plot([], [], linewidth=1.8, label='Score')
    line_thr, = ax_score.plot([], [], linestyle='--', linewidth=1.4, label='Threshold')
    detect_scatter = ax_score.scatter([], [], s=55, marker='o', label='Object')
    ax_score.set_title('Score vs threshold', color='white', fontsize=11, fontweight='bold')
    ax_score.set_xlabel('Frame', color='white')
    ax_score.set_ylabel('Score', color='white')
    ax_score.tick_params(colors='white')
    ax_score.grid(True, alpha=0.25)
    ax_score.legend(loc='upper right', fontsize=8)
    for spine in ax_score.spines.values():
        spine.set_color('#aaaaaa')

    fig.tight_layout(rect=[0, 0, 1, 0.94])
    plt.show(block=False)
    plt.pause(0.1)

    return {
        'fig': fig,
        'ax_status': ax_status, 'ax_radar': ax_radar, 'ax_meter': ax_meter, 'ax_score': ax_score,
        'status_box': status_box, 'status_text': status_text, 'subtitle_text': subtitle_text,
        'object_text': object_text, 'metrics_text': metrics_text,
        'comm_text': comm_text, 'comm_history_text': comm_history_text,
        'radar_status_text': radar_status_text,
        'sweep_line': sweep_line, 'sweep_tail_lines': sweep_tail_lines, 'blip_scatter': blip_scatter,
        'meter_fill': meter_fill, 'meter_text': meter_text, 'meter_line': meter_line,
        'line_score': line_score, 'line_thr': line_thr, 'detect_scatter': detect_scatter,
        'max_history': max_history, 'sweep_frames': sweep_frames,
    }


def update_jcas_comm_plot(handles, state, last_n_displayed):
    """Update one-look object/no-object monitor."""
    if state['n'] == last_n_displayed:
        return last_n_displayed

    cur_n = state['n']
    recs_all = state.get('jcas_records', [])
    latest = recs_all[-1] if recs_all else {}

    status = latest.get('jcas_status', 'KALIBRASI')
    score_latest = float(latest.get('jcas_score', 0.0))
    thr_latest = float(latest.get('jcas_threshold', 0.0))
    ratio_latest = float(latest.get('jcas_score_ratio', score_latest / max(thr_latest, 1e-9))) if thr_latest > 0 else 0.0
    baseline_ready = bool(latest.get('jcas_baseline_ready', False))
    baseline_progress = float(latest.get('jcas_baseline_progress', 0.0))
    ber_latest = state['bers'][-1] if state.get('bers') else float('nan')
    evm_latest = state['evms'][-1] if state.get('evms') else float('nan')
    ovf = state.get('ovf', 0)

    # Confidence: below threshold <100%, above threshold >100%.
    confidence = float(np.clip(ratio_latest * 100.0, 0.0, 180.0)) if baseline_ready else baseline_progress * 100.0
    fill_width = min(confidence, 100.0)
    handles['meter_fill'].set_width(fill_width)

    if status == 'ADA_OBJEK':
        handles['status_box'].set_facecolor('#650000')
        handles['status_box'].set_edgecolor('#ff3333')
        handles['status_text'].set_text('ADA\nOBJEK')
        handles['status_text'].set_color('#ffdddd')
        handles['subtitle_text'].set_text('Objek sedang mengganggu jalur LoS TX–RX')
        handles['subtitle_text'].set_color('#fff0f0')
        handles['object_text'].set_text('TX  ━━━━━━━━  ⚠ OBJEK  ⚠  ━━━━━━━━  RX')
        handles['object_text'].set_color('#fff0aa')
        handles['meter_fill'].set_facecolor('#ff3333')
        handles['meter_text'].set_text(f'DETECTED  | confidence {confidence:.0f}%')
        handles['meter_text'].set_color('#fff0f0')
        handles['radar_status_text'].set_text('BLIP MERAH = OBJEK TERDETEKSI')
        handles['radar_status_text'].set_color('#ffdddd')
    elif status == 'CLEAR':
        handles['status_box'].set_facecolor('#003d16')
        handles['status_box'].set_edgecolor('#00ff66')
        handles['status_text'].set_text('CLEAR')
        handles['status_text'].set_color('#8cffb2')
        handles['subtitle_text'].set_text('Tidak ada benda terdeteksi di antara TX dan RX')
        handles['subtitle_text'].set_color('#d8ffe4')
        handles['object_text'].set_text('TX  ━━━━━━━━━━━━━  LoS  ━━━━━━━━━━━━━  RX')
        handles['object_text'].set_color('#eaffef')
        handles['meter_fill'].set_facecolor('#00cc55')
        handles['meter_text'].set_text(f'CLEAR  | confidence {confidence:.0f}% dari threshold')
        handles['meter_text'].set_color('#ffffff')
        handles['radar_status_text'].set_text('Tidak ada blip = tidak ada objek')
        handles['radar_status_text'].set_color('#d8ffe4')
    else:
        handles['status_box'].set_facecolor('#3a2b00')
        handles['status_box'].set_edgecolor('#ffcc00')
        handles['status_text'].set_text('KALIBRASI')
        handles['status_text'].set_color('#fff0aa')
        handles['subtitle_text'].set_text('Jangan lewatkan objek dulu — sistem sedang ambil baseline kosong')
        handles['subtitle_text'].set_color('#fff0aa')
        handles['object_text'].set_text('TX  ━━━━━━━━━━━━━  LoS  ━━━━━━━━━━━━━  RX')
        handles['object_text'].set_color('#ffffff')
        handles['meter_fill'].set_facecolor('#ffcc00')
        handles['meter_text'].set_text(f'Baseline {baseline_progress*100:.0f}%')
        handles['meter_text'].set_color('#ffffff')
        handles['radar_status_text'].set_text('Kalibrasi baseline kosong...')
        handles['radar_status_text'].set_color('#fff0aa')

    handles['metrics_text'].set_text(
        f"Frame: {cur_n}   Score: {score_latest:.2f}/{thr_latest:.2f}   "
        f"BER: {ber_latest:.4f}   EVM: {evm_latest:.2f} dB   OVF: {ovf}"
    )

    # Radar sweep + held blip.
    sweep_frames = handles.get('sweep_frames', 90)
    theta_now = 2 * np.pi * ((cur_n % sweep_frames) / float(sweep_frames))
    handles['sweep_line'].set_data([theta_now, theta_now], [0, 10])
    for i, ln in enumerate(handles['sweep_tail_lines'], start=1):
        theta_tail = theta_now - 2 * np.pi * (i / float(sweep_frames))
        ln.set_data([theta_tail, theta_tail], [0, 10])

    hold_frames = 35
    recent_det = [r for r in recs_all[-hold_frames:] if r.get('jcas_object_detected', False)]
    if recent_det:
        r_last = recent_det[-1]
        rel = float(r_last.get('jcas_score_ratio', 1.0))
        rel = float(np.clip(rel, 0.8, 4.0))
        theta_blip = 2*np.pi*((r_last['frame_idx'] % sweep_frames)/float(sweep_frames))
        radius_blip = float(np.clip(4.0 + 1.4 * rel, 4.8, 9.5))
        size = float(np.clip(650 * rel, 600, 1900))
        handles['blip_scatter'].set_offsets(np.array([[theta_blip, radius_blip]]))
        handles['blip_scatter'].set_sizes([size])
    else:
        handles['blip_scatter'].set_offsets(np.empty((0, 2)))
        handles['blip_scatter'].set_sizes([])

    # Score plot.
    n_show = min(handles['max_history'], len(recs_all))
    recs = recs_all[-n_show:]
    if recs:
        xr = [r['frame_idx'] for r in recs]
        score = np.asarray([r.get('jcas_score', 0.0) for r in recs], dtype=float)
        thr = np.asarray([r.get('jcas_threshold', 0.0) for r in recs], dtype=float)
        handles['line_score'].set_data(xr, score)
        handles['line_thr'].set_data(xr, thr)
        det_x = [r['frame_idx'] for r in recs if r.get('jcas_object_detected', False)]
        det_y = [r.get('jcas_score', 0.0) for r in recs if r.get('jcas_object_detected', False)]
        if det_x:
            handles['detect_scatter'].set_offsets(np.column_stack([det_x, det_y]))
        else:
            handles['detect_scatter'].set_offsets(np.empty((0, 2)))
        handles['ax_score'].set_xlim(max(1, xr[0]), max(2, xr[-1] + 1))
        max_y = max(float(np.nanmax(score)) if score.size else 1.0,
                    float(np.nanmax(thr)) if thr.size else 1.0,
                    2.0)
        handles['ax_score'].set_ylim(0, max(5.0, max_y * 1.25))

    # ── Realtime TX → RX communication display ─────────────
    # Tampil bila program dijalankan dengan --phase1b --text "STEI".
    # Data berasal dari packet yang sudah didecode di RX: text, counter, CRC.
    if 'comm_text' in handles and 'comm_history_text' in handles:
        p1b_records = state.get('p1b_records', [])
        if p1b_records:
            last_pkt = p1b_records[-1]
            rx_text = str(last_pkt.get('text', '????')).replace('\x00', '').strip()
            rx_ctr = last_pkt.get('counter', '-')
            crc_ok = bool(last_pkt.get('crc_ok', False))
            crc_label = 'CRC OK' if crc_ok else 'CRC FAIL'

            if rx_text == '':
                rx_text = '????'

            handles['comm_text'].set_text(
                f'TX → RX MESSAGE: "{rx_text}"   | CTR={rx_ctr} | {crc_label}'
            )
            handles['comm_text'].set_color('#00ffcc' if crc_ok else '#ff5555')

            recent = p1b_records[-6:]
            recent_lines = []
            for r in recent:
                t = str(r.get('text', '????')).replace('\x00', '').strip() or '????'
                c = r.get('counter', '-')
                ok = 'OK' if r.get('crc_ok', False) else 'FAIL'
                recent_lines.append(f'{t}:{c}:{ok}')
            handles['comm_history_text'].set_text(
                'Recent RX packets: ' + ' | '.join(recent_lines)
            )
        else:
            handles['comm_text'].set_text('TX → RX MESSAGE: waiting for packet...')
            handles['comm_text'].set_color('#cccccc')
            handles['comm_history_text'].set_text('Recent RX packets: -')

    handles['fig'].canvas.draw_idle()
    handles['fig'].canvas.flush_events()
    return cur_n


# ═══════════════════════════════════════════════════════════════════
# HARDWARE LOOP — single FS run, returns metrics dict
# ═══════════════════════════════════════════════════════════════════
def run_hardware(fs, n_frames=100, tx_gain=80.0, rx_gain=70.0,
                 probe_wait=10.0, frame_delay=0.0,
                 dc_offset_auto=True, probe_only=False,
                 warmup_frames=50, verbose=True,
                 phase1b=False, text="STEI", log_csv=None,
                 plot=False, echo_threshold_db=6.0, cir_history_len=30,
                 jcas=True, sense_ma_len=30, sense_cfar_len=160,
                 sense_threshold_k=5.0, sense_min_score=2.2):
    """
    Single-FS hardware run. Returns metrics dict for sweep aggregation.

    phase1b=True: aktifkan V2I "STEI" comm + Phase 1C range estimation.
    plot=True   : aktifkan Phase 1D matplotlib live plot (requires phase1b).
    log_csv     : path untuk per-frame Phase 1B CSV log.
    echo_threshold_db: SNR threshold di atas noise floor (default 6, lower→more sensitive).
    cir_history_len  : length of CIR history buffer untuk background subtraction.
    jcas=True        : aktifkan Forward Scattering / LoS disruption detector.
    sense_ma_len     : panjang moving median untuk DC/clutter removal.
    sense_cfar_len   : panjang rolling reference window untuk threshold adaptif.
    """
    init_params(fs, phase1b=phase1b, text=text)
    import uhd

    tx_gain_sh = mp.Value('d', tx_gain)
    rx_gain_sh = mp.Value('d', rx_gain)
    stop_ev = mp.Event()

    if verbose:
        print(f"\n{'█' * 60}")
        print(f"  RUN @ FS = {fs/1e6:.2f} MS/s  |  TX={tx_gain} dB  RX={rx_gain} dB")
        print(f"  BW eff ≈ {N_ACTIVE * (fs/NSC) / 1e6:.2f} MHz  |  "
              f"δR theoretical ≈ {3e8/(2*N_ACTIVE*(fs/NSC)):.1f} m")
        if phase1b:
            print(f"  PHASE 1B+1C aktif | text='{text}' | "
                  f"comm SC={N_COMM} (QPSK) | sense SC={N_SENSE} (BPSK)")
        print(f"{'█' * 60}")

    print("[INIT] Start TX...")
    tx_proc = mp.Process(target=tx_worker,
                         args=(stop_ev, tx_gain_sh, fs, frame_delay, phase1b, text))
    tx_proc.start()
    time.sleep(3)

    print("[INIT] Start RX...")
    try:
        rx_usrp = init_rx(fs, rx_gain_sh.value)
    except Exception as e:
        print(f"[RX] FATAL: {e}")
        stop_ev.set()
        tx_proc.join(timeout=5)
        return None

    if dc_offset_auto:
        try:
            rx_usrp.set_rx_dc_offset(True, 0)
            rx_usrp.set_rx_iq_balance(True, 0)
            print("[INIT] RX DC + IQ auto-correction: ON")
        except Exception as e:
            print(f"[INIT] DC/IQ not supported: {e}")

    if probe_wait > 0:
        print(f"[INIT] Probe wait {probe_wait}s...")
        time.sleep(probe_wait)

    probe_stats = amplitude_probe(rx_usrp, fs, duration_s=2.0)
    if probe_stats["mean_amp"] < 0.001:
        print("✗ ABORT: signal too low.")
        stop_ev.set()
        tx_proc.join(timeout=5)
        return {
            "fs_mhz": fs/1e6, "tx_gain": tx_gain, "rx_gain": rx_gain,
            "status": "ABORT_LOW_SIGNAL", "n_decoded": 0, "n_target": n_frames,
            "valid_rate": 0.0, "mean_ber": float('nan'), "mean_evm_db": float('nan'),
            "probe_mean_amp": probe_stats["mean_amp"],
            "probe_peak": probe_stats["peak"],
            "ovf_total": 0, "ovf_steady": 0, "ovf_rate_steady": 0.0,
        }
    if probe_only:
        stop_ev.set()
        tx_proc.join(timeout=5)
        return {"fs_mhz": fs/1e6, "status": "PROBE_ONLY",
                "probe_mean_amp": probe_stats["mean_amp"],
                "probe_peak": probe_stats["peak"]}

    state = {"bers": [], "cfos": [], "amps": [], "evms": [],
             "n": 0, "ovf": 0, "ovf_at_warmup_end": None, "running": True,
             # Phase 1B/1C accumulators
             "p1b_records": [],     # per-frame dict
             "p1b_crc_ok_count": 0,
             "p1b_crc_fail_count": 0,
             # Phase 1D plot buffers
             "cir_history": deque(maxlen=cir_history_len),
             # JCAS records are filled in both Phase 0 and Phase 1B modes
             "jcas_records": [],
             }

    # Phase 1F detector: frame-rate based Doppler proxy.
    # Satu frame berdurasi FRAME_LEN/fs detik, sehingga fs_frame = fs/FRAME_LEN.
    jcas_detector = ForwardScatterDetector(
        fs_frame=fs / FRAME_LEN,
        ma_len=sense_ma_len,
        cfar_len=sense_cfar_len,
        threshold_k=sense_threshold_k,
        min_score=sense_min_score,
    ) if jcas else None

    # Live plot: Phase 1B memakai plot ISAC lengkap, mode awal memakai plot komunikasi + JCAS radar.
    plot_handles = None
    plot_mode = None
    if plot:
        if "DISPLAY" not in os.environ and sys.platform == "linux":
            print("[WARN] No DISPLAY env detected (SSH tanpa -X?).")
            print("       Plot membutuhkan X11 display. Disabled.")
            print("       Jika via SSH, jalankan: ssh -X user@host atau jalankan langsung di desktop Ubuntu.")
            print("       Data tetap bisa disimpan dengan --log-csv.")
        else:
            try:
                print("[INIT] Membuka jendela live plot matplotlib...")
                print("       Tkinter canvas init bisa 3-10 detik pertama kali. MOHON SABAR.")
                t0 = time.time()
                # Pakai satu tampilan demo: radar CLEAR/ADA OBJEK + teks komunikasi TX→RX.
                # Dengan --phase1b, payload teks seperti "STEI" tetap didecode dan ditampilkan realtime.
                plot_handles = setup_jcas_comm_plot(fs)
                plot_mode = "jcas_comm"
                print(f"[INIT] Plot window ready ✓ ({time.time()-t0:.1f}s)")
            except Exception as e:
                print(f"[WARN] Plot setup gagal: {e}")
                print("       Continuing without plot. Data tetap tersimpan di log CSV.")
                plot_handles = None
                plot_mode = None

    # CSV writer setup (Phase 1B only)
    csv_file = None
    csv_writer = None
    if phase1b and log_csv:
        csv_file = open(log_csv, "w", newline="")
        csv_writer = csv.writer(csv_file)
        csv_writer.writerow([
            "frame_idx", "text", "counter", "crc_ok",
            "comm_evm_db", "sense_evm_db", "sense_ber",
            "cfo_hz", "amp", "ovf",
            "direct_db", "noise_floor_db", "n_echoes", "echoes_m",
            "jcas_status", "jcas_amp_hp_db", "jcas_phase_hp_rad",
            "jcas_doppler_hz", "jcas_doppler_hp_hz",
            "jcas_score", "jcas_threshold", "jcas_object_detected",
        ])

    def rx_thread():
        st_args = uhd.usrp.StreamArgs("fc32", "sc16")
        st_args.args = "num_recv_frames=1000"
        st  = rx_usrp.get_rx_stream(st_args)
        cmd = uhd.types.StreamCMD(uhd.types.StreamMode.start_cont)
        cmd.stream_now = True
        st.issue_stream_cmd(cmd)

        chunk = np.zeros(FRAME_LEN, dtype=np.complex64)
        md_rx = uhd.types.RXMetadata()
        buf   = np.zeros(0, dtype=np.complex64)
        last_g, ctr = rx_gain_sh.value, 0

        while state["running"]:
            if ctr % 50 == 0:
                g = rx_gain_sh.value
                if g != last_g:
                    rx_usrp.set_rx_gain(g, 0)
                    last_g = g

            nsamp = st.recv(chunk, md_rx)
            if md_rx.error_code == uhd.types.RXMetadataErrorCode.overflow:
                state["ovf"] += 1
                buf = np.zeros(0, dtype=np.complex64)
                ctr += 1
                continue

            buf = np.concatenate([buf, chunk[:nsamp]])

            while len(buf) >= FRAME_LEN + 256:
                res = sync_and_demod(buf, sc_threshold=0.7,
                                     echo_threshold_db=echo_threshold_db)
                if res["bits"] is not None:
                    b = calc_ber(res["bits"])
                    state["bers"].append(b)
                    state["cfos"].append(res["cfo_hz"])
                    state["amps"].append(res["avg_amp"])
                    state["evms"].append(res["evm_db"])
                    state["n"] += 1

                    if state["n"] == warmup_frames:
                        state["ovf_at_warmup_end"] = state["ovf"]

                    # Phase 1F JCAS is calculated for BOTH modes:
                    # - Phase 0/original communication mode: BER/CFO/EVM tetap seperti awal
                    # - Phase 1B mode: V2I packet + CIR/range + JCAS
                    sensing = None
                    if jcas_detector is not None:
                        sensing = jcas_detector.update(
                            avg_amp=res.get("avg_amp"),
                            H_est=res.get("H_est"),
                            cir_mag=res.get("p1c_cir_mag"),
                        )
                    if sensing is not None:
                        state["jcas_records"].append({
                            "frame_idx": state["n"],
                            "jcas_status": sensing["status"],
                            "jcas_amp_db": sensing["amp_db"],
                            "jcas_amp_hp_db": sensing["amp_hp_db"],
                            "jcas_phase_rad": sensing["phase_rad"],
                            "jcas_phase_hp_rad": sensing["phase_hp_rad"],
                            "jcas_doppler_hz": sensing["doppler_hz"],
                            "jcas_doppler_hp_hz": sensing["doppler_hp_hz"],
                            "jcas_score": sensing["score"],
                            "jcas_threshold": sensing["threshold"],
                            "jcas_score_ratio": sensing.get("score_ratio", 0.0),
                            "jcas_baseline_ready": sensing.get("baseline_ready", True),
                            "jcas_baseline_progress": sensing.get("baseline_progress", 1.0),
                            "jcas_object_detected": sensing["object_detected"],
                        })
                    else:
                        state["jcas_records"].append({
                            "frame_idx": state["n"],
                            "jcas_status": "OFF",
                            "jcas_amp_hp_db": 0.0,
                            "jcas_phase_hp_rad": 0.0,
                            "jcas_doppler_hz": 0.0,
                            "jcas_doppler_hp_hz": 0.0,
                            "jcas_score": 0.0,
                            "jcas_threshold": 0.0,
                            "jcas_score_ratio": 0.0,
                            "jcas_baseline_ready": False,
                            "jcas_baseline_progress": 0.0,
                            "jcas_object_detected": False,
                        })

                    # Phase 1B/1C record
                    if phase1b and "p1b_text" in res:
                        rec = {
                            "frame_idx": state["n"],
                            "text": res["p1b_text"],
                            "counter": res["p1b_counter"],
                            "crc_ok": res["p1b_crc_ok"],
                            "comm_evm_db": res["p1b_comm_evm_db"],
                            "sense_evm_db": res["p1b_sense_evm_db"],
                            "sense_ber": res["p1b_sense_ber"],
                            "direct_db": res["p1c_direct_db"],
                            "noise_floor_db": res["p1c_noise_floor_db"],
                            "echoes": res["p1c_echoes"],
                            "n_echoes": res["p1c_n_echoes"],
                            "comm_sample_eq": res.get("p1b_comm_sample_eq"),
                        }
                        if sensing is not None:
                            rec.update({
                                "jcas_status": sensing["status"],
                                "jcas_amp_db": sensing["amp_db"],
                                "jcas_amp_hp_db": sensing["amp_hp_db"],
                                "jcas_phase_rad": sensing["phase_rad"],
                                "jcas_phase_hp_rad": sensing["phase_hp_rad"],
                                "jcas_doppler_hz": sensing["doppler_hz"],
                                "jcas_doppler_hp_hz": sensing["doppler_hp_hz"],
                                "jcas_score": sensing["score"],
                                "jcas_threshold": sensing["threshold"],
                                "jcas_object_detected": sensing["object_detected"],
                            })
                        else:
                            rec.update({
                                "jcas_status": "OFF",
                                "jcas_amp_hp_db": 0.0,
                                "jcas_phase_hp_rad": 0.0,
                                "jcas_doppler_hz": 0.0,
                                "jcas_doppler_hp_hz": 0.0,
                                "jcas_score": 0.0,
                                "jcas_threshold": 0.0,
                                "jcas_object_detected": False,
                            })
                        state["p1b_records"].append(rec)
                        if rec["crc_ok"]:
                            state["p1b_crc_ok_count"] += 1
                        else:
                            state["p1b_crc_fail_count"] += 1
                        # Store CIR untuk Phase 1D plot/background subtraction
                        state["cir_history"].append({
                            "cir_mag": res["p1c_cir_mag"],
                            "bin_to_meter": res["p1c_bin_to_meter"],
                            "frame_idx": state["n"],
                        })
                        if csv_writer is not None:
                            echoes_str = ";".join(
                                f"{r:.2f}@{db:.1f}dB" for r, db in rec["echoes"]
                            )
                            csv_writer.writerow([
                                rec["frame_idx"], rec["text"], rec["counter"],
                                int(rec["crc_ok"]),
                                f"{rec['comm_evm_db']:.2f}",
                                f"{rec['sense_evm_db']:.2f}",
                                f"{rec['sense_ber']:.4f}",
                                f"{res['cfo_hz']:.0f}",
                                f"{res['avg_amp']:.4f}",
                                state["ovf"],
                                f"{rec['direct_db']:.1f}",
                                f"{rec['noise_floor_db']:.1f}",
                                rec["n_echoes"],
                                echoes_str,
                                rec.get("jcas_status", "OFF"),
                                f"{rec.get('jcas_amp_hp_db', 0.0):.3f}",
                                f"{rec.get('jcas_phase_hp_rad', 0.0):.5f}",
                                f"{rec.get('jcas_doppler_hz', 0.0):.3f}",
                                f"{rec.get('jcas_doppler_hp_hz', 0.0):.3f}",
                                f"{rec.get('jcas_score', 0.0):.3f}",
                                f"{rec.get('jcas_threshold', 0.0):.3f}",
                                int(rec.get("jcas_object_detected", False)),
                            ])
                            csv_file.flush()

                    if state["n"] <= 2 and res["sample_eq"] is not None and verbose:
                        print(f"\n--- Frame {state['n']} sample ---")
                        for i, v in enumerate(res["sample_eq"][:5]):
                            print(f"  [{i}] {v.real:+.3f} {v.imag:+.3f}j")
                    buf = buf[res["consume"]:]
                else:
                    advance = res["consume"] if res["consume"] > 0 else 1
                    buf = buf[advance:]
                    if advance == 0:
                        break
            ctr += 1

        st.issue_stream_cmd(uhd.types.StreamCMD(uhd.types.StreamMode.stop_cont))

    threading.Thread(target=rx_thread, daemon=True).start()

    if verbose:
        print(f"\n[RUN] Target {n_frames} frames (warmup={warmup_frames}). Ctrl+C to stop.\n")
        if phase1b:
            print(f"  {'Frm':>4} | {'text':>4} {'ctr':>3} {'CRC':>4} | "
                  f"{'commEVM':>7} {'snsBER':>7} | {'JCAS':>9} {'Score':>5} | "
                  f"{'echoes (m@dB)':<24} | {'OVF':>4}")
        else:
            print(f"  {'Frm':>4} | {'BER':>8} | {'CFO(Hz)':>8} | {'Amp':>6} | "
                  f"{'EVM(dB)':>7} | {'OVF':>4} | Verdict")

    try:
        last_n = 0
        last_plot_n = 0
        plot_update_interval = 0.15  # ~6.7 Hz
        last_plot_time = 0.0
        while state["n"] < n_frames:
            time.sleep(0.01)
            cur = state["n"]
            if verbose:
                for i in range(last_n, cur):
                    if phase1b and i < len(state["p1b_records"]):
                        rec = state["p1b_records"][i]
                        crc_str = "OK " if rec["crc_ok"] else "FAIL"
                        echoes_str = ", ".join(
                            f"{r:.1f}@{db:.0f}" for r, db in rec["echoes"][:3]
                        ) if rec["echoes"] else "-"
                        marker = " [w]" if (i+1) <= warmup_frames else ""
                        print(f"  {i+1:>4} | {rec['text']:>4} {rec['counter']:>3} {crc_str:>4} | "
                              f"{rec['comm_evm_db']:>7.2f} {rec['sense_ber']:>7.4f} | "
                              f"{rec.get('jcas_status','OFF'):>9} {rec.get('jcas_score',0.0):>5.2f} | "
                              f"{echoes_str:<24} | {state['ovf']:>4}{marker}")
                    else:
                        b = state["bers"][i]; cfo = state["cfos"][i]
                        amp = state["amps"][i]; evm = state["evms"][i]
                        v = "✓ OK" if b < 0.05 else ("⚠ MARG" if b < 0.20 else "✗ BAD")
                        marker = " [warmup]" if (i+1) <= warmup_frames else ""
                        print(f"  {i+1:>4} | {b:>8.4f} | {cfo:>8.0f} | "
                              f"{amp:>6.4f} | {evm:>7.2f} | {state['ovf']:>4} | {v}{marker}")
            last_n = cur

            # Live plot update (rate-limited)
            if plot_handles is not None and cur > last_plot_n:
                now = time.time()
                if now - last_plot_time > plot_update_interval:
                    try:
                        last_plot_n = update_jcas_comm_plot(plot_handles, state, last_plot_n)
                    except Exception as e:
                        print(f"[WARN] Plot update error: {e}")
                        plot_handles = None
                    last_plot_time = now
    except KeyboardInterrupt:
        pass
    finally:
        state["running"] = False
        stop_ev.set()
        tx_proc.join(timeout=5)
        if csv_file is not None:
            csv_file.close()
        if plot_handles is not None:
            try:
                # Final plot update + leave window open until user closes
                update_jcas_comm_plot(plot_handles, state, 0)
                import matplotlib.pyplot as plt
                plt.ioff()
                print("\n[INFO] Plot window tetap terbuka. Tutup window untuk exit.")
                plt.show(block=True)
            except Exception:
                pass

    # ── Aggregate metrics ─────────────────────────────────────────
    if not state["bers"]:
        return {
            "fs_mhz": fs/1e6, "tx_gain": tx_gain, "rx_gain": rx_gain,
            "status": "NO_FRAMES", "n_decoded": 0, "n_target": n_frames,
            "valid_rate": 0.0, "mean_ber": float('nan'), "mean_evm_db": float('nan'),
            "ovf_total": state["ovf"], "ovf_steady": 0, "ovf_rate_steady": 0.0,
            "probe_mean_amp": probe_stats["mean_amp"],
            "probe_peak": probe_stats["peak"],
        }

    n_decoded = state["n"]
    bers_arr = np.array(state["bers"])
    evms_arr = np.array(state["evms"])
    valid_mask = bers_arr < 0.05
    valid_rate = float(np.mean(valid_mask))

    ovf_total = state["ovf"]
    ovf_at_warmup = state["ovf_at_warmup_end"] if state["ovf_at_warmup_end"] is not None else ovf_total
    ovf_steady = max(0, ovf_total - ovf_at_warmup)
    n_steady = max(1, n_decoded - warmup_frames)
    ovf_rate_steady = ovf_steady / n_steady

    metrics = {
        "fs_mhz": fs / 1e6,
        "tx_gain": tx_gain, "rx_gain": rx_gain,
        "status": "OK",
        "n_decoded": n_decoded, "n_target": n_frames,
        "valid_rate": valid_rate,
        "mean_ber": float(np.mean(bers_arr)),
        "median_ber": float(np.median(bers_arr)),
        "mean_evm_db": float(np.mean(evms_arr)),
        "mean_amp": float(np.mean(state["amps"])),
        "mean_cfo_hz": float(np.mean(np.abs(state["cfos"]))),
        "ovf_total": ovf_total,
        "ovf_warmup": ovf_at_warmup,
        "ovf_steady": ovf_steady,
        "ovf_rate_steady": ovf_rate_steady,
        "probe_mean_amp": probe_stats["mean_amp"],
        "probe_peak": probe_stats["peak"],
        "probe_spur_dr": probe_stats["spur_dr"],
    }

    # ── Phase 1B summary ──────────────────────────────────────────
    if phase1b and state["p1b_records"]:
        recs = state["p1b_records"]
        n_p1b = len(recs)
        crc_ok = state["p1b_crc_ok_count"]
        crc_fail = state["p1b_crc_fail_count"]
        prr = crc_ok / max(n_p1b, 1)
        comm_evms = np.array([r["comm_evm_db"] for r in recs])
        sense_evms = np.array([r["sense_evm_db"] for r in recs])
        sense_bers = np.array([r["sense_ber"] for r in recs])
        n_echoes_arr = np.array([r["n_echoes"] for r in recs])
        # Histogram top echo range
        all_echoes = [e[0] for r in recs for e in r["echoes"]]
        metrics.update({
            "phase1b": True,
            "p1b_n_packets": n_p1b,
            "p1b_crc_ok": crc_ok,
            "p1b_crc_fail": crc_fail,
            "p1b_prr": prr,
            "p1b_mean_comm_evm_db": float(np.mean(comm_evms)),
            "p1b_mean_sense_evm_db": float(np.mean(sense_evms)),
            "p1b_mean_sense_ber": float(np.mean(sense_bers)),
            "p1b_mean_n_echoes": float(np.mean(n_echoes_arr)),
            "p1b_total_echoes": len(all_echoes),
        })

        # Recovered text histogram (which texts decoded)
        text_counts = {}
        for r in recs:
            if r["crc_ok"]:
                text_counts[r["text"]] = text_counts.get(r["text"], 0) + 1
        metrics["p1b_text_counts"] = text_counts

        # JCAS detection summary
        jcas_recs = [r for r in recs if "jcas_score" in r]
        if jcas_recs:
            det_count = sum(1 for r in jcas_recs if r.get("jcas_object_detected", False))
            metrics.update({
                "jcas_enabled": True,
                "jcas_detection_count": det_count,
                "jcas_detection_rate": det_count / max(len(jcas_recs), 1),
                "jcas_max_score": float(np.max([r.get("jcas_score", 0.0) for r in jcas_recs])),
                "jcas_mean_amp_hp_db": float(np.mean([abs(r.get("jcas_amp_hp_db", 0.0)) for r in jcas_recs])),
                "jcas_mean_abs_doppler_hp_hz": float(np.mean([abs(r.get("jcas_doppler_hp_hz", 0.0)) for r in jcas_recs])),
            })

    if verbose:
        print(f"\n{'═' * 60}")
        if phase1b and state["p1b_records"]:
            # Phase 1B: tampilkan metrics yang relevan (PRR, comm/sense EVM)
            print(f"  FS = {fs/1e6:.2f} MS/s  | Decoded: {n_decoded}/{n_frames}  "
                  f"| PRR: {metrics['p1b_prr']*100:.2f}%")
            print(f"  ─ Phase 1B (V2I '{text}') ────────────────────────")
            print(f"  Packets        : {metrics['p1b_n_packets']}")
            print(f"  PRR            : {metrics['p1b_prr']*100:.2f}%  "
                  f"(OK={metrics['p1b_crc_ok']}, FAIL={metrics['p1b_crc_fail']})")
            print(f"  Mean comm EVM  : {metrics['p1b_mean_comm_evm_db']:.2f} dB (QPSK)")
            print(f"  Mean sense EVM : {metrics['p1b_mean_sense_evm_db']:.2f} dB (BPSK)")
            print(f"  Mean sense BER : {metrics['p1b_mean_sense_ber']:.4f}")
            print(f"  ─ Phase 1C (range) ──────────────────────────────")
            print(f"  Mean #echoes/frame   : {metrics['p1b_mean_n_echoes']:.2f}")
            print(f"  Total echo detections: {metrics['p1b_total_echoes']}")
            if metrics.get("jcas_enabled"):
                print(f"  ─ Phase 1F (forward scattering JCAS) ───────────")
                print(f"  Object detections : {metrics['jcas_detection_count']} "
                      f"({metrics['jcas_detection_rate']*100:.1f}% frames)")
                print(f"  Max JCAS score    : {metrics['jcas_max_score']:.2f}")
                print(f"  Mean |amp HP|     : {metrics['jcas_mean_amp_hp_db']:.3f} dB")
                print(f"  Mean |Doppler HP| : {metrics['jcas_mean_abs_doppler_hp_hz']:.3f} Hz")
            if metrics["p1b_text_counts"]:
                tc = metrics["p1b_text_counts"]
                top = sorted(tc.items(), key=lambda x: -x[1])[:5]
                print(f"  Recovered texts      : {top}")
            print(f"  ─ Hardware ──────────────────────────────────────")
            print(f"  OVF total      : {ovf_total}  (warmup: {ovf_at_warmup}, "
                  f"steady: {ovf_steady} → {ovf_rate_steady*100:.1f}%/frame)")
            print(f"  NOTE: legacy BER/EVM (Phase 0 path) tidak applicable di Phase 1B mode.")
        else:
            # Phase 0/1A: legacy reporting
            valid_mask_legacy = bers_arr < 0.05
            print(f"  FS = {fs/1e6:.2f} MS/s  | Decoded: {n_decoded}/{n_frames}  "
                  f"| Valid: {float(np.mean(valid_mask_legacy))*100:.1f}%")
            print(f"  Mean BER       : {metrics['mean_ber']:.4f}")
            print(f"  Mean EVM       : {metrics['mean_evm_db']:.2f} dB")
            print(f"  OVF total      : {ovf_total}  (warmup: {ovf_at_warmup}, "
                  f"steady: {ovf_steady} → {ovf_rate_steady*100:.1f}%/frame)")
        print(f"{'═' * 60}\n")

    return metrics


# ═══════════════════════════════════════════════════════════════════
# FS SWEEP — Phase 1A main feature
# ═══════════════════════════════════════════════════════════════════
def evaluate_fs(metrics):
    """Pass criteria sesuai CONTEXT_TRANSFER_v2:
       - OVF rate steady < 5%
       - Valid frame rate >= 90%
       - Mean BER < 0.05
    """
    if metrics is None or metrics.get("status") != "OK":
        return False, "no_data"
    reasons = []
    if metrics["ovf_rate_steady"] >= 0.05:
        reasons.append(f"OVF_high({metrics['ovf_rate_steady']*100:.1f}%)")
    if metrics["valid_rate"] < 0.90:
        reasons.append(f"valid_low({metrics['valid_rate']*100:.1f}%)")
    if metrics["mean_ber"] >= 0.05:
        reasons.append(f"BER_high({metrics['mean_ber']:.3f})")
    if reasons:
        return False, ",".join(reasons)
    return True, "pass"


def run_fs_sweep(fs_candidates, frames_per_fs, tx_gain, rx_gain,
                 probe_wait, output_csv, warmup_frames):
    """Loop over FS candidates, write CSV, identify fs_winner."""
    print(f"\n{'#' * 60}")
    print(f"  FS SWEEP — Phase 1A")
    print(f"  Candidates : {[f/1e6 for f in fs_candidates]} MS/s")
    print(f"  Frames/FS  : {frames_per_fs}  (warmup: {warmup_frames})")
    print(f"  Gain       : TX={tx_gain} dB  RX={rx_gain} dB")
    print(f"  Output CSV : {output_csv}")
    print(f"{'#' * 60}\n")

    results = []
    for i, fs in enumerate(fs_candidates):
        print(f"\n>>> [{i+1}/{len(fs_candidates)}] Testing FS = {fs/1e6:.2f} MS/s ...")
        metrics = run_hardware(
            fs=fs, n_frames=frames_per_fs,
            tx_gain=tx_gain, rx_gain=rx_gain,
            probe_wait=probe_wait, warmup_frames=warmup_frames,
            verbose=False,  # less spam during sweep
        )
        if metrics is None:
            metrics = {"fs_mhz": fs/1e6, "tx_gain": tx_gain, "rx_gain": rx_gain,
                       "status": "INIT_FAIL", "n_decoded": 0, "n_target": frames_per_fs,
                       "valid_rate": 0.0, "mean_ber": float('nan'),
                       "mean_evm_db": float('nan'), "ovf_total": 0,
                       "ovf_steady": 0, "ovf_rate_steady": 0.0,
                       "probe_mean_amp": 0.0, "probe_peak": 0.0}
        passed, reason = evaluate_fs(metrics)
        metrics["pass"] = passed
        metrics["reason"] = reason
        results.append(metrics)

        # Quick summary
        if metrics["status"] == "OK":
            print(f"    Decoded: {metrics['n_decoded']}/{metrics['n_target']} "
                  f"| Valid: {metrics['valid_rate']*100:.1f}% "
                  f"| BER: {metrics['mean_ber']:.4f} "
                  f"| EVM: {metrics['mean_evm_db']:.2f} dB "
                  f"| OVF steady: {metrics['ovf_rate_steady']*100:.1f}%/frame "
                  f"| {'✓ PASS' if passed else '✗ FAIL ('+reason+')'}")
        else:
            print(f"    Status: {metrics['status']}")

        time.sleep(2)  # Cooldown between FS

    # ── Write CSV ─────────────────────────────────────────────────
    fieldnames = ["fs_mhz", "tx_gain", "rx_gain", "status", "pass", "reason",
                  "n_decoded", "n_target", "valid_rate", "mean_ber", "median_ber",
                  "mean_evm_db", "mean_amp", "mean_cfo_hz",
                  "ovf_total", "ovf_warmup", "ovf_steady", "ovf_rate_steady",
                  "probe_mean_amp", "probe_peak", "probe_spur_dr"]
    with open(output_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        for r in results:
            w.writerow({k: r.get(k, "") for k in fieldnames})

    # ── Final summary ────────────────────────────────────────────
    print(f"\n{'═' * 60}")
    print(f"  FS SWEEP SUMMARY")
    print(f"{'═' * 60}")
    print(f"  {'FS(MHz)':>8} | {'Valid%':>6} | {'BER':>8} | {'EVM(dB)':>7} | "
          f"{'OVF%':>6} | Verdict")
    print(f"  {'-'*60}")
    for r in results:
        if r["status"] == "OK":
            verdict = "✓ PASS" if r["pass"] else f"✗ {r['reason']}"
            print(f"  {r['fs_mhz']:>8.2f} | "
                  f"{r['valid_rate']*100:>6.1f} | "
                  f"{r['mean_ber']:>8.4f} | "
                  f"{r['mean_evm_db']:>7.2f} | "
                  f"{r['ovf_rate_steady']*100:>6.1f} | {verdict}")
        else:
            print(f"  {r['fs_mhz']:>8.2f} | {r['status']:^45}")

    # Identify fs_winner: PASS dengan FS tertinggi
    passed = [r for r in results if r.get("pass")]
    if passed:
        winner = max(passed, key=lambda r: r["fs_mhz"])
        print(f"\n  🏆 FS WINNER: {winner['fs_mhz']:.2f} MS/s "
              f"(EVM {winner['mean_evm_db']:.2f} dB, BER {winner['mean_ber']:.4f})")
        print(f"     → Gunakan ini untuk Phase 1B-D.")
    else:
        print(f"\n  ⚠ NO FS PASSED. Check gain/antenna, atau turunkan kriteria.")

    print(f"\n  CSV: {output_csv}")
    print(f"{'═' * 60}\n")

    return results


# ═══════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════
if __name__ == "__main__":
    ap = argparse.ArgumentParser(
        description="ofdm_isac_bistatic — Phase 1A FS sweep + single FS run")

    ap.add_argument("--simulate", action="store_true",
                    help="AWGN self-test (no hardware)")
    ap.add_argument("--probe-only", action="store_true",
                    help="HW: capture 2s, print stats, exit")
    ap.add_argument("--fs-sweep", action="store_true",
                    help="Phase 1A: sweep FS candidates")

    ap.add_argument("--fs", type=float, default=20e6,
                    help="Sample rate Hz (single run mode, default 20e6)")
    ap.add_argument("--fs-candidates", type=str,
                    default="5e6,10e6,15e6,20e6,30e6",
                    help="Comma-separated FS list for sweep")
    ap.add_argument("--frames", type=int, default=100,
                    help="Target frames (single run)")
    ap.add_argument("--frames-per-fs", type=int, default=200,
                    help="Frames per FS in sweep mode")
    ap.add_argument("--warmup-frames", type=int, default=50,
                    help="Frames excluded from steady-state OVF rate")

    ap.add_argument("--tx-gain", type=float, default=80.0)
    ap.add_argument("--rx-gain", type=float, default=70.0)
    ap.add_argument("--probe-wait", type=float, default=10.0)
    ap.add_argument("--frame-delay", type=float, default=0.0)
    ap.add_argument("--no-dc-fix", action="store_true")
    ap.add_argument("--output-csv", type=str, default="fs_sweep_results.csv")

    # ── Phase 1B/1C/1D flags ──────────────────────────────────────
    ap.add_argument("--phase1b", action="store_true",
                    help="Enable Phase 1B (V2I 'STEI' comm) + Phase 1C (range est)")
    ap.add_argument("--text", type=str, default="STEI",
                    help="V2I packet text (4 chars max, padded/truncated to 4)")
    ap.add_argument("--log-csv", type=str, default=None,
                    help="Per-frame Phase 1B log CSV path "
                         "(e.g. isac_log_$(date +%%s).csv)")
    ap.add_argument("--plot", action="store_true",
                    help="Matplotlib live plot. Tanpa --phase1b: komunikasi awal + radar JCAS. "
                         "Dengan --phase1b: plot ISAC lengkap range/const/EVM/log.")
    ap.add_argument("--echo-threshold-db", type=float, default=6.0,
                    help="Echo detection SNR threshold above noise floor (default 6 dB). "
                         "Lower → more sensitive (more false alarms).")

    # ── Phase 1F / JCAS forward-scattering detector flags ─────────
    ap.add_argument("--no-jcas", action="store_true",
                    help="Disable Phase 1F forward-scattering detector")
    ap.add_argument("--sense-ma-len", type=int, default=30,
                    help="Moving median window for LoS amplitude/phase clutter removal")
    ap.add_argument("--sense-cfar-len", type=int, default=160,
                    help="Rolling reference window for simple CFAR threshold")
    ap.add_argument("--sense-threshold-k", type=float, default=5.0,
                    help="CFAR K multiplier; lower is more sensitive")
    ap.add_argument("--sense-min-score", type=float, default=2.2,
                    help="Minimum JCAS detection threshold floor")

    args = ap.parse_args()

    if args.simulate:
        init_params(args.fs, phase1b=args.phase1b, text=args.text)
        run_simulation()
        sys.exit(0)

    mp.set_start_method('spawn')

    if args.fs_sweep:
        if args.phase1b:
            print("⚠ --phase1b not supported in --fs-sweep mode. "
                  "Run single FS with --phase1b instead.")
            sys.exit(1)
        fs_list = [float(x) for x in args.fs_candidates.split(",")]
        run_fs_sweep(
            fs_candidates=fs_list,
            frames_per_fs=args.frames_per_fs,
            tx_gain=args.tx_gain, rx_gain=args.rx_gain,
            probe_wait=args.probe_wait,
            output_csv=args.output_csv,
            warmup_frames=args.warmup_frames,
        )
    else:
        # Default log_csv name kalau phase1b aktif tapi log-csv tidak diisi
        log_csv = args.log_csv
        if args.phase1b and log_csv is None:
            log_csv = f"isac_log_{int(time.time())}.csv"
            print(f"[INFO] Phase 1B log CSV: {log_csv}")
        run_hardware(
            fs=args.fs,
            n_frames=args.frames,
            tx_gain=args.tx_gain, rx_gain=args.rx_gain,
            probe_wait=args.probe_wait,
            frame_delay=args.frame_delay,
            dc_offset_auto=not args.no_dc_fix,
            probe_only=args.probe_only,
            warmup_frames=args.warmup_frames,
            phase1b=args.phase1b,
            text=args.text,
            log_csv=log_csv,
            plot=args.plot,
            echo_threshold_db=args.echo_threshold_db,
            jcas=not args.no_jcas,
            sense_ma_len=args.sense_ma_len,
            sense_cfar_len=args.sense_cfar_len,
            sense_threshold_k=args.sense_threshold_k,
            sense_min_score=args.sense_min_score,
        )


