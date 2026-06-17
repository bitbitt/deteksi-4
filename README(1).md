# JCAS / ISAC USRP OFDM Bistatic System

This repository contains the latest development of a **Joint Communication and Sensing (JCAS)** / **Integrated Sensing and Communication (ISAC)** experiment using two USRP devices in a bistatic forward-scattering configuration.

The system is designed to transmit an OFDM communication payload while simultaneously observing channel changes caused by an object moving between the transmitter and receiver.

> Current note: the latest uploaded file reviewed for this README was `baru4.pdf`, which contains the Python script `ofdm_isac_bistatic.py`. If the final GitHub submission must use MATLAB, place the final `.m` code in the repository and update the usage section accordingly.

---

## 1. Project Overview

The main goal is to build a simple and demonstrable JCAS system where:

- one USRP acts as the transmitter,
- one USRP acts as the receiver,
- the transmitter sends an OFDM signal carrying a small text payload such as `STEI`,
- the receiver decodes the communication payload,
- the receiver also detects channel changes when an object passes through the Line-of-Sight path between TX and RX.

This project uses the idea of **bistatic forward scattering**, meaning the object is detected mainly from the disruption it causes in the direct TX-RX signal path.

---

## 2. Current Features

| Feature | Status | Description |
|---|---|---|
| OFDM frame generation | Implemented | Uses STF, LTF, and data symbols |
| Communication payload | Implemented | Sends text such as `STEI` with counter and CRC16 |
| QPSK communication subcarriers | Implemented | Center subcarriers carry the communication packet |
| Sensing subcarriers | Implemented | BPSK known symbols used for channel estimation |
| CIR/range profile | Implemented | Uses channel estimate and zero-padded IFFT |
| FS sweep | Implemented | Tests multiple sample rates to find the most stable setting |
| Live plot | Implemented | Matplotlib visualization for communication and sensing |
| Forward-scattering detector | Partially implemented | Uses amplitude, phase, Doppler proxy, and adaptive threshold |
| Moving target lab validation | Pending | Needs controlled experimental proof |

---

## 3. Hardware Configuration

The current known hardware role assignment is:

```text
Lutetia = TX
Libre   = RX
```

The setup should follow this geometry:

```text
[USRP TX]  ---- Line-of-Sight / Fresnel Zone ----  [USRP RX]
                         ↑
                  object moves here
```

Important setup rules:

1. TX and RX antennas should face each other.
2. The antennas should be placed at approximately the same height.
3. The object must pass between TX and RX, not far to the side.
4. Start with a short and stable TX-RX distance before increasing the range.
5. Keep the surrounding area as still as possible during baseline measurement.

This system should not be claimed as a full 2D radar tracker yet. It is more accurately described as **object presence detection from Line-of-Sight channel disruption**.

---

## 4. Repository Structure

Recommended structure:

```text
JCAS-USRP/
├── src/
│   ├── ofdm_isac_bistatic.py
│   └── main_jcas_ofdm.m              # optional, if MATLAB version is required
├── docs/
│   ├── setup_photo.jpg
│   ├── no_object_plot.png
│   └── object_detected_plot.png
├── logs/
│   ├── fs_sweep_results.csv
│   └── isac_log_lab_test.csv
├── CONTEXT_TRANSFER.md
└── README.md
```

Minimum files to upload to GitHub:

```text
source code file
CONTEXT_TRANSFER.md
README.md
```

If the final requirement says MATLAB, include the latest `.m` file. If the project continues from the current reviewed code, include the `.py` source file rather than a PDF copy of the code.

---

## 5. Software Requirements

For the current Python version:

- Python 3
- NumPy
- Matplotlib
- UHD Python API
- USRP UHD driver and FPGA images

Example Python packages:

```bash
pip install numpy matplotlib
```

UHD installation depends on the operating system and USRP model.

---

## 6. Running the Code

### 6.1 Simulation without hardware

```bash
python3 ofdm_isac_bistatic.py --simulate --fs 20e6
python3 ofdm_isac_bistatic.py --simulate --fs 40e6 --phase1b
```

### 6.2 Sample-rate sweep

```bash
python3 ofdm_isac_bistatic.py --fs-sweep \
  --fs-candidates 5e6,10e6,20e6,30e6,40e6 \
  --frames-per-fs 200 \
  --tx-gain 80 \
  --rx-gain 70
```

The FS sweep should produce a CSV file such as:

```text
fs_sweep_results.csv
```

This file is used to justify which sample rate is stable enough for the hardware.

### 6.3 Communication and sensing run

```bash
python3 ofdm_isac_bistatic.py --phase1b \
  --fs 40e6 \
  --frames 500 \
  --tx-gain 80 \
  --rx-gain 70 \
  --text "STEI"
```

### 6.4 Run with live visualization

```bash
python3 ofdm_isac_bistatic.py --phase1b --plot \
  --fs 40e6 \
  --frames 500 \
  --tx-gain 80 \
  --rx-gain 70
```

### 6.5 Run with CSV logging

```bash
python3 ofdm_isac_bistatic.py --phase1b \
  --fs 40e6 \
  --frames 1000 \
  --tx-gain 80 \
  --rx-gain 70 \
  --log-csv isac_log_lab_test.csv
```

---

## 7. Expected Outputs

The system should produce:

1. decoded communication text, for example `STEI`,
2. CRC status for each packet,
3. communication EVM,
4. sensing EVM or sensing BER,
5. CFO estimate,
6. CIR/range profile,
7. object detection status,
8. JCAS score and adaptive threshold,
9. CSV logs for post-analysis.

Example fields expected in the CSV log:

```text
frame_idx
text
counter
crc_ok
comm_evm_db
sense_evm_db
sense_ber
cfo_hz
amp
n_echoes
echoes_m
jcas_status
jcas_score
jcas_threshold
jcas_object_detected
```

---

## 8. Lesson Learned

### 8.1 Communication and sensing can share one OFDM frame

The project shows that the same OFDM signal can be divided into communication and sensing resources. Communication subcarriers can carry QPSK data, while sensing subcarriers can carry known BPSK symbols for channel estimation.

### 8.2 Forward scattering depends heavily on object position

The most important experimental lesson is that the object must pass through the TX-RX Line-of-Sight region. If the object is placed outside the main path, the detector may not respond clearly.

### 8.3 The radar display is a visualization, not yet true radar localization

The live plot is useful for demonstration, but it should not be interpreted as accurate aircraft-style radar tracking. At the current stage, it is a user-friendly object detection monitor.

### 8.4 Threshold tuning is critical

Detection depends on parameters such as:

```text
echo_threshold_db
sense_threshold_k
sense_min_score
sense_ma_len
sense_cfar_len
```

A lower threshold increases sensitivity but may create more false alarms. A higher threshold reduces false alarms but may miss weak object movement.

### 8.5 Hardware validation is as important as code implementation

A working script is not enough. The project needs repeatable proof using logs, screenshots, and photos of the physical setup.

---

## 9. Proof to Include in the Repository

Recommended proof files:

```text
docs/setup_photo.jpg
docs/no_object_plot.png
docs/object_detected_plot.png
logs/fs_sweep_results.csv
logs/isac_log_lab_test.csv
```

Recommended proof descriptions:

1. TX-RX physical placement.
2. No-object baseline condition.
3. Object crossing condition.
4. Communication payload successfully decoded.
5. JCAS score rising above threshold during object movement.
6. FS sweep showing stable sample rate selection.

---

## 10. Current Limitations

- Absolute target position estimation is not yet validated.
- CIR-based range should be treated as relative delta-range, not final object coordinates.
- Moving target validation is still pending.
- Detector threshold still requires calibration for different rooms and antenna positions.
- Hardware clock synchronization limitations may affect precise ranging.
- If using a public repository, hardware serial numbers should be replaced with placeholders.

---

## 11. Workplan

### Phase A - Repository cleanup

- Put the actual source code in `src/`.
- Avoid uploading only a PDF version of the code.
- Add `README.md` and `CONTEXT_TRANSFER.md`.
- Add logs and screenshots in separate folders.

### Phase B - Hardware setup validation

- Place TX and RX facing each other.
- Measure no-object baseline.
- Move an object through the middle of the TX-RX path.
- Record screenshots and CSV logs.

### Phase C - Detector calibration

Test several parameter combinations:

```bash
--echo-threshold-db 4
--echo-threshold-db 6
--sense-threshold-k 3.5
--sense-threshold-k 5.0
--sense-min-score 1.8
--sense-min-score 2.2
```

Choose the setting that gives the clearest separation between no-object and object-present conditions.

### Phase D - Documentation and report

Document:

- system architecture,
- OFDM frame design,
- communication payload,
- sensing principle,
- hardware placement,
- experimental results,
- limitations.

### Phase E - Final demo

The final demo should show:

1. successful communication text decoding,
2. clear baseline condition,
3. object crossing event,
4. detection status changing from clear to detected,
5. CSV proof that the JCAS score crosses the threshold.

---

## 12. Recommended Wording for Reports

Use these safe technical claims:

- `The system demonstrates OFDM-based communication and channel-change sensing using two USRP devices.`
- `Object presence is detected through Line-of-Sight disruption in a bistatic forward-scattering setup.`
- `The range profile is used as a relative channel indicator, not as absolute target localization.`
- `The live radar-style visualization is designed for intuitive monitoring of object detection status.`

Avoid overclaiming:

- `accurate 2D target tracking`
- `aircraft-like radar localization`
- `precise absolute range estimation`
- `fully validated moving-target radar`

---

## 13. What Was Learned

From this project, the main learning outcomes are:

1. understanding how OFDM can support both communication and sensing,
2. learning how USRP hardware constraints affect sample rate, gain, and overflow,
3. learning how channel estimation can be used for sensing,
4. learning why physical antenna and object position strongly affect detection,
5. learning how to build a real-time visualization for engineering demonstration,
6. learning the importance of honest technical claims and experimental proof.

---

## 14. Final Status

The project has reached a functional prototype stage. Communication and sensing features are implemented in code, but the system still needs controlled validation with correct TX-RX-object positioning. The next priority is to collect strong proof from experiments and clean the repository for GitHub submission.
