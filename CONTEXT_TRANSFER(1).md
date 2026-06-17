# CONTEXT_TRANSFER.md

## 1. Tujuan Dokumen

Dokumen ini dibuat sebagai **context transfer** untuk membuka chat baru tanpa kehilangan arah proyek. Isi dokumen merangkum posisi terakhir pekerjaan, lesson learned, bukti/proof yang perlu dilampirkan, masalah yang masih ada, serta workplan lanjutan untuk pengembangan sistem **JCAS / ISAC berbasis USRP**.

Repo GitHub minimal yang direncanakan berisi:

```text
.
├── kode_terbaru/ atau src/
│   └── kode utama MATLAB/Python terbaru
├── CONTEXT_TRANSFER.md
└── README.md
```

Catatan penting: file ZIP terakhir yang diberikan berisi `baru4.pdf`. Dari isi PDF, kode yang terbaca adalah **Python script `ofdm_isac_bistatic.py`**, bukan file `.m` MATLAB langsung. Jika target akhir GitHub wajib menyertakan **kode MATLAB**, maka kode MATLAB final tetap perlu dilampirkan terpisah di chat/repo baru.

---

## 2. Ringkasan Proyek

Proyek ini bertujuan mengembangkan sistem **Joint Communication and Sensing (JCAS)** atau **Integrated Sensing and Communication (ISAC)** menggunakan dua USRP dalam konfigurasi bistatic / forward scattering.

Konsep sistem:

- **USRP A / Lutetia** berperan sebagai transmitter (TX).
- **USRP B / Libre** berperan sebagai receiver (RX).
- TX mengirimkan sinyal OFDM.
- RX menerima sinyal komunikasi sekaligus memantau perubahan sinyal akibat objek yang melintas di antara TX dan RX.
- Deteksi objek dilakukan dari perubahan jalur Line-of-Sight (LoS), amplitudo, phase, Doppler proxy, EVM, dan estimasi CIR/range profile.

Target demo yang diinginkan:

- Ada komunikasi teks, misalnya `STEI`, dari TX ke RX.
- Ada visualisasi real-time yang mudah dipahami seperti radar monitor.
- Tampilan harus jelas menunjukkan apakah ada objek atau tidak.
- Sistem harus bisa digunakan untuk pengujian objek yang bergerak di antara dua USRP.

---

## 3. Posisi Terakhir Kode

Berdasarkan PDF terakhir, kode utama bernama:

```text
ofdm_isac_bistatic.py
```

Status fitur pada kode:

| Phase | Status | Keterangan |
|---|---|---|
| Phase 1A | Implemented | FS sweep mode untuk mencari sample rate maksimum yang feasible di hardware |
| Phase 1B | Implemented | Payload komunikasi V2I berupa teks `STEI` menggunakan QPSK di subcarrier tengah |
| Phase 1C | Implemented | Range estimation berbasis CIR + delta-range |
| Phase 1D | Implemented | Real-time matplotlib live plot dengan flag `--plot` |
| Phase 1E | Pending | Validasi laboratorium untuk skenario moving target |
| Phase 1F | Implemented sebagian | Forward-scattering / LoS disruption detector berbasis amplitude, phase, Doppler proxy, dan adaptive threshold |

Konfigurasi penting yang terlihat dari kode:

```text
TX_SERIAL = "000000037"
RX_SERIAL = "HQHGTFH"
TX role   = Lutetia / TX
RX role   = Libre / RX
FC        = 5.9 GHz
NFFT      = 64
NCP       = 16
Frame     = STF | LTF1 | LTF2 | DATA x 26
Payload   = "STEI" + counter + CRC16
```

Catatan sebelum publish GitHub: serial number USRP bisa dianggap informasi perangkat. Jika repository publik, sebaiknya serial diganti placeholder seperti `TX_SERIAL = "YOUR_TX_SERIAL"` dan `RX_SERIAL = "YOUR_RX_SERIAL"`.

---

## 4. Posisi Setup Hardware yang Benar

Pertanyaan penting: **apakah posisi setup sekarang masih salah?**

Jawaban ringkas: **kemungkinan iya, jika objek tidak benar-benar melintas di jalur Line-of-Sight antara TX dan RX.**

Untuk konsep forward scattering / bistatic JCAS, posisi yang benar bukan seperti radar monostatic yang memantulkan sinyal ke arah sensor yang sama. Sistem ini mengandalkan perubahan sinyal langsung antara dua USRP.

Setup yang seharusnya:

```text
[USRP TX]  ---> jalur LoS / Fresnel zone --->  [USRP RX]
                  objek lewat di sini
```

Aturan posisi:

1. TX dan RX harus saling berhadapan.
2. Antena TX dan RX harus berada pada tinggi yang relatif sejajar.
3. Objek harus lewat di tengah atau dekat jalur garis lurus TX-RX.
4. Jangan menaruh objek terlalu jauh di samping jalur TX-RX, karena perubahan sinyal bisa sangat kecil.
5. Untuk demo awal, jarak TX-RX jangan terlalu jauh. Mulai dari jarak pendek dan stabil, lalu naik bertahap.
6. Area sekitar sebaiknya tidak terlalu banyak pergerakan orang lain, karena semua gerakan di sekitar LoS bisa memengaruhi sinyal.

Hal yang perlu dipahami:

- Sistem ini lebih cocok disebut **LoS obstruction / forward-scattering object detection**, bukan radar posisi absolut.
- Plot radar yang dibuat adalah visualisasi agar mudah dilihat manusia, tetapi secara fisik sistem belum tentu benar-benar mengetahui koordinat 2D objek seperti radar pesawat.
- Estimasi range dari CIR adalah **delta-range relatif**, bukan posisi objek absolut yang presisi.
- Dengan dua USRP yang tidak tersinkronisasi clock/time secara ketat, absolute ranging akan sulit dan perlu validasi tambahan.

---

## 5. Lesson Learned

### 5.1 Dari sisi konsep JCAS/ISAC

- JCAS bukan hanya mengirim data, tetapi memanfaatkan sinyal komunikasi yang sama untuk sensing.
- Untuk dua USRP, pendekatan yang paling realistis saat ini adalah **bistatic sensing** atau **forward scattering**, bukan radar monostatic penuh.
- Deteksi objek lebih mudah dilihat dari perubahan amplitudo, phase, EVM, dan pola CIR dibanding langsung mencari koordinat objek.
- Jika hanya menggunakan dua node TX-RX, hasil sensing harus dijelaskan sebagai deteksi perubahan kanal, bukan tracking posisi absolut.

### 5.2 Dari sisi komunikasi OFDM

- OFDM cocok karena subcarrier bisa dibagi antara komunikasi dan sensing.
- Kode saat ini membagi subcarrier menjadi:
  - subcarrier komunikasi QPSK untuk payload teks,
  - subcarrier sensing BPSK untuk estimasi kanal,
  - pilot untuk referensi CPE.
- Payload `STEI` sudah menjadi proof bahwa komunikasi TX-RX dapat dimasukkan ke frame OFDM.
- CRC16 penting untuk membedakan payload yang benar dan payload yang error.

### 5.3 Dari sisi sensing

- CIR/range profile dapat digunakan untuk melihat perubahan kanal.
- Echo threshold perlu dituning; threshold terlalu tinggi membuat objek tidak terdeteksi, threshold terlalu rendah membuat false alarm.
- Background subtraction penting agar perubahan kecil akibat objek bisa terlihat lebih jelas.
- Detector berbasis moving median dan CFAR sederhana lebih stabil dibanding hanya melihat amplitudo mentah.

### 5.4 Dari sisi visualisasi

- Visualisasi awal yang terlalu teknis sulit dipahami orang awam.
- Untuk demo, perlu tampilan sederhana:
  - status `CLEAR` atau `OBJECT DETECTED`,
  - radar-style monitor,
  - teks komunikasi TX -> RX,
  - grafik score vs threshold.
- Visualisasi radar harus diberi disclaimer: itu adalah **monitor deteksi**, bukan radar posisi absolut.

### 5.5 Dari sisi hardware USRP

- Dua USRP dapat membutuhkan image FPGA yang berbeda.
- Lutetia dan Libre tidak boleh tertukar role-nya.
- Sample rate besar seperti 40 MS/s dapat feasible, tetapi harus dibuktikan dengan FS sweep, valid rate, BER, EVM, dan overflow rate.
- Gain TX/RX perlu tuning. Sinyal terlalu kecil membuat decoding gagal, sinyal terlalu besar dapat menyebabkan distorsi/saturasi.
- Perlu memastikan RX tidak overflow, karena overflow membuat frame decoding dan sensing tidak stabil.

---

## 6. Proof / Bukti yang Perlu Dilampirkan di Chat Baru atau GitHub

Agar proyek terlihat kuat dan mudah diverifikasi, lampirkan bukti berikut:

### 6.1 Bukti kode

- Kode utama terbaru:
  - jika versi final Python: `ofdm_isac_bistatic.py`
  - jika versi final MATLAB: file `.m` terbaru
- Screenshot struktur folder repo.
- Screenshot bagian kode yang menunjukkan:
  - konfigurasi TX/RX,
  - frame OFDM,
  - payload `STEI`,
  - detector JCAS,
  - live plot.

### 6.2 Bukti komunikasi

- Screenshot terminal saat RX berhasil decode teks `STEI`.
- Log CSV yang menunjukkan:
  - `text`,
  - `counter`,
  - `crc_ok`,
  - `comm_evm_db`,
  - `sense_evm_db`,
  - `sense_ber`,
  - `cfo_hz`.

### 6.3 Bukti sensing

- Screenshot live plot saat kondisi tidak ada objek.
- Screenshot live plot saat objek melintas di antara TX dan RX.
- CSV log sebelum, saat, dan sesudah objek lewat.
- Bukti bahwa `jcas_score` naik melewati `jcas_threshold` saat objek lewat.
- Foto setup fisik TX-RX dan posisi objek.

### 6.4 Bukti FS sweep

- File `fs_sweep_results.csv`.
- Screenshot terminal FS sweep summary.
- Bukti sample rate winner, misalnya 40 MS/s jika valid.
- Catatan parameter:
  - TX gain,
  - RX gain,
  - jarak TX-RX,
  - jenis antena,
  - sample rate,
  - frekuensi carrier,
  - lokasi eksperimen.

### 6.5 Bukti masalah yang masih ada

- Screenshot saat false alarm terjadi.
- Screenshot saat objek tidak terdeteksi walaupun ada objek.
- Catatan posisi objek saat gagal deteksi.
- Catatan apakah antena sudah sejajar atau belum.

---

## 7. Masalah yang Masih Perlu Diperbaiki

1. **Posisi fisik setup kemungkinan belum ideal.** Objek harus lewat di jalur LoS/Fresnel zone, bukan di samping USRP.
2. **Visualisasi radar bisa menyesatkan** jika tidak dijelaskan bahwa radar-style plot hanya representasi deteksi.
3. **Range estimation belum boleh diklaim sebagai posisi absolut.** Saat ini lebih aman disebut delta-range / channel change indicator.
4. **Validasi moving target belum selesai.** Phase 1E masih pending.
5. **Threshold detector perlu dikalibrasi.** Parameter `sense_threshold_k`, `sense_min_score`, dan `echo_threshold_db` perlu dicoba untuk kondisi ruangan berbeda.
6. **Repo GitHub belum clean.** Perlu memastikan file yang diunggah tidak hanya PDF kode, tetapi source code asli.
7. **Kode MATLAB belum terlihat dari upload terakhir.** Jika GitHub harus berisi kode MATLAB, file `.m` final harus dilampirkan lagi.

---

## 8. Workplan Lanjutan

### Step 1 - Rapikan repository

Buat struktur repo yang jelas:

```text
JCAS-USRP/
├── src/
│   ├── ofdm_isac_bistatic.py atau main_jcas.m
│   └── helper files jika ada
├── docs/
│   ├── setup_photo.jpg
│   ├── no_object_plot.png
│   └── object_detected_plot.png
├── logs/
│   ├── fs_sweep_results.csv
│   └── isac_log_test.csv
├── CONTEXT_TRANSFER.md
└── README.md
```

### Step 2 - Pastikan file kode asli tersedia

Jangan hanya upload PDF berisi kode. Upload source code asli:

- Python: `.py`
- MATLAB: `.m`

Jika masih menggunakan Python, nama yang disarankan:

```text
src/ofdm_isac_bistatic.py
```

Jika harus MATLAB, nama yang disarankan:

```text
src/main_jcas_ofdm.m
```

### Step 3 - Validasi setup posisi

Lakukan pengujian dengan tiga kondisi:

1. Tidak ada objek di antara TX-RX.
2. Objek diam di tengah jalur TX-RX.
3. Objek bergerak melintas dari kiri ke kanan melewati jalur TX-RX.

Catat hasilnya di CSV dan screenshot live plot.

### Step 4 - Kalibrasi threshold

Coba beberapa parameter:

```bash
--echo-threshold-db 4
--echo-threshold-db 6
--sense-threshold-k 3.5
--sense-threshold-k 5.0
--sense-min-score 1.8
--sense-min-score 2.2
```

Pilih parameter yang paling stabil untuk demo.

### Step 5 - Pisahkan klaim teknis

Di laporan/GitHub, gunakan wording yang aman:

- Benar: `object presence detection using LoS disruption`
- Benar: `forward-scattering based channel-change detection`
- Benar: `relative CIR/range profile monitoring`
- Hindari: `accurate 2D radar tracking`
- Hindari: `absolute target position estimation`

### Step 6 - Siapkan demo final

Demo final yang disarankan:

1. Jalankan kode.
2. Tampilkan komunikasi `STEI` berhasil diterima.
3. Tunjukkan status `CLEAR` saat tidak ada objek.
4. Lewatkan objek/orang di antara USRP.
5. Tunjukkan status berubah menjadi `OBJECT DETECTED`.
6. Simpan screenshot dan CSV log.

---

## 9. Command Penting dari Kode Saat Ini

### Simulasi tanpa hardware

```bash
python3 ofdm_isac_bistatic.py --simulate --fs 20e6
python3 ofdm_isac_bistatic.py --simulate --fs 40e6 --phase1b
```

### FS sweep

```bash
python3 ofdm_isac_bistatic.py --fs-sweep \
  --fs-candidates 5e6,10e6,20e6,30e6,40e6 \
  --frames-per-fs 200 \
  --tx-gain 80 \
  --rx-gain 70
```

### Single run komunikasi + sensing

```bash
python3 ofdm_isac_bistatic.py --phase1b \
  --fs 40e6 \
  --frames 500 \
  --tx-gain 80 \
  --rx-gain 70 \
  --text "STEI"
```

### Run dengan live plot

```bash
python3 ofdm_isac_bistatic.py --phase1b --plot \
  --fs 40e6 \
  --frames 500 \
  --tx-gain 80 \
  --rx-gain 70
```

### Run dengan log CSV

```bash
python3 ofdm_isac_bistatic.py --phase1b \
  --fs 40e6 \
  --frames 1000 \
  --tx-gain 80 \
  --rx-gain 70 \
  --log-csv isac_log_lab_test.csv
```

---

## 10. Hal yang Perlu Dibawa ke Chat Baru

Saat membuka chat baru, lampirkan minimal:

1. `CONTEXT_TRANSFER.md`
2. `README.md`
3. Source code terbaru asli, bukan PDF:
   - `.py` jika lanjut dari kode sekarang, atau
   - `.m` jika wajib MATLAB.
4. Screenshot live plot terbaru.
5. CSV log terbaru.
6. Foto posisi USRP TX-RX dan posisi objek.
7. Jelaskan role hardware:
   - Lutetia = TX
   - Libre = RX
8. Jelaskan OS dan environment:
   - Linux/Windows yang dipakai,
   - UHD version,
   - Python/MATLAB version,
   - apakah pakai GNU Radio atau tidak.

Prompt pembuka chat baru yang disarankan:

```text
Aku mau lanjut proyek JCAS/ISAC USRP. Aku lampirkan CONTEXT_TRANSFER.md, README.md, kode terbaru, screenshot plot, dan log CSV. Tolong lanjutkan dari posisi terakhir ini. Fokus utama sekarang: validasi posisi TX-RX-object, rapikan kode untuk GitHub, dan buat demo object detection yang jelas.
```

---

## 11. Kesimpulan

Proyek sudah berkembang dari komunikasi USRP biasa menjadi sistem OFDM JCAS/ISAC dengan komunikasi teks, estimasi kanal, range profile, live plot, dan detector berbasis forward scattering. Namun, validasi fisik masih sangat penting. Posisi objek harus benar-benar berada di antara TX dan RX agar konsep forward scattering bekerja. Untuk GitHub, repo perlu dirapikan dengan source code asli, README, context transfer, bukti eksperimen, dan workplan yang jelas.
