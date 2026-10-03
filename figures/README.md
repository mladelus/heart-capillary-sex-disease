# Figures of the combined manuscript

The figures are drawn with Python (matplotlib) from the result files that the R scripts save. They display results; they compute no new statistics.

```
pip install numpy pandas matplotlib pillow
python make_cache.py      # reads the .rds result files of both projects into D.pkl
python figs_main.py       # Figures 2-5 and 7
python figs_extra.py      # Figure 1 and Figures S1-S6
python figs_state.py      # Figure 6 (after the scripts in ../state)
```

`make_cache.py` looks for the result files in `~/Documents/heart_capillary_project` (healthy hearts) and `~/Documents/heart_capillary_disease` (diseased hearts and exploratory analyses); set `HCAP_DATA` and `HCAP2_DATA` to change this. Run R scripts 00-20 first. Output goes to `figures/figures/` as PDF, TIFF (300 dpi) and PNG. Each run checks that no label is cut off or overlaps another panel.
