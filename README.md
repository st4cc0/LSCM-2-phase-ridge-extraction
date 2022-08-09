# Extraction of dynamic, 2-phase wetting ridges from laser scanning confocal (LSCM) images

This code processes images of wetting ridges, captured by laser scanning confocal microscopy. The inputs are 2-channel .tiff files (8bit grayscale). The codes were utilized to process data in the publication - Hauer, Cai, Vollmer and Pham (2022) Phase Separated Wetting Ridges of Sliding Drops on Soft and Swollen Surfaces

-Input files are placed in ./01_Files to Process
-Output files are generated in ./02_Results

The repo holds 3 codes:

1) ExtractRidges                  ;this is a simple profile extraction of both channels
2) ExtractRidges_TempAv_Profile   ;this aligns the PROFILES to the max. ridge point and averages it over all time frames
3) ExtractRidges_TempAv_Frame     ;this aligns the LSCM FRAMES to the max. ridge point and averages it over all time frames

(CC-2022 by Lukas Hauer)
