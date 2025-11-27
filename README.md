# AUDIO AMPLIFIER — QUALITY CONTROL AND DIGITAL EQUALIZER

Description
-----------
Digital amplification and equalizer system, both developed in MATLAB along with the application that runs them.

This repository contains MATLAB code, manuals, and project documents for measuring, validating, and improving audio amplifier performance and for implementing a configurable digital equalizer. The repository's codebase is written entirely in MATLAB.

Repository layout
-----------------
- Manuals/
  - Operator, hardware and user manuals (setup guides, wiring diagrams, safety notes).
- Matlab Code versions/
  - Versioned MATLAB code and scripts. Open the subfolder for the specific version you want to run.
- Project documents/
  - Reports, design documents, datasheets and project deliverables.

What you will find here
-----------------------
- MATLAB implementations of measurement routines (frequency response, THD, SNR, noise floor, gain/offset checks).
- MATLAB scripts or functions implementing a digital equalizer (IIR/FIR implementations, band definitions, presets).
- Documentation to operate the system (in Manuals/) and supporting project reports (in Project documents/).

Prerequisites
-------------
- MATLAB (recommended recent release such as R2018b or later).
- Typical toolboxes likely required:
  - Signal Processing Toolbox
  - DSP System Toolbox
  - (Optional) Audio Toolbox
- Basic audio test hardware: signal source (sound card or DAC), DUT (amplifier), measurement device (sound card, ADC, or audio analyzer), load resistor (e.g., 8 Ω).

Quick start (MATLAB)
--------------------
1. Clone the repository:
   git clone https://github.com/aitorpitarchfontcuberta/AUDIO-AMPLIFIER-QUALITY-CONTROL-AND-DIGITAL-EQUALIZER.git

2. Open MATLAB and add the repository to your path:
   - In MATLAB:
     addpath(genpath(fullfile(pwd, 'Matlab Code versions')));
     savepath;

3. Choose a MATLAB code version:
   - Open the folder Matlab Code versions/ and select the version subfolder you want to run.
   - Look for an entry script or README inside that version folder. Common entry scripts are named `main.m`, `run_all.m`, or a version-specific script.

4. Run the entry script (example pattern — replace with the actual entry script name found in the version folder):
   cd(fullfile(pwd, 'Matlab Code versions', '<version_folder>'));
   run('main.m');    % or run('run_all.m') or the file indicated in that folder's README

General usage notes
-------------------
- Configuration: many measurement scripts expect a configuration structure or file (e.g., sample rate, FFT size, input level, output attenuation). Inspect the chosen version folder for config files or an introductory script that creates a configuration.
- Data I/O: measurement outputs are typically saved to MAT files and/or exported to CSV/plots. Look for `save`, `writetable`, or plotting scripts in the code.
- EQ presets: equalizer presets (if included) will be stored as .mat/.json/.m files inside the version folder or a dedicated presets subfolder. Load and apply them via the provided filter-design and processing functions.

Suggested measurement workflow
------------------------------
1. Prepare hardware: connect signal source → DUT input, DUT output → measurement input. Ensure proper grounding and safe load resistors.
2. Open the appropriate MATLAB version folder and load the measurement configuration.
3. Run the sweep/measurement script to obtain frequency response, THD+N, SNR.
4. Use the equalizer design scripts to compute corrective filters from measured response.
5. Apply the filters in real-time (if supported by your hardware) or process captured audio files to verify correction.
6. Save the test report, plots and configuration used.

What to check if something fails
-------------------------------
- Missing toolboxes: MATLAB will error if required toolboxes are not installed. Use ver or the Add-Ons manager to check.
- No entry script found: open the version folder and inspect filenames; look for README files inside version folders that point to the main script.
- Hardware issues: ensure wiring, power rails and measurement levels are correct and not clipping.

Contributing
------------
- If you want to improve code, add versions, or add documentation: fork the repository, create a branch, and open a pull request with a clear description of changes.
- Open issues for bug reports or feature requests and include reproduction steps and any relevant MATLAB errors or screenshots.

License and attribution
-----------------------
- Add a LICENSE file to the repository if a license is intended (e.g., MIT, Apache-2.0).
- If you reuse code or external libraries, include attribution in the Project documents/ or a dedicated CREDITS file.

Contact / Maintainer
--------------------
Maintainer: aitorpitarchfontcuberta

If more details are needed (for example: which entry script to run, specific toolbox requirements for a particular version, or concrete instructions for your test hardware), inspect the subfolders inside Matlab Code versions/ and Manuals/ for version-specific instructions and scripts.
