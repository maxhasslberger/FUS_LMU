# FUS_LMU
Repository used for tFUS simulations for clinical studies to stimulate the human motor cortex at the LMU hospital in Munich. Some functions originate from other repositories, like [PRESTUS](https://github.com/achetverikov/PRESTUS). Currently, only the [NeuroFUS PRO CTX-500 4-element transducer](https://brainbox-neuro.com/products/neurofus) is supported. Thermal simulations are also theoretically possible but haven't been tested yet.

## Dependencies
- [k-wave](http://www.k-wave.org)
- [BRIC TUS Simulation Tools](https://github.com/sitiny/BRIC_TUS_Simulation_Tools)
- [kArray tools](http://www.k-wave.org/downloads/kWaveArray_alpha_0.3.zip)

## Running acoustic simulations
0. Generate Pseudo CT scans for your T1w scans (1 mm) based on this [Repository](https://github.com/sitiny/mr-to-pct).
1. Type 'simulationApp' in your Matlab cmd window to launch the GUI
2. `Init`__
   2.1. Load settings from past simulations (if applicable) to skip the following `Init` and `Subject Param` steps__
   2.2. Enter complete paths to T1w scan and Pseudo CT__
   2.3. Enter sonication protocol and subject name
3. `Subject Param`__
   3.1. Enter subject (spatial T1w) offset - Open subject T1w scan and enter the (outermost negative) starting index as a positive integer for each dimension (May be exact or close to (Nx, Ny, Nz)  / 2)__
   3.2. Enter Focus coordinates (matching the T1w scan), minimum offset of gel pad (might be desired for an inactive control condition), additional offset (this offset parameter is not considered when computing the focus depth), and ISPPA in water__
   3.3. The Bowl axis is always given relative to a focus point at (128, 128, 128). It solely determines the alignment of the transducer axis, not the focus depth! It can be tuned arbitrarily based on the visual output in `Prepare simulation`.
4. Enter complete Results path
5. Push `Prepare simulation` to observe the focus point on the T1w scan from different views as well as the position of the transducer face center relative to both focus point and Pseudo CT slice.
6. Push `Run simulation` to carry out acoustic simulation and store the results in the specified path

