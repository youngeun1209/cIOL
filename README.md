# cIOL
cIOL (constrained ICA with Online Learning) algorithm

### citation
Y.-E. Lee, N.-S. Kwak, S.-W. Lee, "A Real-Time Movement Artifact Removal Method for Ambulatory Brain-Computer Interfaces," IEEE Trans. on Neural Systems & Rehabilitation Engineering, Vol. 28, No.12, 2020, pp. 2660-2670.

The dataset is available in figshare repository (https://doi.org/10.6084/m9.figshare.13604078) under the terms of Attribution 4.0 International Creative Commons License (http://creativecommons.org/licenses/by/4.0/). The Data description is under review.

## Description
Recently, practical brain-computer interfaces (BCIs) have been widely investigated for detecting human intentions in real world. However, performance differences still exist between the laboratory and the real world environments. One of the main reasons for such differences comes from the user's unstable physical states (e.g., human movements are not strictly controlled), which produce unexpected signal artifacts. Hence, to minimize the performance degradation of electroencephalography (EEG)-based BCIs, we present a novel artifact removal method named constrained independent component analysis with online learning (cIOL). The cIOL can find and reject the noise-like components related to human body movements (i.e., movement artifacts) in the EEG signals. To obtain movement information, isolated electrodes are used to block electrical signals from the brain using high-resistance materials. We estimate artifacts with movement information using constrained independent component analysis from EEG signals and then extract artifact-free signals using online learning in each sample. In addition, the cIOL is evaluated by signal processing under 16 different experimental conditions (two types of EEG devices × two BCI paradigms × four different walking speeds). The experimental results show that the cIOL has the highest accuracy in both scalp- and ear-EEG, and has the highest signal-to-noise ratio in scalp-EEG among the state-of-the-art methods, except for the case of steady-state visual evoked potential at 2.0 m/s with superposition problem.

## Code List
- method_execute: main script of the cIOL. You need to input noisy signals and reference signals of noise. 
- cIOL: main function of cIOL. In this function, it will execute cICA and Adaptive Filter in sequence.
- cICA_multiCh: cICA algorithm managing multi-channel.
- plot_each_channel: plotting each channels separately in a figure.

### Developing Environment
- Matlab 2019b
