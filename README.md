# SampEnMF_Covid
Implementation in Matlab for method SampEnMF apllied in radiographic images for diagnosis of covid.

The file that must be called is main.m. From this file the others will be executed automatically.

The main parameters used were:
m=1:4
e (tolerance constant)=0.06:0.02:0.40
limit(number of windows)=180

General observations:
The algorithm is adapted to color images, which is not the case in this database.
Each base image is divided into sub-images of 64x64 pixels. These sub-pictures are centered according to the original image for best use. Entropy is calculated for each sub-image that contains a proportion of at most 90% black pixels (background).
Comparisons between windows (which are drawn) are performed on each sub-image and are following the rule of at least 98% non-black pixels.
The windows draw is carried out in order to guarantee that the windows are obtained from several parts of the sub-image, for this the sub-images were divided into 9 regions and the windows are drawn from these.
Pixels with 16 bits and signals (int16 type) were used for comparisons, with intensities ranging from -32768 to 32767. These values provided the best results in tests with the colorectal base and the fuzzy function is adapted to these values.
An entropy value is obtained for each sub-image and at the end an average of the sub-images considered (which were not discarded in the process) is obtained.
