# Explanation for the 2D CA-CFAR algorithm

## Initialization

Three matrices with the same size of the Range Doppler map matrix *RDM* are initialized:
- *noise_level* which will store the average noise in dB on the training cells around the cell under test
- *threshold_cfar* which is equal to the noise level offset by the Signal-Noise-Ratio in dB
- *sig_cfar* which contains the output of the CA-CFAR algorithm: 0 if no detection, 1 if detection in cell under test

By initializing the matrices at zero, we automatically ensure that the non-thresholded cells at the edges are suppressed.

## Loop definition

The CA-CFAR loop will run on the cells in the *RDM* matrix:
- along the range axis we will discard from both edges the Training and Guard cell sizes for range: *Tr* and *Gr*
- along the Doppler axis we will discard from both edges the Training and Guard cell sizes for speed: *Td* and *Gd*

By cutting the edges, we will only use cells that are surrounded by enough Training and Guard Cells.

## Selection of Training and Guard cells

For each cell under test (*i,j*), we select two grids centered in (*i,j*):
- a first grid containing the training, guard and test cells, of sizes (2*Tr* + 2*Gr* + 1) on the range axis and (2*Td* + 2*Gd* + 1) on the Doppler axis.
- a second grid containing just the guard and test cells, of sizes (2*Gr* + 1) on the range axis and (2*Gd* + 1) on the Doppler axis.
Then we subtract the second grid from the first one to select only the training cells.

To perform the subtraction, we convert grid indices to linear indices using the function *sub2ind*, then perform a set difference between the linear indices set using the function *setdiff*.

## Calculation of CA-CFAR output

First we average the noise across training cells selected at the previous step and convert it in dB.

Then we add the SNR offset and put the output into the *threshold_cfar* matrix.

Finally, we check whether the Range-Doppler Map output at the cell under test, RDM, is greater than the threshold: if so, we assign a value 1 to the CA-CFAR output *sig_cfar*.

## Suppression of non-thresholded cell at the edges

As mentioned previously, the cells at the edges that were not thresholded by CA-CFAR are already set to zero by the initialization of the CFAR output matrix *sig_cfar*.
