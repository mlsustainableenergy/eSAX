## Energy Time Series Motif Discovery using Symbolic Aggregated Approximation (eSAX)

The scripts were used to perform eSAX for the paper:
*Karrari et al (2019) "A method for sizing centralised energy storage systems using standard patterns"*. All code in R. A python version of the code can be found [here](https://github.com/KIT-IAI/eSAX).

Note: The same time series data is here used with different aggregations, i.e. data recorded in 1s, 10min etc intervals

Assumed structure of files: all input time series are in a folder data with one folder for each aggregation e.g. `data/data_10 min`. The code is
on the same level as the data folder.

### Files
`02_get_subsequences`: Extract subsequences from a time series and calculate ecdf of the data. Needs a data frame `data` as input where the first column is a time stamp and the second column the variable of interest e.g. power. Extracts the sequences from this data frame. There's a daily window implemented and no event detection, if you want to change that, change the input to the `minimum.search` function in lines 66ff. The result of this part is a list of all subsequences and the ecdf of the time series. Additionally the code plots the whole time series, the ecdf function, as well as all subsequences.
There's also an option to calculate the peak power of the sequences (currently this is set to FALSE).

`03_get_motif`: Extracts eSAX based motifs from the subsequences. Needs the ecdf function and the time series subsequences as input, if the previous script was used they are called ecdf_[aggregation] and tssubs_[aggregation] and in the folder data. There are several parameters in this file which can change the outcome of the motif search. There are three categories of parameters: parameters for the alphabet distribution in eSAX, parameters for random projection and parameters for motif discovery. The defaults I mostly worked with are listed in the following table.

| Parameter       | Default       | Description  |
| --------------- |:------------- | ------------ |
| breaks          | 10            | number of breakpoints in alphabet, 10 = all quantiles of the ecdf |
| w               | median(length(subsequences))      | word size (always the same if sequences are of equal length) |
| iter            | min(max(length(subsequences)*0.1), w/10)     | number of iterations of the random projection algorithm (note: motif candidate search depends on it together with count.ratio.1) |
| mask.size       | 2     | mask size for random projection |
| max.dist.ratio  | 2.5      | final distance allowed between occurrences in one motif|
| count.ratio.1   | 2.5      | controls when entries in the collision matrix become candidate motifs  |
| count.ratio.2   | 1.5      | controls whether a candidate motif becomes a motif  |

Results in a list (found.motifs) with all raw subsequences (Subs), all SAX representations of the sequences (Subs.SAX), the raw motifs (Motif.raw), the SAX representations of the motifs (Motif.sax), the collision matrix (Collision.matrix), indices of the time steps where the motifs occur (Indices) and the piecewise approximation (pieces, equal to Subs if the length of the subsequences is equal).


`04_get_plots`: plots all occurrences of found motifs. Needs the list found.motifs as input.

`05_get_profiles`: Calculates the "standard shape" of a motif. currently q20, mean, median and q80 are plotted and stored to a csv file.

© Nicole Ludwig 2021
