# ocSVM
Outlier Detection with One-Class SVM Method

## Prerequisites
Python 2.7 or 3: https://www.python.org/downloads

## Installation and Usage
Clone the repository
```
git clone https://github.com/BooneAndrewsLab/ocSVM.git
cd ocSVM
```

Create a new conda environment with the dependencies (recommended):
```
conda env create -f environment.yml
conda activate ocSVM
```

or if using virtualenv or the system wide Python, do
```
pip install -r requirements.txt
```

Usage
```
python OutlierDetection.py -f <input_file> <data_type> -i <identifier> 
-s <features_file> -l <locations_file> -n <negcontrol_file> -p <poscontrol_file>
-v <variance> -t <threshold> <save_pca> <heatmap> 
```

## Input Options

**--input-file (-f)**: This is the CellProfiler features file for each cell that 
also includes the strain information.

_Examples:_
* _cp_example/input/CDC11_Rep1.csv_
* _raw_example/input/DMA_TSA.csv_

**--data-type (-r)**: Use this flag if input data is raw*. Do not process input 
as CombineCellProfiler* output file.

_*raw = raw_example/input/DMA_TSA.csv_  
_*CombineCellProfiler = cp_example/input/CDC11_Rep1.csv_


**--identifier (-i)**: The information from the mapping sheet that describes a genotype in the 
most unique way. SGA mapping sheets identify strains with ORF - Name - Allele - StrainID format.
Multiple alleles can have the same ORF so an ORF may not be a unique strain identifier. 
Recommendation: try to pick Genotype or StrainID

**--features-file (-s)**: CellProfiler feature sets to include. This file is required 
if input data is raw.

_Example: raw_example/input/features.txt_

**--locations-file (-l)**: Location features. This file is required if input data is raw.

_Example: raw_example/input/features_location.txt_

**--negcontrol-file (-n)**: Negative control file. A txt file that contains all the strain information with the 
same unique strain identifier option -i


_Examples:_
* _cp_example/input/negative_controls.txt_
* _raw_example/input/negative_controls.txt_

**--poscontrol-file (-n)**: Positive control file. This file is not required but only to 
produce the performance of outlier detection. A csv file that contains a list of strain 
information. The strain information should be present in the column as any subset of the 
following strain information: ORF, Name, Allele, StrainID, Gene, Genotype. 

May also include "Penetrance_bins" column if the expected penetrance is estimated for 
positive control strains. 

_Examples:_
* _cp_example/input/positive_controls.csv_
* _cp_example/input/positive_controls_bins.csv_

**--variance (-v)**: Variance explained by PCA, min = 0, max = 1. Default is 0.80

**--threshold (-t)**: Fraction of WT outliers in a population, min = 0, max = 100. Default is 10

**--savepca (-o)**: Use this flag to save PCA results

**--heatmap (-m)**: Use this flag if input data has Row and Column information. This will 
generate a heatmap representing the penetrance values for each plate's wells

## Examples

There are two sets of available example dataset:

_cp_example_: CellProfiler data file <br />
_raw_example_: raw input data file

