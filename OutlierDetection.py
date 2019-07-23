# Nil Sahin, December 2017

from OutlierDetection_Functions import *
from optparse import OptionParser

# Ignore warnings
os.system('export PYTHONWARNINGS="ignore"')

# Parameters to specify
parser = OptionParser()
parser.add_option('-f', '--input-file', dest='input_file', default='',
                  help='Input cell data.')
parser.add_option('-r', '--data-type', dest='data_type', action='store_true',
                  help='Use this flag if input data is raw. Do not process as a CombineCellProfiler output file.')
parser.add_option('-i', '--identifier', type = 'str', dest = 'identifier', default = 'ORF',
                  help = 'Unique strain identifier: ORF - Name - Allele - StrainID - Gene - Genotype')
parser.add_option('-s', '--features_file', default='',
                  help='CellProfiler feature sets to include. This file is required if input data is raw.')
parser.add_option('-l', '--locations_file', default='',
                  help='Location features. This file is required if input data is raw.')
parser.add_option('-o', '--output_pca', action='store_true', 
                  help='Save PCA results')
parser.add_option('-n', '--negcontrol-file', type = 'str', dest = 'negcontrol_file', default = '',
                  help = 'Negative control file')
parser.add_option('-p', '--poscontrol-file', type = 'str', dest = 'poscontrol_file', default = '',
                  help = 'Positive control file')
parser.add_option('-v', '--variance', type = 'float', dest = 'var',
                  default = 0.80, help = 'Variance explained by PCA, min = 0, max = 1')
parser.add_option('-t', '--threshold', type = 'float', dest = 'thres',
                  default = 10, help = 'Fraction of WT outliers in a population, min = 0, max = 100')
(options, args) = parser.parse_args()

# Options
filename = options.input_file
rawdata = options.data_type
identifier = options.identifier
features_file = options.features_file
locations_file = options.locations_file
save_pca = options.output_pca
neg_control_file = options.negcontrol_file
pos_control_file = options.poscontrol_file
variance = options.var
outlier_threshold = options.thres

# Screen name from the input filename
screen_name = filename.split('/')[-1][:-4]

# Get negative control strains
wt_strains = [x.strip() for x in open(neg_control_file, 'r').readlines()]

if __name__ == '__main__':
    # Scale data for each plate and reduce dimensions
    df, plates, CP_features, identifier_features, location_features = extract_plate_information(filename, screen_name,
                                                                                                wt_strains, identifier,
                                                                                                rawdata, features_file,
                                                                                                locations_file)

    # Name output files
    output_files = prepare_output_filenames(screen_name)

    # Dimensionality Reduction
    df = do_PCA(df, output_files, variance, CP_features, save_pca, identifier_features, location_features)

    # Outlier detection
    df = OneClassSVM_method(df, output_files, outlier_threshold, identifier_features, location_features)

    # Prepare penetrance files
    df_OUT = prepare_output_well(df, plates, output_files, rawdata, identifier_features, location_features)
    df_OUT_strain = prepare_output_strain(df, identifier, output_files, rawdata, identifier_features, location_features)
    plot_heatmaps(df_OUT, plates, location_features, output_files)

    # Plot performance results
    plot_performance(pos_control_file, df_OUT_strain, wt_strains, identifier, identifier_features, output_files)
