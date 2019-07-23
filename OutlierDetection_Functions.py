import os
import re
import sys
import time
import itertools
import matplotlib
import numpy as np
import pandas as pd
import seaborn as sns
from ppca import PPCA
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats, spatial
from sklearn.decomposition import PCA
from sklearn import svm, metrics, mixture


def prepare_output_filenames(screen_name):
    """ Prepare output filenames specific for each timepoint.

            Args:
                screen_name:    Extracted from the input filename

            Return:
                output:         Output filenames
    """

    folder_name = '%s_OD_results' % (screen_name)
    if not os.path.isdir(folder_name):
        os.system('mkdir %s' % folder_name)
    log = open('%s/%s_log.txt' % (folder_name, screen_name), 'w')
    log.write('Reading input files\n')
    log.close()

    output = {'HeatmapPrefix': '%s/%s_plate' % (folder_name, screen_name),
              'ODresultsWell': '%s/%s_OD_results_well.csv' % (folder_name, screen_name),
              'ODresultsStrain': '%s/%s_OD_results_strain.csv' % (folder_name, screen_name),
              'PCAResult': '%s/%s_PCA.csv' % (folder_name, screen_name),
              'PCACellsPNG': '%s/%s_PCA_cells.png' % (folder_name, screen_name),
              'PCACellsSVG': '%s/%s_PCA_cells.svg' % (folder_name, screen_name),
              'PCAExplainedVariance': '%s/%s_PCA_explained_variance.txt' % (folder_name, screen_name),
              'PCAExplainedVariancePlot': '%s/%s_PCA_explained_variance_plot.png' % (folder_name, screen_name),
              'PCAFeatureCorrelations': '%s/%s_PCA_feature_correlations.csv' % (folder_name, screen_name),
              'ScoreHistogramPNG': '%s/%s_OD_results_score_histogram.png' % (folder_name, screen_name),
              'ScoreHistogramSVG': '%s/%s_OD_results_score_histogram.svg' % (folder_name, screen_name),
              'ScoreCells': '%s/%s_OD_results_score.csv' % (folder_name, screen_name),
              'ROCCurve': '%s/%s_OD_results_ROC_curve.png' % (folder_name, screen_name),
              'PRCurve': '%s/%s_OD_results_PR_curve.png' % (folder_name, screen_name),
              'CurveNumbers': '%s/%s_OD_results_ROC_PR_curve_numbers.csv' % (folder_name, screen_name),
              'PenetranceBins': '%s/%s_OD_results_penetrance_bins.csv' % (folder_name, screen_name),
              'PenetranceAgreement': '%s/%s_OD_results_penetrance_agreement.png' % (folder_name, screen_name),
              'ConfusionMatrix': '%s/%s_OD_results_confusion_matrix.png' % (folder_name, screen_name),
              'log': '%s/%s_log.txt' % (folder_name, screen_name),
              }

    return output


def log_write(log_f, text):
    """ Open and write to the log file

            Args:
                log_f:  Log file
                text:   Something to write to log file
            """
    
    f = open(log_f, 'a')
    f.write(text)
    f.close()


def add_phenix_filenames(df, fields=4, timepoints=1):
    """ Add phenix filenames for each image

        Well from Plate - 1, Row - 1, Column - 1, Timepoint - 1 will generate these files:
        r01c01f01p01-ch1sk1fk1fl1.tiff
        r01c01f02p01-ch1sk1fk1fl1.tiff
        r01c01f03p01-ch1sk1fk1fl1.tiff
        r01c01f04p01-ch1sk1fk1fl1.tiff

            Args:
                df:         Combined dictionary
                timepoints: Number of timepoints to analyze

            Return:
                df:         Combined dictionary with OperaPhenix filename structure
            """

    num_fields = fields
    num_timepoints = timepoints

    df['Plate'] = df['Plate'].astype(str)
    df['Frame'] = 1
    df['Time'] = 1
    add_df = pd.DataFrame(columns=df.columns.values)
    for i in range(num_timepoints):
        for j in range(num_fields):
            add_this = df.copy()
            add_this['Frame'] = j + 1
            add_this['Time'] = i + 1
            add_df = pd.concat([add_df, add_this])
    df = add_df.sort_values(['Plate', 'Row', 'Column'], ascending=True)
    df = df.reset_index(drop=True)

    df['Filename'] = ''
    for i in range(len(df)):
        # Get well information
        r = '0' + str(int(df.iloc[i, df.columns.get_loc('Row')]))
        c = '0' + str(int(df.iloc[i, df.columns.get_loc('Column')]))
        f = str(int(df.iloc[i, df.columns.get_loc('Frame')])).zfill(2)
        #p = '0' + str(int(df.iloc[i, df.columns.get_loc('Plate')]))
        p = '01'
        sk = str(int(df.iloc[i, df.columns.get_loc('Time')]))

        # Correct names
        if len(r) == 3:
            r = r[1:]
        if len(c) == 3:
            c = c[1:]
        if len(p) == 3:
            p = p[1:]

        # Construct OperaPhenix filenames
        filename = 'r%sc%sf%sp%s-ch1sk%sfk1fl1.tiff' % (r, c, f, p, sk)
        df.iloc[i, df.columns.get_loc('Filename')] = filename

    return df


def combine_sql(path, mapping_sheet, fields, timepoints):
    """ Load up the SQL_ and support files, merge the image and object and write the final table

            Args:
                path:           The folder that contains CellProfiler SQL files
                mapping_sheet:  Mapping sheet file
                fields:         Number of fields to analyze
                timepoints:     Number of timepoints to analyze

            Return:
                CP_out:         SQL and mapping sheet merged for all SQL files
    """

    print('\nCombining CellProfiler SQL files...')

    # Headers
    ImageHeader, ObjectHeader = parse_sql_setup(path + 'SQL_SETUP.SQL')

    # Genes
    Df_GeneNames = pd.read_csv(mapping_sheet, index_col=None)
    Df_GeneNames = add_phenix_filenames(Df_GeneNames, fields=fields, timepoints=timepoints)

    # MetaDataToKeep
    MetaDataKeep = ['Plate', 'ImageNumber', 'Image_Count_Cells', 'Image_Count_Nuclei',
                    'Image_FileName_GFP', 'Image_FileName_RFP']
    MetaDataKeep += Df_GeneNames.columns.values.tolist()

    # Create a list of the SQL files to look through
    files = os.listdir(path)
    files_sql_prefixs = [m.group(1) for x in files for m in [re.search('(SQL_[0-9]*_[0-9]*)_', x)] if m]
    files_sql_prefixs = sorted(list(set(files_sql_prefixs)), key=lambda item: int(item.split('_')[1]))

    # Combine SQL files together in a dataframe
    CP_OUT = pd.DataFrame()
    for f in files_sql_prefixs:
        df_sql = merge_image_object_file(path, f, Df_GeneNames, MetaDataKeep, ImageHeader, ObjectHeader)
        CP_OUT = pd.concat([CP_OUT, df_sql])

    # Add keywords to feature names
    columns = CP_OUT.columns.values.copy()
    cpp = ['AreaShape', 'Granularity', 'Intensity', 'RadialDistribution', 'IntensityDistribution', 'Texture']
    for i in range(len(columns)):
        for feature in cpp:
            if (columns[i] not in MetaDataKeep) and ('Center' not in columns[i]) and (feature in columns[i]):
                columns[i] = 'CP_' + columns[i]
                break

        if columns[i] in Df_GeneNames.columns.values:
            columns[i] = 'MS_' + columns[i]

    CP_OUT.columns = columns
    CP_OUT.reset_index(drop=True)

    return CP_OUT


def merge_image_object_file(path, f, Df_GeneNames, MetaDataKeep, ImageHeader, ObjectHeader):
    """ Merges an image and object file from the SQL prefix f

            Args:
                path:           The folder that contains CellProfiler SQL files
                f:              The specific SQL to read
                Df_GeneNames:   Mapping sheet data
                MetaDataKeep:   The columns to keep for the combined data
                image_header:   Image file header
                object_header:  Object file header

            Return:
                merged:         SQL and mapping sheet merged
    """

    # Load image file
    infile_image = path + f + '_Image.CSV'
    Df_Image = pd.read_csv(infile_image, header=None, index_col=None)
    
    # Determine if this is a large or small file and attach header
    if Df_Image.shape[1] == len(ImageHeader):
        Df_Image.columns = ImageHeader
    else:
        sys.exit('File does not match either the standard or alt headers')

    # Merge
    Df_Image['Filename'] = Df_Image['Image_FileName_GFP'].values.tolist()
    Df_Image_Genes = pd.merge(Df_GeneNames, Df_Image, on=['Filename'])
    assert Df_Image_Genes.shape[0] > 0, 'Failed to merge GeneName information on'
    
    # Reduce to features of interest
    Df_Image_Genes_Sub = Df_Image_Genes[[x for x in Df_Image_Genes.columns if x in MetaDataKeep]]

    # Load object file
    infile_object = path + f + '_Object.CSV'
    Df_Object = pd.read_csv(infile_object, header = None,index_col = None)

    # Determine if this is a large or small file and attach header
    if Df_Object.shape[1] == len(ObjectHeader):
        Df_Object.columns = ObjectHeader
    else:
        sys.exit('File does not match either the standard or alt headers')

    # Merge image and object file
    return pd.merge(Df_Image_Genes_Sub, Df_Object, on='ImageNumber')


def parse_sql_setup(infile):
    """ Reads the SQL_SETUP.SQL file for Image Object file headers

            Args:
                infile:         SQL_SETUP.SQL

            Return:
                image_header:   Image file header
                object_header:  Object file header
    """

    f = open(infile)
    r = {}
    r['image'] = []
    r['object'] = []
    switch = 0
    type_ = ''
    for line in f:
        line = line.rstrip('\n')
        if 'CREATE TABLE' in line:
            if 'Image' in line:
                switch = 1
                type_ = 'image'
            elif 'Object' in line:
                switch = 1
                type_ = 'object'
        elif 'PRIMARY KEY' in line:
            switch = 0
        elif switch == 1:
            val = line.split(' ')[0]
            val = val.replace(',','')
            r[type_].append(val)

    return r['image'], r['object']


def extract_plate_information(filename, screen_name, wt, identifier, rawdata, features_file, locations_file):
    """ Extract information from the input file
        Combines single cell information in a dictionary

            Args:
                filename:               CP output as input file
                screen_name:            Extracted from the input filename
                wt:                     WT names
                identifier:             Unique identifier

            Return:
                df:                     Combined dictionary
                plates:                 List of plate names from the screen
                CP_features:            List of CP feature names
                identifier_features:    List of strain identifiers
                location_features:      List of Plate - Row - Column - Filename
    """

    # Read input file
    print('\nAnalyzing screen: %s' % screen_name)

    # Load CP features data
    input_df = pd.read_csv(filename)

    CP_features = []
    identifier_features = []
    location_features = []
    if not rawdata:
        # Extract features with keywords
        for column in input_df.columns.values:
            if column[:3] == 'CP_':
                CP_features.append(column)
            if column[:3] == 'MS_':
                if column[3:].lower() in ['plate', 'row', 'column', 'time']:
                    location_features.append(column[3:])
                elif column[3:].lower() not in ['frame', 'filename']:
                    identifier_features.append(column[3:])
    else:
        # Get a feature list either from a file or the data itself (method 2)
        if features_file == '':
            raise ValueError('Features file is required if input data is raw.')
        if locations_file == '':
            raise ValueError('Locations file is required if input data is raw.')
        else:
            f = open(features_file, 'r')
            CP_features = list(filter(None, [x.strip() for x in f.readlines()]))
            f.close()

            l = open(locations_file, 'r')
            other_features = list(filter(None, [x.strip() for x in l.readlines()]))
            l.close()
        for c in other_features:
            if c.lower() != identifier.lower():
                location_features.append(c)
            else:
                identifier_features.append(c)

            if 'plate' in c.lower():
                plate_identifier = c

    # Fix cells with empty rows and corrupt nan entries
    input_df[CP_features] = input_df[CP_features].replace(to_replace=r'[a-zA-Z]+', value=np.nan, regex=True)
    nan_count = np.asarray(input_df[CP_features].isnull().sum(axis=1))
    input_df = input_df.iloc[nan_count != len(CP_features), :]
    input_df = input_df.reset_index(drop=True)

    # Get plate and timepoint names
    if not rawdata:
        plates = input_df.MS_Plate.unique()
    else:
        plates = input_df[plate_identifier].unique()
 
   # Create combined dictionary
    df = {}
    for f in identifier_features:
        if not rawdata:
            df[f] = np.array(input_df['MS_' + f])
        else:
            df[f] = np.array(input_df[f])
    for f in location_features:
        if not rawdata:
           if f.lower() in ['row', 'column']:
               df[f] = np.array(input_df['MS_' + f], dtype='int32')
           else:
               df[f] = np.array(input_df['MS_' + f])
        else:
            if 'plate' in f:
                df['Plate'] = np.array(input_df[f])
            else:
                df[f] = np.array(input_df[f])
    
    df['Data'] = np.array(input_df[CP_features], dtype='float64') # (input_df[CP_features], dtype='float64')
    df['DataScaled'] = np.array(input_df[CP_features], dtype='float64')

    # Mark cells that are WT
    df['Mask_WT'] = np.array([x in wt for x in df[identifier]])

    for p in plates:
        mask_plate = np.array([x == p for x in df['Plate']])
        # Mask the plate and wt rows to calculate mean and sd
        data_WT = df['Data'][(mask_plate == 1) & (df['Mask_WT'] == 1)]
        if len(data_WT) == 0:
            data_WT = df['Data'][mask_plate == 1]
        mu = np.nanmean(data_WT, axis=0)
        sd = np.nanstd(data_WT, axis=0)
        df['DataScaled'][mask_plate == 1] = (df['Data'][mask_plate == 1] - mu) / sd
    return df, plates, CP_features, identifier_features, location_features


def do_PCA(df, output, var, feature_set, save_pca, identifier_features, location_features):
    """ Perform PCA or probabilistic PCA depending on nan values.

            Args:
                df:             Existing combined dictionary
                output:         Output filenames
                var:            Minimum explained variance required
                feature_set:    List of CP feature names

            Return:
                df:             Updated with added 'DataPCA' values
            """

    print('Do PCA...')
    log_write(output['log'], 'Do PCA...\n')

    # Check whether there are nan values and choose PCA method
    check_nan = np.isnan(df['DataScaled']).any()
    if check_nan:
        df, exp_var, num_PCs, PCA_feature_loadings = do_probabilistic_PCA(df, var, output)
    else:
        df, exp_var, num_PCs, PCA_feature_loadings = do_regular_PCA(df, var, output)
    
    if save_pca:
        pca_df = pd.DataFrame()
        for f in identifier_features:
            pca_df[f] = df[f]
        for f in location_features:
            if 'plate' not in f.lower():
                pca_df[f] = df[f]
            else:
                pca_df['Plate'] = df['Plate']
        for i in range(num_PCs):
            pca_df['PC%d' %(i+1)] = df['DataPCA'][:,i]
        pca_df.to_csv(output['PCAResult'], index=False)

    # Save correlation matrix of PCs and features
    save_PCA_feature_correlation(df, num_PCs, feature_set, output)

    # Plot total explained variance with each added PC
    plt.plot(exp_var)
    plt.xlabel('Number of PCs')
    plt.ylabel('Total % of variance explained')
    plt.title('Number of PCs to be used = %d / %d features' % (num_PCs, df['Data'].shape[1]))
    fig = plt.gcf()
    fig.savefig(output['PCAExplainedVariancePlot'])
    fig.clf()
    plt.close(fig)

    return df


def do_regular_PCA(df, var, output):
    """ Perform regular PCA on scaled values for the whole screen

            Args:
                df:             Existing combined dictionary
                var:            Minimum explained variance required
                output:         Output filenames

            Return:
                df:             Updated with added 'DataPCA' values
                exp_var:        List of explained variance with each added PC
                num_PCs:        Number of PCs to explain var
                PCA_loadings:   Principal axes in feature space (n_components, n_features)
    """

    print('Feature selection using regular PCA...')
    log_write(output['log'], 'Feature selection using regular PCA...\n')

    exp_var = []
    num_PCs = 0

    # Do PCA with number of components iteratively (max is the number of features)
    for i in range(df['DataScaled'].shape[1]):
        pca = PCA(n_components=i)
        pca.fit(df['DataScaled'])
        total_var = sum(pca.explained_variance_ratio_)
        exp_var.append(total_var)
        # End PCA if the total variance passes the minimum variance required
        if total_var > var:
            num_PCs = i
            np.savetxt(output['PCAExplainedVariance'], pca.explained_variance_ratio_, fmt='%0.4f')
            break

    # Do the final PCA with num_PCs
    pca = PCA(n_components=num_PCs)
    df['DataPCA'] = pca.fit_transform(df['DataScaled'])
    PCA_loadings = pca.components_

    return df, exp_var, num_PCs, PCA_loadings


def do_probabilistic_PCA(df, var, output):
    """ Perform probabilistic PCA (PPCA) on scaled values for the whole screen

            Args:
                df:             Existing combined dictionary
                var:            Minimum explained variance required
                output:         Output filenames

            Return:
                df:             Updated with added 'DataPCA' values
                exp_var:        List of explained variance with each added PC
                num_PCs:        Number of PCs to explain var
                PCA_loadings:   Principal axes in feature space (n_components, n_features)
    """

    print('Feature selection using probabilistic PCA...')
    log_write(output['log'], 'Feature selection using probabilistic PCA...\n')

    # Initialize parameters
    exp_var = [0]
    exp_var_ratio = []
    num_PCs = 0
    ppca = PPCA()
    ppca.fit(df['DataScaled'], d=2)
    exp_var.append(ppca.var_exp[0])
    exp_var_ratio.append(ppca.var_exp[0])

    # Do PPCA with number of components iteratively (max is the number of features, min is 2)
    for i in range(2, df['DataScaled'].shape[1]):
        num_PCs = i
        ppca = PPCA()
        ppca.fit(df['DataScaled'], d=i)
        total_var = ppca.var_exp[i-1]
        exp_var.append(total_var)
        exp_var_ratio.append(ppca.var_exp[i-1] - ppca.var_exp[i-2])
        # End PCA if the total variance passes the minimum variance required
        if total_var > var:
            num_PCs = i
            np.savetxt(output['PCAExplainedVariance'], exp_var_ratio, fmt='%0.4f')
            break

    # Do the final PCA with num_PCs
    ppca = PPCA()
    ppca.fit(df['Data'], d=num_PCs)
    df['DataPCA'] = ppca.transform()
    PPCA_loadings = np.transpose(ppca.C)

    return df, exp_var, num_PCs, PPCA_loadings


def save_PCA_feature_correlation(df, num_PCs, feature_set, output):
    """ Calculate and save correlation between PCs and raw feature data.

            Args:
                df:             Existing combined dictionary
                num_PCs:        Number of PCs to explain var
                feature_set:    List of CP feature names
                output:         Output filenames
    """

    log_write(output['log'], 'Saving PCA feature correlation...\n')

    # Initialize dataframe
    pca_feat_corr = pd.DataFrame(columns=feature_set)
    PCA_columns = []
    for i in range(num_PCs):
        PCA_columns.append('PC' + str(i + 1))
        corr = []
        for j in range(len(feature_set)):
            data_raw = np.copy(df['Data'][:, j])
            data_pca = np.copy(df['DataPCA'][:, i])
            # Remove nan values from the data
            where_nan = np.isnan(data_raw)
            data_raw = data_raw[~where_nan]
            data_pca = data_pca[~where_nan]
            corr.append(stats.pearsonr(data_raw, data_pca)[0])
        pca_feat_corr.loc[i,] = corr

    pca_feat_corr = pca_feat_corr.set_index([PCA_columns])
    pca_feat_corr.to_csv(path_or_buf=output['PCAFeatureCorrelations'])


def OneClassSVM_method(df, output, out_threshold, identifier_features, location_features):
    """ Outlier Detection with One-Class SVM Method.

            Args:
                df:             Existing combined dictionary
                output:         Output filenames
                out_threshold:  WT threshold on the right tail to decide on outlier boundary

            Return:
                df:             With in-outlier information
    """

    log_write(output['log'], 'Outlier detection...\n')
    start_time = time.time()

    # Create a subset with only WT cells and fit the model
    ocsvm = svm.OneClassSVM(kernel='rbf', nu=out_threshold/100.0)
    ocsvm.fit(df['DataPCA'][df['Mask_WT'] == 1])
    dist_to_border = - ocsvm.decision_function(df['DataPCA']).ravel()

    # Threshold and plot data
    threshold = 0
    df['Is_Inlier'] = dist_to_border <= threshold
    plot_title = 'Outlier Detection results after thresholding'
    plot_in_outliers(df['DataPCA'], df['Is_Inlier'], plot_title, output)

    # Plot the distances as an histogram
    plot_title = 'OC-SVM: Distance to boundary'
    plot_score_histogram(df['Is_Inlier'], dist_to_border, threshold, plot_title, output)

    # Save scores
    scores_df = pd.DataFrame()
    for f in identifier_features:
        scores_df[f] = df[f]
    for f in location_features:
        if 'plate' not in f.lower():
            scores_df[f] = df[f]
        else:
            scores_df['Plate'] = df['Plate']
    scores_df['score'] = dist_to_border
    scores_df['is_inlier'] = df['Is_Inlier']
    scores_df.to_csv(output['ScoreCells'], index=False)

    # Print OD runtime
    text = 'Outlier detection method: One-Class SVM\n'
    text += 'WT outlier threshold: %.2f%%\n' % out_threshold
    text += 'Outlier detection runtime: %.2f minutes\n' % ((time.time()-start_time)/60.0)
    text += 'Number of samples: %d\n' % dist_to_border.shape[0]
    text += 'Number of negative samples: %d\n' % dist_to_border[df['Mask_WT'] == 1].shape[0]
    log_write(output['log'], text)

    return df


def plot_in_outliers(data, mask, title, output):
    """ Plot data with in-outlier information using the first 2 PCs

            Args:
                data:           PCA Data to plot
                mask:           Mask the in-outliers
                title:          Title for the plot
                output:         Output filenames
            """

    oc = 'lightskyblue'
    ic = 'navy'
    plt.figure(figsize=(15, 18))
    sns.set_style('white')
    x_all = pd.DataFrame({'PC1': data[:, 0], 'PC2': data[:, 1]})
    x_inliers = pd.DataFrame({'PC1': data[mask == 1, 0], 'PC2': data[mask == 1, 1]})
    x_outliers = pd.DataFrame({'PC1': data[mask == 0, 0], 'PC2': data[mask == 0, 1]})
    # Plot everything first
    g = sns.JointGrid(x='PC1', y='PC2', data=x_all)
    # Plot points
    sns.scatterplot(x_outliers.PC1, x_outliers.PC2, color=oc, ax=g.ax_joint,
                    s=10, linewidth=0, label='Mutant morphology')
    sns.scatterplot(x_inliers.PC1, x_inliers.PC2, color=ic, ax=g.ax_joint,
                    s=10, linewidth=0, label='Normal morphology')
    # Plot kernel density estimates
    sns.distplot(x_outliers.PC1, kde=True, hist=False, color=oc, ax=g.ax_marg_x, axlabel=False)
    sns.distplot(x_inliers.PC1, kde=True, hist=False, color=ic, ax=g.ax_marg_x, axlabel=False)
    sns.distplot(x_outliers.PC2, kde=True, hist=False, color=oc, ax=g.ax_marg_y,
                 vertical=True, axlabel=False)
    sns.distplot(x_inliers.PC2, kde=True, hist=False, color=ic, ax=g.ax_marg_y,
                 vertical=True, axlabel=False)
    fig = plt.gcf()
    plt.title(title, y=1.2)
    fig.savefig(output['PCACellsPNG'], dpi=150, bbox_inches='tight')
    fig.savefig(output['PCACellsSVG'])
    fig.clf()
    plt.close(fig)


def plot_score_histogram(mask, dist, out_threshold, xlabel, output):
    """ Plot distances as an histogram from the complete data

            Args:
                mask:           Mask the in-outliers
                dist:           Distance calculated by the outlier detection algorithm
                out_threshold:  To decide in-outliers
                xlabel:         Distance label
                output:         Output filenames
            """

    oc = 'lightskyblue'
    ic = 'navy'
    plt.figure(figsize=(8, 6))
    plt.axvline(x=out_threshold, color='r', linestyle='-', label='Mutant threshold (%.2f)' % out_threshold)
    sns.set_style('white')
    percent_mut = sum(mask == 0) / float(mask.shape[0]) * 100
    sns.distplot(dist[mask == 0], kde=True, hist=False, color=oc, label='Mutant morphology (%.2f%%)' % percent_mut)
    sns.distplot(dist[mask == 1], kde=True, hist=False, color=ic, label='Normal morphology', axlabel=xlabel)
    fig = plt.gcf()
    plt.ylabel('Density')
    plt.xlim(min(dist), stats.scoreatpercentile(dist, 99.99))
    plt.legend(loc='upper right')
    plt.savefig(output['ScoreHistogramPNG'], dpi=150, bbox_inches='tight')
    plt.savefig(output['ScoreHistogramSVG'])
    fig.clf()
    plt.close(fig)
    

def p_value(df):
    """ Return WT cell numbers to calculate p-value.

            Args:
                df:                 Existing combined dictionary

            Return:
                WT_cells:           Number of cells in WT populations
                WT_cells_outliers:  Number of outlier cells in WT populations
            """
    WT_cells = len(df['Mask_WT'][df['Mask_WT'] == 1])
    WT_cells_outliers = len(df['Mask_WT'][(df['Mask_WT'] == 1) & (df['Is_Inlier'] == 0)])

    return WT_cells, WT_cells_outliers


def dataframe_from_dict(df, features):
    """ Create a new dataframe from the existing combined dictionary

            Args:
                df:         Existing combined dictionary
                features:   Features needed for the dataframe from the dictionary

            Return:
                new_df:     Combined dataframe with features
            """

    new_df = pd.DataFrame()
    for f in features:
        if 'plate' not in f.lower():
            new_df[f] = df[f]
        else:
            new_df['Plate'] = df['Plate']

    return new_df


def prepare_output_well(df, plates, output, rawdata, identifier_features, location_features):
    """ Prepare the output file with plate, row and column information
        Calculate penetrance and p-value

            Args:
                df:                     Existing combined dictionary
                plates:                 Plates in this screen
                output:                 Output filenames
                identifier_features:    List of strain identifiers
                location_features:      List of Plate - Row - Column - Filename

            Return:
                final_df_output:        Combined outlier detection results
            """

    print('Preparing the output values by well...')
    log_write(output['log'], 'Preparing penetrance results by well...\n')

    # Create new dataframe from dict
    append_list = identifier_features + location_features + ['Is_Inlier']
    final_df = dataframe_from_dict(df, append_list)
    if not rawdata:
        well_identifier = 'Row_Col'
    else:
        for f in location_features:
            if 'well' in f.lower():
                well_identifier = f

    if not rawdata:
        final_df[well_identifier] = final_df.Row.map(int).map(str) + '_' + final_df.Column.map(int).map(str)
    else:
        final_df[well_identifier] = final_df.well_number.map(str)

    # Initialize output folder
    final_df_output = pd.DataFrame(columns = identifier_features + location_features +
                                             ['Num_cells', 'Penetrance', 'P-value'])
    this_row = 0

    # Regroup this dataframes by plates then row column info
    WT_cells, WT_cells_outliers = p_value(df)
    plate_column = 'Plate'
    for p in plates:
        final_df_plate = final_df[final_df[plate_column] == p]
        # Regroup this dataframes by Row and Column
        row_col = final_df_plate[well_identifier].unique().tolist()
        for rc in row_col:
            df_rc = final_df_plate[final_df_plate[well_identifier] == rc]
            is_inlier_rc = np.asarray(df_rc['Is_Inlier'])
            num_cells = df_rc.shape[0]
            num_outliers = sum(is_inlier_rc == 0)
            pene = float(num_outliers) / num_cells * 100
            pval = 1 - stats.hypergeom.cdf(num_outliers, WT_cells, WT_cells_outliers, num_cells)

            # Append them to corresponding variables
            line = []
            for i in identifier_features + location_features:
                if 'plate' in i.lower():
                    i = 'Plate'
                line.append(df_rc[i].unique()[0])
            line.append(num_cells)
            line.append(pene)
            line.append(pval)
            final_df_output.loc[this_row, ] = line
            this_row += 1

    # Save into a dataframe
    final_df_output = final_df_output.sort_values('Penetrance', ascending=False)
    final_df_output = final_df_output.reset_index(drop=True)
    final_df_output.to_csv(path_or_buf=output['ODresultsWell'], index=False)

    return final_df_output


def prepare_output_strain(df, identifier, output, rawdata, identifier_features, location_features):
    """ Prepare the output file with strain information
        Calculate penetrance and p-value

            Args:
                df:                     Existing combined dictionary
                identifier:             Unique identifier
                output:                 Output filenames
                identifier_features:    List of strain identifiers
                location_features:      List of Plate - Row - Column - Filename

            Return:
                final_df_output:        Combined outlier detection results
            """

    print('Preparing the output values by strain...')
    log_write(output['log'], 'Preparing penetrance results by strain...\n')

    # Create new dataframe from dict
    append_list = identifier_features + location_features + ['Is_Inlier']
    final_df = dataframe_from_dict(df, append_list)
    if not rawdata:
        final_df['Well'] = final_df.Plate.map(int).map(str) + '_' + \
                       final_df.Row.map(int).map(str) + '_' + \
                       final_df.Column.map(int).map(str)
    else:
        final_df['Well'] = final_df.Plate.map(int).map(str) + '_' + \
                       final_df.well_number.map(str)

    # Initialize output folder
    final_df_output = pd.DataFrame(columns = identifier_features + location_features  +
                                             ['Num_cells', 'Num_wells', 'Penetrance', 'P-value'])
    this_row = 0

    # Regroup this dataframes by strainIDs
    WT_cells, WT_cells_outliers = p_value(df)
    strains = final_df[identifier].unique().tolist()
    strains = [s if not isinstance(s, float) else '' for s in strains]
    final_df[identifier] = final_df[identifier].fillna('')

    for s in strains:
        df_strain = final_df[final_df[identifier] == s]
        is_inlier_strain = np.asarray(df_strain['Is_Inlier'])
        num_cells = df_strain.shape[0]
        num_wells = len(df_strain['Well'].unique())
        num_outliers = sum(is_inlier_strain == 0)
        pene = float(num_outliers) / num_cells * 100
        pval = 1 - stats.hypergeom.cdf(num_outliers, WT_cells, WT_cells_outliers, num_cells)

        # Append them to corresponding variables
        line = []
        for i in identifier_features + location_features:
            if 'plate' in i.lower():
                i = 'Plate'
            line.append(df_strain[i].unique()[0])
        line.append(num_cells)
        line.append(num_wells)
        line.append(pene)
        line.append(pval)
        final_df_output.loc[this_row,] = line
        this_row += 1

    # Save into a dataframe
    final_df_output = final_df_output.sort_values('Penetrance', ascending=False)
    final_df_output = final_df_output.reset_index(drop=True)
    final_df_output.to_csv(path_or_buf=output['ODresultsStrain'], index=False)

    return final_df_output


def plot_heatmaps(df, plates, location_features, output):
    """ Plots penetrance values for each plate's wells if input data has row and column information

            Args:
                df:                 Existing combined dictionary
                plates:             List of plate names from the screen
                screen_name:        Extracted from the input filename
                location_features:  List of Plate - Row - Column - Filename
            """
    log_write(output['log'], 'Plotting plate penetrance heatmaps\n')

    row_column = ''
    col_column = ''
    for f in location_features:
        if f.lower() == 'row':
            row_column = f
        elif f.lower() == 'column':
            col_column = f

    for p in plates:
        penetrance_df = np.ndarray(shape=(16, 24), buffer=np.repeat(np.nan, 384))
        plate_df = df[df.Plate == p]
        for i in range(len(plate_df)):
            row = int(plate_df.iloc[i, plate_df.columns.get_loc(row_column)] - 1)
            column = int(plate_df.iloc[i, plate_df.columns.get_loc(col_column)] - 1)
            penetrance_df[row][column] = plate_df.iloc[i, plate_df.columns.get_loc('Penetrance')]
        sns.set(font_scale=1.5)
        sns.set_style()
        plt.figure(figsize=(14, 8))
        # Mask wells that have no info
        mask = penetrance_df == np.nan
        cg = sns.heatmap(penetrance_df, linewidth=.4, mask=mask, vmin=0, vmax=100, cmap='YlGnBu',
                         cbar_kws={'ticks': [0, 25, 50, 75, 100], 'label': '% penetrance'})
        cg.set_title('Penetrance')
        cg.set_xticklabels(range(1, 25))
        cg.set_yticklabels(range(1, 17))
        fig = plt.gcf()
        fig.savefig('%s%s_penetrance.png' % (output['HeatmapPrefix'], p), dpi=150, bbox_inches='tight')
        fig.clf()
        plt.close(fig)


def plot_ROC_and_PR(PC, NC, output):
    """ Plot ROC and PR curves based on Positive and Negative Control penetrance values
        Threshold changes from 100% penetrance to 0%

            Args:
                PC:       Penetrance values of positive controls
                NC:       Penetrance values of negative controls
                output:   Output filenames
            """

    # Calculate true positive rate (TPR), false positive rate (FPR), and precision
    # Recall is another name for TPR so it is calculated once
    tpr = []
    fpr = []
    prec = []
    for threshold in list(reversed(range(101))):
        # Count True and False Positives with each threshold
        tp = len(PC[PC >= threshold])
        fp = len(NC[NC >= threshold])
        tpr.append(tp / float(len(PC)))
        fpr.append(fp / float(len(NC)))
        if fp == 0:
            prec.append(1)
        else:
            prec.append(tp / float(tp + fp))
    auc = metrics.auc(fpr, tpr)

    # Save plot for ROC
    plt.figure(figsize=(6, 6))
    sns.set(font_scale=1.5)
    sns.set_style('white')
    plt.plot(fpr, tpr, color='darkorange', lw=2, label='AUROC = %0.2f' % auc)
    plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
    plt.xlim([0.0, 1.02])
    plt.ylim([0.0, 1.02])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver operating characteristic')
    plt.legend(loc='lower right')
    fig = plt.gcf()
    plt.savefig(output[0])
    fig.clf()
    plt.close(fig)

    # Save plot for PR
    plt.figure(figsize=(6, 6))
    sns.set(font_scale=1.5)
    sns.set_style('white')
    plt.plot(tpr, prec, color='darkorange', lw=2)
    plt.xlim([0.0, 1.02])
    plt.ylim([0.0, 1.02])
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title('Precision-Recall Curve')
    fig = plt.gcf()
    plt.savefig(output[1])
    fig.clf()
    plt.close(fig)

    # Save numbers
    roc_numbers = pd.DataFrame({'Penetrance_Cutoff': list(reversed(range(101))),
                                'TPR (Recall)': np.asarray(tpr),
                                'FPR': np.asarray(fpr),
                                'Precision': np.asarray(prec)})
    roc_numbers = roc_numbers[['Penetrance_Cutoff', 'TPR (Recall)', 'FPR', 'Precision']]
    roc_numbers = roc_numbers.sort_values('Penetrance_Cutoff', ascending=False)
    roc_numbers.to_csv(path_or_buf=output[2], index=False)


def confusion_matrix(actual, predicted, classes, output):
    """ Produce and plot confusion matrix

            Args:
                actual:     Actual labels
                predicted:  Predicted labels
                classes:    Class information on the confusion matrix
                output:     Output filenames
            """

    acc = metrics.accuracy_score(actual, predicted) * 100
    cm = metrics.confusion_matrix(actual, predicted)

    plt.figure(figsize=(6, 6))
    sns.set(font_scale=1.5)
    sns.set_style()
    for i in range(len(cm)):
        cm[i] = np.array(list(reversed(cm[i])))

    # Normalize confusion matrix
    cm = np.around(cm.astype('float') / cm.sum(axis=1)[:, np.newaxis], decimals=2)
    plt.imshow(cm, interpolation='nearest', cmap=plt.cm.Blues)
    plt.title('Acc %.2f%%' % acc)
    plt.colorbar()

    # Diagonal values are on y=x direction
    tick_marks = np.arange(len(classes))
    plt.xticks(list(reversed(tick_marks)), classes, rotation=45, ha='right')
    plt.yticks(tick_marks, classes)

    # Plot values on confusion matrix, change colors depending on the value
    thresh = cm.max() / 2.
    for i, j in itertools.product(range(cm.shape[0]), range(cm.shape[1])):
        plt.text(j, i, cm[i, j], horizontalalignment='center', color='white' if cm[i, j] > thresh else 'black')

    plt.ylabel('True Label')
    plt.xlabel('Predicted Label')
    plt.grid(False)
    fig = plt.gcf()
    plt.savefig(output, bbox_inches='tight')
    fig.clf()
    plt.close(fig)


def plot_penetrance_agreement(df, wt, output):
    """ Plot penetrance agreement with given penetrance bins and predicted penetrance values on a scatter plot

            Args:
                df:         Existing combined dictionary
                wt:         WT penetrance values
                output:     Output filenames
            """

    # Turn penetrance bins to average penetrance values
    # Bin-0: 0-25% penetrance
    # Bin-1: 25-50% penetrance
    # Bin-2: 50-75% penetrance
    # Bin-3: 75-100% penetrance

    pene = df.Penetrance.tolist()
    bins = df['Penetrance_bin'].tolist()
    true_bins = []
    for b in bins:
        if b == 0:
            true_bins.append(12.5)
        elif b == 1:
            true_bins.append(37.5)
        elif b == 2:
            true_bins.append(62.5)
        elif b == 3:
            true_bins.append(87.5)

    # Add WT penetrance values to penetrance bin-0 block
    pene_wt = []
    true_bins_wt = []
    for w in wt:
        pene_wt.append(w)
        true_bins_wt.append(12.5)

    plt.figure(figsize=(6, 6))
    sns.set(font_scale=1.5)
    sns.set_style()
    plt.scatter(pene, true_bins, marker='.', c='black', alpha=1, s=80)
    plt.scatter(pene_wt, true_bins_wt, marker='.', c='green', alpha=1, s=80)
    plt.xlim([-2, 102])
    plt.ylim([-2, 102])
    plt.xticks([0, 25, 50, 75, 100])
    plt.yticks([0, 25, 50, 75, 100])
    plt.xlabel('Calculated Penetrance')
    plt.ylabel('Penetrance Bin')
    plt.title('Penetrance Agreement')
    fig = plt.gcf()
    plt.savefig(output)
    fig.clf()
    plt.close(fig)


def get_pos_control(df, pos_controls_df, bin, identifier, output, identifier_features):
    """ Read positive control file and append penetrance, p-value, number of cells, and predicted penetrance bin info

            Args:
                df:                     Existing combined dictionary
                pos_controls_file:      Positive control file
                bin:                    True if there is Penetrance_bin info in pos_controls_file
                identifier:             Unique strain identifier
                output:                 Output filenames
                identifier_features:    List of strain identifiers

            Return:
                pos_controls_df:        Dataframe that contains positive controls and penetrance
            """

    # Read and initialize dataframe
    pos_controls_df['Penetrance'] = np.zeros(len(pos_controls_df))
    pos_controls_df['P-value'] = np.zeros(len(pos_controls_df))
    pos_controls_df['Num_cells'] = np.zeros(len(pos_controls_df))
    pos_controls_df['Predicted_Penetrance_bin'] = np.zeros(len(pos_controls_df))

    for i in range(len(pos_controls_df)):
        # Search for the identifier entry from the penetrance output file
        g = pos_controls_df.iloc[i, pos_controls_df.columns.get_loc(identifier)]
        index = df[df[identifier] == g].index.tolist()
        if g in df[identifier].tolist():
            # Place penetrance, p-value, number of cells and predicted penetrance bins
            p = df.iloc[index[0], df.columns.get_loc('Penetrance')]
            pos_controls_df.iloc[i, pos_controls_df.columns.get_loc('Penetrance')] = p
            pval = df.iloc[index[0], df.columns.get_loc('P-value')]
            pos_controls_df.iloc[i, pos_controls_df.columns.get_loc('P-value')] = pval
            cellnum = df.iloc[index[0], df.columns.get_loc('Num_cells')]
            pos_controls_df.iloc[i, pos_controls_df.columns.get_loc('Num_cells')] = cellnum

            # Turn penetrance values to bins
            # Bin-0: 0-25% penetrance
            # Bin-1: 25-50% penetrance
            # Bin-2: 50-75% penetrance
            # Bin-3: 75-100% penetrance
            if p < 25:
                pos_controls_df.iloc[i, pos_controls_df.columns.get_loc('Predicted_Penetrance_bin')] = 0
            elif (p >= 25) and (p < 50):
                pos_controls_df.iloc[i, pos_controls_df.columns.get_loc('Predicted_Penetrance_bin')] = 1
            elif (p >= 50) and (p < 75):
                pos_controls_df.iloc[i, pos_controls_df.columns.get_loc('Predicted_Penetrance_bin')] = 2
            else:
                pos_controls_df.iloc[i, pos_controls_df.columns.get_loc('Predicted_Penetrance_bin')] = 3

        else:
            # If the identifier entry is not present
            pos_controls_df.iloc[i, pos_controls_df.columns.get_loc('Penetrance')] = np.nan
            pos_controls_df.iloc[i, pos_controls_df.columns.get_loc('P-value')] = np.nan
            pos_controls_df.iloc[i, pos_controls_df.columns.get_loc('Num_cells')] = np.nan
            pos_controls_df.iloc[i, pos_controls_df.columns.get_loc('Predicted_Penetrance_bin')] = np.nan

    # Save the penetrance and predicted bins
    pos_controls_df = pos_controls_df.reset_index(drop=True)
    if bin:
        pos_controls_df = pos_controls_df[identifier_features + ['Penetrance_bin', 'Predicted_Penetrance_bin',
                                                                 'Penetrance', 'P-value', 'Num_cells']]
    else:
        pos_controls_df = pos_controls_df[identifier_features + ['Predicted_Penetrance_bin',
                                                                 'Penetrance', 'P-value', 'Num_cells']]
    pos_controls_df.to_csv(path_or_buf=output, index=False)

    return pos_controls_df


def plot_performance(pos_control_file, df, wt, identifier, identifier_features, output):
    """ Plot ROC and PR curves, penetrance agreement and confusion matrices if the positive control file is available

            Args:
                pos_control:            Positive control file
                df:                     Existing combined dictionary
                wt:                     WT ORF names
                identifier:             Unique identifier
                identifier_features:    List of strain identifiers
                output:                 Output filenames
            """

    if pos_control_file != '':
        print('Plotting performance results...')
        log_write(output['log'], 'Plotting performance results...\n')

        output_ROC_PR = [output['ROCCurve'], output['PRCurve'], output['CurveNumbers']]
        pos_control = pd.read_csv(pos_control_file)

        if 'Penetrance_bin' not in pos_control.columns.values:
            # If the penetrance bin information is not available, plot only ROC and PR curves

            # Read positive controls file and add penetrance
            pos_controls_df = get_pos_control(df, pos_control, False, identifier,
                                              output['PenetranceBins'], identifier_features)

            # Remove genes that are not screened
            pos_controls_df = pos_controls_df.dropna(axis=0,
                                                     subset=['Penetrance', 'P-value', 'Predicted_Penetrance_bin'])
            pos_controls_df = pos_controls_df.reset_index(drop=True)

            # Positive and negative controls to plot ROC and PR curves
            PC = np.array(pos_controls_df.Penetrance)
            NC = np.array([])
            for w in wt:
                NC = np.append(NC, np.asarray(df[df.ORF == w].Penetrance))
            plot_ROC_and_PR(PC, NC, output_ROC_PR)

        else:
            # If the penetrance bin information is available, plot all

            # Read positive controls file and add penetrance
            pos_controls_df = get_pos_control(df, pos_control, True, identifier,
                                              output['PenetranceBins'], identifier_features)

            # Remove genes that are not screened
            pos_controls_df = pos_controls_df.dropna(axis=0,
                                                     subset=['Penetrance', 'P-value', 'Predicted_Penetrance_bin'])
            pos_controls_df = pos_controls_df.reset_index(drop=True)

            # Positive and negative controls to plot ROC and PR curves
            PC = np.array(pos_controls_df[(pos_controls_df['Penetrance_bin'] == 2) |
                                          (pos_controls_df['Penetrance_bin'] == 3)].Penetrance)

            NC = np.array(pos_controls_df[pos_controls_df['Penetrance_bin'] == 0].Penetrance)
            for w in wt:
                NC = np.append(NC, np.asarray(df[df.ORF == w].Penetrance))
            plot_ROC_and_PR(PC, NC, output_ROC_PR)

            # Plot confusion matrix
            actual = pos_controls_df.Penetrance_bin.tolist()
            predicted = pos_controls_df.Predicted_Penetrance_bin.tolist()
            actual[:] = [3 - x for x in actual]
            predicted[:] = [3 - x for x in predicted]
            classes = np.asarray(['75-100', '50-75', '25-50', '0-25'])
            confusion_matrix(actual, predicted, classes, output['ConfusionMatrix'])

            # Plot penetrance agreement
            plot_penetrance_agreement(pos_controls_df, NC, output['PenetranceAgreement'])
