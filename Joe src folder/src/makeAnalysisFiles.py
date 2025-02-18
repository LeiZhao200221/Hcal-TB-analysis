import pandas as pd
import numpy as np
import uproot
import os
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.pyplot import get_cmap
import mplhep as hep

hep.style.use(hep.style.ATLAS)

from apply_calibrations import *


# Main class to make analysis files (output is csv format)
class makeAnalysisFiles:
    def __init__(self, root_file_name, out_directory='../analysis_files/',
                 calibration_file='../calibrations/toa_calibration_phase3.csv', do_one_bar=False, do_alignment=True,
                 alignment_threshold=20, output_pulse_shapes=False):
        '''
        Initalization
        @param root_file_name: str or list of str pointing to input ROOT files
        @param out_directory: output directory for analysis csv files
        @param calibration_file: location of TOA calibration csv
        @param do_one_bar: debug flag, performs chain on only one bar
        @param do_alignment: set True if alignment of FPGAs should be performed- if True, then root_file_name should be a list of strs
        @param alignment_threshold: tolerance when considering how different pf_ticks can be between the two FPGAs to consider that event aligned
        @param output_pulse_shapes: save all 8 ADC samples to the output csv
        '''

        # Checks to see if the type of object passed as root_file_name is compatible with a request for alignment
        if (type(root_file_name) is not str and type(root_file_name) is not list):
            raise ValueError(
                'Input file format should be a string of a single file name, or a list of strings of multiple file names!')
        if (do_alignment is True and type(root_file_name) is not list):
            raise ValueError('Alignment procedure was requested, but only one file was input!')

        # Set variables
        if (type(root_file_name) == str):
            self.root_file_name = [root_file_name]
        else:
            self.root_file_name = root_file_name
        try:
            self.run_number = self.root_file_name[0].split('/')[-1].split('_')[5]
        except:
            self.run_number = self.root_file_name[0].split('_')[5]
        self.fpgas = []
        for i in range(len(self.root_file_name)):
            self.root_file_name[i] = self.root_file_name[i] + ':ntuplizehgcroc/hgcroc'
            try:
                self.fpgas.append(self.root_file_name[i].split('/')[-1].split('_')[3])
            except:
                self.fpgas.append(self.root_file_name[i].split('_')[3])
        self.out_directory = out_directory
        self.do_one_bar = do_one_bar
        self.toa_cal_file = pd.read_csv(calibration_file)
        self.do_alignment = do_alignment
        self.alignment_threshold = alignment_threshold
        self.output_pulse_shapes = output_pulse_shapes

        # Create the directory if it doesn't exist
        if not os.path.exists(self.out_directory):
            os.makedirs(self.out_directory)

    # Find events that have consistent pf_ticks in the two FPGAs
    def __align_ticks(self, data, layer, align_file_nbr):
        print('Performing alignment')

        data = data["pf_event", "pf_spill", "pf_ticks"]
        # Open file
        with uproot.open(self.root_file_name[align_file_nbr]) as in_file:
            cut = "layer == " + str(layer)
            in_data = pd.DataFrame()
            for batch in in_file.iterate(["pf_event", "pf_spill", "pf_ticks"], cut, library="pd", stepsize="10 MB"):
                in_data = pd.concat([in_data, batch])

        # Find unique values
        data = data.drop_duplicates().reset_index()
        in_data = in_data.drop_duplicates().reset_index()

        spills_data1 = data['pf_spill'].drop_duplicates().to_list()
        spills_data2 = in_data['pf_spill'].drop_duplicates().to_list()

        spills = list(set(spills_data1).intersection(spills_data2))

        zero_pf_events = []
        one_pf_events = []

        # Loop through pf_spills
        for spill in spills:
            print('Doing spill: ', spill)

            data1_ = data[data['pf_spill'] == spill]
            data2_ = in_data[in_data['pf_spill'] == spill]

            num_rows_1 = len(data1_)
            num_rows_2 = len(data2_)

            # Initialize indices
            index_df1 = 1
            index_df2 = 1

            # Iterate through events in each FPGA, and find pf_events that have consistent pf_ticks
            while index_df1 < num_rows_1 and index_df2 < num_rows_2:
                row_df1 = data1_.iloc[index_df1 - 1]
                row_df2 = data2_.iloc[index_df2 - 1]

                ticks_difference = abs(row_df1['pf_ticks'] - row_df2['pf_ticks'])

                # If the events are consistent, save the pf_event value
                if ticks_difference <= self.alignment_threshold:
                    zero_pf_events.append(row_df1['pf_event'])
                    one_pf_events.append(row_df2['pf_event'])

                    index_df1 += 1
                    index_df2 += 1

                # If not, see if the next event if consistent
                else:
                    temp_row_df1 = data1_.iloc[index_df1]
                    temp_row_df2 = data2_.iloc[index_df2]
                    if ((abs(temp_row_df1['pf_ticks'] - row_df2['pf_ticks'])) <= self.alignment_threshold):
                        index_df1 += 1
                    if ((abs(row_df1['pf_ticks'] - temp_row_df2['pf_ticks'])) <= self.alignment_threshold):
                        index_df2 += 1
                    else:
                        index_df1 += 1
                        index_df2 += 1

        return zero_pf_events, one_pf_events

    # Rework this FPGA file such that only aligned events are kept and reindex pf_events such that they are consistent between the FPGAs
    def __process_half(self, data, pf_events, half):
        print('Processing half: ', half)
        """ with uproot.open(data_file) as in_file:
            filtered_df = pd.DataFrame()
            for batch in in_file.iterate(
                    ["layer","end","strip","raw_id","adc","tot","toa","pf_event","pf_spill","pf_ticks"],
                    library="pd", step_size="10 MB"):

                filtered_df = pd.concat(filtered_df, batch[batch['pf_event'].isin(pf_events)])"""
        filtered_df = data[data['pf_event'].isin(pf_events)]
        pf_events_ = {}

        for i in range(len(pf_events)):
            pf_events_[pf_events[i]] = i

        filtered_df['pf_event'] = filtered_df['pf_event'].map(pf_events_)

        return filtered_df

    # Perform aggregation, where we extract one value for TOT, calibrated TOA and sum, max, mean of all 8 ADC values per event
    def __get_each_end(self, data_frame, end, toa_list):

        if self.output_pulse_shapes:
            aggregated_end = data_frame[data_frame['end'] == end].groupby('pf_event').agg({
                'layer': 'first',
                'strip': 'first',
                'pf_spill': 'first',
                'pf_ticks': 'first',
                'tot': tot_calib,
                'toa': lambda x: toa_calib(x, toa_list, end),
                'adc': ['sum', 'mean', 'max',
                        lambda x: x.iloc[0], lambda x: x.iloc[1], lambda x: x.iloc[2], lambda x: x.iloc[3],
                        lambda x: x.iloc[4], lambda x: x.iloc[5], lambda x: x.iloc[6], lambda x: x.iloc[7]]
            }).reset_index()
            aggregated_end.rename(
                columns={'adc_<lambda_0>': 'adc_0', 'adc_<lambda_1>': 'adc_1', 'adc_<lambda_2>': 'adc_2',
                         'adc_<lambda_3>': 'adc_3',
                         'adc_<lambda_4>': 'adc_4', 'adc_<lambda_5>': 'adc_5', 'adc_<lambda_6>': 'adc_6',
                         'adc_<lambda_7>': 'adc_7'})

            # aggregated_end.columns = ['pf_event','layer','strip','pf_spill','pf_ticks','tot','toa','adc_sum','adc_mean','adc_max', 'adc_first', 'adc_<lambda_0>', 'adc_<lambda_1>']
            aggregated_end.columns = ['pf_event', 'layer', 'strip', 'pf_spill', 'pf_ticks', 'tot', 'toa', 'adc_sum',
                                      'adc_mean', 'adc_max',
                                      'adc_0', 'adc_1', 'adc_2', 'adc_3', 'adc_4', 'adc_5', 'adc_6', 'adc_7']

            return aggregated_end
        else:
            aggregated_end = data_frame[data_frame['end'] == end].groupby('pf_event').agg({
                'layer': 'first',
                'strip': 'first',
                'pf_spill': 'first',
                'pf_ticks': 'first',
                'tot': tot_calib,
                'toa': lambda x: toa_calib(x, toa_list, end),
                'adc': ['sum', 'mean', 'max']
            }).reset_index()

            aggregated_end.columns = ['pf_event', 'layer', 'strip', 'pf_spill', 'pf_ticks', 'tot', 'toa', 'adc_sum',
                                      'adc_mean', 'adc_max']

            return aggregated_end

    # Clean per-end DataFrames for useable formats
    def __clean_frame(self, aggregated_end0, aggregated_end1):
        df_ = pd.merge(aggregated_end0, aggregated_end1, on='pf_event', suffixes=('_end0', '_end1'))

        df_ = df_.drop('layer_end0', axis=1)
        df_ = df_.drop('strip_end0', axis=1)
        df_ = df_.drop('pf_spill_end0', axis=1)
        df_ = df_.drop('pf_ticks_end0', axis=1)

        df_.rename(columns={'layer_end1': 'layer'}, inplace=True)
        df_.rename(columns={'strip_end1': 'strip'}, inplace=True)
        df_.rename(columns={'pf_spill_end1': 'pf_spill'}, inplace=True)
        df_.rename(columns={'pf_ticks_end1': 'pf_ticks'}, inplace=True)

        return df_

    # Process each layer and bar independently
    def __process_group(self, group):
        layer, bar = group.name

        print('layer: ', layer, ', bar: ', bar)
        toa_cal = (self.toa_cal_file['layer'] == layer - 1) & (self.toa_cal_file['bar'] == bar)
        toa_list = self.toa_cal_file[toa_cal].values.tolist()

        try:
            aggregated_end0 = self.__get_each_end(group, 0, toa_list)
            aggregated_end1 = self.__get_each_end(group, 1, toa_list)

            result_df = self.__clean_frame(aggregated_end0, aggregated_end1)
            return result_df
        except:
            print('empty bar')
            return None

    # Main function
    def create_dataframes(self):

        # If we are not doing alignment
        if (self.do_alignment is False):

            # Loop through provided ROOT files
            result_df = pd.DataFrame()
            for j in range(len(self.root_file_name)):
                # number of layers in the LDMX HCAL prototype
                layers = 19
                for i in range(layers):
                    in_data = pd.DataFrame()
                    with uproot.open(self.root_file_name[j]) as in_file:
                        print('reading file for layer ', i + 1, '...')
                        cut = "layer == " + str(i + 1)
                        batchnbr = 1
                        for batch in in_file.iterate(
                                ["layer", "end", "strip", "raw_id", "adc", "tot", "toa", "pf_event", "pf_spill",
                                 "pf_ticks"],
                                cut, library="pd", step_size="10 MB"):
                            print("batch: ", batchnbr)
                            in_data = pd.concat([in_data, batch])
                            batchnbr += 1

                    # If we only want to look at one bar, define this here
                    if self.do_one_bar is True and self.fpgas[j] == 0:
                        in_data = in_data[(in_data['layer'] == 1) & (in_data['strip'] == 3)]

                    if self.do_one_bar is True and self.fpgas[j] != 0:
                        continue

                    # Process each layer and bar independently
                    grouped_data = in_data.groupby(['layer', 'strip'], group_keys=False)

                    del in_data

                    loop_result_df = grouped_data.apply(self.__process_group)

                    del grouped_data

                    result_df = pd.concat([result_df, loop_result_df])

                    del loop_result_df

                # Save to csv file
                result_df.to_csv(self.out_directory + '/run_' + self.run_number + '_fpga_' + self.fpgas[j] + '.csv',
                                 index=False)

        # If we are doing alignment
        if (self.do_alignment is True):
            if self.do_one_bar is True:
                raise ValueError(
                    'Alignment is requested and the debugging feature is also requested. Does not make sense to align '
                    'but only look at one bar!')

            # read data
            layers = 19
            result_df0 = pd.DataFrame()
            for i in range(layers):
                in_data = pd.DataFrame()
                with uproot.open(self.root_file_name[0]) as in_file:
                    print('reading file for layer ', i + 1, '...')
                    cut = "layer == " + str(i + 1)
                    batchnbr = 1
                    for batch in in_file.iterate(
                            ["layer", "end", "strip", "raw_id", "adc", "tot", "toa", "pf_event", "pf_spill",
                             "pf_ticks"],
                            cut, library="pd", step_size="10 MB"):
                        print("batch: ", batchnbr)
                        in_data = pd.concat([in_data, batch])
                        batchnbr += 1

                # Align
                zero_pf_events, one_pf_events = self.__align_ticks(in_data, i, 1)

                # Process FPGA 0
                aligned_data0 = self.__process_half(in_data, zero_pf_events, 0)

                del in_data, zero_pf_events, one_pf_events

                # Process each layer and bar independently
                grouped_data0 = aligned_data0.groupby(['layer', 'strip'], group_keys=False)

                del aligned_data0

                loop_result_df0 = grouped_data0.apply(self.__process_group)

                del grouped_data0

                result_df0 = pd.concat([result_df0, loop_result_df0])

                del loop_result_df0

            # Save to csv file
            result_df0.to_csv(self.out_directory + '/run_' + self.run_number + '_fpga_' + self.fpgas[0] + '.csv',
                              index=False)

            # TODO: note this should probably bee a for loop over the two files to be more readable, but I'm lazy.
            # read data
            layers = 19
            result_df1 = pd.DataFrame()
            for i in range(layers):
                in_data = pd.DataFrame()
                with uproot.open(self.root_file_name[1]) as in_file:
                    print('reading file for layer ', i + 1, '...')
                    cut = "layer == " + str(i + 1)
                    batchnbr = 1
                    for batch in in_file.iterate(
                            ["layer", "end", "strip", "raw_id", "adc", "tot", "toa", "pf_event", "pf_spill",
                             "pf_ticks"],
                            cut, library="pd", step_size="10 MB"):
                        print("batch: ", batchnbr)
                        in_data = pd.concat([in_data, batch])
                        batchnbr += 1

                # Align
                zero_pf_events, one_pf_events = self.__align_ticks(in_data, i, 0)

                # Process FPGA 1
                aligned_data1 = self.__process_half(in_data, one_pf_events, 1)

                del in_data, zero_pf_events, one_pf_events

                # Process each layer and bar independently
                grouped_data1 = aligned_data1.groupby(['layer', 'strip'], group_keys=False)

                del aligned_data1

                loop_result_df1 = grouped_data1.apply(self.__process_group)

                del grouped_data1

                result_df1 = pd.concat([result_df1, loop_result_df1])

                del loop_result_df1

            # Save to csv file
            result_df1.to_csv(
                self.out_directory + '/run_' + self.run_number + '_fpga_' + self.fpgas[1] + '.csv', index=False)

            del result_df

    # Admittedly dirty function that makes many Simple_MIPs
    def make_plots(self, in_csv_file=None, plots_directory=None):
        '''
        Initalization
        @param in_csv_file (optional): csv containing events processed through create_dataframes(), if None, process files in memory
        @param plots_directory (optional): output directory for Simple_MIPs, if None, use default out_directory from makeAnalysisFiles class
        '''

        # Check Simple_MIPs directory and create if needed
        if plots_directory is None:
            plots_directory = self.out_directory
        else:
            if not os.path.exists(plots_directory):
                # Create the directory if it doesn't exist
                os.makedirs(plots_directory)

        # Open csv file
        in_files = None
        if in_csv_file is None:
            for fpga in range(len(self.fpgas)):
                in_files.append(self.out_directory + '/run_' + self.run_number + '_fpga_' + self.fpgas[fpga] + '.csv')
        elif type(in_csv_file) is str:
            in_files = [in_csv_file]
        elif type(in_csv_file) is list:
            in_files = in_csv_file
        else:
            raise ValueError('in_csv_file should be a list, str, or None type!')

        if (len(in_files) != len(self.fpgas)):
            raise ValueError('The length of input csv files must match the number of input files in the previous step!')

        # Loop through FPGAs
        for fpga in range(len(self.fpgas)):
            result_df = pd.read_csv(in_files[fpga])

            # If we only want to look at one bar, define this here
            if self.do_one_bar is True and self.fpgas[fpga] == 0:
                result_df = result_df[(result_df['layer'] == 1) & (result_df['strip'] == 3)]

            if self.do_one_bar is True and self.fpgas[fpga] != 0:
                continue

            # Get layers and bars
            layers = result_df['layer'].unique()
            bars = result_df['strip'].unique()

            # Loop through Simple_MIPs you want to make and plot for each layer/bar (I really don't know why I did it this way)
            for k in range(18):
                print(k)

                # Loop through layers
                for i in range(len(layers)):
                    if self.do_one_bar is True and i > 0:
                        continue

                    # Declare layer-level figures
                    fig1 = None
                    ax1 = None
                    if (k < 6):
                        fig1, ax1 = plt.subplots(figsize=(8, 8))

                    in_data_temp_ = result_df[result_df['layer'] == layers[i]]

                    # Loop through bars
                    for j in range(len(bars)):
                        if self.do_one_bar is True and j > 0:
                            continue

                        # Declare bar-level figures
                        selection = None
                        fig2 = None
                        ax2 = None
                        if (k > 5):
                            fig2, ax2 = plt.subplots(figsize=(8, 8))

                        # Select non-zero TOA events
                        if (k < 3):
                            selection = (in_data_temp_['strip'] == bars[j]) & (in_data_temp_['toa_end0'] > 0) & (
                                    in_data_temp_['toa_end1'] > 0)
                        elif (k < 12 and k > 5):
                            selection = (in_data_temp_['strip'] == bars[j]) & (in_data_temp_['toa_end0'] > 0) & (
                                    in_data_temp_['toa_end1'] > 0)

                        # Select exactly-zero TOA events
                        else:
                            selection = (in_data_temp_['strip'] == bars[j]) & (in_data_temp_['toa_end0'] == 0) & (
                                    in_data_temp_['toa_end1'] == 0)

                        # Impose selection
                        in_data_temp = in_data_temp_[selection]

                        # Non-zero TOA- sum ADC vs TOT (end 0)
                        if (k == 6):
                            ax2.hist2d(in_data_temp['adc_sum_end0'], in_data_temp['tot_end0'], bins=100, cmin=0.01,
                                       cmap='jet', norm=mcolors.LogNorm())
                            ax2.set_xlabel('ADC 0')
                            ax2.set_ylabel('TOT 0')
                            ax2.set_title('TOA Fired, Layer: ' + str(layers[i]) + ' Bar: ' + str(bars[j]) + ' Side: 0')
                            plt.savefig(
                                plots_directory + '/run_' + self.run_number + '_toa_fired_adc_tot_side_0_layer_' + str(
                                    layers[i]) + '_bar_' + str(bars[j]) + '.pdf')
                            plt.close()

                        # Non-zero TOA- sum ADC vs TOA (end 0)
                        if (k == 7):
                            ax2.hist2d(in_data_temp['adc_sum_end0'], in_data_temp['toa_end0'], bins=100, cmin=0.01,
                                       cmap='jet', norm=mcolors.LogNorm())
                            ax2.set_xlabel('ADC 0')
                            ax2.set_ylabel('TOA 0')
                            ax2.set_title('TOA Fired, Layer: ' + str(layers[i]) + ' Bar: ' + str(bars[j]) + ' Side: 0')
                            plt.savefig(
                                plots_directory + '/run_' + self.run_number + '_toa_fired_adc_toa_side_0_layer_' + str(
                                    layers[i]) + '_bar_' + str(bars[j]) + '.pdf')
                            plt.close()

                        # Non-zero TOA- TOT vs TOA (end 0)
                        if (k == 8):
                            ax2.hist2d(in_data_temp['tot_end0'], in_data_temp['toa_end0'], bins=100, cmin=0.01,
                                       cmap='jet', norm=mcolors.LogNorm())
                            ax2.set_xlabel('TOT 0')
                            ax2.set_ylabel('TOA 0')
                            ax2.set_title('TOA Fired, Layer: ' + str(layers[i]) + ' Bar: ' + str(bars[j]) + ' Side: 0')
                            plt.savefig(
                                plots_directory + '/run_' + self.run_number + '_toa_fired_tot_toa_side_0_layer_' + str(
                                    layers[i]) + '_bar_' + str(bars[j]) + '.pdf')
                            plt.close()

                        # Non-zero TOA- sum ADC (end 0) vs sum ADC (end 1)
                        if (k == 9):
                            ax2.hist2d(in_data_temp['adc_sum_end0'], in_data_temp['adc_sum_end1'], bins=100, cmin=0.01,
                                       cmap='jet', norm=mcolors.LogNorm())
                            ax2.set_xlabel('ADC 0')
                            ax2.set_ylabel('ADC 1')
                            ax2.set_title('TOA Fired, Layer: ' + str(layers[i]) + ' Bar: ' + str(bars[j]))
                            plt.savefig(plots_directory + '/run_' + self.run_number + '_toa_fired_adc_adc_layer_' + str(
                                layers[i]) + '_bar_' + str(bars[j]) + '.pdf')
                            plt.close()

                        # Non-zero TOA- TOA (end 0) vs TOA (end 1)
                        if (k == 10):
                            ax2.hist2d(in_data_temp['toa_end0'], in_data_temp['toa_end1'], bins=100, cmin=0.01,
                                       cmap='jet', norm=mcolors.LogNorm())
                            ax2.set_xlabel('TOA 0')
                            ax2.set_ylabel('TOA 1')
                            ax2.set_title('TOA Fired, Layer: ' + str(layers[i]) + ' Bar: ' + str(bars[j]))
                            plt.savefig(plots_directory + '/run_' + self.run_number + '_toa_fired_toa_toa_layer_' + str(
                                layers[i]) + '_bar_' + str(bars[j]) + '.pdf')
                            plt.close()

                        # Non-zero TOA- TOT (end 0) vs TOT (end 1)
                        if (k == 11):
                            ax2.hist2d(in_data_temp['tot_end0'], in_data_temp['tot_end1'], bins=100, cmin=0.01,
                                       cmap='jet', norm=mcolors.LogNorm())
                            ax2.set_xlabel('TOT 0')
                            ax2.set_ylabel('TOT 1')
                            ax2.set_title('TOA Fired, Layer: ' + str(layers[i]) + ' Bar: ' + str(bars[j]))
                            plt.savefig(plots_directory + '/run_' + self.run_number + '_toa_fired_tot_tot_layer_' + str(
                                layers[i]) + '_bar_' + str(bars[j]) + '.pdf')
                            plt.close()

                        # Exactly-zero TOA- sum ADC vs TOT (end 0)
                        if (k == 12):
                            ax2.hist2d(in_data_temp['adc_sum_end0'], in_data_temp['tot_end0'], bins=100, cmin=0.01,
                                       cmap='jet', norm=mcolors.LogNorm())
                            ax2.set_xlabel('ADC 0')
                            ax2.set_ylabel('TOT 0')
                            ax2.set_title(
                                'TOA Not Fired, Layer: ' + str(layers[i]) + ' Bar: ' + str(bars[j]) + ' Side: 0')
                            plt.savefig(
                                plots_directory + '/run_' + self.run_number + '_toa_not_fired_adc_tot_side_0_layer_' + str(
                                    layers[i]) + '_bar_' + str(bars[j]) + '.pdf')
                            plt.close()

                        # Exactly-zero TOA- sum ADC vs TOA (end 0)
                        if (k == 13):
                            ax2.hist2d(in_data_temp['adc_sum_end0'], in_data_temp['toa_end0'], bins=100, cmin=0.01,
                                       cmap='jet', norm=mcolors.LogNorm())
                            ax2.set_xlabel('ADC 0')
                            ax2.set_ylabel('TOA 0')
                            ax2.set_title(
                                'TOA Not Fired, Layer: ' + str(layers[i]) + ' Bar: ' + str(bars[j]) + ' Side: 0')
                            plt.savefig(
                                plots_directory + '/run_' + self.run_number + '_toa_not_fired_adc_toa_side_0_layer_' + str(
                                    layers[i]) + '_bar_' + str(bars[j]) + '.pdf')
                            plt.close()

                        # Exactly-zero TOA- TOT vs TOA (end 0)
                        if (k == 14):
                            ax2.hist2d(in_data_temp['tot_end0'], in_data_temp['toa_end0'], bins=100, cmin=0.01,
                                       cmap='jet', norm=mcolors.LogNorm())
                            ax2.set_xlabel('TOT 0')
                            ax2.set_ylabel('TOA 0')
                            ax2.set_title(
                                'TOA Not Fired, Layer: ' + str(layers[i]) + ' Bar: ' + str(bars[j]) + ' Side: 0')
                            plt.savefig(
                                plots_directory + '/run_' + self.run_number + '_toa_not_fired_tot_toa_side_0_layer_' + str(
                                    layers[i]) + '_bar_' + str(bars[j]) + '.pdf')
                            plt.close()

                        # Exactly-zero TOA- sum ADC (end 0) vs sum ADC (end 1)
                        if (k == 15):
                            ax2.hist2d(in_data_temp['adc_sum_end0'], in_data_temp['adc_sum_end1'], bins=100, cmin=0.01,
                                       cmap='jet', norm=mcolors.LogNorm())
                            ax2.set_xlabel('ADC 0')
                            ax2.set_ylabel('ADC 1')
                            ax2.set_title('TOA Not Fired, Layer: ' + str(layers[i]) + ' Bar: ' + str(bars[j]))
                            plt.savefig(
                                plots_directory + '/run_' + self.run_number + '_toa_not_fired_adc_adc_layer_' + str(
                                    layers[i]) + '_bar_' + str(bars[j]) + '.pdf')
                            plt.close()

                        # Exactly-zero TOA- TOA (end 0) vs TOA (end 1)
                        if (k == 16):
                            ax2.hist2d(in_data_temp['toa_end0'], in_data_temp['toa_end1'], bins=100, cmin=0.01,
                                       cmap='jet', norm=mcolors.LogNorm())
                            ax2.set_xlabel('TOA 0')
                            ax2.set_ylabel('TOA 1')
                            ax2.set_title('TOA Not Fired, Layer: ' + str(layers[i]) + ' Bar: ' + str(bars[j]))
                            plt.savefig(
                                plots_directory + '/run_' + self.run_number + '_toa_not_fired_toa_toa_layer_' + str(
                                    layers[i]) + '_bar_' + str(bars[j]) + '.pdf')
                            plt.close()

                        # Exactly-zero TOA- TOT (end 0) vs TOT (end 1)
                        if (k == 17):
                            ax2.hist2d(in_data_temp['tot_end0'], in_data_temp['tot_end1'], bins=100, cmin=0.01,
                                       cmap='jet', norm=mcolors.LogNorm())
                            ax2.set_xlabel('TOT 0')
                            ax2.set_ylabel('TOT 1')
                            ax2.set_title('TOA Not Fired, Layer: ' + str(layers[i]) + ' Bar: ' + str(bars[j]))
                            plt.savefig(
                                plots_directory + '/run_' + self.run_number + '_toa_not_fired_tot_tot_layer_' + str(
                                    layers[i]) + '_bar_' + str(bars[j]) + '.pdf')
                            plt.close()

                        # Declare a histogram for each layer-level Simple_MIPs (one for each bar)- only considering end 0 since ends are correlated
                        if (k == 0):
                            ax1.hist(in_data_temp['adc_sum_end0'], bins=100, histtype='step',
                                     label='Bar ' + str(bars[j]))
                        if (k == 1):
                            ax1.hist(in_data_temp['toa_end0'], bins=100, histtype='step', label='Bar ' + str(bars[j]))
                        if (k == 2):
                            ax1.hist(in_data_temp['tot_end0'], bins=100, histtype='step', label='Bar ' + str(bars[j]))
                        if (k == 3):
                            ax1.hist(in_data_temp['adc_sum_end0'], bins=100, histtype='step',
                                     label='Bar ' + str(bars[j]))
                        if (k == 4):
                            ax1.hist(in_data_temp['toa_end0'], bins=100, histtype='step', label='Bar ' + str(bars[j]))
                        if (k == 5):
                            ax1.hist(in_data_temp['tot_end0'], bins=100, histtype='step', label='Bar ' + str(bars[j]))

                    # Non-zero TOA- sum ADC all bars (end 0)
                    if (k == 0):
                        ax1.set_yscale('log')
                        ax1.set_xlabel('ADC 0')
                        ax1.set_ylabel('Events')
                        ax1.legend()
                        ax1.set_title('TOA Fired, Layer: ' + str(layers[i]))
                        plt.savefig(plots_directory + '/run_' + self.run_number + '_toa_fired_adc_layer_' + str(
                            layers[i]) + '_side_0.pdf')
                        plt.close()

                    # Non-zero TOA- TOA all bars (end 0)
                    if (k == 1):
                        ax1.set_yscale('log')
                        ax1.set_xlabel('TOA 0')
                        ax1.set_ylabel('Events')
                        ax1.legend()
                        ax1.set_title('TOA Fired, Layer: ' + str(layers[i]))
                        plt.savefig(plots_directory + '/run_' + self.run_number + '_toa_fired_toa_layer_' + str(
                            layers[i]) + '_side_0.pdf')
                        plt.close()

                    # Non-zero TOA- TOT all bars (end 0)
                    if (k == 2):
                        ax1.set_yscale('log')
                        ax1.set_xlabel('TOT 0')
                        ax1.set_ylabel('Events')
                        ax1.legend()
                        ax1.set_title('TOA Fired, Layer: ' + str(layers[i]))
                        plt.savefig(plots_directory + '/run_' + self.run_number + '_toa_fired_tot_layer_' + str(
                            layers[i]) + '_side_0.pdf')
                        plt.close()

                    # Exactly-zero TOA- sum ADC all bars (end 0)
                    if (k == 3):
                        ax1.set_yscale('log')
                        ax1.set_xlabel('ADC 0')
                        ax1.set_ylabel('Events')
                        ax1.legend()
                        ax1.set_title('TOA Not Fired, Layer: ' + str(layers[i]))
                        plt.savefig(plots_directory + '/run_' + self.run_number + '_toa_not_fired_adc_layer_' + str(
                            layers[i]) + '_side_0.pdf')
                        plt.close()

                    # Exactly-zero TOA- TOA all bars (end 0)
                    if (k == 4):
                        ax1.set_yscale('log')
                        ax1.set_xlabel('TOA 0')
                        ax1.set_ylabel('Events')
                        ax1.legend()
                        ax1.set_title('TOA Not Fired, Layer: ' + str(layers[i]))
                        plt.savefig(plots_directory + '/run_' + self.run_number + '_toa_not_fired_toa_layer_' + str(
                            layers[i]) + '_side_0.pdf')
                        plt.close()

                    # Exactly-zero TOA- TOT all bars (end 0)
                    if (k == 5):
                        ax1.set_yscale('log')
                        ax1.set_xlabel('TOT 0')
                        ax1.set_ylabel('Events')
                        ax1.legend()
                        ax1.set_title('TOA Not Fired, Layer: ' + str(layers[i]))
                        plt.savefig(plots_directory + '/run_' + self.run_number + '_toa_not_fired_tot_layer_' + str(
                            layers[i]) + '_side_0.pdf')
                        plt.close()
