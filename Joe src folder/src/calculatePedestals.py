import pandas as pd
import numpy as np
import uproot
import os
import scipy.stats as stats
import scipy.signal as signal
from statistics import mean, stdev, median
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.pyplot import get_cmap
import mplhep as hep
import gc

hep.style.use(hep.style.ATLAS)


# Main class to calculate and save pedestals
class calculatePedestals:
    def __init__(self, root_file_name, out_directory='../calibrations/', plot_pedestals=True,
                 plots_directory='../plots/pedestals', do_one_bar=False, layer_wise=False, in_batches=False, time_trend=False):
        '''
        Initalization
        @param root_file_name: str or list of str pointing to input ROOT files
        @param out_directory: output directory for pedestal calibration csv files
        @param plot_pedestals: flag for plotting pedestal results- if True, then plot
        @param plots_directory: output directory for pedestal plots
        @param do_one_bar: debug flag, performs chain on only one bar
        @param layer_wise: flag to save on memory usage, process one layer at a time
        @param in_batches: flag to save on memory usage, load in data from in file in batches and process one layer at a time
        @param time_trend: DO NOT USE! DOES NOT WORK! flag for plotting the adc values over time.
        @param noise_threshold: defines the distance from the pedestal from which data will be cut when doing a gaussian fit
        '''

        # Checks to see if the type of object passed as root_file_name is compatible with a request for alignment
        if (type(root_file_name) is not str and type(root_file_name) is not list):
            raise ValueError(
                'Input file format should be a string of a single file name, or a list of strings of multiple file names!')

        # Set variables
        if (type(root_file_name) == str):
            self.root_file_name = [root_file_name]
        else:
            self.root_file_name = root_file_name
        try:
            self.run_number = self.root_file_name[0].split('/')[-1].split('_')[5]
        except:
            self.run_number = self.root_file_name[0].split('_')[5]
        # if(self.run_number != '287'):
        # raise ValueError('Expecting run 287 for 4 GeV defocused muons!!')
        self.fpgas = []
        for i in range(len(self.root_file_name)):
            self.root_file_name[i] = self.root_file_name[i] + ':ntuplizehgcroc/hgcroc'
            try:
                self.fpgas.append(self.root_file_name[i].split('/')[-1].split('_')[3])
            except:
                self.fpgas.append(self.root_file_name[i].split('_')[3])
        self.out_directory = out_directory
        self.do_one_bar = do_one_bar
        self.plot_pedestals = plot_pedestals
        self.time_trend = time_trend
        self.layer_wise = layer_wise
        self.in_batches = in_batches
        if (self.in_batches):
            self.layer_wise = True
        self.out_ped_individ = {}
        self.out_ped_sum = {'layer': [],
                            'strip': [],
                            'end': [],
                            'pedestal': [],
                            'mean': [],
                            'std_dev': [],
                            'pedestal_per_time_sample': [],
                            'pedestal_per_time_sample_mean': [],
                            'pedestal_per_time_sample_std_dev': []}
        self.adc_sum_dict = {}

        # Create the directory if it doesn't exist
        if not os.path.exists(self.out_directory):
            os.makedirs(self.out_directory)

        # Check plots directory and create if needed
        if plot_pedestals or time_trend:
            if plots_directory is None:
                self.plots_directory = self.out_directory
            else:
                self.plots_directory = plots_directory
                if not os.path.exists(self.plots_directory):
                    # Create the directory if it doesn't exist
                    os.makedirs(self.plots_directory)

    # Properly index events in DataFrame manipulation step
    def __pivot_dataframe(self, dataframe, end):
        # Do some trickery to properly index the events
        multiindex_df = dataframe[dataframe['end'] == end].set_index(
            ['pf_event', dataframe[dataframe['end'] == end].groupby('pf_event').cumcount()])

        # Pivot the DataFrame to make 'adc' values into separate columns for each end
        pivoted_df = multiindex_df['adc'].unstack().add_prefix('adc_end_' + str(end)).reset_index()

        return pivoted_df

    # Properly format DataFrame in DataFrame manipulation step
    def __format_dataframe(self, dataframe, end):
        dataframe.to_csv(self.out_directory + 'non_agg.csv')
        # Perform aggregation, where we extract one value for TOT, TOA and all 8 values of ADC per event
        aggregated_end = dataframe[dataframe['end'] == end].groupby('pf_event').agg({
            'layer': 'first',
            'strip': 'first',
            'tot': 'sum',
            'toa': 'sum',
            'adc': 'sum'
        }).reset_index()

        # Do some manipulation to merge both ends of the bar
        aggregated_end.columns = ['pf_event', 'layer', 'strip', 'tot', 'toa', 'adc_sum']

        return aggregated_end

    # Clean DataFrame for easy handling
    def __clean_dataframes(self, aggregated_end0, aggregated_end1):
        df_ = pd.merge(aggregated_end0, aggregated_end1, on='pf_event', suffixes=('_end0', '_end1'))

        df_ = df_.drop('layer_end0', axis=1)
        df_ = df_.drop('strip_end0', axis=1)

        df_.rename(columns={'layer_end1': 'layer'}, inplace=True)
        df_.rename(columns={'strip_end1': 'strip'}, inplace=True)

        return df_

    # Calculate per-time sample pedestals (to be subtracted off each individual time sample ADC)
    def __get_individual_pedestals(self, group):
        layer, bar = group.name

        index0 = group.columns.str.contains(r'adc_end_0')
        end_0 = group.iloc[:, index0].to_numpy().flatten()
        end_0 = np.nan_to_num(end_0)

        index1 = group.columns.str.contains(r'adc_end_1')
        end_1 = group.iloc[:, index1].to_numpy().flatten()
        end_1 = np.nan_to_num(end_1)

        pedestal0 = stats.mode(end_0)[0]
        pedestal1 = stats.mode(end_1)[0]
        self.out_ped_individ['layer_' + str(layer) + '_bar_' + str(bar) + '_end_0'] = pedestal0
        self.out_ped_individ['layer_' + str(layer) + '_bar_' + str(bar) + '_end_1'] = pedestal1

        # Identify peak widths
        bins_0 = int(max(end_0) - min(end_0))
        bins_1 = int(max(end_1) - min(end_1))
        hist_0 = np.histogram(end_0, bins=bins_0)[0]
        hist_1 = np.histogram(end_1, bins=bins_1)[0]
        # Under the assumption that the location in hist of max(hist) = pedestal
        hist_0_max = max(hist_0)
        hist_1_max = max(hist_1)
        peak_width_0 = signal.peak_widths(hist_0, np.where(hist_0 == hist_0_max)[0])[0]
        peak_width_1 = signal.peak_widths(hist_1, np.where(hist_1 == hist_1_max)[0])[0]

        # cut away high and low noise
        end_0_cut = end_0[(end_0 <= pedestal0 + peak_width_0 * 3) & (end_0 >= pedestal0 - peak_width_0 * 3)]
        end_1_cut = end_1[(end_1 <= pedestal1 + peak_width_1 * 3) & (end_1 >= pedestal1 - peak_width_1 * 3)]

        # Fit a Gaussian to the pedestals
        mean0, std_dev0 = stats.norm.fit(end_0_cut)
        mean1, std_dev1 = stats.norm.fit(end_1_cut)

        self.out_ped_individ['layer_' + str(layer) + '_bar_' + str(bar) + '_mean_end_0'] = mean0
        self.out_ped_individ['layer_' + str(layer) + '_bar_' + str(bar) + '_std_dev_end_0'] = std_dev0
        self.out_ped_individ['layer_' + str(layer) + '_bar_' + str(bar) + '_mean_end_1'] = mean1
        self.out_ped_individ['layer_' + str(layer) + '_bar_' + str(bar) + '_std_dev_end_1'] = std_dev1
        if self.plot_pedestals:
            self.__plot_pedestal(end_0, layer, bar, 0, mean0, std_dev0, sum_pedestal=False, log_scale=True)
            self.__plot_pedestal(end_1, layer, bar, 1, mean1, std_dev1, sum_pedestal=False, log_scale=True)

    def __plot_pedestal(self, dataframe, layer, bar, end, mean, std_dev, sum_pedestal=True, log_scale=True):
        print("creating plots: layer: ", layer, "bar: ", bar)
        # Define figure for end and plot
        fig = plt.figure(num=1, clear=True)
        ax = fig.add_subplot()
        # fig,ax = plt.subplots(figsize=(8, 8))
        if (sum_pedestal):
            ax.hist(dataframe, bins=600, density=True, histtype='step', range=[-200, 400])
        else:
            ax.hist(dataframe, bins=50, density=True, histtype='step', range=[0, 200])

        xmin, xmax = ax.get_xlim()
        x = np.linspace(xmin, xmax, 1000)
        p = stats.norm.pdf(x, mean, std_dev)

        ax.plot(x, p, 'k', linewidth=1, linestyle='--')

        if (log_scale):
            ax.set_yscale('log')
        if (sum_pedestal):
            ax.set_xlabel('Sum ADC ')
        else:
            ax.set_xlabel('mean ADC')
        ax.set_ylabel('Events (Normalized to density)')
        if (log_scale):
            ax.set_ylim(1e-4, 1)
        ax.legend()

        if (sum_pedestal):
            ax.set_title(
                'Sum pedestal layer: ' + str(layer) + ' Bar: ' + str(bar) + ' Side: ' + str(end) + ' mean: ' + str(
                    round(mean, 2)) + ' std_dev: ' + str(round(std_dev, 2)))
        else:
            ax.set_title('Pedestal of adc mean')
            # ax.set_title('Pedestal per time step layer: ' + str(layer) + ' Bar: ' + str(bar) + ' Side: ' + str(
                # end) + ' mean: ' + str(
                # round(mean, 2)) + ' std_dev: ' + str(round(std_dev, 2)))

        if (sum_pedestal):
            plt.savefig(
                self.plots_directory + '/sum_ped_subtracted_layer_' + str(layer) + '_bar_' + str(
                    bar) + '_side_' + str(end) + '.pdf')
        else:
            plt.savefig(self.plots_directory + '/ped_per_time_step_layer_' + str(layer) + '_bar_' + str(
                bar) + '_side_' + str(end) + '.pdf')
        fig.clear()
        plt.close(fig)

    # TODO: WORK IN PROGRESS
    def __plot_time_trend(self, group):
        layer, bar = group.name
        print("creating plots: layer: ", layer, "bar: ", bar)
        ax1 = group.plot.scatter(x='pf_ticks', y='adc_sum_end0')
        ax1.set_xlabel('pf_ticks')
        ax1.set_ylabel('adc')
        ax1.set_title('adc as function of pf_ticks in layer ' + str(layer) + ' bar ' + str(bar) + ' end ' + ' 0')
        plt.show()
        plt.savefig(self.plots_directory + '/time_trends' + '/time_trend_layer_' + str(layer) + '_bar_' + str(bar) + '_end_' + '0' + '.pdf')
        ax1 = group.plot.scatter(x='pf_ticks', y='adc_sum_end1')
        ax1.set_xlabel('pf_ticks')
        ax1.set_ylabel('adc')
        ax1.set_title('adc as function of pf_ticks in layer ' + str(layer) + ' bar ' + str(bar) + ' end ' + ' 1')
        plt.show()
        plt.savefig(self.plots_directory + '/time_trends' + '/time_trend_layer_' + str(layer) + '_bar_' + str(bar) + '_end_' + '1' + '.pdf')


    # Calculate sum of ADC pedestals (to be subtracted off the case of all 8 time samples added)
    def __get_sum_pedestals(self, group):
        layer, bar = group.name

        # Obtain pedestal appropriate for summation of all 8 time samples
        pedestal_temp0 = 8 * self.out_ped_individ['layer_' + str(layer) + '_bar_' + str(bar) + '_end_0']
        pedestal_temp1 = 8 * self.out_ped_individ['layer_' + str(layer) + '_bar_' + str(bar) + '_end_1']

        # Apply pedestal to sum of ADC events
        group['adc_sum_end0'] = group['adc_sum_end0'] - (pedestal_temp0)
        group['adc_sum_end1'] = group['adc_sum_end1'] - (pedestal_temp1)

        # I don't really understand the point of this? Why are we doing this?
        # Select sum pedestal region for plotting
        fit0_ = (group['adc_sum_end0'] < 200) & (group['adc_sum_end0'] > -200)
        fit1_ = (group['adc_sum_end1'] < 200) & (group['adc_sum_end1'] > -200)
        fit0 = group[fit0_]
        fit1 = group[fit1_]

        # Identify peak widths
        end_0 = fit0['adc_sum_end0'].to_numpy().flatten()
        end_1 = fit1['adc_sum_end1'].to_numpy().flatten()
        print('end_0 in sum: ', end_0)
        bins_0 = int(max(end_0) - min(end_0))
        bins_1 = int(max(end_1) - min(end_1))
        hist_0 = np.histogram(end_0, bins=bins_0)[0]
        print('hist_0 in sum: ', hist_0)
        hist_1 = np.histogram(end_1, bins=bins_1)[0]
        # Under the assumption that the location in hist of max(hist) + min(end) = pedestal
        hist_0_max_loc = np.where(hist_0 == max(hist_0))[0]
        hist_1_max_loc = np.where(hist_1 == max(hist_1))[0]
        peak_width_0 = signal.peak_widths(hist_0, hist_0_max_loc)[0]
        print('peak_width in sum: ', peak_width_0)
        peak_width_1 = signal.peak_widths(hist_1, hist_1_max_loc)[0]
        if len(hist_0_max_loc) > 1:
            end_0_cut = end_0[(end_0 <= min(hist_0_max_loc) + min(end_0) + peak_width_0[0] * 3) & (
                        end_0 >= max(hist_0_max_loc) + min(end_0) - peak_width_0[-1] * 3)]
        # cut away high and low noise
        else:
            end_0_cut = end_0[(end_0 <= hist_0_max_loc + min(end_0) + peak_width_0 * 3) & (end_0 >= hist_0_max_loc + min(end_0) - peak_width_0 * 3)]
        print('end_0_cut in sum: ', end_0_cut)
        if len(hist_1_max_loc) > 1:
            end_1_cut = end_1[(end_1 <= min(hist_1_max_loc) + min(end_1) + peak_width_1[0] * 3) & (
                        end_1 >= max(hist_1_max_loc) + min(end_1) - peak_width_1[-1] * 3)]
        else:
            end_1_cut = end_1[(end_1 <= hist_1_max_loc + min(end_1) + peak_width_1 * 3) & (end_1 >= hist_1_max_loc + min(end_1) - peak_width_1 * 3)]

        del end_0, end_1

        # Fit a Gaussian to the pedestals
        mean0, std_dev0 = stats.norm.fit(end_0_cut)
        print('mean0 and std_dev0 in sum: ', mean0, std_dev0)
        mean1, std_dev1 = stats.norm.fit(end_1_cut)

        del end_0_cut, end_1_cut

        self.out_ped_sum['layer'].append(layer)
        self.out_ped_sum['strip'].append(bar)
        self.out_ped_sum['end'].append(0)
        self.out_ped_sum['pedestal'].append(pedestal_temp0)
        self.out_ped_sum['mean'].append(mean0)
        self.out_ped_sum['std_dev'].append(std_dev0)
        self.out_ped_sum['pedestal_per_time_sample'].append(
            self.out_ped_individ['layer_' + str(layer) + '_bar_' + str(bar) + '_end_0'])
        self.out_ped_sum['pedestal_per_time_sample_mean'].append(
            self.out_ped_individ['layer_' + str(layer) + '_bar_' + str(bar) + '_mean_end_0'])
        self.out_ped_sum['pedestal_per_time_sample_std_dev'].append(
            self.out_ped_individ['layer_' + str(layer) + '_bar_' + str(bar) + '_std_dev_end_0'])

        self.out_ped_sum['layer'].append(layer)
        self.out_ped_sum['strip'].append(bar)
        self.out_ped_sum['end'].append(1)
        self.out_ped_sum['pedestal'].append(pedestal_temp1)
        self.out_ped_sum['mean'].append(mean1)
        self.out_ped_sum['std_dev'].append(std_dev1)
        self.out_ped_sum['pedestal_per_time_sample'].append(
            self.out_ped_individ['layer_' + str(layer) + '_bar_' + str(bar) + '_end_1'])
        self.out_ped_sum['pedestal_per_time_sample_mean'].append(
            self.out_ped_individ['layer_' + str(layer) + '_bar_' + str(bar) + '_mean_end_1'])
        self.out_ped_sum['pedestal_per_time_sample_std_dev'].append(
            self.out_ped_individ['layer_' + str(layer) + '_bar_' + str(bar) + '_std_dev_end_1'])

        print(self.out_ped_sum)

        # Make plots
        if self.plot_pedestals is True:
            self.__plot_pedestal(group['adc_sum_end0'], layer, bar, 0, mean0, std_dev0)
            self.__plot_pedestal(group['adc_sum_end1'], layer, bar, 1, mean1, std_dev1)


    def __process_group(self, group):
        layer, bar = group.name

        print('layer: ', layer, ', bar: ', bar)

        # Manipulation
        pivoted_df_0 = self.__pivot_dataframe(group, 0)
        pivoted_df_1 = self.__pivot_dataframe(group, 1)

        pivoted_df = pd.merge(pivoted_df_0, pivoted_df_1, on='pf_event').reset_index()

        del pivoted_df_0
        del pivoted_df_1
        gc.collect()

        aggregated_end0 = self.__format_dataframe(group, 0)
        aggregated_end1 = self.__format_dataframe(group, 1)

        result_df = self.__clean_dataframes(aggregated_end0, aggregated_end1)

        del aggregated_end0
        del aggregated_end1
        gc.collect()

        result_df = pd.merge(result_df, pivoted_df, on='pf_event').reset_index()

        result_df.to_csv(self.out_directory + 'result_df')

        return result_df

    def __pedestal_calc(self, in_data):

        # First, we have to create a temporary DataFrame that keeps all 8 time samples to calculate a pedestal to subtract off individual time samples
        # Process each layer and bar independently
        print('grouping data...')
        grouped_data = in_data.groupby(['layer', 'strip'], group_keys=False)
        print('finished grouping data')

        del in_data
        gc.collect()

        print('processing groups...')
        final_result_df = grouped_data.apply(self.__process_group)
        print('finished processing groups')

        # Select likely pedestal events where the TOA and TOT are both zero on each end of the bar
        # selection = (final_result_df['tot_end0'] == 0) & (final_result_df['tot_end1'] == 0) & (
        #         final_result_df['toa_end0'] == 0) & (final_result_df['toa_end1'] == 0)

        # final_result_df = final_result_df[selection]

        # Get individual pedestals
        grouped_data_individual = final_result_df.groupby(['layer', 'strip'], group_keys=False)

        print('getting individual pedestals')
        grouped_data_individual.apply(self.__get_individual_pedestals)

        # Now we will get the sum of ADC-appropriate pedestals
        grouped_data_sum = final_result_df.groupby(['layer', 'strip'], group_keys=False)

        grouped_data_sum.apply(self.__get_sum_pedestals)

        if(self.time_trend):
            grouped_data_sum.apply(self.__plot_time_trend)

    # Main function to calculate pedestals. Creates a csv that has both pedestals appropriate for individual time samples as well as for the sum of ADC
    def get_pedestals(self):

        # Loop through provided ROOT files
        for j in range(len(self.root_file_name)):

            if (self.layer_wise):
                layers = 19
                if (self.in_batches):
                    for i in range(layers):
                        in_data = pd.DataFrame()
                        with uproot.open(self.root_file_name[j]) as in_file:
                            print('reading file for layer ', i + 1, '...')
                            cut = "layer == " + str(i + 1)
                            batchnbr = 1
                            for batch in in_file.iterate(
                                    ["layer", "end", "strip", "raw_id", "adc", "tot", "toa", "pf_event", "pf_spill",
                                     "pf_ticks"], cut, library="pd", step_size="10 MB"):
                                print("batch: ", batchnbr)
                                in_data = pd.concat([in_data, batch])
                                batchnbr += 1

                            print('finished reading file')

                        self.__pedestal_calc(in_data)

                        del in_data
                        gc.collect()

                    # Create our final DataFrame with all pedestals
                    print(self.out_ped_sum)
                    ped_df = pd.DataFrame(self.out_ped_sum)

                    # Save our pedestals to a csv file
                    print('saving to file')
                    ped_df.to_csv(self.out_directory + '/pedestals_MIP.csv', index=False)



                else:
                    for i in range(layers):
                        with uproot.open(self.root_file_name[j]) as in_file:
                            print('reading file for layer ', i, '...')
                            cut = "layer == " + str(i + 1)
                            in_data = in_file.arrays(
                                ["layer", "end", "strip", "raw_id", "adc", "tot", "toa", "pf_event", "pf_spill",
                                 "pf_ticks"], cut,
                                library="pd")
                            print('finished reading file')

                        self.__pedestal_calc(in_data)

                        # Create our final DataFrame with all pedestals
                        ped_df = pd.DataFrame(self.out_ped_sum)

                        # Save our pedestals to a csv file
                        ped_df.to_csv(self.out_directory + '/pedestals_MIP' +'_run_' + str(self.run_number) + '.csv', index=False)

                        del in_data
                        gc.collect()
            else:

                with uproot.open(self.root_file_name[j]) as in_file:
                    print('reading file...')
                    in_data = in_file.arrays(
                        ["layer", "end", "strip", "raw_id", "adc", "tot", "toa", "pf_event", "pf_spill", "pf_ticks"],
                        library="pd")
                    print('finished reading file')

                # If we only want to look at one bar, define this here
                if self.do_one_bar is True and self.fpgas[j] == 0:
                    in_data = in_data[(in_data['layer'] == 1) & (in_data['strip'] == 3)]

                if self.do_one_bar is True and self.fpgas[j] != 0:
                    continue

                self.__pedestal_calc(in_data)

                # Create our final DataFrame with all pedestals
                ped_df = pd.DataFrame(self.out_ped_sum)

                # Save our pedestals to a csv file
                ped_df.to_csv(self.out_directory + '/pedestals_MIP' +'_run_' + str(self.run_number) + '.csv', index=False)


    def __plot_time_trend_no_beam(self, group):
        end, strip, layer = group.name
        print("creating plots: layer: ", layer, "bar: ", strip)
        ax = group.plot.scatter(x='pf_ticks', y='adc')
        ax.set_xlabel('pf_ticks')
        ax.set_ylabel('adc')
        ax.set_title('adc as function of pf_ticks in layer ' + str(layer) + ' bar ' + str(strip) + ' end ' + str(end))
        plt.savefig(self.plots_directory + '/time_trends' + '/time_trend_layer_' + str(layer) + '_bar_' + str(strip) + '_end_' + str(end) + '.pdf')
        plt.close()


    def __get_individual_pedestals_no_beam(self, group):
        layer, bar, end = group.name

        index = group.columns.str.contains('adc')
        data = group.iloc[:, index].to_numpy().flatten()

        pedestal = stats.mode(data)[0]
        self.out_ped_individ['layer_' + str(layer) + '_bar_' + str(bar) + '_end_' + str(end)] = pedestal

        # Identify peak widths
        # TODO: Note that number of bins are hard coded to 1023. i.e. how high the ADC goes, not very generic!
        bins = max(data) - min(data)
        hist = np.histogram(data, bins=bins)[0]
        # Under the assumption that the location in hist of max(hist) = pedestal
        hist_max = max(hist)
        peak_width = signal.peak_widths(hist, np.where(hist == hist_max)[0])[0]

        # cut away high and low noise
        data_cut = data[(data <= pedestal + peak_width * 3) & (data >= pedestal - peak_width * 3)]

        # Fit a Gaussian to the pedestals
        mean, std_dev = stats.norm.fit(data_cut)

        self.out_ped_individ['layer_' + str(layer) + '_bar_' + str(bar) + '_mean_end_' + str(end)] = mean
        self.out_ped_individ['layer_' + str(layer) + '_bar_' + str(bar) + '_std_dev_end_' + str(end)] = std_dev
        if(self.plot_pedestals):
            self.__plot_pedestal(data, layer, bar, end,
                                 self.out_ped_individ['layer_' + str(layer) + '_bar_' + str(bar) + '_mean_end_' + str(end)],
                                 self.out_ped_individ[
                                     'layer_' + str(layer) + '_bar_' + str(bar) + '_std_dev_end_' + str(end)],
                                 sum_pedestal=False, log_scale=True)

    def __get_sum_pedestal_no_beam(self, group):
        layer, bar, end = group.name

        # Obtain pedestal appropriate for summation of all 8 time samples
        pedestal_temp0 = 8 * self.out_ped_individ['layer_' + str(layer) + '_bar_' + str(bar) + '_end_' + str(end)]

        # Apply pedestal to sum of ADC events
        # note that this is just 0 for now since there is no sum of adc over any single event
        # TODO: find a reasonable substitute for how we calculate the sum in data with a beam for this data
        # TODO: investigate the use of the "raw_id" to group the datapoints into 8 size series. UPDATE: does not work! Id repeats!
        group['adc_sum_end' + str(end)] = 8 * group['adc'] - pedestal_temp0

        # Select sum pedestal region for plotting
        # Still unsure why we do this, but for the sake of consistency it's here
        fit_ = (group['adc_sum_end' + str(end)] < 200) & (group['adc_sum_end' + str(end)] > -200)
        fit = group[fit_]

        # Identify peak widths
        data = fit['adc_sum_end' + str(end)].to_numpy().flatten()
        # TODO: Note that the bins (and rage) are hard coded to 1023 * 8. i.e. how high the sumADC max goes, not very generic!
        bins = max(data) - min(data)
        hist = np.histogram(data, bins=bins)[0]
        # Under the assumption that the location in hist of max(hist) + min(data) = pedestal
        hist_max_loc = np.where(hist == max(hist))[0]
        peak_width = signal.peak_widths(hist, hist_max_loc)[0]

        # cut away high and low noise
        data_cut = data[(data <= hist_max_loc + min(data) + peak_width * 3) & (data >= hist_max_loc + min(data) - peak_width * 3)]

        del data

        # Fit a Gaussian to the pedestals
        mean, std_dev = stats.norm.fit(data_cut)

        del data_cut

        self.out_ped_sum['layer'].append(layer)
        self.out_ped_sum['strip'].append(bar)
        self.out_ped_sum['end'].append(end)
        self.out_ped_sum['pedestal'].append(pedestal_temp0)
        self.out_ped_sum['mean'].append(mean)
        self.out_ped_sum['std_dev'].append(std_dev)
        self.out_ped_sum['pedestal_per_time_sample'].append(
            self.out_ped_individ['layer_' + str(layer) + '_bar_' + str(bar) + '_end_' + str(end)])
        self.out_ped_sum['pedestal_per_time_sample_mean'].append(
            self.out_ped_individ['layer_' + str(layer) + '_bar_' + str(bar) + '_mean_end_' + str(end)])
        self.out_ped_sum['pedestal_per_time_sample_std_dev'].append(
            self.out_ped_individ['layer_' + str(layer) + '_bar_' + str(bar) + '_std_dev_end_' + str(end)])

        # Make plots
        if self.plot_pedestals:
            self.__plot_pedestal(group['adc_sum_end' + str(end)], layer, bar, end, mean, std_dev)

    #TODO: find a solution that works for calculating adc sums
    def __adc_sum_no_beam(self, group):
        layer, bar, end, raw_id = group.name
        self.adc_sum_dict['layer_' + str(layer) + '_bar_' + str(bar) + '_end_' + str(end)] = group['adc'].sum()

    def __pedestal_calc_no_beam(self, in_data):

        result_df = in_data[['layer', 'strip', 'end', 'raw_id', 'adc', 'tot', 'toa', 'pf_ticks']]
        del in_data
        gc.collect()

        # Select likely pedestal events where the TOA and TOT are both zero on each end of the bar
        # selection = (result_df['tot'] == 0) & (result_df['tot'] == 0)

        # result_df = result_df[selection]

        grouped_data_raw_id = result_df.groupby(['layer', 'strip', 'end', 'raw_id'])

        grouped_data_raw_id.apply(self.__adc_sum_no_beam)

        # Get individual pedestals
        grouped_data_individual = result_df.groupby(['layer', 'strip', 'end'], group_keys=False)

        grouped_data_individual.apply(self.__get_individual_pedestals_no_beam)

        # Now we will get the sum of ADC-appropriate pedestals
        grouped_data_sum = result_df.groupby(['layer', 'strip', 'end'], group_keys=False)

        grouped_data_sum.apply(self.__get_sum_pedestal_no_beam)

        if (self.time_trend):
            grouped_data_sum.apply(self.__plot_time_trend_no_beam)


    def get_pedestals_no_beam(self):
        print('running')
        # Loop through provided ROOT files
        for j in range(len(self.root_file_name)):

            if (self.layer_wise):
                layers = 19
                if (self.in_batches):
                    print('in batch wise')
                    for i in range(layers):
                        in_data = pd.DataFrame()
                        with uproot.open(self.root_file_name[j]) as in_file:
                            print('reading file for layer ', i + 1, '...')
                            cut = "layer == " + str(i + 1)
                            batchnbr = 1
                            for batch in in_file.iterate(
                                    ["layer", "end", "strip", "raw_id", "adc", "tot", "toa", "pf_event", "pf_spill",
                                     "pf_ticks"], cut, library="pd", step_size="10 MB"):
                                print("batch: ", batchnbr)
                                in_data = pd.concat([in_data, batch])
                                batchnbr += 1

                            print('finished reading file')

                        self.__pedestal_calc_no_beam(in_data)

                        # Create our final DataFrame with all pedestals
                        ped_df = pd.DataFrame(self.out_ped_sum)

                        # Save our pedestals to a csv file
                        print(self.out_directory + '/pedestals_no_beam' + '_run_' + str(self.run_number) + '.csv')
                        ped_df.to_csv(self.out_directory + '/pedestals_no_beam' + '_run_' + str(self.run_number) + '.csv', index=False)

                        del in_data
                        gc.collect()

                else:
                    print('in layers')
                    for i in range(layers):
                        with uproot.open(self.root_file_name[j]) as in_file:
                            print('reading file for layer ', i, '...')
                            cut = "layer == " + str(i + 1)
                            in_data = in_file.arrays(
                                ["layer", "end", "strip", "raw_id", "adc", "tot", "toa", "pf_event", "pf_spill",
                                 "pf_ticks"], cut,
                                library="pd")
                            print('finished reading file')

                        self.__pedestal_calc_no_beam(in_data)

                        # Create our final DataFrame with all pedestals
                        ped_df = pd.DataFrame(self.out_ped_sum)

                        # Save our pedestals to a csv file
                        print(self.out_directory + '/pedestals_no_beam' + '_run_' + str(self.run_number) + '.csv')
                        ped_df.to_csv(self.out_directory + '/pedestals_no_beam' + '_run_' + str(self.run_number) + '.csv', index=False)

                        del in_data
                        gc.collect()
            else:
                print('no layers or batches')
                with uproot.open(self.root_file_name[j]) as in_file:
                    print('reading file...')
                    in_data = in_file.arrays(
                        ["layer", "end", "strip", "raw_id", "adc", "tot", "toa", "pf_event", "pf_spill", "pf_ticks"],
                        library="pd")
                    print('finished reading file')

                # If we only want to look at one bar, define this here
                if self.do_one_bar is True and self.fpgas[j] == 0:
                    in_data = in_data[(in_data['layer'] == 1) & (in_data['strip'] == 3)]

                if self.do_one_bar is True and self.fpgas[j] != 0:
                    continue

                self.__pedestal_calc_no_beam(in_data)
                print('finished pedestal clalc')
                # Create our final DataFrame with all pedestals
                ped_df = pd.DataFrame(self.out_ped_sum)
                # Save our pedestals to a csv file
                ped_df.to_csv(self.out_directory + '/pedestals_no_beam' +'_run_' + str(self.run_number) + '.csv', index=False)
