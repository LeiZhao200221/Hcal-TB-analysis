import pandas as pd
import numpy as np
import uproot
import os
import scipy.stats as stats
from statistics import mean,stdev,median
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.pyplot import get_cmap
import math
import mplhep as hep
from scipy.optimize import curve_fit
import pylandau
import csv

hep.style.use(hep.style.ATLAS)



# Main class to calculate and save MIPs
class calculateMIPs:
    def __init__(self, data_file_name, pedestal_file_name='../calibrations/pedestals.csv', mip_fit_cut_file_name=None, 
                 out_directory='.', plot_mips=False, plots_directory='../plots/mips', do_one_bar=False, calc_multipliers=True):
        '''
        Initialization
        @param data_file_name: str or list of str pointing to input analysis csv files
        @param pedestal_file_name: input pedestal csv calibration file
        @param mip_fit_cut_file_name: if using a csv containing how many standard deviations above pedestal to cut to fit MIPs
        @param out_directory: output directory for MIP calibration csv files
        @param plot_mips: flag indicating whether to plot the MIP fit results- if True, then plot
        @param plots_directory: where to put MIP fit plots
        @param do_one_bar: debug flag, performs chain on only one bar
        @param calc_multipliers: multipliers to try for by-hand cuts of pedestal
        '''
        self.out_directory = out_directory
        self.do_one_bar = do_one_bar
        self.plot_mips = plot_mips
        self.data_file_name = None

        if not isinstance(data_file_name, (str, list)):
            raise ValueError('Input file format should be a string of a single file name, or a list of strings of file names')

        if isinstance(data_file_name, str):
            self.data_file_name = [data_file_name]

        if isinstance(data_file_name, list):
            self.data_file_name = data_file_name

        frames = []
        for i in range(len(self.data_file_name)):
            frames.append(pd.read_csv(self.data_file_name[i]))

        
        try:
            self.run_number = self.data_file_name[0].split('/')[-1].split('.')[0].split('_')[1]
        except:
            self.run_number = self.data_file_name[0].split('.')[0].split('_')[1]
        # if(self.run_number != '287'):
          #  raise ValueError('Expecting run 287 for 4 GeV defocused muons!!')
        self.in_data = pd.concat(frames)
        self.in_peds = pd.read_csv(pedestal_file_name)
        self.calc_multipliers = calc_multipliers

        if os.path.exists(self.out_directory+'/mip.csv'):
            self.in_mips = pd.read_csv(self.out_directory+'/mip.csv')
        else:
            with open(self.out_directory+'/mip.csv', 'w', newline='') as csvfile:
                csv_writer = csv.writer(csvfile)
                csv_writer.writerow(['layer','strip','end','multiplier','mpv','eta','sigma','A'])
            self.in_mips = pd.read_csv(self.out_directory+'/mip.csv')

        if calc_multipliers is True:
            print(f'You have elected to calculate the pre-fit cuts by hand. The strategy is to select out enough of the pedestal distribution for a convergent MIP fit.')
            print(f'Be sure you are passing a list of multipliers to the ad_hoc_layer_multipler parameter. These are in units of standard deviations.')
            print(f'When calling the get_mips function, pass the layer (which_layer) and side (which_side) you want and a list of multipliers (ad_hoc_layer_multipler) corresponding to the bars associated to that layer.')
            print(f'Plotting will be done automatically and saved into plots/mips unless the plots_directory field is set explicitly.')
            print(f'A new csv file will automatically be created and populated with the multipliers you set in get_mips and will be overridden each time you iterate.')
            print(f'You can find your csv file here: calibrations/temp_mip_fit_cut.csv')
            self.plot_mips = True
            self.do_one_bar = False
            if os.path.exists('../calibrations/temp_mip_fit_cut.csv'):
                self.in_mip_fit_cut = pd.read_csv('../calibrations/temp_mip_fit_cut.csv')
            else:
                with open('../calibrations/temp_mip_fit_cut.csv', 'w', newline='') as csvfile:
                    csv_writer = csv.writer(csvfile)
                    csv_writer.writerow(['layer','strip','end','multiplier'])
                self.in_mip_fit_cut = pd.read_csv('../calibrations/temp_mip_fit_cut.csv')

        if calc_multipliers is False and mip_fit_cut_file_name is None:
            raise ValueError('You have elected to use per-calculated cuts for proper fit convergence. Please pass the proper csv file into the mip_fit_cut_file_name field.')

        if calc_multipliers is False:
            self.in_mip_fit_cut = pd.read_csv(mip_fit_cut_file_name)

        if not os.path.exists(self.out_directory):
            # Create the directory if it doesn't exist
            os.makedirs(self.out_directory)

        # Check plots directory and create if needed
        if plot_mips is True:
            if plots_directory is None:
                self.plots_directory = self.out_directory
            else:
                self.plots_directory = plots_directory
                if not os.path.exists(self.plots_directory):
                    # Create the directory if it doesn't exist
                    os.makedirs(self.plots_directory)

    def __do_one_end(self, dataframe, layer, bar, end, ad_hoc_layer_multipler):

        # Define selection
        selection_ped = (self.in_peds['strip']==bar) & (self.in_peds['end']==end) & (self.in_peds['layer']==layer)
        selection_mul = (self.in_mip_fit_cut['layer']==layer) & (self.in_mip_fit_cut['strip']==bar) & (self.in_mip_fit_cut['end']==end)
        selection_mips = (self.in_mips['layer']==layer) & (self.in_mips['strip']==bar) & (self.in_mips['end']==end)

        # Select appropriate pedestal
        pedestal_ = self.in_peds[selection_ped]['pedestal'].iloc[0]


        # Get our multiplier
        multiplier = None

        # If we are not deriving it, pull it from our csv
        if self.calc_multipliers is False:
            multiplier = self.in_mip_fit_cut[selection_mul]['multiplier'].iloc[0]

        # If we are deriving it, check to see if we already wrote it to our temp file and either write it or override it
        if self.calc_multipliers is True:
            multiplier = ad_hoc_layer_multipler[bar]
            matching_rows = self.in_mip_fit_cut[selection_mul].all(axis=1)
            if matching_rows.empty is False:
                self.in_mip_fit_cut.loc[selection_mul, 'multiplier'] = multiplier
            else:
                temp_cut = {'layer': [layer], 'strip': [bar], 'end': [end], 'multiplier': [multiplier]}
                new_row = pd.DataFrame(temp_cut)

                self.in_mip_fit_cut = pd.concat([self.in_mip_fit_cut, new_row], ignore_index=True)

        # Rescale ADC with pedestal subtraction
        dataframe['adc_sum_end'+str(end)] = dataframe['adc_sum_end'+str(end)]-(pedestal_)

        # Select data to fit the MIP peak according to our multiplier selection
        selection_ = (dataframe['adc_sum_end'+str(end)] > multiplier * self.in_peds[selection_ped]['std_dev'].iloc[0])
                
        dataframe = dataframe[selection_]

        # Do MIP fits
        # Make a histogram to fit to
        counts,bins = np.histogram(dataframe['adc_sum_end'+str(end)],bins=50,range=[0,2000])
        counts = np.array(counts)
        y = np.array(counts)
        x = (bins[:-1] + bins[1:]) / 2

        new_y = []
        new_x = []
        for kk in range(len(x)):
            if(y[kk]!=0):
                new_y.append(y[kk])
                new_x.append(x[kk])
        new_y = np.array(new_y)
        new_x = np.array(new_x)

        # Define outputs
        coeff = None
        pop_mips = {'layer': [layer],
                    'strip': [bar],
                    'end': [end],
                    'multiplier': [multiplier],
                    'mpv': [],
                    'eta': [],
                    'sigma': [],
                    'A': []}

        # Account for the case where the fit may fail and populate output lists
        try: # If fit converges, fill output with values
            mpv, eta, sigma, A = new_x[np.argmax(new_y)], 1, 1, max(new_y)

            # Fit with initial guesses
            coeff, pcov = curve_fit(pylandau.langau, new_x, new_y,
                        sigma = np.sqrt(new_y),
                        absolute_sigma = True,
                        p0 = (mpv, eta, sigma, A),
                        bounds=(1, 10000))

        except: # If fit does not converge, fill output with garbage
            print('Bad fit')
            coeff = [-999,-999,-999,-999]
            # Populate output with garbage
            
        # See if we already wrote this MIP calibration to our MIP csv file
        matching_rows = self.in_mips[selection_mips].all(axis=1)

        # If we did, just override existing values
        if matching_rows.empty is False:
            self.in_mips.loc[selection_mips, 'multiplier'] = multiplier
            self.in_mips.loc[selection_mips, 'mpv'] = coeff[0]
            self.in_mips.loc[selection_mips, 'eta'] = coeff[1]
            self.in_mips.loc[selection_mips, 'sigma'] = coeff[2]
            self.in_mips.loc[selection_mips, 'A'] = coeff[3]
            
        # If we didn't, append to exisiting MIP file
        else:
            pop_mips['mpv'].append(coeff[0])
            pop_mips['eta'].append(coeff[1])
            pop_mips['sigma'].append(coeff[2])
            pop_mips['A'].append(coeff[3])

            new_row = pd.DataFrame(pop_mips)

            self.in_mips = pd.concat([self.in_mips, new_row], ignore_index=True)

        # Make plots if requested
        if self.plot_mips is True and coeff is not None:
            fig,ax = plt.subplots(figsize=(8, 8))
            try:
                ax.plot(new_x, pylandau.langau(new_x, *coeff), "-")
            except:
                print('Cannot plot bad fit')
            ax.errorbar(x,y,yerr=np.sqrt(y),drawstyle='steps-mid')
            ax.set_yscale('log')
            ax.set_xlabel('Sum ADC '+str(end))
            ax.set_ylabel('Events')
            ax.set_ylim(1e0,1e4)
            ax.legend()
            ax.set_title('Layer: '+str(layer)+' Bar: '+str(bar)+' Side: 0')
            plt.savefig(self.plots_directory+'/mip_ped_subtracted_side_'+str(end)+'_layer_'+str(layer)+'_bar_'+str(bar)+'.pdf')
            plt.close()

    def __process_group(self, group, end, ad_hoc_layer_multipler):
        layer, bar = group.name

        print('layer: ', layer, ', bar: ', bar)

        # If we are trying to define multipliers, we will only look at the side requested
        if self.calc_multipliers is True:
            self.__do_one_end(group, layer, bar, end, ad_hoc_layer_multipler)

        # Else, we will look at both sides
        if self.calc_multipliers is False:
            self.__do_one_end(group, layer, bar, 0, ad_hoc_layer_multipler)
            self.__do_one_end(group, layer, bar, 1, ad_hoc_layer_multipler)

    # Main function to calculate MIPs
    def get_mips(self, which_layer=1, which_end=0, ad_hoc_layer_multipler=None):

        # Copy the data since this code allows the flexibility of real-time updating
        in_data_ = self.in_data.copy()
        
        # If we only want to look at one bar, define this here
        if self.do_one_bar is True:
            in_data_ = in_data_[(in_data_['layer']==1) & (in_data_['strip']==3)]

        if self.calc_multipliers is True:
            in_data_ = in_data_[(in_data_['layer']==which_layer)]

        grouped_data = in_data_.groupby(['layer', 'strip'], group_keys=False)
        grouped_data.apply(self.__process_group, end=which_end, ad_hoc_layer_multipler=ad_hoc_layer_multipler)

        # Write our MIP calibration file
        self.in_mips.to_csv(self.out_directory+'/mip.csv', index=False)

        # Write our multiplier definition file, if requested
        if self.calc_multipliers is True:
            self.in_mip_fit_cut.to_csv('../calibrations/temp_mip_fit_cut.csv', index=False)

        self.in_mips.to_csv(self.out_directory+'/mip.csv',index=False)