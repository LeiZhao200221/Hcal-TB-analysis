import pandas as pd
import numpy as np
import os

# Main class to perform baseline selection on events (output is a .csv file)
class selectCleanEvents:
    def __init__(self, data_file_name, out_directory='../analysis_files/', pedestal_file_name='../calibrations/pedestals.csv', mip_file_name='../calibrations/mip.csv', do_one_bar=False):
        '''
        Initalization
        @param data_file_name: str or list of str pointing to input analysis csv files
        @param out_directory: output directory for cleaned analysis csv files
        @param pedestal_file_name: location of pedestal calibration csv
        @param mip_file_name: location of MIP calibration csv
        @param do_one_bar: debug flag, performs chain on only one bar
        '''

        # Checks to see if the type of object passed as data_file_name is compatible with a request for alignment
        if(type(data_file_name) is not str and type(data_file_name) is not list):
            raise ValueError('Input file format should be a string of a single file name, or a list of strings of multiple file names!')
        self.data_file_name = None
        if(type(data_file_name)==str):
            self.data_file_name = [data_file_name]
        if(type(data_file_name)==list):
            self.data_file_name = data_file_name
        frames = []
        for i in range(len(self.data_file_name)):
            frames.append(pd.read_csv(self.data_file_name[i]))
        try:
            self.run_number = self.data_file_name[0].split('/')[-1].split('.')[0].split('_')[1]
        except:
            self.run_number = self.data_file_name[0].split('.')[0].split('_')[1]
        self.in_data = pd.concat(frames)
        self.out_directory = out_directory
        self.in_peds = pd.read_csv(pedestal_file_name)
        self.in_mips = pd.read_csv(mip_file_name)
        self.do_one_bar = do_one_bar

    # Clean DataFrame for easy handling
    def __clean_dataframes(self, result_df):
        # Make a new column which is the sum of MIPs on both ends of the bars
        result_df['mips']=result_df['adc_sum_end0'] + result_df['adc_sum_end1']

        # Make a new column which is the difference of TOAs on both ends of the bars
        result_df['toa']=result_df['toa_end0'] - result_df['toa_end1']

        # Select events where the TOA is non-zero on both ends of the bar
        result_df = result_df[(result_df['toa_end0']>0)& (result_df['toa_end1']>0)]

        # Group by 'event_number' and calculate the sum of ADC values for each event in layer 1
        layer_1_events = result_df[(result_df['layer'] == 1)]

        adc_sum_per_event = layer_1_events.groupby('pf_event')['mips'].sum()

        # Select events where the sum of ADC values falls within the specified range for layer 1
        selected_events = adc_sum_per_event[adc_sum_per_event.between(0.4, 4)].index

        # Filter the original DataFrame based on the selected events
        final_result_df = result_df[result_df['pf_event'].isin(selected_events)]

        # Remove irrelevant columns
        final_result_df = final_result_df.drop('toa_end0', axis=1)
        final_result_df = final_result_df.drop('toa_end1', axis=1)
        final_result_df = final_result_df.drop('adc_sum_end0', axis=1)
        final_result_df = final_result_df.drop('adc_sum_end1', axis=1)
        final_result_df = final_result_df.drop('adc_mean_end0', axis=1)
        final_result_df = final_result_df.drop('adc_mean_end1', axis=1)
        final_result_df = final_result_df.drop('adc_max_end0', axis=1)
        final_result_df = final_result_df.drop('adc_max_end1', axis=1)

        return final_result_df

    def __process_group(self, group):
        layer, bar = group.name

        print('layer: ', layer, ', bar: ', bar)

        # Select appropriate pedestal
        selection_ped0 = (self.in_peds['strip']==bar) & (self.in_peds['end']==0) & (self.in_peds['layer']==layer)
        selection_ped1 = (self.in_peds['strip']==bar) & (self.in_peds['end']==1) & (self.in_peds['layer']==layer)

        # Select appropriate MIP calibration
        selection_mips0 = (self.in_mips['layer']==layer) & (self.in_mips['strip']==bar) & (self.in_mips['end']==0)
        selection_mips1 = (self.in_mips['layer']==layer) & (self.in_mips['strip']==bar) & (self.in_mips['end']==1)

        ped0 = self.in_peds[selection_ped0]['pedestal'].iloc[0]
        ped1 = self.in_peds[selection_ped1]['pedestal'].iloc[0]

        mip0 = self.in_mips[selection_mips0]['mpv'].iloc[0]
        mip1 = self.in_mips[selection_mips1]['mpv'].iloc[0]

        # Convert sum of ADC into MIPs with pedestal subtracted
        group['adc_sum_end0'] = (group['adc_sum_end0'] - ped0) / mip0
        group['adc_sum_end1'] = (group['adc_sum_end1'] - ped1) / mip1

        return group

    # Main function to select clean events. Creates a csv that contains events with non-zero TOA on both ends of a bar and where the total energy
    # Deposited into the first layer is consistent with a MIP. Also aggregates ADC and TOA information into one column.
    def clean_events(self):
        
        os.makedirs(self.out_directory, exist_ok=True)

        # If we only want to look at one bar, define this here
        if self.do_one_bar is True:
            in_data = in_data[(in_data['layer']==1) & (in_data['strip']==3)]

        # Process each layer and bar independently
        grouped_data = self.in_data.groupby(['layer', 'strip'], group_keys=False)

        result_df = grouped_data.apply(self.__process_group)

        final_result_df = self.__clean_dataframes(result_df)

         # Save our pedestals to a csv file
        final_result_df.to_csv(self.out_directory+ '/cleaned.csv', index=False)
        