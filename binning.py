import sys
import os.path
import numpy as np

class Binning:
    def __init__(self, path_binaries, path_indicators, bin_size, pat_id):
        '''
        path_binaries : str, file path to where activity during film binaries are
        path_indicators: str, path to indicator function files
        bin_size: int, size of bin in msec (e.g. 100)
        pat_id = int, patient number
        '''
        self.path_binaries = path_binaries
        self.path_indicators = path_indicators
        self.bin_size = bin_size
        self.pat_id = pat_id
        self.pat_str = str(pat_id) + '/'

        self.rec_start, self.rec_stop = self.get_rec_refs()
        self.bins = self.make_bins()

    def load_binary(self, binary_name):
        binary = np.load(self.path_binaries + self.pat_str + binary_name)
        return binary

    def load_function(self, function_name):
        function = np.load(self.path_indicators + function_name)
        return function

    def get_rec_refs(self):
        rec_name = "rec_refs{}.npy".format(self.pat_id)
        recs = np.load(self.path_binaries + self.pat_str + rec_name)
        rec_start = recs[0]
        rec_stop = recs[1]
        return rec_start, rec_stop

    def make_bins(self):
        total_msec = self.rec_stop - self.rec_start
        total_bins = int(total_msec / self.bin_size)
        bins = np.linspace(self.rec_start, self.rec_stop, total_bins)
        return bins

    def save_bins(self):
        bin_dir = self.path_binaries + self.pat_str + str(bin_size)
        if not os.path.exists(bin_dir):
            os.makedirs(bin_dir)
        name = "edges_bin{}".format(self.bin_size)
        np.save(bin_dir + "/" + name, self.bins)


    def bin_spks(self):
        """
        binner for spike times

        outputs:
        saves a binary version of the histogram
            histogram[0]: values on the bins
            histogram[1]: bin edges
        """
        bin_dir = self.path_binaries + self.pat_str + str(self.bin_size)
        if not os.path.exists(bin_dir):
            os.makedirs(bin_dir)

        append_name = "_bin{}.npy".format(self.bin_size)
        for filename in os.listdir(self.path_binaries + self.pat_str):
            if filename.startswith('CSC'):
                csc_spikes = self.load_binary(filename)
                csc_binner, _ = np.histogram(csc_spikes, self.bins)
                print(len(csc_binner))
                print(csc_binner[0])
                print(csc_binner[1])
                name = filename[:-4] + append_name
                assert(os.path.isfile(bin_dir + "/" + name) is False), "Binary already exists."
                np.save(bin_dir + "/" + name, csc_binner)

    def bin_fxn(self):
        """
        binner for indicator functions
        """
        bin_dir = self.path_binaries + self.pat_str + str(self.bin_size) + '/indicators'
        if not os.path.exists(bin_dir):
            os.makedirs(bin_dir)

        append_name = "_bin{}.npy".format(self.bin_size)

        for filename in os.listdir(self.path_indicators):
            if filename.endswith('fxn.npy'):
                indicator_fxn = self.load_function(filename)
                indicator_binner, _ = np.histogram(indicator_fxn, self.bins)
                name = filename[:-4] + append_name
                assert(os.path.isfile(bin_dir + "/" + name) is False), "Binary already exists."
                np.save(bin_dir + "/" + name, indicator_binner)



if __name__ == '__main__':
    path_binaries = ('/home/al/Desktop/macke/DeepHumanVision_pilot/data_preprocessing/binaries/')
    path_indicators = ('/home/al/Desktop/macke/DeepHumanVision_pilot/affect_annotation/binaries/indicator/')
    bin_size = 100
    pat_id = 46

    binning = Binning(path_binaries, path_indicators, bin_size, pat_id)
    binning.bin_spks()
    binning.bin_fxn()
    binning.save_bins()
