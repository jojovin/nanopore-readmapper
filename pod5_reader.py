import pod5 as p5
import numpy as np
import matplotlib.pyplot as plt
import SAM_creater
import argparse

def main():
    argparser = argparse.ArgumentParser(description='Analysis of insertions in SAM file')
    argparser.add_argument('pod5', help='pod5 file')
    argparser.add_argument('sam', help='pickled sam object to load')
    args = argparser.parse_args()

    pod5 = args.pod5
    sam_file = args.sam

    SAM = SAM_creater.loadSAM(sam_file)

    SAMfile_readnames = [SAMline.READ_NAME for SAMline in SAM.SAMlines]

    with p5.Reader(pod5) as reader:
        for read in reader.reads(selection=SAMfile_readnames):
            assert str(read.read_id) in SAMfile_readnames

    showSingleReadPlot(pod5, SAMfile_readnames[0])
    
def showSingleReadPlot(pod5: str, selected_read_id: str):
    with p5.Reader(pod5) as reader:

        # Read the selected read from the pod5 file
        # next() is required here as Reader.reads() returns a Generator
        read = next(reader.reads(selection=[selected_read_id]))

        # Get the signal data and sample rate
        sample_rate = read.run_info.sample_rate
        signal = read.signal

        # Compute the time steps over the sampling period
        time = np.arange(len(signal)) / sample_rate # in seconds

        # Plot using matplotlib
        plt.plot(time, signal)
        plt.show()

if __name__ == "__main__":
    main()