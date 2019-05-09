import os
import csv
import sys

class consolidate_metrics(object):

    def __init__(self, mash_file, read_metrics_file, quast_file):
        self.sample_dict = {}
        self.mash_file = mash_file
        self.read_metrics_file = read_metrics_file
        self.quast_file = quast_file
    def parse_mash(self):
        """parse_mash reads the top hits from the mash output and stores the information within a dict. The sample file names are is the key and the value is a list containing the genus and serotype information"""
        with open(self.mash_file, newline='') as mash_file:
            mash_reader = csv.reader(mash_file, delimiter = ',')
            # skips the first row, which is the header of the file
            next(mash_reader)
            # parses through the csv file
            for row in mash_reader:
                # split the path of the sample file by '/'
                sample_path = (row[1].split('/'))
                # sample name
                sample, ext = (sample_path[-1].split('.'))
                # append the sample to the dict as the key
                self.sample_dict[sample] = []
                # mash id (i.e. organism info, genus, species name)
                organism_info_list = row[0]
                # split the genus and species
                genus_species_list = organism_info_list.split('-')
                # split again by the '.'
                final_genus_species_list = genus_species_list[7].split('.')
                # genus information
                organism = (final_genus_species_list[0].split('_'))
                # if the list size is greater than two, pop off the last element of the list and join the genus and species word by '_' 
                if len(organism) > 2:
                    organism.pop()
                    genus_species = '_'.join(organism)
                # if the list size is less than two, join the list and append it to the dict
                else:
                    genus_species = '_'.join(organism)
                # append the genus and species for the sample to a list
                self.sample_dict[sample].append(genus_species)
                
    def parse_read_metrics(self):
        """parses the master read metrics file to calculate the coverage of read1 and read2."""
        # forward sample dict
        forward_sample_dict = {}
        # reverse sample dict
        reverse_sample_dict = {}
        # open to read the master read metrics file
        with open(self.read_metrics_file, newline='') as read_metrics_file:
            read_metrics = csv.reader(read_metrics_file, delimiter = ',')
            # skips the first row, which is the header of the file
            next(read_metrics)
            # parses thorugh the read metrics file
            for row in read_metrics:
                # sample ID
                sample, ext = (((row[0].split('/'))[-1]).split('.'))
                # avgReadLength value
                avg_read_len = float(row[1])
                # totalBases value
                total_bases = int(row[2])
                # avgQuality vlaue
                average_quality = float(row[5])
                # number of reads value
                num_reads = int(row[6])
                # coverage vlaue
                estimated_coverage = float(row[8])
                # forward file stored based on the average quality, average read length, total bases, and number of reads
                if "_1" in sample:
                    sample, num = sample.split('_')
                    self.sample_dict[sample].append(average_quality)
                    self.sample_dict[sample].append(avg_read_len)
                    self.sample_dict[sample].append(total_bases) 
                    self.sample_dict[sample].append(num_reads)
                    estimated_coverage_read_1 = estimated_coverage
                # reverse file stored based on the average quality, average read length, total bases, and number of reads 
                elif "_2" in sample:
                    sample, num = sample.split('_')
                    self.sample_dict[sample].append(average_quality)
                    self.sample_dict[sample].append(avg_read_len)
                    self.sample_dict[sample].append(total_bases)
                    self.sample_dict[sample].append(num_reads)
                    # if there is an estimated coverage (i.e. for the forward sample)
                    if estimated_coverage_read_1:
                        # sum the coverage for the forward and reverse and append it to the list
                        total_coverage = estimated_coverage_read_1+estimated_coverage
                        total_coverage = round(total_coverage, 2)
                        self.sample_dict[sample].append(total_coverage)
                    else:
                        # append the coverage
                        self.sample_dict[sample].append(estimated_coverage)

    def parse_quast(self):
        """parses master quast output csv file, which has the sample, N50, N50 value for the sample"""
        with open(self.quast_file, newline='') as quast_file:
            # initialize reader object, which will read the csv file
            quast = csv.reader(quast_file, delimiter = ',')
            # parse the master csv file
            for row in quast:
                # first item in the list is the sample ID
                sample = row[0]
                # third itme in the list is the N50 vlaue
                N50_value = row[2]
                # N50 value is stored if the sample ID is in the dict
                if sample in self.sample_dict:
                    self.sample_dict[sample].append(N50_value)
                # if the sample ID is not within the dict, create a new key and store the N50 value
                else:
                     self.sample_dict[sample] = []
                     self.sample_dict[sample].append(N50_value)

    def write_report(self):
        """writes out the information harvested from mash output, read metrics output, and quast output. The file output will be in the following format:
        Sample, Predicted Genus, Predicted Species, R1_Q, R1_AvgReadLength, R1_TotalBases, R1_NumReads, R2_Q, R2_AvgReadLength, R2_TotalBases, R2_NumReads, Est_Cvg, N50   
        """
        with(open('AMR_COMP_Report.csv', 'a')) as comp_report:
            comp_report.write("Sample,Predicted_Species,R1_Q,R1_AvgReadLength,R1_TotalBases,R1_NumReads,R2_Q,R2_AvgReadLength,R2_TotalBases,R2_NumReads,Est_Cvg,N50")
            comp_report.write('\n')
            for sample, sample_info in self.sample_dict.items():
                # write the sample name
                comp_report.write(sample)
                # comma after the sample info
                comp_report.write(',')
                # stringify the list of sample information
                sample_info = [str(info_val) for info_val in sample_info]
                # join the sample information by ',' and write the info out to the file
                comp_report.write(','.join(sample_info))
                # write a new line for the next sample information
                comp_report.write('\n')
        
def main():
    try:
        # first input will be the master mash input file
        mash_file = os.path.basename(sys.argv[1])
        # second input will be the master read metrics file
        read_metrics_file = os.path.basename(sys.argv[2])
        # third input will be the master quast file
        quast_file = os.path.basename(sys.argv[3])
    except:
        quit("Please first input the Mash output file, Read Metrics output file, and Quast output file!")
    test = consolidate_metrics(mash_file, read_metrics_file, quast_file)     
    test.parse_mash()
    test.parse_read_metrics()
    test.parse_quast()
    test.write_report()

if __name__ == "__main__":
    main()
