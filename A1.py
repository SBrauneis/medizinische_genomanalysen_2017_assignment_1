import mysql.connector
import sys
import pysam
#import pybedtools

__author__ = 'Shelley Brauneis'
filename = "HG00096.chrom11.ILLUMINA.bwa.GBR.low_coverage.20120522.bam"

class Assignment1:
    def __init__(self):
        ## gene of interest is HBB, I assigned it as my gene
        self.gene = "WT1"
        self.header = []
        self.name = ""
        self.chrom = ""
        self.Start = ""
        self.End = ""
        self.strand = ""
        self.ExonCount = ""
        self.ExonStarts = ""
        self.ExonEnds = ""

    def fetch_gene_coordinates(self, genome_reference, file_name):
        print ("Connecting to UCSC to fetch data")

        ## Open connection
        cnx = mysql.connector.connect(host='genome-mysql.cse.ucsc.edu', user='genomep', passwd='password', db=genome_reference)

        ## Get cursor
        cursor = cnx.cursor()

        ## Build query fields
        query_fields = ["refGene.name2",
                        "refGene.name",
                        "refGene.chrom",
                        "refGene.txStart",
                        "refGene.txEnd",
                        "refGene.strand",
                        "refGene.exonCount",
                        "refGene.exonStarts",
                        "refGene.exonEnds"]

        ## Build query
        query = "SELECT DISTINCT %s from refGene" % ",".join(query_fields)

        ## Execute query
        cursor.execute(query)

        ## Write to file
        with open(file_name, "w") as fh:
            for row in cursor:
                fh.write(str(row) + "\n")
                if row[0] == self.gene:
                    header = row
                    #I will need this information later, so I save it in variables
                    self.gene = header[0]
                    self.name = header[1]
                    self.chrom = header[2]
                    self.Start = header[3]
                    self.End = header[4]
                    self.strand = header[5]
                    self.ExonCount = header[6]
                    self.ExonStarts = header[7]
                    self.ExonEnds = header[8]


        ## Close cursor & connection
        cursor.close()
        cnx.close()

        print ("Done fetching data")

    def get_sam_header(self): #print header(HD) and the info stored under it (version, sorting order, grouping)
        samfile = pysam.AlignmentFile(filename, "rb") #open samfile according to template: pysam.AlignmentFile("ex1.bam", "rb")
        for header, info in samfile.header['HD'].items():
            if header == 'VN':
                a=("Version: {}".format(info))
                return a
            elif header == 'SO':
                b=("Sorting order of alignments: {}".format(info))
                return b
            elif header == 'GO':
                c=("Grouping of alignments: {}".format(info))
                return c

        samfile.close()

    def get_properly_paired_reads_of_gene(self):
        # the file was indexed in the console for this step: samtools index HG00096.chrom11.ILLUMINA.bwa.GBR.low_coverage.20120522.bam
        samfile = pysam.AlignmentFile(filename, "rb")
        ProperPairNr = 0
        for read in samfile.fetch("11", self.Start, self.End):  # for the reads on chromosome 11 between my gene's start and end:
            if read.is_proper_pair:                             #call to check if read has proper pairs
                ProperPairNr = ProperPairNr + 1                 #counting number of proper pairs
                #print(read)                                    #printing reads with proper pairs, remove # if you want to see them printed
        samfile.close()
        return ProperPairNr

    def get_gene_reads_with_indels(self): #indels are D/I's in the CIGAR string
        samfile = pysam.AlignmentFile(filename, "rb")
        indels = 0
        for read in samfile.fetch("11", self.Start, self.End):
            if not read.is_unmapped:            #if the read is not unmapped
                cigar=read.cigar
                for type, nr in cigar:
                    if type == 1 or type == 2:
                        indels = indels + 1    #count number of indels
                        #print(read)            #print the read with an indel
        samfile.close()
        return indels

    def calculate_total_average_coverage(self):
        a = pybedtools.BedTool(filename)        #file was opened using pybedtools
        coverage = a.genome_coverage(bg=True)
        nr=0
        cov=0
        for line in coverage:
            ind_cov=float(line[3])
            cov += ind_cov
            nr += 1
            av_cov = cov/nr
        return av_cov

    def calculate_gene_average_coverage(self):
        a = pybedtools.BedTool(self.header)   #open in pybedtools again
        genecoverage=a.coverage(bg=True)    #calculate gene specific coverage
        nr = 0
        cov = 0
        avCov=0
        for line in coverage:
            st=int(line[1])
            end=int(line[2])
            if self.Start <= st:
                if self.End >= end:
                    ind_cov = float(line[3])
                    cov += ind_cov
                    nr += 1
        av_cov = cov / nr
        return av_cov


    def get_number_mapped_reads(self):
        samfile = pysam.AlignmentFile(filename, "rb")
        mapped=0
        for read in samfile:
            if not read.is_unmapped:    #for all "not unmapped" or better said: all mapped reads:
                mapped = mapped + 1     #count mapped reads
        samfile.close()
        return mapped

    def get_gene_symbol(self):
        return self.gene    #return the gene symbol from the header

    def get_region_of_gene(self):
        return("Chromosome: {}., Start: {}, End: {}".format(self.chrom, self.Start, self.End)) #print the chromosome the gene is on, then the start and end
        #this info was in the header

    def get_number_of_exons(self):
        return(self.ExonCount)  #return the exon count from the header

    def print_summary(self):
        print ("Results: (total of 9, this may take a while)")
        print("1. Sam Header:", self.get_sam_header())
        print("2. Proper Pairs:", self.get_properly_paired_reads_of_gene())
        print("3. Indels:", self.get_gene_reads_with_indels())
        #print("4. Total average coverage:", self.calculate_total_average_coverage())
        #print("5. Gene average coverage:", self.calculate_gene_average_coverage())
        print("6. Number of mapped reads:", self.get_number_mapped_reads())
        print("7. Gene symbol:", self.get_gene_symbol())
        print("8. Gene region:", self.get_region_of_gene())
        print("9. Number of exons:", self.get_number_of_exons())

if __name__ == '__main__':
    print ("Assignment 1", __author__)
    assignment1 = Assignment1()
    print("Fetching gene coordinates:")
    assignment1.fetch_gene_coordinates("hg19", "claudia.txt")
    assignment1.print_summary()