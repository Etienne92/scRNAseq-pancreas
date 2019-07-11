import os
import sys
import json
import gzip


class ID2name_converter:
    """Converts gene IDs to gene names. A gtf file is required to create the mapping between gene IDs and gene names."""
    def __init__(self, gtfFile):
        filelocation=os.path.dirname(sys.argv[0]) +"/../data/id2name.json"

        #Only parse the gtf file if a json file has not been created yet.
        if os.path.isfile(filelocation):
            with open(filelocation,"r") as infile:
                self.id_to_name=json.load(infile)
        else:
            self.id_to_name={}
            with gzip.open(gtfFile,'rt') as gtf:
                for line in gtf:
                    if line[0]!="#":
                        gene_name=""
                        gene_id=""
                        info= line.split("\t")[8].split(";")
                        for x in info:
                            if x[0:7]=="gene_id":
                                gene_id = x.split("\"")[1]
                            if x[1:10]=="gene_name":
                                gene_name = x.split("\"")[1]
                        self.id_to_name[gene_id]=gene_name
            with open(filelocation, 'w') as outfile:  
                json.dump(self.id_to_name, outfile)
    def id2name(self,id):
        return self.id_to_name[id]

