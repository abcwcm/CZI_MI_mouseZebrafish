# the map genes scrip tis from https://github.com/atarashansky/SAMap 

./map_genes.sh --tr1 Mus_musculus.GRCm38.cds.all.fa --tr2 Danio_rerio.GRCz11.cds.all.fa --t1 nucl --t2 nucl --n1 mm --n2 dr

#
#
# Building a new DB, current time: 01/24/2022 18:22:55
# New DB name:   /athena/abc/scratch/paz2005/samap_directory/Mus_musculus.GRCm38.cds.all.fa
# New DB title:  Mus_musculus.GRCm38.cds.all.fa
# Sequence type: Nucleotide
# Keep MBits: T
# Maximum file size: 1000000000B
# Adding sequences from FASTA; added 68381 sequences in 3.56188 seconds.
#
#
# Building a new DB, current time: 01/24/2022 18:22:59
# New DB name:   /athena/abc/scratch/paz2005/samap_directory/Danio_rerio.GRCz11.cds.all.fa
# New DB title:  Danio_rerio.GRCz11.cds.all.fa
# Sequence type: Nucleotide
# Keep MBits: T
# Maximum file size: 1000000000B
# Adding sequences from FASTA; added 52089 sequences in 2.58909 seconds.
# Running tblastx in both directions



## the above code results in maps with transcript ids.  we need a map of genes to transcripts
# see https://github.com/atarashansky/SAMap/issues/54

# names : dict of list of 2D tuples or Nx2 numpy.ndarray, optional, default None
#           If BLAST was run on a transcriptome with Fasta headers that do not match
#           the gene symbols used in the dataset, you can pass a list of tuples mapping
#           the Fasta header name to the Dataset gene symbol:
#           (Fasta header name , Dataset gene symbol). Transcripts with the same gene
#           symbol will be collapsed into a single node in the gene homology graph.
#           By default, the Fasta header IDs are assumed to be equivalent to the
#           gene symbols used in the dataset.
#           The above mapping should be contained in a dicitonary keyed by the corresponding species.
#           For example, if we have `hu` and `mo` species and the `hu` BLAST results need to be translated,
#           then `names = {'hu' : mapping}, where `mapping = [(Fasta header 1, Gene symbol 1), ... , (Fasta header n, Gene symbol n)]`.


## 
grep ">" Danio_rerio.GRCz11.cds.all.fa | cut -d " " -f 1,4 | sed 's/>//g' | sed 's/gene://g' |  awk '{print $1,"\t"$2}' |  sed  's/\.[^.]*$//' | sed 's/[[:space:]]//' >  dr.names.tsv


grep ">"  Mus_musculus.GRCm38.cds.all.fa | cut -d " " -f 1,4 | sed 's/>//g' | sed 's/gene://g' | awk '{print $1,"\t"$2}' |  sed  's/\.[^.]*$//' | sed 's/[[:space:]]//'  >  mm.names.tsv


