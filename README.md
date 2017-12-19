# StatSEG
Stastical implementation of the SEG algorithm for the
masking of low-complexity amino/nucleic acids sequences.

## Installation

1. Clone the repo (`git clone https://github.com/jszym/statseg`)
2. Install requirements (`pip install -r requirements.txt`)

It's as easy as that.

## CLI Usage

Using StatSEG is easy, just specify a FASTA file with sequence that you want to mask using the `--infile` flag.

```
$ python -m statseg --infile prion.fasta

>sp|P04156|PRIO_HUMAN
MANLGCWMLVLFVATWSDLGLCKKRPKPGGxxxxxxxxPxxxSPGGNRYPPQGxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxQWNKPSKPKTNMKHMxxxxxxxxxxxxxxxYMLGSAMSRPIIHFGSDYEDRYYRENMHRYPNQVYYRPMDEYSNQNNFVHDCVNITIKQHTVTTTTKGExxxETDVKMMERVVEQMCITQYERESQAYYQRGSSMVLFSxxPVILLISFLIFLxxG
``` 

You can also output the masked sequence to a new FASTA file instead of just dumping it into the console.

```
$ python -m statseg --infile prion.fasta --outfile prion.masked.fasta
```

## Documentation

API & CLI documentation is [available here](http://bit.ly/2BbgtNt). An explanatory blog post is [available here](https://jszym.com/blog/quick-n-dirty-protein-dna-low-complexity-masking-implementation.html).