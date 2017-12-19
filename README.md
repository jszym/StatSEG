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

StatSEG exposes an easy to use API and a CLI for quick console access. More in-depth documentation on both these systems exist in this repos Sphinx Documentation. 