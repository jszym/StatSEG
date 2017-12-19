# -*- coding: utf-8 -*-
#
# __main__.py
# ========================
# An interactive command-line interface to the
# statseg nucleic/amino acid low-entropy
# masking library.
# ========================
#
# Copyright 2017 Joseph Szymborski
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import click
import statseg
import Bio.SeqIO


@click.command()
@click.option('--infile', help='FASTA file with your sequence.', required=True)
@click.option('--outfile', default=False,
              help='File to output masked result. Defaults to print to console.')
def statseg_ui(infile, outfile):
    """Interactive tool for the masking of low-entropy nucleic/amino acids."""

    records = []

    for seq_record in Bio.SeqIO.parse(infile, "fasta"):
        records.append(statseg.mask_low_complexity(seq_record)[0])
    
    if outfile:
        Bio.SeqIO.write(records, outfile, "fasta")
    else:
        for record in records:
            click.echo(">{}".format(record.id))
            click.echo(record.seq)

if __name__ == '__main__':
    statseg_ui()