# bio-gag

bio-gag is a biogem for detecting and correcting a particular type of error that occurs in (at least) these particular versions of the IonTorrent sequencing kit:

* Ion Xpress Template 100 Kit
* Ion Xpress Template 200 Kit
* Ion Sequencing 100 Kit
* Ion Sequencing 200 Kit

*Gag error* is the term we've coined to describe an error that various people have observed on 
certain sequencing kits with IonTorrent, though it has not previously been characterised very well 
that we know of (we first noticed that the errors seemed to occur at GAG positions in the reads that 
were supposed to be GAAG). This biogem tries to find and fix these errors.

Errors that appear to be of this type were recently refered to in a benchtop sequencing platform 
comparison (Supplementary figure 4):

* http://www.ncbi.nlm.nih.gov/pubmed?term=22522955

There are also some more in-depth discussions about this on the (closed access) Ion Torrent forum: 

* http://ioncommunity.lifetechnologies.com/message/18047
* http://lifetech-it.hosted.jivesoftware.com/message/6233
* http://lifetech-it.hosted.jivesoftware.com/message/7893
* http://lifetech-it.hosted.jivesoftware.com/message/7792

To search for these errors, a pileup format file of aligned sequences is required. These can be generated
either from an assembly or by aligning to a reference, although it has only been tested on de-novo assemblies
assembled with newbler. Note that it has not been entirely optimised due to regular time constraints combined
with the fact they appear to have been fixed in newer kits.

## Installation

Firsts, you'll need to install the Ruby programming language. Instructions on how to do this are available at http://www.ruby-lang.org/

Then, install the bio-gag 'gem':

        gem install bio-gag

## Example usage

First test that installation worked correctly.

        gag --help

This should print out information on how to use this script. If that worked then you can predict errors:

        gag --pileup original_genome.pileup >gags.tsv

Where my.pileup is your de-novo assembly. This will predict gag errors in your assembly. You can also fix these errors

        gag --fix --pileup original_genome.pileup >genome_after_initial_gag_fix.fasta

However, the gag algorithm is not perfect, so spurious fixes may be applied. 
One way to reduce the incidence of spurious prediction is to only fix those errors that result in longer open reading
frames being predicted. Disadvantages of doing this include the fact that only errors within open reading frames
are fixed, and it requires that you are sequencing either a bacteria or an archaeon. The method currently requires
using the gene predictor [Prodigal](http://compbio.ornl.gov/prodigal/), for implementation reasons.

First, predict genes before and after the initial gag fixing by running Prodigal

        prodigal -q -i original_genome.fasta -a original_genome.proteins.fa -o /dev/null
        prodigal -q -i genome_after_initial_gag_fix.fasta -a genome_after_initial_gag_fix.proteins.fa -o /dev/null

Then filter the original gag errors, leaving only those that extend open reading frames

        gag --lookahead original_genome.proteins.fa genome_after_initial_gag_fix.proteins.fa gags.tsv >gags_filtered.tsv

There may be some ORFs that are more complicated than simply being lengthened or not.
These are output to the terminal - you'll have to decide for yourself
whether to keep these or not (sorry).

Then finally, do the final fix using the filtered gags file

        gag --fix original_genome.fasta -g gags_filtered.tsv >genome_after_filtered_gag_fix.fasta

That output file should hopefully have less gag errors in them than the original fasta file.

## Help

If you experience any problems using bio-gag, open an [issue](https://github.com/wwood/bioruby-gag/issues) on GitHub and tell us about it.

## Cite

Currently, this bio-gem is unpublished, but a relevant manuscript is in the works.

## Biogems.info

This Biogem is published at http://biogems.info/index.html#bio-gag

## Copyright

Copyright (c) 2012 Ben J Woodcroft. See LICENSE.txt for further details.

