# phyloscanner
Analysing pathogen genetic diversity and relationships between and within hosts at once, in windows along the genome.  

<p align="center"><img src="InfoAndInputs/PhyloscannerDiagram_big4.jpg" alt="Phyloscanner" width="820", height="314"/></p>

phyloscanner's input is bam files: reads (fragments of nucleotide sequence) that have been mapped (aligned) to the correct part of some reference genome.
We wrote phyloscanner to analyse bam files that each represent a pathogen population in one host, exhibiting within-host and between-host diversity; in general use each bam file should be a sample representing some subpopulation, and we analyse within- and between-sample diversity.  

phyloscanner is freely available under the GNU General Public License version 3, described [here](LICENSE).  
phyloscanner runs natively on Linux and Mac OS, but not Windows.
However on any operating system (including Windows), if you have [VirtualBox](https://www.virtualbox.org/wiki/Downloads) installed, you can run [this](https://www.dropbox.com/sh/j3pmmunhxlc7g1w/AABddPfc5dN9oVnP9vQfAZOta?dl=0) image of Ubuntu Linux 16.04.3 which contains phyloscanner and our separate tool [shiver](https://github.com/ChrisHIV/shiver) (which allows you to process raw sequence data into mapped reads suitable as input for phyloscanner), with all of their dependencies.  
To make phylogenies from mapped reads, phyloscanner requires [samtools](http://www.htslib.org/), [pysam](https://github.com/pysam-developers/pysam) (0.8.1 or later), [biopython](http://biopython.org/wiki/Download), [mafft](http://mafft.cbrc.jp/alignment/software/) and [RAxML](http://sco.h-its.org/exelixis/web/software/raxml/index.html); notes on installing these are [here](InfoAndInputs/InstallationNotesForMakingTrees.sh).
To set up the part of phyloscanner that analyses these phylogenies (or others), follow [these instructions](InfoAndInputs/InstallationNotesForAnalysingTrees.sh). 

Info and help:
* The phyloscanner preprint, discussing the method and its scientific context, is [here](http://www.biorxiv.org/content/early/2017/06/30/157768).  
* The code's manual is [here](InfoAndInputs/PhyloscannerManual.pdf).  
* Instructions and example data for a practical on using phyloscanner are [here](https://tinyurl.com/PhyloscannerPractical).
* Problem with the code? Create a [New Issue](https://github.com/BDI-pathogens/phyloscanner/issues).  
* Query? Ask it publicly at the [google group](https://groups.google.com/forum/#!forum/phyloscanner-users).  

If you use phyloscanner for published work, please cite it and the tools it uses, details [here](InfoAndInputs/CitationDetails.bib).  

### An Example

The simulated bam files in [ExampleInputData](ExampleInputData) illustrate interesting within- and between-host diversity.
Let's analyse them with phyloscanner.
For this example I'll assume you've downloaded the phyloscanner code to your home directory, i.e. that its found in `~/phyloscanner/`.
First we need to make a file listing the files we want analysed; that file should be comma-separated variable format, containing lines like this: `BamFile1,ReferenceFile1,ID1` where the third field is optional (if present it is used as an ID for that bam file).
You can make that csv file manually if you like, but here is an example of making it automatically from the command line for these bam files:
```bash
for i in ~/phyloscanner/ExampleInputData/*.bam; do
  RefFile=${i%.bam}_ref.fasta;
  echo $i,$RefFile
done > InputFileList.csv
```
To make some within- & between-host phylogenies, run this command from the directory where your local copy of this code repostory lives (or adjusting paths appropriately if run from another directory):
```bash
~/phyloscanner/phyloscanner_make_trees.py InputFileList.csv --auto-window-params 300,-700,1000,8300 --alignment-of-other-refs ~/phyloscanner/InfoAndInputs/2refs_HXB2_C.BW.fasta --pairwise-align-to B.FR.83.HXB2_LAI_IIIB_BRU.K03455 
```
(Those window parameters make best use of this simulated data, the `-A` option includes an alignment of extra reference sequences along with the reads, see the manual for the `--pairwise-align-to` option.)  
Now let's analyse those phylogenies:
```bash
~/phyloscanner/phyloscanner_analyse_trees.R RAxMLfiles/RAxML_bestTree. MyOutput s,20 --outgroupName C.BW.00.00BW07621.AF443088 --multifurcationThreshold g
```
In the output you'll see trees and summary information indicating that these samples constitute:
* a straightforward, singly infected individual, 
* a dually infected individual i.e. infected by two distinct strains of virus,
* a contamination pair (one contaminating the other), and
* a pair exhibiting *ancestry*, i.e. one having evolved from the other. For populations of pathogens, this implies transmission, either indirectly (via unsampled intermediate patients) or directly.

Here is one of the trees from the output.
Each patient has many sequences (reads), with one colour per patient.
Extra reference sequences are coloured black.
Can you see why each of the patients merits their label?

<p align="center"><img src="InfoAndInputs/ProcessedTree_MyOutput_InWindow_4000_to_4300.jpg" alt="ExampleTree" width="800", height="741"/></p>

(These bam files were generated by taking HIV sequences from the [LANL HIV database](https://www.hiv.lanl.gov/content/sequence/NEWALIGN/align.html), using them as starting points for  evolution simulated with [SeqGen](https://github.com/rambaut/Seq-Gen) and calibrated to our own real data on within-host diversity.
Reads were then generated in eight windows of the genome, instead of all along the genome, purely to keep file sizes nice and small while remaining an interesting example.)
