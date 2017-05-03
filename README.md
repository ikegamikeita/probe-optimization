# Citation

TBA


# Requirement

This is the authors' development environment. 

* Git (1.7.1)

* Perl (5.24.0)

* Bowtie2 (2.2.9)

* NCBI BLAST (2.2.31+)

* Firfox (50.0) recommended to browse the Agilent eArray.


# Usage and demonstration

1. Clone the repository to your working space.

		$ cd your/local
		$ git clone git@github.com:ikegamikeita/probe-optimization.git
		$ cd GIT_REPO


1. Following demonstration will be peformed at the example directory.

		$ cd example


1. Prepare oligomers from target sequence.

		$ perl ../lib/make_oligomer.pl apmerge.mrna.sample > apmerge.mrna.sample.60mer
		
		
1. Prepare the bowtie2 alignment result. Map oligomer sequence onto target sequence.
   	   
		$ bowtie2-build --threads 1 -f apmerge.mrna.sample apmerge.mrna.sample.index
		$ bowtie2 --threads 1 -f --very-sensitive --norc -x apmerge.mrna.sample.index -U apmerge.mrna.sample.60mer -S result.bwt2
		

1. Prepare the oligomer quality control using the Agilent eArray. First zip oligomer FASTA file (_e.g._ ` zip fasta.zip apmerge.mrna.sample.60mer `). Then manually upload the zipped FASTA file to eArray ([https://earray.chem.agilent.com/earray/](https://earray.chem.agilent.com/earray/)). Click "Probe" then "GE Probe Check" to upload the zipped FASTA. 

   ![image](https://i.gyazo.com/5644fca802f6070877ca7ff9ea272e34.png)

   You will receive an email from the Agilent eArray when process finished. Download the quality control result from the Agilent eArray and upload it to the exmple directory.

1. (Option - constitutional house keeping genes) If to set constitutional genes onto the custo array, blast your target sequence against the list of candidate house keeping genes (`housekeeping/housekeeping.prot`).

		$ makeblastdb -in housekeeping.prot -out housekeeping.prot.db -dbtype prot -hash_index -parse_seqids
		$ blastx -query ../apmerge.mrna.sample -db housekeeping_160627.prot.db -outfmt 6 -evalue 1e-1 -num_threads 1 > result.blastx

1. Set the environment variable PERL5LIB.

		#bash
		$ export PERL5LIB=your/local/probe-optimization/lib:$PERL5LIB


1. Calculate optimized probe set. Create a new file and open it with an editor.
   	     
		$ touch custom.pl
		$ emacs custom.pl

	Perl script:

	    use 5.024;
	    use warnings;
	    use Mymodule::Custom;

		my $data = make_data_structure( \{
											"oligomer" => parse_fasta( 'apmerge.mrna.sample.60mer' ),
											"bc"       => parse_MOST_tdt_simple( 'result_eArray/MOST_tdt.sample' ),
											"xhyb"     => parse_bowtie2( \{ 'path' => 'result.bwt2', 'max_mismatch' => 6 })
									});

		my $probe = select_probe( \{
									'data'  => $data,
									'model' => 'apmerge.mrna.sample',
									'length_3end' => 1000,   # ~~~~~~~~
									'class_size'  => 60,     # ~~~~~~~~
									'max_probe'   => 3,      # ~~~~~~~~
									'max_xhyb'    => 500,    # ~~~~~~~~
									'uniq_id'     => 'myID_' # ~~~~~~~~
									});

		#option
		#my $hk = select_house_keeping_probe( \{
		#										'tdt'      => $probe,
		#										'hk_annot' => 'housekeeping/housekeeping_160627.annot',
		#										'blast'    => 'housekeeping/result.blastx'
		#										});

	Finally, execute the script to obtain a set of optimized probes and summary. 
	
		$ perl custom.pl
	
	

