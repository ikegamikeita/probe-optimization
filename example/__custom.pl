use 5.024;
use warnings;
use Mymodule::Custom;

my $data = make_data_structure( \{
				"oligomer" => parse_fasta( 'apmerge.mrna.sample.60mer' ),
				"bc"       => parse_MOST_tdt_simple( 'MOST_tdt.sample' ),
				"xhyb"     => parse_bowtie2( \{ 'path' => 'result.bwt2', 'max_mismatch' => 6 })
});

my $probe = select_probe( \{
			  'data'  => $data,
			  'model' => 'apmerge.mrna.sample',
			  'length_3end' => 1000,
			  'class_size'  => 60,
			  'max_probe'   => 3,
			  'max_xhyb'    => 500,
			  'uniq_id'     => 'myID_'
		     });

#option
my $hk = select_house_keeping_probe( \{
				     'tdt'      => $probe,
				     'hk_annot' => 'housekeeping/housekeeping_160627.annot',
				     'blast'    => 'housekeeping/result.blastx'
				});

