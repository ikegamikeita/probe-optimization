package Mymodule::Custom;

use base 'Exporter';
our @EXPORT = qw{
		    parse_MOST_tdt_simple
		    parse_bowtie2
		    make_data_structure
		    select_probe
		    select_house_keeping_probe
	    };


use 5.024;
use warnings;



sub select_house_keeping_probe {
    my $option = shift;
    my $tdt    = $$option->{ tdt };
    my $annot  = $$option->{ hk_annot };
    my $blast  = $$option->{ blast };
    my $evalue = $$option->{ evalue };
    
    # annotation file
    my $annot2;
    open  IN, "<", $annot or die $!;
    while ( my $line = <IN> ) {
	chomp $line;
	my @hoge = split /\t/, $line;
	$annot2->{ $hoge[2] }->{ ko } = $hoge[0];
	$annot2->{ $hoge[2] }->{ sp } = $hoge[1];
	$annot2->{ $hoge[2] }->{ description } = $hoge[3];
    }
    close IN;
    
    # get BLAST result and house keeping gene list
    my $hk;
    open  IN, "<", $$option->{ blast } or die $!;
    while ( my $line = <IN> ) {
	chomp $line;
	my @hoge = split /\t/, $line;
	
	# [Criteria 1] remove result if E-value > 1e-5.
	next if $hoge[10] > 1e-5;
	
	# [Criteria 2] remove result if E-value > 1e-40 or identity% < 50.
	next unless ( $hoge[10] < 1e-40 || $hoge[2] > 50 );
	
	next if $hoge[0] =~ m{apau\|gbr};
	$hk->{ $hoge[0] } = 1;
    }
    close IN;
    
    # Output in Agilent COMPLETE format //HOUSE KEEPING//
    open  OUT, ">", "probe." . $$tdt->{ note }->{ uniq_id } . ".probe=" . $$tdt->{ note }->{ max_probe } . ".xhyb=" . $$tdt->{ note }->{ max_xhyb } . ".hk.tdt"or die $!;
    for my $gm ( keys %{ $$tdt->{ probe } } ) {
	next unless $hk->{ $gm };
	my %c;
	@{ $$tdt->{ probe }->{ $gm } } = sort { $a cmp $b } @{ $$tdt->{ probe }->{ $gm } };
	@{ $$tdt->{ probe }->{ $gm } } = grep { !$c{$_}++ } @{ $$tdt->{ probe }->{ $gm } };
	for my $hoge ( @{ $$tdt->{ probe }->{ $gm } } ) {
	    for my $probe ( keys %{ $hoge } ) {
		say OUT join "\t", $hoge->{ $probe }->{ id } . "hk", $hoge->{ $probe }->{ seq }, $hoge->{ $probe }->{ target }, $hoge->{ $probe }->{ accession }, $hoge->{ $probe }->{ genesymbol }, $hoge->{ $probe }->{ description }, "";
	    }
	}
    }
    close OUT;
}


sub select_probe {
    #
    # usage:
    # input:
    #   data: make_data_structure(), data structure
    #   class_size: specify class window size
    #   length_3end: specify range of sequence from 3'-end to make oligomer
    #   max_xhyb: specify maximun value for cross-hybridization allowance
    #   max_probe: specify the number of maximum probes per target sequence
    #
    my $option = shift;
    my $model             = $$option->{ model };
    my $data              = $$option->{ data };
    my $class_length      = $$option->{ class_size };
    my $length_from3prime = $$option->{ length_3end };
    my $max_xhyb          = $$option->{ max_xhyb };
    my $max_probe         = $$option->{ max_probe };
    my $uniq_id           = $$option->{ uniq_id };
    
    # 1: make dictionary
    my $pos2class = pos2class( $length_from3prime, $class_length );
    
    # 2: make coordinates
    my $count;
    for my $gm ( keys %{ $$data } ) {
	for my $oligo ( keys %{ $$data->{ $gm } } ) {
	    $oligo =~ m{.pos=(.*)};
	    my $pos   = $1;
	    my $class = $$pos2class->{ $pos };
	    my $bc    = $$data->{ $gm }->{ $oligo }->{ bc };
	    my $xhyb  = $$data->{ $gm }->{ $oligo }->{ xhyb };
	    $count->{ $gm }->{ $xhyb }->{ $bc }->{ $class }->{ $pos } = $oligo;
	}
    }
    
    # 3:Select probes
    my $all_probe;
    my $probe;
    $probe->{ note }->{ class_size }  = $class_length;
    $probe->{ note }->{ length_3end } = $length_from3prime;
    $probe->{ note }->{ max_xhyb }    = $max_xhyb;
    $probe->{ note }->{ max_probe }   = $max_probe;
    $probe->{ note }->{ uniq_id }     = $uniq_id;
    for my $gm ( keys %{ $count } ) {
	my $check_class;
	my $probe_per_gm = 0;
	my $info;
	for my $xhyb ( sort { $a <=> $b } keys %{ $count->{ $gm } } ) {
	    
	    # Condition (1): when cross-hybridization == 0
	    if ( $xhyb == 0 ) {
		for my $bc ( sort { $a cmp $b } keys %{ $count->{ $gm }->{ $xhyb } } ) {
		    # Skip "POOR_BC" (or "NA")
		    next if $bc eq "POOR_BC";
		    next if $bc eq "NA";
		    for my $class ( sort { $a <=> $b } keys %{ $count->{ $gm }->{ $xhyb }->{ $bc } } ) {
			# Check each window only once per gene model per condition
			$check_class->{ $class } ++;
			next if $check_class->{ $class } > 1;
			
			# Check the number of the probes per gene model
			$probe_per_gm ++;
			next if $probe_per_gm > $$option->{ max_probe };
			
			# Count
			$all_probe ++;
			
			my @oligo_pos = sort { $a <=> $b } keys %{ $count->{ $gm }->{ $xhyb }->{ $bc }->{ $class } };
			my $oligo_id  = $count->{ $gm }->{ $xhyb }->{ $bc }->{ $class }->{ $oligo_pos[0] };
			#my $info;
			$info->{ $oligo_id }->{ id }          = sprintf( $uniq_id . "_%06s", $all_probe );
			$info->{ $oligo_id }->{ seq }         = $$data->{ $gm }->{ $oligo_id }->{ seq };
			$info->{ $oligo_id }->{ target }      = $gm;
			$info->{ $oligo_id }->{ accession }   = "";
			$info->{ $oligo_id }->{ genesymbol }  = "";
			$info->{ $oligo_id }->{ description } = $oligo_id . ",," . $bc . ",,Bin_" . $class_length . "=$class" . ",,xhyb=$xhyb";
			push @{ $probe->{ probe }->{ $gm } }, $info;
		    }
		}
	    }
	    
	    # Condition (2): when cross-hybridization > 0
	    elsif ( $xhyb <= $$option->{ max_xhyb } ) {
		for my $bc ( sort { $a cmp $b } keys %{ $count->{ $gm }->{ $xhyb } } ) {
		    # Skip "POOR_BC" (or "NA")
		    next if $bc eq "POOR_BC";
		    for my $class ( sort { $a <=> $b } keys %{ $count->{ $gm }->{ $xhyb }->{ $bc } } ) {
			# Check each window only once per gene model per condition
			$check_class->{ $class } ++;
			next if $check_class->{ $class } > 1;
			
			# Check the number of the probes per gene model
			$probe_per_gm ++;
			next if $probe_per_gm > $$option->{ max_probe };
			
			# Count
			$all_probe ++;
			
			my @oligo_pos = sort { $a <=> $b } keys %{ $count->{ $gm }->{ $xhyb }->{ $bc }->{ $class } };
			my $oligo_id  = $count->{ $gm }->{ $xhyb }->{ $bc }->{ $class }->{ $oligo_pos[0] };
			$info->{ $oligo_id }->{ id }          = sprintf( $uniq_id . "_%06s", $all_probe );
			$info->{ $oligo_id }->{ seq }         = $$data->{ $gm }->{ $oligo_id }->{ seq };
			$info->{ $oligo_id }->{ target }      = $gm;
			$info->{ $oligo_id }->{ accession }   = "";
			$info->{ $oligo_id }->{ genesymbol }  = "";
			$info->{ $oligo_id }->{ description } = $oligo_id . ",," . $bc . ",,Bin_" . $class_length . "=$class" . ",,xhyb=$xhyb";
			push @{ $probe->{ probe }->{ $gm } }, $info;
		    }
		}
	    }
	}
    }
    # Output in the Agilent COMPLETE format
    open  OUT, ">", "probe." . $probe->{ note }->{ uniq_id } . ".probe=" . $probe->{ note }->{ max_probe } . ".xhyb=" . $probe->{ note }->{ max_xhyb } . ".tdt"or die $!;
    for my $gm ( keys %{ $probe->{ probe } } ) {
	my %c;
	@{ $probe->{ probe }->{ $gm } } = sort { $a cmp $b } @{ $probe->{ probe }->{ $gm } };
	@{ $probe->{ probe }->{ $gm } } = grep { !$c{$_}++ } @{ $probe->{ probe }->{ $gm } };
	for my $hoge ( @{ $probe->{ probe }->{ $gm } } ) {
	    for my $probe ( keys %{ $hoge } ) {
		say OUT join "\t", $hoge->{ $probe }->{ id }, $hoge->{ $probe }->{ seq }, $hoge->{ $probe }->{ target }, $hoge->{ $probe }->{ accession }, $hoge->{ $probe }->{ genesymbol }, $hoge->{ $probe }->{ description }, "";
	    }
	}
    }
    close OUT;
    return \$probe;
}


sub pos2class {
    my $length_from3prime = shift;
    my $class_length = shift;
    my $class;
    my $data;
    for my $pos ( 0 .. $length_from3prime - $class_length ) {
	$class ++ if $pos % $class_length == 0;
	$data->{ $pos } = $class;
    }
    return \$data;
}


sub make_data_structure {
    #
    # usage:
    # input: keys => oligomer, bc_score, xhyb
    # return: reference
    #
    my $option   = shift;
    my $oligomer = $$option->{ oligomer };
    my $bc       = $$option->{ bc };
    my $xhyb     = $$option->{ xhyb };
    my $count;
    for my $oligo ( keys %{ $$oligomer } ) {
	$oligo =~ m{(.*).pos};
	my $gm = $1;
	$count->{ $gm }->{ $oligo }->{ xhyb } = $$xhyb->{ $oligo } // "NA";
	$count->{ $gm }->{ $oligo }->{ bc }   = $$bc->{ $oligo }   // "NA";
	$count->{ $gm }->{ $oligo }->{ seq }  = $$oligomer->{ $oligo }->{ seq };
    }
    return \$count;
}


sub parse_bowtie2 {
    #
    # usage:
    # path, max_mismatch, 
    #
    my $argument = shift;
    my $hit;
    open  IN, "<", $$argument->{ path } or die $!;
    while ( my $line = <IN> ) {
	next if $line =~ m{^@};
	chomp $line;
	my @option = split /\t/, $line;
	my $hit_evidence = $option[1];
	my $flag = 0;
	for my $i ( 11 .. $#option ) {
	    # The number of mismatches
	    if ( $option[ $i ] =~ m{XM:i:(.*)} ) {
		if ( $1 <= $$argument->{ max_mismatch } ) {
		    $flag ++;
		}
	    }
	    # The number of gap opens
	    $flag ++ if $option[ $i ] eq "XO:i:0";
	    # The number of gap extensions
	    $flag ++ if $option[ $i ] eq "XG:i:0";
	}
	if ( $flag == 3 ) {
	    #
	    # Mismatch within threshold && no gap open && no gap extention;
	    #
	    my $probe  = $option[0];
	    my $target = $option[2];
	    my $seq    = $option[9];
	    $probe =~ m{(.*).pos=.*};
	    my $genemodel = $1;
	    push @{ $hit->{ $probe } }, 1;
	    push @{ $hit->{ $probe } }, $target if $genemodel ne $target;
	}
    }
    close IN;

	# 以下、fileter_by_xhyb な処理。その２
    my $data;
    for my $probe ( sort { $a cmp $b } keys %{ $hit } ) {
	my %c;
	@{ $hit->{ $probe } } = sort { $a cmp $b } @{ $hit->{ $probe } };
	@{ $hit->{ $probe } } = grep { !$c{$_}++ } @{ $hit->{ $probe } };
	$data->{ $probe } = scalar @{ $hit->{ $probe } } -1;
    }

	return \$data;
}

sub parse_MOST_tdt_simple {
    #
    # from eArray::ProbeCheck;
    # //BC score only//
    #
    my $tdt    = shift; # MOST_tdt
    my $cutoff = shift || 0;
    my $data;
    my @header;
    open  IN, "<", $tdt or die $!;
    while ( my $line = <IN> ) {
	chomp $line;
	my @element = split /\t/, $line;
	$data->{ $element[0] } = $element[10];
	last if $. == $cutoff;#1_000;
    }
    close IN;
    return \$data;
}



1;
