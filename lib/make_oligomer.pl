use 5.024;
use warnings;

#
# usage: perl make_oligomer.pl FASTA_file
#

my $genemodel = _parse_fasta( $ARGV[0] );

for my $gm ( sort { $a cmp $b } keys %{ $$genemodel } ) {
    my $pos_from_3prime;
    my $window = 60;
    my $limit  = 1000 - $window;
    for my $pos_from_3prime ( 0 .. length( $$genemodel->{ $gm }->{ seq } ) - $window ) {
	last if $pos_from_3prime > $limit;
	my $oligo60 = substr( $$genemodel->{ $gm }->{ seq }, ( length( $$genemodel->{ $gm }->{ seq } ) - $window - $pos_from_3prime ), $window );
	#
	# Remove oligomers containing non-A, -T, -G and -C.
	#
	if ( $oligo60 =~ m{[^ATGCatgc]} ) {
	    # non-A, -T, -G and -C
	}
	else {
	    say $$genemodel->{ $gm }->{ header } . ".pos=" . $pos_from_3prime;
	    say $oligo60;
	}
    }
}


sub _parse_fasta {
    #
    #
    #
    my @file = @_;
    my $data;
    local $/ = "\n>";
    for my $file ( @file ) {
	open  IN, "<", $file or die $!;
	while ( my $line = <IN> ) {
	    chomp $line;
	    $line =~ s{[>]*(.*?)\n}{};
	    my $header = $1;
	    $line =~ s{\n}{}g;
	    my $seq = $line;
	    $data->{ $header }->{ header } = ">" . $header;
	    $data->{ $header }->{ seq }    = $seq;
	}
	close IN;
    }
    return \$data;
}

