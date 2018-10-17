#!/usr/bin/perl
use Getopt::Std;

my %options = ();
getopts("t:", \%options);
$tag=$options{t};

@lines=<>;

$print=0;
$vec=0;
for($i=0; $i<=$#lines; $i++){
    if( $lines[$i] =~ $tag ){
	$print=1;
    }
    if( $print==1 ){
	print $lines[$i];
	if($lines[$i] =~ /^ \$VEC/ ){
	    $vec =1;
	}
	if($lines[$i] =~ /^ \$END/ ){
	    $print =0;
	}
    }
}
if( $vec != 1) {
    die "VEC missing !!!";
}


