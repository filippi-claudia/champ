#!/usr/bin/perl
# create sr p2nmcheck for checks of predefined 
# name-lists

# Input : LIST listname varname1 ... varnameN

$f77_indent='      ';
$f77_splt=  '     $ ';
$f77_scol=70;

%listvar=();
$lists=0;
$listname_len=0;
$vars=0;
$varname_len=0;

sub sep_by_comma{
    $s = @_[0];
    my @w=split(/,/,$s);
    for(my $j=0;$j<$#w;$j++){
	$w[$j] .= ',';
    }
    if(substr($s,length($s)-1,1) eq ','){
	$w[$#w] .= ',';
    }
    
    my @z = ();
    for(my $i=0; $i<=$#w; $i++){
	my @t=split(/\s+/,$w[$i]);
	for($j=0; $j<$#t; $j++){
	$t[$j] .= ' ';
    }
	push( @z, @t );
    }
    
    for(my $k =0; $k<=$#z; $k++){
	if($z[$k]  =~ /^\s*$/){
	    splice(@z,$k,1);
	    $k--;
	}
    }
    @z;
}



sub f77_print_line{
    ($lne_tmp) = @_[0];
    my @line_parts = sep_by_comma($lne_tmp);

    my $clen=length($f77_indent);
    
    $cline=$f77_indent;
    foreach my $q (@line_parts){
	if( ($clen+length($q)) <= $f77_scol ){
	    $clen += length($q);
	    $cline .= $q;
	} else {
	    print "$cline\n";
	    $cline=$f77_splt . $q;
	    $clen=length( $cline );
	    if($clen > $f77_scol){
		edit_error("cant split fortran line to fit");
	    }
	}
    }
    print "$cline\n";
}



sub read_input{
    $current_list='';
    INPUTLOOP: while(<>){
	if( /^s*#/) { next INPUTLOOP; }
	@w = split;
        if( $#w == -1) {next INPUTLOOP; }
	if($w[0] eq 'LIST'){
	    $current_list=$w[1];
	    $lists++;
	    if(length($current_list) > $listname_len){
                $listname_len=length($current_list);
            };
	    $i0=2;
	} else {
	    $i0=0;
	}
	for($i=$i0; $i<=$#w; $i++){
	    push @{$listvar{$current_list}}, $w[$i];
            $vars++;
            if(length($w[$i])>$varname_len){
                $varname_len=length($w[$i]);
            }
	}
    }
}

#foreach $i (keys %listvar){
#    printf("%s %d\n",$i,$#{$listvar{$i}});
#}

sub edit_error{
    $msg = @_[0];
    # printf STDERR ("Edit Error: %s : %d : %s\n",$ARGV,$.,$msg);
    # now ok in emacs compilation mode ?
     printf STDERR ("%s:%d: Error: %s\n",$ARGV,$.,$msg);
    $edit_errors++;
}  

sub generate_f77{
    print "C ---------- predefined namlist check -------\n";
    print "C this file is auto generated, do not edit   \n";
    f77_print_line("subroutine p2nmcheck(p,v,ierr)\n");
    f77_print_line("character p*(*), v*(*)\n");

    if($lists > 0){
	f77_print_line("character lists($lists)*($listname_len)\n");
	f77_print_line("character vars($vars)*($varname_len)\n");
	f77_print_line("dimension iaptr($lists),ieptr($lists)\n");

	$l="data lists/\'";
	$l .= join('\',\'', keys(%listvar));
	$l .= '\'/';
	f77_print_line($l);
	
	
	@pa=();
	$pe=();
	$ii=1;
	@varnames=();
	foreach $i (keys(%listvar)){
	    $s = $#{$listvar{$i}};
	    push @pa, $ii;
	    push @pe, $ii+$s;
	    $ii += ($s+1);
	    push @varnames,@{$listvar{$i}};
	}
	$l="data vars/\'";
	$l .= join('\',\'', @varnames);
	$l .= '\'/';
	f77_print_line($l);
	$l='data iaptr/';
	$l .= join(',',@pa);
	$l .= '/';
	f77_print_line($l);
	$l='data ieptr/';
	$l .= join(',',@pe);
	$l .= '/';
	f77_print_line($l);
	
	f77_print_line("nlist=$lists\n");
	
	@body=<DATA>;
	foreach $l (@body){
	    printf("%s",$l);
	}
    } else {
	f77_print_line("ierr=0\n");
	f77_print_line("end\n");
    }
}


read_input();
generate_f77();

__DATA__
      ierr=0
      do i=1,nlist
       if(lists(i).eq.p) then
        do iv=iaptr(i),ieptr(i)
         if(vars(iv).eq.v) then
          return
         endif
        enddo
        ierr=1
        return
       endif
      enddo
      return
      end
