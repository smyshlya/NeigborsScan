#!/usr/bin/perl  -w 

use strict;
use Bio::DB::EUtilities;
use Bio::Seq;
use List::MoreUtils qw(uniq);
my $debug=0;

#INPUT FILES/PARAMETERS

#this file is 'uniprot to gi' mapping table of the proteins of interest retrieved from http://www.uniprot.org/uploadlists/         
my $file='example.mapping_table'; 

# your NCBI API key (see README)
my $api_key = ''; 

# your email
my $email = ''; 

#the folder where you want to keep your results
my $results_folder='results';

#the length of adjacent regions (in nucleotides) to be scanned (in both directions)
my $range=20000; 

#the hmm database to be queried against                                                                                                                                                                                                                                                                          
my $hmm_folder='hmm';
my $hmm_profile='Conjscan.hmm';

#threshold e-value for hmmscan hits                                                                                                                                                                                                                                                                              
my $dom_set=1e-10;

#the hmm database to be queried against
my $hmm_folder='hmm';
my $hmm_profile='Conjscan.hmm'; 

#threshold e-value for hmmscan hits
my $dom_set=1e-10; 

#hex color definition for hmmscan hits
my $color="#000ff0"; #this is a default color for all the hits
my %color;

# below you can set the colors for the profile HMM hits $color{'profile NAME'} = 'hex color'
$color{'Phage_integrase'} = "#32a852"; 


my $start_run=time();
if (!-e $results_folder)
    {
        system ("mkdir $results_folder");
        print STDERR "all results will be saved into $results_folder\n"
    }
my $gis_hash=get_gis_from_mapping_table($file); 
my $uniprot_to_acc_num=acc_num_using_gis($gis_hash, $file);
my $nuc_hash=elink_nuc($gis_hash);
fetch_annotation_table_using_gi($results_folder, $nuc_hash);
my ($uniprot_to_nuc_hashref, $length_hash) = get_the_longest ($gis_hash, $nuc_hash, $results_folder, $file);
my ($uniprot_to_single_an)=find_an_for_uniprot($uniprot_to_nuc_hashref, $results_folder, $uniprot_to_acc_num, $file);
my ($global_start, $global_end, $global_strand)=download_proteins_from_range ($range,$uniprot_to_nuc_hashref, $results_folder, $uniprot_to_single_an);

run_hmmscan($results_folder, $hmm_profile, $gis_hash);
parse_hmmres($results_folder, $gis_hash, $nuc_hash, $file, $uniprot_to_single_an, $global_start, $global_end, $global_strand, $length_hash, $hmm_profile);

my $end_run=time();
my $run_time=$end_run - $start_run;
print STDERR "Finished. Job took $run_time seconds\n";

=head
run_hmmscan sub runs hmmscan ¯\_(ツ)_/¯ 
Creates .hmm_res file for each of the input proteins. 
=cut
sub run_hmmscan
{
    my ($results_folder, $hmm_profile, $gis_hash) = @_;
    if (!scalar(@_)==3)
    {
	die "number of arguments is incorrect for run_hmmscan\n";
    }
    print "\nhmmscan running..\n";
    if (!-e "$results_folder/$hmm_profile")
    {
	system ("mkdir $results_folder/$hmm_profile");
	print STDERR "results of hmmscan will be saved into $results_folder/$hmm_profile\n"
    }
    foreach my $keys (keys %$gis_hash)
    {
	my $protein_file="$results_folder/$keys.proteins";
	if (! -e "$results_folder/$hmm_profile/$keys.hmm_res")
	{
	    system ("hmmscan --domtblout $results_folder/$hmm_profile/$keys.hmm_res --notextw --domE $dom_set $hmm_folder/$hmm_profile $protein_file")
		and warn "hmmscan failed for $protein_file\n";
	}
    }
    print STDERR "done.\n\n";
    return;
}

=head
get_the_longest sub gets the longest nucleotide sequence corresponding to each of the uniprot acc num. 
The sub returns two hashes:
1.'uniprot number' -> 'GI of the longest nucleotide sequence' 
2.'uniprot number' -> 'the length of the longest sequence'
It also creates a file .long where it keeps the hashes for next runs of the sub on the same data. 
=cut
sub get_the_longest
{
    print STDERR "getting the longest sequence for each of the uniprot accession\n";
    my ($gis_hash, $nuc_hash, $results_folder, $file)=@_;
    if (!scalar( @_)==4)
    {
	die "number of arguments is incorrect for get_the_longest\n";
    }
    my @nuc;
    my %hash;
    my %length;
    if (!-e "$file.long")
    {
	open (my $out_fh, '>', "$file.long" );
	foreach my $key (keys %$nuc_hash)
	{
	    my $array=${$nuc_hash}{$key};
	    push @nuc, @$array;
	}    
	@nuc=uniq (@nuc);
	foreach (@nuc)
	{
	    my $nuc_length;
	    my $nuc=$_;
	    open (my $fh, '<', "$results_folder/$nuc.ft");
	    if ($debug){print "opening $results_folder/$nuc.ft\n";}
	    while ( my $row = <$fh> ) {
		chomp $row;
		my @split_string=split (/\t/, $row);	    
		if (defined $split_string[1] && length $split_string[1]>0 && defined $split_string[0] && length $split_string[0]>0)
		{
		    $split_string[1]=~s/(>|<|\^)//;
		    $split_string[0]=~s/(>|<|\^)//;
		    if ($split_string[1]>$split_string[0])
		    {
			$nuc_length=$split_string[1];
		    }
		    else 
		    {
			$nuc_length=$split_string[0];
		    }
		    
		}
	    }
	    if (defined $nuc_length && length $nuc_length>0)
	    {$length{$nuc}=$nuc_length;}
	    else {$length{$nuc}=0;}
	    if ($debug){print "length of $nuc is $length{$nuc}\n";}
	}
	my %longest_nuc;


	foreach my $key (keys %$nuc_hash)
	{
	    my $array=${$nuc_hash}{$key};
	    my $longest_length=0;
	    my $longest_nuc;
	    foreach (@$array)
	    { 
		my $length=$length{$_};
		if ($length>$longest_length)
		{
		    $longest_length=$length;
		    $longest_nuc=$_;
		}
	    }
	    if (defined $longest_nuc && length $longest_nuc > 0)
	    {
		$longest_nuc{$key}=$longest_nuc;
	    }
	    else
	    {
		$longest_nuc{$key}=0;
	    }
	    if ($debug){print "the longest nuc for protein $key is $longest_nuc{$key} - its length is $longest_length\n";}
	}

	foreach my $key (keys %$gis_hash)
	{
	    if ($debug){print "search for longest for $key\n";}
	    my $array=${$gis_hash}{$key};
	    my $longest_length=0;
	    my $longest_prot_gi;
	    foreach (@$array)
	    {
		my $prot_gi=$_;
		my $length=$length{$longest_nuc{$prot_gi}};
		if ($debug){print "prot gi is $prot_gi, corresponding longest nuc is $longest_nuc{$prot_gi} and its length is $length\n";}
		if ($length>$longest_length)
		{
		    $longest_length=$length;
		    $longest_prot_gi=$prot_gi;
		}
	    }	
	    if ($debug){print "the longest nuc for protein $key is $longest_prot_gi - its length is $longest_length\n";}
	    print $out_fh $key."\t".$longest_nuc{$longest_prot_gi}."\t".$longest_length."\n"; 
	    $hash{$key}=$longest_nuc{$longest_prot_gi};
	}
    }
    open (my $in_fh, '<', "$file.long");
    while (my $row=<$in_fh>) {
	chomp $row;
	my @split_string=split (/\t/, $row);
	$hash{$split_string[0]}=$split_string[1];
	$length{$split_string[0]}=$split_string[2];
    }
    print STDERR "done.\n\n";
    return \%hash, \%length;
}
=head
parse_hmmres sub creates tab-delimited file $file.out, which contains coordinates of hmmscan hits for each of the input protein. 
this file is designed to be used at http://itol.embl.de.
=cut
sub parse_hmmres
{
    print STDERR "\nparsing hmmscan results..\n";
    my ($results_folder, $gis_ref, $nuc_ref, $file, $uniprot_to_single_an, $start_hash, $end_hash, $strand_hash, $length_hash, $hmm_profile)= @_;
    if (!scalar( @_)==10)
    {
	die "number of arguments is incorrect for parse_hmmres\n";
    }
    $file=$file.".out";
    unlink "$file";
    open (my $res_fh, '>>', "$file");
    foreach my $uniprot (keys %$gis_ref)
    {
	my ($uniprot_acc_num, $uniprot_start,  $uniprot_end, $uniprot_strand, $length, $size, $start, $end)=get_uniprot_info($uniprot, $length_hash, $uniprot_to_single_an, $start_hash, $end_hash, $strand_hash);
	my $hits_ref=get_hits ($results_folder, $hmm_profile, $uniprot);

	my @hits = uniq (keys %$hits_ref);
	if ($debug) {print"keys are @hits\n";}
	if (!defined $uniprot_strand)
	{
	    print STDERR "failed to parse $uniprot\n";
	}
	else
	{
	    print $res_fh $uniprot."\t";
  	    if ($uniprot_strand=~m/plus/)
	    {
		print $res_fh "$size\tTR|$start|$end|#ff0000|int\t";
		foreach (@hits)
		{
		    my $hit=$_;
		    if (exists($color{${$hits_ref}{$hit}}))
		    {}
		    else
		    {
			$color{${$hits_ref}{$hit}}=$color;
		    }
			
#		    my @keys333 = uniq (keys %$start_hash);
#		    if($debug){print"keys to compare are: @keys333\n";}
		    if (!defined ${$start_hash}{$hit})
		    {
			print STDERR "failed to find coordinates for hmmscan hit $hit from $uniprot\n";
		    }
		    else
		    {
			if ($uniprot_start+$range>${$start_hash}{$hit} && $uniprot_start-$range<${$start_hash}{$hit})
			{
			    if (${$end_hash}{$hit}<$uniprot_start-$range)
			    {
				${$end_hash}{$hit}=$uniprot_start-$range+1;
			    }
			    if (${$end_hash}{$hit}>$uniprot_start+$range)
			    {
				${$end_hash}{$hit}=$uniprot_start+$range-1;
			    }
			    if ($debug) {print "acc number $hit has start ${$start_hash}{$hit}\n";}
			    my ($hmm_start, $hmm_end);
			    if ($uniprot_start>$range)
			    {
				if ($uniprot_start>${$start_hash}{$hit})
				{
				    $hmm_start=${$start_hash}{$hit}-($uniprot_start-$start);
				    $hmm_end=${$end_hash}{$hit}-($uniprot_start-$start);
				}
				else 
				{
				    $hmm_start=$start-($uniprot_start-${$start_hash}{$hit});
				    $hmm_end=$start-($uniprot_start-${$end_hash}{$hit});
				}
				if (${$strand_hash}{$hit}=~m/plus/)
				{
				    print $res_fh "PR|$hmm_start|$hmm_end|$color{${$hits_ref}{$hit}}|${$hits_ref}{$hit}\t";
				}
				if (${$strand_hash}{$hit}=~m/minus/)
				{
				    print $res_fh "PL|$hmm_end|$hmm_start|$color{${$hits_ref}{$hit}}|${$hits_ref}{$hit}\t";
				}
			    }
			    else 
			    {
				$hmm_start=${$start_hash}{$hit};
				$hmm_end=${$end_hash}{$hit};
				if (${$strand_hash}{$hit}=~m/plus/)
				{
				    print $res_fh "PR|$hmm_start|$hmm_end|$color{${$hits_ref}{$hit}}|${$hits_ref}{$hit}\t";
				}
				if (${$strand_hash}{$hit}=~m/minus/)
				{
				    print $res_fh "PL|$hmm_end|$hmm_start|$color{${$hits_ref}{$hit}}|${$hits_ref}{$hit}\t";
				}
			    }
			}
		    }
		}
	    }
	    else
	    {
		print $res_fh "$size\tTR|".($size-$start)."|".($size-$start+$length)."|#ff0000|int\t";
		foreach (@hits)
		{
		    my $hit=$_;
                    if (exists($color{${$hits_ref}{$hit}}))
                    {}
                    else
                    {
			$color{${$hits_ref}{$hit}}=$color;
                    }
		    if (!defined ${$start_hash}{$hit}) {print STDERR "failed to find coordinates for hmmscan hit $hit from $uniprot\n";}
		    else
		    {
			if ($uniprot_start+$range>${$start_hash}{$hit} && $uniprot_start-$range<${$start_hash}{$hit})
			{
			    if (${$end_hash}{$hit}<$uniprot_start-$range)
			    {
				${$end_hash}{$hit}=$uniprot_start-$range+1;
			    }
			    if (${$end_hash}{$hit}>$uniprot_start+$range)
			    {
				${$end_hash}{$hit}=$uniprot_start+$range-1;
			    }
			    if ($debug) {print "acc number $hit has start ${$start_hash}{$hit}\n";}
			    my ($hmm_start, $hmm_end);
			    $hmm_start=($size-$start)+$uniprot_start-${$start_hash}{$hit};
			    if ($hmm_start<0)
			    {
				print STDERR "something wrong with $uniprot, ";
			    }
			    $hmm_end=($size-$start)+$uniprot_start-${$end_hash}{$hit};
			    if (${$strand_hash}{$hit}=~m/plus/)
			    {
				print $res_fh "PL|$hmm_end|$hmm_start|$color{${$hits_ref}{$hit}}|${$hits_ref}{$hit}\t";
			    }
			    if (${$strand_hash}{$hit}=~m/minus/)
			    {
				print $res_fh "PR|$hmm_start|$hmm_end|$color{${$hits_ref}{$hit}}|${$hits_ref}{$hit}\t";
			    }
			}
		    }
		}
	    }
	print $res_fh "\n";
	}
    }
    print STDERR "done.\n\n";
}

#get_hits sub returns hash 'accession number' -> 'hmmscan hit'. It is used in parse_hmmres sub. 
sub get_hits
{
    my ($results_folder, $hmm_profile, $uniprot)=@_;
    if (!scalar( @_)==3)
    {
	die "number of arguments is incorrect for get_hits\n";
    }
    my $hmm_res="$results_folder/$hmm_profile/$uniprot.hmm_res";
    open (my $fh, "<", $hmm_res);
    my %hits;
    while (my $row=<$fh>) 
    {
	chomp $row;
	if ($row!~m/#.*/)
	{
	    my @split_string=split (/\s+/, $row);
	    if ($debug){print "$uniprot domain i-Evalue is $split_string[12]\n";}
	    if ($split_string[12]<$dom_set)
	    {
		my $name=$split_string[3];
		if($debug){print "corresponding name is $name and its hit is $split_string[0]\n";}
#		my @split=split (/\|/, $name);
		$name=~s/\.[1-9]//;
#		$hits{$split[3]}=$split_string[0];
		$hits{$name}=$split_string[0];
	    }
	}
    }
    return (\%hits);
}
=head
get_uniprot_info sub returns all necessary information for corresponding uniprot number,
such as accession number, length of corresponding nucleotide file etc.
It is used in parse_hmmres subroutine.
=cut
sub get_uniprot_info
{
    my ($uniprot, $length_hash, $uniprot_to_single_an, $start_hash, $end_hash, $strand_hash)=@_;
    if (!scalar (@_)==6)
    {
	die "number of arguments is incorrect for get_uniprot_info\n";
    }
    my $uniprot_nuc_length=${$length_hash}{$uniprot};
    my $uniprot_acc_num=${$uniprot_to_single_an}{$uniprot};
    my ($uniprot_start,  $uniprot_end, $uniprot_strand);
    if (!defined $uniprot_acc_num)
    {
	print STDERR "problems with $uniprot, can't find corresponding coordinates in feature file\n";
    }  
    else
    {
	$uniprot_start=${$start_hash}{$uniprot_acc_num};
	$uniprot_end=${$end_hash}{$uniprot_acc_num};
	$uniprot_strand=${$strand_hash}{$uniprot_acc_num};
	my ($size, $length, $start, $end);
	if ($uniprot_start>$range && $uniprot_start+$range<$uniprot_nuc_length)
	{
	    $start=$range;
	    $end=$range+$uniprot_end-$uniprot_start;
	    $size=2*$range;
	}
	elsif ($uniprot_start<$range && $uniprot_start+$range<$uniprot_nuc_length)
	{
	    $size=$uniprot_start+$range;
	    $start=$uniprot_start;
	    $end=$uniprot_end;
	}
	elsif ($uniprot_start<$range && $uniprot_start+$range>$uniprot_nuc_length)
	{
	    $size=$uniprot_nuc_length;
	    $start=$uniprot_start;
	    $end=$uniprot_end;
	    if ($start>$size) 
	    {
		print STDERR "something wrong, start is bigger than size for uniprot $uniprot (start is $start, and the size is $size)\n";
	    }
	}
	elsif ($uniprot_start>$range && $uniprot_start+$range>$uniprot_nuc_length)
	{       
	    $size=$range+$uniprot_nuc_length-$uniprot_start;
	    $start=$range;
	    $end=$range+$uniprot_end-$uniprot_start;
	    if ($start>$size) 
	    {
		print STDERR "something wrong, start is bigger than size for uniprot $uniprot (start is $start, and the size is $size)\n";
	    }
	}
	if ($uniprot_strand=~m/minus/)
	{
	    $length=$start-$end;
	}
	else
	{
	    $length=$end-$start;
	}
	return ($uniprot_acc_num, $uniprot_start,  $uniprot_end, $uniprot_strand, $length, $size, $start, $end);
    }
}

=head
find_an_for_uniprot sub retrieves unique genbank accession number for each of the uniprot accession number, using correspondin nuc file,
and returns a hash 'protein uniprot_number' -> 'protein genbank accession number'
=cut
sub find_an_for_uniprot
{
    print STDERR "getting unique accession number corresponding to each of the uniprot accession\n";
    my ($uniprot_to_nuc_hash, $results_folder, $uniprot_to_acc_num, $file)= @_;
    if (!scalar(@_)==4)
    {
	die "number of arguments is incorrect for find_an_for_uniprot\n";
    }
    my %uniprot_to_single_an;
    if (!-e "$file.uniprot_to_an")
    {
	open (my $out_fh, '>>', "$file.uniprot_to_an");
	foreach my $uniprot (keys %$uniprot_to_nuc_hash)
	{
	    my $fp="$results_folder/${$uniprot_to_nuc_hash}{$uniprot}.ft";
	    if (-e $fp)
	    {	    
		open my $fh, '<', $fp or die "can't read open '$fp'";
		my $stop;
		my ($start_hash, $end_hash, $strand_hash)=get_start_and_end_of_proteins($fh);
		foreach my $input_array (@{$uniprot_to_acc_num}{$uniprot})
		{
		    foreach (@$input_array)
		    {
			my $input_an=$_;
			my $stop=0;
			if (exists ${$start_hash}{$input_an} and $stop!=1)
			{
			    $stop=1;
			    $uniprot_to_single_an{$uniprot}=$input_an;
			}
		    }
		}
		print $out_fh "$uniprot\t$uniprot_to_single_an{$uniprot}\n";
	    }
	    else {print STDERR "failed to find nuc file corresponding to uniprot $uniprot, skipping..\n";}
	}
    }
    else
    {
	open (my $fh, '<', "$file.uniprot_to_an");
	while (my $row=<$fh>)
	{
	    chomp $row;
	    my @split_string=split (/\t/, $row);
	    $uniprot_to_single_an{$split_string[0]}=$split_string[1];
	}
    }
    print STDERR "done.\n\n";
    return (\%uniprot_to_single_an);
}

=head
download_proteins_from_range sub creates a .prot file with downloaded proteins found within $range from the beginning of input uniprot protein. 
It also creates a .coordinate file with positions of all downloaded proteins and returns them as hashes.
=cut
sub download_proteins_from_range
{
    my ($range, $uniprot_to_nuc_hash, $results_folder, $uniprot_to_single_an)=@_;
    if (!scalar (@_)==4)
    {
	die "number of arguments is incorrect for download_proteins_from_range\n";
    }
    print STDERR "downloading proteins..\n";
    if (!-e "$file.coordinates")
    {
	open (my $out_fh, '>>', $file.".coordinates");
	foreach my $uniprot (keys %$uniprot_to_nuc_hash)
	{
	    my $fp="$results_folder/${$uniprot_to_nuc_hash}{$uniprot}.ft";
	    if (-e $fp)
	    {
		open (my $fh, '<', $fp) or die "can't read open '$fp'";
		my ($start_hash, $end_hash, $strand_hash)=get_start_and_end_of_proteins($fh);
		my $count=0;
		foreach my $key (keys %$start_hash)
		{
		    print $out_fh "$key\t${$start_hash}{$key}\t${$end_hash}{$key}\t${$strand_hash}{$key}\n";
		}
	    }
	    else {print STDERR "problems with $fp corresponding to $uniprot\n";}
	}
    }
    foreach my $uniprot (keys %$uniprot_to_nuc_hash)
    {
	my (%tyrosine_rec_start, %tyrosine_rec_end);
	my $fp="$results_folder/${$uniprot_to_nuc_hash}{$uniprot}.ft";
	if ($debug) {print "found nuc file for $uniprot: $fp\n";}
	if (-e "$results_folder/$uniprot.proteins")
	{
	    if ($debug) {print "file $results_folder/$uniprot.proteins already exisits\n";}
	}
	else
	{
	    if (-e $fp)
	    {
		open my $fh, '<', $fp or die "can't read open '$fp'";
		my $stop;
		my ($start_hash, $end_hash, $strand_hash)=get_start_and_end_of_proteins($fh);
		open ( my $protein_fh, ">>", "$results_folder/$uniprot.proteins");
		my $count=0;
		my (%tyrosine_rec_start, %tyrosine_rec_end);
		my $input_an=${$uniprot_to_single_an}{$uniprot};
		if ($debug) {print "check if tyrosine $input_an exists in the file $fp\n";}
		if (exists ${$start_hash}{$input_an})
		{
		    open ( my $protein_fh, ">>", "$results_folder/$uniprot.proteins");
		    if ($debug) {print "$fp contains tyrosine rec $input_an: start - ${$start_hash}{$input_an}, end - ${$end_hash}{$input_an}\n";}
		    $tyrosine_rec_start{$fp}=${$start_hash}{$input_an};
		    $tyrosine_rec_end{$fp}=${$end_hash}{$input_an};
		    my @keys;
		    foreach my $key (keys %$start_hash) {
			if ($tyrosine_rec_start{$fp}-$range <= ${$start_hash}{$key} && ${$start_hash}{$key} <=$tyrosine_rec_start{$fp}+$range)
			{
			    $count++;
			    push(@keys, $key);
			    if ($debug){print "$key passed! (start is ${$start_hash}{$key}, end is ${$end_hash}{$key}), strand is ${$strand_hash}{$key}\n";}
			    my $temp_file=$fp.".".$key.".tmp";
#			    download_proteins_using_acc_num($temp_file, $key);
#			    open (my $temp_fh, '<', $temp_file);
#			    while ( my $line = <$temp_fh> ) {
#				print $protein_fh $line;
#			    }
			}
		    }
		    if ($debug) {print "total number of proteins is @keys\nsaving to $results_folder/$uniprot.proteins\n";}
		    download_proteins_using_acc_num($protein_fh, @keys);
		}
		print STDERR "$count proteins were downloaded from +/-$range bp from $uniprot\n";
		unlink glob "$results_folder/*.tmp";
	    }
	    else {print STDERR "fail to find table corresponding to uniprot accession $uniprot, skipping..\n"}
	}
    }
    open (my $coord_fh, '<', $file.".coordinates");
    my (%global_start, %global_end, %global_strand);
    while (my $row=<$coord_fh>)
    {
        chomp $row;
        my @split_string=split (/\t/, $row);
        $global_start{$split_string[0]}=$split_string[1];
        $global_end{$split_string[0]}=$split_string[2];
        $global_strand{$split_string[0]}=$split_string[3];
    }
    print STDERR "done.\n";
    return (\%global_start, \%global_end, \%global_strand);
}

=head
this sub employs EUtilities to download proteins, using their Genbank accession numbers.
it is used in download_proteins_from_range subroutine.
=cut
sub download_proteins_using_acc_num
{ 
    my ($temp_file, @id)=@_;
    if (!scalar(@_)==2)
    {
	die "number of arguments is incorrect for download_proteins_using_acc_num\n";
    }
    my $factory = Bio::DB::EUtilities->new(-eutil   => 'epost',
					   -db      => 'protein',
					   -id      => \@id,
					   -email   => $email,
					   -rettype => 'fasta',
					   -api_key => $api_key,
					   -keep_histories => 1
);
# not sure exactly what happens here, took it from https://bioperl.org/howtos/EUtilities_Cookbook_HOWTO.html 
# "How do I retrieve a long list of sequences using a query?"
    if (my $hist = $factory->next_History) {
	$factory->set_parameters(-eutil => 'efetch',
				 -rettype => 'fasta',
				 -history => $hist);
	my $retry = 0; my ($retmax, $retstart) = (500,0);
	$factory->set_parameters(-retmax => $retmax,
				 -retstart => $retstart);
	eval{
	    $factory->get_Response(-cb =>
				   sub {my ($data) = @_; print $temp_file $data} );
	};
	if ($@) {
	    die "Server error: $@.  Try again later" if $retry == 5;
	    print STDERR "Server error, redo #$retry\n";
	    $retry++ && redo RETRIEVE_SEQS;
	}
	print "Retrieved $retstart";
	$retstart += $retmax;
	
    }
    return;
}

=head
this sub returns three hashes from the .ft feature table:
1.'genbank accession number' -> 'start poisition of the protein in the nuc sequence'
2.'genbank accession number' -> 'end position in the nuc sequence'
3.'genbank accession number' -> 'strand in the nuc sequence'
It is used in two other subs: 'find_an_for_uniprot', and 'download_proteins_from_range'
=cut
sub get_start_and_end_of_proteins 
{
    if (!scalar( @_)==1)
    {
        die "number of arguments is incorrect for get_start_and_end_of_proteins\n";
    }
    my $fh=$_[0];
    my $start;
    my $end;
    my $CDS_start;
    my $CDS_end;
    my $acc_num;
    my $mark;
    my $strand="unknown";
    my %start_hash;
    my %end_hash;
    my %strand_hash;
    my @coordinates_array;
    while (my $row=<$fh>) {
	chomp $row;
	my @split_string=split (/\t/, $row);
	if (defined $split_string[0] && length $split_string[0]>0 && defined $split_string[1] && length $split_string[1]>0)
	{
	    $split_string[0]=~s/(>|<|\^)//;
	    $split_string[1]=~s/(>|<|\^)//;
	    if ($split_string[0]>$split_string[1])
	    {
		$start=$split_string[1];
		$end=$split_string[0];
		$strand="minus";
	    }
	    else 
	    {
		$start=$split_string[0];
		$end=$split_string[1];
		$strand="plus";
	    }
	}
	if (defined $split_string[2] && length $split_string[2]>0)
	{
	    $mark=$split_string[2];
	}
	if (defined $split_string[3] && length $split_string[3]>0 && $split_string[3]=~m/protein_id/ && $mark=~/CDS/)
	{
	    if ($strand =~m/plus/)
	    {
		$CDS_start=$start;
		$CDS_end=$end;
	    }
	    if ($strand=~/minus/)
	    {
		    $CDS_end=$start;
		    $CDS_start=$end;
		}
	    if ($split_string[4]=~m/(gb|ref|emb|dbj)\|(.*)\.[1-9]\|/)
    {
		$acc_num=$2;
		$start_hash{$acc_num}=$CDS_start;
		$end_hash{$acc_num}=$CDS_end;
		$strand_hash{$acc_num}=$strand;
	    }
	}
    }
    return (\%start_hash, \%end_hash, \%strand_hash);
    close $fh or die "can't read close $fh";
}

=head
acc_num_using_gis sub returns a hash 'uniprot number' -> 'array of genbank accession number' 
and saves it into tabulated file $file.acc_num
=cut
sub acc_num_using_gis
{
    if (!scalar(@_)==2)
    {
	die "number of arguments is incorrect for acc_num_using_gis\n";
    }
    my ($gis_hash, $file)=@_;
    my @ids;
    foreach my $keys (keys %$gis_hash)
    {
	my $array=${$gis_hash}{$keys};
	push @ids, @$array;
    }
    if ($debug) {print "@ids\n";}
    my %acc_num;
    print "retrieving accession numbers usin GI numbers\n";
    if (!-e $file.".acc_num")
    {
        open (my $acc_num_fh, '>', $file.".acc_num");
        my $factory = Bio::DB::EUtilities->new(-eutil => 'esummary',
                                               -email => $email,
                                               -db    => 'protein',
                                               -id    => \@ids);
        while (my $ds = $factory->next_DocSum) {
            if ($debug){print "ID: ",$ds->get_id,"\n";}
            my $id=$ds->get_id;
            print $acc_num_fh "$id\t";
            while (my $item = $ds->next_Item('flattened'))  {
                my $name=$item->get_name;
                my $content=$item->get_content;
                if ($debug){printf("%-20s:%s\n",$item->get_name,$item->get_content) if $item->get_content;}
                if ($name=~m/Caption/)
                {
                    print $acc_num_fh $content."\n";
                    push @{$acc_num{$id}}, $content;
                }
            }
        }
    }
    open (my $acc_num_fh, '<',$file.".acc_num");
    while (my $row = <$acc_num_fh>)
        {
            chomp $row;
            my @split_string=split (/\t/, $row);
            if ($debug){print "gi is $split_string[0], uniprot AccNum is $split_string[1]\n";}
            push @{$acc_num{$split_string[0]}}, $split_string[1];
        }
    my %uniprot_to_acc_num;
    foreach my $keys (keys %$gis_hash)
    {
	my $gis_array=${$gis_hash}{$keys};
	foreach (@$gis_array)
	{
	    my $array=$acc_num{$_};
	    foreach (@$array)
	    {
		push @{$uniprot_to_acc_num{$keys}}, $_;
	    }
	}
    }
    print STDERR "done.\n\n";
    return \%uniprot_to_acc_num;
}

#get_gis_from_mapping_table sub returns a hash 'protein unirpot number' -> 'array of protein gi numbers'
sub get_gis_from_mapping_table
{
    if (!scalar(@_)==1)
    {
	die "number of arguments is incorrect for get_gis_from_mapping_table\n";
    }
    print STDERR "retrieving GI numbers from file..\n";
    my %gis;
    my $file=$_[0]; 
    open (my $file_fh, '<', $file);
    while (my $row = <$file_fh>)
    {
	chomp $row;
	my @split_string=split (/\t/, $row);
	push @{$gis{$split_string[0]}}, $split_string[1]; 
    }
    print STDERR "done.\n\n";
    return (\%gis);
}

#here we return a hash 'uniprot number' -> 'array of corresponding nucleotide gis'
sub elink_nuc
{
    if (!scalar(@_)==1)
    {
        die "number of arguments is incorrect for elink_nuc\n";
    }
    my $gis_hash=$_[0];
    my @prot_gis;
    foreach my $keys (keys %$gis_hash)
    {
	my $array=${$gis_hash}{$keys};
	push @prot_gis, @$array;
    }
    my $factory = Bio::DB::EUtilities->new(-eutil  => 'elink',
					   -email  => $email,
					   -db     => 'nucleotide',
					   -correspondence => 1,
					   -dbfrom => 'protein',
					   -id     => \@prot_gis); 
    print STDERR "retrieving GIs of corresponding nucleotide sequences...\n";
    my %nucs;
    if (!-e $file.".elink")
    {  open (my $elink_fh, '>',$file.".elink");
       while (my $ds = $factory->next_LinkSet)
       {
	   if ($debug) {
	       print "Processing link name: ",$ds->get_link_name,"\n";
	   }
	   if ($debug) {
	       print "Protein ID: ",join(',',$ds->get_submitted_ids),"\n";
	   }
	   my @submitted_ids=$ds->get_submitted_ids;
	   my @nuc_id=$ds->get_ids;
	   if ($debug) {print "    First nuc IDs: $nuc_id[0]\n";}
	   push @{$nucs{$submitted_ids[0]}}, $nuc_id[0];
	   print "submitted protein is $submitted_ids[0]}, elinked nuc is  $nuc_id[0]\n";
	   print $elink_fh "@submitted_ids\t$nuc_id[0]\n";
       }
    }
    open (my $elink_fh, '<',$file.".elink");
        while (my $row = <$elink_fh>)
        {
            chomp $row;
            my @split_string=split (/\t/, $row);
            if ($debug){print "protein gi is $split_string[0], nuc gi is $split_string[1]\n";}
            push @{$nucs{$split_string[0]}}, $split_string[1];
        }
    print STDERR "done.\n\n";
    return (\%nucs);
}

#fetch_annotation_table_using_gi sub downloads .ft annotation tables for input nucleotide GIs
sub fetch_annotation_table_using_gi
{
    if (!scalar(@_)==2)
    {
        die "number of arguments is incorrect for fetch_annotation_table_using_gi\n";
    }
    print STDERR "fetching annotation files..\n";
    my ($results_folder, $nuc_hash)=@_;
    my @nuc;
    foreach my $key (keys %$nuc_hash)
    {
	my $array=${$nuc_hash}{$key};
	push @nuc, @$array;
    }
    @nuc=uniq(@nuc);
    my $c=0;
    foreach (@nuc)
    {
	my $nuc=$_;
	if ($debug){print "processing $nuc\n"}
	$c++;
	if (! -e "$results_folder/$nuc.ft")
	{
	    if ($debug) {print "sequence $results_folder/$nuc.ft does not exists\ndownloading the feature table";}
	    my $factory = Bio::DB::EUtilities->new(-eutil   => 'efetch',
						   -db      => 'nucleotide',
						   -id      => $nuc,
						   -email   => $email,
						   -rettype => 'ft',
						   -api_key => $api_key
						   
		);
	    if ($debug) {print "processing sequence $nuc ($c out of ". @nuc.")\ncheck if results already exists in results folder $results_folder...\n";}
            $factory->get_Response(-file => "$results_folder/$nuc.ft");
	}
	else {if ($debug) {print "sequence exists\n";}}
	open (my $fh, '<', "$results_folder/$nuc.ft");
	
	
    }
    print STDERR "done.\n\n";
    return;
}
__END__
    
