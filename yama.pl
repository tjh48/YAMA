# Copyright (C) 2009 Thomas J. Hardcastle <tjh48@cam.ac.uk>

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# any later version.

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  You should have received a copy of 
#    the GNU General Public License along with this program. If not, see 
#    <http://www.gnu.org/licenses/>.


#!/usr/bin/perl
use List::Util qw(min);
use File::Path 'rmtree';

use strict ;
use warnings;
use Getopt::Long;
#!/usr/bin/env perl
$| = 1;

my $usage = "\nUsage: $0 -x GENOME_FASTA_FILE -Q SHORT_READ_FASTQ_FILE -o OUTPUT_DIR -n RESULT_NAME -t TEMPDIR --bowtie2=BOWTIE2_OPTIONS --bowtie=BOWTIE_OPTIONS --useBowtie --recordDuplication -i -d -r --no-split";
my $genome_file;
my $fastq_file ;
my $result_dir ;
my $result_name;
my $fastq_file_p1;
my $fastq_file_p2;
my $bowtie_options;
my $bowtie2_options;
my $noninverted;
my $useOnlyBest = 0;
my $deleteTemp;
my $temp_dir;
my $reuse_flag;
my $nosplit_flag;
my $useBowtie;
my $recordDup_flag;

# Get ARGS
GetOptions(
    'x:s'    => \$genome_file,
    'Q:s'    => \$fastq_file,
    '1:s'    => \$fastq_file_p1,
    '2:s'    => \$fastq_file_p2,
    'o:s'    => \$result_dir,
    'n:s'    => \$result_name,
    't:s'    => \$temp_dir,
    'bowtie2:s'    => \$bowtie2_options,
    'bowtie:s'    => \$bowtie_options,
    'useBowtie' => \$useBowtie,
    'recordDuplication' => \$recordDup_flag,
    'i'    => \$noninverted,
    'd'    => \$deleteTemp,
    'r'   => \$reuse_flag,
    'no-split' => \$nosplit_flag
) ;

$result_dir = "." unless $result_dir;
die "\n No output name given\n" unless $result_name;
$temp_dir = $result_dir.'/'."yama_tmp_$result_name" unless $temp_dir;
rmtree($temp_dir) if(-e $temp_dir && !$reuse_flag);
mkdir $temp_dir unless(-e $temp_dir);

my $logfile = $temp_dir.'/'."yama_log_$result_name";
unlink $logfile;

sub printlog {
    my $message = pop(@_);
    print $message;
    open(my $logFH, ">>", "$logfile") or die "cannot open > $logfile: $!";
    print {$logFH} $message;
    close($logFH);
}

sub dielog {
    my $message = pop(@_);
    open(my $logFH, ">>", "$logfile") or die "cannot open > $logfile: $!";
    print {$logFH} $message;
    close($logFH);
    die $message;
}


printlog "\nParameter choices:\n";
printlog "Ok, I won't invert the second element of each pair. This could be a horrible mistake...\n" if $noninverted;
printlog "Re-using temporary files.\n" if $reuse_flag;
printlog "Deleting temporary files afterwards.\n" if $deleteTemp;
if(!$useBowtie && $bowtie2_options && $bowtie_options) {
    printlog "Can't use both bowtie and bowtie2 - using bowtie2.";
    $bowtie_options = 0;
}
$bowtie_options = "" if($useBowtie && !$bowtie_options);
printlog "Using bowtie2 defaults." if(!$bowtie2_options && !$bowtie_options && !$useBowtie);	
printlog "Bowtie2 options: $bowtie2_options\n" if($bowtie2_options && !$useBowtie);
printlog "Bowtie options: $bowtie_options\n" if($bowtie_options || $useBowtie);
printlog "\n";

dielog "\n Bowtie can't handle paired-end alignments. Use the --bowtie2 option" if(($useBowtie || $bowtie_options) && ($fastq_file_p1 || $fastq_file_p2));
dielog "\n No genome file supplied.\n" unless ($genome_file);
dielog "\n No fastq file (or pair of fastq files) supplied.\n" unless ($fastq_file || ($fastq_file_p1 && $fastq_file_p2));

# make C2T genome file

my $genomeMethFile = "$genome_file"."_methConvert";

unless(-e $genomeMethFile) {

    printlog "Converting genome file to $genomeMethFile; cytosines->thymines\n";
    my $C2Tgenome_systemcall = join "", "sed '/^>/ s/^>.*/&_C2T/ ; /^>/! s/C/T/g' <",$genome_file," >", $genome_file,"_C2T";#
    system($C2Tgenome_systemcall);
    my $G2Agenome_systemcall = join "", "sed '/^>/ s/^>.*/&_G2A/ ; /^>/! s/G/A/g' <",$genome_file," >", $genome_file,"_G2A";#
    system($G2Agenome_systemcall);    
    system(join "", "cat $genome_file", "_C2T ", "$genome_file", "_G2A", ">$genomeMethFile") == 0 or dielog "Creation of genome file failed at concatenation stage: $?";
} else {
    printlog "Using existing converted genome file $genomeMethFile\n";
}

if($useBowtie || $bowtie_options) {
    unless(-e $genomeMethFile.".1.ebwt") {
	system("bowtie-build $genomeMethFile $genomeMethFile") == 0 or dielog "bowtie index creation failed: $?";
    }
} else {
    unless(-e $genomeMethFile.".1.bt2") {
	system("bowtie2-build $genomeMethFile $genomeMethFile") == 0 or dielog "bowtie2 index creation failed: $?";
    }
}

# subroutine for reading fastq files - faster (ha!) than BioPerl's native functions.

sub reverse_complement_IUPAC {
    my $dna = shift;
    
    # reverse the DNA sequence
    my $revcomp = reverse($dna);
    
    # complement the reversed DNA sequence
    $revcomp =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
    return $revcomp;
}


sub C2Tconversion {

    my $keepsame = pop(@_);
    my @fastq_files = @_;

    sub readfq {
	my ($fh, $aux) = @_;
	@$aux = [undef, 0] if (!defined(@$aux));
	return if ($aux->[1]);
	if (!defined($aux->[0])) {
	    while (<$fh>) {
		chomp;
		if (substr($_, 0, 1) eq '>' || substr($_, 0, 1) eq '@') {
		    $aux->[0] = $_;
		    last;
		}
	    }
	    if (!defined($aux->[0])) {
		$aux->[1] = 1;
		return;
	    }
	}
	my $name = /^.(\S+)/? $1 : '';
	my $seq = '';
	my $c;
	$aux->[0] = undef;
	while (<$fh>) {
	    chomp;
	    $c = substr($_, 0, 1);
	    last if ($c eq '>' || $c eq '@' || $c eq '+');
	    $seq .= $_;
	}
	$aux->[0] = $_;
	$aux->[1] = 1 if (!defined($aux->[0]));
	return ($name, $seq) if ($c ne '+');
	my $qual = '';
	while (<$fh>) {
	    chomp;
	    $qual .= $_;
	    if (length($qual) >= length($seq)) {
		$aux->[0] = undef;
		return ($name, $seq, $qual);
	    }
	}
	$aux->[1] = 1;
	return ($name, $seq);
    }

    my @fastqfiles_out;

    foreach my $fastfile (@fastq_files) {	
	my @aux = undef;
	my $slen = 0;
	my $fastfile_out = $fastfile;
	$fastfile_out =~ s/.*\///; 
	$fastfile_out = $temp_dir.'/C2T_'.$fastfile_out ;
	push(@fastqfiles_out, $fastfile_out);
	
# If the modified C2T file already exists, we'll just use that. Otherwise, make it!
	
	unless(-e $fastfile_out) {
	    
	    my $n = 0;
	    	    
	    open(my $fhfq, "<", "$fastfile") or dielog "cannot open < $fastfile: $!";
	    open(my $c2tfq, ">", "$fastfile_out") or dielog "cannot open > $fastfile_out: $!";
	    
	    printlog "Creating C2T fastq file: $fastfile_out...";
	    
	    my ($name, $seq, $qual);
	    while (($name, $seq, $qual) = readfq($fhfq, \@aux)) {    
		++$n;
		my $slen = length($seq);
		
#		 IDs in modified fastq file will tell us where the Cs were (minimally padded hexadecimal) before we did C2T conversion
#		my $newname = sprintf("%x", $n);	      		
#		$newname .= "_C";
#		my $padlength = length(sprintf("%x", $slen));
#		$newname .= join("", findletter("C", $seq, $padlength));
		
		if(($n % 1000000) == 0) {
		    print(".");
		}
		
		unless($keepsame) {
		    $seq = reverse_complement_IUPAC($seq);
		    $qual = reverse($qual);
		}

		my $newname = $seq;
		if(length($newname) > 255) {
		    printlog "This sequence ($seq) is too long for YAMA's encoding to work. I'm going to have to leave this one out.\n";
		    next;
		}

		$seq =~tr/C/T/ ;
		print {$c2tfq} "\@$newname\n$seq\n\+$newname\n$qual\n";
	    }
	    
	    printlog "done!\n";
	    
	    close($c2tfq);
	    close($fhfq);
	}
    }
    return(@fastqfiles_out);
}
    
sub processBase {
   
    my($genome_ref, $item_chr, $item_strand, $item_start, $seqLen, $Clocs_ref) = @_;

    my $baseCountFH;
    my $baseCountDupFH;
    my @Cpos;
    my @Ctype;
    my @referenceSeq;  
  
    @referenceSeq = split("", substr($genome_ref->{$item_chr}, $item_start - 1, $seqLen)) if ($item_strand eq "+");
    @referenceSeq = split("", reverse_complement_IUPAC(substr($genome_ref->{$item_chr}, $item_start - 1, $seqLen))) if ($item_strand eq "-");

    my @Cref = @referenceSeq[@$Clocs_ref];
    @Cpos = @$Clocs_ref[grep{$Cref[$_] eq "C"} 0..$#Cref];	            
    @Cpos = map {$_ = $seqLen - $_ -1} @Cpos if($item_strand eq "-");

    return(@Cpos);    
}

sub identifyType {
    
    my($genome_ref, $item_chr, $item_strand, $basePos) = @_;    

    my @Cpos;
    my @Ctype;
    my @bases;
    my $cytType = "";
  
    @bases = split("", substr($genome_ref->{$item_chr}, $basePos - 1, 3)) if ($item_strand eq "+");
    @bases = split("", reverse_complement_IUPAC(substr($genome_ref->{$item_chr}, $basePos - 3, 3))) if ($item_strand eq "-");

    if(defined $bases[1]) {
	$cytType = "CG" if($bases[1] eq "G");
	if(defined $bases[2]) {
	    $cytType = "CHG" if($bases[1] ne "G" && $bases[2] eq "G");
	    $cytType = "CHH" if($bases[1] ne "G" && $bases[2] ne "G");
	}
    }
#    print("$item_chr\t$basePos\t$item_strand\t$cytType\t", join " ", @bases, "\n");
    
    return($cytType);    
}


# This is the subroutine for printing out the best of bowtie alignments    

sub onlyTheBest {
    my ($baseCountFH, $baseCountDupFH, $paired, $genome_ref, $alignList_ref) = @_;    
    my $minqual;
    my $outline = "";

    if($useOnlyBest) {
	my @keep = ();
	foreach my $arrayIndex (0..(scalar(@$alignList_ref) - 1)) {	    
	    my $qual = $alignList_ref->[$arrayIndex][1];
	    $minqual //= $qual;	    
	    if($qual < $minqual) {
		@keep = ($arrayIndex);
		$minqual = $qual;
	    } elsif($qual == $minqual) {push @keep, $arrayIndex; }
	}
	@$alignList_ref = (@$alignList_ref)[@keep];
    }
	
    my $numAligns = scalar(@$alignList_ref);
    if($paired eq "=") {$numAligns = $numAligns / 2 ; }
    
    while(my $alignItem = shift(@$alignList_ref))
    {
	my ($trimID, $mapqual, $item_chr, $item_strand, $item_start, $cigar) = @$alignItem;
	my $itempos = $item_chr.":".$item_strand.":".$item_start;
	    
	my @cigarAction = grep length, split(/\d+/, $cigar);
	my @cigarLength = split(/[MID]/, $cigar);
	
	my $seq = $trimID;

	if(("I" ~~ @cigarAction) || "D" ~~ @cigarAction) {
	    if($item_strand eq "-") {$trimID = reverse($trimID); }
	    
	    my $seqpos = 0;
	    $seq = "";
	    my $cigarLen = @cigarAction;
	    for (my $cc = 0; $cc < $cigarLen; $cc++) {
		if($cigarAction[$cc] eq "M") {
		    $seq .= substr($trimID, $seqpos, $cigarLength[$cc]);			
		    $seqpos = $seqpos + $cigarLength[$cc];
		}
		$seq .= "-" x $cigarLength[$cc] if $cigarAction[$cc] eq "D";
		$seqpos = $seqpos + $cigarLength[$cc] if $cigarAction[$cc] eq "I";
	    }

	    if($item_strand eq "-") {$seq = reverse($seq); }
	}

	sub findletter {
	    my ($letter, $string) = @_;	    
	    my @stringArray = split "", $string;	    
	    my @letterpos = grep{$stringArray[$_] eq $letter} 0..$#stringArray;
	    return(@letterpos);
	}
	    
	my @Clocs = findletter("C", $seq);
	my @Tlocs = findletter("T", $seq);
	my $seqLen = length($seq);

	@Clocs = processBase($genome_ref, $item_chr, $item_strand, $item_start, $seqLen, \@Clocs);
	@Tlocs = processBase($genome_ref, $item_chr, $item_strand, $item_start, $seqLen, \@Tlocs);

	foreach my $Cloc (@Clocs) {
	    my $basePos = $Cloc + $item_start;
	    my $baseType = identifyType($genome_ref, $item_chr, $item_strand, $basePos);
	    print {$baseCountFH} "$item_chr\t$basePos\t$item_strand\t$baseType\t1\t0\t$numAligns\n";
	}
	foreach my $Tloc (@Tlocs) {
	    my $basePos = $Tloc + $item_start;
	    my $baseType = identifyType($genome_ref, $item_chr, $item_strand, $basePos);
	    print {$baseCountFH} "$item_chr\t$basePos\t$item_strand\t$baseType\t0\t1\t$numAligns\n";
	}
	if($numAligns > 1 && $recordDup_flag) { print {$baseCountDupFH} "$trimID\t$item_chr\t$item_start\t$item_strand\t$numAligns\n"; }	

	$outline .= join("\t", $itempos, join("_", @Clocs), join("_", @Tlocs), $seqLen, $numAligns);
	$outline .= "\n";
    }
    return $outline;
}
    
sub dec2bin {
    my $str = unpack("B32", pack("N", shift));
    # $str =~ s/^0+(?=\d)//;   # otherwise you'll get leading zeros
    return $str;
}

sub processLine {
    my($line, $mapqual, $alignList_ref) = @_;	   
    
    my ($id, $flag, $chr, $start, undef, $cigar, $paired, 
	undef, undef, undef, undef, @optionals) = split("\t", $line);
    my $strand = "\+";
    $strand = "-" if(substr(dec2bin($flag), -5, 1) == "1");
    
# reads must align to the right strand
    
    if(($chr =~ m/G2A/ && $strand eq "-") || ($chr =~ m/C2T/ && $strand eq "\+")) {	
	
# select only multiple matching reads with mismatches equal to the best performing (top) reported match	    	    	    
	
	$chr =~ s/_G2A//;
	$chr =~ s/_C2T//;		
	push @$alignList_ref, [$id, $mapqual, $chr, $strand, $start, $cigar];		
    }
    return($alignList_ref);
}



my @fastqfiles_out;
my @fastqfiles_out_p1;
my @fastqfiles_out_p2;
my $pairedend;

if($fastq_file) {
    $pairedend = 0;
    if($fastq_file_p1 || $fastq_file_p2) {
	printlog "YAMA doesn't currently handle simultaneous alignment of paired and single end reads. Only the single end reads will be used in this run.\n";
    }
    my @fastqfiles = ($fastq_file);
    @fastqfiles_out = C2Tconversion(@fastqfiles, 1);
} elsif($fastq_file_p1 && $fastq_file_p2) {
    $pairedend = 1;
    my @fastqfiles_p1 = ($fastq_file_p1);
    @fastqfiles_out_p1 = C2Tconversion(@fastqfiles_p1, 1);
    
    my @fastqfiles_p2 = ($fastq_file_p2);
    @fastqfiles_out_p2 = C2Tconversion(@fastqfiles_p2, $noninverted);
}


# again, if bowtie ouput file exists we'll just use that, otherwise, we'll use bowtie to align the reads

my $bowalignFile =  $temp_dir.'/'.$result_name.'_meth.bowalign' ;

unless(-e $bowalignFile) {
    my $bowtie_system_call;

    if($useBowtie || $bowtie_options) {
	my $bowtieinfiles = join ",", @fastqfiles_out;
	$bowtie_system_call = "bowtie $bowtie_options -q $genomeMethFile $bowtieinfiles --sam $bowalignFile";	
    } else {
	if($pairedend) {	
	    my $bowtieinfiles_1 = join ",", @fastqfiles_out_p1;
	    my $bowtieinfiles_2 = join ",", @fastqfiles_out_p2;
	    
	    $bowtie_system_call = "bowtie2 --no-unal --ff $bowtie2_options -q $genomeMethFile -1 $bowtieinfiles_1 -2 $bowtieinfiles_2 -S $bowalignFile";
	} else {
	    my $bowtieinfiles = join ",", @fastqfiles_out;
	    $bowtie_system_call = "bowtie2 --no-unal --ff $bowtie2_options -q $genomeMethFile -U $bowtieinfiles -S $bowalignFile";
	    unlink @fastqfiles_out if $deleteTemp;
	}
    }    
    printlog $bowtie_system_call, "\n";    
    system($bowtie_system_call) == 0 or dielog "bowtie(2) call failed: $?";

    unlink @fastqfiles_out if $deleteTemp && !$pairedend;
    unlink @fastqfiles_out_p1 if $deleteTemp && $pairedend;
    unlink @fastqfiles_out_p2 if $deleteTemp && $pairedend;

} else {
    printlog "Using existing bowtie alignment file: $bowalignFile\n";
}



open(my $bowalignFH, "<", "$bowalignFile") or dielog "cannot open < $bowalignFile: $!";

opendir(DIR, $temp_dir) || dielog "can't open directory : $!";
closedir(DIR);

my $baseCountFile = $temp_dir.'/'.$result_name."_basecount";
my $baseCountDupFile = $temp_dir.'/'.$result_name."_basecount_duplication";

# same-old same-old - if the file exists, we won't bother to create new ones.

# my $bestFile = $temp_dir.'/'.$result_name."_bestalignments";

unless(-e $baseCountFile)
{       
    open(my $baseCountFH , ">", $baseCountFile) or dielog "cannot open > $baseCountFile: $!";
    open(my $baseCountDupFH , ">", $baseCountDupFile) or dielog "cannot open > $baseCountDupFile: $!";

    my %chrSeqs;
    my $header;
    open(my $genomeFH, "<", $genome_file) or dielog "cannot open < $genome_file: $!";    
    printlog "Reading genome file...\n";

    while (<$genomeFH>) {
	chomp;
	if (substr($_, 0, 1) eq '>') {	       
	    $header = $_;
	    $header =~ s/\s.*$//;
	    $header =~ s/^>//;
	    dielog "Duplicate sequence for $header found in $genome_file!" if(exists($chrSeqs{$header}));
	    $chrSeqs{"$header"} = "";
	    print $header, "\n";	   	
	} else {
	    $chrSeqs{"$header"} .= $_;       
	}
    }

#    open(my $bestFileFH, ">", $bestFile) or dielog "cannot open > $bestFile: $!";    

    my @alignList = ();
    my $m = 0;    
    my $currid = "";
    my $currpairid = "";
    my $paired;

    while () {
	my $line = <$bowalignFH>;
	last if not defined $line;
	
	chomp $line;
	next if $line =~ m/^@/;

	my $pairline;	
	my $id; 
	my $qual;

	($id, $qual, $paired) = (split("\t", $line))[0, 1, 6];
	
	my $pairid = "";
	if($paired eq "=") {
	    $pairline = <$bowalignFH>;
	    my ($pairid, $pairqual) = (split("\t", $pairline))[0, 1];
	    $qual = $qual + $pairqual
	}

	if(scalar(@alignList) > 0) { $currid = $alignList[0][0]; }	    

	if($id ne $currid || $pairid ne $currpairid) {
	    onlyTheBest($baseCountFH, $baseCountDupFH, $paired, \%chrSeqs, \@alignList) if(scalar(@alignList) > 0);
	    @alignList = ();
	    $m++;
	    if(($m % 1000000) == 0) { print("."); }
	    $currid = $id;
	    $currpairid = $pairid;
	}
	
	processLine($line, $qual, \@alignList);
	if($paired eq "=") { processLine($pairline, $qual, \@alignList); }
	
    }
    
    onlyTheBest($baseCountFH, $baseCountDupFH, $paired, \%chrSeqs, \@alignList) if(scalar(@alignList) > 0);
    
    close($bowalignFH);
    close($baseCountFH);
    close($baseCountDupFH);
#    close($bestFileFH);

    printlog("done!\n");
} else {
    printlog "Basecount files already exist for this output name; using these:\n";
    printlog $baseCountFile, "\n";
}

unlink $bowalignFile if $deleteTemp;

my $baseCountSortFile = $baseCountFile;
$baseCountSortFile =~ s/basecount/basecount_sort/;

my $baseCountDupSortFile = $baseCountDupFile;
$baseCountDupSortFile =~ s/basecount/basecount_sort/;

unless(-e $baseCountSortFile) {
    printlog "Sorting basecount file...";
    system("sort -T $temp_dir -k1,1 -k2,2n -k3,3 -k7,7 $baseCountFile > $baseCountSortFile") == 0 or dielog "Sorting failed: $?";
    printlog "done!\n";
    
    printlog "Sorting duplicate basecount file...";
    system("sort -T $temp_dir -k1,1 -k2,2n -k3,3 $baseCountDupFile > $baseCountDupSortFile") == 0 or dielog "Sorting duplicate basecount file failed: $?";
    printlog "done!\n";
    
} else {
    printlog "Sorted basecount file already exists for this output name; using this:\n";
    printlog $baseCountSortFile, "\n";
}

unlink $baseCountFile if $deleteTemp;
unlink $baseCountDupFile if $deleteTemp;

#unless(-e $baseCountNRFile) {
printlog "Making sorted basecount non-redundant...";

open(my $SRfh, "<", $baseCountSortFile) or dielog "cannot open < $baseCountSortFile: $!";

my $NRfh;
my $CGNRfh;
my $CHGNRfh;
my $CHHNRfh;

if($nosplit_flag) {
    my $baseCountNRFile = $result_dir.'/'.$result_name."_methCalls";
    open($NRfh, ">", $baseCountNRFile) or dielog "cannot open > $baseCountNRFile: $!";
} else {
    my $baseCountNRCGFile = $result_dir.'/'.$result_name."_CG_methCalls";
    open($CGNRfh, ">", $baseCountNRCGFile) or dielog "cannot open > $baseCountNRCGFile: $!";
    my $baseCountNRCHGFile = $result_dir.'/'.$result_name."_CHG_methCalls";
    open($CHGNRfh, ">", $baseCountNRCHGFile) or dielog "cannot open > $baseCountNRCHGFile: $!";
    my $baseCountNRCHHFile = $result_dir.'/'.$result_name."_CHH_methCalls";
    open($CHHNRfh, ">", $baseCountNRCHHFile) or dielog "cannot open > $baseCountNRCHHFile: $!";
}

my $chr = "";
my $pos = -1;
my $strand = "";
my $context = "";
my $methCount = 0;
my $umethCount = 0;
my $nu = 0;
    
while (<$SRfh>) {
    chomp;
    my ($lchr, $lpos, $lstrand, $lcontext, $methc, $umethc, $lnu) = split("\t");
    if($lchr eq $chr && $lpos == $pos && $lstrand eq $strand)
    {
	$methCount += ($methc / $lnu);
	$umethCount += ($umethc / $lnu);	    	    
    } else {
	if($methCount > 0 || $umethCount > 0) {
	    if($nosplit_flag) {
		print {$NRfh} "$chr\t$pos\t$strand\t$methCount\t$umethCount\n";
	    } else {		
		if($context eq "CG") {
		    print {$CGNRfh} "$chr\t$pos\t$strand\t$methCount\t$umethCount\n";
		} elsif($context eq "CHG") {
		    print {$CHGNRfh} "$chr\t$pos\t$strand\t$methCount\t$umethCount\n";
		} elsif($context eq "CHH") {
		    print {$CHHNRfh} "$chr\t$pos\t$strand\t$methCount\t$umethCount\n";
		}
	    }
	}
	$chr = $lchr; $pos = $lpos; $strand = $lstrand; $methCount = $methc; $umethCount = $umethc; $nu = $lnu; $context = $lcontext;
    }	    
}
    
if($nosplit_flag) {
    print {$NRfh} "$chr\t$pos\t$strand\t$methCount\t$umethCount\n";
} else {		
    if($context eq "CG") {
	print {$CGNRfh} "$chr\t$pos\t$strand\t$methCount\t$umethCount\n";
    } elsif($context eq "CHG") {
	print {$CHGNRfh} "$chr\t$pos\t$strand\t$methCount\t$umethCount\n";
    } elsif($context eq "CHH") {
	print {$CHHNRfh} "$chr\t$pos\t$strand\t$methCount\t$umethCount\n";
    }
}


close($SRfh);

if($nosplit_flag) {
    close($NRfh);
} else {
    close($CGNRfh);
    close($CHGNRfh);
    close($CHHNRfh);
}


printlog "done!\n";
unlink $baseCountSortFile if $deleteTemp;


#} else {
#    print "Non-redundant basecount file already exists for this output name; using this:\n";
#    print $baseCountNRFile, "\n";
#}

#rmtree $temp_dir if $deleteTemp;
