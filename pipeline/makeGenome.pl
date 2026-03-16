#!/bin/env perl

##
## ecLego3 v3.2404.01
## Copyright (C) 2019-2026 Chee-Hong WONG (Dr Chia-Lin Wei Laboratory)
##
## ecLego component
##

use strict;
use Getopt::Long;
use File::Basename;
use FindBin;
use lib $FindBin::Bin;
use Data::Dumper;

my $G_USAGE = "
$0 --seq <seqfile> --path <path> --vcf <vcffile>
  --seq STR           sequence file (.fasta)
  --vcf STR           variant call file(.vcf)
  --path STR          genome path
  --var STR           variant ids
  --oprefix STR       prefix to output filename
";

# TODO: handle nested approach in future... path{variant}
my $seqFile = undef;
my $vcfFile = undef;
my $pathid = undef;
my $paths = undef;
my $prefix = undef;
my $verbose = 0;
my $help = 0;
my $rotate = 0;

GetOptions (
"seq=s"        => \$seqFile,
"vcf=s"        => \$vcfFile,
"id=s"         => \$pathid,
"path=s"       => \$paths,
"oprefix=s"    => \$prefix,
"rotate=i"     => \$rotate,
"verbose!"     => \$verbose,
"help!"        => \$help)
or die("Error in command line arguments\n$G_USAGE");

die "$G_USAGE" if ($help);

if (!defined $pathid && defined $prefix) {
    $pathid=basename($prefix);
}


# extract variants and path components
my %seqs = ();
my @paths = ();
my %variants = ();
# contig_2:Sniffles2.INV.387S1
foreach my $path (split(/;/, $paths)) {
    my ($seqId, $muts) = split(/:/, $path, 2);

    my $id = $seqId;
    $id =~ s/^[\+\-]//;
    if (!exists $seqs{$id}) {
        my %item = (seqid=>$id, seq=>'');
        $seqs{$id} = \%item;
    }

    $muts = '' if (!defined $muts);
    my @muts = split(/,/, $muts);
    foreach my $mut (@muts) {
        if (!exists $variants{$mut}) {
            my %item = (variant=>$mut);
            $variants{$mut} = \%item;
        }
    }

    my %path = (seqid=>$seqId, id=>$id, mutations=>\@muts);
    push @paths, \%path;
}

# read seqfile for components and tally
open INFILE, "$seqFile" || die "Fail to open $seqFile\n$!\n";
my $id = '';
my $seqRef = undef;
while (<INFILE>) {
    chomp();
    if (/^>/) {
        $id = substr($_, 1);
        ($id) = split(/\s+/, $id);
        $seqRef = (exists $seqs{$id}) ? $seqs{$id} : undef;
    } else {
        if (defined $seqRef) {
            $seqRef->{seq} .= $_;
        }
    }
}
close INFILE;

while (my ($id, $seqRef) = each %seqs) {
    die "$id does not have sequence info in $seqFile!\n" if (0==length($seqRef->{seq}));
}

# read vcffile for variants and tally
if (defined $vcfFile) {
    open INFILE, "$vcfFile" || die "Fail to open $vcfFile\n$!\n";
    while (<INFILE>) {
        next if (/^#/);
        chomp();
        my @bits = split(/\t/);
        if (exists $variants{$bits[2]}) {
            my $varRef = $variants{$bits[2]};
            $varRef->{CHROM} = $bits[0];
            $varRef->{POS} = $bits[1];
            $varRef->{REF} = $bits[3];
            $varRef->{ALT} = $bits[4];
            $varRef->{QUAL} = $bits[5];
            $varRef->{FILTER} = $bits[6];
            my %info = ();
            foreach my $infoBit (split(/;/, $bits[7])) {
                my ($key, $value) = split(/\=/, $infoBit);
                $info{$key} = $value;
            }
            $varRef->{INFO} = \%info;
        }
    }
    close INFILE;
}

while (my ($id, $varRef) = each %variants) {
    next if ($id =~ /del\.\d+\-\d+/); # skip operation that is non-vcf related
    next if ($id =~ /rotate\.[+\-]*\d+/); # skip operation that is non-vcf related
    next if ($id =~ /ins\.\d+\..+/); # skip operation that is non-vcf related
    die "$id does not have info in $vcfFile!\n" if (!exists $varRef->{CHROM});
    die "Variant $id $varRef->{CHROM} is not on your path!\n" if (!exists $seqs{$varRef->{CHROM}});
}
foreach my $pathRef (@paths) {
    foreach my $variant (@{$pathRef->{mutations}}) {
        next if ($variant =~ /del\.\d+\-\d+/); # skip operation that is non-vcf related
        next if ($variant =~ /rotate\.[+\-]*\d+/); # skip operation that is non-vcf related
        next if ($variant =~ /ins\.\d+\..+/); # skip operation that is non-vcf related
        my $varRef = $variants{$variant};
        die "Variant $id ($varRef->{CHROM}) on path $pathRef->{seqid} ?!\n" if ($pathRef->{id} ne $varRef->{CHROM});
    }
}


# generate output
my $finalSeq = '';
foreach my $pathRef (@paths) {
    my $pathSeq = '';
    my $toRCSeq = 0;
    my $state = substr($pathRef->{seqid},0,1);
    my $pid = $pathRef->{seqid};
    if ("+" eq $state) {
        $pid = substr($pid, 1);
    } elsif ("-" eq $state) {
        $pid = substr($pid, 1);
        $toRCSeq = 1;
    }
    $pathSeq = $seqs{$pid}->{seq};
    # TODO: check!
    # apply the variant
    my $nonVCFINS = 0;
    my $VCFINS = 0;
    my $cratedRotation = 0;
    my $VCFDEL = 0;
    my $nonVCFDEL = 0;
    my $numMuts = scalar(@{$pathRef->{mutations}});
    if ($numMuts>0) {
        my %muttypes = ();
        foreach my $mut (@{$pathRef->{mutations}}) {
            if ($mut =~ /del\.\d+\-\d+/) {
                # non-vcf related
                $muttypes{DEL}++;
                $nonVCFDEL = 1;
            } elsif ($mut =~ /rotate\.[+\-]*\d+/) {
                # non-vcf related
                $muttypes{ROT}++;
                $cratedRotation = 1;
            } elsif ($mut =~ /ins\.\d+\..+/) {
                # non-vcf related
                $muttypes{INS}++;
                $nonVCFINS = 1;
            } else {
                $muttypes{$variants{$mut}->{INFO}->{SVTYPE}}++;
                $VCFDEL=1 if ('DEL' eq $variants{$mut}->{INFO}->{SVTYPE});
                $VCFINS=1 if ('INS' eq $variants{$mut}->{INFO}->{SVTYPE});
            }
        }
        my $numMutTypes = scalar(keys %muttypes);
        if ($numMutTypes<=1) {
            # okie, proceed
        } else {
            # TODO: check INS!
            if (2==$numMutTypes) {
                # Q: composite deletion supported should have been 1==$numMutTypes !
                if (exists $muttypes{DEL} && exists $muttypes{ROT}) {
                    # deletion and rotation combination is fine
                } elsif (exists $muttypes{INS} && exists $muttypes{ROT}) {
                    # insertion and rotation combination is fine
                } elsif (exists $muttypes{INS} && exists $muttypes{DEL}) {
                    # insertion and deletion combination is fine
                } else {
                    die "Sorry, but multiple mutations on a single contig has not been implemented.\n";
                }
            } elsif (3==$numMutTypes) {
                if (exists $muttypes{DEL} && exists $muttypes{INS} && exists $muttypes{ROT}) {
                    # deletion, insertion and rotation combination is fine
                } else {
                    die "Sorry, but multiple mutations on a single contig has not been implemented.\n";
                }
            } else {
                die "Sorry, but multiple mutations on a single contig has not been implemented.\n";
            }
        }
        # if (scalar(keys %muttypes)>1) {
        #     die "Sorry, but multiple mutations on a single contig has not been implemented.\n";
        # } elsif (exists $muttypes{DEL}) {
        #     # composite deletion supported
        # } else {
        #     die "Sorry, but multiple mutations on a single contig has not been implemented.\n";
        # }
    }
    if (0!=$nonVCFDEL || 0!=$VCFDEL || 0!=$cratedRotation || 0!=$VCFINS || 0!=$nonVCFINS) {
        # handle nonVCF deletion and/or multiple deletion alone
        # + possibly rotation
        my $seenRotation = 0;
        my @rotations = ();
        my @bases = split(//, $pathSeq);
        foreach my $mut (@{$pathRef->{mutations}}) {
            if ($mut =~ /del\.\d+\-\d+/) {
                die "See mutation $mut after rotation. Unsupported feature!\n" if (0!=$seenRotation);
                my ($start, $end) = $mut =~ /del\.(\d+)\-(\d+)/;
                for(my $i=$start-1; $i<$end; ++$i) {
                    $bases[$i] = '';
                }
            } elsif ($mut =~ /ins\.\d+\..+/) {
                # TODO: perform insertion
                my ($start, $seqOrFile) = $mut =~ /ins\.(\d+)\.(.+)/;
                if (-f $seqOrFile) {
                    my $fileSeq = '';
                    open SEQFILE, $seqOrFile || die "Fail to open $seqOrFile\n$!\n";
                    while (<SEQFILE>) {
                        chomp();
                        if (/^>/) {
                            die "Multiple sequences in $seqOrFile!\n" if ('' ne $fileSeq);
                            $fileSeq = '';
                        } else {
                            $fileSeq .= $_;
                        }
                    }
                    close SEQFILE;
                    $seqOrFile = $fileSeq;
                }
                $bases[$start-1] .= $seqOrFile;
            } elsif ($mut =~ /rotate\.[+\-]*\d+/) {
                $seenRotation = 1;
                push @rotations, $mut;
            } else {
                die "See variant $mut after rotation. Unsupported feature!\n" if (0!=$seenRotation);
                my $vcfRef = $variants{$mut};
                my $svType = $vcfRef->{INFO}->{SVTYPE};
                my $svLen = $vcfRef->{INFO}->{SVLEN};
                if ('DEL' eq $svType) {
                    for(my $i=$vcfRef->{POS}-1; $i<($vcfRef->{INFO}->{END}-1); ++$i) {
                        $bases[$i] = '';
                    }
                } elsif ('INS' eq $svType) {
                    die "INS expects $vcfRef->{REF} at $vcfRef->{POS} but sees $bases[$vcfRef->{POS}-1]\n" if ($bases[$vcfRef->{POS}-1] ne $vcfRef->{REF});
                    $bases[$vcfRef->{POS}-1] = $vcfRef->{ALT};
                }
            }
        }
        $pathSeq = join('', @bases);
        # TODO: to rotate if necessary!
        foreach my $mut (@rotations) {
            my ($rotateBp) = $mut =~ /rotate\.([+\-]*\d+)/;
            $rotateBp = int($rotateBp);
            if ($rotateBp>0) {
                if ($rotateBp>length($pathSeq)) {
                    $rotateBp = $rotateBp % length($pathSeq);
                }
                my $segL = substr($pathSeq, 0, $rotateBp);
                my $segR = substr($pathSeq, $rotateBp);
                $pathSeq = $segR . $segL;
            } else {
                die "Unsupported negative rotation $mut\n";
            }
        }
    } else {
        foreach my $mut (@{$pathRef->{mutations}}) {
            my $vcfRef = $variants{$mut};
            my $svType = $vcfRef->{INFO}->{SVTYPE};
            # CHROM POS INFO/END
            # sniffles/2.2 has 1-based start, its END is 1 bp after SV
            my $svLen = $vcfRef->{INFO}->{SVLEN};
            # my $svLen = $vcfRef->{INFO}->{END} - $vcfRef->{POS} + 1;
            if ('INV' eq $svType) {
                # 3 segments, rc, stitch back
                my $segL = substr($pathSeq, 0, $vcfRef->{POS}-1);
                my $segM = substr($pathSeq, $vcfRef->{POS}-1, $svLen);
                my $segR = substr($pathSeq, $vcfRef->{INFO}->{END}-1);
                ## print $vcfRef->{INFO}->{SVLEN}, " vs ", ($vcfRef->{INFO}->{END} - $vcfRef->{POS} + 1), "\n";
                ##### FIXME: use alignment to figure out!
                ## die "wrong implementation for inv! new(",length($segL.$segM.$segR),") vs original(",length($pathSeq),")" if ($pathSeq ne ($segL.$segM.$segR));
                printf "segL(%d) + segM(%d) + segR(%d) = delta(%d) vs SVLEN = %d\n", 
                    length($segL), length($segM), length($segR), 
                    (length($pathSeq)-length($segL.$segR)), abs($svLen);
                $segM = reverse $segM;
                $segM =~ tr/ACGTNacgtn/TGCANtgcan/;
                $pathSeq = $segL . $segM . $segR;
            } elsif ('DEL' eq $svType) {
                # 2 segments, stitch together
                my $segL = substr($pathSeq, 0, $vcfRef->{POS}-1);
                my $segR = substr($pathSeq, $vcfRef->{INFO}->{END}-1);
                printf "| segL + segR | = %d vs SVLEN = %d\n", (length($pathSeq)-length($segL.$segR)), abs($svLen);
                ##### FIXME: use alignment to figure out!
                ## die "wrong implementation for del! new(",(length($pathSeq)-length($segL.$segR)),") vs original(",abs($vcfRef->{INFO}->{SVLEN}),")" if ((length($pathSeq)-length($segL.$segR))!=(abs($svLen)+1));
                $pathSeq = $segL . $segR;
            } elsif ('DUP' eq $svType) {
                # 3 segments, rc, stitch back
                my $segL = substr($pathSeq, 0, $vcfRef->{POS}-1);
                my $segM = substr($pathSeq, $vcfRef->{POS}-1, $svLen);
                my $segR = substr($pathSeq, $vcfRef->{INFO}->{END}-1);
                ## print $vcfRef->{INFO}->{SVLEN}, " vs ", ($vcfRef->{INFO}->{END} - $vcfRef->{POS} + 1), "\n";
                ##### FIXME: figure out how to check
                ## die "wrong implementation for inv! new(",length($segL.$segM.$segR),") vs original(",length($pathSeq),")" if ($pathSeq ne ($segL.$segM.$segR));
                printf "segL(%d) + segM(%d) + segM(%d) + segR(%d) = delta(%d) vs SVLEN = %d\n", 
                    length($segL), length($segM), length($segM), length($segR), 
                    abs(length($pathSeq)-length($segL.$segM.$segM.$segR)), abs($svLen);
                $pathSeq = $segL . $segM . $segM . $segR;
            } else {
                die "Sorry, but svType $svType has not been implemented.\n";
            }
        }
    }

    if (0!=$toRCSeq) {
        $pathSeq = reverse $pathSeq;
        $pathSeq =~ tr/ACGTNacgtn/TGCANtgcan/;
    }
    $finalSeq .= $pathSeq;
}

my @ids = (); grep { push @ids, $_->{seqid} } @paths;
my $comments = sprintf("path=%s", join(",",@ids));
if (scalar(keys %variants)>0) {
    $comments .= sprintf(" variants=%s", join(",", sort keys %variants));
}
$comments .= sprintf(" length=%d", length($finalSeq));
if ($rotate>0) {
    if ($rotate>length($finalSeq)) {
        $rotate = $rotate % length($finalSeq);
    }
    $comments .= sprintf(" rotate=%d", $rotate);

    my $segL = substr($finalSeq, 0, $rotate);
    my $segR = substr($finalSeq, $rotate);
    $finalSeq = $segR . $segL;
}

open OUTFILE, ">$prefix.fa" || die "fail to open $prefix.fa\n$!\n";
printf OUTFILE ">%s %s\n", $pathid, $comments;
print OUTFILE $finalSeq, "\n";
close OUTFILE;

my $rcSeq = reverse $finalSeq;
$rcSeq =~ tr/ACGTNacgtn/TGCANtgcan/;
open OUTFILE, ">$prefix.rc.fa" || die "fail to open $prefix.rc.fa\n$!\n";
printf OUTFILE ">%s/rc %s\n", $pathid, $comments;
print OUTFILE $rcSeq, "\n";
close OUTFILE;

exit 0;


