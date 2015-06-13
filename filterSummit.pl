#! /usr/bin/perl
#
use strict;
use warnings;

my $summitInstance = shift;
my $filterInstance = shift;
my $window = shift;

&main ();

sub main
{
    $window = 10 if ( not defined $window );
    my $nrSummitInstance = $summitInstance . ".nr";
    rmRdp ( $summitInstance, $nrSummitInstance );

    filter ( $nrSummitInstance, $filterInstance, $window );
}

sub rmRdp 
{
    my $summitInst = shift;
    my $nrSummitInstance = shift;

    my $summitInstTmp = $summitInst . ".tmp";
    print STDERR `sort -k1,1 -k2,3n $summitInst -o $summitInstTmp`;

    my $count = 0;
    my $lastLine = "NULL";
    my $lastCount = 0;
    open ( SMT, $summitInstTmp );
    open ( NR, ">$nrSummitInstance" );
    while ( my $line = <SMT> ) {
        next if ( $line =~ /^#/ );
        chomp $line;

        $count++;
        #last if ( $count > 50 );

        ## new site
        if ( $lastLine ne $line ) {
            if ( $lastLine ne "NULL" ) {
                my ( $chr, $start, $end, $ensembl, $tag, $strand ) = split ( /\t/, $lastLine );
                print NR join "\t", $chr, $start, $end, $ensembl, $lastCount, $strand;
                print NR "\n";
            }

            $lastLine = $line;
            $lastCount = 1;
        }
        else { $lastCount++; }
    }
    close SMT;
    my ( $chr, $start, $end, $ensembl, $tag, $strand ) = split ( /\t/, $lastLine );
    print NR join "\t", $chr, $start, $end, $ensembl, $lastCount, $strand;
    print NR "\n";
    close NR;

    print STDERR `/bin/rm $summitInstTmp`;

    1;
}

sub filter
{
    my $nrSummitInstance = shift;
    my $outSummit = shift;
    my $win = shift;

    my $count = 0;
    my %cluster = ();

    my $lastchr = "NULL";
    my $lastEnd = 0;
    my $lastStrand = "NULL";
    open ( SMT, $nrSummitInstance );
    open ( OUT, ">$outSummit" );
    while ( my $line = <SMT> ) {
        next if ( $line =~ /^#/ );
        chomp $line;

        $count++;
        #last if ( $count > 50 );

        my ( $chr, $start, $end, $ensembl, $height, $strand ) = split ( /\t/, $line );

        ## new cluster
        if ( $lastchr ne $chr ) {
            if ( $lastchr ne "NULL" ) {
                foreach my $strand ( "+", "-" ) {
                    if ( $cluster{$strand}{chrom} ne "0" ) {
                        print OUT $cluster{$strand}{chrom}, "\t", $cluster{$strand}{start}, "\t", $cluster{$strand}{end};
                        print OUT "\t", $strand, "\t", $cluster{$strand}{height}, "\n";
                    }
                }
            }

            %cluster = ();
            $cluster{"+"}{chrom} = 0;
            $cluster{"-"}{chrom} = 0;

            $cluster{$strand}{chrom} = $chr;
            $cluster{$strand}{start} = $start;
            $cluster{$strand}{end} = $end;
            $cluster{$strand}{height} = $height;

            $lastchr = $chr;
            $lastStrand = $strand;
            $lastEnd = $end;
        }
        elsif ( ( $strand ne $lastStrand ) or ( $start - $lastEnd ) > $win ) {
            foreach my $strand ( "+", "-" ) {
                if ( $cluster{$strand}{chrom} ne "0" ) {
                    print OUT $cluster{$strand}{chrom}, "\t", $cluster{$strand}{start}, "\t", $cluster{$strand}{end};
                    print OUT "\t.\t", $cluster{$strand}{height}, "\t", $strand, "\n";
                }
            }

            %cluster = ();
            $cluster{"+"}{chrom} = 0;
            $cluster{"-"}{chrom} = 0;

            $cluster{$strand}{chrom} = $chr;
            $cluster{$strand}{start} = $start;
            $cluster{$strand}{end} = $end;
            $cluster{$strand}{height} = $height;

            $lastchr = $chr;
            $lastStrand = $strand;
            $lastEnd = $end;
        }
        else {
            if ( $height > $cluster{$strand}{height} ) {
                $cluster{$strand}{chrom} = $chr;
                $cluster{$strand}{start} = $start;
                $cluster{$strand}{end} = $end;
                $cluster{$strand}{height} = $height;
            }

            $lastEnd = $end;
        }
    }
    close SMT;

    foreach my $strand ( "+", "-" ) {
        if ( $cluster{$strand}{chrom} ne "0" ) {
            print OUT $cluster{$strand}{chrom}, "\t", $cluster{$strand}{start}, "\t", $cluster{$strand}{end};
            print OUT "\t.\t", $cluster{$strand}{height}, "\t", $strand, "\n";
        }
    }
    close OUT;

    1;
}

