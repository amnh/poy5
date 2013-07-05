use HTML::HTML5::Parser;

sub getContent {
    my ($secname,$node) = @_;

    #skip unhelpful nodes
    while( $node->textContent !~ $secname ){
        $node2 = $node->nextNonBlankSibling();
        if ( defined $node2 ){
            $node = $node2;
        } else {
            return $node;
        }
    }
    $node = $node->nextNonBlankSibling();
    $node = $node->nextNonBlankSibling();
    $node = $node->nextNonBlankSibling();
    return $node;
}

sub clean_name {
    my ($mess) = @_;
    $mess =~ /([A-Za-z_]+)/s;
    return $1;
}

#return a Document
my $parser = HTML::HTML5::Parser->new;
my $dom    = $parser->parse_file( shift );

#deal with output file
my $file = shift;
my $outfile;
open $outfile, '>>', $file;

#get elements by id -- h3 are the tag names.
my @elts   = $dom->getElementsByTagName("h3");
shift @elts;
shift @elts;

#iterate through each section gathering description and examples.
print {$outfile} "let index = [\n";
foreach $elt (@elts){
    $n = clean_name($elt->textContent);

    $d = getContent("Description",$elt)->textContent;
    $d =~ s/^\s+//; #remove leading spaces
    $d =~ s/\s+$//; #remove trailing spaces
    $d =~ s/"/\\"/g; #replace quotes
    $d =~ s/\\n/\n/g; #replace newlines

    $en= getContent("Example",$elt);
    $e = $en->textContent;
    $e =~ s/^\s+//; #remove leading spaces
    $e =~ s/\s+$//; #remove trailing spaces
    $e =~ s/"/\\"/g; #replace quotes

    print {$outfile} "\t(\"$n\",(\"$d\",\"$e\"));\n";
}
print {$outfile} "]";
