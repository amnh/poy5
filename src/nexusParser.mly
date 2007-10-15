/* Parser for nexus files */
%{
let parse_error s = 
    try
        let b = (Parsing.symbol_start_pos ()) 
        and e = (Parsing.symbol_end_pos ()) in
        let b = string_of_int (b.Lexing.pos_cnum)
        and e = string_of_int (e.Lexing.pos_cnum) in
        Status.user_message Status.Error 
        (s ^ "@ between@ characters@ " ^ b ^ 
        "@ and@ " ^ e)
    with
    | _ -> Status.user_message Status.Error s

let report_error text b e =
    let b = string_of_int (b.Lexing.pos_cnum)
    and e = string_of_int (e.Lexing.pos_cnum) in
    Status.user_message Status.Error 
    ("There@ was@ a@ parsing@ error@ between@ characters@ " ^ b ^ "@ and@ "
    ^ e ^ " in the " ^ text)
%}
%token EOF
%token <string> ANCSTATES
%token <string> ASSUMPTIONS
%token <string> AVERAGE
%token <string> BEGIN
%token <string> BINHEX
%token <string> BOTH
%token <string> CHANGESET
%token <char> CHAR
%token <string> CHARACTER
%token <string> CHARACTERS
%token <string> CHARLABELS
%token <string> CHARPARTITION
%token <string> CHARSET
%token <string> CHARSTATELABELS
%token <string> CODEORDER
%token <string> CODESET
%token <string> CODONS
%token <string> CONDONPOSSET
%token <string> CONTINUOUS
%token <string> COUNT
%token <string> CSTREE
%token <string> DATA
%token <string> DATA_CSTREE
%token <string> QUOTED
%token <string> DATATYPE
%token <string> DEFTYPE
%token <string> DIAGONAL
%token <string> DIMENSIONS
%token <string> DISTANCES
%token <string> DNA
%token <string> DROS
%token <string> ELIMINATE
%token <string> ENCODE
%token <string> END
%token <string> EPS
%token <string> EQUATE
%token <string> EXSET
%token <string> EXT
%token <string> EXTENSIONS
%token <string> FILE
%token <string> FORMAT
%token <string> FREQUENCY
%token <string> GAP
%token <string> GAPMODE
%token <string> GENETICCODE
%token <string> GIF
%token <string> INDIVIDUALS
%token <string> INLINE
%token <string> INTERLEAVE
%token <string> ITEMS
%token <string> JPEG
%token <string> LABELS
%token <string> LOWER
%token <string> MAM
%token <string> MATCHCHAR
%token <string> MATRIX
%token <string> MAX
%token <string> MAXSTEPS
%token <string> MEDIAN
%token <string> MIN
%token <string> MINSTEPS
%token <string> MISSING
%token <string> MTDNA
%token <string> NCHAR
%token <string> NEWSTATE
%token <string> NEWTAXA
%token <string> NOLABELS
%token <string> NONE
%token <string> NOTES
%token <string> NTAX
%token <string> NOTOKENS
%token <string> NUCLEOTIDE
%token <string> NUCORDER
%token <string> OPTIONS
%token <string> PICT
%token <string> PICTURE
%token <string> POLYTCOUNT
%token <string> PROTEIN
%token <string> RESOURCE
%token <string> RESPECTCASE
%token <string> RNA
%token <string> SAMPLESIZE
%token <string> SETS
%token <string> SOURCE
%token <string> STANDARD
%token <string> STATE
%token <string> STATELABELS
%token <string> STATES
%token <string> STATESET
%token <string> STATESFORMAT
%token <string> STATESPRESENT
%token <string> STDERROR
%token <string> STEPMATRIX
%token <string> SYMBOLS
%token <string> TAXON
%token <string> TAXA
%token <string> TAXLABELS
%token <string> TAXPARTITION
%token <string> TAXSET
%token <string> TEXT
%token <string> TIFF
%token <string> TOKENS
%token <string> TRANSLATE
%token <string> TRANSPOSE
%token <string> TREE
%token <string> TREEPARTITION
%token <string> TREES
%token <string> TREESET
%token <string> TRIANGLE
%token <string> TYPESET
%token <string> UNALIGNED
%token <string> UNIVERSAL
%token <string> UPPER
%token <string> USERTYPE
%token <string> UUENCODE
%token <string> VARIANCE
%token <string> VECTOR
%token <string> WTSET
%token <string> YEAST
%token <string> EIDENT
%token NEXUS SEMICOLON EQUAL COMMA QUOTE BACKSLASH DASH LPARENT RPARENT STAR
COLON
%token <string> IDENT
%token <string> FLOAT
%token <string> INTEGER
%start tree
%type <Nexus.tree> tree
%start header
%type <unit> header
%start block
%type <Nexus.block> block
%start symbol_pair
%type <(string * string list)> symbol_pair
%start symbol_list
%type <string list> symbol_list
%%

header:
    | NEXUS { () }
    ;
block:
    | BEGIN TAXA SEMICOLON taxa END SEMICOLON  
        { Nexus.Taxa $4 }
    | BEGIN CHARACTERS SEMICOLON characters END SEMICOLON 
        { Nexus.Characters $4 }
    | BEGIN DATA SEMICOLON characters END SEMICOLON
        { Nexus.Characters $4 }
    | BEGIN UNALIGNED SEMICOLON unaligned END SEMICOLON
        { Nexus.Unaligned $4 }
    | BEGIN TREES SEMICOLON trees END SEMICOLON
        { Nexus.Trees $4 }
    | BEGIN NOTES SEMICOLON notes END SEMICOLON 
        { Nexus.Notes $4 }
    | BEGIN DISTANCES SEMICOLON distances END SEMICOLON 
        { Nexus.Distances $4 }
    | BEGIN ASSUMPTIONS SEMICOLON assumptions END SEMICOLON
        { Nexus.Assumptions $4 }
    | BEGIN SETS SEMICOLON list_of_anything END SEMICOLON
        { Nexus.Ignore $2 }
    | BEGIN IDENT error END SEMICOLON
        { Nexus.Error $2 }
    | BEGIN TAXA SEMICOLON error END SEMICOLON  
        { Nexus.Error $2 }
    | BEGIN UNALIGNED SEMICOLON error END SEMICOLON
        { Nexus.Error $2 }
    | BEGIN NOTES SEMICOLON error END SEMICOLON 
        { Nexus.Error $2 }
    | BEGIN DISTANCES SEMICOLON error END SEMICOLON 
        { Nexus.Error $2 }
    | BEGIN SETS SEMICOLON error END SEMICOLON
        { Nexus.Error $2 }
    ;
assumptions:
    | assumption_items assumptions { $1 :: $2 }
    | { [] }
    ;
assumption_items:
    | optional_assumption_options   { Nexus.Options $1 }
    | optional_user_type            { Nexus.UserType $1 }
    | optional_type_set             { Nexus.TypeDef $1 }
    | optional_wtset                { Nexus.WeightDef $1 }
    | optional_exset                { Nexus.ExcludeSet $1 }
    | optional_ancstates            { Nexus.AncestralDef $1 }
    | error                         
    { report_error "Assumption Block" (Parsing.symbol_start_pos ()) (Parsing.symbol_end_pos ());
        raise Parsing.Parse_error }
    ;
optional_assumption_options:
    | OPTIONS deftype polytcount gapmode SEMICOLON { ($2, $3, $4) }
    ;
deftype:
    | DEFTYPE EQUAL IDENT { Some $3 }
    | { None }
    ;
polytcount:
    | POLYTCOUNT EQUAL MINSTEPS { Nexus.MinSteps }
    | POLYTCOUNT EQUAL MAXSTEPS { Nexus.MaxSteps }
    | { Nexus.MinSteps }
    ;
gapmode:
    | GAPMODE EQUAL MISSING { Nexus.Missing }
    | GAPMODE EQUAL NEWSTATE { Nexus.NewState }
    | { Nexus.NewState }
    ;
optional_user_type:
    | USERTYPE IDENT user_type_definition SEMICOLON { ($2, $3) }
    ;
user_type_definition:
    | EQUAL INTEGER numbers_and_chars { Nexus.StepMatrix ($2, $3) }
    | STEPMATRIX EQUAL INTEGER numbers_and_chars 
        { Nexus.StepMatrix ($3, $4) }
    | CSTREE DATA_CSTREE { Nexus.CSTree $2 }
    ;
numbers_and_chars:
    | number_and_char numbers_and_chars { $1 :: $2 }
    | number_and_char { [$1] }
    ;
number_and_char:
    | INTEGER   { $1 }
    | FLOAT     { $1 }
    | CHAR      { Char.escaped $1 }
    | COLON     { ":" }
    | EQUAL     { "=" }
    | COMMA     { "," }
    | BACKSLASH { "/" }
    | DASH      { "-" }
    | STAR      { "*" }
    | IDENT     { $1 }
    ;
optional_type_set:
    | TYPESET optional_set_for_assumptions { $2 }
    ;
optional_wtset:
    | WTSET optional_set_for_assumptions { $2 }
    ;
do_token:
    | NOTOKENS { false }
    | TOKENS { true }
    | { false }
do_star:
    | STAR { false }
    | { true }
optional_set_for_assumptions:
    | do_star IDENT do_token EQUAL standard_type_set SEMICOLON 
        { ($1, $2, $3, Nexus.Standard $5) }
    | do_star IDENT STANDARD do_token EQUAL standard_type_set SEMICOLON 
        { ($1, $2, $4, Nexus.Standard $6) }
    | do_star IDENT VECTOR do_token EQUAL vector_type_set SEMICOLON 
        { ($1, $2, $4,Nexus.Vector $6) }
    ;
standard_type_set_item:
    | INTEGER COLON characterset_list { Nexus.Code ($1, $3) }
    | FLOAT COLON characterset_list   { Nexus.Code ($1, $3) }
    | IDENT COLON characterset_list   { Nexus.IName ($1, $3) }
    ;
standard_type_set:
    | standard_type_set_item COMMA standard_type_set { ($1 :: $3) }
    | standard_type_set_item { [$1] }
    ;
vector_type_set:
    | INTEGER vector_type_set { ($1 :: $2) }
    | INTEGER { [$1] }
    ;
optional_exset:
    | EXSET optional_set_for_assumptions { $2 }
    ;
optional_ancstates:
    | ANCSTATES optional_set_for_assumptions { $2 }
    ;
distances:
    | optional_distances_dimensions optional_format optional_taxlabels DATA SEMICOLON
        { ($1, $2, $3, $4) }
    ;
optional_distances_dimensions:
    | DIMENSIONS NEWTAXA NTAX EQUAL INTEGER NCHAR EQUAL INTEGER SEMICOLON
        { Some (true, $5, $8) }
    | DIMENSIONS NTAX EQUAL INTEGER NCHAR EQUAL INTEGER SEMICOLON
        { Some (false, $4, $7) }
    | { None }
    ;
notes:
    | optional_text optional_picture 
        { $1, $2 }
    ;
optional_text:
    | TEXT optional_set_pair_list SOURCE EQUAL source DATA SEMICOLON
        { Some ($2, $5, $6) }
    | { None }
    ;
optional_picture:
    | PICTURE optional_set_pair_list optional_pictureformat optional_encode SOURCE EQUAL source 
    DATA SEMICOLON { Some ($2, $3, $4, $7, $8) }
    | { None } 
    ;
optional_set_pair_list:
    | set_pair optional_set_pair_list { $1 :: $2 }
    | { [] }
    ;
set_pair:
    | TAXON EQUAL characterset     { Nexus.TaxonSet $3 }
    | CHARACTER EQUAL characterset { Nexus.CharacterSet $3 }
    | STATE EQUAL characterset      { Nexus.StateSet $3 }
    | TREE EQUAL characterset        { Nexus.TreeSet $3 }
    ;
source:
    | INLINE    {Nexus.Inline }
    | FILE      { Nexus.File }
    | RESOURCE  { Nexus.Resource }
    ;
optional_pictureformat:
    | FORMAT EQUAL pictureformat { Some $3 }
    | { None }
    ;
pictureformat:
    | PICT      { Nexus.Pict }
    | TIFF      { Nexus.Tiff }
    | EPS       { Nexus.Eps }
    | JPEG      { Nexus.Jpeg }
    | GIF       { Nexus.Gif }
    ;
optional_encode:
    | ENCODE EQUAL pictureencoding { Some $3 }
    | { None }
    ;
pictureencoding:
    | NONE  { Nexus.None }
    | UUENCODE { Nexus.UUEncode }
    | BINHEX    { Nexus.BinHex }
    ;
trees:
    | optional_translate tree_list { ($1, $2) }
    ;
optional_translate:
    | TRANSLATE pairs_list SEMICOLON { $2 }
    | { [] }
    ;
tree_list:
    | DATA SEMICOLON tree_list { $1 :: $3 }
    | DATA SEMICOLON { [ $1 ] }
    ;
pairs_list:
    | IDENT IDENT COMMA pairs_list { ($1, $2) :: $4 }
    | IDENT IDENT { [$1, $2] }
    ;
characters:
    | DIMENSIONS optional_taxa_dimensions NCHAR EQUAL INTEGER SEMICOLON 
     optional_format optional_eliminate optional_taxlabels
     optional_charstatelabels optional_charlabels optional_statelabels 
     DATA SEMICOLON 
     { 
         { Nexus.char_taxon_dimensions = $2;
         Nexus.char_char_dimensions = $5;
         Nexus.char_format = $7;
         Nexus.char_eliminate = $8;
         Nexus.char_taxlabels = $9;
         Nexus.char_statelabels = $10;
         Nexus.char_charlabels = $11;
         Nexus.char_charstates = $12;
         Nexus.chars = $13;}
     }
unaligned:
    | optional_unaligned_dimensions optional_format DATA SEMICOLON
        { { Nexus.unal_taxon_dimensions = $1; unal_format = $2; unal = $3 } }
     ;
optional_unaligned_dimensions:
    | DIMENSIONS optional_taxa_dimensions { $2 }
    | { None }
    ;
optional_charlabels:
    | CHARLABELS taxonlist SEMICOLON    { $2 }
    |                                   { [] }
    ;
optional_charstatelabels:
    | CHARSTATELABELS charstatelables SEMICOLON { $2 }
    |                                           { [] }
    ;
charstatelables:
    | INTEGER IDENT BACKSLASH taxonlist COMMA charstatelables { ($1, $2, $4) :: $6 }
    | INTEGER IDENT BACKSLASH taxonlist { [] }
    ;
optional_statelabels:
    | STATELABELS statelabels SEMICOLON { $2 }
    |                           { [] }
    ;
statelabels:
    | INTEGER taxonlist COMMA statelabels { ($1, $2) :: $4 }
    |                                       { [] }
    ;
optional_taxlabels:
    | TAXLABELS taxonlist SEMICOLON { $2 }
    |                           { [] }
    ;
optional_eliminate:
    | ELIMINATE characterset SEMICOLON { Some $2 }
    |                           { None }
    ;
optional_format:
    | FORMAT format_items_list SEMICOLON { $2 }
    |           { [] }
    ;
format_items_list:
    | format_items format_items_list { $1 :: $2 }
    | error format_items_list { 
        report_error "Format list" (Parsing.symbol_start_pos ()) (Parsing.symbol_end_pos ());
        $2 }
    |           { [] }
    ;
format_items:
    | DATATYPE EQUAL datatype { Nexus.Datatype $3 }
    | DATATYPE EQUAL error {
        Status.user_message Status.Error
        ("Many@ \"NEXUS\"@ files@ include@ in@ the@ CHARACTERS@ or@ DATA@ " ^
        "block@ the@ dataype@ @[<b>mixed@]@ or@ @[<b>restriction@]@ or@ something@ else.@ " ^
        "@ Those@ are@ most@ likely@ Mr.Bayes@ files.@ No@ matter@ that@ " ^
        "they@ say@ #NEXUS@ at@ the@ beginning,@ this@ is@ not@ a@ valid@ " ^
        "NEXUS@ file.@ I can@ not@ process@ this.@ I@ am@ sorry. " ^
        "You'll@ have@ to@ convert@ it@ to@ a@ valid@ nexus@ file.");
        raise Parsing.Parse_error }
    | RESPECTCASE { Nexus.RespectCase }
    | MISSING EQUAL symbol { Nexus.FMissing $3 }
    | GAP EQUAL symbol { Nexus.Gap $3 }
    | SYMBOLS EQUAL QUOTED { 
        let len = String.length $3 in
        Nexus.Symbols (String.sub $3 0 (len - 1))}
    | EQUATE EQUAL QUOTED { Nexus.Equate $3 }
    | MATCHCHAR EQUAL symbol { Nexus.MatchChar $3 }
    | NOLABELS     { Nexus.Labels false }
    | LABELS        { Nexus.Labels true }
    | TRANSPOSE     { Nexus.Transpose }
    | INTERLEAVE    EQUAL { 
        Status.user_message Status.Error
        ("Many@ \"NEXUS\"@ files@ include@ in@ the@ CHARACTERS@ or@ DATA@ " ^
        "block@ the@ directives@ interleave=yes@ or@ interleave=no.@ " ^
        "@[<b>That@ is@ not@ a@ valid@ NEXUS@ directive@].@ To@ fix@ your@ file@ " ^
        "replace@ interleave=yes@ with@ INTERLEAVE,@ or@ just@ remove@ it@ " ^
        "if@ you@ see@ interleave=no.");
        raise Parsing.Parse_error }
    | INTERLEAVE    { Nexus.Interleave }
    | ITEMS EQUAL item { Nexus.Items $3 }
    | STATESFORMAT EQUAL states_format { Nexus.StatesFormat $3 }
    | NOTOKENS     { Nexus.Tokens false }
    | TOKENS        { Nexus.Tokens true }
    | TRIANGLE EQUAL triangle_format { Nexus.Triangle $3 }
    ;
triangle_format:
    | LOWER { Nexus.Lower }
    | UPPER { Nexus.Upper }
    | BOTH { Nexus.Both }
    ;
datatype:
    | STANDARD { Nexus.DStandard }
    | DNA { Nexus.Dna }
    | RNA { Nexus.Rna }
    | NUCLEOTIDE { Nexus.Nucleotide }
    | PROTEIN { Nexus.Protein }
    | CONTINUOUS { Nexus.Continuous }
    ;
symbol:
    | IDENT { $1 }
    | DASH { "-" }
    | CHAR  { Char.escaped $1 }
    ;
symbol_list:
    | symbol symbol_list { $1 :: $2 }
    |       { [] }
    ;
symbol_pair:
    | symbol EQUAL symbol {($1, [$3])}
    | symbol EQUAL LPARENT symbol_list RPARENT { ($1, $4) }
    ;
item:
    | MIN { Nexus.Min }
    | MAX { Nexus.Max }
    | MEDIAN  { Nexus.Median }
    | AVERAGE { Nexus.Average }
    | VARIANCE { Nexus.Variance }
    | STDERROR { Nexus.Stderror }
    | SAMPLESIZE { Nexus.SampleSize }
    | STATES     { Nexus.States }
    ;
states_format:
    | STATESPRESENT { Nexus.StatesPresent }
    | INDIVIDUALS   { Nexus.Individuals }
    | COUNT         { Nexus.Count }
    | FREQUENCY     { Nexus.Frequency }
    ;
optional_taxa_dimensions:
    | NTAX EQUAL INTEGER { Some $3 }
    |               { None }
    ;
taxa:
    | DIMENSIONS NTAX EQUAL INTEGER SEMICOLON TAXLABELS taxonlist SEMICOLON 
        { ($4, $7) }
    ;
taxonlist:
    | IDENT taxonlist   { $1 :: $2 }
    |                   { [] }
    ;
characterset_list:
    | characterset characterset_list { $1 :: $2 }
    | characterset      { [$1] }
    ;
characterset:
    | INTEGER DASH CHAR    { Nexus.Range ($1, None) }
    | INTEGER DASH INTEGER { Nexus.Range ($1, Some $3) }
    | INTEGER              { Nexus.Single ($1) }
    | IDENT                { Nexus.Name $1 }
    | CHAR                 { Nexus.Name (String.make 1 $1) }
    ;
list_of_anything:
    | any_thing_minus_end list_of_anything { () }
    | any_thing_minus_end { () }
    ;
any_thing_minus_end:
    | TAXA  { () }
    | DATA  { () }
    | UNALIGNED  { () }
    | TREES  { () }
    | NOTES  { () }
    | DISTANCES  { () }
    | ASSUMPTIONS  { () }
    | SETS  { () }
    | CHARACTERS { () }
    | ANCSTATES { () }
    | AVERAGE { () }
    | BINHEX { () }
    | BOTH { () }
    | CHANGESET { () }
    | CHAR { () }
    | CHARACTER { () }
    | CHARLABELS { () }
    | CHARPARTITION { () }
    | CHARSET { () }
    | CHARSTATELABELS { () }
    | CODEORDER { () }
    | CODESET { () }
    | CODONS { () }
    | CONDONPOSSET { () }
    | CONTINUOUS { () }
    | COUNT { () }
    | CSTREE { () }
    | DATA_CSTREE { () }
    | QUOTED { () }
    | DATATYPE { () }
    | DEFTYPE { () }
    | DIAGONAL { () }
    | DIMENSIONS { () }
    | DNA { () }
    | DROS { () }
    | ELIMINATE { () }
    | ENCODE { () }
    | EPS { () }
    | EQUATE { () }
    | EXSET { () }
    | EXT { () }
    | EXTENSIONS { () }
    | FILE { () }
    | FORMAT { () }
    | FREQUENCY { () }
    | GAP { () }
    | GAPMODE { () }
    | GENETICCODE { () }
    | GIF { () }
    | INDIVIDUALS { () }
    | INLINE { () }
    | INTERLEAVE { () }
    | ITEMS { () }
    | JPEG { () }
    | LABELS { () }
    | LOWER { () }
    | MAM { () }
    | MATCHCHAR { () }
    | MATRIX { () }
    | MAX { () }
    | MAXSTEPS { () }
    | MEDIAN { () }
    | MIN { () }
    | MINSTEPS { () }
    | MISSING { () }
    | MTDNA { () }
    | NCHAR { () }
    | NEWSTATE { () }
    | NEWTAXA { () }
    | NONE { () }
    | NTAX { () }
    | NUCLEOTIDE { () }
    | NUCORDER { () }
    | OPTIONS { () }
    | PICT { () }
    | PICTURE { () }
    | POLYTCOUNT { () }
    | PROTEIN { () }
    | RESOURCE { () }
    | RESPECTCASE { () }
    | RNA { () }
    | SAMPLESIZE { () }
    | SOURCE { () }
    | STANDARD { () }
    | STATE { () }
    | STATELABELS { () }
    | STATES { () }
    | STATESET { () }
    | STATESFORMAT { () }
    | STATESPRESENT { () }
    | STDERROR { () }
    | STEPMATRIX { () }
    | SYMBOLS { () }
    | TAXON { () }
    | TAXLABELS { () }
    | TAXPARTITION { () }
    | TAXSET { () }
    | TEXT { () }
    | TIFF { () }
    | TOKENS { () }
    | TRANSLATE { () }
    | TRANSPOSE { () }
    | TREE { () }
    | TREEPARTITION { () }
    | TREESET { () }
    | TRIANGLE { () }
    | TYPESET { () }
    | UNIVERSAL { () }
    | UPPER { () }
    | USERTYPE { () }
    | UUENCODE { () }
    | VARIANCE { () }
    | VECTOR { () }
    | WTSET { () }
    | YEAST { () }
    | EIDENT { () }
    | STAR { () }
    | COLON { () }
    | IDENT { () }
    | FLOAT { () }
    | INTEGER { () }
    | SEMICOLON  { () }
    | EQUAL { () }
    | COMMA { () }
    | QUOTE { () }
    | BACKSLASH { () }
    | DASH { () }
    | LPARENT { () }
    | RPARENT { () }
    ;
tree:
    | do_star IDENT EQUAL single_tree EOF { $4 }
    ;
single_tree:
    | IDENT optional_length { Nexus.Leaf ($1, $2) }
    | LPARENT single_tree_list RPARENT optional_label optional_length 
                            { Nexus.Node ($2, $4, $5) }
    ;
single_tree_list:
    | single_tree COMMA single_tree_list { $1 :: $3 }
    | single_tree { [$1] }
    ;
optional_length:
    | COLON INTEGER { Some (float_of_string $2) }
    | COLON FLOAT { Some (float_of_string $2) }
    | { None }
    ;
optional_label:
    | IDENT     { Some $1 }
    | { None }
    ;
