/* Parser for hennig files */
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
let report_error b e =
    let b = string_of_int (b.Lexing.pos_cnum)
    and e = string_of_int (e.Lexing.pos_cnum) in
    Status.user_message Status.Error 
    ("Unrecognized@ command@ between@ characters@ " ^ b ^ "@ and@ "
    ^ e)
%}
%token <string> DATA TREES
%token <string> TREAD
%token <string> WORD
%token <string> INT
%token <string> CHARNAME
%token <char> CHAR
%token CCODE COST PROCESS OPTCODE CHARNAMECMD NSTATES DNA PROTEINS NUMBER GAP NOGAP
%token LPARENT RPARENT GT EQUAL QUESTION SEMICOLON DASH LSQ RSQ PLUS STAR BACKSLASH LBRACKET RBRACKET DOT 
%type <Hennig.command> command
%start command
%type <(int * int * string)> xread
%start xread
%%

command:
    | DATA SEMICOLON { Hennig.Xread $1 }
    | CCODE character_change_list SEMICOLON { Hennig.Ccode $2 }
    | TREES SEMICOLON { Hennig.Tread $1 }
    | COST cost_change_list SEMICOLON { Hennig.Cost $2 }
    | PROCESS BACKSLASH SEMICOLON { Hennig.Ignore }
    | OPTCODE INT DOT INT SEMICOLON { Hennig.Ignore }
    | CHARNAMECMD char_names_list SEMICOLON { Hennig.Charname $2 }
    | NSTATES number_of_states { Hennig.Nstates (Some $2) }
    | error SEMICOLON { 
        report_error (Parsing.symbol_start_pos ()) (Parsing.symbol_end_pos ());
        Hennig.Ignore 
    }

gap:
    |           { None }
    | GAP       { Some `Gap }
    | NOGAP     { Some `Nogap }

number_of_states:
    | STAR          { `Number 8 }
    | DNA gap       { `Dna $2 }
    | PROTEINS gap  { `Protein $2 }
    | NUMBER INT    { `Number (int_of_string $2) }
char_names_list:
    | CHARNAME char_names_list { 
        let res = 
            (* Get rif of the closing semicolon *)
            let res = $1 in
            String.sub res 0 ((String.length res) - 1)
        in
        res :: $2 }
    | { [] }
character_change_list:
    | character_change character_change_list { $1 :: $2 }
    | { [] }
character_change:
    | PLUS character_list   { Hennig.Additive $2 }
    | DASH character_list   { Hennig.NonAdditive $2 }
    | LSQ character_list    { Hennig.Active $2 }
    | RSQ character_list    { Hennig.Inactive $2 }
    | LPARENT character_list { Hennig.Sankoff $2 }
    | RPARENT character_list { Hennig.NonAdditive $2 }
    | BACKSLASH INT character_list { Hennig.Weight (int_of_string $2, $3) }
character_list:
    | DOT { [Hennig.All] }
    | aux_character_list { $1 }
aux_character_list:
    | INT DOT INT aux_character_list { (Hennig.Range (int_of_string $1,
    int_of_string $3)) :: $4 }
    | INT aux_character_list { (Hennig.Single (int_of_string $1)) :: $2 }
    | { [] }
cost_change_list:
    | cost_change cost_change_list { $1 :: $2 }
    | { [] }
cost_change:
    | INT EQUAL chars GT chars INT 
        { (false, int_of_string $1, ($3, $5) , int_of_string $6) }
    | INT EQUAL chars BACKSLASH chars INT
        { (true, int_of_string $1, ($3, $5) , int_of_string $6) }
chars:
    | INT { $1 }
    | LSQ int_list RSQ { String.concat "" $2 }
int_list:
    | INT int_list { $1 :: $2 }
    | { [] }
xread:
    | INT INT DATA { (int_of_string $1), (int_of_string $2), $3 }
